/* BlueMolecule — file-format readers and format detection.
 *
 * Public API (all on the global BM object):
 *   BM.detectFormat(fileName, text)  -> format id ('gaussian-log', 'orca-out', ...) or null
 *   BM.parseFile(fileName, text)     -> Dataset (throws Error with a user-readable message)
 *   BM.parseLines(fileName, lines)   -> Dataset from a line array (big files, see bigfile.js); readers take text or lines
 *   BM.parseSource(fileName, blob|chunks, opts) -> Promise<Dataset>  DCD/TRR/XTC read in slices (stride, frame range)
 *   BM.parsers.{gaussianLog, gaussianInput, orca, xyz, cube, molden, fchk, molsdf, pdb, gro}(text, fileName)
 *   BM.zmatrixToCartesian(rows)      -> Float64Array(3N)   (see the comment above the function)
 *   BM.formats                       -> list of {id, name, extensions}
 *
 * Dataset shape: see newDataset(). Lengths are Angstrom, energies Hartree, dipoles Debye,
 * frequencies cm-1. Written from public descriptions of the file formats only.
 *
 * QM/MM and ONIOM outputs (Gaussian ONIOM, ORCA %qmmm / QM/XTB / QM/QM2, CP2K QMMM, xtb --oniom) set
 * ds.meta.qmmm = true and
 *   ds.layers = { scheme: 'oniom' | 'qmmm' | 'qm-qm2',
 *     atomLayer: 'H' (high / QM) | 'M' (medium / QM2) | 'L' (low / MM) per full-system atom (null if the output never says),
 *     frozen, active: boolean[] | null,  linkAtoms: [{atom, bondedTo, element, type}] (0-based; atom = the atom the
 *     link atom replaces, bondedTo = the QM-side atom), mmTypes: string[] | null, mmCharges: number[] | null,
 *     labels: {H, M, L} method strings,  extras: nLayers, embedding, chargeMult [{level, system, charge, multiplicity}],
 *     freezeCodes (Gaussian), qmAtoms / qm2Atoms / subsystemOnly (ORCA: true when only the QM region is printed) }.
 *   Frames are full-system geometries, one per optimisation step; frame.energy is the total only (never a
 *   sub-calculation): Gaussian 'ONIOM' with layer energies 'high-model', 'low-model', 'low-real', 'med-mid', ...;
 *   ORCA 'QM/MM' or 'QM/QM2' (else 'Final') with 'QM', 'MM', 'QM-MM', 'QM SCF' or 'high-model', 'low-model', 'low-real'.
 */
(function (BM) {
  'use strict';

  const BOHR = 0.529177210903;          // Angstrom per Bohr
  const AU2DEBYE = 2.541746;            // Debye per e*Bohr
  const HARTREE2EV = 27.211386245988;
  const EV2NM = 1239.84198;             // nm * eV
  const CM2EV = 1 / 8065.543937;

  // ------------------------------------------------------------------ utilities
  /** Number from a Fortran-style field: D exponents, '******' overflow -> NaN. */
  function num(s) {
    if (typeof s === 'number') return s;
    if (s === undefined || s === null) return NaN;
    s = String(s).trim();
    if (!s || s.indexOf('*') >= 0) return NaN;
    if (/[dD]/.test(s)) s = s.replace(/[dD]/, 'e');
    const v = Number(s);
    return Number.isNaN(v) ? parseFloat(s) : v;
  }
  function int(s) { const v = parseInt(s, 10); return Number.isNaN(v) ? NaN : v; }
  function isNumTok(s) { return /^[-+]?(\d+\.?\d*|\.\d+)([eEdD][-+]?\d+)?$/.test(s); }
  function isIntTok(s) { return /^[-+]?\d+$/.test(s); }
  function tok(line) { const t = line.trim(); return t ? t.split(/\s+/) : []; }
  /** Lines of a text; a line array (BM.bigfile.readLines, BM.parseLines) is returned unchanged. */
  function splitLines(text) { return Array.isArray(text) ? text : String(text).split(/\r\n|\r|\n/); }
  /** Total characters of a line array joined with '\n'. */
  function joinedLength(L) { let n = L.length ? L.length - 1 : 0; for (let k = 0; k < L.length; k++) n += L[k].length; return n; }
  /** Whole text for readers that need it (whole-text regexes, JSON): joins a line array when that is safe
   *  (BM.maxJoinChars, 300 million characters; V8 strings end near 536 million), otherwise a clear error. */
  function asText(x, what) {
    if (!Array.isArray(x)) return String(x);
    const n = joinedLength(x), lim = BM.maxJoinChars || 3e8;
    if (n > lim) {
      throw err('This file is too large for the ' + (what || 'reader') + ' reader (' + Math.round(n / 1e6) + ' million characters; it reads the ' +
        'whole file as one text, at most ' + Math.round(lim / 1e6) + ' million).');
    }
    return x.join('\n');
  }
  /** The first `max` characters of a text or line array (format sniffing). */
  function headOf(x, max) {
    if (!Array.isArray(x)) return String(x).slice(0, max);
    const parts = [];
    let n = 0;
    for (let k = 0; k < x.length && n < max; k++) { parts.push(x[k]); n += x[k].length + 1; }
    return parts.join('\n').slice(0, max);
  }
  function baseName(f) { return f ? String(f).split(/[\\/]/).pop() : ''; }
  function stripExt(f) { return baseName(f).replace(/\.[^.]*$/, ''); }
  /** All decimal numbers in a string, also when glued ("-10.20000-10.19999"). */
  function floatsIn(s) {
    const m = String(s).match(/[-+]?(?:\d+\.\d*|\.\d+|\d+)(?:[eEdD][-+]?\d+)?/g);
    return m ? m.map(num) : [];
  }
  /** Fixed-width fields (Fortran Fw.d); empty fields are skipped. */
  function fixedFields(s, w) {
    const out = [];
    for (let p = 0; p < s.length; p += w) { const f = s.substr(p, w); if (f.trim()) out.push(num(f)); }
    return out;
  }
  function err(msg) { return new Error(msg); }
  function last(a) { return a.length ? a[a.length - 1] : null; }

  // ------------------------------------------------------------------ condensed huge outputs (condense.js)
  // A condensed line array (BM.condense) has a .bmSide store; a line starting with U+0001 stands for a compacted
  // table: PH + 'G' + k -> bmSide.geoms[k] (geometry rows), PH + 'M' + k -> bmSide.modes[k] (normal-mode rows).
  const isPH = s => s !== undefined && s.charCodeAt(0) === 1;
  const phIndex = s => +s.slice(2);
  const DEC10 = [1, 10, 100, 1e3, 1e4, 1e5, 1e6, 1e7, 1e8, 1e9, 1e10, 1e11, 1e12, 1e13, 1e14, 1e15];
  const ZERO_PACK = Object.freeze({ vals: new Int8Array(0), d: 0 });   // a mode printed without rows: zeros
  /**
   * Mode vectors kept as integers with d decimals (huge outputs: 1 byte per value for Gaussian's 2-decimal modes)
   * and decoded on access: modes[k] is a Float64Array of the printed values, or null when the mode's
   * displacements were not loaded. packs[k] = {vals: Int8Array | Int16Array | Int32Array | Float64Array, d} | null.
   */
  function packedModes(packs, N) {
    const out = new Array(packs.length), cache = new Map();
    const decode = k => {
      const p = packs[k];
      if (!p) return null;
      let v = cache.get(k);
      if (v) return v;
      v = new Float64Array(3 * N);
      const src = p.vals, n = Math.min(src.length, v.length);
      if (p.d < 0) for (let i = 0; i < n; i++) v[i] = src[i];
      else { const s = DEC10[p.d]; for (let i = 0; i < n; i++) v[i] = src[i] / s; }
      if (cache.size >= 8) cache.delete(cache.keys().next().value);
      cache.set(k, v);
      return v;
    };
    for (let k = 0; k < packs.length; k++) Object.defineProperty(out, k, { enumerable: true, configurable: true, get: () => decode(k) });
    Object.defineProperty(out, 'packed', { value: true });
    Object.defineProperty(out, 'has', { value: k => !!packs[k] });   // loaded, without decoding it
    return out;
  }

  function newDataset(format, fileName) {
    return {
      format, program: null, title: baseName(fileName) || '',
      meta: { route: null, keywords: null, method: null, basis: null, charge: null, multiplicity: null,
        jobTypes: [], normalTermination: null, solvent: null },
      frames: [], vibrations: null, excitedStates: null, charges: null, spins: null, dipole: null,
      orbitals: null, thermo: null, volumes: null, bonds: null, cell: null, warnings: [],
      // basis: {source: 'fchk' | 'molden', shells: [{atom, center (Bohr), l, pure, exps, coefs}], nbf} and
      // mos: {restricted, alpha: {energies, occupations, symmetries, coeffs (orbital after orbital)}, beta}: BM.orbitals
      basis: null, mos: null,
      // jobs: [{index, label, kinds: subset of opt/freq/sp/td/scan/irc/nmr, frameStart, frameEnd, normalTermination}]
      // nmr: {iso: ppm per atom (NaN = not computed), aniso: number[]|null, method, couplings?: [{i, j, J (Hz)}] (0-based)}
      // optSteps: [{step, energy, opt}] when an output reports optimisation cycles without their geometries
      // ircPath: [{step, energy, dE (kcal/mol), maxGrad, rmsGrad, ts}] (ORCA "IRC PATH SUMMARY")
      jobs: null, nmr: null, optSteps: null, ircPath: null
    };
  }
  /** xyz null: a step whose geometry was not loaded (huge outputs, BM.condense); it becomes an entry of
   *  ds.optSteps with its energies, and no frame (splitPhantoms). */
  function newFrame(Z, xyz) {
    return { Z: Array.from(Z), xyz: xyz === null ? null : xyz instanceof Float64Array ? xyz : Float64Array.from(xyz),
      energy: null, energies: null, label: null, converged: false, opt: null,
      atomNames: null, resNames: null, resIds: null };
  }
  function newThermo() {
    return { zpe: null, eCorr: null, hCorr: null, gCorr: null, eZPE: null, eTotal: null, H: null, G: null,
      temperature: null };
  }
  function sameGeom(f, Z, xyz, tol) {
    tol = tol || 2e-5;
    if (!f || !f.xyz || !xyz || f.Z.length !== Z.length || f.xyz.length !== xyz.length) return false;
    for (let i = 0; i < Z.length; i++) if (f.Z[i] !== Z[i]) return false;
    for (let i = 0; i < xyz.length; i++) if (!(Math.abs(f.xyz[i] - xyz[i]) <= tol)) return false;
    return true;
  }
  function padCell(vecs) {
    if (!vecs || !vecs.length) return null;
    const c = vecs.slice(0, 3).map(v => [v[0], v[1], v[2]]);
    while (c.length < 3) c.push([0, 0, 0]);
    return c;
  }
  function addJobType(meta, t) { if (t && meta.jobTypes.indexOf(t) < 0) meta.jobTypes.push(t); }
  function multName(m) {
    return ['', 'Singlet', 'Doublet', 'Triplet', 'Quartet', 'Quintet', 'Sextet', 'Septet', 'Octet'][m] ||
      (m > 0 ? 'Mult ' + m : '');
  }

  // Energy priority used to pick frame.energy from frame.energies ('MM': Gaussian molecular-mechanics "Energy=").
  const ENERGY_PRIORITY = ['Final', 'CP', 'ONIOM', 'CCSD(T)', 'QCISD(T)', 'CCSD', 'QCISD', 'CCD', 'MP5', 'MP4', 'MP4(SDQ)', 'MP4(DQ)', 'MP3', 'MP2', 'DH', 'SCF', 'MM'];
  const sameZ = (a, b) => { if (!a || a.length !== b.length) return false; for (let i = 0; i < b.length; i++) if (a[i] !== b[i]) return false; return true; };

  /** Collects geometry frames from program outputs, merging repeated identical geometries
   * (Gaussian re-prints the final structure; ORCA's final evaluation; Link1 freq jobs). */
  function FrameSink(ds) {
    this.ds = ds;
    this.pending = null;   // {frame, rank, locked}
    this.dup = false;      // the current frame was re-used for an identical geometry
    this.seen = new Set(); // energy kinds stored since the last geometry block
    this.firstE = new WeakMap();   // frame -> {kind: first value stored for that frame}
    this.priority = ENERGY_PRIORITY;   // kinds that may become frame.energy, best first (QM/MM: only the total)
    this.splitDup = true;              // same geometry re-printed with a different first energy -> separate frame
  }
  FrameSink.prototype.current = function () { return last(this.ds.frames); };
  FrameSink.prototype.geometry = function (Z, xyz, rank) {
    this.seen = new Set();
    const fr = this.ds.frames;
    const cur = last(fr);
    const p = this.pending;
    if (p && p.frame === cur && !p.locked && rank > p.rank) {
      // A preferred orientation of the same structure (Z-matrix < Input < Standard orientation).
      if (p.fresh) {
        const prev = fr.length > 1 ? fr[fr.length - 2] : null;
        if (prev && sameGeom(prev, Z, xyz)) {
          fr.pop(); this.pending = { frame: prev, rank, locked: false, fresh: false }; this.dup = true; return;
        }
        // (a geometry that was not loaded never replaces one that was)
        if (xyz || !cur.xyz) { cur.Z = Array.from(Z); cur.xyz = xyz ? Float64Array.from(xyz) : null; }
        p.rank = rank; return;
      }
      if (sameGeom(cur, Z, xyz)) { p.rank = rank; return; }
    }
    if (cur && sameGeom(cur, Z, xyz)) {
      this.pending = { frame: cur, rank, locked: false, fresh: false }; this.dup = true; return;
    }
    const f = newFrame(Z, xyz);
    if (cur && sameZ(cur.Z, f.Z)) f.Z = cur.Z;   // consecutive frames share one element list (big outputs)
    fr.push(f);
    this.pending = { frame: f, rank, locked: false, fresh: true };
    this.dup = false;
  };
  FrameSink.prototype.energy = function (kind, value) {
    let f = this.current();
    if (!f || !Number.isFinite(value)) return;
    // A re-printed geometry is the same calculation only if its first energy of each kind matches the
    // first one stored for the frame (ONIOM / counterpoise print several "SCF Done" per geometry).
    const first = this.firstE.get(f);
    if (this.splitDup && this.dup && !this.seen.has(kind) && first && first[kind] !== undefined && Math.abs(first[kind] - value) > 2e-6) {
      const g = newFrame(f.Z, f.xyz);   // same geometry, different calculation: separate frame
      this.ds.frames.push(g); f = g; this.dup = false;
      this.pending = { frame: g, rank: this.pending ? this.pending.rank : 0, locked: true, fresh: false };
    }
    this.seen.add(kind);
    let fe = this.firstE.get(f);
    if (!fe) { fe = {}; this.firstE.set(f, fe); }
    if (fe[kind] === undefined) fe[kind] = value;
    if (this.pending && this.pending.frame === f) this.pending.locked = true;
    f.energies = f.energies || {};
    f.energies[kind] = value;
    for (const k of this.priority) if (f.energies[k] !== undefined) { f.energy = f.energies[k]; break; }
  };

  // ------------------------------------------------------------- job list
  const JOB_KINDS = ['opt', 'freq', 'sp', 'td', 'scan', 'irc', 'nmr'];
  function jobKinds(types) {
    const k = JOB_KINDS.filter(t => t !== 'sp' && types.indexOf(t) >= 0);
    return k.length ? k : ['sp'];
  }
  /** ds.jobs bookkeeping (Gaussian Link1 / ORCA $new_job). Call geometry() after every FrameSink.geometry(). */
  function JobTracker(ds) { this.ds = ds; this.list = []; this.cur = null; }
  JobTracker.prototype.open = function (label, kinds) {
    this.close();
    this.cur = { index: this.list.length, label: label || null, kinds: kinds || ['sp'], frameStart: -1, frameEnd: -1,
      normalTermination: null, n0: this.ds.frames.length, merged: -1, geom: false };
    this.list.push(this.cur);
  };
  JobTracker.prototype.geometry = function () {
    const j = this.cur;
    if (!j || j.geom) return;
    j.geom = true;
    if (this.ds.frames.length === j.n0) j.merged = j.n0 - 1;   // first geometry = previous job's last frame
  };
  JobTracker.prototype.close = function () {
    const j = this.cur;
    if (!j) return;
    this.cur = null;
    const fr = this.ds.frames, n = fr.length;
    if (n > j.n0) {
      let s = j.merged >= 0 ? j.merged : j.n0;
      // the merged frame was split off again (same geometry, different calculation): it is not this job's
      if (j.merged >= 0 && sameGeom(fr[j.n0], fr[j.merged].Z, fr[j.merged].xyz)) s = j.n0;
      j.frameStart = s; j.frameEnd = n - 1;
    } else {
      j.frameStart = j.frameEnd = j.merged >= 0 ? j.merged : Math.max(0, n - 1);
    }
    delete j.n0; delete j.merged; delete j.geom;
  };
  JobTracker.prototype.finish = function () { this.close(); return this.list.length ? this.list : null; };

  /** Value of 'R(1,2)' (Angstrom), 'A(1,2,3)' or 'D(1,2,3,4)' (degrees) for a frame; atoms 1-based. */
  function geomValue(coord, xyz) {
    const m = /^([RAD])\(([\d,]+)\)$/.exec(coord || '');
    if (!m || !xyz) return NaN;
    const ix = m[2].split(',').map(s => int(s) - 1);
    if (ix.some(k => k < 0 || 3 * k + 2 >= xyz.length)) return NaN;
    const p = k => [xyz[3 * k], xyz[3 * k + 1], xyz[3 * k + 2]];
    if (m[1] === 'R' && ix.length === 2) return Math.hypot(...vsub(p(ix[0]), p(ix[1])));
    if (m[1] === 'A' && ix.length === 3) {
      const a = vnorm(vsub(p(ix[0]), p(ix[1]))), b = vnorm(vsub(p(ix[2]), p(ix[1])));
      return Math.acos(Math.max(-1, Math.min(1, vdot(a, b)))) / RAD;
    }
    if (m[1] === 'D' && ix.length === 4) {
      const b0 = vsub(p(ix[0]), p(ix[1])), b1 = vnorm(vsub(p(ix[2]), p(ix[1]))), b2 = vsub(p(ix[3]), p(ix[2]));
      const v = vsub(b0, vscale(b1, vdot(b0, b1))), w = vsub(b2, vscale(b1, vdot(b2, b1)));
      return Math.atan2(vdot(vcross(b1, v), w), vdot(v, w)) / RAD;
    }
    return NaN;
  }
  /** Gaussian ModRedundant line -> 'R(1,2)' etc. for scanned coordinates, null for others, undefined if not such a line. */
  function modRedCoord(tk) {
    let k = 0, type = null;
    if (/^[BADLX]$/i.test(tk[0])) { type = tk[0].toUpperCase(); k = 1; }
    const at = [];
    while (k < tk.length && isIntTok(tk[k])) at.push(int(tk[k++]));
    if (!at.length && !type) return undefined;
    if ((tk[k] || '').toUpperCase() !== 'S') return null;
    if (!type) type = ['', '', 'B', 'A', 'D'][at.length] || null;
    const c = { B: 'R', A: 'A', D: 'D' }[type];
    return c && at.length === { R: 2, A: 3, D: 4 }[c] ? c + '(' + at.join(',') + ')' : null;
  }
  /** Lower-triangular matrix printed in column blocks (Gaussian spin-spin couplings). Returns [{i, j, J}] (0-based, i > j). */
  function readLowerTriangle(L, i0) {
    const out = [];
    let cols = null, j = i0;
    for (; j < L.length; j++) {
      const tk = tok(L[j]);
      if (!tk.length) break;
      if (tk.every(isIntTok)) { cols = tk.map(int); continue; }
      if (!cols || !isIntTok(tk[0]) || !tk.slice(1).every(isNumTok)) break;
      const r = int(tk[0]);
      for (let k = 1; k < tk.length && k <= cols.length; k++) if (cols[k - 1] < r) out.push({ i: r - 1, j: cols[k - 1] - 1, J: num(tk[k]) });
    }
    return { pairs: out, end: j };
  }
  function setScan(f, coords, k, value) {
    if (!f || k < 0) return;
    if (!f.scan) f.scan = { coord: coords[0], value: NaN };
    if (k === 0) f.scan.value = value;
    if (coords.length > 1) {
      if (!f.scan.coords) { f.scan.coords = coords.slice(); f.scan.values = coords.map(() => NaN); }
      f.scan.values[k] = value;
    }
  }

  // ------------------------------------------------------- method / basis hints
  const METHOD_RE = /^(?:r|u|ro)?(hf|dft|b3lyp|b3pw91|b3p86|pbe0|pbe1pbe|pbepbe|pbe|revpbe|blyp|bp86|bpw91|b97(?:-?d3?|-?1|-?3c)?|b98|tpssh?|tpsstpss|r2scan(?:-3c)?|scan|m06(?:-?2x|l|hf)?|m05(?:-?2x)?|m08hx|m11l?|mn1[25]l?|wb97x?(?:-?d3?(?:bj)?|-?v)?|wb97m-?v|cam-b3lyp|lc-wpbe|lc-blyp|hse06|hseh1pbe|apfd|x3lyp|o3lyp|b2plypd?3?|b2gp-plyp|mpw2plyp|dsd-?pbep86|pwpb95|mp2|ri-mp2|scs-mp2|mp3|mp4(?:\(sdq\)|sdq|sdtq)?|dlpno-ccsd\(t\)|ccsd(?:\(t\))?|cisd|qcisd(?:\(t\))?|cis|cbs-qb3|g[34](?:mp2)?|pm[367]|am1|pddg|zindo|dftb|gfn[012]?-?xtb|xtb[012]?|hf-3c|pbeh-3c|uff|amber|dreiding|svwn5?|lsda|casscf|nevpt2|caspt2|eom-ccsd)$/i;
  const BASIS_RE = /^(sto-\d+g|\d-\d+\+*g.*|6-31.*|6-311.*|3-21.*|(?:aug-|d-aug-|jun-|may-)?cc-p(?:w?c)?v[dtq56]z.*|(?:ma-|ma-zora-|zora-|dkh-)?def2-?\S+|def-\S+|lanl2(?:dz|tz|mb)\S*|lanl08\S*|sdd\S*|cep-\S+|dgdzvp2?|dgtzvp|tzvpp?|svp|qzvpp?|ugbs|gen|genecp|ano-\S+|pc(?:seg)?-\d|aug-pc-\d|minix|midix|sv\(p\)|x2c-\S+)$/i;

  /** Splits a Gaussian route into keywords, keeping parenthesised option lists together. */
  function routeTokens(route) {
    const out = [];
    let cur = '', depth = 0;
    for (const ch of route) {
      if (ch === '(') depth++;
      else if (ch === ')') depth = Math.max(0, depth - 1);
      if (/\s/.test(ch) && depth === 0) { if (cur) out.push(cur); cur = ''; } else cur += ch;
    }
    if (cur) out.push(cur);
    const res = [];
    for (const t of out) {
      if (res.length && (t[0] === '=' || /=$/.test(res[res.length - 1]))) res[res.length - 1] += t; else res.push(t);
    }
    return res;
  }

  function applyGaussianRoute(route, meta, first) {
    const toks = routeTokens(route.replace(/^\s*#[pPnNtT]?/, '').trim());
    if (first) meta.keywords = toks.slice();
    let method = null, basis = null;
    for (const t of toks) {
      const lo = t.toLowerCase();
      const name = lo.replace(/[=(].*$/, '');
      if (!method && /^oniom/.test(lo)) method = t;
      else if (!method && t.indexOf('/') > 0 && !/^(iop|scrf|geom|guess|pop|density|output|test|ext)\b/.test(lo)) {
        const k = t.indexOf('/');
        method = t.slice(0, k); basis = t.slice(k + 1) || null;
      } else if (!method && METHOD_RE.test(name)) method = t;
      else if (!basis && BASIS_RE.test(lo)) basis = t;
      if (name === 'opt' || name === 'popt') {
        addJobType(meta, 'opt');
        if (/\bts\b|qst[23]/.test(lo)) addJobType(meta, 'ts');
      } else if (name === 'freq' || name === 'frequency') addJobType(meta, 'freq');
      else if (/^(td|tda|tddft|cis|zindo|eomccsd)$/.test(name)) addJobType(meta, 'td');
      else if (name === 'scan') addJobType(meta, 'scan');
      else if (name === 'irc' || name === 'ircmax') addJobType(meta, 'irc');
      else if (name === 'nmr') addJobType(meta, 'nmr');
      else if (name === 'force') addJobType(meta, 'force');
      else if (name === 'stable') addJobType(meta, 'stable');
      else if (name === 'polar') addJobType(meta, 'polar');
      if (name === 'scrf' && !/check/i.test(lo)) {
        const m = /solvent\s*=\s*([^,)\s]+)/i.exec(t);
        if (!meta.solvent) meta.solvent = m ? m[1] : 'Water';
      }
    }
    return { method, basis };
  }
  function finishJobTypes(meta) {
    const main = ['opt', 'freq', 'irc', 'scan', 'td', 'nmr', 'force', 'stable', 'polar'];
    if (!meta.jobTypes.some(t => main.indexOf(t) >= 0)) addJobType(meta, 'sp');
  }

  /** Rows "idx Sym v1 v2 ..." (Gaussian atom tables). Returns {cols, end}. */
  function readIndexedTable(L, start, ncols) {
    const cols = [];
    for (let k = 0; k < ncols; k++) cols.push([]);
    let j = start;
    for (; j < L.length; j++) {
      const t = tok(L[j]);
      if (t.length < 2 + ncols || !isIntTok(t[0]) || isNumTok(t[1])) break;
      for (let k = 0; k < ncols; k++) cols[k].push(num(t[2 + k]));
    }
    return { cols, end: j };
  }

  // --------------------------------------------------------- Gaussian log
  /** A label line of a frequency block: "Frequencies --", hpmodes "Frequencies ---", ROA "Raman1 Fr= 1--" (no blank
   *  before the dashes). A value too wide for its column is printed as asterisks right after the dashes
   *  ("IR Inten    --************", e.g. the huge IR intensities of a charged QM/MM system): it still is a label line. */
  const FREQ_LABEL = /^([A-Za-z].*?)\s*(---|--)(?![-\d.])\s*(.*)$/;
  const STARS = /^\*+$/;
  function parseGaussFreqSection(L, i0) {
    const rec = { hp: false, freqs: [], redMass: [], forceConst: [], ir: [], raman: [], depol: [], sym: [],
      modes: [], nAtoms: 0 };
    let block = null, expectSym = false, depolSeen = false, rowsIn = false, odd = 0, lastA = -1;
    let j = i0 + 1;
    for (; j < L.length; j++) {
      if (block && isPH(L[j]) && L.bmSide) {   // the atom rows of this block, condensed (vecs[k] null: not loaded)
        const run = L.bmSide.modes[phIndex(L[j])];
        if (!rec.packed) rec.packed = [];
        block.forEach((m, k) => { rec.packed[m] = run.vecs[k] || null; });
        rec.nAtoms = Math.max(rec.nAtoms, run.nAtoms);
        rowsIn = true;
        continue;
      }
      const t = L[j].trim();
      if (!t) { if (rec.freqs.length) break; continue; }
      if (!block && /^(activities|incident light|and normal coordinates|depolarization|reduced masses)/i.test(t)) continue;
      const tk = t.split(/\s+/);
      if (tk.every(isIntTok)) {
        block = tk.map(x => int(x) - 1);
        for (const m of block) if (!rec.modes[m]) rec.modes[m] = [];
        expectSym = true; depolSeen = false; rowsIn = false; odd = 0; lastA = -1;
        continue;
      }
      const lm = FREQ_LABEL.exec(t);
      if (lm && block) {
        expectSym = false;
        if (lm[2] === '---') rec.hp = true;
        const label = lm[1].toLowerCase();
        const vals = tok(lm[3] || '').map(num);
        let arr = null;
        if (label.startsWith('frequencies')) arr = rec.freqs;
        else if (label.startsWith('red. masses') || label.startsWith('reduced masses')) arr = rec.redMass;
        else if (label.startsWith('frc consts') || label.startsWith('force constants')) arr = rec.forceConst;
        else if (label.startsWith('ir inten')) arr = rec.ir;
        else if (label.startsWith('raman activ')) arr = rec.raman;
        else if (label.startsWith('depolar') && !depolSeen && label.indexOf('(u)') < 0) { arr = rec.depol; depolSeen = true; }
        if (arr) block.forEach((m, k) => { arr[m] = k < vals.length ? vals[k] : NaN; });
        continue;
      }
      if (expectSym && block && tk.length === block.length) {
        block.forEach((m, k) => { rec.sym[m] = tk[k]; });
        expectSym = false;
        continue;
      }
      if (/^Atom\s+AN\b/.test(t) || /^Coord\s+Atom\s+Element/.test(t)) continue;
      // an atom number too wide for its column (asterisks) is the one after the previous row's
      if (block && tk.length >= 3 && (isIntTok(tk[0]) || (STARS.test(tk[0]) && lastA >= 0)) && isIntTok(tk[1])) {
        if (rec.hp && isIntTok(tk[0]) && isIntTok(tk[2]) && tk.length === 3 + block.length) {
          const c = int(tk[0]) - 1, a = int(tk[1]) - 1;
          if (c >= 0 && c < 3 && a >= 0) {
            block.forEach((m, k) => { rec.modes[m][3 * a + c] = num(tk[3 + k]); });
            rec.nAtoms = Math.max(rec.nAtoms, a + 1);
            rowsIn = true;
            continue;
          }
        }
        if (tk.length === 2 + 3 * block.length) {
          const a = isIntTok(tk[0]) ? int(tk[0]) - 1 : lastA + 1;
          block.forEach((m, k) => { for (let c = 0; c < 3; c++) rec.modes[m][3 * a + c] = num(tk[2 + 3 * k + c]); });
          rec.nAtoms = Math.max(rec.nAtoms, a + 1);
          rowsIn = true; lastA = a;
          continue;
        }
      }
      if (/^Harmonic frequencies/.test(t)) break;   // the next section (hpmodes prints both)
      if (rowsIn && /^-{5,}$/.test(t)) break;   // the box of " - Thermochemistry -" right after the last block
      // a line this reader does not know among the labels of a block, before its atom rows: skipped (a few at most)
      if (block && !rowsIn && ++odd <= 8) continue;
      if (rec.freqs.length) rec.stopped = { line: j, text: /^\d/.test(t) ? null : t.slice(0, 70) };
      break;
    }
    return { rec, end: j };
  }

  function buildVibrations(rec, frameIndex, ds) {
    const n = rec.freqs.length;
    if (!n) return null;
    const frame = ds.frames[frameIndex];
    const N = frame ? frame.Z.length : rec.nAtoms;
    if (rec.nAtoms && rec.nAtoms !== N) {
      ds.warnings.push('Normal modes list ' + rec.nAtoms + ' atoms but the geometry has ' + N + '.');
    }
    const fill = a => {
      if (!a || !a.length) return null;
      const out = new Array(n);
      for (let i = 0; i < n; i++) out[i] = a[i] === undefined ? NaN : a[i];
      return out;
    };
    let modes = [];
    // (packed[m] null: not loaded; undefined: no rows were printed, which the full reading turns into zeros)
    if (rec.packed) modes = packedModes(Array.from({ length: n }, (_, m) => (rec.packed[m] === undefined ? ZERO_PACK : rec.packed[m])), N);
    else {
      for (let m = 0; m < n; m++) {
        const v = new Float64Array(3 * N);
        const src = rec.modes[m] || [];
        for (let k = 0; k < 3 * N; k++) { const x = src[k]; v[k] = x === undefined || Number.isNaN(x) ? 0 : x; }
        modes.push(v);
      }
    }
    const out = { frameIndex: Math.max(0, frameIndex), freqs: fill(rec.freqs), ir: fill(rec.ir), raman: fill(rec.raman),
      depol: fill(rec.depol), redMass: fill(rec.redMass), forceConst: fill(rec.forceConst),
      sym: rec.sym.length ? fill(rec.sym).map(s => (typeof s === 'string' ? s : '')) : null, modes };
    if (rec.packed) out.modesLoaded = modes.reduce((s, m, k) => s + (rec.packed[k] === null ? 0 : 1), 0);
    return out;
  }

  function readGaussOrientation(L, i) {
    if (isPH(L[i + 1]) && L.bmSide) {   // condensed table (xyz null: a step whose geometry was not loaded)
      const g = L.bmSide.geoms[phIndex(L[i + 1])];
      return { Z: g.Z, xyz: g.xyz, tv: g.tv, end: i + 1, complete: g.complete };
    }
    let j = i + 1, dashes = 0;
    while (j < L.length && dashes < 2) { if (/^\s*-{10,}\s*$/.test(L[j])) dashes++; j++; }
    const Z = [], xyz = [], tv = [];
    for (; j < L.length; j++) {
      const t = tok(L[j]);
      if (t.length < 3 || !isIntTok(t[0]) || !isIntTok(t[1])) break;
      // 3F12.6 fields run together when |coordinate| >= 1000 ("-1234.567890-1000.100000")
      const f = floatsIn(L[j]);
      if (f.length < 5) break;
      const an = int(t[1]);
      const x = f[f.length - 3], y = f[f.length - 2], z = f[f.length - 1];
      if (an === -2) tv.push([x, y, z]);
      else if (an >= 0) { Z.push(an); xyz.push(x, y, z); }
    }
    // A table without its closing dashes is a file cut off while being written.
    const complete = j < L.length && /^\s*-{10,}\s*$/.test(L[j]);
    return { Z, xyz, tv, end: j, complete };
  }

  // ------------------------------------------------------------------ ONIOM (Gaussian)
  /** "C-CT--0.1234(PDBName=CA,ResName=ALA,ResNum=1)" -> {el, type, charge, props}; null if not an atom spec. */
  function oniomAtomSpec(s) {
    const m = /^([A-Za-z][a-z]?)(\d*)((?:-[^\s(]*)?)(?:\((.*)\))?$/.exec(s);
    if (!m || BM.el.toZ(m[1]) < 0) return null;
    let type = null, charge = null;
    const rest = m[3] ? m[3].slice(1) : '';
    if (m[3]) {
      const k = rest.indexOf('-');
      type = (k < 0 ? rest : rest.slice(0, k)) || null;
      if (k >= 0) { const q = num(rest.slice(k + 1)); if (Number.isFinite(q)) charge = q; }
    }
    const props = {};
    if (m[4]) for (const kv of m[4].split(',')) { const e = kv.indexOf('='); if (e > 0) props[kv.slice(0, e).trim()] = kv.slice(e + 1).trim(); }
    return { el: m[1], Z: BM.el.toZ(m[1]), type, charge, props };
  }
  /** The fields of an atom row, split on spaces and on commas outside parentheses:
   *  "N-N3-0.1592(PDBName=N,ResName=SER,ResNum=1),0,1.2,0.3,0.0,L" -> [spec, '0', '1.2', '0.3', '0.0', 'L']. */
  function rowTokens(t) {
    const out = [];
    let cur = '', d = 0;
    for (const ch of String(t)) {
      if (ch === '(') d++; else if (ch === ')') d = Math.max(0, d - 1);
      if (d === 0 && (ch === ',' || ch === ' ' || ch === '\t')) { if (cur) out.push(cur); cur = ''; continue; }
      cur += ch;
    }
    if (cur) out.push(cur);
    return out;
  }
  /** One atom row of an ONIOM input echo: spec [freeze] x y z [H|M|L [link-spec bondedTo [factors...]]]. The same
   *  fields separated by commas are the rows of a structure taken from a checkpoint file ("C,0,x,y,z"). */
  function oniomAtomRow(t) {
    const tk = rowTokens(t);
    if (!tk.length) return null;
    const a = oniomAtomSpec(tk[0]);
    if (!a) return null;
    let k = 1, freeze = null;
    if (tk.length >= 5 && isIntTok(tk[1]) && isNumTok(tk[2]) && isNumTok(tk[3]) && isNumTok(tk[4])) { freeze = int(tk[1]); k = 2; }
    if (!(isNumTok(tk[k]) && isNumTok(tk[k + 1]) && isNumTok(tk[k + 2]))) return null;
    a.freeze = freeze;
    a.xyz = [num(tk[k]), num(tk[k + 1]), num(tk[k + 2])];
    k += 3;
    a.layer = k < tk.length && /^[HMLhml]$/.test(tk[k]) ? tk[k++].toUpperCase() : null;
    a.link = null;
    if (a.layer && k < tk.length) {
      const l = oniomAtomSpec(tk[k]);
      if (l) {
        a.link = { el: l.el, type: l.type, charge: l.charge, bondedTo: isIntTok(tk[k + 1] || '') ? int(tk[k + 1]) - 1 : -1, factors: [] };
        for (let q = k + 2; q < tk.length && isNumTok(tk[q]); q++) a.link.factors.push(num(tk[q]));
      }
    }
    return a;
  }
  const ONIOM_LEVEL = { high: 'high', hi: 'high', med: 'med', medium: 'med', mid: 'med', middle: 'med', low: 'low', lo: 'low' };
  const ONIOM_SYSTEM = { model: 'model', small: 'model', mid: 'mid', middle: 'mid', med: 'mid', medium: 'mid', inter: 'mid', intermediate: 'mid', int: 'mid', real: 'real', full: 'real' };
  const oniomKey = (lev, sys) => (ONIOM_LEVEL[String(lev).toLowerCase()] || String(lev).toLowerCase()) + '-' +
    (ONIOM_SYSTEM[String(sys).toLowerCase()] || String(sys).toLowerCase());
  /** The "Symbolic Z-matrix:" echo: per-layer charge / multiplicity and the atom rows (Cartesian ONIOM input). */
  function readOniomEcho(L, j) {
    const cm = [], atoms = [];
    for (; j < L.length; j++) {
      const t = L[j].trim();
      if (!t) { if (atoms.length) break; continue; }
      const m = /^Charge\s*=\s*(-?\d+)\s+Multiplicity\s*=\s*(\d+)(?:\s+for\s+(\w+)\s+level\s+calculation\s+on\s+(\w+)\s+system)?/i.exec(t);
      if (m) {
        cm.push({ charge: int(m[1]), multiplicity: int(m[2]), level: m[3] ? (ONIOM_LEVEL[m[3].toLowerCase()] || m[3].toLowerCase()) : null,
          system: m[4] ? (ONIOM_SYSTEM[m[4].toLowerCase()] || m[4].toLowerCase()) : null });
        continue;
      }
      const a = oniomAtomRow(t);
      if (!a) break;
      atoms.push(a);
    }
    return { cm, atoms };
  }
  /** "Structure from the checkpoint file:" (Geom=Check / AllCheck): a title, the charge / multiplicity lines and rows
   *  "C,0,x,y,z" (element, freeze code, coordinates). Gaussian prints no ONIOM layers here; rows that carry them anyway
   *  are read like the input echo. */
  function readChkStructure(L, j) {
    const cm = [], atoms = [];
    const end = Math.min(L.length, j + 40);
    let first = -1;
    for (; j < end; j++) {          // up to the first row: title, charges, "Redundant internal coordinates found in file."
      const t = L[j].trim();
      const m = /^Charge\s*=\s*(-?\d+)\s+Multiplicity\s*=\s*(\d+)(?:\s+for\s+(\w+)\s+level\s+calculation\s+on\s+(\w+)\s+system)?/i.exec(t);
      if (m) {
        cm.push({ charge: int(m[1]), multiplicity: int(m[2]), level: m[3] ? (ONIOM_LEVEL[m[3].toLowerCase()] || m[3].toLowerCase()) : null,
          system: m[4] ? (ONIOM_SYSTEM[m[4].toLowerCase()] || m[4].toLowerCase()) : null });
        continue;
      }
      if (t.indexOf(',') > 0 && oniomAtomRow(t)) { first = j; break; }
    }
    if (first < 0) return { cm, atoms };
    for (j = first; j < L.length; j++) {
      const a = oniomAtomRow(L[j].trim());
      if (!a) break;
      atoms.push(a);
    }
    return { cm, atoms };
  }
  /**
   * BM.qmRegionFromLinks(Z, xyz, links, nQM) -> sorted QM atom indices | null
   * The QM region of a QM/MM partition that the output does not list atom by atom: the atoms covalently bonded (the
   * viewer's own rule, radii + 0.4 Å) to the QM ends of the link atoms (links: [{atom: replaced MM atom, bondedTo: QM
   * atom}]), never crossing a link bond or entering a replaced atom. A hydrogen ends the walk (one bond: short contacts
   * to other atoms are no path). Accepted only when it has exactly nQM atoms; QM atoms bonded to no link atom (waters,
   * ions, a ligand on its own) cannot be found this way, and the result is null.
   */
  function qmRegionFromLinks(Z, X, links, nQM) {
    const N = Z ? Z.length : 0;
    if (!N || !X || X.length < 3 * N || !links || !links.length || !(nQM > 0) || nQM >= N) return null;
    const host = new Uint8Array(N), blocked = new Set();
    for (const l of links) {
      if (!(l.atom >= 0 && l.atom < N && l.bondedTo >= 0 && l.bondedTo < N)) return null;
      host[l.atom] = 1;
      blocked.add(l.atom * N + l.bondedTo); blocked.add(l.bondedTo * N + l.atom);
    }
    const cov = z => { const r = BM.el.covalent[z]; return r > 0 ? r : 1.5; };
    let maxR = 0;
    for (let a = 0; a < N; a++) if (Z[a] > 0) maxR = Math.max(maxR, cov(Z[a]));
    const cell = 2 * maxR + 0.4, grid = new Map(), key = (i, j, k) => i + ',' + j + ',' + k;
    const ci = a => [Math.floor(X[3 * a] / cell), Math.floor(X[3 * a + 1] / cell), Math.floor(X[3 * a + 2] / cell)];
    for (let a = 0; a < N; a++) {
      if (!(Z[a] > 0)) continue;
      const c = ci(a), k = key(c[0], c[1], c[2]);
      const l = grid.get(k); if (l) l.push(a); else grid.set(k, [a]);
    }
    const inQM = new Uint8Array(N), stack = [];
    for (const l of links) if (!host[l.bondedTo] && !inQM[l.bondedTo]) { inQM[l.bondedTo] = 1; stack.push(l.bondedTo); }
    let count = stack.length;
    while (stack.length) {
      const a = stack.pop();
      if (Z[a] === 1) continue;
      const c = ci(a), ra = cov(Z[a]);
      for (let di = -1; di <= 1; di++) for (let dj = -1; dj <= 1; dj++) for (let dk = -1; dk <= 1; dk++) {
        const l = grid.get(key(c[0] + di, c[1] + dj, c[2] + dk));
        if (!l) continue;
        for (const b of l) {
          if (inQM[b] || host[b] || blocked.has(a * N + b)) continue;
          const dx = X[3 * b] - X[3 * a], dy = X[3 * b + 1] - X[3 * a + 1], dz = X[3 * b + 2] - X[3 * a + 2], d2 = dx * dx + dy * dy + dz * dz;
          const lim = ra + cov(Z[b]) + 0.4;
          if (d2 < lim * lim && d2 > 0.16) {
            inQM[b] = 1; stack.push(b);
            if (++count > nQM) return null;   // the region leaks past the count: some boundary bond has no link atom
          }
        }
      }
    }
    if (count !== nQM) return null;
    const out = [];
    for (let a = 0; a < N; a++) if (inQM[a]) out.push(a);
    return out;
  }
  /**
   * BM.oniomLayersFromPopulation(Z, pop, opts) -> {atomLayer, hosts: [{atom, element, system}], qm, mid, how} | null
   * The layers of a Gaussian ONIOM job that does not list them (Geom=Check / AllCheck), from the population analysis
   * of its QM sub-calculations. Gaussian prints a Mulliken charge for every atom of the real system there:
   *  - the atoms of the sub-calculation's system (its QM atoms, and the atoms that link atoms replace, printed with
   *    the link atom's element, usually H) have charges;
   *  - the other atoms have 0 or, with electronic embedding, their fixed MM charges, the same in every table, while the
   *    charges of the QM atoms change from one geometry to the next.
   * pop: {'high-model': P, 'med-mid': P} with P = {el: element symbol per atom, q: the first table's charges,
   * changed: per atom 1 when a later table differs, tables: count, nActive: Gaussian's atom count of that system}.
   * opts: {cutHosts: atoms replaced by link atoms ("ONIOM: Cut between"), three: 3-layer job, nModel: model-system
   * atoms when NActive is missing (NQM, with Gaussian's own MM)}. The atoms found must be exactly as many as Gaussian
   * counts, or the result is null (nothing is guessed).
   */
  function oniomLayersFromPopulation(Z, pop, opts) {
    opts = opts || {};
    const N = Z ? Z.length : 0;
    if (!N || !pop) return null;
    const cut = new Set((opts.cutHosts || []).filter(a => a >= 0 && a < N));
    const system = (p, alt) => {
      if (!p || !p.q || p.q.length !== N || !p.el || p.el.length !== N) return null;
      const target = p.nActive > 0 ? p.nActive : alt > 0 ? alt : 0;
      if (!target || target > N) return null;
      const elDiff = new Uint8Array(N);
      for (let a = 0; a < N; a++) {
        const z = BM.el.toZ(p.el[a]);
        if (z < 0) return null;
        if (z !== Z[a]) elDiff[a] = 1;
      }
      let act = new Uint8Array(N), n = 0, how = 'charges';
      for (let a = 0; a < N; a++) if (p.q[a] !== 0) { act[a] = 1; n++; }
      if (n !== target) {     // electronic embedding: MM charges on the other atoms, the same in every table
        if (!(p.changed && p.tables > 1)) return null;
        act = new Uint8Array(N); n = 0; how = 'changed';
        for (let a = 0; a < N; a++) if (p.changed[a] || elDiff[a]) { act[a] = 1; n++; }
        if (n !== target) return null;
      }
      return { act, elDiff, how, tables: p.tables || 1, el: p.el };
    };
    const m = system(pop['high-model'], opts.three ? 0 : opts.nModel);
    if (!m) return null;
    const hostM = new Uint8Array(N), layer = new Array(N).fill('L');
    let qm = 0, mid = 0;
    for (let a = 0; a < N; a++) {
      if (!m.act[a]) continue;
      if (m.elDiff[a] || cut.has(a)) hostM[a] = 1; else { layer[a] = 'H'; qm++; }
    }
    if (!qm) return null;
    const hosts = [];
    for (let a = 0; a < N; a++) if (hostM[a]) hosts.push({ atom: a, element: m.el[a], system: 'model' });
    let how = m.how;
    if (opts.three) {
      // the medium layer: the mid system's own atoms that are not high; replaced atoms there are low-layer atoms
      const i = system(pop['med-mid'], 0);
      if (!i) return null;
      for (let a = 0; a < N; a++) {
        if (layer[a] === 'H') { if (!i.act[a]) return null; continue; }
        if (!i.act[a]) continue;
        if (i.elDiff[a] || (cut.has(a) && !hostM[a])) { if (!hostM[a]) hosts.push({ atom: a, element: i.el[a], system: 'mid' }); }
        else { layer[a] = 'M'; mid++; }
      }
      if (!mid) return null;
      if (i.how !== how) how = 'changed';
    } else {
      for (const a of cut) if (!m.act[a]) return null;   // every link atom belongs to the model system
    }
    return { atomLayer: layer, hosts, qm, mid, how, tables: m.tables };
  }
  /** Method labels {H, M, L} from "ONIOM(B3LYP/6-31G(d):Amber=SoftFirst)=EmbedCharge". */
  function oniomLabels(route) {
    const i = route.search(/oniom\s*\(/i);
    const out = { H: null, M: null, L: null };
    if (i < 0) return out;
    let k = route.indexOf('(', i) + 1, depth = 1, cur = '';
    const parts = [];
    for (; k < route.length && depth > 0; k++) {
      const ch = route[k];
      if (ch === '(') depth++;
      else if (ch === ')') { depth--; if (!depth) break; }
      if (ch === ':' && depth === 1) { parts.push(cur); cur = ''; } else cur += ch;
    }
    parts.push(cur);
    const tidy = p => {
      let d = 0, o = '';
      for (const ch of p.trim()) {           // drop "=options", upper-case outside parentheses: 6-31G(d)
        if (ch === '(') d++; else if (ch === ')') d--;
        if (ch === '=' && d === 0) break;
        o += d > 0 || ch === ')' ? ch : ch.toUpperCase();
      }
      return o || null;
    };
    const lab = parts.map(tidy);
    out.H = lab[0] || null;
    if (lab.length >= 3) { out.M = lab[1]; out.L = lab[2]; } else out.L = lab[1] || null;
    return out;
  }
  /**
   * The layers of a Gaussian ONIOM job (output or input). Per-atom layers come from, in this order:
   *  - the molecule echo ("Symbolic Z-matrix:", or the rows of an input file), source 'input';
   *  - a structure taken from the checkpoint file whose rows carry layers, source 'checkpoint';
   *  - the population analysis of the high-level (and medium-level) sub-calculations, source 'population': with
   *    Geom=Check / AllCheck Gaussian lists no layers, but its Mulliken charges there cover exactly the atoms of the
   *    model (mid) system (oniomLayersFromPopulation);
   *  - rebuilt (source 'rebuilt'): with Geom=Check / AllCheck Gaussian lists no layers, but prints the link atoms
   *    ("ONIOM: Cut between") and NQM, the atoms of the model system other than link atoms (high-layer atoms plus the
   *    atoms the link atoms replace, with Gaussian's own MM). The atoms bonded to the QM ends of the link atoms make the
   *    high layer when their number matches (2-layer QM:MM only);
   *  - otherwise atomLayer stays null (unknown).
   * more: {chk: readChkStructure result, counts: {nAtoms, nqm}, frame: a full-system frame, pop: population analyses}.
   */
  function buildOniomLayers(route, echo, nReal, cuts, ee, more) {
    more = more || {};
    let atoms = echo && echo.atoms && echo.atoms.length === nReal ? echo.atoms : null, source = atoms ? 'input' : null;
    const chk = more.chk && more.chk.atoms.length === nReal ? more.chk.atoms : null;
    if (!atoms && chk && chk.some(a => a.layer)) { atoms = chk; source = 'checkpoint'; }
    const typed = atoms || chk;      // freeze codes, MM types and charges also from a structure without layers
    const labels = oniomLabels(route);
    const anyFreeze = typed && typed.some(a => a.freeze !== null && a.freeze !== 0);
    const frozen = anyFreeze ? typed.map(a => a.freeze === -1) : null;
    let linkAtoms = [];
    if (atoms) atoms.forEach((a, k) => { if (a.link) linkAtoms.push({ atom: k, bondedTo: a.link.bondedTo, element: a.link.el, type: a.link.type || null, factors: a.link.factors.length ? a.link.factors : null }); });
    if (!linkAtoms.length && cuts.length) linkAtoms = cuts.map(c => ({ atom: c.atom, bondedTo: c.bondedTo, element: c.element, type: null, factors: c.factor ? [c.factor] : null }));
    let atomLayer = atoms ? atoms.map(a => a.layer || 'H') : null, rebuilt = null, population = null;
    const c = more.counts, fr = more.frame;
    if (!atomLayer && more.pop && fr && fr.Z.length === nReal) {
      const r = oniomLayersFromPopulation(fr.Z, more.pop, { cutHosts: cuts.map(x => x.atom), three: !!labels.M,
        nModel: c && c.nAtoms === nReal && c.nqm > 0 && c.nqm < nReal ? c.nqm : 0 });
      if (r) {
        atomLayer = r.atomLayer; source = 'population';
        population = { how: r.how, qm: r.qm, mid: r.mid, links: r.hosts.length, tables: r.tables };
        if (!linkAtoms.length && r.hosts.length) {
          // no "ONIOM: Cut between" lines: each replaced atom is bonded to the nearest atom of the layer above it
          const X = fr.xyz, lay = r.atomLayer;
          linkAtoms = r.hosts.map(h => {
            let best = -1, bd = Infinity;
            for (let b = 0; b < nReal; b++) {
              if (b === h.atom || !(h.system === 'model' ? lay[b] === 'H' : lay[b] !== 'L')) continue;
              const dx = X[3 * b] - X[3 * h.atom], dy = X[3 * b + 1] - X[3 * h.atom + 1], dz = X[3 * b + 2] - X[3 * h.atom + 2], d2 = dx * dx + dy * dy + dz * dz;
              if (d2 < bd) { bd = d2; best = b; }
            }
            return { atom: h.atom, bondedTo: best, element: h.element, type: null, factors: null };
          });
        }
      }
    }
    if (!atomLayer && c && fr && fr.Z.length === nReal && !labels.M && /^(amber|uff|dreiding)/i.test(labels.L || '') &&
        c.nAtoms === nReal && c.nqm > 0 && c.nqm < nReal && linkAtoms.length) {
      const hosts = new Set(linkAtoms.map(l => l.atom));
      const qm = qmRegionFromLinks(fr.Z, fr.xyz, linkAtoms, c.nqm - hosts.size);
      if (qm) {
        const s = new Set(qm);
        atomLayer = Array.from({ length: nReal }, (_, a) => (s.has(a) ? 'H' : 'L'));
        source = 'rebuilt'; rebuilt = { qm: qm.length, links: linkAtoms.length };
      }
    }
    const layers = {
      scheme: 'oniom',
      atomLayer,
      frozen, active: frozen ? frozen.map(f => !f) : null,
      linkAtoms,
      mmTypes: typed && typed.some(a => a.type) ? typed.map(a => a.type || null) : null,
      mmCharges: typed && typed.some(a => a.charge !== null) ? typed.map(a => (a.charge === null ? NaN : a.charge)) : null,
      labels,
      nLayers: labels.M ? 3 : 2,
      embedding: ee ? 'electronic' : 'mechanical',
      chargeMult: echo && echo.cm ? echo.cm.filter(x => x.level) : more.chk ? more.chk.cm.filter(x => x.level) : [],
      freezeCodes: anyFreeze ? typed.map(a => (a.freeze === null ? 0 : a.freeze)) : null,
      source, rebuilt, population
    };
    return layers;
  }

  function parseGaussianLog(text, fileName) {
    const ds = newDataset('gaussian-log', fileName);
    const meta = ds.meta;
    const L = splitLines(text);
    const n = L.length;
    const sink = new FrameSink(ds);
    const charges = {}, spins = {};
    const freqSets = [];
    const jobStatus = [];
    let job = 0, routes = 0, routeEnd = -1, lastDash = -1, prevDash = -1, titleSet = false;
    let states = null, inState = false, temperature = null;
    let eig = null, eigLast = null, prevEig = false, nboSpin = false;
    let scfMethod = null, stdBasis = null, solventOut = null, routeMB = null, cellVecs = null;
    // Excited-state optimisations / scans / IRCs: the TD total energy is the frame energy of that job.
    let tdGeomJob = false;
    let negEig = null;   // negative eigenvalues of the Hessian in this optimisation step (1 for a TS search)
    let ccKind = 'CCSD', lastCorr = null;   // coupled-cluster / QCI: the kind of this job, the last E(Corr) printed
    const tdFrames = new Set();
    const jobsT = new JobTracker(ds);
    let scanCoords = [], paramTable = null, rotLen = null, rotVel = null, nmr = null;
    const ircJobs = new Set();
    // ONIOM: frames are real-system geometries, one per optimisation step; model-system prints (QM region + link
    // atoms) are skipped and only "ONIOM: extrapolated energy" may become frame.energy (sub-calculation SCF energies
    // never stand in for it). Layer energies (gridpoints) go to frame.energies as 'high-model', 'low-real', ...
    let oniomJob = false, oniomOpt = false, oniomRoute = null, oniomEcho = null, nReal = 0, stepOpen = false, oniomEDone = false, oniomEE = false;
    let oniomSub = null;   // layer key of the sub-calculation being printed ("ONIOM: generating point 2 -- high level on model system.")
    let oniomChk = null, oniomCounts = null;   // structure taken from the checkpoint (Geom=Check); "NAtoms= NQM=" counts
    // population analyses of the QM sub-calculations (see oniomLayersFromPopulation), by key ('high-model', 'med-mid')
    const oniomPop = {};
    const popOf = k => oniomPop[k] || (oniomPop[k] = { nActive: null, el: null, q: null, changed: null, tables: 0 });
    const oniomFrames = new Map(), oniomCuts = [], oniomMM = new Map();   // ONIOM frame -> rank of its orientation; frame -> MM "Energy=" values
    const E = (kind, v) => { if (!oniomJob) sink.energy(kind, v); };   // sub-calculation energies: not in ONIOM jobs
    const closeStep = () => { stepOpen = false; oniomEDone = false; oniomSub = null; };
    const oniomGeometry = (g, rank) => {
      const cur = sink.current();
      if (oniomOpt && stepOpen && cur && oniomFrames.has(cur)) {
        // another real-system print within the same step (micro-iterations, a second orientation): it updates the
        // step's coordinates until the step's ONIOM energy is known (Standard orientation preferred), then is ignored
        if (!oniomEDone && sameZ(cur.Z, g.Z) && rank >= (oniomFrames.get(cur) || 0)) {
          if (g.xyz) cur.xyz = Float64Array.from(g.xyz);   // (a geometry that was not loaded changes nothing)
          oniomFrames.set(cur, rank);
        }
        return;
      }
      const n0 = ds.frames.length;
      sink.geometry(g.Z, g.xyz, rank); jobsT.geometry();
      if (ds.frames.length > n0) { oniomEDone = false; oniomSub = null; }
      oniomFrames.set(sink.current(), Math.max(rank, oniomFrames.get(sink.current()) || 0));
      stepOpen = true;
    };
    // Every orientation block is passed on: within one step "Standard" (rank 2) replaces "Input" (1) or
    // "Z-Matrix" (0) orientation; steps printed without a Standard orientation (NoSymm Link1 job,
    // symmetry turned off) keep their Input orientation instead of being dropped.

    for (let i = 0; i < n; i++) {
      const line = L[i];
      // Rows starting with a digit (MO coefficients, density and distance matrices, SCF tables, population
      // tables) are handled by no rule here except excited-state transitions: skip them without allocating.
      let p0 = 0, c0 = line.charCodeAt(0);
      while (c0 === 32) c0 = line.charCodeAt(++p0);
      if (c0 >= 48 && c0 <= 57 && !inState) { prevEig = false; continue; }
      const t = line.trimStart();
      if (!t) { prevEig = false; inState = false; continue; }
      const c = t[0];
      if ((c === 'A' || c === 'B') && /^(Alpha|Beta)\s+(occ|virt)\.\s+eigenvalues\s*--/.test(t)) {
        const m = /^(Alpha|Beta)\s+(occ|virt)\.\s+eigenvalues\s*--(.*)$/.exec(t);
        if (!prevEig) { eig = { aOcc: [], aVirt: [], bOcc: null, bVirt: null }; eigLast = eig; }
        const vals = m[3].indexOf('*') >= 0 ? fixedFields(m[3], 10) : floatsIn(m[3]);
        if (m[1] === 'Alpha') (m[2] === 'occ' ? eig.aOcc : eig.aVirt).push(...vals);
        else {
          if (!eig.bOcc) { eig.bOcc = []; eig.bVirt = []; }
          (m[2] === 'occ' ? eig.bOcc : eig.bVirt).push(...vals);
        }
        prevEig = true;
        continue;
      }
      prevEig = false;
      if (inState && /^\d+[AB]?\s*(->|<-)\s*\d+[AB]?\s+\S+/.test(t) && states && states.length) {
        const m = /^(\d+[AB]?)\s*(->|<-)\s*(\d+[AB]?)\s+(\S+)/.exec(t);
        const tr = { from: m[1], to: m[3], coeff: num(m[4]) };
        if (m[2] === '<-') tr.deexcitation = true;
        last(states).transitions.push(tr);
        continue;
      }
      inState = false;
      if (c === '%' && routes === 0) {   // Link 0 commands of the first job: %chk, %oldchk (where Geom=Check reads from)
        const m = /^%(old)?chk\s*=\s*(\S+)/i.exec(t);
        if (m) { if (m[1]) meta.oldchk = m[2]; else meta.chk = m[2]; }
        continue;
      }
      if (c === '!') {
        // "! Initial / Optimized / Non-Optimized Parameters !" tables: scanned coordinates and their values
        const h = /^!\s*(Initial|Optimized|Non-Optimized)\s+Parameters\s*!/.exec(t);
        if (h) { paramTable = h[1]; continue; }
        const r = paramTable && /^!\s*\S+\s+([RAD]\([\d,]+\))\s+(-?\d+\.?\d*)\s+(.*?)\s*!\s*$/.exec(t);
        if (r) {
          if (paramTable === 'Initial') { if (/\bScan\b/.test(r[3]) && scanCoords.indexOf(r[1]) < 0) scanCoords.push(r[1]); }
          else setScan(sink.current(), scanCoords, scanCoords.indexOf(r[1]), num(r[2]));
        }
        continue;
      }
      if (c === 'T' && t.startsWith('The following ModRedundant input section has been read')) {
        for (let j = i + 1; j < n && j < i + 200; j++) {
          const tk = tok(L[j]);
          const sc = tk.length ? modRedCoord(tk) : undefined;
          if (sc === undefined) break;
          if (sc && scanCoords.indexOf(sc) < 0) scanCoords.push(sc);
          i = j;
        }
        continue;
      }
      if ((c === 'S' || c === 'C' || c === 'I' || c === 'M' || c === 'G') && t.indexOf('Magnetic shielding tensor (ppm)') >= 0) {
        const method = /^(.*?)\s*Magnetic shielding tensor/.exec(t)[1].replace(/^SCF\s+/, '').trim() || 'GIAO';
        const iso = [], aniso = [];
        let j = i + 1;
        for (; j < n; j++) {
          const s = L[j];
          const m = /^\s*(\d+)\s+\S+\s+Isotropic\s*=\s*(\S+)\s+Anisotropy\s*=\s*(\S+)/.exec(s);
          if (m) { iso[int(m[1]) - 1] = num(m[2]); aniso[int(m[1]) - 1] = num(m[3]); continue; }
          if (/^\s*(XX=|XY=|XZ=|Eigenvalues:)/.test(s)) continue;
          break;
        }
        if (iso.length) nmr = { iso: Array.from(iso, v => (v === undefined ? NaN : v)),
          aniso: Array.from(aniso, v => (v === undefined ? NaN : v)), method, couplings: nmr ? nmr.couplings : undefined };
        i = j - 1;
        continue;
      }
      if (c === 'T' && t.startsWith('Total nuclear spin-spin coupling J (Hz)')) {
        const r = readLowerTriangle(L, i + 1);
        if (r.pairs.length) { if (!nmr) nmr = { iso: null, aniso: null, method: null }; nmr.couplings = r.pairs; }
        i = r.end - 1;
        continue;
      }
      if (c === 'R' && t.startsWith('Rotatory Strengths (R) in cgs')) {
        const head = tok(L[i + 1] || '');
        const k = head.findIndex(s => /^R\((length|velocity)\)$/i.test(s));
        if (k > 0) {
          const map = new Map();
          let j = i + 2;
          for (; j < n; j++) { const r = tok(L[j]); if (r.length <= k || !isIntTok(r[0])) break; map.set(int(r[0]), num(r[k])); }
          if (/velocity/i.test(head[k])) rotVel = map; else rotLen = map;
          i = j - 1;
        }
        continue;
      }
      if (c === '-') {
        if (/^-{3,}\s*$/.test(t)) {
          const nx = i + 1 < n ? L[i + 1].trimStart() : '';
          if (nx[0] === '#') {
            const dashLen = t.trimEnd().length;
            const parts = [];
            let j = i + 1;
            for (; j < n && !/^\s*-{3,}\s*$/.test(L[j]); j++) parts.push(L[j].replace(/^ /, ''));
            let route = parts[0] || '';
            for (let k = 1; k < parts.length; k++) route += (parts[k - 1].length < dashLen ? ' ' : '') + parts[k];
            route = route.replace(/\s+/g, ' ').trim();
            routes++;
            if (routes > 1) job++;
            const mb = applyGaussianRoute(route, meta, routes === 1);
            if (routes === 1) { meta.route = route; routeMB = mb; }
            const jt = { jobTypes: [], solvent: 'x' };
            applyGaussianRoute(route, jt, false);
            tdGeomJob = jt.jobTypes.indexOf('td') >= 0 && jt.jobTypes.some(x => x === 'opt' || x === 'scan' || x === 'irc');
            ccKind = /\bqcisd/i.test(route) ? 'QCISD' : /\bccd\b/i.test(route) ? 'CCD' : 'CCSD';
            lastCorr = null;
            oniomJob = /\boniom\s*\(/i.test(route);
            oniomOpt = oniomJob && jt.jobTypes.indexOf('opt') >= 0;
            if (oniomJob && !oniomRoute) oniomRoute = route;
            if (oniomJob && /embed/i.test(route) && !/noembed/i.test(route)) oniomEE = true;
            sink.priority = oniomJob ? ['ONIOM'] : ENERGY_PRIORITY;
            sink.splitDup = !oniomJob;
            nReal = 0; closeStep();
            jobsT.open(route.length > 100 ? route.slice(0, 97) + '...' : route, jobKinds(jt.jobTypes));
            if (jt.jobTypes.indexOf('irc') >= 0) ircJobs.add(jobsT.cur);
            scanCoords = [];
            routeEnd = j; prevDash = lastDash; lastDash = j; i = j;
          } else { prevDash = lastDash; lastDash = i; }
        }
        continue;
      }
      if (oniomJob && c === 'S' && t.startsWith('Structure from the checkpoint')) {
        const s = readChkStructure(L, i + 1);    // rows are left to the main loop (charge lines are read there too)
        if (s.atoms.length && !oniomChk) oniomChk = s;
      }
      if (oniomJob && c === 'S' && t.startsWith('Symbolic Z-matrix')) {
        const echo = readOniomEcho(L, i + 1);   // rows are left to the main loop (charge lines are read there too)
        if (echo.atoms.length) nReal = echo.atoms.length;
        if ((echo.atoms.length || echo.cm.length) && (!oniomEcho || (!oniomEcho.atoms.length && echo.atoms.length))) oniomEcho = echo;
      }
      if (!titleSet && (t.startsWith('Symbolic Z-matrix') || t.startsWith('Charge =') ||
          t.startsWith('Structure from the checkpoint') || t.startsWith('Redundant internal coordinates') ||
          t.startsWith('Z-Matrix taken from the checkpoint'))) {
        if (lastDash === i - 1 && prevDash > routeEnd && lastDash - prevDash > 1 && lastDash - prevDash <= 6) {
          ds.title = L.slice(prevDash + 1, lastDash).map(s => s.trim()).join(' ').trim() || ds.title;
        }
        titleSet = true;
      }
      if (t.startsWith('Charge =')) {
        const m = /^Charge\s*=\s*(-?\d+)\s+Multiplicity\s*=\s*(\d+)/.exec(t);
        if (m && meta.charge === null) { meta.charge = int(m[1]); meta.multiplicity = int(m[2]); }
        continue;
      }
      if (t.startsWith('Standard orientation:') || t.startsWith('Input orientation:') ||
          t.startsWith('Z-Matrix orientation:')) {
        const rank = c === 'S' ? 2 : c === 'I' ? 1 : 0;
        const g = readGaussOrientation(L, i);
        if (g.complete) {
          if (g.Z.length && oniomJob) {
            if (!nReal) nReal = g.Z.length;   // the real system is printed first
            if (g.Z.length === nReal) oniomGeometry(g, rank);
          } else if (g.Z.length) { sink.geometry(g.Z, g.xyz, rank); jobsT.geometry(); }
          if (g.tv.length) cellVecs = g.tv;
        }
        i = g.end;
        continue;
      }
      if (t.startsWith('SCF Done:')) {
        const m = /E\(([^)]+)\)\s*=\s*(\S+)/.exec(t);
        if (m) {
          scfMethod = m[1];
          if (!oniomJob) sink.energy('SCF', num(m[2]));
          else if (oniomSub && !(oniomOpt && oniomEDone)) sink.energy(oniomSub, num(m[2]));
        }
        continue;
      }
      if (t.indexOf('EUMP2') >= 0) { const m = /EUMP2\s*=\s*(\S+)/.exec(t); if (m) E('MP2', num(m[1])); }
      if (t.indexOf('EUMP3') >= 0) { const m = /EUMP3\s*=\s*(\S+)/.exec(t); if (m) E('MP3', num(m[1])); }
      // " E4(SDTQ)= -0.3267D-02   UMP4(SDTQ)= -0.75016068D+02": DQ, SDQ and SDTQ in that order, each kept ('MP4' is the full
      // SDTQ); the best one printed is the energy of an MP4 job
      if (t.indexOf('UMP4(') >= 0) {
        const m = /UMP4\((DQ|SDQ|SDTQ)\)\s*=\s*(\S+)/.exec(t);
        if (m) { E(m[1] === 'SDTQ' ? 'MP4' : 'MP4(' + m[1] + ')', num(m[2])); continue; }
      }
      if (t.startsWith('DEMP5')) { const m = /\bMP5\s*=\s*(\S+)/.exec(t); if (m) E('MP5', num(m[1])); continue; }
      if (t.startsWith('CCSD(T)=')) { E('CCSD(T)', num(t.split('=')[1])); continue; }
      if (t.startsWith('QCISD(T)=')) { E('QCISD(T)', num(t.split('=')[1])); continue; }
      if (t.startsWith('DE(Corr)=')) { const m = /\sE\(CORR\)\s*=\s*(\S+)/i.exec(t); if (m) lastCorr = num(m[1]); continue; }
      // the converged CCSD / CCD / QCISD energy: G16 "Wavefunction amplitudes converged. E(Corr)= ...", G09 "CI/CC converged
      // in N iterations" after the last "DE(Corr)= ... E(CORR)= ..." line
      if (t.startsWith('Wavefunction amplitudes converged.')) {
        const m = /E\(Corr\)\s*=\s*(\S+)/i.exec(t);
        if (m) E(ccKind, num(m[1]));
        lastCorr = null;
        continue;
      }
      if (t.startsWith('CI/CC converged')) { if (lastCorr !== null) E(ccKind, lastCorr); lastCorr = null; continue; }
      if (t.startsWith('E2(')) {
        const m = /^E2\(([^)]+)\)\s*=\s*\S+\s+E\(\1\)\s*=\s*(\S+)/.exec(t);
        if (m) E('DH', num(m[2]));
        continue;
      }
      if (c === 'O' && t.startsWith('ONIOM:')) {
        if (t.startsWith('ONIOM: extrapolated energy')) {
          if (!(oniomOpt && oniomEDone)) { sink.energy('ONIOM', num(t.split('=')[1])); if (oniomJob) oniomEDone = true; }
          continue;
        }
        let m = /^ONIOM: gridpoint\s+\d+\s+method:\s*(\S+)\s+system:\s*(\S+)\s+energy:\s*(\S+)/.exec(t);
        if (m) { if (!(oniomOpt && oniomEDone)) sink.energy(oniomKey(m[1], m[2]), num(m[3])); continue; }
        m = /^ONIOM: generating point\s+\d+\s*--\s*(\w+)\s+level\s+on\s+(\w+)\s+system/i.exec(t);
        if (m) { oniomSub = oniomKey(m[1], m[2]); continue; }
        if (t.startsWith('ONIOM: calculating energy')) { oniomSub = null; continue; }
        m = /^ONIOM: Cut between\s+([A-Za-z]+)\s*\/\s*([A-Za-z]+)\s+(\d+)\s+and\s+([A-Za-z]+)\s+(\d+)(?:\s+factor=\s*(\S+))?/.exec(t);
        if (m && oniomRoute && !(oniomEcho && oniomEcho.atoms.length)) { if (!oniomCuts.some(x => x.atom === int(m[3]) - 1)) oniomCuts.push({ atom: int(m[3]) - 1, bondedTo: int(m[5]) - 1, element: m[2], factor: m[6] ? num(m[6]) : null }); continue; }
        if (t.startsWith('ONIOM: generating point charges')) oniomEE = true;
        continue;
      }
      if (c === 'E' && t.startsWith('Energy=') && /^Energy=\s*\S+\s+NIter=/.test(t)) {
        // molecular mechanics (Amber, UFF, Dreiding): the total of an MM job, a layer energy in ONIOM
        const v = num(tok(t)[1]);
        const f = sink.current();
        if (!oniomJob) sink.energy('MM', v);
        else if (oniomSub) { if (!(oniomOpt && oniomEDone)) sink.energy(oniomSub, v); }
        else if (f && !(oniomOpt && oniomEDone) && Number.isFinite(v)) (oniomMM.get(f) || oniomMM.set(f, []).get(f)).push(v);
        continue;
      }
      if (c === 'N' && t.startsWith('NAtoms=') && t.indexOf('NQM=') > 0) {
        if (oniomJob && !nReal) { const m = /^NAtoms=\s*(\d+)/.exec(t); if (m) nReal = int(m[1]); }
        const q = oniomJob && !oniomCounts && /^NAtoms=\s*(\d+)\s+NQM=\s*(\d+)/.exec(t);
        if (q) oniomCounts = { nAtoms: int(q[1]), nqm: int(q[2]) };
        continue;
      }
      if (c === 'N' && oniomSub && t.startsWith('NAtoms=') && t.indexOf('NActive=') > 0) {
        // "NAtoms=  127 NActive=   25 NUniq= ..." in a sub-calculation: the atoms of its system (QM and link atoms)
        const m = /NActive=\s*(\d+)/.exec(t), p = popOf(oniomSub);
        if (m && p.nActive === null) p.nActive = int(m[1]);
        continue;
      }
      if (t.startsWith('Counterpoise') && /corrected energy\s*=/.test(t)) { E('CP', num(t.split('=')[1])); continue; }
      if (t.startsWith('Total Energy, E(')) {
        // G16 "E(TD-HF/TD-DFT)", G09 "E(TD-HF/TD-KS)", CIS variants
        const m = /^Total Energy, E\([^)]*\)\s*=\s*(\S+)/.exec(t);
        if (m && !oniomJob) { sink.energy('TD', num(m[1])); if (tdGeomJob && sink.current()) tdFrames.add(sink.current()); }
        continue;
      }
      if (t.startsWith('ITU=')) {
        // "Eigenvalues ---" rows of the step's Hessian (Berny optimiser): count the negative ones; a line such as
        // "Use linear search instead of GDIIS." may come first
        let neg = 0, seen = false, j = i + 1;
        while (j < n && j < i + 4 && !L[j].trim().startsWith('Eigenvalues ---')) j++;
        for (; j < n && j < i + 2000; j++) {
          const s = L[j].trim();
          if (!s.startsWith('Eigenvalues ---')) break;
          seen = true;
          for (const x of s.slice(15).match(/-?\d*\.\d+(?:[DdEe][-+]?\d+)?/g) || []) if (num(x) < 0) neg++;
        }
        if (seen) { negEig = neg; i = j - 1; }
        continue;
      }
      if (t.startsWith('Item') && /Value\s+Threshold\s+Converged/.test(t)) {
        const opt = {};
        if (negEig !== null) { opt.negEig = negEig; negEig = null; }
        for (let j = i + 1; j < Math.min(n, i + 8); j++) {
          const s = L[j].trim();
          const m = /^(Maximum|RMS)\s+(Force|Displacement)\s+(\S+)\s+(\S+)\s+(YES|NO)/.exec(s);
          if (m) {
            const key = (m[1] === 'RMS' ? 'rms' : 'max') + (m[2] === 'Force' ? 'Force' : 'Disp');
            opt[key] = num(m[3]); opt[key + 'Thr'] = num(m[4]); opt[key + 'Conv'] = m[5] === 'YES';
            continue;
          }
          const p = /^Predicted change in Energy\s*=\s*(\S+)/.exec(s);
          if (p) opt.predictedChange = num(p[1]);
          break;
        }
        const f = sink.current();
        if (f) f.opt = opt;
        closeStep();
        continue;
      }
      if (t.startsWith('Step number')) {
        closeStep();
        const m = /^Step number\s+(\d+)\s+out of a maximum of\s+(\d+)(?:\s+on scan point\s+(\d+)\s+out of\s+(\d+))?/.exec(t);
        const f = sink.current();
        if (m && f) f.label = m[3] ? 'Scan point ' + m[3] + ', step ' + m[1] : 'Step ' + m[1];
        continue;
      }
      if (t.startsWith('Point Number:')) {
        closeStep();
        const m = /^Point Number:\s*(\d+)\s+Path Number:\s*(\d+)/.exec(t);
        const f = sink.current();
        if (m && f) {
          f.label = 'IRC path ' + m[2] + ', point ' + m[1];
          let coord = null;
          for (let j = i + 1; j < Math.min(n, i + 6); j++) {
            const q = /NET REACTION COORDINATE UP TO THIS POINT\s*=\s*(\S+)/i.exec(L[j]);
            if (q) { coord = Math.abs(num(q[1])); break; }
          }
          // signed reaction coordinate: path 2 (reverse direction) is negative, the TS is 0
          f.irc = { path: int(m[2]), point: int(m[1]), coord: coord === null ? null : (int(m[2]) === 2 ? -coord : coord) };
        }
        continue;
      }
      if (t.startsWith('Optimization completed')) {
        closeStep();
        const f = sink.current();
        if (f) {
          f.converged = true;
          // provisional scan value from the geometry; the "Optimized Parameters" table that follows replaces it
          scanCoords.forEach((sc, k) => setScan(f, scanCoords, k, geomValue(sc, f.xyz)));
        }
        continue;
      }
      if (t.startsWith('Mulliken')) {
        const m = /^Mulliken (atomic )?(charges|spin densities)( and spin densities)?:\s*$/.exec(t.trimEnd());
        if (m) {
          const two = !!m[3];
          const r = readIndexedTable(L, i + 2, two ? 2 : 1);
          if (r.cols[0].length) {
            if (m[2] === 'spin densities') spins.Mulliken = r.cols[0];
            else { charges.Mulliken = r.cols[0]; if (two) spins.Mulliken = r.cols[1]; }
            if (oniomSub && m[2] === 'charges') {   // the population analysis of an ONIOM sub-calculation: its atoms
              const p = popOf(oniomSub), q = r.cols[0];
              if (!p.q) {
                p.q = Float64Array.from(q); p.el = []; p.changed = new Uint8Array(q.length); p.tables = 1;
                for (let j = i + 2; j < r.end; j++) p.el.push(tok(L[j])[1]);
              } else if (q.length === p.q.length) {
                for (let k = 0; k < q.length; k++) if (q[k] !== p.q[k]) p.changed[k] = 1;
                p.tables++;
              }
            }
          }
          i = r.end - 1;
        }
        continue;
      }
      if (t.startsWith('ESP charges:')) {
        const r = readIndexedTable(L, i + 2, 1);
        if (r.cols[0].length) charges.ESP = r.cols[0];
        i = r.end - 1;
        continue;
      }
      if (t.startsWith('Charges from ESP fit')) {
        let j = i + 1;
        while (j < n && /^\s*Charge\s*=/.test(L[j])) j++;
        if (j < n && /^\s*1(\s+\d+)*\s*$/.test(L[j])) {
          const r = readIndexedTable(L, j + 1, 1);
          if (r.cols[0].length) charges.ESP = r.cols[0];
          i = r.end - 1;
        }
        continue;
      }
      if (t.startsWith('APT ')) {   // G16 "APT charges:", G09 "APT atomic charges:"
        if (/^APT (atomic )?charges:\s*$/.test(t)) {
          const r = readIndexedTable(L, i + 2, 1);
          if (r.cols[0].length) charges.APT = r.cols[0];
          i = r.end - 1;
        }
        continue;
      }
      if (t.startsWith('Hirshfeld charges') || t.startsWith('Hirshfeld spin densities')) {
        // "... with hydrogens summed into heavy atoms:" (Q-H / Q-CM5 per heavy atom) must not replace the table
        if (t.indexOf('summed') >= 0) continue;
        const head = tok(L[i + 1] || '');
        const iq = head.indexOf('Q-H'), is = head.indexOf('S-H'), ic = head.indexOf('Q-CM5');
        const q = [], s = [], cm5 = [];
        let j = i + 2;
        for (; j < n; j++) {
          const r = tok(L[j]);
          if (r.length < 3 || !isIntTok(r[0])) break;
          if (iq >= 0) q.push(num(r[2 + iq]));
          if (is >= 0) s.push(num(r[2 + is]));
          if (ic >= 0) cm5.push(num(r[2 + ic]));
        }
        if (q.length) charges.Hirshfeld = q;
        if (cm5.length) charges.CM5 = cm5;
        if (s.some(v => Math.abs(v) > 1e-6)) spins.Hirshfeld = s;
        i = j - 1;
        continue;
      }
      if (t.indexOf('N A T U R A L   A T O M I C') >= 0) { nboSpin = false; continue; }
      if (c === '*' && /^\*+\s+(Alpha|Beta)\s+spin orbitals\s+\*+/.test(t)) { nboSpin = true; continue; }
      if (t.startsWith('Summary of Natural Population Analysis')) {
        let j = i + 1;
        while (j < n && !/^\s*-{10,}\s*$/.test(L[j]) && j < i + 8) j++;
        const q = [];
        for (j++; j < n; j++) {
          const r = tok(L[j]);
          if (r.length < 3 || !isIntTok(r[1])) break;
          q[int(r[1]) - 1] = num(r[2]);
        }
        if (q.length && !nboSpin) charges.NPA = Array.from(q, v => (v === undefined ? NaN : v));
        i = j - 1;
        continue;
      }
      if (t.startsWith('Dipole moment (field-independent basis, Debye)')) {
        const m = /X=\s*(\S+)\s+Y=\s*(\S+)\s+Z=\s*(\S+)\s+Tot=\s*(\S+)/.exec(L[i + 1] || '');
        if (m) ds.dipole = { x: num(m[1]), y: num(m[2]), z: num(m[3]), total: num(m[4]) };
        continue;
      }
      if (t.startsWith('Excitation energies and oscillator strengths')) { states = []; continue; }
      if (t.startsWith('Excited State')) {
        // F8.2 wavelength: >= 10000 nm runs into " eV", >= 100000 nm overflows to "********"
        const m = /^Excited State\s+(\d+):\s+(\S+)\s+(\S+?)\s*eV\s*(\S+?)\s*nm\s+f=\s*(\S+)(?:\s+<S\*\*2>=\s*(\S+))?/.exec(t);
        if (m) {
          if (!states) states = [];
          const ev = num(m[3]);
          let nm = num(m[4]);
          if (!Number.isFinite(nm)) nm = ev > 0 ? EV2NM / ev : NaN;
          states.push({ index: int(m[1]), energyEV: ev, wavelengthNM: nm, f: num(m[5]), label: m[2],
            s2: m[6] !== undefined ? num(m[6]) : null, transitions: [] });
          inState = true;
        }
        continue;
      }
      if (t.startsWith('Harmonic frequencies')) {
        const r = parseGaussFreqSection(L, i);
        if (r.rec.freqs.length) freqSets.push({ rec: r.rec, frameIndex: ds.frames.length - 1, job });
        if (r.rec.stopped) {
          ds.warnings.push('The vibrations section stopped early, after ' + r.rec.freqs.length + ' vibrations, at a line it did not understand' +
            (r.rec.stopped.text ? ': "' + r.rec.stopped.text + '"' : ' (an atom row of an unexpected form)') + '. Vibrations printed after it were not read.');
        }
        i = r.end - 1;
        continue;
      }
      if (t.startsWith('Temperature')) {
        const m = /^Temperature\s+(\S+)\s+Kelvin/.exec(t);
        if (m) temperature = num(m[1]);
        continue;
      }
      if (t.startsWith('Zero-point correction=')) {
        ds.thermo = newThermo(); ds.thermo.temperature = temperature; ds.thermo.zpe = num(t.split('=')[1]);
        continue;
      }
      if (ds.thermo && (t.startsWith('Thermal correction to') || t.startsWith('Sum of electronic and'))) {
        const v = num(t.split('=')[1]);
        const keys = [['Thermal correction to Energy', 'eCorr'], ['Thermal correction to Enthalpy', 'hCorr'],
          ['Thermal correction to Gibbs Free Energy', 'gCorr'], ['Sum of electronic and zero-point Energies', 'eZPE'],
          ['Sum of electronic and thermal Energies', 'eTotal'], ['Sum of electronic and thermal Enthalpies', 'H'],
          ['Sum of electronic and thermal Free Energies', 'G']];
        for (const [p, k] of keys) if (t.startsWith(p + '=')) ds.thermo[k] = v;
        continue;
      }
      if (t.startsWith('Normal termination')) { jobStatus[job] = true; continue; }
      if (t.startsWith('Error termination')) { jobStatus[job] = false; continue; }
      if (t.startsWith('Gaussian ') && !ds.program) {
        const m = /^Gaussian (\d+), Revision ([A-Za-z]\.\d+)/.exec(t) || /^Gaussian (\d+):\s+\S*?Rev([A-Za-z]\.\d+)/.exec(t);
        if (m) ds.program = 'Gaussian ' + m[1] + ', Revision ' + m[2];
        continue;
      }
      if (t.startsWith('Standard basis:')) { const m = /^Standard basis:\s+(\S+)/.exec(t); if (m) stdBasis = m[1]; continue; }
      if (t.startsWith('Solvent') && !solventOut) {
        const m = /^Solvent\s+:\s*([^,]+),/.exec(t);
        if (m) solventOut = m[1].trim();
      }
    }

    if (!ds.frames.length) throw err('No molecular geometry found (no "Standard orientation" or "Input orientation" block).');
    if (oniomRoute) {
      // MM layer energies printed only as "Energy=" lines (no gridpoint summary): assigned in calculation order
      const lab = oniomLabels(oniomRoute), three = !!lab.M;
      const mmLow = /^(amber|uff|dreiding)/i.test(lab.L || '');
      const order = mmLow ? (three ? ['low-mid', 'low-real'] : ['low-model', 'low-real']) : [];
      for (const [f, vals] of oniomMM) {
        const has = f.energies && Object.keys(f.energies).some(k => /^(high|med|low)-/.test(k));
        if (!has && order.length && vals.length === order.length) {
          f.energies = f.energies || {};
          order.forEach((k, q) => { f.energies[k] = vals[q]; });
        }
      }
      if (!nReal) nReal = ds.frames[0].Z.length;
      const full = ds.frames.find(f => f.Z.length === nReal && f.xyz) || null;
      ds.layers = buildOniomLayers(oniomRoute, oniomEcho, nReal, oniomCuts, oniomEE, { chk: oniomChk, counts: oniomCounts, frame: full, pop: oniomPop });
      ds.meta.qmmm = true;
      const at = oniomEcho && oniomEcho.atoms.length === nReal ? oniomEcho.atoms : oniomChk && oniomChk.atoms.length === nReal ? oniomChk.atoms : null;
      if (at && at.some(a => a.props.PDBName || a.props.ResName)) {
        const names = at.map(a => a.props.PDBName || null), res = at.map(a => a.props.ResName || null), ids = at.map(a => (a.props.ResNum ? int(a.props.ResNum) : null));
        for (const f of ds.frames) if (f.Z.length === nReal) { f.atomNames = names; f.resNames = res; f.resIds = ids; }
      }
      const fromChk = !(oniomEcho && oniomEcho.atoms.length) && (oniomChk || /geom\s*=\s*\(?[\w,]*\b(all)?check/i.test(oniomRoute));
      const series = meta.oldchk ? ' (this job continued from ' + meta.oldchk + ')' : '';
      const why = fromChk ? 'Gaussian took this molecule from the checkpoint file (Geom=Check), so the output does not list which atoms are in which ONIOM layer'
        : 'The molecule echo is not in Cartesian form, so the ONIOM layers could not be read from it';
      const nl = ds.layers.linkAtoms.length, links = nl + (nl === 1 ? ' link atom' : ' link atoms');
      if (ds.layers.source === 'population') {
        const P = ds.layers.population, s = k => (k === 1 ? '' : 's');
        const who = (P.mid ? `${P.qm} high-layer and ${P.mid} medium-layer atoms` : `${P.qm} QM atom${s(P.qm)}`) +
          (P.links === 1 ? ' and the atom replaced by a link atom' : P.links ? ` and the ${P.links} atoms replaced by link atoms` : '');
        const how = P.how === 'changed'
          ? `with electronic embedding the other atoms keep their MM charges there, and the Mulliken charges that change from one geometry to the next (${P.tables} tables) are those of the ${who}`
          : `only the ${who} ${P.qm + P.mid + P.links === 1 ? 'has a Mulliken charge' : 'have Mulliken charges'} there`;
        ds.warnings.push(`${why}. BlueMolecule read them from the ${P.mid ? 'high- and medium-level calculations' : 'high-level calculation on the model system'}: ${how}, ` +
          `as many atoms as Gaussian counted. To use the layers of another file of the series${series} instead, open it and take its layers in the Info tab (ONIOM layers).`);
      } else if (ds.layers.source === 'rebuilt') {
        ds.warnings.push(why + '. BlueMolecule rebuilt them from the ' + links + ': the ' + ds.layers.rebuilt.qm +
          ' atoms bonded to their QM side are as many as the QM atoms Gaussian counted. To use the original assignment instead, open the input or output that started the series' +
          series + ' and take its layers in the Info tab (ONIOM layers).');
      } else if (!ds.layers.atomLayer) {
        ds.warnings.push(why + (nl ? ', and the ' + links + (nl === 1 ? ' alone does' : ' alone do') + ' not give the whole QM region' : '') +
          '. Open the input or output that started the series' + series + ', then take its layers in the Info tab (ONIOM layers).');
      }
      const bad = ds.frames.filter(f => oniomFrames.has(f) && f.Z.length !== nReal).length;
      if (bad) ds.warnings.push(bad + ' ONIOM geometries do not have the ' + nReal + ' atoms of the real system.');
    }
    for (const f of tdFrames) if (f.energies && f.energies.TD !== undefined) f.energy = f.energies.TD;
    // ---- job list, scan / IRC bookkeeping, energy change between optimisation steps
    const jobList = jobsT.finish() || [{ index: 0, label: meta.route, kinds: jobKinds(meta.jobTypes), frameStart: 0,
      frameEnd: ds.frames.length - 1, normalTermination: null }];
    jobList.forEach((jb, k) => {
      jb.normalTermination = jobStatus[k] === undefined ? null : jobStatus[k];
      const fr = ds.frames.slice(jb.frameStart, jb.frameEnd + 1);
      if (fr.some(f => f.scan) && jb.kinds.indexOf('scan') < 0) {
        jb.kinds = jb.kinds.filter(x => x !== 'sp').concat(['scan']);
        addJobType(meta, 'scan');
      }
      if (ircJobs.has(jb) && fr.length && !fr[0].irc) fr[0].irc = { path: 0, point: 0, coord: 0 };
      for (let q = 1; q < fr.length; q++) {
        const a = fr[q - 1], b = fr[q];
        if (b.opt && a.opt && !a.converged && b.opt.energyChange === undefined && Number.isFinite(a.energy) &&
            Number.isFinite(b.energy)) b.opt.energyChange = b.energy - a.energy;
      }
    });
    ds.jobs = jobList;
    if (states && (rotLen || rotVel)) {
      for (const st of states) {
        if (rotLen && rotLen.has(st.index)) st.R = rotLen.get(st.index);
        if (rotVel && rotVel.has(st.index)) st.Rvel = rotVel.get(st.index);
      }
    }
    if (nmr) { ds.nmr = nmr; if (nmr.couplings === undefined) delete nmr.couplings; addJobType(meta, 'nmr'); }
    if (cellVecs) ds.cell = padCell(cellVecs);
    if (freqSets.length) {
      // the last geometry's frequencies; of several sections for it (hpmodes prints two), the one with the most
      // vibrations, then the high-precision one
      const end = last(freqSets);
      let pick = end;
      for (let k = freqSets.length - 2; k >= 0; k--) {
        const s = freqSets[k];
        if (s.job !== end.job || s.frameIndex !== end.frameIndex) break;
        const n = s.rec.freqs.length, np = pick.rec.freqs.length;
        if (n > np || (n === np && s.rec.hp && !pick.rec.hp)) pick = s;
      }
      ds.vibrations = buildVibrations(pick.rec, pick.frameIndex, ds);
    }
    if (states && states.length) ds.excitedStates = states;
    if (Object.keys(charges).length) ds.charges = charges;
    if (Object.keys(spins).length) ds.spins = spins;
    if (eigLast && (eigLast.aOcc.length || eigLast.aVirt.length)) {
      ds.orbitals = { alpha: { energies: eigLast.aOcc.concat(eigLast.aVirt), nOcc: eigLast.aOcc.length },
        beta: eigLast.bOcc ? { energies: eigLast.bOcc.concat(eigLast.bVirt), nOcc: eigLast.bOcc.length } : null };
    }
    if (routeMB) {
      meta.method = routeMB.method ? routeMB.method.toUpperCase() : null;
      meta.basis = routeMB.basis;
    }
    if (!meta.method && scfMethod) meta.method = scfMethod.replace(/^(RO|R|U)(?=[A-Z0-9])/, '');
    if (stdBasis && (!meta.basis || meta.basis.toLowerCase() === stdBasis.toLowerCase())) meta.basis = stdBasis;
    if (solventOut) meta.solvent = solventOut;
    finishJobTypes(meta);
    if (jobStatus.some(s => s === false)) meta.normalTermination = false;
    else if (jobStatus[job] === true) meta.normalTermination = true;
    if (!ds.program) ds.program = 'Gaussian';
    return ds;
  }

  // -------------------------------------------------------------- Z-matrix
  const vsub = (a, b) => [a[0] - b[0], a[1] - b[1], a[2] - b[2]];
  const vadd = (a, b) => [a[0] + b[0], a[1] + b[1], a[2] + b[2]];
  const vscale = (a, s) => [a[0] * s, a[1] * s, a[2] * s];
  const vdot = (a, b) => a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
  const vcross = (a, b) => [a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]];
  const vnorm = a => { const l = Math.hypot(a[0], a[1], a[2]); return l > 1e-12 ? vscale(a, 1 / l) : [0, 0, 0]; };
  const RAD = Math.PI / 180;

  /**
   * BM.zmatrixToCartesian(rows) -> Float64Array(3 * rows.length), Angstrom.
   * rows[i] describes atom i; references are 0-based indices of EARLIER rows:
   *   { xyz: [x, y, z] }                          Cartesian position, used as given
   *   { }                                         first atom: origin
   *   { a, r }                                    bonded to atom a at distance r (second atom: +z axis)
   *   { a, r, b, angle }                          bond angle i-a-b in degrees (third atom: xz plane)
   *   { a, r, b, angle, c, dihedral }             dihedral i-a-b-c in degrees (IUPAC sign convention)
   *   { a, r, b, angle, c, dihedral, mode: ±1 }   Gaussian format code 1: 'dihedral' is a second bond
   *                                               angle i-a-c; the sign of mode selects the side
   * Orientation follows Gaussian: atom 1 at the origin, atom 2 on +z, atom 3 in the xz plane.
   */
  function zmatrixToCartesian(rows) {
    const n = rows.length;
    const X = new Float64Array(3 * n);
    const P = i => {
      if (!(i >= 0 && i < n)) throw err('Z-matrix reference to an undefined atom.');
      return [X[3 * i], X[3 * i + 1], X[3 * i + 2]];
    };
    for (let i = 0; i < n; i++) {
      const r = rows[i] || {};
      let p;
      if (r.xyz) p = [r.xyz[0], r.xyz[1], r.xyz[2]];
      else if (r.a === undefined || r.a === null) p = [0, 0, 0];
      else if (r.b === undefined || r.b === null) { const A = P(r.a); p = [A[0], A[1], A[2] + r.r]; }
      else if (r.c === undefined || r.c === null) {
        const A = P(r.a), B = P(r.b);
        const u = vnorm(vsub(B, A));
        let w = vcross([0, 1, 0], u);
        if (Math.hypot(w[0], w[1], w[2]) < 1e-8) w = vcross([1, 0, 0], u);
        w = vnorm(w);
        const th = r.angle * RAD;
        p = vadd(A, vscale(vadd(vscale(u, Math.cos(th)), vscale(w, Math.sin(th))), r.r));
      } else if (!r.mode) {
        const A = P(r.a), B = P(r.b), C = P(r.c);
        const bc = vnorm(vsub(A, B));
        let nn = vcross(vsub(B, C), bc);
        if (Math.hypot(nn[0], nn[1], nn[2]) < 1e-8) {   // collinear reference atoms
          nn = vcross(bc, Math.abs(bc[0]) < 0.9 ? [1, 0, 0] : [0, 1, 0]);
        }
        nn = vnorm(nn);
        const m = vcross(nn, bc);
        const th = r.angle * RAD, ph = r.dihedral * RAD;
        p = vadd(A, vadd(vscale(bc, -r.r * Math.cos(th)),
          vadd(vscale(m, r.r * Math.sin(th) * Math.cos(ph)), vscale(nn, r.r * Math.sin(th) * Math.sin(ph)))));
      } else {
        const A = P(r.a), u = vnorm(vsub(P(r.b), A)), v = vnorm(vsub(P(r.c), A));
        const g = vdot(u, v), ca = Math.cos(r.angle * RAD), cb = Math.cos(r.dihedral * RAD);
        const den = 1 - g * g || 1e-12;
        const x = (ca - g * cb) / den, y = (cb - g * ca) / den;
        const z2 = 1 - (x * x + y * y + 2 * x * y * g);
        const w = vnorm(vcross(u, v));
        const d = vadd(vadd(vscale(u, x), vscale(v, y)), vscale(w, (r.mode > 0 ? 1 : -1) * Math.sqrt(Math.max(0, z2))));
        p = vadd(A, vscale(d, r.r));
      }
      X[3 * i] = p[0]; X[3 * i + 1] = p[1]; X[3 * i + 2] = p[2];
    }
    return X;
  }

  // ------------------------------------------------------ Gaussian input
  const DUMMY_LABEL = /^(x|xx|du|dum|da|q)(?![a-z])/i;
  const GHOST_LABEL = /^(bq|gh)(?![a-z])|-bq\b|\(.*ghost/i;
  const VARDEF = /^\s*([A-Za-z][\w.]*)\s*(?:=\s*|\s+)([-+]?(?:\d+\.?\d*|\.\d+)(?:[eEdD][-+]?\d+)?)(?:\s+.*)?$/;

  function gaussLabelZ(label, ds) {
    if (GHOST_LABEL.test(label)) return 0;
    const z = BM.el.toZ(label.replace(/\(.*$/, ''));
    if (z < 0) { ds.warnings.push('Unknown element label "' + label + '" read as a dummy atom.'); return 0; }
    return z;
  }

  function parseGaussianInput(text, fileName) {
    const ds = newDataset('gaussian-input', fileName);
    const meta = ds.meta;
    let L = splitLines(text).map(s => s.replace(/\t/g, ' ').replace(/!.*$/, ''));
    const li = L.findIndex(s => /^\s*--link1--\s*$/i.test(s));
    if (li >= 0) {
      if (L.slice(li + 1).some(s => s.trim())) ds.warnings.push('Only the first job of this --Link1-- input is shown.');
      L = L.slice(0, li);
    }
    let i = 0;
    const link0 = [];
    while (i < L.length && !L[i].trim()) i++;
    while (i < L.length && L[i].trim().startsWith('%')) { link0.push(L[i].trim()); i++; }
    while (i < L.length && !L[i].trim()) i++;
    if (i >= L.length || !L[i].trim().startsWith('#')) throw err('No route section (a line starting with "#") was found.');
    const rl = [];
    while (i < L.length && L[i].trim()) rl.push(L[i++].trim());
    const route = rl.join(' ').replace(/\s+/g, ' ').trim();
    meta.route = route;
    meta.link0 = link0;
    const mb = applyGaussianRoute(route, meta, true);
    meta.method = mb.method ? mb.method.toUpperCase() : null;
    meta.basis = mb.basis;
    finishJobTypes(meta);
    if (/geom\s*=\s*\(?[\w,]*\b(all)?check/i.test(route)) {
      throw err('This input takes its geometry from a checkpoint file (Geom=Check), so it contains no coordinates.');
    }
    i++;
    const title = [];
    while (i < L.length && L[i].trim()) title.push(L[i++].trim());
    if (title.length) ds.title = title.join(' ');
    i++;
    while (i < L.length && !L[i].trim()) i++;
    const cm = tok(L[i] || '').map(s => s.replace(/,$/, ''));
    const cmv = cm.join(' ').split(/[\s,]+/).filter(Boolean);
    if (cmv.length < 2 || !cmv.every(isIntTok)) throw err('Expected the charge and multiplicity line after the title.');
    meta.charge = int(cmv[0]); meta.multiplicity = int(cmv[1]);
    i++;
    const spec = [];
    while (i < L.length && L[i].trim()) spec.push(L[i++]);
    if (!spec.length) throw err('The molecule specification is empty.');
    const first = rowTokens(spec[0]);
    let cart = first.length >= 4 && ((isNumTok(first[1]) && isNumTok(first[2]) && isNumTok(first[3])));
    const labels = [], rows = [], specOf = [];
    const tv = [];
    if (cart) {
      for (const s of spec) {
        const t = rowTokens(s);   // "N-N3-0.1592(PDBName=N,ResName=SER,ResNum=1) 0 x y z L": commas inside () kept
        if (!t.length) continue;
        let xyz = null;
        if (t.length >= 5 && isIntTok(t[1]) && int(t[1]) <= 0 && isNumTok(t[2]) && isNumTok(t[3]) && isNumTok(t[4])) {
          xyz = [num(t[2]), num(t[3]), num(t[4])];
        } else if (t.length >= 4 && isNumTok(t[1]) && isNumTok(t[2]) && isNumTok(t[3]) &&
                   !(t.length >= 5 && isIntTok(t[1]) && int(t[1]) > 0 && isIntTok(t[3]))) {
          xyz = [num(t[1]), num(t[2]), num(t[3])];
        } else { cart = false; break; }   // mixed Cartesian / Z-matrix input: use the Z-matrix reader
        if (/^tv$/i.test(t[0])) { tv.push(xyz); continue; }
        labels.push(t[0]); rows.push({ xyz }); specOf.push(s);
      }
      if (!cart) { labels.length = 0; rows.length = 0; tv.length = 0; specOf.length = 0; }
    }
    if (!cart) {
      const zl = [], vl = [];
      let inVars = false;
      for (const s of spec) {
        if (/^\s*(variables|constants)\s*:?\s*$/i.test(s)) { inVars = true; continue; }
        (inVars ? vl : zl).push(s);
      }
      let j = i + 1;
      for (let blk = 0; blk < 2 && j < L.length; blk++) {
        const b = [];
        let k = j;
        while (k < L.length && L[k].trim()) b.push(L[k++]);
        const body = b.filter(s => !/^\s*(variables|constants)\s*:?\s*$/i.test(s));
        if (!b.length || !body.every(s => VARDEF.test(s))) break;
        vl.push(...body);
        j = k + 1;
      }
      const vars = Object.create(null);
      for (const s of vl) {
        const m = VARDEF.exec(s);
        if (m) { vars[m[1]] = num(m[2]); vars[m[1].toLowerCase()] = num(m[2]); }
      }
      const value = tk => {
        if (isNumTok(tk)) return num(tk);
        const m = /^([-+]?)(.+)$/.exec(tk);
        const v = vars[m[2]] !== undefined ? vars[m[2]] : vars[m[2].toLowerCase()];
        if (v === undefined) throw err('Z-matrix variable "' + m[2] + '" is not defined.');
        return m[1] === '-' ? -v : v;
      };
      const ref = (tk, row) => {
        let k;
        if (isIntTok(tk)) k = int(tk) - 1;
        else k = labels.findIndex(l => l === tk);
        if (k < 0 || k >= row) {
          const lo = labels.findIndex(l => l.toLowerCase() === String(tk).toLowerCase());
          if (lo >= 0 && lo < row) k = lo;
          else throw err('Z-matrix line ' + (row + 1) + ' refers to unknown atom "' + tk + '".');
        }
        return k;
      };
      zl.forEach((s, row) => {
        const t = rowTokens(s);
        labels.push(t[0]); specOf.push(null);
        if (t.length >= 5 && t[1] === '0' && isNumTok(t[2]) && isNumTok(t[3]) && isNumTok(t[4])) {
          rows.push({ xyz: [num(t[2]), num(t[3]), num(t[4])] }); return;
        }
        const r = {};
        if (t.length >= 3) { r.a = ref(t[1], row); r.r = value(t[2]); }
        if (t.length >= 5) { r.b = ref(t[3], row); r.angle = value(t[4]); }
        if (t.length >= 7) {
          r.c = ref(t[5], row); r.dihedral = value(t[6]);
          if (t.length >= 8 && isIntTok(t[7]) && int(t[7]) !== 0) r.mode = int(t[7]);
        }
        rows.push(r);
      });
    }
    const X = zmatrixToCartesian(rows);
    const Z = [], xyz = [], kept = [];
    labels.forEach((lab, k) => {
      if (DUMMY_LABEL.test(lab)) return;
      if (/^tv$/i.test(lab)) { tv.push([X[3 * k], X[3 * k + 1], X[3 * k + 2]]); return; }
      Z.push(gaussLabelZ(lab, ds)); xyz.push(X[3 * k], X[3 * k + 1], X[3 * k + 2]); kept.push(k);
    });
    if (!Z.length) throw err('No atoms found in the molecule specification.');
    const f = newFrame(Z, xyz);
    f.label = cart ? 'Input geometry' : 'Input geometry (from Z-matrix)';
    ds.frames.push(f);
    // residue information written by GaussView: "C-CT--0.1(PDBName=CA,ResName=ALA,ResNum=12)"
    const props = kept.map(k => { const a = oniomAtomSpec(labels[k]); return a ? a.props : {}; });
    if (props.some(p => p.PDBName || p.ResName)) {
      f.atomNames = props.map(p => p.PDBName || null); f.resNames = props.map(p => p.ResName || null);
      f.resIds = props.map(p => (p.ResNum && isIntTok(p.ResNum) ? int(p.ResNum) : null));
    }
    if (/\boniom\s*\(/i.test(route)) {
      // ONIOM: layers, link atoms, freeze codes, MM types and charges from the rows (Cartesian input)
      const atoms = cart ? kept.map(k => oniomAtomRow(specOf[k])) : [];
      const three = !!oniomLabels(route).M;
      const order = three ? [['low', 'real'], ['med', 'mid'], ['low', 'mid'], ['high', 'model'], ['med', 'model'], ['low', 'model']] : [['low', 'real'], ['high', 'model'], ['low', 'model']];
      const cmPairs = [];
      for (let q = 0; q + 1 < cmv.length && q / 2 < order.length; q += 2) cmPairs.push({ charge: int(cmv[q]), multiplicity: int(cmv[q + 1]), level: order[q / 2][0], system: order[q / 2][1] });
      const ok = atoms.length === Z.length && atoms.every(Boolean);
      ds.layers = buildOniomLayers(route, { cm: cmPairs, atoms: ok ? atoms : [] }, Z.length, [], /embed/i.test(route) && !/noembed/i.test(route));
      meta.qmmm = true;
      if (!ok) ds.warnings.push('The ONIOM layers of this input could not be read' + (cart ? ' (some atom rows are not in the form "atom [freeze] x y z layer")' : ' (the molecule is a Z-matrix)') + '.');
    }
    if (tv.length) ds.cell = padCell(tv);
    ds.program = null;
    return ds;
  }

  // ------------------------------------------------------------------ ORCA
  function orcaAtomRows(L, start, re) {
    const out = [];
    let j = start;
    for (; j < L.length; j++) {
      const m = re.exec(L[j]);
      if (!m) break;
      out.push(m);
    }
    return { rows: out, end: j };
  }
  const ORCA_CHG = /^\s*(\d+)\s*([A-Za-z]{1,2}[^\s:]*)\s*:\s*(\S+)(?:\s+(\S+))?\s*$/;

  // ------------------------------------------------------------------ ORCA QM/MM and QM/QM2 (multiscale)
  /** "{0:5 10 12:20}" (0-based, inclusive ranges) -> sorted unique indices. */
  function orcaRanges(s) {
    const out = new Set();
    for (const p of String(s).replace(/[{}]/g, ' ').split(/[\s,]+/)) {
      if (!p) continue;
      const m = /^(\d+)(?::(\d+))?$/.exec(p);
      if (!m) continue;
      const a = int(m[1]), b = m[2] !== undefined ? int(m[2]) : a;
      for (let k = Math.min(a, b); k <= Math.max(a, b) && out.size < 5e6; k++) out.add(k);
    }
    return Array.from(out).sort((x, y) => x - y);
  }
  /** QM/MM settings from the input echo: scheme, %qmmm items, the size of an inline *xyz geometry. */
  function orcaQmmmInput(input) {
    const bang = input.filter(x => /^\s*!/.test(x)).map(x => x.replace(/^\s*!/, '').replace(/#.*$/, '')).join(' ');
    const kw = bang.split(/\s+/).filter(Boolean);
    let scheme = null, qm2 = null, mm = false, extra = null;
    for (const k of kw) {
      const u = k.toUpperCase();
      if (u === 'QMMM' || u === 'QM/MM' || u === 'L-QM/MM' || u === 'IONIC-CRYSTAL-QMMM' || u === 'MOL-CRYS-QMMM') { scheme = 'qmmm'; mm = true; if (/CRYS/.test(u)) extra = u; }
      else {
        const m = /^QM\/([A-Z0-9-]+?)(\/MM)?$/.exec(u);
        if (m && m[1] !== 'MM') { qm2 = m[1]; mm = !!m[2]; scheme = mm ? 'qmmm' : 'qm-qm2'; }
      }
    }
    const inp = input.join('\n');
    const at = inp.search(/%qmmm\b/i);
    const items = {};
    if (at >= 0) {
      if (!scheme) scheme = 'qmmm';
      const re = /\{[^}]*\}|"[^"]*"|#[^\n]*|[^\s{}"#]+/g;
      re.lastIndex = at + 5;
      let m, key = null;
      while ((m = re.exec(inp))) {
        const tkn = m[0];
        if (tkn[0] === '#') continue;
        if (key === null) { if (/^end$/i.test(tkn)) break; key = tkn.toLowerCase(); continue; }
        items[key] = tkn[0] === '"' ? tkn.slice(1, -1) : tkn;
        if (tkn[0] === '{') { const q = /\s*end\b/iy; q.lastIndex = re.lastIndex; if (q.exec(inp)) re.lastIndex = q.lastIndex; }
        key = null;
      }
    }
    if (!scheme && /%pointcharges\b/i.test(inp)) scheme = 'external';
    if (!scheme) return null;
    let nInput = 0;
    const g = /^\s*\*\s*xyz\s+-?\d+\s+\d+\s*$/im.exec(inp);
    if (g) {
      for (const line of inp.slice(g.index + g[0].length).split('\n').slice(1)) {
        const tt = line.trim();
        if (!tt) continue;
        if (tt[0] === '*') break;
        if (/^[A-Za-z]/.test(tt)) nInput++;
      }
    }
    const list = k => (items[k] && items[k][0] === '{' ? orcaRanges(items[k]) : null);
    const charge = k => (items[k] !== undefined && isIntTok(items[k]) ? int(items[k]) : null);
    return {
      scheme, qm2, mm, extra, nInput,
      qmAtoms: list('qmatoms'), qm2Atoms: list('qm2atoms'), activeAtoms: list('activeatoms'), fixedAtoms: list('optregion_fixedatoms'),
      ff: items.orcafffilename || null, qm2Method: items.qm2custommethod || null, qm2Basis: items.qm2custombasis || null,
      qm2File: items.qm2customfile || null, embedding: items.embedding ? String(items.embedding).toLowerCase() : null,
      chargeMult: [['total', charge('charge_total'), charge('mult_total')], ['medium', charge('charge_medium'), charge('mult_medium')]]
        .filter(c => c[1] !== null || c[2] !== null).map(c => ({ system: c[0], charge: c[1], multiplicity: c[2] })),
      summary: { nQM: 0, nLink: 0, nActive: 0, nFull: 0, links: [] }
    };
  }
  /** A QM/MM energy line ("E(QM/MM) ... -228.58 Eh", "QM/MM energy: ...", "E(QM2, full system) ...") -> {key, value}. */
  function orcaQmmmEnergy(t) {
    const m = /^(E\s*\([^)]*\)|[A-Za-z][A-Za-z0-9 ,/()+_-]*?)\s*(?:\.{2,}|:|=)\s*(-?\d+\.\d+(?:[eEdD][-+]?\d+)?)\s*(?:Eh|a\.u\.|au|hartrees?)?\s*(?:\S+\s*eV)?\s*$/.exec(t);
    if (!m) return null;
    const lab = m[1].toLowerCase().replace(/\s+/g, ' ').trim();
    if (!/\bqm|\bmm\b|qmmm|multiscale/.test(lab) || /grad|time|charge|dipole|number|thresh|conv|cutoff|radius|scal/.test(lab)) return null;
    const inner = lab.replace(/^e\s*\(\s*/, '').replace(/\s*\)$/, '');
    let key = null;
    if (/qm\/qm2|qm\/xtb|multiscale|qm\/hf|qm\/pbeh|qm\/r2scan/.test(inner)) key = 'QM/QM2';
    else if (/qm2|\blow\b/.test(inner)) key = /full|total system|real|whole|entire/.test(inner) ? 'low-real' : /qm region|qm part|qm subsystem|model|inner|qm atoms|small/.test(inner) ? 'low-model' : 'QM2';
    else if (/qm\/mm|qmmm/.test(inner) && !/interaction|electrostatic|vdw|van der waals|coupling/.test(inner)) key = 'QM/MM';
    else if (/qm.?mm|interaction|electrostatic|vdw|coupling/.test(inner) && /qm/.test(inner) && /mm/.test(inner)) key = 'QM-MM';
    else if (/^(total )?mm( energy| part| total)*$|^mm\b/.test(inner)) key = 'MM';
    else if (/^qm\b/.test(inner) && /region|model|part|subsystem|high/.test(inner)) key = 'high-model';
    else if (/^qm( energy| total)*$/.test(inner)) key = 'QM';
    else if (/total/.test(inner) && /qm/.test(inner) && /mm/.test(inner)) key = 'QM/MM';
    return key ? { key, value: num(m[2]) } : null;
  }
  /** Summary counts / link atom rows of the QM/MM setup print; true when the line was one of them. */
  function orcaQmmmSummary(t, q) {
    const s = q.summary;
    let m = /number of qm atoms\s*(?:\.{2,}|:|=)\s*(\d+)/i.exec(t) || /^qm region size\s*(?:\.{2,}|:|=)\s*(\d+)/i.exec(t) ||
      /^size of qm subsystem\s*(?:\.{2,}|:|=)\s*(\d+)$/i.exec(t);   // ORCA 6 ("... plus link atoms" is another line)
    if (m) { s.nQM = int(m[1]); return true; }
    if ((m = /number of link atoms\s*(?:\.{2,}|:|=)\s*(\d+)/i.exec(t))) { s.nLink = int(m[1]); return true; }
    if ((m = /(?:number of active atoms|optimized atoms\s*=\s*activeregion)\s*(?:\.{2,}|:|=)\s*(\d+)/i.exec(t))) { s.nActive = int(m[1]); return true; }
    if ((m = /(?:total system size|total number of atoms|number of atoms in (?:the )?(?:full|total|whole|real) system|^size of qm\/?mm system)\s*(?:\.{2,}|:|=)\s*(\d+)/i.exec(t))) { s.nFull = int(m[1]); return true; }
    if ((m = /link atom.*?\bqm atom\s*(\d+)\b.*?\bmm atom\s*(\d+)/i.exec(t))) { s.links.push({ bondedTo: int(m[1]), atom: int(m[2]) }); return true; }
    return false;
  }
  /**
   * BM.qmmmMatch(fullZ, fullXyz, subZ, subXyz, {shift}) -> {qm, links, shift} | null
   * Finds the atoms of a QM-region print (QM atoms + link atoms) in the full system: QM atoms coincide with
   * full-system atoms (after one common translation when shift is true: CP2K centres the QM box), the others are
   * link atoms, each placed on a bond between a QM atom (bondedTo) and the MM atom it replaces (atom).
   */
  function qmmmMatch(fZ, fX, sZ, sX, opts) {
    const N = fZ.length, n = sZ.length, tol = (opts && opts.tol) || 2e-3, cell = 0.05;
    if (!n || !N) return null;
    const grid = new Map();
    const key = (x, y, z) => Math.floor(x / cell) + ',' + Math.floor(y / cell) + ',' + Math.floor(z / cell);
    for (let a = 0; a < N; a++) {
      const k = key(fX[3 * a], fX[3 * a + 1], fX[3 * a + 2]);
      const l = grid.get(k); if (l) l.push(a); else grid.set(k, [a]);
    }
    const find = (z, x, y, w) => {
      const i = Math.floor(x / cell), j = Math.floor(y / cell), k = Math.floor(w / cell);
      for (let di = -1; di <= 1; di++) for (let dj = -1; dj <= 1; dj++) for (let dk = -1; dk <= 1; dk++) {
        const l = grid.get((i + di) + ',' + (j + dj) + ',' + (k + dk));
        if (l) for (const a of l) if (fZ[a] === z && Math.abs(fX[3 * a] - x) < tol && Math.abs(fX[3 * a + 1] - y) < tol && Math.abs(fX[3 * a + 2] - w) < tol) return a;
      }
      return -1;
    };
    const tryShift = (tx, ty, tz) => {
      const hit = new Array(n);
      let c = 0;
      for (let b = 0; b < n; b++) { hit[b] = find(sZ[b], sX[3 * b] + tx, sX[3 * b + 1] + ty, sX[3 * b + 2] + tz); if (hit[b] >= 0) c++; }
      return { hit, c };
    };
    let best = tryShift(0, 0, 0), t = [0, 0, 0];
    if (opts && opts.shift && best.c < n / 2) {
      // anchor: the heaviest atom of the QM print against every full-system atom of the same element
      let b0 = 0;
      for (let b = 1; b < n; b++) if (sZ[b] > sZ[b0]) b0 = b;
      for (let a = 0; a < N && best.c < n; a++) {
        if (fZ[a] !== sZ[b0]) continue;
        const tx = fX[3 * a] - sX[3 * b0], ty = fX[3 * a + 1] - sX[3 * b0 + 1], tz = fX[3 * a + 2] - sX[3 * b0 + 2];
        const r = tryShift(tx, ty, tz);
        if (r.c > best.c) { best = r; t = [tx, ty, tz]; }
      }
    }
    if (best.c < Math.max(1, n / 2)) return null;
    const qm = best.hit.filter(a => a >= 0).sort((x, y) => x - y), inQM = new Set(qm), links = [];
    for (let b = 0; b < n; b++) {
      if (best.hit[b] >= 0) continue;
      const x = sX[3 * b] + t[0], y = sX[3 * b + 1] + t[1], z = sX[3 * b + 2] + t[2];
      let host = -1, hd = Infinity, q = -1, qd = Infinity;
      for (let a = 0; a < N; a++) {
        const dx = fX[3 * a] - x, dy = fX[3 * a + 1] - y, dz = fX[3 * a + 2] - z, d = dx * dx + dy * dy + dz * dz;
        if (inQM.has(a)) { if (d < qd) { qd = d; q = a; } } else if (d < hd) { hd = d; host = a; }
      }
      if (host >= 0) links.push({ atom: host, bondedTo: q, element: BM.el.symbol(sZ[b]) || 'H', type: null });
    }
    return { qm, links, shift: t };
  }
  /** Layers of an ORCA multiscale run on the full-system frame f (sub: the QM-region print of the same step). */
  function orcaLayers(q, f, sub, meta, subsystemOnly) {
    const N = f.Z.length, X = f.xyz;
    const d2 = (a, xyz, b) => { const dx = X[3 * a] - xyz[3 * b], dy = X[3 * a + 1] - xyz[3 * b + 1], dz = X[3 * a + 2] - xyz[3 * b + 2]; return dx * dx + dy * dy + dz * dz; };
    let qm = q.qmAtoms ? q.qmAtoms.filter(k => k < N) : null;
    let links = [];
    if (!subsystemOnly && sub && sub.Z.length < N) {
      // QM atoms are printed at their full-system positions; the atoms that match none are the link atoms
      const r = qmmmMatch(f.Z, X, sub.Z, sub.xyz, { shift: false });
      if (r) { if (!qm) qm = r.qm; links = r.links; }
    }
    if (!links.length && q.summary.links.length) links = q.summary.links.filter(l => l.atom < N).map(l => ({ atom: l.atom, bondedTo: l.bondedTo, element: 'H', type: null }));
    let rebuilt = false;
    if (!qm && !subsystemOnly && X && links.length && q.summary.nQM > 0) {
      // no QMAtoms in the input echo and no QM-region print to match: the atoms bonded to the QM ends of the link atoms,
      // when they are as many as the QM atoms ORCA counts
      const r = qmRegionFromLinks(f.Z, X, links, q.summary.nQM);
      if (r) { qm = r; rebuilt = true; }
    }
    if (!links.length && subsystemOnly && X && q.summary.nQM > 0 && q.summary.nLink > 0 && N === q.summary.nQM + q.summary.nLink) {
      // only the QM region is printed: its QM atoms, then its link atoms, each bonded to the nearest QM atom
      for (let a = q.summary.nQM; a < N; a++) {
        let best = -1, bd = Infinity;
        for (let b = 0; b < q.summary.nQM; b++) { const d = d2(a, X, b); if (d < bd) { bd = d; best = b; } }
        links.push({ atom: a, bondedTo: best, element: 'H', type: null });
      }
    }
    if (!links.length && qm && !subsystemOnly && X) {
      // covalent bonds that cross the QM boundary carry the link atoms
      const inQM = new Set(qm);
      for (const a of qm) for (let b = 0; b < N; b++) {
        if (inQM.has(b)) continue;
        const rc = (BM.el.covalent[f.Z[a]] || 0.8) + (BM.el.covalent[f.Z[b]] || 0.8);
        if (d2(a, X, b) < 1.15 * 1.15 * rc * rc && f.Z[a] > 1 && f.Z[b] > 1) links.push({ atom: b, bondedTo: a, element: 'H', type: null });
      }
    }
    const inQM = new Set(qm || []), inQM2 = new Set(q.qm2Atoms || []);
    const outer = q.scheme === 'qm-qm2' ? 'M' : 'L';
    // QM atoms not known: the layers stay unknown (null) rather than every atom shown as MM
    const atomLayer = subsystemOnly ? new Array(N).fill('H') : qm ? Array.from({ length: N }, (_, a) => (inQM.has(a) ? 'H' : inQM2.has(a) ? 'M' : outer)) : null;
    let active = null, frozen = null;
    if (q.activeAtoms && !subsystemOnly) { const s = new Set(q.activeAtoms); active = Array.from({ length: N }, (_, a) => s.has(a)); frozen = active.map(x => !x); }
    if (q.fixedAtoms && !subsystemOnly) { frozen = frozen || new Array(N).fill(false); for (const a of q.fixedAtoms) if (a < N) frozen[a] = true; }
    const H = [meta.method, meta.basis].filter(Boolean).join('/') || null;
    const qm2Label = q.qm2Method ? q.qm2Method + (q.qm2Basis ? '/' + q.qm2Basis : '') : q.qm2 ? ({ XTB: 'GFN2-xTB', XTB0: 'GFN0-xTB', XTB1: 'GFN1-xTB', XTB2: 'GFN2-xTB', 'HF-3C': 'HF-3c', 'PBEH-3C': 'PBEh-3c', 'R2SCAN-3C': 'r2SCAN-3c', QM2: 'QM2' }[q.qm2] || q.qm2) : null;
    const mmLabel = q.mm || q.scheme === 'qmmm' ? 'MM' + (q.ff ? ' (' + baseName(q.ff) + ')' : '') : null;
    return {
      scheme: q.scheme === 'external' ? 'qmmm' : q.scheme,
      atomLayer, frozen, active, linkAtoms: links, mmTypes: null, mmCharges: null,
      labels: { H, M: qm2Label, L: mmLabel },
      nLayers: [H, qm2Label, mmLabel].filter(Boolean).length,
      embedding: q.embedding || (q.scheme === 'qm-qm2' ? null : q.scheme === 'external' ? 'point charges' : 'electronic'),
      chargeMult: q.chargeMult, subsystemOnly: !!subsystemOnly, qmAtoms: qm, qm2Atoms: q.qm2Atoms || null,
      source: subsystemOnly ? null : rebuilt ? 'rebuilt' : qm ? 'input' : null, rebuilt: rebuilt ? { qm: qm.length, links: links.length } : null
    };
  }

  function parseOrca(text, fileName) {
    const ds = newDataset('orca-out', fileName);
    const meta = ds.meta;
    const L = splitLines(text);
    const n = L.length;
    const sink = new FrameSink(ds);
    const input = [], charges = {}, spins = {};
    let cycle = 0, scan = 0, status = null, orbitals = null, eel = null, basisOut = null;
    let vib = null;                      // {frameIndex, freqs: {}, modes: [], ir: {}, raman: {}, depol: {}}
    let absRows = null, tdBlocks = [], tdCur = null, dipAU = null, dipTot = null;
    const jobsT = new JobTracker(ds);
    jobsT.open(null, null);
    let scanCur = null, cdLen = null, cdVel = null, shield = null, ccHead = -1;
    // QM/MM, QM/QM2 (multiscale): every geometry print is tagged with its step ("epoch"); after the loop each step
    // becomes one frame with the full-system geometry and all of the step's energies (QM-region prints are no steps).
    let qmmm = null, qmmmChecked = false, epoch = 0, qmmmTotalLine = false;
    const epochOf = new Map(), zMemo = new WeakMap();
    // ORCA 6 prints a step's MM energy ("FINAL SINGLE POINT ENERGY (MM)") before the step's geometry: an energy read
    // after a new step began and before its geometry waits for that geometry (it was attached to the previous step)
    let pendE = [];
    const qE = (kind, value) => {
      const f = sink.current();
      if (qmmm && (!f || (epochOf.has(f) && epochOf.get(f) !== epoch))) pendE.push([kind, value, epoch]);
      else sink.energy(kind, value);
    };
    const checkQmmm = () => {
      if (qmmmChecked) return;
      qmmmChecked = true;
      qmmm = orcaQmmmInput(input);
      if (qmmm) { sink.priority = ['QM/MM', 'QM/QM2', 'Final']; sink.splitDup = false; }
    };

    for (let i = 0; i < n; i++) {
      const raw = L[i];
      // rows starting with a digit (SCF iterations, MO coefficients, charges, gradients, A.U. coordinates) are
      // consumed by the block readers below or not needed: skip them without allocating
      let p0 = 0, c0 = raw.charCodeAt(0);
      while (c0 === 32 || c0 === 9) c0 = raw.charCodeAt(++p0);
      if (c0 >= 48 && c0 <= 57) continue;
      const t = raw.trim();
      if (!t) continue;
      const c = t[0];
      if (c === '|') {
        const m = /^\|\s*\d+>\s?(.*)$/.exec(t);
        if (m) { input.push(m[1]); if (m[1].indexOf('****END OF INPUT****') >= 0) checkQmmm(); }
        continue;
      }
      if (c === '$') {
        const m = /JOB NUMBER\s+(\d+)/.exec(t);
        if (m && int(m[1]) > 1) { jobsT.open(null, null); cycle = 0; scan = 0; scanCur = null; epoch++; }
        continue;
      }
      if (c === 'C' && (t === 'CD SPECTRUM' || t.startsWith('CD SPECTRUM VIA TRANSITION'))) {
        // ORCA 5: "State cm-1 nm R MX MY MZ"; ORCA 6: "0-1A -> 1-1A eV cm-1 nm R MX MY MZ" (R in 1e-40 cgs)
        let j = i + 2;
        while (j < n && !/^-{10,}$/.test(L[j].trim())) j++;
        const map = new Map();
        for (j++; j < n; j++) {
          const s = L[j].trim();
          if (!s || /^-{5,}$/.test(s)) break;
          const m6 = /^(\d+)-(\d+)\S*\s*->\s*(\d+)-(\d+)\S*\s+(.*)$/.exec(s);
          if (m6) { map.set(int(m6[4]) + ':' + int(m6[3]), num(tok(m6[5])[3])); continue; }   // multiplicity:state, as the states
          const v = tok(s);
          if (v.length >= 4 && isIntTok(v[0])) map.set(int(v[0]), num(v[3])); else break;
        }
        if (/VELOCITY/.test(t)) cdVel = map; else cdLen = map;
        i = j - 1;
        continue;
      }
      if (c === 'C' && t === 'CHEMICAL SHIELDING SUMMARY (ppm)') {
        let j = i + 1;
        while (j < n && j < i + 10 && !/^\s*Nucleus\s+Element/.test(L[j])) j++;
        const rows = [];
        for (j++; j < n; j++) {
          const s = L[j].trim();
          if (!s) { if (!rows.length) continue; break; }
          if (/^-+(\s+-+)*$/.test(s)) continue;
          const q = /^(\d+)\s+([A-Za-z]{1,2})\s+(\S+)\s+(\S+)/.exec(s);
          if (!q) break;
          rows.push([int(q[1]), num(q[3]), num(q[4])]);
        }
        if (rows.length) shield = rows;
        i = j - 1;
        continue;
      }
      if (t === 'CARTESIAN COORDINATES (ANGSTROEM)') {
        let j = i + 1;
        if (j < n && /^-+$/.test(L[j].trim())) j++;
        const Z = [];
        let xyz = [];
        checkQmmm();
        for (; j < n; j++) {
          if (isPH(L[j]) && L.bmSide) {   // condensed rows: labels (interned) and coordinates (null: not loaded)
            const g = L.bmSide.geoms[phIndex(L[j])];
            let zz = zMemo.get(g.codes);
            if (!zz) {
              const zOf = g.labels.map(t0 => {
                const lab = t0.replace(/[:>].*$/, ''), z = qmmm && /^H(LA|L|\*)$/i.test(lab) ? 1 : BM.el.toZ(lab);
                return z < 0 || t0.indexOf(':') > 0 ? 0 : z;
              });
              zz = Array.from(g.codes, c => zOf[c]);
              zMemo.set(g.codes, zz);
            }
            for (let a = 0; a < zz.length; a++) Z.push(zz[a]);
            xyz = g.xyz;
            continue;
          }
          let tk = tok(L[j]);
          // ghost atoms "H:" / "H :", embedding "Zr>", point charges "Q  -0.834  x y z"
          if (tk.length > 1 && (tk[1] === ':' || tk[1] === '>')) tk = [tk[0] + tk[1]].concat(tk.slice(2));
          const k = tk.length;
          if (k < 4 || k > 5 || !/^[A-Za-z]/.test(tk[0]) || !tk.slice(1).every(isNumTok)) break;
          const lab = tk[0].replace(/[:>].*$/, '');
          const z = qmmm && /^H(LA|L|\*)$/i.test(lab) ? 1 : BM.el.toZ(lab);   // QM/MM link-atom labels
          Z.push(z < 0 || tk[0].indexOf(':') > 0 ? 0 : z);
          xyz.push(num(tk[k - 3]), num(tk[k - 2]), num(tk[k - 1]));
        }
        // a block running into the end of the file is an output still being written
        if (j >= n - 1) Z.length = 0;
        if (Z.length) {
          sink.geometry(Z, xyz, 1);
          jobsT.geometry();
          const f = sink.current();
          if (qmmm) {
            epochOf.set(f, epoch);
            for (const [k, v, e] of pendE) if (e === epoch) sink.energy(k, v);
            pendE = [];
          }
          if (!f.label && (cycle || scan)) f.label = scan ? 'Scan step ' + scan + ', cycle ' + cycle : 'Cycle ' + cycle;
        }
        i = j - 1;
        continue;
      }
      if (c === '*') {
        let m = /GEOMETRY OPTIMIZATION CYCLE\s+(\d+)/.exec(t);
        if (m) { cycle = int(m[1]); epoch++; continue; }
        if (t.indexOf('FINAL ENERGY EVALUATION AT THE STATIONARY POINT') >= 0) { epoch++; continue; }
        m = /RELAXED SURFACE SCAN STEP\s+(\d+)/.exec(t);
        if (m) {
          scan = int(m[1]); cycle = 0; addJobType(meta, 'scan'); epoch++;
          // "*   Angle (  1,   0,   2)  :   100.00000000   *" (0-based atoms), one line per scanned coordinate
          const cs = [], vs = [];
          for (let j = i + 1; j < Math.min(n, i + 12); j++) {
            const q = /\b(Bond|Angle|Dihedral)\s*\(([\d,\s]+)\)\s*:\s*(-?\d+\.?\d*)/.exec(L[j]);
            if (q) cs.push({ Bond: 'R', Angle: 'A', Dihedral: 'D' }[q[1]] + '(' + q[2].split(',').map(s => int(s) + 1).join(',') + ')'), vs.push(num(q[3]));
            else if (/^\s*\*{20,}\s*$/.test(L[j])) break;
          }
          scanCur = cs.length ? { coords: cs, values: vs } : null;
          continue;
        }
        if (t.indexOf('THE OPTIMIZATION HAS CONVERGED') >= 0) {
          const f = sink.current();
          if (f) { f.converged = true; if (scanCur) scanCur.coords.forEach((c2, k) => setScan(f, scanCur.coords, k, scanCur.values[k])); }
          continue;
        }
        if (t.indexOf('IRC PATH SUMMARY') >= 0) {
          const pts = [];
          let j = i + 1;
          while (j < n && j < i + 15 && !/^\s*Step\s+E\(Eh\)/.test(L[j])) j++;
          for (j++; j < n; j++) {
            const q = /^\s*(\d+)\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)(.*)$/.exec(L[j]);
            if (!q) { if (!L[j].trim() && !pts.length) continue; break; }
            pts.push({ step: int(q[1]), energy: num(q[2]), dE: num(q[3]), maxGrad: num(q[4]), rmsGrad: num(q[5]), ts: /<=\s*TS/.test(q[6]) });
          }
          if (pts.length) { ds.ircPath = pts; addJobType(meta, 'irc'); }
          i = j - 1;
          continue;
        }
        if (t.indexOf('ORCA TERMINATED NORMALLY') >= 0) { status = true; continue; }
        continue;
      }
      if (qmmm && (t.indexOf('QM') >= 0 || t.indexOf('MM') >= 0 || t.indexOf('ink') >= 0 || t.indexOf('system size') >= 0 || t.indexOf('activeRegion') >= 0)) {
        if (orcaQmmmSummary(t, qmmm)) continue;
        const q = orcaQmmmEnergy(t);
        if (q) { qE(q.key, q.value); continue; }
      }
      if (t.startsWith('Program Version')) {
        const m = /^Program Version\s+(\d+\.\d+(?:\.\d+)?)/.exec(t);
        if (m && !ds.program) ds.program = 'ORCA ' + m[1];
        continue;
      }
      if (t.startsWith('Total Energy') && /^Total Energy\s*:\s*\S+\s+Eh/.test(t)) {
        qE(qmmm ? 'QM SCF' : 'SCF', num(/^Total Energy\s*:\s*(\S+)/.exec(t)[1]));
        continue;
      }
      if (t.startsWith('FINAL SINGLE POINT ENERGY')) {
        // ORCA 6 QM/MM prints three: "(MM)" (the full system with the force field), plain (the QM region) and
        // "(QM/MM)" (the total, the energy of the step)
        const v = num(last(tok(t))), sfx = /^FINAL SINGLE POINT ENERGY\s*\(([^)]+)\)/.exec(t);
        const q = qmmm && sfx ? orcaQmmmEnergy('E(' + sfx[1] + ') ... ' + v) : null;
        if (q) { qE(q.key, q.value); if (q.key === 'QM/MM' || q.key === 'QM/QM2') qmmmTotalLine = true; }
        else qE('Final', v);
        continue;
      }
      if (t.startsWith('MP2 TOTAL ENERGY:') || t.startsWith('RI-MP2 TOTAL ENERGY:')) {
        sink.energy('MP2', num(tok(t.split(':')[1])[0])); continue;
      }
      // ADC(2): the ground state (MP2) that the excitation energies are added to
      if (t.startsWith('Second-order corrected state energy')) { const v = num(last(tok(t))); if (Number.isFinite(v)) sink.energy('MP2', v); continue; }
      if (t.startsWith('E(CCSD')) {
        const m = /^E\((CCSD(?:\(T\))?)\)\s*\.\.\.\s*(\S+)/.exec(t);
        if (m) sink.energy(m[1], num(m[2]));
        continue;
      }
      // "COUPLED CLUSTER ENERGY" ... "E(TOT) ... -380.500397409": the ground-state CCSD (QCISD) energy, which EOM- and
      // STEOM-CCSD runs print without an "E(CCSD)" line
      if (t === 'COUPLED CLUSTER ENERGY') { ccHead = i; continue; }
      if (t.startsWith('E(TOT)') && ccHead >= 0 && i - ccHead < 10) {
        const m = /^E\(TOT\)\s*\.\.\.\s*(\S+)/.exec(t), inp = input.join('\n');
        const kind = /\bqcisd/i.test(inp) ? 'QCISD' : /ccsd/i.test(inp) ? 'CCSD' : null;
        if (m && kind) sink.energy(kind, num(m[1]));
        ccHead = -1;
        continue;
      }
      if (t.indexOf('|Geometry convergence|') >= 0) {
        const opt = {};
        const keys = { 'Energy change': 'energyChange', 'RMS gradient': 'rmsForce', 'MAX gradient': 'maxForce',
          'RMS step': 'rmsDisp', 'MAX step': 'maxDisp' };
        let j = i + 1;
        for (; j < Math.min(n, i + 12); j++) {
          const m = /^\s*(Energy change|RMS gradient|MAX gradient|RMS step|MAX step)\s+(\S+)\s+(\S+)\s+(YES|NO)/.exec(L[j]);
          if (m) { const k = keys[m[1]]; opt[k] = num(m[2]); opt[k + 'Thr'] = num(m[3]); opt[k + 'Conv'] = m[4] === 'YES'; }
          else if (/^\s*\.{5,}/.test(L[j])) break;
        }
        const f = sink.current();
        if (f) f.opt = opt;
        i = j;
        continue;
      }
      if (t.startsWith('Total Charge') || t.startsWith('Multiplicity')) {
        const m = /^(Total Charge|Multiplicity)\s+\w+\s+\.+\s+(-?\d+)/.exec(t);
        // first job of a $new_job chain, as for Gaussian Link1 files (route/method/basis are the first job's too)
        if (m) {
          if (m[1] === 'Multiplicity') { if (meta.multiplicity === null) meta.multiplicity = int(m[2]); }
          else if (meta.charge === null) meta.charge = int(m[2]);
        }
        continue;
      }
      if (t.startsWith('Your calculation utilizes the basis:')) { basisOut = t.split(':')[1].trim() || null; continue; }
      if (t === 'VIBRATIONAL FREQUENCIES') {
        vib = { frameIndex: ds.frames.length - 1, freqs: [], modes: [], ir: null, raman: null, depol: null };
        let j = i + 1;
        for (; j < n; j++) {
          const s = L[j].trim();
          const m = /^(\d+):\s+(-?\d+\.\d*)\s+cm\*\*-1/.exec(s);
          if (m) { vib.freqs[int(m[1])] = num(m[2]); continue; }
          if ((vib.freqs.length && s) || (!vib.freqs.length && j > i + 12)) break;
        }
        i = j - 1;
        continue;
      }
      if (t === 'NORMAL MODES' && vib) {
        let j = i + 1;
        const cols = vib.modes;
        let head = null;
        for (; j < n; j++) {
          if (head && isPH(L[j]) && L.bmSide) {   // the rows of this block, condensed (vecs[q] null: not loaded)
            const run = L.bmSide.modes[phIndex(L[j])];
            if (!vib.packed) vib.packed = [];
            head.forEach((k, q) => { vib.packed[k] = run.vecs[q] || null; });
            continue;
          }
          const s = L[j].trim();
          if (!s) { if (head && cols.length) { const nx = (L[j + 1] || '').trim(); if (!nx || !/^\d/.test(nx)) break; } continue; }
          const tk = s.split(/\s+/);
          if (tk.every(isIntTok)) { head = tk.map(int); for (const k of head) if (!cols[k]) cols[k] = []; continue; }
          if (head && isIntTok(tk[0]) && tk.length === head.length + 1 && tk.slice(1).every(isNumTok)) {
            const r = int(tk[0]);
            head.forEach((k, q) => { cols[k][r] = num(tk[q + 1]); });
            continue;
          }
          if (head) break;
        }
        i = j - 1;
        continue;
      }
      if (t === 'IR SPECTRUM' && vib) {
        vib.ir = [];
        let j = i + 2;
        for (; j < Math.min(n, i + 9); j++) if (/^-{10,}$/.test(L[j].trim())) break;
        for (j++; j < n; j++) {
          const m = /^\s*(\d+):\s+([^(]*)/.exec(L[j]);
          if (!m) break;
          const v = tok(m[2]).map(num);
          vib.ir[int(m[1])] = v.length >= 4 ? v[2] : v[1];
        }
        i = j - 1;
        continue;
      }
      if (t === 'RAMAN SPECTRUM' && vib) {
        vib.raman = []; vib.depol = [];
        let j = i + 2;
        for (; j < Math.min(n, i + 9); j++) if (/^-{10,}$/.test(L[j].trim())) break;
        for (j++; j < n; j++) {
          const m = /^\s*(\d+):\s+(\S+)\s+(\S+)(?:\s+(\S+))?/.exec(L[j]);
          if (!m) break;
          vib.raman[int(m[1])] = num(m[3]);
          vib.depol[int(m[1])] = m[4] !== undefined ? num(m[4]) : NaN;
        }
        i = j - 1;
        continue;
      }
      if (/EXCITED STATES/.test(t)) {
        const m = /^(TD-DFT(?:\/TDA)?|CIS|RPA|TDA|TD-DFT\/RPA)\s+EXCITED STATES(?:\s+\((SINGLETS|TRIPLETS)\))?/.exec(t);
        if (m) {
          const kind = m[2] || null;
          if (tdBlocks.some(b => b.kind === kind)) { tdBlocks = []; absRows = null; cdLen = cdVel = null; }   // a new TD-DFT run
          tdCur = { kind, states: [] }; tdBlocks.push(tdCur); addJobType(meta, 'td');
        }
        continue;
      }
      if (t.startsWith('STATE') && tdCur) {
        const m = /^STATE\s+(\d+):\s+E=\s*(\S+)\s+au\s+(\S+)\s+eV\s+(\S+)\s+cm\*\*-1(?:\s+<S\*\*2>\s*=\s*(\S+))?(?:\s+Mult\s+(\d+))?/.exec(t);
        if (m) {
          const st = { index: int(m[1]), eV: num(m[3]), cm: num(m[4]), s2: m[5] !== undefined ? num(m[5]) : null,
            mult: m[6] ? int(m[6]) : null, transitions: [] };
          tdCur.states.push(st);
          for (let j = i + 1; j < n; j++) {
            const q = /^\s*(\d+[ab]?)\s*(->|<-)\s*(\d+[ab]?)\s*:\s*(\S+)(?:\s+\(c=\s*(\S+)\))?/.exec(L[j]);
            if (!q) { i = j - 1; break; }
            const tr = { from: q[1], to: q[3], coeff: q[5] !== undefined ? num(q[5]) : Math.sqrt(Math.abs(num(q[4]))) };
            if (q[2] === '<-') tr.deexcitation = true;   // full TD-DFT (RPA) de-excitation
            st.transitions.push(tr);
          }
        }
        continue;
      }
      if (t === 'ABSORPTION SPECTRUM VIA TRANSITION ELECTRIC DIPOLE MOMENTS') {
        let j = i + 2;
        while (j < n && !/^-{10,}$/.test(L[j].trim())) j++;
        const rows = [];
        for (j++; j < n; j++) {
          const s = L[j].trim();
          if (!s || /^-{5,}$/.test(s)) break;
          const m6 = /^(\d+)-(\d+)([A-Za-z0-9'"_]*)\s*->\s*(\d+)-(\d+)([A-Za-z0-9'"_]*)\s+(.*)$/.exec(s);
          if (m6) {
            const v = tok(m6[7]).map(num);
            rows.push({ index: int(m6[4]), mult: int(m6[5]), irrep: m6[6] || '', eV: v[0], cm: v[1], nm: v[2], f: v[3] });
            continue;
          }
          const v = tok(s);
          // ORCA 5 numbers the triplets after the singlets: "6   25241.0   396.2   spin forbidden (mult=3)"
          const sf = /spin forbidden \(mult=(\d+)\)/.exec(s);
          if (sf && v.length >= 3 && isIntTok(v[0])) rows.push({ index: int(v[0]), cm: num(v[1]), nm: num(v[2]), f: 0, mult: int(sf[1]), forbidden: true });
          else if (v.length >= 4 && isIntTok(v[0])) rows.push({ index: int(v[0]), cm: num(v[1]), nm: num(v[2]), f: num(v[3]) });
          else break;
        }
        absRows = rows;
        i = j - 1;
        continue;
      }
      if (t.startsWith('MULLIKEN ATOMIC CHARGES') || t.startsWith('LOEWDIN ATOMIC CHARGES')) {
        const key = c === 'M' ? 'Mulliken' : 'Loewdin';
        const r = orcaAtomRows(L, i + 2, ORCA_CHG);
        if (r.rows.length) {
          charges[key] = r.rows.map(m => num(m[3]));
          if (r.rows[0][4] !== undefined && /SPIN/.test(t)) spins[key] = r.rows.map(m => num(m[4]));
        }
        i = r.end - 1;
        continue;
      }
      if (t === 'HIRSHFELD ANALYSIS') {
        let j = i + 1;
        while (j < n && j < i + 12 && !/^\s*ATOM\s+CHARGE/.test(L[j])) j++;
        const q = [], s = [];
        for (j++; j < n; j++) {
          const tk = tok(L[j]);
          if (tk.length < 3 || !isIntTok(tk[0])) break;
          q.push(num(tk[2])); if (tk.length > 3) s.push(num(tk[3]));
        }
        if (q.length) charges.Hirshfeld = q;
        if (s.some(v => Math.abs(v) > 1e-6)) spins.Hirshfeld = s;
        i = j - 1;
        continue;
      }
      if (t.startsWith('Total Dipole Moment')) {
        const v = floatsIn(t.split(':')[1] || '');
        if (v.length >= 3) dipAU = v.slice(0, 3);
        continue;
      }
      if (t.startsWith('Magnitude (Debye)')) { dipTot = num(t.split(':')[1]); continue; }
      if (t === 'ORBITAL ENERGIES') {
        const sets = { alpha: [], beta: [] };
        let cur = sets.alpha, j = i + 1;
        for (; j < n; j++) {
          const s = L[j].trim();
          if (/^-+$/.test(s) || /^NO\s+OCC/.test(s)) continue;
          if (/SPIN UP ORBITALS/.test(s)) { cur = sets.alpha; continue; }
          if (/SPIN DOWN ORBITALS/.test(s)) { cur = sets.beta; continue; }
          if (!s) {
            if (!cur.length) continue;
            let k = j + 1;
            while (k < n && !L[k].trim()) k++;
            if (k < n && /SPIN DOWN ORBITALS/.test(L[k])) { j = k - 1; continue; }
            break;
          }
          const tk = s.split(/\s+/);
          if (tk.length >= 3 && isIntTok(tk[0]) && isNumTok(tk[1]) && isNumTok(tk[2])) cur.push([num(tk[1]), num(tk[2])]);
          else break;
        }
        const mk = a => {
          let h = -1;
          a.forEach((p, k) => { if (p[0] > 1e-3) h = k; });
          return { energies: a.map(p => p[1]), nOcc: h + 1 };
        };
        if (sets.alpha.length) orbitals = { alpha: mk(sets.alpha), beta: sets.beta.length ? mk(sets.beta) : null };
        i = j - 1;
        continue;
      }
      if (t.startsWith('THERMOCHEMISTRY AT')) { ds.thermo = newThermo(); continue; }
      if (ds.thermo) {
        let m;
        if ((m = /^Temperature\s+\.\.\.\s+(\S+)\s+K/.exec(t))) ds.thermo.temperature = num(m[1]);
        else if ((m = /^Electronic energy\s+\.\.\.\s+(\S+)/.exec(t))) eel = num(m[1]);
        else if ((m = /^Zero point energy\s+\.\.\.\s+(\S+)/.exec(t))) ds.thermo.zpe = num(m[1]);
        else if ((m = /^Total thermal energy\s+(\S+)/.exec(t))) ds.thermo.eTotal = num(m[1]);
        else if ((m = /^Total correction\s+(\S+)/.exec(t))) ds.thermo.eCorr = num(m[1]);
        else if ((m = /^Total Enthalpy\s+\.\.\.\s+(\S+)/i.exec(t))) ds.thermo.H = num(m[1]);
        else if ((m = /^Final Gibbs free energy\s+\.\.\.\s+(\S+)/.exec(t))) ds.thermo.G = num(m[1]);
        else if ((m = /^G-E\(el\)\s+\.\.\.\s+(\S+)/.exec(t))) ds.thermo.gCorr = num(m[1]);
        if (m) continue;
      }
      if (/ORCA finished by error termination|aborting the run|ORCA TERMINATED ABNORMALLY/i.test(t)) status = false;
    }

    let qmmmSub = null, qmmmOnlySub = false, qmmmMap = null;
    if (qmmm && qmmmTotalLine) {
      for (const f of ds.frames) {
        const e = f.energies;
        if (!e || e.Final === undefined) continue;
        if (e.QM === undefined) e.QM = e.Final;
        delete e.Final;
        f.energy = null;
        for (const k of sink.priority) if (e[k] !== undefined) { f.energy = e[k]; break; }
      }
    }
    if (qmmm && ds.frames.length) {
      // one frame per step: its full-system geometry (the largest print) with the energies, convergence and labels
      // of all prints of the step; a final evaluation at the stationary point merges into the last cycle
      jobsT.close();
      const old = ds.frames;
      let big = 0;
      for (const f of old) big = Math.max(big, f.Z.length);
      const nFull = Math.max(qmmm.nInput || 0, qmmm.summary.nFull || 0, big);
      qmmmOnlySub = big < nFull;
      const target = big;
      const groups = [];
      for (const f of old) {
        const e = epochOf.has(f) ? epochOf.get(f) : -1;
        if (!groups.length || last(groups).e !== e) groups.push({ e, list: [] });
        last(groups).list.push(f);
      }
      const out = [], idx = new Map();
      const pick = f => { f.energy = null; if (f.energies) for (const k of sink.priority) if (f.energies[k] !== undefined) { f.energy = f.energies[k]; break; } };
      let skipped = 0;
      for (const gr of groups) {
        const geo = gr.list.filter(f => f.Z.length === target);
        const keep = geo.length ? geo[geo.length - 1] : null;
        if (!keep) { skipped++; for (const f of gr.list) idx.set(f, Math.max(0, out.length - 1)); continue; }
        const en = {};
        let opt = null, conv = false, label = null, sc = null;
        for (const f of gr.list) {
          if (f.energies) Object.assign(en, f.energies);
          if (f.opt) opt = f.opt;
          if (f.converged) conv = true;
          if (!label && f.label) label = f.label;
          if (f.scan) sc = f.scan;
        }
        keep.energies = Object.keys(en).length ? en : null; keep.opt = opt; keep.converged = conv; keep.label = keep.label || label;
        if (sc) keep.scan = sc;
        pick(keep);
        if (!qmmmSub && !qmmmOnlySub && keep.xyz) { const sub = gr.list.find(f => f.Z.length < target && f.xyz); if (sub) qmmmSub = { full: keep, sub }; }
        const prev = last(out);
        if (prev && sameGeom(prev, keep.Z, keep.xyz)) {
          prev.energies = Object.assign(prev.energies || {}, keep.energies || {});
          pick(prev);
          prev.converged = prev.converged || keep.converged;
          if (!prev.opt) prev.opt = keep.opt;
          for (const f of gr.list) idx.set(f, out.length - 1);
          continue;
        }
        out.push(keep);
        for (const f of gr.list) idx.set(f, out.length - 1);
      }
      for (let k = 1; k < out.length; k++) if (sameZ(out[k - 1].Z, out[k].Z)) out[k].Z = out[k - 1].Z;
      qmmmMap = k => (k >= 0 && k < old.length && idx.has(old[k]) ? idx.get(old[k]) : Math.max(0, out.length - 1));
      if (vib) vib.frameIndex = qmmmMap(vib.frameIndex);
      ds.frames = out;
      if (skipped) ds.warnings.push(skipped + ' QM/MM step(s) printed no full-system geometry and were left out.');
    }
    if (!ds.frames.length) throw err('No "CARTESIAN COORDINATES (ANGSTROEM)" block found in this ORCA output.');
    // ---- input echo -> meta
    const bang = input.filter(s => /^\s*!/.test(s)).map(s => s.replace(/^\s*!/, '').replace(/#.*$/, '').trim());
    if (bang.length) {
      meta.route = '! ' + bang.join(' ');
      meta.keywords = bang.join(' ').split(/\s+/).filter(Boolean);
    }
    const kw = (meta.keywords || []).map(s => s.toLowerCase());
    for (const k of kw) {
      if (/^(opt|copt|zopt|gdiis-opt|tightopt|verytightopt|looseopt|optts|scants|l-opt|lopt)$/.test(k)) addJobType(meta, 'opt');
      if (/^(optts|scants|neb-ts)$/.test(k)) addJobType(meta, 'ts');
      if (/^(freq|numfreq|anfreq)$/.test(k)) addJobType(meta, 'freq');
      if (k === 'irc') addJobType(meta, 'irc');
      if (k === 'nmr') addJobType(meta, 'nmr');
      if (!meta.method && METHOD_RE.test(k)) meta.method = (meta.keywords[kw.indexOf(k)]);
      if (!meta.basis && BASIS_RE.test(k) && !/\/(j|jk|c)$/i.test(k)) meta.basis = meta.keywords[kw.indexOf(k)];
      const sm = /^(?:c?pcm|smd|alpb|ddcosmo|cosmo|cpcmc)\(([^)]+)\)$/.exec(k);
      if (sm && !meta.solvent) meta.solvent = meta.keywords[kw.indexOf(k)].replace(/^[^(]*\(|\)$/g, '');
    }
    const inp = input.join('\n');
    if (/%tddft|%cis/i.test(inp)) addJobType(meta, 'td');
    if (/%geom[\s\S]*?\bscan\b/i.test(inp)) addJobType(meta, 'scan');
    if (!meta.solvent) { const m = /(?:smdsolvent|solvent)\s+"([^"]+)"/i.exec(inp); if (m) meta.solvent = m[1]; }
    if (meta.charge === null) {
      const m = /^\s*\*\s*(?:xyz|int|gzmt|xyzfile|pdbfile)\s+(-?\d+)\s+(\d+)/im.exec(inp);
      if (m) { meta.charge = int(m[1]); meta.multiplicity = int(m[2]); }
    }
    if (!meta.basis && basisOut) meta.basis = basisOut;
    finishJobTypes(meta);
    meta.normalTermination = status;
    if (!ds.program) ds.program = 'ORCA';
    if (dipAU) ds.dipole = { x: dipAU[0] * AU2DEBYE, y: dipAU[1] * AU2DEBYE, z: dipAU[2] * AU2DEBYE,
      total: dipTot !== null ? dipTot : Math.hypot(dipAU[0], dipAU[1], dipAU[2]) * AU2DEBYE };
    if (Object.keys(charges).length) ds.charges = charges;
    if (Object.keys(spins).length) ds.spins = spins;
    ds.orbitals = orbitals;
    if (ds.thermo && eel !== null) {
      if (ds.thermo.zpe !== null) ds.thermo.eZPE = eel + ds.thermo.zpe;
      if (ds.thermo.H !== null) ds.thermo.hCorr = ds.thermo.H - eel;
      if (ds.thermo.gCorr === null && ds.thermo.G !== null) ds.thermo.gCorr = ds.thermo.G - eel;
    }
    // ---- vibrations: drop the exactly-zero translation/rotation modes
    if (vib && vib.freqs.length) {
      const fi = Math.max(0, vib.frameIndex);
      const N = ds.frames[fi].Z.length;
      const keep = [];
      vib.freqs.forEach((f, k) => { if (f !== undefined && Math.abs(f) > 1e-3) keep.push(k); });
      const pick = a => (a ? keep.map(k => (a[k] === undefined ? NaN : a[k])) : null);
      const modes = vib.packed ? packedModes(keep.map(k => (vib.packed[k] === undefined ? ZERO_PACK : vib.packed[k])), N) : keep.map(k => {
        const v = new Float64Array(3 * N);
        const col = vib.modes[k];
        if (col) for (let r = 0; r < 3 * N; r++) v[r] = col[r] === undefined ? 0 : col[r];
        return v;
      });
      if (!vib.modes.length) ds.warnings.push('Frequencies were found but no NORMAL MODES block; modes cannot be animated.');
      ds.vibrations = { frameIndex: fi, freqs: pick(vib.freqs), ir: pick(vib.ir), raman: pick(vib.raman),
        depol: pick(vib.depol), redMass: null, forceConst: null, sym: null, modes };
      if (vib.packed) ds.vibrations.modesLoaded = keep.filter(k => vib.packed[k] !== null).length;
    }
    ds.excitedStates = orcaStates(absRows, tdBlocks, meta.multiplicity);
    // excited-state methods (TD-DFT, CIS, EOM-CCSD...): ORCA's "FINAL SINGLE POINT ENERGY" is then the total energy of
    // the state it follows (IROOT). It is the frame's energy only in a job that follows that state (optimisation,
    // scan, IRC); otherwise it is kept as energies['Excited state'] and the frame's energy is the ground state's, as
    // for Gaussian TD jobs. Recognised by being a ground-state energy of the frame plus an excitation energy.
    if (ds.excitedStates && ds.excitedStates.length && !meta.jobTypes.some(x => x === 'opt' || x === 'scan' || x === 'irc')) {
      const ex = ds.excitedStates.map(st => st.energyEV / HARTREE2EV).filter(Number.isFinite);
      for (const f of ds.frames) {
        const en = f.energies;
        if (!en || en.Final === undefined) continue;
        const ground = Object.keys(en).filter(k => k !== 'Final' && sink.priority.indexOf(k) >= 0).map(k => en[k]);
        if (!ground.some(g => ex.some(w => Math.abs(en.Final - g - w) < 5e-5))) continue;
        en['Excited state'] = en.Final;
        delete en.Final;
        f.energy = null;
        for (const k of sink.priority) if (en[k] !== undefined) { f.energy = en[k]; break; }
      }
    }
    if (ds.excitedStates && (cdLen || cdVel)) {
      // CD rows and states share a key: "multiplicity:state" (ORCA 6 rows) or the row number (ORCA 5)
      for (const st of ds.excitedStates) {
        const k = st.cdKey;
        if (k === null || k === undefined) continue;
        if (cdLen && Number.isFinite(cdLen.get(k))) st.R = cdLen.get(k);
        if (cdVel && Number.isFinite(cdVel.get(k))) st.Rvel = cdVel.get(k);
      }
    }
    if (ds.excitedStates) for (const st of ds.excitedStates) delete st.cdKey;
    if (shield) {
      const N = last(ds.frames).Z.length;
      const iso = new Array(N).fill(NaN), aniso = new Array(N).fill(NaN);
      for (const [k, a, b] of shield) { iso[k] = a; aniso[k] = b; }
      ds.nmr = { iso: Array.from(iso, v => (v === undefined ? NaN : v)), aniso: Array.from(aniso, v => (v === undefined ? NaN : v)), method: 'GIAO' };
      addJobType(meta, 'nmr');
    }
    // ---- jobs ($new_job): label and kinds from each job's part of the input echo
    const segs = [[]];
    for (const s of input) { if (/^\s*\$new_job\b/i.test(s)) segs.push([]); else last(segs).push(s); }
    const jl = jobsT.finish();
    if (qmmmMap) jl.forEach(jb => { jb.frameStart = qmmmMap(jb.frameStart); jb.frameEnd = qmmmMap(jb.frameEnd); });
    jl.forEach((jb, k) => {
      const seg = segs[k] || [];
      const bang = seg.filter(s => /^\s*!/.test(s)).map(s => s.replace(/^\s*!/, '').replace(/#.*$/, '').trim()).join(' ').trim();
      jb.label = bang ? '! ' + bang : null;
      jb.kinds = jobKinds(seg.length ? orcaJobTypes(bang, seg.join('\n')) : meta.jobTypes);
      if (ds.frames.slice(jb.frameStart, jb.frameEnd + 1).some(f => f.scan) && jb.kinds.indexOf('scan') < 0) jb.kinds.push('scan');
      if (k === jl.length - 1 && ds.ircPath && jb.kinds.indexOf('irc') < 0) jb.kinds = jb.kinds.filter(x => x !== 'sp').concat(['irc']);
      if (k === jl.length - 1 && shield && jb.kinds.indexOf('nmr') < 0) jb.kinds = jb.kinds.filter(x => x !== 'sp').concat(['nmr']);
      jb.normalTermination = k < jl.length - 1 ? true : status;
    });
    ds.jobs = jl;
    if (qmmm) {
      ds.layers = orcaLayers(qmmm, qmmmSub ? qmmmSub.full : ds.frames.find(f => f.xyz) || ds.frames[0], qmmmSub ? qmmmSub.sub : null, meta, qmmmOnlySub);
      meta.qmmm = true;
      if (qmmmOnlySub) ds.warnings.push('This output prints only the QM region (with link atoms); ORCA writes the full-system geometries to its trajectory and PDB files.');
      else if (ds.layers.source === 'rebuilt') ds.warnings.push('The QM atoms are not listed in this output (no QMAtoms in the input echo). BlueMolecule rebuilt them from the ' +
        ds.layers.rebuilt.links + ' link atom' + (ds.layers.rebuilt.links === 1 ? '' : 's') + ': the ' + ds.layers.rebuilt.qm + ' atoms bonded to their QM side are as many as the QM atoms ORCA counted.');
      else if (!ds.layers.atomLayer) ds.warnings.push('The QM atoms are not listed in this output (no QMAtoms in the input echo, and no QM-region geometry to match), so the layers are not shown. ' +
        'Open the input or another output of the same system that lists them, then take its layers in the Info tab (QM/MM layers).');
    }
    return ds;
  }
  function orcaJobTypes(kwStr, inp) {
    const t = [];
    for (const k of kwStr.toLowerCase().split(/\s+/)) {
      if (/^(opt|copt|zopt|gdiis-opt|tightopt|verytightopt|looseopt|optts|scants|l-opt|lopt)$/.test(k)) t.push('opt');
      if (/^(freq|numfreq|anfreq)$/.test(k)) t.push('freq');
      if (k === 'irc') t.push('irc');
      if (k === 'nmr') t.push('nmr');
    }
    if (/%tddft|%cis\b/i.test(inp)) t.push('td');
    if (/%geom[\s\S]*?\bscan\b/i.test(inp)) t.push('scan');
    return t;
  }

  function orcaStateLabel(kind, st, row, groundMult) {
    if (row && row.mult) return multName(row.mult) + (row.irrep ? '-' + row.irrep : '');
    if (st && st.mult) return multName(st.mult);
    if (kind === 'SINGLETS') return 'Singlet';
    if (kind === 'TRIPLETS') return 'Triplet';
    if ((groundMult || 1) > 1 && st && st.s2 !== null) {
      const S = (-1 + Math.sqrt(1 + 4 * Math.max(0, st.s2))) / 2;
      return multName(Math.round(2 * S + 1));
    }
    return multName(groundMult || 1);
  }

  function orcaStates(absRows, blocks, groundMult) {
    if (!blocks.length && !absRows) return null;
    const main = blocks.find(b => b.kind === 'SINGLETS' || b.kind === null) || null;
    const ofMult = m => (m === 3 ? blocks.find(b => b.kind === 'TRIPLETS') : m === 1 ? main : null) || null;
    const out = [];
    const used = new Set();
    if (absRows && absRows.length) {
      for (const r of absRows) {
        const ev = Number.isFinite(r.eV) ? r.eV : r.cm * CM2EV;
        // the row's state in an excited-states block: by its number within its multiplicity (ORCA 6: "0-1A -> 2-3B"),
        // or by energy for ORCA 5's spin-forbidden rows, which are numbered after the singlets
        const blk = r.mult ? ofMult(r.mult) : main;
        let st = null;
        if (blk) {
          if (r.mult) {   // (ORCA 6 numbers a triplet block's states 6, 7, ... after the singlets, its table 1-3, 2-3, ...)
            const near = blk.states.filter(x => !used.has(x) && Math.abs(x.eV - ev) < 2e-3);
            st = near.find(x => x.index === r.index) || near[0] || null;
          } else st = blk.states.find(x => x.index === r.index) || null;
          if (st) used.add(st);
        }
        out.push({ index: r.forbidden && st ? st.index : r.index, energyEV: ev, wavelengthNM: r.cm > 0 ? 1e7 / r.cm : EV2NM / ev, f: r.f,
          label: orcaStateLabel(blk ? blk.kind : null, st, r, groundMult), s2: st ? st.s2 : null,
          transitions: st ? st.transitions : [], cdKey: r.mult && !r.forbidden ? r.mult + ':' + r.index : r.index });
      }
    }
    for (const b of blocks) {
      for (const st of b.states) {
        if (used.has(st) || (b === main && absRows && absRows.length)) continue;
        out.push({ index: st.index, energyEV: st.eV, wavelengthNM: st.eV > 0 ? EV2NM / st.eV : NaN,
          f: b.kind === 'TRIPLETS' ? 0 : null, label: orcaStateLabel(b.kind, st, null, groundMult), s2: st.s2,
          transitions: st.transitions, cdKey: b === main ? st.index : null });
      }
    }
    return out.length ? out : null;
  }

  // ------------------------------------------------------------ ORCA input
  /** ORCA atom line "O 0 0 0", "O(1) ...", "H : ..." (ghost), "DA ..." (dummy); null for point charges "Q". */
  function orcaInpAtom(t) {
    let tk = tok(t);
    if (tk.length > 1 && tk[1] === ':') tk = [tk[0] + ':'].concat(tk.slice(2));
    if (!tk.length || /^q$/i.test(tk[0])) return null;
    const lab = tk[0].replace(/[({:].*$/, '');
    let z = /^(da|x)$/i.test(lab) ? 0 : BM.el.toZ(lab);
    if (z < 0) return undefined;
    if (tk[0].indexOf(':') > 0) z = 0;
    return { z, v: tk.slice(1) };
  }
  function looksLikeOrcaInp(head) {
    if (/&GLOBAL\b|&FORCE_EVAL\b/i.test(head) || /^\s*\$(CONTRL|DATA)\b/im.test(head)) return false;
    return /^\s*!\s*\S/m.test(head) &&
      (/^\s*\*\s*(xyz|int|internal|gzmt|xyzfile|gzmtfile|pdbfile)\s+-?\d+\s+\d+/im.test(head) || /^\s*%coords\b/im.test(head));
  }
  /** ORCA %paras: NAME = value | start, end, steps | [v1 v2 ...] -> {name (lower case): first value}; scans lists
   *  the parameters that take several values (the geometry shown uses the first). */
  function orcaParas(L) {
    const vals = {}, scans = [];
    for (let i = 0; i < L.length; i++) {
      if (!/^\s*%paras\b/i.test(L[i])) continue;
      for (let j = i + 1; j < L.length && !/^\s*end\s*$/i.test(L[j]); j++) {
        const m = /^\s*(\w+)\s*=\s*(.+?)\s*$/.exec(L[j]);
        if (!m) continue;
        const parts = m[2].replace(/[[\]]/g, ' ').split(/[\s,]+/).filter(Boolean);
        if (!parts.length || !isNumTok(parts[0])) continue;
        vals[m[1].toLowerCase()] = num(parts[0]);
        if (parts.length > 1) scans.push(m[1]);
        i = j;
      }
    }
    return { vals, scans };
  }
  function parseOrcaInput(text, fileName) {
    const ds = newDataset('orca-input', fileName);
    const meta = ds.meta;
    const paras = orcaParas(splitLines(text).map(s => s.replace(/#.*$/, '')));
    const unknown = new Set();
    // coordinates may name %paras values in braces: "C 1 0 0 {B1} 0 0"
    const sub = s => s.replace(/\{\s*(\w+)\s*\}/g, (all, k) => {
      const v = paras.vals[k.toLowerCase()];
      if (v === undefined) { unknown.add(k); return all; }
      return String(v);
    });
    const L = splitLines(text).map(s => s.replace(/#.*$/, '')).map(s => (s.indexOf('{') >= 0 ? sub(s) : s));
    const kw = [];
    let job = 0, bohr = false;
    const addGeom = (Z, xyz, label) => {
      if (!Z.length) return;
      if (bohr) for (let k = 0; k < xyz.length; k++) xyz[k] *= BOHR;
      const f = newFrame(Z, xyz);
      f.label = label;
      ds.frames.push(f);
    };
    const setCM = (c, m) => { if (meta.charge === null) { meta.charge = c; meta.multiplicity = m; } };
    for (let i = 0; i < L.length; i++) {
      const t = L[i].trim();
      if (!t) continue;
      if (/^\$new_job\b/i.test(t)) { job++; continue; }
      if (t[0] === '!') { if (!job) kw.push(...t.slice(1).trim().split(/\s+/).filter(Boolean)); continue; }
      if (t[0] === '%') {
        const name = t.slice(1).split(/\s+/)[0].toLowerCase();
        if (name === 'coords') {
          let typ = 'xyz', ch = 0, mu = 1, units = 'angs';
          const Z = [], xyz = [];
          let j = i + 1;
          for (; j < L.length; j++) {
            const s = L[j].trim();
            if (/^end$/i.test(s)) break;
            const kv = /^(ctyp|charge|mult|units)\s+(\S+)/i.exec(s);
            if (kv) {
              const k = kv[1].toLowerCase();
              if (k === 'ctyp') typ = kv[2].toLowerCase(); else if (k === 'charge') ch = int(kv[2]);
              else if (k === 'mult') mu = int(kv[2]); else units = kv[2].toLowerCase();
              continue;
            }
            if (/^coords$/i.test(s)) {
              for (j++; j < L.length && !/^end$/i.test(L[j].trim()); j++) {
                const a = orcaInpAtom(L[j]);
                if (a && a.v.length >= 3 && a.v.slice(0, 3).every(isNumTok)) { Z.push(a.z); xyz.push(num(a.v[0]), num(a.v[1]), num(a.v[2])); }
              }
            }
          }
          i = j;
          setCM(ch, mu);
          if (typ !== 'xyz') ds.warnings.push('%coords with CTyp ' + typ + ' is not read; only Cartesian %coords blocks are.');
          else {
            if (/^bohr/.test(units)) for (let k = 0; k < xyz.length; k++) xyz[k] *= BOHR;
            const f = newFrame(Z, xyz); f.label = job ? 'Job ' + (job + 1) : null;
            if (Z.length) ds.frames.push(f);
          }
          continue;
        }
        // skip "%block ... end" (one line, or up to an "end" line; a new directive also ends it, e.g. %maxcore 2000)
        if (/\bend$/i.test(t) || /^%(maxcore|moinp|base|id)\b/i.test(t)) continue;
        let j = i + 1;
        for (; j < L.length; j++) {
          const s = L[j].trim();
          if (/^end$/i.test(s)) break;
          if (/^[!%*$]/.test(s)) { j--; break; }
        }
        i = j;
        continue;
      }
      if (t[0] === '*') {
        const m = /^\*\s*(xyz|int|internal|gzmt|xyzfile|gzmtfile|pdbfile)\s+(-?\d+)\s+(\d+)(?:\s+(\S+))?/i.exec(t);
        if (!m) continue;
        const typ = m[1].toLowerCase();
        setCM(int(m[2]), int(m[3]));
        bohr = kw.some(k => /^bohrs?$/i.test(k));
        const label = job ? 'Job ' + (job + 1) : null;
        if (/file$/.test(typ)) {
          ds.warnings.push('Coordinates of ' + (job ? 'job ' + (job + 1) : 'this job') + ' are read from "' + (m[4] || '?') +
            '" (* ' + typ + '); open that file to see them.');
          continue;
        }
        const rows = [], Z = [], xyz = [];
        let j = i + 1;
        for (; j < L.length && !/^\s*\*\s*$/.test(L[j]); j++) {
          const s = L[j].trim();
          if (!s) continue;
          const a = orcaInpAtom(s);
          if (a === null) continue;
          if (a === undefined) break;
          if (typ === 'xyz') {
            if (a.v.length < 3 || !a.v.slice(0, 3).every(isNumTok)) break;
            Z.push(a.z); xyz.push(num(a.v[0]), num(a.v[1]), num(a.v[2]));
          } else if (typ === 'gzmt') {   // Gaussian-style: Sym [a r [b angle [c dihedral]]], 1-based references
            const v = a.v;
            const r = {};
            if (v.length >= 2) { r.a = int(v[0]) - 1; r.r = num(v[1]); }
            if (v.length >= 4) { r.b = int(v[2]) - 1; r.angle = num(v[3]); }
            if (v.length >= 6) { r.c = int(v[4]) - 1; r.dihedral = num(v[5]); }
            rows.push(r); Z.push(a.z);
          } else {                          // * int: Sym NA NB NC R A D (1-based references, 0 = none)
            const v = a.v.map(num);
            if (v.length < 6) break;
            const r = {};
            if (rows.length >= 1) { r.a = v[0] - 1; r.r = v[3]; }
            if (rows.length >= 2) { r.b = v[1] - 1; r.angle = v[4]; }
            if (rows.length >= 3) { r.c = v[2] - 1; r.dihedral = v[5]; }
            rows.push(r); Z.push(a.z);
          }
        }
        i = j;
        if (typ === 'xyz') addGeom(Z, xyz, label);
        else if (Z.length) addGeom(Z, Array.from(zmatrixToCartesian(rows)), label);   // lengths scale with the units
        continue;
      }
    }
    if (unknown.size) {
      throw err('The coordinates use ' + Array.from(unknown, k => '{' + k + '}').join(', ') + ', which the input does not define in %paras.');
    }
    if (!ds.frames.length) {
      throw err(ds.warnings.length ? ds.warnings.join(' ') : 'No "* xyz", "* int" or %coords geometry found in this ORCA input.');
    }
    if (paras.scans.length) {
      ds.warnings.push('This input scans ' + paras.scans.join(', ') + ' (%paras); the geometry shown uses the first value' +
        (paras.scans.length > 1 ? 's.' : '.'));
      addJobType(meta, 'scan');
    }
    if (kw.length) { meta.route = '! ' + kw.join(' '); meta.keywords = kw.slice(); }
    const lo = kw.map(s => s.toLowerCase());
    lo.forEach((k, q) => {
      if (!meta.method && METHOD_RE.test(k)) meta.method = kw[q];
      if (!meta.basis && BASIS_RE.test(k) && !/\/(j|jk|c)$/i.test(k)) meta.basis = kw[q];
    });
    for (const t of orcaJobTypes(kw.join(' '), L.join('\n'))) addJobType(meta, t);
    if (lo.some(k => /^(optts|scants|neb-ts)$/.test(k))) addJobType(meta, 'ts');
    finishJobTypes(meta);
    return ds;
  }

  // ------------------------------------------------------------------- XYZ
  function xyzCommentEnergy(c) {
    if (!c) return null;
    let m;
    // extended XYZ (ASE convention): energy= is in eV
    if (/(?:^|\s)(?:Properties|Lattice)\s*=/i.test(c)) {
      m = /(?:^|\s)energy\s*=\s*"?([-+]?(?:\d+\.?\d*|\.\d+)(?:[eEdD][-+]?\d+)?)/i.exec(c);
      return m ? num(m[1]) / HARTREE2EV : null;
    }
    m = /(?:^|[\s,;(])(?:e|energy|etot|e_tot|total_energy)\s*[=:]\s*([-+]?(?:\d+\.?\d*|\.\d+)(?:[eEdD][-+]?\d+)?)/i.exec(c);
    if (m) return num(m[1]);
    m = /(?:^|\s)E\s+([-+]?\d+\.\d+(?:[eEdD][-+]?\d+)?)(?=\s|$)/.exec(c);
    if (m) return num(m[1]);
    m = /^\s*([-+]?\d+\.\d+(?:[eEdD][-+]?\d+)?)\s*$/.exec(c);
    return m ? num(m[1]) : null;
  }

  function xyzAtom(t, ds, where) {
    let z = BM.el.toZ(t[0]);
    if (z < 0) { ds.warnings.push('Unknown element "' + t[0] + '" ' + where + ' read as a dummy atom.'); z = 0; }
    return z;
  }

  function parseXYZ(text, fileName) {
    const ds = newDataset('xyz', fileName);
    const L = splitLines(text);
    const n = L.length;
    const zMemo = new WeakMap();
    let i = 0;
    while (i < n && !L[i].trim()) i++;
    const t0 = tok(L[i] || '');
    if (t0.length >= 1 && isIntTok(t0[0]) && (t0.length === 1 || !isNumTok(t0[1]))) {
      while (i < n) {
        while (i < n && !L[i].trim()) i++;
        if (i >= n) break;
        const h = tok(L[i]);
        if (!isIntTok(h[0])) {
          if (ds.frames.length) { ds.warnings.push('Stopped reading at line ' + (i + 1) + ' (expected an atom count).'); break; }
          throw err('Line ' + (i + 1) + ' should contain the number of atoms.');
        }
        const na = int(h[0]);
        const comment = i + 1 < n ? L[i + 1] : '';
        let f;
        if (na > 0 && isPH(L[i + 2]) && L.bmSide) {   // condensed frame (condense.js): labels and coordinates
          const g = L.bmSide.geoms[phIndex(L[i + 2])];
          let zz = zMemo.get(g.codes);
          if (!zz) {
            const zOf = g.labels.map(lb => xyzAtom([lb], ds, 'in frame ' + (ds.frames.length + 1)));
            zz = Array.from(g.codes, c => zOf[c]);
            zMemo.set(g.codes, zz);
          }
          f = newFrame([], g.xyz);
          f.Z = zz;   // frames with the same atoms share one element list
          i -= na - 1;   // the placeholder stands for the na atom lines
        } else {
          const Z = [], xyz = [];
          let ok = true;
          for (let k = 0; k < na; k++) {
            const t = tok(L[i + 2 + k] || '');
            if (t.length < 4 || !isNumTok(t[1]) || !isNumTok(t[2]) || !isNumTok(t[3])) { ok = false; break; }
            Z.push(xyzAtom(t, ds, 'on line ' + (i + 3 + k)));
            xyz.push(num(t[1]), num(t[2]), num(t[3]));
          }
          if (!ok) {
            if (ds.frames.length) { ds.warnings.push('Frame ' + (ds.frames.length + 1) + ' is incomplete and was skipped.'); break; }
            throw err('Expected ' + na + ' atom lines ("Element x y z") after line ' + (i + 2) + '.');
          }
          f = newFrame(Z, xyz);
        }
        f.energy = xyzCommentEnergy(comment);
        f.label = comment.trim() || null;
        ds.frames.push(f);
        if (ds.frames.length === 1) {
          const lat = /Lattice\s*=\s*"([^"]+)"/i.exec(comment);
          if (lat) { const v = floatsIn(lat[1]); if (v.length >= 9) ds.cell = [v.slice(0, 3), v.slice(3, 6), v.slice(6, 9)]; }
        }
        i += 2 + na;
      }
    } else {
      const Z = [], xyz = [];
      for (; i < n; i++) {
        const t = tok(L[i]);
        if (!t.length) continue;
        if (t.length < 4 || !isNumTok(t[1]) || !isNumTok(t[2]) || !isNumTok(t[3]) || BM.el.toZ(t[0]) < 0) {
          if (Z.length) break;
          throw err('Not an XYZ file: line ' + (i + 1) + ' is not "Element x y z".');
        }
        Z.push(xyzAtom(t, ds, 'on line ' + (i + 1)));
        xyz.push(num(t[1]), num(t[2]), num(t[3]));
      }
      ds.frames.push(newFrame(Z, xyz));
    }
    if (!ds.frames.length) throw err('No atoms found.');
    if (ds.frames.length === 1 && ds.frames[0].label) ds.title = ds.frames[0].label;
    return ds;
  }

  // ------------------------------------------------------------------ cube
  const POW10 = [];
  for (let k = 0; k <= 340; k++) POW10.push(Number('1e' + k));

  function parseCube(text, fileName) {
    const ds = newDataset('cube', fileName);
    // text, or a line array (big files): the header is read line by line, the values by the scanner below
    const arr = Array.isArray(text) ? text : null;
    const N = arr ? arr.length : text.length;
    let pos = 0;
    const line = arr ? () => (pos < N ? arr[pos++] : '') : () => {
      if (pos >= N) return '';
      let e = text.indexOf('\n', pos);
      if (e < 0) e = N;
      const s = text.slice(pos, e);
      pos = e + 1;
      return s;
    };
    const c1 = line(), c2 = line();
    let t = tok(line());
    if (t.length < 4 || !isIntTok(t[0])) throw err('Line 3 of a cube file must hold the atom count and the grid origin.');
    let nat = int(t[0]);
    const hasMO = nat < 0;
    nat = Math.abs(nat);
    let nval = t.length >= 5 && isIntTok(t[4]) ? Math.max(1, int(t[4])) : 1;
    const o = [num(t[1]), num(t[2]), num(t[3])];
    const dims = [], axes = [];
    let ang = false;
    for (let a = 0; a < 3; a++) {
      t = tok(line());
      if (t.length < 4 || !isIntTok(t[0])) throw err('Cube grid axis line ' + (a + 1) + ' is malformed.');
      const nn = int(t[0]);
      if (a === 0) ang = nn < 0;
      dims.push(Math.abs(nn));
      const f = nn < 0 ? 1 : BOHR;
      axes.push([num(t[1]) * f, num(t[2]) * f, num(t[3]) * f]);
    }
    const fo = ang ? 1 : BOHR;
    const Z = [], xyz = [];
    for (let k = 0; k < nat; k++) {
      t = tok(line());
      if (t.length < 5) throw err('Cube atom line ' + (k + 1) + ' is incomplete.');
      Z.push(Math.max(0, int(t[0])));
      xyz.push(num(t[2]) * fo, num(t[3]) * fo, num(t[4]) * fo);
    }
    let ids = null;
    if (hasMO) {
      ids = [];
      let count = -1;
      while (pos < N && (count < 0 || ids.length < count)) {
        for (const x of tok(line())) { if (count < 0) count = int(x); else ids.push(int(x)); }
      }
      nval = Math.max(1, count);
      ids = ids.slice(0, nval);
    }
    const npts = dims[0] * dims[1] * dims[2];
    if (!(npts > 0)) throw err('Cube grid dimensions are invalid.');
    const total = npts * nval;
    const vols = [];
    for (let v = 0; v < nval; v++) vols.push(new Float32Array(npts));
    // Hand-written float scanner over the raw text (fast for 10^6+ values); over each line for a line array.
    let cnt = 0;
    const single = nval === 1 ? vols[0] : null;
    const scan = (text, p, N) => {
      while (cnt < total) {
        let ch = text.charCodeAt(p);
        while (ch === 32 || ch === 10 || ch === 13 || ch === 9 || ch === 44) ch = text.charCodeAt(++p);
        if (p >= N) break;
        let neg = false;
        if (ch === 45) { neg = true; ch = text.charCodeAt(++p); } else if (ch === 43) ch = text.charCodeAt(++p);
        let mant = 0, ex = 0, nd = 0;
        while (ch >= 48 && ch <= 57) { mant = mant * 10 + (ch - 48); nd++; ch = text.charCodeAt(++p); }
        if (ch === 46) {
          ch = text.charCodeAt(++p);
          while (ch >= 48 && ch <= 57) { mant = mant * 10 + (ch - 48); ex--; nd++; ch = text.charCodeAt(++p); }
        }
        let val;
        if (nd === 0) {
          while (p < N && ch !== 32 && ch !== 10 && ch !== 13 && ch !== 9) ch = text.charCodeAt(++p);
          val = NaN;
        } else {
          if (ch === 69 || ch === 101 || ch === 68 || ch === 100) {
            ch = text.charCodeAt(++p);
            let en = false, e = 0;
            if (ch === 45) { en = true; ch = text.charCodeAt(++p); } else if (ch === 43) ch = text.charCodeAt(++p);
            while (ch >= 48 && ch <= 57) { e = e * 10 + (ch - 48); ch = text.charCodeAt(++p); }
            ex += en ? -e : e;
          }
          if (ex < 0) val = -ex <= 340 ? mant / POW10[-ex] : 0;
          else val = ex <= 308 ? mant * POW10[ex] : Infinity;
          if (neg) val = -val;
        }
        if (single) single[cnt] = val;
        else vols[cnt % nval][(cnt / nval) | 0] = val;
        cnt++;
      }
    };
    if (arr) { for (let k = pos; k < N && cnt < total; k++) scan(arr[k], 0, arr[k].length); }
    else scan(text, pos, N);
    if (cnt < total) throw err('Cube data are truncated: expected ' + total + ' values, found ' + cnt + '.');
    const head = (c1 + ' ' + c2).toLowerCase();
    let kind = 'other';
    if (hasMO || /\bmo\b|orbital/.test(head)) kind = 'orbital';
    else if (/potential|\besp\b/.test(head)) kind = 'potential';   // ("Electrostatic potential from Total SCF Density")
    else if (/density|rho/.test(head)) kind = 'density';
    const base = c1.trim() || stripExt(fileName) || 'volume';
    ds.volumes = vols.map((data, v) => ({
      name: ids ? 'MO ' + ids[v] : (nval > 1 ? base + ' #' + (v + 1) : base),
      kind, origin: [o[0] * fo, o[1] * fo, o[2] * fo], axes: axes.map(a => a.slice()), dims: dims.slice(), data
    }));
    if (ids) ds.volumes.forEach((vol, v) => { vol.mo = ids[v]; });
    ds.title = c1.trim() || ds.title;
    const f = newFrame(Z, xyz);
    ds.frames.push(f);
    ds.program = null;
    return ds;
  }

  // ---------------------------------------------------------------- Molden
  function parseMolden(text, fileName) {
    const ds = newDataset('molden', fileName);
    const sections = [];
    let cur = null;
    for (const line of splitLines(text)) {
      const m = /^\s*\[([^\]]+)\]\s*(.*)$/.exec(line);
      if (m) { cur = { name: m[1].trim().toUpperCase(), arg: m[2].trim(), lines: [] }; sections.push(cur); }
      else if (cur) cur.lines.push(line);
    }
    if (!sections.length) throw err('No [Section] headers found; this is not a Molden file.');
    const sec = name => sections.find(s => s.name === name) || null;
    const title = sec('TITLE');
    if (title) { const tl = title.lines.map(s => s.trim()).filter(Boolean); if (tl.length) ds.title = tl.join(' '); }
    const atoms = sec('ATOMS');
    if (atoms) {
      const f = /au|bohr/i.test(atoms.arg) ? BOHR : 1;
      const Z = [], xyz = [];
      for (const s of atoms.lines) {
        const t = tok(s);
        if (t.length < 6) continue;
        let z = isIntTok(t[2]) ? int(t[2]) : BM.el.toZ(t[0]);
        if (z < 0) z = 0;
        Z.push(z); xyz.push(num(t[3]) * f, num(t[4]) * f, num(t[5]) * f);
      }
      if (Z.length) ds.frames.push(newFrame(Z, xyz));
    }
    const geoms = sec('GEOMETRIES');
    if (geoms) {
      try {
        const g = parseXYZ(geoms.lines.join('\n'), fileName);
        ds.frames = g.frames;
        ds.frames.forEach((fr, k) => { fr.energy = null; fr.label = 'Step ' + (k + 1); });
      } catch (e) { ds.warnings.push('Could not read [GEOMETRIES]: ' + e.message); }
    }
    const conv = sec('GEOCONV');
    if (conv) {
      const series = {};
      let key = null;
      for (const s of conv.lines) {
        const t = s.trim();
        if (!t) continue;
        if (isNumTok(t.split(/\s+/)[0])) { if (key) series[key].push(num(t.split(/\s+/)[0])); }
        else { key = t.toLowerCase(); series[key] = []; }
      }
      const en = series.energy || [];
      const off = Math.max(0, ds.frames.length - en.length);
      en.forEach((e, k) => { const fr = ds.frames[k + off]; if (fr) { fr.energy = e; fr.energies = { SCF: e }; } });
      ds.frames.forEach((fr, k) => {
        const o = {};
        [['max-force', 'maxForce'], ['rms-force', 'rmsForce'], ['max-step', 'maxDisp'], ['rms-step', 'rmsDisp']]
          .forEach(([a, b]) => { if (series[a] && series[a][k - off] !== undefined) o[b] = series[a][k - off]; });
        if (Object.keys(o).length) fr.opt = o;
      });
    }
    const freq = sec('FREQ'), frc = sec('FR-COORD'), nrm = sec('FR-NORM-COORD');
    if (freq && nrm) {
      const freqs = freq.lines.map(s => tok(s)).filter(t => t.length && isNumTok(t[0])).map(t => num(t[0]));
      let fi = ds.frames.length - 1;
      if (frc) {
        const Z = [], xyz = [];
        for (const s of frc.lines) {
          const t = tok(s);
          if (t.length < 4) continue;
          const z = BM.el.toZ(t[0]);
          Z.push(z < 0 ? 0 : z); xyz.push(num(t[1]) * BOHR, num(t[2]) * BOHR, num(t[3]) * BOHR);
        }
        if (Z.length && !(fi >= 0 && sameGeom(ds.frames[fi], Z, xyz, 2e-3))) {
          const f = newFrame(Z, xyz); f.label = 'Vibrational analysis geometry';
          ds.frames.push(f); fi = ds.frames.length - 1;
        }
      }
      if (fi < 0) throw err('Molden file has frequencies but no geometry.');
      const N = ds.frames[fi].Z.length;
      const modes = [];
      let curm = null;
      for (const s of nrm.lines) {
        const t = tok(s);
        if (!t.length) continue;
        if (/^vibration$/i.test(t[0])) { curm = []; modes.push(curm); continue; }
        if (curm && t.length >= 3) curm.push(num(t[0]), num(t[1]), num(t[2]));
      }
      if (modes.length !== freqs.length) ds.warnings.push('Molden file lists ' + freqs.length + ' frequencies but ' + modes.length + ' normal modes.');
      const nm = Math.min(modes.length, freqs.length);
      const ints = sec('INT');
      let ir = null, raman = null;
      if (ints) {
        const rows = ints.lines.map(s => tok(s)).filter(t => t.length && isNumTok(t[0]));
        ir = rows.slice(0, nm).map(t => num(t[0]));
        if (rows.length && rows.every(t => t.length >= 2 && isNumTok(t[1]))) raman = rows.slice(0, nm).map(t => num(t[1]));
        if (ir.length !== nm) ir = null;
        if (raman && raman.length !== nm) raman = null;
      }
      const mv = [];
      for (let k = 0; k < nm; k++) {
        const v = new Float64Array(3 * N);
        if (modes[k].length !== 3 * N) ds.warnings.push('Normal mode ' + (k + 1) + ' does not match the atom count.');
        for (let q = 0; q < Math.min(3 * N, modes[k].length); q++) v[q] = modes[k][q];
        mv.push(v);
      }
      ds.vibrations = { frameIndex: fi, freqs: freqs.slice(0, nm), ir, raman, depol: null, redMass: null,
        forceConst: null, sym: null, modes: mv };
    }
    if (sec('MO')) moldenOrbitals(ds, sections, sec);
    if (!ds.frames.length) throw err('Molden file contains no [Atoms], [GEOMETRIES] or [FR-COORD] geometry.');
    return ds;
  }

  /** The [GTO] basis and the [MO] orbitals of a Molden file as ds.basis and ds.mos (the layout of BM.orbitals:
   *  sp shells split into s + p, centres in Bohr, MO coefficients one orbital after the other), and the orbital
   *  energies as ds.orbitals. [5D] and [5D7F] make d and f pure, [7F] f, [5D10F] d, [9G] g and higher. A set whose
   *  energies are not in ascending order (symmetry blocks: Molpro, Molcas) is sorted, with each orbital's place in
   *  the file kept as fileIndex. */
  function moldenOrbitals(ds, sections, sec) {
    const bad = why => { ds.warnings.push('The orbitals of this Molden file cannot be computed: ' + why + '.'); };
    const gto = sec('GTO'), mo = sec('MO'), atoms = sec('ATOMS');
    if (!gto) return bad(sec('STO') ? 'Slater-type orbitals ([STO]) are not supported' : 'it has no [GTO] basis set');
    if (!atoms) return bad('it has no [Atoms] section');
    const fa = /au|bohr/i.test(atoms.arg) ? 1 : 1 / BOHR;
    const centre = new Map();
    let na = 0;
    for (const s of atoms.lines) {
      const t = tok(s);
      if (t.length < 6) continue;
      const id = isIntTok(t[1]) ? int(t[1]) : na + 1;
      if (!centre.has(id)) centre.set(id, { index: na, xyz: [num(t[3]) * fa, num(t[4]) * fa, num(t[5]) * fa] });
      na++;
    }
    const has = n => sections.some(s => s.name === n);
    const pure = [false, false, has('5D') || has('5D7F') || has('5D10F'), has('5D') || has('5D7F') || has('7F')];
    const pureHigh = has('9G');
    const shells = [];
    const L = gto.lines;
    let atom = null, scaled = false;
    for (let k = 0; k < L.length; k++) {
      const t = tok(L[k]);
      if (!t.length) continue;
      if (isIntTok(t[0]) && (t.length === 1 || (t.length === 2 && isIntTok(t[1])))) {   // "n 0": the atom of the shells below
        atom = centre.get(int(t[0]));
        if (!atom) return bad('[GTO] names atom ' + t[0] + ', which [Atoms] does not list');
        continue;
      }
      const lab = t[0].toLowerCase(), sp = lab === 'sp' || lab === 'l';
      const l = sp ? 0 : 'spdfghik'.indexOf(lab), np = t.length > 1 && isIntTok(t[1]) ? int(t[1]) : -1;
      if (lab.length > 2 || l < 0 || np < 1) return bad('line ' + (k + 1) + ' of [GTO] is not a shell ("' + L[k].trim().slice(0, 40) + '")');
      if (!atom) return bad('a shell comes before its atom in [GTO]');
      // the third number scales the exponents (by its square, as in Gaussian basis input); it is almost always 1
      const sc = t.length > 2 && isNumTok(t[2]) ? num(t[2]) : 1, f2 = sc > 0 && Math.abs(sc - 1) > 1e-12 ? sc * sc : 1;
      if (f2 !== 1) scaled = true;
      const exps = [], c1 = [], c2 = [];
      for (let q = 1; q <= np; q++) {
        const u = tok(L[k + q] || '');
        if (u.length < (sp ? 3 : 2) || !isNumTok(u[0]) || !isNumTok(u[1])) return bad('a shell of [GTO] has fewer primitives than it declares');
        exps.push(num(u[0]) * f2); c1.push(num(u[1])); if (sp) c2.push(num(u[2]));
      }
      k += np;
      if (sp) {
        shells.push({ atom: atom.index, center: atom.xyz.slice(), l: 0, pure: false, exps, coefs: c1 });
        shells.push({ atom: atom.index, center: atom.xyz.slice(), l: 1, pure: false, exps: exps.slice(), coefs: c2 });
      } else shells.push({ atom: atom.index, center: atom.xyz.slice(), l, pure: l >= 2 && (l >= 4 ? pureHigh : pure[l]), exps, coefs: c1 });
    }
    if (!shells.length) return bad('[GTO] lists no shells');
    if (scaled) ds.warnings.push('The [GTO] shells of this Molden file scale their exponents (a factor other than 1.00); the exponents were scaled by its square.');
    if (shells.some(s => s.l >= 5 && !s.pure)) ds.warnings.push('Cartesian h (and higher) functions have no defined order in the Molden format; the orbitals may be wrong.');
    let nbf = 0;
    for (const s of shells) nbf += s.pure ? 2 * s.l + 1 : (s.l + 1) * (s.l + 2) / 2;
    // [MO]: "Key= value" lines (Sym, Ene, Spin, Occup), then "index coefficient" lines
    const list = [];
    let cur = null, over = 0;
    for (const s of mo.lines) {
      const m = /^\s*([A-Za-z]+)\s*=\s*(.*?)\s*$/.exec(s);
      if (m) {
        if (!cur || cur.n) { cur = { sym: null, e: NaN, beta: false, occ: NaN, c: new Float64Array(nbf), n: 0 }; list.push(cur); }
        const key = m[1].toLowerCase(), v = m[2];
        if (key === 'sym') cur.sym = v;
        else if (key === 'ene') cur.e = num(v.split(/\s+/)[0]);
        else if (key === 'spin') cur.beta = /^beta/i.test(v);
        else if (key === 'occup') cur.occ = num(v.split(/\s+/)[0]);
        continue;
      }
      const t = tok(s);
      if (!cur || t.length < 2 || !isIntTok(t[0])) continue;
      for (let q = 0; q + 1 < t.length; q += 2) {
        const i = int(t[q]);
        if (i >= 1 && i <= nbf) cur.c[i - 1] = num(t[q + 1]); else over = Math.max(over, i);
        cur.n++;
      }
    }
    if (over) return bad('[MO] has coefficients for ' + over + ' basis functions, but [GTO] gives ' + nbf);
    const orbs = list.filter(o => o.n);
    if (!orbs.length) return bad('[MO] lists no orbital coefficients');
    const set = arr => {
      if (!arr.length) return null;
      let order = arr.map((o, i) => i);
      const asc = arr.every((o, i) => i === 0 || !(o.e < arr[i - 1].e));
      if (!asc && arr.every(o => Number.isFinite(o.e))) order = order.sort((a, b) => arr[a].e - arr[b].e || a - b);
      else order = null;
      const pick = order ? order.map(i => arr[i]) : arr;
      const coeffs = new Float64Array(pick.length * nbf);
      pick.forEach((o, k) => coeffs.set(o.c, k * nbf));
      const out = { energies: pick.map(o => o.e), occupations: pick.map(o => (Number.isFinite(o.occ) ? o.occ : 0)),
        symmetries: pick.map(o => o.sym), coeffs };
      if (order) out.fileIndex = order.map(i => i + 1);
      return out;
    };
    const alpha = set(orbs.filter(o => !o.beta)), beta = set(orbs.filter(o => o.beta));
    if (!alpha) return bad('[MO] has beta orbitals only');
    ds.basis = { source: 'molden', shells, nbf, program: /orca_2mkl|\bORCA\b/i.test(ds.title || '') ? 'orca' : null };
    ds.mos = { restricted: !beta, alpha, beta };
    const lev = s => s && { energies: s.energies.slice(), nOcc: s.occupations.filter(x => x > 0.5).length, occupations: s.occupations.slice() };
    if (alpha.energies.some(Number.isFinite)) ds.orbitals = { alpha: lev(alpha), beta: lev(beta) };
  }

  // ------------------------------------------------------------------ fchk
  const FCHK_NUM = /[-+]?(?:\d+\.?\d*|\.\d+)(?:[eEdD][-+]?\d+)?/g;
  function parseFchk(text, fileName) {
    const ds = newDataset('fchk', fileName);
    const L = splitLines(text);
    ds.title = (L[0] || '').trim() || ds.title;
    // line 2: the job type in 10 columns, then the method and the basis set (Gaussian 16 gives the method 60 columns,
    // Q-Chem and older versions 30): the basis is the last word after a run of spaces
    const l2 = L[1] || '';
    let jobType = l2.slice(0, 10).trim(), method = '', basis = '';
    const mb = /^\s*(.*?\S)\s{2,}(\S+)$/.exec(l2.slice(10).replace(/\s+$/, ''));
    if (mb) { method = mb[1]; basis = mb[2]; } else method = l2.slice(10).trim();
    if (!method) { const t = tok(l2); jobType = t[0] || ''; method = t[1] || ''; basis = t[2] || ''; }
    const D = Object.create(null), coordSets = [];
    let moGeom = null;
    for (let i = 2; i < L.length; i++) {
      const s = L[i];
      if (!s || s[0] === ' ') continue;
      const m = /^(.+?)\s+([IRCLH])\s+(N=\s*)?(\S+)\s*$/.exec(s.length > 43 && /[IRCLH]/.test(s[43]) ?
        s.slice(0, 40).trim() + '   ' + s.slice(43) : s);
      if (!m) continue;
      const name = m[1].trim(), type = m[2];
      if (m[3]) {
        const cnt = int(m[4]);
        const vals = [];
        let j = i + 1;
        if (type === 'I' || type === 'R') {
          while (j < L.length && vals.length < cnt && (L[j][0] === ' ' || L[j][0] === '-')) {
            // (Psi4 writes 16-column numbers with three-digit exponents that touch: 1.1e+00-3.4e-118)
            const t = /\d[-+]/.test(L[j]) ? L[j].match(FCHK_NUM) || [] : tok(L[j]);
            for (const x of t) vals.push(type === 'I' ? int(x) : num(x));
            j++;
          }
        } else {
          while (j < L.length && L[j][0] === ' ') j++;
        }
        if ((type === 'I' || type === 'R') && vals.length !== cnt) {
          ds.warnings.push('fchk array "' + name + '" holds ' + vals.length + ' values but declares ' + cnt + '.');
        }
        if (type === 'C') {   // text in 12-character words, five to a line (the route, the title): kept as one string
          const nl = Math.ceil(cnt / 5);
          let txt = '';
          for (let k = 0; k < nl && i + 1 + k < L.length; k++) txt += (L[i + 1 + k] || '').padEnd(60).slice(0, 60);
          D[name] = txt.replace(/\s+/g, ' ').trim();
          i += nl;
          continue;
        }
        D[name] = vals.length > cnt ? vals.slice(0, cnt) : vals;
        if (name === 'Current cartesian coordinates' && vals.length === cnt) coordSets.push(D[name]);   // (Q-Chem repeats it per step)
        if (name === 'Alpha MO coefficients') moGeom = coordSets[coordSets.length - 1] || null;   // (the orbitals' geometry)
        i = j - 1;
      } else D[name] = type === 'I' ? int(m[4]) : type === 'R' ? num(m[4]) : m[4];
    }
    const Zs = D['Atomic numbers'], C = D['Current cartesian coordinates'];
    if (!Zs || !C) throw err('fchk file lacks "Atomic numbers" or "Current cartesian coordinates".');
    const nat = Zs.length;
    if (C.length !== 3 * nat) throw err('fchk coordinate array length does not match the atom count.');
    const f = newFrame(Zs.map(z => Math.max(0, z)), C.map(v => v * BOHR));
    const E = D['Total Energy'];
    if (typeof E === 'number') { f.energy = E; f.energies = { Total: E }; }
    if (typeof D['SCF Energy'] === 'number') { f.energies = f.energies || {}; f.energies.SCF = D['SCF Energy']; if (f.energy === null) f.energy = D['SCF Energy']; }
    // a TD / CIS job: "Total Energy" is the excited state's ("CIS Energy"). As in the log reader, it is the frame's
    // energy only in a job that follows the state (optimisation, scan, IRC); otherwise the ground state's is.
    if (typeof D['CIS Energy'] === 'number' && typeof D['SCF Energy'] === 'number') {
      f.energies = { SCF: D['SCF Energy'], TD: D['CIS Energy'] };
      f.energy = /opt|scan|irc/i.test(jobType) ? D['CIS Energy'] : D['SCF Energy'];
    }
    // the optimisation history (Gaussian "Opt point N" sections: every geometry of an optimisation, or of each point
    // of a relaxed scan, with its energy). A checkpoint file keeps the history of an earlier optimisation for later
    // jobs that read their geometry from it, so it is used only for an optimisation or scan job, and only when its
    // last geometry is the file's current one.
    const route = typeof D.Route === 'string' ? D.Route : '';
    const hist = fchkHistory(D, Zs, nat, C);
    if (hist && (/opt|scan/i.test(jobType) || /(^|[\s#])opt\b/i.test(route))) {
      const scan = /scan/i.test(jobType) || hist.points > 1;
      hist.frames.forEach((h, k) => {
        const fr = k === hist.frames.length - 1 ? f : h.frame;
        if (fr !== f && Number.isFinite(h.energy)) { fr.energy = h.energy; fr.energies = { Total: h.energy }; }
        fr.converged = h.last;
        fr.label = scan ? 'Scan point ' + h.point + ', step ' + h.step : 'Step ' + h.step;
        ds.frames.push(fr);
      });
    } else if (coordSets.length > 1 && !/freq/i.test(jobType)) {
      // Q-Chem writes the coordinates again at each step of an optimisation or of a molecular-dynamics run (a
      // frequency job's sets are displacements, ending with the reference): each set is a frame, the last one the file's
      const opt = /opt|scan/i.test(jobType);
      coordSets.forEach((c, k) => {
        if (c.length !== 3 * nat) return;
        const fr = k === coordSets.length - 1 ? f : newFrame(Zs.map(z => Math.max(0, z)), c.map(v => v * BOHR));
        if (k > 0 && k < coordSets.length - 1 && c.every((v, q) => v === coordSets[k - 1][q])) return;
        fr.label = (opt ? 'Step ' : 'Frame ') + (k + 1);
        ds.frames.push(fr);
      });
    } else ds.frames.push(f);
    const meta = ds.meta;
    meta.route = route || null;
    meta.charge = typeof D.Charge === 'number' ? D.Charge : null;
    meta.multiplicity = typeof D.Multiplicity === 'number' ? D.Multiplicity : null;
    meta.method = method && !/^(RO|R|U)$/.test(method) ? method.replace(/^(RO|R|U)(?=[A-Z0-9])/, '') : null;   // (Q-Chem writes only "R")
    meta.basis = basis || null;
    const jt = jobType.toLowerCase();
    if (/opt/.test(jt)) addJobType(meta, 'opt');
    if (/freq/.test(jt)) addJobType(meta, 'freq');
    if (/scan/.test(jt)) addJobType(meta, 'scan');
    if (/irc/.test(jt)) addJobType(meta, 'irc');
    // excited states: Gaussian's "ETran" arrays (energy, transition dipole... of each state) or Q-Chem's
    // "Excitation Energies" (Eh) and "Oscillator Strengths"
    const states = fchkStates(D);
    if (states) { ds.excitedStates = states; addJobType(meta, 'td'); }
    if (typeof D['CIS Energy'] === 'number') addJobType(meta, 'td');
    finishJobTypes(meta);
    // Q-Chem and Psi4 write fchk files too (Q-Chem titles them "Jobname.Temp" and adds Dipole_Data)
    const title1 = (L[0] || '').trim();
    ds.program = title1 === 'Generated by Psi4' ? 'Psi4 (formatted checkpoint)'
      : title1 === 'Jobname.Temp' || D.Dipole_Data !== undefined ? 'Q-Chem (formatted checkpoint)' : 'Gaussian (formatted checkpoint)';
    const ch = {};
    [['Mulliken Charges', 'Mulliken'], ['ESP Charges', 'ESP'], ['NPA Charges', 'NPA']].forEach(([k, n]) => {
      if (Array.isArray(D[k]) && D[k].length === nat) ch[n] = D[k];
    });
    if (Object.keys(ch).length) ds.charges = ch;
    const dm = D['Dipole Moment'], dq = D['Dipole_Data'];   // Gaussian: atomic units; Q-Chem's fchk: Debye
    if (Array.isArray(dm) && dm.length === 3) {
      ds.dipole = { x: dm[0] * AU2DEBYE, y: dm[1] * AU2DEBYE, z: dm[2] * AU2DEBYE, total: Math.hypot(dm[0], dm[1], dm[2]) * AU2DEBYE };
    } else if (Array.isArray(dq) && dq.length === 3 && dq.every(Number.isFinite)) {
      ds.dipole = { x: dq[0], y: dq[1], z: dq[2], total: Math.hypot(dq[0], dq[1], dq[2]) };
    }
    const ae = D['Alpha Orbital Energies'], be = D['Beta Orbital Energies'];
    if (Array.isArray(ae) && ae.length) {
      const na = D['Number of alpha electrons'], nb = D['Number of beta electrons'];
      ds.orbitals = { alpha: { energies: ae, nOcc: typeof na === 'number' ? na : 0 },
        beta: Array.isArray(be) && be.length ? { energies: be, nOcc: typeof nb === 'number' ? nb : 0 } : null };
    }
    fchkBasis(D, ds, nat, moGeom && moGeom.length === C.length ? moGeom : C, coordSets[0]);
    const E2 = D['Vib-E2'], VM = D['Vib-Modes'];
    if (Array.isArray(E2) || Array.isArray(VM)) {
      const nm = Array.isArray(VM) ? VM.length / (3 * nat) : 0;
      if (!Array.isArray(E2) || !Array.isArray(VM) || !Number.isInteger(nm) || nm < 1 ||
          E2.length < 4 * nm || E2.length % nm !== 0 || E2.slice(0, nm).some(v => !Number.isFinite(v))) {
        ds.warnings.push('Vibrational data in the fchk file are inconsistent (Vib-E2 / Vib-Modes sizes) and were skipped.');
      } else {
        const modes = [];
        for (let k = 0; k < nm; k++) modes.push(Float64Array.from(VM.slice(3 * nat * k, 3 * nat * (k + 1))));
        ds.vibrations = { frameIndex: 0, freqs: E2.slice(0, nm), ir: E2.slice(3 * nm, 4 * nm), raman: null, depol: null,
          redMass: E2.slice(nm, 2 * nm), forceConst: E2.slice(2 * nm, 3 * nm), sym: null, modes };
      }
    }
    return ds;
  }

  /** The "Opt point N" sections of a Gaussian fchk file: [{frame, energy, point, step, last}] in order, or null when
   *  there are none or their last geometry is not the file's current one (cur, Bohr). */
  function fchkHistory(D, Zs, nat, cur) {
    const pts = new Map();
    for (const k of Object.keys(D)) {
      const m = /^Opt point\s+(\d+)\s+(Geometries|Results for each geome)/.exec(k);
      if (!m || !Array.isArray(D[k])) continue;
      const p = int(m[1]);
      if (!pts.has(p)) pts.set(p, {});
      pts.get(p)[m[2] === 'Geometries' ? 'G' : 'R'] = D[k];
    }
    if (!pts.size) return null;
    const per = typeof D['Optimization Num results per geometry'] === 'number' && D['Optimization Num results per geometry'] > 0 ? D['Optimization Num results per geometry'] : 2;
    const out = [], nums = Array.from(pts.keys()).sort((a, b) => a - b);
    let points = 0;
    for (const p of nums) {
      const { G, R } = pts.get(p), ng = G ? Math.floor(G.length / (3 * nat)) : 0;
      if (ng) points++;
      for (let g = 0; g < ng; g++) {
        const xyz = G.slice(3 * nat * g, 3 * nat * (g + 1));
        out.push({ frame: newFrame(Zs.map(z => Math.max(0, z)), xyz.map(v => v * BOHR)), xyzBohr: xyz,
          energy: R && Number.isFinite(R[per * g]) ? R[per * g] : null, point: points, step: g + 1, last: g === ng - 1 });
      }
    }
    if (out.length < 2) return null;
    const lastG = out[out.length - 1].xyzBohr;
    for (let q = 0; q < cur.length; q++) if (Math.abs(lastG[q] - cur[q]) > 1e-4) return null;
    return { frames: out, points };
  }
  /** The basis set and MO coefficients of an fchk file as ds.basis and ds.mos (the layout of BM.orbitals: SP shells
   *  split into s + p, negative shell types pure, centres in Bohr). Occupations come from the electron counts. */
  function fchkBasis(D, ds, nat, cur, first) {
    const st = D['Shell types'], npr = D['Number of primitives per shell'], sa = D['Shell to atom map'],
      ex = D['Primitive exponents'], cc = D['Contraction coefficients'], pc = D['P(S=P) Contraction coefficients'],
      sc = D['Coordinates of each shell'], ca = D['Alpha MO coefficients'], cb = D['Beta MO coefficients'];
    if (!Array.isArray(st) || !Array.isArray(ca)) return;
    const bad = why => { ds.warnings.push('The orbitals of this fchk file cannot be computed: ' + why + '.'); };
    const ns = st.length;
    if (![npr, sa, ex, cc].every(Array.isArray) || npr.length !== ns || sa.length !== ns) return bad('its shell arrays are incomplete');
    const shells = [];
    let p = 0, nbf = 0;
    for (let k = 0; k < ns; k++) {
      const t = st[k], n = npr[k], atom = sa[k] - 1;
      if (!Number.isInteger(t) || !(n >= 1) || p + n > ex.length || p + n > cc.length || !(atom >= 0 && atom < nat)) return bad('shell ' + (k + 1) + ' is malformed');
      let center = Array.isArray(sc) && sc.length === 3 * ns ? [sc[3 * k], sc[3 * k + 1], sc[3 * k + 2]]
        : [cur[3 * atom], cur[3 * atom + 1], cur[3 * atom + 2]];
      // Q-Chem writes the shells once, at the first geometry, and the orbitals after a later one (optimisations)
      if (first && first !== cur && first.length === cur.length && [0, 1, 2].every(a => Math.abs(center[a] - first[3 * atom + a]) < 1e-5)) {
        center = [cur[3 * atom], cur[3 * atom + 1], cur[3 * atom + 2]];
      }
      const exps = ex.slice(p, p + n), coefs = cc.slice(p, p + n);
      if (t === -1) {
        if (!Array.isArray(pc) || pc.length < p + n) return bad('an SP shell has no P(S=P) coefficients');
        shells.push({ atom, center, l: 0, pure: false, exps, coefs });
        shells.push({ atom, center: center.slice(), l: 1, pure: false, exps: exps.slice(), coefs: pc.slice(p, p + n) });
        nbf += 4;
      } else {
        const l = Math.abs(t), pure = t < -1;
        shells.push({ atom, center, l, pure, exps, coefs });
        nbf += pure ? 2 * l + 1 : (l + 1) * (l + 2) / 2;
      }
      p += n;
    }
    if (typeof D['Number of basis functions'] === 'number' && D['Number of basis functions'] !== nbf) return bad('its shells make ' + nbf + ' basis functions, but it declares ' + D['Number of basis functions']);
    const nmo = ca.length / nbf;
    if (!Number.isInteger(nmo) || nmo < 1) return bad('the MO coefficients do not fit the basis set');
    const na = typeof D['Number of alpha electrons'] === 'number' ? D['Number of alpha electrons'] : 0;
    const nb = typeof D['Number of beta electrons'] === 'number' ? D['Number of beta electrons'] : na;
    let unres = Array.isArray(cb) && cb.length === ca.length;
    if (unres && cb.every((v, i) => v === ca[i])) unres = false;   // (some writers repeat restricted orbitals as beta)
    const ea = D['Alpha Orbital Energies'], eb = D['Beta Orbital Energies'];
    const occ = f => Array.from({ length: nmo }, (_, i) => f(i));
    const spin = (C, E, o) => ({ energies: Array.isArray(E) && E.length === nmo ? E.slice() : new Array(nmo).fill(NaN),
      occupations: o, symmetries: null, coeffs: Float64Array.from(C) });
    ds.basis = { source: 'fchk', shells, nbf };
    ds.mos = unres
      ? { restricted: false, alpha: spin(ca, ea, occ(i => (i < na ? 1 : 0))), beta: spin(cb, eb, occ(i => (i < nb ? 1 : 0))) }
      : { restricted: true, alpha: spin(ca, ea, occ(i => (i < nb ? 2 : i < na ? 1 : 0))), beta: null };
  }
  /** Excited states of an fchk file (Gaussian ETran arrays, Q-Chem excitation arrays), or null. */
  function fchkStates(D) {
    const mk = (k, de, f, label) => {
      const ev = de * HARTREE2EV;
      return { index: k + 1, energyEV: ev, wavelengthNM: ev > 0 ? EV2NM / ev : NaN, f, label, s2: null, transitions: [] };
    };
    const ets = D['ETran scalars'], etv = D['ETran state values'], ref = D['SCF Energy'];
    if (Array.isArray(ets) && Array.isArray(etv) && ets.length >= 2 && typeof ref === 'number') {
      // per state: total energy, transition dipole (length), (velocity), magnetic dipole, ... (ets[1] values)
      const n = ets[0], per = ets[1], spin = D['ETran spin'];
      if (n > 0 && per >= 4 && etv.length >= n * per) {
        const out = [];
        for (let k = 0; k < n; k++) {
          const de = etv[k * per] - ref, mx = etv[k * per + 1], my = etv[k * per + 2], mz = etv[k * per + 3];
          const sp = Array.isArray(spin) ? spin[k] : null;   // 0: singlet, 2: triplet
          out.push(mk(k, de, 2 / 3 * de * (mx * mx + my * my + mz * mz), sp === 0 ? 'Singlet' : sp === 2 ? 'Triplet' : ''));
        }
        return out;
      }
    }
    const xe = D['Excitation Energies'], xf = D['Oscillator Strengths'];
    if (Array.isArray(xe) && xe.length) {
      return xe.map((e, k) => mk(k, e, Array.isArray(xf) && Number.isFinite(xf[k]) ? xf[k] : null, ''));
    }
    return null;
  }

  // --------------------------------------------------------------- MOL/SDF
  function bondOrder(t) { return t === 1 || t === 2 || t === 3 ? t : t === 4 ? 1.5 : 1; }

  function parseMolV3000(rec, ds) {
    const lines = [];
    for (let k = 0; k < rec.length; k++) {
      let s = rec[k];
      if (!s.startsWith('M  V30')) continue;
      s = s.slice(6).trim();
      while (s.endsWith('-') && k + 1 < rec.length) { s = s.slice(0, -1) + rec[++k].replace(/^M  V30\s?/, '').trim(); }
      lines.push(s);
    }
    let mode = null;
    const Z = [], xyz = [], bonds = [], idmap = {};
    for (const s of lines) {
      if (/^BEGIN ATOM/i.test(s)) { mode = 'a'; continue; }
      if (/^BEGIN BOND/i.test(s)) { mode = 'b'; continue; }
      if (/^END/i.test(s)) { mode = null; continue; }
      const t = s.split(/\s+/);
      if (mode === 'a' && t.length >= 5) {
        idmap[t[0]] = Z.length;
        const z = BM.el.toZ(t[1]);
        Z.push(z < 0 ? 0 : z); xyz.push(num(t[2]), num(t[3]), num(t[4]));
      } else if (mode === 'b' && t.length >= 4 && idmap[t[2]] !== undefined && idmap[t[3]] !== undefined) {
        bonds.push([idmap[t[2]], idmap[t[3]], bondOrder(int(t[1]))]);
      }
    }
    if (!Z.length) throw err('No atoms found in the V3000 connection table.');
    ds.frames.push(newFrame(Z, xyz));
    ds.bonds = bonds;
    return ds;
  }

  /** One V2000 record: {Z, xyz, bonds, energy}. Errors name the record when the file holds several. */
  function molRecord(rec, ds, nRec) {
    const where = nRec > 1 ? ' (record ' + nRec + ')' : '';
    const counts = rec[3];
    let na = int(counts.slice(0, 3)), nb = int(counts.slice(3, 6));
    if (!(na >= 0) || !(nb >= 0)) { const t = tok(counts); na = int(t[0]); nb = int(t[1]); }
    if (!(na >= 0) || !(nb >= 0)) throw err('Cannot read the atom/bond counts line' + where + '.');
    const Z = [], xyz = [];
    for (let k = 0; k < na; k++) {
      const s = rec[4 + k];
      if (s === undefined) throw err('The atom block is truncated' + where + '.');
      let x = num(s.slice(0, 10)), y = num(s.slice(10, 20)), z = num(s.slice(20, 30));
      let sym = s.slice(31, 34).trim();
      if (!Number.isFinite(x) || !Number.isFinite(y) || !Number.isFinite(z) || !sym) {
        const t = tok(s); x = num(t[0]); y = num(t[1]); z = num(t[2]); sym = t[3] || '';
      }
      if (!Number.isFinite(x) || !Number.isFinite(y) || !Number.isFinite(z)) {
        throw err('Atom line ' + (5 + k) + where + ' has no readable coordinates (the atom block may be truncated).');
      }
      let zz = sym === 'D' || sym === 'T' ? 1 : BM.el.toZ(sym);
      if (zz < 0) { ds.warnings.push('Unknown atom symbol "' + sym + '" read as a dummy atom.'); zz = 0; }
      Z.push(zz); xyz.push(x, y, z);
    }
    const bonds = [];
    for (let k = 0; k < nb; k++) {
      const s = rec[4 + na + k];
      if (s === undefined) { ds.warnings.push('The bond block is truncated.'); break; }
      let a = int(s.slice(0, 3)), b = int(s.slice(3, 6)), o = int(s.slice(6, 9));
      if (!(a > 0 && b > 0)) { const t = tok(s); a = int(t[0]); b = int(t[1]); o = int(t[2]); }
      if (a > 0 && b > 0 && a <= na && b <= na) bonds.push([a - 1, b - 1, bondOrder(o)]);
    }
    // the energy BlueMolecule's own "SDF (all geometries)" stores, in hartree ("> <ENERGY>")
    let energy = null;
    if (/BlueMol/.test(rec[1] || '')) {
      const k = rec.findIndex(s => /^>\s+.*<ENERGY>/.test(s));
      if (k >= 0 && isNumTok((rec[k + 1] || '').trim())) energy = num(rec[k + 1].trim());
    }
    return { Z, xyz, bonds, energy };
  }
  function parseMol(text, fileName) {
    const ds = newDataset('mol', fileName);
    const L = splitLines(text);
    // records, separated by "$$$$" lines
    const recs = [];
    let st = 0;
    for (let k = 0; k < L.length; k++) if (L[k].startsWith('$$$$')) { recs.push(L.slice(st, k)); st = k + 1; }
    const tail = L.slice(st);
    if (!recs.length || tail.some(s => s.trim())) recs.push(tail);
    const rec = recs[0];
    if (rec.length < 4) throw err('The molfile header (3 lines) and counts line are missing.');
    ds.title = rec[0].trim() || stripExt(fileName);
    ds.program = rec[1].slice(2, 10).trim() || null;
    if (/V3000/.test(rec[3])) {
      parseMolV3000(rec, ds);
      const Z0 = ds.frames[0].Z;
      let same = true;
      for (let r = 1; r < recs.length && same; r++) {
        let d2 = null;
        if (recs[r].length >= 4 && /V3000/.test(recs[r][3])) { try { d2 = parseMolV3000(recs[r], newDataset('mol', fileName)); } catch (e) { d2 = null; } }
        const f2 = d2 && d2.frames[0];
        if (!f2 || f2.Z.length !== Z0.length || f2.Z.some((z, i) => z !== Z0[i])) { same = false; break; }
        f2.label = recs[r][0].trim() || null;
        ds.frames.push(f2);
      }
      if (!same) {
        ds.frames.length = 1;
        ds.warnings.push('This SD file holds ' + recs.length + ' records of different molecules; only the first is shown.');
      } else if (ds.frames.length > 1) ds.frames[0].label = rec[0].trim() || null;
      return ds;
    }
    const first = molRecord(rec, ds, recs.length > 1 ? 1 : 0);
    // several records of the same atoms (conformers, a trajectory, BlueMolecule's "SDF (all geometries)") are the
    // geometries of one structure; records of different molecules: the first is shown
    const more = [];
    let same = recs.length > 1;
    for (let r = 1; r < recs.length && same; r++) {
      if (recs[r].length < 4 || /V3000/.test(recs[r][3])) { same = false; break; }
      let m;
      try { m = molRecord(recs[r], ds, r + 1); } catch (e) { same = false; break; }
      if (m.Z.length !== first.Z.length || m.Z.some((z, i) => z !== first.Z[i])) { same = false; break; }
      m.title = recs[r][0].trim();
      more.push(m);
    }
    const add = (m, label) => { const f = newFrame(m.Z, m.xyz); if (m.energy !== null) f.energy = m.energy; f.label = label || null; ds.frames.push(f); return f; };
    if (same) {
      add(first, rec[0].trim() || null);
      for (const m of more) add(m, m.title || null);
    } else {
      if (recs.length > 1) ds.warnings.push('This SD file holds ' + recs.length + ' records of different molecules; only the first is shown.');
      add(first, null);
    }
    ds.bonds = first.bonds;
    return ds;
  }

  // ------------------------------------------------------------- PDB / GRO
  const STD_RES = new Set(('ALA ARG ASN ASP CYS GLN GLU GLY HIS HID HIE HIP HSD HSE HSP ILE LEU LYS LYN MET PHE ' +
    'PRO SER THR TRP TYR VAL ASH GLH CYX CYM ACE NME NMA NH2 DA DC DG DT DU A C G U T RA RC RG RU DA5 DA3 DC5 DC3 ' +
    'DG5 DG3 DT5 DT3 RA5 RA3 RC5 RC3 RG5 RG3 RU5 RU3 HOH WAT SOL TIP3 TIP4 TIP5 TP3 T3P T4P SPC SPCE H2O DOD ' +
    'POPC POPE POPS POPG DPPC DOPC DMPC CHL1 CHOL PA PC PE OL ST').split(' '));
  const ION_RES = { NA: 'Na', 'NA+': 'Na', SOD: 'Na', NIO: 'Na', K: 'K', 'K+': 'K', POT: 'K', KIO: 'K',
    CL: 'Cl', 'CL-': 'Cl', CLA: 'Cl', CIO: 'Cl', BR: 'Br', 'BR-': 'Br', IOD: 'I', 'I-': 'I', F: 'F', 'F-': 'F',
    LI: 'Li', 'LI+': 'Li', LIT: 'Li', RB: 'Rb', CS: 'Cs', CES: 'Cs', MG: 'Mg', MG2: 'Mg', 'MG2+': 'Mg', MGI: 'Mg',
    CA: 'Ca', CA2: 'Ca', 'CA2+': 'Ca', CAL: 'Ca', ZN: 'Zn', ZN2: 'Zn', 'ZN2+': 'Zn', ZIN: 'Zn', FE: 'Fe', FE2: 'Fe',
    FE3: 'Fe', CU: 'Cu', CU1: 'Cu', CU2: 'Cu', MN: 'Mn', MN3: 'Mn', CO: 'Co', NI: 'Ni', CD: 'Cd', CD2: 'Cd',
    HG: 'Hg', PT: 'Pt', PD: 'Pd', AU: 'Au', AG: 'Ag', BA: 'Ba', SR: 'Sr', AL: 'Al', CR: 'Cr', PB: 'Pb' };
  const TWO_LETTER = new Set('CL BR FE ZN MG MN CU PT PD AU AG LI SE SI RU RH IR OS GD'.split(' '));

  /** Element symbol guessed from an atom name and residue name (PDB/GRO without element columns). */
  function elementFromName(name, resName) {
    const res = String(resName || '').trim().toUpperCase();
    const nm = String(name || '').trim();
    if (ION_RES[res]) return ION_RES[res];
    const core = nm.replace(/^[\d\s']+/, '');
    if (/[a-z]/.test(core.slice(0, 2))) {
      const z = BM.el.toZ(core);
      if (z > 0) return BM.el.symbol(z);
    }
    const letters = core.replace(/[^A-Za-z].*$/, '').toUpperCase();
    if (!letters) return null;
    if (!STD_RES.has(res) && letters.length >= 2 && TWO_LETTER.has(letters.slice(0, 2)) && !/^[A-Za-z]/.test(core.slice(2))) {
      return letters[0] + letters[1].toLowerCase();
    }
    if (letters[0] === 'D' && STD_RES.has(res)) return 'H';    // deuterium in a standard residue
    return letters[0];
  }

  function cellFromParams(a, b, c, al, be, ga) {
    const ca = Math.cos(al * RAD), cb = Math.cos(be * RAD), cg = Math.cos(ga * RAD), sg = Math.sin(ga * RAD);
    const cy = (ca - cb * cg) / sg;
    return [[a, 0, 0], [b * cg, b * sg, 0], [c * cb, c * cy, c * Math.sqrt(Math.max(0, 1 - cb * cb - cy * cy))]];
  }

  function parsePDB(text, fileName) {
    const ds = newDataset('pdb', fileName);
    const L = splitLines(text);
    const models = [];
    const fresh = () => ({ Z: [], xyz: [], names: [], res: [], ids: [], ch: [], bf: [], oc: [], anyCh: false });
    let cur = null, title = [], header = '', unknown = 0;
    const ss = []; // HELIX / SHEET records -> ds.secondaryStructure
    for (let raw of L) {
      if (raw.indexOf('\t') >= 0) raw = raw.replace(/\t/g, ' ');
      const rec = raw.slice(0, 6);
      if (rec.startsWith('ATOM') || rec.startsWith('HETATM')) {
        if (!cur) cur = fresh();
        const name4 = raw.slice(12, 16), resName = raw.slice(17, 21).trim();
        let x = num(raw.slice(30, 38)), y = num(raw.slice(38, 46)), z = num(raw.slice(46, 54));
        if (!Number.isFinite(x) || !Number.isFinite(y) || !Number.isFinite(z)) {
          const dec = (raw.slice(26).match(/-?\d+\.\d+/g) || []).map(num);
          if (dec.length < 3) throw err('Cannot read the coordinates of: "' + raw.trim() + '".');
          x = dec[0]; y = dec[1]; z = dec[2];
        }
        const el = raw.slice(76, 78).trim();
        let zz = -1;
        if (el && /^[A-Za-z]{1,2}$/.test(el)) zz = el.toUpperCase() === 'D' ? 1 : BM.el.toZ(el);
        if (zz < 0) { const g = elementFromName(name4, resName); zz = g ? (g === 'D' ? 1 : BM.el.toZ(g)) : -1; }
        if (zz < 0) { unknown++; zz = 0; }
        cur.Z.push(zz); cur.xyz.push(x, y, z);
        cur.names.push(name4.trim()); cur.res.push(resName); cur.ids.push(int(raw.slice(22, 26)));
        // chain id (column 22), else the segment id (columns 73-76, CHARMM); occupancy and B-factor
        const cid = (raw[21] || ' ').trim() || raw.slice(72, 76).trim();
        if (cid) cur.anyCh = true;
        cur.ch.push(cid); cur.oc.push(num(raw.slice(54, 60))); cur.bf.push(num(raw.slice(60, 66)));
        continue;
      }
      if (rec === 'HELIX ' || rec === 'SHEET ') {
        const helix = rec === 'HELIX ';
        const a = int(raw.slice(helix ? 21 : 22, helix ? 25 : 26)), b = int(raw.slice(33, 37));
        if (Number.isFinite(a) && Number.isFinite(b)) ss.push({ type: helix ? 'helix' : 'sheet', chain: (raw[helix ? 19 : 21] || '').trim(), start: { resId: a }, end: { resId: b } });
        continue;
      }
      if (rec.startsWith('MODEL')) { if (cur && cur.Z.length) models.push(cur); cur = fresh(); continue; }
      if (rec.startsWith('ENDMDL')) { if (cur && cur.Z.length) models.push(cur); cur = null; continue; }
      if (rec === 'END' || rec.startsWith('END ') || raw.trim() === 'END') { if (cur && cur.Z.length) models.push(cur); cur = null; break; }
      if (rec.startsWith('CRYST1')) {
        const p = [raw.slice(6, 15), raw.slice(15, 24), raw.slice(24, 33), raw.slice(33, 40), raw.slice(40, 47), raw.slice(47, 54)].map(num);
        if (p.every(Number.isFinite) && !(p[0] === 1 && p[1] === 1 && p[2] === 1)) ds.cell = cellFromParams(...p);
        continue;
      }
      if (rec.startsWith('TITLE')) title.push(raw.slice(10).trim());
      else if (rec.startsWith('HEADER')) header = raw.slice(10, 50).trim();
    }
    if (cur && cur.Z.length) models.push(cur);
    if (!models.length) throw err('No ATOM or HETATM records found.');
    if (unknown) ds.warnings.push(unknown + ' atom(s) with an unrecognised element were read as dummy atoms.');
    models.forEach((m, k) => {
      const f = newFrame(m.Z, m.xyz);
      f.atomNames = m.names; f.resNames = m.res; f.resIds = m.ids;
      if (m.anyCh) f.chains = m.ch;
      if (m.bf.some(v => Number.isFinite(v) && v !== 0)) f.bfactors = Float32Array.from(m.bf);
      if (m.oc.some(v => Number.isFinite(v) && v !== 1)) f.occupancy = Float32Array.from(m.oc);
      if (models.length > 1) f.label = 'Model ' + (k + 1);
      ds.frames.push(f);
    });
    if (ss.length) ds.secondaryStructure = ss;
    ds.title = title.join(' ').trim() || header || ds.title;
    return ds;
  }

  function parseGRO(text, fileName) {
    const ds = newDataset('gro', fileName);
    const L = splitLines(text);
    let i = 0, unknown = 0;
    while (i < L.length) {
      while (i < L.length && !L[i].trim()) i++;
      if (i >= L.length) break;
      const title = L[i];
      const cnt = (L[i + 1] || '').trim();
      if (!isIntTok(cnt)) {
        if (ds.frames.length) { ds.warnings.push('Stopped reading at line ' + (i + 2) + '.'); break; }
        throw err('Line 2 of a .gro file must hold the number of atoms.');
      }
      const na = int(cnt);
      const first = L[i + 2] || '';
      let w = 8;
      const d1 = first.indexOf('.', 20), d2 = d1 >= 0 ? first.indexOf('.', d1 + 1) : -1;
      if (d1 >= 0 && d2 > d1) w = d2 - d1;
      const Z = [], xyz = [], names = [], res = [], ids = [];
      for (let k = 0; k < na; k++) {
        const s = L[i + 2 + k];
        if (s === undefined) throw err('The .gro file is truncated (expected ' + na + ' atoms).');
        let x = num(s.slice(20, 20 + w)), y = num(s.slice(20 + w, 20 + 2 * w)), z = num(s.slice(20 + 2 * w, 20 + 3 * w));
        const rid = int(s.slice(0, 5)), rn = s.slice(5, 10).trim(), an = s.slice(10, 15).trim();
        if (!Number.isFinite(x) || !Number.isFinite(y) || !Number.isFinite(z)) {
          const t = tok(s);
          x = num(t[3]); y = num(t[4]); z = num(t[5]);
          if (!Number.isFinite(x + y + z)) throw err('Cannot read the coordinates on line ' + (i + 3 + k) + '.');
        }
        const g = elementFromName(an, rn);
        let zz = g ? BM.el.toZ(g) : -1;
        if (zz < 0) { if (!/^(MW|LP|EP|M|V)\d*$/i.test(an)) unknown++; zz = 0; }
        Z.push(zz); xyz.push(x * 10, y * 10, z * 10);
        names.push(an); res.push(rn); ids.push(rid);
      }
      const box = floatsIn(L[i + 2 + na] || '');
      if (!ds.cell && box.length >= 3) {
        const b = box.map(v => v * 10);
        ds.cell = box.length >= 9 ? [[b[0], b[3], b[4]], [b[5], b[1], b[6]], [b[7], b[8], b[2]]] :
          [[b[0], 0, 0], [0, b[1], 0], [0, 0, b[2]]];
      }
      const f = newFrame(Z, xyz);
      f.atomNames = names; f.resNames = res; f.resIds = ids;
      f.label = title.trim() || null;
      const tm = /\bt=\s*(\S+)/.exec(title);
      if (tm) f.time = num(tm[1]);
      ds.frames.push(f);
      i += 3 + na;
    }
    if (!ds.frames.length) throw err('No atoms found.');
    if (unknown) ds.warnings.push(unknown + ' atom name(s) could not be mapped to an element and are shown as dummies.');
    ds.title = (ds.frames[0].label || '').replace(/\s*\bt=.*$/, '').trim() || ds.title;
    return ds;
  }

  // ------------------------------------------------------ detection & entry
  const FORMATS = [
    { id: 'gaussian-log', name: 'Gaussian output', ext: ['log', 'out'], parse: parseGaussianLog },
    { id: 'orca-out', name: 'ORCA output', ext: ['out', 'log', 'orca'], parse: parseOrca },
    { id: 'gaussian-input', name: 'Gaussian input', ext: ['gjf', 'com', 'gau', 'gjc'], parse: parseGaussianInput },
    { id: 'orca-input', name: 'ORCA input', ext: ['inp'], parse: parseOrcaInput },
    { id: 'xyz', name: 'XYZ', ext: ['xyz', 'txt', 'trj'], parse: parseXYZ },
    { id: 'cube', name: 'Gaussian cube', ext: ['cube', 'cub'], parse: parseCube },
    { id: 'molden', name: 'Molden', ext: ['molden', 'mold', 'molf'], parse: parseMolden },
    { id: 'fchk', name: 'Gaussian formatted checkpoint', ext: ['fchk', 'fch', 'fck'], parse: parseFchk },
    { id: 'mol', name: 'MDL Molfile / SD file', ext: ['mol', 'sdf', 'sd', 'mdl'], parse: parseMol },
    { id: 'pdb', name: 'PDB', ext: ['pdb', 'ent'], parse: parsePDB },
    { id: 'gro', name: 'GROMACS .gro', ext: ['gro'], parse: parseGRO }
  ];

  function looksLikeCube(L) {
    if (L.length < 7) return false;
    const a = tok(L[2]), b = [tok(L[3]), tok(L[4]), tok(L[5])];
    return a.length >= 4 && isIntTok(a[0]) && a.slice(1, 4).every(isNumTok) &&
      b.every(t => t.length === 4 && isIntTok(t[0]) && t.slice(1).every(isNumTok));
  }
  function looksLikeGro(L) {
    if (L.length < 4 || !isIntTok((L[1] || '').trim())) return false;
    const s = L[2] || '';
    return s.length >= 44 && isIntTok(s.slice(0, 5).trim()) && /^\s*-?\d+\.\d+$/.test(s.slice(20, 28)) &&
      /^\s*-?\d+\.\d+$/.test(s.slice(28, 36));
  }
  /** Link 0 lines, a "#" route, a blank line, the title, a blank line and the charge and multiplicity (when the head
   *  reaches them): a comment line starting with "#" in another format (Psi4, xtb, shell) is not enough. */
  function looksLikeGjf(L) {
    let k = 0;
    while (k < L.length && (!L[k].trim() || L[k].trim()[0] === '%')) k++;
    if (!(k < L.length && /^\s*#/.test(L[k]))) return false;
    while (k < L.length && L[k].trim()) k++;          // the route (one or more lines)
    while (k < L.length && !L[k].trim()) k++;
    const t0 = k;
    while (k < L.length && L[k].trim()) k++;          // the title (a few lines at most)
    if (k - t0 > 6) return false;
    while (k < L.length && !L[k].trim()) k++;
    if (k >= L.length) return true;                    // (the head ends before: as before, a "#" line decides)
    return /^\s*[-+]?\d+[\s,]+\d+/.test(L[k]);        // charge and multiplicity (ONIOM: several pairs)
  }
  function looksLikeXyz(L) {
    const t0 = tok(L[0] || '');
    if (t0.length === 1 && isIntTok(t0[0])) {
      const t = tok(L[2] || '');
      return t.length >= 4 && BM.el.toZ(t[0]) >= 0 && t.slice(1, 4).every(isNumTok);
    }
    const t = tok(L.find(s => s.trim()) || '');
    return t.length >= 4 && BM.el.toZ(t[0]) >= 0 && t.slice(1, 4).every(isNumTok);
  }

  function sniff(text) {
    const head = text.slice(0, 65536);
    if (/^Number of atoms\s+I\s+\d+/m.test(head)) return 'fchk';
    if (/Entering Gaussian System|Gaussian, Inc\.|This is part of the Gaussian|^ Gaussian \d\d:/m.test(head)) return 'gaussian-log';
    if (/\*\s*O\s{3}R\s{3}C\s{3}A\s*\*|Program Version \d+\.\d+\.\d+\s+-\s+RELEASE/.test(head)) return 'orca-out';
    if (/^\s*\[Molden Format\]/im.test(head) || /^\s*\[(Atoms|FREQ|FR-COORD|GEOMETRIES)\]/im.test(head)) return 'molden';
    if (looksLikeOrcaInp(head)) return 'orca-input';
    const L = head.split(/\r\n|\r|\n/, 12);
    if (L.length >= 4 && /V[23]000\s*$/.test(L[3])) return 'mol';
    if (/^(ATOM  |HETATM).{24}[ \d-]{4}\.\d{3}[ \d-]{4}\.\d{3}/m.test(head)) return 'pdb';   // x and y in columns 31-46 (F8.3)
    if (looksLikeCube(L)) return 'cube';
    if (looksLikeGro(L)) return 'gro';
    if (looksLikeGjf(head.split(/\r\n|\r|\n/, 80))) return 'gaussian-input';
    if (looksLikeXyz(L)) return 'xyz';
    if (/^M {2}END/m.test(head)) return 'mol';
    return null;
  }

  // ------------------------------------------------------ extensible format registry
  // Other modules add readers with BM.registerFormat({ id, name, extensions: ['cif'], binary?, kind?,
  //   sniff?(head, lowerCaseFileName) -> 0..1, parse?(text, fileName) -> Dataset,
  //   parseBinary?(bytes: Uint8Array, fileName) -> Dataset }).
  // sniff() gets the first 64 KB of text; return >= 0.9 only for unmistakable signatures (program banners).
  const EXTRA = [];
  function registerFormat(f) {
    if (!f || !f.id || !f.name || !(f.parse || f.parseBinary)) throw new Error('registerFormat: id, name and parse/parseBinary are required');
    const entry = Object.assign({ extensions: [], binary: false, kind: 'structure' }, f);
    entry.extensions = (entry.extensions || []).map(e => String(e).toLowerCase().replace(/^\./, ''));
    const k = EXTRA.findIndex(x => x.id === entry.id);
    if (k >= 0) EXTRA[k] = entry; else EXTRA.push(entry);
    const listed = BM.formats.findIndex(x => x.id === entry.id);
    const info = { id: entry.id, name: entry.name, extensions: entry.extensions.slice(), binary: !!entry.binary, kind: entry.kind };
    if (listed >= 0) BM.formats[listed] = info; else BM.formats.push(info);
    return entry;
  }
  function extOf(fileName) {
    const name = baseName(fileName || '').toLowerCase();
    if (/\.molden\.input$/.test(name)) return 'molden';
    return name.indexOf('.') >= 0 ? name.split('.').pop() : '';
  }
  function builtinByExt(ext, s) {
    switch (ext) {
      case 'log': case 'out': case 'orca': return s;
      case 'gjf': case 'gjc': case 'gau': return 'gaussian-input';
      case 'inp': return s === 'orca-input' ? s : null;
      case 'com': return s === 'gaussian-log' ? s : (s === 'gaussian-input' || !s ? 'gaussian-input' : s);
      case 'xyz': case 'trj': return 'xyz';
      case 'cube': case 'cub': return 'cube';
      case 'molden': case 'mold': case 'molf': return 'molden';
      case 'fchk': case 'fch': case 'fck': return 'fchk';
      case 'mol': case 'sdf': case 'sd': case 'mdl': return 'mol';
      case 'pdb': case 'ent': return 'pdb';
      case 'gro': return 'gro';
      default: return null;
    }
  }
  function detectFormat(fileName, text) {
    const name = baseName(fileName || '').toLowerCase();
    const ext = extOf(fileName);
    const head = text === undefined || text === null ? null : String(text).slice(0, 65536);
    const s = head === null ? null : sniff(head);
    let best = null, bestScore = 0;
    for (const f of EXTRA) {
      if (f.binary && !f.parse) continue;
      let sc = 0;
      if (head !== null && f.sniff) { try { sc = +f.sniff(head, name) || 0; } catch (e) { sc = 0; } }
      if (f.extensions.indexOf(ext) >= 0) sc += f.sniff ? 0.25 : 0.6;
      if (sc > bestScore) { bestScore = sc; best = f.id; }
    }
    const strongBuiltin = s === 'fchk' || s === 'gaussian-log' || s === 'orca-out' || s === 'molden';
    if (strongBuiltin && bestScore < 1.2) return s;
    if (bestScore >= 0.9) return best;
    const byExt = builtinByExt(ext, s);
    if (byExt) return byExt;
    if (bestScore >= 0.5) return best;
    return s || (bestScore > 0.3 ? best : null);
  }
  function isBinaryFormat(fileName) {
    const ext = extOf(fileName);
    return EXTRA.some(f => f.binary && f.extensions.indexOf(ext) >= 0);
  }
  /** A file name whose contents may be binary although its extension is a text one (Amber .rst7: ASCII or NetCDF). */
  function maybeBinary(fileName) {
    const ext = extOf(fileName);
    return EXTRA.some(f => f.binary && f.magic && (f.magicExtensions || []).indexOf(ext) >= 0);
  }
  /** The first bytes of such a file belong to a binary format (its signature): the binary reader takes it. */
  function sniffBinary(fileName, u8) {
    if (!u8 || !u8.length) return false;
    const ext = extOf(fileName);
    return EXTRA.some(f => f.binary && f.magic && (f.magicExtensions || []).indexOf(ext) >= 0 && f.magic(u8));
  }
  function supportedList() {
    return 'Gaussian output/input/fchk, ORCA output/input, XYZ, cube, Molden, MOL/SDF, PDB, GRO' +
      (EXTRA.length ? ', ' + EXTRA.map(f => f.name).join(', ') : '');
  }

  function validate(ds) {
    if (!ds.frames.length && !(ds.volumes && ds.volumes.length)) throw err('The file contains no atoms.');
    // an energy that could not be read is "not known" (null), never NaN: every consumer tests for numbers
    for (const f of ds.frames) {
      if (typeof f.energy === 'number' && !Number.isFinite(f.energy)) f.energy = null;
      if (f.energies) for (const k of Object.keys(f.energies)) if (typeof f.energies[k] === 'number' && !Number.isFinite(f.energies[k])) delete f.energies[k];
    }
    if (Array.isArray(ds.optSteps)) for (const st of ds.optSteps) if (typeof st.energy === 'number' && !Number.isFinite(st.energy)) st.energy = null;
    const trajOnly = !!(ds.meta && ds.meta.needsTopology);
    let badZ = 0;
    const seenZ = new WeakSet();   // frames of a trajectory often share one element list
    for (const f of ds.frames) {
      if (trajOnly && !f.Z) { if (f.xyz.length % 3) throw err('Internal error: coordinate array length is not a multiple of 3.'); continue; }
      if (f.xyz.length !== 3 * f.Z.length) throw err('Internal error: coordinate array length mismatch.');
      // an element that could not be read (a damaged number) is a dummy atom, never NaN
      if (seenZ.has(f.Z)) continue;
      seenZ.add(f.Z);
      for (let i = 0; i < f.Z.length; i++) if (!(Number.isInteger(f.Z[i]) && f.Z[i] >= 0 && f.Z[i] <= 118)) { f.Z[i] = 0; badZ++; }
    }
    if (badZ) ds.warnings.push(badZ + ' atom' + (badZ > 1 ? 's have an element' : ' has an element') + ' that could not be read; shown as dummy atoms (X).');
    // overflowing ("1e999") or NaN coordinates never reach the viewer: keep the geometries before the first such
    // one (e.g. an MD run that blew up), as the XYZ reader does for an unreadable frame
    for (let k = 0; k < ds.frames.length; k++) {
      const x = ds.frames[k].xyz;
      let i = 0;
      while (i < x.length && Number.isFinite(x[i])) i++;
      if (i === x.length) continue;
      const what = 'geometry ' + (k + 1) + ' has a coordinate that is not a finite number (atom ' + (Math.floor(i / 3) + 1) + ')';
      if (k === 0) throw err(what[0].toUpperCase() + what.slice(1) + '.');
      ds.warnings.push('Stopped at ' + what + '; ' + k + ' of ' + ds.frames.length + ' geometries are shown.');
      ds.frames.length = k;
      // the jobs keep to the geometries that are left
      if (Array.isArray(ds.jobs)) {
        ds.jobs = ds.jobs.filter(j => !(j.frameStart >= k));
        ds.jobs.forEach(j => { if (j.frameEnd >= k) j.frameEnd = k - 1; });
        if (!ds.jobs.length) ds.jobs = null;
      }
      break;
    }
    if (ds.vibrations) {
      const v = ds.vibrations;
      if (!ds.frames[v.frameIndex]) v.frameIndex = ds.frames.length - 1;
    }
    // a charge table cut while it was written (a job still running) fits no geometry: it is left out until it is whole
    if (ds.charges && typeof ds.charges === 'object') {
      const counts = new Set(ds.frames.map(f => (f.Z ? f.Z.length : -1)));
      for (const k of Object.keys(ds.charges)) {
        const q = ds.charges[k];
        if (Array.isArray(q) && !counts.has(q.length)) delete ds.charges[k];
      }
      if (!Object.keys(ds.charges).length) ds.charges = null;
    }
    // orbital energies that stop before the HOMO (a table cut, or a layout a reader missed): the HOMO, the LUMO and the
    // level diagram would be wrong, so they are left out
    if (ds.orbitals && typeof ds.orbitals === 'object') {
      const short = s => !!s && (!s.energies || !s.energies.length || (Number.isFinite(s.nOcc) && s.nOcc > s.energies.length));
      if (short(ds.orbitals.alpha) || short(ds.orbitals.beta)) {
        ds.orbitals = null;
        ds.warnings.push('The orbital energies in this file are incomplete (fewer than the occupied orbitals), so they are not shown.');
      }
    }
    return ds;
  }

  /**
   * Steps whose geometry was not loaded (frames with xyz null: huge outputs read with a step selection, see
   * condense.js) leave ds.frames. ds.steps then lists every step in order, {energy, energies, opt, label,
   * converged, scan, irc, job, frame: index in ds.frames or null}, so the energy and convergence plots keep all
   * of them; job ranges and the vibrations' frame are renumbered.
   */
  function splitPhantoms(ds) {
    const fr = ds.frames;
    if (!fr.some(f => !f.xyz)) return ds;
    const jobOf = new Int32Array(fr.length).fill(-1);
    (ds.jobs || []).forEach((jb, k) => { for (let i = Math.max(0, jb.frameStart); i <= jb.frameEnd && i < fr.length; i++) if (jobOf[i] < 0) jobOf[i] = k; });
    const map = new Int32Array(fr.length), keep = [], steps = [];
    for (let i = 0; i < fr.length; i++) {
      const f = fr[i];
      map[i] = f.xyz ? keep.length : -1;
      steps.push({ energy: f.energy, energies: f.energies, opt: f.opt, label: f.label, converged: f.converged,
        scan: f.scan || null, irc: f.irc || null, job: jobOf[i], frame: f.xyz ? keep.length : null });
      if (f.xyz) keep.push(f);
    }
    if (!keep.length) throw err('None of the geometries of this file was loaded.');
    const before = i => { for (let k = Math.min(i, fr.length - 1); k >= 0; k--) if (map[k] >= 0) return map[k]; return -1; };
    const after = i => { for (let k = Math.max(0, i); k < fr.length; k++) if (map[k] >= 0) return map[k]; return -1; };
    const near = i => { const b = before(i); return b >= 0 ? b : Math.max(0, after(i)); };
    for (const jb of ds.jobs || []) {
      const s = after(jb.frameStart), e = before(jb.frameEnd);
      if (s >= 0 && e >= s) { jb.frameStart = s; jb.frameEnd = e; } else jb.frameStart = jb.frameEnd = near(jb.frameEnd);
    }
    if (ds.vibrations) ds.vibrations.frameIndex = near(ds.vibrations.frameIndex);
    ds.frames = keep;
    ds.steps = steps;
    return ds;
  }

  // ------------------------------------------------------------------ jobs that are still running
  /** Program outputs (as opposed to structure files): Gaussian and ORCA logs and the 'output' readers. */
  function isOutputFormat(fmt) { return !!fmt && (fmt.kind === 'output' || fmt.id === 'gaussian-log' || fmt.id === 'orca-out'); }
  /** Readers that report whether the job terminated; for them, no termination message means "not finished". */
  const TRACKS_TERMINATION = new Set(['gaussian-log', 'orca-out', 'gamess-out', 'nwchem-out', 'qchem-out', 'psi4-out',
    'xtb-out', 'mopac-out', 'cp2k-out']);
  const PARTIAL_LINE_NOTE = 'The file ends in the middle of a line, as the output of a job that is still running does; ' +
    'that incomplete line was not read.';
  /** A running job's output usually ends part-way through a line, e.g. "SCF Done:  E(RB3LYP) =  -0.", and a number
   *  on that line is incomplete (-0. would read as an energy of 0). The text up to the last line break, or null. */
  function withoutPartialLine(text) {
    const c = text.charCodeAt(text.length - 1);
    if (c === 10 || c === 13) return null;
    const k = Math.max(text.lastIndexOf('\n'), text.lastIndexOf('\r'));
    return k > 0 ? text.slice(0, k + 1) : null;
  }
  /** Marks unfinished jobs (meta.running) and notes a dropped partial line. */
  function finishOutput(ds, fmt, partial) {
    if (!ds || !isOutputFormat(fmt)) return ds;
    const m = ds.meta || (ds.meta = {});
    if (partial) {
      if (!Array.isArray(ds.warnings)) ds.warnings = [];
      ds.warnings.push(PARTIAL_LINE_NOTE);
    }
    if (m.normalTermination !== true && m.normalTermination !== false && (partial || TRACKS_TERMINATION.has(fmt.id))) m.running = true;
    return ds;
  }

  /** validate()'s message ("Geometry 1 has a coordinate that is not a finite number") names the file and the reader,
   *  as the readers' own messages do. */
  function checked(fn, shown, fmt) {
    try { return fn(); } catch (e) {
      const m = (e && e.message) || String(e);
      throw err(/^Could not read /.test(m) ? m : 'Could not read "' + shown + '" as ' + fmt.name + ': ' + m);
    }
  }
  function parseFile(fileName, text) {
    const shown = baseName(fileName) || 'file';
    if (text === null || text === undefined) throw err('"' + shown + '" could not be read (no content).');
    text = String(text);
    if (text.charCodeAt(0) === 0xFEFF) text = text.slice(1);
    if (!text.trim()) throw err('"' + shown + '" is empty.');
    if (/[\u0000-\u0008]/.test(text.slice(0, 4096))) {
      throw err('"' + shown + '" looks like a binary file that BlueMolecule cannot read. Supported: ' + supportedList() + '.');
    }
    const id = detectFormat(fileName, text);
    if (!id) throw err('The format of "' + shown + '" was not recognised. Supported: ' + supportedList() + '.');
    const fmt = FORMATS.find(f => f.id === id) || EXTRA.find(f => f.id === id);
    if (!fmt.parse) throw err('"' + shown + '" must be opened as a binary file.');
    let partial = false;
    if (isOutputFormat(fmt)) {
      const cut = withoutPartialLine(text);
      if (cut !== null) { text = cut; partial = true; }
    }
    let ds;
    try { ds = fmt.parse(text, fileName); } catch (e) {
      throw err('Could not read "' + shown + '" as ' + fmt.name + ': ' + (e && e.message ? e.message : e));
    }
    if (ds && !ds.format) ds.format = fmt.id;
    return checked(() => validate(finishOutput(ds, fmt, partial)), shown, fmt);
  }
  /**
   * BM.parseLines(fileName, lines) -> Dataset: like parseFile for a file already split into lines
   * (BM.bigfile.readLines; files above V8's ~512 MiB string limit). The format is detected from the first 64 KB;
   * readers that need the whole text join the lines when that is safe and otherwise report the file as too large.
   * A string argument is passed to parseFile.
   */
  function parseLines(fileName, lines) {
    if (!Array.isArray(lines)) return parseFile(fileName, lines);
    const shown = baseName(fileName) || 'file';
    const side = lines.bmSide || null;   // a condensed huge output (condense.js): its placeholders are U+0001 lines
    if (lines.length && lines[0].charCodeAt(0) === 0xFEFF) lines[0] = lines[0].slice(1);
    const head = headOf(lines, 65536);
    if (!head.trim()) {
      let k = 0;
      while (k < lines.length && !lines[k].trim()) k++;
      if (k === lines.length) throw err('"' + shown + '" is empty.');
    }
    if (!side && /[\u0000-\u0008]/.test(head.slice(0, 4096))) {
      throw err('"' + shown + '" looks like a binary file that BlueMolecule cannot read. Supported: ' + supportedList() + '.');
    }
    const id = side && side.format ? side.format : detectFormat(fileName, head);
    if (!id) throw err('The format of "' + shown + '" was not recognised. Supported: ' + supportedList() + '.');
    const fmt = FORMATS.find(f => f.id === id) || EXTRA.find(f => f.id === id);
    if (!fmt.parse) throw err('"' + shown + '" must be opened as a binary file.');
    // a complete file ends with a line break, so its last element is ''; anything else is a partial last line
    // (a condenser has already set it aside: lines.bmPartial)
    let partial = !!lines.bmPartial;
    if (!side && isOutputFormat(fmt) && lines.length > 1 && lines[lines.length - 1] !== '') { lines[lines.length - 1] = ''; partial = true; }
    let ds;
    try { ds = fmt.parse(lines, fileName); } catch (e) {
      throw err('Could not read "' + shown + '" as ' + fmt.name + ': ' + (e && e.message ? e.message : e));
    }
    if (ds && !ds.format) ds.format = fmt.id;
    if (side) {
      splitPhantoms(ds);
      const ld = side.loaded;
      if (ld) {   // what was loaded, counted on the dataset (the Info tab; "Reload with other choices")
        const v = ds.vibrations;
        ds.meta.loaded = { steps: Object.assign({}, ld.steps, { total: ds.steps ? ds.steps.length : ld.steps.frames ? ld.steps.total : ds.frames.length, loaded: ds.frames.length }),
          modes: v ? Object.assign({}, ld.modes, { total: v.freqs.length, loaded: v.modesLoaded === undefined ? v.freqs.length : v.modesLoaded }) : null, bytes: ld.bytes };
      }
    }
    return checked(() => validate(finishOutput(ds, fmt, partial)), shown, fmt);
  }
  /** Binary files (trajectories etc.). Falls back to text parsing when no binary reader claims the extension. */
  function parseFileBinary(fileName, bytes) {
    const shown = baseName(fileName) || 'file';
    const u8 = bytes instanceof Uint8Array ? bytes : new Uint8Array(bytes);
    const ext = extOf(fileName);
    let cands = EXTRA.filter(f => f.parseBinary && f.extensions.indexOf(ext) >= 0);
    if (!cands.length) cands = EXTRA.filter(f => f.parseBinary && f.magic && f.magic(u8));   // e.g. a NetCDF restart named .rst7
    if (!cands.length) return parseFile(fileName, new TextDecoder('utf-8').decode(u8));
    let lastErr = null;
    for (const fmt of cands) {
      try { const ds = fmt.parseBinary(u8, fileName); if (ds && !ds.format) ds.format = fmt.id; return validate(ds); }
      catch (e) { lastErr = e; }
    }
    throw err('Could not read "' + shown + '" as ' + cands.map(f => f.name).join(' / ') + ': ' + (lastErr && lastErr.message ? lastErr.message : lastErr));
  }

  /** Readers that read a Blob / chunk source themselves (big binary trajectories: every Nth frame). */
  function sourceReaderFor(fileName) { const ext = extOf(fileName); return EXTRA.find(f => f.parseSource && f.extensions.indexOf(ext) >= 0) || null; }
  /** BM.parseSource(fileName, source, opts) -> Promise<Dataset> for formats with such a reader (DCD, TRR, XTC, NetCDF). */
  async function parseSource(fileName, source, opts) {
    const fmt = sourceReaderFor(fileName);
    if (!fmt) throw err('"' + (baseName(fileName) || 'file') + '" has no reader for big files.');
    let ds;
    try { ds = await fmt.parseSource(source, fileName, opts || {}); } catch (e) {
      if (e && e.name === 'AbortError') throw e;
      throw err('Could not read "' + (baseName(fileName) || 'file') + '" as ' + fmt.name + ': ' + (e && e.message ? e.message : e));
    }
    if (ds && !ds.format) ds.format = fmt.id;
    const reload = ds.meta && ds.meta.reload;
    validate(ds);
    if (reload) Object.defineProperty(ds.meta, 'reload', { enumerable: false, configurable: true, value: o => reload(o).then(validate) });
    return ds;
  }
  BM.detectFormat = detectFormat;
  // shared with condense.js (the condensers must classify lines exactly as the readers do)
  Object.defineProperty(BM, '_parse', { value: Object.freeze({ tok, isIntTok, isNumTok, num, int, floatsIn, modRedCoord, DEC10, FREQ_LABEL, STARS }), configurable: true });
  BM.parseSource = parseSource;
  BM.hasSourceReader = fileName => !!sourceReaderFor(fileName);
  /** BM.probeSource(fileName, source) -> Promise<{frames, atoms, perFrame, exact}>: a big trajectory's header only. */
  BM.probeSource = async (fileName, source) => {
    const fmt = sourceReaderFor(fileName);
    return fmt ? fmt.parseSource(source, fileName, { probe: true }) : null;
  };
  BM.parseFile = parseFile;
  BM.parseLines = parseLines;
  BM.qmmmMatch = qmmmMatch;
  BM.qmRegionFromLinks = qmRegionFromLinks;
  BM.oniomLayersFromPopulation = oniomLayersFromPopulation;
  BM.oniomAtomRow = oniomAtomRow;
  BM.linesToText = asText;
  BM.maxJoinChars = 3e8;
  BM.zmatrixToCartesian = zmatrixToCartesian;
  BM.parsers = { gaussianLog: parseGaussianLog, gaussianInput: parseGaussianInput, orca: parseOrca, xyz: parseXYZ,
    cube: parseCube, molden: parseMolden, fchk: parseFchk, molsdf: parseMol, pdb: parsePDB, gro: parseGRO };
  // non-enumerable: the ten classic parsers stay the enumerable set
  Object.defineProperty(BM.parsers, 'orcaInput', { value: parseOrcaInput, enumerable: false });
  BM.formats = FORMATS.map(f => ({ id: f.id, name: f.name, extensions: f.ext.slice(), binary: false, kind: 'builtin' }));
  BM.registerFormat = registerFormat;
  BM.parseFileBinary = parseFileBinary;
  BM.isBinaryFormat = isBinaryFormat;
  BM.maybeBinary = maybeBinary;
  BM.sniffBinary = sniffBinary;
  BM.formatRegistry = EXTRA;
  BM.emptyDataset = (format, fileName) => newDataset(format, fileName);
  BM.units = { BOHR, AU2DEBYE, HARTREE2EV, EV2NM };
  BM.elementFromAtomName = elementFromName;
})(typeof window !== 'undefined' ? (window.BM = window.BM || {}) : (globalThis.BM = globalThis.BM || {}));
