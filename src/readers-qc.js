/* BlueMolecule — readers for quantum-chemistry programs besides Gaussian and ORCA.
 *
 * Registered with BM.registerFormat (ids):
 *   gamess-out, gamess-input        GAMESS-US / Firefly output (.log/.out) and input ($DATA, C1)
 *   nwchem-out, nwchem-input        NWChem output and .nw input
 *   qchem-out, qchem-input          Q-Chem output and .in input
 *   psi4-out, psi4-input            Psi4 output and .dat/.in input
 *   xtb-out                         xtb output
 *   mopac-out, mopac-input          MOPAC .out / .arc and .mop
 *   cp2k-out, cp2k-input            CP2K output and .inp (&COORD, &CELL)
 *   turbomole-coord, turbomole-gradient   Turbomole coord ($coord, Bohr) and gradient files
 *
 * Datasets follow parsers.js (Angstrom, Hartree, Debye, cm-1). frame.opt uses the shared keys energyChange,
 * maxForce, rmsForce, maxDisp, rmsDisp plus gradNorm (gradient norm: xtb, MOPAC GNORM, Turbomole |dE/dxyz|),
 * each optionally with <key>Thr and <key>Conv. Outputs that report optimisation cycles without printing their
 * geometries (xtb, MOPAC, CP2K) put them in ds.optSteps = [{step, energy, opt}].
 * Written from public descriptions of the output layouts; no code from the programs themselves.
 */
(function (BM) {
  'use strict';

  const BOHR = 0.529177210903, AU2DEBYE = 2.541746, EV2NM = 1239.84198;
  const HARTREE2EV = 27.211386245988, KCAL = 627.509474, KB_KCAL = 0.0019872043;   // kcal/(mol K)
  const DEBYE2_AMU_A2_TO_KMMOL = 42.2561;

  // ------------------------------------------------------------------ helpers
  function num(s) {
    if (typeof s === 'number') return s;
    if (s === undefined || s === null) return NaN;
    s = String(s).trim();
    if (!s || s.indexOf('*') >= 0) return NaN;
    if (/[dD]/.test(s)) s = s.replace(/[dD]/, 'e');
    const v = Number(s);
    return Number.isNaN(v) ? parseFloat(s) : v;
  }
  const int = s => { const v = parseInt(s, 10); return Number.isNaN(v) ? NaN : v; };
  const isNum = s => /^[-+]?(\d+\.?\d*|\.\d+)([eEdD][-+]?\d+)?$/.test(s);
  const isInt = s => /^[-+]?\d+$/.test(s);
  const tok = s => { const t = String(s === undefined ? '' : s).trim(); return t ? t.split(/\s+/) : []; };
  const lines = text => (Array.isArray(text) ? text : String(text).split(/\r\n|\r|\n/));   // text or a line array (big files)
  const floats = s => (String(s).match(/[-+]?(?:\d+\.\d*|\.\d+|\d+)(?:[eEdD][-+]?\d+)?/g) || []).map(num);
  const last = a => (a && a.length ? a[a.length - 1] : null);
  const cap = s => (s ? s[0].toUpperCase() + s.slice(1).toLowerCase() : s);
  function elZ(label) {
    const s = String(label || '');
    if (/^(gh|bq|x|xx|du|dum)\b|^gh\(/i.test(s)) return 0;
    const z = BM.el.toZ(s);
    return z < 0 ? 0 : z;
  }
  function addJT(meta, t) { if (t && meta.jobTypes.indexOf(t) < 0) meta.jobTypes.push(t); }
  const JOB_KINDS = ['opt', 'freq', 'td', 'scan', 'irc', 'nmr'];
  const jobKinds = types => { const k = JOB_KINDS.filter(t => types.indexOf(t) >= 0); return k.length ? k : ['sp']; };
  function newThermo() {
    return { zpe: null, eCorr: null, hCorr: null, gCorr: null, eZPE: null, eTotal: null, H: null, G: null, temperature: null };
  }
  /** Completes totals from corrections (or the reverse) given the electronic energy. */
  function thermoTotals(th, eel) {
    if (!th || !Number.isFinite(eel)) return th;
    const pairs = [['zpe', 'eZPE'], ['eCorr', 'eTotal'], ['hCorr', 'H'], ['gCorr', 'G']];
    for (const [c, t] of pairs) {
      if (th[c] !== null && th[t] === null) th[t] = eel + th[c];
      else if (th[t] !== null && th[c] === null) th[c] = th[t] - eel;
    }
    return th;
  }
  // Gaussian / ORCA outputs are never claimed (the ORCA banner also wraps xtb runs started by ORCA).
  const FOREIGN = /Entering Gaussian System|Gaussian, Inc\.|This is part of the Gaussian|^ Gaussian \d\d:|\*\s*O\s{3}R\s{3}C\s{3}A\s*\*|Program Version \d+\.\d+\.\d+\s+-\s+RELEASE/m;

  function mkFrame(Z, xyz) {
    return { Z: Array.from(Z), xyz: Float64Array.from(xyz), energy: null, energies: null, label: null, converged: false,
      opt: null, atomNames: null, resNames: null, resIds: null };
  }
  function sameGeom(f, Z, xyz, tol) {
    if (!f || f.Z.length !== Z.length) return false;
    for (let i = 0; i < Z.length; i++) if (f.Z[i] !== Z[i]) return false;
    for (let i = 0; i < xyz.length; i++) if (!(Math.abs(f.xyz[i] - xyz[i]) <= tol)) return false;
    return true;
  }
  const PRIORITY = ['CCSD(T)', 'CCSD', 'MP2', 'SCF', 'xTB', 'HOF'];
  /** Frame collector: identical consecutive geometries (re-printed final structures) are merged. */
  function Frames(ds) { this.ds = ds; this.merged = false; }
  Frames.prototype.cur = function () { return last(this.ds.frames); };
  Frames.prototype.add = function (Z, xyz, label, tol) {
    const fr = this.ds.frames, c = last(fr);
    if (!Z.length) return c;
    if (c && sameGeom(c, Z, xyz, tol || 2e-5)) { this.merged = true; if (label && !c.label) c.label = label; return c; }
    const f = mkFrame(Z, xyz);
    f.label = label || null;
    fr.push(f);
    this.merged = false;
    return f;
  };
  Frames.prototype.energy = function (kind, v, f) {
    f = f || this.cur();
    if (!f || !Number.isFinite(v)) return;
    f.energies = f.energies || {};
    f.energies[kind] = v;
    for (const k of PRIORITY) if (f.energies[k] !== undefined) { f.energy = f.energies[k]; return; }
    f.energy = v;
  };
  function setOpt(o, key, v, thr, conv) {
    if (!o || !Number.isFinite(v)) return;
    o[key] = v;
    if (Number.isFinite(thr)) { o[key + 'Thr'] = thr; o[key + 'Conv'] = conv === undefined || conv === null ? Math.abs(v) <= thr : !!conv; }
    else if (conv !== undefined && conv !== null) o[key + 'Conv'] = !!conv;
  }
  /** Multi-job bookkeeping (Q-Chem "@@@" jobs): frame ranges per job; call geometry() after each Frames.add. */
  function Jobs(ds) { this.ds = ds; this.list = []; this.cur = null; }
  Jobs.prototype.open = function (label, kinds) {
    this.close();
    this.cur = { index: this.list.length, label: label || null, kinds: kinds || ['sp'], frameStart: -1, frameEnd: -1,
      normalTermination: null, n0: this.ds.frames.length, m: -1, g: false };
    this.list.push(this.cur);
  };
  Jobs.prototype.geometry = function () {
    const j = this.cur;
    if (j && !j.g) { j.g = true; if (this.ds.frames.length === j.n0) j.m = j.n0 - 1; }
  };
  Jobs.prototype.close = function () {
    const j = this.cur;
    if (!j) return;
    this.cur = null;
    const n = this.ds.frames.length;
    if (n > j.n0) { j.frameStart = j.m >= 0 ? j.m : j.n0; j.frameEnd = n - 1; }
    else j.frameStart = j.frameEnd = j.m >= 0 ? j.m : Math.max(0, n - 1);
    delete j.n0; delete j.m; delete j.g;
  };
  function state(index, eV, f, label) {
    return { index, energyEV: eV, wavelengthNM: eV > 0 ? EV2NM / eV : NaN, f: Number.isFinite(f) ? f : null,
      label: label || '', s2: null, transitions: [] };
  }
  function isLinear(xyz) {
    const N = xyz.length / 3;
    if (N <= 2) return true;
    const a = [xyz[3] - xyz[0], xyz[4] - xyz[1], xyz[5] - xyz[2]], la = Math.hypot(...a);
    for (let i = 2; i < N; i++) {
      const b = [xyz[3 * i] - xyz[0], xyz[3 * i + 1] - xyz[1], xyz[3 * i + 2] - xyz[2]];
      const c = [a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]];
      if (Math.hypot(...c) > 1e-3 * la * Math.max(1e-9, Math.hypot(...b))) return false;
    }
    return true;
  }
  /**
   * ds.vibrations from per-mode arrays. modes[k] = displacement list (3N). With dropTR, the 6 (5 if linear)
   * modes with the smallest |frequency| (translations/rotations) are removed.
   */
  function buildVib(ds, fi, freqs, modes, extra, dropTR) {
    extra = extra || {};
    let keep = freqs.map((f, k) => k).filter(k => Number.isFinite(freqs[k]));
    fi = Math.max(0, Math.min(fi, ds.frames.length - 1));
    const fr = ds.frames[fi];
    const N = fr ? fr.Z.length : 0;
    if (dropTR && fr && keep.length >= 3 * N - 1) {
      const nTR = isLinear(fr.xyz) ? 5 : 6;
      const drop = new Set(keep.slice().sort((a, b) => Math.abs(freqs[a]) - Math.abs(freqs[b])).slice(0, nTR));
      keep = keep.filter(k => !drop.has(k));
    }
    if (!keep.length) return null;
    const pick = a => (a && keep.some(k => Number.isFinite(a[k])) ? keep.map(k => (Number.isFinite(a[k]) ? a[k] : NaN)) : null);
    const mv = keep.map(k => {
      const v = new Float64Array(3 * N), s = modes && modes[k];
      if (s) for (let q = 0; q < 3 * N; q++) v[q] = Number.isFinite(s[q]) ? s[q] : 0;
      return v;
    });
    return { frameIndex: fi, freqs: keep.map(k => freqs[k]), ir: pick(extra.ir), raman: pick(extra.raman), depol: null,
      redMass: pick(extra.redMass), forceConst: pick(extra.forceConst),
      sym: extra.sym && keep.some(k => extra.sym[k]) ? keep.map(k => extra.sym[k] || '') : null, modes: mv };
  }
  function finish(ds, noGeom) {
    const meta = ds.meta;
    if (!ds.frames.length) throw new Error(noGeom || 'No geometry found in this file.');
    if (!meta.jobTypes.some(t => ['opt', 'freq', 'irc', 'scan', 'td', 'nmr', 'force'].indexOf(t) >= 0)) addJT(meta, 'sp');
    if (!ds.jobs) {
      ds.jobs = [{ index: 0, label: meta.route || (meta.keywords ? meta.keywords.join(' ') : null) || null,
        kinds: jobKinds(meta.jobTypes), frameStart: 0, frameEnd: ds.frames.length - 1, normalTermination: meta.normalTermination }];
    }
    if (ds.vibrations && !ds.frames[ds.vibrations.frameIndex]) ds.vibrations.frameIndex = ds.frames.length - 1;
    return ds;
  }
  /** Reads a coordinate table: rows parsed by rowFn(tokens) -> [Z, x, y, z] or null (stop). */
  function readTable(L, j, rowFn, scale) {
    const Z = [], xyz = [];
    scale = scale || 1;
    for (; j < L.length; j++) {
      const r = rowFn(tok(L[j]), L[j]);
      if (!r) break;
      Z.push(r[0]); xyz.push(r[1] * scale, r[2] * scale, r[3] * scale);
    }
    return { Z, xyz, end: j };
  }
  const allNum = (a, i, k) => a.length >= i + k && a.slice(i, i + k).every(isNum);
  function reg(id, name, extensions, kind, sniff, parse) {
    BM.registerFormat({ id, name, extensions, kind, sniff: (h, nm) => (FOREIGN.test(h) ? 0 : sniff(h, nm)), parse });
  }

  // =================================================================== GAMESS-US / Firefly
  const GAMESS_SIG = /GAMESS VERSION\s*=|\bFirefly version\b|PC GAMESS version|FROM IOWA STATE UNIVERSITY|EXECUTION OF GAMESS BEGUN/i;
  function kvs(s, o) {
    const re = /([A-Z][A-Z0-9]*)\s*=\s*([^\s=$]+)/gi;
    let m;
    while ((m = re.exec(s))) o[m[1].toUpperCase()] = m[2].toUpperCase();
  }
  function gamessBasis(o) {
    const g = o.GBASIS || '', ng = o.IGAUSS || o.NGAUSS || '';
    if (!g) return null;
    if (g === 'STO') return 'STO-' + (ng || 3) + 'G';
    if (g === 'N21') return (ng || 3) + '-21G';
    if (g === 'N31' || g === 'N311') {
      const tf = v => /^(T|\.T|\.TRUE)/.test(v || '');
      let b = (ng || 6) + '-' + g.slice(1) + (tf(o.DIFFSP) ? (tf(o.DIFFS) ? '++' : '+') : '') + 'G';
      const d = +o.NDFUNC || 0, p = +o.NPFUNC || 0, f = +o.NFFUNC || 0;
      const heavy = d ? (d > 1 ? d : '') + 'd' + (f ? 'f' : '') : '';
      const pol = [heavy, p ? (p > 1 ? p : '') + 'p' : ''].filter(Boolean).join(',');
      return pol ? b + '(' + pol + ')' : b;
    }
    const cc = /^(A?)CC([DTQ56])(C?)$/.exec(g);
    if (cc) return (cc[1] ? 'aug-' : '') + 'cc-p' + (cc[3] ? 'C' : '') + 'V' + cc[2] + 'Z';
    return g;
  }
  function gamessRow(t) {
    if (t.length !== 5 || !allNum(t, 1, 4)) return null;
    let z = Math.round(num(t[1]));
    if (!(z >= 0 && z <= 118)) z = elZ(t[0]);
    return [z, num(t[2]), num(t[3]), num(t[4])];
  }
  /** Orbital energies from "EIGENVECTORS" / "MOLECULAR ORBITALS" blocks. */
  function gamessEig(L, j) {
    const e = [];
    let blank = 0;
    for (; j < L.length; j++) {
      const t = L[j].trim();
      if (!t || /^-+$/.test(t)) { if (++blank > 4) break; continue; }
      const tk = t.split(/\s+/);
      if (!tk.every(isInt)) break;
      let k = j + 1;
      while (k < L.length && !L[k].trim()) k++;
      const ev = tok(L[k]);
      if (ev.length !== tk.length || !ev.every(isNum)) break;
      e.push(...ev.map(num));
      for (j = k + 1; j < L.length && L[j].trim(); j++);
      blank = 0;
    }
    return { e, end: j };
  }
  function gamessFreq(L, i, vib) {
    let h = i - 1;
    while (h > 0 && !L[h].trim()) h--;
    const idx = tok(L[h]).map(int);
    const vals = [];
    for (const x of tok(L[i].trim().replace(/^FREQUENCY:/, ''))) { if (x === 'I') vals[vals.length - 1] *= -1; else vals.push(num(x)); }
    if (!idx.length || idx.length !== vals.length || !idx.every(Number.isFinite)) return i;
    idx.forEach((m, k) => { vib.freqs[m - 1] = vals[k]; vib.modes[m - 1] = vib.modes[m - 1] || []; });
    const put = (arr, rest, f) => tok(rest).forEach((x, k) => { if (idx[k] !== undefined) arr[idx[k] - 1] = f(x); });
    let j = i + 1, atom = -1;
    const coord = (c, rest) => tok(rest).forEach((x, k) => { if (idx[k] !== undefined) vib.modes[idx[k] - 1][3 * atom + c] = num(x); });
    for (; j < L.length; j++) {
      const t = L[j].trim();
      if (!t) continue;
      let m;
      if ((m = /^SYMMETRY:\s*(.*)$/.exec(t))) { put(vib.sym, m[1], x => x); continue; }
      if ((m = /^REDUCED MASS:\s*(.*)$/.exec(t))) { put(vib.redMass, m[1], num); continue; }
      if ((m = /^IR INTENSITY:\s*(.*)$/.exec(t))) { put(vib.ir, m[1], x => num(x) * DEBYE2_AMU_A2_TO_KMMOL); continue; }
      if ((m = /^(\d+)\s+\S+\s+X\s+(.*)$/.exec(t))) { atom = int(m[1]) - 1; coord(0, m[2]); continue; }
      if (atom >= 0 && (m = /^([YZ])\s+(.*)$/.exec(t))) { coord(m[1] === 'Y' ? 1 : 2, m[2]); continue; }
      if (/^[A-Z][A-Z .]*:\s/.test(t) && atom < 0) continue;
      break;
    }
    return j - 1;
  }

  function parseGamess(text, name) {
    const ds = BM.emptyDataset('gamess-out', name), meta = ds.meta, F = new Frames(ds), L = lines(text), n = L.length;
    const opts = {}, charges = {};
    let opttol = 1e-4, nOccA = null, nOccB = null, orbA = null, orbB = null, beta = false, eqPending = false;
    let vib = null, trSkip = null, thermo = null, states = null, normal = null, nserch = null;
    for (let i = 0; i < n; i++) {
      const t = L[i].trim();
      if (!t) continue;
      let m;
      if (t.startsWith('INPUT CARD>')) { kvs(t.slice(11), opts); continue; }
      if ((m = /GAMESS VERSION\s*=\s*(.+?)\s*\*?$/.exec(t))) { ds.program = 'GAMESS (US) ' + m[1]; continue; }
      if ((m = /Firefly version\s+([\d.]+)/i.exec(t))) { ds.program = 'Firefly ' + m[1]; continue; }
      if (t === '$CONTRL OPTIONS' || t === 'BASIS OPTIONS') { for (let j = i + 2; j < n && L[j].trim(); j++) kvs(L[j], opts); continue; }
      if ((m = /^OPTTOL\s*=\s*(\S+)/.exec(t))) { opttol = num(m[1]); continue; }
      if ((m = /^CHARGE OF MOLECULE\s*=\s*(-?\d+)/.exec(t))) { meta.charge = int(m[1]); continue; }
      if ((m = /^SPIN MULTIPLICITY\s*=\s*(\d+)/.exec(t))) { meta.multiplicity = int(m[1]); continue; }
      if ((m = /^NUMBER OF OCCUPIED ORBITALS \(ALPHA\)\s*=\s*(\d+)/.exec(t))) { nOccA = int(m[1]); continue; }
      if ((m = /^NUMBER OF OCCUPIED ORBITALS \(BETA ?\)\s*=\s*(\d+)/.exec(t))) { nOccB = int(m[1]); continue; }
      if (/^ATOM\s+ATOMIC\s+COORDINATES \(BOHR\)/.test(t)) {
        const g = readTable(L, i + 2, gamessRow, BOHR);
        if (g.Z.length && !ds.frames.length) F.add(g.Z, g.xyz, 'Input geometry');
        i = g.end - 1;
        continue;
      }
      if ((m = /BEGINNING GEOMETRY SEARCH POINT NSERCH=\s*(\d+)/.exec(t))) { nserch = int(m[1]); continue; }
      if (/^COORDINATES OF ALL ATOMS ARE \(ANGS\)/.test(t)) {
        let j = i + 1;
        while (j < n && j < i + 4 && !/^\s*-{10,}/.test(L[j])) j++;
        const g = readTable(L, j + 1, gamessRow);
        if (g.Z.length) {
          const f = F.add(g.Z, g.xyz, nserch !== null ? 'NSERCH ' + nserch : null, 5e-5);
          if (eqPending) { f.converged = true; eqPending = false; }
        }
        i = g.end - 1;
        continue;
      }
      if ((m = /^FINAL\s+(\S+)\s+ENERGY IS\s+(\S+)/.exec(t))) { F.energy('SCF', num(m[2])); continue; }
      if ((m = /^E\(MP2\)\s*=\s*(\S+)/.exec(t))) { F.energy('MP2', num(m[1])); continue; }
      if ((m = /NSERCH:?\s*(\d+)\s+E=\s*(\S+)\s+GRAD\.\s*MAX=\s*(\S+)\s+R\.M\.S\.=\s*(\S+)/.exec(t))) {
        const f = F.cur();
        if (f) {
          if (f.energy === null) F.energy('SCF', num(m[2]), f);
          f.opt = f.opt || {};
          setOpt(f.opt, 'maxForce', num(m[3]), opttol);
          setOpt(f.opt, 'rmsForce', num(m[4]), opttol / 3);
          if (!f.label) f.label = 'NSERCH ' + m[1];
        }
        continue;
      }
      if (t.indexOf('EQUILIBRIUM GEOMETRY LOCATED') >= 0) { eqPending = true; continue; }
      if (/^TOTAL MULLIKEN AND LOWDIN ATOMIC POPULATIONS/.test(t)) {
        let j = i + 1;
        while (j < n && j < i + 4 && !/MULL\.POP\./.test(L[j])) j++;
        const mu = [], lo = [];
        for (j++; j < n; j++) { const r = tok(L[j]); if (r.length < 6 || !isInt(r[0])) break; mu.push(num(r[3])); lo.push(num(r[5])); }
        if (mu.length) { charges.Mulliken = mu; charges.Loewdin = lo; }
        i = j - 1;
        continue;
      }
      if (/^DX\s+DY\s+DZ\s+\/D\/\s+\(DEBYE\)/.test(t)) {
        const v = floats(L[i + 1] || '');
        if (v.length >= 4) ds.dipole = { x: v[0], y: v[1], z: v[2], total: v[3] };
        continue;
      }
      if (/-+\s*ALPHA SET\s*-+/.test(t)) { beta = false; continue; }
      if (/-+\s*BETA SET\s*-+/.test(t)) { beta = true; continue; }
      if (t === 'EIGENVECTORS' || t === 'MOLECULAR ORBITALS') {
        const r = gamessEig(L, i + 1);
        if (r.e.length) { if (beta) orbB = r.e; else { orbA = r.e; orbB = null; } }
        i = r.end - 1;
        continue;
      }
      if ((m = /MODES\s+(\d+)\s+TO\s+(\d+)\s+ARE TAKEN AS ROTATIONS AND TRANSLATIONS/.exec(t))) { trSkip = [int(m[1]), int(m[2])]; continue; }
      if (t.startsWith('FREQUENCY:')) {
        if (!vib) vib = { fi: ds.frames.length - 1, freqs: [], modes: [], ir: [], redMass: [], sym: [] };
        i = gamessFreq(L, i, vib);
        continue;
      }
      if ((m = /^THERMOCHEMISTRY AT T=\s*(\S+)\s*K/.exec(t))) { thermo = newThermo(); thermo.temperature = num(m[1]); continue; }
      if (thermo && /^THE HARMONIC ZERO POINT ENERGY IS/.test(t)) { const v = floats(L[i + 1] || ''); if (v.length) thermo.zpe = v[0]; continue; }
      if (thermo && /^KCAL\/MOL\s+KCAL\/MOL\s+KCAL\/MOL/.test(t)) {
        for (let j = i + 1; j < Math.min(n, i + 8); j++) {
          const r = tok(L[j]);
          if (r[0] === 'TOTAL') { thermo.eCorr = num(r[1]) / KCAL; thermo.hCorr = num(r[2]) / KCAL; thermo.gCorr = num(r[3]) / KCAL; break; }
        }
        continue;
      }
      if (/^SUMMARY OF (TDDFT|CIS) RESULTS/.test(t)) {
        states = [];
        for (let j = i + 1; j < n; j++) {
          const r = tok(L[j]);
          if (!r.length) { if (states.length) { i = j; break; } continue; }
          if (r.length >= 8 && isInt(r[0]) && int(r[0]) > 0 && allNum(r, 2, 2)) states.push(state(int(r[0]), num(r[3]), num(last(r)), r[1]));
          else if (states.length) { i = j - 1; break; }
        }
        addJT(meta, 'td');
        continue;
      }
      if (/EXECUTION OF (GAMESS|FIREFLY) TERMINATED NORMALLY/i.test(t)) { normal = true; continue; }
      if (/TERMINATED\s+-?ABNORMALLY/i.test(t)) normal = false;
    }
    if (eqPending && F.cur()) F.cur().converged = true;
    // ---- metadata
    const run = opts.RUNTYP || 'ENERGY';
    if (run === 'OPTIMIZE') addJT(meta, 'opt');
    if (run === 'SADPOINT') { addJT(meta, 'opt'); addJT(meta, 'ts'); }
    if (run === 'HESSIAN' || /^(\.?T|\.TRUE)/.test(opts.HSSEND || '') || vib) addJT(meta, 'freq');
    if (run === 'IRC') addJT(meta, 'irc');
    if ((opts.TDDFT && opts.TDDFT !== 'NONE') || opts.CITYP === 'CIS') addJT(meta, 'td');
    const scf = opts.SCFTYP || null;
    meta.method = opts.DFTTYP && opts.DFTTYP !== 'NONE' ? opts.DFTTYP : opts.MPLEVL === '2' ? 'MP2'
      : opts.CCTYP && opts.CCTYP !== 'NONE' ? opts.CCTYP : scf;
    meta.basis = gamessBasis(opts);
    meta.keywords = Object.keys(opts).filter(k => /^(SCFTYP|RUNTYP|DFTTYP|MPLEVL|CCTYP|TDDFT|GBASIS|NGAUSS|NDFUNC|NPFUNC)$/.test(k) && !/^(NONE|0)$/.test(opts[k])).map(k => k + '=' + opts[k]);
    meta.route = meta.keywords.join(' ') || null;
    if (meta.charge === null && opts.ICHARG !== undefined) meta.charge = int(opts.ICHARG);
    if (meta.multiplicity === null && opts.MULT !== undefined) meta.multiplicity = int(opts.MULT);
    meta.normalTermination = normal;
    if (!ds.program) ds.program = 'GAMESS';
    if (Object.keys(charges).length) ds.charges = charges;
    if (orbA) ds.orbitals = { alpha: { energies: orbA, nOcc: nOccA === null ? 0 : nOccA },
      beta: orbB ? { energies: orbB, nOcc: nOccB === null ? 0 : nOccB } : null };
    if (states && states.length) ds.excitedStates = states;
    if (vib && vib.freqs.length) {
      if (trSkip) for (let k = trSkip[0] - 1; k < trSkip[1]; k++) vib.freqs[k] = NaN;
      ds.vibrations = buildVib(ds, vib.fi, vib.freqs, vib.modes, vib, !trSkip);
    }
    if (thermo) {
      const vf = ds.frames[ds.vibrations ? ds.vibrations.frameIndex : ds.frames.length - 1];
      ds.thermo = thermoTotals(thermo, vf ? vf.energy : NaN);
    }
    return finish(ds, 'No coordinates found in this GAMESS output.');
  }

  /** GAMESS input: $DATA with C1 symmetry (Cartesian, Angstrom unless UNITS=BOHR). */
  function parseGamessInput(text, name) {
    const ds = BM.emptyDataset('gamess-input', name), meta = ds.meta, L = lines(text);
    const opts = {};
    let data = -1;
    for (let i = 0; i < L.length; i++) {
      const t = L[i].replace(/!.*$/, '');
      if (/^\s*\$DATA\b/i.test(t)) { data = i; break; }
      if (/^\s*\$/.test(t) || data < 0) kvs(t, opts);
    }
    if (data < 0) throw new Error('No $DATA group found in this GAMESS input.');
    const title = (L[data + 1] || '').trim();
    const sym = tok(L[data + 2] || '')[0] || '';
    if (!/^C1$/i.test(sym)) throw new Error('Only $DATA with C1 symmetry is read (this input uses ' + (sym || '?') + '); use C1 or open the output.');
    const f = /^BOHR/.test(opts.UNITS || '') ? BOHR : 1;
    const Z = [], xyz = [];
    let j = data + 3;
    for (; j < L.length && !/^\s*\$END\b/i.test(L[j]); j++) {
      const r = gamessRow(tok(L[j].replace(/!.*$/, '')));
      if (r) { Z.push(r[0]); xyz.push(r[1] * f, r[2] * f, r[3] * f); }
    }
    if (!Z.length) throw new Error('The $DATA group has no Cartesian atoms (Z-matrix $DATA is not read).');
    const fr = mkFrame(Z, xyz);
    ds.frames.push(fr);
    ds.title = title || ds.title;
    ds.program = null;
    if (opts.ICHARG !== undefined) meta.charge = int(opts.ICHARG);
    if (opts.MULT !== undefined) meta.multiplicity = int(opts.MULT);
    meta.method = opts.DFTTYP && opts.DFTTYP !== 'NONE' ? opts.DFTTYP : opts.SCFTYP || null;
    meta.basis = gamessBasis(opts);
    const run = opts.RUNTYP || 'ENERGY';
    if (run === 'OPTIMIZE' || run === 'SADPOINT') addJT(meta, 'opt');
    if (run === 'HESSIAN' || /^(\.?T|\.TRUE)/.test(opts.HSSEND || '')) addJT(meta, 'freq');
    meta.keywords = Object.keys(opts).map(k => k + '=' + opts[k]);
    return finish(ds);
  }

  // =================================================================== NWChem
  const NWCHEM_SIG = /Northwest Computational Chemistry Package \(NWChem\)/;
  function nwModes(L, j) {
    const freqs = [], modes = [];
    let cols = null;
    for (; j < L.length; j++) {
      const t = L[j].trim();
      if (!t || /^\(.*\)$/.test(t) || (!cols && /^-+$/.test(t))) continue;
      const tk = t.split(/\s+/);
      if (tk.every(isInt)) { cols = tk.map(x => int(x) - 1); continue; }
      if (cols && /^(P\.)?Frequency$/.test(tk[0])) { cols.forEach((m, k) => { freqs[m] = num(tk[k + 1]); modes[m] = modes[m] || []; }); continue; }
      if (cols && isInt(tk[0]) && tk.length === cols.length + 1 && allNum(tk, 1, cols.length)) {
        const r = int(tk[0]) - 1;
        cols.forEach((m, k) => { modes[m][r] = num(tk[k + 1]); });
        continue;
      }
      break;
    }
    return { freqs, modes, end: j };
  }
  function parseNwchem(text, name) {
    const ds = BM.emptyDataset('nwchem-out', name), meta = ds.meta, F = new Frames(ds), L = lines(text), n = L.length;
    const echo = [], charges = {}, thr = { gmax: 4.5e-4, grms: 3e-4, xmax: 1.8e-3, xrms: 1.2e-3 };
    let inEcho = false, convPending = false, normal = null, nm = null, ir = null, thermo = null, S = null, T = null;
    let states = null, cur = null, orbA = null, orbB = null;
    for (let i = 0; i < n; i++) {
      const s = L[i], t = s.trim();
      if (!t) continue;
      let m;
      if (/^=+\s*echo of input deck\s*=+$/.test(t)) { inEcho = true; continue; }
      if (inEcho) { if (/^=+$/.test(t)) inEcho = false; else echo.push(t); continue; }
      if ((m = /Northwest Computational Chemistry Package \(NWChem\)\s+(\S+)/.exec(t))) { ds.program = 'NWChem ' + m[1]; continue; }
      if (/^Geometry "[^"]*" -> "[^"]*"/.test(t)) {
        let j = i + 1, scale = 1;
        for (; j < n && j < i + 10 && !/^\s*No\.\s+Tag\s+Charge/.test(L[j]); j++) if (/Output coordinates in a\.u\./.test(L[j])) scale = BOHR;
        while (j < n && !/^\s*-{3,}/.test(L[j])) j++;
        const g = readTable(L, j + 1, r => (r.length >= 6 && isInt(r[0]) && allNum(r, 2, 4)
          ? [Math.round(num(r[2])) > 0 ? Math.round(num(r[2])) : elZ(r[1]), num(r[3]), num(r[4]), num(r[5])] : null), scale);
        if (g.Z.length) F.add(g.Z, g.xyz);
        i = g.end - 1;
        continue;
      }
      if ((m = /^Total (DFT|SCF) energy\s*=\s*(\S+)/.exec(t))) { F.energy('SCF', num(m[2])); continue; }
      if ((m = /^Total MP2 energy\s+(\S+)/.exec(t))) { F.energy('MP2', num(m[1])); continue; }
      if ((m = /^Charge\s*:\s*(-?\d+)/.exec(t))) { meta.charge = int(m[1]); continue; }
      if ((m = /^Spin multiplicity\s*:\s*(\d+)/.exec(t))) { meta.multiplicity = int(m[1]); continue; }
      if ((m = /\((gmax|grms|xmax|xrms)\)\s*=\s*(\S+)/.exec(t))) { thr[m[1]] = num(m[2]); continue; }
      if (t === 'Optimization converged') { convPending = true; continue; }
      if (t[0] === '@' && (m = /^@\s+(\d+)\s+(-?\d+\.\d+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)/.exec(t))) {
        const f = F.cur(), step = int(m[1]);
        if (!f) continue;
        if (f.energy === null) F.energy('SCF', num(m[2]), f);
        const o = f.opt = f.opt || {};
        if (step > 0) setOpt(o, 'energyChange', num(m[3]));
        setOpt(o, 'maxForce', num(m[4]), thr.gmax);
        setOpt(o, 'rmsForce', num(m[5]), thr.grms);
        if (step > 0) { setOpt(o, 'rmsDisp', num(m[6]), thr.xrms); setOpt(o, 'maxDisp', num(m[7]), thr.xmax); }
        if (!f.label) f.label = 'Step ' + step;
        if (convPending) { f.converged = true; convPending = false; }
        continue;
      }
      if (/^Mulliken analysis of the total density/.test(t)) {
        let j = i + 1;
        while (j < n && j < i + 8 && !/^\s*-{5,}\s+-{3,}/.test(L[j])) j++;
        const q = [];
        for (j++; j < n; j++) { const r = tok(L[j]); if (r.length < 4 || !isInt(r[0]) || !isInt(r[2]) || !isNum(r[3])) break; q.push(int(r[2]) - num(r[3])); }
        if (q.length) charges.Mulliken = q;
        i = j - 1;
        continue;
      }
      if (/^Multipole analysis of the density/.test(t)) {
        const d = [NaN, NaN, NaN];
        for (let j = i + 1; j < Math.min(n, i + 16); j++) {
          const r = tok(L[j]);
          if (r.length >= 5 && r[0] === '1' && r.slice(1, 4).every(isInt) && isNum(r[4])) { const k = r.slice(1, 4).indexOf('1'); if (k >= 0) d[k] = num(r[4]); }
        }
        if (d.every(Number.isFinite)) ds.dipole = { x: d[0] * AU2DEBYE, y: d[1] * AU2DEBYE, z: d[2] * AU2DEBYE, total: Math.hypot(...d) * AU2DEBYE };
        continue;
      }
      if (/Final (Alpha |Beta )?Molecular Orbital Analysis/.test(t)) {
        const e = [], occ = [];
        let j = i + 1;
        for (; j < n; j++) {
          const q = /^\s*Vector\s+(\d+)\s+Occ=\s*(\S+)\s+E=\s*(\S+)/.exec(L[j]);
          if (q) { e[int(q[1]) - 1] = num(q[3]); occ[int(q[1]) - 1] = num(q[2]); continue; }
          if (/center of mass|Task\s+times|Mulliken analysis|Multipole analysis|Molecular Orbital Analysis|Parallel integral/.test(L[j])) break;
        }
        let h = -1;
        occ.forEach((o, k) => { if (o > 1e-3) h = k; });
        const set = { energies: Array.from(e, v => (v === undefined ? NaN : v)), nOcc: h + 1 };
        if (/Beta/.test(t)) orbB = set; else { orbA = set; orbB = null; }
        i = j - 1;
        continue;
      }
      if (/^NORMAL MODE EIGENVECTORS IN CARTESIAN COORDINATES/.test(t)) { nm = nwModes(L, i + 1); nm.fi = ds.frames.length - 1; i = nm.end - 1; continue; }
      if (/\|\|\s+(Projected\s+)?Infra Red Intensities/.test(t)) {
        const a = [];
        let j = i + 1;
        for (; j < n; j++) {
          const q = /^\s*(\d+)\s+(-?\d+\.\d+)\s+\|\|\s+(\S+)\s+(\S+)\s+(\S+)/.exec(L[j]);
          if (q) { a[int(q[1]) - 1] = num(q[5]); continue; }
          if (a.length || j > i + 4) break;
        }
        ir = a;
        i = j - 1;
        continue;
      }
      if ((m = /^Temperature\s*=\s*(\S+?)K/.exec(t))) { T = num(m[1]); thermo = thermo || newThermo(); thermo.temperature = T; continue; }
      if ((m = /^(Zero-Point correction to Energy|Thermal correction to Energy|Thermal correction to Enthalpy)\s*=.*\(\s*(\S+)\s+au\)/.exec(t))) {
        thermo = thermo || newThermo();
        thermo[{ 'Zero-Point correction to Energy': 'zpe', 'Thermal correction to Energy': 'eCorr', 'Thermal correction to Enthalpy': 'hCorr' }[m[1]]] = num(m[2]);
        continue;
      }
      if ((m = /^Total Entropy\s*=\s*(\S+)\s*cal\/mol-K/.exec(t))) { S = num(m[1]); continue; }
      if ((m = /^Root\s+(\d+)\s+(singlet|triplet|doublet|quartet|quintet)?\s*(\S*)\s+(-?\d+\.\d+)\s+a\.u\.\s*\(?\s*(-?\d+\.\d+)\s*eV/i.exec(t))) {
        if (!states || (states.length && int(m[1]) <= last(states).index && int(m[1]) === 1)) states = [];
        cur = state(int(m[1]), num(m[5]), NaN, cap(m[2] || '') + (m[3] ? (m[2] ? '-' : '') + m[3] : ''));
        states.push(cur);
        continue;
      }
      if (cur && (m = /^(Dipole|Total) Oscillator Strength\s+(\S+)/.exec(t))) { if (m[1] === 'Dipole' || cur.f === null) cur.f = num(m[2]); continue; }
      if (cur && (m = /^Occ\.\s+(\d+)\s+\S+\s+---\s+Virt\.\s+(\d+)\s+\S+\s+(\S+)/.exec(t))) { cur.transitions.push({ from: m[1], to: m[2], coeff: num(m[3]) }); continue; }
      if (/^Total times\s+cpu:/.test(t)) { normal = true; continue; }
      if (/For further details see manual section|There is an error in the input file/.test(t)) normal = false;
    }
    // ---- metadata from the input echo
    const tasks = [];
    let basis = null, xc = null;
    for (const e of echo) {
      let m;
      if ((m = /^task\s+(\w+)(?:\s+(\w+))?/i.exec(e))) tasks.push([m[1].toLowerCase(), (m[2] || 'energy').toLowerCase()]);
      if (!basis && (m = /^\*\s+library\s+(\S+)/i.exec(e))) basis = m[1];
      if (!xc && (m = /^xc\s+(.+)$/i.exec(e))) xc = m[1].trim();
      if (meta.charge === null && (m = /^charge\s+(-?\d+)/i.exec(e))) meta.charge = int(m[1]);
    }
    for (const [th, op] of tasks) {
      if (op === 'optimize' || op === 'saddle') addJT(meta, 'opt');
      if (op === 'saddle') addJT(meta, 'ts');
      if (op === 'freq' || op === 'frequencies' || op === 'hessian') addJT(meta, 'freq');
      if (th === 'tddft') addJT(meta, 'td');
    }
    const th = tasks.length ? tasks[0][0] : null;
    meta.method = th === 'dft' || th === 'tddft' ? (xc ? xc.toUpperCase() : 'DFT') : th === 'scf' ? 'HF' : th ? th.toUpperCase() : null;
    meta.basis = basis;
    meta.route = tasks.map(x => 'task ' + x.join(' ')).join('; ') || null;
    meta.normalTermination = normal;
    if (!ds.program) ds.program = 'NWChem';
    if (Object.keys(charges).length) ds.charges = charges;
    if (orbA) ds.orbitals = { alpha: orbA, beta: orbB };
    if (states && states.length) ds.excitedStates = states;
    if (nm && nm.freqs.length) ds.vibrations = buildVib(ds, nm.fi, nm.freqs, nm.modes, { ir }, true);
    if (thermo) {
      if (S !== null && thermo.hCorr !== null && T) thermo.gCorr = thermo.hCorr - T * S / 1000 / KCAL;
      const vf = ds.frames[ds.vibrations ? ds.vibrations.frameIndex : ds.frames.length - 1];
      ds.thermo = thermoTotals(thermo, vf ? vf.energy : NaN);
    }
    return finish(ds, 'No "Geometry" block found in this NWChem output.');
  }
  /** NWChem .nw input: first geometry block (Angstrom unless "units bohr|au"). */
  function parseNwInput(text, name) {
    const ds = BM.emptyDataset('nwchem-input', name), meta = ds.meta, L = lines(text).map(s => s.replace(/#.*$/, ''));
    const Z = [], xyz = [];
    let i = L.findIndex(s => /^\s*geometry\b/i.test(s));
    if (i < 0) throw new Error('No geometry block found in this NWChem input.');
    const f = /\bunits\s+(bohr|au|a\.u\.)/i.test(L[i]) ? BOHR : 1;
    for (i++; i < L.length && !/^\s*end\s*$/i.test(L[i]); i++) {
      const r = tok(L[i]);
      if (r.length >= 4 && !/^(symmetry|system|zmatrix|zcoord|units|autosym|noautoz|center|nocenter|load)$/i.test(r[0]) && allNum(r, 1, 3)) {
        Z.push(elZ(r[0].replace(/^bq/i, '') && /^bq/i.test(r[0]) ? 'X' : r[0])); xyz.push(num(r[1]) * f, num(r[2]) * f, num(r[3]) * f);
      }
    }
    if (!Z.length) throw new Error('The NWChem geometry block has no Cartesian atoms.');
    ds.frames.push(mkFrame(Z, xyz));
    ds.program = null;
    for (const s of L) {
      let m;
      if ((m = /^\s*charge\s+(-?\d+)/i.exec(s))) meta.charge = int(m[1]);
      if ((m = /^\s*mult\s+(\d+)/i.exec(s))) meta.multiplicity = int(m[1]);
      if ((m = /^\s*\*\s+library\s+(\S+)/i.exec(s)) && !meta.basis) meta.basis = m[1];
      if ((m = /^\s*xc\s+(.+)$/i.exec(s)) && !meta.method) meta.method = m[1].trim().toUpperCase();
      if ((m = /^\s*title\s+"([^"]*)"/i.exec(s))) ds.title = m[1];
      if ((m = /^\s*task\s+(\w+)(?:\s+(\w+))?/i.exec(s))) {
        const op = (m[2] || 'energy').toLowerCase();
        if (op === 'optimize' || op === 'saddle') addJT(meta, 'opt');
        if (op === 'freq' || op === 'frequencies' || op === 'hessian') addJT(meta, 'freq');
        if (/^tddft$/i.test(m[1])) addJT(meta, 'td');
        if (!meta.method && /^scf$/i.test(m[1])) meta.method = 'HF';
        else if (!meta.method && !/^dft$/i.test(m[1])) meta.method = m[1].toUpperCase();
      }
    }
    return finish(ds);
  }

  // =================================================================== Q-Chem
  const QCHEM_SIG = /Welcome to Q-Chem|Q-Chem, Inc\., Pleasanton/;
  /** $molecule / $rem from Q-Chem input text (also the "User input:" echo of an output). */
  function qchemSections(L) {
    const out = { charge: null, mult: null, atoms: [], rem: {}, read: false };
    for (let i = 0; i < L.length; i++) {
      const t = L[i].replace(/!.*$/, '').trim();
      if (/^\$molecule$/i.test(t)) {
        let j = i + 1;
        const cm = tok(L[j]);
        if (cm.length >= 2 && isInt(cm[0]) && isInt(cm[1])) { out.charge = int(cm[0]); out.mult = int(cm[1]); j++; }
        for (; j < L.length && !/^\s*\$end/i.test(L[j]); j++) {
          const r = tok(L[j]);
          if (/^read$/i.test(r[0] || '')) out.read = true;
          if (r.length >= 4 && allNum(r, 1, 3)) out.atoms.push([elZ(r[0]), num(r[1]), num(r[2]), num(r[3])]);
        }
        i = j;
      } else if (/^\$rem$/i.test(t)) {
        let j = i + 1;
        for (; j < L.length && !/^\s*\$end/i.test(L[j]); j++) {
          const m = /^\s*(\w+)\s*(?:=\s*|\s+)(\S+)/.exec(L[j].replace(/!.*$/, ''));
          if (m) out.rem[m[1].toUpperCase()] = m[2];
        }
        i = j;
      }
    }
    return out;
  }
  function qchemJobTypes(rem, meta) {
    const jt = (rem.JOBTYPE || 'sp').toLowerCase();
    const types = [];
    if (jt === 'opt' || jt === 'optimization' || jt === 'ts') types.push('opt');
    if (jt === 'freq' || jt === 'frequency') types.push('freq');
    if (jt === 'rpath') types.push('irc');
    if (jt === 'nmr') types.push('nmr');
    if (rem.CIS_N_ROOTS && +rem.CIS_N_ROOTS > 0) types.push('td');
    if (meta) { types.forEach(t => addJT(meta, t)); if (jt === 'ts') addJT(meta, 'ts'); }
    return types;
  }
  function parseQchem(text, name) {
    const ds = BM.emptyDataset('qchem-out', name), meta = ds.meta, F = new Frames(ds), L = lines(text), n = L.length;
    const jobs = new Jobs(ds), charges = {};
    let normal = null, thermo = null, S = null, states = null, cur = null, nOccA = null, nOccB = null, orb = null;
    let vib = null, firstInput = null;
    for (let i = 0; i < n; i++) {
      const s = L[i], t = s.trim();
      if (!t) continue;
      let m;
      if (t === 'User input:') {
        let j = i + 1;
        while (j < n && /^-+$/.test(L[j].trim())) j++;
        const start = j;
        while (j < n && !/^-{20,}$/.test(L[j].trim())) j++;
        const sec = qchemSections(L.slice(start, j));
        if (!firstInput) firstInput = sec;
        jobs.open('JOBTYPE ' + (sec.rem.JOBTYPE || 'sp') + (sec.rem.METHOD ? ' ' + sec.rem.METHOD : '') + (sec.rem.BASIS ? '/' + sec.rem.BASIS : ''),
          jobKinds(qchemJobTypes(sec.rem, meta)));
        if (meta.charge === null && sec.charge !== null) { meta.charge = sec.charge; meta.multiplicity = sec.mult; }
        i = j;
        continue;
      }
      if ((m = /^Q-Chem\s+(\d+\.\d+(?:\.\d+)?)/.exec(t)) && !ds.program) { ds.program = 'Q-Chem ' + m[1]; continue; }
      if (/^Standard Nuclear Orientation/.test(t)) {
        const scale = /Bohr/i.test(t) ? BOHR : 1;
        let j = i + 1;
        while (j < n && j < i + 4 && !/^\s*-{10,}/.test(L[j])) j++;
        const g = readTable(L, j + 1, r => (r.length === 5 && isInt(r[0]) && allNum(r, 2, 3) ? [elZ(r[1]), num(r[2]), num(r[3]), num(r[4])] : null), scale);
        if (g.Z.length) { F.add(g.Z, g.xyz); jobs.geometry(); }
        i = g.end - 1;
        continue;
      }
      if ((m = /^There are\s+(\d+)\s+alpha and\s+(\d+)\s+beta electrons/.exec(t))) { nOccA = int(m[1]); nOccB = int(m[2]); continue; }
      if ((m = /^Total energy in the final basis set\s*=\s*(\S+)/.exec(t))) { F.energy('SCF', num(m[1])); continue; }
      if ((m = /^(?:RI-)?MP2\s+total energy\s*=\s*(\S+)/i.exec(t))) { F.energy('MP2', num(m[1])); continue; }
      if ((m = /^Optimization Cycle:\s*(\d+)/.exec(t))) { const f = F.cur(); if (f && !f.label) f.label = 'Cycle ' + m[1]; continue; }
      if (/^Maximum\s+Tolerance\s+Cnvgd\?/.test(t)) {
        const f = F.cur(), o = f ? (f.opt = f.opt || {}) : null;
        for (let j = i + 1; j < Math.min(n, i + 5); j++) {
          const q = /^\s*(Gradient|Displacement|Energy change)\s+(\S+)\s+(\S+)\s+(YES|NO)/.exec(L[j]);
          if (!q) break;
          setOpt(o, { Gradient: 'maxForce', Displacement: 'maxDisp', 'Energy change': 'energyChange' }[q[1]], num(q[2]), num(q[3]), q[4] === 'YES');
          i = j;
        }
        continue;
      }
      if (/\*\*\s*OPTIMIZATION CONVERGED\s*\*\*/.test(t)) { const f = F.cur(); if (f) f.converged = true; continue; }
      if (/^Mode:\s+\d/.test(t)) {
        if (!vib) vib = { fi: ds.frames.length - 1, freqs: [], modes: [], ir: [], redMass: [], forceConst: [] };
        const idx = tok(t).slice(1).map(x => int(x) - 1);
        let j = i + 1;
        const put = (arr, rest) => tok(rest).forEach((x, k) => { if (idx[k] !== undefined) arr[idx[k]] = num(x); });
        for (; j < n; j++) {
          const u = L[j].trim();
          let q;
          if (!u) continue;
          if ((q = /^Frequency:\s+(.*)$/.exec(u))) { put(vib.freqs, q[1]); continue; }
          if ((q = /^Force Cnst:\s+(.*)$/.exec(u))) { put(vib.forceConst, q[1]); continue; }
          if ((q = /^Red\. Mass:\s+(.*)$/.exec(u))) { put(vib.redMass, q[1]); continue; }
          if ((q = /^IR Intens:\s+(.*)$/.exec(u))) { put(vib.ir, q[1]); continue; }
          if (/^(IR Active|Raman Active|Raman Intens|Depolar)\s*:/.test(u) || /^X\s+Y\s+Z/.test(u)) continue;
          if (/^TransDip/.test(u)) break;
          const r = u.split(/\s+/);
          if (r.length === 1 + 3 * idx.length && allNum(r, 1, 3 * idx.length)) {
            idx.forEach((md, k) => { (vib.modes[md] = vib.modes[md] || []).push(num(r[1 + 3 * k]), num(r[2 + 3 * k]), num(r[3 + 3 * k])); });
            continue;
          }
          j--;
          break;
        }
        i = j;
        continue;
      }
      if ((m = /^STANDARD THERMODYNAMIC QUANTITIES AT\s+(\S+)\s*K/.exec(t))) { thermo = newThermo(); thermo.temperature = num(m[1]); continue; }
      if (thermo && (m = /^Zero point vibrational energy:\s+(\S+)\s+kcal\/mol/.exec(t))) { thermo.zpe = num(m[1]) / KCAL; continue; }
      if (thermo && (m = /^Total Enthalpy:\s+(\S+)\s+kcal\/mol/.exec(t))) { thermo.hCorr = num(m[1]) / KCAL; continue; }
      if (thermo && (m = /^Total Entropy:\s+(\S+)\s+cal\/mol\.K/.exec(t))) { S = num(m[1]); continue; }
      if (/^Ground-State Mulliken Net Atomic Charges/.test(t)) {
        let j = i + 1;
        while (j < n && j < i + 6 && !/^\s*-{10,}/.test(L[j])) j++;
        const q = [];
        for (j++; j < n; j++) { const r = tok(L[j]); if (r.length < 3 || !isInt(r[0]) || !isNum(r[2])) break; q.push(num(r[2])); }
        if (q.length) charges.Mulliken = q;
        i = j - 1;
        continue;
      }
      if (/^Dipole Moment \(Debye\)/.test(t)) {
        const a = /X\s+(\S+)\s+Y\s+(\S+)\s+Z\s+(\S+)/.exec(L[i + 1] || ''), b = /Tot\s+(\S+)/.exec(L[i + 2] || '');
        if (a) ds.dipole = { x: num(a[1]), y: num(a[2]), z: num(a[3]), total: b ? num(b[1]) : Math.hypot(num(a[1]), num(a[2]), num(a[3])) };
        continue;
      }
      if (/^Orbital Energies \(a\.u\.\)/.test(t)) {
        const sets = { Alpha: { occ: [], virt: [] }, Beta: { occ: [], virt: [] } };
        let set = null, part = null, j = i + 1;
        for (; j < n; j++) {
          const u = L[j].trim();
          if (!u) continue;
          if (/^-{10,}$/.test(u)) { if (set) break; continue; }
          let q;
          if ((q = /^(Alpha|Beta) MOs/.exec(u))) { set = sets[q[1]]; continue; }
          if (/^--\s*Occupied\s*--/.test(u)) { part = 'occ'; continue; }
          if (/^--\s*Virtual\s*--/.test(u)) { part = 'virt'; continue; }
          if (set && part && /^[-\d.\s*]+$/.test(u)) { set[part].push(...tok(u).map(num)); continue; }
          break;
        }
        const mk = x => ({ energies: x.occ.concat(x.virt), nOcc: x.occ.length });
        if (sets.Alpha.occ.length) orb = { alpha: mk(sets.Alpha), beta: sets.Beta.occ.length ? mk(sets.Beta) : null };
        i = j;
        continue;
      }
      if (/^(TDDFT\/TDA|TDDFT|CIS|RPA|SF-TDDFT|SF-CIS)\s+Excitation Energies/.test(t)) { states = []; cur = null; continue; }
      if (states && (m = /^Excited state\s+(\d+):\s+excitation energy \(eV\)\s*=\s*(\S+)/.exec(t))) {
        cur = state(int(m[1]), num(m[2]), NaN, '');
        states.push(cur);
        continue;
      }
      if (cur && (m = /^Multiplicity:\s*(\S+)/.exec(t))) { cur.label = m[1]; continue; }
      if (cur && (m = /^<S\*\*2>\s*:\s*(\S+)/.exec(t))) { cur.s2 = num(m[1]); continue; }
      if (cur && (m = /^Strength\s*:\s*(\S+)/.exec(t))) { cur.f = num(m[1]); continue; }
      if (cur && (m = /^([DS])\(\s*(\d+)\)\s*-->\s*([VS])\(\s*(\d+)\)\s+amplitude\s*=\s*(\S+)(?:\s+(alpha|beta))?/.exec(t))) {
        const occN = m[6] === 'beta' ? nOccB : nOccA;
        const to = m[3] === 'V' && occN !== null ? int(m[4]) + occN : int(m[4]);
        cur.transitions.push({ from: m[2], to: String(to), coeff: num(m[5]) });
        continue;
      }
      if (/Thank you very much for using Q-Chem/.test(t)) { normal = true; if (jobs.cur) jobs.cur.normalTermination = true; continue; }
      if (/Q-Chem fatal error occurred/i.test(t)) { normal = false; if (jobs.cur) jobs.cur.normalTermination = false; }
    }
    jobs.close();
    const rem = firstInput ? firstInput.rem : {};
    meta.method = rem.METHOD || rem.EXCHANGE || null;
    if (meta.method) meta.method = meta.method.toUpperCase();
    meta.basis = rem.BASIS || null;
    meta.keywords = Object.keys(rem).map(k => k + ' ' + rem[k]);
    meta.route = meta.keywords.join('; ') || null;
    if (states && states.length) addJT(meta, 'td');
    meta.normalTermination = normal;
    if (!ds.program) ds.program = 'Q-Chem';
    if (jobs.list.length) ds.jobs = jobs.list;
    if (Object.keys(charges).length) ds.charges = charges;
    if (orb) { orb.alpha.nOcc = nOccA !== null ? nOccA : orb.alpha.nOcc; ds.orbitals = orb; }
    if (states && states.length) ds.excitedStates = states;
    if (vib && vib.freqs.length) ds.vibrations = buildVib(ds, vib.fi, vib.freqs, vib.modes, vib, false);
    if (thermo) {
      const T = thermo.temperature || 298.15;
      if (thermo.hCorr !== null) {
        thermo.eCorr = thermo.hCorr - KB_KCAL * T / KCAL;
        if (S !== null) thermo.gCorr = thermo.hCorr - T * S / 1000 / KCAL;
      }
      const vf = ds.frames[ds.vibrations ? ds.vibrations.frameIndex : ds.frames.length - 1];
      ds.thermo = thermoTotals(thermo, vf ? vf.energy : NaN);
    }
    if (!ds.frames.length && firstInput && firstInput.atoms.length) {
      ds.frames.push(mkFrame(firstInput.atoms.map(a => a[0]), [].concat(...firstInput.atoms.map(a => a.slice(1)))));
    }
    return finish(ds, 'No "Standard Nuclear Orientation" block found in this Q-Chem output.');
  }
  function parseQchemInput(text, name) {
    const ds = BM.emptyDataset('qchem-input', name), meta = ds.meta;
    const sec = qchemSections(lines(text));
    if (!sec.atoms.length) throw new Error(sec.read ? 'The $molecule section reads its geometry from a previous job ("read").' : 'No Cartesian $molecule section found in this Q-Chem input.');
    ds.frames.push(mkFrame(sec.atoms.map(a => a[0]), [].concat(...sec.atoms.map(a => a.slice(1)))));
    ds.program = null;
    meta.charge = sec.charge; meta.multiplicity = sec.mult;
    meta.method = sec.rem.METHOD ? sec.rem.METHOD.toUpperCase() : (sec.rem.EXCHANGE ? sec.rem.EXCHANGE.toUpperCase() : null);
    meta.basis = sec.rem.BASIS || null;
    meta.keywords = Object.keys(sec.rem).map(k => k + ' ' + sec.rem[k]);
    qchemJobTypes(sec.rem, meta);
    return finish(ds);
  }

  // =================================================================== Psi4
  const PSI4_SIG = /Psi4: An Open-Source Ab Initio Electronic Structure Package/;
  function psi4Row(r) {
    if (r.length < 4 || !/^[A-Za-z]/.test(r[0]) || !allNum(r, 1, 3) || (r.length > 5)) return null;
    return [elZ(r[0].replace(/^Gh\((\w+)\)$/i, 'X')), num(r[1]), num(r[2]), num(r[3])];
  }
  function parsePsi4(text, name) {
    const ds = BM.emptyDataset('psi4-out', name), meta = ds.meta, F = new Frames(ds), L = lines(text), n = L.length;
    const input = [], charges = {};
    let normal = null, thermo = null, orb = null, vib = null, basisOut = null, funcOut = null, dip = null;
    for (let i = 0; i < n; i++) {
      const s = L[i], t = s.trim();
      if (!t) continue;
      let m;
      if (t === '==> Input File <==') {
        let j = i + 1;
        while (j < n && !/^-{20,}$/.test(L[j].trim())) j++;
        for (j++; j < n && !/^-{20,}$/.test(L[j].trim()); j++) input.push(L[j]);
        i = j;
        continue;
      }
      if ((m = /^Psi4\s+(\d[\w.+-]*)/.exec(t)) && !ds.program) { ds.program = 'Psi4 ' + m[1]; continue; }
      if ((m = /^Geometry \(in (Angstrom|Bohr)\), charge = (-?\d+), multiplicity = (\d+):/.exec(t))) {
        meta.charge = int(m[2]); meta.multiplicity = int(m[3]);
        let j = i + 1;
        while (j < n && (!L[j].trim() || /^\s*(Center\s+X|-{5,})/.test(L[j]))) j++;
        const g = readTable(L, j, psi4Row, m[1] === 'Bohr' ? BOHR : 1);
        if (g.Z.length) F.add(g.Z, g.xyz, null, 2e-6);
        i = g.end - 1;
        continue;
      }
      if ((m = /^@\S*?(?:RHF|UHF|ROHF|CUHF|RKS|UKS) Final Energy:\s+(\S+)/.exec(t))) { F.energy('SCF', num(m[1])); continue; }
      if ((m = /^\*?\s*(CCSD\(T\)|CCSD|MP2) total energy\s*=\s*(\S+)/i.exec(t))) { F.energy(m[1].toUpperCase(), num(m[2])); continue; }
      if ((m = /^Basis Set:\s+(\S+)/.exec(t)) && !basisOut) { basisOut = m[1]; continue; }
      if ((m = /^=>\s*Composite Functional:\s*(\S+)\s*<=/.exec(t)) && !funcOut) { funcOut = m[1]; continue; }
      if (t === '==> Convergence Check <==') {
        let j = i + 1, crit = null;
        const cols = ['energyChange', 'maxForce', 'rmsForce', 'maxDisp', 'rmsDisp'];
        for (; j < Math.min(n, i + 14); j++) {
          const r = tok(L[j]);
          if (r[0] === 'Convergence' && r[1] === 'Criteria') {
            crit = [];
            for (let k = 2; k < r.length && crit.length < 5; k++) {
              if (isNum(r[k])) { crit.push(num(r[k])); if (r[k + 1] === '*') k++; } else if (r[k] === 'o') crit.push(null);
            }
            continue;
          }
          if (crit && r.length > 6 && isInt(r[0]) && isNum(r[1])) {
            const f = F.cur(), o = f ? (f.opt = f.opt || {}) : null;
            const E = num(r[1]);
            let c = 0;
            for (let k = 2; k < r.length && c < 5; k++) {
              if (!isNum(r[k])) continue;
              const mark = r[k + 1] === '*' ? '*' : r[k + 1] === 'o' ? 'o' : '';
              const v = num(r[k]), key = cols[c], thr = crit[c];
              if (!(c === 0 && Math.abs(v - E) <= 0.01 * Math.abs(E))) setOpt(o, key, v, mark === 'o' || thr === null ? NaN : thr, mark === 'o' ? undefined : mark === '*');
              c++;
            }
            if (f && !f.label) f.label = 'Step ' + r[0];
            break;
          }
        }
        i = j;
        continue;
      }
      if (/Optimization is complete!/.test(t)) { const f = F.cur(); if (f) f.converged = true; addJT(meta, 'opt'); continue; }
      if (/^Vibration\s+\d/.test(t)) {
        if (!vib) vib = { fi: ds.frames.length - 1, freqs: [], modes: [], ir: [], redMass: [], forceConst: [], sym: [] };
        const idx = tok(t).slice(1).map(x => int(x) - 1);
        const put = (arr, rest, f) => tok(rest).forEach((x, k) => { if (idx[k] !== undefined) arr[idx[k]] = f(x); });
        let j = i + 1;
        for (; j < n; j++) {
          const u = L[j].trim();
          let q;
          if (!u || /^-{10,}$/.test(u)) continue;
          if ((q = /^Freq \[cm\^-1\]\s+(.*)$/.exec(u))) { put(vib.freqs, q[1], x => (/i$/.test(x) ? -num(x.slice(0, -1)) : num(x))); continue; }
          if ((q = /^Irrep\s+(.*)$/.exec(u))) { put(vib.sym, q[1], x => x); continue; }
          if ((q = /^Reduced mass \[u\]\s+(.*)$/.exec(u))) { put(vib.redMass, q[1], num); continue; }
          if ((q = /^Force const \[mDyne\/A\]\s+(.*)$/.exec(u))) { put(vib.forceConst, q[1], num); continue; }
          if ((q = /^IR activ \[km\/mol\]\s+(.*)$/.exec(u))) { put(vib.ir, q[1], num); continue; }
          if (/^(Turning point|RMS dev|Char temp)/.test(u)) continue;
          const r = u.split(/\s+/);
          if (r.length === 2 + 3 * idx.length && isInt(r[0]) && allNum(r, 2, 3 * idx.length)) {
            const a = int(r[0]) - 1;
            idx.forEach((md, k) => { const v = vib.modes[md] = vib.modes[md] || []; for (let c = 0; c < 3; c++) v[3 * a + c] = num(r[2 + 3 * k + c]); });
            continue;
          }
          break;
        }
        i = j - 1;
        continue;
      }
      if ((m = /^Correction (ZPE|E|H|G)\s+.*?(-?\d+\.\d+)\s+\[Eh\]\s*$/.exec(t))) {
        thermo = thermo || newThermo();
        thermo[{ ZPE: 'zpe', E: 'eCorr', H: 'hCorr', G: 'gCorr' }[m[1]]] = num(m[2]);
        continue;
      }
      if ((m = /^Total (ZPE|E|H|G),.*?(-?\d+\.\d+)\s+\[Eh\]\s*$/.exec(t))) {
        thermo = thermo || newThermo();
        thermo[{ ZPE: 'eZPE', E: 'eTotal', H: 'H', G: 'G' }[m[1]]] = num(m[2]);
        const tt = /at\s+(\d+\.\d+)\s+\[K\]/.exec(t);
        if (tt && m[1] !== 'ZPE') thermo.temperature = num(tt[1]);
        continue;
      }
      if ((m = /^(Mulliken|Lowdin) Charges: \(a\.u\.\)/.exec(t))) {
        let j = i + 1;
        if (/Center\s+Symbol/.test(L[j] || '')) j++;
        const q = [];
        for (; j < n; j++) { const r = tok(L[j]); if (r.length < 6 || !isInt(r[0]) || !isNum(r[5])) break; q.push(num(r[5])); }
        if (q.length) charges[m[1] === 'Lowdin' ? 'Loewdin' : 'Mulliken'] = q;
        i = j - 1;
        continue;
      }
      if ((m = /^Dipole ([XYZ])\s+:\s+(\S+)\s+(\S+)\s+(\S+)/.exec(t))) { dip = dip || [0, 0, 0]; dip['XYZ'.indexOf(m[1])] = num(m[4]) * AU2DEBYE; continue; }
      if ((m = /^X:\s+(\S+)\s+Y:\s+(\S+)\s+Z:\s+(\S+)\s+Total:\s+(\S+)/.exec(t)) && /Dipole Moment: \[D\]/.test(L[i - 1] || '')) {
        ds.dipole = { x: num(m[1]), y: num(m[2]), z: num(m[3]), total: num(m[4]) }; dip = null;
        continue;
      }
      if (t === 'Orbital Energies [Eh]') {
        const sets = { a: { occ: [], virt: [] }, b: { occ: [], virt: [] } };
        let cur = null, j = i + 1;
        for (; j < n; j++) {
          const u = L[j].trim();
          if (!u || /^-+$/.test(u)) continue;
          if (/^(Doubly Occupied|Singly Occupied|Alpha Occupied):/.test(u)) { cur = sets.a.occ; continue; }
          if (/^(Virtual|Alpha Virtual):/.test(u)) { cur = sets.a.virt; continue; }
          if (/^Beta Occupied:/.test(u)) { cur = sets.b.occ; continue; }
          if (/^Beta Virtual:/.test(u)) { cur = sets.b.virt; continue; }
          const pairs = u.match(/\d+[A-Za-z"'0-9]*\s+-?\d+\.\d+/g);
          if (cur && pairs && pairs.length) { for (const p of pairs) cur.push(num(tok(p)[1])); continue; }
          break;
        }
        const srt = a => a.slice().sort((x, y) => x - y);
        const mk = x => ({ energies: srt(x.occ).concat(srt(x.virt)), nOcc: x.occ.length });
        if (sets.a.occ.length) orb = { alpha: mk(sets.a), beta: sets.b.occ.length ? mk(sets.b) : null };
        i = j - 1;
        continue;
      }
      if (/Psi4 exiting successfully/.test(t)) { normal = true; continue; }
      if (/Psi4 encountered an error/.test(t)) normal = false;
    }
    // ---- metadata from the input echo
    const inp = input.join('\n');
    let m;
    const calls = [];
    const re = /\b(energy|optimize|opt|frequency|frequencies|freq|hessian|gradient|properties|prop)\s*\(\s*['"]([^'"]+)['"]/g;
    while ((m = re.exec(inp))) calls.push([m[1].toLowerCase(), m[2]]);
    for (const [c] of calls) {
      if (c === 'optimize' || c === 'opt') addJT(meta, 'opt');
      if (/^freq|^hessian/.test(c)) addJT(meta, 'freq');
    }
    if (/\btdscf_states\b/i.test(inp)) addJT(meta, 'td');
    const bm = /^\s*set\s+basis\s+(\S+)/im.exec(inp) || /^\s*basis\s+(\S+)/im.exec(inp);
    meta.method = calls.length ? calls[0][1].toUpperCase() : funcOut;
    meta.basis = bm ? bm[1] : basisOut;
    meta.route = calls.map(c => c[0] + "('" + c[1] + "')").join(' ') || null;
    meta.normalTermination = normal;
    if (!ds.program) ds.program = 'Psi4';
    if (dip) ds.dipole = { x: dip[0], y: dip[1], z: dip[2], total: Math.hypot(...dip) };
    if (Object.keys(charges).length) ds.charges = charges;
    ds.orbitals = orb;
    if (vib && vib.freqs.length) ds.vibrations = buildVib(ds, vib.fi, vib.freqs, vib.modes, vib, false);
    if (thermo) {
      const vf = ds.frames[ds.vibrations ? ds.vibrations.frameIndex : ds.frames.length - 1];
      ds.thermo = thermoTotals(thermo, vf ? vf.energy : NaN);
    }
    return finish(ds, 'No "Geometry (in Angstrom)" block found in this Psi4 output.');
  }
  /** Psi4 input: the first molecule { } block (Cartesian lines; charge/multiplicity line; units). */
  function parsePsi4Input(text, name) {
    const ds = BM.emptyDataset('psi4-input', name), meta = ds.meta, L = lines(text).map(s => s.replace(/#.*$/, ''));
    let i = L.findIndex(s => /^\s*molecule\s*\w*\s*\{/.test(s));
    if (i < 0) throw new Error('No "molecule { ... }" block found in this Psi4 input.');
    const Z = [], xyz = [];
    let bohr = false;
    const body = [];
    const first = L[i].replace(/^[^{]*\{/, '');
    if (first.trim()) body.push(first);
    for (i++; i < L.length && L[i].indexOf('}') < 0; i++) body.push(L[i]);
    if (i < L.length) body.push(L[i].slice(0, L[i].indexOf('}')));
    for (const s of body) {
      const r = tok(s);
      if (!r.length) continue;
      if (r.length === 2 && isInt(r[0]) && isInt(r[1]) && !Z.length) { meta.charge = int(r[0]); meta.multiplicity = int(r[1]); continue; }
      if (/^units?$/i.test(r[0])) { bohr = /^(bohr|au|a\.u\.)$/i.test(r[1] || ''); continue; }
      if (r.length >= 4 && /^[A-Za-z]/.test(r[0]) && allNum(r, 1, 3)) { Z.push(elZ(r[0].replace(/^@/, 'X'))); xyz.push(num(r[1]), num(r[2]), num(r[3])); }
    }
    if (!Z.length) throw new Error('The Psi4 molecule block has no Cartesian atoms (Z-matrix input is not read).');
    if (bohr) for (let k = 0; k < xyz.length; k++) xyz[k] *= BOHR;
    ds.frames.push(mkFrame(Z, xyz));
    ds.program = null;
    const inp = L.join('\n');
    const c = /\b(energy|optimize|opt|frequency|frequencies|freq)\s*\(\s*['"]([^'"]+)['"]/.exec(inp);
    if (c) meta.method = c[2].toUpperCase();
    if (/\b(optimize|opt)\s*\(/.test(inp)) addJT(meta, 'opt');
    if (/\b(frequency|frequencies|freq|hessian)\s*\(/.test(inp)) addJT(meta, 'freq');
    const b = /^\s*set\s+basis\s+(\S+)/im.exec(inp);
    if (b) meta.basis = b[1];
    return finish(ds);
  }

  // =================================================================== xtb
  const XTB_SIG = /\|\s+x\s?T\s?B\s+\||\*\s*xtb version\s+\d/;
  const pairsIn = s => { const out = []; const re = /(\d+):\s*(-?\d+\.?\d*)/g; let m; while ((m = re.exec(s))) out.push([int(m[1]) - 1, num(m[2])]); return out; };
  function parseXtb(text, name) {
    const ds = BM.emptyDataset('xtb-out', name), meta = ds.meta, F = new Frames(ds), L = lines(text), n = L.length;
    const steps = [], charges = {}, freqs = [], red = [], ir = [];
    let ethr = NaN, gthr = NaN, conv = null, normal = null, cyc = null, call = '', orb = null, thermo = null, fE = null, list = null, temp = null;
    for (let i = 0; i < n; i++) {
      const t = L[i].trim();
      if (!t) continue;
      let m;
      if ((m = /xtb version\s+(\S+)/.exec(t)) && !ds.program) { ds.program = 'xtb ' + m[1]; continue; }
      if ((m = /^program call\s*:\s*(.*)$/.exec(t))) { call = m[1]; continue; }
      if ((m = /^:\s*Hamiltonian\s+(\S+)\s*:/.exec(t))) { meta.method = m[1]; continue; }
      if ((m = /^:\s*net charge\s+(-?\d+)\s*:/.exec(t))) { meta.charge = int(m[1]); continue; }
      if ((m = /^:\s*unpaired electrons\s+(\d+)\s*:/.exec(t))) { meta.multiplicity = int(m[1]) + 1; continue; }
      if ((m = /^:\s*energy convergence\s+(\S+)\s+Eh\s*:/.exec(t))) { ethr = num(m[1]); continue; }
      if ((m = /^:\s*grad\. convergence\s+(\S+)\s+Eh/.exec(t))) { gthr = num(m[1]); continue; }
      if ((m = /^\.+\s*CYCLE\s+(\d+)\s*\.+$/.exec(t))) { cyc = { step: int(m[1]), energy: null, opt: {} }; steps.push(cyc); list = null; continue; }
      if (cyc && (m = /^\*\s*total energy\s*:\s*(\S+)\s+Eh\s+change\s+(\S+)\s+Eh/.exec(t))) { cyc.energy = num(m[1]); setOpt(cyc.opt, 'energyChange', num(m[2]), ethr); continue; }
      if (cyc && (m = /^gradient norm\s*:\s*(\S+)\s+Eh/.exec(t))) { setOpt(cyc.opt, 'gradNorm', num(m[1]), gthr); continue; }
      if (cyc && (m = /^maximum displ\.\s*:\s*(\S+)/.exec(t))) { setOpt(cyc.opt, 'maxDisp', num(m[1])); continue; }
      if (/GEOMETRY OPTIMIZATION CONVERGED/.test(t)) { conv = true; continue; }
      if (/FAILED TO CONVERGE GEOMETRY OPTIMIZATION/.test(t)) { conv = false; continue; }
      if (t === 'final structure:') {
        let j = i + 1;
        while (j < n && (/^=+$/.test(L[j].trim()) || !L[j].trim())) j++;
        let g = null;
        if (/^\$coord/.test(L[j].trim())) g = tmRows(L, j + 1);
        else if (isInt(L[j].trim())) {
          const nat = int(L[j].trim());
          const em = /energy:\s*(\S+)/.exec(L[j + 1] || '');
          if (em) fE = num(em[1]);
          g = readTable(L, j + 2, r => (r.length >= 4 && /^[A-Za-z]/.test(r[0]) && allNum(r, 1, 3) ? [elZ(r[0]), num(r[1]), num(r[2]), num(r[3])] : null));
          if (g.Z.length > nat) { g.Z.length = nat; g.xyz.length = 3 * nat; }
        }
        if (g && g.Z.length) F.add(g.Z, g.xyz, 'Final structure');
        i = g ? g.end - 1 : j;
        continue;
      }
      if ((m = /^\|\s*TOTAL ENERGY\s+(\S+)\s+Eh/.exec(t))) { fE = num(m[1]); continue; }
      if ((m = /^\|\s*TOTAL (ENTHALPY|FREE ENERGY)\s+(\S+)\s+Eh/.exec(t))) { thermo = thermo || newThermo(); thermo[m[1] === 'ENTHALPY' ? 'H' : 'G'] = num(m[2]); continue; }
      if ((m = /^::\s*(zero point energy|G\(RRHO\) contrib\.)\s+(\S+)\s+Eh/.exec(t))) { thermo = thermo || newThermo(); thermo[m[1][0] === 'z' ? 'zpe' : 'gCorr'] = num(m[2]); continue; }
      if (temp === null && (m = /^(\d+\.\d+)\s+VIB\b/.exec(t))) { temp = num(m[1]); continue; }
      if (/^#\s+Occupation\s+Energy\/Eh\s+Energy\/eV/.test(t)) {
        const e = [], occ = [];
        let j = i + 1;
        for (; j < n; j++) {
          const u = L[j].trim();
          if (/^-{10,}$/.test(u) || /\.\.\./.test(u)) continue;
          const r = tok(u.replace(/\((HOMO|LUMO)\)/, ''));
          if (r.length < 3 || !isInt(r[0]) || !r.every(isNum)) break;
          const k = int(r[0]) - 1;
          if (r.length >= 4) { occ[k] = num(r[1]); e[k] = num(r[2]); } else { occ[k] = 0; e[k] = num(r[1]); }
        }
        let h = -1;
        occ.forEach((o, k) => { if (o > 1e-3) h = k; });
        if (e.length) orb = { alpha: { energies: Array.from(e, v => (v === undefined ? NaN : v)), nOcc: h + 1 }, beta: null };
        i = j - 1;
        continue;
      }
      if (/^#\s+Z\s+covCN\s+q\s+C6AA/.test(t)) {
        const q = [];
        let j = i + 1;
        for (; j < n; j++) { const r = tok(L[j]); if (r.length < 5 || !isInt(r[0]) || !isInt(r[1]) || !isNum(r[4])) break; q.push(num(r[4])); }
        if (q.length) charges.Mulliken = q;
        i = j - 1;
        continue;
      }
      if (t === 'molecular dipole:') {
        for (let j = i + 1; j < Math.min(n, i + 5); j++) {
          const q = /^\s*full:\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)/.exec(L[j]);
          if (q) { ds.dipole = { x: num(q[1]) * AU2DEBYE, y: num(q[2]) * AU2DEBYE, z: num(q[3]) * AU2DEBYE, total: num(q[4]) }; break; }
        }
        continue;
      }
      if ((m = /^eigval\s*:\s*(.*)$/.exec(t))) { freqs.push(...tok(m[1]).map(num)); continue; }
      if (/^reduced masses \(amu\)/.test(t)) { list = red; continue; }
      if (/^IR intensities/.test(t)) { list = ir; continue; }
      if (/^(Raman intensities|output can be read|writing <)/.test(t)) { list = null; continue; }
      if (list && /^\d+:/.test(t)) { for (const [k, v] of pairsIn(t)) list[k] = v; continue; }
      if (/^normal termination of xtb/.test(t)) { normal = true; continue; }
      if (/abnormal termination of xtb/.test(t)) normal = false;
    }
    if (!ds.frames.length) {
      throw new Error('This xtb output has no coordinates (xtb prints them only after an optimisation); open the input structure or xtbopt.xyz.');
    }
    const f = F.cur();
    if (fE !== null) F.energy('xTB', fE, f);
    if (steps.length) {
      ds.optSteps = steps;
      f.opt = Object.assign({}, last(steps).opt);
      f.converged = conv === true;
    }
    if (/--(opt|ohess)\b/.test(call) || steps.length) addJT(meta, 'opt');
    if (/--(hess|ohess|bhess)\b/.test(call) || freqs.length) addJT(meta, 'freq');
    if (!meta.method) meta.method = /--gfnff\b/.test(call) ? 'GFN-FF' : 'GFN' + ((/--gfn\s+([012])/.exec(call) || [])[1] || '2') + '-xTB';
    if (meta.charge === null) { const c = /--(?:chrg|charge)\s+(-?\d+)/.exec(call); meta.charge = c ? int(c[1]) : 0; }
    if (meta.multiplicity === null) { const u = /--uhf\s+(\d+)/.exec(call); meta.multiplicity = u ? int(u[1]) + 1 : 1; }
    meta.route = call || null;
    const om = /--oniom\s+(\S+)(?:\s+([\d,-]+)(?=\s|$))?/.exec(call);
    if (om && ds.frames.length) {
      // xtb ONIOM: "--oniom gfn2:gfnff 1-5,8" (inner region, 1-based); the energies printed are the ONIOM totals
      const nm = x => ({ gfn2: 'GFN2-xTB', gfn1: 'GFN1-xTB', gfn0: 'GFN0-xTB', gfnff: 'GFN-FF', orca: 'ORCA', turbomole: 'Turbomole' }[String(x).toLowerCase()] || x || null);
      const lv = om[1].split(':');
      const f0 = ds.frames[ds.frames.length - 1], N = f0.Z.length;
      let inner = null;
      if (om[2]) {
        inner = new Set();
        for (const p of om[2].split(',')) { const r = /^(\d+)(?:-(\d+))?$/.exec(p); if (r) for (let k = int(r[1]); k <= int(r[2] || r[1]) && k <= N; k++) inner.add(k - 1); }
      }
      const links = [];
      if (inner) for (const a of inner) for (let b = 0; b < N; b++) {
        if (inner.has(b) || f0.Z[a] < 2 || f0.Z[b] < 2) continue;
        const d = Math.hypot(f0.xyz[3 * a] - f0.xyz[3 * b], f0.xyz[3 * a + 1] - f0.xyz[3 * b + 1], f0.xyz[3 * a + 2] - f0.xyz[3 * b + 2]);
        if (d < 1.15 * ((BM.el.covalent[f0.Z[a]] || 0.8) + (BM.el.covalent[f0.Z[b]] || 0.8))) links.push({ atom: b, bondedTo: a, element: 'H', type: null });
      }
      ds.layers = { scheme: 'oniom', atomLayer: inner ? Array.from({ length: N }, (_, a) => (inner.has(a) ? 'H' : 'L')) : null, frozen: null, active: null,
        linkAtoms: links, mmTypes: null, mmCharges: null, labels: { H: nm(lv[0]), M: null, L: nm(lv[1]) }, nLayers: 2, embedding: null };
      meta.qmmm = true;
      if (!inner) ds.warnings.push('The ONIOM inner region is given in a file (' + (om[2] || 'not shown') + '); layers are not shown.');
    }
    meta.normalTermination = normal;
    if (!ds.program) ds.program = 'xtb';
    if (Object.keys(charges).length) ds.charges = charges;
    ds.orbitals = orb;
    if (freqs.length) {
      ds.vibrations = buildVib(ds, ds.frames.length - 1, freqs, null, { ir, redMass: red }, true);
      if (ds.vibrations) ds.warnings.push('xtb prints no normal-mode vectors; open the g98.out file written by xtb to animate the modes.');
    }
    if (thermo) { if (thermo.temperature === null) thermo.temperature = temp; ds.thermo = thermoTotals(thermo, f.energy); }
    return finish(ds);
  }

  // =================================================================== MOPAC
  const MOPAC_SIG = /MOPAC20\d\d|\bMOPAC v\d|\*\*\s+MOPAC\b|Stewart Computational Chemistry|\bMOPAC\s+\d+\.\d+/;
  const MOPAC_METHOD = /^(PM7(-TS)?|PM6(-D3H4X?|-DH\+|-DH2X?|-D3)?|PM3|AM1|RM1|MNDOD?|PM5)$/i;
  const MULTS = ['', 'SINGLET', 'DOUBLET', 'TRIPLET', 'QUARTET', 'QUINTET', 'SEXTET', 'SEPTET', 'OCTET', 'NONET'];
  /** Input-style geometry lines: "Sym x f y f z f" (Cartesian) or "Sym r f a f d f na nb nc" (internal). */
  function mopacGeom(L, j) {
    const Z = [], vals = [];
    for (; j < L.length; j++) {
      const r = tok(L[j]);
      if (!r.length || !/^[A-Za-z]/.test(r[0])) break;
      const z = /^(XX|X|Tv)$/i.test(r[0]) ? 0 : elZ(r[0].replace(/\(.*$/, ''));
      const v = r.slice(1).filter(isNum).map(num);
      if (v.length < 3) break;
      Z.push(z); vals.push(v);
    }
    const internal = vals.some((v, k) => k > 0 && v.length >= 9 && v.slice(6, 9).some(x => x > 0));
    let xyz;
    if (internal) {
      const rows = vals.map((v, k) => {
        if (k === 0) return {};
        const a = v[6] > 0 ? v[6] - 1 : k - 1;
        if (k === 1) return { a, r: v[0] };
        const b = v[7] > 0 ? v[7] - 1 : (a === 0 ? 1 : 0);
        if (k === 2) return { a, r: v[0], b, angle: v[2] };
        return { a, r: v[0], b, angle: v[2], c: v[8] > 0 ? v[8] - 1 : 0, dihedral: v[4] };
      });
      xyz = Array.from(BM.zmatrixToCartesian(rows));
    } else xyz = [].concat(...vals.map(v => (v.length >= 5 ? [v[0], v[2], v[4]] : v.slice(0, 3))));
    return { Z, xyz, end: j };
  }
  function mopacKeywords(kw, meta) {
    const K = tok(kw.toUpperCase());
    meta.keywords = K;
    meta.route = K.join(' ') || null;
    const ch = K.map(k => /^CHARGE=([-+]?\d+)/.exec(k)).find(Boolean);
    if (ch) meta.charge = int(ch[1]);
    const mu = K.map(k => MULTS.indexOf(k)).find(k => k > 0);
    if (mu) meta.multiplicity = mu;
    meta.method = K.find(k => MOPAC_METHOD.test(k)) || meta.method || 'PM7';
    const has = re => K.some(k => re.test(k));
    if (has(/^(FORCE|FREQCY|THERMO)/)) addJT(meta, 'freq');
    if (has(/^IRC/)) addJT(meta, 'irc');
    if (has(/^TS$/)) { addJT(meta, 'opt'); addJT(meta, 'ts'); }
    if (has(/^OPT/) || !(has(/^1SCF$/) || has(/^(FORCE|FREQCY|IRC|DRC|THERMO)/))) addJT(meta, 'opt');
    return K;
  }
  function parseMopac(text, name) {
    const ds = BM.emptyDataset('mopac-out', name), meta = ds.meta, F = new Frames(ds), L = lines(text), n = L.length;
    const steps = [], charges = {}, vf = [], vt = [], vr = [], vs = [];
    let kwPending = false, kw = null, hof = null, hofFrame = null, gthr = 1, normal = null, conv = false, nFilled = null;
    let eig = null, modes = null, fi = -1;
    for (let i = 0; i < n; i++) {
      const t = L[i].trim();
      if (!t) continue;
      let m;
      if ((m = /(MOPAC20\d\d)|Version:\s*([\w.]+)|\bMOPAC\s+v?(\d[\w.]*)/.exec(t)) && (!ds.program || /Version/.test(t))) {
        ds.program = 'MOPAC' + (m[1] ? ' ' + m[1].slice(5) : '') + (m[2] ? ' ' + m[2] : m[3] ? ' ' + m[3] : '');
      }
      if (/^\*\s+CALCULATION DONE:/.test(t)) { kwPending = true; continue; }
      if (kwPending && /^\*{40,}$/.test(t)) { kwPending = false; if (kw === null) { kw = L[i + 1] || ''; mopacKeywords(kw, meta); ds.title = (L[i + 2] || '').trim() || ds.title; } continue; }
      if (t === 'FINAL GEOMETRY OBTAINED') {
        if (kw === null) { kw = L[i + 1] || ''; mopacKeywords(kw, meta); }
        const g = mopacGeom(L, i + 4);
        if (g.Z.length) { const f = F.add(g.Z, g.xyz, 'Final geometry', 1e-4); if (hof !== null && hofFrame === null) { F.energy('HOF', hof / KCAL, f); hofFrame = f; } }
        i = g.end - 1;
        continue;
      }
      if (t === 'CARTESIAN COORDINATES') {
        let j = i + 1;
        while (j < n && j < i + 6 && !/^\s*\d+\s+[A-Za-z]{1,2}\s+-?\d/.test(L[j])) j++;
        const g = readTable(L, j, r => (r.length >= 5 && isInt(r[0]) && /^[A-Za-z]/.test(r[1]) && allNum(r, 2, 3) ? [elZ(r[1]), num(r[2]), num(r[3]), num(r[4])] : null));
        if (g.Z.length) {
          const f = F.add(g.Z, g.xyz, ds.frames.length ? 'Final geometry' : 'Input geometry', 1e-4);
          if (hof !== null && hofFrame === null) { F.energy('HOF', hof / KCAL, f); hofFrame = f; }
        }
        i = g.end - 1;
        continue;
      }
      if ((m = /^GNORM=\s*(\S+)/i.exec(t))) { gthr = num(m[1]); continue; }
      if ((m = /CYCLE:\s*(\d+).*?GRAD\.:\s*(\S+)\s+HEAT:\s*(\S+)/.exec(t))) {
        const o = {};
        setOpt(o, 'gradNorm', num(m[2]), gthr);
        steps.push({ step: int(m[1]), energy: num(m[3]) / KCAL, opt: o });
        continue;
      }
      if (/GRADIENT\s+=\s*\S+\s+IS LESS THAN CUTOFF|GEOMETRY OPTIMISED|HEAT OF FORMATION TEST SATISFIED|GRADIENT TEST PASSED/.test(t)) { conv = true; continue; }
      if ((m = /^(?:FINAL )?HEAT OF FORMATION\s*=\s*(\S+)\s+KCAL/.exec(t))) { hof = num(m[1]); meta.heatOfFormation = hof; hofFrame = null; continue; }
      if ((m = /^NO\. OF FILLED LEVELS\s*=\s*(\d+)/.exec(t))) { nFilled = int(m[1]); continue; }
      if ((m = /^DIPOLE\s*=\s*(\S+)\s+DEBYE/.exec(t))) { if (!ds.dipole) ds.dipole = { x: null, y: null, z: null, total: num(m[1]) }; continue; }
      if ((m = /^SUM\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)/.exec(t))) { ds.dipole = { x: num(m[1]), y: num(m[2]), z: num(m[3]), total: num(m[4]) }; continue; }
      if (t === 'EIGENVALUES') {
        const e = [];
        let j = i + 1;
        for (; j < n; j++) { const r = tok(L[j]); if (!r.length) { if (e.length) break; continue; } if (!r.every(isNum)) break; e.push(...r.map(x => num(x) / HARTREE2EV)); }
        if (e.length) eig = e;
        i = j - 1;
        continue;
      }
      if (/^NET ATOMIC CHARGES/.test(t)) {
        let j = i + 1;
        while (j < n && j < i + 4 && !/ATOM NO\./.test(L[j])) j++;
        const q = [];
        for (j++; j < n; j++) { const r = tok(L[j]); if (r.length < 3 || !isInt(r[0]) || !isNum(r[2])) break; q.push(num(r[2])); }
        if (q.length) charges.Coulson = q;
        i = j - 1;
        continue;
      }
      if ((m = /^VIBRATION\s+(\d+)\s*(\S*)/.exec(t))) { if (fi < 0) fi = ds.frames.length - 1; vs[int(m[1]) - 1] = /^\d*[A-Za-z]/.test(m[2]) ? m[2].replace(/^\d+/, '') : ''; continue; }
      if ((m = /^FREQ\.\s+(\S+)/.exec(t))) { vf.push(num(m[1])); continue; }
      if ((m = /^T-DIPOLE\s+(\S+)/.exec(t))) { vt[vf.length - 1] = num(m[1]); continue; }
      if ((m = /^RED\. MASS\s+(\S+)/.exec(t))) { vr[vf.length - 1] = num(m[1]); continue; }
      if (/^NORMAL COORDINATE ANALYSIS/.test(t)) {
        modes = [];
        let cols = null, j = i + 1, gotF = false;
        for (; j < n; j++) {
          const u = L[j].trim();
          if (!u) continue;
          const r = u.split(/\s+/);
          if (/^Root No\./i.test(u)) { cols = r.slice(2).map(x => int(x) - 1); gotF = false; continue; }
          if (!cols) break;
          if (!gotF) { if (r.length === cols.length && r.every(isNum)) gotF = true; continue; }
          if (r.length === cols.length + 1 && isInt(r[0]) && allNum(r, 1, cols.length)) {
            cols.forEach((md, k) => { (modes[md] = modes[md] || [])[int(r[0]) - 1] = num(r[k + 1]); });
            continue;
          }
          break;
        }
        i = j - 1;
        continue;
      }
      if (/==\s*MOPAC DONE\s*==/.test(t)) { normal = true; continue; }
      if (/^(ERROR|JOB ENDED ABNORMALLY)/.test(t)) normal = false;
    }
    if (hof !== null && hofFrame === null && F.cur()) F.energy('HOF', hof / KCAL);
    if (steps.length) {
      ds.optSteps = steps;
      const f = F.cur();
      if (f) { f.opt = Object.assign({}, last(steps).opt); f.converged = conv; }
    }
    meta.normalTermination = normal;
    if (!ds.program) ds.program = 'MOPAC';
    if (Object.keys(charges).length) ds.charges = charges;
    if (eig) ds.orbitals = { alpha: { energies: eig, nOcc: nFilled === null ? 0 : nFilled }, beta: null };
    if (vf.length) {
      ds.vibrations = buildVib(ds, fi < 0 ? ds.frames.length - 1 : fi, vf, modes,
        { ir: vt.map(x => x * x), redMass: vr, sym: vs }, false);
      addJT(meta, 'freq');
      if (ds.vibrations && ds.vibrations.ir) ds.warnings.push('MOPAC IR intensities are relative (squared transition dipoles, Debye^2).');
    }
    return finish(ds, 'No CARTESIAN COORDINATES or FINAL GEOMETRY OBTAINED block found in this MOPAC file.');
  }
  function parseMopacInput(text, name) {
    const ds = BM.emptyDataset('mopac-input', name), meta = ds.meta, L = lines(text);
    let i = 0;
    while (i < L.length && /^\s*\*/.test(L[i])) i++;   // comment lines
    let kw = L[i] || '';
    let titles = 2;
    if (/(^|\s)\+\s*$/.test(kw)) { kw = kw.replace(/\+\s*$/, ' ') + (L[++i] || ''); }
    else if (/(^|\s)&\s*$/.test(kw)) { kw = kw.replace(/&\s*$/, ' ') + (L[++i] || ''); titles = 1; }
    mopacKeywords(kw, meta);
    ds.title = (L[i + 1] || '').trim() || ds.title;
    const g = mopacGeom(L, i + 1 + titles);
    if (!g.Z.length) throw new Error('No geometry found after the keyword and title lines of this MOPAC input.');
    const keep = g.Z.map((z, k) => k).filter(k => g.Z[k] > 0 || !/^Tv/i.test(tok(L[i + 1 + titles + k])[0] || ''));
    ds.frames.push(mkFrame(keep.map(k => g.Z[k]), [].concat(...keep.map(k => g.xyz.slice(3 * k, 3 * k + 3)))));
    ds.program = null;
    return finish(ds);
  }

  // =================================================================== CP2K
  const CP2K_SIG = /CP2K\|\s+version string:|CP2K\|\s+source code revision/;
  function parseCp2k(text, name) {
    const ds = BM.emptyDataset('cp2k-out', name), meta = ds.meta, F = new Frames(ds), L = lines(text), n = L.length;
    const steps = [], charges = {}, cell = [], vf = [], vir = [], vred = [], vmodes = [];
    let run = null, normal = null, conv = false, e0 = null, vfi = -1, qmmm = false;
    const modOf = new Map();   // frame -> module that printed it (FIST: whole system, QUICKSTEP: QM region in QM/MM)
    for (let i = 0; i < n; i++) {
      const t = L[i].trim();
      if (!t) continue;
      let m;
      if ((m = /^CP2K\|\s+version string:\s+(.*)$/.exec(t))) { ds.program = m[1].trim().replace(/^CP2K version/i, 'CP2K'); continue; }
      if ((m = /^GLOBAL\|\s+Run type\s+(\S+)/.exec(t))) { run = m[1].toUpperCase(); continue; }
      if ((m = /^DFT\|\s+Charge\s+(-?\d+)/.exec(t))) { meta.charge = int(m[1]); continue; }
      if ((m = /^DFT\|\s+Multiplicity\s+(\d+)/.exec(t))) { meta.multiplicity = int(m[1]); continue; }
      if ((m = /^FUNCTIONAL\|\s+([^:\s][^:]*?):/.exec(t)) && !meta.method) { meta.method = m[1].trim(); continue; }
      if ((m = /^Orbital Basis Set\s+(\S+)/.exec(t)) && !meta.basis) { meta.basis = m[1]; continue; }
      if ((m = /^CELL\|\s+Vector ([abc]) \[angstrom\]:\s+(\S+)\s+(\S+)\s+(\S+)/.exec(t))) { cell['abc'.indexOf(m[1])] = [num(m[2]), num(m[3]), num(m[4])]; continue; }
      if (t.startsWith('QMMM|')) { qmmm = true; continue; }
      if ((m = /ATOMIC COORDINATES IN\s+(angstrom|bohr)/i.exec(t))) {
        const scale = /bohr/i.test(m[1]) ? BOHR : 1;
        const mod = ((/MODULE\s+([A-Z_]+)/i.exec(t) || [])[1] || '').toUpperCase();
        let j = i + 1;
        while (j < n && j < i + 6 && !/^\s*\d+\s+\d+\s+[A-Za-z]/.test(L[j])) j++;
        const g = readTable(L, j, r => (r.length >= 7 && isInt(r[0]) && isInt(r[1]) && isInt(r[3]) && allNum(r, 4, 3)
          ? [int(r[3]) > 0 ? int(r[3]) : elZ(r[2]), num(r[4]), num(r[5]), num(r[6])] : null), scale);
        if (g.Z.length) modOf.set(F.add(g.Z, g.xyz, ds.frames.length ? null : 'Input geometry'), mod);
        i = g.end - 1;
        continue;
      }
      if ((m = /^ENERGY\|\s+Total FORCE_EVAL \(\s*(\w+)\s*\) energy \[(?:a\.u\.|hartree)\]:\s+(\S+)/.exec(t))) {
        if (/^QMMM$/i.test(m[1])) qmmm = true;
        if (e0 === null) { e0 = num(m[2]); F.energy(qmmm ? 'QM/MM' : 'SCF', e0); }
        continue;
      }
      if ((m = /^-+\s*Informations at step\s*=\s*(\d+)/.exec(t))) {
        const st = { step: int(m[1]), energy: null, opt: {} }, v = {};
        let j = i + 1;
        for (; j < n && !/^\s*-{20,}\s*$/.test(L[j]); j++) {
          const q = /^\s*(.+?)\s*=\s*(\S+)\s*$/.exec(L[j]);
          if (q) v[q[1]] = q[2];
        }
        const yes = k => (v[k] === undefined ? undefined : /YES/i.test(v[k]));
        st.energy = num(v['Total Energy']);
        setOpt(st.opt, 'energyChange', num(v['Real energy change']));
        setOpt(st.opt, 'maxDisp', num(v['Max. step size']), num(v['Conv. limit for step size']), yes('Convergence in step size'));
        setOpt(st.opt, 'rmsDisp', num(v['RMS step size']), num(v['Conv. limit for RMS step']), yes('Convergence in RMS step'));
        setOpt(st.opt, 'maxForce', num(v['Max. gradient']), num(v['Conv. limit for gradients']), yes('Conv. for gradients') !== undefined ? yes('Conv. for gradients') : yes('Conv. in gradients'));
        setOpt(st.opt, 'rmsForce', num(v['RMS gradient']), num(v['Conv. limit for RMS grad.']), yes('Conv. in RMS gradients'));
        steps.push(st);
        i = j;
        continue;
      }
      if (/GEOMETRY OPTIMIZATION COMPLETED/.test(t)) { conv = true; continue; }
      if (t === 'Mulliken Population Analysis') {
        let j = i + 1;
        while (j < n && j < i + 6 && !/^\s*#\s*Atom\s+Element/.test(L[j])) j++;
        const q = [];
        for (j++; j < n; j++) { const r = tok(L[j]); if (r.length < 5 || !isInt(r[0]) || !r.slice(3).every(isNum)) break; q.push(num(r.length >= 7 ? r[5] : r[4])); }
        if (q.length) charges.Mulliken = q;
        i = j - 1;
        continue;
      }
      if ((m = /^VIB\|\s*Frequency \(cm\^-1\)\s+(.*)$/.exec(t))) { if (vfi < 0) vfi = ds.frames.length - 1; vf.push(...tok(m[1]).map(num)); continue; }
      if ((m = /^VIB\|\s*(?:IR int \(KM\/Mole\)|Intensities)\s+(.*)$/.exec(t))) { vir.push(...tok(m[1]).map(num)); continue; }
      if ((m = /^VIB\|\s*Red\.\s*Masses \(a\.u\.\)\s+(.*)$/.exec(t))) { vred.push(...tok(m[1]).map(num)); continue; }
      if (/^ATOM\s+EL\s+X\s+Y\s+Z/.test(t)) {
        const nb = tok(t).filter(x => x === 'X').length, base = vf.length - nb;
        let j = i + 1;
        for (; j < n; j++) {
          const r = tok(L[j]);
          if (r.length !== 2 + 3 * nb || !isInt(r[0]) || !allNum(r, 2, 3 * nb)) break;
          const a = int(r[0]) - 1;
          for (let k = 0; k < nb; k++) { const v = vmodes[base + k] = vmodes[base + k] || []; for (let c = 0; c < 3; c++) v[3 * a + c] = num(r[2 + 3 * k + c]); }
        }
        i = j - 1;
        continue;
      }
      if (/PROGRAM ENDED AT/.test(t)) { normal = true; continue; }
      if (/\[ABORT\]|CP2K\|\s+ABORT/.test(t)) normal = false;
    }
    const mods = new Set(Array.from(modOf.values()));
    if (qmmm || (mods.has('FIST') && (mods.has('QUICKSTEP') || mods.has('QS')))) {
      // QM/MM: FIST prints the whole system, QUICKSTEP only the QM region (+ link atoms, shifted into the QM box)
      const full = ds.frames.filter(f => modOf.get(f) === 'FIST'), qs = ds.frames.filter(f => modOf.get(f) !== 'FIST');
      if (full.length && qs.length) {
        for (const f of qs) if (f.energies && !full[0].energies) { full[0].energies = f.energies; full[0].energy = f.energy; }
        full[0].label = full[0].label || 'Input geometry';
        ds.frames = full;
        if (vfi >= 0) vfi = 0;
        const r = BM.qmmmMatch ? BM.qmmmMatch(full[0].Z, full[0].xyz, qs[0].Z, qs[0].xyz, { shift: true }) : null;
        const N = full[0].Z.length, inQM = new Set(r ? r.qm : []);
        ds.layers = { scheme: 'qmmm', atomLayer: r ? Array.from({ length: N }, (_, a) => (inQM.has(a) ? 'H' : 'L')) : null, frozen: null, active: null,
          linkAtoms: r ? r.links : [], mmTypes: null, mmCharges: null,
          labels: { H: [meta.method, meta.basis].filter(Boolean).join('/') || 'QM', M: null, L: 'MM (FIST)' }, nLayers: 2, embedding: null };
        if (!r) ds.warnings.push('The QM atoms could not be located in the MM coordinates.');
      }
      meta.qmmm = true;
    }
    if (run === 'GEO_OPT' || run === 'GEOMETRY_OPTIMIZATION' || run === 'CELL_OPT') addJT(meta, 'opt');
    if (run === 'VIBRATIONAL_ANALYSIS' || vf.length) addJT(meta, 'freq');
    if (steps.length) {
      ds.optSteps = steps;
      if (conv) steps[steps.length - 1].converged = true;
      ds.warnings.push('CP2K prints only the starting coordinates; the optimised geometries are in the *-pos-1.xyz trajectory.');
    }
    meta.route = run;
    meta.normalTermination = normal;
    if (!ds.program) ds.program = 'CP2K';
    if (cell.length === 3 && cell.every(Boolean)) ds.cell = cell;
    if (Object.keys(charges).length) ds.charges = charges;
    if (vf.length) ds.vibrations = buildVib(ds, vfi < 0 ? 0 : vfi, vf, vmodes, { ir: vir, redMass: vred }, false);
    return finish(ds, 'No "ATOMIC COORDINATES" block found in this CP2K output.');
  }
  /** CP2K input: &COORD (Angstrom unless UNIT bohr; SCALED with &CELL) and &CELL (ABC or A/B/C vectors). */
  function parseCp2kInput(text, name) {
    const ds = BM.emptyDataset('cp2k-input', name), meta = ds.meta, L = lines(text).map(s => s.replace(/[!#].*$/, ''));
    const Z = [], X = [];
    let cell = null, unit = 1, scaled = false, inCoord = false, inCell = false, abc = null, ang = [90, 90, 90];
    const vec = {};
    for (const s of L) {
      const t = s.trim();
      if (!t) continue;
      if (/^&COORD\b/i.test(t)) { inCoord = true; continue; }
      if (/^&END\s+COORD\b/i.test(t) || (inCoord && /^&END\s*$/i.test(t))) { inCoord = false; continue; }
      if (/^&CELL\b/i.test(t)) { inCell = true; continue; }
      if (/^&END\s+CELL\b/i.test(t) || (inCell && /^&END\s*$/i.test(t))) { inCell = false; continue; }
      const r = t.split(/\s+/);
      if (inCoord) {
        if (/^UNIT$/i.test(r[0])) { unit = /^bohr$/i.test(r[1] || '') ? BOHR : 1; continue; }
        if (/^SCALED$/i.test(r[0])) { scaled = !/^(\.?F|\.FALSE\.|FALSE|NO)/i.test(r[1] || 'T'); continue; }
        if (r.length >= 4 && /^[A-Za-z]/.test(r[0]) && allNum(r, 1, 3)) { Z.push(elZ(r[0].replace(/[_\d].*$/, ''))); X.push(num(r[1]), num(r[2]), num(r[3])); }
        continue;
      }
      if (inCell) {
        if (/^ABC$/i.test(r[0]) && allNum(r, r.length > 4 ? 2 : 1, 3)) { const o = r.length > 4 && /^\[/.test(r[1]) ? 2 : 1; abc = r.slice(o, o + 3).map(num); }
        else if (/^ALPHA_BETA_GAMMA$/i.test(r[0])) ang = r.slice(1, 4).map(num);
        else if (/^[ABC]$/i.test(r[0]) && allNum(r, 1, 3)) vec[r[0].toUpperCase()] = r.slice(1, 4).map(num);
      }
    }
    if (vec.A && vec.B && vec.C) cell = [vec.A, vec.B, vec.C];
    else if (abc) {
      const [al, be, ga] = ang.map(x => x * Math.PI / 180);
      const cx = Math.cos(be), cy = (Math.cos(al) - Math.cos(be) * Math.cos(ga)) / Math.sin(ga);
      cell = [[abc[0], 0, 0], [abc[1] * Math.cos(ga), abc[1] * Math.sin(ga), 0], [abc[2] * cx, abc[2] * cy, abc[2] * Math.sqrt(Math.max(0, 1 - cx * cx - cy * cy))]]
        .map(v => v.map(x => (Math.abs(x) < 1e-9 ? 0 : x)));
    }
    if (!Z.length) throw new Error('No &COORD section found in this CP2K input (coordinates may come from an external file).');
    const xyz = new Array(X.length);
    for (let a = 0; a < Z.length; a++) {
      const p = X.slice(3 * a, 3 * a + 3);
      for (let c = 0; c < 3; c++) xyz[3 * a + c] = scaled && cell ? p[0] * cell[0][c] + p[1] * cell[1][c] + p[2] * cell[2][c] : p[c] * unit;
    }
    ds.frames.push(mkFrame(Z, xyz));
    ds.cell = cell;
    ds.program = null;
    const g = /^\s*RUN_TYPE\s+(\S+)/im.exec(L.join('\n'));
    if (g && /GEO_OPT|CELL_OPT/i.test(g[1])) addJT(meta, 'opt');
    if (g && /VIBRATIONAL/i.test(g[1])) addJT(meta, 'freq');
    const c = /^\s*CHARGE\s+(-?\d+)/im.exec(L.join('\n')), mu = /^\s*MULTIPLICITY\s+(\d+)/im.exec(L.join('\n'));
    if (c) meta.charge = int(c[1]);
    if (mu) meta.multiplicity = int(mu[1]);
    const xc = /&XC_FUNCTIONAL\s+(\w+)/i.exec(L.join('\n'));
    if (xc) meta.method = xc[1].toUpperCase();
    return finish(ds);
  }

  // =================================================================== Turbomole
  const tmEl = s => elZ(s ? s[0].toUpperCase() + s.slice(1).toLowerCase() : s);
  /** $coord rows "x y z el [f]" in Bohr, up to the next $ line. */
  function tmRows(L, j) {
    const Z = [], xyz = [];
    for (; j < L.length; j++) {
      const t = L[j].trim();
      if (!t) continue;
      if (t[0] === '$') break;
      const r = t.split(/\s+/);
      if (r.length < 4 || !allNum(r, 0, 3) || !/^[A-Za-z]/.test(r[3])) break;
      Z.push(tmEl(r[3])); xyz.push(num(r[0]) * BOHR, num(r[1]) * BOHR, num(r[2]) * BOHR);
    }
    return { Z, xyz, end: j };
  }
  function parseTmCoord(text, name) {
    const ds = BM.emptyDataset('turbomole-coord', name), L = lines(text);
    const i = L.findIndex(s => /^\s*\$coord\b/.test(s));
    if (i < 0) throw new Error('No $coord block found.');
    const fm = /file\s*=\s*(\S+)/.exec(L[i]);
    if (fm) throw new Error('$coord refers to the file "' + fm[1] + '"; open that file instead.');
    if (/\bfrac\b/i.test(L[i])) throw new Error('Fractional $coord (riper) is not read.');
    const g = tmRows(L, i + 1);
    if (!g.Z.length) throw new Error('The $coord block is empty.');
    ds.frames.push(mkFrame(g.Z, g.xyz));
    ds.program = null;
    return finish(ds);
  }
  function parseTmGrad(text, name) {
    const ds = BM.emptyDataset('turbomole-gradient', name), meta = ds.meta, F = new Frames(ds), L = lines(text), n = L.length;
    for (let i = 0; i < n; i++) {
      const m = /cycle\s*=\s*(\d+)\s+(.*?)energy\s*=\s*(\S+)(?:\s+\|dE\/dxyz\|\s*=\s*(\S+))?/i.exec(L[i]);
      if (!m) continue;
      const g = tmRows(L, i + 1);
      const N = g.Z.length, grad = [];
      let j = g.end;
      for (; j < n && grad.length < 3 * N; j++) {
        const r = tok(L[j]);
        if (r.length !== 3 || !r.every(isNum)) break;
        grad.push(...r.map(num));
      }
      if (!N || grad.length < 3 * N) { i = j - 1; continue; }
      const prev = F.cur();
      const f = mkFrame(g.Z, g.xyz);
      ds.frames.push(f);
      f.label = 'Cycle ' + m[1];
      const kind = /mp2/i.test(m[2]) ? 'MP2' : /cc/i.test(m[2]) ? 'CCSD' : 'SCF';
      F.energy(kind, num(m[3]), f);
      const o = f.opt = {};
      let mx = 0, ss = 0;
      for (const v of grad) { mx = Math.max(mx, Math.abs(v)); ss += v * v; }
      setOpt(o, 'maxForce', mx);
      setOpt(o, 'rmsForce', Math.sqrt(ss / grad.length));
      setOpt(o, 'gradNorm', m[4] !== undefined ? num(m[4]) : Math.sqrt(ss));
      if (prev && Number.isFinite(prev.energy)) setOpt(o, 'energyChange', f.energy - prev.energy);
      i = j - 1;
    }
    if (ds.frames.length) addJT(meta, 'opt');
    ds.program = 'Turbomole';
    return finish(ds, 'No "cycle = ..." blocks with coordinates and gradients found in this Turbomole gradient file.');
  }

  // =================================================================== registration
  const mopacInpSniff = h => {
    const first = (lines(h).find(s => s.trim() && !/^\s*\*/.test(s)) || '').toUpperCase();
    return tok(first).some(k => MOPAC_METHOD.test(k) || /^(1SCF|CHARGE=|XYZ|PRECISE|GNORM=)/.test(k)) ? 0.5 : 0;
  };
  reg('gamess-out', 'GAMESS-US / Firefly output', ['log', 'out', 'gamout'], 'output', h => (GAMESS_SIG.test(h) ? 0.95 : 0), parseGamess);
  reg('gamess-input', 'GAMESS input', ['inp'], 'structure',
    h => (!GAMESS_SIG.test(h) && /^\s*\$(CONTRL|BASIS|SYSTEM|SCF)\b/im.test(h) && /^\s*\$DATA\b/im.test(h) ? 0.6 : 0), parseGamessInput);
  reg('nwchem-out', 'NWChem output', ['out', 'log', 'nwout'], 'output', h => (NWCHEM_SIG.test(h) ? 0.95 : 0), parseNwchem);
  reg('nwchem-input', 'NWChem input', ['nw', 'nwi'], 'structure',
    h => (!NWCHEM_SIG.test(h) && /^\s*geometry\b/im.test(h) && /^\s*task\s+\w+/im.test(h) ? 0.6 : 0), parseNwInput);
  reg('qchem-out', 'Q-Chem output', ['out', 'log', 'qcout'], 'output', h => (QCHEM_SIG.test(h) ? 0.95 : 0), parseQchem);
  reg('qchem-input', 'Q-Chem input', ['in', 'qcin', 'inp'], 'structure',
    h => (!QCHEM_SIG.test(h) && /^\s*\$molecule\b/im.test(h) && /^\s*\$rem\b/im.test(h) ? 0.6 : 0), parseQchemInput);
  reg('psi4-out', 'Psi4 output', ['out', 'log'], 'output', h => (PSI4_SIG.test(h) ? 0.95 : 0), parsePsi4);
  reg('psi4-input', 'Psi4 input', ['dat', 'in', 'psi4'], 'structure',
    h => (!PSI4_SIG.test(h) && /^\s*molecule\s*\w*\s*\{/m.test(h) && /\b(energy|optimize|opt|frequency|frequencies|freq|properties|gradient|hessian)\s*\(/.test(h) ? 0.6 : 0), parsePsi4Input);
  reg('xtb-out', 'xtb output', ['out', 'log', 'xtbout'], 'output', h => (XTB_SIG.test(h) ? 0.95 : 0), parseXtb);
  reg('mopac-out', 'MOPAC output / archive', ['out', 'arc', 'mopout'], 'output', h => (MOPAC_SIG.test(h) ? 0.95 : 0), parseMopac);
  reg('mopac-input', 'MOPAC input', ['mop'], 'structure', h => (MOPAC_SIG.test(h) ? 0 : mopacInpSniff(h)), parseMopacInput);
  reg('cp2k-out', 'CP2K output', ['out', 'log'], 'output', h => (CP2K_SIG.test(h) ? 0.95 : 0), parseCp2k);
  reg('cp2k-input', 'CP2K input', ['inp', 'restart'], 'structure',
    h => (!CP2K_SIG.test(h) && /^\s*&(FORCE_EVAL|GLOBAL)\b/im.test(h) ? 0.6 : 0), parseCp2kInput);
  reg('turbomole-coord', 'Turbomole coord (Bohr)', ['coord', 'tmol'], 'structure',
    (h, nm) => (/^\s*\$coord\b/m.test(h) && !/^\s*\$grad\b/m.test(h) ? (nm === 'coord' ? 0.9 : 0.8) : 0), parseTmCoord);
  reg('turbomole-gradient', 'Turbomole gradient', ['gradient', 'grad'], 'output',
    h => (/^\s*\$grad\b/m.test(h) && /cycle\s*=/.test(h) ? 0.9 : 0), parseTmGrad);

  BM.qcReaders = { gamess: parseGamess, gamessInput: parseGamessInput, nwchem: parseNwchem, nwchemInput: parseNwInput,
    qchem: parseQchem, qchemInput: parseQchemInput, psi4: parsePsi4, psi4Input: parsePsi4Input, xtb: parseXtb,
    mopac: parseMopac, mopacInput: parseMopacInput, cp2k: parseCp2k, cp2kInput: parseCp2kInput,
    turbomoleCoord: parseTmCoord, turbomoleGradient: parseTmGrad };
})(typeof window !== 'undefined' ? (window.BM = window.BM || {}) : (globalThis.BM = globalThis.BM || {}));
