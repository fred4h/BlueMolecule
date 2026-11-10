/* BlueMolecule — crystal, materials-science and molecular-dynamics structure readers.
 *
 * Registers with BM.registerFormat (Dataset shape: see parsers.js):
 *   cif, mmcif                CIF (cell, fractional sites, symmetry expansion) and mmCIF/PDBx atom_site
 *   mol2, cml, cjson          Tripos Mol2, Chemical Markup Language, Avogadro CJSON
 *   poscar, xdatcar, outcar   VASP
 *   qe-in, qe-out             Quantum ESPRESSO pw.x input and output
 *   castep-cell, xsf, cp2k-in CASTEP .cell, XCrySDen XSF/AXSF (with DATAGRID_3D), CP2K input
 *   prmtop, inpcrd, psf, dcd, lammps-data, lammps-dump, trr, xtc   MD topologies and trajectories
 * Units: Å and Hartree. frame.opt.maxForce / rmsForce are in Hartree/Bohr (max = largest atomic force,
 * rms over all Cartesian components). Extras beyond the core Dataset: frame.chains, frame.cell (only for
 * variable-cell runs), frame.occupancy (CIF), frame.formalCharges (CML, CJSON), frame.time (ps),
 * ds.cellOrigin (LAMMPS box origin), ds.meta.spaceGroup (CIF).
 * Trajectory-only files: meta.needsTopology + meta.natoms, frames with Z = null. DCD / TRR / XTC also read from Blob slices
 * (parseSource, BM.structReaders.sources): every Nth frame of multi-GB files, meta.totalFrames / stride / reload(opts).
 * Topology-only files: meta.topologyOnly, one frame with Z and zero coordinates, ds.bonds.
 * Written from public descriptions of the file formats only.
 */
(function (BM) {
  'use strict';
  const BOHR = 0.529177210903, HARTREE_EV = 27.211386245988, RAD = Math.PI / 180;
  const EVA_AU = BOHR / HARTREE_EV;            // eV/Å -> Hartree/Bohr
  const AMBER_Q = 18.2223;                     // prmtop charge units per elementary charge

  // ------------------------------------------------------------------ helpers
  function num(s) {
    if (typeof s === 'number') return s;
    if (s === undefined || s === null) return NaN;
    s = String(s).trim();
    return s ? Number(s.replace(/[dD]/, 'e')) : NaN;
  }
  const splitLines = t => (Array.isArray(t) ? t : String(t).split(/\r\n|\r|\n/));   // text or a line array (big files)
  /** Whole text for readers that need it; a line array is joined when that is safe (BM.linesToText). */
  const wholeText = (t, what) => (!Array.isArray(t) ? String(t) : BM.linesToText ? BM.linesToText(t, what) : t.join('\n'));
  const tok = s => { s = String(s === undefined || s === null ? '' : s).trim(); return s ? s.split(/\s+/) : []; };
  const baseName = f => (f ? String(f).split(/[\\/]/).pop() : '');
  const stripExt = f => baseName(f).replace(/\.[^.]*$/, '');
  const err = m => new Error(m);
  const isNum = s => /^[-+]?(\d+\.?\d*|\.\d+)([eEdD][-+]?\d+)?$/.test(s);
  const isInt = s => /^[-+]?\d+$/.test(s);
  const fin = a => a.every(Number.isFinite);
  const sum = a => a.reduce((s, x) => s + x, 0);

  function newFrame(Z, xyz) {
    return { Z: Z ? Array.from(Z) : null, xyz: xyz instanceof Float64Array ? xyz : Float64Array.from(xyz),
      energy: null, energies: null, label: null, converged: false, opt: null,
      atomNames: null, resNames: null, resIds: null, chains: null };
  }
  function newDs(id, fileName, program) {
    const ds = BM.emptyDataset(id, fileName);
    if (program) ds.program = program;
    return ds;
  }
  /** Z from an element symbol or label ("Na1+", "Fe_pv", "C12"); D/T -> 1; -1 if unknown. */
  function zOf(label) {
    if (label === undefined || label === null) return -1;
    const s = String(label).trim();
    if (/^[DT](?:$|[^A-Za-z])/.test(s)) return 1;
    return BM.el.toZ(s);
  }
  /** Z from a PDB-style atom name; -1 if unknown. */
  function zFromName(name, resName) {
    const g = BM.elementFromAtomName ? BM.elementFromAtomName(String(name || ''), resName) : null;
    if (!g) return -1;
    return g === 'D' || g === 'T' ? 1 : BM.el.toZ(g);
  }
  const UNITED = [[13.019, 6], [14.027, 6], [15.035, 6]];   // united-atom CH, CH2, CH3
  /** Element with the nearest standard atomic mass (H for < 3.5, covering D and repartitioned H). */
  function zFromMass(m) {
    if (!(m > 0.5)) return 0;
    if (m < 3.5) return 1;
    const M = BM.el.mass;
    let best = 0, bd = Infinity;
    for (let z = 2; z < M.length; z++) { if (!(M[z] > 0)) continue; const d = Math.abs(M[z] - m); if (d < bd) { bd = d; best = z; } }
    for (const u of UNITED) if (Math.abs(u[0] - m) < bd) { bd = Math.abs(u[0] - m); best = u[1]; }
    return best;
  }
  /** Mass first; the atom name decides when the mass is not a standard one (mass repartitioning). */
  function zFromMassName(m, name, res) {
    const z = zFromMass(m);
    if (z > 1 && Math.abs(BM.el.mass[z] - m) > 0.1 && !UNITED.some(u => Math.abs(u[0] - m) < 0.01)) {
      const g = zFromName(name, res);
      if (g > 0) return g;
    }
    return z;
  }

  // ------------------------------------------------------------------ cell math (row vectors a, b, c)
  function cosd(x) { const c = Math.cos(x * RAD); return Math.abs(c) < 1e-12 ? 0 : c; }
  function cellFromParams(a, b, c, al, be, ga) {
    const ca = cosd(al), cb = cosd(be), cg = cosd(ga), sg = Math.sin(ga * RAD);
    const cy = (ca - cb * cg) / sg;
    return [[a, 0, 0], [b * cg, b * sg, 0], [c * cb, c * cy, c * Math.sqrt(Math.max(0, 1 - cb * cb - cy * cy))]];
  }
  function det3(m) {
    return m[0][0] * (m[1][1] * m[2][2] - m[1][2] * m[2][1]) - m[0][1] * (m[1][0] * m[2][2] - m[1][2] * m[2][0]) +
      m[0][2] * (m[1][0] * m[2][1] - m[1][1] * m[2][0]);
  }
  function inv3(M) {
    const [[a, b, c], [d, e, f], [g, h, i]] = M;
    const A = e * i - f * h, B = -(d * i - f * g), C = d * h - e * g, det = a * A + b * B + c * C;
    if (!det) return null;
    return [[A / det, -(b * i - c * h) / det, (b * f - c * e) / det],
      [B / det, (a * i - c * g) / det, -(a * f - c * d) / det],
      [C / det, -(a * h - b * g) / det, (a * e - b * d) / det]];
  }
  const f2c = (C, u, v, w) => [u * C[0][0] + v * C[1][0] + w * C[2][0], u * C[0][1] + v * C[1][1] + w * C[2][1], u * C[0][2] + v * C[1][2] + w * C[2][2]];
  const c2f = (I, x, y, z) => [x * I[0][0] + y * I[1][0] + z * I[2][0], x * I[0][1] + y * I[1][1] + z * I[2][1], x * I[0][2] + y * I[1][2] + z * I[2][2]];
  const scaleCell = (C, s) => C.map(v => v.map(x => x * s));
  const sameCell = (a, b) => a.every((v, i) => v.every((x, k) => Math.abs(x - b[i][k]) < 1e-6));
  const goodCell = C => !!C && C.length === 3 && C.every(v => v.length === 3 && fin(v)) && Math.abs(det3(C)) > 1e-8;
  /** Per-frame cells: one shared cell goes to ds.cell; varying cells stay on the frames as well. */
  function settleCells(ds, useLast) {
    const known = ds.frames.map(f => f.cell).filter(Boolean);
    if (!known.length) return;
    ds.cell = useLast ? known[known.length - 1] : known[0];
    if (known.every(c => sameCell(c, known[0]))) for (const f of ds.frames) delete f.cell;
  }
  /** Cell from DCD/LAMMPS-style parameters given as lengths and angles (degrees). */
  function boxCell(a, b, c, al, be, ga) { return cellFromParams(a, b, c, al, be, ga); }

  // ================================================================== CIF / mmCIF
  function cifTokens(text) {
    const v = [], q = [];
    const L = splitLines(text);
    for (let li = 0; li < L.length; li++) {
      const s = L[li];
      if (s.charCodeAt(0) === 59) {                         // ';' text field
        const buf = [s.slice(1)];
        while (++li < L.length && L[li].charCodeAt(0) !== 59) buf.push(L[li]);
        v.push(buf.join('\n').trim()); q.push(1);
        continue;
      }
      const n = s.length;
      let i = 0;
      while (i < n) {
        const c = s.charCodeAt(i);
        if (c <= 32) { i++; continue; }
        if (c === 35) break;                                // '#' comment
        let j = i + 1;
        if (c === 39 || c === 34) {                         // a quote ends only when followed by a blank
          while (j < n && !(s.charCodeAt(j) === c && (j + 1 === n || s.charCodeAt(j + 1) <= 32))) j++;
          v.push(s.slice(i + 1, j)); q.push(1); i = j + 1; continue;
        }
        while (j < n && s.charCodeAt(j) > 32) j++;
        v.push(s.slice(i, j)); q.push(0); i = j;
      }
    }
    return { v, q };
  }
  function cifBlocks(text) {
    const { v, q } = cifTokens(text);
    const blocks = [];
    let blk = null;
    const kw = k => !q[k] && /^(data_|save_|loop_$|global_$|stop_$)/i.test(v[k]);
    const isTag = k => !q[k] && v[k][0] === '_';
    const need = () => { if (!blk) { blk = { name: '', tags: new Map(), loops: [] }; blocks.push(blk); } return blk; };
    for (let i = 0; i < v.length;) {
      const t = v[i];
      if (!q[i]) {
        const lo = t.toLowerCase();
        if (lo.startsWith('data_')) { blk = { name: t.slice(5), tags: new Map(), loops: [] }; blocks.push(blk); i++; continue; }
        if (lo === 'loop_') {
          i++;
          const tags = [], vals = [];
          while (i < v.length && isTag(i)) tags.push(v[i++].toLowerCase());
          while (i < v.length && !isTag(i) && !kw(i)) vals.push(v[i++]);
          if (tags.length) need().loops.push({ tags, vals });
          continue;
        }
        if (t[0] === '_') {
          if (i + 1 < v.length && !isTag(i + 1) && !kw(i + 1)) { need().tags.set(lo, v[i + 1]); i += 2; } else i++;
          continue;
        }
      }
      i++;
    }
    return blocks;
  }
  const unk = s => s === undefined || s === null || s === '?' || s === '.';
  const cifNum = s => (unk(s) ? NaN : num(String(s).replace(/\(\d*\)\s*$/, '')));
  function cifLoop(blk, names) { for (const lp of blk.loops) for (const n of names) if (lp.tags.indexOf(n) >= 0) return lp; return null; }
  function cifCol(lp, names) {
    if (!lp) return null;
    for (const n of names) {
      const k = lp.tags.indexOf(n);
      if (k < 0) continue;
      const w = lp.tags.length, rows = Math.floor(lp.vals.length / w), out = new Array(rows);
      for (let r = 0; r < rows; r++) out[r] = lp.vals[r * w + k];
      return out;
    }
    return null;
  }
  function cifTag(blk, names) {
    for (const n of names) if (blk.tags.has(n)) return blk.tags.get(n);
    for (const n of names) { const c = cifCol(cifLoop(blk, [n]), [n]); if (c && c.length) return c[0]; }
    return undefined;
  }
  function fracVal(s) { const k = s.indexOf('/'); return k >= 0 ? num(s.slice(0, k)) / num(s.slice(k + 1)) : num(s); }
  /** "-x+1/2, y, z-y" -> {R (rows act on fractional coordinates), t}; null if unreadable. */
  function parseSymop(s) {
    const parts = String(s).toLowerCase().replace(/\s+/g, '').split(',');
    if (parts.length !== 3) return null;
    const R = [[0, 0, 0], [0, 0, 0], [0, 0, 0]], t = [0, 0, 0];
    for (let r = 0; r < 3; r++) {
      let p = parts[r];
      if (!p) return null;
      if (p[0] !== '+' && p[0] !== '-') p = '+' + p;
      const terms = p.match(/[+-][^+-]+/g);
      if (!terms || terms.join('') !== p) return null;
      for (const term of terms) {
        const sg = term[0] === '-' ? -1 : 1, body = term.slice(1), m = /[xyz]/.exec(body);
        if (m) {
          let c = body.replace(m[0], '').replace('*', ''), val = 1;
          if (c) { if (c[0] === '/') c = '1' + c; val = fracVal(c); }
          if (!Number.isFinite(val)) return null;
          R[r]['xyz'.indexOf(m[0])] += sg * val;
        } else {
          const val = fracVal(body);
          if (!Number.isFinite(val)) return null;
          t[r] += sg * val;
        }
      }
    }
    return { R, t };
  }
  const isIdentity = op => op.R.every((row, i) => row.every((x, k) => Math.abs(x - (i === k ? 1 : 0)) < 1e-9)) && op.t.every(x => Math.abs(x - Math.round(x)) < 1e-9);
  /** Applies the operations, wraps into [0,1) and merges positions closer than tol (Å). */
  function expandSites(cell, sites, ops, tol) {
    const nb = cell.map(v => Math.max(1, Math.min(64, Math.floor(Math.hypot(v[0], v[1], v[2]) / 1.5))));
    const grid = new Map(), out = [], t2 = tol * tol;
    const md = (a, n) => ((a % n) + n) % n;
    let mixed = 0;
    for (const s of sites) for (const op of ops) {
      const f = [0, 1, 2].map(r => {
        const x = op.R[r][0] * s.f[0] + op.R[r][1] * s.f[1] + op.R[r][2] * s.f[2] + op.t[r];
        const w = x - Math.floor(x + 1e-6);
        return w < 0 ? 0 : w;
      });
      const b = f.map((x, k) => Math.min(nb[k] - 1, Math.floor(x * nb[k])));
      let dup = null;
      for (let da = -1; da <= 1 && !dup; da++) for (let db = -1; db <= 1 && !dup; db++) for (let dc = -1; dc <= 1 && !dup; dc++) {
        const list = grid.get((md(b[0] + da, nb[0]) * nb[1] + md(b[1] + db, nb[1])) * nb[2] + md(b[2] + dc, nb[2]));
        if (!list) continue;
        for (const o of list) {
          let d0 = f[0] - o.f[0], d1 = f[1] - o.f[1], d2 = f[2] - o.f[2];
          d0 -= Math.round(d0); d1 -= Math.round(d1); d2 -= Math.round(d2);
          const r = f2c(cell, d0, d1, d2);
          if (r[0] * r[0] + r[1] * r[1] + r[2] * r[2] < t2) { dup = o; break; }
        }
      }
      if (dup) { if (dup.z !== s.z) mixed++; continue; }
      const site = { z: s.z, f, label: s.label, occ: s.occ };
      out.push(site);
      const key = (b[0] * nb[1] + b[1]) * nb[2] + b[2];
      if (!grid.has(key)) grid.set(key, []);
      grid.get(key).push(site);
    }
    return { sites: out, mixed };
  }
  const SYMOP_TAGS = ['_symmetry_equiv_pos_as_xyz', '_space_group_symop_operation_xyz', '_space_group_symop.operation_xyz', '_symmetry_equiv.pos_as_xyz'];

  function readCifSmall(blk, ds) {
    const len = n => cifNum(cifTag(blk, ['_cell_length_' + n, '_cell.length_' + n]));
    const ang = n => { const v = cifNum(cifTag(blk, ['_cell_angle_' + n, '_cell.angle_' + n])); return Number.isFinite(v) ? v : 90; };
    const a = len('a'), b = len('b'), c = len('c');
    const lp = cifLoop(blk, ['_atom_site_fract_x', '_atom_site.fract_x', '_atom_site_cartn_x']);
    const C = n => cifCol(lp, ['_atom_site_' + n, '_atom_site.' + n]);
    const label = C('label'), type = C('type_symbol'), occ = C('occupancy'), calc = C('calc_flag');
    let fx = C('fract_x'), fy = C('fract_y'), fz = C('fract_z');
    const cart = !fx;
    if (cart) { fx = C('cartn_x'); fy = C('cartn_y'); fz = C('cartn_z'); }
    if (!fx || !fy || !fz) throw err('The atom_site loop has no coordinates.');
    const cell = [a, b, c].every(v => v > 0) ? cellFromParams(a, b, c, ang('alpha'), ang('beta'), ang('gamma')) : null;
    if (!cell && !cart) throw err('Fractional coordinates need _cell_length_a, _b and _c.');
    const inv = cell ? inv3(cell) : null;
    const sites = [];
    let unknown = 0, partial = 0;
    for (let r = 0; r < fx.length; r++) {
      if (calc && /^dum/i.test(calc[r] || '')) continue;
      let p = [cifNum(fx[r]), cifNum(fy[r]), cifNum(fz[r])];
      if (!fin(p)) continue;
      if (cart && inv) p = c2f(inv, p[0], p[1], p[2]);
      let z = type && !unk(type[r]) ? zOf(type[r]) : -1;
      if (z < 0 && label) z = zOf(label[r]);
      if (z < 0) { unknown++; z = 0; }
      const o = occ ? cifNum(occ[r]) : 1;
      if (Number.isFinite(o) && o < 0.999) partial++;
      sites.push({ z, f: p, label: label ? label[r] : (type ? type[r] : ''), occ: Number.isFinite(o) ? o : 1 });
    }
    if (!sites.length) throw err('No atom sites with valid coordinates.');
    const sg = cifTag(blk, ['_symmetry_space_group_name_h-m', '_space_group_name_h-m_alt', '_symmetry.space_group_name_h-m', '_space_group.name_h-m_alt']);
    if (!unk(sg)) ds.meta.spaceGroup = String(sg).trim();
    const opStr = cifCol(cifLoop(blk, SYMOP_TAGS), SYMOP_TAGS) || [];
    let ops = opStr.map(parseSymop).filter(Boolean);
    if (ops.length < opStr.length) ds.warnings.push((opStr.length - ops.length) + ' symmetry operation(s) could not be read and were skipped.');
    if (!ops.length && !unk(sg) && !/^p\s*1$/i.test(String(sg).trim())) {
      ds.warnings.push('Space group ' + String(sg).trim() + ' is named but no symmetry operations are listed; only the asymmetric unit is shown.');
    }
    let list = sites;
    if (cell && ops.length && !(ops.length === 1 && isIdentity(ops[0]))) {
      const ex = expandSites(cell, sites, ops, 0.01);
      list = ex.sites;
      if (ex.mixed) ds.warnings.push(ex.mixed + ' symmetry-generated position(s) coincide with an atom of another element (mixed site); the first atom was kept.');
    }
    const Z = [], xyz = [], names = [], occs = [];
    for (const s of list) {
      const r = cell ? f2c(cell, s.f[0], s.f[1], s.f[2]) : s.f;
      Z.push(s.z); xyz.push(r[0], r[1], r[2]); names.push(s.label); occs.push(s.occ);
    }
    const f = newFrame(Z, xyz);
    f.atomNames = names;
    if (partial) { f.occupancy = occs; ds.warnings.push(partial + ' site(s) are partially occupied (disorder); all positions are shown.'); }
    if (unknown) ds.warnings.push(unknown + ' site(s) with an unrecognised element are shown as dummy atoms.');
    ds.frames.push(f);
    ds.cell = cell;
  }

  function readMmcif(blk, ds) {
    const lp = cifLoop(blk, ['_atom_site.cartn_x']);
    const C = n => cifCol(lp, ['_atom_site.' + n]);
    const X = C('cartn_x'), Y = C('cartn_y'), W = C('cartn_z');
    if (!X || !Y || !W) throw err('The atom_site loop has no Cartn_x/y/z columns.');
    const type = C('type_symbol'), latom = C('label_atom_id'), aatom = C('auth_atom_id'), lcomp = C('label_comp_id'),
      acomp = C('auth_comp_id'), lasym = C('label_asym_id'), aasym = C('auth_asym_id'), lseq = C('label_seq_id'),
      aseq = C('auth_seq_id'), alt = C('label_alt_id'), model = C('pdbx_pdb_model_num'), ins = C('pdbx_pdb_ins_code');
    const pick = (p, q, r) => (p && !unk(p[r]) ? p[r] : (q && !unk(q[r]) ? q[r] : ''));
    const models = [], byId = new Map();
    let unknown = 0, altSkipped = 0;
    for (let r = 0; r < X.length; r++) {
      const mid = model ? model[r] : '1';
      let M = byId.get(mid);
      if (!M) { M = { Z: [], xyz: [], names: [], res: [], ids: [], ch: [], firstAlt: new Map() }; byId.set(mid, M); models.push(M); }
      const an = pick(aatom, latom, r), rn = pick(acomp, lcomp, r), ch = pick(aasym, lasym, r), rs = pick(aseq, lseq, r);
      if (alt && !unk(alt[r])) {                       // keep the first alternate location of each residue
        const key = ch + '|' + rs + '|' + (ins && !unk(ins[r]) ? ins[r] : '');
        const fa = M.firstAlt.get(key);
        if (fa === undefined) M.firstAlt.set(key, alt[r]);
        else if (fa !== alt[r]) { altSkipped++; continue; }
      }
      const x = num(X[r]), y = num(Y[r]), z = num(W[r]);
      if (!Number.isFinite(x + y + z)) continue;
      let el = type && !unk(type[r]) ? zOf(type[r]) : -1;
      if (el < 0) el = zFromName(an, rn);
      if (el < 0) { unknown++; el = 0; }
      M.Z.push(el); M.xyz.push(x, y, z); M.names.push(an); M.res.push(rn); M.ids.push(parseInt(rs, 10) || 0); M.ch.push(ch);
    }
    if (!models.length || !models[0].Z.length) throw err('The atom_site loop holds no atoms.');
    models.forEach((M, k) => {
      const f = newFrame(M.Z, M.xyz);
      f.atomNames = M.names; f.resNames = M.res; f.resIds = M.ids; f.chains = M.ch;
      if (models.length > 1) f.label = 'Model ' + (k + 1);
      ds.frames.push(f);
    });
    const p = ['length_a', 'length_b', 'length_c', 'angle_alpha', 'angle_beta', 'angle_gamma'].map(n => cifNum(cifTag(blk, ['_cell.' + n])));
    if (fin(p.slice(0, 3)) && !(p[0] === 1 && p[1] === 1 && p[2] === 1)) ds.cell = cellFromParams(p[0], p[1], p[2], p[3] || 90, p[4] || 90, p[5] || 90);
    const t = cifTag(blk, ['_struct.title']);
    if (!unk(t) && String(t).trim()) ds.title = String(t).trim();
    // secondary structure: struct_conf (HELX_*) and struct_sheet_range, in author numbering like the atoms
    const ssl = [], range = (lp, pre, type, keep) => {
      if (!lp) return;
      const col = n => cifCol(lp, [pre + n]), ty = col('conf_type_id');
      const ch = col('beg_auth_asym_id') || col('beg_label_asym_id'), b = col('beg_auth_seq_id') || col('beg_label_seq_id'), e = col('end_auth_seq_id') || col('end_label_seq_id');
      if (!b || !e) return;
      for (let r = 0; r < b.length; r++) {
        if (keep && !(ty && keep.test(String(ty[r])))) continue;
        const a = parseInt(b[r], 10), z = parseInt(e[r], 10);
        if (Number.isFinite(a) && Number.isFinite(z)) ssl.push({ type, chain: ch && !unk(ch[r]) ? String(ch[r]) : '', start: { resId: a }, end: { resId: z } });
      }
    };
    range(cifLoop(blk, ['_struct_conf.conf_type_id']), '_struct_conf.', 'helix', /^HELX/i);
    range(cifLoop(blk, ['_struct_sheet_range.sheet_id', '_struct_sheet_range.beg_auth_seq_id']), '_struct_sheet_range.', 'sheet', null);
    if (ssl.length) ds.secondaryStructure = ssl;
    if (altSkipped) ds.warnings.push('Alternate locations: the first one of each residue is shown (' + altSkipped + ' atom record(s) skipped).');
    if (unknown) ds.warnings.push(unknown + ' atom(s) with an unrecognised element are shown as dummy atoms.');
  }

  function readCif(text, fileName) {
    const blocks = cifBlocks(text);
    if (!blocks.length) throw err('No data_ block found.');
    const isMm = b => !!cifLoop(b, ['_atom_site.cartn_x']);
    const hasAtoms = b => isMm(b) || !!cifLoop(b, ['_atom_site_fract_x', '_atom_site.fract_x', '_atom_site_cartn_x']);
    const k = blocks.findIndex(hasAtoms);
    if (k < 0) throw err('No atom sites (_atom_site_fract_x or _atom_site.Cartn_x) were found.');
    const blk = blocks[k], mm = isMm(blk);
    const ds = newDs(mm ? 'mmcif' : 'cif', fileName);
    if (blk.name) ds.title = blk.name;
    if (mm) readMmcif(blk, ds); else readCifSmall(blk, ds);
    const others = blocks.filter((b, i) => i !== k).map(b => b.name || '(unnamed)');
    if (others.length) {
      ds.warnings.push('The file holds ' + blocks.length + ' data blocks; "' + (blk.name || '(unnamed)') + '" is shown. Other blocks: ' +
        others.slice(0, 20).join(', ') + (others.length > 20 ? ', …' : '') + '.');
    }
    return ds;
  }
  function sniffCif(h) {
    if (!/^\s*data_/m.test(h)) return 0;
    if (/^\s*_(atom_site_fract_x|cell_length_a|symmetry_equiv_pos_as_xyz|space_group_symop_operation_xyz)\b/im.test(h)) return 0.9;
    return /^\s*loop_/m.test(h) ? 0.5 : 0.3;
  }
  function sniffMmcif(h) {
    if (/^\s*_atom_site\.(group_PDB|Cartn_x|label_atom_id)\b/im.test(h)) return 0.95;
    return /^\s*data_/m.test(h) && /^\s*_entry\.id/m.test(h) ? 0.6 : 0;
  }

  // ================================================================== Tripos Mol2
  const MOL2_BOND = { ar: 1.5, am: 1, du: 1, un: 1, 1: 1, 2: 2, 3: 3 };
  function mol2Z(type, name) {
    if (/^(LP|Du)(\.|$)/i.test(type)) return 0;
    const m = /^([A-Z][a-z]?)(?:\.|$)/.exec(type);           // Tripos types: element symbol [.hybridisation]
    if (m && !/^(Any|Hal|Het|Hev)$/.test(type)) { const z = BM.el.toZ(m[1]); if (z > 0) return z; }
    let z = zFromName(name);                                 // GAFF / Amber force-field types: use the atom name
    if (z > 0) return z;
    const g = /^(cl|br|[a-z])/i.exec(type);
    if (g) { z = BM.el.toZ(g[1].length === 2 ? g[1][0].toUpperCase() + g[1][1].toLowerCase() : g[1].toUpperCase()); if (z > 0) return z; }
    return -1;
  }
  function readMol2(text, fileName) {
    const L = splitLines(text);
    const mols = [];
    let cur = null, sec = null;
    for (let i = 0; i < L.length; i++) {
      const t = L[i].trim();
      if (!t || t[0] === '#') continue;
      if (t[0] === '@') {
        const m = /^@<TRIPOS>(\w+)/i.exec(t);
        sec = m ? m[1].toUpperCase() : null;
        if (sec === 'MOLECULE') {
          cur = { name: (L[i + 1] || '').trim(), charge: (L[i + 4] || '').trim().toUpperCase(), atoms: [], bonds: [], subs: new Map(), cell: null };
          mols.push(cur);
          i += 2;
        }
        continue;
      }
      if (!cur) continue;
      const f = tok(t);
      if (sec === 'ATOM' && f.length >= 6) cur.atoms.push(f);
      else if (sec === 'BOND' && f.length >= 4) cur.bonds.push(f);
      else if (sec === 'SUBSTRUCTURE' && f.length >= 2) cur.subs.set(f[0], { name: f[1], chain: f[5] });
      else if (sec === 'CRYSIN' && f.length >= 6) { const p = f.slice(0, 6).map(num); if (fin(p)) cur.cell = p; }
    }
    if (!mols.length || !mols[0].atoms.length) throw err('No @<TRIPOS>MOLECULE record with an ATOM section was found.');
    const ds = newDs('mol2', fileName);
    let unknown = 0;
    const build = mol => {
      const Z = [], xyz = [], names = [], res = [], ids = [], ch = [], q = [], idx = new Map();
      let hasQ = false, hasChain = false;
      mol.atoms.forEach((f, i) => {
        idx.set(f[0], i);
        let z = mol2Z(f[5], f[1]);
        if (z < 0) { unknown++; z = 0; }
        const p = [num(f[2]), num(f[3]), num(f[4])];
        if (!fin(p)) throw err('Bad coordinates for atom ' + f[0] + '.');
        Z.push(z); xyz.push(p[0], p[1], p[2]); names.push(f[1]);
        const sub = f[6] !== undefined ? mol.subs.get(f[6]) : null;
        const sname = f[7] !== undefined ? f[7] : (sub ? sub.name : '');
        let rn = sname || 'UNL', rid = parseInt(f[6], 10) || 1;
        const mm = sname ? /^(.*?[^\d-])(-?\d+)$/.exec(sname) : null;
        if (mm && /[A-Za-z]/.test(mm[1])) { rn = mm[1]; rid = parseInt(mm[2], 10); }
        res.push(rn); ids.push(rid);
        const c = sub && sub.chain && !/^\*+$/.test(sub.chain) ? sub.chain : '';
        if (c) hasChain = true;
        ch.push(c);
        const qq = num(f[8]);
        if (Number.isFinite(qq)) hasQ = true;
        q.push(Number.isFinite(qq) ? qq : 0);
      });
      const bonds = [];
      for (const b of mol.bonds) {
        const i = idx.get(b[1]), j = idx.get(b[2]), ty = String(b[3]).toLowerCase();
        if (i === undefined || j === undefined || ty === 'nc') continue;
        bonds.push([i, j, MOL2_BOND[ty] !== undefined ? MOL2_BOND[ty] : (num(ty) > 0 ? num(ty) : 1)]);
      }
      const fr = newFrame(Z, xyz);
      fr.atomNames = names; fr.resNames = res; fr.resIds = ids; fr.chains = hasChain ? ch : null;
      fr.label = mol.name && !/^\*+$/.test(mol.name) ? mol.name : null;
      return { fr, bonds, q: hasQ && mol.charge !== 'NO_CHARGES' ? q : null };
    };
    const first = build(mols[0]);
    ds.frames.push(first.fr);
    ds.bonds = first.bonds;
    if (first.q) ds.charges = { Mol2: first.q };
    for (let k = 1; k < mols.length; k++) {
      const m = mols[k];
      if (m.atoms.length !== mols[0].atoms.length) { ds.warnings.push('The file holds ' + mols.length + ' molecules; molecules that differ from the first are not shown.'); break; }
      const b = build(m);
      if (b.fr.Z.some((z, i) => z !== first.fr.Z[i])) { ds.warnings.push('The file holds ' + mols.length + ' molecules; molecules that differ from the first are not shown.'); break; }
      ds.frames.push(b.fr);
    }
    if (first.fr.label) ds.title = first.fr.label;
    const c = mols[0].cell;
    if (c && c[0] > 0 && c[1] > 0 && c[2] > 0 && !(c[0] === 1 && c[1] === 1 && c[2] === 1)) ds.cell = cellFromParams(c[0], c[1], c[2], c[3], c[4], c[5]);
    if (unknown) ds.warnings.push(unknown + ' atom(s) with an unrecognised type are shown as dummy atoms.');
    return ds;
  }

  // ================================================================== CML
  const xmlUnesc = s => s.replace(/&(lt|gt|quot|apos|amp|#\d+|#x[0-9a-f]+);/gi, (m, e) => {
    const k = e.toLowerCase();
    if (k[0] !== '#') return { lt: '<', gt: '>', quot: '"', apos: "'", amp: '&' }[k];
    return String.fromCharCode(k[1] === 'x' ? parseInt(k.slice(2), 16) : parseInt(k.slice(1), 10));
  });
  function xmlAttrs(s) {
    const o = {}, re = /([\w:.-]+)\s*=\s*(?:"([^"]*)"|'([^']*)')/g;
    let m;
    while ((m = re.exec(s))) o[m[1].replace(/^[\w.-]*:/, '')] = xmlUnesc(m[2] !== undefined ? m[2] : m[3]);
    return o;
  }
  const CML_ORDER = { 1: 1, 2: 2, 3: 3, s: 1, d: 2, t: 3, a: 1.5 };
  const CML_LEAF = new Set(['string', 'float', 'integer', 'stringArray', 'floatArray', 'integerArray', 'coordinate2', 'coordinate3', 'scalar']);
  function readCML(text, fileName) {
    text = wholeText(text, 'CML');
    const ds = newDs('cml', fileName);
    const re = /<!--[\s\S]*?-->|<\?[\s\S]*?\?>|<!\[CDATA\[[\s\S]*?\]\]>|<![^>]*>|<(\/?)([A-Za-z_][\w:.-]*)((?:[^>"']|"[^"]*"|'[^']*')*)>/g;
    const stack = [], atoms = [], bonds = [], cp = {};
    let molDepth = -1, molDone = false, extra = 0, title = '', m;
    const SPLIT = s => String(s).trim().split(/\s+/);
    const CELLKEY = { a: 'a', acell: 'a', b: 'b', bcell: 'b', c: 'c', ccell: 'c', alpha: 'alpha', beta: 'beta', gamma: 'gamma' };
    function expandAtoms(o) {
      if (o.elementType === undefined) return;
      const el = SPLIT(o.elementType), ids = o.atomID || o.atomId || o.id, idl = ids ? SPLIT(ids) : null;
      const keys = ['x3', 'y3', 'z3', 'x2', 'y2', 'xFract', 'yFract', 'zFract', 'formalCharge'];
      const cols = keys.map(k => (o[k] !== undefined ? SPLIT(o[k]) : null));
      for (let i = 0; i < el.length; i++) {
        const a = { id: idl ? idl[i] : 'a' + (atoms.length + 1), elementType: el[i] };
        keys.forEach((k, j) => { if (cols[j]) a[k] = cols[j][i]; });
        atoms.push(a);
      }
    }
    function onEnd(e, content) {
      const n = e.name, p = e.parent;
      if (n === 'molecule') { if (stack.length === molDepth) molDone = true; return; }
      if (molDone) return;
      if (CML_LEAF.has(n) && p) {
        const key = (e.a.builtin || e.a.dictRef || e.a.title || '').replace(/^.*:/, '');
        const v = content.trim();
        if (p.name === 'crystal') { const ck = CELLKEY[key.toLowerCase()]; if (ck) cp[ck] = num(v); return; }
        if (key === 'atomRef') (p.kids.atomRef = p.kids.atomRef || []).push(v);
        else if (key) p.kids[key] = v;
        return;
      }
      if (n === 'atom') { atoms.push(Object.assign({}, e.kids, e.a)); return; }
      if (n === 'atomArray') { expandAtoms(Object.assign({}, e.kids, e.a)); return; }
      if (n === 'bond') {
        const o = Object.assign({}, e.kids, e.a);
        const refs = o.atomRefs2 || o.atomRefs || (Array.isArray(o.atomRef) ? o.atomRef.join(' ') : null);
        if (refs) { const r = SPLIT(refs); if (r.length >= 2) bonds.push([r[0], r[1], o.order]); }
        return;
      }
      if (n === 'bondArray') {
        const o = Object.assign({}, e.kids, e.a);
        if (o.atomRef1 && o.atomRef2) {
          const a = SPLIT(o.atomRef1), b = SPLIT(o.atomRef2), ord = o.order ? SPLIT(o.order) : [];
          for (let k = 0; k < Math.min(a.length, b.length); k++) bonds.push([a[k], b[k], ord[k]]);
        }
      }
    }
    while ((m = re.exec(text))) {
      if (m[2] === undefined) continue;
      const name = m[2].replace(/^[\w.-]*:/, '');
      if (m[1] === '/') {
        let e = null;
        while (stack.length) { const x = stack.pop(); if (x.name === name) { e = x; break; } }
        if (e) onEnd(e, text.slice(e.cs, m.index));
        continue;
      }
      let raw = m[3];
      const self = /\/\s*$/.test(raw);
      if (self) raw = raw.replace(/\/\s*$/, '');
      const e = { name, a: xmlAttrs(raw), kids: {}, cs: re.lastIndex, parent: stack[stack.length - 1] || null };
      if (name === 'molecule') {
        if (molDepth < 0) { molDepth = stack.length; title = e.a.title || e.a.id || ''; } else if (molDone) extra++;
      }
      if (self) onEnd(e, ''); else stack.push(e);
    }
    if (!atoms.length) throw err('No <atom> or <atomArray> elements were found.');
    const cell = cp.a > 0 && cp.b > 0 && cp.c > 0 ? cellFromParams(cp.a, cp.b, cp.c, cp.alpha || 90, cp.beta || 90, cp.gamma || 90) : null;
    const Z = [], xyz = [], fc = [], idx = new Map();
    let unknown = 0, missing = 0, twoD = 0, hasFc = false;
    atoms.forEach((a, i) => {
      if (a.id !== undefined) idx.set(String(a.id), i);
      let z = zOf(a.elementType);
      if (z < 0) { unknown++; z = 0; }
      let p = null;
      if (a.xyz3) p = SPLIT(a.xyz3).map(num);
      else if (a.x3 !== undefined) p = [num(a.x3), num(a.y3), num(a.z3)];
      else if (a.xFract !== undefined && cell) p = f2c(cell, num(a.xFract), num(a.yFract), num(a.zFract));
      else if (a.xy2) { const q = SPLIT(a.xy2).map(num); p = [q[0], q[1], 0]; twoD++; }
      else if (a.x2 !== undefined) { p = [num(a.x2), num(a.y2), 0]; twoD++; }
      if (!p || p.length < 3 || !fin(p.slice(0, 3))) { missing++; p = [0, 0, 0]; }
      Z.push(z); xyz.push(p[0], p[1], p[2]);
      const q = a.formalCharge !== undefined ? parseInt(a.formalCharge, 10) || 0 : 0;
      if (q) hasFc = true;
      fc.push(q);
    });
    const out = [];
    for (const [r1, r2, o] of bonds) {
      const i = idx.get(r1), j = idx.get(r2);
      if (i === undefined || j === undefined) continue;
      const key = String(o === undefined ? '1' : o).trim().toLowerCase();
      let ord = CML_ORDER[key];
      if (ord === undefined) { ord = num(key); if (!(ord > 0)) ord = 1; }
      out.push([i, j, ord]);
    }
    const f = newFrame(Z, xyz);
    if (hasFc) f.formalCharges = fc;
    ds.frames.push(f);
    ds.bonds = out;
    ds.cell = cell;
    if (title) ds.title = title;
    if (extra) ds.warnings.push('The file holds ' + (extra + 1) + ' molecules; only the first is shown.');
    if (twoD) ds.warnings.push(twoD + ' atom(s) have only 2D coordinates (z = 0).');
    if (missing) ds.warnings.push(missing + ' atom(s) have no coordinates and were placed at the origin.');
    if (unknown) ds.warnings.push(unknown + ' atom(s) with an unrecognised element are shown as dummy atoms.');
    return ds;
  }

  // ================================================================== Avogadro CJSON
  function readCJSON(text, fileName) {
    text = wholeText(text, 'CJSON');
    let J;
    try { J = JSON.parse(text); } catch (e) { throw err('Invalid JSON: ' + e.message); }
    const A = J.atoms || {};
    const E = A.elements || {};
    let Z = E.number || E.atomicNumbers;
    if (!Array.isArray(Z) && Array.isArray(E.symbol)) Z = E.symbol.map(s => Math.max(0, zOf(s)));
    if (!Array.isArray(Z) || !Z.length) throw err('No atoms.elements.number array.');
    Z = Z.map(z => (Number.isInteger(+z) && +z >= 0 ? +z : 0));
    const n = Z.length;
    const ds = newDs('cjson', fileName, 'Avogadro');
    const uc = J.unitCell || J.unitcell;
    let cell = null;
    if (uc) {
      if (Array.isArray(uc.cellVectors) && uc.cellVectors.length >= 9) cell = [0, 1, 2].map(k => uc.cellVectors.slice(3 * k, 3 * k + 3).map(Number));
      else if (uc.a > 0 && uc.b > 0 && uc.c > 0) cell = cellFromParams(+uc.a, +uc.b, +uc.c, +(uc.alpha || 90), +(uc.beta || 90), +(uc.gamma || 90));
      if (cell && !goodCell(cell)) cell = null;
    }
    const co = A.coords || {};
    const toCart = fr => { const out = new Float64Array(3 * n); for (let i = 0; i < n; i++) out.set(f2c(cell, +fr[3 * i], +fr[3 * i + 1], +fr[3 * i + 2]), 3 * i); return out; };
    const sets = [];
    if (Array.isArray(co['3dSets']) && co['3dSets'].length) for (const s of co['3dSets']) sets.push(Float64Array.from(s, Number));
    else if (Array.isArray(co['3d'])) sets.push(Float64Array.from(co['3d'], Number));
    else if (Array.isArray(co['3dFractional']) && cell) sets.push(toCart(co['3dFractional']));
    else if (Array.isArray(co['2d'])) { const p = new Float64Array(3 * n); for (let i = 0; i < n; i++) { p[3 * i] = +co['2d'][2 * i]; p[3 * i + 1] = +co['2d'][2 * i + 1]; } sets.push(p); ds.warnings.push('Only 2D coordinates are present (z = 0).'); }
    if (!sets.length) throw err('No atoms.coords.3d (or 3dFractional with a unitCell) array.');
    const labels = Array.isArray(A.labels) && A.labels.length === n ? A.labels.map(String) : null;
    const fcs = Array.isArray(A.formalCharges) && A.formalCharges.length === n ? A.formalCharges.map(q => +q || 0) : null;
    let resN = null, resI = null, resC = null, names = labels;
    if (Array.isArray(J.residues)) {
      resN = new Array(n).fill(''); resI = new Array(n).fill(0); resC = new Array(n).fill('');
      names = names || new Array(n).fill('');
      for (const r of J.residues) {
        if (!r || !r.atoms) continue;
        for (const [an, ia] of Object.entries(r.atoms)) if (ia >= 0 && ia < n) { names[ia] = an; resN[ia] = r.name || ''; resI[ia] = +r.id || 0; resC[ia] = r.chainId || ''; }
      }
    }
    sets.forEach((p, k) => {
      if (p.length < 3 * n || !fin(Array.from(p.subarray(0, 3 * n)))) throw err('Coordinate set ' + (k + 1) + ' does not match the atom count.');
      const f = newFrame(Z, p.subarray(0, 3 * n));
      f.atomNames = names; f.resNames = resN; f.resIds = resI; f.chains = resC;
      if (fcs) f.formalCharges = fcs;
      ds.frames.push(f);
    });
    const B = J.bonds || {};
    const ix = B.connections && (B.connections.index || B.connections.indices);
    if (Array.isArray(ix)) {
      ds.bonds = [];
      for (let k = 0; k + 1 < ix.length; k += 2) {
        const i = +ix[k], j = +ix[k + 1];
        if (i >= 0 && j >= 0 && i < n && j < n) ds.bonds.push([i, j, B.order && +B.order[k / 2] > 0 ? +B.order[k / 2] : 1]);
      }
    }
    const pc = J.partialCharges || A.partialCharges;
    if (pc && typeof pc === 'object') {
      for (const [k, v] of Object.entries(pc)) if (Array.isArray(v) && v.length === n) { ds.charges = ds.charges || {}; ds.charges[k] = v.map(Number); }
    }
    const V = J.vibrations;
    if (V && Array.isArray(V.frequencies) && V.frequencies.length) {
      const nm = V.frequencies.length, ev = V.eigenVectors || V.eigenvectors || [];
      const modes = [];
      for (let k = 0; k < nm; k++) {
        const src = Array.isArray(ev[k]) ? ev[k] : (ev.length === nm * 3 * n ? ev.slice(k * 3 * n, (k + 1) * 3 * n) : []);
        const v = new Float64Array(3 * n);
        for (let q = 0; q < Math.min(3 * n, src.length); q++) v[q] = +src[q] || 0;
        modes.push(v);
      }
      const arr = x => (Array.isArray(x) && x.length === nm ? x.map(Number) : null);
      ds.vibrations = { frameIndex: 0, freqs: V.frequencies.map(Number), ir: arr(V.intensities), raman: arr(V.ramanIntensities),
        depol: null, redMass: null, forceConst: null, sym: null, modes };
    }
    if (cell) ds.cell = cell;
    if (J.name) ds.title = String(J.name);
    return ds;
  }

  // ================================================================== VASP
  function vaspHeader(L, k) {
    const comment = L[k] || '';
    const sc = tok(L[k + 1]).map(num);
    if (!sc.length || !Number.isFinite(sc[0])) throw err('Line ' + (k + 2) + ' must hold the scaling factor.');
    let cell = [0, 1, 2].map(i => tok(L[k + 2 + i]).slice(0, 3).map(num));
    if (!cell.every(v => v.length === 3 && fin(v))) throw err('Lines ' + (k + 3) + '-' + (k + 5) + ' must hold the three lattice vectors.');
    let scale = 1, scale3 = null;
    if (sc.length >= 3 && fin(sc.slice(0, 3))) scale3 = sc.slice(0, 3);                 // per Cartesian axis
    else if (sc[0] < 0) scale = Math.cbrt(-sc[0] / Math.abs(det3(cell)));              // negative = target volume
    else scale = sc[0];
    cell = scale3 ? cell.map(v => v.map((x, i) => x * scale3[i])) : scaleCell(cell, scale);
    let j = k + 5, t = tok(L[j]), names = null;
    if (t.length && !isInt(t[0])) { names = t; j++; t = tok(L[j]); }
    const counts = [];
    for (const x of t) { if (!isInt(x)) break; counts.push(parseInt(x, 10)); }
    if (!counts.length) throw err('Line ' + (j + 1) + ' must hold the number of atoms of each species.');
    return { comment, cell, scale, scale3, names, counts, next: j + 1 };
  }
  function vaspSpecies(h, ds) {
    let zs = null;
    if (h.names) zs = h.names.slice(0, h.counts.length).map(zOf);
    else {
      const t = tok(h.comment).map(zOf);
      if (t.length >= h.counts.length && t.slice(0, h.counts.length).every(z => z > 0)) zs = t.slice(0, h.counts.length);
    }
    if (!zs || zs.length < h.counts.length) {
      ds.warnings.push('The species are not named (VASP 4 format); atoms are shown as dummies. Add the element line above the atom counts.');
      zs = h.counts.map(() => 0);
    }
    if (zs.some(z => z < 0)) { ds.warnings.push('Unknown species name(s) are shown as dummy atoms.'); zs = zs.map(z => Math.max(0, z)); }
    const Z = [], names = [];
    h.counts.forEach((c, s) => { for (let i = 0; i < c; i++) { Z.push(zs[s]); names.push(h.names ? h.names[s] : 'T' + (s + 1)); } });
    return { Z, names };
  }
  function readPoscar(text, fileName) {
    const L = splitLines(text);
    const ds = newDs('poscar', fileName, 'VASP');
    const h = vaspHeader(L, 0);
    const sp = vaspSpecies(h, ds), N = sp.Z.length;
    let j = h.next, sel = false;
    if (/^\s*s/i.test(L[j] || '')) { sel = true; j++; }
    const mode = (L[j] || '').trim();
    if (!mode) throw err('Line ' + (j + 1) + ' must say Direct or Cartesian.');
    const cart = /^[ck]/i.test(mode);
    j++;
    const xyz = new Float64Array(3 * N), flags = sel ? [] : null;
    for (let i = 0; i < N; i++) {
      const t = tok(L[j + i]);
      const x = num(t[0]), y = num(t[1]), z = num(t[2]);
      if (!Number.isFinite(x + y + z)) throw err('Expected ' + N + ' positions; line ' + (j + i + 1) + ' is not a position.');
      let r;
      if (cart) r = h.scale3 ? [x * h.scale3[0], y * h.scale3[1], z * h.scale3[2]] : [x * h.scale, y * h.scale, z * h.scale];
      else r = f2c(h.cell, x, y, z);
      xyz.set(r, 3 * i);
      if (flags) flags.push((t[3] || 'T') + (t[4] || 'T') + (t[5] || 'T'));
    }
    const f = newFrame(sp.Z, xyz);
    f.atomNames = sp.names;
    ds.frames.push(f);
    ds.cell = h.cell;
    ds.title = h.comment.trim() || ds.title;
    if (flags) ds.meta.selectiveDynamics = flags;
    return ds;
  }
  function sniffPoscar(h, name) {
    name = name || '';
    if (/xdatcar|outcar|oszicar|incar|kpoints|potcar|doscar|chgcar|locpot|eigenval|procar|vasprun|elfcar|parchg/.test(name)) return 0;
    let content = 0;
    if (!/configuration=|Konfig=/.test(h)) {
      const L = h.split(/\r\n|\r|\n/, 10);
      const sc = tok(L[1]);
      let ok = L.length >= 8 && (sc.length === 1 || sc.length === 3) && sc.every(isNum);
      for (let k = 2; ok && k <= 4; k++) { const t = tok(L[k]); ok = t.length >= 3 && t.slice(0, 3).every(isNum); }
      if (ok) {
        let k = 5, t = tok(L[k]), named = true;
        // unknown species labels (e.g. "X" for dummy atoms) are read as dummies by readPoscar: lower confidence only
        if (t.length && !isInt(t[0])) { named = t.every(x => BM.el.toZ(x) > 0); ok = t.every(x => /^[A-Za-z]/.test(x)); k++; t = tok(L[k]); }
        ok = ok && t.length > 0 && t.every(isInt);
        let md = (L[k + 1] || '').trim();
        if (/^s/i.test(md)) md = (L[k + 2] || '').trim();
        if (ok && /^[dck]/i.test(md)) content = named ? 0.7 : 0.4;
      }
    }
    if (/^(poscar|contcar)|[._-](poscar|contcar)$/.test(name)) return content ? 0.95 : 0.5;
    return content;
  }
  function readXdatcar(text, fileName) {
    const L = splitLines(text);
    const ds = newDs('xdatcar', fileName, 'VASP');
    let h = null, sp = null, N = 0, k = 0;
    outer: while (k < L.length) {
      const s = L[k];
      if (!s.trim()) { k++; continue; }
      const m = /^\s*(Direct|Cartesian)\s*configuration\s*=?\s*(\d*)/i.exec(s);
      if (m) {
        if (!h) throw err('A configuration appears before the header.');
        const cart = /^c/i.test(m[1]), xyz = new Float64Array(3 * N);
        for (let i = 0; i < N; i++) {
          const t = tok(L[k + 1 + i]);
          const x = num(t[0]), y = num(t[1]), z = num(t[2]);
          if (!Number.isFinite(x + y + z)) { ds.warnings.push('Configuration ' + (ds.frames.length + 1) + ' is incomplete and was skipped.'); break outer; }
          xyz.set(cart ? [x * h.scale, y * h.scale, z * h.scale] : f2c(h.cell, x, y, z), 3 * i);
        }
        const f = newFrame(sp.Z, xyz);
        f.label = 'Configuration ' + (m[2] || ds.frames.length + 1);
        f.cell = h.cell;
        ds.frames.push(f);
        k += N + 1;
        continue;
      }
      h = vaspHeader(L, k);                     // first header, or a repeated one (variable cell)
      if (!sp) { sp = vaspSpecies(h, ds); N = sp.Z.length; ds.title = h.comment.trim() || ds.title; }
      else if (sum(h.counts) !== N) throw err('The atom count changes between configurations.');
      k = h.next;
    }
    if (!ds.frames.length) throw err('No "Direct configuration=" blocks were found.');
    for (const f of ds.frames) f.atomNames = sp.names;
    settleCells(ds, false);
    return ds;
  }
  function readOutcar(text, fileName) {
    const L = splitLines(text);
    const ds = newDs('outcar', fileName, 'VASP');
    const titel = [], vrh = [], pot = [];
    let counts = null, nions = 0, cell = null, Z = null, finalPending = false, lastToten = null;
    const m0 = /vasp\.(\S+)/.exec(L[0] || '');
    if (m0) ds.meta.version = m0[1];
    const species = () => {
      const names = titel.length ? titel : vrh.length ? vrh : pot.slice(0, counts ? counts.length : pot.length);
      if (!counts) counts = nions ? [nions] : [];
      const out = [];
      counts.forEach((c, s) => { const z = Math.max(0, zOf(names[s])); for (let i = 0; i < c; i++) out.push(z); });
      if (!names.length) ds.warnings.push('No POTCAR species names were found; atoms are shown as dummies.');
      return out;
    };
    for (let i = 0; i < L.length; i++) {
      const s = L[i];
      if (s.indexOf('TITEL') >= 0) { const t = tok(s.split('=')[1]); if (t.length) titel.push(t[1] || t[0]); continue; }
      if (s.indexOf('VRHFIN') >= 0) { const m = /VRHFIN\s*=\s*([A-Za-z]+)/.exec(s); if (m) vrh.push(m[1]); continue; }
      if (/^\s*POTCAR:/.test(s)) { const t = tok(s.split(':')[1]); if (t.length) pot.push(t[1] || t[0]); continue; }
      if (s.indexOf('ions per type') >= 0) { counts = tok(s.split('=')[1]).map(x => parseInt(x, 10)).filter(x => x >= 0); continue; }
      if (s.indexOf('NIONS') >= 0) { const m = /NIONS\s*=\s*(\d+)/.exec(s); if (m) nions = parseInt(m[1], 10); continue; }
      if (s.indexOf('direct lattice vectors') >= 0) {
        const c = [1, 2, 3].map(k => tok(L[i + k]).slice(0, 3).map(num));
        if (c.every(v => v.length === 3 && fin(v))) cell = c;
        i += 3; continue;
      }
      if (/^\s*POSITION\s+TOTAL-FORCE/.test(s)) {
        if (!Z) { Z = species(); if (!nions) nions = Z.length; if (Z.length !== nions) throw err('"ions per type" does not add up to NIONS.'); }
        const xyz = new Float64Array(3 * nions);
        let fmax = 0, fsq = 0, j = i + 2;
        for (let a = 0; a < nions; a++, j++) {
          const v = tok(L[j]).map(num);
          if (v.length < 6 || !fin(v.slice(0, 6))) { ds.warnings.push('The last ionic step is incomplete and was skipped.'); i = L.length; break; }
          xyz[3 * a] = v[0]; xyz[3 * a + 1] = v[1]; xyz[3 * a + 2] = v[2];
          fmax = Math.max(fmax, Math.hypot(v[3], v[4], v[5])); fsq += v[3] * v[3] + v[4] * v[4] + v[5] * v[5];
        }
        if (i >= L.length) break;
        const f = newFrame(Z, xyz);
        f.label = 'Ionic step ' + (ds.frames.length + 1);
        f.opt = { maxForce: fmax * EVA_AU, rmsForce: Math.sqrt(fsq / (3 * nions)) * EVA_AU };
        if (cell) f.cell = cell;
        ds.frames.push(f);
        finalPending = false; lastToten = null;
        i = j; continue;
      }
      if (s.indexOf('FREE ENERGIE OF THE ION-ELECTRON SYSTEM') >= 0) { finalPending = true; continue; }
      const te = /free\s+energy\s+TOTEN\s*=\s*(\S+)/.exec(s);
      if (te) {
        const e = num(te[1]) / HARTREE_EV;
        const f = ds.frames[ds.frames.length - 1];
        if (finalPending && f) { f.energy = e; finalPending = false; } else lastToten = e;
      }
    }
    if (!ds.frames.length) throw err('No "POSITION ... TOTAL-FORCE" blocks were found.');
    const lastF = ds.frames[ds.frames.length - 1];
    if (lastF.energy === null && lastToten !== null) lastF.energy = lastToten;
    settleCells(ds, true);
    return ds;
  }

  // ================================================================== Quantum ESPRESSO
  function stripBang(s) {
    let q = null;
    for (let i = 0; i < s.length; i++) {
      const c = s[i];
      if (q) { if (c === q) q = null; } else if (c === "'" || c === '"') q = c; else if (c === '!' || c === '#') return s.slice(0, i);
    }
    return s;
  }
  function parseNamelists(text) {
    const src = splitLines(text).map(stripBang).join('\n');
    const nl = {}, re = /&(\w+)/g;
    let m, end = 0;
    while ((m = re.exec(src))) {
      const ls = src.lastIndexOf('\n', m.index - 1) + 1;
      if (src.slice(ls, m.index).trim()) continue;
      let i = re.lastIndex, q = null;
      for (; i < src.length; i++) { const c = src[i]; if (q) { if (c === q) q = null; } else if (c === "'" || c === '"') q = c; else if (c === '/') break; }
      const body = src.slice(re.lastIndex, i), vals = {};
      const ar = /([A-Za-z_][\w%]*(?:\s*\([^)]*\))?)\s*=\s*('[^']*'|"[^"]*"|[^,\s]+)/g;
      let a;
      while ((a = ar.exec(body))) {
        let v = a[2];
        if (/^['"]/.test(v)) v = v.slice(1, -1);
        else if (/^\.?(t|true)\.?$/i.test(v)) v = true;
        else if (/^\.?(f|false)\.?$/i.test(v)) v = false;
        else if (Number.isFinite(num(v))) v = num(v);
        vals[a[1].replace(/\s+/g, '').toLowerCase()] = v;
      }
      nl[m[1].toLowerCase()] = vals;
      end = i + 1;
      re.lastIndex = i + 1;
    }
    return { nl, rest: src.slice(end) };
  }
  const QE_CARD = /^\s*(ATOMIC_SPECIES|ATOMIC_POSITIONS|CELL_PARAMETERS|K_POINTS|CONSTRAINTS|OCCUPATIONS|ATOMIC_VELOCITIES|ATOMIC_FORCES|ADDITIONAL_K_POINTS|SOLVENTS|HUBBARD)\b\s*[({]?\s*([A-Za-z_]*)/i;
  /** Bravais lattice vectors of pw.x (ibrav), in units of alat; c[2..6] = celldm(2..6). */
  function qeLattice(ibrav, c) {
    const s3 = Math.sqrt(3), b = c[2], cc = c[3];
    switch (ibrav) {
      case 1: return [[1, 0, 0], [0, 1, 0], [0, 0, 1]];
      case 2: return [[-0.5, 0, 0.5], [0, 0.5, 0.5], [-0.5, 0.5, 0]];
      case 3: return [[0.5, 0.5, 0.5], [-0.5, 0.5, 0.5], [-0.5, -0.5, 0.5]];
      case -3: return [[-0.5, 0.5, 0.5], [0.5, -0.5, 0.5], [0.5, 0.5, -0.5]];
      case 4: return [[1, 0, 0], [-0.5, s3 / 2, 0], [0, 0, cc]];
      case 5: { const g = c[4], tx = Math.sqrt((1 - g) / 2), ty = Math.sqrt((1 - g) / 6), tz = Math.sqrt((1 + 2 * g) / 3); return [[tx, -ty, tz], [0, 2 * ty, tz], [-tx, -ty, tz]]; }
      case 6: return [[1, 0, 0], [0, 1, 0], [0, 0, cc]];
      case 7: return [[0.5, -0.5, cc / 2], [0.5, 0.5, cc / 2], [-0.5, -0.5, cc / 2]];
      case 8: return [[1, 0, 0], [0, b, 0], [0, 0, cc]];
      case 9: return [[0.5, b / 2, 0], [-0.5, b / 2, 0], [0, 0, cc]];
      case -9: return [[0.5, -b / 2, 0], [0.5, b / 2, 0], [0, 0, cc]];
      case 10: return [[0.5, 0, cc / 2], [0.5, b / 2, 0], [0, b / 2, cc / 2]];
      case 11: return [[0.5, b / 2, cc / 2], [-0.5, b / 2, cc / 2], [-0.5, -b / 2, cc / 2]];
      case 12: { const g = c[4]; return [[1, 0, 0], [b * g, b * Math.sqrt(1 - g * g), 0], [0, 0, cc]]; }
      case -12: { const g = c[5]; return [[1, 0, 0], [0, b, 0], [cc * g, 0, cc * Math.sqrt(1 - g * g)]]; }
      case 14: {
        const ca = c[4], cb = c[5], cg = c[6], sg = Math.sqrt(1 - cg * cg);
        return [[1, 0, 0], [b * cg, b * sg, 0], [cc * cb, cc * (ca - cb * cg) / sg, cc * Math.sqrt(1 + 2 * ca * cb * cg - ca * ca - cb * cb - cg * cg) / sg]];
      }
      default: return null;
    }
  }
  function zQe(label, sp) {
    let z = zOf(label);
    if (z > 0) return z;
    const s = sp && sp[label];
    if (s) {
      const m = /^([A-Za-z]{1,2})(?=[._\-\d]|$)/.exec(s.pseudo || '');
      if (m) { z = BM.el.toZ(m[1]); if (z > 0) return z; }
      if (s.mass > 0) return zFromMass(s.mass);
    }
    return 0;
  }
  function readQeIn(text, fileName) {
    const { nl, rest } = parseNamelists(text);
    const sys = nl.system;
    if (!sys) throw err('No &SYSTEM namelist found.');
    const ds = newDs('qe-in', fileName, 'Quantum ESPRESSO');
    const cards = {};
    let cur = null;
    for (const raw of splitLines(rest)) {
      const s = raw.trim();
      const m = QE_CARD.exec(s);
      if (m) { cur = cards[m[1].toUpperCase()] = { unit: m[2].toLowerCase(), head: s, lines: [] }; continue; }
      if (cur && s) cur.lines.push(s);
    }
    const ibrav = Math.round(num(sys.ibrav === undefined ? 0 : sys.ibrav));
    const cd = [0, 1, 2, 3, 4, 5, 6].map(k => num(sys['celldm(' + k + ')']));
    let alat = Number.isFinite(cd[1]) ? cd[1] * BOHR : (num(sys.a) > 0 ? num(sys.a) : NaN);   // Å
    let cell;
    if (ibrav === 0) {
      const cp = cards.CELL_PARAMETERS;
      if (!cp) throw err('ibrav = 0 needs a CELL_PARAMETERS card.');
      const v = cp.lines.slice(0, 3).map(s => tok(s).slice(0, 3).map(num));
      if (v.length < 3 || !v.every(r => r.length === 3 && fin(r))) throw err('CELL_PARAMETERS needs three lines of three numbers.');
      const am = /alat\s*=\s*([-+\d.eEdD]+)/i.exec(cp.head);
      if (am) alat = num(am[1]) * BOHR;
      const u = cp.unit;
      if (u === 'alat' || (!u && Number.isFinite(alat))) {
        if (!Number.isFinite(alat)) throw err('CELL_PARAMETERS in alat units needs celldm(1) or A.');
        cell = scaleCell(v, alat);
      } else if (u === 'bohr' || !u) cell = scaleCell(v, BOHR);
      else if (u === 'angstrom') cell = v;
      else throw err('Unknown CELL_PARAMETERS unit "' + u + '".');
      if (!Number.isFinite(alat)) alat = Math.hypot(cell[0][0], cell[0][1], cell[0][2]);
    } else {
      const c = cd.slice();
      if (!Number.isFinite(cd[1]) && num(sys.a) > 0) {
        const A = num(sys.a);
        c[2] = num(sys.b) / A; c[3] = num(sys.c) / A;
        if (ibrav === 14) { c[4] = num(sys.cosbc); c[5] = num(sys.cosac); c[6] = num(sys.cosab); }
        else if (ibrav === -12 || ibrav === -13) c[5] = num(sys.cosac);
        else c[4] = num(sys.cosab);
      }
      if (!Number.isFinite(alat)) throw err('ibrav = ' + ibrav + ' needs celldm(1) or A.');
      const lat = qeLattice(ibrav, c);
      if (!lat) throw err('ibrav = ' + ibrav + ' is not supported; use ibrav = 0 with CELL_PARAMETERS.');
      cell = scaleCell(lat, alat);
      if (!goodCell(cell)) throw err('The celldm / A, B, C values do not define a valid cell for ibrav = ' + ibrav + '.');
    }
    const sp = {};
    for (const s of (cards.ATOMIC_SPECIES || { lines: [] }).lines) { const t = tok(s); if (t.length >= 2) sp[t[0]] = { mass: num(t[1]), pseudo: t[2] || '' }; }
    const ap = cards.ATOMIC_POSITIONS;
    if (!ap) throw err('No ATOMIC_POSITIONS card found.');
    const nat = Math.round(num(sys.nat)) || ap.lines.length;
    if (ap.lines.length < nat) throw err('ATOMIC_POSITIONS lists ' + ap.lines.length + ' atoms, nat = ' + nat + '.');
    const unit = ap.unit || 'alat';
    if (unit === 'crystal_sg') ds.warnings.push('crystal_sg positions are shown without the space-group expansion.');
    const Z = [], xyz = [], names = [];
    for (let k = 0; k < nat; k++) {
      const t = tok(ap.lines[k]);
      const x = num(t[1]), y = num(t[2]), z = num(t[3]);
      if (!Number.isFinite(x + y + z)) throw err('Bad position line: "' + ap.lines[k] + '".');
      let r;
      if (unit.startsWith('crystal')) r = f2c(cell, x, y, z);
      else if (unit === 'alat') r = [x * alat, y * alat, z * alat];
      else if (unit === 'bohr') r = [x * BOHR, y * BOHR, z * BOHR];
      else if (unit === 'angstrom') r = [x, y, z];
      else throw err('Unknown ATOMIC_POSITIONS unit "' + unit + '".');
      Z.push(zQe(t[0], sp)); xyz.push(r[0], r[1], r[2]); names.push(t[0]);
    }
    const f = newFrame(Z, xyz);
    f.atomNames = names;
    ds.frames.push(f);
    ds.cell = cell;
    const ctl = nl.control || {};
    ds.title = ctl.title || ctl.prefix || ds.title;
    if (Number.isFinite(num(sys.tot_charge))) ds.meta.charge = num(sys.tot_charge);
    if (ctl.calculation) ds.meta.jobTypes.push(String(ctl.calculation));
    return ds;
  }
  function readQeOut(text, fileName) {
    const L = splitLines(text);
    const ds = newDs('qe-out', fileName, 'Quantum ESPRESSO');
    let alat = NaN, nat = 0, cell = null;
    const sp = {};
    const zLab = lab => { const z = zOf(lab); return z > 0 ? z : (sp[lab] || 0); };
    const push = (Z, xyz, label) => {
      const prev = ds.frames[ds.frames.length - 1];
      const pc = prev ? prev.cell || null : null;
      if (prev && prev.xyz.length === xyz.length && prev.xyz.every((v, i) => Math.abs(v - xyz[i]) < 1e-5) && ((!pc && !cell) || (pc && cell && sameCell(pc, cell)))) {
        prev.cell = cell; if (label) prev.label = label; return prev;
      }
      const f = newFrame(Z, xyz);
      f.cell = cell; f.label = label || 'Step ' + (ds.frames.length + 1);
      ds.frames.push(f);
      return f;
    };
    for (let i = 0; i < L.length; i++) {
      const s = L[i];
      let m;
      if ((m = /lattice parameter \(alat\)\s*=\s*([-\d.eE+]+)/.exec(s))) { alat = num(m[1]) * BOHR; continue; }
      if ((m = /number of atoms\/cell\s*=\s*(\d+)/.exec(s))) { nat = parseInt(m[1], 10); continue; }
      if (/crystal axes: \(cart\. coord\. in units of alat\)/.test(s)) {
        const c = [1, 2, 3].map(k => { const q = /=\s*\(\s*([^)]*)\)/.exec(L[i + k] || ''); return q ? tok(q[1]).map(num) : []; });
        if (c.every(v => v.length === 3 && fin(v)) && Number.isFinite(alat)) cell = scaleCell(c, alat);
        i += 3; continue;
      }
      if (/atomic species\s+valence\s+mass\s+pseudopotential/.test(s)) {
        for (let j = i + 1; j < L.length && L[j].trim(); j++) {
          const t = tok(L[j]), pm = /\s([A-Za-z]{1,2})\s*\(/.exec(L[j]);
          if (t.length >= 3) sp[t[0]] = pm && BM.el.toZ(pm[1]) > 0 ? BM.el.toZ(pm[1]) : Math.max(0, zOf(t[0]));
        }
        continue;
      }
      if (/site n\.\s+atom\s+positions \(alat units\)/.test(s)) {
        const Z = [], xyz = [];
        for (let j = i + 1; j <= i + nat; j++) {
          const q = /^\s*\d+\s+(\S+)\s+tau\(\s*\d+\)\s*=\s*\(\s*([^)]*)\)/.exec(L[j] || '');
          if (!q) break;
          const v = tok(q[2]).map(num);
          Z.push(zLab(q[1])); xyz.push(v[0] * alat, v[1] * alat, v[2] * alat);
        }
        if (Z.length === nat && nat && fin(xyz)) push(Z, xyz, null);
        i += nat; continue;
      }
      if ((m = /^\s*CELL_PARAMETERS\s*\(?\s*(\w+)\s*=?\s*([-\d.eE+]*)/.exec(s))) {
        const v = [1, 2, 3].map(k => tok(L[i + k]).slice(0, 3).map(num));
        if (v.every(r => r.length === 3 && fin(r))) {
          const u = m[1].toLowerCase();
          cell = u === 'alat' ? scaleCell(v, (m[2] ? num(m[2]) * BOHR : alat)) : u === 'bohr' ? scaleCell(v, BOHR) : v;
        }
        i += 3; continue;
      }
      if ((m = /^\s*ATOMIC_POSITIONS\s*[({]?\s*(\w*)/.exec(s))) {
        const u = (m[1] || 'alat').toLowerCase(), Z = [], xyz = [];
        if (u.startsWith('crystal') && !cell) throw err('Crystal coordinates appear before the cell.');
        for (let j = i + 1; j <= i + nat; j++) {
          const t = tok(L[j]);
          const x = num(t[1]), y = num(t[2]), z = num(t[3]);
          if (!Number.isFinite(x + y + z)) break;
          const r = u.startsWith('crystal') ? f2c(cell, x, y, z) : u === 'bohr' ? [x * BOHR, y * BOHR, z * BOHR] : u === 'angstrom' ? [x, y, z] : [x * alat, y * alat, z * alat];
          Z.push(zLab(t[0])); xyz.push(r[0], r[1], r[2]);
        }
        if (Z.length === nat && nat) push(Z, xyz, null);
        i += nat; continue;
      }
      if ((m = /^\s*!\s+total energy\s*=\s*([-\d.eEdD+]+)\s*Ry/.exec(s))) {
        const f = ds.frames[ds.frames.length - 1];
        if (f) f.energy = num(m[1]) / 2;
        continue;
      }
      if (/Forces acting on atoms/.test(s)) {
        const f = ds.frames[ds.frames.length - 1];
        let got = 0, fmax = 0, fsq = 0;
        for (let j = i + 1; j < L.length && got < nat && j < i + nat + 4; j++) {
          const q = /atom\s+\d+\s+type\s+\d+\s+force\s*=\s*(\S+)\s+(\S+)\s+(\S+)/.exec(L[j]);
          if (!q) continue;
          const v = [num(q[1]), num(q[2]), num(q[3])];
          fmax = Math.max(fmax, Math.hypot(v[0], v[1], v[2])); fsq += v[0] * v[0] + v[1] * v[1] + v[2] * v[2]; got++;
        }
        if (f && got === nat && nat) f.opt = { maxForce: fmax / 2, rmsForce: Math.sqrt(fsq / (3 * nat)) / 2 };   // Ry/Bohr -> Ha/Bohr
        continue;
      }
      if (/End final coordinates/.test(s)) { const f = ds.frames[ds.frames.length - 1]; if (f) { f.converged = true; f.label = 'Final'; } continue; }
      if (/bfgs converged/i.test(s)) { const f = ds.frames[ds.frames.length - 1]; if (f) f.converged = true; }
    }
    if (!ds.frames.length) throw err('No atomic positions were found.');
    settleCells(ds, true);
    const v = /Program PWSCF\s+(\S+)/.exec((Array.isArray(text) ? text.slice(0, 200).join('\n') : text).slice(0, 4000));
    if (v) ds.meta.version = v[1];
    return ds;
  }

  // ================================================================== CASTEP .cell
  const LEN_UNIT = { ang: 1, angstrom: 1, bohr: BOHR, a0: BOHR, nm: 10, m: 1e10, cm: 1e8 };
  function takeUnit(lines) {
    if (lines.length && /^[a-z]+0?$/i.test(lines[0])) { const u = lines.shift().toLowerCase(); return LEN_UNIT[u] !== undefined ? LEN_UNIT[u] : 1; }
    return 1;
  }
  function readCastepCell(text, fileName) {
    const blocks = {};
    let cur = null;
    for (const raw of splitLines(text)) {
      const s = raw.replace(/[!#].*$/, '').trim();
      if (!s) continue;
      let m = /^%block\s+(\w+)/i.exec(s);
      if (m) { cur = blocks[m[1].toLowerCase()] = []; continue; }
      if (/^%endblock/i.test(s)) { cur = null; continue; }
      if (cur) cur.push(s);
    }
    const ds = newDs('castep-cell', fileName, 'CASTEP');
    let cell = null;
    if (blocks.lattice_cart) {
      const b = blocks.lattice_cart.slice(), u = takeUnit(b);
      cell = b.slice(0, 3).map(s => tok(s).slice(0, 3).map(x => num(x) * u));
    } else if (blocks.lattice_abc) {
      const b = blocks.lattice_abc.slice(), u = takeUnit(b);
      const p = tok(b[0]).map(num), q = tok(b[1]).map(num);
      if (p.length >= 3 && q.length >= 3) cell = cellFromParams(p[0] * u, p[1] * u, p[2] * u, q[0], q[1], q[2]);
    }
    if (!cell || !goodCell(cell)) throw err('No valid LATTICE_CART or LATTICE_ABC block.');
    const frac = blocks.positions_frac;
    const pos = frac ? frac.slice() : blocks.positions_abs ? blocks.positions_abs.slice() : null;
    if (!pos) throw err('No POSITIONS_FRAC or POSITIONS_ABS block.');
    const u = frac ? 1 : takeUnit(pos);
    const Z = [], xyz = [], names = [];
    let unknown = 0;
    for (const s of pos) {
      const t = tok(s);
      const x = num(t[1]), y = num(t[2]), z = num(t[3]);
      if (!Number.isFinite(x + y + z)) throw err('Bad position line: "' + s + '".');
      let el = isInt(t[0]) ? parseInt(t[0], 10) : zOf(t[0]);
      if (el < 0) { unknown++; el = 0; }
      const r = frac ? f2c(cell, x, y, z) : [x * u, y * u, z * u];
      Z.push(el); xyz.push(r[0], r[1], r[2]); names.push(t[0]);
    }
    const f = newFrame(Z, xyz);
    f.atomNames = names;
    ds.frames.push(f);
    ds.cell = cell;
    if (unknown) ds.warnings.push(unknown + ' atom(s) with an unrecognised species are shown as dummy atoms.');
    return ds;
  }

  // ================================================================== XCrySDen XSF / AXSF
  function readXsf(text, fileName) {
    const L = splitLines(text);
    const ds = newDs('xsf', fileName, 'XCrySDen');
    let cell = null, periodic = false;
    const vols = [];
    const zTok = s => { const z = isInt(s) ? parseInt(s, 10) : zOf(s); return z >= 0 ? z : 0; };
    const atomLine = t => t.length >= 4 && (isInt(t[0]) || /^[A-Za-z]{1,2}\d*$/.test(t[0])) && isNum(t[1]) && isNum(t[2]) && isNum(t[3]);
    const addFrame = (rows, label) => {
      const Z = [], xyz = [];
      for (const t of rows) { Z.push(zTok(t[0])); xyz.push(num(t[1]), num(t[2]), num(t[3])); }
      const f = newFrame(Z, xyz);
      if (periodic && cell) f.cell = cell;
      if (label) f.label = label;
      ds.frames.push(f);
    };
    let i = 0;
    while (i < L.length) {
      const t = tok(L[i]);
      if (!t.length || t[0][0] === '#') { i++; continue; }
      const kw = t[0].toUpperCase();
      if (kw === 'CRYSTAL' || kw === 'SLAB' || kw === 'POLYMER') { periodic = true; i++; continue; }
      if (kw === 'MOLECULE') { periodic = false; i++; continue; }
      if (kw === 'PRIMVEC') {
        const c = [1, 2, 3].map(k => tok(L[i + k]).slice(0, 3).map(num));
        if (!c.every(v => v.length === 3 && fin(v))) throw err('PRIMVEC needs three lattice vectors.');
        cell = c; i += 4; continue;
      }
      if (kw === 'CONVVEC') { i += 4; continue; }
      if (kw === 'PRIMCOORD' || kw === 'CONVCOORD') {
        const n = parseInt(tok(L[i + 1])[0], 10);
        if (!(n > 0)) throw err(kw + ' must be followed by the number of atoms.');
        if (kw === 'PRIMCOORD') {
          const rows = [];
          for (let k = 0; k < n; k++) { const r = tok(L[i + 2 + k]); if (!atomLine(r)) throw err('PRIMCOORD block is truncated.'); rows.push(r); }
          addFrame(rows, t[1] ? 'Step ' + t[1] : null);
        }
        i += 2 + n; continue;
      }
      if (kw === 'ATOMS') {
        const rows = [];
        i++;
        while (i < L.length && atomLine(tok(L[i]))) rows.push(tok(L[i++]));
        if (rows.length) addFrame(rows, t[1] ? 'Step ' + t[1] : null);
        continue;
      }
      if (/^BEGIN_BLOCK_DATAGRID_?3D/.test(kw)) { i++; if (i < L.length && !/^\s*(BEGIN_)?DATAGRID_3D/i.test(L[i])) i++; continue; }
      if (/^BEGIN_BLOCK_DATAGRID_?2D/.test(kw)) { while (i < L.length && !/^\s*END_BLOCK_DATAGRID_?2D/i.test(L[i])) i++; i++; continue; }
      if (/^(BEGIN_)?DATAGRID_3D/.test(kw)) {
        const name = t[0].replace(/^(BEGIN_)?DATAGRID_3D_?/i, '') || 'grid';
        const dims = tok(L[i + 1]).slice(0, 3).map(x => parseInt(x, 10));
        const o = tok(L[i + 2]).slice(0, 3).map(num), sp = [3, 4, 5].map(k => tok(L[i + k]).slice(0, 3).map(num));
        if (!(dims.length === 3 && dims.every(d => d > 0)) || o.length < 3 || !fin(o) || !sp.every(v => v.length === 3 && fin(v))) throw err('Bad DATAGRID_3D header.');
        const [n1, n2, n3] = dims, total = n1 * n2 * n3, raw = new Float32Array(total);
        let cnt = 0, j = i + 6;
        for (; j < L.length; j++) {
          if (/^\s*END_DATAGRID/i.test(L[j])) break;
          for (const x of tok(L[j])) { if (cnt < total) raw[cnt] = num(x); cnt++; }
        }
        if (cnt < total) throw err('DATAGRID_3D "' + name + '" is truncated: ' + cnt + ' of ' + total + ' values.');
        const data = new Float32Array(total);                 // XSF: x fastest -> BlueMolecule: k (axis 3) fastest
        for (let c = 0, q = 0; c < n3; c++) for (let b = 0; b < n2; b++) for (let a = 0; a < n1; a++, q++) data[(a * n2 + b) * n3 + c] = raw[q];
        const axes = sp.map((v, k) => v.map(x => x / Math.max(1, dims[k] - 1)));
        vols.push({ name, kind: /dens|rho|charge/i.test(name) ? 'density' : /orb|wfn|psi|mo/i.test(name) ? 'orbital' : 'other', origin: o, axes, dims, data });
        i = j + 1; continue;
      }
      i++;
    }
    if (!ds.frames.length && !vols.length) throw err('No PRIMCOORD or ATOMS block found.');
    if (vols.length) ds.volumes = vols;
    if (periodic) settleCells(ds, false);
    if (periodic && !ds.cell && cell) ds.cell = cell;
    return ds;
  }

  // ================================================================== CP2K input (&SUBSYS &CELL/&COORD)
  const CP2K_UNIT = { angstrom: 1, ang: 1, bohr: BOHR, nm: 10, pm: 0.01 };
  function readCp2kIn(text, fileName) {
    const ds = newDs('cp2k-in', fileName, 'CP2K');
    const stack = [], cellK = {}, rows = [], kinds = {};
    let unitC = 1, scaled = false, kind = null, topo = null;
    const unitTok = (v) => { if (v.length && /^\[.*\]$/.test(v[0])) { const u = v[0].slice(1, -1).toLowerCase(); v.shift(); return CP2K_UNIT[u] !== undefined ? CP2K_UNIT[u] : 1; } return 1; };
    for (const raw of splitLines(text)) {
      const s = raw.replace(/[!#].*$/, '').trim();
      if (!s) continue;
      if (/^&END\b/i.test(s)) { stack.pop(); continue; }
      const m = /^&(\w+)\s*(.*)$/.exec(s);
      if (m) { const nm = m[1].toUpperCase(); stack.push(nm); if (nm === 'KIND') kind = m[2].trim(); continue; }
      const top = stack[stack.length - 1], parent = stack[stack.length - 2], t = tok(s), key = t[0].toUpperCase();
      if (top === 'CELL' && parent === 'SUBSYS') {
        const v = t.slice(1), u = unitTok(v);
        if (key === 'A' || key === 'B' || key === 'C' || key === 'ABC') cellK[key] = v.slice(0, 3).map(x => num(x) * u);
        else if (key === 'ALPHA_BETA_GAMMA') cellK.ang = v.slice(0, 3).map(num);
      } else if (top === 'COORD' && parent === 'SUBSYS') {
        if (key === 'UNIT') unitC = CP2K_UNIT[(t[1] || '').toLowerCase()] || 1;
        else if (key === 'SCALED') scaled = !t[1] || /^(t|true|\.true\.|yes|on)$/i.test(t[1]);
        else if (t.length >= 4) rows.push(t);
      } else if (top === 'KIND' && key === 'ELEMENT' && kind) kinds[kind.toUpperCase()] = t[1];
      else if (top === 'TOPOLOGY' && key === 'COORD_FILE_NAME') topo = t[1];
      else if (top === 'GLOBAL' && (key === 'PROJECT' || key === 'PROJECT_NAME')) ds.title = t[1] || ds.title;
    }
    let cell = null;
    if (cellK.A && cellK.B && cellK.C) cell = [cellK.A, cellK.B, cellK.C];
    else if (cellK.ABC) { const g = cellK.ang || [90, 90, 90]; cell = cellFromParams(cellK.ABC[0], cellK.ABC[1], cellK.ABC[2], g[0], g[1], g[2]); }
    if (cell && !goodCell(cell)) cell = null;
    if (!rows.length) {
      if (topo) throw err('The coordinates are in the external file "' + topo + '"; open that file instead.');
      throw err('No &COORD section with atoms found in &SUBSYS.');
    }
    if (scaled && !cell) throw err('SCALED coordinates need a &CELL.');
    const Z = [], xyz = [], names = [];
    let unknown = 0;
    for (const t of rows) {
      const x = num(t[1]), y = num(t[2]), z = num(t[3]);
      if (!Number.isFinite(x + y + z)) throw err('Bad &COORD line: "' + t.join(' ') + '".');
      let el = zOf(kinds[t[0].toUpperCase()] || t[0]);
      if (el < 0) { unknown++; el = 0; }
      const r = scaled ? f2c(cell, x, y, z) : [x * unitC, y * unitC, z * unitC];
      Z.push(el); xyz.push(r[0], r[1], r[2]); names.push(t[0]);
    }
    const f = newFrame(Z, xyz);
    f.atomNames = names;
    ds.frames.push(f);
    ds.cell = cell;
    if (unknown) ds.warnings.push(unknown + ' atom(s) with an unrecognised kind are shown as dummy atoms.');
    return ds;
  }

  // ================================================================== Amber prmtop / inpcrd
  function readPrmtop(text, fileName) {
    const flags = {};
    let cur = null;
    for (const s of splitLines(text)) {
      if (s.startsWith('%FLAG')) { cur = flags[tok(s.slice(5))[0]] = { fmt: null, lines: [] }; continue; }
      if (s.startsWith('%FORMAT')) { if (cur) cur.fmt = /\((\d+)([aAiIeEfF])(\d+)/.exec(s); continue; }
      if (s.startsWith('%')) continue;
      if (cur) cur.lines.push(s);
    }
    const get = (name, count) => {
      const f = flags[name];
      if (!f || !f.fmt) return null;
      const w = parseInt(f.fmt[3], 10), str = /a/i.test(f.fmt[2]), out = [];
      for (const s of f.lines) for (let p = 0; p < s.length; p += w) {
        const v = s.substr(p, w);
        if (str) out.push(v.trim()); else if (v.trim()) out.push(num(v));
      }
      return count !== undefined ? out.slice(0, count) : out;
    };
    const P = get('POINTERS');
    if (!P || P.length < 12) throw err('No %FLAG POINTERS section.');
    const n = P[0], nres = P[11], ifbox = P[27] || 0;
    const ds = newDs('prmtop', fileName, 'Amber');
    const names = get('ATOM_NAME', n) || [];
    const anum = get('ATOMIC_NUMBER', n), mass = get('MASS', n), charge = get('CHARGE', n);
    const rlab = get('RESIDUE_LABEL', nres) || [], rptr = get('RESIDUE_POINTER', nres) || [];
    const rnum = get('RESIDUE_NUMBER', nres), rch = get('RESIDUE_CHAINID', nres);
    const res = new Array(n).fill(''), rid = new Array(n).fill(0), ch = rch ? new Array(n).fill('') : null;
    for (let r = 0; r < nres; r++) {
      const a0 = rptr[r] - 1, a1 = r + 1 < nres ? rptr[r + 1] - 1 : n;
      for (let a = Math.max(0, a0); a < Math.min(n, a1); a++) { res[a] = rlab[r] || ''; rid[a] = rnum ? rnum[r] : r + 1; if (ch) ch[a] = rch[r] || ''; }
    }
    const useNum = anum && anum.length === n && anum.some(z => z > 0);
    const Z = [];
    for (let a = 0; a < n; a++) Z.push(useNum ? Math.max(0, anum[a] | 0) : zFromMassName(mass ? mass[a] : 0, names[a], res[a]));
    if (!useNum) ds.warnings.push('No ATOMIC_NUMBER section: elements were assigned from the atomic masses.');
    const bonds = [];
    for (const key of ['BONDS_INC_HYDROGEN', 'BONDS_WITHOUT_HYDROGEN']) {
      const b = get(key) || [];
      for (let k = 0; k + 2 < b.length; k += 3) { const i = Math.abs(b[k]) / 3, j = Math.abs(b[k + 1]) / 3; if (i < n && j < n) bonds.push([i, j, 1]); }
    }
    const f = newFrame(Z, new Float64Array(3 * n));
    f.atomNames = names; f.resNames = res; f.resIds = rid; f.chains = ch;
    ds.frames.push(f);
    ds.bonds = bonds;
    if (charge && charge.length === n) ds.charges = { Amber: charge.map(q => q / AMBER_Q) };
    const box = get('BOX_DIMENSIONS');
    if (ifbox > 0 && box && box.length >= 4 && box[1] > 0) {
      const be = box[0];
      ds.cell = Math.abs(be - 109.4712206) < 0.01 ? cellFromParams(box[1], box[2], box[3], be, be, be) : cellFromParams(box[1], box[2], box[3], 90, be, 90);
    }
    const tf = flags.TITLE || flags.CTITLE;
    const title = tf ? tf.lines.join(' ').trim() : '';
    if (title) ds.title = title;
    ds.meta.topologyOnly = true;
    ds.meta.natoms = n;
    return ds;
  }
  function readInpcrd(text, fileName) {
    const L = splitLines(text);
    const h = tok(L[1]);
    const n = parseInt(h[0], 10);
    if (!(n > 0)) throw err('Line 2 must hold the number of atoms.');
    let vals = [];
    for (let k = 2; k < L.length; k++) {
      const s = L[k];
      if (!s.trim()) continue;
      for (let p = 0; p < s.length; p += 12) { const f = s.substr(p, 12); if (f.trim()) vals.push(num(f)); }
    }
    if (!fin(vals)) { vals = []; for (let k = 2; k < L.length; k++) for (const x of tok(L[k])) vals.push(num(x)); }
    if (vals.length < 3 * n || !fin(vals)) throw err('Expected ' + 3 * n + ' coordinates, found ' + vals.length + '.');
    const ds = newDs('inpcrd', fileName, 'Amber');
    const extra = vals.length - 3 * n;
    const boxN = extra >= 3 * n ? extra - 3 * n : extra;
    const f = newFrame(null, vals.slice(0, 3 * n));
    if (Number.isFinite(num(h[1]))) f.time = num(h[1]);
    ds.frames.push(f);
    if (boxN === 6 || boxN === 3) {
      const b = vals.slice(vals.length - boxN);
      ds.cell = boxN === 6 ? cellFromParams(b[0], b[1], b[2], b[3], b[4], b[5]) : [[b[0], 0, 0], [0, b[1], 0], [0, 0, b[2]]];
    }
    if (extra >= 3 * n) ds.meta.hasVelocities = true;
    ds.title = L[0].trim() || ds.title;
    ds.meta.needsTopology = true;
    ds.meta.natoms = n;
    return ds;
  }
  function sniffInpcrd(h) {
    const L = h.split(/\r\n|\r|\n/, 5);
    if (L.length < 3 || !/^\s*\d+(\s+[-+\d.eE]+)?\s*$/.test(L[1])) return 0;
    return /^(\s*-?\d+\.\d{7}){1,6}\s*$/.test(L[2]) ? 0.45 : 0;
  }

  // ================================================================== CHARMM / NAMD PSF
  function readPsf(text, fileName) {
    const L = splitLines(text);
    if (!/^\s*PSF\b/.test(L[0] || '')) throw err('The first line must start with "PSF".');
    const ds = newDs('psf', fileName, 'CHARMM');
    let i = 1, n = -1;
    const remarks = [];
    for (; i < L.length; i++) {
      if (/^\s*REMARKS/i.test(L[i])) remarks.push(L[i].replace(/^\s*REMARKS\s*/i, '').trim());
      if (/!NATOM/.test(L[i])) { n = parseInt(tok(L[i])[0], 10); i++; break; }
    }
    if (!(n > 0)) throw err('No !NATOM section.');
    const Z = [], names = [], res = [], rid = [], seg = [], q = [];
    for (let k = 0; k < n; i++) {
      if (i >= L.length) throw err('The atom section is truncated.');
      const t = tok(L[i]);
      if (!t.length) continue;
      if (t.length < 8) throw err('Bad atom line ' + (i + 1) + '.');
      seg.push(t[1]); rid.push(parseInt(t[2], 10) || 0); res.push(t[3]); names.push(t[4]); q.push(num(t[6]));
      Z.push(zFromMassName(num(t[7]), t[4], t[3]));
      k++;
    }
    const bonds = [];
    for (; i < L.length; i++) if (/!NBOND/.test(L[i])) break;
    if (i < L.length) {
      const nb = parseInt(tok(L[i])[0], 10), ints = [];
      for (i++; i < L.length && ints.length < 2 * nb; i++) for (const x of tok(L[i])) { if (!isInt(x)) break; ints.push(parseInt(x, 10)); }
      for (let k = 0; k + 1 < ints.length; k += 2) if (ints[k] >= 1 && ints[k + 1] >= 1 && ints[k] <= n && ints[k + 1] <= n) bonds.push([ints[k] - 1, ints[k + 1] - 1, 1]);
    }
    const f = newFrame(Z, new Float64Array(3 * n));
    f.atomNames = names; f.resNames = res; f.resIds = rid; f.chains = seg;
    ds.frames.push(f);
    ds.bonds = bonds;
    if (fin(q)) ds.charges = { CHARMM: q };
    const title = remarks.find(r => r && !/^(original|topology|segment|patch|psfgen|VMD)/i.test(r));
    if (title) ds.title = title;
    ds.meta.topologyOnly = true;
    ds.meta.natoms = n;
    return ds;
  }

  // ================================================================== DCD (CHARMM / NAMD / X-PLOR)
  function readDcd(u8, fileName) {
    const dv = new DataView(u8.buffer, u8.byteOffset, u8.byteLength), N = u8.length;
    const cord = p => u8[p] === 67 && u8[p + 1] === 79 && u8[p + 2] === 82 && u8[p + 3] === 68;
    let le, M;
    if (N >= 8 && dv.getInt32(0, true) === 84 && cord(4)) { le = true; M = 4; }
    else if (N >= 8 && dv.getInt32(0, false) === 84 && cord(4)) { le = false; M = 4; }
    else if (N >= 12 && dv.getInt32(0, true) === 84 && dv.getInt32(4, true) === 0 && cord(8)) { le = true; M = 8; }
    else if (N >= 12 && dv.getInt32(0, false) === 0 && dv.getInt32(4, false) === 84 && cord(8)) { le = false; M = 8; }
    else throw err('Not a DCD file (no CORD header record).');
    const i32 = p => dv.getInt32(p, le);
    const mark = p => (M === 4 ? i32(p) : (le ? i32(p) : i32(p + 4)));
    let pos = 0;
    const rec = () => {
      if (pos + M > N) return null;
      const len = mark(pos), start = pos + M, end = start + len;
      if (len < 0 || end + M > N || mark(end) !== len) return null;
      pos = end + M;
      return { start, len };
    };
    const h = rec();
    if (!h || h.len < 84) throw err('The DCD header record is damaged.');
    const ic = []; for (let k = 0; k < 20; k++) ic.push(i32(h.start + 4 + 4 * k));
    const charmm = ic[19] !== 0, namnf = ic[8], hasCell = charmm && ic[10] !== 0, has4D = charmm && ic[11] !== 0;
    const tr = rec();
    let title = '';
    if (tr && tr.len >= 4) {
      const nt = i32(tr.start);
      for (let k = 0; k < nt && 4 + 80 * (k + 1) <= tr.len; k++) title += String.fromCharCode.apply(null, u8.subarray(tr.start + 4 + 80 * k, tr.start + 4 + 80 * (k + 1))).replace(/\0/g, ' ').trim() + ' ';
    }
    const nr = rec();
    if (!nr || nr.len < 4) throw err('The DCD atom-count record is missing.');
    const n = i32(nr.start);
    if (!(n > 0)) throw err('The DCD file declares no atoms.');
    let free = null;
    if (namnf > 0) {
      const fr = rec();
      if (!fr || fr.len !== 4 * (n - namnf)) throw err('The DCD free-atom list is damaged.');
      free = new Int32Array(n - namnf);
      for (let k = 0; k < free.length; k++) free[k] = i32(fr.start + 4 * k) - 1;
    }
    const ds = newDs('dcd', fileName, charmm ? 'CHARMM/NAMD' : 'X-PLOR');
    let last = null, truncated = false;
    while (pos < N) {
      let cell = null;
      if (hasCell) {
        const c = rec();
        if (!c || c.len !== 48) { truncated = true; break; }
        const u = []; for (let k = 0; k < 6; k++) u.push(dv.getFloat64(c.start + 8 * k, le));
        // unitcell = [A, gamma, B, beta, alpha, C]; angles are cosines (CHARMM, NAMD >= 2.5) or degrees
        const cosines = [u[1], u[3], u[4]].every(x => x >= -1 && x <= 1);
        const ang = x => (cosines ? Math.acos(x) / RAD : x);
        if (u[0] > 0 && u[2] > 0 && u[5] > 0) cell = boxCell(u[0], u[2], u[5], ang(u[4]), ang(u[3]), ang(u[1]));
      }
      const cnt = free && last ? free.length : n;
      const xs = [rec(), rec(), rec()];
      if (xs.some(r => !r || r.len !== 4 * cnt)) { truncated = true; break; }
      if (has4D && !rec()) { truncated = true; break; }
      const xyz = last ? Float64Array.from(last) : new Float64Array(3 * n);
      for (let d = 0; d < 3; d++) {
        const st = xs[d].start;
        for (let k = 0; k < cnt; k++) xyz[3 * (cnt === n ? k : free[k]) + d] = dv.getFloat32(st + 4 * k, le);
      }
      const f = newFrame(null, xyz);
      if (cell) f.cell = cell;
      ds.frames.push(f);
      last = xyz;
    }
    if (!ds.frames.length) throw err('The DCD file holds no complete frame.');
    if (truncated) ds.warnings.push('The last frame is incomplete and was skipped.');
    if (ic[0] > 0 && ic[0] !== ds.frames.length) ds.warnings.push('The header announces ' + ic[0] + ' frames; ' + ds.frames.length + ' were read.');
    settleCells(ds, false);
    if (title.trim()) ds.title = title.trim().slice(0, 120);
    ds.meta.needsTopology = true;
    ds.meta.natoms = n;
    return ds;
  }

  // ================================================================== LAMMPS data / dump
  const LMP_STYLE = {
    atomic: { type: 1, x: 2 }, charge: { type: 1, q: 2, x: 3 }, molecular: { mol: 1, type: 2, x: 3 },
    bond: { mol: 1, type: 2, x: 3 }, angle: { mol: 1, type: 2, x: 3 }, dihedral: { mol: 1, type: 2, x: 3 },
    full: { mol: 1, type: 2, q: 3, x: 4 }, sphere: { type: 1, x: 4 }, ellipsoid: { type: 1, x: 4 }, body: { type: 1, x: 4 },
    peri: { type: 1, x: 4 }, dipole: { type: 1, q: 2, x: 3 }, hybrid: { type: 1, x: 2 }
  };
  function readLammpsData(text, fileName) {
    const L = splitLines(text);
    const ds = newDs('lammps-data', fileName, 'LAMMPS');
    const box = { xlo: 0, xhi: 0, ylo: 0, yhi: 0, zlo: 0, zhi: 0, xy: 0, xz: 0, yz: 0 }, gen = {};
    let natoms = 0, i = 1, hasBox = false;
    for (; i < L.length; i++) {
      const s = L[i].replace(/#.*$/, '').trim();
      if (!s) continue;
      const t = tok(s);
      if (!isNum(t[0])) break;
      const words = t.filter(x => !isNum(x)).join(' ').toLowerCase(), v = t.filter(isNum).map(num);
      if (words === 'atoms') natoms = v[0];
      else if (words === 'xlo xhi') { box.xlo = v[0]; box.xhi = v[1]; hasBox = true; }
      else if (words === 'ylo yhi') { box.ylo = v[0]; box.yhi = v[1]; }
      else if (words === 'zlo zhi') { box.zlo = v[0]; box.zhi = v[1]; }
      else if (words === 'xy xz yz') { box.xy = v[0]; box.xz = v[1]; box.yz = v[2]; }
      else if (words === 'avec' || words === 'bvec' || words === 'cvec') gen[words[0]] = v.slice(0, 3);
      else if (words === 'abc origin') gen.o = v.slice(0, 3);
    }
    const sec = {};
    let cur = null;
    for (; i < L.length; i++) {
      const raw = L[i], s = raw.replace(/#.*$/, '').trim();
      if (!s) continue;
      if (/^[A-Za-z]/.test(s)) { cur = sec[s.toLowerCase()] = { head: raw, rows: [] }; continue; }
      if (cur) cur.rows.push(raw);
    }
    const atoms = sec.atoms;
    if (!atoms || !atoms.rows.length) throw err('No Atoms section found.');
    const masses = {}, labels = {}, labelId = {};
    for (const r of (sec['atom type labels'] || { rows: [] }).rows) { const t = tok(r.replace(/#.*$/, '')); if (t.length >= 2) { labels[t[0]] = t[1]; labelId[t[1]] = t[0]; } }
    for (const r of (sec.masses || { rows: [] }).rows) {
      const t = tok(r.replace(/#.*$/, '')), c = /#\s*(\S+)/.exec(r);
      if (t.length >= 2) { const id = labelId[t[0]] || t[0]; masses[id] = { m: num(t[1]), label: c ? c[1] : labels[id] || null }; }
    }
    const rows = atoms.rows.map(r => tok(r.replace(/#.*$/, '')));
    const hint = (/#\s*(\w+)/.exec(atoms.head) || [])[1];
    const w = rows[0].length;
    let style = hint && LMP_STYLE[hint.toLowerCase()] ? hint.toLowerCase() : null;
    if (!style) {
      const base = w >= 8 && rows.every(t => t.slice(-3).every(isInt)) ? w - 3 : w;
      if (base === 5) style = 'atomic';
      else if (base === 7) style = 'full';
      else if (base === 6) style = rows.some(t => /[.eE]/.test(t[2])) ? 'charge' : 'molecular';
      else throw err('Cannot tell the atom style of the Atoms section (' + w + ' columns); add a "# style" comment.');
    }
    const lay = LMP_STYLE[style];
    if (natoms && rows.length !== natoms) ds.warnings.push('The header announces ' + natoms + ' atoms; ' + rows.length + ' were read.');
    const recs = rows.map(t => {
      if (t.length < lay.x + 3) throw err('Atom line "' + t.join(' ') + '" is too short for atom style ' + style + '.');
      const ty = labelId[t[lay.type]] || t[lay.type];
      return { id: parseInt(t[0], 10), type: ty, mol: lay.mol ? parseInt(t[lay.mol], 10) : 0, q: lay.q ? num(t[lay.q]) : 0,
        x: num(t[lay.x]), y: num(t[lay.x + 1]), z: num(t[lay.x + 2]) };
    }).sort((a, b) => a.id - b.id);
    const n = recs.length, Z = [], xyz = [], names = [], mol = [], q = [], idx = new Map();
    let noMass = 0;
    recs.forEach((r, k) => {
      if (!Number.isFinite(r.x + r.y + r.z)) throw err('Bad coordinates for atom ' + r.id + '.');
      idx.set(r.id, k);
      const ms = masses[r.type];
      let z = ms ? zFromMass(ms.m) : -1;
      if (z < 0) { const lab = labels[r.type]; z = lab ? Math.max(0, zOf(lab)) : 0; noMass++; }
      Z.push(z); xyz.push(r.x, r.y, r.z); names.push(ms && ms.label ? ms.label : labels[r.type] || String(r.type)); mol.push(r.mol); q.push(r.q);
    });
    const f = newFrame(Z, xyz);
    f.atomNames = names;
    if (lay.mol) { f.resIds = mol; f.resNames = new Array(n).fill('MOL'); }
    ds.frames.push(f);
    if (lay.q) ds.charges = { LAMMPS: q };
    const bs = sec.bonds;
    if (bs) {
      ds.bonds = [];
      for (const r of bs.rows) { const t = tok(r.replace(/#.*$/, '')); const a = idx.get(parseInt(t[2], 10)), b = idx.get(parseInt(t[3], 10)); if (a !== undefined && b !== undefined) ds.bonds.push([a, b, 1]); }
    }
    if (gen.a && gen.b && gen.c) { ds.cell = [gen.a, gen.b, gen.c]; ds.cellOrigin = gen.o || [0, 0, 0]; }
    else if (hasBox) {
      ds.cell = [[box.xhi - box.xlo, 0, 0], [box.xy, box.yhi - box.ylo, 0], [box.xz, box.yz, box.zhi - box.zlo]];
      ds.cellOrigin = [box.xlo, box.ylo, box.zlo];
    }
    if (noMass) ds.warnings.push(noMass + ' atom(s) have no Masses entry; their elements come from type labels or are dummies.');
    ds.title = (L[0] || '').trim() || ds.title;
    ds.meta.atomStyle = style;
    return ds;
  }
  function readLammpsDump(text, fileName) {
    const L = splitLines(text);
    const ds = newDs('lammps-dump', fileName, 'LAMMPS');
    let i = 0, step = null, n = 0, box = null, needTop = false;
    while (i < L.length) {
      const s = L[i];
      if (s.startsWith('ITEM: TIMESTEP')) { step = parseInt(L[i + 1], 10); i += 2; continue; }
      if (s.startsWith('ITEM: NUMBER OF ATOMS')) { n = parseInt(L[i + 1], 10); i += 2; continue; }
      if (s.startsWith('ITEM: BOX BOUNDS')) {
        const tri = /xy\s+xz\s+yz/.test(s), b = [1, 2, 3].map(k => tok(L[i + k]).map(num));
        if (!b.every(r => r.length >= 2 && fin(r))) throw err('Bad BOX BOUNDS block.');
        let [xlo, xhi] = b[0], [ylo, yhi] = b[1];
        const [zlo, zhi] = b[2];
        let xy = 0, xz = 0, yz = 0;
        if (tri) {
          xy = b[0][2] || 0; xz = b[1][2] || 0; yz = b[2][2] || 0;
          xlo -= Math.min(0, xy, xz, xy + xz); xhi -= Math.max(0, xy, xz, xy + xz);
          ylo -= Math.min(0, yz); yhi -= Math.max(0, yz);
        }
        box = { lo: [xlo, ylo, zlo], cell: [[xhi - xlo, 0, 0], [xy, yhi - ylo, 0], [xz, yz, zhi - zlo]] };
        i += 4; continue;
      }
      if (s.startsWith('ITEM: ATOMS')) {
        const cols = tok(s.slice(11)), c = k => cols.indexOf(k);
        let px = ['x', 'y', 'z'].map(c), mode = 'abs';
        if (px[0] < 0) { px = ['xu', 'yu', 'zu'].map(c); }
        if (px[0] < 0) { px = ['xs', 'ys', 'zs'].map(c); mode = 'scaled'; }
        if (px[0] < 0) { px = ['xsu', 'ysu', 'zsu'].map(c); mode = 'scaled'; }
        if (px.some(k => k < 0)) throw err('The ATOMS columns hold no x/y/z, xu, xs or xsu coordinates.');
        if (mode === 'scaled' && !box) throw err('Scaled coordinates need a BOX BOUNDS block.');
        const ci = c('id'), ce = c('element'), cm = c('mass'), ct = c('type');
        const rows = [];
        for (let k = 0; k < n; k++) {
          const t = tok(L[i + 1 + k]);
          if (t.length < cols.length) { ds.warnings.push('Timestep ' + step + ' is incomplete and was skipped.'); i = L.length; break; }
          rows.push(t);
        }
        if (i >= L.length) break;
        if (ci >= 0) { const ids = new Map(rows.map(t => [t, parseInt(t[ci], 10)])); rows.sort((a, b) => ids.get(a) - ids.get(b)); }
        const xyz = new Float64Array(3 * n);
        rows.forEach((t, k) => {
          let x = num(t[px[0]]), y = num(t[px[1]]), z = num(t[px[2]]);
          if (mode === 'scaled') { const r = f2c(box.cell, x, y, z); x = r[0] + box.lo[0]; y = r[1] + box.lo[1]; z = r[2] + box.lo[2]; }
          xyz[3 * k] = x; xyz[3 * k + 1] = y; xyz[3 * k + 2] = z;
        });
        let Z = null;
        if (ce >= 0) Z = rows.map(t => Math.max(0, zOf(t[ce])));
        else if (cm >= 0) Z = rows.map(t => zFromMass(num(t[cm])));
        else needTop = true;
        const f = newFrame(needTop ? null : Z, xyz);
        f.label = 'Timestep ' + step;
        f.step = step;
        if (box) f.cell = box.cell;
        if (ct >= 0) f.types = rows.map(t => t[ct]);
        ds.frames.push(f);
        i += 1 + n; continue;
      }
      i++;
    }
    if (!ds.frames.length) throw err('No "ITEM: ATOMS" blocks were found.');
    if (needTop) {
      for (const f of ds.frames) f.Z = null;
      ds.meta.needsTopology = true;
      ds.meta.natoms = ds.frames[0].xyz.length / 3;
      ds.meta.atomTypes = ds.frames[0].types || null;
    }
    if (box) ds.cellOrigin = box.lo;
    settleCells(ds, false);
    return ds;
  }

  // ================================================================== GROMACS TRR / XTC (XDR, big-endian)
  function readTrr(u8, fileName) {
    const dv = new DataView(u8.buffer, u8.byteOffset, u8.byteLength), N = u8.length;
    const ds = newDs('trr', fileName, 'GROMACS');
    let p = 0, natoms = 0, truncated = false;
    while (p + 4 <= N) {
      if (dv.getInt32(p) !== 1993) { if (!ds.frames.length) throw err('Not a TRR file (bad magic number).'); truncated = true; break; }
      if (p + 12 > N) { truncated = true; break; }
      const slen = dv.getInt32(p + 8);
      let q = p + 12 + ((slen + 3) & ~3);
      if (q + 13 * 4 > N) { truncated = true; break; }
      const H = []; for (let k = 0; k < 13; k++) H.push(dv.getInt32(q + 4 * k));
      q += 52;
      const [irS, eS, boxS, virS, presS, topS, symS, xS, vS, fS, nat, step] = H;
      const ref = boxS ? boxS / 9 : xS ? xS / (3 * nat) : vS ? vS / (3 * nat) : fS ? fS / (3 * nat) : 4;
      const R = ref === 8 ? 8 : 4, rd = o => (R === 8 ? dv.getFloat64(o) : dv.getFloat32(o));
      if (slen < 0 || H.slice(0, 11).some(v => v < 0) || (xS && xS !== 3 * nat * R) || (boxS && boxS !== 9 * R)) {   // corrupt header
        if (!ds.frames.length) throw err('The TRR frame header is corrupt (negative or inconsistent block sizes).');
        truncated = true; break;
      }
      const end = q + 2 * R + irS + eS + boxS + virS + presS + topS + symS + xS + vS + fS;
      if (end > N) { truncated = true; break; }
      const t = rd(q);
      q += 2 * R + irS + eS;
      let cell = null;
      if (boxS) { const b = []; for (let k = 0; k < 9; k++) b.push(rd(q + R * k) * 10); cell = [b.slice(0, 3), b.slice(3, 6), b.slice(6, 9)]; if (!goodCell(cell)) cell = null; }
      q += boxS + virS + presS + topS + symS;
      if (xS) {
        if (natoms && nat !== natoms) throw err('The atom count changes between TRR frames.');
        natoms = nat;
        const xyz = new Float64Array(3 * nat);
        for (let k = 0; k < 3 * nat; k++) xyz[k] = rd(q + R * k) * 10;
        const f = newFrame(null, xyz);
        f.time = t; f.step = step; f.label = 't = ' + (+t.toFixed(4)) + ' ps';
        if (cell) f.cell = cell;
        ds.frames.push(f);
      }
      p = end;
    }
    if (!ds.frames.length) throw err('The TRR file holds no frame with coordinates.');
    if (truncated) ds.warnings.push('The file ends with an incomplete frame, which was skipped.');
    settleCells(ds, false);
    ds.meta.needsTopology = true;
    ds.meta.natoms = natoms;
    return ds;
  }

  // XTC coordinate decompression (the xdr3dfcoord scheme of xdrfile / GROMACS).
  const MAGICINTS = [0, 0, 0, 0, 0, 0, 0, 0, 0, 8, 10, 12, 16, 20, 25, 32, 40, 50, 64, 80, 101, 128, 161, 203, 256, 322, 406, 512,
    645, 812, 1024, 1290, 1625, 2048, 2580, 3250, 4096, 5060, 6501, 8192, 10321, 13003, 16384, 20642, 26007, 32768, 41285, 52015,
    65536, 82570, 104031, 131072, 165140, 208063, 262144, 330280, 416127, 524287, 660561, 832255, 1048576, 1321122, 1664510,
    2097152, 2642245, 3329021, 4194304, 5284491, 6658042, 8388607, 10568983, 13316085, 16777216];
  const FIRSTIDX = 9;
  function sizeofint(size) { let n = 1, bits = 0; while (size >= n && bits < 32) { bits++; n *= 2; } return bits; }
  function sizeofints(sizes) {
    const bytes = [1];
    for (const s of sizes) {
      let tmp = 0, k = 0;
      for (; k < bytes.length; k++) { tmp = bytes[k] * s + tmp; bytes[k] = tmp % 256; tmp = Math.floor(tmp / 256); }
      while (tmp !== 0) { bytes[k++] = tmp % 256; tmp = Math.floor(tmp / 256); }
    }
    let n = 1, bits = 0;
    const top = bytes[bytes.length - 1];
    while (top >= n) { bits++; n *= 2; }
    return bits + (bytes.length - 1) * 8;
  }
  function receivebits(st, nbits) {
    const mask = nbits >= 32 ? 0xffffffff : 2 ** nbits - 1, b = st.buf;
    let cnt = st.cnt, lastbits = st.lastbits, lastbyte = st.lastbyte, n = 0;
    while (nbits >= 8) {
      lastbyte = ((lastbyte << 8) | b[cnt++]) & 0xffff;
      n = (n | ((lastbyte >> lastbits) << (nbits - 8))) >>> 0;
      nbits -= 8;
    }
    if (nbits > 0) {
      if (lastbits < nbits) { lastbits += 8; lastbyte = ((lastbyte << 8) | b[cnt++]) & 0xffff; }
      lastbits -= nbits;
      n = (n | ((lastbyte >> lastbits) & ((1 << nbits) - 1))) >>> 0;
    }
    st.cnt = cnt; st.lastbits = lastbits; st.lastbyte = lastbyte;
    return (n & mask) >>> 0;
  }
  function receiveints(st, nbits, sizes, out, o) {
    const bytes = [0, 0, 0, 0];
    let nb = 0;
    while (nbits > 8) { bytes[nb++] = receivebits(st, 8); nbits -= 8; }
    if (nbits > 0) bytes[nb++] = receivebits(st, nbits);
    for (let i = 2; i > 0; i--) {
      let v = 0;
      for (let j = nb - 1; j >= 0; j--) { v = v * 256 + bytes[j]; const q = Math.floor(v / sizes[i]); bytes[j] = q; v -= q * sizes[i]; }
      out[o + i] = v;
    }
    out[o] = (bytes[0] | (bytes[1] << 8) | (bytes[2] << 16) | ((bytes[3] || 0) << 24));
  }
  function xtcDecode(buf, n, minint, maxint, smallidx0, precision) {
    const sizeint = [0, 1, 2].map(k => maxint[k] - minint[k] + 1);
    let bitsize = 0, bitsizeint = null;
    if ((sizeint[0] | sizeint[1] | sizeint[2]) > 0xffffff || sizeint.some(s => s > 0xffffff)) bitsizeint = sizeint.map(sizeofint);
    else bitsize = sizeofints(sizeint);
    let smallidx = smallidx0;
    if (smallidx < FIRSTIDX || smallidx >= MAGICINTS.length) throw err('Bad XTC compression index.');
    let smaller = Math.floor(MAGICINTS[Math.max(FIRSTIDX, smallidx - 1)] / 2), smallnum = Math.floor(MAGICINTS[smallidx] / 2);
    let sizesmall = [MAGICINTS[smallidx], MAGICINTS[smallidx], MAGICINTS[smallidx]];
    const st = { buf, cnt: 0, lastbits: 0, lastbyte: 0 }, inv = 1 / precision;
    const out = new Float64Array(3 * n), cur = [0, 0, 0], prev = [0, 0, 0], tc = [0, 0, 0];
    let i = 0, run = 0, w = 0;
    while (i < n) {
      if (bitsizeint) { for (let k = 0; k < 3; k++) cur[k] = receivebits(st, bitsizeint[k]); }
      else receiveints(st, bitsize, sizeint, cur, 0);
      i++;
      for (let k = 0; k < 3; k++) { cur[k] += minint[k]; prev[k] = cur[k]; }
      const flag = receivebits(st, 1);
      let isSmaller = 0;
      if (flag === 1) { run = receivebits(st, 5); isSmaller = run % 3; run -= isSmaller; isSmaller--; }
      if (run > 0) {
        if (w + run + 3 > 3 * n) throw err('Corrupt XTC frame (run exceeds the atom count).');
        for (let k = 0; k < run; k += 3) {
          receiveints(st, smallidx, sizesmall, tc, 0);
          i++;
          for (let d = 0; d < 3; d++) tc[d] += prev[d] - smallnum;
          if (k === 0) {            // the first small atom was stored after the large one (water: O,H,H)
            for (let d = 0; d < 3; d++) { const t = tc[d]; tc[d] = prev[d]; prev[d] = t; }
            out[w++] = prev[0] * inv; out[w++] = prev[1] * inv; out[w++] = prev[2] * inv;
          } else { prev[0] = tc[0]; prev[1] = tc[1]; prev[2] = tc[2]; }
          out[w++] = tc[0] * inv; out[w++] = tc[1] * inv; out[w++] = tc[2] * inv;
        }
      } else { out[w++] = cur[0] * inv; out[w++] = cur[1] * inv; out[w++] = cur[2] * inv; }
      smallidx += isSmaller;
      if (isSmaller < 0) { smallnum = smaller; smaller = smallidx > FIRSTIDX ? Math.floor(MAGICINTS[smallidx - 1] / 2) : 0; }
      else if (isSmaller > 0) { smaller = smallnum; smallnum = Math.floor(MAGICINTS[smallidx] / 2); }
      if (smallidx < FIRSTIDX || smallidx >= MAGICINTS.length) throw err('Corrupt XTC frame (compression index out of range).');
      sizesmall = [MAGICINTS[smallidx], MAGICINTS[smallidx], MAGICINTS[smallidx]];
    }
    return out;
  }
  function readXtc(u8, fileName) {
    const dv = new DataView(u8.buffer, u8.byteOffset, u8.byteLength), N = u8.length;
    const ds = newDs('xtc', fileName, 'GROMACS');
    let p = 0, natoms = 0, truncated = false;
    while (p + 4 <= N) {
      try {
        if (dv.getInt32(p) !== 1995) { if (!ds.frames.length) throw err('Not an XTC file (bad magic number).'); truncated = true; break; }
        if (p + 56 > N) { truncated = true; break; }
        const nat = dv.getInt32(p + 4), step = dv.getInt32(p + 8), t = dv.getFloat32(p + 12);
        const b = []; for (let k = 0; k < 9; k++) b.push(dv.getFloat32(p + 16 + 4 * k) * 10);
        let q = p + 52;
        const lsize = dv.getInt32(q); q += 4;
        if (lsize !== nat || nat <= 0) throw err('Corrupt XTC frame header.');
        if (natoms && nat !== natoms) throw err('The atom count changes between XTC frames.');
        natoms = nat;
        let xyz;
        if (nat <= 9) {
          if (q + 12 * nat > N) { truncated = true; break; }
          xyz = new Float64Array(3 * nat);
          for (let k = 0; k < 3 * nat; k++) xyz[k] = dv.getFloat32(q + 4 * k) * 10;
          q += 12 * nat;
        } else {
          if (q + 36 > N) { truncated = true; break; }
          const prec = dv.getFloat32(q);
          const minint = [dv.getInt32(q + 4), dv.getInt32(q + 8), dv.getInt32(q + 12)], maxint = [dv.getInt32(q + 16), dv.getInt32(q + 20), dv.getInt32(q + 24)];
          const smallidx = dv.getInt32(q + 28), nbytes = dv.getInt32(q + 32);
          q += 36;
          if (nbytes < 0 || q + nbytes > N) { truncated = true; break; }
          const buf = new Uint8Array(nbytes + 8);
          buf.set(u8.subarray(q, q + nbytes));
          q += (nbytes + 3) & ~3;
          xyz = xtcDecode(buf, nat, minint, maxint, smallidx, prec);
          for (let k = 0; k < xyz.length; k++) xyz[k] *= 10;
        }
        const f = newFrame(null, xyz);
        f.time = t; f.step = step; f.label = 't = ' + (+t.toFixed(4)) + ' ps';
        const cell = [b.slice(0, 3), b.slice(3, 6), b.slice(6, 9)];
        if (goodCell(cell)) f.cell = cell;
        ds.frames.push(f);
        p = q;
      } catch (e) {
        if (!ds.frames.length) throw e;
        ds.warnings.push('Reading stopped at a damaged frame: ' + e.message);
        break;
      }
    }
    if (!ds.frames.length) throw err('The XTC file holds no complete frame.');
    if (truncated) ds.warnings.push('The file ends with an incomplete frame, which was skipped.');
    settleCells(ds, false);
    ds.meta.needsTopology = true;
    ds.meta.natoms = natoms;
    return ds;
  }

  // ================================================================== trajectories from Blob slices (big files)
  // DCD frames have a fixed size: only the selected frames are read. TRR / XTC are scanned frame by frame in large
  // windows; frames that are not selected are skipped without decoding. Options for all three:
  //   {stride, first, last (0-based frame indices, inclusive), budget (bytes of coordinates to keep, 512 MB),
  //    onProgress({bytes, total, frames}), signal}
  // Without a stride, one is chosen so that the kept coordinates fit the budget. The dataset carries
  // meta.totalFrames, meta.stride, meta.first, meta.last and meta.reload(opts) -> Promise<Dataset> (same source).
  const TRAJ_BUDGET = 512 * 1048576, WINDOW = 16 * 1048576;
  function trajReader(src) {
    const R = BM.bigfile && BM.bigfile.byteReader(src);
    if (!R) throw err('This trajectory source cannot be read in pieces.');
    let w0 = 0, w = new Uint8Array(0);
    return {
      size: R.size,
      /** bytes [off, off + len) (the window grows to `ahead` bytes when that is cheaper than many small reads) */
      async get(off, len, ahead) {
        if (off >= w0 && off + len <= w0 + w.byteLength) return w.subarray(off - w0, off - w0 + len);
        const want = Math.min(R.size - off, Math.max(len, ahead || 0));
        w = await R.read(off, want); w0 = off;
        if (w.byteLength < len) throw err('The file ended unexpectedly at byte ' + (off + w.byteLength) + '.');
        return w.subarray(0, len);
      }
    };
  }
  function trajSelection(total, perFrame, opts) {
    const first = Math.max(0, Math.floor(+opts.first || 0));
    const last = Math.min(total - 1, Number.isFinite(+opts.last) && opts.last !== null && opts.last !== undefined ? Math.floor(+opts.last) : total - 1);
    const count = Math.max(0, last - first + 1);
    let stride = Math.max(1, Math.floor(+opts.stride || 0));
    if (!opts.stride) {
      const max = Math.max(1, Math.floor((opts.budget || TRAJ_BUDGET) / Math.max(1, perFrame)));
      stride = Math.max(1, Math.ceil(count / max));
    }
    return { first, last, stride };
  }
  function finishTraj(ds, sel, total, natoms, src, fileName, opts, reader) {
    ds.meta.needsTopology = true;
    ds.meta.natoms = natoms;
    ds.meta.totalFrames = total;
    ds.meta.stride = sel.stride; ds.meta.first = sel.first; ds.meta.last = sel.last;
    Object.defineProperty(ds.meta, 'reload', { enumerable: false, configurable: true,
      value: o => reader(src, fileName, Object.assign({}, opts, { stride: null, first: null, last: null }, o || {})) });
    if (sel.stride > 1 || sel.first > 0 || sel.last < total - 1) {
      ds.warnings.push('Loaded ' + ds.frames.length + ' of ' + total + ' frames' + (sel.stride > 1 ? ' (every ' + sel.stride + 'th)' : '') +
        (sel.first > 0 || sel.last < total - 1 ? ', frames ' + (sel.first + 1) + '-' + (sel.last + 1) : '') + '.');
    }
    settleCells(ds, false);
    return ds;
  }
  const aborted = o => o && o.signal && o.signal.aborted;
  const abortErr = () => { const e = new Error('Reading was cancelled.'); e.name = 'AbortError'; return e; };

  async function readDcdSource(src, fileName, opts) {
    opts = opts || {};
    const rd = trajReader(src), N = rd.size;
    if (N < 12) throw err('Not a DCD file (too short).');
    let h = await rd.get(0, Math.min(N, 65536));
    let dv = new DataView(h.buffer, h.byteOffset, h.byteLength);
    const cord = p => h[p] === 67 && h[p + 1] === 79 && h[p + 2] === 82 && h[p + 3] === 68;
    let le, M;
    if (dv.getInt32(0, true) === 84 && cord(4)) { le = true; M = 4; }
    else if (dv.getInt32(0, false) === 84 && cord(4)) { le = false; M = 4; }
    else if (dv.getInt32(0, true) === 84 && dv.getInt32(4, true) === 0 && cord(8)) { le = true; M = 8; }
    else if (dv.getInt32(0, false) === 0 && dv.getInt32(4, false) === 84 && cord(8)) { le = false; M = 8; }
    else throw err('Not a DCD file (no CORD header record).');
    let pos = 0;
    const mark = p => (M === 4 ? dv.getInt32(p, le) : (le ? dv.getInt32(p, le) : dv.getInt32(p + 4, le)));
    const rec = async () => {
      if (pos + M > h.byteLength) return null;
      const len = mark(pos), start = pos + M, end = start + len;
      if (len < 0 || end + M > N) return null;
      if (end + M > h.byteLength) { h = await rd.get(0, end + M); dv = new DataView(h.buffer, h.byteOffset, h.byteLength); }
      if (mark(end) !== len) return null;
      pos = end + M;
      return { start, len };
    };
    const hd = await rec();
    if (!hd || hd.len < 84) throw err('The DCD header record is damaged.');
    const ic = []; for (let k = 0; k < 20; k++) ic.push(dv.getInt32(hd.start + 4 + 4 * k, le));
    const charmm = ic[19] !== 0, namnf = ic[8], hasCell = charmm && ic[10] !== 0, has4D = charmm && ic[11] !== 0;
    const tr = await rec();
    let title = '';
    if (tr && tr.len >= 4) {
      const nt = dv.getInt32(tr.start, le);
      for (let k = 0; k < nt && 4 + 80 * (k + 1) <= tr.len; k++) title += String.fromCharCode.apply(null, h.subarray(tr.start + 4 + 80 * k, tr.start + 4 + 80 * (k + 1))).replace(/\0/g, ' ').trim() + ' ';
    }
    const nr = await rec();
    if (!nr || nr.len < 4) throw err('The DCD atom-count record is missing.');
    const n = dv.getInt32(nr.start, le);
    if (!(n > 0)) throw err('The DCD file declares no atoms.');
    let free = null;
    if (namnf > 0) {
      const fr = await rec();
      if (!fr || fr.len !== 4 * (n - namnf)) throw err('The DCD free-atom list is damaged.');
      free = new Int32Array(n - namnf);
      for (let k = 0; k < free.length; k++) free[k] = dv.getInt32(fr.start + 4 * k, le) - 1;
    }
    const cellLen = hasCell ? 2 * M + 48 : 0, recLen = c => 2 * M + 4 * c;
    const size0 = cellLen + (has4D ? 4 : 3) * recLen(n), cnt = free ? free.length : n, sizeK = cellLen + (has4D ? 4 : 3) * recLen(cnt);
    const total = N - pos >= size0 ? 1 + Math.floor((N - pos - size0) / sizeK) : 0;
    if (!total) throw err('The DCD file holds no complete frame.');
    const sel = trajSelection(total, 24 * n, opts);
    const ds = newDs('dcd', fileName, charmm ? 'CHARMM/NAMD' : 'X-PLOR');
    const offset = k => pos + (k === 0 ? 0 : size0 + (k - 1) * sizeK);
    /** One frame from its bytes; null when the record markers do not fit (damaged file). */
    const decode = (u, k, base) => {
      const v = new DataView(u.buffer, u.byteOffset, u.byteLength);
      const mk = p => (M === 4 ? v.getInt32(p, le) : (le ? v.getInt32(p, le) : v.getInt32(p + 4, le)));
      let p = 0, cell = null;
      if (hasCell) {
        if (mk(p) !== 48 || mk(p + M + 48) !== 48) return null;
        const c = []; for (let q = 0; q < 6; q++) c.push(v.getFloat64(p + M + 8 * q, le));
        const cosines = [c[1], c[3], c[4]].every(x => x >= -1 && x <= 1), ang = x => (cosines ? Math.acos(x) / RAD : x);
        if (c[0] > 0 && c[2] > 0 && c[5] > 0) cell = boxCell(c[0], c[2], c[5], ang(c[4]), ang(c[3]), ang(c[1]));
        p += cellLen;
      }
      const m = k === 0 ? n : cnt;
      const xyz = base ? Float64Array.from(base) : new Float64Array(3 * n);
      for (let d = 0; d < 3; d++) {
        if (mk(p) !== 4 * m || mk(p + M + 4 * m) !== 4 * m) return null;
        const st = p + M;
        for (let q = 0; q < m; q++) xyz[3 * (m === n ? q : free[q]) + d] = v.getFloat32(st + 4 * q, le);
        p += recLen(m);
      }
      const f = newFrame(null, xyz);
      f.label = 'Frame ' + (k + 1);
      if (cell) f.cell = cell;
      return f;
    };
    let base = null, damaged = false;
    if (free) { base = decode(await rd.get(offset(0), size0), 0, null); if (!base) throw err('The first DCD frame is damaged.'); base = base.xyz; }
    const ahead = sel.stride * sizeK < WINDOW ? WINDOW : 0;
    for (let k = sel.first; k <= sel.last; k += sel.stride) {
      if (aborted(opts)) throw abortErr();
      const f = k === 0 && free ? newFrame(null, base) : decode(await rd.get(offset(k), k === 0 ? size0 : sizeK, ahead), k, base);
      if (!f) { damaged = true; break; }
      if (k === 0 && free) f.label = 'Frame 1';
      ds.frames.push(f);
      if (opts.onProgress && ds.frames.length % 64 === 0) opts.onProgress({ bytes: offset(k), total: N, frames: ds.frames.length });
    }
    if (!ds.frames.length) throw err('The DCD file holds no readable frame.');
    if (damaged) ds.warnings.push('A damaged frame stopped the reading.');
    if (pos + size0 + (total - 1) * sizeK < N) ds.warnings.push('The last frame is incomplete and was skipped.');
    if (ic[0] > 0 && ic[0] !== total) ds.warnings.push('The header announces ' + ic[0] + ' frames; ' + total + ' are in the file.');
    if (title.trim()) ds.title = title.trim().slice(0, 120);
    return finishTraj(ds, sel, total, n, src, fileName, opts, readDcdSource);
  }

  /** Sequential XDR trajectories: frame(u8 window, p) -> {end, nat, make()} or null (end of data / damaged). */
  async function scanXdr(src, fileName, opts, id, magic, headLen, frameAt) {
    opts = opts || {};
    const rd = trajReader(src), N = rd.size;
    const ds = newDs(id, fileName, 'GROMACS');
    let p = 0, k = 0, natoms = 0, sel = null, truncated = false;
    while (p + 4 <= N) {
      if (aborted(opts)) throw abortErr();
      const hl = Math.min(N - p, headLen);
      const hb = await rd.get(p, hl, WINDOW);
      const hv = new DataView(hb.buffer, hb.byteOffset, hb.byteLength);
      if (hv.getInt32(0) !== magic) { if (!k) throw err('Not a ' + id.toUpperCase() + ' file (bad magic number).'); truncated = true; break; }
      const info = frameAt(hv, hb.byteLength);
      if (!info || p + info.len > N) { if (!k) throw err('The first ' + id.toUpperCase() + ' frame is damaged or incomplete.'); truncated = true; break; }
      if (natoms && info.nat !== natoms) throw err('The atom count changes between ' + id.toUpperCase() + ' frames.');
      natoms = info.nat;
      if (!sel) sel = trajSelection(Math.max(1, Math.floor(N / info.len)), 24 * natoms, opts);
      if (k >= sel.first && k <= sel.last && (k - sel.first) % sel.stride === 0 && info.hasX) {
        const fb = await rd.get(p, info.len, WINDOW);
        try { ds.frames.push(info.make(new DataView(fb.buffer, fb.byteOffset, fb.byteLength), fb)); }
        catch (e) { if (!ds.frames.length) throw e; ds.warnings.push('Reading stopped at a damaged frame: ' + e.message); break; }
      }
      if (info.hasX || id === 'xtc') k++;
      p += info.len;
      if (opts.onProgress && k % 256 === 0) opts.onProgress({ bytes: p, total: N, frames: ds.frames.length });
    }
    if (!ds.frames.length) throw err('The ' + id.toUpperCase() + ' file holds no frame with coordinates' + (sel && sel.first ? ' in the chosen range' : '') + '.');
    if (truncated) ds.warnings.push('The file ends with an incomplete frame, which was skipped.');
    sel.last = Math.min(sel.last, k - 1);
    return finishTraj(ds, sel, k, natoms, src, fileName, opts, id === 'trr' ? readTrrSource : readXtcSource);
  }
  function readTrrSource(src, fileName, opts) {
    return scanXdr(src, fileName, opts, 'trr', 1993, 1024, (dv, L) => {
      if (L < 12) return null;
      const slen = dv.getInt32(8);
      let q = 12 + ((slen + 3) & ~3);
      if (slen < 0 || q + 52 > L) return null;
      const H = []; for (let k = 0; k < 13; k++) H.push(dv.getInt32(q + 4 * k));
      q += 52;
      const [irS, eS, boxS, virS, presS, topS, symS, xS, vS, fS, nat, step] = H;
      const ref = boxS ? boxS / 9 : xS ? xS / (3 * nat) : vS ? vS / (3 * nat) : fS ? fS / (3 * nat) : 4;
      const R = ref === 8 ? 8 : 4;
      if (H.slice(0, 11).some(v => v < 0) || (xS && xS !== 3 * nat * R) || (boxS && boxS !== 9 * R)) return null;
      const len = q + 2 * R + irS + eS + boxS + virS + presS + topS + symS + xS + vS + fS;
      return { len, nat, hasX: xS > 0, make: (v) => {
        const rdv = o => (R === 8 ? v.getFloat64(o) : v.getFloat32(o));
        const t = rdv(q);
        let o = q + 2 * R + irS + eS, cell = null;
        if (boxS) { const b = []; for (let k = 0; k < 9; k++) b.push(rdv(o + R * k) * 10); cell = [b.slice(0, 3), b.slice(3, 6), b.slice(6, 9)]; if (!goodCell(cell)) cell = null; }
        o += boxS + virS + presS + topS + symS;
        const xyz = new Float64Array(3 * nat);
        for (let k = 0; k < 3 * nat; k++) xyz[k] = rdv(o + R * k) * 10;
        const f = newFrame(null, xyz);
        f.time = t; f.step = step; f.label = 't = ' + (+t.toFixed(4)) + ' ps';
        if (cell) f.cell = cell;
        return f;
      } };
    });
  }
  function readXtcSource(src, fileName, opts) {
    return scanXdr(src, fileName, opts, 'xtc', 1995, 92, (dv, L) => {
      if (L < 56) return null;
      const nat = dv.getInt32(4), step = dv.getInt32(8), t = dv.getFloat32(12);
      if (dv.getInt32(52) !== nat || nat <= 0) return null;
      let len, nbytes = 0;
      if (nat <= 9) len = 56 + 12 * nat;
      else { if (L < 92) return null; nbytes = dv.getInt32(88); if (nbytes < 0) return null; len = 92 + ((nbytes + 3) & ~3); }
      return { len, nat, hasX: true, make: (v, u) => {
        const b = []; for (let k = 0; k < 9; k++) b.push(v.getFloat32(16 + 4 * k) * 10);
        let xyz;
        if (nat <= 9) { xyz = new Float64Array(3 * nat); for (let k = 0; k < 3 * nat; k++) xyz[k] = v.getFloat32(56 + 4 * k) * 10; }
        else {
          const prec = v.getFloat32(56);
          const minint = [v.getInt32(60), v.getInt32(64), v.getInt32(68)], maxint = [v.getInt32(72), v.getInt32(76), v.getInt32(80)];
          const buf = new Uint8Array(nbytes + 8);
          buf.set(u.subarray(92, 92 + nbytes));
          xyz = xtcDecode(buf, nat, minint, maxint, v.getInt32(84), prec);
          for (let k = 0; k < xyz.length; k++) xyz[k] *= 10;
        }
        const f = newFrame(null, xyz);
        f.time = t; f.step = step; f.label = 't = ' + (+t.toFixed(4)) + ' ps';
        const cell = [b.slice(0, 3), b.slice(3, 6), b.slice(6, 9)];
        if (goodCell(cell)) f.cell = cell;
        return f;
      } };
    });
  }

  // ================================================================== registration
  const REG = o => BM.registerFormat(o);
  REG({ id: 'cif', name: 'CIF', extensions: ['cif'], kind: 'structure', sniff: sniffCif, parse: readCif });
  REG({ id: 'mmcif', name: 'mmCIF / PDBx', extensions: ['cif', 'mmcif', 'mcif'], kind: 'structure', sniff: sniffMmcif, parse: readCif });
  REG({ id: 'mol2', name: 'Tripos Mol2', extensions: ['mol2', 'ml2', 'sy2'], kind: 'structure',
    sniff: h => (/^\s*@<TRIPOS>MOLECULE/m.test(h) ? 0.95 : /^\s*@<TRIPOS>ATOM/m.test(h) ? 0.9 : 0), parse: readMol2 });
  REG({ id: 'cml', name: 'Chemical Markup Language', extensions: ['cml'], kind: 'structure',
    sniff: h => (/xml-cml\.org/.test(h) ? 0.95 : /<([\w.-]+:)?molecule[\s>]/.test(h) && /<([\w.-]+:)?atom(Array)?[\s>]/.test(h) ? 0.9 : 0), parse: readCML });
  REG({ id: 'cjson', name: 'Avogadro CJSON', extensions: ['cjson'], kind: 'structure',
    sniff: h => (!/^\s*\{/.test(h) ? 0 : /"chemical[ _]?json"/i.test(h) ? 0.95 : /"atoms"\s*:\s*\{/.test(h) && /"elements"\s*:/.test(h) ? 0.6 : 0), parse: readCJSON });
  REG({ id: 'poscar', name: 'VASP POSCAR/CONTCAR', extensions: ['vasp', 'poscar', 'contcar'], kind: 'structure', sniff: sniffPoscar, parse: readPoscar });
  REG({ id: 'xdatcar', name: 'VASP XDATCAR', extensions: ['xdatcar'], kind: 'trajectory',
    sniff: (h, n) => (/xdatcar/.test(n || '') || /^\s*(Direct|Cartesian) configuration=/m.test(h) ? 0.95 : 0), parse: readXdatcar });
  REG({ id: 'outcar', name: 'VASP OUTCAR', extensions: ['outcar'], kind: 'output',
    sniff: (h, n) => (/^\s*vasp\.\d/.test(h) || (/outcar/.test(n || '') && /POSITION\s+TOTAL-FORCE|TITEL/.test(h)) ? 0.95 : 0), parse: readOutcar });
  REG({ id: 'qe-in', name: 'Quantum ESPRESSO pw.x input', extensions: ['pwi', 'in'], kind: 'structure',
    sniff: h => (/^\s*&(control|system)\b/im.test(h) ? (/ATOMIC_POSITIONS/i.test(h) ? 0.95 : 0.6) : 0), parse: readQeIn });
  REG({ id: 'qe-out', name: 'Quantum ESPRESSO output', extensions: ['pwo', 'out'], kind: 'output',
    sniff: h => (/Program PWSCF/.test(h) ? 0.95 : 0), parse: readQeOut });
  REG({ id: 'castep-cell', name: 'CASTEP cell', extensions: ['cell'], kind: 'structure',
    sniff: h => (/^\s*%block\s+(lattice_cart|lattice_abc|positions_frac|positions_abs)\b/im.test(h) ? 0.9 : 0), parse: readCastepCell });
  REG({ id: 'xsf', name: 'XCrySDen XSF/AXSF', extensions: ['xsf', 'axsf'], kind: 'structure',
    sniff: h => (/^\s*(PRIMCOORD|PRIMVEC|BEGIN_BLOCK_DATAGRID_?3D|ANIMSTEPS)\b/m.test(h) ? 0.9 : /^\s*ATOMS\s*$/m.test(h) ? 0.4 : 0), parse: readXsf });
  REG({ id: 'cp2k-in', name: 'CP2K input', extensions: ['inp', 'cp2k'], kind: 'structure',
    sniff: h => (/^\s*&SUBSYS\b/im.test(h) ? (/^\s*&FORCE_EVAL\b/im.test(h) ? 0.9 : 0.6) : 0), parse: readCp2kIn });
  REG({ id: 'prmtop', name: 'Amber topology (prmtop)', extensions: ['prmtop', 'parm7'], kind: 'topology',
    sniff: h => (/^%VERSION/.test(h) && /%FLAG POINTERS/.test(h) ? 0.95 : /^%FLAG\s+(TITLE|POINTERS)/m.test(h) ? 0.8 : 0), parse: readPrmtop });
  REG({ id: 'inpcrd', name: 'Amber coordinates (inpcrd/rst7)', extensions: ['inpcrd', 'rst7', 'restrt'], kind: 'trajectory', sniff: sniffInpcrd, parse: readInpcrd });
  REG({ id: 'psf', name: 'CHARMM/NAMD PSF', extensions: ['psf'], kind: 'topology',
    sniff: h => (/^\s*PSF\b/.test(h) && /!NATOM/.test(h) ? 0.95 : 0), parse: readPsf });
  REG({ id: 'dcd', name: 'CHARMM/NAMD DCD', extensions: ['dcd'], kind: 'trajectory', binary: true, parseBinary: readDcd, parseSource: readDcdSource });
  REG({ id: 'lammps-data', name: 'LAMMPS data', extensions: ['data', 'lmp', 'lammps'], kind: 'structure',
    sniff: h => (/^\s*[-+\d.eE]+\s+[-+\d.eE]+\s+xlo\s+xhi\b/m.test(h) ? (/^\s*\d+\s+atoms\b/m.test(h) ? 0.95 : 0.7) : 0), parse: readLammpsData });
  REG({ id: 'lammps-dump', name: 'LAMMPS dump', extensions: ['lammpstrj', 'dump'], kind: 'trajectory',
    sniff: h => (/^ITEM: TIMESTEP/m.test(h) ? 0.95 : 0), parse: readLammpsDump });
  REG({ id: 'trr', name: 'GROMACS TRR', extensions: ['trr'], kind: 'trajectory', binary: true, parseBinary: readTrr, parseSource: readTrrSource });
  REG({ id: 'xtc', name: 'GROMACS XTC', extensions: ['xtc'], kind: 'trajectory', binary: true, parseBinary: readXtc, parseSource: readXtcSource });

  BM.structReaders = { cifBlocks, parseSymop, cellFromParams, zFromMass, xtcDecode,
    parsers: { cif: readCif, mol2: readMol2, cml: readCML, cjson: readCJSON, poscar: readPoscar, xdatcar: readXdatcar, outcar: readOutcar,
      qeIn: readQeIn, qeOut: readQeOut, castepCell: readCastepCell, xsf: readXsf, cp2kIn: readCp2kIn, prmtop: readPrmtop,
      inpcrd: readInpcrd, psf: readPsf, dcd: readDcd, lammpsData: readLammpsData, lammpsDump: readLammpsDump, trr: readTrr, xtc: readXtc },
    sources: { dcd: readDcdSource, trr: readTrrSource, xtc: readXtcSource } };
})(typeof window !== 'undefined' ? (window.BM = window.BM || {}) : (globalThis.BM = globalThis.BM || {}));
