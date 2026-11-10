/* BlueMolecule — molecule building and editing interface.
 *
 *   Tool palette and tools (navigate, select, draw, template, bond-centric, manipulate), the context options
 *   bar, the Build tab, force-field clean-up runs (BM.ff), symmetry (BM.sym), the coordinate editor and the
 *   element / template pickers.  Small chemistry helpers shared with app.js live in BM.buildUI.kit
 *   (bond-order perception for display, benzenoid rings, suggested total charge, hydrogen bonds).
 *
 * Original code.  Loaded after the engines and before app.js; app.js calls BM.buildUI.create(api).
 */
(function (BM) {
  'use strict';
  const el = BM.el, geom = BM.geom;
  const LIMIT = { perceive: 2000, auto: 1500, clean: 3000, sym: 1000, editor: 20000 };
  const QUICK = [6, 1, 7, 8, 16, 15, 9, 17, 35, 53, 78];
  const TOOLS = [
    { id: 'navigate', name: 'Navigate', key: 'N', icon: 'i-t-nav', tip: 'rotate, pan and zoom; click atoms to measure' },
    { id: 'select', name: 'Select', key: 'S', icon: 'i-t-select', tip: 'box-select and pick atoms' },
    { id: 'draw', name: 'Draw', key: 'D', icon: 'i-t-draw', tip: 'add atoms and bonds' },
    { id: 'template', name: 'Template', key: 'T', icon: 'i-t-template', tip: 'molecules, groups, metal centres, ligands' },
    { id: 'bond', name: 'Bond-centric', key: 'K', icon: 'i-t-bond', tip: 'rotate about or stretch a bond' },
    { id: 'manipulate', name: 'Manipulate', key: 'M', icon: 'i-t-move', tip: 'move and rotate fragments' }
  ];
  const TPL_MODES = [['groups', 'Groups', 'Attach a group to an atom'], ['molecules', 'Molecules', 'Insert a molecule'], ['metal', 'Metal', 'Build a metal centre'], ['ligands', 'Ligands', 'Replace ligands at selected sites']];
  const ORDER_NAME = { 1: 'single', 1.5: 'aromatic', 2: 'double', 3: 'triple' };
  const SYM_TOLS = [[0.02, 'tight · 0.02 Å'], [0.1, 'normal · 0.1 Å'], [0.3, 'loose · 0.3 Å']];

  // ------------------------------------------------------------------------------------------- vectors
  const vsub = (a, b) => [a[0] - b[0], a[1] - b[1], a[2] - b[2]];
  const vadd = (a, b) => [a[0] + b[0], a[1] + b[1], a[2] + b[2]];
  const vmul = (a, s) => [a[0] * s, a[1] * s, a[2] * s];
  const vdot = (a, b) => a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
  const vcross = (a, b) => [a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]];
  const vlen = a => Math.sqrt(vdot(a, a));
  const vunit = a => { const l = vlen(a); return l > 1e-12 ? vmul(a, 1 / l) : [0, 0, 0]; };
  const P = (xyz, i) => [xyz[3 * i], xyz[3 * i + 1], xyz[3 * i + 2]];
  const uniq = list => [...new Set(list)];
  const isMetal = z => BM.chem.isMetal(z);

  // ============================================================================ kit (no DOM)
  function pairsToList(pairs, orders) {
    const out = [];
    for (let k = 0; k + 1 < pairs.length; k += 2) out.push([pairs[k], pairs[k + 1], orders ? +orders[k >> 1] || 1 : 1]);
    return out;
  }
  /** Kekulé six-rings whose atoms each carry one double bond to an atom of such a ring (benzene, pyridine,
   *  fused arenes) get order 1.5, so the clean-up force field and the renderer treat them as aromatic. */
  function aromatize(Z, bonds) {
    const n = Z.length, out = bonds.map(b => [b[0], b[1], b[2] === undefined ? 1 : +b[2]]);
    if (n < 6 || out.length < 6) return out;
    const adj = Array.from({ length: n }, () => []), dbl = new Int32Array(n), mate = new Int32Array(n).fill(-1);
    out.forEach((b, k) => {
      const [i, j, o] = b;
      if (i === j || i < 0 || j < 0 || i >= n || j >= n) return;
      adj[i].push([j, k]); adj[j].push([i, k]);
      if (o === 2) { dbl[i]++; dbl[j]++; mate[i] = j; mate[j] = i; }
    });
    const RING = new Set([5, 6, 7, 8, 15, 16]), cand = new Uint8Array(n);
    for (let i = 0; i < n; i++) cand[i] = RING.has(Z[i]) && dbl[i] === 1 && adj[i].length <= 3 ? 1 : 0;
    const rings = new Map(), seen = new Int32Array(n).fill(-1), dist = new Int32Array(n), prev = new Int32Array(n);
    out.forEach((b, k) => {
      const a = b[0], c = b[1];
      if (a === c || !cand[a] || !cand[c]) return;
      const queue = [a];
      seen[a] = k; dist[a] = 0; prev[a] = -1;
      let found = false;
      for (let h = 0; h < queue.length && !found; h++) {
        const u = queue[h];
        if (dist[u] >= 5) continue;
        for (const [v, kk] of adj[u]) {
          if (kk === k || !cand[v] || seen[v] === k) continue;
          seen[v] = k; dist[v] = dist[u] + 1; prev[v] = u;
          if (v === c) { found = true; break; }
          queue.push(v);
        }
      }
      if (!found || dist[c] !== 5) return;
      const ring = [];
      for (let x = c; x !== -1; x = prev[x]) ring.push(x);
      const key = ring.slice().sort((p, q) => p - q).join(',');
      if (!rings.has(key)) rings.set(key, ring);
    });
    if (!rings.size) return out;
    const inRing = new Uint8Array(n);
    for (const r of rings.values()) for (const a of r) inRing[a] = 1;
    const index = new Map();
    out.forEach((b, k) => index.set(Math.min(b[0], b[1]) * n + Math.max(b[0], b[1]), k));
    for (const r of rings.values()) {
      if (!r.every(a => mate[a] >= 0 && inRing[mate[a]])) continue;
      for (let t = 0; t < 6; t++) {
        const a = r[t], c = r[(t + 1) % 6], k = index.get(Math.min(a, c) * n + Math.max(a, c));
        if (k !== undefined) out[k][2] = 1.5;
      }
    }
    return out;
  }
  /** Bond orders and formal charges for a fixed connectivity (pairs: flat Uint32Array or [[i, j, …]]). */
  function perceive(Z, xyz, pairs, fc) {
    const list = ArrayBuffer.isView(pairs) ? pairsToList(pairs) : pairs.map(b => [b[0], b[1], 1]);
    const n = Z.length;
    if (!list.length) return { bonds: [], charges: fc && fc.length === n ? Array.from(fc, q => q || 0) : new Array(n).fill(0) };
    const mol = { Z: Array.from(Z), xyz, bonds: list };
    if (fc && fc.length === n && Array.from(fc).some(q => q)) mol.formalCharges = Array.from(fc, q => q || 0);
    const r = BM.chem.perceiveBondOrders(mol, { withCharges: true });
    return { bonds: aromatize(mol.Z, r.bonds), charges: Array.from(r.formalCharges, q => q || 0) };
  }

  // Common oxidation states, used only to suggest a total charge for metal complexes.
  const COMMON_OX = (() => {
    const o = {};
    [3, 11, 19, 37, 55, 47].forEach(z => { o[z] = 1; });
    [4, 12, 20, 38, 56, 25, 26, 28, 29, 30, 44, 46, 48, 78, 80, 82].forEach(z => { o[z] = 2; });
    [13, 21, 24, 27, 31, 39, 45, 49, 77, 83].forEach(z => { o[z] = 3; });
    for (let z = 57; z <= 71; z++) o[z] = 3;
    [22, 40, 72].forEach(z => { o[z] = 4; });
    return o;
  })();
  const VALENCES = { 1: [1], 5: [3], 6: [4], 7: [3], 8: [2], 9: [1], 14: [4], 15: [3, 5], 16: [2, 4, 6], 17: [1, 3, 5, 7],
    32: [4], 33: [3, 5], 34: [2, 4, 6], 35: [1, 3, 5, 7], 52: [2, 4, 6], 53: [1, 3, 5, 7] };
  /** Ionic-model charge of a ligand atom from its bond-order sum without metal bonds (Cl− 0 bonds, CO 0 net, …). */
  function valenceCharge(z, s) {
    const v = VALENCES[z];
    if (!v) return 0;
    s = Math.round(s);
    if (z === 6 && s === 2) return 0; // carbene donor
    if (s < v[0]) return s - v[0];
    if (v.indexOf(s) >= 0) return 0;
    const max = v[v.length - 1];
    if (s > max) return z === 5 ? max - s : s - max;
    return 0;
  }
  const ROMAN = ['0', 'I', 'II', 'III', 'IV', 'V', 'VI', 'VII', 'VIII'];
  /** Suggested total charge of a structure {Z, bonds, formalCharges}: the sum of formal charges, plus — for
   *  metal complexes — ligand charges (ionic model) and each metal's oxidation state (its formal charge if set,
   *  otherwise a common value).  determinable is false if a metal has no default. */
  function suggestCharge(f) {
    const Z = f.Z, n = Z.length, fc = f.formalCharges || [], bonds = f.bonds || [];
    let any = false;
    for (let i = 0; i < n && !any; i++) if (isMetal(Z[i])) any = true;
    if (!any) { let s = 0; for (let i = 0; i < n; i++) s += fc[i] || 0; return { total: s, determinable: true, metals: [], ligand: 0 }; }
    const adj = Array.from({ length: n }, () => []);
    for (const b of bonds) { if (b[0] < n && b[1] < n) { adj[b[0]].push([b[1], b[2] === undefined ? 1 : +b[2]]); adj[b[1]].push([b[0], b[2] === undefined ? 1 : +b[2]]); } }
    const inCx = new Uint8Array(n), stack = [];
    for (let i = 0; i < n; i++) if (isMetal(Z[i]) && !inCx[i]) { inCx[i] = 1; stack.push(i); while (stack.length) { const u = stack.pop(); for (const [v] of adj[u]) if (!inCx[v]) { inCx[v] = 1; stack.push(v); } } }
    let total = 0, ligand = 0, ok = true;
    const metals = [];
    for (let i = 0; i < n; i++) {
      if (!(Z[i] > 0)) continue;
      if (isMetal(Z[i])) {
        const ox = fc[i] ? fc[i] : COMMON_OX[Z[i]];
        if (ox === undefined) { ok = false; metals.push({ z: Z[i], ox: null, assumed: true }); } else { total += ox; metals.push({ z: Z[i], ox, assumed: !fc[i] }); }
        continue;
      }
      if (!inCx[i]) { total += fc[i] || 0; continue; }
      let q = fc[i] || 0;
      if (!q) { let s = 0; for (const [j, o] of adj[i]) if (!isMetal(Z[j])) s += o; q = valenceCharge(Z[i], s); }
      total += q; ligand += q;
    }
    return { total, determinable: ok, metals, ligand };
  }
  function chargeText(q) { return q === 0 ? '0' : (q > 0 ? '+' : '−') + Math.abs(q); }
  function oxText(metals) {
    const seen = new Set(), parts = [];
    for (const m of metals) {
      const key = m.z + ':' + m.ox;
      if (!m.assumed || seen.has(key)) continue;
      seen.add(key);
      parts.push(m.ox === null ? el.symbol(m.z) + ' (oxidation state unknown)' : el.symbol(m.z) + '(' + (m.ox < 0 ? '−' + ROMAN[-m.ox] : ROMAN[m.ox] || m.ox) + ')');
    }
    return parts.join(', ');
  }

  /** D–H…A hydrogen bonds (D, A in N/O/F, H…A ≤ 2.5 Å, D–H…A ≥ 120°) for bonded pairs (flat array). */
  function hydrogenBonds(Z, xyz, pairs) {
    const n = Z.length, NOF = z => z === 7 || z === 8 || z === 9;
    const donor = new Int32Array(n).fill(-1);
    let nd = 0;
    for (let k = 0; k + 1 < pairs.length; k += 2) {
      const i = pairs[k], j = pairs[k + 1];
      if (i >= n || j >= n) continue;
      if (Z[i] === 1 && NOF(Z[j])) { donor[i] = j; nd++; } else if (Z[j] === 1 && NOF(Z[i])) { donor[j] = i; nd++; }
    }
    if (!nd) return [];
    const CUT = 2.5, grid = new Map(), cellOf = v => Math.floor(v / CUT);
    const key = (a, b, c) => ((a + 4096) * 8192 + (b + 4096)) * 8192 + (c + 4096);
    for (let i = 0; i < n; i++) {
      if (!NOF(Z[i])) continue;
      const k = key(cellOf(xyz[3 * i]), cellOf(xyz[3 * i + 1]), cellOf(xyz[3 * i + 2]));
      const l = grid.get(k);
      if (l) l.push(i); else grid.set(k, [i]);
    }
    const out = [], cosMax = -0.5; // 120°
    for (let h = 0; h < n; h++) {
      const D = donor[h];
      if (D < 0) continue;
      const hx = xyz[3 * h], hy = xyz[3 * h + 1], hz = xyz[3 * h + 2];
      const dx = xyz[3 * D] - hx, dy = xyz[3 * D + 1] - hy, dz = xyz[3 * D + 2] - hz, dl = Math.hypot(dx, dy, dz) || 1;
      const cx = cellOf(hx), cy = cellOf(hy), cz = cellOf(hz);
      for (let a = -1; a <= 1; a++) for (let b = -1; b <= 1; b++) for (let c = -1; c <= 1; c++) {
        const l = grid.get(key(cx + a, cy + b, cz + c));
        if (!l) continue;
        for (const A of l) {
          if (A === D) continue;
          const ax = xyz[3 * A] - hx, ay = xyz[3 * A + 1] - hy, az = xyz[3 * A + 2] - hz, r = Math.hypot(ax, ay, az);
          if (r > CUT || r < 1.1) continue;
          if ((dx * ax + dy * ay + dz * az) / (dl * r) <= cosMax) out.push([h, A]);
        }
      }
    }
    return out;
  }

  /** "Symbol x y z" lines with 8 decimals. */
  function coordText(Z, xyz) {
    const L = [];
    const f8 = v => ' ' + (Math.abs(v) < 5e-9 ? 0 : v).toFixed(8).padStart(12);
    for (let i = 0; i < Z.length; i++) L.push((Z[i] > 0 ? el.symbol(Z[i]) : 'X').padEnd(2) + f8(xyz[3 * i]) + f8(xyz[3 * i + 1]) + f8(xyz[3 * i + 2]));
    return L.join('\n');
  }
  const NUM = /^[-+]?(\d+\.?\d*|\.\d+)([eEdD][-+]?\d+)?$/;
  const num = s => Number(String(s).replace(/[dD]/, 'e'));
  function atomLine(s) {
    const t = s.trim().split(/[\s,]+/);
    if (t.length < 4 || !NUM.test(t[1]) || !NUM.test(t[2]) || !NUM.test(t[3])) return null;
    const z = el.toZ(t[0].replace(/[:(].*$/, ''));
    if (z < 0) return null;
    return { z, p: [num(t[1]), num(t[2]), num(t[3])] };
  }
  /** Coordinates from the editor: "Symbol x y z" lines, an ORCA "* xyz" block, a whole XYZ file, a Gaussian
   *  input or any other format BM.parseFile recognises.  Returns {Z, xyz, source}. */
  function parseCoords(text) {
    const src = String(text || '').replace(/\r\n?/g, '\n');
    const lines = src.split('\n').map(s => s.trim()).filter(s => s && !/^(#|!|\/\/)\s/.test(s) && !/^#\s*$/.test(s));
    if (!lines.length) throw new Error('The coordinate box is empty.');
    const plain = lines.map(s => (s.split(/[\s,]+/).length === 4 ? atomLine(s) : null));
    if (plain.every(Boolean)) return { Z: plain.map(a => a.z), xyz: Float64Array.from([].concat(...plain.map(a => a.p))), source: 'coordinates' };
    const m = /^[ \t]*\*[ \t]*xyz[ \t]+[-+]?\d+[ \t]+\d+[ \t]*$/im.exec(src);
    if (m) {
      const rest = src.slice(m.index + m[0].length).split('\n'), atoms = [];
      for (const s of rest) {
        if (/^\s*\*\s*$/.test(s)) break;
        if (!s.trim()) continue;
        const a = atomLine(s);
        if (!a) throw new Error('Cannot read this line of the ORCA geometry block: "' + s.trim() + '".');
        atoms.push(a);
      }
      if (!atoms.length) throw new Error('The ORCA "* xyz" block has no atoms.');
      return { Z: atoms.map(a => a.z), xyz: Float64Array.from([].concat(...atoms.map(a => a.p))), source: 'ORCA input' };
    }
    const first = src.split('\n').find(s => s.trim()) || '';
    let ds, source;
    if (/^\s*\d+\s*$/.test(first)) { ds = BM.parseFile('pasted.xyz', src); source = 'XYZ'; }
    else if (/^\s*#/m.test(src)) { ds = BM.parseFile('pasted.gjf', src); source = 'Gaussian input'; }
    else {
      try { ds = BM.parseFile('pasted', src); } catch (e) {
        throw new Error('Could not read the coordinates: use one "Symbol x y z" line per atom, or paste a whole XYZ file, Gaussian input or ORCA input.');
      }
      source = 'pasted ' + ds.format + ' text';
    }
    const fr = ds.frames[ds.frames.length - 1];
    if (!fr || !fr.Z.length) throw new Error('No atoms found in the pasted text.');
    if (ds.frames.length > 1) source += ' (last of ' + ds.frames.length + ' geometries)';
    return { Z: Array.from(fr.Z), xyz: Float64Array.from(fr.xyz), source };
  }

  // Periodic-table cell [row, column]; rows 9 and 10 are the f block.
  function ptCell(z) {
    if (z === 1) return [1, 1];
    if (z === 2) return [1, 18];
    if (z <= 10) return [2, z <= 4 ? z - 2 : z + 8];
    if (z <= 18) return [3, z <= 12 ? z - 10 : z];
    if (z <= 36) return [4, z - 18];
    if (z <= 54) return [5, z - 36];
    if (z <= 57) return [6, z - 54];
    if (z <= 71) return [9, z - 54];
    if (z <= 86) return [6, z - 68];
    if (z <= 89) return [7, z - 86];
    if (z <= 103) return [10, z - 86];
    return [7, z - 100];
  }

  const kit = { aromatize, perceive, suggestCharge, chargeText, oxText, hydrogenBonds, coordText, parseCoords, pairsToList, LIMIT };

  // ============================================================================ user interface
  function create(app) {
    const S = app.state, B = S.build, esc = app.esc;
    const $ = (s, r) => (r || document).querySelector(s), $$ = (s, r) => Array.from((r || document).querySelectorAll(s));
    const R = () => app.R();
    const doc = () => S.active;
    const cur = () => (S.active ? app.curFrame(S.active) : null);
    let drag = null, relax = null, chip = null, popup = null;

    // ---------------------------------------------------------------- small helpers
    const label = (f, i) => (f.Z[i] > 0 ? el.symbol(f.Z[i]) : 'X') + (i + 1);
    const orderName = o => ORDER_NAME[o] || String(o);
    const rgb = c => 'rgb(' + c.map(v => Math.round(v * 255)).join(',') + ')';
    function localXY(e) { const r = R().canvas.getBoundingClientRect(); return [e.clientX - r.left, e.clientY - r.top]; }
    function adjacency(d, fi) {
      const f = d.data.frames[fi], n = f.Z.length, pairs = app.bondsFor(d, fi), adj = Array.from({ length: n }, () => []);
      for (let k = 0; k + 1 < pairs.length; k += 2) { const i = pairs[k], j = pairs[k + 1]; if (i < n && j < n && i !== j) { adj[i].push(j); adj[j].push(i); } }
      return adj;
    }
    const curAdj = () => adjacency(doc(), doc().frameIndex);
    function neighbours(i) { return curAdj()[i] || []; }
    function bondAtoms(k) { const b = R().bonds; return k >= 0 && 2 * k + 1 < b.length ? [b[2 * k], b[2 * k + 1]] : null; }
    function bondOrder(a, b) {
      const d = doc(), t = app.topoFor(d, d.frameIndex), pairs = t.pairs, orders = app.ordersFor(d, d.frameIndex, true);
      for (let k = 0; k + 1 < pairs.length; k += 2) if ((pairs[k] === a && pairs[k + 1] === b) || (pairs[k] === b && pairs[k + 1] === a)) return orders ? +orders[k >> 1] || 1 : 1;
      return 0;
    }
    function fragmentOf(i) {
      const adj = curAdj(), seen = new Set([i]), stack = [i];
      while (stack.length) { const u = stack.pop(); for (const v of adj[u]) if (!seen.has(v)) { seen.add(v); stack.push(v); } }
      return [...seen].sort((a, b) => a - b);
    }
    function centroid(xyz, list) {
      const c = [0, 0, 0];
      for (const i of list) { c[0] += xyz[3 * i]; c[1] += xyz[3 * i + 1]; c[2] += xyz[3 * i + 2]; }
      return list.length ? vmul(c, 1 / list.length) : c;
    }
    function composeMap(a, b, n0) {
      if (!a && !b) return null;
      const out = new Int32Array(n0);
      for (let i = 0; i < n0; i++) { const j = a ? a[i] : i; out[i] = j < 0 ? -1 : b ? (j < b.length ? b[j] : -1) : j; }
      return out;
    }
    function status(msg) { app.status(msg); }
    function toast(msg, err) { app.toast(msg, err); }
    function screenRight() { const m = BM.quat.mat(R().q); return [m[0], m[1], m[2]]; }
    function molRadius(xyz, n) {
      const c = centroid(xyz, Array.from({ length: n }, (_, i) => i));
      let r = 0;
      for (let i = 0; i < n; i++) r = Math.max(r, vlen(vsub(P(xyz, i), c)));
      return { c, r };
    }
    /** Where to put a new fragment of radius rNew next to the current structure (origin if empty). */
    function besidePoint(rNew) {
      const f = cur();
      if (!f || !f.Z.length) return [0, 0, 0];
      const { c, r } = molRadius(f.xyz, f.Z.length);
      return vadd(c, vmul(screenRight(), r + rNew + 2.6));
    }
    function setSelection(list) {
      const d = doc(); if (!d) return;
      const n = cur().Z.length;
      d.selection = uniq(list.filter(i => i >= 0 && i < n));
      R().setSelection(d.selection);
      app.updateMeasure();
    }
    function combineSel(base, list, mode) {
      if (mode === 'replace') return list.slice();
      if (mode === 'add') return uniq(base.concat(list));
      if (mode === 'remove') { const s = new Set(list); return base.filter(i => !s.has(i)); }
      const s = new Set(base);
      for (const i of list) { if (s.has(i)) s.delete(i); else s.add(i); }
      return [...s];
    }
    const selMode = e => (e.shiftKey ? 'add' : e.ctrlKey || e.metaKey ? 'remove' : 'replace');
    function ensureDoc() { const d = doc(); return d && !d.pending ? d : newMolecule({ keepTool: true }); }

    // ---------------------------------------------------------------- edits
    /** Run fn(mol, f, d) → {mol, map?, select?} | null as one undoable edit; msg may be a function of the result. */
    function edit(fn, msg, opts) {
      opts = opts || {};
      const d = doc(); if (!d) return false;
      const wasEmpty = !cur().Z.length;
      let res = null;
      const ok = app.doEdit((f, dd) => {
        res = fn(app.frameMol(f), f, dd);
        if (!res) return false;
        app.commitMol(f, dd, res.mol, res.map || null);
        if (res.select) dd.selection = res.select.filter(i => i >= 0 && i < f.Z.length);
        return true;
      }, typeof msg === 'function' ? () => msg(res) : msg, { fit: wasEmpty || opts.fit, auto: opts.auto });
      return ok;
    }
    function adjusted(m, map, affected, n0) {
      if (!B.adjustH) return { mol: m, map };
      const atoms = uniq(affected.filter(i => i >= 0 && i < m.Z.length));
      if (!atoms.length) return { mol: m, map };
      const r = BM.chem.adjustHydrogens(m, { atoms });
      return { mol: r.mol, map: composeMap(map, r.map, n0) };
    }
    function hParent(i) {
      const f = cur();
      if (f.Z[i] !== 1) return -1;
      const nb = neighbours(i);
      return nb.length === 1 ? nb[0] : -1;
    }

    // ---- draw operations
    function drawAddAt(px, py) {
      const f = cur(), Zn = B.element;
      const pos = f.Z.length ? R().screenToWorld(px, py) : [0, 0, 0];
      edit(mol => {
        const r = BM.chem.addAtom(mol, Zn, pos);
        return adjusted(r.mol, null, [r.index], mol.Z.length);
      }, 'Added ' + el.symbol(Zn), { auto: true });
    }
    function drawGrow(start, px, py, alsoDelete) {
      const f = cur(), Zn = B.element, order = B.order;
      let anchor = start;
      const del = [];
      if (B.adjustH) { const p = hParent(start); if (p >= 0) { anchor = p; del.push(start); } }
      if (alsoDelete >= 0 && del.indexOf(alsoDelete) < 0) del.push(alsoDelete);
      const pa = P(f.xyz, anchor), end = R().screenToWorld(px, py, pa), dir = vsub(end, pa);
      const L = BM.chem.idealBondLength(f.Z[anchor], Zn, order);
      const name = label(f, anchor);
      edit(mol => {
        let m = mol, map = null;
        if (del.length) { const r = BM.chem.deleteAtoms(m, del); m = r.mol; map = r.map; }
        const a = map ? map[anchor] : anchor;
        const pos = vlen(dir) > 0.05 ? vadd(pa, vmul(vunit(dir), L)) : BM.chem.newAtomPosition(m, a, Zn, order);
        const r1 = BM.chem.addAtom(m, Zn, pos);
        m = BM.chem.setBond(r1.mol, a, r1.index, order);
        return adjusted(m, map, [a, r1.index], mol.Z.length);
      }, 'Added ' + el.symbol(Zn) + ' bonded to ' + name, { auto: true });
    }
    function drawConnect(start, target, px, py) {
      const f = cur(), order = B.order;
      let a = start, b = target;
      const del = [];
      if (B.adjustH) {
        const pa = hParent(a);
        if (pa >= 0) { del.push(a); a = pa; }
        const pb = hParent(b);
        if (pb >= 0) { if (pb === a) return drawGrow(start, px, py, b); del.push(b); b = pb; }
      }
      if (a === b) return;
      const names = label(f, a) + '–' + label(f, b);
      edit(mol => {
        let m = mol, map = null;
        if (del.length) { const r = BM.chem.deleteAtoms(m, del); m = r.mol; map = r.map; }
        const A = map ? map[a] : a, C = map ? map[b] : b;
        m = BM.chem.setBond(m, A, C, order);
        return adjusted(m, map, [A, C], mol.Z.length);
      }, names + ': ' + orderName(order) + ' bond', { auto: true });
    }
    function drawSetElement(i) {
      const f = cur(), Zn = B.element;
      if (f.Z[i] === Zn) { status(label(f, i) + ' is already ' + el.symbol(Zn) + '.'); return; }
      const nb = neighbours(i), name = label(f, i);
      edit(mol => {
        const m = Object.assign({}, mol, { Z: mol.Z.slice(), xyz: Float64Array.from(mol.xyz) });
        m.Z[i] = Zn;
        if (nb.length === 1) { // keep a terminal atom at a sensible bond length
          const j = nb[0], pj = P(m.xyz, j), u = vunit(vsub(P(m.xyz, i), pj));
          const o = (m.bonds.find(b => (b[0] === i && b[1] === j) || (b[0] === j && b[1] === i)) || [0, 0, 1])[2];
          if (vlen(u) > 0) { const p = vadd(pj, vmul(u, BM.chem.idealBondLength(m.Z[j], Zn, o))); m.xyz[3 * i] = p[0]; m.xyz[3 * i + 1] = p[1]; m.xyz[3 * i + 2] = p[2]; }
        }
        return adjusted(m, null, [i].concat(nb.filter(j => mol.Z[j] > 1)), mol.Z.length);
      }, name + ' → ' + el.symbol(Zn), { auto: true });
    }
    function drawBond(a, b, order) {
      const f = cur(), names = label(f, a) + '–' + label(f, b);
      edit(mol => adjusted(BM.chem.setBond(mol, a, b, order), null, [a, b], mol.Z.length),
        order ? names + ': ' + orderName(order) + ' bond' : 'Deleted bond ' + names, { auto: true });
    }
    function drawCycleBond(a, b) {
      const o = bondOrder(a, b);
      drawBond(a, b, o >= 2.5 ? 1 : o >= 1.75 ? 3 : 2);
    }
    function drawDeleteAtom(i) {
      const f = cur(), nb = neighbours(i), name = label(f, i);
      const del = [i], aff = [];
      if (f.Z[i] !== 1) {
        for (const j of nb) {
          if (B.adjustH && f.Z[j] === 1 && neighbours(j).length === 1) del.push(j);
          else aff.push(j);
        }
      }
      edit(mol => {
        const r = BM.chem.deleteAtoms(mol, del);
        return adjusted(r.mol, r.map, aff.map(j => r.map[j]).filter(j => j >= 0), mol.Z.length);
      }, 'Deleted ' + name + (del.length > 1 ? ' and ' + (del.length - 1) + ' H' : ''), { auto: true });
    }

    // ---- template operations
    const tplEntry = (kind, id) => BM.templates.get(kind, id);
    function siteOf(i) {
      const f = cur();
      if (isMetal(f.Z[i])) return -1;
      const adj = curAdj(), toMetal = v => adj[v].some(w => isMetal(f.Z[w]));
      if (toMetal(i)) return i;
      const seen = new Set([i]);
      let front = [i];
      while (front.length) {
        const next = [];
        for (const u of front) for (const v of adj[u]) {
          if (seen.has(v) || isMetal(f.Z[v])) continue;
          if (toMetal(v)) return v;
          seen.add(v); next.push(v);
        }
        front = next;
      }
      return -1;
    }
    function tplInsertMolecule(id, at) {
      const e = tplEntry('molecule', id);
      if (!e) return;
      ensureDoc();
      edit(mol => ({ mol: BM.templates.insertMolecule(mol, id, { at }).mol, map: null }), 'Inserted ' + e.name, { auto: true });
    }
    function tplBuildCenter(at) {
      ensureDoc();
      const t = B.tpl, lig = tplEntry('ligand', t.metalLigand);
      edit(mol => ({ mol: BM.templates.buildCenter(mol, { Z: t.metalZ, geometry: t.geometry, ligand: t.metalLigand, at }).mol, map: null }),
        'Built ' + el.symbol(t.metalZ) + ' centre (' + (BM.ff.GEOMETRIES[t.geometry] || {}).name + ', ' + (lig ? lig.name : t.metalLigand) + ')', { auto: true });
    }
    function tplAttach(i, e, groupId) {
      const f = cur(), id = groupId || B.tpl.group, g = tplEntry('group', id);
      if (!g) return;
      let target = i;
      if (f.Z[i] > 1 && !isMetal(f.Z[i])) {
        let full = false;
        try { full = BM.chem.valenceInfo(app.frameMol(f), i).missingH === 0; } catch (err) { full = false; }
        const hs = neighbours(i).filter(j => f.Z[j] === 1 && neighbours(j).length === 1);
        if (full && hs.length) {
          target = hs[0];
          if (e && e.clientX !== undefined) {
            const [px, py] = localXY(e);
            let best = Infinity;
            for (const h of hs) { const s = R().worldToScreen(P(f.xyz, h)), dd = Math.hypot(s.x - px, s.y - py); if (dd < best) { best = dd; target = h; } }
          }
        }
      }
      const name = label(f, i);
      edit(mol => { const r = BM.templates.attachGroup(mol, target, id); return { mol: r.mol, map: r.map || null }; }, 'Attached ' + g.name + ' to ' + name, { auto: true });
    }
    function tplReplace(atoms, ligandId) {
      const id = ligandId || B.tpl.ligand, lig = tplEntry('ligand', id);
      if (!lig) return;
      if (lig.denticity === 2) return tplBidentate(atoms, id);
      const sites = uniq(atoms.map(siteOf));
      if (!sites.length || sites.some(s => s < 0)) { toast('Pick ligand atoms bonded to a metal (not the metal itself).', true); return; }
      edit(mol => {
        let m = mol, map = null;
        for (const s0 of sites) {
          const s = map ? map[s0] : s0;
          if (s < 0) continue;
          const r = BM.templates.replaceLigand(m, s, id);
          map = composeMap(map, r.map, mol.Z.length); m = r.mol;
        }
        return { mol: m, map, select: [] };
      }, (sites.length > 1 ? sites.length + ' ligands' : 'Ligand') + ' replaced by ' + lig.name, { auto: true });
    }
    function tplBidentate(atoms, ligandId) {
      const id = ligandId || B.tpl.ligand, lig = tplEntry('ligand', id);
      const sites = uniq(atoms.map(siteOf));
      if (sites.length !== 2 || sites.some(s => s < 0)) { toast('Select two cis ligand atoms on the same metal for ' + lig.name + '.', true); return; }
      edit(mol => { const r = BM.templates.addBidentate(mol, sites[0], sites[1], id); return { mol: r.mol, map: r.map, select: [] }; },
        'Placed ' + lig.name + ' on two cis sites', { auto: true });
    }
    function templateAtomClick(i, e) {
      const m = B.tpl.mode;
      if (m === 'groups') return tplAttach(i, e);
      if (m === 'ligands') {
        const lig = tplEntry('ligand', B.tpl.ligand);
        if ((e && e.shiftKey) || (lig && lig.denticity === 2)) {
          const d = doc();
          let sel = d.selection.indexOf(i) >= 0 ? d.selection.filter(x => x !== i) : d.selection.concat([i]);
          if (lig && lig.denticity === 2 && sel.length > 2) sel = sel.slice(-2);
          setSelection(sel); renderOptions();
          return;
        }
        return tplReplace([i]);
      }
      status('Click on empty space to place the ' + (m === 'metal' ? 'metal centre' : 'molecule') + '.');
    }
    function templateBackgroundClick(e) {
      const m = B.tpl.mode;
      if (m === 'groups' || m === 'ligands') { status(m === 'groups' ? 'Click an atom (an H or a heavy atom) to attach the group.' : 'Click a ligand atom bonded to a metal.'); return; }
      const f = cur(), [px, py] = localXY(e);
      const at = f && f.Z.length ? R().screenToWorld(px, py) : [0, 0, 0];
      if (m === 'molecules') tplInsertMolecule(B.tpl.molecule, at); else tplBuildCenter(at);
    }

    // ---------------------------------------------------------------- tools
    const T = { navigate: null };
    T.select = {
      cursor: hit => (hit.kind === 'atom' ? 'pointer' : 'crosshair'),
      onDown(ev, hit) {
        // plain presses on atoms stay with the renderer (click toggles, double-click picks the fragment, drag rotates);
        // with Shift or Ctrl/⌘ a box can start anywhere
        if (!doc() || ev.altKey || (hit.kind === 'atom' && selMode(ev) === 'replace')) return false;
        drag = { kind: 'band', hit, x0: hit.px, y0: hit.py, x1: hit.px, y1: hit.py, moved: false, mode: selMode(ev), base: doc().selection.slice() };
        return true;
      },
      onMove(ev, hit) {
        const g = drag; if (!g || g.kind !== 'band') return;
        g.x1 = hit.px; g.y1 = hit.py;
        if (!g.moved && Math.hypot(g.x1 - g.x0, g.y1 - g.y0) > 3) g.moved = true;
        if (g.moved) R().setSelection(combineSel(g.base, R().atomsInRect(g.x0, g.y0, g.x1, g.y1), g.mode));
      },
      onUp(ev, hit) {
        const g = drag; drag = null; if (!g || g.kind !== 'band') return;
        if (ev.type === 'pointercancel') { R().setSelection(g.base); return; }
        g.x1 = hit.px; g.y1 = hit.py;
        if (!g.moved) {
          if (g.hit.kind === 'bond') { const ab = bondAtoms(g.hit.index); if (ab) setSelection(combineSel(g.base, ab, g.mode)); }
          else if (g.hit.kind === 'atom') setSelection(combineSel(g.base, [g.hit.index], g.mode));
          else if (g.mode === 'replace') setSelection([]);
          return;
        }
        setSelection(combineSel(g.base, R().atomsInRect(g.x0, g.y0, g.x1, g.y1), g.mode));
        status(doc().selection.length + ' atoms selected');
      }
    };
    T.draw = {
      cursor: hit => (hit.kind === 'none' ? 'crosshair' : 'pointer'),
      onHover(ev, hit) { const k = hit.kind === 'bond' ? hit.index : -1; if (k !== R().highlightedBond) R().highlightBond(k); },
      onDown(ev, hit) {
        if (!doc() || hit.kind === 'none') return false;
        R().highlightBond(-1);
        drag = { kind: 'draw', hit, start: hit.kind === 'atom' ? hit.index : -1, x0: hit.px, y0: hit.py, px: hit.px, py: hit.py, moved: false, alt: ev.altKey, target: -1 };
        return true;
      },
      onMove(ev, hit) {
        const g = drag; if (!g || g.kind !== 'draw') return;
        g.px = hit.px; g.py = hit.py;
        if (!g.moved && Math.hypot(g.px - g.x0, g.py - g.y0) > 4) g.moved = true;
        g.target = hit.kind === 'atom' && hit.index !== g.start ? hit.index : -1;
      },
      onUp(ev, hit) {
        const g = drag; drag = null; if (!g || g.kind !== 'draw' || ev.type === 'pointercancel') return;
        g.px = hit.px; g.py = hit.py;
        if (!g.moved) {
          const del = g.alt || ev.altKey;
          if (g.hit.kind === 'atom') { if (del) drawDeleteAtom(g.start); else drawSetElement(g.start); }
          else if (g.hit.kind === 'bond') { const ab = bondAtoms(g.hit.index); if (ab) { if (del) drawBond(ab[0], ab[1], 0); else drawCycleBond(ab[0], ab[1]); } }
          return;
        }
        if (g.start < 0) return;
        const t = hit.kind === 'atom' && hit.index !== g.start ? hit.index : -1;
        if (t >= 0) drawConnect(g.start, t, g.px, g.py); else drawGrow(g.start, g.px, g.py, -1);
      }
    };
    T.template = {
      cursor: hit => { const m = B.tpl.mode; return m === 'groups' || m === 'ligands' ? (hit.kind === 'atom' ? 'pointer' : '') : (hit.kind === 'none' ? 'copy' : ''); },
      onDown() { return false; }
    };
    T.bond = {
      cursor: hit => (drag && drag.kind === 'bond' ? 'grabbing' : hit.kind === 'bond' ? 'grab' : ''),
      onHover(ev, hit) { const k = hit.kind === 'bond' ? hit.index : -1; if (k !== R().highlightedBond) R().highlightBond(k); },
      onDown(ev, hit) {
        if (!doc() || hit.kind !== 'bond') return false;
        const ab = bondAtoms(hit.index); if (!ab) return false;
        const adj = curAdj(), sa = geom.side(adj, ab[0], ab[1]), sb = geom.side(adj, ab[1], ab[0]);
        if (sa.ring || sb.ring) { toast('This bond is part of a ring, so it cannot be rotated or stretched rigidly.'); return true; }
        const moveA = sa.atoms.length < sb.atoms.length;
        const p = moveA ? ab[0] : ab[1], q = moveA ? ab[1] : ab[0];
        drag = { kind: 'bond', k: hit.index, p, q, atoms: moveA ? sa.atoms : sb.atoms, x0: hit.px, y0: hit.py, lx: hit.px, ly: hit.py, moved: false, stretch: ev.shiftKey, turned: 0, ref: dihedralRef(adj, q, p), f: null };
        R().highlightBond(hit.index);
        bondChip(drag);
        return true;
      },
      onMove(ev, hit) {
        const g = drag; if (!g || g.kind !== 'bond') return;
        if (!g.moved) {
          if (Math.hypot(hit.px - g.x0, hit.py - g.y0) < 2) return;
          g.moved = true; g.f = app.beginContinuousEdit(doc());
        }
        const xyz = g.f.xyz, pp = P(xyz, g.p), pq = P(xyz, g.q), axis = vsub(pp, pq);
        if (g.stretch) {
          const w0 = R().screenToWorld(g.lx, g.ly, pp), w1 = R().screenToWorld(hit.px, hit.py, pp), u = vunit(axis), L = vlen(axis);
          let dl = vdot(vsub(w1, w0), u);
          if (L + dl < 0.5) dl = 0.5 - L;
          for (const a of g.atoms) { xyz[3 * a] += u[0] * dl; xyz[3 * a + 1] += u[1] * dl; xyz[3 * a + 2] += u[2] * dl; }
        } else {
          const sp = R().worldToScreen(pp), sq = R().worldToScreen(pq), dx = hit.px - g.lx, dy = hit.py - g.ly;
          let ux = sp.x - sq.x, uy = sp.y - sq.y;
          const ul = Math.hypot(ux, uy);
          let amt = dx;
          if (ul > 8) { ux /= ul; uy /= ul; amt = -dx * uy + dy * ux; }
          const th = amt * 0.012;
          if (th) { geom.rotateAtoms(xyz, g.atoms, pq, axis, th); g.turned += th; }
        }
        g.lx = hit.px; g.ly = hit.py;
        R().updatePositions(xyz);
        bondChip(g);
      },
      onUp(ev, hit) {
        const g = drag; drag = null; if (!g || g.kind !== 'bond') return;
        chip = null;
        if (!g.moved) { setSelection([g.q, g.p]); return; }
        const f = cur();
        app.endContinuousEdit(doc(), g.stretch ? 'Bond ' + label(f, g.q) + '–' + label(f, g.p) + ' stretched to ' + geom.dist(f.xyz, g.q, g.p).toFixed(3) + ' Å'
          : 'Rotated about ' + label(f, g.q) + '–' + label(f, g.p) + ' by ' + (g.turned * 180 / Math.PI).toFixed(1) + '°');
      }
    };
    T.manipulate = {
      cursor: hit => (drag && drag.kind === 'move' ? 'grabbing' : hit.kind === 'atom' ? 'grab' : ''),
      onDown(ev, hit) {
        if (!doc() || hit.kind !== 'atom') return false;
        const d = doc(), i = hit.index;
        const set = d.selection.indexOf(i) >= 0 ? d.selection.slice() : fragmentOf(i);
        drag = { kind: 'move', i, set, rotate: ev.shiftKey, x0: hit.px, y0: hit.py, lx: hit.px, ly: hit.py, moved: false, f: null, plane: P(R().xyz, i) };
        return true;
      },
      onMove(ev, hit) {
        const g = drag; if (!g || g.kind !== 'move') return;
        if (!g.moved) {
          if (Math.hypot(hit.px - g.x0, hit.py - g.y0) < 3) return;
          g.moved = true; g.f = app.beginContinuousEdit(doc());
          if (B.autoOpt) relaxStart('auto', { frozen: g.set, drag: true });
        }
        const xyz = g.f.xyz;
        if (g.rotate) {
          const dx = hit.px - g.lx, dy = hit.py - g.ly, ang = Math.hypot(dx, dy) * 0.01;
          if (ang > 0) {
            const m = BM.quat.mat(R().q), v = [dy, dx, 0];
            const axis = [m[0] * v[0] + m[3] * v[1] + m[6] * v[2], m[1] * v[0] + m[4] * v[1] + m[7] * v[2], m[2] * v[0] + m[5] * v[1] + m[8] * v[2]];
            geom.rotateAtoms(xyz, g.set, centroid(xyz, g.set), axis, ang);
          }
        } else {
          const w0 = R().screenToWorld(g.lx, g.ly, g.plane), w1 = R().screenToWorld(hit.px, hit.py, g.plane), dv = vsub(w1, w0);
          for (const a of g.set) { xyz[3 * a] += dv[0]; xyz[3 * a + 1] += dv[1]; xyz[3 * a + 2] += dv[2]; }
          g.plane = vadd(g.plane, dv);
        }
        g.lx = hit.px; g.ly = hit.py;
        R().updatePositions(xyz);
        app.updateMeasure(true);
        if (relax && relax.mode === 'auto') relaxPositions(xyz);
      },
      onUp(ev, hit) {
        const g = drag; drag = null; if (!g || g.kind !== 'move') return;
        if (!g.moved) { if (ev.type !== 'pointercancel') { const d = doc(); setSelection(d.selection.indexOf(g.i) >= 0 ? d.selection.filter(x => x !== g.i) : d.selection.concat([g.i])); } return; }
        app.endContinuousEdit(doc(), (g.rotate ? 'Rotated ' : 'Moved ') + g.set.length + ' atom' + (g.set.length > 1 ? 's' : ''));
        if (relax && relax.mode === 'auto') relaxRelease();
      }
    };
    function dihedralRef(adj, q, p) {
      const f = cur(), pick = (a, skip) => { const l = adj[a].filter(x => x !== skip); return l.find(x => f.Z[x] > 1) ?? l[0]; };
      const a = pick(q, p), b = pick(p, q);
      return a === undefined || b === undefined ? null : [a, q, p, b];
    }
    function bondChip(g) {
      const f = cur(), xyz = R().xyz, names = label(f, g.q) + '–' + label(f, g.p);
      if (g.stretch) chip = `<span class="atoms">${esc(names)}</span> <span>length</span><span class="val">${geom.dist(xyz, g.q, g.p).toFixed(4)} Å</span>`;
      else if (g.ref) chip = `<span class="atoms">${esc(g.ref.map(i => label(f, i)).join(' – '))}</span> <span>dihedral</span><span class="val">${geom.dihedral(xyz, ...g.ref).toFixed(2)}°</span>`;
      else chip = `<span class="atoms">${esc(names)}</span> <span>rotated</span><span class="val">${(g.turned * 180 / Math.PI).toFixed(1)}°</span>`;
      app.updateMeasure(true);
    }

    // ---------------------------------------------------------------- clean-up (force field) runs
    function relaxStart(mode, opts) {
      opts = opts || {};
      const d = doc(), f = cur();
      if (!d || !f || !f.Z.length) return false;
      const n = f.Z.length;
      if (mode === 'auto' && n > LIMIT.auto) { status('Auto-optimize is off for structures above ' + LIMIT.auto + ' atoms.'); return false; }
      if (mode === 'optimize' && n > LIMIT.clean) { toast('Geometry clean-up is limited to ' + LIMIT.clean + ' atoms (this structure has ' + n + ').', true); return false; }
      relaxStop(true);
      const frozen = uniq(d.frozen.concat(opts.frozen || [])).filter(i => i < n);
      let opt;
      try { opt = BM.ff.createOptimizer(app.frameMol(f), { frozen }); } catch (e) { toast('Clean-up failed: ' + e.message, true); return false; }
      const s0 = opt.step(0);
      relax = { d, f, opt, mode, steps: 0, e0: s0.energy, e: s0.energy, rms: s0.rmsGrad, stall: 0, max: mode === 'optimize' ? 2000 : 800, dragging: !!opts.drag, raf: 0, t0: performance.now() };
      relax.raf = requestAnimationFrame(relaxTick);
      if (mode === 'optimize') refreshCleanup();
      return true;
    }
    function relaxTick() {
      const r = relax;
      if (!r) return;
      r.raf = 0;
      if (doc() !== r.d || cur() !== r.f) { relaxStop(true); return; }
      if ((r.stall >= 3 || r.converged) && r.dragging) return; // idle until the drag moves atoms again
      const t0 = performance.now();
      let s = null;
      do {
        s = r.opt.step(1);
        r.steps++;
        r.stall = Math.abs(s.energy - r.e) < 1e-9 ? r.stall + 1 : 0;
        r.e = s.energy; r.rms = s.rmsGrad; r.converged = s.converged;
      } while (!s.converged && r.stall < 3 && r.steps < r.max && performance.now() - t0 < 9);
      r.f.xyz = s.xyz;
      R().updatePositions(r.f.xyz);
      app.updateMeasure(true);
      if (r.mode === 'optimize') refreshCleanup();
      const done = s.converged || r.stall >= 3 || r.steps >= r.max;
      if (done && !r.dragging) relaxFinish();
      else if (!done || !r.dragging) r.raf = requestAnimationFrame(relaxTick);
    }
    function relaxPositions(xyz) {
      const r = relax;
      if (!r) return;
      try { r.opt.setPositions(xyz); } catch (e) { relaxStop(true); return; }
      r.stall = 0; r.converged = false;
      if (!r.raf) r.raf = requestAnimationFrame(relaxTick);
    }
    function relaxRelease() {
      const r = relax;
      if (!r) return;
      r.dragging = false;
      try { r.opt.setFrozen(r.d.frozen.filter(i => i < r.f.Z.length)); } catch (e) { relaxStop(true); return; }
      r.stall = 0; r.converged = false; r.steps = 0;
      if (!r.raf) r.raf = requestAnimationFrame(relaxTick);
    }
    function relaxFinish() {
      const r = relax;
      if (!r) return;
      relax = null;
      if (r.raf) cancelAnimationFrame(r.raf);
      if (r.mode === 'optimize') {
        r.d._ff = { e0: r.e0, e: r.e, steps: r.steps, converged: !!r.converged, stopped: false };
        status('Clean-up ' + (r.converged ? 'converged' : 'finished') + ' after ' + r.steps + ' step' + (r.steps === 1 ? '' : 's') + ': E = ' + r.e.toFixed(3) + ' kcal/mol (clean-up force field, not a QM energy)');
      }
      app.afterPositions(r.d);
    }
    /** Stop a running clean-up (keeps the geometry reached so far). */
    function relaxStop(quiet) {
      const r = relax;
      if (!r) return false;
      relax = null;
      if (r.raf) cancelAnimationFrame(r.raf);
      if (r.mode === 'optimize') {
        r.d._ff = { e0: r.e0, e: r.e, steps: r.steps, converged: !!r.converged, stopped: true };
        if (!quiet) status('Clean-up stopped after ' + r.steps + ' steps: E = ' + r.e.toFixed(3) + ' kcal/mol');
      }
      app.afterPositions(r.d);
      return true;
    }
    function optimizeGeometry() {
      const d = doc(), f0 = cur();
      if (!d || !f0) return;
      if (relax && relax.mode === 'optimize') return;
      if (!f0.Z.length) { toast('There is nothing to optimize yet.', true); return; }
      if (f0.Z.length > LIMIT.clean) { toast('Geometry clean-up is limited to ' + LIMIT.clean + ' atoms (this structure has ' + f0.Z.length + ').', true); return; }
      app.beginContinuousEdit(d);
      relaxStart('optimize');
    }
    function refreshCleanup() {
      const d = doc(); if (!d) return;
      const parts = { clean: cleanBtnDyn, cinfo: cleanInfoDyn, freeze: freezeDyn };
      for (const k in parts) { const box = $('#panel [data-bp="' + k + '"]'); if (box) box.innerHTML = parts[k](d); }
    }

    // ---------------------------------------------------------------- tool switching, palette, options bar
    function setTool(id, quiet) {
      if (!TOOLS.some(t => t.id === id)) return;
      if (relax) relaxStop(true);
      drag = null; chip = null;
      closePopup();
      B.tool = id;
      const r = R();
      if (r) { r.highlightBond(-1); r.tool = T[id] || null; }
      renderToolbox(); renderOptions(); showHint(!quiet);
      if (!quiet) app.updateMeasure(true);
      app.savePrefs();
    }
    function renderToolbox() {
      const box = $('#toolbox');
      if (!box) return;
      if (!box.children.length) {
        box.innerHTML = TOOLS.map(t => `<button data-tool="${t.id}" title="${t.name} (${t.key}) — ${t.tip}" aria-label="${t.name}"><svg class="ic"><use href="#${t.icon}"/></svg></button>`).join('');
        box.addEventListener('click', e => { const b = e.target.closest('[data-tool]'); if (b) setTool(b.dataset.tool); });
      }
      $$('[data-tool]', box).forEach(b => b.classList.toggle('on', b.dataset.tool === B.tool));
      box.classList.toggle('hidden', !doc() || !!doc().pending);
    }
    const autoChk = () => `<label class="to-chk" title="Auto clean-up: relax with the clean-up force field after each change (up to ${LIMIT.auto} atoms)"><input type="checkbox" data-opt="autoOpt" ${B.autoOpt ? 'checked' : ''}>Auto-clean</label>`;
    function renderOptions() {
      const bar = $('#toolopts');
      if (!bar) return;
      const d = doc();
      let h = '';
      if (d && B.tool === 'draw') {
        const inQuick = QUICK.indexOf(B.element) >= 0;
        h = `<span class="to-group">${QUICK.map(z => `<button class="to-el ${z === B.element ? 'on' : ''}" data-el="${z}" title="${el.name(z)}">${el.symbol(z)}</button>`).join('')}<button class="to-el more ${inQuick ? '' : 'on'}" data-pt="draw" title="Other element (periodic table)">${inQuick ? '…' : el.symbol(B.element)}</button></span>
          <span class="to-sep"></span>
          <span class="to-group" title="Order of new bonds">${[1, 2, 3].map(o => `<button class="to-ord ${o === B.order ? 'on' : ''}" data-order="${o}" title="${ORDER_NAME[o][0].toUpperCase() + ORDER_NAME[o].slice(1)} bond for new bonds">${['–', '=', '≡'][o - 1]}</button>`).join('')}</span>
          <span class="to-sep"></span>
          <label class="to-chk" title="Adjust hydrogens after each change to keep normal valences"><input type="checkbox" data-opt="adjustH" ${B.adjustH ? 'checked' : ''}>±H</label>${autoChk()}`;
      } else if (d && B.tool === 'template') {
        const t = B.tpl;
        h = `<span class="to-seg">${TPL_MODES.map(([id, n, tip]) => `<button class="${t.mode === id ? 'on' : ''}" data-mode="${id}" title="${tip}">${n}</button>`).join('')}</span><span class="to-sep"></span>`;
        if (t.mode === 'molecules') h += pickBtn('molecule', t.molecule);
        else if (t.mode === 'groups') h += pickBtn('group', t.group);
        else if (t.mode === 'metal') {
          h += `<button class="to-el on" data-pt="metal" title="Metal element (periodic table)">${el.symbol(t.metalZ)}</button>
            <select class="to-select" data-tpl="geometry" title="Coordination geometry">${BM.templates.geometries.map(g => `<option value="${g.id}" ${g.id === t.geometry ? 'selected' : ''}>${esc(g.name)} (${g.cn})</option>`).join('')}</select>
            <select class="to-select" data-tpl="metalLigand" title="Ligand on every site">${ligandOptions(t.metalLigand, true)}</select>`;
        } else {
          const lig = tplEntry('ligand', t.ligand), nSel = d.selection.length;
          h += pickBtn('ligand', t.ligand) + `<span class="to-lbl" title="Ligand charge">q ${lig ? chargeText(lig.charge) : '?'}</span>
            <button class="to-btn" data-tact="apply-sites" ${nSel ? '' : 'disabled'} title="${lig && lig.denticity === 2 ? 'Place the ligand on the two selected cis sites' : 'Replace the ligand at every selected site'}">Apply${nSel ? ' (' + nSel + ')' : ''}</button>`;
        }
      } else if (d && B.tool === 'select') {
        const f = cur(), els = uniq(f.Z.filter(z => z >= 0)).sort((a, b) => a - b);
        h = `<span class="to-group">${[['all', 'All', 'Select all atoms (Ctrl/⌘+A)'], ['none', 'None', 'Clear the selection'], ['invert', 'Invert', 'Invert the selection'], ['fragment', 'Fragment', 'Extend the selection to whole molecules']].map(([k, n, tip]) => `<button class="to-btn" data-sel="${k}" title="${tip}">${n}</button>`).join('')}</span>
          <span class="to-sep"></span>
          <select class="to-select" data-sel-el title="Select every atom of one element"><option value="">Element…</option>${els.map(z => `<option value="${z}">${z > 0 ? el.symbol(z) + ' · ' + el.name(z) : 'X · dummy'}</option>`).join('')}</select>
          <span class="to-sep"></span>
          <input class="to-num" type="number" data-expand min="0.1" step="0.5" value="${esc(B.inputs.expand)}" title="Distance (Å) for Expand"><span class="to-lbl">Å</span><button class="to-btn" data-sel="expand" title="Add atoms within this distance of the selection">Expand</button>`;
      } else if (d && B.tool === 'bond') {
        h = '<span class="to-lbl" title="Drag a bond to rotate its smaller side; Shift-drag to stretch it"><b>Drag</b> rotate · <b>Shift</b> stretch</span>';
      } else if (d && B.tool === 'manipulate') {
        h = `<span class="to-lbl" title="Drag an atom to move its fragment (or the selection); Shift-drag to rotate it"><b>Drag</b> move · <b>Shift</b> rotate</span><span class="to-sep"></span>${autoChk()}`;
      }
      bar.innerHTML = h;
      bar.classList.toggle('hidden', !h);
    }
    function pickBtn(kind, id) {
      const e = tplEntry(kind, id);
      return `<button class="to-pick" data-pick="${kind}" title="Choose a ${kind} (searchable list)"><span>${esc(e ? e.name : id)}</span><svg class="ic sm"><use href="#i-chev"/></svg></button>`;
    }
    function ligandOptions(current, withEmpty) {
      const L = BM.templates.ligands;
      const opt = l => `<option value="${l.id}" ${l.id === current ? 'selected' : ''}>${esc(l.name)} ${chargeText(l.charge) === '0' ? '' : '(' + chargeText(l.charge) + ')'}</option>`;
      let h = `<optgroup label="Monodentate">${L.filter(l => l.denticity === 1 && l.id !== 'X').map(opt).join('')}</optgroup><optgroup label="Bidentate">${L.filter(l => l.denticity === 2).map(opt).join('')}</optgroup>`;
      if (withEmpty) h += `<optgroup label="Other"><option value="X" ${current === 'X' ? 'selected' : ''}>Empty sites (X)</option></optgroup>`;
      return h;
    }
    function bindOptions() {
      const bar = $('#toolopts');
      bar.addEventListener('click', e => {
        const t = e.target;
        const b = t.closest('button');
        if (!b || b.disabled) return;
        if (b.dataset.el) { B.element = +b.dataset.el; renderOptions(); showHint(true); app.savePrefs(); return; }
        if (b.dataset.order) { B.order = +b.dataset.order; renderOptions(); app.savePrefs(); return; }
        if (b.dataset.pt === 'draw') { openPeriodic(b, B.element, z => { B.element = z; renderOptions(); showHint(true); app.savePrefs(); }); return; }
        if (b.dataset.pt === 'metal') { openPeriodic(b, B.tpl.metalZ, z => { B.tpl.metalZ = z; renderOptions(); showHint(true); app.savePrefs(); }); return; }
        if (b.dataset.mode) { B.tpl.mode = b.dataset.mode; renderOptions(); showHint(true); app.savePrefs(); return; }
        if (b.dataset.pick) { openPicker(b, b.dataset.pick); return; }
        if (b.dataset.sel) { selectAction(b.dataset.sel); return; }
        if (b.dataset.tact === 'apply-sites') { tplReplace(doc().selection.slice()); }
      });
      bar.addEventListener('change', e => {
        const t = e.target;
        if (t.dataset.opt) { B[t.dataset.opt] = t.checked; app.savePrefs(); if (t.dataset.opt === 'autoOpt') { if (!t.checked && relax && relax.mode === 'auto') relaxStop(true); if (doc() && doc().tab === 'build') refreshCleanup(); } return; }
        if (t.dataset.tpl) { B.tpl[t.dataset.tpl] = t.value; showHint(true); app.savePrefs(); return; }
        if (t.dataset.selEl !== undefined && t.value !== '') { const z = +t.value, f = cur(); setSelection(f.Z.map((x, i) => (x === z ? i : -1)).filter(i => i >= 0)); t.value = ''; }
      });
      bar.addEventListener('input', e => { if (e.target.dataset.expand !== undefined) B.inputs.expand = e.target.value; });
    }
    function selectAction(k) {
      const d = doc(); if (!d) return;
      const f = cur(), n = f.Z.length, all = Array.from({ length: n }, (_, i) => i);
      if (k === 'all') setSelection(all);
      else if (k === 'none') setSelection([]);
      else if (k === 'invert') { const s = new Set(d.selection); setSelection(all.filter(i => !s.has(i))); }
      else if (k === 'fragment') {
        if (!d.selection.length) { status('Select an atom first; Fragment extends the selection to whole molecules.'); return; }
        const adj = curAdj(), seen = new Set(d.selection), stack = d.selection.slice();
        while (stack.length) { const u = stack.pop(); for (const v of adj[u]) if (!seen.has(v)) { seen.add(v); stack.push(v); } }
        setSelection([...seen]);
      } else if (k === 'expand') {
        const r = parseFloat(B.inputs.expand);
        if (!(r > 0)) { toast('Enter a distance in Å.', true); return; }
        if (!d.selection.length) { status('Select atoms first.'); return; }
        const xyz = f.xyz, cell = Math.max(r, 0.5), grid = new Map(), key = (a, b, c) => a + ',' + b + ',' + c, cf = v => Math.floor(v / cell);
        for (const i of d.selection) { const k2 = key(cf(xyz[3 * i]), cf(xyz[3 * i + 1]), cf(xyz[3 * i + 2])); (grid.get(k2) || grid.set(k2, []).get(k2)).push(i); }
        const out = new Set(d.selection);
        for (let i = 0; i < n; i++) {
          if (out.has(i)) continue;
          const cx = cf(xyz[3 * i]), cy = cf(xyz[3 * i + 1]), cz = cf(xyz[3 * i + 2]);
          let hit = false;
          for (let a = -1; a <= 1 && !hit; a++) for (let b = -1; b <= 1 && !hit; b++) for (let c = -1; c <= 1 && !hit; c++) {
            const l = grid.get(key(cx + a, cy + b, cz + c));
            if (l) for (const j of l) if (geom.dist(xyz, i, j) <= r) { hit = true; break; }
          }
          if (hit) out.add(i);
        }
        setSelection([...out]);
        status(out.size + ' atoms within ' + r + ' Å of the selection');
      }
    }

    // ---------------------------------------------------------------- hints
    function hintText() {
      const d = doc();
      const tool = TOOLS.find(t => t.id === B.tool) || TOOLS[0];
      const head = tool.name + ' — ';
      switch (B.tool) {
        case 'select': return head + 'drag on empty space: box select (Shift adds, Ctrl/⌘ removes) · click an atom: toggle · double-click: whole fragment · Ctrl/⌘+A: all';
        case 'draw': return head + 'click empty space: add ' + el.symbol(B.element) + ' · drag from an atom: new bonded atom · click an atom: make it ' + el.symbol(B.element) + ' · click a bond: cycle 1→2→3 · Alt/Option-click: delete';
        case 'template': {
          const t = B.tpl;
          if (t.mode === 'molecules') { const e = tplEntry('molecule', t.molecule); return head + 'click empty space to insert ' + (e ? e.name : t.molecule); }
          if (t.mode === 'groups') { const e = tplEntry('group', t.group); return head + 'click an atom (an H or a heavy atom) to attach ' + (e ? e.name : t.group); }
          const id = t.mode === 'metal' ? t.metalLigand : t.ligand, lig = tplEntry('ligand', id);
          let s = head + (t.mode === 'metal' ? 'click empty space to build ' + el.symbol(t.metalZ) + ' (' + ((BM.ff.GEOMETRIES[t.geometry] || {}).name || t.geometry) + ')'
            : 'click a ligand atom to replace it' + (lig && lig.denticity === 2 ? ' — pick two cis sites, then Apply' : ' · Shift-click picks sites'));
          if (lig) s += ' · ligand ' + lig.id + ' charge ' + chargeText(lig.charge);
          if (d && cur().Z.length) {
            const sc = suggestCharge(app.topoFrame(d));
            const ox = oxText(sc.metals);
            s += sc.determinable ? ' · suggested total charge ' + chargeText(sc.total) + (ox ? ' (assuming ' + ox + ')' : '') : ' · total charge: set the metal oxidation state as its formal charge';
          }
          return s;
        }
        case 'bond': return head + 'drag a bond: rotate the smaller side · Shift-drag: stretch it · ring bonds cannot be rotated';
        case 'manipulate': return head + 'drag an atom: move its fragment, or the selection if the atom is selected · Shift-drag: rotate · each drag is one undo step';
        default: return head + 'drag: rotate · right-drag or Shift-drag: pan · wheel: zoom · click atoms: measure (2–4) · double-click: centre';
      }
    }
    function showHint(now) { app.setHint(doc() ? hintText() : '', now); }

    // ---------------------------------------------------------------- popups (periodic table, template picker)
    function closePopup() {
      if (!popup) return;
      popup.el.remove();
      document.removeEventListener('pointerdown', popup.outside, true);
      popup = null;
    }
    function openPopup(anchor, html, cls) {
      closePopup();
      const box = document.createElement('div');
      box.className = 'popup ' + (cls || '');
      box.innerHTML = html;
      document.body.appendChild(box);
      const r = anchor.getBoundingClientRect(), w = box.offsetWidth, h = box.offsetHeight;
      let left = Math.min(r.left, window.innerWidth - w - 8), top = r.bottom + 6;
      if (top + h > window.innerHeight - 8) top = Math.max(8, r.top - h - 6);
      box.style.left = Math.max(8, left) + 'px'; box.style.top = top + 'px';
      const outside = e => { if (!box.contains(e.target) && !anchor.contains(e.target)) closePopup(); };
      document.addEventListener('pointerdown', outside, true);
      popup = { el: box, outside, anchor };
      return box;
    }
    function openPeriodic(anchor, current, onPick) {
      if (popup && popup.anchor === anchor) { closePopup(); return; }
      let cells = '';
      for (let z = 1; z <= 118; z++) {
        const [row, col] = ptCell(z), c = el.color[z] || [0.6, 0.6, 0.6];
        cells += `<button class="pt-el ${z === current ? 'on' : ''} ${isMetal(z) ? 'metal' : ''}" data-z="${z}" style="grid-row:${row};grid-column:${col};--elc:${rgb(c)}" title="${z} · ${el.name(z)}">${el.symbol(z)}</button>`;
      }
      const box = openPopup(anchor, `<div class="pt-grid">${cells}<div class="pt-gap" style="grid-row:8;grid-column:1 / span 18"></div></div><div class="pt-foot"><span data-pt-info>${current} · ${el.symbol(current)} · ${el.name(current)}</span><span class="dim">Esc to close</span></div>`, 'pt');
      box.addEventListener('mouseover', e => { const b = e.target.closest('[data-z]'); if (b) $('[data-pt-info]', box).textContent = b.dataset.z + ' · ' + el.symbol(+b.dataset.z) + ' · ' + el.name(+b.dataset.z); });
      box.addEventListener('click', e => { const b = e.target.closest('[data-z]'); if (!b) return; closePopup(); onPick(+b.dataset.z); });
    }
    function libraryEntries(kinds) {
      const T2 = BM.templates, out = [];
      if (kinds.indexOf('molecule') >= 0) for (const m of T2.molecules) out.push({ kind: 'molecule', id: m.id, name: m.name, cat: m.category });
      if (kinds.indexOf('group') >= 0) for (const g of T2.groups) out.push({ kind: 'group', id: g.id, name: g.name, cat: g.category + ' groups' });
      if (kinds.indexOf('ligand') >= 0) for (const l of T2.ligands) out.push({ kind: 'ligand', id: l.id, name: l.name, cat: l.category, extra: 'charge ' + chargeText(l.charge) });
      return out;
    }
    function filterEntries(list, q) {
      const words = String(q || '').toLowerCase().split(/\s+/).filter(Boolean);
      if (!words.length) return list;
      return list.filter(e => { const s = (e.name + ' ' + e.id + ' ' + e.cat + ' ' + e.kind).toLowerCase(); return words.every(w => s.indexOf(w) >= 0); });
    }
    function entryRows(list, current, picker) {
      if (!list.length) return '<div class="lib-empty">No matches</div>';
      let h = '', cat = null;
      for (const e of list) {
        if (e.cat !== cat) { cat = e.cat; h += `<div class="lib-cat">${esc(cat)}</div>`; }
        const attr = picker ? `data-pickid="${esc(e.id)}"` : `data-act="b-lib" data-lib="${e.kind}:${esc(e.id)}"`;
        h += `<button class="lib-item ${e.id === current ? 'on' : ''}" ${attr} title="${esc(e.name)}"><span>${esc(e.name)}</span>${e.extra ? `<span class="dim">${esc(e.extra)}</span>` : ''}</button>`;
      }
      return h;
    }
    function openPicker(anchor, kind) {
      if (popup && popup.anchor === anchor) { closePopup(); return; }
      const all = libraryEntries([kind]), current = B.tpl[kind === 'molecule' ? 'molecule' : kind === 'group' ? 'group' : 'ligand'];
      const box = openPopup(anchor, `<input type="search" class="pick-q" placeholder="Search ${kind}s…" spellcheck="false"><div class="pick-list">${entryRows(all, current, true)}</div>`, 'pick');
      const q = $('.pick-q', box), listEl = $('.pick-list', box);
      q.addEventListener('input', () => { listEl.innerHTML = entryRows(filterEntries(all, q.value), current, true); });
      q.addEventListener('keydown', e => {
        if (e.key === 'Enter') { const first = $('.lib-item', listEl); if (first) first.click(); }
        if (e.key === 'Escape') { e.stopPropagation(); closePopup(); }
      });
      box.addEventListener('click', e => {
        const b = e.target.closest('[data-pickid]'); if (!b) return;
        e.stopPropagation();
        B.tpl[kind] = b.dataset.pickid;
        closePopup(); renderOptions(); showHint(true); app.savePrefs();
      });
      setTimeout(() => q.focus(), 0);
    }

    // ---------------------------------------------------------------- library (Build tab) and SMILES
    function useLibrary(kind, id) {
      const d = doc();
      if (kind === 'molecule') {
        const e = tplEntry('molecule', id); if (!e) return;
        B.tpl.molecule = id;
        const m = e.mol, rNew = molRadius(m.xyz, m.Z.length).r;
        tplInsertMolecule(id, d ? besidePoint(rNew) : [0, 0, 0]);
        return;
      }
      if (!d || !cur().Z.length) { toast('Build or open a structure first, then attach the ' + kind + ' to an atom.', true); return; }
      const sel = d.selection;
      if (kind === 'group') {
        B.tpl.group = id;
        if (sel.length === 1) { tplAttach(sel[0], null, id); return; }
        B.tpl.mode = 'groups'; setTool('template');
        status('Click an atom (an H or a heavy atom) to attach ' + tplEntry('group', id).name + '.');
        return;
      }
      const lig = tplEntry('ligand', id); if (!lig) return;
      B.tpl.ligand = id;
      if ((lig.denticity === 1 && sel.length >= 1) || (lig.denticity === 2 && sel.length === 2)) { tplReplace(sel.slice(), id); return; }
      B.tpl.mode = 'ligands'; setTool('template');
      status(lig.denticity === 2 ? 'Pick two cis ligand atoms, then press Apply to selected sites.' : 'Click a ligand atom to replace it by ' + lig.name + '.');
    }
    function insertSmiles(text) {
      const s = String(text || '').trim();
      if (!s) { toast('Type a SMILES string first, e.g. CC(=O)Oc1ccccc1C(=O)O', true); return false; }
      let r;
      try { r = BM.chem.smilesTo3D(s, { seed: 1 }); } catch (e) { toast('SMILES: ' + e.message, true); return false; }
      const add = r.mol, arom = new Set(r.aromaticBonds || []), n = add.Z.length;
      const c = centroid(add.xyz, Array.from({ length: n }, (_, i) => i));
      const rNew = molRadius(add.xyz, n).r;
      ensureDoc();
      const at = besidePoint(rNew);
      const formula = geom.formula(add.Z);
      const ok = edit(mol => {
        const n0 = mol.Z.length;
        const xyz = new Float64Array(3 * (n0 + n));
        xyz.set(mol.xyz);
        for (let i = 0; i < n; i++) for (let k = 0; k < 3; k++) xyz[3 * (n0 + i) + k] = add.xyz[3 * i + k] - c[k] + at[k];
        const bonds = mol.bonds.concat(add.bonds.map((b, k) => [b[0] + n0, b[1] + n0, arom.has(k) ? 1.5 : b[2]]));
        return { mol: { Z: mol.Z.concat(Array.from(add.Z)), xyz, bonds, formalCharges: mol.formalCharges.concat(Array.from(add.formalCharges || new Array(n).fill(0), q => q || 0)) }, map: null };
      }, 'Inserted ' + formula + ' from SMILES' + (r.warnings && r.warnings.length ? ' (' + r.warnings.join('; ') + ')' : ''), {});
      if (ok && r.warnings && r.warnings.length) toast('SMILES: ' + r.warnings.join('; '));
      if (ok) { B.inputs.smiles = ''; const inp = $('#bp-smiles'); if (inp) inp.value = ''; }
      return ok;
    }
    function currentSmiles() {
      const d = doc(), f = cur();
      if (!f || !f.Z.length) return '';
      const t = app.topoFrame(d);
      return BM.chem.toSmiles({ Z: t.Z, bonds: t.bonds, formalCharges: t.formalCharges });
    }
    function newMolecule(opts) {
      const d = app.addBuiltDoc();
      if (!opts || !opts.keepTool) setTool('draw');
      return d;
    }

    // ---------------------------------------------------------------- Build tab
    const fold = id => (B.fold[id] ? 'folded' : '');
    function sec(id, title, body, extra) {
      return `<div class="sec bsec ${fold(id)}" data-sec="${id}"><h3 class="fold" data-act="b-fold" data-sec="${id}"><span>${title}</span>${extra || ''}<svg class="ic sm chev"><use href="#i-chev"/></svg></h3><div class="bbody">${body}</div></div>`;
    }
    const btn = (act, text, enabled, title, cls) => `<button class="btn small ${cls || ''}" data-act="${act}" ${enabled === false ? 'disabled' : ''} ${title ? `title="${esc(title)}"` : ''}>${text}</button>`;
    const ibtn = (act, icon, enabled, title, cls) => `<button class="btn small icon-only ${cls || ''}" data-act="${act}" ${enabled === false ? 'disabled' : ''} title="${esc(title)}" aria-label="${esc(title)}">${icon}</button>`;
    const svg = id => `<svg class="ic sm"><use href="#${id}"/></svg>`;
    const row = (label, body, title) => `<div class="brow"${title ? ` title="${esc(title)}"` : ''}><span class="blbl">${label}</span><span class="bctl">${body}</span></div>`;
    function selDyn(d) {
      const f = app.curFrame(d), s = d.selection, m = app.measurement(d);
      let txt;
      if (!s.length) txt = '<span class="dim" title="Click atoms (Navigate or Select tool) or press Ctrl/⌘+A">No selection</span>';
      else if (s.length <= 4) txt = esc(s.map(i => label(f, i)).join(' – ')) + (m ? `<b class="mval">${m.text}</b>` : '');
      else txt = `<b>${s.length}</b> atoms selected`;
      return `<div class="route bsel">${txt}</div>`;
    }
    function cleanBtnDyn(d) {
      const n = app.curFrame(d).Z.length, running = relax && relax.mode === 'optimize' && relax.d === d;
      return running ? btn('b-opt-stop', svg('i-stop') + 'Stop', true, 'Stop the clean-up (keeps the geometry reached)', 'primary')
        : btn('b-opt', 'Clean', n > 0 && n <= LIMIT.clean, n > LIMIT.clean ? 'Limited to ' + LIMIT.clean + ' atoms' : 'Relax the geometry with the clean-up force field (not a QM method)', 'primary');
    }
    function cleanInfoDyn(d) {
      const running = relax && relax.mode === 'optimize' && relax.d === d;
      if (running) return `<div class="route mono binfo">step ${relax.steps} · E ${relax.e.toFixed(3)} kcal/mol · rms grad ${relax.rms.toFixed(3)}</div>`;
      if (d._ff) return `<div class="route mono binfo" title="Clean-up force field energy (kcal/mol), not a quantum-chemical energy">E ${d._ff.e.toFixed(3)} (from ${d._ff.e0.toFixed(3)}) kcal/mol · ${d._ff.steps} step${d._ff.steps === 1 ? '' : 's'} · ${d._ff.stopped ? 'stopped' : d._ff.converged ? 'converged' : 'not converged'}</div>`;
      return '';
    }
    function freezeDyn(d) {
      const nf = d.frozen.length;
      return `${btn('b-freeze', 'Freeze sel.', d.selection.length > 0, 'Selected atoms keep their positions during clean-up')}${btn('b-unfreeze', 'Unfreeze', nf > 0, 'Unfreeze all atoms')}<span class="dim">${nf ? nf + ' frozen' : ''}</span>`;
    }
    function editDyn(d) {
      const f = app.curFrame(d), s = d.selection, has = s.length > 0, m = app.measurement(d);
      let bond = '';
      if (s.length === 2) {
        const o = bondOrder(s[0], s[1]);
        bond = `<span class="bseg">${[1, 2, 3].map(k => btn('b-bond-' + k, ['–', '=', '≡'][k - 1], o !== k, ['Single', 'Double', 'Triple'][k - 1] + ' bond ' + label(f, s[0]) + '–' + label(f, s[1]), o === k ? 'on' : '')).join('')}${btn('b-bond-0', '×', !!o, 'Delete the bond')}</span>`;
      } else bond = '<span class="dim" title="Select two atoms to set, create or delete their bond">2 atoms</span>';
      const unit = m ? { distance: 'Å', angle: '°', dihedral: '°' }[m.kind] : '';
      return row('Atom', `<input type="text" id="edit-el" placeholder="El" title="New element for the selected atoms, e.g. N, Cl, Pt" value="${esc(B.inputs.el)}" ${has ? '' : 'disabled'}>
          ${ibtn('b-pt', '⋯', has, 'Periodic table')}${btn('edit-element', 'Set', has, 'Change the element of the selected atoms')}
          <span class="bsep"></span>${btn('b-charge--1', 'q−', has, 'Lower the formal charge of the selected atoms')}${btn('b-charge-1', 'q+', has, 'Raise the formal charge of the selected atoms')}
          ${ibtn('edit-delete', '<svg class="ic sm"><use href="#i-close"/></svg>', has, 'Delete the selected atoms (Delete)')}`) +
        row('Bond', bond + `<span class="bsep"></span>${btn('b-perceive', 'Orders', f.Z.length > 1, 'Perceive single/double/triple/aromatic bond orders')}${btn('b-rebond', 'Rebond', f.Z.length > 1, 'Recompute the connectivity from distances')}`) +
        row(m ? { distance: 'Dist', angle: 'Angle', dihedral: 'Dihed' }[m.kind] : 'Value', `<input type="number" id="edit-val" step="any" value="${m ? m.value.toFixed(m.kind === 'distance' ? 4 : 2) : ''}" ${m ? '' : 'disabled'} title="Select 2, 3 or 4 atoms; the smaller fragment moves (inside rings only the last atom)"><span class="dim">${unit}</span>${btn('edit-set', 'Apply', !!m, 'Set the distance, angle or dihedral', 'primary')}`);
    }
    function alignDyn(d) {
      const s = d.selection;
      return `<div class="btnrow">${btn('b-align-origin', 'Atom → 0', s.length === 1, 'Translate so the selected atom is at the origin')}
          ${['X', 'Y', 'Z'].map(a => btn('b-align-' + a.toLowerCase(), '→' + a, s.length === 2, 'Bond along +' + a + ': first atom at the origin, second on the +' + a + ' axis')).join('')}
          ${btn('b-align-plane', '→XY', s.length === 3, 'Plane of three atoms into XY (first at the origin, second on +X)')}
          ${btn('edit-center', 'Centre', true, 'Move the centroid to the origin')}</div>`;
    }
    function symDyn(d) {
      const r = d._sym;
      if (!r) return '';
      const grp = pgHtml(r.group), ord = r.order === Infinity ? '∞' : r.order;
      if (r.kind === 'detect') return `<div class="route">Point group <b>${grp}</b> · order ${ord}${r.tol && Math.abs(r.tol - B.symTol) > 1e-9 ? ` · tolerance used ${r.tol.toFixed(3)} Å` : ''}</div>`;
      return `<div class="route">Symmetrized to <b>${grp}</b> · order ${ord} · max shift ${r.maxShift.toFixed(4)} Å</div>`;
    }
    function pgHtml(g) {
      const s = String(g || '').replace('inf', '∞');
      const m = /^([A-Z])(.*)$/.exec(s);
      return m ? esc(m[1]) + (m[2] ? '<sub>' + esc(m[2]) + '</sub>' : '') : esc(s);
    }
    function panelBuild(d) {
      const f = app.curFrame(d), n = f.Z.length, coordsOpen = !B.fold.coords;
      const lib = filterEntries(libraryEntries(['molecule', 'group', 'ligand']), B.inputs.libq);
      let coords = '';
      if (coordsOpen) {
        if (n > LIMIT.editor) coords = `<p class="note">Limited to ${LIMIT.editor} atoms.</p>`;
        else coords = `<textarea id="bp-coords" class="coords" spellcheck="false" wrap="off" title="One “Symbol x y z” line per atom (Å), or paste an XYZ file, Gaussian input or ORCA input">${esc(d._coordDirty ? d._coordText : kit.coordText(f.Z, f.xyz))}</textarea>
          <div class="btnrow" style="margin-top:6px">${btn('b-coords-apply', 'Apply', true, 'Replace the geometry with these coordinates (recomputes bonds and orders)', 'primary')}${btn('b-coords-copy', 'Copy')}${btn('b-coords-revert', 'Revert', !!d._coordDirty)}<span class="dim">${d._coordDirty ? 'modified' : n + ' atoms'}</span></div>`;
      }
      const editing = d.built ? '' : f.edited ? `<span class="dim" title="Edits apply to this copy; the file geometries are unchanged">editing “${esc(f.label || 'Edited')}”</span>` : '';
      return `<div class="sec bhead"><div data-bp="sel">${selDyn(d)}</div>
          <div class="bq">
            <span class="bgrp">${ibtn('undo', svg('i-undo'), d.undo.length > 0, 'Undo (Ctrl/⌘+Z)')}${ibtn('redo', svg('i-redo'), d.redo.length > 0, 'Redo (Ctrl/⌘+Shift+Z)')}</span>
            <span class="bgrp">${btn('b-h-add', '+H', n > 0, 'Add hydrogens to satisfy valences (selected atoms, or all)')}${btn('b-h-remove', '−H', n > 0, 'Remove hydrogens (selected atoms, or all)')}${btn('b-h-adjust', '±H', n > 0, 'Adjust hydrogens: remove wrong counts and re-add')}</span>
            <span class="bgrp" data-bp="clean">${cleanBtnDyn(d)}</span>
            <span class="bgrp">${btn('b-copy-smiles', 'SMILES', n > 0, 'Copy the SMILES of the structure')}</span>
          </div>
          <div data-bp="cinfo">${cleanInfoDyn(d)}</div>
          ${d._smiles ? `<div class="route mono smiles-out" title="SMILES of the current structure">${esc(d._smiles)}</div>` : ''}${editing}</div>` +
        sec('insert', 'Insert', `
          <div class="inline"><input type="text" id="bp-smiles" class="grow" placeholder="SMILES, e.g. c1ccccc1O" spellcheck="false" value="${esc(B.inputs.smiles)}" title="Build 3D coordinates from SMILES (Enter)">${btn('b-smiles', 'Insert', true, 'Build 3D coordinates from SMILES')}</div>
          <input type="search" id="bp-lib-q" class="libq" placeholder="Library: molecules, groups, ligands…" spellcheck="false" value="${esc(B.inputs.libq)}" title="Molecules are inserted beside the structure; a group or ligand goes onto the selected atom(s), otherwise the Template tool is armed">
          <div class="lib" id="bp-lib">${entryRows(lib, null)}</div>`) +
        sec('edit', 'Edit', `<div data-bp="edit">${editDyn(d)}</div>`) +
        sec('align', 'Align', `<div data-bp="align">${alignDyn(d)}</div>`) +
        sec('clean', 'Clean-up', `<div class="btnrow"><label class="inline chk" title="Relax after each drawing, template or manipulation step (up to ${LIMIT.auto} atoms)"><input type="checkbox" data-bopt="autoOpt" ${B.autoOpt ? 'checked' : ''}> Auto</label><span data-bp="freeze" class="inline">${freezeDyn(d)}</span></div>`) +
        sec('sym', 'Symmetry', `<div class="inline"><select id="bp-symtol" title="Tolerance for symmetry detection">${SYM_TOLS.map(([v, t]) => `<option value="${v}" ${Math.abs(v - B.symTol) < 1e-9 ? 'selected' : ''}>${t}</option>`).join('')}</select>
            ${btn('b-sym-detect', 'Detect', n > 0 && n <= LIMIT.sym, 'Detect the point group')}${btn('b-sym-apply', 'Symmetrize', n > 1 && n <= LIMIT.sym, 'Make the structure exactly symmetric (new edited geometry)')}</div>
          <div data-bp="sym" style="margin-top:6px">${symDyn(d)}</div>${n > LIMIT.sym ? `<p class="note">Limited to ${LIMIT.sym} atoms.</p>` : ''}`) +
        sec('coords', 'Coordinates', coords);
    }
    function refreshSelection(d) {
      if (!d) return;
      if (B.tool === 'template' && B.tpl.mode === 'ligands') renderOptions();
      const p = $('#panel');
      if (!p || d.tab !== 'build') return;
      const parts = { sel: selDyn, edit: editDyn, align: alignDyn, freeze: freezeDyn };
      for (const k in parts) { const box = $(`[data-bp="${k}"]`, p); if (box) box.innerHTML = parts[k](d); }
    }

    // ---- Build tab actions
    function hTargets() { const d = doc(); return d.selection.length ? { atoms: d.selection.slice() } : {}; }
    function warnUnperceived() {
      const d = doc(), f = cur();
      if (!f.edited && f.Z.length > LIMIT.perceive) toast('Bond orders are not perceived above ' + LIMIT.perceive + ' atoms, so hydrogens assume single bonds.');
      return d;
    }
    function action(act, elm) {
      const d = doc();
      if (act === 'b-new') { newMolecule(); return true; }
      if (act === 'b-fold') { const id = elm.dataset.sec; B.fold[id] = !B.fold[id]; app.savePrefs(); if (id === 'coords') app.renderPanel(); else { const s = elm.closest('.bsec'); if (s) s.classList.toggle('folded', !!B.fold[id]); } return true; }
      if (act === 'b-lib') { const [kind, ...rest] = elm.dataset.lib.split(':'); useLibrary(kind, rest.join(':')); return true; }
      if (act === 'b-smiles') { const inp = $('#bp-smiles'); insertSmiles(inp ? inp.value : B.inputs.smiles); return true; }
      if (!d) return false;
      const f = cur();
      switch (act) {
        case 'b-copy-smiles': {
          let s;
          try { s = currentSmiles(); } catch (e) { toast('SMILES: ' + e.message, true); return true; }
          d._smiles = s; app.renderPanel();
          copyText(s, 'SMILES copied to the clipboard');
          return true;
        }
        case 'b-h-add': {
          let none = false;
          warnUnperceived();
          edit(mol => { const r = BM.chem.addHydrogens(mol, hTargets()); if (!r.added.length) { none = true; return null; } return { mol: r.mol, map: null, n: r.added.length }; }, r => 'Added ' + r.n + ' hydrogen' + (r.n > 1 ? 's' : ''));
          if (none) status('No hydrogens were missing.');
          return true;
        }
        case 'b-h-remove': {
          let none = false;
          edit(mol => { const r = BM.chem.removeHydrogens(mol, hTargets()); const nn = mol.Z.length - r.mol.Z.length; if (!nn) { none = true; return null; } return { mol: r.mol, map: r.map, n: nn }; }, r => 'Removed ' + r.n + ' hydrogen' + (r.n > 1 ? 's' : ''));
          if (none) status('There were no hydrogens to remove.');
          return true;
        }
        case 'b-h-adjust': warnUnperceived(); edit(mol => { const r = BM.chem.adjustHydrogens(mol, hTargets()); return { mol: r.mol, map: r.map, n: r.mol.Z.filter(z => z === 1).length - mol.Z.filter(z => z === 1).length }; }, r => 'Hydrogens adjusted (' + (r.n >= 0 ? '+' : '−') + Math.abs(r.n) + ')'); return true;
        case 'b-perceive': app.doEdit(ff => { const p = perceive(ff.Z, ff.xyz, ff.bonds, ff.formalCharges); ff.bonds = p.bonds; ff.formalCharges = p.charges; }, 'Bond orders perceived'); return true;
        case 'b-rebond': app.doEdit(ff => { const pairs = geom.computeBonds(ff.Z, ff.xyz, { tolerance: S.view.tol }); const p = ff.Z.length <= LIMIT.perceive * 10 ? perceive(ff.Z, ff.xyz, pairs, ff.formalCharges) : { bonds: pairsToList(pairs), charges: ff.formalCharges }; ff.bonds = p.bonds; ff.formalCharges = p.charges; }, 'Bonds recomputed from distances'); return true;
        case 'b-bond-0': case 'b-bond-1': case 'b-bond-2': case 'b-bond-3': {
          const s = d.selection; if (s.length !== 2) return true;
          const o = +act.slice(-1), names = label(f, s[0]) + '–' + label(f, s[1]);
          edit(mol => ({ mol: BM.chem.setBond(mol, s[0], s[1], o), map: null }), o ? names + ': ' + orderName(o) + ' bond' : 'Deleted bond ' + names);
          return true;
        }
        case 'b-pt': { const inp = $('#edit-el'); openPeriodic(elm, Math.max(1, el.toZ(inp && inp.value) > 0 ? el.toZ(inp.value) : 6), z => { B.inputs.el = el.symbol(z); const i2 = $('#edit-el'); if (i2) i2.value = B.inputs.el; }); return true; }
        case 'b-charge--1': case 'b-charge-1': {
          const dq = act === 'b-charge-1' ? 1 : -1, s = d.selection.slice(); if (!s.length) return true;
          app.doEdit(ff => { if (!ff.formalCharges || ff.formalCharges.length !== ff.Z.length) ff.formalCharges = new Array(ff.Z.length).fill(0); for (const i of s) ff.formalCharges[i] = (ff.formalCharges[i] || 0) + dq; },
            s.length === 1 ? 'Formal charge of ' + label(f, s[0]) + ': ' + chargeText((f.formalCharges ? f.formalCharges[s[0]] || 0 : 0) + dq) : 'Formal charges changed by ' + chargeText(dq));
          return true;
        }
        case 'b-align-origin': case 'b-align-x': case 'b-align-y': case 'b-align-z': case 'b-align-plane': align(act.slice(8)); return true;
        case 'b-opt': optimizeGeometry(); return true;
        case 'b-opt-stop': relaxStop(false); app.renderPanel(); return true;
        case 'b-freeze': d.frozen = uniq(d.frozen.concat(d.selection)); R().setAtomMarkers(d.frozen); refreshCleanup(); status(d.frozen.length + ' atoms frozen for clean-up'); return true;
        case 'b-unfreeze': d.frozen = []; R().setAtomMarkers([]); refreshCleanup(); status('All atoms unfrozen'); return true;
        case 'b-sym-detect': symDetect(); return true;
        case 'b-sym-apply': symApply(); return true;
        case 'b-coords-apply': coordsApply(); return true;
        case 'b-coords-copy': { const t = $('#bp-coords'); copyText(t ? t.value : kit.coordText(f.Z, f.xyz), 'Coordinates copied to the clipboard'); return true; }
        case 'b-coords-revert': d._coordDirty = false; d._coordText = ''; app.renderPanel(); return true;
      }
      return false;
    }
    function copyText(text, okMsg) { app.copyText(text, okMsg); }
    function align(kind) {
      const d = doc(), s = d.selection.slice();
      const need = { origin: 1, x: 2, y: 2, z: 2, plane: 3 }[kind];
      if (s.length !== need) { toast('Select ' + need + ' atom' + (need > 1 ? 's' : '') + ' first.', true); return; }
      app.doEdit(f => {
        const xyz = Float64Array.from(f.xyz), n = f.Z.length, o = P(xyz, s[0]);
        for (let i = 0; i < n; i++) for (let k = 0; k < 3; k++) xyz[3 * i + k] -= o[k];
        if (kind !== 'origin') {
          let M;
          const u = vunit(P(xyz, s[1]));
          if (vlen(u) < 1e-9) throw new Error('The selected atoms coincide.');
          if (kind === 'plane') {
            const w = vunit(vcross(u, P(xyz, s[2])));
            if (vlen(w) < 1e-9) throw new Error('The three atoms are collinear, so they do not define a plane.');
            M = [u, vcross(w, u), w];
          } else {
            const e = { x: [1, 0, 0], y: [0, 1, 0], z: [0, 0, 1] }[kind];
            let ax = vcross(u, e);
            const c = Math.max(-1, Math.min(1, vdot(u, e)));
            if (vlen(ax) < 1e-9) ax = c > 0 ? null : vunit(vcross(u, Math.abs(u[0]) < 0.9 ? [1, 0, 0] : [0, 1, 0]));
            if (ax) geom.rotateAtoms(xyz, Array.from({ length: n }, (_, i) => i), [0, 0, 0], ax, Math.acos(c));
          }
          if (M) for (let i = 0; i < n; i++) { const p = P(xyz, i); for (let k = 0; k < 3; k++) xyz[3 * i + k] = vdot(M[k], p); }
        }
        f.xyz = xyz;
      }, { origin: label(cur(), s[0]) + ' moved to the origin', x: 'Bond aligned with X', y: 'Bond aligned with Y', z: 'Bond aligned with Z', plane: 'Plane aligned with XY' }[kind]);
      R().fit();
    }
    function symDetect() {
      const d = doc(), f = cur();
      if (f.Z.length > LIMIT.sym) { toast('Symmetry tools are limited to ' + LIMIT.sym + ' atoms.', true); return; }
      try {
        const r = BM.sym.detect(f.Z, f.xyz, { tol: B.symTol });
        d._sym = { kind: 'detect', group: r.pointGroup, order: r.order, tol: r.tol };
        status('Point group ' + r.pointGroup + ' (order ' + (r.order === Infinity ? '∞' : r.order) + ')');
      } catch (e) { toast('Symmetry: ' + e.message, true); return; }
      const box = $('#panel [data-bp="sym"]'); if (box) box.innerHTML = symDyn(d);
    }
    function symApply() {
      const d = doc(), f0 = cur();
      if (f0.Z.length > LIMIT.sym) { toast('Symmetry tools are limited to ' + LIMIT.sym + ' atoms.', true); return; }
      let info = null;
      const ok = app.doEdit(f => { info = BM.sym.symmetrize(f.Z, f.xyz, { tol: B.symTol }); f.xyz = Float64Array.from(info.xyz); },
        () => 'Symmetrized to ' + info.pointGroup + '; max shift ' + info.maxShift.toFixed(4) + ' Å');
      if (ok && info) { d._sym = { kind: 'sym', group: info.pointGroup, order: info.order, maxShift: info.maxShift }; const box = $('#panel [data-bp="sym"]'); if (box) box.innerHTML = symDyn(d); }
    }
    function coordsApply() {
      const d = doc(), t = $('#bp-coords');
      const text = t ? t.value : d._coordText;
      let r;
      try { r = kit.parseCoords(text); } catch (e) { toast(e.message, true); return; }
      const ok = app.doEdit(f => {
        const same = r.Z.length === f.Z.length && r.Z.every((z, i) => z === f.Z[i]);
        if (!same) { f.atomNames = null; f.resNames = null; f.resIds = null; d.selection = []; d.frozen = []; }
        f.Z = r.Z.slice(); f.xyz = Float64Array.from(r.xyz);
        const pairs = geom.computeBonds(f.Z, f.xyz, { tolerance: S.view.tol });
        const p = f.Z.length <= LIMIT.perceive ? perceive(f.Z, f.xyz, pairs, null) : { bonds: pairsToList(pairs), charges: new Array(f.Z.length).fill(0) };
        f.bonds = p.bonds; f.formalCharges = p.charges;
      }, 'Applied ' + r.Z.length + ' atoms from the ' + r.source);
      if (ok) { d._coordDirty = false; d._coordText = ''; app.renderPanel(); R().fit(); }
    }
    function panelInput(e) {
      const t = e.target, d = doc();
      if (t.id === 'bp-smiles') { B.inputs.smiles = t.value; return true; }
      if (t.id === 'edit-el') { B.inputs.el = t.value; return true; }
      if (t.id === 'bp-lib-q') { B.inputs.libq = t.value; const box = $('#bp-lib'); if (box) box.innerHTML = entryRows(filterEntries(libraryEntries(['molecule', 'group', 'ligand']), t.value), null); return true; }
      if (t.id === 'bp-coords' && d) { d._coordText = t.value; if (!d._coordDirty) { d._coordDirty = true; const span = t.parentElement.querySelector('.btnrow .dim'); if (span) span.textContent = 'modified'; const rv = t.parentElement.querySelector('[data-act="b-coords-revert"]'); if (rv) rv.disabled = false; } return true; }
      if (t.id === 'bp-symtol') { B.symTol = parseFloat(t.value) || 0.1; app.savePrefs(); return true; }
      if (t.dataset.bopt === 'autoOpt') { B.autoOpt = t.checked; app.savePrefs(); if (!t.checked && relax && relax.mode === 'auto') relaxStop(true); renderOptions(); return true; }
      return false;
    }
    function panelKey(e) {
      const t = e.target;
      if (e.key === 'Enter' && t.id === 'bp-smiles') { e.preventDefault(); insertSmiles(t.value); return true; }
      if (e.key === 'Enter' && t.id === 'edit-el') { e.preventDefault(); app.editElement(); return true; }
      if (e.key === 'Enter' && t.id === 'edit-val') { e.preventDefault(); app.editSet(); return true; }
      if (e.key === 'Enter' && t.id === 'bp-lib-q') { e.preventDefault(); const b = $('#bp-lib .lib-item'); if (b) b.click(); return true; }
      return false;
    }

    // ---------------------------------------------------------------- dispatch from app.js
    function atomClick(i, e) {
      if (!doc()) return false;
      if (B.tool === 'select') { const d = doc(); setSelection(combineSel(d.selection, [i], e && e.shiftKey ? 'add' : e && (e.ctrlKey || e.metaKey) ? 'remove' : 'toggle')); return true; }
      if (B.tool === 'template') { templateAtomClick(i, e); return true; }
      if (B.tool === 'draw') return true;
      return false;
    }
    function backgroundClick(e) {
      if (!doc()) return false;
      if (B.tool === 'draw') { if (!e || !e.altKey) { const [px, py] = e ? localXY(e) : [R().canvas.clientWidth / 2, R().canvas.clientHeight / 2]; drawAddAt(px, py); } return true; }
      if (B.tool === 'template') { templateBackgroundClick(e); return true; }
      return false;
    }
    function atomDblClick(i, e) {
      if (!doc()) return false;
      if (B.tool === 'select') { if (i >= 0) { const frag = fragmentOf(i); setSelection(combineSel(doc().selection, frag, 'add')); status(frag.length + '-atom fragment selected'); } return true; }
      if (B.tool === 'draw' || B.tool === 'template') { if (i < 0) backgroundClick(e); else if (B.tool === 'template') templateAtomClick(i, e); return true; }
      return false;
    }
    /** After every undoable edit (app.doEdit). */
    function afterEdit(d, opts) {
      opts = opts || {};
      d._sym = null; d._ff = null; d._smiles = null;
      if (opts.fit) fitForBuilding();
      if (opts.auto && B.autoOpt && B.tool !== 'navigate') relaxStart('auto');
      showHint();
      if (B.tool === 'select') renderOptions();
    }
    function fitForBuilding() {
      const r = R();
      r.fit();
      const minR = 6;
      if (r.radius < minR) {
        const half = (r.opt.fov / 2) * Math.PI / 180, aspect = Math.max(0.2, r.canvas.clientWidth / Math.max(1, r.canvas.clientHeight)), fr = minR * 1.05;
        r.dist = Math.max(fr / Math.sin(half), fr / Math.sin(Math.atan(Math.tan(half) * aspect)));
        r.requestRender();
      }
    }
    function onActivate() {
      if (relax) relaxStop(true);
      drag = null; chip = null; closePopup();
      const r = R(); if (r) r.highlightBond(-1);
      renderToolbox(); renderOptions(); showHint();
    }
    function overlay(ctx, project) {
      const g = drag;
      if (!g) return;
      if (g.kind === 'band' && g.moved) {
        const x = Math.min(g.x0, g.x1), y = Math.min(g.y0, g.y1), w = Math.abs(g.x1 - g.x0), h = Math.abs(g.y1 - g.y0);
        ctx.fillStyle = 'rgba(79,140,255,0.12)'; ctx.fillRect(x, y, w, h);
        ctx.setLineDash([5, 4]); ctx.lineWidth = 1.2; ctx.strokeStyle = 'rgba(138,180,255,0.95)'; ctx.strokeRect(x + 0.5, y + 0.5, w, h); ctx.setLineDash([]);
      } else if (g.kind === 'draw' && g.moved && g.start >= 0) {
        const xyz = R().xyz, a = project(P(xyz, g.start)), b = g.target >= 0 ? project(P(xyz, g.target)) : { x: g.px, y: g.py };
        ctx.setLineDash([6, 4]); ctx.lineWidth = 2; ctx.strokeStyle = 'rgba(255,212,94,0.95)';
        ctx.beginPath(); ctx.moveTo(a.x, a.y); ctx.lineTo(b.x, b.y); ctx.stroke(); ctx.setLineDash([]);
        ctx.beginPath(); ctx.arc(b.x, b.y, g.target >= 0 ? 15 : 11, 0, 2 * Math.PI);
        if (g.target >= 0) { ctx.lineWidth = 2; ctx.strokeStyle = 'rgba(255,212,94,0.95)'; ctx.stroke(); return; }
        ctx.fillStyle = 'rgba(10,16,26,0.8)'; ctx.fill(); ctx.lineWidth = 1.6; ctx.strokeStyle = rgb(el.color[B.element] || [1, 1, 1]); ctx.stroke();
        ctx.fillStyle = '#f3f6fb'; ctx.font = '600 11px ui-sans-serif, -apple-system, "Segoe UI", Roboto, sans-serif'; ctx.textAlign = 'center'; ctx.textBaseline = 'middle';
        ctx.fillText(el.symbol(B.element), b.x, b.y + 0.5);
      }
    }
    function keydown(e) {
      if (e.key === 'Escape') {
        if (popup) { closePopup(); return true; }
        if (relax) { relaxStop(false); if (doc() && doc().tab === 'build') app.renderPanel(); return true; }
      }
      return false;
    }

    bindOptions();
    renderToolbox();
    return {
      setTool, panelBuild, refreshSelection, action, panelInput, panelKey, atomClick, backgroundClick, atomDblClick,
      afterEdit, onActivate, overlay, keydown, insertSmiles, newMolecule, useLibrary, relaxStop, optimizeGeometry,
      showHint, renderOptions, closePopup, chipHtml: () => chip, running: () => (relax ? relax.mode : null),
      tools: T, selectAction, currentSmiles, fitForBuilding, hintText,
      cancelDrag: () => { drag = null; chip = null; }
    };
  }

  BM.buildUI = { create, kit, TOOLS };
})(typeof window !== 'undefined' ? (window.BM = window.BM || {}) : (globalThis.BM = globalThis.BM || {}));
