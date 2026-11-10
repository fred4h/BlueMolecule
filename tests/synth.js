/* Synthetic test structures for the representation tests (Node: require('./synth.js'); browser: add as a script).
 * Ideal backbones from phi/psi (NeRF placement); a helix-bundle + antiparallel-sheet protein with arc-length loops and a
 * ligand pocket; a POPC bilayer with hydrogens, TIP3 water and SOD/CLA; a heparan-sulfate chain on a serine (sugars
 * from BM.chem.smilesTo3D, linked by superposing the anomeric O1 onto the acceptor O). Clashing side chains are
 * trimmed so that distance-based bonding finds no spurious bonds.
 * Writers: toPDB(sys, {ss}) with HELIX/SHEET records, toGRO(sys), toDataset(sys).
 * sys = {Z, xyz: number[], names, res, ids, chains, bonds: [[i, j, order]], ss?: [{type, chain, start, end}]} */
(function (root) {
  'use strict';
  const D = Math.PI / 180;
  const sub = (a, b) => [a[0] - b[0], a[1] - b[1], a[2] - b[2]], add = (a, b) => [a[0] + b[0], a[1] + b[1], a[2] + b[2]];
  const mul = (a, s) => [a[0] * s, a[1] * s, a[2] * s], dot = (a, b) => a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
  const cross = (a, b) => [a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]];
  const len = a => Math.hypot(a[0], a[1], a[2]), norm = a => mul(a, 1 / (len(a) || 1));
  const perp = a => norm(Math.abs(a[0]) < 0.9 ? cross(a, [1, 0, 0]) : cross(a, [0, 1, 0]));
  /** D with |CD| = r, angle BCD = theta (deg), dihedral ABCD = phi (deg). */
  function place(A, B, C, r, theta, phi) {
    const bc = norm(sub(C, B)), n = norm(cross(sub(B, A), bc)), m = cross(n, bc), t = theta * D, p = phi * D;
    return add(C, add(add(mul(bc, -r * Math.cos(t)), mul(m, r * Math.sin(t) * Math.cos(p))), mul(n, r * Math.sin(t) * Math.sin(p))));
  }
  function dihedral(a, b, c, d) {
    const b0 = sub(a, b), b1 = norm(sub(c, b)), b2 = sub(d, c);
    const v = sub(b0, mul(b1, dot(b0, b1))), w = sub(b2, mul(b1, dot(b2, b1)));
    return Math.atan2(dot(cross(b1, v), w), dot(v, w)) / D;
  }
  /** Backbone N, CA, C, O of residues with the given phi/psi (deg), trans peptides. */
  function backbone(phis, psis) {
    const out = [];
    let N = [0, 0, 0], CA = [1.458, 0, 0], C = add(CA, [1.525 * Math.cos(68.8 * D), 1.525 * Math.sin(68.8 * D), 0]);
    for (let k = 0; k < phis.length; k++) {
      if (k > 0) {
        const p = out[k - 1];
        N = place(p.N, p.CA, p.C, 1.329, 116.2, psis[k - 1]);
        CA = place(p.CA, p.C, N, 1.458, 121.7, 180);
        C = place(p.C, N, CA, 1.525, 111.2, phis[k]);
      }
      out.push({ N, CA, C });
    }
    for (let k = 0; k < out.length; k++) out[k].O = place(out[k].N, out[k].CA, out[k].C, 1.231, 120.5, psis[k] + 180);
    return out;
  }
  // side chains: [name, Z, parentName (default: previous atom), dihedral offset]
  const SIDE = {
    ALA: [], SER: [['OG', 8]], CYS: [['SG', 16]], VAL: [['CG1', 6], ['CG2', 6, 'CB', 120]], THR: [['OG1', 8], ['CG2', 6, 'CB', 120]],
    LEU: [['CG', 6], ['CD1', 6], ['CD2', 6, 'CG', 120]], ILE: [['CG1', 6], ['CG2', 6, 'CB', 120], ['CD1', 6, 'CG1']], MET: [['CG', 6], ['SD', 16], ['CE', 6]],
    LYS: [['CG', 6], ['CD', 6], ['CE', 6], ['NZ', 7]], GLU: [['CG', 6], ['CD', 6], ['OE1', 8], ['OE2', 8, 'CD', 120]], ASP: [['CG', 6], ['OD1', 8], ['OD2', 8, 'CG', 120]],
    ASN: [['CG', 6], ['OD1', 8], ['ND2', 7, 'CG', 120]], GLN: [['CG', 6], ['CD', 6], ['OE1', 8], ['NE2', 7, 'CD', 120]],
    ARG: [['CG', 6], ['CD', 6], ['NE', 7], ['CZ', 6], ['NH1', 7], ['NH2', 7, 'CZ', 120]], GLY: null
  };
  function residueAtoms(bb, name) {
    const at = [['N', 7, bb.N], ['CA', 6, bb.CA], ['C', 6, bb.C], ['O', 8, bb.O]];
    const sc = SIDE[name];
    if (sc === null || bb.loop) return at;
    const CB = place(bb.C, bb.N, bb.CA, 1.53, 110.5, -122.6);
    at.push(['CB', 6, CB]);
    const P = { N: bb.N, CA: bb.CA, CB }, parent = { CA: 'N', CB: 'CA' };
    let last = 'CB';
    for (const [nm, z, par, off] of sc || []) {
      const p = par || last, pp = parent[p], ppp = parent[pp];
      const pos = place(P[ppp], P[pp], P[p], z === 6 ? 1.52 : z === 16 ? 1.81 : 1.43, 112, 180 + (off || 0));
      P[nm] = pos; parent[nm] = p; at.push([nm, z, pos]);
      if (!par) last = nm;
    }
    return at;
  }
  function newSys() { return { Z: [], xyz: [], names: [], res: [], ids: [], chains: [], bonds: [] }; }
  function push(sys, name, z, p, res, id, chain) { sys.Z.push(z); sys.xyz.push(p[0], p[1], p[2]); sys.names.push(name); sys.res.push(res); sys.ids.push(id); sys.chains.push(chain); return sys.Z.length - 1; }
  const pos = (sys, i) => [sys.xyz[3 * i], sys.xyz[3 * i + 1], sys.xyz[3 * i + 2]];
  const xf = (R, t) => p => [R[0] * p[0] + R[1] * p[1] + R[2] * p[2] + t[0], R[3] * p[0] + R[4] * p[1] + R[5] * p[2] + t[1], R[6] * p[0] + R[7] * p[1] + R[8] * p[2] + t[2]];
  function rotTo(a, b) {
    a = norm(a); b = norm(b);
    const v = cross(a, b), c = dot(a, b);
    if (c < -0.9999) { const p = perp(a); return [2 * p[0] * p[0] - 1, 2 * p[0] * p[1], 2 * p[0] * p[2], 2 * p[1] * p[0], 2 * p[1] * p[1] - 1, 2 * p[1] * p[2], 2 * p[2] * p[0], 2 * p[2] * p[1], 2 * p[2] * p[2] - 1]; }
    const k = 1 / (1 + c);
    return [v[0] * v[0] * k + c, v[0] * v[1] * k - v[2], v[0] * v[2] * k + v[1], v[1] * v[0] * k + v[2], v[1] * v[1] * k + c, v[1] * v[2] * k - v[0], v[2] * v[0] * k - v[1], v[2] * v[1] * k + v[0], v[2] * v[2] * k + c];
  }
  function rotAxis(ax, ang) {
    const [x, y, z] = norm(ax), c = Math.cos(ang), s = Math.sin(ang), t = 1 - c;
    return [t * x * x + c, t * x * y - s * z, t * x * z + s * y, t * x * y + s * z, t * y * y + c, t * y * z - s * x, t * x * z - s * y, t * y * z + s * x, t * z * z + c];
  }
  const matMul = (A, B) => [0, 1, 2].flatMap(r => [0, 1, 2].map(c => A[3 * r] * B[c] + A[3 * r + 1] * B[3 + c] + A[3 * r + 2] * B[6 + c]));
  const centroid = pts => mul(pts.reduce((s, p) => add(s, p), [0, 0, 0]), 1 / pts.length);
  function axisOf(pts) {
    const c = centroid(pts), M = [0, 0, 0, 0, 0, 0, 0, 0, 0];
    for (const p of pts) { const v = sub(p, c); for (let a = 0; a < 3; a++) for (let b = 0; b < 3; b++) M[3 * a + b] += v[a] * v[b]; }
    let d = norm(sub(pts[pts.length - 1], pts[0]));
    for (let i = 0; i < 60; i++) d = norm([M[0] * d[0] + M[1] * d[1] + M[2] * d[2], M[3] * d[0] + M[4] * d[1] + M[5] * d[2], M[6] * d[0] + M[7] * d[1] + M[8] * d[2]]);
    if (dot(d, sub(pts[pts.length - 1], pts[0])) < 0) d = mul(d, -1);
    return { c, d };
  }
  const moveBB = (bb, T) => bb.map(r => Object.assign({}, r, { N: T(r.N), CA: T(r.CA), C: T(r.C), O: T(r.O) }));
  /** Kabsch–Sander energy (the formula of BM.cartoon.hbEnergy). */
  function hbE(c, o, n, h) { const d = (a, b) => Math.max(0.5, len(sub(a, b))); return 27.888 * (1 / d(o, n) + 1 / d(c, h) - 1 / d(o, h) - 1 / d(c, n)); }
  const H_of = (bb, k) => (k > 0 ? add(bb[k].N, norm(sub(bb[k - 1].C, bb[k - 1].O))) : null);

  // ---------------------------------------------------------------- secondary-structure elements
  const idealHelix = n => backbone(Array(n).fill(-57), Array(n).fill(-47));
  const idealStrand = n => backbone(Array(n).fill(-139), Array(n).fill(135));
  /** Helix with its axis along dir, centred at c. */
  function helixAt(n, c, dir, spin) {
    const h = idealHelix(n), ax = axisOf(h.map(r => r.CA)), R = matMul(rotAxis(dir, spin || 0), rotTo(ax.d, dir)), cc = xf(R, [0, 0, 0])(ax.c);
    return moveBB(h, xf(R, sub(c, cc)));
  }
  /** Two antiparallel strands in the xy plane (strand 1 along +x), strand 2 placed to maximise backbone H-bonds. */
  function antiparallelPair(n) {
    const A = idealStrand(n), ax = axisOf(A.map(r => r.CA));
    const R = rotTo(ax.d, [1, 0, 0]);
    let a = moveBB(A, xf(R, mul(xf(R, [0, 0, 0])(ax.c), -1)));
    const ov = norm(a.reduce((s, r, k) => add(s, mul(sub(r.O, r.C), k % 2 ? -1 : 1)), [0, 0, 0]));
    a = moveBB(a, xf(rotAxis([1, 0, 0], -Math.atan2(ov[2], ov[1])), [0, 0, 0]));
    const Rz = rotAxis([0, 0, 1], Math.PI);
    let best = null;
    for (let dy = 4.2; dy <= 5.6; dy += 0.1) for (let dx = -3.6; dx <= 3.6; dx += 0.15) for (const sy of [1, -1]) {
      const b = moveBB(a, xf(Rz, [dx, sy * dy, 0]));
      let e = 0, clash = 0;
      for (let i = 0; i < n; i++) for (let j = 0; j < n; j++) {
        const hj = H_of(b, j), hi = H_of(a, i);
        if (hj) { const v = hbE(a[i].C, a[i].O, b[j].N, hj); if (v < -0.5) e += v; }
        if (hi) { const v = hbE(b[j].C, b[j].O, a[i].N, hi); if (v < -0.5) e += v; }
        if (len(sub(a[i].CA, b[j].CA)) < 3.8) clash++;
      }
      if (!clash && (!best || e < best.e)) best = { e, b };
    }
    return [a, best.b];
  }
  /** Loop residues between two elements: CAs evenly spaced (≈3.75 Å) on a quadratic curve whose bulge direction is
   *  chosen to keep away from `avoid` (CA positions already placed); rough but bonded backbone (the flanking residues'
   *  C/O and N are adjusted to join it). */
  function loop(from, to, avoid) {
    const mid = mul(add(from.CA, to.CA), 0.5), axis = norm(sub(to.CA, from.CA)), d0 = len(sub(to.CA, from.CA));
    const u = perp(axis), v = cross(axis, u);
    let best = null;
    for (let q = 0; q < 16; q++) for (const mag of [3, 5, 7.5, 10]) {
      const dir = add(mul(u, Math.cos(q * Math.PI / 8)), mul(v, Math.sin(q * Math.PI / 8))), bulge = add(mid, mul(dir, mag + 0.2 * d0));
      const curve = t => add(add(mul(from.CA, (1 - t) * (1 - t)), mul(bulge, 2 * t * (1 - t))), mul(to.CA, t * t));
      const S = [], cum = [0];
      for (let k = 0; k <= 300; k++) S.push(curve(k / 300));
      for (let k = 1; k <= 300; k++) cum.push(cum[k - 1] + len(sub(S[k], S[k - 1])));
      const L = cum[300], nres = Math.max(1, Math.round(L / 3.75) - 1), step = L / (nres + 1), cas = [];
      for (let k = 1; k <= nres; k++) { const target = k * step; let w = 1; while (w < 300 && cum[w] < target) w++; cas.push(S[w]); }
      const chain = [from.CA].concat(cas, [to.CA]);
      if (chain.some((c, k) => k && (len(sub(c, chain[k - 1])) < 3.4 || len(sub(c, chain[k - 1])) > 4.2))) continue;
      let mind = Infinity;
      for (const c of cas) for (const a of avoid) { const dd = len(sub(c, a)); if (dd < mind) mind = dd; }
      const score = Math.min(mind, 7) - 0.05 * nres;
      if (!best || score > best.score) best = { score, cas };
    }
    const res = best.cas.map(ca => ({ CA: ca, loop: true }));
    const all = [from].concat(res, [to]), ups = [];
    for (let k = 1; k + 1 < all.length; k++) ups[k] = norm(cross(sub(all[k + 1].CA, all[k].CA), sub(all[k - 1].CA, all[k].CA)));
    ups[0] = ups[1] || perp(sub(all[1].CA, all[0].CA));
    for (let k = 1; k + 1 < all.length; k++) if (dot(ups[k], ups[k - 1]) < 0) ups[k] = mul(ups[k], -1);
    // loop residues get N, C, O; the flanking element residues keep their own backbone
    for (let k = 0; k + 1 < all.length; k++) {
      const a = all[k].CA, b = all[k + 1].CA, d = sub(b, a), up = ups[k];
      if (k > 0) { all[k].C = add(a, add(mul(d, 0.32), mul(up, 0.75))); all[k].O = add(all[k].C, mul(up, 1.23)); }
      if (k + 1 < all.length - 1) all[k + 1].N = add(a, add(mul(d, 0.68), mul(up, 0.75)));
    }
    for (const r of res) avoid.push(r.CA);
    return res;
  }

  // ---------------------------------------------------------------- systems
  function addChain(sys, bbs, seq, chain, id0) {
    bbs.forEach((bb, k) => { const nm = bb.loop ? 'GLY' : seq[k % seq.length]; for (const [an, z, p] of residueAtoms(bb, nm)) push(sys, an, z, p, nm, id0 + k, chain); });
  }
  /** Trim side chains (to CB, then to Gly) and drop waters that come closer than 2.5 Å to another residue. */
  function declash(sys) {
    const n = sys.Z.length, cell = 3, grid = new Map(), key = p => p.map(v => Math.floor(v / cell)).join(',');
    for (let i = 0; i < n; i++) { if (sys.Z[i] <= 1) continue; const k = key(pos(sys, i)); if (!grid.has(k)) grid.set(k, []); grid.get(k).push(i); }
    const resKey = i => sys.chains[i] + ':' + sys.ids[i], bbName = new Set(['N', 'CA', 'C', 'O']);
    const kill = new Set(), trimmed = new Map();
    for (let i = 0; i < n; i++) {
      if (sys.Z[i] <= 1) continue;
      const p = pos(sys, i), g = p.map(v => Math.floor(v / cell));
      for (let dx = -1; dx <= 1; dx++) for (let dy = -1; dy <= 1; dy++) for (let dz = -1; dz <= 1; dz++) {
        for (const j of grid.get([g[0] + dx, g[1] + dy, g[2] + dz].join(',')) || []) {
          if (j <= i || resKey(i) === resKey(j) || len(sub(p, pos(sys, j))) >= 2.5) continue;
          // peptide C-N (adjacent residues) and explicit bonds are fine
          if (sys.chains[i] === sys.chains[j] && Math.abs(sys.ids[i] - sys.ids[j]) === 1 && ((sys.names[i] === 'C' && sys.names[j] === 'N') || (sys.names[i] === 'N' && sys.names[j] === 'C'))) continue;
          if (sys.bonded && (sys.bonded.has(i + ':' + j) || sys.bonded.has(j + ':' + i))) continue;
          if (bbName.has(sys.names[i]) && bbName.has(sys.names[j])) {
            // two backbone atoms (rough loop junctions): drop a carbonyl O, else the loop atom
            const o = sys.names[i] === 'O' ? i : sys.names[j] === 'O' ? j : sys.res[i] === 'GLY' ? i : j;
            if (len(sub(p, pos(sys, j))) < 2.0) sys.Z[o] = -1;
            continue;
          }
          if (sys.res[i] === 'POPC' && sys.res[j] === 'POPC') continue; // tight lipid packing is fine
          for (const a of [i, j]) {
            if (sys.res[a] === 'SOL' || sys.res[a] === 'TIP3' || sys.res[a] === 'POPC') { kill.add(resKey(a)); continue; }
            if (!AA.test(sys.res[a]) || bbName.has(sys.names[a]) || (sys.fixed && sys.fixed.has(a))) continue;
            trimmed.set(resKey(a), sys.names[a] === 'CB' ? 'all' : trimmed.get(resKey(a)) === 'all' ? 'all' : 'cb');
          }
        }
      }
    }
    for (let i = 0; i < n; i++) {
      const k = resKey(i), t = trimmed.get(k);
      if (sys.Z[i] < 0) continue;
      if (kill.has(k) || (t === 'all' && !bbName.has(sys.names[i])) || (t === 'cb' && !bbName.has(sys.names[i]) && sys.names[i] !== 'CB')) sys.Z[i] = -1;
      else if (t === 'all') sys.res[i] = 'GLY';
      else if (t === 'cb' && !(sys.fixed && sys.fixed.has(i))) sys.res[i] = 'ALA';
    }
    compact(sys);
  }
  function compact(sys) {
    const keep = [], remap = new Map();
    for (let i = 0; i < sys.Z.length; i++) if (sys.Z[i] >= 0) { remap.set(i, keep.length); keep.push(i); }
    if (keep.length === sys.Z.length) return remap;
    const xyz = [];
    for (const i of keep) xyz.push(sys.xyz[3 * i], sys.xyz[3 * i + 1], sys.xyz[3 * i + 2]);
    sys.xyz = xyz;
    for (const k of ['Z', 'names', 'res', 'ids', 'chains']) sys[k] = keep.map(i => sys[k][i]);
    sys.bonds = sys.bonds.filter(b => remap.has(b[0]) && remap.has(b[1])).map(b => [remap.get(b[0]), remap.get(b[1]), b[2]]);
    if (sys.fixed) sys.fixed = new Set([...sys.fixed].filter(i => remap.has(i)).map(i => remap.get(i)));
    return remap;
  }
  /** Three-helix bundle + two-stranded antiparallel sheet in one chain A, ligand LIG in the bundle, waters and ions. */
  function protein(BM, opts) {
    opts = opts || {};
    const sys = newSys(), H = 16;
    const h1 = helixAt(H, [-6.6, -3, 0], [0, 0, 1]), h2 = helixAt(H, [6.6, -3, 0], [0, 0, -1], 0.8), h3 = helixAt(H, [0, 8, 0], [0, 0, 1], 1.9);
    const [sa, sb] = antiparallelPair(7), Rs = rotAxis([1, 0, 0], Math.PI / 2);
    const s1 = moveBB(sa, xf(Rs, [0, -15.5, -2])), s2 = moveBB(sb, xf(Rs, [0, -15.5, -2]));
    const avoid = [].concat(h1, h2, h3, s1, s2).map(r => r.CA);
    const L1 = loop(h1[H - 1], h2[0], avoid), L2 = loop(h2[H - 1], h3[0], avoid), L3 = loop(h3[H - 1], s1[0], avoid), L4 = loop(s1[s1.length - 1], s2[0], avoid);
    const elems = [h1, L1, h2, L2, h3, L3, s1, L4, s2];
    addChain(sys, [].concat(...elems), ['LEU', 'ALA', 'GLU', 'LYS', 'LEU', 'ILE', 'SER', 'ALA', 'VAL', 'GLN', 'LEU', 'ASP', 'ARG', 'THR', 'MET', 'ASN'], 'A', 1);
    const ss = [];
    let k0 = 1;
    for (const e of elems) { const t = e === h1 || e === h2 || e === h3 ? 'helix' : e === s1 || e === s2 ? 'sheet' : null; if (t) ss.push({ type: t, chain: 'A', start: { resId: k0 }, end: { resId: k0 + e.length - 1 } }); k0 += e.length; }
    sys.ss = ss; sys.loops = { L1: [h1.length + 1, h1.length + L1.length], L2: [h1.length + L1.length + h2.length + 1, h1.length + L1.length + h2.length + L2.length] };
    if (BM && BM.chem && opts.ligand !== false) {
      const lig = BM.chem.smilesTo3D(opts.ligandSmiles || 'Cc1ccc(cc1)C(=O)Nc1ccc(O)cc1', { seed: 3 }).mol, nl = lig.Z.length;
      const lp = i => [lig.xyz[3 * i], lig.xyz[3 * i + 1], lig.xyz[3 * i + 2]], c = centroid(Array.from({ length: nl }, (_, i) => lp(i)));
      const ax = axisOf(Array.from({ length: nl }, (_, i) => lp(i))), T = xf(rotTo(ax.d, [0.2, 0.25, 1]), [0, 0.6, 0.5]);
      const s0 = sys.Z.length;
      sys.fixed = new Set();
      for (let i = 0; i < nl; i++) sys.fixed.add(push(sys, (BM.el.symbol(lig.Z[i]) + (i + 1)).slice(0, 4), lig.Z[i], T(sub(lp(i), c)), 'LIG', 900, 'L'));
      sys.bonded = new Set(lig.bonds.map(b => (s0 + b[0]) + ':' + (s0 + b[1])));
      for (const b of lig.bonds) sys.bonds.push([s0 + b[0], s0 + b[1], b[2]]);
    }
    if (opts.water !== false) {
      let w = 0;
      const ions = [[16.6, 15.4, 14.2], [-16.8, -20.6, 13.4]];
      for (let x = -18; x <= 18; x += 6.0) for (let y = -24; y <= 18; y += 6.0) for (let z = -16; z <= 16; z += 6.0) {
        const p = [x + 0.7 * Math.sin(w), y + 0.6 * Math.cos(1.3 * w), z];
        w++;
        if (ions.some(q => len(sub(q, p)) < 4)) continue;
        push(sys, 'OW', 8, p, 'SOL', 2000 + w, 'W'); push(sys, 'HW1', 1, add(p, [0.757, 0.586, 0]), 'SOL', 2000 + w, 'W'); push(sys, 'HW2', 1, add(p, [-0.757, 0.586, 0]), 'SOL', 2000 + w, 'W');
      }
      push(sys, 'NA', 11, ions[0], 'NA', 3001, 'I'); push(sys, 'CL', 17, ions[1], 'CL', 3002, 'I');
    }
    declash(sys);
    return sys;
  }

  /** POPC-like lipid (CHARMM-style names) with hydrogens; head at the origin, tails along -z. [name, Z, xyz][] + heavy bonds. */
  function lipidTemplate() {
    const at = [], bonds = [], idx = {};
    const A = (n, z, p) => { idx[n] = at.length; at.push([n, z, p]); };
    const B = (a, b) => bonds.push([idx[a], idx[b]]);
    A('N', 7, [0, 0, 0]); A('C13', 6, [1.0, 0.6, 0.9]); A('C14', 6, [-1.0, 0.6, 0.9]); A('C15', 6, [0, -1.2, 0.9]); A('C12', 6, [0, 0.3, -1.47]);
    A('C11', 6, [0.3, -0.5, -2.68]); A('O12', 8, [0.1, 0.2, -3.88]); A('P', 15, [0.4, -0.4, -5.3]); A('O13', 8, [1.85, -0.65, -5.4]); A('O14', 8, [-0.4, -1.6, -5.5]);
    A('O11', 8, [0.0, 0.7, -6.4]); A('C1', 6, [0.3, 0.6, -7.8]); A('C2', 6, [-0.2, 1.7, -8.65]); A('C3', 6, [0.5, 1.6, -10.0]);
    A('O21', 8, [-1.6, 1.8, -8.85]); A('O31', 8, [1.9, 1.5, -10.2]);
    for (const [a, b] of [['N', 'C13'], ['N', 'C14'], ['N', 'C15'], ['N', 'C12'], ['C12', 'C11'], ['C11', 'O12'], ['O12', 'P'], ['P', 'O13'], ['P', 'O14'], ['P', 'O11'], ['O11', 'C1'],
      ['C1', 'C2'], ['C2', 'C3'], ['C2', 'O21'], ['C3', 'O31']]) B(a, b);
    // all-trans tails (C-C 1.53 Å, C-C-C 112°) tilted by 25° towards +x; carbonyl O along the trigonal direction
    const tilt = 25 * D, ax = [Math.sin(tilt), 0, -Math.cos(tilt)], zz = [Math.cos(tilt), 0, Math.sin(tilt)];
    const tail = (pref, start, nC, carbonylO, esterO) => {
      for (let k = 1; k <= nC; k++) { A(pref + k, 6, add(start, add(mul(ax, 1.27 * (k - 1)), mul(zz, k % 2 ? 0.43 : -0.43)))); if (k > 1) B(pref + (k - 1), pref + k); }
      const c1 = at[idx[pref + 1]][2], u1 = norm(sub(at[idx[esterO]][2], c1)), u2 = norm(sub(at[idx[pref + 2]][2], c1));
      A(carbonylO, 8, add(c1, mul(norm(mul(add(u1, u2), -1)), 1.23))); B(pref + 1, carbonylO); B(esterO, pref + 1);
    };
    tail('C2', [-2.8, 1.8, -10.19], 18, 'O22', 'O21');
    tail('C3', [2.6, 1.5, -11.3], 16, 'O32', 'O31');
    // hydrogens from the heavy-atom bonds (tetrahedral directions); carbonyl carbons and N get none
    const nb = at.map(() => []);
    for (const [a, b] of bonds) { nb[a].push(b); nb[b].push(a); }
    const heavy = at.length;
    for (let i = 0; i < heavy; i++) {
      const [name, z, p] = at[i];
      if (z !== 6 || name === 'C21' || name === 'C31') continue;
      const nH = 4 - nb[i].length;
      if (nH <= 0) continue;
      const u = nb[i].map(j => norm(sub(at[j][2], p))), base = norm(mul(u.reduce((s, v) => add(s, v), [0, 0, 0]), -1));
      let dirs;
      if (nH === 1) dirs = [base];
      else if (nH === 2) { const w = norm(cross(u[0], u[1])); dirs = [norm(add(mul(base, 0.58), mul(w, 0.82))), norm(add(mul(base, 0.58), mul(w, -0.82)))]; }
      else { const a = mul(u[0], -1), q = perp(a), r = cross(a, q); dirs = [0, 1, 2].map(k => norm(add(mul(a, 0.33), mul(add(mul(q, Math.cos(k * 2.094)), mul(r, Math.sin(k * 2.094))), 0.94)))); }
      dirs.forEach((d, k) => A('H' + name.slice(1) + 'RST'[k], 1, add(p, mul(d, 1.0))));
    }
    return { atoms: at, bonds };
  }
  /** Bilayer: n x n lipids per leaflet (choline N at z = ±32 Å), TIP3 water slabs, SOD/CLA; optional 4-helix TM bundle. */
  function membrane(opts) {
    opts = opts || {};
    const n = opts.n || 8, sp = opts.spacing || 9.0, tpl = lipidTemplate().atoms, sys = newSys();
    const half = (n - 1) * sp / 2, protR = opts.protein ? 13 : 0;
    let rid = 1;
    if (opts.protein) {
      for (let q = 0; q < 4; q++) {
        const ang = q * Math.PI / 2 + 0.4;
        addChain(sys, helixAt(36, [7 * Math.cos(ang), 7 * Math.sin(ang), 0], [0.12 * Math.sin(ang), -0.12 * Math.cos(ang), q % 2 ? -1 : 1], q), ['LEU', 'ALA', 'VAL', 'ILE', 'LEU', 'GLY', 'SER', 'LEU', 'ALA', 'LYS'], 'P' + 'ABCD'[q], rid);
        rid += 36;
      }
    }
    for (const side of [1, -1]) {
      for (let a = 0; a < n; a++) for (let b = 0; b < n; b++) {
        const x = -half + a * sp + (b % 2 ? sp / 4 : -sp / 4), y = -half + b * sp;
        if (protR && Math.hypot(x, y) < protR) continue;
        const ang = Math.sin(a * 7.1 + b * 13.7 + (side > 0 ? 0 : 5)) * 0.12 + (side > 0 ? 0 : 2.2), c = Math.cos(ang), s = Math.sin(ang);
        for (const [nm, z, p] of tpl) push(sys, nm, z, [x + c * p[0] - s * p[1], y + s * p[0] + c * p[1], side * (32 + p[2])], 'POPC', rid, 'MEMB');
        rid++;
      }
    }
    if (opts.water !== false) {
      const L = half + sp / 2;
      let w = 0;
      for (const zc of [40, -40]) for (let x = -L; x <= L; x += 3.4) for (let y = -L; y <= L; y += 3.4) for (const dz of [-3.3, 0, 3.3]) {
        const p = [x, y, zc + dz];
        w++;
        push(sys, 'OH2', 8, p, 'TIP3', rid, 'TIP3'); push(sys, 'H1', 1, add(p, [0.757, 0.586, 0]), 'TIP3', rid, 'TIP3'); push(sys, 'H2', 1, add(p, [-0.757, 0.586, 0]), 'TIP3', rid, 'TIP3');
        rid++;
      }
      const ni = Math.max(4, Math.round(w / 55));
      for (let q = 0; q < ni; q++) { const t = q * 2.399, r = L * 0.85 * Math.sqrt((q + 0.5) / ni); push(sys, q % 2 ? 'CLA' : 'SOD', q % 2 ? 17 : 11, [Math.cos(t) * r, Math.sin(t) * r, (q % 4 < 2 ? 1 : -1) * 46.5], q % 2 ? 'CLA' : 'SOD', rid++, q % 2 ? 'CLA' : 'SOD'); }
    }
    if (opts.protein) declash(sys);
    return sys;
  }

  /** Crude B-DNA-like single strand (names DA/DT/DG/DC, chain N): backbone atoms on a 36°/3.38 Å helix and planar
   *  base rings stacked near the axis (geometry for cartoon tests, not chemistry). */
  function nucleic(n, opts) {
    opts = opts || {};
    const sys = newSys(), seq = ['DA', 'DT', 'DG', 'DC'], chain = opts.chain || 'N';
    const cyl = (r, deg, z) => [r * Math.cos(deg * D), r * Math.sin(deg * D), z];
    for (let k = 0; k < n; k++) {
      const th = k * 36 + (opts.phase || 0), z0 = k * 3.38, rn = seq[k % 4], pur = rn === 'DA' || rn === 'DG';
      const at = [['P', 15, cyl(8.9, th, z0)], ['OP1', 8, cyl(10.2, th - 4, z0 + 0.8)], ['OP2', 8, cyl(10.0, th + 5, z0 - 0.9)], ["O5'", 8, cyl(8.3, th + 8, z0 + 0.6)],
        ["C5'", 6, cyl(8.6, th + 16, z0 + 1.3)], ["C4'", 6, cyl(7.8, th + 22, z0 + 1.8)], ["O4'", 8, cyl(6.7, th + 25, z0 + 1.2)], ["C3'", 6, cyl(8.4, th + 30, z0 + 2.6)],
        ["O3'", 8, cyl(8.8, th + 31, z0 + 3.0)], ["C2'", 6, cyl(7.4, th + 36, z0 + 2.4)], ["C1'", 6, cyl(5.9, th + 30, z0 + 1.5)]];
      // base: hexagon centred 2 Å from the axis in the plane z = z0 + 1.5, five-ring / N1 towards the backbone
      const ang = (th + 30) * D, u = [Math.cos(ang), Math.sin(ang), 0], v = [-Math.sin(ang), Math.cos(ang), 0], zc = z0 + 1.5;
      const P2 = (a, b) => [u[0] * a + v[0] * b, u[1] * a + v[1] * b, zc];
      const hc = 2.0, hex = (deg) => P2(hc + 1.4 * Math.cos(deg * D), 1.4 * Math.sin(deg * D));
      if (pur) {
        [['C5', 6, 30], ['C4', 6, -30], ['N3', 7, -90], ['C2', 6, -150], ['N1', 7, 150], ['C6', 6, 90]].forEach(([nm, z, d]) => at.push([nm, z, hex(d)]));
        const pc = hc + 1.21 + 0.963, R5 = 1.191;
        [['N9', 7, -72], ['C8', 6, 0], ['N7', 7, 72]].forEach(([nm, z, d]) => at.push([nm, z, P2(pc + R5 * Math.cos(d * D), R5 * Math.sin(d * D))]));
      } else [['N1', 7, 0], ['C2', 6, -60], ['N3', 7, -120], ['C4', 6, 180], ['C5', 6, 120], ['C6', 6, 60]].forEach(([nm, z, d]) => at.push([nm, z, P2(hc + 0.5 + 1.4 * Math.cos(d * D), 1.4 * Math.sin(d * D))]));
      for (const [nm, z, p] of at) push(sys, nm, z, p, rn, k + 1, chain);
    }
    return sys;
  }

  // ---------------------------------------------------------------- glycans
  // monosaccharides as SMILES with heavy atoms in a known order (names), C1 and O1 (removed when linking)
  const MONO = {
    Xyl: { smi: 'C1[C@H]([C@@H]([C@H]([C@@H](O1)O)O)O)O', names: ['C5', 'C4', 'C3', 'C2', 'C1', 'O5', 'O1', 'O2', 'O3', 'O4'] },
    Gal: { smi: 'C([C@@H]1[C@@H]([C@@H]([C@H]([C@@H](O1)O)O)O)O)O', names: ['C6', 'C5', 'C4', 'C3', 'C2', 'C1', 'O5', 'O1', 'O2', 'O3', 'O4', 'O6'] },
    Glc: { smi: 'C([C@@H]1[C@H]([C@@H]([C@H]([C@@H](O1)O)O)O)O)O', names: ['C6', 'C5', 'C4', 'C3', 'C2', 'C1', 'O5', 'O1', 'O2', 'O3', 'O4', 'O6'] },
    Man: { smi: 'C([C@@H]1[C@H]([C@@H]([C@@H]([C@H](O1)O)O)O)O)O', names: ['C6', 'C5', 'C4', 'C3', 'C2', 'C1', 'O5', 'O1', 'O2', 'O3', 'O4', 'O6'] },
    Fuc: { smi: 'C[C@H]1[C@H]([C@H]([C@@H]([C@@H](O1)O)O)O)O', names: ['C6', 'C5', 'C4', 'C3', 'C2', 'C1', 'O5', 'O1', 'O2', 'O3', 'O4'] },
    GlcA: { smi: 'C(=O)([C@@H]1[C@H]([C@@H]([C@H]([C@@H](O1)O)O)O)O)O', names: ['C6', 'O6A', 'C5', 'C4', 'C3', 'C2', 'C1', 'O5', 'O1', 'O2', 'O3', 'O4', 'O6B'] },
    IdoA: { smi: 'C(=O)([C@H]1[C@H]([C@@H]([C@H]([C@@H](O1)O)O)O)O)O', names: ['C6', 'O6A', 'C5', 'C4', 'C3', 'C2', 'C1', 'O5', 'O1', 'O2', 'O3', 'O4', 'O6B'] },
    GlcNAc: { smi: 'C([C@@H]1[C@H]([C@@H]([C@H]([C@H](O1)O)NC(C)=O)O)O)O', names: ['C6', 'C5', 'C4', 'C3', 'C2', 'C1', 'O5', 'O1', 'N2', 'C7', 'C8', 'O7', 'O3', 'O4', 'O6'] },
    GlcN: { smi: 'C([C@@H]1[C@H]([C@@H]([C@H]([C@H](O1)O)N)O)O)O', names: ['C6', 'C5', 'C4', 'C3', 'C2', 'C1', 'O5', 'O1', 'N2', 'O3', 'O4', 'O6'] },
    GlcNS6S: { smi: 'C([C@@H]1[C@H]([C@@H]([C@H]([C@H](O1)O)NS(=O)(=O)O)O)O)OS(=O)(=O)O', names: ['C6', 'C5', 'C4', 'C3', 'C2', 'C1', 'O5', 'O1', 'N2', 'S2', 'O2S1', 'O2S2', 'O2S3', 'O3', 'O4', 'O6', 'S6', 'O6S1', 'O6S2', 'O6S3'] },
    IdoA2S: { smi: 'C(=O)([C@H]1[C@H]([C@@H]([C@H]([C@@H](O1)O)OS(=O)(=O)O)O)O)O', names: ['C6', 'O6A', 'C5', 'C4', 'C3', 'C2', 'C1', 'O5', 'O1', 'O2', 'S2', 'O2S1', 'O2S2', 'O2S3', 'O3', 'O4', 'O6B'] },
    Neu5Ac: { smi: 'CC(=O)N[C@@H]1[C@@H](O)C[C@@](O)(O[C@H]1[C@H](O)[C@H](O)CO)C(O)=O', names: ['C11', 'C10', 'O10', 'N5', 'C5', 'C4', 'O4', 'C3', 'C2', 'O2', 'O6', 'C6', 'C7', 'O7', 'C8', 'O8', 'C9', 'O9', 'C1', 'O1B', 'O1A'] }
  };
  /** Heparan-sulfate proteoglycan fragment: the bundle protein (no ligand) with Xyl-Gal-Gal-GlcA linker on a loop serine and
   *  GlcNAc/GlcNS6S-uronate repeats. Residue names mix PDB ids, CHARMM-GUI truncations (BGLC, AGLC, AIDO) and unknown names
   *  (HS1, HS2) that only the ring topology identifies. sys.expected = [{resn, rid, label}]. */
  function glycan(BM, opts) {
    opts = opts || {};
    const sys = protein(BM, { ligand: false, water: false });
    const sid = sys.loops.L2[0] + 1;
    const at = nm => sys.names.findIndex((x, i) => x === nm && sys.ids[i] === sid && sys.chains[i] === 'A');
    for (let i = 0; i < sys.Z.length; i++) if (sys.ids[i] === sid && sys.chains[i] === 'A') sys.res[i] = 'SER';
    const ca = at('CA'), CB = place(pos(sys, at('C')), pos(sys, at('N')), pos(sys, ca), 1.53, 110.5, -122.6);
    const OG = place(pos(sys, at('N')), pos(sys, ca), CB, 1.43, 111, 180);
    let accC = push(sys, 'CB', 6, CB, 'SER', sid, 'A'), accO = push(sys, 'OG', 8, OG, 'SER', sid, 'A');
    const seq = opts.seq || [['Xyl', 'XYS'], ['Gal', 'GAL', 'O4'], ['Gal', 'GAL', 'O4'], ['GlcA', 'BGLC', 'O3'], ['GlcNAc', 'NAG', 'O4'], ['GlcA', 'BDP', 'O4'],
      ['GlcNS6S', 'AGLC', 'O4'], ['IdoA2S', 'AIDO', 'O4'], ['GlcNS6S', 'HS1', 'O4'], ['IdoA2S', 'HS2', 'O4'], ['GlcNS6S', 'SGN', 'O4'], ['IdoA2S', 'IDS', 'O4']];
    let rid = 501, prev = null;
    const expected = [];
    for (const [kind, resn, accName] of seq) {
      if (prev) { accO = prev[accName]; accC = prev['C' + accName.slice(1)]; }
      const m = MONO[kind], mol = BM.chem.smilesTo3D(m.smi, { seed: 7 }).mol, nh = m.names.length;
      const P = i => [mol.xyz[3 * i], mol.xyz[3 * i + 1], mol.xyz[3 * i + 2]];
      const iC1 = m.names.indexOf('C1'), iO1 = m.names.indexOf('O1');
      const drop = new Set([iO1]);
      for (const b of mol.bonds) { if (b[0] === iO1 && mol.Z[b[1]] === 1) drop.add(b[1]); if (b[1] === iO1 && mol.Z[b[0]] === 1) drop.add(b[0]); }
      // the acceptor hydroxyl loses its H (atoms of the previous residue keep their indices: H atoms come last)
      for (let i = sys.Z.length - 1; i >= 0; i--) if (sys.Z[i] === 1 && len(sub(pos(sys, i), pos(sys, accO))) < 1.15) sys.Z[i] = -1;
      const remap = compact(sys);
      accO = remap.get(accO); accC = remap.get(accC);
      if (prev) for (const k of Object.keys(prev)) prev[k] = remap.get(prev[k]);
      const target = pos(sys, accO), back = norm(sub(pos(sys, accC), target)), u = perp(back), v = cross(back, u);
      let best = null;
      for (let q = 0; q < 18; q++) for (let w = 0; w < 18; w++) {
        const ph = q * Math.PI / 9, dir = norm(add(mul(back, Math.cos(117 * D)), mul(add(mul(u, Math.cos(ph)), mul(v, Math.sin(ph))), Math.sin(117 * D))));
        const R = matMul(rotAxis(dir, w * Math.PI / 9), rotTo(sub(P(iC1), P(iO1)), dir)), Rx = xf(R, [0, 0, 0]);
        const T = p => add(target, Rx(sub(p, P(iO1))));
        let clash = 0;
        for (let i = 0; i < mol.Z.length; i++) if (!drop.has(i) && i !== iC1) clash += clashScore(sys, T(P(i)), accO, mol.Z[i] === 1 ? 2.2 : 3.2);
        if (!best || clash < best.clash) best = { clash, T };
      }
      if (root.SYNTH_DEBUG) console.log('glycan', resn, 'clash', best.clash.toFixed(2));
      const atoms = {}, map = new Map();
      for (let i = 0; i < mol.Z.length; i++) {
        if (drop.has(i)) continue;
        const nm = i < nh ? m.names[i] : 'H' + (i - nh + 1);
        map.set(i, push(sys, nm, mol.Z[i], best.T(P(i)), resn, rid, 'G'));
        atoms[nm] = map.get(i);
      }
      for (const b of mol.bonds) if (map.has(b[0]) && map.has(b[1])) sys.bonds.push([map.get(b[0]), map.get(b[1]), b[2]]);
      sys.bonds.push([accO, atoms.C1, 1]);
      expected.push({ resn, rid, label: kind });
      prev = atoms; rid++;
    }
    sys.expected = expected;
    if (opts.ions !== false) { push(sys, 'SOD', 11, [20, 20, 20], 'SOD', 3001, 'I'); push(sys, 'SOD', 11, [-22, 18, -20], 'SOD', 3002, 'I'); push(sys, 'CLA', 17, [22, -24, 18], 'CLA', 3003, 'I'); }
    return sys;
  }
  function clashScore(sys, p, skip, lim) {
    let s = 0;
    for (let i = 0; i < sys.Z.length; i++) {
      if (i === skip || (sys.Z[i] <= 1 && lim > 2.5)) continue;
      const d = Math.hypot(sys.xyz[3 * i] - p[0], sys.xyz[3 * i + 1] - p[1], sys.xyz[3 * i + 2] - p[2]);
      if (d < lim) s += (lim - d) * (lim - d);
    }
    return s;
  }

  // ---------------------------------------------------------------- writers
  const SYM = ['X', 'H', 'He', 'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne', 'Na', 'Mg', 'Al', 'Si', 'P', 'S', 'Cl'];
  const sym = z => (root.BM && root.BM.el ? root.BM.el.symbol(z) : SYM[z] || 'X');
  const AA = /^(ALA|ARG|ASN|ASP|CYS|GLN|GLU|GLY|HIS|ILE|LEU|LYS|MET|PHE|PRO|SER|THR|TRP|TYR|VAL)$/;
  function toPDB(sys, opts) {
    opts = opts || {};
    const L = ['HEADER    SYNTHETIC TEST STRUCTURE'];
    if (opts.ss && sys.ss) {
      let h = 0, s = 0;
      for (const r of sys.ss) {
        if (r.type === 'helix') { h++; L.push(('HELIX  ' + String(h).padStart(3) + ' ' + String(h).padStart(3) + ' GLY ' + r.chain + ' ' + String(r.start.resId).padStart(4) + '  GLY ' + r.chain + ' ' + String(r.end.resId).padStart(4) + '  1').padEnd(80)); }
        else { s++; L.push(('SHEET  ' + String(s).padStart(3) + '  S1 2 GLY ' + r.chain + String(r.start.resId).padStart(4) + '  GLY ' + r.chain + String(r.end.resId).padStart(4) + '  0').padEnd(80)); }
      }
    }
    const f = v => v.toFixed(3).padStart(8);
    for (let i = 0; i < sys.Z.length; i++) {
      const nm = sys.names[i], name4 = nm.length >= 4 ? nm.slice(0, 4) : ' ' + nm.padEnd(3);
      L.push((AA.test(sys.res[i]) ? 'ATOM  ' : 'HETATM') + String((i + 1) % 100000).padStart(5) + ' ' + name4 + ' ' + sys.res[i].slice(0, 4).padEnd(4) + (sys.chains[i] || ' ').slice(0, 1) +
        String(sys.ids[i] % 10000).padStart(4) + '    ' + f(sys.xyz[3 * i]) + f(sys.xyz[3 * i + 1]) + f(sys.xyz[3 * i + 2]) + '  1.00' + (opts.bfac ? opts.bfac(i).toFixed(2).padStart(6) : '  0.00') + '          ' + sym(sys.Z[i]).padStart(2));
    }
    L.push('END');
    return L.join('\n') + '\n';
  }
  function toGRO(sys, title) {
    const L = [title || 'synthetic system', String(sys.Z.length)], mx = [-Infinity, -Infinity, -Infinity], mn = [Infinity, Infinity, Infinity];
    for (let i = 0; i < sys.Z.length; i++) {
      const p = [sys.xyz[3 * i] / 10, sys.xyz[3 * i + 1] / 10, sys.xyz[3 * i + 2] / 10];
      for (let k = 0; k < 3; k++) { if (p[k] > mx[k]) mx[k] = p[k]; if (p[k] < mn[k]) mn[k] = p[k]; }
      L.push(String(sys.ids[i] % 100000).padStart(5) + sys.res[i].slice(0, 5).padEnd(5) + sys.names[i].slice(0, 5).padStart(5) + String((i + 1) % 100000).padStart(5) + p.map(v => v.toFixed(3).padStart(8)).join(''));
    }
    L.push(mx.map((v, k) => (v - mn[k] + 0.5).toFixed(5).padStart(10)).join(''));
    return L.join('\n') + '\n';
  }
  /** Dataset for BM.app.openData / BM.reps.context (no explicit bonds: the app bonds by distance). */
  function toDataset(sys, name) {
    return { format: 'pdb', program: null, title: name || 'synthetic', meta: { route: null, keywords: null, method: null, basis: null, charge: null, multiplicity: null, jobTypes: [], normalTermination: null, solvent: null },
      frames: [{ Z: sys.Z.slice(), xyz: Float64Array.from(sys.xyz), energy: null, energies: null, label: null, converged: false, opt: null, atomNames: sys.names.slice(), resNames: sys.res.slice(), resIds: sys.ids.slice(), chains: sys.chains.slice() }],
      vibrations: null, excitedStates: null, charges: null, spins: null, dipole: null, orbitals: null, thermo: null, volumes: null, bonds: null, cell: null, warnings: [], secondaryStructure: sys.ss || null };
  }

  root.BMSynth = { place, dihedral, backbone, idealHelix, idealStrand, helixAt, antiparallelPair, residueAtoms, protein, membrane, lipidTemplate, glycan, nucleic, MONO, toPDB, toGRO, toDataset, hbE, declash };
  if (typeof module !== 'undefined' && module.exports) module.exports = root.BMSynth;
})(typeof window !== 'undefined' ? window : globalThis);
