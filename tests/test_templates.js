#!/usr/bin/env node
/* BlueMolecule building-library tests.  Run:  node tests/test_templates.js */
'use strict';
const path = require('path');
const ROOT = path.join(__dirname, '..');
for (const f of ['elements', 'geom', 'ff', 'templates']) require(path.join(ROOT, 'src', f + '.js'));
const BM = globalThis.BM, T = BM.templates, FF = BM.ff, G = BM.geom;

let pass = 0, fail = 0, current = '';
const failures = [];
function ok(cond, msg) { if (cond) pass++; else { fail++; failures.push('[' + current + '] ' + msg); } }
function near(a, b, tol, msg) { ok(typeof a === 'number' && Math.abs(a - b) <= tol, msg + ': got ' + (typeof a === 'number' ? a.toFixed(4) : a) + ', expected ' + b + ' ± ' + tol); }
function test(name, fn) {
  current = name;
  const t0 = Date.now();
  try { fn(); } catch (e) { fail++; failures.push('[' + name + '] threw: ' + (e && e.stack || e)); }
  console.log('  ' + name + ' (' + (Date.now() - t0) + ' ms)');
}

// ---------------------------------------------------------------------------------------------- helpers
const EMPTY = { Z: [], xyz: [], bonds: [] };
const ang = (x, i, j, k) => G.angle(x, i, j, k), dist = (x, i, j) => G.dist(x, i, j);
const sym = z => BM.el.symbol(z);
/** Sane bonds: every bond within tol of the shared ideal length. Returns the worst deviation. */
function worstBond(m) {
  let w = 0, what = '';
  for (const [i, j, o] of m.bonds) {
    if (!(m.Z[i] > 0) || !(m.Z[j] > 0)) continue;
    const d = Math.abs(dist(m.xyz, i, j) - FF.idealBondLength(m.Z[i], m.Z[j], o));
    if (d > w) { w = d; what = sym(m.Z[i]) + (i + 1) + '–' + sym(m.Z[j]) + (j + 1); }
  }
  return { w, what };
}
/** Clash criterion: atoms ≥ 3 bonds apart (or unbonded) closer than 0.7 × the sum of their vdW radii,
 * or any two atoms closer than 0.6 Å. */
function clashes(m) {
  const list = T.findClashes(m, { scale: 0.7 });
  for (let i = 0; i < m.Z.length; i++) for (let j = i + 1; j < m.Z.length; j++) if (dist(m.xyz, i, j) < 0.6) list.push([i, j, dist(m.xyz, i, j)]);
  return list;
}
function sane(m, label, tol) {
  const wb = worstBond(m), cl = clashes(m);
  ok(wb.w <= (tol || 0.08), label + ': bond lengths sane (worst ' + wb.what + ' off by ' + wb.w.toFixed(3) + ' Å)');
  ok(cl.length === 0, label + ': no clashes' + (cl.length ? ' (' + cl.slice(0, 3).map(c => sym(m.Z[c[0]]) + (c[0] + 1) + '…' + sym(m.Z[c[1]]) + (c[1] + 1) + ' ' + c[2].toFixed(2)).join(', ') + ')' : ''));
  ok(Array.from(m.xyz).every(Number.isFinite), label + ': finite coordinates');
}
const donorsOf = (m, M) => m.bonds.filter(b => b[0] === M || b[1] === M).map(b => (b[0] === M ? b[1] : b[0]));
function ligandAngles(m, M) {
  const d = donorsOf(m, M), a = [];
  for (let p = 0; p < d.length; p++) for (let q = p + 1; q < d.length; q++) a.push(ang(m.xyz, d[p], M, d[q]));
  return a;
}
const neighbours = (m, i) => m.bonds.filter(b => b[0] === i || b[1] === i).map(b => (b[0] === i ? b[1] : b[0]));
function planeDev(X, idx) { // via the smallest-eigenvalue normal (3x3 Jacobi)
  const n = idx.length, c = [0, 0, 0];
  for (const i of idx) for (let k = 0; k < 3; k++) c[k] += X[3 * i + k] / n;
  const a = [[0, 0, 0], [0, 0, 0], [0, 0, 0]], v = [[1, 0, 0], [0, 1, 0], [0, 0, 1]];
  for (const i of idx) { const d = [0, 1, 2].map(k => X[3 * i + k] - c[k]); for (let p = 0; p < 3; p++) for (let q = 0; q < 3; q++) a[p][q] += d[p] * d[q]; }
  for (let s = 0; s < 60; s++) for (let p = 0; p < 3; p++) for (let q = p + 1; q < 3; q++) {
    if (Math.abs(a[p][q]) < 1e-300) continue;
    const th = (a[q][q] - a[p][p]) / (2 * a[p][q]), t = (th >= 0 ? 1 : -1) / (Math.abs(th) + Math.sqrt(th * th + 1)), cs = 1 / Math.sqrt(t * t + 1), sn = t * cs;
    for (let k = 0; k < 3; k++) { const x1 = a[k][p], x2 = a[k][q]; a[k][p] = cs * x1 - sn * x2; a[k][q] = sn * x1 + cs * x2; }
    for (let k = 0; k < 3; k++) { const x1 = a[p][k], x2 = a[q][k]; a[p][k] = cs * x1 - sn * x2; a[q][k] = sn * x1 + cs * x2; }
    for (let k = 0; k < 3; k++) { const x1 = v[k][p], x2 = v[k][q]; v[k][p] = cs * x1 - sn * x2; v[k][q] = sn * x1 + cs * x2; }
  }
  let m = 0;
  for (let k = 1; k < 3; k++) if (a[k][k] < a[m][m]) m = k;
  let mx = 0;
  for (const i of idx) mx = Math.max(mx, Math.abs((X[3 * i] - c[0]) * v[0][m] + (X[3 * i + 1] - c[1]) * v[1][m] + (X[3 * i + 2] - c[2]) * v[2][m]));
  return mx;
}
const snapshot = m => JSON.stringify({ Z: Array.from(m.Z), xyz: Array.from(m.xyz), bonds: m.bonds, fc: m.formalCharges || null });

console.log('BlueMolecule template tests');

// ---------------------------------------------------------------------------------------------- tests
test('library contents and metadata', () => {
  const ids = list => list.map(e => e.id);
  const want = {
    geometries: ['linear', 'trigonal-planar', 'tetrahedral', 'square-planar', 'seesaw', 't-shaped', 'trigonal-bipyramidal', 'square-pyramidal', 'octahedral', 'pentagonal-bipyramidal'],
    ligands: ['H', 'F', 'Cl', 'Br', 'I', 'OH', 'H2O', 'NH3', 'CO', 'CN', 'NCMe', 'PH3', 'PMe3', 'PPh3', 'py', 'DMSO', 'en', 'bipy', 'ox', 'acac', 'phen'],
    groups: ['methyl', 'ethyl', 'isopropyl', 'tert-butyl', 'phenyl', 'benzyl', 'vinyl', 'ethynyl', 'OH', 'OMe', 'NH2', 'NMe2', 'NO2', 'CN', 'CF3', 'CHO', 'COOH', 'COOMe', 'CONH2', 'acetyl', 'SH', 'SO3H', 'F', 'Cl', 'Br', 'I'],
    molecules: ['water', 'ammonia', 'methane', 'ethane', 'ethylene', 'acetylene', 'benzene', 'toluene', 'cyclohexane', 'cyclopentane', 'pyridine', 'pyrrole', 'furan', 'thiophene', 'imidazole', 'naphthalene', 'methanol', 'ethanol', 'acetone', 'acetonitrile', 'dmso', 'dichloromethane', 'chloroform', 'thf', 'co2', 'formaldehyde', 'acetic-acid', 'cisplatin', 'transplatin']
  };
  for (const k of Object.keys(want)) for (const id of want[k]) ok(ids(T[k]).includes(id), k + ' contains ' + id);
  for (const k of Object.keys(want)) for (const e of T[k]) ok(typeof e.id === 'string' && typeof e.name === 'string' && e.name && typeof e.category === 'string' && e.category, k + ' entry ' + e.id + ' has id/name/category');
  for (const g of T.geometries) {
    ok(g.cn === g.vectors.length, g.id + ': cn matches vectors');
    for (const v of g.vectors) near(Math.hypot(...v), 1, 1e-12, g.id + ' unit vectors');
  }
  const charges = { H: -1, F: -1, Cl: -1, Br: -1, I: -1, OH: -1, H2O: 0, NH3: 0, CO: 0, CN: -1, NCMe: 0, PH3: 0, PMe3: 0, PPh3: 0, py: 0, DMSO: 0, en: 0, bipy: 0, ox: -2, acac: -1, phen: 0 };
  const donorEl = { CO: 'C', CN: 'C', NCMe: 'N', py: 'N', DMSO: 'S', PPh3: 'P', H2O: 'O', en: 'N', bipy: 'N', ox: 'O', acac: 'O', phen: 'N' };
  for (const [id, q] of Object.entries(charges)) {
    const L = T.get('ligand', id);
    ok(L.charge === q, id + ' charge ' + q);
    ok(Array.isArray(L.donors) && L.donors.length === L.denticity, id + ' donor indices');
    if (donorEl[id]) for (const d of L.donors) ok(sym(L.mol.Z[d]) === donorEl[id], id + ' binds through ' + donorEl[id]);
  }
  for (const g of T.groups) ok(g.attach === 0 && g.direction.length === 3, g.id + ' attachment atom and direction');
  ok(T.get('group', 'nope') === null && T.get('molecule', 'Benzene').id === 'benzene', 'get() lookup');
  const a = T.get('molecule', 'benzene').mol; a.xyz[0] = 99;
  ok(T.get('molecule', 'benzene').mol.xyz[0] !== 99, 'entries hand out copies');
});

test('every molecule template inserts cleanly', () => {
  let host = T.insertMolecule(EMPTY, 'water', { at: [-30, 0, 0] }).mol;
  for (const e of T.molecules) {
    const r = T.insertMolecule(EMPTY, e.id, { at: [1, -2, 3] });
    sane(r.mol, e.id);
    ok(r.added.length === r.mol.Z.length && r.added[0] === 0, e.id + ': added indices');
    const c = [0, 1, 2].map(k => r.added.reduce((s, i) => s + r.mol.xyz[3 * i + k], 0) / r.added.length);
    near(Math.hypot(c[0] - 1, c[1] + 2, c[2] - 3), 0, 1e-6, e.id + ': centred at `at`');
    const before = snapshot(host);
    const r2 = T.insertMolecule(host, e.id, { at: [0, 0, 0] });
    ok(snapshot(host) === before, e.id + ': input not mutated');
    ok(r2.added[0] === host.Z.length && r2.mol.Z.length === host.Z.length + r2.added.length, e.id + ': appended at the end');
    for (let k = 0; k < host.xyz.length; k++) if (r2.mol.xyz[k] !== host.xyz[k]) { ok(false, e.id + ': existing atoms moved'); break; }
    host = r2.mol.Z.length < 400 ? r2.mol : host;
  }
  for (const e of [...T.groups, ...T.ligands]) sane(e.mol, e.kind + ' ' + e.id);
});

test('molecule template geometries', () => {
  const m = id => T.get('molecule', id).mol;
  let x = m('water'); near(ang(x.xyz, 1, 0, 2), 104.5, 0.5, 'water H–O–H');
  x = m('ammonia'); near(ang(x.xyz, 1, 0, 2), 107, 0.5, 'ammonia H–N–H');
  x = m('methane'); near(ang(x.xyz, 1, 0, 2), 109.47, 0.3, 'methane');
  x = m('co2'); near(ang(x.xyz, 1, 0, 2), 180, 0.5, 'CO2 linear');
  x = m('acetylene'); near(ang(x.xyz, 1, 0, 2), 180, 0.5, 'acetylene linear');
  x = m('ethane'); near(Math.abs(G.dihedral(x.xyz, 1, 0, 4, 5)) % 120, 60, 3, 'ethane staggered');
  x = m('ethylene'); ok(planeDev(x.xyz, [0, 1, 2, 3, 4, 5]) < 0.01, 'ethylene planar');
  x = m('benzene'); ok(planeDev(x.xyz, [...Array(12).keys()]) < 0.02, 'benzene planar');
  for (let k = 0; k < 6; k++) { const d = dist(x.xyz, k, (k + 1) % 6); ok(d > 1.38 && d < 1.42, 'benzene C–C ' + d.toFixed(3)); }
  x = m('naphthalene'); ok(planeDev(x.xyz, [...Array(18).keys()]) < 0.02, 'naphthalene planar'); ok(x.Z.filter(z => z === 1).length === 8, 'C10H8');
  x = m('cyclohexane');
  for (let k = 0; k < 6; k++) { const d = G.dihedral(x.xyz, k, (k + 1) % 6, (k + 2) % 6, (k + 3) % 6); ok(Math.abs(Math.abs(d) - 55) <= 8, 'cyclohexane chair dihedral ' + d.toFixed(1)); }
  for (const id of ['pyridine', 'pyrrole', 'furan', 'thiophene', 'imidazole']) { x = m(id); ok(planeDev(x.xyz, [...Array(x.Z.length).keys()]) < 0.02, id + ' planar'); }
  x = m('dmso'); { const S = 0, nb = neighbours(x, S); let sum = 0; for (let p = 0; p < 3; p++) for (let q = p + 1; q < 3; q++) sum += ang(x.xyz, nb[p], S, nb[q]); ok(sum < 340, 'DMSO sulfur pyramidal (angle sum ' + sum.toFixed(1) + ')'); }
  x = m('acetic-acid'); ok(x.Z.filter(z => z === 8).length === 2 && x.Z.length === 8, 'acetic acid C2H4O2');
  for (const [id, cisT] of [['cisplatin', true], ['transplatin', false]]) {
    x = m(id);
    const N = x.Z.map((z, i) => (z === 7 ? i : -1)).filter(i => i >= 0), Cl = x.Z.map((z, i) => (z === 17 ? i : -1)).filter(i => i >= 0);
    near(ang(x.xyz, N[0], 0, N[1]), cisT ? 90 : 180, 3, id + ' N–Pt–N'); near(ang(x.xyz, Cl[0], 0, Cl[1]), cisT ? 90 : 180, 3, id + ' Cl–Pt–Cl');
    near(dist(x.xyz, 0, N[0]), 2.07, 0.02, id + ' Pt–N'); near(dist(x.xyz, 0, Cl[0]), 2.38, 0.02, id + ' Pt–Cl');
  }
});

test('every group attaches to a benzene H and to a methane H', () => {
  for (const baseId of ['benzene', 'methane']) {
    const base = T.get('molecule', baseId).mol;
    const h = base.Z.indexOf(1), p = neighbours(base, h)[0];
    const dirH = [0, 1, 2].map(k => base.xyz[3 * h + k] - base.xyz[3 * p + k]);
    for (const g of T.groups) {
      const tag = g.id + ' on ' + baseId;
      const before = snapshot(base);
      const r = T.attachGroup(base, h, g.id), m = r.mol, gm = g.mol;
      ok(snapshot(base) === before, tag + ': input not mutated');
      ok(m.Z.length === base.Z.length - 1 + gm.Z.length, tag + ': atom count');
      ok(r.removed && r.removed.length === 1 && r.removed[0] === h, tag + ': H reported removed');
      ok(r.map && r.map[h] === -1 && r.map[p] === p - (p > h ? 1 : 0), tag + ': index map');
      ok(r.added.length === gm.Z.length && r.added[0] === base.Z.length - 1, tag + ': new atoms appended at the end');
      const A = r.added[0], pn = r.map[p];
      ok(neighbours(m, pn).includes(A), tag + ': bonded to the parent');
      near(dist(m.xyz, pn, A), FF.idealBondLength(m.Z[pn], m.Z[A], 1), 0.03, tag + ': parent–group bond');
      const dirA = [0, 1, 2].map(k => m.xyz[3 * A + k] - m.xyz[3 * pn + k]);
      const cos = dirA.reduce((s, v, k) => s + v * dirH[k], 0) / (Math.hypot(...dirA) * Math.hypot(...dirH));
      ok(Math.acos(Math.min(1, cos)) * 180 / Math.PI < 6, tag + ': along the old C–H direction');
      sane(m, tag);
      for (let i = 0; i < base.Z.length; i++) {
        if (i === h) continue;
        const j = r.map[i];
        if (m.xyz[3 * j] !== base.xyz[3 * i] || m.Z[j] !== base.Z[i]) { ok(false, tag + ': existing atoms untouched'); break; }
      }
    }
  }
});

test('groups on heavy atoms and dummies', () => {
  const eth = T.get('molecule', 'ethane').mol;
  const c2h5 = { Z: eth.Z.slice(0, 7), xyz: eth.xyz.slice(0, 21), bonds: eth.bonds.filter(b => b[0] < 7 && b[1] < 7) };
  let r = T.attachGroup(c2h5, 4, 'methyl');
  const A = r.added[0];
  ok(!r.map && !r.removed, 'no atoms removed for a heavy-atom attachment');
  for (const j of neighbours(r.mol, 4).filter(j => j !== A)) near(ang(r.mol.xyz, j, 4, A), 109.47, 3, 'open tetrahedral site');
  sane(r.mol, 'propane from C2H5 + methyl');
  r = T.attachGroup(T.get('molecule', 'pyridine').mol, 0, 'methyl');
  near(ang(r.mol.xyz, 1, 0, r.added[0]), 120, 3, 'N-methylpyridinium: in-plane site');
  const pt = T.buildCenter(EMPTY, { Z: 78, geometry: 'square-planar' }).mol;
  const pt3 = { Z: pt.Z.slice(0, 4), xyz: pt.xyz.slice(0, 12), bonds: pt.bonds.filter(b => b[1] < 4) };
  r = T.attachGroup(pt3, 0, 'methyl');
  const angs = [1, 2, 3].map(j => ang(r.mol.xyz, j, 0, r.added[0])).sort((a, b) => a - b);
  ok(Math.abs(angs[0] - 90) < 3 && Math.abs(angs[1] - 90) < 3 && angs[2] > 177, 'methyl fills the free square-planar site: ' + angs.map(a => a.toFixed(1)).join(', '));
  near(dist(r.mol.xyz, 0, r.added[0]), 2.12, 0.02, 'Pt–C');
  // a dummy (Z = 0) marker is replaced like an H
  const withX = { Z: [6, 1, 1, 1, 0], xyz: Float64Array.from([...T.get('molecule', 'methane').mol.xyz]), bonds: [[0, 1, 1], [0, 2, 1], [0, 3, 1], [0, 4, 1]] };
  r = T.attachGroup(withX, 4, 'OH');
  ok(r.mol.Z.length === 6 && r.mol.Z.indexOf(0) < 0, 'dummy replaced by the group');
  sane(r.mol, 'methanol from dummy');
  // chain growth stays zig-zag: butane by attaching ethyl to an ethane H
  r = T.attachGroup(eth, 5, 'ethyl');
  near(Math.abs(G.dihedral(r.mol.xyz, 0, 4, r.added[0], r.added[0] + 3)), 180, 15, 'butane anti');
  // nitro carries its formal charges along
  r = T.attachGroup(T.get('molecule', 'benzene').mol, 6, 'NO2');
  ok(r.mol.formalCharges && r.mol.formalCharges.reduce((s, q) => s + q, 0) === 0 && r.mol.formalCharges.some(q => q === 1), 'nitro formal charges N+ / O−');
});

test('cisplatin from square-planar PtCl4', () => {
  const pt = T.buildCenter(EMPTY, { Z: 'Pt', geometry: 'square-planar', ligand: 'Cl' });
  ok(pt.mol.Z.length === 5 && pt.added.length === 5, 'PtCl4 atoms');
  for (const a of ligandAngles(pt.mol, 0)) ok(Math.abs(a - 90) < 1 || a > 179, 'PtCl4 angle ' + a.toFixed(1));
  for (const d of donorsOf(pt.mol, 0)) near(dist(pt.mol.xyz, 0, d), 2.38, 0.01, 'Pt–Cl');
  const cl = donorsOf(pt.mol, 0), a = cl[0];
  const b = cl.find(j => j !== a && Math.abs(ang(pt.mol.xyz, a, 0, j) - 90) < 5);
  const r1 = T.replaceLigand(pt.mol, a, 'NH3');
  ok(r1.removed.length === 1 && r1.map[a] === -1, 'old Cl removed with map');
  const r2 = T.replaceLigand(r1.mol, r1.map[b], 'NH3');
  const m = r2.mol, N = m.Z.map((z, i) => (z === 7 ? i : -1)).filter(i => i >= 0), Cl = m.Z.map((z, i) => (z === 17 ? i : -1)).filter(i => i >= 0);
  ok(N.length === 2 && Cl.length === 2 && m.Z.length === 11, 'PtCl2(NH3)2');
  near(ang(m.xyz, N[0], 0, N[1]), 90, 3, 'cis N–Pt–N'); near(ang(m.xyz, Cl[0], 0, Cl[1]), 90, 3, 'cis Cl–Pt–Cl');
  for (const n of N) for (const c of Cl) { const t = ang(m.xyz, n, 0, c); ok(Math.abs(t - 90) < 3 || t > 177, 'N–Pt–Cl ' + t.toFixed(1)); }
  for (const n of N) {
    near(dist(m.xyz, 0, n), 2.07, 0.02, 'Pt–N');
    const hs = neighbours(m, n).filter(j => m.Z[j] === 1);
    ok(hs.length === 3, 'NH3 has three H');
    for (const h of hs) { near(dist(m.xyz, n, h), 1.01, 0.02, 'N–H'); near(ang(m.xyz, 0, n, h), 109.47, 6, 'Pt–N–H'); }
  }
  for (const c of Cl) near(dist(m.xyz, 0, c), 2.38, 0.02, 'Pt–Cl');
  ok(planeDev(m.xyz, [0, ...N, ...Cl]) < 0.05, 'square-planar core');
  sane(m, 'cisplatin');
});

test('ethylenediamine and other chelates on cis sites', () => {
  for (const [Zm, geo] of [['Pt', 'square-planar'], ['Co', 'octahedral']]) {
    const c = T.buildCenter(EMPTY, { Z: Zm, geometry: geo, ligand: 'Cl' }).mol;
    const cl = donorsOf(c, 0), a = cl[0], b = cl.find(j => j !== a && Math.abs(ang(c.xyz, a, 0, j) - 90) < 5);
    for (const L of ['en', 'bipy', 'ox', 'acac', 'phen']) {
      const tag = L + ' on ' + Zm + ' (' + geo + ')';
      const r = T.addBidentate(c, a, b, L), m = r.mol;
      const lig = T.get('ligand', L).mol;
      ok(r.removed.length === 2 && r.removed.includes(a) && r.removed.includes(b), tag + ': both old ligands removed');
      ok(m.Z.length === c.Z.length - 2 + lig.Z.length, tag + ': atom count');
      const D = [r.added[0], r.added[1]];
      ok(donorsOf(m, 0).includes(D[0]) && donorsOf(m, 0).includes(D[1]), tag + ': both donors bonded to the metal');
      for (const d of D) near(dist(m.xyz, 0, d), FF.idealBondLength(m.Z[0], m.Z[d], 1), 0.06, tag + ': M–donor');
      const bite = ang(m.xyz, D[0], 0, D[1]);
      ok(bite > 75 && bite < 100, tag + ': bite angle ' + bite.toFixed(1));
      // the donors sit on the two chosen (cis) sites
      const ua = [0, 1, 2].map(k => c.xyz[3 * a + k] - c.xyz[3 * 0 + k]), ub = [0, 1, 2].map(k => c.xyz[3 * b + k] - c.xyz[3 * 0 + k]);
      const dirOf = i => [0, 1, 2].map(k => m.xyz[3 * i + k] - m.xyz[k]);
      const angTo = (u, v) => Math.acos(Math.max(-1, Math.min(1, u.reduce((s, x, k) => s + x * v[k], 0) / (Math.hypot(...u) * Math.hypot(...v))))) * 180 / Math.PI;
      const fit = Math.min(Math.max(angTo(dirOf(D[0]), ua), angTo(dirOf(D[1]), ub)), Math.max(angTo(dirOf(D[0]), ub), angTo(dirOf(D[1]), ua)));
      ok(fit < 12, tag + ': donors on the chosen sites (within ' + fit.toFixed(1) + '°)');
      sane(m, tag, 0.1);
    }
  }
  const pt = T.buildCenter(EMPTY, { Z: 'Pt', geometry: 'square-planar' }).mol;
  let threw = '';
  try { T.addBidentate(pt, 1, 3, 'en'); } catch (e) { threw = e.message; }
  ok(/cis/.test(threw), 'trans sites rejected');
  threw = '';
  try { T.replaceLigand(pt, 1, 'en'); } catch (e) { threw = e.message; }
  ok(/addBidentate/.test(threw), 'replaceLigand refuses a bidentate');
  const en = T.addBidentate(pt, 1, 2, 'en').mol;
  const N = donorsOf(en, 0).filter(i => en.Z[i] === 7);
  const C = N.map(n => neighbours(en, n).find(j => en.Z[j] === 6));
  const tor = Math.abs(G.dihedral(en.xyz, N[0], C[0], C[1], N[1]));
  ok(tor > 35 && tor < 70, 'en chelate ring puckered (N–C–C–N ' + tor.toFixed(1) + '°)');
  // replacing one end of a chelate removes the whole ligand
  const r = T.replaceLigand(en, N[0], 'Cl');
  ok(r.mol.Z.filter(z => z === 7).length === 0 && r.removed.length === T.get('ligand', 'en').mol.Z.length, 'whole chelate removed');
  // tris-chelates
  for (const L of ['en', 'bipy']) {
    const m = T.buildCenter(EMPTY, { Z: 'Ru', geometry: 'octahedral', ligand: L }).mol;
    ok(donorsOf(m, 0).length === 6 && m.Z.indexOf(0) < 0, L + ': three chelates fill the octahedron');
    sane(m, '[Ru(' + L + ')3]', 0.1);
  }
});

test('octahedral Co(NH3)6, tetrahedral ZnCl4 and other centres', () => {
  let m = T.buildCenter(EMPTY, { Z: 27, geometry: 'octahedral', ligand: 'NH3', at: [2, 0, -1] }).mol;
  ok(m.Z.length === 25, 'Co(NH3)6 atoms');
  near(m.xyz[0], 2, 1e-9, 'metal placed at `at`');
  let cis = 0, trans = 0;
  for (const a of ligandAngles(m, 0)) { if (Math.abs(a - 90) < 3) cis++; else if (a > 177) trans++; }
  ok(cis === 12 && trans === 3, 'octahedral angles (' + cis + ' cis, ' + trans + ' trans)');
  for (const d of donorsOf(m, 0)) { ok(m.Z[d] === 7, 'N donors'); near(dist(m.xyz, 0, d), 1.97, 0.02, 'Co–N'); }
  sane(m, 'Co(NH3)6');
  m = T.buildCenter(EMPTY, { Z: 'Zn', geometry: 'tetrahedral', ligand: 'Cl' }).mol;
  for (const a of ligandAngles(m, 0)) near(a, 109.47, 2, 'Zn tetrahedral');
  for (const d of donorsOf(m, 0)) near(dist(m.xyz, 0, d), 2.24, 0.02, 'Zn–Cl');
  sane(m, 'ZnCl4');
  // every geometry with dummy markers; every monodentate ligand on an octahedron
  for (const g of T.geometries) {
    m = T.buildCenter(EMPTY, { Z: 'Fe', geometry: g.id, ligand: 'X' }).mol;
    const d = donorsOf(m, 0);
    ok(d.length === g.cn && d.every(i => m.Z[i] === 0), g.id + ': ' + g.cn + ' dummy markers');
    const sorted = ligandAngles(m, 0).sort((a, b) => a - b);
    const ideal = [];
    for (let p = 0; p < g.cn; p++) for (let q = p + 1; q < g.cn; q++) ideal.push(Math.acos(Math.max(-1, Math.min(1, g.vectors[p].reduce((s, x, k) => s + x * g.vectors[q][k], 0)))) * 180 / Math.PI);
    ideal.sort((a, b) => a - b);
    ok(sorted.every((a, k) => Math.abs(a - ideal[k]) < 0.5), g.id + ': marker angles ideal');
  }
  // every monodentate ligand: small ones six times on an octahedron, bulky ones (cone angle too large for
  // six around Co) four times on a tetrahedron
  const bulky = new Set(['PMe3', 'PPh3', 'py', 'DMSO']);
  for (const L of T.ligands.filter(l => l.denticity === 1 && l.id !== 'X')) {
    if (bulky.has(L.id)) {
      m = T.buildCenter(EMPTY, { Z: 'Pd', geometry: 'tetrahedral', ligand: L.id }).mol;
      for (const a of ligandAngles(m, 0)) near(a, 109.47, 6, 'Pd(' + L.id + ')4 tetrahedral');
      sane(m, 'Pd(' + L.id + ')4', 0.1);
      continue;
    }
    m = T.buildCenter(EMPTY, { Z: 'Co', geometry: 'octahedral', ligand: L.id }).mol;
    let c2 = 0, t2 = 0;
    for (const a of ligandAngles(m, 0)) { if (Math.abs(a - 90) < 5) c2++; else if (a > 175) t2++; }
    ok(c2 === 12 && t2 === 3, 'Co(' + L.id + ')6 stays octahedral (' + c2 + '/' + t2 + ')');
    sane(m, 'Co(' + L.id + ')6', 0.05);
  }
  // dummy marker → real ligand along the same direction
  const fe = T.buildCenter(EMPTY, { Z: 'Fe', geometry: 'octahedral', ligand: 'X' }).mol;
  const r = T.replaceLigand(fe, 3, 'CO');
  const C = r.added[0], u = [0, 1, 2].map(k => fe.xyz[9 + k] - fe.xyz[k]), v = [0, 1, 2].map(k => r.mol.xyz[3 * C + k] - r.mol.xyz[k]);
  near(Math.acos(u.reduce((s, x, k) => s + x * v[k], 0) / (Math.hypot(...u) * Math.hypot(...v))) * 180 / Math.PI, 0, 1, 'CO along the marker direction');
  near(dist(r.mol.xyz, 0, C), FF.idealBondLength(26, 6, 1), 0.02, 'Fe–C ideal');
  near(ang(r.mol.xyz, 0, C, C + 1), 180, 2, 'Fe–C–O linear');
  // existing atoms are kept and frozen when a centre is built next to them
  const host = T.insertMolecule(EMPTY, 'benzene', { at: [8, 0, 0] }).mol;
  const rc = T.buildCenter(host, { Z: 'Pd', geometry: 'square-planar', ligand: 'PH3' });
  ok(rc.added[0] === 12 && rc.mol.Z[12] === 46, 'metal appended after existing atoms');
  ok(Array.from(host.xyz).every((x, k) => rc.mol.xyz[k] === x), 'existing atoms untouched');
});

// ---------------------------------------------------------------------------------------------- summary
console.log('\nBlueMolecule template tests: ' + pass + ' passed, ' + fail + ' failed');
if (fail) {
  for (const f of failures) console.log('  FAIL ' + f);
  process.exit(1);
}
