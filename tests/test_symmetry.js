#!/usr/bin/env node
/* BlueMolecule point-group tests.  Run:  node tests/test_symmetry.js
 * All structures are built in code: named molecules, orbits of every axial and cubic group from
 * generic seed atoms (randomly rotated and translated), noisy copies, and 150-200 atom clusters. */
'use strict';
const path = require('path');
const ROOT = path.join(__dirname, '..');
require(path.join(ROOT, 'src', 'elements.js'));
require(path.join(ROOT, 'src', 'symmetry.js'));
const BM = globalThis.BM;

let pass = 0, fail = 0, current = '';
const failures = [];
function ok(cond, msg) { if (cond) pass++; else { fail++; failures.push('[' + current + '] ' + msg); } }
function eq(a, b, msg) { ok(a === b, msg + ': got ' + JSON.stringify(a) + ', expected ' + JSON.stringify(b)); }
function test(name, fn) {
  current = name;
  try { fn(); } catch (e) { fail++; failures.push('[' + name + '] threw: ' + (e && e.stack || e)); }
}
function throwsMsg(fn, re, msg) {
  try { fn(); ok(false, msg + ': no error thrown'); } catch (e) { ok(e instanceof Error && re.test(e.message), msg + ': message "' + e.message + '"'); }
}

// ------------------------------------------------------------------ helpers
function rng(seed) {   // mulberry32
  return function () {
    seed = (seed + 0x6D2B79F5) | 0;
    let t = Math.imul(seed ^ (seed >>> 15), 1 | seed);
    t = (t + Math.imul(t ^ (t >>> 7), 61 | t)) ^ t;
    return ((t ^ (t >>> 14)) >>> 0) / 4294967296;
  };
}
const R = rng(20261105);
const mol = atoms => ({ Z: atoms.map(a => BM.el.toZ(a[0])), xyz: [].concat(...atoms.map(a => a.slice(1))) });
const mm = (A, B) => [0, 1, 2].flatMap(i => [0, 1, 2].map(j => A[3 * i] * B[j] + A[3 * i + 1] * B[3 + j] + A[3 * i + 2] * B[6 + j]));
const mv = (M, v) => [M[0] * v[0] + M[1] * v[1] + M[2] * v[2], M[3] * v[0] + M[4] * v[1] + M[5] * v[2], M[6] * v[0] + M[7] * v[1] + M[8] * v[2]];
const unit = v => { const l = Math.hypot(...v); return v.map(x => x / l); };
const dot = (a, b) => a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
const crossv = (a, b) => [a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]];
function rotAxis(a, t) {
  const [x, y, z] = unit(a), c = Math.cos(t), s = Math.sin(t), C = 1 - c;
  return [c + x * x * C, x * y * C - z * s, x * z * C + y * s, y * x * C + z * s, c + y * y * C, y * z * C - x * s,
    z * x * C - y * s, z * y * C + x * s, c + z * z * C];
}
const refl = v => { const [x, y, z] = unit(v); return [1 - 2 * x * x, -2 * x * y, -2 * x * z, -2 * x * y, 1 - 2 * y * y, -2 * y * z, -2 * x * z, -2 * y * z, 1 - 2 * z * z]; };
const rotZ = t => rotAxis([0, 0, 1], t);
function randomRotation(r) {
  const u1 = r(), u2 = r() * 2 * Math.PI, u3 = r() * 2 * Math.PI;
  const a = Math.sqrt(1 - u1), b = Math.sqrt(u1);
  const [w, x, y, z] = [a * Math.sin(u2), a * Math.cos(u2), b * Math.sin(u3), b * Math.cos(u3)];
  return [1 - 2 * (y * y + z * z), 2 * (x * y - z * w), 2 * (x * z + y * w), 2 * (x * y + z * w), 1 - 2 * (x * x + z * z),
    2 * (y * z - x * w), 2 * (x * z - y * w), 2 * (y * z + x * w), 1 - 2 * (x * x + y * y)];
}
function moved(m, r) {
  const Rm = randomRotation(r), t = [10 * r() - 5, 10 * r() - 5, 10 * r() - 5];
  const xyz = [];
  for (let i = 0; i < m.Z.length; i++) {
    const p = mv(Rm, m.xyz.slice(3 * i, 3 * i + 3));
    xyz.push(p[0] + t[0], p[1] + t[1], p[2] + t[2]);
  }
  return { Z: m.Z.slice(), xyz };
}
const noisy = (m, amp, r) => ({ Z: m.Z.slice(), xyz: m.xyz.map(v => v + amp * (2 * r() - 1)) });
function com(Z, xyz) {
  let W = 0; const c = [0, 0, 0];
  Z.forEach((z, i) => { const w = BM.el.mass[z]; W += w; for (let k = 0; k < 3; k++) c[k] += w * xyz[3 * i + k]; });
  return c.map(v => v / W);
}
function minDist(xyz) {
  let d = Infinity;
  for (let i = 0; i < xyz.length / 3; i++) for (let j = i + 1; j < xyz.length / 3; j++) {
    d = Math.min(d, Math.hypot(xyz[3 * i] - xyz[3 * j], xyz[3 * i + 1] - xyz[3 * j + 1], xyz[3 * i + 2] - xyz[3 * j + 2]));
  }
  return d;
}
function closure(gens) {
  const ops = [[1, 0, 0, 0, 1, 0, 0, 0, 1]];
  for (let i = 0; i < ops.length; i++) for (const g of gens) {
    const p = mm(g, ops[i]);
    if (!ops.some(q => q.every((v, k) => Math.abs(v - p[k]) < 1e-8))) ops.push(p);
    if (ops.length > 500) throw new Error('closure runaway');
  }
  return ops;
}
function orbitMol(G, seeds) {
  const atoms = [];
  for (const [sym, p] of seeds) {
    const pts = [];
    for (const M of G) {
      const q = mv(M, p);
      if (!pts.some(s => Math.hypot(s[0] - q[0], s[1] - q[1], s[2] - q[2]) < 1e-6)) pts.push(q);
    }
    for (const q of pts) atoms.push([sym, ...q]);
  }
  return mol(atoms);
}
/** Operations reported for a structure must be orthogonal, map every atom onto its image and close. */
function checkOps(res, m, msg) {
  const c = res.center, n = m.Z.length;
  let worst = 0, orth = 0, perms = true;
  for (const op of res.operations) {
    const M = op.matrix;
    const MtM = mm([M[0], M[3], M[6], M[1], M[4], M[7], M[2], M[5], M[8]], M);
    orth = Math.max(orth, ...MtM.map((v, k) => Math.abs(v - (k % 4 === 0 ? 1 : 0))));
    const seen = new Set(op.perm);
    if (seen.size !== n) perms = false;
    for (let i = 0; i < n; i++) {
      const j = op.perm[i];
      if (m.Z[j] !== m.Z[i]) perms = false;
      const p = mv(M, [m.xyz[3 * i] - c[0], m.xyz[3 * i + 1] - c[1], m.xyz[3 * i + 2] - c[2]]);
      worst = Math.max(worst, Math.hypot(p[0] - m.xyz[3 * j] + c[0], p[1] - m.xyz[3 * j + 1] + c[1], p[2] - m.xyz[3 * j + 2] + c[2]));
    }
  }
  ok(orth < 1e-9, msg + ': orthogonal matrices (' + orth.toExponential(1) + ')');
  ok(perms, msg + ': permutations are bijections between like atoms');
  ok(worst <= 2 * res.tol + 1e-9, msg + ': operations map atoms onto atoms (worst ' + worst.toFixed(4) + ' Å)');
  if (Number.isFinite(res.order) && res.order <= 48) {
    const key = M => M.map(v => Math.round(v * 1e6) + 0).join(',');
    const set = new Set(res.operations.map(o => key(o.matrix)));
    let closed = true;
    for (const a of res.operations) for (const b of res.operations) if (!set.has(key(mm(a.matrix, b.matrix)))) closed = false;
    ok(closed, msg + ': operation set is closed under multiplication');
  }
  ok(res.operations[0].type === 'E', msg + ': identity listed first');
}
function detectIs(m, group, msg, opts) {
  const res = BM.sym.detect(m.Z, m.xyz, opts);
  eq(res.pointGroup, group, msg + ' point group');
  return res;
}
const counts = res => res.operations.reduce((a, o) => { const k = o.type + (o.n ? o.n : ''); a[k] = (a[k] || 0) + 1; return a; }, {});

// ------------------------------------------------------------------ molecules
const TH = Math.acos(-1 / 3);
const tetra = (d) => { const a = d / Math.sqrt(3); return [[a, a, a], [-a, -a, a], [-a, a, -a], [a, -a, -a]]; };
const M = {
  water: mol([['O', 0, 0, 0.1173], ['H', 0, 0.7572, -0.4692], ['H', 0, -0.7572, -0.4692]]),
  ammonia: (() => {
    const at = [['N', 0, 0, 0.1162]];
    for (let k = 0; k < 3; k++) at.push(['H', 0.9397 * Math.cos(2 * Math.PI * k / 3), 0.9397 * Math.sin(2 * Math.PI * k / 3), -0.2711]);
    return mol(at);
  })(),
  methane: mol([['C', 0, 0, 0], ...tetra(1.087).map(p => ['H', ...p])]),
  ch2cl2: mol([['C', 0, 0, 0.77], ['Cl', 1.47, 0, -0.24], ['Cl', -1.47, 0, -0.24], ['H', 0, 0.89, 1.37], ['H', 0, -0.89, 1.37]]),
  benzene: (() => {
    const at = [];
    for (let k = 0; k < 6; k++) { const t = k * Math.PI / 3; at.push(['C', 1.39 * Math.cos(t), 1.39 * Math.sin(t), 0], ['H', 2.47 * Math.cos(t), 2.47 * Math.sin(t), 0]); }
    return mol(at);
  })(),
  ethane: stag => {
    const at = [['C', 0, 0, 0.765], ['C', 0, 0, -0.765]];
    for (let k = 0; k < 3; k++) {
      const t1 = 2 * Math.PI * k / 3, t2 = t1 + (stag ? Math.PI / 3 : 0);
      at.push(['H', 1.02 * Math.cos(t1), 1.02 * Math.sin(t1), 1.16], ['H', 1.02 * Math.cos(t2), 1.02 * Math.sin(t2), -1.16]);
    }
    return mol(at);
  },
  ethylene: mol([['C', 0.667, 0, 0], ['C', -0.667, 0, 0], ['H', 1.232, 0.923, 0], ['H', 1.232, -0.923, 0], ['H', -1.232, 0.923, 0], ['H', -1.232, -0.923, 0]]),
  allene: mol([['C', 0, 0, 0], ['C', 0, 0, 1.308], ['C', 0, 0, -1.308], ['H', 0.928, 0, 1.869], ['H', -0.928, 0, 1.869],
    ['H', 0, 0.928, -1.869], ['H', 0, -0.928, -1.869]]),
  transDCE: mol([['C', -0.665, 0.0, 0], ['C', 0.665, 0, 0], ['Cl', -1.62, 1.45, 0], ['Cl', 1.62, -1.45, 0], ['H', -1.25, -0.92, 0], ['H', 1.25, 0.92, 0]]),
  h2o2: mol([['O', 0, 0.7375, 0], ['O', 0, -0.7375, 0], ['H', 0.8563, 0.9231, 0.2864], ['H', -0.8563, -0.9231, 0.2864]]),
  chfclbr: (() => {
    const t = tetra(1);
    return mol([['C', 0, 0, 0], ['H', ...t[0].map(v => v * 1.09)], ['F', ...t[1].map(v => v * 1.35)],
      ['Cl', ...t[2].map(v => v * 1.77)], ['Br', ...t[3].map(v => v * 1.94)]]);
  })(),
  co2: mol([['C', 0, 0, 0], ['O', 0, 0, 1.16], ['O', 0, 0, -1.16]]),
  hcn: mol([['H', 0, 0, -1.066], ['C', 0, 0, 0], ['N', 0, 0, 1.156]]),
  sf6: mol([['S', 0, 0, 0], ['F', 1.564, 0, 0], ['F', -1.564, 0, 0], ['F', 0, 1.564, 0], ['F', 0, -1.564, 0], ['F', 0, 0, 1.564], ['F', 0, 0, -1.564]]),
  cisPt: mol([['Pt', 0, 0, 0], ['Cl', 2.31, 0, 0], ['Cl', 0, 2.31, 0], ['N', -2.05, 0, 0], ['N', 0, -2.05, 0]]),
  transPt: mol([['Pt', 0, 0, 0], ['Cl', 2.31, 0, 0], ['Cl', -2.31, 0, 0], ['N', 0, 2.05, 0], ['N', 0, -2.05, 0]]),
  b12: (() => {
    const phi = (1 + Math.sqrt(5)) / 2, s = 1.7 / 2, at = [];
    for (const [x, y, z] of [[0, 1, phi], [0, -1, phi], [0, 1, -phi], [0, -1, -phi]]) at.push(['B', s * x, s * y, s * z], ['B', s * y, s * z, s * x], ['B', s * z, s * x, s * y]);
    return mol(at);
  })(),
  // 1,3,5,7-tetramethylcyclooctatetraene-like S4 skeleton: tub-shaped C8 ring with alternating substituents
  s4: (() => {
    const S4 = [0, 1, 0, -1, 0, 0, 0, 0, -1];   // (x, y, z) -> (y, -x, -z)
    return orbitMol(closure([S4]), [['C', [1.52, 0.48, 0.37]], ['C', [0.47, 1.49, -0.36]], ['H', [2.31, 1.02, 0.95]], ['C', [2.62, -0.71, 0.93]]]);
  })(),
  c4h: (() => {
    const gens = [rotZ(Math.PI / 2), [1, 0, 0, 0, 1, 0, 0, 0, -1]];
    return orbitMol(closure(gens), [['C', [1.4, 0.35, 0.72]], ['O', [2.6, 0.9, 1.1]], ['Fe', [0, 0, 0]]]);
  })(),
  // C60: no atom lies on a C5 or C3 axis
  c60: (() => {
    const phi = (1 + Math.sqrt(5)) / 2, at = [];
    for (const b of [[0, 1, 3 * phi], [1, 2 + phi, 2 * phi], [phi, 2, 2 * phi + 1]]) {
      for (const sx of [1, -1]) for (const sy of [1, -1]) for (const sz of [1, -1]) {
        const p = [b[0] * sx, b[1] * sy, b[2] * sz];
        if ((b[0] === 0 && sx < 0)) continue;
        for (const q of [p, [p[1], p[2], p[0]], [p[2], p[0], p[1]]]) at.push(['C', ...q.map(v => 0.7 * v)]);
      }
    }
    return mol(at);
  })(),
  cubane: (() => {
    const at = [];
    for (const sx of [1, -1]) for (const sy of [1, -1]) for (const sz of [1, -1]) at.push(['C', 0.78 * sx, 0.78 * sy, 0.78 * sz], ['H', 1.41 * sx, 1.41 * sy, 1.41 * sz]);
    return mol(at);
  })(),
  ferrocene: stag => {
    const at = [['Fe', 0, 0, 0]];
    for (const z of [1.66, -1.66]) for (let k = 0; k < 5; k++) {
      const t = 2 * Math.PI * k / 5 + (z < 0 && stag ? Math.PI / 5 : 0);
      at.push(['C', 1.21 * Math.cos(t), 1.21 * Math.sin(t), z], ['H', 2.28 * Math.cos(t), 2.28 * Math.sin(t), z * 1.03]);
    }
    return mol(at);
  },
  hc8h: (() => {
    const at = [['H', 0, 0, 0]];
    let z = 1.06, last = 0;
    for (let k = 0; k < 8; k++) { at.push(['C', 0, 0, z]); last = z; z += k % 2 ? 1.36 : 1.21; }
    at.push(['H', 0, 0, last + 1.06]);
    return mol(at);
  })(),
  cyclo18: (() => {   // polyynic cyclo[18]carbon: alternating bonds, D9h
    const at = [];
    for (let k = 0; k < 18; k++) { const t = 2 * Math.PI * k / 18 + (k % 2 ? 0.06 : 0); at.push(['C', 3.8 * Math.cos(t), 3.8 * Math.sin(t), 0]); }
    return mol(at);
  })()
};

test('named molecules', () => {
  const cases = [
    ['water', M.water, 'C2v', 4], ['ammonia', M.ammonia, 'C3v', 6], ['methane', M.methane, 'Td', 24],
    ['CH2Cl2', M.ch2cl2, 'C2v', 4], ['benzene', M.benzene, 'D6h', 24], ['staggered ethane', M.ethane(true), 'D3d', 12],
    ['eclipsed ethane', M.ethane(false), 'D3h', 12], ['ethylene', M.ethylene, 'D2h', 8], ['allene', M.allene, 'D2d', 8],
    ['trans-ClHC=CHCl', M.transDCE, 'C2h', 4], ['H2O2', M.h2o2, 'C2', 2], ['CHFClBr', M.chfclbr, 'C1', 1],
    ['CO2', M.co2, 'Dinfh', Infinity], ['HCN', M.hcn, 'Cinfv', Infinity], ['SF6', M.sf6, 'Oh', 48],
    ['cis-PtCl2N2', M.cisPt, 'C2v', 4], ['trans-PtCl2N2', M.transPt, 'D2h', 8], ['B12 icosahedron', M.b12, 'Ih', 120],
    ['S4 tub', M.s4, 'S4', 4], ['C4h', M.c4h, 'C4h', 8], ['C60', M.c60, 'Ih', 120], ['cubane', M.cubane, 'Oh', 48],
    ['eclipsed ferrocene', M.ferrocene(false), 'D5h', 20], ['staggered ferrocene', M.ferrocene(true), 'D5d', 20],
    ['HC8H', M.hc8h, 'Dinfh', Infinity], ['cyclo[18]carbon', M.cyclo18, 'D9h', 36]
  ];
  for (const [name, m, group, order] of cases) {
    const res = detectIs(m, group, name);
    eq(res.order, order, name + ' order');
    if (Number.isFinite(order)) eq(res.operations.length, order, name + ' operation count');
    checkOps(res, m, name);
    const moved1 = moved(m, R);
    eq(BM.sym.detect(moved1.Z, moved1.xyz).pointGroup, group, name + ' after random rotation + translation');
  }
});

test('operations, axes and classes', () => {
  const nh3 = BM.sym.detect(M.ammonia.Z, M.ammonia.xyz);
  const c = com(M.ammonia.Z, M.ammonia.xyz);
  ok(nh3.center.every((v, k) => Math.abs(v - c[k]) < 1e-12), 'centre of mass');
  ok(Math.abs(Math.abs(nh3.principalAxis[2]) - 1) < 1e-9, 'NH3 principal axis along the C3 axis: ' + nh3.principalAxis);
  eq(nh3.topType, 'symmetric', 'NH3 symmetric top');
  eq(JSON.stringify(nh3.classes.map(s => s.length).sort()), '[1,3]', 'NH3 equivalent atoms N | H3');
  const c3 = nh3.operations.filter(o => o.type === 'C');
  ok(c3.length === 2 && c3.every(o => o.n === 3 && o.k === 1 && Math.abs(Math.abs(o.axis[2]) - 1) < 1e-9), 'two C3 about z');
  eq(nh3.operations.filter(o => o.type === 'sigma').length, 3, 'three sigma_v');
  const bz = BM.sym.detect(M.benzene.Z, M.benzene.xyz);
  eq(JSON.stringify(counts(bz)), JSON.stringify({ E: 1, C6: 2, C3: 2, C2: 7, i: 1, S6: 2, S3: 2, sigma: 7 }), 'D6h class sizes');
  eq(JSON.stringify(bz.operations.filter(o => o.type === 'sigma').map(o => o.label)), JSON.stringify(['σh', 'σv', 'σv', 'σv', 'σd', 'σd', 'σd']), 'D6h mirror labels');
  eq(bz.operations.filter(o => /^C2[′″]$/.test(o.label)).length, 6, "D6h C2' and C2'' labels");
  const ch4 = BM.sym.detect(M.methane.Z, M.methane.xyz);
  eq(JSON.stringify(counts(ch4)), JSON.stringify({ E: 1, C3: 8, C2: 3, S4: 6, sigma: 6 }), 'Td class sizes');
  eq(ch4.topType, 'spherical', 'CH4 spherical top');
  const sf6 = BM.sym.detect(M.sf6.Z, M.sf6.xyz);
  eq(JSON.stringify(counts(sf6)), JSON.stringify({ E: 1, C4: 6, C3: 8, C2: 9, i: 1, S6: 8, S4: 6, sigma: 9 }), 'Oh class sizes');
  const b12 = BM.sym.detect(M.b12.Z, M.b12.xyz);
  eq(JSON.stringify(counts(b12)), JSON.stringify({ E: 1, C5: 24, C3: 20, C2: 15, i: 1, S10: 24, S6: 20, sigma: 15 }), 'Ih class sizes');
  ok(b12.operations.filter(o => o.type === 'C' && o.n === 5).every(o => o.k === 1 || o.k === 2), 'C5 and C5^2 powers');
  const s4 = BM.sym.detect(M.s4.Z, M.s4.xyz);
  ok(s4.operations.some(o => o.type === 'S' && o.n === 4) && !s4.operations.some(o => o.type === 'sigma' || o.type === 'i'), 'S4: S4 axis, no planes, no i');
  const w = BM.sym.detect(M.water.Z, M.water.xyz);
  eq(w.topType, 'asymmetric', 'water asymmetric top');
  ok(Math.abs(Math.abs(w.principalAxis[2]) - 1) < 1e-9, 'water C2 axis');
  eq(BM.sym.detect(M.co2.Z, M.co2.xyz).topType, 'linear', 'CO2 linear');
  const one = BM.sym.detect([26], [1, 2, 3]);
  ok(one.pointGroup === 'Kh' && one.order === Infinity && one.topType === 'atom', 'single atom: Kh');
  eq(BM.sym.detect([1, 1], [0, 0, 0, 0.74, 0, 0]).pointGroup, 'Dinfh', 'H2');
  eq(BM.sym.detect([], []).pointGroup, 'C1', 'empty structure');
  eq(BM.sym.detect([6, 0, 0], [0, 0, 0, 1, 0, 0, -1, 0, 0]).pointGroup, 'Dinfh', 'dummy atoms are matched with each other');
  throwsMsg(() => BM.sym.detect([1, 1], [0, 0, 0]), /matching size/, 'short coordinate array');
});

test('names and orders', () => {
  const nm = BM.sym.normalizeName;
  eq(nm('c2v'), 'C2v', 'c2v'); eq(nm('D*h'), 'Dinfh', 'D*h'); eq(nm('C∞v'), 'Cinfv', 'C∞v'); eq(nm('d3D'), 'D3d', 'd3D');
  eq(nm('S2'), 'Ci', 'S2 = Ci'); eq(nm('C3i'), 'S6', 'C3i = S6'); eq(nm('S3'), 'C3h', 'S3 = C3h'); eq(nm('TD'), 'Td', 'TD');
  eq(nm('C2d'), null, 'C2d is not a group'); eq(nm('X7'), null, 'garbage');
  const go = BM.sym.groupOrder;
  eq(go('D6h'), 24, 'D6h'); eq(go('D2d'), 8, 'D2d'); eq(go('S8'), 8, 'S8'); eq(go('C5v'), 10, 'C5v'); eq(go('Ih'), 120, 'Ih');
  eq(go('Dinfh'), Infinity, 'Dinfh');
});

// ------------------------------------------------------------------ every point group from orbits
const PHI = (1 + Math.sqrt(5)) / 2;
const SH = [1, 0, 0, 0, 1, 0, 0, 0, -1], INV = [-1, 0, 0, 0, -1, 0, 0, 0, -1];
const C2X = rotAxis([1, 0, 0], Math.PI), C3D = rotAxis([1, 1, 1], 2 * Math.PI / 3);
const T = [rotZ(Math.PI), C2X, C3D];
const GROUPS = {
  Cs: [SH], Ci: [INV], S4: [mm(SH, rotZ(Math.PI / 2))], S6: [mm(SH, rotZ(Math.PI / 3))], S8: [mm(SH, rotZ(Math.PI / 4))],
  T, Td: T.concat([refl([1, -1, 0])]), Th: T.concat([INV]), O: [rotZ(Math.PI / 2), C3D], Oh: [rotZ(Math.PI / 2), C3D, INV],
  I: T.concat([rotAxis([0, 1, PHI], 2 * Math.PI / 5)]), Ih: T.concat([rotAxis([0, 1, PHI], 2 * Math.PI / 5), INV])
};
for (const n of [2, 3, 4, 5, 6]) {
  const Cn = rotZ(2 * Math.PI / n);
  GROUPS['C' + n] = [Cn]; GROUPS['C' + n + 'v'] = [Cn, refl([0, 1, 0])]; GROUPS['C' + n + 'h'] = [Cn, SH];
  GROUPS['D' + n] = [Cn, C2X]; GROUPS['D' + n + 'h'] = [Cn, C2X, SH];
  GROUPS['D' + n + 'd'] = [Cn, C2X, refl([-Math.sin(Math.PI / (2 * n)), Math.cos(Math.PI / (2 * n)), 0])];
}
GROUPS.C8v = [rotZ(Math.PI / 4), refl([0, 1, 0])];

/** Orbits of `nseeds` generic atoms (different elements) under the group. With few seeds a structure can
 *  have extra symmetry by accident (azimuths aligned within tol), so clean orbits are checked at tol 0.02. */
function groupStructure(name, r, nseeds) {
  const G = closure(GROUPS[name]);
  const scale = Math.max(1, Math.sqrt(G.length / 5));
  for (let attempt = 0; attempt < 1000; attempt++) {
    const seeds = ['C', 'N', 'H', 'O', 'S'].slice(0, nseeds || 4).map((s, k) => {
      const rad = scale * (1.2 + 0.9 * k + 0.5 * r()), u = unit([r() - 0.5, r() - 0.5, r() - 0.5]);
      return [s, u.map(v => v * rad)];
    });
    const m = orbitMol(G, seeds);
    if (minDist(m.xyz) > 0.75) return { m, order: G.length };
  }
  throw new Error('could not place generic atoms for ' + name);
}

test('all point groups from generic orbits', () => {
  for (const name of Object.keys(GROUPS)) {
    const { m, order } = groupStructure(name, R);
    eq(order, BM.sym.groupOrder(name), name + ' generator check');
    const mo = moved(m, R);
    const res = detectIs(mo, name, name + ' (' + m.Z.length + ' atoms)', { tol: 0.02 });
    eq(res.order, order, name + ' order');
    checkOps(res, mo, name);
    eq(res.classes.length, 4, name + ': four orbits of equivalent atoms');
  }
});

// ------------------------------------------------------------------ noise and symmetrization
function checkSymmetrized(m, group, amp, msg) {
  const nz = noisy(moved(m, R), amp, R);
  const res = detectIs(nz, group, msg + ' with ±' + amp + ' Å noise');
  const s = BM.sym.symmetrize(nz.Z, nz.xyz);
  eq(s.pointGroup, group, msg + ' symmetrize group');
  ok(s.xyz instanceof Float64Array && s.xyz.length === nz.xyz.length, msg + ' output array');
  const tight = BM.sym.detect(nz.Z, s.xyz, { tol: 1e-6 });
  eq(tight.pointGroup, group, msg + ' symmetrized structure is exactly ' + group + ' (tol 1e-6)');
  ok(tight.maxDeviation < 1e-9, msg + ' residual asymmetry ' + tight.maxDeviation.toExponential(2));
  const c0 = com(nz.Z, nz.xyz), c1 = com(nz.Z, Array.from(s.xyz));
  ok(c0.every((v, k) => Math.abs(v - c1[k]) < 1e-9), msg + ' centre of mass kept');
  ok(s.maxShift < 3 * amp && s.maxShift > 0 && s.rmsShift <= s.maxShift, msg + ' shifts ' + s.maxShift.toFixed(4) + ' / ' + s.rmsShift.toFixed(4));
  let far = 0;
  for (let i = 0; i < nz.xyz.length; i++) far = Math.max(far, Math.abs(s.xyz[i] - nz.xyz[i]));
  ok(far < 3 * amp, msg + ' orientation kept (largest coordinate change ' + far.toFixed(4) + ')');
  return res;
}

test('0.02 Å noise: detection and symmetrization', () => {
  const cases = [['water', M.water, 'C2v'], ['ammonia', M.ammonia, 'C3v'], ['methane', M.methane, 'Td'], ['CH2Cl2', M.ch2cl2, 'C2v'],
    ['benzene', M.benzene, 'D6h'], ['staggered ethane', M.ethane(true), 'D3d'], ['eclipsed ethane', M.ethane(false), 'D3h'],
    ['ethylene', M.ethylene, 'D2h'], ['allene', M.allene, 'D2d'], ['trans-ClHC=CHCl', M.transDCE, 'C2h'], ['H2O2', M.h2o2, 'C2'],
    ['CO2', M.co2, 'Dinfh'], ['HCN', M.hcn, 'Cinfv'], ['SF6', M.sf6, 'Oh'], ['cis-PtCl2N2', M.cisPt, 'C2v'],
    ['trans-PtCl2N2', M.transPt, 'D2h'], ['B12', M.b12, 'Ih'], ['S4 tub', M.s4, 'S4'], ['C4h', M.c4h, 'C4h'],
    ['C60', M.c60, 'Ih'], ['staggered ferrocene', M.ferrocene(true), 'D5d'], ['HC8H', M.hc8h, 'Dinfh']];
  for (const [name, m, g] of cases) checkSymmetrized(m, g, 0.02, name);
  for (const g of ['C3h', 'D5d', 'S6', 'T', 'Th', 'O', 'I', 'C5v', 'D4h']) checkSymmetrized(groupStructure(g, R, 5).m, g, 0.02, g + ' orbit');
  // Linear: the symmetrized molecule is exactly straight.
  const nz = noisy(moved(M.co2, R), 0.02, R);
  const s = BM.sym.symmetrize(nz.Z, nz.xyz);
  const p = i => Array.from(s.xyz.slice(3 * i, 3 * i + 3));
  const u = unit(p(1).map((v, k) => v - p(2)[k])), w = p(0).map((v, k) => v - p(2)[k]);
  ok(Math.hypot(...crossv(u, w)) < 1e-12, 'CO2 symmetrized onto a straight line');
  ok(Math.abs(Math.hypot(...p(1).map((v, k) => v - p(0)[k])) - Math.hypot(...p(2).map((v, k) => v - p(0)[k]))) < 1e-12, 'CO2 equal C=O bonds');
  eq(BM.sym.detect(nz.Z, s.xyz, { tol: 1e-9 }).pointGroup, 'Dinfh', 'symmetrized CO2 is linear even at tol 1e-9');
});

test('symmetrize with a requested point group', () => {
  const nz = noisy(moved(M.benzene, R), 0.02, R);
  const sub = BM.sym.symmetrize(nz.Z, nz.xyz, { pointGroup: 'c2v' });
  eq(sub.pointGroup, 'C2v', 'subgroup C2v of D6h');
  eq(sub.order, 4, 'subgroup order');
  const g1 = BM.sym.detect(nz.Z, sub.xyz, { tol: 1e-6 }).pointGroup;
  ok(g1 === 'C2v', 'C2v-symmetrized noisy benzene keeps only C2v exactly: ' + g1);
  const d2h = BM.sym.symmetrize(nz.Z, nz.xyz, { pointGroup: 'D2h' });
  eq(BM.sym.detect(nz.Z, d2h.xyz, { tol: 1e-6 }).pointGroup, 'D2h', 'D2h subgroup (three generators)');
  const d3h = BM.sym.symmetrize(nz.Z, nz.xyz, { pointGroup: 'D3h' });
  eq(BM.sym.detect(nz.Z, d3h.xyz, { tol: 1e-6 }).pointGroup, 'D3h', 'D3h subgroup');
  // Larger distortions than tol: the requested group is found by relaxing the tolerance.
  const bad = noisy(moved(M.benzene, R), 0.15, R);
  ok(BM.sym.detect(bad.Z, bad.xyz).pointGroup !== 'D6h', 'distorted benzene is not D6h at tol 0.1');
  const fix = BM.sym.symmetrize(bad.Z, bad.xyz, { pointGroup: 'D6h' });
  eq(fix.pointGroup, 'D6h', 'forced D6h');
  eq(BM.sym.detect(bad.Z, fix.xyz, { tol: 1e-6 }).pointGroup, 'D6h', 'forced D6h is exact');
  ok(fix.maxShift > 0.1, 'large shift reported: ' + fix.maxShift.toFixed(3));
  const lin = noisy(M.co2, 0.02, R);
  eq(BM.sym.symmetrize(lin.Z, lin.xyz, { pointGroup: 'Cinfv' }).pointGroup, 'Cinfv', 'Cinfv subgroup of Dinfh');
  throwsMsg(() => BM.sym.symmetrize(M.water.Z, M.water.xyz, { pointGroup: 'Oh' }), /not close to Oh/, 'impossible group');
  throwsMsg(() => BM.sym.symmetrize(M.water.Z, M.water.xyz, { pointGroup: 'Q7' }), /Unknown point group/, 'unknown group');
  const same = BM.sym.symmetrize(M.water.Z, M.water.xyz, { pointGroup: 'C2v' });
  ok(same.maxShift < 1e-12, 'already symmetric water is unchanged');
});

test('tolerance behaviour', () => {
  const w = noisy(M.water, 0.0, R);
  w.xyz[4] += 0.05;   // one O-H longer by 0.05 Å
  eq(BM.sym.detect(w.Z, w.xyz).pointGroup, 'C2v', 'within tol 0.1');
  eq(BM.sym.detect(w.Z, w.xyz, { tol: 0.01 }).pointGroup, 'Cs', 'tight tol 0.01 sees only the molecular plane');
  const bent = mol([['O', 0, 0, 0], ['C', 1.16, 0, 0], ['O', 2.32, 0.03, 0]]);
  eq(BM.sym.detect(bent.Z, bent.xyz).pointGroup, 'Dinfh', 'near-linear within tol');
  eq(BM.sym.detect(bent.Z, bent.xyz, { tol: 0.005 }).pointGroup, 'C2v', 'bent at tight tol');
  // A tolerance comparable to the spacing of like atoms is clamped (and the search stays fast).
  const cloud = [], r3 = rng(9);
  while (cloud.length < 150) {
    const p = [12 * r3(), 12 * r3(), 12 * r3()];
    if (cloud.every(a => Math.hypot(a[1] - p[0], a[2] - p[1], a[3] - p[2]) > 1.0)) cloud.push(['C', ...p]);
  }
  const cm = mol(cloud), t0 = Date.now();
  const big = BM.sym.detect(cm.Z, cm.xyz, { tol: 0.8 });
  ok(big.tol < 0.8 && big.tol > 0.3, 'tolerance clamped below half the like-atom spacing: ' + big.tol.toFixed(3));
  eq(big.pointGroup, 'C1', 'random carbon cloud');
  throwsMsg(() => BM.sym.symmetrize(cm.Z, cm.xyz, { pointGroup: 'Oh' }), /not close to Oh/, 'impossible group on a 150-atom cloud');
  ok(Date.now() - t0 < 5000, 'large-tolerance searches finish quickly (' + (Date.now() - t0) + ' ms)');
});

// ------------------------------------------------------------------ performance
test('performance (150-200 atoms)', () => {
  const big = [];
  // I: 60 + 60 + 30 = 150 atoms.  Td: 6 x 24 + 6 = 150.  D3h: 12 x 12 + 2 + 3 + 1 = 150.  Oh: 4 x 48 + 8 = 200.
  const gI = closure(GROUPS.I);
  big.push(['I (150)', orbitMol(gI, [['C', [1.1, 2.3, 4.1]], ['C', [2.9, 0.7, 4.4]], ['H', [0, 0, 6.2]]]), 'I']);
  const gTd = closure(GROUPS.Td);
  const tdSeeds = [];
  for (let k = 0; k < 6; k++) tdSeeds.push([k % 2 ? 'C' : 'H', [1.5 + 0.9 * k, 0.6 + 0.37 * k, 0.21 + 0.53 * k]]);
  tdSeeds.push(['O', [0, 0, 1.3]]);
  big.push(['Td (150)', orbitMol(gTd, tdSeeds), 'Td']);
  const gD3h = closure(GROUPS.D3h);
  const d3Seeds = [['Fe', [0, 0, 0]], ['C', [0, 0, 2.6]], ['N', [3.1, 0, 0]]];
  for (let k = 0; k < 12; k++) d3Seeds.push([['C', 'H', 'O'][k % 3], [2.0 + 0.8 * k, 0.9 + 0.11 * k, 0.7 + 0.43 * k]]);
  big.push(['D3h (150)', orbitMol(gD3h, d3Seeds), 'D3h']);
  const gOh = closure(GROUPS.Oh);
  const ohSeeds = [['N', [1.9, 1.9, 1.9]]];
  for (let k = 0; k < 4; k++) ohSeeds.push([k % 2 ? 'C' : 'H', [2.2 + 1.1 * k, 0.8 + 0.29 * k, 0.35 + 0.61 * k]]);
  big.push(['Oh (200)', orbitMol(gOh, ohSeeds), 'Oh']);
  const rnd = [];
  const r2 = rng(7);
  while (rnd.length < 150) {
    const p = [12 * r2(), 12 * r2(), 12 * r2()];
    if (rnd.every(a => Math.hypot(a[1] - p[0], a[2] - p[1], a[3] - p[2]) > 1.0)) rnd.push([['C', 'H', 'N', 'O'][rnd.length % 4], ...p]);
  }
  big.push(['C1 random (150)', mol(rnd), 'C1']);
  big.push(['I (150) noisy', noisy(moved(big[0][1], R), 0.02, R), 'I']);
  const timing = [];
  for (const [name, m, g] of big) {
    const t0 = process.hrtime.bigint();
    const res = BM.sym.detect(m.Z, m.xyz);
    const ms = Number(process.hrtime.bigint() - t0) / 1e6;
    eq(res.pointGroup, g, name);
    ok(ms < 1000, name + ' detected in ' + ms.toFixed(0) + ' ms (< 1 s)');
    const t1 = process.hrtime.bigint();
    BM.sym.symmetrize(m.Z, m.xyz);
    const ms2 = Number(process.hrtime.bigint() - t1) / 1e6;
    ok(ms2 < 2000, name + ' symmetrized in ' + ms2.toFixed(0) + ' ms');
    timing.push(name + ' ' + m.Z.length + ' atoms: ' + ms.toFixed(0) + ' ms');
  }
  console.log('  timings: ' + timing.join('; '));
});

// ---------------------------------------------------------------------- summary
console.log('\nBlueMolecule symmetry tests: ' + pass + ' passed, ' + fail + ' failed');
if (fail) {
  for (const f of failures) console.log('  FAIL ' + f);
  process.exit(1);
}
