#!/usr/bin/env node
/* BlueMolecule clean-up force-field tests.  Run:  node tests/test_ff.js */
'use strict';
const path = require('path');
const ROOT = path.join(__dirname, '..');
require(path.join(ROOT, 'src', 'elements.js'));
require(path.join(ROOT, 'src', 'geom.js'));
require(path.join(ROOT, 'src', 'ff.js'));
const BM = globalThis.BM, FF = BM.ff, G = BM.geom;

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
let seed = 12345;
const rnd = () => { seed = (seed * 16807) % 2147483647; return seed / 2147483647; };
function jitter(xyz, amp) { const x = Float64Array.from(xyz); for (let i = 0; i < x.length; i++) x[i] += (rnd() - 0.5) * 2 * amp; return x; }
const dist = (x, i, j) => G.dist(x, i, j), ang = (x, i, j, k) => G.angle(x, i, j, k), dih = (x, i, j, k, l) => G.dihedral(x, i, j, k, l);
/** Largest distance of the listed atoms from their least-squares plane (3x3 Jacobi eigen-solver). */
function planeDev(X, idx) {
  const n = idx.length, c = [0, 0, 0];
  for (const i of idx) for (let k = 0; k < 3; k++) c[k] += X[3 * i + k] / n;
  const a = [[0, 0, 0], [0, 0, 0], [0, 0, 0]], v = [[1, 0, 0], [0, 1, 0], [0, 0, 1]];
  for (const i of idx) { const d = [0, 1, 2].map(k => X[3 * i + k] - c[k]); for (let p = 0; p < 3; p++) for (let q = 0; q < 3; q++) a[p][q] += d[p] * d[q]; }
  for (let sweep = 0; sweep < 60; sweep++) {
    for (let p = 0; p < 3; p++) for (let q = p + 1; q < 3; q++) {
      if (Math.abs(a[p][q]) < 1e-300) continue;
      const th = (a[q][q] - a[p][p]) / (2 * a[p][q]);
      const t = (th >= 0 ? 1 : -1) / (Math.abs(th) + Math.sqrt(th * th + 1)), cs = 1 / Math.sqrt(t * t + 1), sn = t * cs;
      for (let k = 0; k < 3; k++) { const x1 = a[k][p], x2 = a[k][q]; a[k][p] = cs * x1 - sn * x2; a[k][q] = sn * x1 + cs * x2; }
      for (let k = 0; k < 3; k++) { const x1 = a[p][k], x2 = a[q][k]; a[p][k] = cs * x1 - sn * x2; a[q][k] = sn * x1 + cs * x2; }
      for (let k = 0; k < 3; k++) { const x1 = v[k][p], x2 = v[k][q]; v[k][p] = cs * x1 - sn * x2; v[k][q] = sn * x1 + cs * x2; }
    }
  }
  let m = 0;
  for (let k = 1; k < 3; k++) if (a[k][k] < a[m][m]) m = k;
  const nv = [v[0][m], v[1][m], v[2][m]];
  let mx = 0;
  for (const i of idx) mx = Math.max(mx, Math.abs((X[3 * i] - c[0]) * nv[0] + (X[3 * i + 1] - c[1]) * nv[1] + (X[3 * i + 2] - c[2]) * nv[2]));
  return mx;
}
const range = n => [...Array(n).keys()];
/** Central finite-difference check of the analytic gradient. Returns the largest absolute error. */
function gradError(mol, opts) {
  const ctx = FF.setup(mol, opts);
  const x = Float64Array.from(mol.xyz), ev = FF.evaluate(ctx, x);
  let err = 0;
  for (let k = 0; k < x.length; k++) {
    if (ev.grad[k] === 0 && !ctx.free[Math.floor(k / 3)]) continue;
    const h = 1e-5, xp = Float64Array.from(x), xm = Float64Array.from(x);
    xp[k] += h; xm[k] -= h;
    const num = (FF.evaluate(ctx, xp).energy - FF.evaluate(ctx, xm).energy) / (2 * h);
    err = Math.max(err, Math.abs(num - ev.grad[k]) / Math.max(1, Math.abs(num)));
  }
  return { err, terms: ev.terms };
}

// ---------------------------------------------------------------------------------------------- molecules
function water() { return { Z: [8, 1, 1], xyz: Float64Array.from([0, 0, 0, 1.12, 0.05, 0, -0.62, 1.0, 0.15]), bonds: [[0, 1, 1], [0, 2, 1]] }; }
function benzene() {
  const Z = [], xyz = [], bonds = [];
  for (let k = 0; k < 6; k++) { const a = k * Math.PI / 3; Z.push(6); xyz.push(1.39 * Math.cos(a), 1.39 * Math.sin(a), 0); }
  for (let k = 0; k < 6; k++) { const a = k * Math.PI / 3; Z.push(1); xyz.push(2.48 * Math.cos(a), 2.48 * Math.sin(a), 0); bonds.push([k, (k + 1) % 6, 1.5], [k, k + 6, 1]); }
  return { Z, xyz: Float64Array.from(xyz), bonds };
}
function cyclohexaneFlat(flatten) {
  const b = 1.52, th = 109.47 * Math.PI / 180, R = 2 * b * Math.sin(th / 2) / Math.sqrt(3), h = Math.sqrt((b * b - R * R) / 4);
  const Z = [], xyz = [], bonds = [];
  for (let k = 0; k < 6; k++) { const a = k * Math.PI / 3; Z.push(6); xyz.push(R * Math.cos(a), R * Math.sin(a), (k % 2 ? h : -h) * flatten); }
  for (let k = 0; k < 6; k++) {
    const a = k * Math.PI / 3, s = k % 2 ? 1 : -1;
    Z.push(1); xyz.push((R + 0.4) * Math.cos(a), (R + 0.4) * Math.sin(a), s * 1.0); bonds.push([k, Z.length - 1, 1]);
    Z.push(1); xyz.push((R + 1.0) * Math.cos(a), (R + 1.0) * Math.sin(a), -s * 0.35); bonds.push([k, Z.length - 1, 1]);
    bonds.push([k, (k + 1) % 6, 1]);
  }
  return { Z, xyz: Float64Array.from(xyz), bonds };
}
/** Ethane with an H–C–C–H dihedral of `phi0` degrees (0 = eclipsed). */
function ethane(phi0) {
  const Z = [6, 6], xyz = [0, 0, 0, 1.53, 0, 0], bonds = [[0, 1, 1]];
  for (let c = 0; c < 2; c++) for (let k = 0; k < 3; k++) {
    const a = (k * 120 + (c ? phi0 : 0)) * Math.PI / 180, sx = c ? 1.53 + 0.36 : -0.36;
    Z.push(1); xyz.push(sx, 1.03 * Math.cos(a), 1.03 * Math.sin(a)); bonds.push([c, Z.length - 1, 1]);
  }
  return { Z, xyz: Float64Array.from(xyz), bonds };
}
function ethylene(twist) {
  const t = twist * Math.PI / 180;
  return {
    Z: [6, 6, 1, 1, 1, 1],
    xyz: Float64Array.from([0, 0, 0, 1.36, 0, 0, -0.55, 0.92, 0, -0.55, -0.92, 0, 1.91, 0.92 * Math.cos(t), 0.92 * Math.sin(t), 1.91, -0.92 * Math.cos(t), -0.92 * Math.sin(t)]),
    bonds: [[0, 1, 2], [0, 2, 1], [0, 3, 1], [1, 4, 1], [1, 5, 1]]
  };
}
function acetyleneBent() {
  return { Z: [1, 6, 6, 1], xyz: Float64Array.from([-1.0, 0.55, 0, 0, 0, 0, 1.22, 0, 0.1, 2.2, -0.4, -0.2]), bonds: [[0, 1, 1], [1, 2, 3], [2, 3, 1]] };
}
function cisplatinDistorted() {
  const Z = [78, 17, 17, 7, 1, 1, 1, 7, 1, 1, 1];
  const xyz = [0, 0, 0, -1.8, 1.4, 0.3, -1.5, -1.9, -0.2, 1.3, 1.6, 0.2, 2.0, 1.4, 0.9, 2.1, 1.3, -0.8, 1.0, 2.4, 0.1,
    1.6, -1.2, -0.3, 1.3, -2.1, 0.8, 1.4, -2.0, -0.8, 2.4, -1.0, 0.0];
  const bonds = [[0, 1, 1], [0, 2, 1], [0, 3, 1], [0, 7, 1], [3, 4, 1], [3, 5, 1], [3, 6, 1], [7, 8, 1], [7, 9, 1], [7, 10, 1]];
  return { Z, xyz: Float64Array.from(xyz), bonds };
}
/** A metal (Z) with ligand atoms (Zl, distance d) along the given unit directions, plus H3 on N donors. */
function complex(Zm, dirs, ligZ, d) {
  const Z = [Zm], xyz = [0, 0, 0], bonds = [];
  dirs.forEach((v, k) => {
    const L = Z.length, zl = ligZ[k];
    Z.push(zl); xyz.push(v[0] * d[k], v[1] * d[k], v[2] * d[k]); bonds.push([0, L, 1]);
    if (zl === 7) {
      const u = v, p = Math.abs(u[0]) < 0.9 ? [1, 0, 0] : [0, 1, 0];
      const cr = (a, b) => [a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]];
      const nrm = a => { const l = Math.hypot(...a); return a.map(c => c / l); };
      const w1 = nrm(cr(u, p)), w2 = cr(u, w1);
      for (let h = 0; h < 3; h++) {
        const a = h * 2 * Math.PI / 3, cs = Math.cos(a), sn = Math.sin(a);
        const dir = [0, 1, 2].map(c => 0.34 * u[c] + 0.94 * (cs * w1[c] + sn * w2[c]));
        Z.push(1); xyz.push(xyz[3 * L] + dir[0] * 1.01, xyz[3 * L + 1] + dir[1] * 1.01, xyz[3 * L + 2] + dir[2] * 1.01); bonds.push([L, Z.length - 1, 1]);
      }
    }
  });
  return { Z, xyz: Float64Array.from(xyz), bonds };
}
const TET = [[1, 1, 1], [-1, -1, 1], [-1, 1, -1], [1, -1, -1]].map(v => v.map(c => c / Math.sqrt(3)));
const SQP = [[1, 0, 0], [0, 1, 0], [-1, 0, 0], [0, -1, 0]];
const OCT = [[1, 0, 0], [0, 1, 0], [-1, 0, 0], [0, -1, 0], [0, 0, 1], [0, 0, -1]];
/** All L–M–L angles of the metal (atom 0) in a complex() molecule. */
function ligandAngles(x, mol) {
  const d = mol.bonds.filter(b => b[0] === 0).map(b => b[1]), a = [];
  for (let p = 0; p < d.length; p++) for (let q = p + 1; q < d.length; q++) a.push(ang(x, d[p], 0, d[q]));
  return a;
}
/** N-methylacetamide-like amide plus a phenyl ring, an ether and a close nonbonded contact. */
function mixedMolecule() {
  // CH3-C(=O)-NH-CH3 (amide), built roughly, then distorted
  const Z = [6, 6, 8, 7, 6, 1, 1, 1, 1, 1, 1, 1, 1];
  const xyz = [0, 0, 0, 1.5, 0, 0, 2.1, 1.05, 0, 2.2, -1.2, 0.1, 3.65, -1.25, 0, -0.4, 1.0, 0.2, -0.35, -0.5, 0.9, -0.4, -0.55, -0.85,
    1.7, -2.05, 0.3, 4.0, -2.2, 0.35, 4.05, -0.5, 0.7, 4.0, -1.0, -1.0, 2.0, 1.9, 0.6];
  const bonds = [[0, 1, 1], [1, 2, 2], [1, 3, 1], [3, 4, 1], [0, 5, 1], [0, 6, 1], [0, 7, 1], [3, 8, 1], [4, 9, 1], [4, 10, 1], [4, 11, 1]];
  // atom 12: an isolated H placed close to the carbonyl O (nonbonded contact)
  return { Z, xyz: jitter(xyz, 0.08), bonds };
}

console.log('BlueMolecule force-field tests');

// ---------------------------------------------------------------------------------------------- tests
test('shared ideal bond lengths', () => {
  near(FF.idealBondLength(6, 6, 1), 1.52, 1e-6, 'C–C');
  near(FF.idealBondLength(6, 6, 1.5), 1.52 * 0.915, 1e-6, 'C–C aromatic');
  near(FF.idealBondLength(6, 6, 2), 1.52 * 0.875, 1e-6, 'C=C');
  near(FF.idealBondLength(6, 6, 3), 1.52 * 0.79, 1e-6, 'C≡C');
  near(FF.idealBondLength(6, 1, 1), 1.09, 1e-12, 'C–H'); near(FF.idealBondLength(1, 8, 1), 0.96, 1e-12, 'O–H');
  near(FF.idealBondLength(7, 1, 1), 1.01, 1e-12, 'N–H'); near(FF.idealBondLength(1, 1, 1), 0.74, 1e-12, 'H–H');
  near(FF.idealBondLength(16, 1, 1), 1.34, 1e-12, 'S–H'); near(FF.idealBondLength(15, 1, 1), 1.42, 1e-12, 'P–H');
  near(FF.idealBondLength(14, 1, 1), 1.48, 1e-12, 'Si–H'); near(FF.idealBondLength(5, 1, 1), 1.19, 1e-12, 'B–H');
  near(FF.idealBondLength(78, 1, 1), BM.el.covalent[78] + 0.31, 1e-6, 'Pt–H (other X–H)');
  near(FF.idealBondLength(78, 7, 1), 2.07, 1e-6, 'Pt–N'); near(FF.idealBondLength(78, 17, 1), 2.38, 1e-6, 'Pt–Cl');
});

test('water', () => {
  const r = FF.optimize(water());
  ok(r.converged, 'converged');
  near(ang(r.xyz, 1, 0, 2), 104.5, 0.5, 'H–O–H'); near(dist(r.xyz, 0, 1), 0.96, 0.005, 'O–H'); near(dist(r.xyz, 0, 2), 0.96, 0.005, 'O–H\'');
});

test('benzene', () => {
  const m = benzene(); m.xyz = jitter(m.xyz, 0.2);
  ok(planeDev(m.xyz, range(12)) > 0.1, 'start is distorted');
  const r = FF.optimize(m);
  ok(r.converged, 'converged');
  ok(planeDev(r.xyz, range(12)) < 0.02, 'planar: max deviation ' + planeDev(r.xyz, range(12)).toFixed(4) + ' Å');
  for (let k = 0; k < 6; k++) { const d = dist(r.xyz, k, (k + 1) % 6); ok(d >= 1.38 && d <= 1.42, 'C–C ' + d.toFixed(4)); }
  for (let k = 0; k < 6; k++) near(ang(r.xyz, k, (k + 1) % 6, (k + 2) % 6), 120, 1, 'C–C–C');
});

test('cyclohexane chair', () => {
  const m = cyclohexaneFlat(0.25);
  const r = FF.optimize(m);
  ok(r.converged, 'converged');
  for (let k = 0; k < 6; k++) {
    const d = dih(r.xyz, k, (k + 1) % 6, (k + 2) % 6, (k + 3) % 6);
    ok(Math.abs(Math.abs(d) - 55) <= 8, 'ring dihedral ' + d.toFixed(1));
    ok(k === 0 || Math.sign(d) !== Math.sign(dih(r.xyz, k - 1, k, (k + 1) % 6, (k + 2) % 6)), 'alternating signs (chair)');
  }
});

test('ethane staggered, ethylene planar, acetylene linear', () => {
  let r = FF.optimize(ethane(20));
  const phi = Math.abs(dih(r.xyz, 2, 0, 1, 5));
  ok(r.converged, 'ethane converged');
  ok([60, 180].some(t => Math.abs(phi - t) < 3), 'ethane H–C–C–H staggered: ' + phi.toFixed(1));
  near(dist(r.xyz, 0, 1), 1.52, 0.01, 'ethane C–C');
  r = FF.optimize(ethylene(25));
  ok(r.converged, 'ethylene converged');
  ok(planeDev(r.xyz, range(6)) < 0.01, 'ethylene planar');
  near(Math.abs(dih(r.xyz, 2, 0, 1, 4)), 0, 1, 'cis H–C=C–H'); near(Math.abs(dih(r.xyz, 2, 0, 1, 5)), 180, 1, 'trans H–C=C–H');
  near(dist(r.xyz, 0, 1), 1.33, 0.01, 'C=C');
  r = FF.optimize(acetyleneBent());
  ok(r.converged, 'acetylene converged');
  ok(ang(r.xyz, 0, 1, 2) > 178.5 && ang(r.xyz, 1, 2, 3) > 178.5, 'linear: ' + ang(r.xyz, 0, 1, 2).toFixed(2) + ', ' + ang(r.xyz, 1, 2, 3).toFixed(2));
  near(dist(r.xyz, 1, 2), 1.2008, 0.01, 'C≡C');
});

test('cisplatin', () => {
  const r = FF.optimize(cisplatinDistorted());
  const x = r.xyz;
  ok(r.converged, 'converged');
  for (const [a, b] of [[1, 2], [3, 7], [1, 3], [2, 7]]) near(ang(x, a, 0, b), 90, 3, 'cis angle ' + a + '-Pt-' + b);
  for (const [a, b] of [[1, 7], [2, 3]]) near(ang(x, a, 0, b), 180, 3, 'trans angle ' + a + '-Pt-' + b);
  near(dist(x, 0, 3), 2.07, 0.02, 'Pt–N'); near(dist(x, 0, 7), 2.07, 0.02, 'Pt–N\'');
  near(dist(x, 0, 1), 2.38, 0.02, 'Pt–Cl'); near(dist(x, 0, 2), 2.38, 0.02, 'Pt–Cl\'');
  ok(planeDev(x, [0, 1, 2, 3, 7]) < 0.05, 'PtCl2N2 core planar');
  for (const N of [3, 7]) for (let h = 1; h <= 3; h++) near(dist(x, N, N + h), 1.01, 0.01, 'N–H');
});

test('metal centres: Zn tetrahedral, Ni planar/tetrahedral, Co octahedral, CN5', () => {
  // Zn(II): a tetrahedral complex stays tetrahedral
  const zn = complex(30, TET, [17, 17, 7, 7], [2.2, 2.3, 2.0, 1.9]);
  zn.xyz = jitter(zn.xyz, 0.12);
  let r = FF.optimize(zn);
  ok(r.converged, 'Zn converged');
  for (const a of ligandAngles(r.xyz, zn)) near(a, 109.47, 3, 'Zn L–M–L');
  near(dist(r.xyz, 0, 1), 1.22 + 1.02, 0.02, 'Zn–Cl');
  // Ni: near-planar input → square planar; tetrahedral input → tetrahedral
  const flat = complex(28, SQP.map(v => [v[0], v[1], v[0] * 0.15]).map(v => { const l = Math.hypot(...v); return v.map(c => c / l); }), [17, 17, 17, 17], [2.2, 2.2, 2.2, 2.2]);
  r = FF.optimize(flat);
  const sa = ligandAngles(r.xyz, flat).sort((a, b) => a - b);
  ok(sa[4] > 177 && sa[5] > 177 && Math.abs(sa[0] - 90) < 3, 'Ni near-planar → square planar: ' + sa.map(a => a.toFixed(0)).join(','));
  const niT = complex(28, TET, [17, 17, 17, 17], [2.2, 2.2, 2.2, 2.2]);
  r = FF.optimize(niT);
  for (const a of ligandAngles(r.xyz, niT)) near(a, 109.47, 2, 'Ni tetrahedral L–M–L');
  // Pt with a tetrahedral start still becomes square planar
  const ptT = complex(78, TET.map((v, k) => (k === 0 ? [v[0] * 1.3, v[1], v[2]] : v)).map(v => { const l = Math.hypot(...v); return v.map(c => c / l); }), [17, 17, 17, 17], [2.4, 2.4, 2.4, 2.4]);
  r = FF.optimize(ptT);
  const pa = ligandAngles(r.xyz, ptT).sort((a, b) => a - b);
  ok(pa[4] > 175 && pa[5] > 175 && Math.abs(pa[0] - 90) < 3, 'Pt tetrahedral start → square planar: ' + pa.map(a => a.toFixed(0)).join(','));
  // Co(NH3)6 octahedral
  const co = complex(27, OCT, [7, 7, 7, 7, 7, 7], [2.1, 1.9, 2.0, 2.05, 1.95, 2.0]);
  co.xyz = jitter(co.xyz, 0.1);
  r = FF.optimize(co);
  ok(r.converged, 'Co converged');
  const donors = co.Z.map((z, i) => (z === 7 ? i : -1)).filter(i => i >= 0);
  let cis = 0, trans = 0;
  for (let p = 0; p < 6; p++) for (let q = p + 1; q < 6; q++) {
    const a = ang(r.xyz, donors[p], 0, donors[q]);
    if (Math.abs(a - 90) < 4) cis++; else if (a > 176) trans++;
  }
  ok(cis === 12 && trans === 3, 'octahedral: ' + cis + ' cis, ' + trans + ' trans');
  for (const N of donors) near(dist(r.xyz, 0, N), 1.97, 0.02, 'Co–N');
  // CN5: trigonal-bipyramidal and square-pyramidal inputs keep their shape
  const tbp = [[0, 0, 1], [0, 0, -1], [1, 0, 0], [-0.5, Math.sqrt(3) / 2, 0], [-0.5, -Math.sqrt(3) / 2, 0]];
  r = FF.optimize({ ...complex(26, tbp, [6, 6, 6, 6, 6].map(() => 17), [2.3, 2.3, 2.3, 2.3, 2.3]), xyz: jitter(complex(26, tbp, [17, 17, 17, 17, 17], [2.3, 2.3, 2.3, 2.3, 2.3]).xyz, 0.1) });
  near(ang(r.xyz, 1, 0, 2), 180, 3, 'TBP axial–axial'); near(ang(r.xyz, 3, 0, 4), 120, 3, 'TBP equatorial');
  const spy = [[0, 0, 1], [1, 0, 0], [0, 1, 0], [-1, 0, 0], [0, -1, 0]];
  r = FF.optimize(complex(26, spy, [17, 17, 17, 17, 17], [2.3, 2.3, 2.3, 2.3, 2.3]));
  near(ang(r.xyz, 2, 0, 4), 180, 3, 'SP basal trans'); near(ang(r.xyz, 1, 0, 2), 90, 3, 'SP apex–basal');
  ok(FF.perceive(complex(26, spy, [17, 17, 17, 17, 17], [2.3, 2.3, 2.3, 2.3, 2.3])).atoms[0].geom === 'square-pyramidal', 'CN5 picks square pyramid');
});

test('hybridization perception', () => {
  const P = FF.perceive(mixedMolecule()).atoms;
  ok(P[1].hyb === 2 && P[1].theta0 === 120, 'carbonyl C sp2');
  ok(P[3].planarN && P[3].hyb === 2, 'amide N planar');
  ok(P[0].hyb === 3 && Math.abs(P[0].theta0 - 109.47) < 1e-9, 'methyl C sp3');
  const nh3 = FF.perceive({ Z: [7, 1, 1, 1], xyz: Float64Array.from([0, 0, 0, 1, 0, -0.3, -0.5, 0.87, -0.3, -0.5, -0.87, -0.3]), bonds: [[0, 1, 1], [0, 2, 1], [0, 3, 1]] }).atoms[0];
  ok(nh3.hyb === 3 && nh3.theta0 === 107, 'NH3 N sp3 107°');
  const co2 = FF.perceive({ Z: [6, 8, 8], xyz: Float64Array.from([0, 0, 0, 1.2, 0, 0, -1.2, 0, 0]), bonds: [[0, 1, 2], [0, 2, 2]] }).atoms[0];
  ok(co2.hyb === 1 && co2.theta0 === 180, 'CO2 carbon sp');
  const sf4 = FF.perceive({ Z: [16, 9, 9, 9, 9], xyz: Float64Array.from([0, 0, 0, 0, 0, 1.6, 0, 0, -1.6, 1.6, 0, 0, -0.8, 1.4, 0]), bonds: [[0, 1, 1], [0, 2, 1], [0, 3, 1], [0, 4, 1]] }).atoms[0];
  ok(sf4.geom === 'seesaw', 'SF4 seesaw (got ' + sf4.geom + ')');
  const so4 = FF.perceive({ Z: [16, 8, 8, 8, 8], xyz: Float64Array.from([0, 0, 0, 1, 1, 1, -1, -1, 1, -1, 1, -1, 1, -1, -1]), bonds: [[0, 1, 1], [0, 2, 1], [0, 3, 1], [0, 4, 1]] }).atoms[0];
  ok(so4.sn === 4 && Math.abs(so4.theta0 - 109.47) < 1e-9, 'sulfate drawn with single bonds stays tetrahedral');
});

test('frozen atoms do not move; dummies are ignored and fixed', () => {
  const m = cisplatinDistorted();
  const frozen = [0, 1, 4];
  const r = FF.optimize(m, { frozen });
  for (const i of frozen) for (let c = 0; c < 3; c++) ok(r.xyz[3 * i + c] === m.xyz[3 * i + c], 'atom ' + i + ' unchanged');
  ok(dist(r.xyz, 0, 3) !== dist(m.xyz, 0, 3), 'free atoms moved');
  // a dummy bonded to an atom changes nothing and never moves
  const w = water(), wd = { Z: [8, 1, 1, 0], xyz: Float64Array.from([...w.xyz, 0.3, 0.3, 0.3]), bonds: [...w.bonds, [0, 3, 1]] };
  const r1 = FF.optimize(w), r2 = FF.optimize(wd);
  for (let k = 0; k < 9; k++) near(r2.xyz[k], r1.xyz[k], 1e-9, 'dummy has no effect');
  ok(r2.xyz[9] === 0.3 && r2.xyz[10] === 0.3 && r2.xyz[11] === 0.3, 'dummy fixed');
  // isolated atoms only feel repulsion
  const iso = { Z: [11, 17], xyz: Float64Array.from([0, 0, 0, 1.0, 0, 0]), bonds: [] };
  const ri = FF.optimize(iso);
  ok(dist(ri.xyz, 0, 1) > 1.0 + 1e-3, 'overlapping ions pushed apart');
  const ev = FF.evaluate(FF.setup(iso), ri.xyz);
  ok(ev.terms.bond === 0 && ev.terms.angle === 0 && ev.terms.torsion === 0, 'no bonded terms for unbonded atoms');
});

test('gradient check (finite differences)', () => {
  const cases = [
    ['water', water()],
    ['benzene', { ...benzene(), xyz: jitter(benzene().xyz, 0.15) }],
    ['ethylene', ethylene(30)],
    ['cyclohexane', cyclohexaneFlat(0.5)],
    ['cisplatin', cisplatinDistorted()],
    ['Co(NH3)6', { ...complex(27, OCT, [7, 7, 7, 7, 7, 7], [2.1, 1.9, 2.0, 2.05, 1.95, 2.0]), xyz: jitter(complex(27, OCT, [7, 7, 7, 7, 7, 7], [2.1, 1.9, 2.0, 2.05, 1.95, 2.0]).xyz, 0.1) }],
    ['amide + contact', mixedMolecule()],
    ['acetylene', acetyleneBent()]
  ];
  // a crowded case: two waters squeezed together to exercise the repulsion term (incl. a 1-4 pair)
  const w2 = { Z: [8, 1, 1, 8, 1, 1], xyz: Float64Array.from([0, 0, 0, 0.96, 0, 0, -0.3, 0.9, 0, 1.6, 0.4, 0.3, 2.3, 0.9, 0.2, 1.4, -0.5, 0.6]), bonds: [[0, 1, 1], [0, 2, 1], [3, 4, 1], [3, 5, 1]] };
  cases.push(['crowded waters', w2]);
  const seen = { bond: 0, angle: 0, torsion: 0, oop: 0, vdw: 0 };
  for (const [name, mol] of cases) {
    const { err, terms } = gradError(mol);
    ok(err < 1e-5, name + ': max relative gradient error ' + err.toExponential(2));
    for (const k in seen) if (terms[k] > 1e-6) seen[k]++;
  }
  for (const k in seen) ok(seen[k] > 0, 'term "' + k + '" exercised by the gradient check');
  const fr = gradError(cisplatinDistorted(), { frozen: [1, 2] });
  ok(fr.err < 1e-5, 'gradient with frozen atoms');
});

test('robustness: degenerate starts, small rings, hypervalent centres', () => {
  // two H atoms exactly on the carbon and the other two collinear: still ends tetrahedral
  let r = FF.optimize({ Z: [6, 1, 1, 1, 1], xyz: Float64Array.from([0, 0, 0, 0, 0, 0, 1, 0, 0, -1, 0, 0, 0, 0, 0]), bonds: [[0, 1, 1], [0, 2, 1], [0, 3, 1], [0, 4, 1]] });
  ok(r.converged, 'coincident start converges');
  for (const [a, b] of [[1, 2], [1, 3], [1, 4], [2, 3], [2, 4], [3, 4]]) near(ang(r.xyz, a, 0, b), 109.47, 1, 'H–C–H from a coincident start');
  r = FF.optimize({ Z: [8, 1, 1], xyz: Float64Array.from([0, 0, 0, 0.96, 0, 0, -0.96, 0, 0]), bonds: [[0, 1, 1], [0, 2, 1]] });
  near(ang(r.xyz, 1, 0, 2), 104.5, 0.5, 'exactly linear water start bends');
  r = FF.optimize({ Z: [8, 1, 1], xyz: Float64Array.from([0, 0, NaN, 0.96, 0, 0, -0.5, 0.8, 0]), bonds: [[0, 1, 1], [0, 2, 1]] });
  ok(r.converged === false, 'non-finite input reports not converged (no throw)');
  r = FF.optimize({ Z: [], xyz: new Float64Array(0), bonds: [] });
  ok(r.converged && r.steps === 0, 'empty molecule');
  // cyclopropane (no ring-angle term) and puckered cyclobutane (90° ring angles)
  const ringMol = (n, h) => {
    const Z = [], xyz = [], bonds = [], R = 1.52 / (2 * Math.sin(Math.PI / n));
    for (let k = 0; k < n; k++) { const a = 2 * Math.PI * k / n; Z.push(6); xyz.push(R * Math.cos(a), R * Math.sin(a), k % 2 ? h : -h); bonds.push([k, (k + 1) % n, 1]); }
    for (let k = 0; k < n; k++) {
      const a = 2 * Math.PI * k / n, z0 = k % 2 ? h : -h;
      for (const s of [1, -1]) { Z.push(1); xyz.push((R + 0.6) * Math.cos(a), (R + 0.6) * Math.sin(a), z0 + s * 0.9); bonds.push([k, Z.length - 1, 1]); }
    }
    return { Z, xyz: Float64Array.from(xyz), bonds };
  };
  r = FF.optimize(ringMol(3, 0));
  near(dist(r.xyz, 0, 1), 1.52, 0.02, 'cyclopropane C–C'); near(ang(r.xyz, 3, 0, 4), 117, 5, 'cyclopropane H–C–H opens');
  r = FF.optimize(ringMol(4, 0.1));
  near(ang(r.xyz, 0, 1, 2), 88, 3, 'cyclobutane C–C–C');
  const pk = Math.abs(dih(r.xyz, 0, 1, 2, 3));
  ok(pk > 15 && pk < 35, 'cyclobutane puckered (' + pk.toFixed(1) + '°)');
  // hypervalent main-group centres follow VSEPR polyhedra
  const star = (Zc, g) => {
    const Z = [Zc], xyz = [0, 0, 0], bonds = [];
    FF.GEOMETRIES[g].vectors.forEach((v, k) => { Z.push(9); xyz.push(1.6 * v[0] + 0.12 * Math.sin(7 * k), 1.6 * v[1] + 0.12 * Math.cos(3 * k), 1.6 * v[2] + 0.1 * Math.sin(k)); bonds.push([0, k + 1, 1]); });
    return { Z, xyz: Float64Array.from(xyz), bonds };
  };
  for (const [name, Zc, g] of [['SF6', 16, 'octahedral'], ['PF5', 15, 'trigonal-bipyramidal'], ['SF4', 16, 'seesaw'], ['ClF3', 17, 't-shaped'], ['XeF4', 54, 'square-planar']]) {
    const m = star(Zc, g), o = FF.optimize(m);
    ok(FF.perceive(m).atoms[0].geom === g, name + ' perceived as ' + g);
    const got = [], want = [], V = FF.GEOMETRIES[g].vectors;
    for (let p = 1; p < m.Z.length; p++) for (let q = p + 1; q < m.Z.length; q++) {
      got.push(ang(o.xyz, p, 0, q));
      want.push(Math.acos(Math.max(-1, Math.min(1, V[p - 1][0] * V[q - 1][0] + V[p - 1][1] * V[q - 1][1] + V[p - 1][2] * V[q - 1][2]))) * 180 / Math.PI);
    }
    got.sort((a, b) => a - b); want.sort((a, b) => a - b);
    ok(got.every((a, k) => Math.abs(a - want[k]) < 3), name + ' angles ideal');
  }
});

test('createOptimizer (live)', () => {
  const m = cisplatinDistorted();
  const opt = FF.createOptimizer(m, { frozen: [] });
  let s = opt.step(5);
  const e5 = s.energy;
  ok(s.xyz.length === 33 && Number.isFinite(e5), 'step returns state');
  s = opt.step(10);
  ok(s.energy <= e5 + 1e-9, 'energy decreases');
  // the user drags atom 1 somewhere else while it is held
  const x = Float64Array.from(s.xyz); x[3] -= 0.8; x[4] += 0.5;
  opt.setFrozen([1]); opt.setPositions(x);
  for (let k = 0; k < 40; k++) s = opt.step(5);
  ok(s.xyz[3] === x[3] && s.xyz[4] === x[4], 'held atom stays where it was dragged');
  ok(s.converged || s.rmsGrad < 0.05, 'converges while atom is held (rms ' + s.rmsGrad.toFixed(3) + ')');
  opt.setFrozen([]);
  for (let k = 0; k < 60 && !s.converged; k++) s = opt.step(10);
  ok(s.converged, 'converges after release');
  near(ang(s.xyz, 1, 0, 2), 90, 3, 'cis angle restored');
  // rebuild with a new topology (an extra isolated water)
  const w = water();
  const bigger = { Z: [...m.Z, ...w.Z], xyz: Float64Array.from([...s.xyz, ...w.xyz.map((v, k) => v + (k % 3 === 0 ? 6 : 0))]), bonds: [...m.bonds, ...w.bonds.map(b => [b[0] + 11, b[1] + 11, b[2]])] };
  opt.rebuild(bigger);
  for (let k = 0; k < 50 && !s.converged; k++) s = opt.step(10);
  s = opt.step(200);
  ok(s.xyz.length === 42, 'rebuild changes the atom count');
  near(ang(s.xyz, 12, 11, 13), 104.5, 1, 'new water optimized after rebuild');
  const cheap = Date.now();
  for (let k = 0; k < 100; k++) opt.step(1);
  ok(Date.now() - cheap < 500, 'step calls are cheap');
  let threw = false;
  try { opt.setPositions(new Float64Array(5)); } catch (e) { threw = true; }
  ok(threw, 'setPositions rejects wrong length');
});

test('performance: 1000-atom water box, 50 steps', () => {
  const Z = [], xyz = [], bonds = [];
  let k = 0;
  for (let a = 0; a < 7 && k < 333; a++) for (let b = 0; b < 7 && k < 333; b++) for (let c = 0; c < 7 && k < 333; c++, k++) {
    const o = Z.length, cx = a * 3.1, cy = b * 3.1, cz = c * 3.1, th = rnd() * 6.28, ph = rnd() * 3.14;
    Z.push(8, 1, 1);
    xyz.push(cx, cy, cz, cx + 1.05 * Math.cos(th), cy + 1.05 * Math.sin(th), cz + 0.2 * rnd(),
      cx + 0.9 * Math.cos(th + 1.7), cy + 0.9 * Math.sin(th + 1.7), cz + 0.4 * Math.cos(ph));
    bonds.push([o, o + 1, 1], [o, o + 2, 1]);
  }
  Z.push(11); xyz.push(-2.5, -2.5, -2.5);
  ok(Z.length === 1000, '1000 atoms');
  const t0 = Date.now();
  const r = FF.optimize({ Z, xyz: Float64Array.from(xyz), bonds }, { maxSteps: 50, tol: 1e-9 });
  const ms = Date.now() - t0;
  ok(r.steps === 50, '50 steps taken (' + r.steps + ')');
  ok(ms < 3000, '50 steps in ' + ms + ' ms (< 3000)');
  let worst = 0;
  for (let w = 0; w < 333; w++) worst = Math.max(worst, Math.abs(ang(r.xyz, 3 * w + 1, 3 * w, 3 * w + 2) - 104.5));
  ok(worst < 5, 'waters relaxed (worst H–O–H error ' + worst.toFixed(2) + '°)');
});

// ---------------------------------------------------------------------------------------------- summary
console.log('\nBlueMolecule force-field tests: ' + pass + ' passed, ' + fail + ' failed');
if (fail) {
  for (const f of failures) console.log('  FAIL ' + f);
  process.exit(1);
}
