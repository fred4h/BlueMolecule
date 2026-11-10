#!/usr/bin/env node
/* BlueMolecule tests for src/analysis.js (superposition, RMSD, convergence series, scan/IRC profiles, NMR sticks,
 * orbital levels) and the ECD simulation in src/spectrum.js.
 * Run:  node tests/test_analysis.js */
'use strict';
const fs = require('fs');
const path = require('path');
const ROOT = path.join(__dirname, '..');
for (const f of ['elements.js', 'parsers.js', 'readers-struct.js', 'writers.js', 'readers-qc.js', 'spectrum.js', 'analysis.js']) require(path.join(ROOT, 'src', f));
const BM = globalThis.BM, A = BM.analysis;
const FIX = path.join(__dirname, 'fixtures');
const load = f => BM.parseFile(path.basename(f), fs.readFileSync(path.join(FIX, f), 'utf8'));

let pass = 0, fail = 0, current = '';
const failures = [];
function ok(c, msg) { if (c) pass++; else { fail++; failures.push('[' + current + '] ' + msg); } }
function eq(a, b, msg) { ok(a === b, msg + ': got ' + JSON.stringify(a) + ', expected ' + JSON.stringify(b)); }
function near(a, b, tol, msg) { ok(typeof a === 'number' && Math.abs(a - b) <= tol, msg + ': got ' + a + ', expected ' + b); }
function test(name, fn) { current = name; try { fn(); } catch (e) { fail++; failures.push('[' + name + '] threw: ' + (e && e.stack || e)); } }

// deterministic pseudo-random numbers
let seed = 12345;
const rnd = () => { seed = (seed * 1103515245 + 12345) % 2147483648; return seed / 2147483648; };
function rotMat(axis, th) {
  const [x, y, z] = axis.map(v => v / Math.hypot(...axis)), c = Math.cos(th), s = Math.sin(th), C = 1 - c;
  return [[c + x * x * C, x * y * C - z * s, x * z * C + y * s], [y * x * C + z * s, c + y * y * C, y * z * C - x * s], [z * x * C - y * s, z * y * C + x * s, c + z * z * C]];
}
function transform(X, R, t) {
  const out = new Float64Array(X.length);
  for (let i = 0; i < X.length / 3; i++) for (let p = 0; p < 3; p++) out[3 * i + p] = R[p][0] * X[3 * i] + R[p][1] * X[3 * i + 1] + R[p][2] * X[3 * i + 2] + t[p];
  return out;
}

test('Kabsch superposition and RMSD', () => {
  const n = 12, X = new Float64Array(3 * n);
  for (let i = 0; i < 3 * n; i++) X[i] = (rnd() - 0.5) * 6;
  for (const [axis, th] of [[[1, 2, 3], 0.7], [[0, 0, 1], Math.PI], [[1, -1, 0.2], 2.9], [[0.3, 0.1, -1], 1e-3]]) {
    const Y = transform(X, rotMat(axis, th), [1.5, -2, 0.25]);
    near(A.rmsd(X, Y), 0, 1e-6, 'rotated + translated copy: RMSD 0 (θ = ' + th + ')');
    const S = A.superpose(X, Y);
    let m = 0; for (let i = 0; i < 3 * n; i++) m = Math.max(m, Math.abs(S[i] - X[i]));
    near(m, 0, 1e-6, 'superposed copy matches the reference');
    near(A.maxDeviation(X, Y), 0, 1e-6, 'largest deviation 0 after alignment');
    ok(A.rmsd(X, Y, { align: false }) > 0.5, 'unaligned RMSD is large');
  }
  // perturbed copy: the aligned RMSD equals the RMSD of the perturbation when no rigid motion is involved
  const P = Float64Array.from(X, v => v + (rnd() - 0.5) * 0.02);
  const r0 = A.rmsd(X, P, { align: false }), r1 = A.rmsd(X, P);
  ok(r1 <= r0 + 1e-12 && r1 > 0.25 * r0, 'aligned RMSD ≤ unaligned for a small perturbation: ' + r1 + ' vs ' + r0);
  near(A.rmsd(X, transform(P, rotMat([2, 1, 0], 1.1), [3, 3, 3])), r1, 1e-9, 'RMSD invariant to rigid motion');
  // two atoms: distance change only
  const a = Float64Array.from([0, 0, 0, 1, 0, 0]), b = Float64Array.from([5, 5, 5, 5, 6.2, 5]);
  near(A.rmsd(a, b), 0.1, 1e-9, 'diatomic RMSD = |Δr|/2');
  near(A.maxDeviation(a, b), 0.1, 1e-9, 'diatomic largest deviation');
});

test('convergence series: Gaussian optimisation', () => {
  const ds = load('g16_water_optfreq.log');
  const c = A.convergenceSeries(ds, {});
  eq(c.source, 'frames', 'source');
  eq(c.steps.length, 4, 'four steps');
  ok(c.isOpt, 'is an optimisation');
  const get = k => c.series.find(s => s.key === k);
  for (const k of ['energy', 'dE', 'maxForce', 'rmsForce', 'maxDisp', 'rmsDisp', 'energyChange', 'rmsd', 'maxStep']) ok(get(k), 'series ' + k);
  near(get('maxForce').threshold, 0.00045, 1e-12, 'max force threshold');
  near(get('rmsDisp').threshold, 0.0012, 1e-12, 'RMS displacement threshold');
  ok(get('maxForce').log && !get('energy').log, 'log flags');
  near(get('maxForce').values[3], 0.000027, 1e-12, 'last max force');
  eq(get('rmsd').values[3], 0, 'RMSD to final ends at 0');
  ok(get('rmsd').values[0] > 0.001, 'first RMSD > 0');
  eq(get('maxStep').values[0], null, 'no step before the first geometry');
  near(get('dE').values[3], 0, 1e-9, 'ΔE vs final ends at 0');
  ok(get('dE').values[0] > 0.3, 'first ΔE positive (kcal/mol)');
  eq(c.ratios.length, 4, 'four criteria with thresholds');
  ok(c.ratios.every(r => r.values[3] < 1), 'converged below 1');
  eq(c.converged, 3, 'converged at the last step');
  const j1 = A.convergenceSeries(ds, { job: 1 });
  eq(j1.steps.length, 1, 'freq job: one geometry');
  eq(j1.series.length, 0, 'no series for one geometry');
  eq(A.jobRange(ds, 1).join(','), '3,3', 'job range');
});

test('convergence series: ORCA, xtb, MOPAC, CP2K', () => {
  const o = A.convergenceSeries(load('orca5_water_optfreq.out'), {});
  const eg = o.series.find(s => s.key === 'energyChange');
  near(eg.threshold, 5e-6, 1e-15, 'ORCA energy change threshold');
  ok(o.series.find(s => s.key === 'rmsd').values[2] === 0, 'ORCA RMSD ends at 0');
  const x = A.convergenceSeries(load('qc/xtb_water_ohess.out'), {});
  eq(x.source, 'optSteps', 'xtb cycles without geometries');
  ok(x.steps.length === 3 && x.series.some(s => s.key === 'gradNorm' && s.threshold === 0.001), 'xtb gradient norm with threshold');
  ok(!x.series.some(s => s.key === 'rmsd'), 'no RMSD without geometries');
  eq(x.steps[2].frame, 0, 'last cycle linked to the final structure');
  const m = A.convergenceSeries(load('qc/mopac_water_opt_force.out'), {});
  ok(m.source === 'optSteps' && m.series.find(s => s.key === 'gradNorm').unit === 'kcal/mol/Å', 'MOPAC GNORM unit');
  const p = A.convergenceSeries(load('qc/cp2k_water_geoopt.out'), {});
  ok(p.source === 'optSteps' && p.series.filter(s => s.threshold).length === 4, 'CP2K criteria with thresholds');
  const t = A.convergenceSeries(load('water_traj.xyz'), {});
  ok(!t.isOpt, 'an XYZ trajectory is not an optimisation');
});

test('scan and IRC profiles', () => {
  const s = A.profiles(load('qc/g16_water_modred_scan.log'), {});
  eq(s.length, 1, 'one profile');
  eq(s[0].kind, 'scan', 'scan');
  eq(s[0].points.map(p => p.x).join(','), '100,105', 'scan values');
  eq(s[0].xLabel, 'A(2,1,3) (°)', 'scan axis label');
  const g = A.profiles(load('g16_water_scan.log'), {});
  ok(g.length === 1 && g[0].kind === 'scan' && g[0].points.length === 2, 'Gaussian relaxed scan from converged points');
  const i = A.profiles(load('qc/g16_hcn_irc.log'), {});
  ok(i.length === 1 && i[0].kind === 'irc' && i[0].points.length === 5, 'Gaussian IRC points');
  eq(i[0].points.map(p => p.x).join(','), '-0.2,-0.1,0,0.1,0.2', 'IRC sorted by reaction coordinate');
  const oi = A.profiles(load('qc/orca5_hcn_irc.out'), {});
  ok(oi.length === 1 && oi[0].xLabel === 'IRC step' && oi[0].points.length >= 3 && oi[0].points[0].frame === null, 'ORCA IRC path without geometries');
});

test('NMR sticks and orbital levels', () => {
  const ds = load('qc/g16_water_nmr.log'), Z = ds.frames[0].Z;
  const nuc = A.nmrNuclei(ds, Z);
  eq(nuc.map(n => n.label + ':' + n.count).join(' '), '1H:2 17O:1', 'nuclei');
  let st = A.nmrSticks(ds, Z, 1, {});
  ok(st.length === 1 && st[0].count === 2, 'equivalent protons merged');
  near(st[0].x, 31.2341, 1e-9, 'shielding plotted when ref = 0');
  st = A.nmrSticks(ds, Z, 1, { ref: 31.88 });
  near(st[0].x, 31.88 - 31.2341, 1e-9, 'δ = σref − σ');
  eq(A.nmrSticks(ds, Z, 1, { merge: false }).length, 2, 'unmerged sticks');
  const lv = A.moLevels(load('review/g16_oh_radical_pop.log'));
  eq(lv.map(c => c.label).join(','), 'α,β', 'open-shell columns');
  eq(lv[0].homo, 4, 'alpha HOMO index');
  eq(lv[1].homo, 3, 'beta HOMO index');
  ok(lv[0].levels[4].occ && !lv[0].levels[5].occ, 'occupations');
  eq(A.moLevels(load('g16_water_optfreq.log')).length, 1, 'closed shell: one column');
});

test('ECD spectrum', () => {
  const S = load('qc/orca5_water_cd.out').excitedStates;
  const r = BM.spectrum.ecd(S, { xmin: 2, xmax: 12, unit: 'eV', sigmaEV: 0.3, n: 2001 });
  eq(r.sticks.length, S.filter(s => isFinite(s.R)).length, 'one stick per state with R');
  const s0 = S[0], i0 = Math.round((s0.energyEV - 2) / 10 * 2000);
  ok(Math.sign(r.curve.y[i0]) === Math.sign(s0.R) || Math.abs(r.curve.y[i0]) < 1e-9, 'sign of Δε follows R near the first state');
  // single state: peak value = E R / (22.96 √π σ)
  const one = BM.spectrum.ecd([{ energyEV: 5, R: 10, index: 1 }], { xmin: 5, xmax: 5.0001, unit: 'eV', sigmaEV: 0.25, n: 2 });
  near(one.curve.y[0], 5 * 10 / (22.96 * Math.sqrt(Math.PI) * 0.25), 1e-9, 'peak height');
  const nm = BM.spectrum.ecd([{ energyEV: 5, R: -10, index: 1 }], { xmin: 200, xmax: 300, unit: 'nm', sigmaEV: 0.25, n: 101 });
  ok(nm.curve.y.every(v => v <= 0) && nm.sticks[0].y === -10, 'negative R gives a negative band');
});

console.log('\nBlueMolecule analysis tests: ' + pass + ' passed, ' + fail + ' failed');
if (fail) { for (const f of failures) console.log('  FAIL ' + f); process.exit(1); }
