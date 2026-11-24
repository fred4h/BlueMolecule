#!/usr/bin/env node
/* BlueMolecule tests for src/vibtools.js (issue #14): the character of a stationary point from its frequencies, the
 * two ends of a mode, the bonds that form and break along it, the meshes of both ends, and the routes / keywords of
 * the forward and reverse optimisations; with the negative Hessian eigenvalues of transition-state searches
 * (Gaussian and Q-Chem readers, BM.analysis.convergenceSeries).
 * Run:  node tests/test_vibtools.js */
'use strict';
const fs = require('fs');
const path = require('path');
const ROOT = path.join(__dirname, '..');
for (const f of ['elements.js', 'bigfile.js', 'parsers.js', 'condense.js', 'readers-struct.js', 'writers.js', 'readers-qc.js', 'geom.js', 'chem.js', 'analysis.js', 'cartoon.js', 'vibtools.js'])
  require(path.join(ROOT, 'src', f));
const BM = globalThis.BM, V = BM.vibtools;
const FIX = path.join(__dirname, 'fixtures');
const read = f => fs.readFileSync(path.join(FIX, f), 'utf8');
const load = f => BM.parseFile(path.basename(f), read(f));
const expected = JSON.parse(read('expected.json'));

let pass = 0, fail = 0, current = '';
const failures = [];
function ok(c, msg) { if (c) pass++; else { fail++; failures.push('[' + current + '] ' + msg); } }
function eq(a, b, msg) { ok(a === b, msg + ': got ' + JSON.stringify(a) + ', expected ' + JSON.stringify(b)); }
function near(a, b, tol, msg) { ok(typeof a === 'number' && Math.abs(a - b) <= tol, msg + ': got ' + a + ', expected ' + b); }
function test(name, fn) { current = name; try { fn(); } catch (e) { fail++; failures.push('[' + name + '] threw: ' + (e && e.stack || e)); } }

test('stationary points', () => {
  const st = f => V.stationary({ freqs: f });
  eq(st([120, 340, 1600]).kind, 'minimum', 'no imaginary frequency');
  eq(st([-0.4, 0, 0.3, 120]).kind, 'minimum', 'translations and rotations printed as -0.4 are not imaginary');
  const ts = st([210, -1532.4, 980]);
  ok(ts.kind === 'ts' && ts.imag.join() === '1' && ts.large.join() === '1' && !ts.small.length, 'one imaginary frequency: a transition state, mode 2');
  const tsn = st([-12.3, 300, -1532.4]);
  ok(tsn.kind === 'ts' && tsn.imag.join() === '2,0' && tsn.small.join() === '0', 'a TS with a small imaginary mode: most negative first: ' + tsn.imag);
  eq(st([-1532.4, -420.1, -15.2, 200]).kind, 'saddle', 'two large imaginary frequencies');
  eq(st([-15.2, 200]).kind, 'noise', 'a small imaginary frequency only');
  eq(V.stationary({ freqs: [-80, 200] }, { small: 100 }).kind, 'noise', 'the threshold can be changed');
  eq(st([]).kind, 'minimum', 'no frequencies');
});

test('descriptions', () => {
  const d = (f, ts) => V.describe({ freqs: f }, null, ts);
  eq(d([100, 200], false).text, 'No imaginary frequencies: a minimum.', 'minimum');
  ok(d([100, 200], false).ok && !d([100, 200], true).ok, 'a minimum is fine, except for a TS job');
  eq(d([100, 200], true).text, 'No imaginary frequencies: this is a minimum, not a transition state.', 'TS job without imaginary frequency');
  eq(d([-1532.4, 200], true).text, 'One imaginary frequency (1532.4i cm⁻¹): a transition state.', 'transition state');
  eq(d([-1532.4, 200], true).short, 'one imaginary frequency (1532.4i cm⁻¹)', 'short form');
  ok(/^2 imaginary frequencies \(1532\.4i, 12\.3i cm⁻¹\): one transition-state mode and one small one below 50i cm⁻¹, often numerical noise/.test(d([-1532.4, -12.3, 200]).text), 'TS with noise: ' + d([-1532.4, -12.3, 200]).text);
  ok(/higher-order saddle point/.test(d([-1532.4, -420.1, 200]).text) && !d([-1532.4, -420.1, 200]).ok, 'saddle point');
  ok(/^One imaginary frequency \(15\.2i cm⁻¹\), below 50i cm⁻¹: often numerical noise in a minimum .*unless this is a very flat transition state\.$/.test(d([-15.2, 200]).text), 'noise: ' + d([-15.2, 200]).text);
  ok(/a transition state needs one clearly imaginary frequency\.$/.test(d([-15.2, 200], true).text), 'noise in a TS job');
  ok(/^6 imaginary frequencies \(90\.0i, 80\.0i, 70\.0i, 60\.0i, … cm⁻¹\)/.test(d([-90, -80, -70, -60, -55, -51]).text), 'long lists are cut: ' + d([-90, -80, -70, -60, -55, -51]).text);
});

test('ends of a mode', () => {
  const eq0 = Float64Array.from([0, 0, 0, 0, 0, 1, 1, 0, 0]), mode = [0, 0, 2, 0, 0, -1, 0, 0, 0];
  const e = V.ends(eq0, mode, 0.3);
  near(e.u[2], 0.3, 1e-12, 'the largest shift is the amplitude');
  near(e.u[5], -0.15, 1e-12, 'the others in proportion');
  ok(e.minus[2] === -0.3 && e.plus[2] === 0.3 && e.minus[5] === 1.15 && e.plus[5] === 0.85, 'eq - u and eq + u');
});

const sn2 = load('g16_sn2_ts_optfreq.log'), sv = sn2.vibrations, seq = sn2.frames[sv.frameIndex];
test('the transition-state fixture', () => {
  const ex = expected.g16_sn2_ts;
  eq(sv.freqs.length, 12, 'twelve vibrations');
  ok(sv.freqs.every((f, i) => Math.abs(f - ex.freqs[i]) < 1e-9), 'frequencies as written: ' + sv.freqs.slice(0, 3));
  eq(JSON.stringify(sv.sym), JSON.stringify(ex.sym), 'D3h labels, A2" for the imaginary mode');
  ok(sn2.meta.jobTypes.indexOf('ts') >= 0, 'the route makes it a TS job');
  eq(sn2.warnings.length, 0, 'no reader notes (the thermochemistry box right after the last block ends the section)');
  eq(V.stationary(sv).kind, 'ts', 'a transition state');
});

test('bonds that form and break', () => {
  const u = V.ends(seq.xyz, sv.modes[0], 0.35).u;
  const rows = V.bondChanges(seq.Z, seq.xyz, u, { imaginary: true });
  eq(rows.length, 2, 'the two C–F bonds and nothing else');
  const [a, b] = rows;
  ok(a.i === 0 && a.j === 1 && a.change === 'forms' && b.i === 0 && b.j === 2 && b.change === 'breaks', 'C1–F2 forms and C1–F3 breaks towards +: ' + JSON.stringify(rows.map(r => [r.i, r.j, r.change])));
  near(a.d0, 1.85, 1e-4, 'C–F at the TS');
  near(a.dMinus - a.d0, a.d0 - a.dPlus, 1e-9, 'symmetric distances at the two ends');
  near(a.delta, -b.delta, 1e-9, 'opposite changes');
  ok(!a.bonded0 && a.bondedPlus && !a.bondedMinus, 'the C–F bond is drawn at the + end only');
  near(a.ratio, a.d0 / (BM.el.covalent[6] + BM.el.covalent[9]), 1e-12, 'ratio to the covalent radii');
  // a real mode: no "forms" or "breaks"
  const sym = sv.freqs.findIndex(f => f > 3000);
  const r2 = V.bondChanges(seq.Z, seq.xyz, V.ends(seq.xyz, sv.modes[sym], 0.35).u, { imaginary: false });
  ok(r2.length === 3 && r2.every(r => r.i === 0 && seq.Z[r.j] === 1 && (r.change === 'shorter' || r.change === 'longer')), 'the symmetric C–H stretch: three C–H, shorter or longer: ' + JSON.stringify(r2.map(r => [r.i, r.j, r.change])));
  ok(!V.bondChanges(seq.Z, seq.xyz, V.ends(seq.xyz, sv.modes[sym], 0.35).u, { imaginary: true }).some(r => r.change === 'forms' || r.change === 'breaks'), 'C–H bonds of normal length do not form or break');
  // an angle (H–O–H) is no bond: the water bend changes only angles... and O–H a little
  const w = load('g16_water_optfreq.log'), wf = w.frames[w.vibrations.frameIndex];
  const bend = V.bondChanges(wf.Z, wf.xyz, V.ends(wf.xyz, w.vibrations.modes[0], 0.35).u, { fraction: 0 });
  ok(bend.length && bend.every(r => r.i === 0), 'the H···H pair of water spans an angle and is left out: ' + JSON.stringify(bend.map(r => [r.i, r.j])));
  // atoms that do not move, far away, add nothing; no displacement, no rows
  const Z = Array.from(seq.Z).concat([8]), X = Float64Array.from(Array.from(seq.xyz).concat([20, 20, 20])), U = Float64Array.from(Array.from(u).concat([0, 0, 0]));
  eq(V.bondChanges(Z, X, U, { imaginary: true }).length, 2, 'a far, still atom changes nothing');
  eq(V.bondChanges(seq.Z, seq.xyz, new Float64Array(18)).length, 0, 'a zero displacement gives no rows');
  eq(V.bondChanges(seq.Z, seq.xyz, new Float64Array(3)).length, 0, 'a displacement of the wrong length gives no rows');
  eq(V.bondChanges(seq.Z, seq.xyz, u, { max: 1 }).length, 1, 'max rows');
});

test('meshes of both ends', () => {
  const u = V.ends(seq.xyz, sv.modes[0], 0.35).u;
  const m = V.endMeshes(seq.Z, seq.xyz, u);
  ok(m.length === 2 && m[0].sign === -1 && m[1].sign === 1 && m.every(x => x.opacity === 0.5), 'two translucent meshes, - then +');
  ok(m.every(x => x.positions.length === x.normals.length && x.positions.length === x.colors.length && x.indices.length % 3 === 0), 'complete meshes');
  ok(m[0].colors[2] > m[0].colors[0] && m[1].colors[0] > m[1].colors[2], '- is blue, + is orange');
  let maxI = 0; for (const i of m[1].indices) maxI = Math.max(maxI, i);
  ok(maxI < m[1].positions.length / 3, 'indices within the vertices');
  // the + end: 6 spheres (42 vertices) and 4 bonds (3 C–H and the C–F that formed; 18 vertices each)
  eq(m[1].positions.length / 3, 6 * 42 + 4 * 18, 'vertices of the + end');
  const noH = V.endMeshes(seq.Z, seq.xyz, u, { radius: i => (seq.Z[i] === 1 ? 0 : 0.3) });
  eq(noH[1].positions.length / 3, 3 * 42 + 18, 'radius 0 hides an atom and its bonds (hydrogens hidden)');
  eq(V.endMeshes(seq.Z, seq.xyz, new Float64Array(18)).length, 0, 'no displacement, no meshes');
  // a large structure: only the atoms that move, at most maxAtoms
  const n = 2000, Z = new Array(n).fill(6), X = new Float64Array(3 * n), U = new Float64Array(3 * n);
  for (let i = 0; i < n; i++) { X[3 * i] = (i % 20) * 1.5; X[3 * i + 1] = Math.floor(i / 20) % 10 * 1.5; X[3 * i + 2] = Math.floor(i / 200) * 1.5; }
  for (const i of [0, 1, 2, 3, 4]) U[3 * i + 2] = 0.1 * (i + 1);
  const big = V.endMeshes(Z, X, U);
  eq(big[0].positions.length / 3, 5 * 42 + 4 * 18, 'a large structure: the five moving atoms (and their four bonds) only');
  eq(V.endMeshes(Z, X, U, { maxAtoms: 2 })[0].positions.length / 3, 2 * 42 + 18, 'at most maxAtoms, the largest shifts first');
});

test('routes and keywords of the forward and reverse optimisations', () => {
  const r = V.minimumRoute;
  eq(r('#p opt=(ts,calcfc,noeigentest) freq b3lyp/6-31g(d)'), '#p opt freq b3lyp/6-31g(d)', 'TS options removed');
  eq(r('# Opt=(TS, CalcFC, MaxCycles=100) Freq B3LYP/6-31G(d) SCRF=(SMD,Solvent=Water)'), '# Opt=(MaxCycles=100) Freq B3LYP/6-31G(d) SCRF=(SMD,Solvent=Water)', 'other options and keywords kept');
  eq(r('#p opt=(ts,calcfc,tight) b3lyp/6-31g(d)'), '#p opt=tight b3lyp/6-31g(d)', 'one option left: opt=tight');
  eq(r('#p freq b3lyp/6-31g(d) geom=check guess=read'), '#p opt freq b3lyp/6-31g(d)', 'opt added; geom= and guess= dropped');
  eq(r('#P B3LYP/6-31G* Opt=QST3'), '#P B3LYP/6-31G* Opt', 'QST3');
  eq(r('#opt=ts b3lyp/sto-3g'), '# opt b3lyp/sto-3g', 'a # written with the first keyword');
  eq(r('#p opt(ts,modredundant) b3lyp/6-31g(d)'), '#p opt b3lyp/6-31g(d)', 'modredundant needs input that is not written');
  eq(r('#p IRC=(calcfc,maxpoints=20) b3lyp/6-31g(d)'), '#p opt b3lyp/6-31g(d)', 'an IRC route becomes an optimisation');
  eq(r('#n M062X/def2TZVP opt=(calcfc,ts,saddle=1) int=ultrafine'), '#n M062X/def2TZVP opt int=ultrafine', 'saddle=1');
  eq(V.routeTokens('#p opt = ( ts , calcfc ) b3lyp/6-31g(d,p)').join('|'), '#p|opt=(ts,calcfc)|b3lyp/6-31g(d,p)', 'route tokens: spaces inside parentheses');
  const k = V.minimumKeywords;
  eq(k('! B3LYP D3BJ def2-SVP OptTS Freq'), 'B3LYP D3BJ def2-SVP Opt Freq', 'OptTS becomes Opt');
  eq(k(['B3LYP', 'def2-SVP', 'NumFreq']), 'B3LYP def2-SVP NumFreq Opt', 'Opt added (keywords as an array)');
  eq(k('! r2SCAN-3c ScanTS'), 'r2SCAN-3c Opt', 'ScanTS');
  eq(k('! B3LYP def2-SVP TightOpt Freq'), 'B3LYP def2-SVP TightOpt Freq', 'an optimisation already: kept');
  eq(k('! B3LYP def2-SVP NEB-TS Freq'), 'B3LYP def2-SVP Opt Freq', 'NEB-TS');
});

test('negative Hessian eigenvalues', () => {
  eq(sn2.frames.map(f => f.opt && f.opt.negEig).join(), expected.g16_sn2_ts.negEig.join(), 'Gaussian: one negative eigenvalue at each step of the TS search (one step prints a line before them)');
  const c = BM.analysis.convergenceSeries(sn2, {}), s = c.series.find(x => x.key === 'negEig');
  ok(s && s.integer && !s.log && s.values.join() === '1,1,1,1', 'a plot series of whole numbers: ' + JSON.stringify(s && s.values));
  const w = load('g16_water_optfreq.log');
  ok(!w.frames.some(f => f.opt && f.opt.negEig !== undefined) && !BM.analysis.convergenceSeries(w, {}).series.some(x => x.key === 'negEig'), 'no Hessian printed, no series');
  // a minimisation whose Hessian has only positive eigenvalues: counted, but no plot (all 0)
  const pos = read('g16_water_optfreq.log').replace(/( Item {15}Value {5}Threshold {2}Converged\?)/g, ' ITU=  0\n     Eigenvalues ---    0.16000   0.47012   0.56000\n$1');
  const dp = BM.parseFile('w.log', pos);
  ok(dp.frames.every(f => !f.opt || f.opt.negEig === 0) && !BM.analysis.convergenceSeries(dp, {}).series.some(x => x.key === 'negEig'), 'all zero in a minimisation: no series');
  // ... and one with a negative eigenvalue: the series appears (a minimisation that is not on a minimum)
  const neg = read('g16_water_optfreq.log').replace(/( Item {15}Value {5}Threshold {2}Converged\?)/g, ' ITU=  0\n Use linear search instead of GDIIS.\n     Eigenvalues ---   -0.01000   0.47012   0.56000\n$1');
  const dn = BM.parseFile('w.log', neg), sn = BM.analysis.convergenceSeries(dn, {}).series.find(x => x.key === 'negEig');
  ok(sn && sn.values.every(v => v === 1), 'a negative eigenvalue in a minimisation is plotted: ' + JSON.stringify(sn && sn.values));
  // condensed reading (huge outputs) keeps the rows
  const lines = read('g16_sn2_ts_optfreq.log').split('\n'), id = BM.detectFormat('g16_sn2_ts_optfreq.log', lines.slice(0, 400).join('\n'));
  const cond = BM.parseLines('g16_sn2_ts_optfreq.log', BM.condense.lines(id, lines));
  eq(cond.frames.map(f => f.opt && f.opt.negEig).join(), '1,1,1,1', 'the same in condensed reading');
  // Q-Chem prints "Hessian Eigenvalues:" before the convergence table of each step
  const qc = read('qc/qchem_water_optfreq.out');
  const block = vals => '  Hessian Eigenvalues:\n' + vals.map(v => v.toFixed(6).padStart(13)).join('') + '\n\n Transition State search - taking P-RFO step\n\n';
  const qts = qc.replace('JOBTYPE         opt', 'JOBTYPE         ts').replace(/(   Maximum     Tolerance    Cnvgd\?)/g, (m0, m1, off) => block([-0.041234, 0.160000, 0.452011]) + m1);
  const dq = BM.parseFile('ts.out', qts), qs = BM.analysis.convergenceSeries(dq, {}).series.find(x => x.key === 'negEig');
  ok(dq.meta.jobTypes.indexOf('ts') >= 0 && qs && qs.values.join() === '1,1,1', 'Q-Chem TS search: one negative eigenvalue at each step: ' + JSON.stringify(qs && qs.values));
  const qmin = qc.replace(/(   Maximum     Tolerance    Cnvgd\?)/g, m1 => block([0.160000, 0.452011, 0.604808]) + m1);
  const dm = BM.parseFile('min.out', qmin);
  ok(dm.frames.filter(f => f.opt && f.opt.negEig === 0).length === 3 && !BM.analysis.convergenceSeries(dm, {}).series.some(x => x.key === 'negEig'), 'Q-Chem minimisation: counted (0), not plotted');
});

console.log('\nBlueMolecule normal-mode tool tests: ' + pass + ' passed, ' + fail + ' failed');
if (fail) { for (const f of failures) console.log('  FAIL ' + f); process.exit(1); }
