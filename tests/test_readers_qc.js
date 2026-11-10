#!/usr/bin/env node
/* BlueMolecule tests for src/readers-qc.js (GAMESS, NWChem, Q-Chem, Psi4, xtb, MOPAC, CP2K, Turbomole).
 * Run:  node tests/test_readers_qc.js */
'use strict';
const fs = require('fs');
const path = require('path');
const ROOT = path.join(__dirname, '..');
for (const f of ['elements.js', 'parsers.js', 'writers.js', 'readers-qc.js']) require(path.join(ROOT, 'src', f));
const BM = globalThis.BM;
const QC = path.join(__dirname, 'fixtures', 'qc');
const read = f => fs.readFileSync(path.join(QC, f), 'utf8');
const load = f => BM.parseFile(f, read(f));

let pass = 0, fail = 0, current = '';
const failures = [];
function ok(c, msg) { if (c) pass++; else { fail++; failures.push('[' + current + '] ' + msg); } }
function eq(a, b, msg) { ok(a === b, msg + ': got ' + JSON.stringify(a) + ', expected ' + JSON.stringify(b)); }
function near(a, b, tol, msg) { ok(typeof a === 'number' && Math.abs(a - b) <= tol, msg + ': got ' + a + ', expected ' + b); }
function arrNear(a, b, tol, msg) {
  if (!a || a.length !== b.length) { ok(false, msg + ': length ' + (a && a.length) + ' vs ' + b.length); return; }
  for (let i = 0; i < b.length; i++) if (!(Math.abs(a[i] - b[i]) <= tol)) { ok(false, msg + '[' + i + ']: got ' + a[i] + ', expected ' + b[i]); return; }
  pass++;
}
function test(name, fn) { current = name; try { fn(); } catch (e) { fail++; failures.push('[' + name + '] threw: ' + (e && e.stack || e)); } }

// Reference water optimisation used by the fixtures (3 steps, converged at the last one).
const E3 = [-76.4089074962, -76.4098631177, -76.409871723];
const GEO = [[0, 0, 0.1173, 0, 0.7572, -0.4692, 0, -0.7572, -0.4692], [0, 0, 0.1189, 0, 0.7598, -0.4756, 0, -0.7598, -0.4756],
  [0, 0, 0.119018, 0, 0.759242, -0.475927, 0, -0.759242, -0.475927]];
const FREQ = [1626.8, 3773.58, 3884.37], IR = [69.334, 4.187, 36.912];
const KEYS = ['program', 'meta', 'frames', 'vibrations', 'excitedStates', 'charges', 'dipole', 'orbitals', 'thermo', 'jobs', 'warnings'];

function common(ds, id, prog) {
  eq(ds.format, id, 'format id');
  ok(KEYS.every(k => k in ds), 'dataset keys');
  ok(new RegExp(prog).test(ds.program || ''), 'program "' + ds.program + '"');
  ok(Array.isArray(ds.jobs) && ds.jobs.length >= 1 && ds.jobs.every(j => j.frameStart >= 0 && j.frameEnd < ds.frames.length), 'job list');
  ds.frames.forEach(f => ok(f.xyz instanceof Float64Array && f.xyz.length === 3 * f.Z.length, 'frame arrays'));
}
function optSeries(ds, keys, thrKeys) {
  ds.frames.forEach((f, k) => {
    ok(f.opt && keys.every(q => Number.isFinite(f.opt[q])), 'step ' + (k + 1) + ' has ' + keys.join('/'));
    ok(thrKeys.every(q => Number.isFinite(f.opt[q + 'Thr']) && typeof f.opt[q + 'Conv'] === 'boolean'), 'step ' + (k + 1) + ' thresholds/flags');
  });
}
function vibCheck(ds, freqs, ir, tolF, tolI) {
  ok(ds.vibrations, 'vibrations present');
  arrNear(ds.vibrations.freqs, freqs, tolF, 'frequencies');
  if (ir) arrNear(ds.vibrations.ir, ir, tolI, 'IR intensities (km/mol)');
  const N = ds.frames[ds.vibrations.frameIndex].Z.length;
  ok(ds.vibrations.modes.length === freqs.length && ds.vibrations.modes.every(m => m instanceof Float64Array && m.length === 3 * N), 'mode shapes (3N each)');
}

// ===================================================================== GAMESS
test('GAMESS-US optimisation + Hessian', () => {
  const t = read('gamess_water_opt_hess.log');
  eq(BM.detectFormat('water.log', t), 'gamess-out', 'detect .log');
  eq(BM.detectFormat('water.out', t), 'gamess-out', 'detect .out');
  const ds = BM.parseFile('water.log', t);
  common(ds, 'gamess-out', '^GAMESS \\(US\\) 30 SEP 2021');
  eq(ds.frames.length, 3, 'NSERCH 0..2 (input and final re-print merged)');
  arrNear(ds.frames.map(f => f.energy), E3, 1e-10, 'energies');
  ds.frames.forEach((f, k) => arrNear(Array.from(f.xyz), GEO[k], 1e-6, 'geometry ' + k));
  optSeries(ds, ['maxForce', 'rmsForce'], ['maxForce', 'rmsForce']);
  near(ds.frames[0].opt.maxForceThr, 1e-4, 1e-12, 'OPTTOL');
  near(ds.frames[0].opt.rmsForceThr, 1e-4 / 3, 1e-12, 'OPTTOL/3');
  eq(JSON.stringify(ds.frames.map(f => f.converged)), '[false,false,true]', 'EQUILIBRIUM GEOMETRY LOCATED');
  eq(ds.frames[2].opt.maxForceConv && ds.frames[2].opt.rmsForceConv, true, 'converged flags');
  eq(ds.meta.method + '/' + ds.meta.basis, 'B3LYP/6-31G(d)', 'method/basis');
  eq(ds.meta.charge + ' ' + ds.meta.multiplicity, '0 1', 'charge/multiplicity');
  eq(JSON.stringify(ds.meta.jobTypes), '["opt","freq"]', 'job types');
  eq(ds.meta.normalTermination, true, 'terminated normally');
  vibCheck(ds, FREQ, IR, 1e-9, 0.01);
  eq(ds.vibrations.frameIndex, 2, 'frequencies at the optimised geometry');
  arrNear(Array.from(ds.vibrations.modes[2]), [0, -0.0681, 0, 0, 0.5382, -0.4531, 0, 0.5382, 0.4531], 1e-9, 'mode 9 vector');
  arrNear(ds.charges.Mulliken, [-0.684245, 0.342122, 0.342122], 1e-9, 'Mulliken');
  arrNear(ds.charges.Loewdin, [-0.453412, 0.226706, 0.226706], 1e-9, 'Lowdin');
  near(ds.dipole.total, 2.093745, 1e-9, 'dipole');
  eq(ds.orbitals.alpha.nOcc, 5, 'occupied orbitals');
  arrNear(ds.orbitals.alpha.energies, [-19.1351, -0.9971, -0.5203, -0.3662, -0.2921, 0.0654, 0.1491], 1e-9, 'orbital energies');
  near(ds.thermo.zpe, 0.021245, 1e-9, 'ZPE');
  near(ds.thermo.temperature, 298.15, 1e-9, 'temperature');
  near(ds.thermo.hCorr, 15.703 / 627.509474, 1e-9, 'H correction (kcal/mol table)');
  near(ds.thermo.G, E3[2] + 1.873 / 627.509474, 1e-9, 'G = E + G correction');
});

test('GAMESS input ($DATA C1)', () => {
  const t = read('gamess_water.inp');
  eq(BM.detectFormat('water.inp', t), 'gamess-input', 'detect');
  const ds = BM.parseFile('water.inp', t);
  arrNear(Array.from(ds.frames[0].xyz), GEO[0], 1e-9, 'coordinates');
  eq(ds.title, 'water', 'title');
  eq(ds.meta.method + '/' + ds.meta.basis, 'B3LYP/6-31G(d)', 'method/basis');
  let msg = '';
  try { BM.parseFile('w.inp', t.replace('\nC1\n', '\nCNV 2\n\n')); } catch (e) { msg = e.message; }
  ok(/C1/.test(msg), 'non-C1 symmetry explained: ' + msg);
  const w = BM.inputPrograms.find(p => p.id === 'gamess').write({ Z: [8, 1, 1], xyz: GEO[0], title: 'rt', task: 'opt' }).text;
  arrNear(Array.from(BM.parseFile('rt.inp', w).frames[0].xyz), GEO[0], 1e-8, 'writer round trip');
});

// ===================================================================== NWChem
test('NWChem optimisation + frequencies', () => {
  const t = read('nwchem_water_optfreq.out');
  eq(BM.detectFormat('water.out', t), 'nwchem-out', 'detect');
  const ds = BM.parseFile('water.out', t);
  common(ds, 'nwchem-out', '^NWChem 7\\.0\\.2');
  eq(ds.frames.length, 3, 'steps (final geometry merged)');
  arrNear(ds.frames.map(f => f.energy), E3, 1e-9, 'energies');
  optSeries(ds, ['maxForce', 'rmsForce'], ['maxForce', 'rmsForce']);
  ok(ds.frames.slice(1).every(f => ['energyChange', 'rmsDisp', 'maxDisp'].every(q => Number.isFinite(f.opt[q]))), 'Delta E, Xrms, Xmax from step 1');
  near(ds.frames[1].opt.maxDisp, 0.02367, 1e-12, 'Xmax');
  near(ds.frames[1].opt.maxDispThr, 0.0018, 1e-12, 'xmax threshold');
  near(ds.frames[1].opt.energyChange, -9.6e-4, 1e-12, 'Delta E (D exponent)');
  eq(JSON.stringify(ds.frames.map(f => f.converged)), '[false,false,true]', 'converged');
  ok(['maxForce', 'rmsForce', 'maxDisp', 'rmsDisp'].every(q => ds.frames[2].opt[q + 'Conv'] === true), 'all criteria met');
  eq(ds.meta.method + '/' + ds.meta.basis, 'B3LYP/6-31g*', 'method/basis');
  eq(JSON.stringify(ds.meta.jobTypes), '["opt","freq"]', 'job types');
  vibCheck(ds, FREQ, IR, 1e-9, 1e-9);
  arrNear(Array.from(ds.vibrations.modes[1]), [0, 0, 0.0478, 0, 0.5835, -0.3795, 0, -0.5835, -0.3795], 1e-9, 'mode 8');
  arrNear(ds.charges.Mulliken, [-0.68, 0.34, 0.34], 1e-9, 'Mulliken (Z - population)');
  near(ds.dipole.z, -0.823745 * 2.541746, 1e-6, 'dipole (a.u. -> Debye)');
  eq(ds.orbitals.alpha.nOcc, 5, 'nOcc');
  near(ds.orbitals.alpha.energies[0], -19.13513, 1e-9, 'orbital energy');
  near(ds.thermo.zpe, 0.021245, 1e-9, 'ZPE');
  near(ds.thermo.gCorr, 0.025024 - 298.15 * 45.096 / 1000 / 627.509474, 1e-9, 'G correction from H - TS');
  eq(ds.meta.normalTermination, true, 'normal termination');
});

test('NWChem TDDFT and input', () => {
  const ds = load('nwchem_water_tddft.out');
  eq(ds.excitedStates.length, 3, 'roots');
  arrNear(ds.excitedStates.map(s => s.energyEV), [7.2961, 9.0971, 9.8418], 1e-9, 'energies (eV)');
  arrNear(ds.excitedStates.map(s => s.f), [0.01292, 0, 0.04502], 1e-9, 'oscillator strengths');
  eq(ds.excitedStates[0].label, 'Singlet-b2', 'label');
  eq(ds.excitedStates[2].transitions.length, 2, 'transitions');
  ok(ds.meta.jobTypes.indexOf('td') >= 0, 'td job');
  const nw = BM.inputPrograms.find(p => p.id === 'nwchem').write({ Z: [8, 1, 1], xyz: GEO[0], title: 'rt', task: 'opt', charge: 0, multiplicity: 1 }).text;
  eq(BM.detectFormat('rt.nw', nw), 'nwchem-input', 'detect .nw');
  const back = BM.parseFile('rt.nw', nw);
  arrNear(Array.from(back.frames[0].xyz), GEO[0], 1e-8, 'writer round trip');
  ok(back.meta.jobTypes.indexOf('opt') >= 0 && back.meta.charge === 0, 'task and charge');
});

// ===================================================================== Q-Chem
test('Q-Chem optimisation + frequency jobs', () => {
  const t = read('qchem_water_optfreq.out');
  eq(BM.detectFormat('water.out', t), 'qchem-out', 'detect');
  const ds = BM.parseFile('water.out', t);
  common(ds, 'qchem-out', '^Q-Chem 6\\.0\\.2');
  eq(ds.frames.length, 3, 'cycles');
  arrNear(ds.frames.map(f => f.energy), E3, 1e-10, 'energies');
  optSeries(ds, ['maxForce', 'maxDisp'], ['maxForce', 'maxDisp']);
  ok(ds.frames[0].opt.energyChange === undefined, 'cycle 1 energy change "*********" skipped');
  near(ds.frames[1].opt.energyChange, -0.000955621, 1e-12, 'energy change');
  near(ds.frames[1].opt.energyChangeThr, 1e-6, 1e-15, 'energy threshold');
  eq(JSON.stringify(ds.frames.map(f => f.converged)), '[false,false,true]', 'OPTIMIZATION CONVERGED');
  eq(ds.jobs.length, 2, 'two jobs');
  eq(JSON.stringify(ds.jobs.map(j => j.kinds)), '[["opt"],["freq"]]', 'job kinds');
  eq(ds.jobs[1].frameStart, 2, 'freq job at the optimised geometry');
  vibCheck(ds, FREQ, IR, 1e-9, 1e-9);
  arrNear(ds.vibrations.redMass, [1.0828, 1.0447, 1.0817], 1e-9, 'reduced masses');
  arrNear(ds.charges.Mulliken, [-0.851498, 0.425749, 0.425749], 1e-9, 'Mulliken');
  near(ds.dipole.total, 2.0937, 1e-9, 'dipole');
  eq(ds.orbitals.alpha.nOcc, 5, 'nOcc');
  eq(ds.orbitals.alpha.energies.length, 7, 'orbital energies');
  near(ds.thermo.hCorr, 15.703 / 627.509474, 1e-9, 'enthalpy correction');
  eq(ds.meta.charge + ' ' + ds.meta.multiplicity + ' ' + ds.meta.method + ' ' + ds.meta.basis, '0 1 B3LYP 6-31G*', 'meta');
  eq(ds.meta.normalTermination, true, 'normal termination');
});

test('Q-Chem TDDFT and input', () => {
  const ds = load('qchem_water_tddft.out');
  arrNear(ds.excitedStates.map(s => s.energyEV), [7.2961, 9.0971, 9.8418], 1e-9, 'excitation energies');
  arrNear(ds.excitedStates.map(s => s.f), [0.0129, 0, 0.045], 1e-9, 'strengths');
  eq(ds.excitedStates[0].label, 'Singlet', 'multiplicity');
  eq(JSON.stringify(ds.excitedStates[0].transitions[0]), '{"from":"5","to":"6","coeff":0.9963}', 'V(1) -> orbital 6');
  ok(ds.meta.jobTypes.indexOf('td') >= 0, 'td job');
  const w = BM.inputPrograms.find(p => p.id === 'qchem').write({ Z: [8, 1, 1], xyz: GEO[0], title: 'rt', task: 'opt' }).text;
  eq(BM.detectFormat('rt.in', w), 'qchem-input', 'detect .in');
  const back = BM.parseFile('rt.in', w);
  arrNear(Array.from(back.frames[0].xyz), GEO[0], 1e-8, 'writer round trip');
  ok(back.meta.jobTypes.indexOf('opt') >= 0 && back.meta.method === 'B3LYP', 'jobtype/method');
});

// ===================================================================== Psi4
test('Psi4 optimisation + frequencies', () => {
  const t = read('psi4_water_optfreq.out');
  eq(BM.detectFormat('output.dat', t), 'psi4-out', 'detect (any extension)');
  const ds = BM.parseFile('water.out', t);
  common(ds, 'psi4-out', '^Psi4 1\\.7');
  eq(ds.frames.length, 3, 'steps (final optimized geometry merged)');
  arrNear(ds.frames.map(f => f.energy), E3, 1e-10, 'energies');
  optSeries(ds, ['maxForce', 'rmsForce', 'maxDisp', 'rmsDisp'], ['maxForce', 'maxDisp']);
  ok(ds.frames[0].opt.energyChange === undefined, 'step 1 Delta E (= total energy) skipped');
  near(ds.frames[1].opt.energyChange, -9.56e-4, 1e-12, 'Delta E');
  ok(ds.frames[0].opt.rmsForceThr === undefined, 'inactive criterion has no threshold');
  eq(ds.frames[2].opt.maxForceConv, true, 'met (*)');
  eq(ds.frames[1].opt.maxForceConv, false, 'unmet');
  eq(ds.frames[2].converged, true, 'Optimization is complete');
  vibCheck(ds, FREQ, IR, 1e-9, 1e-9);
  eq(JSON.stringify(ds.vibrations.sym), '["A1","A1","B2"]', 'irreps');
  arrNear(ds.charges.Mulliken, [-0.86158, 0.43079, 0.43079], 1e-9, 'Mulliken');
  arrNear(ds.charges.Loewdin, [-0.63754, 0.31877, 0.31877], 1e-9, 'Lowdin');
  near(ds.dipole.z, 0.8237763 * 2.541746, 1e-6, 'dipole');
  eq(ds.orbitals.alpha.nOcc, 5, 'nOcc');
  near(ds.thermo.G, -76.40627396, 1e-9, 'G');
  near(ds.thermo.zpe, 0.0212443, 1e-9, 'ZPE');
  eq(ds.meta.charge + ' ' + ds.meta.multiplicity + ' ' + ds.meta.method + ' ' + ds.meta.basis, '0 1 B3LYP 6-31G*', 'meta');
  eq(ds.meta.normalTermination, true, 'exited successfully');
});

test('Psi4 input', () => {
  const w = BM.inputPrograms.find(p => p.id === 'psi4').write({ Z: [8, 1, 1], xyz: GEO[0], title: 'rt', task: 'optfreq' }).text;
  eq(BM.detectFormat('rt.dat', w), 'psi4-input', 'detect .dat');
  const back = BM.parseFile('rt.dat', w);
  arrNear(Array.from(back.frames[0].xyz), GEO[0], 1e-8, 'writer round trip');
  ok(back.meta.jobTypes.indexOf('opt') >= 0 && back.meta.jobTypes.indexOf('freq') >= 0, 'optimize + frequency');
});

// ===================================================================== xtb
test('xtb --ohess', () => {
  const t = read('xtb_water_ohess.out');
  eq(BM.detectFormat('xtb.out', t), 'xtb-out', 'detect');
  const ds = BM.parseFile('xtb.out', t);
  common(ds, 'xtb-out', '^xtb 6\\.5\\.1');
  eq(ds.frames.length, 1, 'final structure only (xtb prints no per-cycle geometry)');
  near(ds.frames[0].energy, -5.070544440612, 1e-12, 'total energy');
  eq(ds.optSteps.length, 3, 'optimisation cycles in ds.optSteps');
  arrNear(ds.optSteps.map(s => s.energy), [-5.0697251, -5.0705308, -5.0705444], 1e-9, 'cycle energies');
  ok(ds.optSteps.every(s => ['energyChange', 'gradNorm', 'maxDisp'].every(q => Number.isFinite(s.opt[q]))), 'cycle criteria');
  near(ds.optSteps[0].opt.gradNormThr, 1e-3, 1e-12, 'gradient threshold');
  near(ds.optSteps[0].opt.energyChangeThr, 5e-6, 1e-12, 'energy threshold');
  eq(ds.frames[0].converged, true, 'converged');
  eq(ds.frames[0].opt.gradNormConv, true, 'final gradient criterion met');
  eq(ds.meta.method, 'GFN2-xTB', 'method');
  eq(ds.meta.charge + ' ' + ds.meta.multiplicity, '0 1', 'charge/multiplicity');
  vibCheck(ds, [1539.85, 3643.12, 3651.51], [98.5, 6.17, 12.58], 1e-9, 1e-9);
  ok(ds.warnings.some(w => /g98\.out/.test(w)), 'no-mode-vectors note');
  arrNear(ds.charges.Mulliken, [-0.565, 0.282, 0.282], 1e-9, 'charges');
  near(ds.dipole.total, 1.925, 1e-9, 'dipole');
  eq(ds.orbitals.alpha.nOcc, 4, 'HOMO index');
  near(ds.thermo.G, -5.061290115523, 1e-12, 'free energy');
  near(ds.thermo.temperature, 298.15, 1e-9, 'temperature');
  eq(ds.meta.normalTermination, true, 'normal termination');
  let msg = '';
  try { BM.parseFile('sp.out', t.slice(0, t.indexOf(' ================\n final structure'))); } catch (e) { msg = e.message; }
  ok(/no coordinates/.test(msg), 'output without coordinates explained');
});

// ===================================================================== MOPAC
test('MOPAC output, archive and input', () => {
  const t = read('mopac_water_opt_force.out');
  eq(BM.detectFormat('water.out', t), 'mopac-out', 'detect');
  const ds = BM.parseFile('water.out', t);
  common(ds, 'mopac-out', '^MOPAC 2016');
  eq(ds.frames.length, 2, 'input and final geometry');
  arrNear(Array.from(ds.frames[1].xyz), [0, 0, 0.11901813, 0, 0.75924211, -0.4759269, 0, -0.75924211, -0.4759269], 1e-8, 'final geometry');
  near(ds.frames[1].energy, -53.42174 / 627.509474, 1e-12, 'heat of formation (Hartree)');
  near(ds.meta.heatOfFormation, -53.42174, 1e-9, 'heat of formation (kcal/mol)');
  eq(ds.optSteps.length, 3, 'CYCLE lines');
  arrNear(ds.optSteps.map(s => s.opt.gradNorm), [26.917, 6.241, 0.358], 1e-9, 'GNORM per cycle');
  eq(ds.frames[1].converged && ds.frames[1].opt.gradNormConv, true, 'converged');
  eq(ds.meta.method + ' ' + ds.meta.charge + ' ' + ds.meta.multiplicity, 'PM7 0 1', 'keywords');
  ok(ds.meta.jobTypes.indexOf('opt') >= 0 && ds.meta.jobTypes.indexOf('freq') >= 0, 'OPT FORCE');
  vibCheck(ds, [1626.75, 3773.58, 3884.37], null, 1e-9);
  arrNear(Array.from(ds.vibrations.modes[0]), [0, 0, -0.068, 0, -0.4165, 0.5401, 0, 0.4165, 0.5401], 1e-9, 'mode 1');
  eq(JSON.stringify(ds.vibrations.sym), '["A1","A1","B2"]', 'symmetry');
  arrNear(ds.charges.Coulson, [-0.63772, 0.31886, 0.31886], 1e-9, 'net atomic charges');
  near(ds.dipole.z, 2.123, 1e-9, 'dipole vector');
  eq(ds.orbitals.alpha.nOcc, 4, 'filled levels');
  near(ds.orbitals.alpha.energies[3] * 27.211386245988, -12.18217, 1e-6, 'eigenvalues (eV -> Hartree)');
  eq(ds.meta.normalTermination, true, 'MOPAC DONE');
  const arc = load('mopac_water.arc');
  eq(arc.format, 'mopac-out', 'archive');
  arrNear(Array.from(arc.frames[0].xyz), Array.from(ds.frames[1].xyz), 1e-8, 'FINAL GEOMETRY OBTAINED');
  near(arc.frames[0].energy, -53.42174 / 627.509474, 1e-12, 'archive heat of formation');
  const mop = load('mopac_water_zmat.mop');
  eq(mop.format, 'mopac-input', '.mop');
  const x = mop.frames[0].xyz;
  near(Math.hypot(x[6] - x[0], x[7] - x[1], x[8] - x[2]), 0.959476, 1e-6, 'Z-matrix bond');
  eq(mop.meta.jobTypes[0], 'sp', '1SCF');
  const w = BM.inputPrograms.find(p => p.id === 'mopac').write({ Z: [8, 1, 1], xyz: GEO[0], title: 'rt', task: 'opt' }).text;
  arrNear(Array.from(BM.parseFile('rt.mop', w).frames[0].xyz), GEO[0], 1e-8, 'writer round trip');
});

// ===================================================================== CP2K
test('CP2K geometry optimisation and input', () => {
  const t = read('cp2k_water_geoopt.out');
  eq(BM.detectFormat('water.out', t), 'cp2k-out', 'detect');
  const ds = BM.parseFile('water.out', t);
  common(ds, 'cp2k-out', '^CP2K 9\\.1');
  eq(ds.frames.length, 1, 'starting coordinates');
  near(ds.frames[0].xyz[2], 5.1173, 1e-9, 'coordinates');
  eq(ds.optSteps.length, 3, 'Informations at step blocks');
  const o = ds.optSteps[1].opt;
  ok(['energyChange', 'maxDisp', 'rmsDisp', 'maxForce', 'rmsForce'].every(q => Number.isFinite(o[q])), 'all criteria');
  near(o.maxForceThr, 4.5e-4, 1e-12, 'gradient limit');
  eq(o.maxDispConv, true, 'step converged (YES)');
  eq(o.rmsDispConv, false, 'RMS step not converged (NO)');
  ok(['maxDisp', 'rmsDisp', 'maxForce', 'rmsForce'].every(q => ds.optSteps[2].opt[q + 'Conv'] === true), 'last step converged');
  near(ds.optSteps[2].energy, -17.154996221, 1e-10, 'step energy');
  eq(ds.meta.method + ' ' + ds.meta.basis + ' ' + ds.meta.charge + ' ' + ds.meta.multiplicity, 'PBE DZVP-MOLOPT-SR-GTH 0 1', 'meta');
  eq(JSON.stringify(ds.cell), '[[10,0,0],[0,10,0],[0,0,10]]', 'cell');
  arrNear(ds.charges.Mulliken, [-0.665234, 0.332617, 0.332617], 1e-9, 'Mulliken');
  eq(ds.meta.normalTermination, true, 'PROGRAM ENDED');
  const inp = load('cp2k_water.inp');
  eq(inp.format, 'cp2k-input', 'input');
  arrNear(Array.from(inp.frames[0].xyz), Array.from(ds.frames[0].xyz), 1e-9, '&COORD');
  eq(JSON.stringify(inp.cell), '[[10,0,0],[0,10,0],[0,0,10]]', '&CELL ABC');
  ok(inp.meta.jobTypes.indexOf('opt') >= 0 && inp.meta.method === 'PBE', 'RUN_TYPE / functional');
});

// ===================================================================== Turbomole
test('Turbomole coord and gradient', () => {
  const c = read('coord');
  eq(BM.detectFormat('coord', c), 'turbomole-coord', 'detect "coord"');
  const ds = BM.parseFile('coord', c);
  arrNear(Array.from(ds.frames[0].xyz), GEO[0], 1e-8, 'Bohr -> Angstrom');
  eq(JSON.stringify(ds.frames[0].Z), '[8,1,1]', 'lower-case elements, frozen flag');
  const g = read('gradient');
  eq(BM.detectFormat('gradient', g), 'turbomole-gradient', 'detect "gradient"');
  const gd = BM.parseFile('gradient', g);
  eq(gd.frames.length, 3, 'cycles');
  arrNear(gd.frames.map(f => f.energy), E3, 1e-10, 'energies');
  gd.frames.forEach((f, k) => arrNear(Array.from(f.xyz), GEO[k], 1e-8, 'geometry ' + k));
  arrNear(gd.frames.map(f => f.opt.maxForce), [0.012346, 0.001634, 0.000012], 1e-12, 'max gradient');
  near(gd.frames[0].opt.rmsForce, Math.sqrt((0.012346 ** 2 + 2 * (0.008765432 ** 2 + 0.006173 ** 2)) / 9), 1e-12, 'RMS gradient');
  arrNear(gd.frames.map(f => f.opt.gradNorm), [0.0208, 0.0027, 0.00002], 1e-12, '|dE/dxyz|');
  near(gd.frames[2].opt.energyChange, E3[2] - E3[1], 1e-12, 'energy change');
  const w = BM.inputPrograms.find(p => p.id === 'turbomole').write({ Z: [8, 1, 1], xyz: GEO[0] }).text;
  arrNear(Array.from(BM.parseFile('coord', w).frames[0].xyz), GEO[0], 1e-10, 'writer round trip');
});

// ===================================================================== detection
test('detection: Gaussian / ORCA never claimed, inputs weaker than outputs', () => {
  const FX = path.join(__dirname, 'fixtures');
  const gauss = fs.readFileSync(path.join(FX, 'g16_water_optfreq.log'), 'utf8');
  const orca = fs.readFileSync(path.join(FX, 'orca5_water_optfreq.out'), 'utf8');
  for (const ext of ['log', 'out', 'inp', 'in', 'dat']) {
    eq(BM.detectFormat('a.' + ext, gauss), 'gaussian-log', 'Gaussian .' + ext);
    eq(BM.detectFormat('a.' + ext, orca), 'orca-out', 'ORCA .' + ext);
  }
  // an ORCA output that ran xtb (contains the xtb banner) stays ORCA
  const xb = read('xtb_water_ohess.out').split('\n').slice(0, 12).join('\n');
  eq(BM.detectFormat('a.out', orca.replace('CARTESIAN COORDINATES (ANGSTROEM)', xb + '\nCARTESIAN COORDINATES (ANGSTROEM)')), 'orca-out', 'ORCA with embedded xtb banner');
  for (const f of BM.formatRegistry) ok(f.sniff(gauss.slice(0, 65536), 'a.log') === 0 && f.sniff(orca.slice(0, 65536), 'a.out') === 0, f.id + ' sniff is 0 for Gaussian/ORCA');
  const outs = { 'gamess_water_opt_hess.log': 'gamess-out', 'nwchem_water_optfreq.out': 'nwchem-out', 'qchem_water_optfreq.out': 'qchem-out',
    'psi4_water_optfreq.out': 'psi4-out', 'xtb_water_ohess.out': 'xtb-out', 'mopac_water_opt_force.out': 'mopac-out', 'cp2k_water_geoopt.out': 'cp2k-out' };
  for (const [f, id] of Object.entries(outs)) {
    const s = read(f).slice(0, 65536);
    const fmt = BM.formatRegistry.find(x => x.id === id);
    ok(fmt.sniff(s, f) >= 0.9, id + ' banner sniff >= 0.9');
    ok(BM.formatRegistry.filter(x => x.kind === 'structure').every(x => x.sniff(s, f) < 0.9), f + ': input sniffs stay below 0.9');
    eq(BM.detectFormat('noext', read(f)), id, f + ' without extension');
  }
  eq(BM.detectFormat('notes.txt', 'just some words\n'), null, 'unknown text');
  eq(BM.detectFormat('x.in', 'hello\n'), null, 'unknown .in');
});

console.log('\nBlueMolecule QC reader tests: ' + pass + ' passed, ' + fail + ' failed');
if (fail) { for (const f of failures) console.log('  FAIL ' + f); process.exit(1); }
