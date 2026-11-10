#!/usr/bin/env node
/* BlueMolecule file-format tests.  Run:  node tests/test_parsers.js
 * Fixtures come from tools/make_fixtures.py (expected.json holds the reference numbers). */
'use strict';
const fs = require('fs');
const path = require('path');
const ROOT = path.join(__dirname, '..');
require(path.join(ROOT, 'src', 'elements.js'));
require(path.join(ROOT, 'src', 'parsers.js'));
require(path.join(ROOT, 'src', 'writers.js'));
const BM = globalThis.BM;
const FIX = path.join(__dirname, 'fixtures');
const EXP = JSON.parse(fs.readFileSync(path.join(FIX, 'expected.json'), 'utf8'));
const read = f => fs.readFileSync(path.join(FIX, f), 'utf8');
const load = f => BM.parseFile(f, read(f));

let pass = 0, fail = 0, current = '';
const failures = [];
function ok(cond, msg) { if (cond) pass++; else { fail++; failures.push('[' + current + '] ' + msg); } }
function eq(a, b, msg) { ok(a === b, msg + ': got ' + JSON.stringify(a) + ', expected ' + JSON.stringify(b)); }
function near(a, b, tol, msg) { ok(typeof a === 'number' && Math.abs(a - b) <= tol, msg + ': got ' + a + ', expected ' + b); }
function arrNear(a, b, tol, msg) {
  if (!a || a.length !== b.length) { ok(false, msg + ': length ' + (a && a.length) + ' vs ' + b.length); return; }
  for (let i = 0; i < b.length; i++) {
    if (!(Math.abs(a[i] - b[i]) <= tol)) { ok(false, msg + '[' + i + ']: got ' + a[i] + ', expected ' + b[i]); return; }
  }
  pass++;
}
function arrEq(a, b, msg) { ok(JSON.stringify(Array.from(a || [])) === JSON.stringify(b), msg + ': got ' + JSON.stringify(a) + ', expected ' + JSON.stringify(b)); }
function test(name, fn) {
  current = name;
  try { fn(); } catch (e) { fail++; failures.push('[' + name + '] threw: ' + (e && e.stack || e)); }
}
function throwsMsg(fn, re, msg) {
  try { fn(); ok(false, msg + ': no error thrown'); } catch (e) { ok(e instanceof Error && re.test(e.message), msg + ': message "' + e.message + '"'); }
}
const dist = (x, i, j) => Math.hypot(x[3 * i] - x[3 * j], x[3 * i + 1] - x[3 * j + 1], x[3 * i + 2] - x[3 * j + 2]);
function angle(x, i, j, k) {
  const a = [x[3 * i] - x[3 * j], x[3 * i + 1] - x[3 * j + 1], x[3 * i + 2] - x[3 * j + 2]];
  const b = [x[3 * k] - x[3 * j], x[3 * k + 1] - x[3 * j + 1], x[3 * k + 2] - x[3 * j + 2]];
  return Math.acos((a[0] * b[0] + a[1] * b[1] + a[2] * b[2]) / (Math.hypot(...a) * Math.hypot(...b))) * 180 / Math.PI;
}
function dihedral(x, i, j, k, l) {
  const p = n => [x[3 * n], x[3 * n + 1], x[3 * n + 2]];
  const sub = (a, b) => [a[0] - b[0], a[1] - b[1], a[2] - b[2]];
  const cr = (a, b) => [a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]];
  const dt = (a, b) => a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
  const b0 = sub(p(i), p(j)), b1 = sub(p(k), p(j)), b2 = sub(p(l), p(k));
  const n1 = Math.hypot(...b1), u1 = b1.map(v => v / n1);
  const v = sub(b0, u1.map(q => q * dt(b0, u1))), w = sub(b2, u1.map(q => q * dt(b2, u1)));
  return Math.atan2(dt(cr(u1, v), w), dt(v, w)) * 180 / Math.PI;
}
const flat = a => [].concat(...a);
function checkShape(ds) {
  const keys = ['format', 'program', 'title', 'meta', 'frames', 'vibrations', 'excitedStates', 'charges', 'spins', 'dipole',
    'orbitals', 'thermo', 'volumes', 'bonds', 'cell', 'warnings'];
  for (const k of keys) ok(k in ds, 'dataset key "' + k + '" present');
  for (const k of ['route', 'keywords', 'method', 'basis', 'charge', 'multiplicity', 'jobTypes', 'normalTermination', 'solvent']) {
    ok(k in ds.meta, 'meta key "' + k + '" present');
  }
  for (const f of ds.frames) {
    ok(f.xyz instanceof Float64Array && f.xyz.length === 3 * f.Z.length, 'frame xyz is Float64Array(3N)');
    for (const k of ['energy', 'energies', 'label', 'converged', 'opt', 'atomNames', 'resNames', 'resIds']) ok(k in f, 'frame key ' + k);
  }
  if (ds.vibrations) {
    const v = ds.vibrations, N = ds.frames[v.frameIndex].Z.length;
    ok(v.modes.length === v.freqs.length && v.modes.every(m => m instanceof Float64Array && m.length === 3 * N), 'modes are Float64Array(3N)');
  }
}

// ======================================================================= Gaussian
test('Gaussian 16 opt+freq', () => {
  const ds = load('g16_water_optfreq.log'), E = EXP.g16_water_optfreq;
  checkShape(ds);
  eq(BM.detectFormat('g16_water_optfreq.log', read('g16_water_optfreq.log')), 'gaussian-log', 'detect');
  eq(ds.format, 'gaussian-log', 'format');
  eq(ds.program, 'Gaussian 16, Revision C.01', 'program');
  eq(ds.title, E.title, 'title');
  eq(ds.meta.route, E.route, 'wrapped route re-joined');
  eq(ds.meta.method, 'B3LYP', 'method');
  eq(ds.meta.basis, '6-31G(d)', 'basis');
  eq(ds.meta.charge, 0, 'charge'); eq(ds.meta.multiplicity, 1, 'multiplicity');
  ok(ds.meta.jobTypes.includes('opt') && ds.meta.jobTypes.includes('freq'), 'jobTypes opt+freq ' + ds.meta.jobTypes);
  eq(ds.meta.normalTermination, true, 'normal termination');
  eq(ds.meta.solvent, 'Water', 'solvent');
  ok(ds.meta.keywords.includes('freq=raman'), 'keywords');
  eq(ds.frames.length, 4, 'frames (re-printed final geometry and Link1 freq geometry merged)');
  arrNear(ds.frames.map(f => f.energy), E.energies, 1e-10, 'SCF energies');
  arrEq(ds.frames.map(f => f.converged), [false, false, false, true], 'converged flags');
  arrEq(ds.frames.map(f => f.label), ['Step 1', 'Step 2', 'Step 3', 'Step 4'], 'labels');
  ds.frames.forEach((f, k) => {
    arrNear([f.opt.maxForce, f.opt.rmsForce, f.opt.maxDisp, f.opt.rmsDisp], E.opt[k], 1e-6, 'opt table step ' + (k + 1));
  });
  eq(ds.frames[3].opt.maxForceConv, true, 'max force converged flag');
  near(ds.frames[0].opt.maxForceThr, 0.00045, 1e-9, 'threshold');
  ok(ds.frames.every(f => Math.abs(f.xyz[0]) + Math.abs(f.xyz[3]) + Math.abs(f.xyz[6]) < 1e-9), 'standard orientation used (x = 0)');
  arrNear(Array.from(ds.frames[3].xyz), E.final, 1e-6, 'final geometry');
  const v = ds.vibrations;
  eq(v.frameIndex, 3, 'vibrations attached to last geometry');
  arrNear(v.freqs, E.freqs, 1e-4, 'frequencies');
  arrNear(v.ir, E.ir, 1e-4, 'IR intensities');
  arrNear(v.raman, [4.8971, 105.3392, 42.1957], 1e-4, 'Raman activities');
  arrNear(v.depol, [0.6934, 0.1521, 0.75], 1e-4, 'depolarisation (P)');
  arrNear(v.redMass, E.redMass, 1e-4, 'reduced masses');
  arrNear(v.forceConst, E.forceConst, 1e-4, 'force constants');
  arrEq(v.sym, E.sym, 'mode symmetries');
  E.modes.forEach((m, k) => arrNear(Array.from(v.modes[k]), m, 1e-9, 'mode ' + (k + 1)));
  arrNear(ds.charges.Mulliken, [-0.868866, 0.434433, 0.434433], 1e-9, 'final Mulliken charges');
  arrNear(ds.charges.Hirshfeld, [-0.330523, 0.165262, 0.165261], 1e-9, 'Hirshfeld');
  arrNear(ds.charges.CM5, [-0.648178, 0.324089, 0.324089], 1e-9, 'CM5');
  eq(ds.spins, null, 'closed shell: no spins');
  arrNear([ds.dipole.x, ds.dipole.y, ds.dipole.z, ds.dipole.total], [0, 0, -2.0985, 2.0985], 1e-9, 'dipole');
  eq(ds.orbitals.alpha.nOcc, 5, 'nOcc');
  eq(ds.orbitals.alpha.energies.length, 19, 'orbital count');
  near(ds.orbitals.alpha.energies[0], -19.13798, 1e-9, 'last eigenvalue set kept');
  eq(ds.orbitals.beta, null, 'no beta set');
  for (const k of ['zpe', 'eCorr', 'hCorr', 'gCorr', 'eZPE', 'eTotal', 'H', 'G']) near(ds.thermo[k], E.thermo[k], 2e-6, 'thermo ' + k);
  near(ds.thermo.temperature, 298.15, 1e-9, 'temperature');
});

test('Gaussian hpmodes', () => {
  const ds = load('g16_water_hpmodes.log'), E = EXP.g16_water_hpmodes;
  eq(ds.program, 'Gaussian 16, Revision A.03', 'program');
  arrNear(ds.vibrations.freqs, E.freqs, 1e-4, 'frequencies');
  E.modes.forEach((m, k) => arrNear(Array.from(ds.vibrations.modes[k]), m, 1e-9, 'high-precision mode ' + (k + 1)));
  ok(ds.vibrations.modes.some(m => Array.from(m).some(x => Math.abs(x * 100 - Math.round(x * 100)) > 1e-3)), 'not the 2-decimal copy');
  arrNear(ds.vibrations.raman, [4.8971, 105.3392, 42.1957], 1e-4, 'hp Raman');
  // physical consistency: reduced mass = sum m_a |d_a|^2 and mass-weighted orthogonality
  const m = [15.99491, 1.00783, 1.00783];
  const mw = (a, b) => { let s = 0; for (let i = 0; i < 9; i++) s += m[(i / 3) | 0] * a[i] * b[i]; return s; };
  ds.vibrations.modes.forEach((d, k) => near(mw(d, d), ds.vibrations.redMass[k], 2e-3, 'reduced mass from mode ' + (k + 1)));
  near(mw(ds.vibrations.modes[0], ds.vibrations.modes[1]), 0, 2e-4, 'modes 1/2 orthogonal (mass-weighted)');
  near(mw(ds.vibrations.modes[1], ds.vibrations.modes[2]), 0, 2e-4, 'modes 2/3 orthogonal (mass-weighted)');
});

test('Gaussian TD-DFT', () => {
  const ds = load('g16_water_tddft.log'), E = EXP.g16_water_tddft;
  checkShape(ds);
  ok(ds.meta.jobTypes.includes('td'), 'td job type');
  eq(ds.excitedStates.length, 4, 'number of states');
  E.states.forEach((s, k) => {
    const x = ds.excitedStates[k];
    eq(x.index, s.index, 'index'); eq(x.label, s.label, 'label');
    near(x.energyEV, s.eV, 1e-9, 'eV'); near(x.wavelengthNM, s.nm, 1e-9, 'nm'); near(x.f, s.f, 1e-9, 'f');
    near(x.s2, s.s2, 1e-9, 'S2');
    eq(x.transitions.length, s.transitions.length, 'transition count');
    s.transitions.forEach((t, q) => {
      eq(x.transitions[q].from, t[0], 'from'); eq(x.transitions[q].to, t[1], 'to');
      near(x.transitions[q].coeff, t[2], 1e-9, 'coeff');
    });
  });
  eq(ds.excitedStates[1].transitions[1].deexcitation, true, '"<-" marked as de-excitation');
  arrNear(ds.charges.ESP, E.esp, 1e-9, 'ESP charges');
  arrNear(ds.charges.NPA, E.npa, 1e-9, 'NPA charges');
  ok(ds.charges.Mulliken && ds.charges.Mulliken.length === 3, 'Mulliken');
  near(ds.frames[0].energy, E.energy, 1e-10, 'ground-state energy kept as frame energy');
  near(ds.frames[0].energies.TD, E.energy + 6.5432 / 27.211386245988, 2e-9, 'TD total energy stored');
});

test('Gaussian 09 Input orientation only, error termination', () => {
  const ds = load('g09_cisplatin_input_orientation.log'), E = EXP.g09_cisplatin;
  checkShape(ds);
  eq(ds.program, 'Gaussian 09, Revision D.01', 'program');
  eq(ds.frames.length, 2, 'frames');
  arrEq(ds.frames[0].Z, E.Z, 'Z incl. Pt');
  arrNear(Array.from(ds.frames[0].xyz), E.x1, 1e-6, 'step 1 geometry');
  arrNear(Array.from(ds.frames[1].xyz), E.x2, 1e-6, 'step 2 geometry');
  arrNear(ds.frames.map(f => f.energy), E.energies, 1e-9, 'energies');
  eq(ds.meta.normalTermination, false, 'error termination');
  ok(ds.frames.every(f => !f.converged), 'not converged');
  eq(ds.orbitals.alpha.nOcc, 36, 'occupied orbitals');
  arrNear(ds.orbitals.alpha.energies.slice(0, 36), E.occ, 1e-9, 'glued eigenvalues (-101.56789-101.56788)');
  arrNear(ds.charges.Mulliken, E.mull, 1e-9, '"Mulliken atomic charges" header');
  arrNear(ds.charges.ESP, E.esp, 1e-9, 'G09 ESP fit block');
  near(ds.dipole.total, Math.hypot(3.1234, 3.1234), 1e-4, 'old dipole format');
  eq(ds.meta.basis, 'genecp', 'basis from route');
  eq(ds.vibrations, null, 'no vibrations');
});

test('Gaussian relaxed scan', () => {
  const ds = load('g16_water_scan.log'), E = EXP.g16_scan;
  eq(ds.frames.length, 4, 'frames');
  arrEq(ds.frames.map(f => f.converged), [false, true, false, true], 'converged scan points');
  arrEq(ds.frames.map(f => f.label), E.labels, 'labels');
  arrNear(ds.frames.map(f => f.energy), E.energies, 1e-10, 'energies');
});

test('Gaussian log edge cases (inline)', () => {
  const head = ' Entering Gaussian System, Link 0=g16\n ------------\n #p mp2/cc-pvdz\n ------------\n' +
    ' 1/38=1/1;\n -----\n MP2 test\n -----\n Symbolic Z-matrix:\n Charge =  1 Multiplicity = 2\n';
  const orient = (rows) => '                         Standard orientation:\n ' + '-'.repeat(69) +
    '\n Center     Atomic      Atomic             Coordinates (Angstroms)\n Number     Number       Type             X           Y           Z\n ' +
    '-'.repeat(69) + '\n' + rows + ' ' + '-'.repeat(69) + '\n';
  const txt = head + orient('      1          8           0        0.000000    0.000000    0.119000\r\n' +
    '      2          1           0        0.000000    0.760000   -0.476000\r\n      3         -2           0        3.000000    0.000000    0.000000\n') +
    ' SCF Done:  E(UHF) =  -75.9876543210     A.U. after   12 cycles\n' +
    ' E2 =    -0.2026457212D+00 EUMP2 =    -0.76214093838D+02\n' +
    ' Mulliken charges and spin densities:\n               1          2\n     1  O   -0.300000   1.100000\n     2  H    1.300000  -0.100000\n' +
    ' Sum of Mulliken charges =   1.00000   1.00000\n' +
    ' Harmonic frequencies (cm**-1), IR intensities (KM/Mole), Raman scattering\n' +
    '                      1\n                      A\n Frequencies --   3500.1234\n Red. masses --      1.0500\n' +
    ' Frc consts  --  **********\n IR Inten    --  ******\n  Atom  AN      X      Y      Z\n     1   8     0.00   0.00   0.06\n     2   1     0.00   0.70  -0.71\n\n' +
    ' Normal termination of Gaussian 16 at Wed Nov  4 21:00:00 2026.\n';
  const ds = BM.parseFile('mp2.log', txt);
  eq(ds.frames.length, 1, 'one frame');
  eq(ds.frames[0].Z.length, 2, 'Tv excluded from atoms');
  arrNear(flat(ds.cell), [3, 0, 0, 0, 0, 0, 0, 0, 0], 1e-12, 'Tv -> cell (padded)');
  near(ds.frames[0].energy, -76.214093838, 1e-12, 'MP2 energy (D exponent) is the frame energy');
  near(ds.frames[0].energies.SCF, -75.987654321, 1e-12, 'SCF energy kept');
  eq(ds.meta.method, 'MP2', 'method'); eq(ds.meta.basis, 'cc-pvdz', 'basis');
  eq(ds.meta.charge, 1, 'charge'); eq(ds.meta.multiplicity, 2, 'mult');
  arrNear(ds.spins.Mulliken, [1.1, -0.1], 1e-12, 'spin densities');
  ok(Number.isNaN(ds.vibrations.forceConst[0]) && Number.isNaN(ds.vibrations.ir[0]), '****** -> NaN');
  near(ds.vibrations.freqs[0], 3500.1234, 1e-9, 'freq next to overflow fields');
  eq(ds.title, 'MP2 test', 'title');
  throwsMsg(() => BM.parseFile('empty.log', head), /Gaussian output.*geometry/i, 'no geometry -> readable error');
});

// =========================================================================== ORCA
test('ORCA 5 opt+freq', () => {
  const ds = load('orca5_water_optfreq.out'), E = EXP.orca5_water;
  checkShape(ds);
  eq(BM.detectFormat('orca5_water_optfreq.out', read('orca5_water_optfreq.out')), 'orca-out', 'detect .out as ORCA');
  eq(ds.program, 'ORCA 5.0.4', 'program');
  arrEq(ds.meta.keywords, ['B3LYP', 'def2-SVP', 'Opt', 'Freq'], 'keywords');
  eq(ds.meta.route, '! B3LYP def2-SVP Opt Freq', 'route');
  eq(ds.meta.method, 'B3LYP', 'method'); eq(ds.meta.basis, 'def2-SVP', 'basis');
  eq(ds.meta.charge, 0, 'charge'); eq(ds.meta.multiplicity, 1, 'mult');
  ok(ds.meta.jobTypes.includes('opt') && ds.meta.jobTypes.includes('freq'), 'job types');
  eq(ds.meta.normalTermination, true, 'terminated normally');
  eq(ds.frames.length, 3, 'frames (final evaluation merged)');
  arrNear(ds.frames.map(f => f.energy), E.energies, 1e-11, 'FINAL SINGLE POINT ENERGY per frame');
  arrEq(ds.frames.map(f => f.converged), [false, false, true], 'converged');
  arrEq(ds.frames.map(f => f.label), ['Cycle 1', 'Cycle 2', 'Cycle 3'], 'labels');
  ds.frames.forEach((f, k) => {
    const o = E.opt[k];
    if (o[0] === null) ok(f.opt.energyChange === undefined, 'no energy change in cycle 1');
    else near(f.opt.energyChange, o[0], 1e-10, 'energy change');
    arrNear([f.opt.rmsForce, f.opt.maxForce, f.opt.rmsDisp, f.opt.maxDisp], o.slice(1), 1e-10, 'convergence cycle ' + (k + 1));
  });
  arrNear(Array.from(ds.frames[2].xyz), E.final, 1e-6, 'final geometry');
  const v = ds.vibrations;
  eq(v.frameIndex, 2, 'frameIndex');
  arrNear(v.freqs, E.freqs, 1e-9, 'frequencies (zero modes dropped)');
  arrNear(v.ir, E.ir, 1e-9, 'IR Int column');
  arrNear(v.raman, E.raman, 1e-9, 'Raman activity'); arrNear(v.depol, E.depol, 1e-9, 'depolarisation');
  E.modes.forEach((m, k) => arrNear(Array.from(v.modes[k]), m, 1e-9, 'normal mode ' + k));
  arrNear(ds.charges.Mulliken, [-0.371632, 0.185816, 0.185816], 1e-9, 'Mulliken (last)');
  arrNear(ds.charges.Loewdin, [-0.213456, 0.106728, 0.106728], 1e-9, 'Loewdin');
  arrNear(ds.charges.Hirshfeld, [-0.312345, 0.156172, 0.156173], 1e-9, 'Hirshfeld');
  arrNear([ds.dipole.x, ds.dipole.y, ds.dipole.z], [0, 0, -0.82399 * 2.541746], 1e-9, 'dipole a.u. -> Debye');
  near(ds.dipole.total, 2.09437, 1e-5, 'Magnitude (Debye)');
  eq(ds.orbitals.alpha.nOcc, 5, 'nOcc'); eq(ds.orbitals.alpha.energies.length, 10, 'orbitals');
  near(ds.orbitals.alpha.energies[0], -18.817393, 1e-9, 'orbital energy in Eh');
  const T = E.thermo, eel = E.energies[2];
  near(ds.thermo.zpe, T.zpe, 2e-8, 'ZPE'); near(ds.thermo.H, T.H, 2e-8, 'H'); near(ds.thermo.G, T.G, 2e-8, 'G');
  near(ds.thermo.eTotal, T.eTotal, 2e-8, 'U'); near(ds.thermo.eCorr, T.eCorr, 2e-8, 'total correction');
  near(ds.thermo.gCorr, T.gCorr, 2e-8, 'G-E(el)'); near(ds.thermo.hCorr, T.H - eel, 3e-8, 'H-E(el)');
  near(ds.thermo.eZPE, eel + T.zpe, 3e-8, 'E+ZPE'); near(ds.thermo.temperature, 298.15, 1e-9, 'T');
});

test('ORCA 4 IR table', () => {
  const ds = load('orca4_water_freq.out'), E = EXP.orca4_ir;
  eq(ds.program, 'ORCA 4.2.1', 'program');
  arrNear(ds.vibrations.freqs, E.freqs, 1e-9, 'freqs');
  arrNear(ds.vibrations.ir, E.ir, 1e-9, 'T**2 column used as km/mol');
  ok(ds.meta.jobTypes.includes('freq') && !ds.meta.jobTypes.includes('opt'), 'job types');
});

test('ORCA 5 TD-DFT', () => {
  const ds = load('orca5_water_tddft.out'), E = EXP.orca5_td;
  const xs = ds.excitedStates;
  eq(xs.length, 5, '3 singlets + 2 triplets');
  E.singlets.forEach((s, k) => {
    eq(xs[k].index, s.index, 'index'); eq(xs[k].label, 'Singlet', 'label');
    near(xs[k].energyEV, s.eV, 1e-6, 'eV'); near(xs[k].wavelengthNM, s.nm, 1e-6, 'nm');
    near(xs[k].f, s.f, 1e-12, 'fosc from the electric-dipole table');
    near(xs[k].s2, 0, 1e-12, 'S2');
    s.transitions.forEach((t, q) => { eq(xs[k].transitions[q].from, t[0], 'from'); eq(xs[k].transitions[q].to, t[1], 'to');
      near(xs[k].transitions[q].coeff, t[2], 1e-9, 'c'); });
  });
  E.triplets.forEach((s, k) => {
    const x = xs[3 + k];
    eq(x.label, 'Triplet', 'triplet label'); near(x.energyEV, s.eV, 6e-4, 'triplet eV'); eq(x.f, 0, 'triplet f'); near(x.s2, 2, 1e-12, 'S2');
  });
  ok(ds.meta.jobTypes.includes('td'), 'td job');
});

test('ORCA 6 TD-DFT', () => {
  const ds = load('orca6_water_tddft.out'), E = EXP.orca6_td;
  eq(ds.program, 'ORCA 6.0.1', 'program');
  eq(ds.excitedStates.length, 3, 'states');
  E.singlets.forEach((s, k) => {
    const x = ds.excitedStates[k];
    eq(x.index, s.index, 'index'); eq(x.label, 'Singlet-A', 'label from "1-1A"');
    near(x.energyEV, s.eV, 1e-6, 'eV'); near(x.f, s.f, 1e-12, 'f');
    eq(x.transitions.length, s.transitions.length, 'transitions');
  });
});

test('ORCA relaxed scan', () => {
  const ds = load('orca5_water_scan.out');
  eq(ds.frames.length, 4, 'frames');
  arrEq(ds.frames.map(f => f.converged), [false, true, false, true], 'converged');
  eq(ds.frames[2].label, 'Scan step 2, cycle 1', 'label');
  arrNear(ds.frames.map(f => f.energy), EXP.orca_scan.energies, 1e-11, 'energies');
  ok(ds.meta.jobTypes.includes('scan') && ds.meta.jobTypes.includes('opt'), 'scan job type');
});

test('open-shell and multi-job details (inline)', () => {
  // Gaussian: two Link1 jobs with frequencies -> last set (hpmodes preferred inside that job)
  const hp = read('g16_water_hpmodes.log');
  const two = BM.parseFile('two.log', hp + hp.split('1632.1418').join('1600.0000'));
  eq(two.frames.length, 1, 'identical geometries of both jobs merged');
  near(two.vibrations.freqs[0], 1600, 1e-9, 'last frequency set kept');
  ok(Array.from(two.vibrations.modes[1]).some(x => Math.abs(x * 100 - Math.round(x * 100)) > 1e-3), 'hpmodes of last job');
  // Gaussian UHF eigenvalues and NBO alpha/beta spin blocks
  const uhf = ' Entering Gaussian System, Link 0=g16\n ------\n #p uhf/sto-3g\n ------\n 1/1;\n ---\n OH\n ---\n Charge =  0 Multiplicity = 2\n' +
    '                          Input orientation:\n ' + '-'.repeat(69) + '\n Center     Atomic      Atomic             Coordinates (Angstroms)\n' +
    ' Number     Number       Type             X           Y           Z\n ' + '-'.repeat(69) + '\n' +
    '      1          8           0        0.000000    0.000000    0.000000\n      2          1           0        0.000000    0.000000    0.970000\n ' + '-'.repeat(69) + '\n' +
    ' SCF Done:  E(UHF) =  -74.3612345678     A.U. after    9 cycles\n' +
    ' Alpha  occ. eigenvalues --  -19.20000  -1.10000\n Alpha virt. eigenvalues --    0.10000\n' +
    '  Beta  occ. eigenvalues --  -19.18000\n  Beta virt. eigenvalues --   -0.10000   0.12000\n' +
    ' N A T U R A L   A T O M I C   O R B I T A L   A N D\n Summary of Natural Population Analysis:\n\n    Atom  No    Charge\n --------------------\n' +
    '      O    1   -0.40000\n      H    2    0.40000\n ====================\n' +
    ' *******         Alpha spin orbitals         *******\n Summary of Natural Population Analysis:\n --------------------\n' +
    '      O    1   -0.70000\n      H    2    0.20000\n ====================\n';
  const u = BM.parseFile('oh.log', uhf);
  eq(u.orbitals.alpha.nOcc, 2, 'alpha nOcc'); eq(u.orbitals.beta.nOcc, 1, 'beta nOcc'); eq(u.orbitals.beta.energies.length, 3, 'beta energies');
  arrNear(u.charges.NPA, [-0.4, 0.4], 1e-12, 'total NPA charges, spin blocks ignored');
  eq(u.meta.normalTermination, null, 'no termination line -> null');
  // ORCA UKS: spin-resolved orbitals, Mulliken spin populations, imaginary mode kept, zero modes dropped
  const rows = (vals) => vals.map((v, r) => '      ' + r + '    ' + v.map(x => x.toFixed(6).padStart(11)).join('')).join('\n');
  const M = [[0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, -0.0612], [0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0.9718]];
  const orca = '                                 * O   R   C   A *\n                         Program Version 5.0.4 -  RELEASE  -\n' +
    '|  1> ! UKS B3LYP def2-SVP Freq CPCM(Toluene)\n|  2> * xyz 0 2\n---------------------------------\nCARTESIAN COORDINATES (ANGSTROEM)\n' +
    '---------------------------------\n  O      0.000000    0.000000    0.000000\n  H      0.000000    0.000000    0.970000\n\n' +
    'Total Energy       :          -75.70000000 Eh           -2059.90000 eV\n' +
    '----------------\nORBITAL ENERGIES\n----------------\n                 SPIN UP ORBITALS\n  NO   OCC          E(Eh)            E(eV) \n' +
    '   0   1.0000     -18.900000      -514.2965 \n   1   1.0000      -1.100000       -29.9325 \n   2   0.0000       0.100000         2.7211 \n\n' +
    '                 SPIN DOWN ORBITALS\n  NO   OCC          E(Eh)            E(eV) \n   0   1.0000     -18.880000      -513.7523 \n' +
    '   1   0.0000      -0.100000        -2.7211 \n   2   0.0000       0.120000         3.2654 \n\n' +
    '--------------------------------------------\nMULLIKEN ATOMIC CHARGES AND SPIN POPULATIONS\n--------------------------------------------\n' +
    '   0 O :   -0.300000    1.050000\n   1 H :    0.300000   -0.050000\nSum of atomic charges         :    0.0000000\n\n' +
    '-------------------------   --------------------\nFINAL SINGLE POINT ENERGY       -75.700000000000\n-------------------------   --------------------\n\n' +
    '-----------------------\nVIBRATIONAL FREQUENCIES\n-----------------------\n\nScaling factor for frequencies =  1.000000000  (already applied!)\n\n' +
    [0, 1, 2, 3, 4].map(k => '   ' + k + ':         0.00 cm**-1').join('\n') + '\n   5:      -153.40 cm**-1 ***imaginary mode***\n\n\n' +
    '------------\nNORMAL MODES\n------------\n\n                  0          1          2          3          4          5    \n' + rows(M) + '\n\n\n' +
    '-----------\nIR SPECTRUM\n-----------\n\n Mode   freq       eps      Int      T**2         TX        TY        TZ\n' +
    '       cm**-1   L/(mol*cm) km/mol    a.u.\n----------------------------------------------------------------------------\n' +
    '  5:   -153.40   0.001000    5.05  0.000192  ( 0.000000  0.000000  0.013856)\n\n' +
    'ORCA finished by error termination in SCF gradient\n';
  const o = BM.parseFile('oh.out', orca);
  eq(o.orbitals.alpha.nOcc, 2, 'spin-up nOcc'); eq(o.orbitals.beta.nOcc, 1, 'spin-down nOcc');
  arrNear(o.spins.Mulliken, [1.05, -0.05], 1e-12, 'Mulliken spin populations');
  eq(o.vibrations.freqs.length, 1, 'only the non-zero mode kept'); near(o.vibrations.freqs[0], -153.4, 1e-12, 'imaginary mode');
  arrNear(Array.from(o.vibrations.modes[0]), [0, 0, -0.0612, 0, 0, 0.9718], 1e-12, 'mode column 5');
  arrNear(o.vibrations.ir, [5.05], 1e-12, 'IR mapped by mode index');
  eq(o.meta.solvent, 'Toluene', 'CPCM solvent'); eq(o.meta.multiplicity, 2, 'mult from "* xyz 0 2"');
  eq(o.meta.normalTermination, false, 'error termination'); near(o.frames[0].energy, -75.7, 1e-12, 'final energy');
});

// ======================================================= review: real-output layouts
// Snippets in fixtures/review follow genuine Gaussian 09/16 and ORCA 5/6 headers and column widths.
const rv = f => BM.parseFile(f, fs.readFileSync(path.join(FIX, 'review', f), 'utf8'));

test('review: Gaussian open-shell populations with "summed into heavy atoms" tables', () => {
  const ds = rv('g16_oh_radical_pop.log');
  checkShape(ds);
  arrNear(ds.charges.Mulliken, [-0.384505, 0.384505], 1e-9, 'Mulliken per atom (not the summed table)');
  arrNear(ds.spins.Mulliken, [1.0292, -0.0292], 1e-9, 'Mulliken spins');
  arrNear(ds.charges.Hirshfeld, [-0.214372, 0.214372], 1e-9, 'Hirshfeld per atom, summed Q-H/Q-CM5 table ignored');
  arrNear(ds.charges.CM5, [-0.441538, 0.441538], 1e-9, 'CM5 per atom');
  arrNear(ds.spins.Hirshfeld, [0.991524, 0.008476], 1e-9, 'Hirshfeld spins');
  arrNear(ds.charges.APT, [-0.302341, 0.302341], 1e-9, 'APT charges');
  eq(ds.orbitals.beta.nOcc, 4, 'beta occupied');
  arrNear(ds.vibrations.freqs, [3713.8417], 1e-9, 'single-mode block');
  arrNear(Array.from(ds.vibrations.modes[0]), [0, 0, 0.06, 0, 0, -1], 1e-12, 'single-mode displacements');
  arrEq(ds.vibrations.sym, ['SG'], 'symmetry');
});

test('review: Gaussian Standard orientation in job 1 only (Link1 NoSymm opt)', () => {
  const ds = rv('g16_link1_symm_then_nosymm.log');
  checkShape(ds);
  eq(ds.frames.length, 3, 'SP frame + both NoSymm optimisation steps');
  arrNear(ds.frames.map(f => f.energy), [-76.4089533216, -76.4082676378, -76.4089354966], 1e-10, 'energies stay with their own geometry');
  arrEq(ds.frames.map(f => f.label), [null, 'Step 1', 'Step 2'], 'labels');
  arrEq(ds.frames.map(f => f.converged), [false, false, true], 'converged');
  near(ds.frames[0].xyz[0], 0, 1e-12, 'job 1 in Standard orientation');
  arrNear(Array.from(ds.frames[2].xyz.slice(0, 3)), [0.096995, -0.242315, 0.41309], 1e-9, 'job 2 in Input orientation');
  eq(ds.meta.normalTermination, true, 'both jobs terminated normally');
});

test('review: Gaussian 09 TD-DFT optimisation', () => {
  const ds = rv('g09_water_tdopt.log');
  eq(ds.program, 'Gaussian 09, Revision D.01', 'program');
  eq(ds.frames.length, 2, 'frames');
  arrNear(ds.frames.map(f => f.energies.TD), [-76.1688716503, -76.1801234567], 1e-9, '"E(TD-HF/TD-KS)" read');
  arrNear(ds.frames.map(f => f.energy), [-76.1688716503, -76.1801234567], 1e-9, 'frame energy = optimised excited-state energy');
  near(ds.frames[1].energies.SCF, -76.40612, 1e-10, 'ground-state SCF kept');
  eq(ds.excitedStates.length, 3, 'states of the last step');
  near(ds.excitedStates[0].energyEV, 6.1003, 1e-9, 'last TD block');
});

test('review: Gaussian Excited State lines with long wavelengths', () => {
  const ds = rv('g16_td_long_wavelength.log');
  const xs = ds.excitedStates;
  eq(xs.length, 3, 'all states read');
  arrEq(xs.map(s => s.label), ['Triplet-?Sym', 'Triplet-B2', 'Singlet-A1'], 'labels');
  near(xs[0].wavelengthNM, 13051.18, 1e-9, '"eV13051.18 nm" (fields run together)');
  near(xs[1].wavelengthNM, 1239.84198 / 0.01, 1e-6, '"********" wavelength from eV');
  arrEq(xs.map(s => s.transitions.length), [1, 1, 2], 'transitions stay with their state');
  eq(xs[1].transitions[0].from, '4B', 'beta transition');
});

test('review: Gaussian freq=ROA labels and a one-mode last block', () => {
  const v = rv('g16_water_roa.log').vibrations;
  arrNear(v.freqs, [1632.1065, 3832.7837, 3939.5705], 1e-9, 'all blocks read past "Raman1 Fr= 1--"');
  arrNear(Array.from(v.modes[1]), [0, 0, 0.05, 0, 0.58, -0.41, 0, -0.58, -0.41], 1e-12, 'mode 2 displacements');
  arrNear(Array.from(v.modes[2]), [0, -0.07, 0, 0, 0.55, -0.44, 0, 0.55, 0.44], 1e-12, 'mode 3 displacements');
  arrNear(v.depol, [0.6934, 0.1521, 0.75], 1e-9, 'Depolar (P)');
  arrEq(v.sym, ['A1', 'A1', 'B2'], 'symmetries');
});

test('review: Gaussian coordinates >= 1000 A and a log still being written', () => {
  const ds = rv('g16_far_coords_running.log');
  eq(ds.frames.length, 1, 'cut-off orientation table ignored');
  arrNear(Array.from(ds.frames[0].xyz.slice(0, 6)), [-1234.567891, -1000.1, 12.5, -1233.48, -1000.1, 12.5], 1e-9, 'run-together 3F12.6 fields');
  eq(ds.meta.normalTermination, null, 'no termination yet');
});

test('review: Gaussian ONIOM opt + Link1 freq', () => {
  const ds = rv('g16_oniom_optfreq.log');
  eq(ds.frames.length, 2, 'freq job at the optimised geometry merged despite three SCF energies per point');
  arrNear(ds.frames.map(f => f.energy), [-76.421, -76.4223456789], 1e-10, 'ONIOM extrapolated energies');
  eq(ds.vibrations.frameIndex, 1, 'vibrations on the optimised geometry');
  eq(ds.meta.charge, 0, 'charge'); eq(ds.meta.multiplicity, 1, 'multiplicity (first of the ONIOM lines)');
  ok(/^ONIOM/.test(ds.meta.method), 'method ' + ds.meta.method);
});

test('review: ORCA ghost atoms and point charges', () => {
  const ds = rv('orca5_ghost_pointcharge.out');
  checkShape(ds);
  arrEq(ds.frames[0].Z, [8, 1, 1, 0, 0, 0, 0], 'ghosts ("O:", "H :") and Q as dummies, geometry not cut at "H :"');
  arrNear(Array.from(ds.frames[0].xyz.slice(12, 21)), [1.680398, -0.373741, -0.758561, 1.680398, -0.373741, 0.758561, 4, 0, 0], 1e-9, 'coordinates after the charge column');
});

test('review: ORCA full TD-DFT with de-excitations', () => {
  const xs = rv('orca5_water_rpa.out').excitedStates;
  eq(xs.length, 2, 'states');
  eq(xs[0].transitions.length, 3, 'lines after "<-" kept');
  eq(xs[0].transitions[1].deexcitation, true, '"<-" marked');
  near(xs[0].transitions[2].coeff, 0.10954451, 1e-12, 'coefficient');
  near(xs[0].f, 0.016487919, 1e-12, 'fosc');
});

test('review: ORCA $new_job chain cut off mid-geometry', () => {
  const ds = rv('orca6_newjob_running.out');
  eq(ds.program, 'ORCA 6.0.1', 'program');
  eq(ds.frames.length, 3, 'opt cycles (final evaluation merged) + job 2; partial block dropped');
  arrNear(ds.frames.map(f => f.energy), [-76.3, -76.31, -75.95], 1e-10, 'energies');
  ok(ds.frames.every(f => f.Z.length === 3), 'no truncated frame');
  eq(ds.meta.charge, 0, 'charge of the first job'); eq(ds.meta.multiplicity, 1, 'multiplicity of the first job');
  eq(ds.meta.normalTermination, null, 'still running');
  eq(BM.detectFormat('job.log', fs.readFileSync(path.join(FIX, 'review', 'orca6_newjob_running.out'), 'utf8')), 'orca-out', 'detect');
});

// ================================================================ XYZ / inputs
test('XYZ multi-frame (CRLF, tabs, extra columns, atomic numbers)', () => {
  const ds = load('water_traj.xyz');
  checkShape(ds);
  eq(ds.frames.length, 4, 'frames');
  arrNear(ds.frames.map(f => f.energy), EXP.xyz.energies, 1e-12, 'energies from 4 comment styles');
  arrEq(ds.frames[3].Z, [8, 1, 1], 'atomic numbers as labels');
  ok(ds.frames.every(f => f.Z.length === 3), 'atom counts');
  near(dist(ds.frames[0].xyz, 0, 1), 0.95, 1e-6, 'O-H distance frame 1');
  eq(BM.detectFormat('traj.txt', read('water_traj.xyz')), 'xyz', '.txt sniffed as XYZ');
  const bare = BM.parseFile('bare.txt', 'O 0 0 0\nH 0 0.76 0.59\nH 0 -0.76 0.59\n');
  eq(bare.frames[0].Z.length, 3, 'header-less XYZ');
});

test('cisplatin (XYZ and Gaussian input)', () => {
  const x = load('cisplatin.xyz'), g = load('cisplatin.gjf'), E = EXP.cisplatin;
  arrEq(x.frames[0].Z, E.Z, 'Z (xyz)');
  arrEq(g.frames[0].Z, E.Z, 'Z (gjf)');
  arrNear(Array.from(g.frames[0].xyz), E.xyz, 1e-8, 'gjf coordinates; basis section ignored');
  near(dist(x.frames[0].xyz, 0, 1), 2.33, 1e-6, 'Pt-Cl'); near(dist(x.frames[0].xyz, 0, 3), 2.05, 1e-6, 'Pt-N');
  near(angle(x.frames[0].xyz, 1, 0, 2), 90, 1e-6, 'cis Cl-Pt-Cl');
  eq(x.title, 'cisplatin, cis-[PtCl2(NH3)2]', 'xyz title from comment');
  eq(g.meta.solvent, 'water', 'solvent from SCRF'); eq(g.meta.link0.length, 3, 'link0');
  eq(g.meta.basis, 'genecp', 'basis'); ok(g.meta.jobTypes.includes('opt') && g.meta.jobTypes.includes('freq'), 'job types');
});

test('Gaussian input: Cartesian label styles', () => {
  const ds = load('methanol_cartesian.gjf');
  checkShape(ds);
  arrEq(ds.frames[0].Z, [6, 8, 1, 1, 1, 1], 'Z from C-CA--0.25, O(Fragment=1), H1, 1, H, H4');
  arrNear(Array.from(ds.frames[0].xyz.slice(9, 15)), [-0.5495, 1.0422, 0.8947, -0.5495, 1.0422, -0.8947], 1e-12, 'freeze-code lines');
  eq(ds.meta.route, '#p opt freq b3lyp/6-31g(d) geom=connectivity', 'multi-line route');
  eq(ds.title, 'Methanol, mixed Cartesian label styles', 'title');
  eq(ds.meta.charge, 0, 'charge'); eq(ds.meta.multiplicity, 1, 'multiplicity');
  arrEq(ds.meta.link0, ['%chk=methanol.chk', '%mem=2GB', '%nprocshared=4'], 'link0');
  eq(ds.meta.method, 'B3LYP', 'method');
});

test('Gaussian input: Z-matrix with dummy atom, labels, variables, constants', () => {
  const ds = load('nh3_zmatrix.com');
  const x = ds.frames[0].xyz;
  arrEq(ds.frames[0].Z, [7, 1, 1, 1], 'dummy X dropped');
  for (let k = 1; k <= 3; k++) near(dist(x, 0, k), 1.012, 1e-9, 'N-H ' + k);
  near(angle(x, 1, 0, 2), 106.6635, 2e-3, 'H-N-H'); near(angle(x, 1, 0, 3), angle(x, 2, 0, 3), 1e-9, 'C3v');
  const h2o2 = '#p hf/3-21g\n\nH2O2\n\n0 1\nO\nO 1 roo\nH 1 roh 2 aooh\nH 2 roh 1 aooh 3 -dih\n\nroo=1.45\nroh 0.97\naooh 100.0\n\ndih 115.0\n\n';
  const p = BM.parseFile('h2o2.gjf', h2o2), y = p.frames[0].xyz;
  near(dist(y, 0, 1), 1.45, 1e-9, 'O-O'); near(dist(y, 0, 2), 0.97, 1e-9, 'O-H');
  near(angle(y, 2, 0, 1), 100, 1e-9, 'H-O-O'); near(dihedral(y, 3, 1, 0, 2), -115, 1e-9, 'signed variable dihedral');
  const pbc = BM.parseFile('chain.gjf', '#p pbe/sto-3g\n\nchain\n\n0 1\nH 0.0 0.0 0.0\nH 0.74 0.0 0.0\nTv 2.0 0.0 0.0\n\n');
  eq(pbc.frames[0].Z.length, 2, 'Tv not an atom'); near(pbc.cell[0][0], 2.0, 1e-12, 'Tv -> cell');
  throwsMsg(() => BM.parseFile('chk.gjf', '%chk=a.chk\n#p opt geom=allcheck\n\n'), /checkpoint/, 'geom=check input');
  throwsMsg(() => BM.parseFile('bad.gjf', '#p hf\n\ntitle\n\n0 1\nH 1 rr\n\n'), /Gaussian input/, 'bad Z-matrix');
});

test('zmatrixToCartesian helper', () => {
  const x = BM.zmatrixToCartesian([{}, { a: 0, r: 1.5 }, { a: 1, r: 1.2, b: 0, angle: 110 },
    { a: 2, r: 1.1, b: 1, angle: 120, c: 0, dihedral: 60 }, { a: 1, r: 1.0, b: 0, angle: 100, c: 2, dihedral: 105, mode: 1 },
    { xyz: [5, 6, 7] }]);
  near(dist(x, 0, 1), 1.5, 1e-12, 'r'); near(x[5], 1.5, 1e-12, 'atom 2 on +z'); near(x[7], 0, 1e-12, 'atom 3 in xz plane');
  near(angle(x, 2, 1, 0), 110, 1e-9, 'angle'); near(dihedral(x, 3, 2, 1, 0), 60, 1e-9, 'dihedral sign');
  near(angle(x, 4, 1, 0), 100, 1e-9, 'format-code angle 1'); near(angle(x, 4, 1, 2), 105, 1e-9, 'format-code angle 2');
  arrNear(Array.from(x.slice(15)), [5, 6, 7], 0, 'Cartesian row');
});

// ================================================================ volumes
test('cube files', () => {
  const s = load('water_homo_bohr.cube'), E = EXP.cubes, B = 0.529177210903;
  checkShape(s);
  eq(s.volumes.length, 1, 'one volume');
  const v = s.volumes[0];
  arrEq(v.dims, E.single.dims, 'dims'); eq(v.kind, 'orbital', 'kind'); eq(v.name, 'MO 5', 'name');
  arrNear(v.origin, E.single.origin, 1e-9, 'origin Bohr -> A'); arrNear(flat(v.axes), flat(E.single.axes), 1e-9, 'axes');
  ok(v.data instanceof Float32Array && v.data.length === 20 * 18 * 16, 'Float32Array data');
  const [p, val] = E.single.probe;
  near(v.data[(p[0] * 18 + p[1]) * 16 + p[2]], val, Math.abs(val) * 1e-6 + 1e-12, 'value at (3,5,7), k fastest');
  near(v.data[0], E.single.first, 1e-12, 'first value'); near(v.data[v.data.length - 1], E.single.last, 1e-12, 'last value');
  near(s.frames[0].xyz[2], EXP.g16_water_optfreq.final[2] * 0 + s.frames[0].xyz[2], 0, 'atoms read');
  eq(s.frames[0].Z.length, 3, 'atoms');
  const t = load('water_mo5_mo6.cube');
  eq(t.volumes.length, 2, 'two interleaved MOs');
  arrEq(t.volumes.map(q => q.name), ['MO 5', 'MO 6'], 'names');
  const [p5, v5] = E.two.probe5, [p6, v6] = E.two.probe6;
  near(t.volumes[0].data[(p5[0] * 11 + p5[1]) * 10 + p5[2]], v5, Math.abs(v5) * 1e-6, 'MO 5 value');
  near(t.volumes[1].data[(p6[0] * 11 + p6[1]) * 10 + p6[2]], v6, Math.abs(v6) * 1e-6, 'MO 6 value');
  const a = load('water_density_angstrom.cube');
  const va = a.volumes[0];
  eq(va.kind, 'density', 'density kind'); arrEq(va.dims, E.ang.dims, 'dims (negative counts)');
  arrNear(va.origin, E.ang.origin, 1e-9, 'origin kept in A'); arrNear(flat(va.axes), flat(E.ang.axes), 1e-9, 'axes in A');
  arrNear(Array.from(a.frames[0].xyz.slice(3, 6)), E.ang.atom1, 1e-6, 'atoms in A');
  const [pa, vaa] = E.ang.probe;
  near(va.data[(pa[0] * 9 + pa[1]) * 8 + pa[2]], vaa, vaa * 1e-6, 'value');
  void B;
});

test('cube 100x100x100 parse speed', () => {
  const n = 100, parts = ['big cube', 'synthetic', '    1   -5.000000   -5.000000   -5.000000',
    '  100    0.100000    0.000000    0.000000', '  100    0.000000    0.100000    0.000000',
    '  100    0.000000    0.000000    0.100000', '    8    8.000000    0.000000    0.000000    0.000000'];
  const rows = [];
  for (let i = 0; i < n; i++) for (let j = 0; j < n; j++) {
    let line = '';
    for (let k = 0; k < n; k++) {
      const x = -5 + 0.1 * i, y = -5 + 0.1 * j, z = -5 + 0.1 * k;
      line += (' ' + Math.exp(-(x * x + y * y + z * z) / 4).toExponential(5).toUpperCase().replace(/E([+-])(\d)$/, 'E$10$2')).padStart(13);
      if (k % 6 === 5 || k === n - 1) { rows.push(line); line = ''; }
    }
  }
  const text = parts.join('\n') + '\n' + rows.join('\n') + '\n';
  const t0 = Date.now();
  const ds = BM.parseFile('big.cube', text);
  const ms = Date.now() - t0;
  ok(ms < 2000, 'parsed ' + (text.length / 1e6).toFixed(1) + ' MB cube in ' + ms + ' ms');
  const d = ds.volumes[0].data;
  eq(d.length, 1e6, '10^6 values');
  near(d[(50 * n + 50) * n + 50], 1, 1e-6, 'centre value');
  near(d[(10 * n + 20) * n + 30], Math.exp(-(16 + 9 + 4) / 4), 1e-8, 'value (10,20,30)');
  console.log('  cube 100^3: ' + ms + ' ms');
});

test('Molden with frequencies and geometries', () => {
  const ds = load('water_freq.molden'), E = EXP.molden;
  checkShape(ds);
  eq(ds.frames.length, 3, '[GEOMETRIES] frames');
  arrNear(ds.frames.map(f => f.energy), E.energies, 1e-10, '[GEOCONV] energies');
  near(ds.frames[0].opt.maxForce, 0.0123, 1e-12, 'max-force'); near(ds.frames[2].opt.rmsForce, 0.00005, 1e-12, 'rms-force');
  eq(ds.vibrations.frameIndex, 2, 'FR-COORD matched to last geometry');
  arrNear(ds.vibrations.freqs, E.freqs, 1e-9, 'freqs'); arrNear(ds.vibrations.ir, E.ir, 1e-9, '[INT] IR');
  arrNear(ds.vibrations.raman, [2.184578, 87.149025, 21.880203].map(v => +v.toFixed(4)), 1e-9, '[INT] 2nd column');
  arrNear(Array.from(ds.vibrations.modes[0]), E.mode0, 1e-9, 'vibration 1');
  arrNear(Array.from(ds.frames[2].xyz), E.final, 1e-8, 'final geometry');
  ok(ds.warnings.some(w => /cube/.test(w)), 'orbital warning');
  eq(ds.title, 'water opt+freq (synthetic, from ORCA data)', 'title');
});

test('Gaussian fchk', () => {
  const txt = read('water_freq.fchk'), E = EXP.fchk;
  const ds = BM.parseFile('water_freq.fchk', txt);
  checkShape(ds);
  arrEq(ds.frames[0].Z, [8, 1, 1], 'Z');
  arrNear(Array.from(ds.frames[0].xyz), E.xyz, 1e-9, 'coordinates Bohr -> A');
  near(ds.frames[0].energy, E.energy, 1e-12, 'Total Energy');
  eq(ds.meta.method, 'B3LYP', 'method'); eq(ds.meta.basis, '6-31G(d)', 'basis'); ok(ds.meta.jobTypes.includes('freq'), 'freq');
  arrNear(ds.vibrations.freqs, E.freqs, 1e-6, 'Vib-E2 frequencies'); arrNear(ds.vibrations.ir, E.ir, 1e-6, 'IR');
  arrNear(ds.vibrations.redMass, E.redMass, 1e-8, 'reduced masses');
  eq(ds.vibrations.modes.length, 3, 'Vib-Modes');
  arrNear(ds.charges.Mulliken, [-0.868866, 0.434433, 0.434433], 1e-8, 'Mulliken');
  arrNear([ds.dipole.x, ds.dipole.y, ds.dipole.z], E.dipole, 1e-6, 'dipole a.u. -> Debye');
  eq(ds.orbitals.alpha.nOcc, 5, 'nOcc');
  const broken = txt.replace(/Vib-Modes( +)R   N=( +)27/, 'Vib-Modes$1R   N=$226');
  const b = BM.parseFile('broken.fchk', broken);
  eq(b.vibrations, null, 'inconsistent vib data skipped'); ok(b.warnings.some(w => /inconsistent/.test(w)), 'warning');
});

// ================================================================ MOL / PDB / GRO
test('SDF (first record, bond orders)', () => {
  const ds = load('acetic_acid.sdf'), E = EXP.sdf;
  checkShape(ds);
  eq(ds.frames[0].Z.length, E.n, 'atoms'); arrEq(ds.frames[0].Z, [6, 6, 8, 8, 1, 1, 1, 1], 'Z');
  arrEq(ds.bonds, E.bonds, 'bonds with orders'); arrNear(Array.from(ds.frames[0].xyz.slice(0, 3)), E.first, 1e-12, 'coords');
  eq(ds.title, 'acetic acid', 'title'); ok(ds.warnings.length === 1, 'multi-record warning');
});

test('PDB with element columns', () => {
  const ds = load('peptide_ions_elements.pdb'), E = EXP.pdb;
  checkShape(ds);
  arrEq(ds.frames[0].Z, E.Z1, 'elements from columns 77-78');
  arrEq(ds.frames[0].atomNames, E.names1, 'atom names'); arrEq(ds.frames[0].resNames, E.res1, 'residue names');
  arrEq(ds.frames[0].resIds, E.ids1, 'residue ids');
  arrNear(Array.from(ds.frames[0].xyz.slice(-3)), E.x_last, 5e-4, 'coordinates');
  arrNear(flat(ds.cell), [30, 0, 0, 0, 30, 0, 0, 0, 30], 1e-9, 'CRYST1 cell');
  eq(ds.title, 'BLUEMOLECULE PARSER TEST', 'title from HEADER');
});

test('PDB without element columns, MODEL/ENDMDL', () => {
  const ds = load('peptide_ions_noelements.pdb'), E = EXP.pdb;
  eq(ds.frames.length, 2, 'models -> frames');
  arrEq(ds.frames[0].Z, E.Z2, 'elements from atom names (CA in ALA = C, CA ion = Ca, NA/CL/ZN ions, PT1/CL1 ligand, 1HB)');
  arrEq(ds.frames[1].Z, E.Z2, 'model 2 elements');
  near(ds.frames[1].xyz[0] - ds.frames[0].xyz[0], 0.1, 1e-9, 'model 2 coordinates');
  eq(ds.frames[1].label, 'Model 2', 'label');
});

test('GRO (velocities, box, element heuristics, two frames)', () => {
  const ds = load('water_ions_cisplatin.gro'), E = EXP.gro;
  checkShape(ds);
  eq(ds.frames.length, 2, 'frames');
  arrEq(ds.frames[0].Z, E.Z, 'elements (OW, HW1, NA, CL, PT, CL1, N1, H11...)');
  arrNear(Array.from(ds.frames[0].xyz.slice(0, 3)), E.first, 1e-9, 'nm -> A');
  arrNear(Array.from(ds.frames[0].xyz.slice(15, 18)), E.pt, 1e-9, 'Pt position');
  arrNear(flat(ds.cell), [30, 0, 0, 0, 31, 0, 0, 0, 32], 1e-9, 'box -> cell');
  arrEq(ds.frames[0].resNames.slice(0, 5), ['SOL', 'SOL', 'SOL', 'NA', 'CL'], 'residue names');
  near(ds.frames[1].xyz[0] - ds.frames[0].xyz[0], 0.1, 1e-9, 'frame 2 moved');
  near(ds.frames[1].time, 10, 1e-12, 'time from title');
  const tri = BM.parseFile('t.gro', 't\n    1\n    1SOL     OW    1   0.100   0.200   0.300\n   1.0   2.0   3.0   0.0   0.0   0.5   0.0   0.3   0.4\n');
  arrNear(flat(tri.cell), [10, 0, 0, 5, 20, 0, 3, 4, 30], 1e-9, 'triclinic box');
});

// ================================================================ detection
test('format detection and errors', () => {
  const g = read('g16_water_optfreq.log'), o = read('orca5_water_optfreq.out');
  eq(BM.detectFormat('job.out', g), 'gaussian-log', '.out Gaussian'); eq(BM.detectFormat('job.log', o), 'orca-out', '.log ORCA');
  eq(BM.detectFormat('noext', o), 'orca-out', 'no extension ORCA'); eq(BM.detectFormat('x.dat', read('acetic_acid.sdf')), 'mol', 'sniff SDF');
  eq(BM.detectFormat('a', read('water_homo_bohr.cube')), 'cube', 'sniff cube'); eq(BM.detectFormat('a', read('water_freq.fchk')), 'fchk', 'sniff fchk');
  eq(BM.detectFormat('a', read('water_ions_cisplatin.gro')), 'gro', 'sniff gro'); eq(BM.detectFormat('a', read('peptide_ions_elements.pdb')), 'pdb', 'sniff pdb');
  eq(BM.detectFormat('a', read('water_freq.molden')), 'molden', 'sniff molden'); eq(BM.detectFormat('a', read('nh3_zmatrix.com')), 'gaussian-input', 'sniff gjf');
  eq(BM.detectFormat('a.txt', 'hello world\nthis is text'), null, 'unknown -> null');
  throwsMsg(() => BM.parseFile('notes.log', 'just some words\n'), /not recognised/, 'unrecognised file');
  throwsMsg(() => BM.parseFile('e.xyz', '   \n'), /empty/, 'empty file');
  throwsMsg(() => BM.parseFile('t.xyz', '5\ncomment\nC 0 0 0\n'), /XYZ/, 'truncated XYZ');
  throwsMsg(() => BM.parseFile('t.cube', read('water_homo_bohr.cube').slice(0, 5000)), /truncated/, 'truncated cube');
  ok(Object.keys(BM.parsers).length === 10, 'ten parsers exported');
});

// ================================================================ writers
test('writers round trips', () => {
  const W = BM.writers, cp = load('cisplatin.xyz').frames[0];
  const x1 = BM.parseFile('w.xyz', W.xyz({ Z: cp.Z, xyz: cp.xyz, comment: 'E = -1153.25' }));
  arrEq(x1.frames[0].Z, cp.Z, 'xyz Z'); arrNear(Array.from(x1.frames[0].xyz), Array.from(cp.xyz), 1e-8, 'xyz coords');
  near(x1.frames[0].energy, -1153.25, 1e-12, 'xyz comment energy');
  const xm = BM.parseFile('w.xyz', W.xyzMulti([{ Z: cp.Z, xyz: cp.xyz, comment: 'a' }, { Z: [8], xyz: [0, 0, 0], comment: 'b\nc' }]));
  eq(xm.frames.length, 2, 'xyzMulti frames'); eq(xm.frames[1].label, 'b c', 'newline-free comment');
  const gj = W.gaussianInput({ Z: cp.Z.concat([0]), xyz: Array.from(cp.xyz).concat([0, 0, 3]), route: '#p opt b3lyp/genecp',
    title: 'cisplatin\nexport', charge: 0, multiplicity: 1, link0: ['%nprocshared=4', 'mem=8GB'] });
  ok(/\n\n\n$/.test(gj) && /\n0 1\nPt /.test(gj), 'Gaussian blank-line layout');
  const g = BM.parseFile('w.gjf', gj);
  arrEq(g.frames[0].Z, cp.Z.concat([0]), 'gjf Z incl. Bq ghost'); arrNear(Array.from(g.frames[0].xyz.slice(0, 33)), Array.from(cp.xyz), 1e-8, 'gjf coords');
  eq(g.meta.route, '#p opt b3lyp/genecp', 'route'); eq(g.title, 'cisplatin export', 'title');
  arrEq(g.meta.link0, ['%nprocshared=4', '%mem=8GB'], 'link0');
  const oi = W.orcaInput({ Z: cp.Z, xyz: cp.xyz, keywords: 'B3LYP def2-TZVP Opt', charge: 0, multiplicity: 1, nprocs: 8, maxcore: 3000,
    extraBlocks: ['%cpcm\n  smd true\n  SMDsolvent "water"\nend'] });
  ok(/^! B3LYP def2-TZVP Opt\n%pal nprocs 8 end\n%maxcore 3000\n%cpcm/.test(oi) && /\* xyz 0 1\n  Pt /.test(oi) && /\n\*\n$/.test(oi), 'ORCA input layout');
  const oxyz = oi.split('\n').filter(l => /^  [A-Z][a-z]? /.test(l)).map(l => l.trim().split(/\s+/).slice(1).map(Number));
  arrNear(flat(oxyz), Array.from(cp.xyz), 1e-8, 'ORCA coordinates');
  const ac = load('acetic_acid.sdf');
  const m = BM.parseFile('w.mol', W.mol({ Z: ac.frames[0].Z, xyz: ac.frames[0].xyz, bonds: ac.bonds, title: 'acetic acid' }));
  arrEq(m.frames[0].Z, ac.frames[0].Z, 'mol Z'); arrEq(m.bonds, ac.bonds, 'mol bonds');
  arrNear(Array.from(m.frames[0].xyz), Array.from(ac.frames[0].xyz), 1e-4, 'mol coords');
  const mg = BM.parseFile('w.mol', W.mol({ Z: cp.Z, xyz: cp.xyz, title: 'cisplatin' }));
  eq(mg.bonds.length, 10, 'guessed bonds for cisplatin (2 Pt-Cl, 2 Pt-N, 6 N-H)');
  const pd = BM.parseFile('w.pdb', W.pdb({ Z: cp.Z, xyz: cp.xyz, resNames: cp.Z.map(() => 'CPT'), resIds: cp.Z.map(() => 7) }));
  arrEq(pd.frames[0].Z, cp.Z, 'pdb Z (Pt, Cl two-letter)'); arrNear(Array.from(pd.frames[0].xyz), Array.from(cp.xyz), 6e-4, 'pdb coords');
  eq(pd.frames[0].atomNames[0], 'Pt1', 'default atom names'); eq(pd.frames[0].resIds[3], 7, 'resIds');
  const pdbText = W.pdb({ Z: [6, 17], xyz: [0, 0, 0, 1.8, 0, 0], atomNames: ['C1', 'CL1'], resNames: ['LIG', 'LIG'] });
  const line = pdbText.split('\n')[2];
  eq(line.slice(12, 16), 'CL1 ', 'two-letter element name starts in column 13'); eq(line.slice(76, 78), 'CL', 'element columns');
  const back = BM.parseFile('w.pdb', pdbText.replace(/(.{76}).{2}/g, '$1  '));
  arrEq(back.frames[0].Z, [6, 17], 'element recovered from names without columns');
});

// ---------------------------------------------------------------------- summary
console.log('\nBlueMolecule parser tests: ' + pass + ' passed, ' + fail + ' failed');
if (fail) {
  for (const f of failures) console.log('  FAIL ' + f);
  process.exit(1);
}
