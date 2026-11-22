#!/usr/bin/env node
/* BlueMolecule tests for the 2.3 readers (src/readers-qc2.js and the fchk additions in src/parsers.js):
 *   Molpro outputs and logs, Dalton outputs, Turbomole program outputs (dscf, escf, aoforce, ricc2), GAMESS punch
 *   files (with point-group expansion), and fchk excited states and optimisation history.
 * The fixtures come from tools/make_fixtures_qc2.py, in the layouts of real outputs (the cclib test files); their
 * values are in tests/fixtures/qc/expected_qc2.json.
 * Run:  node tests/test_readers_qc2.js */
'use strict';
const fs = require('fs');
const path = require('path');
const ROOT = path.join(__dirname, '..');
for (const m of ['elements', 'bigfile', 'parsers', 'condense', 'readers-struct', 'writers', 'readers-qc', 'readers-qc2', 'geom', 'chem', 'analysis']) {
  require(path.join(ROOT, 'src', m + '.js'));
}
const BM = globalThis.BM;
const QC = path.join(__dirname, 'fixtures', 'qc');
const read = f => fs.readFileSync(path.join(QC, f), 'utf8');
const load = f => BM.parseFile(f, read(f));
const X = JSON.parse(read('expected_qc2.json'));

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
const xyzOf = f => Array.from(f.xyz);
const conv = ds => BM.analysis.convergenceSeries(ds);

// ------------------------------------------------------------------ Molpro
test('Molpro optimisation and frequencies', () => {
  const t = read('molpro_water_optfreq.out'), ds = BM.parseFile('water.out', t), E = X.molpro;
  eq(BM.detectFormat('water.out', t), 'molpro-out', 'detected by its banner');
  eq(ds.program, 'Molpro 2018.2', 'program and version');
  eq(ds.title, 'water optimisation and frequencies', 'title from the input echo');
  eq(ds.meta.method + '/' + ds.meta.basis, 'B3LYP/cc-pVDZ', 'method (functional of the ks command) and basis (as written in the summary)');
  eq(ds.meta.charge + '/' + ds.meta.multiplicity, '0/1', 'charge and multiplicity from the nuclear charge and the electrons');
  eq(ds.frames.length, 2, 'the first and the last geometry (the output has no others)');
  arrNear(xyzOf(ds.frames[0]), E.first, 1e-8, 'first geometry (Bohr in the output)');
  arrNear(xyzOf(ds.frames[1]), E.last, 1e-9, 'optimised geometry (the xyz block)');
  near(ds.frames[0].energy, E.E[0], 1e-11, 'energy of the first geometry');
  near(ds.frames[1].energy, E.E[2], 1e-11, 'energy of the optimised geometry');
  ok(ds.frames[1].converged, 'the optimised geometry is converged');
  eq(ds.optSteps.length, 3, 'the table: two cycles and the final geometry');
  arrNear(ds.optSteps.map(s => s.energy), E.cycles, 1e-8, 'cycle energies');
  eq(ds.optSteps[0].opt.maxForceThr, 3e-4, 'gradient threshold from THRESHOLDS (OPTGRAD)');
  near(ds.optSteps[0].opt.stepNorm, 1.9e-2, 1e-12, 'step length (STEPLEN)');
  const c = conv(ds);
  eq(c.source, 'optSteps', 'the convergence plot uses the cycles');
  eq(c.converged, 2, 'converged at the last cycle (Molpro: gradient, and energy or step)');
  arrNear(ds.vibrations.freqs, E.freqs, 1e-9, 'frequencies');
  arrNear(ds.vibrations.ir, E.ir, 1e-9, 'IR intensities');
  eq(ds.vibrations.sym.join(' '), 'A1 A1 B2', 'symmetries');
  eq(ds.vibrations.frameIndex, 1, 'vibrations of the last geometry');
  ok(ds.vibrations.modes.every(m => m.length === 9 && m.some(v => v !== 0)), 'displacements of every mode');
  near(ds.thermo.zpe, E.zpe, 1e-9, 'zero-point energy');
  eq(ds.orbitals.alpha.nOcc, 5, 'occupied orbitals');
  arrNear(ds.orbitals.alpha.energies, E.occ.slice().sort((a, b) => a - b), 1e-5, 'orbital energies (sorted over symmetries)');
  arrNear(ds.charges.Mulliken, E.mulliken, 1e-9, 'Mulliken charges ("- 0.44900")');
  near(ds.dipole.total, 0.8256 * 2.541746, 1e-6, 'dipole (a.u. to Debye)');
  ok(ds.meta.normalTermination === true, 'normal termination');
  ok(ds.jobs && ds.jobs[0].kinds.join() === 'opt,freq', 'job kinds');
});
test('Molpro 2006: the table after END OF GEOMETRY OPTIMIZATION and no final geometry', () => {
  let t = read('molpro_water_optfreq.out');
  const a = t.indexOf(' Current geometry (xyz format'), b = t.indexOf(' Normal Modes\n');
  t = t.slice(0, a) + t.slice(b);   // no xyz block: the final energy belongs to a geometry the output does not print
  const ds = BM.parseFile('water.out', t);
  eq(ds.frames.length, 1, 'one geometry');
  near(ds.frames[0].energy, X.molpro.E[0], 1e-11, 'the first geometry keeps its own energy');
  ok(ds.warnings.some(w => /Z-matrix variables/.test(w)), 'a warning says the structure is the first one');
  eq(ds.optSteps.length, 3, 'the cycles are still read');
});
test('Molpro logs', () => {
  const atoms = (z) => [' ATOMIC COORDINATES', '', ' NR  ATOM    CHARGE       X              Y              Z', '',
    '   1  O       8.00    0.000000000    0.000000000    ' + z, '   2  H1      1.00    0.000000000    1.430000000   -0.880000000',
    '   3  H2      1.00    0.000000000   -1.430000000   -0.880000000', ''];
  const log = ['1PROGRAM * ALASKA (Gradient of the energy)     Author: R. Lindh', '', ' Optimization point 1', '', ' Geometry written to block  2 of record 700', '',
    ' PROGRAM * SEWARD (Integral evaluation for generally contracted gaussian basis sets)     Author: Roland Lindh, 1990', '']
    .concat(atoms('0.220000000'), [' NUMBER OF ELECTRONS:       5+    5-    SPACE SYMMETRY=1', ' !RKS STATE  1.1 Energy              -76.408800000000', '',
      ' Optimization point 2', ''], atoms('0.225000000'), [' !RKS STATE  1.1 Energy              -76.408950000000', '']).join('\n');
  eq(BM.detectFormat('water.log', log), 'molpro-out', 'an optimisation log (no banner) is a Molpro file');
  const ds = BM.parseFile('water.log', log);
  eq(ds.frames.length, 2, 'one geometry per optimisation point');
  arrNear(ds.frames.map(f => f.energy), [-76.4088, -76.40895], 1e-9, 'their energies');
  ok(ds.meta.jobTypes.indexOf('opt') >= 0, 'an optimisation');
  eq(ds.meta.method, 'DFT', 'an RKS energy without the input: DFT');
  const fq = [' Cartesian Coordinates are in [bohr]', ' Creating new Variables:'].concat(atoms('0.220000000'), [' !RKS STATE  1.1 Energy  -76.408953000000'],
    atoms('0.230000000'), [' !RKS STATE  1.1 Energy  -76.408900000000'], atoms('0.210000000'), [' !RKS STATE  1.1 Energy  -76.408901000000']).join('\n');
  const f2 = BM.parseFile('water.log', fq);
  eq(f2.frames.length, 1, 'a numerical-frequency log: the reference geometry only');
  near(f2.frames[0].energy, -76.408953, 1e-9, 'with its energy');
  ok(f2.warnings.some(w => /displaced/.test(w)), 'and a warning about the displaced geometries');
});

// ------------------------------------------------------------------ Dalton
test('Dalton optimisation and frequencies', () => {
  const t = read('dalton_water_optfreq.out'), ds = BM.parseFile('water.out', t), E = X.dalton;
  eq(BM.detectFormat('water.out', t), 'dalton-out', 'detected by its banner');
  eq(ds.program, 'Dalton 2015.0', 'program and version');
  eq(ds.meta.method + '/' + ds.meta.basis, 'B3LYP/cc-pVDZ', 'method (.DFT B3LYP) and basis');
  eq(ds.meta.charge + '/' + ds.meta.multiplicity, '0/1', 'charge and multiplicity');
  eq(ds.frames.length, 3, 'three geometries (repeated ones merged)');
  E.frames.forEach((f, k) => { arrNear(xyzOf(ds.frames[k]), f.xyz, 2e-9, 'geometry ' + (k + 1)); near(ds.frames[k].energy, f.energy, 1e-11, 'energy ' + (k + 1)); });
  ok(ds.frames[2].converged && !ds.frames[1].converged, 'the last geometry is converged');
  near(ds.frames[1].opt.gradNorm, 0.000458, 1e-12, 'gradient norm of the @ Iter table');
  near(ds.frames[1].opt.stepNorm, 0.00061, 1e-12, 'step length of the @ Iter table');
  eq(conv(ds).converged, 2, 'convergence plot: converged at the last geometry');
  arrNear(ds.vibrations.freqs, E.freqs, 1e-9, 'frequencies (Dalton lists the highest first)');
  arrNear(ds.vibrations.ir, E.ir, 1e-9, 'IR intensities');
  eq(ds.vibrations.sym.join(' '), 'A1 A1 B2', 'symmetries');
  near(ds.thermo.zpe, E.zpe, 1e-9, 'zero-point energy (Total Molecular Energy)');
  near(ds.dipole.total, E.dipole, 1e-5, 'dipole');
  eq(ds.orbitals.alpha.nOcc, 5, 'occupied orbitals from the electron count');
  eq(ds.orbitals.alpha.energies.length, 8, 'orbital energies');
  ok(ds.meta.normalTermination === true, 'normal termination');
});
test('Dalton response (excited states)', () => {
  const ds = load('dalton_water_tddft.out'), E = X.dalton_td;
  eq(ds.frames.length, 1, 'the molecule of the .mol echo (the output prints no coordinates)');
  arrNear(xyzOf(ds.frames[0]), E.xyz, 1e-7, 'its coordinates, as far as the 50 columns of the echo go');
  ok(ds.warnings.some(w => /50 columns/.test(w)), 'with a warning about the cut echo');
  eq(ds.excitedStates.length, 4, 'four states');
  arrNear(ds.excitedStates.map(s => s.energyEV), E.states.map(s => s[0]), 1e-6, 'sorted by energy over the symmetries');
  arrNear(ds.excitedStates.map(s => s.f), E.states.map(s => s[1]), 1e-8, 'oscillator strengths: the sum of the components, 0 when forbidden');
  eq(ds.excitedStates[0].label, 'B1, state 1', 'symmetry label');
  ok(ds.meta.jobTypes.indexOf('td') >= 0, 'a TD job');
});

// ------------------------------------------------------------------ Turbomole
test('Turbomole dscf', () => {
  const t = read('turbomole_water_dscf.out'), ds = BM.parseFile('dscf.out', t), E = X.turbomole;
  eq(BM.detectFormat('dscf.out', t), 'turbomole-out', 'detected by its program line');
  eq(ds.program, 'Turbomole 7.4.1 (dscf)', 'program, version and module');
  arrNear(xyzOf(ds.frames[0]), E.xyz, 1e-8, 'coordinates (Bohr)');
  near(ds.frames[0].energy, E.energy, 1e-10, 'total energy');
  eq(ds.meta.method + '/' + ds.meta.basis, 'B3LYP/def2-SVP', 'functional and basis');
  eq(ds.meta.charge + '/' + ds.meta.multiplicity, '0/1', 'charge from the electrostatic moments');
  near(ds.dipole.total, E.dipole, 1e-5, 'dipole');
  ok(ds.meta.normalTermination === true, 'ended normally');
});
test('Turbomole escf', () => {
  const ds = load('turbomole_water_escf.out'), E = X.turbomole_escf;
  eq(ds.excitedStates.length, 4, 'four states');
  arrNear(ds.excitedStates.map(s => s.energyEV), E.states.map(s => s[0]), 1e-9, 'energies (eV)');
  arrNear(ds.excitedStates.map(s => s.f), E.states.map(s => s[1]), 1e-12, 'oscillator strengths (length)');
  arrNear(ds.excitedStates.map(s => s.R), E.R, 1e-12, 'rotatory strengths (length, 10^-40 cgs) for the ECD spectrum');
  eq(ds.excitedStates[0].label, 'Singlet b1', 'multiplicity and symmetry');
  const tr = ds.excitedStates[0].transitions[0];
  eq(tr.from + ' -> ' + tr.to, '5a1 -> 6a1', 'dominant contribution');
  near(tr.coeff * tr.coeff, E.weight[0], 1e-9, 'its weight (|coeff|^2)');
});
test('Turbomole aoforce', () => {
  const ds = load('turbomole_water_aoforce.out'), E = X.turbomole_force;
  arrNear(ds.vibrations.freqs, E.freqs, 1e-9, 'frequencies (the six zero modes dropped)');
  arrNear(ds.vibrations.ir, E.ir, 1e-9, 'IR intensities');
  eq(ds.vibrations.sym.join(' '), 'a1 a1 b2', 'symmetries');
  ok(ds.vibrations.modes.every(m => m.length === 9 && Math.abs(Math.hypot(...m) - 1) < 1e-4), 'displacements (normalised to unity)');
  near(ds.thermo.zpe, E.zpe, 1e-9, 'zero-point energy');
});
test('Turbomole ricc2', () => {
  const ds = load('turbomole_water_ricc2.out'), E = X.turbomole_ricc2;
  near(ds.frames[0].energies.MP2, E.mp2, 1e-10, 'MP2 energy');
  near(ds.frames[0].energy, E.cc2, 1e-10, 'the frame energy is the CC2 one (not the MP2 printed first)');
  eq(ds.meta.method, 'CC2', 'method');
  arrNear(ds.excitedStates.map(s => s.energyEV), E.states.map(s => s[0]), 1e-5, 'the CC2 states (the last table), not the CCS ones nor the summary');
  arrNear(ds.excitedStates.map(s => s.f), E.states.map(s => s[1]), 1e-8, 'oscillator strengths from the transition blocks');
  eq(ds.excitedStates[0].label, 'Singlet a (CC2)', 'label with the model');
});
test('Turbomole outputs without coordinates', () => {
  let msg = '';
  try { BM.parseFile('define.out', ' define (host) : TURBOMOLE rev. V7.4.1 compiled\n\n nothing here\n'); } catch (e) { msg = e.message; }
  ok(/prints no atomic coordinates.*coord file/.test(msg), 'explained: ' + msg);
});

// ------------------------------------------------------------------ GAMESS punch files
test('GAMESS punch: optimisation in C2v (symmetry-unique atoms)', () => {
  const t = read('gamess_water_opt.dat'), ds = BM.parseFile('water.dat', t), E = X.gamess_opt;
  eq(BM.detectFormat('water.dat', t), 'gamess-punch', 'a .dat file with $DATA and punched results');
  eq(ds.frames.length, 3, 'one geometry per NSERCH (the final one merged)');
  E.frames.forEach((f, k) => { arrNear(xyzOf(ds.frames[k]), f.xyz, 1e-9, 'geometry ' + (k + 1) + ': O and one H expanded to water'); near(ds.frames[k].energy, f.energy, 1e-10, 'energy ' + (k + 1)); });
  ok(ds.frames[2].converged, 'the RESULTS FROM SUCCESSFUL ... GEOMETRY SEARCH geometry is converged');
  eq(ds.meta.method, 'B3LYP', 'method from E(R-B3LYP)');
  ok(ds.warnings.some(w => /generated from its symmetry-unique/.test(w)), 'a warning about the atom order');
});
test('GAMESS punch: C1 Hessian with normal modes and the AIMPAC section', () => {
  const t = read('gamess_water_hess.dat'), ds = BM.parseFile('water.dat', t), E = X.gamess_hess;
  arrNear(xyzOf(ds.frames[0]), E.xyz, 1e-9, 'geometry');
  near(ds.frames[0].energy, E.energy, 1e-10, 'energy');
  arrNear(ds.vibrations.freqs, E.freqs, 1e-9, 'frequencies (translations and rotations dropped)');
  ok(ds.vibrations.modes.every(m => m.length === 9), 'displacements');
  eq(ds.orbitals.alpha.nOcc, 5, 'occupied orbitals (AIMPAC)');
  arrNear(ds.orbitals.alpha.energies, E.occ, 1e-8, 'their energies');
  near(ds.dipole.total, 2.098512, 1e-9, 'dipole (Debye)');
  eq(ds.warnings.length, 0, 'no warning');
  // GAMESS orders the modes by eigenvalue: a mode before the six of translation and rotation is imaginary
  const lines = t.split('\n'), m0 = lines.findIndex(l => /^MODE\s+1\s/.test(l)), m7 = lines.findIndex(l => /^MODE\s+7\s/.test(l));
  const ts = lines.slice(0, m0).concat(lines.slice(m7, m7 + 4).map(l => l.replace(/MODE\s+7\s+FREQUENCY=\s*\S+/, 'MODE     1   FREQUENCY= 812.34000')), lines.slice(m0)).join('\n');
  const d2 = BM.parseFile('ts.dat', ts);
  eq(d2.vibrations.freqs.length, 4, 'four vibrations with the extra mode');
  near(d2.vibrations.freqs[0], -812.34, 1e-9, 'the first mode is imaginary');
});
test('GAMESS punch: symmetric file without AIMPAC keeps no modes', () => {
  const t = read('gamess_water_hess.dat').replace('C1       0', 'CNV      2\n')
    .replace(/(H {10}1\.0[^\n]*\n {3}N31 {5}6\n {8}\n)H {10}1\.0[^\n]*\n {3}N31 {5}6\n {8}\n/, '$1')
    .replace(/----- TOP OF INPUT FILE FOR BADER[\s\S]*?----- END OF INPUT FILE FOR BADER'S AIMPAC PROGRAM -----/, '');
  const ds = BM.parseFile('w.dat', t);
  eq(ds.frames[0].Z.length, 3, 'water expanded from O and one H');
  eq(ds.vibrations, null, 'the normal modes follow an atom order the file does not give: not shown');
  ok(ds.warnings.some(w => /normal modes are not shown/.test(w)), 'explained');
});
test('GAMESS punch detection', () => {
  const psi4 = 'molecule {\n0 1\nO\nH 1 0.96\nH 1 0.96 2 104.5\n}\nset basis cc-pvdz\nenergy(\'scf\')\n';
  eq(BM.detectFormat('input.dat', psi4), 'psi4-input', 'a Psi4 input .dat is still a Psi4 input');
  const inp = ' $CONTRL SCFTYP=RHF RUNTYP=ENERGY $END\n $DATA\nwater\nC1\nO 8.0 0.0 0.0 0.0\n $END\n';
  eq(BM.detectFormat('water.inp', inp), 'gamess-input', 'an input with $CONTRL is still an input');
});
test('GAMESS input with a point group', () => {
  const inp = ' $CONTRL SCFTYP=RHF RUNTYP=ENERGY $END\n $DATA\nmethane\nTD\n\nC 6.0 0.0 0.0 0.0\nH 1.0 0.629 0.629 0.629\n $END\n';
  const ds = BM.parseFile('ch4.inp', inp);
  eq(ds.frames[0].Z.join(','), '6,1,1,1,1', 'TD: C and four H from one');
  const d = [1, 2, 3, 4].map(k => Math.hypot(...Array.from(ds.frames[0].xyz.slice(3 * k, 3 * k + 3))));
  arrNear(d, [1.0895, 1.0895, 1.0895, 1.0895].map(() => Math.sqrt(3) * 0.629), 1e-9, 'all C-H at the same distance');
  const benz = ' $DATA\nbenzene\nDNH 6\n\nC 6.0 1.39 0.0 0.0\nH 1.0 2.47 0.0 0.0\n $END\n $CONTRL SCFTYP=RHF $END\n';
  eq(BM.parseFile('bz.inp', benz).frames[0].Z.length, 12, 'DNH 6: benzene from one C and one H');
  const sf6 = ' $CONTRL SCFTYP=RHF $END\n $DATA\nSF6\nOH\n\nS 16.0 0.0 0.0 0.0\nF 9.0 1.56 0.0 0.0\n $END\n';
  eq(BM.parseFile('sf6.inp', sf6).frames[0].Z.length, 7, 'OH: SF6 from S and one F');
  const ops = BM.qcReaders._h.gamessOps;
  eq([['CS', 0], ['CI', 0], ['CN', 3], ['S2N', 2], ['CNH', 2], ['CNV', 3], ['DN', 3], ['DNH', 6], ['DND', 2], ['T', 0], ['TH', 0], ['TD', 0], ['O', 0], ['OH', 0]]
    .map(([g, n]) => ops(g, n).length).join(' '), '2 2 3 4 4 6 6 24 8 12 24 24 24 48', 'group orders');
});

// ------------------------------------------------------------------ fchk files
test('fchk: optimisation history', () => {
  const t = read('g16_water_opt.fchk'), ds = BM.parseFile('opt.fchk', t), E = X.fchk_opt;
  eq(ds.frames.length, 3, 'every geometry of "Opt point 1"');
  arrNear(ds.frames.map(f => f.energy), E.E, 1e-9, 'their energies (the last one at full precision)');
  arrNear(xyzOf(ds.frames[2]), E.last, 1e-9, 'the last geometry is the current one');
  ok(ds.frames[2].converged && !ds.frames[0].converged, 'the last geometry is converged');
  eq(ds.frames.map(f => f.label).join(', '), 'Step 1, Step 2, Step 3', 'labels');
  eq(ds.meta.route, '#p b3lyp/6-31g(d) opt', 'route (a text array)');
  eq(ds.meta.basis, '6-31G(d)', 'basis from line 2 (Gaussian 16 gives the method 60 columns)');
  const sp = t.replace(/^FOpt/m, 'SP  ').replace(/#p b3lyp\/6-31g\(d\) opt/, '#p b3lyp/6-31g(d) geom=check');
  eq(BM.parseFile('sp.fchk', sp).frames.length, 1, 'a single point that read its geometry from an optimisation\'s checkpoint: no history');
  const moved = t.replace(/(Current cartesian coordinates\s+R\s+N=\s+9\n\s+)(\S+)/, (m0, a, b) => a + (parseFloat(b) + 0.01).toExponential(8).toUpperCase());
  eq(BM.parseFile('m.fchk', moved).frames.length, 1, 'a history that does not end at the current geometry is not used');
});
test('fchk: excited states', () => {
  const ds = load('g16_water_td.fchk'), E = X.fchk_td;
  eq(ds.excitedStates.length, 4, 'four states (ETran scalars)');
  arrNear(ds.excitedStates.map(s => s.energyEV), E.states.map(s => s[0]), 2e-6, 'energies from the state energies minus the SCF energy');
  arrNear(ds.excitedStates.map(s => s.f), E.states.map(s => s[1]), 2e-6, 'oscillator strengths from the transition dipoles');
  eq(ds.excitedStates.map(s => s.label).join(' '), E.states.map(s => s[2]).join(' '), 'singlets and triplets (ETran spin)');
  ok(ds.meta.jobTypes.indexOf('td') >= 0, 'a TD job');
  near(ds.frames[0].energy, ds.frames[0].energies.SCF, 0, 'a TD single point shows the ground-state energy');
  const q = load('qchem_water_cis.fchk');
  arrNear(q.excitedStates.map(s => s.energyEV), X.fchk_qchem.states.map(s => s[0]), 1e-6, 'Q-Chem: excitation energies (Eh)');
  arrNear(q.excitedStates.map(s => s.f), X.fchk_qchem.states.map(s => s[1]), 1e-9, 'Q-Chem: oscillator strengths');
});

test('fchk: Q-Chem coordinates of every step', () => {
  // Q-Chem writes "Current cartesian coordinates" again at each step of an optimisation or MD run
  const t = read('qchem_water_cis.fchk').replace(/^SP        /m, 'FOPT      ');
  const blk = t.slice(t.indexOf('Current cartesian coordinates'), t.indexOf('Number of Excited States'));
  const step = blk.replace(/(R   N=\s+9\n\s+)(\S+)/, (m0, a, b) => a + (parseFloat(b) + 0.05).toExponential(8).toUpperCase());
  const two = t.replace(blk, step + blk);
  const ds = BM.parseFile('q.fchk', two);
  eq(ds.frames.length, 2, 'two steps');
  eq(ds.frames.map(f => f.label).join(', '), 'Step 1, Step 2', 'labelled as steps');
  near(ds.frames[0].xyz[0] - ds.frames[1].xyz[0], 0.05 * 0.529177210903, 1e-9, 'the first step is the first coordinate set');
  eq(BM.parseFile('q.fchk', two.replace(/^FOPT      /m, 'FREQ      ')).frames.length, 1, 'a frequency job keeps its last (reference) geometry');
});

// ------------------------------------------------------------------ robustness: cut while running, damaged
for (const f of ['molpro_water_optfreq.out', 'dalton_water_optfreq.out', 'dalton_water_tddft.out', 'turbomole_water_dscf.out', 'turbomole_water_escf.out',
  'turbomole_water_aoforce.out', 'turbomole_water_ricc2.out', 'gamess_water_opt.dat', 'gamess_water_hess.dat', 'g16_water_opt.fchk', 'g16_water_td.fchk']) {
  test('cut ' + f, () => {
    const text = read(f), n = BM.parseFile(f, text).frames[0].Z.length;
    let bad = [], read2 = 0;
    for (let k = 1; k < 60; k++) {
      const c = Math.floor(text.length * k / 60);
      let ds;
      try { ds = BM.parseFile(f, text.slice(0, c)); } catch (e) { if (/is not a function|Cannot read prop|undefined|NaN/.test(e.message)) bad.push('@' + c + ' ' + e.message); continue; }
      read2++;
      ds.frames.forEach((fr, i) => { if (fr.Z.length !== n || fr.xyz.some(v => !Number.isFinite(v))) bad.push('@' + c + ' frame ' + (i + 1) + ': ' + fr.Z.length + ' atoms'); });
      if (ds.vibrations) ds.vibrations.modes.forEach((m, i) => { if (m.length !== 3 * n) bad.push('@' + c + ' mode ' + i); });
    }
    ok(!bad.length, f + ' cut at 59 places (' + read2 + ' read): ' + bad.slice(0, 3).join('; '));
  });
}

console.log(failures.length ? failures.map(f => '  FAIL ' + f).join('\n') + '\n' : '');
console.log('BlueMolecule 2.3 reader tests: ' + pass + ' passed, ' + fail + ' failed');
process.exit(fail ? 1 : 0);
