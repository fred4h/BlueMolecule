#!/usr/bin/env node
/* BlueMolecule tests for QM/MM and ONIOM outputs: Gaussian ONIOM (2 and 3 layers), ORCA QM/MM and QM/XTB, CP2K QM/MM,
 * xtb ONIOM. Fixtures: tests/fixtures/qmmm (written by `python3 tools/make_big_fixture.py qmmm tests/fixtures/qmmm`,
 * expected values in expected.json), plus variants made here (no gridpoint summary, micro-iteration prints, a file cut
 * off mid-step, ORCA prints without the full system / without the input QM atom list / without energy summaries).
 * Run:  node tests/test_qmmm.js */
'use strict';
const fs = require('fs');
const path = require('path');
const ROOT = path.join(__dirname, '..');
for (const f of ['elements.js', 'bigfile.js', 'parsers.js', 'readers-struct.js', 'writers.js', 'readers-qc.js', 'analysis.js']) require(path.join(ROOT, 'src', f));
const BM = globalThis.BM;
const DIR = path.join(__dirname, 'fixtures', 'qmmm');
const EXP = JSON.parse(fs.readFileSync(path.join(DIR, 'expected.json'), 'utf8'));
const read = f => fs.readFileSync(path.join(DIR, f), 'utf8');

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
function same(a, b, msg) { eq(JSON.stringify(a), JSON.stringify(b), msg); }
function test(name, fn) { current = name; try { fn(); } catch (e) { fail++; failures.push('[' + name + '] threw: ' + (e && e.stack || e)); } }
const load = (name, text) => BM.parseFile(name, text === undefined ? read(name) : text);

/** Invariants for every QM/MM dataset: full-system frames, one total per step, no sub-energy as the total. */
function common(ds, exp) {
  eq(ds.meta.qmmm, true, 'meta.qmmm');
  const n = exp.atoms;
  ok(ds.frames.every(f => f.Z.length === n && f.xyz.length === 3 * n), 'every frame has the ' + n + ' atoms of the full system');
  const L = ds.layers;
  ok(L && ['oniom', 'qmmm', 'qm-qm2'].indexOf(L.scheme) >= 0, 'layers.scheme ' + (L && L.scheme));
  for (const k of ['atomLayer', 'frozen', 'active', 'linkAtoms', 'mmTypes', 'mmCharges', 'labels']) ok(L && k in L, 'layers.' + k + ' present');
  if (L.atomLayer) { eq(L.atomLayer.length, n, 'one layer per atom'); ok(L.atomLayer.every(x => /^[HML]$/.test(x)), 'layer letters'); }
  for (const k of ['frozen', 'active', 'mmTypes', 'mmCharges']) if (L[k]) eq(L[k].length, n, k + ' per atom');
  for (const la of L.linkAtoms) ok(la.atom >= 0 && la.atom < n && la.bondedTo >= 0 && la.bondedTo < n && typeof la.element === 'string', 'link atom indices ' + JSON.stringify(la));
  ok(['H', 'M', 'L'].every(k => k in L.labels), 'labels H/M/L');
}
/** The total-energy series: one value per step, strictly decreasing (the fixtures converge monotonically). */
function series(ds, exp) {
  const E = ds.frames.map(f => f.energy);
  arrNear(E, exp.energies, 1e-8, 'frame.energy = total of every step');
  ok(E.every((e, k) => k === 0 || e < E[k - 1] + 1e-12), 'no alternation: energies only go down');
  const cs = BM.analysis.convergenceSeries(ds, { job: 0 });
  eq(cs.steps.length, exp.frames, 'convergenceSeries: one step per optimisation cycle');
  ok(cs.isOpt, 'convergenceSeries: optimisation');
  const s = cs.series.find(x => x.key === 'energy');
  ok(s && s.values.every((v, k) => Math.abs(v - exp.energies[k]) < 1e-8), 'convergenceSeries energy = totals');
  const d = cs.series.find(x => x.key === 'dE');
  ok(d && d.values.every(v => v >= -1e-9) && Math.abs(d.values[d.values.length - 1]) < 1e-9, 'ΔE vs final: positive, zero at the end');
  const span = Math.max(...E) - Math.min(...E);
  ok(span < 0.05, 'energy span ' + span.toFixed(4) + ' Eh (sub-calculation energies would be hundreds of Eh away)');
  return cs;
}

// ------------------------------------------------------------------ Gaussian ONIOM(B3LYP/6-31G(d):Amber)
test('Gaussian ONIOM QM/MM optimisation', () => {
  const f = 'g16_oniom_amber_opt.log', exp = EXP[f], ds = load(f);
  eq(ds.format, 'gaussian-log', 'format');
  eq(ds.frames.length, exp.frames, 'one frame per step (model-system orientation tables are not steps)');
  common(ds, exp);
  series(ds, exp);
  ds.frames.forEach((fr, k) => {
    const want = exp.layerEnergies[k];
    for (const key of Object.keys(want)) near(fr.energies[key], want[key], 1e-9, 'step ' + (k + 1) + ' ' + key);
    near(fr.energies.ONIOM, exp.energies[k], 1e-9, 'step ' + (k + 1) + ' ONIOM');
    eq(fr.energies.SCF, undefined, 'no bare SCF key (sub-calculation)');
  });
  const L = ds.layers;
  eq(L.scheme, 'oniom', 'scheme'); eq(L.atomLayer.join(''), exp.atomLayer.join(''), 'H/L flags');
  same(L.frozen, exp.frozen, 'freeze codes -1');
  same(L.active, exp.frozen.map(x => !x), 'active = not frozen');
  same(L.mmTypes, exp.mmTypes, 'MM types');
  arrNear(L.mmCharges, exp.mmCharges, 1e-12, 'MM charges ("C-CT--0.0932")');
  eq(L.linkAtoms.length, 1, 'one link atom');
  same([L.linkAtoms[0].atom, L.linkAtoms[0].bondedTo, L.linkAtoms[0].element, L.linkAtoms[0].type], [0, 4, 'H', 'HC'], 'link atom H-HC on C1 bonded to C5');
  same(L.labels, exp.labels, 'method labels');
  eq(L.embedding, 'electronic', 'EmbedCharge');
  eq(L.chargeMult.length, 3, 'per-layer charge / multiplicity');
  same(L.chargeMult[1], { charge: 0, multiplicity: 1, level: 'high', system: 'model' }, 'high level on model system');
  ok(ds.frames[ds.frames.length - 1].converged, 'converged at the last step');
  ok(ds.frames.every(fr => fr.opt && fr.opt.maxForce !== undefined), 'convergence table on every step');
  eq(ds.meta.charge, 0, 'charge'); eq(ds.warnings.length, 0, 'no warnings');
  // frozen atoms never move
  const x0 = ds.frames[0].xyz, x5 = ds.frames[5].xyz;
  ok([14, 15, 16].every(a => [0, 1, 2].every(d => x0[3 * a + d] === x5[3 * a + d])), 'frozen waters keep their coordinates');
});

test('Gaussian ONIOM variants: no gridpoint summary, micro-iteration prints, cut-off file', () => {
  const f = 'g16_oniom_amber_opt.log', exp = EXP[f], text = read(f);
  // 1) without the "ONIOM: gridpoint" lines the MM "Energy=" lines give the MM layer energies (calculation order)
  let ds = load(f, text.replace(/^ ONIOM: gridpoint.*\n/gm, ''));
  arrNear(ds.frames.map(x => x.energy), exp.energies, 1e-9, 'totals unchanged');
  ds.frames.forEach((fr, k) => {
    near(fr.energies['low-model'], exp.layerEnergies[k]['low-model'], 1e-9, 'MM model energy from Energy= #1');
    near(fr.energies['low-real'], exp.layerEnergies[k]['low-real'], 1e-9, 'MM real energy from Energy= #2');
  });
  // 2) a real-system print with moved MM atoms after the step's ONIOM energy (micro-iterations) is not a new step,
  //    nor is one before it a separate frame
  const blocks = text.split(/(?= ONIOM: Dipole moment)/);
  const moved = b => b.replace(/^(\s+1[2-4]\s+[18]\s+0\s+)(-?\d+\.\d+)/gm, (m, a, x) => a + (+x + 0.05).toFixed(6).padStart(x.length));
  const orient = text.match(/ {26}Input orientation:[\s\S]*?-{69}\n[\s\S]*?-{69}\n[\s\S]*?-{69}\n/g);
  ok(orient && orient.length >= 12, 'orientation tables found');
  let micro = blocks[0];
  for (let k = 1; k < blocks.length; k++) micro += moved(orient[2 * (k - 1)]) + g_mm('    0.0123') + blocks[k];
  ds = load(f, micro);
  eq(ds.frames.length, exp.frames, 'micro-iteration prints: still one frame per step');
  arrNear(ds.frames.map(x => x.energy), exp.energies, 1e-9, 'micro-iteration prints: totals unchanged');
  const ref = load(f);
  ok(ds.frames.every((fr, k) => fr.xyz.every((v, i) => v === ref.frames[k].xyz[i])), 'the step keeps the geometry of its energy');
  // 3) a file cut off after the QM sub-calculation of the last step: that step has no total, not the QM energy
  const cut = text.slice(0, text.lastIndexOf(' SCF Done:')) + text.slice(text.lastIndexOf(' SCF Done:')).split('\n').slice(0, 3).join('\n') + '\n';
  ds = load(f, cut);
  eq(ds.frames.length, exp.frames, 'cut-off file: all steps');
  eq(ds.frames[ds.frames.length - 1].energy, null, 'cut-off step: energy null (the sub-calculation SCF energy is no total)');
  arrNear(ds.frames.slice(0, -1).map(x => x.energy), exp.energies.slice(0, -1), 1e-9, 'earlier steps unchanged');
});
function g_mm(e) { return ' Energy=   ' + e + ' NIter=   0.\n'; }

test('Gaussian 3-layer ONIOM(B3LYP:PM6:UFF)', () => {
  const f = 'g16_oniom3_opt.log', exp = EXP[f], ds = load(f);
  eq(ds.frames.length, exp.frames, 'one frame per step (model, mid and re-printed real geometries are not steps)');
  common(ds, exp);
  series(ds, exp);
  ds.frames.forEach((fr, k) => {
    for (const key of Object.keys(exp.layerEnergies[k])) near(fr.energies[key], exp.layerEnergies[k][key], 1e-9, 'step ' + (k + 1) + ' ' + key);
  });
  const L = ds.layers;
  eq(L.atomLayer.join(''), exp.atomLayer.join(''), 'H/M/L flags'); eq(L.nLayers, 3, 'three layers');
  same(L.labels, exp.labels, 'labels high/medium/low');
  same(L.linkAtoms.map(l => [l.atom, l.bondedTo, l.element, l.type]), exp.linkAtoms.map(l => [l.atom, l.bondedTo, l.element, l.type]), 'link atoms H->M and M->L');
  same(L.mmTypes, exp.mmTypes, 'UFF types');
  eq(L.mmCharges, null, 'no MM charges given');
  eq(L.frozen, null, 'no freeze codes');
  eq(L.chargeMult.length, 5, 'five charge / multiplicity pairs');
  same(L.chargeMult.map(c => c.level + '-' + c.system), ['low-real', 'med-mid', 'low-mid', 'high-model', 'med-model'], 'levels and systems');
  eq(L.embedding, 'mechanical', 'mechanical embedding');
  // the standard orientation (rotated) replaces the input orientation of the same step
  const x = ds.frames[0].xyz;
  ok(Math.abs(x[3 * 3] - 0) < 0.1 && Math.abs(x[3 * 3 + 1] - 1.51) < 0.1, 'frames use the Standard orientation');
});

test('Gaussian ONIOM opt + Link1 freq with sub-calculation markers (review fixture)', () => {
  const ds = BM.parseFile('g16_oniom_optfreq.log', fs.readFileSync(path.join(__dirname, 'fixtures', 'review', 'g16_oniom_optfreq.log'), 'utf8'));
  eq(ds.frames.length, 2, 'frames');
  arrNear(ds.frames.map(f => f.energy), [-76.421, -76.4223456789], 1e-10, 'ONIOM totals');
  near(ds.frames[0].energies['high-model'], -76.4, 1e-10, '"generating point 2 -- high level on model system" -> high-model');
  near(ds.frames[0].energies['low-real'], -0.051, 1e-10, 'low-real from its marker');
  near(ds.frames[0].energies['low-model'], -0.03, 1e-10, 'low-model (marker and gridpoint)');
  eq(ds.layers.scheme, 'oniom', 'layers');
  same(ds.layers.labels, { H: 'B3LYP/6-31G(d)', M: null, L: 'PM6' }, 'labels');
  eq(ds.layers.atomLayer, null, 'no Cartesian echo: layers unknown (not guessed)');
  ok(ds.warnings.some(w => /layers could not be read/.test(w)), 'warning about the missing layer flags');
  eq(ds.vibrations.frameIndex, 1, 'vibrations on the optimised geometry');
});

// ------------------------------------------------------------------ ORCA
test('ORCA QM/MM optimisation (%qmmm, ORCAFF)', () => {
  const f = 'orca6_qmmm_opt.out', exp = EXP[f], ds = load(f);
  eq(ds.format, 'orca-out', 'format');
  eq(ds.frames.length, exp.frames, 'one frame per cycle (QM-region prints and the final evaluation are not steps)');
  common(ds, exp);
  series(ds, exp);
  ds.frames.forEach((fr, k) => {
    for (const key of ['QM', 'MM', 'QM-MM']) near(fr.energies[key], exp.components[k][key], 1e-9, 'cycle ' + (k + 1) + ' ' + key);
    near(fr.energies['QM/MM'], exp.energies[k], 1e-9, 'QM/MM total');
    near(fr.energies.Final, exp.energies[k], 1e-9, 'FINAL SINGLE POINT ENERGY');
    ok(Math.abs(fr.energies['QM SCF'] - fr.energy) > 0.01, 'QM SCF energy kept apart from the total');
  });
  const L = ds.layers;
  eq(L.scheme, 'qmmm', 'scheme');
  eq(L.atomLayer.join(''), exp.atomLayer.join(''), 'QMAtoms {4:10}');
  same(L.active, exp.active, 'ActiveAtoms {0:13}');
  same(L.frozen, exp.active.map(x => !x), 'frozen = outside the active region');
  same(L.linkAtoms.map(l => [l.atom, l.bondedTo, l.element]), [[0, 4, 'H']], 'link atom from the QM-region print');
  same(L.labels, exp.labels, 'labels');
  eq(L.subsystemOnly, false, 'full system printed');
  ok(ds.frames[ds.frames.length - 1].converged, 'converged');
  same(ds.jobs.map(j => [j.frameStart, j.frameEnd]), [[0, exp.frames - 1]], 'job range after merging the prints');
  eq(ds.warnings.length, 0, 'no warnings');
});

test('ORCA QM/MM variants: no QMAtoms in the echo, no summaries, only QM-region prints', () => {
  const f = 'orca6_qmmm_opt.out', exp = EXP[f], text = read(f);
  // layers from matching the QM-region print against the full system
  let ds = load(f, text.replace(/^\|\s*\d+>\s*QMAtoms \{4:10\} end\n/m, ''));
  eq(ds.layers.atomLayer.join(''), exp.atomLayer.join(''), 'QM atoms found by coordinates');
  same(ds.layers.linkAtoms.map(l => [l.atom, l.bondedTo]), [[0, 4]], 'link atom by coordinates');
  // without the QM/MM energy block the FINAL SINGLE POINT ENERGY is the total; QM SCF never is
  ds = load(f, text.replace(/^E\((QM|MM|QM-MM|QM\/MM)\).*\n/gm, '').replace(/^(Total system size|Number of|Link atom).*\n/gm, ''));
  arrNear(ds.frames.map(x => x.energy), exp.energies, 1e-9, 'totals from FINAL SINGLE POINT ENERGY');
  eq(ds.frames.length, exp.frames, 'frames');
  same(ds.layers.linkAtoms.map(l => [l.atom, l.bondedTo]), [[0, 4]], 'link atom without the setup summary');
  // only the QM region printed: frames are the QM region, flagged
  const noFull = text.replace(/-{33}\nCARTESIAN COORDINATES \(ANGSTROEM\)\n-{33}\n((?:  [A-Z][a-z]?\s+\S+\s+\S+\s+\S+\n){17})\n/g, '');
  ds = load(f, noFull);
  eq(ds.frames.length, exp.frames, 'QM-region only: one frame per cycle');
  ok(ds.frames.every(fr => fr.Z.length === 8), 'QM region + link atom in every frame');
  eq(ds.layers.subsystemOnly, true, 'subsystemOnly flag');
  ok(ds.warnings.some(w => /only the QM region/.test(w)), 'warning');
  arrNear(ds.frames.map(x => x.energy), exp.energies, 1e-9, 'totals');
});

test('ORCA QM/XTB (QM/QM2) optimisation', () => {
  const f = 'orca5_qmxtb_opt.out', exp = EXP[f], ds = load(f);
  eq(ds.frames.length, exp.frames, 'frames');
  common(ds, exp);
  series(ds, exp);
  ds.frames.forEach((fr, k) => {
    for (const key of ['low-real', 'low-model', 'high-model']) near(fr.energies[key], exp.components[k][key], 1e-9, 'cycle ' + (k + 1) + ' ' + key);
    near(fr.energies['QM/QM2'], exp.energies[k], 1e-9, 'QM/QM2 total');
  });
  const L = ds.layers;
  eq(L.scheme, 'qm-qm2', 'scheme qm-qm2');
  eq(L.atomLayer.join(''), exp.atomLayer.join(''), 'QM atoms H, the rest QM2 (M)');
  same(L.labels, exp.labels, 'labels');
  same(L.linkAtoms.map(l => [l.atom, l.bondedTo, l.element]), [[0, 4, 'H']], 'link atom');
  eq(L.active, null, 'no active region given');
});

// ------------------------------------------------------------------ CP2K, xtb
test('CP2K QM/MM geometry optimisation', () => {
  const f = 'cp2k_qmmm_geoopt.out', exp = EXP[f], ds = load(f);
  eq(ds.format, 'cp2k-out', 'format');
  eq(ds.frames.length, 1, 'only the whole system (FIST) is a geometry; the QUICKSTEP QM region is not');
  common(ds, exp);
  eq(ds.layers.atomLayer.join(''), exp.atomLayer.join(''), 'QM atoms located through the shifted QM box');
  same(ds.layers.linkAtoms.map(l => [l.atom, l.bondedTo, l.element]), [[0, 4, 'H']], 'link atom');
  arrNear(ds.optSteps.map(s => s.energy), exp.energies, 1e-8, 'QMMM totals per step');
  near(ds.frames[0].energies['QM/MM'], -45.1234, 1e-9, 'first ENERGY| Total FORCE_EVAL ( QMMM ) energy');
  const cs = BM.analysis.convergenceSeries(ds, { job: 0 });
  eq(cs.source, 'optSteps', 'convergence from the optimisation steps');
  ok(cs.series.find(s => s.key === 'energy').values.every((v, k) => Math.abs(v - exp.energies[k]) < 1e-8), 'series');
});

test('xtb ONIOM', () => {
  const f = 'xtb_oniom_water.out', exp = EXP[f], ds = load(f);
  common(ds, exp);
  eq(ds.layers.atomLayer.join(''), exp.atomLayer.join(''), '--oniom inner atoms 1,2');
  same(ds.layers.labels, exp.labels, 'labels');
});

// ------------------------------------------------------------------ big-file path
test('QM/MM outputs read from lines give the same datasets', () => {
  const ser = ds => JSON.stringify(ds, (k, v) => (ArrayBuffer.isView(v) ? Array.from(v) : v));
  for (const f of Object.keys(EXP)) {
    const t = read(f);
    eq(ser(BM.parseLines(f, t.split(/\r\n|\r|\n/))), ser(BM.parseFile(f, t)), f + ': parseLines = parseFile');
  }
});

test('non-QM/MM outputs are not flagged', () => {
  for (const f of ['g16_water_optfreq.log', 'orca5_water_optfreq.out']) {
    const ds = BM.parseFile(f, fs.readFileSync(path.join(__dirname, 'fixtures', f), 'utf8'));
    ok(!ds.meta.qmmm && !ds.layers, f + ': no layers');
  }
});

console.log('');
for (const m of failures) console.log('FAIL ' + m);
console.log('BlueMolecule QM/MM tests: ' + pass + ' passed, ' + fail + ' failed');
process.exit(fail ? 1 : 0);
