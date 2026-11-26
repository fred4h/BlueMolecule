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
for (const f of ['elements.js', 'bigfile.js', 'parsers.js', 'readers-struct.js', 'writers.js', 'readers-qc.js', 'analysis.js', 'condense.js']) require(path.join(ROOT, 'src', f));
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

// ------------------------------------------------------------------ layers the output does not list (2.6.1)
/**
 * The Amber fixture as a job that takes its molecule from the checkpoint (Geom=AllCheck): Gaussian then prints
 * "Structure from the checkpoint file" rows "El,freeze,x,y,z" without layers, as in real outputs, and NQM counts the
 * model-system atoms other than link atoms (high layer + replaced atoms, with Gaussian's own MM). opts: nqm, route,
 * rows (a function of the echo row tokens).
 */
function allcheck(text, opts) {
  opts = opts || {};
  const head = text.indexOf(' Symbolic Z-matrix:'), end = text.indexOf('\n \n', head);
  const rows = text.slice(head, end).split('\n').filter(l => /^ [A-Z][a-z]?-/.test(l)).map(l => l.trim().split(/\s+/));
  const row = opts.rows || (t => t[0].split('-')[0] + ',' + t[1] + ',' + t[2] + ',' + t[3] + ',' + t[4]);
  const block = ' Structure from the checkpoint file:  "job0.chk"\n' +
    ' Charge =  0 Multiplicity = 1 for low   level calculation on real  system.\n' +
    ' Charge =  0 Multiplicity = 1 for high  level calculation on model system.\n' +
    ' Charge =  0 Multiplicity = 1 for low   level calculation on model system.\n' +
    ' Redundant internal coordinates found in file.\n' + rows.map(t => ' ' + row(t)).join('\n') + '\n Recover connectivity data from disk.';
  let out = (text.slice(0, head) + block + text.slice(end)).replace(/ %chk=job\.chk\n/, ' %oldchk=job0.chk\n %chk=job.chk\n')
    .replace(/geom=connectivity/, 'geom=allcheck')
    .replace(/NAtoms=\s+17 NQM=\s+\d+ NQMF=\s+\d+ NMMI=\s+\d+ NMMIF=\s+\d+/, 'NAtoms=     17 NQM=        ' + String(opts.nqm || 8).padStart(2) + ' NQMF=       0 NMMI=      9 NMMIF=      3');
  if (opts.route) out = out.replace(/oniom\(b3lyp\/6-31g\(d\):amber\)/, opts.route);
  return out;
}

test('Gaussian ONIOM from the checkpoint (Geom=AllCheck): layers rebuilt from the link atoms', () => {
  const f = 'g16_oniom_amber_opt.log', exp = EXP[f], text = allcheck(read(f));
  ok(!/Symbolic Z-matrix/.test(text) && /Structure from the checkpoint file/.test(text), 'variant has no molecule echo');
  const ds = load('restart.log', text);
  common(ds, exp);
  series(ds, exp);
  const L = ds.layers;
  eq(L.source, 'rebuilt', 'source');
  eq(L.atomLayer.join(''), exp.atomLayer.join(''), 'high layer = the atoms bonded to the QM side of the link atom (NQM 8 = 7 + 1 replaced)');
  same(L.rebuilt, { qm: 7, links: 1 }, 'rebuilt counts');
  same(L.frozen, exp.frozen, 'frozen atoms from the freeze codes of the checkpoint rows');
  same(L.linkAtoms.map(l => [l.atom, l.bondedTo, l.element]), [[0, 4, 'H']], 'link atom from "ONIOM: Cut between"');
  eq(L.mmTypes, null, 'no MM types (not printed)');
  same(L.labels, exp.labels, 'method labels from the route');
  eq(L.chargeMult.length, 3, 'per-layer charge / multiplicity from the checkpoint block');
  eq(ds.meta.oldchk, 'job0.chk', '%oldchk'); eq(ds.meta.chk, 'job.chk', '%chk');
  eq(ds.warnings.length, 1, 'one reader note');
  ok(/checkpoint file \(Geom=Check\).*rebuilt them from the 1 link atom: the 7 atoms.*continued from job0\.chk.*take its layers in the Info tab/.test(ds.warnings[0]), 'note: ' + ds.warnings[0]);
  // read from lines (big files) the same
  const ser = x => JSON.stringify(x, (k, v) => (ArrayBuffer.isView(v) ? Array.from(v) : v));
  eq(ser(BM.parseLines('restart.log', text.split('\n'))), ser(ds), 'parseLines = parseFile');
});

test('Gaussian ONIOM from the checkpoint: layers left unknown when the link atoms cannot account for them', () => {
  const f = 'g16_oniom_amber_opt.log', exp = EXP[f];
  // a water in the QM region: 10 QM atoms, but only 7 are bonded to the link atom's QM end
  let ds = load('restart.log', allcheck(read(f), { nqm: 11 }));
  eq(ds.layers.atomLayer, null, 'not rebuilt: layers unknown, not guessed');
  eq(ds.layers.source, null, 'no source');
  same(ds.layers.frozen, exp.frozen, 'freeze codes still read');
  eq(ds.layers.linkAtoms.length, 1, 'link atom still known');
  ok(ds.warnings.length === 1 && /does not list which atoms.*1 link atom alone does not give the whole QM region.*continued from job0\.chk.*take its layers in the Info tab/.test(ds.warnings[0]), 'note: ' + ds.warnings[0]);
  arrNear(ds.frames.map(x => x.energy), exp.energies, 1e-9, 'energies unchanged');
  // QM:QM (no Gaussian MM low level): NQM is not the model-system count, so nothing is rebuilt
  ds = load('restart.log', allcheck(read(f), { route: 'oniom(b3lyp/6-31g(d):pm6)' }));
  eq(ds.layers.atomLayer, null, 'QM:QM: not rebuilt');
  eq(ds.layers.labels.L, 'PM6', 'low level PM6');
  // too many counted atoms in reach: a boundary bond without a link atom (the walk leaks): not rebuilt
  ds = load('restart.log', allcheck(read(f), { nqm: 6 }));
  eq(ds.layers.atomLayer, null, 'fewer QM atoms counted than bonded: not rebuilt');
});

test('Gaussian ONIOM: checkpoint rows that carry layers, and molecule echoes in Z-matrix form', () => {
  const f = 'g16_oniom_amber_opt.log', exp = EXP[f], text = read(f);
  // rows "El-Type-Charge,freeze,x,y,z,Layer[,link,host]" are read like the echo
  let ds = load('restart.log', allcheck(text, { rows: t => t.slice(0, 5).concat(t.slice(5)).join(',') }));
  eq(ds.layers.source, 'checkpoint', 'layers from the checkpoint rows');
  eq(ds.layers.atomLayer.join(''), exp.atomLayer.join(''), 'H/L flags');
  same(ds.layers.mmTypes, exp.mmTypes, 'MM types'); same(ds.layers.frozen, exp.frozen, 'frozen');
  same(ds.layers.linkAtoms.map(l => [l.atom, l.bondedTo, l.type]), [[0, 4, 'HC']], 'link atom with its type');
  eq(ds.warnings.length, 0, 'no note');
  // an echo whose rows are a Z-matrix (no Cartesian coordinates): rebuilt from the link atoms as well
  const zm = text.replace(/^( [A-Z][a-z]?-\S+\s+-?\d)\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)( [HL])/gm, '$1   R1   A1   D1 $5')
    .replace(/NAtoms=\s+17 NQM=\s+\d+/, 'NAtoms=     17 NQM=        8');
  ds = load(f, zm);
  eq(ds.layers.source, 'rebuilt', 'Z-matrix echo: rebuilt');
  eq(ds.layers.atomLayer.join(''), exp.atomLayer.join(''), 'Z-matrix echo: H/L flags');
  ok(/not in Cartesian form.*rebuilt them/.test(ds.warnings[0] || ''), 'note: ' + ds.warnings[0]);
});

/** Population analyses of QM sub-calculations added to a fixture: Gaussian prints a Mulliken charge for every atom of the
 *  real system there. tables: [{step (0-based), sys ('high level on model' by default), q, el, nActive (null: no
 *  "NActive=" line)}], each inserted before that step's "ONIOM: calculating energy.". */
function withPop(text, tables) {
  const parts = text.split(' ONIOM: calculating energy.\n');
  for (const t of tables) {
    const n = t.q.length, w = v => String(v).padStart(5);
    parts[t.step] += ' ONIOM: generating point  2 -- ' + (t.sys || 'high level on model') + ' system.\n' +
      (t.nActive === null ? '' : ' NAtoms=' + w(n) + ' NActive=' + w(t.nActive) + ' NUniq=' + w(t.nActive) + ' SFac= 1.00D+00 NAtFMM=   60 NAOKFM=F Big=F\n') +
      ' Mulliken charges:\n               1\n' + t.q.map((q, a) => String(a + 1).padStart(6) + '  ' + t.el[a].padEnd(2) + q.toFixed(6).padStart(11)).join('\n') +
      '\n Sum of Mulliken charges =   0.00000\n';
  }
  return parts.join(' ONIOM: calculating energy.\n');
}
/** Condensed reading (the big-file path) of text, then BM.parseLines. */
function condensed(name, text) {
  const c = BM.condense.create('gaussian-log', {});
  const lines = text.split('\n');
  for (let k = 0; k < lines.length - 1; k++) c.push(lines[k]);
  return BM.parseLines(name, c.finish(lines[lines.length - 1]));
}

test('Gaussian ONIOM from the checkpoint: layers from the population analysis of the high-level calculation', () => {
  const f = 'g16_oniom_amber_opt.log', exp = EXP[f], base = allcheck(read(f)), want = exp.atomLayer.join('');
  const el = Array.from(load(f).frames[0].Z, z => BM.el.symbol(z));
  const mmq = [0, 0.04, 0.04, 0.04, 0, 0, 0, 0, 0, 0, 0, -0.834, 0.417, 0.417, -0.834, 0.417, 0.417];
  const qmq = [0.07, 0, 0, 0, -0.31, 0.15, 0.15, 0.62, -0.52, -0.58, 0.43];
  // the model system's table: the QM atoms and the link atom (printed as H, atom 1) have charges, which change from
  // step to step; the others have 0 (mechanical embedding) or their fixed MM charges (electronic embedding)
  const table = (step, o) => {
    o = o || {};
    const e = el.slice(); if (!o.sameEl) e[0] = 'H';
    return { step, el: e, nActive: o.nActive === undefined ? 8 : o.nActive,
      q: el.map((_, a) => (want[a] === 'H' || a === 0 ? qmq[a] + step * 0.0013 * (a + 1) : o.ee ? mmq[a] : 0)) };
  };
  let text = withPop(base, [table(0)]);
  let ds = load('restart.log', text);
  common(ds, exp);
  series(ds, exp);
  let L = ds.layers;
  eq(L.source, 'population', 'source');
  eq(L.atomLayer.join(''), want, 'layers: the atoms with Mulliken charges in the model system, not the link atom (printed as H)');
  same(L.population, { how: 'charges', qm: 7, mid: 0, links: 1, tables: 1 }, 'population summary');
  same(L.linkAtoms.map(l => [l.atom, l.bondedTo, l.element]), [[0, 4, 'H']], 'link atom from "ONIOM: Cut between"');
  same(L.frozen, exp.frozen, 'frozen atoms from the checkpoint rows');
  eq(ds.warnings.length, 1, 'one reader note');
  ok(/checkpoint file \(Geom=Check\).*read them from the high-level calculation on the model system: only the 7 QM atoms and the atom replaced by a link atom have Mulliken charges there, as many atoms as Gaussian counted.*continued from job0\.chk.*take its layers in the Info tab/.test(ds.warnings[0]), 'note: ' + ds.warnings[0]);
  const ser = x => JSON.stringify(x, (k, v) => (ArrayBuffer.isView(v) ? Array.from(v) : v));
  eq(ser(BM.parseLines('restart.log', text.split('\n'))), ser(ds), 'parseLines = parseFile');
  eq(condensed('restart.log', text).layers.atomLayer.join(''), want, 'condensed reading (big files): the same layers');

  // electronic embedding: the other atoms keep their MM charges, the QM charges change between the two tables
  text = withPop(base, [table(0, { ee: true }), table(5, { ee: true })]);
  ds = load('restart.log', text);
  L = ds.layers;
  eq(L.source, 'population', 'embedding: source');
  eq(L.atomLayer.join(''), want, 'embedding: layers');
  same(L.population, { how: 'changed', qm: 7, mid: 0, links: 1, tables: 2 }, 'embedding: the charges that change');
  ok(/electronic embedding the other atoms keep their MM charges there, and the Mulliken charges that change from one geometry to the next \(2 tables\) are those of the 7 QM atoms and the atom replaced by a link atom, as many atoms as Gaussian counted/.test(ds.warnings[0]), 'embedding note: ' + ds.warnings[0]);
  eq(condensed('restart.log', text).layers.population.how, 'changed', 'condensed reading: embedding');
  // one table with embedding cannot tell MM charges from QM charges: the link-atom rebuild (NQM 8) applies instead
  eq(load('restart.log', withPop(base, [table(0, { ee: true })])).layers.source, 'rebuilt', 'embedding, one table: link-atom rebuild');
  // two tables of the same geometry: nothing changes, nothing is read (and the rebuild cannot account for NQM 11)
  ds = load('restart.log', withPop(allcheck(read(f), { nqm: 11 }), [table(0, { ee: true }), Object.assign(table(0, { ee: true }), { step: 5 })]));
  eq(ds.layers.atomLayer, null, 'embedding, unchanged charges: unknown');

  // a water in the QM region, bonded to no link atom: the link-atom rebuild cannot find it, the population analysis does
  const water = step => { const t = table(step, { nActive: 11 }); t.q = t.q.map((v, a) => (a >= 11 && a <= 13 ? [-0.8, 0.4, 0.4][a - 11] + 0.001 * step : v)); return t; };
  ds = load('restart.log', withPop(allcheck(read(f), { nqm: 11 }), [water(0)]));
  eq(ds.layers.source, 'population', 'water in the QM region: read from the population analysis');
  eq(ds.layers.atomLayer.join(''), want.slice(0, 11) + 'HHH' + want.slice(14), 'the water is in the high layer');
  // a link atom printed with its host's element is known from "ONIOM: Cut between"
  eq(load('restart.log', withPop(base, [table(0, { sameEl: true })])).layers.atomLayer.join(''), want, 'same-element link atom');
  // no "ONIOM: Cut between" lines: the link atom from the element printed, bonded to the nearest QM atom
  ds = load('restart.log', withPop(base.replace(/^ ONIOM: Cut between.*\n/gm, ''), [table(0)]));
  eq(ds.layers.atomLayer.join(''), want, 'no cut lines: layers');
  same(ds.layers.linkAtoms.map(l => [l.atom, l.bondedTo, l.element]), [[0, 4, 'H']], 'no cut lines: link atom');
  // the count must match Gaussian's: NActive, or NQM (Gaussian's own MM) when no NActive line is printed
  eq(load('restart.log', withPop(base, [table(0, { nActive: 9 })])).layers.source, 'rebuilt', 'NActive differs: not used (link-atom rebuild)');
  eq(load('restart.log', withPop(base, [table(0, { nActive: null })])).layers.source, 'population', 'no NActive line: NQM counts the model system');
  // an input echo still comes first
  eq(load(f, withPop(read(f), [table(0)])).layers.source, 'input', 'with the echo: source input');
});

test('Gaussian 3-layer ONIOM from the checkpoint: high and medium layers from the population analyses', () => {
  const f = 'g16_oniom3_opt.log', exp = EXP[f], text = read(f), want = exp.atomLayer.join('');
  const head = text.indexOf(' Symbolic Z-matrix:'), end = text.indexOf('\n \n', head), echo = text.slice(head, end).split('\n');
  const rows = echo.filter(l => /^ [A-Z][a-z]?-/.test(l)).map(l => l.trim().split(/\s+/));
  const base = text.slice(0, head) + ' Structure from the checkpoint file:  "job0.chk"\n' + echo.filter(l => /^ Charge =/.test(l)).join('\n') +
    '\n Redundant internal coordinates found in file.\n' + rows.map(t => ' ' + t[0].split('-')[0] + ',0,' + t.slice(1, 4).join(',')).join('\n') +
    '\n Recover connectivity data from disk.' + text.slice(end);
  const el = Array.from(load(f).frames[0].Z, z => BM.el.symbol(z));
  // mid system: H and M atoms, and atom 7 replaced by a link atom; model system: H atoms, and atom 4 replaced
  const tab = (sys, nActive, inside, host) => ({ step: 0, sys, nActive, el: el.map((e, a) => (a === host ? 'H' : e)),
    q: el.map((_, a) => (inside(a) || a === host ? 0.1 + 0.01 * a : 0)) });
  let ds = load('restart3.log', withPop(base, [tab('med level on mid', 7, a => a < 6, 6), tab('high level on model', 4, a => a < 3, 3)]));
  const L = ds.layers;
  eq(L.source, 'population', 'source');
  eq(L.atomLayer.join(''), want, 'high, medium and low layers');
  same(L.population, { how: 'charges', qm: 3, mid: 3, links: 2, tables: 1 }, 'counts');
  same(L.linkAtoms.map(l => [l.atom, l.bondedTo, l.element]), exp.linkAtoms.map(l => [l.atom, l.bondedTo, l.element]), 'link atoms, bonded to the nearest atom of the layer above');
  ok(/high- and medium-level calculations: only the 3 high-layer and 3 medium-layer atoms and the 2 atoms replaced by link atoms have Mulliken charges there/.test(ds.warnings[0] || ''), 'note: ' + ds.warnings[0]);
  ds = load('restart3.log', withPop(base, [tab('high level on model', 4, a => a < 3, 3)]));
  eq(ds.layers.atomLayer, null, 'no medium-level table: unknown (medium cannot be told from low)');
});

test('BM.oniomLayersFromPopulation', () => {
  const Z = [6, 1, 1, 6, 8, 1, 1], el = ['C', 'H', 'H', 'H', 'O', 'H', 'H'];   // atom 3 (C) replaced by a link atom (H)
  const q0 = [-0.3, 0.1, 0.1, 0.1, 0, 0, 0];
  const P = (q, e, nActive, extra) => ({ 'high-model': Object.assign({ q: Float64Array.from(q), el: e, nActive, tables: 1, changed: new Uint8Array(q.length) }, extra || {}) });
  const F = BM.oniomLayersFromPopulation;
  let r = F(Z, P(q0, el, 4));
  eq(r && r.atomLayer.join(''), 'HHHLLLL', 'QM atoms: those with charges, without the link atom (element H)');
  same(r && r.hosts, [{ atom: 3, element: 'H', system: 'model' }], 'replaced atom');
  eq(F(Z, P(q0, el, 5)), null, 'count differs from NActive: null');
  eq(F(Z, P(q0.slice(0, 6), el.slice(0, 6), 4)), null, 'table shorter than the system: null');
  eq(F(Z, P(q0, ['??'].concat(el.slice(1)), 4)), null, 'unreadable element: null');
  eq(F(Z, {}), null, 'no table: null');
  eq(F(Z, P([0, 0, 0, 0.1, 0, 0, 0], el, 1)), null, 'only a link atom: null');
  const ee = [-0.3, 0.1, 0.1, 0.1, -0.834, 0.417, 0.417];
  r = F(Z, P(ee, el, 4, { tables: 2, changed: Uint8Array.from([1, 1, 1, 1, 0, 0, 0]) }));
  eq(r && r.how + ' ' + r.atomLayer.join(''), 'changed HHHLLLL', 'embedding: the atoms whose charges changed');
  eq(F(Z, P(ee, el, 4)), null, 'embedding, one table: null');
  eq(F(Z, P(q0, el, null), { nModel: 4 }).qm, 3, 'no NActive: the model-system count (NQM) instead');
  eq(F(Z, P(q0, el, null)), null, 'no count at all: null');
  eq(F(Z, P(q0, ['C', 'H', 'H', 'C', 'O', 'H', 'H'], 4), { cutHosts: [3] }).atomLayer.join(''), 'HHHLLLL', 'replaced atom printed with its own element: from the cut lines');
  eq(F(Z, P(q0, el, 4), { cutHosts: [4] }), null, 'a cut atom outside the model system: null');
  eq(F(Z, P(q0, el, 4), { three: true }), null, '3 layers without the medium-level table: null');
});

test('BM.qmRegionFromLinks', () => {
  // ethanol-like chain C0-C1-O2-H3 with H on C0 (4..6) and C1 (7, 8); link between C0 (MM) and C1 (QM)
  const Z = [6, 6, 8, 1, 1, 1, 1, 1, 1];
  const X = [0, 0, 0, 1.52, 0, 0, 2.0, 1.35, 0, 2.95, 1.35, 0, -0.36, 1.02, 0, -0.36, -0.51, 0.89, -0.36, -0.51, -0.89, 1.88, -0.51, 0.89, 1.88, -0.51, -0.89];
  const qm = BM.qmRegionFromLinks(Z, X, [{ atom: 0, bondedTo: 1 }], 5);
  same(qm, [1, 2, 3, 7, 8], 'QM side of the link: C1, O2, H3 and the two H on C1');
  eq(BM.qmRegionFromLinks(Z, X, [{ atom: 0, bondedTo: 1 }], 6), null, 'count mismatch: null');
  eq(BM.qmRegionFromLinks(Z, X, [], 5), null, 'no link atoms: null');
  eq(BM.qmRegionFromLinks(Z, X, [{ atom: 0, bondedTo: 99 }], 5), null, 'bad index: null');
  // a hydrogen 1.0 Å from another hydrogen (short contact) is no path
  const Z2 = Z.concat([1, 8, 1]), X2 = X.concat([3.4, 2.2, 0, 4.0, 2.9, 0, 4.9, 2.6, 0]);
  same(BM.qmRegionFromLinks(Z2, X2, [{ atom: 0, bondedTo: 1 }], 5), [1, 2, 3, 7, 8], 'the walk stops at hydrogens');
});

test('ORCA QM/MM without the QM atom list: rebuilt from the link atoms, else unknown (not all MM)', () => {
  const f = 'orca6_qmmm_opt.out', exp = EXP[f], text = read(f);
  const noList = text.replace(/^\|\s*\d+>\s*QMAtoms \{4:10\} end\n/m, '')
    .replace(/-{33}\nCARTESIAN COORDINATES \(ANGSTROEM\)\n-{33}\n((?:  [A-Z][a-z]?\s+\S+\s+\S+\s+\S+\n){8})\n/g, '');
  let ds = load(f, noList);
  eq(ds.layers.source, 'rebuilt', 'rebuilt from the link atom rows and "Number of QM atoms"');
  eq(ds.layers.atomLayer.join(''), exp.atomLayer.join(''), 'QM atoms');
  ok(/not listed in this output.*rebuilt them from the 1 link atom: the 7 atoms/.test(ds.warnings.join(' ')), 'note');
  ds = load(f, noList.replace(/^Link atom.*\n/gm, ''));
  eq(ds.layers.atomLayer, null, 'without link atoms: unknown, not every atom MM');
  ok(/not listed in this output.*take its layers in the Info tab/.test(ds.warnings.join(' ')), 'note: take them from another document');
  eq(load(f).layers.source, 'input', 'with QMAtoms: source input');
});

test('Gaussian ONIOM inputs: layers, link atoms, freeze codes, residues (GaussView style)', () => {
  const gjf = ['%chk=test.chk', '#p oniom(b3lyp/6-31g(d):amber)=embedcharge opt geom=connectivity', '', 'GaussView-style ONIOM input', '',
    '0 1 0 1 0 1',
    ' N-N3-0.1592(PDBName=N,ResName=SER,ResNum=1)     0   -1.200000    0.300000    0.000000 L',
    ' H-H-0.1849(PDBName=H1,ResName=SER,ResNum=1)     0   -1.900000    1.000000    0.000000 L',
    ' C-CX-0.0567(PDBName=CA,ResName=SER,ResNum=1)    0    0.000000    0.000000    0.000000 L H-HC 4',
    ' C-2C-0.2596(PDBName=CB,ResName=SER,ResNum=1)    0    0.800000    1.300000    0.000000 H',
    ' O-OH--0.6714(PDBName=OG,ResName=SER,ResNum=1)  -1    2.100000    1.100000    0.000000 H',
    ' H-HO-0.4239(PDBName=HG,ResName=SER,ResNum=1)    0    2.500000    1.950000    0.000000 H',
    '', ' 1 2 1.0 3 1.0', ' 2', ' 3 4 1.0', ' 4 5 1.0', ' 5 6 1.0', ' 6', '', ''].join('\n');
  let ds = BM.parseFile('serine.gjf', gjf);
  eq(ds.format, 'gaussian-input', 'opens (residue information with commas inside parentheses)');
  same(Array.from(ds.frames[0].Z), [7, 1, 6, 6, 8, 1], 'elements');
  const L = ds.layers;
  eq(ds.meta.qmmm, true, 'qmmm'); eq(L.source, 'input', 'source');
  eq(L.atomLayer.join(''), 'LLLHHH', 'layers');
  same(L.linkAtoms.map(l => [l.atom, l.bondedTo, l.element, l.type]), [[2, 3, 'H', 'HC']], 'link atom');
  same(L.frozen, [false, false, false, false, true, false], 'freeze code -1');
  same(L.mmTypes, ['N3', 'H', 'CX', '2C', 'OH', 'HO'], 'MM types');
  arrNear(L.mmCharges, [0.1592, 0.1849, 0.0567, 0.2596, -0.6714, 0.4239], 1e-12, 'MM charges');
  same(L.chargeMult.map(c => c.level + '-' + c.system), ['low-real', 'high-model', 'low-model'], 'charge / multiplicity pairs');
  eq(L.embedding, 'electronic', 'EmbedCharge');
  same(ds.frames[0].resNames, ['SER', 'SER', 'SER', 'SER', 'SER', 'SER'], 'residue names');
  same(ds.frames[0].atomNames, ['N', 'H1', 'CA', 'CB', 'OG', 'HG'], 'atom names');
  eq(ds.warnings.length, 0, 'no warnings');
  // the same rows separated by commas, and three layers with six charge / multiplicity pairs
  const three = gjf.replace('amber)=embedcharge', 'pm6:amber)').replace('0 1 0 1 0 1', '0 1 0 1 0 1 0 1 0 1 0 1')
    .replace(/(PDBName=N,ResName=SER,ResNum=1\))     0   -1.200000    0.300000    0.000000 L/, '$1,0,-1.2,0.3,0.0,M');
  ds = BM.parseFile('three.gjf', three);
  eq(ds.layers.atomLayer.join(''), 'MLLHHH', 'comma-separated row and a medium layer');
  eq(ds.layers.nLayers, 3, 'three layers'); eq(ds.layers.labels.M, 'PM6', 'medium level');
  same(ds.layers.chargeMult.map(c => c.level + '-' + c.system), ['low-real', 'med-mid', 'low-mid', 'high-model', 'med-model', 'low-model'], 'six pairs');
  // a plain (non-ONIOM) input with residue information also opens, with its residues
  ds = BM.parseFile('plain.gjf', ['#p b3lyp/6-31g(d) opt', '', 'plain', '', '0 1',
    ' O(PDBName=O,ResName=HOH,ResNum=5)  0.0 0.0 0.117', ' H(PDBName=H1,ResName=HOH,ResNum=5) 0.0 0.757 -0.469', ' H(PDBName=H2,ResName=HOH,ResNum=5) 0.0 -0.757 -0.469', '', ''].join('\n'));
  eq(ds.frames[0].Z.length, 3, 'plain input opens'); same(ds.frames[0].resNames, ['HOH', 'HOH', 'HOH'], 'plain input residues'); ok(!ds.layers && !ds.meta.qmmm, 'no layers');
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
