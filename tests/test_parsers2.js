#!/usr/bin/env node
/* BlueMolecule tests for the Gaussian / ORCA reader additions (convergence data, job list, scans, IRC, NMR, ECD),
 * the ORCA input reader and the Molden / cube / Gaussian-input / ORCA-input save formats.
 * Run:  node tests/test_parsers2.js */
'use strict';
const fs = require('fs');
const path = require('path');
const ROOT = path.join(__dirname, '..');
require(path.join(ROOT, 'src', 'elements.js'));
require(path.join(ROOT, 'src', 'parsers.js'));
require(path.join(ROOT, 'src', 'writers.js'));
const BM = globalThis.BM;
const FIX = path.join(__dirname, 'fixtures');
const read = f => fs.readFileSync(path.join(FIX, f), 'utf8');
const load = f => BM.parseFile(path.basename(f), read(f));

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
const dist = (x, i, j) => Math.hypot(x[3 * i] - x[3 * j], x[3 * i + 1] - x[3 * j + 1], x[3 * i + 2] - x[3 * j + 2]);
const OPT4 = ['maxForce', 'rmsForce', 'maxDisp', 'rmsDisp'];

// Detection of every existing Gaussian / ORCA / other fixture, before the extra readers are loaded.
const ALL = [];
for (const d of ['', 'review']) for (const f of fs.readdirSync(path.join(FIX, d))) if (/\.(log|out|gjf|com|xyz|cube|molden|fchk|sdf|pdb|gro)$/.test(f)) ALL.push(path.join(d, f));
const before = ALL.map(f => BM.detectFormat(path.basename(f), read(f)));

// ===================================================================== Gaussian
test('Gaussian opt+freq: convergence table on every step, job list', () => {
  const ds = load('g16_water_optfreq.log');
  eq(ds.frames.length, 4, 'frames');
  ds.frames.forEach((f, k) => {
    ok(f.opt && OPT4.every(q => Number.isFinite(f.opt[q]) && Number.isFinite(f.opt[q + 'Thr']) && typeof f.opt[q + 'Conv'] === 'boolean'), 'step ' + (k + 1) + ' has all four criteria with thresholds and flags');
    ok(Number.isFinite(f.opt.predictedChange), 'step ' + (k + 1) + ' predicted change');
  });
  ok(ds.frames[0].opt.energyChange === undefined, 'no energy change for the first step');
  near(ds.frames[2].opt.energyChange, ds.frames[2].energy - ds.frames[1].energy, 1e-12, 'energy change = E(k) - E(k-1)');
  eq(ds.jobs.length, 2, 'opt freq runs as two Link1 jobs');
  eq(JSON.stringify(ds.jobs[0].kinds), '["opt","freq"]', 'job 1 kinds');
  eq(JSON.stringify(ds.jobs[1].kinds), '["freq"]', 'job 2 kinds');
  eq(ds.jobs[0].frameStart + ':' + ds.jobs[0].frameEnd, '0:3', 'job 1 frames');
  eq(ds.jobs[1].frameStart + ':' + ds.jobs[1].frameEnd, '3:3', 'job 2 = final geometry');
  ok(ds.jobs.every(j => j.normalTermination === true), 'both jobs terminated normally');
  ok(/^#p opt freq/.test(ds.jobs[0].label) && /Freq/.test(ds.jobs[1].label), 'labels are the routes');
});

test('Gaussian Link1: single point then optimisation', () => {
  const ds = load('review/g16_link1_symm_then_nosymm.log');
  eq(ds.jobs.length, 2, 'two jobs');
  eq(JSON.stringify(ds.jobs.map(j => j.kinds)), '[["sp"],["opt"]]', 'kinds');
  eq(ds.jobs[1].frameStart, 1, 'job 2 starts after the single point');
  eq(ds.jobs[1].frameEnd, ds.frames.length - 1, 'job 2 ends at the last frame');
  const one = load('g16_water_tddft.log');
  eq(one.jobs.length, 1, 'single-job file -> one entry');
  eq(JSON.stringify(one.jobs[0].kinds), '["td"]', 'TD job kind');
});

test('Gaussian ModRedundant relaxed scan', () => {
  const ds = load('qc/g16_water_modred_scan.log');
  const pts = ds.frames.filter(f => f.scan);
  eq(pts.length, 2, 'two converged scan points');
  ok(pts.every(f => f.converged && f.scan.coord === 'A(2,1,3)'), 'coordinate label and converged flag');
  arrNear(pts.map(f => f.scan.value), [100, 105], 1e-9, 'values from the Optimized Parameters tables');
  ok(ds.jobs[0].kinds.indexOf('scan') >= 0 && ds.meta.jobTypes.indexOf('scan') >= 0, 'scan job kind');
  ok(ds.frames.every(f => f.opt && OPT4.every(q => q in f.opt)), 'every scan step has convergence data');
  const plain = load('g16_water_scan.log');   // no ModRedundant echo: values come from geometry only if a scan coordinate is known
  ok(plain.frames.every(f => !f.scan), 'no scan coordinate known -> no frame.scan');
});

test('Gaussian IRC', () => {
  const ds = load('qc/g16_hcn_irc.log');
  eq(ds.frames.length, 5, 'TS + 2 forward + 2 reverse points');
  eq(JSON.stringify(ds.frames[0].irc), '{"path":0,"point":0,"coord":0}', 'TS');
  eq(JSON.stringify(ds.frames.map(f => f.irc.path)), '[0,1,1,2,2]', 'paths');
  arrNear(ds.frames.map(f => f.irc.coord), [0, 0.1, 0.2, -0.1, -0.2], 1e-9, 'signed reaction coordinate');
  eq(JSON.stringify(ds.jobs[0].kinds), '["irc"]', 'IRC job');
});

test('Gaussian NMR shieldings and spin-spin couplings', () => {
  const ds = load('qc/g16_water_nmr.log');
  arrNear(ds.nmr.iso, [325.5673, 31.2341, 31.2341], 1e-9, 'isotropic shielding');
  arrNear(ds.nmr.aniso, [51.2234, 19.5672, 19.5672], 1e-9, 'anisotropy');
  eq(ds.nmr.method, 'GIAO', 'method');
  eq(ds.nmr.couplings.length, 3, 'three pairs (J, not K)');
  near(ds.nmr.couplings.find(c => c.i === 2 && c.j === 1).J, -7.12345, 1e-9, 'J(H,H)');
  near(ds.nmr.couplings.find(c => c.i === 1 && c.j === 0).J, -75.4321, 1e-9, 'J(O,H)');
  ok(ds.meta.jobTypes.indexOf('nmr') >= 0 && ds.jobs[0].kinds[0] === 'nmr', 'NMR job');
});

test('Gaussian ECD rotatory strengths', () => {
  const ds = load('qc/g16_water_ecd.log');
  arrNear(ds.excitedStates.map(s => s.R), [0, -0.2552, 0, 1.1522], 1e-9, 'R(length)');
  arrNear(ds.excitedStates.map(s => s.Rvel), [0, -0.2552, 0, 1.1111], 1e-9, 'R(velocity)');
  eq(ds.excitedStates[1].label, 'Singlet-B1', 'states unchanged');
});

// ===================================================================== ORCA
test('ORCA convergence, scans and jobs', () => {
  const ds = load('orca5_water_optfreq.out');
  const keys = ['energyChange', 'rmsForce', 'maxForce', 'rmsDisp', 'maxDisp'];
  ok(ds.frames.every(f => f.opt && keys.slice(1).every(q => Number.isFinite(f.opt[q]) && Number.isFinite(f.opt[q + 'Thr']) && typeof f.opt[q + 'Conv'] === 'boolean')), 'every cycle has the gradient/step criteria');
  ok(ds.frames.slice(1).every(f => Number.isFinite(f.opt.energyChange) && Number.isFinite(f.opt.energyChangeThr)), 'energy change from cycle 2');
  eq(ds.jobs.length, 1, 'one job');
  eq(JSON.stringify(ds.jobs[0].kinds), '["opt","freq"]', 'kinds');
  const sc = load('orca5_water_scan.out');
  const pts = sc.frames.filter(f => f.scan);
  eq(pts.length, 2, 'scan points');
  ok(pts.every(f => f.converged && f.scan.coord === 'A(2,1,3)'), '0-based ORCA atoms -> 1-based label');
  arrNear(pts.map(f => f.scan.value), [100, 105], 1e-9, 'scan values');
  eq(JSON.stringify(sc.jobs[0].kinds), '["opt","scan"]', 'scan job kinds');
  const nj = load('review/orca6_newjob_running.out');
  eq(nj.jobs.length, 2, '$new_job -> two jobs');
  eq(nj.jobs[0].label, '! B3LYP def2-SVP Opt', 'job 1 label');
  eq(nj.jobs[1].label, '! PBE0 def2-TZVP', 'job 2 label');
  eq(JSON.stringify(nj.jobs.map(j => j.kinds)), '[["opt"],["sp"]]', 'kinds per job');
  eq(nj.jobs[1].frameStart, nj.frames.length - 1, 'job 2 frame');
  eq(JSON.stringify(nj.jobs.map(j => j.normalTermination)), '[true,null]', 'still-running second job');
});

test('ORCA NMR, CD (ORCA 5 and 6) and IRC summary', () => {
  const n = load('qc/orca5_water_nmr.out');
  arrNear(n.nmr.iso, [325.567, 31.234, 31.234], 1e-9, 'shielding');
  arrNear(n.nmr.aniso, [51.223, 19.567, 19.567], 1e-9, 'anisotropy');
  ok(n.jobs[0].kinds.indexOf('nmr') >= 0, 'nmr kind');
  const c5 = load('qc/orca5_water_cd.out');
  arrNear(c5.excitedStates.slice(0, 3).map(s => s.R), [-12.34567, 0, 3.21], 1e-9, 'ORCA 5 R');
  ok(c5.excitedStates.slice(3).every(s => s.R === undefined), 'triplets get no singlet R');
  const c6 = load('qc/orca6_water_cd.out');
  arrNear(c6.excitedStates.map(s => s.R), [-12.34567, 0, 3.21], 1e-9, 'ORCA 6 R (length)');
  arrNear(c6.excitedStates.map(s => s.Rvel), [-11, 0, 2.9], 1e-9, 'ORCA 6 R (velocity)');
  const irc = load('qc/orca5_hcn_irc.out');
  eq(irc.ircPath.length, 5, 'IRC points');
  eq(irc.ircPath.findIndex(p => p.ts), 2, 'TS marker');
  near(irc.ircPath[0].energy, -93.360140123, 1e-12, 'energy');
  near(irc.ircPath[1].dE, -14.760307, 1e-9, 'dE');
  eq(JSON.stringify(irc.jobs[0].kinds), '["irc"]', 'IRC kind');
});

// ===================================================================== ORCA input
test('ORCA input reader', () => {
  const t = read('qc/orca_water_opt.inp');
  eq(BM.detectFormat('job.inp', t), 'orca-input', 'detect .inp');
  eq(BM.detectFormat('noext', t), 'orca-input', 'detect without extension');
  const ds = BM.parseFile('job.inp', t);
  eq(JSON.stringify(ds.frames[0].Z), '[8,1,1]', 'elements');
  arrNear(Array.from(ds.frames[0].xyz), [0, 0, 0.1173, 0, 0.7572, -0.4692, 0, -0.7572, -0.4692], 1e-9, 'coordinates');
  eq(ds.meta.route, '! B3LYP def2-SVP Opt TightSCF', 'keywords (comments removed)');
  eq(ds.meta.method + '/' + ds.meta.basis, 'B3LYP/def2-SVP', 'method/basis');
  eq(ds.meta.charge + ' ' + ds.meta.multiplicity, '0 1', 'charge/multiplicity');
  ok(ds.meta.jobTypes.indexOf('opt') >= 0, 'opt job');
  const zi = load('qc/orca_water_int.inp');
  const x = zi.frames[0].xyz;
  near(dist(x, 0, 1), 0.959476, 1e-6, '* int bond');
  near(dist(x, 0, 2), 0.959476, 1e-6, '* int bond 2');
  near(Math.acos(((x[3] - x[0]) * (x[6] - x[0]) + (x[4] - x[1]) * (x[7] - x[1]) + (x[5] - x[2]) * (x[8] - x[2])) / 0.959476 ** 2) * 180 / Math.PI, 104.25, 1e-4, '* int angle');
  const nj = load('qc/orca_newjob_xyzfile.inp');
  eq(nj.frames.length, 1, 'first job geometry');
  eq(nj.meta.charge + ' ' + nj.meta.multiplicity, '-1 2', 'charge/multiplicity');
  ok(nj.warnings.some(w => /orca_newjob_xyzfile\.xyz/.test(w) && /xyzfile/.test(w)), '* xyzfile warning');
  let msg = '';
  try { load('qc/orca_xyzfile_only.inp'); } catch (e) { msg = e.message; }
  ok(/water\.xyz/.test(msg), '* xyzfile only -> error naming the file: ' + msg);
  ok(BM.formats.some(f => f.id === 'orca-input' && f.extensions.indexOf('inp') >= 0), 'listed in BM.formats');
  ok(typeof BM.parsers.orcaInput === 'function', 'BM.parsers.orcaInput');
  // not claimed: CP2K / GAMESS inputs and Gaussian inputs with "!" comments
  ok(BM.detectFormat('a.inp', read('qc/cp2k_water.inp')) !== 'orca-input', 'CP2K input not ORCA');
  ok(BM.detectFormat('a.inp', read('qc/gamess_water.inp')) !== 'orca-input', 'GAMESS input not ORCA');
  eq(BM.detectFormat('a.gjf', '! comment\n#p opt b3lyp/6-31g(d)\n\ntitle\n\n0 1\nH 0 0 0\nH 0 0 0.74\n\n'), 'gaussian-input', 'Gaussian input with ! comment');
  eq(BM.detectFormat('a', '# comment\n! HF def2-SVP\n* xyz 0 1\nH 0 0 0\nH 0 0 0.74\n*\n'), 'orca-input', 'ORCA input starting with # comment');
});

// ===================================================================== save formats
test('Molden export round trip (atoms + vibrations)', () => {
  const src = load('g16_water_optfreq.log');
  const fr = src.frames[src.frames.length - 1];
  const ctx = { Z: fr.Z, xyz: fr.xyz, title: 'water', frames: src.frames, frameIndex: src.frames.length - 1, vibrations: src.vibrations, meta: src.meta };
  const fmt = BM.saveFormats.find(f => f.id === 'molden');
  const text = fmt.write(ctx);
  eq(BM.detectFormat('w.molden', text), 'molden', 'detected');
  const back = BM.parseFile('w.molden', text);
  eq(JSON.stringify(back.frames[0].Z), JSON.stringify(fr.Z), 'atoms');
  arrNear(Array.from(back.frames[0].xyz), Array.from(fr.xyz), 1e-9, 'coordinates');
  eq(back.frames.length, 1, 'FR-COORD matches [Atoms]');
  arrNear(back.vibrations.freqs, src.vibrations.freqs, 1e-4, 'frequencies');
  arrNear(back.vibrations.ir, src.vibrations.ir, 1e-4, 'IR intensities');
  arrNear(Array.from(back.vibrations.modes[2]), Array.from(src.vibrations.modes[2]), 1e-8, 'mode 3');
  const plain = fmt.write({ Z: [8, 1], xyz: [0, 0, 0, 0, 0, 0.97], title: 'OH' });
  ok(plain.indexOf('[FREQ]') < 0 && BM.parseFile('oh.molden', plain).frames[0].Z.length === 2, 'atoms only without vibrations');
});

test('cube export round trip', () => {
  const src = load('water_homo_bohr.cube');
  const v = src.volumes[0];
  const fmt = BM.saveFormats.find(f => f.id === 'cube');
  ok(fmt.needs.indexOf('volume') >= 0 && fmt.kind === 'volume', 'needs a volume');
  const text = fmt.write({ Z: src.frames[0].Z, xyz: src.frames[0].xyz, title: 'w', volume: v });
  const back = BM.parseFile('w.cube', text);
  const b = back.volumes[0];
  eq(JSON.stringify(b.dims), JSON.stringify(v.dims), 'dims');
  arrNear(b.origin, v.origin, 1e-6, 'origin (Angstrom)');
  arrNear([].concat(...b.axes), [].concat(...v.axes), 1e-6, 'axes');
  arrNear(Array.from(back.frames[0].xyz), Array.from(src.frames[0].xyz), 1e-6, 'atoms');
  let worst = 0;
  for (let i = 0; i < v.data.length; i++) worst = Math.max(worst, Math.abs(b.data[i] - v.data[i]) / Math.max(1e-30, Math.abs(v.data[i])));
  ok(worst < 1e-5, 'data to 6 significant digits (worst relative error ' + worst + ')');
  eq(b.mo, v.mo, 'MO index kept');
  const two = load('water_mo5_mo6.cube');
  const t2 = fmt.write({ Z: two.frames[0].Z, xyz: two.frames[0].xyz, volume: two.volumes[1] });
  const b2 = BM.parseFile('mo.cube', t2).volumes[0];
  eq(b2.mo, two.volumes[1].mo, 'MO cube: orbital number');
  let threw = false;
  try { fmt.write({ Z: [1], xyz: [0, 0, 0] }); } catch (e) { threw = /volume/i.test(e.message); }
  ok(threw, 'error without volume data');
});

test('Gaussian / ORCA input export', () => {
  const g = load('g16_water_optfreq.log');
  const f = g.frames[3];
  const gi = BM.saveFormats.find(x => x.id === 'gaussian-input');
  const t = gi.write({ Z: f.Z, xyz: f.xyz, title: 'water', meta: g.meta });
  const back = BM.parseFile('w.gjf', t);
  eq(back.meta.route, g.meta.route, 'route from meta.route');
  arrNear(Array.from(back.frames[0].xyz), Array.from(f.xyz), 1e-7, 'coordinates');
  const meta2 = Object.assign({}, g.meta, { route: '#p opt b3lyp/6-31g(d) geom=connectivity guess=read' });
  eq(BM.parseFile('w.gjf', gi.write({ Z: f.Z, xyz: f.xyz, meta: meta2 })).meta.route, '#p opt b3lyp/6-31g(d)', 'geom=/guess= dropped');
  const o = load('orca5_water_optfreq.out');
  eq(BM.parseFile('w.gjf', gi.write({ Z: o.frames[0].Z, xyz: o.frames[0].xyz, meta: o.meta })).meta.route, '#p', 'non-Gaussian route -> #p');
  const oi = BM.saveFormats.find(x => x.id === 'orca-input');
  const ot = oi.write({ Z: o.frames[0].Z, xyz: o.frames[0].xyz, meta: o.meta });
  const ob = BM.parseFile('w.inp', ot);
  eq(ob.format, 'orca-input', 'ORCA input read back');
  eq(ob.meta.route, o.meta.route, 'ORCA keywords kept');
  arrNear(Array.from(ob.frames[0].xyz), Array.from(o.frames[0].xyz), 1e-7, 'coordinates');
  const fromG = BM.parseFile('w.inp', oi.write({ Z: f.Z, xyz: f.xyz, meta: g.meta }));
  eq(fromG.meta.route, '! B3LYP 6-31G(d)', 'Gaussian method/basis -> ORCA keywords');
});

// ===================================================================== detection unchanged
test('detection of existing fixtures unchanged with the extra readers loaded', () => {
  require(path.join(ROOT, 'src', 'readers-qc.js'));
  ok(BM.formatRegistry.length >= 10, 'extra readers registered');
  const after = ALL.map(f => BM.detectFormat(path.basename(f), read(f)));
  const diff = ALL.filter((f, k) => after[k] !== before[k]).map((f, k) => f);
  eq(diff.length, 0, 'changed detections: ' + diff.join(', '));
  ok(before.filter(s => s === 'gaussian-log').length >= 8 && before.filter(s => s === 'orca-out').length >= 6, 'Gaussian/ORCA fixtures counted');
  for (const f of ALL) if (/\.(log|out)$/.test(f)) { const s = BM.detectFormat(path.basename(f).replace(/\.\w+$/, '.out'), read(f)); ok(s === 'gaussian-log' || s === 'orca-out', f + ' as .out -> ' + s); }
});

console.log('\nBlueMolecule parser additions tests: ' + pass + ' passed, ' + fail + ' failed');
if (fail) { for (const f of failures) console.log('  FAIL ' + f); process.exit(1); }
