#!/usr/bin/env node
/* BlueMolecule input-generator tests.  Run:  node tests/test_writers.js
 * Checks every program of BM.inputPrograms: required sections, charge/multiplicity, task keywords,
 * full-precision coordinates (bohr for Turbomole) and read-back of Gaussian, ORCA and xtb files. */
'use strict';
const path = require('path');
const ROOT = path.join(__dirname, '..');
require(path.join(ROOT, 'src', 'elements.js'));
require(path.join(ROOT, 'src', 'parsers.js'));
require(path.join(ROOT, 'src', 'writers.js'));
const BM = globalThis.BM;
const W = BM.writers;

let pass = 0, fail = 0, current = '';
const failures = [], notes = [];
function ok(cond, msg) { if (cond) pass++; else { fail++; failures.push('[' + current + '] ' + msg); } }
function eq(a, b, msg) { ok(a === b, msg + ': got ' + JSON.stringify(a) + ', expected ' + JSON.stringify(b)); }
function has(text, s, msg) { ok(typeof text === 'string' && (s instanceof RegExp ? s.test(text) : text.indexOf(s) >= 0), (msg || 'contains') + ' ' + s); }
function hasNot(text, s, msg) { ok(typeof text === 'string' && !(s instanceof RegExp ? s.test(text) : text.indexOf(s) >= 0), (msg || 'does not contain') + ' ' + s); }
function before(text, a, b, msg) { const i = text.indexOf(a), k = text.indexOf(b); ok(i >= 0 && k > i, (msg || 'order') + ': "' + a + '" before "' + b + '"'); }
function test(name, fn) {
  current = name;
  try { fn(); } catch (e) { fail++; failures.push('[' + name + '] threw: ' + (e && e.stack || e)); }
}
function throwsMsg(fn, re, msg) {
  try { fn(); ok(false, msg + ': no error thrown'); } catch (e) { ok(e instanceof Error && re.test(e.message), msg + ': message "' + e.message + '"'); }
}
function arrNear(a, b, tol, msg) {
  if (!a || a.length !== b.length) { ok(false, msg + ': length ' + (a && a.length) + ' vs ' + b.length); return; }
  for (let i = 0; i < b.length; i++) if (!(Math.abs(a[i] - b[i]) <= tol)) { ok(false, msg + '[' + i + ']: got ' + a[i] + ', expected ' + b[i]); return; }
  pass++;
}
const BOHR = 0.529177210903;
const num = s => Number(String(s).replace(/[dD]/, 'e'));

// Coordinates with 10 significant decimals: a full-precision writer must reproduce them to 5e-11.
const WATER = { Z: [8, 1, 1], xyz: [0.0000000000, 0.0000000000, 0.1173456789, 0.0000000000, 0.7572345678, -0.4692123456,
  0.0000000000, -0.7572345678, -0.4692123456] };
const PT = { Z: [78, 17, 17, 7, 1, 1, 1, 7, 1, 1, 1], xyz: [
  0.0123456789, -0.0234567891, 0.0012345678, 2.3123456789, 0.0011223344, -0.0123456789, 0.0011223344, 2.3198765432, 0.0101010101,
  -2.0512345678, 0.0021212121, 0.0123123123, -2.4012345678, 0.9512345678, 0.0222222222, -2.4098765432, -0.4812345678, 0.8312345678,
  -2.4087654321, -0.4723456789, -0.8156789012, 0.0034343434, -2.0598765432, -0.0111111111, 0.9612345678, -2.4012345678, 0.0198765432,
  -0.4834567891, -2.4123456789, 0.8234567891, -0.4767891234, -2.4156789123, -0.8212345678] };
const base = extra => Object.assign({ Z: WATER.Z, xyz: WATER.xyz, charge: 0, multiplicity: 1, task: 'opt', nprocs: 4, memoryMB: 4000,
  title: 'water test' }, extra || {});
const prog = id => BM.inputPrograms.find(p => p.id === id);

/** Coordinates from atom lines "Sym x y z" between two marker lines (exclusive). */
function coordBlock(text, startRe, endRe, cols) {
  const L = text.split('\n');
  let i = L.findIndex(s => startRe.test(s));
  const out = [];
  if (i < 0) return out;
  for (i++; i < L.length && !endRe.test(L[i]); i++) {
    const t = L[i].trim().split(/\s+/);
    if (t.length >= 4 && /^[A-Za-z]{1,2}$/.test(t[0]) && cols.every(c => /^-?\d+\.\d+$/.test(t[c]))) out.push(...cols.map(c => num(t[c])));
  }
  return out;
}

// ------------------------------------------------------------------ registry
test('program registry', () => {
  const ids = BM.inputPrograms.map(p => p.id);
  eq(JSON.stringify(ids), JSON.stringify(['gaussian', 'orca', 'psi4', 'qchem', 'nwchem', 'gamess', 'mopac', 'xtb', 'molpro', 'turbomole']), 'ids');
  for (const p of BM.inputPrograms) {
    ok(typeof p.name === 'string' && p.name && typeof p.ext === 'string', p.id + ' name/ext');
    ok(Array.isArray(p.methods) && Array.isArray(p.bases) && Array.isArray(p.tasks) && typeof p.write === 'function', p.id + ' lists and write()');
    ok(p.tasks.join() === 'sp,opt,freq,optfreq', p.id + ' tasks');
    ok(p.defaultMethod === null || p.methods.includes(p.defaultMethod), p.id + ' default method is listed');
    ok(p.defaultBasis === null || p.bases.includes(p.defaultBasis), p.id + ' default basis is listed');
    for (const task of p.tasks) {
      for (const mol of [WATER, PT]) {
        const r = p.write(Object.assign({}, mol, { task, charge: 0, multiplicity: 1, title: p.id + ' ' + task }));
        ok(r && typeof r.text === 'string' && r.text.length > 20 && /\n$/.test(r.text), p.id + ' ' + task + ' text');
        ok(r.note === undefined || typeof r.note === 'string', p.id + ' note type');
      }
    }
  }
  eq(prog('mopac').defaultBasis, null, 'MOPAC has no basis'); eq(prog('xtb').defaultBasis, null, 'xtb has no basis');
  ok(BM.inputPrograms.every(p => p.id === 'turbomole' || typeof p.write(base()).command === 'string'), 'run commands');
  // Notes: impossible spin state, dummies left out.
  has(prog('psi4').write(base({ multiplicity: 2 })).note, 'cannot have multiplicity 2', 'parity warning');
  ok(prog('psi4').write(base({ multiplicity: 3 })).note === undefined, 'triplet water is fine');
  has(prog('qchem').write(base({ Z: [8, 1, 1, 0], xyz: WATER.xyz.concat([0, 0, 2]) })).note, '1 dummy/ghost atom(s) were left out', 'dummy note');
  throwsMsg(() => prog('orca').write(base({ task: 'irc' })), /Unknown task/, 'unknown task');
  throwsMsg(() => W.psi4(base({ charge: 0.5 })), /charge must be an integer/, 'fractional charge');
  throwsMsg(() => W.nwchem(base({ multiplicity: 0 })), /multiplicity must be 1 or larger/, 'zero multiplicity');
  throwsMsg(() => W.qchem({ Z: [0], xyz: [0, 0, 0] }), /no real atoms/, 'only dummies');
  eq(W.psi4(base({ task: 'Opt+Freq' })), W.psi4(base({ task: 'optfreq' })), 'task aliases');
  ok(typeof W.xyz === 'function' && typeof W.gaussianInput === 'function' && typeof W.orcaInput === 'function' &&
    typeof W.mol === 'function' && typeof W.pdb === 'function' && typeof W.xyzMulti === 'function', 'existing writers kept');
});

// ------------------------------------------------------------------ Gaussian
test('Gaussian', () => {
  const r = prog('gaussian').write(base({ task: 'optfreq', solvent: 'water' }));
  has(r.text, /^%nprocshared=4\n%mem=4000MB\n#p B3LYP\/6-31G\(d\) opt freq scrf=\(smd,solvent=water\)\n\nwater test\n\n0 1\nO /, 'link0, route, title, charge line');
  ok(/\n\n\n$/.test(r.text), 'blank-line terminated');
  const ds = BM.parseFile('job.gjf', r.text);
  arrNear(Array.from(ds.frames[0].xyz), WATER.xyz, 1e-8, 'read back coordinates');
  eq(JSON.stringify(ds.frames[0].Z), '[8,1,1]', 'read back elements');
  eq(ds.meta.charge, 0, 'charge'); eq(ds.meta.multiplicity, 1, 'multiplicity');
  ok(ds.meta.jobTypes.includes('opt') && ds.meta.jobTypes.includes('freq'), 'job types ' + ds.meta.jobTypes);
  eq(ds.meta.solvent, 'water', 'solvent'); eq(ds.meta.method, 'B3LYP', 'method'); eq(ds.meta.basis, '6-31G(d)', 'basis');
  const sp = BM.parseFile('a.gjf', prog('gaussian').write(base({ task: 'sp', nprocs: 1, memoryMB: 0 })).text);
  eq(JSON.stringify(sp.meta.jobTypes), '["sp"]', 'sp has no task keyword'); eq(sp.meta.link0.length, 0, 'no link0 by default');
  const fr = prog('gaussian').write(base({ task: 'freq', method: 'PBE1PBE', basis: 'def2-SVP' })).text;
  has(fr, '#p PBE1PBE/Def2SVP freq', 'def2 basis in Gaussian spelling');
  const pt = prog('gaussian').write(Object.assign({}, PT, { task: 'opt', charge: 1, multiplicity: 2, method: 'B3LYP EmpiricalDispersion=GD3BJ',
    basis: 'GenECP', extra: 'Pt 0\nLANL2DZ\n****\nH N Cl 0\n6-31G(d)\n****\n\nPt 0\nLANL2DZ' }));
  has(pt.text, '#p B3LYP/GenECP EmpiricalDispersion=GD3BJ opt', 'extra route words after method/basis');
  has(pt.text, /\n1 2\nPt /, 'charge/multiplicity'); has(pt.text, /\n\nPt 0\nLANL2DZ\n\*\*\*\*\n/, 'extra input after the geometry');
  const pds = BM.parseFile('pt.gjf', pt.text);
  arrNear(Array.from(pds.frames[0].xyz), PT.xyz, 1e-8, 'Pt complex coordinates'); eq(pds.meta.multiplicity, 2, 'doublet read back');
  has(prog('gaussian').write(base({ method: 'PM7', basis: 'def2-SVP' })).text, '#p PM7 opt\n', 'semi-empirical: no basis');
});

// ------------------------------------------------------------------ ORCA
test('ORCA', () => {
  const r = prog('orca').write(base({ task: 'optfreq', solvent: 'water', extra: '%scf\n  maxiter 300\nend' }));
  has(r.text, /^! B3LYP def2-SVP Opt Freq CPCM\(water\)\n# water test\n%pal nprocs 4 end\n%maxcore 1000\n%scf\n  maxiter 300\nend\n\n\* xyz 0 1\n/, 'keywords, title comment, %pal, %maxcore per core, blocks');
  ok(/\n\*\n$/.test(r.text), 'geometry closed with *');
  ok(BM.detectFormat('job.inp', r.text) !== 'gaussian-input', 'not mistaken for a Gaussian input');
  const keywords = t => t.split('\n').filter(s => s[0] === '!').join(' ');
  eq(keywords(prog('orca').write(base({ task: 'sp' })).text), '! B3LYP def2-SVP', 'sp keywords');
  eq(keywords(prog('orca').write(base({ task: 'opt' })).text), '! B3LYP def2-SVP Opt', 'opt keywords');
  eq(keywords(prog('orca').write(base({ task: 'freq' })).text), '! B3LYP def2-SVP Freq', 'freq keywords');
  eq(keywords(prog('orca').write(base({ task: 'opt', method: 'r2SCAN-3c', basis: 'def2-TZVP' })).text), '! r2SCAN-3c Opt', 'composite method: no basis');
  has(prog('orca').write(base({ charge: -1, multiplicity: 2 })).text, '* xyz -1 2', 'charge/multiplicity');
  hasNot(prog('orca').write(base({ nprocs: 1, memoryMB: 0 })).text, '%', 'serial job without memory: no % blocks');
  // Read back: through BM.parseFile when it knows ORCA inputs, otherwise with a minimal "* xyz" reader.
  const text = prog('orca').write(Object.assign({}, PT, { task: 'opt', nprocs: 2 })).text;
  let xyz = null;
  if (BM.detectFormat('job.inp', text)) {
    const ds = BM.parseFile('job.inp', text);
    xyz = Array.from(ds.frames[0].xyz);
    eq(JSON.stringify(ds.frames[0].Z), JSON.stringify(PT.Z), 'read back elements (BM.parseFile)');
  } else {
    notes.push('BM.parseFile has no ORCA-input reader yet; the ORCA geometry was read back with the test\'s own "* xyz" reader.');
    xyz = coordBlock(text, /^\* xyz /, /^\*$/, [1, 2, 3]);
  }
  arrNear(xyz, PT.xyz, 1e-8, 'read back coordinates');
});

// ------------------------------------------------------------------ Psi4
test('Psi4', () => {
  const t = W.psi4(base({ task: 'optfreq' }));
  has(t, /^# water test\nmemory 4000 MB\nset_num_threads\(4\)\n\nmolecule \{\n0 1\nO /, 'header and molecule block');
  has(t, 'units angstrom\n}', 'units'); has(t, '\nset basis def2-SVP\n', 'set basis');
  before(t, "optimize('b3lyp')", "frequency('b3lyp')", 'opt then freq');
  hasNot(t, 'energy(', 'no separate energy call');
  arrNear(coordBlock(t, /^molecule \{/, /^units/, [1, 2, 3]), WATER.xyz, 5e-11, 'coordinates (10 decimals)');
  has(W.psi4(base({ task: 'sp' })), "\nenergy('b3lyp')\n", 'sp'); hasNot(W.psi4(base({ task: 'sp' })), 'optimize(', 'sp has no optimize');
  const fq = W.psi4(base({ task: 'freq' }));
  has(fq, "\nfrequency('b3lyp')\n", 'freq'); hasNot(fq, 'optimize(', 'freq has no optimize');
  const open = W.psi4(base({ charge: 1, multiplicity: 2 }));
  has(open, 'molecule {\n1 2\n', 'cation doublet'); has(open, 'set reference uks', 'UKS for DFT');
  has(W.psi4(base({ method: 'HF', multiplicity: 3 })), 'set reference uhf', 'UHF for HF');
  const sol = W.psi4(base({ solvent: 'dmso', extra: 'set scf_type df' }));
  has(sol, 'set pcm true', 'PCM on'); has(sol, 'Solvent = Dimethylsulfoxide', 'PCMSolver solvent name');
  before(sol, 'set scf_type df', "optimize('b3lyp')", 'extra before the call');
  hasNot(W.psi4(base({ nprocs: 1, memoryMB: 0 })), /memory|set_num_threads/, 'no memory/threads by default');
});

// ------------------------------------------------------------------ Q-Chem
test('Q-Chem', () => {
  const t = W.qchem(base({ task: 'opt', solvent: 'water' }));
  has(t, /^\$comment\nwater test\n\$end\n\n\$molecule\n0 1\n   O /, '$comment and $molecule');
  has(t, /\$rem\n   JOBTYPE +opt\n   METHOD +b3lyp\n   BASIS +def2-SVP\n   MEM_TOTAL +4000\n   SOLVENT_METHOD +PCM\n\$end/, '$rem');
  has(t, '$solvent\n   SolventName      water\n$end', '$solvent');
  arrNear(coordBlock(t, /^\$molecule/, /^\$end/, [1, 2, 3]), WATER.xyz, 5e-11, 'coordinates');
  for (const [task, jt] of [['sp', 'sp'], ['freq', 'freq']]) {
    const s = W.qchem(base({ task }));
    has(s, new RegExp('JOBTYPE +' + jt + '\\n'), task + ' JOBTYPE'); hasNot(s, '@@@', task + ' single job');
  }
  const of = W.qchem(base({ task: 'optfreq', multiplicity: 3, extra: 'SCF_CONVERGENCE 8' }));
  const jobs = of.split('\n@@@\n');
  eq(jobs.length, 2, 'opt+freq is two jobs joined by @@@');
  has(jobs[0], /JOBTYPE +opt/, 'first job optimises'); has(jobs[0], '$molecule\n0 3\n', 'triplet');
  has(jobs[1], '$molecule\n   read\n$end', 'second job reads the geometry'); has(jobs[1], /JOBTYPE +freq/, 'second job frequencies');
  has(jobs[1], /SCF_GUESS +read/, 'reads the orbitals'); has(jobs[0], '   SCF_CONVERGENCE 8\n$end', 'extra $rem line');
  has(jobs[1], '   SCF_CONVERGENCE 8\n$end', 'extra $rem line in both jobs');
  const sec = W.qchem(base({ extra: '$pcm\n   Theory CPCM\n$end' }));
  has(sec, '\n$pcm\n   Theory CPCM\n$end\n', 'extra section appended'); hasNot(sec, '   $pcm', 'section not put inside $rem');
});

// ------------------------------------------------------------------ NWChem
test('NWChem', () => {
  const t = W.nwchem(base({ task: 'optfreq', baseName: 'water opt' }));
  has(t, /^echo\nstart water_opt\ntitle "water test"\nmemory total 1000 mb\ncharge 0\n\ngeometry units angstrom\n/, 'header');
  has(t, 'basis\n  * library def2-svp\nend', 'basis block'); has(t, 'dft\n  xc b3lyp\n  mult 1\nend', 'dft block');
  before(t, 'task dft optimize', 'task dft freq', 'optimize then freq');
  arrNear(coordBlock(t, /^geometry/, /^end/, [1, 2, 3]), WATER.xyz, 5e-11, 'coordinates');
  has(W.nwchem(base({ task: 'sp' })), /\ntask dft energy\n$/, 'sp task'); has(W.nwchem(base({ task: 'freq' })), /\ntask dft freq\n$/, 'freq task');
  const open = W.nwchem(base({ charge: -1, multiplicity: 2, basis: '6-31+G(d,p)' }));
  has(open, 'charge -1', 'charge'); has(open, '  mult 2\n  odft', 'open-shell DFT'); has(open, '* library 6-31+g**', 'Pople star notation');
  const hf = W.nwchem(base({ method: 'HF', task: 'sp', multiplicity: 3 }));
  has(hf, 'scf\n  uhf\n  nopen 2\nend', 'UHF triplet'); has(hf, 'task scf energy', 'scf task');
  has(W.nwchem(base({ method: 'MP2' })), 'task mp2 optimize', 'MP2'); has(W.nwchem(base({ method: 'CCSD(T)', task: 'sp' })), 'task ccsd(t) energy', 'CCSD(T)');
  has(W.nwchem(base({ method: 'PBE' })), 'xc xpbe96 cpbe96', 'PBE functional');
  has(W.nwchem(base({ method: 'B3LYP-D3BJ' })), '  xc b3lyp\n  disp vdw 4', 'D3(BJ) dispersion');
  const sol = W.nwchem(base({ solvent: 'water', extra: 'driver\n  maxiter 50\nend' }));
  has(sol, 'cosmo\n  solvent water\nend', 'COSMO'); before(sol, 'driver', 'task dft optimize', 'extra before task');
});

// ------------------------------------------------------------------ GAMESS
/** Keywords of a GAMESS group, joined across continuation lines. */
function gamessGroup(text, name) {
  const m = new RegExp('(^|\\n) \\$' + name + '\\b([\\s\\S]*?)\\$END').exec(text);
  return m ? m[2].trim().split(/\s+/).join(' ') : null;
}

test('GAMESS (US)', () => {
  const t = W.gamess(base({ task: 'optfreq', solvent: 'water' }));
  ok(/^ \$CONTRL /.test(t), 'starts with $CONTRL');
  eq(gamessGroup(t, 'CONTRL'), 'SCFTYP=RHF RUNTYP=OPTIMIZE ICHARG=0 MULT=1 DFTTYP=B3LYP MAXIT=100', '$CONTRL');
  eq(gamessGroup(t, 'SYSTEM'), 'MWORDS=131', '$SYSTEM (4000 MB over 4 processes)');
  eq(gamessGroup(t, 'BASIS'), 'GBASIS=N31 NGAUSS=6 NDFUNC=1', '6-31G(d)');
  eq(gamessGroup(t, 'STATPT'), 'NSTEP=100 HSSEND=.TRUE.', 'Hessian after the optimisation');
  eq(gamessGroup(t, 'FORCE'), 'METHOD=SEMINUM', 'semi-numerical DFT Hessian'); eq(gamessGroup(t, 'PCM'), 'SOLVNT=WATER', '$PCM');
  has(t, /\n \$DATA\nwater test\nC1\nO +8\.0 /, '$DATA: title, C1, nuclear charge');
  ok(/\n \$END\n$/.test(t), '$DATA closed');
  ok(t.split('\n').every(s => s.length <= 80), 'lines within 80 columns');
  ok(t.split('\n').filter(s => s.indexOf('$') >= 0 && !/^ \$/.test(s) && !/^  /.test(s)).length === 0, 'groups start in column 2');
  const L = t.split('\n'), k = L.indexOf('C1'), xyz = [];
  for (let i = k + 1; i < L.length && !/\$END/.test(L[i]); i++) { const s = L[i].trim().split(/\s+/); xyz.push(num(s[2]), num(s[3]), num(s[4])); }
  arrNear(xyz, WATER.xyz, 5e-11, 'coordinates');
  has(W.gamess(base({ task: 'sp' })), 'RUNTYP=ENERGY', 'sp'); has(W.gamess(base({ task: 'freq' })), 'RUNTYP=HESSIAN', 'freq');
  hasNot(W.gamess(base({ task: 'opt' })), 'HSSEND', 'plain opt');
  const open = W.gamess(base({ charge: 1, multiplicity: 2, method: 'MP2', basis: 'aug-cc-pVTZ' }));
  has(open, 'SCFTYP=UHF', 'UHF'); has(open, 'ICHARG=1 MULT=2', 'charge/multiplicity'); has(open, 'MPLEVL=2', 'MP2');
  has(open, 'GBASIS=ACCT', 'aug-cc-pVTZ'); hasNot(open, 'DFTTYP', 'no DFT for MP2');
  eq(gamessGroup(W.gamess(base({ basis: '6-311++G(2d,p)' })), 'BASIS'), 'GBASIS=N311 NGAUSS=6 NDFUNC=2 NPFUNC=1 DIFFSP=.TRUE. DIFFS=.TRUE.', '6-311++G(2d,p)');
  eq(gamessGroup(W.gamess(base({ basis: '6-31G**' })), 'BASIS'), 'GBASIS=N31 NGAUSS=6 NDFUNC=1 NPFUNC=1', '6-31G**');
  eq(gamessGroup(W.gamess(base({ basis: 'STO-3G' })), 'BASIS'), 'GBASIS=STO NGAUSS=3', 'STO-3G');
  eq(gamessGroup(W.gamess(base({ basis: 'cc-pCVTZ' })), 'BASIS'), 'GBASIS=CCTC', 'cc-pCVTZ');
  eq(gamessGroup(W.gamess(base({ method: 'PM3' })), 'BASIS'), 'GBASIS=PM3', 'semi-empirical via $BASIS');
  eq(gamessGroup(W.gamess(base({ method: 'HF', task: 'freq' })), 'FORCE'), null, 'analytic RHF Hessian: no $FORCE');
  const many = W.gamess(base({ method: 'CAM-B3LYP', multiplicity: 3, charge: -2, task: 'optfreq', basis: '6-311++G(2df,2pd)' }));
  ok(many.split('\n').every(s => s.length <= 79), 'continued groups stay within 79 columns');
  eq(gamessGroup(many, 'CONTRL'), 'SCFTYP=UHF RUNTYP=OPTIMIZE ICHARG=-2 MULT=3 DFTTYP=CAMB3LYP MAXIT=100', 'continued $CONTRL');
  has(W.gamess(Object.assign({}, PT, { task: 'sp' })), /\nPt +78\.0 /, 'Pt nuclear charge');
  has(prog('gamess').write(base({ basis: 'def2-SVP' })).note, 'GBASIS=DEF2-SVP', 'unknown basis keyword is flagged');
});

// ------------------------------------------------------------------ MOPAC
test('MOPAC', () => {
  const t = W.mopac(base({ task: 'opt', method: undefined, basis: 'def2-SVP' }));
  const L = t.split('\n');
  eq(L[0], 'PM7 CHARGE=0 SINGLET THREADS=4', 'keyword line (basis ignored)');
  eq(L[1], 'water test', 'first title line'); ok(/^BlueMolecule: /.test(L[2]), 'second title line');
  has(L[3], /^O +0\.0000000000 1 +0\.0000000000 1 +0\.1173456789 1$/, 'Cartesian line with optimisation flags');
  const xyz = [];
  for (const s of L.slice(3).filter(Boolean)) { const q = s.trim().split(/\s+/); xyz.push(num(q[1]), num(q[3]), num(q[5])); }
  arrNear(xyz, WATER.xyz, 5e-11, 'coordinates');
  eq(W.mopac(base({ task: 'sp', nprocs: 1 })).split('\n')[0], 'PM7 CHARGE=0 SINGLET 1SCF', 'sp keywords');
  ok(/ 0$/.test(W.mopac(base({ task: 'sp' })).split('\n')[3]), 'sp: flags 0');
  eq(W.mopac(base({ task: 'freq', nprocs: 1 })).split('\n')[0], 'PM7 CHARGE=0 SINGLET FORCE', 'freq keywords');
  eq(W.mopac(base({ task: 'optfreq', nprocs: 1 })).split('\n')[0], 'PM7 CHARGE=0 SINGLET OPT FORCE', 'opt+freq keywords');
  eq(W.mopac(base({ charge: 1, multiplicity: 2, nprocs: 1, method: 'PM6-D3H4', solvent: 'water' })).split('\n')[0],
    'PM6-D3H4 CHARGE=1 DOUBLET EPS=78.36', 'cation doublet, COSMO dielectric');
  eq(W.mopac(base({ multiplicity: 3, nprocs: 1, solvent: '36.6' })).split('\n')[0], 'PM7 CHARGE=0 TRIPLET EPS=36.6', 'numeric dielectric');
  throwsMsg(() => W.mopac(base({ solvent: 'unobtainium' })), /dielectric constant/, 'unknown solvent');
  const long = W.mopac(base({ extra: 'GNORM=0.01 PRECISE LARGE AUX(0,PRECISION=9) GEO-OK LET DISP SUPER BONDS LOCALIZE', nprocs: 1 })).split('\n');
  ok(/ \+$/.test(long[0]) && long[0].length <= 80 && !/\+$/.test(long[1]) && long[2] === 'water test', 'long keyword line continued with +');
});

// ------------------------------------------------------------------ xtb
test('xtb', () => {
  const r = W.xtb(base({ task: 'opt', nprocs: 1 }));
  ok(Array.isArray(r.files) && r.files.length === 1 && r.files[0].name === 'mol.xyz', 'one mol.xyz file');
  eq(r.command, 'xtb mol.xyz --opt --chrg 0 --uhf 0 --gfn 2', 'opt command');
  const ds = BM.parseFile(r.files[0].name, r.files[0].text);
  arrNear(Array.from(ds.frames[0].xyz), WATER.xyz, 5e-11, 'XYZ reads back with BM.parseFile');
  eq(JSON.stringify(ds.frames[0].Z), '[8,1,1]', 'elements'); eq(ds.frames[0].label, 'water test', 'title as comment');
  eq(W.xtb(base({ task: 'sp', nprocs: 1 })).command, 'xtb mol.xyz --chrg 0 --uhf 0 --gfn 2', 'sp');
  eq(W.xtb(base({ task: 'freq', nprocs: 1 })).command, 'xtb mol.xyz --hess --chrg 0 --uhf 0 --gfn 2', 'freq');
  eq(W.xtb(base({ task: 'optfreq', nprocs: 1 })).command, 'xtb mol.xyz --ohess --chrg 0 --uhf 0 --gfn 2', 'opt+freq');
  eq(W.xtb(base({ charge: -1, multiplicity: 2, solvent: 'Dichloromethane', method: 'GFN1-xTB', baseName: 'anion.xyz' })).command,
    'xtb anion.xyz --opt --chrg -1 --uhf 1 --gfn 1 --alpb ch2cl2 -P 4', 'anion doublet, ALPB, threads');
  eq(W.xtb(base({ method: 'GFN-FF', nprocs: 1 })).command, 'xtb mol.xyz --opt --chrg 0 --uhf 0 --gfnff', 'GFN-FF');
  const xc = W.xtb(base({ nprocs: 1, extra: '$constrain\n  force constant=1.0\n  distance: 1, 2, auto' }));
  eq(xc.files.length, 2, 'xcontrol file'); eq(xc.files[1].name, 'xcontrol', 'xcontrol name');
  ok(/\$end\n$/.test(xc.files[1].text), 'xcontrol terminated'); ok(/--input xcontrol$/.test(xc.command), 'xcontrol on the command line');
  eq(W.xtb(base({ nprocs: 1, extra: '--etemp 500' })).command, 'xtb mol.xyz --opt --chrg 0 --uhf 0 --gfn 2 --etemp 500', 'extra flags');
  throwsMsg(() => W.xtb(base({ method: 'B3LYP' })), /GFN2-xTB/, 'non-xtb method');
  const reg = prog('xtb').write(Object.assign({}, PT, { task: 'optfreq', multiplicity: 1, nprocs: 1 }));
  eq(reg.text, reg.files[0].text, 'registry text = xyz file'); has(reg.command, '--ohess', 'registry command');
  arrNear(Array.from(BM.parseFile('mol.xyz', reg.text).frames[0].xyz), PT.xyz, 5e-11, 'Pt complex reads back');
});

// ------------------------------------------------------------------ Molpro
test('Molpro', () => {
  const t = W.molpro(base({ task: 'optfreq' }));
  has(t, /^\*\*\*,water test\nmemory,125,m {3}! megawords per process\ngeomtyp=xyz\ngeometry=\{\n3\nwater test\nO /, 'header and geometry block');
  has(t, '}\nbasis=cc-pVDZ\nset,charge=0\nset,spin=0\n{ks,b3lyp}\noptg\n{frequencies}\n---\n', 'basis, charge/spin, method, optg, frequencies');
  arrNear(coordBlock(t, /^geometry=\{/, /^\}/, [1, 2, 3]), WATER.xyz, 5e-11, 'coordinates');
  const sp = W.molpro(base({ task: 'sp' }));
  hasNot(sp, 'optg', 'sp'); hasNot(sp, 'frequencies', 'sp');
  const fq = W.molpro(base({ task: 'freq' }));
  hasNot(fq, 'optg', 'freq without optg'); has(fq, '{frequencies}', 'freq');
  const open = W.molpro(base({ charge: 1, multiplicity: 2 }));
  has(open, 'set,charge=1\nset,spin=1\n{uks,b3lyp}', 'open-shell UKS');
  has(W.molpro(base({ method: 'CCSD(T)', task: 'sp' })), '{hf}\n{ccsd(t)}', 'closed-shell CCSD(T)');
  has(W.molpro(base({ method: 'CCSD(T)', multiplicity: 3 })), '{rhf}\n{uccsd(t)}', 'open-shell RHF-UCCSD(T)');
  has(W.molpro(base({ method: 'MP2' })), '{hf}\n{mp2}', 'MP2');
  has(W.molpro(base({ title: 'a, b; c!' })), '***,a b c\n', 'title separators removed');
  has(prog('molpro').write(base({ solvent: 'water' })).note, 'Implicit solvent is not written', 'solvent note');
});

// ------------------------------------------------------------------ Turbomole
test('Turbomole coord', () => {
  const t = W.turbomole(Object.assign({}, PT));
  const L = t.trim().split('\n');
  eq(L[0], '$coord', 'starts with $coord'); eq(L[L.length - 1], '$end', 'ends with $end');
  const rows = L.slice(1, -1).map(s => s.trim().split(/\s+/));
  eq(JSON.stringify(rows.map(q => q[3])), JSON.stringify(PT.Z.map(z => BM.el.symbol(z).toLowerCase())), 'lower-case element symbols last');
  arrNear([].concat(...rows.map(q => q.slice(0, 3).map(num))), PT.xyz.map(v => v / BOHR), 5e-14, 'coordinates in bohr, 14 decimals');
  ok(rows.every(q => q.slice(0, 3).every(s => /\.\d{14}$/.test(s))), '14 decimals');
  const r = prog('turbomole').write(base({ task: 'optfreq', solvent: 'water' }));
  eq(r.files[0].name, 'coord', 'file name coord'); eq(r.text, r.files[0].text, 'text = coord');
  has(r.note, 'define', 'note: use define'); has(r.note, 'cosmoprep', 'note: solvent via cosmoprep'); has(r.note, 'aoforce', 'note: commands');
  has(W.turbomole(base({ extra: '$user-defined bonds\n1-2, 1-3' })), '$user-defined bonds\n1-2, 1-3\n$end\n', 'extra data groups before $end');
});

// ------------------------------------------------------------------ cross-program checks
test('charge, multiplicity and precision everywhere', () => {
  const opts = Object.assign({}, PT, { charge: 2, multiplicity: 3, task: 'opt', nprocs: 1, title: 'Pt test' });
  const texts = {
    psi4: [W.psi4(opts), /molecule \{\n2 3\n/], qchem: [W.qchem(opts), /\$molecule\n2 3\n/],
    nwchem: [W.nwchem(opts), /charge 2[\s\S]*mult 3/], gamess: [W.gamess(opts), /ICHARG=2 MULT=3/],
    mopac: [W.mopac(Object.assign({}, opts, { method: 'PM7' })), /^PM7 CHARGE=2 TRIPLET\n/],
    molpro: [W.molpro(opts), /set,charge=2\nset,spin=2\n/], xtb: [W.xtb(Object.assign({}, opts, { method: 'GFN2-xTB' })).command, /--chrg 2 --uhf 2/],
    gaussian: [prog('gaussian').write(opts).text, /\n2 3\nPt /], orca: [prog('orca').write(opts).text, /\* xyz 2 3\n/]
  };
  for (const [id, [text, re]] of Object.entries(texts)) has(text, re, id + ' charge/multiplicity');
  // Every coordinate value appears in each Ångström input with at least 10 decimals (8 for the Gaussian/ORCA writers).
  for (const id of ['psi4', 'qchem', 'nwchem', 'gamess', 'mopac', 'molpro']) {
    const text = texts[id][0];
    ok(PT.xyz.every(v => text.indexOf(v.toFixed(10)) >= 0), id + ': all coordinates at full precision');
  }
  for (const id of ['gaussian', 'orca']) ok(PT.xyz.every(v => texts[id][0].indexOf(v.toFixed(8)) >= 0), id + ': coordinates to 1e-8 Å');
});

// ---------------------------------------------------------------------- summary
for (const n of notes) console.log('  note: ' + n);
console.log('\nBlueMolecule writer tests: ' + pass + ' passed, ' + fail + ' failed');
if (fail) {
  for (const f of failures) console.log('  FAIL ' + f);
  process.exit(1);
}
