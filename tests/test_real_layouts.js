#!/usr/bin/env node
/* Reader fixes from the 2.0 check against real output files (cclib's and MDAnalysis's test collections): each case
 * follows the layout of the real files, which the imitation fixtures had missed.
 *   - a job still running: a geometry table cut while it is written makes no frame with part of the atoms
 *     (GAMESS, NWChem, Psi4, Q-Chem, xtb, MOPAC, CP2K), nor does a cut "User input" echo (Q-Chem)
 *   - Q-Chem inputs of several jobs (@@@): one molecule, every job type
 *   - ORCA inputs with %paras values in braces; ORCA / xtb .engrad files; MOPAC 2012 archives ("CHARGE" column) and
 *     .mop inputs without a method keyword; a Turbomole mpshift.out is no PDB file
 *   - Gaussian MP4 ("E4(SDTQ)= ... UMP4(SDTQ)= ..."), MP5, CCD and Gaussian 09 CCSD energies
 *   - ORCA TD-DFT: the "FINAL SINGLE POINT ENERGY" of the followed excited state; ORCA 5 "spin forbidden" triplet rows
 *     and ORCA 6 per-multiplicity rows, with their CD
 *   - fchk files of TD jobs and of Q-Chem; GAMESS-US doubled population tables; Firefly "FINAL ENERGY IS"; NWChem's
 *     5-decimal gross populations; an unreadable element; validation messages that name the reader
 *   Run:  node tests/test_real_layouts.js */
'use strict';
const fs = require('fs');
const path = require('path');
const ROOT = path.resolve(__dirname, '..');
for (const m of ['elements', 'bigfile', 'parsers', 'condense', 'readers-struct', 'writers', 'readers-qc', 'spectrum', 'geom', 'chem', 'analysis']) require(path.join(ROOT, 'src', m + '.js'));
const BM = globalThis.BM;
const FIX = path.join(ROOT, 'tests', 'fixtures');
const last = a => a[a.length - 1];
const read = f => fs.readFileSync(path.join(FIX, f), 'utf8');

let passes = 0, failures = 0;
const fails = [];
const check = (c, msg) => { if (c) passes++; else { failures++; fails.push(msg); } };
const near = (a, b, tol) => typeof a === 'number' && Math.abs(a - b) <= tol;
const dist = (f, i, j) => Math.hypot(f.xyz[3 * i] - f.xyz[3 * j], f.xyz[3 * i + 1] - f.xyz[3 * j + 1], f.xyz[3 * i + 2] - f.xyz[3 * j + 2]);
const throws = (fn, re) => { try { fn(); return false; } catch (e) { return re.test(e.message); } };
/** A group of checks; an exception fails the group instead of ending the run. */
function section(name, fn) { try { fn(); } catch (e) { check(false, name + ' threw: ' + (e && e.message)); } }

// ------------------------------------------------------------------ a job still running: whole geometries only
for (const f of ['qc/gamess_water_opt_hess.log', 'qc/nwchem_water_optfreq.out', 'qc/nwchem_water_tddft.out', 'qc/psi4_water_optfreq.out',
  'qc/qchem_water_optfreq.out', 'qc/qchem_water_tddft.out', 'qc/xtb_water_ohess.out', 'qc/mopac_water_opt_force.out', 'qc/cp2k_water_geoopt.out']) section(f, () => {
  const text = read(f), name = path.basename(f);
  const n = BM.parseFile(name, text).frames[0].Z.length;
  let cuts = 0, bad = [];
  for (let k = 1; k < 80; k++) {
    const c = Math.floor(text.length * k / 80);
    let ds;
    try { ds = BM.parseFile(name, text.slice(0, c)); } catch (e) { continue; }
    cuts++;
    ds.frames.forEach((fr, i) => { if (fr.Z.length !== n) bad.push('@' + c + ' frame ' + (i + 1) + ': ' + fr.Z.length + ' of ' + n + ' atoms'); });
  }
  check(cuts > 20 && !bad.length, name + ' cut at 79 places (' + cuts + ' read): every frame has all ' + n + ' atoms ' + bad.slice(0, 3).join('; '));
});

// ------------------------------------------------------------------ Q-Chem inputs of several jobs
section('Q-Chem inputs of several jobs', () => {
  const mol = (z) => '$molecule\n0 1\nO   0.000   0.000   ' + z + '\nH   0.000   0.757  -0.469\nH   0.000  -0.757  -0.469\n$end\n';
  const rem = jt => '\n$rem\nJOBTYPE ' + jt + '\nMETHOD b3lyp\nBASIS 6-31G*\n$end\n';
  const two = mol('0.117') + rem('opt') + '\n@@@\n\n' + mol('0.120') + rem('freq');
  const ds = BM.parseFile('two_jobs.in', two);
  check(ds.format === 'qchem-input' && ds.frames.length === 1 && ds.frames[0].Z.length === 3, 'Q-Chem two jobs: one water, not ' + ds.frames.map(f => f.Z.length) + ' atoms');
  check(near(ds.frames[0].xyz[2], 0.117, 1e-9), 'Q-Chem two jobs: the first job\'s geometry');
  check(ds.meta.jobTypes.indexOf('opt') >= 0 && ds.meta.jobTypes.indexOf('freq') >= 0, 'Q-Chem two jobs: both job types: ' + ds.meta.jobTypes);
  const read2 = mol('0.117') + rem('opt') + '\n@@@\n\n$molecule\nread\n$end\n' + rem('freq');
  check(BM.parseFile('r.in', read2).frames[0].Z.length === 3, 'Q-Chem: a second job that reads the molecule');
  // the output's echo of the input, cut inside $molecule: no frame of 1 or 2 atoms
  const out = read('qc/qchem_water_optfreq.out'), at = out.indexOf('$molecule');
  let bad = 0;
  for (let c = at + 10; c < at + 140; c += 7) {
    try { const d = BM.parseFile('q.out', out.slice(0, c)); if (d.frames.some(f => f.Z.length !== 3)) bad++; } catch (e) { /* no geometry yet */ }
  }
  check(!bad, 'Q-Chem output cut inside the echoed $molecule: no partial molecule (' + bad + ')');
});

// ------------------------------------------------------------------ ORCA input with %paras
section('ORCA input with %paras', () => {
  const inp = '! B3LYP def2-SVP\n%paras\n  R = 0.96\n  A = 104.5\n  RS = 0.90, 1.10, 5\nend\n\n* int 0 1\nO 0 0 0 0.0 0.0 0.0\nH 1 0 0 {R} 0.0 0.0\nH 1 2 0 {RS} {A} 0.0\n*\n';
  const ds = BM.parseFile('scan.inp', inp), f = ds.frames[0];
  check(ds.format === 'orca-input' && near(dist(f, 0, 1), 0.96, 1e-9) && near(dist(f, 0, 2), 0.90, 1e-9), 'ORCA %paras: the bond lengths ' + dist(f, 0, 1) + ', ' + dist(f, 0, 2));
  const ang = Math.acos(((f.xyz[3] - f.xyz[0]) * (f.xyz[6] - f.xyz[0]) + (f.xyz[4] - f.xyz[1]) * (f.xyz[7] - f.xyz[1]) + (f.xyz[5] - f.xyz[2]) * (f.xyz[8] - f.xyz[2])) / (dist(f, 0, 1) * dist(f, 0, 2))) * 180 / Math.PI;
  check(near(ang, 104.5, 1e-6), 'ORCA %paras: the angle ' + ang);
  check(ds.warnings.some(w => /scans RS/.test(w)) && ds.meta.jobTypes.indexOf('scan') >= 0, 'ORCA %paras: the scanned parameter is named: ' + ds.warnings);
  check(throws(() => BM.parseFile('x.inp', inp.replace('{A}', '{ANGLE}')), /\{ANGLE\}.*%paras/), 'ORCA %paras: an undefined name is named in the message');
});

// ------------------------------------------------------------------ .engrad, MOPAC, a Turbomole file that is no PDB
section('.engrad, MOPAC, a Turbomole file that is no PDB', () => {
  const B = 0.529177210903;
  const eng = '#\n# Number of atoms\n#\n 3\n#\n# The current total energy in Eh\n#\n   -76.400000000000\n#\n# The current gradient in Eh/bohr\n#\n' +
    ['0.0', '0.0', '0.012', '0.0', '0.005', '-0.006', '0.0', '-0.005', '-0.006'].map(v => '      ' + v).join('\n') +
    '\n#\n# The atomic numbers and current coordinates in Bohr\n#\n   8     0.0000000    0.0000000    0.2217000\n   1     0.0000000    1.4309000   -0.8867000\n   1     0.0000000   -1.4309000   -0.8867000\n';
  const ds = BM.parseFile('water.engrad', eng), f = ds.frames[0];
  check(ds.format === 'engrad' && f.Z.join() === '8,1,1' && near(f.xyz[4], 1.4309 * B, 1e-9) && f.energy === -76.4, '.engrad: atoms, Bohr coordinates and energy');
  check(f.opt && near(f.opt.maxForce, 0.012, 1e-12) && near(f.opt.rmsForce, Math.sqrt((0.012 ** 2 * 1 + 0.005 ** 2 * 2 + 0.006 ** 2 * 2) / 9), 1e-12), '.engrad: gradient max and RMS');
  check(BM.detectFormat('x.engrad', eng) === 'engrad', '.engrad is not taken for a Gaussian input (it starts with "#")');
  const arc = read('qc/mopac_water.arc').replace(/^(\s*FINAL GEOMETRY OBTAINED)\s*$/m, '$1                                    CHARGE');
  check(/CHARGE/.test(arc.split('FINAL GEOMETRY OBTAINED')[1].split('\n')[0]), '(the archive heading has its CHARGE column)');
  check(BM.parseFile('w.arc', arc).frames.length >= 1, 'MOPAC 2012 archive with a CHARGE column after FINAL GEOMETRY OBTAINED');
  const mop = ' UHF TRIPLET\nO2\n\n   O       -5.40592             1        -0.47619             1        -0.00000             1 \n   O       -4.36197             1        -0.31061             1        -0.00000             1 \n\n';
  const dm = BM.parseFile('o2_trip.mop', mop);
  check(dm.format === 'mopac-input' && dm.frames[0].Z.join() === '8,8' && near(dist(dm.frames[0], 0, 1), Math.hypot(1.04395, 0.16558), 1e-6), '.mop without a method keyword is a MOPAC input');
  // MOPAC UHF: ALPHA and BETA EIGENVALUES, occupied by the alpha and beta electrons (no orbitals were read before)
  const uhf = read('qc/mopac_water_opt_force.out').replace(/^(\s*)NO\. OF FILLED LEVELS\s*=\s*4$/m, '$1NO. OF ALPHA ELECTRONS  =          5\n$1NO. OF BETA  ELECTRONS  =          4')
    .replace(/^(\s*)EIGENVALUES\n\n(.*)\n/m, '$1ALPHA EIGENVALUES\n$2\n\n\n$1 BETA EIGENVALUES\n  -40.1 -19.5 -14.2 -11.9   3.1   4.4\n');
  check(/BETA EIGENVALUES/.test(uhf) && /ALPHA ELECTRONS/.test(uhf), '(MOPAC UHF layout made)');
  const ou = BM.parseFile('uhf.out', uhf).orbitals;
  check(ou && ou.alpha.nOcc === 5 && ou.beta && ou.beta.nOcc === 4 && ou.beta.energies.length === 6 && near(ou.beta.energies[3] * 27.211386245988, -11.9, 1e-6),
    'MOPAC UHF: alpha and beta orbital energies with their occupations: ' + JSON.stringify(ou && { a: ou.alpha.nOcc, b: ou.beta && ou.beta.nOcc }));
  const dv = BM.parseFile('w.out', read('qc/mopac_water_opt_force.out')).vibrations;
  check(dv && dv.ir && dv.irUnit === 'D²', 'MOPAC IR intensities carry their unit (D², not km/mol): ' + (dv && dv.irUnit));
  if (globalThis.BM.spectrum) check(BM.spectrum.ir(dv, { yMode: 'epsilon', scale: 1, shape: 'lorentzian', fwhm: 12, xmin: 400, xmax: 4000 }).yLabel === 'Relative absorbance', 'MOPAC IR is never shown as ε');
  check(BM.qcReaders.cp2kInput(read('struct/h2o.inp'), 'h2o.inp').frames[0].Z.join() === '8,1,1', 'CP2K input: a kind named "Wat" with ELEMENT O is oxygen');
  // the two CP2K input readers (structure files: .inp; calculation files: .restart) must read the same
  const cp = {
    fixture: read('struct/h2o.inp'),
    calc: '&GLOBAL\n RUN_TYPE GEO_OPT\n&END GLOBAL\n&FORCE_EVAL\n &DFT\n  CHARGE -1\n  MULTIPLICITY 2\n  &XC\n   &XC_FUNCTIONAL PBE\n   &END XC_FUNCTIONAL\n  &END XC\n &END DFT\n' +
      ' &SUBSYS\n  &CELL\n   ABC [bohr] 18.8973 18.8973 18.8973\n  &END CELL\n  &COORD\n   SCALED\n   O 0.5 0.5 0.5\n   H 0.5 0.575 0.45\n  &END COORD\n &END SUBSYS\n&END FORCE_EVAL\n',
    nm: '&FORCE_EVAL\n &SUBSYS\n  &CELL\n   A [nm] 1 0 0\n   B [nm] 0 1 0\n   C [nm] 0 0 1\n  &END CELL\n  &COORD\n   UNIT nm\n   O1 0 0 0\n   H 0 0.076 -0.047\n  &END COORD\n  &KIND O1\n   ELEMENT O\n  &END KIND\n &END SUBSYS\n&END FORCE_EVAL\n',
  };
  for (const [k, t] of Object.entries(cp)) {
    const a = BM.structReaders.parsers.cp2kIn(t, k + '.inp'), b = BM.qcReaders.cp2kInput(t, k + '.restart');
    const sig = d => [d.frames[0].Z.join(), Array.from(d.frames[0].xyz).map(v => v.toFixed(4)).join(), d.cell ? d.cell.flat().map(v => v.toFixed(4)).join() : '-',
      d.meta.charge, d.meta.multiplicity, d.meta.method, d.meta.jobTypes.join('/')].join(' | ');
    check(sig(a) === sig(b), 'the two CP2K input readers agree (' + k + '): ' + sig(a) + ' vs ' + sig(b));
  }
  check(BM.parseFile('c.inp', cp.calc).meta.charge === -1 && BM.parseFile('c.inp', cp.calc).meta.method === 'PBE', 'a CP2K input shows its charge and functional');
  check(BM.detectFormat('mpshift.out', ' ATOM  c    1      ISOTROPIC:      114.0801043       ANISOTROPIC:      150.2683398\n') !== 'pdb', 'a Turbomole mpshift.out line is no PDB atom');
  check(BM.detectFormat('x.out', 'ATOM      1  N   ALA A   1      11.104   6.134  -6.504  1.00  0.00           N\n') === 'pdb', 'a real PDB line is still found in a .out file');
});

// ------------------------------------------------------------------ Gaussian correlated energies (real layouts)
section('Gaussian correlated energies (real layouts)', () => {
  const g9 = BM.parseFile('g09_water_ccsd_t.log', read('review/g09_water_ccsd_t.log')).frames[0];
  check(near(g9.energies.CCSD, -75.017683639, 1e-9) && near(g9.energies['MP4(SDQ)'], -75.016013656, 1e-9) && near(g9.energy, -75.017760423, 1e-9),
    'Gaussian 09 CCSD(T): CCSD from "CI/CC converged" and MP4(SDQ) from mid-line UMP4: ' + JSON.stringify(g9.energies));
  const m5 = BM.parseFile('g16_water_mp5.log', read('review/g16_water_mp5.log')).frames[0];
  check(near(m5.energies.MP4, -75.016068053, 1e-9) && near(m5.energies.MP5, -75.017172934, 1e-9) && near(m5.energy, -75.017172934, 1e-9),
    'Gaussian 16 MP5 job: MP4(SDTQ) and MP5, and MP5 is the energy (was MP3): ' + JSON.stringify(m5.energies));
  const ccd = read('review/g09_water_ccsd_t.log').replace('#p ccsd(t)/sto-3g', '#p ccd/sto-3g  ').replace(/ CCSD\(T\)= .*\n/, '');
  const cd = BM.parseFile('ccd.log', ccd).frames[0];
  check(near(cd.energies.CCD, -75.017683639, 1e-9) && cd.energy === cd.energies.CCD, 'a CCD job keeps its energy as CCD: ' + JSON.stringify(cd.energies));
});

// ------------------------------------------------------------------ ORCA TD-DFT with triplets, ORCA 5 and 6 layouts
for (const [f, v] of [['review/orca5_water_tddft_triplets.out', 5], ['review/orca6_water_tddft_triplets.out', 6]]) section('ORCA ' + v + ' TD-DFT', () => {
  const ds = BM.parseFile(path.basename(f), read(f)), fr = ds.frames[ds.frames.length - 1], es = ds.excitedStates;
  check(near(fr.energy, -76.31794234, 1e-8) && near(fr.energies['Excited state'], -76.0357033435, 1e-9) && fr.energies.Final === undefined,
    'ORCA ' + v + ' TD-DFT single point: the energy is the ground state, the followed state is kept apart: ' + JSON.stringify(fr.energies));
  const sing = es.filter(s => /^Singlet/.test(s.label)), trip = es.filter(s => /^Triplet/.test(s.label));
  check(es.length === 5 && sing.length === 3 && trip.length === 2, 'ORCA ' + v + ': 3 singlets and 2 triplets, each once: ' + es.map(s => s.index + s.label).join(' '));
  check(trip.every(s => s.f === 0 && s.transitions.length === 1) && near(trip[0].energyEV, 6.987, 2e-3), 'ORCA ' + v + ': triplets with f = 0 and their transitions');
  const s1 = sing.find(s => s.index === 1);
  check(s1 && near(s1.energyEV, 7.680, 2e-3) && near(s1.R, -12.34567, 1e-9) && near(s1.f, 0.016487919, 1e-12), 'ORCA ' + v + ': singlet 1 keeps its CD and f');
  // ADC(2): the excitations are added to the MP2 ground state ("Second-order corrected state energy"), not to the SCF
  const adc = read(f).replace(/^(Total Energy\s*:\s*)-76\.31794234/m, '$1-76.01234567')
    .replace(/^(FINAL SINGLE POINT ENERGY)/m, 'Second-order corrected state energy     ...   -76.317942340\n\n$1');
  check(/Second-order/.test(adc) && /-76\.01234567/.test(adc), '(ADC(2) layout made)');
  const fa = BM.parseFile('adc2.out', adc).frames.slice(-1)[0];
  check(near(fa.energy, -76.31794234, 1e-8) && fa.energies.MP2 === fa.energy && near(fa.energies['Excited state'], -76.0357033435, 1e-9),
    'ORCA ' + v + ' ADC(2): the energy is the MP2 ground state, the followed state is kept apart: ' + JSON.stringify(fa.energies));
  check(trip.every(s => !(Math.abs(s.R || 0) > 0)), 'ORCA ' + v + ': no singlet CD on a triplet');
  if (v === 6) check(near(s1.Rvel, -11, 1e-9), 'ORCA 6: velocity CD of singlet 1');
});

// ------------------------------------------------------------------ fchk: TD jobs, Q-Chem
section('fchk: TD jobs, Q-Chem', () => {
  const base = read('water_freq.fchk');
  const td = base.replace(/^(Total Energy\s+R\s+)(\S+)/m, '$1-7.459032049491202E+01').replace(/^(SCF Energy\s+R\s+\S+\n)/m, '$1CIS Energy                                 R     -7.459032049491202E+01\n');
  check(/CIS Energy/.test(td), '(the TD fchk has a CIS Energy)');
  const d = BM.parseFile('td.fchk', td), f = d.frames[0];
  check(near(f.energies.TD, -74.59032049491202, 1e-9) && f.energy === f.energies.SCF && d.meta.jobTypes.indexOf('td') >= 0,
    'fchk of a TD job: the energy is the ground state, the excited state is energies.TD: ' + JSON.stringify(f.energies));
  const qc = base.replace(/^[^\n]*/, 'Jobname.Temp').replace(/^Dipole Moment[\s\S]*?\n(?=\S)/m, '') + 'Dipole_Data                                R   N=           3\n  1.00000000E-15 -2.00000000E-15 -1.67349843E+00\n';
  const q = BM.parseFile('qchem.fchk', qc);
  check(/Q-Chem/.test(q.program) && q.dipole && near(q.dipole.total, 1.67349843, 1e-8), 'Q-Chem fchk: program and dipole in Debye: ' + q.program + ' ' + (q.dipole && q.dipole.total));
});

// ------------------------------------------------------------------ GAMESS-US / Firefly, NWChem
section('GAMESS-US / Firefly, NWChem', () => {
  const g = read('qc/gamess_water_opt_hess.log');
  const m = /( *TOTAL MULLIKEN AND LOWDIN ATOMIC POPULATIONS\n(?: *-+\n)? *ATOM +MULL\.POP\. +CHARGE +LOW\.POP\. +CHARGE\n)((?: +\d+ +\S+ +\S+ +\S+ +\S+ +\S+\n)+)/.exec(g);
  check(!!m, '(the GAMESS fixture has a population table)');
  if (m) {
    // some GAMESS-US builds: the heading twice (first without CHARGE columns) and a populations-only row before each full row
    const rows = m[2].split('\n').filter(Boolean).map(r => r.replace(/^( +\d+ +\S+ +)(\S+) +(\S+) +(\S+) +(\S+)$/, '$1$2                     $4\n$1$2   $3         $4    $5')).join('\n') + '\n';
    const dbl = g.replace(m[0], '          TOTAL MULLIKEN AND LOWDIN ATOMIC POPULATIONS\n       ATOM         MULL.POP.                    LOW.POP.\n\n' + m[1] + rows);
    const want = BM.parseFile('g.log', g).charges, got = BM.parseFile('g.log', dbl).charges;
    check(got && got.Mulliken && got.Mulliken.join() === want.Mulliken.join() && got.Loewdin.join() === want.Loewdin.join(),
      'GAMESS-US doubled population table: the same charges ' + (got && got.Mulliken));
  }
  // (without the optimisation's "NSERCH ... E=" lines, which give energies too: a Firefly single point has none)
  const g0 = g.replace(/^.*NSERCH:?\s*\d+\s+E=.*\n/gm, ''), ff = g0.replace(/FINAL R-B3LYP ENERGY IS/g, 'FINAL ENERGY IS');
  const e1 = BM.parseFile('g.log', g0).frames.map(f => f.energy).join(), e2 = BM.parseFile('ff.log', ff).frames.map(f => f.energy).join();
  check(e1 === e2 && /-76\.4/.test(e2), 'Firefly "FINAL ENERGY IS" lines give the same energies: ' + e2);
  const nw = read('qc/nwchem_water_optfreq.out');
  const at = nw.lastIndexOf('Mulliken analysis of the total density');
  const end = nw.lastIndexOf('\n', nw.indexOf('Multipole analysis of the density', at));
  const gross = '\n\n          ----- Total      gross population on atoms ----\n\n              1  O    8.0     8.66212\n              2  H    1.0     0.66894\n              3  H    1.0     0.66894\n';
  const dsN = BM.parseFile('nw.out', nw.slice(0, end) + gross + nw.slice(end));
  check(dsN.charges && near(dsN.charges.Mulliken[0], -0.66212, 1e-9) && near(dsN.charges.Mulliken[1], 0.33106, 1e-9), 'NWChem: Mulliken charges to 5 decimals from the gross populations: ' + (dsN.charges && dsN.charges.Mulliken));
  const nwf = read('qc/nwchem_water_optfreq.out').replace(/^ *P?\.?Frequency .*\n/m, '');
  check((() => { try { BM.parseFile('nw.out', nwf); return true; } catch (e) { return !/Cannot (set|read) prop/.test(e.message); } })(), 'NWChem normal modes without their Frequency line do not break the reader');
  // Q-Chem "Orbital Energies (a.u.) and Symmetries": a line of labels under each line of energies (only the first
  // line of energies was read, fewer than the occupied orbitals, and the level diagram broke)
  const qs = read('qc/qchem_water_optfreq.out').replace('Orbital Energies (a.u.)\n', 'Orbital Energies (a.u.) and Symmetries\n')
    .replace(' Alpha MOs\n', ' Alpha MOs, Restricted\n')
    .replace('-19.1351  -0.9971  -0.5203  -0.3662  -0.2921\n', '-19.1351  -0.9971  -0.5203\n  1 A1    2 A1    1 B2\n -0.3662  -0.2921\n  3 A1    1 B1\n')
    .replace('  0.0654   0.1491\n', '  0.0654   0.1491\n  4 A1    2 B2\n');
  check(/1 B2/.test(qs), '(Q-Chem layout with symmetries made)');
  const oq = BM.parseFile('q.out', qs).orbitals;
  check(oq && oq.alpha.energies.length === 7 && oq.alpha.nOcc === 5 && oq.alpha.energies[4] === -0.2921, 'Q-Chem orbital energies with symmetry labels: all 7 read, HOMO -0.2921: ' + (oq && oq.alpha.energies));
});

// ------------------------------------------------------------------ NWChem ab initio MD (task dft qmd)
section('NWChem QMD', () => {
  // each step prints its energy, then its geometry in a.u. in the gradient table, then the QMD summary (only the first
  // geometry was read before)
  const B = 0.529177210903;
  const nw = read('qc/nwchem_water_optfreq.out');
  const head = nw.slice(0, nw.indexOf('      Atomic Mass')) + '\n                                 NWChem QMD Module\n                                 -----------------\n';
  const grad = (dz, E) => `\n         Total DFT energy =     ${E}\n\n                         DFT ENERGY GRADIENTS\n\n    atom               coordinates                        gradient\n` +
    '                 x          y          z           x          y          z\n' +
    [['O', 0, 0, 0.117300 + dz], ['H', 0, 0.757200, -0.469200], ['H', 0, -0.757200, -0.469200]].map((a, k) =>
      `   ${k + 1} ${a[0]}    ${a.slice(1).map(v => (v / B).toFixed(6).padStart(11)).join('')}    0.000000   0.000000  -0.004108`).join('\n') + '\n';
  const info = (k, t, E) => `\n            QMD Run Information\n            -------------------\n  Time elapsed (fs) :                     ${t}\n  Pot. energy (a.u.):        ${k}          ${E}\n`;
  const qmd = head + grad(0, '-76.400234495670') + grad(0.01, '-76.399603733417') + info(1, '0.604721', '-76.399604') + grad(0.02, '-76.396794145899') + info(2, '1.209442', '-76.396794');
  const dq = BM.parseFile('qmd.out', qmd);
  check(dq.frames.length === 3 && dq.frames.map(f => f.energy).join() === '-76.40023449567,-76.399603733417,-76.396794145899',
    'NWChem QMD: every MD step is a geometry with its energy: ' + dq.frames.map(f => f.energy).join());
  check(near(dq.frames[2].xyz[2], 0.1373, 2e-6) && dq.frames.map(f => f.label).join() === 't = 0.00 fs,t = 0.60 fs,t = 1.21 fs',
    'NWChem QMD: geometries from a.u., labelled with their time: ' + dq.frames.map(f => f.label).join());
  // a step cut in its gradient table (a job still running) adds no partial geometry
  const cut = qmd.slice(0, qmd.lastIndexOf('   2 H'));
  check(BM.parseFile('qmd.out', cut).frames.length === 2, 'NWChem QMD cut in a gradient table: no partial geometry');
});

// ------------------------------------------------------------------ ORCA 6 QM/MM: the layout of real outputs
section('ORCA 6 QM/MM', () => {
  // Each step prints FINAL SINGLE POINT ENERGY (MM) before its geometry, then the QM region's energy and (QM/MM).
  // All three were read as one "final" energy, so each step showed the next step's MM energy (about +2 Eh in a real
  // output) and only the last step the QM/MM total, which made ΔE and the energy plots meaningless.
  const t = read('review/orca6_qmmm_qmregion_opt.out');
  const d = BM.parseFile('qmmm.out', t);
  const E = d.frames.map(f => f.energies || {});
  check(d.frames.length === 3 && d.frames.every((f, k) => f.energy === E[k]['QM/MM']) && E.every(e => e.Final === undefined),
    'ORCA 6 QM/MM: every step\'s energy is its QM/MM total: ' + d.frames.map(f => f.energy).join());
  check(E.map(e => e.MM).join() === '0.0123456789,0.0119876543,0.0117654321' && near(E[0].QM, -437.999999321, 1e-9) && near(E[0]['QM/MM'], E[0].QM + E[0].MM, 1e-9),
    'ORCA 6 QM/MM: each step keeps its own MM energy (printed before its geometry) and its QM energy: ' + JSON.stringify(E[0]));
  const c = BM.analysis.convergenceSeries(d), dE = c.series.find(s => s.key === 'dE');
  check(dE && dE.values.every(v => Math.abs(v) < 2) && near(dE.values[0], (E[0]['QM/MM'] - E[2]['QM/MM']) * 627.5094740631, 1e-6), 'ORCA 6 QM/MM: ΔE vs final in kcal/mol: ' + (dE && dE.values.map(v => v.toFixed(3))));
  // only the QM region is printed (the full system is 120 atoms): its QM atoms, then the link atom, bonded to S
  const L = d.layers;
  check(L && L.subsystemOnly && L.linkAtoms.length === 1 && L.linkAtoms[0].atom === 5 && L.linkAtoms[0].bondedTo === 0 && L.labels.L === 'MM (system.ORCAFF.prms)',
    'ORCA 6 QM/MM: QM region alone, with its link atom: ' + JSON.stringify(L && { only: L.subsystemOnly, links: L.linkAtoms }));
  // a job still running, cut after the next step's MM energy and before its geometry: that energy waits for the step
  const cut = t.slice(0, t.indexOf('FINAL SINGLE POINT ENERGY (MM)        0.011987') + 60);
  const dc = BM.parseFile('qmmm.out', cut);
  check(dc.frames.length === 1 && dc.frames[0].energies.MM === 0.0123456789, 'ORCA 6 QM/MM cut before a step\'s geometry: the previous step keeps its own MM energy: ' + JSON.stringify(dc.frames[0].energies));
  // cut after a step's QM energy, before its (QM/MM) line: the step has no total yet, so no energy
  const cut2 = t.slice(0, t.indexOf('FINAL SINGLE POINT ENERGY (QM/MM)    -437.98859') - 2);
  const d2 = BM.parseFile('qmmm.out', cut2);
  check(d2.frames.length === 2 && d2.frames[1].energy === null && d2.frames[0].energy === E[0]['QM/MM'], 'ORCA 6 QM/MM cut before a step\'s total: that step has no energy yet: ' + d2.frames.map(f => f.energy).join());
});

// ------------------------------------------------------------------ frequencies by finite differences (#32)
section('finite-difference frequencies', () => {
  // Psi4 FINDIF: one computation per displaced geometry, the reference first (Psi4 1.4+); they were all steps, and the
  // vibrations went to the last displaced geometry
  const ps = read('qc/psi4_water_optfreq.out');
  const hv = ps.indexOf('  ==> Harmonic Vibrational Analysis <==');
  const B = 0.529177210903, ref = [[8, 0, 0, 0.119018], [1, 0, 0.759242, -0.475927], [1, 0, -0.759242, -0.475927]];
  const comp = (dz, dy, E, mu) => '    Geometry (in Bohr), charge = 0, multiplicity = 1:\n\n       Center              X                  Y                   Z       \n    ------------   -----------------  -----------------  -----------------\n' +
    ref.map((a, k) => '         ' + ['O', 'H', 'H'][k] + '   ' + [a[1], a[2] + (k === 1 ? dy : 0), a[3] + (k === 0 ? dz : 0)].map(v => (v / B).toFixed(12).padStart(19)).join('')).join('\n') +
    '\n\n  @DF-RKS Final Energy:   ' + E + '\n\n Dipole X            :          0.0000000            0.0000000            0.0000000\n' +
    ' Dipole Y            :          0.0000000            0.0000000            0.0000000\n Dipole Z            :         -0.1594632            0.9832395            ' + mu + '\n\n';
  const findif = ps.slice(0, ps.indexOf('    Geometry (in Angstrom)')) + '    Number of geometries (including reference) is 5.\n\n  //    FiniteDifference Computations  //\n\n' +
    comp(0, 0, '-76.40987172299999', '0.8237763') + comp(0.004, 0, '-76.40986501234567', '0.8251000') + comp(-0.004, 0, '-76.40986412345678', '0.8224000') +
    comp(0, 0.004, '-76.40986623456789', '0.8239000') + comp(0, -0.004, '-76.40986634567891', '0.8236000') +
    '  //      FiniteDifference Results     //\n\n  5 gradients passed in, including the reference geometry.\n\n' + ps.slice(hv);
  const df = BM.parseFile('fd.out', findif);
  check(df.frames.length === 1 && df.frames[0].energy === -76.40987172299999 && df.vibrations && df.vibrations.frameIndex === 0,
    'Psi4 FINDIF: the displaced geometries are no steps; the reference keeps its energy and the vibrations: ' + df.frames.length + ' frames, energy ' + df.frames.map(f => f.energy));
  check(df.dipole && near(df.dipole.z, 0.8237763 * 2.541746473, 1e-5), 'Psi4 FINDIF: the dipole is the reference geometry\'s, not the last displacement\'s: ' + JSON.stringify(df.dipole));
  // an optimisation, then finite differences from the optimised geometry: the steps stay, the displacements go
  const opt = ps.slice(0, hv) + '    Number of geometries (including reference) is 5.\n\n' + comp(0, 0, '-76.40987172299999', '0.8237763') +
    comp(0.004, 0, '-76.40986501234567', '0.8251000') + comp(-0.004, 0, '-76.40986412345678', '0.8224000') + '  3 gradients passed in, including the reference geometry.\n\n' + ps.slice(hv);
  const n0 = BM.parseFile('w.out', ps).frames.length, dop = BM.parseFile('opt.out', opt);
  check(dop.frames.length === n0 && dop.vibrations.frameIndex === n0 - 1, 'Psi4 optimisation + FINDIF: the optimisation steps stay (' + n0 + '), no displaced geometry is added: ' + dop.frames.length);
  // a job still running, in its displacements: the reference only
  check(BM.parseFile('fd.out', findif.slice(0, findif.indexOf('  //      FiniteDifference Results'))).frames.length === 1, 'Psi4 FINDIF still running: only the reference geometry');

  // Q-Chem anharmonic (VPT2): displaced geometries after "Starting finite difference calculation", then the reference
  // printed again before the anharmonic analysis
  const qc = read('qc/qchem_water_optfreq.out');
  const geo = (dz, E) => ' ----------------------------------------------------------------\n             Standard Nuclear Orientation (Angstroms)\n    I     Atom           X                Y                Z\n ----------------------------------------------------------------\n' +
    ref.map((a, k) => '    ' + (k + 1) + '      ' + ['O', 'H', 'H'][k] + '   ' + [a[1], a[2], a[3] + (k === 0 ? dz : 0)].map(v => v.toFixed(10).padStart(17)).join('')).join('\n') +
    '\n ----------------------------------------------------------------\n' + (E ? ' Total energy in the final basis set =      ' + E + '\n' : '');
  const end = qc.lastIndexOf(' --------------------------------------------------------------\n \n        *');
  const anh = qc.slice(0, end) + ' =====================================================\n = Starting finite difference calculation for ANHAR  =\n Finite diff step size:  9.99854e-02 a.u.\n\n' +
    geo(0.05, '-76.4091234567') + geo(-0.05, '-76.4090123456') + geo(0.1, '-76.4070123456') + geo(-0.1, '-76.4068123456') +
    '      ANHARMONIC CUBIC & QUARTIC FORCE FIELD (CM**-1)\n Eta[1,1,1]=-12.3456\n' + geo(0, null) + '                    VIBRATIONAL ANHARMONIC ANALYSIS\n' + qc.slice(end);
  const qa = BM.parseFile('qchem.out', qc), da = BM.parseFile('anh.out', anh);
  check(da.frames.length === qa.frames.length && last(da.frames).energy === last(qa.frames).energy && da.vibrations.frameIndex === qa.vibrations.frameIndex,
    'Q-Chem anharmonic: the displaced geometries and the reference printed again are no steps: ' + da.frames.length + ' vs ' + qa.frames.length + ' frames');
});

// ------------------------------------------------------------------ what the level diagram and the Info tab are given
section('incomplete orbital energies', () => {
  // orbital energies that stop before the HOMO are left out, with a note, instead of breaking the level diagram
  const cut = read('qc/qchem_water_optfreq.out').replace('-19.1351  -0.9971  -0.5203  -0.3662  -0.2921\n', '-19.1351  -0.9971\n');
  const dc = BM.parseFile('cut.out', cut);
  check(dc.orbitals === null && dc.warnings.some(w => /orbital energies in this file are incomplete/.test(w)), 'incomplete orbital energies are left out with a note: ' + JSON.stringify(dc.orbitals));
});

// ------------------------------------------------------------------ validation
section('validation', () => {
  const xyz = '3\nwater\nO 0 0 0.117\nH 0 0.757 -0.469\nH 0 -0.757 -0.469\n';
  check(throws(() => BM.parseFile('big.xyz', '1\n\nH 1e999 0 0\n'), /^Could not read "big\.xyz" as XYZ: .*not a finite number/), 'a validation message names the file and the reader');
  const ds = BM.parseFile('w.xyz', xyz);
  check(ds.frames[0].Z.join() === '8,1,1', '(water)');
  const fch = read('water_freq.fchk').replace(/^(Atomic numbers\s+I\s+N=\s+3\n)\s*8/m, '$1 ****');
  const df = BM.parseFile('bad.fchk', fch);
  check(df.frames[0].Z[0] === 0 && df.warnings.some(w => /could not be read; shown as dummy atoms/.test(w)), 'an unreadable element becomes a dummy atom, with a warning: ' + df.frames[0].Z + ' ' + df.warnings);
});

for (const f of fails) console.log('  FAIL:', f);
console.log(`Real-layout reader tests: ${passes} passed, ${failures} failed`);
process.exit(failures ? 1 : 0);
