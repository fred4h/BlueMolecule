/* Synthetic huge QM/MM outputs for the big-file tests, produced as a stream of bytes: nothing is written to disk, so
 * files of 5 or 10 GB can be read in a sandbox with little disk space. The layout follows the public appearance of a
 * Gaussian ONIOM(QM:MM) optimisation + frequency job (tools/make_big_fixture.py writes the small versions):
 *   - a QM region of a few atoms in a box of MM water, "Symbolic Z-matrix" with layers and MM types;
 *   - per step: the real-system "Input orientation", the three ONIOM sub-calculations (with the model system's own
 *     orientation table), "ONIOM: extrapolated energy", the forces table and the Berny convergence table;
 *   - "Optimization completed", the structure printed again, then a Link1 frequency job whose "Harmonic frequencies"
 *     section lists all 3N - 6 normal modes (the bulk of the file: about N^2 rows).
 *
 *   gaussianOniom({atoms, steps, qm, modes}) -> {stream(): ReadableStream<Uint8Array>, size (exact bytes),
 *                                               expected: {frames, atoms, energies, nModes, freq(m), mode(m, a, c)}}
 *   atomsForSize(bytes) -> the atom count that gives a file of about that size (with steps = 20)
 */
'use strict';

const enc = new TextEncoder();
const pad = (s, n) => { s = String(s); return s.length >= n ? s : ' '.repeat(n - s.length) + s; };
const f = (v, w, d) => pad(v.toFixed(d), w);
const GLINE = ' ' + '-'.repeat(69);
const QM = [[6, 'CT', -0.18, 0.000, 0.000, 0.000], [6, 'C', 0.58, 1.520, 0.000, 0.000], [8, 'O', -0.55, 2.150, 1.050, 0.000],
  [8, 'OH', -0.61, 2.150, -1.150, 0.000], [1, 'HO', 0.45, 3.080, -1.020, 0.000], [1, 'HC', 0.06, -0.380, 1.020, 0.000],
  [1, 'HC', 0.06, -0.380, -0.510, 0.880], [1, 'HC', 0.06, -0.380, -0.510, -0.880], [7, 'N', -0.40, 0.000, 0.000, 1.450]];
const SYM = { 1: 'H', 6: 'C', 7: 'N', 8: 'O' };

function system(N, nq) {
  const Z = new Int32Array(N), X = new Float64Array(3 * N), type = [], charge = [], layer = [];
  for (let a = 0; a < N; a++) {
    if (a < nq) {
      const q = QM[a % QM.length], shift = Math.floor(a / QM.length) * 4.0;
      Z[a] = q[0]; X[3 * a] = q[3] + shift; X[3 * a + 1] = q[4]; X[3 * a + 2] = q[5];
      type.push(q[1]); charge.push(q[2]); layer.push('H');
    } else {   // MM water on a lattice
      const w = Math.floor((a - nq) / 3), r = (a - nq) % 3, side = 40;
      const ox = 6 + 3.1 * (w % side), oy = -60 + 3.1 * (Math.floor(w / side) % side), oz = -60 + 3.1 * Math.floor(w / (side * side));
      const d = [[0, 0, 0], [0.757, 0.586, 0], [-0.757, 0.586, 0]][r];
      Z[a] = r ? 1 : 8; X[3 * a] = ox + d[0]; X[3 * a + 1] = oy + d[1]; X[3 * a + 2] = oz + d[2];
      type.push(r ? 'HW' : 'OW'); charge.push(r ? 0.417 : -0.834); layer.push('L');
    }
  }
  return { Z, X, type, charge, layer };
}

/** Values in hundredths for mode m, atom a, coordinate c: deterministic, in [-100, 100] (Gaussian prints F7.2). */
const modeValue = (m, a, c) => ((m * 7919 + a * 104729 + c * 31337) % 201 + 201) % 201 - 100;
const freqOf = m => (m === 0 ? -35.1234 : 12.5 + m * 0.137);

function gaussianOniom(opts) {
  const N = opts.atoms, S = Math.max(2, opts.steps || 20), nq = Math.min(opts.qm || 9, N);
  const nModes = opts.modes === undefined ? 3 * N - 6 : opts.modes;
  const sys = system(N, nq);
  const E = k => (k < S - 1 ? -268.4 + 0.012 * Math.exp(-0.25 * k) : -268.4);
  const route = '#p oniom(b3lyp/6-31g(d):amber)=embedcharge opt freq nosymm geom=connectivity';
  const F72 = [];
  for (let v = -100; v <= 100; v++) F72.push(pad((v / 100).toFixed(2), 7));
  const PFX = Array.from({ length: N }, (_, a) => pad(a + 1, 6) + pad(sys.Z[a], 4));
  const step = k => { const s = 0.03 * Math.exp(-0.3 * k) * (k < S - 1 ? 1 : 0); return a => (a < nq ? s * Math.sin(1.3 * a + 0.5 * k) : 0); };

  function* header() {
    yield [' Entering Gaussian System, Link 0=g16', ' Input=job.gjf', ' Output=job.log',
      ' Copyright (c) 1988-2019, Gaussian, Inc.  All Rights Reserved.', ' This is part of the Gaussian(R) 16 program.',
      ' ******************************************', ' Gaussian 16:  ES64L-G16RevC.01  3-Jul-2019', '                16-Nov-2026 ',
      ' ******************************************', ' %chk=job.chk', ' %mem=64GB',
      ' ' + '-'.repeat(route.length + 1), ' ' + route, ' ' + '-'.repeat(route.length + 1),
      ' 1/14=-1,18=20,19=15,38=1/1,3;', ' 2/9=110,12=2,17=6,18=5,40=1/2;',
      ' -------------------------', ' synthetic QM/MM test box', ' -------------------------', ' Symbolic Z-matrix:',
      ' Charge =  0 Multiplicity = 1 for low   level calculation on real  system.',
      ' Charge =  0 Multiplicity = 1 for high  level calculation on model system.',
      ' Charge =  0 Multiplicity = 1 for low   level calculation on model system.'].join('\n') + '\n';
    let out = [];
    for (let a = 0; a < N; a++) {
      const spec = SYM[sys.Z[a]] + '-' + sys.type[a] + '-' + sys.charge[a].toFixed(4);
      out.push(' ' + spec.padEnd(28) + pad(0, 3) + f(sys.X[3 * a], 14, 6) + f(sys.X[3 * a + 1], 14, 6) + f(sys.X[3 * a + 2], 14, 6) + ' ' + sys.layer[a]);
      if (out.length >= 20000) { yield out.join('\n') + '\n'; out = []; }
    }
    if (out.length) yield out.join('\n') + '\n';
    yield ' \n NAtoms=' + pad(N, 8) + ' NQM=' + pad(nq, 9) + ' NQMF=       0 NMMI=      0 NMMIF=      0\n' +
      '                NMic=' + pad(N - nq, 9) + ' NMicF=      0.\n';
  }
  function* orient(k, real) {
    const d = step(k);
    yield ['                          Input orientation:                          ', GLINE,
      ' Center     Atomic      Atomic             Coordinates (Angstroms)', ' Number     Number       Type             X           Y           Z', GLINE].join('\n') + '\n';
    const n = real ? N : nq;
    let out = [];
    for (let a = 0; a < n; a++) {
      const s = d(a);
      out.push(pad(a + 1, 7) + pad(sys.Z[a], 11) + pad(0, 12) + '    ' + f(sys.X[3 * a] + s, 12, 6) + f(sys.X[3 * a + 1] - s, 12, 6) + f(sys.X[3 * a + 2] + 0.5 * s, 12, 6));
      if (out.length >= 20000) { yield out.join('\n') + '\n'; out = []; }
    }
    out.push(GLINE);
    yield out.join('\n') + '\n';
  }
  function* forces(k) {
    yield [' ***** Axes restored to original set *****', ' ' + '-'.repeat(67), ' Center     Atomic                   Forces (Hartrees/Bohr)',
      ' Number     Number              X              Y              Z', ' ' + '-'.repeat(67)].join('\n') + '\n';
    const s = Math.exp(-k / 6);
    let out = [];
    for (let a = 0; a < N; a++) {
      out.push(pad(a + 1, 7) + pad(sys.Z[a], 9) + '       ' + f(0.01 * s * Math.sin(a), 15, 9) + f(-0.01 * s * Math.cos(a), 15, 9) + f(0.002 * s, 15, 9));
      if (out.length >= 20000) { yield out.join('\n') + '\n'; out = []; }
    }
    out.push(' ' + '-'.repeat(67), ' Cartesian Forces:  Max     ' + (0.01 * s).toFixed(9) + ' RMS     ' + (0.004 * s).toFixed(9));
    yield out.join('\n') + '\n';
  }
  function oniomEnergies(k) {
    const e = E(k), lr = -0.0281 + 0.002 * Math.exp(-k), lm = 0.0135 + 0.001 * Math.exp(-k), hm = e - lr + lm;
    return { e, lr, lm, hm };
  }
  function* subcalcs(k) {
    const { e, lr, lm, hm } = oniomEnergies(k);
    yield ' ONIOM: generating point  1 -- low level on model system.\n Energy= ' + f(lm, 20, 12) + ' NIter=   0.\n' +
      ' ONIOM: generating point  2 -- high level on model system.\n';
    yield* orient(k, false);
    yield ' SCF Done:  E(RB3LYP) =  ' + hm.toFixed(12) + '     A.U. after  11 cycles\n' +
      ' ONIOM: generating point  3 -- low level on real system.\n Energy= ' + f(lr, 20, 12) + ' NIter=   0.\n' +
      ' ONIOM: calculating energy.\n' +
      ' ONIOM: gridpoint  1 method:  low   system:  model energy: ' + f(lm, 20, 12) + '\n' +
      ' ONIOM: gridpoint  2 method:  high  system:  model energy: ' + f(hm, 20, 12) + '\n' +
      ' ONIOM: gridpoint  3 method:  low   system:  real  energy: ' + f(lr, 20, 12) + '\n' +
      ' ONIOM: extrapolated energy = ' + f(e, 20, 12) + '\n';
  }
  function conv(k) {
    const s = Math.exp(-k / 6), last = k === S - 1;
    const vals = last ? [[0.000012, 0.00045], [0.000004, 0.0003], [0.00015, 0.0018], [0.00006, 0.0012]]
      : [[0.02 * s, 0.00045], [0.008 * s, 0.0003], [0.05 * s, 0.0018], [0.02 * s, 0.0012]];
    const names = ['Maximum Force', 'RMS     Force', 'Maximum Displacement', 'RMS     Displacement'];
    return [' GradGradGradGradGradGradGradGradGradGradGradGradGradGradGradGradGradGrad', ' Berny optimization.',
      ' Step number' + pad(k + 1, 4) + ' out of a maximum of' + pad(100, 5), '         Item               Value     Threshold  Converged?']
      .concat(names.map((nm, q) => ' ' + nm.padEnd(20) + f(vals[q][0], 11, 6) + f(vals[q][1], 13, 6) + '     ' + (vals[q][0] < vals[q][1] ? 'YES' : 'NO ')),
        [' Predicted change in Energy=-' + (1.234 * s).toFixed(6) + 'D-04'], last ? [' Optimization completed.', '    -- Stationary point found.'] : [],
        [' GradGradGradGradGradGradGradGradGradGradGradGradGradGradGradGradGradGrad']).join('\n') + '\n';
  }
  function* freq() {
    yield ' Harmonic frequencies (cm**-1), IR intensities (KM/Mole), Raman scattering\n activities (A**4/AMU), depolarization ratios for plane and unpolarized\n' +
      ' incident light, reduced masses (AMU), force constants (mDyne/A),\n and normal coordinates:\n';
    for (let m0 = 0; m0 < nModes; m0 += 3) {
      const ms = [m0, m0 + 1, m0 + 2].filter(m => m < nModes);
      const col = (fn, w, first) => ms.map((m, q) => pad(fn(m), q ? 23 : first)).join('');
      const head = [col(m => m + 1, 0, 22), col(() => 'A', 0, 22), ' Frequencies --' + col(m => freqOf(m).toFixed(4), 0, 12),
        ' Red. masses --' + col(m => (1 + (m % 7) * 0.5).toFixed(4), 0, 12), ' Frc consts  --' + col(m => (0.001 * m).toFixed(4), 0, 12),
        ' IR Inten    --' + col(m => (m % 13 * 1.7).toFixed(4), 0, 12),
        '  Atom  AN' + ms.map(() => '      X      Y      Z  ').join('').replace(/\s+$/, '')];
      let out = head.join('\n') + '\n';
      const rows = [];
      for (let a = 0; a < N; a++) {
        let r = PFX[a];
        for (const m of ms) r += '  ' + F72[modeValue(m, a, 0) + 100] + F72[modeValue(m, a, 1) + 100] + F72[modeValue(m, a, 2) + 100];
        rows.push(r);
      }
      out += rows.join('\n') + '\n';
      yield out;
    }
    yield [' -------------------', ' - Thermochemistry -', ' -------------------', ' Temperature   298.150 Kelvin.  Pressure   1.00000 Atm.',
      ' Zero-point correction=                           1.234567 (Hartree/Particle)',
      ' Thermal correction to Energy=                    1.345678', ' Thermal correction to Enthalpy=                  1.346622',
      ' Thermal correction to Gibbs Free Energy=         1.100000', ' Sum of electronic and zero-point Energies=           -267.165433',
      ' Sum of electronic and thermal Energies=              -267.054322', ' Sum of electronic and thermal Enthalpies=            -267.053378',
      ' Sum of electronic and thermal Free Energies=         -267.300000'].join('\n') + '\n';
  }
  function* all() {
    yield* header();
    for (let k = 0; k < S; k++) {
      yield* orient(k, true);
      yield* subcalcs(k);
      yield* forces(k);
      yield conv(k);
    }
    yield* orient(S - 1, true);   // the optimised structure printed again
    const r2 = '#p oniom(b3lyp/6-31g(d):amber)=embedcharge freq geom=allcheck guess=tcheck';
    yield ' Normal termination of Gaussian 16 at Mon Nov 16 12:00:00 2026.\n Link1:  Proceeding to internal job step number  2.\n' +
      ' ' + '-'.repeat(r2.length + 1) + '\n ' + r2 + '\n ' + '-'.repeat(r2.length + 1) + '\n';
    yield* orient(S - 1, true);
    yield* subcalcs(S - 1);
    if (nModes) yield* freq();
    yield ' Normal termination of Gaussian 16 at Mon Nov 16 13:00:00 2026.\n';
  }
  let size = null;
  const measure = () => { let n = 0; for (const s of all()) n += s.length; return n; };   // (ASCII: one byte per character)
  return {
    get size() { if (size === null) size = measure(); return size; },
    text: () => { const parts = []; for (const s of all()) parts.push(s); return parts.join(''); },
    stream() {
      const it = all();
      let buf = [], n = 0;
      return new ReadableStream({
        pull(ctl) {
          for (;;) {
            const r = it.next();
            if (r.done) { if (buf.length) { ctl.enqueue(enc.encode(buf.join(''))); buf = []; } ctl.close(); return; }
            buf.push(r.value); n += r.value.length;
            if (n >= 4 * 1048576) { ctl.enqueue(enc.encode(buf.join(''))); buf = []; n = 0; return; }
          }
        }
      });
    },
    expected: { frames: S, atoms: N, nModes, qm: nq, energies: Array.from({ length: S }, (_, k) => E(k)), freq: freqOf, mode: (m, a, c) => modeValue(m, a, c) / 100,
      xyz0: [0, 1, 2].flatMap(a => { const d = step(0)(a); return [+(sys.X[3 * a] + d).toFixed(6), +(sys.X[3 * a + 1] - d).toFixed(6), +(sys.X[3 * a + 2] + 0.5 * d).toFixed(6)]; }) }
  };
}
/** Atoms for a file of about `bytes` (20 optimisation steps; the frequency rows are ~80 bytes each, N^2 of them). */
function atomsForSize(bytes, steps) {
  const S = steps || 20;
  let N = Math.floor(Math.sqrt(bytes / 81));
  for (let it = 0; it < 30; it++) {   // N^2 * 80.3 (modes) + S * N * ~230 (orientation, forces, model system) + header
    const est = N * N * 80.3 + (S + 2) * N * 235 + N * 74;
    N = Math.max(10, Math.round(N * Math.sqrt(bytes / est)));
  }
  return N;
}
module.exports = { gaussianOniom, atomsForSize, modeValue };
