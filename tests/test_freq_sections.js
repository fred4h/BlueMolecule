#!/usr/bin/env node
/* Gaussian frequency sections of many blocks, as real outputs print them, read in full and condensed.
 *   node tests/test_freq_sections.js
 * The column layout is that of real Gaussian 09/16 outputs (Atom AN rows: I6, I4, then 2X, 3F7.2 per mode; labels in
 * 15 columns and values in F12.4 + 11X). Variants: imaginary modes, a value too wide for its column printed as asterisks
 * right after the dashes (huge IR intensities of charged QM/MM systems), an unknown line among the labels, atom numbers
 * too wide for their column, a damaged row (the reader says where it stopped), and a second, shorter section. */
'use strict';
const path = require('path');
const ROOT = path.join(__dirname, '..');
for (const f of ['elements', 'parsers', 'readers-qc', 'condense', 'bigfile']) require(path.join(ROOT, 'src', f + '.js'));
const BM = globalThis.BM;

let pass = 0, fail = 0;
const failures = [];
const ok = (c, m) => { if (c) pass++; else { fail++; failures.push(m); } };

const N = 7, NM = 3 * N - 6;   // 15 vibrations: 5 blocks of 3
const Z = [6, 6, 8, 1, 1, 1, 1];
const pad = (s, w) => String(s).padStart(w);
const disp = (m, a, c) => Math.round(100 * Math.sin(1.3 * m + 0.7 * a + 2.1 * c) * 0.6) / 100;   // 2 decimals, |x| <= 0.6
const freq = m => (m < 2 ? -(54.0897 / (m + 1)) : 100 + 37.5 * m);

function log(opts) {
  opts = opts || {};
  const L = [' Entering Gaussian System, Link 0=g16', ' Gaussian 16:  ES64L-G16RevC.01  3-Jul-2019', ' ******************************************',
    ' ---------------------------------', ' #p b3lyp/6-31g(d) freq', ' ---------------------------------', ' ',
    ' title', ' ', '                          Input orientation:',
    ' ---------------------------------------------------------------------',
    ' Center     Atomic      Atomic             Coordinates (Angstroms)',
    ' Number     Number       Type             X           Y           Z',
    ' ---------------------------------------------------------------------'];
  for (let a = 0; a < N; a++) L.push(pad(a + 1, 7) + pad(Z[a], 11) + pad(0, 12) + '    ' + [0, 1, 2].map(c => pad((a * 0.9 + c * 0.31).toFixed(6), 12)).join(''));
  L.push(' ---------------------------------------------------------------------',
    ' SCF Done:  E(RB3LYP) =  -154.123456789     A.U. after   10 cycles',
    ' Low frequencies ---  -54.0897  -27.0449    0.0011    0.0012    0.0013    0.0014',
    ' Harmonic frequencies (cm**-1), IR intensities (KM/Mole), Raman scattering',
    ' activities (A**4/AMU), depolarization ratios for plane and unpolarized',
    ' incident light, reduced masses (AMU), force constants (mDyne/A),',
    ' and normal coordinates:');
  const nm = opts.modes || NM;
  for (let m0 = 0; m0 < nm; m0 += 3) {
    const ms = [m0, m0 + 1, m0 + 2].filter(m => m < nm), b = m0 / 3;
    const col = (f, first) => ms.map((m, q) => pad(f(m), q ? 23 : first)).join('');
    const val = f => ' ' + ms.map((m, q) => (q ? ' '.repeat(11) : '') + pad(f(m), 12)).join('').trimEnd().slice(1);
    L.push(col(m => m + 1, 22), col(() => 'A', 22));
    L.push(' Frequencies --' + val(m => freq(m).toFixed(4)));
    L.push(' Red. masses --' + val(m => (1 + m * 0.25).toFixed(4)));
    L.push(' Frc consts  --' + val(m => (0.01 * m).toFixed(4)));
    if (opts.stars && b === 0) L.push(' IR Inten    --' + ms.map((m, q) => (q ? ' '.repeat(11) : '') + (q === 0 ? '************' : pad((m * 1.5).toFixed(4), 12))).join(''));
    else L.push(' IR Inten    --' + val(m => (m * 1.5).toFixed(4)));
    if (opts.unknown && b === 1) L.push(' (some line this reader does not know)');
    if (opts.noRowsLast && m0 + 3 >= nm) continue;   // a block printed without its atom rows
    L.push('  Atom  AN' + ms.map(() => '      X      Y      Z  ').join('').trimEnd());
    for (let a = 0; a < N; a++) {
      let r = opts.bigAtoms && a >= 4 ? '******' + pad(Z[a], 4) : pad(a + 1, 6) + pad(Z[a], 4);
      for (const m of ms) r += '  ' + [0, 1, 2].map(c => pad(disp(m, a, c).toFixed(2), 7)).join('');
      if (opts.damaged && b === 2 && a === 3) r = r.slice(0, 30) + ' x' + r.slice(32);
      L.push(r);
    }
  }
  if (opts.again) {   // a second, shorter section for the same geometry
    L.push(' Harmonic frequencies (cm**-1), IR intensities (KM/Mole), Raman scattering', ' and normal coordinates:',
      pad(1, 22), pad('A', 22), ' Frequencies --' + pad(freq(0).toFixed(4), 12), '  Atom  AN      X      Y      Z');
    for (let a = 0; a < N; a++) L.push(pad(a + 1, 6) + pad(Z[a], 4) + '  ' + [0, 1, 2].map(c => pad(disp(0, a, c).toFixed(2), 7)).join(''));
  }
  L.push('', ' -------------------', ' - Thermochemistry -', ' -------------------', ' Temperature   298.150 Kelvin.  Pressure   1.00000 Atm.',
    ' Zero-point correction=                           0.081234 (Hartree/Particle)', ' Normal termination of Gaussian 16 at Tue Nov 17 12:00:00 2026.', '');
  return L.join('\n');
}

async function both(name, text) {
  const full = BM.parseFile(name, text);
  const lines = await BM.bigfile.readCondensed(new Blob([text]), () => BM.condense.create('gaussian-log', { load: { steps: { mode: 'all' }, modes: { mode: 'all' } } }));
  return { full, cond: BM.parseLines(name, lines) };
}
const vec = (v, m) => Array.from(v.modes[m] || []);
function sameVibs(r, label) {
  const a = r.full.vibrations, b = r.cond.vibrations;
  ok(a && b && a.freqs.length === b.freqs.length, label + ': same number of vibrations read in full and condensed (' + (a && a.freqs.length) + ', ' + (b && b.freqs.length) + ')');
  if (!a || !b) return;
  let same = true;
  for (let m = 0; m < a.freqs.length; m++) if (a.freqs[m] !== b.freqs[m] || vec(a, m).join() !== vec(b, m).join()) same = false;
  ok(same, label + ': same frequencies and displacements in full and condensed');
  ok(JSON.stringify(r.full.warnings) === JSON.stringify(r.cond.warnings), label + ': same reader notes: ' + JSON.stringify([r.full.warnings, r.cond.warnings]));
}
const expectDisp = (v, m) => { const out = []; for (let a = 0; a < N; a++) for (let c = 0; c < 3; c++) out.push(disp(m, a, c)); return out; };
const allModes = (v, upto) => { for (let m = 0; m < upto; m++) if (vec(v, m).join() !== expectDisp(v, m).join()) return false; return true; };

module.exports = { log, NM, disp };
if (require.main === module) (async () => {
  const plain = await both('plain.log', log());
  for (const [k, r] of [['full', plain.full], ['condensed', plain.cond]]) {
    const v = r.vibrations;
    ok(v && v.freqs.length === NM, `plain (${k}): all ${NM} vibrations, 5 blocks: ${v && v.freqs.length}`);
    ok(v && v.freqs[0] === -54.0897 && v.freqs[1] === -27.0449, `plain (${k}): imaginary modes as negative frequencies`);
    ok(v && allModes(v, NM), `plain (${k}): every displacement of every block`);
    ok(r.warnings.length === 0, `plain (${k}): no reader notes: ` + r.warnings.join(' | '));
  }
  sameVibs(plain, 'plain');

  const stars = await both('stars.log', log({ stars: true }));
  ok(stars.full.vibrations.freqs.length === NM && allModes(stars.full.vibrations, NM), 'an IR intensity printed as asterisks after the dashes: every block is still read');
  ok(Number.isNaN(stars.full.vibrations.ir[0]) && stars.full.vibrations.ir[1] === 1.5, 'the asterisks are no value; the next column is read');
  sameVibs(stars, 'asterisks');

  const unknown = await both('unknown.log', log({ unknown: true }));
  ok(unknown.full.vibrations.freqs.length === NM && allModes(unknown.full.vibrations, NM), 'an unknown line among the labels of a block is skipped');
  sameVibs(unknown, 'unknown line');

  const big = await both('bigatoms.log', log({ bigAtoms: true }));
  ok(big.full.vibrations.freqs.length === NM && allModes(big.full.vibrations, NM), 'atom numbers printed as asterisks: the next atom after the previous row');
  sameVibs(big, 'atom numbers as asterisks');

  const damaged = await both('damaged.log', log({ damaged: true }));
  const dv = damaged.full.vibrations;
  ok(dv.freqs.length === 9, 'a damaged row in the third block: the three first blocks are read (' + dv.freqs.length + ')');
  ok(allModes(dv, 6), 'the complete blocks keep their displacements');
  ok(damaged.full.warnings.some(w => /stopped early, after 9 vibrations/.test(w) && /atom row of an unexpected form/.test(w)), 'the reader notes say where it stopped: ' + damaged.full.warnings.join(' | '));
  sameVibs(damaged, 'damaged row');

  const again = await both('again.log', log({ again: true }));
  ok(again.full.vibrations.freqs.length === NM && allModes(again.full.vibrations, NM), 'a second, shorter section for the same geometry does not replace the full one');
  sameVibs(again, 'second section');

  const norows = await both('norows.log', log({ noRowsLast: true }));
  ok(norows.full.vibrations.freqs.length === NM && allModes(norows.full.vibrations, NM - 3) && vec(norows.full.vibrations, NM - 1).every(x => x === 0), 'a last block without atom rows: its vibrations are listed, without displacements');
  sameVibs(norows, 'block without rows');

  // a cut file (a job still running): the blocks before the cut, no note about a stop
  const text = log(), cutAt = text.indexOf(' Frequencies --', text.indexOf(' Frequencies --') + 1) + 5;
  const cut = await both('cut.log', text.slice(0, cutAt));
  ok(cut.full.vibrations && cut.full.vibrations.freqs.length === 3 && !cut.full.warnings.some(w => /stopped early/.test(w)), 'a file cut inside the second block: the first block, no stop note');
  sameVibs(cut, 'cut');

  console.log(`\nFrequency-section tests: ${pass} passed, ${fail} failed`);
  if (fail) { for (const f of failures) console.log('  FAIL ' + f); process.exit(1); }
})().catch(e => { console.error(e); process.exit(1); });
