#!/usr/bin/env node
/* Outputs of jobs that are still running (issue #1): the file ends part-way through the job, often part-way through a
 * line. Every program output in tests/fixtures is cut at many places (line ends, random characters, and right after
 * energy keywords, where "SCF Done:  E(RB3LYP) =  -0." would read as an energy of 0) and parsed as text and as a
 * line array. A cut file must never produce an energy that the complete file does not contain, and an unfinished
 * job is marked as such.   Run:  node tests/test_running_jobs.js */
'use strict';
const path = require('path');
const fs = require('fs');

const ROOT = path.resolve(__dirname, '..');
for (const m of ['elements', 'bigfile', 'parsers', 'readers-struct', 'writers', 'readers-qc', 'geom', 'chem', 'analysis']) {
  require(path.join(ROOT, 'src', m + '.js'));
}
const BM = globalThis.BM;

let passes = 0, failures = 0;
const fails = [];
const check = (cond, msg) => { if (cond) passes++; else { failures++; fails.push(msg); } };

const OUTPUT_IDS = new Set(['gaussian-log', 'orca-out']);
for (const f of BM.formats) if (f.kind === 'output') OUTPUT_IDS.add(f.id);

/** Every energy a dataset holds: frames, their components, and optimisation cycles. */
function energies(ds) {
  const out = [];
  for (const f of ds.frames) {
    if (f.energy !== null && f.energy !== undefined) out.push(f.energy);
    if (f.energies) for (const v of Object.values(f.energies)) if (typeof v === 'number') out.push(v);
  }
  for (const s of ds.optSteps || []) if (s.energy !== null && s.energy !== undefined) out.push(s.energy);
  for (const p of ds.ircPath || []) if (p.energy !== null && p.energy !== undefined) out.push(p.energy);
  return out;
}
const near = (a, b) => Math.abs(a - b) <= 1e-9 * Math.max(1, Math.abs(b));

let seed = 12345;
const rnd = () => (seed = (seed * 1103515245 + 12345) % 2147483648) / 2147483648;
const ENERGY_KEYS = [/SCF Done:\s+E\([^)]*\)\s+=\s*/g, /FINAL SINGLE POINT ENERGY[^\n\d-]*/g, /extrapolated energy\s*=\s*/g,
  /TOTAL ENERGY\s*[:=]?\s*/gi, /Total Energy\s*[:=]\s*/g, /ENERGY\s*=\s*/g, /energy[:=]\s*/gi, /HEAT OF FORMATION\s*=\s*/g];

const files = [];
for (const d of ['tests/fixtures', 'tests/fixtures/qc', 'tests/fixtures/qmmm', 'tests/fixtures/review']) {
  for (const f of fs.readdirSync(path.join(ROOT, d))) if (/\.(log|out)$/.test(f)) files.push(path.join(ROOT, d, f));
}
let outputs = 0, cuts = 0, parsed = 0, running = 0, noted = 0, plottedCuts = 0;
for (const p of files) {
  const name = path.basename(p);
  const text = fs.readFileSync(p, 'utf8');
  let full;
  try { full = BM.parseFile(name, text); } catch (e) { continue; }
  if (!OUTPUT_IDS.has(full.format)) continue;
  outputs++;
  const known = energies(full).filter(Number.isFinite);
  if (!/\r?\n$/.test(text)) continue;   // the complete file itself must end with a line break for the comparison
  check(!(full.warnings || []).some(w => /middle of a line/.test(w)), name + ': a complete file has no partial-line note');
  // cut positions: line ends, positions just after energy keywords (and a few characters into the number), random
  const pos = new Set();
  const lineEnds = [];
  for (let i = text.indexOf('\n'); i >= 0; i = text.indexOf('\n', i + 1)) lineEnds.push(i + 1);
  for (let k = 1; k <= 25; k++) pos.add(lineEnds[Math.floor(lineEnds.length * k / 26)]);
  for (const re of ENERGY_KEYS) {
    for (const m of text.matchAll(re)) for (const extra of [0, 2, 4, 7]) pos.add(m.index + m[0].length + extra);
  }
  for (let k = 0; k < 30; k++) pos.add(1 + Math.floor(rnd() * (text.length - 1)));
  for (const c of pos) {
    if (!(c > 0 && c < text.length)) continue;
    const cut = text.slice(0, c);
    cuts++;
    for (const asLines of [false, true]) {
      let ds;
      try { ds = asLines ? BM.parseLines(name, cut.split(/\r\n|\r|\n/)) : BM.parseFile(name, cut); }
      catch (e) { continue; }   // too early in the file: no geometry yet, which is reported as an error
      parsed++;
      const got = energies(ds);
      const bad = got.filter(e => !Number.isFinite(e) || !known.some(k => near(e, k)));
      check(bad.length === 0, `${name} cut at ${c}${asLines ? ' (lines)' : ''}: energies not in the complete file: ${bad.slice(0, 3).join(', ')} ` +
        `after ...${JSON.stringify(cut.slice(-60))}`);
      // the plots (convergence of every job, scan and IRC profiles) use only printed energies, and do not crash
      let plotted = [], crash = '';
      try {
        const jobs = [null].concat(Array.isArray(ds.jobs) ? ds.jobs.map((_, k) => k) : []);
        for (const job of jobs) {
          const cs = BM.analysis.convergenceSeries(ds, { job });
          for (const s of cs.series) {
            if (s.key === 'energy') plotted = plotted.concat(s.values.filter(v => v !== null));
            if (s.values.some(v => v !== null && !Number.isFinite(v))) crash = crash || s.key + ' series has a non-number';
          }
        }
        for (const p of BM.analysis.profiles(ds, {})) plotted = plotted.concat(p.points.map(q => q.energy));
      } catch (e) { crash = e.message; }
      check(!crash, `${name} cut at ${c}${asLines ? ' (lines)' : ''}: plots: ${crash}`);
      const badPlot = plotted.filter(e => !Number.isFinite(e) || !known.some(k => near(e, k)));
      check(badPlot.length === 0, `${name} cut at ${c}${asLines ? ' (lines)' : ''}: plotted energies not in the complete file: ${badPlot.slice(0, 3).join(', ')}`);
      if (plotted.length) plottedCuts++;
      const unfinished = ds.meta && ds.meta.normalTermination !== true && ds.meta.normalTermination !== false;
      if (unfinished && ds.meta.running) running++;
      if (!/\n$/.test(cut) && (ds.warnings || []).some(w => /middle of a line/.test(w))) noted++;
      // (an early cut can look like another format, e.g. the input echoed at the top of an NWChem output)
      check(!OUTPUT_IDS.has(ds.format) || /\n$/.test(cut) || (ds.warnings || []).some(w => /middle of a line/.test(w)), `${name} cut at ${c}: the partial last line is noted`);
    }
  }
}
check(outputs >= 35, 'program outputs among the fixtures: ' + outputs);
check(parsed > 2000, 'cut files that could be read: ' + parsed);
check(running > 1000, 'cut files marked as not finished: ' + running);
check(noted > 1000, 'cut files with the partial-line note: ' + noted);
check(plottedCuts > 1000, 'cut files with an energy plot or profile: ' + plottedCuts);

// the example the issue describes: Gaussian, last SCF energy cut after its first characters
const g = fs.readFileSync(path.join(ROOT, 'tests', 'fixtures', 'g16_water_optfreq.log'), 'utf8');
const at = g.lastIndexOf('SCF Done:');
const head = g.slice(0, g.indexOf('=', at) + 4);   // "... SCF Done:  E(RB3LYP) =  -7"
const ds = BM.parseFile('running.log', head);
check(ds.frames.every(f => f.energy === null || f.energy < -70), 'Gaussian cut inside the last "SCF Done" number: no bogus energy (' +
  ds.frames.map(f => f.energy).join(', ') + ')');
check(ds.meta.running === true, 'Gaussian job without "Normal termination" is marked as running');
const fullG = BM.parseFile('done.log', g);
check(!fullG.meta.running && fullG.meta.normalTermination === true, 'a finished Gaussian job is not marked as running');
const c = BM.analysis.convergenceSeries(ds, {});
const eSer = c.series.find(s => s.key === 'energy');
check(!eSer || eSer.values.every(v => v === null || v < -70), 'the energy series of the running job has no bogus point');
// structure files are left alone: an XYZ file without a final line break keeps its last atom
const xyz = BM.parseFile('w.xyz', '3\nwater\nO 0 0 0.117\nH 0 0.757 -0.469\nH 0 -0.757 -0.469');
check(xyz.frames[0].Z.length === 3, 'an XYZ file without a final line break keeps its last atom');

for (const f of fails) console.log('  FAIL:', f);
console.log(`Running-job tests: ${passes} passed, ${failures} failed (${outputs} outputs, ${cuts} cuts, ${parsed} parsed)`);
process.exit(failures ? 1 : 0);
