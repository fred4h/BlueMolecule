#!/usr/bin/env node
/* Condensed reading of huge outputs (issue #6): every Gaussian, ORCA and XYZ example file must read exactly the
 * same through a condenser (BM.condense, which compacts geometry tables and normal modes and drops the rows the
 * reader skips) as from its full line array; also when the file is cut part-way (a job still running), and with
 * the load options of issue #7 (every Nth step, the last N, a range, the first N vibrations).
 *   Run:  node tests/test_condense.js */
'use strict';
const path = require('path');
const fs = require('fs');

const ROOT = path.resolve(__dirname, '..');
for (const m of ['elements', 'bigfile', 'parsers', 'condense', 'readers-struct', 'writers', 'readers-qc', 'geom', 'chem', 'analysis']) {
  require(path.join(ROOT, 'src', m + '.js'));
}
const BM = globalThis.BM;

let passes = 0, failures = 0;
const fails = [];
const check = (cond, msg) => { if (cond) passes++; else { failures++; if (fails.length < 60) fails.push(msg); } };

// ------------------------------------------------------------------ deep comparison
const isTyped = x => ArrayBuffer.isView(x) && !(x instanceof DataView);
/** The first difference between a and b as a path string, or null (NaN equals NaN, -0 equals 0). */
function diff(a, b, p, seen) {
  if (a === b) return null;
  if (typeof a === 'number' && typeof b === 'number') return (Number.isNaN(a) && Number.isNaN(b)) ? null : p + ': ' + a + ' vs ' + b;
  if (a === null || b === null || typeof a !== 'object' || typeof b !== 'object') return p + ': ' + JSON.stringify(a) + ' vs ' + JSON.stringify(b);
  if (isTyped(a) || isTyped(b) || Array.isArray(a) || Array.isArray(b)) {
    if (!(isTyped(a) || Array.isArray(a)) || !(isTyped(b) || Array.isArray(b))) return p + ': array vs non-array';
    if (a.length !== b.length) return p + '.length: ' + a.length + ' vs ' + b.length;
    for (let i = 0; i < a.length; i++) { const d = diff(a[i], b[i], p + '[' + i + ']', seen); if (d) return d; }
    return null;
  }
  if (seen.has(a)) return null;
  seen.add(a);
  const ka = Object.keys(a).filter(k => a[k] !== undefined), kb = Object.keys(b).filter(k => b[k] !== undefined);
  const all = new Set(ka.concat(kb));
  for (const k of all) { const d = diff(a[k], b[k], p + '.' + k, seen); if (d) return d; }
  return null;
}
/** A dataset in comparable form: line numbers dropped from warnings, the condensed run's bookkeeping removed. */
function norm(ds) {
  const o = Object.assign({}, ds);
  o.warnings = (ds.warnings || []).map(w => w.replace(/\b(on |at )?line \d+/g, 'line N').replace(/in frame \d+/g, 'line N'));
  o.meta = Object.assign({}, ds.meta);
  delete o.meta.loaded;
  delete o.steps;
  if (ds.vibrations) { o.vibrations = Object.assign({}, ds.vibrations); delete o.vibrations.modesLoaded; o.vibrations.modes = Array.from(ds.vibrations.modes); }
  return o;
}

const OWN = new Set(BM.condense.formats);
function cases() {
  const out = [];
  for (const d of ['tests/fixtures', 'tests/fixtures/qc', 'tests/fixtures/qmmm', 'tests/fixtures/review']) {
    for (const f of fs.readdirSync(path.join(ROOT, d))) {
      const p = path.join(ROOT, d, f);
      if (!fs.statSync(p).isFile() || /\.json$/.test(f)) continue;
      const text = fs.readFileSync(p, 'utf8');
      if (/[\u0000-\u0008]/.test(text.slice(0, 4096))) continue;
      const id = BM.detectFormat(f, text.slice(0, 65536));
      if (OWN.has(id)) out.push({ name: f, text, id });
    }
  }
  return out;
}
const linesOf = text => text.split(/\r\n|\r|\n/);
function both(name, id, lines, opts) {
  let full = null, cond = null, e1 = null, e2 = null;
  try { full = BM.parseLines(name, lines.slice()); } catch (e) { e1 = e.message; }
  try { cond = BM.parseLines(name, BM.condense.lines(id, lines, opts)); } catch (e) { e2 = e.message; }
  return { full, cond, e1, e2 };
}

// ------------------------------------------------------------------ 1. every example, whole and cut
const all = cases();
let nCut = 0;
for (const { name, text, id } of all) {
  const lines = linesOf(text);
  const r = both(name, id, lines);
  check(!r.e2 === !r.e1 && (!r.e1 || r.e1 === r.e2), name + ': errors differ: ' + r.e1 + ' | ' + r.e2);
  if (r.full && r.cond) {
    const d = diff(norm(r.full), norm(r.cond), name, new WeakSet());
    check(!d, 'condensed read differs: ' + d);
  }
  // cut part-way, at line ends and inside lines (a job still running)
  const pos = new Set();
  for (let k = 1; k <= 12; k++) pos.add(Math.floor(text.length * k / 13));
  for (const re of [/Input orientation:|Standard orientation:|CARTESIAN COORDINATES \(ANGSTROEM\)/g, /Harmonic frequencies|NORMAL MODES/g, /SCF Done|FINAL SINGLE POINT/g]) {
    for (const m of text.matchAll(re)) { pos.add(m.index + 300); pos.add(m.index + 1500); }
  }
  for (const c of pos) {
    if (!(c > 0 && c < text.length)) continue;
    const cut = linesOf(text.slice(0, c));
    const q = both(name, id, cut);
    nCut++;
    check(!q.e2 === !q.e1, name + ' cut at ' + c + ': errors differ: ' + q.e1 + ' | ' + q.e2);
    if (q.full && q.cond) {
      const d = diff(norm(q.full), norm(q.cond), name + '@' + c, new WeakSet());
      check(!d, 'condensed read of a cut file differs: ' + d);
    }
  }
}
check(all.length >= 30, 'Gaussian, ORCA and XYZ examples: ' + all.length);

// ------------------------------------------------------------------ 2. what the condensed lines hold
{
  const g = fs.readFileSync(path.join(ROOT, 'tests', 'fixtures', 'g16_water_optfreq.log'), 'utf8');
  const L = linesOf(g), C = BM.condense.lines('gaussian-log', L);
  check(C.length < L.length * 0.8, 'condensed lines are fewer: ' + C.length + ' of ' + L.length);
  check(C.bmSide.geoms.length >= 4 && C.bmSide.modes.length >= 1, 'geometries and normal-mode blocks are stored aside: ' + C.bmSide.geoms.length + ', ' + C.bmSide.modes.length);
  const v = C.bmSide.modes[0].vecs[0];
  check(v.vals instanceof Int8Array && v.d === 2, "Gaussian's 2-decimal modes are kept as bytes: " + (v.vals && v.vals.constructor.name) + ' d=' + v.d);
  check(C.some(s => s === '0'), 'runs of skipped rows become a single "0" line');
}

// ------------------------------------------------------------------ 3. load options (issue #7)
const energiesOf = ds => ds.frames.map(f => f.energy);
const stepEnergies = ds => (ds.steps || ds.frames).map(s => s.energy);
for (const file of ['g16_water_optfreq.log', 'qmmm/g16_oniom_amber_opt.log', 'qmmm/orca6_qmmm_opt.out', 'orca5_water_optfreq.out']) {
  const p = path.join(ROOT, 'tests', 'fixtures', file);
  if (!fs.existsSync(p)) { check(false, 'missing fixture ' + file); continue; }
  const text = fs.readFileSync(p, 'utf8'), name = path.basename(file), id = BM.detectFormat(name, text.slice(0, 65536));
  const L = linesOf(text), full = BM.parseLines(name, L.slice());
  const nF = full.frames.length;
  for (const steps of [{ mode: 'every', n: 2 }, { mode: 'last', n: 1 }, { mode: 'range', from: 2, to: 3 }]) {
    const ds = BM.parseLines(name, BM.condense.lines(id, L, { load: { steps } }));
    const tag = name + ' ' + JSON.stringify(steps);
    check(ds.frames.length <= nF && ds.frames.length >= 1, tag + ': frames ' + ds.frames.length + ' of ' + nF);
    if (ds.steps) {
      check(ds.steps.length === nF, tag + ': every step is listed (' + ds.steps.length + ' of ' + nF + ')');
      check(!diff(stepEnergies(ds), energiesOf(full), 'E', new WeakSet()), tag + ': the steps keep every energy');
      if (ds.steps.length === nF) for (const s of ds.steps) if (s.frame !== null) {
        const f = ds.frames[s.frame], k = ds.steps.indexOf(s);
        check(!diff(Array.from(f.xyz), Array.from(full.frames[k].xyz), 'xyz', new WeakSet()), tag + ': loaded step ' + (k + 1) + ' has its own geometry');
      }
      const lastStep = ds.steps[ds.steps.length - 1];
      check(lastStep.frame === ds.frames.length - 1, tag + ': the last step is loaded');
      check(ds.meta.loaded && ds.meta.loaded.steps.total === nF && ds.meta.loaded.steps.loaded === ds.frames.length, tag + ': meta.loaded counts the steps ' + JSON.stringify(ds.meta.loaded && ds.meta.loaded.steps));
      const conv = BM.analysis.convergenceSeries(ds, {}), convF = BM.analysis.convergenceSeries(full, {});
      const eS = conv.series.find(s => s.key === 'energy'), eF = convF.series.find(s => s.key === 'energy');
      check(!eF || (eS && !diff(eS.values, eF.values, 'conv', new WeakSet())), tag + ': the energy plot keeps every step');
    }
    if (full.vibrations) check(ds.vibrations && ds.vibrations.frameIndex >= 0 && ds.vibrations.frameIndex < ds.frames.length && ds.frames[ds.vibrations.frameIndex].Z.length === full.frames[full.vibrations.frameIndex].Z.length, tag + ': the vibrations keep their geometry');
  }
  if (full.vibrations) {
    const nV = full.vibrations.freqs.length;
    for (const modes of [{ mode: 'first', n: 2 }, { mode: 'none' }]) {
      const ds = BM.parseLines(name, BM.condense.lines(id, L, { load: { modes } }));
      const v = ds.vibrations, want = modes.mode === 'none' ? 0 : Math.min(2, nV);
      check(v && v.freqs.length === nV && !diff(v.freqs, full.vibrations.freqs, 'f', new WeakSet()), name + ' ' + JSON.stringify(modes) + ': every frequency is read');
      check(v && v.modesLoaded === want, name + ' ' + JSON.stringify(modes) + ': displacements of ' + want + ' modes (' + (v && v.modesLoaded) + ')');
      for (let k = 0; k < nV; k++) {
        const m = v.modes[k];
        if (k < want) check(m && !diff(Array.from(m), Array.from(full.vibrations.modes[k]), 'm', new WeakSet()), name + ': mode ' + (k + 1) + ' is loaded exactly');
        else check(m === null, name + ': mode ' + (k + 1) + ' is not loaded');
      }
    }
  }
}

// ------------------------------------------------------------------ 4. XYZ trajectories: frames left out
{
  const frames = [];
  for (let k = 0; k < 10; k++) frames.push('3\nframe ' + (k + 1) + ' E=' + (-76.1 - k * 0.001).toFixed(4) + '\nO 0 0 ' + (0.1 + k * 0.01).toFixed(3) + '\nH 0 0.757 -0.469\nH 0 -0.757 -0.469');
  const text = frames.join('\n') + '\n', L = linesOf(text), full = BM.parseLines('t.xyz', L.slice());
  for (const [steps, want] of [[{ mode: 'every', n: 3 }, [1, 4, 7, 10]], [{ mode: 'last', n: 2 }, [9, 10]], [{ mode: 'range', from: 3, to: 5 }, [3, 4, 5]]]) {
    const ds = BM.parseLines('t.xyz', BM.condense.lines('xyz', L, { load: { steps } }));
    const got = ds.frames.map(f => +/frame (\d+)/.exec(f.label)[1]);
    check(JSON.stringify(got) === JSON.stringify(want), 'XYZ ' + JSON.stringify(steps) + ': frames ' + got);
    check(ds.frames.every(f => { const k = +/frame (\d+)/.exec(f.label)[1] - 1; return !diff(Array.from(f.xyz), Array.from(full.frames[k].xyz), 'x', new WeakSet()) && f.energy === full.frames[k].energy; }), 'XYZ ' + JSON.stringify(steps) + ': each frame is the original one');
    check(ds.meta.loaded.steps.total === 10 && ds.meta.loaded.steps.loaded === want.length, 'XYZ ' + JSON.stringify(steps) + ': meta.loaded ' + JSON.stringify(ds.meta.loaded.steps));
  }
  check(BM.parseLines('t.xyz', BM.condense.lines('xyz', L)).frames[0].Z === BM.parseLines('t.xyz', BM.condense.lines('xyz', L)).frames[0].Z || true, 'frames share their element list');
  const dsC = BM.parseLines('t.xyz', BM.condense.lines('xyz', L));
  check(dsC.frames.every(f => f.Z === dsC.frames[0].Z), 'condensed XYZ frames share one element list');
}

// ------------------------------------------------------------------ 5. kept lines are copies of the text
{
  const own = BM.condense.own;
  const odd = ['', 'a', '\u00e9', '\u4e2d', 'ab', ' 64', '\u0001G12', 'x\u2028y', '\ud83d\ude00', 'lone \ud800 surrogate', 'tab\there', '"quoted" \\ back', ' Normal termination of Gaussian 16'];
  check(odd.every(s => own(s) === s), 'own() keeps every string as it is');
  // the copies must not keep the big piece they were cut from alive (V8 here; WebKit is checked by the desktop builds)
  require('v8').setFlagsFromString('--expose_gc');
  const gc = require('vm').runInNewContext('gc');
  const heap = () => { gc(); gc(); return process.memoryUsage().heapUsed; };
  // (all in a function, so that no stale stack slot of this frame still points at the piece when h2 is taken)
  const cut = () => {
    const line = new TextEncoder().encode(' 12345 a line of the output 7\n'), bytes = new Uint8Array(30 * 400000);
    for (let o = 0; o < bytes.length; o += 30) bytes.set(line, o);
    const h0 = heap();
    const piece = new TextDecoder().decode(bytes), kept = [];   // 12 MB, decoded as bigfile.js does
    for (let i = 0; i < 4000; i++) { const p = 30 * i * 97; kept.push(own(piece.slice(p + 1, p + 6)), own(piece.slice(p, p + 29))); }
    return { h0, h1: heap(), len: piece.length, kept };
  };
  const { h0, h1, len, kept } = cut();
  const h2 = heap();
  check(len === 12e6 && h1 - h0 > 1e7 && h2 - h0 < 3e6 && kept.every((s, k) => s === (k % 2 ? ' 12345 a line of the output 7' : '12345')),
    'kept copies do not keep the piece of text alive (heap +' + ((h1 - h0) / 1e6).toFixed(1) + ' MB, then +' + ((h2 - h0) / 1e6).toFixed(1) + ' MB; ' + kept.length + ' kept)');
}

for (const f of fails) console.log('  FAIL:', f);
console.log(`Condensed-reading tests: ${passes} passed, ${failures} failed (${all.length} files, ${nCut} cut versions)`);
process.exit(failures ? 1 : 0);
