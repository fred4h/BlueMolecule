#!/usr/bin/env node
/* BlueMolecule tests for big files: BM.bigfile.readLines (Blob, stream and chunk sources), BM.parseLines, readers on
 * line arrays, and the slice readers for DCD / TRR / XTC (stride, frame range, reload).
 *
 *   node tests/test_bigfiles.js            fast checks (a few seconds; ~40 MB generated logs, deleted afterwards)
 *   node tests/test_bigfiles.js --big      also ~600 MB Gaussian and ORCA logs, a ~1.2 GB log and a ~1.1 GB DCD:
 *                                          Node (each in a child process: time, peak RSS) and headless Chromium
 *                                          through the real UI (dist/BlueMolecule.html, setInputFiles)
 *   options: --keep (keep generated files), --no-browser; BM_BIG_DIR=dir for the generated files (default: a
 *   temporary directory). Big logs come from tools/make_big_fixture.py.
 */
'use strict';
const fs = require('fs');
const os = require('os');
const path = require('path');
const { spawnSync, execSync } = require('child_process');
const ROOT = path.join(__dirname, '..');
for (const f of ['elements.js', 'bigfile.js', 'parsers.js', 'readers-struct.js', 'writers.js', 'readers-qc.js']) require(path.join(ROOT, 'src', f));
const BM = globalThis.BM;
const BF = BM.bigfile;

// ------------------------------------------------------------------ child mode: one big file in a fresh process
if (process.argv[2] === '--child') {
  (async () => {
    const file = process.argv[3], kind = process.argv[4] || 'text';
    const out = { file: path.basename(file), bytes: fs.statSync(file).size };
    const t0 = Date.now();
    try {
      const blob = await fs.openAsBlob(file);
      if (kind === 'dcd') {
        const ds = await BM.parseSource(path.basename(file), blob, {});
        Object.assign(out, { readS: (Date.now() - t0) / 1000, frames: ds.frames.length, totalFrames: ds.meta.totalFrames, stride: ds.meta.stride,
          natoms: ds.meta.natoms, sample: Array.from(ds.frames[ds.frames.length - 1].xyz.slice(0, 3)) });
      } else {
        const lines = await BF.readLines(blob, {});
        out.readS = (Date.now() - t0) / 1000; out.lines = lines.length;
        const t1 = Date.now();
        const ds = BM.parseLines(path.basename(file), lines);
        out.parseS = (Date.now() - t1) / 1000;
        out.frames = ds.frames.length; out.atoms = ds.frames[0].Z.length;
        out.energies = ds.frames.map(f => f.energy);
      }
    } catch (e) { out.error = e.message; out.errorName = e.name; out.oom = BF.isOutOfMemory(e); }
    out.totalS = (Date.now() - t0) / 1000;
    out.maxRssMB = Math.round(process.resourceUsage().maxRSS / 1024);
    process.stdout.write(JSON.stringify(out));
  })();
  return;
}

const BIGRUN = process.argv.includes('--big'), KEEP = process.argv.includes('--keep'), NOBROWSER = process.argv.includes('--no-browser');
let pass = 0, fail = 0, current = '';
const failures = [];
function ok(c, msg) { if (c) pass++; else { fail++; failures.push('[' + current + '] ' + msg); } }
function eq(a, b, msg) { ok(a === b, msg + ': got ' + JSON.stringify(a) + ', expected ' + JSON.stringify(b)); }
function arrNear(a, b, tol, msg) {
  if (!a || a.length !== b.length) { ok(false, msg + ': length ' + (a && a.length) + ' vs ' + b.length); return; }
  for (let i = 0; i < b.length; i++) if (!(Math.abs(a[i] - b[i]) <= tol)) { ok(false, msg + '[' + i + ']: got ' + a[i] + ', expected ' + b[i]); return; }
  pass++;
}
async function test(name, fn) { current = name; try { await fn(); } catch (e) { fail++; failures.push('[' + name + '] threw: ' + (e && e.stack || e)); } }
const ser = ds => JSON.stringify(ds, (k, v) => (ArrayBuffer.isView(v) ? Array.from(v) : v));
const split = s => s.replace(/^﻿/, '').split(/\r\n|\r|\n/);
const sources = u8 => ({
  bytes: () => u8,
  blob: () => new Blob([u8]),
  chunks: () => ({ name: 'x', size: u8.length, readChunk: async (o, n) => u8.slice(o, o + n) }),
  stream: () => { let o = 0; return new ReadableStream({ pull(c) { if (o >= u8.length) { c.close(); return; } const n = 1 + ((Math.random() * 9000) | 0); c.enqueue(u8.slice(o, o + n)); o += n; } }); }
});
const WORK = process.env.BM_BIG_DIR || fs.mkdtempSync(path.join(os.tmpdir(), 'bm-big-'));
fs.mkdirSync(WORK, { recursive: true });
const made = [];
function gen(kind, name, mb) {
  const out = path.join(WORK, name);
  const r = spawnSync('python3', [path.join(ROOT, 'tools', 'make_big_fixture.py'), kind, out, '--mb', String(mb)], { encoding: 'utf8' });
  if (r.status !== 0) throw new Error('make_big_fixture.py failed: ' + r.stderr);
  made.push(out, out.replace(/\.(log|out)$/, '.json'));
  return { file: out, exp: JSON.parse(fs.readFileSync(out.replace(/\.(log|out)$/, '.json'), 'utf8')) };
}

// ------------------------------------------------------------------ trajectory writers (independent of the readers)
function dcdBytes(o) {
  const n = o.natoms, le = true;
  const i32 = v => { const b = Buffer.alloc(4); b.writeInt32LE(v); return b; };
  const rec = body => Buffer.concat([i32(body.length), body, i32(body.length)]);
  const f32 = arr => { const b = Buffer.alloc(4 * arr.length); arr.forEach((v, k) => b.writeFloatLE(v, 4 * k)); return b; };
  const f64 = arr => { const b = Buffer.alloc(8 * arr.length); arr.forEach((v, k) => b.writeDoubleLE(v, 8 * k)); return b; };
  const hdr = Buffer.alloc(84); hdr.write('CORD', 0, 'latin1');
  const ic = new Array(20).fill(0);
  ic[0] = o.frames.length; ic[2] = 1; ic[8] = o.fixed ? o.fixed.length : 0; ic[10] = o.cell ? 1 : 0; ic[19] = 24;
  ic.forEach((v, k) => { if (k !== 9) hdr.writeInt32LE(v, 4 + 4 * k); });
  hdr.writeFloatLE(0.002, 40);
  const t = Buffer.alloc(84, 32); t.writeInt32LE(1, 0); t.write('big test', 4, 'latin1');
  const parts = [rec(hdr), rec(t), rec(i32(n))];
  const free = o.fixed ? [...Array(n).keys()].filter(k => o.fixed.indexOf(k) < 0) : null;
  if (free) parts.push(rec(Buffer.concat(free.map(k => i32(k + 1)))));
  o.frames.forEach((fr, fi) => {
    if (o.cell) parts.push(rec(f64([o.cell[0], 0, o.cell[1], 0, 0, o.cell[2]])));
    const idx = free && fi > 0 ? free : [...Array(n).keys()];
    for (let d = 0; d < 3; d++) parts.push(rec(f32(idx.map(k => fr[3 * k + d]))));
  });
  return new Uint8Array(Buffer.concat(parts));
}
function xdr() {
  const parts = [];
  const api = {
    int: v => { const b = Buffer.alloc(4); b.writeInt32BE(v); parts.push(b); return api; },
    f32: v => { const b = Buffer.alloc(4); b.writeFloatBE(v); parts.push(b); return api; },
    raw: b => { parts.push(Buffer.from(b)); return api; },
    bytes: () => new Uint8Array(Buffer.concat(parts))
  };
  return api;
}
function trrBytes(frames, n) {
  const w = xdr();
  frames.forEach((fr, k) => {
    w.int(1993).int(13).int(12).raw('GMX_trn_file');
    [0, 0, 36, 0, 0, 0, 0, 12 * n, 0, 0, n, 10 * k, 0].forEach(v => w.int(v));
    w.f32(0.5 * k); w.f32(0);
    fr.box.forEach(w.f32); fr.x.forEach(w.f32);
  });
  return w.bytes();
}
function xtcSmall(frames, n) {           // <= 9 atoms: uncompressed XTC frames
  const w = xdr();
  frames.forEach((fr, k) => { w.int(1995).int(n).int(5 * k).f32(0.1 * k); fr.box.forEach(w.f32); w.int(n); fr.x.forEach(w.f32); });
  return w.bytes();
}

(async () => {
  // ================================================================ readLines
  await test('readLines equals String.split for every source and piece size', async () => {
    const al = ['a', 'b', '\r', '\n', '\r\n', ' ', 'é', '日', '😀'];
    let bad = 0;
    for (let t = 0; t < 400; t++) {
      let s = (t % 7 === 0 ? '﻿' : '');
      const len = t % 5 === 0 ? 20000 : (Math.random() * 80) | 0;
      for (let i = 0; i < len; i++) s += t % 5 === 0 && i % 997 ? 'x\n'[i % 13 ? 0 : 1] : al[(Math.random() * al.length) | 0];
      const u8 = new TextEncoder().encode(s), exp = JSON.stringify(split(s));
      const cs = t % 5 === 0 ? 4096 + ((Math.random() * 3000) | 0) : 1 + ((Math.random() * 9) | 0);
      for (const [kind, mk] of Object.entries(sources(u8))) {
        const got = await BF.readLines(mk(), { chunkSize: cs });
        if (JSON.stringify(got) !== exp && bad++ < 3) failures.push('[readLines] ' + kind + ' chunk ' + cs + ' differs for ' + JSON.stringify(s.slice(0, 60)));
      }
    }
    eq(bad, 0, 'mismatching splits');
    eq(JSON.stringify(await BF.readLines(new Uint8Array(0))), '[""]', 'empty file -> [""]');
    eq(JSON.stringify(await BF.readLines(new TextEncoder().encode('a\r'), { chunkSize: 1 })), '["a",""]', 'CR at the very end');
    eq(JSON.stringify(await BF.readLines(new TextEncoder().encode('a\r\nb'), { chunkSize: 2 })), '["a","b"]', 'CR LF across pieces');
    const bad8 = new Uint8Array([0x61, 0xE6, 0x97, 0x0A, 0x62, 0xC3]);
    eq(JSON.stringify(await BF.readLines(bad8, { chunkSize: 2 })), JSON.stringify(new TextDecoder().decode(bad8).split('\n')), 'invalid UTF-8 -> U+FFFD as TextDecoder');
  });

  await test('readLines progress, cancellation and read errors', async () => {
    const u8 = new TextEncoder().encode('line\n'.repeat(400000));     // 2 MB
    const seen = [];
    await BF.readLines(new Blob([u8]), { chunkSize: 65536, onProgress: p => seen.push(p) });
    ok(seen.length >= 1 && seen[seen.length - 1].done && seen[seen.length - 1].bytes === u8.length, 'final progress report');
    ok(seen.every((p, k) => !k || p.bytes >= seen[k - 1].bytes), 'bytes grow');
    const ctl = new AbortController();
    let n = 0;
    const src = { name: 'c', size: u8.length, readChunk: async (o, len) => { if (++n === 3) ctl.abort(); return u8.slice(o, o + len); } };
    let err = null;
    try { await BF.readLines(src, { chunkSize: 65536, signal: ctl.signal }); } catch (e) { err = e; }
    ok(err && BF.isAbort(err), 'aborted read rejects with AbortError');
    ok(n <= 4, 'reading stops right after the abort (' + n + ' pieces)');
    let cancelled = false;
    const st = new ReadableStream({ pull(c) { c.enqueue(u8.slice(0, 1000)); }, cancel() { cancelled = true; } });
    const c2 = new AbortController();
    setTimeout(() => c2.abort(), 20);
    err = null;
    try { await BF.readLines(st, { signal: c2.signal }); } catch (e) { err = e; }
    ok(err && BF.isAbort(err) && cancelled, 'endless stream: aborted and cancelled');
    err = null;
    try { await BF.readLines({ name: 'e', size: 100, readChunk: async () => { throw new Error('disk gone'); } }); } catch (e) { err = e; }
    ok(err && /disk gone/.test(err.message), 'read errors propagate');
    err = null;
    try { await BF.readLines({ name: 's', size: 100, readChunk: async () => new Uint8Array(0) }); } catch (e) { err = e; }
    ok(err && /past byte 0/.test(err.message), 'a short source is reported');
  });

  await test('memory guard and error messages', async () => {
    let err = null;
    try { await BF.readLines(new Blob([new Uint8Array(40 * 1048576).fill(97)]), { memoryLimit: 32 * 1048576 }); } catch (e) { err = e; }
    ok(err && err.name === 'OutOfMemoryError' && BF.isOutOfMemory(err), 'projected heap above the limit: OutOfMemoryError');
    const msg = BF.describeError(err, 'huge.log', 3e9);
    ok(/Not enough memory to open "huge\.log" \(2\.79 GB\)/.test(msg), 'clear message: ' + msg);
    ok(BF.isOutOfMemory(new RangeError('Invalid string length')) && BF.isOutOfMemory(new RangeError('Array buffer allocation failed')), 'engine limits');
    ok(!BF.isOutOfMemory(new Error('No geometry found')), 'other errors are not memory errors');
    ok(/cancelled/.test(BF.describeError(Object.assign(new Error('x'), { name: 'AbortError' }), 'a.log')), 'cancel message');
  });

  await test('toSource normalises host objects', async () => {
    const u8 = new TextEncoder().encode('3\nw\nO 0 0 0\nH 0 0 1\nH 0 1 0\n');
    const kinds = [[new Blob([u8]), 'blob'], [{ name: 'a.xyz', bytes: u8 }, 'bytes'], [{ name: 'a.xyz', base64: Buffer.from(u8).toString('base64') }, 'bytes'],
      [{ name: 'a.xyz', size: u8.length, readChunk: async (o, n) => u8.slice(o, o + n) }, 'chunks'], [{ name: 'a.xyz', text: 'x' }, 'text'], [u8.buffer, 'bytes']];
    for (const [x, k] of kinds) eq(BF.toSource(x, 'a.xyz').kind, k, 'kind ' + k);
    for (const [x] of kinds.slice(0, 4)) eq((await BF.readLines(BF.toSource(x, 'a.xyz'))).length, 6, 'lines through toSource');
    let err = null; try { BF.toSource(42, 'n'); } catch (e) { err = e; }
    ok(err && /Unsupported/.test(err.message), 'unsupported object');
    eq(await BF.readHead(new Blob([u8]), 5), '3\nw\nO', 'readHead');
  });

  // ================================================================ parseLines and readers on line arrays
  await test('parseLines gives the same dataset as parseFile for every fixture', async () => {
    const walk = d => fs.readdirSync(d, { withFileTypes: true }).flatMap(e => (e.isDirectory() ? walk(path.join(d, e.name)) : [path.join(d, e.name)]));
    let n = 0;
    for (const f of walk(path.join(__dirname, 'fixtures'))) {
      const name = path.basename(f);
      if (BM.isBinaryFormat(name) || /\.(png|json|gz)$/.test(name)) continue;
      const text = fs.readFileSync(f, 'utf8');
      let a, b;
      try { a = ser(BM.parseFile(name, text)); } catch (e) { a = 'ERR ' + e.message; }
      const lines = await BF.readLines(new Blob([fs.readFileSync(f)]), { chunkSize: 4099 });
      try { b = ser(BM.parseLines(name, lines)); } catch (e) { b = 'ERR ' + e.message; }
      ok(a === b, name + ': parseLines = parseFile');
      n++;
    }
    ok(n > 100, n + ' fixtures compared');
  });

  await test('parseLines edge cases', async () => {
    const err = f => { try { f(); return null; } catch (e) { return e.message; } };
    ok(/is empty/.test(err(() => BM.parseLines('e.xyz', ['', '   ', '']))), 'empty');
    ok(/binary file/.test(err(() => BM.parseLines('b.xyz', ['\u0001\u0002', 'x']))), 'binary');
    ok(/not recognised/.test(err(() => BM.parseLines('q.zzz', ['hello', 'world']))), 'unknown format');
    eq(BM.parseLines('w.xyz', '3\nw\nO 0 0 0\nH 0 0 1\nH 0 1 0\n').frames.length, 1, 'a string goes to parseFile');
    eq(BM.parseLines('w.xyz', ['﻿3', 'w', 'O 0 0 0', 'H 0 0 1', 'H 0 1 0']).frames[0].Z.length, 3, 'BOM on the first line');
    const cml = fs.readFileSync(path.join(__dirname, 'fixtures', 'struct', 'array.cml'), 'utf8').split('\n');
    eq(BM.parseLines('a.cml', cml).frames[0].Z.length > 0, true, 'CML (whole-text reader) from lines');
    const keep = BM.maxJoinChars;
    BM.maxJoinChars = 100;
    ok(/too large for the CML reader/.test(err(() => BM.parseLines('a.cml', cml))), 'too large for a whole-text reader: clear error');
    eq(BM.parseLines('w.xyz', ['3', 'w', 'O 0 0 0', 'H 0 0 1', 'H 0 1 0']).frames.length, 1, 'line-based readers have no join limit');
    BM.maxJoinChars = keep;
    const cube = fs.readFileSync(path.join(__dirname, 'fixtures', 'water_mo5_mo6.cube'), 'utf8');
    const a = BM.parseFile('w.cube', cube), b = BM.parseLines('w.cube', split(cube));
    ok(a.volumes.length === 2 && a.volumes.every((v, k) => v.data.every((x, i) => x === b.volumes[k].data[i])), 'cube values from lines');
  });

  await test('generated 40 MB Gaussian / ORCA logs through the big-file path', async () => {
    for (const [kind, name] of [['gaussian', 'midg.log'], ['orca', 'mido.out']]) {
      const { file, exp } = gen(kind, name, 40);
      const t0 = Date.now();
      const lines = await BF.readLines(await fs.openAsBlob(file), { chunkSize: 8 * 1048576 });
      const ds = BM.parseLines(name, lines);
      const dt = Date.now() - t0;
      eq(ds.frames.length, exp.frames, name + ' frames');
      arrNear(ds.frames.map(f => f.energy), exp.energies, 1e-7, name + ' energies');
      ok(ds.frames.every(f => f.Z.length === exp.atoms), name + ' atoms');
      ok(ds.frames[ds.frames.length - 1].converged, name + ' converged');
      eq(ser(ds), ser(BM.parseFile(name, fs.readFileSync(file, 'utf8'))), name + ': same dataset as from one string');
      ok(dt < 15000, name + ' read + parse in ' + dt + ' ms');
      console.log('  ' + name + ': ' + (exp.bytes / 1e6).toFixed(0) + ' MB, ' + lines.length + ' lines, read + parse ' + dt + ' ms');
    }
  });

  // ================================================================ trajectories from slices
  await test('DCD from slices: stride, range, fixed atoms, cell, reload, budget', async () => {
    const n = 5, F = 23;
    const frames = Array.from({ length: F }, (_, k) => Array.from({ length: 3 * n }, (_, i) => +(0.1 * k + 0.01 * i).toFixed(3)));
    for (const o of [{}, { cell: [20, 21, 22] }, { fixed: [1, 3] }]) {
      const u8 = dcdBytes(Object.assign({ natoms: n, frames }, o));
      const ref = BM.parseFileBinary('t.dcd', u8);
      const tag = JSON.stringify(o);
      for (const src of [new Blob([u8]), { name: 't.dcd', size: u8.length, readChunk: async (p, l) => u8.slice(p, p + l) }]) {
        const all = await BM.parseSource('t.dcd', src);
        eq(all.frames.length, F, tag + ' all frames'); eq(all.meta.totalFrames, F, tag + ' totalFrames'); eq(all.meta.stride, 1, tag + ' stride 1');
        ok(all.frames.every((f, k) => f.xyz.every((v, i) => v === ref.frames[k].xyz[i])), tag + ' same coordinates as the whole-file reader');
        if (o.cell) ok(JSON.stringify(all.cell) === JSON.stringify(ref.cell), tag + ' cell');
        const st = await BM.parseSource('t.dcd', src, { stride: 5, first: 2, last: 20 });
        eq(st.frames.map(f => f.label).join(','), 'Frame 3,Frame 8,Frame 13,Frame 18', tag + ' stride 5 from frame 3 to 21');
        ok(st.frames.every((f, k) => f.xyz.every((v, i) => v === ref.frames[2 + 5 * k].xyz[i])), tag + ' strided coordinates');
        eq(st.meta.stride, 5, tag + ' meta.stride'); eq(st.meta.first, 2, tag + ' meta.first'); eq(st.meta.last, 20, tag + ' meta.last');
        ok(st.warnings.some(w => /Loaded 4 of 23 frames/.test(w)), tag + ' warning says what was loaded');
        const re = await st.meta.reload({ stride: 11 });
        eq(re.frames.length, 3, tag + ' reload every 11th'); eq(re.meta.first, 0, tag + ' reload resets the range');
        const auto = await BM.parseSource('t.dcd', src, { budget: 24 * n * 6 });
        eq(auto.meta.stride, 4, tag + ' budget for 6 frames: stride 4'); eq(auto.frames.length, 6, tag + ' 6 frames kept');
        ok(!Object.keys(st.meta).includes('reload'), tag + ' reload is not enumerable');
      }
    }
    const cut = dcdBytes({ natoms: n, frames }).slice(0, -10);
    const ds = await BM.parseSource('c.dcd', new Blob([cut]));
    eq(ds.frames.length, F - 1, 'truncated file: complete frames only');
    ok(ds.warnings.some(w => /incomplete/.test(w)), 'truncation warning');
    let err = null; try { await BM.parseSource('x.dcd', new Blob([new Uint8Array(100)])); } catch (e) { err = e; }
    ok(err && /Not a DCD file/.test(err.message), 'not a DCD');
    const ctl = new AbortController(); ctl.abort(); err = null;
    try { await BM.parseSource('t.dcd', new Blob([dcdBytes({ natoms: n, frames })]), { signal: ctl.signal }); } catch (e) { err = e; }
    ok(err && err.name === 'AbortError', 'cancel');
  });

  await test('TRR / XTC from slices: stride and range', async () => {
    const n = 3, F = 10, box = [3, 0, 0, 0, 3, 0, 0, 0, 3];
    const frames = Array.from({ length: F }, (_, k) => ({ box, x: Array.from({ length: 3 * n }, (_, i) => +(0.05 * k + 0.1 * i).toFixed(3)) }));
    for (const [nm, u8] of [['t.trr', trrBytes(frames, n)], ['t.xtc', xtcSmall(frames, n)]]) {
      const ref = BM.parseFileBinary(nm, u8);
      const all = await BM.parseSource(nm, new Blob([u8]));
      eq(all.frames.length, F, nm + ' frames'); eq(all.meta.totalFrames, F, nm + ' totalFrames');
      ok(all.frames.every((f, k) => f.xyz.every((v, i) => v === ref.frames[k].xyz[i]) && f.time === ref.frames[k].time), nm + ' = whole-file reader');
      const st = await BM.parseSource(nm, new Blob([u8]), { stride: 3, first: 1 });
      eq(st.frames.map(f => +f.time.toFixed(3)).join(','), nm === 't.trr' ? '0.5,2,3.5' : '0.1,0.4,0.7', nm + ' stride 3 from the 2nd frame');
      eq(st.meta.totalFrames, F, nm + ' all frames counted');
      eq((await st.meta.reload({ stride: 5 })).frames.length, 2, nm + ' reload');
    }
  });

  // ================================================================ the app's reading paths in Chromium
  if (!NOBROWSER && fs.existsSync(path.join(ROOT, 'dist', 'BlueMolecule.html'))) await hostTests();

  // ================================================================ big files (opt-in)
  if (BIGRUN) await bigTests();

  if (!KEEP) for (const f of made) try { fs.unlinkSync(f); } catch (e) { /* gone */ }
  if (!KEEP && !process.env.BM_BIG_DIR) try { fs.rmdirSync(WORK); } catch (e) { /* not empty */ }
  console.log('');
  for (const m of failures) console.log('FAIL ' + m);
  console.log('BlueMolecule big-file tests: ' + pass + ' passed, ' + fail + ' failed' + (BIGRUN ? '' : ' (run with --big for the 600 MB / 1.2 GB files)'));
  process.exit(fail ? 1 : 0);
})();

// ------------------------------------------------------------------ host paths (default run)
/** A ~70 MB log (above BM.bigfile.BIG) opened as a chunk object pushed by a native shell, through a native
 *  open dialog, streamed by the launcher (fetch body), and cancelled while reading. */
async function hostTests() {
  const { chromium } = require('playwright');
  const { file, exp } = gen('gaussian', 'host70.log', 70);
  const browser = await chromium.launch({ args: ['--enable-unsafe-swiftshader'] });
  const page = await browser.newPage();
  const errors = [];
  page.on('pageerror', e => errors.push(e.message));
  // capture the callback the app registers with BM.native.onOpenFile (what a shell's file delivery calls)
  await page.addInitScript(() => {
    let nat = null;
    const B = window.BM = window.BM || {};
    Object.defineProperty(B, 'native', { configurable: true, get() { return nat; }, set(v) {
      nat = v;
      const orig = v.onOpenFile;
      v.onOpenFile = cb => { window.__hostCb = cb; return orig.call(v, cb); };
    } });
  });
  await page.goto('file://' + path.join(ROOT, 'dist', 'BlueMolecule.html'));
  await page.waitForFunction(() => window.BM && window.BM.app && window.__hostCb);
  await page.evaluate(() => { const i = document.createElement('input'); i.type = 'file'; i.id = 'tmp-in'; document.body.appendChild(i); });
  await page.setInputFiles('#tmp-in', file);
  // in-page helpers: a chunk object over the File (optionally slow), and a document-count waiter
  await page.evaluate(() => {
    const f = document.querySelector('#tmp-in').files[0];
    window.__chunks = (delay, name) => ({ name: name || 'host70.log', size: f.size, reads: 0,
      readChunk(o, n) { this.reads++; return new Promise(r => setTimeout(r, delay || 0)).then(() => f.slice(o, o + n).arrayBuffer()).then(b => new Uint8Array(b)); } });
    window.__wait = (n, ms) => new Promise((res, rej) => { const t0 = Date.now(); const iv = setInterval(() => {
      const s = window.BM.app.state, bar = document.querySelector('#bm-progress');
      if (s.docs.length >= n && (!bar || bar.style.display === 'none')) { clearInterval(iv); res(Date.now() - t0); }
      else if (Date.now() - t0 > ms) { clearInterval(iv); rej(new Error('timeout; toast: ' + document.querySelector('#toast').textContent)); }
    }, 50); });
    window.__doc = () => { const d = window.BM.app.state.active; return { name: d.name, frames: d.data.frames.length, E: d.data.frames.map(x => x.energy) }; };
  });
  await test('host: a big file pushed as a chunk object (native shell)', async () => {
    await page.evaluate(() => { window.__c = window.__chunks(0, 'pushed.log'); window.__hostCb('pushed.log', window.__c); });
    await page.evaluate(() => window.__wait(1, 60000));
    const d = await page.evaluate(() => Object.assign(window.__doc(), { reads: window.__c.reads }));
    eq(d.name, 'pushed.log', 'document'); eq(d.frames, exp.frames, 'frames');
    arrNear(d.E, exp.energies, 1e-7, 'energies');
    ok(d.reads >= 3, 'read in slices (' + d.reads + ' readChunk calls)');
    const small = await page.evaluate(() => { window.__hostCb('w.xyz', new TextEncoder().encode('3\nw\nO 0 0 0\nH 0 0 0.96\nH 0.93 0 -0.24\n')); return window.__wait(2, 10000).then(() => window.__doc()); });
    eq(small.frames, 1, 'small file pushed as (name, Uint8Array)');
  });
  await test('host: Open button with BM.native.openDialog', async () => {
    await page.evaluate(() => { window.BM.native.openDialog = async () => [window.__chunks(0, 'dialog.log')]; });
    await page.click('#btn-open');
    await page.evaluate(() => window.__wait(3, 60000));
    const d = await page.evaluate(() => window.__doc());
    eq(d.name, 'dialog.log', 'opened from the dialog'); eq(d.frames, exp.frames, 'frames');
    await page.evaluate(() => { window.BM.native.openDialog = async () => undefined; });
    await page.click('#btn-open');
    await page.waitForTimeout(200);
    eq(await page.evaluate(() => window.BM.app.state.docs.length), 3, 'cancelled dialog opens nothing');
    await page.evaluate(() => { delete window.BM.native.openDialog; });
  });
  await test('host: launcher files are streamed (fetch body)', async () => {
    await page.evaluate(() => { window.BM_LAUNCH = { files: [{ name: 'launched.log', url: URL.createObjectURL(document.querySelector('#tmp-in').files[0]) }] }; return window.BM.app.loadLaunchFiles(); });
    await page.evaluate(() => window.__wait(4, 60000));
    const d = await page.evaluate(() => window.__doc());
    eq(d.name, 'launched.log', 'launcher document'); eq(d.frames, exp.frames, 'frames');
    arrNear(d.E, exp.energies, 1e-7, 'energies');
  });
  await test('host: progress bar and Cancel', async () => {
    await page.evaluate(() => { window.__c2 = window.__chunks(400, 'slow.log'); window.__hostCb('slow.log', window.__c2); });
    await page.waitForFunction(() => { const b = document.querySelector('#bm-progress'); return b && b.style.display !== 'none'; }, null, { timeout: 10000 });
    const txt = await page.textContent('#bm-progress');
    ok(/Reading slow\.log/.test(txt) && /Cancel/.test(txt), 'progress card: ' + txt.trim());
    await page.click('#bm-progress .bmp-x');
    await page.waitForFunction(() => /cancelled/.test(document.querySelector('#toast').textContent), null, { timeout: 10000 });
    const st = await page.evaluate(() => ({ docs: window.BM.app.state.docs.length, reads: window.__c2.reads, bar: document.querySelector('#bm-progress').style.display }));
    eq(st.docs, 4, 'no document after Cancel'); ok(st.reads < 3, 'reading stopped (' + st.reads + ' slices)'); eq(st.bar, 'none', 'bar hidden');
  });
  await test('host: .gz file decompressed as a stream', async () => {
    const gz = file + '.gz';
    fs.writeFileSync(gz, require('zlib').gzipSync(fs.readFileSync(file), { level: 1 }));
    made.push(gz);
    await page.setInputFiles('#file-input', gz);
    await page.evaluate(() => window.__wait(5, 60000));
    const d = await page.evaluate(() => window.__doc());
    eq(d.name, 'host70.log', 'name without .gz'); eq(d.frames, exp.frames, 'frames from the decompressed stream');
  });
  await test('host: a file that cannot fit the heap is refused with a clear message', async () => {
    await page.evaluate(() => { Object.defineProperty(performance, 'memory', { configurable: true, get: () => ({ jsHeapSizeLimit: 60e6, usedJSHeapSize: 20e6, totalJSHeapSize: 30e6 }) }); });
    await page.evaluate(() => { window.__hostCb('huge.log', window.__chunks(0, 'huge.log')); });
    await page.waitForFunction(() => /Not enough memory/.test(document.querySelector('#toast').textContent), null, { timeout: 20000 });
    const t = await page.textContent('#toast');
    ok(/Not enough memory to open "huge\.log" \(\d+ MB\)/.test(t), 'toast: ' + t);
    eq(await page.evaluate(() => window.BM.app.state.docs.length), 5, 'nothing opened');
  });
  ok(!errors.length, 'no page errors: ' + errors.join('; '));
  await browser.close();
}

// ------------------------------------------------------------------ --big
function child(file, kind) {
  const r = spawnSync(process.execPath, [__filename, '--child', file, kind || 'text'], { encoding: 'utf8', maxBuffer: 1 << 26, timeout: 600000 });
  try { return JSON.parse(r.stdout); } catch (e) { return { error: 'child failed (status ' + r.status + ', signal ' + r.signal + '): ' + (r.stderr || '').slice(-400) }; }
}
/** Atom i of the big test system on a 3.5 Å grid (30 x 30 x n): no bonds, cheap to draw. */
const gridPos = i => [3.5 * (i % 30), 3.5 * (Math.floor(i / 30) % 30), 3.5 * Math.floor(i / 900)];
function writeBigDcd(file, natoms, nframes) {
  // frame k: atom i at its grid position, y shifted by k * 0.001 Å; written in 64 MB pieces
  const fd = fs.openSync(file, 'w');
  const hb = Buffer.from(dcdBytes({ natoms, frames: [] })); hb.writeInt32LE(nframes, 8);
  fs.writeSync(fd, hb);
  const rec = 4 * natoms, fb = Buffer.alloc(3 * (rec + 8));
  for (let d = 0; d < 3; d++) {
    const o = d * (rec + 8);
    fb.writeInt32LE(rec, o); fb.writeInt32LE(rec, o + 4 + rec);
    if (d !== 1) for (let i = 0; i < natoms; i++) fb.writeFloatLE(gridPos(i)[d], o + 4 + 4 * i);
  }
  const buf = [];
  let size = 0;
  for (let k = 0; k < nframes; k++) {
    const o = rec + 8;
    for (let i = 0; i < natoms; i++) fb.writeFloatLE(gridPos(i)[1] + k * 0.001, o + 4 + 4 * i);
    buf.push(Buffer.from(fb)); size += fb.length;
    if (size > 64 * 1048576) { fs.writeSync(fd, Buffer.concat(buf)); buf.length = 0; size = 0; }
  }
  if (buf.length) fs.writeSync(fd, Buffer.concat(buf));
  fs.closeSync(fd);
  made.push(file);
}
const rendererRss = () => { try { return +execSync("ps -eo rss,args | grep -- '--type=renderer' | grep -v grep | awk '{s+=$1} END {print s}'").toString().trim() / 1024 | 0; } catch (e) { return -1; } };

async function bigTests() {
  const rows = [];
  const files = {};
  current = 'big: generate';
  for (const [kind, name, mb] of [['gaussian', 'big600g.log', 600], ['orca', 'big600o.out', 600], ['gaussian', 'big1200g.log', 1200]]) {
    const t0 = Date.now();
    files[name] = gen(kind, name, mb);
    console.log('  generated ' + name + ': ' + (files[name].exp.bytes / 1e6).toFixed(0) + ' MB in ' + ((Date.now() - t0) / 1000).toFixed(1) + ' s');
  }
  const dcd = path.join(WORK, 'big.dcd');
  writeBigDcd(dcd, 20000, 4600);
  console.log('  generated big.dcd: ' + (fs.statSync(dcd).size / 1e6).toFixed(0) + ' MB (20000 atoms, 4600 frames)');

  await test('big: Node, ~600 MB Gaussian and ORCA logs', async () => {
    for (const name of ['big600g.log', 'big600o.out']) {
      const { file, exp } = files[name], r = child(file);
      ok(!r.error, name + ': ' + (r.error || 'ok'));
      eq(r.frames, exp.frames, name + ' frames'); eq(r.atoms, exp.atoms, name + ' atoms');
      if (r.energies) arrNear(r.energies, exp.energies, 1e-7, name + ' energies');
      ok(r.totalS < 60, name + ' well under a minute: ' + r.totalS + ' s');
      rows.push(['Node', name, exp.bytes, r.lines, r.readS, r.parseS, r.totalS, r.maxRssMB, r.frames + ' frames']);
    }
  });
  await test('big: Node, ~1.2 GB log (works, or a clear out-of-memory error)', async () => {
    const { file, exp } = files['big1200g.log'], r = child(file);
    if (r.error) ok(r.oom && /memory/i.test(r.error), '1.2 GB: failure must be a clear memory error: ' + r.error);
    else { eq(r.frames, exp.frames, '1.2 GB frames'); arrNear(r.energies, exp.energies, 1e-7, '1.2 GB energies'); }
    rows.push(['Node', 'big1200g.log', exp.bytes, r.lines, r.readS, r.parseS, r.totalS, r.maxRssMB, r.error ? 'ERROR ' + r.error : r.frames + ' frames']);
  });
  await test('big: Node, 1.1 GB DCD as every Nth frame', async () => {
    const r = child(dcd, 'dcd');
    ok(!r.error, 'dcd: ' + (r.error || 'ok'));
    eq(r.totalFrames, 4600, 'totalFrames'); eq(r.natoms, 20000, 'natoms');
    ok(r.stride > 1 && r.frames === Math.ceil(4600 / r.stride), 'auto stride ' + r.stride + ' -> ' + r.frames + ' frames');
    ok(r.sample && Math.abs(r.sample[1] - 0.001 * (r.stride * (r.frames - 1))) < 1e-5, 'coordinates of the last frame read');
    rows.push(['Node', 'big.dcd', fs.statSync(dcd).size, null, r.readS, null, r.totalS, r.maxRssMB, r.frames + '/' + r.totalFrames + ' frames, stride ' + r.stride]);
  });

  if (!NOBROWSER && fs.existsSync(path.join(ROOT, 'dist', 'BlueMolecule.html'))) {
    const { chromium } = require('playwright');
    const browser = await chromium.launch({ args: ['--enable-unsafe-swiftshader'] });
    const page = await browser.newPage();
    const errors = [];
    page.on('pageerror', e => errors.push(e.message));
    page.on('crash', () => errors.push('renderer crashed'));
    await page.goto('file://' + path.join(ROOT, 'dist', 'BlueMolecule.html'));
    await page.waitForFunction(() => window.BM && window.BM.app && window.BM.bigfile);
    const openBig = async file => {
      const before = await page.evaluate(() => ({ n: window.BM.app.state.docs.length, f: window.BM.app.state.active ? window.BM.app.state.active.data.frames.length : 0 }));
      let peak = 0, bar = false;
      const t0 = Date.now();
      await page.setInputFiles('#file-input', file);
      let st = null;
      for (let i = 0; i < 1500; i++) {
        await page.waitForTimeout(200);
        peak = Math.max(peak, rendererRss());
        st = await page.evaluate(b => {
          const s = window.BM.app.state, el = document.querySelector('#bm-progress');
          return { busy: !!el && el.style.display !== 'none', changed: s.docs.length !== b.n || (s.active && s.active.data.frames.length !== b.f), toast: document.querySelector('#toast').textContent };
        }, before);
        if (st.busy) bar = true;
        if (!st.busy && (st.changed || /could not|not enough|memory|cancel/i.test(st.toast))) break;
      }
      const d = await page.evaluate(() => { const a = window.BM.app.state.active; return a ? { frames: a.data.frames.length, atoms: a.data.frames[0].Z.length, energies: a.data.frames.map(f => f.energy), meta: { totalFrames: a.data.meta.totalFrames, stride: a.data.meta.stride } } : null; });
      return { s: (Date.now() - t0) / 1000, peak, bar, toast: st && st.toast, d };
    };
    await test('big: Chromium opens the ~600 MB logs through the file input', async () => {
      for (const name of ['big600g.log', 'big600o.out']) {
        const { file, exp } = files[name], r = await openBig(file);
        ok(r.bar, name + ': progress bar shown');
        ok(r.d && r.d.frames === exp.frames, name + ' frames ' + (r.d && r.d.frames) + ' / ' + exp.frames + ' (' + r.toast + ')');
        ok(r.d && r.d.energies.filter(Number.isFinite).length === exp.frames, name + ' every frame has its energy');
        if (r.d) arrNear(r.d.energies, exp.energies, 1e-7, name + ' energies');
        rows.push(['Chromium', name, exp.bytes, null, null, null, r.s, r.peak, r.d ? r.d.frames + ' frames' : r.toast]);
        await page.evaluate(() => { const s = window.BM.app.state; s.docs.length = 0; s.active = null; });
      }
    });
    await test('big: Chromium, ~1.2 GB log (works, or a clear out-of-memory message)', async () => {
      const { file, exp } = files['big1200g.log'], r = await openBig(file);
      if (errors.includes('renderer crashed')) ok(false, 'the renderer crashed');
      else if (r.d && r.d.frames === exp.frames) ok(true, '1.2 GB opened');
      else ok(/not enough memory/i.test(r.toast || ''), '1.2 GB: clear memory message, got: ' + r.toast);
      rows.push(['Chromium', 'big1200g.log', exp.bytes, null, null, null, r.s, r.peak, r.d && r.d.frames === exp.frames ? r.d.frames + ' frames' : r.toast]);
      await page.evaluate(() => { const s = window.BM.app.state; s.docs.length = 0; s.active = null; });
    });
    await test('big: Chromium, 20000-atom topology + 1.1 GB DCD', async () => {
      const xyz = path.join(WORK, 'topo20k.xyz');
      fs.writeFileSync(xyz, '20000\nbig topology\n' + Array.from({ length: 20000 }, (_, i) => 'Ar ' + gridPos(i).map(v => v.toFixed(3)).join(' ')).join('\n') + '\n');
      made.push(xyz);
      await openBig(xyz);
      const r = await openBig(dcd);
      ok(r.d && r.d.atoms === 20000 && r.d.frames > 100 && r.d.frames < 4601, 'DCD frames attached: ' + (r.d && r.d.frames) + ' (' + r.toast + ')');
      ok(/of 4600 frames loaded/.test(r.toast || ''), 'toast tells the stride: ' + r.toast);
      rows.push(['Chromium', 'big.dcd', fs.statSync(dcd).size, null, null, null, r.s, r.peak, r.toast]);
    });
    ok(!errors.length, 'no page errors: ' + errors.join('; '));
    await browser.close();
  }
  console.log('\n  where     file          MB      lines        read s  parse s  total s  peak RSS MB  result');
  for (const r of rows) {
    console.log('  ' + r[0].padEnd(9) + ' ' + r[1].padEnd(13) + String(Math.round(r[2] / 1e6)).padStart(5) + (r[3] ? String(r[3]).padStart(12) : ''.padStart(12)) +
      (r[4] !== null && r[4] !== undefined ? r[4].toFixed(1).padStart(10) : ''.padStart(10)) + (r[5] !== null && r[5] !== undefined ? r[5].toFixed(1).padStart(9) : ''.padStart(9)) +
      (r[6] !== undefined ? r[6].toFixed(1).padStart(9) : '') + String(r[7]).padStart(13) + '  ' + r[8]);
  }
}
