#!/usr/bin/env node
/* BlueMolecule Amber NetCDF tests (issue #8).  Run:  node tests/test_netcdf.js
 * The files in tests/fixtures/netcdf were written by SciPy's NetCDF 3 writer (classic and 64-bit offset) and, for the
 * CDF-5 variant that SciPy does not write, by a separate byte writer in Python (tests/make_netcdf.py, which also
 * writes expected.json). Damaged, streaming and NetCDF 4 files are made here by editing those bytes. */
'use strict';
const fs = require('fs');
const path = require('path');
const ROOT = path.join(__dirname, '..');
for (const f of ['elements', 'parsers', 'writers', 'readers-struct', 'bigfile']) require(path.join(ROOT, 'src', f + '.js'));
const BM = globalThis.BM;
const DIR = path.join(__dirname, 'fixtures', 'netcdf');
const EXP = JSON.parse(fs.readFileSync(path.join(DIR, 'expected.json'), 'utf8'));
const bytes = f => new Uint8Array(fs.readFileSync(path.join(DIR, f)));

let pass = 0, fail = 0, current = '';
const failures = [];
function ok(c, msg) { if (c) pass++; else { fail++; failures.push('[' + current + '] ' + msg); } }
function eq(a, b, msg) { ok(a === b, msg + ': got ' + JSON.stringify(a) + ', expected ' + JSON.stringify(b)); }
function same(xyz, want, msg) {   // float32 values stored by the writer, read back exactly
  if (!xyz || xyz.length !== want.length) { ok(false, msg + ': length ' + (xyz && xyz.length) + ' vs ' + want.length); return; }
  for (let i = 0; i < want.length; i++) if (xyz[i] !== want[i]) { ok(false, msg + '[' + i + ']: got ' + xyz[i] + ', expected ' + want[i]); return; }
  pass++;
}
async function test(name, fn) { current = name; try { await fn(); } catch (e) { fail++; failures.push('[' + name + '] threw: ' + (e && e.stack || e)); } }
const cellOf = c => { const n = v => Math.hypot(...v), d = (u, v) => Math.acos((u[0] * v[0] + u[1] * v[1] + u[2] * v[2]) / (n(u) * n(v))) * 180 / Math.PI; return [n(c[0]), n(c[1]), n(c[2]), d(c[1], c[2]), d(c[0], c[2]), d(c[0], c[1])]; };
const closeAll = (a, b, tol) => a.length === b.length && a.every((v, i) => Math.abs(v - b[i]) <= tol);
/** A Blob-like source for parseSource (what the app passes for a big file). */
const blob = u8 => new Blob([u8]);
const chunks = u8 => ({ name: 'x', size: u8.length, readChunk: async (p, l) => u8.slice(p, p + l) });
const withNumrecs = (u8, n) => { const c = u8.slice(); new DataView(c.buffer).setUint32(4, n, false); return c; };   // CDF-1/2: numrecs after the magic

(async () => {
  await test('every fixture, read whole', () => {
    for (const f of Object.keys(EXP)) {
      const e = EXP[f], ds = BM.parseFileBinary(f, bytes(f));
      eq(ds.format, 'netcdf', f + ': format');
      eq(ds.frames.length, e.frames, f + ': frames');
      ok(ds.meta.needsTopology && ds.meta.natoms === e.atoms, f + ': coordinates waiting for a topology of ' + e.atoms + ' atoms');
      ds.frames.forEach((fr, k) => { same(fr.xyz, e.xyz[k], f + ': frame ' + (k + 1)); ok(fr.Z === null, f + ': no elements in a trajectory'); });
      if (e.time.length) ok(ds.frames.every((fr, k) => fr.time === e.time[k]), f + ': times ' + ds.frames.map(x => x.time).join());
      if (e.cell.length) {
        const cells = ds.frames.map(fr => fr.cell || ds.cell);
        ok(cells.every((c, k) => c && closeAll(cellOf(c), e.cell[k], 1e-6)), f + ': box lengths and angles of every frame');
      } else ok(!ds.cell && ds.frames.every(fr => !fr.cell), f + ': no box');
    }
  });
  await test('Amber attributes and variables', () => {
    const ds = BM.parseFileBinary('traj_64bit.nc', bytes('traj_64bit.nc'));
    eq(ds.program, 'sander 18.0', 'program and version');
    eq(ds.title, 'BlueMolecule NetCDF test (written by SciPy)', 'title');
    ok(ds.meta.hasVelocities, 'velocities are noticed (and skipped: they sit between the time and the coordinates)');
    ok(ds.frames[0].cell && ds.frames[1].cell, 'a box that changes stays on each frame');
    ok(/^Frame 3 · t = 6 ps$/.test(ds.frames[2].label), 'frame label with the time: ' + ds.frames[2].label);
    const r = BM.parseFileBinary('restart.ncrst', bytes('restart.ncrst'));
    eq(r.frames[0].label, 'Restart', 'a restart is one geometry');
    eq(r.frames[0].time, 1234.5, 'restart time');
    ok(r.meta.hasVelocities, 'restart velocities');
    const h = BM.structReaders.netcdfHeader(bytes('traj_cdf5.nc'));
    eq(h.version, 5, 'CDF-5 header');
    eq(h.gatts.big_number, Math.pow(2, 40) + 5, 'a 64-bit integer attribute');
    eq(h.gatts.unsigned, 4000000000, 'an unsigned attribute');
    eq(h.gatts.ratio, 0.25, 'a double attribute');
    eq(BM.structReaders.netcdfHeader(bytes('traj_classic.nc')).version, 1, 'classic header');
  });
  await test('an Amber restart named .rst7 (pmemd writes NetCDF restarts by default)', () => {
    const u8 = bytes('restart_netcdf.rst7');
    ok(BM.maybeBinary('md.rst7') && BM.maybeBinary('md.restrt') && !BM.maybeBinary('md.log'), 'rst7 and restrt may hold NetCDF');
    ok(BM.sniffBinary('md.rst7', u8) && !BM.sniffBinary('md.rst7', new TextEncoder().encode('title\n    4\n')), 'the NetCDF signature decides');
    ok(!BM.isBinaryFormat('md.rst7'), 'the extension stays a text one');
    const ds = BM.parseFileBinary('md.rst7', u8);
    eq(ds.format, 'netcdf', 'read by the NetCDF reader');
    same(ds.frames[0].xyz, EXP['restart_netcdf.rst7'].xyz[0], 'coordinates');
    const txt = BM.parseFile('w.inpcrd', fs.readFileSync(path.join(__dirname, 'fixtures', 'struct', 'water.inpcrd'), 'utf8'));
    eq(txt.format, 'inpcrd', 'an ASCII restart still reads as text');
  });
  await test('slices of a big file: stride, range, last N, probe, reload', async () => {
    const u8 = bytes('traj_64bit.nc'), e = EXP['traj_64bit.nc'];
    for (const src of [blob(u8), chunks(u8)]) {
      const all = await BM.parseSource('t.nc', src);
      eq(all.frames.length, 12, 'all frames');
      all.frames.forEach((fr, k) => same(fr.xyz, e.xyz[k], 'frame ' + (k + 1)));
      ok(all.frames.every((fr, k) => fr.time === e.time[k]), 'times');
      const pr = await BM.probeSource('t.nc', src);
      ok(pr.frames === 12 && pr.atoms === 4 && pr.exact === true && pr.perFrame === 96, 'probe: ' + JSON.stringify(pr));
      const st = await BM.parseSource('t.nc', src, { stride: 5 });
      eq(st.frames.map(f => f.label.split(' ·')[0]).join(), 'Frame 1,Frame 6,Frame 11', 'every 5th');
      same(st.frames[2].xyz, e.xyz[10], 'every 5th: frame 11');
      eq(st.meta.totalFrames, 12, 'total frames');
      eq(st.meta.stride, 5, 'stride');
      ok(st.warnings.some(w => /Loaded 3 of 12 frames \(every 5th\)/.test(w)), 'says what was loaded: ' + st.warnings.join(' | '));
      const rg = await BM.parseSource('t.nc', src, { first: 3, last: 6 });
      eq(rg.frames.length, 4, 'range 4-7');
      same(rg.frames[0].xyz, e.xyz[3], 'range: first frame');
      const tl = await BM.parseSource('t.nc', src, { tail: 3 });
      eq(tl.frames.length, 3, 'the last 3');
      same(tl.frames[2].xyz, e.xyz[11], 'the last frame');
      eq(tl.meta.first, 9, 'last 3: first index');
      const rl = await st.meta.reload({ stride: 1, first: 10 });
      eq(rl.frames.length, 2, 'reload with another choice');
      same(rl.frames[1].xyz, e.xyz[11], 'reload: the last frame');
      const auto = await BM.parseSource('t.nc', src, { budget: 96 * 4 });
      ok(auto.frames.length <= 4 && auto.meta.stride === 3, 'a stride from the memory budget: ' + auto.meta.stride);
      const r = await BM.parseSource('r.ncrst', blob(bytes('restart.ncrst')));
      same(r.frames[0].xyz, EXP['restart.ncrst'].xyz[0], 'a restart from slices');
      ok(r.cell && closeAll(cellOf(r.cell), [40, 41, 42, 90, 90, 90], 1e-9), 'restart box from slices');
    }
    const ctl = new AbortController(); ctl.abort();
    let err = null;
    try { await BM.parseSource('t.nc', blob(u8), { signal: ctl.signal }); } catch (x) { err = x; }
    ok(err && err.name === 'AbortError', 'cancel');
  });
  await test('cut short, streaming, a header that lags, NetCDF 4, text', async () => {
    const u8 = bytes('traj_64bit.nc');
    const cut = u8.slice(0, u8.length - 30);
    const c1 = BM.parseFileBinary('c.nc', cut);
    eq(c1.frames.length, 11, 'a cut last frame is skipped');
    ok(c1.warnings.some(w => /counts 12 frames; the file ends after 11/.test(w)), 'and the header count is reported: ' + c1.warnings.join(' | '));
    eq((await BM.parseSource('c.nc', blob(cut))).frames.length, 11, 'the same from slices');
    const lag = withNumrecs(u8, 5);
    const c2 = BM.parseFileBinary('l.nc', lag);
    eq(c2.frames.length, 12, 'a header that counts fewer frames than the file holds (a run that was stopped)');
    ok(c2.warnings.some(w => /counts 5 frames; the file holds 12 complete ones/.test(w)), 'is reported');
    const stream = withNumrecs(u8, 0xFFFFFFFF);
    eq(BM.parseFileBinary('s.nc', stream).frames.length, 12, 'streaming: frames from the file size');
    eq((await BM.probeSource('s.nc', blob(stream))).frames, 12, 'streaming probe');
    let err = null;
    try { BM.parseFileBinary('h.nc', u8.slice(0, 200)); } catch (x) { err = x; }
    ok(err && /header is incomplete/.test(err.message), 'a cut header: ' + (err && err.message));
    const hdf = new Uint8Array(64); hdf.set([0x89, 0x48, 0x44, 0x46, 0x0d, 0x0a, 0x1a, 0x0a]);
    err = null; try { BM.parseFileBinary('n4.nc', hdf); } catch (x) { err = x; }
    ok(err && /NetCDF 4/.test(err.message) && /nccopy/.test(err.message), 'NetCDF 4 is refused with advice: ' + (err && err.message));
    err = null; try { BM.parseFileBinary('old.mdcrd', new TextEncoder().encode('title\n   1.000   2.000   3.000\n')); } catch (x) { err = x; }
    ok(err && /text file/.test(err.message) && /cpptraj/.test(err.message), 'an ASCII mdcrd: ' + (err && err.message));
    const rst = bytes('restart.ncrst');
    for (const cutAt of [rst.length - 8, rst.length - 60]) {
      err = null; try { BM.parseFileBinary('r.ncrst', rst.slice(0, cutAt)); } catch (x) { err = x; }
      ok(err && /restart is cut short/.test(err.message), 'a cut restart (' + cutAt + ' bytes): ' + (err && err.message));
      err = null; try { await BM.parseSource('r.ncrst', blob(rst.slice(0, cutAt))); } catch (x) { err = x; }
      ok(err && /restart is cut short/.test(err.message), 'a cut restart from slices: ' + (err && err.message));
    }
    err = null; try { BM.parseFileBinary('x.nc', new Uint8Array(100)); } catch (x) { err = x; }
    ok(err && /Not a NetCDF file/.test(err.message), 'not NetCDF');
    // a NetCDF file that is not an Amber trajectory
    const other = u8.slice();
    const at = Buffer.from(other).indexOf('coordinates');
    other[at] = 'k'.charCodeAt(0);
    err = null; try { BM.parseFileBinary('o.nc', other); } catch (x) { err = x; }
    ok(err && /No "coordinates" variable/.test(err.message), 'no coordinates variable: ' + (err && err.message));
  });
  await test('registration', () => {
    const f = BM.formats.find(x => x.id === 'netcdf');
    ok(f && f.binary && ['nc', 'ncdf', 'ncrst', 'mdcrd'].every(x => f.extensions.indexOf(x) >= 0), 'listed as a binary trajectory format');
    ok(BM.isBinaryFormat('run.nc') && BM.hasSourceReader('run.nc') && BM.hasSourceReader('mdcrd.mdcrd'), 'read from slices when big');
  });

  console.log('\nBlueMolecule NetCDF tests: ' + pass + ' passed, ' + fail + ' failed');
  if (fail) { for (const f of failures) console.log('  FAIL ' + f); process.exit(1); }
})();
