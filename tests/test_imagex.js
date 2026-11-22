#!/usr/bin/env node
/* BlueMolecule tests for src/imagex.js (issue #16), BlueMolecule's own image encoders: deflate and zlib (checked with
 * Node's zlib), PNG, TIFF and PDF decoded here, JPEG and WebP checked by structure, ICC profiles, sRGB -> Display
 * P3. When Python 3 with Pillow is installed, every format is also decoded by Pillow (exact pixels for the lossless
 * formats, PSNR for JPEG); the browser tests decode PNG, JPEG and WebP with Chromium's own decoders.
 * Run:  node tests/test_imagex.js */
'use strict';
const path = require('path');
const fs = require('fs');
const os = require('os');
const zlibN = require('zlib');
const { spawnSync } = require('child_process');
require(path.join(__dirname, '..', 'src', 'imagex.js'));
const X = globalThis.BM.imagex;

let pass = 0, fail = 0, current = '';
const failures = [];
function ok(c, msg) { if (c) pass++; else { fail++; failures.push('[' + current + '] ' + msg); } }
function eq(a, b, msg) { ok(a === b, msg + ': got ' + JSON.stringify(a) + ', expected ' + JSON.stringify(b)); }
function test(name, fn) { current = name; try { fn(); } catch (e) { fail++; failures.push('[' + name + '] threw: ' + (e && e.stack || e)); } }
const same = (a, b) => a.length === b.length && Buffer.compare(Buffer.from(a), Buffer.from(b)) === 0;
let seed = 7;
const rnd = () => { seed = (Math.imul(seed, 1103515245) + 12345) >>> 0; return seed / 4294967296; };
const u32be = (b, p) => ((b[p] << 24) | (b[p + 1] << 16) | (b[p + 2] << 8) | b[p + 3]) >>> 0;
const u16le = (b, p) => b[p] | (b[p + 1] << 8), u32le = (b, p) => (b[p] | (b[p + 1] << 8) | (b[p + 2] << 16) | (b[p + 3] << 24)) >>> 0;
const ascii = (b, p, n) => String.fromCharCode(...b.subarray(p, p + n));

/** A test picture: transparent background, shaded discs with soft edges, a gradient band, a noisy corner. */
function picture(W, H, opaque) {
  const d = new Uint8Array(W * H * 4);
  for (let y = 0; y < H; y++) for (let x = 0; x < W; x++) {
    let r = 0, g = 0, b = 0, a = 0;
    if (y < H * 0.15) { r = x * 255 / W; g = 128; b = 255 - r; a = 255; }
    for (const [cx, cy, rad, col] of [[0.3, 0.55, 0.25, [230, 40, 40]], [0.66, 0.6, 0.3, [60, 120, 250]], [0.5, 0.4, 0.12, [250, 250, 250]]]) {
      const dd = Math.hypot(x - cx * W, y - cy * H), R = rad * H;
      if (dd < R) { const s = Math.sqrt(1 - (dd / R) ** 2); r = col[0] * (0.3 + 0.7 * s); g = col[1] * (0.3 + 0.7 * s); b = col[2] * (0.3 + 0.7 * s); a = 255; }
      else if (dd < R + 1.5) { a = Math.max(a, Math.round(255 * (R + 1.5 - dd) / 1.5)); r = col[0] * 0.3; g = col[1] * 0.3; b = col[2] * 0.3; }
    }
    if (x > W * 0.86 && y > H * 0.74) { r = rnd() * 255; g = rnd() * 255; b = rnd() * 255; a = 255; }
    const i = (y * W + x) * 4;
    d[i] = r; d[i + 1] = g; d[i + 2] = b; d[i + 3] = opaque ? 255 : a;
  }
  return { width: W, height: H, data: d };
}
const IMG = picture(203, 131, false), OPAQUE = picture(203, 131, true);
const rgbOf = img => { const o = new Uint8Array(img.width * img.height * 3); for (let i = 0, j = 0; i < img.data.length; i += 4) { o[j++] = img.data[i]; o[j++] = img.data[i + 1]; o[j++] = img.data[i + 2]; } return o; };

test('checksums', () => {
  const t = new TextEncoder().encode('123456789');
  eq(X.crc32(t), 0xcbf43926, 'CRC-32 of "123456789"');
  eq(X.adler32(new TextEncoder().encode('Wikipedia')), 0x11e60398, 'Adler-32 of "Wikipedia"');
  eq(X.adler32(new Uint8Array(0)), 1, 'Adler-32 of nothing');
});

test('Huffman code lengths', () => {
  const fib = [1, 1];
  while (fib.length < 30) fib.push(fib[fib.length - 1] + fib[fib.length - 2]);   // the worst case for code length
  for (const limit of [7, 15, 16]) {
    const l = X.codeLengths(fib, limit);
    let kraft = 0, max = 0;
    for (const v of l) if (v) { kraft += Math.pow(2, -v); max = Math.max(max, v); }
    ok(max <= limit && Math.abs(kraft - 1) < 1e-12, 'Fibonacci frequencies, limit ' + limit + ': max ' + max + ', Kraft sum ' + kraft);
  }
  const l = X.codeLengths([5, 0, 0, 9, 1], 15);
  ok(l[1] === 0 && l[2] === 0 && l[3] <= l[0] && l[0] <= l[4], 'unused symbols get no code; frequent ones the shortest: ' + Array.from(l));
  eq(Array.from(X.codeLengths([0, 3, 0], 15)).join(), '0,1,0', 'one symbol: length 1');
  const c = X.canonical(Uint8Array.from([2, 1, 3, 3]));
  eq(Array.from(c).join(), '2,0,6,7', 'canonical codes (RFC 1951 example order)');
});

test('deflate and zlib', () => {
  const cases = {
    empty: new Uint8Array(0), one: new Uint8Array([42]), zeros: new Uint8Array(300000),
    text: new TextEncoder().encode(Array.from({ length: 20000 }, (_, i) => 'atom ' + (i % 97) + ' x=' + (rnd() * 10).toFixed(3)).join('\n')),
    random: Uint8Array.from({ length: 200000 }, () => (rnd() * 256) | 0),
    mixed: Uint8Array.from({ length: 600000 }, (_, i) => (i % 9973 < 4000 ? 0 : i % 3000 < 1500 ? (rnd() * 4) | 0 : (rnd() * 256) | 0)),
    image: (() => { const W = 3000, H = 300, d = new Uint8Array((W * 4 + 1) * H); for (let y = 0; y < H; y++) for (let x = 0; x < W; x++) { const i = y * (W * 4 + 1) + 1 + x * 4, inDisc = Math.hypot(x - 1500, y - 150) < 120 + (y % 7); d[i] = inDisc ? (x * 7) & 255 : 16; d[i + 1] = inDisc ? (y * 3) & 255 : 23; d[i + 2] = inDisc ? 200 : 34; d[i + 3] = 255; } return d; })()
  };
  for (const [k, d] of Object.entries(cases)) for (const lv of [0, 1, 4, 6, 9]) {
    let back = null;
    try { back = zlibN.inflateSync(Buffer.from(X.zlib(d, lv))); } catch (e) { ok(false, k + ' level ' + lv + ': ' + e.message); continue; }
    ok(same(back, d), k + ' level ' + lv + ': the zlib stream inflates to the data');
  }
  ok(X.zlib(cases.zeros, 6).length < 1000 && X.zlib(cases.text, 9).length < cases.text.length / 3, 'compresses: zeros and text');
  ok(X.zlib(cases.random, 6).length <= cases.random.length + 64, 'random data: stored blocks, barely larger');
  ok(same(X.zlib(cases.image, 6), X.zlib(cases.image, 6)), 'the same bytes every time');
  const raw = X.deflate(cases.text, 6);
  ok(same(zlibN.inflateRawSync(Buffer.from(raw)), cases.text), 'raw deflate');
});

/** PNG: chunks with their CRC, then the pixels inflated and unfiltered. */
function readPng(b) {
  ok(same(b.subarray(0, 8), [137, 80, 78, 71, 13, 10, 26, 10]), 'PNG signature');
  const chunks = [];
  let idat = [];
  for (let p = 8; p < b.length;) {
    const len = u32be(b, p), type = ascii(b, p + 4, 4), data = b.subarray(p + 8, p + 8 + len);
    ok(X.crc32(b, p + 4, p + 8 + len) === u32be(b, p + 8 + len), type + ' chunk CRC');
    chunks.push({ type, data });
    if (type === 'IDAT') idat.push(data);
    p += 12 + len;
  }
  const ih = chunks[0].data, W = u32be(ih, 0), H = u32be(ih, 4), ct = ih[9], ch = ct === 6 ? 4 : 3, stride = W * ch;
  const raw = zlibN.inflateSync(Buffer.concat(idat.map(x => Buffer.from(x))));
  const out = new Uint8Array(W * H * 4), prev = new Uint8Array(stride), row = new Uint8Array(stride);
  for (let y = 0; y < H; y++) {
    const f = raw[y * (stride + 1)];
    for (let i = 0; i < stride; i++) {
      const x = raw[y * (stride + 1) + 1 + i], a = i >= ch ? row[i - ch] : 0, up = prev[i], c = i >= ch ? prev[i - ch] : 0;
      let v;
      if (f === 0) v = x; else if (f === 1) v = x + a; else if (f === 2) v = x + up; else if (f === 3) v = x + ((a + up) >> 1);
      else { const p = a + up - c, pa = Math.abs(p - a), pb = Math.abs(p - up), pc = Math.abs(p - c); v = x + (pa <= pb && pa <= pc ? a : pb <= pc ? up : c); }
      row[i] = v & 255;
    }
    for (let x = 0; x < W; x++) for (let k = 0; k < 4; k++) out[(y * W + x) * 4 + k] = k < ch ? row[x * ch + k] : 255;
    prev.set(row);
  }
  return { W, H, ct, chunks, data: out };
}
test('PNG', () => {
  for (const level of [0, 1, 6, 9]) {
    const r = readPng(X.png(IMG, { level }));
    ok(r.W === 203 && r.H === 131 && r.ct === 6 && same(r.data, IMG.data), 'level ' + level + ': RGBA pixels read back exactly');
  }
  const o = readPng(X.png(OPAQUE, { dpi: 300, profile: 'srgb' }));
  ok(o.ct === 2 && same(o.data, OPAQUE.data), 'an opaque image is written as RGB');
  const phys = o.chunks.find(c => c.type === 'pHYs');
  ok(phys && u32be(phys.data, 0) === 11811 && phys.data[8] === 1, '300 DPI = 11811 pixels per metre');
  ok(o.chunks.some(c => c.type === 'sRGB') && o.chunks.some(c => c.type === 'gAMA') && !o.chunks.some(c => c.type === 'iCCP'), 'sRGB: the sRGB chunk');
  const p3 = readPng(X.png(IMG, { profile: 'p3' })), iccp = p3.chunks.find(c => c.type === 'iCCP');
  const z = iccp.data.indexOf(0);
  ok(ascii(iccp.data, 0, z) === 'Display P3' && same(zlibN.inflateSync(Buffer.from(iccp.data.subarray(z + 2))), X.icc('p3')), 'Display P3: the profile in iCCP');
  const t = o.chunks.find(c => c.type === 'tEXt');
  ok(t && ascii(t.data, 0, 8) === 'Software', 'a Software text chunk');
  ok(same(X.png(IMG, {}), X.png(IMG, {})), 'the same bytes every time');
  ok(same(readPng(X.png(IMG, { alpha: false })).data, OPAQUE.data.map((v, i) => (i % 4 === 3 ? 255 : IMG.data[i]))), 'alpha: false drops the alpha channel');
});

test('ICC profiles', () => {
  for (const name of ['srgb', 'p3']) {
    const p = X.icc(name), s15 = o => ((u32be(p, o) | 0) / 65536);
    eq(u32be(p, 0), p.length, name + ': size in the header');
    ok(ascii(p, 36, 4) === 'acsp' && ascii(p, 12, 4) === 'mntr' && ascii(p, 16, 4) === 'RGB ' && ascii(p, 20, 4) === 'XYZ ' && u32be(p, 8) === 0x02100000, name + ': a v2.1 RGB display profile');
    const n = u32be(p, 128), tags = {};
    for (let k = 0; k < n; k++) { const q = 132 + 12 * k; tags[ascii(p, q, 4)] = [u32be(p, q + 4), u32be(p, q + 8)]; }
    ok(['desc', 'cprt', 'wtpt', 'rXYZ', 'gXYZ', 'bXYZ', 'rTRC', 'gTRC', 'bTRC', 'chad'].every(t => tags[t] && tags[t][0] + tags[t][1] <= p.length && tags[t][0] % 4 === 0), name + ': every tag, inside the profile, aligned');
    const xyz = t => [0, 1, 2].map(k => s15(tags[t][0] + 8 + 4 * k));
    const sum = [0, 1, 2].map(k => xyz('rXYZ')[k] + xyz('gXYZ')[k] + xyz('bXYZ')[k]);
    ok(sum.every((v, k) => Math.abs(v - [0.9642, 1, 0.8249][k]) < 2e-4), name + ': the colorants add up to the D50 white: ' + sum.map(v => v.toFixed(4)));
    ok(tags.rTRC[0] === tags.gTRC[0] && tags.gTRC[0] === tags.bTRC[0], name + ': one tone curve shared');
  }
  const red = p => { const n = u32be(p, 128); for (let k = 0; k < n; k++) if (ascii(p, 132 + 12 * k, 4) === 'rXYZ') return (u32be(p, u32be(p, 136 + 12 * k) + 8) | 0) / 65536; return 0; };
  ok(Math.abs(red(X.icc('srgb')) - 0.4361) < 2e-4 && Math.abs(red(X.icc('p3')) - 0.5151) < 2e-4, 'red colorants of sRGB (0.4361) and Display P3 (0.5151)');
  eq(X.icc('cmyk'), null, 'unknown profile: null');
});

test('sRGB to Display P3', () => {
  const px = { width: 4, height: 1, data: Uint8Array.from([255, 255, 255, 255, 0, 0, 0, 128, 255, 0, 0, 255, 0, 255, 0, 255]) };
  const q = X.toP3(px).data;
  ok(q[0] === 255 && q[1] === 255 && q[2] === 255 && q[4] === 0 && q[5] === 0 && q[6] === 0 && q[7] === 128, 'white and black stay; alpha is kept');
  ok(Math.abs(q[8] - 234) <= 1 && Math.abs(q[9] - 51) <= 1 && Math.abs(q[10] - 35) <= 1, 'sRGB red is (234, 51, 35) in Display P3: ' + Array.from(q.subarray(8, 11)));
  ok(Math.abs(q[12] - 117) <= 1 && Math.abs(q[13] - 251) <= 1 && Math.abs(q[14] - 76) <= 1, 'sRGB green is (117, 251, 76): ' + Array.from(q.subarray(12, 15)));
});

/** JPEG markers: [{m, data}] up to the start of scan. */
function jpegMarkers(b) {
  const out = [];
  ok(b[0] === 0xff && b[1] === 0xd8 && b[b.length - 2] === 0xff && b[b.length - 1] === 0xd9, 'JPEG: SOI ... EOI');
  for (let p = 2; p < b.length;) {
    const m = b[p + 1], len = (b[p + 2] << 8) | b[p + 3];
    out.push({ m, data: b.subarray(p + 4, p + 2 + len) });
    if (m === 0xda) break;
    p += 2 + len;
  }
  return out;
}
test('JPEG', () => {
  const b = X.jpeg(OPAQUE, { quality: 50, dpi: 300, profile: 'srgb' }), mk = jpegMarkers(b), get = m => mk.filter(x => x.m === m);
  const app0 = get(0xe0)[0];
  ok(app0 && ascii(app0.data, 0, 4) === 'JFIF' && app0.data[7] === 1 && ((app0.data[8] << 8) | app0.data[9]) === 300, 'JFIF with 300 DPI');
  const app2 = get(0xe2)[0];
  ok(app2 && ascii(app2.data, 0, 11) === 'ICC_PROFILE' && same(app2.data.subarray(14), X.icc('srgb')), 'the ICC profile in APP2');
  const dqt = get(0xdb)[0];
  ok(dqt && dqt.data.length === 130 && dqt.data[1] === 16 && dqt.data[2] === 11 && dqt.data[66] === 17, 'quality 50: the standard tables (16, 11, ... and 17, ...)');
  const q100 = jpegMarkers(X.jpeg(OPAQUE, { quality: 100 })).find(x => x.m === 0xdb);
  ok(q100.data.subarray(1, 65).every(v => v === 1), 'quality 100: every quantiser 1');
  const sof = get(0xc0)[0];
  ok(sof && ((sof.data[1] << 8) | sof.data[2]) === 131 && ((sof.data[3] << 8) | sof.data[4]) === 203 && sof.data[5] === 3 && sof.data[7] === 0x11, 'baseline, 203 x 131, three components, no subsampling');
  const dht = get(0xc4)[0];
  let p = 0, tables = 0, okT = true;
  while (dht && p < dht.data.length) {
    const counts = dht.data.subarray(p + 1, p + 17), n = counts.reduce((s, v) => s + v, 0);
    let kraft = 0; counts.forEach((c, k) => { kraft += c * Math.pow(2, -(k + 1)); });
    if (!(n > 0 && kraft < 1)) okT = false;   // < 1: no code of all ones
    p += 17 + n; tables++;
  }
  ok(tables === 4 && okT, 'four Huffman tables made for the image, none with a code of all ones');
  ok(same(X.jpeg(OPAQUE, { quality: 80 }), X.jpeg(OPAQUE, { quality: 80 })), 'the same bytes every time');
  ok(X.jpeg(OPAQUE, { quality: 30 }).length < X.jpeg(OPAQUE, { quality: 95 }).length, 'lower quality, smaller file');
  let thrown = '';
  try { X.jpeg({ width: 70000, height: 1, data: new Uint8Array(280000) }); } catch (e) { thrown = e.message; }
  ok(/at most 65,535/.test(thrown), 'too large for JPEG: ' + thrown);
});

test('WebP', () => {
  const b = X.webp(IMG, {});
  ok(ascii(b, 0, 4) === 'RIFF' && u32le(b, 4) === b.length - 8 && ascii(b, 8, 4) === 'WEBP' && ascii(b, 12, 4) === 'VP8L', 'RIFF / WEBP / VP8L');
  const v = 20;
  const bits = u32le(b, v + 1);
  ok(b[v] === 0x2f && (bits & 0x3fff) + 1 === 203 && ((bits >>> 14) & 0x3fff) + 1 === 131 && ((bits >>> 28) & 1) === 1, 'lossless signature, 203 x 131, alpha used');
  const o = X.webp(OPAQUE, {});
  ok(((u32le(o, v + 1) >>> 28) & 1) === 0, 'no alpha in an opaque image');
  const ic = X.webp(IMG, { profile: 'srgb' });
  ok(ascii(ic, 12, 4) === 'VP8X' && (ic[20] & 0x20) && (ic[20] & 0x10) && ascii(ic, 30, 4) === 'ICCP' && same(ic.subarray(38, 38 + u32le(ic, 34)), X.icc('srgb')), 'VP8X with the ICC and alpha flags, then ICCP');
  ok(same(X.webp(IMG, {}), b), 'the same bytes every time');
  let thrown = '';
  try { X.webp({ width: 16385, height: 1, data: new Uint8Array(16385 * 4) }); } catch (e) { thrown = e.message; }
  ok(/at most 16,384/.test(thrown), 'too large for WebP: ' + thrown);
});

/** Baseline TIFF: the tags and the pixels (strips inflated, the predictor undone). */
function readTiff(b) {
  ok(ascii(b, 0, 2) === 'II' && u16le(b, 2) === 42, 'TIFF: little-endian, 42');
  const ifd = u32le(b, 4), n = u16le(b, ifd), tags = {};
  const size = { 2: 1, 3: 2, 4: 4, 5: 8, 7: 1 };
  for (let k = 0; k < n; k++) {
    const q = ifd + 2 + 12 * k, tag = u16le(b, q), type = u16le(b, q + 2), count = u32le(b, q + 4), bytes = size[type] * count;
    const at = bytes > 4 ? u32le(b, q + 8) : q + 8, vals = [];
    for (let i = 0; i < count; i++) vals.push(type === 3 ? u16le(b, at + 2 * i) : type === 4 ? u32le(b, at + 4 * i) : type === 5 ? u32le(b, at + 8 * i) / u32le(b, at + 8 * i + 4) : b[at + i]);
    tags[tag] = { type, vals, at, count };
  }
  const W = tags[256].vals[0], H = tags[257].vals[0], spp = tags[277].vals[0], rows = tags[278].vals[0], comp = tags[259].vals[0], pred = tags[317] ? tags[317].vals[0] : 1;
  const out = new Uint8Array(W * H * 4);
  tags[273].vals.forEach((off, k) => {
    let s = b.subarray(off, off + tags[279].vals[k]);
    if (comp === 8) s = new Uint8Array(zlibN.inflateSync(Buffer.from(s)));
    const y0 = k * rows;
    for (let y = y0, t = 0; y < Math.min(H, y0 + rows); y++) {
      const row = s.subarray((y - y0) * W * spp, (y - y0 + 1) * W * spp);
      if (pred === 2) for (let i = spp; i < row.length; i++) row[i] = (row[i] + row[i - spp]) & 255;
      for (let x = 0; x < W; x++) for (let c = 0; c < 4; c++) out[(y * W + x) * 4 + c] = c < spp ? row[x * spp + c] : 255;
      t++;
    }
  });
  return { W, H, spp, comp, pred, tags, data: out };
}
test('TIFF', () => {
  const t = readTiff(X.tiff(IMG, { dpi: 300, profile: 'srgb' }));
  ok(t.W === 203 && t.H === 131 && t.spp === 4 && t.comp === 8 && t.pred === 2 && same(t.data, IMG.data), 'RGBA, Deflate with the horizontal predictor: pixels read back exactly');
  ok(t.tags[282].vals[0] === 300 && t.tags[296].vals[0] === 2 && t.tags[338].vals[0] === 2, '300 DPI (inches); unassociated alpha');
  const icc = t.tags[34675];
  ok(icc && same(X.tiff(IMG, { profile: 'srgb' }).subarray(icc.at, icc.at + icc.count), X.icc('srgb')), 'the ICC profile tag');
  const r = readTiff(X.tiff(OPAQUE, { level: 0 }));
  ok(r.spp === 3 && r.comp === 1 && same(r.data, OPAQUE.data), 'uncompressed RGB');
  const big = picture(700, 500, true), rb = readTiff(X.tiff(big, {}));
  ok(rb.tags[273].vals.length > 1 && same(rb.data, big.data), 'several strips: ' + rb.tags[273].vals.length);
});

test('PDF', () => {
  const b = X.pdf(IMG, { dpi: 300, profile: 'srgb', title: 'water (test)' }), s = Buffer.from(b).toString('latin1');
  ok(s.startsWith('%PDF-1.4') && s.trimEnd().endsWith('%%EOF'), 'a PDF 1.4 file');
  const sx = +/startxref\n(\d+)/.exec(s)[1], xr = s.slice(sx);
  ok(xr.startsWith('xref'), 'startxref points at the xref table');
  const offs = [...xr.matchAll(/^(\d{10}) 00000 n $/gm)].map(m => +m[1]);
  ok(offs.length >= 6 && offs.every((o, k) => s.startsWith((k + 1) + ' 0 obj', o)), 'every object at its xref offset (' + offs.length + ')');
  const mb = /\/MediaBox \[0 0 ([\d.]+) ([\d.]+)\]/.exec(s);
  ok(mb && Math.abs(+mb[1] - 203 / 300 * 72) < 0.01 && Math.abs(+mb[2] - 131 / 300 * 72) < 0.01, 'page size from 300 DPI: ' + (mb && mb.slice(1)));
  const stream = re => { const m = re.exec(s); const len = +/\/Length (\d+)/.exec(m[0])[1], at = s.indexOf('stream\n', m.index) + 7; return zlibN.inflateSync(Buffer.from(b.subarray(at, at + len))); };
  ok(same(stream(/\/Subtype \/Image \/Width 203 \/Height 131 \/ColorSpace \[\/ICCBased[^>]*/), rgbOf(IMG)), 'the image: RGB, ICCBased');
  ok(same(stream(/\/ColorSpace \/DeviceGray[^>]*/), IMG.data.filter((_, i) => i % 4 === 3)), 'the transparency as an SMask');
  ok(same(stream(/\/N 3 \/Alternate \/DeviceRGB[^>]*/), X.icc('srgb')), 'the ICC profile');
  ok(/\/Title \(water \\\(test\\\)\)/.test(s), 'the title, escaped');
  const o = Buffer.from(X.pdf(OPAQUE, {})).toString('latin1');
  ok(!/SMask/.test(o) && /\/DeviceRGB/.test(o) && /\/MediaBox \[0 0 203 131\]/.test(o), 'opaque, no profile, 72 DPI: no SMask, DeviceRGB, 1 pt per pixel');
});

test('encode', () => {
  ok(same(X.encode('png', IMG, {}), X.png(IMG, {})), 'encode("png") is png()');
  ok(same(X.encode('tiff', IMG, { profile: 'p3' }), X.tiff(X.toP3(IMG), { profile: 'p3' })), 'the P3 profile converts the pixels first');
  let thrown = '';
  try { X.encode('bmp', IMG, {}); } catch (e) { thrown = e.message; }
  ok(/Unknown image format/.test(thrown), 'unknown format');
  ok(X.FORMATS.map(f => f.id).join() === 'png,jpeg,webp,tiff,pdf' && !X.FORMATS.find(f => f.id === 'jpeg').alpha, 'five formats; JPEG has no transparency');
});

// ------------------------------------------------------------------ Pillow, when it is installed
test('decoded by Pillow', () => {
  const probe = spawnSync('python3', ['-c', 'import PIL, PIL.features; print(PIL.features.check("webp"))'], { encoding: 'utf8' });
  if (probe.status !== 0) { console.log('  (Pillow not installed: the decoding by Pillow is skipped)'); return; }
  const dir = fs.mkdtempSync(path.join(os.tmpdir(), 'bm-imagex-'));
  const files = [['a.png', X.png(IMG, { dpi: 300 })], ['o.png', X.png(OPAQUE, {})], ['a.tif', X.tiff(IMG, { profile: 'srgb' })], ['o.tif', X.tiff(OPAQUE, { level: 0 })],
    ['a.jpg', X.jpeg(OPAQUE, { quality: 92, profile: 'srgb' })], ['q50.jpg', X.jpeg(OPAQUE, { quality: 50 })]];
  const big = picture(640, 480, false);
  if (probe.stdout.trim() === 'True') files.push(['a.webp', X.webp(IMG, {})], ['o.webp', X.webp(OPAQUE, { profile: 'srgb' })], ['big.webp', X.webp(big, {})]);
  for (const [n, data] of files) fs.writeFileSync(path.join(dir, n), data);
  fs.writeFileSync(path.join(dir, 'a.rgba'), IMG.data); fs.writeFileSync(path.join(dir, 'o.rgba'), OPAQUE.data); fs.writeFileSync(path.join(dir, 'big.rgba'), big.data);
  const py = `
import sys, math, os
from PIL import Image
d = sys.argv[1]
ref = {k: open(os.path.join(d, k + '.rgba'), 'rb').read() for k in ('a', 'o', 'big')}
for n in sorted(os.listdir(d)):
    if n.endswith('.rgba'): continue
    im = Image.open(os.path.join(d, n)); im.load()
    key = 'big' if n.startswith('big') else ('o' if n.startswith('o') or n.endswith('.jpg') else 'a')
    got = im.convert('RGBA').tobytes(); r = ref[key]
    if n.endswith('.jpg'):
        se = sum((got[i] - r[i]) ** 2 for i in range(len(r)) if i % 4 != 3); mse = se / (len(r) * 3 / 4)
        print(n, 'psnr', round(10 * math.log10(255 ** 2 / mse), 1) if mse else 99, int('icc_profile' in im.info))
    else:
        print(n, 'diff', sum(1 for a, b in zip(got, r) if a != b), int('icc_profile' in im.info))
`;
  const r = spawnSync('python3', ['-c', py, dir], { encoding: 'utf8' });
  if (r.status !== 0) { ok(false, 'Pillow failed: ' + r.stderr.slice(-400)); return; }
  for (const line of r.stdout.trim().split('\n')) {
    const [n, kind, v, icc] = line.split(' ');
    if (kind === 'psnr') ok(+v > (n === 'q50.jpg' ? 28 : 38), n + ': Pillow decodes it, PSNR ' + v + ' dB');
    else ok(v === '0', n + ': Pillow reads exactly the same pixels (' + v + ' bytes differ)');
    if (/^a\.(tif|jpg)$|^o\.webp$/.test(n)) ok(icc === '1', n + ': Pillow finds the ICC profile');
  }
  fs.rmSync(dir, { recursive: true, force: true });
});

console.log('\nBlueMolecule image encoder tests: ' + pass + ' passed, ' + fail + ' failed');
if (fail) { for (const f of failures) console.log('  FAIL ' + f); process.exit(1); }
