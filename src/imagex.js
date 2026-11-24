/* BlueMolecule — image encoders (issue #16), BM.imagex. BlueMolecule writes every format itself, so an image
 * is the same file in every browser engine (canvas.toBlob differs between engines, and WebKit cannot write WebP).
 *   png(img, opts)    RGBA or RGB, adaptive filters, zlib level 0-9, pHYs (DPI), sRGB chunk or an ICC profile
 *   jpeg(img, opts)   baseline, no chroma subsampling, quality 1-100, Huffman tables made for the image, JFIF
 *                     density (DPI), ICC profile in APP2
 *   webp(img, opts)   lossless (VP8L): subtract-green and predictor transforms, runs copied from the left or above;
 *                     an ICC profile goes in a VP8X / ICCP container. At most 16384 × 16384 pixels
 *   tiff(img, opts)   baseline RGB(A), uncompressed or Deflate with the horizontal predictor, resolution, ICC profile
 *   pdf(img, opts)    one page of the print size; the image FlateDecode, an SMask for transparency, ICCBased colour
 *   icc(name)         ICC v2 display profiles: 'srgb' (IEC 61966-2.1 primaries and curve) and 'p3' (Display P3)
 *   toP3(img)         the same colours as Display P3 values (for a file tagged with the P3 profile)
 *   deflate(bytes, level), zlib(bytes, level), crc32(bytes), adler32(bytes)
 * img = {width, height, data: RGBA bytes, not premultiplied}. opts: dpi, profile ('srgb' | 'p3' | null; the pixels
 * must already be in that space), alpha (false drops the alpha channel), level (0-9), quality (JPEG). Each encoder
 * returns a Uint8Array. */
(function (BM) {
  'use strict';

  // ---------------------------------------------------------------- checksums
  const CRC = new Int32Array(256);
  for (let n = 0; n < 256; n++) { let c = n; for (let k = 0; k < 8; k++) c = c & 1 ? 0xedb88320 ^ (c >>> 1) : c >>> 1; CRC[n] = c; }
  function crc32(b, start, end, crc) {
    let c = crc === undefined ? -1 : ~crc;
    for (let i = start || 0, e = end === undefined ? b.length : end; i < e; i++) c = CRC[(c ^ b[i]) & 255] ^ (c >>> 8);
    return ~c >>> 0;
  }
  function adler32(b) {
    let a = 1, s = 0;
    for (let i = 0, n = b.length; i < n;) {
      const e = Math.min(n, i + 3800);
      for (; i < e; i++) { a += b[i]; s += a; }
      a %= 65521; s %= 65521;
    }
    return ((s << 16) | a) >>> 0;
  }

  // ---------------------------------------------------------------- byte and bit writers
  class Bytes {
    constructor(cap) { this.buf = new Uint8Array(cap || 4096); this.pos = 0; }
    need(k) { if (this.pos + k > this.buf.length) { const b = new Uint8Array(Math.max(this.buf.length * 2, this.pos + k + 4096)); b.set(this.buf.subarray(0, this.pos)); this.buf = b; } }
    u8(v) { this.need(1); this.buf[this.pos++] = v & 255; }
    u16be(v) { this.need(2); this.buf[this.pos++] = (v >>> 8) & 255; this.buf[this.pos++] = v & 255; }
    u32be(v) { this.need(4); for (let s = 24; s >= 0; s -= 8) this.buf[this.pos++] = (v >>> s) & 255; }
    u16le(v) { this.need(2); this.buf[this.pos++] = v & 255; this.buf[this.pos++] = (v >>> 8) & 255; }
    u32le(v) { this.need(4); for (let s = 0; s < 32; s += 8) this.buf[this.pos++] = (v >>> s) & 255; }
    bytes(a) { this.need(a.length); this.buf.set(a, this.pos); this.pos += a.length; }
    ascii(s) { this.need(s.length); for (let i = 0; i < s.length; i++) this.buf[this.pos++] = s.charCodeAt(i) & 255; }
    done() { return this.buf.slice(0, this.pos); }
  }
  /** Bits packed from the least significant bit (deflate, WebP lossless). */
  class BitsLE {
    constructor(cap) { this.buf = new Uint8Array(cap || 65536); this.pos = 0; this.acc = 0; this.n = 0; }
    put(v, n) {
      this.acc |= v << this.n; this.n += n;
      if (this.n >= 8) {
        if (this.pos + 4 > this.buf.length) { const b = new Uint8Array(this.buf.length * 2 + 4096); b.set(this.buf.subarray(0, this.pos)); this.buf = b; }
        do { this.buf[this.pos++] = this.acc & 255; this.acc >>>= 8; this.n -= 8; } while (this.n >= 8);
      }
    }
    align() { if (this.n > 0) this.put(0, 8 - this.n); }
    bytes(a) {
      this.align();
      if (this.pos + a.length > this.buf.length) { const b = new Uint8Array(Math.max(this.buf.length * 2, this.pos + a.length + 4096)); b.set(this.buf.subarray(0, this.pos)); this.buf = b; }
      this.buf.set(a, this.pos); this.pos += a.length;
    }
    done() { this.align(); return this.buf.slice(0, this.pos); }
  }

  // ---------------------------------------------------------------- Huffman codes
  /** Code lengths of at most `limit` bits for the frequencies: a Huffman code, then lengths over the limit moved down
   *  (JPEG Annex K.3). Symbols are ranked by frequency, then by number, so the result is the same everywhere. */
  function codeLengths(freq, limit) {
    const n = freq.length, len = new Uint8Array(n), syms = [];
    for (let i = 0; i < n; i++) if (freq[i] > 0) syms.push(i);
    if (!syms.length) return len;
    if (syms.length === 1) { len[syms[0]] = 1; return len; }
    syms.sort((a, b) => freq[a] - freq[b] || a - b);
    const m = syms.length, w = new Float64Array(2 * m), parent = new Int32Array(2 * m);
    for (let k = 0; k < m; k++) w[k] = freq[syms[k]];
    let leaf = 0, node = m;
    for (let next = m; next < 2 * m - 1; next++) {
      const pick = () => (leaf < m && (node >= next || w[leaf] <= w[node]) ? leaf++ : node++);
      const a = pick(), b = pick();
      w[next] = w[a] + w[b]; parent[a] = parent[b] = next;
    }
    const depth = new Int32Array(2 * m);
    let maxD = 0;
    for (let k = 2 * m - 3; k >= 0; k--) { depth[k] = depth[parent[k]] + 1; if (k < m && depth[k] > maxD) maxD = depth[k]; }
    const count = new Int32Array(Math.max(maxD, limit) + 2);
    for (let k = 0; k < m; k++) count[depth[k]]++;
    for (let i = maxD; i > limit; i--) {
      while (count[i] > 0) {
        let j = i - 2;
        while (count[j] === 0) j--;
        count[i] -= 2; count[i - 1]++; count[j + 1] += 2; count[j]--;
      }
    }
    // the most frequent symbols get the shortest codes
    let k = m - 1;
    for (let l = 1; l <= limit; l++) for (let c = count[l]; c > 0; c--) len[syms[k--]] = l;
    return len;
  }
  /** Canonical codes (most significant bit first) for code lengths (RFC 1951, 3.2.2). */
  function canonical(len) {
    let maxL = 0;
    for (let i = 0; i < len.length; i++) if (len[i] > maxL) maxL = len[i];
    const bl = new Int32Array(maxL + 1), next = new Int32Array(maxL + 1), codes = new Int32Array(len.length);
    for (let i = 0; i < len.length; i++) if (len[i]) bl[len[i]]++;
    let code = 0;
    for (let b = 1; b <= maxL; b++) { code = (code + bl[b - 1]) << 1; next[b] = code; }
    for (let i = 0; i < len.length; i++) if (len[i]) codes[i] = next[len[i]]++;
    return codes;
  }
  function reverse(c, l) { let r = 0; for (let i = 0; i < l; i++) { r = (r << 1) | (c & 1); c >>>= 1; } return r; }
  /** Codes for a writer that packs from the least significant bit: the canonical codes bit-reversed. */
  const lsbCodes = len => { const c = canonical(len); for (let i = 0; i < c.length; i++) if (len[i]) c[i] = reverse(c[i], len[i]); return c; };

  // ---------------------------------------------------------------- deflate (RFC 1951) and zlib (RFC 1950)
  const LBASE = [3, 4, 5, 6, 7, 8, 9, 10, 11, 13, 15, 17, 19, 23, 27, 31, 35, 43, 51, 59, 67, 83, 99, 115, 131, 163, 195, 227, 258];
  const LEXT = [0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 2, 2, 2, 2, 3, 3, 3, 3, 4, 4, 4, 4, 5, 5, 5, 5, 0];
  const DBASE = [1, 2, 3, 4, 5, 7, 9, 13, 17, 25, 33, 49, 65, 97, 129, 193, 257, 385, 513, 769, 1025, 1537, 2049, 3073, 4097, 6145, 8193, 12289, 16385, 24577];
  const DEXT = [0, 0, 0, 0, 1, 1, 2, 2, 3, 3, 4, 4, 5, 5, 6, 6, 7, 7, 8, 8, 9, 9, 10, 10, 11, 11, 12, 12, 13, 13];
  const LCODE = new Uint8Array(259), DCODE = new Uint8Array(32769);
  for (let c = 0; c < 29; c++) for (let l = LBASE[c], e = c < 28 ? LBASE[c + 1] : 259; l < e; l++) LCODE[l] = c;
  for (let c = 0; c < 30; c++) for (let d = DBASE[c], e = c < 29 ? DBASE[c + 1] : 32769; d < e; d++) DCODE[d] = c;
  const CL_ORDER = [16, 17, 18, 0, 8, 7, 9, 6, 10, 5, 11, 4, 12, 3, 13, 2, 14, 1, 15];
  // per level: hash-chain length, lazy matching, a match long enough to stop looking
  const LEVELS = [null, [4, 0, 8], [8, 0, 16], [16, 0, 32], [16, 1, 16], [32, 1, 32], [64, 1, 128], [128, 1, 128], [256, 1, 258], [1024, 1, 258]];

  /** Run-length coded code lengths: [[symbol, extra bits value], ...] with 16 (repeat the previous length 3-6 times),
   *  17 and 18 (3-10 and 11-138 zeros). For deflate "previous" is the last length, zeros included (zeroResets); for
   *  WebP it is the last non-zero length, 8 at first. */
  function rleLengths(lens, prevInit, zeroResets) {
    const out = [];
    let prev = prevInit;
    for (let i = 0; i < lens.length;) {
      const v = lens[i];
      let r = 1;
      while (i + r < lens.length && lens[i + r] === v) r++;
      if (v === 0) {
        let left = r;
        while (left >= 11) { const k = Math.min(138, left); out.push([18, k - 11]); left -= k; }
        if (left >= 3) { out.push([17, left - 3]); left = 0; }
        while (left-- > 0) out.push([0, 0]);
        if (zeroResets) prev = 0;
      } else {
        let left = r;
        if (v !== prev) { out.push([v, 0]); left--; prev = v; }
        while (left >= 3) { const k = Math.min(6, left); out.push([16, k - 3]); left -= k; }
        while (left-- > 0) out.push([v, 0]);
      }
      i += r;
    }
    return out;
  }
  const EXTRA_CL = { 16: 2, 17: 3, 18: 7 };

  /** Raw deflate stream of `data` (Uint8Array) at level 0 (stored) to 9. */
  function deflate(data, level) {
    level = level === undefined ? 6 : Math.max(0, Math.min(9, level | 0));
    const n = data.length, out = new BitsLE(Math.max(4096, (n >> 1) + 1024));
    const stored = (from, to, final) => {
      if (from === to) { out.put(final ? 1 : 0, 1); out.put(0, 2); out.align(); out.bytes([0, 0, 255, 255]); return; }
      for (let p = from; p < to;) {
        const k = Math.min(65535, to - p), last = final && p + k === to;
        out.put(last ? 1 : 0, 1); out.put(0, 2); out.align();
        out.bytes([k & 255, k >>> 8, ~k & 255, (~k >>> 8) & 255]);
        out.bytes(data.subarray(p, p + k));
        p += k;
      }
    };
    if (level === 0 || n === 0) { stored(0, n, true); return out.done(); }
    const [chainMax, lazy, nice] = LEVELS[level];
    const WSIZE = 32768, WMASK = WSIZE - 1, HSIZE = 1 << 15, HMASK = HSIZE - 1;
    const head = new Int32Array(HSIZE).fill(-1), prev = new Int32Array(WSIZE);
    const BLOCK = 1 << 15;
    const sLit = new Uint16Array(BLOCK + 2), sDist = new Uint16Array(BLOCK + 2);
    let nsym = 0, blockStart = 0;
    const hashAt = i => ((data[i] << 10) ^ (data[i + 1] << 5) ^ data[i + 2]) & HMASK;
    const insert = i => { const h = hashAt(i); prev[i & WMASK] = head[h]; head[h] = i; };
    const search = i => {
      const maxLen = Math.min(258, n - i);
      if (maxLen < 3) return 0;
      let p = head[hashAt(i)], best = 2, bestD = 0, left = chainMax;
      while (p >= 0 && i - p <= WSIZE && left-- > 0) {
        if (data[p + best] === data[i + best] && data[p] === data[i] && data[p + 1] === data[i + 1]) {
          let l = 2;
          while (l < maxLen && data[p + l] === data[i + l]) l++;
          if (l > best) { best = l; bestD = i - p; if (l >= nice) break; }
        }
        p = prev[p & WMASK];
      }
      return bestD ? best * 65536 + bestD : 0;
    };
    const flush = (end, final) => {
      const lf = new Float64Array(286), df = new Float64Array(30);
      for (let k = 0; k < nsym; k++) {
        if (sDist[k] === 0) lf[sLit[k]]++;
        else { lf[257 + LCODE[sLit[k]]]++; df[DCODE[sDist[k]]]++; }
      }
      lf[256]++;
      const two = f => { let c = 0; for (let k = 0; k < f.length; k++) if (f[k]) c++; for (let k = 0; c < 2 && k < f.length; k++) if (!f[k]) { f[k] = 1; c++; } };
      two(lf); two(df);
      const ll = codeLengths(lf, 15), dl = codeLengths(df, 15);
      let hlit = 286, hdist = 30;
      while (hlit > 257 && !ll[hlit - 1]) hlit--;
      while (hdist > 1 && !dl[hdist - 1]) hdist--;
      const all = new Uint8Array(hlit + hdist);
      all.set(ll.subarray(0, hlit)); all.set(dl.subarray(0, hdist), hlit);
      const rle = rleLengths(all, -1, true), cf = new Float64Array(19);
      for (const [s] of rle) cf[s]++;
      two(cf);
      const cl = codeLengths(cf, 7);
      let hclen = 19;
      while (hclen > 4 && !cl[CL_ORDER[hclen - 1]]) hclen--;
      // size in bits: dynamic block against stored blocks
      let bits = 3 + 14 + 3 * hclen;
      for (const [s] of rle) bits += cl[s] + (EXTRA_CL[s] || 0);
      for (let k = 0; k < nsym; k++) {
        if (sDist[k] === 0) bits += ll[sLit[k]];
        else { const c = LCODE[sLit[k]], d = DCODE[sDist[k]]; bits += ll[257 + c] + LEXT[c] + dl[d] + DEXT[d]; }
      }
      bits += ll[256];
      const raw = end - blockStart, storedBits = (raw + 5 * Math.max(1, Math.ceil(raw / 65535))) * 8 + 7;
      if (storedBits < bits) stored(blockStart, end, final);
      else {
        const lc = lsbCodes(ll), dc = lsbCodes(dl), cc = lsbCodes(cl);
        out.put(final ? 1 : 0, 1); out.put(2, 2);
        out.put(hlit - 257, 5); out.put(hdist - 1, 5); out.put(hclen - 4, 4);
        for (let k = 0; k < hclen; k++) out.put(cl[CL_ORDER[k]], 3);
        for (const [s, x] of rle) { out.put(cc[s], cl[s]); if (EXTRA_CL[s]) out.put(x, EXTRA_CL[s]); }
        for (let k = 0; k < nsym; k++) {
          if (sDist[k] === 0) out.put(lc[sLit[k]], ll[sLit[k]]);
          else {
            const len = sLit[k], c = LCODE[len], d = sDist[k], dcd = DCODE[d];
            out.put(lc[257 + c], ll[257 + c]); if (LEXT[c]) out.put(len - LBASE[c], LEXT[c]);
            out.put(dc[dcd], dl[dcd]); if (DEXT[dcd]) out.put(d - DBASE[dcd], DEXT[dcd]);
          }
        }
        out.put(lc[256], ll[256]);
      }
      nsym = 0; blockStart = end;
    };
    let i = 0, m = search(0);
    while (i < n) {
      if (i + 2 < n) insert(i);
      const len = (m / 65536) | 0;
      if (len >= 3) {
        if (lazy && len < nice && i + 3 < n) {
          const m2 = search(i + 1);
          if (((m2 / 65536) | 0) > len) {
            sLit[nsym] = data[i]; sDist[nsym++] = 0; i++; m = m2;
            if (nsym >= BLOCK) flush(i, false);
            continue;
          }
        }
        sLit[nsym] = len; sDist[nsym++] = m & 0xffff;
        const end = i + len;
        for (let k = i + 1; k < end && k + 2 < n; k++) insert(k);
        i = end;
      } else { sLit[nsym] = data[i]; sDist[nsym++] = 0; i++; }
      if (nsym >= BLOCK) flush(i, false);
      m = i < n ? search(i) : 0;
    }
    if (nsym || blockStart < n) flush(n, true);
    else { out.put(1, 1); out.put(1, 2); out.put(0, 7); }   // an empty final block (fixed codes, end of block only)
    return out.done();
  }
  function zlib(data, level) {
    level = level === undefined ? 6 : level;
    const d = deflate(data, level), out = new Uint8Array(d.length + 6);
    out[0] = 0x78; out[1] = level <= 1 ? 0x01 : level < 6 ? 0x5e : level === 6 ? 0x9c : 0xda;
    out.set(d, 2);
    const a = adler32(data);
    out[d.length + 2] = a >>> 24; out[d.length + 3] = (a >>> 16) & 255; out[d.length + 4] = (a >>> 8) & 255; out[d.length + 5] = a & 255;
    return out;
  }

  // ---------------------------------------------------------------- colour: ICC profiles, sRGB -> Display P3
  const PRIMARIES = {
    srgb: { name: 'sRGB IEC61966-2.1', r: [0.64, 0.33], g: [0.30, 0.60], b: [0.15, 0.06], w: [0.3127, 0.3290] },
    p3: { name: 'Display P3', r: [0.680, 0.320], g: [0.265, 0.690], b: [0.150, 0.060], w: [0.3127, 0.3290] }
  };
  const xyz = ([x, y]) => [x / y, 1, (1 - x - y) / y];
  function inv3(m) {
    const [a, b, c, d, e, f, g, h, i] = m, A = e * i - f * h, B = f * g - d * i, C = d * h - e * g, det = a * A + b * B + c * C;
    return [A / det, (c * h - b * i) / det, (b * f - c * e) / det, B / det, (a * i - c * g) / det, (c * d - a * f) / det, C / det, (b * g - a * h) / det, (a * e - b * d) / det];
  }
  const mul3 = (a, b) => { const o = new Array(9); for (let r = 0; r < 3; r++) for (let c = 0; c < 3; c++) o[3 * r + c] = a[3 * r] * b[c] + a[3 * r + 1] * b[3 + c] + a[3 * r + 2] * b[6 + c]; return o; };
  const apply3 = (m, v) => [m[0] * v[0] + m[1] * v[1] + m[2] * v[2], m[3] * v[0] + m[4] * v[1] + m[5] * v[2], m[6] * v[0] + m[7] * v[1] + m[8] * v[2]];
  /** Linear RGB -> XYZ (white of the space) for chromaticities. */
  function rgbToXyz(p) {
    const R = xyz(p.r), G = xyz(p.g), B = xyz(p.b), W = xyz(p.w);
    const M = [R[0], G[0], B[0], R[1], G[1], B[1], R[2], G[2], B[2]], S = apply3(inv3(M), W);
    return [M[0] * S[0], M[1] * S[1], M[2] * S[2], M[3] * S[0], M[4] * S[1], M[5] * S[2], M[6] * S[0], M[7] * S[1], M[8] * S[2]];
  }
  const D50 = [0.9642, 1.0, 0.8249];   // the ICC profile connection space white
  const BRADFORD = [0.8951, 0.2664, -0.1614, -0.7502, 1.7135, 0.0367, 0.0389, -0.0685, 1.0296];
  function adapt(from, to) {   // Bradford chromatic adaptation matrix (XYZ -> XYZ)
    const a = apply3(BRADFORD, from), b = apply3(BRADFORD, to);
    return mul3(inv3(BRADFORD), mul3([b[0] / a[0], 0, 0, 0, b[1] / a[1], 0, 0, 0, b[2] / a[2]], BRADFORD));
  }
  const s15 = v => Math.round(v * 65536) | 0;
  const ICC_CACHE = {};
  /** An ICC v2.1 display profile (matrix and tone curves) for 'srgb' or 'p3'. */
  function icc(name) {
    const p = PRIMARIES[name];
    if (!p) return null;
    if (ICC_CACHE[name]) return ICC_CACHE[name];
    const M = rgbToXyz(p), chad = adapt(xyz(p.w), D50), C = mul3(chad, M);   // columns: the D50 XYZ of R, G, B
    const tags = [], ASCII = s => { const b = new Bytes(); b.ascii(s); return b.done(); };
    const desc = (() => {
      const b = new Bytes(), s = p.name;
      b.ascii('desc'); b.u32be(0); b.u32be(s.length + 1); b.ascii(s); b.u8(0);
      b.u32be(0); b.u32be(0); b.u16be(0); b.u8(0); for (let k = 0; k < 67; k++) b.u8(0);
      return b.done();
    })();
    const text = s => { const b = new Bytes(); b.ascii('text'); b.u32be(0); b.ascii(s); b.u8(0); return b.done(); };
    const XYZ = v => { const b = new Bytes(); b.ascii('XYZ '); b.u32be(0); for (const x of v) b.u32be(s15(x)); return b.done(); };
    const curve = (() => {   // the sRGB transfer function, 1024 points
      const b = new Bytes(), N = 1024;
      b.ascii('curv'); b.u32be(0); b.u32be(N);
      for (let k = 0; k < N; k++) { const v = k / (N - 1), l = v <= 0.04045 ? v / 12.92 : Math.pow((v + 0.055) / 1.055, 2.4); b.u16be(Math.round(l * 65535)); }
      return b.done();
    })();
    const sf32 = m => { const b = new Bytes(); b.ascii('sf32'); b.u32be(0); for (const x of m) b.u32be(s15(x)); return b.done(); };
    tags.push(['desc', desc], ['cprt', text('No copyright, use freely')], ['wtpt', XYZ(D50)],
      ['rXYZ', XYZ([C[0], C[3], C[6]])], ['gXYZ', XYZ([C[1], C[4], C[7]])], ['bXYZ', XYZ([C[2], C[5], C[8]])],
      ['rTRC', curve], ['gTRC', curve], ['bTRC', curve], ['chad', sf32(chad)]);
    // layout: header (128), tag table, data (4-byte aligned; the three curves share one copy)
    const table = 4 + 12 * tags.length, at = new Map();
    let off = 128 + table;
    const blobs = [];
    for (const [, data] of tags) {
      if (at.has(data)) continue;
      at.set(data, off); blobs.push([off, data]);
      off += (data.length + 3) & ~3;
    }
    const size = off, b = new Bytes(size);
    b.u32be(size); b.u32be(0); b.u32be(0x02100000); b.ascii('mntrRGB XYZ ');
    for (const v of [2026, 11, 22, 0, 0, 0]) b.u16be(v);
    b.ascii('acsp'); b.u32be(0); b.u32be(0); b.u32be(0); b.u32be(0); b.u32be(0); b.u32be(0); b.u32be(0);
    for (const v of D50) b.u32be(s15(v));
    b.ascii('BMol'); for (let k = 0; k < 44; k++) b.u8(0);
    b.u32be(tags.length);
    for (const [sig, data] of tags) { b.ascii(sig); b.u32be(at.get(data)); b.u32be(data.length); }
    for (const [o, data] of blobs) { while (b.pos < o) b.u8(0); b.bytes(data); }
    while (b.pos < size) b.u8(0);
    return (ICC_CACHE[name] = b.done());
  }
  /** sRGB pixels as Display P3 values of the same colours (exact matrix in linear light; alpha kept). */
  function toP3(img) {
    const toLin = new Float64Array(256);
    for (let k = 0; k < 256; k++) { const v = k / 255; toLin[k] = v <= 0.04045 ? v / 12.92 : Math.pow((v + 0.055) / 1.055, 2.4); }
    const N = 4096, enc = new Uint8Array(N + 1);
    for (let k = 0; k <= N; k++) { const l = k / N, v = l <= 0.0031308 ? 12.92 * l : 1.055 * Math.pow(l, 1 / 2.4) - 0.055; enc[k] = Math.max(0, Math.min(255, Math.round(v * 255))); }
    const M = mul3(inv3(rgbToXyz(PRIMARIES.p3)), rgbToXyz(PRIMARIES.srgb)), src = img.data, out = new Uint8Array(src.length);
    const q = l => enc[Math.max(0, Math.min(N, Math.round(l * N)))];
    for (let i = 0; i < src.length; i += 4) {
      const r = toLin[src[i]], g = toLin[src[i + 1]], b = toLin[src[i + 2]];
      out[i] = q(M[0] * r + M[1] * g + M[2] * b); out[i + 1] = q(M[3] * r + M[4] * g + M[5] * b); out[i + 2] = q(M[6] * r + M[7] * g + M[8] * b); out[i + 3] = src[i + 3];
    }
    return { width: img.width, height: img.height, data: out };
  }

  const hasAlpha = d => { for (let i = 3; i < d.length; i += 4) if (d[i] !== 255) return true; return false; };
  const VERSION = () => (BM.version || '');

  // ---------------------------------------------------------------- PNG
  function png(img, o) {
    o = o || {};
    const W = img.width, H = img.height, src = img.data, level = o.level === undefined ? 6 : o.level;
    const alpha = o.alpha !== false && hasAlpha(src), ch = alpha ? 4 : 3, stride = W * ch;
    const raw = new Uint8Array((stride + 1) * H);
    let prev = new Uint8Array(stride), row = new Uint8Array(stride);
    const cand = [0, 1, 2, 3, 4].map(() => new Uint8Array(stride));
    for (let y = 0; y < H; y++) {
      for (let x = 0, s = y * W * 4, t = 0; x < W; x++, s += 4) { row[t++] = src[s]; row[t++] = src[s + 1]; row[t++] = src[s + 2]; if (alpha) row[t++] = src[s + 3]; }
      let best = 0;
      if (level > 0) {
        let bestSum = Infinity;
        for (let f = 0; f < 5; f++) {
          const c = cand[f];
          let sum = 0;
          for (let i = 0; i < stride; i++) {
            const a = i >= ch ? row[i - ch] : 0, b = prev[i], cc = i >= ch ? prev[i - ch] : 0;
            let v;
            if (f === 0) v = row[i];
            else if (f === 1) v = row[i] - a;
            else if (f === 2) v = row[i] - b;
            else if (f === 3) v = row[i] - ((a + b) >> 1);
            else { const p = a + b - cc, pa = Math.abs(p - a), pb = Math.abs(p - b), pc = Math.abs(p - cc); v = row[i] - (pa <= pb && pa <= pc ? a : pb <= pc ? b : cc); }
            v &= 255; c[i] = v; sum += v < 128 ? v : 256 - v;
          }
          if (sum < bestSum) { bestSum = sum; best = f; }
        }
      } else cand[0].set(row);
      raw[y * (stride + 1)] = best;
      raw.set(cand[best], y * (stride + 1) + 1);
      const t = prev; prev = row; row = t;
    }
    const b = new Bytes(raw.length / 2 + 1024);
    const chunk = (type, data) => {
      b.u32be(data.length); const s = b.pos; b.ascii(type); b.bytes(data);
      b.u32be(crc32(b.buf, s, b.pos));
    };
    b.bytes([137, 80, 78, 71, 13, 10, 26, 10]);
    const ihdr = new Bytes(13); ihdr.u32be(W); ihdr.u32be(H); ihdr.bytes([8, alpha ? 6 : 2, 0, 0, 0]); chunk('IHDR', ihdr.done());
    if (o.profile === 'srgb') {
      chunk('sRGB', new Uint8Array([0]));
      const g = new Bytes(); g.u32be(45455); chunk('gAMA', g.done());
      const c = new Bytes(); for (const v of [31270, 32900, 64000, 33000, 30000, 60000, 15000, 6000]) c.u32be(v); chunk('cHRM', c.done());
    } else if (o.profile) {
      const p = icc(o.profile), c = new Bytes(); c.ascii(PRIMARIES[o.profile].name.slice(0, 79)); c.u8(0); c.u8(0); c.bytes(zlib(p, 9)); chunk('iCCP', c.done());
    }
    if (o.dpi > 0) { const ppm = Math.round(o.dpi / 0.0254), p = new Bytes(9); p.u32be(ppm); p.u32be(ppm); p.u8(1); chunk('pHYs', p.done()); }
    const sw = new Bytes(); sw.ascii('Software'); sw.u8(0); sw.ascii(('BlueMolecule ' + VERSION()).trim()); chunk('tEXt', sw.done());
    const z = zlib(raw, level);
    for (let p = 0; p < z.length || p === 0; p += 1 << 18) chunk('IDAT', z.subarray(p, Math.min(z.length, p + (1 << 18))));
    chunk('IEND', new Uint8Array(0));
    return b.done();
  }

  // ---------------------------------------------------------------- JPEG (baseline, 4:4:4)
  const ZIGZAG = new Uint8Array([0, 1, 8, 16, 9, 2, 3, 10, 17, 24, 32, 25, 18, 11, 4, 5, 12, 19, 26, 33, 40, 48, 41, 34, 27, 20, 13, 6, 7, 14, 21, 28,
    35, 42, 49, 56, 57, 50, 43, 36, 29, 22, 15, 23, 30, 37, 44, 51, 58, 59, 52, 45, 38, 31, 39, 46, 53, 60, 61, 54, 47, 55, 62, 63]);   // zigzag position -> natural index
  const Q_LUMA = [16, 11, 10, 16, 24, 40, 51, 61, 12, 12, 14, 19, 26, 58, 60, 55, 14, 13, 16, 24, 40, 57, 69, 56, 14, 17, 22, 29, 51, 87, 80, 62,
    18, 22, 37, 56, 68, 109, 103, 77, 24, 35, 55, 64, 81, 104, 113, 92, 49, 64, 78, 87, 103, 121, 120, 101, 72, 92, 95, 98, 112, 100, 103, 99];
  const Q_CHROMA = [17, 18, 24, 47, 99, 99, 99, 99, 18, 21, 26, 66, 99, 99, 99, 99, 24, 26, 56, 99, 99, 99, 99, 99, 47, 66, 99, 99, 99, 99, 99, 99,
    99, 99, 99, 99, 99, 99, 99, 99, 99, 99, 99, 99, 99, 99, 99, 99, 99, 99, 99, 99, 99, 99, 99, 99, 99, 99, 99, 99, 99, 99, 99, 99];   // natural order (ITU T.81 Annex K)
  const AASF = [1.0, 1.387039845, 1.306562965, 1.175875602, 1.0, 0.785694958, 0.541196100, 0.275899379];
  function fdct(d) {   // AAN floating-point forward DCT, in place (IJG jfdctflt)
    for (let p = 0; p < 2; p++) {
      const step = p === 0 ? 1 : 8, line = p === 0 ? 8 : 1;
      for (let k = 0; k < 8; k++) {
        const o = k * line, i0 = o, i1 = o + step, i2 = o + 2 * step, i3 = o + 3 * step, i4 = o + 4 * step, i5 = o + 5 * step, i6 = o + 6 * step, i7 = o + 7 * step;
        const t0 = d[i0] + d[i7], t7 = d[i0] - d[i7], t1 = d[i1] + d[i6], t6 = d[i1] - d[i6], t2 = d[i2] + d[i5], t5 = d[i2] - d[i5], t3 = d[i3] + d[i4], t4 = d[i3] - d[i4];
        let t10 = t0 + t3, t13 = t0 - t3, t11 = t1 + t2, t12 = t1 - t2;
        d[i0] = t10 + t11; d[i4] = t10 - t11;
        const z1 = (t12 + t13) * 0.707106781;
        d[i2] = t13 + z1; d[i6] = t13 - z1;
        t10 = t4 + t5; t11 = t5 + t6; t12 = t6 + t7;
        const z5 = (t10 - t12) * 0.382683433, z2 = 0.541196100 * t10 + z5, z4 = 1.306562965 * t12 + z5, z3 = t11 * 0.707106781;
        const z11 = t7 + z3, z13 = t7 - z3;
        d[i5] = z13 + z2; d[i3] = z13 - z2; d[i1] = z11 + z4; d[i7] = z11 - z4;
      }
    }
  }
  const nbits = v => { v = v < 0 ? -v : v; let n = 0; while (v) { n++; v >>= 1; } return n; };
  function jpeg(img, o) {
    o = o || {};
    const W = img.width, H = img.height, src = img.data;
    if (W > 65535 || H > 65535) throw new Error('JPEG images are at most 65,535 pixels wide and high.');
    const quality = Math.max(1, Math.min(100, Math.round(o.quality === undefined ? 92 : o.quality)));
    const scale = quality < 50 ? 5000 / quality : 200 - 2 * quality;
    const qt = base => base.map(v => Math.max(1, Math.min(255, Math.floor((v * scale + 50) / 100))));
    const qY = qt(Q_LUMA), qC = qt(Q_CHROMA);
    const ft = q => Float64Array.from(q, (v, i) => 1 / (v * AASF[i >> 3] * AASF[i & 7] * 8));
    const fY = ft(qY), fC = ft(qC);
    // a transparent image goes on white (JPEG has no alpha)
    const bx = Math.ceil(W / 8), by = Math.ceil(H / 8);
    const blk = [new Float64Array(64), new Float64Array(64), new Float64Array(64)], coef = new Int32Array(64);
    const blockCoefs = (bxi, byi, comp) => {   // quantised coefficients of one 8x8 block, zigzag order, into coef
      const d = blk[comp];
      for (let yy = 0; yy < 8; yy++) {
        const y = Math.min(H - 1, byi * 8 + yy);
        for (let xx = 0; xx < 8; xx++) {
          const x = Math.min(W - 1, bxi * 8 + xx), s = (y * W + x) * 4, a = src[s + 3] / 255;
          const r = src[s] * a + 255 * (1 - a), g = src[s + 1] * a + 255 * (1 - a), b = src[s + 2] * a + 255 * (1 - a);
          d[yy * 8 + xx] = comp === 0 ? 0.299 * r + 0.587 * g + 0.114 * b - 128 : comp === 1 ? -0.168735892 * r - 0.331264108 * g + 0.5 * b : 0.5 * r - 0.418687589 * g - 0.081312411 * b;
        }
      }
      fdct(d);
      const f = comp === 0 ? fY : fC;
      for (let k = 0; k < 64; k++) { const i = ZIGZAG[k], v = d[i] * f[i]; coef[k] = v > 0 ? Math.floor(v + 0.5) : -Math.floor(-v + 0.5); }
    };
    // two passes: count the symbols for Huffman tables made for this image, then write
    const freq = [new Float64Array(257), new Float64Array(257), new Float64Array(257), new Float64Array(257)];   // DC Y, AC Y, DC C, AC C
    let codes = null, lens = null, bw = null;
    const emit = (t, sym) => { if (!codes) freq[t][sym]++; else bw.bits(codes[t][sym], lens[t][sym]); };
    const value = (v, n) => { if (codes && n) bw.bits(v < 0 ? v + (1 << n) - 1 : v, n); };
    const scan = () => {
      const pred = [0, 0, 0];
      for (let byi = 0; byi < by; byi++) for (let bxi = 0; bxi < bx; bxi++) for (let c = 0; c < 3; c++) {
        blockCoefs(bxi, byi, c);
        const tDC = c === 0 ? 0 : 2, tAC = tDC + 1, diff = coef[0] - pred[c];
        pred[c] = coef[0];
        const s = nbits(diff);
        emit(tDC, s); value(diff, s);
        let run = 0;
        for (let k = 1; k < 64; k++) {
          const v = coef[k];
          if (v === 0) { run++; continue; }
          while (run > 15) { emit(tAC, 0xf0); run -= 16; }
          const n = nbits(v);
          emit(tAC, (run << 4) | n); value(v, n);
          run = 0;
        }
        if (run > 0) emit(tAC, 0);
      }
    };
    scan();
    // Huffman tables: at most 16 bits, no code of all ones (a reserved symbol 256 takes it, then is removed)
    lens = []; codes = [];
    const dht = [];
    for (let t = 0; t < 4; t++) {
      const f = freq[t];
      f[256] = 1;
      let used = 0;
      for (let s = 0; s < 256; s++) if (f[s]) used++;
      if (!used) f[0] = 1;
      const l = codeLengths(f, 16);
      // the reserved symbol must hold the longest length (then its code is the one of all ones): swap it with a
      // symbol of the same frequency that got it
      let maxL = 0;
      for (let s = 0; s <= 256; s++) if (l[s] > maxL) maxL = l[s];
      if (l[256] < maxL) for (let s = 255; s >= 0; s--) if (l[s] === maxL && f[s] === 1) { l[s] = l[256]; l[256] = maxL; break; }
      l[256] = 0;
      const c = canonical(l);
      lens.push(l); codes.push(c);
      const counts = new Array(16).fill(0), vals = [];
      for (let L = 1; L <= 16; L++) for (let s = 0; s < 256; s++) if (l[s] === L) { counts[L - 1]++; vals.push(s); }
      dht.push([t === 0 ? 0x00 : t === 1 ? 0x10 : t === 2 ? 0x01 : 0x11, counts, vals]);
    }
    const out = new Bytes(W * H + 4096);
    const marker = (m, body) => { out.u8(0xff); out.u8(m); out.u16be(body.length + 2); out.bytes(body); };
    out.u8(0xff); out.u8(0xd8);
    const dpi = o.dpi > 0 ? Math.min(65535, Math.round(o.dpi)) : 0, j = new Bytes();
    j.ascii('JFIF'); j.u8(0); j.u8(1); j.u8(1); j.u8(dpi ? 1 : 0); j.u16be(dpi || 1); j.u16be(dpi || 1); j.u8(0); j.u8(0);
    marker(0xe0, j.done());
    if (o.profile) {
      const p = icc(o.profile), max = 65519, n = Math.ceil(p.length / max);
      for (let k = 0; k < n; k++) { const a = new Bytes(); a.ascii('ICC_PROFILE'); a.u8(0); a.u8(k + 1); a.u8(n); a.bytes(p.subarray(k * max, (k + 1) * max)); marker(0xe2, a.done()); }
    }
    const q = new Bytes();
    q.u8(0); for (let k = 0; k < 64; k++) q.u8(qY[ZIGZAG[k]]);
    q.u8(1); for (let k = 0; k < 64; k++) q.u8(qC[ZIGZAG[k]]);
    marker(0xdb, q.done());
    const sof = new Bytes(); sof.u8(8); sof.u16be(H); sof.u16be(W); sof.u8(3);
    sof.bytes([1, 0x11, 0, 2, 0x11, 1, 3, 0x11, 1]);
    marker(0xc0, sof.done());
    const h = new Bytes();
    for (const [id, counts, vals] of dht) { h.u8(id); for (const c of counts) h.u8(c); h.bytes(vals); }
    marker(0xc4, h.done());
    marker(0xda, new Uint8Array([3, 1, 0x00, 2, 0x11, 3, 0x11, 0, 63, 0]));
    // entropy-coded data, most significant bit first, 0xFF followed by 0x00
    let acc = 0, nb = 0;
    bw = {
      bits(v, n) {
        acc = (acc << n) | (v & ((1 << n) - 1)); nb += n;
        while (nb >= 8) { const byte = (acc >>> (nb - 8)) & 255; out.u8(byte); if (byte === 255) out.u8(0); nb -= 8; acc &= (1 << nb) - 1; }
      }
    };
    scan();
    if (nb > 0) bw.bits((1 << (8 - nb)) - 1, 8 - nb);   // pad with ones
    out.u8(0xff); out.u8(0xd9);
    return out.done();
  }

  // ---------------------------------------------------------------- WebP lossless (VP8L)
  const CL_ORDER_WEBP = [17, 18, 0, 1, 2, 3, 4, 5, 16, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15];
  function prefix(v) {   // value >= 1 -> [prefix code, extra bits, extra value]
    const d = v - 1;
    if (d < 4) return [d, 0, 0];
    const hb = 31 - Math.clz32(d), second = (d >>> (hb - 1)) & 1, extra = hb - 1;
    return [2 * hb + second, extra, d & ((1 << extra) - 1)];
  }
  /** Write one prefix code (VP8L): a simple code for 1 or 2 symbols below 256, otherwise code lengths, which
   *  themselves are coded. Returns the codes for writing the data (bit-reversed) and their lengths. */
  function writeWebpCode(bw, freq) {
    const n = freq.length, used = [];
    for (let s = 0; s < n; s++) if (freq[s] > 0) used.push(s);
    if (used.length <= 2 && used.every(s => s < 256)) {
      const len = new Uint8Array(n), codes = new Int32Array(n);
      if (!used.length) used.push(0);
      bw.put(1, 1); bw.put(used.length - 1, 1);
      if (used[0] < 2) { bw.put(0, 1); bw.put(used[0], 1); } else { bw.put(1, 1); bw.put(used[0], 8); }
      if (used.length === 2) { bw.put(used[1], 8); len[used[0]] = len[used[1]] = 1; codes[used[1]] = 1; }
      return { len, codes };
    }
    const len = codeLengths(freq, 15);   // one symbol: its length is written, and it is then coded with no bits
    const rle = rleLengths(len, 8, false), cf = new Float64Array(19);
    for (const [s] of rle) cf[s]++;
    let nz = 0;
    for (let s = 0; s < 19; s++) if (cf[s]) nz++;
    for (let s = 0; nz < 2 && s < 19; s++) if (!cf[s]) { cf[s] = 1; nz++; }
    const cl = codeLengths(cf, 7), cc = lsbCodes(cl);
    let ncl = 19;
    while (ncl > 4 && !cl[CL_ORDER_WEBP[ncl - 1]]) ncl--;
    bw.put(0, 1); bw.put(ncl - 4, 4);
    for (let k = 0; k < ncl; k++) bw.put(cl[CL_ORDER_WEBP[k]], 3);
    bw.put(0, 1);   // every symbol has a length
    for (const [s, x] of rle) { bw.put(cc[s], cl[s]); if (EXTRA_CL[s]) bw.put(x, EXTRA_CL[s]); }
    if (used.length === 1) return { len: new Uint8Array(n), codes: new Int32Array(n) };
    return { len, codes: lsbCodes(len) };
  }
  /** An entropy-coded image: pixels (Uint32 ARGB) as literals or runs copied from the pixel on the left or above. */
  function writeWebpImage(bw, px, width, main) {
    bw.put(0, 1);              // no colour cache
    if (main) bw.put(0, 1);    // one group of prefix codes
    const n = px.length, tok = new Int32Array(n), tl = new Int32Array(n);   // tl > 0: run of tl pixels, tok = distance code
    let nt = 0;
    for (let i = 0; i < n;) {
      let left = 0, up = 0;
      if (i >= 1) while (i + left < n && left < 4096 && px[i + left] === px[i + left - 1]) left++;
      if (i >= width) while (i + up < n && up < 4096 && px[i + up] === px[i + up - width]) up++;
      if (Math.max(left, up) >= 3) {
        const useUp = up > left;
        tok[nt] = useUp ? 1 : 2; tl[nt++] = useUp ? up : left; i += useUp ? up : left;
      } else { tok[nt] = px[i]; tl[nt++] = 0; i++; }
    }
    const fG = new Float64Array(280), fR = new Float64Array(256), fB = new Float64Array(256), fA = new Float64Array(256), fD = new Float64Array(40);
    for (let k = 0; k < nt; k++) {
      if (tl[k]) { fG[256 + prefix(tl[k])[0]]++; fD[prefix(tok[k])[0]]++; }
      else { const p = tok[k]; fG[(p >>> 8) & 255]++; fR[(p >>> 16) & 255]++; fB[p & 255]++; fA[p >>> 24]++; }
    }
    const G = writeWebpCode(bw, fG), Rc = writeWebpCode(bw, fR), Bc = writeWebpCode(bw, fB), A = writeWebpCode(bw, fA), D = writeWebpCode(bw, fD);
    for (let k = 0; k < nt; k++) {
      if (tl[k]) {
        const [lc, le, lx] = prefix(tl[k]), [dc, de, dx] = prefix(tok[k]);
        bw.put(G.codes[256 + lc], G.len[256 + lc]); if (le) bw.put(lx, le);
        bw.put(D.codes[dc], D.len[dc]); if (de) bw.put(dx, de);
      } else {
        const p = tok[k], g = (p >>> 8) & 255, r = (p >>> 16) & 255, b = p & 255, a = p >>> 24;
        bw.put(G.codes[g], G.len[g]); bw.put(Rc.codes[r], Rc.len[r]); bw.put(Bc.codes[b], Bc.len[b]); bw.put(A.codes[a], A.len[a]);
      }
    }
  }
  function webp(img, o) {
    o = o || {};
    const W = img.width, H = img.height, src = img.data;
    if (W > 16384 || H > 16384) throw new Error('WebP images are at most 16,384 pixels wide and high.');
    const alpha = o.alpha !== false && hasAlpha(src), n = W * H, px = new Uint32Array(n);
    for (let i = 0, s = 0; i < n; i++, s += 4) {   // ARGB, then subtract green
      const r = src[s], g = src[s + 1], b = src[s + 2], a = alpha ? src[s + 3] : 255;
      px[i] = ((a << 24) | (((r - g) & 255) << 16) | (g << 8) | ((b - g) & 255)) >>> 0;
    }
    // predictor transform: one mode per 32 x 32 block, the one with the smallest residuals
    const bits = 5, bs = 1 << bits, tw = Math.ceil(W / bs), th = Math.ceil(H / bs), modes = new Uint32Array(tw * th), res = new Uint32Array(n);
    const ch = (p, s) => (p >>> s) & 255;
    const predict = (mode, L, T, TL) => {
      if (mode === 1) return L;
      if (mode === 2) return T;
      let out = 0;
      if (mode === 11) {
        let pl = 0, pt = 0;
        for (let s = 0; s < 32; s += 8) { const e = ch(L, s) + ch(T, s) - ch(TL, s); pl += Math.abs(e - ch(L, s)); pt += Math.abs(e - ch(T, s)); }
        return pl < pt ? L : T;
      }
      for (let s = 0; s < 32; s += 8) out |= Math.max(0, Math.min(255, ch(L, s) + ch(T, s) - ch(TL, s))) << s;   // mode 12
      return out >>> 0;
    };
    const sub = (p, q) => {   // channel by channel, modulo 256
      const a = ((p >>> 24) - (q >>> 24)) & 255, r = (ch(p, 16) - ch(q, 16)) & 255, g = (ch(p, 8) - ch(q, 8)) & 255, b = ((p & 255) - (q & 255)) & 255;
      return ((a << 24) | (r << 16) | (g << 8) | b) >>> 0;
    };
    const cost = r => { let c = 0; for (let s = 0; s < 32; s += 8) { const v = ch(r, s); c += v < 128 ? v : 256 - v; } return c; };
    const residual = (x, y, i, mode) => {
      if (i === 0) return sub(px[0], 0xff000000);
      if (y === 0) return sub(px[i], px[i - 1]);
      if (x === 0) return sub(px[i], px[i - W]);
      return sub(px[i], predict(mode, px[i - 1], px[i - W], px[i - W - 1]));
    };
    for (let ty = 0; ty < th; ty++) for (let tx = 0; tx < tw; tx++) {
      let best = 1, bestCost = Infinity;
      for (const mode of [1, 2, 11, 12]) {
        let c = 0;
        for (let y = ty * bs; y < Math.min(H, ty * bs + bs) && c < bestCost; y++) for (let x = tx * bs; x < Math.min(W, tx * bs + bs); x++) c += cost(residual(x, y, y * W + x, mode));
        if (c < bestCost) { bestCost = c; best = mode; }
      }
      modes[ty * tw + tx] = best << 8;
      for (let y = ty * bs; y < Math.min(H, ty * bs + bs); y++) for (let x = tx * bs; x < Math.min(W, tx * bs + bs); x++) res[y * W + x] = residual(x, y, y * W + x, best);
    }
    const bw = new BitsLE(n + 1024);
    bw.put(0x2f, 8); bw.put(W - 1, 14); bw.put(H - 1, 14); bw.put(alpha ? 1 : 0, 1); bw.put(0, 3);
    bw.put(1, 1); bw.put(2, 2);                                       // subtract green
    bw.put(1, 1); bw.put(0, 2); bw.put(bits - 2, 3); writeWebpImage(bw, modes, tw, false);   // predictor
    bw.put(0, 1);                                                     // no more transforms
    writeWebpImage(bw, res, W, true);
    const vp8l = bw.done(), out = new Bytes(vp8l.length + 4096);
    const chunk = (type, data) => { out.ascii(type); out.u32le(data.length); out.bytes(data); if (data.length & 1) out.u8(0); };
    out.ascii('RIFF'); out.u32le(0); out.ascii('WEBP');
    if (o.profile) {
      const x = new Bytes(10);
      x.u8(0x20 | (alpha ? 0x10 : 0)); x.u8(0); x.u8(0); x.u8(0);
      for (const v of [W - 1, H - 1]) { x.u8(v & 255); x.u8((v >>> 8) & 255); x.u8((v >>> 16) & 255); }
      chunk('VP8X', x.done());
      chunk('ICCP', icc(o.profile));
    }
    chunk('VP8L', vp8l);
    const r = out.done(), size = r.length - 8;
    r[4] = size & 255; r[5] = (size >>> 8) & 255; r[6] = (size >>> 16) & 255; r[7] = (size >>> 24) & 255;
    return r;
  }

  // ---------------------------------------------------------------- TIFF
  function tiff(img, o) {
    o = o || {};
    const W = img.width, H = img.height, src = img.data, alpha = o.alpha !== false && hasAlpha(src), spp = alpha ? 4 : 3;
    const deflated = (o.level === undefined ? 6 : o.level) > 0, stride = W * spp;
    const rows = Math.max(1, Math.min(H, Math.floor((1 << 18) / Math.max(1, stride)))), nStrips = Math.ceil(H / rows);
    const strips = [];
    for (let k = 0; k < nStrips; k++) {
      const y0 = k * rows, y1 = Math.min(H, y0 + rows), s = new Uint8Array((y1 - y0) * stride);
      for (let y = y0, t = 0; y < y1; y++) {
        for (let x = 0, p = y * W * 4; x < W; x++, p += 4) { s[t++] = src[p]; s[t++] = src[p + 1]; s[t++] = src[p + 2]; if (alpha) s[t++] = src[p + 3]; }
        if (deflated) for (let i = (y - y0 + 1) * stride - 1; i >= (y - y0) * stride + spp; i--) s[i] = (s[i] - s[i - spp]) & 255;   // horizontal predictor
      }
      strips.push(deflated ? zlib(s, o.level === undefined ? 6 : o.level) : s);
    }
    const profile = o.profile ? icc(o.profile) : null, dpi = o.dpi > 0 ? o.dpi : 72;
    const software = ('BlueMolecule ' + VERSION()).trim() + '\0';
    // tags: [tag, type (3 SHORT, 4 LONG, 5 RATIONAL, 2 ASCII, 7 UNDEFINED), values or bytes]
    const tags = [[256, 4, [W]], [257, 4, [H]], [258, 3, new Array(spp).fill(8)], [259, 3, [deflated ? 8 : 1]], [262, 3, [2]],
      [273, 4, strips.map(() => 0)], [277, 3, [spp]], [278, 4, [rows]], [279, 4, strips.map(s => s.length)],
      [282, 5, [[Math.round(dpi * 1000), 1000]]], [283, 5, [[Math.round(dpi * 1000), 1000]]], [284, 3, [1]], [296, 3, [2]], [305, 2, software]];
    if (deflated) tags.push([317, 3, [2]]);
    if (alpha) tags.push([338, 3, [2]]);
    if (profile) tags.push([34675, 7, profile]);
    tags.sort((a, b) => a[0] - b[0]);
    const size = t => (t[1] === 3 ? 2 : t[1] === 5 ? 8 : t[1] === 4 ? 4 : 1) * (t[1] === 2 ? t[2].length : t[2].length);
    const ifdSize = 2 + 12 * tags.length + 4;
    let off = 8 + ifdSize;
    const big = new Map();
    for (const t of tags) if (size(t) > 4) { big.set(t, off); off += (size(t) + 1) & ~1; }
    const stripAt = [];
    for (const s of strips) { stripAt.push(off); off += (s.length + 1) & ~1; }
    tags.find(t => t[0] === 273)[2] = stripAt;
    const b = new Bytes(off + 16);
    b.ascii('II'); b.u16le(42); b.u32le(8);
    b.u16le(tags.length);
    const putVals = t => {
      if (t[1] === 3) for (const v of t[2]) b.u16le(v);
      else if (t[1] === 4) for (const v of t[2]) b.u32le(v);
      else if (t[1] === 5) for (const [nu, de] of t[2]) { b.u32le(nu); b.u32le(de); }
      else if (t[1] === 2) b.ascii(t[2]);
      else b.bytes(t[2]);
    };
    for (const t of tags) {
      b.u16le(t[0]); b.u16le(t[1]); b.u32le(t[1] === 2 ? t[2].length : t[2].length);
      if (big.has(t)) b.u32le(big.get(t));
      else { const s = b.pos; putVals(t); while (b.pos < s + 4) b.u8(0); }
    }
    b.u32le(0);
    for (const t of tags) if (big.has(t)) { while (b.pos < big.get(t)) b.u8(0); putVals(t); }
    strips.forEach((s, k) => { while (b.pos < stripAt[k]) b.u8(0); b.bytes(s); });
    return b.done();
  }

  // ---------------------------------------------------------------- PDF
  function pdf(img, o) {
    o = o || {};
    const W = img.width, H = img.height, src = img.data, alpha = o.alpha !== false && hasAlpha(src), level = o.level === undefined ? 6 : Math.max(1, o.level);
    const dpi = o.dpi > 0 ? o.dpi : 72, wPt = W / dpi * 72, hPt = H / dpi * 72, fmt = v => (Math.round(v * 1000) / 1000).toString();
    const rgb = new Uint8Array(W * H * 3), a = alpha ? new Uint8Array(W * H) : null;
    for (let i = 0, s = 0; i < W * H; i++, s += 4) { rgb[3 * i] = src[s]; rgb[3 * i + 1] = src[s + 1]; rgb[3 * i + 2] = src[s + 2]; if (a) a[i] = src[s + 3]; }
    const objs = [];   // [dict text, stream bytes | null]
    const add = (dict, stream) => { objs.push([dict, stream || null]); return objs.length; };
    const catalog = add(null), pages = add(null), page = add(null);
    const prof = o.profile ? icc(o.profile) : null;
    const iccObj = prof ? add('/N 3 /Alternate /DeviceRGB /Filter /FlateDecode', zlib(prof, 9)) : 0;
    const smask = a ? add(`/Type /XObject /Subtype /Image /Width ${W} /Height ${H} /ColorSpace /DeviceGray /BitsPerComponent 8 /Filter /FlateDecode`, zlib(a, level)) : 0;
    const image = add(`/Type /XObject /Subtype /Image /Width ${W} /Height ${H} /ColorSpace ${iccObj ? `[/ICCBased ${iccObj} 0 R]` : '/DeviceRGB'} /BitsPerComponent 8 /Filter /FlateDecode${smask ? ` /SMask ${smask} 0 R` : ''}`, zlib(rgb, level));
    const content = add('', new TextEncoder().encode(`q ${fmt(wPt)} 0 0 ${fmt(hPt)} 0 0 cm /Im0 Do Q\n`));
    const info = add(`/Producer (BlueMolecule ${VERSION()}) /Title (${String(o.title || 'BlueMolecule image').replace(/[()\\]/g, '\\$&').replace(/[^\x20-\x7e]/g, '?')})`);
    objs[catalog - 1][0] = `/Type /Catalog /Pages ${pages} 0 R`;
    objs[pages - 1][0] = `/Type /Pages /Kids [${page} 0 R] /Count 1`;
    objs[page - 1][0] = `/Type /Page /Parent ${pages} 0 R /MediaBox [0 0 ${fmt(wPt)} ${fmt(hPt)}] /Resources << /XObject << /Im0 ${image} 0 R >> >> /Contents ${content} 0 R`;
    const b = new Bytes(rgb.length / 2 + 4096), at = [];
    b.ascii('%PDF-1.4\n'); b.bytes([37, 0xe2, 0xe3, 0xcf, 0xd3, 10]);
    objs.forEach(([dict, stream], k) => {
      at.push(b.pos);
      b.ascii(`${k + 1} 0 obj\n<< ${dict}${stream ? ` /Length ${stream.length}` : ''} >>\n`);
      if (stream) { b.ascii('stream\n'); b.bytes(stream); b.ascii('\nendstream\n'); }
      b.ascii('endobj\n');
    });
    const xref = b.pos;
    b.ascii(`xref\n0 ${objs.length + 1}\n0000000000 65535 f \n`);
    for (const p of at) b.ascii(String(p).padStart(10, '0') + ' 00000 n \n');
    b.ascii(`trailer\n<< /Size ${objs.length + 1} /Root ${catalog} 0 R /Info ${info} 0 R >>\nstartxref\n${xref}\n%%EOF\n`);
    return b.done();
  }

  const FORMATS = [
    { id: 'png', name: 'PNG', ext: 'png', mime: 'image/png', alpha: true, max: 1e9 },
    { id: 'jpeg', name: 'JPEG', ext: 'jpg', mime: 'image/jpeg', alpha: false, max: 65535 },
    { id: 'webp', name: 'WebP (lossless)', ext: 'webp', mime: 'image/webp', alpha: true, max: 16384 },
    { id: 'tiff', name: 'TIFF', ext: 'tif', mime: 'image/tiff', alpha: true, max: 1e9 },
    { id: 'pdf', name: 'PDF', ext: 'pdf', mime: 'application/pdf', alpha: true, max: 1e9 }
  ];
  /** Encode img in a format of FORMATS (the pixels converted to Display P3 first for that profile). */
  function encode(format, img, o) {
    o = o || {};
    const im = o.profile === 'p3' ? toP3(img) : img;
    switch (format) {
      case 'png': return png(im, o);
      case 'jpeg': return jpeg(im, o);
      case 'webp': return webp(im, o);
      case 'tiff': return tiff(im, o);
      case 'pdf': return pdf(im, o);
    }
    throw new Error('Unknown image format: ' + format);
  }

  BM.imagex = { FORMATS, encode, png, jpeg, webp, tiff, pdf, icc, toP3, deflate, zlib, crc32, adler32, codeLengths, canonical, hasAlpha };
})(typeof window !== 'undefined' ? (window.BM = window.BM || {}) : (globalThis.BM = globalThis.BM || {}));
