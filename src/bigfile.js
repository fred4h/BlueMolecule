/* BlueMolecule — reading files too large for one JavaScript string (BM.bigfile).
 *
 * V8 caps a string at about 512 MiB (2^29 - 24 UTF-16 units), so Blob.text() of a larger output fails and the
 * file looked "empty". Here the bytes are decoded piece by piece with a streaming TextDecoder and split into
 * lines; the readers (parsers.js, readers-qc.js, readers-struct.js) accept such a line array wherever they
 * accept text, and BM.parseLines(fileName, lines) detects the format and parses it.
 *
 *   BM.bigfile.readLines(source, {onProgress, signal, chunkSize, total, memoryLimit}) -> Promise<string[]>
 *       source: a Blob/File (read in 32 MB slices), a ReadableStream of bytes (fetch body, DecompressionStream),
 *       a {name, size, readChunk(offset, length) -> Promise<Uint8Array>} object (native shells), a Uint8Array /
 *       ArrayBuffer, or a string. The lines are exactly those of text.split(/\r\n|\r|\n/) (a leading BOM is
 *       dropped, invalid UTF-8 becomes U+FFFD), also when a CR LF pair or a multi-byte character straddles two
 *       pieces. onProgress({bytes, total, lines}) runs after every piece; total is null when unknown (a stream
 *       without opts.total). signal is an AbortSignal: the promise then rejects with an Error named 'AbortError'.
 *       When the JavaScript heap can be measured (Node; Chromium with precise memory info) reading stops with an
 *       Error named 'OutOfMemoryError' before the heap would run out, instead of crashing the page.
 *   BM.bigfile.readBytes(source, {onProgress, signal}) -> Promise<Uint8Array>     the whole content
 *   BM.bigfile.readHead(source, n = 65536) -> Promise<string>                      the first n bytes as text
 *   BM.bigfile.byteReader(source) -> {size, read(offset, length) -> Promise<Uint8Array>} | null   random access
 *   BM.bigfile.toStream(source) -> ReadableStream;  BM.bigfile.gunzip(source, onBytes) -> ReadableStream
 *       (onBytes(compressedBytesSoFar) gives progress on the compressed side, whose size is known)
 *   BM.bigfile.toSource(item, name) -> {name, size, kind: 'blob'|'chunks'|'bytes'|'stream'|'text', data}
 *       normalises what a host hands over (File/Blob, {name, bytes|data|base64|text}, chunk object, stream)
 *   BM.bigfile.isAbort(e), isOutOfMemory(e), describeError(e, name, size)   error helpers for the UI
 *   BM.bigfile.estimateMemory(size) -> bytes of memory a text file of `size` bytes needs once split into lines
 *   BM.bigfile.CHUNK (32 MB), BIG (64 MB: from this size the app reads files with readLines)
 */
(function (BM) {
  'use strict';
  const MB = 1048576;
  const CHUNK = 32 * MB;
  const BIG = 64 * MB;

  // ------------------------------------------------------------------ errors
  function abortError() { const e = new Error('Reading was cancelled.'); e.name = 'AbortError'; return e; }
  function isAbort(e) { return !!e && (e.name === 'AbortError' || e.code === 20); }
  function oomError(msg) { const e = new Error(msg); e.name = 'OutOfMemoryError'; return e; }
  /** True for the errors a JavaScript engine raises when memory or string/array limits are exhausted. */
  function isOutOfMemory(e) {
    if (!e) return false;
    if (e.name === 'OutOfMemoryError') return true;
    const m = String(e.message || e);
    return /Invalid string length|Invalid array length|Invalid typed array length|Array buffer allocation failed|allocation failed|out of memory/i.test(m) ||
      (e.name === 'RangeError' && /alloc/i.test(m));
  }
  const fmtMB = b => (b >= 1e9 ? (b / 1073741824).toFixed(2) + ' GB' : Math.round(b / MB) + ' MB');
  /** A message for the user: cancellation, memory exhaustion (with what to do), or the error itself. */
  function describeError(e, name, size) {
    const what = '"' + (name || 'the file') + '"' + (size ? ' (' + fmtMB(size) + ')' : '');
    if (isAbort(e)) return 'Opening ' + what + ' was cancelled.';
    if (isOutOfMemory(e)) {
      return 'Not enough memory to open ' + what + '. ' + (e && e.name === 'OutOfMemoryError' ? e.message + ' ' : '') +
        'Close other documents or tabs and try again, or open a smaller part of the file (for example the last ' +
        'optimisation cycles, cut out with a text tool).';
    }
    return (e && e.message) || String(e);
  }

  // ------------------------------------------------------------------ sources
  const isBlob = x => !!x && typeof x === 'object' && typeof x.slice === 'function' && typeof x.arrayBuffer === 'function' && typeof x.size === 'number';
  const isChunks = x => !!x && typeof x === 'object' && typeof x.readChunk === 'function';
  const isStream = x => !!x && typeof x === 'object' && typeof x.getReader === 'function';
  const isBytes = x => !!x && (x instanceof Uint8Array || Object.prototype.toString.call(x) === '[object ArrayBuffer]' || ArrayBuffer.isView(x));
  const asU8 = x => (x instanceof Uint8Array ? x : ArrayBuffer.isView(x) ? new Uint8Array(x.buffer, x.byteOffset, x.byteLength) : new Uint8Array(x));
  function sizeOf(src) {
    if (isBlob(src)) return src.size;
    if (isChunks(src)) return Number.isFinite(+src.size) ? +src.size : null;
    if (isBytes(src)) return asU8(src).byteLength;
    return null;
  }
  /** Random access to a Blob, chunk object or byte array; null for streams. */
  function byteReader(src) {
    if (src && src.kind && src.data !== undefined && !isBlob(src) && !isChunks(src)) src = src.data;   // toSource() record
    if (isBlob(src)) return { size: src.size, read: async (o, n) => new Uint8Array(await src.slice(o, o + n).arrayBuffer()) };
    if (isChunks(src)) {
      return { size: sizeOf(src), read: async (o, n) => { const r = await src.readChunk(o, n); return r == null ? new Uint8Array(0) : asU8(r); } };
    }
    if (isBytes(src)) { const u8 = asU8(src); return { size: u8.byteLength, read: async (o, n) => u8.subarray(o, Math.min(u8.byteLength, o + n)) }; }
    return null;
  }
  /** A ReadableStream of the source's bytes (Blob.stream(), or pieces of a chunk object). */
  function toStream(src, chunkSize) {
    if (isStream(src)) return src;
    if (isBlob(src) && typeof src.stream === 'function') return src.stream();
    const R = byteReader(src);
    if (!R) throw new TypeError('BM.bigfile: unsupported source');
    let off = 0;
    const step = chunkSize || 4 * MB;
    return new ReadableStream({
      async pull(ctl) {
        if (R.size !== null && off >= R.size) { ctl.close(); return; }
        const u8 = await R.read(off, R.size === null ? step : Math.min(step, R.size - off));
        if (!u8.byteLength) { ctl.close(); return; }
        off += u8.byteLength;
        ctl.enqueue(u8);
      }
    });
  }
  /** Decompressed .gz content as a stream; onBytes(n) reports the compressed bytes consumed so far. */
  function gunzip(src, onBytes) {
    if (typeof DecompressionStream === 'undefined') throw new Error('This browser cannot decompress .gz files; please unzip the file first.');
    let n = 0;
    const count = new TransformStream({ transform(chunk, ctl) { n += chunk.byteLength; if (onBytes) onBytes(n); ctl.enqueue(chunk); } });
    return toStream(src).pipeThrough(count).pipeThrough(new DecompressionStream('gzip'));
  }
  function base64ToBytes(b64) {
    const s = atob(String(b64 || '').replace(/\s+/g, ''));
    const out = new Uint8Array(s.length);
    for (let i = 0; i < s.length; i++) out[i] = s.charCodeAt(i);
    return out;
  }
  /**
   * Normalises a file handed over by a host or a test: File/Blob, chunk object, ReadableStream, bytes, or a plain
   * record {name, bytes | data | buffer | base64 | text | blob | file | stream, size}. Returns
   * {name, size (null if unknown), kind: 'blob'|'chunks'|'bytes'|'stream'|'text', data}.
   */
  function toSource(item, name) {
    if (item && item.kind && 'data' in item && ['blob', 'chunks', 'bytes', 'stream', 'text'].indexOf(item.kind) >= 0 && !isBlob(item) && !isChunks(item)) return item;
    const nm = String(name || (item && (item.name || item.fileName || item.filename)) || 'untitled');
    if (typeof item === 'string') return { name: nm, size: item.length, kind: 'text', data: item };
    if (isBlob(item)) return { name: nm, size: item.size, kind: 'blob', data: item };
    if (isChunks(item)) return { name: nm, size: sizeOf(item), kind: 'chunks', data: item };
    if (isStream(item)) return { name: nm, size: null, kind: 'stream', data: item };
    if (isBytes(item)) { const u8 = asU8(item); return { name: nm, size: u8.byteLength, kind: 'bytes', data: u8 }; }
    if (item && typeof item === 'object') {
      for (const k of ['file', 'blob']) if (isBlob(item[k])) return { name: nm, size: item[k].size, kind: 'blob', data: item[k] };
      for (const k of ['bytes', 'data', 'buffer', 'content']) {
        const v = item[k];
        if (isBytes(v)) { const u8 = asU8(v); return { name: nm, size: u8.byteLength, kind: 'bytes', data: u8 }; }
        if (isBlob(v)) return { name: nm, size: v.size, kind: 'blob', data: v };
      }
      if (typeof item.base64 === 'string') { const u8 = base64ToBytes(item.base64); return { name: nm, size: u8.byteLength, kind: 'bytes', data: u8 }; }
      if (typeof item.text === 'string') return { name: nm, size: item.text.length, kind: 'text', data: item.text };
      if (isStream(item.stream)) return { name: nm, size: Number.isFinite(+item.size) ? +item.size : null, kind: 'stream', data: item.stream };
    }
    throw new TypeError('Unsupported file object for "' + nm + '".');
  }

  // ------------------------------------------------------------------ memory
  function heapStats() {
    try {
      const m = typeof performance !== 'undefined' && performance.memory;
      if (m && m.jsHeapSizeLimit > 0 && m.usedJSHeapSize > 0) return { used: m.usedJSHeapSize, limit: m.jsHeapSizeLimit, precise: !!BM.bigfile.preciseMemory };
    } catch (e) { /* not available */ }
    try {
      if (typeof process !== 'undefined' && process.versions && process.versions.node && typeof require === 'function') {
        const h = require('v8').getHeapStatistics();
        return { used: h.used_heap_size, limit: h.heap_size_limit, precise: true };
      }
    } catch (e) { /* not available */ }
    return null;
  }
  /** Approximate memory (bytes) a text file of `size` bytes takes as a line array: the decoded text (1 byte per
   *  character for ASCII) plus one string header and array slot per line (~60-character lines). */
  function estimateMemory(size) { return Math.round(size * (1 + 40 / 60)); }
  /**
   * Watches the heap while lines accumulate and stops early (clear message) when the projected heap at the end of
   * the file, plus room for the reader's own data, would exceed the limit. Without precise numbers (Chromium
   * quantises performance.memory and refreshes it rarely) only the static limit is checked once.
   */
  function MemoryGuard(total, opts) {
    this.total = total; this.h0 = heapStats(); this.limitOverride = opts && opts.memoryLimit;
    this.reserve = total ? Math.min(total * 0.25, 512 * MB) : 64 * MB;   // the reader's frames, tables, temporaries
  }
  MemoryGuard.prototype.check = function (bytes, final) {
    const h0 = this.h0;
    if (this.limitOverride && this.total && estimateMemory(this.total) > this.limitOverride) {
      throw oomError('The file needs about ' + fmtMB(estimateMemory(this.total)) + ' of memory and ' + fmtMB(this.limitOverride) + ' is available.');
    }
    if (!h0 || !h0.precise) return;
    const h = heapStats();
    if (!h) return;
    const limit = this.limitOverride || h.limit;
    const grown = Math.max(0, h.used - h0.used);
    const rest = this.total && !final ? Math.max(0, this.total - bytes) : 0;
    const perByte = bytes > 8 * MB ? grown / bytes : 1.7;
    const projected = h.used + perByte * rest + this.reserve;
    if (h.used > 0.9 * limit || (projected > 0.97 * limit && bytes > 16 * MB)) {
      throw oomError('After ' + fmtMB(bytes) + (this.total ? ' of ' + fmtMB(this.total) : '') + ' the memory in use is ' + fmtMB(h.used) +
        ' of at most ' + fmtMB(limit) + '; the whole file would need about ' + fmtMB(projected) + '.');
    }
  };

  // ------------------------------------------------------------------ line splitting
  /** Splits decoded text pieces into lines; the result equals text.split(/\r\n|\r|\n/) of the whole text. */
  function LineSplitter(out) { this.out = out; this.carry = ''; this.cr = false; }
  LineSplitter.prototype.feed = function (text) {
    if (!text) return;
    if (this.cr) {                       // the previous piece ended with CR: a LF right here belongs to it
      this.cr = false;
      if (text.charCodeAt(0) === 10) { text = text.slice(1); if (!text) return; }
    }
    let endsCR = false;
    if (text.charCodeAt(text.length - 1) === 13) { endsCR = true; text = text.slice(0, -1); }
    let parts;
    if (text.indexOf('\r') < 0) parts = text.split('\n');
    else {
      parts = text.split('\r\n');
      let fixed = null;                  // lone CR or LF left inside a part (mixed or classic-Mac line ends)
      for (let k = 0; k < parts.length; k++) {
        const p = parts[k];
        if (p.indexOf('\n') >= 0 || p.indexOf('\r') >= 0) {
          if (!fixed) fixed = parts.slice(0, k);
          const q = p.split(/\r|\n/);
          for (let m = 0; m < q.length; m++) fixed.push(q[m]);
        } else if (fixed) fixed.push(p);
      }
      if (fixed) parts = fixed;
    }
    const out = this.out;
    if (this.carry) parts[0] = this.carry + parts[0];
    const end = endsCR ? parts.length : parts.length - 1;
    for (let k = 0; k < end; k++) out.push(parts[k]);
    this.carry = endsCR ? '' : parts[parts.length - 1];
    if (endsCR) this.cr = true;
  };
  LineSplitter.prototype.finish = function () { this.out.push(this.carry); this.carry = ''; return this.out; };

  // Node's TextDecoder returns two-byte strings (2 bytes per character); pure-ASCII pieces are decoded as Latin-1
  // there, which gives one-byte strings (browsers already produce those for Latin-1 text).
  const NODE = typeof process !== 'undefined' && !!(process.versions && process.versions.node) && typeof Buffer === 'function';
  function isAscii(u8) {
    const n = u8.byteLength;
    let i = 0;
    if ((u8.byteOffset & 3) === 0 && n >= 64) {
      const w = new Uint32Array(u8.buffer, u8.byteOffset, n >> 2);
      for (let k = 0; k < w.length; k++) if (w[k] & 0x80808080) return false;
      i = w.length << 2;
    }
    for (; i < n; i++) if (u8[i] & 0x80) return false;
    return true;
  }
  /** Streaming UTF-8 decoding; a leading BOM is dropped (once, at the very start). */
  function Decoder() { this.dec = new TextDecoder('utf-8', { ignoreBOM: true }); this.used = false; this.start = true; }
  Decoder.prototype.piece = function (u8, final) {
    let t;
    if (final) t = this.used ? this.dec.decode() : '';
    else if (NODE && u8.byteLength >= 4096 && isAscii(u8)) {
      t = (this.used ? this.dec.decode() : '') + Buffer.from(u8.buffer, u8.byteOffset, u8.byteLength).toString('latin1');
      this.used = false;
    } else { t = this.dec.decode(u8, { stream: true }); this.used = true; }
    if (this.start && t) { this.start = false; if (t.charCodeAt(0) === 0xFEFF) t = t.slice(1); }
    return t;
  };

  const now = () => (typeof performance !== 'undefined' && performance.now ? performance.now() : Date.now());
  /** Lets timers and UI events (a Cancel click) run when a source answers without real I/O (in-memory streams). */
  function Breather() { this.t = now(); }
  Breather.prototype.maybe = function () {
    if (now() - this.t < 50) return null;
    return new Promise(r => setTimeout(r, 0)).then(() => { this.t = now(); });
  };

  /** Pieces of text from any source, decoded with one streaming decoder. */
  async function eachPiece(source, opts, fn) {
    const signal = opts.signal;
    const dec = new Decoder(), breathe = new Breather();
    let bytes = 0;
    const src = source && source.kind && 'data' in source && !isBlob(source) && !isChunks(source) ? source.data : source;
    const total = sizeOf(src) !== null ? sizeOf(src) : (Number.isFinite(+opts.total) && opts.total > 0 ? +opts.total : null);
    if (signal && signal.aborted) throw abortError();
    if (isStream(src)) {
      const reader = src.getReader();
      let done = false;
      try {
        for (;;) {
          const b = breathe.maybe();
          if (b) await b;
          if (signal && signal.aborted) throw abortError();
          const r = await reader.read();
          if (r.done) { done = true; break; }
          const u8 = asU8(r.value);
          bytes += u8.byteLength;
          fn(dec.piece(u8, false), bytes, total, false);
        }
      } finally {
        if (!done) { try { await reader.cancel(); } catch (e) { /* already closed */ } }
        try { reader.releaseLock(); } catch (e) { /* ignore */ }
      }
    } else {
      const R = byteReader(src);
      if (!R) throw new TypeError('BM.bigfile.readLines: unsupported source (use a Blob/File, ReadableStream, chunk object or bytes).');
      const step = Math.max(1, opts.chunkSize || CHUNK);
      for (let off = 0; R.size === null || off < R.size;) {
        const b = breathe.maybe();
        if (b) await b;
        if (signal && signal.aborted) throw abortError();
        const want = R.size === null ? step : Math.min(step, R.size - off);
        const u8 = await R.read(off, want);
        if (!u8.byteLength) {
          if (R.size === null) break;
          throw new Error('The file could not be read past byte ' + off + ' of ' + R.size + '.');
        }
        off += u8.byteLength; bytes = off;
        fn(dec.piece(u8, false), bytes, total, false);
      }
    }
    if (signal && signal.aborted) throw abortError();
    fn(dec.piece(null, true), bytes, total, true);
    return { bytes, total };
  }

  async function readLines(source, opts) {
    opts = opts || {};
    if (typeof source === 'string') return source.split(/\r\n|\r|\n/);
    if (source && source.kind === 'text') return String(source.data).split(/\r\n|\r|\n/);
    const lines = [];
    const sp = new LineSplitter(lines);
    const src0 = source && source.kind && 'data' in source && !isBlob(source) && !isChunks(source) ? source.data : source;
    const guard = new MemoryGuard(sizeOf(src0) !== null ? sizeOf(src0) : (+opts.total || 0), opts);
    let lastReport = 0;
    const r = await eachPiece(source, opts, (text, bytes, total, final) => {
      sp.feed(text);
      if (final) return;
      guard.check(bytes, false);
      if (opts.onProgress && (bytes - lastReport >= 4 * MB || (total && bytes >= total))) {
        lastReport = bytes;
        opts.onProgress({ bytes, total, lines: lines.length });
      }
    });
    sp.finish();
    if (opts.onProgress) opts.onProgress({ bytes: r.bytes, total: r.total === null ? r.bytes : r.total, lines: lines.length, done: true });
    return lines;
  }

  async function readBytes(source, opts) {
    opts = opts || {};
    const src = source && source.kind && 'data' in source && !isBlob(source) && !isChunks(source) ? source.data : source;
    if (typeof src === 'string') return new TextEncoder().encode(src);
    if (isBytes(src)) return asU8(src);
    if (isBlob(src) && !opts.onProgress) return new Uint8Array(await src.arrayBuffer());
    const R = isStream(src) ? null : byteReader(src);
    if (R && R.size !== null) {
      const out = new Uint8Array(R.size);
      for (let off = 0; off < R.size;) {
        if (opts.signal && opts.signal.aborted) throw abortError();
        const u8 = await R.read(off, Math.min(opts.chunkSize || CHUNK, R.size - off));
        if (!u8.byteLength) throw new Error('The file could not be read past byte ' + off + '.');
        out.set(u8, off); off += u8.byteLength;
        if (opts.onProgress) opts.onProgress({ bytes: off, total: R.size });
      }
      return out;
    }
    const parts = [];
    let n = 0;
    const reader = toStream(src).getReader();
    for (;;) {
      if (opts.signal && opts.signal.aborted) { try { await reader.cancel(); } catch (e) { /* ignore */ } throw abortError(); }
      const r = await reader.read();
      if (r.done) break;
      const u8 = asU8(r.value);
      parts.push(u8); n += u8.byteLength;
      if (opts.onProgress) opts.onProgress({ bytes: n, total: Number.isFinite(+opts.total) ? +opts.total : null });
    }
    const out = new Uint8Array(n);
    let o = 0;
    for (const p of parts) { out.set(p, o); o += p.byteLength; }
    return out;
  }

  async function readHead(source, n) {
    n = n || 65536;
    const src = source && source.kind && 'data' in source && !isBlob(source) && !isChunks(source) ? source.data : source;
    if (typeof src === 'string') return src.slice(0, n);
    const R = byteReader(src);
    if (!R) throw new TypeError('BM.bigfile.readHead needs a Blob, chunk object or bytes.');
    const u8 = await R.read(0, R.size === null ? n : Math.min(n, R.size));
    return new TextDecoder('utf-8').decode(u8);   // a character cut at the end becomes U+FFFD: harmless for sniffing
  }

  BM.bigfile = {
    CHUNK, BIG, readLines, readBytes, readHead, byteReader, toStream, gunzip, toSource, sizeOf,
    isAbort, isOutOfMemory, describeError, estimateMemory, heapStats, LineSplitter,
    preciseMemory: false   // set by a host that runs Chromium with --enable-precise-memory-info
  };
})(typeof window !== 'undefined' ? (window.BM = window.BM || {}) : (globalThis.BM = globalThis.BM || {}));
