/* BlueMolecule — native bridge (BM.native).            BM:native-bridge
 *
 * Talks to whatever hosts the page:
 *   'macos-app' — the native macOS app (WKWebView message handlers bmReady, bmSave, bmObabel, bmCopy;
 *                 replies come back through BM.native._resolve(id, result) and documents opened from
 *                 Finder/the Dock through BM.native._receiveFile(name, base64));
 *   'linux-app' — the Linux Qt WebEngine window (native/linux/bluemolecule_qt.py, which defines
 *                 window.BM_QT before the page runs). The same messages plus bmRead go through
 *                 QWebChannel as channel.objects.bmHost.post(name, JSON text); replies again come through
 *                 _resolve. Documents named on the command line arrive as _receiveFile(name, base64) or,
 *                 from 32 MB on, as _receiveLargeFile(name, size, fileId), read on demand with
 *                 bmRead {id, fileId, offset, length} → {ok, base64, offset, length, eof} (≤ 8 MB of base64);
 *   'launcher'  — the Python launcher (window.BM_LAUNCH; HTTP endpoints on 127.0.0.1);
 *   'browser'   — a plain browser: save/copy report "not handled" (null) so the UI uses its own paths.
 *
 * Open Babel uses one contract on all hosts: request {op: 'formats'|'convert', input (base64),
 * inFormat, outFormat, options[]}, reply {ok, output (base64), error}. 'formats' returns the text of
 * `obabel -L formats`, which is parsed here.
 */
(function (BM) {
  'use strict';
  if (BM.native) return;
  const G = typeof window !== 'undefined' ? window : globalThis;
  const VERSION = '1.2.0';
  const FORMAT_RE = /^[A-Za-z0-9_-]{1,16}$/;

  const handlers = () => (G.webkit && G.webkit.messageHandlers) || null;
  const QT = G.BM_QT && typeof G.BM_QT === 'object' ? G.BM_QT : null;
  const kind = (() => {
    const h = handlers();
    if (h && h.bmSave) return 'macos-app';
    if (QT) return 'linux-app';
    if (G.BM_LAUNCH) return 'launcher';
    return 'browser';
  })();
  /** Bytes per bmRead reply: 6 MiB, i.e. 8 MiB of base64 (the Linux app may lower it). */
  const READ_LIMIT = QT && Number.isSafeInteger(QT.readLimit) && QT.readLimit > 0 ? Math.min(QT.readLimit, 6291456) : 6291456;

  // ---------------- bytes, base64, text ----------------
  const isBuffer = x => Object.prototype.toString.call(x) === '[object ArrayBuffer]';
  async function toBytes(data) {
    if (data == null) return new Uint8Array(0);
    if (typeof data === 'string') return new TextEncoder().encode(data);
    if (ArrayBuffer.isView(data)) return new Uint8Array(data.buffer, data.byteOffset, data.byteLength);
    if (isBuffer(data)) return new Uint8Array(data);
    if (typeof data.arrayBuffer === 'function') return new Uint8Array(await data.arrayBuffer());
    throw new TypeError('BM.native: unsupported data (use a Blob, ArrayBuffer, typed array or string)');
  }
  function bytesToBase64(bytes) {
    let s = '';
    for (let i = 0; i < bytes.length; i += 0x8000) s += String.fromCharCode.apply(null, bytes.subarray(i, i + 0x8000));
    return btoa(s);
  }
  function base64ToBytes(b64) {
    const s = atob(String(b64 || '').replace(/\s+/g, ''));
    const out = new Uint8Array(s.length);
    for (let i = 0; i < s.length; i++) out[i] = s.charCodeAt(i);
    return out;
  }
  const bytesToText = bytes => new TextDecoder('utf-8').decode(bytes);

  // ---------------- the Linux app's QWebChannel ----------------
  // The channel opens asynchronously; messages wait in qtQueue until it is open. If it never opens,
  // waiting requests are answered with an error and save/copy fall back to the page's own paths.
  let qtHost = null, qtState = 'idle';          // idle → connecting → open | failed
  const qtQueue = [];
  function qtFail(why) {
    qtState = 'failed';
    console.error('BlueMolecule: ' + why);
    while (qtQueue.length) { const m = qtQueue.shift(); if (m.id) _resolve(m.id, { ok: false, error: why }); }
  }
  function qtConnect() {
    if (qtState !== 'idle') return;
    qtState = 'connecting';
    const limit = Number(QT.connectTimeout) > 0 ? Number(QT.connectTimeout) : 10000;
    const t0 = Date.now();
    const attempt = () => {
      const transport = G.qt && G.qt.webChannelTransport;
      if (transport && typeof G.QWebChannel === 'function') {
        try {
          new G.QWebChannel(transport, ch => {
            const host = ch && ch.objects && ch.objects.bmHost;
            if (!host || typeof host.post !== 'function') { qtFail('the BlueMolecule window offers no bmHost object'); return; }
            qtHost = host;
            qtState = 'open';
            while (qtQueue.length) { const m = qtQueue.shift(); host.post(m.name, m.text); }
          });
        } catch (e) { qtFail('could not open the channel to the BlueMolecule window: ' + ((e && e.message) || e)); }
        return;
      }
      if (Date.now() - t0 > limit) { qtFail('the channel to the BlueMolecule window did not appear'); return; }
      setTimeout(attempt, 20);
    };
    attempt();
  }
  /** Whether the host's own Save dialog, clipboard and messages can be used. */
  const appBridge = () => kind === 'macos-app' || (kind === 'linux-app' && qtState !== 'failed');

  // ---------------- requests to the app ----------------
  let seq = 0;
  const pending = new Map();
  function post(name, msg) {
    if (kind === 'linux-app') {
      if (qtState === 'failed') throw new Error('the channel to the BlueMolecule window is not available');
      const text = JSON.stringify(msg);
      if (qtHost) qtHost.post(name, text);
      else { qtQueue.push({ name, text, id: msg && msg.id }); qtConnect(); }
      return;
    }
    const h = handlers();
    if (!h || !h[name]) throw new Error('the BlueMolecule app has no "' + name + '" handler');
    h[name].postMessage(msg);
  }
  /** Post {id, ...payload} to the app; resolves with the result the app passes to _resolve(id, result). */
  function call(name, payload, timeoutMs) {
    return new Promise(resolve => {
      const id = 'bm' + (++seq);
      let timer = 0;
      pending.set(id, r => { if (timer) clearTimeout(timer); resolve(r); });
      if (timeoutMs) timer = setTimeout(() => { if (pending.delete(id)) resolve({ ok: false, error: 'No reply from the BlueMolecule app.' }); }, timeoutMs);
      try { post(name, Object.assign({ id }, payload)); }
      catch (e) { pending.delete(id); if (timer) clearTimeout(timer); resolve({ ok: false, error: String((e && e.message) || e) }); }
    });
  }
  function _resolve(id, result) {
    const fn = pending.get(String(id));
    if (!fn) return false;
    pending.delete(String(id));
    fn(result && typeof result === 'object' ? result : { ok: !!result });
    return true;
  }

  // ---------------- save / copy / log ----------------
  /** Save through the host's own Save dialog. true = saved, false = cancelled or failed, null = not handled here. */
  async function save(name, data) {
    if (!appBridge()) return null;
    const bytes = await toBytes(data);
    const mime = (data && typeof data.type === 'string' && data.type) || 'application/octet-stream';
    const r = await call('bmSave', { name: String(name || 'untitled'), mime, base64: bytesToBase64(bytes) });
    native.lastSave = r;
    return !!(r && r.ok);
  }
  /** Plain text to the system clipboard. true/false, or null when the host has no native clipboard. */
  function copyText(text) {
    if (!appBridge()) return Promise.resolve(null);
    try { post('bmCopy', { text: String(text == null ? '' : text) }); return Promise.resolve(true); }
    catch (e) { return Promise.resolve(false); }
  }
  /** A line for the app's log file (bmLog {level, message}). true/false, or null without an app. */
  function log(level, message) {
    if (!appBridge()) return null;
    try { post('bmLog', { level: String(level || 'log'), message: String(message == null ? '' : message).slice(0, 8000) }); return true; }
    catch (e) { return false; }
  }

  // ---------------- Open Babel ----------------
  async function obabelRequest(payload) {
    let r;
    if (kind === 'macos-app' || kind === 'linux-app') r = await call('bmObabel', payload, 90000);
    else if (kind === 'launcher' && G.BM_LAUNCH.obabel) {
      const resp = await G.fetch(G.BM_LAUNCH.obabel, { method: 'POST', headers: { 'Content-Type': 'application/json' }, body: JSON.stringify(payload) });
      if (!resp.ok) throw new Error('Open Babel request failed (HTTP ' + resp.status + ')');
      r = await resp.json();
    } else throw new Error('Open Babel needs the BlueMolecule app or the bluemolecule launcher.');
    if (!r || !r.ok) throw new Error((r && r.error) || 'Open Babel failed.');
    return r;
  }
  /** Parse `obabel -L formats` ("xyz -- XYZ cartesian coordinates format [Read-only]") into records. */
  function parseFormats(text) {
    const out = [], seen = new Set();
    for (const line of String(text || '').split(/\r?\n/)) {
      const m = /^\s*([A-Za-z0-9_+-]+)\s+--\s+(.*?)\s*$/.exec(line);
      if (!m || seen.has(m[1].toLowerCase())) continue;
      let description = m[2], read = true, write = true;
      if (/\[read-only\]$/i.test(description)) write = false;
      else if (/\[write-only\]$/i.test(description)) read = false;
      description = description.replace(/\s*\[(read|write)-only\]$/i, '');
      seen.add(m[1].toLowerCase());
      out.push({ code: m[1], description, read, write });
    }
    return out;
  }
  let formatsPromise = null;
  const obabel = {
    /** Whether Open Babel can be used from this page (the host found an obabel executable). */
    available() {
      if (kind === 'browser') return Promise.resolve(false);
      return obabel.formats().then(f => f.length > 0, () => false);
    },
    /** [{code, description, read, write}], cached after the first success. */
    formats() {
      if (!formatsPromise) {
        formatsPromise = obabelRequest({ op: 'formats' }).then(r => parseFormats(bytesToText(base64ToBytes(r.output))));
        formatsPromise.catch(() => { formatsPromise = null; });
      }
      return formatsPromise;
    },
    /** Convert data (string, bytes or Blob) from inFmt to outFmt; resolves with the output text. */
    async convert(data, inFmt, outFmt, options) {
      if (!FORMAT_RE.test(String(inFmt)) || !FORMAT_RE.test(String(outFmt))) throw new Error('Invalid Open Babel format code: ' + inFmt + ' → ' + outFmt);
      const opts = options == null ? [] : Array.from(options, String);
      const input = bytesToBase64(await toBytes(data));
      const r = await obabelRequest({ op: 'convert', input, inFormat: String(inFmt), outFormat: String(outFmt), options: opts });
      return bytesToText(base64ToBytes(r.output));
    }
  };

  // ---------------- documents opened by the host ----------------
  /** A large document the Linux app keeps on disk: {name, size, readChunk(offset, length) → Promise<Uint8Array>}.
   *  readChunk returns the bytes [offset, min(offset + length, size)), fetched in bmRead pieces of ≤ READ_LIMIT. */
  function fileHandle(name, size, fileId) {
    async function readChunk(offset, length) {
      if (!Number.isSafeInteger(offset) || !Number.isSafeInteger(length) || offset < 0 || length < 0) {
        throw new RangeError('readChunk(offset, length) needs whole numbers ≥ 0');
      }
      const end = Math.min(size, offset + length);
      if (end <= offset) return new Uint8Array(0);
      const out = new Uint8Array(end - offset);
      let pos = offset;
      while (pos < end) {
        const r = await call('bmRead', { fileId, offset: pos, length: Math.min(READ_LIMIT, end - pos) }, 120000);
        if (!r || !r.ok) throw new Error((r && r.error) || 'Could not read ' + name);
        const bytes = base64ToBytes(r.base64);
        if (!bytes.length) break;                       // the file became shorter since it was opened
        const take = Math.min(bytes.length, end - pos);
        out.set(take === bytes.length ? bytes : bytes.subarray(0, take), pos - offset);
        pos += take;
      }
      return pos === end ? out : out.slice(0, pos - offset);
    }
    return Object.freeze({ name, size, readChunk });
  }
  const isHandle = x => !!x && typeof x === 'object' && typeof x.readChunk === 'function';

  const openCallbacks = [];
  const inbox = [];
  /** Small documents go to cb(name, Uint8Array), large ones to cb({name, size, readChunk}). */
  function deliver(name, data) {
    if (openCallbacks.length) {
      for (const cb of openCallbacks.slice()) { try { if (isHandle(data)) cb(data); else cb(name, data); } catch (e) { console.error(e); } }
      return true;
    }
    const app = BM.app;   // nobody registered: hand a small file to the viewer directly (large ones wait)
    if (!isHandle(data) && app && typeof app.openFile === 'function' && typeof G.File === 'function') {
      Promise.resolve(app.openFile(new G.File([data], name))).catch(e => console.error(e));
      return true;
    }
    return false;
  }
  function flushInbox() { while (inbox.length && deliver(inbox[0][0], inbox[0][1])) inbox.shift(); }
  /** cb for every document the host opens (Finder, Dock, command line): cb(name, Uint8Array), or, for a
   *  large file of the Linux app, cb({name, size, readChunk}). Returns an unsubscribe function. */
  function onOpenFile(cb) {
    if (typeof cb !== 'function') throw new TypeError('BM.native.onOpenFile needs a function');
    openCallbacks.push(cb);
    flushInbox();
    return () => { const i = openCallbacks.indexOf(cb); if (i >= 0) openCallbacks.splice(i, 1); };
  }
  function _receiveFile(name, base64) {
    inbox.push([String(name || 'untitled'), base64ToBytes(base64)]);
    flushInbox();
    return true;
  }
  function _receiveLargeFile(name, size, fileId) {
    const n = Number(size);
    if (!Number.isSafeInteger(n) || n < 0 || fileId == null || fileId === '') return false;
    const nm = String(name || 'untitled');
    inbox.push([nm, fileHandle(nm, n, String(fileId))]);
    flushInbox();
    return true;
  }

  // ---------------- host-specific start-up ----------------
  function announceReady() {
    let tries = 0;
    const tick = () => {
      if (BM.app || ++tries > 300) {   // poll for up to 15 s, then let the app deliver anyway
        flushInbox();
        try { post('bmReady', { version: VERSION, app: !!BM.app }); } catch (e) { /* nothing to tell */ }
        return;
      }
      setTimeout(tick, 50);
    };
    tick();
  }
  /** App safety nets for code that does not ask BM.native first: an <a download> click (WKWebView has no
   *  download UI) goes to the Save dialog and navigator.clipboard.writeText to the system clipboard. */
  function installAppShims() {
    const blobs = new Map();
    const U = G.URL;
    if (U && typeof U.createObjectURL === 'function') {
      const create = U.createObjectURL, revoke = U.revokeObjectURL;
      U.createObjectURL = function (obj) {
        const url = create.apply(U, arguments);
        if (obj && typeof obj.arrayBuffer === 'function') { blobs.set(url, obj); if (blobs.size > 24) blobs.delete(blobs.keys().next().value); }
        return url;
      };
      if (typeof revoke === 'function') U.revokeObjectURL = function (url) { blobs.delete(url); return revoke.apply(U, arguments); };
    }
    const A = G.HTMLAnchorElement;
    if (A && A.prototype && typeof A.prototype.click === 'function') {
      const click = A.prototype.click;
      A.prototype.click = function () {
        const href = String(this.href || '');
        if (appBridge() && this.hasAttribute && this.hasAttribute('download') && /^(blob|data):/i.test(href)) {
          const name = this.getAttribute('download') || 'download';
          const blob = blobs.get(href);
          (blob ? Promise.resolve(blob) : G.fetch(href).then(r => r.blob()))
            .then(b => save(name, b))
            .catch(e => console.error('BlueMolecule: could not save ' + name + ': ' + ((e && e.message) || e)));
          return undefined;
        }
        return click.apply(this, arguments);
      };
    }
    const nav = G.navigator;
    if (nav) {
      const original = nav.clipboard && typeof nav.clipboard.writeText === 'function' ? nav.clipboard.writeText.bind(nav.clipboard) : null;
      const writeText = text => copyText(text).then(ok => {
        if (ok === null && original) return original(text);
        if (!ok) throw new Error('Copy failed');
      });
      try {
        if (nav.clipboard) nav.clipboard.writeText = writeText;
        else Object.defineProperty(nav, 'clipboard', { value: { writeText }, configurable: true });
      } catch (e) { /* keep the page's own fallbacks */ }
    }
  }
  /** Launcher: tell the server when the window is hidden (browsers throttle timers in hidden pages, so
   *  heartbeats may come only once a minute) and when it closes (the server then quits within seconds). */
  function installLauncherHeartbeat() {
    const L = G.BM_LAUNCH, D = G.document;
    if (!L || !L.ping || !D || typeof D.addEventListener !== 'function') return;
    const send = state => {
      const url = L.ping + (String(L.ping).indexOf('?') < 0 ? '?' : '&') + 'state=' + state;
      try { if (G.navigator && G.navigator.sendBeacon && G.navigator.sendBeacon(url, '')) return; } catch (e) { /* fall through */ }
      try { G.fetch(url, { method: 'POST', keepalive: true }).catch(() => {}); } catch (e) { /* ignore */ }
    };
    let hideTimer = 0;
    D.addEventListener('visibilitychange', () => {
      clearTimeout(hideTimer);
      // delayed so that closing the window (hidden, then pagehide) reports only 'closed'
      if (D.visibilityState === 'hidden') hideTimer = setTimeout(() => send('hidden'), 400);
      else send('visible');
    });
    G.addEventListener('pagehide', () => { clearTimeout(hideTimer); send('closed'); });
    G.addEventListener('pageshow', e => { if (e && e.persisted) send('visible'); });
  }

  const native = {
    kind, version: VERSION,
    save, copyText, log, obabel, onOpenFile,
    _resolve, _receiveFile, _receiveLargeFile, _parseFormats: parseFormats,
    lastSave: null
  };
  BM.native = native;
  if (kind === 'linux-app') qtConnect();
  if (kind === 'macos-app' || kind === 'linux-app') { installAppShims(); announceReady(); }
  else if (kind === 'launcher') installLauncherHeartbeat();
})(typeof window !== 'undefined' ? (window.BM = window.BM || {}) : (globalThis.BM = globalThis.BM || {}));
