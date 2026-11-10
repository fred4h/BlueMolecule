// Unit tests for src/native-bridge.js (BM.native) against mocked hosts: plain browser, macOS app
// (WKWebView message handlers) and the Python launcher (fetch endpoints).
// Usage: node tests/test_native_bridge.js
const fs = require('fs');
const path = require('path');
const vm = require('vm');

const SRC = fs.readFileSync(path.join(__dirname, '..', 'src', 'native-bridge.js'), 'utf8');
let passes = 0, failures = 0;
const check = (cond, msg) => { if (cond) passes++; else { failures++; console.log('  FAIL:', msg); } };
const sleep = ms => new Promise(r => setTimeout(r, ms));
const b64 = s => Buffer.from(s).toString('base64');
const unb64 = s => Buffer.from(s, 'base64');
const rejects = async (p, re) => { try { await p; return false; } catch (e) { return re ? re.test(e.message) : true; } };

/** Run native-bridge.js in a fresh global with the given extra globals; returns that global. */
function load(extra) {
  const g = { console, setTimeout, clearTimeout, atob, btoa, TextEncoder, TextDecoder, Blob, File, URL };
  Object.assign(g, extra || {});
  g.window = g;
  vm.createContext(g);
  vm.runInContext(SRC, g, { filename: 'native-bridge.js' });
  return g;
}

/** A mocked WKWebView host: handler name -> function(msg, global) that may reply via _resolve. */
function macHost(reply) {
  const posted = [];
  const messageHandlers = {};
  const ref = {};
  for (const name of ['bmReady', 'bmSave', 'bmObabel', 'bmCopy', 'bmLog']) {
    messageHandlers[name] = {
      postMessage(msg) {
        posted.push({ name, msg: JSON.parse(JSON.stringify(msg)) });
        if (reply[name]) setTimeout(() => reply[name](msg, ref.g), 1);
      }
    };
  }
  return { posted, ref, globals: { webkit: { messageHandlers } } };
}

const FORMATS_TEXT = 'acesin -- ACES input format [Write-only]\nmol2 -- Sybyl Mol2 format\n' +
  'xyz -- XYZ cartesian coordinates format\nlog -- Generic output file format [Read-only]\n\nnot a format line\n';

(async () => {
  // ---------- plain browser ----------
  console.log('Plain browser');
  {
    const g = load();
    const N = g.BM.native;
    check(N.kind === 'browser', 'kind browser: ' + N.kind);
    check(await N.save('a.xyz', new Blob(['x'])) === null, 'save not handled in a browser');
    check(await N.copyText('x') === null, 'copyText not handled in a browser');
    check(await N.obabel.available() === false, 'obabel unavailable in a browser');
    check(await rejects(N.obabel.formats(), /app or the bluemolecule launcher/), 'formats rejects in a browser');
    check(await rejects(N.obabel.convert('x', 'xyz', 'sdf')), 'convert rejects in a browser');
    check(N._resolve('nope', {}) === false, '_resolve ignores unknown ids');
    check(typeof N.onOpenFile(() => {}) === 'function', 'onOpenFile returns an unsubscribe function');
  }

  // ---------- macOS app ----------
  console.log('macOS app');
  {
    let saveReply = { ok: true, path: '/Users/me/water.xyz' };
    const host = macHost({
      bmSave: (m, g) => g.BM.native._resolve(m.id, saveReply),
      bmObabel: (m, g) => {
        if (m.op === 'formats') g.BM.native._resolve(m.id, { ok: true, output: b64(FORMATS_TEXT), error: '' });
        else if (m.inFormat === 'bad') g.BM.native._resolve(m.id, { ok: false, output: '', error: 'boom: 0 molecules converted' });
        else g.BM.native._resolve(m.id, { ok: true, output: b64('converted ' + m.inFormat + '->' + m.outFormat + ' ' + unb64(m.input).toString() + ' ' + m.options.join(',')), error: '' });
      }
    });
    const g = load(host.globals);
    host.ref.g = g;
    const N = g.BM.native;
    check(N.kind === 'macos-app', 'kind macos-app: ' + N.kind);

    // save: text blob, cancel, binary bytes
    const text = '3\nwater\nO 0 0 0\nH 0 0 1\nH 0 1 0\n';
    check(await N.save('water.xyz', new Blob([text], { type: 'chemical/x-xyz' })) === true, 'save resolves true');
    const s1 = host.posted.filter(p => p.name === 'bmSave').pop().msg;
    check(s1.name === 'water.xyz' && s1.mime === 'chemical/x-xyz' && unb64(s1.base64).toString() === text && /^bm\d+$/.test(s1.id), 'bmSave payload ' + JSON.stringify(s1).slice(0, 120));
    check(N.lastSave && N.lastSave.path === '/Users/me/water.xyz', 'lastSave records the path');
    saveReply = { ok: false };
    check(await N.save('x.png', new Blob(['x'])) === false, 'cancelled save resolves false');
    saveReply = { ok: true, path: '/tmp/b.bin' };
    const bytes = new Uint8Array(256).map((_, i) => i);
    check(await N.save('b.bin', bytes) === true, 'save accepts a Uint8Array');
    const s3 = host.posted.filter(p => p.name === 'bmSave').pop().msg;
    check(Buffer.compare(unb64(s3.base64), Buffer.from(bytes)) === 0 && s3.mime === 'application/octet-stream', 'binary bytes round-trip through base64');
    check(await N.save('s.txt', 'héllo') === true && unb64(host.posted.filter(p => p.name === 'bmSave').pop().msg.base64).toString() === 'héllo', 'save accepts a string (UTF-8)');

    // copy
    check(await N.copyText('C1=CC=CC=C1') === true, 'copyText resolves true');
    const c = host.posted.filter(p => p.name === 'bmCopy').pop();
    check(c && c.msg.text === 'C1=CC=CC=C1' && !('id' in c.msg), 'bmCopy payload {text}');

    // obabel
    const nObabelBefore = host.posted.filter(p => p.name === 'bmObabel').length;
    const f = await N.obabel.formats();
    check(f.length === 4, 'four formats parsed: ' + JSON.stringify(f));
    const byCode = Object.fromEntries(f.map(x => [x.code, x]));
    check(byCode.acesin && byCode.acesin.read === false && byCode.acesin.write === true && byCode.acesin.description === 'ACES input format', 'write-only format parsed');
    check(byCode.log && byCode.log.read === true && byCode.log.write === false, 'read-only format parsed');
    check(byCode.xyz && byCode.xyz.read && byCode.xyz.write && byCode.xyz.description === 'XYZ cartesian coordinates format', 'read/write format parsed');
    check(await N.obabel.available() === true, 'obabel available');
    await N.obabel.formats();
    check(host.posted.filter(p => p.name === 'bmObabel').length === nObabelBefore + 1, 'formats are cached after the first success');
    const fm = host.posted.find(p => p.name === 'bmObabel').msg;
    check(fm.op === 'formats', 'formats request op');
    const out = await N.obabel.convert('CCO', 'smi', 'xyz', ['--gen3d', '-h']);
    check(out === 'converted smi->xyz CCO --gen3d,-h', 'convert result: ' + out);
    const cm = host.posted.filter(p => p.name === 'bmObabel').pop().msg;
    check(cm.op === 'convert' && cm.inFormat === 'smi' && cm.outFormat === 'xyz' && unb64(cm.input).toString() === 'CCO' && cm.options.length === 2, 'convert request payload');
    check(await rejects(N.obabel.convert('x', 'bad', 'xyz'), /0 molecules converted/), 'convert error message is passed on');
    const nPosted = host.posted.length;
    check(await rejects(N.obabel.convert('x', 'x y', 'xyz'), /Invalid Open Babel format/), 'bad format code rejected');
    check(host.posted.length === nPosted, 'bad format code is not sent to the app');

    // documents from Finder/Dock: queued until someone listens, then in order
    check(N._receiveFile('first.xyz', b64('3\nfirst\n')) === true, '_receiveFile accepts');
    N._receiveFile('second.cube', b64('cube'));
    const got = [];
    const off = N.onOpenFile((name, data) => got.push([name, Buffer.from(data).toString(), Object.prototype.toString.call(data)]));
    check(got.length === 2 && got[0][0] === 'first.xyz' && got[0][1] === '3\nfirst\n' && got[1][0] === 'second.cube', 'queued files delivered in order: ' + JSON.stringify(got));
    check(got[0][2] === '[object Uint8Array]', 'callback receives a Uint8Array');
    N._receiveFile('third.pdb', b64('ATOM'));
    check(got.length === 3 && got[2][0] === 'third.pdb', 'later files delivered at once');
    off();
    N._receiveFile('fourth.xyz', b64('x'));
    check(got.length === 3, 'unsubscribed callback no longer called');

    // bmReady waits for BM.app
    await sleep(120);
    check(!host.posted.some(p => p.name === 'bmReady'), 'bmReady not posted before BM.app exists');
    const opened = [];
    g.BM.app = { openFile: file => { opened.push(file.name); return Promise.resolve(); } };
    await sleep(150);
    const ready = host.posted.filter(p => p.name === 'bmReady');
    check(ready.length === 1 && ready[0].msg.app === true && ready[0].msg.version === '1.2.0', 'bmReady posted once BM.app exists: ' + JSON.stringify(ready));
    check(opened.join() === 'fourth.xyz', 'with no listener, queued files go to BM.app.openFile: ' + opened);
    N._receiveFile('fifth.mol', b64('mol'));
    check(opened.join() === 'fourth.xyz,fifth.mol', 'BM.app.openFile gets later files too');

    // a reply that never comes
    check(N._resolve('bm9999', { ok: true }) === false, 'unknown reply id ignored');
  }

  // ---------- macOS safety nets: <a download> and navigator.clipboard ----------
  console.log('macOS shims');
  {
    const host = macHost({ bmSave: (m, g) => g.BM.native._resolve(m.id, { ok: true, path: '/x/' + m.name }) });
    let clickedThrough = 0;
    class FakeAnchor {
      constructor(href, download) { this.href = href; this.attrs = download ? { download } : {}; }
      hasAttribute(n) { return n in this.attrs; }
      getAttribute(n) { return this.attrs[n]; }
    }
    FakeAnchor.prototype.click = function () { clickedThrough++; };
    let created = 0;
    const FakeURL = { createObjectURL: () => 'blob:fake/' + (++created), revokeObjectURL: () => {} };
    const nav = { clipboard: { writeText: () => Promise.reject(new Error('NotAllowedError')) } };
    const g = load(Object.assign({}, host.globals, { HTMLAnchorElement: FakeAnchor, URL: FakeURL, navigator: nav, fetch: async () => { throw new Error('fetch should not be needed'); } }));
    host.ref.g = g;
    const blob = new Blob(['spectrum'], { type: 'image/svg+xml' });
    const url = g.URL.createObjectURL(blob);
    new FakeAnchor(url, 'ir.svg').click();
    await sleep(30);
    const s = host.posted.filter(p => p.name === 'bmSave').pop();
    check(s && s.msg.name === 'ir.svg' && unb64(s.msg.base64).toString() === 'spectrum', '<a download> click goes to the native Save dialog');
    new FakeAnchor('https://example.org/', null).click();
    check(clickedThrough === 1, 'ordinary links still click through');
    await nav.clipboard.writeText('copied text');
    const c = host.posted.filter(p => p.name === 'bmCopy').pop();
    check(c && c.msg.text === 'copied text', 'navigator.clipboard.writeText goes to the pasteboard');
  }

  // ---------- launcher ----------
  console.log('Launcher');
  {
    const requests = [];
    let status = 200;
    const fetch = async (url, init) => {
      requests.push({ url, init });
      const body = init && init.body ? JSON.parse(init.body) : {};
      let res;
      if (body.op === 'formats') res = { ok: true, output: b64(FORMATS_TEXT), error: '' };
      else if (body.op === 'convert') res = body.inFormat === 'bad' ? { ok: false, error: 'not found' } : { ok: true, output: b64('OUT:' + unb64(body.input).toString()), error: '' };
      else res = {};
      return { ok: status === 200, status, json: async () => res };
    };
    const docL = {}, winL = {};
    const document = { visibilityState: 'visible', addEventListener: (t, f) => { docL[t] = f; } };
    const beacons = [];
    const navigator = { sendBeacon: (u, b) => { beacons.push(u); return true; } };
    const BM_LAUNCH = { files: [], ping: '/TOK/ping', obabel: '/TOK/obabel', version: '1.2.0' };
    const g = load({ BM_LAUNCH, fetch, document, navigator, addEventListener: (t, f) => { winL[t] = f; } });
    const N = g.BM.native;
    check(N.kind === 'launcher', 'kind launcher: ' + N.kind);
    check(await N.save('a.xyz', new Blob(['x'])) === null, 'launcher save not handled (the browser has its own Save dialog)');
    check(await N.copyText('x') === null, 'launcher copy not handled');
    const f = await N.obabel.formats();
    check(f.length === 4 && requests[0].url === '/TOK/obabel' && requests[0].init.method === 'POST', 'formats via POST /TOKEN/obabel');
    check(JSON.parse(requests[0].init.body).op === 'formats' && requests[0].init.headers['Content-Type'] === 'application/json', 'formats request body');
    check(await N.obabel.convert(new Uint8Array([65, 66]), 'xyz', 'pdb', []) === 'OUT:AB', 'convert via the launcher');
    const body = JSON.parse(requests[requests.length - 1].init.body);
    check(body.op === 'convert' && body.inFormat === 'xyz' && body.outFormat === 'pdb' && Array.isArray(body.options), 'convert request body');
    check(await rejects(N.obabel.convert('x', 'bad', 'pdb'), /not found/), 'launcher error passed on');
    status = 500;
    check(await rejects(N.obabel.convert('x', 'xyz', 'pdb'), /HTTP 500/), 'HTTP error rejects');

    // heartbeat states
    check(typeof docL.visibilitychange === 'function' && typeof winL.pagehide === 'function', 'visibility/pagehide listeners installed');
    document.visibilityState = 'hidden'; docL.visibilitychange();
    await sleep(500);
    check(beacons.pop() === '/TOK/ping?state=hidden', 'hidden window reported');
    document.visibilityState = 'visible'; docL.visibilitychange();
    check(beacons.pop() === '/TOK/ping?state=visible', 'visible window reported');
    document.visibilityState = 'hidden'; docL.visibilitychange(); winL.pagehide({ persisted: false });
    await sleep(500);
    check(beacons.length === 1 && beacons[0] === '/TOK/ping?state=closed', 'closing reports only "closed": ' + beacons);
  }

  // ---------- Linux app (Qt WebEngine window, QWebChannel) ----------
  console.log('Linux app');
  {
    // a mocked channel object with the contract of native/linux/bluemolecule_qt.py (class Bridge)
    const FILE = Buffer.alloc(13 * 1024 * 1024 + 4321);
    for (let i = 0; i < FILE.length; i++) FILE[i] = (i * 7 + (i >> 13)) & 255;
    const LIMIT = 6291456;
    const posted = [], ref = {};
    let saveReply = { ok: true, path: '/home/u/water.xyz' };
    const reply = {
      bmSave: m => ref.g.BM.native._resolve(m.id, saveReply),
      bmObabel: m => ref.g.BM.native._resolve(m.id, m.op === 'formats' ? { ok: true, output: b64(FORMATS_TEXT), error: '' }
        : { ok: true, output: b64('OB ' + m.inFormat + '>' + m.outFormat + ' ' + unb64(m.input).toString()), error: '' }),
      bmRead: m => {
        if (m.fileId !== 'f1') return ref.g.BM.native._resolve(m.id, { ok: false, error: 'Unknown file.' });
        const data = FILE.subarray(m.offset, m.offset + Math.min(m.length, LIMIT));
        ref.g.BM.native._resolve(m.id, { ok: true, offset: m.offset, length: data.length, eof: m.offset + data.length >= FILE.length, base64: data.toString('base64') });
      }
    };
    const bmHost = { post(name, text) { const msg = JSON.parse(text); posted.push({ name, msg }); if (reply[name]) setTimeout(() => reply[name](msg), 1); } };
    let channels = 0;
    function QWebChannel(transport, cb) { channels++; setTimeout(() => cb({ objects: { bmHost } }), 30); }
    const clicked = [];
    class HTMLAnchorElement {
      constructor(href, download) { this.href = href; this.dl = download; }
      hasAttribute(n) { return n === 'download' && this.dl != null; }
      getAttribute(n) { return n === 'download' ? this.dl : null; }
      click() { clicked.push(this.href); }
    }
    const navigator = { clipboard: { writeText: async () => 'browser clipboard' } };
    const g = load({ BM_QT: { shell: 'linux-app', version: '1.2.0', readLimit: LIMIT }, qt: { webChannelTransport: {} }, QWebChannel,
                     HTMLAnchorElement, navigator, fetch: async () => ({ blob: async () => new Blob(['fetched']) }) });
    ref.g = g;
    const N = g.BM.native;
    check(N.kind === 'linux-app', 'kind linux-app: ' + N.kind);
    check(channels === 1, 'the channel is opened at start-up');

    // requests made before the channel is open wait, then go out in order
    const early = N.save('water.xyz', new Blob(['3\nw\n'], { type: 'chemical/x-xyz' }));
    await sleep(5);
    check(posted.length === 0, 'nothing is posted before the channel is open');
    check(typeof N.copyText('first').then === 'function', 'copyText returns a promise');
    check(await early === true && N.lastSave.path === '/home/u/water.xyz', 'save through the Qt Save dialog');
    check(posted[0].name === 'bmSave' && posted[1].name === 'bmCopy' && posted[1].msg.text === 'first', 'queued messages keep their order');
    const s = posted[0].msg;
    check(/^bm\d+$/.test(s.id) && s.name === 'water.xyz' && s.mime === 'chemical/x-xyz' && unb64(s.base64).toString() === '3\nw\n', 'bmSave payload ' + JSON.stringify(s));
    saveReply = { ok: false, cancelled: true };
    check(await N.save('x.png', new Uint8Array([1, 2])) === false, 'cancelled save resolves false');
    check(await N.copyText('C1=CC=CC=C1') === true && posted[posted.length - 1].msg.text === 'C1=CC=CC=C1' && !('id' in posted[posted.length - 1].msg), 'bmCopy {text}');
    check(N.log('error', 'boom') === true && posted[posted.length - 1].name === 'bmLog' && posted[posted.length - 1].msg.level === 'error'
          && posted[posted.length - 1].msg.message === 'boom', 'bmLog {level, message}');
    const f = await N.obabel.formats();
    check(f.length === 4 && posted.some(p => p.name === 'bmObabel' && p.msg.op === 'formats'), 'Open Babel formats through bmObabel');
    check(await N.obabel.convert('CCO', 'smi', 'xyz', ['--gen3d']) === 'OB smi>xyz CCO', 'Open Babel convert through bmObabel');

    // app shims: <a download> goes to the Save dialog, navigator.clipboard.writeText to bmCopy
    saveReply = { ok: true, path: '/home/u/pic.png' };
    const url = g.URL.createObjectURL(new Blob(['PNGDATA'], { type: 'image/png' }));
    new HTMLAnchorElement(url, 'pic.png').click();
    await sleep(20);
    const sv = posted.filter(p => p.name === 'bmSave').pop().msg;
    check(sv.name === 'pic.png' && unb64(sv.base64).toString() === 'PNGDATA' && sv.mime === 'image/png' && clicked.length === 0, '<a download> -> bmSave');
    new HTMLAnchorElement('https://example.org/', null).click();
    check(clicked.length === 1, 'other links click normally');
    await navigator.clipboard.writeText('from clipboard API');
    check(posted[posted.length - 1].name === 'bmCopy' && posted[posted.length - 1].msg.text === 'from clipboard API', 'navigator.clipboard.writeText -> bmCopy');

    // the page announces itself once the viewer exists
    g.BM.app = {};
    await sleep(120);
    const ready = posted.filter(p => p.name === 'bmReady');
    check(ready.length === 1 && ready[0].msg.version === '1.2.0' && ready[0].msg.app === true, 'bmReady {version, app}');

    // documents from the command line: small ones as bytes, large ones as a handle that reads on demand
    const got = [];
    N.onOpenFile((a, b) => got.push([a, b]));
    N._receiveFile('small.xyz', b64('1\n\nH 0 0 0\n'));
    check(got.length === 1 && got[0][0] === 'small.xyz' && Buffer.from(got[0][1]).toString() === '1\n\nH 0 0 0\n', 'small file -> (name, Uint8Array)');
    check(N._receiveLargeFile('big.xyz', FILE.length, 'f1') === true, '_receiveLargeFile accepts a descriptor');
    const h = got[1] && got[1][0];
    check(got.length === 2 && got[1][1] === undefined && h && h.name === 'big.xyz' && h.size === FILE.length && typeof h.readChunk === 'function'
          && Object.isFrozen(h), 'large file -> one {name, size, readChunk} argument');
    let reads = posted.filter(p => p.name === 'bmRead').length;
    const head = await h.readChunk(0, 100);
    check(Object.prototype.toString.call(head) === '[object Uint8Array]' && Buffer.compare(Buffer.from(head), FILE.subarray(0, 100)) === 0, 'readChunk(0, 100)');
    const span = await h.readChunk(1000, 13 * 1024 * 1024);
    const spanReads = posted.filter(p => p.name === 'bmRead').slice(reads + 1);
    check(span.length === 13 * 1024 * 1024 && Buffer.compare(Buffer.from(span), FILE.subarray(1000, 1000 + 13 * 1024 * 1024)) === 0, 'readChunk across several replies');
    check(spanReads.length === 3 && spanReads.every(p => p.msg.length <= LIMIT && p.msg.fileId === 'f1') && spanReads[1].msg.offset === 1000 + LIMIT,
          'bmRead requests of at most 6 MiB (8 MiB of base64): ' + JSON.stringify(spanReads.map(p => [p.msg.offset, p.msg.length])));
    const tail = await h.readChunk(FILE.length - 10, 100);
    check(tail.length === 10 && Buffer.compare(Buffer.from(tail), FILE.subarray(FILE.length - 10)) === 0, 'readChunk at the end is cut to the file size');
    reads = posted.filter(p => p.name === 'bmRead').length;
    check((await h.readChunk(FILE.length, 10)).length === 0 && (await h.readChunk(5, 0)).length === 0
          && posted.filter(p => p.name === 'bmRead').length === reads, 'empty ranges need no request');
    check(await rejects(h.readChunk(-1, 10), /whole numbers/) && await rejects(h.readChunk(0, 1.5), /whole numbers/) && await rejects(h.readChunk('0', 1)),
          'bad ranges are refused');
    N._receiveLargeFile('gone.xyz', 50, 'f9');
    check(await rejects(got[2][0].readChunk(0, 10), /Unknown file/), 'a read error rejects with the host message');
    check(N._receiveLargeFile('x', -1, 'f1') === false && N._receiveLargeFile('x', 10, '') === false && got.length === 3, 'bad descriptors are refused');
  }
  {
    // no channel ever appears: requests are answered with an error, save/copy fall back to the page's own paths
    const g = load({ BM_QT: { shell: 'linux-app', connectTimeout: 60 } });
    const N = g.BM.native;
    const early = N.save('a.xyz', 'x');
    const r = await early;
    check(N.kind === 'linux-app' && r === false && /did not appear/.test(N.lastSave.error), 'a request made while waiting gets an error: ' + JSON.stringify(N.lastSave));
    check(await N.save('b.xyz', 'x') === null && await N.copyText('x') === null && N.log('info', 'x') === null, 'afterwards save/copy/log are not handled (null)');
    check(await rejects(N.obabel.formats(), /not available/), 'Open Babel reports that the window is not reachable');
  }

  console.log(`\nNative bridge tests: ${passes} passed, ${failures} failed`);
  process.exit(failures ? 1 : 0);
})().catch(e => { console.error(e); process.exit(2); });
