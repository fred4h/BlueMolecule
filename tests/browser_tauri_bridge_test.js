// The desktop app's side of BM.native (kind 'tauri-app'), in headless Chromium with a stand-in for Tauri's
// window.__TAURI__: documents from the command line and later launches, reading them in pieces, the Open and Save
// dialogs, Open Babel, the clipboard, web links, the updater dialog and the self-test script of the cloud builds
// (desktop/src-tauri/src/selftest.js). The Rust side is built and started on each system by .github/workflows/desktop.yml.
// Usage: node tests/browser_tauri_bridge_test.js
const { chromium } = require('playwright');
const path = require('path');
const fs = require('fs');

const ROOT = path.resolve(__dirname, '..');
const HTML = 'file://' + path.join(ROOT, 'dist', 'BlueMolecule.html');
const FIX = path.join(ROOT, 'tests', 'fixtures');
let passes = 0, failures = 0;
const check = (c, m) => { if (c) passes++; else { failures++; console.log('  FAIL:', m); } };
const sleep = ms => new Promise(r => setTimeout(r, ms));

/** window.__TAURI__ and window.BM_TAURI as the desktop app gives them, with its commands answered in the page. */
function mockTauri(files, opts) {
  return `(() => {
    const files = ${JSON.stringify(files)};
    const opts = ${JSON.stringify(opts || {})};
    const bytes = b64 => Uint8Array.from(atob(b64), c => c.charCodeAt(0));
    const store = {}; let next = 0;
    const doc = f => { const id = 'f' + (++next); store[id] = bytes(f.b64); return { name: f.name, size: store[id].length, fileId: id }; };
    const listeners = {};
    window.__calls = [];
    window.__tauriMock = { listeners, later: f => listeners['bm-open'] && listeners['bm-open']({ payload: doc(f) }) };
    window.BM_TAURI = Object.freeze({ version: '2.0.0', os: 'linux', memTotal: 16e9, memAvailable: 8e9, updates: !!opts.updates,
      selfTest: !!opts.selfTest, selfTestFiles: files.map(f => f.name) });
    window.__TAURI__ = {
      core: { invoke: async (cmd, args) => {
        window.__calls.push([cmd, args]);
        switch (cmd) {
          case 'bm_ready': return files.map(doc);
          case 'bm_read': { const b = store[args.fileId]; if (!b) throw 'This file is not open in BlueMolecule.'; return b.slice(args.offset, Math.min(b.length, args.offset + args.length)).buffer; }
          case 'bm_open_dialog': return opts.dialogCancel ? null : [doc(opts.dialogFile)];
          case 'bm_save': window.__saved = args; return { ok: true, cancelled: false, path: '/home/u/' + args.name, error: null };
          case 'bm_copy': window.__copied = args.text; return true;
          case 'bm_log': return null;
          case 'bm_open_url': window.__opened = args.url; return null;
          case 'bm_obabel': return args.request.op === 'formats' ? { ok: true, output: btoa('xyz -- XYZ cartesian coordinates format\\nmol -- MDL MOL format\\n'), error: '' } : { ok: true, output: btoa('converted'), error: '' };
          case 'bm_update_check': return opts.update || { supported: true, available: false, current: '2.0.0' };
          case 'bm_update_install': window.__installed = true; return null;
          case 'bm_self_test_done': window.__report = args; return null;
          default: throw 'unknown command ' + cmd;
        }
      } },
      event: { listen: async (name, cb) => { listeners[name] = cb; return () => {}; } }
    };
  })();`;
}

(async () => {
  const browser = await chromium.launch({ args: ['--enable-unsafe-swiftshader'] });
  const b64 = f => fs.readFileSync(path.join(FIX, f)).toString('base64');
  const files = [{ name: 'g16_water_optfreq.log', b64: b64('g16_water_optfreq.log') }, { name: 'cisplatin.xyz', b64: b64('cisplatin.xyz') }];
  const errors = [];
  const open = async (init, extra) => {
    const ctx = await browser.newContext({ viewport: { width: 1300, height: 900 } });
    await ctx.addInitScript(init);
    if (extra) await ctx.addInitScript(extra);
    const page = await ctx.newPage();
    page.on('pageerror', e => errors.push(e.message));
    await page.goto(HTML);
    await page.waitForFunction(() => window.BM && BM.app && BM.native);
    return page;
  };

  console.log('documents, reading, dialogs');
  let page = await open(mockTauri(files, { dialogFile: { name: 'water_traj.xyz', b64: b64('water_traj.xyz') }, updates: true,
    update: { supported: true, available: true, current: '2.0.0', version: '2.0.1', date: '2026-12-01T00:00:00Z', notes: 'Fixes <b>things</b>.\n- one\n- two' } }));
  check(await page.evaluate(() => BM.native.kind) === 'tauri-app', 'the bridge recognises the desktop app');
  await page.waitForFunction(() => BM.app.state.docs.length === 2 && BM.app.state.docs.every(d => !d.pending), null, { timeout: 15000 });
  let st = await page.evaluate(() => BM.app.state.docs.map(d => [d.name, d.data.frames.length, d.data.frames[0].Z.length]));
  check(JSON.stringify(st.map(x => x[0])) === '["g16_water_optfreq.log","cisplatin.xyz"]', 'documents of the command line open: ' + JSON.stringify(st));
  check(st[0][1] === 4 && st[1][2] === 11, 'and are read whole through bm_read');
  const reads = await page.evaluate(() => window.__calls.filter(c => c[0] === 'bm_read').length);
  check(reads >= 2, 'bytes come from bm_read (' + reads + ' calls)');
  check(await page.evaluate(() => window.__calls[0] && window.__calls[0][0] === 'bm_ready'), 'bm_ready is asked once the "bm-open" listener is in place');
  // a later launch with a document ("Open with" while running)
  await page.evaluate(f => window.__tauriMock.later(f), { name: 'g16_water_scan.log', b64: b64('g16_water_scan.log') });
  await page.waitForFunction(() => BM.app.state.docs.length === 3, null, { timeout: 10000 });
  check(await page.evaluate(() => BM.app.state.active.name) === 'g16_water_scan.log', 'a document sent later ("bm-open") opens too');
  // Open dialog
  await page.evaluate(() => BM.app.openDialog());
  await page.waitForFunction(() => BM.app.state.docs.length === 4, null, { timeout: 10000 });
  check(await page.evaluate(() => BM.app.state.active.name) === 'water_traj.xyz', 'the system Open dialog gives documents read the same way');
  // Save
  const saved = await page.evaluate(async () => { const ok = await BM.native.save('x.xyz', new Blob(['hello'])); return { ok, args: window.__saved, last: BM.native.lastSave }; });
  check(saved.ok === true && saved.args.name === 'x.xyz' && atob(saved.args.data) === 'hello' && saved.last.path === '/home/u/x.xyz', 'saving goes to bm_save with the bytes');
  // clipboard, links, Open Babel
  check(await page.evaluate(() => BM.native.copyText('abc').then(ok => ok && window.__copied === 'abc')), 'copying goes to the system clipboard');
  await page.evaluate(() => { const a = document.createElement('a'); a.href = 'https://example.org/docs'; a.textContent = 'x'; document.body.appendChild(a); a.click(); a.remove(); });
  check(await page.evaluate(() => window.__opened) === 'https://example.org/docs', 'web links open in the default browser');
  check(await page.evaluate(() => location.href.startsWith('file:')), 'and the app window stays on the app');
  check(await page.evaluate(() => BM.native.obabel.formats().then(f => f.length === 2 && f[0].code === 'xyz')), 'Open Babel formats through bm_obabel');
  check(await page.evaluate(() => BM.native.obabel.convert('x', 'xyz', 'mol').then(t => t === 'converted')), 'Open Babel conversion through bm_obabel');
  // memory for big files
  check(await page.evaluate(() => window.BM_TAURI.memAvailable === 8e9), 'the free memory comes from BM_TAURI');
  // updates
  await page.keyboard.press('Escape');
  await page.click('#btn-help');
  await page.waitForSelector('#dlg[open]');
  const helpButtons = await page.$$eval('#dlg .dlg-foot button', b => b.map(x => x.textContent));
  check(helpButtons.includes('Check for updates…'), 'Help offers Check for updates: ' + helpButtons.join(', '));
  await page.click('#dlg .dlg-foot button:has-text("Check for updates")');
  await page.waitForFunction(() => /2\.0\.1 is available/.test(document.querySelector('#dlg .dlg-head').textContent));
  const upd = await page.$eval('#dlg', e => ({ text: e.textContent.replace(/\s+/g, ' '), notesHtml: e.querySelector('.upd-notes').innerHTML }));
  check(/You have version 2\.0\.0; this one was released on 2026-12-01/.test(upd.text) && /Check for updates when BlueMolecule starts/.test(upd.text), 'the update dialog: versions, date, the opt-in box');
  check(upd.notesHtml.includes('&lt;b&gt;things&lt;/b&gt;') && !upd.notesHtml.includes('<b>'), 'release notes are shown as text, never as markup');
  await page.check('#dlg [data-upd="auto"]');
  await page.click('#dlg .dlg-foot button:has-text("Install and restart")');
  await sleep(150);
  check(await page.evaluate(() => window.__installed === true && localStorage.getItem('bluemolecule.updates.auto') === '1'), 'Install and restart installs; the start-up check is now on');
  await page.context().close();

  console.log('start-up check and cancelled dialog');
  page = await open(mockTauri([], { dialogCancel: true, updates: true, update: { supported: true, available: true, current: '2.0.0', version: '2.1.0' } }),
    `try { localStorage.setItem('bluemolecule.updates.auto', '1'); } catch (e) {}`);
  const nDocs = await page.evaluate(() => BM.app.state.docs.length);
  await page.evaluate(() => BM.app.openDialog());
  await sleep(200);
  check(await page.evaluate(n => BM.app.state.docs.length === n, nDocs), 'a cancelled Open dialog opens nothing');
  await page.waitForFunction(() => /2\.1\.0 is available/.test(document.querySelector('#toast').textContent), null, { timeout: 15000 });
  check(true, 'with the box ticked, start-up tells about a new version');
  await page.context().close();
  page = await open(mockTauri([], { updates: true, update: { supported: true, available: true, current: '2.0.0', version: '2.1.0' } }));
  await sleep(9000);
  check(!(await page.evaluate(() => window.__calls.some(c => c[0] === 'bm_update_check'))), 'without it, BlueMolecule never checks by itself');
  await page.context().close();

  console.log('self-test script');
  const selfTest = fs.readFileSync(path.join(ROOT, 'desktop', 'src-tauri', 'src', 'selftest.js'), 'utf8');
  page = await open(mockTauri(files, { selfTest: true }), selfTest);
  await page.waitForFunction(() => window.__report, null, { timeout: 60000 });
  const rep = await page.evaluate(() => window.__report);
  check(rep.report.ok === true, 'the self-test passes in Chromium: ' + JSON.stringify(rep.report.checks) + ' ' + (rep.report.error || '') + ' ' + rep.report.errors.join(' | '));
  check(rep.report.files.length === 2 && rep.report.files[1].atoms === 11 && rep.report.gpu, 'it reports the files and the GPU: ' + JSON.stringify(rep.report.files.map(f => f.name)) + ' ' + rep.report.gpu);
  check(typeof rep.image === 'string' && rep.image.length > 4000 && atob(rep.image.slice(0, 12)).startsWith('\x89PNG'), 'and a picture of the view');
  await page.context().close();

  const errs = errors.filter(e => !/WebGL|GPU stall/i.test(e));
  check(errs.length === 0, 'no page errors: ' + errs.join(' | '));
  await browser.close();
  console.log(`Desktop bridge tests: ${passes} passed, ${failures} failed`);
  process.exit(failures ? 1 : 0);
})().catch(e => { console.error(e); process.exit(1); });
