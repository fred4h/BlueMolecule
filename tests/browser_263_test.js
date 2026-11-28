// Browser checks of 2.6.3 in dist/BlueMolecule.html (headless Chromium):
//   - exports, calculation inputs and Save As use the stored geometry while a vibration is animated, never the
//     momentary picture of the animation (a random point of the oscillation)
//   - every side-panel tab stays visible (a second row instead of tabs scrolled out of sight)
//   - the smallest desktop window (820 x 560): the plot keeps its room, the start screen fits and scrolls
//   - plain messages for the updater's usual failures
// Usage: node tests/browser_263_test.js
const { chromium } = require('playwright');
const path = require('path');
const fs = require('fs');

const ROOT = path.resolve(__dirname, '..');
const HTML = 'file://' + path.join(ROOT, 'dist', 'BlueMolecule.html');
const FIX = path.join(ROOT, 'tests', 'fixtures');

let failures = 0, passes = 0;
const check = (cond, msg) => { if (cond) passes++; else { failures++; console.log('  FAIL:', msg); } };

async function newPage(browser, w, h, errors) {
  const page = await (await browser.newContext({ viewport: { width: w, height: h } })).newPage();
  page.on('pageerror', e => errors.push('pageerror: ' + e.message));
  page.on('console', m => { if (m.type() === 'error' && !/Failed to load resource/.test(m.text())) errors.push('console: ' + m.text()); });
  await page.goto(HTML);
  await page.waitForFunction(() => window.BM && window.BM.app);
  return page;
}
const open = async (page, f) => {
  const n = await page.evaluate(() => BM.app.state.docs.length);
  await page.setInputFiles('#file-input', path.join(FIX, f));
  await page.waitForFunction(n => BM.app.state.docs.length > n, n);
  await page.waitForTimeout(250);
};

(async () => {
  const browser = await chromium.launch({ args: ['--enable-unsafe-swiftshader'] });
  const errors = [];
  let page = await newPage(browser, 1440, 900, errors);
  await page.evaluate(() => { window.__saved = []; window.__bmSaveHook = async (name, blob) => { window.__saved.push({ name, text: await blob.text() }); }; });

  // ------------------------------------------------------------ exports while the imaginary mode of a TS is animated
  await open(page, 'g16_sn2_ts_optfreq.log');
  await page.waitForFunction(() => { const d = BM.app.state.active, f = d.data.frames[d.frameIndex].xyz, x = d._animXyz; if (BM.app.state.anim.kind !== 'vib' || !x) return false; for (let i = 0; i < f.length; i++) if (Math.abs(f[i] - x[i]) > 0.02) return true; return false; }, null, { timeout: 5000 });
  const r = await page.evaluate(async () => {
    const a = BM.app, d = a.state.active, f = d.data.frames[d.frameIndex].xyz;
    const diff = x => { let m = 0; for (let i = 0; i < f.length; i++) m = Math.max(m, Math.abs(f[i] - x[i])); return m; };
    const shown = diff(d._animXyz);
    window.__saved = [];
    a.exportSimple('xyz');
    await new Promise(res => setTimeout(res, 150));
    const quick = diff(BM.parseFile('q.xyz', window.__saved[0].text).frames[0].xyz);
    const ctx = diff(a.saveContext(d, false).xyz);
    return { anim: a.state.anim.kind, shown, quick, ctx };
  });
  check(r.anim === 'vib' && r.shown > 0.02, 'the imaginary mode is animated, and the picture is displaced: ' + JSON.stringify(r));
  check(r.quick < 1e-6, 'Export → XYZ writes the stored geometry, not the animated one: ' + r.quick);
  check(r.ctx < 1e-6, 'Save As uses the stored geometry: ' + r.ctx);
  // the calculation-input dialog: its preview holds the stored geometry
  await page.evaluate(() => BM.app.exportCalc('gaussian'));
  await page.waitForSelector('#dlg[open] textarea[data-c="preview"]');
  const inp = await page.evaluate(() => {
    const d = BM.app.state.active, f = d.data.frames[d.frameIndex].xyz, t = document.querySelector('#dlg textarea[data-c="preview"]').value;
    const x = BM.parseFile('in.gjf', t).frames[0].xyz;
    let m = 0; for (let i = 0; i < f.length; i++) m = Math.max(m, Math.abs(f[i] - x[i]));
    return m;
  });
  check(inp < 1e-5, 'the calculation input uses the stored geometry: ' + inp);
  await page.keyboard.press('Escape');

  // ------------------------------------------------------------ every tab visible
  await open(page, 'g16_water_optfreq.log');
  const tabs = await page.evaluate(() => {
    const box = document.querySelector('#tabs').getBoundingClientRect();
    return [...document.querySelectorAll('#tabs button')].map(b => { const r = b.getBoundingClientRect(); return { t: b.dataset.tab, inside: r.width > 0 && r.left >= box.left - 1 && r.right <= box.right + 1 && r.bottom <= box.bottom + 1 }; });
  });
  check(tabs.length === 6 && tabs.every(t => t.inside), 'all six tabs are inside the panel: ' + JSON.stringify(tabs));
  await page.click('#tabs button[data-tab="build"]');
  check(await page.evaluate(() => BM.app.state.active.tab === 'build'), 'the Build tab can be clicked');

  // ------------------------------------------------------------ updater messages
  const updater = async (error) => {
    await page.evaluate(err => { window.BM.native = { kind: 'test', updates: { supported: true, check: async () => ({ supported: true, available: false, current: '2.6.3', error: err }) } }; }, error);
    await page.click('#btn-help');
    await page.waitForSelector('#dlg[open]');
    await page.click('#dlg button:has-text("Check for updates")');
    await page.waitForFunction(() => document.querySelector('#dlg').open && document.querySelector('#dlg .dlg-head').textContent === 'Updates', null, { timeout: 5000 });
    const t = await page.evaluate(() => document.querySelector('#dlg').textContent);
    await page.keyboard.press('Escape');
    return t;
  };
  let t = await updater('Could not fetch a valid release JSON from the remote');
  check(/still being prepared/.test(t) && /The updater said: Could not fetch a valid release JSON/.test(t), 'missing latest.json: ' + t.slice(0, 300));
  t = await updater('the platform `windows-x86_64` was not found on the response `platforms` object');
  check(/no update for this system yet/.test(t), 'platform missing: ' + t.slice(0, 300));
  t = await updater('error sending request for url (https://github.com/fred4h/BlueMolecule/releases/latest/download/latest.json)');
  check(/could not reach GitHub/.test(t), 'offline: ' + t.slice(0, 300));
  await page.evaluate(() => { delete window.BM.native; });
  await page.context().close();

  // ------------------------------------------------------------ the smallest desktop window
  page = await newPage(browser, 820, 560, errors);
  const empty = await page.evaluate(() => {
    const e = document.querySelector('#empty'), box = e.getBoundingClientRect(), c = document.querySelector('.empty-card').getBoundingClientRect(), f = document.querySelector('#formats').getBoundingClientRect();
    return { cardInside: c.left >= box.left - 1 && c.right <= box.right + 1, formatsInside: f.left >= box.left - 1 && f.right <= box.right + 1, scrolls: e.scrollHeight > e.clientHeight && getComputedStyle(e).overflowY === 'auto' };
  });
  check(empty.cardInside && empty.formatsInside, 'start screen: the card and the list of formats fit the width: ' + JSON.stringify(empty));
  check(empty.scrolls, 'start screen: scrolls when the window is short: ' + JSON.stringify(empty));
  await open(page, 'g16_water_optfreq.log');
  await page.click('#tabs button[data-tab="vib"]');
  await page.waitForTimeout(200);
  const lay = await page.evaluate(() => {
    const h = s => document.querySelector(s).getBoundingClientRect().height;
    const c = document.querySelector('#plot-controls'), pb = document.querySelector('#playbar').getBoundingClientRect(), tb = document.querySelector('#toolbox').getBoundingClientRect();
    return { plot: Math.round(h('.plot-body')), view: Math.round(h('#viewport')), head: Math.round(h('.plot-head')), controls: Math.round(c.getBoundingClientRect().width), scrollable: getComputedStyle(c).overflowX,
      playbarClear: pb.width === 0 || pb.left >= tb.right - 1 || pb.top >= tb.bottom - 1 };
  });
  check(lay.plot >= 80 && lay.view >= 150, 'plot and viewer both keep room: ' + JSON.stringify(lay));
  check(lay.head <= 80 && lay.controls >= 400 && lay.scrollable === 'auto', 'narrow plot area: the controls get their own scrolling row: ' + JSON.stringify(lay));
  check(lay.playbarClear, 'the playbar does not cover the tool box: ' + JSON.stringify(lay));

  check(!errors.length, 'no page errors: ' + errors.join(' | '));
  await browser.close();
  console.log(`2.6.3 browser tests: ${passes} passed, ${failures} failed`);
  process.exit(failures ? 1 : 0);
})().catch(e => { console.error(e); process.exit(1); });
