#!/usr/bin/env node
/* The Windows launcher's page contract (native/windows/BlueMolecule.cs): files opened with "Open with" are
 * embedded as data: URLs in window.BM_LAUNCH, written into a copy of BlueMolecule.html just before </head>,
 * with a UTF-8 byte-order mark. This test writes such a page exactly as the C# code does (mirrored below; the
 * literal pieces are checked against the C# source) and opens it from file:// in Chromium, like Edge does.
 * Run after the build:   node tests/browser_windows_launch_test.js */
'use strict';
const { chromium } = require('playwright');
const path = require('path');
const fs = require('fs');
const os = require('os');
const zlib = require('zlib');

const ROOT = path.resolve(__dirname, '..');
const PAGE = path.join(ROOT, 'dist', 'BlueMolecule.html');
const CS = path.join(ROOT, 'native', 'windows', 'BlueMolecule.cs');
const FIX = path.join(ROOT, 'tests', 'fixtures');

let passes = 0, failures = 0;
const check = (cond, msg) => { if (cond) passes++; else { failures++; console.log('  FAIL:', msg); } };

// ---- the C# code's pieces, mirrored; the source must contain each literal
const PIECES = {
  open: '<script>window.BM_LAUNCH={"files":[',
  name: '{"name":',
  url: ',"url":"data:application/octet-stream;base64,',
  close: '"}',
  version: '],"version":',
  end: '};</script>',
};
const csEsc = s => s.replace(/\\/g, '\\\\').replace(/"/g, '\\"');
/** JsonString() of the C# code. */
function jsonString(s) {
  let out = '"';
  for (const ch of s) {
    const c = ch.codePointAt(0);
    if (ch === '"') out += '\\"';
    else if (ch === '\\') out += '\\\\';
    else if (c < 0x20 || ch === '<' || ch === '>' || ch === '&' || c === 0x2028 || c === 0x2029) out += '\\u' + c.toString(16).padStart(4, '0');
    else out += ch;
  }
  return out + '"';
}
/** WriteLaunchPage() of the C# code. */
function launchPage(html, files, version) {
  const at = html.indexOf('</head>');
  let s = html.slice(0, at) + PIECES.open;
  files.forEach((f, i) => { if (i) s += ','; s += PIECES.name + jsonString(f.name) + PIECES.url + f.bytes.toString('base64') + PIECES.close; });
  s += PIECES.version + jsonString(version) + PIECES.end + html.slice(at);
  return Buffer.concat([Buffer.from([0xef, 0xbb, 0xbf]), Buffer.from(s, 'utf8')]);
}

(async () => {
  if (!fs.existsSync(PAGE)) { console.log('dist/BlueMolecule.html not found: run python3 tools/build.py first'); process.exit(2); }
  const cs = fs.readFileSync(CS, 'utf8');
  for (const [k, v] of Object.entries(PIECES)) check(cs.includes('"' + csEsc(v) + '"') || cs.includes('w.Write("' + csEsc(v)), 'C# source writes the "' + k + '" piece ' + JSON.stringify(v));
  check(/IndexOf\("<\/head>", StringComparison\.Ordinal\)/.test(cs), 'C# inserts before the first </head>');
  check(/new UTF8Encoding\(true\)/.test(cs), 'C# writes a UTF-8 byte-order mark');
  check(/c == '<' \|\| c == '>' \|\| c == '&'/.test(cs), 'C# escapes < > & in JSON strings');
  const version = (cs.match(/const string Version = "([^"]+)"/) || [])[1];
  check(version === '@VERSION@', 'C# version placeholder (filled in by the build): ' + version);

  const html = fs.readFileSync(PAGE, 'utf8');
  check(html.indexOf('</head>') > html.indexOf('<meta charset="utf-8">'), 'the charset meta tag comes before </head>');
  const tricky = 'wa<t>er & "opt" \\ ñ — test.log';
  const xyz = '3\nwater\nO 0 0 0.117\nH 0 0.757 -0.469\nH 0 -0.757 -0.469\n';
  const files = [
    { name: 'water.inpcrd', bytes: fs.readFileSync(path.join(FIX, 'struct', 'water.inpcrd')) },   // coordinates before their topology:
    { name: 'water.prmtop', bytes: fs.readFileSync(path.join(FIX, 'struct', 'water.prmtop')) },   // the page opens the topology first
    { name: tricky, bytes: fs.readFileSync(path.join(FIX, 'g16_water_optfreq.log')) },
    { name: 'small.xyz.gz', bytes: zlib.gzipSync(Buffer.from(xyz)) },
  ];
  const dir = fs.mkdtempSync(path.join(os.tmpdir(), 'bm-winlaunch-'));
  const out = path.join(dir, 'launch-test.html');
  fs.writeFileSync(out, launchPage(html, files, '1.4.1-test'));

  const browser = await chromium.launch({ args: ['--enable-unsafe-swiftshader'] });
  const page = await browser.newPage({ viewport: { width: 1280, height: 800 } });
  const errors = [];
  page.on('pageerror', e => errors.push(String(e)));
  page.on('console', m => { if (m.type() === 'error') errors.push('console: ' + m.text()); });
  await page.goto('file://' + out);
  await page.waitForFunction(() => window.BM && window.BM.app && window.BM.app.state.docs.length >= 3, null, { timeout: 20000 }).catch(() => {});
  await page.waitForTimeout(500);
  const st = await page.evaluate(() => {
    const a = window.BM.app, docs = a.state.docs;
    return {
      kind: window.BM.native.kind, charset: document.characterSet, launchVersion: window.BM_LAUNCH && window.BM_LAUNCH.version,
      docs: docs.map(d => ({ name: d.name, frames: d.data.frames.length, atoms: d.data.frames.length ? d.data.frames[0].Z.length : 0, vib: !!d.data.vibrations })),
      angstrom: document.querySelector('#tg-bl').textContent, injected: document.querySelectorAll('t, [data-pwn]').length,
      toast: document.querySelector('#toast').textContent
    };
  });
  check(st.kind === 'launcher', 'the page sees a launcher host: ' + st.kind);
  check(st.charset === 'UTF-8', 'the page is decoded as UTF-8: ' + st.charset);
  check(st.angstrom === 'Å', 'non-ASCII interface text intact: ' + JSON.stringify(st.angstrom));
  check(st.launchVersion === '1.4.1-test', 'BM_LAUNCH parsed: version ' + st.launchVersion);
  const byName = n => st.docs.find(d => d.name === n);
  check(byName(tricky) && byName(tricky).frames === 4 && byName(tricky).vib, 'Gaussian log with a tricky name opened intact: ' + JSON.stringify(st.docs.map(d => d.name)));
  check(byName('small.xyz') && byName('small.xyz').atoms === 3, 'gzip-compressed file opened (named without .gz)');
  const top = st.docs.find(d => /water\.prmtop/.test(d.name));
  check(top && top.frames === 1 && top.atoms === 4 && /water\.inpcrd: 1 geometry loaded into water\.prmtop/.test(st.toast),
    'topology opened first and the coordinates attached to it: ' + JSON.stringify(top) + ' ' + st.toast);
  check(st.docs.length === 3, 'three documents (log, xyz, prmtop with its coordinates): ' + st.docs.length);
  check(st.injected === 0, 'the file name did not become markup');
  // Open Babel is not provided by this host: asking must fail quietly
  const ob = await page.evaluate(() => window.BM.native.obabel.available());
  check(ob === false, 'Open Babel reported as unavailable');
  await page.evaluate(() => { try { window.BM.app.saveAs(); } catch (e) { /* ignore */ } });
  await page.waitForTimeout(400);
  const saveAs = await page.evaluate(() => (document.querySelector('[data-sa="ob"]') || {}).textContent || '');
  check(!/checking/.test(saveAs), 'Save as does not wait forever for Open Babel: ' + JSON.stringify(saveAs));
  check(errors.length === 0, 'no page errors: ' + errors.join(' | '));
  await browser.close();
  fs.rmSync(dir, { recursive: true, force: true });
  console.log(`Windows launch-page tests: ${passes} passed, ${failures} failed`);
  process.exit(failures ? 1 : 0);
})().catch(e => { console.error(e); process.exit(1); });
