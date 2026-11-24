// Browser checks of 2.6 in dist/BlueMolecule.html (headless Chromium): Download from PubChem and the PDB (#22), with
// every request to the network intercepted (nothing leaves this computer):
//   - no request at start-up, on opening the dialog, or while "Allow downloads" is unticked
//   - PubChem by name (3D conformer as SDF), PubChem without a 3D conformer (built from the SMILES), the PDB by ID
//     and its biological assembly; what is sent (GET, no cookies) and to which two hosts
//   - errors shown in the dialog, a bad ID sends nothing, Cancel stops a download, unticking turns downloads off
// Usage: node tests/browser_fetch_test.js
const { chromium } = require('playwright');
const path = require('path');
const fs = require('fs');

const ROOT = path.resolve(__dirname, '..');
const HTML = 'file://' + path.join(ROOT, 'dist', 'BlueMolecule.html');
const MMCIF = fs.readFileSync(path.join(ROOT, 'tests', 'fixtures', 'struct', 'small.mmcif'), 'utf8');
const PUG = 'https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/', RCSB = 'https://files.rcsb.org/download/';

let failures = 0, passes = 0;
const check = (cond, msg) => { if (cond) passes++; else { failures++; console.log('  FAIL:', msg); } };

(async () => {
  const browser = await chromium.launch({ args: ['--enable-unsafe-swiftshader'] });
  const context = await browser.newContext({ viewport: { width: 1440, height: 1000 } });
  const page = await context.newPage();
  const errors = [];
  page.on('pageerror', e => errors.push('pageerror: ' + e.message));
  page.on('console', m => { if (m.type() === 'error' && !/Failed to load resource/.test(m.text())) errors.push('console: ' + m.text()); });
  page.on('dialog', d => { errors.push('dialog: ' + d.message()); d.dismiss(); });

  // every http(s) request is recorded and answered here; only the two databases get an answer
  const requests = [];
  let aspirin = '', hold = null;
  const cors = { 'access-control-allow-origin': '*' };
  await context.route(/^https?:/, async route => {
    const req = route.request(), url = req.url();
    requests.push({ url, method: req.method(), headers: req.headers(), body: req.postData() });
    const send = (status, body, type) => route.fulfill({ status, body, headers: Object.assign({ 'content-type': type || 'text/plain' }, cors) });
    if (url === PUG + 'name/aspirin/SDF?record_type=3d') return send(200, aspirin, 'chemical/x-mdl-sdfile');
    if (url === PUG + 'name/sodium%20chloride/SDF?record_type=3d') return send(404, '{"Fault":{"Code":"PUGREST.NotFound","Message":"No records found"}}', 'application/json');
    if (url === PUG + 'name/sodium%20chloride/property/SMILES/JSON') return send(200, '{"PropertyTable":{"Properties":[{"CID":5234,"SMILES":"[Na+].[Cl-]"}]}}', 'application/json');
    if (url === RCSB + '1ABC.cif' || url === RCSB + '1ABC-assembly1.cif') return send(200, MMCIF);
    if (url === RCSB + '2ABC.cif') { await new Promise(r => { hold = r; setTimeout(r, 4000); }); try { await send(200, MMCIF); } catch (e) { /* the page gave up */ } return; }
    if (url.startsWith(PUG) || url.startsWith(RCSB)) return send(404, 'Not Found');
    return route.abort();
  });

  await page.goto(HTML);
  await page.waitForFunction(() => window.BM && window.BM.app && window.BM.fetchdb);
  aspirin = await page.evaluate(() => {
    const m = BM.chem.smilesTo3D('CC(=O)Oc1ccccc1C(=O)O', { seed: 3 }).mol;
    return BM.writers.mol({ Z: m.Z, xyz: m.xyz, bonds: m.bonds.map(b => [b[0], b[1], Math.round(b[2]) || 1]), title: '2244' }).replace(/\n?$/, '\n') + '> <PUBCHEM_COMPOUND_CID>\n2244\n\n$$$$\n';
  });
  check(requests.length === 0, 'no request at start-up: ' + requests.map(r => r.url).join(' '));
  check(await page.evaluate(() => { try { return localStorage.getItem('bluemolecule.downloads.allow'); } catch (e) { return 'x'; } }) === null, 'downloads not allowed yet');

  const dlgOpen = () => page.evaluate(() => document.querySelector('#dlg').open);
  const msg = () => page.$eval('[data-fdb=msg]', e => ({ text: e.textContent, err: e.classList.contains('err') }));
  const goBtn = '#dlg [data-dbtn="1"]', cancelBtn = '#dlg [data-dbtn="0"]';
  const docOpened = async name => { await page.waitForFunction(n => BM.app.state.active && BM.app.state.active.name === n, name, { timeout: 15000 }); };

  // ------------------------------------------------------------ the empty screen's button, then the toolbar's
  check(await page.isVisible('#btn-fetch-2'), 'a PubChem or PDB button on the empty screen');
  await page.click('#btn-fetch-2');
  check(await dlgOpen(), 'the dialog opens from the empty screen');
  await page.click(cancelBtn);
  check(!(await dlgOpen()), 'Cancel closes it');
  await page.click('#btn-fetch');
  check(await dlgOpen() && (await page.$eval('#dlg .dlg-head', e => e.textContent)) === 'Download a structure', 'the dialog opens from the toolbar');
  const note = await page.$eval('#dlg .fdb', e => e.innerText);
  check(/pubchem\.ncbi\.nlm\.nih\.gov/.test(note) && /files\.rcsb\.org/.test(note) && /never leave your computer/.test(note), 'the dialog names both hosts and what is sent');
  check(await page.$eval(goBtn, b => b.disabled), 'Download is disabled until downloads are allowed');
  check(await page.evaluate(() => document.activeElement && document.activeElement.dataset.fdb === 'q'), 'the query field has the focus');
  await page.fill('[data-fdb=q]', 'aspirin');
  await page.press('[data-fdb=q]', 'Enter');
  check(await dlgOpen(), 'Enter does not close the dialog');
  let m = await msg();
  check(m.err && /Tick/.test(m.text), 'Enter while not allowed asks for the tick: ' + m.text);
  check(requests.length === 0, 'no request while downloads are not allowed');

  // ------------------------------------------------------------ PubChem by name
  await page.check('[data-fdb=allow]');
  check(await page.evaluate(() => localStorage.getItem('bluemolecule.downloads.allow')) === '1', 'the permission is remembered');
  check(!(await page.$eval(goBtn, b => b.disabled)), 'Download enabled once allowed');
  check(requests.length === 0, 'ticking the box sends nothing');
  await page.click(goBtn);
  await docOpened('aspirin_CID2244.sdf');
  check(!(await dlgOpen()), 'the dialog closes after the download');
  check(requests.length === 1 && requests[0].url === PUG + 'name/aspirin/SDF?record_type=3d', 'one request, to PubChem: ' + requests.map(r => r.url).join(' '));
  check(requests[0].method === 'GET' && !requests[0].body && !requests[0].headers.cookie && !requests[0].headers.referer, 'a GET without a body, cookies or referrer: ' + JSON.stringify(requests[0].headers));
  let d = await page.evaluate(() => { const d = BM.app.state.active; return { n: BM.app.curFrame(d).Z.length, title: d.data.title, fmt: d.data.format, dl: d.data.downloaded }; });
  check(d.n === 21 && d.fmt === 'mol' && d.title === 'aspirin · PubChem CID 2244', 'aspirin opened as an SD file: ' + JSON.stringify(d));
  check(d.dl && d.dl.from === 'PubChem CID 2244' && d.dl.url === requests[0].url && /^\d{4}-\d\d-\d\d$/.test(d.dl.date), 'where it came from is kept: ' + JSON.stringify(d.dl));
  await page.click('#tabs button[data-tab="info"]');
  const infoTxt = await page.$eval('#panel', e => e.innerText);
  check(/Downloaded\s*PubChem CID 2244 · \d{4}-\d\d-\d\d/.test(infoTxt), 'the Info tab says where it came from');

  // ------------------------------------------------------------ PubChem without a 3D conformer: the SMILES
  await page.click('#btn-fetch');
  check(await page.isChecked('[data-fdb=allow]'), 'the box stays ticked');
  await page.fill('[data-fdb=q]', 'sodium chloride');
  await page.press('[data-fdb=q]', 'Enter');
  await docOpened('sodium_chloride_CID5234');
  d = await page.evaluate(() => { const d = BM.app.state.active; return { Z: Array.from(BM.app.curFrame(d).Z), warn: d.data.warnings.join(' '), tab: d.tab, title: d.data.title }; });
  check(d.Z.length === 2 && d.Z.includes(11) && d.Z.includes(17), 'built from the SMILES: ' + JSON.stringify(d.Z));
  check(/built one from its SMILES, \[Na\+\]\.\[Cl-\]/.test(d.warn) && d.tab === 'info', 'and says so, on the Info tab: ' + d.warn.slice(0, 120));
  check(requests.slice(1).map(r => r.url).join(' ') === PUG + 'name/sodium%20chloride/SDF?record_type=3d ' + PUG + 'name/sodium%20chloride/property/SMILES/JSON', 'the conformer, then the SMILES');

  // ------------------------------------------------------------ the PDB, an assembly, errors
  await page.click('#btn-fetch');
  check(!(await page.isVisible('[data-fdb=asmrow]')), 'no assembly box for PubChem');
  await page.check('input[name=fdb-src][value=pdb]');
  check(await page.isVisible('[data-fdb=asmrow]') && (await page.$eval('[data-fdb=ql]', e => e.textContent)) === 'PDB ID', 'the PDB fields');
  let n0 = requests.length;
  await page.fill('[data-fdb=q]', 'abc');
  await page.click(goBtn);
  m = await msg();
  check(await dlgOpen() && m.err && /four characters/.test(m.text) && requests.length === n0, 'a bad ID is explained and sends nothing: ' + m.text);
  await page.fill('[data-fdb=q]', '9zzz');
  await page.click(goBtn);
  await page.waitForFunction(() => /no entry/.test(document.querySelector('[data-fdb=msg]').textContent), null, { timeout: 10000 });
  m = await msg();
  check(await dlgOpen() && m.err && /The PDB has no entry 9ZZZ/.test(m.text), 'a missing entry is said in the dialog: ' + m.text);
  check(!(await page.$eval(goBtn, b => b.disabled)) && !(await page.$eval('[data-fdb=q]', e => e.disabled)), 'and the dialog can be used again');
  await page.fill('[data-fdb=q]', '1abc');
  await page.check('[data-fdb=asm]');
  await page.click(goBtn);
  await docOpened('1ABC-assembly1.cif');
  d = await page.evaluate(() => { const d = BM.app.state.active; return { fmt: d.data.format, title: d.data.title, dl: d.data.downloaded, n: BM.app.curFrame(d).Z.length }; });
  check(d.fmt === 'mmcif' && d.title === 'Tiny test peptide' && d.n > 0, 'the assembly opened as mmCIF: ' + JSON.stringify(d));
  check(d.dl && d.dl.from === 'RCSB Protein Data Bank, biological assembly 1' && requests[requests.length - 1].url === RCSB + '1ABC-assembly1.cif', 'assembly 1 asked for: ' + requests[requests.length - 1].url);
  const prefs = await page.evaluate(() => JSON.parse(localStorage.getItem('bluemolecule.downloads')));
  check(prefs.source === 'pdb' && prefs.assembly === true, 'the choice of database is remembered');

  // ------------------------------------------------------------ Cancel stops a download
  await page.click('#btn-fetch');
  check(await page.isChecked('input[name=fdb-src][value=pdb]'), 'the dialog opens on the PDB again');
  await page.uncheck('[data-fdb=asm]');
  await page.fill('[data-fdb=q]', '2abc');
  const docs0 = await page.evaluate(() => BM.app.state.docs.length);
  await page.click(goBtn);
  await page.waitForFunction(() => /Downloading/.test(document.querySelector('[data-fdb=msg]').textContent));
  check(await page.$eval(goBtn, b => b.disabled) && await page.$eval('[data-fdb=q]', e => e.disabled), 'inputs wait during the download');
  await page.click(cancelBtn);
  check(!(await dlgOpen()), 'Cancel closes the dialog');
  if (hold) hold();
  await page.waitForTimeout(600);
  check((await page.evaluate(() => BM.app.state.docs.length)) === docs0, 'a cancelled download opens nothing');

  // ------------------------------------------------------------ unticking turns downloads off
  await page.click('#btn-fetch');
  await page.uncheck('[data-fdb=allow]');
  check(await page.evaluate(() => localStorage.getItem('bluemolecule.downloads.allow')) === '0', 'unticking is remembered');
  n0 = requests.length;
  await page.fill('[data-fdb=q]', '1abc');
  await page.press('[data-fdb=q]', 'Enter');
  check(await page.$eval(goBtn, b => b.disabled) && requests.length === n0, 'nothing is sent once unticked');
  await page.keyboard.press('Escape');
  check(!(await dlgOpen()), 'Esc closes the dialog');

  const hosts = new Set(requests.map(r => new URL(r.url).host));
  check([...hosts].every(h => h === 'pubchem.ncbi.nlm.nih.gov' || h === 'files.rcsb.org'), 'only the two databases were asked: ' + [...hosts].join(', '));
  check(!errors.length, 'no page errors: ' + errors.join(' | '));
  await browser.close();
  console.log(`Download browser tests: ${passes} passed, ${failures} failed`);
  process.exit(failures ? 1 : 0);
})().catch(e => { console.error(e); process.exit(1); });
