#!/usr/bin/env node
/* Injection test: text read from files (titles, names, labels, comments) and file names must never become markup.
 * Opens crafted files in the built viewer, visits every panel tab, hovers the 3D view, opens Save as and Help,
 * and fails if any injected element or event handler made it into the page.
 * Run after the build:   node tests/browser_xss_test.js */
'use strict';
const { chromium } = require('playwright');
const path = require('path');
const fs = require('fs');

const ROOT = path.resolve(__dirname, '..');
const HTML = 'file://' + path.join(ROOT, 'dist', 'BlueMolecule.html');
const FIX = path.join(ROOT, 'tests', 'fixtures');
const HIT = 'window.__xss=(window.__xss||0)+1';
const SP = `<b data-pwn=1>B</b><img src="x" onerror="${HIT}">`;      // for free-text fields
const NS = `<b/data-pwn=1>B</b><img/src="x"/onerror="${HIT}">`;      // for whitespace-separated fields
const xmlEsc = s => s.replace(/&/g, '&amp;').replace(/</g, '&lt;').replace(/>/g, '&gt;').replace(/"/g, '&quot;');
const rd = f => fs.readFileSync(path.join(FIX, f), 'utf8');
const lines = f => rd(f).split('\n');

let passes = 0, failures = 0;
const check = (cond, msg) => { if (cond) passes++; else { failures++; console.log('  FAIL:', msg); } };

function cases() {
  const out = [];
  const add = (name, text) => out.push({ name, text });
  const water = 'O 0 0 0.117\nH 0 0.757 -0.469\nH 0 -0.757 -0.469\n';
  add('w' + SP + '.xyz', '3\n' + SP + '\n' + water);
  add('traj.xyz', [1, 2, 3].map(k => '3\n' + SP + ' step ' + k + '\n' + water).join(''));
  add('badel.xyz', '1\nc\n' + NS + ' 0 0 0\n');
  let L = lines('g16_water_optfreq.log');
  L[38] = ' ' + SP;
  L = L.map(s => (/^ +A1 +A1 +B2$/.test(s) ? '   ' + [NS, NS, NS].join('   ') : s));
  add('g16_title.log', L.join('\n'));
  add('g16_states.log', rd('g16_water_tddft.log').replace(/(Singlet|Triplet)-B1/g, '$1-' + NS));
  add('orca_kw.out', rd('orca5_water_optfreq.out').replace('! B3LYP def2-SVP Opt Freq', '! B3LYP def2-SVP Opt Freq ' + SP));
  L = lines('water_homo_bohr.cube'); L[0] = ' ' + SP; L[1] = ' ' + SP;
  add('title.cube', L.join('\n'));
  L = lines('water_freq.molden'); L[2] = SP;
  add('title.molden', L.join('\n').replace(/Sym=\s*\S+/g, 'Sym= ' + NS));
  L = lines('water_freq.fchk'); L[0] = SP;
  add('title.fchk', L.join('\n'));
  L = lines('struct/ligand.mol2'); L[2] = SP;
  add('names.mol2', L.join('\n').replace(/BEN1/g, NS).replace(' C1 ', ' ' + NS + ' '));
  add('title.cml', rd('struct/attr.cml').replace('title="methylammonium"', 'title="' + xmlEsc(SP) + '"'));
  add('name.cjson', rd('struct/water.cjson').replace('"name": "water"', '"name": ' + JSON.stringify(SP)));
  L = lines('acetic_acid.sdf'); L[0] = SP; L[2] = SP;
  const tag = L.findIndex(s => s.startsWith('> <'));
  if (tag >= 0) L[tag + 1] = SP;
  add('title.sdf', L.join('\n'));
  L = lines('water_ions_cisplatin.gro'); L[0] = SP;
  add('title.gro', L.join('\n'));
  add('title.pdb', 'TITLE     ' + SP + '\nCOMPND    ' + SP + '\n' + rd('peptide_ions_elements.pdb'));
  add('oniom_types.log', rd('qmmm/g16_oniom_amber_opt.log').replace('-CT-', '-' + NS + '-'));
  return out;
}

(async () => {
  if (!fs.existsSync(HTML.slice(7))) { console.log('dist/BlueMolecule.html not found: run python3 tools/build.py first'); process.exit(2); }
  const browser = await chromium.launch({ args: ['--enable-unsafe-swiftshader'] });
  const page = await browser.newPage({ viewport: { width: 1400, height: 900 } });
  const errors = [];
  page.on('pageerror', e => errors.push(String(e)));
  page.on('dialog', d => { errors.push('dialog: ' + d.message()); d.dismiss(); });
  await page.goto(HTML);
  await page.waitForFunction(() => window.BM && window.BM.app && window.BM.app.state);
  const injected = () => page.evaluate(() => ({ xss: window.__xss || 0, nodes: document.querySelectorAll('[data-pwn], img[src="x"]').length }));
  const shown = new Set();
  // self-check: the detector must see markup that really is injected
  await page.evaluate(sp => { const d = document.createElement('div'); d.id = 'selfcheck'; d.innerHTML = sp; document.body.appendChild(d); }, SP);
  await page.waitForTimeout(150);
  const s0 = await injected();
  check(s0.xss >= 1 && s0.nodes >= 1, 'self-check: the detector sees injected markup');
  await page.evaluate(() => { document.getElementById('selfcheck').remove(); window.__xss = 0; });

  for (const c of cases()) {
    const before = await page.evaluate(() => { const s = window.BM.app.state; return { n: s.docs.length, t: document.querySelector('#toast').textContent }; });
    await page.setInputFiles('#file-input', { name: c.name, mimeType: 'application/octet-stream', buffer: Buffer.from(c.text) });
    await page.waitForFunction(b => { const s = window.BM.app.state; return s.docs.length !== b.n || document.querySelector('#toast').textContent !== b.t; }, before, { timeout: 15000 }).catch(() => {});
    await page.waitForTimeout(150);
    const tabs = await page.evaluate(() => [...document.querySelectorAll('#tabs [data-tab]')].map(b => b.dataset.tab));
    for (const t of tabs) {
      await page.click(`#tabs [data-tab="${t}"]`).catch(() => {});
      await page.waitForTimeout(60);
      if (await page.evaluate(() => document.body.innerText.includes('data-pwn'))) shown.add(c.name + ':' + t);
    }
    await page.mouse.move(700, 450); await page.waitForTimeout(80);
    await page.mouse.move(705, 455); await page.waitForTimeout(80);
    const r = await injected();
    check(r.xss === 0 && r.nodes === 0, `${c.name}: markup from the file reached the page (handlers ${r.xss}, nodes ${r.nodes})`);
  }
  // dialogs that show names taken from the active document
  for (const fn of ['saveAs', 'exportCalc']) {
    await page.evaluate(f => { try { const a = window.BM.app; a.activate(a.state.docs[0]); a[f](); } catch (e) { /* not available for this document */ } }, fn);
    await page.waitForTimeout(120);
    const r = await injected();
    check(r.xss === 0 && r.nodes === 0, `${fn} dialog: injected markup (handlers ${r.xss}, nodes ${r.nodes})`);
    await page.keyboard.press('Escape');
  }
  await page.keyboard.press('?'); await page.waitForTimeout(80);
  const docs = await page.evaluate(() => window.BM.app.state.docs.length);
  check(docs >= 12, 'most crafted files opened (' + docs + ' documents)');
  check(shown.size >= 5, 'the payload text was displayed (as text) in several places: ' + shown.size);
  check(errors.filter(e => !/^dialog/.test(e)).length === 0, 'no page errors: ' + errors.join(' | '));
  check(!errors.some(e => /^dialog/.test(e)), 'no alert/confirm dialogs were opened');
  await browser.close();
  console.log(`browser_xss_test: ${passes} passed, ${failures} failed (payload shown as text in ${shown.size} places)`);
  process.exit(failures ? 1 : 0);
})().catch(e => { console.error(e); process.exit(1); });
