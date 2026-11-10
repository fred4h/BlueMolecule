#!/usr/bin/env node
/* Regenerates the screenshots in docs/images from the built viewer, the test fixtures and the synthetic
 * macromolecules of tests/synth.js (needs Playwright with Chromium):
 *     python3 tools/build.py && node tools/docs_screenshots.js
 */
'use strict';
const { chromium } = require('playwright');
const path = require('path');
const fs = require('fs');

const ROOT = path.resolve(__dirname, '..');
const HTML = 'file://' + path.join(ROOT, 'dist', 'BlueMolecule.html');
const FIX = path.join(ROOT, 'tests', 'fixtures');
const OUT = path.join(ROOT, 'docs', 'images');

(async () => {
  if (!fs.existsSync(HTML.slice(7))) { console.error('dist/BlueMolecule.html not found: run python3 tools/build.py first'); process.exit(2); }
  fs.mkdirSync(OUT, { recursive: true });
  const browser = await chromium.launch({ args: ['--enable-unsafe-swiftshader'] });
  const ctx = await browser.newContext({ viewport: { width: 1440, height: 900 }, deviceScaleFactor: 1 });
  const page = await ctx.newPage();
  await page.goto(HTML);
  await page.waitForFunction(() => window.BM && window.BM.app);
  await page.addStyleTag({ content: '#toast { display: none !important; }' });
  await page.addScriptTag({ path: path.join(ROOT, 'tests', 'synth.js') });
  const ev = (fn, arg) => page.evaluate(fn, arg);
  const settle = (ms = 350) => page.waitForTimeout(ms);
  const shot = async name => { await settle(); await page.mouse.move(1, 899); await settle(120); await page.screenshot({ path: path.join(OUT, name + '.png') }); console.log('docs/images/' + name + '.png'); };
  const closeAll = () => ev(() => { for (const b of [...document.querySelectorAll('#doc-list [data-close]')]) b.click(); });
  const open = async file => {
    await page.setInputFiles('#file-input', path.join(FIX, file));
    await page.waitForFunction(n => { const s = window.BM.app.state; return s.active && s.active.name === n; }, path.basename(file), { timeout: 15000 });
    await settle();
  };
  const openSynth = async (name, code) => {
    await ev(([nm, c]) => { const S = window.BMSynth, sys = new Function('S', 'BM', 'return ' + c)(S, BM); BM.app.openText(nm, /\.gro$/.test(nm) ? S.toGRO(sys) : S.toPDB(sys, { ss: true })); }, [name, code]);
    await settle(600);
  };
  const tab = async id => { await page.click(`#tabs [data-tab="${id}"]`); await settle(150); };
  const plot = async id => { await page.selectOption('#plot-pick', id); await settle(200); };
  const pause = () => ev(() => { const b = document.querySelector('#playbar [data-act="play"].on, #playbar .playing'); if (b) b.click(); });

  // 1. building from SMILES, Build tab
  await page.fill('#empty-smiles', 'CC(=O)Oc1ccccc1C(=O)O');
  await page.click('#btn-smiles-2');
  await settle(800);
  await tab('build');
  await ev(() => BM.app.setTool('draw'));
  await shot('build');
  await ev(() => BM.app.setTool('navigate'));
  await closeAll();

  // 2. optimisation + frequencies: Info contents and the convergence plot
  await open('g16_water_optfreq.log');
  await tab('info');
  await plot('opt:all');
  await shot('optimisation');

  // 3. the same file: vibrations and the IR spectrum
  await tab('vib');
  await page.click('#vib-table tr.row[data-i="2"]');
  await ev(() => { const c = document.querySelector('[data-vib="arrows"]'); if (c && !c.checked) c.click(); });
  await plot('ir');
  await settle(700);
  await shot('vibrations');
  await closeAll();

  // 4. excited states and the UV-Vis spectrum
  await open('orca5_water_tddft.out');
  await tab('states');
  await page.click('#states-table tr.row[data-i="0"]');
  await plot('uv');
  await shot('uv-vis');
  await closeAll();

  // 5. three-layer ONIOM: QM/MM layer preset with legend, Display tab
  await open('qmmm/g16_oniom3_opt.log');
  await ev(() => { const d = BM.app.state.active, R = BM.app.R(); BM.app.setFrame(d, d.data.frames.length - 1); R.dist *= 0.7; R.requestRender(); });
  await tab('display');
  await shot('qmmm-layers');
  await closeAll();

  // 6. macromolecules (synthetic structures from tests/synth.js)
  await openSynth('protein.pdb', 'S.protein(BM)');
  await tab('display');
  await shot('protein-ligand');
  await closeAll();
  await openSynth('hspg.pdb', 'S.glycan(BM)');
  await tab('info');
  await shot('glycan-snfg');
  await closeAll();
  await openSynth('membrane.gro', 'S.membrane({ n: 10, protein: true })');
  await tab('info');
  await shot('membrane');

  await browser.close();
})().catch(e => { console.error(e); process.exit(1); });
