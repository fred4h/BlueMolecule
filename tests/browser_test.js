// End-to-end checks of dist/BlueMolecule.html in headless Chromium.
// Usage: node tests/browser_test.js [screenshotDir]
const { chromium } = require('playwright');
const path = require('path');
const fs = require('fs');
const { spawn } = require('child_process');

const ROOT = path.resolve(__dirname, '..');
const HTML = 'file://' + path.join(ROOT, 'dist', 'BlueMolecule.html');
const FIX = path.join(ROOT, 'tests', 'fixtures');
const SHOTS = process.argv[2] || path.join(ROOT, 'tests', 'shots');
fs.mkdirSync(SHOTS, { recursive: true });

let failures = 0, passes = 0;
const check = (cond, msg) => { if (cond) passes++; else { failures++; console.log('  FAIL:', msg); } };

(async () => {
  const browser = await chromium.launch({ args: ['--enable-unsafe-swiftshader'] });
  const ctx = await browser.newContext({ viewport: { width: 1440, height: 900 }, deviceScaleFactor: 1 });
  const page = await ctx.newPage();
  const errors = [];
  page.on('pageerror', e => errors.push('pageerror: ' + e.message));
  page.on('console', m => { if (m.type() === 'error') errors.push('console: ' + m.text()); });
  await page.goto(HTML);
  await page.waitForFunction(() => window.BM && window.BM.app);
  await page.evaluate(() => {
    window.__saved = [];
    window.__bmSaveHook = async (name, blob) => { window.__saved.push({ name, text: blob.type === 'image/png' ? null : await blob.text(), size: blob.size, type: blob.type }); };
  });
  await page.screenshot({ path: path.join(SHOTS, '00_empty.png') });

  const open = async (file) => {
    await page.setInputFiles('#file-input', path.join(FIX, file));
    await page.waitForFunction(n => { const s = window.BM.app.state; return s.active && s.active.name === n; }, file, { timeout: 15000 });
    await page.waitForTimeout(300);
  };
  const st = () => page.evaluate(() => {
    const d = window.BM.app.state.active;
    return { name: d.name, frames: d.data.frames.length, fi: d.frameIndex, atoms: d.data.frames[d.frameIndex].Z.length, vib: d.data.vibrations ? d.data.vibrations.freqs.length : 0, plotType: d.plotType, tab: d.tab, docs: window.BM.app.state.docs.length, nSphere: window.BM.app.R().nSphere, nCyl: window.BM.app.R().nCyl, surfaces: window.BM.app.R().surfaces.length };
  });
  const shot = async name => page.screenshot({ path: path.join(SHOTS, name + '.png') });

  // 1. Gaussian opt+freq
  console.log('Gaussian opt+freq');
  await open('g16_water_optfreq.log');
  let s = await st();
  check(s.frames === 4 && s.atoms === 3 && s.vib === 3, 'water opt+freq state ' + JSON.stringify(s));
  check(s.nSphere === 3 && s.nCyl === 4, 'spheres/bonds drawn ' + s.nSphere + '/' + s.nCyl);
  check(s.plotType === 'ir' && s.tab === 'vib', 'default IR plot + vib tab');
  await shot('01_gaussian_vib_tab');
  await page.click('#vib-table tr.row[data-i="0"]');
  await page.waitForTimeout(500);
  const animating = await page.evaluate(() => window.BM.app.state.anim.kind);
  check(animating === 'vib', 'vibration animates after click: ' + animating);
  await page.evaluate(() => { window.BM.app.state.vib.arrows = true; });
  await page.click('#playbar button[data-pb="arrows"]');
  await page.waitForTimeout(300);
  await shot('02_gaussian_vib_anim');
  await page.selectOption('#plot-pick', 'raman');
  await page.waitForTimeout(200);
  await shot('03_raman');
  await page.selectOption('#plot-pick', 'opt:energy');
  await page.waitForTimeout(200);
  await page.keyboard.press('Escape');
  await page.click('#tabs button[data-tab="steps"]');
  await page.click('#steps-table tr.row[data-i="0"]');
  s = await st();
  check(s.fi === 0, 'step click selects geometry 1');
  await shot('04_steps_energy');
  await page.click('#tabs button[data-tab="info"]');
  await page.waitForTimeout(100);
  await shot('05_info');

  // measurement + edit
  console.log('Measure and edit');
  await page.evaluate(() => { const a = window.BM.app; a.setFrame(a.state.active, a.state.active.data.frames.length - 1); a.onAtomClick(1); a.onAtomClick(0); a.onAtomClick(2); });
  const chip = await page.textContent('#measure');
  check(/angle/.test(chip), 'angle measurement chip: ' + chip);
  await page.click('#tabs button[data-tab="build"]');
  await page.fill('#edit-val', '109.5');
  await page.click('[data-act="edit-set"]');
  await page.waitForTimeout(150);
  const ang = await page.evaluate(() => { const a = window.BM.app, d = a.state.active, f = d.data.frames[d.frameIndex]; return [f.edited, BM.geom.angle(f.xyz, 1, 0, 2)]; });
  check(ang[0] === true && Math.abs(ang[1] - 109.5) < 1e-6, 'angle edited to 109.5: ' + ang);
  await shot('06_edit');
  await page.click('[data-act="undo"]');
  const angU = await page.evaluate(() => { const d = window.BM.app.state.active, f = d.data.frames[d.frameIndex]; return BM.geom.angle(f.xyz, 1, 0, 2); });
  check(Math.abs(angU - 109.5) > 1, 'undo restores angle: ' + angU);

  // exports
  console.log('Exports');
  await page.evaluate(() => window.BM.app.exportSimple('xyz'));
  await page.evaluate(() => window.BM.app.exportSimple('mol'));
  await page.evaluate(() => window.BM.app.exportSimple('pdb'));
  await page.evaluate(() => window.BM.app.exportSimple('xyz-all'));
  await page.evaluate(() => window.BM.app.exportGaussian());
  await page.waitForSelector('#dlg[open]');
  await shot('07_gjf_dialog');
  await page.click('#dlg .btn.primary');
  await page.evaluate(() => window.BM.app.exportOrca());
  await page.click('#dlg .btn.primary');
  await page.evaluate(() => { const d = window.BM.app.state.active; d.plotType = 'ir'; window.BM.app.plotCSV(); });
  await page.click('[data-act="plot-svg"]');
  await page.waitForTimeout(300);
  const saved = await page.evaluate(() => window.__saved.map(x => ({ name: x.name, text: x.text, size: x.size })));
  const names = saved.map(x => x.name);
  console.log('  saved:', names.join(', '));
  check(names.length >= 7, 'seven exports saved');
  const reparsed = await page.evaluate(files => files.filter(f => /\.(xyz|mol|pdb|gjf|inp)$/.test(f.name) && !/\.inp$/.test(f.name)).map(f => { try { const d = BM.parseFile(f.name, f.text); return [f.name, d.frames.length, d.frames[0].Z.length]; } catch (e) { return [f.name, 'ERR ' + e.message]; } }), saved);
  console.log('  reparsed:', JSON.stringify(reparsed));
  check(reparsed.every(r => typeof r[1] === 'number' && r[2] === 3), 'exports parse back');
  const inp = saved.find(x => x.name.endsWith('.inp'));
  check(inp && /\* xyz 0 1/.test(inp.text) && /^! /m.test(inp.text), 'ORCA input looks right');
  const svg = saved.find(x => x.name.endsWith('.svg'));
  check(svg && svg.text.startsWith('<?xml') && svg.text.includes('<polyline'), 'SVG spectrum export');

  // 2. ORCA TD-DFT (UV-Vis)
  console.log('ORCA TD-DFT');
  await open('orca5_water_tddft.out');
  s = await st();
  check(s.plotType === 'uv' && s.docs === 2, 'ORCA TD-DFT → UV plot, two docs ' + JSON.stringify(s));
  await page.click('#states-table tr.row[data-i="1"]');
  await page.waitForTimeout(150);
  await shot('08_orca_uv');

  // 3. cube
  console.log('Cube');
  await open('water_homo_bohr.cube');
  s = await st();
  check(s.surfaces === 2, 'cube: two lobes drawn (' + s.surfaces + ')');
  await shot('09_cube');
  await page.fill('.surf input[data-s="iso"]', '0.05');
  await page.dispatchEvent('.surf input[data-s="iso"]', 'input');
  await page.waitForTimeout(200);
  // built-in example cube (finer grid)
  await page.evaluate(() => { const ex = BM.examples.find(e => e.id === 'benzene-cube'); window.BM.app.openText(ex.name, ex.text); });
  await page.waitForTimeout(400);
  await page.evaluate(() => { const d = window.BM.app.state.active; d.surfaces[0].opacity = 0.8; });
  await page.dispatchEvent('.surf input[data-s="opacity"]', 'input');
  await page.waitForTimeout(300);
  await shot('10_benzene_orbital');

  // 4. cisplatin, labels, styles
  console.log('Cisplatin and styles');
  await open('cisplatin.xyz');
  s = await st();
  check(s.atoms === 11 && s.nCyl === 20, 'cisplatin 11 atoms 10 bonds ' + JSON.stringify(s));
  await page.selectOption('#sel-labels', 'symnum');
  await page.click('#tg-bl');
  await page.selectOption('#sel-bg', 'white');
  await page.waitForTimeout(200);
  await shot('11_cisplatin_labels_white');
  for (const [k, n] of [['2', 'sticks'], ['3', 'spacefill'], ['4', 'wireframe'], ['1', 'ballstick']]) {
    await page.mouse.move(700, 400);
    await page.keyboard.press(k);
    await page.waitForTimeout(120);
    if (n !== 'ballstick') await shot('12_style_' + n);
  }
  await page.selectOption('#sel-labels', 'none');
  await page.click('#tg-bl');
  await page.selectOption('#sel-bg', 'navy');

  // 5. other formats
  console.log('Other formats');
  const others = ['g16_water_hpmodes.log', 'g16_water_tddft.log', 'g09_cisplatin_input_orientation.log', 'g16_water_scan.log', 'methanol_cartesian.gjf', 'nh3_zmatrix.com', 'cisplatin.gjf',
    'orca5_water_optfreq.out', 'orca4_water_freq.out', 'orca5_water_scan.out', 'orca6_water_tddft.out', 'water_traj.xyz', 'water_mo5_mo6.cube', 'water_density_angstrom.cube',
    'water_freq.molden', 'water_freq.fchk', 'acetic_acid.sdf', 'peptide_ions_elements.pdb', 'peptide_ions_noelements.pdb', 'water_ions_cisplatin.gro'];
  for (const f of others) {
    const before = errors.length;
    await open(f);
    const r = await st();
    check(r.atoms > 0 && r.nSphere > 0, f + ' renders ' + JSON.stringify(r));
    check(errors.length === before, f + ' no errors: ' + errors.slice(before).join(' | '));
  }
  await shot('13_scan_orca');
  await page.evaluate(() => { const a = window.BM.app; a.activate(a.state.docs.find(d => d.name === 'g16_water_scan.log')); });
  await page.waitForTimeout(200);
  await shot('14_gaussian_scan');
  await page.evaluate(() => { const a = window.BM.app; a.activate(a.state.docs.find(d => d.name === 'water_ions_cisplatin.gro')); });
  await page.waitForTimeout(200);
  await shot('15_gro');

  // image export
  await page.evaluate(() => window.BM.app.R().snapshot(2, true).then(b => window.__saved.push({ name: 'snap.png', size: b.size, type: b.type })));
  await page.waitForTimeout(500);
  const snap = await page.evaluate(() => window.__saved.find(x => x.name === 'snap.png'));
  check(snap && snap.size > 5000, 'PNG snapshot ' + (snap && snap.size));

  // gz file
  const zlib = require('zlib');
  const gzPath = path.join(SHOTS, 'water_freq.log.gz');
  fs.writeFileSync(gzPath, zlib.gzipSync(fs.readFileSync(path.join(FIX, 'g16_water_optfreq.log'))));
  await page.setInputFiles('#file-input', gzPath);
  await page.waitForFunction(() => window.BM.app.state.active.name === 'water_freq.log', null, { timeout: 15000 });
  check(true, 'gz opens');

  // light theme
  await page.click('#btn-theme');
  await page.waitForTimeout(200);
  await shot('16_light_theme');
  await page.click('#btn-theme');

  // help dialog
  await page.keyboard.press('?');
  await page.waitForSelector('#dlg[open]');
  await shot('17_help');
  await page.keyboard.press('Escape');

  console.log('Console/page errors:', errors.length ? errors : 'none');
  check(errors.length === 0, 'no console errors');
  await browser.close();

  // 6. launcher (serves files to the page)
  console.log('Launcher');
  const proc = spawn('python3', [path.join(ROOT, 'dist', 'build', require('fs').readdirSync(path.join(ROOT, 'dist', 'build')).find(d => /^BlueMolecule-[\d.]+-linux$/.test(d)), 'bluemolecule'), '--no-open', path.join(FIX, 'orca5_water_optfreq.out'), path.join(FIX, 'water_homo_bohr.cube')]);
  const url = await new Promise((res, rej) => { proc.stdout.on('data', d => res(String(d).trim())); setTimeout(() => rej(new Error('launcher timeout')), 8000); });
  const b2 = await chromium.launch();
  const p2 = await b2.newPage();
  await p2.goto(url);
  await p2.waitForFunction(() => window.BM && window.BM.app && window.BM.app.state.docs.length === 2, null, { timeout: 15000 });
  const names2 = await p2.evaluate(() => window.BM.app.state.docs.map(d => d.name));
  check(names2.join(',') === 'orca5_water_optfreq.out,water_homo_bohr.cube', 'launcher served both files: ' + names2);
  const bad = await p2.evaluate(async u => (await fetch(u.replace(/\/[^/]+\/$/, '/wrongtoken/'))).status, url);
  check(bad === 404, 'launcher rejects wrong token: ' + bad);
  await p2.screenshot({ path: path.join(SHOTS, '18_launcher.png') });
  await b2.close();
  const exited = await new Promise(res => { const t = setTimeout(() => res(false), 40000); proc.on('exit', () => { clearTimeout(t); res(true); }); });
  check(exited, 'launcher exits after the window closes');

  console.log(`\nBrowser tests: ${passes} passed, ${failures} failed`);
  process.exit(failures ? 1 : 0);
})().catch(e => { console.error(e); process.exit(2); });
