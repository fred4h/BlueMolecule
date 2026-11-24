// Browser checks of 2.4 in dist/BlueMolecule.html (headless Chromium): orbitals and densities computed from fchk and
// Molden files in workers (#17), molecular surfaces and the electrostatic potential on surfaces (#18).
//   - the Surfaces tab: orbital menu with HOMO and LUMO, Compute, a progress bar, one new surface shown at a time
//   - the result matches BM.orbitals.evaluate on the same grid (the workers compute what the page would)
//   - density coloured by the ESP of the charges; a potential cube coloring a density cube of another document
//   - unrestricted files: α/β menus and the spin density; the MO level plot's "Show orbital"
//   - molecular surfaces (vdW, accessible, excluded); Save as Molden writes the orbitals; text from files is escaped
// Usage: node tests/browser_orbitals_test.js
const { chromium } = require('playwright');
const path = require('path');
const fs = require('fs');

const ROOT = path.resolve(__dirname, '..');
const HTML = 'file://' + path.join(ROOT, 'dist', 'BlueMolecule.html');
const ORB = path.join(ROOT, 'tests', 'fixtures', 'orbitals');

let failures = 0, passes = 0;
const check = (cond, msg) => { if (cond) passes++; else { failures++; console.log('  FAIL:', msg); } };

(async () => {
  const browser = await chromium.launch({ args: ['--enable-unsafe-swiftshader'] });
  const ctx = await browser.newContext({ viewport: { width: 1440, height: 1000 }, deviceScaleFactor: 1 });
  const page = await ctx.newPage();
  const errors = [];
  page.on('pageerror', e => errors.push('pageerror: ' + e.message));
  page.on('console', m => { if (m.type() === 'error') errors.push('console: ' + m.text()); });
  page.on('dialog', d => { errors.push('dialog: ' + d.message()); d.dismiss(); });
  await page.goto(HTML);
  await page.waitForFunction(() => window.BM && window.BM.app);
  const open = async (name, text) => {
    await page.evaluate(([n, t]) => BM.app.openText(n, t), [name, text]);
    await page.waitForFunction(n => BM.app.state.active && BM.app.state.active.name === n, name);
  };
  const tab = async id => { await page.click(`#tabs button[data-tab="${id}"]`); await page.waitForTimeout(100); };
  const idle = sel => page.waitForFunction(s => !document.querySelector(s + ' .sf-job'), sel, { timeout: 60000 });
  const cards = () => page.$$eval('#panel .surf', cs => cs.map(c => ({ name: c.querySelector('.t').textContent, kind: c.querySelector('.badge').textContent, on: c.querySelector('[data-s=visible]').checked })));
  // the view's pixels: how many are far from the background colour
  const painted = () => page.evaluate(() => {
    const R = BM.app.R(), gl = R.gl;
    R.render();   // (read in the same task: the drawing buffer is not preserved)
    const w = gl.drawingBufferWidth, h = gl.drawingBufferHeight, px = new Uint8Array(w * h * 4);
    gl.bindFramebuffer(gl.FRAMEBUFFER, null);
    gl.readPixels(0, 0, w, h, gl.RGBA, gl.UNSIGNED_BYTE, px);
    const bg = [px[0], px[1], px[2]];
    let n = 0, red = 0, blue = 0;
    for (let i = 0; i < px.length; i += 16) {
      if (Math.abs(px[i] - bg[0]) + Math.abs(px[i + 1] - bg[1]) + Math.abs(px[i + 2] - bg[2]) > 60) n++;
      if (px[i] > 150 && px[i + 1] < 110 && px[i + 2] < 110) red++;
      if (px[i + 2] > 150 && px[i] < 110) blue++;
    }
    return { n, red, blue };
  });

  // ------------------------------------------------------------ a restricted fchk file
  await open('h2o_631gd_pure.fchk', fs.readFileSync(path.join(ORB, 'h2o_631gd_pure.fchk'), 'utf8'));
  const tabs = await page.$$eval('#tabs button', b => b.map(x => x.dataset.tab));
  check(tabs.includes('surf'), 'Surfaces tab for an fchk file: ' + tabs.join(','));
  const info = await page.$eval('#panel', p => p.innerText);
  check(/Basis set and orbitals/.test(info), 'Info tab lists the basis set and orbitals');
  await tab('surf');
  const sel = await page.$eval('[data-orb=index]', s => ({ v: s.value, text: s.options[s.selectedIndex].text, n: s.options.length }));
  check(sel.v === '4' && /HOMO/.test(sel.text) && sel.n === 18, 'orbital menu: HOMO selected, 18 orbitals: ' + JSON.stringify(sel));
  await page.waitForFunction(() => /orthonormal/.test(document.querySelector('#sf-orb').innerText), null, { timeout: 5000 }).catch(() => {});
  check(/orthonormal to \de-\d+ \(the standard convention\)/.test(await page.$eval('#sf-orb', e => e.innerText)), 'the orthonormality check is reported');
  await page.click('[data-act=sf-orb]');
  await idle('#sf-orb');
  let cs = await cards();
  check(cs.length === 1 && cs[0].name === 'MO 5 (HOMO)' && cs[0].kind === 'orbital' && cs[0].on, 'HOMO surface: ' + JSON.stringify(cs));
  const st = await page.$eval('#status-msg', e => e.textContent);
  check(/MO 5 \(HOMO\): \d+ × \d+ × \d+ points in [\d.]+ s \(\d+ threads?\)/.test(st), 'computed in workers: ' + st);
  // the workers computed what the page computes
  const cmp = await page.evaluate(() => {
    const d = BM.app.state.active, s = d.surfaces[0], M = BM.orbitals.model(d.data), v = s.vol;
    const B = 0.529177210903, grid = { origin: v.origin.map(x => x / B), axes: v.axes.map(a => a.map(x => x / B)), dims: v.dims };
    const ref = BM.orbitals.evaluate(M, { kind: 'mo', spin: 'alpha', index: 4 }, grid);
    let m = 0; for (let i = 0; i < ref.length; i++) m = Math.max(m, Math.abs(ref[i] - v.data[i]));
    return { m, n: ref.length };
  });
  check(cmp.n > 1000 && cmp.m < 1e-6, 'worker grid = page grid: ' + JSON.stringify(cmp));
  await page.waitForTimeout(300);
  const p1 = await painted();
  check(p1.n > 2000 && p1.blue > 200, 'the orbital is drawn: ' + JSON.stringify(p1));
  // LUMO: a second card, the HOMO hidden; HOMO again: reused, not recomputed
  await page.click('[data-act=sf-lumo]');
  await idle('#sf-orb');
  cs = await cards();
  check(cs.length === 2 && cs[1].name === 'MO 6 (LUMO)' && cs[1].on && !cs[0].on, 'LUMO shown alone: ' + JSON.stringify(cs));
  await page.click('[data-act=sf-homo]');
  await idle('#sf-orb');
  cs = await cards();
  check(cs.length === 2 && cs[0].on && !cs[1].on, 'HOMO again from the list: ' + JSON.stringify(cs));
  // density, coloured by the ESP of the Mulliken charges
  await page.selectOption('[data-orb=kind]', 'density');
  check(!(await page.$('[data-orb=index]')), 'no orbital menu for the density');
  await page.click('[data-act=sf-orb]');
  await idle('#sf-orb');
  cs = await cards();
  check(cs.length === 3 && cs[2].name === 'Electron density' && cs[2].kind === 'density', 'density card: ' + JSON.stringify(cs));
  const opts = await page.$$eval('.surf[data-surf="2"] [data-s=esp] option', o => o.map(x => x.value));
  check(opts.join() === 'none,q:Mulliken', 'colour sources: ' + opts.join());
  await page.selectOption('.surf[data-surf="2"] [data-s=esp]', 'q:Mulliken');
  await page.waitForTimeout(400);
  const esp = await page.$eval('.surf[data-surf="2"]', c => ({ range: c.querySelector('[data-s=espRange]') && c.querySelector('[data-s=espRange]').placeholder, scale: !!c.querySelector('.esp-scale'), note: c.innerText }));
  check(esp.scale && +esp.range > 0 && /kcal\/mol/.test(esp.note), 'ESP scale and range: ' + JSON.stringify(esp.range));
  const p2 = await painted();
  check(p2.red > 100 && p2.blue > 100, 'ESP colours drawn (red and blue): ' + JSON.stringify(p2));
  await page.fill('.surf[data-surf="2"] [data-s=espRange]', '0.2');
  await page.dispatchEvent('.surf[data-surf="2"] [data-s=espRange]', 'input');
  await page.waitForTimeout(300);
  check(await page.evaluate(() => BM.app.state.active.surfaces[2].esp.range === 0.2), 'ESP range set');
  // molecular surfaces
  for (const [type, re] of [['vdw', /van der Waals/], ['sas', /accessible/], ['ses', /excluded/]]) {
    await page.selectOption('[data-ms=type]', type);
    await page.click('[data-act=sf-ms]');
    await idle('#sf-ms');
    cs = await cards();
    const last = cs[cs.length - 1];
    check(re.test(last.name) && last.kind === 'surface' && last.on, type + ' surface: ' + JSON.stringify(last));
  }
  check(!(await page.$('.surf:last-child [data-s=isoSlider]')), 'no isovalue slider on a molecular surface');
  // remove a computed surface
  const n0 = (await cards()).length;
  await page.click('.surf:last-child [data-act=sf-del]');
  check((await cards()).length === n0 - 1, 'surface removed');
  // Save as Molden: the orbitals are written
  const saved = await page.evaluate(() => {
    const d = BM.app.state.active, c = BM.app.saveContext(d, false);
    return BM.qcWriters.molden(c);
  });
  check(/\[GTO\]/.test(saved) && /\[5D\]/.test(saved) && (saved.match(/ Sym=/g) || []).length === 18, 'Molden export with [GTO], [5D] and 18 orbitals');

  // ------------------------------------------------------------ an unrestricted Molden file, with a hostile label
  const uhf = fs.readFileSync(path.join(ORB, 'nh2_uhf_631g.molden'), 'utf8').replace(' Sym=       2a', ' Sym= <img src=x onerror=alert(1)>');
  await open('nh2_uhf.molden', uhf);
  await tab('surf');
  check(!!(await page.$('[data-orb=spin]')), 'α/β menu for unrestricted orbitals');
  const kinds = await page.$$eval('[data-orb=kind] option', o => o.map(x => x.value));
  check(kinds.join() === 'mo,density,spin', 'spin density offered: ' + kinds.join());
  const optText = await page.$$eval('[data-orb=index] option', o => o.map(x => x.textContent));
  check(optText.some(t => /<img src=x onerror=alert\(1\)>/.test(t)), 'a file\'s symmetry label shown as text');
  check(!(await page.$('#panel img')), 'no element made from a file\'s text');
  await page.selectOption('[data-orb=spin]', 'beta');
  const bsel = await page.$eval('[data-orb=index]', s => s.options[s.selectedIndex].text);
  await page.selectOption('[data-orb=kind]', 'spin');
  await page.click('[data-act=sf-orb]');
  await idle('#sf-orb');
  cs = await cards();
  check(cs.length === 1 && cs[0].kind === 'spin density', 'spin density card: ' + JSON.stringify(cs) + ' (beta menu: ' + bsel + ')');
  const both = await page.evaluate(() => BM.app.state.active.surfaces[0].both);
  check(both === true, 'spin density shows both signs');
  // MO level plot: pick a level, Show orbital
  await page.evaluate(() => { const d = BM.app.state.active; d.plotType = 'mo'; BM.app.renderPlot(); });
  await page.waitForSelector('#plot-controls [data-show-mo]');
  check(await page.$eval('#plot-controls [data-show-mo]', b => b.disabled), '"Show orbital" waits for a level');
  await page.evaluate(() => { BM.app.plot().onPick(100000 + 3); });
  await page.click('#plot-controls [data-show-mo]');
  await idle('#sf-orb');
  cs = await cards();
  check(cs.length === 2 && /^MO 3 β/.test(cs[1].name) && cs[1].on, 'Show orbital computes the picked level: ' + JSON.stringify(cs));

  // ------------------------------------------------------------ a density cube coloured by a potential cube (two documents)
  const cube = (title, f) => {
    const n = 12, h = 0.6, L = [' ' + title, ' made for the test', '    1   -3.300000   -3.300000   -3.300000'];
    for (let a = 0; a < 3; a++) L.push('   ' + n + '   ' + [0, 1, 2].map(b => (a === b ? h : 0).toFixed(6)).join('   '));
    L.push('    8    8.000000    0.000000    0.000000    0.000000');
    const v = [];
    for (let i = 0; i < n; i++) for (let j = 0; j < n; j++) for (let k = 0; k < n; k++) v.push(f(-3.3 + i * h, -3.3 + j * h, -3.3 + k * h).toExponential(5));
    for (let k = 0; k < v.length; k += 6) L.push(v.slice(k, k + 6).join(' '));
    return L.join('\n') + '\n';
  };
  await open('pot.cube', cube('Electrostatic potential from Total SCF Density', (x, y, z) => 0.05 * x));
  await open('dens.cube', cube('Electron density from Total SCF Density', (x, y, z) => Math.exp(-(x * x + y * y + z * z) / 2)));
  await tab('surf');
  const srcs = await page.$$eval('.surf [data-s=esp] option', o => o.map(x => x.textContent));
  check(srcs.some(t => /Potential: .*\(pot\.cube\)/.test(t)), 'the potential of the other document is offered: ' + srcs.join(' | '));
  const v = await page.$$eval('.surf [data-s=esp] option', o => o.map(x => x.value).find(x => x.startsWith('v:')));
  await page.selectOption('.surf [data-s=esp]', v);
  await page.waitForTimeout(400);
  const p3 = await painted();
  check(p3.red > 100 && p3.blue > 100, 'density coloured by the potential cube: ' + JSON.stringify(p3));
  const msOnly = await page.evaluate(() => !BM.orbitals.model(BM.app.state.active.data));
  check(msOnly && !(await page.$('#sf-orb')) && !!(await page.$('#sf-ms')), 'a cube file has molecular surfaces but no orbital section');

  // ------------------------------------------------------------ a structure file: only molecular surfaces
  await open('w.xyz', '3\nwater\nO 0 0 0.12\nH 0 0.76 -0.47\nH 0 -0.76 -0.47\n');
  await tab('surf');
  check(!!(await page.$('#sf-ms')) && !(await page.$('#sf-orb')), 'XYZ: molecular surface section only');
  check(/fchk or Molden/.test(await page.$eval('#panel', p => p.innerText)), 'a hint where orbitals come from');

  check(!errors.length, 'no page errors: ' + errors.join(' | '));
  await browser.close();
  console.log(`Orbital and surface browser tests: ${passes} passed, ${failures} failed`);
  process.exit(failures ? 1 : 0);
})().catch(e => { console.error(e); process.exit(1); });
