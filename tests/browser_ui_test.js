// Browser checks of the 1.4.2 interface fixes in dist/BlueMolecule.html (headless Chromium):
//   #1 an output that is still being written: no bogus energy, "Not finished" and "in progress" shown
//   #2 drop-down menus stay open after a click, and work by keyboard and by press-drag-release
//   #3 the Job name of the calculation-input dialog is the saved file name
//   #4 ↑ / ↓ (and Page Up / Page Down, Home / End) move through the Steps, Vibrations and States lists
//   #5 right-click deletes an atom (and, in the Draw tool, a bond) while building; right-drag still pans
// Usage: node tests/browser_ui_test.js
const { chromium } = require('playwright');
const path = require('path');
const fs = require('fs');

const ROOT = path.resolve(__dirname, '..');
const HTML = 'file://' + path.join(ROOT, 'dist', 'BlueMolecule.html');
const FIX = path.join(ROOT, 'tests', 'fixtures');

let failures = 0, passes = 0;
const check = (cond, msg) => { if (cond) passes++; else { failures++; console.log('  FAIL:', msg); } };
const sleep = ms => new Promise(r => setTimeout(r, ms));

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
    window.__bmSaveHook = async (name, blob) => { window.__saved.push({ name, text: await blob.text() }); };
    window.__ctx = [];   // was the browser's context menu suppressed over the 3D view?
    document.addEventListener('contextmenu', e => window.__ctx.push({ prevented: e.defaultPrevented, canvas: e.target.id === 'gl' || !!(e.target.closest && e.target.closest('#viewport, #gl')) }));
  });
  const open = async file => {
    await page.setInputFiles('#file-input', path.join(FIX, file));
    await page.waitForFunction(n => BM.app.state.active && BM.app.state.active.name === n, file);
    await sleep(150);
  };
  const openText = async (name, text) => {
    await page.evaluate(([n, t]) => BM.app.openText(n, t), [name, text]);
    await page.waitForFunction(n => BM.app.state.active && BM.app.state.active.name === n, name);
    await sleep(150);
  };
  const tab = async id => { await page.click('#tabs [data-tab="' + id + '"]'); await sleep(80); };
  const menu = () => page.evaluate(() => {
    const m = document.querySelector('.bm-menu');
    if (!m) return null;
    const items = [...m.querySelectorAll('.bm-opt')];
    return { n: items.length, groups: m.querySelectorAll('.bm-grp').length, hi: items.findIndex(x => x.classList.contains('hi')),
      labels: items.map(x => x.textContent), off: items.map(x => x.classList.contains('off')), sel: BM.selectMenu.current() && BM.selectMenu.current().id };
  });
  const valueOf = sel => page.$eval(sel, e => e.value);
  const focused = () => page.evaluate(() => document.activeElement && (document.activeElement.id || document.activeElement.tagName));

  // ================================================================ #2 drop-down menus
  console.log('#2 drop-down menus');
  await open('g16_water_optfreq.log');
  /** Centre of the plot picker now (its width follows the chosen plot's name). */
  const ppc = async () => { const b = await page.locator('#plot-pick').boundingBox(); return [b.x + b.width / 2, b.y + b.height / 2]; };
  const nOpts = await page.$eval('#plot-pick', s => [...s.options].length);
  const v0 = await valueOf('#plot-pick');
  // a plain click opens the menu and it stays open
  await page.mouse.click(...(await ppc()));
  await sleep(400);
  let m = await menu();
  check(m && m.sel === 'plot-pick' && m.n === nOpts, 'a click opens the plot menu with every option: ' + JSON.stringify(m && { n: m.n, sel: m.sel }) + ' of ' + nOpts);
  check(m && m.labels[m.hi] === (await page.$eval('#plot-pick', s => s.options[s.selectedIndex].textContent)), 'the current plot is highlighted');
  // click an option that is in view and differs from the current one
  const target = await page.evaluate(() => {
    const menuEl = document.querySelector('.bm-menu'), r = menuEl.getBoundingClientRect(), sel = document.querySelector('#plot-pick');
    const items = [...menuEl.querySelectorAll('.bm-opt')];
    const k = items.findIndex(x => !x.classList.contains('off') && x._opt.value !== sel.value && x.getBoundingClientRect().bottom <= r.bottom && x.getBoundingClientRect().top >= r.top);
    const b = items[k].getBoundingClientRect();
    return { k, value: items[k]._opt.value, x: b.left + b.width / 2, y: b.top + b.height / 2 };
  });
  await page.mouse.click(target.x, target.y);
  await sleep(150);
  check(!(await menu()), 'the menu closes after a choice');
  check((await valueOf('#plot-pick')) === target.value && target.value !== v0, 'clicking an option chooses it: ' + target.value);
  check((await focused()) !== 'plot-pick', 'after a mouse choice the keyboard goes back to the window (focus: ' + (await focused()) + ')');
  // keyboard: Enter opens, ↓ moves, Enter chooses; focus stays on the drop-down
  await page.focus('#plot-pick');
  const before = await page.$eval('#plot-pick', s => s.selectedIndex);
  await page.keyboard.press('Enter');
  await sleep(80);
  m = await menu();
  check(!!m, 'Enter opens the menu of a focused drop-down');
  await page.keyboard.press('ArrowDown');
  await page.keyboard.press('Enter');
  await sleep(120);
  const afterKey = await page.$eval('#plot-pick', s => ({ i: s.selectedIndex, n: s.options.length, off: [...s.options].map(o => o.disabled) }));
  let expect = before;
  for (let k = before + 1; k < afterKey.n; k++) if (!afterKey.off[k]) { expect = k; break; }
  check(afterKey.i === expect, 'ArrowDown + Enter chooses the next option: ' + before + ' → ' + afterKey.i + ' (expected ' + expect + ')');
  check(!(await menu()), 'after a keyboard choice the menu closes');
  // (the plot picker gives the keyboard back to the window after every change; other drop-downs keep it)
  await page.focus('#sel-labels');
  const lab0 = await valueOf('#sel-labels');
  await page.keyboard.press('Enter');
  await page.keyboard.press('ArrowDown');
  await page.keyboard.press('Enter');
  await sleep(100);
  check((await valueOf('#sel-labels')) !== lab0 && (await focused()) === 'sel-labels', 'a keyboard choice in the Labels drop-down keeps the focus there');
  await page.selectOption('#sel-labels', 'none');
  // Escape closes without changing
  const vEsc = await valueOf('#plot-pick');
  await page.focus('#plot-pick');
  await page.keyboard.press('Space');
  await sleep(60);
  check(!!(await menu()), 'Space opens the menu');
  await page.keyboard.press('ArrowDown');
  await page.keyboard.press('Escape');
  await sleep(60);
  check(!(await menu()) && (await valueOf('#plot-pick')) === vEsc, 'Escape closes the menu and keeps the value');
  // type-ahead
  await page.keyboard.press('Enter');
  await sleep(60);
  m = await menu();
  const letterAt = m.labels.findIndex((t, k) => k !== m.hi && !m.off[k] && /^[a-z]/i.test(t.trim()));
  const letter = m.labels[letterAt].trim()[0].toLowerCase();
  await page.keyboard.press(letter);
  await sleep(40);
  m = await menu();
  check(m && m.labels[m.hi].trim().toLowerCase().startsWith(letter), 'typing "' + letter + '" jumps to an option that starts with it: ' + (m && m.labels[m.hi]));
  await page.keyboard.press('Escape');
  // press, drag to an option and release, as with a system menu
  const vDrag = await valueOf('#plot-pick');
  await page.mouse.move(...(await ppc()));
  await page.mouse.down();
  await sleep(60);
  const dragTo = await page.evaluate(() => {
    const menuEl = document.querySelector('.bm-menu'), r = menuEl.getBoundingClientRect(), sel = document.querySelector('#plot-pick');
    const items = [...menuEl.querySelectorAll('.bm-opt')];
    const ok = items.filter(x => !x.classList.contains('off') && x._opt.value !== sel.value && x.getBoundingClientRect().bottom <= r.bottom && x.getBoundingClientRect().top >= r.top);
    const t = ok[ok.length - 1], b = t.getBoundingClientRect();
    return { value: t._opt.value, x: b.left + b.width / 2, y: b.top + b.height / 2 };
  });
  await page.mouse.move(dragTo.x, dragTo.y, { steps: 8 });
  await page.mouse.up();
  await sleep(150);
  check(!(await menu()) && (await valueOf('#plot-pick')) === dragTo.value && dragTo.value !== vDrag, 'press, drag and release chooses the option under the pointer: ' + dragTo.value);
  // a click elsewhere closes the menu without a change
  await page.mouse.click(...(await ppc()));
  await sleep(120);
  check(!!(await menu()), 'the menu opens again');
  const vOut = await valueOf('#plot-pick');
  const gl = await page.locator('#gl').boundingBox();
  await page.mouse.click(gl.x + 30, gl.y + gl.height - 30);
  await sleep(120);
  check(!(await menu()) && (await valueOf('#plot-pick')) === vOut, 'a click outside closes the menu and keeps the value');
  // a second click on the drop-down closes its menu
  await page.mouse.click(...(await ppc()));
  await sleep(380);
  await page.mouse.click(...(await ppc()));
  await sleep(120);
  check(!(await menu()), 'a second click on the drop-down closes its menu');
  // a toolbar drop-down
  const lb = await page.locator('#sel-labels').boundingBox();
  await page.mouse.click(lb.x + lb.width / 2, lb.y + lb.height / 2);
  await sleep(380);
  m = await menu();
  check(m && m.sel === 'sel-labels' && m.n === 5, 'the Labels drop-down of the toolbar opens its menu: ' + JSON.stringify(m && m.labels));
  await page.locator('.bm-menu .bm-opt', { hasText: 'Number' }).first().click();
  await sleep(120);
  check((await page.evaluate(() => BM.app.state.view.labels)) === 'number', 'choosing "Number" sets the atom labels');
  await page.keyboard.press('l');   // shortcuts work right after a mouse choice (L cycles the labels)
  await sleep(80);
  check((await page.evaluate(() => BM.app.state.view.labels)) === 'symnum', 'the L shortcut works right after a mouse choice in a drop-down');
  await page.selectOption('#sel-labels', 'none');   // a value set by code (as tests and scripts do) still applies
  check((await page.evaluate(() => BM.app.state.view.labels)) === 'none', 'setting the value of the <select> directly still works');
  // grouped options (the calculation dialog: Program has no groups, Basis may; the toolbar has none) and dialogs
  await page.evaluate(() => BM.app.exportCalc('orca'));
  await page.waitForSelector('#dlg[open]');
  const prog = await page.locator('#dlg [data-c="program"]').boundingBox();
  await page.mouse.click(prog.x + prog.width / 2, prog.y + prog.height / 2);
  await sleep(380);
  m = await menu();
  const nProg = await page.evaluate(() => BM.inputPrograms.length);
  check(m && m.n === nProg, 'the Program drop-down in the dialog opens its menu with ' + nProg + ' programs: ' + (m && m.n));
  await page.locator('.bm-menu .bm-opt', { hasText: /^Gaussian/ }).first().click();
  await sleep(120);
  check((await valueOf('#dlg [data-c="program"]')) === 'gaussian', 'choosing a program in the dialog works');
  check((await focused()) !== 'BODY', 'inside a dialog the focus stays in the dialog after a mouse choice');
  await page.keyboard.press('Escape');   // closes the dialog
  await sleep(120);

  // ================================================================ #3 Job name = file name
  console.log('#3 Job name');
  await page.evaluate(() => BM.app.exportCalc('orca'));
  await page.waitForSelector('#dlg[open]');
  const jn = await page.inputValue('#dlg [data-c="title"]');
  check(jn === 'g16_water_optfreq', 'the Job name of a file document defaults to the file name without extension: ' + jn);
  let note = await page.textContent('#dlg [data-c="note"]');
  check(/Saves g16_water_optfreq\.inp\./.test(note), 'the note names the file to be saved: ' + note);
  await page.fill('#dlg [data-c="title"]', 'Água: opt 1');
  await sleep(80);
  note = await page.textContent('#dlg [data-c="note"]');
  check(/Saves Agua_opt_1\.inp\./.test(note), 'the note follows the Job name: ' + note);
  await page.evaluate(() => { window.__saved.length = 0; });
  await page.click('#dlg .btn.primary');
  await sleep(150);
  let saved = await page.evaluate(() => window.__saved.map(x => x.name));
  check(saved.length === 1 && saved[0] === 'Agua_opt_1.inp', 'ORCA input saved as <Job name>.inp: ' + saved);
  await page.evaluate(() => BM.app.exportCalc('gaussian'));
  await page.waitForSelector('#dlg[open]');
  await page.fill('#dlg [data-c="title"]', 'water scan.gjf');
  await sleep(80);
  const gpre = await page.inputValue('#dlg [data-c="preview"]');
  check(/%chk=water_scan\.chk/i.test(gpre), 'the Gaussian checkpoint takes the Job name: ' + (gpre.match(/%chk=\S+/i) || ['none'])[0]);
  await page.evaluate(() => { window.__saved.length = 0; });
  await page.click('#dlg .btn.primary');
  await sleep(150);
  saved = await page.evaluate(() => window.__saved.map(x => x.name));
  check(saved.length === 1 && saved[0] === 'water_scan.gjf', 'a typed extension is not doubled: ' + saved);
  await page.evaluate(() => BM.app.exportCalc('orca'));
  await page.waitForSelector('#dlg[open]');
  await page.fill('#dlg [data-c="title"]', '   ');
  await sleep(80);
  await page.evaluate(() => { window.__saved.length = 0; });
  await page.click('#dlg .btn.primary');
  await sleep(150);
  saved = await page.evaluate(() => window.__saved.map(x => x.name));
  check(saved.length === 1 && saved[0] === 'job.inp', 'an empty Job name saves job.inp: ' + saved);
  // every program, with spaces and special characters in the Job name
  const progs = await page.evaluate(() => BM.inputPrograms.map(p => ({ id: p.id, ext: p.ext })));
  const wrong = [];
  for (const p of progs) {
    await page.evaluate(id => BM.app.exportCalc(id), p.id);
    await page.waitForSelector('#dlg[open]');
    await page.fill('#dlg [data-c="title"]', 'My job (test) #1 – é');
    await sleep(60);
    await page.evaluate(() => { window.__saved.length = 0; });
    await page.click('#dlg .btn.primary');
    await sleep(200);
    const names = await page.evaluate(() => window.__saved.map(x => x.name));
    const want = p.id === 'turbomole' ? 'coord' : 'My_job_test_1_e.' + (p.ext || 'inp');   // Turbomole reads a file named coord
    if (names[0] !== want) wrong.push(p.id + ': ' + names.join(',') + ' (expected ' + want + ')');
    if (await page.$('#dlg[open]')) { await page.keyboard.press('Escape'); await sleep(60); }
  }
  check(wrong.length === 0 && progs.length === 10, 'all 10 programs save <Job name>.<extension>: ' + (wrong.join(' | ') || progs.map(p => p.id).join(', ')));

  // ================================================================ #4 arrow keys in the lists
  console.log('#4 arrow keys in the lists');
  const fi = () => page.evaluate(() => BM.app.state.active.frameIndex);
  await tab('steps');
  const nFr = await page.evaluate(() => BM.app.state.active.data.frames.length);
  await page.keyboard.press('Home');
  check((await fi()) === 0, 'Home: first geometry');
  await page.keyboard.press('ArrowDown');
  check((await fi()) === 1, 'ArrowDown: next geometry');
  await page.keyboard.press('ArrowDown');
  await page.keyboard.press('ArrowUp');
  check((await fi()) === 1, 'ArrowUp: previous geometry');
  check(await page.evaluate(() => [...document.querySelectorAll('#panel tr.row.on')].map(r => r.dataset.i).join() === '1'), 'the Steps list highlights the current row');
  await page.keyboard.press('ArrowRight');
  check((await fi()) === 2, 'ArrowRight still steps');
  await page.keyboard.press('End');
  check((await fi()) === nFr - 1, 'End: last geometry');
  await page.keyboard.press('ArrowDown');
  check((await fi()) === 0, 'ArrowDown on the last geometry wraps to the first, as → does');
  await page.keyboard.press('PageDown');
  check((await fi()) === Math.min(nFr - 1, 10), 'PageDown: 10 geometries on (or the last)');
  await page.keyboard.press('PageUp');
  check((await fi()) === 0, 'PageUp: 10 geometries back (or the first)');
  // Vibrations
  await tab('vib');
  const vi = () => page.evaluate(() => BM.app.state.active.vibIndex);
  const nVib = await page.evaluate(() => BM.app.state.active.data.vibrations.freqs.length);
  await page.keyboard.press('ArrowDown');
  const v1 = await vi();
  await page.keyboard.press('ArrowDown');
  const v2 = await vi();
  check(v1 >= 0 && v2 === (v1 + 1) % nVib, 'ArrowDown moves to the next vibration: ' + v1 + ' → ' + v2);
  await page.keyboard.press('ArrowUp');
  check((await vi()) === v1, 'ArrowUp moves back: ' + (await vi()));
  await page.keyboard.press('End');
  check((await vi()) === nVib - 1, 'End: last vibration');
  await page.keyboard.press('Home');
  check((await vi()) === 0, 'Home: first vibration');
  await page.keyboard.press('ArrowUp');
  check((await vi()) === nVib - 1, 'ArrowUp on the first vibration wraps to the last');
  check(await page.evaluate(() => [...document.querySelectorAll('#vib-table tr.row.on')].map(r => +r.dataset.i).join() === String(BM.app.state.active.vibIndex)), 'the current vibration row is highlighted');
  await page.keyboard.press('Space');   // stop the animation that a selected vibration starts
  await page.keyboard.press('Escape');
  // States
  await open('g16_water_tddft.log');
  await tab('states');
  const si = () => page.evaluate(() => BM.app.state.active.stateIndex);
  const nSt = await page.evaluate(() => BM.app.state.active.data.excitedStates.length);
  await page.keyboard.press('ArrowDown');
  const s1 = await si();
  await page.keyboard.press('ArrowDown');
  check(s1 >= 0 && (await si()) === (s1 + 1) % nSt, 'ArrowDown moves to the next excited state: ' + s1 + ' → ' + (await si()));
  await page.keyboard.press('End');
  check((await si()) === nSt - 1, 'End: last excited state');
  await page.keyboard.press('ArrowDown');
  check((await si()) === 0, 'ArrowDown on the last state wraps to the first');
  // a long list keeps the current row in view
  const traj = [];
  for (let k = 0; k < 80; k++) traj.push('3\nE = ' + (-76.4 - k * 1e-4).toFixed(6) + '\nO 0 0 ' + (0.117 + k * 1e-3).toFixed(4) + '\nH 0 0.757 -0.469\nH 0 -0.757 -0.469');
  await openText('long.xyz', traj.join('\n') + '\n');
  await tab('steps');
  const rowVisible = i => page.evaluate(i => {
    const p = document.querySelector('#panel').getBoundingClientRect(), r = document.querySelector('#panel tr.row[data-i="' + i + '"]');
    if (!r) return false;
    const b = r.getBoundingClientRect();
    return b.top >= p.top - 1 && b.bottom <= p.bottom + 1;
  }, i);
  await page.keyboard.press('End');
  await sleep(60);
  check((await fi()) === 79 && await rowVisible(79), 'End in a long Steps list shows the last row');
  await page.keyboard.press('PageUp');
  await sleep(60);
  check((await fi()) === 69 && await rowVisible(69), 'PageUp moves 10 rows back and keeps the row in view');
  await page.keyboard.press('Home');
  await sleep(60);
  check((await fi()) === 0 && await rowVisible(0), 'Home shows the first row again');
  for (let k = 0; k < 40; k++) await page.keyboard.press('ArrowDown');
  await sleep(60);
  check((await fi()) === 40 && await rowVisible(40), 'ArrowDown 40 times: row 41 is current and in view');
  // the keys are left to text fields
  await tab('info');
  const fBefore = await fi();
  await page.evaluate(() => { const i = document.createElement('input'); i.id = 'tmp-input'; document.body.appendChild(i); i.focus(); });
  await page.keyboard.press('ArrowDown');
  check((await fi()) === fBefore, 'arrow keys typed in a text field do not move');
  await page.evaluate(() => document.getElementById('tmp-input').remove());
  await page.click('#btn-help');
  await page.waitForSelector('#dlg[open]');
  const help = await page.textContent('#dlg');
  check(/↑ \/ ↓/.test(help) && /Right-click atom/.test(help), 'the Help dialog lists ↑ / ↓ and right-click');
  await page.keyboard.press('Escape');
  await sleep(80);

  // ================================================================ #1 a job that is still running
  console.log('#1 running job');
  const g = fs.readFileSync(path.join(FIX, 'g16_water_optfreq.log'), 'utf8');
  const at = g.lastIndexOf('SCF Done:');
  await openText('running.log', g.slice(0, g.indexOf('=', at) + 4));   // "... SCF Done:  E(RB3LYP) =  -7"
  await tab('info');
  const run = await page.evaluate(() => {
    const d = BM.app.state.active, fr = d.data.frames;
    const c = BM.analysis.convergenceSeries(d.data, {}), e = c.series.find(s => s.key === 'energy');
    return { n: fr.length, energies: fr.map(f => f.energy), series: e ? e.values : null, running: !!d.data.meta.running, info: document.querySelector('#panel').textContent };
  });
  check(run.running && /Not finished/.test(run.info), 'the Info tab says the job is not finished');
  check(run.energies.every(e => e === null || e < -70), 'no geometry gets a bogus energy: ' + run.energies.join(', '));
  check(!run.series || run.series.every(v => v === null || v < -70), 'the energy plot has no bogus point');
  await tab('steps');
  const busy = await page.evaluate(() => [...document.querySelectorAll('#panel .lab.busy')].map(x => x.closest('tr').dataset.i));
  check(run.energies[run.n - 1] === null ? busy.length === 1 && +busy[0] === run.n - 1 : busy.length === 0, 'the step without an energy yet is labelled "in progress": ' + busy);

  // ================================================================ #5 right-click deletes while building
  console.log('#5 right-click delete');
  await page.click('#btn-new');
  await page.waitForFunction(() => BM.app.state.active && /^New molecule/.test(BM.app.state.active.name));
  await page.keyboard.press('d');
  const box = await page.locator('#gl').boundingBox();
  const info = () => page.evaluate(() => { const d = BM.app.state.active, f = d.data.frames[d.frameIndex]; return { n: f.Z.length, formula: BM.geom.formula(f.Z), bonds: f.bonds ? f.bonds.length : 0, tool: BM.app.state.build.tool }; });
  const atomXY = i => page.evaluate(i => { const R = BM.app.R(), x = R.xyz; const s = R.worldToScreen([x[3 * i], x[3 * i + 1], x[3 * i + 2]]); return [s.x, s.y]; }, i).then(([x, y]) => [box.x + x, box.y + y]);
  const pickAt = (x, y) => page.evaluate(([x, y]) => BM.app.R().pickAny(x, y), [x - box.x, y - box.y]);
  /** Page position of atom i where a click picks that atom. */
  const atomSpot = async i => { const [x, y] = await atomXY(i); const h = await pickAt(x, y); return h.kind === 'atom' && h.index === i ? [x, y] : null; };
  const rightClick = async (x, y) => { await page.mouse.click(x, y, { button: 'right' }); await sleep(160); };
  const hIndex = () => page.evaluate(() => { const f = BM.app.curFrame(BM.app.state.active); return f.Z.indexOf(1); });
  await page.mouse.click(box.x + box.width / 2, box.y + box.height / 2);
  await sleep(150);
  let s = await info();
  check(s.formula === 'CH4', 'Draw: a click places methane: ' + s.formula);
  await page.keyboard.press('d');   // shows the tool hint again at once
  check(/right-click/.test(await page.textContent('#status-msg')), 'the Draw hint mentions right-click: ' + (await page.textContent('#status-msg')));
  let hs = await atomSpot(await hIndex());
  check(!!hs, 'an H atom can be clicked');
  await page.evaluate(() => { window.__ctx.length = 0; });
  if (hs) await rightClick(hs[0], hs[1]);
  s = await info();
  check(s.formula === 'CH3' && s.n === 4, 'Draw: right-click on an H deletes it: ' + s.formula);
  const ctxEv = await page.evaluate(() => window.__ctx.slice());
  check(ctxEv.length >= 1 && ctxEv.every(c => c.prevented && c.canvas), 'the browser context menu is suppressed over the 3D view: ' + JSON.stringify(ctxEv));
  await page.keyboard.press('Control+z');
  await sleep(120);
  check((await info()).formula === 'CH4', 'Ctrl+Z restores the deleted atom');
  // right-drag pans and deletes nothing
  const c0 = await page.evaluate(() => BM.app.R().center.slice());
  const [cx, cy] = await atomXY(0);
  await page.mouse.move(cx, cy);
  await page.mouse.down({ button: 'right' });
  await page.mouse.move(cx + 80, cy + 40, { steps: 8 });
  await page.mouse.up({ button: 'right' });
  await sleep(150);
  const c1 = await page.evaluate(() => BM.app.R().center.slice());
  s = await info();
  check(Math.hypot(c1[0] - c0[0], c1[1] - c0[1], c1[2] - c0[2]) > 0.1 && s.formula === 'CH4', 'right-drag from an atom pans the view and deletes nothing: ' + s.formula);
  // a bond: grow a C–C bond, then right-click the bond
  const [ax, ay] = await atomXY(0);
  await page.mouse.move(ax, ay); await page.mouse.down(); await page.mouse.move(ax + 150, ay - 60, { steps: 8 }); await page.mouse.up();
  await sleep(150);
  s = await info();
  check(s.formula === 'C2H6', 'Draw: dragging from C builds ethane: ' + s.formula);
  const cc = await page.evaluate(() => {
    const R = BM.app.R(), b = R.bonds, x = R.xyz, Z = R.Z;
    for (let k = 0; 2 * k + 1 < b.length; k++) {
      const i = b[2 * k], j = b[2 * k + 1];
      if (Z[i] !== 6 || Z[j] !== 6) continue;
      const p = R.worldToScreen([(x[3 * i] + x[3 * j]) / 2, (x[3 * i + 1] + x[3 * j + 1]) / 2, (x[3 * i + 2] + x[3 * j + 2]) / 2]);
      const h = R.pickAny(p.x, p.y);
      if (h.kind === 'bond' && h.index === k) return { x: p.x, y: p.y };
    }
    return null;
  });
  check(!!cc, 'the C–C bond can be clicked');
  if (cc) {
    await rightClick(box.x + cc.x, box.y + cc.y);
    s = await info();
    const ccLeft = await page.evaluate(() => { const f = BM.app.curFrame(BM.app.state.active); return f.bonds.some(b => f.Z[b[0]] === 6 && f.Z[b[1]] === 6); });
    check(!ccLeft && s.formula === 'C2H8', 'Draw: right-click on a bond deletes the bond, and ±H caps both carbons: ' + s.formula + ', ' + s.bonds + ' bonds');
    await page.keyboard.press('Control+z');
    await sleep(120);
  }
  // a carbon: deletes it with its hydrogens (±H on)
  s = await info();
  const cSpot = await atomSpot(0);
  if (cSpot) await rightClick(cSpot[0], cSpot[1]);
  const s2 = await info();
  check(s.formula === 'C2H6' && s2.formula === 'CH4', 'Draw: right-click on a carbon deletes it and its H; the other carbon gets back its H: ' + s.formula + ' → ' + s2.formula);
  // the other tools
  const tryTool = async (key, deletes) => {
    await page.keyboard.press(key);
    await sleep(60);
    const t = (await info()).tool;
    const n0 = (await info()).n, h = await hIndex(), spot = h >= 0 ? await atomSpot(h) : null;
    if (!spot) { check(false, t + ': no clickable H'); return; }
    await rightClick(spot[0], spot[1]);
    const n1 = (await info()).n;
    check(deletes ? n1 === n0 - 1 : n1 === n0, t + (deletes ? ': right-click deletes the atom' : ': right-click deletes nothing') + ' (' + n0 + ' → ' + n1 + ')');
    if (deletes && n1 !== n0) { await page.keyboard.press('Control+z'); await sleep(100); }
  };
  await tryTool('t', true);
  await tryTool('k', true);
  await tryTool('m', true);
  await tryTool('n', false);
  await tryTool('s', false);
  // a file geometry: deleting every atom is refused
  await openText('water.xyz', '3\nwater\nO 0 0 0.117\nH 0 0.757 -0.469\nH 0 -0.757 -0.469\n');
  await page.keyboard.press('d');
  await sleep(80);
  const oSpot = await atomSpot(0);
  if (oSpot) await rightClick(oSpot[0], oSpot[1]);
  s = await info();
  const toast = await page.textContent('#toast');
  check(s.n === 3 && /Cannot delete every atom/.test(toast), 'a file geometry keeps its last atoms: ' + s.n + ' atoms, ' + JSON.stringify(toast));
  const h1 = await atomSpot(1);
  if (h1) await rightClick(h1[0], h1[1]);
  check((await info()).n === 2, 'right-click deletes an H of a file geometry: ' + (await info()).n);

  // 2.0.1: the Info tab names the level of the geometry's energy and lists the job's other levels apart
  await openText('triplets.out', fs.readFileSync(path.join(FIX, 'review/orca6_water_tddft_triplets.out'), 'utf8'));
  await tab('info');
  const en = await page.$eval('#panel', e => e.textContent.replace(/\s+/g, ' '));
  check(/This geometry\s*-76\.31794234 Eh \(SCF\)/.test(en) && /Excited state\s*-76\.03570334 Eh/.test(en),
    'the Info tab names the energy level and shows the excited state apart: ' + (en.match(/Energy.{0,140}/) || [''])[0]);
  await openText('mp5.log', fs.readFileSync(path.join(FIX, 'review/g16_water_mp5.log'), 'utf8'));
  await tab('info');
  const e5 = await page.$eval('#panel', e => e.textContent.replace(/\s+/g, ' '));
  check(/This geometry\s*-\d+\.\d{8} Eh \(MP5\)/.test(e5) && /SCF\s*-\d+\.\d{8} Eh/.test(e5) && /MP4\s*-\d+\.\d{8} Eh/.test(e5),
    'an MP5 job: the energy is MP5, with SCF, MP2 to MP4 listed: ' + (e5.match(/Energy.{0,240}/) || [''])[0]);

  // ================================================================ 2.1: normal-mode tools for transition states (#14)
  console.log('#14 transition states');
  await page.evaluate(() => { BM.app.state.vib.playing = true; });   // the arrow-key checks above paused the animations
  await open('g16_sn2_ts_optfreq.log');
  await sleep(250);
  let ts = await page.evaluate(() => ({ vib: BM.app.state.active.vibIndex, anim: BM.app.state.anim.kind, tab: BM.app.state.active.tab, toast: document.querySelector('#toast').textContent,
    sum: document.querySelector('#vib-summary') && document.querySelector('#vib-summary').textContent, cls: document.querySelector('#vib-summary') && document.querySelector('#vib-summary').className }));
  check(ts.vib === 0 && ts.anim === 'vib' && ts.tab === 'vib', 'opening a TS frequency job selects and animates its imaginary mode: ' + JSON.stringify(ts));
  check(/one imaginary frequency \(492\.2i cm⁻¹\), selected/.test(ts.toast), 'the message on opening says there is one imaginary frequency: ' + ts.toast);
  check(ts.sum === 'One imaginary frequency (492.2i cm⁻¹): a transition state.' && / ok\b/.test(ts.cls), 'the Vibrations tab says it is a transition state: ' + ts.sum);
  const bondRows = () => page.$$eval('#vib-bonds-table tr.row', rs => rs.map(r => [...r.cells].map(c => c.textContent)));
  let br = await bondRows();
  check(br.length === 2 && br[0].join('|') === 'C1–F2|2.336|1.850|1.364|forms' && br[1].join('|') === 'C1–F3|1.364|1.850|2.336|breaks',
    'the bond table shows the C–F bond that forms and the one that breaks, with the distances at both ends: ' + JSON.stringify(br));
  // amplitude: the distances at the ends follow it
  await page.$eval('#panel [data-vib="amplitude"]', e => { e.value = '0.5'; e.dispatchEvent(new Event('input', { bubbles: true })); });
  await sleep(80);
  br = await bondRows();
  const dm = +br[0][1], d0 = +br[0][2], dp = +br[0][3];
  check(dm > 2.336 && dp < 1.364 && Math.abs((dm - d0) - (d0 - dp)) < 0.002, 'a larger amplitude moves the ends further apart: ' + JSON.stringify(br[0]));
  // arrow length is its own control: the longest arrow is that long, whatever the amplitude
  await page.$eval('#panel [data-vib="arrows"]', e => { e.checked = true; e.dispatchEvent(new Event('input', { bubbles: true })); });
  await page.$eval('#panel [data-vib="arrowLen"]', e => { e.value = '1.5'; e.dispatchEvent(new Event('input', { bubbles: true })); });
  await sleep(60);
  const longest = await page.evaluate(() => { const a = BM.app.R().arrows; let m = 0; for (let i = 0; a && i < a.length; i += 3) m = Math.max(m, Math.hypot(a[i], a[i + 1], a[i + 2])); return m; });
  check(Math.abs(longest - 1.5) < 1e-6, 'the Arrow length slider sets the longest arrow: ' + longest);
  await page.$eval('#panel [data-vib="arrows"]', e => { e.checked = false; e.dispatchEvent(new Event('input', { bubbles: true })); });
  // both ends: a still view with two translucent structures at the amplitude
  await page.check('#panel [data-vib="ends"]');
  await sleep(120);
  ts = await page.evaluate(() => { const d = BM.app.state.active, R = BM.app.R(), eq = d.data.frames[d.data.vibrations.frameIndex].xyz;
    return { anim: BM.app.state.anim.kind, playing: BM.app.state.vib.playing, ends: R.vibEnds.length, op: R.vibEnds.map(m => m.opacity), still: Array.from(R.xyz).every((x, i) => Math.abs(x - eq[i]) < 1e-9) }; });
  check(ts.anim === null && !ts.playing && ts.ends === 2 && ts.op.every(o => o > 0 && o < 1) && ts.still, 'Both ends stops the animation at the centre and draws both ends, translucent: ' + JSON.stringify(ts));
  await page.click('#vib-table tr.row[data-i="1"]');
  await sleep(120);
  ts = await page.evaluate(() => ({ vib: BM.app.state.active.vibIndex, anim: BM.app.state.anim.kind, ends: BM.app.R().vibEnds.length, link: !!document.querySelector('#vib-summary a') }));
  check(ts.vib === 1 && ts.anim === null && ts.ends === 2 && ts.link, 'another mode keeps the still view with its own ends, and the summary offers to show the imaginary one: ' + JSON.stringify(ts));
  await page.click('#vib-summary a');
  await sleep(120);
  check((await page.evaluate(() => BM.app.state.active.vibIndex)) === 0, '"Show it" selects the imaginary mode again');
  await page.uncheck('#panel [data-vib="ends"]');
  await sleep(60);
  check((await page.evaluate(() => BM.app.R().vibEnds.length)) === 0, 'unticking Both ends removes them');
  // a row of the bond table measures the pair
  await page.click('#vib-bonds-table tr.row');
  await sleep(80);
  ts = await page.evaluate(() => ({ sel: BM.app.state.active.selection.join(), chip: document.querySelector('#measure').textContent }));
  check(ts.sel === '0,1' && /C1 – F2/.test(ts.chip) && /distance/.test(ts.chip), 'clicking a bond row measures that distance: ' + JSON.stringify(ts));
  // forward and reverse inputs: two optimisation inputs displaced by -/+ the Displace amount
  await page.fill('#disp-amt', '0.2');
  const taskBefore = await page.evaluate(() => JSON.parse(localStorage.getItem('bluemolecule.calc') || '{}').task);
  await page.click('[data-act="vib-inputs"]');
  await page.waitForSelector('#dlg[open]');
  await sleep(120);
  ts = await page.evaluate(() => ({ title: document.querySelector('#dlg').textContent.includes('Forward and reverse inputs'), route: document.querySelector('#dlg [data-c="route"]').value,
    task: document.querySelector('#dlg [data-c="task"]').value, ro: document.querySelector('#dlg [data-c="preview"]').readOnly, note: document.querySelector('#dlg [data-c="note"]').textContent }));
  check(ts.title && ts.route === '#p opt freq b3lyp/6-31+g(d)' && ts.task === 'opt' && ts.ro, 'the dialog turns the TS route into an optimisation: ' + JSON.stringify(ts));
  check(/_minus\.gjf, g16_sn2_ts_optfreq_plus\.gjf: the geometry displaced by −0\.2 and \+0\.2 Å/.test(ts.note), 'the note names both files and the displacement: ' + ts.note);
  await page.evaluate(() => { window.__saved.length = 0; });
  await page.click('#dlg .btn.primary');
  await sleep(250);
  saved = await page.evaluate(() => window.__saved.map(x => ({ name: x.name, text: x.text })));
  const cz = t => { const l = t.split('\n').find(s => /^C\s/.test(s)); return l ? +l.trim().split(/\s+/)[3] : NaN; };
  check(saved.length === 2 && saved[0].name === 'g16_sn2_ts_optfreq_minus.gjf' && saved[1].name === 'g16_sn2_ts_optfreq_plus.gjf', 'two inputs saved: ' + saved.map(x => x.name));
  check(saved.length === 2 && Math.abs(cz(saved[0].text) + 0.2) < 1e-6 && Math.abs(cz(saved[1].text) - 0.2) < 1e-6 && /%chk=g16_sn2_ts_optfreq_plus\.chk/.test(saved[1].text) && !/\bts\b/i.test(saved[1].text.split('\n').find(s => s[0] === '#')),
    'the carbon (the largest shift) is displaced by -0.2 and +0.2 Å; own checkpoint files; no TS keywords: ' + saved.map(x => cz(x.text)));
  check((await page.evaluate(() => JSON.parse(localStorage.getItem('bluemolecule.calc') || '{}').task)) === taskBefore, 'saving the pair keeps the default Task of the dialog: ' + taskBefore);
  // the number of negative Hessian eigenvalues of the TS search, as a plot and in the Steps tab
  await page.selectOption('#plot-pick', 'opt:negEig');
  await sleep(120);
  ts = await page.evaluate(() => { const s = BM.app.plot().spec; return { int: s.integerY, ys: s.points.map(p => p.y), label: s.yLabel }; });
  check(ts.int && ts.ys.join() === '1,1,1,1' && /Negative Hessian eigenvalues/.test(ts.label), 'the TS search plots one negative eigenvalue at each step: ' + JSON.stringify(ts));
  await tab('steps');
  const stepsText = await page.$eval('#panel', e => e.textContent.replace(/\s+/g, ' '));
  check(/Negative Hessian eigenvalues\s*1\s*✓/.test(stepsText), 'the Steps tab shows the count for the current step: ' + (stepsText.match(/Negative Hessian.{0,30}/) || [''])[0]);
  // a minimum: said so, nothing selected
  await open('g16_water_optfreq.log');
  await sleep(150);
  ts = await page.evaluate(() => ({ vib: BM.app.state.active.vibIndex, sum: document.querySelector('#vib-summary').textContent }));
  check(ts.vib === -1 && ts.sum === 'No imaginary frequencies: a minimum.', 'a minimum is called one and no mode is selected: ' + JSON.stringify(ts));
  // the ORCA keywords of a file are written with spaces
  await open('orca5_water_optfreq.out');
  await page.evaluate(() => BM.app.exportCalc('orca'));
  await page.waitForSelector('#dlg[open]');
  const kw = await page.inputValue('#dlg [data-c="kw"]');
  check(kw && !/,/.test(kw) && /\bOpt\b/i.test(kw), 'the ORCA keywords of the file are separated by spaces: ' + kw);
  await page.keyboard.press('Escape');
  await sleep(80);

  check(errors.length === 0, 'no page errors: ' + errors.slice(0, 5).join(' | '));
  await browser.close();
  console.log(`UI tests (1.4.2 fixes): ${passes} passed, ${failures} failed`);
  process.exit(failures ? 1 : 0);
})().catch(e => { console.error(e); process.exit(1); });
