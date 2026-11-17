// Browser checks of milestone 1.5 in dist/BlueMolecule.html (headless Chromium):
//   #9  zoomable plots: drag a box, the magnifier buttons, double-click, a click still selects, every plot,
//       and the CSV, PNG and SVG exports follow the zoom
//   #10 the QM/MM and ONIOM component energies plot
//   #7  the loading dialog for a large output and a trajectory, "Load other parts…" and the Info tab
//   #8  Amber NetCDF trajectories and restarts after a prmtop (files written by SciPy, when it is installed)
// Usage: node tests/browser_15_test.js
const { chromium } = require('playwright');
const path = require('path');
const fs = require('fs');
const os = require('os');
const { execFileSync } = require('child_process');

const ROOT = path.resolve(__dirname, '..');
const HTML = 'file://' + path.join(ROOT, 'dist', 'BlueMolecule.html');
const FIX = path.join(ROOT, 'tests', 'fixtures');

let failures = 0, passes = 0;
const check = (cond, msg) => { if (cond) passes++; else { failures++; console.log('  FAIL:', msg); } };
const sleep = ms => new Promise(r => setTimeout(r, ms));
const near = (a, b, tol) => Math.abs(a - b) <= tol;

(async () => {
  const browser = await chromium.launch({ args: ['--enable-unsafe-swiftshader'] });
  const ctx = await browser.newContext({ viewport: { width: 1440, height: 1000 }, deviceScaleFactor: 1 });
  const page = await ctx.newPage();
  const errors = [];
  page.on('pageerror', e => errors.push('pageerror: ' + e.message));
  page.on('console', m => { if (m.type() === 'error') errors.push('console: ' + m.text()); });
  await page.goto(HTML);
  await page.waitForFunction(() => window.BM && window.BM.app);
  await page.evaluate(() => {
    window.__saved = [];
    window.__bmSaveHook = async (name, blob) => {
      const u8 = new Uint8Array(await blob.arrayBuffer());
      window.__saved.push({ name, text: /\.png$/.test(name) ? '' : new TextDecoder().decode(u8), bytes: Array.from(u8.slice(0, 64)), size: u8.length, sum: u8.reduce((a, v, i) => (a + v * (i % 251 + 1)) % 1e9, 0) });
    };
    /** The plot's geometry: data coordinates to client pixels, and the visible items that can be picked. */
    window.__plot = () => {
      const p = BM.app.plot(), s = p.spec, c = p.canvas, r = c.getBoundingClientRect(), L = p._layout(c.clientWidth, c.clientHeight), m = p.m;
      const inside = (X, Y) => X > m.l + 3 && X < m.l + L.pw - 3 && Y > m.t + 3 && Y < m.t + L.ph - 3;
      let items = [];
      if (s.type === 'levels') items = p._levelGeom(L).map(g => ({ X: (g.x0 + g.x1) / 2, Y: g.y, it: g.item }));
      else if (s.bars && s.bars.length) items = s.bars.map(b => ({ X: L.X(b.x), Y: m.t + L.ph / 2, it: b }));
      else if (s.lines || (s.points && s.points.length)) items = (s.lines ? [].concat(...s.lines.map(l => l.pts)) : s.points).map(q => ({ X: L.X(q.x), Y: L.Y(q.y), it: q }));
      else items = (s.sticks || []).map(q => ({ X: L.X(q.x), Y: m.t + L.ph / 2, it: q }));
      const vis = items.filter(q => q.it.id !== undefined && q.it.id !== null && isFinite(q.Y) && inside(q.X, q.Y));
      return {
        rect: { x: r.x, y: r.y }, m, pw: L.pw, ph: L.ph, L: { x0: L.x0, x1: L.x1, y0: L.y0, y1: L.y1 }, zoom: p.zoomed(), type: s.type || null, logY: !!s.logY,
        items: vis.map(q => ({ cx: r.x + q.X, cy: r.y + q.Y, id: q.it.id, x: q.it.x, y: q.it.y, sel: s.selected === q.it.id })),
        hidden: items.length - vis.length, xLabels: (s.xTickLabels || []).length
      };
    };
    /** Client pixels of a data point in the plot shown. */
    window.__px = (x, y) => { const p = BM.app.plot(), c = p.canvas, r = c.getBoundingClientRect(), L = p._layout(c.clientWidth, c.clientHeight); return [r.x + L.X(x), r.y + L.Y(y)]; };
  });
  const geom = () => page.evaluate(() => window.__plot());
  const zoom = () => page.evaluate(() => BM.app.plot().zoomed());
  const lastSaved = () => page.evaluate(() => window.__saved[window.__saved.length - 1] || null);
  const save = async act => { const n = await page.evaluate(() => window.__saved.length); await page.click(`[data-act="${act}"]`); await page.waitForFunction(k => window.__saved.length > k, n); return lastSaved(); };
  /** Drags a box between two points given as fractions of the plot area (0,0 top left). */
  const dragF = async (fx0, fy0, fx1, fy1) => {
    const g = await geom(), X = f => g.rect.x + g.m.l + f * g.pw, Y = f => g.rect.y + g.m.t + f * g.ph;
    await page.mouse.move(X(fx0), Y(fy0)); await page.mouse.down();
    await page.mouse.move((X(fx0) + X(fx1)) / 2, (Y(fy0) + Y(fy1)) / 2, { steps: 4 }); await page.mouse.move(X(fx1), Y(fy1), { steps: 4 });
    await page.mouse.up(); await sleep(60);
  };
  const dragPx = async (a, b) => { await page.mouse.move(a[0], a[1]); await page.mouse.down(); await page.mouse.move((a[0] + b[0]) / 2, (a[1] + b[1]) / 2, { steps: 4 }); await page.mouse.move(b[0], b[1], { steps: 4 }); await page.mouse.up(); await sleep(60); };
  const open = async file => {
    await page.setInputFiles('#file-input', path.join(FIX, file));
    const nm = path.basename(file);
    await page.waitForFunction(n => BM.app.state.active && BM.app.state.active.name === n && !document.querySelector('#dlg[open]'), nm);
    await sleep(150);
  };
  const pickPlot = async t => { await page.selectOption('#plot-pick', t); await sleep(80); };
  const doc = f => page.evaluate(f);

  // ================================================================ #9 zoomable plots
  console.log('#9 zoomable plots: IR spectrum');
  await open('g16_water_optfreq.log');
  await pickPlot('ir');
  let g = await geom();
  check(!g.zoom, 'the plot starts unzoomed');
  check(await page.$eval('#btn-plot-fit', b => b.disabled && !b.classList.contains('on')) && await page.$eval('#btn-plot-zoom-out', b => b.disabled), 'Fit and Zoom out are off while the whole plot is shown');
  // a box from 3990 to 3780 cm-1 (IR runs right to left), over the O-H stretches; released on a stick
  const sticks = await doc(() => BM.app.plot().spec.sticks.map(s => ({ x: s.x, id: s.id })));
  const s3833 = sticks.find(s => Math.abs(s.x - 3832.8) < 1), s3940 = sticks.find(s => Math.abs(s.x - 3939.6) < 1);
  const a = await page.evaluate(() => window.__px(3990, 1.1)), b = await page.evaluate(x => window.__px(x, 0), s3940.x);
  await dragPx(a, [await page.evaluate(() => window.__px(3780, 0)).then(p => p[0]), b[1] - 1]);
  let z = await zoom();
  check(z && near(z.x0, 3780, 6) && near(z.x1, 3990, 6), 'drag a box: the x range becomes the box: ' + JSON.stringify(z));
  check(z && z.y0 > -0.05 && z.y1 < 1.15 && z.y1 > 0.9, 'drag a box: the y range becomes the box: ' + JSON.stringify(z && [z.y0, z.y1]));
  check((await doc(() => BM.app.state.active.vibIndex)) === -1, 'the end of a box selects nothing');
  check(await page.$eval('#btn-plot-fit', b => !b.disabled && b.classList.contains('on')), 'the Fit button shows that the plot is zoomed');
  g = await geom();
  check(g.items.length === 2 && g.items.every(it => it.x > 3780 && it.x < 3990), 'only the two sticks in the box are shown: ' + JSON.stringify(g.items.map(i => i.x)));
  // the tick labels are those of the zoomed range
  const svg = await save('plot-svg');
  const nums = [...svg.text.matchAll(/<text[^>]*>(-?[\d.]+)<\/text>/g)].map(m => +m[1]).filter(v => v >= 100);
  check(nums.length >= 3 && nums.every(v => v >= z.x0 - 1 && v <= z.x1 + 1), 'SVG export follows the zoom (x ticks ' + nums.join(', ') + ')');
  const csv = (await save('plot-csv')).text;
  const curve = csv.split('\n\n')[0].split('\n').filter(l => /^\d/.test(l)).map(l => +l.split(',')[0]);
  const stk = csv.split('# sticks\n')[1].split('\n').filter(l => /^\d/.test(l)).map(l => +l.split(',')[0]);
  check(/^# zoomed in: only wavenumbers \(cm-1\) from 37\d\d(\.\d+)? to 39\d\d/m.test(csv), 'the CSV says it holds the zoomed range: ' + (csv.split('\n')[1] || ''));
  check(curve.length > 50 && curve.every(x => x >= z.x0 - 1e-6 && x <= z.x1 + 1e-6), 'CSV curve rows only inside the zoom (' + curve.length + ' rows, ' + curve[0] + '…' + curve[curve.length - 1] + ')');
  check(stk.join() === [s3833.id + 1, s3940.id + 1].sort((p, q) => p - q).join(), 'CSV sticks only inside the zoom: modes ' + stk.join());
  const pngZoomed = await save('plot-png');
  // a click on a stick still selects it, and the zoom stays
  const it = (await geom()).items.find(q => q.id === s3833.id);
  await page.mouse.click(it.cx, it.cy); await sleep(120);
  check((await doc(() => BM.app.state.active.vibIndex)) === s3833.id, 'a click on a stick in the zoomed plot selects its mode');
  const z2 = await zoom();
  check(z2 && near(z2.x0, z.x0, 1e-6) && near(z2.x1, z.x1, 1e-6), 'the zoom is kept when the plot is redrawn after the selection');
  // plot settings: the line width keeps the zoom, the y quantity starts again
  await page.fill('#plot-controls input[data-pc="fwhm"]', '6'); await page.dispatchEvent('#plot-controls input[data-pc="fwhm"]', 'change'); await sleep(80);
  check(!!(await zoom()), 'a new FWHM keeps the zoom');
  // the buttons: in, out, fit
  const w0 = z.x1 - z.x0, c0 = (z.x0 + z.x1) / 2, h0 = z.y1 - z.y0;
  await page.click('#btn-plot-zoom-in'); await sleep(60);
  z = await zoom();
  check(z && near(z.x1 - z.x0, w0 / 2, 1e-6) && near((z.x0 + z.x1) / 2, c0, 1e-6) && near(z.y1 - z.y0, h0 / 2, 1e-9), 'Zoom in halves both ranges around the centre');
  await page.click('#btn-plot-zoom-out'); await sleep(60);
  z = await zoom();
  check(z && near(z.x1 - z.x0, w0, 1e-6), 'Zoom out doubles them again');
  for (let k = 0; k < 8 && (await zoom()); k++) { await page.click('#btn-plot-zoom-out'); await sleep(40); }
  check(!(await zoom()), 'zooming out stops at the whole plot');
  g = await geom();
  check(near(Math.min(g.L.x0, g.L.x1), 400, 1e-6) && near(Math.max(g.L.x0, g.L.x1), 4100, 1e-6), 'and the whole range is shown again: ' + JSON.stringify(g.L));
  check(await page.$eval('#btn-plot-fit', b => b.disabled && !b.classList.contains('on')), 'the Fit button is off again');
  const pngWhole = await save('plot-png');
  check(pngZoomed.size > 1000 && pngWhole.size > 1000 && pngZoomed.sum !== pngWhole.sum, 'PNG export follows the zoom (zoomed and whole images differ)');
  // the zoom button works on an unzoomed plot, then Fit
  await page.click('#btn-plot-zoom-in'); await sleep(60);
  z = await zoom();
  check(z && near(z.x1 - z.x0, 1850, 1e-6), 'Zoom in on the whole plot: half of the range: ' + JSON.stringify(z));
  await page.click('#btn-plot-fit'); await sleep(60);
  check(!(await zoom()), 'Fit shows the whole plot');
  // a flat box zooms x only; double-click resets
  await dragF(0.2, 0.5, 0.6, 0.5);
  z = await zoom(); g = await geom();
  const full = await page.evaluate(() => { const p = BM.app.plot(); return p._ranges(p._full().L); });
  check(z && near(z.y0, full.y0, 1e-9) && near(z.y1, full.y1, 1e-9) && z.x1 - z.x0 < 0.5 * (full.x1 - full.x0), 'a flat box zooms the x axis only');
  await page.mouse.dblclick(g.rect.x + g.m.l + g.pw / 2, g.rect.y + g.m.t + g.ph / 2); await sleep(80);
  check(!(await zoom()), 'a double-click shows the whole plot again');
  // another y quantity starts again
  await dragF(0.2, 0.2, 0.6, 0.8);
  check(!!(await zoom()), 'zoomed again');
  await page.selectOption('#plot-controls select[data-pc="yMode"]', 'epsilon'); await sleep(80);
  check(!(await zoom()), 'another y quantity (ε) starts unzoomed');
  // a zoom is kept per plot: Raman then IR again
  await dragF(0.1, 0.1, 0.5, 0.9);
  const zIr = await zoom();
  await pickPlot('raman');
  check(!(await zoom()), 'another plot starts unzoomed');
  await pickPlot('ir');
  z = await zoom();
  check(z && zIr && near(z.x0, zIr.x0, 1e-9) && near(z.x1, zIr.x1, 1e-9), 'coming back to a plot shows its zoom again');
  await page.click('#btn-plot-fit');

  // tick labels of many bars (charges of a 3,000-atom system) never overlap, whole or zoomed
  const labelGaps = await page.evaluate(() => {
    const p = BM.app.plot(), keep = p.spec, n = 3000, c = p.canvas;
    p.set({ xRange: [0.4, n + 0.6], bars: Array.from({ length: n }, (_, i) => ({ x: i + 1, y: Math.sin(i), id: i })), barWidth: 0.72, xTickLabels: Array.from({ length: n }, (_, i) => ({ x: i + 1, label: 'H' + (i + 1) })), zoomKey: 'labels-test' });
    const gaps = () => {
      const t = []; const P = { line() {}, rect() {}, textV() {}, text(s, x, y) { if (/^H\d+$/.test(s)) t.push({ s, x }); } };
      p._axes(P, p._layout(c.clientWidth, c.clientHeight), p._theme(), c.clientWidth, c.clientHeight);
      let worst = Infinity;
      for (let k = 1; k < t.length; k++) worst = Math.min(worst, t[k].x - t[k - 1].x - (t[k].s.length + t[k - 1].s.length) / 2 * 10.5 * 0.56);
      return { n: t.length, worst };
    };
    const whole = gaps();
    p.zoomIn(0.05);
    const zoomed = gaps();
    p.fit(); p.set(keep);
    return { whole, zoomed };
  });
  check(labelGaps.whole.n > 5 && labelGaps.whole.worst >= 0, '3,000 bar labels do not overlap: ' + JSON.stringify(labelGaps.whole));
  check(labelGaps.zoomed.n > labelGaps.whole.n / 3 && labelGaps.zoomed.worst >= 0, 'nor when zoomed in: ' + JSON.stringify(labelGaps.zoomed));

  console.log('#9 zoomable plots: every plot type');
  /** Zooms into the middle of the plot shown, clicks a visible item, checks the selection, the CSV and the reset. */
  const everyPlot = async (label, t, box, selected) => {
    await pickPlot(t);
    const before = await geom();
    await dragF(box[0], box[1], box[2], box[3]);
    const zz = await zoom(), gg = await geom();
    check(!!zz, `${label}: a box zooms (${t})`);
    if (!zz) return;
    if (gg.type === 'levels') check(zz.x0 === 0 && zz.x1 === 1 && zz.y1 - zz.y0 < 0.7 * (before.L.y1 - before.L.y0), `${label}: the levels zoom in energy only: ` + JSON.stringify(zz));
    if (gg.logY) check(zz.y0 > 0 && zz.y1 > zz.y0 && Math.log10(zz.y1 / zz.y0) < 0.7 * Math.log10(before.L.y1 / before.L.y0), `${label}: a log axis zooms in log space: ` + JSON.stringify(zz));
    check(gg.items.length > 0, `${label}: something is shown in the zoom`);
    const target = gg.items.find(q => !q.sel) || gg.items[0];
    if (target) {
      await page.mouse.click(target.cx, target.cy); await sleep(120);
      check(await page.evaluate(selected, target.id), `${label}: a click in the zoomed plot selects item ${target.id}`);
      check(!!(await zoom()), `${label}: the zoom is kept after the selection`);
    }
    const c = await save('plot-csv');
    check(/^# zoomed in: only /m.test(c.text), `${label}: the CSV says it is zoomed (${(c.text.split('\n')[1] || '').slice(0, 60)})`);
    const g2 = await geom();
    await page.mouse.dblclick(g2.rect.x + g2.m.l + 4, g2.rect.y + g2.m.t + 4); await sleep(80);
    check(!(await zoom()), `${label}: double-click resets`);
    return c.text;
  };
  const frameSel = id => BM.app.state.active.frameIndex === id;
  let text = await everyPlot('energy (opt)', 'opt:energy', [0.3, 0.02, 0.98, 0.98], frameSel);
  check(text && text.split('\n').filter(l => /^\d/.test(l)).length < 4, 'opt energy CSV: only the steps in the zoom');
  await everyPlot('convergence (log)', 'opt:all', [0.3, 0.1, 0.98, 0.7], frameSel);
  await everyPlot('RMSD', 'opt:rmsd', [0.3, 0.02, 0.98, 0.98], frameSel);
  await everyPlot('orbital levels', 'mo', [0.1, 0.3, 0.9, 0.7], id => BM.app.state.active.moSel === id);
  text = await everyPlot('charges', 'charges', [0.4, 0.02, 0.98, 0.98], id => BM.app.state.active.selection.indexOf(id) >= 0);
  check(text && text.split('\n').filter(l => /^\d/.test(l)).length === 2, 'charges CSV: only the atoms in the zoom (' + (text || '').split('\n').filter(l => /^\d/.test(l)).map(l => l.split(',')[0]).join() + ')');
  await open('g16_water_scan.log');
  await everyPlot('scan', 'scan', [0.3, 0.02, 0.98, 0.98], frameSel);
  await open('qc/g16_hcn_irc.log');
  await everyPlot('IRC', 'irc', [0.2, 0.02, 0.98, 0.98], frameSel);
  await open('water_traj.xyz');
  await everyPlot('energy profile', 'energy', [0.3, 0.02, 0.98, 0.98], frameSel);
  await open('qc/g16_water_ecd.log');
  const stateSel = id => { const d = BM.app.state.active; return d.stateIndex >= 0 && d.data.excitedStates[d.stateIndex].index === id; };
  await everyPlot('UV-Vis', 'uv', [0.02, 0.02, 0.7, 0.98], stateSel);
  await everyPlot('ECD', 'ecd', [0.02, 0.02, 0.7, 0.98], stateSel);
  await open('qc/g16_water_nmr.log');
  await everyPlot('NMR', 'nmr', [0.02, 0.02, 0.6, 0.98], id => BM.app.state.active.selection.indexOf(id) >= 0);

  // ================================================================ #10 component energies
  console.log('#10 QM/MM and ONIOM component energies');
  for (const [file, names, steps] of [
    ['qmmm/g16_oniom_amber_opt.log', ['Total', 'High level, model system', 'Low level, model system', 'Low level, real system'], 6],
    ['qmmm/g16_oniom3_opt.log', ['Total', 'High level, model system', 'Medium level, model system', 'Medium level, intermediate system', 'Low level, intermediate system', 'Low level, real system'], 5],
    ['qmmm/orca6_qmmm_opt.out', ['Total', 'QM', 'MM', 'QM–MM interaction'], 5]]) {
    await open(file);
    const opts = await page.$eval('#plot-pick', sl => [...sl.options].map(o => [o.value, o.textContent]));
    const parts = opts.find(o => o[0] === 'opt:parts');
    check(!!parts, `${file}: an "Energy components" plot is offered (${parts && parts[1]})`);
    if (!parts) continue;
    await pickPlot('opt:parts');
    const sp = await doc(() => { const s = BM.app.plot().spec; return { names: s.lines.map(l => l.name), n: s.lines.map(l => l.pts.length), last: s.lines.map(l => l.pts[l.pts.length - 1].y), yLabel: s.yLabel }; });
    check(JSON.stringify(sp.names) === JSON.stringify(names), `${file}: the total and each component: ${sp.names.join(' / ')}`);
    check(sp.n.every(k => k === steps), `${file}: every step of each (${sp.n.join()})`);
    check(sp.last.every(v => Math.abs(v) < 1e-9) && /vs final/.test(sp.yLabel), `${file}: each relative to its own final value`);
    const c = (await save('plot-csv')).text.split('\n');
    const head = c.find(l => /^step,frame,/.test(l)) || '';
    check(head.split(',').length === 2 + names.length && /total_hartree/.test(head), `${file}: CSV columns: ${head}`);
    const rows = c.filter(l => /^\d/.test(l));
    const want = await doc(() => { const d = BM.app.state.active, f = d.data.frames[d.data.frames.length - 1]; return f.energy; });
    const lastRow = rows[rows.length - 1].split(',');
    check(rows.length === steps && Math.abs(+lastRow[2] - want) < 1e-9, `${file}: CSV rows hold the energies in hartree (${lastRow[2]} vs ${want})`);
    await page.click('#tabs [data-tab="info"]'); await sleep(80);
    const info = await page.$eval('#panel', e => e.textContent);
    check(names.slice(1).every(n => info.indexOf(n) >= 0), `${file}: the Info tab lists the components`);
  }
  // a plain optimisation has no component plot
  await open('g16_water_optfreq.log');
  check(!(await page.$eval('#plot-pick', sl => [...sl.options].some(o => o.value === 'opt:parts'))), 'no component plot without components');

  // ================================================================ #7 the loading dialog
  console.log('#7 choose what to load');
  const dialogOpen = () => page.evaluate(() => !!document.querySelector('#dlg[open]'));
  const openAsk = async file => { await page.setInputFiles('#file-input', path.join(FIX, file)); await page.waitForSelector('#dlg[open]', { timeout: 10000 }); await sleep(60); };
  const choose = async (steps, modes) => {
    await page.check(`#dlg input[name="ld-steps"][value="${steps[0]}"]`);
    if (steps[1] !== undefined) { const k = steps[0] === 'range' ? 'from' : steps[0]; await page.fill(`#dlg [data-ld="${k}"]`, String(steps[1])); }
    if (steps[2] !== undefined) await page.fill('#dlg [data-ld="to"]', String(steps[2]));
    if (modes) {
      await page.check(`#dlg input[name="ld-modes"][value="${modes[0]}"]`);
      if (modes[1] !== undefined) await page.fill('#dlg [data-ld="first"]', String(modes[1]));
    }
    await page.click('#dlg .dlg-foot .btn.primary');
    await page.waitForFunction(() => !document.querySelector('#dlg[open]'));
    await sleep(250);
  };
  const loadedInfo = async () => { await page.click('#tabs [data-tab="info"]'); await sleep(80); return page.$eval('#panel', e => { const s = e.querySelector('#sec-loaded'); return s ? s.innerText.replace(/\s+/g, ' ') : ''; }); };
  const reloadAsk = async () => { await page.click('#tabs [data-tab="info"]'); await sleep(60); await page.click('[data-act="reload-part"]'); await page.waitForSelector('#dlg[open]'); await sleep(40); };
  await page.evaluate(() => { window.__bmAskSize = 1; window.__bmBigSize = 1; });   // ask for every file, read the big way
  await openAsk('g16_water_optfreq.log');
  const dlgText = await page.$eval('#dlg', e => e.textContent.replace(/\s+/g, ' '));
  check(/g16_water_optfreq\.log/.test(dlgText) && /3 atoms/.test(dlgText) && /opt, freq/.test(dlgText), 'the dialog names the file, its atoms and jobs: ' + dlgText.slice(0, 160));
  check(/every/.test(dlgText) && /only the last/.test(dlgText) && /to/.test(dlgText) && /Vibrations/.test(dlgText) && /frequencies and intensities only/.test(dlgText), 'it offers every Nth, the last N, a range, and the vibrations');
  await choose(['last', 1], ['none']);
  let st = await doc(() => { const d = BM.app.state.active, L = d.data.meta.loaded; return { name: d.name, frames: d.data.frames.length, loaded: L, vib: d.data.vibrations && d.data.vibrations.freqs.length, has: d.data.vibrations && d.data.vibrations.modesLoaded }; });
  check(st.name === 'g16_water_optfreq.log' && st.loaded && st.loaded.steps.mode === 'last' && st.loaded.modes.mode === 'none', 'opened with the choice: ' + JSON.stringify(st.loaded && { steps: st.loaded.steps, modes: st.loaded.modes }));
  check(st.vib === 3 && st.has === 0, 'frequencies without displacements: ' + st.vib + ' modes, ' + st.has + ' with displacements');
  let li = await loadedInfo();
  check(/Loaded part of the file/i.test(li) && /Displacements 0 of 3 vibrations/.test(li) && /Load other parts/.test(li), 'the Info tab says what was loaded: ' + li);
  check(/Geometries \d of 4 \(the last 1\); the plots show the energies of every step/.test(li), 'the Info tab: geometries and how: ' + li);
  await page.click('#tabs [data-tab="vib"]'); await sleep(80);
  check(await page.$$eval('#vib-table tr.row', r => r.length === 3 && r.every(x => x.classList.contains('nodisp'))), 'the Vibrations tab marks modes without displacements');
  // Load other parts: the dialog again with the previous choice, then the whole file
  await reloadAsk();
  check(await page.$eval('#dlg input[name="ld-steps"][value="last"]', r => r.checked) && await page.$eval('#dlg input[name="ld-modes"][value="none"]', r => r.checked), 'Load other parts starts from the previous choice');
  await choose(['all'], ['all']);
  st = await doc(() => { const d = BM.app.state.active; return { n: BM.app.state.docs.length, frames: d.data.frames.length, has: d.data.vibrations.modesLoaded, mode: d.data.meta.loaded.steps.mode }; });
  check(st.frames === 4 && st.has === 3 && st.mode === 'all', 'reloaded with everything: ' + JSON.stringify(st));
  li = await loadedInfo();
  check(/Geometries 4 of 4/.test(li) && /all 3 vibrations/.test(li), 'the Info tab follows: ' + li);
  // every Nth and a range on a scan (the first and last geometry of each job are kept)
  await openAsk('g16_water_scan.log');
  await choose(['every', 2]);
  st = await doc(() => { const d = BM.app.state.active, L = d.data.meta.loaded.steps; return { frames: d.data.frames.length, total: L.total, loaded: L.loaded, mode: L.mode, n: L.n, steps: d.data.steps ? d.data.steps.length : 0 }; });
  check(st.mode === 'every' && st.n === 2 && st.loaded === st.frames && st.total >= st.frames, 'every 2nd step: ' + JSON.stringify(st));
  await reloadAsk();
  await choose(['range', 2, 3]);
  st = await doc(() => { const d = BM.app.state.active, L = d.data.meta.loaded.steps; return { frames: d.data.frames.length, mode: L.mode, from: L.from, to: L.to }; });
  check(st.mode === 'range' && st.from === 2 && st.to === 3, 'a range of steps: ' + JSON.stringify(st));
  // a multi-frame XYZ: every Nth frame
  await openAsk('water_traj.xyz');
  check(/Frames/.test(await page.$eval('#dlg', e => e.textContent)), 'an XYZ trajectory is asked about in frames');
  await choose(['every', 2]);
  st = await doc(() => { const d = BM.app.state.active, L = d.data.meta.loaded; return { frames: d.data.frames.length, L: L && L.steps }; });
  check(st.frames === 2 && st.L && st.L.mode === 'every' && st.L.total === 4, 'XYZ every 2nd frame: ' + JSON.stringify(st));
  await reloadAsk();
  await choose(['last', 3]);
  check((await doc(() => BM.app.state.active.data.frames.length)) === 3, 'XYZ reloaded as the last 3 frames');
  // a DCD after a matching structure: a range of frames
  const dcd = (() => {   // CHARMM DCD, 3 atoms, 9 frames (independent of the reader)
    const n = 3, frames = Array.from({ length: 9 }, (_, k) => [0, 0, 0.12 + 0.01 * k, 0, 0.76, -0.47, 0, -0.76, -0.47]);
    const i32 = v => { const b = Buffer.alloc(4); b.writeInt32LE(v); return b; }, rec = body => Buffer.concat([i32(body.length), body, i32(body.length)]);
    const f32 = a => { const b = Buffer.alloc(4 * a.length); a.forEach((v, k) => b.writeFloatLE(v, 4 * k)); return b; };
    const hdr = Buffer.alloc(84); hdr.write('CORD', 0, 'latin1');
    const ic = new Array(20).fill(0); ic[0] = frames.length; ic[2] = 1; ic[19] = 24;
    ic.forEach((v, k) => { if (k !== 9) hdr.writeInt32LE(v, 4 + 4 * k); });
    hdr.writeFloatLE(0.002, 40);
    const t = Buffer.alloc(84, 32); t.writeInt32LE(1, 0); t.write('dialog test', 4, 'latin1');
    const parts = [rec(hdr), rec(t), rec(i32(n))];
    for (const fr of frames) for (let d = 0; d < 3; d++) parts.push(rec(f32([0, 1, 2].map(k => fr[3 * k + d]))));
    return Buffer.concat(parts);
  })();
  await page.evaluate(() => { window.__bmAskSize = 0; window.__bmBigSize = 0; });
  await open('g16_water_tddft.log');   // 3 atoms, like the DCD
  await page.evaluate(() => { window.__bmAskSize = 1; window.__bmBigSize = 1; });
  await page.setInputFiles('#file-input', { name: 'water.dcd', mimeType: 'application/octet-stream', buffer: dcd });
  await page.waitForSelector('#dlg[open]');
  const ddlg = await page.$eval('#dlg', e => e.textContent.replace(/\s+/g, ' '));
  check(/3 atoms/.test(ddlg) && /9 frames/.test(ddlg), 'the DCD dialog: ' + ddlg.slice(0, 120));
  await choose(['range', 3, 5]);
  st = await doc(() => { const d = BM.app.state.active; return { n: d.data.frames.length, z: Array.from(d.data.frames.slice(-3).map(f => +f.xyz[2].toFixed(3))), t: d.trajLoad && d.trajLoad.count }; });
  check(st.t === 3 && JSON.stringify(st.z) === '[0.14,0.15,0.16]', 'DCD frames 3–5 added to the matching structure: ' + JSON.stringify(st));
  // cancel opens nothing
  const nDocs = await doc(() => BM.app.state.docs.length);
  await openAsk('g16_water_tddft.log');
  await page.click('#dlg .dlg-foot .btn:not(.primary)');
  await sleep(200);
  check(!(await dialogOpen()) && (await doc(() => BM.app.state.docs.length)) === nDocs, 'Cancel opens nothing');

  // ================================================================ #8 + #7 Amber NetCDF after a prmtop
  console.log('#8 Amber NetCDF trajectories and restarts');
  await page.evaluate(() => { window.__bmAskSize = 0; window.__bmBigSize = 0; });
  await open('struct/water.prmtop');
  check(await doc(() => BM.app.state.active.pending), 'the topology waits for coordinates');
  await page.setInputFiles('#file-input', path.join(FIX, 'netcdf', 'traj_64bit.nc'));
  await page.waitForFunction(() => !BM.app.state.active.pending);
  await sleep(150);
  const exp = JSON.parse(fs.readFileSync(path.join(FIX, 'netcdf', 'expected.json'), 'utf8'));
  st = await doc(() => { const d = BM.app.state.active, f = d.data.frames; return { n: f.length, x0: Array.from(f[0].xyz), x11: Array.from(f[11].xyz), Z: f[0].Z, cells: f.filter(x => x.cell).length, label: f[2].label, time: f[2].time }; });
  check(st.n === 12 && JSON.stringify(st.Z) === '[8,1,1,11]', 'NetCDF frames on the prmtop atoms: ' + st.n + ' frames, Z ' + JSON.stringify(st.Z));
  check(JSON.stringify(st.x0) === JSON.stringify(exp['traj_64bit.nc'].xyz[0]) && JSON.stringify(st.x11) === JSON.stringify(exp['traj_64bit.nc'].xyz[11]), 'coordinates as written by SciPy');
  check(st.cells === 12 && st.time === 6 && /Frame 3/.test(st.label), 'box per frame and time: ' + JSON.stringify({ cells: st.cells, time: st.time, label: st.label }));
  await page.setInputFiles('#file-input', path.join(FIX, 'netcdf', 'restart_netcdf.rst7'));
  await page.waitForFunction(() => BM.app.state.active.data.frames.length === 13);
  st = await doc(() => { const d = BM.app.state.active, f = d.data.frames[12]; return { label: f.label, x: Array.from(f.xyz) }; });
  check(st.label === 'Restart' && JSON.stringify(st.x) === JSON.stringify(exp['restart_netcdf.rst7'].xyz[0]), 'a NetCDF restart named .rst7 is added: ' + st.label);
  // a trajectory big enough for the dialog: every 3rd frame, then the last 5
  await page.evaluate(() => { window.__bmAskSize = 1; window.__bmBigSize = 1; });
  await open('struct/water.prmtop');
  await openAsk('netcdf/traj_64bit.nc');
  const tdlg = await page.$eval('#dlg', e => e.textContent.replace(/\s+/g, ' '));
  check(/Amber NetCDF/.test(tdlg) && /4 atoms/.test(tdlg) && /12 frames/.test(tdlg) && !/about 12 frames/.test(tdlg), 'the trajectory dialog: format, atoms, exact frame count: ' + tdlg.slice(0, 140));
  await choose(['every', 3]);
  st = await doc(() => { const d = BM.app.state.active; return { n: d.data.frames.length, labels: d.data.frames.map(f => f.label.split(' ·')[0]), t: d.trajLoad && { total: d.trajLoad.total, stride: d.trajLoad.stride, count: d.trajLoad.count } }; });
  check(st.n === 4 && st.labels.join() === 'Frame 1,Frame 4,Frame 7,Frame 10', 'every 3rd frame: ' + st.labels.join());
  check(st.t && st.t.total === 12 && st.t.stride === 3 && st.t.count === 4, 'the document remembers the part loaded: ' + JSON.stringify(st.t));
  li = await loadedInfo();
  check(/Frames of traj_64bit\.nc 4 of 12 \(every 3rd\)/.test(li) && /Load other parts/.test(li), 'Info tab: ' + li);
  await reloadAsk();
  check(await page.$eval('#dlg input[name="ld-steps"][value="every"]', r => r.checked) && (await page.$eval('#dlg [data-ld="every"]', e => e.value)) === '3', 'reload starts from every 3rd');
  await choose(['last', 5]);
  st = await doc(() => { const d = BM.app.state.active; return { n: d.data.frames.length, labels: d.data.frames.map(f => f.label.split(' ·')[0]) }; });
  check(st.n === 5 && st.labels.join() === 'Frame 8,Frame 9,Frame 10,Frame 11,Frame 12', 'reloaded as the last 5: ' + st.labels.join());
  li = await loadedInfo();
  check(/Frames of traj_64bit\.nc 5 of 12 \(8–12\)/.test(li), 'Info tab after the reload: ' + li);
  await reloadAsk();
  await choose(['range', 2, 4]);
  st = await doc(() => BM.app.state.active.data.frames.map(f => f.label.split(' ·')[0]).join());
  check(st === 'Frame 2,Frame 3,Frame 4', 'a frame range: ' + st);
  await page.evaluate(() => { window.__bmAskSize = 0; window.__bmBigSize = 0; });

  // ================================================================ the vibrations of a big frequency job
  console.log('Vibrations: every block, also after an IR intensity printed as asterisks');
  const FS = require('./test_freq_sections.js');
  for (const big of [false, true]) {
    await page.evaluate(b => { window.__bmAskSize = 0; window.__bmBigSize = b ? 1 : 0; }, big);
    await page.setInputFiles('#file-input', { name: 'stars-freq.log', mimeType: 'text/plain', buffer: Buffer.from(FS.log({ stars: true })) });
    await page.waitForFunction(() => BM.app.state.active && BM.app.state.active.name === 'stars-freq.log' && BM.app.state.active.data.vibrations);
    await sleep(150);
    const how = big ? 'condensed' : 'whole';
    st = await doc(() => { const v = BM.app.state.active.data.vibrations; return { n: v.freqs.length, ir: !!v.ir }; });
    check(st.n === FS.NM && st.ir, `${how}: all ${FS.NM} vibrations and their IR intensities (${st.n})`);
    await page.click('#tabs [data-tab="vib"]'); await sleep(80);
    check((await page.$$('#vib-table tr.row')).length === FS.NM && (await page.$eval('#tabs [data-tab="vib"]', e => e.textContent)).includes(String(FS.NM)), `${how}: the Vibrations tab lists them all`);
    await page.click('#vib-table tr.row[data-i="7"]'); await sleep(400);
    st = await doc(() => { const d = BM.app.state.active, eq = d.data.frames[d.data.vibrations.frameIndex].xyz, x = d._animXyz; let m = 0; for (let i = 0; i < eq.length; i++) m = Math.max(m, Math.abs(x[i] - eq[i])); return { vib: d.vibIndex, moved: m, anim: !!d._vibVec }; });
    check(st.vib === 7 && st.anim && st.moved > 0.05, `${how}: a click on mode 8 animates it (largest shift so far ${st.moved.toFixed(3)} Å)`);
  }
  await page.evaluate(() => { window.__bmAskSize = 0; window.__bmBigSize = 0; });
  // a vibration without displacements in the file says so instead of standing still silently
  await page.setInputFiles('#file-input', { name: 'nodisp-freq.log', mimeType: 'text/plain', buffer: Buffer.from(FS.log({ damaged: true })) });
  await page.waitForFunction(() => BM.app.state.active && BM.app.state.active.name === 'nodisp-freq.log');
  await sleep(150);
  st = await doc(() => { const d = BM.app.state.active; return { n: d.data.vibrations.freqs.length, notes: d.data.warnings.join(' | ') }; });
  check(st.n === 9 && /stopped early, after 9 vibrations/.test(st.notes), 'a damaged row: the blocks before it, and a reader note: ' + st.notes.slice(0, 120));
  await page.setInputFiles('#file-input', { name: 'norows-freq.log', mimeType: 'text/plain', buffer: Buffer.from(FS.log({ noRowsLast: true })) });
  await page.waitForFunction(() => BM.app.state.active && BM.app.state.active.name === 'norows-freq.log');
  await sleep(150);
  await page.click('#tabs [data-tab="vib"]'); await sleep(80);
  await page.click(`#vib-table tr.row[data-i="${FS.NM - 1}"]`); await sleep(200);
  const stText = await page.$eval('#status-msg', e => e.textContent);
  check(/no atom displacements for this vibration/.test(stText) && (await doc(() => BM.app.state.active._vibVec)) === null, 'a vibration printed without displacements says so when clicked: ' + stText);

  // ================================================================ summary
  const errs = errors.filter(e => !/favicon|WebGL|GPU stall|swiftshader/i.test(e));
  check(errs.length === 0, 'no page errors: ' + errs.join(' | '));
  await browser.close();
  console.log(`Milestone 1.5 browser tests: ${passes} passed, ${failures} failed`);
  process.exit(failures ? 1 : 0);
})().catch(e => { console.error(e); process.exit(1); });
