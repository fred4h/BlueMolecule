// End-to-end checks of the building and editing tools in dist/BlueMolecule.html (headless Chromium).
// Real mouse events are sent at canvas positions computed with R().worldToScreen.
// Usage: node tests/browser_build_test.js [screenshotDir]
const { chromium } = require('playwright');
const path = require('path');
const fs = require('fs');

const ROOT = path.resolve(__dirname, '..');
const HTML = 'file://' + path.join(ROOT, 'dist', 'BlueMolecule.html');
const SHOTS = process.argv[2] || path.join(ROOT, 'tests', 'shots', 'build');
fs.mkdirSync(SHOTS, { recursive: true });

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
    window.__bmSaveHook = async (name, blob) => { window.__saved.push({ name, text: blob.type === 'image/png' ? null : await blob.text(), size: blob.size }); };
  });
  const shot = name => page.screenshot({ path: path.join(SHOTS, name + '.png') });
  /** Expand a Build-tab section (sections other than Insert and Edit start collapsed). */
  const openSec = id => page.evaluate(id => { const s = document.querySelector('.bsec[data-sec="' + id + '"]'); if (s && s.classList.contains('folded')) s.querySelector('h3').click(); }, id);
  const box = await page.locator('#gl').boundingBox();
  const info = () => page.evaluate(() => {
    const a = BM.app, d = a.state.active, f = d.data.frames[d.frameIndex];
    return { name: d.name, n: f.Z.length, formula: BM.geom.formula(f.Z), bonds: f.bonds ? f.bonds.length : -1, sel: d.selection.slice(), tool: a.state.build.tool, frames: d.data.frames.length, fi: d.frameIndex };
  });
  /** Page coordinates of atom i (or a world point). */
  const atomXY = i => page.evaluate(i => { const R = BM.app.R(), x = R.xyz; const s = R.worldToScreen([x[3 * i], x[3 * i + 1], x[3 * i + 2]]); return [s.x, s.y]; }, i).then(([x, y]) => [box.x + x, box.y + y]);
  const pick = (x, y) => page.evaluate(([x, y]) => BM.app.R().pickAny(x, y), [x - box.x, y - box.y]);
  const drag = async (x0, y0, x1, y1, opts) => {
    opts = opts || {};
    if (opts.shift) await page.keyboard.down('Shift');
    await page.mouse.move(x0, y0); await page.mouse.down();
    await page.mouse.move(x1, y1, { steps: opts.steps || 8 });
    await page.mouse.up();
    if (opts.shift) await page.keyboard.up('Shift');
    await sleep(opts.wait || 120);
  };
  const waitIdle = () => page.waitForFunction(() => !BM.app.ui().running(), null, { timeout: 20000 });
  const heavy = () => page.evaluate(() => { const d = BM.app.state.active, f = d.data.frames[d.frameIndex]; return f.Z.map((z, i) => (z > 1 ? i : -1)).filter(i => i >= 0); });
  /** A C–C (or given) bond whose midpoint picks as that bond, as {k, i, j, x, y} in page coordinates. */
  const clickableBond = (filter) => page.evaluate(filter => {
    const R = BM.app.R(), b = R.bonds, x = R.xyz, Z = R.Z, out = [];
    for (let k = 0; 2 * k + 1 < b.length; k++) {
      const i = b[2 * k], j = b[2 * k + 1];
      if (filter === 'heavy' && (Z[i] === 1 || Z[j] === 1)) continue;
      const s = R.worldToScreen([(x[3 * i] + x[3 * j]) / 2, (x[3 * i + 1] + x[3 * j + 1]) / 2, (x[3 * i + 2] + x[3 * j + 2]) / 2]);
      const h = R.pickAny(s.x, s.y);
      if (h.kind === 'bond' && h.index === k) out.push({ k, i, j, x: s.x, y: s.y });
    }
    return out;
  }, filter).then(list => list.map(o => Object.assign(o, { x: o.x + box.x, y: o.y + box.y })));
  /** A point near (x, y) where nothing is picked. */
  const emptyNear = async (x, y) => {
    for (const r of [0, 20, 40, 60, 90]) for (let a = 0; a < 8; a++) {
      const px = x + r * Math.cos(a * Math.PI / 4), py = y + r * Math.sin(a * Math.PI / 4);
      if ((await pick(px, py)).kind === 'none') return [px, py];
    }
    return [x, y];
  };

  // ---------------------------------------------------------------- 1. draw tool
  console.log('Draw tool');
  await shot('00_start_screen');
  await page.click('#btn-new-2');
  await page.waitForFunction(() => BM.app.state.active && /^New molecule/.test(BM.app.state.active.name));
  await page.keyboard.press('d');
  let s = await info();
  check(s.tool === 'draw' && s.n === 0, 'new molecule opens empty in the Draw tool ' + JSON.stringify(s));
  await page.mouse.click(box.x + box.width / 2, box.y + box.height / 2);
  await sleep(150);
  s = await info();
  check(s.n === 5 && s.formula === 'CH4', 'click on empty space places C, Adjust H gives methane: ' + s.formula);
  const origin = await page.evaluate(() => { const f = BM.app.curFrame(BM.app.state.active); return [f.xyz[0], f.xyz[1], f.xyz[2]]; });
  check(origin.every(v => Math.abs(v) < 1e-9), 'first atom at the origin ' + origin);
  const dirs = [[150, -85], [150, 85], [150, -85]];
  for (const [dx, dy] of dirs) {
    const hv = await heavy();
    const xs = await Promise.all(hv.map(atomXY));
    let best = 0; xs.forEach((p, k) => { if (p[0] > xs[best][0]) best = k; });
    const [x0, y0] = xs[best];
    await drag(x0, y0, x0 + dx, y0 + dy);
  }
  s = await info();
  check(s.formula === 'C4H10' && s.n === 14, 'three drags give butane C4H10: ' + s.formula + ' / ' + s.n);
  const cc = await page.evaluate(() => { const d = BM.app.state.active, f = BM.app.curFrame(d); return f.bonds.filter(b => f.Z[b[0]] === 6 && f.Z[b[1]] === 6).map(b => BM.geom.dist(f.xyz, b[0], b[1])); });
  check(cc.length === 3 && cc.every(r => Math.abs(r - 1.52) < 0.03), 'C–C bonds at the ideal length: ' + cc.map(r => r.toFixed(3)));
  let bonds = await clickableBond('heavy');
  check(bonds.length > 0, 'a C–C bond is clickable');
  if (bonds.length) {
    const b0 = bonds[0];
    await page.mouse.click(b0.x, b0.y);
    await sleep(150);
    s = await info();
    const ord = await page.evaluate(([i, j]) => { const f = BM.app.curFrame(BM.app.state.active); const b = f.bonds.find(q => (q[0] === i && q[1] === j) || (q[0] === j && q[1] === i)); return b ? b[2] : null; }, [b0.i, b0.j]);
    check(ord === 2, 'click on a bond sets order 2: ' + ord);
    check(s.formula === 'C4H8' && s.n === 12, 'Adjust H after the double bond: C4H8, got ' + s.formula);
    const multi = await page.evaluate(() => BM.app.R().nMultiCyl);
    check(multi >= 4, 'double bond drawn as two thin cylinders (' + multi + ')');
  }
  await page.mouse.move(box.x + 300, box.y + 300);
  await shot('01_draw_butene');
  // element change: a terminal, singly bonded carbon becomes O
  const term = await page.evaluate(() => {
    const f = BM.app.curFrame(BM.app.state.active), Z = f.Z;
    const heavyNb = i => f.bonds.filter(b => (b[0] === i || b[1] === i) && Z[b[0] + b[1] - i] > 1);
    return Z.map((z, i) => i).find(i => Z[i] === 6 && heavyNb(i).length === 1 && heavyNb(i)[0][2] === 1);
  });
  await page.click('.toolopts [data-el="8"]');
  const [ox, oy] = await atomXY(term);
  if ((await pick(ox, oy)).kind === 'atom') {
    await page.mouse.click(ox, oy); await sleep(120);
    s = await info();
    check(s.formula === 'C3H6O', 'click on an atom with O selected: ' + s.formula);
    await page.keyboard.press('Control+z'); await sleep(120);
    s = await info();
    check(s.formula === 'C4H8', 'undo element change: ' + s.formula);
  } else check(false, 'terminal carbon not pickable');
  await page.click('.toolopts [data-el="6"]');
  const butene = await page.evaluate(() => BM.app.state.active.id);

  // ---------------------------------------------------------------- 2. SMILES
  console.log('SMILES and hydrogens');
  await page.click('#btn-new');
  await page.waitForFunction(() => BM.app.state.active.data.frames[0].Z.length === 0);
  await page.fill('#bp-smiles', 'CC(=O)Oc1ccccc1C(=O)O');
  await page.click('[data-act="b-smiles"]');
  await sleep(200);
  s = await info();
  check(s.n === 21 && s.formula === 'C9H8O4', 'aspirin from SMILES: ' + s.n + ' atoms, ' + s.formula);
  const arom = await page.evaluate(() => BM.app.curFrame(BM.app.state.active).bonds.filter(b => b[2] === 1.5).length);
  check(arom === 6, 'aromatic ring stored with order 1.5 (' + arom + ')');
  await page.click('[data-act="b-copy-smiles"]');
  await sleep(100);
  const smi = await page.textContent('.smiles-out');
  const smiBack = await page.evaluate(t => { const m = BM.chem.smilesTo3D(t).mol; return BM.geom.formula(m.Z); }, smi);
  check(smiBack === 'C9H8O4', 'Copy SMILES round-trips: ' + smi + ' → ' + smiBack);
  await shot('02_smiles_aspirin');
  await page.click('[data-act="b-h-remove"]');
  await sleep(150);
  s = await info();
  check(s.formula === 'C9O4' && s.n === 13, 'remove all H: ' + s.formula);
  await page.click('[data-act="b-h-add"]');
  await sleep(150);
  s = await info();
  check(s.formula === 'C9H8O4' && s.n === 21, 'add H back gives the same formula: ' + s.formula);

  // ---------------------------------------------------------------- 3. templates
  console.log('Templates');
  await page.click('#btn-new');
  await page.waitForFunction(() => BM.app.state.active.data.frames[0].Z.length === 0);
  await page.keyboard.press('t');
  await page.click('.toolopts [data-mode="metal"]');
  check(await page.evaluate(() => { const t = BM.app.state.build.tpl; return t.metalZ === 78 && t.geometry === 'square-planar' && t.metalLigand === 'Cl'; }), 'metal centre defaults Pt / square planar / Cl');
  await page.mouse.click(box.x + box.width / 2, box.y + box.height / 2);
  await sleep(200);
  s = await info();
  check(s.formula === 'Cl4Pt' && s.n === 5, 'buildCenter PtCl4: ' + s.formula);
  await page.click('.toolopts [data-mode="ligands"]');
  await page.click('.toolopts [data-pick="ligand"]');
  await page.fill('.popup .pick-q', 'ammine');
  await page.waitForTimeout(50);
  await shot('04a_ligand_picker');
  await page.keyboard.press('Enter');
  check(await page.evaluate(() => BM.app.state.build.tpl.ligand) === 'NH3', 'ligand picker selects NH3');
  const hintL = await page.evaluate(() => BM.app.ui().hintText());
  check(hintL === await page.textContent('#status-msg'), 'changing the template option shows its hint at once');
  check(/charge 0/.test(hintL) && /suggested total charge −2/.test(hintL) && /Pt\(II\)/.test(hintL), 'ligand hint shows charges: ' + hintL);
  for (let rep = 0; rep < 2; rep++) {
    const target = await page.evaluate(() => {
      const f = BM.app.curFrame(BM.app.state.active), Z = f.Z, pt = Z.indexOf(78), n = Z.indexOf(7);
      const cls = Z.map((z, i) => (z === 17 ? i : -1)).filter(i => i >= 0);
      if (n < 0) return cls[0];
      return cls.find(i => Math.abs(BM.geom.angle(f.xyz, n, pt, i) - 90) < 10);
    });
    const [x, y] = await atomXY(target);
    const h = await pick(x, y);
    check(h.kind === 'atom' && h.index === target, 'Cl ' + target + ' is clickable');
    await page.mouse.click(x, y);
    await sleep(250);
  }
  s = await info();
  check(s.formula === 'Cl2H6N2Pt', 'two cis Cl replaced by NH3: ' + s.formula);
  const cisAng = await page.evaluate(() => { const f = BM.app.curFrame(BM.app.state.active), Z = f.Z, pt = Z.indexOf(78), ns = Z.map((z, i) => (z === 7 ? i : -1)).filter(i => i >= 0); return BM.geom.angle(f.xyz, ns[0], pt, ns[1]); });
  check(Math.abs(cisAng - 90) < 8, 'the two N are cis: ' + cisAng.toFixed(1) + '°');
  const hintC = await page.evaluate(() => BM.app.ui().hintText());
  check(/suggested total charge 0/.test(hintC), 'cisplatin: suggested total charge 0 (' + hintC + ')');
  await shot('04_template_cisplatin');
  // groups: methyl on an ammine H
  await page.click('.toolopts [data-mode="groups"]');
  const hN = await page.evaluate(() => { const f = BM.app.curFrame(BM.app.state.active); return f.Z.indexOf(1); });
  const [hx, hy] = await atomXY(hN);
  check((await pick(hx, hy)).kind === 'atom', 'ammine H pickable');
  await page.mouse.click(hx, hy); await sleep(200);
  s = await info();
  check(s.formula === 'CH8Cl2N2Pt', 'methyl attached to an H: ' + s.formula);
  await page.keyboard.press('Control+z'); await sleep(150);
  // molecules: click inserts benzene
  await page.click('.toolopts [data-mode="molecules"]');
  const [ex, ey] = await emptyNear(box.x + 120 + 150, box.y + box.height - 160);
  await page.mouse.click(ex, ey); await sleep(250);
  s = await info();
  check(s.formula === 'C6H12Cl2N2Pt', 'benzene inserted where clicked: ' + s.formula);
  await page.keyboard.press('Control+z'); await sleep(150);
  // bidentate ligand: pick two cis Cl sites on a fresh PtCl4, then "Apply to selected sites"
  await page.click('#btn-new');
  await page.keyboard.press('t');
  await page.click('.toolopts [data-mode="metal"]');
  await page.mouse.click(box.x + box.width / 2, box.y + box.height / 2); await sleep(200);
  await page.click('.toolopts [data-mode="ligands"]');
  await page.click('.toolopts [data-pick="ligand"]');
  await page.fill('.popup .pick-q', 'ethylenediamine'); await page.keyboard.press('Enter');
  const cisCl = await page.evaluate(() => { const f = BM.app.curFrame(BM.app.state.active), c = f.Z.map((z, i) => (z === 17 ? i : -1)).filter(i => i >= 0); return [c[0], c.find(i => i !== c[0] && Math.abs(BM.geom.angle(f.xyz, c[0], 0, i) - 90) < 10)]; });
  for (const i of cisCl) { const [a, b] = await atomXY(i); await page.mouse.click(a, b); await sleep(80); }
  check((await info()).sel.length === 2 && !(await page.$eval('.toolopts [data-tact="apply-sites"]', b => b.disabled)), 'two cis sites picked for the bidentate ligand');
  await page.click('.toolopts [data-tact="apply-sites"]'); await sleep(250);
  s = await info();
  check(s.formula === 'C2H8Cl2N2Pt', 'en placed on the two selected sites: ' + s.formula);
  await shot('04b_template_bidentate_en');

  // ---------------------------------------------------------------- 4. optimize
  console.log('Clean-up and symmetry');
  await page.evaluate(id => { const a = BM.app; a.activate(a.state.docs.find(d => d.id === id)); }, butene);
  await page.click('#tabs [data-tab="build"]');
  await page.click('[data-act="b-opt"]');
  await waitIdle();
  const ffr = await page.evaluate(() => BM.app.state.active._ff);
  check(ffr && ffr.e < ffr.e0 - 1e-6, 'Optimize geometry lowers the clean-up energy: ' + (ffr && ffr.e0.toFixed(3) + ' → ' + ffr.e.toFixed(3)));
  check(ffr && ffr.converged, 'clean-up converged');
  await shot('05_optimized');

  // ---------------------------------------------------------------- 5. benzene symmetry
  await page.click('#btn-new');
  await page.fill('#bp-smiles', 'c1ccccc1');
  await page.press('#bp-smiles', 'Enter');
  await sleep(200);
  await page.click('[data-act="b-opt"]');
  await waitIdle();
  await openSec('sym');
  await page.selectOption('#bp-symtol', '0.1');
  await page.click('[data-act="b-sym-apply"]');
  await sleep(150);
  let sym = await page.evaluate(() => BM.app.state.active._sym);
  if (!sym || sym.group !== 'D6h') {
    await page.selectOption('#bp-symtol', '0.3');
    await page.click('[data-act="b-sym-detect"]');
    await sleep(100);
    sym = await page.evaluate(() => BM.app.state.active._sym);
  }
  check(sym && sym.group === 'D6h', 'benzene symmetrized/detected as D6h: ' + JSON.stringify(sym));
  await page.click('[data-act="b-sym-detect"]');
  await sleep(100);
  const det = await page.evaluate(() => BM.app.state.active._sym);
  check(det && det.group === 'D6h' && det.order === 24, 'detect point group after symmetrize: ' + JSON.stringify(det));
  await shot('06_benzene_D6h');

  // ---------------------------------------------------------------- 6. bond-centric
  console.log('Bond-centric, manipulate, select');
  await page.click('#btn-new');
  await page.fill('#bp-smiles', 'CCCC');
  await page.click('[data-act="b-smiles"]');
  await sleep(200);
  await page.keyboard.press('k');
  const mid = await page.evaluate(() => {
    const f = BM.app.curFrame(BM.app.state.active), Z = f.Z, adj = Z.map(() => []);
    for (const b of f.bonds) { adj[b[0]].push(b[1]); adj[b[1]].push(b[0]); }
    const c = Z.map((z, i) => (z === 6 ? i : -1)).filter(i => i >= 0);
    const inner = c.filter(i => adj[i].filter(j => Z[j] === 6).length === 2);
    const ends = inner.map(i => adj[i].find(j => Z[j] === 6 && inner.indexOf(j) < 0));
    return { a: ends[0], b: inner[0], c: inner[1], d: ends[1] };
  });
  const dih = () => page.evaluate(m => { const f = BM.app.curFrame(BM.app.state.active); return BM.geom.dihedral(f.xyz, m.a, m.b, m.c, m.d); }, mid);
  const d0 = await dih();
  let cb = (await clickableBond('heavy')).find(o => (o.i === mid.b && o.j === mid.c) || (o.i === mid.c && o.j === mid.b));
  if (!cb) { await page.evaluate(() => BM.app.R().rotateBy([1, 0, 0], 0.5)); await sleep(80); cb = (await clickableBond('heavy')).find(o => (o.i === mid.b && o.j === mid.c) || (o.i === mid.c && o.j === mid.b)); }
  check(!!cb, 'central C–C bond is pickable');
  if (cb) {
    await page.mouse.move(cb.x, cb.y); await sleep(80);
    const hl = await page.evaluate(() => BM.app.R().highlightedBond);
    check(hl === cb.k, 'hover highlights the bond (' + hl + ')');
    const undoN = await page.evaluate(() => BM.app.state.active.undo.length);
    const [p1, p2] = [await atomXY(mid.b), await atomXY(mid.c)];
    const ul = Math.hypot(p2[0] - p1[0], p2[1] - p1[1]) || 1, perp = [-(p2[1] - p1[1]) / ul, (p2[0] - p1[0]) / ul];
    await page.mouse.down(); await page.mouse.move(cb.x + 80 * perp[0], cb.y + 80 * perp[1], { steps: 10 });
    const chipLive = await page.textContent('#measure');
    await page.mouse.up(); await sleep(120);
    const d1 = await dih();
    const ddiff = Math.abs(((d1 - d0 + 540) % 360) - 180);
    check(ddiff > 15, 'bond-centric drag changes the dihedral: ' + d0.toFixed(1) + ' → ' + d1.toFixed(1));
    check(/dihedral/.test(chipLive), 'chip shows the live dihedral: ' + chipLive);
    const undoN2 = await page.evaluate(() => BM.app.state.active.undo.length);
    check(undoN2 === undoN + 1, 'one undo step per drag (' + undoN + ' → ' + undoN2 + ')');
    const lens = await page.evaluate(() => { const f = BM.app.curFrame(BM.app.state.active); return f.bonds.map(b => BM.geom.dist(f.xyz, b[0], b[1])); });
    check(lens.every(r => r > 0.9 && r < 1.7), 'bond lengths unchanged by the rotation');
    await page.mouse.move(cb.x, cb.y); await sleep(50);
    await shot('07_bond_centric');
  }

  // ---------------------------------------------------------------- 7. coordinate editor + H-bonds + manipulate + undo
  await page.click('#btn-new');
  await page.click('.bsec[data-sec="coords"] h3');
  await page.waitForSelector('#bp-coords');
  const dimer = ['O  -1.5510  -0.1145   0.0000', 'H  -1.9343   0.7625   0.0000', 'H  -0.5997   0.0407   0.0000',
    'O   1.3505   0.1112   0.0000', 'H   1.6803  -0.3740  -0.7563', 'H   1.6803  -0.3740   0.7563'].join('\n');
  await page.fill('#bp-coords', dimer);
  await page.click('[data-act="b-coords-apply"]');
  await sleep(200);
  s = await info();
  check(s.n === 6 && s.formula === 'H4O2' && s.bonds === 4, 'coordinate editor Apply: ' + JSON.stringify(s));
  const back = await page.inputValue('#bp-coords');
  check(back.split('\n').length === 6 && /^O\s+-1\.55100000\s+-0\.11450000\s+0\.00000000$/.test(back.split('\n')[0]), 'editor shows Symbol x y z with 8 decimals: ' + back.split('\n')[0]);
  await page.click('#dd-settings [data-dd]');
  await page.check('#set-hbonds');
  await sleep(100);
  const hb = await page.evaluate(() => BM.app.R().dashedBonds.map(p => p.join('-')));
  check(hb.length === 1 && hb[0] === '2-3', 'water dimer: one H-bond O–H…O drawn (' + hb + ')');
  await page.click('#status');
  await shot('08_water_dimer_hbond');
  // manipulate: drag the second water away
  await page.keyboard.press('m');
  const before = await page.evaluate(() => Array.from(BM.app.curFrame(BM.app.state.active).xyz));
  const [mx, my] = await atomXY(3);
  check((await pick(mx, my)).index === 3, 'second water O pickable');
  await drag(mx, my, mx + 160, my - 60);
  const after = await page.evaluate(() => Array.from(BM.app.curFrame(BM.app.state.active).xyz));
  const moved = k => Math.hypot(after[3 * k] - before[3 * k], after[3 * k + 1] - before[3 * k + 1], after[3 * k + 2] - before[3 * k + 2]);
  check([3, 4, 5].every(k => moved(k) > 1.5) && [0, 1, 2].every(k => moved(k) < 1e-9), 'manipulate moves only the dragged fragment: ' + [0, 1, 2, 3, 4, 5].map(k => moved(k).toFixed(2)));
  const rigid = Math.abs(Math.hypot(after[9] - after[12], after[10] - after[13], after[11] - after[14]) - Math.hypot(before[9] - before[12], before[10] - before[13], before[11] - before[14]));
  check(rigid < 1e-9, 'fragment moved rigidly');
  const hb2 = await page.evaluate(() => BM.app.R().dashedBonds.length);
  check(hb2 === 0, 'H-bond disappears when the waters are pulled apart');
  await shot('09_manipulate');
  // shift-drag rotates the fragment about its centroid
  const [rx, ry] = await atomXY(3);
  await drag(rx, ry, rx + 50, ry + 30, { shift: true });
  const rot = await page.evaluate(() => Array.from(BM.app.curFrame(BM.app.state.active).xyz));
  const cen = a => [0, 1, 2].map(k => (a[9 + k] + a[12 + k] + a[15 + k]) / 3);
  const c1 = cen(after), c2 = cen(rot);
  check(Math.hypot(c1[0] - c2[0], c1[1] - c2[1], c1[2] - c2[2]) < 1e-6 && Math.hypot(rot[12] - after[12], rot[13] - after[13], rot[14] - after[14]) > 0.05, 'Shift-drag rotates the fragment about its centroid');
  // undo twice restores the original geometry
  await page.keyboard.press('Control+z'); await sleep(100);
  await page.keyboard.press('Control+z'); await sleep(120);
  const undone = await page.evaluate(() => Array.from(BM.app.curFrame(BM.app.state.active).xyz));
  check(undone.every((v, k) => Math.abs(v - before[k]) < 1e-9), 'undo restores the previous geometry');
  check((await page.evaluate(() => BM.app.R().dashedBonds.length)) === 1, 'H-bond back after undo');
  await page.keyboard.press('Control+Shift+z'); await sleep(100);
  const redone = await page.evaluate(() => Array.from(BM.app.curFrame(BM.app.state.active).xyz));
  check(redone.every((v, k) => Math.abs(v - after[k]) < 1e-9), 'redo re-applies the move');

  // ---------------------------------------------------------------- 8. select tool rubber band
  await page.evaluate(id => { const a = BM.app; a.activate(a.state.docs.find(d => d.id === id)); }, butene);
  await page.keyboard.press('s');
  await page.keyboard.press('Escape');
  const pts = await Promise.all((await page.evaluate(() => BM.app.curFrame(BM.app.state.active).Z.map((_, i) => i))).map(atomXY));
  const xsorted = pts.map(p => p[0]).sort((a, b) => a - b), ys = pts.map(p => p[1]);
  const cut = (xsorted[Math.floor(xsorted.length / 2) - 1] + xsorted[Math.floor(xsorted.length / 2)]) / 2;
  const x0 = xsorted[0] - 30, y0 = Math.min(...ys) - 30, y1 = Math.max(...ys) + 30;
  const [bx0, by0] = await emptyNear(x0, y0);
  await page.mouse.move(bx0, by0); await page.mouse.down();
  await page.mouse.move(cut, y1, { steps: 8 });
  await shot('10_select_band');
  await page.mouse.up(); await sleep(100);
  const expect = await page.evaluate(([a, b, c, d]) => BM.app.R().atomsInRect(a, b, c, d).sort((p, q) => p - q), [bx0 - box.x, by0 - box.y, cut - box.x, y1 - box.y]);
  s = await info();
  check(s.sel.length > 0 && s.sel.length < s.n && JSON.stringify(s.sel.slice().sort((p, q) => p - q)) === JSON.stringify(expect), 'rubber band selects the atoms inside: ' + s.sel.length + ' of ' + s.n);
  const chipSel = await page.textContent('#measure');
  check(s.sel.length <= 4 || /atoms selected/.test(chipSel), 'chip shows "N atoms selected": ' + chipSel);
  // shift adds, ctrl removes, double-click selects the whole fragment
  const firstSel = s.sel.length;
  await page.click('.toolopts [data-sel="invert"]');
  s = await info();
  check(s.sel.length === s.n - firstSel, 'Invert selection');
  await page.keyboard.press('Control+a');
  s = await info();
  check(s.sel.length === s.n, 'Ctrl+A selects all');
  await page.click('.toolopts [data-sel="none"]');
  const hvs = await heavy();
  const [sx1, sy1] = await atomXY(hvs[0]), [sx2, sy2] = await atomXY(hvs[1]);
  await page.mouse.click(sx1, sy1); await sleep(400);
  await page.keyboard.down('Shift'); await page.mouse.click(sx2, sy2); await page.keyboard.up('Shift'); await sleep(100);
  s = await info();
  check(s.sel.length === 2 && s.sel.indexOf(hvs[0]) >= 0 && s.sel.indexOf(hvs[1]) >= 0, 'click then Shift-click selects two atoms: ' + s.sel);
  await page.click('.toolopts [data-sel="none"]');
  await sleep(400);
  const [ax, ay] = await atomXY((await heavy())[0]);
  await page.mouse.dblclick(ax, ay); await sleep(120);
  s = await info();
  check(s.sel.length === s.n, 'double-click selects the fragment (' + s.sel.length + ')');
  await shot('11_select_all');

  // auto-optimize while drawing relaxes the new atoms (one undo step per change)
  await page.click('#btn-new');
  await page.check('.toolopts [data-opt="autoOpt"]');
  await page.mouse.click(box.x + box.width / 2, box.y + box.height / 2); await sleep(120);
  const [cx0, cy0] = await atomXY(0);
  await drag(cx0, cy0, cx0 + 40, cy0 - 140, { wait: 40 });
  await waitIdle();
  const auto = await page.evaluate(() => { const d = BM.app.state.active, f = BM.app.curFrame(d), c = f.Z.map((z, i) => (z === 6 ? i : -1)).filter(i => i >= 0); return { formula: BM.geom.formula(f.Z), cc: BM.geom.dist(f.xyz, c[0], c[1]), undo: d.undo.length, e: BM.ff.optimize({ Z: f.Z, xyz: f.xyz, bonds: f.bonds }, { maxSteps: 0 }).energy, e2: BM.ff.optimize({ Z: f.Z, xyz: f.xyz, bonds: f.bonds }).energy }; });
  check(auto.formula === 'C2H6' && Math.abs(auto.cc - 1.52) < 0.03 && auto.undo === 2 && auto.e - auto.e2 < 0.05, 'auto-optimize relaxes ethane while drawing: ' + JSON.stringify(auto));
  await page.uncheck('.toolopts [data-opt="autoOpt"]');
  // perceived aromatic rings for a file without bond orders: benzene.xyz keeps six equal C–C bonds after clean-up
  const bzx = ['12', 'benzene', 'C 0.000 1.397 0', 'C 1.2098 0.6985 0', 'C 1.2098 -0.6985 0', 'C 0 -1.397 0', 'C -1.2098 -0.6985 0', 'C -1.2098 0.6985 0',
    'H 0 2.484 0', 'H 2.1512 1.242 0', 'H 2.1512 -1.242 0', 'H 0 -2.484 0', 'H -2.1512 -1.242 0', 'H -2.1512 1.242 0'].join('\n');
  await page.evaluate(t => BM.app.openText('benzene.xyz', t), bzx);
  const bzOrders = await page.evaluate(() => Array.from(BM.app.ordersFor(BM.app.state.active, 0)).filter(o => o === 1.5).length);
  check(bzOrders === 6, 'benzene.xyz: six aromatic bonds perceived for display (' + bzOrders + ')');
  await page.click('#tabs [data-tab="build"]');
  await page.click('[data-act="b-opt"]');
  await waitIdle();
  const bzLens = await page.evaluate(() => { const f = BM.app.curFrame(BM.app.state.active); return f.bonds.filter(b => f.Z[b[0]] === 6 && f.Z[b[1]] === 6).map(b => BM.geom.dist(f.xyz, b[0], b[1])); });
  check(bzLens.length === 6 && Math.max(...bzLens) - Math.min(...bzLens) < 0.01, 'clean-up of benzene.xyz keeps equal C–C bonds: ' + bzLens.map(r => r.toFixed(3)));
  check(await page.evaluate(() => BM.app.curFrame(BM.app.state.active).label) === 'Edited', 'clean-up of a file geometry works on an "Edited" copy');

  // ---------------------------------------------------------------- 9. calculation input
  console.log('Calculation input');
  await page.click('#dd-export [data-dd]');
  await page.click('[data-act="export-calc"]');
  await page.waitForSelector('#dlg[open]');
  const progs = await page.evaluate(() => BM.inputPrograms.map(p => p.id));
  const bad = [];
  for (const id of progs) {
    await page.selectOption('[data-c="program"]', id);
    await sleep(40);
    const r = await page.evaluate(() => ({ text: document.querySelector('[data-c="preview"]').value, err: document.querySelector('[data-c="note"]').classList.contains('err'), note: document.querySelector('[data-c="note"]').textContent, save: document.querySelector('[data-dbtn="2"]').textContent }));
    if (!r.text.trim() || r.err) bad.push(id + ': ' + r.note);
    if (id === 'orca') await shot('12_calc_dialog_orca');
  }
  check(bad.length === 0, 'every program writes a non-empty input: ' + (bad.join(' | ') || progs.join(', ')));
  await page.selectOption('[data-c="program"]', 'xtb');
  await page.selectOption('[data-c="task"]', 'opt');
  const xtbNote = await page.textContent('[data-c="note"]');
  check(/Run: xtb .*--opt/.test(xtbNote), 'xtb shows the command: ' + xtbNote);
  await page.evaluate(() => { window.__saved.length = 0; });
  await page.click('#dlg .btn.primary');
  await sleep(150);
  let saved = await page.evaluate(() => window.__saved.map(x => x.name));
  check(saved.length === 1 && /\.xyz$/.test(saved[0]), 'xtb saves its coordinate file: ' + saved);
  // charge default for built molecules = sum of formal charges; Gaussian route-from-file is off for built docs
  await page.click('#btn-new');
  await page.fill('#bp-smiles', 'C[N+](C)(C)C');
  await page.click('[data-act="b-smiles"]');
  await sleep(150);
  await page.evaluate(() => BM.app.exportCalc('gaussian'));
  await page.waitForSelector('#dlg[open]');
  const ch = await page.inputValue('[data-c="charge"]');
  check(ch === '1', 'built tetramethylammonium: charge defaults to +1 (' + ch + ')');
  check(!(await page.$('[data-c="useRoute"]')), 'no route checkbox for a built molecule');
  await page.click('#dlg .btn.primary'); await sleep(100);
  saved = await page.evaluate(() => window.__saved.map(x => x.name));
  check(saved.some(n => /^New_molecule.*\.gjf$/.test(n)), 'Gaussian input saved for the built molecule: ' + saved);

  // ---------------------------------------------------------------- 10. file documents: Build tab and extras
  console.log('File documents');
  await page.setInputFiles('#file-input', path.join(ROOT, 'tests', 'fixtures', 'g16_water_optfreq.log'));
  await page.waitForFunction(() => BM.app.state.active.name === 'g16_water_optfreq.log');
  await page.click('#tabs [data-tab="build"]');
  await page.evaluate(() => BM.app.exportCalc());
  await page.waitForSelector('#dlg[open]');
  const useRoute = await page.$eval('[data-c="useRoute"]', e => e.checked);
  const prevG = await page.inputValue('[data-c="preview"]');
  check(useRoute && /^#/m.test(prevG) && /0 1/.test(prevG), 'Gaussian file: route from the file is used by default');
  await page.click('#dlg .btn:not(.primary):nth-child(2)');
  await page.click('#dd-settings [data-dd]');
  await page.check('#set-dipole');
  await page.click('#status');
  await sleep(80);
  const arrows = await page.evaluate(() => BM.app.R().extraArrows.length);
  check(arrows === 1, 'dipole arrow shown for a file with a dipole');
  await page.keyboard.press('s');
  await page.keyboard.press('Control+a');
  await page.keyboard.press('Delete');
  s = await info();
  check(s.n === 3, 'Delete on every atom of a file geometry is refused');
  await page.keyboard.press('Escape');
  await shot('13_water_file_build_tab');
  // light theme
  await page.click('#btn-theme');
  await page.evaluate(id => { const a = BM.app; a.activate(a.state.docs.find(d => d.id === id)); }, butene);
  await page.keyboard.press('d');
  await page.click('.toolopts [data-pt="draw"]');
  await sleep(60);
  await shot('14_light_periodic_table');
  await page.click('.popup [data-z="7"]');
  check(await page.evaluate(() => BM.app.state.build.element) === 7, 'periodic table picks N');
  await page.click('.toolopts [data-el="6"]');
  await page.click('#btn-theme');

  // ---------------------------------------------------------------- 11. large structure stays responsive
  console.log('Large structure');
  const gro = (() => {
    const L = ['water box', '5001'], n = 17;
    let k = 0;
    for (let a = 0; a < n && k < 1667; a++) for (let b = 0; b < n && k < 1667; b++) for (let c = 0; c < 6 && k < 1667; c++, k++) {
      const x = a * 0.31, y = b * 0.31, z = c * 0.31, r = k + 1;
      L.push(String(r % 100000).padStart(5) + 'SOL'.padEnd(5) + '   OW' + String((3 * k + 1) % 100000).padStart(5) + x.toFixed(3).padStart(8) + y.toFixed(3).padStart(8) + z.toFixed(3).padStart(8));
      L.push(String(r % 100000).padStart(5) + 'SOL'.padEnd(5) + '  HW1' + String((3 * k + 2) % 100000).padStart(5) + (x + 0.0957).toFixed(3).padStart(8) + y.toFixed(3).padStart(8) + z.toFixed(3).padStart(8));
      L.push(String(r % 100000).padStart(5) + 'SOL'.padEnd(5) + '  HW2' + String((3 * k + 3) % 100000).padStart(5) + (x - 0.024).toFixed(3).padStart(8) + (y + 0.0927).toFixed(3).padStart(8) + z.toFixed(3).padStart(8));
    }
    L.push('   5.27000   5.27000   1.86000');
    return L.join('\n') + '\n';
  })();
  const tOpen = await page.evaluate(t => { const t0 = performance.now(); BM.app.openText('waterbox.gro', t); return performance.now() - t0; }, gro);
  s = await info();
  check(s.n === 5001, 'water box with 5001 atoms opens (' + s.n + ')');
  check(tOpen < 4000, 'opening 5001 atoms takes ' + tOpen.toFixed(0) + ' ms');
  const tUi = await page.evaluate(() => {
    const a = BM.app, t0 = performance.now();
    a.state.active.tab = 'build'; document.querySelector('#tabs [data-tab="build"]') && document.querySelector('#tabs [data-tab="build"]').click();
    for (const t of ['select', 'draw', 'template', 'bond', 'manipulate', 'navigate']) a.setTool(t);
    a.selectAll(); a.ui().selectAction('none');
    return performance.now() - t0;
  });
  check(tUi < 1500, 'Build tab, tool switching and select-all on 5001 atoms: ' + tUi.toFixed(0) + ' ms');
  await page.keyboard.press('k');
  const tHover = await page.evaluate(() => {
    const cv = BM.app.R().canvas, r = cv.getBoundingClientRect(), t0 = performance.now();
    for (let k = 0; k < 20; k++) cv.dispatchEvent(new PointerEvent('pointermove', { clientX: r.left + 300 + 15 * k, clientY: r.top + 300 + 8 * k, bubbles: true }));
    return (performance.now() - t0) / 20;
  });
  check(tHover < 60, 'bond-tool hover (atom + bond picking, highlight) on 5001 atoms: ' + tHover.toFixed(1) + ' ms per move');
  const tHb = await page.evaluate(() => { const t0 = performance.now(); BM.app.setView({ hbonds: true }); return performance.now() - t0; });
  check(tHb < 1500, 'H-bond search on 5001 atoms: ' + tHb.toFixed(0) + ' ms');
  console.log(`  5001 atoms: open ${tOpen.toFixed(0)} ms · Build tab + tools + select-all ${tUi.toFixed(0)} ms · hover ${tHover.toFixed(1)} ms/move · H-bonds ${tHb.toFixed(0)} ms`);
  await page.keyboard.press('n');
  await shot('15_large_waterbox');

  console.log('Console/page errors:', errors.length ? errors : 'none');
  check(errors.length === 0, 'no console errors');
  await browser.close();
  console.log(`\nBuild UI tests: ${passes} passed, ${failures} failed`);
  process.exit(failures ? 1 : 0);
})().catch(e => { console.error(e); process.exit(2); });
