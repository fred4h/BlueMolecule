// Browser checks of 2.5 in dist/BlueMolecule.html (headless Chromium): the Build tab's Crystal section (#21) and
// conformer search (#20).
//   - supercell, slab, cell lengths and angles (scaled or kept), wrap, remove the cell, a box around a molecule;
//     undo restores atoms and cell; the view shows the new cell; Save as writes it (CIF, POSCAR)
//   - conformer search on a built molecule: a new document, lowest first, cancel
// Usage: node tests/browser_build25_test.js
const { chromium } = require('playwright');
const path = require('path');
const fs = require('fs');

const ROOT = path.resolve(__dirname, '..');
const HTML = 'file://' + path.join(ROOT, 'dist', 'BlueMolecule.html');
const STRUCT = path.join(ROOT, 'tests', 'fixtures', 'struct');

let failures = 0, passes = 0;
const check = (cond, msg) => { if (cond) passes++; else { failures++; console.log('  FAIL:', msg); } };
const near = (a, b, t) => Math.abs(a - b) <= t;

(async () => {
  const browser = await chromium.launch({ args: ['--enable-unsafe-swiftshader'] });
  const page = await (await browser.newContext({ viewport: { width: 1440, height: 1000 } })).newPage();
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
  const info = () => page.evaluate(() => {
    const d = BM.app.state.active, f = BM.app.curFrame(d), c = BM.app.cellOf(d), R = BM.app.R();
    return { n: f.Z.length, edited: !!f.edited, cell: c, p: c ? BM.crystal.params(c) : null, shown: !!R.cell, status: document.querySelector('#status-msg').textContent, frames: d.data.frames.length };
  });
  const unfold = async id => { const folded = await page.$eval(`.bsec[data-sec="${id}"]`, e => e.classList.contains('folded')); if (folded) await page.click(`.bsec[data-sec="${id}"] h3`); };
  const setNum = async (id, v) => { await page.fill('#' + id, String(v)); };

  // ------------------------------------------------------------ rock salt
  await open('nacl_fm3m.cif', fs.readFileSync(path.join(STRUCT, 'nacl_fm3m.cif'), 'utf8'));
  await page.evaluate(() => BM.app.setView({ cell: true }));
  await page.click('#tabs button[data-tab="build"]');
  await unfold('crystal');
  const txt = await page.$eval('.bsec[data-sec="crystal"]', e => e.innerText);
  check(/a 5\.640\d · b 5\.640\d · c 5\.640\d Å/.test(txt) && /V 179\.4/.test(txt), 'lattice parameters shown: ' + txt.slice(0, 120));
  let s0 = await info();
  check(s0.n === 8 && s0.shown, 'eight atoms and the cell drawn');
  await setNum('bp-sc-1', 2); await setNum('bp-sc-2', 2); await setNum('bp-sc-3', 3);
  await page.click('[data-act="b-cell-super"]');
  let s1 = await info();
  check(s1.n === 96 && s1.edited && near(s1.p.c, 3 * s0.p.c, 1e-9) && near(s1.p.a, 2 * s0.p.a, 1e-9), 'supercell 2 × 2 × 3: ' + JSON.stringify({ n: s1.n, a: s1.p && s1.p.a }));
  check(s1.shown && /Supercell 2 × 2 × 3/.test(s1.status), 'the new cell drawn; status: ' + s1.status);
  // Save as: the writers get the new cell (written as POSCAR and CIF, read back)
  const written = await page.evaluate(() => {
    const d = BM.app.state.active, c = BM.app.saveContext(d, false), out = {};
    for (const id of ['poscar', 'cif']) {
      const f = BM.saveFormats.find(x => x.id === id);
      const text = f.write(c), ds = BM.parseFile(id === 'cif' ? 'x.cif' : 'POSCAR', typeof text === 'string' ? text : text.files[0].text);
      const cell = ds.cell || ds.frames[0].cell;
      out[id] = { n: ds.frames[0].Z.length, c: Math.hypot(...cell[2]) };
    }
    return out;
  });
  check(written.poscar.n === 96 && near(written.poscar.c, 3 * s0.p.c, 1e-6), 'POSCAR with the supercell: ' + JSON.stringify(written.poscar));
  check(written.cif.n === 96 && near(written.cif.c, 3 * s0.p.c, 1e-3), 'CIF with the supercell: ' + JSON.stringify(written.cif));
  // undo: atoms and cell back
  await page.keyboard.press('Control+z');
  await page.waitForTimeout(200);
  let s2 = await info();
  check(s2.n === 8 && near(s2.p.c, s0.p.c, 1e-9), 'undo restores the atoms and the cell: ' + JSON.stringify({ n: s2.n, c: s2.p && s2.p.c }));
  await page.keyboard.press('Control+Shift+z');
  await page.waitForTimeout(200);
  check((await info()).n === 96, 'redo');
  await page.keyboard.press('Control+z');
  await page.waitForTimeout(200);
  // slab (1 1 1)
  await unfold('crystal');
  await setNum('bp-hkl-h', 1); await setNum('bp-hkl-k', 1); await setNum('bp-hkl-l', 1); await setNum('bp-layers', 2); await setNum('bp-vac', 10);
  await page.click('[data-act="b-cell-slab"]');
  const s3 = await info();
  check(s3.n === 16 && near(s3.p.gamma, 120, 1e-6) && near(s3.p.c, 2 * s0.p.a / Math.sqrt(3) + 10, 1e-6), 'slab (111): ' + JSON.stringify({ n: s3.n, g: s3.p && s3.p.gamma, c: s3.p && s3.p.c }));
  check(/Slab \(1 1 1\): 2 layers, 10 Å vacuum/.test(s3.status), 'slab status: ' + s3.status);
  // the cell: lengths and angles, scaled or kept
  await unfold('crystal');
  const a0 = await page.$eval('#bp-cell-a', e => +e.value);
  await setNum('bp-cell-c', 30);
  if (await page.$eval('#bp-cell-scale', e => e.checked)) await page.click('#bp-cell-scale');
  const before = await page.evaluate(() => Array.from(BM.app.curFrame(BM.app.state.active).xyz));
  await page.click('[data-act="b-cell-set"]');
  const s4 = await info(), after = await page.evaluate(() => Array.from(BM.app.curFrame(BM.app.state.active).xyz));
  check(near(s4.p.c, 30, 1e-6) && near(s4.p.a, a0, 1e-3) && before.every((v, i) => near(v, after[i], 1e-9)), 'c = 30 Å, atoms kept in place');
  await unfold('crystal');
  await setNum('bp-cell-a', (a0 * 1.1).toFixed(4));
  await page.click('#bp-cell-scale');
  await page.click('[data-act="b-cell-set"]');
  const s5 = await info();
  check(near(s5.p.a, a0 * 1.1, 1e-3) && near(s5.p.c, 30, 1e-6), 'a lengthened by 10%, c kept');
  const span = await page.evaluate(() => { const f = BM.app.curFrame(BM.app.state.active); let lo = Infinity, hi = -Infinity; for (let i = 0; i < f.Z.length; i++) { lo = Math.min(lo, f.xyz[3 * i]); hi = Math.max(hi, f.xyz[3 * i]); } return hi - lo; });
  const span0 = (() => { let lo = Infinity, hi = -Infinity; for (let i = 0; i < after.length; i += 3) { lo = Math.min(lo, after[i]); hi = Math.max(hi, after[i]); } return hi - lo; })();
  check(near(span, span0 * 1.1, 0.05 * span0 + 1e-6), 'atoms scaled with the cell along x: ' + span.toFixed(3) + ' vs ' + span0.toFixed(3));
  // wrap, then remove the cell
  await unfold('crystal');
  await page.click('[data-act="b-cell-wrap"]');
  const fw = await page.evaluate(() => { const d = BM.app.state.active; return Array.from(BM.crystal.toFrac(BM.app.cellOf(d), BM.app.curFrame(d).xyz)); });
  check(fw.every(v => v >= -1e-9 && v < 1), 'wrapped into the cell');
  await unfold('crystal');
  await page.click('[data-act="b-cell-del"]');
  const s6 = await info();
  check(!s6.cell && !s6.shown, 'cell removed (not drawn, none for Save as)');
  check(await page.$eval('.bsec[data-sec="crystal"]', e => /no unit cell/.test(e.innerText) && !!e.querySelector('[data-act="b-cell-box"]')), 'the section offers a box');
  await page.keyboard.press('Control+z');
  await page.waitForTimeout(200);
  check(!!(await info()).cell, 'undo brings the cell back');

  // ------------------------------------------------------------ a molecule: a box, and conformers
  await page.evaluate(() => BM.app.newMolecule());
  await page.waitForFunction(() => document.querySelector('#bp-smiles'));
  await page.fill('#bp-smiles', 'CCCCO');
  await page.click('[data-act="b-smiles"]');
  await page.waitForFunction(() => BM.app.curFrame(BM.app.state.active).Z.length === 15);
  await unfold('crystal');
  await setNum('bp-box-m', 6);
  await page.click('[data-act="b-cell-box"]');
  const s7 = await info();
  check(s7.cell && s7.shown && s7.p.alpha === 90 && s7.p.gamma === 90, 'box around the molecule: ' + JSON.stringify(s7.p));
  const inside = await page.evaluate(() => { const d = BM.app.state.active; return Array.from(BM.crystal.toFrac(BM.app.cellOf(d), BM.app.curFrame(d).xyz)).every(v => v > 0.05 && v < 0.95); });
  check(inside, 'the molecule inside, away from the faces');
  await page.keyboard.press('Control+z');
  await page.waitForTimeout(200);
  await unfold('conf');
  const ctext = await page.$eval('.bsec[data-sec="conf"]', e => e.innerText);
  check(/2 rotatable bonds/.test(ctext), 'rotatable bonds counted: ' + ctext.slice(-160));
  check(!(await page.$('#panel .note')), 'no prose in the Build tab');
  const parent = await page.evaluate(() => BM.app.state.active.name);
  await page.click('[data-act="b-conf"]');
  await page.waitForFunction(p => BM.app.state.active && BM.app.state.active.name === p + ' conformers', parent, { timeout: 60000 });
  const doc = await page.evaluate(() => { const d = BM.app.state.active; return { frames: d.data.frames.length, labels: d.data.frames.map(f => f.label), e: d.data.frames.map(f => f.energy), warn: d.data.warnings.join(' '), tabs: Array.from(document.querySelectorAll('#tabs button')).map(b => b.dataset.tab) }; });
  check(doc.frames >= 2 && /^Conformer 1 · ΔE 0\.00 kcal\/mol/.test(doc.labels[0]), 'conformers as a new document: ' + doc.labels.slice(0, 2).join(' | '));
  check(doc.e.every((v, i) => i === 0 || v >= doc.e[i - 1]), 'lowest first');
  check(/clean-up force field/.test(doc.warn) && doc.tabs.includes('steps'), 'reader note on the energies; a Steps tab');
  // cancel a search
  await page.evaluate(() => BM.app.newMolecule());
  await page.fill('#bp-smiles', 'CCCCCCCCCCCCCC');
  await page.click('[data-act="b-smiles"]');
  await page.waitForFunction(() => BM.app.curFrame(BM.app.state.active).Z.length > 40);
  await unfold('conf');
  const docs0 = await page.evaluate(() => BM.app.state.docs.length);
  await page.click('[data-act="b-conf"]');
  await page.waitForSelector('[data-act="b-conf-cancel"]', { timeout: 5000 });
  await page.click('[data-act="b-conf-cancel"]');
  await page.waitForTimeout(300);
  check((await page.evaluate(() => BM.app.state.docs.length)) === docs0 && !!(await page.$('[data-act="b-conf"]')), 'cancelled: no new document, Search again');

  check(!errors.length, 'no page errors: ' + errors.join(' | '));
  await browser.close();
  console.log(`Building tools browser tests: ${passes} passed, ${failures} failed`);
  process.exit(failures ? 1 : 0);
})().catch(e => { console.error(e); process.exit(1); });
