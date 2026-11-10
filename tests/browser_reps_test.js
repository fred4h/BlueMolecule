// End-to-end checks of representations in dist/BlueMolecule.html (headless Chromium): presets chosen on load for the
// QM/MM fixtures (layer styles, legend, link markers) and for synthetic protein, membrane, glycan and nucleic-acid
// systems built by tests/synth.js; the Display tab (rows, live count, errors, autocomplete, presets); the toolbar on
// the "All atoms" row; Info-tab chips; frame changes and edits; timings on a ~290k-atom membrane; screenshots.
// Usage: node tests/browser_reps_test.js [screenshotDir]
const { chromium } = require('playwright');
const path = require('path');
const fs = require('fs');

const ROOT = path.resolve(__dirname, '..');
const HTML = 'file://' + path.join(ROOT, 'dist', 'BlueMolecule.html');
const FIX = path.join(ROOT, 'tests', 'fixtures');
const SHOTS = process.argv[2] || path.join(ROOT, 'tests', 'shots', 'reps');
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
  await page.addScriptTag({ path: path.join(ROOT, 'tests', 'synth.js') });
  const ev = (fn, arg) => page.evaluate(fn, arg);
  const shot = name => page.screenshot({ path: path.join(SHOTS, name + '.png') });
  const settle = () => page.waitForTimeout(250);
  const state = () => ev(() => {
    const a = BM.app, d = a.state.active, R = a.R();
    return { name: d.name, preset: d.repPreset, rows: d.reps.map(r => r.selection + ':' + r.style), counts: d._repRes.counts, errors: d._repRes.errors, rep: !!R.rep, nSphere: R.nSphere, nCyl: R.nCyl,
      nLine: R.nLine, nGhost: R.nGhostSph + R.nGhostCyl, nMarker: R.nMarker, overlays: R.overlays.length, verts: R.overlays.reduce((s, o) => s + o.vertices, 0), proxies: R._proxy ? R._proxy.list.length : 0,
      legend: R.legend ? [R.legend.title].concat(R.legend.items.map(i => i.label + '|' + (i.sub || '') + '|' + i.glyph)) : null, sel: d.selection.length, n: a.curFrame(d).Z.length, tab: d.tab };
  });
  const openFixture = async file => {
    await page.setInputFiles('#file-input', path.join(FIX, file));
    await page.waitForFunction(n => { const s = window.BM.app.state; return s.active && s.active.name === n; }, path.basename(file), { timeout: 15000 });
    await settle();
  };
  const openSynth = async (name, code) => { await ev(([nm, c]) => { const S = window.BMSynth, sys = new Function('S', 'BM', 'return ' + c)(S, BM); BM.app.openText(nm, /\.gro$/.test(nm) ? S.toGRO(sys) : S.toPDB(sys, { ss: true })); }, [name, code]); await settle(); };
  const tab = async id => { await page.click(`#tabs [data-tab="${id}"]`); await page.waitForTimeout(120); };
  const chips = () => ev(() => [...document.querySelectorAll('#panel .ct .ct-n')].map(e => e.textContent));
  /** Frame a selection for a screenshot (centre on it, distance from its radius). */
  const zoomOn = (sel, k) => ev(([sel, k]) => {
    const a = BM.app, R = a.R(), m = BM.reps.evaluate(sel, a.repCtx(a.state.active)), x = R.xyz, c = [0, 0, 0];
    let n = 0, r = 0;
    for (let i = 0; i < m.length; i++) if (m[i]) { c[0] += x[3 * i]; c[1] += x[3 * i + 1]; c[2] += x[3 * i + 2]; n++; }
    c[0] /= n; c[1] /= n; c[2] /= n;
    for (let i = 0; i < m.length; i++) if (m[i]) r = Math.max(r, Math.hypot(x[3 * i] - c[0], x[3 * i + 1] - c[1], x[3 * i + 2] - c[2]));
    R.center = c; R.dist = Math.max(12, r * k); R.requestRender();
  }, [sel, k]);

  // ---------------------------------------------------------------- 1. QM/MM and ONIOM outputs
  console.log('QM/MM layers preset');
  const qmmm = ['g16_oniom3_opt.log', 'g16_oniom_amber_opt.log', 'orca6_qmmm_opt.out', 'orca5_qmxtb_opt.out', 'cp2k_qmmm_geoopt.out', 'xtb_oniom_water.out'];
  for (const f of qmmm) {
    const before = errors.length;
    await openFixture('qmmm/' + f);
    const s = await state(), L = await ev(() => { const d = BM.app.state.active, L = d.data.layers, c = {}; for (const x of L.atomLayer) c[x] = (c[x] || 0) + 1; return { counts: c, labels: L.labels, links: L.linkAtoms.length, frozen: L.frozen ? L.frozen.filter(Boolean).length : 0, scheme: L.scheme }; });
    const nLayers = Object.keys(L.counts).length;
    check(s.preset === 'qmmm' && s.rep, f + ': QM/MM preset with per-atom styles ' + JSON.stringify(s.rows));
    check(s.rows[0] === 'layer H:ballstick' && (!L.counts.M || s.rows.includes('layer M:sticks')) && (!L.counts.L || s.rows.includes('layer L:lines')), f + ': H ball-and-stick, M sticks, L lines');
    check(s.nSphere > 0 && (!L.counts.L || s.nLine > 0), f + ': spheres for the QM atoms, lines for MM ' + [s.nSphere, s.nLine]);
    check(s.legend && s.legend.length >= 1 + nLayers, f + ': legend ' + JSON.stringify(s.legend));
    for (const k of Object.keys(L.counts)) if (L.labels[k]) check(s.legend.some(x => x.includes(L.labels[k]) && x.includes(L.counts[k] + ' atom')), f + ': legend names layer ' + k + ' (' + L.labels[k] + ')');
    if (L.links) check(s.nMarker >= 1 && s.legend.some(x => /^Link atoms/.test(x)), f + ': link atoms marked ' + s.nMarker);
    if (L.frozen) check(s.rows.some(r => /^frozen/.test(r)) && s.legend.some(x => /^Frozen/.test(x)), f + ': frozen atoms dimmed and in the legend');
    await tab('info');
    const c = await chips(), want = L.scheme === 'oniom' ? ['High layer', 'Low layer'] : ['QM region'];
    check(want.every(w => c.includes(w)), f + ': layer chips in Contents ' + JSON.stringify(c));
    check(errors.length === before, f + ': no errors ' + errors.slice(before).join(' | '));
  }
  await ev(() => { const a = BM.app; a.activate(a.state.docs.find(d => d.name === 'g16_oniom3_opt.log')); });
  await settle();
  await shot('oniom_layers');
  // a layer chip selects that layer (and keeps the QM/MM preset)
  await tab('info');
  await page.click('#panel .ct:has-text("Medium layer")');
  await settle();
  let s = await state();
  check(s.sel === 3 && s.preset === 'qmmm', 'Medium-layer chip selects its 3 atoms ' + s.sel);
  await ev(() => { const a = BM.app; a.state.active.selection = []; a.R().setSelection([]); });

  // ---------------------------------------------------------------- 2. protein with a ligand
  console.log('Protein cartoon + ligand');
  await openSynth('protein.pdb', 'S.protein(BM)');
  s = await state();
  check(s.preset === 'protein', 'protein-ligand preset chosen on load ' + s.preset);
  check(s.rows.join(' | ') === 'protein:cartoon | same residue as (protein and within 5 of ligand):sticks | ligand:ballstick | ions:spacefill | water:hidden', 'protein-ligand rows ' + s.rows.join(' | '));
  check(s.overlays >= 1 && s.verts > 2000 && s.proxies > 60, 'cartoon mesh and CA pick proxies ' + [s.overlays, s.verts, s.proxies]);
  check(s.counts[2] === 30 && s.counts[4] > 500 && s.errors.every(e => !e), 'per-row counts ' + JSON.stringify(s.counts));
  await zoomOn('protein or ligand', 3.9);
  await settle();
  await shot('protein_ligand');
  // hovering the cartoon reports the residue under it
  const hov = await ev(() => {
    const a = BM.app, R = a.R(), f = a.curFrame(a.state.active);
    const ca = f.atomNames.findIndex((n, i) => n === 'CA' && f.resIds[i] === 8), p = R.worldToScreen([R.xyz[3 * ca], R.xyz[3 * ca + 1], R.xyz[3 * ca + 2]]);
    const box = R.canvas.getBoundingClientRect();
    const k = R.pick(p.x, p.y);
    return { x: box.left + p.x, y: box.top + p.y, ca, pick: k, name: k >= 0 ? f.atomNames[k] : null, res: k >= 0 ? f.resNames[k] + ' ' + f.resIds[k] : null };
  });
  await page.mouse.move(hov.x, hov.y);
  await page.waitForTimeout(150);
  const hoverText = await ev(() => document.querySelector('#status-hover').textContent);
  check(hov.name === 'CA' && hoverText.includes(hov.res) && /helix/.test(hoverText), 'a CA under the cartoon is pickable; status shows its residue ' + JSON.stringify([hov.name, hov.res, hoverText]));
  await tab('info');
  let c = await chips();
  check(['Protein', 'Ligands', 'Water', 'Ions'].every(w => c.includes(w)), 'residue-class chips ' + JSON.stringify(c));
  await page.click('#panel .ct:has-text("Water")');
  await settle();
  s = await state();
  check(s.sel === s.counts[4], 'Water chip selects every water atom ' + s.sel);
  await ev(() => { const a = BM.app; a.state.active.selection = []; a.R().setSelection([]); a.renderPanel(); });

  // ---------------------------------------------------------------- 3. Display tab
  console.log('Display tab');
  await tab('display');
  check(await page.locator('#rp-list .rp').count() === 5, 'five rows listed');
  await page.click('#rp-list .rp:nth-child(3)');
  await page.waitForTimeout(100);
  check(await page.inputValue('#rp-q') === 'ligand' && /30 atoms/.test(await page.textContent('#rp-qn')), 'row selected: editor shows its selection and live count');
  await page.fill('#rp-q', 'ligand and');
  await page.waitForTimeout(300);
  check(/"and" needs a selection after it/.test(await page.textContent('#rp-err')) && (await ev(() => BM.app.state.active.reps[2].selection)) === 'ligand', 'invalid selection: error shown, row unchanged');
  await page.fill('#rp-q', 'resname L');
  await page.waitForTimeout(80);
  const ac = await ev(() => [...document.querySelectorAll('#rp-ac .ac-it .mono')].map(e => e.textContent));
  await shot('display_tab');
  check(ac.includes('LIG') && ac.includes('LEU'), 'autocomplete offers residue names present ' + JSON.stringify(ac));
  await page.keyboard.press('ArrowDown');
  const pick = await ev(() => document.querySelector('#rp-ac .ac-it.on .mono').textContent);
  await page.keyboard.press('Enter');
  await page.waitForTimeout(300);
  check((await page.inputValue('#rp-q')).startsWith('resname ' + pick), 'Enter accepts the highlighted suggestion ' + await page.inputValue('#rp-q'));
  await page.fill('#rp-q', 'resname LIG or ions');
  await page.waitForTimeout(350);
  s = await state();
  check(s.rows[2] === 'resname LIG or ions:ballstick' && s.counts[2] === 32 && /32 atoms/.test(await page.textContent('#rp-qn')), 'valid edit applied live ' + s.rows[2] + ' ' + s.counts[2]);
  await page.selectOption('#rp-style', 'spacefill');
  await settle();
  const sp = await state();
  check(sp.rows[2] === 'resname LIG or ions:spacefill' && sp.nSphere >= s.nSphere, 'style change applied');
  await page.selectOption('#rp-color', 'uniform');
  await settle();
  check(/^#[0-9a-f]{6}$/.test(await ev(() => BM.app.state.active.reps[2].color)) && await page.locator('#rp-uni').count() === 1, 'uniform colour with a colour picker');
  await page.click('[data-act="rp-add"]');
  await page.waitForTimeout(150);
  check((await state()).rows.length === 6 && await ev(() => document.activeElement && document.activeElement.id) === 'rp-q', '+ Add inserts a row and focuses its selection');
  await page.click('[data-act="rp-dup"]'); await page.waitForTimeout(100);
  check((await state()).rows.length === 7, 'duplicate');
  await page.click('[data-act="rp-del"]'); await page.waitForTimeout(100);
  await page.click('[data-act="rp-del"]'); await page.waitForTimeout(100);
  check((await state()).rows.length === 5, 'delete');
  await page.click('#rp-list .rp:nth-child(2)'); await page.waitForTimeout(80);
  await page.click('[data-act="rp-up"]'); await page.waitForTimeout(100);
  s = await state();
  check(s.rows[0].startsWith('same residue as') && s.rows[1] === 'protein:cartoon', 'move up reorders ' + s.rows.slice(0, 2));
  await page.click('[data-act="rp-down"]'); await page.waitForTimeout(100);
  await page.click('#rp-list .rp:nth-child(1) .rp-eye'); await settle();
  check((await state()).overlays === 0, 'hiding the cartoon row removes the cartoon');
  await page.click('#rp-list .rp:nth-child(1) .rp-eye'); await settle();
  check((await state()).overlays >= 1, 'showing it again');
  // presets menu
  await page.click('#dd-presets [data-dd]');
  await page.click('#dd-presets [data-preset="all"]');
  await settle();
  s = await state();
  check(s.preset === 'all' && !s.rep && s.overlays === 0 && s.nSphere === s.n, 'All atoms preset = plain toolbar style ' + JSON.stringify([s.rep, s.nSphere, s.n]));
  await page.click('#dd-presets [data-dd]');
  await page.click('#dd-presets [data-preset="protein"]');
  await settle();
  check((await state()).overlays >= 1, 'Protein-ligand preset back');
  // toolbar style buttons act on the "All atoms" representation
  await page.click('#style-seg button[data-style="sticks"]');
  await settle();
  s = await state();
  check(s.preset === 'all' && s.rows.join() === 'all:sticks' && !s.rep && s.nCyl > 0, 'toolbar Sticks: all atoms as sticks ' + s.rows);
  await page.click('#style-seg button[data-style="ballstick"]');
  await settle();
  check((await state()).rows.join() === 'all:ballstick', 'toolbar Ball & stick updates the All-atoms row');
  // a QM/MM document keeps its own rows when another document changes style
  await ev(() => { const a = BM.app; a.activate(a.state.docs.find(d => d.name === 'g16_oniom3_opt.log')); });
  check((await state()).preset === 'qmmm', 'other documents keep their representations');

  // ---------------------------------------------------------------- 4. membrane
  console.log('Membrane (lipid style)');
  await openSynth('membrane.gro', 'S.membrane({ n: 10, protein: true })');
  s = await state();
  check(s.preset === 'membrane', 'membrane preset chosen on load ' + s.preset);
  check(['protein:cartoon', 'lipid:lipid', 'ions:spacefill', 'water:hidden'].every(r => s.rows.includes(r)), 'membrane rows ' + s.rows.join(' | '));
  const lip = await ev(() => {
    const a = BM.app, d = a.state.active, R = a.R(), f = a.curFrame(d), st = R.rep.style, out = { P: 0, Psph: 0, tail: 0, tailLine: 0, H: 0, Hhid: 0, wat: 0, watHid: 0 };
    for (let i = 0; i < f.Z.length; i++) {
      const r = f.resNames[i];
      if (r === 'POPC') {
        if (f.Z[i] === 15) { out.P++; if (st[i] === 3) out.Psph++; }
        else if (f.Z[i] === 6) { out.tail++; if (st[i] === 5) out.tailLine++; }
        else if (f.Z[i] === 1) { out.H++; if (st[i] === 0) out.Hhid++; }
      } else if (r === 'TIP3') { out.wat++; if (st[i] === 0) out.watHid++; }
    }
    return out;
  });
  check(lip.P > 100 && lip.Psph === lip.P && lip.tailLine === lip.tail && lip.Hhid === lip.H && lip.watHid === lip.wat, 'lipid style: P spheres, carbon tails as lines, hydrogens and water hidden ' + JSON.stringify(lip));
  check(s.nLine > 5000 && s.overlays >= 1, 'lines and the protein cartoon ' + [s.nLine, s.overlays]);
  await tab('info');
  check((await chips()).includes('Lipids'), 'Lipids chip');
  await shot('membrane');
  // trajectory frame: the rows persist and the overlays follow the new coordinates
  const fr = await ev(() => {
    const a = BM.app, d = a.state.active, f = d.data.frames[0], R = a.R(), ov0 = R.overlays[0];
    d.data.frames.push(Object.assign({}, f, { xyz: f.xyz.map((v, i) => v + (i % 3 === 2 ? 1.5 : 0)), label: 'frame 2' }));
    a.setFrame(d, 1);
    return { frame: d.frameIndex, rows: d.reps.map(r => r.selection + ':' + r.style).join(' | '), rebuilt: R.overlays[0] !== ov0, verts: R.overlays[0].vertices, verts0: ov0 ? ov0.vertices : -1, counts: d._repRes.counts.join() };
  });
  check(fr.frame === 1 && fr.rows === s.rows.join(' | ') && fr.rebuilt && fr.verts === fr.verts0 && fr.counts === s.counts.join(), 'frame change keeps the rows; cartoon rebuilt for the new frame ' + JSON.stringify(fr));
  // timings on a ~290k-atom membrane (SwiftShader renders on the CPU: draw times are pessimistic)
  const perf = await ev(() => {
    const S = window.BMSynth, a = BM.app, sys = S.membrane({ n: 27, protein: true }), txt = S.toGRO(sys);
    let t = performance.now();
    a.openText('big_membrane.gro', txt);
    const open = performance.now() - t, d = a.state.active, R = a.R();
    t = performance.now(); a.refreshReps(d); const reps = performance.now() - t;
    t = performance.now(); R.updatePositions(R.xyz); const rebuild = performance.now() - t;
    t = performance.now(); R.render(); R.gl.finish(); const frame = performance.now() - t;
    const out = { atoms: sys.Z.length, open: Math.round(open), reapply: Math.round(reps), resolve: Math.round(d._repRes.ms), rebuild: Math.round(rebuild), draw: Math.round(frame), lines: R.nLine, spheres: R.nSphere, cartoonVerts: R.overlays.reduce((q, o) => q + o.vertices, 0) };
    document.querySelector(`[data-close="${d.id}"]`).click();
    return out;
  });
  console.log('  ~290k-atom membrane:', JSON.stringify(perf));
  check(perf.atoms > 250000 && perf.lines > 60000 && perf.open < 20000 && perf.reapply < 3000 && perf.rebuild < 1500, 'big membrane opens and re-renders within budget ' + JSON.stringify(perf));

  // ---------------------------------------------------------------- 5. heparan sulfate proteoglycan fragment
  console.log('Glycan: rings + SNFG');
  await openSynth('hspg.pdb', 'S.glycan(BM)');
  s = await state();
  check(s.preset === 'glycan', 'glycan preset chosen on load ' + s.preset);
  check(['protein:cartoon', 'carbohydrate:sticks', 'carbohydrate:rings', 'carbohydrate:snfg'].every(r => s.rows.includes(r)), 'glycan rows ' + s.rows.join(' | '));
  check(s.overlays >= 3, 'cartoon, filled rings and symbols ' + s.overlays);
  const sg = await ev(() => {
    const a = BM.app, R = a.R(), f = a.curFrame(a.state.active), i = f.atomNames.findIndex((n, k) => n === 'C1' && f.resNames[k] === 'HS2');
    const p = R.worldToScreen([R.xyz[3 * i], R.xyz[3 * i + 1], R.xyz[3 * i + 2]]), box = R.canvas.getBoundingClientRect();
    return { x: box.left + p.x, y: box.top + p.y };
  });
  await page.mouse.move(sg.x, sg.y);
  await page.waitForTimeout(150);
  const sgText = await ev(() => document.querySelector('#status-hover').textContent);
  check(/IdoA2S|GlcNS6S/.test(sgText), 'hovering a sugar names the monosaccharide (unknown residue name HS2) ' + sgText);
  await tab('info');
  check((await chips()).includes('Glycans'), 'Glycans chip');
  // look across the chain: screen x along the glycan, screen y from the protein towards it, then fit both
  await ev(() => {
    const a = BM.app, R = a.R(), c = a.repCtx(a.state.active), x = R.xyz;
    const cen = sel => { const m = BM.reps.evaluate(sel, c), p = [0, 0, 0]; let n = 0; for (let i = 0; i < m.length; i++) if (m[i]) { p[0] += x[3 * i]; p[1] += x[3 * i + 1]; p[2] += x[3 * i + 2]; n++; } return p.map(v => v / n); };
    const sub = (p, q) => [p[0] - q[0], p[1] - q[1], p[2] - q[2]], nrm = v => { const l = Math.hypot(...v); return v.map(t => t / l); };
    const u = nrm(sub(cen('resname IDS'), cen('resname XYS'))), g = cen('carbohydrate'), w0 = sub(g, cen('protein')), d = w0[0] * u[0] + w0[1] * u[1] + w0[2] * u[2];
    const w = nrm([w0[0] - d * u[0], w0[1] - d * u[1], w0[2] - d * u[2]]), n = [u[1] * w[2] - u[2] * w[1], u[2] * w[0] - u[0] * w[2], u[0] * w[1] - u[1] * w[0]];
    const m = [u, w, n], tr = m[0][0] + m[1][1] + m[2][2];
    let q;
    if (tr > 0) { const S = Math.sqrt(tr + 1) * 2; q = [(m[2][1] - m[1][2]) / S, (m[0][2] - m[2][0]) / S, (m[1][0] - m[0][1]) / S, 0.25 * S]; }
    else if (m[0][0] > m[1][1] && m[0][0] > m[2][2]) { const S = Math.sqrt(1 + m[0][0] - m[1][1] - m[2][2]) * 2; q = [0.25 * S, (m[0][1] + m[1][0]) / S, (m[0][2] + m[2][0]) / S, (m[2][1] - m[1][2]) / S]; }
    else if (m[1][1] > m[2][2]) { const S = Math.sqrt(1 + m[1][1] - m[0][0] - m[2][2]) * 2; q = [(m[0][1] + m[1][0]) / S, 0.25 * S, (m[1][2] + m[2][1]) / S, (m[0][2] - m[2][0]) / S]; }
    else { const S = Math.sqrt(1 + m[2][2] - m[0][0] - m[1][1]) * 2; q = [(m[0][2] + m[2][0]) / S, (m[1][2] + m[2][1]) / S, 0.25 * S, (m[1][0] - m[0][1]) / S]; }
    R.setOrientation(q);
  });
  await zoomOn('carbohydrate or protein', 3.9);
  await settle();
  await shot('glycan_snfg');
  // an edit (delete an ion): the rows survive and are re-evaluated on the new atom indices
  const ed = await ev(() => {
    const a = BM.app, d = a.state.active, f = a.curFrame(d), n0 = f.Z.length, ions0 = d._repRes.counts[d.reps.findIndex(r => r.selection === 'ions')];
    d.selection = [f.resNames.indexOf('SOD')];
    a.editDelete();
    const g = a.curFrame(d);
    return { n0, n1: g.Z.length, rows: d.reps.length, ions0, ions1: d._repRes.counts[d.reps.findIndex(r => r.selection === 'ions')], sugars: d._repRes.counts[d.reps.findIndex(r => r.style === 'snfg')] };
  });
  check(ed.n1 === ed.n0 - 1 && ed.ions1 === ed.ions0 - 1 && ed.sugars > 200 && ed.rows === s.rows.length, 'after deleting an ion the rows are re-evaluated ' + JSON.stringify(ed));

  // ---------------------------------------------------------------- 6. nucleic acid
  console.log('Nucleic acid');
  await openSynth('dna.pdb', 'S.nucleic(16)');
  s = await state();
  check(s.preset === 'nucleic' && s.rows[0] === 'nucleic:cartoon', 'nucleic preset ' + s.rows);
  check(s.overlays === 2 && s.proxies === 16, 'backbone tube + base slabs, one proxy per nucleotide ' + [s.overlays, s.proxies]);
  await shot('nucleic');

  console.log('Console/page errors:', errors.length ? errors : 'none');
  check(errors.length === 0, 'no console errors');
  await browser.close();
  console.log(`\nRepresentation browser tests: ${passes} passed, ${failures} failed  (screenshots in ${SHOTS})`);
  process.exit(failures ? 1 : 0);
})().catch(e => { console.error(e); process.exit(2); });
