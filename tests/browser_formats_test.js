// End-to-end checks of file formats, the Contents list, the plot picker (convergence, scan/IRC, ECD, NMR, MO levels,
// charges), Save As (with a mocked Open Babel) and MD topology + trajectory loading in dist/BlueMolecule.html.
// Usage: node tests/browser_formats_test.js [screenshotDir]
const { chromium } = require('playwright');
const path = require('path');
const fs = require('fs');

const ROOT = path.resolve(__dirname, '..');
const HTML = 'file://' + path.join(ROOT, 'dist', 'BlueMolecule.html');
const FIX = path.join(ROOT, 'tests', 'fixtures');
const SHOTS = process.argv[2] || path.join(ROOT, 'tests', 'shots', 'formats');
fs.mkdirSync(SHOTS, { recursive: true });

let failures = 0, passes = 0;
const check = (cond, msg) => { if (cond) passes++; else { failures++; console.log('  FAIL:', msg); } };

/** A small little-endian CHARMM DCD file (no unit cell) written for this test. */
function dcd(natoms, frames) {
  const i32 = v => { const b = Buffer.alloc(4); b.writeInt32LE(v); return b; };
  const rec = body => Buffer.concat([i32(body.length), body, i32(body.length)]);
  const head = Buffer.alloc(84);
  head.write('CORD', 0, 'latin1');
  head.writeInt32LE(frames.length, 4); head.writeInt32LE(1, 12); head.writeInt32LE(24, 80);
  head.writeFloatLE(0.002, 40);
  const title = Buffer.alloc(84, 32); title.writeInt32LE(1, 0); title.write('BlueMolecule test trajectory', 4, 'latin1');
  const parts = [rec(head), rec(title), rec(i32(natoms))];
  for (const f of frames) for (let d = 0; d < 3; d++) {
    const b = Buffer.alloc(4 * natoms);
    for (let k = 0; k < natoms; k++) b.writeFloatLE(f[3 * k + d], 4 * k);
    parts.push(rec(b));
  }
  return Buffer.concat(parts);
}

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
  const shot = name => page.screenshot({ path: path.join(SHOTS, name + '.png') });
  const toast = () => page.textContent('#toast');
  const docs = () => page.evaluate(() => window.BM.app.state.docs.length);
  /** Open a fixture (or {name, buffer}) and wait until the document count or the active document changes. */
  const open = async (f, buffer) => {
    const before = await page.evaluate(() => { const s = window.BM.app.state; return { n: s.docs.length, frames: s.active ? s.active.data.frames.length : 0, t: document.querySelector('#toast').textContent }; });
    await page.setInputFiles('#file-input', buffer ? { name: f, mimeType: 'application/octet-stream', buffer } : path.join(FIX, f));
    await page.waitForFunction(b => { const s = window.BM.app.state; return s.docs.length !== b.n || (s.active && s.active.data.frames.length !== b.frames) || document.querySelector('#toast').textContent !== b.t; }, before, { timeout: 15000 });
    await page.waitForTimeout(150);
  };
  const st = () => page.evaluate(() => {
    const a = window.BM.app, d = a.state.active, f = d.data.frames[d.frameIndex];
    return { name: d.name, format: d.data.format, frames: d.data.frames.length, fi: d.frameIndex, atoms: f.Z.length, pending: !!d.pending, plot: d.plotType, tab: d.tab,
      types: a.plotTypes(d).map(t => t[0]), nSphere: a.R().nSphere, tabs: [...document.querySelectorAll('#tabs button')].map(b => b.dataset.tab) };
  });
  const contents = () => page.evaluate(() => { const a = window.BM.app, d = a.state.active; d.tab = 'info'; document.querySelector('#tabs [data-tab="info"]').click(); return [...document.querySelectorAll('#panel .ct .ct-n')].map(e => e.textContent); });
  const spec = () => page.evaluate(() => { const s = window.BM.app.plot().spec; return s && JSON.parse(JSON.stringify(s, (k, v) => (v instanceof Float64Array ? Array.from(v) : typeof v === 'function' ? undefined : v))); });
  const pick = async id => { await page.selectOption('#plot-pick', id); await page.waitForTimeout(80); return spec(); };

  // ---------------------------------------------------------------- 1. quantum-chemistry outputs
  console.log('QC formats');
  const qc = [
    ['qc/gamess_water_opt_hess.log', 'gamess-out', ['Geometries', 'Vibrations', 'Charges', 'Orbital energies']],
    ['qc/nwchem_water_optfreq.out', 'nwchem-out', ['Geometries', 'Vibrations']],
    ['qc/qchem_water_optfreq.out', 'qchem-out', ['Geometries', 'Jobs', 'Vibrations']],
    ['qc/psi4_water_optfreq.out', 'psi4-out', ['Geometries', 'Vibrations']],
    ['qc/xtb_water_ohess.out', 'xtb-out', ['Geometries', 'Vibrations', 'Warnings']],
    ['qc/mopac_water_opt_force.out', 'mopac-out', ['Geometries', 'Vibrations']],
    ['qc/cp2k_water_geoopt.out', 'cp2k-out', ['Geometries', 'Unit cell']],
    ['qc/nwchem_water_tddft.out', 'nwchem-out', ['Excited states']],
    ['qc/coord', 'turbomole-coord', []]
  ];
  for (const [f, fmt, want] of qc) {
    const before = errors.length;
    await open(f);
    const s = await st();
    check(s.format === fmt && s.atoms > 0 && s.nSphere > 0, f + ' opens as ' + fmt + ': ' + JSON.stringify(s));
    const c = await contents();
    check(want.every(w => c.includes(w)), f + ' Contents ' + JSON.stringify(c) + ' has ' + want.join(', '));
    if (/opt/.test(f)) check(s.types.includes('opt:energy') && s.types.includes('opt:dE'), f + ' optimisation plots: ' + s.types.join(' '));
    check(errors.length === before, f + ' no errors: ' + errors.slice(before).join(' | '));
  }
  // cycles without geometries (xtb): gradient norm with threshold, no RMSD
  await page.evaluate(() => { const a = window.BM.app; a.activate(a.state.docs.find(d => d.name === 'xtb_water_ohess.out')); });
  let p = await pick('opt:gradNorm');
  check(p && p.logY && p.hlines.length === 1 && Math.abs(p.hlines[0].y - 0.001) < 1e-12 && p.points.length === 3, 'xtb gradient norm vs threshold ' + JSON.stringify(p && p.hlines));
  check(!(await st()).types.includes('opt:rmsd'), 'no RMSD plot without geometries');

  // ---------------------------------------------------------------- 2. structure formats
  console.log('Structure formats');
  const sf = [
    ['struct/nacl_fm3m.cif', 'cif', ['Unit cell']], ['struct/POSCAR', 'poscar', ['Unit cell']], ['struct/ligand.mol2', 'mol2', ['Charges']],
    ['struct/water.cjson', 'cjson', ['Vibrations', 'Charges']], ['struct/grid.xsf', 'xsf', ['Volumes', 'Unit cell']], ['struct/small.mmcif', 'mmcif', []],
    ['struct/traj.lammpstrj', 'lammps-dump', ['Geometries']], ['struct/full.data', 'lammps-data', ['Charges', 'Unit cell']], ['struct/relax.pwo', 'qe-out', ['Geometries']],
    ['struct/si.cell', 'castep-cell', ['Unit cell']], ['struct/array.cml', 'cml', []], ['struct/XDATCAR', 'xdatcar', ['Geometries']]
  ];
  for (const [f, fmt, want] of sf) {
    const before = errors.length;
    await open(f);
    const s = await st();
    check(s.format === fmt && s.atoms > 0 && s.nSphere > 0, f + ' opens as ' + fmt + ': ' + JSON.stringify(s));
    const c = await contents();
    check(want.every(w => c.includes(w)), f + ' Contents ' + JSON.stringify(c) + ' has ' + want.join(', '));
    check(errors.length === before, f + ' no errors: ' + errors.slice(before).join(' | '));
  }
  await page.evaluate(() => { const a = window.BM.app; a.activate(a.state.docs.find(d => d.name === 'full.data')); });
  const origin = await page.evaluate(() => window.BM.app.R().cellOrigin);
  check(origin && origin[1] === -1, 'LAMMPS box origin passed to the renderer: ' + JSON.stringify(origin));

  // ---------------------------------------------------------------- 3. topology, then coordinates and trajectories
  console.log('Topology + coordinates');
  await open('struct/water.prmtop');
  let s = await st();
  check(s.pending && s.atoms === 4 && !s.tabs.includes('build'), 'prmtop opens with coordinates pending ' + JSON.stringify(s));
  check(await page.isVisible('#pending'), 'pending hint shown in the view');
  const nDocs = await docs();
  await open('struct/water.inpcrd');
  s = await st();
  check(s.name === 'water.prmtop' && !s.pending && s.frames === 1 && s.nSphere === 4 && (await docs()) === nDocs, 'inpcrd fills the topology ' + JSON.stringify(s));
  check(!(await page.isVisible('#pending')), 'pending hint hidden');
  const base = await page.evaluate(() => Array.from(window.BM.app.state.active.data.frames[0].xyz));
  const traj = [0, 1, 2].map(k => base.map((v, i) => v + (i % 3 === 0 ? 0.05 * k : 0)));
  await open('water_md.dcd', dcd(4, traj));
  s = await st();
  check(s.name === 'water.prmtop' && s.frames === 4 && s.fi === 1 && (await docs()) === nDocs, 'DCD frames attached to the active document ' + JSON.stringify(s));
  const zOk = await page.evaluate(() => { const fr = window.BM.app.state.active.data.frames; return fr.every(f => f.Z.length === 4 && f.xyz.length === 12); });
  check(zOk, 'attached frames carry the topology atoms');
  await open('other.dcd', dcd(7, [new Array(21).fill(0)]));
  check(/open the topology\/structure first/.test(await toast()) && (await docs()) === nDocs, 'mismatched trajectory: clear message: ' + await toast());
  // regression: files handed over together by the macOS app (one message each) or the launcher (command-line
  // order) open topology first, like a multi-file drop
  const dcdB64 = dcd(4, traj).toString('base64'), topB64 = fs.readFileSync(path.join(FIX, 'struct/water.prmtop')).toString('base64');
  await page.evaluate(([a, b]) => { window.BM.native._receiveFile('host_md.dcd', a); window.BM.native._receiveFile('host.prmtop', b); }, [dcdB64, topB64]);
  await page.waitForFunction(n => window.BM.app.state.docs.length === n + 1 && window.BM.app.state.active.data.frames.length === 3, nDocs, { timeout: 5000 }).catch(() => {});
  s = await st();
  check(s.name === 'host.prmtop' && !s.pending && s.frames === 3 && (await docs()) === nDocs + 1, 'macOS app: a trajectory sent before its topology still attaches ' + JSON.stringify(s));
  const lp = await ctx.newPage();
  lp.on('pageerror', e => errors.push('launcher page: ' + e.message));
  await lp.addInitScript(files => { window.BM_LAUNCH = { files, version: 'test' }; },
    [{ name: 'l_md.dcd', url: 'data:application/octet-stream;base64,' + dcdB64 }, { name: 'l.prmtop', url: 'data:application/octet-stream;base64,' + topB64 }]);
  await lp.goto(HTML);
  await lp.waitForFunction(() => window.BM && window.BM.app && window.BM.app.state.active && window.BM.app.state.active.data.frames.length === 3, null, { timeout: 15000 }).catch(() => {});
  const ls = await lp.evaluate(() => { const d = window.BM.app && window.BM.app.state.active; return d ? { name: d.name, frames: d.data.frames.length, pending: !!d.pending, docs: window.BM.app.state.docs.length } : null; });
  check(ls && ls.name === 'l.prmtop' && ls.frames === 3 && !ls.pending && ls.docs === 1, 'launcher: files named trajectory first still open topology first ' + JSON.stringify(ls));
  await lp.close();
  await open('struct/two_seg.psf');
  s = await st();
  check(s.pending && s.format === 'psf', 'PSF topology pending');
  await shot('03_pending_topology');

  // ---------------------------------------------------------------- 4. convergence plots (Gaussian, ORCA)
  console.log('Convergence plots');
  await open('g16_water_optfreq.log');
  s = await st();
  check(s.plot === 'ir', 'opt+freq keeps the IR plot by default');
  check(await page.isVisible('[data-plot-go]'), 'Convergence shortcut shown next to the spectrum');
  const crit = ['opt:maxForce', 'opt:rmsForce', 'opt:maxDisp', 'opt:rmsDisp', 'opt:energyChange', 'opt:rmsd', 'opt:maxStep', 'opt:all'];
  check(crit.every(c => s.types.includes(c)), 'Gaussian criteria in the picker: ' + s.types.join(' '));
  p = await pick('opt:maxForce');
  check(p.logY && p.hlines.length === 1 && Math.abs(p.hlines[0].y - 0.00045) < 1e-12 && p.points.length === 4, 'max force on a log axis with its threshold');
  p = await pick('opt:rmsd');
  check(p.points.length === 4 && p.points[3].y === 0 && p.points[0].y > 0, 'RMSD to final ends at 0: ' + p.points.map(q => q.y.toFixed(4)));
  p = await pick('opt:all');
  check(p.logY && p.lines.length === 4 && p.hlines[0].y === 1 && p.lines.every(l => l.pts[l.pts.length - 1].y < 1), 'all criteria / threshold with a line at 1');
  await page.evaluate(() => window.BM.app.plot().onPick(1));
  check((await st()).fi === 1, 'clicking a convergence point selects that geometry');
  await page.click('#tabs [data-tab="info"]');
  await page.waitForTimeout(100);
  await shot('01_contents_convergence');
  // job selector
  await page.selectOption('#panel [data-job]', '1');
  s = await st();
  check(!s.types.some(t => t.startsWith('opt:')) && s.fi === 3, 'job 2 (freq) hides optimisation plots ' + JSON.stringify(s.types));
  await page.selectOption('#panel [data-job]', '-1');
  check((await st()).types.includes('opt:all'), 'all jobs again');
  // regression: switching to a job without the shown optimisation criterion must not try to draw it
  await pick('opt:maxForce');
  await page.selectOption('#panel [data-job]', '1');
  s = await st();
  check(s.plot !== 'opt:maxForce' && !/could not be drawn/.test(JSON.stringify((await spec()) || {})), 'job switch replaces a criterion plot the job lacks: ' + s.plot);
  await page.selectOption('#panel [data-job]', '-1');
  // regression: the calculation-input dialog keeps no listeners from an earlier opening (this file has a route)
  await page.evaluate(() => window.BM.app.exportCalc());
  await page.keyboard.press('Escape');
  // Contents entry jumps to the steps tab and the optimisation plot
  await page.click('#panel .ct[data-go-tab="steps"]');
  s = await st();
  check(s.tab === 'steps' && s.plot === 'opt:dE', 'Geometries entry opens Steps + ΔE plot ' + s.tab + ' ' + s.plot);
  await open('orca5_water_optfreq.out');
  await page.evaluate(() => window.BM.app.exportCalc());
  for (const prog of ['gaussian', 'orca', 'psi4', 'gaussian']) await page.selectOption('#dlg [data-c="program"]', prog);
  check(/^%|#/m.test(await page.inputValue('#dlg [data-c="preview"]')), 'calculation dialog of a second file switches programs cleanly');
  await page.keyboard.press('Escape');
  p = await pick('opt:energyChange');
  check(p.hlines.length === 1 && Math.abs(p.hlines[0].y - 5e-6) < 1e-15, 'ORCA energy change threshold');
  p = await pick('opt:rmsd');
  check(p.points[p.points.length - 1].y === 0, 'ORCA RMSD ends at 0');
  for (const id of ['opt:all', 'opt:maxForce']) {
    await pick(id);
    await page.click('[data-act="plot-svg"]'); await page.click('[data-act="plot-csv"]'); await page.click('[data-act="plot-png"]');
  }

  // ---------------------------------------------------------------- 5. scan / IRC, NMR, ECD, orbitals, charges
  console.log('Scan, IRC, NMR, ECD, MO, charges');
  await open('qc/g16_water_modred_scan.log');
  s = await st();
  check(s.plot === 'scan', 'scan file defaults to the scan profile');
  p = await spec();
  check(p.points.length === 2 && p.xLabel === 'A(2,1,3) (°)', 'scan energy vs coordinate');
  await open('qc/g16_hcn_irc.log');
  p = await pick('irc');
  check(p.points.length === 5 && p.points[0].x < 0, 'IRC energy vs reaction coordinate');
  await open('qc/orca5_hcn_irc.out');
  check((await st()).types.includes('irc'), 'ORCA IRC path plot');
  await open('qc/g16_water_nmr.log');
  s = await st();
  check(s.plot === 'nmr', 'NMR file defaults to the NMR plot');
  p = await spec();
  check(p.sticks.length === 1 && p.sticks[0].y === 2 && !p.reverseX, '1H shieldings merged into one stick');
  await page.fill('#plot-controls [data-pc="ref"]', '31.88');
  await page.dispatchEvent('#plot-controls [data-pc="ref"]', 'change');
  p = await spec();
  check(Math.abs(p.sticks[0].x - (31.88 - 31.2341)) < 1e-6 && p.reverseX, 'δ = σref − σ with a reversed axis');
  await page.selectOption('#plot-controls [data-pc="z"]', '8');
  p = await spec();
  check(p.sticks.length === 1 && /17O/.test(p.xLabel), 'nucleus select: 17O');
  check((await contents()).includes('NMR shieldings'), 'Contents lists NMR shieldings');
  await page.click('[data-act="plot-svg"]');
  await open('qc/orca5_water_cd.out');
  p = await pick('ecd');
  const ecdY = p.series[0].y;
  check(p.sticks.length >= 2 && ecdY.some(v => v > 0) && ecdY.some(v => v < 0), 'ECD curve with both signs');
  await page.click('[data-act="plot-csv"]');
  await open('qc/g16_water_ecd.log');
  check((await st()).types.includes('ecd'), 'Gaussian ECD plot available');
  p = await pick('mo');
  check(p.type === 'levels' && p.columns.length >= 1 && p.columns[0].levels.some(l => /HOMO/.test(l.tag)) && p.gaps.length >= 1, 'MO level diagram with HOMO and gap');
  const lumo = p.columns[0].levels.find(l => /LUMO/.test(l.tag));
  await page.evaluate(id => window.BM.app.plot().onPick(id, null), lumo.id);
  check(/LUMO/.test(await page.textContent('#mo-read')), 'clicking a level shows its energy: ' + await page.textContent('#mo-read'));
  await page.click('[data-act="plot-svg"]');
  await open('review/g16_oh_radical_pop.log');
  p = await pick('mo');
  check(p.columns.length === 2, 'open shell: alpha and beta columns');
  p = await pick('charges');
  check(p.bars.length === 2 && p.xTickLabels.length === 2, 'charges bar chart');
  await page.selectOption('#plot-controls [data-pc="scheme"]', 'Hirshfeld');
  check(await page.evaluate(() => window.BM.app.state.active.chargeScheme) === 'Hirshfeld', 'charge scheme select');
  await page.evaluate(() => window.BM.app.plot().onPick(1, null));
  check(await page.evaluate(() => window.BM.app.state.active.selection.join()) === '1', 'clicking a bar selects the atom');
  await page.click('[data-act="plot-svg"]'); await page.click('[data-act="plot-csv"]');
  await page.waitForTimeout(300);
  const plotFiles = await page.evaluate(() => window.__saved.filter(x => /_(opt_all|opt_maxForce|nmr|ecd|mo|charges)\.(svg|csv|png)$/.test(x.name)).map(x => ({ name: x.name, ok: x.text === null ? x.size > 1000 : x.text.length > 40 && (!/\.svg$/.test(x.name) || x.text.startsWith('<?xml')) })));
  check(plotFiles.length >= 11 && plotFiles.every(x => x.ok), 'SVG/PNG/CSV exports of the new plots: ' + plotFiles.map(x => x.name + (x.ok ? '' : ' (bad)')).join(', '));

  // ---------------------------------------------------------------- 6. Save As round trips
  console.log('Save As');
  await open('struct/ligand.mol2');
  const nAt = (await st()).atoms;
  await page.keyboard.press('Control+Shift+S');
  await page.waitForSelector('#dlg[open] [data-sa="fmt"]');
  await page.waitForTimeout(150);
  await shot('02_save_as');
  const disabled = await page.evaluate(() => [...document.querySelectorAll('#dlg [data-sa="fmt"] option:disabled')].map(o => o.value));
  check(disabled.includes('cube'), 'cube disabled without volume data: ' + disabled.join(','));
  await page.click('#dlg .btn:not(.primary)');
  const roundTrip = ['xyz', 'extxyz', 'mol', 'sdf', 'mol2', 'cml', 'cjson', 'pdb', 'cif', 'poscar', 'qe-in', 'xsf', 'gro', 'molden', 'gaussian-input', 'orca-input', 'turbomole-coord'];
  const rt = [];
  for (const id of roundTrip) {
    await page.keyboard.press('Control+Shift+S');
    await page.waitForSelector('#dlg[open] [data-sa="fmt"]');
    await page.selectOption('#dlg [data-sa="fmt"]', id);
    const preview = await page.inputValue('#dlg [data-sa="preview"]');
    const n0 = await page.evaluate(() => window.__saved.length);
    await page.click('#dlg .btn.primary');
    await page.waitForFunction(n => window.__saved.length > n, n0);
    const r = await page.evaluate(() => { const f = window.__saved[window.__saved.length - 1]; try { const d = BM.parseFile(f.name, f.text); return { name: f.name, n: d.frames[0].Z.length }; } catch (e) { return { name: f.name, err: e.message }; } });
    rt.push(r);
    check(preview.length > 20 && r.n === nAt, `Save As ${id} → ${r.name} reopens with ${r.n} atoms (expected ${nAt})${r.err ? ': ' + r.err : ''}`);
  }
  // all geometries
  await open('g16_water_optfreq.log');
  await page.keyboard.press('Control+Shift+S');
  await page.waitForSelector('#dlg[open] [data-sa="fmt"]');
  await page.selectOption('#dlg [data-sa="fmt"]', 'xyz');
  await page.check('#dlg input[name="sa-fr"][value="all"]');
  check(/_all\.xyz$/.test(await page.inputValue('#dlg [data-sa="name"]')), 'all-geometries file name');
  await page.click('#dlg .btn.primary');
  await page.waitForTimeout(200);
  const allFr = await page.evaluate(() => { const f = window.__saved[window.__saved.length - 1]; return BM.parseFile(f.name, f.text).frames.length; });
  check(allFr === 4, 'Save As XYZ with all geometries: ' + allFr + ' frames');

  // ---------------------------------------------------------------- 7. Open Babel (mocked host)
  console.log('Open Babel (mocked)');
  await page.evaluate(() => {
    window.__ob = [];
    const sdf = BM.saveFormats.find(f => f.id === 'sdf').write({ Z: [6, 1, 1, 1, 1], xyz: Float64Array.from([0, 0, 0, 0.63, 0.63, 0.63, -0.63, -0.63, 0.63, -0.63, 0.63, -0.63, 0.63, -0.63, -0.63]), bonds: [[0, 1, 1], [0, 2, 1], [0, 3, 1], [0, 4, 1]], title: 'methane' });
    BM.native.obabel = {
      available: async () => true,
      formats: async () => [{ code: 'pdbqt', description: 'AutoDock PDBQT format', read: true, write: true }, { code: 'mdl', description: 'MDL MOL format', read: true, write: true }, { code: 'smy', description: 'SMILES (read only)', read: true, write: false }],
      convert: async (data, inF, outF) => { window.__ob.push(inF + '>' + outF); return outF === 'sdf' ? sdf : 'CONVERTED ' + outF + '\n' + (typeof data === 'string' ? data : '').split('\n')[0] + '\n'; }
    };
  });
  await open('mystery.cdx', Buffer.from('\u0001\u0002binary-ish chemdraw content\n'));
  s = await st();
  check(s.name === 'mystery.cdx' && s.atoms === 5 && /Open Babel/.test(await toast()), 'unknown format converted with Open Babel: ' + JSON.stringify(s) + ' ' + await toast());
  check(await page.evaluate(() => window.__ob.join()) === 'cdx>sdf', 'Open Babel asked for SDF');
  await page.keyboard.press('Control+Shift+S');
  await page.waitForSelector('#dlg[open] [data-sa="fmt"]');
  await page.waitForFunction(() => document.querySelector('#dlg [data-sa="fmt"] option[value="ob:pdbqt"]'));
  check(!(await page.$('#dlg [data-sa="fmt"] option[value="ob:smy"]')), 'read-only Open Babel formats not offered');
  await page.selectOption('#dlg [data-sa="fmt"]', 'ob:pdbqt');
  check(/\.pdbqt$/.test(await page.inputValue('#dlg [data-sa="name"]')), 'Open Babel file name extension');
  await page.click('#dlg .btn.primary');
  await page.waitForTimeout(200);
  const ob = await page.evaluate(() => ({ calls: window.__ob.join(), last: window.__saved[window.__saved.length - 1] }));
  check(ob.calls === 'cdx>sdf,sdf>pdbqt' && /^CONVERTED pdbqt/.test(ob.last.text) && /\.pdbqt$/.test(ob.last.name), 'Save via Open Babel: SDF converted to pdbqt ' + JSON.stringify(ob.calls));

  // ---------------------------------------------------------------- 8. compact Build tab
  await page.click('#btn-new');
  await page.fill('#bp-smiles', 'CC(=O)Oc1ccccc1C(=O)O');
  await page.click('[data-act="b-smiles"]');
  await page.waitForTimeout(500);
  await page.evaluate(() => { const a = BM.app; a.setTool('select'); a.onAtomClick(1); a.onAtomClick(2); });
  await page.waitForTimeout(150);
  const folded = await page.evaluate(() => [...document.querySelectorAll('.bsec')].map(s => s.dataset.sec + (s.classList.contains('folded') ? '-' : '+')).join(' '));
  check(folded === 'insert+ edit+ align- clean- sym- coords-', 'Build sections: first two open, rest collapsed: ' + folded);
  check(await page.evaluate(() => !document.querySelector('#panel .note')), 'no prose paragraphs in the Build tab');
  await shot('04_build_compact');

  console.log('Console/page errors:', errors.length ? errors : 'none');
  check(errors.length === 0, 'no console errors');
  await browser.close();
  console.log(`\nBrowser format tests: ${passes} passed, ${failures} failed`);
  process.exit(failures ? 1 : 0);
})().catch(e => { console.error(e); process.exit(2); });
