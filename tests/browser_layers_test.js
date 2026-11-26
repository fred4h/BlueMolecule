// Browser checks of 2.6.1 and 2.6.2 in dist/BlueMolecule.html (headless Chromium): QM/MM layers that an output does not list.
//   - a Gaussian ONIOM job that took its molecule from the checkpoint (Geom=AllCheck): layers rebuilt from the link
//     atoms and shown with the QM/MM display; when they cannot be rebuilt, the Info tab says so
//   - Info tab "ONIOM layers": take the layers (and residue names) from an open document with the same atoms
//   - a GaussView-style ONIOM input with residue information opens with its layers
//   - ORCA QM/MM without the QM atom list: layers unknown (not every atom MM), then taken from another document
//   - 2.6.2: layers read from the population analysis of the high-level calculation; "Open file with the layers…";
//     the Contents item of a file without layers
// Usage: node tests/browser_layers_test.js
const { chromium } = require('playwright');
const path = require('path');
const fs = require('fs');

const ROOT = path.resolve(__dirname, '..');
const HTML = 'file://' + path.join(ROOT, 'dist', 'BlueMolecule.html');
const FIX = path.join(ROOT, 'tests', 'fixtures');
const rd = f => fs.readFileSync(path.join(FIX, f), 'utf8');

let failures = 0, passes = 0;
const check = (cond, msg) => { if (cond) passes++; else { failures++; console.log('  FAIL:', msg); } };

/** The Amber ONIOM fixture as a Geom=AllCheck job (see tests/test_qmmm.js): checkpoint rows without layers, NQM. */
function allcheck(text, nqm) {
  const head = text.indexOf(' Symbolic Z-matrix:'), end = text.indexOf('\n \n', head);
  const rows = text.slice(head, end).split('\n').filter(l => /^ [A-Z][a-z]?-/.test(l)).map(l => l.trim().split(/\s+/));
  const block = ' Structure from the checkpoint file:  "job0.chk"\n' +
    ' Charge =  0 Multiplicity = 1 for low   level calculation on real  system.\n' +
    ' Charge =  0 Multiplicity = 1 for high  level calculation on model system.\n' +
    ' Charge =  0 Multiplicity = 1 for low   level calculation on model system.\n' +
    rows.map(t => ' ' + t[0].split('-')[0] + ',' + t[1] + ',' + t[2] + ',' + t[3] + ',' + t[4]).join('\n') + '\n Recover connectivity data from disk.';
  return (text.slice(0, head) + block + text.slice(end)).replace(/ %chk=job\.chk\n/, ' %oldchk=job0.chk\n %chk=job.chk\n')
    .replace(/geom=connectivity/, 'geom=allcheck')
    .replace(/NAtoms=\s+17 NQM=\s+\d+ NQMF=\s+\d+ NMMI=\s+\d+ NMMIF=\s+\d+/, 'NAtoms=     17 NQM=        ' + String(nqm).padStart(2) + ' NQMF=       0 NMMI=      9 NMMIF=      3');
}

/** The high-level model calculation's Mulliken table (every atom: QM atoms and the link atom, printed as H, with charges;
 *  the others 0) added to the first step (see withPop in tests/test_qmmm.js). */
function withPop(text, layers) {
  const Z = [6, 1, 1, 1, 6, 1, 1, 6, 8, 8, 1, 8, 1, 1, 8, 1, 1], sym = ['', 'H', '', '', '', '', 'C', '', 'O'];
  const rows = Z.map((z, a) => String(a + 1).padStart(6) + '  ' + (a === 0 ? 'H' : sym[z]).padEnd(2) + (layers[a] === 'H' || a === 0 ? 0.05 + 0.01 * a : 0).toFixed(6).padStart(11));
  const k = text.indexOf(' ONIOM: calculating energy.');
  return text.slice(0, k) + ' ONIOM: generating point  2 -- high level on model system.\n NAtoms=   17 NActive=    8 NUniq=    8 SFac= 1.00D+00 NAtFMM=   60 NAOKFM=F Big=F\n' +
    ' Mulliken charges:\n               1\n' + rows.join('\n') + '\n Sum of Mulliken charges =   0.00000\n' + text.slice(k);
}

(async () => {
  const browser = await chromium.launch({ args: ['--enable-unsafe-swiftshader'] });
  const page = await (await browser.newContext({ viewport: { width: 1440, height: 950 } })).newPage();
  const errors = [];
  page.on('pageerror', e => errors.push('pageerror: ' + e.message));
  page.on('console', m => { if (m.type() === 'error' && !/Failed to load resource/.test(m.text())) errors.push('console: ' + m.text()); });
  await page.goto(HTML);
  await page.waitForFunction(() => window.BM && window.BM.app);
  const open = (name, text) => page.evaluate(([n, t]) => { BM.app.openText(n, t); }, [name, text]);
  const act = name => page.evaluate(n => BM.app.activate(BM.app.state.docs.find(d => d.name === n)), name);
  const info = async () => { await page.click('#tabs button[data-tab="info"]'); return page.$eval('#sec-layers', e => e.innerText).catch(() => null); };
  const state = () => page.evaluate(() => { const d = BM.app.state.active, L = d.data.layers;
    return { layers: L && L.atomLayer ? L.atomLayer.join('') : null, source: L ? L.source : null, from: L ? L.from || null : null, preset: d.repPreset,
      rows: d.reps ? d.reps.map(r => r.selection) : [], notes: d.data.warnings.slice(), res: BM.app.curFrame(d).resNames || null }; });
  const amber = rd('qmmm/g16_oniom_amber_opt.log');
  const want = 'LLLLHHHHHHHLLLLLL';

  // ------------------------------------------------------------ rebuilt automatically
  await open('restart.log', allcheck(amber, 8));
  let s = await state();
  check(s.layers === want && s.source === 'rebuilt', 'Geom=AllCheck job: layers rebuilt from the link atom: ' + JSON.stringify(s));
  check(s.preset === 'qmmm' && s.rows.includes('layer H') && s.rows.includes('layer L'), 'the QM/MM display is used: ' + s.rows.join(', '));
  const lg = await page.evaluate(() => { const r = BM.app.state.active._repRes; return r && r.legend ? r.legend.title + ': ' + r.legend.items.map(i => i.label + ' ' + i.sub).join('; ') : null; });
  check(/^ONIOM layers: High B3LYP\/6-31G\(d\) · 7 atoms; Low AMBER · 10 atoms; Frozen 3 atoms; Link atoms 1$/.test(lg || ''), 'legend: ' + lg);
  let sec = await info();
  check(sec && /ONIOM LAYERS/i.test(sec) && /Rebuilt from the link atoms/.test(sec), 'Info: rebuilt: ' + sec);
  check(/started the series \(this job continued from job0\.chk\)/.test(sec), 'Info: which file to open: ' + sec);

  // ------------------------------------------------------------ not rebuildable: unknown, then taken from another document
  await open('restart-water.log', allcheck(amber, 11));
  s = await state();
  check(s.layers === null && s.preset !== 'qmmm', 'water in the QM region: layers unknown, no QM/MM display: ' + JSON.stringify(s));
  check(s.notes.some(n => /does not list which atoms are in which ONIOM layer/.test(n)), 'reader note');
  sec = await info();
  check(/Not listed in this file/.test(sec) && /From\s+restart\.log/.test(sec), 'Info: not listed; the rebuilt document (same atoms) is offered: ' + sec);
  check(await page.$eval('#panel .contents', e => { const b = [...e.querySelectorAll('.ct.warn')].find(x => /ONIOM layers/.test(x.textContent)); return !!b && b.dataset.goAnchor === 'sec-layers' && /not listed/.test(b.textContent); }),
    'Contents: an "ONIOM layers: not listed" item that leads to the section');
  await open('water.log', rd('g16_water_optfreq.log'));        // other atoms: never offered
  await open('job0.log', amber);                               // the job that started the series
  await act('restart-water.log');
  await info();
  const offered = await page.$$eval('[data-take-layers] option', o => o.map(x => x.textContent));
  check(offered.join() === 'restart.log,job0.log', 'documents with the same atoms are offered (not water.log): ' + offered.join());
  await page.selectOption('[data-take-layers]', { label: 'job0.log' });
  await page.click('[data-take-layers-go]');
  s = await state();
  check(s.layers === want && s.source === 'document' && s.from === 'job0.log', 'layers taken: ' + JSON.stringify(s));
  check(s.preset === 'qmmm' && s.rows.includes('link'), 'QM/MM display after taking them: ' + s.rows.join(', '));
  check(!s.notes.some(n => /take its layers/.test(n)), 'the reader note about missing layers is gone');
  check(/Layers taken from job0\.log: 7 high, 10 low, 1 link atom/.test(await page.$eval('#status', e => e.textContent)), 'status message');
  check(/Taken from job0\.log/.test(await info()), 'Info: taken from');
  check(await page.evaluate(() => { const d = BM.app.state.active; return d.data.layers.frozen && d.data.layers.frozen.filter(Boolean).length === 3; }), 'frozen atoms kept from this job');

  // ------------------------------------------------------------ read from the population analysis of the high-level calculation
  await open('restart-pop.log', withPop(allcheck(amber, 11), want));
  s = await state();
  check(s.layers === want && s.source === 'population' && s.preset === 'qmmm', 'population analysis: layers and the QM/MM display: ' + JSON.stringify(s));
  sec = await info();
  check(/Read from the population analysis of the high-level calculation/.test(sec) && /Open file with the layers/.test(sec), 'Info: read from the population analysis: ' + sec);

  // ------------------------------------------------------------ "Open file with the layers…": the file opened gives its layers
  await open('restart-open.log', allcheck(amber, 11));
  s = await state();
  check(s.layers === null, 'another job without layers: unknown');
  await info();
  let chooser = page.waitForEvent('filechooser');
  await page.click('[data-act="layers-open"]');
  await (await chooser).setFiles({ name: 'water.log', mimeType: 'text/plain', buffer: Buffer.from(rd('g16_water_optfreq.log')) });
  await page.waitForFunction(() => BM.app.state.docs.filter(d => d.name === 'water.log').length === 2);
  check(/water\.log does not list QM\/MM layers either: nothing was copied to restart-open\.log/.test(await page.$eval('#status', e => e.textContent)), 'a file without layers: nothing copied, and the status says so');
  await act('restart-open.log');
  check((await state()).layers === null, 'still unknown');
  await info();
  chooser = page.waitForEvent('filechooser');
  await page.click('[data-act="layers-open"]');
  await (await chooser).setFiles({ name: 'first-job.log', mimeType: 'text/plain', buffer: Buffer.from(amber) });
  await page.waitForFunction(() => BM.app.state.active && BM.app.state.active.name === 'restart-open.log' && BM.app.state.active.data.layers.source === 'document');
  s = await state();
  check(s.layers === want && s.from === 'first-job.log' && s.preset === 'qmmm', 'the file opened gives its layers to the job: ' + JSON.stringify(s));
  check(await page.evaluate(() => BM.app.state.docs.some(d => d.name === 'first-job.log')), 'the file stays open as its own document');

  // ------------------------------------------------------------ a GaussView-style input with residue information
  const gjf = ['%chk=test.chk', '#p oniom(b3lyp/6-31g(d):amber)=embedcharge opt geom=connectivity', '', 'GaussView-style ONIOM input', '', '0 1 0 1 0 1',
    ' N-N3-0.1592(PDBName=N,ResName=SER,ResNum=1)     0   -1.200000    0.300000    0.000000 L',
    ' H-H-0.1849(PDBName=H1,ResName=SER,ResNum=1)     0   -1.900000    1.000000    0.000000 L',
    ' C-CX-0.0567(PDBName=CA,ResName=SER,ResNum=1)    0    0.000000    0.000000    0.000000 L H-HC 4',
    ' C-2C-0.2596(PDBName=CB,ResName=SER,ResNum=1)    0    0.800000    1.300000    0.000000 H',
    ' O-OH--0.6714(PDBName=OG,ResName=SER,ResNum=1)  -1    2.100000    1.100000    0.000000 H',
    ' H-HO-0.4239(PDBName=HG,ResName=SER,ResNum=1)    0    2.500000    1.950000    0.000000 H', '', ''].join('\n');
  await open('serine.gjf', gjf);
  s = await state();
  check(s.layers === 'LLLHHH' && s.preset === 'qmmm' && s.res && s.res[0] === 'SER', 'GaussView input: layers, QM/MM display, residues: ' + JSON.stringify(s));
  check((await info()) === null, 'Info: no layers section when the file lists them');

  // ------------------------------------------------------------ ORCA without the QM atom list
  const orca = rd('qmmm/orca6_qmmm_opt.out');
  const bare = orca.replace(/^\|\s*\d+>\s*QMAtoms \{4:10\} end\n/m, '')
    .replace(/-{33}\nCARTESIAN COORDINATES \(ANGSTROEM\)\n-{33}\n((?:  [A-Z][a-z]?\s+\S+\s+\S+\s+\S+\n){8})\n/g, '').replace(/^Link atom.*\n/gm, '');
  await open('qmmm-bare.out', bare);
  s = await state();
  check(s.layers === null && s.preset !== 'qmmm', 'ORCA without QM atoms: unknown, not all MM: ' + JSON.stringify(s));
  sec = await info();
  check(/QM\/MM LAYERS/i.test(sec) && /Not listed in this file/.test(sec), 'Info: QM/MM layers not listed: ' + sec);
  await open('qmmm.out', orca);
  await act('qmmm-bare.out');
  await info();
  await page.click('[data-take-layers-go]');
  s = await state();
  check(s.layers === 'LLLLHHHHHHHLLLLLL' && s.preset === 'qmmm', 'ORCA layers taken from the full output: ' + JSON.stringify(s));

  check(!errors.length, 'no page errors: ' + errors.join(' | '));
  await browser.close();
  console.log(`Layers browser tests: ${passes} passed, ${failures} failed`);
  process.exit(failures ? 1 : 0);
})().catch(e => { console.error(e); process.exit(1); });
