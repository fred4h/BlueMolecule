#!/usr/bin/env node
/* BlueMolecule tests for src/reps.js: the selection language (every keyword, precedence, errors), residue
 * classification (name tables and the topology fallbacks), monosaccharide identification for SNFG (names, CHARMM-GUI
 * truncations, ring topology + stereochemistry for unknown names, sulfation), QM/MM layer keywords on the fixtures in
 * tests/fixtures/qmmm, colour schemes, the representation resolver and presets.
 * Run: node tests/test_reps.js */
'use strict';
const fs = require('fs');
const path = require('path');
const ROOT = path.join(__dirname, '..');
for (const f of ['elements.js', 'geom.js', 'bigfile.js', 'parsers.js', 'readers-struct.js', 'writers.js', 'readers-qc.js', 'chem.js', 'reps.js', 'cartoon.js', 'renderer.js']) require(path.join(ROOT, 'src', f));
const BM = globalThis.BM, R = BM.reps, S = require('./synth.js');

let pass = 0, fail = 0, current = '';
const failures = [];
const ok = (c, msg) => { if (c) pass++; else { fail++; failures.push('[' + current + '] ' + msg); } };
const eq = (a, b, msg) => ok(a === b, msg + ': got ' + JSON.stringify(a) + ', expected ' + JSON.stringify(b));
function test(name, fn) { current = name; try { fn(); } catch (e) { fail++; failures.push('[' + name + '] threw: ' + (e && e.stack || e)); } }

/** Context from a synthetic system (distance bonds like the app). */
function ctxOf(sys, extra) {
  const f = { Z: sys.Z.slice(), xyz: Float64Array.from(sys.xyz), atomNames: sys.names ? sys.names.slice() : null, resNames: sys.res ? sys.res.slice() : null,
    resIds: sys.ids ? sys.ids.slice() : null, chains: sys.chains ? sys.chains.slice() : null };
  if (sys.bfactors) f.bfactors = sys.bfactors;
  return R.context(Object.assign({ frame: f, bonds: BM.geom.computeBonds(f.Z, f.xyz, { tolerance: 0.4 }) }, extra || {}));
}
const idx = (text, ctx) => R.indices(R.evaluate(text, ctx));
const cnt = (text, ctx) => R.count(R.evaluate(text, ctx));
const err = text => { try { R.parse(text); return null; } catch (e) { return e; } };
/** A small hand-made system: ALA-GLY peptide (chain A), HOH, NA, SOD, CLA, a ligand, POPC-named and CHL1-named residues. */
function mini() {
  const sys = { Z: [], xyz: [], names: [], res: [], ids: [], chains: [] };
  const add = (n, z, p, r, id, ch) => { sys.Z.push(z); sys.xyz.push(...p); sys.names.push(n); sys.res.push(r); sys.ids.push(id); sys.chains.push(ch); };
  const bb = S.backbone([-60, -60], [-45, -45]);
  S.residueAtoms(bb[0], 'ALA').forEach(([n, z, p]) => add(n, z, p, 'ALA', 10, 'A'));
  S.residueAtoms(bb[1], 'GLY').forEach(([n, z, p]) => add(n, z, p, 'GLY', 11, 'A'));
  add('H', 1, [bb[1].N[0], bb[1].N[1] + 1.0, bb[1].N[2]], 'GLY', 11, 'A');
  add('OW', 8, [10, 0, 0], 'HOH', 20, 'W'); add('HW1', 1, [10.76, 0.59, 0], 'HOH', 20, 'W'); add('HW2', 1, [9.24, 0.59, 0], 'HOH', 20, 'W');
  add('NA', 11, [14, 0, 0], 'NA', 21, 'I'); add('SOD', 11, [18, 0, 0], 'SOD', 22, 'I'); add('CLA', 17, [22, 0, 0], 'CLA', 23, 'I');
  add('C1', 6, [0, 10, 0], 'LIG', 30, 'L'); add('O1', 8, [1.2, 10, 0], 'LIG', 30, 'L'); add('PT1', 78, [0, 12.1, 0], 'LIG', 30, 'L');
  add('P', 15, [0, 20, 0], 'POPC', 40, 'M'); add('C2', 6, [0, 25, 0], 'CHL1', 41, 'M');
  return sys;
}

// ---------------------------------------------------------------- selection language
test('parse: precedence and parentheses', () => {
  const ast = R.parse('protein or water and not hydrogen');
  eq(ast.op, 'or', 'or is the outermost operator');
  eq(ast.b.op, 'and', 'and binds tighter than or');
  eq(ast.b.b.op, 'not', 'not binds tighter than and');
  eq(R.parse('(protein or water) and heavy').op, 'and', 'parentheses group');
  eq(R.parse('not protein and water').op, 'and', 'not applies to the next term only');
  eq(R.parse('within 5 of ligand and protein').op, 'and', 'within takes one term (use parentheses for more)');
  eq(R.parse('same residue as within 3 of ligand').a.op, 'within', 'same residue as within …');
  eq(R.parse('PROTEIN AND Water').op, 'and', 'keywords are case-insensitive');
});
test('evaluate: keywords on a small system', () => {
  const c = ctxOf(mini());
  const n = c.n;
  eq(cnt('all', c), n, 'all'); eq(cnt('none', c), 0, 'none');
  eq(cnt('protein', c), 10, 'protein = ALA (5 heavy) + GLY (4) + H');
  eq(cnt('backbone', c), 9, 'backbone N CA C O x2 + amide H');
  eq(idx('sidechain', c).map(i => c.frame.atomNames[i]).join(), 'CB', 'sidechain = CB');
  eq(cnt('water', c), 3, 'water'); eq(cnt('ions', c), 3, 'ions NA SOD CLA');
  eq(cnt('lipid', c), 2, 'lipid POPC + CHL1'); eq(cnt('ligand', c), 3, 'ligand LIG');
  eq(cnt('hydrogen', c), 3, 'hydrogen'); eq(cnt('heavy', c), n - 3, 'heavy');
  eq(idx('metal', c).map(i => c.frame.atomNames[i]).join(), 'NA,SOD,PT1', 'metal: Na, Na, Pt');
  eq(cnt('chain A', c), 10, 'chain A'); eq(cnt('chain A W', c), 13, 'chain A W');
  eq(cnt('resname POPC CHL1', c), 2, 'resname list'); eq(cnt('resname po*', c), 1, 'resname wildcard, case-insensitive');
  eq(cnt('resid 10', c), 5, 'resid single'); eq(cnt('resid 10-11', c), 10, 'resid range'); eq(cnt('resid 10 to 11 20', c), 13, 'resid "to" range + single');
  eq(cnt('name CA N C O', c), 8, 'name list'); eq(cnt('name H?1', c), 1, 'name with ?');
  eq(cnt('element C N', c), idx('all', c).filter(i => c.Z[i] === 6 || c.Z[i] === 7).length, 'element');
  eq(cnt('element pt', c), 1, 'element case-insensitive');
  eq(idx('index 0-2', c).join(), '0,1,2', 'index 0-based'); eq(idx('serial 1-3', c).join(), '0,1,2', 'serial 1-based');
  eq(cnt('protein and not hydrogen', c), 9, 'and not'); eq(cnt('water or ions', c), 6, 'or');
  eq(cnt('not (protein or water)', c), n - 13, 'not ( … )');
  eq(cnt('protein or water and hydrogen', c), 12, 'precedence: protein or (water and hydrogen)');
  eq(cnt('(protein or water) and hydrogen', c), 3, 'parentheses override precedence');
  eq(cnt('same residue as name CB', c), 5, 'same residue as');
  eq(cnt('within 1.5 of name OW', c), 3, 'within (includes the source)');
  eq(cnt('sugar', c), 0, 'sugar alias'); eq(cnt('glycan', c), 0, 'glycan alias'); eq(cnt('carbohydrate', c), 0, 'carbohydrate');
  eq(cnt('nucleic', c), 0, 'nucleic');
  eq(cnt('helix or sheet', c), 0, 'too short for secondary structure: no helix or sheet'); eq(cnt('coil', c), 10, 'coil = protein outside helices and sheets');
  void 0;
});
test('evaluate: within and same residue as match a brute-force reference', () => {
  const sys = S.protein(BM), c = ctxOf(sys), n = c.n, xyz = c.xyz;
  const lig = R.evaluate('ligand', c), mask = R.evaluate('within 4.5 of ligand', c);
  let bad = 0;
  for (let i = 0; i < n; i++) {
    let near = false;
    for (let j = 0; j < n && !near; j++) if (lig[j] && Math.hypot(xyz[3 * i] - xyz[3 * j], xyz[3 * i + 1] - xyz[3 * j + 1], xyz[3 * i + 2] - xyz[3 * j + 2]) <= 4.5) near = true;
    if (!!mask[i] !== near) bad++;
  }
  eq(bad, 0, 'within 4.5 of ligand: atoms differing from brute force');
  const site = R.evaluate('same residue as (protein and within 5 of ligand)', c), cls = R.classify(c);
  const hit = new Set();
  for (let i = 0; i < n; i++) if (site[i]) hit.add(cls.resOf[i]);
  ok(hit.size >= 4 && [...hit].every(r => cls.resClass[r] === 1), 'binding site = whole protein residues: ' + hit.size);
  ok(cnt('within 0 of ligand', c) === cnt('ligand', c), 'within 0 = the source');
});
test('evaluate: selected, frozen, active, layer, link (context data)', () => {
  const sys = mini(), n = sys.Z.length;
  const layers = { scheme: 'oniom', atomLayer: sys.Z.map((z, i) => (i < 5 ? 'H' : i < 10 ? 'M' : 'L')), frozen: sys.Z.map((z, i) => i >= 12), active: null, linkAtoms: [{ atom: 5, bondedTo: 2 }], labels: { H: 'B3LYP', M: 'PM6', L: 'UFF' } };
  const c = ctxOf(sys, { layers, selection: [1, 3], frozen: [0] });
  eq(idx('selected', c).join(), '1,3', 'selected'); eq(cnt('layer H', c), 5, 'layer H'); eq(cnt('layer M L', c), n - 5, 'layer M L');
  eq(cnt('qm', c), 5, 'qm = layer H'); eq(cnt('high', c), 5, 'high'); eq(cnt('medium', c), 5, 'medium'); eq(cnt('low', c), n - 10, 'low'); eq(cnt('mm', c), n - 10, 'mm');
  eq(cnt('layer qm2', c), 5, 'layer qm2 alias');
  eq(cnt('frozen', c), n - 12 + 1, 'frozen = layers.frozen + Build-tab frozen');
  eq(cnt('active', c), n - (n - 12 + 1), 'active = not frozen');
  eq(idx('link', c).join(), '5', 'link atoms');
  const c2 = ctxOf(sys);
  eq(cnt('layer H', c2), 0, 'no layers: layer selects nothing (no error)');
  eq(cnt('link', c2), 0, 'no layers: link selects nothing');
});
test('parse errors are clear', () => {
  const cases = [
    ['', /Empty selection/], ['protien', /Unknown keyword "protien" — did you mean "protein"\?/], ['protein and', /"and" needs a selection after it/],
    ['and water', /"and" needs a selection before it/], ['(protein or water', /Missing "\)"/], ['protein)', /Unexpected "\)"/], ['()', /Empty parentheses/],
    ['within of ligand', /"within" needs a distance/], ['within 5 ligand', /Expected "of"/], ['within 5 of', /needs a selection after it/],
    ['same as protein', /Expected "same residue as/], ['resid', /"resid" needs numbers or ranges/], ['resid ten', /"resid" expects numbers or ranges like 10-20, not "ten"/],
    ['layer X', /Layer must be H, M or L/], ['element Xx', /Unknown element "Xx"/], ['protein water', /Unexpected "water": join terms with and \/ or/],
    ['10-20', /put it after a keyword such as resid/], ['name "CA', /Missing closing quote/], ['not', /"not" needs a selection after it/]
  ];
  for (const [t, re] of cases) { const e = err(t); ok(e && e.name === 'SelectionError' && re.test(e.message), JSON.stringify(t) + ' -> ' + (e ? e.message : 'no error')); }
  const e = err('protein and wter');
  ok(e && e.pos === 12 && e.end === 16, 'error position points at the bad word: ' + (e && [e.pos, e.end]));
  const chk = R.check('resname POPC or', ctxOf(mini()));
  ok(!chk.ok && /or/.test(chk.error), 'check() reports errors without throwing');
  ok(R.check('lipid', ctxOf(mini())).count === 2, 'check() counts');
});
test('autocomplete', () => {
  const c = ctxOf(mini());
  let r = R.complete('prot', 4, c);
  ok(r.items.some(i => i.text === 'protein') && r.from === 0 && r.to === 4, 'keyword prefix');
  r = R.complete('resname ', 8, c);
  ok(['ALA', 'GLY', 'HOH', 'POPC', 'CHL1', 'LIG'].every(x => r.items.some(i => i.text === x)), 'resname values: ' + r.items.map(i => i.text));
  r = R.complete('resname PO', 10, c);
  eq(r.items.map(i => i.text).join(), 'POPC', 'resname value prefix');
  r = R.complete('chain ', 6, c);
  ok(r.items.some(i => i.text === 'A') && r.items.some(i => i.text === 'W'), 'chain values');
  r = R.complete('element ', 8, c);
  ok(r.items.some(i => i.text === 'Pt'), 'element values');
  r = R.complete('within ', 7, c);
  ok(r.items.some(i => i.text === '5'), 'distance after within');
  r = R.complete('protein and w', 13, c);
  ok(r.items.some(i => i.text === 'water') && r.items.some(i => i.text === 'within') && r.from === 12, 'keyword after an operator');
  r = R.complete('same ', 5, c);
  eq(r.items.map(i => i.text).join(), 'residue', 'same -> residue');
});

// ---------------------------------------------------------------- classification
test('residue tables', () => {
  const T = R.tables;
  for (const x of ['HID', 'HIE', 'HIP', 'HSD', 'HSE', 'HSP', 'CYX', 'ASH', 'GLH', 'LYN', 'NALA', 'CGLY', 'ACE', 'NME', 'MSE']) ok(T.AMINO.has(x), 'amino variant ' + x);
  for (const x of ['HOH', 'WAT', 'SOL', 'TIP3', 'T3P', 'SPC', 'TIP4']) ok(T.WATER.has(x), 'water ' + x);
  for (const x of ['SOD', 'CLA', 'POT', 'CAL', 'NA', 'CL', 'K', 'MG', 'ZN', 'Na+'.toUpperCase()]) ok(T.IONS.has(x), 'ion ' + x);
  for (const x of ['POPC', 'POPE', 'POPS', 'POPG', 'DPPC', 'DOPC', 'DMPC', 'CHL1', 'PSM', 'PC', 'PE', 'PS', 'PGR', 'PA', 'OL', 'ST', 'CHL']) ok(T.LIPID.has(x), 'lipid ' + x);
  for (const x of ['GLC', 'BGC', 'NAG', 'NDG', 'MAN', 'BMA', 'GAL', 'GLA', 'FUC', 'FUL', 'SIA', 'XYS', 'IDR', 'BDP', 'SGN', 'IDS']) ok(x in T.CCD, 'sugar CCD ' + x);
  for (const x of ['DA', 'DT', 'A', 'U', 'ADE', 'RA5']) ok(T.NUCLEIC.has(x), 'nucleotide ' + x);
});
test('classification by names, CHARMM sugar names and fallbacks', () => {
  const sys = { Z: [], xyz: [], names: [], res: [], ids: [], chains: [] };
  const add = (n, z, p, r, id) => { sys.Z.push(z); sys.xyz.push(...p); sys.names.push(n); sys.res.push(r); sys.ids.push(id); sys.chains.push('X'); };
  // an unknown amino acid (atoms N CA C) and an unknown nucleotide (P C4' C1')
  add('N', 7, [0, 0, 0], 'QZX', 1); add('CA', 6, [1.46, 0, 0], 'QZX', 1); add('C', 6, [2.0, 1.4, 0], 'QZX', 1);
  add('P', 15, [0, 20, 0], 'NUC', 2); add("C4'", 6, [3, 20, 0], 'NUC', 2); add("C1'", 6, [5, 20, 0], 'NUC', 2);
  // lone metal (unknown name) and lone O (water without H)
  add('ZN', 30, [0, -20, 0], 'ZZZ', 3); add('O', 8, [30, 0, 0], 'WWW', 4);
  // CHARMM-GUI style sugar names (no geometry needed for the name rule)
  add('C1', 6, [40, 0, 0], 'BGLCNA', 5); add('C1', 6, [44, 0, 0], 'AIDO', 6); add('C1', 6, [48, 0, 0], 'AGLC', 7);
  const c = ctxOf(sys), cls = R.classify(c), name = r => R.CLASS_IDS[cls.resClass[r]];
  eq([0, 1, 2, 3, 4, 5, 6].map(name).join(), 'protein,nucleic,ions,water,carbohydrate,carbohydrate,carbohydrate', 'fallback classes');
});
test('ring-topology sugar fallback and lipid fallback', () => {
  // glucose from SMILES with an unknown residue name; cyclohexanol-like ring without ring O stays a ligand
  const mk = (smi, resn, shift) => { const m = BM.chem.smilesTo3D(smi, { seed: 2 }).mol; return { Z: m.Z, xyz: Array.from(m.xyz, (v, k) => v + (k % 3 === 0 ? shift : 0)), names: m.Z.map((z, i) => BM.el.symbol(z) + i), res: m.Z.map(() => resn), ids: m.Z.map(() => shift), chains: m.Z.map(() => 'S') }; };
  const join = list => list.reduce((a, b) => ({ Z: a.Z.concat(b.Z), xyz: a.xyz.concat(b.xyz), names: a.names.concat(b.names), res: a.res.concat(b.res), ids: a.ids.concat(b.ids), chains: a.chains.concat(b.chains) }));
  const parts = [mk(S.MONO.Glc.smi, 'QQ1', 0), mk('OC1CCCCC1', 'QQ2', 20), mk('C1CCOCC1', 'QQ3', 40)];
  // three copies of a long-chain amphiphile with an unknown name, and one single copy
  for (let k = 0; k < 3; k++) parts.push(mk('CCCCCCCCCCCCCCCC(=O)OCC(O)CO', 'MYC', 60 + 25 * k));
  parts.push(mk('CCCCCCCCCCCCCCCC(=O)OCC(O)CO', 'ONE', 160));
  const c = ctxOf(join(parts)), cls = R.classify(c), cl = nm => R.CLASS_IDS[cls.resClass[cls.resName.indexOf(nm)]];
  eq(cl('QQ1'), 'carbohydrate', 'glucose with an unknown name -> carbohydrate (ring O + exocyclic O)');
  eq(cl('QQ2'), 'ligand', 'cyclohexanol (no ring O) -> ligand');
  eq(cl('QQ3'), 'ligand', 'tetrahydropyran without exocyclic O -> ligand');
  eq(cl('MYC'), 'lipid', 'long-chain residue occurring 3 times -> lipid');
  eq(cl('ONE'), 'ligand', 'a single long-chain molecule stays a ligand');
});
test('no residue names: molecules are the residues', () => {
  const m = BM.chem.smilesTo3D('O.O.CCO', { seed: 1 }).mol;
  const c = R.context({ frame: { Z: m.Z, xyz: m.xyz }, bonds: BM.geom.computeBonds(m.Z, m.xyz, { tolerance: 0.4 }) });
  const cls = R.classify(c);
  eq(cls.nRes, 3, 'three molecules'); eq(R.count(R.evaluate('water', c)), 6, 'two waters by topology');
  eq(R.count(R.evaluate('same residue as index 0', c)), 3, 'same residue = same molecule');
});

// ---------------------------------------------------------------- monosaccharides
test('SNFG identification from ring topology and stereochemistry', () => {
  const want = { Glc: ['Glc', 'sphere', 'D'], Gal: ['Gal', 'sphere', 'D'], Man: ['Man', 'sphere', 'D'], Fuc: ['Fuc', 'cone', 'L'], Xyl: ['Xyl', 'star', 'D'], GlcA: ['GlcA', 'ddiamond', 'D'],
    IdoA: ['IdoA', 'ddiamond', 'L'], GlcNAc: ['GlcNAc', 'cube', 'D'], GlcN: ['GlcN', 'dcube', 'D'], GlcNS6S: ['GlcNS6S', 'dcube', 'D'], IdoA2S: ['IdoA2S', 'ddiamond', 'L'], Neu5Ac: ['Neu5Ac', 'diamond', ''] };
  for (const [k, [label, shape, dl]] of Object.entries(want)) {
    const m = BM.chem.smilesTo3D(S.MONO[k].smi, { seed: 5 }).mol;
    const c = R.context({ frame: { Z: m.Z, xyz: m.xyz, resNames: m.Z.map(() => 'UNK'), resIds: m.Z.map(() => 1), atomNames: m.Z.map((z, i) => 'X' + i) }, bonds: BM.geom.computeBonds(m.Z, m.xyz, { tolerance: 0.4 }) });
    const s = R.sugars(c);
    ok(s.length === 1 && s[0].label === label && s[0].shape === shape && (!dl || s[0].dl === dl), k + ' -> ' + (s[0] ? [s[0].label, s[0].shape, s[0].dl].join(' ') : 'none'));
  }
  const col = (k, c) => { const m = BM.chem.smilesTo3D(S.MONO[k].smi, { seed: 5 }).mol; const cc = R.context({ frame: { Z: m.Z, xyz: m.xyz }, bonds: BM.geom.computeBonds(m.Z, m.xyz, { tolerance: 0.4 }) }); return R.toHex(R.sugars(cc)[0].color) === c; };
  ok(col('Glc', '#0090bc') && col('Gal', '#ffd400') && col('Man', '#00a651') && col('Fuc', '#ed1c24') && col('Xyl', '#f47920') && col('Neu5Ac', '#a54399') && col('IdoA', '#a17a4d'), 'SNFG colours (blue, yellow, green, red, orange, purple, brown)');
});
test('heparan sulfate chain: names, truncated CHARMM names and unknown names', () => {
  const sys = S.glycan(BM), c = ctxOf(sys), cls = R.classify(c), s = R.sugars(c);
  eq(s.length, sys.expected.length, 'one ring per residue');
  for (const e of sys.expected) {
    const r = cls.resName.findIndex((nm, k) => nm === e.resn && cls.resId[k] === e.rid), hit = s.find(x => x.res === r);
    ok(hit && hit.label === e.label, e.resn + ' ' + e.rid + ' -> ' + (hit && hit.label) + ' (expected ' + e.label + ')');
  }
  eq(R.count(R.evaluate('carbohydrate', c)), sys.res.filter((r, i) => sys.chains[i] === 'G').length, 'carbohydrate covers every glycan atom (unknown names via ring topology)');
  const i = sys.names.findIndex((n, k) => n === 'C1' && sys.res[k] === 'HS2');
  ok(/HS2 510 G · IdoA2S/.test(R.describe(c, i)), 'describe() for the status bar: ' + R.describe(c, i));
  eq(R.autoPreset(c), 'glycan', 'glycan preset for a proteoglycan');
});

// ---------------------------------------------------------------- QM/MM fixtures
test('layer keywords and presets on the QM/MM fixtures', () => {
  const DIR = path.join(__dirname, 'fixtures', 'qmmm'), EXP = JSON.parse(fs.readFileSync(path.join(DIR, 'expected.json'), 'utf8'));
  let n = 0;
  for (const file of fs.readdirSync(DIR).filter(f => !f.endsWith('.json'))) {
    const ds = BM.parseFile(file, fs.readFileSync(path.join(DIR, file), 'utf8')), f = ds.frames[ds.frames.length - 1], L = ds.layers;
    const c = R.context({ frame: f, bonds: BM.geom.computeBonds(f.Z, f.xyz, { tolerance: 0.4 }), layers: L && L.atomLayer && L.atomLayer.length === f.Z.length ? L : null });
    if (!c.layers) continue;
    n++;
    for (const k of ['H', 'M', 'L']) eq(R.count(R.evaluate('layer ' + k, c)), L.atomLayer.filter(x => x === k).length, file + ' layer ' + k);
    eq(R.count(R.evaluate('link', c)), new Set(L.linkAtoms.map(a => a.atom)).size, file + ' link atoms');
    if (L.frozen) eq(R.count(R.evaluate('frozen', c)), L.frozen.filter(Boolean).length, file + ' frozen');
    eq(R.autoPreset(c), 'qmmm', file + ' auto preset');
    const rows = R.presetRows('qmmm', c, {}), res = R.resolve(rows, c);
    ok(rows[0].selection === 'layer H' && rows[0].style === 'ballstick', file + ' H as ball-and-stick');
    if (L.atomLayer.includes('M')) ok(rows.some(r => r.selection === 'layer M' && r.style === 'sticks'), file + ' M as sticks');
    if (L.atomLayer.includes('L')) ok(rows.some(r => r.selection === 'layer L' && r.style === 'lines'), file + ' L as lines');
    if (L.frozen && L.frozen.some(Boolean)) ok(rows.some(r => /^frozen/.test(r.selection) && (r.opacity < 1 || r.style === 'ghost')), file + ' frozen atoms dimmed or ghost');
    if (L.linkAtoms.length) ok(res.markers.length === 1 && res.markers[0].indices.length === new Set(L.linkAtoms.map(a => a.atom)).size, file + ' link markers');
    const H = L.atomLayer.indexOf('H');
    eq(res.style[H], R.CODE.ballstick, file + ' per-atom style of a high-layer atom');
    void EXP;
  }
  ok(n >= 5, 'fixtures with layers: ' + n);
});

// ---------------------------------------------------------------- colours, resolver, presets
test('colour schemes', () => {
  const sys = mini(), layers = { scheme: 'oniom', atomLayer: sys.Z.map((z, i) => (i < 5 ? 'H' : 'L')), frozen: null, active: null, linkAtoms: [], labels: {} };
  sys.bfactors = Float32Array.from(sys.Z, (z, i) => i);
  const c = ctxOf(sys, { layers, charges: sys.Z.map((z, i) => (i % 2 ? 0.5 : -0.5)) });
  const at = (cols, i) => R.toHex([cols[3 * i], cols[3 * i + 1], cols[3 * i + 2]]);
  eq(at(R.colorize('layer', c), 0), R.toHex(R.LAYER_COL.H), 'layer colour H'); eq(at(R.colorize('layer', c), 7), R.toHex(R.LAYER_COL.L), 'layer colour L');
  eq(at(R.colorize('#ff8800', c), 3), '#ff8800', 'uniform colour');
  eq(at(R.colorize('element', c), 1), R.toHex(BM.el.color[6]), 'element colour of CA');
  ok(at(R.colorize('chain', c), 0) !== at(R.colorize('chain', c), 11), 'chains differ');
  ok(at(R.colorize('ss', c), 1) === R.toHex(R.SS_COL.coil), 'coil colour for a short peptide');
  const q = R.colorize('charge', c);
  ok(q[2] > q[0] && q[3] > q[5], 'charge: negative blue, positive red');
  const b = R.colorize('bfactor', c);
  ok(b[2] > b[0] && b[3 * (sys.Z.length - 1)] > b[3 * (sys.Z.length - 1) + 2], 'B-factor ramp blue -> red');
  const r = R.colorize('restype', c);
  ok(at(r, 10) === R.toHex([0.4980392156862745, 0.8274509803921568, 0.9411764705882353]) || at(r, 10) !== at(r, 0), 'residue type colours');
});
test('resolver: later rows override, overlays add, extent, markers, errors', () => {
  const sys = S.protein(BM), c = ctxOf(sys);
  const rows = [R.row('all', 'all', 'lines'), R.row('protein', 'protein', 'cartoon', 'ss'), R.row('lig', 'ligand', 'ballstick'), R.row('bad', 'protien', 'sticks'),
    R.row('water', 'water', 'hidden'), R.row('ca', 'name CA', 'markers', '#00ff00'), R.row('off', 'ions', 'spacefill', 'element', { visible: false }), R.row('ghost', 'resid 1-3', 'ghost', 'element', { opacity: 0.4 })];
  const res = R.resolve(rows, c), cl = R.evaluate('ligand', c), wa = R.evaluate('water', c), io = R.evaluate('ions', c), first = R.evaluate('resid 1-3 and chain A', c);
  const li = R.indices(cl)[0], wi = R.indices(wa)[0], ii = R.indices(io)[0], fi = R.indices(first)[0];
  eq(res.style[li], R.CODE.ballstick, 'ligand row overrides the lines row');
  eq(res.style[wi], R.CODE.hidden, 'water hidden overrides lines');
  eq(res.style[ii], R.CODE.lines, 'a hidden (invisible) row changes nothing');
  eq(res.style[fi], R.CODE.ghost, 'ghost row'); ok(Math.abs(res.alpha[fi] - 0.4) < 1e-6, 'ghost opacity');
  eq(res.overlays.length, 1, 'cartoon overlay'); eq(res.overlays[0].kind, 'cartoon', 'overlay kind');
  ok(res.errors[3] && /did you mean "protein"/.test(res.errors[3]), 'error kept per row');
  eq(res.counts[2], R.count(cl), 'per-row atom counts');
  eq(res.markers.length, 1, 'markers row'); eq(R.toHex(res.markers[0].color), '#00ff00', 'marker colour');
  eq(res.extent[wi], 0, 'hidden water not in the fitting extent');
  const prot = R.indices(R.evaluate('protein', c));
  ok(prot.every(i => res.extent[i] === 1), 'cartoon atoms count for fitting');
  const lip = R.resolve([R.row('l', 'all', 'lipid')], c);
  eq(lip.style[R.indices(R.evaluate('hydrogen', c))[0]], R.CODE.hidden, 'lipid style hides hydrogens');
  ok(R.isDefault([R.row('All atoms', 'all', 'sticks', 'element', { base: true })]) && !R.isDefault([R.row('All atoms', 'all', 'lines', 'element', { base: true })]) &&
    !R.isDefault([R.row('All atoms', 'all', 'sticks', 'element', { base: true, opacity: 0.5 })]), 'isDefault: only the plain toolbar-style row');
  eq(JSON.stringify(R.CODE), JSON.stringify(BM.Renderer.STYLE), 'reps style codes = renderer style codes');
});
test('presets for each kind of system', () => {
  const prot = ctxOf(S.protein(BM)), mem = ctxOf(S.membrane({ n: 6, protein: true })), gly = ctxOf(S.glycan(BM));
  eq(R.autoPreset(prot), 'protein', 'protein auto preset'); eq(R.autoPreset(mem), 'membrane', 'membrane auto preset'); eq(R.autoPreset(gly), 'glycan', 'glycan auto preset');
  const small = R.context({ frame: { Z: [8, 1, 1], xyz: Float64Array.from([0, 0, 0, 0.96, 0, 0, -0.24, 0.93, 0]) } });
  eq(R.autoPreset(small), 'all', 'small molecule keeps the default');
  const p = R.presetRows('protein', prot, {});
  ok(p[0].style === 'cartoon' && p[0].color === 'ss' && p.some(r => r.style === 'sticks' && /within 5 of ligand/.test(r.selection)) && p.some(r => r.selection === 'ligand' && r.style === 'ballstick') &&
    p.some(r => r.selection === 'water' && r.style === 'hidden') && p.some(r => r.selection === 'ions' && r.style === 'spacefill' && r.scale < 1), 'protein-ligand rows: ' + p.map(r => r.selection + ':' + r.style).join(' | '));
  const m = R.presetRows('membrane', mem, {});
  ok(m.some(r => r.selection === 'protein' && r.style === 'cartoon') && m.some(r => r.selection === 'lipid' && r.style === 'lipid') && m.some(r => r.selection === 'water' && r.style === 'hidden'), 'membrane rows');
  const mr = R.resolve(m, mem), P = R.indices(R.evaluate('lipid and element P', mem))[0], C = R.indices(R.evaluate('lipid and name C216', mem))[0];
  ok(mr.style[P] === R.CODE.spacefill && mr.style[C] === R.CODE.lines, 'lipid style: P sphere, tail carbon line');
  const g = R.presetRows('glycan', gly, {});
  ok(['sticks', 'rings', 'snfg'].every(st => g.some(r => r.selection === 'carbohydrate' && r.style === st)) && g.some(r => r.style === 'cartoon'), 'glycan rows');
  const a = R.presetRows('all', null, { style: 'spacefill', color: 'charge' });
  ok(a.length === 1 && a[0].base && a[0].style === 'spacefill' && a[0].color === 'charge', 'all-atoms row follows the toolbar');
  const cont = R.contents(mem);
  ok(cont.classes.lipid.res > 20 && cont.classes.protein.res === 144 && cont.classes.water.res > 100 && cont.classes.ions.res > 0, 'contents counts: ' + JSON.stringify(Object.fromEntries(Object.entries(cont.classes).map(([k, v]) => [k, v.res]))));
});
test('large system timings (230k atoms)', () => {
  const sys = S.membrane({ n: 26, protein: true }), f = { Z: sys.Z, xyz: Float64Array.from(sys.xyz), atomNames: sys.names, resNames: sys.res, resIds: sys.ids, chains: sys.chains };
  const b = BM.geom.computeBonds(f.Z, f.xyz, { tolerance: 0.4 }), c = R.context({ frame: f, bonds: b });
  let t = Date.now();
  R.classify(c);
  const tc = Date.now() - t;
  t = Date.now();
  const res = R.resolve(R.presetRows('membrane', c, {}), c);
  const tr = Date.now() - t;
  t = Date.now();
  const w = R.count(R.evaluate('same residue as (lipid and within 8 of protein)', c));
  const tw = Date.now() - t;
  console.log(`  ${f.Z.length} atoms: classify ${tc} ms, membrane preset ${tr} ms, within-query ${tw} ms (${w} atoms)`);
  ok(tc < 1500 && tr < 1500 && tw < 1500 && res.counts.every(x => x > 0), 'classification, resolution and within stay fast');
});

console.log(`\nBlueMolecule representation tests: ${pass} passed, ${fail} failed`);
for (const f of failures) console.log('  FAIL ' + f);
process.exit(fail ? 1 : 0);
