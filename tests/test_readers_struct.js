#!/usr/bin/env node
/* BlueMolecule structure-reader tests (src/readers-struct.js).  Run:  node tests/test_readers_struct.js
 * Text fixtures live in tests/fixtures/struct/. DCD, TRR and XTC files are generated here with small
 * independent byte writers (the XTC encoder packs bits with BigInt arithmetic, unlike the reader). */
'use strict';
const fs = require('fs');
const path = require('path');
const vm = require('vm');
const ROOT = path.join(__dirname, '..');
for (const f of ['elements', 'parsers', 'writers', 'readers-struct']) require(path.join(ROOT, 'src', f + '.js'));
const BM = globalThis.BM;
const FIX = path.join(__dirname, 'fixtures', 'struct');
const read = f => fs.readFileSync(path.join(FIX, f), 'utf8');
const load = f => BM.parseFile(f, read(f));
const BOHR = 0.529177210903, HARTREE_EV = 27.211386245988, EVA_AU = BOHR / HARTREE_EV, D2R = Math.PI / 180;

let pass = 0, fail = 0, current = '';
const failures = [];
function ok(c, msg) { if (c) pass++; else { fail++; failures.push('[' + current + '] ' + msg); } }
function eq(a, b, msg) { ok(a === b, msg + ': got ' + JSON.stringify(a) + ', expected ' + JSON.stringify(b)); }
function near(a, b, tol, msg) { ok(typeof a === 'number' && Math.abs(a - b) <= tol, msg + ': got ' + a + ', expected ' + b); }
function arrNear(a, b, tol, msg) {
  if (!a || a.length !== b.length) { ok(false, msg + ': length ' + (a && a.length) + ' vs ' + b.length); return; }
  for (let i = 0; i < b.length; i++) if (!(Math.abs(a[i] - b[i]) <= tol)) { ok(false, msg + '[' + i + ']: got ' + a[i] + ', expected ' + b[i]); return; }
  pass++;
}
function arrEq(a, b, msg) { ok(JSON.stringify(Array.from(a || [])) === JSON.stringify(b), msg + ': got ' + JSON.stringify(a && Array.from(a)) + ', expected ' + JSON.stringify(b)); }
const flat = c => c.reduce((a, v) => a.concat(v), []);
function test(name, fn) { current = name; try { fn(); } catch (e) { fail++; failures.push('[' + name + '] threw: ' + (e && e.stack || e)); } }
const atom = (f, i) => [f.xyz[3 * i], f.xyz[3 * i + 1], f.xyz[3 * i + 2]];
const dist = (f, i, j) => Math.hypot(f.xyz[3 * i] - f.xyz[3 * j], f.xyz[3 * i + 1] - f.xyz[3 * j + 1], f.xyz[3 * i + 2] - f.xyz[3 * j + 2]);
const vadd = (...v) => v.reduce((s, x) => s.map((c, k) => c + x[k]), [0, 0, 0]);
const vs = (v, s) => v.map(x => x * s);

// ------------------------------------------------------------------ detection
test('format detection', () => {
  const want = { 'nacl_fm3m.cif': 'cif', 'p1bar_labels.cif': 'cif', 'small.mmcif': 'mmcif', 'ligand.mol2': 'mol2', 'gaff2.mol2': 'mol2',
    'attr.cml': 'cml', 'array.cml': 'cml', 'water.cjson': 'cjson', 'nacl_frac.cjson': 'cjson', POSCAR: 'poscar', CONTCAR: 'poscar',
    'si_v4.vasp': 'poscar', XDATCAR: 'xdatcar', OUTCAR: 'outcar', 'si_ibrav0.pwi': 'qe-in', 'si_ibrav2.in': 'qe-in', 'h2_bohr.pwi': 'qe-in',
    'relax.pwo': 'qe-out', 'si.cell': 'castep-cell', 'h2o_abs.cell': 'castep-cell', 'grid.xsf': 'xsf', 'anim.axsf': 'xsf',
    'h2o.inp': 'cp2k-in', 'water.prmtop': 'prmtop', 'water.inpcrd': 'inpcrd', 'two_seg.psf': 'psf', 'full.data': 'lammps-data',
    'atomic.data': 'lammps-data', 'traj.lammpstrj': 'lammps-dump' };
  for (const [f, id] of Object.entries(want)) eq(BM.detectFormat(f, read(f)), id, f);
  // content-only detection (file renamed)
  eq(BM.detectFormat('structure', read('POSCAR')), 'poscar', 'POSCAR without its name');
  eq(BM.detectFormat('run.out', read('relax.pwo')), 'qe-out', 'QE output named .out');
  eq(BM.detectFormat('vasp.log', read('OUTCAR')), 'outcar', 'OUTCAR renamed');
  for (const id of ['dcd', 'trr', 'xtc']) ok(BM.isBinaryFormat('traj.' + id), id + ' is a binary format');
});

test('builtin detection unchanged', () => {
  const ctx = { console };
  ctx.globalThis = ctx;
  vm.createContext(ctx);
  for (const f of ['elements.js', 'parsers.js']) vm.runInContext(fs.readFileSync(path.join(ROOT, 'src', f), 'utf8'), ctx, { filename: f });
  const base = ctx.BM;
  const dirs = [path.join(__dirname, 'fixtures'), path.join(__dirname, 'fixtures', 'review'), path.join(__dirname, 'fixtures', 'orbitals')];
  let n = 0;
  for (const d of dirs) {
    if (!fs.existsSync(d)) continue;
    for (const f of fs.readdirSync(d)) {
      const p = path.join(d, f);
      if (!fs.statSync(p).isFile() || /\.json$/.test(f)) continue;
      const text = fs.readFileSync(p, 'utf8');
      eq(BM.detectFormat(f, text), base.detectFormat(f, text), 'detection of ' + f);
      n++;
    }
  }
  ok(n >= 20, 'checked ' + n + ' builtin fixtures');
});

// ------------------------------------------------------------------ CIF / mmCIF
test('CIF NaCl Fm-3m', () => {
  const ds = load('nacl_fm3m.cif');
  const f = ds.frames[0];
  eq(f.Z.length, 8, 'conventional cell atoms');
  eq(f.Z.filter(z => z === 11).length, 4, 'Na count');
  eq(f.Z.filter(z => z === 17).length, 4, 'Cl count');
  arrNear(flat(ds.cell), [5.6402, 0, 0, 0, 5.6402, 0, 0, 0, 5.6402], 1e-9, 'cell');
  let dmin = Infinity;
  for (let i = 0; i < 8; i++) for (let j = 0; j < 8; j++) if (f.Z[i] === 11 && f.Z[j] === 17) dmin = Math.min(dmin, dist(f, i, j));
  near(dmin, 5.6402 / 2, 1e-9, 'shortest Na-Cl');
  ok(f.xyz.every(v => Math.abs(v) < 1e-9 || Math.abs(v - 2.8201) < 1e-9), 'positions on the 0 / a/2 lattice');
  eq(ds.meta.spaceGroup, 'F m -3 m', 'space group');
  ok(ds.warnings.some(w => /extra_block/.test(w)), 'other data blocks listed');
  eq(BM.structReaders.parseSymop('-y+1/2, x-y, 2/3+z').t[2], 2 / 3, 'symop translation');
  arrEq(BM.structReaders.parseSymop('-y+1/2, x-y, 2/3+z').R[1], [1, -1, 0], 'symop rotation row');
});
test('CIF P-1, labels, occupancy, text fields', () => {
  const ds = load('p1bar_labels.cif');
  const f = ds.frames[0];
  eq(ds.title, 'test_p1bar', 'first block with atoms');
  arrEq(f.Z, [26, 8, 8, 6, 6, 1, 1], 'Z from labels, Fe on the inversion centre merged');
  const c = ds.cell;
  near(Math.hypot(...c[2]), 7, 1e-9, '|c|');
  near(Math.acos((c[0][0] * c[2][0] + c[0][1] * c[2][1] + c[0][2] * c[2][2]) / (5 * 7)) / D2R, 100, 1e-9, 'beta');
  const inv = (x, y, z) => vadd(vs(c[0], x), vs(c[1], y), vs(c[2], z));
  arrNear(atom(f, 2), inv(0.9, 0.8, 0.7), 1e-9, 'inverted O wrapped into the cell');
  ok(f.occupancy && f.occupancy[3] === 0.5, 'occupancy kept');
  ok(ds.warnings.some(w => /global/.test(w)), 'global block listed');
  arrEq(f.atomNames.slice(0, 3), ['Fe1', 'O1', 'O1'], 'labels');
});
test('mmCIF models, chains, alt locs', () => {
  const ds = load('small.mmcif');
  eq(ds.format, 'mmcif', 'format');
  eq(ds.frames.length, 2, 'models');
  const f = ds.frames[0];
  arrEq(f.Z, [7, 6, 6, 8, 8, 8, 30], 'Z (first alt loc kept)');
  arrEq(f.chains, ['X', 'X', 'X', 'X', 'Y', 'W', 'W'], 'auth chains');
  arrEq(f.resIds, [10, 10, 11, 11, 20, 101, 201], 'auth residue numbers');
  arrEq(f.resNames, ['GLY', 'GLY', 'SER', 'SER', 'DA', 'HOH', 'ZN'], 'residue names');
  eq(f.atomNames[4], "O5'", 'quoted atom name');
  arrNear(atom(f, 2), [3, 2, 3], 1e-9, 'alt A position');
  near(ds.frames[1].xyz[0], 1.5, 1e-9, 'model 2 coordinates');
  eq(ds.frames[1].label, 'Model 2', 'model label');
  eq(ds.cell, null, 'unit cell 1x1x1 ignored');
  eq(ds.title, 'Tiny test peptide', 'title');
});

// ------------------------------------------------------------------ Mol2, CML, CJSON
test('Mol2 Tripos', () => {
  const ds = load('ligand.mol2');
  const f = ds.frames[0];
  eq(f.Z.length, 18, 'atoms');
  arrEq(f.Z.slice(12), [6, 8, 7, 1, 1, 1], 'formamide elements');
  eq(ds.bonds.length, 17, 'bonds (nc skipped)');
  arrEq(ds.bonds[0], [0, 1, 1.5], 'aromatic bond');
  arrEq(ds.bonds[12], [12, 13, 2], 'double bond');
  arrEq(ds.bonds[13], [12, 14, 1], 'amide bond = 1');
  near(ds.charges.Mol2[0], -0.115, 1e-12, 'partial charge');
  arrEq([f.resNames[0], f.resIds[0], f.resNames[12], f.resIds[12]], ['BEN', 1, 'FOR', 2], 'substructures as residues');
  arrEq([f.chains[0], f.chains[17]], ['A', 'B'], 'substructure chains');
  arrNear(flat(ds.cell), [20, 0, 0, 0, 21, 0, 0, 0, 22], 1e-9, 'CRYSIN cell');
  eq(ds.title, 'bz_formamide', 'title');
});
test('Mol2 GAFF types, several molecules', () => {
  const ds = load('gaff2.mol2');
  eq(ds.frames.length, 2, 'molecules as frames');
  arrEq(ds.frames[0].Z, [6, 8, 1, 1, 1, 1], 'elements from atom names');
  near(ds.frames[1].xyz[0], 0.1, 1e-9, 'second molecule');
  near(ds.charges.Mol2[1], -0.5988, 1e-12, 'charges');
  eq(ds.frames[0].chains, null, 'no chains');
});
test('CML attribute and array forms', () => {
  let ds = load('attr.cml');
  arrEq(ds.frames[0].Z, [6, 7, 1, 1, 1, 1, 1, 1], 'atoms (comment ignored)');
  eq(ds.bonds.length, 7, 'bonds');
  eq(ds.frames[0].formalCharges[1], 1, 'formal charge');
  arrNear(flat(ds.cell), [10, 0, 0, 0, 11, 0, 0, 0, 12], 1e-9, 'crystal');
  eq(ds.title, 'methylammonium', 'title');
  ds = load('array.cml');
  arrEq(ds.frames[0].Z, [6, 6, 1, 1, 1, 1], 'array-form atoms');
  arrEq(ds.bonds[0], [0, 1, 2], 'array-form double bond');
  eq(ds.bonds.length, 5, 'array-form bonds');
  near(ds.frames[0].xyz[3], 1.33, 1e-12, 'array-form x3');
  const cml1 = '<molecule id="m"><atomArray><atom id="a1"><string builtin="elementType">O</string><float builtin="x3">0</float>' +
    '<float builtin="y3">0</float><float builtin="z3">0</float></atom><atom id="a2"><string builtin="elementType">H</string>' +
    '<float builtin="x3">0.96</float><float builtin="y3">0</float><float builtin="z3">0</float></atom></atomArray><bondArray><bond>' +
    '<string builtin="atomRef">a1</string><string builtin="atomRef">a2</string><string builtin="order">1</string></bond></bondArray></molecule>';
  ds = BM.parseFile('old.cml', cml1);
  arrEq(ds.frames[0].Z, [8, 1], 'CML1 builtin children');
  near(ds.frames[0].xyz[3], 0.96, 1e-12, 'CML1 x3');
  arrEq(ds.bonds[0], [0, 1, 1], 'CML1 bond');
});
test('CJSON', () => {
  let ds = load('water.cjson');
  arrEq(ds.frames[0].Z, [8, 1, 1], 'elements');
  eq(ds.bonds.length, 2, 'bonds');
  arrEq(ds.charges.Mulliken, [-0.66, 0.33, 0.33], 'partial charges');
  eq(ds.vibrations.freqs.length, 3, 'frequencies');
  near(ds.vibrations.modes[0][2], 0.07, 1e-12, 'eigenvector');
  arrEq(ds.vibrations.ir, [70, 3, 45], 'intensities');
  eq(ds.title, 'water', 'name');
  ds = load('nacl_frac.cjson');
  const c = ds.cell;
  near(Math.hypot(...c[1]), 3.988, 1e-9, 'unitCell b');
  arrNear(atom(ds.frames[0], 1), vs(vadd(c[0], c[1], c[2]), 0.5), 1e-9, '3dFractional');
});

// ------------------------------------------------------------------ VASP
test('POSCAR / CONTCAR', () => {
  let ds = load('POSCAR');
  const f = ds.frames[0];
  arrEq(f.Z, [11, 11, 11, 11, 17, 17, 17, 17], 'species');
  arrNear(atom(f, 4), [2.82, 2.82, 2.82], 1e-9, 'Direct coordinates');
  eq(ds.meta.selectiveDynamics[0], 'FFF', 'selective dynamics flags');
  eq(ds.title, 'NaCl rocksalt conventional', 'comment line');
  ds = load('CONTCAR');
  arrNear(flat(ds.cell), [2, 0, 0, 0, 2, 0, 0, 0, 4], 1e-9, 'negative scale = volume');
  arrNear(ds.frames[0].xyz, [0, 0, 0, 1, 1, 2, 0.5, 0, 1], 1e-9, 'Cartesian scaled, velocities ignored');
  ds = load('si_v4.vasp');
  arrEq(ds.frames[0].Z, [14, 14], 'VASP 4 species from the comment');
  arrNear(atom(ds.frames[0], 1), [1.3575, 1.3575, 1.3575], 1e-9, 'VASP 4 position');
});
test('XDATCAR', () => {
  const ds = load('XDATCAR');
  eq(ds.frames.length, 2, 'configurations');
  near(ds.frames[1].xyz[0], 0.0543, 1e-9, 'frame 2');
  near(ds.cell[0][0], 5.43, 1e-12, 'cell');
  ok(!ds.frames[0].cell, 'constant cell not repeated per frame');
});
test('OUTCAR', () => {
  const ds = load('OUTCAR');
  eq(ds.frames.length, 2, 'ionic steps');
  arrEq(ds.frames[0].Z, [11, 17], 'species from TITEL');
  near(ds.frames[0].energy, -10.84651153 / HARTREE_EV, 1e-12, 'free energy TOTEN step 1');
  near(ds.frames[1].energy, -10.95 / HARTREE_EV, 1e-12, 'free energy TOTEN step 2');
  near(ds.frames[0].opt.maxForce, 0.1 * EVA_AU, 1e-12, 'max force');
  near(ds.frames[1].opt.maxForce, 0.05 * EVA_AU, 1e-12, 'max force step 2');
  near(ds.frames[0].opt.rmsForce, Math.sqrt(0.02 / 6) * EVA_AU, 1e-12, 'rms force');
  near(ds.frames[0].cell[0][1], 2.82, 1e-12, 'per-step cell 1');
  near(ds.cell[0][1], 2.83, 1e-12, 'final cell');
  near(ds.frames[1].xyz[0], 0.01, 1e-12, 'positions');
});

// ------------------------------------------------------------------ QE, CASTEP, XSF, CP2K
test('Quantum ESPRESSO input', () => {
  const alat = 10.2 * BOHR, fcc = [[-0.5, 0, 0.5], [0, 0.5, 0.5], [-0.5, 0.5, 0]];
  let ds = load('si_ibrav0.pwi');
  arrNear(flat(ds.cell), flat(fcc).map(v => v * alat), 1e-9, 'CELL_PARAMETERS alat');
  arrNear(atom(ds.frames[0], 1), [-0.25 * alat, 0.25 * alat, 0.25 * alat], 1e-9, 'crystal positions');
  arrEq(ds.frames[0].Z, [14, 14], 'species');
  eq(ds.title, 'si', 'prefix as title');
  ds = load('si_ibrav2.in');
  arrNear(flat(ds.cell), flat(fcc).map(v => v * alat), 1e-9, 'ibrav = 2');
  arrNear(atom(ds.frames[0], 1), [0.25 * alat, 0.25 * alat, 0.25 * alat], 1e-9, 'alat positions');
  ds = load('h2_bohr.pwi');
  near(ds.cell[0][0], 20 * BOHR, 1e-9, 'CELL_PARAMETERS bohr');
  near(dist(ds.frames[0], 0, 1), 1.4 * BOHR, 1e-9, 'bohr positions');
  arrEq(ds.frames[0].Z, [1, 1], 'labels H1');
});
test('Quantum ESPRESSO output', () => {
  const ds = load('relax.pwo'), alat = 10.2 * BOHR;
  eq(ds.frames.length, 2, 'ionic steps (final coordinates merged)');
  near(ds.frames[0].energy, -7.9, 1e-12, 'energy Ry -> Hartree');
  near(ds.frames[1].energy, -7.92, 1e-12, 'energy step 2');
  arrNear(atom(ds.frames[0], 1), [0.26 * alat, 0.25 * alat, 0.25 * alat], 1e-9, 'initial positions');
  arrNear(atom(ds.frames[1], 1), [0.255 * alat, 0.25 * alat, 0.25 * alat], 1e-9, 'final positions');
  ok(ds.frames[1].converged && ds.frames[1].label === 'Final', 'final coordinates flagged');
  near(ds.frames[0].opt.maxForce, 0.005, 1e-12, 'forces Ry/Bohr -> Ha/Bohr');
  near(ds.cell[0][0], -0.5 * alat, 1e-9, 'cell');
  arrEq(ds.frames[0].Z, [14, 14], 'species');
});
test('CASTEP .cell', () => {
  let ds = load('si.cell');
  const c = ds.cell;
  near(Math.hypot(...c[2]), 3.8396, 1e-9, 'LATTICE_ABC');
  arrNear(atom(ds.frames[0], 1), vs(vadd(c[0], c[1], c[2]), 0.25), 1e-9, 'POSITIONS_FRAC');
  arrEq(ds.frames[0].Z, [14, 14], 'labels Si:1');
  ds = load('h2o_abs.cell');
  near(ds.cell[1][1], 18 * BOHR, 1e-9, 'LATTICE_CART bohr');
  near(ds.frames[0].xyz[2], 0.1173, 1e-12, 'POSITIONS_ABS ang');
});
test('XSF / AXSF', () => {
  let ds = load('grid.xsf');
  arrEq(ds.frames[0].Z, [8, 1], 'atoms (forces ignored)');
  arrNear(flat(ds.cell), [3, 0, 0, 0, 4, 0, 0, 0, 5], 0, 'PRIMVEC');
  const v = ds.volumes[0];
  arrEq(v.dims, [3, 2, 2], 'grid dims');
  eq(v.data[(1 * 2 + 1) * 2 + 0], 11, 'grid value (1,1,0)');
  eq(v.data[(2 * 2 + 1) * 2 + 1], 112, 'grid value (2,1,1)');
  arrNear(flat(v.axes), [1.5, 0, 0, 0, 4, 0, 0, 0, 5], 1e-12, 'step vectors');
  eq(v.kind, 'density', 'kind');
  ds = load('anim.axsf');
  eq(ds.frames.length, 2, 'animation steps');
  near(ds.frames[1].xyz[3], 1.3, 1e-12, 'step 2');
  eq(ds.cell, null, 'molecule: no cell');
});
test('CP2K input', () => {
  const ds = load('h2o.inp');
  arrEq(ds.frames[0].Z, [8, 1, 1], 'KIND ELEMENT mapping');
  arrNear(flat(ds.cell), [10, 0, 0, 0, 11, 0, 0, 0, 12], 1e-9, 'ABC cell');
  eq(ds.title, 'water_test', 'project name');
});

// ------------------------------------------------------------------ Amber, PSF, LAMMPS
test('Amber prmtop + inpcrd', () => {
  let ds = load('water.prmtop');
  ok(ds.meta.topologyOnly, 'topologyOnly');
  const f = ds.frames[0];
  arrEq(f.Z, [8, 1, 1, 11], 'ATOMIC_NUMBER');
  arrEq(f.atomNames, ['O', 'H1', 'H2', 'Na+'], 'atom names');
  arrEq(f.resNames, ['WAT', 'WAT', 'WAT', 'Na+'], 'residue labels');
  arrEq(f.resIds, [1, 1, 1, 2], 'residue pointers');
  arrEq(ds.bonds, [[0, 1, 1], [0, 2, 1]], 'bonds');
  arrNear(ds.charges.Amber, [-0.834, 0.417, 0.417, 1], 1e-7, 'charges / 18.2223');
  arrNear(flat(ds.cell), [20, 0, 0, 0, 21, 0, 0, 0, 22], 1e-9, 'box');
  ok(f.xyz.every(v => v === 0), 'zero coordinates');
  eq(ds.title, 'water and sodium', 'title');
  ds = BM.parseFile('w.prmtop', read('water.prmtop').replace(/%FLAG ATOMIC_NUMBER\n%FORMAT\(10I8\)\n.*\n/, ''));
  arrEq(ds.frames[0].Z, [8, 1, 1, 11], 'elements from masses');
  ds = load('water.inpcrd');
  ok(ds.meta.needsTopology && ds.meta.natoms === 4 && ds.frames[0].Z === null, 'needsTopology');
  near(ds.frames[0].xyz[2], 0.1173, 1e-12, 'coordinates');
  arrNear(ds.frames[0].xyz.slice(9), [5, 5, 5], 1e-12, 'last atom');
  arrNear(flat(ds.cell), [20, 0, 0, 0, 21, 0, 0, 0, 22], 1e-9, 'box line');
});
test('CHARMM PSF', () => {
  const ds = load('two_seg.psf');
  const f = ds.frames[0];
  ok(ds.meta.topologyOnly, 'topologyOnly');
  arrEq(f.Z, [7, 1, 8, 11], 'elements from masses');
  arrEq(f.chains, ['PROA', 'PROA', 'WAT', 'IONS'], 'segments as chains');
  arrEq(f.resIds, [1, 1, 5, 7], 'residue ids');
  arrEq(f.resNames, ['GLY', 'GLY', 'TIP3', 'SOD'], 'residue names');
  arrEq(ds.bonds, [[0, 1, 1]], 'bonds');
  near(ds.charges.CHARMM[3], 1, 1e-12, 'charges');
  eq(ds.title, 'tiny test system', 'title from REMARKS');
});
test('LAMMPS data', () => {
  let ds = load('full.data');
  const f = ds.frames[0];
  arrEq(f.Z, [8, 1, 1, 8, 1, 1], 'elements from masses, sorted by id');
  arrNear(atom(f, 3), [3, 3, 3], 0, 'atom id 4');
  arrEq(ds.bonds, [[0, 1, 1], [0, 2, 1], [3, 4, 1], [3, 5, 1]], 'bonds by atom id');
  near(ds.charges.LAMMPS[0], -0.834, 1e-12, 'charges');
  arrEq(f.resIds, [1, 1, 1, 2, 2, 2], 'molecule ids');
  arrNear(flat(ds.cell), [10, 0, 0, 1, 10, 0, 0.5, -0.5, 12], 1e-12, 'tilted box');
  arrEq(ds.cellOrigin, [0, -1, 0], 'box origin');
  arrEq(f.atomNames.slice(0, 2), ['O', 'H'], 'type labels from Masses comments');
  eq(ds.meta.atomStyle, 'full', 'style');
  ds = load('atomic.data');
  arrEq(ds.frames[0].Z, [11, 17], 'atomic style without a hint');
  eq(ds.meta.atomStyle, 'atomic', 'inferred style');
  ds = BM.parseFile('c.data', read('atomic.data').replace('Atoms', 'Atoms # charge').replace('1 1 0.0', '1 1 1.0 0.0').replace('2 2 2.82', '2 2 -1.0 2.82'));
  arrNear(ds.charges.LAMMPS, [1, -1], 0, 'charge style');
});
test('LAMMPS dump', () => {
  const ds = load('traj.lammpstrj');
  eq(ds.frames.length, 2, 'timesteps');
  arrEq(ds.frames[0].Z, [8, 1, 1], 'element column, sorted by id');
  arrNear(ds.frames[0].xyz, [0, 0, 0, 1, 0, 0, 0.1, 1, 0], 1e-12, 'scaled triclinic coordinates');
  arrNear(flat(ds.cell), [10, 0, 0, 1, 10, 0, 0, 0, 10], 1e-12, 'box from bounds');
  near(ds.frames[1].xyz[2], 0.1, 1e-12, 'timestep 100');
  const noel = read('traj.lammpstrj').replace(/ element/g, '').replace(/^(\d) (\d) [OH] /gm, '$1 $2 ');
  const d2 = BM.parseFile('t.lammpstrj', noel);
  ok(d2.meta.needsTopology && d2.frames[0].Z === null && d2.meta.natoms === 3, 'type-only dump needs a topology');
});

// ------------------------------------------------------------------ DCD (independent byte writer)
function dcdBytes(o) {
  const le = o.le !== false, m8 = !!o.marker8, n = o.natoms;
  const i32 = v => { const b = Buffer.alloc(4); if (le) b.writeInt32LE(v); else b.writeInt32BE(v); return b; };
  const mark = v => (m8 ? (le ? Buffer.concat([i32(v), i32(0)]) : Buffer.concat([i32(0), i32(v)])) : i32(v));
  const rec = body => Buffer.concat([mark(body.length), body, mark(body.length)]);
  const f32 = arr => { const b = Buffer.alloc(4 * arr.length); arr.forEach((v, k) => (le ? b.writeFloatLE(v, 4 * k) : b.writeFloatBE(v, 4 * k))); return b; };
  const f64 = arr => { const b = Buffer.alloc(8 * arr.length); arr.forEach((v, k) => (le ? b.writeDoubleLE(v, 8 * k) : b.writeDoubleBE(v, 8 * k))); return b; };
  const hdr = Buffer.alloc(84);
  hdr.write('CORD', 0, 'latin1');
  const ic = new Array(20).fill(0);
  ic[0] = o.frames.length; ic[2] = 1; ic[8] = o.fixed ? o.fixed.length : 0; ic[10] = o.cell ? 1 : 0; ic[19] = 24;
  ic.forEach((v, k) => { if (k !== 9) { if (le) hdr.writeInt32LE(v, 4 + 4 * k); else hdr.writeInt32BE(v, 4 + 4 * k); } });
  if (le) hdr.writeFloatLE(0.002, 40); else hdr.writeFloatBE(0.002, 40);
  const parts = [rec(hdr)];
  const t = Buffer.alloc(84, 32);
  if (le) t.writeInt32LE(1, 0); else t.writeInt32BE(1, 0);
  t.write(o.title || 'test dcd', 4, 'latin1');
  parts.push(rec(t), rec(i32(n)));
  const free = o.fixed ? [...Array(n).keys()].filter(k => o.fixed.indexOf(k) < 0) : null;
  if (free) parts.push(rec(Buffer.concat(free.map(k => i32(k + 1)))));
  o.frames.forEach((fr, fi) => {
    if (o.cell) {
      const [a, b, c, al, be, ga] = o.cell;
      const ang = x => (o.degrees ? x : Math.cos(x * D2R));
      parts.push(rec(f64([a, ang(ga), b, ang(be), ang(al), c])));
    }
    const idx = free && fi > 0 ? free : [...Array(n).keys()];
    for (let d = 0; d < 3; d++) parts.push(rec(f32(idx.map(k => fr[3 * k + d]))));
  });
  return new Uint8Array(Buffer.concat(parts));
}
const F32 = a => Array.from(Float32Array.from(a));
test('DCD', () => {
  const frames = [[0, 0, 0, 1, 0, 0, 0, 1.5, 0], [0.1, 0, 0, 1.1, 0, 0, 0, 1.6, 0], [0.2, 0, 0, 1.2, 0, 0, 0, 1.7, 0.25]];
  const cell = [30, 31, 32, 90, 100, 90];
  const expCell = [30, 0, 0, 0, 31, 0, 32 * Math.cos(100 * D2R), 0, 32 * Math.sin(100 * D2R)];
  let ds = BM.parseFileBinary('a.dcd', dcdBytes({ le: true, natoms: 3, frames, cell }));
  eq(ds.frames.length, 3, 'LE frames');
  ok(ds.meta.needsTopology && ds.meta.natoms === 3 && ds.frames[0].Z === null, 'needsTopology');
  arrNear(ds.frames[2].xyz, F32(frames[2]), 0, 'LE coordinates');
  arrNear(flat(ds.cell), expCell, 1e-9, 'unit cell from cosines');
  eq(ds.title, 'test dcd', 'title');
  ds = BM.parseFileBinary('b.dcd', dcdBytes({ le: false, natoms: 3, frames }));
  eq(ds.frames.length, 3, 'BE frames');
  arrNear(ds.frames[1].xyz, F32(frames[1]), 0, 'BE coordinates');
  eq(ds.cell, null, 'no unit cell');
  ds = BM.parseFileBinary('c.dcd', dcdBytes({ le: true, natoms: 3, frames, cell, degrees: true }));
  arrNear(flat(ds.cell), expCell, 1e-9, 'unit cell in degrees');
  ds = BM.parseFileBinary('d.dcd', dcdBytes({ le: true, natoms: 3, frames, marker8: true }));
  arrNear(ds.frames[2].xyz, F32(frames[2]), 0, '64-bit record markers');
  ds = BM.parseFileBinary('e.dcd', dcdBytes({ le: true, natoms: 3, frames, fixed: [1] }));
  arrNear(ds.frames[2].xyz, F32([0.2, 0, 0, 1, 0, 0, 0, 1.7, 0.25]), 0, 'fixed atom keeps its first-frame position');
  const full = dcdBytes({ le: true, natoms: 3, frames, cell });
  ds = BM.parseFileBinary('f.dcd', full.subarray(0, full.length - 6));
  eq(ds.frames.length, 2, 'truncated file');
  ok(ds.warnings.some(w => /incomplete/.test(w)), 'truncation warning');
});

// ------------------------------------------------------------------ GROMACS TRR / XTC
function xdr() {
  const parts = [];
  const api = {
    int: v => { const b = Buffer.alloc(4); b.writeInt32BE(v); parts.push(b); return api; },
    f32: v => { const b = Buffer.alloc(4); b.writeFloatBE(v); parts.push(b); return api; },
    f64: v => { const b = Buffer.alloc(8); b.writeDoubleBE(v); parts.push(b); return api; },
    raw: b => { parts.push(Buffer.from(b)); return api; },
    bytes: () => new Uint8Array(Buffer.concat(parts))
  };
  return api;
}
function trrBytes(frames, n, dbl) {
  const w = xdr(), R = dbl ? 8 : 4, real = v => (dbl ? w.f64(v) : w.f32(v));
  frames.forEach((fr, k) => {
    w.int(1993).int(13).int(12).raw('GMX_trn_file');
    [0, 0, 9 * R, 0, 0, 0, 0, 3 * n * R, 0, 0, n, 10 * k, 0].forEach(v => w.int(v));
    real(0.5 * k); real(0);
    fr.box.forEach(real);
    fr.x.forEach(real);
  });
  return w.bytes();
}
test('TRR', () => {
  const box = [3, 0, 0, 0, 3.1, 0, 0.5, 0, 3.2];
  const frames = [{ box, x: [0.1, 0.2, 0.3, 1.1, 1.2, 1.3] }, { box, x: [0.15, 0.25, 0.35, 1.15, 1.25, 1.35] }];
  for (const dbl of [false, true]) {
    const ds = BM.parseFileBinary('t.trr', trrBytes(frames, 2, dbl));
    const p = dbl ? 'double: ' : 'single: ';
    eq(ds.frames.length, 2, p + 'frames');
    arrNear(ds.frames[1].xyz, frames[1].x.map(v => v * 10), dbl ? 1e-12 : 1e-5, p + 'coordinates nm -> Å');
    arrNear(flat(ds.cell), box.map(v => v * 10), dbl ? 1e-12 : 1e-5, p + 'box');
    near(ds.frames[1].time, 0.5, 1e-7, p + 'time');
    ok(ds.meta.needsTopology && ds.meta.natoms === 2, p + 'needsTopology');
  }
});
const MAGIC = [0, 0, 0, 0, 0, 0, 0, 0, 0, 8, 10, 12, 16, 20, 25, 32, 40, 50, 64, 80, 101, 128, 161, 203, 256, 322, 406, 512, 645, 812, 1024,
  1290, 1625, 2048, 2580, 3250, 4096, 5060, 6501, 8192, 10321, 13003, 16384, 20642, 26007, 32768, 41285, 52015, 65536, 82570, 104031,
  131072, 165140, 208063, 262144, 330280, 416127, 524287, 660561, 832255, 1048576, 1321122, 1664510, 2097152, 2642245, 3329021,
  4194304, 5284491, 6658042, 8388607, 10568983, 13316085, 16777216];
/** One XTC frame (coordinates in nm). Runs of nearby atoms use the small-integer code, with the
 * first two atoms of a run swapped as the format prescribes, and the small-code size is stepped up
 * and down so the decoder's adaptive index is exercised. */
function xtcFrame(x, n, prec, step, time, box) {
  const w = xdr().int(1995).int(n).int(step).f32(time);
  box.forEach(v => w.f32(v));
  w.int(n);
  if (n <= 9) { x.forEach(v => w.f32(v)); return w.bytes(); }
  const I = x.map(v => Math.round(v * prec));
  const mn = [0, 1, 2].map(d => Math.min(...I.filter((_, k) => k % 3 === d)));
  const mx = [0, 1, 2].map(d => Math.max(...I.filter((_, k) => k % 3 === d)));
  const size = [0, 1, 2].map(d => mx[d] - mn[d] + 1);
  const bitlen = v => (v === 0n ? 0 : v.toString(2).length);
  const large = size.some(s => s > 0xffffff);
  const bitsize = large ? 0 : bitlen(size.reduce((p, s) => p * BigInt(s), 1n));
  const bits = [];
  const put = (v, nb) => { const V = BigInt(v); for (let b = nb - 1; b >= 0; b--) bits.push(Number((V >> BigInt(b)) & 1n)); };
  const putInts = (nb, sz, nums) => {
    const V = (BigInt(nums[0]) * BigInt(sz[1]) + BigInt(nums[1])) * BigInt(sz[2]) + BigInt(nums[2]);
    const full = Math.floor(nb / 8), rem = nb % 8;
    for (let k = 0; k < full; k++) put((V >> BigInt(8 * k)) & 255n, 8);
    if (rem) put((V >> BigInt(8 * full)) & ((1n << BigInt(rem)) - 1n), rem);
  };
  const smallidx0 = 12, cycle = [0, 1, -1, 0, -1, 1];
  let smallidx = smallidx0, prevrun = 0, i = 0, g = 0;
  while (i < n) {
    const sm = Math.floor(MAGIC[smallidx] / 2), sz = MAGIC[smallidx];
    const fits = (a, b) => [0, 1, 2].every(d => { const v = I[3 * a + d] - I[3 * b + d] + sm; return v >= 0 && v < sz; });
    let k = 0;
    if (i + 1 < n && fits(i, i + 1)) { k = 1; while (k < 8 && i + k + 1 < n && fits(i + k + 1, k === 1 ? i : i + k)) k++; }
    let isSmaller = cycle[g++ % cycle.length];
    if (smallidx + isSmaller < 9 || smallidx + isSmaller > 20) isSmaller = 0;
    const L = k ? i + 1 : i;
    const lv = [0, 1, 2].map(d => I[3 * L + d] - mn[d]);
    if (large) [0, 1, 2].forEach(d => put(lv[d], bitlen(BigInt(size[d])))); else putInts(bitsize, size, lv);
    const run = 3 * k;
    if (run !== prevrun || isSmaller !== 0) { put(1, 1); put(run + isSmaller + 1, 5); prevrun = run; } else put(0, 1);
    for (let j = 0; j < k; j++) {
      const a = j === 0 ? i : i + 1 + j, ref = j === 0 ? i + 1 : (j === 1 ? i : i + j);
      putInts(smallidx, [sz, sz, sz], [0, 1, 2].map(d => I[3 * a + d] - I[3 * ref + d] + sm));
    }
    smallidx += isSmaller;
    i += k ? k + 1 : 1;
  }
  const nbytes = Math.ceil(bits.length / 8), data = Buffer.alloc((nbytes + 3) & ~3);
  bits.forEach((b, k) => { if (b) data[k >> 3] |= 0x80 >> (k & 7); });
  w.f32(prec);
  mn.forEach(v => w.int(v)); mx.forEach(v => w.int(v));
  w.int(smallidx0).int(nbytes).raw(data);
  return w.bytes();
}
test('XTC', () => {
  const box = [3, 0, 0, 0, 3, 0, 0, 0, 3];
  const mk = shift => {
    const x = [];
    for (let m = 0; m < 10; m++) {                            // water-like triplets: runs of small differences
      const o = [0.3 * (m % 4) + shift, 0.31 * Math.floor(m / 4) + 0.01 * m, 0.5 + 0.02 * m];
      x.push(...o, o[0] + 0.0957, o[1], o[2], o[0] - 0.024, o[1] + 0.0927, o[2]);
    }
    x.push(2.5, 2.6, 2.7, 0.05, 2.9, 0.4, 1.234567, 0.000123, 2.999);   // isolated atoms
    return x;
  };
  const f1 = mk(0), f2 = mk(0.0123), n = f1.length / 3;
  const bytes = Buffer.concat([xtcFrame(f1, n, 1000, 0, 0, box), xtcFrame(f2, n, 1000, 100, 0.2, box)]);
  const ds = BM.parseFileBinary('t.xtc', new Uint8Array(bytes));
  const q = (x, p) => Math.round(x * p) / p * 10;
  eq(ds.frames.length, 2, 'frames');
  arrNear(ds.frames[0].xyz, f1.map(v => q(v, 1000)), 1e-5, 'frame 1 decompressed');
  arrNear(ds.frames[1].xyz, f2.map(v => q(v, 1000)), 1e-5, 'frame 2 decompressed');
  near(ds.frames[1].time, 0.2, 1e-6, 'time');
  arrNear(flat(ds.cell), box.map(v => v * 10), 1e-5, 'box');
  ok(ds.meta.needsTopology && ds.meta.natoms === n, 'needsTopology');
  const small = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6];
  const d2 = BM.parseFileBinary('s.xtc', xtcFrame(small, 2, 1000, 0, 0, box));
  arrNear(d2.frames[0].xyz, small.map(v => v * 10), 1e-5, 'uncompressed frame (<= 9 atoms)');
  const wide = f1.slice();
  wide[0] = 20000.5;                                           // > 2^24 integer range: per-axis bit sizes
  const d3 = BM.parseFileBinary('w.xtc', xtcFrame(wide, n, 1000, 0, 0, box));
  arrNear(d3.frames[0].xyz, wide.map(v => q(v, 1000)), 1e-3, 'large coordinate range');
  const cut = BM.parseFileBinary('c.xtc', new Uint8Array(bytes.subarray(0, bytes.length - 20)));
  eq(cut.frames.length, 1, 'truncated file keeps complete frames');
});

// ---------------------------------------------------------------------- summary
console.log('\nBlueMolecule structure-reader tests: ' + pass + ' passed, ' + fail + ' failed');
if (fail) { for (const f of failures) console.log('  FAIL ' + f); process.exit(1); }
