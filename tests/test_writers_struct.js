#!/usr/bin/env node
/* BlueMolecule structure-writer tests (src/writers-struct.js).  Run:  node tests/test_writers_struct.js
 * Every writer is read back with its reader (builtin MOL/XYZ/GRO readers or src/readers-struct.js):
 * geometry within 1e-4 Å, cell (vectors, or parameters for CIF/CML/Mol2) and bonds preserved. */
'use strict';
const path = require('path');
const ROOT = path.join(__dirname, '..');
for (const f of ['elements', 'parsers', 'writers', 'chem', 'readers-struct', 'writers-struct']) require(path.join(ROOT, 'src', f + '.js'));
const BM = globalThis.BM;
const D2R = Math.PI / 180, HARTREE_EV = 27.211386245988;

let pass = 0, fail = 0, current = '';
const failures = [];
function ok(c, msg) { if (c) pass++; else { fail++; failures.push('[' + current + '] ' + msg); } }
function eq(a, b, msg) { ok(a === b, msg + ': got ' + JSON.stringify(a) + ', expected ' + JSON.stringify(b)); }
function has(text, s, msg) { ok(typeof text === 'string' && (s instanceof RegExp ? s.test(text) : text.indexOf(s) >= 0), (msg || 'contains') + ' ' + s); }
function arrNear(a, b, tol, msg) {
  if (!a || a.length !== b.length) { ok(false, msg + ': length ' + (a && a.length) + ' vs ' + b.length); return; }
  for (let i = 0; i < b.length; i++) if (!(Math.abs(a[i] - b[i]) <= tol)) { ok(false, msg + '[' + i + ']: got ' + a[i] + ', expected ' + b[i]); return; }
  pass++;
}
function arrEq(a, b, msg) { ok(JSON.stringify(Array.from(a || [])) === JSON.stringify(Array.from(b)), msg + ': got ' + JSON.stringify(a && Array.from(a)) + ', expected ' + JSON.stringify(Array.from(b))); }
function test(name, fn) { current = name; try { fn(); } catch (e) { fail++; failures.push('[' + name + '] threw: ' + (e && e.stack || e)); } }
const flat = c => c.reduce((a, v) => a.concat(Array.from(v)), []);
const normBonds = bs => bs.map(b => [Math.min(b[0], b[1]), Math.max(b[0], b[1]), b[2] || 1]).sort((p, q) => p[0] - q[0] || p[1] - q[1]);
const params = c => BM.structWriters.cellParams(c);
function cellOf(a, b, c, al, be, ga) {            // a along x, b in the xy plane (independent of the readers)
  const ca = Math.cos(al * D2R), cb = Math.cos(be * D2R), cg = Math.cos(ga * D2R), sg = Math.sin(ga * D2R);
  const cy = (ca - cb * cg) / sg;
  return [[a, 0, 0], [b * cg, b * sg, 0], [c * cb, c * cy, c * Math.sqrt(1 - cb * cb - cy * cy)]];
}

// ------------------------------------------------------------------ test structures
const PHENOL = (() => {
  const Z = [6, 6, 6, 6, 6, 6, 8, 1, 1, 1, 1, 1, 1], x = [];
  for (let k = 0; k < 6; k++) { const a = k * Math.PI / 3; x.push(1.3912345 * Math.cos(a), 1.3912345 * Math.sin(a), 0.0012345 * k); }
  x.push(2.7512345, 0.0001234, 0.0123456);
  for (let k = 1; k < 6; k++) { const a = k * Math.PI / 3; x.push(2.4712345 * Math.cos(a), 2.4712345 * Math.sin(a), -0.0011111 * k); }
  x.push(3.1012345, 0.9101234, 0.0234567);
  const bonds = [];
  for (let k = 0; k < 6; k++) bonds.push([k, (k + 1) % 6, 1.5]);
  bonds.push([0, 6, 1]);
  for (let k = 1; k < 6; k++) bonds.push([k, 6 + k, 1]);
  bonds.push([6, 12, 1]);
  return { Z, xyz: Float64Array.from(x), bonds, title: 'phenol test', atomNames: ['C1', 'C2', 'C3', 'C4', 'C5', 'C6', 'O1', 'H2', 'H3', 'H4', 'H5', 'H6', 'HO'],
    resNames: new Array(13).fill('PHL'), resIds: new Array(13).fill(7), chains: new Array(13).fill('B'),
    charges: { Mulliken: [0.12, -0.1, -0.11, -0.09, -0.1, -0.12, -0.55, 0.1, 0.1, 0.1, 0.1, 0.1, 0.35] } };
})();
const CELL = cellOf(5.1, 6.2, 7.3, 80, 95, 105);
const CRYST = (() => {
  const fr = [[0.01, 0.02, 0.03], [0.25, 0.3, 0.35], [0.5, 0.1, 0.9], [0.75, 0.7, 0.65]];
  const x = [];
  for (const f of fr) for (let d = 0; d < 3; d++) x.push(f[0] * CELL[0][d] + f[1] * CELL[1][d] + f[2] * CELL[2][d]);
  return { Z: [26, 8, 14, 8], xyz: Float64Array.from(x), cell: CELL, title: 'triclinic test', bonds: [] };
})();
const NH4 = { Z: [7, 1, 1, 1, 1], xyz: Float64Array.from([0, 0, 0, 0.59, 0.59, 0.59, -0.59, -0.59, 0.59, -0.59, 0.59, -0.59, 0.59, -0.59, -0.59]),
  bonds: [[0, 1, 1], [0, 2, 1], [0, 3, 1], [0, 4, 1]], fc: [1, 0, 0, 0, 0], title: 'ammonium', meta: { charge: 1 } };
const VOL = (() => {
  const dims = [3, 4, 5], data = new Float32Array(60);
  for (let a = 0; a < 3; a++) for (let b = 0; b < 4; b++) for (let c = 0; c < 5; c++) data[(a * 4 + b) * 5 + c] = (a * 7 + b * 3 + c) * 0.137 - 1;
  return { name: 'density', kind: 'density', origin: [-1.5, -2, -2.5], axes: [[0.5, 0, 0], [0, 0.45, 0], [0.1, 0, 0.4]], dims, data };
})();
function ctxOf(m, extra) {
  return Object.assign({ Z: m.Z, xyz: m.xyz, bonds: m.bonds || null, formalCharges: m.fc || null, cell: m.cell || null, title: m.title || '',
    atomNames: m.atomNames || null, resNames: m.resNames || null, resIds: m.resIds || null, chains: m.chains || null,
    frames: [{ Z: m.Z, xyz: m.xyz, energy: null }], frameIndex: 0, volume: null, charges: m.charges || null, meta: m.meta || {}, baseName: 'test' }, extra || {});
}
function write(id, ctx) {
  const f = BM.saveFormats.find(x => x.id === id);
  if (!f) throw new Error('save format ' + id + ' is not registered');
  const out = f.write(ctx);
  return typeof out === 'string' ? out : out.files[0].text;
}
const rt = (id, ctx, name) => { const text = write(id, ctx); return { text, ds: BM.parseFile(name, text) }; };
/** Shift applied by the writers to centre a cell-less structure in its box. */
function centring(xyz, L) {
  const mn = [Infinity, Infinity, Infinity], mx = [-Infinity, -Infinity, -Infinity];
  for (let i = 0; i < xyz.length; i++) { mn[i % 3] = Math.min(mn[i % 3], xyz[i]); mx[i % 3] = Math.max(mx[i % 3], xyz[i]); }
  return [0, 1, 2].map(d => L / 2 - (mn[d] + mx[d]) / 2);
}
const shifted = (xyz, s) => Array.from(xyz, (v, i) => v + s[i % 3]);

// ------------------------------------------------------------------ tests
test('registry', () => {
  for (const id of ['sdf', 'sdf-frames', 'mol2', 'cml', 'cjson', 'cif', 'extxyz', 'extxyz-frames', 'poscar', 'qe-in', 'cp2k-in', 'xsf', 'gro', 'lammps-data', 'smiles']) {
    const f = BM.saveFormats.find(x => x.id === id);
    ok(f && f.ext && typeof f.write === 'function', id + ' registered');
  }
  eq(BM.saveFormats.find(x => x.id === 'sdf-frames').needs[0], 'frames', 'sdf-frames needs frames');
});
test('SDF', () => {
  let { text, ds } = rt('sdf', ctxOf(PHENOL), 'p.sdf');
  arrNear(ds.frames[0].xyz, PHENOL.xyz, 1e-4, 'geometry');
  arrEq(normBonds(ds.bonds), normBonds(PHENOL.bonds), 'bonds with aromatic orders');
  eq(ds.title, 'phenol test', 'title');
  ({ text } = rt('sdf', ctxOf(NH4), 'n.sdf'));
  has(text, 'M  CHG  1   1   1', 'formal charge');
  const frames = [0, 1, 2].map(k => ({ Z: PHENOL.Z, xyz: PHENOL.xyz.map(v => v + 0.1 * k), energy: -307.1 - 0.01 * k }));
  text = write('sdf-frames', ctxOf(PHENOL, { frames }));
  const recs = text.split('$$$$\n').filter(s => s.trim());
  eq(recs.length, 3, 'one record per frame');
  recs.forEach((r, k) => arrNear(BM.parseFile('f.sdf', r + '$$$$\n').frames[0].xyz, frames[k].xyz, 1e-4, 'frame ' + (k + 1)));
  has(text, '>  <ENERGY>\n-307.1200000000', 'energy property');
  const big = { Z: new Array(1000).fill(2), xyz: Float64Array.from({ length: 3000 }, (_, k) => (k % 3 === 0 ? 3.1 * Math.floor(k / 3) : 0.5)), bonds: [] };
  ({ text, ds } = rt('sdf', ctxOf(big), 'big.sdf'));
  has(text, 'V3000', 'V3000 for more than 999 atoms');
  eq(ds.frames[0].Z.length, 1000, 'V3000 atoms read back');
  arrNear(ds.frames[0].xyz, big.xyz, 1e-4, 'V3000 geometry');
  ds = rt('sdf', ctxOf({ Z: [8, 1, 1], xyz: Float64Array.from([0, 0, 0.1173, 0, 0.7572, -0.4692, 0, -0.7572, -0.4692]) }), 'w.sdf').ds;
  eq(ds.bonds.length, 2, 'bonds guessed when the structure has none');
});
test('Mol2', () => {
  const { text, ds } = rt('mol2', ctxOf(PHENOL), 'p.mol2');
  arrNear(ds.frames[0].xyz, PHENOL.xyz, 1e-4, 'geometry');
  arrEq(ds.frames[0].Z, PHENOL.Z, 'elements');
  arrEq(normBonds(ds.bonds), normBonds(PHENOL.bonds), 'bonds (ar = 1.5)');
  arrNear(ds.charges.Mol2, PHENOL.charges.Mulliken, 1e-6, 'partial charges');
  arrEq([ds.frames[0].resNames[0], ds.frames[0].resIds[0], ds.frames[0].chains[0]], ['PHL', 7, 'B'], 'residue and chain');
  arrEq(ds.frames[0].atomNames, PHENOL.atomNames, 'atom names');
  has(text, ' C.ar ', 'aromatic carbon type');
  has(text, ' O.3 ', 'hydroxyl oxygen type');
  const c = rt('mol2', ctxOf(CRYST), 'c.mol2').ds;
  arrNear(params(c.cell), params(CELL), 1e-5, 'CRYSIN cell parameters');
  arrNear(c.frames[0].xyz, CRYST.xyz, 1e-4, 'crystal geometry');
});
test('CML', () => {
  let ds = rt('cml', ctxOf(PHENOL), 'p.cml').ds;
  arrNear(ds.frames[0].xyz, PHENOL.xyz, 1e-6, 'geometry');
  arrEq(normBonds(ds.bonds), normBonds(PHENOL.bonds), 'bonds (A = aromatic)');
  eq(ds.title, 'phenol test', 'title');
  ds = rt('cml', ctxOf(NH4), 'n.cml').ds;
  arrEq(ds.frames[0].formalCharges, [1, 0, 0, 0, 0], 'formal charges');
  ds = rt('cml', ctxOf(CRYST), 'c.cml').ds;
  arrNear(params(ds.cell), params(CELL), 1e-6, 'crystal parameters');
  arrNear(ds.frames[0].xyz, CRYST.xyz, 1e-6, 'crystal geometry');
});
test('CJSON', () => {
  const vib = { freqs: [1000, 1100], ir: [5, 6], modes: [new Float64Array(39).fill(0.1), new Float64Array(39).fill(-0.2)] };
  let ds = rt('cjson', ctxOf(PHENOL, { vibrations: vib }), 'p.cjson').ds;
  arrNear(ds.frames[0].xyz, PHENOL.xyz, 1e-8, 'geometry');
  arrEq(normBonds(ds.bonds).map(b => [b[0], b[1]]), normBonds(PHENOL.bonds).map(b => [b[0], b[1]]), 'bond pairs');
  const ring = ds.bonds.filter(b => b[0] < 6 && b[1] < 6).map(b => b[2]);
  ok(ring.every(o => o === 1 || o === 2) && ring.reduce((s, o) => s + o, 0) === 9, 'aromatic ring written as Kekule orders: ' + ring);
  arrNear(ds.charges.Mulliken, PHENOL.charges.Mulliken, 1e-12, 'partial charges');
  arrNear(ds.vibrations.freqs, [1000, 1100], 0, 'frequencies');
  near(ds.vibrations.modes[1][5], -0.2, 1e-12, 'eigenvectors');
  ds = rt('cjson', ctxOf(CRYST), 'c.cjson').ds;
  arrNear(flat(ds.cell), flat(CELL), 1e-9, 'cellVectors');
  arrNear(ds.frames[0].xyz, CRYST.xyz, 1e-8, 'crystal geometry');
});
function near(a, b, tol, msg) { ok(typeof a === 'number' && Math.abs(a - b) <= tol, msg + ': got ' + a + ', expected ' + b); }
test('CIF', () => {
  let { text, ds } = rt('cif', ctxOf(CRYST), 'c.cif');
  arrNear(ds.frames[0].xyz, CRYST.xyz, 1e-6, 'geometry');
  arrNear(params(ds.cell), params(CELL), 1e-7, 'cell parameters');
  arrEq(ds.frames[0].Z, CRYST.Z, 'elements');
  has(text, "'P 1'", 'P1');
  ({ text, ds } = rt('cif', ctxOf(PHENOL), 'p.cif'));
  has(text, /# No unit cell.*20 Angstrom cubic box/, 'box comment');
  arrNear(flat(ds.cell), [20, 0, 0, 0, 20, 0, 0, 0, 20], 1e-9, '20 Å box');
  arrNear(ds.frames[0].xyz, shifted(PHENOL.xyz, centring(PHENOL.xyz, 20)), 1e-6, 'centred geometry');
  const long = { Z: [6, 6], xyz: Float64Array.from([0, 0, 0, 25, 0, 0]) };
  has(write('cif', ctxOf(long)), /35 Angstrom cubic box/, 'larger box for a long structure (5 Angstrom vacuum)');
});
test('Extended XYZ', () => {
  let { text, ds } = rt('extxyz', ctxOf(CRYST), 'c.xyz');
  arrNear(ds.frames[0].xyz, CRYST.xyz, 1e-8, 'geometry');
  arrNear(flat(ds.cell), flat(CELL), 1e-8, 'Lattice');
  has(text, 'Properties=species:S:1:pos:R:3', 'Properties');
  has(text, 'pbc="T T T"', 'pbc');
  const frames = [0, 1].map(k => ({ Z: CRYST.Z, xyz: CRYST.xyz.map(v => v + 0.05 * k), energy: -10 - k }));
  ({ text, ds } = rt('extxyz-frames', ctxOf(CRYST, { frames }), 't.xyz'));
  eq(ds.frames.length, 2, 'frames');
  arrNear(ds.frames[1].xyz, frames[1].xyz, 1e-8, 'frame 2');
  has(text, 'energy=' + (-11 * HARTREE_EV).toFixed(8), 'energy in eV');
  near(ds.frames[1].energy, -11, 1e-7, 'extxyz energy (eV) read back in Hartree');
  near(BM.parseFile('p.xyz', '1\nenergy=-0.5 comment\nH 0 0 0\n').frames[0].energy, -0.5, 1e-12, 'plain XYZ energy stays Hartree');
  near(BM.parseFile('q.xyz', '1\nProperties=species:S:1:pos:R:3 energy=-27.211386245988 pbc="F F F"\nH 0 0 0\n').frames[0].energy, -1, 1e-9, 'extxyz energy= in eV');
  text = write('extxyz', ctxOf(PHENOL));
  ok(text.indexOf('Lattice') < 0 && text.indexOf('pbc="F F F"') > 0, 'molecule: no Lattice, pbc F');
});
test('POSCAR', () => {
  const mixed = { Z: [26, 8, 14, 8], xyz: CRYST.xyz, cell: CELL, title: 'mixed' };
  let { text, ds } = rt('poscar', ctxOf(mixed), 'POSCAR');
  arrEq(ds.frames[0].Z, [26, 8, 8, 14], 'species grouped');
  const order = [0, 1, 3, 2];
  arrNear(ds.frames[0].xyz, flat(order.map(i => Array.from(CRYST.xyz.slice(3 * i, 3 * i + 3)))), 1e-8, 'geometry (grouped order)');
  arrNear(flat(ds.cell), flat(CELL), 1e-10, 'lattice vectors');
  has(text, /^\s+Fe\s+O\s+Si\s*$/m, 'VASP 5 species line');
  ({ text, ds } = rt('poscar', ctxOf(PHENOL), 'CONTCAR'));
  has(text.split('\n')[0], /No unit cell/, 'box comment');
  arrNear(flat(ds.cell), [20, 0, 0, 0, 20, 0, 0, 0, 20], 1e-10, 'box');
});
test('Quantum ESPRESSO input', () => {
  let { text, ds } = rt('qe-in', ctxOf(CRYST), 'c.pwi');
  arrNear(ds.frames[0].xyz, CRYST.xyz, 1e-8, 'geometry');
  arrNear(flat(ds.cell), flat(CELL), 1e-8, 'CELL_PARAMETERS angstrom');
  arrEq(ds.frames[0].Z, CRYST.Z, 'elements');
  for (const s of ['ibrav = 0', 'nat = 4', 'ntyp = 3', 'ecutwfc', "occupations = 'smearing'", 'conv_thr', 'ATOMIC_SPECIES', 'Fe.UPF', 'ATOMIC_POSITIONS angstrom', 'K_POINTS automatic']) has(text, s);
  has(text, /\n\s+Fe\s+55\.845/, 'mass in ATOMIC_SPECIES');
  has(text, /K_POINTS automatic\n\s+6 5 5 0 0 0/, 'k-point grid from the cell');
  ({ text, ds } = rt('qe-in', ctxOf(PHENOL), 'p.pwi'));
  has(text, 'K_POINTS gamma', 'gamma point for a boxed molecule');
  has(text, /! No unit cell/, 'box comment');
  arrNear(ds.frames[0].xyz, shifted(PHENOL.xyz, centring(PHENOL.xyz, 20)), 1e-8, 'centred molecule');
  text = write('qe-in', ctxOf(NH4, { meta: { charge: 1, multiplicity: 3 } }));
  has(text, 'tot_charge = 1'); has(text, 'nspin = 2'); has(text, 'tot_magnetization = 2');
});
test('CP2K input', () => {
  let { text, ds } = rt('cp2k-in', ctxOf(CRYST), 'c.inp');
  arrNear(ds.frames[0].xyz, CRYST.xyz, 1e-8, 'geometry');
  arrNear(flat(ds.cell), flat(CELL), 1e-8, 'cell');
  for (const s of ['&FORCE_EVAL', '&DFT', '&SUBSYS', '&CELL', '&COORD', '&KIND Fe', 'BASIS_SET DZVP-MOLOPT-SR-GTH', 'POTENTIAL GTH-PBE', 'PERIODIC XYZ', '&END FORCE_EVAL']) has(text, s);
  ({ text, ds } = rt('cp2k-in', ctxOf(PHENOL), 'p.inp'));
  has(text, 'PERIODIC NONE'); has(text, 'POISSON_SOLVER WAVELET');
  arrNear(ds.frames[0].xyz, shifted(PHENOL.xyz, centring(PHENOL.xyz, 20)), 1e-8, 'centred molecule');
});
test('XSF', () => {
  let ds = rt('xsf', ctxOf(CRYST, { volume: VOL }), 'c.xsf').ds;
  arrNear(ds.frames[0].xyz, CRYST.xyz, 1e-8, 'geometry');
  arrNear(flat(ds.cell), flat(CELL), 1e-8, 'PRIMVEC');
  const v = ds.volumes[0];
  arrEq(v.dims, VOL.dims, 'grid dims');
  arrNear(v.origin, VOL.origin, 1e-9, 'origin');
  arrNear(flat(v.axes), flat(VOL.axes), 1e-9, 'step vectors');
  arrNear(v.data, VOL.data, 1e-6, 'grid values');
  ds = rt('xsf', ctxOf(PHENOL), 'p.xsf').ds;
  arrNear(ds.frames[0].xyz, PHENOL.xyz, 1e-8, 'molecule geometry');
  eq(ds.cell, null, 'molecule: no cell');
});
test('GRO', () => {
  let { text, ds } = rt('gro', ctxOf(Object.assign({}, PHENOL, { cell: CELL })), 'p.gro');
  arrNear(ds.frames[0].xyz, PHENOL.xyz, 1e-4, 'geometry (nm, 7 decimals)');
  arrNear(flat(ds.cell), flat(CELL), 1e-5, 'triclinic box');
  arrEq(ds.frames[0].Z, PHENOL.Z, 'elements from atom names');
  arrEq(ds.frames[0].atomNames, PHENOL.atomNames, 'atom names');
  arrEq([ds.frames[0].resNames[0], ds.frames[0].resIds[0]], ['PHL', 7], 'residue');
  ds = rt('gro', ctxOf(CRYST), 'c.gro').ds;
  arrEq(ds.frames[0].Z, CRYST.Z, 'generated atom names keep the elements');
  ({ text, ds } = rt('gro', ctxOf(PHENOL), 'm.gro'));
  has(text.split('\n')[0], /No unit cell/, 'box comment in the title');
  arrNear(ds.frames[0].xyz, shifted(PHENOL.xyz, centring(PHENOL.xyz, 20)), 1e-4, 'centred molecule');
  eq(text.trim().split('\n').pop().trim(), '2.0000000   2.0000000   2.0000000', 'cubic box line');
});
test('LAMMPS data', () => {
  let ds = rt('lammps-data', ctxOf(CRYST), 'c.data').ds;
  arrNear(ds.frames[0].xyz, CRYST.xyz, 1e-8, 'geometry');
  arrNear(flat(ds.cell), flat(CELL), 1e-8, 'restricted triclinic box');
  arrEq(ds.frames[0].Z, CRYST.Z, 'elements from masses');
  const R = (() => { const t = 0.3, c = Math.cos(t), s = Math.sin(t); return [[c, -s, 0], [s, c, 0], [0, 0, 1]]; })();
  const rot = v => [0, 1, 2].map(i => R[i][0] * v[0] + R[i][1] * v[1] + R[i][2] * v[2]);
  const rc = CELL.map(rot), rx = [];
  for (let i = 0; i < 4; i++) rx.push(...rot(Array.from(CRYST.xyz.slice(3 * i, 3 * i + 3))));
  ds = rt('lammps-data', ctxOf({ Z: CRYST.Z, xyz: Float64Array.from(rx), cell: rc }), 'r.data').ds;
  ok(Math.abs(ds.cell[0][1]) < 1e-12 && Math.abs(ds.cell[0][2]) < 1e-12 && Math.abs(ds.cell[1][2]) < 1e-12, 'rotated cell written in LAMMPS orientation');
  arrNear(params(ds.cell), params(CELL), 1e-8, 'cell parameters kept');
  arrNear(ds.frames[0].xyz, CRYST.xyz, 1e-8, 'atoms rotated with the cell');
  ds = rt('lammps-data', ctxOf(PHENOL), 'p.data').ds;
  arrNear(ds.frames[0].xyz, PHENOL.xyz, 1e-8, 'molecule coordinates unchanged');
  const o = ds.cellOrigin;
  ok(Array.from(PHENOL.xyz).every((v, i) => v > o[i % 3] && v < o[i % 3] + 20), 'box around the molecule');
});
test('SMILES', () => {
  const s = write('smiles', ctxOf(PHENOL)).split('\t')[0].trim();
  const p = BM.chem.parseSmiles(s);
  eq(p.atoms.length, 7, 'heavy atoms in ' + s);
  eq(p.bonds.length, 7, 'heavy-atom bonds');
  eq(p.atoms.filter(a => a.Z === 8).length, 1, 'one oxygen');
  eq(write('smiles', ctxOf(NH4)).trim(), '[NH4+]\tammonium', 'charged atom and title');
  const guessed = write('smiles', ctxOf({ Z: PHENOL.Z, xyz: PHENOL.xyz })).trim();
  eq(BM.chem.parseSmiles(guessed).atoms.length, 7, 'bonds perceived when none are given: ' + guessed);
});

// ---------------------------------------------------------------------- summary
console.log('\nBlueMolecule structure-writer tests: ' + pass + ' passed, ' + fail + ' failed');
if (fail) { for (const f of failures) console.log('  FAIL ' + f); process.exit(1); }
