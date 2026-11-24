#!/usr/bin/env node
/* BlueMolecule tests for the 2.5 building tools: crystal tools (src/crystal.js: lattice parameters, supercells,
 * surface slabs, cell edits, boxes) and the conformer search (src/conformers.js).
 * Run:  node tests/test_crystal_conformers.js */
'use strict';
const fs = require('fs');
const path = require('path');
const ROOT = path.join(__dirname, '..');
for (const m of ['elements', 'bigfile', 'parsers', 'condense', 'readers-struct', 'writers', 'readers-qc', 'geom', 'chem', 'ff', 'analysis', 'crystal', 'conformers']) {
  require(path.join(ROOT, 'src', m + '.js'));
}
const BM = globalThis.BM, C = BM.crystal, K = BM.conformers;

let pass = 0, fail = 0, current = '';
const failures = [];
function ok(c, msg) { if (c) pass++; else { fail++; failures.push('[' + current + '] ' + msg); } }
function eq(a, b, msg) { ok(a === b, msg + ': got ' + JSON.stringify(a) + ', expected ' + JSON.stringify(b)); }
function near(a, b, tol, msg) { ok(typeof a === 'number' && Math.abs(a - b) <= tol, msg + ': got ' + a + ', expected ' + b); }
function test(name, fn) { current = name; try { fn(); } catch (e) { fail++; failures.push('[' + name + '] threw: ' + (e && e.stack || e)); } }
async function testAsync(name, fn) { current = name; try { await fn(); } catch (e) { fail++; failures.push('[' + name + '] threw: ' + (e && e.stack || e)); } }
let seed = 4242;
const rnd = () => (seed = (seed * 16807) % 2147483647) / 2147483647;

// rock salt from the CIF fixture (the symmetry-expanded conventional cell)
const nacl = BM.parseFile('nacl_fm3m.cif', fs.readFileSync(path.join(__dirname, 'fixtures', 'struct', 'nacl_fm3m.cif'), 'utf8'));
const NA = { Z: nacl.frames[0].Z, xyz: nacl.frames[0].xyz }, CELL = nacl.cell || nacl.frames[0].cell;
const A = C.params(CELL).a;

/** The shortest distance between atoms, with the images of the in-plane (or all) cell vectors. */
function minDist(xyz, cell, dims) {
  const n = xyz.length / 3, sh = [];
  const r = k => (dims.indexOf(k) >= 0 ? [-1, 0, 1] : [0]);
  for (const i of r(0)) for (const j of r(1)) for (const k of r(2)) sh.push([0, 1, 2].map(q => i * cell[0][q] + j * cell[1][q] + k * cell[2][q]));
  let m = Infinity;
  for (let p = 0; p < n; p++) for (let q = 0; q < n; q++) for (const s of sh) {
    if (p === q && !s.some(Boolean)) continue;
    const d = Math.hypot(xyz[3 * p] - xyz[3 * q] - s[0], xyz[3 * p + 1] - xyz[3 * q + 1] - s[1], xyz[3 * p + 2] - xyz[3 * q + 2] - s[2]);
    if (d < m) m = d;
  }
  return m;
}

test('lattice parameters', () => {
  for (let t = 0; t < 20; t++) {
    const p = { a: 2 + 8 * rnd(), b: 2 + 8 * rnd(), c: 2 + 8 * rnd(), alpha: 60 + 60 * rnd(), beta: 60 + 60 * rnd(), gamma: 60 + 60 * rnd() };
    let cell;
    try { cell = C.fromParams(p.a, p.b, p.c, p.alpha, p.beta, p.gamma); } catch (e) { continue; }
    const q = C.params(cell);
    ok(['a', 'b', 'c', 'alpha', 'beta', 'gamma'].every(k => Math.abs(q[k] - p[k]) < 1e-9), 'round trip ' + JSON.stringify(p));
    ok(cell[0][1] === 0 && cell[0][2] === 0 && cell[1][2] === 0, 'a along x, b in the xy plane');
    const xyz = Float64Array.from({ length: 30 }, () => 20 * rnd() - 10), back = C.toCart(cell, C.toFrac(cell, xyz));
    ok(xyz.every((v, i) => Math.abs(v - back[i]) < 1e-9), 'fractional and Cartesian coordinates');
    const w = C.toFrac(cell, C.wrap(cell, xyz));
    ok(w.every(v => v >= -1e-9 && v < 1), 'wrapped into the cell');
  }
  let err = ''; try { C.fromParams(5, 5, 5, 10, 10, 170); } catch (e) { err = e.message; }
  ok(/incompatible/.test(err), 'impossible angles: ' + err);
  err = ''; try { C.params([[1, 0, 0], [2, 0, 0], [0, 0, 1]]); } catch (e) { err = e.message; }
  ok(/no volume/.test(err), 'flat cell: ' + err);
  near(A, 5.64, 0.01, 'rock salt from the CIF fixture'); eq(NA.Z.length, 8, 'eight atoms');
});

test('supercells', () => {
  const s = C.supercell(NA, CELL, [2, 3, 1]);
  eq(s.Z.length, 48, 'atoms'); near(C.params(s.cell).b, 3 * A, 1e-9, 'b tripled'); near(C.params(s.cell).volume, 6 * A ** 3, 1e-6, 'volume');
  eq(s.map.slice(0, 8).join(), '0,1,2,3,4,5,6,7', 'the first copy is the original');
  near(s.xyz[3 * 8 + 1] - NA.xyz[1], A, 1e-9, 'second copy shifted by b');
  near(s.xyz[3 * 24] - NA.xyz[0], A, 1e-9, 'fourth copy shifted by a');
  near(minDist(s.xyz, s.cell, [0, 1, 2]), A / 2, 1e-6, 'no overlapping atoms (shortest Na–Cl = a/2)');
  const r = C.supercell(NA, CELL, [[1, 1, 0], [-1, 1, 0], [0, 0, 1]]);
  eq(r.Z.length, 16, '√2 × √2 cell: twice the atoms'); near(C.params(r.cell).a, A * Math.SQRT2, 1e-9, 'its a');
  near(minDist(r.xyz, r.cell, [0, 1, 2]), A / 2, 1e-6, 'no overlapping atoms');
  ok(C.toFrac(r.cell, r.xyz).every(v => v > -1e-6 && v < 1), 'atoms inside the new cell');
  let err = ''; try { C.supercell(NA, CELL, [0, 1, 1]); } catch (e) { err = e.message; }
  ok(/positive/.test(err), 'zero repeats: ' + err);
  err = ''; try { C.supercell(NA, CELL, [1.5, 1, 1]); } catch (e) { err = e.message; }
  ok(/whole numbers/.test(err), 'fractions: ' + err);
});

test('slabs of rock salt', () => {
  for (const hkl of [[1, 0, 0], [0, 1, 0], [0, 0, 1], [1, 1, 0], [1, -1, 0], [1, 1, 1], [-1, 1, 1], [2, 1, 0], [3, 1, 1], [2, 2, 0]]) {
    const s = C.slab(NA, CELL, hkl, { layers: 3, vacuum: 12 }), p = C.params(s.cell);
    const red = hkl.map(x => x / (hkl[0] === 2 && hkl[1] === 2 ? 2 : 1));
    const d = A / Math.hypot(...red);
    const tag = hkl.join(' ');
    eq(s.Z.length, 24, tag + ': three layers of eight atoms');
    near(s.d, d, 1e-9, tag + ': layer spacing a/√(h²+k²+l²)');
    near(p.alpha, 90, 1e-9, tag + ': α'); near(p.beta, 90, 1e-9, tag + ': β');
    ok(Math.abs(s.cell[2][0]) < 1e-12 && Math.abs(s.cell[2][1]) < 1e-12, tag + ': third vector along z');
    near(p.c, 3 * d + 12, 1e-9, tag + ': c = layers × spacing + vacuum');
    near(Math.abs(s.cell[0][1]) + Math.abs(s.cell[0][2]), 0, 1e-12, tag + ': a along x');
    const f = C.toFrac(s.cell, s.xyz);
    ok(f.every(v => v > -1e-6 && v < 1 + 1e-6), tag + ': atoms inside the cell');
    let zlo = Infinity, zhi = -Infinity;
    for (let i = 0; i < s.Z.length; i++) { zlo = Math.min(zlo, s.xyz[3 * i + 2]); zhi = Math.max(zhi, s.xyz[3 * i + 2]); }
    ok(zhi - zlo < 3 * d, tag + ': thinner than three layers');
    near((zlo + zhi) / 2, p.c / 2, 1e-6, tag + ': centred in the cell');
    near(minDist(s.xyz, s.cell, [0, 1]), A / 2, 1e-6, tag + ': the bulk distances (shortest Na–Cl = a/2)');
    eq(s.Z.filter(z => z === 11).length, 12, tag + ': as many Na as Cl');
    // a slab of the whole stoichiometry: every layer has the same Na/Cl ratio as the crystal
  }
  // the (111) surface cell is hexagonal
  const h = C.params(C.slab(NA, CELL, [1, 1, 1], { layers: 1 }).cell);
  near(h.gamma, 120, 1e-6, '(111): γ = 120°'); near(h.a, h.b, 1e-9, '(111): a = b');
  // termination shift: the other layer at the bottom
  const s0 = C.slab(NA, CELL, [1, 1, 1], { layers: 2, vacuum: 10 }), s1 = C.slab(NA, CELL, [1, 1, 1], { layers: 2, vacuum: 10, shift: 0.5 });
  const bottom = s => { let z = Infinity, k = -1; for (let i = 0; i < s.Z.length; i++) if (s.xyz[3 * i + 2] < z - 1e-6) { z = s.xyz[3 * i + 2]; k = i; } return s.Z[k]; };
  ok(bottom(s0) !== bottom(s1), '(111): the shift changes the terminating element (' + bottom(s0) + ' → ' + bottom(s1) + ')');
  let err = ''; try { C.slab(NA, CELL, [0, 0, 0], {}); } catch (e) { err = e.message; }
  ok(/zero/.test(err), 'Miller indices all zero: ' + err);
});

test('slabs of a triclinic crystal', () => {
  // one atom in a triclinic cell: the layer spacing is 1 / |h b1 + k b2 + l b3| (reciprocal vectors without 2π)
  const cell = C.fromParams(3.1, 4.3, 5.2, 77, 98, 106), inv = C.inverse(cell);
  const atoms = { Z: [6, 8], xyz: Float64Array.from([0.3, 0.2, 0.1, 1.2, 1.1, 0.9]) };
  for (const hkl of [[1, 0, 0], [0, 1, 1], [1, -1, 2], [2, 1, -1], [1, 2, 3]]) {
    const s = C.slab(atoms, cell, hkl, { layers: 2, vacuum: 5 });
    const G = [0, 1, 2].map(q => hkl[0] * inv[q][0] + hkl[1] * inv[q][1] + hkl[2] * inv[q][2]);
    near(s.d, 1 / Math.hypot(...G), 1e-9, hkl.join(' ') + ': layer spacing');
    eq(s.Z.length, 4, hkl.join(' ') + ': two layers of the two atoms');
    near(C.params(s.cell).volume, 5 * C.params(s.cell).a * C.params(s.cell).b * Math.sin(C.params(s.cell).gamma * Math.PI / 180) + 2 * C.params(cell).volume, 1e-6, hkl.join(' ') + ': volume = two cells + the vacuum');
    // the C–O distance of the crystal is kept (the pair is not split across the slab's side faces)
    const B = C.planeBasis(cell, ...hkl), det = B[0][0] * (B[1][1] * B[2][2] - B[1][2] * B[2][1]) - B[0][1] * (B[1][0] * B[2][2] - B[1][2] * B[2][0]) + B[0][2] * (B[1][0] * B[2][1] - B[1][1] * B[2][0]);
    eq(Math.abs(det), 1, hkl.join(' ') + ': unimodular change of basis');
    ok(B.flat().every(Number.isInteger), hkl.join(' ') + ': integer vectors');
    eq(hkl[0] * B[2][0] + hkl[1] * B[2][1] + hkl[2] * B[2][2], 1, hkl.join(' ') + ': third vector one layer out of the plane');
    ok([0, 1].every(k => hkl[0] * B[k][0] + hkl[1] * B[k][1] + hkl[2] * B[k][2] === 0), hkl.join(' ') + ': two vectors in the plane');
  }
});

test('cell edits and boxes', () => {
  const nc = C.fromParams(6, 6, 6, 90, 90, 90);
  const scaled = C.setCell(NA, CELL, nc, true), kept = C.setCell(NA, CELL, nc, false);
  const f0 = C.toFrac(CELL, NA.xyz), f1 = C.toFrac(nc, scaled);
  ok(f0.every((v, i) => Math.abs(v - f1[i]) < 1e-9), 'scaled: the same fractional coordinates');
  ok(Array.from(kept).every((v, i) => v === NA.xyz[i]), 'kept: the same positions');
  // a turned cell keeps its orientation when its parameters change
  const th = 0.3, R = [[Math.cos(th), Math.sin(th), 0], [-Math.sin(th), Math.cos(th), 0], [0, 0, 1]];
  const turned = C.fromParams(4, 5, 6, 80, 90, 100).map(v => [0, 1, 2].map(j => v[0] * R[0][j] + v[1] * R[1][j] + v[2] * R[2][j]));
  const same = C.oriented(turned, C.fromParams(4, 5, 6, 80, 90, 100));
  ok(same.every((v, i) => v.every((x, j) => Math.abs(x - turned[i][j]) < 1e-9)), 'unchanged parameters: the same cell');
  const longer = C.oriented(turned, C.fromParams(8, 5, 6, 80, 90, 100));
  ok(longer[0].every((x, j) => Math.abs(x - 2 * turned[0][j]) < 1e-9), 'a doubled along the old a');
  const mol = { Z: [8, 1, 1], xyz: Float64Array.from([0, 0, 0.12, 0, 0.76, -0.47, 0, -0.76, -0.47]) };
  const b = C.box(mol, 5), p = C.params(b.cell);
  near(p.a, 10, 1e-9, 'box x'); near(p.b, 11.52, 1e-9, 'box y'); near(p.c, 10.59, 1e-9, 'box z');
  const cx = (b.xyz[1] + b.xyz[4] + b.xyz[7]) / 3;
  ok(Math.abs((Math.max(b.xyz[1], b.xyz[4], b.xyz[7]) + Math.min(b.xyz[1], b.xyz[4], b.xyz[7])) / 2 - p.b / 2) < 1e-9 && cx > 0, 'molecule centred');
});

test('rotatable bonds', () => {
  const R = smi => K.rotors(BM.chem.smilesTo3D(smi).mol || BM.chem.smilesTo3D(smi));
  for (const [smi, n, why] of [['CC', 0, 'ethane'], ['CCCC', 1, 'butane'], ['CCCCCC', 3, 'hexane'], ['C1CCCCC1', 0, 'cyclohexane: a ring'],
    ['CC(=O)NC', 0, 'N-methylacetamide: the amide bond is not turned'], ['CC#CC', 0, '2-butyne: linear'], ['c1ccccc1-c1ccccc1', 1, 'biphenyl'],
    ['OCCO', 1, 'ethylene glycol'], ['CCOC(=O)C', 2, 'ethyl acetate (C–O–C and O–C(=O))']]) {
    const r = R(smi);
    eq(r.length, n, why);
  }
  eq(R('c1ccccc1-c1ccccc1')[0].step, 180, 'biphenyl: two planar ends, turned by 180°');
  eq(R('CCCC')[0].step, 120, 'butane: 120°');
});

(async () => {
  await testAsync('butane: anti below gauche', async () => {
    const m = BM.chem.smilesTo3D('CCCC'), mol = m.mol || m;
    const r = await K.search(mol, {}).promise;
    ok(r.conformers.length >= 2, 'two conformers at least: ' + r.conformers.length);
    eq(r.mode, 'systematic', 'systematic'); eq(r.tried, 3, 'three structures relaxed');
    const C4 = mol.Z.map((z, i) => i).filter(i => mol.Z[i] === 6);
    const dih = x => Math.abs(K.dihedral(x, C4[0], C4[1], C4[2], C4[3])) * 180 / Math.PI;
    near(dih(r.conformers[0].xyz), 180, 5, 'the lowest is anti');
    near(dih(r.conformers[1].xyz), 65, 15, 'the next is gauche');
    ok(r.conformers[1].rel > 0 && r.conformers[0].rel === 0, 'energies relative to the lowest');
    ok(r.conformers.every((c, i) => i === 0 || c.energy >= r.conformers[i - 1].energy), 'lowest first');
  });
  await testAsync('hexane: distinct, within the window, mirror images merged', async () => {
    const m = BM.chem.smilesTo3D('CCCCCC'), mol = m.mol || m;
    const r = await K.search(mol, { rmsd: 0.5, window: 3 }).promise;
    eq(r.tried, 27, 'every combination of three rotors');
    ok(r.conformers.length >= 4 && r.conformers.length <= 15, 'kept: ' + r.conformers.length);
    ok(r.conformers.every(c => c.rel <= 3 + 1e-9), 'within 3 kcal/mol');
    let close = 0;
    for (let i = 0; i < r.conformers.length; i++) for (let j = i + 1; j < r.conformers.length; j++) if (K.heavyRmsd(mol.Z, r.conformers[i].xyz, r.conformers[j].xyz) < 0.5) close++;
    eq(close, 0, 'pairwise RMSD above the threshold');
    const few = await K.search(mol, { maxConformers: 2 }).promise;
    eq(few.conformers.length, 2, 'at most maxConformers');
    const tight = await K.search(mol, { window: 0.01 }).promise;
    ok(tight.conformers.every(c => c.rel <= 0.01), 'a narrow window');
  });
  await testAsync('many rotors: random combinations; cancel', async () => {
    const m = BM.chem.smilesTo3D('CCCCCCCCCC'), mol = m.mol || m;
    const r = await K.search(mol, { samples: 40, maxSystematic: 100 }).promise;
    eq(r.mode, 'random', 'random combinations above maxSystematic');
    ok(r.tried <= 41 && r.tried >= 2, 'tried ' + r.tried);
    const j = K.search(mol, {});
    j.cancel();
    let msg = ''; try { await j.promise; } catch (e) { msg = e.message; }
    eq(msg, 'cancelled', 'cancelled');
  });
  console.log('BlueMolecule crystal and conformer tests: ' + pass + ' passed, ' + fail + ' failed');
  if (fail) { for (const f of failures.slice(0, 40)) console.log('  FAIL ' + f); process.exit(1); }
})();
