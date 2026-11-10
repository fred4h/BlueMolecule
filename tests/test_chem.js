#!/usr/bin/env node
/* BlueMolecule chemistry tests (src/chem.js).  Run:  node tests/test_chem.js
 * Independent of src/ff.js: BM.ff is never loaded here (a stub is installed only to test the hand-off). */
'use strict';
const path = require('path');
const ROOT = path.join(__dirname, '..');
require(path.join(ROOT, 'src', 'elements.js'));
require(path.join(ROOT, 'src', 'geom.js'));
require(path.join(ROOT, 'src', 'chem.js'));
const BM = globalThis.BM;
delete BM.ff;
const C = BM.chem, G = BM.geom;

let pass = 0, fail = 0, current = '';
const failures = [];
function ok(cond, msg) { if (cond) pass++; else { fail++; failures.push('[' + current + '] ' + msg); } }
function eq(a, b, msg) { ok(a === b, msg + ': got ' + JSON.stringify(a) + ', expected ' + JSON.stringify(b)); }
function near(a, b, tol, msg) { ok(typeof a === 'number' && Math.abs(a - b) <= tol, msg + ': got ' + a + ', expected ' + b + ' +/- ' + tol); }
function test(name, fn) {
  current = name;
  try { fn(); } catch (e) { fail++; failures.push('[' + name + '] threw: ' + (e && e.stack || e)); }
}
function throwsMsg(fn, re, msg) {
  try { fn(); ok(false, msg + ': no error thrown'); } catch (e) { ok(e instanceof Error && re.test(e.message), msg + ': message "' + (e && e.message) + '"'); }
}

// ------------------------------------------------------------------ helpers
const formulaOfParse = p => { const Z = []; p.atoms.forEach(a => { Z.push(a.Z); for (let k = 0; k < a.hCount; k++) Z.push(1); }); return G.formula(Z); };
const orderCounts = bonds => { const c = {}; for (const b of bonds) c[b[2]] = (c[b[2]] || 0) + 1; return JSON.stringify(Object.keys(c).sort().map(k => [k, c[k]])); };
const nOrder = (bonds, o) => bonds.filter(b => b[2] === o).length;
const touches = (b, i) => b[0] === i || b[1] === i;
const deepCopy = m => JSON.parse(JSON.stringify(m, (k, v) => (ArrayBuffer.isView(v) ? Array.from(v) : v)));
const same = (a, b) => JSON.stringify(deepCopy(a)) === JSON.stringify(deepCopy(b));
function signedVolume(x, c, a, b, d) {
  const p = i => [x[3 * i] - x[3 * c], x[3 * i + 1] - x[3 * c + 1], x[3 * i + 2] - x[3 * c + 2]];
  const u = p(a), v = p(b), w = p(d);
  return u[0] * (v[1] * w[2] - v[2] * w[1]) + u[1] * (v[2] * w[0] - v[0] * w[2]) + u[2] * (v[0] * w[1] - v[1] * w[0]);
}
/** Largest distance of atoms idx from their least-squares plane (Jacobi eigen-solver on the 3x3 covariance). */
function planeDeviation(x, idx) {
  const c = [0, 0, 0];
  idx.forEach(i => { for (let k = 0; k < 3; k++) c[k] += x[3 * i + k] / idx.length; });
  const A = [[0, 0, 0], [0, 0, 0], [0, 0, 0]], V = [[1, 0, 0], [0, 1, 0], [0, 0, 1]];
  idx.forEach(i => { const d = [0, 1, 2].map(k => x[3 * i + k] - c[k]); for (let a = 0; a < 3; a++) for (let b = 0; b < 3; b++) A[a][b] += d[a] * d[b]; });
  for (let sweep = 0; sweep < 30; sweep++) {
    for (let p = 0; p < 3; p++) for (let q = p + 1; q < 3; q++) {
      if (Math.abs(A[p][q]) < 1e-15) continue;
      const th = 0.5 * Math.atan2(2 * A[p][q], A[q][q] - A[p][p]), cs = Math.cos(th), sn = Math.sin(th);
      for (let k = 0; k < 3; k++) { const kp = A[k][p], kq = A[k][q]; A[k][p] = cs * kp - sn * kq; A[k][q] = sn * kp + cs * kq; }
      for (let k = 0; k < 3; k++) { const pk = A[p][k], qk = A[q][k]; A[p][k] = cs * pk - sn * qk; A[q][k] = sn * pk + cs * qk; }
      for (let k = 0; k < 3; k++) { const vp = V[k][p], vq = V[k][q]; V[k][p] = cs * vp - sn * vq; V[k][q] = sn * vp + cs * vq; }
    }
  }
  let mi = 0;
  for (let k = 1; k < 3; k++) if (A[k][k] < A[mi][mi]) mi = k;
  const nrm = [V[0][mi], V[1][mi], V[2][mi]];
  let mx = 0;
  idx.forEach(i => { const d = Math.abs([0, 1, 2].reduce((s, k) => s + (x[3 * i + k] - c[k]) * nrm[k], 0)); if (d > mx) mx = d; });
  return mx;
}
function neighbours(mol) {
  const adj = mol.Z.map(() => []);
  for (const b of mol.bonds) { adj[b[0]].push(b[1]); adj[b[1]].push(b[0]); }
  return adj;
}
/** Topological distance matrix (BFS), capped at 4. */
function topoDist(mol) {
  const adj = neighbours(mol), n = mol.Z.length;
  return mol.Z.map((_, s) => {
    const d = new Array(n).fill(4); d[s] = 0;
    const q = [s];
    while (q.length) { const u = q.shift(); if (d[u] >= 3) continue; for (const v of adj[u]) if (d[v] > d[u] + 1) { d[v] = d[u] + 1; q.push(v); } }
    return d;
  });
}
/** Every atom's bond-order sum plus implicit H must be an allowed valence (Kekulé check). */
function valencesOk(p) {
  const sums = p.atoms.map(a => a.hCount);
  for (const b of p.bonds) { sums[b[0]] += b[2]; sums[b[1]] += b[2]; }
  return p.atoms.every((a, i) => {
    if (C.isMetal(a.Z) || a.Z === 0 || p.bonds.some(b => touches(b, i) && C.isMetal(p.atoms[b[0] === i ? b[1] : b[0]].Z))) return true;
    const list = C._internal.valenceList(a.Z, a.charge);
    const smilesDefault = { 7: [3, 5], 15: [3, 5], 16: [2, 4, 6] }[a.Z];
    return list.includes(sums[i]) || (!a.bracket && smilesDefault && smilesDefault.includes(sums[i]));
  });
}

// ================================================================== parsing, kekulization, implicit H
const PARSE_CASES = [
  // [smiles, formula, double bonds, triple bonds]
  ['C', 'CH4', 0, 0], ['CCO', 'C2H6O', 0, 0], ['c1ccccc1', 'C6H6', 3, 0], ['Cc1ccccc1', 'C7H8', 3, 0],
  ['c1ccncc1', 'C5H5N', 3, 0], ['c1cc[nH]c1', 'C4H5N', 2, 0], ['c1ccoc1', 'C4H4O', 2, 0], ['c1ccsc1', 'C4H4S', 2, 0],
  ['c1cnc[nH]1', 'C3H4N2', 2, 0], ['c1ccc2[nH]ccc2c1', 'C8H7N', 4, 0], ['c1ccc2ccccc2c1', 'C10H8', 5, 0],
  ['Cn1cccc1', 'C5H7N', 2, 0], ['c1ncc2[nH]cnc2n1', 'C5H4N4', 4, 0],
  ['CN1C=NC2=C1C(=O)N(C(=O)N2C)C', 'C8H10N4O2', 4, 0], ['Cn1cnc2c1c(=O)n(C)c(=O)n2C', 'C8H10N4O2', 4, 0],
  ['CC(=O)Oc1ccccc1C(=O)O', 'C9H8O4', 5, 0], ['CC(=O)O', 'C2H4O2', 1, 0], ['[O-]C(=O)C', 'C2H3O2', 1, 0],
  ['[NH4+]', 'H4N', 0, 0], ['c1ccc(cc1)[N+](=O)[O-]', 'C6H5NO2', 4, 0], ['O=N(=O)c1ccccc1', 'C6H5NO2', 5, 0],
  ['Cl[Pt](Cl)([NH3])[NH3]', 'Cl2H6N2Pt', 0, 0], ['CC#N', 'C2H3N', 0, 1], ['O=c1cc[nH]cc1', 'C5H5NO', 3, 0],
  ['C%10CCCCC%10', 'C6H12', 0, 0], ['c%10ccc2c(c%10)ccc1ccccc12', 'C14H10', 7, 0], ['C%(123)CC%(123)', 'C3H6', 0, 0],
  ['[Na+].[Cl-]', 'ClNa', 0, 0], ['CCO.O', 'C2H8O2', 0, 0], ['OS(=O)(=O)O', 'H2O4S', 2, 0], ['CS(C)=O', 'C2H6OS', 1, 0],
  ['ClP(Cl)(Cl)(Cl)Cl', 'Cl5P', 0, 0], ['c1cc[cH-]c1', 'C5H5', 2, 0], ['C[n+]1ccccc1', 'C6H8N', 3, 0], ['[2H]C([2H])([2H])O', 'CH4O', 0, 0]
];

test('parse: formulas, implicit H and Kekulé forms', () => {
  for (const [s, formula, dbl, tpl] of PARSE_CASES) {
    const p = C.parseSmiles(s);
    eq(formulaOfParse(p), formula, s + ' formula');
    eq(nOrder(p.bonds, 2), dbl, s + ' double bonds');
    eq(nOrder(p.bonds, 3), tpl, s + ' triple bonds');
    ok(p.bonds.every(b => b[2] === 1 || b[2] === 2 || b[2] === 3), s + ' bond orders are Kekulé (1/2/3)');
    ok(valencesOk(p), s + ' every atom has a valid valence after kekulization');
    ok(Array.isArray(p.stereo.doubleBonds), s + ' stereo.doubleBonds present');
  }
});

test('parse: atom details', () => {
  const pyr = C.parseSmiles('c1cc[nH]c1');
  const n = pyr.atoms.findIndex(a => a.Z === 7);
  eq(pyr.atoms[n].hCount, 1, 'pyrrole [nH] keeps its H');
  ok(pyr.bonds.filter(b => touches(b, n)).every(b => b[2] === 1), 'pyrrole N has only single bonds');
  ok(pyr.atoms.every(a => a.aromatic), 'pyrrole atoms flagged aromatic');
  const nmp = C.parseSmiles('Cn1cccc1');
  ok(nmp.bonds.filter(b => touches(b, 1)).every(b => b[2] === 1), 'N-methylpyrrole n has no double bond');
  eq(nmp.atoms[1].hCount, 0, 'N-methylpyrrole n has no H');
  const py = C.parseSmiles('c1ccncc1');
  eq(py.bonds.filter(b => touches(b, 3)).map(b => b[2]).sort().join(','), '1,2', 'pyridine N has one double bond');
  const fur = C.parseSmiles('c1ccoc1');
  ok(fur.bonds.filter(b => touches(b, 3)).every(b => b[2] === 1), 'furan O single-bonded');
  const caf = C.parseSmiles('Cn1cnc2c1c(=O)n(C)c(=O)n2C');
  ok(caf.atoms.every((a, i) => a.Z !== 7 || caf.bonds.filter(b => touches(b, i)).reduce((s, b) => s + b[2], 0) === 3), 'caffeine N valences 3');
  const iso = C.parseSmiles('[13CH4]').atoms[0];
  eq(iso.isotope, 13, 'isotope'); eq(iso.hCount, 4, 'bracket H count'); eq(iso.Z, 6, 'bracket C');
  const fe = C.parseSmiles('[Fe+2]').atoms[0];
  eq(fe.Z, 26, 'Fe'); eq(fe.charge, 2, '+2 charge'); eq(fe.hCount, 0, 'bracket atom without H');
  eq(C.parseSmiles('[Cu++]').atoms[0].charge, 2, '++ charge');
  eq(C.parseSmiles('[O-2]').atoms[0].charge, -2, '-2 charge');
  eq(C.parseSmiles('[N---]').atoms[0].charge, -3, '--- charge');
  eq(C.parseSmiles('[NH4+:7]').atoms[0].atomClass, 7, 'atom class');
  eq(C.parseSmiles('[Pt]').atoms[0].Z, 78, '[Pt]');
  eq(C.parseSmiles('[Sc]').atoms[0].Z, 21, '[Sc] is scandium');
  eq(C.parseSmiles('[se]1cccc1').atoms[0].Z, 34, 'aromatic [se]');
  eq(C.parseSmiles('[CH3]').atoms[0].hCount, 3, 'radical keeps stated H');
  eq(C.parseSmiles('C[C@TH2H](N)O').atoms[1].chirality, '@@', '@TH2 read as @@');
  eq(C.parseSmiles('C[C@SP1H](N)O').atoms[1].chirality, null, 'non-tetrahedral class ignored');
  eq(C.parseSmiles('*C').atoms[0].Z, 0, 'wildcard atom');
  const cis = C.parseSmiles('Cl[Pt](Cl)([NH3])[NH3]');
  eq(cis.atoms[1].hCount, 0, 'Pt has no implicit H'); eq(cis.atoms[3].hCount, 3, '[NH3] has exactly 3 H');
  eq(cis.bonds.filter(b => touches(b, 1)).length, 4, 'Pt has four bonds');
  eq(C.parseSmiles('[Na+].[Cl-]').bonds.length, 0, 'dot-separated parts are not bonded');
  eq(C.parseSmiles('  CCO ethanol').atoms.length, 3, 'whitespace ends the SMILES');
  eq(C.parseSmiles('C=1CCCCC=1').bonds.filter(b => b[2] === 2).length, 1, 'ring-bond symbol on both sides');
  eq(C.parseSmiles('C1CCCCC=1').bonds.filter(b => b[2] === 2).length, 1, 'ring-bond symbol on closing side');
  eq(C.parseSmiles('c1ccccc1-c1ccccc1').bonds.filter(b => b[2] === 2).length, 6, 'biphenyl single link');
  eq(C.parseSmiles('c1ccccc1c1ccccc1').bonds.filter(b => b[2] === 2).length, 6, 'implicit link between rings is single');
});

test('parse: stereo', () => {
  const ala = C.parseSmiles('N[C@@H](C)C(=O)O');
  eq(ala.atoms[1].chirality, '@@', 'chirality @@');
  eq(JSON.stringify(ala.stereo.tetrahedral), JSON.stringify([{ atom: 1, chirality: '@@', neighbors: [0, -1, 2, 3] }]), 'neighbour order with implicit H');
  const ring = C.parseSmiles('C[C@@H]1CCCCN1');
  eq(JSON.stringify(ring.stereo.tetrahedral[0].neighbors), JSON.stringify([0, -1, 6, 2]), 'ring-closure neighbour takes the digit position');
  const first = C.parseSmiles('[C@@H](F)(Cl)Br');
  eq(JSON.stringify(first.stereo.tetrahedral[0].neighbors), JSON.stringify([-1, 1, 2, 3]), 'leading atom: implicit H first');
  const so = C.parseSmiles('C[S@](=O)CC');
  eq(JSON.stringify(so.stereo.tetrahedral[0].neighbors), JSON.stringify([0, -1, 2, 3]), 'lone pair takes the implicit-H slot');
  const cfg = s => C.parseSmiles(s).stereo.doubleBonds.map(d => d.config).join(',');
  eq(cfg('F/C=C/F'), 'trans', 'F/C=C/F'); eq(cfg('F/C=C\\F'), 'cis', 'F/C=C\\F');
  eq(cfg('C(\\F)=C/F'), 'trans', 'C(\\F)=C/F'); eq(cfg('C(/F)=C/F'), 'cis', 'C(/F)=C/F'); eq(cfg('F\\C=C\\F'), 'trans', 'F\\C=C\\F');
  eq(cfg('C/C=C/C=C/C'), 'trans,trans', 'conjugated diene'); eq(cfg('CC=CC'), '', 'no marks, no stereo');
  eq(cfg('F/C=C/1.F1'), 'trans', 'direction on a ring-closure bond');
});

test('parse: errors', () => {
  const bad = [['C1CC', 2], ['C(C', 2], ['C)C', 2], ['[Xx]', 2], ['C==C', 3], ['[C', 1], ['C%1', 2], ['CC(=)C', 4], ['1CC', 1], ['C.', 2],
    ['C1C1', 4], ['HC', 1], ['C()C', 3], ['.C', 1], ['CC..C', 4], ['C$', 2], ['[C@X]', 4], ['C[NH4+:]', 8], ['C11', 3], ['C=(O)C', 2], ['CC]', 3], ['C^C', 2]];
  for (const [s, pos] of bad) throwsMsg(() => C.parseSmiles(s), new RegExp('position ' + pos + '\\b'), 'syntax error in "' + s + '" at position ' + pos);
  throwsMsg(() => C.parseSmiles(''), /empty/i, 'empty SMILES');
  throwsMsg(() => C.parseSmiles('c1cccc1'), /kekul/i, 'five aromatic carbons cannot be kekulized');
  throwsMsg(() => C.parseSmiles('n1cccc1'), /kekul.*\[nH\]/i, 'pyrrole without [nH] hints at the fix');
  throwsMsg(() => C.parseSmiles('c1ccccc1c'), /kekul/i, 'acyclic aromatic atom');
});

// ================================================================== SMILES writer
test('toSmiles round trips', () => {
  const cases = PARSE_CASES.map(c => c[0]).concat(['C1CC2CCC1C2', '[H][H]', 'C#N', 'C1=CC=CC=C1C1CCCCC1', '[Fe+2].[O-]C(=O)C.[O-]C(=O)C',
    'C[S+](C)C', '[CH3]', 'F[P-](F)(F)(F)(F)F', 'C=C=C', 'CC(C)(C)c1ccc2ccccc2c1', 'C1CC2CC1C1CCCCC21', 'N#Cc1ccccc1C#N']);
  for (const s of cases) {
    const p = C.parseSmiles(s), out = C.toSmiles(p), back = C.parseSmiles(out);
    eq(formulaOfParse(back), formulaOfParse(p), s + ' -> ' + out + ' formula');
    eq(orderCounts(back.bonds), orderCounts(p.bonds), s + ' -> ' + out + ' bond orders');
    eq(back.atoms.reduce((q, a) => q + a.charge, 0), p.atoms.reduce((q, a) => q + a.charge, 0), s + ' total charge');
  }
  for (const s of ['CCO', 'CN1C=NC2=C1C(=O)N(C(=O)N2C)C', 'c1ccc(cc1)[N+](=O)[O-]', 'Cl[Pt](Cl)([NH3])[NH3]', '[2H]C([2H])([2H])O', '[NH4+].[Cl-]']) {
    const m = C.smilesTo3D(s).mol, out = C.toSmiles(m), back = C.parseSmiles(out);
    const plainH = i => m.Z[i] === 1 && !(m.isotopes && m.isotopes[i]);
    const heavy = m.bonds.filter(b => !plainH(b[0]) && !plainH(b[1]));
    eq(formulaOfParse(back), G.formula(m.Z), s + ' (explicit-H molecule) -> ' + out + ' formula');
    eq(orderCounts(back.bonds), orderCounts(heavy), s + ' (explicit-H molecule) bond orders');
  }
  eq(C.toSmiles({ Z: [8, 1, 1], bonds: [[0, 1, 1], [0, 2, 1]] }), 'O', 'hydrogens folded into the organic subset');
  eq(C.toSmiles({ Z: [6, 1, 1, 1], bonds: [[0, 1, 1], [0, 2, 1], [0, 3, 1]] }), '[CH3]', 'unusual H count bracketed');
  eq(C.toSmiles({ Z: [7, 1, 1, 1, 1], bonds: [[0, 1, 1], [0, 2, 1], [0, 3, 1], [0, 4, 1]], formalCharges: [1, 0, 0, 0, 0] }), '[NH4+]', 'charge bracketed');
  ok(/\[Pt\]/.test(C.toSmiles({ Z: [17, 78, 17], bonds: [[0, 1, 1], [1, 2, 1]] })), 'metal bracketed');
  eq(C.toSmiles({ Z: [11, 17], bonds: [], formalCharges: [1, -1] }), '[Na+].[Cl-]', 'components joined with dots');
  const big = C.parseSmiles(C.toSmiles({ Z: new Array(24).fill(6), bonds: Array.from({ length: 24 }, (_, i) => [i, (i + 1) % 24, 1]).concat(
    Array.from({ length: 11 }, (_, i) => [2 * i, 2 * i + 13 > 23 ? 2 * i + 13 - 24 : 2 * i + 13, 1])) }));
  eq(big.atoms.length, 24, 'many simultaneous ring closures (two-digit ring numbers) parse back');
  const arom = C.toSmiles({ Z: [6, 6, 6, 6, 6, 6], bonds: [[0, 1, 1.5], [1, 2, 1.5], [2, 3, 1.5], [3, 4, 1.5], [4, 5, 1.5], [5, 0, 1.5]], hCount: [1, 1, 1, 1, 1, 1] });
  eq(nOrder(C.parseSmiles(arom).bonds, 2), 3, 'transient 1.5 orders are written as a Kekulé structure');
  const inp = { Z: [6, 8, 1], bonds: [[0, 1, 2], [0, 2, 1]], formalCharges: [0, 0, 0] }, before = deepCopy(inp);
  C.toSmiles(inp);
  ok(same(inp, before), 'toSmiles does not mutate its input');
});

// ================================================================== 3D generation
const GEOM_SET = ['CCO', 'Cc1ccccc1', 'C1CCCCC1', 'c1ccc2ccccc2c1', 'CN1C=NC2=C1C(=O)N(C(=O)N2C)C', 'CC(=O)Oc1ccccc1C(=O)O',
  'N[C@@H](C)C(=O)O', 'C/C=C/C', 'C/C=C\\C', 'Cl[Pt](Cl)([NH3])[NH3]', 'c1ccc2[nH]ccc2c1', 'c1ccc(cc1)[N+](=O)[O-]', 'CC#N',
  'CC(C)Cc1ccc(cc1)C(C)C(=O)O', 'C1=CC=CC=C1', 'C[C@]12CC[C@H]3[C@@H](CC=C4C[C@@H](O)CC[C@]34C)[C@@H]1CC[C@@H]2[C@H](C)CCCC(C)C', 'C1CC2CCC1C2', '[O-]C(=O)C.[NH4+]'];

test('smilesTo3D: bond lengths, clashes, components', () => {
  for (const s of GEOM_SET) {
    const res = C.smilesTo3D(s, { seed: 11 });
    const m = res.mol, x = m.xyz, n = m.Z.length;
    ok(Array.isArray(res.warnings), s + ' warnings array');
    ok(x instanceof Float64Array && x.length === 3 * n && x.every(Number.isFinite), s + ' coordinates finite Float64Array(3N)');
    const p = C.parseSmiles(s);
    eq(G.formula(m.Z), formulaOfParse(p), s + ' explicit-H formula');
    const arom = new Set(res.aromaticBonds);
    let worst = 0;
    m.bonds.forEach((b, k) => {
      const ideal = C.idealBondLength(m.Z[b[0]], m.Z[b[1]], arom.has(k) ? 1.5 : b[2]);
      worst = Math.max(worst, Math.abs(G.dist(x, b[0], b[1]) - ideal));
    });
    ok(worst <= 0.08, s + ' bond lengths within 0.08 Å of ideal (worst ' + worst.toFixed(3) + ')');
    // non-bonded = separated by three or more bonds (or in different components)
    const td = topoDist(m);
    let minHeavy = Infinity, minH = Infinity;
    for (let i = 0; i < n; i++) for (let j = i + 1; j < n; j++) {
      const d = G.dist(x, i, j);
      if (m.Z[i] > 1 && m.Z[j] > 1 && td[i][j] >= 3) minHeavy = Math.min(minHeavy, d);
      if (m.Z[i] === 1 && m.Z[j] === 1) minH = Math.min(minH, d); // includes geminal pairs
    }
    ok(minHeavy >= 2.2, s + ' no non-bonded heavy-atom pair closer than 2.2 Å (min ' + minHeavy.toFixed(2) + ')');
    ok(minH >= 1.6, s + ' no H pair closer than 1.6 Å (min ' + minH.toFixed(2) + ')');
  }
  const salt = C.smilesTo3D('[Na+].[Cl-]').mol;
  ok(G.dist(salt.xyz, 0, 1) >= 3.0, 'ions placed at least 3 Å apart');
  const mix = C.smilesTo3D('CCO.O.[NH4+]').mol, fr = C.fragments(mix);
  eq(fr.length, 3, 'three components');
  let minSep = Infinity;
  for (let a = 0; a < fr.length; a++) for (let b = a + 1; b < fr.length; b++) for (const i of fr[a]) for (const j of fr[b]) minSep = Math.min(minSep, G.dist(mix.xyz, i, j));
  ok(minSep >= 3.0, 'components (including H) at least 3 Å apart: ' + minSep.toFixed(2));
  const a1 = C.smilesTo3D('CC(C)Cc1ccc(cc1)C(C)C(=O)O', { seed: 5 }).mol, a2 = C.smilesTo3D('CC(C)Cc1ccc(cc1)C(C)C(=O)O', { seed: 5 }).mol;
  ok(same(a1.xyz, a2.xyz), 'same seed gives identical coordinates');
  const hydr = C.smilesTo3D('CCO').mol;
  ok(hydr.Z.slice(0, 3).every(z => z > 1) && hydr.Z.slice(3).every(z => z === 1), 'SMILES atoms first, added H appended');
  throwsMsg(() => C.smilesTo3D('C1CC'), /position/, 'smilesTo3D reports parse errors');
});

test('smilesTo3D: planarity, rings, metals', () => {
  for (const s of ['c1ccccc1', 'c1ccc2ccccc2c1']) {
    const m = C.smilesTo3D(s, { seed: 3 }).mol;
    const all = m.Z.map((_, i) => i);
    ok(planeDeviation(m.xyz, all) <= 0.05, s + ' planar within 0.05 Å (' + planeDeviation(m.xyz, all).toFixed(4) + ')');
  }
  const caf = C.smilesTo3D('CN1C=NC2=C1C(=O)N(C(=O)N2C)C').mol;
  const ring = [1, 2, 3, 4, 5, 6, 7, 9, 10, 12];
  ok(planeDeviation(caf.xyz, ring) <= 0.05, 'caffeine ring system planar');
  const cy = C.smilesTo3D('C1CCCCC1', { seed: 3 }).mol;
  ok(planeDeviation(cy.xyz, [0, 1, 2, 3, 4, 5]) > 0.2, 'cyclohexane is not planar');
  const tors = [0, 1, 2, 3, 4, 5].map(k => G.dihedral(cy.xyz, k, (k + 1) % 6, (k + 2) % 6, (k + 3) % 6));
  ok(tors.every((t, k) => Math.abs(Math.abs(t) - 56) < 12 && Math.sign(t) !== Math.sign(tors[(k + 1) % 6])), 'cyclohexane is a chair: ' + tors.map(t => t.toFixed(0)).join(' '));
  const cp = C.smilesTo3D('Cl[Pt](Cl)([NH3])[NH3]').mol;
  near(G.angle(cp.xyz, 0, 1, 2), 90, 4, 'cis Cl-Pt-Cl'); near(G.angle(cp.xyz, 0, 1, 3), 180, 4, 'trans Cl-Pt-N');
  near(G.angle(cp.xyz, 3, 1, 4), 90, 4, 'cis N-Pt-N');
  near(G.dist(cp.xyz, 0, 1), C.idealBondLength(78, 17, 1), 0.03, 'Pt-Cl length');
  const oct = C.smilesTo3D('F[S](F)(F)(F)(F)F').mol;
  near(G.angle(oct.xyz, 0, 1, 6), 180, 5, 'SF6 trans F-S-F'); near(G.angle(oct.xyz, 0, 1, 2), 90, 5, 'SF6 cis F-S-F');
  const nit = C.smilesTo3D('CC#N').mol;
  near(G.angle(nit.xyz, 0, 1, 2), 180, 3, 'nitrile linear');
  const cyc = C.smilesTo3D('C1CC1').mol;
  near(G.angle(cyc.xyz, 0, 1, 2), 60, 3, 'cyclopropane 60°');
  const ind = C.smilesTo3D('c1ccc2[nH]ccc2c1').mol;
  ok(planeDeviation(ind.xyz, ind.Z.map((_, i) => i)) <= 0.05, 'indole planar including H');
});

test('smilesTo3D: stereochemistry', () => {
  // L-alanine (S): heavy neighbours N, CH3, COOH in SMILES order; '@@' with H second => V(N, CH3, COOH) < 0
  for (const seed of [1, 2, 3, 4]) {
    const L = C.smilesTo3D('N[C@@H](C)C(=O)O', { seed }).mol, D = C.smilesTo3D('N[C@H](C)C(=O)O', { seed }).mol;
    const vL = signedVolume(L.xyz, 1, 0, 2, 3), vD = signedVolume(D.xyz, 1, 0, 2, 3);
    ok(vL < -0.5 && vD > 0.5, 'L/D-alanine opposite signed volumes (seed ' + seed + '): ' + vL.toFixed(2) + ' / ' + vD.toFixed(2));
  }
  const L2 = C.smilesTo3D('C[C@@H](C(=O)O)N').mol; // the same L-alanine written differently
  ok(signedVolume(L2.xyz, 1, 5, 0, 2) < 0, 'L-alanine from a second spelling has the same handedness');
  // CIP check for L-alanine = (S): with priorities N > COOH > CH3, V(N, COOH, CH3) > 0 means counter-clockwise = S
  const L = C.smilesTo3D('N[C@@H](C)C(=O)O').mol;
  ok(signedVolume(L.xyz, 1, 0, 3, 2) > 0, 'L-alanine is (S)');
  const E = C.smilesTo3D('C/C=C/C').mol, Zb = C.smilesTo3D('C/C=C\\C').mol;
  ok(Math.abs(G.dihedral(E.xyz, 0, 1, 2, 3)) > 170, '(E)-2-butene trans: ' + G.dihedral(E.xyz, 0, 1, 2, 3).toFixed(1));
  ok(Math.abs(G.dihedral(Zb.xyz, 0, 1, 2, 3)) < 10, '(Z)-2-butene cis: ' + G.dihedral(Zb.xyz, 0, 1, 2, 3).toFixed(1));
  const st = C.smilesTo3D('F/C=C(/Cl)\\Br').mol; // F trans to Cl
  ok(Math.abs(G.dihedral(st.xyz, 0, 1, 2, 3)) > 170 && Math.abs(G.dihedral(st.xyz, 0, 1, 2, 4)) < 10, 'trisubstituted alkene F/Cl trans, F/Br cis');
  // a stereocentre written with a ring-closure digit: four equivalent spellings of one enantiomer of 2-methylpiperidine
  const forms = [['C[C@@H]1CCCCN1', 1, 6, 2, 0], ['N1CCCC[C@H]1C', 5, 0, 4, 6], ['[C@@H]1(C)CCCCN1', 0, 6, 2, 1], ['C[C@H]1NCCCC1', 1, 2, 6, 0]];
  const signs = forms.map(([s, c, nN, r, me]) => Math.sign(signedVolume(C.smilesTo3D(s).mol.xyz, c, nN, r, me)));
  ok(signs.every(v => v === signs[0]), 'ring-closure chirality consistent across spellings: ' + signs.join(' '));
  const camphor = C.smilesTo3D('CC1(C)[C@@H]2CC[C@@]1(C)C(=O)C2');
  eq(camphor.warnings.length, 0, 'camphor embeds without strain warnings');
});

test('smilesTo3D: BM.ff hand-off', () => {
  let calls = 0;
  BM.ff = { optimize(mol, opts) { calls++; ok(mol.bonds.some(b => b[2] === 1.5), 'aromatic bonds passed to the force field as 1.5'); ok(opts && opts.maxSteps > 0, 'maxSteps passed'); return { xyz: Float64Array.from(mol.xyz), energy: 0, steps: 1, converged: true }; } };
  const r1 = C.smilesTo3D('Cc1ccccc1');
  ok(calls === 1 && r1.warnings.length === 0, 'force field used when present');
  BM.ff = { optimize(mol) { return { xyz: new Float64Array(mol.xyz.length).fill(NaN), energy: 0, steps: 1, converged: false }; } };
  const r2 = C.smilesTo3D('N[C@@H](C)C(=O)O');
  ok(r2.mol.xyz.every(Number.isFinite) && r2.warnings.some(w => /rejected/.test(w)), 'bad force-field result rejected');
  BM.ff = { optimize() { throw new Error('boom'); } };
  const r3 = C.smilesTo3D('CCO');
  ok(r3.mol.xyz.every(Number.isFinite) && r3.warnings.some(w => /boom/.test(w)), 'force-field exception reported as a warning');
  delete BM.ff;
});

// ================================================================== bond-order perception
test('perceiveBondOrders', () => {
  const cases = ['c1ccccc1', 'CN1C=NC2=C1C(=O)N(C(=O)N2C)C', 'CC(=O)O', 'CC#N', 'c1ccncc1', 'c1ccc2ccccc2c1', 'CC(=O)Oc1ccccc1C(=O)O',
    '[O-]C(=O)C', 'c1ccc(cc1)[N+](=O)[O-]', 'O=C=O', 'C=C=C', 'CS(=O)(=O)C', 'OP(=O)(O)O', 'CN=[N+]=[N-]', 'c1cc[nH]c1', 'CC(=O)Nc1ccc(O)cc1',
    'Cl[Pt](Cl)([NH3])[NH3]', 'C#CC#C', '[O-][n+]1ccccc1'];
  for (const s of cases) {
    const m = C.smilesTo3D(s).mol;
    const ring = C._internal.ringBondFlags(m.Z.length, m.bonds);
    const reset = { Z: m.Z, xyz: m.xyz, bonds: m.bonds.map(b => [b[0], b[1], 1]) }, before = deepCopy(reset);
    for (const withQ of [false, true]) {
      const input = withQ ? Object.assign({}, reset, { formalCharges: m.formalCharges }) : reset;
      const got = C.perceiveBondOrders(input);
      const tag = s + (withQ ? ' (with charges)' : '');
      eq(got.length, m.bonds.length, tag + ' bond count');
      ok(got.every((b, k) => b[0] === m.bonds[k][0] && b[1] === m.bonds[k][1]), tag + ' keeps bond order/indices');
      ok(got.every((b, k) => ring[k] || b[2] === m.bonds[k][2]), tag + ' acyclic bonds match the original orders');
      const sums = bonds => { const v = new Array(m.Z.length).fill(0); bonds.forEach(b => { v[b[0]] += b[2]; v[b[1]] += b[2]; }); return v.join(','); };
      eq(sums(got), sums(m.bonds), tag + ' every atom regains its original valence (valid Kekulé structure)');
      eq(orderCounts(got), orderCounts(m.bonds), tag + ' same multiset of bond orders');
    }
    ok(same(reset, before), s + ' input not mutated');
  }
  const cp = C.smilesTo3D('Cl[Pt](Cl)([NH3])[NH3]').mol;
  ok(C.perceiveBondOrders(cp).filter(b => b[0] === 1 || b[1] === 1).every(b => b[2] === 1), 'bonds to Pt stay single');
  // an idealised QM-like benzene (C-C 1.397, C-H 1.084) with connectivity only from coordinates
  const Z = [], xyz = [];
  for (let k = 0; k < 6; k++) { const a = k * Math.PI / 3; Z.push(6); xyz.push(1.397 * Math.cos(a), 1.397 * Math.sin(a), 0); }
  for (let k = 0; k < 6; k++) { const a = k * Math.PI / 3; Z.push(1); xyz.push(2.481 * Math.cos(a), 2.481 * Math.sin(a), 0); }
  const bz = C.perceiveBondOrders({ Z, xyz: Float64Array.from(xyz) });
  eq(bz.length, 12, 'benzene connectivity from coordinates');
  eq(nOrder(bz, 2), 3, 'benzene from coordinates: three double bonds');
  ok([0, 1, 2, 3, 4, 5].every(i => bz.filter(b => touches(b, i) && b[2] === 2).length === 1), 'each benzene carbon has one double bond');
  const nitro = C.perceiveBondOrders({ Z: C.smilesTo3D('C[N+](=O)[O-]').mol.Z, xyz: C.smilesTo3D('C[N+](=O)[O-]').mol.xyz }, { withCharges: true });
  eq(nOrder(nitro.bonds, 2), 1, 'nitromethane without charges: one N=O');
  eq(nitro.formalCharges.filter(q => q).sort().join(','), '-1,1', 'implied charge separation N+ / O-');
});

// ================================================================== valence and hydrogens
test('valenceInfo', () => {
  const one = (Z, q) => C.valenceInfo({ Z: [Z], xyz: new Float64Array(3), bonds: [], formalCharges: [q] }, 0);
  eq(one(6, 0).missingH, 4, 'C needs 4 H'); eq(one(7, 0).missingH, 3, 'N needs 3 H'); eq(one(8, 0).missingH, 2, 'O needs 2 H');
  eq(one(7, 1).missingH, 4, 'N+ target 4'); eq(one(8, 1).missingH, 3, 'O+ target 3'); eq(one(8, -1).missingH, 1, 'O- target 1');
  eq(one(6, -1).missingH, 3, 'C- target 3'); eq(one(6, 1).missingH, 3, 'C+ target 3'); eq(one(5, -1).missingH, 4, 'B- target 4');
  eq(one(78, 0).missingH, 0, 'no H on Pt'); eq(one(78, 0).hybridization, 'metal', 'Pt hybridization');
  eq(one(18, 0).missingH, 0, 'no H on Ar'); eq(one(1, 0).missingH, 0, 'no H on H');
  const ac = C.smilesTo3D('CC(=O)O').mol;
  const v1 = C.valenceInfo(ac, 1);
  eq(v1.element, 'C', 'element'); eq(v1.bondOrderSum, 4, 'carboxyl C bond-order sum'); eq(v1.targetValence, 4, 'target'); eq(v1.hybridization, 'sp2', 'carboxyl C sp2');
  eq(C.valenceInfo(ac, 0).hybridization, 'sp3', 'methyl sp3'); eq(C.valenceInfo(ac, 2).hybridization, 'sp2', 'carbonyl O sp2');
  eq(C.valenceInfo(C.smilesTo3D('CC#N').mol, 1).hybridization, 'sp', 'nitrile C sp');
  eq(C.valenceInfo(C.smilesTo3D('CC(N)=O').mol, 2).hybridization, 'sp2', 'amide N conjugated sp2');
  eq(C.valenceInfo(C.smilesTo3D('CCN').mol, 2).hybridization, 'sp3', 'amine N sp3');
  const so2 = C.valenceInfo(C.smilesTo3D('CS(=O)(=O)C').mol, 1);
  eq(so2.missingH, 0, 'sulfone S saturated'); eq(so2.hybridization, 'sp3', 'sulfone S tetrahedral');
  const nitro = C.smilesTo3D('O=N(=O)C').mol;
  eq(C.valenceInfo(nitro, 1).missingH, 0, 'pentavalent nitro N needs no H');
  const quat = { Z: [7, 6, 6, 6, 6], xyz: new Float64Array(15), bonds: [[0, 1, 1], [0, 2, 1], [0, 3, 1], [0, 4, 1]] };
  eq(C.valenceInfo(quat, 0).missingH, 0, 'quaternary N without charge gets no H');
  throwsMsg(() => C.valenceInfo(ac, 99), /range/, 'index check');
});

test('add / remove / adjust hydrogens', () => {
  const caf = C.smilesTo3D('CN1C=NC2=C1C(=O)N(C(=O)N2C)C').mol, before = deepCopy(caf);
  const stripped = C.removeHydrogens(caf);
  eq(G.formula(stripped.mol.Z), 'C8N4O2', 'caffeine stripped');
  eq(stripped.map.length, caf.Z.length, 'map covers old atoms');
  ok(Array.from(stripped.map).every((v, i) => (caf.Z[i] === 1 ? v === -1 : v >= 0)), 'map marks deleted H with -1');
  const back = C.addHydrogens(stripped.mol);
  eq(G.formula(back.mol.Z), G.formula(caf.Z), 'caffeine H restored');
  eq(back.added.length, 10, 'ten H added');
  ok(same(caf, before), 'inputs not mutated');
  const bm = back.mol, adj = neighbours(bm);
  const ch = [], hch = [];
  bm.bonds.forEach(b => { if (bm.Z[b[0]] === 6 && bm.Z[b[1]] === 1) ch.push(G.dist(bm.xyz, b[0], b[1])); });
  ok(ch.length === 10 && ch.every(d => Math.abs(d - 1.09) < 0.005), 'C-H = 1.09 Å');
  adj.forEach((nb, i) => { const hs = nb.filter(j => bm.Z[j] === 1); if (bm.Z[i] === 6 && hs.length === 3) for (let a = 0; a < 3; a++) for (let b = a + 1; b < 3; b++) hch.push(G.angle(bm.xyz, hs[a], i, hs[b])); });
  ok(hch.length === 9 && hch.every(a => Math.abs(a - 109.47) < 0.5), 'methyl H-C-H = 109.5°');
  const ringH = adj[2].find(j => bm.Z[j] === 1);
  ok(planeDeviation(bm.xyz, [1, 2, 3, 4, 5, ringH]) < 0.05, 'aromatic C-H added in the ring plane');
  // staggering: re-added ethane hydrogens
  const eth = C.addHydrogens(C.removeHydrogens(C.smilesTo3D('CC').mol).mol).mol, ea = neighbours(eth);
  const d = G.dihedral(eth.xyz, ea[0].find(j => eth.Z[j] === 1), 0, 1, ea[1].find(j => eth.Z[j] === 1));
  ok([60, 180].some(t => Math.abs(Math.abs(d) - t) < 5), 'ethane H staggered: ' + d.toFixed(1));
  // conjugated sp2: amide NH2 and vinyl CH2 in plane
  const acm = C.addHydrogens(C.removeHydrogens(C.smilesTo3D('CC(N)=O').mol).mol).mol, am = neighbours(acm);
  ok(planeDeviation(acm.xyz, [0, 1, 2, 3].concat(am[2].filter(j => acm.Z[j] === 1))) < 0.05, 'amide NH2 planar with the carbonyl');
  // only selected atoms
  const ethanol = C.removeHydrogens(C.smilesTo3D('CCO').mol).mol;
  const partial = C.addHydrogens(ethanol, { atoms: [2] });
  eq(partial.added.length, 1, 'H added only to the selected O'); eq(G.formula(partial.mol.Z), 'C2HO', 'partial formula');
  near(G.dist(partial.mol.xyz, 2, partial.added[0]), 0.96, 1e-9, 'O-H length');
  eq(C.addHydrogens({ Z: [78, 17, 17], xyz: Float64Array.from([0, 0, 0, 2.3, 0, 0, 0, 2.3, 0]), bonds: [[0, 1, 1], [0, 2, 1]] }).added.length, 0, 'no H on metals or terminal halides');
  eq(C.addHydrogens({ Z: [54], xyz: new Float64Array(3), bonds: [] }).added.length, 0, 'no H on noble gases');
  const nh4 = C.addHydrogens({ Z: [7], xyz: new Float64Array(3), bonds: [], formalCharges: [1] }).mol;
  eq(G.formula(nh4.Z), 'H4N', 'N+ gets four H'); eq(nh4.formalCharges.length, 5, 'formalCharges extended');
  near(G.angle(nh4.xyz, 1, 0, 2), 109.47, 0.01, 'tetrahedral NH4+');
  // remove selected
  const etoh = C.smilesTo3D('CCO').mol;
  const rmO = C.removeHydrogens(etoh, { atoms: [2] });
  eq(G.formula(rmO.mol.Z), 'C2H5O', 'removeHydrogens on a heavy atom removes its H');
  const firstH = etoh.Z.indexOf(1);
  eq(C.removeHydrogens(etoh, { atoms: [firstH] }).mol.Z.length, etoh.Z.length - 1, 'removeHydrogens on an H removes it');
  // adjust after an edit: ethane C-C made double, then back
  const ethane = C.smilesTo3D('CC').mol;
  const ethene = C.adjustHydrogens(C.setBond(ethane, 0, 1, 2));
  eq(G.formula(ethene.mol.Z), 'C2H4', 'adjust: ethane -> ethene');
  eq(ethene.map.length, ethane.Z.length, 'adjust map length');
  const en = ethene.mol, ena = neighbours(en), h0 = ena[0].filter(j => en.Z[j] === 1);
  near(G.angle(en.xyz, h0[0], 0, h0[1]), 120, 0.5, 'sp2 H-C-H');
  ok(planeDeviation(en.xyz, en.Z.map((_, i) => i)) < 0.05, 'ethene H in the plane of the double bond');
  const untouched = C.adjustHydrogens(ethane);
  eq(untouched.added.length, 0, 'adjust leaves correct atoms alone'); ok(same(untouched.mol.xyz, ethane.xyz), 'coordinates unchanged');
  const charged = C.adjustHydrogens(Object.assign({}, C.removeHydrogens(C.smilesTo3D('C[O-]').mol).mol));
  eq(G.formula(charged.mol.Z), 'CH3O', 'formal charge respected (methoxide keeps O without H)');
});

// ================================================================== editing helpers
test('editing helpers', () => {
  near(C.idealBondLength(6, 6, 1), 1.52, 1e-6, 'C-C'); near(C.idealBondLength(6, 6, 2), 1.33, 1e-6, 'C=C');
  near(C.idealBondLength(6, 6, 3), 1.2008, 1e-6, 'C#C'); near(C.idealBondLength(6, 6, 1.5), 1.3908, 1e-6, 'aromatic C:C');
  near(C.idealBondLength(6, 1, 1), 1.09, 1e-9, 'C-H'); near(C.idealBondLength(1, 8, 1), 0.96, 1e-9, 'O-H');
  near(C.idealBondLength(1, 1, 1), 0.74, 1e-9, 'H-H'); near(C.idealBondLength(17, 1, 1), 1.33, 1e-6, 'Cl-H uses rc + 0.31');
  near(C.idealBondLength(78, 17, 1), 2.38, 1e-6, 'Pt-Cl');
  const lone = { Z: [6], xyz: Float64Array.from([1, 2, 3]), bonds: [] };
  const p0 = C.newAtomPosition(lone, 0, 8, 2);
  near(Math.hypot(p0[0] - 1, p0[1] - 2, p0[2] - 3), C.idealBondLength(6, 8, 2), 1e-9, 'isolated atom: ideal distance');
  ok(same(p0, C.newAtomPosition(lone, 0, 8, 2)), 'isolated atom: deterministic direction');
  const ethane = C.removeHydrogens(C.smilesTo3D('CC').mol).mol;
  const p1 = C.newAtomPosition(ethane, 1, 6, 1), prop = C.addAtom(ethane, 6, p1);
  eq(prop.index, 2, 'addAtom index'); near(G.angle(prop.mol.xyz, 0, 1, 2), 109.47, 0.1, 'new atom at the tetrahedral site');
  near(G.dist(prop.mol.xyz, 1, 2), 1.52, 1e-6, 'new atom at the ideal length');
  const ethyl = C.setBond(prop.mol, 1, 2, 1);
  const p2 = C.newAtomPosition(ethyl, 1, 8, 1), iso = C.addAtom(ethyl, 8, p2).mol;
  near(G.angle(iso.xyz, 0, 1, 3), 109.47, 1, 'third substituent tetrahedral'); near(G.angle(iso.xyz, 2, 1, 3), 109.47, 1, 'third substituent tetrahedral (2)');
  const ald = C.removeHydrogens(C.smilesTo3D('CC=O').mol).mol;
  const p3 = C.newAtomPosition(ald, 1, 1, 1), withH = C.addAtom(ald, 1, p3).mol;
  near(G.angle(withH.xyz, 0, 1, 3), 180 - G.angle(withH.xyz, 0, 1, 2) / 2, 0.5, 'sp2 site bisects the exterior angle'); ok(planeDeviation(withH.xyz, [0, 1, 2, 3]) < 0.02, 'sp2 site in plane');
  const pt = C.removeHydrogens(C.smilesTo3D('Cl[Pt](Cl)[NH3]').mol).mol;
  const p4 = C.newAtomPosition(pt, 1, 7, 1), pt4 = C.addAtom(pt, 7, p4).mol;
  const angs = [0, 2, 3].map(j => G.angle(pt4.xyz, j, 1, 4)).sort((a, b) => a - b);
  ok(Math.abs(angs[0] - 90) < 6 && Math.abs(angs[1] - 90) < 6 && Math.abs(angs[2] - 180) < 6, 'fourth Pt ligand completes the square: ' + angs.map(a => a.toFixed(0)).join('/'));
  const m = C.smilesTo3D('CCO').mol, before = deepCopy(m);
  const del = C.deleteAtoms(m, [1]);
  eq(del.mol.Z.length, m.Z.length - 1, 'deleteAtoms count'); eq(del.map[1], -1, 'deleted atom maps to -1'); eq(del.map[2], 1, 'later atoms shift');
  ok(del.mol.bonds.every(b => b[0] < del.mol.Z.length && b[1] < del.mol.Z.length), 'bonds remapped');
  eq(del.mol.bonds.length, m.bonds.length - 4, 'bonds to the deleted atom removed');
  ok(del.map instanceof Int32Array, 'map is Int32Array');
  const sb = C.setBond(m, 0, 2, 1);
  eq(sb.bonds.length, m.bonds.length + 1, 'setBond adds'); eq(C.setBond(sb, 2, 0, 0).bonds.length, m.bonds.length, 'order 0 deletes');
  eq(C.setBond(m, 0, 1, 2).bonds.find(b => b[0] === 0 && b[1] === 1)[2], 2, 'setBond changes order');
  ok(same(m, before), 'editing helpers do not mutate input');
  throwsMsg(() => C.setBond(m, 0, 0, 1), /itself/, 'self bond rejected');
  eq(JSON.stringify(C.fragments({ Z: [6, 6, 8, 11], bonds: [[0, 1, 1]] })), '[[0,1],[2],[3]]', 'fragments');
  const named = C.addAtom({ Z: [6], xyz: new Float64Array(3), bonds: [], atomNames: ['C1'] }, 1, [1, 0, 0]).mol;
  eq(named.atomNames.length, 2, 'per-atom arrays extended');
});

// ================================================================== performance
test('performance', () => {
  const peptide = 'CC(C)C[C@H](NC(=O)[C@H](Cc1ccccc1)NC(=O)[C@@H](N)CO)C(=O)N[C@@H](Cc1c[nH]c2ccccc12)C(=O)N[C@@H](CCCCN)C(=O)N[C@@H](C)C(=O)NCC(=O)O';
  const steroid = 'C[C@]12CC[C@H]3[C@@H](CC=C4C[C@@H](O)CC[C@]34C)[C@@H]1CC[C@@H]2[C@H](C)CCCC(C)C';
  const heavy = C.parseSmiles(peptide).atoms.length;
  ok(heavy >= 55, 'peptide has about 60 heavy atoms (' + heavy + ')');
  C.smilesTo3D('CCO'); // warm-up
  const times = [];
  for (const seed of [1, 2, 3]) { const t0 = process.hrtime.bigint(); C.smilesTo3D(peptide, { seed }); times.push(Number(process.hrtime.bigint() - t0) / 1e6); }
  times.sort((a, b) => a - b);
  ok(times[1] <= 300, 'median 3D time for ' + heavy + ' heavy atoms <= 0.3 s: ' + times.map(t => t.toFixed(0)).join('/') + ' ms');
  console.log('  smilesTo3D, ' + heavy + ' heavy atoms: ' + times.map(t => t.toFixed(0)).join(' / ') + ' ms');
  const st = C.smilesTo3D(steroid, { seed: 1 });
  eq(st.warnings.length, 0, 'cholesterol embeds cleanly');
  const ch = st.mol, sv = C._internal; // ring B double bond present and stereo kept
  ok(ch.bonds.some(b => b[2] === 2), 'cholesterol C5=C6 kept');
  ok(sv && true, 'internal helpers exposed');
  const t1 = process.hrtime.bigint();
  for (let k = 0; k < 20; k++) C.parseSmiles(peptide);
  ok(Number(process.hrtime.bigint() - t1) / 1e6 < 200, 'parsing is fast');
});

// ------------------------------------------------------------------ summary
console.log('\nBlueMolecule chemistry tests: ' + pass + ' passed, ' + fail + ' failed');
if (fail) {
  for (const f of failures) console.log('  FAIL ' + f);
  process.exit(1);
}
