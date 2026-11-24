#!/usr/bin/env node
/* BlueMolecule tests for orbitals and densities computed from fchk and Molden files (src/orbitals.js, the basis-set
 * readers in src/parsers.js) and for molecular surfaces and the electrostatic potential (src/molsurf.js).
 *
 * The references are independent of the JavaScript code: tools/mini_scf.py (its own integrals and solid harmonics)
 * wrote the fixtures in tests/fixtures/orbitals/ with, for each file, the basis and orbitals a reader should return,
 * orbital values, densities and spin densities at 24 points, and the overlap matrix. The Molden files include the
 * conventions of other programs (ORCA's normalisation and signs, Cartesian functions normalised like x^l).
 * Run:  node tests/test_orbitals.js */
'use strict';
const fs = require('fs');
const path = require('path');
const ROOT = path.join(__dirname, '..');
for (const m of ['elements', 'bigfile', 'parsers', 'condense', 'readers-struct', 'writers', 'readers-qc', 'readers-qc2', 'geom', 'chem', 'analysis', 'isosurface', 'orbitals', 'molsurf']) {
  require(path.join(ROOT, 'src', m + '.js'));
}
const BM = globalThis.BM, O = BM.orbitals, S = BM.molsurf;
const DIR = path.join(__dirname, 'fixtures', 'orbitals');
const read = f => fs.readFileSync(path.join(DIR, f), 'utf8');

let pass = 0, fail = 0, current = '';
const failures = [];
function ok(c, msg) { if (c) pass++; else { fail++; failures.push('[' + current + '] ' + msg); } }
function eq(a, b, msg) { ok(a === b, msg + ': got ' + JSON.stringify(a) + ', expected ' + JSON.stringify(b)); }
function near(a, b, tol, msg) { ok(typeof a === 'number' && Math.abs(a - b) <= tol, msg + ': got ' + a + ', expected ' + b); }
function maxDiff(a, b) { let m = 0; for (let i = 0; i < b.length; i++) m = Math.max(m, Math.abs(a[i] - b[i])); return m; }
function test(name, fn) { current = name; try { fn(); } catch (e) { fail++; failures.push('[' + name + '] threw: ' + (e && e.stack || e)); } }
async function testAsync(name, fn) { current = name; try { await fn(); } catch (e) { fail++; failures.push('[' + name + '] threw: ' + (e && e.stack || e)); } }

const FILES = fs.readdirSync(DIR).filter(f => /\.(fchk|molden)$/.test(f)).sort();
const loaded = {};
function load(f) {
  if (!loaded[f]) {
    const ds = BM.parseFile(f, read(f)), J = JSON.parse(read(f + '.json')), M = O.model(ds);
    loaded[f] = { ds, J, M };
  }
  return loaded[f];
}
// printed precision of the coefficients: ORCA-style files have six decimals
const tolFor = f => (/_orca/.test(f) ? 3e-6 : 1e-7);

// ------------------------------------------------------------------ readers
function same(a, b, p, out) {
  if (b === null || typeof b !== 'object') {
    if (typeof b === 'number' ? !(typeof a === 'number' && Math.abs(a - b) <= 1e-12 * Math.max(1, Math.abs(b))) : a !== b) out.push(p + ': ' + a + ' vs ' + b);
    return;
  }
  if (Array.isArray(b)) {
    if (!a || a.length !== b.length) { out.push(p + ': length ' + (a && a.length) + ' vs ' + b.length); return; }
    for (let i = 0; i < b.length && out.length < 5; i++) same(a[i], b[i], p + '[' + i + ']', out);
    return;
  }
  for (const k of Object.keys(b)) same(a ? a[k] : undefined, b[k], p + '.' + k, out);
}
for (const f of FILES) {
  test('reader: ' + f, () => {
    const { ds, J } = load(f), out = [];
    ok(ds.basis && ds.mos, 'basis and orbitals read');
    const got = { basis: ds.basis, mos: { restricted: ds.mos.restricted,
      alpha: Object.assign({}, ds.mos.alpha, { coeffs: Array.from(ds.mos.alpha.coeffs) }),
      beta: ds.mos.beta && Object.assign({}, ds.mos.beta, { coeffs: Array.from(ds.mos.beta.coeffs) }) } };
    same(got, J.ds, 'ds', out);
    ok(!out.length, 'the basis set and orbitals the fixture expects: ' + out.join('; '));
    ok(!ds.warnings.some(w => /orbitals/i.test(w)), 'no orbital warnings: ' + ds.warnings.join(' | '));
    ok(ds.orbitals && ds.orbitals.alpha.energies.length === ds.mos.alpha.energies.length, 'orbital energies for the level plot');
  });
}

// ------------------------------------------------------------------ conventions, values, densities, overlap
for (const f of FILES) {
  test('orbitals: ' + f, () => {
    const { J, M } = load(f), E = J.expect;
    ok(M, 'model built');
    eq(M.nbf, E.nbf, 'basis functions');
    const det = O.detect(M);
    ok(det.error < 1e-5, 'orbitals orthonormal (deviation ' + det.error.toExponential(1) + ')');
    for (const k of ['prim', 'cart', 'sign']) {
      if (E.convention[k] && det.relevant[k]) eq(det.conv[k], E.convention[k], 'convention ' + k);
    }
    const pts = E.points_bohr.flat(), tol = tolFor(f);
    for (const spin of Object.keys(E.values)) {
      for (const k of Object.keys(E.values[spin])) {
        const v = O.values(M, { kind: 'mo', spin, index: +k }, pts);
        ok(maxDiff(v, E.values[spin][k]) < tol, spin + ' MO ' + (+k + 1) + ' at 24 points: max error ' + maxDiff(v, E.values[spin][k]).toExponential(2));
      }
    }
    const rho = O.values(M, { kind: 'density' }, pts);
    ok(maxDiff(rho, E.density) < 4 * tol * Math.max(1, ...E.density), 'electron density at 24 points: max error ' + maxDiff(rho, E.density).toExponential(2));
    if (E.spin_density) {
      ok(O.hasSpin(M), 'spin density offered');
      const sd = O.values(M, { kind: 'spin' }, pts);
      ok(maxDiff(sd, E.spin_density) < 4 * tol, 'spin density at 24 points: max error ' + maxDiff(sd, E.spin_density).toExponential(2));
    } else ok(!O.hasSpin(M), 'no spin density for a closed shell');
    const Sm = O.overlap(M, det.conv.prim), n = M.nbf, ref = E.overlap;
    // the fixture's overlap is in Gaussian's order; Molden lists Cartesian g functions in another order
    const perm = [];
    for (const sh of M.shells) {
      if (M.source === 'molden' && !sh.pure && sh.l >= 4) { const g = O.cartOrder(sh.l, 'fchk').map(c => c.join()); O.cartOrder(sh.l, 'molden').forEach(c => perm.push(sh.off + g.indexOf(c.join()))); }
      else for (let q = 0; q < sh.nf; q++) perm.push(sh.off + q);
    }
    const R = (i, j) => { const a = Math.max(perm[i], perm[j]), b = Math.min(perm[i], perm[j]); return ref[a * (a + 1) / 2 + b]; };
    let worst = 0;
    for (let i = 0; i < n; i++) for (let j = 0; j < n; j++) worst = Math.max(worst, Math.abs(Sm[i * n + j] - R(i, j)));
    // the Molden files print exponents and coefficients with fewer digits than the fixture's own basis
    ok(worst < (f.endsWith('.fchk') ? 1e-7 : /_orca/.test(f) ? 1e-6 : 1e-7), 'overlap matrix: max error ' + worst.toExponential(2));
  });
}

test('the conventions really differ: a wrong one is not orthonormal', () => {
  const right = f => { const c = load(f).M.detected.conv; return c; };
  for (const [f, wrong] of [['fg_cart_axis.molden', { cart: 'component' }], ['fg_pure_orca.molden', { sign: 'standard' }], ['h2o_631gd_pure_orca.molden', { prim: 'normalized' }],
    ['fg_cart.molden', { cart: 'axis' }], ['fg_cart.molden', { cart: 'radial' }], ['fg_cart.molden', { cart: 'scaled' }], ['fg_cart.fchk', { prim: 'raw' }]]) {
    const { M } = load(f);
    O.detect(M);
    const e = O.check(M, Object.assign({}, right(f), wrong));
    ok(e > 1e-2, f + ': ' + JSON.stringify(wrong) + ' gives ' + e.toExponential(1));
  }
  const { M } = load('fg_pure_orca.molden');
  ok(/ORCA/.test(O.describe(M.detected)) && /primitive normalisation/.test(O.describe(M.detected)), 'described: ' + O.describe(M.detected));
});

// ------------------------------------------------------------------ the pieces
test('solid harmonics: every pure shell orthonormal (l = 2 to 6), Gaussian signs', () => {
  for (let l = 2; l <= 6; l++) {
    const ds = { basis: { source: 'fchk', nbf: 2 * l + 1, shells: [{ atom: 0, center: [0.1, -0.2, 0.3], l, pure: true, exps: [0.9, 0.3], coefs: [0.6, 0.5] }] },
      mos: { restricted: true, alpha: { energies: [], occupations: [], symmetries: null, coeffs: new Float64Array(0) }, beta: null } };
    ds.mos.alpha.coeffs = new Float64Array((2 * l + 1) * (2 * l + 1));
    const M = O.model(ds), Sm = O.overlap(M, 'normalized'), n = 2 * l + 1;
    let worst = 0;
    for (let i = 0; i < n; i++) for (let j = 0; j < n; j++) worst = Math.max(worst, Math.abs(Sm[i * n + j] - (i === j ? 1 : 0)));
    ok(worst < 1e-12, 'l = ' + l + ': |S - 1| = ' + worst.toExponential(1));
  }
  // signs and shapes: d(+2) = x^2 - y^2, f(-3) = y(3x^2 - y^2)
  const d2 = O.solidHarmonic(2, 2), f3 = O.solidHarmonic(3, -3);
  ok(d2.get('2,0,0') > 0 && d2.get('0,2,0') < 0 && d2.size === 2, 'd+2 is x² − y²');
  ok(f3.get('2,1,0') > 0 && f3.get('0,3,0') < 0 && Math.abs(f3.get('2,1,0') / f3.get('0,3,0') + 3) < 1e-12, 'f−3 is 3x²y − y³');
  eq(O.mOrder(3).join(','), '0,1,-1,2,-2,3,-3', 'pure order');
  eq(O.cartOrder(4, 'molden').map(c => c.join('')).slice(3, 6).join(' '), '310 301 130', 'Molden g order starts xxxx yyyy zzzz xxxy xxxz yyyx');
  eq(O.cartOrder(4, 'fchk').slice(0, 2).map(c => c.join('')).join(' '), '004 013', 'Gaussian g order starts zzzz yzzz');
});

test('overlap of two s Gaussians, analytic', () => {
  const a = 0.8, b = 0.35, R = 1.7;
  const ds = { basis: { source: 'fchk', nbf: 2, shells: [{ atom: 0, center: [0, 0, 0], l: 0, pure: false, exps: [a], coefs: [1] }, { atom: 1, center: [0, 0, R], l: 0, pure: false, exps: [b], coefs: [1] }] },
    mos: { restricted: true, alpha: { energies: [0, 0], occupations: [2, 0], symmetries: null, coeffs: Float64Array.from([1, 0, 0, 1]) }, beta: null } };
  const M = O.model(ds), Sm = O.overlap(M, 'normalized');
  near(Sm[1], Math.pow(2 * Math.sqrt(a * b) / (a + b), 1.5) * Math.exp(-a * b / (a + b) * R * R), 1e-14, 'S12');
  near(Sm[0], 1, 1e-14, 'S11');
  // a normalised s function: (2a/pi)^(3/4) at its centre
  near(O.values(M, { kind: 'mo', spin: 'alpha', index: 0 }, [0, 0, 0])[0], Math.pow(2 * a / Math.PI, 0.75) + Math.pow(2 * b / Math.PI, 0) * 0, 1e-6 + Math.pow(2 * b / Math.PI, 0.75) * Math.exp(-b * R * R) * 0, 'value at the centre');
});

test('grid evaluation matches point evaluation (line screening)', () => {
  for (const f of ['h2o_631gd_pure.fchk', 'fg_cart.molden', 'nh2_uhf_631g.fchk']) {
    const { M } = load(f), g = O.gridFor(M, { spacing: 0.9, margin: 4 });
    const fr = O.frontier(M, 'alpha');
    for (const what of [{ kind: 'mo', spin: 'alpha', index: Math.max(0, fr.homo) }, { kind: 'density' }]) {
      const v = O.evaluate(M, what, g), [n1, n2, n3] = g.dims, pts = [], ref = [];
      for (let t = 0; t < 40; t++) {
        const i = (t * 7) % n1, j = (t * 11) % n2, k = (t * 13) % n3;
        pts.push(g.origin[0] + i * g.axes[0][0], g.origin[1] + j * g.axes[1][1], g.origin[2] + k * g.axes[2][2]);
        ref.push(v[(i * n2 + j) * n3 + k]);
      }
      const pv = O.values(M, what, pts);
      ok(maxDiff(pv, ref) < 1e-6 * Math.max(1, ...ref.map(Math.abs)), f + ' ' + what.kind + ': grid and points agree (' + maxDiff(pv, ref).toExponential(1) + ')');
    }
  }
});

test('integrals on a grid: the HOMO is normalised, one unpaired electron in NH2', () => {
  const h = 0.2;
  { const { M } = load('h2o_631gd_pure.fchk'), g = O.gridFor(M, { spacing: h, margin: 6 }), fr = O.frontier(M, 'alpha');
    const v = O.evaluate(M, { kind: 'mo', spin: 'alpha', index: fr.homo }, g);
    let s = 0; for (let i = 0; i < v.length; i++) s += v[i] * v[i];
    near(s * h * h * h, 1, 0.01, '∫ψ² of the HOMO'); }
  { const { M } = load('nh2_uhf_631g.fchk'), g = O.gridFor(M, { spacing: h, margin: 6 });
    const v = O.evaluate(M, { kind: 'spin' }, g);
    let s = 0; for (let i = 0; i < v.length; i++) s += v[i];
    near(s * h * h * h, 1, 0.03, '∫ spin density = one unpaired electron'); }
});

test('menus: HOMO and LUMO, labels, grid box', () => {
  const { M } = load('nh2_uhf_631g.molden');
  const fa = O.frontier(M, 'alpha'), fb = O.frontier(M, 'beta');
  eq(fa.homo, 4, 'alpha HOMO (5 alpha electrons)'); eq(fb.homo, 3, 'beta HOMO (4 beta electrons)');
  const L = O.list(M);
  eq(L.filter(x => x.tag === 'HOMO').length, 2, 'a HOMO per spin');
  eq(L.find(x => x.spin === 'beta' && x.index === 4).tag, 'LUMO', 'beta LUMO');
  eq(O.label(M, { kind: 'mo', spin: 'beta', index: 3 }), 'MO 4 β (HOMO)', 'label');
  const g = O.gridFor(M, { spacing: 0.5, margin: 5 });
  ok(g.dims.every(n => n >= 21), 'box covers the molecule and the margin: ' + g.dims.join('x'));
  const big = O.gridFor(M, { spacing: 0.05, margin: 5, maxPoints: 2e5 });
  ok(big.dims[0] * big.dims[1] * big.dims[2] <= 2e5 && big.spacing > 0.05, 'grid widened to stay under maxPoints');
  const { M: R } = load('h2o_sto3g.fchk');
  eq(O.list(R).find(x => x.index === 4).tag, 'HOMO', 'restricted HOMO');
  eq(O.label(R, { kind: 'density' }), 'Electron density', 'density label');
  let threw = false; try { O.job(R, { kind: 'spin' }); } catch (e) { threw = /unpaired/.test(e.message); }
  ok(threw, 'no spin density for a closed shell');
});

// ------------------------------------------------------------------ file variants
test('Molden: orbitals in symmetry blocks are sorted by energy', () => {
  const t = read('h2o_sto3g.molden');
  const head = t.slice(0, t.indexOf(' Sym='));
  const blocks = t.slice(t.indexOf(' Sym=')).split(/(?= Sym=)/);
  const shuffled = head + [blocks[4], blocks[0], blocks[5], blocks[1], blocks[2], blocks[3], blocks[6]].join('');
  const ds = BM.parseFile('w.molden', shuffled), ref = load('h2o_sto3g.molden').ds;
  ok(ds.mos.alpha.energies.every((e, i) => i === 0 || e >= ds.mos.alpha.energies[i - 1]), 'ascending');
  eq(ds.mos.alpha.fileIndex.join(','), '2,4,5,6,1,3,7', 'each orbital\'s place in the file');
  ok(maxDiff(ds.mos.alpha.coeffs, ref.mos.alpha.coeffs) === 0, 'the same orbitals as the sorted file');
  eq(ds.orbitals.alpha.nOcc, 5, 'occupied orbitals for the level plot');
});
test('Molden: unusable orbital sections give a warning, not an error', () => {
  const t = read('h2o_sto3g.molden');
  const cases = [
    [t.replace('[GTO]', '[STO]'), /Slater/],
    [t.replace(/\[GTO\][\s\S]*?(?=\[MO\])/, ''), /no \[GTO\]/],
    [t.replace('   2 0\n', '   9 0\n'), /atom 9/],
    [t + '    8   0.5\n', /coefficients for 8 basis functions/],
    [t.replace(' sp    3 1.00', ' sp    4 1.00'), /fewer primitives|not a shell/],
  ];
  for (const [txt, re] of cases) {
    const ds = BM.parseFile('w.molden', txt);
    ok(!ds.basis && ds.frames.length === 1 && ds.warnings.some(w => re.test(w)), re + ': ' + ds.warnings.join(' | '));
  }
  // Fortran D exponents, an Angstrom [Atoms] section and an ORCA title
  const orca = read('h2o_631gd_pure_orca.molden').replace(/(\d)E([-+])/g, '$1D$2').replace(' h2o_631gd_pure', ' Molden file created by orca_2mkl for BaseName=x');
  const ds = BM.parseFile('o.molden', orca), M = O.model(ds);
  eq(ds.basis.program, 'orca', 'ORCA recognised from the title');
  eq(O.detect(M).tried[0].conv.sign, 'orca', 'ORCA\'s convention tried first');
  ok(O.detect(M).error < 1e-5, 'and it fits');
});
test('fchk: touching numbers (Psi4) and shells written at the first geometry (Q-Chem)', () => {
  const t = read('h2o_sto3g.fchk'), ref = load('h2o_sto3g.fchk').ds;
  // squeeze the real arrays so that negative numbers touch (as Psi4 writes three-digit exponents)
  const squeezed = t.split('\n').map(l => (/^\s+-?\d\.\d{8}E[-+]\d\d/.test(l) ? l.replace(/ +(-\d\.\d{8}E[-+]\d\d)/g, '$1') : l)).join('\n');
  ok(squeezed !== t && /\d-\d\.\d/.test(squeezed), 'the test file has touching numbers');
  const ds = BM.parseFile('p.fchk', squeezed);
  ok(maxDiff(ds.mos.alpha.coeffs, ref.mos.alpha.coeffs) === 0 && ds.basis.nbf === 7, 'read as the original');
  ok(maxDiff(ds.frames[0].xyz, ref.frames[0].xyz) === 0, 'coordinates too');
  // a later geometry (moved by 0.3 Bohr in x) before the orbitals: the shells follow it
  const L = t.split('\n'), c = L.findIndex(l => l.startsWith('Current cartesian coordinates')), m = L.findIndex(l => l.startsWith('Alpha MO coefficients'));
  const block = [L[c]], moved = [];
  for (let k = c + 1; L[k][0] === ' '; k++) block.push(L[k]);
  const vals = block.slice(1).join(' ').trim().split(/\s+/).map(Number).map((v, i) => (i % 3 === 0 ? v + 0.3 : v));
  moved.push(block[0]);
  for (let k = 0; k < vals.length; k += 5) moved.push(vals.slice(k, k + 5).map(v => (v < 0 ? '' : ' ') + v.toExponential(8).toUpperCase().replace(/E([-+])(\d)$/, 'E$10$2')).map(s => s.padStart(16)).join(''));
  const q = L.slice(0, m).concat(moved, L.slice(m)).join('\n');
  const dq = BM.parseFile('q.fchk', q);
  ok(dq.basis.shells.every((s, i) => Math.abs(s.center[0] - ref.basis.shells[i].center[0] - 0.3) < 1e-7 && s.center[1] === ref.basis.shells[i].center[1]), 'shells moved to the orbitals\' geometry');
  eq(dq.frames.length, 2, 'two geometries');
});

test('Molden export: every fixture written in the Molden standard and read back', () => {
  for (const f of FILES) {
    const { ds, J, M } = load(f), fr = ds.frames[ds.frames.length - 1];
    const text = BM.qcWriters.molden({ Z: fr.Z, xyz: fr.xyz, title: f, orbitalModel: M });
    const back = BM.parseFile('back.molden', text), M2 = O.model(back);
    ok(M2 && M2.nbf === M.nbf, f + ': basis read back');
    ok(/\[Atoms\] AU/.test(text) && M2.shells.every((s2, i) => [0, 1, 2].every(a => Math.abs(s2.center[a] - M.shells[i].center[a]) < 1e-9)), f + ': the atoms at the basis centres (Bohr)');
    if (!M2) continue;
    const det = O.detect(M2);
    ok(det.error < 1e-5 && det.conv.prim === 'normalized' && det.conv.cart === 'component' && det.conv.sign === 'standard', f + ': written in the standard convention (' + JSON.stringify(det.conv) + ', ' + det.error.toExponential(1) + ')');
    const pts = J.expect.points_bohr.flat(), fa = O.frontier(M, 'alpha');
    for (const what of [{ kind: 'mo', spin: 'alpha', index: fa.homo }, { kind: 'mo', spin: 'alpha', index: M.nbf - 1 }, { kind: 'density' }]) {
      const a = O.values(M, what, pts), b2 = O.values(M2, what, pts);
      ok(maxDiff(a, b2) < 1e-8 * Math.max(1, ...Array.from(a, Math.abs)), f + ': ' + what.kind + ' ' + (what.index !== undefined ? what.index + 1 : '') + ' the same after the round trip (' + maxDiff(a, b2).toExponential(1) + ')');
    }
    eq(back.mos.restricted, ds.mos.restricted, f + ': restricted or not');
  }
  // another geometry: no orbitals written
  const { ds, M } = load('h2o_sto3g.fchk'), fr = ds.frames[0], moved = Float64Array.from(fr.xyz, (v, i) => v + (i === 0 ? 0.1 : 0));
  ok(!/\[GTO\]/.test(BM.qcWriters.molden({ Z: fr.Z, xyz: moved, title: 'x', orbitalModel: M })), 'no orbitals for a moved geometry');
});

test('cut files: no crash, and a basis only when it is complete', () => {
  let n = 0;
  for (const f of FILES) {
    const t = read(f);
    for (let k = 1; k <= 24; k++) {
      const cut = t.slice(0, Math.floor(t.length * k / 25));
      let ds = null;
      try { ds = BM.parseFile(f, cut); } catch (e) { ok(e instanceof Error && e.message.length > 0, f + ' cut ' + k + ': clean error'); continue; }
      n++;
      if (ds.basis) {
        const M = O.model(ds);
        if (M) {
          let threw = null;
          try { O.values(M, { kind: 'mo', spin: 'alpha', index: 0 }, [0.1, 0.2, 0.3]); } catch (e) { threw = e; }
          ok(!threw, f + ' cut ' + k + ': evaluates (' + (threw && threw.message) + ')');
        }
      }
    }
  }
  ok(n > 100, 'cut files parsed: ' + n);
});

test('level plot uses the file\'s occupations', () => {
  const ds = { orbitals: { alpha: { energies: [-2, -1, -0.5, 0.2], nOcc: 2, occupations: [2, 0, 2, 0] }, beta: null } };
  const c = BM.analysis.moLevels(ds)[0];
  eq(c.levels.map(l => l.occ ? 1 : 0).join(''), '1010', 'occupied flags'); eq(c.homo, 2, 'HOMO'); eq(c.lumo, 3, 'LUMO');
  const c2 = BM.analysis.moLevels({ orbitals: { alpha: { energies: [-2, -1, 0.3], nOcc: 2 }, beta: null } })[0];
  eq(c2.homo, 1, 'HOMO from nOcc'); eq(c2.lumo, 2, 'LUMO from nOcc');
});
test('cube files: an electrostatic potential is a potential', () => {
  const cube = ' Electrostatic potential from Total SCF Density\n x\n    1    0.0 0.0 0.0\n    2  0.5 0 0\n    2  0 0.5 0\n    2  0 0 0.5\n    1  1.0  0 0 0\n 1 2 3 4 5 6 7 8\n';
  eq(BM.parseFile('p.cube', cube).volumes[0].kind, 'potential', 'kind');
  eq(BM.parseFile('d.cube', cube.replace('Electrostatic potential', 'Electron density')).volumes[0].kind, 'density', 'a density stays a density');
});

// ------------------------------------------------------------------ background computation
(async () => {
  await testAsync('compute(): on the page when there are no workers, cancel', async () => {
    const { M } = load('nh3_631g.fchk'), g = O.gridFor(M, { spacing: 0.6 }), what = { kind: 'density' };
    const fr = [];
    const job = O.compute(M, what, g, { workers: 0, onProgress: f => fr.push(f) });
    const vol = await job.promise;
    const ref = O.evaluate(M, what, g);
    ok(maxDiff(vol.data, ref) === 0, 'the same values as evaluate()');
    ok(fr.length >= 1 && fr[fr.length - 1] === 1, 'progress reaches 1');
    eq(vol.kind, 'density', 'kind'); eq(vol.name, 'Electron density', 'name');
    near(vol.axes[0][0], g.axes[0][0] * 0.529177210903, 1e-12, 'axes in Å');
    near(vol.origin[2], g.origin[2] * 0.529177210903, 1e-12, 'origin in Å');
    const j2 = O.compute(M, what, O.gridFor(M, { spacing: 0.1 }), { workers: 0 });
    j2.cancel();
    let msg = '';
    try { await j2.promise; } catch (e) { msg = e.message; }
    eq(msg, 'cancelled', 'cancelled');
  });
  await testAsync('the worker code runs on its own', async () => {
    const { M } = load('fg_pure.molden'), g = O.gridFor(M, { spacing: 0.8 }), what = { kind: 'mo', spin: 'alpha', index: 11 };
    const posted = [];
    const self = { postMessage: m => posted.push(m) };
    new Function('self', O.workerSource())(self);
    self.onmessage({ data: { job: O.job(M, what) } });
    const ax = [].concat(g.axes[0], g.axes[1], g.axes[2]);
    self.onmessage({ data: { o: g.origin, ax, dims: g.dims, i0: 2, i1: 5 } });
    const ref = O.evaluate(M, what, g), plane = g.dims[1] * g.dims[2];
    ok(posted.length === 2 && posted[0].ready, 'the worker says it runs');
    ok(posted[1] && posted[1].i0 === 2 && maxDiff(posted[1].out, ref.subarray(2 * plane, 5 * plane)) === 0, 'slab from the worker code');
    const msPosted = [];
    const core = new Function('return ' + S.fieldCore.toString())();
    const r = core([0, 0, 0], [3.1], 'ses', 1.4, 0.3, 1e6);
    msPosted.push(r);
    ok(r.dims[0] > 10 && r.data.length === r.dims[0] * r.dims[1] * r.dims[2], 'molecular surface core on its own');
  });
  await testAsync('workers that fail or never answer: the page takes over', async () => {
    const { M } = load('nh3_631g.fchk'), g = O.gridFor(M, { spacing: 0.7 }), what = { kind: 'density' }, ref = O.evaluate(M, what, g);
    const saved = { Worker: globalThis.Worker, Blob: globalThis.Blob, c: URL.createObjectURL };
    let made = 0;
    // a stand-in for Worker: 'silent' never answers, 'broken' reports an error, 'good' runs the worker source
    const kinds = ['silent', 'broken', 'good'];
    let kind = 'silent';
    globalThis.Blob = class { constructor(parts) { this.text = parts.join(''); } };
    URL.createObjectURL = b => { made++; return 'blob:test-' + made; };
    globalThis.Worker = class {
      constructor() {
        this.kind = kind;
        const self = { postMessage: (m) => setTimeout(() => this.onmessage && this.onmessage({ data: m }), 0) };
        if (this.kind === 'good') new Function('self', O.workerSource())(self);
        this.self = self;
      }
      postMessage(m) {
        if (this.kind === 'silent') return;
        if (this.kind === 'broken') { setTimeout(() => this.onerror && this.onerror({ message: 'cannot start', preventDefault() {} }), 0); return; }
        setTimeout(() => this.self.onmessage({ data: m }), 0);
      }
      terminate() { this.dead = true; }
    };
    try {
      for (kind of kinds) {
        const t0 = Date.now(), vol = await O.compute(M, what, g, { workers: 2, watchdogMs: 200 }).promise;
        ok(maxDiff(vol.data, ref) === 0, kind + ' workers: the right values');
        eq(vol.computed.workers, kind === 'good' ? 2 : 0, kind + ' workers: who computed');
        if (kind === 'silent') ok(Date.now() - t0 >= 150, 'the page waited for the watchdog');
      }
      kind = 'silent';
      const ms = await S.fieldAsync([6], [0, 0, 0], { type: 'vdw', watchdogMs: 200 }).promise;
      ok(ms.kind === 'molsurf' && ms.data.length > 100, 'molecular surface: the page takes over from a silent worker');
    } finally {
      globalThis.Worker = saved.Worker; globalThis.Blob = saved.Blob; URL.createObjectURL = saved.c;
    }
  });
  await testAsync('molecular surfaces', async () => {
    const at = (vol, x, y, z) => S.sample(vol, [x, y, z])[0];
    const crossing = (vol, u, r0, r1) => {
      let a = r0, b = r1;
      const fa0 = at(vol, u[0] * a, u[1] * a, u[2] * a);
      for (let it = 0; it < 60; it++) { const m = (a + b) / 2, fm = at(vol, u[0] * m, u[1] * m, u[2] * m); if ((fm > 0) === (fa0 > 0)) a = m; else b = m; }
      return (a + b) / 2;
    };
    const R = BM.el.vdw[6], dir = [0.6, 0.48, 0.64];
    near(crossing(S.field([6], [0, 0, 0], { type: 'vdw', spacing: 0.2 }), dir, 0.5, 4), R, 0.01, 'one atom: van der Waals radius');
    near(crossing(S.field([6], [0, 0, 0], { type: 'sas', spacing: 0.2 }), dir, 0.5, 4), R + 1.4, 0.01, 'one atom: accessible radius');
    near(crossing(S.field([6], [0, 0, 0], { type: 'ses', spacing: 0.2 }), dir, 0.5, 4), R, 0.01, 'one atom: excluded surface = van der Waals');
    near(crossing(S.field([6], [0, 0, 0], { type: 'sas', probe: 2, spacing: 0.2 }), dir, 0.5, 5), R + 2, 0.01, 'probe radius');
    // two atoms 3 Å apart: the probe fills the groove; at the midplane the excluded surface is a torus
    const two = S.field([6, 6], [-1.5, 0, 0, 1.5, 0, 0], { type: 'ses', spacing: 0.15 });
    near(crossing(two, [0, 0.6, 0.8], 0.1, 4), Math.sqrt((R + 1.4) ** 2 - 1.5 ** 2) - 1.4, 0.01, 're-entrant surface at the midplane');
    const mesh = BM.iso.extract(two, 0);
    ok(mesh.indices.length > 1000, 'surface mesh: ' + mesh.indices.length / 3 + ' triangles');
    const vol = await S.fieldAsync([8, 1, 1], [0, 0, 0.12, 0, 0.76, -0.47, 0, -0.76, -0.47], { type: 'vdw' }).promise;
    eq(vol.kind, 'molsurf', 'kind'); eq(vol.surface.type, 'vdw', 'type'); ok(/van der Waals/.test(vol.name), 'name: ' + vol.name);
    let err = ''; try { S.field([], [], {}); } catch (e) { err = e.message; }
    ok(/no atoms/.test(err), 'no atoms: ' + err);
  });
  await testAsync('electrostatic potential on surfaces', async () => {
    near(S.espCharges([1], [0, 0, 0], [1], [0.529177210903, 0, 0])[0], 1, 1e-12, '+1 at 1 Bohr: 1 Hartree/e');
    near(S.espCharges([8, 1, 1], [0, 0, 0, 1, 0, 0, -1, 0, 0], [-0.8, 0.4, 0.4], [0, 3, 0])[0], (-0.8 / 3 + 0.8 / Math.sqrt(10)) * 0.529177210903, 1e-12, 'sum over atoms');
    // trilinear sampling is exact for a linear field
    const vol = { origin: [-1, -2, 0.5], axes: [[0.5, 0, 0], [0, 0.25, 0], [0.1, 0, 0.4]], dims: [6, 9, 5], data: new Float32Array(6 * 9 * 5) };
    for (let i = 0; i < 6; i++) for (let j = 0; j < 9; j++) for (let k = 0; k < 5; k++) {
      const x = -1 + 0.5 * i + 0.1 * k, y = -2 + 0.25 * j, z = 0.5 + 0.4 * k;
      vol.data[(i * 9 + j) * 5 + k] = 2 * x - y + 3 * z;
    }
    const p = [0.3, -1.1, 1.3], v = S.sample(vol, p)[0];
    near(v, 2 * 0.3 + 1.1 + 3 * 1.3, 1e-5, 'trilinear on a skewed grid');
    ok(Number.isNaN(S.sample(vol, [9, 9, 9])[0]), 'NaN outside');
    const c = S.colors([-1, 0, 1, NaN], -1, 1);
    ok(c[0] > 0.8 && c[2] < 0.2, 'negative is red'); ok(c[3] > 0.9 && c[4] > 0.9 && c[5] > 0.9, 'zero is white'); ok(c[8] > 0.8 && c[6] < 0.2, 'positive is blue'); ok(Math.abs(c[9] - 0.6) < 1e-6, 'NaN is grey');
  });
  console.log('BlueMolecule orbital and surface tests: ' + pass + ' passed, ' + fail + ' failed');
  if (fail) { for (const f of failures.slice(0, 40)) console.log('  FAIL ' + f); process.exit(1); }
})();
