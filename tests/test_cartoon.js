#!/usr/bin/env node
/* BlueMolecule tests for src/cartoon.js: simplified DSSP on generated ideal structures (α-helix, antiparallel sheet,
 * isolated strand, 3-10 helix), file secondary structure (PDB HELIX/SHEET, mmCIF struct_conf) taking precedence, and
 * mesh sanity for cartoons (ribbons, arrows, tubes, helix cylinders), nucleic tubes with base slabs, traces, filled
 * rings and SNFG symbols with link sticks; vertex budgets for 5000 residues.
 * Run: node tests/test_cartoon.js */
'use strict';
const path = require('path');
const ROOT = path.join(__dirname, '..');
for (const f of ['elements.js', 'geom.js', 'bigfile.js', 'parsers.js', 'readers-struct.js', 'chem.js', 'reps.js', 'cartoon.js']) require(path.join(ROOT, 'src', f));
const BM = globalThis.BM, R = BM.reps, CT = BM.cartoon, S = require('./synth.js');

let pass = 0, fail = 0, current = '';
const failures = [];
const ok = (c, msg) => { if (c) pass++; else { fail++; failures.push('[' + current + '] ' + msg); } };
const eq = (a, b, msg) => ok(a === b, msg + ': got ' + JSON.stringify(a) + ', expected ' + JSON.stringify(b));
function test(name, fn) { current = name; try { fn(); } catch (e) { fail++; failures.push('[' + name + '] threw: ' + (e && e.stack || e)); } }

function sysOf(chains) {
  const sys = { Z: [], xyz: [], names: [], res: [], ids: [], chains: [] };
  chains.forEach(([bb, ch, name]) => bb.forEach((r, k) => { for (const [nm, z, p] of S.residueAtoms(r, name || 'ALA')) { sys.Z.push(z); sys.xyz.push(...p); sys.names.push(nm); sys.res.push(name || 'ALA'); sys.ids.push(k + 1); sys.chains.push(ch); } }));
  return sys;
}
function ctxOf(sys, extra) {
  const f = { Z: sys.Z, xyz: Float64Array.from(sys.xyz), atomNames: sys.names, resNames: sys.res, resIds: sys.ids, chains: sys.chains };
  return R.context(Object.assign({ frame: f, bonds: BM.geom.computeBonds(f.Z, f.xyz, { tolerance: 0.4 }) }, extra || {}));
}
const ssString = (ctx, ss) => { const cls = R.classify(ctx); let s = ''; for (let r = 0; r < cls.nRes; r++) if (cls.resClass[r] === 1) s += '-HE'[ss[r]]; return s; };
/** Sanity of every mesh: finite values, indices in range, unit normals; returns {verts, tris, meshes}. */
function meshCheck(g, label) {
  let verts = 0, tris = 0, bad = 0;
  for (const m of g.meshes) {
    const nv = m.positions.length / 3;
    verts += nv; tris += m.indices.length / 3;
    if (m.normals.length !== m.positions.length || m.colors.length !== m.positions.length || m.indices.length % 3) bad++;
    for (let i = 0; i < m.positions.length; i++) if (!Number.isFinite(m.positions[i]) || !Number.isFinite(m.normals[i]) || !(m.colors[i] >= 0 && m.colors[i] <= 1.0001)) { bad++; break; }
    for (let i = 0; i < m.indices.length; i++) if (m.indices[i] >= nv) { bad++; break; }
    for (let v = 0; v < nv; v++) { const l = Math.hypot(m.normals[3 * v], m.normals[3 * v + 1], m.normals[3 * v + 2]); if (Math.abs(l - 1) > 2e-3) { bad++; break; } }
  }
  ok(bad === 0 && tris > 0, label + ': mesh sane (' + verts + ' vertices, ' + tris + ' triangles)');
  return { verts, tris, meshes: g.meshes };
}
/** Largest distance from a mesh vertex to the nearest of the points. */
function maxDistTo(meshes, pts) {
  let worst = 0;
  for (const m of meshes) for (let v = 0; v < m.positions.length / 3; v++) {
    let best = Infinity;
    for (const p of pts) { const d = Math.hypot(m.positions[3 * v] - p[0], m.positions[3 * v + 1] - p[1], m.positions[3 * v + 2] - p[2]); if (d < best) best = d; }
    if (best > worst) worst = best;
  }
  return worst;
}
/** Largest distance from a mesh vertex to the CA trace (segments between consecutive CAs < 4.3 Å apart). */
function maxDistToTrace(meshes, pts) {
  const segs = [];
  for (let k = 0; k + 1 < pts.length; k++) if (Math.hypot(pts[k + 1][0] - pts[k][0], pts[k + 1][1] - pts[k][1], pts[k + 1][2] - pts[k][2]) < 4.3) segs.push([pts[k], pts[k + 1]]);
  let worst = 0;
  for (const m of meshes) for (let v = 0; v < m.positions.length / 3; v++) {
    const p = [m.positions[3 * v], m.positions[3 * v + 1], m.positions[3 * v + 2]];
    let best = Infinity;
    for (const [a, b] of segs) {
      const ab = [b[0] - a[0], b[1] - a[1], b[2] - a[2]], ap = [p[0] - a[0], p[1] - a[1], p[2] - a[2]];
      const t = Math.max(0, Math.min(1, (ab[0] * ap[0] + ab[1] * ap[1] + ab[2] * ap[2]) / (ab[0] * ab[0] + ab[1] * ab[1] + ab[2] * ab[2])));
      const d = Math.hypot(ap[0] - t * ab[0], ap[1] - t * ab[1], ap[2] - t * ab[2]);
      if (d < best) best = d;
    }
    if (best > worst) worst = best;
  }
  return worst;
}
const caPts = ctx => { const out = []; const f = ctx.frame; f.atomNames.forEach((n, i) => { if (n === 'CA') out.push([f.xyz[3 * i], f.xyz[3 * i + 1], f.xyz[3 * i + 2]]); }); return out; };

// ---------------------------------------------------------------- secondary structure
test('DSSP: ideal alpha helix', () => {
  const c = ctxOf(sysOf([[S.idealHelix(20), 'A']])), s = ssString(c, CT.dssp(c));
  eq(s.length, 20, 'residues');
  ok(/^-?H{17,20}-?$/.test(s), 'helix over the interior residues: ' + s);
  eq(R.count(R.evaluate('helix', c)), R.count(R.evaluate('protein and resid 2-19', c)) + (s[0] === 'H' ? 5 : 0) + (s[19] === 'H' ? 5 : 0), '"helix" keyword follows the assignment');
});
test('DSSP: antiparallel two-stranded sheet', () => {
  const [a, b] = S.antiparallelPair(8), c = ctxOf(sysOf([[a, 'A'], [b, 'B']])), s = ssString(c, CT.dssp(c));
  const sa = s.slice(0, 8), sb = s.slice(8);
  ok((sa.match(/E/g) || []).length >= 5 && (sb.match(/E/g) || []).length >= 5 && !/H/.test(s), 'both strands E: ' + sa + ' / ' + sb);
  // the energy function itself: an ideal linear N-H...O=C bond is strongly attractive
  const e = CT.hbEnergy([0, 0, 0], [1.23, 0, 0], [4.13, 0, 0], [3.13, 0, 0]);
  ok(e < -2.5 && e > -4.5, 'linear H-bond energy ' + e.toFixed(2) + ' kcal/mol');
});
test('DSSP: isolated strand, 3-10 helix and a mixed protein', () => {
  let c = ctxOf(sysOf([[S.idealStrand(10), 'A']]));
  eq(ssString(c, CT.dssp(c)), '----------', 'isolated strand is coil (no partner)');
  c = ctxOf(sysOf([[S.backbone(Array(12).fill(-49), Array(12).fill(-26)), 'A']]));
  ok(!/E/.test(ssString(c, CT.dssp(c))), '3-10 helix is never a strand: ' + ssString(c, CT.dssp(c)));
  const sys = S.protein(BM, { ligand: false, water: false });
  c = ctxOf(sys);
  const file = R.secondary(R.context(Object.assign({}, c, { cache: new Map(), _cls: undefined, ssRecords: sys.ss })));
  const d = CT.dssp(c), cls = R.classify(c);
  let hel = 0, helOk = 0, sh = 0, shOk = 0;
  for (let r = 0; r < cls.nRes; r++) { if (file[r] === 1) { hel++; if (d[r] === 1) helOk++; } if (file[r] === 2) { sh++; if (d[r] === 2) shOk++; } }
  ok(helOk / hel > 0.8, 'DSSP finds the built helices: ' + helOk + '/' + hel);
  ok(shOk / sh >= 0.5, 'DSSP finds the built strands (edges next to rough loops may be coil): ' + shOk + '/' + sh);
});
test('file secondary structure takes precedence (PDB HELIX/SHEET, mmCIF)', () => {
  const sys = S.protein(BM, { ligand: false, water: false });
  sys.ss.push({ type: 'helix', chain: 'A', start: { resId: 17 }, end: { resId: 18 } }); // a loop the file calls helix
  const ds = BM.parseFile('t.pdb', S.toPDB(sys, { ss: true })), f = ds.frames[0];
  ok(ds.secondaryStructure && ds.secondaryStructure.length === sys.ss.length, 'HELIX/SHEET records read: ' + (ds.secondaryStructure || []).length);
  ok(f.chains && f.chains.includes('A'), 'chain ids read');
  const c = R.context({ frame: f, bonds: BM.geom.computeBonds(f.Z, f.xyz, { tolerance: 0.4 }), ssRecords: ds.secondaryStructure });
  const ss = R.secondary(c), cls = R.classify(c);
  const at = id => ss[cls.resId.findIndex((x, r) => x === id && cls.resClass[r] === 1)];
  ok(at(17) === 1 && at(18) === 1, 'loop residues marked helix by the file stay helix');
  ok(at(64) === 2 && at(5) === 1, 'sheet and helix from records');
  eq(cls.ssSource, 'file', 'source recorded');
  // other chains are not affected by a record for chain A
  sys.chains = sys.chains.map(x => (x === 'A' ? 'B' : x));
  const ds2 = BM.parseFile('t2.pdb', S.toPDB(sys, { ss: true }).replace(/(HELIX.{12})B/g, '$1A'));
  const f2 = ds2.frames[0], c2 = R.context({ frame: f2, bonds: BM.geom.computeBonds(f2.Z, f2.xyz, { tolerance: 0.4 }), ssRecords: [{ type: 'helix', chain: 'Z', start: { resId: 1 }, end: { resId: 200 } }] });
  const s2 = R.secondary(c2);
  ok(!Array.from(s2).some(v => v === 1), 'a record for another chain marks nothing');
});

// ---------------------------------------------------------------- meshes
test('protein cartoon: ribbons, arrows, tubes', () => {
  const sys = S.protein(BM), c = ctxOf(sys, { ssRecords: sys.ss }), mask = R.evaluate('protein', c), cols = R.colorize('ss', c, mask);
  const g = CT.build('cartoon', c, mask, cols, {}), m = meshCheck(g, 'cartoon');
  eq(g.proxies.length, R.count(R.evaluate('protein and name CA', c)), 'one pick proxy per CA');
  const all = caPts(c);
  const md = maxDistToTrace(m.meshes, all);
  ok(md > 1.3 && md < 2.4, 'widest point = sheet arrowheads: ' + md.toFixed(2));
  // cross-section profile: strand body 1.1 Å half-width, arrowhead 1.75 Å narrowing to the coil, helix 1.25 Å, coil 0.3 Å
  const prof = [];
  CT.build('cartoon', c, mask, cols, { profile: prof });
  const maxA = f => Math.max(...prof.filter(f).map(p => p.a));
  const nearly = (a, b) => Math.abs(a - b) < 0.02;
  ok(nearly(maxA(p => p.ss === 2 && !p.arrow), 1.1), 'strand half-width ' + maxA(p => p.ss === 2 && !p.arrow));
  ok(nearly(maxA(p => p.arrow), 1.75) && prof.filter(p => p.arrow && p.u === 0).length >= 2, 'arrowheads start at 1.75 Å (one per strand)');
  ok(nearly(maxA(p => p.ss === 1), 1.25) && nearly(Math.min(...prof.map(p => p.a)), 0.3), 'helix ribbon 1.25 Å, coil tube 0.3 Å');
  const prof2 = [];
  CT.build('cartoon', c, mask, cols, { profile: prof2, scale: 1.5 });
  ok(nearly(Math.max(...prof2.map(p => p.a)), 1.75 * 1.5), 'size 1.5x scales the cross-sections');
  const prof3 = [];
  CT.build('cartoon', c, mask, cols, { profile: prof3, cylinders: true });
  ok(prof3.filter(p => p.ss === 1).every(p => p.a < 0.31), 'cylinder mode: no helix ribbons');
  // helix ribbons stay within ~1.35 Å of the CA spline, coil tubes within 0.35 Å
  const cls = R.classify(c), ss = R.secondary(c), bb = CT.proteinChains(c), helixCA = [], coilCA = [];
  for (let k = 0; k < bb.count; k++) { const p = [c.xyz[3 * bb.ca[k]], c.xyz[3 * bb.ca[k] + 1], c.xyz[3 * bb.ca[k] + 2]]; (ss[bb.res[k]] === 1 ? helixCA : coilCA).push(p); }
  void cls;
  const hm = CT.build('cartoon', c, R.evaluate('helix', c), cols, {});
  const hd = maxDistToTrace(hm.meshes, helixCA);
  ok(hd > 1.0 && hd < 1.9, 'helix ribbons (1.25 Å half-width) hug the CA trace: ' + hd.toFixed(2));
  const colours = new Set();
  for (const mm of m.meshes) for (let v = 0; v < mm.colors.length; v += 3) colours.add(R.toHex([mm.colors[v], mm.colors[v + 1], mm.colors[v + 2]]));
  ok(colours.has(R.toHex(R.SS_COL.helix)) && colours.has(R.toHex(R.SS_COL.sheet)) && colours.has(R.toHex(R.SS_COL.coil)), 'secondary-structure colours on the mesh');
  const cy = CT.build('cartoon', c, mask, cols, { cylinders: true }), mc = meshCheck(cy, 'cartoon with helix cylinders');
  ok(maxDistTo(cy.meshes.slice(0, 1), all) > 2.0, 'cylinders reach out to the helix radius');
  void coilCA;
  ok(mc.verts !== m.verts, 'cylinder mode builds different geometry');
  void md;
});
test('trace, nucleic tube with base slabs, rings, SNFG', () => {
  const sys = S.protein(BM), c = ctxOf(sys);
  const tr = CT.build('trace', c, R.evaluate('protein', c), R.colorize('chain', c), {});
  meshCheck(tr, 'trace');
  const td = maxDistToTrace(tr.meshes, caPts(c));
  ok(td > 0.3 && td < 0.9, 'trace tube (0.38 Å) on the CA trace: ' + td.toFixed(2));
  const nuc = S.nucleic(12), cn = ctxOf(nuc), mk = R.evaluate('nucleic', cn);
  eq(R.count(mk), nuc.Z.length, 'synthetic strand classified as nucleic');
  eq(R.autoPreset(cn), 'nucleic', 'nucleic preset');
  const gn = CT.build('cartoon', cn, mk, R.colorize('residue', cn, mk), {});
  meshCheck(gn, 'nucleic cartoon');
  ok(gn.meshes.length === 2 && gn.meshes[1].twoSided, 'tube + two-sided base slabs');
  const slabTris = gn.meshes[1].indices.length / 3;
  eq(slabTris, 6 * ((6 + 5) * 2 + 9 * 2) + 6 * (6 * 2 + 6 * 2), 'slab triangles: purines (two fans of 6+5, 9 side quads) and pyrimidines (6, 6 quads)');
  eq(gn.proxies.length, 12, 'one proxy per nucleotide');
  const colours = new Set();
  const sl = gn.meshes[1];
  for (let v = 0; v < sl.colors.length; v += 3) colours.add(R.toHex([sl.colors[v], sl.colors[v + 1], sl.colors[v + 2]]));
  eq(colours.size, 4, 'slabs coloured by base (A, T, G, C)');
  // rings: the ligand's two benzene rings
  const rg = CT.build('rings', c, R.evaluate('ligand', c), R.colorize('element', c), {});
  meshCheck(rg, 'ligand rings');
  eq(rg.meshes[0].indices.length / 3, 12, 'two six-membered rings -> 12 fan triangles');
  // glycan rings and symbols
  const gs = S.glycan(BM), cg = ctxOf(gs), carb = R.evaluate('carbohydrate', cg);
  const gr = CT.build('rings', cg, carb, R.colorize('sugar', cg, carb), { scheme: 'sugar' });
  meshCheck(gr, 'sugar rings');
  eq(gr.meshes[0].indices.length / 3, 12 * 6, 'twelve pyranose rings (6 triangles each)');
  const sn = CT.build('snfg', cg, carb, R.colorize('element', cg, carb), {});
  const ms = meshCheck(sn, 'SNFG symbols');
  eq(sn.proxies.length, 12, 'one symbol per sugar');
  let links = 0, white = 0;
  for (const mm of ms.meshes) for (let v = 0; v < mm.colors.length; v += 3) {
    const h = R.toHex([mm.colors[v], mm.colors[v + 1], mm.colors[v + 2]]);
    if (h === '#ccd1db') links++; if (h === '#f7f7f7') white++;
  }
  eq(links / 22, 12, 'link sticks: 11 glycosidic + 1 to the serine');
  ok(white > 0, 'divided symbols have a white half');
  const seen = new Set(R.sugars(cg).map(s => s.shape));
  ok(['star', 'sphere', 'ddiamond', 'cube', 'dcube'].every(x => seen.has(x)), 'shapes: ' + [...seen].join(', '));
});
test('vertex budget and timing for 5000 residues', () => {
  const chains = [];
  for (let k = 0; k < 20; k++) {
    const h = S.helixAt(250, [(k % 5) * 26, Math.floor(k / 5) * 26, 0], [0, 0.2, 1], k);
    chains.push([h, String.fromCharCode(65 + k), 'ALA']);
  }
  const sys = sysOf(chains), c = ctxOf(sys);
  let t = Date.now();
  const ss = R.secondary(c);
  const tss = Date.now() - t;
  ok(Array.from(ss).filter(v => v === 1).length > 4500, 'DSSP on 5000 residues finds the helices');
  t = Date.now();
  const mask = R.evaluate('protein', c), g = CT.build('cartoon', c, mask, R.colorize('ss', c, mask), {});
  const tb = Date.now() - t, m = meshCheck(g, '5000-residue cartoon');
  console.log(`  5000 residues (${c.n} atoms): DSSP ${tss} ms, cartoon ${tb} ms, ${m.verts} vertices, ${m.tris} triangles`);
  ok(m.verts < 400000 && m.tris < 800000, 'reasonable vertex count');
  ok(tss < 3000 && tb < 3000, 'fast enough');
});

console.log(`\nBlueMolecule cartoon tests: ${pass} passed, ${fail} failed`);
for (const f of failures) console.log('  FAIL ' + f);
process.exit(fail ? 1 : 0);
