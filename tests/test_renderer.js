// Renderer checks in headless Chromium: bond orders, picking, tool input hooks, overlay hook and extra geometry.
// Drives tests/renderer_harness.html (loads src/elements.js, geom.js, renderer.js only).
// Usage: node tests/test_renderer.js [screenshotDir]
const { chromium } = require('playwright');
const path = require('path');
const fs = require('fs');

const ROOT = path.resolve(__dirname, '..');
const URL = 'file://' + path.join(ROOT, 'tests', 'renderer_harness.html');
const SHOTS = process.argv[2] || path.join(ROOT, 'tests', 'shots', 'renderer');
fs.mkdirSync(SHOTS, { recursive: true });

let failures = 0, passes = 0;
const check = (cond, msg) => { if (cond) passes++; else { failures++; console.log('  FAIL:', msg); } };
const near = (a, b, tol) => Math.abs(a - b) <= tol;

(async () => {
  const browser = await chromium.launch({ args: ['--enable-unsafe-swiftshader'] });
  const ctx = await browser.newContext({ viewport: { width: 820, height: 660 }, deviceScaleFactor: 1 });
  const page = await ctx.newPage();
  const errors = [];
  page.on('pageerror', e => errors.push('pageerror: ' + e.message));
  page.on('console', m => { if (m.type() === 'error') errors.push('console: ' + m.text()); });
  await page.goto(URL);
  await page.waitForFunction(() => window.harnessReady);
  const ev = (fn, arg) => page.evaluate(fn, arg);
  const shot = async name => { await ev(() => frame()); await page.locator('#vp').screenshot({ path: path.join(SHOTS, name + '.png') }); };

  // ---------------- bond orders ----------------
  console.log('Bond orders: benzene (Kekulé)');
  let r = await ev(() => {
    // built under a tilted view: an offset taken from the view would leave the ring plane
    load(benzene(false), { rot: RING_ROT, view: BM.quat.mul(BM.quat.axis([1, 0, 0], -1.0), conj(RING_ROT)) });
    const M = BM.quat.mat(RING_ROT), N = [M[2], M[5], M[8]];
    const dot = (a, b) => a[0] * b[0] + a[1] * b[1] + a[2] * b[2], sub = (a, b) => [a[0] - b[0], a[1] - b[1], a[2] - b[2]];
    const cen = [0, 0, 0]; for (let i = 0; i < 6; i++) for (let k = 0; k < 3; k++) cen[k] += R.xyz[3 * i + k] / 6;
    const lines = readCyl(R.multiInst, R.nMultiCyl);
    let maxOut = 0;
    for (const L of lines) maxOut = Math.max(maxOut, Math.abs(dot(sub(L.a, cen), N)), Math.abs(dot(sub(L.b, cen), N)));
    const axisDist = p => Math.min(...[0, 2, 4].map(k => {
      const a = atomPos(k), u = sub(atomPos((k + 1) % 6), a), t = dot(sub(p, a), u) / dot(u, u);
      return Math.hypot(...sub(p, [a[0] + u[0] * t, a[1] + u[1] * t, a[2] + u[2] * t]));
    }));
    const offs = lines.map(L => axisDist(L.a));
    return { nCyl: R.nCyl, nMain: R.nMainCyl, nMulti: R.nMultiCyl, maxOut, offMin: Math.min(...offs), offMax: Math.max(...offs), rads: [...new Set(lines.map(L => +L.r.toFixed(4)))] };
  });
  check(r.nMain === 18 && r.nMulti === 12 && r.nCyl === 30, 'Kekulé counts (main 18, multi 12, nCyl 30): ' + JSON.stringify(r));
  check(r.maxOut < 1e-4, 'double-bond lines lie in the ring plane (max out-of-plane ' + r.maxOut + ' Å)');
  check(near(r.offMin, 0.1045, 1e-4) && near(r.offMax, 0.1045, 1e-4), 'double-bond lines offset ±0.95·r from the axis: ' + r.offMin + '..' + r.offMax);
  check(r.rads.length === 1 && near(r.rads[0], 0.066, 1e-4), 'double-bond lines thinner (0.6·r): ' + r.rads);
  await ev(() => R.setOrientation(conj(RING_ROT)));
  await shot('01_benzene_kekule_faceon');
  await ev(() => R.setOrientation(BM.quat.mul(BM.quat.axis([1, 0, 0], -1.0), conj(RING_ROT))));
  await shot('02_benzene_kekule_tilted');
  await ev(() => R.setOrientation(BM.quat.mul(BM.quat.axis([1, 0, 0], -88 * Math.PI / 180), conj(RING_ROT))));
  await shot('03_benzene_kekule_edgeon');
  await ev(() => { R.setOrientation(conj(RING_ROT)); R.setOptions({ style: 'wireframe' }); });
  await shot('04_benzene_kekule_wireframe');
  await ev(() => R.setOptions({ style: 'sticks' }));
  await shot('05_benzene_kekule_sticks');

  console.log('Bond orders: aromatic benzene');
  r = await ev(() => {
    load(benzene(true), { rot: RING_ROT, view: BM.quat.mul(BM.quat.axis([1, 0, 0], -1.0), conj(RING_ROT)) });
    const M = BM.quat.mat(RING_ROT), N = [M[2], M[5], M[8]];
    const cen = [0, 0, 0]; for (let i = 0; i < 6; i++) for (let k = 0; k < 3; k++) cen[k] += R.xyz[3 * i + k] / 6;
    const dashes = readCyl(R.multiInst, R.nMultiCyl);
    let maxOut = 0, maxR = 0;
    for (const D of dashes) {
      const m = [(D.a[0] + D.b[0]) / 2 - cen[0], (D.a[1] + D.b[1]) / 2 - cen[1], (D.a[2] + D.b[2]) / 2 - cen[2]];
      maxOut = Math.max(maxOut, Math.abs(m[0] * N[0] + m[1] * N[1] + m[2] * N[2])); maxR = Math.max(maxR, Math.hypot(...m));
    }
    return { nMain: R.nMainCyl, nMulti: R.nMultiCyl, maxOut, maxR };
  });
  check(r.nMain === 24 && r.nMulti >= 12, 'aromatic: normal cylinder for every bond plus dashes ' + JSON.stringify(r));
  check(r.maxOut < 1e-4, 'aromatic dashes lie in the ring plane');
  check(r.maxR < 1.15, 'aromatic dashes sit inside the ring (max centroid distance ' + r.maxR.toFixed(3) + ' Å < apothem 1.20)');
  await ev(() => R.setOrientation(conj(RING_ROT)));
  await shot('06_benzene_aromatic');
  await ev(() => R.setOptions({ style: 'wireframe' }));
  await shot('07_benzene_aromatic_wireframe');

  console.log('Bond orders: persistence');
  r = await ev(() => {
    const { xyz, bonds } = load(benzene(false), { rot: RING_ROT, view: conj(RING_ROT) });
    const orders = () => Float32Array.from(benzene(false).orders), out = { initial: R.nMultiCyl };
    const moved = Float64Array.from(xyz, (v, i) => v + 0.01 * Math.sin(i));
    R.updatePositions(moved); out.afterUpdate = R.nMultiCyl;
    R.setScene(R.Z, moved, bonds, true); out.sameBonds = R.nMultiCyl;
    R.setBonds(bonds); out.setBondsSame = R.nMultiCyl;
    R.setScene(R.Z, moved, Uint32Array.from(bonds), true); out.newBonds = [R.nMultiCyl, R.bondOrders, R.nCyl];
    R.setScene(R.Z, moved, bonds, true, orders()); out.withOrders = R.nMultiCyl;
    R.setBonds(Uint32Array.from(bonds)); out.setBondsNew = R.nMultiCyl;
    R.setBonds(bonds, orders()); out.setBondsOrders = R.nMultiCyl;
    R.setOptions({ style: 'sticks' }); out.sticks = [R.nMultiCyl, R.nCyl];
    R.setOptions({ style: 'spacefill' }); out.spacefill = [R.nMultiCyl, R.nCyl];
    R.setOptions({ style: 'wireframe' }); out.wire = [R.nMultiCyl, R.nCyl];
    R.setOptions({ style: 'ballstick' }); out.back = R.nMultiCyl;
    R.setBondOrders(null); out.cleared = [R.nMultiCyl, R.nCyl];
    return out;
  });
  check(r.initial === 12 && r.afterUpdate === 12, 'orders persist across updatePositions');
  check(r.sameBonds === 12 && r.setBondsSame === 12, 'orders kept when setScene/setBonds get the same bonds');
  check(r.newBonds[0] === 0 && r.newBonds[1] === null && r.newBonds[2] === 24, 'orders reset when setScene gets new bonds without orders');
  check(r.withOrders === 12 && r.setBondsNew === 0 && r.setBondsOrders === 12, 'setScene/setBonds orders argument');
  check(r.sticks.join() === '0,24' && r.spacefill.join() === '0,0' && r.wire.join() === '12,30' && r.back === 12, 'sticks/spacefill stay single: ' + JSON.stringify(r));
  check(r.cleared.join() === '0,24', 'setBondOrders(null) back to singles');

  console.log('Bond orders: CO2, acetylene');
  r = await ev(() => {
    const view = BM.quat.axis([0.3, 1, 0.2], 0.7);
    load(co2(), { view });
    const perpTo = (lines, q) => { const V = BM.quat.mat(q); return Math.max(...lines.map(L => Math.abs(L.a[1] * V[7] + L.a[2] * V[8]) / Math.hypot(L.a[1], L.a[2]))); };
    const lines = readCyl(R.multiInst, R.nMultiCyl), red = BM.el.color[8], grey = BM.el.color[6];
    const res = { n: R.nMultiCyl, main: R.nMainCyl, viewDot: perpTo(lines, view), offs: lines.map(L => Math.hypot(L.a[1], L.a[2])) };
    const atO = x => Math.abs(Math.abs(x) - 1.16) < 1e-4; // half-line ending at an oxygen
    res.split = lines.every(L => { const c = atO(L.a[0]) || atO(L.b[0]) ? red : grey; return [0, 1, 2].every(k => Math.abs(L.c[k] - c[k]) < 1e-5); });
    const view2 = BM.quat.axis([1, 0, 0], 1.1);
    R.setOrientation(view2); R.setBondOrders(Float32Array.from([2, 2]));
    res.viewDot2 = perpTo(readCyl(R.multiInst, R.nMultiCyl), view2);
    R.setOrientation(view); R.setBondOrders(Float32Array.from([2, 2]));
    return res;
  });
  check(r.n === 8 && r.main === 0, 'CO2: two double bonds = 8 half-lines, no single cylinders');
  check(r.viewDot < 1e-4 && r.viewDot2 < 1e-4, 'CO2 (no usable neighbour): offset ⟂ view direction at build time');
  check(r.offs.every(o => near(o, 0.1045, 1e-4)), 'CO2 offsets');
  check(r.split, 'multi-bond colour split at the midpoint (O red / C grey)');
  await shot('08_co2');
  await ev(() => R.setOptions({ style: 'wireframe' }));
  await shot('09_co2_wireframe');
  r = await ev(() => {
    const view = BM.quat.axis([0.3, 1, 0.2], 0.7), V = BM.quat.mat(view);
    load(acetylene(), { view });
    const lines = readCyl(R.multiInst, R.nMultiCyl), offs = lines.map(L => Math.hypot(L.a[1], L.a[2]));
    const perp = lines.filter((L, k) => offs[k] > 1e-6).every(L => Math.abs(L.a[1] * V[7] + L.a[2] * V[8]) / Math.hypot(L.a[1], L.a[2]) < 1e-4);
    return { n: R.nMultiCyl, main: R.nMainCyl, offs: offs.map(o => +o.toFixed(4)).sort(), perp };
  });
  check(r.n === 6 && r.main === 4, 'acetylene: triple = 3 split lines, two C–H singles ' + JSON.stringify(r));
  check(r.perp && r.offs.join() === '0,0,0.154,0.154,0.154,0.154', 'acetylene line offsets 0 and ±1.4·r, ⟂ view: ' + r.offs);
  await shot('10_acetylene');
  await ev(() => R.setOptions({ style: 'wireframe' }));
  await shot('11_acetylene_wireframe');

  // ---------------- extra geometry ----------------
  console.log('Dashed H-bond');
  r = await ev(() => {
    const m = waterDimer();
    load(m, { view: BM.quat.axis([1, 0.2, 0], -0.5) });
    R.setDashedBonds([m.hbond], [0.4, 0.75, 1.0]);
    const onSeg = (i, j) => readCyl(R.decoInst, R.nDeco).every(D => [D.a, D.b].every(p => {
      const a = atomPos(i), b = atomPos(j), u = [b[0] - a[0], b[1] - a[1], b[2] - a[2]], L2 = u[0] ** 2 + u[1] ** 2 + u[2] ** 2;
      const t = ((p[0] - a[0]) * u[0] + (p[1] - a[1]) * u[1] + (p[2] - a[2]) * u[2]) / L2;
      return t > 0 && t < 1 && Math.hypot(a[0] + u[0] * t - p[0], a[1] + u[1] * t - p[1], a[2] + u[2] * t - p[2]) < 1e-4;
    }));
    const res = { n: R.nDeco, onSeg: onSeg(1, 3), main: R.nCyl };
    const x = Float64Array.from(R.xyz); for (const i of [3, 4, 5]) { x[3 * i] += 0.25; x[3 * i + 1] += 0.35; }
    R.updatePositions(x);
    res.n2 = R.nDeco; res.onSeg2 = onSeg(1, 3);
    R.setOptions({ showH: false }); res.hidden = R.nDeco; R.setOptions({ showH: true });
    R.setDashedBonds([[1, 3], [0, 99], [2, 2]]); res.filtered = R.nDeco === res.n2;
    R.setDashedBonds([1, 3]); res.flat = R.nDeco === res.n2;
    load(m, { view: BM.quat.axis([1, 0.2, 0], -0.5) }); R.setDashedBonds([m.hbond], [0.4, 0.75, 1.0]);
    return res;
  });
  check(r.n >= 3 && r.onSeg && r.main === 8, 'H-bond dashes on the H···O segment: ' + JSON.stringify(r));
  check(r.n2 >= 3 && r.onSeg2, 'dashes follow the atoms after updatePositions');
  check(r.hidden === 0 && r.filtered && r.flat, 'hidden / out-of-range pairs skipped; flat pair list accepted');
  await shot('12_water_dimer_hbond');

  console.log('Cell, dipole arrow, markers, bond highlight');
  r = await ev(() => {
    load(formaldehyde(), { view: BM.quat.axis([1, 0.35, 0], -0.75) });
    R.setCell([[5.2, 0, 0], [1.1, 4.8, 0], [0.7, 0.9, 4.6]], [-3.0, -2.9, -2.3]);
    R.setExtraArrows([{ from: [1.4, 0, 1.2], to: [-1.6, 0, 1.2], color: [1.0, 0.55, 0.15], radius: 0.07 }]);
    R.setAtomMarkers([2, 3, 3, 99]);
    R.setSelection([1]);
    R.highlightBond(0);
    R.dist *= 1.9;
    const lines = readCyl(R.multiInst, R.nMultiCyl), cone = readCyl(R.xConeInst, R.nXCone)[0];
    return { deco: R.nDeco, cones: R.nXCone, tip: cone.b, markers: R.nMarker, halo: R.nBondHalo, multi: R.nMultiCyl, maxZ: Math.max(...lines.map(L => Math.max(Math.abs(L.a[2]), Math.abs(L.b[2])))) };
  });
  check(r.deco === 13 && r.cones === 1 && near(r.tip[0], -1.6, 1e-5) && near(r.tip[2], 1.2, 1e-5), 'cell = 12 edges, arrow = shaft + head ending at `to`: ' + JSON.stringify(r));
  check(r.markers === 2 && r.halo === 1, 'markers deduplicated/filtered, bond halo drawn');
  check(r.multi === 4 && r.maxZ < 1e-5, 'C=O lines lie in the molecular plane (neighbour H defines it)');
  await shot('13_cell_dipole_markers');
  await ev(() => R.setOptions({ bg: [1, 1, 1], fog: false }));
  await shot('14_cell_dipole_markers_white');
  r = await ev(() => {
    R.setOptions({ bg: [0.063, 0.090, 0.133], fog: true });
    const out = {};
    R.setScene(R.Z, R.xyz, Uint32Array.from(R.bonds), true); out.haloAfterNewBonds = [R.nBondHalo, R.highlightedBond];
    R.setExtraArrows([]); R.setCell([[5, 0, 0], [0, 5, 0], [0, 0, 0]]); out.cell2d = R.nDeco;
    R.setCell([5, 0, 0, 0, 5, 0, 0, 0, 5], [-2.5, -2.5, -2.5]); out.cellFlat = R.nDeco;
    R.setCell(null); out.none = R.nDeco;
    R.setAtomMarkers(null); out.markers = R.nMarker;
    return out;
  });
  check(r.haloAfterNewBonds.join() === '0,-1', 'bond highlight resets when the bonds change');
  check(r.cell2d === 4 && r.cellFlat === 12 && r.none === 0 && r.markers === 0, 'cell variants and clearing: ' + JSON.stringify(r));

  // ---------------- picking ----------------
  console.log('Picking');
  r = await ev(() => {
    load(benzene(false), { rot: RING_ROT, view: conj(RING_ROT) });
    const scr = p => R.worldToScreen(p);
    const mid = k => { const a = atomPos(R.bonds[2 * k]), b = atomPos(R.bonds[2 * k + 1]); return [(a[0] + b[0]) / 2, (a[1] + b[1]) / 2, (a[2] + b[2]) / 2]; };
    const res = {};
    res.atoms = [...Array(12).keys()].every(i => { const s = scr(atomPos(i)), h = R.pickAny(s.x, s.y); return h.kind === 'atom' && h.index === i && R.pick(s.x, s.y) === i; });
    res.bonds = [...Array(12).keys()].every(k => { const s = scr(mid(k)), h = R.pickAny(s.x, s.y); return R.pickBond(s.x, s.y) === k && h.kind === 'bond' && h.index === k; });
    const V = BM.quat.mat(R.q), vz = [V[6], V[7], V[8]];
    const side = (k, d) => {
      const a = atomPos(R.bonds[2 * k]), b = atomPos(R.bonds[2 * k + 1]), u = [b[0] - a[0], b[1] - a[1], b[2] - a[2]];
      const w = [u[1] * vz[2] - u[2] * vz[1], u[2] * vz[0] - u[0] * vz[2], u[0] * vz[1] - u[1] * vz[0]], l = Math.hypot(...w), m = mid(k);
      const s = scr([m[0] + w[0] / l * d, m[1] + w[1] / l * d, m[2] + w[2] / l * d]);
      return R.pickBond(s.x, s.y);
    };
    res.sides = [side(1, 0.13), side(1, 0.2), side(0, 0.16), side(0, 0.22), side(6, -0.14)];
    res.empty = R.pickAny(5, 5);
    res.emptyBond = R.pickBond(5, 5);
    return res;
  });
  check(r.atoms, 'pickAny/pick hit every atom at its centre');
  check(r.bonds, 'pickBond/pickAny hit every bond at its midpoint');
  check(r.sides.join() === '1,-1,0,-1,6', 'pick radius = max(drawn, 0.15 Å): single 0.13 hit / 0.20 miss, double 0.16 hit / 0.22 miss: ' + r.sides);
  check(r.empty.kind === 'none' && r.empty.index === -1 && r.emptyBond === -1, 'empty space → none');
  r = await ev(() => {
    clearExtras();
    R.setOptions({ style: 'ballstick', atomScale: 0.1, bondRadius: 0.28 }); R.setOrientation([0, 0, 0, 1]);
    R.setScene([1, 1], Float64Array.from([0, 0, 0, 0, 0, -2]), Uint32Array.from([0, 1]), false);
    const w = R.canvas.clientWidth, h = R.canvas.clientHeight, out = { tie: R.pickAny(w / 2, h / 2), bond: R.pickBond(w / 2, h / 2) };
    R.setOptions({ atomScale: 0.25, bondRadius: 0.11 });
    // a bond 4 Å in front of atom 0, and atom 3 2 Å in front of a bond, along the pick rays
    const P = [[0, 0, -4], [-1.5, 0, 0], [1.5, 0, 0], [0, 2.5, 2], [-1.5, 2.5, 0], [1.5, 2.5, 0]], x = Float64Array.from(P.flat());
    R.setScene([6, 6, 6, 1, 6, 6], x, Uint32Array.from([1, 2, 4, 5]), false);
    const along = (m, dt) => { const s = R.worldToScreen(m), { o, d } = R.ray(s.x, s.y), t = Math.hypot(m[0] - o[0], m[1] - o[1], m[2] - o[2]) + dt; return { s, p: [o[0] + d[0] * t, o[1] + d[1] * t, o[2] + d[2] * t] }; };
    const A = along([0, 0, 0], 4), B = along([0, 2.5, 0], -2);
    x.set(A.p, 0); x.set(B.p, 9); R.updatePositions(x);
    Object.assign(out, { behind: R.pickAny(A.s.x, A.s.y), atomOnly: R.pick(A.s.x, A.s.y), bondOnly: R.pickBond(A.s.x, A.s.y), front: R.pickAny(B.s.x, B.s.y), bondBehind: R.pickBond(B.s.x, B.s.y) });
    return out;
  });
  check(r.tie.kind === 'atom' && r.tie.index === 0 && r.bond === 0, 'atoms win ties (identical atom sphere and bond capsule entry): ' + JSON.stringify(r.tie));
  check(r.behind.kind === 'bond' && r.behind.index === 0 && r.atomOnly === 0 && r.bondOnly === 0, 'nearer bond beats the atom behind it; pick() still finds the atom');
  check(r.front.kind === 'atom' && r.front.index === 3 && r.bondBehind === 1, 'nearer atom beats the bond behind it');

  console.log('atomsInRect, screen/world conversion');
  r = await ev(() => {
    load(benzene(false), { rot: RING_ROT, view: conj(RING_ROT) });
    const w = R.canvas.clientWidth, h = R.canvas.clientHeight, cx = R.worldToScreen(R.center).x;
    const out = { all: R.atomsInRect(0, 0, w, h).length, right: R.atomsInRect(w, h, cx, 0) };
    out.expRight = [...Array(12).keys()].filter(i => R.worldToScreen(atomPos(i)).x >= cx);
    const one = R.worldToScreen(atomPos(4)); out.tiny = R.atomsInRect(one.x - 2, one.y - 2, one.x + 2, one.y + 2);
    R.setOptions({ showH: false }); out.noH = R.atomsInRect(0, 0, w, h); R.setOptions({ showH: true });
    const V = BM.quat.mat(R.q), p = atomPos(2), s = R.worldToScreen(p), back = R.screenToWorld(s.x, s.y, p);
    out.err1 = Math.hypot(back[0] - p[0], back[1] - p[1], back[2] - p[2]); out.z = s.z;
    const w0 = R.screenToWorld(123, 456), c = R.center, s2 = R.worldToScreen(w0);
    out.plane = Math.abs((w0[0] - c[0]) * V[6] + (w0[1] - c[1]) * V[7] + (w0[2] - c[2]) * V[8]); out.err2 = Math.hypot(s2.x - 123, s2.y - 456);
    R.setOptions({ perspective: false });
    const s3 = R.worldToScreen(p), b3 = R.screenToWorld(s3.x, s3.y, p); out.err3 = Math.hypot(b3[0] - p[0], b3[1] - p[1], b3[2] - p[2]);
    out.orthoAll = R.atomsInRect(0, 0, w, h).length;
    R.setOptions({ perspective: true });
    return out;
  });
  check(r.all === 12 && r.right.length > 0 && r.right.length < 12 && r.right.join() === r.expRight.join(), 'atomsInRect (any corner order): ' + r.right);
  check(r.tiny.join() === '4' && r.noH.join() === '0,1,2,3,4,5' && r.orthoAll === 12, 'atomsInRect: small box, hidden H excluded, orthographic');
  check(r.err1 < 1e-9 && r.err2 < 1e-6 && r.plane < 1e-9 && r.err3 < 1e-9 && r.z < 0, 'screenToWorld/worldToScreen round trips: ' + [r.err1, r.err2, r.plane, r.err3].map(v => v.toExponential(1)));

  // ---------------- tool hooks (real pointer input) ----------------
  console.log('Tool hooks');
  await ev(() => {
    load(benzene(false), { rot: RING_ROT, view: conj(RING_ROT) });
    window.log = [];
    R.onAtomClick = i => log.push(['click', i]); R.onAtomDblClick = i => log.push(['dbl', i]); R.onBackgroundClick = () => log.push(['bg']);
    window.owner = {
      cursor: 'crosshair',
      onDown(e, hit) { log.push(['down', hit.kind, hit.index, hit.world, hit.px, hit.py]); return true; },
      onMove(e, hit) { log.push(['move', hit.kind, hit.index]); },
      onUp(e, hit) { log.push(['up', hit.kind, hit.index, e.type]); }
    };
    R.tool = owner;
  });
  const box = await page.locator('#gl').boundingBox();
  const at = async i => { const s = await ev(i => R.worldToScreen(atomPos(i)), i); return [box.x + s.x, box.y + s.y]; };
  const state = () => ev(() => ({ q: R.q.slice(), c: R.center.slice(), dist: R.dist, cursor: R.canvas.style.cursor, log: window.log.slice() }));
  const same = (a, b) => a.q.every((v, k) => v === b.q[k]) && a.c.every((v, k) => v === b.c[k]);
  const kinds = (log, k) => log.filter(e => e[0] === k);
  const drag = async (x, y, dx, dy, button) => { await page.mouse.move(x, y); await page.mouse.down({ button }); await page.mouse.move(x + dx, y + dy, { steps: 6 }); await page.mouse.up({ button }); };

  let s0 = await state(), s1;
  check(s0.cursor === 'crosshair', 'tool cursor applied');
  let [x0, y0] = await at(0);
  await drag(x0, y0, 80, 40);
  s1 = await state();
  const d0 = kinds(s1.log, 'down')[0];
  check(kinds(s1.log, 'down').length === 1 && d0[1] === 'atom' && d0[2] === 0 && near(d0[4], x0 - box.x, 1) && near(d0[5], y0 - box.y, 1), 'onDown gets the atom hit and CSS-pixel position');
  check(kinds(s1.log, 'move').length >= 3 && kinds(s1.log, 'up').length === 1 && kinds(s1.log, 'up')[0][3] === 'pointerup', 'onMove/onUp keep firing for the owned drag');
  check(same(s0, s1) && !s1.log.some(e => e[0] === 'click' || e[0] === 'bg'), 'owned drag: no rotate, no pan, no click');
  r = await ev(w => { const V = BM.quat.mat(R.q), p = atomPos(0); return { plane: Math.abs((w[0] - p[0]) * V[6] + (w[1] - p[1]) * V[7] + (w[2] - p[2]) * V[8]), dist: Math.hypot(w[0] - p[0], w[1] - p[1], w[2] - p[2]) }; }, d0[3]);
  check(r.plane < 1e-9 && r.dist < 0.05, 'hit.world lies on the screen-parallel plane through the hit atom: ' + JSON.stringify(r));

  await ev(() => { log.length = 0; });
  [x0, y0] = await at(3);
  await page.mouse.click(x0, y0); await page.mouse.click(x0, y0);
  s1 = await state();
  check(kinds(s1.log, 'down').length === 2 && kinds(s1.log, 'up').length === 2 && !s1.log.some(e => e[0] === 'click' || e[0] === 'dbl'), 'owned clicks never reach onAtomClick/onAtomDblClick');

  await ev(() => { log.length = 0; });
  s0 = await state();
  await drag(box.x + 400, box.y + 300, 60, 30, 'right');
  s1 = await state();
  check(s1.log.length === 0 && s0.q.every((v, k) => v === s1.q[k]) && !s0.c.every((v, k) => v === s1.c[k]), 'right-drag pans and never reaches the tool');
  await page.mouse.move(box.x + 400, box.y + 300); await page.mouse.wheel(0, 240); await page.waitForTimeout(80);
  s1 = await state();
  check(s1.dist > s0.dist * 1.2, 'wheel zoom with a tool active: ' + s0.dist.toFixed(2) + ' → ' + s1.dist.toFixed(2));

  await ev(() => { log.length = 0; R.tool = { onDown(e, hit) { log.push(['down', hit.kind]); return false; }, onMove() { log.push(['move']); }, onUp() { log.push(['up']); } }; });
  s0 = await state();
  [x0, y0] = await at(1);
  await drag(x0, y0, 60, 20);
  s1 = await state();
  check(!same(s0, s1) && kinds(s1.log, 'down').length === 1 && !kinds(s1.log, 'move').length && !kinds(s1.log, 'up').length, 'onDown → false: default rotation, no onMove/onUp');
  await ev(() => { R.setOrientation(conj(RING_ROT)); log.length = 0; });
  await page.waitForTimeout(400);
  [x0, y0] = await at(2);
  await page.mouse.click(x0, y0);
  s1 = await state();
  check(JSON.stringify(s1.log.filter(e => e[0] !== 'down')) === '[["click",2]]', 'onDown → false: click still reaches onAtomClick ' + JSON.stringify(s1.log));

  await ev(() => {
    log.length = 0;
    R.tool = { cursor: hit => (hit.kind === 'bond' ? 'pointer' : 'crosshair'), onHover(e, hit) { log.push(['hover', hit.kind, hit.index, hit.world]); R.highlightBond(hit.kind === 'bond' ? hit.index : -1); } };
  });
  const m0 = await ev(() => { const a = atomPos(R.bonds[0]), b = atomPos(R.bonds[1]); return R.worldToScreen([(a[0] + b[0]) / 2, (a[1] + b[1]) / 2, (a[2] + b[2]) / 2]); });
  await page.mouse.move(box.x + m0.x, box.y + m0.y);
  s1 = await state();
  let hv = kinds(s1.log, 'hover').pop();
  check(hv && hv[1] === 'bond' && hv[2] === 0 && s1.cursor === 'pointer' && (await ev(() => R.nBondHalo)) === 1, 'onHover reports the bond; cursor function; highlightBond');
  await shot('15_hover_bond_highlight');
  await page.mouse.move(box.x + 30, box.y + 30);
  s1 = await state(); hv = kinds(s1.log, 'hover').pop();
  r = await ev(w => { const V = BM.quat.mat(R.q), c = R.center; return Math.abs((w[0] - c[0]) * V[6] + (w[1] - c[1]) * V[7] + (w[2] - c[2]) * V[8]); }, hv[3]);
  check(hv[1] === 'none' && hv[2] === -1 && s1.cursor === 'crosshair' && r < 1e-9, 'hover on empty space: none, world on the plane through the view centre');
  await page.mouse.move(box.x + m0.x, box.y + m0.y); await ev(() => { log.length = 0; });
  await page.mouse.move(box.x + 810, box.y + 300);
  s1 = await state(); hv = kinds(s1.log, 'hover').pop();
  check(hv && hv[1] === 'none', 'pointerleave sends onHover(none)');

  await ev(() => { R.highlightBond(-1); R.tool = owner; log.length = 0; });
  s0 = await state();
  [x0, y0] = await at(5);
  await page.mouse.move(x0, y0); await page.mouse.down();
  await ev(() => { R.tool = null; });
  await page.mouse.move(x0 + 70, y0 + 30, { steps: 5 }); await page.mouse.up();
  s1 = await state();
  check(s1.log.length === 1 && s1.log[0][0] === 'down' && same(s0, s1) && s1.cursor === '', 'tool removed mid-drag: no further events, no rotation, cursor reset');
  await page.waitForTimeout(400);
  await ev(() => { log.length = 0; });
  [x0, y0] = await at(4);
  await page.mouse.click(x0, y0); await page.waitForTimeout(400); await page.mouse.click(box.x + 20, box.y + 20);
  s1 = await state();
  check(JSON.stringify(s1.log) === '[["click",4],["bg"]]', 'no tool: atom click and background click as before ' + JSON.stringify(s1.log));

  // ---------------- overlay hook ----------------
  console.log('Overlay hook');
  r = await ev(async () => {
    window.hookLog = [];
    const ov = document.getElementById('overlay');
    R.labelProvider = () => { hookLog.push('labels'); return [{ atom: 0, text: 'C1' }, { atom: 3, text: 'C4' }]; };
    R.overlayHook = (c, project) => {
      hookLog.push(c.canvas === ov ? 'hook' : 'hook-elsewhere');
      const a = project(atomPos(0)), b = project(atomPos(3));
      window.hookProj = a;
      c.fillStyle = 'rgba(124,178,255,0.14)'; c.strokeStyle = '#7cb2ff'; c.lineWidth = 1.5; c.setLineDash([5, 4]);
      const x = Math.min(a.x, b.x) - 30, y = Math.min(a.y, b.y) - 30, w = Math.abs(b.x - a.x) + 60, h = Math.abs(b.y - a.y) + 60;
      c.fillRect(x, y, w, h); c.strokeRect(x, y, w, h);
      c.setLineDash([]); c.strokeStyle = '#ffd04a'; c.beginPath(); c.moveTo(a.x, a.y); c.lineTo(b.x, b.y); c.stroke();
    };
    R.render();
    const order = hookLog.slice(), exp = R.worldToScreen(atomPos(0));
    hookLog.length = 0;
    const pending = R.snapshot(1), snap = hookLog.slice(), blob = await pending;
    return { order, projErr: Math.hypot(hookProj.x - exp.x, hookProj.y - exp.y), snap, blob: blob.size };
  });
  check(r.order.join() === 'labels,hook', 'overlay hook runs after labels, on the overlay canvas: ' + r.order);
  check(r.projErr < 1e-9, 'hook project() matches worldToScreen');
  check(r.snap.join() === 'labels,labels,hook' && r.blob > 1000, 'snapshot image skips the hook (only the follow-up screen render calls it): ' + r.snap);
  await shot('16_overlay_hook');
  await ev(() => clearExtras());

  // ---------------- per-atom representations ----------------
  console.log('Per-atom representations');
  r = await ev(() => {
    const S = BM.Renderer.STYLE, out = {};
    const readBuf = (buf, count, stride) => { const f = new Float32Array(count * stride); if (count) { R.gl.bindBuffer(R.gl.ARRAY_BUFFER, buf); R.gl.getBufferSubData(R.gl.ARRAY_BUFFER, 0, f); } return f; };
    const rep = (fn, extra) => { const n = R.Z.length, style = new Uint8Array(n); for (let i = 0; i < n; i++) style[i] = fn(i); return Object.assign({ style }, extra || {}); };
    load(benzene(false), { rot: RING_ROT, view: conj(RING_ROT) });
    const base = [R.nSphere, R.nMainCyl, R.nMultiCyl, R.nCyl];
    // carbons ball-and-stick, hydrogens lines: each bond half in the style of its atom
    R.setRepresentation(rep(i => (R.Z[i] === 6 ? S.ballstick : S.lines)));
    out.mixed = [R.nSphere, R.nMainCyl, R.nMultiCyl, R.nLine, R.nGhostSph];
    out.visH = R.visible(6); out.radH = R.atomRadius(6);
    // everything as lines: one instance per bond with both colours
    R.setRepresentation(rep(() => S.lines));
    out.lines = [R.nSphere, R.nMainCyl, R.nMultiCyl, R.nLine];
    const L = readBuf(R.lineInst, R.nLine, 12);
    const ch = R.bonds.length / 2 - 1; // last bond: C5-H11
    out.lineCols = [Array.from(L.subarray(12 * ch + 6, 12 * ch + 9)).map(v => +v.toFixed(3)), Array.from(L.subarray(12 * ch + 9, 12 * ch + 12)).map(v => +v.toFixed(3))];
    out.elemC = BM.el.color[6].map(v => +v.toFixed(3)); out.elemH = BM.el.color[1].map(v => +v.toFixed(3));
    // translucent lines fade towards the background colour
    R.setRepresentation(rep(() => S.lines, { alpha: new Float32Array(12).fill(0.4) }));
    const L2 = readBuf(R.lineInst, R.nLine, 12), bg = R.opt.bg;
    out.dim = Math.abs(L2[6] - (bg[0] + (BM.el.color[6][0] - bg[0]) * 0.4)) < 1e-5;
    // ghosts: separate translucent buffers with per-instance alpha (default 0.3)
    R.setRepresentation(rep(i => (R.Z[i] === 1 ? S.ghost : S.sticks)));
    const G = readBuf(R.ghostSphInst, R.nGhostSph, 8), GC = readBuf(R.ghostCylInst, R.nGhostCyl, 11);
    out.ghost = [R.nSphere, R.nGhostSph, R.nMainCyl, R.nGhostCyl, +G[7].toFixed(3), +GC[10].toFixed(3)];
    R.setRepresentation(rep(i => (R.Z[i] === 1 ? S.ghost : S.sticks), { alpha: new Float32Array(12).fill(0.55) }));
    out.ghostAlpha = +readBuf(R.ghostSphInst, 1, 8)[7].toFixed(3);
    // radius scale
    R.setRepresentation(rep(() => S.ballstick, { scale: new Float32Array(12).fill(2) }));
    out.scaled = +(R.atomRadius(0) / (BM.el.vdw[6] * R.opt.atomScale)).toFixed(4);
    // hidden atoms are neither drawn nor picked; proxies make them pickable again
    R.setRepresentation(rep(i => (R.Z[i] === 1 ? S.hidden : S.ballstick)));
    const hp = R.worldToScreen(atomPos(6)), cp = R.worldToScreen(atomPos(0));
    out.hidden = [R.nSphere, R.visible(6), R.pick(hp.x, hp.y), R.pick(cp.x, cp.y)];
    R.setRepresentation(rep(() => S.hidden));
    out.allHidden = [R.nSphere, R.pick(cp.x, cp.y)];
    R.setPickProxies([0], 1.0); R.setSelection([0]);
    out.proxy = [R.pick(cp.x, cp.y), R.pick(hp.x, hp.y), R.nHalo];
    R.setSelection([]); R.setPickProxies(null);
    // extent: hidden atoms that count for fitting (cartoon atoms)
    R._updateRadius(); out.radNone = R.radius;
    R.setRepresentation(rep(() => S.hidden, { extent: new Uint8Array(12).fill(1) })); R._updateRadius(); out.radExt = +R.radius.toFixed(3);
    R.setRepresentation(null); R._updateRadius(); out.radAll = +R.radius.toFixed(3);
    // null restores the toolbar style exactly
    out.restored = [R.nSphere, R.nMainCyl, R.nMultiCyl, R.nCyl].join() === base.join();
    // markers from representations add to the frozen-atom markers
    R.setAtomMarkers([2]); R.setRepMarkers([{ indices: [0, 1], color: [1, 0, 1] }]); out.markers = R.nMarker; R.setAtomMarkers([]); R.setRepMarkers([]);
    return out;
  });
  check(r.mixed.join() === '6,12,12,6,0', 'C ball-and-stick + H lines: 6 spheres, 12 cylinder halves, 12 double-bond lines, 6 line halves ' + JSON.stringify(r.mixed));
  check(r.visH && r.radH === 0, 'line-style atoms are visible with no sphere');
  check(r.lines.join() === '0,0,0,12', 'all lines: one line instance per bond, no spheres or cylinders ' + JSON.stringify(r.lines));
  check(JSON.stringify(r.lineCols) === JSON.stringify([r.elemC, r.elemH]), 'line colour split at the bond middle (C grey / H white)');
  check(r.dim, 'translucent lines fade towards the background');
  check(r.ghost.join() === '6,6,18,6,0.3,0.3', 'ghost H: translucent spheres and bond halves (alpha 0.3), C sticks 18 halves ' + JSON.stringify(r.ghost));
  check(r.ghostAlpha === 0.55, 'ghost alpha from the representation');
  check(r.scaled === 2, 'radius scale');
  check(r.hidden.join() === '6,false,-1,0', 'hidden H: not drawn, not pickable; C picked ' + JSON.stringify(r.hidden));
  check(r.allHidden.join() === '0,-1' && r.proxy.join() === '0,-1,1', 'proxy makes a hidden atom pickable with a halo ' + JSON.stringify(r.proxy));
  check(r.radNone === 1 && Math.abs(r.radExt - r.radAll) < 1e-6, 'extent atoms count for the view radius (' + r.radExt + ' vs ' + r.radAll + ')');
  check(r.restored, 'setRepresentation(null) restores the toolbar-style buffers');
  check(r.markers === 3, 'representation markers + frozen markers');
  await shot('17_reps_mixed_lines');

  // overlays, legend, level of detail
  r = await ev(() => {
    const S = BM.Renderer.STYLE, out = {};
    load(benzene(false), { rot: RING_ROT, view: conj(RING_ROT) });
    const n = R.Z.length, style = new Uint8Array(n);
    R.setRepresentation({ style });
    const px = (x, y) => { R.render(); const gl = R.gl, d = new Uint8Array(4), dpr = R.dpr; gl.readPixels(Math.round(x * dpr), Math.round(R.canvas.height - y * dpr), 1, 1, gl.RGBA, gl.UNSIGNED_BYTE, d); return Array.from(d); };
    const c = R.worldToScreen([0, 0, 0]);
    out.bg = px(c.x, c.y);
    // a red two-sided triangle through the ring centre
    const tri = s => ({ positions: new Float32Array([-s, -s, 0, s, -s, 0, 0, s, 0]), normals: new Float32Array([0, 0, 1, 0, 0, 1, 0, 0, 1]), colors: new Float32Array([1, 0, 0, 1, 0, 0, 1, 0, 0]), indices: new Uint32Array([0, 1, 2]), twoSided: true });
    const M = BM.quat.mat(RING_ROT), rot = p => [M[0] * p[0] + M[1] * p[1] + M[2] * p[2], M[3] * p[0] + M[4] * p[1] + M[5] * p[2], M[6] * p[0] + M[7] * p[1] + M[8] * p[2]];
    const m = tri(2.5); for (let v = 0; v < 3; v++) { const q = rot([m.positions[3 * v], m.positions[3 * v + 1], m.positions[3 * v + 2]]); m.positions.set(q, 3 * v); const nn = rot([0, 0, 1]); m.normals.set(nn, 3 * v); }
    R.setOverlays([m]);
    out.opaque = px(c.x, c.y); out.nOverlay = R.overlays.length;
    R.setOverlays([Object.assign({}, m, { opacity: 0.5 })]);
    out.half = px(c.x, c.y);
    R.setOverlays([]);
    // legend box in the top-right corner of the overlay canvas, also in snapshots
    R.setLegend({ title: 'ONIOM layers', items: [{ label: 'High', sub: 'B3LYP/6-31G(d) · 3 atoms', glyph: 'ballstick' }, { label: 'Low', sub: 'UFF', glyph: 'lines' }] });
    R.render();
    const ctx = R.overlay.getContext('2d'), w = R.overlay.width, a = ctx.getImageData(w - 30, 20, 1, 1).data[3], b = ctx.getImageData(10, 20, 1, 1).data[3];
    out.legend = [a > 0, b === 0];
    R.setLegend(null);
    // level of detail follows the drawn instances; wire tubes become lines for huge bond counts
    const many = k => { const Z = new Array(k).fill(6), xyz = new Float64Array(3 * k); for (let i = 0; i < k; i++) { xyz[3 * i] = (i % 100) * 3; xyz[3 * i + 1] = Math.floor(i / 100) % 100 * 3; xyz[3 * i + 2] = Math.floor(i / 10000) * 3; } return { Z, xyz }; };
    R.setRepresentation(null, true);
    let s = many(25000); R.setScene(s.Z, s.xyz, new Uint32Array(0)); out.lod25k = R._lod;
    s = many(4000); R.setScene(s.Z, s.xyz, new Uint32Array(0)); out.lod4k = R._lod;
    s = many(25000); R.setRepresentation({ style: new Uint8Array(25000).fill(S.lines) }, true); R.setScene(s.Z, s.xyz, new Uint32Array(0)); out.lodLines = [R._lod, R.nSphere, R.nLine];
    const k = 70000, Z = new Array(k).fill(6), xyz = new Float64Array(3 * k), bonds = new Uint32Array(2 * (k - 1));
    for (let i = 0; i < k; i++) { xyz[3 * i] = (i % 2) * 0.9; xyz[3 * i + 1] = i * 1.2; }
    for (let i = 0; i + 1 < k; i++) { bonds[2 * i] = i; bonds[2 * i + 1] = i + 1; }
    R.setRepresentation(null, true); R.setOptions({ style: 'wireframe' }); R.setScene(Z, xyz, bonds);
    out.wireBig = [R.nMainCyl, R.nLine, R._wireLines];
    R.setOptions({ style: 'ballstick' }); out.cylLod = R._cylLod;
    return out;
  });
  check(r.opaque[0] > 200 && r.opaque[1] < 60 && r.nOverlay === 1, 'opaque overlay mesh drawn (red at the ring centre): ' + r.opaque + ' over ' + r.bg);
  check(r.half[0] > r.bg[0] + 40 && r.half[0] < r.opaque[0] - 40, 'translucent overlay blended: ' + r.half);
  check(r.legend[0] && r.legend[1], 'legend drawn in the top-right corner only');
  check(r.lod25k === 1 && r.lod4k === 2 && r.lodLines.join() === '3,0,75000', 'sphere detail by drawn spheres (lines need none): ' + [r.lod25k, r.lod4k, r.lodLines]);
  check(r.wireBig[0] === 0 && r.wireBig[1] === 69999 && r.wireBig[2] === true && r.cylLod === 10, 'huge wire models use lines; cylinder detail drops above 40k halves ' + JSON.stringify(r.wireBig) + ' ' + r.cylLod);
  await ev(() => { clearExtras(); R.setOptions({ style: 'ballstick' }); });

  console.log('Console/page errors:', errors.length ? errors : 'none');
  check(errors.length === 0, 'no console errors');
  await browser.close();
  console.log(`\nRenderer tests: ${passes} passed, ${failures} failed  (screenshots in ${SHOTS})`);
  process.exit(failures ? 1 : 0);
})().catch(e => { console.error(e); process.exit(2); });
