/* BlueMolecule — crystal tools: lattice parameters, supercells, surface slabs, cell edits (BM.crystal).
 *
 *   BM.crystal.params(cell)                     -> {a, b, c, alpha, beta, gamma, volume} (Å, degrees, Å³)
 *   BM.crystal.fromParams(a, b, c, al, be, ga)   -> cell [[ax ay az], [bx by bz], [cx cy cz]], a along x, b in xy
 *   BM.crystal.toFrac(cell, xyz) / toCart(cell, frac)
 *   BM.crystal.wrap(cell, xyz)                   -> Float64Array: every atom moved into the cell (fractions in [0, 1))
 *   BM.crystal.supercell(atoms, cell, n)         -> {Z, xyz, cell, map}: n = [n1, n2, n3] or a 3×3 integer matrix
 *   BM.crystal.slab(atoms, cell, hkl, opts)      -> {Z, xyz, cell, map, d, normal}: opts {layers, vacuum, shift}
 *   BM.crystal.setCell(atoms, cell, newCell, scale) -> Float64Array: positions in the new cell (scaled or kept)
 *   BM.crystal.box(atoms, margin)                -> {xyz, cell}: an orthorhombic box around a molecule, centred
 *   atoms: {Z, xyz (Å, flat)}; map[k]: the source atom of new atom k (to carry names along)
 *
 * The slab: the Miller indices (h k l) of the cell's lattice give two lattice vectors in the plane and a third with
 * h u + k v + l w = 1, a unimodular change of basis (same volume, same atoms). The oriented cell is repeated
 * `layers` times along the third vector; the slab's cell keeps the two in-plane vectors and has its third vector along
 * the normal, one layer spacing per layer plus the vacuum. The result is rotated so that the first vector lies along
 * x and the normal along z.
 */
(function (BM) {
  'use strict';
  const RAD = Math.PI / 180;
  const dot = (a, b) => a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
  const cross = (a, b) => [a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]];
  const norm = a => Math.hypot(a[0], a[1], a[2]);
  const scale = (a, s) => [a[0] * s, a[1] * s, a[2] * s];
  const add = (a, b) => [a[0] + b[0], a[1] + b[1], a[2] + b[2]];
  const det3 = m => dot(m[0], cross(m[1], m[2]));

  function check(cell) {
    if (!Array.isArray(cell) || cell.length !== 3 || cell.some(v => !v || v.length !== 3 || !v.every(Number.isFinite))) throw new Error('There is no unit cell.');
    if (Math.abs(det3(cell)) < 1e-8) throw new Error('The unit cell has no volume.');
  }
  function params(cell) {
    check(cell);
    const [a, b, c] = cell, la = norm(a), lb = norm(b), lc = norm(c);
    const ang = (u, v, lu, lv) => Math.acos(Math.max(-1, Math.min(1, dot(u, v) / (lu * lv)))) / RAD;
    return { a: la, b: lb, c: lc, alpha: ang(b, c, lb, lc), beta: ang(a, c, la, lc), gamma: ang(a, b, la, lb), volume: Math.abs(det3(cell)) };
  }
  function fromParams(a, b, c, alpha, beta, gamma) {
    if (![a, b, c].every(x => x > 0)) throw new Error('Cell lengths must be positive.');
    const ca = Math.cos(alpha * RAD), cb = Math.cos(beta * RAD), cg = Math.cos(gamma * RAD), sg = Math.sin(gamma * RAD);
    if (!(sg > 1e-8)) throw new Error('γ must be between 0° and 180°.');
    const cx = cb, cy = (ca - cb * cg) / sg, cz2 = 1 - cx * cx - cy * cy;
    if (!(cz2 > 1e-10)) throw new Error('These angles do not make a cell (α, β and γ are incompatible).');
    return [[a, 0, 0], [b * cg, b * sg, 0], [c * cx, c * cy, c * Math.sqrt(cz2)]];
  }
  /** Inverse of the matrix whose rows are the cell vectors (r = f · cell, so f = r · inv). */
  function inverse(cell) {
    const [a, b, c] = cell, d = det3(cell), bc = cross(b, c), ca = cross(c, a), ab = cross(a, b);
    // columns of the inverse are (b × c, c × a, a × b) / det
    return [[bc[0] / d, ca[0] / d, ab[0] / d], [bc[1] / d, ca[1] / d, ab[1] / d], [bc[2] / d, ca[2] / d, ab[2] / d]];
  }
  function toFrac(cell, xyz) {
    check(cell);
    const M = inverse(cell), n = xyz.length / 3, out = new Float64Array(3 * n);
    for (let i = 0; i < n; i++) {
      const x = xyz[3 * i], y = xyz[3 * i + 1], z = xyz[3 * i + 2];
      for (let k = 0; k < 3; k++) out[3 * i + k] = x * M[0][k] + y * M[1][k] + z * M[2][k];
    }
    return out;
  }
  function toCart(cell, f) {
    const n = f.length / 3, out = new Float64Array(3 * n);
    for (let i = 0; i < n; i++) for (let k = 0; k < 3; k++) out[3 * i + k] = f[3 * i] * cell[0][k] + f[3 * i + 1] * cell[1][k] + f[3 * i + 2] * cell[2][k];
    return out;
  }
  const into01 = v => { let w = v - Math.floor(v); if (w >= 1 - 1e-9) w = 0; return w; };
  function wrap(cell, xyz) {
    const f = toFrac(cell, xyz);
    for (let i = 0; i < f.length; i++) f[i] = into01(f[i]);
    return toCart(cell, f);
  }

  function supercell(atoms, cell, n) {
    check(cell);
    let M;
    if (Array.isArray(n) && n.length === 3 && n.every(x => typeof x === 'number')) M = [[n[0], 0, 0], [0, n[1], 0], [0, 0, n[2]]];
    else M = n;
    if (!Array.isArray(M) || M.length !== 3 || M.some(r => !r || r.length !== 3 || r.some(x => !Number.isInteger(x)))) throw new Error('A supercell needs whole numbers.');
    const dM = det3(M);
    if (!(dM >= 1)) throw new Error('The supercell matrix must have a positive determinant.');
    if (dM > 1e5) throw new Error('That supercell is too large.');
    const nat = atoms.Z.length;
    if (nat * dM > 2e6) throw new Error('That supercell would have ' + Math.round(nat * dM).toLocaleString('en') + ' atoms (at most 2,000,000).');
    const newCell = M.map(r => add(add(scale(cell[0], r[0]), scale(cell[1], r[1])), scale(cell[2], r[2])));
    const Z = [], xyz = [], map = [];
    if (M.every((r, i) => r.every((x, j) => (i === j ? x >= 1 : x === 0)))) {
      // diagonal: every translation of the box, atoms as they are (molecules stay whole)
      for (let i = 0; i < M[0][0]; i++) for (let j = 0; j < M[1][1]; j++) for (let k = 0; k < M[2][2]; k++) {
        const t = add(add(scale(cell[0], i), scale(cell[1], j)), scale(cell[2], k));
        for (let p = 0; p < nat; p++) { Z.push(atoms.Z[p]); xyz.push(atoms.xyz[3 * p] + t[0], atoms.xyz[3 * p + 1] + t[1], atoms.xyz[3 * p + 2] + t[2]); map.push(p); }
      }
      return { Z, xyz: Float64Array.from(xyz), cell: newCell, map };
    }
    // general matrix: the lattice translations inside the new cell
    return fill(atoms, cell, newCell, dM);
  }
  /** Atoms of the lattice (cell + atoms) that fall inside newCell, wrapped; expects |newCell| / |cell| = count. */
  function fill(atoms, cell, newCell, count) {
    const nat = atoms.Z.length, f0 = toFrac(cell, atoms.xyz), Minv = inverse(newCell);
    // the integer translations to try: the corners of newCell in old fractional coordinates
    const corners = [];
    for (const a of [0, 1]) for (const b of [0, 1]) for (const c of [0, 1]) corners.push(toFrac(cell, add(add(scale(newCell[0], a), scale(newCell[1], b)), scale(newCell[2], c))));
    const lo = [0, 1, 2].map(k => Math.floor(Math.min(...corners.map(v => v[k]))) - 1), hi = [0, 1, 2].map(k => Math.ceil(Math.max(...corners.map(v => v[k]))) + 1);
    const Z = [], xyz = [], map = [], seen = new Set();
    for (let i = lo[0]; i <= hi[0]; i++) for (let j = lo[1]; j <= hi[1]; j++) for (let k = lo[2]; k <= hi[2]; k++) {
      for (let p = 0; p < nat; p++) {
        const fx = f0[3 * p] + i, fy = f0[3 * p + 1] + j, fz = f0[3 * p + 2] + k;
        const r = [fx * cell[0][0] + fy * cell[1][0] + fz * cell[2][0], fx * cell[0][1] + fy * cell[1][1] + fz * cell[2][1], fx * cell[0][2] + fy * cell[1][2] + fz * cell[2][2]];
        const g = [0, 1, 2].map(q => r[0] * Minv[0][q] + r[1] * Minv[1][q] + r[2] * Minv[2][q]);
        if (g.some(v => v < -1e-7 || v >= 1 - 1e-7)) continue;
        const key = p + ':' + g.map(v => Math.round(into01(v) * 1e5)).join(',');
        if (seen.has(key)) continue;
        seen.add(key);
        Z.push(atoms.Z[p]); xyz.push(r[0], r[1], r[2]); map.push(p);
      }
    }
    if (Math.abs(Z.length - nat * count) > 0) throw new Error('The new cell holds ' + Z.length + ' atoms instead of ' + Math.round(nat * count) + '; the structure may not be periodic in this cell.');
    return { Z, xyz: Float64Array.from(xyz), cell: newCell, map };
  }

  function gcd(a, b) { a = Math.abs(a); b = Math.abs(b); while (b) { [a, b] = [b, a % b]; } return a; }
  /** [g, x, y] with a x + b y = g = gcd(a, b). */
  function egcd(a, b) {
    if (b === 0) return [Math.abs(a), a < 0 ? -1 : 1, 0];
    const [g, x, y] = egcd(b, a % b);
    return [g, y, x - Math.trunc(a / b) * y];
  }
  /** Integer vectors (in lattice coordinates) c1, c2 spanning the (hkl) plane and c3 with h·c3 = 1. */
  function planeBasis(cell, h, k, l) {
    let c1, c2, c3;
    const g = gcd(k, l);
    if (g === 0) {   // (h 0 0)
      c1 = [0, 1, 0]; c2 = [0, 0, 1]; c3 = [h > 0 ? 1 : -1, 0, 0];
    } else {
      const [, p, q] = egcd(k, l);                 // p k + q l = g
      c1 = [g, -p * h, -q * h];
      c2 = [0, l / g, -k / g];
      const [, a, b] = egcd(g, h);                 // a g + b h = 1 (gcd(h, k, l) = 1)
      c3 = [b, a * p, a * q];
    }
    const cart = v => add(add(scale(cell[0], v[0]), scale(cell[1], v[1])), scale(cell[2], v[2]));
    // the most orthogonal in-plane pair (reduce c1 against c2), then c3 as close to the normal as possible
    const v2 = cart(c2), t = Math.round(dot(cart(c1), v2) / dot(v2, v2));
    if (t) c1 = c1.map((x, i) => x - t * c2[i]);
    const V1 = cart(c1), V2 = cart(c2), nrm = cross(V1, V2), n2 = dot(nrm, nrm);
    // remove the in-plane part of c3 in units of c1 and c2 (least squares in the plane)
    const w = cart(c3), wp = add(w, scale(nrm, -dot(w, nrm) / n2));
    const g11 = dot(V1, V1), g12 = dot(V1, V2), g22 = dot(V2, V2), b1 = dot(wp, V1), b2 = dot(wp, V2), D = g11 * g22 - g12 * g12;
    const m1 = Math.round((b1 * g22 - b2 * g12) / D), m2 = Math.round((b2 * g11 - b1 * g12) / D);
    c3 = c3.map((x, i) => x - m1 * c1[i] - m2 * c2[i]);
    // right-handed
    if (det3([cart(c1), cart(c2), cart(c3)]) < 0) { const tmp = c1; c1 = c2; c2 = tmp; }
    return [c1, c2, c3];
  }
  function slab(atoms, cell, hkl, opts) {
    check(cell);
    opts = opts || {};
    let [h, k, l] = hkl.map(x => Math.round(+x));
    if (![h, k, l].every(Number.isFinite) || (!h && !k && !l)) throw new Error('Miller indices cannot all be zero.');
    const g0 = gcd(gcd(h, k), l);
    h /= g0; k /= g0; l /= g0;
    const layers = Math.max(1, Math.round(opts.layers || 1)), vacuum = Math.max(0, +opts.vacuum || 0);
    const C = planeBasis(cell, h, k, l);
    const oriented = C.map(v => add(add(scale(cell[0], v[0]), scale(cell[1], v[1])), scale(cell[2], v[2])));
    // the atoms in the oriented cell (same volume), with the termination shift along the third vector
    const one = fill(atoms, cell, oriented, 1);
    const nat = one.Z.length, f = toFrac(oriented, one.xyz), shift = +opts.shift || 0;
    for (let p = 0; p < nat; p++) { f[3 * p] = into01(f[3 * p]); f[3 * p + 1] = into01(f[3 * p + 1]); f[3 * p + 2] = into01(f[3 * p + 2] - shift); }
    const [A1, A2, A3] = oriented, nrm0 = cross(A1, A2), nhat = scale(nrm0, 1 / norm(nrm0)), d = dot(A3, nhat);
    // rotate: A1 along x, the normal along z
    const e1 = scale(A1, 1 / norm(A1)), e3 = nhat, e2 = cross(e3, e1);
    const rot = v => [dot(v, e1), dot(v, e2), dot(v, e3)];
    const R1 = rot(A1), R2 = rot(A2), R3 = rot(A3);
    const Z = [], xyz = [], map = [];
    for (let L = 0; L < layers; L++) {
      for (let p = 0; p < nat; p++) {
        const fz = f[3 * p + 2] + L;
        let x = f[3 * p] * R1[0] + f[3 * p + 1] * R2[0] + fz * R3[0];
        let y = f[3 * p] * R1[1] + f[3 * p + 1] * R2[1] + fz * R3[1];
        const z = f[3 * p] * R1[2] + f[3 * p + 1] * R2[2] + fz * R3[2];
        // back into the in-plane cell (the third vector's slant moves each layer sideways)
        const u = R2[1] ? into01((x - y * R2[0] / R2[1]) / R1[0]) : into01(x / R1[0]), v = R2[1] ? into01(y / R2[1]) : 0;
        x = u * R1[0] + v * R2[0]; y = v * R2[1];
        Z.push(one.Z[p]); xyz.push(x, y, z); map.push(one.map[p]);
      }
    }
    const thick = layers * d;
    const out = { Z, xyz: Float64Array.from(xyz), cell: [R1, R2, [0, 0, thick + vacuum]], map, d, normal: [0, 0, 1], hkl: [h, k, l] };
    if (vacuum > 0) {   // the slab in the middle of the cell
      let zlo = Infinity, zhi = -Infinity;
      for (let p = 0; p < Z.length; p++) { zlo = Math.min(zlo, out.xyz[3 * p + 2]); zhi = Math.max(zhi, out.xyz[3 * p + 2]); }
      const dz = (thick + vacuum) / 2 - (zlo + zhi) / 2;
      for (let p = 0; p < Z.length; p++) out.xyz[3 * p + 2] += dz;
    }
    return out;
  }

  /** Positions in a new cell: scale = true keeps the fractional coordinates (the atoms move with the cell),
   *  false the Cartesian ones. The new cell is turned to the old cell's orientation. */
  function oriented(cell, newCell) {
    // the rotation that takes the standard orientation of the old parameters onto the old cell, applied to the new
    const p = params(cell), S0 = fromParams(p.a, p.b, p.c, p.alpha, p.beta, p.gamma);
    const Sinv = inverse(S0);   // r = f · S0 -> f = r · Sinv; the old cell = S0 · T with T = Sinv · cell
    const T = [0, 1, 2].map(i => [0, 1, 2].map(j => Sinv[i][0] * cell[0][j] + Sinv[i][1] * cell[1][j] + Sinv[i][2] * cell[2][j]));
    return newCell.map(v => [0, 1, 2].map(j => v[0] * T[0][j] + v[1] * T[1][j] + v[2] * T[2][j]));
  }
  function setCell(atoms, cell, newCell, scaleAtoms) {
    check(cell); check(newCell);
    if (!scaleAtoms) return Float64Array.from(atoms.xyz);
    return toCart(newCell, toFrac(cell, atoms.xyz));
  }
  function box(atoms, margin) {
    const n = atoms.Z.length, m = margin > 0 ? margin : 10;
    if (!n) throw new Error('There are no atoms.');
    const lo = [Infinity, Infinity, Infinity], hi = [-Infinity, -Infinity, -Infinity];
    for (let i = 0; i < n; i++) for (let a = 0; a < 3; a++) { lo[a] = Math.min(lo[a], atoms.xyz[3 * i + a]); hi[a] = Math.max(hi[a], atoms.xyz[3 * i + a]); }
    const L = [0, 1, 2].map(a => hi[a] - lo[a] + 2 * m), xyz = new Float64Array(3 * n);
    for (let i = 0; i < n; i++) for (let a = 0; a < 3; a++) xyz[3 * i + a] = atoms.xyz[3 * i + a] - (lo[a] + hi[a]) / 2 + L[a] / 2;
    return { xyz, cell: [[L[0], 0, 0], [0, L[1], 0], [0, 0, L[2]]] };
  }

  BM.crystal = { params, fromParams, toFrac, toCart, wrap, supercell, slab, setCell, oriented, box, planeBasis, inverse };
})(typeof window !== 'undefined' ? (window.BM = window.BM || {}) : (globalThis.BM = globalThis.BM || {}));
