/* BlueMolecule — molecular point-group symmetry: detection, operations and symmetrization.
 *
 * BM.sym.detect(Z, xyz, {tol = 0.1}) -> {
 *   pointGroup     Schoenflies symbol: C1 Cs Ci Cn Cnv Cnh Dn Dnh Dnd S2n T Td Th O Oh I Ih Cinfv Dinfh
 *                  ('Kh' for a single atom)
 *   order          number of operations (Infinity for linear molecules and single atoms)
 *   operations     [{type: 'E'|'C'|'S'|'sigma'|'i', n, k, axis, normal, label, matrix, perm}]
 *                    C/S: rotation by 2*pi*k/n about `axis` (k <= n/2; S = rotation then reflection
 *                    through the perpendicular plane); sigma: mirror plane with unit `normal`.
 *                    matrix: row-major 3x3 acting on positions relative to `center`.
 *                    perm[i] = j: the operation carries atom i onto atom j.
 *                  Linear molecules and atoms list a finite subgroup (C2v / D2h / Ci), which is
 *                  enough to symmetrize them.
 *   center         centre of mass [x, y, z] (Å)
 *   principalAxis  unit vector of the highest-order proper axis (Cs: mirror normal; C1, Ci: null)
 *   frame          {x, y, z}: orthonormal symmetry frame (z = principal axis for axial groups,
 *                  the C2/C4 axes for cubic groups)
 *   topType        'atom' | 'linear' | 'spherical' | 'symmetric' | 'asymmetric'
 *   moments        principal moments of inertia (amu Å², ascending)
 *   classes        sets of symmetry-equivalent atoms (arrays of atom indices)
 *   tol            tolerance actually used: at most 0.45 x the median distance between neighbouring
 *                  atoms of the same element, and halved when the group found at `tol` is inconsistent
 *   maxDeviation   largest distance between an operated atom and its image atom (Å)
 * }
 * `tol` is the largest allowed distance between an operated atom and the atom it lands on, so
 * random displacements up to about tol / 2 per atom keep the group.
 *
 * BM.sym.symmetrize(Z, xyz, {tol = 0.1, pointGroup}) -> {xyz, pointGroup, order, maxShift, rmsShift}
 *   Averages every atom over the operations of the group so that the result is exactly symmetric,
 *   keeping the centre of mass and the orientation. With `pointGroup`, the structure is made to have
 *   that group: a subgroup of the detected group, or a larger group found by relaxing the tolerance
 *   (up to 8 x tol). Throws if the structure is not close to the requested group.
 *
 * BM.sym.normalizeName(s) -> canonical Schoenflies symbol or null ('c2v' -> 'C2v', 'D*h' -> 'Dinfh')
 * BM.sym.groupOrder(name) -> number of operations of a point group
 *
 * Method: centre of mass, inertia tensor (linear test), symmetry-equivalent atom sets from
 * distance signatures, candidate axes (principal axes, atom directions, midpoints of equivalent
 * pairs, normals of planes through the centre and atom pairs, normals of planes through three
 * equivalent atoms) and candidate mirror normals (principal axes, equivalent-pair differences,
 * plane normals). Each candidate is tested by mapping every atom onto an atom of the same element;
 * accepted elements are refined by least squares from the atom permutation, and products of found
 * rotations supply further axis candidates. The group is assigned with the standard flowchart,
 * rebuilt exactly in a symmetry frame from its generators and closed; each operation's atom
 * permutation is the composition of its generators' permutations, and every operation is verified
 * against the structure. Typical cost: < 0.1 s for 200 atoms.
 */
(function (BM) {
  'use strict';

  const TWO_PI = 2 * Math.PI;
  const PHI = (1 + Math.sqrt(5)) / 2;
  const SIN_T = Math.sin(0.1), COS_T = Math.cos(0.1);   // angular slack when relating symmetry elements
  const COS_SKIP = Math.cos(0.1);                        // candidate already covered by a found element
  const CELL = 0.03;                                     // hash cell for candidate-direction de-duplication
  const WORK_LIMIT = 6000;                               // pair candidates generated before stopping early
  const TEST_BUDGET = 4e7;                               // atom-distance evaluations per element search
  const TOL_CAP = 0.45;                                  // tol <= TOL_CAP x closest same-element distance

  // ------------------------------------------------------------------ vectors, 3x3 matrices (row-major)
  const dot = (a, b) => a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
  const cross = (a, b) => [a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]];
  const vnorm = a => Math.sqrt(dot(a, a));
  const vscale = (a, s) => [a[0] * s, a[1] * s, a[2] * s];
  const vadd = (a, b) => [a[0] + b[0], a[1] + b[1], a[2] + b[2]];
  const vsub = (a, b) => [a[0] - b[0], a[1] - b[1], a[2] - b[2]];
  function unit(a) { const l = vnorm(a); return l > 0 ? [a[0] / l, a[1] / l, a[2] / l] : [0, 0, 1]; }
  /** Sign convention for directions: the first clearly non-zero component is positive. */
  function canon(u) {
    for (let k = 0; k < 3; k++) {
      if (Math.abs(u[k]) > 1e-6) return u[k] < 0 ? [-u[0], -u[1], -u[2]] : u.slice();
    }
    return u.slice();
  }
  function perpendicular(a) {
    const t = Math.abs(a[0]) < 0.6 ? [1, 0, 0] : Math.abs(a[1]) < 0.6 ? [0, 1, 0] : [0, 0, 1];
    return unit(cross(a, t));
  }
  const IDENT = [1, 0, 0, 0, 1, 0, 0, 0, 1];
  const INVERSION = [-1, 0, 0, 0, -1, 0, 0, 0, -1];
  const diag = (a, b, c) => [a, 0, 0, 0, b, 0, 0, 0, c];
  function mmul(A, B) {
    const C = new Array(9);
    for (let i = 0; i < 3; i++) {
      for (let j = 0; j < 3; j++) C[3 * i + j] = A[3 * i] * B[j] + A[3 * i + 1] * B[3 + j] + A[3 * i + 2] * B[6 + j];
    }
    return C;
  }
  const mvec = (M, v) => [M[0] * v[0] + M[1] * v[1] + M[2] * v[2], M[3] * v[0] + M[4] * v[1] + M[5] * v[2],
    M[6] * v[0] + M[7] * v[1] + M[8] * v[2]];
  const mT = M => [M[0], M[3], M[6], M[1], M[4], M[7], M[2], M[5], M[8]];
  const det3 = M => M[0] * (M[4] * M[8] - M[5] * M[7]) - M[1] * (M[3] * M[8] - M[5] * M[6]) + M[2] * (M[3] * M[7] - M[4] * M[6]);
  /** Rotation by angle t (radians, counter-clockwise) about the unit vector a. */
  function rotation(a, t) {
    const c = Math.cos(t), s = Math.sin(t), C = 1 - c, x = a[0], y = a[1], z = a[2];
    return [c + x * x * C, x * y * C - z * s, x * z * C + y * s,
      y * x * C + z * s, c + y * y * C, y * z * C - x * s,
      z * x * C - y * s, z * y * C + x * s, c + z * z * C];
  }
  /** Reflection through the plane with unit normal v. */
  function reflection(v) {
    const x = v[0], y = v[1], z = v[2];
    return [1 - 2 * x * x, -2 * x * y, -2 * x * z, -2 * x * y, 1 - 2 * y * y, -2 * y * z, -2 * x * z, -2 * y * z, 1 - 2 * z * z];
  }
  const improper = (a, t) => mmul(reflection(a), rotation(a, t));
  function maxDiff(A, B) { let m = 0; for (let i = 0; i < 9; i++) { const d = Math.abs(A[i] - B[i]); if (d > m) m = d; } return m; }

  /** Eigen-decomposition of a symmetric 3x3 matrix by Jacobi rotations; values ascending. */
  function eigSym(A0) {
    const a = A0.slice(), V = IDENT.slice();
    for (let sweep = 0; sweep < 60; sweep++) {
      const off = a[1] * a[1] + a[2] * a[2] + a[5] * a[5];
      const dg = a[0] * a[0] + a[4] * a[4] + a[8] * a[8];
      if (off <= 1e-30 * dg || off < 1e-300) break;
      for (let pq = 0; pq < 3; pq++) {
        const p = pq === 2 ? 1 : 0, q = pq === 0 ? 1 : 2;
        const apq = a[3 * p + q];
        if (apq === 0) continue;
        const th = (a[3 * q + q] - a[3 * p + p]) / (2 * apq);
        const t = (th >= 0 ? 1 : -1) / (Math.abs(th) + Math.sqrt(th * th + 1));
        const c = 1 / Math.sqrt(t * t + 1), s = t * c;
        for (let k = 0; k < 3; k++) {
          const akp = a[3 * k + p], akq = a[3 * k + q];
          a[3 * k + p] = c * akp - s * akq; a[3 * k + q] = s * akp + c * akq;
        }
        for (let k = 0; k < 3; k++) {
          const apk = a[3 * p + k], aqk = a[3 * q + k];
          a[3 * p + k] = c * apk - s * aqk; a[3 * q + k] = s * apk + c * aqk;
        }
        for (let k = 0; k < 3; k++) {
          const vkp = V[3 * k + p], vkq = V[3 * k + q];
          V[3 * k + p] = c * vkp - s * vkq; V[3 * k + q] = s * vkp + c * vkq;
        }
      }
    }
    const idx = [0, 1, 2].sort((i, j) => a[4 * i] - a[4 * j]);
    return { values: idx.map(i => a[4 * i]), vectors: idx.map(i => [V[i], V[3 + i], V[6 + i]]) };
  }

  // ------------------------------------------------------------------ preparation
  function massOf(z) {
    const m = z > 0 && BM.el && BM.el.mass ? BM.el.mass[z] : 0;
    return m > 0 ? m : 0;
  }

  /** Centre of mass, centred coordinates, inertia tensor. */
  function prepare(Z, xyz) {
    if (!Z || !xyz || xyz.length < 3 * Z.length) throw new Error('Symmetry needs Z and xyz arrays of matching size.');
    const n = Z.length;
    const type = new Int32Array(n);
    const w = new Float64Array(n);
    let W = 0;
    for (let i = 0; i < n; i++) {
      const z = Number(Z[i]);
      type[i] = z > 0 ? Math.round(z) : 0;
      w[i] = massOf(type[i]); W += w[i];
      if (!Number.isFinite(xyz[3 * i]) || !Number.isFinite(xyz[3 * i + 1]) || !Number.isFinite(xyz[3 * i + 2])) {
        throw new Error('Atom ' + (i + 1) + ' has a non-numeric coordinate.');
      }
    }
    if (!(W > 0)) { w.fill(1); W = n; }
    const c = [0, 0, 0];
    for (let i = 0; i < n; i++) for (let k = 0; k < 3; k++) c[k] += w[i] * xyz[3 * i + k];
    if (n) for (let k = 0; k < 3; k++) c[k] /= W;
    const r = new Float64Array(3 * n), rad = new Float64Array(n);
    const I = [0, 0, 0, 0, 0, 0, 0, 0, 0];
    let wr = 0;
    for (let i = 0; i < n; i++) {
      const x = xyz[3 * i] - c[0], y = xyz[3 * i + 1] - c[1], z = xyz[3 * i + 2] - c[2];
      r[3 * i] = x; r[3 * i + 1] = y; r[3 * i + 2] = z;
      rad[i] = Math.sqrt(x * x + y * y + z * z);
      const m = w[i];
      wr += m * rad[i];
      I[0] += m * (y * y + z * z); I[4] += m * (x * x + z * z); I[8] += m * (x * x + y * y);
      I[1] -= m * x * y; I[2] -= m * x * z; I[5] -= m * y * z;
    }
    I[3] = I[1]; I[6] = I[2]; I[7] = I[5];
    // Typical (median) distance from an atom to the nearest atom of the same element. A tolerance near half
    // of it cannot tell which atom an operation maps onto, so detect() clamps tol below that. The median
    // ignores the odd duplicated atom.
    const nn = new Float64Array(n).fill(Infinity);
    for (let i = 0; i < n; i++) {
      for (let j = i + 1; j < n; j++) {
        if (type[i] !== type[j]) continue;
        const d = Math.hypot(r[3 * i] - r[3 * j], r[3 * i + 1] - r[3 * j + 1], r[3 * i + 2] - r[3 * j + 2]);
        if (d < nn[i]) nn[i] = d;
        if (d < nn[j]) nn[j] = d;
      }
    }
    const fin = Array.from(nn).filter(Number.isFinite).sort((a, b) => a - b);
    const dSame = fin.length ? fin[fin.length >> 1] : Infinity;
    return { n, type, w, center: c, r, rad, inertia: eigSym(I), wr, dSame };
  }

  /** Sets of possibly equivalent atoms: same element and the same sorted distances (within 2 tol)
   *  to the atoms of every element. Any operation that passes the tol test maps an atom inside its set. */
  function equivalenceSets(P, tol) {
    const { n, r, rad, type } = P;
    const D = new Float64Array(n * n);
    for (let i = 0; i < n; i++) {
      for (let j = i + 1; j < n; j++) {
        const dx = r[3 * i] - r[3 * j], dy = r[3 * i + 1] - r[3 * j + 1], dz = r[3 * i + 2] - r[3 * j + 2];
        D[i * n + j] = D[j * n + i] = Math.sqrt(dx * dx + dy * dy + dz * dz);
      }
    }
    const types = Array.from(new Set(type)).sort((a, b) => a - b);
    const sig = new Array(n);
    for (let i = 0; i < n; i++) {
      const parts = [];
      for (const t of types) {
        const g = [];
        for (let j = 0; j < n; j++) if (j !== i && type[j] === t) g.push(D[i * n + j]);
        g.sort((a, b) => a - b);
        for (const d of g) parts.push(d);
      }
      sig[i] = Float64Array.from(parts);
    }
    const parent = new Int32Array(n);
    for (let i = 0; i < n; i++) parent[i] = i;
    const find = i => { while (parent[i] !== i) { parent[i] = parent[parent[i]]; i = parent[i]; } return i; };
    const thr = 2 * tol + 1e-9;
    for (let i = 0; i < n; i++) {
      for (let j = i + 1; j < n; j++) {
        if (type[i] !== type[j] || Math.abs(rad[i] - rad[j]) > thr) continue;
        const fi = find(i), fj = find(j);
        if (fi === fj) continue;
        const a = sig[i], b = sig[j];
        let same = true;
        for (let k = 0; k < a.length; k++) if (Math.abs(a[k] - b[k]) > thr) { same = false; break; }
        if (same) parent[fi] = fj;
      }
    }
    return groupsFrom(n, find);
  }

  function groupsFrom(n, find) {
    const classOf = new Int32Array(n), members = [], map = new Map();
    for (let i = 0; i < n; i++) {
      const root = find(i);
      let k = map.get(root);
      if (k === undefined) { k = members.length; map.set(root, k); members.push([]); }
      classOf[i] = k; members[k].push(i);
    }
    return { classOf, members };
  }

  function elementSets(P) {
    const map = new Map();
    for (let i = 0; i < P.n; i++) { if (!map.has(P.type[i])) map.set(P.type[i], i); }
    return groupsFrom(P.n, i => map.get(P.type[i]));
  }

  /** Returns test(M, base, ang) -> {perm, err} | null. Every atom i must have an unused atom j of its
   *  set with |M r_i - r_j| <= base + 2 ang |r_i| (ang: angular uncertainty of an approximate element).
   *  The widened tolerance is capped at half the distance the operation moves the atom, so a small
   *  rotation cannot pass by leaving every atom on itself. */
  function makeMatcher(P, sets) {
    const { n, r, rad } = P;
    const order = Array.from({ length: n }, (_, i) => i).sort((a, b) => rad[b] - rad[a]);
    const mark = new Int32Array(n);
    const perm = new Int32Array(n);
    const members = sets.members, classOf = sets.classOf;
    let stamp = 0;
    const test = function (M, base, ang) {
      stamp++;
      let worst = 0;
      for (let o = 0; o < n; o++) {
        const i = order[o];
        const x = r[3 * i], y = r[3 * i + 1], z = r[3 * i + 2];
        const px = M[0] * x + M[1] * y + M[2] * z, py = M[3] * x + M[4] * y + M[5] * z, pz = M[6] * x + M[7] * y + M[8] * z;
        let t = base + 2 * ang * rad[i];
        if (t > base) {
          const ex = px - x, ey = py - y, ez = pz - z;
          const half = 0.5 * Math.sqrt(ex * ex + ey * ey + ez * ez);
          const cap = half > base ? half : base;
          if (t > cap) t = cap;
        }
        let bd = t * t, best = -1;
        const S = members[classOf[i]];
        test.work += S.length;
        for (let s = 0; s < S.length; s++) {
          const j = S[s];
          if (mark[j] === stamp) continue;
          const dx = r[3 * j] - px, dy = r[3 * j + 1] - py, dz = r[3 * j + 2] - pz;
          const d2 = dx * dx + dy * dy + dz * dz;
          if (d2 <= bd) { bd = d2; best = j; }
        }
        if (best < 0) return null;
        mark[best] = stamp; perm[i] = best;
        if (bd > worst) worst = bd;
      }
      return { perm: Int32Array.from(perm), err: Math.sqrt(worst) };
    };
    test.work = 0;                  // atom-distance evaluations so far (bounds the search on huge inputs)
    return test;
  }

  // ------------------------------------------------------------------ element refinement (least squares)
  function topVector(M, guess) {
    M[3] = M[1]; M[6] = M[2]; M[7] = M[5];
    const e = eigSym(M);
    const span = Math.abs(e.values[2]) + Math.abs(e.values[0]) + 1e-300;
    if (e.values[2] - e.values[1] < 1e-9 * span) return guess;
    const v = e.vectors[2];
    return dot(v, guess) < 0 ? vscale(v, -1) : v;
  }
  /** Axis of a rotation from its permutation: cycle centroids lie on the axis and the displacements
   *  r_p(i) - r_i are perpendicular to it. */
  function refineAxis(P, perm, guess) {
    const { n, r } = P;
    const M = [0, 0, 0, 0, 0, 0, 0, 0, 0];
    const seen = new Uint8Array(n);
    for (let i = 0; i < n; i++) {
      if (seen[i]) continue;
      const cyc = [];
      for (let j = i; !seen[j]; j = perm[j]) { seen[j] = 1; cyc.push(j); }
      let cx = 0, cy = 0, cz = 0;
      for (const k of cyc) { cx += r[3 * k]; cy += r[3 * k + 1]; cz += r[3 * k + 2]; }
      const L = cyc.length;
      cx /= L; cy /= L; cz /= L;
      M[0] += L * cx * cx; M[1] += L * cx * cy; M[2] += L * cx * cz; M[4] += L * cy * cy; M[5] += L * cy * cz; M[8] += L * cz * cz;
    }
    for (let i = 0; i < n; i++) {
      const j = perm[i];
      const dx = r[3 * j] - r[3 * i], dy = r[3 * j + 1] - r[3 * i + 1], dz = r[3 * j + 2] - r[3 * i + 2];
      M[0] -= dx * dx; M[1] -= dx * dy; M[2] -= dx * dz; M[4] -= dy * dy; M[5] -= dy * dz; M[8] -= dz * dz;
    }
    return topVector(M, guess);
  }
  /** Mirror normal from its permutation: displacements are along the normal, midpoints in the plane. */
  function refinePlane(P, perm, guess) {
    const { n, r } = P;
    const M = [0, 0, 0, 0, 0, 0, 0, 0, 0];
    for (let i = 0; i < n; i++) {
      const j = perm[i];
      const dx = r[3 * j] - r[3 * i], dy = r[3 * j + 1] - r[3 * i + 1], dz = r[3 * j + 2] - r[3 * i + 2];
      const mx = (r[3 * j] + r[3 * i]) / 2, my = (r[3 * j + 1] + r[3 * i + 1]) / 2, mz = (r[3 * j + 2] + r[3 * i + 2]) / 2;
      M[0] += dx * dx - mx * mx; M[1] += dx * dy - mx * my; M[2] += dx * dz - mx * mz;
      M[4] += dy * dy - my * my; M[5] += dy * dz - my * mz; M[8] += dz * dz - mz * mz;
    }
    return topVector(M, guess);
  }

  // ------------------------------------------------------------------ candidate directions
  /** Candidate directions closer than mergeAngle are merged, keeping the one with the larger lever. The
   *  angle must stay below what the kept candidate's loose test absorbs (tol / rMax), or an exact
   *  candidate could be replaced by a nearby inexact one (e.g. an arbitrary degenerate inertia axis). */
  function CandSet(minLever, mergeAngle) {
    this.minLever = minLever; this.mergeCos = Math.cos(mergeAngle); this.list = []; this.grid = new Map();
  }
  CandSet.prototype.cellKey = (a, b, c) => (a + 64) * 16384 + (b + 64) * 128 + (c + 64);
  CandSet.prototype.near = function (u) {
    const bx = Math.floor(u[0] / CELL), by = Math.floor(u[1] / CELL), bz = Math.floor(u[2] / CELL);
    for (let dx = -1; dx <= 1; dx++) {
      for (let dy = -1; dy <= 1; dy++) {
        for (let dz = -1; dz <= 1; dz++) {
          const cell = this.grid.get(this.cellKey(bx + dx, by + dy, bz + dz));
          if (cell) for (const c of cell) if (dot(c.u, u) > this.mergeCos) return c;
        }
      }
    }
    return null;
  };
  /** lever: length scale of the construction; the direction error is about tol / lever. */
  CandSet.prototype.add = function (v, lever) {
    if (!(lever >= this.minLever)) return;
    const l = vnorm(v);
    if (!(l > 1e-12)) return;
    const u = [v[0] / l, v[1] / l, v[2] / l], neg = [-u[0], -u[1], -u[2]];
    const c = this.near(u) || this.near(neg);
    if (c) {
      if (lever > c.lever) { c.lever = lever; c.u = dot(c.u, u) >= 0 ? u : neg; }
      return;
    }
    const e = { u, lever };
    this.list.push(e);
    const key = this.cellKey(Math.floor(u[0] / CELL), Math.floor(u[1] / CELL), Math.floor(u[2] / CELL));
    let cell = this.grid.get(key);
    if (!cell) this.grid.set(key, cell = []);
    cell.push(e);
  };
  CandSet.prototype.sorted = function () { return this.list.slice().sort((a, b) => b.lever - a.lever); };

  /** Finds proper axes (with their highest order), mirror planes, inversion and S2n axes. */
  function findElements(P, sets, test, tol) {
    const { n, r, rad } = P;
    const pos = i => [r[3 * i], r[3 * i + 1], r[3 * i + 2]];
    let nMax = 1, rMax = 0;
    for (const S of sets.members) if (S.length > nMax) nMax = S.length;
    for (let i = 0; i < n; i++) if (rad[i] > rMax) rMax = rad[i];
    const mergeAngle = Math.max(1e-7, Math.min(0.02, 0.25 * tol / Math.max(rMax, 1e-9)));
    const axisC = new CandSet(tol, mergeAngle), planeC = new CandSet(tol, mergeAngle);
    for (const v of P.inertia.vectors) { axisC.add(v, rMax); planeC.add(v, rMax); }
    for (let i = 0; i < n; i++) if (rad[i] > tol) axisC.add(pos(i), rad[i]);
    const bySize = sets.members.filter(S => S.length > 1).sort((a, b) => a.length - b.length);
    let work = 0, good = false;
    for (const S of bySize) {
      if (good && work > WORK_LIMIT) break;
      const m = S.length, i0 = S[0];
      if (rad[i0] <= tol) continue;
      const Q = S.map(pos), p0 = Q[0];
      for (let a = 0; a < m; a++) {
        for (let b = a + 1; b < m; b++) {
          const mid = vscale(vadd(Q[a], Q[b]), 0.5);
          axisC.add(mid, vnorm(mid));                       // C2 through the midpoint of an equivalent pair
          const d = vsub(Q[a], Q[b]);
          planeC.add(d, vnorm(d) / 2);                      // mirror bisecting an equivalent pair
        }
      }
      for (let b = 1; b < m; b++) {
        const c = cross(p0, Q[b]), lv = vnorm(c) / (rad[i0] + rad[S[b]]);
        axisC.add(c, lv); planeC.add(c, lv);                // normal of the plane through the centre and two atoms
      }
      for (let b = 1; b < m; b++) {
        const e1 = vsub(Q[b], p0), l1 = vnorm(e1);
        for (let c = b + 1; c < m; c++) {
          const e2 = vsub(Q[c], p0), t = cross(e1, e2);
          axisC.add(t, vnorm(t) / (2 * (l1 + vnorm(e2))));  // Cn (n >= 3) through a ring of equivalent atoms
        }
      }
      work += 2 * m * m;
      if (m >= 3 && !good) {
        const u = unit(p0);
        for (let b = 1; b < m; b++) if (vnorm(cross(Q[b], u)) > tol) { good = true; break; }
      }
    }

    /** Possible orders k of a rotation about u, largest first: the atom farthest from the axis must have an
     *  equivalent atom at the same height and radius, rotated by 2 pi / k (within base + 2 ang |r|). */
    const ordersFor = (u, base, ang) => {
      let p = -1, pr = 0, q0 = null;
      for (let i = 0; i < n; i++) {
        const v = pos(i), q = vsub(v, vscale(u, dot(v, u))), l = vnorm(q);
        if (l > pr) { pr = l; p = i; q0 = q; }
      }
      if (p < 0) return [];
      const tp = base + 2 * ang * rad[p];
      if (pr <= tp) {                      // every atom is close to this axis: no filtering possible
        const all = [];
        for (let k = nMax; k >= 2; k--) all.push(k);
        return all;
      }
      const e1 = vscale(q0, 1 / pr), e2 = cross(u, e1), h = dot(pos(p), u);
      const dphi = Math.min(Math.PI, 2.4 * Math.asin(Math.min(1, tp / (2 * pr))) + 1e-9);
      const ks = new Set();
      for (const j of sets.members[sets.classOf[p]]) {
        if (j === p) continue;
        const v = pos(j), a = dot(v, u);
        if (Math.abs(a - h) > tp) continue;
        const q = vsub(v, vscale(u, a));
        if (Math.abs(vnorm(q) - pr) > tp) continue;
        let phi = Math.atan2(dot(q, e2), dot(q, e1));
        if (phi <= 0) phi += TWO_PI;
        const kHi = phi - dphi > 0 ? Math.min(nMax, Math.floor(TWO_PI / (phi - dphi))) : nMax;
        for (let k = Math.max(2, Math.ceil(TWO_PI / (phi + dphi))); k <= kHi; k++) ks.add(k);
      }
      return Array.from(ks).sort((a, b) => b - a);
    };
    const axes = [];
    const maxOrder = ax => {
      for (const k of ordersFor(ax, tol, 0)) {
        const res = test(rotation(ax, TWO_PI / k), tol, 0);
        if (res) return { n: k, perm: res.perm };
      }
      return null;
    };
    const covered = u => axes.some(a => Math.abs(dot(a.axis, u)) > COS_SKIP);
    /** Approximate axis u -> loose test for each order, refine, strict test; true if a new axis was added. */
    const tryAxis = (u, ang) => {
      if (covered(u)) return false;
      for (const k of ordersFor(u, 2 * tol, ang)) {
        const loose = test(rotation(u, TWO_PI / k), 2 * tol, ang);
        if (!loose) continue;
        let ax = refineAxis(P, loose.perm, u);
        let best = maxOrder(ax);
        if (!best) continue;
        const ax2 = refineAxis(P, best.perm, ax), best2 = maxOrder(ax2);
        if (best2 && best2.n >= best.n) { ax = ax2; best = best2; }
        if (covered(ax)) return false;
        axes.push({ axis: ax, n: best.n, perm: best.perm, s2n: false });
        return true;
      }
      return false;
    };
    const work0 = test.work, spent = () => test.work - work0 > TEST_BUDGET;
    if (nMax >= 2) {
      for (const c of axisC.sorted()) { if (spent()) break; tryAxis(c.u, Math.min(0.2, tol / c.lever)); }
      // Products of found rotations are symmetry operations too, and their axes are accurate:
      // this recovers axes whose own candidates were too inexact (e.g. C2 x C2' -> C3 or C5 in I).
      for (let round = 0; round < 3 && !spent(); round++) {
        let added = false;
        const list = axes.slice();
        for (let a = 0; a < list.length && !spent(); a++) {
          const Ra = rotation(list[a].axis, TWO_PI / list[a].n);
          for (let b = a + 1; b < list.length; b++) {
            for (const sgn of list[b].n > 2 ? [1, -1] : [1]) {
              const ra = rotAxisAngle(mmul(Ra, rotation(list[b].axis, sgn * TWO_PI / list[b].n)));
              if (ra.axis && tryAxis(ra.axis, 0.02)) added = true;
            }
          }
        }
        if (!added) break;
      }
    }
    const planes = [];
    const work1 = test.work;
    for (const c of planeC.sorted()) {
      if (test.work - work1 > TEST_BUDGET / 2) break;
      if (planes.some(p => Math.abs(dot(p.normal, c.u)) > COS_SKIP)) continue;
      const ang = Math.min(0.2, tol / c.lever);
      const loose = test(reflection(c.u), 2 * tol, ang);
      if (!loose) continue;
      let nv = refinePlane(P, loose.perm, c.u);
      let st = test(reflection(nv), tol, 0);
      if (!st) continue;
      const nv2 = refinePlane(P, st.perm, nv), st2 = test(reflection(nv2), tol, 0);
      if (st2) { nv = nv2; st = st2; }
      planes.push({ normal: nv, perm: st.perm });
    }
    const inversion = !!test(INVERSION, tol, 0);
    for (const a of axes) a.s2n = !!test(improper(a.axis, Math.PI / a.n), tol, 0);
    return { axes, planes, inversion };
  }

  // ------------------------------------------------------------------ flowchart
  const CUBIC_ORDER = { T: 12, Td: 24, Th: 24, O: 24, Oh: 48, I: 60, Ih: 120 };

  /** Point group from symmetry elements {axes: [{axis, n, s2n}], planes: [{normal}], inversion}. */
  function classify(el) {
    const { axes, planes, inversion } = el;
    if (!axes.length) {
      if (planes.length) return { name: 'Cs', order: 2, kind: 'Cs', sigma: planes[0] };
      if (inversion) return { name: 'Ci', order: 2, kind: 'Ci' };
      return { name: 'C1', order: 1, kind: 'C1' };
    }
    if (axes.filter(a => a.n >= 3).length >= 2) {
      const family = axes.some(a => a.n % 5 === 0) ? 'I' : axes.some(a => a.n % 4 === 0) ? 'O' : 'T';
      let name;
      if (family === 'T') name = inversion ? 'Th' : planes.length ? 'Td' : 'T';
      else name = inversion || planes.length ? family + 'h' : family;
      return { name, order: CUBIC_ORDER[name], kind: 'cubic', family };
    }
    let nmax = 0;
    for (const a of axes) if (a.n > nmax) nmax = a.n;
    let best = null;
    for (const main of axes) {
      if (main.n !== nmax) continue;
      const c = classifyAxial(main, el);
      if (!best || c.order > best.order) best = c;
    }
    return best;
  }
  function classifyAxial(main, el) {
    const n = main.n, z = main.axis;
    const perp = el.axes.filter(a => a !== main && a.n % 2 === 0 && Math.abs(dot(a.axis, z)) < SIN_T);
    const sh = el.planes.find(p => Math.abs(dot(p.normal, z)) > COS_T);
    const sv = el.planes.filter(p => Math.abs(dot(p.normal, z)) < SIN_T);
    if (perp.length) {
      if (sh) return { name: 'D' + n + 'h', order: 4 * n, kind: 'Dnh', n, main, sec: perp[0] };
      if (sv.length) return { name: 'D' + n + 'd', order: 4 * n, kind: 'Dnd', n, main, sec: perp[0] };
      return { name: 'D' + n, order: 2 * n, kind: 'Dn', n, main, sec: perp[0] };
    }
    if (sh) return { name: 'C' + n + 'h', order: 2 * n, kind: 'Cnh', n, main };
    if (sv.length) return { name: 'C' + n + 'v', order: 2 * n, kind: 'Cnv', n, main, sigma: sv[0] };
    if (main.s2n) return { name: 'S' + 2 * n, order: 2 * n, kind: 'S2n', n, main };
    return { name: 'C' + n, order: n, kind: 'Cn', n, main };
  }

  // ------------------------------------------------------------------ exact groups from generators
  function orthoFrame(z, xh) {
    z = unit(z);
    let x = vsub(xh, vscale(z, dot(xh, z)));
    x = vnorm(x) > 1e-9 ? unit(x) : perpendicular(z);
    return { x, y: cross(z, x), z };
  }
  function frameFrom2(a, b) {
    const x = unit(a);
    const y = unit(vsub(b, vscale(x, dot(b, x))));
    return { x, y, z: cross(x, y) };
  }
  function toWorld(fr, T) {
    const F = [fr.x[0], fr.y[0], fr.z[0], fr.x[1], fr.y[1], fr.z[1], fr.x[2], fr.y[2], fr.z[2]];
    return mmul(mmul(F, T), mT(F));
  }
  const fromFrame = (fr, v) => [fr.x[0] * v[0] + fr.y[0] * v[1] + fr.z[0] * v[2],
    fr.x[1] * v[0] + fr.y[1] * v[1] + fr.z[1] * v[2], fr.x[2] * v[0] + fr.y[2] * v[1] + fr.z[2] * v[2]];
  const Rz = t => rotation([0, 0, 1], t);
  const C2Z = diag(-1, -1, 1), C2X = diag(1, -1, -1), SIGMA_H = diag(1, 1, -1), SIGMA_Y = diag(1, -1, 1);
  const C3D = [0, 0, 1, 1, 0, 0, 0, 1, 0];           // +120 degrees about (1,1,1): x -> y -> z -> x
  const SWAP_XY = [0, 1, 0, 1, 0, 0, 0, 0, 1];       // mirror x = y (sigma_d of Td)

  /** Closes the exact generator matrices together with their atom permutations. Each new operation's
   *  permutation is the composition of its factors', so the permutations form a group homomorphism by
   *  construction (which makes symmetrize exact). Returns null if two products of the same operation
   *  disagree on the permutation, if more than `limit` operations appear, or if any operation leaves an
   *  atom farther than 2 tol from its image atom. */
  function closeGroup(P, gens, testEl, tol, limit) {
    const n = P.n, r = P.r;
    const gp = [];
    for (const g of gens) {
      const res = testEl(g, 2 * tol, 0);
      if (!res) return null;
      gp.push(res.perm);
    }
    const id = new Int32Array(n);
    for (let i = 0; i < n; i++) id[i] = i;
    const mats = [IDENT.slice()], perms = [id];
    for (let i = 0; i < mats.length; i++) {
      for (let t = 0; t < gens.length; t++) {
        const M = mmul(gens[t], mats[i]), pi = perms[i], s = gp[t];
        const q = new Int32Array(n);
        for (let a = 0; a < n; a++) q[a] = s[pi[a]];
        let k = -1;
        for (let m = 0; m < mats.length; m++) if (maxDiff(M, mats[m]) < 1e-6) { k = m; break; }
        if (k >= 0) {
          const pk = perms[k];
          for (let a = 0; a < n; a++) if (pk[a] !== q[a]) return null;
          continue;
        }
        mats.push(M); perms.push(q);
        if (mats.length > limit) return null;
      }
    }
    let dev = 0;
    const lim = 2 * tol;
    for (let m = 0; m < mats.length; m++) {
      const M = mats[m], q = perms[m];
      for (let a = 0; a < n; a++) {
        const x = r[3 * a], y = r[3 * a + 1], z = r[3 * a + 2], b = q[a];
        const d = Math.hypot(M[0] * x + M[1] * y + M[2] * z - r[3 * b], M[3] * x + M[4] * y + M[5] * z - r[3 * b + 1],
          M[6] * x + M[7] * y + M[8] * z - r[3 * b + 2]);
        if (d > lim) return null;
        if (d > dev) dev = d;
      }
    }
    return { mats, perms, dev };
  }

  function inertiaFrame(P) {
    const v = P.inertia.vectors;
    return { x: v[0].slice(), y: v[1].slice(), z: cross(v[0], v[1]) };
  }
  /** Direction perpendicular to z towards the atom farthest from the axis (a tidy frame x axis). */
  function perpHint(P, z) {
    let best = null, bl = 1e-6;
    for (let i = 0; i < P.n; i++) {
      const p = [P.r[3 * i], P.r[3 * i + 1], P.r[3 * i + 2]];
      const q = vsub(p, vscale(z, dot(p, z))), l = vnorm(q);
      if (l > bl + 1e-6) { bl = l; best = q; }
    }
    return best || perpendicular(z);
  }

  /** Symmetry frame, exact generator matrices and principal axis for a classified group. */
  function groupSetup(P, cl, el, tol, testEl) {
    let fr, gens, pa = null;
    if (cl.kind === 'C1' || cl.kind === 'Ci') {
      fr = inertiaFrame(P);
      gens = cl.kind === 'Ci' ? [INVERSION] : [];
    } else if (cl.kind === 'Cs') {
      const z = unit(cl.sigma.normal);
      fr = orthoFrame(z, perpHint(P, z));
      gens = [SIGMA_H]; pa = fr.z;
    } else if (cl.kind === 'cubic') {
      const fam = cl.family;
      const ax = el.axes.filter(a => (fam === 'O' ? a.n % 4 === 0 : a.n === 2));
      fr = null;
      for (let a = 0; a < ax.length && !fr; a++) {
        for (let b = a + 1; b < ax.length; b++) {
          if (Math.abs(dot(ax[a].axis, ax[b].axis)) < SIN_T) { fr = frameFrom2(ax[a].axis, ax[b].axis); break; }
        }
      }
      if (!fr && fam !== 'I') {
        const c3 = el.axes.filter(a => a.n % 3 === 0);
        if (c3.length >= 2) {
          const a = c3[0].axis;
          let b = c3[1].axis;
          if (dot(a, b) > 0) b = vscale(b, -1);
          const x = unit(vadd(a, b));
          fr = frameFrom2(x, mvec(rotation(a, TWO_PI / 3), x));
        }
      }
      if (!fr) return null;
      if (fam === 'T') {
        gens = [C2Z, C2X, C3D];
        if (cl.name === 'Td') gens.push(SWAP_XY);
        if (cl.name === 'Th') gens.push(INVERSION);
        pa = unit(fromFrame(fr, [1, 1, 1]));
      } else if (fam === 'O') {
        gens = [Rz(Math.PI / 2), C3D];
        if (cl.name === 'Oh') gens.push(INVERSION);
        pa = fr.z;
      } else {
        // Two icosahedra share a set of three perpendicular C2 axes; keep the one that fits.
        let c5 = null;
        for (const v of [[0, 1, PHI], [1, 0, PHI]]) {
          const R = rotation(unit(v), TWO_PI / 5);
          if (testEl(toWorld(fr, R), 2 * tol, 0)) { c5 = { v, R }; break; }
        }
        if (!c5) return null;
        gens = [C2Z, C2X, C3D, c5.R];
        if (cl.name === 'Ih') gens.push(INVERSION);
        pa = unit(fromFrame(fr, c5.v));
      }
    } else {
      const n = cl.n, z = unit(cl.main.axis);
      const xh = cl.sec ? cl.sec.axis : cl.sigma ? cross(cl.sigma.normal, z) : perpHint(P, z);
      fr = orthoFrame(z, xh);
      pa = fr.z;
      const Cn = Rz(TWO_PI / n);
      switch (cl.kind) {
        case 'Cn': gens = [Cn]; break;
        case 'Cnv': gens = [Cn, SIGMA_Y]; break;
        case 'Cnh': gens = [Cn, SIGMA_H]; break;
        case 'S2n': gens = [mmul(SIGMA_H, Rz(Math.PI / n))]; break;
        case 'Dn': gens = [Cn, C2X]; break;
        case 'Dnh': gens = [Cn, C2X, SIGMA_H]; break;
        case 'Dnd': gens = [Cn, C2X, reflection([-Math.sin(Math.PI / (2 * n)), Math.cos(Math.PI / (2 * n)), 0])]; break;
        default: return null;
      }
    }
    return { frame: fr, gens: gens.map(T => toWorld(fr, T)), principalAxis: pa ? canon(pa) : null };
  }

  // ------------------------------------------------------------------ describing operations
  function rotAxisAngle(R) {
    const c = Math.max(-1, Math.min(1, (R[0] + R[4] + R[8] - 1) / 2));
    const theta = Math.acos(c);
    if (theta < 1e-7) return { axis: null, theta: 0 };
    const w = [R[7] - R[5], R[2] - R[6], R[3] - R[1]];          // 2 sin(theta) * axis
    if (theta < Math.PI - 1e-5 && vnorm(w) > 1e-9) return { axis: unit(w), theta };
    const B = [R[0] + 1, R[1], R[2], R[3], R[4] + 1, R[5], R[6], R[7], R[8] + 1];  // = 2 a a^T at theta = pi
    let j = 0;
    if (B[4] > B[4 * j]) j = 1;
    if (B[8] > B[4 * j]) j = 2;
    return { axis: canon(unit([B[j], B[3 + j], B[6 + j]])), theta: Math.PI };
  }
  function angleFraction(theta) {
    const f = theta / TWO_PI;
    for (let n = 1; n <= 360; n++) {
      const k = Math.round(f * n);
      if (k > 0 && Math.abs(f * n - k) < 1e-5) return { n, k };
    }
    return { n: Math.max(1, Math.round(1 / f)), k: 1 };
  }
  function describe(M) {
    if (det3(M) > 0) {
      const { axis, theta } = rotAxisAngle(M);
      if (!axis) return { type: 'E' };
      const { n, k } = angleFraction(theta);
      return { type: 'C', n, k, axis };
    }
    const { axis, theta } = rotAxisAngle(M.map(v => -v));        // -M = C(theta) about axis
    if (!axis) return { type: 'i' };
    if (theta >= Math.PI - 1e-5) return { type: 'sigma', normal: canon(axis) };
    const { n, k } = angleFraction(Math.PI - theta);              // M = S(pi - theta) about -axis
    return { type: 'S', n, k, axis: vscale(axis, -1) };
  }
  const TYPE_RANK = { E: 0, C: 1, i: 2, S: 3, sigma: 4 };
  const SIGMA_RANK = { 'σh': 0, 'σv': 1, 'σd': 2, 'σ': 3 };

  function sigmaLabel(normal, group, fr, kind, n) {
    if (!fr || group === 'Cs' || group === 'Ih') return 'σ';
    if (group === 'Th') return 'σh';
    if (group === 'Oh') return [fr.x, fr.y, fr.z].some(a => Math.abs(dot(a, normal)) > COS_T) ? 'σh' : 'σd';
    if (group === 'Td') return 'σd';
    const z = fr.z, c = Math.abs(dot(normal, z));
    if (c > COS_T) return 'σh';
    if (c > SIN_T) return 'σ';
    if (kind === 'Dnd') return 'σd';
    if ((kind === 'Dnh' || kind === 'Cnv') && n % 2 === 0 && n > 2) return evenClass(cross(normal, z), fr, n) ? 'σv' : 'σd';
    return 'σv';
  }
  /** For D2n/C2nv-type groups: is the in-plane direction d in the class of the frame x axis? */
  function evenClass(d, fr, n) {
    let a = Math.atan2(dot(d, fr.y), dot(d, fr.x));
    if (a < 0) a += Math.PI;
    return Math.round(a / (Math.PI / n)) % 2 === 0;
  }
  function rotationLabel(d, fr, kind, n) {
    let s = d.type + d.n + (d.k > 1 ? '^' + d.k : '');
    if (d.type === 'C' && d.n === 2 && fr && (kind === 'Dn' || kind === 'Dnh' || kind === 'Dnd') &&
        Math.abs(dot(d.axis, fr.z)) < SIN_T) {
      s += kind !== 'Dnd' && n % 2 === 0 && !evenClass(d.axis, fr, n) ? '″' : '′';
    }
    return s;
  }

  function makeOperations(mats, perms, group, fr, kind, n) {
    const ops = mats.map((M, idx) => {
      const d = describe(M);
      const op = { type: d.type };
      if (d.type === 'C' || d.type === 'S') { op.n = d.n; op.k = d.k; op.axis = d.axis; }
      if (d.type === 'sigma') op.normal = d.normal;
      op.label = d.type === 'sigma' ? sigmaLabel(d.normal, group, fr, kind, n)
        : d.type === 'C' || d.type === 'S' ? rotationLabel(d, fr, kind, n) : d.type;
      op.matrix = M.map(v => (Math.abs(v) < 1e-15 ? 0 : v));
      op.perm = Array.from(perms[idx]);
      return op;
    });
    ops.sort((a, b) => (TYPE_RANK[a.type] - TYPE_RANK[b.type]) ||
      (a.type === 'sigma' ? SIGMA_RANK[a.label] - SIGMA_RANK[b.label] : 0) ||
      ((b.n || 0) - (a.n || 0)) || ((a.k || 0) - (b.k || 0)));
    return ops;
  }

  function orbits(n, perms) {
    const parent = new Int32Array(n);
    for (let i = 0; i < n; i++) parent[i] = i;
    const find = i => { while (parent[i] !== i) { parent[i] = parent[parent[i]]; i = parent[i]; } return i; };
    for (const p of perms) for (let i = 0; i < n; i++) { const a = find(i), b = find(p[i]); if (a !== b) parent[a] = b; }
    return groupsFrom(n, find).members;
  }

  function topTypeOf(P, tol, group) {
    if (P.n === 1) return 'atom';
    if (group === 'Cinfv' || group === 'Dinfh') return 'linear';
    if (CUBIC_ORDER[group]) return 'spherical';
    const m = /^[CDS](\d+)/.exec(group);
    if (m && +m[1] >= 3 && !(group[0] === 'S' && +m[1] === 4)) return 'symmetric';
    if (group === 'S4' || group === 'D2d') return 'symmetric';
    const [a, b, c] = P.inertia.values;
    const eps = Math.max(1e-9 * c, tol * P.wr);
    const ab = Math.abs(b - a) < eps, bc = Math.abs(c - b) < eps;
    return ab && bc ? 'spherical' : ab || bc ? 'symmetric' : 'asymmetric';
  }

  // ------------------------------------------------------------------ detection
  function result(P, tol, group, order, mats, perms, frame, pa, kind, n, maxDev) {
    return {
      pointGroup: group, order,
      operations: makeOperations(mats, perms, group, frame, kind, n),
      center: P.center.slice(),
      principalAxis: pa,
      frame: frame ? { x: frame.x.slice(), y: frame.y.slice(), z: frame.z.slice() } : null,
      topType: topTypeOf(P, tol, group),
      moments: P.inertia.values.slice(),
      classes: orbits(P.n, perms),
      tol, maxDeviation: maxDev
    };
  }

  function detectLinear(P, tol, testEl) {
    const u = P.inertia.vectors[0];
    for (let i = 0; i < P.n; i++) {
      const p = [P.r[3 * i], P.r[3 * i + 1], P.r[3 * i + 2]];
      if (vnorm(vsub(p, vscale(u, dot(p, u)))) > tol) return null;   // no |p|^2 - (p.u)^2 cancellation
    }
    const inv = !!testEl(INVERSION, tol, 0);
    const fr = orthoFrame(u, perpendicular(u));
    const gens = (inv ? [C2Z, C2X, SIGMA_H] : [C2Z, SIGMA_Y]).map(T => toWorld(fr, T));
    const g = closeGroup(P, gens, testEl, tol, 8);
    if (!g) return null;
    const group = inv ? 'Dinfh' : 'Cinfv';
    return result(P, tol, group, Infinity, g.mats, g.perms, fr, canon(u), 'linear', 0, g.dev);
  }

  function detectAt(P, tol) {
    const testEl = makeMatcher(P, elementSets(P));
    const lin = detectLinear(P, tol, testEl);
    if (lin) return lin;
    const sets = equivalenceSets(P, tol);
    const el = findElements(P, sets, makeMatcher(P, sets), tol);
    const cl = classify(el);
    const setup = groupSetup(P, cl, el, tol, testEl);
    if (!setup) return null;
    const g = closeGroup(P, setup.gens, testEl, tol, cl.order);
    if (!g || g.mats.length !== cl.order) return null;
    return result(P, tol, cl.name, cl.order, g.mats, g.perms, setup.frame, setup.principalAxis, cl.kind, cl.n || 0, g.dev);
  }

  function detect(Z, xyz, opts) {
    opts = opts || {};
    const tol0 = opts.tol > 0 ? +opts.tol : 0.1;
    const P = prepare(Z, xyz);
    if (P.n === 0) {
      return { pointGroup: 'C1', order: 1, operations: makeOperations([IDENT.slice()], [[]], 'C1', null, 'C1', 0),
        center: [0, 0, 0], principalAxis: null, frame: null, topType: 'atom', moments: [0, 0, 0], classes: [],
        tol: tol0, maxDeviation: 0 };
    }
    if (P.n === 1) {
      const fr = { x: [1, 0, 0], y: [0, 1, 0], z: [0, 0, 1] };
      return result(P, tol0, 'Kh', Infinity, [IDENT.slice(), INVERSION.slice()], [[0], [0]], fr, null, 'atom', 0, 0);
    }
    let tol = Number.isFinite(P.dSame) ? Math.min(tol0, TOL_CAP * P.dSame) : tol0;
    for (let attempt = 0; attempt < 5; attempt++) {
      const res = detectAt(P, tol);
      if (res) return res;
      tol /= 2;
    }
    const ident = new Int32Array(P.n);
    for (let i = 0; i < P.n; i++) ident[i] = i;
    return result(P, tol, 'C1', 1, [IDENT.slice()], [ident], inertiaFrame(P), null, 'C1', 0, 0);
  }

  // ------------------------------------------------------------------ names, orders, subgroups
  function normalizeName(s) {
    let t = String(s === undefined || s === null ? '' : s).trim().replace(/[\s_{}$\\]/g, '');
    if (!t) return null;
    const lin = /^([cd])(?:inf|infinity|∞|oo|\*)([vh])$/i.exec(t);
    if (lin) {
      const g = lin[1].toUpperCase() + 'inf' + lin[2].toLowerCase();
      return g === 'Cinfv' || g === 'Dinfh' ? g : null;
    }
    t = t[0].toUpperCase() + t.slice(1).toLowerCase();
    const alias = { C1h: 'Cs', C1v: 'Cs', S1: 'Cs', S2: 'Ci', D1: 'C2', D1h: 'C2v', D1d: 'C2h', C3i: 'S6', C5i: 'S10' };
    if (alias[t]) return alias[t];
    let m = /^S(\d+)$/.exec(t);
    if (m && +m[1] % 2 === 1 && +m[1] >= 3) return 'C' + m[1] + 'h';
    if (/^(C1|Cs|Ci|T|Td|Th|O|Oh|I|Ih|Kh)$/.test(t)) return t;
    m = /^([CDS])(\d+)([vhd]?)$/.exec(t);
    if (!m) return null;
    const k = +m[2];
    if (m[1] === 'C' && k >= 2 && m[3] !== 'd') return t;
    if (m[1] === 'D' && k >= 2) return t;
    if (m[1] === 'S' && k >= 4 && k % 2 === 0 && !m[3]) return t;
    return null;
  }

  function groupOrder(name) {
    const g = normalizeName(name);
    if (!g) return NaN;
    if (g === 'C1') return 1;
    if (g === 'Cs' || g === 'Ci') return 2;
    if (CUBIC_ORDER[g]) return CUBIC_ORDER[g];
    if (g === 'Cinfv' || g === 'Dinfh' || g === 'Kh') return Infinity;
    const m = /^([CDS])(\d+)([vhd]?)$/.exec(g);
    const n = +m[2];
    if (m[1] === 'S') return n;
    if (m[1] === 'C') return m[3] ? 2 * n : n;
    return m[3] ? 4 * n : 2 * n;
  }

  /** Point-group name of an exact set of described operations. */
  function nameOfOps(ops) {
    const axes = [], planes = [];
    let inversion = false;
    for (const o of ops) {
      if (o.type === 'C') {
        const a = axes.find(x => Math.abs(dot(x.axis, o.axis)) > 1 - 1e-6);
        if (a) a.n = Math.max(a.n, o.n); else axes.push({ axis: o.axis, n: o.n, s2n: false });
      } else if (o.type === 'sigma') planes.push({ normal: o.normal });
      else if (o.type === 'i') inversion = true;
    }
    for (const o of ops) {
      if (o.type !== 'S') continue;
      const a = axes.find(x => Math.abs(dot(x.axis, o.axis)) > 1 - 1e-6);
      if (a && o.n === 2 * a.n) a.s2n = true;
    }
    return classify({ axes, planes, inversion }).name;
  }

  /** Operations of a subgroup called `want` inside the detected group, or null. */
  function subgroupOps(res, want) {
    const ops = res.operations, N = ops.length;
    const E = ops.find(o => o.type === 'E');
    if (want === 'C1') return [E];
    if (res.pointGroup === 'Kh') return want === 'Ci' ? ops.slice() : null;
    if (res.pointGroup === 'Dinfh' && want === 'Cinfv') {
      const z = res.principalAxis;
      return ops.filter(o => o.type === 'E' || (o.type === 'C' && Math.abs(dot(o.axis, z)) > COS_T) ||
        (o.type === 'sigma' && Math.abs(dot(o.normal, z)) < SIN_T));
    }
    const target = groupOrder(want);
    if (!(target >= 1) || !Number.isFinite(target) || target > N || N % target !== 0) return null;
    const key = M => M.map(v => Math.round(v * 1e5)).join(',');
    const index = new Map();
    ops.forEach((o, i) => index.set(key(o.matrix), i));
    const lookup = M => {
      const k = index.get(key(M));
      if (k !== undefined) return k;
      for (let i = 0; i < N; i++) if (maxDiff(M, ops[i].matrix) < 1e-6) return i;
      return -1;
    };
    const tab = new Int32Array(N * N);
    for (let a = 0; a < N; a++) for (let b = 0; b < N; b++) tab[a * N + b] = lookup(mmul(ops[a].matrix, ops[b].matrix));
    const e = ops.indexOf(E);
    const gen = seeds => {
      const inH = new Uint8Array(N), H = [e];
      inH[e] = 1;
      for (let i = 0; i < H.length; i++) {
        for (const g of seeds) {
          const p = tab[g * N + H[i]];
          if (p < 0) return null;
          if (!inH[p]) { inH[p] = 1; H.push(p); if (H.length > target) return null; }
        }
      }
      return H;
    };
    const fits = H => H && H.length === target && nameOfOps(H.map(i => ops[i])) === want;
    const pick = H => H.slice().sort((a, b) => a - b).map(i => ops[i]);
    for (let g = 0; g < N; g++) { const H = gen([g]); if (fits(H)) return pick(H); }
    const seen = new Set(), halves = [];
    for (let a = 0; a < N; a++) {
      for (let b = a + 1; b < N; b++) {
        const H = gen([a, b]);
        if (!H) continue;
        const k = H.slice().sort((x, y) => x - y).join(',');
        if (seen.has(k)) continue;
        seen.add(k);
        if (fits(H)) return pick(H);
        if (2 * H.length === target) halves.push(H);
      }
    }
    for (const H of halves) {
      const inH = new Set(H);
      for (let g = 0; g < N; g++) {
        if (inH.has(g)) continue;
        const H2 = gen(H.concat([g]));
        if (fits(H2)) return pick(H2);
      }
    }
    return null;
  }

  // ------------------------------------------------------------------ symmetrization
  function symmetrize(Z, xyz, opts) {
    opts = opts || {};
    const tol = opts.tol > 0 ? +opts.tol : 0.1;
    const n = Z.length;
    if (!xyz || xyz.length < 3 * n) throw new Error('Symmetry needs Z and xyz arrays of matching size.');
    let res = detect(Z, xyz, { tol });
    let ops = res.operations, group = res.pointGroup;
    if (opts.pointGroup !== undefined && opts.pointGroup !== null && opts.pointGroup !== '') {
      const want = normalizeName(opts.pointGroup);
      if (!want) throw new Error('Unknown point group "' + opts.pointGroup + '".');
      if (want !== group) {
        let sub = subgroupOps(res, want);
        const dSame = prepare(Z, xyz).dSame;
        const cap = Number.isFinite(dSame) ? TOL_CAP * dSame : Infinity;
        let tried = Math.min(tol, cap);
        for (const f of [1.5, 2, 3, 4, 6, 8]) {
          if (sub) break;
          const t = Math.min(tol * f, cap);
          if (t <= tried * (1 + 1e-9)) break;           // clamped: a larger tolerance would not change anything
          tried = t;
          const r2 = detect(Z, xyz, { tol: t });
          const s2 = r2.pointGroup === want ? r2.operations : subgroupOps(r2, want);
          if (s2) { res = r2; sub = s2; }
        }
        if (!sub) {
          throw new Error('The structure is not close to ' + want + ' symmetry (detected ' + group +
            '; tolerances up to ' + +tried.toFixed(3) + ' Å were tried).');
        }
        ops = sub; group = want;
      }
    }
    const c = res.center, G = ops.length;
    const r = new Float64Array(3 * n);
    for (let i = 0; i < n; i++) for (let k = 0; k < 3; k++) r[3 * i + k] = xyz[3 * i + k] - c[k];
    let s = average(r, ops.map(o => o.perm), ops, n);
    // The permutations come from independent nearest-atom matches; if they were not consistent the
    // average would not be exactly symmetric, so re-match on the averaged structure and repeat.
    for (let iter = 0; iter < 3 && !isExact(s, ops, ops.map(o => o.perm), n); iter++) {
      const perms = rematch(s, Z, ops, n);
      if (!perms) break;
      s = average(r, perms, ops, n);
      ops = ops.map((o, k) => Object.assign({}, o, { perm: Array.from(perms[k]) }));
    }
    const out = new Float64Array(3 * n);
    let maxShift = 0, sum = 0;
    for (let i = 0; i < n; i++) {
      let d2 = 0;
      for (let k = 0; k < 3; k++) {
        out[3 * i + k] = s[3 * i + k] + c[k];
        const d = out[3 * i + k] - xyz[3 * i + k];
        d2 += d * d;
      }
      sum += d2;
      if (d2 > maxShift) maxShift = d2;
    }
    const order = group === res.pointGroup ? res.order : group === 'Cinfv' || group === 'Dinfh' ? Infinity : G;
    return { xyz: out, pointGroup: group, order, maxShift: Math.sqrt(maxShift), rmsShift: n ? Math.sqrt(sum / n) : 0 };
  }

  function average(r, perms, ops, n) {
    const s = new Float64Array(3 * n), G = ops.length;
    for (let i = 0; i < n; i++) {
      let ax = 0, ay = 0, az = 0;
      for (let g = 0; g < G; g++) {
        const M = ops[g].matrix, j = perms[g][i];
        const x = r[3 * j], y = r[3 * j + 1], z = r[3 * j + 2];
        ax += M[0] * x + M[3] * y + M[6] * z;                    // M^T r_j  (M^-1 = M^T)
        ay += M[1] * x + M[4] * y + M[7] * z;
        az += M[2] * x + M[5] * y + M[8] * z;
      }
      s[3 * i] = ax / G; s[3 * i + 1] = ay / G; s[3 * i + 2] = az / G;
    }
    return s;
  }
  function isExact(s, ops, perms, n) {
    let scale = 1;
    for (let i = 0; i < 3 * n; i++) scale = Math.max(scale, Math.abs(s[i]));
    for (let g = 0; g < ops.length; g++) {
      const M = ops[g].matrix, p = perms[g];
      for (let i = 0; i < n; i++) {
        const q = mvec(M, [s[3 * i], s[3 * i + 1], s[3 * i + 2]]), j = p[i];
        if (Math.abs(q[0] - s[3 * j]) + Math.abs(q[1] - s[3 * j + 1]) + Math.abs(q[2] - s[3 * j + 2]) > 1e-10 * scale) return false;
      }
    }
    return true;
  }
  function rematch(s, Z, ops, n) {
    const out = [];
    for (const o of ops) {
      const M = o.matrix, used = new Uint8Array(n), p = new Int32Array(n);
      for (let i = 0; i < n; i++) {
        const q = mvec(M, [s[3 * i], s[3 * i + 1], s[3 * i + 2]]);
        let best = -1, bd = Infinity;
        for (let j = 0; j < n; j++) {
          if (used[j] || (Z[j] > 0 ? Z[j] : 0) !== (Z[i] > 0 ? Z[i] : 0)) continue;
          const d = (q[0] - s[3 * j]) ** 2 + (q[1] - s[3 * j + 1]) ** 2 + (q[2] - s[3 * j + 2]) ** 2;
          if (d < bd) { bd = d; best = j; }
        }
        if (best < 0) return null;
        used[best] = 1; p[i] = best;
      }
      out.push(p);
    }
    return out;
  }

  BM.sym = { detect, symmetrize, normalizeName, groupOrder };
})(typeof window !== 'undefined' ? (window.BM = window.BM || {}) : (globalThis.BM = globalThis.BM || {}));
