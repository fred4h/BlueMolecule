/* BlueMolecule — geometry utilities: bonds, measurements, editing, formula. */
(function (BM) {
  'use strict';
  const el = BM.el;
  const DEG = 180 / Math.PI;

  /** Distance-based bonds. Returns Uint32Array of pairs [i0, j0, i1, j1, ...]. */
  function computeBonds(Z, xyz, opts) {
    opts = opts || {};
    const tol = opts.tolerance !== undefined ? opts.tolerance : 0.4;
    const n = Z.length;
    let maxR = 0;
    for (let i = 0; i < n; i++) if (Z[i] > 0 && el.covalent[Z[i]] > maxR) maxR = el.covalent[Z[i]];
    if (n < 2 || maxR === 0) return new Uint32Array(0);
    const cell = 2 * maxR + tol;
    let minx = Infinity, miny = Infinity, minz = Infinity;
    for (let i = 0; i < n; i++) {
      if (xyz[3 * i] < minx) minx = xyz[3 * i];
      if (xyz[3 * i + 1] < miny) miny = xyz[3 * i + 1];
      if (xyz[3 * i + 2] < minz) minz = xyz[3 * i + 2];
    }
    const grid = new Map();
    const keyOf = (ix, iy, iz) => (ix * 1048576 + iy) * 1048576 + iz; // safe below 2^53 for reasonable grids
    const cx = new Int32Array(n), cy = new Int32Array(n), cz = new Int32Array(n);
    for (let i = 0; i < n; i++) {
      if (Z[i] <= 0) continue;
      const ix = Math.floor((xyz[3 * i] - minx) / cell), iy = Math.floor((xyz[3 * i + 1] - miny) / cell), iz = Math.floor((xyz[3 * i + 2] - minz) / cell);
      cx[i] = ix; cy[i] = iy; cz[i] = iz;
      const k = keyOf(ix, iy, iz);
      let list = grid.get(k);
      if (!list) { list = []; grid.set(k, list); }
      list.push(i);
    }
    const out = [];
    for (let i = 0; i < n; i++) {
      if (Z[i] <= 0) continue;
      const ri = el.covalent[Z[i]];
      const xi = xyz[3 * i], yi = xyz[3 * i + 1], zi = xyz[3 * i + 2];
      for (let dx = -1; dx <= 1; dx++) for (let dy = -1; dy <= 1; dy++) for (let dz = -1; dz <= 1; dz++) {
        const list = grid.get(keyOf(cx[i] + dx, cy[i] + dy, cz[i] + dz));
        if (!list) continue;
        for (let t = 0; t < list.length; t++) {
          const j = list[t];
          if (j <= i) continue;
          const ddx = xyz[3 * j] - xi, ddy = xyz[3 * j + 1] - yi, ddz = xyz[3 * j + 2] - zi;
          const d2 = ddx * ddx + ddy * ddy + ddz * ddz;
          const lim = ri + el.covalent[Z[j]] + tol;
          if (d2 < lim * lim && d2 > 0.16) { out.push(i, j); }
        }
      }
    }
    return Uint32Array.from(out);
  }

  function adjacency(n, bonds) {
    const adj = Array.from({ length: n }, () => []);
    for (let b = 0; b < bonds.length; b += 2) { adj[bonds[b]].push(bonds[b + 1]); adj[bonds[b + 1]].push(bonds[b]); }
    return adj;
  }

  function dist(xyz, i, j) {
    const dx = xyz[3 * j] - xyz[3 * i], dy = xyz[3 * j + 1] - xyz[3 * i + 1], dz = xyz[3 * j + 2] - xyz[3 * i + 2];
    return Math.sqrt(dx * dx + dy * dy + dz * dz);
  }
  function vec(xyz, i, j) { return [xyz[3 * j] - xyz[3 * i], xyz[3 * j + 1] - xyz[3 * i + 1], xyz[3 * j + 2] - xyz[3 * i + 2]]; }
  const dot = (a, b) => a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
  const cross = (a, b) => [a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]];
  const norm = a => Math.sqrt(dot(a, a));
  const scale = (a, s) => [a[0] * s, a[1] * s, a[2] * s];
  const sub = (a, b) => [a[0] - b[0], a[1] - b[1], a[2] - b[2]];
  function unit(a) { const l = norm(a); return l > 1e-12 ? scale(a, 1 / l) : [0, 0, 0]; }

  function angle(xyz, i, j, k) {
    const a = vec(xyz, j, i), b = vec(xyz, j, k);
    const c = dot(a, b) / (norm(a) * norm(b));
    return Math.acos(Math.max(-1, Math.min(1, c))) * DEG;
  }
  /** IUPAC-signed dihedral i-j-k-l in degrees (-180, 180]. */
  function dihedral(xyz, i, j, k, l) {
    const b0 = vec(xyz, j, i), b1 = unit(vec(xyz, j, k)), b2 = vec(xyz, k, l);
    const v = sub(b0, scale(b1, dot(b0, b1))), w = sub(b2, scale(b1, dot(b2, b1)));
    const x = dot(v, w), y = dot(cross(b1, v), w);
    return Math.atan2(y, x) * DEG;
  }

  /** Atoms reachable from `start` without crossing the direct edge start->`cut`. Returns {set, ring}. */
  function side(adj, start, cut) {
    const seen = new Uint8Array(adj.length);
    const stack = [start]; seen[start] = 1;
    let ring = false;
    const list = [];
    while (stack.length) {
      const u = stack.pop(); list.push(u);
      for (const v of adj[u]) {
        if (u === start && v === cut) continue;
        if (v === cut) { ring = true; continue; }
        if (!seen[v]) { seen[v] = 1; stack.push(v); }
      }
    }
    return { atoms: list, ring };
  }

  function rotateAtoms(xyz, atoms, origin, axis, theta) {
    const a = unit(axis), c = Math.cos(theta), s = Math.sin(theta);
    for (const p of atoms) {
      const v = [xyz[3 * p] - origin[0], xyz[3 * p + 1] - origin[1], xyz[3 * p + 2] - origin[2]];
      const axv = cross(a, v), ad = dot(a, v) * (1 - c);
      xyz[3 * p] = origin[0] + v[0] * c + axv[0] * s + a[0] * ad;
      xyz[3 * p + 1] = origin[1] + v[1] * c + axv[1] * s + a[1] * ad;
      xyz[3 * p + 2] = origin[2] + v[2] * c + axv[2] * s + a[2] * ad;
    }
  }
  const pos = (xyz, i) => [xyz[3 * i], xyz[3 * i + 1], xyz[3 * i + 2]];

  /** Choose which side of a bond to move: the smaller one, falling back to a single atom inside rings. */
  function chooseSide(adj, a, b) {
    const sb = side(adj, b, a);
    const sa = side(adj, a, b);
    if (sb.ring || sa.ring) return { atoms: [b], moveB: true, ring: true };
    return sb.atoms.length <= sa.atoms.length ? { atoms: sb.atoms, moveB: true, ring: false } : { atoms: sa.atoms, moveB: false, ring: false };
  }

  function setDistance(xyz0, bonds, i, j, target) {
    if (!(target > 0.05)) throw new Error('Distance must be a positive number of Ångströms.');
    const xyz = Float64Array.from(xyz0);
    const adj = adjacency(xyz.length / 3, bonds);
    const pick = chooseSide(adj, i, j);
    const d = dist(xyz, i, j);
    if (d < 1e-8) throw new Error('The two atoms coincide.');
    const u = unit(vec(xyz, i, j));
    const shift = (target - d) * (pick.moveB ? 1 : -1);
    for (const p of pick.atoms) { xyz[3 * p] += u[0] * shift; xyz[3 * p + 1] += u[1] * shift; xyz[3 * p + 2] += u[2] * shift; }
    return { xyz, moved: pick.atoms.length, ring: pick.ring };
  }

  function setAngle(xyz0, bonds, i, j, k, targetDeg) {
    const xyz = Float64Array.from(xyz0);
    const adj = adjacency(xyz.length / 3, bonds);
    const cur = angle(xyz, i, j, k);
    let n = cross(vec(xyz, j, i), vec(xyz, j, k));
    if (norm(n) < 1e-8) { const a = vec(xyz, j, i); n = cross(a, Math.abs(a[0]) < 0.9 ? [1, 0, 0] : [0, 1, 0]); }
    const delta = (targetDeg - cur) / DEG;
    const sk = side(adj, k, j), si = side(adj, i, j);
    let atoms, theta, ring = false;
    if (sk.ring || si.ring) { atoms = [k]; theta = delta; ring = true; }
    else if (sk.atoms.length <= si.atoms.length) { atoms = sk.atoms; theta = delta; }
    else { atoms = si.atoms; theta = -delta; }
    rotateAtoms(xyz, atoms, pos(xyz, j), n, theta);
    return { xyz, moved: atoms.length, ring };
  }

  function setDihedral(xyz0, bonds, i, j, k, l, targetDeg) {
    const xyz = Float64Array.from(xyz0);
    const adj = adjacency(xyz.length / 3, bonds);
    const sk = side(adj, k, j), sj = side(adj, j, k);
    if (sk.ring || sj.ring) throw new Error('The central bond is part of a ring, so the dihedral cannot be changed by a rigid rotation.');
    let delta = targetDeg - dihedral(xyz, i, j, k, l);
    while (delta > 180) delta -= 360;
    while (delta <= -180) delta += 360;
    const axis = vec(xyz, j, k);
    if (sk.atoms.length <= sj.atoms.length) rotateAtoms(xyz, sk.atoms, pos(xyz, k), axis, delta / DEG);
    else rotateAtoms(xyz, sj.atoms, pos(xyz, j), axis, -delta / DEG);
    return { xyz, moved: Math.min(sk.atoms.length, sj.atoms.length), ring: false };
  }

  function formula(Z) {
    const counts = new Map();
    for (const z of Z) if (z > 0) counts.set(z, (counts.get(z) || 0) + 1);
    const syms = [...counts.keys()].map(z => [el.symbol(z), counts.get(z)]);
    const hasC = counts.has(6);
    syms.sort((a, b) => {
      if (hasC) {
        const rank = s => (s === 'C' ? 0 : s === 'H' ? 1 : 2);
        if (rank(a[0]) !== rank(b[0])) return rank(a[0]) - rank(b[0]);
      }
      return a[0] < b[0] ? -1 : a[0] > b[0] ? 1 : 0;
    });
    return syms.map(([s, c]) => s + (c > 1 ? c : '')).join('');
  }
  function molWeight(Z) { let m = 0; for (const z of Z) if (z > 0) m += el.mass[z]; return m; }

  function bounds(xyz, Z, onlyVisible) {
    let cx = 0, cy = 0, cz = 0, cnt = 0;
    const n = xyz.length / 3;
    for (let i = 0; i < n; i++) { if (onlyVisible && !onlyVisible(i)) continue; cx += xyz[3 * i]; cy += xyz[3 * i + 1]; cz += xyz[3 * i + 2]; cnt++; }
    if (!cnt) return { center: [0, 0, 0], radius: 1 };
    cx /= cnt; cy /= cnt; cz /= cnt;
    let r = 0;
    for (let i = 0; i < n; i++) {
      if (onlyVisible && !onlyVisible(i)) continue;
      const dx = xyz[3 * i] - cx, dy = xyz[3 * i + 1] - cy, dz = xyz[3 * i + 2] - cz;
      const rr = Math.sqrt(dx * dx + dy * dy + dz * dz) + (Z ? (el.vdw[Z[i]] || 1.5) : 1.5);
      if (rr > r) r = rr;
    }
    return { center: [cx, cy, cz], radius: Math.max(r, 1.5) };
  }

  /** Displacement vector scaled so the largest atomic displacement equals `amp` Å. */
  function normalizedMode(mode, amp) {
    let m = 0;
    for (let i = 0; i < mode.length; i += 3) {
      const d = Math.hypot(mode[i], mode[i + 1], mode[i + 2]);
      if (d > m) m = d;
    }
    const s = m > 0 ? amp / m : 0;
    const out = new Float64Array(mode.length);
    for (let i = 0; i < mode.length; i++) out[i] = mode[i] * s;
    return out;
  }

  function atomLabel(Z, i) { return el.symbol(Z[i]) + (i + 1); }

  BM.geom = { computeBonds, adjacency, dist, angle, dihedral, side, setDistance, setAngle, setDihedral, formula, molWeight, bounds, normalizedMode, atomLabel, rotateAtoms };
})(typeof window !== 'undefined' ? (window.BM = window.BM || {}) : (globalThis.BM = globalThis.BM || {}));
