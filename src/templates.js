/* BlueMolecule — building library (BM.templates).
 *
 * Coordination geometries, ligands, substituent groups and small molecules with explicit 3D coordinates.
 * Coordinates are generated here in code from ideal internal coordinates — the shared bond-length rule
 * (BM.ff.idealBondLength), VSEPR / coordination polyhedra for directions, rigid ring units — and then
 * polished with BM.ff.optimize.  Entries are built lazily on first use and cached.
 *
 * Editing operations (insertMolecule, attachGroup, buildCenter, replaceLigand, addBidentate) never mutate
 * their input; they return {mol, added, removed?, map?}.  New atoms are appended at the end and relaxed with
 * a short BM.ff.optimize in which every other atom is frozen.
 */
(function (BM) {
  'use strict';
  const el = BM.el, ff = BM.ff;
  const r0 = ff.idealBondLength;
  const RAD = Math.PI / 180;

  // ------------------------------------------------------------------------------------------- vectors
  const add = (a, b) => [a[0] + b[0], a[1] + b[1], a[2] + b[2]];
  const sub = (a, b) => [a[0] - b[0], a[1] - b[1], a[2] - b[2]];
  const mul = (a, s) => [a[0] * s, a[1] * s, a[2] * s];
  const dot = (a, b) => a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
  const cross = (a, b) => [a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]];
  const len = a => Math.sqrt(dot(a, a));
  const unit = a => { const l = len(a); return l > 1e-12 ? mul(a, 1 / l) : [0, 0, 0]; };
  const perp = a => unit(Math.abs(a[0]) < 0.9 * len(a) ? cross(a, [1, 0, 0]) : cross(a, [0, 1, 0]));
  const perpPart = (a, n) => sub(a, mul(n, dot(a, n)));
  const I3 = [1, 0, 0, 0, 1, 0, 0, 0, 1];
  function rotAxis(axis, th) {
    const [x, y, z] = unit(axis), c = Math.cos(th), s = Math.sin(th), t = 1 - c;
    return [t * x * x + c, t * x * y - s * z, t * x * z + s * y,
      t * x * y + s * z, t * y * y + c, t * y * z - s * x,
      t * x * z - s * y, t * y * z + s * x, t * z * z + c];
  }
  const apply = (R, v) => [R[0] * v[0] + R[1] * v[1] + R[2] * v[2], R[3] * v[0] + R[4] * v[1] + R[5] * v[2], R[6] * v[0] + R[7] * v[1] + R[8] * v[2]];
  function matMul(A, B) {
    const C = new Array(9);
    for (let r = 0; r < 3; r++) for (let c = 0; c < 3; c++) C[3 * r + c] = A[3 * r] * B[c] + A[3 * r + 1] * B[3 + c] + A[3 * r + 2] * B[6 + c];
    return C;
  }
  /** Rotation taking direction a onto direction b. */
  function rotFromTo(a, b) {
    a = unit(a); b = unit(b);
    const c = dot(a, b);
    if (c > 1 - 1e-12) return I3.slice();
    if (c < -1 + 1e-12) return rotAxis(perp(a), Math.PI);
    return rotAxis(cross(a, b), Math.acos(c));
  }
  /** Rotation whose columns map the orthonormal frame (a1, a2, a1×a2) onto (b1, b2, b1×b2). */
  function frameRotation(a1, a2, b1, b2) {
    const a3 = cross(a1, a2), b3 = cross(b1, b2), R = new Array(9);
    for (let r = 0; r < 3; r++) for (let c = 0; c < 3; c++) R[3 * r + c] = b1[r] * a1[c] + b2[r] * a2[c] + b3[r] * a3[c];
    return R;
  }
  const twistAngle = (n, a, b) => { const ap = perpPart(a, n), bp = perpPart(b, n); return Math.atan2(dot(cross(ap, bp), n), dot(ap, bp)); };

  // ------------------------------------------------------------------------------------------- site directions
  const GEO = ff.GEOMETRIES;
  /** Ideal site set for an atom from BM.ff.perceive information (null if it has no defined polyhedron). */
  function siteSet(info) {
    if (info.geom && GEO[info.geom]) return GEO[info.geom].vectors;
    if (!info.metal) {
      if (info.sn === 2) return GEO.linear.vectors;
      if (info.sn === 3) return GEO['trigonal-planar'].vectors;
      if (info.sn === 4) return GEO.tetrahedral.vectors;
    }
    return null;
  }
  /** Free sites of an ideal polyhedron after matching the existing neighbour directions to it.  With a single
   * neighbour the twist about that bond puts the first free site anti to `ref` (staggered / zig-zag). */
  function openDirections(sites, dirs, ref) {
    if (!dirs.length) return sites.map(v => v.slice());
    if (dirs.length >= sites.length) return [];
    const d0 = unit(dirs[0]);
    if (dirs.length === 1) {
      const R = rotFromTo(sites[0], d0);
      let free = sites.slice(1).map(v => apply(R, v));
      const want = ref ? mul(ref, -1) : null;
      if (want && len(perpPart(want, d0)) > 1e-6 && len(perpPart(free[0], d0)) > 1e-6) {
        const R2 = rotAxis(d0, twistAngle(d0, free[0], want));
        free = free.map(v => apply(R2, v));
      }
      return free;
    }
    const asg = ff.assignSites(dirs.map(unit), sites);
    const R1 = rotFromTo(sites[asg.perm[0]], d0);
    let sy = 0, sx = 0;
    for (let k = 1; k < dirs.length; k++) {
      const ap = perpPart(apply(R1, sites[asg.perm[k]]), d0), bp = perpPart(unit(dirs[k]), d0);
      sy += dot(cross(ap, bp), d0); sx += dot(ap, bp);
    }
    const R = matMul(rotAxis(d0, Math.atan2(sy, sx)), R1);
    const used = new Set(asg.perm);
    return sites.filter((_, p) => !used.has(p)).map(v => apply(R, v));
  }
  /** Fallback: `count` directions spread as far as possible from the existing ones. */
  function spreadDirections(dirs, count) {
    const out = [], taken = dirs.map(unit), N = 240, cand = [];
    for (let k = 0; k < N; k++) {
      const z = 1 - 2 * (k + 0.5) / N, r = Math.sqrt(1 - z * z), ph = k * 2.399963229728653;
      cand.push([r * Math.cos(ph), r * Math.sin(ph), z]);
    }
    for (let c = 0; c < count; c++) {
      if (!taken.length) { out.push([1, 0, 0]); taken.push([1, 0, 0]); continue; }
      let best = cand[0], bestS = -Infinity;
      for (const v of cand) {
        let m = Infinity;
        for (const t of taken) m = Math.min(m, 1 - dot(v, t));
        if (m > bestS) { bestS = m; best = v; }
      }
      out.push(best); taken.push(best);
    }
    return out;
  }

  // ------------------------------------------------------------------------------------------- fragment builder
  /** Internal builder: atoms, bonds and rigid units (rings with local coordinates); build() places the rest
   * breadth-first on ideal open sites at ideal bond lengths, then polishes with the force field. */
  class Frag {
    constructor() { this.Z = []; this.pos = []; this.bonds = []; this.fc = []; this.units = []; }
    atom(sym, p) {
      this.Z.push(typeof sym === 'number' ? sym : el.toZ(sym));
      this.pos.push(p ? p.slice() : null); this.fc.push(0);
      return this.Z.length - 1;
    }
    bond(i, j, o) { this.bonds.push([i, j, o || 1]); return this; }
    to(i, sym, o, nH) { const j = this.atom(sym); this.bond(i, j, o || 1); if (nH) this.H(j, nH); return j; }
    H(i, n) { for (let k = 0; k < (n === undefined ? 1 : n); k++) this.to(i, 'H'); return this; }
    /** A ring of new atoms as one rigid unit: regular polygon, 'chair' (6) or 'envelope' (5). */
    ring(syms, order, shape) {
      const n = syms.length, ids = syms.map(s => this.atom(s));
      const orders = Array.isArray(order) ? order : ids.map(() => order || 1);
      let side = 0;
      for (let k = 0; k < n; k++) {
        this.bond(ids[k], ids[(k + 1) % n], orders[k]);
        side += r0(this.Z[ids[k]], this.Z[ids[(k + 1) % n]], orders[k]);
      }
      side /= n;
      let local;
      if (shape === 'chair' && n === 6) {
        const th = 109.47 * RAD, R = 2 * side * Math.sin(th / 2) / Math.sqrt(3), h = Math.sqrt(Math.max(0, side * side - R * R) / 4);
        local = ids.map((_, k) => [R * Math.cos(k * Math.PI / 3), R * Math.sin(k * Math.PI / 3), k % 2 ? h : -h]);
      } else {
        const Rc = side / (2 * Math.sin(Math.PI / n));
        local = ids.map((_, k) => [Rc * Math.cos(2 * Math.PI * k / n), Rc * Math.sin(2 * Math.PI * k / n), 0]);
        if (shape === 'envelope') local[2 % n][2] = 0.55;
      }
      this.units.push({ atoms: ids, local });
      return ids;
    }
    /** Existing atoms as one rigid unit with explicit local coordinates. */
    rigid(ids, coords) { this.units.push({ atoms: ids.slice(), local: coords.map(c => [c[0], c[1], c[2] || 0]) }); }

    build(opts) {
      opts = opts || {};
      const n = this.Z.length, Z = this.Z, pos = this.pos.map(p => (p ? p.slice() : null));
      const nbr = Array.from({ length: n }, () => []), ordOf = new Map();
      for (const [i, j, o] of this.bonds) { nbr[i].push(j); nbr[j].push(i); ordOf.set(i * n + j, o); ordOf.set(j * n + i, o); }
      const info = ff.perceive({ Z, xyz: null, bonds: this.bonds, formalCharges: this.fc }, { geometries: opts.geometries }).atoms;
      const unitOf = new Int32Array(n).fill(-1);
      this.units.forEach((u, k) => u.atoms.forEach(a => { unitOf[a] = k; }));
      const placed = pos.map(p => !!p), queue = [];
      for (let i = 0; i < n; i++) if (placed[i]) queue.push(i);
      const placeUnitLocal = k => this.units[k].atoms.forEach((a, t) => { pos[a] = this.units[k].local[t].slice(); placed[a] = true; queue.push(a); });
      for (let root = 0; root < n; root++) {
        if (!queue.length && !placed[root]) {
          if (unitOf[root] >= 0) placeUnitLocal(unitOf[root]);
          else { pos[root] = [0, 0, 0]; placed[root] = true; queue.push(root); }
          if (root > 0) { const off = [8 * root, 0, 0]; for (const a of queue) pos[a] = add(pos[a], off); }
        }
        while (queue.length) {
          const a = queue.shift();
          const un = nbr[a].filter(j => !placed[j]);
          if (!un.length) continue;
          un.sort((p, q) => (Z[q] > 1 ? 1 : 0) - (Z[p] > 1 ? 1 : 0));
          const have = nbr[a].filter(j => placed[j]);
          const dirs = have.map(j => unit(sub(pos[j], pos[a])));
          let ref = null;
          if (have.length === 1) {
            const p = have[0];
            const r = nbr[p].find(j => j !== a && placed[j] && Z[j] > 1) ?? nbr[p].find(j => j !== a && placed[j]);
            if (r !== undefined) ref = sub(pos[r], pos[p]);
          }
          const sites = siteSet(info[a]);
          let free = sites && sites.length >= nbr[a].length ? openDirections(sites, dirs, ref) : [];
          if (free.length < un.length) free = free.concat(spreadDirections(dirs.concat(free), un.length - free.length));
          let s = 0;
          for (const u of un) {
            if (placed[u]) continue;
            const dir = free[s++], L = r0(Z[a], Z[u], ordOf.get(a * n + u));
            if (unitOf[u] >= 0) this._placeUnit(unitOf[u], u, a, dir, L, pos, placed, queue, nbr);
            else { pos[u] = add(pos[a], mul(dir, L)); placed[u] = true; queue.push(u); }
          }
        }
      }
      const xyz = new Float64Array(3 * n);
      pos.forEach((p, i) => { xyz[3 * i] = p[0]; xyz[3 * i + 1] = p[1]; xyz[3 * i + 2] = p[2]; });
      const mol = { Z: Z.slice(), xyz, bonds: this.bonds.map(b => b.slice()) };
      if (this.fc.some(q => q)) mol.formalCharges = this.fc.slice();
      if (opts.polish !== false) {
        mol.xyz = ff.optimize(mol, { tol: 0.005, maxSteps: 4000, frozen: opts.frozen, geometries: opts.geometries }).xyz;
      }
      return mol;
    }
    _placeUnit(k, u, a, dir, L, pos, placed, queue, nbr) {
      const U = this.units[k], iu = U.atoms.indexOf(u);
      let ext = [0, 0, 0];
      for (const v of nbr[u]) { const iv = U.atoms.indexOf(v); if (iv >= 0) ext = sub(ext, unit(sub(U.local[iv], U.local[iu]))); }
      ext = len(ext) > 1e-9 ? unit(ext) : [1, 0, 0];
      const base = add(pos[a], mul(dir, L));
      const R1 = rotFromTo(ext, mul(dir, -1));
      const rel = U.local.map(p => apply(R1, sub(p, U.local[iu])));
      const others = [];
      for (let i = 0; i < pos.length; i++) if (placed[i] && i !== a) others.push(pos[i]);
      let bestR = I3, best = Infinity;
      for (let s = 0; s < 36; s++) {
        const R = rotAxis(dir, s * Math.PI / 18);
        let score = 0;
        for (const p of rel) {
          const w = add(base, apply(R, p));
          for (const o of others) { const d = sub(w, o), d2 = dot(d, d); score += 1 / (d2 * d2 * d2 + 1e-6); }
        }
        if (score < best - 1e-12) { best = score; bestR = R; }
      }
      U.atoms.forEach((v, t) => { pos[v] = add(base, apply(bestR, rel[t])); placed[v] = true; queue.push(v); });
    }
  }

  // ------------------------------------------------------------------------------------------- library helpers
  const C6 = ['C', 'C', 'C', 'C', 'C', 'C'];
  const PM = 78;               // a Pt atom stands in for "the metal" while ligand templates are built
  const posOf = (xyz, i) => [xyz[3 * i], xyz[3 * i + 1], xyz[3 * i + 2]];
  const range = n => Array.from({ length: n }, (_, i) => i);

  /** Pointy-top hexagons of side s centred at `centers`: unique vertices and ring edges. */
  function honeycomb(centers, s) {
    const pts = [], edges = new Set();
    const idx = p => {
      for (let i = 0; i < pts.length; i++) if (Math.hypot(pts[i][0] - p[0], pts[i][1] - p[1]) < 1e-6) return i;
      pts.push(p); return pts.length - 1;
    };
    for (const [cx, cy] of centers) {
      const ring = range(6).map(k => idx([cx + s * Math.cos((30 + 60 * k) * RAD), cy + s * Math.sin((30 + 60 * k) * RAD), 0]));
      for (let k = 0; k < 6; k++) { const a = ring[k], b = ring[(k + 1) % 6]; edges.add(a < b ? a + ',' + b : b + ',' + a); }
    }
    return { pts, edges: [...edges].map(e => e.split(',').map(Number)) };
  }
  const nearestPt = (pts, p) => {
    let b = 0, bd = Infinity;
    pts.forEach((q, i) => { const d = Math.hypot(q[0] - p[0], q[1] - p[1]); if (d < bd) { bd = d; b = i; } });
    return b;
  };
  /** Planar fused aromatic rings; `hetero` = [[symbol, [x, y]], ...]; H on every ring-degree-2 carbon. */
  function aromaticPolycycle(f, centers, hetero) {
    const { pts, edges } = honeycomb(centers, r0(6, 6, 1.5));
    const syms = pts.map(() => 'C');
    for (const [sym, p] of hetero) syms[nearestPt(pts, p)] = sym;
    const ids = syms.map(s => f.atom(s)), deg = new Array(pts.length).fill(0);
    for (const [a, b] of edges) { f.bond(ids[a], ids[b], 1.5); deg[a]++; deg[b]++; }
    f.rigid(ids, pts);
    ids.forEach((id, k) => { if (deg[k] === 2 && syms[k] === 'C') f.H(id); });
    return { ids, pts, unit: f.units[f.units.length - 1] };
  }
  /** Position at distance d from both donors, on the side away from `centroid`. */
  function metalBetween(p1, p2, centroid, d) {
    const mid = mul(add(p1, p2), 0.5), half = len(sub(p2, p1)) / 2;
    const w = unit(perpPart(sub(mid, centroid), unit(sub(p2, p1))));
    return add(mid, mul(w, Math.sqrt(Math.max(0.25, d * d - half * half))));
  }

  // ------------------------------------------------------------------------------------------- molecules
  const MOLECULE_SPECS = [
    ['water', 'Water', 'Small molecules', f => { f.H(f.atom('O'), 2); }],
    ['ammonia', 'Ammonia', 'Small molecules', f => { f.H(f.atom('N'), 3); }],
    ['co2', 'Carbon dioxide', 'Small molecules', f => { const c = f.atom('C'); f.to(c, 'O', 2); f.to(c, 'O', 2); }],
    ['formaldehyde', 'Formaldehyde', 'Small molecules', f => { const c = f.atom('C'); f.H(c, 2); f.to(c, 'O', 2); }],
    ['acetic-acid', 'Acetic acid', 'Small molecules', f => { const c = f.atom('C'); f.H(c, 3); const k = f.to(c, 'C'); f.to(k, 'O', 2); f.to(k, 'O', 1, 1); }],
    ['methane', 'Methane', 'Hydrocarbons', f => { f.H(f.atom('C'), 4); }],
    ['ethane', 'Ethane', 'Hydrocarbons', f => { const c = f.atom('C'); f.H(c, 3); f.to(c, 'C', 1, 3); }],
    ['ethylene', 'Ethylene', 'Hydrocarbons', f => { const c = f.atom('C'); f.H(c, 2); f.to(c, 'C', 2, 2); }],
    ['acetylene', 'Acetylene', 'Hydrocarbons', f => { const c = f.atom('C'); f.H(c, 1); f.to(c, 'C', 3, 1); }],
    ['cyclohexane', 'Cyclohexane (chair)', 'Hydrocarbons', f => { f.ring(C6, 1, 'chair').forEach(c => f.H(c, 2)); }],
    ['cyclopentane', 'Cyclopentane', 'Hydrocarbons', f => { f.ring(C6.slice(1), 1, 'envelope').forEach(c => f.H(c, 2)); }],
    ['benzene', 'Benzene', 'Aromatics', f => { f.ring(C6, 1.5).forEach(c => f.H(c)); }],
    ['toluene', 'Toluene', 'Aromatics', f => { const r = f.ring(C6, 1.5); f.to(r[0], 'C', 1, 3); r.slice(1).forEach(c => f.H(c)); }],
    ['naphthalene', 'Naphthalene', 'Aromatics', f => { const h = r0(6, 6, 1.5) * Math.sqrt(3) / 2; aromaticPolycycle(f, [[-h, 0], [h, 0]], []); }],
    ['pyridine', 'Pyridine', 'Heterocycles', f => { f.ring(['N', 'C', 'C', 'C', 'C', 'C'], 1.5).slice(1).forEach(c => f.H(c)); }],
    ['pyrrole', 'Pyrrole', 'Heterocycles', f => { f.ring(['N', 'C', 'C', 'C', 'C'], 1.5).forEach(c => f.H(c)); }],
    ['furan', 'Furan', 'Heterocycles', f => { f.ring(['O', 'C', 'C', 'C', 'C'], 1.5).slice(1).forEach(c => f.H(c)); }],
    ['thiophene', 'Thiophene', 'Heterocycles', f => { f.ring(['S', 'C', 'C', 'C', 'C'], 1.5).slice(1).forEach(c => f.H(c)); }],
    ['imidazole', 'Imidazole', 'Heterocycles', f => { const r = f.ring(['N', 'C', 'N', 'C', 'C'], 1.5); [0, 1, 3, 4].forEach(k => f.H(r[k])); }],
    ['methanol', 'Methanol', 'Solvents', f => { const c = f.atom('C'); f.H(c, 3); f.to(c, 'O', 1, 1); }],
    ['ethanol', 'Ethanol', 'Solvents', f => { const c = f.atom('C'); f.H(c, 3); f.to(f.to(c, 'C', 1, 2), 'O', 1, 1); }],
    ['acetone', 'Acetone', 'Solvents', f => { const c = f.atom('C'); f.to(c, 'O', 2); f.to(c, 'C', 1, 3); f.to(c, 'C', 1, 3); }],
    ['acetonitrile', 'Acetonitrile', 'Solvents', f => { const c = f.atom('C'); f.H(c, 3); f.to(f.to(c, 'C'), 'N', 3); }],
    ['dmso', 'Dimethyl sulfoxide (DMSO)', 'Solvents', f => { const s = f.atom('S'); f.to(s, 'O', 2); f.to(s, 'C', 1, 3); f.to(s, 'C', 1, 3); }],
    ['dichloromethane', 'Dichloromethane', 'Solvents', f => { const c = f.atom('C'); f.H(c, 2); f.to(c, 'Cl'); f.to(c, 'Cl'); }],
    ['chloroform', 'Chloroform', 'Solvents', f => { const c = f.atom('C'); f.H(c, 1); for (let k = 0; k < 3; k++) f.to(c, 'Cl'); }],
    ['thf', 'Tetrahydrofuran (THF)', 'Solvents', f => { f.ring(['O', 'C', 'C', 'C', 'C'], 1, 'envelope').slice(1).forEach(c => f.H(c, 2)); }],
    ['cisplatin', 'Cisplatin, cis-[PtCl2(NH3)2]', 'Metal complexes', () => platin(true)],
    ['transplatin', 'Transplatin, trans-[PtCl2(NH3)2]', 'Metal complexes', () => platin(false)]
  ];
  function buildMolecule(fn) {
    const f = new Frag();
    const res = fn(f);
    const mol = res && res.Z ? res : f.build();
    const n = mol.Z.length, c = [0, 0, 0];
    for (let i = 0; i < n; i++) for (let k = 0; k < 3; k++) c[k] += mol.xyz[3 * i + k] / n;
    const xyz = Float64Array.from(mol.xyz, (v, k) => v - c[k % 3]);
    const out = { Z: mol.Z.slice(), xyz, bonds: mol.bonds.map(b => b.slice()) };
    if (mol.formalCharges && mol.formalCharges.some(q => q)) out.formalCharges = mol.formalCharges.slice();
    return out;
  }
  function platin(cis) {
    let m = buildCenter({ Z: [], xyz: [], bonds: [] }, { Z: 78, geometry: 'square-planar', ligand: 'Cl' }).mol;
    let r = replaceLigand(m, 1, 'NH3');          // Cl atoms 1..4 sit on sites x, y, −x, −y
    r = replaceLigand(r.mol, r.map[cis ? 2 : 3], 'NH3');
    m = r.mol;
    return { Z: m.Z, xyz: ff.optimize(m, { tol: 0.005, maxSteps: 4000 }).xyz, bonds: m.bonds };
  }

  // ------------------------------------------------------------------------------------------- groups
  // The attachment atom is atom 0; in the stored frame it sits at the origin with its parent along −x.
  const GROUP_SPECS = [
    ['methyl', 'Methyl', 'Alkyl', f => { f.H(f.atom('C'), 3); }],
    ['ethyl', 'Ethyl', 'Alkyl', f => { const c = f.atom('C'); f.H(c, 2); f.to(c, 'C', 1, 3); }],
    ['isopropyl', 'Isopropyl', 'Alkyl', f => { const c = f.atom('C'); f.H(c, 1); f.to(c, 'C', 1, 3); f.to(c, 'C', 1, 3); }],
    ['tert-butyl', 'tert-Butyl', 'Alkyl', f => { const c = f.atom('C'); for (let k = 0; k < 3; k++) f.to(c, 'C', 1, 3); }],
    ['phenyl', 'Phenyl', 'Aryl', f => { f.ring(C6, 1.5).slice(1).forEach(c => f.H(c)); }],
    ['benzyl', 'Benzyl', 'Aryl', f => { const c = f.atom('C'); f.H(c, 2); const r = f.ring(C6, 1.5); f.bond(c, r[0]); r.slice(1).forEach(x => f.H(x)); }],
    ['vinyl', 'Vinyl', 'Unsaturated', f => { const c = f.atom('C'); f.H(c, 1); f.to(c, 'C', 2, 2); }],
    ['ethynyl', 'Ethynyl', 'Unsaturated', f => { f.to(f.atom('C'), 'C', 3, 1); }],
    ['OH', 'Hydroxy', 'Oxygen', f => { f.H(f.atom('O'), 1); }],
    ['OMe', 'Methoxy', 'Oxygen', f => { f.to(f.atom('O'), 'C', 1, 3); }],
    ['NH2', 'Amino', 'Nitrogen', f => { f.H(f.atom('N'), 2); }],
    ['NMe2', 'Dimethylamino', 'Nitrogen', f => { const n = f.atom('N'); f.to(n, 'C', 1, 3); f.to(n, 'C', 1, 3); }],
    ['NO2', 'Nitro', 'Nitrogen', f => { const n = f.atom('N'); f.to(n, 'O', 2); const o = f.to(n, 'O', 1); f.fc[n] = 1; f.fc[o] = -1; }],
    ['CN', 'Cyano', 'Nitrogen', f => { f.to(f.atom('C'), 'N', 3); }],
    ['CF3', 'Trifluoromethyl', 'Halogens', f => { const c = f.atom('C'); for (let k = 0; k < 3; k++) f.to(c, 'F'); }],
    ['CHO', 'Formyl', 'Carbonyl', f => { const c = f.atom('C'); f.H(c, 1); f.to(c, 'O', 2); }],
    ['COOH', 'Carboxyl', 'Carbonyl', f => { const c = f.atom('C'); f.to(c, 'O', 2); f.to(c, 'O', 1, 1); }],
    ['COOMe', 'Methoxycarbonyl', 'Carbonyl', f => { const c = f.atom('C'); f.to(c, 'O', 2); f.to(f.to(c, 'O'), 'C', 1, 3); }],
    ['CONH2', 'Carbamoyl', 'Carbonyl', f => { const c = f.atom('C'); f.to(c, 'O', 2); f.to(c, 'N', 1, 2); }],
    ['acetyl', 'Acetyl', 'Carbonyl', f => { const c = f.atom('C'); f.to(c, 'O', 2); f.to(c, 'C', 1, 3); }],
    ['SH', 'Sulfanyl (thiol)', 'Sulfur', f => { f.H(f.atom('S'), 1); }],
    ['SO3H', 'Sulfo', 'Sulfur', f => { const s = f.atom('S'); f.to(s, 'O', 2); f.to(s, 'O', 2); f.to(s, 'O', 1, 1); }],
    ['F', 'Fluoro', 'Halogens', f => { f.atom('F'); }],
    ['Cl', 'Chloro', 'Halogens', f => { f.atom('Cl'); }],
    ['Br', 'Bromo', 'Halogens', f => { f.atom('Br'); }],
    ['I', 'Iodo', 'Halogens', f => { f.atom('I'); }]
  ];
  /** Build a group with a placeholder H as its parent, then store it with atom 0 at the origin, parent along −x. */
  function buildGroup(fn) {
    const f = new Frag();
    fn(f);
    const X = f.atom('H', [0, 0, 0]);
    f.bond(0, X);
    const mol = f.build({ frozen: [X] });
    const A = posOf(mol.xyz, 0), R = rotFromTo(sub(posOf(mol.xyz, X), A), [-1, 0, 0]);
    const keep = range(mol.Z.length).filter(i => i !== X);
    return {
      Z: keep.map(i => mol.Z[i]), local: keep.map(i => apply(R, sub(posOf(mol.xyz, i), A))),
      bonds: mol.bonds.filter(b => b[0] !== X && b[1] !== X).map(b => b.slice()),
      fc: mol.formalCharges ? keep.map(i => mol.formalCharges[i]) : null, attach: 0, dir: [-1, 0, 0]
    };
  }

  // ------------------------------------------------------------------------------------------- ligands
  // Donor atoms come first (0, or 0 and 1); the stored frame has the metal at the origin.
  function buildEn(f) {
    const d = r0(PM, 7, 1), b = 42.5 * RAD;
    const M = f.atom(PM), n1 = f.atom('N'), n2 = f.atom('N'), c1 = f.atom('C'), c2 = f.atom('C');
    f.bond(M, n1); f.bond(M, n2); f.bond(n1, c1); f.bond(c1, c2); f.bond(c2, n2);
    f.rigid([M, n1, c1, c2, n2], [[0, 0, 0], [d * Math.cos(b), d * Math.sin(b), 0], [2.85, 0.76, 0.32], [2.85, -0.76, -0.32], [d * Math.cos(b), -d * Math.sin(b), 0]]);
    f.H(n1, 2); f.H(n2, 2); f.H(c1, 2); f.H(c2, 2);
    return { M, donors: [n1, n2] };
  }
  function buildOx(f) {
    const d = r0(PM, 8, 1), b = 42 * RAD;
    const M = f.atom(PM), o1 = f.atom('O'), o2 = f.atom('O'), c1 = f.atom('C'), c2 = f.atom('C');
    f.bond(M, o1); f.bond(M, o2); f.bond(o1, c1); f.bond(c1, c2); f.bond(c2, o2);
    f.rigid([M, o1, c1, c2, o2], [[0, 0, 0], [d * Math.cos(b), d * Math.sin(b), 0], [2.79, 0.76, 0], [2.79, -0.76, 0], [d * Math.cos(b), -d * Math.sin(b), 0]]);
    f.to(c1, 'O', 2); f.to(c2, 'O', 2);
    return { M, donors: [o1, o2] };
  }
  function buildAcac(f) {
    const M = f.atom(PM), o1 = f.atom('O'), o2 = f.atom('O'), c1 = f.atom('C'), c3 = f.atom('C'), c2 = f.atom('C');
    f.bond(M, o1); f.bond(o1, c1, 2); f.bond(c1, c3); f.bond(c3, c2, 2); f.bond(c2, o2); f.bond(o2, M);
    f.rigid([M, o1, c1, c3, c2, o2], [[0, 0, 0], [1.43, 1.43, 0], [2.75, 1.25, 0], [3.45, 0, 0], [2.75, -1.25, 0], [1.43, -1.43, 0]]);
    f.to(c1, 'C', 1, 3); f.to(c2, 'C', 1, 3); f.H(c3, 1);
    return { M, donors: [o1, o2] };
  }
  function buildBipy(f) {
    const s = r0(6, 6, 1), a = r0(6, 6, 1.5);
    const vtx = (cx, ang) => [cx + a * Math.cos(ang * RAD), a * Math.sin(ang * RAD), 0];
    const sym = ['C', 'N', 'C', 'C', 'C', 'C'];
    const A = sym.map(x => f.atom(x)), B = sym.map(x => f.atom(x));
    for (let k = 0; k < 6; k++) { f.bond(A[k], A[(k + 1) % 6], 1.5); f.bond(B[k], B[(k + 1) % 6], 1.5); }
    f.bond(A[0], B[0], 1);
    // ring A: C2 at 0°, N1 at 60°, ...; ring B is its mirror image, so both N point to +y (cis, chelating)
    const pts = [...[0, 60, 120, 180, 240, 300].map(t => vtx(-s / 2 - a, t)), ...[180, 120, 60, 0, 300, 240].map(t => vtx(s / 2 + a, t))];
    const M = f.atom(PM);
    f.bond(M, A[1]); f.bond(M, B[1]);
    f.rigid([...A, ...B, M], [...pts, metalBetween(pts[1], pts[7], [0, 0, 0], r0(PM, 7, 1))]);
    [2, 3, 4, 5].forEach(k => { f.H(A[k]); f.H(B[k]); });
    return { M, donors: [A[1], B[1]] };
  }
  function buildPhen(f) {
    const s = r0(6, 6, 1.5), h = s * Math.sqrt(3) / 2, p1 = [2 * h, s], p10 = [0, 2 * s];
    const { ids, pts, unit: U } = aromaticPolycycle(f, [[0, 0], [2 * h, 0], [-h, Math.sqrt(3) * h]], [['N', p1], ['N', p10]]);
    const k1 = nearestPt(pts, p1), k10 = nearestPt(pts, p10);
    const cen = mul(pts.reduce((acc, p) => add(acc, p), [0, 0, 0]), 1 / pts.length);
    const M = f.atom(PM);
    f.bond(M, ids[k1]); f.bond(M, ids[k10]);
    U.atoms.push(M); U.local.push(metalBetween(pts[k1], pts[k10], cen, r0(PM, 7, 1)));
    return { M, donors: [ids[k1], ids[k10]] };
  }
  const LIGAND_SPECS = [
    ['H', 'Hydride', -1, f => { f.atom('H'); }],
    ['F', 'Fluoride', -1, f => { f.atom('F'); }],
    ['Cl', 'Chloride', -1, f => { f.atom('Cl'); }],
    ['Br', 'Bromide', -1, f => { f.atom('Br'); }],
    ['I', 'Iodide', -1, f => { f.atom('I'); }],
    ['OH', 'Hydroxide', -1, f => { f.H(f.atom('O'), 1); }],
    ['H2O', 'Aqua (H2O)', 0, f => { f.H(f.atom('O'), 2); }],
    ['NH3', 'Ammine (NH3)', 0, f => { f.H(f.atom('N'), 3); }],
    ['CO', 'Carbonyl (C-bound)', 0, f => { f.to(f.atom('C'), 'O', 3); }],
    ['CN', 'Cyanide (C-bound)', -1, f => { f.to(f.atom('C'), 'N', 3); }],
    ['NCMe', 'Acetonitrile (N-bound)', 0, f => { f.to(f.to(f.atom('N'), 'C', 3), 'C', 1, 3); }],
    ['PH3', 'Phosphine (PH3)', 0, f => { f.H(f.atom('P'), 3); }],
    ['PMe3', 'Trimethylphosphine', 0, f => { const p = f.atom('P'); for (let k = 0; k < 3; k++) f.to(p, 'C', 1, 3); }],
    ['PPh3', 'Triphenylphosphine', 0, f => {
      const p = f.atom('P');
      for (let k = 0; k < 3; k++) { const r = f.ring(C6, 1.5); f.bond(p, r[0]); r.slice(1).forEach(c => f.H(c)); }
    }],
    ['py', 'Pyridine (N-bound)', 0, f => { f.ring(['N', 'C', 'C', 'C', 'C', 'C'], 1.5).slice(1).forEach(c => f.H(c)); }],
    ['DMSO', 'DMSO (S-bound)', 0, f => { const s = f.atom('S'); f.to(s, 'O', 2); f.to(s, 'C', 1, 3); f.to(s, 'C', 1, 3); }],
    ['X', 'Dummy site marker (X)', 0, null],
    ['en', 'Ethylenediamine (en)', 0, buildEn, 2],
    ['bipy', '2,2′-Bipyridine (bipy)', 0, buildBipy, 2],
    ['ox', 'Oxalate (ox)', -2, buildOx, 2],
    ['acac', 'Acetylacetonate (acac)', -1, buildAcac, 2],
    ['phen', '1,10-Phenanthroline (phen)', 0, buildPhen, 2]
  ];
  /** Build a ligand around a placeholder metal; store it with donors first and the metal at the origin. */
  function buildLigand(fn, dent) {
    if (!fn) return { Z: [0], local: [[1, 0, 0]], bonds: [], donors: [0], fc: null };
    const f = new Frag();
    let M, donors;
    if (dent === 2) ({ M, donors } = fn(f));
    else { fn(f); donors = [0]; M = f.atom(PM, [0, 0, 0]); f.bond(M, 0); }
    const mol = f.build(dent === 2 ? { frozen: [M], geometries: { [M]: 'square-planar' } } : { frozen: [M] });
    const order = [...donors, ...range(mol.Z.length).filter(i => i !== M && donors.indexOf(i) < 0)];
    const inv = new Int32Array(mol.Z.length).fill(-1);
    order.forEach((o, k) => { inv[o] = k; });
    const m0 = posOf(mol.xyz, M);
    let pts = order.map(i => sub(posOf(mol.xyz, i), m0));
    let R;
    if (dent === 2) {
      const e1 = unit(add(unit(pts[0]), unit(pts[1]))), e2 = unit(perpPart(sub(pts[0], pts[1]), e1));
      R = frameRotation(e1, e2, [1, 0, 0], [0, 1, 0]);
    } else R = rotFromTo(pts[0], [1, 0, 0]);
    pts = pts.map(p => apply(R, p));
    const bonds = mol.bonds.filter(b => b[0] !== M && b[1] !== M).map(b => [inv[b[0]], inv[b[1]], b[2]]);
    return { Z: order.map(i => mol.Z[i]), local: pts, bonds, donors: dent === 2 ? [0, 1] : [0], fc: null };
  }

  // ------------------------------------------------------------------------------------------- public lists
  function lazyGetter(entry, build, publicView) {
    let cache = null;
    const get = () => (cache || (cache = build()));
    Object.defineProperty(entry, '_data', { enumerable: false, get });
    Object.defineProperty(entry, 'mol', { enumerable: false, get: () => publicView(get()) });
  }
  const copyMol = m => {
    const o = { Z: m.Z.slice(), xyz: Float64Array.from(m.xyz), bonds: m.bonds.map(b => b.slice()) };
    if (m.formalCharges) o.formalCharges = m.formalCharges.slice();
    return o;
  };
  const fragView = d => {
    const o = { Z: d.Z.slice(), xyz: Float64Array.from([].concat(...d.local)), bonds: d.bonds.map(b => b.slice()) };
    if (d.fc && d.fc.some(q => q)) o.formalCharges = d.fc.slice();
    return o;
  };
  const geometries = ff.GEOMETRY_IDS.map(id => ({
    id, name: GEO[id].name, category: 'Coordination geometries', kind: 'geometry', cn: GEO[id].cn, vectors: GEO[id].vectors.map(v => v.slice())
  }));
  const molecules = MOLECULE_SPECS.map(([id, name, category, fn]) => {
    const e = { id, name, category, kind: 'molecule' };
    lazyGetter(e, () => buildMolecule(fn), copyMol);
    return e;
  });
  const groups = GROUP_SPECS.map(([id, name, category, fn]) => {
    const e = { id, name, category, kind: 'group', attach: 0, direction: [-1, 0, 0] };
    lazyGetter(e, () => buildGroup(fn), fragView);
    return e;
  });
  const ligands = LIGAND_SPECS.map(([id, name, charge, fn, dent]) => {
    dent = dent || 1;
    const e = { id, name, category: dent === 2 ? 'Bidentate ligands' : 'Monodentate ligands', kind: 'ligand', denticity: dent, charge, donors: dent === 2 ? [0, 1] : [0] };
    lazyGetter(e, () => buildLigand(fn, dent), d => Object.assign(fragView(d), { donors: d.donors.slice(), charge }));
    return e;
  });
  const INDEX = { geometry: geometries, molecule: molecules, group: groups, ligand: ligands };
  /** Look up a library entry by kind ('geometry' | 'molecule' | 'group' | 'ligand') and id (case-insensitive). */
  function get(kind, id) {
    const list = INDEX[kind];
    if (!list) throw new Error('templates: unknown kind "' + kind + '".');
    const key = String(id);
    return list.find(e => e.id === key) || list.find(e => e.id.toLowerCase() === key.toLowerCase()) || null;
  }
  function need(kind, id) {
    const e = get(kind, id);
    if (!e) throw new Error('templates: unknown ' + kind + ' "' + id + '".');
    return e;
  }

  // ------------------------------------------------------------------------------------------- editing helpers
  const DUMMY_EXTRA = 0.70;    // a dummy site marker sits at covalent(M) + 0.70 Å, like a first-row donor
  function normMol(mol) {
    mol = mol || {};
    const Z = Array.from(mol.Z || [], z => (z > 0 ? z | 0 : 0)), n = Z.length;
    const xyz = Float64Array.from(mol.xyz || []);
    if (xyz.length !== 3 * n) throw new Error('templates: mol.xyz must hold 3 coordinates per atom.');
    const bonds = (mol.bonds || []).map(b => [b[0] | 0, b[1] | 0, b[2] === undefined ? 1 : +b[2]]);
    for (const b of bonds) if (b[0] < 0 || b[1] < 0 || b[0] >= n || b[1] >= n) throw new Error('templates: bond refers to a missing atom.');
    const fc = mol.formalCharges ? Array.from(mol.formalCharges, q => +q || 0) : null;
    return { Z, xyz, bonds, fc };
  }
  function result(M, extra) {
    const mol = { Z: M.Z, xyz: M.xyz, bonds: M.bonds };
    if (M.fc) mol.formalCharges = M.fc;
    return Object.assign({ mol }, extra);
  }
  function adjacencyAll(M) {
    const nbr = M.Z.map(() => []);
    for (const [i, j] of M.bonds) if (i !== j) { nbr[i].push(j); nbr[j].push(i); }
    return nbr;
  }
  function removeAtoms(M, set) {
    const n = M.Z.length, map = new Int32Array(n);
    let k = 0;
    for (let i = 0; i < n; i++) map[i] = set.has(i) ? -1 : k++;
    const Z = [], xyz = new Float64Array(3 * k), fc = M.fc ? [] : null;
    for (let i = 0; i < n; i++) {
      if (map[i] < 0) continue;
      Z.push(M.Z[i]);
      xyz[3 * map[i]] = M.xyz[3 * i]; xyz[3 * map[i] + 1] = M.xyz[3 * i + 1]; xyz[3 * map[i] + 2] = M.xyz[3 * i + 2];
      if (fc) fc.push(M.fc[i]);
    }
    const bonds = M.bonds.filter(b => map[b[0]] >= 0 && map[b[1]] >= 0).map(b => [map[b[0]], map[b[1]], b[2]]);
    return { M: { Z, xyz, bonds, fc }, map };
  }
  /** Append atoms (in place on the internal copy M); returns their new indices. */
  function appendAtoms(M, Zs, pts, bonds, fcs) {
    const start = M.Z.length, m = Zs.length, xyz = new Float64Array(3 * (start + m));
    xyz.set(M.xyz);
    pts.forEach((p, k) => { xyz[3 * (start + k)] = p[0]; xyz[3 * (start + k) + 1] = p[1]; xyz[3 * (start + k) + 2] = p[2]; });
    M.Z = M.Z.concat(Zs); M.xyz = xyz;
    M.bonds = M.bonds.concat(bonds.map(b => [b[0] + start, b[1] + start, b[2]]));
    if (M.fc || (fcs && fcs.some(q => q))) M.fc = (M.fc || new Array(start).fill(0)).concat(fcs ? fcs.slice() : new Array(m).fill(0));
    return range(m).map(k => start + k);
  }
  /** Short force-field relaxation of the new atoms only; everything else is frozen. */
  function relax(M, added, geometries) {
    if (!added.length) return;
    const n = M.Z.length, isNew = new Uint8Array(n), frozen = [];
    for (const i of added) isNew[i] = 1;
    for (let i = 0; i < n; i++) if (!isNew[i]) frozen.push(i);
    const mol = { Z: M.Z, xyz: M.xyz, bonds: M.bonds };
    if (M.fc) mol.formalCharges = M.fc;
    M.xyz = ff.optimize(mol, { frozen, geometries, maxSteps: 500, tol: 0.05 }).xyz;
  }
  function nearAtoms(M, center, radius, skip) {
    const out = [], r2 = radius * radius;
    for (let i = 0; i < M.Z.length; i++) {
      if (!(M.Z[i] > 0) || (skip && skip.has(i))) continue;
      const dx = M.xyz[3 * i] - center[0], dy = M.xyz[3 * i + 1] - center[1], dz = M.xyz[3 * i + 2] - center[2];
      if (dx * dx + dy * dy + dz * dz < r2) out.push(i);
    }
    return out;
  }
  const fragRadius = d => d.local.reduce((m, p) => Math.max(m, len(p)), 0);
  /** Smooth clash measure Σ (σ/r)^8 between new points and nearby existing atoms. */
  function clashScore(pts, Zs, M, near) {
    let s = 0;
    for (let k = 0; k < pts.length; k++) {
      if (!(Zs[k] > 0)) continue;
      const rk = el.vdw[Zs[k]], p = pts[k];
      for (const j of near) {
        const dx = p[0] - M.xyz[3 * j], dy = p[1] - M.xyz[3 * j + 1], dz = p[2] - M.xyz[3 * j + 2];
        const sig = 0.8 * (rk + el.vdw[M.Z[j]]), t = sig * sig / Math.max(dx * dx + dy * dy + dz * dz, 1e-6);
        s += t * t * t * t;
      }
    }
    return s;
  }
  /** Rigid placement: anchor atom at `at`, local direction `ldir` turned onto `wdir`, best twist about wdir. */
  function orient(local, anchor, ldir, at, wdir, Zs, M, near) {
    const R1 = rotFromTo(ldir, wdir), rel = local.map(p => apply(R1, sub(p, local[anchor])));
    let best = null;
    for (let s = 0; s < (local.length > 1 ? 36 : 1); s++) {
      const R = rotAxis(wdir, s * Math.PI / 18), pts = rel.map(p => add(at, apply(R, p)));
      const sc = clashScore(pts, Zs, M, near);
      if (!best || sc < best.score - 1e-12) best = { pts, score: sc };
    }
    return best;
  }
  /** Best-fitting polyhedron with exactly as many sites as the metal has neighbours now. */
  function currentGeometry(M, m, nbr) {
    const pm = posOf(M.xyz, m), dirs = nbr[m].map(j => unit(sub(posOf(M.xyz, j), pm)));
    if (dirs.length < 2 || dirs.some(d => len(d) === 0)) return null;
    let best = null, bs = Infinity;
    for (const id of ff.GEOMETRY_IDS) {
      if (GEO[id].cn !== dirs.length) continue;
      const r = ff.assignSites(dirs, GEO[id].vectors);
      if (r && r.score < bs) { bs = r.score; best = id; }
    }
    return best;
  }
  /** Directions of open sites on atom a for one more substituent (VSEPR or coordination polyhedron). */
  function openSites(M, a, nbr) {
    const n = M.Z.length, have = nbr[a], pa = posOf(M.xyz, a);
    const dirs = have.map(j => unit(sub(posOf(M.xyz, j), pa))).filter(d => len(d) > 0);
    let sites = null;
    if (ff.usesPolyhedron(M.Z[a])) {
      const g = ff.metalGeometry(M.Z[a], have.length + 1, dirs.length ? dirs : null);
      sites = g ? GEO[g].vectors : null;
    } else {
      const tmp = { Z: M.Z.concat([1]), xyz: null, bonds: M.bonds.concat([[a, n, 1]]), formalCharges: M.fc ? M.fc.concat([0]) : undefined };
      sites = siteSet(ff.perceive(tmp).atoms[a]);
    }
    let ref = null;
    if (have.length === 1) {
      const p = have[0], r = nbr[p].find(j => j !== a && M.Z[j] > 1) ?? nbr[p].find(j => j !== a);
      if (r !== undefined) ref = sub(posOf(M.xyz, r), posOf(M.xyz, p));
    }
    let free = sites && sites.length > dirs.length ? openDirections(sites, dirs, ref) : [];
    if (!free.length) free = spreadDirections(dirs, 1);
    return free;
  }
  /** The ligand attached through site atom s: everything reachable from s without passing the metal m
   * (other metals are never entered). */
  function ligandAtoms(M, nbr, s, m) {
    const seen = new Set([s]), stack = [s];
    while (stack.length) {
      const u = stack.pop();
      for (const v of nbr[u]) {
        if (v === m || seen.has(v) || ff.isMetal(M.Z[v])) continue;
        seen.add(v); stack.push(v);
      }
    }
    return seen;
  }
  function metalOf(M, nbr, s, fn) {
    if (!Number.isInteger(s) || s < 0 || s >= M.Z.length) throw new Error(fn + ': atom index ' + s + ' is out of range.');
    if (ff.isMetal(M.Z[s])) throw new Error(fn + ': atom ' + (s + 1) + ' is a metal; pick a ligand atom bonded to it.');
    const ps = posOf(M.xyz, s), metals = nbr[s].filter(j => ff.isMetal(M.Z[j]));
    if (!metals.length) throw new Error(fn + ': atom ' + (s + 1) + ' is not bonded to a metal.');
    return metals.reduce((b, j) => (len(sub(posOf(M.xyz, j), ps)) < len(sub(posOf(M.xyz, b), ps)) ? j : b));
  }
  function placeMono(M, m, u, lig) {
    const d = lig._data, Zm = M.Z[m], ZD = d.Z[0];
    const L = ZD > 0 ? r0(Zm, ZD, 1) : el.covalent[Zm] + DUMMY_EXTRA;
    const pm = posOf(M.xyz, m), at = add(pm, mul(u, L));
    const near = nearAtoms(M, pm, L + fragRadius(d) + 4, null);
    const o = orient(d.local, 0, mul(unit(d.local[0]), -1), at, mul(u, -1), d.Z, M, near);
    const added = appendAtoms(M, d.Z, o.pts, d.bonds, d.fc);
    M.bonds.push([m, added[0], 1]);
    return added;
  }
  function placeBidentate(M, m, uA, uB, lig) {
    const d = lig._data, Zm = M.Z[m], pm = posOf(M.xyz, m);
    const bis = add(unit(uA), unit(uB));
    if (len(bis) < 0.2) throw new Error('addBidentate: the two sites must be cis to each other.');
    const PA = add(pm, mul(uA, r0(Zm, d.Z[0], 1))), PB = add(pm, mul(uB, r0(Zm, d.Z[1], 1)));
    const e1t = unit(bis), e2t = unit(perpPart(sub(PA, PB), e1t));
    const L1 = d.local[0], L2 = d.local[1];
    const e1l = unit(add(unit(L1), unit(L2))), e2l = unit(perpPart(sub(L1, L2), e1l));
    const shift = dot(sub(mul(add(PA, PB), 0.5), pm), e1t) - dot(mul(add(L1, L2), 0.5), e1l);
    const near = nearAtoms(M, pm, fragRadius(d) + 4, null);
    let best = null;
    for (const flip of [1, -1]) {         // a half turn about the bisector swaps which donor takes which site
      const R = frameRotation(e1l, e2l, e1t, mul(e2t, flip));
      const pts = d.local.map(p => add(add(pm, apply(R, p)), mul(e1t, shift)));
      const sc = clashScore(pts, d.Z, M, near);
      if (!best || sc < best.score) best = { pts, score: sc };
    }
    const added = appendAtoms(M, d.Z, best.pts, d.bonds, d.fc);
    M.bonds.push([m, added[0], 1], [m, added[1], 1]);
    return added;
  }
  /** Largest set of disjoint cis site pairs (|cos| < 0.6). */
  function cisPairs(vecs) {
    const n = vecs.length, used = new Uint8Array(n), cur = [];
    let best = [];
    (function rec(i) {
      if (best.length === Math.floor(n / 2)) return;
      while (i < n && used[i]) i++;
      if (i >= n) { if (cur.length > best.length) best = cur.slice(); return; }
      used[i] = 1;
      for (let j = i + 1; j < n; j++) {
        if (used[j] || Math.abs(dot(vecs[i], vecs[j])) >= 0.6) continue;
        used[j] = 1; cur.push([i, j]); rec(i + 1); cur.pop(); used[j] = 0;
      }
      rec(i + 1);
      used[i] = 0;
    })(0);
    return best;
  }

  // ------------------------------------------------------------------------------------------- operations
  /** Insert a library molecule centred (centroid) at opts.at. */
  function insertMolecule(mol, id, opts) {
    const t = need('molecule', id)._data, M = normMol(mol);
    const at = opts && opts.at ? [+opts.at[0] || 0, +opts.at[1] || 0, +opts.at[2] || 0] : [0, 0, 0];
    const pts = range(t.Z.length).map(i => add(posOf(t.xyz, i), at));
    const added = appendAtoms(M, t.Z, pts, t.bonds, t.formalCharges || null);
    relax(M, added);
    return result(M, { added });
  }

  /** Attach a substituent: an H (or dummy) is replaced along its bond; a heavy atom gets it on an open site. */
  function attachGroup(mol, atomIndex, groupId) {
    const g = need('group', groupId)._data, M0 = normMol(mol), n = M0.Z.length, a = atomIndex;
    if (!Number.isInteger(a) || a < 0 || a >= n) throw new Error('attachGroup: atom index ' + atomIndex + ' is out of range.');
    const nbr = adjacencyAll(M0), pa = posOf(M0.xyz, a);
    let parent = -1, cands = null;
    const removed = [];
    if (M0.Z[a] <= 1) {
      removed.push(a);
      if (nbr[a].length) {
        parent = nbr[a].reduce((b, j) => (len(sub(posOf(M0.xyz, j), pa)) < len(sub(posOf(M0.xyz, b), pa)) ? j : b));
        const u = unit(sub(pa, posOf(M0.xyz, parent)));
        cands = [len(u) > 0 ? u : [1, 0, 0]];
      }
    } else { parent = a; cands = openSites(M0, a, nbr); }
    const skip = new Set(removed);
    let best = null;
    if (parent < 0) best = { pts: g.local.map(p => add(pa, p)) };
    else {
      const pp = posOf(M0.xyz, parent), L = r0(M0.Z[parent], g.Z[0], 1);
      const near = nearAtoms(M0, pp, L + fragRadius(g) + 4, skip);
      for (const u of cands) {
        const o = orient(g.local, 0, g.dir, add(pp, mul(u, L)), mul(u, -1), g.Z, M0, near);
        if (!best || o.score < best.score) best = o;
      }
    }
    let M = M0, map = null;
    if (removed.length) ({ M, map } = removeAtoms(M0, skip));
    const p2 = parent >= 0 ? (map ? map[parent] : parent) : -1;
    const added = appendAtoms(M, g.Z, best.pts, g.bonds, g.fc);
    if (p2 >= 0) M.bonds.push([p2, added[0], 1]);
    let hint;
    if (p2 >= 0 && ff.usesPolyhedron(M.Z[p2])) { const gid = currentGeometry(M, p2, adjacencyAll(M)); if (gid) hint = { [p2]: gid }; }
    relax(M, added, hint);
    const extra = { added };
    if (removed.length) { extra.removed = removed; extra.map = map; }
    return result(M, extra);
  }

  /** A metal atom (opts.Z: number or symbol) at opts.at with one ligand per site of opts.geometry.
   * opts.ligand defaults to 'Cl'; 'X' puts dummy site markers; a bidentate ligand fills cis pairs. */
  function buildCenter(mol, opts) {
    opts = opts || {};
    const M = normMol(mol);
    const Zm = typeof opts.Z === 'string' ? el.toZ(opts.Z) : Math.floor(+opts.Z);
    if (!(Zm > 0) || Zm >= el.count) throw new Error('buildCenter: unknown element "' + opts.Z + '".');
    const gid = opts.geometry || (ff.metalGeometry(Zm, 4, null) === 'square-planar' ? 'square-planar' : 'octahedral');
    const geom = GEO[gid];
    if (!geom) throw new Error('buildCenter: unknown geometry "' + gid + '".');
    const lig = need('ligand', opts.ligand === undefined ? 'Cl' : opts.ligand);
    const at = opts.at ? [+opts.at[0] || 0, +opts.at[1] || 0, +opts.at[2] || 0] : [0, 0, 0];
    const [m] = appendAtoms(M, [Zm], [at], [], null);
    const added = [m];
    if (lig.denticity === 2) {
      const pairs = cisPairs(geom.vectors), paired = new Set([].concat(...pairs));
      for (const [p, q] of pairs) added.push(...placeBidentate(M, m, geom.vectors[p], geom.vectors[q], lig));
      geom.vectors.forEach((v, k) => { if (!paired.has(k)) added.push(...placeMono(M, m, v, need('ligand', 'X'))); });
    } else for (const v of geom.vectors) added.push(...placeMono(M, m, v, lig));
    relax(M, added.slice(1), { [m]: gid });       // the metal stays exactly at `at`
    return result(M, { added });
  }

  /** Replace the ligand bound through site atom siteAtomIndex by a monodentate library ligand. */
  function replaceLigand(mol, siteAtomIndex, ligandId) {
    const lig = need('ligand', ligandId);
    if (lig.denticity !== 1) throw new Error('replaceLigand: "' + ligandId + '" is bidentate; use addBidentate.');
    const M0 = normMol(mol), nbr = adjacencyAll(M0), s = siteAtomIndex;
    const m = metalOf(M0, nbr, s, 'replaceLigand');
    const hint = currentGeometry(M0, m, nbr);
    const u = unit(sub(posOf(M0.xyz, s), posOf(M0.xyz, m)));
    const gone = ligandAtoms(M0, nbr, s, m);
    const { M, map } = removeAtoms(M0, gone);
    const m2 = map[m];
    const added = placeMono(M, m2, u, lig);
    relax(M, added, hint ? { [m2]: hint } : undefined);
    return result(M, { added, removed: [...gone].sort((x, y) => x - y), map });
  }

  /** Put a bidentate ligand on two cis sites (site atoms bonded to the same metal). */
  function addBidentate(mol, siteA, siteB, ligandId) {
    const lig = need('ligand', ligandId);
    if (lig.denticity !== 2) throw new Error('addBidentate: "' + ligandId + '" is not a bidentate ligand.');
    const M0 = normMol(mol), nbr = adjacencyAll(M0);
    if (siteA === siteB) throw new Error('addBidentate: pick two different site atoms.');
    const m = metalOf(M0, nbr, siteA, 'addBidentate');
    if (nbr[siteB].indexOf(m) < 0) throw new Error('addBidentate: the two site atoms are not bonded to the same metal.');
    const pm = posOf(M0.xyz, m);
    const uA = unit(sub(posOf(M0.xyz, siteA), pm)), uB = unit(sub(posOf(M0.xyz, siteB), pm));
    if (Math.acos(Math.max(-1, Math.min(1, dot(uA, uB)))) / RAD > 135) throw new Error('addBidentate: the two sites must be cis to each other.');
    const hint = currentGeometry(M0, m, nbr);
    const gone = ligandAtoms(M0, nbr, siteA, m);
    for (const i of ligandAtoms(M0, nbr, siteB, m)) gone.add(i);
    const { M, map } = removeAtoms(M0, gone);
    const m2 = map[m];
    const added = placeBidentate(M, m2, uA, uB, lig);
    relax(M, added, hint ? { [m2]: hint } : undefined);
    return result(M, { added, removed: [...gone].sort((x, y) => x - y), map });
  }

  /** Pairs of atoms at least three bonds apart (or unconnected) that are closer than scale·(Rvdw,i + Rvdw,j).
   * Returns [[i, j, distance], ...]; dummies are ignored.  opts.scale defaults to 0.7. */
  function findClashes(mol, opts) {
    const M = normMol(mol), n = M.Z.length, scale = opts && opts.scale > 0 ? opts.scale : 0.7;
    const nbr = adjacencyAll(M), close = new Set();
    for (let i = 0; i < n; i++) for (const j of nbr[i]) { close.add(i * n + j); for (const k of nbr[j]) close.add(i * n + k); }
    const out = [];
    for (let i = 0; i < n; i++) {
      if (!(M.Z[i] > 0)) continue;
      for (let j = i + 1; j < n; j++) {
        if (!(M.Z[j] > 0) || close.has(i * n + j)) continue;
        const dx = M.xyz[3 * j] - M.xyz[3 * i], dy = M.xyz[3 * j + 1] - M.xyz[3 * i + 1], dz = M.xyz[3 * j + 2] - M.xyz[3 * i + 2];
        const d = Math.sqrt(dx * dx + dy * dy + dz * dz), lim = scale * (el.vdw[M.Z[i]] + el.vdw[M.Z[j]]);
        if (d < lim) out.push([i, j, d]);
      }
    }
    return out;
  }

  BM.templates = {
    geometries, ligands, groups, molecules, get,
    insertMolecule, attachGroup, buildCenter, replaceLigand, addBidentate, findClashes
  };
})(typeof window !== 'undefined' ? (window.BM = window.BM || {}) : (globalThis.BM = globalThis.BM || {}));
