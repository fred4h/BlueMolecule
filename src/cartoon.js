/* BlueMolecule — macromolecule geometry.
 *   Protein cartoon: Catmull-Rom spline through CA, ribbon width along the carbonyl C→O direction; helices as ribbons
 *   (or cylinders), sheets as arrows, coil as tubes. Secondary structure from ds.secondaryStructure (PDB HELIX/SHEET,
 *   mmCIF struct_conf) or a simplified DSSP (Kabsch–Sander backbone H-bond energies, amide H from the backbone).
 *   Nucleic acids: backbone tube through P (C4′ where P is missing) with base slabs on rods. Trace: CA/P tube.
 *   Rings: filled ring polygons (sugars coloured by monosaccharide). SNFG: 3D symbols at sugar ring centroids, in the
 *   ring plane, with link sticks between connected residues.
 *
 *   BM.cartoon.secondary(ctx)             -> Uint8Array per residue: 0 coil, 1 helix, 2 sheet
 *   BM.cartoon.dssp(ctx)                  -> same, always computed from the coordinates
 *   BM.cartoon.build(kind, ctx, mask, colors, opts) -> {meshes: [{positions, normals, colors, indices, twoSided}], proxies, proxyRadius}
 *       kind: 'cartoon' | 'trace' | 'rings' | 'snfg'; colors: Float32Array(3n) per atom;
 *       opts: {cylinders, scheme, scale (size factor), profile (array: receives the cartoon cross-sections, for tests)}
 * Pure geometry code (browser and Node). */
(function (BM) {
  'use strict';
  const reps = () => BM.reps;
  const C_PROTEIN = 1, C_NUCLEIC = 2, C_SUGAR = 6;

  // ------------------------------------------------------------------ small vector helpers
  const sub = (a, b) => [a[0] - b[0], a[1] - b[1], a[2] - b[2]];
  const add = (a, b) => [a[0] + b[0], a[1] + b[1], a[2] + b[2]];
  const mul = (a, s) => [a[0] * s, a[1] * s, a[2] * s];
  const dot = (a, b) => a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
  const cross = (a, b) => [a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]];
  const len = a => Math.hypot(a[0], a[1], a[2]);
  const norm = a => { const l = len(a) || 1; return [a[0] / l, a[1] / l, a[2] / l]; };
  const lerp = (a, b, t) => [a[0] + (b[0] - a[0]) * t, a[1] + (b[1] - a[1]) * t, a[2] + (b[2] - a[2]) * t];
  const perp = a => norm(Math.abs(a[0]) < 0.9 ? cross(a, [1, 0, 0]) : cross(a, [0, 1, 0]));
  const P = (xyz, i) => [xyz[3 * i], xyz[3 * i + 1], xyz[3 * i + 2]];
  const colAt = (colors, i) => [colors[3 * i], colors[3 * i + 1], colors[3 * i + 2]];

  // ------------------------------------------------------------------ mesh builder (growable typed arrays)
  class Mesh {
    constructor(twoSided) { this.p = new Float32Array(3 * 4096); this.nr = new Float32Array(3 * 4096); this.c = new Float32Array(3 * 4096); this.ix = new Uint32Array(3 * 4096); this.nv = 0; this.ni = 0; this.twoSided = !!twoSided; }
    _grow(kind, need) {
      if (kind === 'v') {
        const cap = Math.max(this.p.length * 2, need);
        for (const k of ['p', 'nr', 'c']) { const a = new Float32Array(cap); a.set(this[k]); this[k] = a; }
      } else { const a = new Uint32Array(Math.max(this.ix.length * 2, need)); a.set(this.ix); this.ix = a; }
    }
    v(p, n, c) {
      const k = 3 * this.nv;
      if (k + 3 > this.p.length) this._grow('v', k + 3);
      this.p[k] = p[0]; this.p[k + 1] = p[1]; this.p[k + 2] = p[2];
      this.nr[k] = n[0]; this.nr[k + 1] = n[1]; this.nr[k + 2] = n[2];
      this.c[k] = c[0]; this.c[k + 1] = c[1]; this.c[k + 2] = c[2];
      return this.nv++;
    }
    t(a, b, c) { if (this.ni + 3 > this.ix.length) this._grow('i', this.ni + 3); this.ix[this.ni++] = a; this.ix[this.ni++] = b; this.ix[this.ni++] = c; }
    flatTri(a, b, c, col) { const n = norm(cross(sub(b, a), sub(c, a))); const i = this.v(a, n, col), j = this.v(b, n, col), k = this.v(c, n, col); this.t(i, j, k); }
    flatQuad(a, b, c, d, col) { const n = norm(cross(sub(b, a), sub(c, a))); const i = this.v(a, n, col), j = this.v(b, n, col), k = this.v(c, n, col), l = this.v(d, n, col); this.t(i, j, k); this.t(i, k, l); }
    done() {
      if (!this.ni) return null;
      return { positions: this.p.slice(0, 3 * this.nv), normals: this.nr.slice(0, 3 * this.nv), colors: this.c.slice(0, 3 * this.nv), indices: this.ix.slice(0, this.ni), twoSided: this.twoSided };
    }
  }
  /** Capped cylinder a→b (radius r, seg sides) with colour c. */
  function cylinder(m, a, b, r, c, seg, capped) {
    const ax = sub(b, a), L = len(ax);
    if (!(L > 1e-6)) return;
    const T = mul(ax, 1 / L), U = perp(T), V = cross(T, U), ring = [];
    for (let k = 0; k <= seg; k++) {
      const th = 2 * Math.PI * k / seg, d = add(mul(U, Math.cos(th)), mul(V, Math.sin(th)));
      ring.push([m.v(add(a, mul(d, r)), d, c), m.v(add(b, mul(d, r)), d, c)]);
    }
    for (let k = 0; k < seg; k++) { const [a0, b0] = ring[k], [a1, b1] = ring[k + 1]; m.t(a0, a1, b1); m.t(a0, b1, b0); }
    if (capped) for (const [p, s] of [[a, -1], [b, 1]]) {
      const nn = mul(T, s), ctr = m.v(p, nn, c), idx = [];
      for (let k = 0; k <= seg; k++) { const th = 2 * Math.PI * k / seg; idx.push(m.v(add(p, mul(add(mul(U, Math.cos(th)), mul(V, Math.sin(th))), r)), nn, c)); }
      for (let k = 0; k < seg; k++) { if (s > 0) m.t(ctr, idx[k], idx[k + 1]); else m.t(ctr, idx[k + 1], idx[k]); }
    }
  }

  // ------------------------------------------------------------------ backbone extraction (cached per topology)
  const O_NAMES = ['O', 'OT1', 'OC1', 'O1', 'OXT'];
  function proteinChains(ctx) {
    const cls = reps().classify(ctx);
    if (cls.bbProt) return cls.bbProt;
    const names = ctx.frame.atomNames, xyz = ctx.xyz, Z = ctx.Z;
    const res = [], ca = [], nn = [], cc = [], oo = [], seg = [];
    if (names) {
      for (let r = 0; r < cls.nRes; r++) {
        if (cls.resClass[r] !== C_PROTEIN) continue;
        let a = -1, n = -1, c = -1, o = -1, oRank = 9;
        for (let t = cls.resOff[r]; t < cls.resOff[r + 1]; t++) {
          const i = cls.resAtoms[t], nm = String(names[i] || '').trim().toUpperCase();
          if (nm === 'CA' && Z[i] === 6) a = i; else if (nm === 'N') n = i; else if (nm === 'C') c = i;
          else { const k = O_NAMES.indexOf(nm); if (k >= 0 && k < oRank) { o = i; oRank = k; } }
        }
        if (a < 0) continue;
        res.push(r); ca.push(a); nn.push(n); cc.push(c); oo.push(o);
      }
    }
    let s = -1;
    for (let k = 0; k < res.length; k++) {
      let brk = k === 0 || cls.chain[res[k]] !== cls.chain[res[k - 1]];
      if (!brk) {
        const d = len(sub(P(xyz, ca[k]), P(xyz, ca[k - 1])));
        if (d > 4.3) brk = true;
        else if (cc[k - 1] >= 0 && nn[k] >= 0 && len(sub(P(xyz, nn[k]), P(xyz, cc[k - 1]))) > 2.2) brk = true;
      }
      if (brk) s++;
      seg.push(s);
    }
    cls.bbProt = { count: res.length, res: Int32Array.from(res), ca: Int32Array.from(ca), n: Int32Array.from(nn), c: Int32Array.from(cc), o: Int32Array.from(oo), seg: Int32Array.from(seg) };
    return cls.bbProt;
  }
  const PUR6 = ['N1', 'C2', 'N3', 'C4', 'C5', 'C6'], PUR5 = ['C4', 'C5', 'N7', 'C8', 'N9'], PURP = ['N1', 'C2', 'N3', 'C4', 'N9', 'C8', 'N7', 'C5', 'C6'];
  function nucleicChains(ctx) {
    const cls = reps().classify(ctx);
    if (cls.bbNuc) return cls.bbNuc;
    const names = ctx.frame.atomNames, xyz = ctx.xyz;
    const res = [], ctl = [], seg = [], bases = [];
    if (names) {
      for (let r = 0; r < cls.nRes; r++) {
        if (cls.resClass[r] !== C_NUCLEIC) continue;
        const at = new Map();
        for (let t = cls.resOff[r]; t < cls.resOff[r + 1]; t++) { const i = cls.resAtoms[t]; at.set(String(names[i] || '').trim().toUpperCase().replace(/\*/g, "'"), i); }
        const c = at.has('P') ? at.get('P') : at.has("C4'") ? at.get("C4'") : at.has("C3'") ? at.get("C3'") : -1;
        if (c < 0) continue;
        const pur = at.has('N9') && at.has('C8');
        const ring = (pur ? PURP : PUR6).map(x => at.get(x));
        const b = ring.every(x => x !== undefined) ? { pur, perim: ring, r6: PUR6.map(x => at.get(x)), r5: pur ? PUR5.map(x => at.get(x)) : null, attach: at.get(pur ? 'N9' : 'N1') } : null;
        res.push(r); ctl.push(c); bases.push(b);
      }
    }
    let s = -1;
    for (let k = 0; k < res.length; k++) {
      let brk = k === 0 || cls.chain[res[k]] !== cls.chain[res[k - 1]];
      if (!brk && len(sub(P(xyz, ctl[k]), P(xyz, ctl[k - 1]))) > 8.5) brk = true;
      if (brk) s++;
      seg.push(s);
    }
    cls.bbNuc = { count: res.length, res: Int32Array.from(res), ctl: Int32Array.from(ctl), seg: Int32Array.from(seg), bases };
    return cls.bbNuc;
  }

  // ------------------------------------------------------------------ secondary structure
  /** Per-residue secondary structure: file records when present, otherwise the simplified DSSP. */
  function secondary(ctx) {
    const cls = reps().classify(ctx), ss = new Uint8Array(cls.nRes), bb = proteinChains(ctx);
    if (!bb.count) return ss;
    const rec = ctx.ssRecords;
    if (rec && rec.length) {
      const hasChains = cls.chain.some(c => c);
      for (const R of rec) {
        const code = R.type === 'helix' ? 1 : R.type === 'sheet' ? 2 : 0;
        const a = R.start && R.start.resId, b = R.end && R.end.resId;
        if (!code || !Number.isFinite(a) || !Number.isFinite(b)) continue;
        for (let k = 0; k < bb.count; k++) {
          const r = bb.res[k], id = cls.resId[r];
          if (id >= a && id <= b && (!hasChains || !R.chain || cls.chain[r] === String(R.chain).trim())) ss[r] = code;
        }
      }
      cls.ssSource = 'file';
      return ss;
    }
    const codes = dsspCodes(ctx, bb);
    for (let k = 0; k < bb.count; k++) ss[bb.res[k]] = codes[k];
    cls.ssSource = 'dssp';
    return ss;
  }
  function dssp(ctx) {
    const cls = reps().classify(ctx), ss = new Uint8Array(cls.nRes), bb = proteinChains(ctx);
    const codes = dsspCodes(ctx, bb);
    for (let k = 0; k < bb.count; k++) ss[bb.res[k]] = codes[k];
    return ss;
  }
  /** Kabsch–Sander electrostatic H-bond energy (kcal/mol) between C=O (c, o) and N-H (n, h). */
  function hbEnergy(c, o, n, h) {
    const d = (a, b) => Math.max(0.5, len(sub(a, b)));
    return 27.888 * (1 / d(o, n) + 1 / d(c, h) - 1 / d(o, h) - 1 / d(c, n));
  }
  function dsspCodes(ctx, bb) {
    const m = bb.count, xyz = ctx.xyz, codes = new Uint8Array(m);
    if (m < 3) return codes;
    const prev = new Int32Array(m).fill(-1), next = new Int32Array(m).fill(-1);
    for (let k = 1; k < m; k++) if (bb.seg[k] === bb.seg[k - 1]) { prev[k] = k - 1; next[k - 1] = k; }
    const N = [], C = [], O = [], H = [], CA = [];
    for (let k = 0; k < m; k++) {
      CA.push(P(xyz, bb.ca[k]));
      N.push(bb.n[k] >= 0 ? P(xyz, bb.n[k]) : null); C.push(bb.c[k] >= 0 ? P(xyz, bb.c[k]) : null); O.push(bb.o[k] >= 0 ? P(xyz, bb.o[k]) : null);
    }
    const cls = reps().classify(ctx);
    for (let k = 0; k < m; k++) {
      const p = prev[k];
      if (N[k] && p >= 0 && C[p] && O[p] && cls.resName[bb.res[k]] !== 'PRO') H.push(add(N[k], norm(sub(C[p], O[p])))); else H.push(null);
    }
    // neighbour pairs on CA (9 Å)
    const cell = 9, grid = new Map(), key = (x, y, z) => x + ',' + y + ',' + z;
    for (let k = 0; k < m; k++) {
      const g = CA[k].map(v => Math.floor(v / cell)), kk = key(g[0], g[1], g[2]);
      if (!grid.has(kk)) grid.set(kk, []);
      grid.get(kk).push(k);
    }
    const hb = new Set(), pairs = [];
    for (let i = 0; i < m; i++) {
      const g = CA[i].map(v => Math.floor(v / cell));
      for (let dx = -1; dx <= 1; dx++) for (let dy = -1; dy <= 1; dy++) for (let dz = -1; dz <= 1; dz++) {
        const list = grid.get(key(g[0] + dx, g[1] + dy, g[2] + dz));
        if (!list) continue;
        for (const j of list) {
          if (j <= i) continue;
          const d = len(sub(CA[i], CA[j]));
          if (d > 9) continue;
          pairs.push(i, j);
          // acceptor i (C=O) - donor j (N-H) and the reverse; neighbours in sequence never H-bond
          if (Math.abs(i - j) < 2 && bb.seg[i] === bb.seg[j]) continue;
          if (C[i] && O[i] && N[j] && H[j] && hbEnergy(C[i], O[i], N[j], H[j]) < -0.5) hb.add(i * m + j);
          if (C[j] && O[j] && N[i] && H[i] && hbEnergy(C[j], O[j], N[i], H[i]) < -0.5) hb.add(j * m + i);
        }
      }
    }
    const HB = (a, b) => a >= 0 && b >= 0 && hb.has(a * m + b);
    const ahead = (k, s) => { let x = k; for (let q = 0; q < s && x >= 0; q++) x = next[x]; return x; };
    // alpha helices: two consecutive 4-turns at i-1 and i -> residues i..i+3
    const turn4 = new Uint8Array(m), helix = new Uint8Array(m);
    for (let k = 0; k < m; k++) { const j = ahead(k, 4); if (j >= 0 && HB(k, j)) turn4[k] = 1; }
    for (let k = 1; k < m; k++) if (turn4[k] && prev[k] >= 0 && turn4[prev[k]]) { let x = k; for (let q = 0; q < 4 && x >= 0; q++) { helix[x] = 1; x = next[x]; } }
    // beta bridges (parallel and antiparallel), ladders of consecutive bridges -> strands
    const bridge = new Uint8Array(m);
    for (let q = 0; q < pairs.length; q += 2) {
      const i = pairs[q], j = pairs[q + 1];
      if (bb.seg[i] === bb.seg[j] && Math.abs(i - j) < 3) continue;
      const ip = prev[i], in_ = next[i], jp = prev[j], jn = next[j];
      const par = (HB(ip, j) && HB(j, in_)) || (HB(jp, i) && HB(i, jn));
      const anti = (HB(i, j) && HB(j, i)) || (HB(ip, jn) && HB(jp, in_));
      if (par || anti) { bridge[i] = 1; bridge[j] = 1; }
    }
    const strand = new Uint8Array(m);
    for (let k = 0; k < m; k++) if (bridge[k] && ((prev[k] >= 0 && bridge[prev[k]]) || (next[k] >= 0 && bridge[next[k]]))) strand[k] = 1;
    for (let k = 0; k < m; k++) if (!strand[k] && prev[k] >= 0 && next[k] >= 0 && strand[prev[k]] && strand[next[k]]) strand[k] = 3; // one-residue bulges
    for (let k = 0; k < m; k++) codes[k] = helix[k] ? 1 : strand[k] ? 2 : 0;
    return codes;
  }

  // ------------------------------------------------------------------ spline tubes and ribbons
  function catmull(p0, p1, p2, p3, u) {
    const u2 = u * u, u3 = u2 * u, out = [0, 0, 0], tan = [0, 0, 0];
    for (let k = 0; k < 3; k++) {
      const a = p0[k], b = p1[k], c = p2[k], d = p3[k];
      out[k] = 0.5 * (2 * b + (-a + c) * u + (2 * a - 5 * b + 4 * c - d) * u2 + (-a + 3 * b - 3 * c + d) * u3);
      tan[k] = 0.5 * ((-a + c) + 2 * (2 * a - 5 * b + 4 * c - d) * u + 3 * (-a + 3 * b - 3 * c + d) * u2);
    }
    return [out, tan];
  }
  /** Samples along a polyline of control points: {p, t, w, k (control index of the nearer point), u, k0}. */
  function sampleSpline(pts, sides, nsub) {
    const K = pts.length, out = [];
    const at = i => pts[Math.max(0, Math.min(K - 1, i))];
    const ext = i => (i < 0 ? sub(pts[0], sub(pts[1], pts[0])) : i >= K ? add(pts[K - 1], sub(pts[K - 1], pts[K - 2])) : pts[i]);
    for (let k = 0; k < K - 1; k++) {
      for (let s = 0; s < nsub; s++) {
        const u = s / nsub, [p, tg] = catmull(ext(k - 1), at(k), at(k + 1), ext(k + 2), u);
        out.push({ p, t: norm(tg), w: sides ? lerp(sides[k], sides[k + 1], u) : null, k: u < 0.5 ? k : k + 1, k0: k, u });
      }
    }
    const [p, tg] = catmull(ext(K - 2), at(K - 2), at(K - 1), ext(K), 1);
    out.push({ p, t: norm(tg), w: sides ? sides[K - 1] : null, k: K - 1, k0: K - 1, u: 0 });
    // rotation-minimising frame when no side vectors are given
    let prevW = null;
    for (const s of out) {
      let w = s.w ? sub(s.w, mul(s.t, dot(s.w, s.t))) : prevW ? sub(prevW, mul(s.t, dot(prevW, s.t))) : perp(s.t);
      if (len(w) < 1e-4) w = prevW ? sub(prevW, mul(s.t, dot(prevW, s.t))) : perp(s.t);
      if (len(w) < 1e-4) w = perp(s.t);
      w = norm(w);
      if (prevW && dot(w, prevW) < 0 && !s.w) w = mul(w, -1);
      s.w = w; prevW = w;
    }
    return out;
  }
  /** Elliptical tube through samples (a = half-width along w, b = half-thickness), caps at strip ends.
   *  samples[q].dims = [a, b], .col, .brk (start a new strip here), .skip (no geometry). */
  function sweep(m, samples, seg) {
    let ring = null, prevS = null;
    const cs = [], sn = [];
    for (let j = 0; j <= seg; j++) { const th = 2 * Math.PI * j / seg; cs.push(Math.cos(th)); sn.push(Math.sin(th)); }
    const makeRing = s => {
      const [a, b] = s.dims, W = s.w, B = norm(cross(s.t, W)), idx = [];
      for (let j = 0; j <= seg; j++) {
        const pos = add(s.p, add(mul(W, a * cs[j]), mul(B, b * sn[j])));
        const nrm = norm(add(mul(W, cs[j] / a), mul(B, sn[j] / b)));
        idx.push(m.v(pos, nrm, s.col));
      }
      return idx;
    };
    const cap = (s, dir) => {
      const [a, b] = s.dims, W = s.w, B = norm(cross(s.t, W)), nn = mul(s.t, dir), ctr = m.v(s.p, nn, s.col), idx = [];
      for (let j = 0; j <= seg; j++) idx.push(m.v(add(s.p, add(mul(W, a * cs[j]), mul(B, b * sn[j]))), nn, s.col));
      for (let j = 0; j < seg; j++) { if (dir > 0) m.t(ctr, idx[j], idx[j + 1]); else m.t(ctr, idx[j + 1], idx[j]); }
    };
    for (const s of samples) {
      if (s.skip) { if (ring) { cap(prevS, 1); ring = null; } prevS = null; continue; }
      if (s.brk && ring) { cap(prevS, 1); ring = null; }
      const r = makeRing(s);
      if (!ring) cap(s, -1);
      else for (let j = 0; j < seg; j++) { m.t(ring[j], r[j + 1], r[j]); m.t(ring[j], ring[j + 1], r[j + 1]); }
      ring = r; prevS = s;
    }
    if (ring && prevS) cap(prevS, 1);
  }
  const DIMS = { 0: [0.3, 0.3], 1: [1.25, 0.24], 2: [1.1, 0.26] };
  const smooth = t => t * t * (3 - 2 * t);

  function proteinCartoon(ctx, mask, colors, opts, m, proxies) {
    const bb = proteinChains(ctx), cls = reps().classify(ctx), xyz = ctx.xyz;
    if (!bb.count) return;
    const ssr = reps().secondary(ctx) || new Uint8Array(cls.nRes);
    const tot = bb.count, nsub = tot > 3000 ? 3 : tot > 1200 ? 4 : 7, seg = tot > 3000 ? 8 : tot > 1200 ? 10 : 14;
    const cyl = !!(opts && opts.cylinders), f = (opts && opts.scale) || 1;
    // runs of consecutive residues in the mask within one segment
    let k = 0;
    while (k < tot) {
      if (!mask[bb.ca[k]]) { k++; continue; }
      let e = k;
      while (e + 1 < tot && bb.seg[e + 1] === bb.seg[k] && mask[bb.ca[e + 1]]) e++;
      if (e > k) chainRun(k, e);
      else proxies.push(bb.ca[k]);
      k = e + 1;
    }
    function chainRun(s, e) {
      const K = e - s + 1, pts = [], sides = [], ss = [], cols = [];
      for (let q = 0; q < K; q++) {
        const j = s + q;
        pts.push(P(xyz, bb.ca[j])); ss.push(ssr[bb.res[j]]); cols.push(colAt(colors, bb.ca[j])); proxies.push(bb.ca[j]);
      }
      // flatten strands: average pleated CA positions
      const sm = pts.map((p, q) => (ss[q] === 2 && q > 0 && q < K - 1 && ss[q - 1] === 2 && ss[q + 1] === 2 ? mul(add(add(pts[q - 1], mul(p, 2)), pts[q + 1]), 0.25) : p));
      for (let q = 0; q < K; q++) {
        const j = s + q;
        const tg = norm(sub(sm[Math.min(K - 1, q + 1)], sm[Math.max(0, q - 1)]));
        let w = bb.c[j] >= 0 && bb.o[j] >= 0 ? sub(P(xyz, bb.o[j]), P(xyz, bb.c[j])) : (q ? sides[q - 1] : perp(tg));
        w = sub(w, mul(tg, dot(w, tg)));
        if (len(w) < 1e-3) w = q ? sides[q - 1] : perp(tg);
        w = norm(w);
        if (q && dot(w, sides[q - 1]) < 0) w = mul(w, -1);
        sides.push(w);
      }
      // helix runs drawn as cylinders (cylinder mode)
      const inCyl = new Uint8Array(K);
      if (cyl) {
        let q = 0;
        while (q < K) {
          if (ss[q] !== 1) { q++; continue; }
          let r = q;
          while (r + 1 < K && ss[r + 1] === 1) r++;
          if (r - q + 1 >= 4) { helixCylinder(m, pts.slice(q, r + 1), cols.slice(q, r + 1), 2.3 * f); for (let x = q + 1; x < r; x++) inCyl[x] = 1; inCyl[q] = 2; inCyl[r] = 3; }
          q = r + 1;
        }
      }
      const dimsOf = q => { const d0 = DIMS[cyl && ss[q] === 1 ? 0 : ss[q]]; return [d0[0] * f, d0[1] * f]; };
      const samples = sampleSpline(sm, sides, nsub);
      for (const S of samples) {
        const k0 = S.k0, k1 = Math.min(K - 1, k0 + 1), u = S.u;
        let d0 = dimsOf(k0), d1 = dimsOf(k1);
        // arrowhead over the span after a strand's last residue (or over its final span at a chain end)
        const arrow = ss[k0] === 2 && k1 !== k0 && (ss[k1] !== 2 || (k1 === K - 1 && k0 === K - 2));
        if (arrow) { d0 = [1.75 * f, 0.26 * f]; d1 = [DIMS[0][0] * f, DIMS[0][1] * f]; }
        const t = arrow ? u : smooth(u);
        S.dims = [d0[0] + (d1[0] - d0[0]) * t, d0[1] + (d1[1] - d0[1]) * t];
        S.col = cols[S.k];
        if (arrow && u === 0) S.brk = true;
        if (opts && opts.profile) opts.profile.push({ ss: ss[S.k], arrow, a: S.dims[0], b: S.dims[1], u });
        if (cyl && (inCyl[k0] === 1 || inCyl[k0] === 2) && (inCyl[k1] === 1 || inCyl[k1] === 3) && !(inCyl[k0] === 2 && u === 0)) S.skip = true;
      }
      // the arrow base starts a new strip whose first ring is the wide one: the old strip first runs on to the base
      const out = [];
      for (const S of samples) { if (S.brk) out.push(Object.assign({}, S, { brk: false, dims: dimsOf(S.k0) })); out.push(S); }
      sweep(m, out, seg);
    }
  }
  function helixCylinder(m, pts, cols, radius) {
    const K = pts.length, A = [];
    for (let q = 0; q + 3 < K; q++) A.push(mul(add(add(pts[q], pts[q + 1]), add(pts[q + 2], pts[q + 3])), 0.25));
    if (A.length < 1) return;
    const c = mul(A.reduce((s, p) => add(s, p), [0, 0, 0]), 1 / A.length);
    let d = A.length > 1 ? norm(sub(A[A.length - 1], A[0])) : norm(sub(pts[K - 1], pts[0]));
    // refine with the principal axis of the axis points
    const cov = [[0, 0, 0], [0, 0, 0], [0, 0, 0]];
    for (const p of A) { const v = sub(p, c); for (let a = 0; a < 3; a++) for (let b = 0; b < 3; b++) cov[a][b] += v[a] * v[b]; }
    for (let it = 0; it < 30 && A.length > 2; it++) { const v = [dot(cov[0], d), dot(cov[1], d), dot(cov[2], d)]; if (len(v) < 1e-9) break; d = norm(v); }
    if (dot(d, sub(pts[K - 1], pts[0])) < 0) d = mul(d, -1);
    const t0 = dot(sub(pts[0], c), d), t1 = dot(sub(pts[K - 1], c), d);
    const a = add(c, mul(d, t0)), b = add(c, mul(d, t1));
    for (let q = 0; q < K; q++) cylinder(m, lerp(a, b, q / K), lerp(a, b, (q + 1) / K), radius || 2.3, cols[q], 24, q === 0 || q === K - 1);
  }

  function nucleicCartoon(ctx, mask, colors, opts, m, slabs, proxies) {
    const nb = nucleicChains(ctx), xyz = ctx.xyz;
    if (!nb.count) return;
    const tot = nb.count, nsub = tot > 2000 ? 3 : 6, seg = tot > 2000 ? 8 : 12;
    let k = 0;
    while (k < tot) {
      if (!mask[nb.ctl[k]]) { k++; continue; }
      let e = k;
      while (e + 1 < tot && nb.seg[e + 1] === nb.seg[k] && mask[nb.ctl[e + 1]]) e++;
      const pts = [], cols = [];
      for (let q = k; q <= e; q++) { pts.push(P(xyz, nb.ctl[q])); cols.push(colAt(colors, nb.ctl[q])); proxies.push(nb.ctl[q]); }
      let samples = null;
      if (pts.length > 1) {
        samples = sampleSpline(pts, null, nsub);
        const f = (opts && opts.scale) || 1;
        for (const S of samples) { S.dims = [0.9 * f, 0.9 * f]; S.col = cols[S.k]; }
        sweep(m, samples, seg);
      }
      for (let q = k; q <= e; q++) {
        const b = nb.bases[q];
        if (!b) continue;
        const col = colAt(colors, nb.ctl[q]);
        baseSlab(slabs, xyz, b, col);
        // rod from the backbone (spline point between this phosphate and the next) to the glycosidic N
        const j = q - k, sp = samples ? samples[Math.min(samples.length - 1, j * nsub + Math.floor(nsub / 2))].p : pts[j];
        cylinder(m, sp, P(xyz, b.attach), 0.3, col, 8, false);
      }
      k = e + 1;
    }
  }
  function baseSlab(m, xyz, b, col) {
    const per = b.perim.map(i => P(xyz, i)), c = mul(per.reduce((s, p) => add(s, p), [0, 0, 0]), 1 / per.length);
    let nrm = [0, 0, 0];
    for (let q = 0; q < per.length; q++) nrm = add(nrm, cross(per[q], per[(q + 1) % per.length]));
    nrm = norm(nrm);
    const h = 0.22, up = mul(nrm, h), dn = mul(nrm, -h);
    const fan = (ring, off, n) => {
      const pts = ring.map(i => add(P(xyz, i), off)), cc = mul(pts.reduce((s, p) => add(s, p), [0, 0, 0]), 1 / pts.length), ci = m.v(cc, n, col), ids = pts.map(p => m.v(p, n, col));
      for (let q = 0; q < ids.length; q++) m.t(ci, ids[q], ids[(q + 1) % ids.length]);
    };
    for (const ring of [b.r6].concat(b.r5 ? [b.r5] : [])) { fan(ring, up, nrm); fan(ring, dn, mul(nrm, -1)); }
    for (let q = 0; q < per.length; q++) {
      const a = per[q], bq = per[(q + 1) % per.length];
      m.flatQuad(add(a, dn), add(bq, dn), add(bq, up), add(a, up), col);
    }
    void c;
  }

  function trace(ctx, mask, colors, m, proxies, f) {
    f = f || 1;
    const xyz = ctx.xyz;
    for (const bb of [proteinChains(ctx), nucleicChains(ctx)]) {
      const ctl = bb.ca || bb.ctl, tot = bb.count, nsub = tot > 3000 ? 2 : 4, seg = tot > 3000 ? 6 : 10;
      let k = 0;
      while (k < tot) {
        if (!mask[ctl[k]]) { k++; continue; }
        let e = k;
        while (e + 1 < tot && bb.seg[e + 1] === bb.seg[k] && mask[ctl[e + 1]]) e++;
        const pts = [], cols = [];
        for (let q = k; q <= e; q++) { pts.push(P(xyz, ctl[q])); cols.push(colAt(colors, ctl[q])); proxies.push(ctl[q]); }
        if (pts.length > 1) { const S = sampleSpline(pts, null, nsub); for (const s of S) { s.dims = [0.38 * f, 0.38 * f]; s.col = cols[s.k]; } sweep(m, S, seg); }
        k = e + 1;
      }
    }
  }

  // ------------------------------------------------------------------ rings and SNFG symbols
  function ringFrame(xyz, ring, first) {
    const pts = ring.map(i => P(xyz, i)), c = mul(pts.reduce((s, p) => add(s, p), [0, 0, 0]), 1 / pts.length);
    let n = [0, 0, 0];
    for (let q = 0; q < pts.length; q++) n = add(n, cross(sub(pts[q], c), sub(pts[(q + 1) % pts.length], c)));
    n = norm(n);
    let u = sub(P(xyz, first === undefined ? ring[0] : first), c);
    u = norm(sub(u, mul(n, dot(u, n))));
    return { c, n, u, v: cross(n, u), pts };
  }
  function perceiveRings(ctx, mask, maxRing) {
    const A = reps().adjacency(ctx), n = ctx.n, seen = new Set(), out = [];
    let cnt = 0;
    for (let i = 0; i < n; i++) if (mask[i]) cnt++;
    if (cnt > 60000) return out;
    for (let i = 0; i < n; i++) {
      if (!mask[i] || ctx.Z[i] <= 1) continue;
      for (let q = A.off[i]; q < A.off[i + 1]; q++) {
        const j = A.nbr[q];
        if (j <= i || !mask[j] || ctx.Z[j] <= 1) continue;
        // shortest path j -> i avoiding the direct bond (heavy masked atoms only)
        const prev = new Map([[i, -1]]);
        let front = [i], path = null;
        for (let depth = 1; depth < maxRing && front.length && !path; depth++) {
          const nx = [];
          for (const u of front) {
            for (let r = A.off[u]; r < A.off[u + 1] && !path; r++) {
              const v = A.nbr[r];
              if ((u === i && v === j) || prev.has(v) || !mask[v] || ctx.Z[v] <= 1) continue;
              prev.set(v, u);
              if (v === j) { path = []; for (let x = v; x !== -1; x = prev.get(x)) path.push(x); }
              else nx.push(v);
            }
          }
          front = nx;
        }
        if (!path || path.length < 3) continue;
        const key = path.slice().sort((a, b) => a - b).join(',');
        if (seen.has(key)) continue;
        seen.add(key); out.push(path);
      }
    }
    return out;
  }
  function rings(ctx, mask, colors, opts, m) {
    const xyz = ctx.xyz, scheme = opts && opts.scheme, done = new Set();
    for (const s of reps().sugars(ctx)) {
      if (!s.ring.some(a => mask[a])) continue;
      const col = scheme === 'sugar' || !scheme ? s.color : avgColor(colors, s.ring);
      fillRing(m, xyz, s.ring, col);
      done.add(s.ring.slice().sort((a, b) => a - b).join(','));
    }
    for (const r of perceiveRings(ctx, mask, 7)) {
      const key = r.slice().sort((a, b) => a - b).join(',');
      if (done.has(key)) continue;
      const c = avgColor(colors, r);
      fillRing(m, xyz, r, [0.35 + 0.65 * c[0], 0.35 + 0.65 * c[1], 0.35 + 0.65 * c[2]].map(v => Math.min(1, v * 0.92)));
    }
  }
  function avgColor(colors, ring) { const c = [0, 0, 0]; for (const i of ring) for (let k = 0; k < 3; k++) c[k] += colors[3 * i + k] / ring.length; return c; }
  function fillRing(m, xyz, ring, col) {
    const f = ringFrame(xyz, ring), ci = m.v(f.c, f.n, col), ids = f.pts.map(p => m.v(p, f.n, col));
    for (let q = 0; q < ids.length; q++) m.t(ci, ids[q], ids[(q + 1) % ids.length]);
  }

  let ICO = null;
  function unitSphere() {
    if (ICO) return ICO;
    const t = (1 + Math.sqrt(5)) / 2;
    let v = [[-1, t, 0], [1, t, 0], [-1, -t, 0], [1, -t, 0], [0, -1, t], [0, 1, t], [0, -1, -t], [0, 1, -t], [t, 0, -1], [t, 0, 1], [-t, 0, -1], [-t, 0, 1]].map(norm);
    let f = [[0, 11, 5], [0, 5, 1], [0, 1, 7], [0, 7, 10], [0, 10, 11], [1, 5, 9], [5, 11, 4], [11, 10, 2], [10, 7, 6], [7, 1, 8], [3, 9, 4], [3, 4, 2], [3, 2, 6], [3, 6, 8], [3, 8, 9], [4, 9, 5], [2, 4, 11], [6, 2, 10], [8, 6, 7], [9, 8, 1]];
    for (let s = 0; s < 2; s++) {
      const cache = new Map(), nf = [];
      const mid = (a, b) => { const k = a < b ? a * 1e5 + b : b * 1e5 + a; if (cache.has(k)) return cache.get(k); v.push(norm(lerp(v[a], v[b], 0.5))); cache.set(k, v.length - 1); return v.length - 1; };
      for (const [a, b, c] of f) { const ab = mid(a, b), bc = mid(b, c), ca = mid(c, a); nf.push([a, ab, ca], [b, bc, ab], [c, ca, bc], [ab, bc, ca]); }
      f = nf;
    }
    ICO = { v, f };
    return ICO;
  }
  /** One SNFG symbol: local axes u, v in the ring plane, n along its normal; R = size. */
  function symbol(m, shape, F, R, col, white, inv) {
    const W = (x, y, z) => add(F.c, add(add(mul(F.u, x), mul(F.v, y)), mul(F.n, z)));
    const quad = (a, b, c, d, cl) => m.flatQuad(W(...a), W(...b), W(...c), W(...d), cl);
    const tri = (a, b, c, cl) => m.flatTri(W(...a), W(...b), W(...c), cl);
    switch (shape) {
      case 'sphere': {
        const s = unitSphere(), base = m.nv;
        for (const p of s.v) m.v(W(p[0] * R, p[1] * R, p[2] * R), norm(add(add(mul(F.u, p[0]), mul(F.v, p[1])), mul(F.n, p[2]))), col);
        for (const [a, b, c] of s.f) m.t(base + a, base + b, base + c);
        break;
      }
      case 'cube': case 'dcube': {
        const a = R * 0.82, h = R * 0.82;
        const box = (pts, cl) => { // prism over a polygon in the u-v plane (CCW), from z = -h to +h
          const k = pts.length;
          for (let q = 1; q + 1 < k; q++) { tri([pts[0][0], pts[0][1], h], [pts[q][0], pts[q][1], h], [pts[q + 1][0], pts[q + 1][1], h], cl); tri([pts[0][0], pts[0][1], -h], [pts[q + 1][0], pts[q + 1][1], -h], [pts[q][0], pts[q][1], -h], cl); }
          for (let q = 0; q < k; q++) { const p = pts[q], r = pts[(q + 1) % k]; quad([p[0], p[1], -h], [r[0], r[1], -h], [r[0], r[1], h], [p[0], p[1], h], cl); }
        };
        if (shape === 'cube') box([[-a, -a], [a, -a], [a, a], [-a, a]], col);
        else { box([[-a, -a], [a, -a], [a, a]], col); box([[-a, -a], [a, a], [-a, a]], white); }
        break;
      }
      case 'diamond': case 'ddiamond': {
        const a = R * 1.08, h = R * 0.62, X = [[a, 0], [0, a], [-a, 0], [0, -a]];
        for (let q = 0; q < 4; q++) {
          const p = X[q], r = X[(q + 1) % 4];
          const upper = p[1] > 0 || r[1] > 0, cl = shape === 'diamond' ? col : (upper !== !!inv ? col : white);
          tri([p[0], p[1], 0], [r[0], r[1], 0], [0, 0, h], cl); tri([r[0], r[1], 0], [p[0], p[1], 0], [0, 0, -h], cl);
        }
        break;
      }
      case 'cone': {
        const seg = 18, apex = [0, R * 1.05, 0], base = -R * 0.8, rb = R * 0.85;
        for (let q = 0; q < seg; q++) {
          const a0 = 2 * Math.PI * q / seg, a1 = 2 * Math.PI * (q + 1) / seg;
          const p0 = [rb * Math.cos(a0), base, rb * Math.sin(a0)], p1 = [rb * Math.cos(a1), base, rb * Math.sin(a1)];
          tri(apex, p1, p0, col); tri([0, base, 0], p0, p1, col);
        }
        break;
      }
      case 'star': case 'pentagon': case 'hexagon': {
        const h = R * 0.36, pts = [];
        const k = shape === 'star' ? 10 : shape === 'pentagon' ? 5 : 6;
        for (let q = 0; q < k; q++) { const ang = Math.PI / 2 + 2 * Math.PI * q / k, r = shape === 'star' && q % 2 ? R * 0.46 : R * 1.12; pts.push([r * Math.cos(ang), r * Math.sin(ang)]); }
        for (let q = 0; q < k; q++) {
          const p = pts[q], r = pts[(q + 1) % k];
          tri([0, 0, h], [p[0], p[1], h], [r[0], r[1], h], col); tri([0, 0, -h], [r[0], r[1], -h], [p[0], p[1], -h], col);
          quad([p[0], p[1], -h], [r[0], r[1], -h], [r[0], r[1], h], [p[0], p[1], h], col);
        }
        break;
      }
    }
  }
  function snfg(ctx, mask, colors, opts, m, proxies) {
    const xyz = ctx.xyz, cls = reps().classify(ctx), list = reps().sugars(ctx).filter(s => s.ring.some(a => mask[a]));
    const R = 1.25 * ((opts && opts.scale) || 1), white = [0.97, 0.97, 0.97], cen = new Map();
    for (const s of list) {
      const F = ringFrame(xyz, s.ring, s.anomeric);
      symbol(m, s.shape, F, R, s.color, white, s.inv);
      if (!cen.has(s.res)) cen.set(s.res, { c: F.c, ring: s.ring });
      proxies.push(s.anomeric);
    }
    // link sticks: sugar-sugar bonds between residues, and sugar-aglycone/protein attachments
    const A = reps().adjacency(ctx), linkCol = [0.80, 0.82, 0.86], seen = new Set();
    for (const [r, info] of cen) {
      for (let t = cls.resOff[r]; t < cls.resOff[r + 1]; t++) {
        const a = cls.resAtoms[t];
        for (let q = A.off[a]; q < A.off[a + 1]; q++) {
          const b = A.nbr[q], rb = cls.resOf[b];
          if (rb === r) continue;
          if (cls.resClass[rb] === C_SUGAR && cen.has(rb)) {
            const key = r < rb ? r + ':' + rb : rb + ':' + r;
            if (seen.has(key)) continue;
            seen.add(key);
            cylinder(m, info.c, cen.get(rb).c, 0.2, linkCol, 10, false);
          } else if (cls.resClass[rb] !== C_SUGAR && ctx.Z[b] > 1) {
            const key = r + '>' + b;
            if (seen.has(key)) continue;
            seen.add(key);
            cylinder(m, info.c, P(xyz, b), 0.2, linkCol, 10, false);
          }
        }
      }
    }
    return R;
  }

  /** Geometry for one overlay row. */
  function build(kind, ctx, mask, colors, opts) {
    const m = new Mesh(false), two = new Mesh(true), proxies = [];
    let radius = 1.4;
    switch (kind) {
      case 'cartoon': proteinCartoon(ctx, mask, colors, opts, m, proxies); nucleicCartoon(ctx, mask, colors, opts, m, two, proxies); break;
      case 'trace': trace(ctx, mask, colors, m, proxies, opts && opts.scale); radius = 0.9; break;
      case 'rings': rings(ctx, mask, colors, opts, two); break;
      case 'snfg': radius = snfg(ctx, mask, colors, opts, m, proxies); break;
    }
    return { meshes: [m.done(), two.done()].filter(Boolean), proxies, proxyRadius: radius };
  }

  BM.cartoon = { secondary, dssp, hbEnergy, proteinChains, nucleicChains, build, perceiveRings, Mesh };
})(typeof window !== 'undefined' ? (window.BM = window.BM || {}) : (globalThis.BM = globalThis.BM || {}));
