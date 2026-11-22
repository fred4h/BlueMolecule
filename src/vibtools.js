/* BlueMolecule — normal-mode tools for transition states (issue #14), BM.vibtools:
 *   stationary(vib)              character of the structure from its frequencies (minimum, transition state, ...)
 *   describe(vib, st, tsJob)     that character in words: {short, text, ok}
 *   ends(eq, mode, amp)          the two ends of a mode: eq -/+ the displacement scaled to a largest atom shift of amp Å
 *   bondChanges(Z, eq, u)        atom pairs whose distance changes most along a displacement u: the bonds that form and
 *                                break along a transition-state mode, with the distances at the - and + ends
 *   endMeshes(Z, eq, u)          translucent ball-and-stick meshes of both ends (renderer overlay meshes)
 *   minimumRoute(route)          a Gaussian route for optimising a geometry displaced from a transition state
 *   minimumKeywords(keywords)    the same for an ORCA keyword line
 * Frequencies are in cm-1, imaginary ones negative, as in every dataset. */
(function (BM) {
  'use strict';
  const SMALL = 50;   // cm-1: imaginary frequencies below this are often numerical noise

  const isNum = v => typeof v === 'number' && Number.isFinite(v);
  const imagText = f => Math.abs(f).toFixed(1) + 'i';

  /** {kind, imag, large, small, n, threshold}: kind is 'minimum' (no imaginary frequency), 'ts' (exactly one
   *  imaginary frequency of at least `small` cm-1), 'saddle' (more than one) or 'noise' (only small imaginary
   *  frequencies). imag lists the imaginary modes, most negative first; large and small split it at `small`.
   *  Frequencies between -1 and 0 cm-1 (translations and rotations printed as -0.0) do not count. */
  function stationary(vib, opts) {
    const freqs = (vib && vib.freqs) || [];
    const small = opts && opts.small > 0 ? opts.small : SMALL;
    const imag = [];
    for (let i = 0; i < freqs.length; i++) if (isNum(freqs[i]) && freqs[i] < -1) imag.push(i);
    imag.sort((a, b) => freqs[a] - freqs[b] || a - b);
    const large = imag.filter(i => -freqs[i] >= small), tiny = imag.filter(i => -freqs[i] < small);
    const kind = !imag.length ? 'minimum' : large.length === 1 ? 'ts' : large.length > 1 ? 'saddle' : 'noise';
    return { kind, imag, large, small: tiny, n: freqs.length, threshold: small };
  }

  /** Frequencies as "1532.4i, 45.2i cm⁻¹" (at most `max`, then "…"). */
  function listImag(freqs, idx, max) {
    max = max || 4;
    return idx.slice(0, max).map(i => imagText(freqs[i])).join(', ') + (idx.length > max ? ', …' : '') + ' cm⁻¹';
  }

  /** The stationary point in words. short: for a status line ("one imaginary frequency (1532.4i cm⁻¹)");
   *  text: one or two sentences for the Vibrations tab; ok: false when the result needs a second look (imaginary
   *  frequencies in a minimum, none or several in a transition-state job). */
  function describe(vib, st, tsJob) {
    st = st || stationary(vib);
    const f = (vib && vib.freqs) || [], nI = st.imag.length, nS = st.small.length, thr = st.threshold;
    const count = n => (n === 1 ? 'one imaginary frequency' : n + ' imaginary frequencies');
    const short = nI ? count(nI) + ' (' + listImag(f, st.imag) + ')' : 'no imaginary frequencies';
    const cap = s => s.charAt(0).toUpperCase() + s.slice(1);
    const why = '(loose convergence or a coarse integration grid)';
    let text, ok;
    if (st.kind === 'minimum') {
      text = tsJob ? 'No imaginary frequencies: this is a minimum, not a transition state.' : 'No imaginary frequencies: a minimum.';
      ok = !tsJob;
    } else if (st.kind === 'ts') {
      text = nS ? cap(short) + ': one transition-state mode and ' + (nS === 1 ? 'one small one' : nS + ' small ones') + ' below ' + thr + 'i cm⁻¹, often numerical noise ' + why + '.'
        : cap(short) + ': a transition state.';
      ok = true;
    } else if (st.kind === 'saddle') {
      text = cap(short) + ': a higher-order saddle point, neither a minimum nor a transition state.';
      ok = false;
    } else {
      text = cap(short) + (nI === 1 ? ', below ' : ', all below ') + thr + 'i cm⁻¹: often numerical noise in a minimum ' + why +
        (tsJob ? '; a transition state needs one clearly imaginary frequency.' : (nI === 1 ? ', unless this is a very flat transition state.' : '.'));
      ok = false;
    }
    return { short, text, ok };
  }

  /** Displacement u of `mode` scaled so the largest atom shift is amp Å, and the geometries eq - u and eq + u. */
  function ends(eq, mode, amp) {
    const u = BM.geom.normalizedMode(mode, amp);
    const minus = new Float64Array(eq.length), plus = new Float64Array(eq.length);
    for (let k = 0; k < eq.length; k++) { minus[k] = eq[k] - u[k]; plus[k] = eq[k] + u[k]; }
    return { u, minus, plus };
  }

  const dist = (x, i, j) => Math.hypot(x[3 * j] - x[3 * i], x[3 * j + 1] - x[3 * i + 1], x[3 * j + 2] - x[3 * i + 2]);
  const BOND_TOL = 0.4, NEAR = 1.7, PARTIAL = 1.15;

  /** Atom pairs whose distance changes most between eq - u and eq + u, largest change first:
   *  [{i, j, d0, dMinus, dPlus, delta (dPlus - dMinus), ratio (d0 / sum of covalent radii), bonded0, bondedMinus,
   *  bondedPlus, change}]. Pairs: bonded (covalent radii + 0.4 Å, as drawn) at one of the three geometries, or
   *  closer than 1.7 × the sum of their covalent radii here, with at least one atom that moves; pairs that only
   *  span an angle (both bonded to the same atom, not to each other) are left out. change says what happens
   *  towards the + end: 'forms' or 'breaks' for a partial bond along an imaginary mode (opts.imaginary: longer than
   *  1.15 × the covalent radii here, or bonded at one end only), otherwise 'shorter' or 'longer'.
   *  opts: max rows (8), fraction of the largest change a row needs (0.25), imaginary. */
  function bondChanges(Z, eq, u, opts) {
    opts = opts || {};
    const n = Z.length, cov = BM.el.covalent, out = [];
    if (!n || !u || u.length !== 3 * n) return out;
    const um = new Float64Array(n);
    let umax = 0, maxR = 0;
    for (let i = 0; i < n; i++) {
      um[i] = Math.hypot(u[3 * i], u[3 * i + 1], u[3 * i + 2]);
      if (um[i] > umax) umax = um[i];
      if (Z[i] > 0 && cov[Z[i]] > maxR) maxR = cov[Z[i]];
    }
    if (!(umax > 0) || !(maxR > 0)) return out;
    const minus = new Float64Array(3 * n), plus = new Float64Array(3 * n);
    for (let k = 0; k < 3 * n; k++) { minus[k] = eq[k] - u[k]; plus[k] = eq[k] + u[k]; }
    // neighbours within 1.7 x the covalent radii (or bonded at an end) of the atoms that move, on a grid
    const cut = Math.max(2 * maxR * NEAR, 2 * maxR + BOND_TOL) + 2 * umax;
    const grid = new Map(), key = (a, b, c) => a + ',' + b + ',' + c, cell = i => [0, 1, 2].map(k => Math.floor(eq[3 * i + k] / cut));
    for (let i = 0; i < n; i++) {
      if (!(Z[i] > 0)) continue;
      const k = key(...cell(i));
      let l = grid.get(k);
      if (!l) grid.set(k, l = []);
      l.push(i);
    }
    const bonds = BM.geom.computeBonds(Z, eq), adj = BM.geom.adjacency(n, bonds);
    const moving = i => um[i] > 0.02 * umax, seen = new Set();
    for (let i = 0; i < n; i++) {
      if (!(Z[i] > 0) || !moving(i)) continue;
      const [cx, cy, cz] = cell(i);
      for (let dx = -1; dx <= 1; dx++) for (let dy = -1; dy <= 1; dy++) for (let dz = -1; dz <= 1; dz++) {
        for (const j of grid.get(key(cx + dx, cy + dy, cz + dz)) || []) {
          if (j === i || (moving(j) && j < i)) continue;
          const a = Math.min(i, j), b = Math.max(i, j), id = a * n + b;
          if (seen.has(id)) continue;
          seen.add(id);
          const rc = cov[Z[a]] + cov[Z[b]], d0 = dist(eq, a, b), dm = dist(minus, a, b), dp = dist(plus, a, b);
          const bond = d => d < rc + BOND_TOL && d > 0.4;
          const b0 = bond(d0), bm = bond(dm), bp = bond(dp);
          if (!b0 && !bm && !bp) {
            if (d0 >= NEAR * rc) continue;
            if (adj[a].some(c => adj[b].indexOf(c) >= 0)) continue;   // an angle, not a bond
          }
          out.push({ i: a, j: b, d0, dMinus: dm, dPlus: dp, delta: dp - dm, ratio: d0 / rc, bonded0: b0, bondedMinus: bm, bondedPlus: bp });
        }
      }
    }
    out.sort((p, q) => Math.abs(q.delta) - Math.abs(p.delta) || p.i - q.i || p.j - q.j);
    const top = out.length ? Math.abs(out[0].delta) : 0, frac = isNum(opts.fraction) ? opts.fraction : 0.25;
    const rows = out.filter(r => Math.abs(r.delta) >= Math.max(frac * top, 1e-3)).slice(0, opts.max || 8);
    for (const r of rows) {
      const partial = opts.imaginary && (r.ratio > PARTIAL || r.bondedMinus !== r.bondedPlus);
      r.change = partial ? (r.delta < 0 ? 'forms' : 'breaks') : (r.delta < 0 ? 'shorter' : 'longer');
    }
    return rows;
  }

  // ---------------------------------------------------------------- meshes of the two ends
  let SPHERE = null;
  /** Icosphere with one subdivision (42 vertices): ghosts of many atoms stay light. */
  function unitSphere() {
    if (SPHERE) return SPHERE;
    const nrm = p => { const l = Math.hypot(p[0], p[1], p[2]); return [p[0] / l, p[1] / l, p[2] / l]; };
    const t = (1 + Math.sqrt(5)) / 2;
    const v = [[-1, t, 0], [1, t, 0], [-1, -t, 0], [1, -t, 0], [0, -1, t], [0, 1, t], [0, -1, -t], [0, 1, -t], [t, 0, -1], [t, 0, 1], [-t, 0, -1], [-t, 0, 1]].map(nrm);
    const f0 = [[0, 11, 5], [0, 5, 1], [0, 1, 7], [0, 7, 10], [0, 10, 11], [1, 5, 9], [5, 11, 4], [11, 10, 2], [10, 7, 6], [7, 1, 8], [3, 9, 4], [3, 4, 2], [3, 2, 6], [3, 6, 8], [3, 8, 9], [4, 9, 5], [2, 4, 11], [6, 2, 10], [8, 6, 7], [9, 8, 1]];
    const cache = new Map(), f = [];
    const mid = (a, b) => {
      const k = a < b ? a * 64 + b : b * 64 + a;
      if (!cache.has(k)) { v.push(nrm([(v[a][0] + v[b][0]) / 2, (v[a][1] + v[b][1]) / 2, (v[a][2] + v[b][2]) / 2])); cache.set(k, v.length - 1); }
      return cache.get(k);
    };
    for (const [a, b, c] of f0) { const ab = mid(a, b), bc = mid(b, c), ca = mid(c, a); f.push([a, ab, ca], [b, bc, ab], [c, ca, bc], [ab, bc, ca]); }
    SPHERE = { v, f };
    return SPHERE;
  }
  function sphere(m, p, r, col) {
    const s = unitSphere(), base = m.nv;
    for (const q of s.v) m.v([p[0] + q[0] * r, p[1] + q[1] * r, p[2] + q[2] * r], q, col);
    for (const [a, b, c] of s.f) m.t(base + a, base + b, base + c);
  }
  function cylinder(m, a, b, r, col, seg) {
    const ax = [b[0] - a[0], b[1] - a[1], b[2] - a[2]], L = Math.hypot(ax[0], ax[1], ax[2]);
    if (!(L > 1e-6)) return;
    const T = [ax[0] / L, ax[1] / L, ax[2] / L];
    const o = Math.abs(T[0]) < 0.9 ? [1, 0, 0] : [0, 1, 0];
    let U = [T[1] * o[2] - T[2] * o[1], T[2] * o[0] - T[0] * o[2], T[0] * o[1] - T[1] * o[0]];
    const lu = Math.hypot(U[0], U[1], U[2]); U = [U[0] / lu, U[1] / lu, U[2] / lu];
    const V = [T[1] * U[2] - T[2] * U[1], T[2] * U[0] - T[0] * U[2], T[0] * U[1] - T[1] * U[0]], ring = [];
    for (let k = 0; k <= seg; k++) {
      const th = 2 * Math.PI * k / seg, c = Math.cos(th), s = Math.sin(th), d = [U[0] * c + V[0] * s, U[1] * c + V[1] * s, U[2] * c + V[2] * s];
      ring.push([m.v([a[0] + d[0] * r, a[1] + d[1] * r, a[2] + d[2] * r], d, col), m.v([b[0] + d[0] * r, b[1] + d[1] * r, b[2] + d[2] * r], d, col)]);
    }
    for (let k = 0; k < seg; k++) { const [a0, b0] = ring[k], [a1, b1] = ring[k + 1]; m.t(a0, a1, b1); m.t(a0, b1, b0); }
  }
  const MINUS_COLOR = [0.5, 0.72, 1.0], PLUS_COLOR = [1.0, 0.7, 0.38];

  /** Translucent ball-and-stick meshes of eq - u (blue) and eq + u (orange): [minusMesh, plusMesh] for
   *  renderer.setVibEnds. Every atom of small structures; for larger ones the atoms that move (at least 2 % of the
   *  largest shift), at most opts.maxAtoms (1500) of them. opts: radius(i) for the spheres (default 0.85 × the
   *  ball-and-stick radius; 0 hides the atom and its bonds), bondRadius (0.065 Å), opacity (0.5), minusColor,
   *  plusColor. */
  function endMeshes(Z, eq, u, opts) {
    opts = opts || {};
    const n = Z.length, Mesh = BM.cartoon && BM.cartoon.Mesh;
    if (!Mesh || !n || !u || u.length !== 3 * n) return [];
    const um = new Float64Array(n);
    let umax = 0;
    for (let i = 0; i < n; i++) { um[i] = Math.hypot(u[3 * i], u[3 * i + 1], u[3 * i + 2]); if (um[i] > umax) umax = um[i]; }
    if (!(umax > 0)) return [];
    let keep = [];
    for (let i = 0; i < n; i++) if (Z[i] > 0 && (n <= 400 || um[i] >= 0.02 * umax)) keep.push(i);
    const cap = opts.maxAtoms || 1500;
    if (keep.length > cap) keep = keep.sort((a, b) => um[b] - um[a] || a - b).slice(0, cap).sort((a, b) => a - b);
    const inKeep = new Uint8Array(n);
    for (const i of keep) inKeep[i] = 1;
    const vdw = BM.el.vdw, radius = typeof opts.radius === 'function' ? opts.radius : i => 0.85 * 0.25 * vdw[Z[i]];
    const br = opts.bondRadius > 0 ? opts.bondRadius : 0.065, out = [];
    for (const [sign, col] of [[-1, opts.minusColor || MINUS_COLOR], [1, opts.plusColor || PLUS_COLOR]]) {
      const x = new Float64Array(3 * n);
      for (let k = 0; k < 3 * n; k++) x[k] = eq[k] + sign * u[k];
      const m = new Mesh(false), P = i => [x[3 * i], x[3 * i + 1], x[3 * i + 2]];
      for (const i of keep) { const r = radius(i); if (r > 0) sphere(m, P(i), r, col); else inKeep[i] = 0; }   // radius 0: hidden, with its bonds
      const b = BM.geom.computeBonds(Z, x);
      for (let k = 0; k < b.length; k += 2) if (inKeep[b[k]] && inKeep[b[k + 1]]) cylinder(m, P(b[k]), P(b[k + 1]), br, col, 8);
      const mesh = m.done();
      if (mesh) { mesh.opacity = opts.opacity > 0 ? opts.opacity : 0.5; mesh.sign = sign; out.push(mesh); }
    }
    return out;
  }

  // ---------------------------------------------------------------- follow-up optimisations
  /** Split a Gaussian route into keywords (spaces or commas outside parentheses separate them). */
  function routeTokens(route) {
    const s = String(route || '').replace(/\s+/g, ' ').replace(/\s*=\s*/g, '=').replace(/\s*\(\s*/g, '(').replace(/\s*\)/g, ')').replace(/\s*,\s*/g, ',').trim();
    const out = [];
    let cur = '', depth = 0;
    for (const ch of s) {
      if (ch === '(') depth++;
      else if (ch === ')') depth = Math.max(0, depth - 1);
      if ((ch === ' ' || ch === ',') && depth === 0) { if (cur) out.push(cur); cur = ''; continue; }
      cur += ch;
    }
    if (cur) out.push(cur);
    return out;
  }
  // options of opt= that only make sense for a transition-state search, or need input we do not write
  const TS_OPT = /^(ts|calcfc|calcall|rcfc|readfc|noeigentest|eigentest|qst2|qst3|saddle=?\d*|modredundant|addredundant|addgic|readallgic|readgic|readopt|readoptimize|readfreeze)$/i;
  /** Gaussian route for optimising a geometry displaced from a transition state: the TS options of opt= removed
   *  (ts, calcfc, noeigentest, qst2, ...), opt added when the route has none, IRC and scan keywords dropped, and
   *  geom= and guess= dropped (the geometry is in the new input, the checkpoint file is new). */
  function minimumRoute(route) {
    const toks = [], out = [];
    for (const t of routeTokens(route)) {   // "#opt=ts": the # and a keyword written together
      if (t[0] === '#' && !/^#[pnt]?$/i.test(t)) toks.push('#', t.slice(1)); else toks.push(t);
    }
    let hasOpt = false;
    for (const t of toks) {
      const m = /^([a-z]+)(?:=?\((.*)\)|=(.*))?$/i.exec(t), name = m ? m[1].toLowerCase() : '';
      if (/^(geom|geometry|guess|irc|ircmax|scan)$/.test(name)) continue;
      if (name === 'opt' || name === 'optimize') {
        const opts = (m[2] !== undefined ? m[2] : m[3] !== undefined ? m[3] : '').split(',').map(s => s.trim()).filter(s => s && !TS_OPT.test(s));
        out.push(m[1] + (!opts.length ? '' : opts.length === 1 && opts[0].indexOf('=') < 0 ? '=' + opts[0] : '=(' + opts.join(',') + ')'));
        hasOpt = true;
        continue;
      }
      out.push(t);
    }
    if (!hasOpt) out.splice(out.length && out[0][0] === '#' ? 1 : 0, 0, 'opt');
    return out.join(' ');
  }
  /** ORCA keywords for the same purpose: OptTS, ScanTS, NEB-TS, NEB-CI, NEB and IRC become Opt (added when
   *  missing); the rest is kept. Accepts a string ("! B3LYP def2-SVP OptTS") or an array of keywords. */
  function minimumKeywords(keywords) {
    const toks = (Array.isArray(keywords) ? keywords.join(' ') : String(keywords || '')).replace(/!/g, ' ').split(/\s+/).filter(Boolean), out = [];
    let opt = false;
    for (const t of toks) {
      const lo = t.toLowerCase();
      if (/^(optts|scants|neb-ts|neb-ci|zoom-neb-ts|zoom-neb-ci|neb|irc)$/.test(lo)) { if (!opt) { out.push('Opt'); opt = true; } continue; }
      if (/^(opt|copt|zopt|gdiis-opt|tightopt|verytightopt|looseopt|l-opt|lopt)$/.test(lo)) { if (opt) continue; opt = true; }
      out.push(t);
    }
    if (!opt) out.push('Opt');
    return out.join(' ');
  }

  BM.vibtools = { SMALL, stationary, describe, listImag, ends, bondChanges, endMeshes, minimumRoute, minimumKeywords, routeTokens };
})(typeof window !== 'undefined' ? (window.BM = window.BM || {}) : (globalThis.BM = globalThis.BM || {}));
