/* BlueMolecule — molecular orbitals and electron densities computed from the basis set and MO coefficients of an
 * fchk or Molden file (ds.basis and ds.mos, see parsers.js), without cubegen or orca_plot.
 *
 *   BM.orbitals.model(ds)                -> model, or null without a usable basis (cached on the dataset)
 *   BM.orbitals.detect(model)            -> {conv, error, tried}: the file's normalisation convention
 *   BM.orbitals.overlap(model, prim)     -> Float64Array(nbf * nbf): the AO overlap matrix
 *   BM.orbitals.list(model)              -> [{spin, index, energy, occ, sym, tag}]: the orbitals, for menus
 *   BM.orbitals.values(model, what, pts) -> Float64Array: values at points (Bohr, x y z one after the other)
 *   BM.orbitals.gridFor(model, opts)     -> {origin, axes, dims} (Bohr): a box around the molecule
 *   BM.orbitals.compute(model, what, grid, opts) -> {promise (the volume, in Å), cancel()}
 *   BM.orbitals.evaluate(model, what, grid)      -> Float32Array: the same, at once (tests, small grids)
 *   what: {kind: 'mo', spin: 'alpha' | 'beta', index (0-based)} | {kind: 'density'} | {kind: 'spin'}
 *
 * An AO here is a contracted Cartesian component x^a y^b z^c, each component normalised on its own, or a normalised
 * real solid harmonic in the order m = 0, +1, -1, +2, -2, … with Gaussian's signs. Programs write their MOs in other
 * conventions:
 *   prim  'normalized': the contraction coefficients multiply normalised primitives (Gaussian, the Molden standard);
 *         'raw': they include the primitive normalisation (ORCA, Psi4 before 1.0). Contractions are always
 *         renormalised, as Molden does (Turbomole writes unnormalised Cartesian d contractions).
 *   cart  'component': each Cartesian component normalised on its own (Gaussian, Molden);
 *         'axis': every component normalised like x^l (Psi4 1.0 to 1.3.2);
 *         'radial': without the double-factorial factors (CFOUR);
         'scaled': each component times sqrt((2l-1)!!) (Turbomole's tm2molden).
 *   sign  'standard', or 'orca': ORCA's pure functions with |m| = 3 and 4 have the opposite sign.
 * detect() takes the convention in which the file's MOs are orthonormal (C^T S C = 1), which also checks the reading.
 * The grid is evaluated in Web Workers when the page can start them (blob URLs), otherwise in slices on the page.
 */
(function (BM) {
  'use strict';

  const BOHR = 0.529177210903;   // Å per Bohr
  const dfact = n => { let r = 1; for (let k = n; k > 1; k -= 2) r *= k; return r; };
  const binom = (n, k) => { if (k < 0 || k > n) return 0; let r = 1; for (let i = 1; i <= k; i++) r = r * (n - k + i) / i; return r; };
  const count = (s, ch) => s.split(ch).length - 1;
  const MOLDEN_G = ['xxxx', 'yyyy', 'zzzz', 'xxxy', 'xxxz', 'yyyx', 'yyyz', 'zzzx', 'zzzy', 'xxyy', 'xxzz', 'yyzz', 'xxyz', 'yyxz', 'zzxy']
    .map(s => [count(s, 'x'), count(s, 'y'), count(s, 'z')]);
  const DEFAULT = { prim: 'normalized', cart: 'component', sign: 'standard' };

  /** Cartesian components [a, b, c] of a shell in the order its file lists them. Gaussian (and Molden up to f):
   *  x y z; xx yy zz xy xz yz; xxx yyy zzz xyy xxy xxz xzz yzz yyz xyz; g and higher zzzz yzzz yyzz … xxxx.
   *  Molden's g: xxxx yyyy zzzz xxxy xxxz yyyx yyyz zzzx zzzy xxyy xxzz yyzz xxyz yyxz zzxy. */
  function cartOrder(l, source) {
    if (l === 0) return [[0, 0, 0]];
    if (l === 1) return [[1, 0, 0], [0, 1, 0], [0, 0, 1]];
    if (l === 2) return [[2, 0, 0], [0, 2, 0], [0, 0, 2], [1, 1, 0], [1, 0, 1], [0, 1, 1]];
    if (l === 3) return [[3, 0, 0], [0, 3, 0], [0, 0, 3], [1, 2, 0], [2, 1, 0], [2, 0, 1], [1, 0, 2], [0, 1, 2], [0, 2, 1], [1, 1, 1]];
    if (l === 4 && source === 'molden') return MOLDEN_G.map(c => c.slice());
    const out = [];
    for (let a = 0; a <= l; a++) for (let b = 0; b <= l - a; b++) out.push([a, b, l - a - b]);
    return out;
  }
  /** Orders m of a pure shell: 0, +1, -1, +2, -2, … */
  function mOrder(l) { const m = [0]; for (let k = 1; k <= l; k++) m.push(k, -k); return m; }

  /** Real solid harmonic of degree l and order m as a map 'a,b,c' -> coefficient (unnormalised; Helgaker, Jørgensen
   *  and Olsen, Molecular Electronic-Structure Theory, eq. 6.4.48): m > 0 cosine-like, m < 0 sine-like. */
  function solidHarmonic(l, m) {
    const am = Math.abs(m), v0 = m < 0 ? 1 : 0, out = new Map();
    for (let t = 0; t <= (l - am) >> 1; t++) {
      for (let u = 0; u <= t; u++) {
        for (let v2 = v0; v2 <= am; v2 += 2) {   // v2 = 2v
          const sgn = (t + (v2 - v0) / 2) % 2 ? -1 : 1;
          const c = sgn * Math.pow(0.25, t) * binom(l, t) * binom(l - t, am + t) * binom(t, u) * binom(am, v2);
          if (!c) continue;
          const key = (2 * t + am - 2 * u - v2) + ',' + (2 * u + v2) + ',' + (l - 2 * t - am);
          out.set(key, (out.get(key) || 0) + c);
        }
      }
    }
    return out;
  }
  /** <x^a y^b z^c | x^a' y^b' z^c'> for one radial part, up to a factor common to the shell. */
  function gram(a, b) {
    const s0 = a[0] + b[0], s1 = a[1] + b[1], s2 = a[2] + b[2];
    return s0 % 2 || s1 % 2 || s2 % 2 ? 0 : dfact(s0 - 1) * dfact(s1 - 1) * dfact(s2 - 1);
  }
  /** T[f][c]: pure function f (order mOrder) as a combination of the normalised Cartesian components comps[c]. */
  function pureMatrix(l, comps) {
    const idx = new Map(comps.map((c, k) => [c.join(','), k]));
    return mOrder(l).map(m => {
      const P = solidHarmonic(l, m), row = new Float64Array(comps.length);
      let q = 0;
      for (const [k1, v1] of P) for (const [k2, v2] of P) q += v1 * v2 * gram(k1.split(',').map(Number), k2.split(',').map(Number));
      for (const [k, v] of P) { const c = idx.get(k); row[c] = v * Math.sqrt(gram(comps[c], comps[c]) / q); }
      return row;
    });
  }
  const TCACHE = new Map();
  function pureT(l) { if (!TCACHE.has(l)) TCACHE.set(l, pureMatrix(l, cartOrder(l, 'fchk'))); return TCACHE.get(l); }

  // ------------------------------------------------------------------ the model
  function model(ds) {
    if (!ds || !ds.basis || !ds.mos || !ds.mos.alpha || !Array.isArray(ds.basis.shells)) return null;
    if (ds._orbModel !== undefined) return ds._orbModel;
    let m = null;
    try { m = build(ds); } catch (e) { m = null; }
    Object.defineProperty(ds, '_orbModel', { value: m, writable: true, configurable: true, enumerable: false });
    return m;
  }
  function build(ds) {
    const B = ds.basis, src = B.source === 'molden' ? 'molden' : 'fchk';
    const shells = [];
    let off = 0;
    for (const s of B.shells) {
      const l = s.l | 0, pure = !!s.pure && l >= 2;
      if (l < 0 || l > 7 || !s.exps || !s.exps.length || s.exps.length !== s.coefs.length) return null;
      const comps = pure ? cartOrder(l, 'fchk') : cartOrder(l, src);
      const cn = comps.map(c => 1 / Math.sqrt(dfact(2 * c[0] - 1) * dfact(2 * c[1] - 1) * dfact(2 * c[2] - 1)));
      const nf = pure ? 2 * l + 1 : comps.length;
      shells.push({ atom: s.atom, center: Float64Array.from(s.center), l, pure, exps: Float64Array.from(s.exps), coefs: Float64Array.from(s.coefs),
        comps, cn, T: pure ? pureT(l) : null, off, nf });
      off += nf;
    }
    const nbf = off;
    if (nbf !== B.nbf) return null;
    const sets = [['alpha', ds.mos.alpha], ['beta', ds.mos.beta]].filter(x => x[1] && x[1].coeffs && x[1].coeffs.length);
    for (const [, s] of sets) if (s.coeffs.length % nbf) return null;
    const M = { nbf, shells, source: src, program: B.program || null, restricted: !!ds.mos.restricted || sets.length === 1, sets: {}, conv: null, detected: null };
    for (const [name, s] of sets) {
      const n = s.coeffs.length / nbf;
      M.sets[name] = { name, n, coeffs: s.coeffs, energies: s.energies || [], occupations: s.occupations || [], symmetries: s.symmetries || null, fileIndex: s.fileIndex || null };
    }
    return M;
  }
  /** Normalised radial weights of a shell: R(r) = sum_k w_k exp(-a_k r^2), times cn[c] x^a y^b z^c per component. */
  function weights(sh, prim) {
    const l = sh.l, n = sh.exps.length, w = new Float64Array(n);
    for (let k = 0; k < n; k++) {
      const a = sh.exps[k];
      w[k] = prim === 'raw' ? sh.coefs[k] / (Math.pow(2 * a / Math.PI, 0.75) * Math.pow(4 * a, l / 2)) : sh.coefs[k];
    }
    let S = 0;
    for (let i = 0; i < n; i++) for (let j = 0; j < n; j++) {
      const ai = sh.exps[i], aj = sh.exps[j];
      S += w[i] * w[j] * Math.pow(2 * Math.sqrt(ai * aj) / (ai + aj), l + 1.5);
    }
    if (!(S > 0) || !Number.isFinite(S)) throw new Error('a shell of the basis set cannot be normalised');
    const f = 1 / Math.sqrt(S);
    for (let k = 0; k < n; k++) { const a = sh.exps[k]; w[k] *= f * Math.pow(2 * a / Math.PI, 0.75) * Math.pow(4 * a, l / 2); }
    return w;
  }

  // ------------------------------------------------------------------ overlap (Obara–Saika)
  function fill1D(E, la, lb, PA, PB, p) {
    const n = lb + 1, h = 1 / (2 * p);
    E[0] = Math.sqrt(Math.PI / p);
    for (let i = 0; i <= la; i++) {
      for (let j = 0; j <= lb; j++) {
        if (!i && !j) continue;
        let v;
        if (i > 0) v = PA * E[(i - 1) * n + j] + (i > 1 ? h * (i - 1) * E[(i - 2) * n + j] : 0) + (j > 0 ? h * j * E[(i - 1) * n + j - 1] : 0);
        else v = PB * E[j - 1] + (j > 1 ? h * (j - 1) * E[j - 2] : 0);
        E[i * n + j] = v;
      }
    }
  }
  /** Overlaps of the Cartesian components of two shells (component normalisation included), row-major na x nb. */
  function cartOverlap(sa, wa, sb, wb) {
    const la = sa.l, lb = sb.l, A = sa.center, B = sb.center, na = sa.comps.length, nb = sb.comps.length;
    const out = new Float64Array(na * nb), n = lb + 1;
    const Ex = new Float64Array((la + 1) * n), Ey = new Float64Array((la + 1) * n), Ez = new Float64Array((la + 1) * n);
    const dx = A[0] - B[0], dy = A[1] - B[1], dz = A[2] - B[2], ab2 = dx * dx + dy * dy + dz * dz;
    for (let i = 0; i < sa.exps.length; i++) {
      const a = sa.exps[i];
      for (let j = 0; j < sb.exps.length; j++) {
        const b = sb.exps[j], p = a + b, K = Math.exp(-a * b / p * ab2) * wa[i] * wb[j];
        if (Math.abs(K) * Math.pow(Math.PI / p, 1.5) < 1e-18) continue;
        const Px = (a * A[0] + b * B[0]) / p, Py = (a * A[1] + b * B[1]) / p, Pz = (a * A[2] + b * B[2]) / p;
        fill1D(Ex, la, lb, Px - A[0], Px - B[0], p);
        fill1D(Ey, la, lb, Py - A[1], Py - B[1], p);
        fill1D(Ez, la, lb, Pz - A[2], Pz - B[2], p);
        for (let ca = 0; ca < na; ca++) {
          const c1 = sa.comps[ca];
          for (let cb = 0; cb < nb; cb++) {
            const c2 = sb.comps[cb];
            out[ca * nb + cb] += K * Ex[c1[0] * n + c2[0]] * Ey[c1[1] * n + c2[1]] * Ez[c1[2] * n + c2[2]];
          }
        }
      }
    }
    for (let ca = 0; ca < na; ca++) for (let cb = 0; cb < nb; cb++) out[ca * nb + cb] *= sa.cn[ca] * sb.cn[cb];
    return out;
  }
  /** The AO overlap matrix in this module's convention, for the primitive convention prim. */
  function overlap(M, prim) {
    const nbf = M.nbf, S = new Float64Array(nbf * nbf), sh = M.shells;
    const W = sh.map(s => weights(s, prim || 'normalized'));
    // the most diffuse exponent of each shell, to skip pairs that do not overlap
    const amin = sh.map(s => Math.min.apply(null, Array.from(s.exps)));
    for (let A = 0; A < sh.length; A++) {
      for (let B = A; B < sh.length; B++) {
        const a = sh[A], b = sh[B];
        const dx = a.center[0] - b.center[0], dy = a.center[1] - b.center[1], dz = a.center[2] - b.center[2];
        if (amin[A] * amin[B] / (amin[A] + amin[B]) * (dx * dx + dy * dy + dz * dz) > 46) continue;   // e^-46 ~ 1e-20
        const C = cartOverlap(a, W[A], b, W[B]), nca = a.comps.length, ncb = b.comps.length;
        for (let f = 0; f < a.nf; f++) {
          for (let g = 0; g < b.nf; g++) {
            let v = 0;
            if (!a.pure && !b.pure) v = C[f * ncb + g];
            else {
              for (let p = 0; p < nca; p++) {
                const ta = a.pure ? a.T[f][p] : (p === f ? 1 : 0);
                if (!ta) continue;
                for (let q = 0; q < ncb; q++) {
                  const tb = b.pure ? b.T[g][q] : (q === g ? 1 : 0);
                  if (tb) v += ta * tb * C[p * ncb + q];
                }
              }
            }
            S[(a.off + f) * nbf + b.off + g] = v;
            S[(b.off + g) * nbf + a.off + f] = v;
          }
        }
      }
    }
    return S;
  }

  // ------------------------------------------------------------------ conventions
  /** Factors that turn the file's MO coefficients into this module's convention (per AO, in file order). */
  function convFactors(M, conv) {
    const f = new Float64Array(M.nbf).fill(1);
    for (const sh of M.shells) {
      if (!sh.pure && sh.l >= 2 && conv.cart !== 'component') {
        sh.comps.forEach((c, k) => {
          const g = dfact(2 * c[0] - 1) * dfact(2 * c[1] - 1) * dfact(2 * c[2] - 1);
          f[sh.off + k] = conv.cart === 'axis' ? Math.sqrt(g / dfact(2 * sh.l - 1)) : conv.cart === 'radial' ? Math.sqrt(g) : Math.sqrt(dfact(2 * sh.l - 1));
        });
      }
      if (sh.pure && sh.l >= 3 && conv.sign === 'orca') mOrder(sh.l).forEach((m, k) => { if (Math.abs(m) === 3 || Math.abs(m) === 4) f[sh.off + k] = -1; });
    }
    return f;
  }
  /** Indices of the orbitals used to test orthonormality: all of a small set, otherwise the lowest, those around the
   *  frontier and the highest. */
  function sample(set, max) {
    const n = set.n;
    if (n <= max) return Array.from({ length: n }, (_, i) => i);
    let h = -1;
    for (let i = 0; i < n; i++) if (set.occupations[i] > 0.5) h = i;
    const pick = new Set();
    for (let i = 0; i < 6; i++) pick.add(i);
    for (let i = h - 9; i <= h + 10; i++) if (i >= 0 && i < n) pick.add(i);
    for (let i = n - 6; i < n; i++) pick.add(i);
    for (let k = 0; pick.size < max && k < n; k += Math.max(1, Math.floor(n / max))) pick.add(k);
    return Array.from(pick).sort((a, b) => a - b);
  }
  /** Largest |C^T S C - 1| over the sampled orbitals of every set, with the coefficients scaled by fac. */
  function deviation(M, S, fac, picks) {
    const nbf = M.nbf;
    let worst = 0;
    for (const name of Object.keys(M.sets)) {
      const set = M.sets[name], idx = picks[name];
      const V = idx.map(m => { const v = new Float64Array(nbf); for (let u = 0; u < nbf; u++) v[u] = set.coeffs[m * nbf + u] * fac[u]; return v; });
      const SV = V.map(v => {
        const r = new Float64Array(nbf);
        for (let u = 0; u < nbf; u++) { const vu = v[u]; if (!vu) continue; const row = u * nbf; for (let w = 0; w < nbf; w++) r[w] += S[row + w] * vu; }
        return r;
      });
      for (let a = 0; a < V.length; a++) {
        for (let b = a; b < V.length; b++) {
          let d = 0;
          const x = V[a], y = SV[b];
          for (let u = 0; u < nbf; u++) d += x[u] * y[u];
          const e = Math.abs(d - (a === b ? 1 : 0));
          if (!(e <= worst)) worst = Number.isFinite(e) ? Math.max(worst, e) : Infinity;
        }
      }
    }
    return worst;
  }
  const convKey = c => c.prim + '/' + c.cart + '/' + c.sign;
  /** The convention in which the MOs are orthonormal. Tries the standard one first; the others only when it fails. */
  function detect(M) {
    if (M.detected) return M.detected;
    const sh = M.shells;
    const first = M.program === 'orca' ? { prim: 'raw', cart: 'component', sign: 'orca' } : DEFAULT;
    if (M.nbf > 4000) {   // the overlap matrix alone would take hundreds of megabytes: the program's convention, unchecked
      M.conv = Object.assign({}, first);
      M.detected = { conv: M.conv, error: NaN, tried: [], skipped: true, relevant: {} };
      return M.detected;
    }
    const prims = sh.some(s => s.exps.length > 1) ? ['normalized', 'raw'] : ['normalized'];
    const carts = sh.some(s => !s.pure && s.l >= 2) ? ['component', 'axis', 'radial', 'scaled'] : ['component'];
    const signs = sh.some(s => s.pure && s.l >= 3) ? ['standard', 'orca'] : ['standard'];
    const max = M.nbf > 1500 ? 10 : M.nbf > 600 ? 16 : 24, picks = {};
    for (const name of Object.keys(M.sets)) picks[name] = sample(M.sets[name], max);
    const tried = [];
    const Scache = {};
    const S = prim => Scache[prim] || (Scache[prim] = overlap(M, prim));
    const test = c => { const e = deviation(M, S(c.prim), convFactors(M, c), picks); tried.push({ conv: c, error: e }); return e; };
    // the program's own convention first, when the file says which program wrote it
    let best = { conv: Object.assign({}, first), error: test(first) };
    if (!(best.error <= 1e-5)) {
      for (const prim of prims) for (const cart of carts) for (const sign of signs) {
        const c = { prim, cart, sign };
        if (convKey(c) === convKey(first)) continue;
        const e = test(c);
        if (e < best.error * 0.5) best = { conv: c, error: e };
      }
    }
    M.conv = best.conv;
    M.detected = { conv: best.conv, error: best.error, tried, relevant: { prim: prims.length > 1, cart: carts.length > 1, sign: signs.length > 1 } };
    return M.detected;
  }
  /** |C^T S C - 1| of the sampled orbitals in a given convention (any of prim, cart, sign; the others standard). */
  function check(M, conv) {
    const c = Object.assign({}, DEFAULT, conv), max = M.nbf > 1500 ? 10 : M.nbf > 600 ? 16 : 24, picks = {};
    for (const name of Object.keys(M.sets)) picks[name] = sample(M.sets[name], max);
    return deviation(M, overlap(M, c.prim), convFactors(M, c), picks);
  }
  /** Words for a convention (for the interface and warnings). */
  function describe(det) {
    const c = det.conv, r = det.relevant || {}, parts = [];
    if (r.prim && c.prim === 'raw') parts.push('contraction coefficients that include the primitive normalisation');
    if (r.cart && c.cart === 'axis') parts.push('Cartesian functions normalised like x^l');
    if (r.cart && c.cart === 'radial') parts.push('Cartesian functions without their double-factorial normalisation (CFOUR)');
    if (r.cart && c.cart === 'scaled') parts.push('Cartesian functions scaled as Turbomole writes them');
    if (r.sign && c.sign === 'orca') parts.push('ORCA’s signs for f and g functions');
    return parts.length ? parts.join(', ') : 'the standard convention';
  }

  // ------------------------------------------------------------------ orbitals for menus
  function list(M) {
    const out = [];
    for (const name of ['alpha', 'beta']) {
      const s = M.sets[name];
      if (!s) continue;
      let homo = -1;
      for (let i = 0; i < s.n; i++) if (s.occupations[i] > 0.5) homo = i;
      for (let i = 0; i < s.n; i++) {
        const d = i - homo;
        const tag = homo < 0 ? '' : d === 0 ? 'HOMO' : d === 1 ? 'LUMO' : d < 0 ? (d >= -9 ? 'HOMO−' + -d : '') : (d <= 10 ? 'LUMO+' + (d - 1) : '');
        out.push({ spin: name, index: i, energy: s.energies[i], occ: s.occupations[i], sym: s.symmetries ? s.symmetries[i] : null, tag, homo: d === 0, lumo: d === 1 });
      }
    }
    return out;
  }
  function frontier(M, spin) {
    const s = M.sets[spin || 'alpha'];
    if (!s) return { homo: -1, lumo: -1 };
    let homo = -1;
    for (let i = 0; i < s.n; i++) if (s.occupations[i] > 0.5) homo = i;
    return { homo, lumo: homo + 1 < s.n ? homo + 1 : -1 };
  }
  /** Whether the spin density is defined: unrestricted orbitals, or restricted ones with singly occupied orbitals. */
  function hasSpin(M) {
    if (M.sets.beta) return true;
    const o = M.sets.alpha.occupations;
    return Array.from(o).some(x => Math.abs(x - 1) < 1e-6);
  }

  // ------------------------------------------------------------------ jobs for the kernel
  /** Plain data for evalBlock: the shells, and the orbitals to evaluate as Cartesian coefficients. */
  function job(M, what) {
    if (!M.conv) detect(M);
    const conv = M.conv, fac = convFactors(M, conv), nbf = M.nbf;
    const orbs = [];   // [{set, index, weight}]
    const kind = what.kind;
    if (kind === 'mo') {
      const s = M.sets[what.spin || 'alpha'];
      if (!s || !(what.index >= 0 && what.index < s.n)) throw new Error('There is no such orbital.');
      orbs.push({ set: s, index: what.index, weight: 1 });
    } else if (kind === 'density' || kind === 'spin') {
      for (const name of Object.keys(M.sets)) {
        const s = M.sets[name];
        for (let i = 0; i < s.n; i++) {
          const o = s.occupations[i];
          if (!(Math.abs(o) > 1e-10)) continue;
          if (kind === 'density') orbs.push({ set: s, index: i, weight: o });
          else if (M.sets.beta) orbs.push({ set: s, index: i, weight: name === 'beta' ? -o : o });
          else if (Math.abs(o - 1) < 1e-6) orbs.push({ set: s, index: i, weight: 1 });
        }
      }
      if (!orbs.length) throw new Error(kind === 'spin' ? 'These orbitals have no unpaired electrons.' : 'No orbital is occupied.');
    } else throw new Error('Unknown quantity: ' + kind);
    const sh = M.shells, ns = sh.length;
    let nc = 0, np = 0, maxL = 0;
    for (const s of sh) { nc += s.comps.length; np += s.exps.length; maxL = Math.max(maxL, s.l); }
    const J = {
      ns, maxL, nso: orbs.length, nc, density: kind !== 'mo',
      center: new Float64Array(3 * ns), l: new Int32Array(ns), poff: new Int32Array(ns + 1), coff: new Int32Array(ns + 1),
      exp: new Float64Array(np), w: new Float64Array(np), abc: new Int32Array(3 * nc), rcut2: new Float64Array(ns),
      D: new Float64Array(orbs.length * nc), wt: Float64Array.from(orbs.map(o => o.weight))
    };
    let p = 0, c = 0;
    sh.forEach((s, k) => {
      const W = weights(s, conv.prim);
      J.center.set(s.center, 3 * k); J.l[k] = s.l; J.poff[k] = p; J.coff[k] = c;
      let r2 = 0;
      for (let q = 0; q < s.exps.length; q++) {
        J.exp[p + q] = s.exps[q]; J.w[p + q] = W[q];
        // where |w| e^(-a r^2) r^l falls below 1e-9
        const a = s.exps[q], lw = Math.log(Math.abs(W[q]) / 1e-9);
        let rr = Math.max(1, lw / a);
        for (let it = 0; it < 8; it++) rr = Math.max(1e-6, (lw + s.l * 0.5 * Math.log(rr)) / a);
        if (lw > 0) r2 = Math.max(r2, rr);
      }
      J.rcut2[k] = r2;
      s.comps.forEach((cc, q) => { J.abc[3 * (c + q)] = cc[0]; J.abc[3 * (c + q) + 1] = cc[1]; J.abc[3 * (c + q) + 2] = cc[2]; });
      orbs.forEach((o, oi) => {
        const base = o.index * nbf;
        for (let q = 0; q < s.comps.length; q++) {
          let v = 0;
          if (s.pure) { for (let f = 0; f < s.nf; f++) v += s.T[f][q] * o.set.coeffs[base + s.off + f] * fac[s.off + f]; }
          else v = o.set.coeffs[base + s.off + q] * fac[s.off + q];
          J.D[oi * nc + c + q] = v * s.cn[q];
        }
      });
      p += s.exps.length; c += s.comps.length;
    });
    J.poff[ns] = p; J.coff[ns] = c;
    return J;
  }

  /** The kernel (also run in workers, so it uses nothing from outside): values on the planes i0 <= i < i1 of a grid
   *  (Bohr) with origin o, step vectors ax (9 numbers) and dims [n1, n2, n3]; returns Float32Array((i1-i0) n2 n3).
   *  For each line of the grid along its third axis, the shells whose cut-off sphere it crosses are evaluated over
   *  that chord only. */
  function evalBlock(J, o, ax, dims, i0, i1, tick) {
    const n2 = dims[1], n3 = dims[2], nso = J.nso, nc = J.nc, ns = J.ns, maxL = J.maxL;
    const out = new Float32Array((i1 - i0) * n2 * n3);
    const psi = new Float64Array(nso * n3);
    const px = new Float64Array(maxL + 1), py = new Float64Array(maxL + 1), pz = new Float64Array(maxL + 1);
    const a3x = ax[6], a3y = ax[7], a3z = ax[8], a3n = a3x * a3x + a3y * a3y + a3z * a3z;
    const C = J.center, L = J.l, E = J.exp, W = J.w, P = J.poff, CO = J.coff, ABC = J.abc, D = J.D, R2 = J.rcut2, wt = J.wt;
    for (let i = i0; i < i1; i++) {
      for (let j = 0; j < n2; j++) {
        const bx = o[0] + i * ax[0] + j * ax[3], by = o[1] + i * ax[1] + j * ax[4], bz = o[2] + i * ax[2] + j * ax[5];
        psi.fill(0);
        for (let s = 0; s < ns; s++) {
          const ux = C[3 * s] - bx, uy = C[3 * s + 1] - by, uz = C[3 * s + 2] - bz;
          const t = (ux * a3x + uy * a3y + uz * a3z) / a3n;
          const d2 = ux * ux + uy * uy + uz * uz - t * t * a3n;
          const rc = R2[s];
          if (d2 > rc) continue;
          const half = Math.sqrt((rc - d2) / a3n);
          const k0 = Math.max(0, Math.ceil(t - half)), k1 = Math.min(n3 - 1, Math.floor(t + half));
          if (k0 > k1) continue;
          const l = L[s], p0 = P[s], p1 = P[s + 1], c0 = CO[s], c1 = CO[s + 1];
          for (let k = k0; k <= k1; k++) {
            const x = bx + k * a3x - C[3 * s], y = by + k * a3y - C[3 * s + 1], z = bz + k * a3z - C[3 * s + 2];
            const r2 = x * x + y * y + z * z;
            if (r2 > rc) continue;
            let rad = 0;
            for (let q = p0; q < p1; q++) rad += W[q] * Math.exp(-E[q] * r2);
            if (l === 0) {
              for (let so = 0; so < nso; so++) psi[so * n3 + k] += rad * D[so * nc + c0];
            } else if (l === 1) {   // (x, y, z in every format)
              for (let so = 0; so < nso; so++) { const b = so * nc + c0; psi[so * n3 + k] += rad * (D[b] * x + D[b + 1] * y + D[b + 2] * z); }
            } else {
              px[0] = py[0] = pz[0] = 1;
              for (let e = 1; e <= l; e++) { px[e] = px[e - 1] * x; py[e] = py[e - 1] * y; pz[e] = pz[e - 1] * z; }
              for (let so = 0; so < nso; so++) {
                let v = 0;
                const b = so * nc;
                for (let cc = c0; cc < c1; cc++) v += D[b + cc] * px[ABC[3 * cc]] * py[ABC[3 * cc + 1]] * pz[ABC[3 * cc + 2]];
                psi[so * n3 + k] += rad * v;
              }
            }
          }
        }
        const base = ((i - i0) * n2 + j) * n3;
        if (J.density) {
          for (let k = 0; k < n3; k++) { let v = 0; for (let so = 0; so < nso; so++) { const f = psi[so * n3 + k]; v += wt[so] * f * f; } out[base + k] = v; }
        } else {
          for (let k = 0; k < n3; k++) out[base + k] = psi[k];
        }
      }
      if (tick) tick(i);
    }
    return out;
  }

  /** Values at points (Bohr): each point is a grid of one point. */
  function values(M, what, pts) {
    const J = job(M, what), n = pts.length / 3, out = new Float64Array(n);
    const ax = [1, 0, 0, 0, 1, 0, 0, 0, 1];
    for (let k = 0; k < n; k++) out[k] = evalBlock(J, [pts[3 * k], pts[3 * k + 1], pts[3 * k + 2]], ax, [1, 1, 1], 0, 1)[0];
    return out;
  }

  /** A box around the basis centres: margin and spacing in Bohr (defaults 7.5 and 0.35), at most maxPoints points
   *  (spacing widened to fit). */
  function gridFor(M, opts) {
    opts = opts || {};
    const margin = opts.margin > 0 ? opts.margin : 7.5;
    let h = opts.spacing > 0 ? opts.spacing : 0.35;
    const lo = [Infinity, Infinity, Infinity], hi = [-Infinity, -Infinity, -Infinity];
    for (const s of M.shells) for (let a = 0; a < 3; a++) { lo[a] = Math.min(lo[a], s.center[a]); hi[a] = Math.max(hi[a], s.center[a]); }
    const ext = [0, 1, 2].map(a => hi[a] - lo[a] + 2 * margin);
    const maxP = opts.maxPoints > 0 ? opts.maxPoints : 4e6;
    let dims = ext.map(e => Math.max(2, Math.ceil(e / h) + 1));
    while (dims[0] * dims[1] * dims[2] > maxP) { h *= 1.05; dims = ext.map(e => Math.max(2, Math.ceil(e / h) + 1)); }
    const origin = [0, 1, 2].map(a => (lo[a] + hi[a]) / 2 - (dims[a] - 1) * h / 2);
    return { origin, axes: [[h, 0, 0], [0, h, 0], [0, 0, h]], dims, spacing: h };
  }

  function label(M, what) {
    if (what.kind === 'density') return 'Electron density';
    if (what.kind === 'spin') return 'Spin density';
    const s = M.sets[what.spin || 'alpha'], it = list(M).find(x => x.spin === s.name && x.index === what.index);
    const spin = M.sets.beta ? (s.name === 'beta' ? ' β' : ' α') : '';
    return 'MO ' + (what.index + 1) + spin + (it && it.tag ? ' (' + it.tag + ')' : '');
  }
  function volumeOf(M, what, grid, data) {
    const b = BOHR;
    return { name: label(M, what), kind: what.kind === 'mo' ? 'orbital' : what.kind === 'spin' ? 'spin' : 'density',
      origin: grid.origin.map(x => x * b), axes: grid.axes.map(v => v.map(x => x * b)), dims: grid.dims.slice(), data,
      computed: { what: Object.assign({}, what), spacing: grid.spacing } };
  }
  /** The grid at once, on this thread. */
  function evaluate(M, what, grid) {
    const J = job(M, what), ax = [].concat(grid.axes[0], grid.axes[1], grid.axes[2]);
    return evalBlock(J, grid.origin, ax, grid.dims, 0, grid.dims[0]);
  }

  // ------------------------------------------------------------------ workers
  let workerURL;
  function workerSource() {
    return '"use strict";\nconst evalBlock = ' + evalBlock.toString() + ';\nlet J = null;\n' +
      'self.onmessage = e => { const m = e.data; if (m.job) { J = m.job; self.postMessage({ ready: true }); return; }\n' +
      '  const out = evalBlock(J, m.o, m.ax, m.dims, m.i0, m.i1);\n' +
      '  self.postMessage({ i0: m.i0, i1: m.i1, out }, [out.buffer]); };\n';
  }
  function makeWorker() {
    if (typeof Worker === 'undefined' || typeof Blob === 'undefined' || typeof URL === 'undefined' || !URL.createObjectURL) return null;
    try {
      if (!workerURL) workerURL = URL.createObjectURL(new Blob([workerSource()], { type: 'text/javascript' }));
      return new Worker(workerURL);
    } catch (e) { return null; }
  }
  /** The grid in the background: {promise, cancel}. opts.onProgress(fraction); opts.workers: how many (0: on the
   *  page). The promise resolves with the volume (Å) or rejects with an Error ('cancelled' when cancelled). */
  function compute(M, what, grid, opts) {
    opts = opts || {};
    const J = job(M, what), ax = [].concat(grid.axes[0], grid.axes[1], grid.axes[2]), dims = grid.dims;
    const n1 = dims[0], plane = dims[1] * dims[2], data = new Float32Array(n1 * plane);
    const hc = typeof navigator !== 'undefined' && navigator.hardwareConcurrency ? navigator.hardwareConcurrency : 2;
    const nw = opts.workers !== undefined ? opts.workers : Math.max(1, Math.min(8, hc - 1));
    // slabs of a few planes, handed out as workers finish (the work per plane varies)
    const step = Math.max(1, Math.min(8, Math.ceil(n1 / Math.max(1, nw) / 6)));
    const slabs = [];
    for (let i = 0; i < n1; i += step) slabs.push([i, Math.min(n1, i + step)]);
    let done = 0, cancelled = false, workers = [], timer = 0, reject = null;
    const progress = () => { if (opts.onProgress) try { opts.onProgress(done / n1); } catch (e) { /* ignore */ } };
    const finish = () => { workers.forEach(w => w.terminate()); workers = []; };
    const promise = new Promise((resolve, rej) => {
      reject = rej;
      const queue = slabs.slice(), busy = new Map();   // slabs not yet done; slab given to each worker
      let used = 0;
      const ok = () => { finish(); const v = volumeOf(M, what, grid, data); v.computed.workers = used; resolve(v); };
      const store = (i0, i1, out) => { data.set(out, i0 * plane); done += i1 - i0; progress(); };
      // on the page, one slab at a time (no workers, or workers that failed: their slabs go back to the queue)
      const run = () => {
        if (cancelled) return;
        const t0 = Date.now();
        try {
          while (queue.length && Date.now() - t0 < 40) { const [i0, i1] = queue.shift(); store(i0, i1, evalBlock(J, grid.origin, ax, dims, i0, i1)); }
        } catch (e) { rej(e); return; }
        if (done >= n1) ok(); else timer = setTimeout(run, 0);
      };
      for (let k = 0; k < nw; k++) { const w = makeWorker(); if (!w) break; workers.push(w); }
      used = workers.length;
      if (!workers.length) { timer = setTimeout(run, 0); return; }
      let failed = false;
      // a worker that never answers (a web view that starts workers but cannot run them): the page takes over
      const fallBack = () => {
        if (failed || cancelled) return;
        failed = true;
        for (const sl of busy.values()) queue.push(sl);
        busy.clear(); finish(); used = 0;
        timer = setTimeout(run, 0);
      };
      let alive = false;
      const watchdog = setTimeout(() => { if (!alive) fallBack(); }, opts.watchdogMs > 0 ? opts.watchdogMs : 4000);
      const feed = w => {
        if (!queue.length) return;
        const sl = queue.shift();
        busy.set(w, sl);
        w.postMessage({ o: grid.origin, ax, dims, i0: sl[0], i1: sl[1] });
      };
      workers.forEach(w => {
        w.onmessage = e => {
          if (cancelled || failed) return;
          const m = e.data;
          if (m.ready) { if (!alive) { alive = true; clearTimeout(watchdog); } return; }
          busy.delete(w);
          store(m.i0, m.i1, m.out);
          if (done >= n1) ok(); else feed(w);
        };
        w.onerror = e => {
          if (e && e.preventDefault) e.preventDefault();
          clearTimeout(watchdog);
          fallBack();
        };
        w.postMessage({ job: J });
        feed(w);
      });
    });
    return { promise, cancel() { if (cancelled) return; cancelled = true; clearTimeout(timer); finish(); if (reject) reject(new Error('cancelled')); } };
  }

  // ------------------------------------------------------------------ Molden output
  const fE = (x, w) => { if (!Number.isFinite(x)) x = 0; const t = x.toExponential(10).replace(/e([-+])(\d)$/, 'e$10$2').toUpperCase(); return t.padStart(w); };
  /** The [GTO] and [MO] sections (lines) of the orbitals in the Molden standard (normalised primitives, components
   *  normalised on their own, Gaussian's signs, Molden's g order), for atoms Z at xyz (Å) — or null when the
   *  atoms are not those of the basis set, or the shells cannot be written in Molden's format. lines.centres holds
   *  the basis centre of each atom (Bohr; null for an atom without shells), for an exact [Atoms] AU section. */
  function moldenSections(M, Z, xyz) {
    if (!M || !Z || !xyz) return null;
    if (!M.conv) detect(M);
    const byAtom = new Map();
    for (const sh of M.shells) {
      if (!(sh.atom >= 0 && sh.atom < Z.length)) return null;
      for (let a = 0; a < 3; a++) if (Math.abs(sh.center[a] * BOHR - xyz[3 * sh.atom + a]) > 1e-4) return null;
      if (!byAtom.has(sh.atom)) byAtom.set(sh.atom, []);
      byAtom.get(sh.atom).push(sh);
    }
    // one purity per angular momentum (Molden's flags), and no Cartesian h or higher (no order defined)
    const pureL = {};
    for (const sh of M.shells) {
      if (sh.l < 2) continue;
      if (pureL[sh.l] !== undefined && pureL[sh.l] !== sh.pure) return null;
      pureL[sh.l] = sh.pure;
      if (sh.l >= 5 && !sh.pure) return null;
    }
    if (pureL[4] === false && Object.keys(pureL).some(l => l >= 5)) return null;
    const L = ['[GTO]'];
    for (const [atom, list] of Array.from(byAtom).sort((a, b) => a[0] - b[0])) {
      L.push(String(atom + 1).padStart(4) + ' 0');
      for (const sh of list) {
        const W = weights(sh, M.conv.prim);
        L.push(' ' + 'spdfghik'[sh.l] + String(sh.exps.length).padStart(5) + ' 1.00');
        for (let k = 0; k < sh.exps.length; k++) {
          const a = sh.exps[k];
          L.push(fE(a, 20) + fE(W[k] / (Math.pow(2 * a / Math.PI, 0.75) * Math.pow(4 * a, sh.l / 2)), 20));
        }
      }
      L.push('');
    }
    if (pureL[2] && pureL[3] === false) L.push('[5D10F]');
    else if (pureL[2]) L.push('[5D]');
    else if (pureL[3]) L.push('[7F]');
    if (pureL[4]) L.push('[9G]');
    // the rows of each orbital in Molden's order (Cartesian g from Gaussian's order), with the convention applied
    const fac = convFactors(M, M.conv), rows = [];
    for (const sh of M.shells) {
      if (!sh.pure && sh.l === 4 && M.source !== 'molden') {
        const g = cartOrder(4, 'fchk').map(c => c.join());
        for (const c of cartOrder(4, 'molden')) rows.push(sh.off + g.indexOf(c.join()));
      } else for (let q = 0; q < sh.nf; q++) rows.push(sh.off + q);
    }
    L.push('[MO]');
    for (const name of ['alpha', 'beta']) {
      const set = M.sets[name];
      if (!set) continue;
      for (let i = 0; i < set.n; i++) {
        L.push(' Sym= ' + (set.symmetries && set.symmetries[i] ? String(set.symmetries[i]).replace(/\s+/g, '') : (i + 1) + 'a'));
        L.push(' Ene= ' + fE(set.energies[i], 20).trim());
        L.push(' Spin= ' + (name === 'beta' ? 'Beta' : 'Alpha'));
        L.push(' Occup= ' + (Number.isFinite(set.occupations[i]) ? set.occupations[i] : 0).toFixed(6));
        const base = i * M.nbf;
        rows.forEach((r, k) => L.push(String(k + 1).padStart(5) + fE(set.coeffs[base + r] * fac[r], 20)));
      }
    }
    L.centres = Z.map((z, i) => (byAtom.has(i) ? Array.from(byAtom.get(i)[0].center) : null));
    return L;
  }

  BM.orbitals = { model, detect, check, describe, moldenSections, overlap, list, frontier, hasSpin, values, gridFor, evaluate, compute, label, job, evalBlock,
    cartOrder, solidHarmonic, pureMatrix, mOrder, convFactors, workerSource };
})(typeof window !== 'undefined' ? (window.BM = window.BM || {}) : (globalThis.BM = globalThis.BM || {}));
