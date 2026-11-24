/* BlueMolecule — conformer search with the clean-up force field (BM.conformers).
 *
 *   BM.conformers.rotors(mol)          -> [{i, j, a, b, side, step}]: rotatable bonds of mol ({Z, xyz, bonds})
 *   BM.conformers.search(mol, opts)    -> {promise, cancel()}; the promise gives {conformers, rotors, tried, mode}
 *       conformers: [{xyz, energy (kcal/mol, clean-up force field), rel, rmsd (to the lowest)}], lowest first
 *       opts: {maxConformers 20, window 10 (kcal/mol), rmsd 0.5 (Å, heavy atoms after superposition),
 *              maxSystematic 1000 (torsion combinations tried one by one), samples 300 (random combinations
 *              otherwise), seed, onProgress(fraction, found)}
 *
 * A rotatable bond is a single bond outside rings whose two atoms both have other heavy neighbours, and which is not
 * an amide C–N bond or part of a linear chain (sp atoms). Each rotor is turned in steps of 120° (180° between two
 * planar sp2 atoms). Every combination of turns (all of them when there are few, random ones otherwise) is relaxed
 * with the clean-up force field. The structures are sorted by its energy and kept when they differ from every lower
 * one by more than the RMSD threshold (heavy atoms, after superposition), are not one of them numbered another way or
 * mirrored (the same energy and the same sorted heavy-atom distances), and lie within the energy window.
 * The clean-up force field is generic: the ranking is a guide for picking starting points, not a result.
 */
(function (BM) {
  'use strict';

  function adjacency(mol) {
    const n = mol.Z.length, nbr = Array.from({ length: n }, () => []), ord = new Map();
    for (const b of mol.bonds || []) {
      const i = b[0], j = b[1], o = b[2] === undefined ? 1 : +b[2];
      if (i === j || i < 0 || j < 0 || i >= n || j >= n) continue;
      nbr[i].push(j); nbr[j].push(i); ord.set(i < j ? i + ',' + j : j + ',' + i, o);
    }
    return { nbr, order: (i, j) => ord.get(i < j ? i + ',' + j : j + ',' + i) };
  }
  /** Atoms reachable from start without crossing the bond start–block. */
  function sideOf(nbr, start, block) {
    const seen = new Set([start]), stack = [start];
    while (stack.length) {
      const a = stack.pop();
      for (const b of nbr[a]) {
        if (a === start && b === block) continue;
        if (!seen.has(b)) { seen.add(b); stack.push(b); }
      }
    }
    return seen;
  }
  function rotors(mol) {
    const { nbr, order } = adjacency(mol), Z = mol.Z, out = [];
    const heavy = i => Z[i] > 1;
    const planar = i => nbr[i].some(k => order(i, k) >= 1.5) || (Z[i] === 7 && nbr[i].some(k => nbr[k].some(m => order(k, m) >= 1.5)));
    const linear = i => nbr[i].some(k => order(i, k) >= 2.5) || (nbr[i].length === 2 && nbr[i].every(k => order(i, k) >= 1.75));
    const amide = (i, j) => {
      const cn = (c, nn) => Z[c] === 6 && Z[nn] === 7 && nbr[c].some(k => k !== nn && (Z[k] === 8 || Z[k] === 16) && order(c, k) >= 1.75);
      return cn(i, j) || cn(j, i);
    };
    const seen = new Set();
    for (const b of mol.bonds || []) {
      let i = b[0], j = b[1];
      const o = b[2] === undefined ? 1 : +b[2];
      if (o !== 1 || !(Z[i] > 0) || !(Z[j] > 0)) continue;
      const key = i < j ? i + ',' + j : j + ',' + i;
      if (seen.has(key)) continue;
      seen.add(key);
      const hi = nbr[i].filter(k => k !== j && Z[k] > 0), hj = nbr[j].filter(k => k !== i && Z[k] > 0);
      if (!hi.some(heavy) || !hj.some(heavy)) continue;          // a methyl or OH turns nothing that matters
      if (linear(i) || linear(j) || amide(i, j)) continue;
      const side = sideOf(nbr, j, i);
      if (side.has(i)) continue;                                  // in a ring
      // move the smaller side
      let s = side;
      const other = sideOf(nbr, i, j);
      if (other.size < side.size) { s = other; [i, j] = [j, i]; }
      const a = nbr[i].filter(k => k !== j && Z[k] > 0).sort((p, q) => Z[q] - Z[p])[0];
      const bb = nbr[j].filter(k => k !== i && Z[k] > 0).sort((p, q) => Z[q] - Z[p])[0];
      out.push({ i, j, a, b: bb, side: Array.from(s), step: planar(i) && planar(j) ? 180 : 120 });
    }
    return out;
  }

  function dihedral(x, a, b, c, d) {
    const p = k => [x[3 * k], x[3 * k + 1], x[3 * k + 2]];
    const A = p(a), B = p(b), C = p(c), D = p(d);
    const b0 = [A[0] - B[0], A[1] - B[1], A[2] - B[2]], b1 = [C[0] - B[0], C[1] - B[1], C[2] - B[2]], b2 = [D[0] - C[0], D[1] - C[1], D[2] - C[2]];
    const l = Math.hypot(b1[0], b1[1], b1[2]) || 1, u = [b1[0] / l, b1[1] / l, b1[2] / l];
    const pr = v => { const t = v[0] * u[0] + v[1] * u[1] + v[2] * u[2]; return [v[0] - t * u[0], v[1] - t * u[1], v[2] - t * u[2]]; };
    const v = pr(b0), w = pr(b2), cr = [u[1] * v[2] - u[2] * v[1], u[2] * v[0] - u[0] * v[2], u[0] * v[1] - u[1] * v[0]];
    return Math.atan2(cr[0] * w[0] + cr[1] * w[1] + cr[2] * w[2], v[0] * w[0] + v[1] * w[1] + v[2] * w[2]);
  }
  /** Turn the atoms of `side` about the bond i–j by ang (radians), in place. */
  function turn(x, i, j, side, ang) {
    const o = [x[3 * i], x[3 * i + 1], x[3 * i + 2]], ax = [x[3 * j] - o[0], x[3 * j + 1] - o[1], x[3 * j + 2] - o[2]];
    const l = Math.hypot(ax[0], ax[1], ax[2]);
    if (!(l > 1e-8)) return;
    const k = [ax[0] / l, ax[1] / l, ax[2] / l], c = Math.cos(ang), s = Math.sin(ang);
    for (const m of side) {
      if (m === i || m === j) continue;
      const v = [x[3 * m] - o[0], x[3 * m + 1] - o[1], x[3 * m + 2] - o[2]];
      const kv = k[0] * v[0] + k[1] * v[1] + k[2] * v[2], cx = [k[1] * v[2] - k[2] * v[1], k[2] * v[0] - k[0] * v[2], k[0] * v[1] - k[1] * v[0]];
      for (let q = 0; q < 3; q++) x[3 * m + q] = o[q] + v[q] * c + cx[q] * s + k[q] * kv * (1 - c);
    }
  }
  function rng(seed) { let s = (seed >>> 0) || 0x9e3779b9; return () => { s ^= s << 13; s >>>= 0; s ^= s >>> 17; s ^= s << 5; s >>>= 0; return s / 4294967296; }; }
  /** Heavy-atom RMSD after superposition. */
  function heavyRmsd(Z, A, B) {
    const idx = []; for (let i = 0; i < Z.length; i++) if (Z[i] > 1) idx.push(i);
    const take = X => { const o = new Float64Array(3 * idx.length); idx.forEach((i, k) => { o[3 * k] = X[3 * i]; o[3 * k + 1] = X[3 * i + 1]; o[3 * k + 2] = X[3 * i + 2]; }); return o; };
    return idx.length ? BM.analysis.rmsd(take(A), take(B)) : 0;
  }

  /** Sorted heavy-atom distances: the same for structures that differ only by numbering or by a mirror. */
  function shape(Z, X) {
    const idx = []; for (let i = 0; i < Z.length; i++) if (Z[i] > 1) idx.push(i);
    const d = [];
    for (let p = 0; p < idx.length; p++) for (let q = p + 1; q < idx.length; q++) {
      const i = idx[p], j = idx[q];
      d.push(Math.hypot(X[3 * i] - X[3 * j], X[3 * i + 1] - X[3 * j + 1], X[3 * i + 2] - X[3 * j + 2]));
    }
    return Float64Array.from(d).sort();
  }
  function sameShape(a, b) {
    if (a.length !== b.length) return false;
    let s = 0; for (let k = 0; k < a.length; k++) { const t = a[k] - b[k]; s += t * t; }
    return Math.sqrt(s / Math.max(1, a.length)) < 0.03;
  }
  function search(mol, opts) {
    opts = opts || {};
    const maxConf = opts.maxConformers > 0 ? opts.maxConformers : 20, win = opts.window > 0 ? opts.window : 10;
    const thr = opts.rmsd > 0 ? opts.rmsd : 0.5, maxSys = opts.maxSystematic > 0 ? opts.maxSystematic : 1000, samples = opts.samples > 0 ? opts.samples : 300;
    const n = mol.Z.length;
    let cancelled = false, timer = 0, reject = null;
    const promise = new Promise((resolve, rej) => {
      reject = rej;
      if (!n) { rej(new Error('There are no atoms.')); return; }
      const R = rotors(mol);
      const steps = R.map(r => Math.round(360 / r.step));
      let total = 1;
      for (const s of steps) { total *= s; if (total > 1e9) break; }
      const systematic = total <= maxSys, count = systematic ? total : samples;
      const rand = rng(opts.seed || 12345);
      const relax = x => BM.ff.optimize({ Z: mol.Z, xyz: x, bonds: mol.bonds, formalCharges: mol.formalCharges }, { maxSteps: opts.maxSteps || 600, tol: 0.05 });
      const start = relax(Float64Array.from(mol.xyz));
      const found = [];   // {xyz, energy}
      const consider = r => { if (Number.isFinite(r.energy)) found.push({ xyz: r.xyz, energy: r.energy }); };
      consider(start);
      const dih0 = R.map(r => dihedral(start.xyz, r.a, r.i, r.j, r.b));
      let k = 0;
      const one = idx => {
        // the combination: each rotor's turn (multiples of its step)
        const x = Float64Array.from(start.xyz);
        let rest = idx, changed = false;
        R.forEach((r, q) => {
          const m = systematic ? rest % steps[q] : Math.floor(rand() * steps[q]);
          if (systematic) rest = Math.floor(rest / steps[q]);
          if (!m) return;
          const target = dih0[q] + m * r.step * Math.PI / 180, now = dihedral(x, r.a, r.i, r.j, r.b);
          turn(x, r.i, r.j, r.side, target - now);
          changed = true;
        });
        if (changed) consider(relax(x));
      };
      const finish = () => {
        found.sort((p, q) => p.energy - q.energy);
        const kept = [];
        for (const c of found) {
          if (kept.length >= maxConf) break;
          if (c.energy - found[0].energy > win) break;
          if (kept.some(q => heavyRmsd(mol.Z, q.xyz, c.xyz) < thr)) continue;
          // the same conformer numbered another way, or its mirror image (same energy, same distances)
          c.shape = shape(mol.Z, c.xyz);
          if (kept.some(q => Math.abs(q.energy - c.energy) < 0.05 && sameShape(q.shape, c.shape))) continue;
          kept.push(c);
        }
        const e0 = kept.length ? kept[0].energy : 0;
        resolve({ conformers: kept.map(c => ({ xyz: c.xyz, energy: c.energy, rel: c.energy - e0, rmsd: heavyRmsd(mol.Z, kept[0].xyz, c.xyz) })),
          rotors: R, tried: found.length, mode: systematic ? 'systematic' : 'random', combinations: total });
      };
      const run = () => {
        if (cancelled) return;
        const t0 = Date.now();
        try {
          while (k < count && Date.now() - t0 < 40) one(k++);
        } catch (e) { rej(e); return; }
        if (opts.onProgress) try { opts.onProgress(count ? k / count : 1, found.length); } catch (e) { /* ignore */ }
        if (k >= count) finish(); else timer = setTimeout(run, 0);
      };
      timer = setTimeout(run, 0);
    });
    return { promise, cancel() { if (cancelled) return; cancelled = true; clearTimeout(timer); if (reject) reject(new Error('cancelled')); } };
  }

  BM.conformers = { rotors, search, dihedral, turn, heavyRmsd };
})(typeof window !== 'undefined' ? (window.BM = window.BM || {}) : (globalThis.BM = globalThis.BM || {}));
