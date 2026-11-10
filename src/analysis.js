/* BlueMolecule — analysis of geometry series for the plot panel.
 *
 *   BM.analysis.kabsch(A, B)                 -> {rmsd, rot (3x3 rows), ca, cb}: rotation that superposes B onto A
 *                                              (Float64Array 3N each; rot is applied to B - cb, then + ca)
 *   BM.analysis.rmsd(A, B, {align})          -> RMSD in Å (align defaults to true)
 *   BM.analysis.maxDeviation(A, B, {align})  -> largest atomic distance in Å (after alignment by default)
 *   BM.analysis.jobRange(ds, job)            -> [first, last] frame indices (inclusive) of a job, or of all frames
 *   BM.analysis.convergenceSeries(ds, {job}) -> {source, isOpt, steps, series, ratios, converged}
 *        steps  [{n, frame, label}]  (frame = index into ds.frames, or null for cycles without geometry)
 *        series [{key, label, unit, values, threshold, log}]  (values parallel to steps; null = missing)
 *        ratios [{key, label, values}]  criterion / threshold for every criterion that has a threshold
 *   BM.analysis.profiles(ds, {job})          -> [{kind: 'scan'|'irc', xLabel, points: [{x, energy, frame, label}]}]
 *   BM.analysis.nmrNuclei(ds, Z)             -> [{z, label, count}]
 *   BM.analysis.nmrSticks(ds, Z, z, {ref, merge, tol}) -> [{x, sigma, count, atoms}]  (x = ref - sigma when ref != 0)
 *   BM.analysis.moLevels(ds)                 -> [{spin, label, levels: [{index, e (Eh), occ}], homo, lumo}]
 * Energies in Hartree, lengths in Å. Uses the quaternion form of the Kabsch superposition (Horn, 1987).
 */
(function (BM) {
  'use strict';
  const KCAL = 627.509474;
  const isNum = v => typeof v === 'number' && isFinite(v);

  // ------------------------------------------------------------------ superposition
  /** Eigen-decomposition of a symmetric 4x4 matrix (cyclic Jacobi); returns {values, vectors (columns)}. */
  function jacobi4(Min) {
    const a = Min.map(r => r.slice()), v = [[1, 0, 0, 0], [0, 1, 0, 0], [0, 0, 1, 0], [0, 0, 0, 1]];
    for (let sweep = 0; sweep < 60; sweep++) {
      let off = 0;
      for (let p = 0; p < 3; p++) for (let q = p + 1; q < 4; q++) off += a[p][q] * a[p][q];
      if (off < 1e-30) break;
      for (let p = 0; p < 3; p++) for (let q = p + 1; q < 4; q++) {
        if (Math.abs(a[p][q]) < 1e-300) continue;
        const th = (a[q][q] - a[p][p]) / (2 * a[p][q]);
        const t = Math.sign(th || 1) / (Math.abs(th) + Math.sqrt(th * th + 1));
        const c = 1 / Math.sqrt(t * t + 1), s = t * c;
        for (let k = 0; k < 4; k++) { const akp = a[k][p], akq = a[k][q]; a[k][p] = c * akp - s * akq; a[k][q] = s * akp + c * akq; }
        for (let k = 0; k < 4; k++) { const apk = a[p][k], aqk = a[q][k]; a[p][k] = c * apk - s * aqk; a[q][k] = s * apk + c * aqk; }
        for (let k = 0; k < 4; k++) { const vkp = v[k][p], vkq = v[k][q]; v[k][p] = c * vkp - s * vkq; v[k][q] = s * vkp + c * vkq; }
      }
    }
    return { values: [a[0][0], a[1][1], a[2][2], a[3][3]], vectors: v };
  }
  function centroid(X, n) {
    let x = 0, y = 0, z = 0;
    for (let i = 0; i < n; i++) { x += X[3 * i]; y += X[3 * i + 1]; z += X[3 * i + 2]; }
    return n ? [x / n, y / n, z / n] : [0, 0, 0];
  }
  /** Superpose B onto A: the rotation maximising sum a_i . (R b_i) for the centred sets, and the resulting RMSD. */
  function kabsch(A, B) {
    const n = Math.min(A.length, B.length) / 3 | 0;
    if (!n) return { rmsd: 0, rot: [[1, 0, 0], [0, 1, 0], [0, 0, 1]], ca: [0, 0, 0], cb: [0, 0, 0] };
    const ca = centroid(A, n), cb = centroid(B, n);
    const S = [[0, 0, 0], [0, 0, 0], [0, 0, 0]];
    let ga = 0, gb = 0;
    for (let i = 0; i < n; i++) {
      const a = [A[3 * i] - ca[0], A[3 * i + 1] - ca[1], A[3 * i + 2] - ca[2]];
      const b = [B[3 * i] - cb[0], B[3 * i + 1] - cb[1], B[3 * i + 2] - cb[2]];
      for (let p = 0; p < 3; p++) { ga += a[p] * a[p]; gb += b[p] * b[p]; for (let q = 0; q < 3; q++) S[p][q] += b[p] * a[q]; }
    }
    const [[xx, xy, xz], [yx, yy, yz], [zx, zy, zz]] = S;
    const N = [[xx + yy + zz, yz - zy, zx - xz, xy - yx],
      [yz - zy, xx - yy - zz, xy + yx, zx + xz],
      [zx - xz, xy + yx, -xx + yy - zz, yz + zy],
      [xy - yx, zx + xz, yz + zy, -xx - yy + zz]];
    const e = jacobi4(N);
    let k = 0;
    for (let t = 1; t < 4; t++) if (e.values[t] > e.values[k]) k = t;
    const lam = e.values[k];
    let [q0, q1, q2, q3] = [e.vectors[0][k], e.vectors[1][k], e.vectors[2][k], e.vectors[3][k]];
    const qn = Math.hypot(q0, q1, q2, q3) || 1; q0 /= qn; q1 /= qn; q2 /= qn; q3 /= qn;
    const rot = [
      [q0 * q0 + q1 * q1 - q2 * q2 - q3 * q3, 2 * (q1 * q2 - q0 * q3), 2 * (q1 * q3 + q0 * q2)],
      [2 * (q1 * q2 + q0 * q3), q0 * q0 - q1 * q1 + q2 * q2 - q3 * q3, 2 * (q2 * q3 - q0 * q1)],
      [2 * (q1 * q3 - q0 * q2), 2 * (q2 * q3 + q0 * q1), q0 * q0 - q1 * q1 - q2 * q2 + q3 * q3]];
    return { rmsd: Math.sqrt(Math.max(0, (ga + gb - 2 * lam) / n)), rot, ca, cb };
  }
  /** B superposed onto A (new Float64Array). */
  function superpose(A, B) {
    const k = kabsch(A, B), n = Math.min(A.length, B.length) / 3 | 0, out = new Float64Array(3 * n), R = k.rot;
    for (let i = 0; i < n; i++) {
      const b = [B[3 * i] - k.cb[0], B[3 * i + 1] - k.cb[1], B[3 * i + 2] - k.cb[2]];
      for (let p = 0; p < 3; p++) out[3 * i + p] = R[p][0] * b[0] + R[p][1] * b[1] + R[p][2] * b[2] + k.ca[p];
    }
    return out;
  }
  function rmsd(A, B, opts) {
    const n = Math.min(A.length, B.length) / 3 | 0;
    if (!n) return 0;
    if (!opts || opts.align !== false) return kabsch(A, B).rmsd;
    let s = 0;
    for (let i = 0; i < 3 * n; i++) { const d = A[i] - B[i]; s += d * d; }
    return Math.sqrt(s / n);
  }
  function maxDeviation(A, B, opts) {
    const n = Math.min(A.length, B.length) / 3 | 0;
    const C = !opts || opts.align !== false ? superpose(A, B) : B;
    let m = 0;
    for (let i = 0; i < n; i++) m = Math.max(m, Math.hypot(A[3 * i] - C[3 * i], A[3 * i + 1] - C[3 * i + 1], A[3 * i + 2] - C[3 * i + 2]));
    return m;
  }

  // ------------------------------------------------------------------ optimisation convergence
  function jobRange(ds, job) {
    const n = ds.frames.length;
    const j = isNum(job) && job >= 0 && ds.jobs ? ds.jobs[job] : null;
    if (!j) return [0, n - 1];
    return [Math.max(0, Math.min(n - 1, j.frameStart)), Math.max(0, Math.min(n - 1, j.frameEnd))];
  }
  const CRITERIA = [
    ['maxForce', 'Max force', 'Eh/bohr'], ['rmsForce', 'RMS force', 'Eh/bohr'],
    ['maxDisp', 'Max displacement', 'bohr'], ['rmsDisp', 'RMS displacement', 'bohr'],
    ['energyChange', 'Energy change', 'Eh'], ['gradNorm', 'Gradient norm', 'Eh/bohr']];

  function convergenceSeries(ds, opts) {
    const job = opts ? opts.job : null;
    const [a, b] = jobRange(ds, job);
    const fr = [];
    for (let i = a; i <= b; i++) { const f = ds.frames[i]; if (f && f.Z && f.xyz && (isNum(f.energy) || f.opt)) fr.push(i); }
    const nOpt = fr.filter(i => ds.frames[i].opt).length;
    const os = Array.isArray(ds.optSteps) ? ds.optSteps : null;
    const jobKinds = isNum(job) && job >= 0 && ds.jobs && ds.jobs[job] ? ds.jobs[job].kinds || [] : (ds.meta && ds.meta.jobTypes) || [];
    const optJob = jobKinds.some(k => /opt/i.test(k));
    let source = 'frames', steps;
    if (os && os.length >= 2 && os.length > nOpt && (optJob || !ds.jobs || ds.jobs.length < 2)) {
      source = 'optSteps';
      steps = os.map((s, k) => ({ n: isNum(s.step) ? s.step : k + 1, frame: k === os.length - 1 ? b : null, label: 'Cycle ' + (isNum(s.step) ? s.step : k + 1), energy: isNum(s.energy) ? s.energy : null, opt: s.opt || null, converged: false }));
    } else {
      steps = fr.map((i, k) => { const f = ds.frames[i]; return { n: k + 1, frame: i, label: f.label || 'Geometry ' + (i + 1), energy: isNum(f.energy) ? f.energy : null, opt: f.opt || null, converged: !!f.converged }; });
    }
    const series = [], ratios = [];
    const out = { source, isOpt: false, steps: steps.map(s => ({ n: s.n, frame: s.frame, label: s.label })), series, ratios, converged: -1 };
    if (steps.length < 2) return out;
    const E = steps.map(s => s.energy);
    const lastE = [...E].reverse().find(isNum);
    if (E.some(isNum)) {
      series.push({ key: 'energy', label: 'Energy', unit: 'Eh', values: E, threshold: null, log: false });
      series.push({ key: 'dE', label: 'ΔE vs final', unit: 'kcal/mol', values: E.map(e => (isNum(e) ? (e - lastE) * KCAL : null)), threshold: null, log: false });
    }
    const mopac = ds.format === 'mopac-out';
    for (const [key, label, unit] of CRITERIA) {
      const vals = steps.map(s => (s.opt && isNum(s.opt[key]) ? Math.abs(s.opt[key]) : null));
      if (!vals.some(isNum)) continue;
      let thr = null;
      for (let k = steps.length - 1; k >= 0 && thr === null; k--) { const t = steps[k].opt && steps[k].opt[key + 'Thr']; if (isNum(t) && t > 0) thr = t; }
      series.push({ key, label: key === 'energyChange' ? '|Energy change|' : label, unit: mopac && key === 'gradNorm' ? 'kcal/mol/Å' : unit, values: vals, threshold: thr, log: true });
      if (thr) ratios.push({ key, label, values: vals.map(v => (isNum(v) ? v / thr : null)) });
    }
    out.isOpt = optJob || series.some(s => s.log);
    if (source === 'frames') {
      const F = steps.map(s => ds.frames[s.frame]), ref = F[F.length - 1];
      // geometry series only for optimisations (not for long MD trajectories with energies)
      if (out.isOpt && F.length * ref.Z.length <= 2e7 && F.every(f => f.Z.length === ref.Z.length && f.xyz.length === ref.xyz.length) && ref.Z.length) {
        series.push({ key: 'rmsd', label: 'RMSD to final', unit: 'Å', values: F.map((f, k) => (k === F.length - 1 ? 0 : rmsd(ref.xyz, f.xyz))), threshold: null, log: false });
        series.push({ key: 'maxStep', label: 'Largest atomic step', unit: 'Å', values: F.map((f, k) => (k ? maxDeviation(F[k - 1].xyz, f.xyz) : null)), threshold: null, log: false });
      }
      out.converged = steps.findIndex(s => s.converged);
    } else {
      const ok = s => s.opt && Object.keys(s.opt).filter(k => /Conv$/.test(k)).length > 0 && Object.keys(s.opt).filter(k => /Conv$/.test(k)).every(k => s.opt[k] === true);
      for (let k = 0; k < steps.length; k++) if (ok(steps[k])) { out.converged = k; break; }
    }
    return out;
  }

  // ------------------------------------------------------------------ scan / IRC profiles
  function coordUnit(name) {
    const s = String(name || '').trim();
    if (/^[RB]\b|^[RB]\(/i.test(s)) return ' (Å)';
    if (/^[AD]\b|^[AD]\(|^L\(/i.test(s)) return ' (°)';
    return '';
  }
  function profiles(ds, opts) {
    const [a, b] = jobRange(ds, opts ? opts.job : null), out = [];
    const scan = [], irc = [];
    for (let i = a; i <= b; i++) {
      const f = ds.frames[i];
      if (!f || !isNum(f.energy)) continue;
      if (f.scan && isNum(f.scan.value)) scan.push({ x: f.scan.value, energy: f.energy, frame: i, label: (f.scan.coord ? f.scan.coord + ' = ' : '') + f.scan.value, coord: f.scan.coord });
      if (f.irc && isNum(f.irc.coord)) irc.push({ x: f.irc.coord, energy: f.energy, frame: i, label: f.label || 'IRC point' });
    }
    if (scan.length >= 2) {
      const c = scan[scan.length - 1].coord;
      out.push({ kind: 'scan', xLabel: (c || 'Scan coordinate') + coordUnit(c), points: scan });
    } else {
      const kinds = (ds.meta && ds.meta.jobTypes) || [];
      let scanLike = kinds.indexOf('scan') >= 0;
      if (!scanLike) for (let i = a; i <= b && !scanLike; i++) scanLike = /\bscan\b/i.test((ds.frames[i] && ds.frames[i].label) || '');
      if (scanLike) {
        const pts = [];
        for (let i = a; i <= b; i++) { const f = ds.frames[i]; if (f && f.converged && isNum(f.energy)) pts.push({ x: pts.length + 1, energy: f.energy, frame: i, label: f.label || 'Scan point ' + (pts.length + 1) }); }
        if (pts.length >= 2) out.push({ kind: 'scan', xLabel: 'Scan point', points: pts });
      }
    }
    if (irc.length >= 2) out.push({ kind: 'irc', xLabel: 'Reaction coordinate', points: irc.sort((p, q) => p.x - q.x) });
    else if (Array.isArray(ds.ircPath) && ds.ircPath.filter(p => isNum(p.energy)).length >= 2) {
      out.push({ kind: 'irc', xLabel: 'IRC step', points: ds.ircPath.filter(p => isNum(p.energy)).map((p, k) => ({ x: isNum(p.step) ? p.step : k + 1, energy: p.energy, frame: null, label: 'IRC step ' + (isNum(p.step) ? p.step : k + 1) + (p.ts ? ' (TS)' : '') })) });
    }
    return out;
  }

  // ------------------------------------------------------------------ NMR
  const ISOTOPE = { 1: '1H', 3: '7Li', 5: '11B', 6: '13C', 7: '15N', 8: '17O', 9: '19F', 11: '23Na', 13: '27Al', 14: '29Si', 15: '31P', 34: '77Se', 50: '119Sn', 78: '195Pt' };
  const sym = z => (BM.el && z > 0 ? BM.el.symbol(z) : 'X');
  const nucleusLabel = z => ISOTOPE[z] || sym(z);
  const PREF = [1, 6, 7, 9, 15];
  function nmrNuclei(ds, Z) {
    const iso = ds.nmr && ds.nmr.iso;
    if (!iso || !Z) return [];
    const cnt = new Map();
    Z.forEach((z, i) => { if (z > 0 && isNum(iso[i])) cnt.set(z, (cnt.get(z) || 0) + 1); });
    const rank = z => { const k = PREF.indexOf(z); return k >= 0 ? k : 100 + z; };
    return [...cnt.keys()].sort((p, q) => rank(p) - rank(q)).map(z => ({ z, label: nucleusLabel(z), count: cnt.get(z) }));
  }
  function nmrSticks(ds, Z, z, opts) {
    const o = opts || {}, iso = ds.nmr && ds.nmr.iso, ref = isNum(o.ref) ? o.ref : 0, tol = isNum(o.tol) ? o.tol : 0.05;
    if (!iso || !Z) return [];
    const at = [];
    Z.forEach((zz, i) => { if (zz === z && isNum(iso[i])) at.push(i); });
    at.sort((p, q) => iso[p] - iso[q]);
    const groups = [];
    for (const i of at) {
      const g = groups[groups.length - 1];
      if (o.merge !== false && g && Math.abs(iso[i] - g.sigmaLast) <= tol) { g.atoms.push(i); g.sum += iso[i]; g.sigmaLast = iso[i]; }
      else groups.push({ atoms: [i], sum: iso[i], sigmaLast: iso[i] });
    }
    return groups.map(g => { const s = g.sum / g.atoms.length; return { x: ref ? ref - s : s, sigma: s, count: g.atoms.length, atoms: g.atoms }; });
  }

  // ------------------------------------------------------------------ orbital levels
  function moLevels(ds) {
    const o = ds.orbitals, out = [];
    if (!o) return out;
    const col = (s, spin, label) => {
      if (!s || !s.energies || !s.energies.length) return;
      const nOcc = isNum(s.nOcc) ? s.nOcc : 0;
      out.push({ spin, label, levels: Array.from(s.energies, (e, i) => ({ index: i + 1, e, occ: i < nOcc })), homo: nOcc > 0 ? nOcc - 1 : -1, lumo: nOcc < s.energies.length ? nOcc : -1 });
    };
    const beta = o.beta && o.beta.energies && o.beta.energies.length;
    col(o.alpha, 'alpha', beta ? 'α' : 'MO');
    if (beta) col(o.beta, 'beta', 'β');
    return out;
  }

  BM.analysis = { kabsch, superpose, rmsd, maxDeviation, jobRange, convergenceSeries, profiles, nmrNuclei, nmrSticks, nucleusLabel, moLevels, KCAL };
})(typeof window !== 'undefined' ? (window.BM = window.BM || {}) : (globalThis.BM = globalThis.BM || {}));
