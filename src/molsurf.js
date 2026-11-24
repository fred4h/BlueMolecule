/* BlueMolecule — molecular surfaces, and the electrostatic potential mapped onto surfaces.
 *
 *   BM.molsurf.field(Z, xyz, opts) -> volume whose zero level is the surface (positive inside), in Å:
 *       opts.type 'vdw' (van der Waals: the union of the atoms' spheres), 'sas' (solvent accessible: the spheres grown
 *       by the probe radius, traced by the probe's centre), 'ses' (solvent excluded, or Connolly: the surface the
 *       probe's front touches as it rolls over the atoms); opts.probe (Å, 1.4 = water), opts.spacing (Å), opts.maxPoints
 *   BM.molsurf.espCharges(Z, xyz, q, pts) -> Float64Array: the potential of point charges q at pts (Å), atomic units
 *   BM.molsurf.sample(vol, pts)          -> Float64Array: a volume interpolated (trilinear) at pts (Å), NaN outside
 *   BM.molsurf.colors(values, lo, hi)     -> Float32Array RGB: red (negative), white, blue (positive); grey for NaN
 *
 * The solvent-excluded surface: the probe centres that touch no atom form the region A outside the solvent-
 * accessible surface; the excluded region is every point farther than the probe radius from A. The distance to A is
 * measured to points of A's boundary (the grid points next to it, projected onto the accessible sphere they belong
 * to), so the surface is accurate to a small fraction of the grid spacing.
 */
(function (BM) {
  'use strict';
  const BOHR = 0.529177210903;

  function radii(Z, add) {
    const v = (BM.el && BM.el.vdw) || [];
    return Z.map(z => (z > 0 ? (v[z] || 1.8) + add : 0));
  }
  /** The surface's field (self-contained, so that it also runs in a worker): xyz (Å), Rs (grown radii, 0 = no
   *  sphere), type, probe, spacing h, maxPoints -> {origin, h, dims, data}. */
  function fieldCore(xyz, Rs, type, probe, h, maxPoints) {
    /** Grid of spacing h (widened to stay under maxPoints) around spheres of radii R plus pad. */
    function boxFor(xyz, R, h, pad, maxPoints) {
      const lo = [Infinity, Infinity, Infinity], hi = [-Infinity, -Infinity, -Infinity];
      for (let i = 0; i < R.length; i++) {
        if (!(R[i] > 0)) continue;
        for (let a = 0; a < 3; a++) { lo[a] = Math.min(lo[a], xyz[3 * i + a] - R[i] - pad); hi[a] = Math.max(hi[a], xyz[3 * i + a] + R[i] + pad); }
      }
      if (!(lo[0] < hi[0])) return null;
      const ext = [0, 1, 2].map(a => hi[a] - lo[a]);
      let dims = ext.map(e => Math.max(2, Math.ceil(e / h) + 1));
      while (dims[0] * dims[1] * dims[2] > maxPoints) { h *= 1.06; dims = ext.map(e => Math.max(2, Math.ceil(e / h) + 1)); }
      return { origin: lo, h, dims };
    }

    /** Spatial hash of points (flat xyz) in cubic cells of size c. */
    function hash(pts, c) {
      const map = new Map(), key = (i, j, k) => ((i * 73856093) ^ (j * 19349663) ^ (k * 83492791));
      for (let p = 0; p < pts.length / 3; p++) {
        const i = Math.floor(pts[3 * p] / c), j = Math.floor(pts[3 * p + 1] / c), k = Math.floor(pts[3 * p + 2] / c), kk = key(i, j, k);
        let b = map.get(kk);
        if (!b) { b = []; map.set(kk, b); }
        b.push(p);
      }
      // every point within s cells (s * c) of (x, y, z), and perhaps a few more
      return { c, near(x, y, z, fn, s) {
        s = s || 1;
        const i = Math.floor(x / c), j = Math.floor(y / c), k = Math.floor(z / c);
        for (let a = i - s; a <= i + s; a++) for (let b = j - s; b <= j + s; b++) for (let d = k - s; d <= k + s; d++) {
          const l = map.get(key(a, b, d));
          if (l) for (let q = 0; q < l.length; q++) fn(l[q]);
        }
      } };
    }
    const n = Rs.length, maxR = Math.max.apply(null, Rs);
    const box = boxFor(xyz, Rs, h, 3 * h, maxPoints);
    if (!box) throw new Error('There are no atoms with a radius.');
    h = box.h;
    const [n1, n2, n3] = box.dims, o = box.origin, N = n1 * n2 * n3;
    // D: signed distance to the nearest (grown) sphere, min over atoms; near: that atom
    const D = new Float32Array(N).fill(1e9), near = type === 'ses' ? new Int32Array(N).fill(-1) : null;
    for (let i = 0; i < n; i++) {
      const R = Rs[i];
      if (!(R > 0)) continue;
      const reach = R + 3 * h, cx = xyz[3 * i], cy = xyz[3 * i + 1], cz = xyz[3 * i + 2];
      const i0 = Math.max(0, Math.floor((cx - reach - o[0]) / h)), i1 = Math.min(n1 - 1, Math.ceil((cx + reach - o[0]) / h));
      const j0 = Math.max(0, Math.floor((cy - reach - o[1]) / h)), j1 = Math.min(n2 - 1, Math.ceil((cy + reach - o[1]) / h));
      const k0 = Math.max(0, Math.floor((cz - reach - o[2]) / h)), k1 = Math.min(n3 - 1, Math.ceil((cz + reach - o[2]) / h));
      for (let a = i0; a <= i1; a++) {
        const dx = o[0] + a * h - cx;
        for (let b = j0; b <= j1; b++) {
          const dy = o[1] + b * h - cy, dxy = dx * dx + dy * dy, row = (a * n2 + b) * n3;
          if (dxy > reach * reach) continue;
          for (let c = k0; c <= k1; c++) {
            const dz = o[2] + c * h - cz, d = Math.sqrt(dxy + dz * dz) - R;
            if (d < D[row + c]) { D[row + c] = d; if (near) near[row + c] = i; }
          }
        }
      }
    }
    const data = new Float32Array(N);
    if (type !== 'ses') {
      const cap = 3 * h;
      for (let p = 0; p < N; p++) { const v = -D[p]; data[p] = v < -cap ? -cap : v; }
    } else {
      // points of A's boundary: grid points outside every grown sphere next to one inside, projected onto the sphere
      const seeds = [];
      const atomHash = hash(xyz, Math.max(maxR, 1));
      // atoms in buckets of maxR over the box (compressed rows): is a point outside every grown sphere?
      const ca = Math.max(maxR, 1), ad = [0, 1, 2].map(q => Math.max(1, Math.ceil((box.dims[q] - 1) * h / ca) + 1));
      const aCell = (x, q) => Math.max(0, Math.min(ad[q] - 1, Math.floor((x - o[q]) / ca)));
      const aStart = new Int32Array(ad[0] * ad[1] * ad[2] + 1), aIdx = new Int32Array(n), aOf = new Int32Array(n);
      for (let i = 0; i < n; i++) { aOf[i] = (aCell(xyz[3 * i], 0) * ad[1] + aCell(xyz[3 * i + 1], 1)) * ad[2] + aCell(xyz[3 * i + 2], 2); aStart[aOf[i] + 1]++; }
      for (let q = 0; q < aStart.length - 1; q++) aStart[q + 1] += aStart[q];
      { const fill = aStart.slice(0, -1); for (let i = 0; i < n; i++) aIdx[fill[aOf[i]]++] = i; }
      const R2in = Float64Array.from(Rs, r => (r > 0 ? (r - 1e-4) * (r - 1e-4) : -1));
      const inA = (x, y, z) => {
        const ci = aCell(x, 0), cj = aCell(y, 1), ck = aCell(z, 2);
        for (let u = Math.max(0, ci - 1); u <= Math.min(ad[0] - 1, ci + 1); u++) {
          for (let v = Math.max(0, cj - 1); v <= Math.min(ad[1] - 1, cj + 1); v++) {
            const row = (u * ad[1] + v) * ad[2], e1 = aStart[row + Math.min(ad[2] - 1, ck + 1) + 1];
            for (let e = aStart[row + Math.max(0, ck - 1)]; e < e1; e++) {
              const j = aIdx[e], ex = x - xyz[3 * j], ey = y - xyz[3 * j + 1], ez = z - xyz[3 * j + 2];
              if (ex * ex + ey * ey + ez * ez < R2in[j]) return false;
            }
          }
        }
        return true;
      };
      // atoms whose grown sphere comes within two grid steps of A's boundary; the others are buried
      const surfaceAtom = new Uint8Array(n), reachA = Float64Array.from(Rs, r => (r + 2 * h) * (r + 2 * h));
      const markNear = (x, y, z) => {
        const ci = aCell(x, 0), cj = aCell(y, 1), ck = aCell(z, 2);
        for (let u = Math.max(0, ci - 1); u <= Math.min(ad[0] - 1, ci + 1); u++) {
          for (let v = Math.max(0, cj - 1); v <= Math.min(ad[1] - 1, cj + 1); v++) {
            const row = (u * ad[1] + v) * ad[2], e1 = aStart[row + Math.min(ad[2] - 1, ck + 1) + 1];
            for (let e = aStart[row + Math.max(0, ck - 1)]; e < e1; e++) {
              const j = aIdx[e], ex = x - xyz[3 * j], ey = y - xyz[3 * j + 1], ez = z - xyz[3 * j + 2];
              if (Rs[j] > 0 && ex * ex + ey * ey + ez * ez < reachA[j]) surfaceAtom[j] = 1;
            }
          }
        }
      };
      const st = [n2 * n3, n3, 1];
      for (let a = 0; a < n1; a++) for (let b = 0; b < n2; b++) for (let c = 0; c < n3; c++) {
        const p = (a * n2 + b) * n3 + c;
        if (D[p] < 0) continue;
        const nb = (a > 0 && D[p - st[0]] < 0) || (a < n1 - 1 && D[p + st[0]] < 0) || (b > 0 && D[p - st[1]] < 0) || (b < n2 - 1 && D[p + st[1]] < 0) ||
          (c > 0 && D[p - 1] < 0) || (c < n3 - 1 && D[p + 1] < 0);
        if (!nb) continue;
        const x = o[0] + a * h, y = o[1] + b * h, z = o[2] + c * h, i = near[p];
        markNear(x, y, z);
        if (i >= 0) {
          const vx = x - xyz[3 * i], vy = y - xyz[3 * i + 1], vz = z - xyz[3 * i + 2], l = Math.hypot(vx, vy, vz) || 1, f = Rs[i] / l;
          const px = xyz[3 * i] + vx * f, py = xyz[3 * i + 1] + vy * f, pz = xyz[3 * i + 2] + vz * f;
          if (inA(px, py, pz)) { seeds.push(px, py, pz); continue; }
        }
        seeds.push(x, y, z);
      }
      // A's sharp edges, where the probe touches two atoms (circles where two grown spheres meet) or three (points):
      // the nearest points of A to the re-entrant parts of the surface
      const nbrs = [];
      for (let i = 0; i < n; i++) {
        const l = [];
        // (circles and points of buried atoms are buried too)
        if (Rs[i] > 0 && surfaceAtom[i]) atomHash.near(xyz[3 * i], xyz[3 * i + 1], xyz[3 * i + 2], j => {
          if (j <= i || !(Rs[j] > 0) || !surfaceAtom[j]) return;
          const d = Math.hypot(xyz[3 * j] - xyz[3 * i], xyz[3 * j + 1] - xyz[3 * i + 1], xyz[3 * j + 2] - xyz[3 * i + 2]);
          if (d < Rs[i] + Rs[j] && d > Math.abs(Rs[i] - Rs[j])) l.push(j);
        }, 2);
        nbrs.push(l);
      }
      const C = i => [xyz[3 * i], xyz[3 * i + 1], xyz[3 * i + 2]];
      const sub = (a, b) => [a[0] - b[0], a[1] - b[1], a[2] - b[2]], dot = (a, b) => a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
      const cross = (a, b) => [a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]];
      const unit = a => { const l = Math.hypot(a[0], a[1], a[2]) || 1; return [a[0] / l, a[1] / l, a[2] / l]; };
      for (let i = 0; i < n; i++) {
        for (const j of nbrs[i]) {
          const ci = C(i), cj = C(j), v = sub(cj, ci), d = Math.hypot(v[0], v[1], v[2]), u = [v[0] / d, v[1] / d, v[2] / d];
          const a = (d * d + Rs[i] * Rs[i] - Rs[j] * Rs[j]) / (2 * d), rho = Math.sqrt(Math.max(0, Rs[i] * Rs[i] - a * a));
          const e1 = unit(Math.abs(u[0]) < 0.9 ? cross(u, [1, 0, 0]) : cross(u, [0, 1, 0])), e2 = cross(u, e1);
          const m = Math.max(8, Math.ceil(2 * Math.PI * rho / h));
          for (let t = 0; t < m; t++) {
            const th = 2 * Math.PI * t / m, cs = Math.cos(th) * rho, sn = Math.sin(th) * rho;
            const px = ci[0] + a * u[0] + cs * e1[0] + sn * e2[0], py = ci[1] + a * u[1] + cs * e1[1] + sn * e2[1], pz = ci[2] + a * u[2] + cs * e1[2] + sn * e2[2];
            if (inA(px, py, pz)) seeds.push(px, py, pz);
          }
          // three spheres: i, j and a common neighbour k > j
          for (const k of nbrs[i]) {
            if (k <= j) continue;
            const ck = C(k), dk = Math.hypot(ck[0] - cj[0], ck[1] - cj[1], ck[2] - cj[2]);
            if (!(dk < Rs[j] + Rs[k] && dk > Math.abs(Rs[j] - Rs[k]))) continue;
            const w = sub(ck, ci), ii = dot(u, w), ey0 = [w[0] - ii * u[0], w[1] - ii * u[1], w[2] - ii * u[2]], jj = Math.hypot(ey0[0], ey0[1], ey0[2]);
            if (jj < 1e-8) continue;
            const ey = [ey0[0] / jj, ey0[1] / jj, ey0[2] / jj], ez = cross(u, ey);
            const x = (Rs[i] * Rs[i] - Rs[j] * Rs[j] + d * d) / (2 * d), y = (Rs[i] * Rs[i] - Rs[k] * Rs[k] + ii * ii + jj * jj) / (2 * jj) - ii / jj * x;
            const z2 = Rs[i] * Rs[i] - x * x - y * y;
            if (z2 < 0) continue;
            for (const s of [1, -1]) {
              const z = s * Math.sqrt(z2);
              const px = ci[0] + x * u[0] + y * ey[0] + z * ez[0], py = ci[1] + x * u[1] + y * ey[1] + z * ez[1], pz = ci[2] + x * u[2] + y * ey[2] + z * ez[2];
              if (inA(px, py, pz)) seeds.push(px, py, pz);
            }
          }
        }
      }
      // nearest seed of each point within reach of A: seeds in buckets of reach / 2 over the box (compressed rows)
      const reach = probe + 2 * h, cs = reach / 2, ns = seeds.length / 3;
      const bd = [0, 1, 2].map(q => Math.max(1, Math.ceil((box.dims[q] - 1) * h / cs) + 1));
      const cellOf = (x, q) => Math.max(0, Math.min(bd[q] - 1, Math.floor((x - o[q]) / cs)));
      const start = new Int32Array(bd[0] * bd[1] * bd[2] + 1), idx = new Int32Array(ns), cell = new Int32Array(ns);
      for (let q = 0; q < ns; q++) { cell[q] = (cellOf(seeds[3 * q], 0) * bd[1] + cellOf(seeds[3 * q + 1], 1)) * bd[2] + cellOf(seeds[3 * q + 2], 2); start[cell[q] + 1]++; }
      for (let q = 0; q < start.length - 1; q++) start[q + 1] += start[q];
      const fill = start.slice(0, -1);
      for (let q = 0; q < ns; q++) idx[fill[cell[q]]++] = q;
      const S = Float64Array.from(seeds);
      for (let a = 0; a < n1; a++) for (let b = 0; b < n2; b++) for (let c = 0; c < n3; c++) {
        const p = (a * n2 + b) * n3 + c, d = D[p];
        if (d >= 0) { data[p] = -probe; continue; }
        if (d < -reach) { data[p] = probe; continue; }   // deeper than the probe inside a grown sphere: excluded
        const x = o[0] + a * h, y = o[1] + b * h, z = o[2] + c * h;
        const ci = cellOf(x, 0), cj = cellOf(y, 1), ck = cellOf(z, 2);
        let best = reach * reach;
        for (let u = Math.max(0, ci - 2); u <= Math.min(bd[0] - 1, ci + 2); u++) {
          for (let v = Math.max(0, cj - 2); v <= Math.min(bd[1] - 1, cj + 2); v++) {
            const row = (u * bd[1] + v) * bd[2];
            const e0 = start[row + Math.max(0, ck - 2)], e1 = start[row + Math.min(bd[2] - 1, ck + 2) + 1];
            for (let e = e0; e < e1; e++) {
              const q = 3 * idx[e], ex = x - S[q], ey = y - S[q + 1], ez = z - S[q + 2], r2 = ex * ex + ey * ey + ez * ez;
              if (r2 < best) best = r2;
            }
          }
        }
        const v = Math.sqrt(best) - probe;
        data[p] = v > probe ? probe : v;
      }
    }
    return { origin: o.slice(), h, dims: [n1, n2, n3], data };
  }
  function args(Z, opts) {
    opts = opts || {};
    const type = opts.type === 'sas' || opts.type === 'ses' ? opts.type : 'vdw';
    const probe = type === 'vdw' ? 0 : (opts.probe > 0 ? opts.probe : 1.4);
    if (!Z.length) throw new Error('There are no atoms.');
    const n = Z.length, h = opts.spacing > 0 ? opts.spacing : n <= 300 ? 0.25 : n <= 3000 ? 0.4 : 0.6;
    return { type, probe, Rs: radii(Z, probe), h, maxPoints: opts.maxPoints > 0 ? opts.maxPoints : 8e6 };
  }
  function volume(r, a) {
    const name = a.type === 'vdw' ? 'van der Waals surface' : (a.type === 'sas' ? 'Solvent-accessible surface' : 'Solvent-excluded surface') + ' (probe ' + a.probe.toFixed(2) + ' Å)';
    return { name, kind: 'molsurf', origin: r.origin, axes: [[r.h, 0, 0], [0, r.h, 0], [0, 0, r.h]], dims: r.dims, data: r.data, iso: 0,
      surface: { type: a.type, probe: a.probe, spacing: r.h } };
  }
  /** The surface at once (tests, small molecules). */
  function field(Z, xyz, opts) {
    const a = args(Z, opts);
    return volume(fieldCore(Array.from(xyz), a.Rs, a.type, a.probe, a.h, a.maxPoints), a);
  }
  let workerURL = null;
  /** The surface in a worker when the page can start one: {promise (the volume), cancel()}. */
  function fieldAsync(Z, xyz, opts) {
    const a = args(Z, opts), X = Array.from(xyz);
    let w = null, timer = 0, cancelled = false, rej = null;
    const promise = new Promise((resolve, reject) => {
      rej = reject;
      const local = () => { if (cancelled) return; try { resolve(volume(fieldCore(X, a.Rs, a.type, a.probe, a.h, a.maxPoints), a)); } catch (e) { reject(e); } };
      try {
        if (typeof Worker !== 'undefined' && typeof Blob !== 'undefined' && typeof URL !== 'undefined' && URL.createObjectURL) {
          if (!workerURL) workerURL = URL.createObjectURL(new Blob(['"use strict";\nconst fieldCore = ' + fieldCore.toString() + ';\n' +
            'self.onmessage = e => { const m = e.data; self.postMessage({ ready: true }); try { const r = fieldCore(m.xyz, m.Rs, m.type, m.probe, m.h, m.maxPoints); self.postMessage(r, [r.data.buffer]); } catch (err) { self.postMessage({ error: String(err && err.message || err) }); } };\n'], { type: 'text/javascript' }));
          w = new Worker(workerURL);
        }
      } catch (e) { w = null; }
      if (w) {
        // a worker that cannot run (a page that forbids them, or a web view where they never start): the same on the page
        let alive = false;
        const giveUp = () => { if (w) w.terminate(); w = null; if (!cancelled) timer = setTimeout(local, 0); };
        const watchdog = setTimeout(() => { if (!alive) giveUp(); }, opts && opts.watchdogMs > 0 ? opts.watchdogMs : 4000);
        w.onmessage = e => {
          const r = e.data;
          if (r.ready) { alive = true; clearTimeout(watchdog); return; }
          if (w) w.terminate(); w = null;
          if (cancelled) return;
          if (r.error) reject(new Error(r.error)); else resolve(volume(r, a));
        };
        w.onerror = e => { if (e && e.preventDefault) e.preventDefault(); clearTimeout(watchdog); giveUp(); };
        w.postMessage({ xyz: X, Rs: a.Rs, type: a.type, probe: a.probe, h: a.h, maxPoints: a.maxPoints });
      } else timer = setTimeout(local, 0);
    });
    return { promise, cancel() { if (cancelled) return; cancelled = true; clearTimeout(timer); if (w) { w.terminate(); w = null; } if (rej) rej(new Error('cancelled')); } };
  }

  /** Electrostatic potential of point charges (e) at points (Å), in atomic units (Hartree per e). */
  function espCharges(Z, xyz, q, pts) {
    const n = Z.length, m = pts.length / 3, out = new Float64Array(m);
    for (let p = 0; p < m; p++) {
      const x = pts[3 * p], y = pts[3 * p + 1], z = pts[3 * p + 2];
      let v = 0;
      for (let i = 0; i < n; i++) {
        const qi = q[i];
        if (!qi) continue;
        const dx = x - xyz[3 * i], dy = y - xyz[3 * i + 1], dz = z - xyz[3 * i + 2];
        v += qi / Math.max(1e-6, Math.sqrt(dx * dx + dy * dy + dz * dz));
      }
      out[p] = v * BOHR;   // e / Å -> e / Bohr
    }
    return out;
  }

  function inverse3(a1, a2, a3) {
    const m = [a1[0], a2[0], a3[0], a1[1], a2[1], a3[1], a1[2], a2[2], a3[2]];
    const det = m[0] * (m[4] * m[8] - m[5] * m[7]) - m[1] * (m[3] * m[8] - m[5] * m[6]) + m[2] * (m[3] * m[7] - m[4] * m[6]);
    if (!det) return null;
    const id = 1 / det;
    return [(m[4] * m[8] - m[5] * m[7]) * id, (m[2] * m[7] - m[1] * m[8]) * id, (m[1] * m[5] - m[2] * m[4]) * id,
      (m[5] * m[6] - m[3] * m[8]) * id, (m[0] * m[8] - m[2] * m[6]) * id, (m[2] * m[3] - m[0] * m[5]) * id,
      (m[3] * m[7] - m[4] * m[6]) * id, (m[1] * m[6] - m[0] * m[7]) * id, (m[0] * m[4] - m[1] * m[3]) * id];
  }
  /** A volume interpolated at points (Å): trilinear, NaN outside the grid. */
  function sample(vol, pts) {
    const [n1, n2, n3] = vol.dims, d = vol.data, o = vol.origin, M = inverse3(vol.axes[0], vol.axes[1], vol.axes[2]);
    const m = pts.length / 3, out = new Float64Array(m);
    if (!M) { out.fill(NaN); return out; }
    for (let p = 0; p < m; p++) {
      const x = pts[3 * p] - o[0], y = pts[3 * p + 1] - o[1], z = pts[3 * p + 2] - o[2];
      const u = M[0] * x + M[1] * y + M[2] * z, v = M[3] * x + M[4] * y + M[5] * z, w = M[6] * x + M[7] * y + M[8] * z;
      if (!(u >= 0 && v >= 0 && w >= 0 && u <= n1 - 1 && v <= n2 - 1 && w <= n3 - 1)) { out[p] = NaN; continue; }
      const i = Math.min(n1 - 2, Math.floor(u)), j = Math.min(n2 - 2, Math.floor(v)), k = Math.min(n3 - 2, Math.floor(w));
      const fu = u - i, fv = v - j, fw = w - k, at = (a, b, c) => d[((i + a) * n2 + j + b) * n3 + k + c];
      if (n1 < 2 || n2 < 2 || n3 < 2) { out[p] = d[0]; continue; }
      out[p] = (1 - fu) * ((1 - fv) * ((1 - fw) * at(0, 0, 0) + fw * at(0, 0, 1)) + fv * ((1 - fw) * at(0, 1, 0) + fw * at(0, 1, 1))) +
        fu * ((1 - fv) * ((1 - fw) * at(1, 0, 0) + fw * at(1, 0, 1)) + fv * ((1 - fw) * at(1, 1, 0) + fw * at(1, 1, 1)));
    }
    return out;
  }

  /** Diverging colours: red at lo (negative potential, electron rich), white at 0, blue at hi; grey for NaN. */
  function colors(values, lo, hi) {
    const out = new Float32Array(values.length * 3);
    const neg = [0.84, 0.19, 0.15], mid = [0.97, 0.97, 0.97], pos = [0.17, 0.35, 0.86];
    for (let p = 0; p < values.length; p++) {
      const v = values[p];
      let c;
      if (!Number.isFinite(v)) c = [0.6, 0.6, 0.6];
      else if (v < 0) { const t = lo < 0 ? Math.min(1, v / lo) : 1; c = [0, 1, 2].map(a => mid[a] + (neg[a] - mid[a]) * t); }
      else { const t = hi > 0 ? Math.min(1, v / hi) : 1; c = [0, 1, 2].map(a => mid[a] + (pos[a] - mid[a]) * t); }
      out[3 * p] = c[0]; out[3 * p + 1] = c[1]; out[3 * p + 2] = c[2];
    }
    return out;
  }

  BM.molsurf = { field, fieldAsync, fieldCore, espCharges, sample, colors };
})(typeof window !== 'undefined' ? (window.BM = window.BM || {}) : (globalThis.BM = globalThis.BM || {}));
