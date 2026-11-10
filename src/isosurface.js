/* BlueMolecule — isosurface extraction (surface nets) for volumetric data (cube files). */
(function (BM) {
  'use strict';

  const EDGES = [[0, 1], [2, 3], [4, 5], [6, 7], [0, 2], [1, 3], [4, 6], [5, 7], [0, 4], [1, 5], [2, 6], [3, 7]];
  const CORNER = [[0, 0, 0], [1, 0, 0], [0, 1, 0], [1, 1, 0], [0, 0, 1], [1, 0, 1], [0, 1, 1], [1, 1, 1]];

  function inv3T(a1, a2, a3) {
    // M = [a1 a2 a3] (columns). Return inverse-transpose as row-major 3x3.
    const m = [a1[0], a2[0], a3[0], a1[1], a2[1], a3[1], a1[2], a2[2], a3[2]];
    const det = m[0] * (m[4] * m[8] - m[5] * m[7]) - m[1] * (m[3] * m[8] - m[5] * m[6]) + m[2] * (m[3] * m[7] - m[4] * m[6]);
    const id = 1 / det;
    const inv = [
      (m[4] * m[8] - m[5] * m[7]) * id, (m[2] * m[7] - m[1] * m[8]) * id, (m[1] * m[5] - m[2] * m[4]) * id,
      (m[5] * m[6] - m[3] * m[8]) * id, (m[0] * m[8] - m[2] * m[6]) * id, (m[2] * m[3] - m[0] * m[5]) * id,
      (m[3] * m[7] - m[4] * m[6]) * id, (m[1] * m[6] - m[0] * m[7]) * id, (m[0] * m[4] - m[1] * m[3]) * id];
    return [inv[0], inv[3], inv[6], inv[1], inv[4], inv[7], inv[2], inv[5], inv[8]];
  }

  /**
   * Extract the surface where value == iso. For iso < 0 the surface encloses values below iso.
   * Returns { positions: Float32Array, normals: Float32Array, indices: Uint32Array }.
   */
  function extract(vol, iso) {
    const [n1, n2, n3] = vol.dims, data = vol.data;
    const s = iso >= 0 ? 1 : -1, lvl = Math.abs(iso);
    const o = vol.origin, a1 = vol.axes[0], a2 = vol.axes[1], a3 = vol.axes[2];
    const G = inv3T(a1, a2, a3);
    const stride1 = n2 * n3, stride2 = n3;
    const g = (i, j, k) => s * data[i * stride1 + j * stride2 + k] - lvl;
    const gc = (i, j, k) => {
      i = i < 0 ? 0 : i >= n1 ? n1 - 1 : i; j = j < 0 ? 0 : j >= n2 ? n2 - 1 : j; k = k < 0 ? 0 : k >= n3 ? n3 - 1 : k;
      return s * data[i * stride1 + j * stride2 + k] - lvl;
    };
    const c1 = n1 - 1, c2 = n2 - 1, c3 = n3 - 1;
    if (c1 < 1 || c2 < 1 || c3 < 1) return { positions: new Float32Array(0), normals: new Float32Array(0), indices: new Uint32Array(0) };
    const cellIdx = new Int32Array(c1 * c2 * c3).fill(-1);
    const pos = [], nrm = [];
    const val = new Float64Array(8);
    const grad = new Float64Array(24);
    for (let i = 0; i < c1; i++) for (let j = 0; j < c2; j++) for (let k = 0; k < c3; k++) {
      let mask = 0;
      for (let c = 0; c < 8; c++) {
        const v = g(i + CORNER[c][0], j + CORNER[c][1], k + CORNER[c][2]);
        val[c] = v; if (v > 0) mask |= 1 << c;
      }
      if (mask === 0 || mask === 255) continue;
      let px = 0, py = 0, pz = 0, cnt = 0;
      for (const [a, b] of EDGES) {
        const va = val[a], vb = val[b];
        if ((va > 0) === (vb > 0)) continue;
        const t = va / (va - vb);
        px += CORNER[a][0] + t * (CORNER[b][0] - CORNER[a][0]);
        py += CORNER[a][1] + t * (CORNER[b][1] - CORNER[a][1]);
        pz += CORNER[a][2] + t * (CORNER[b][2] - CORNER[a][2]);
        cnt++;
      }
      const u = px / cnt, v = py / cnt, w = pz / cnt;
      // central-difference gradients at the 8 corners, blended trilinearly
      for (let c = 0; c < 8; c++) {
        const x = i + CORNER[c][0], y = j + CORNER[c][1], z = k + CORNER[c][2];
        grad[3 * c] = (gc(x + 1, y, z) - gc(x - 1, y, z)) * 0.5;
        grad[3 * c + 1] = (gc(x, y + 1, z) - gc(x, y - 1, z)) * 0.5;
        grad[3 * c + 2] = (gc(x, y, z + 1) - gc(x, y, z - 1)) * 0.5;
      }
      let gx = 0, gy = 0, gz = 0;
      for (let c = 0; c < 8; c++) {
        const wt = (CORNER[c][0] ? u : 1 - u) * (CORNER[c][1] ? v : 1 - v) * (CORNER[c][2] ? w : 1 - w);
        gx += wt * grad[3 * c]; gy += wt * grad[3 * c + 1]; gz += wt * grad[3 * c + 2];
      }
      // world gradient = M^-T * index gradient; outward normal = -grad
      let wx = G[0] * gx + G[1] * gy + G[2] * gz, wy = G[3] * gx + G[4] * gy + G[5] * gz, wz = G[6] * gx + G[7] * gy + G[8] * gz;
      const l = Math.hypot(wx, wy, wz) || 1;
      wx = -wx / l; wy = -wy / l; wz = -wz / l;
      const fi = i + u, fj = j + v, fk = k + w;
      pos.push(o[0] + fi * a1[0] + fj * a2[0] + fk * a3[0], o[1] + fi * a1[1] + fj * a2[1] + fk * a3[1], o[2] + fi * a1[2] + fj * a2[2] + fk * a3[2]);
      nrm.push(wx, wy, wz);
      cellIdx[(i * c2 + j) * c3 + k] = pos.length / 3 - 1;
    }
    const cell = (i, j, k) => (i < 0 || j < 0 || k < 0 || i >= c1 || j >= c2 || k >= c3) ? -1 : cellIdx[(i * c2 + j) * c3 + k];
    const idx = [];
    const quad = (a, b, c, d) => { if (a < 0 || b < 0 || c < 0 || d < 0) return; idx.push(a, b, c, a, c, d); };
    for (let i = 0; i < n1; i++) for (let j = 0; j < n2; j++) for (let k = 0; k < n3; k++) {
      const p = g(i, j, k) > 0;
      if (i + 1 < n1 && (g(i + 1, j, k) > 0) !== p) quad(cell(i, j - 1, k - 1), cell(i, j, k - 1), cell(i, j, k), cell(i, j - 1, k));
      if (j + 1 < n2 && (g(i, j + 1, k) > 0) !== p) quad(cell(i - 1, j, k - 1), cell(i, j, k - 1), cell(i, j, k), cell(i - 1, j, k));
      if (k + 1 < n3 && (g(i, j, k + 1) > 0) !== p) quad(cell(i - 1, j - 1, k), cell(i, j - 1, k), cell(i, j, k), cell(i - 1, j, k));
    }
    return { positions: new Float32Array(pos), normals: new Float32Array(nrm), indices: Uint32Array.from(idx) };
  }

  function stats(vol) {
    let mn = Infinity, mx = -Infinity, absMax = 0;
    const d = vol.data;
    for (let i = 0; i < d.length; i++) { const v = d[i]; if (v < mn) mn = v; if (v > mx) mx = v; }
    absMax = Math.max(Math.abs(mn), Math.abs(mx));
    return { min: mn, max: mx, absMax };
  }

  BM.iso = { extract, stats };
})(typeof window !== 'undefined' ? (window.BM = window.BM || {}) : (globalThis.BM = globalThis.BM || {}));
