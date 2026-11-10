/* BlueMolecule — structure, crystal and periodic-code exporters for the "Save as" registry.
 *
 * Registers with BM.registerSaveFormat (ctx: see writers.js):
 *   sdf, sdf-frames, mol2, cml, cjson, cif, extxyz, extxyz-frames, poscar, qe-in, cp2k-in, xsf, gro,
 *   lammps-data, smiles
 * Bonds come from ctx.bonds; without them, connectivity is guessed from covalent radii and bond orders
 * are perceived with BM.chem when it is loaded. Periodic formats need a cell: without ctx.cell the
 * structure is centred in a cubic box of 20 Å (larger when the structure would not fit with 5 Å of
 * vacuum) and the file says so in a comment; LAMMPS keeps the coordinates and moves the box instead.
 * SMILES needs BM.chem.toSmiles (chem.js); it is looked up when the writer runs.
 */
(function (BM) {
  'use strict';
  const RAD = Math.PI / 180, HARTREE_EV = 27.211386245988;

  // ------------------------------------------------------------------ helpers
  const sym = z => (z > 0 ? BM.el.symbol(z) : 'X');
  /** Fixed-width field (SDF V2000, GRO columns). */
  function fx(x, w, d) {
    const s = Number.isFinite(x) ? (Math.abs(x) < 0.5 * Math.pow(10, -d) ? 0 : x).toFixed(d) : 'NaN';
    return s.length >= w ? s : ' '.repeat(w - s.length) + s;
  }
  /** Free-format field: like fx, but always separated from the previous field by a blank. */
  const fs = (x, w, d) => { const s = fx(x, w, d); return s[0] === ' ' ? s : ' ' + s; };
  const ix = (v, w) => { const s = String(v); return s.length >= w ? s : ' '.repeat(w - s.length) + s; };
  const oneLine = s => String(s === undefined || s === null ? '' : s).replace(/[\r\n\t]+/g, ' ').trim();
  const isNum = v => typeof v === 'number' && Number.isFinite(v);
  const flat = c => c.reduce((a, v) => a.concat(Array.from(v)), []);
  function check(c) {
    if (!c || !c.Z || !c.xyz || c.xyz.length < 3 * c.Z.length) throw new Error('Writer needs Z and xyz arrays of matching size.');
    if (!c.Z.length) throw new Error('There are no atoms to write.');
  }
  const det3 = m => m[0][0] * (m[1][1] * m[2][2] - m[1][2] * m[2][1]) - m[0][1] * (m[1][0] * m[2][2] - m[1][2] * m[2][0]) +
    m[0][2] * (m[1][0] * m[2][1] - m[1][1] * m[2][0]);
  function inv3(M) {
    const [[a, b, c], [d, e, f], [g, h, i]] = M;
    const A = e * i - f * h, B = -(d * i - f * g), C = d * h - e * g, det = a * A + b * B + c * C;
    return [[A / det, -(b * i - c * h) / det, (b * f - c * e) / det], [B / det, (a * i - c * g) / det, -(a * f - c * d) / det],
      [C / det, -(a * h - b * g) / det, (a * e - b * d) / det]];
  }
  const c2f = (I, x, y, z) => [x * I[0][0] + y * I[1][0] + z * I[2][0], x * I[0][1] + y * I[1][1] + z * I[2][1], x * I[0][2] + y * I[1][2] + z * I[2][2]];
  const f2c = (C, u, v, w) => [u * C[0][0] + v * C[1][0] + w * C[2][0], u * C[0][1] + v * C[1][1] + w * C[2][1], u * C[0][2] + v * C[1][2] + w * C[2][2]];
  const validCell = c => Array.isArray(c) && c.length === 3 && c.every(v => v && v.length === 3 && Array.from(v).every(Number.isFinite)) && Math.abs(det3(c)) > 1e-6;
  /** [a, b, c, alpha, beta, gamma] of a cell given as row vectors. */
  function cellParams(C) {
    const len = v => Math.hypot(v[0], v[1], v[2]);
    const ang = (u, v) => Math.acos(Math.max(-1, Math.min(1, (u[0] * v[0] + u[1] * v[1] + u[2] * v[2]) / (len(u) * len(v))))) / RAD;
    return [len(C[0]), len(C[1]), len(C[2]), ang(C[1], C[2]), ang(C[0], C[2]), ang(C[0], C[1])];
  }
  const fileStem = ctx => String(ctx.baseName || ctx.title || 'structure').replace(/[^\w.-]+/g, '_').replace(/^[._]+/, '') || 'structure';
  const titleOf = ctx => oneLine(ctx.title || ctx.baseName || '');
  /** The cell for periodic formats: ctx.cell, or a cubic box with the structure at its centre. */
  function boxed(ctx, center) {
    const n = ctx.Z.length, xyz = new Float64Array(3 * n);
    for (let k = 0; k < 3 * n; k++) xyz[k] = +ctx.xyz[k];
    if (validCell(ctx.cell)) return { cell: ctx.cell.map(v => Array.from(v, Number)), xyz, note: null, boxed: false, lo: [0, 0, 0] };
    const mn = [Infinity, Infinity, Infinity], mx = [-Infinity, -Infinity, -Infinity];
    for (let i = 0; i < n; i++) for (let d = 0; d < 3; d++) { mn[d] = Math.min(mn[d], xyz[3 * i + d]); mx[d] = Math.max(mx[d], xyz[3 * i + d]); }
    const L = Math.max(20, Math.ceil(Math.max(mx[0] - mn[0], mx[1] - mn[1], mx[2] - mn[2]) + 10));
    const mid = [0, 1, 2].map(d => (mn[d] + mx[d]) / 2);
    if (center) for (let i = 0; i < n; i++) for (let d = 0; d < 3; d++) xyz[3 * i + d] += L / 2 - mid[d];
    return { cell: [[L, 0, 0], [0, L, 0], [0, 0, L]], xyz, boxed: true, lo: mid.map(m => m - L / 2),
      note: 'No unit cell in the source: the structure was ' + (center ? 'centred in' : 'placed in') + ' a ' + L + ' Angstrom cubic box.' };
  }
  function speciesOrder(Z) { const s = []; for (const z of Z) if (s.indexOf(z) < 0) s.push(z); return s; }
  const curEnergy = ctx => { const f = ctx.frames && ctx.frames[ctx.frameIndex]; return f && isNum(f.energy) ? f.energy : null; };

  // ------------------------------------------------------------------ bonds
  function guessPairs(Z, X) {
    const n = Z.length, cov = BM.el.covalent, out = [], grid = new Map();
    let rmax = 0;
    for (const z of Z) if (z > 0) rmax = Math.max(rmax, cov[z] || 0);
    const cut = 2.3 * rmax + 0.1, key = (a, b, c) => a + ',' + b + ',' + c;
    for (let i = 0; i < n; i++) {
      if (!(Z[i] > 0)) continue;
      const g = [0, 1, 2].map(d => Math.floor(X[3 * i + d] / cut));
      for (let a = -1; a <= 1; a++) for (let b = -1; b <= 1; b++) for (let c = -1; c <= 1; c++) {
        const list = grid.get(key(g[0] + a, g[1] + b, g[2] + c));
        if (!list) continue;
        for (const j of list) {
          const d = Math.hypot(X[3 * i] - X[3 * j], X[3 * i + 1] - X[3 * j + 1], X[3 * i + 2] - X[3 * j + 2]);
          if (d > 0.1 && d < 1.15 * (cov[Z[i]] + cov[Z[j]])) out.push([j, i, 1]);
        }
      }
      const k = key(g[0], g[1], g[2]);
      if (!grid.has(k)) grid.set(k, []);
      grid.get(k).push(i);
    }
    return out.sort((p, q) => p[0] - q[0] || p[1] - q[1]);
  }
  /** BM.chem.perceiveBondOrders returns the bond list ([[i, j, order]]); older builds returned {bonds}. */
  const perceived = p => (Array.isArray(p) ? p : p && Array.isArray(p.bonds) ? p.bonds : null);
  /** {bonds: [[i, j, order]], fc: formal charges or null} for formats that store bonds. */
  function bondsOf(ctx) {
    const fc0 = ctx.formalCharges && ctx.formalCharges.length === ctx.Z.length ? Array.from(ctx.formalCharges, q => +q || 0) : null;
    if (ctx.bonds && ctx.bonds.length) return { bonds: ctx.bonds.map(b => [b[0], b[1], b[2] || 1]), fc: fc0 };
    const pairs = guessPairs(ctx.Z, ctx.xyz);
    if (pairs.length && BM.chem && BM.chem.perceiveBondOrders) {
      try {
        const p = BM.chem.perceiveBondOrders({ Z: Array.from(ctx.Z), xyz: ctx.xyz, bonds: pairs });
        const list = perceived(p);
        if (list && list.length === pairs.length) return { bonds: list.map(b => [b[0], b[1], b[2] || 1]), fc: fc0 && fc0.some(Boolean) ? fc0 : (p.formalCharges || fc0) };
      } catch (e) { /* keep single bonds */ }
    }
    return { bonds: pairs, fc: fc0 };
  }
  /** Integer bond orders (aromatic 1.5 bonds kekulized with BM.chem when available). */
  function integerOrders(Z, xyz, bonds) {
    if (bonds.every(b => b[2] === Math.round(b[2]))) return bonds.map(b => b[2]);
    if (BM.chem && BM.chem.perceiveBondOrders) {
      try {
        const p = BM.chem.perceiveBondOrders({ Z: Array.from(Z), xyz, bonds: bonds.map(b => [b[0], b[1], 1]) });
        const m = new Map((perceived(p) || []).map(b => [Math.min(b[0], b[1]) + ',' + Math.max(b[0], b[1]), b[2]]));
        return bonds.map(b => (b[2] === Math.round(b[2]) ? b[2] : m.get(Math.min(b[0], b[1]) + ',' + Math.max(b[0], b[1])) || 1));
      } catch (e) { /* fall through */ }
    }
    return bonds.map(b => Math.max(1, Math.round(b[2])));
  }

  // ------------------------------------------------------------------ SDF
  const sdfOrder = o => (o === 1.5 ? 4 : Math.max(1, Math.min(3, Math.round(o || 1))));
  function stamp() {
    const d = new Date(), p2 = v => String(v).padStart(2, '0');
    return p2(d.getMonth() + 1) + p2(d.getDate()) + p2(d.getFullYear() % 100) + p2(d.getHours()) + p2(d.getMinutes());
  }
  function molBlock(Z, X, bonds, fc, title, props) {
    const n = Z.length, out = [oneLine(title).slice(0, 80), '  BlueMol ' + stamp() + '3D', ''];
    const el = i => (Z[i] > 0 ? sym(Z[i]) : 'Q');
    if (n <= 999 && bonds.length <= 999) {
      out.push(ix(n, 3) + ix(bonds.length, 3) + '  0  0  0  0  0  0  0  0999 V2000');
      for (let i = 0; i < n; i++) out.push(fx(X[3 * i], 10, 4) + fx(X[3 * i + 1], 10, 4) + fx(X[3 * i + 2], 10, 4) + ' ' + el(i).padEnd(3) + ' 0  0  0  0  0  0  0  0  0  0  0  0');
      for (const b of bonds) out.push(ix(b[0] + 1, 3) + ix(b[1] + 1, 3) + ix(sdfOrder(b[2]), 3) + '  0');
      const ch = [];
      if (fc) fc.forEach((q, i) => { if (q) ch.push([i + 1, q]); });
      for (let k = 0; k < ch.length; k += 8) { const part = ch.slice(k, k + 8); out.push('M  CHG' + ix(part.length, 3) + part.map(p => ix(p[0], 4) + ix(p[1], 4)).join('')); }
    } else {
      out.push('  0  0  0     0  0            999 V3000', 'M  V30 BEGIN CTAB', 'M  V30 COUNTS ' + n + ' ' + bonds.length + ' 0 0 0', 'M  V30 BEGIN ATOM');
      for (let i = 0; i < n; i++) out.push('M  V30 ' + (i + 1) + ' ' + el(i) + ' ' + X[3 * i].toFixed(6) + ' ' + X[3 * i + 1].toFixed(6) + ' ' + X[3 * i + 2].toFixed(6) + ' 0' + (fc && fc[i] ? ' CHG=' + fc[i] : ''));
      out.push('M  V30 END ATOM');
      if (bonds.length) {
        out.push('M  V30 BEGIN BOND');
        bonds.forEach((b, k) => out.push('M  V30 ' + (k + 1) + ' ' + sdfOrder(b[2]) + ' ' + (b[0] + 1) + ' ' + (b[1] + 1)));
        out.push('M  V30 END BOND');
      }
      out.push('M  V30 END CTAB');
    }
    out.push('M  END');
    for (const [k, v] of props) out.push('>  <' + k + '>', String(v), '');
    out.push('$$$$');
    return out.join('\n') + '\n';
  }
  function writeSdf(ctx, all) {
    check(ctx);
    const frames = all && ctx.frames && ctx.frames.length ? ctx.frames : [{ Z: ctx.Z, xyz: ctx.xyz, energy: curEnergy(ctx) }];
    const main = bondsOf(ctx), title = titleOf(ctx);
    return frames.map((f, k) => {
      const same = f.Z.length === ctx.Z.length && Array.from(f.Z).every((z, i) => z === ctx.Z[i]);
      const b = same ? main : bondsOf({ Z: f.Z, xyz: f.xyz, bonds: null, formalCharges: null });
      const props = isNum(f.energy) ? [['ENERGY', f.energy.toFixed(10)]] : [];
      return molBlock(f.Z, f.xyz, b.bonds, b.fc, frames.length > 1 ? (title ? title + ' ' : '') + '#' + (k + 1) : title, props);
    }).join('');
  }

  // ------------------------------------------------------------------ Tripos Mol2
  function sybylTypes(Z, bonds) {
    const n = Z.length, adj = Array.from({ length: n }, () => []);
    for (const [i, j, o] of bonds) { adj[i].push([j, o]); adj[j].push([i, o]); }
    const maxO = i => adj[i].reduce((m, e) => (e[1] === 1.5 ? m : Math.max(m, e[1])), 0);
    const arom = i => adj[i].some(e => e[1] === 1.5);
    const nDouble = i => adj[i].filter(e => e[1] === 2).length;
    const carbonyl = j => Z[j] === 6 && adj[j].some(e => e[1] === 2 && (Z[e[0]] === 8 || Z[e[0]] === 16));
    return Z.map((z, i) => {
      switch (z) {
        case 0: return 'Du';
        case 6: return arom(i) ? 'C.ar' : maxO(i) >= 3 || nDouble(i) >= 2 ? 'C.1' : maxO(i) === 2 ? 'C.2' : 'C.3';
        case 7: return arom(i) ? 'N.ar' : maxO(i) >= 3 ? 'N.1' : maxO(i) === 2 ? 'N.2' : adj[i].length >= 4 ? 'N.4' : adj[i].some(e => e[1] === 1 && carbonyl(e[0])) ? 'N.am' : 'N.3';
        case 8: return maxO(i) >= 2 ? 'O.2' : 'O.3';
        case 16: return maxO(i) >= 2 ? 'S.2' : 'S.3';
        case 15: return 'P.3';
        default: return sym(z);
      }
    });
  }
  function atomNamesOf(ctx) {
    const cnt = {};
    return Array.from(ctx.Z, (z, i) => {
      const s = sym(z);
      cnt[s] = (cnt[s] || 0) + 1;
      const a = ctx.atomNames && ctx.atomNames[i] ? String(ctx.atomNames[i]).trim().replace(/\s+/g, '_') : '';
      return a || s + cnt[s];
    });
  }
  function chargeScheme(ctx, prefer) {
    const n = ctx.Z.length, C = ctx.charges;
    if (!C) return null;
    const okArr = k => Array.isArray(C[k]) && C[k].length === n && C[k].every(Number.isFinite);
    const k = okArr(prefer) ? prefer : Object.keys(C).find(okArr);
    return k ? C[k] : null;
  }
  function writeMol2(ctx) {
    check(ctx);
    const n = ctx.Z.length, { bonds } = bondsOf(ctx), types = sybylTypes(ctx.Z, bonds), names = atomNamesOf(ctx);
    const q = chargeScheme(ctx, 'Mol2');
    const subs = [], subOf = new Int32Array(n), map = new Map();
    for (let i = 0; i < n; i++) {
      const rn = ctx.resNames && ctx.resNames[i] ? String(ctx.resNames[i]).trim().replace(/\s+/g, '_') : 'UNL';
      const ri = ctx.resIds && Number.isFinite(+ctx.resIds[i]) ? +ctx.resIds[i] : 1;
      const ch = ctx.chains && ctx.chains[i] ? String(ctx.chains[i]).trim() : '';
      const key = ch + '|' + ri + '|' + rn;
      if (!map.has(key)) { map.set(key, subs.length); subs.push({ name: rn + ri, root: i + 1, chain: ch, res: rn }); }
      subOf[i] = map.get(key);
    }
    const L = ['@<TRIPOS>MOLECULE', titleOf(ctx) || '*****', ix(n, 5) + ix(bonds.length, 6) + ix(subs.length, 6) + '     0     0',
      ctx.resNames && subs.length > 1 ? 'BIOPOLYMER' : 'SMALL', q ? 'USER_CHARGES' : 'NO_CHARGES', '', '@<TRIPOS>ATOM'];
    for (let i = 0; i < n; i++) {
      L.push(ix(i + 1, 7) + ' ' + names[i].padEnd(8) + fs(ctx.xyz[3 * i], 12, 6) + fs(ctx.xyz[3 * i + 1], 12, 6) + fs(ctx.xyz[3 * i + 2], 12, 6) +
        ' ' + types[i].padEnd(6) + ix(subOf[i] + 1, 5) + '  ' + subs[subOf[i]].name.padEnd(8) + fs(q ? q[i] : 0, 11, 6));
    }
    L.push('@<TRIPOS>BOND');
    bonds.forEach((b, k) => {
      let t = b[2] === 1.5 ? 'ar' : String(Math.max(1, Math.round(b[2])));
      if (t === '1' && ((types[b[0]] === 'N.am' && ctx.Z[b[1]] === 6 && types[b[1]] === 'C.2') || (types[b[1]] === 'N.am' && ctx.Z[b[0]] === 6 && types[b[0]] === 'C.2'))) t = 'am';
      L.push(ix(k + 1, 6) + ix(b[0] + 1, 6) + ix(b[1] + 1, 6) + ' ' + t);
    });
    L.push('@<TRIPOS>SUBSTRUCTURE');
    subs.forEach((s, k) => L.push(ix(k + 1, 6) + ' ' + s.name.padEnd(8) + ix(s.root, 7) + ' ' + (ctx.resNames ? 'RESIDUE' : 'GROUP').padEnd(9) +
      '1 ' + (s.chain || '****').padEnd(5) + s.res.padEnd(5) + '0 ROOT'));
    if (validCell(ctx.cell)) L.push('@<TRIPOS>CRYSIN', cellParams(ctx.cell).map(v => fs(v, 12, 6)).join('') + '     1     1');
    return L.join('\n') + '\n';
  }

  // ------------------------------------------------------------------ CML
  const xmlEsc = s => String(s).replace(/&/g, '&amp;').replace(/</g, '&lt;').replace(/>/g, '&gt;').replace(/"/g, '&quot;');
  function writeCml(ctx) {
    check(ctx);
    const n = ctx.Z.length, { bonds, fc } = bondsOf(ctx), X = ctx.xyz;
    const L = ['<?xml version="1.0" encoding="UTF-8"?>', '<molecule xmlns="http://www.xml-cml.org/schema" id="m1" title="' + xmlEsc(titleOf(ctx)) + '">'];
    if (validCell(ctx.cell)) {
      const p = cellParams(ctx.cell);
      L.push('  <crystal>');
      ['a', 'b', 'c', 'alpha', 'beta', 'gamma'].forEach((k, i) => L.push('    <scalar title="' + k + '" dictRef="cml:' + k + '" units="units:' + (i < 3 ? 'angstrom' : 'degree') + '">' + p[i].toFixed(8) + '</scalar>'));
      L.push('  </crystal>');
    }
    L.push('  <atomArray>');
    for (let i = 0; i < n; i++) {
      L.push('    <atom id="a' + (i + 1) + '" elementType="' + (ctx.Z[i] > 0 ? sym(ctx.Z[i]) : 'Du') + '" x3="' + (+X[3 * i]).toFixed(8) + '" y3="' +
        (+X[3 * i + 1]).toFixed(8) + '" z3="' + (+X[3 * i + 2]).toFixed(8) + '"' + (fc && fc[i] ? ' formalCharge="' + fc[i] + '"' : '') + '/>');
    }
    L.push('  </atomArray>');
    if (bonds.length) {
      L.push('  <bondArray>');
      for (const b of bonds) L.push('    <bond atomRefs2="a' + (b[0] + 1) + ' a' + (b[1] + 1) + '" order="' + (b[2] === 1.5 ? 'A' : Math.max(1, Math.round(b[2]))) + '"/>');
      L.push('  </bondArray>');
    }
    L.push('</molecule>');
    return L.join('\n') + '\n';
  }

  // ------------------------------------------------------------------ Avogadro CJSON
  function pretty(v, ind) {
    ind = ind || '';
    if (Array.isArray(v)) {
      if (v.every(x => x === null || typeof x !== 'object')) return JSON.stringify(v);
      return '[\n' + v.map(x => ind + '  ' + pretty(x, ind + '  ')).join(',\n') + '\n' + ind + ']';
    }
    if (v && typeof v === 'object') return '{\n' + Object.entries(v).map(([k, x]) => ind + '  ' + JSON.stringify(k) + ': ' + pretty(x, ind + '  ')).join(',\n') + '\n' + ind + '}';
    return JSON.stringify(v);
  }
  function writeCjson(ctx) {
    check(ctx);
    const n = ctx.Z.length, r = v => +(+v).toFixed(10), X = Array.from({ length: 3 * n }, (_, k) => r(ctx.xyz[k]));
    const { bonds, fc } = bondsOf(ctx);
    const J = { chemicalJson: 1, name: titleOf(ctx), atoms: { elements: { number: Array.from(ctx.Z) }, coords: { '3d': X } } };
    if (validCell(ctx.cell)) {
      const p = cellParams(ctx.cell), I = inv3(ctx.cell), fr = [];
      for (let i = 0; i < n; i++) fr.push(...c2f(I, X[3 * i], X[3 * i + 1], X[3 * i + 2]).map(r));
      J.atoms.coords['3dFractional'] = fr;
      J.unitCell = { a: r(p[0]), b: r(p[1]), c: r(p[2]), alpha: r(p[3]), beta: r(p[4]), gamma: r(p[5]), cellVectors: flat(ctx.cell).map(r) };
    }
    if (fc && fc.some(Boolean)) J.atoms.formalCharges = fc;
    if (ctx.atomNames && ctx.atomNames.length === n) J.atoms.labels = Array.from(ctx.atomNames, s => String(s || ''));
    if (bonds.length) J.bonds = { connections: { index: flat(bonds.map(b => [b[0], b[1]])) }, order: integerOrders(ctx.Z, ctx.xyz, bonds) };
    if (ctx.charges) {
      const pc = {};
      for (const [k, v] of Object.entries(ctx.charges)) if (Array.isArray(v) && v.length === n) pc[k] = v.map(Number);
      if (Object.keys(pc).length) J.partialCharges = pc;
    }
    const V = ctx.vibrations;
    if (V && Array.isArray(V.freqs) && V.freqs.length && V.modes && V.modes.length === V.freqs.length) {
      J.vibrations = { frequencies: V.freqs.map(Number), modes: V.freqs.map((_, k) => k + 1), eigenVectors: V.modes.map(m => Array.from(m, r)) };
      if (V.ir && V.ir.length === V.freqs.length) J.vibrations.intensities = V.ir.map(Number);
      if (V.raman && V.raman.length === V.freqs.length) J.vibrations.ramanIntensities = V.raman.map(Number);
    }
    return pretty(J) + '\n';
  }

  // ------------------------------------------------------------------ CIF (P1)
  function writeCif(ctx) {
    check(ctx);
    const B = boxed(ctx, true), I = inv3(B.cell), p = cellParams(B.cell);
    const L = ['# CIF written by BlueMolecule (space group P1, all atoms listed)'];
    if (B.note) L.push('# ' + B.note);
    L.push('data_' + fileStem(ctx).replace(/[^\w.-]/g, '_'), "_symmetry_space_group_name_H-M   'P 1'", '_symmetry_Int_Tables_number      1');
    ['a', 'b', 'c'].forEach((k, i) => L.push(('_cell_length_' + k).padEnd(25) + p[i].toFixed(8)));
    ['alpha', 'beta', 'gamma'].forEach((k, i) => L.push(('_cell_angle_' + k).padEnd(25) + p[3 + i].toFixed(8)));
    L.push('_cell_volume'.padEnd(25) + Math.abs(det3(B.cell)).toFixed(6), 'loop_', '_symmetry_equiv_pos_as_xyz', "'x, y, z'", 'loop_',
      '_atom_site_label', '_atom_site_type_symbol', '_atom_site_fract_x', '_atom_site_fract_y', '_atom_site_fract_z', '_atom_site_occupancy');
    const cnt = {};
    ctx.Z.forEach((z, i) => {
      const s = sym(z);
      cnt[s] = (cnt[s] || 0) + 1;
      const f = c2f(I, B.xyz[3 * i], B.xyz[3 * i + 1], B.xyz[3 * i + 2]);
      L.push((s + cnt[s]).padEnd(8) + s.padEnd(3) + f.map(v => fs(v, 15, 10)).join('') + ' 1.0');
    });
    return L.join('\n') + '\n';
  }

  // ------------------------------------------------------------------ extended XYZ
  function extxyzBlock(Z, X, cell, energy) {
    const head = [];
    const periodic = validCell(cell);
    if (periodic) head.push('Lattice="' + flat(cell).map(v => (+v).toFixed(10)).join(' ') + '"');
    head.push('Properties=species:S:1:pos:R:3');
    if (isNum(energy)) head.push('energy=' + (energy * HARTREE_EV).toFixed(8));
    head.push('pbc="' + (periodic ? 'T T T' : 'F F F') + '"');
    const L = [String(Z.length), head.join(' ')];
    for (let i = 0; i < Z.length; i++) L.push(sym(Z[i]).padEnd(3) + fs(X[3 * i], 18, 10) + fs(X[3 * i + 1], 18, 10) + fs(X[3 * i + 2], 18, 10));
    return L.join('\n') + '\n';
  }
  function writeExtxyz(ctx, all) {
    check(ctx);
    if (!all || !ctx.frames || !ctx.frames.length) return extxyzBlock(ctx.Z, ctx.xyz, ctx.cell, curEnergy(ctx));
    return ctx.frames.map(f => extxyzBlock(f.Z, f.xyz, ctx.cell, f.energy)).join('');
  }

  // ------------------------------------------------------------------ VASP POSCAR (VASP 5, species grouped)
  function writePoscar(ctx) {
    check(ctx);
    const B = boxed(ctx, true), I = inv3(B.cell), sp = speciesOrder(ctx.Z);
    const L = [(titleOf(ctx) || 'structure') + (B.note ? '   # ' + B.note : ''), '   1.0'];
    for (const v of B.cell) L.push(v.map(x => fs(x, 22, 14)).join(''));
    L.push(sp.map(z => sym(z).padStart(5)).join(''), sp.map(z => ix(ctx.Z.filter(q => q === z).length, 5)).join(''), 'Direct');
    for (const z of sp) ctx.Z.forEach((q, i) => {
      if (q !== z) return;
      L.push(c2f(I, B.xyz[3 * i], B.xyz[3 * i + 1], B.xyz[3 * i + 2]).map(v => fs(v, 20, 14)).join('') + '   ' + sym(z));
    });
    return L.join('\n') + '\n';
  }

  // ------------------------------------------------------------------ Quantum ESPRESSO pw.x input
  function chargeMult(ctx) {
    const m = ctx.meta || {};
    return { charge: isNum(m.charge) ? m.charge : 0, mult: isNum(m.multiplicity) && m.multiplicity >= 1 ? Math.round(m.multiplicity) : 1 };
  }
  function writeQe(ctx) {
    check(ctx);
    const B = boxed(ctx, true), sp = speciesOrder(ctx.Z), { charge, mult } = chargeMult(ctx);
    const kgrid = B.boxed ? null : B.cell.map(v => Math.max(1, Math.ceil(30 / Math.hypot(v[0], v[1], v[2]))));
    const L = ['! pw.x input written by BlueMolecule' + (titleOf(ctx) ? ': ' + titleOf(ctx) : '')];
    if (B.note) L.push('! ' + B.note);
    L.push('&CONTROL', "  calculation = 'scf'", "  prefix = '" + fileStem(ctx).replace(/[^\w-]/g, '_') + "'", "  outdir = './tmp'", "  pseudo_dir = './pseudo'",
      '  tprnfor = .true.', '  tstress = .true.', '/', '&SYSTEM', '  ibrav = 0', '  nat = ' + ctx.Z.length, '  ntyp = ' + sp.length,
      '  ecutwfc = 40.0', '  ecutrho = 320.0', "  occupations = 'smearing'", "  smearing = 'mv'", '  degauss = 0.01');
    if (charge) L.push('  tot_charge = ' + charge);
    if (mult > 1) L.push('  nspin = 2', '  tot_magnetization = ' + (mult - 1));
    L.push('/', '&ELECTRONS', '  conv_thr = 1.0d-8', '  mixing_beta = 0.4', '/', '&IONS', '/', '&CELL', '/',
      '! Placeholder pseudopotentials: replace <El>.UPF with the files in pseudo_dir and adjust ecutwfc/ecutrho.', 'ATOMIC_SPECIES');
    for (const z of sp) L.push('  ' + sym(z).padEnd(3) + fs(BM.el.mass[z] || 1, 12, 6) + '  ' + sym(z) + '.UPF');
    L.push('CELL_PARAMETERS angstrom');
    for (const v of B.cell) L.push(v.map(x => fs(x, 20, 12)).join(''));
    L.push('ATOMIC_POSITIONS angstrom');
    ctx.Z.forEach((z, i) => L.push('  ' + sym(z).padEnd(3) + fs(B.xyz[3 * i], 20, 12) + fs(B.xyz[3 * i + 1], 20, 12) + fs(B.xyz[3 * i + 2], 20, 12)));
    if (kgrid) L.push('K_POINTS automatic', '  ' + kgrid.join(' ') + ' 0 0 0'); else L.push('K_POINTS gamma');
    return L.join('\n') + '\n';
  }

  // ------------------------------------------------------------------ CP2K input
  function writeCp2k(ctx) {
    check(ctx);
    const B = boxed(ctx, true), sp = speciesOrder(ctx.Z), { charge, mult } = chargeMult(ctx);
    const L = ['# CP2K input written by BlueMolecule' + (titleOf(ctx) ? ': ' + titleOf(ctx) : '')];
    if (B.note) L.push('# ' + B.note);
    L.push('# Placeholder basis set and pseudopotential: check BASIS_MOLOPT and GTH_POTENTIALS for your elements.',
      '&GLOBAL', '  PROJECT ' + fileStem(ctx).replace(/[^\w-]/g, '_'), '  RUN_TYPE ENERGY', '  PRINT_LEVEL LOW', '&END GLOBAL',
      '&FORCE_EVAL', '  METHOD QUICKSTEP', '  &DFT', '    BASIS_SET_FILE_NAME BASIS_MOLOPT', '    POTENTIAL_FILE_NAME GTH_POTENTIALS');
    if (charge) L.push('    CHARGE ' + charge);
    if (mult > 1) L.push('    MULTIPLICITY ' + mult, '    UKS .TRUE.');
    L.push('    &MGRID', '      CUTOFF 400', '      REL_CUTOFF 50', '    &END MGRID', '    &QS', '      EPS_DEFAULT 1.0E-12', '    &END QS');
    if (B.boxed) L.push('    &POISSON', '      PERIODIC NONE', '      POISSON_SOLVER WAVELET', '    &END POISSON');
    L.push('    &SCF', '      SCF_GUESS ATOMIC', '      EPS_SCF 1.0E-6', '      MAX_SCF 50', '      &OT ON', '        MINIMIZER DIIS',
      '        PRECONDITIONER FULL_SINGLE_INVERSE', '      &END OT', '    &END SCF', '    &XC', '      &XC_FUNCTIONAL PBE', '      &END XC_FUNCTIONAL',
      '    &END XC', '  &END DFT', '  &SUBSYS', '    &CELL');
    ['A', 'B', 'C'].forEach((k, i) => L.push('      ' + k + B.cell[i].map(x => fs(x, 20, 12)).join('')));
    L.push('      PERIODIC ' + (B.boxed ? 'NONE' : 'XYZ'), '    &END CELL', '    &COORD');
    ctx.Z.forEach((z, i) => L.push('      ' + sym(z).padEnd(3) + fs(B.xyz[3 * i], 20, 12) + fs(B.xyz[3 * i + 1], 20, 12) + fs(B.xyz[3 * i + 2], 20, 12)));
    L.push('    &END COORD');
    for (const z of sp) L.push('    &KIND ' + sym(z), '      BASIS_SET DZVP-MOLOPT-SR-GTH', '      POTENTIAL GTH-PBE', '    &END KIND');
    L.push('  &END SUBSYS', '&END FORCE_EVAL');
    return L.join('\n') + '\n';
  }

  // ------------------------------------------------------------------ XSF (with DATAGRID_3D)
  function writeXsf(ctx) {
    check(ctx);
    const n = ctx.Z.length, X = ctx.xyz;
    const L = ['# XSF written by BlueMolecule' + (titleOf(ctx) ? ': ' + titleOf(ctx) : '')];
    const atomLine = i => ix(ctx.Z[i], 4) + fs(X[3 * i], 18, 10) + fs(X[3 * i + 1], 18, 10) + fs(X[3 * i + 2], 18, 10);
    if (validCell(ctx.cell)) {
      L.push('CRYSTAL', 'PRIMVEC');
      for (const v of ctx.cell) L.push(Array.from(v).map(x => fs(+x, 18, 10)).join(''));
      L.push('PRIMCOORD', ix(n, 6) + ' 1');
    } else L.push('ATOMS');
    for (let i = 0; i < n; i++) L.push(atomLine(i));
    const V = ctx.volume;
    if (V && V.data && V.dims && V.axes && V.origin) {
      const [n1, n2, n3] = V.dims;
      if (V.data.length < n1 * n2 * n3) throw new Error('The volume data do not match its grid dimensions.');
      const nm = String(V.name || 'grid').replace(/\W+/g, '_').replace(/^_+|_+$/g, '') || 'grid';
      L.push('BEGIN_BLOCK_DATAGRID_3D', '  ' + nm, '  BEGIN_DATAGRID_3D_' + nm, '    ' + n1 + ' ' + n2 + ' ' + n3,
        '  ' + V.origin.map(x => fs(+x, 18, 10)).join(''));
      V.axes.forEach((a, k) => L.push('  ' + a.map(x => fs(x * Math.max(1, V.dims[k] - 1), 18, 10)).join('')));
      const row = [];
      for (let c = 0; c < n3; c++) for (let b = 0; b < n2; b++) for (let a = 0; a < n1; a++) {   // x fastest
        row.push((+V.data[(a * n2 + b) * n3 + c]).toExponential(7));
        if (row.length === 6) { L.push('    ' + row.join(' ')); row.length = 0; }
      }
      if (row.length) L.push('    ' + row.join(' '));
      L.push('  END_DATAGRID_3D', 'END_BLOCK_DATAGRID_3D');
    }
    return L.join('\n') + '\n';
  }

  // ------------------------------------------------------------------ GROMACS .gro (nm)
  function writeGro(ctx) {
    check(ctx);
    const B = boxed(ctx, true), n = ctx.Z.length, names = atomNamesOf(ctx), nm = v => fx(v / 10, 12, 7);
    const L = [(titleOf(ctx) || 'structure') + (B.note ? ' (' + B.note + ')' : ''), ix(n, 5)];
    for (let i = 0; i < n; i++) {
      const rid = ctx.resIds && Number.isFinite(+ctx.resIds[i]) ? +ctx.resIds[i] : 1;
      const rn = ctx.resNames && ctx.resNames[i] ? String(ctx.resNames[i]).trim() : 'MOL';
      L.push(ix(((rid % 100000) + 100000) % 100000, 5) + rn.slice(0, 5).padEnd(5) + ix(names[i].slice(0, 5), 5) + ix((i + 1) % 100000, 5) +
        nm(B.xyz[3 * i]) + nm(B.xyz[3 * i + 1]) + nm(B.xyz[3 * i + 2]));
    }
    const c = B.cell;
    const off = [c[0][1], c[0][2], c[1][0], c[1][2], c[2][0], c[2][1]];
    L.push((off.every(v => Math.abs(v) < 1e-10) ? [c[0][0], c[1][1], c[2][2]] : [c[0][0], c[1][1], c[2][2], c[0][1], c[0][2], c[1][0], c[1][2], c[2][0], c[2][1]]).map(nm).join(''));
    return L.join('\n') + '\n';
  }

  // ------------------------------------------------------------------ LAMMPS data (atomic style)
  function writeLammps(ctx) {
    check(ctx);
    const n = ctx.Z.length, sp = speciesOrder(ctx.Z), type = new Map(sp.map((z, k) => [z, k + 1]));
    let xyz = new Float64Array(3 * n), lo = [0, 0, 0], lx, ly, lz, xy = 0, xz = 0, yz = 0, note = null;
    for (let k = 0; k < 3 * n; k++) xyz[k] = +ctx.xyz[k];
    if (validCell(ctx.cell)) {
      const C = ctx.cell.map(v => Array.from(v, Number));
      const restricted = Math.abs(C[0][1]) < 1e-10 && Math.abs(C[0][2]) < 1e-10 && Math.abs(C[1][2]) < 1e-10 && C[0][0] > 0 && C[1][1] > 0 && C[2][2] > 0;
      if (restricted) { lx = C[0][0]; xy = C[1][0]; ly = C[1][1]; xz = C[2][0]; yz = C[2][1]; lz = C[2][2]; }
      else {
        const [a, b, c, al, be, ga] = cellParams(C);
        lx = a; xy = b * Math.cos(ga * RAD); ly = Math.sqrt(b * b - xy * xy); xz = c * Math.cos(be * RAD);
        yz = (b * c * Math.cos(al * RAD) - xy * xz) / ly; lz = Math.sqrt(Math.max(0, c * c - xz * xz - yz * yz));
        const I = inv3(C), N = [[lx, 0, 0], [xy, ly, 0], [xz, yz, lz]];
        for (let i = 0; i < n; i++) xyz.set(f2c(N, ...c2f(I, xyz[3 * i], xyz[3 * i + 1], xyz[3 * i + 2])), 3 * i);
        note = 'cell rotated into the LAMMPS orientation (a along x, b in the xy plane)' + (det3(C) < 0 ? ' and made right-handed' : '');
      }
    } else {
      const B = boxed(ctx, false);
      lo = B.lo; lx = ly = lz = B.cell[0][0];
      note = B.note.replace(/^No unit cell in the source: the structure was placed in/, 'no unit cell: box of').replace(/\.$/, '') + ' around the structure';
    }
    const tilt = Math.abs(xy) > 1e-12 || Math.abs(xz) > 1e-12 || Math.abs(yz) > 1e-12;
    const g = v => (+v).toFixed(10);
    const L = ['LAMMPS data file written by BlueMolecule' + (titleOf(ctx) ? ': ' + titleOf(ctx) : '') + (note ? ' (' + note + ')' : ''), '',
      n + ' atoms', sp.length + ' atom types', '',
      g(lo[0]) + ' ' + g(lo[0] + lx) + ' xlo xhi', g(lo[1]) + ' ' + g(lo[1] + ly) + ' ylo yhi', g(lo[2]) + ' ' + g(lo[2] + lz) + ' zlo zhi'];
    if (tilt) L.push(g(xy) + ' ' + g(xz) + ' ' + g(yz) + ' xy xz yz');
    L.push('', 'Masses', '');
    sp.forEach((z, k) => L.push((k + 1) + ' ' + (BM.el.mass[z] || 1).toFixed(6) + '  # ' + sym(z)));
    L.push('', 'Atoms  # atomic', '');
    for (let i = 0; i < n; i++) L.push((i + 1) + ' ' + type.get(ctx.Z[i]) + ' ' + g(xyz[3 * i]) + ' ' + g(xyz[3 * i + 1]) + ' ' + g(xyz[3 * i + 2]));
    return L.join('\n') + '\n';
  }

  // ------------------------------------------------------------------ SMILES
  function writeSmiles(ctx) {
    check(ctx);
    const chem = BM.chem;
    if (!chem || !chem.toSmiles) throw new Error('SMILES export needs the chemistry module (chem.js).');
    const n = ctx.Z.length, keep = [], map = new Int32Array(n).fill(-1);
    for (let i = 0; i < n; i++) if (ctx.Z[i] > 0) { map[i] = keep.length; keep.push(i); }
    if (!keep.length) throw new Error('There are no real atoms to write.');
    const { bonds, fc } = bondsOf(ctx);
    const b2 = bonds.filter(b => map[b[0]] >= 0 && map[b[1]] >= 0).map(b => [map[b[0]], map[b[1]], b[2]]);
    const s = chem.toSmiles({ Z: keep.map(i => ctx.Z[i]), bonds: b2, formalCharges: fc ? keep.map(i => fc[i] || 0) : undefined });
    const t = titleOf(ctx);
    return s + (t ? '\t' + t : '') + '\n';
  }

  // ------------------------------------------------------------------ registration
  const R = f => BM.registerSaveFormat(f);
  R({ id: 'sdf', name: 'SD file', ext: 'sdf', kind: 'structure', write: c => writeSdf(c, false) });
  R({ id: 'sdf-frames', name: 'SD file (all geometries)', ext: 'sdf', kind: 'trajectory', needs: ['frames'], write: c => writeSdf(c, true) });
  R({ id: 'mol2', name: 'Tripos Mol2', ext: 'mol2', kind: 'structure', write: writeMol2 });
  R({ id: 'cml', name: 'Chemical Markup Language', ext: 'cml', kind: 'structure', write: writeCml });
  R({ id: 'cjson', name: 'Avogadro CJSON', ext: 'cjson', kind: 'structure', write: writeCjson });
  R({ id: 'cif', name: 'CIF (P1)', ext: 'cif', kind: 'structure', write: writeCif });
  R({ id: 'extxyz', name: 'Extended XYZ', ext: 'xyz', kind: 'structure', write: c => writeExtxyz(c, false) });
  R({ id: 'extxyz-frames', name: 'Extended XYZ (all geometries)', ext: 'xyz', kind: 'trajectory', needs: ['frames'], write: c => writeExtxyz(c, true) });
  R({ id: 'poscar', name: 'VASP POSCAR', ext: 'vasp', kind: 'structure', write: writePoscar });
  R({ id: 'qe-in', name: 'Quantum ESPRESSO pw.x input', ext: 'pwi', kind: 'structure', write: writeQe });
  R({ id: 'cp2k-in', name: 'CP2K input', ext: 'inp', kind: 'structure', write: writeCp2k });
  R({ id: 'xsf', name: 'XCrySDen XSF (+ volume)', ext: 'xsf', kind: 'structure', write: writeXsf });
  R({ id: 'gro', name: 'GROMACS .gro', ext: 'gro', kind: 'structure', write: writeGro });
  R({ id: 'lammps-data', name: 'LAMMPS data (atomic)', ext: 'data', kind: 'structure', write: writeLammps });
  R({ id: 'smiles', name: 'SMILES', ext: 'smi', kind: 'structure', write: writeSmiles });

  BM.structWriters = { sdf: writeSdf, mol2: writeMol2, cml: writeCml, cjson: writeCjson, cif: writeCif, extxyz: writeExtxyz, poscar: writePoscar,
    qe: writeQe, cp2k: writeCp2k, xsf: writeXsf, gro: writeGro, lammps: writeLammps, smiles: writeSmiles, cellParams, boxed };
})(typeof window !== 'undefined' ? (window.BM = window.BM || {}) : (globalThis.BM = globalThis.BM || {}));
