/* BlueMolecule — chemistry core.
 *   SMILES reading (with kekulization) and writing, bond-order perception, valence and hydrogen
 *   handling, distance-geometry style 3D structure generation and small editing helpers.
 * Original implementation, no dependencies.  Needs BM.el (elements.js); uses BM.ff when present.
 *
 * Molecule objects: { Z: number[], xyz: Float64Array(3N), bonds: [[i, j, order]], formalCharges?: number[] }
 * Nothing here mutates its inputs; functions that remove or reorder atoms also return `map` (old -> new, -1 deleted).
 */
(function (BM) {
  'use strict';
  const EL = () => BM.el;
  const DEG = 180 / Math.PI, RAD = Math.PI / 180;
  const TET = Math.acos(-1 / 3) * DEG; // 109.4712...

  // =================================================================== element classes and valences
  const SQUARE_PLANAR = new Set([78, 46, 79, 45, 77]); // Pt Pd Au Rh Ir
  const NOBLE = new Set([2, 10, 18, 36, 54, 86, 118]);

  /** Periodic-table group (1-18); lanthanides/actinides report 3. */
  function groupOf(z) {
    if (z === 1) return 1;
    if (z === 2) return 18;
    let g;
    if (z <= 10) g = z - 2;
    else if (z <= 18) g = z - 10;
    else if (z <= 36) return z - 18;
    else if (z <= 54) return z - 36;
    else if (z <= 56) return z - 54;
    else if (z <= 71) return 3;
    else if (z <= 86) return z - 68;
    else if (z <= 88) return z - 86;
    else if (z <= 103) return 3;
    else return Math.min(18, z - 100);
    return g >= 3 ? g + 10 : g;
  }
  function periodOf(z) { return z <= 2 ? 1 : z <= 10 ? 2 : z <= 18 ? 3 : z <= 36 ? 4 : z <= 54 ? 5 : z <= 86 ? 6 : 7; }
  function isMetal(z) {
    if (!(z > 1)) return false;
    const g = groupOf(z);
    if (g <= 12) return true;
    if (g === 13) return z !== 5;
    if (g === 14) return z === 50 || z === 82 || z === 114;
    if (g === 15) return z === 83 || z === 115;
    if (g === 16) return z === 84 || z === 116;
    return false;
  }
  /** Allowed valences of a non-metal for a given formal charge (null for metals / dummies). */
  function valenceList(z, charge) {
    charge = charge || 0;
    if (z === 1) return charge === 0 ? [1] : [0];
    if (!(z > 1) || isMetal(z)) return null;
    const ve = groupOf(z) - 10 - charge;
    if (ve < 0 || ve > 8) return [0];
    const base = ve <= 4 ? ve : 8 - ve;
    const out = [base];
    if (periodOf(z) >= 3 && ve > 4) for (let v = base + 2; v <= ve; v += 2) out.push(v);
    return out;
  }
  function valenceElectrons(z, charge) { return z === 1 ? 1 - (charge || 0) : groupOf(z) - 10 - (charge || 0); }
  // Default valences for implicit hydrogens on organic-subset SMILES atoms.
  const SMILES_VAL = { 5: [3], 6: [4], 7: [3, 5], 8: [2], 15: [3, 5], 16: [2, 4, 6], 9: [1], 17: [1], 35: [1], 53: [1] };
  function smallestAtLeast(list, x) { if (list) for (const v of list) if (v >= x - 1e-9) return v; return undefined; }

  // =================================================================== ideal bond lengths
  const XH_LEN = { 1: 0.74, 5: 1.19, 6: 1.09, 7: 1.01, 8: 0.96, 14: 1.48, 15: 1.42, 16: 1.34 };
  const ORDER_F = [[1, 1], [1.5, 0.915], [2, 0.875], [3, 0.79]];
  function cov(z) { const r = EL().covalent[z]; return r > 0 ? r : 1.5; }
  function orderFactor(o) {
    if (!(o > 1)) return 1;
    if (o >= 3) return 0.79;
    for (let k = 1; k < ORDER_F.length; k++) {
      const [o1, f1] = ORDER_F[k];
      if (o <= o1) { const [o0, f0] = ORDER_F[k - 1]; return f0 + (f1 - f0) * (o - o0) / (o1 - o0); }
    }
    return 0.79;
  }
  /** r0(i, j, order) = (rc_i + rc_j) * f(order); X-H bonds use fixed lengths. */
  function idealBondLength(z1, z2, order) {
    if (z1 === 1 || z2 === 1) {
      const x = z1 === 1 ? z2 : z1;
      return XH_LEN[x] !== undefined ? XH_LEN[x] : cov(x) + 0.31;
    }
    return (cov(z1) + cov(z2)) * orderFactor(order === undefined ? 1 : order);
  }

  // =================================================================== small vector helpers
  const vsub = (a, b) => [a[0] - b[0], a[1] - b[1], a[2] - b[2]];
  const vadd = (a, b) => [a[0] + b[0], a[1] + b[1], a[2] + b[2]];
  const vscale = (a, s) => [a[0] * s, a[1] * s, a[2] * s];
  const vdot = (a, b) => a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
  const vcross = (a, b) => [a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]];
  const vnorm = a => Math.sqrt(vdot(a, a));
  function vunit(a) { const l = vnorm(a); return l > 1e-12 ? [a[0] / l, a[1] / l, a[2] / l] : [0, 0, 0]; }
  const P = (xyz, i) => [xyz[3 * i], xyz[3 * i + 1], xyz[3 * i + 2]];
  /** A deterministic unit vector perpendicular to u. */
  function perpendicular(u) {
    const ax = Math.abs(u[0]), ay = Math.abs(u[1]), az = Math.abs(u[2]);
    const t = ax <= ay && ax <= az ? [1, 0, 0] : ay <= az ? [0, 1, 0] : [0, 0, 1];
    return vunit(vcross(u, t));
  }
  /** Seeded PRNG (xorshift32 with a splitmix-style seed scramble). */
  function makeRng(seed) {
    let x = (seed === undefined ? 12345 : seed) >>> 0;
    x = Math.imul(x ^ (x >>> 16), 0x45d9f3b) >>> 0; x = Math.imul(x ^ (x >>> 16), 0x45d9f3b) >>> 0; x = (x ^ (x >>> 16)) >>> 0;
    if (x === 0) x = 0x9e3779b9;
    return function () { x ^= x << 13; x >>>= 0; x ^= x >>> 17; x ^= x << 5; x >>>= 0; return x / 4294967296; };
  }

  // =================================================================== graph helpers
  /** adjacency with bond info: adj[i] = [[j, order, bondIndex], ...] */
  function adjOf(n, bonds) {
    const adj = Array.from({ length: n }, () => []);
    for (let k = 0; k < bonds.length; k++) {
      const b = bonds[k];
      adj[b[0]].push([b[1], b[2] === undefined ? 1 : b[2], k]);
      adj[b[1]].push([b[0], b[2] === undefined ? 1 : b[2], k]);
    }
    return adj;
  }
  /** Flags bonds that lie on a cycle (i.e. are not bridges). Iterative Tarjan low-link. */
  function ringBondFlags(n, bonds, adj) {
    adj = adj || adjOf(n, bonds);
    const disc = new Int32Array(n).fill(-1), low = new Int32Array(n);
    const ring = new Uint8Array(bonds.length).fill(1);
    let t = 0;
    for (let r = 0; r < n; r++) {
      if (disc[r] >= 0) continue;
      disc[r] = low[r] = t++;
      const stack = [[r, -1, 0]];
      while (stack.length) {
        const top = stack[stack.length - 1], u = top[0];
        if (top[2] < adj[u].length) {
          const [v, , e] = adj[u][top[2]++];
          if (e === top[1]) continue;
          if (disc[v] < 0) { disc[v] = low[v] = t++; stack.push([v, e, 0]); }
          else if (disc[v] < low[u]) low[u] = disc[v];
        } else {
          stack.pop();
          if (stack.length) {
            const p = stack[stack.length - 1][0];
            if (low[u] < low[p]) low[p] = low[u];
            if (low[u] > disc[p]) ring[top[1]] = 0;
          }
        }
      }
    }
    return ring;
  }
  /** Breadth-first shortest paths with banned atoms, reusable for one graph (nbrs[i] = neighbour index list). */
  function pathTool(nbrs) {
    const n = nbrs.length, seen = new Int32Array(n), par = new Int32Array(n), depth = new Int32Array(n), queue = new Int32Array(n);
    let stamp = 0;
    /** Shortest path src..dst (inclusive) avoiding ban1/ban2, at most maxEdges edges; skipDirect forbids the edge src-dst. */
    function path(src, dst, ban1, ban2, maxEdges, skipDirect) {
      stamp++;
      if (src === dst) return [src];
      let h = 0, t = 0;
      queue[t++] = src; seen[src] = stamp; depth[src] = 0; par[src] = -1;
      while (h < t) {
        const u = queue[h++];
        if (depth[u] >= maxEdges) continue;
        for (const v of nbrs[u]) {
          if (v === ban1 || v === ban2 || seen[v] === stamp) continue;
          if (skipDirect && u === src && v === dst) continue;
          seen[v] = stamp; par[v] = u; depth[v] = depth[u] + 1;
          if (v === dst) { const out = [v]; let x = u; while (x >= 0) { out.push(x); x = par[x]; } return out.reverse(); }
          queue[t++] = v;
        }
      }
      return null;
    }
    return { path };
  }
  function fragmentsOf(n, bonds) {
    const parent = new Int32Array(n);
    for (let i = 0; i < n; i++) parent[i] = i;
    const find = x => { while (parent[x] !== x) { parent[x] = parent[parent[x]]; x = parent[x]; } return x; };
    for (const b of bonds) { const a = find(b[0]), c = find(b[1]); if (a !== c) parent[Math.max(a, c)] = Math.min(a, c); }
    const groups = new Map();
    for (let i = 0; i < n; i++) { const r = find(i); if (!groups.has(r)) groups.set(r, []); groups.get(r).push(i); }
    return [...groups.values()];
  }

  // =================================================================== maximum matching (Edmonds' blossom algorithm)
  /** Maximum-cardinality matching of a general graph. adj: neighbour lists; init: starting matching;
   *  order: vertices in the order augmenting paths are searched from. Returns Int32Array mate (-1 = free). */
  function maxMatching(adj, init, order) {
    const n = adj.length;
    const mate = init ? Int32Array.from(init) : new Int32Array(n).fill(-1);
    if (!init) for (let v = 0; v < n; v++) if (mate[v] < 0) for (const u of adj[v]) if (mate[u] < 0) { mate[v] = u; mate[u] = v; break; }
    const parent = new Int32Array(n), base = new Int32Array(n), queue = new Int32Array(n);
    const queued = new Uint8Array(n), inBlossom = new Uint8Array(n), onPath = new Uint8Array(n);
    function lca(a, b) {
      onPath.fill(0);
      for (;;) { a = base[a]; onPath[a] = 1; if (mate[a] < 0) break; a = parent[mate[a]]; }
      for (;;) { b = base[b]; if (onPath[b]) return b; b = parent[mate[b]]; }
    }
    function markPath(v, b, child) {
      while (base[v] !== b) {
        inBlossom[base[v]] = 1; inBlossom[base[mate[v]]] = 1;
        parent[v] = child; child = mate[v]; v = parent[mate[v]];
      }
    }
    function augmentFrom(root) {
      parent.fill(-1); queued.fill(0);
      for (let i = 0; i < n; i++) base[i] = i;
      let h = 0, t = 0;
      queue[t++] = root; queued[root] = 1;
      while (h < t) {
        const v = queue[h++];
        for (const u of adj[v]) {
          if (base[v] === base[u] || mate[v] === u) continue;
          if (u === root || (mate[u] >= 0 && parent[mate[u]] >= 0)) {
            const b = lca(v, u);
            inBlossom.fill(0);
            markPath(v, b, u); markPath(u, b, v);
            for (let i = 0; i < n; i++) {
              if (inBlossom[base[i]]) { base[i] = b; if (!queued[i]) { queued[i] = 1; queue[t++] = i; } }
            }
          } else if (parent[u] < 0) {
            parent[u] = v;
            if (mate[u] < 0) {
              let x = u;
              while (x >= 0) { const pv = parent[x], nx = mate[pv]; mate[x] = pv; mate[pv] = x; x = nx; }
              return true;
            }
            const w = mate[u];
            if (!queued[w]) { queued[w] = 1; queue[t++] = w; }
          }
        }
      }
      return false;
    }
    if (order) { for (const v of order) if (mate[v] < 0 && adj[v].length) augmentFrom(v); }
    else for (let v = 0; v < n; v++) if (mate[v] < 0 && adj[v].length) augmentFrom(v);
    return mate;
  }
  /** Pick double bonds among candidate bonds so that each atom with need[i] gets exactly one.
   *  pairs: [[a, b], ...] for all bonds; cand: candidate bond indices (preferred first). */
  function kekuleMatch(n, pairs, cand, need) {
    const local = new Int32Array(n).fill(-1), verts = [];
    for (let i = 0; i < n; i++) if (need[i]) { local[i] = verts.length; verts.push(i); }
    const adj = verts.map(() => []), edge = new Map();
    for (const k of cand) {
      const a = local[pairs[k][0]], b = local[pairs[k][1]];
      if (a < 0 || b < 0) continue;
      adj[a].push(b); adj[b].push(a);
      edge.set(a < b ? a * verts.length + b : b * verts.length + a, k);
    }
    const mate = maxMatching(adj);
    const chosen = [], unmatched = [];
    for (let v = 0; v < verts.length; v++) {
      if (mate[v] < 0) unmatched.push(verts[v]);
      else if (v < mate[v]) chosen.push(edge.get(v * verts.length + mate[v]));
    }
    return { chosen, unmatched };
  }

  // =================================================================== SMILES parser
  const ORGANIC = { B: 5, C: 6, N: 7, O: 8, P: 15, S: 16, F: 9, I: 53 };
  const AROM_ORGANIC = { b: 5, c: 6, n: 7, o: 8, p: 15, s: 16 };
  const AROM_BRACKET = { b: 5, c: 6, n: 7, o: 8, p: 15, s: 16, se: 34, as: 33, te: 52 };
  const isDigit = ch => ch >= '0' && ch <= '9';

  function smilesError(msg, pos, src) {
    const e = new Error('SMILES syntax error at position ' + (pos + 1) + ': ' + msg + (src ? ' (in "' + src + '")' : ''));
    e.position = pos;
    return e;
  }
  function symbolZ(sym) { const i = EL().SYMBOLS.indexOf(sym); return i > 0 ? i : -1; }

  function parseBracket(s, i, end) {
    let j = i + 1;
    const d0 = j;
    while (j < end && isDigit(s[j])) j++;
    const isotope = j > d0 ? parseInt(s.slice(d0, j), 10) : 0;
    let Z = -1, aromatic = false;
    const ch = s[j];
    if (ch === '*') { Z = 0; j++; }
    else if (ch >= 'A' && ch <= 'Z') {
      const nx = s[j + 1];
      if (nx >= 'a' && nx <= 'z' && symbolZ(ch + nx) > 0) { Z = symbolZ(ch + nx); j += 2; }
      else if (symbolZ(ch) > 0) { Z = symbolZ(ch); j++; }
      else throw smilesError("unknown element '" + ch + (nx >= 'a' && nx <= 'z' ? nx : '') + "'", j, s);
    } else if (ch >= 'a' && ch <= 'z') {
      const two = s.slice(j, j + 2);
      if (two.length === 2 && AROM_BRACKET[two] !== undefined) { Z = AROM_BRACKET[two]; j += 2; }
      else if (AROM_BRACKET[ch] !== undefined) { Z = AROM_BRACKET[ch]; j++; }
      else throw smilesError("unknown aromatic symbol '" + ch + "'", j, s);
      aromatic = true;
    } else if (j >= end) throw smilesError("unterminated bracket atom '['", i, s);
    else throw smilesError("expected an element symbol after '[' but found '" + ch + "'", j, s);
    let chirality = null;
    if (s[j] === '@') {
      j++;
      if (s[j] === '@') { chirality = '@@'; j++; }
      else if (/^(TH|AL|SP|TB|OH)/.test(s.slice(j, j + 2))) {
        const cls = s.slice(j, j + 2);
        j += 2;
        const n0 = j;
        while (j < end && isDigit(s[j])) j++;
        if (j === n0) throw smilesError('chirality class @' + cls + ' needs a number', j, s);
        const num = parseInt(s.slice(n0, j), 10);
        chirality = cls === 'TH' && num === 1 ? '@' : cls === 'TH' && num === 2 ? '@@' : null; // other classes ignored
      } else chirality = '@';
    }
    let hCount = 0;
    if (s[j] === 'H') {
      j++;
      const n0 = j;
      while (j < end && isDigit(s[j])) j++;
      hCount = j > n0 ? parseInt(s.slice(n0, j), 10) : 1;
    }
    let charge = 0;
    if (s[j] === '+' || s[j] === '-') {
      const sym = s[j], sg = sym === '+' ? 1 : -1;
      j++;
      if (isDigit(s[j])) {
        const n0 = j;
        while (j < end && isDigit(s[j])) j++;
        charge = sg * parseInt(s.slice(n0, j), 10);
      } else { let k = 1; while (s[j] === sym) { k++; j++; } charge = sg * k; }
    }
    let atomClass = 0;
    if (s[j] === ':') {
      j++;
      const n0 = j;
      while (j < end && isDigit(s[j])) j++;
      if (j === n0) throw smilesError('atom class after ":" needs a number', j, s);
      atomClass = parseInt(s.slice(n0, j), 10);
    }
    if (j >= end) throw smilesError("unterminated bracket atom '['", i, s);
    if (s[j] !== ']') throw smilesError("unexpected '" + s[j] + "' inside bracket atom", j, s);
    return { atom: { Z, charge, hCount, aromatic, chirality, isotope, atomClass, bracket: true }, next: j + 1 };
  }

  /** Parse a SMILES string. Returns {atoms, bonds: [[i, j, order]] (Kekulé), stereo: {doubleBonds, tetrahedral}, aromaticBonds}. */
  function parseSmiles(smiles) {
    if (typeof smiles !== 'string') throw new Error('parseSmiles: expected a string');
    let start = 0;
    while (start < smiles.length && /\s/.test(smiles[start])) start++;
    let end = start;
    while (end < smiles.length && !/\s/.test(smiles[end])) end++;
    if (end === start) throw new Error('SMILES syntax error: empty string');
    const s = smiles, src = smiles.slice(start, end);
    const atoms = [], bonds = [], nbr = [], hasPrev = [];
    const branches = [], rings = new Map(), bondKeys = new Set();
    let prev = -1, pending = null, afterOpen = false, lastDot = -1;

    function addBond(a, b, sym, pos, dirFrom) {
      const key = Math.min(a, b) + ',' + Math.max(a, b);
      if (bondKeys.has(key)) throw smilesError('atoms ' + (Math.min(a, b) + 1) + ' and ' + (Math.max(a, b) + 1) + ' are bonded twice', pos, src);
      bondKeys.add(key);
      let order = 1, arom = false, dir = null;
      if (sym === '=') order = 2;
      else if (sym === '#') order = 3;
      else if (sym === '$') order = 4;
      else if (sym === ':') arom = true;
      else if (sym === '/' || sym === '\\') dir = sym;
      else if (sym === null) arom = atoms[a].aromatic && atoms[b].aromatic;
      bonds.push({ a, b, order, arom, dir, dirFrom, pos });
    }

    let i = start;
    while (i < end) {
      const c = s[i];
      if (c === '(') {
        if (prev < 0) throw smilesError("branch '(' without a preceding atom", i, src);
        if (afterOpen) throw smilesError("a branch must start with an atom or bond, not '('", i, src);
        if (pending) throw smilesError("bond symbol before '(' (write it inside the branch)", pending.pos, src);
        branches.push({ atom: prev, pos: i });
        afterOpen = true; i++; continue;
      }
      if (c === ')') {
        if (!branches.length) throw smilesError("unmatched ')'", i, src);
        if (pending) throw smilesError('bond symbol without a following atom', pending.pos, src);
        if (afterOpen) throw smilesError('empty branch "()"', i, src);
        if (prev < 0) throw smilesError("branch ends with '.'", i, src);
        prev = branches.pop().atom; i++; continue;
      }
      if (c === '-' || c === '=' || c === '#' || c === '$' || c === ':' || c === '/' || c === '\\') {
        if (pending) throw smilesError('two bond symbols in a row', i, src);
        if (prev < 0) throw smilesError("bond symbol '" + c + "' without a preceding atom", i, src);
        pending = { sym: c, pos: i }; i++; continue;
      }
      if (isDigit(c) || c === '%') {
        if (prev < 0) throw smilesError('ring-closure digit without a preceding atom', i, src);
        if (afterOpen) throw smilesError('ring-closure digit at the start of a branch', i, src);
        const p0 = i;
        let num;
        if (c === '%') {
          if (s[i + 1] === '(') {
            const close = s.indexOf(')', i + 2);
            if (close < 0 || close >= end || !/^\d+$/.test(s.slice(i + 2, close))) throw smilesError("'%(' must be followed by digits and ')'", i, src);
            num = parseInt(s.slice(i + 2, close), 10); i = close + 1;
          } else {
            if (!(isDigit(s[i + 1]) && isDigit(s[i + 2]) && i + 2 < end)) throw smilesError("'%' must be followed by two digits", i, src);
            num = parseInt(s.slice(i + 1, i + 3), 10); i += 3;
          }
        } else { num = c.charCodeAt(0) - 48; i++; }
        const open = rings.get(num);
        if (open) {
          if (open.atom === prev) throw smilesError('ring bond ' + num + ' opens and closes on the same atom', p0, src);
          const a = open.bond, b = pending;
          const isDir = x => x && (x.sym === '/' || x.sym === '\\');
          let use = a || b;
          if (a && b) {
            if (isDir(a) && isDir(b)) use = a;
            else if (isDir(a) && b.sym === '-') use = a;
            else if (isDir(b) && a.sym === '-') use = b;
            else if (a.sym !== b.sym) throw smilesError('conflicting bond symbols for ring bond ' + num, b.pos, src);
          }
          rings.delete(num);
          addBond(open.atom, prev, use ? use.sym : null, p0, use === a ? open.atom : prev);
          nbr[open.atom][open.slot] = prev;
          nbr[prev].push(open.atom);
        } else {
          rings.set(num, { atom: prev, bond: pending, slot: nbr[prev].length, pos: p0 });
          nbr[prev].push(null);
        }
        pending = null; continue;
      }
      if (c === '.') {
        if (prev < 0) throw smilesError("unexpected '.'", i, src);
        if (pending) throw smilesError("bond symbol before '.'", pending.pos, src);
        prev = -1; lastDot = i; i++; continue;
      }
      // ---- atoms
      let atom, p0 = i;
      if (c === '[') { const r = parseBracket(s, i, end); atom = r.atom; i = r.next; }
      else if (c === '*') { atom = { Z: 0, aromatic: false }; i++; }
      else if (c === 'C' && s[i + 1] === 'l') { atom = { Z: 17, aromatic: false }; i += 2; }
      else if (c === 'B' && s[i + 1] === 'r') { atom = { Z: 35, aromatic: false }; i += 2; }
      else if (ORGANIC[c] !== undefined) { atom = { Z: ORGANIC[c], aromatic: false }; i++; }
      else if (AROM_ORGANIC[c] !== undefined) { atom = { Z: AROM_ORGANIC[c], aromatic: true }; i++; }
      else if (c === 'H') throw smilesError("hydrogen must be written in brackets, e.g. [H] or [CH4]", i, src);
      else if (c === ']') throw smilesError("unmatched ']'", i, src);
      else throw smilesError("unexpected character '" + c + "'", i, src);
      if (!atom.bracket) Object.assign(atom, { charge: 0, hCount: 0, chirality: null, isotope: 0, atomClass: 0, bracket: false });
      const idx = atoms.length;
      atoms.push(atom); nbr.push([]); hasPrev.push(prev >= 0);
      if (prev >= 0) {
        addBond(prev, idx, pending ? pending.sym : null, pending ? pending.pos : p0, prev);
        nbr[prev].push(idx); nbr[idx].push(prev);
      }
      if (atom.bracket && atom.hCount > 0) nbr[idx].push(-1); // implicit H follows the preceding atom in stereo order
      pending = null; prev = idx; afterOpen = false;
    }
    if (pending) throw smilesError('bond symbol at the end of the SMILES', pending.pos, src);
    if (branches.length) throw smilesError("unclosed branch '('", branches[branches.length - 1].pos, src);
    if (rings.size) { const [num, r] = rings.entries().next().value; throw smilesError('ring bond ' + num + ' is never closed', r.pos, src); }
    if (prev < 0) throw smilesError("SMILES ends with '.'", lastDot >= 0 ? lastDot : end - 1, src);

    const n = atoms.length;
    // aromatic bonds must join aromatic atoms and lie on a ring
    const pairs = bonds.map(b => [b.a, b.b, b.order]);
    const ring = ringBondFlags(n, pairs);
    bonds.forEach((b, k) => { if (b.arom && !(atoms[b.a].aromatic && atoms[b.b].aromatic && ring[k])) b.arom = false; });
    const aromaticBonds = [];
    bonds.forEach((b, k) => { if (b.arom) aromaticBonds.push(k); });

    // ---- kekulization: which aromatic atoms still need a double bond?
    const need = new Uint8Array(n);
    const incident = Array.from({ length: n }, () => []);
    bonds.forEach(b => { incident[b.a].push(b); incident[b.b].push(b); });
    for (let a = 0; a < n; a++) {
      const at = atoms[a];
      if (!at.aromatic) continue;
      let bsum = 0, multi = false, nArom = 0;
      for (const b of incident[a]) {
        if (b.arom) { bsum += 1; nArom++; } else { bsum += b.order; if (b.order >= 2) multi = true; }
      }
      if (at.bracket) {
        bsum += at.hCount;
        const v = smallestAtLeast(valenceList(at.Z, at.charge), bsum);
        need[a] = !multi && v !== undefined && v - bsum >= 1 ? 1 : 0;
      } else {
        const v = smallestAtLeast(SMILES_VAL[at.Z], bsum);
        const r = v === undefined ? 0 : v - bsum;
        need[a] = !multi && r >= 1 ? 1 : 0;
        at.hCount = Math.max(0, r - need[a]);
      }
      if (need[a] && nArom === 0) {
        throw new Error('Cannot kekulize SMILES "' + src + '": aromatic atom ' + EL().symbol(at.Z).toLowerCase() + (a + 1) + ' is not part of an aromatic ring');
      }
    }
    const km = kekuleMatch(n, pairs, aromaticBonds, need);
    if (km.unmatched.length) {
      const list = km.unmatched.map(a => EL().symbol(atoms[a].Z).toLowerCase() + (a + 1)).join(', ');
      const hint = atoms.some((at, a) => need[a] && at.Z === 7 && !at.bracket && incident[a].length === 2)
        ? ' (a pyrrole-type nitrogen needs an explicit H, e.g. [nH])' : '';
      throw new Error('Cannot kekulize SMILES "' + src + '": no alternating single/double bond pattern exists for aromatic atoms ' + list + hint);
    }
    for (const k of km.chosen) bonds[k].order = 2;
    // ---- implicit hydrogens on non-aromatic organic-subset atoms
    for (let a = 0; a < n; a++) {
      const at = atoms[a];
      if (at.bracket || at.aromatic || at.Z === 0) continue;
      let bsum = 0;
      for (const b of incident[a]) bsum += b.order;
      const v = smallestAtLeast(SMILES_VAL[at.Z], bsum);
      at.hCount = v === undefined ? 0 : v - bsum;
    }
    // ---- stereo: double bonds from '/' and '\'
    const doubleBonds = [];
    for (const db of bonds) {
      if (db.order !== 2) continue;
      const side = center => {
        for (const d of incident[center]) {
          if (d === db || !d.dir) continue;
          const sub = d.a === center ? d.b : d.a;
          const up = d.dirFrom === sub ? (d.dir === '/' ? -1 : 1) : (d.dir === '/' ? 1 : -1);
          return { sub, up };
        }
        return null;
      };
      const sa = side(db.a), sb = side(db.b);
      if (sa && sb) doubleBonds.push({ bond: [db.a, db.b], ref: [sa.sub, sb.sub], config: sa.up === sb.up ? 'cis' : 'trans' });
    }
    // ---- stereo: tetrahedral centres with their neighbour order (-1 = implicit H or lone pair)
    const tetrahedral = [];
    for (let a = 0; a < n; a++) {
      const at = atoms[a];
      if (!at.chirality) continue;
      const list = nbr[a].slice();
      if (at.hCount > 1) continue;
      if (list.length === 3 && at.hCount === 0) list.splice(hasPrev[a] ? 1 : 0, 0, -1);
      if (list.length === 4) tetrahedral.push({ atom: a, chirality: at.chirality, neighbors: list });
    }
    return {
      atoms: atoms.map(at => ({ Z: at.Z, charge: at.charge, hCount: at.hCount, aromatic: at.aromatic, chirality: at.chirality,
        isotope: at.isotope, atomClass: at.atomClass, bracket: at.bracket })),
      bonds: bonds.map(b => [b.a, b.b, b.order]),
      stereo: { doubleBonds, tetrahedral },
      aromaticBonds
    };
  }

  // =================================================================== SMILES writer
  const ORGANIC_Z = new Set([5, 6, 7, 8, 9, 15, 16, 17, 35, 53]);
  function implicitHFor(z, bsum) { const v = smallestAtLeast(SMILES_VAL[z], bsum); return v === undefined ? 0 : Math.round(v - bsum); }

  /** Copy of bonds with 1.5 orders turned into an alternating Kekulé pattern (best effort; unmatched become single). */
  function kekulizeFractional(Z, bonds, charges, implicitH) {
    const n = Z.length;
    const out = bonds.map(b => [b[0], b[1], b[2] === undefined ? 1 : b[2]]);
    const frac = [];
    out.forEach((b, k) => { if (b[2] > 1 && b[2] < 2) frac.push(k); });
    if (!frac.length) return out;
    const need = new Uint8Array(n), sigma = new Float64Array(n), hasFrac = new Uint8Array(n), multi = new Uint8Array(n);
    for (const b of out) {
      for (const a of [b[0], b[1]]) {
        if (b[2] > 1 && b[2] < 2) { sigma[a] += 1; hasFrac[a] = 1; } else { sigma[a] += b[2]; if (b[2] >= 2) multi[a] = 1; }
      }
    }
    for (let a = 0; a < n; a++) {
      if (!hasFrac[a] || multi[a] || isMetal(Z[a])) continue;
      const s = sigma[a] + (implicitH ? implicitH[a] || 0 : 0);
      const v = smallestAtLeast(valenceList(Z[a], charges ? charges[a] || 0 : 0), s);
      if (v !== undefined && v - s >= 1) need[a] = 1;
    }
    const km = kekuleMatch(n, out, frac, need);
    for (const k of frac) out[k][2] = 1;
    for (const k of km.chosen) out[k][2] = 2;
    return out;
  }

  /** Write a (non-canonical) Kekulé SMILES for a molecule {Z, bonds, formalCharges?} or a parseSmiles() result. */
  function toSmiles(mol) {
    if (!mol) throw new Error('toSmiles: no molecule');
    let Z, charge, implicitH, iso, bonds;
    if (Array.isArray(mol.atoms)) {
      Z = mol.atoms.map(a => a.Z);
      charge = mol.atoms.map(a => a.charge || 0);
      implicitH = mol.atoms.map(a => a.hCount || 0);
      iso = mol.atoms.map(a => a.isotope || 0);
      bonds = mol.bonds || [];
    } else {
      Z = Array.from(mol.Z || []);
      const n0 = Z.length;
      charge = mol.formalCharges ? Array.from(mol.formalCharges, c => c || 0) : new Array(n0).fill(0);
      implicitH = mol.hCount ? Array.from(mol.hCount, h => h || 0) : new Array(n0).fill(0);
      iso = mol.isotopes ? Array.from(mol.isotopes, v => v || 0) : new Array(n0).fill(0);
      bonds = mol.bonds || [];
    }
    const n = Z.length;
    if (!n) return '';
    bonds = kekulizeFractional(Z, bonds, charge, implicitH).map(b => [b[0], b[1], Math.max(1, Math.min(4, Math.round(b[2])))]);
    const adj = adjOf(n, bonds);
    const folded = new Uint8Array(n), hTot = implicitH.slice();
    for (let i = 0; i < n; i++) {
      if (Z[i] !== 1 || charge[i] || iso[i] || adj[i].length !== 1) continue;
      const [j, o] = adj[i][0];
      if (o !== 1 || Z[j] === 1) continue;
      folded[i] = 1; hTot[j]++;
    }
    const bsum = new Float64Array(n);
    for (const [a, b, o] of bonds) { if (!folded[a] && !folded[b]) { bsum[a] += o; bsum[b] += o; } }
    const sym = z => (z > 0 && z < EL().SYMBOLS.length ? EL().SYMBOLS[z] : '*');
    function atomStr(i) {
      const z = Z[i], h = hTot[i], q = charge[i];
      if (z === 0 && !h && !q && !iso[i]) return '*';
      if (ORGANIC_Z.has(z) && !q && !iso[i] && implicitHFor(z, bsum[i]) === h) return sym(z);
      let s = '[' + (iso[i] ? iso[i] : '') + sym(z);
      if (h) s += 'H' + (h > 1 ? h : '');
      if (q) s += (q > 0 ? '+' : '-') + (Math.abs(q) > 1 ? Math.abs(q) : '');
      return s + ']';
    }
    // pass 1: DFS tree, children order and ring closures (back edges)
    const state = new Uint8Array(n), parent = new Int32Array(n).fill(-1);
    const children = Array.from({ length: n }, () => []);
    const ringOpen = Array.from({ length: n }, () => []), ringClose = Array.from({ length: n }, () => []);
    const roots = [];
    let ringId = 0;
    for (let r = 0; r < n; r++) {
      if (folded[r] || state[r]) continue;
      roots.push(r);
      state[r] = 1;
      const stack = [[r, 0]];
      while (stack.length) {
        const top = stack[stack.length - 1], u = top[0];
        if (top[1] < adj[u].length) {
          const [v, o] = adj[u][top[1]++];
          if (folded[v]) continue;
          if (!state[v]) { state[v] = 1; parent[v] = u; children[u].push([v, o]); stack.push([v, 0]); }
          else if (state[v] === 1 && v !== parent[u]) { const id = ringId++; ringOpen[v].push({ id, order: o }); ringClose[u].push({ id }); }
        } else { state[u] = 2; stack.pop(); }
      }
    }
    // pass 2: emit in the same pre-order
    const bondSym = o => (o === 2 ? '=' : o === 3 ? '#' : o === 4 ? '$' : '');
    const ringStr = k => (k < 10 ? String(k) : k < 100 ? '%' + k : '%(' + k + ')');
    const free = [], numOf = new Map();
    let next = 1;
    const take = () => { if (free.length) { free.sort((a, b) => a - b); return free.shift(); } return next++; };
    const out = [];
    roots.forEach((root, ci) => {
      if (ci) out.push('.');
      const stack = [{ u: root, s: '', branch: false }];
      while (stack.length) {
        const it = stack.pop();
        if (it.close) { out.push(')'); continue; }
        const u = it.u;
        if (it.branch) out.push('(');
        out.push(it.s + atomStr(u));
        const closing = [];
        for (const rc of ringClose[u]) { const k = numOf.get(rc.id); closing.push(k); out.push(ringStr(k)); }
        for (const ro of ringOpen[u]) { const k = take(); numOf.set(ro.id, k); out.push(bondSym(ro.order) + ringStr(k)); }
        for (const k of closing) free.push(k);
        const ch = children[u];
        if (ch.length) {
          stack.push({ u: ch[ch.length - 1][0], s: bondSym(ch[ch.length - 1][1]), branch: false });
          for (let c = ch.length - 2; c >= 0; c--) { stack.push({ close: true }); stack.push({ u: ch[c][0], s: bondSym(ch[c][1]), branch: true }); }
        }
      }
    });
    return out.join('');
  }

  // =================================================================== molecule context, valence, hybridization
  function ctxOf(mol) {
    const Z = mol.Z, n = Z.length, bonds = mol.bonds || [];
    const fc = mol.formalCharges;
    const ctx = { n, Z, xyz: mol.xyz, bonds, adj: adjOf(n, bonds), charge: i => (fc ? fc[i] || 0 : 0), aromatic: null, _ring: null };
    ctx.ringBond = k => { if (!ctx._ring) ctx._ring = ringBondFlags(n, bonds, ctx.adj); return ctx._ring[k] === 1; };
    return ctx;
  }
  /** Target valence for a non-metal with bond-order sum `sum` (rules described in valenceInfo). */
  function targetValence(z, q, sum, maxOrder) {
    const list = valenceList(z, q);
    if (!list) return sum;
    const lowest = list[0];
    if (sum <= lowest + 1e-9) return lowest;
    if (list.some(v => Math.abs(v - sum) < 1e-9)) return sum;
    const nx = smallestAtLeast(list, sum);
    return nx !== undefined && maxOrder >= 2 ? nx : sum; // hypervalent expansion only next to multiple bonds
  }
  function hasPi(ctx, j) {
    if (ctx.aromatic && ctx.aromatic[j]) return true;
    if (periodOf(ctx.Z[j]) !== 2) return false;
    for (const e of ctx.adj[j]) if (e[1] >= 1.5) return true;
    return false;
  }
  /** Valence bookkeeping and VSEPR hybridization for atom i (extraOrder: a hypothetical extra bond). */
  function valenceCore(ctx, i, extraOrder) {
    const z = ctx.Z[i], q = ctx.charge(i);
    let sum = 0, maxOrder = 0, deg = 0;
    for (const e of ctx.adj[i]) { sum += e[1]; deg++; if (e[1] > maxOrder) maxOrder = e[1]; }
    if (extraOrder) { sum += extraOrder; deg++; if (extraOrder > maxOrder) maxOrder = extraOrder; }
    const info = { element: EL().symbol(z), bondOrderSum: sum, targetValence: sum, missingH: 0, hybridization: 'none',
      charge: q, degree: deg, lonePairs: 0, stericNumber: deg };
    if (z === 1) { info.targetValence = q === 0 ? 1 : 0; return info; }
    if (!(z > 1) || NOBLE.has(z)) return info;
    if (isMetal(z)) { info.hybridization = 'metal'; return info; }
    const target = targetValence(z, q, sum, maxOrder);
    info.targetValence = target;
    info.missingH = Math.max(0, Math.floor(target - sum + 1e-6));
    const lp = Math.max(0, Math.floor((valenceElectrons(z, q) - sum - info.missingH) / 2 + 1e-9));
    const sigma = deg + info.missingH;
    info.lonePairs = lp;
    info.stericNumber = sigma + lp;
    let hyb = info.stericNumber <= 2 ? 'sp' : info.stericNumber === 3 ? 'sp2' : 'sp3';
    if (ctx.aromatic && ctx.aromatic[i]) hyb = 'sp2';
    else if (hyb === 'sp3' && lp > 0 && sigma <= 3 && info.stericNumber === 4) {
      // lone pair conjugated with a neighbouring pi system (amide/aniline/pyrrole N; furan-type ring O/S)
      let piNbrs = 0, ringPi = 0;
      for (const e of ctx.adj[i]) if (e[1] === 1 && hasPi(ctx, e[0])) { piNbrs++; if (ctx.ringBond(e[2])) ringPi++; }
      if ((z === 7 && piNbrs >= 1) || ((z === 8 || z === 16 || z === 34) && ringPi >= 2)) hyb = 'sp2';
    }
    info.hybridization = hyb;
    return info;
  }
  /** Public: {element, bondOrderSum, targetValence, missingH, hybridization, charge, degree, lonePairs, stericNumber}. */
  function valenceInfo(mol, i) {
    if (!(i >= 0 && i < mol.Z.length)) throw new RangeError('valenceInfo: atom index ' + i + ' out of range');
    return valenceCore(ctxOf(mol), i);
  }
  function geometryFor(info, z, cn) {
    if (info.hybridization === 'metal') {
      if (cn <= 4 && SQUARE_PLANAR.has(z)) return 'square';
      if (cn <= 4) return 'tetrahedral';
      if (cn <= 6) return 'octahedral';
      return 'generic';
    }
    if (info.hybridization === 'none') return 'linear';
    if (info.hybridization === 'sp') return 'linear';
    if (info.hybridization === 'sp2') return 'trigonal';
    return info.stericNumber <= 4 ? 'tetrahedral' : info.stericNumber === 5 ? 'generic' : 'octahedral';
  }

  // =================================================================== VSEPR sites
  const fib = (() => { // 240 quasi-uniform directions for generic placement
    const out = [], N = 240, ga = Math.PI * (3 - Math.sqrt(5));
    for (let k = 0; k < N; k++) { const y = 1 - 2 * (k + 0.5) / N, r = Math.sqrt(1 - y * y); out.push([Math.cos(ga * k) * r, y, Math.sin(ga * k) * r]); }
    return out;
  })();
  function genericSites(nb, count) {
    const used = nb.slice(), out = [];
    for (let c = 0; c < count; c++) {
      let best = null, bestScore = -Infinity;
      for (const d of fib) {
        let m = Infinity;
        for (const u of used) { const t = -vdot(d, u); if (t < m) m = t; } // maximise the smallest angle
        if (!used.length) m = d[0];
        if (m > bestScore + 1e-12) { bestScore = m; best = d; }
      }
      out.push(best); used.push(best);
    }
    return out;
  }
  /** Free directions for geometry `geom` around existing neighbour unit vectors `nb`.
   *  ref: direction (from the single neighbour to one of its substituents) used to stagger / pick the trans side. */
  function siteDirections(geom, nb, ref) {
    const k = nb.length, u1 = nb[0];
    const side = () => { // unit vector perpendicular to u1, pointing away from ref
      if (ref) { const p = vsub(ref, vscale(u1, vdot(ref, u1))); if (vnorm(p) > 1e-3) return vscale(vunit(p), -1); }
      return perpendicular(u1);
    };
    if (geom === 'linear') {
      if (k === 0) return [[1, 0, 0], [-1, 0, 0]];
      if (k === 1) return [vscale(u1, -1)];
      return [];
    }
    if (geom === 'trigonal') {
      if (k === 0) return [[1, 0, 0], [-0.5, Math.sqrt(3) / 2, 0], [-0.5, -Math.sqrt(3) / 2, 0]];
      if (k === 1) { const w = side(), c = -0.5, s = Math.sqrt(3) / 2; return [vadd(vscale(u1, c), vscale(w, s)), vadd(vscale(u1, c), vscale(w, -s))]; }
      if (k === 2) { const sm = vadd(nb[0], nb[1]); return vnorm(sm) < 0.1 ? [perpendicular(nb[0])] : [vscale(vunit(sm), -1)]; }
      return [];
    }
    if (geom === 'tetrahedral') {
      if (k === 0) { const a = 1 / Math.sqrt(3); return [[a, a, a], [a, -a, -a], [-a, a, -a], [-a, -a, a]]; }
      if (k === 1) {
        const ea = side(), eb = vcross(u1, ea), c = Math.cos(TET * RAD), s = Math.sin(TET * RAD), out = [];
        for (let m = 0; m < 3; m++) {
          const ph = m * 2 * Math.PI / 3;
          out.push(vadd(vscale(u1, c), vadd(vscale(ea, s * Math.cos(ph)), vscale(eb, s * Math.sin(ph)))));
        }
        return out;
      }
      if (k === 2) {
        const sm = vadd(nb[0], nb[1]);
        const b = vnorm(sm) < 0.1 ? perpendicular(nb[0]) : vscale(vunit(sm), -1);
        let nrm = vcross(nb[0], nb[1]);
        nrm = vnorm(nrm) < 1e-3 ? vunit(vcross(b, perpendicular(b))) : vunit(nrm);
        const h = TET / 2 * RAD;
        return [vadd(vscale(b, Math.cos(h)), vscale(nrm, Math.sin(h))), vadd(vscale(b, Math.cos(h)), vscale(nrm, -Math.sin(h)))];
      }
      if (k === 3) {
        const sm = vadd(vadd(nb[0], nb[1]), nb[2]);
        if (vnorm(sm) > 0.2) return [vscale(vunit(sm), -1)];
        return [vunit(vcross(vsub(nb[1], nb[0]), vsub(nb[2], nb[0])))];
      }
      return [];
    }
    if (geom === 'square' || geom === 'octahedral') {
      let e1 = k ? u1 : [1, 0, 0], e2 = null;
      for (let m = 1; m < k && !e2; m++) { const p = vsub(nb[m], vscale(e1, vdot(nb[m], e1))); if (vnorm(p) > 0.3) e2 = vunit(p); }
      if (!e2) e2 = perpendicular(e1);
      const e3 = vcross(e1, e2);
      const cand = geom === 'square' ? [e2, vscale(e1, -1), vscale(e2, -1), e1]
        : [e2, e3, vscale(e1, -1), vscale(e2, -1), vscale(e3, -1), e1];
      return cand.filter(d => nb.every(u => vdot(u, d) < 0.7));
    }
    return genericSites(nb, Math.max(1, 6 - k));
  }
  /** Reference vector for staggering when atom i has exactly one placed neighbour j. */
  function staggerRef(Z, xyz, nbrsOf, isPlaced, i, j) {
    let best = -1, bestKey = Infinity;
    for (const [s, o] of nbrsOf(j)) {
      if (s === i || !isPlaced(s)) continue;
      const key = (Z[s] === 1 ? 100 : 0) + o * 10 + s / 1e6;
      if (key < bestKey) { bestKey = key; best = s; }
    }
    return best < 0 ? null : vsub(P(xyz, best), P(xyz, j));
  }
  /** Positions for `count` new atoms bonded to i (bond length r), given placed neighbours. */
  function placeAround(Z, xyz, nbrsOf, isPlaced, i, count, geom, r) {
    const ci = P(xyz, i), placed = [], nb = [];
    for (const [j] of nbrsOf(i)) {
      if (!isPlaced(j)) continue;
      const d = vsub(P(xyz, j), ci);
      if (vnorm(d) > 1e-6) { placed.push(j); nb.push(vunit(d)); }
    }
    const ref = placed.length === 1 ? staggerRef(Z, xyz, nbrsOf, isPlaced, i, placed[0]) : null;
    let dirs = siteDirections(geom, nb, ref).slice(0, count);
    if (dirs.length < count) dirs = dirs.concat(genericSites(nb.concat(dirs), count - dirs.length));
    return dirs.map((d, m) => { const rr = Array.isArray(r) ? r[m] : r; return [ci[0] + d[0] * rr, ci[1] + d[1] * rr, ci[2] + d[2] * rr]; });
  }

  // =================================================================== building new molecule objects
  const PER_ATOM = ['formalCharges', 'atomNames', 'resNames', 'resIds', 'chainIds', 'chains', 'labels', 'isotopes', 'masses',
    'charges', 'spins', 'occupancy', 'bfactors', 'hCount', 'aromatic', 'atomTypes'];
  function isArrayLike(v) { return Array.isArray(v) || (ArrayBuffer.isView(v) && !(v instanceof DataView)); }
  /** New molecule from per-atom selection `keep` (old indices, in new order) plus `extra` appended atoms. */
  function rebuild(mol, keep, map, extraZ, extraXyz, bonds) {
    const out = Object.assign({}, mol);
    const n0 = mol.Z.length, nNew = keep.length + extraZ.length;
    out.Z = keep.map(i => mol.Z[i]).concat(extraZ);
    const xyz = new Float64Array(3 * nNew);
    const src = mol.xyz || new Float64Array(3 * n0);
    keep.forEach((i, k) => { xyz[3 * k] = src[3 * i]; xyz[3 * k + 1] = src[3 * i + 1]; xyz[3 * k + 2] = src[3 * i + 2]; });
    for (let m = 0; m < extraZ.length; m++) for (let d = 0; d < 3; d++) xyz[3 * (keep.length + m) + d] = extraXyz[3 * m + d];
    out.xyz = xyz;
    out.bonds = bonds;
    for (const key of PER_ATOM) {
      const v = mol[key];
      if (!isArrayLike(v) || v.length !== n0) continue;
      const fill = key === 'formalCharges' || ArrayBuffer.isView(v) ? 0 : null;
      const arr = keep.map(i => v[i]);
      for (let m = 0; m < extraZ.length; m++) arr.push(fill);
      out[key] = ArrayBuffer.isView(v) ? new v.constructor(arr) : arr;
    }
    return out;
  }
  function copyBonds(bonds) { return (bonds || []).map(b => [b[0], b[1], b[2] === undefined ? 1 : b[2]]); }
  function deleteMask(mol, mask) {
    const n = mol.Z.length, map = new Int32Array(n).fill(-1), keep = [];
    for (let i = 0; i < n; i++) if (!mask[i]) { map[i] = keep.length; keep.push(i); }
    const bonds = [];
    for (const b of mol.bonds || []) if (map[b[0]] >= 0 && map[b[1]] >= 0) bonds.push([map[b[0]], map[b[1]], b[2] === undefined ? 1 : b[2]]);
    return { mol: rebuild(mol, keep, map, [], [], bonds), map };
  }
  const identity = n => Array.from({ length: n }, (_, i) => i);

  // =================================================================== hydrogens
  /** Add hydrogens to satisfy valences (all atoms, or opts.atoms). Returns {mol, added}. */
  function addHydrogens(mol, opts) {
    opts = opts || {};
    const ctx = ctxOf(mol), n = ctx.n;
    const targets = opts.atoms ? [...new Set(opts.atoms)] : identity(n);
    const Z = Array.from(mol.Z), xyz = Array.from(mol.xyz || new Float64Array(3 * n));
    const nbrs = ctx.adj.map(l => l.map(e => [e[0], e[1]]));
    const bonds = copyBonds(mol.bonds), added = [];
    const nbrsOf = j => nbrs[j], isPlaced = () => true;
    for (const i of targets) {
      if (!(i >= 0 && i < n)) throw new RangeError('addHydrogens: atom index ' + i + ' out of range');
      const z = Z[i];
      if (!(z > 1) || isMetal(z) || NOBLE.has(z)) continue;
      const info = valenceCore(ctx, i);
      if (!info.missingH) continue;
      const pos = placeAround(Z, xyz, nbrsOf, isPlaced, i, info.missingH, geometryFor(info, z, info.degree + info.missingH), idealBondLength(z, 1, 1));
      for (const p of pos) {
        const k = Z.length;
        Z.push(1); xyz.push(p[0], p[1], p[2]); nbrs.push([[i, 1]]); nbrs[i].push([k, 1]);
        bonds.push([i, k, 1]); added.push(k);
      }
    }
    const extraZ = Z.slice(n), extraXyz = xyz.slice(3 * n);
    return { mol: rebuild(mol, identity(n), null, extraZ, extraXyz, bonds), added };
  }
  /** Remove hydrogens: all, or opts.atoms (H atoms themselves, or heavy atoms whose H are removed). Returns {mol, map}. */
  function removeHydrogens(mol, opts) {
    opts = opts || {};
    const n = mol.Z.length, mask = new Uint8Array(n);
    if (opts.atoms) {
      const adj = adjOf(n, mol.bonds || []);
      for (const a of opts.atoms) {
        if (!(a >= 0 && a < n)) throw new RangeError('removeHydrogens: atom index ' + a + ' out of range');
        if (mol.Z[a] === 1) mask[a] = 1; else for (const [j] of adj[a]) if (mol.Z[j] === 1) mask[j] = 1;
      }
    } else for (let i = 0; i < n; i++) if (mol.Z[i] === 1) mask[i] = 1;
    return deleteMask(mol, mask);
  }
  /** Make H counts match valences for opts.atoms (default all heavy atoms): wrong counts are removed and re-added. */
  function adjustHydrogens(mol, opts) {
    opts = opts || {};
    const ctx = ctxOf(mol), n = ctx.n;
    let targets = [];
    if (opts.atoms) {
      for (const a of opts.atoms) {
        if (!(a >= 0 && a < n)) throw new RangeError('adjustHydrogens: atom index ' + a + ' out of range');
        if (mol.Z[a] === 1) { for (const e of ctx.adj[a]) if (mol.Z[e[0]] !== 1) targets.push(e[0]); } else targets.push(a);
      }
      targets = [...new Set(targets)];
    } else targets = identity(n).filter(i => mol.Z[i] !== 1);
    const mask = new Uint8Array(n), redo = [];
    for (const i of targets) {
      const z = mol.Z[i];
      if (!(z > 1) || isMetal(z) || NOBLE.has(z)) continue;
      let heavySum = 0, maxOrder = 0;
      const hs = [];
      for (const e of ctx.adj[i]) {
        if (mol.Z[e[0]] === 1 && ctx.adj[e[0]].length === 1) hs.push(e[0]);
        else { heavySum += e[1]; if (e[1] > maxOrder) maxOrder = e[1]; }
      }
      const need = Math.max(0, Math.floor(targetValence(z, ctx.charge(i), heavySum, maxOrder) - heavySum + 1e-6));
      if (need !== hs.length) { for (const h of hs) mask[h] = 1; redo.push(i); }
    }
    const del = deleteMask(mol, mask);
    const res = addHydrogens(del.mol, { atoms: redo.map(i => del.map[i]) });
    return { mol: res.mol, map: del.map, added: res.added };
  }

  // =================================================================== editing helpers
  /** Best open VSEPR site for a new atom Znew bonded to atom i with the given order. */
  function newAtomPosition(mol, i, Znew, order) {
    if (!(i >= 0 && i < mol.Z.length)) throw new RangeError('newAtomPosition: atom index ' + i + ' out of range');
    order = order || 1;
    const ctx = ctxOf(mol), xyz = mol.xyz, ci = P(xyz, i);
    const r = idealBondLength(mol.Z[i], Znew, order);
    const placedNbrs = ctx.adj[i].filter(e => vnorm(vsub(P(xyz, e[0]), ci)) > 1e-6);
    if (!placedNbrs.length) return [ci[0] + r, ci[1], ci[2]];
    const info = valenceCore(ctx, i, order);
    const geom = geometryFor(info, mol.Z[i], info.degree + info.missingH);
    const nbrsOf = j => ctx.adj[j].map(e => [e[0], e[1]]);
    return placeAround(mol.Z, xyz, nbrsOf, () => true, i, 1, geom, r)[0];
  }
  function addAtom(mol, Z, xyz) {
    const n = mol.Z.length;
    const p = xyz ? [xyz[0], xyz[1], xyz[2]] : [0, 0, 0];
    return { mol: rebuild(mol, identity(n), null, [Z], p, copyBonds(mol.bonds)), index: n };
  }
  function deleteAtoms(mol, indices) {
    const n = mol.Z.length, mask = new Uint8Array(n);
    for (const i of indices || []) { if (!(i >= 0 && i < n)) throw new RangeError('deleteAtoms: atom index ' + i + ' out of range'); mask[i] = 1; }
    return deleteMask(mol, mask);
  }
  /** Set (or with order 0 remove) the bond i-j. Returns a new molecule. */
  function setBond(mol, i, j, order) {
    const n = mol.Z.length;
    if (!(i >= 0 && i < n && j >= 0 && j < n)) throw new RangeError('setBond: atom index out of range');
    if (i === j) throw new Error('setBond: an atom cannot be bonded to itself');
    if (!(order >= 0)) throw new Error('setBond: order must be 0 (delete) or positive');
    const bonds = [];
    let found = false;
    for (const b of mol.bonds || []) {
      if ((b[0] === i && b[1] === j) || (b[0] === j && b[1] === i)) { found = true; if (order > 0) bonds.push([b[0], b[1], order]); }
      else bonds.push([b[0], b[1], b[2] === undefined ? 1 : b[2]]);
    }
    if (!found && order > 0) bonds.push([i, j, order]);
    const out = Object.assign({}, mol, { Z: Array.from(mol.Z), bonds });
    if (mol.xyz) out.xyz = Float64Array.from(mol.xyz);
    if (mol.formalCharges) out.formalCharges = Array.from(mol.formalCharges);
    return out;
  }
  function fragments(mol) { return fragmentsOf(mol.Z.length, mol.bonds || []); }

  // =================================================================== bond-order perception
  function angleAt(xyz, a, j, b) {
    const u = vsub(P(xyz, a), P(xyz, j)), v = vsub(P(xyz, b), P(xyz, j));
    const c = vdot(u, v) / ((vnorm(u) * vnorm(v)) || 1);
    return Math.acos(Math.max(-1, Math.min(1, c))) * DEG;
  }
  /** Assign bond orders from connectivity, valence and geometry. Returns [[i, j, order]] in the input bond order.
   *  Bonds to metals stay single. With opts.withCharges returns {bonds, formalCharges} (charges implied by the
   *  chosen structure when mol.formalCharges is absent). */
  function perceiveBondOrders(mol, opts) {
    opts = opts || {};
    const Z = mol.Z, n = Z.length, xyz = mol.xyz;
    let pairs;
    if (mol.bonds && mol.bonds.length) pairs = mol.bonds.map(b => [b[0], b[1]]);
    else if (xyz && BM.geom && BM.geom.computeBonds) {
      const u = BM.geom.computeBonds(Z, xyz);
      pairs = [];
      for (let k = 0; k < u.length; k += 2) pairs.push([u[k], u[k + 1]]);
    } else pairs = [];
    const given = !!mol.formalCharges;
    const q = i => (given ? mol.formalCharges[i] || 0 : 0);
    const adj = adjOf(n, pairs.map(p => [p[0], p[1], 1]));
    const deg = adj.map(l => l.length);
    const noPi = i => Z[i] <= 1 || isMetal(Z[i]) || NOBLE.has(Z[i]);
    // pi capacity of every atom from its valence
    const cap = new Int32Array(n), extra = new Int32Array(n);
    for (let i = 0; i < n; i++) {
      if (noPi(i)) continue;
      const list = valenceList(Z[i], q(i));
      let v = smallestAtLeast(list, deg[i]);
      if (v === undefined) continue;
      let c = v - deg[i];
      if (periodOf(Z[i]) >= 3) { // S(=O)2, P=O, ClO4: expand while terminal O/S/N neighbours remain
        let term = 0;
        for (const [j] of adj[i]) if (deg[j] === 1 && [7, 8, 16, 34].includes(Z[j])) term++;
        for (const w of list) { if (c >= term) break; if (w > v) { v = w; c = v - deg[i]; } }
      }
      if (xyz && c >= 2 && deg[i] === 2 && periodOf(Z[i]) === 2 && angleAt(xyz, adj[i][0][0], i, adj[i][1][0]) < 155) c = 1;
      cap[i] = Math.max(0, c);
      // charge-separated alternatives (N+ / O+), used only if they add pi bonds: nitro, N-oxides, azides, CO, onium rings
      if (!given && (Z[i] === 7 || Z[i] === 8)) extra[i] = Math.max(0, (Z[i] === 7 ? 4 : 3) - deg[i] - cap[i]);
    }
    if (xyz && !opts.noGeometryCaps) {
      for (let i = 0; i < n; i++) if (extra[i] && deg[i] === 2 && angleAt(xyz, adj[i][0][0], i, adj[i][1][0]) < 155 && cap[i] >= 1) extra[i] = 0;
    }
    // candidate pi bonds, most double-like first
    const cand = [];
    pairs.forEach(([a, b], k) => {
      if (noPi(a) || noPi(b) || !(cap[a] + extra[a]) || !(cap[b] + extra[b])) return;
      let score = 0.5;
      if (xyz) {
        const d = vnorm(vsub(P(xyz, a), P(xyz, b)));
        const r1 = idealBondLength(Z[a], Z[b], 1), r2 = idealBondLength(Z[a], Z[b], 2);
        if (d > r1 + 0.12) return;
        score = (d - r2) / (r1 - r2);
      }
      cand.push({ k, a, b, score });
    });
    cand.sort((x, y) => x.score - y.score || x.k - y.k);
    // Split every atom into one vertex per pi bond it may form; promotion copies (N+ / O+) come last.
    function solve(ext, init) {
      const first = new Int32Array(n);
      let nv = 0;
      for (let i = 0; i < n; i++) { first[i] = nv; nv += cap[i] + ext[i]; }
      const owner = new Int32Array(nv), promo = new Uint8Array(nv);
      for (let i = 0; i < n; i++) for (let c = 0; c < cap[i] + ext[i]; c++) { owner[first[i] + c] = i; promo[first[i] + c] = c >= cap[i] ? 1 : 0; }
      const vadj = Array.from({ length: nv }, () => []), edgeBond = new Map();
      for (const { k, a, b } of cand) {
        for (let ca = 0; ca < cap[a] + ext[a]; ca++) for (let cb = 0; cb < cap[b] + ext[b]; cb++) {
          const u = first[a] + ca, v = first[b] + cb;
          vadj[u].push(v); vadj[v].push(u); edgeBond.set(u * nv + v, k); edgeBond.set(v * nv + u, k);
        }
      }
      const mate = new Int32Array(nv).fill(-1);
      if (init) { // carry over pi bonds of a previous solution (normal copies only)
        const used = new Int32Array(n);
        for (const [a, b, times] of init) for (let t = 0; t < times; t++) {
          if (used[a] < cap[a] && used[b] < cap[b]) { const u = first[a] + used[a]++, v = first[b] + used[b]++; mate[u] = v; mate[v] = u; }
        }
      } else { // greedy: most double-like bonds first; a second pass forms triple bonds where lengths say so
        const freeCopy = i => { for (let c = 0; c < cap[i]; c++) if (mate[first[i] + c] < 0) return first[i] + c; return -1; };
        for (let pass = 0; pass < 2; pass++) {
          for (const { a, b, score } of cand) {
            if (pass === 1 && xyz && score > -0.4) continue;
            const u = freeCopy(a), v = freeCopy(b);
            if (u >= 0 && v >= 0) { mate[u] = v; mate[v] = u; }
          }
        }
      }
      // augment from unsatisfied C/N sites first, then promotion copies, so negative charges end on O/S
      const rank = u => (promo[u] ? 1 : [8, 16, 34].includes(Z[owner[u]]) ? 2 : 0);
      const order = Array.from({ length: nv }, (_, u) => u).sort((p, q) => rank(p) - rank(q) || p - q);
      const fin = maxMatching(vadj, mate, order);
      const piCount = new Map(), promoted = new Uint8Array(n), unsat = new Uint8Array(n);
      for (let u = 0; u < nv; u++) {
        const v = fin[u];
        if (v > u) { const k = edgeBond.get(u * nv + v); piCount.set(k, (piCount.get(k) || 0) + 1); }
        if (v >= 0 && promo[u]) promoted[owner[u]] = 1;
        if (v < 0 && !promo[u]) unsat[owner[u]] = 1;
      }
      return { piCount, promoted, unsat };
    }
    let sol = solve(new Int32Array(n), null);
    if (!given) {
      // charge separation only next to an atom that is still unsatisfied, and only if it ends up compensated
      const ext = new Int32Array(n);
      for (let i = 0; i < n; i++) if (extra[i] && adj[i].some(e => sol.unsat[e[0]])) ext[i] = extra[i];
      for (let round = 0; round < 4 && ext.some(v => v); round++) {
        const init = [...sol.piCount].map(([k, t]) => [pairs[k][0], pairs[k][1], t]);
        const trial = solve(ext, init);
        let bad = false;
        for (let i = 0; i < n; i++) if (trial.promoted[i] && !adj[i].some(e => trial.unsat[e[0]])) { ext[i] = 0; bad = true; }
        if (!bad) { sol = trial; break; }
      }
    }
    const orders = pairs.map(() => 1);
    const charges = new Array(n).fill(0);
    for (const [k, t] of sol.piCount) orders[k] += t;
    for (let i = 0; i < n; i++) if (sol.promoted[i]) charges[i] = 1;
    const out = pairs.map(([a, b], k) => [a, b, orders[k]]);
    if (!opts.withCharges) return out;
    if (given) return { bonds: out, formalCharges: Array.from(mol.formalCharges, c => c || 0) };
    for (let i = 0; i < n; i++) {
      if (charges[i] || noPi(i)) continue;
      let s = 0;
      for (const e of adj[i]) s += orders[e[2]];
      if (Z[i] === 8 && s === 1 && deg[i] === 1) charges[i] = -1;
      if (Z[i] === 7 && s === 2 && deg[i] <= 2 && cap[i] > 0) charges[i] = -1;
    }
    return { bonds: out, formalCharges: charges };
  }

  // =================================================================== 3D structure generation
  const SP3_RING = { 3: 60, 4: 88, 5: 104, 6: 111, 7: 114, 8: 116 };
  const SP2_RING = { 3: 60, 4: 90, 5: 108, 6: 120, 7: 124, 8: 126 };
  const W12 = 10, W13 = 3, W14P = 2, W14 = 0.5, W_NB = 0.6, W_CHIRAL = 2, W_PLANAR = 3, W_LIN = 10, W_TOR = 0.3, NB_SCALE = 0.75;

  /** Kekulé-written benzenoid six-rings (each atom has one double bond, to a ring atom) are treated as aromatic. */
  function markBenzenoid(n, bonds, aromB) {
    const adj = adjOf(n, bonds), ring = ringBondFlags(n, bonds, adj);
    const pt = pathTool(adj.map(l => l.map(e => e[0])));
    const dblCount = new Uint8Array(n), dblRing = new Uint8Array(n), bondOf = new Map();
    bonds.forEach((b, k) => {
      bondOf.set(Math.min(b[0], b[1]) * n + Math.max(b[0], b[1]), k);
      if (b[2] === 2) { dblCount[b[0]]++; dblCount[b[1]]++; if (ring[k]) { dblRing[b[0]] = 1; dblRing[b[1]] = 1; } }
    });
    const done = new Set();
    bonds.forEach((b, k) => {
      if (!ring[k]) return;
      const path = pt.path(b[0], b[1], -1, -1, 5, true);
      if (!path || path.length !== 6) return;
      const key = path.slice().sort((x, y) => x - y).join(',');
      if (done.has(key)) return;
      done.add(key);
      if (!path.every(a => dblCount[a] === 1 && dblRing[a])) return;
      const ks = [];
      for (let t = 0; t < 6; t++) {
        const a = path[t], c = path[(t + 1) % 6], kk = bondOf.get(Math.min(a, c) * n + Math.max(a, c));
        if (kk === undefined) return;
        ks.push(kk);
      }
      for (const kk of ks) aromB[kk] = 1;
    });
  }

  function polyAngles(z, m, hyb, info) {
    const metal = hyb === 'metal';
    if (!metal && m < 5 && !(m === 4 && info.stericNumber === 6)) return null;
    // square planar (also partially built: two or three ligands of Pt/Pd/Au/Rh/Ir sit on square sites)
    if ((m === 4 && (metal ? SQUARE_PLANAR.has(z) : info.stericNumber === 6)) || (metal && m < 4 && SQUARE_PLANAR.has(z))) {
      return (p, q) => ((q - p) % 2 === 0 ? 180 : 90);
    }
    if (m === 6) { const opp = [5, 3, 4, 1, 2, 0]; return (p, q) => (opp[p] === q ? 180 : 90); }
    if (m === 5) return (p, q) => { const ax = x => x === 0 || x === 4; return ax(p) && ax(q) ? 180 : ax(p) || ax(q) ? 90 : 120; };
    if (!metal) return null;
    const g = groupOf(z);
    if (m === 2) return () => (g === 11 || g === 12 ? 180 : TET);
    if (m === 3) return () => (g === 11 || g === 12 ? 120 : TET);
    if (m === 4) return () => TET;
    return () => null; // high coordination numbers: only repulsion between ligands
  }

  /** Topology-derived distance bounds, chirality and planarity data for a whole molecule with explicit H. */
  function buildModel(Z, bonds, charges, aromB, aromA, stereo) {
    const n = Z.length, warnings = [];
    const ctx = ctxOf({ Z, bonds, formalCharges: charges });
    ctx.aromatic = aromA;
    const info = [], hyb = [];
    for (let i = 0; i < n; i++) { info.push(valenceCore(ctx, i)); hyb.push(info[i].hybridization); }
    const nbrs = ctx.adj.map(l => l.map(e => e[0])), adjO = ctx.adj.map(l => l.map(e => [e[0], e[1]]));
    const pt = pathTool(nbrs);
    const K = (a, b) => (a < b ? a * n + b : b * n + a);
    const bondOf = new Map();
    const blen = bonds.map((b, k) => { bondOf.set(K(b[0], b[1]), k); return idealBondLength(Z[b[0]], Z[b[1]], aromB[k] ? 1.5 : b[2]); });
    const len = (a, b) => blen[bondOf.get(K(a, b))];
    const termH = Z.map((z, i) => z === 1 && nbrs[i].length === 1 && Z[nbrs[i][0]] !== 1);
    const planarHyb = i => hyb[i] === 'sp2' || hyb[i] === 'sp';
    const ringFlag = ringBondFlags(n, bonds, ctx.adj);
    const pairs = new Map(); // K -> {a, b, lo, hi, w, kind (1: 1-2, 2: 1-3, 3: 1-4)}
    function putPair(a, b, lo, hi, w, kind) {
      const key = K(a, b), old = pairs.get(key);
      if (!old || old.kind > kind) { pairs.set(key, { a, b, lo, hi, w, kind }); return; }
      if (old.kind < kind) return;
      let nlo = Math.max(old.lo, lo), nhi = Math.min(old.hi, hi);
      if (nlo > nhi) { const mid = (nlo + nhi) / 2; nlo = mid - 0.01; nhi = mid + 0.01; }
      old.lo = nlo; old.hi = nhi; old.w = Math.max(old.w, w);
    }
    bonds.forEach((b, k) => putPair(b[0], b[1], blen[k] - 0.01, blen[k] + 0.01, W12, 1));
    // ---- 1-3: bond angles
    const ang = new Map(), A = (j, a, b) => j * n * n + K(a, b), linear = [], pyramidal = new Uint8Array(n);
    for (let j = 0; j < n; j++) {
      const nb = nbrs[j], m = nb.length;
      if (m < 2) continue;
      const poly = polyAngles(Z[j], m, hyb[j], info[j]);
      const list = [];
      for (let p = 0; p < m; p++) {
        for (let q = p + 1; q < m; q++) {
          const a = nb[p], b = nb[q];
          let t, tl = 2.5, ring = false;
          if (poly) { t = poly(p, q); tl = 3; }
          else {
            const path = pt.path(a, b, j, -1, 6, false);
            if (path) {
              ring = true;
              const size = path.length + 1;
              const flat = planarHyb(j) && path.every(planarHyb);
              t = flat ? 180 * (size - 2) / size : hyb[j] === 'sp3' ? SP3_RING[size] : SP2_RING[size];
            } else { t = hyb[j] === 'sp' ? 180 : hyb[j] === 'sp2' ? 120 : TET; tl = hyb[j] === 'sp3' ? 4 : 2.5; }
          }
          list.push({ a, b, t, tl, ring });
        }
      }
      if (!poly && hyb[j] === 'sp2' && m === 3) {
        let s = 0;
        const free = list.filter(e => !e.ring);
        for (const e of list) if (e.ring) s += e.t;
        if (free.length && free.length < 3) { const t = Math.max(100, Math.min(160, (360 - s) / free.length)); for (const e of free) e.t = t; }
        if (!free.length && s < 355) pyramidal[j] = 1; // e.g. fullerene / corannulene: ring angles cannot sum to 360
      } else if (!poly && hyb[j] === 'sp3') {
        let mn = Infinity;
        for (const e of list) if (e.ring) mn = Math.min(mn, e.t);
        if (mn < TET - 1) for (const e of list) if (!e.ring) e.t = TET + (TET - mn) * 0.14;
      }
      for (const e of list) {
        if (e.t === null || e.t === undefined) continue;
        ang.set(A(j, e.a, e.b), e.t);
        if (e.t >= 179) linear.push([e.a, j, e.b]);
        if (bondOf.has(K(e.a, e.b))) continue; // three-membered ring
        const r1 = len(j, e.a), r2 = len(j, e.b);
        const d = t => Math.sqrt(Math.max(0, r1 * r1 + r2 * r2 - 2 * r1 * r2 * Math.cos(Math.min(180, Math.max(0, t)) * RAD)));
        putPair(e.a, e.b, d(e.t - e.tl), d(e.t + e.tl), W13, 2);
      }
    }
    // ---- 1-4: torsions (cis/trans for planar bonds, a range for rotatable ones)
    const dbStereo = new Map();
    for (const s of (stereo && stereo.doubleBonds) || []) dbStereo.set(K(s.bond[0], s.bond[1]), s);
    const planarBonds = [], staggered = [];
    bonds.forEach((bd, k) => {
      const j = bd[0], kk = bd[1];
      const J = nbrs[j].filter(x => x !== kk), KN = nbrs[kk].filter(x => x !== j);
      if (!J.length || !KN.length) return;
      const planar = hyb[j] === 'sp2' && hyb[kk] === 'sp2' && !pyramidal[j] && !pyramidal[kk] &&
        (aromB[k] || bd[2] >= 2 || (ringFlag[k] && pt.path(j, kk, -1, -1, 7, true) !== null));
      let ref = null;
      if (planar) {
        for (const i of J) {
          for (const l of KN) if (i !== l && pt.path(i, l, j, kk, 5, false)) { ref = { x: i, y: l, cis: true, ring: true }; break; }
          if (ref) break;
        }
        const st = dbStereo.get(K(j, kk));
        if (st) {
          const fwd = st.bond[0] === j, x = fwd ? st.ref[0] : st.ref[1], y = fwd ? st.ref[1] : st.ref[0], want = st.config === 'cis';
          if (ref) {
            const implied = (ref.x === x) === (ref.y === y) ? ref.cis : !ref.cis;
            if (implied !== want) warnings.push('Double-bond stereo for atoms ' + (j + 1) + '=' + (kk + 1) + ' conflicts with its small ring and was ignored.');
          } else ref = { x, y, cis: want };
        }
        if (!ref) {
          const hx = J.find(i => !termH[i]), hy = KN.find(l => !termH[l]);
          ref = { x: hx === undefined ? J[0] : hx, y: hy === undefined ? KN[0] : hy, cis: false };
        }
        planarBonds.push({ j, k: kk, ref });
      } else if (hyb[j] === 'sp3' && hyb[kk] === 'sp3' && bd[2] === 1 && !(ringFlag[k] && pt.path(j, kk, -1, -1, 4, true))) {
        staggered.push([j, kk]); // 3-fold torsion preference (chairs, staggered chains); not in 3-5 rings
      }
      for (const i of J) {
        for (const l of KN) {
          if (i === l || bondOf.has(K(i, l))) continue;
          const t1 = ang.get(A(j, i, kk)), t2 = ang.get(A(kk, j, l));
          if (t1 === undefined || t2 === undefined) continue;
          const r1 = len(i, j), r2 = blen[k], r3 = len(kk, l);
          const c1 = Math.cos(t1 * RAD), s1 = Math.sin(t1 * RAD), c2 = Math.cos(t2 * RAD), s2 = Math.sin(t2 * RAD);
          const dAt = cp => { const dx = r2 - r3 * c2 - r1 * c1, dy = r3 * s2 * cp - r1 * s1, dz = r3 * s2 * Math.sqrt(Math.max(0, 1 - cp * cp)); return Math.sqrt(dx * dx + dy * dy + dz * dz); };
          const dc = dAt(1), dt = dAt(-1);
          if (planar) {
            const cis = (i === ref.x) === (l === ref.y) ? ref.cis : !ref.cis;
            const d = cis ? dc : dt;
            putPair(i, l, d - 0.04, d + 0.04, W14P, 3);
          } else putPair(i, l, Math.min(dc, dt), Math.max(dc, dt), W14, 3);
        }
      }
    });
    // ---- tetrahedral chirality as signed volumes over three non-H neighbours
    const chiral = [];
    for (const t of (stereo && stereo.tetrahedral) || []) {
      const c = t.atom;
      const nl = t.neighbors.map(x => (x < 0 || termH[x] ? -1 : x));
      const holes = nl.filter(x => x < 0).length;
      if (holes > 1 || nl.length !== 4) continue;
      const k = holes ? nl.indexOf(-1) : 0;
      const tri = nl.filter((x, idx) => idx !== k);
      const sign = (t.chirality === '@@' ? 1 : -1) * (k % 2 === 0 ? 1 : -1);
      const v0 = 0.7698 * len(c, tri[0]) * len(c, tri[1]) * len(c, tri[2]);
      chiral.push({ c, a: tri[0], b: tri[1], d: tri[2], sign, vlo: 0.45 * v0 });
    }
    return { n, Z, bonds, blen, aromB, hyb, info, nbrs, adjO, termH, pairs, chiral, planarBonds, staggered, linear, pyramidal, warnings, charges };
  }

  /** Packed error-function terms for a subset of atoms (local indices). */
  function buildTerms(model, atoms) {
    const n = model.n, m = atoms.length, loc = new Int32Array(n).fill(-1);
    atoms.forEach((a, k) => { loc[a] = k; });
    const bi = [], bj = [], blo = [], bhi = [], bw = [], kind = [], seen = new Set();
    for (const p of model.pairs.values()) {
      const u = loc[p.a], v = loc[p.b];
      if (u < 0 || v < 0) continue;
      bi.push(u); bj.push(v); blo.push(p.lo); bhi.push(p.hi); bw.push(p.w); kind.push(p.kind);
      seen.add(u < v ? u * m + v : v * m + u);
    }
    const vdw = EL().vdw, ni = [], nj = [], nlo = [];
    for (let u = 0; u < m; u++) {
      const ru = vdw[model.Z[atoms[u]]] || 1.5;
      for (let v = u + 1; v < m; v++) {
        if (seen.has(u * m + v)) continue;
        ni.push(u); nj.push(v); nlo.push(NB_SCALE * (ru + (vdw[model.Z[atoms[v]]] || 1.5)));
      }
    }
    const cc = [], ca = [], cb = [], cd = [], cs = [], cv = [];
    for (const t of model.chiral) {
      if (loc[t.c] < 0 || loc[t.a] < 0 || loc[t.b] < 0 || loc[t.d] < 0) continue;
      cc.push(loc[t.c]); ca.push(loc[t.a]); cb.push(loc[t.b]); cd.push(loc[t.d]); cs.push(t.sign); cv.push(t.vlo);
    }
    const pc = [], pa = [], pb = [], pd = [];
    for (const i of atoms) {
      const nb = model.nbrs[i];
      if (model.hyb[i] !== 'sp2' || model.pyramidal[i] || nb.length !== 3 || nb.some(j => loc[j] < 0)) continue;
      pc.push(loc[i]); pa.push(loc[nb[0]]); pb.push(loc[nb[1]]); pd.push(loc[nb[2]]);
    }
    // coplanarity of i-j-k-l across every planar (double / aromatic / conjugated ring) bond
    for (const pbd of model.planarBonds) {
      const j = pbd.j, k = pbd.k;
      if (loc[j] < 0 || loc[k] < 0) continue;
      for (const i of model.nbrs[j]) {
        if (i === k || loc[i] < 0) continue;
        for (const l of model.nbrs[k]) {
          if (l === j || l === i || loc[l] < 0) continue;
          pc.push(loc[j]); pa.push(loc[i]); pb.push(loc[k]); pd.push(loc[l]);
        }
      }
    }
    const la = [], lj = [], lb = [];
    for (const [a, j, b] of model.linear) if (loc[a] >= 0 && loc[j] >= 0 && loc[b] >= 0) { la.push(loc[a]); lj.push(loc[j]); lb.push(loc[b]); }
    const ti = [], tj = [], tk = [], tl = [], tw = [];
    for (const [j, k] of model.staggered) {
      if (loc[j] < 0 || loc[k] < 0) continue;
      const I1 = model.nbrs[j].filter(i => i !== k && loc[i] >= 0), L1 = model.nbrs[k].filter(l => l !== j && loc[l] >= 0);
      if (!I1.length || !L1.length) continue;
      const w = W_TOR / (I1.length * L1.length) * 3; // about W_TOR per bond with full substitution
      for (const i of I1) for (const l of L1) if (i !== l) { ti.push(loc[i]); tj.push(loc[j]); tk.push(loc[k]); tl.push(loc[l]); tw.push(w); }
    }
    const I = a => Int32Array.from(a), F = a => Float64Array.from(a);
    return { m, atoms, bi: I(bi), bj: I(bj), blo: F(blo), bhi: F(bhi), bw: F(bw), kind: Uint8Array.from(kind),
      ni: I(ni), nj: I(nj), nlo: F(nlo), cc: I(cc), ca: I(ca), cb: I(cb), cd: I(cd), cs: F(cs), cv: F(cv),
      pc: I(pc), pa: I(pa), pb: I(pb), pd: I(pd), la: I(la), lj: I(lj), lb: I(lb),
      ti: I(ti), tj: I(tj), tk: I(tk), tl: I(tl), tw: F(tw) };
  }

  /** Error function and gradient (D = 3 or 4). wPlanar: improper weight (3D only); w4: 4th-dimension squash. */
  function evalTerms(T, x, g, D, wPlanar, w4, wTor) {
    g.fill(0);
    let E = 0;
    const bi = T.bi, bj = T.bj, blo = T.blo, bhi = T.bhi, bw = T.bw;
    for (let t = 0, nb = bi.length; t < nb; t++) {
      const i = bi[t] * D, j = bj[t] * D;
      const dx = x[i] - x[j], dy = x[i + 1] - x[j + 1], dz = x[i + 2] - x[j + 2], dw = D === 4 ? x[i + 3] - x[j + 3] : 0;
      const d2 = dx * dx + dy * dy + dz * dz + dw * dw, lo = blo[t], hi = bhi[t];
      if (d2 >= lo * lo && d2 <= hi * hi) continue;
      const d = Math.sqrt(d2) + 1e-12, v = d < lo ? d - lo : d - hi, f = 2 * bw[t] * v / d;
      E += bw[t] * v * v;
      g[i] += f * dx; g[i + 1] += f * dy; g[i + 2] += f * dz; g[j] -= f * dx; g[j + 1] -= f * dy; g[j + 2] -= f * dz;
      if (D === 4) { g[i + 3] += f * dw; g[j + 3] -= f * dw; }
    }
    const ni = T.ni, nj = T.nj, nlo = T.nlo;
    for (let t = 0, nn = ni.length; t < nn; t++) {
      const i = ni[t] * D, j = nj[t] * D;
      const dx = x[i] - x[j], dy = x[i + 1] - x[j + 1], dz = x[i + 2] - x[j + 2], dw = D === 4 ? x[i + 3] - x[j + 3] : 0;
      const d2 = dx * dx + dy * dy + dz * dz + dw * dw, lo = nlo[t];
      if (d2 >= lo * lo) continue;
      const d = Math.sqrt(d2) + 1e-12, v = d - lo, f = 2 * W_NB * v / d;
      E += W_NB * v * v;
      g[i] += f * dx; g[i + 1] += f * dy; g[i + 2] += f * dz; g[j] -= f * dx; g[j + 1] -= f * dy; g[j + 2] -= f * dz;
      if (D === 4) { g[i + 3] += f * dw; g[j + 3] -= f * dw; }
    }
    const vol = (c, a, b, e, wf) => { // wf(V) -> [energy, dE/dV] or null
      const ax = x[a] - x[c], ay = x[a + 1] - x[c + 1], az = x[a + 2] - x[c + 2];
      const bx = x[b] - x[c], by = x[b + 1] - x[c + 1], bz = x[b + 2] - x[c + 2];
      const ex = x[e] - x[c], ey = x[e + 1] - x[c + 1], ez = x[e + 2] - x[c + 2];
      const kx = by * ez - bz * ey, ky = bz * ex - bx * ez, kz = bx * ey - by * ex;
      const r = wf(ax * kx + ay * ky + az * kz);
      if (!r) return;
      E += r[0];
      const f = r[1];
      const lx = ey * az - ez * ay, ly = ez * ax - ex * az, lz = ex * ay - ey * ax;
      const mx = ay * bz - az * by, my = az * bx - ax * bz, mz = ax * by - ay * bx;
      g[a] += f * kx; g[a + 1] += f * ky; g[a + 2] += f * kz;
      g[b] += f * lx; g[b + 1] += f * ly; g[b + 2] += f * lz;
      g[e] += f * mx; g[e + 1] += f * my; g[e + 2] += f * mz;
      g[c] -= f * (kx + lx + mx); g[c + 1] -= f * (ky + ly + my); g[c + 2] -= f * (kz + lz + mz);
    };
    for (let t = 0; t < T.cc.length; t++) {
      const s = T.cs[t], lo = T.cv[t];
      vol(T.cc[t] * D, T.ca[t] * D, T.cb[t] * D, T.cd[t] * D, V => { const sv = s * V; if (sv >= lo) return null; const v = sv - lo; return [W_CHIRAL * v * v, 2 * W_CHIRAL * v * s]; });
    }
    if (wPlanar > 0) {
      for (let t = 0; t < T.pc.length; t++) vol(T.pc[t] * D, T.pa[t] * D, T.pb[t] * D, T.pd[t] * D, V => [wPlanar * V * V, 2 * wPlanar * V]);
    }
    // linear centres: W_LIN * (1 + cos(theta)), stiff where distance bounds are not
    for (let t = 0; t < T.la.length; t++) {
      const a = T.la[t] * D, j = T.lj[t] * D, b = T.lb[t] * D;
      let uu = 0, vv = 0, uv = 0;
      for (let c = 0; c < D; c++) { const u = x[a + c] - x[j + c], v = x[b + c] - x[j + c]; uu += u * u; vv += v * v; uv += u * v; }
      const lu = Math.sqrt(uu) + 1e-12, lv = Math.sqrt(vv) + 1e-12, cs = uv / (lu * lv);
      E += W_LIN * (1 + cs);
      for (let c = 0; c < D; c++) {
        const u = x[a + c] - x[j + c], v = x[b + c] - x[j + c];
        const ga = W_LIN * (v / (lu * lv) - cs * u / uu), gb = W_LIN * (u / (lu * lv) - cs * v / vv);
        g[a + c] += ga; g[b + c] += gb; g[j + c] -= ga + gb;
      }
    }
    // 3-fold torsion preference on the first three coordinates: w * (1 + cos 3phi)
    if (wTor > 0) {
      for (let t = 0, nt = T.ti.length; t < nt; t++) {
        const i = T.ti[t] * D, j = T.tj[t] * D, k = T.tk[t] * D, l = T.tl[t] * D;
        const Fx = x[i] - x[j], Fy = x[i + 1] - x[j + 1], Fz = x[i + 2] - x[j + 2];
        const Gx = x[j] - x[k], Gy = x[j + 1] - x[k + 1], Gz = x[j + 2] - x[k + 2];
        const Hx = x[l] - x[k], Hy = x[l + 1] - x[k + 1], Hz = x[l + 2] - x[k + 2];
        const Ax = Fy * Gz - Fz * Gy, Ay = Fz * Gx - Fx * Gz, Az = Fx * Gy - Fy * Gx;
        const Bx = Hy * Gz - Hz * Gy, By = Hz * Gx - Hx * Gz, Bz = Hx * Gy - Hy * Gx;
        const AA = Ax * Ax + Ay * Ay + Az * Az, BB = Bx * Bx + By * By + Bz * Bz, GG = Gx * Gx + Gy * Gy + Gz * Gz;
        if (AA < 1e-10 || BB < 1e-10 || GG < 1e-10) continue;
        const gl = Math.sqrt(GG);
        const sy = ((By * Az - Bz * Ay) * Gx + (Bz * Ax - Bx * Az) * Gy + (Bx * Ay - By * Ax) * Gz) / gl;
        const phi = Math.atan2(sy, Ax * Bx + Ay * By + Az * Bz);
        const w = T.tw[t] * wTor;
        E += w * (1 + Math.cos(3 * phi));
        const dE = -3 * w * Math.sin(3 * phi);
        const fa = -gl / AA, fb = gl / BB;
        const fja = (Fx * Gx + Fy * Gy + Fz * Gz) / (AA * gl), fjb = (Hx * Gx + Hy * Gy + Hz * Gz) / (BB * gl);
        g[i] += dE * fa * Ax; g[i + 1] += dE * fa * Ay; g[i + 2] += dE * fa * Az;
        g[l] += dE * fb * Bx; g[l + 1] += dE * fb * By; g[l + 2] += dE * fb * Bz;
        g[j] += dE * ((fja - fa) * Ax - fjb * Bx); g[j + 1] += dE * ((fja - fa) * Ay - fjb * By); g[j + 2] += dE * ((fja - fa) * Az - fjb * Bz);
        g[k] += dE * ((fjb - fb) * Bx - fja * Ax); g[k + 1] += dE * ((fjb - fb) * By - fja * Ay); g[k + 2] += dE * ((fjb - fb) * Bz - fja * Az);
      }
    }
    if (D === 4 && w4 > 0) for (let i = 3; i < x.length; i += 4) { E += w4 * x[i] * x[i]; g[i] += 2 * w4 * x[i]; }
    return E;
  }

  /** Limited-memory BFGS with a backtracking line search and a per-coordinate step cap. Minimises in place. */
  function lbfgs(x, fg, maxIter, gtol, maxStep) {
    const n = x.length, M = 8;
    const g = new Float64Array(n), gN = new Float64Array(n), xN = new Float64Array(n), d = new Float64Array(n);
    const S = [], Y = [], R = [], alpha = new Float64Array(M), pool = [];
    let f = fg(x, g), it = 0;
    for (; it < maxIter; it++) {
      let gmax = 0;
      for (let i = 0; i < n; i++) { const a = g[i] < 0 ? -g[i] : g[i]; if (a > gmax) gmax = a; }
      if (gmax < gtol) break;
      for (let i = 0; i < n; i++) d[i] = -g[i];
      const m = S.length;
      for (let k = m - 1; k >= 0; k--) {
        const Sk = S[k], Yk = Y[k];
        let s = 0;
        for (let i = 0; i < n; i++) s += Sk[i] * d[i];
        alpha[k] = R[k] * s;
        for (let i = 0; i < n; i++) d[i] -= alpha[k] * Yk[i];
      }
      if (m) {
        const Sk = S[m - 1], Yk = Y[m - 1];
        let sy = 0, yy = 0;
        for (let i = 0; i < n; i++) { sy += Sk[i] * Yk[i]; yy += Yk[i] * Yk[i]; }
        const gam = sy / yy;
        for (let i = 0; i < n; i++) d[i] *= gam;
      }
      for (let k = 0; k < m; k++) {
        const Sk = S[k], Yk = Y[k];
        let s = 0;
        for (let i = 0; i < n; i++) s += Yk[i] * d[i];
        const b = R[k] * s;
        for (let i = 0; i < n; i++) d[i] += Sk[i] * (alpha[k] - b);
      }
      let gd = 0;
      for (let i = 0; i < n; i++) gd += g[i] * d[i];
      if (!(gd < 0)) { for (let i = 0; i < n; i++) d[i] = -g[i]; gd = 0; for (let i = 0; i < n; i++) gd -= g[i] * g[i]; S.length = Y.length = R.length = 0; }
      let dmax = 0;
      for (let i = 0; i < n; i++) { const a = d[i] < 0 ? -d[i] : d[i]; if (a > dmax) dmax = a; }
      let t = dmax > maxStep ? maxStep / dmax : 1, fN = Infinity, ok = false;
      for (let ls = 0; ls < 25; ls++) {
        for (let i = 0; i < n; i++) xN[i] = x[i] + t * d[i];
        fN = fg(xN, gN);
        if (fN <= f + 1e-4 * t * gd) { ok = true; break; }
        t *= 0.35;
      }
      if (!ok) { if (!S.length) break; S.length = Y.length = R.length = 0; continue; }
      let sNew, yNew;
      if (S.length === M) { sNew = S.shift(); yNew = Y.shift(); R.shift(); } else if (pool.length) [sNew, yNew] = pool.pop(); else { sNew = new Float64Array(n); yNew = new Float64Array(n); }
      let sy = 0;
      for (let i = 0; i < n; i++) { sNew[i] = xN[i] - x[i]; yNew[i] = gN[i] - g[i]; sy += sNew[i] * yNew[i]; }
      if (sy > 1e-12) { S.push(sNew); Y.push(yNew); R.push(1 / sy); } else pool.push([sNew, yNew]);
      const df = f - fN;
      x.set(xN); g.set(gN); f = fN;
      if (f < 1e-12 || (df >= 0 && df < 1e-11 * Math.max(1, f))) break;
    }
    return { f, iter: it };
  }

  function assessTerms(T, x) {
    let chiralWrong = 0, maxBond = 0;
    for (let t = 0; t < T.cc.length; t++) {
      const c = 3 * T.cc[t], a = 3 * T.ca[t], b = 3 * T.cb[t], e = 3 * T.cd[t];
      const u = [x[a] - x[c], x[a + 1] - x[c + 1], x[a + 2] - x[c + 2]], v = [x[b] - x[c], x[b + 1] - x[c + 1], x[b + 2] - x[c + 2]];
      const w = [x[e] - x[c], x[e + 1] - x[c + 1], x[e + 2] - x[c + 2]];
      if (T.cs[t] * vdot(u, vcross(v, w)) <= 0) chiralWrong++;
    }
    for (let t = 0; t < T.bi.length; t++) {
      if (T.kind[t] !== 1) continue;
      const i = 3 * T.bi[t], j = 3 * T.bj[t];
      const d = Math.hypot(x[i] - x[j], x[i + 1] - x[j + 1], x[i + 2] - x[j + 2]);
      const dev = Math.max(T.blo[t] - d, d - T.bhi[t], 0);
      if (dev > maxBond) maxBond = dev;
    }
    return { chiralWrong, maxBond };
  }

  function embedOnce(T, rng) {
    const m = T.m;
    if (m === 1) return { x3: new Float64Array(3), E: 0 };
    const L = 1 + 1.5 * Math.cbrt(m);
    const x4 = new Float64Array(4 * m);
    for (let i = 0; i < 4 * m; i++) x4[i] = (2 * rng() - 1) * L;
    lbfgs(x4, (x, g) => evalTerms(T, x, g, 4, 0, 0, 0), 400, 1e-3, 0.6);
    lbfgs(x4, (x, g) => evalTerms(T, x, g, 4, 0, 0.3, 1), 150, 1e-3, 0.4);
    lbfgs(x4, (x, g) => evalTerms(T, x, g, 4, 0, 3, 1), 150, 1e-3, 0.4);
    const x3 = new Float64Array(3 * m);
    for (let i = 0; i < m; i++) { x3[3 * i] = x4[4 * i]; x3[3 * i + 1] = x4[4 * i + 1]; x3[3 * i + 2] = x4[4 * i + 2]; }
    const r = lbfgs(x3, (x, g) => evalTerms(T, x, g, 3, W_PLANAR, 0, 1), 400, 1e-4, 0.3);
    return { x3, E: r.f, Ebounds: evalTerms(T, x3, new Float64Array(3 * m), 3, W_PLANAR, 0, 0) };
  }

  function placeTerminalH(model, comp, xyz) {
    const placed = new Uint8Array(model.n);
    for (const i of comp) if (!model.termH[i]) placed[i] = 1;
    const nbrsOf = j => model.adjO[j], isPlaced = i => placed[i] === 1;
    for (const p of comp) {
      if (model.termH[p]) continue;
      const hs = model.nbrs[p].filter(h => model.termH[h]);
      if (!hs.length) continue;
      const geom = geometryFor(model.info[p], model.Z[p], model.nbrs[p].length);
      const pos = placeAround(model.Z, xyz, nbrsOf, isPlaced, p, hs.length, geom, idealBondLength(model.Z[p], 1, 1));
      hs.forEach((h, k) => { xyz[3 * h] = pos[k][0]; xyz[3 * h + 1] = pos[k][1]; xyz[3 * h + 2] = pos[k][2]; placed[h] = 1; });
    }
  }

  /** Accept a force-field result only if it keeps chirality, double-bond geometry and sane bond lengths. */
  function polishAcceptable(model, T, before, after) {
    for (let i = 0; i < after.length; i++) if (!Number.isFinite(after[i])) return false;
    const a = assessTerms(T, after);
    if (a.chiralWrong > assessTerms(T, before).chiralWrong) return false;
    for (let t = 0; t < T.bi.length; t++) {
      if (T.kind[t] !== 1) continue;
      const i = 3 * T.bi[t], j = 3 * T.bj[t];
      const d = Math.hypot(after[i] - after[j], after[i + 1] - after[j + 1], after[i + 2] - after[j + 2]);
      if (Math.abs(d - (T.blo[t] + T.bhi[t]) / 2) > 0.25) return false;
    }
    const loc = new Map(T.atoms.map((v, k) => [v, k]));
    for (const pb of model.planarBonds) { // cis/trans relation of the reference pair must survive
      const r = pb.ref;
      const ids = [r.x, pb.j, pb.k, r.y].map(v => loc.get(v));
      if (ids.some(v => v === undefined)) continue;
      const dh = x => Math.abs(BM.geom ? BM.geom.dihedral(x, ids[0], ids[1], ids[2], ids[3]) : 0);
      if ((dh(before) < 90) !== (dh(after) < 90)) return false;
    }
    return true;
  }

  function embedComponent(model, comp, xyz, rng, warnings, opts) {
    const heavy = comp.filter(i => !model.termH[i]);
    const T = buildTerms(model, heavy);
    const m0 = heavy.length;
    const tries = opts.maxAttempts || (m0 <= 3 ? 2 : m0 <= 30 ? 8 : m0 <= 80 ? 4 : 2);
    let best = null;
    for (let a = 0; a < tries; a++) {
      const r = embedOnce(T, rng);
      const q = assessTerms(T, r.x3);
      const score = r.E + 50 * q.chiralWrong + 20 * Math.max(0, q.maxBond - 0.03);
      if (!best || score < best.score) best = { x3: r.x3, score, q };
      if (q.chiralWrong === 0 && q.maxBond < 0.03 && r.Ebounds < 0.02 + 0.004 * heavy.length && r.E - r.Ebounds < 0.05 * T.ti.length + 0.05) break;
    }
    heavy.forEach((i, k) => { xyz[3 * i] = best.x3[3 * k]; xyz[3 * i + 1] = best.x3[3 * k + 1]; xyz[3 * i + 2] = best.x3[3 * k + 2]; });
    placeTerminalH(model, comp, xyz);
    const TA = buildTerms(model, comp);
    const x = new Float64Array(3 * comp.length);
    comp.forEach((i, k) => { x[3 * k] = xyz[3 * i]; x[3 * k + 1] = xyz[3 * i + 1]; x[3 * k + 2] = xyz[3 * i + 2]; });
    lbfgs(x, (xx, g) => evalTerms(TA, xx, g, 3, W_PLANAR, 0, 1), 500, 1e-5, 0.2);
    if (opts.ff !== false && BM.ff && typeof BM.ff.optimize === 'function' && comp.length > 1) {
      try {
        const loc = new Map(comp.map((v, k) => [v, k])), sb = [];
        model.bonds.forEach((b, k) => { if (loc.has(b[0])) sb.push([loc.get(b[0]), loc.get(b[1]), model.aromB[k] ? 1.5 : b[2]]); });
        const sub = { Z: comp.map(i => model.Z[i]), xyz: Float64Array.from(x), bonds: sb, formalCharges: comp.map(i => model.charges[i] || 0) };
        const r = BM.ff.optimize(sub, { maxSteps: opts.ffSteps || 400 });
        const nx = r && r.xyz;
        if (nx && nx.length === x.length && polishAcceptable(model, TA, x, nx)) x.set(nx);
        else warnings.push('Force-field polish was rejected; kept the distance-geometry structure.');
      } catch (e) { warnings.push('Force-field polish failed: ' + (e && e.message)); }
    }
    const fin = assessTerms(TA, x);
    if (fin.chiralWrong) warnings.push('Could not satisfy ' + fin.chiralWrong + ' tetrahedral stereo centre(s).');
    if (fin.maxBond > 0.08) warnings.push('Strained structure: a bond deviates ' + fin.maxBond.toFixed(2) + ' Å from its ideal length.');
    let cx = 0, cy = 0, cz = 0;
    for (let k = 0; k < comp.length; k++) { cx += x[3 * k]; cy += x[3 * k + 1]; cz += x[3 * k + 2]; }
    cx /= comp.length; cy /= comp.length; cz /= comp.length;
    comp.forEach((i, k) => { xyz[3 * i] = x[3 * k] - cx; xyz[3 * i + 1] = x[3 * k + 1] - cy; xyz[3 * i + 2] = x[3 * k + 2] - cz; });
  }

  /** Lay out components along x with at least `gap` Å between them. */
  function arrangeComponents(comps, xyz, gap) {
    let maxX = null;
    for (const comp of comps) {
      let lo = Infinity, hi = -Infinity;
      for (const i of comp) { lo = Math.min(lo, xyz[3 * i]); hi = Math.max(hi, xyz[3 * i]); }
      if (maxX === null) { maxX = hi; continue; }
      const shift = maxX + gap - lo;
      for (const i of comp) xyz[3 * i] += shift;
      maxX = hi + shift;
    }
  }

  /** SMILES -> {mol (explicit H, 3D coordinates in Å), warnings, aromaticBonds}. opts: {seed, ff, maxAttempts}. */
  function smilesTo3D(smiles, opts) {
    opts = opts || {};
    const parsed = parseSmiles(smiles);
    const Z = parsed.atoms.map(a => a.Z), charges = parsed.atoms.map(a => a.charge || 0);
    const bonds = parsed.bonds.map(b => b.slice());
    const nHeavy = Z.length;
    parsed.atoms.forEach((a, i) => {
      for (let h = 0; h < a.hCount; h++) { const k = Z.length; Z.push(1); charges.push(0); bonds.push([i, k, 1]); }
    });
    const n = Z.length;
    const aromB = new Uint8Array(bonds.length);
    for (const k of parsed.aromaticBonds) aromB[k] = 1;
    markBenzenoid(n, bonds, aromB);
    const aromA = new Uint8Array(n);
    for (let i = 0; i < nHeavy; i++) if (parsed.atoms[i].aromatic) aromA[i] = 1;
    bonds.forEach((b, k) => { if (aromB[k]) { aromA[b[0]] = 1; aromA[b[1]] = 1; } });
    const model = buildModel(Z, bonds, charges, aromB, aromA, parsed.stereo);
    const warnings = model.warnings.slice();
    const rng = makeRng(opts.seed === undefined ? 1 : opts.seed);
    const xyz = new Float64Array(3 * n);
    const comps = fragmentsOf(n, bonds);
    for (const comp of comps) embedComponent(model, comp, xyz, rng, warnings, opts);
    arrangeComponents(comps, xyz, 3.2);
    const aromaticBonds = [];
    aromB.forEach((v, k) => { if (v) aromaticBonds.push(k); });
    const mol = { Z, xyz, bonds, formalCharges: charges };
    if (parsed.atoms.some(a => a.isotope)) mol.isotopes = Z.map((z, i) => (i < nHeavy ? parsed.atoms[i].isotope || 0 : 0));
    return { mol, warnings, aromaticBonds };
  }

  // ---- end of implementation (more sections are inserted above this line)

  BM.chem = {
    parseSmiles, toSmiles, smilesTo3D, perceiveBondOrders,
    valenceInfo, addHydrogens, removeHydrogens, adjustHydrogens,
    idealBondLength, newAtomPosition, addAtom, deleteAtoms, setBond, fragments,
    isMetal,
    _internal: { valenceList, groupOf, maxMatching, ringBondFlags, makeRng, kekulizeFractional, lbfgs }
  };
})(typeof window !== 'undefined' ? (window.BM = window.BM || {}) : (globalThis.BM = globalThis.BM || {}));
