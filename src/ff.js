/* BlueMolecule — clean-up force field (BM.ff).
 *
 * An original, generic molecular-mechanics model for turning sketched or distorted structures into reasonable
 * starting geometries before a quantum-chemistry job.  It is NOT a published force field (not UFF, not MMFF);
 * the parameters are generic and chosen for robustness rather than accuracy.
 *
 *   E = Σ bonds    ½ kb (r − r0)²             r0: the shared BlueMolecule ideal-length rule (idealBondLength)
 *     + Σ angles   ½ kθ (θ − θ0)²             θ0 from VSEPR-style hybridization, or from a coordination polyhedron
 *     + Σ torsions ½ V (1 − d cos nφ)         3-fold sp3–sp3, 2-fold planar for double/aromatic/conjugated bonds,
 *                                              weak 6-fold sp2–sp3
 *     + Σ sp2      K h²                       h: height of an sp2 centre above the plane of its three neighbours
 *     + Σ pairs    w ε (1 − r²/d0²)²  (r<d0)  soft repulsion, d0 = 0.9 (Rvdw,i + Rvdw,j); 1-4 pairs w = ½
 *
 * Units: kcal/mol, Å; angles are radians internally.  Dummy atoms (Z = 0) take part in nothing and never move.
 */
(function (BM) {
  'use strict';
  const el = BM.el;
  const RAD = Math.PI / 180;

  // ------------------------------------------------------------------------------------------- parameters
  const K_BOND = { 1: 700, 1.5: 900, 2: 1100, 3: 1400 }, K_BOND_METAL = 300;
  const K_ANGLE = 100, K_ANGLE_LINEAR = 150, K_ANGLE_DONOR = 50, K_ANGLE_POLY = 70;
  const V_DOUBLE = 40, V_AROMATIC = 25, V_CONJ = 5, V_AMIDE = 12, V_ESTER = 6, V_ARYL_ETHER = 3, V_SP2_SP3 = 0.5;
  const V_SP3 = 2.0, K_OOP = 40, EPS_VDW = 50, VDW_SCALE = 0.9, HBOND_D0 = 1.75, SCALE_14 = 0.5, SKIN = 1.0;

  // ------------------------------------------------------------------------------------------- elements
  const MAIN_VALENCE = {
    1: 1, 2: 2, 3: 1, 4: 2, 5: 3, 6: 4, 7: 5, 8: 6, 9: 7, 10: 8, 11: 1, 12: 2, 13: 3, 14: 4, 15: 5, 16: 6, 17: 7, 18: 8,
    19: 1, 20: 2, 31: 3, 32: 4, 33: 5, 34: 6, 35: 7, 36: 8, 37: 1, 38: 2, 49: 3, 50: 4, 51: 5, 52: 6, 53: 7, 54: 8,
    55: 1, 56: 2, 81: 3, 82: 4, 83: 5, 84: 6, 85: 7, 86: 8, 87: 1, 88: 2
  };
  const P_BLOCK_METALS = new Set([13, 31, 49, 50, 81, 82, 83, 84]);
  const SQUARE_PLANAR_METALS = new Set([78, 46, 79, 45, 77]); // Pt, Pd, Au, Rh, Ir
  const CHALCOGEN = new Set([8, 16, 34, 52]);
  const XH_LENGTH = { 1: 0.74, 5: 1.19, 6: 1.09, 7: 1.01, 8: 0.96, 14: 1.48, 15: 1.42, 16: 1.34 };

  function period(Z) { return Z <= 2 ? 1 : Z <= 10 ? 2 : Z <= 18 ? 3 : Z <= 36 ? 4 : Z <= 54 ? 5 : Z <= 86 ? 6 : 7; }
  /** s-, d- and f-block metals: their angles come from a coordination polyhedron chosen by coordination number. */
  function usesPolyhedron(Z) {
    if (Z === 3 || Z === 4 || Z === 11 || Z === 12 || Z === 19 || Z === 20 || Z === 37 || Z === 38 ||
      Z === 55 || Z === 56 || Z === 87 || Z === 88) return true;
    return (Z >= 21 && Z <= 30) || (Z >= 39 && Z <= 48) || (Z >= 57 && Z <= 80) || (Z >= 89 && Z <= 112);
  }
  function isMetal(Z) { return Z > 0 && (usesPolyhedron(Z) || P_BLOCK_METALS.has(Z)); }
  function valence(Z) { return MAIN_VALENCE[Z] || 0; }

  function orderFactor(o) { return o >= 2.5 ? 0.79 : o >= 1.75 ? 0.875 : o >= 1.25 ? 0.915 : 1; }
  /** The shared BlueMolecule ideal bond length (Å) for elements Zi, Zj and bond order 1 | 1.5 | 2 | 3. */
  function idealBondLength(Zi, Zj, order) {
    if (Zi === 1 || Zj === 1) {
      const X = Zi === 1 ? Zj : Zi;
      return X in XH_LENGTH ? XH_LENGTH[X] : el.covalent[X] + 0.31;
    }
    return (el.covalent[Zi] + el.covalent[Zj]) * orderFactor(order || 1);
  }
  function bondK(o) { return o >= 2.5 ? K_BOND[3] : o >= 1.75 ? K_BOND[2] : o >= 1.25 ? K_BOND[1.5] : K_BOND[1]; }

  // ------------------------------------------------------------------------------------------- vectors
  const dot3 = (a, b) => a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
  function unit3(v) { const l = Math.hypot(v[0], v[1], v[2]); return l > 1e-12 ? [v[0] / l, v[1] / l, v[2] / l] : [0, 0, 0]; }

  // ------------------------------------------------------------------------------------------- geometries
  const C120 = -0.5, S120 = Math.sqrt(3) / 2;
  const pent = k => [Math.cos(2 * Math.PI * k / 5), Math.sin(2 * Math.PI * k / 5), 0];
  /** Ideal coordination polyhedra: unit vectors from the centre to each site. */
  const GEOMETRY_LIST = [
    ['linear', 'Linear', [[0, 0, 1], [0, 0, -1]]],
    ['trigonal-planar', 'Trigonal planar', [[1, 0, 0], [C120, S120, 0], [C120, -S120, 0]]],
    ['tetrahedral', 'Tetrahedral', [[1, 1, 1], [-1, -1, 1], [-1, 1, -1], [1, -1, -1]].map(unit3)],
    ['square-planar', 'Square planar', [[1, 0, 0], [0, 1, 0], [-1, 0, 0], [0, -1, 0]]],
    ['seesaw', 'Seesaw', [[0, 0, 1], [0, 0, -1], [1, 0, 0], [C120, S120, 0]]],
    ['t-shaped', 'T-shaped', [[0, 0, 1], [0, 0, -1], [1, 0, 0]]],
    ['trigonal-bipyramidal', 'Trigonal bipyramidal', [[0, 0, 1], [0, 0, -1], [1, 0, 0], [C120, S120, 0], [C120, -S120, 0]]],
    ['square-pyramidal', 'Square pyramidal', [[0, 0, 1], [1, 0, 0], [0, 1, 0], [-1, 0, 0], [0, -1, 0]]],
    ['octahedral', 'Octahedral', [[1, 0, 0], [0, 1, 0], [-1, 0, 0], [0, -1, 0], [0, 0, 1], [0, 0, -1]]],
    ['pentagonal-bipyramidal', 'Pentagonal bipyramidal', [[0, 0, 1], [0, 0, -1], pent(0), pent(1), pent(2), pent(3), pent(4)]]
  ];
  const GEOMETRIES = {};
  for (const [id, name, vectors] of GEOMETRY_LIST) {
    GEOMETRIES[id] = { id, name, cn: vectors.length, vectors: vectors.map(v => v.map(c => (Math.abs(c) < 1e-15 ? 0 : c))) };
  }

  /** Best injective assignment of direction vectors to polyhedron sites, judged by all pairwise cosines
   * (rotation-invariant).  Returns {perm, score}: perm[a] = site index of direction a. */
  function assignSites(dirs, sites) {
    const n = dirs.length, m = sites.length;
    if (n > m) return null;
    const D = dirs.map(a => dirs.map(b => dot3(a, b)));
    const S = sites.map(a => sites.map(b => dot3(a, b)));
    let best = Infinity, bestPerm = null, budget = 400000;
    const perm = new Array(n).fill(-1), used = new Uint8Array(m);
    (function rec(a, acc) {
      if (acc >= best - 1e-12 || budget-- <= 0) return;
      if (a === n) { best = acc; bestPerm = perm.slice(); return; }
      for (let p = 0; p < m; p++) {
        if (used[p]) continue;
        let add = 0;
        for (let b = 0; b < a; b++) { const t = D[a][b] - S[p][perm[b]]; add += t * t; }
        used[p] = 1; perm[a] = p;
        rec(a + 1, acc + add);
        used[p] = 0;
      }
    })(0, 0);
    return { perm: bestPerm, score: best };
  }
  function fitScore(dirs, id) { const r = assignSites(dirs, GEOMETRIES[id].vectors); return r ? r.score : Infinity; }

  /** Coordination polyhedron for a metal with ligand directions `dirs` (null when positions are unknown). */
  function metalGeometry(Zm, cn, dirs, override) {
    if (override && GEOMETRIES[override] && GEOMETRIES[override].cn >= cn && cn >= 2) return override;
    switch (cn) {
      case 2: return 'linear';
      case 3: return 'trigonal-planar';
      case 4:
        if (SQUARE_PLANAR_METALS.has(Zm)) return 'square-planar';
        if (Zm === 28 && dirs && fitScore(dirs, 'square-planar') < fitScore(dirs, 'tetrahedral')) return 'square-planar';
        return 'tetrahedral';
      case 5:
        if (dirs && fitScore(dirs, 'square-pyramidal') < fitScore(dirs, 'trigonal-bipyramidal')) return 'square-pyramidal';
        return 'trigonal-bipyramidal';
      case 6: return 'octahedral';
      case 7: return 'pentagonal-bipyramidal';
      default: return null;
    }
  }
  /** Hypervalent main-group polyhedra from VSEPR steric number and number of bonded neighbours. */
  function mainGroupGeometry(sn, nb) {
    if (sn === 5) return nb === 5 ? 'trigonal-bipyramidal' : nb === 4 ? 'seesaw' : nb === 3 ? 't-shaped' : nb === 2 ? 'linear' : null;
    if (sn === 6) return nb === 6 ? 'octahedral' : nb === 5 ? 'square-pyramidal' : nb === 4 ? 'square-planar' : nb === 3 ? 't-shaped' : nb === 2 ? 'linear' : null;
    if (sn === 7) return nb === 7 ? 'pentagonal-bipyramidal' : nb === 6 ? 'octahedral' : null;
    return null;
  }

  // ------------------------------------------------------------------------------------------- topology
  /** Neighbour lists of real atoms (dummies and duplicate bonds dropped). */
  function buildAdjacency(mol) {
    const Z = mol.Z, n = Z.length;
    const nbr = Array.from({ length: n }, () => []), ord = Array.from({ length: n }, () => []);
    const seen = new Set(), bonds = [];
    for (const b of (mol.bonds || [])) {
      const i = b[0] | 0, j = b[1] | 0;
      let o = +b[2]; if (!(o > 0)) o = 1;
      if (i === j || i < 0 || j < 0 || i >= n || j >= n || !(Z[i] > 0) || !(Z[j] > 0)) continue;
      const key = i < j ? i * n + j : j * n + i;
      if (seen.has(key)) continue;
      seen.add(key);
      nbr[i].push(j); ord[i].push(o); nbr[j].push(i); ord[j].push(o);
      bonds.push(i < j ? [i, j, o] : [j, i, o]);
    }
    return { nbr, ord, bonds };
  }

  /** Per-atom perception: VSEPR steric number, hybridization, ideal angle or coordination polyhedron.
   * opts: {xyz, adj, geometries: {atomIndex: geometryId}}.  Returns {atoms, nbr, ord, bonds, metal}. */
  function perceive(mol, opts) {
    opts = opts || {};
    const Z = mol.Z, n = Z.length;
    const xyz = opts.xyz || mol.xyz;
    const adj = opts.adj || buildAdjacency(mol);
    const { nbr, ord } = adj;
    const fc = mol.formalCharges;
    const over = opts.geometries || {};
    const metal = new Uint8Array(n);
    for (let i = 0; i < n; i++) metal[i] = isMetal(Z[i]) ? 1 : 0;
    const haveXYZ = xyz && xyz.length >= 3 * n;
    const dirTo = (i, j) => unit3([xyz[3 * j] - xyz[3 * i], xyz[3 * j + 1] - xyz[3 * i + 1], xyz[3 * j + 2] - xyz[3 * i + 2]]);
    const atoms = new Array(n);
    for (let i = 0; i < n; i++) {
      const A = { nb: 0, sn: 0, lp: 0, hyb: 0, theta0: 0, geom: null, sites: null, list: [], metal: !!metal[i], arom: false, pi: false, planarN: false, conjO: false };
      atoms[i] = A;
      if (!(Z[i] > 0)) continue;
      if (usesPolyhedron(Z[i])) { A.list = nbr[i].slice(); A.nb = A.list.length; continue; }
      let S = 0, multi = false;
      for (let t = 0; t < nbr[i].length; t++) {
        const j = nbr[i][t], o = ord[i][t];
        if (metal[j] && !metal[i]) {
          // A metal bonded to another neighbour of i is π-bonded (η² or higher): leave it out of i's geometry.
          let pi = false;
          for (const k of nbr[i]) if (k !== j && nbr[j].indexOf(k) >= 0) { pi = true; break; }
          if (pi) continue;
          A.list.push(j); S += 2; continue;              // dative bond: the ligand supplies the pair
        }
        A.list.push(j);
        if (o > 1.25 && o < 1.75) A.arom = true;
        if (o > 1.25 && !metal[j]) multi = true;
        S += (o === 1 && CHALCOGEN.has(Z[j]) && nbr[j].length === 1) ? 2 : o; // terminal X–O read as X=O
      }
      const nb = A.list.length;
      const q = fc ? (+fc[i] || 0) : 0;
      let lp = Z[i] === 1 ? 0 : Math.max(0, Math.floor((valence(Z[i]) - q - S) / 2 + 1e-9));
      if (period(Z[i]) <= 2) lp = Math.min(lp, Math.max(0, 4 - nb));
      let sn = nb + lp;
      if (A.arom && nb <= 3) { sn = 3; lp = Math.max(0, 3 - nb); }
      A.nb = nb; A.lp = lp; A.sn = sn;
      if (Z[i] === 1 || nb === 0) continue;
      if (sn >= 2 && sn <= 4) {
        A.hyb = sn - 1;
        if (sn === 2) A.theta0 = 180;
        else if (sn === 3) A.theta0 = 120;
        else if (lp === 0) A.theta0 = 109.47;
        else if (period(Z[i]) <= 2) A.theta0 = lp === 1 ? 107 : 104.5;
        else A.theta0 = lp === 1 ? 100 : 98;
        A.pi = multi && sn <= 3;
      } else if (sn >= 5 && period(Z[i]) >= 3) {
        A.geom = (over[i] && GEOMETRIES[over[i]] && GEOMETRIES[over[i]].cn >= nb) ? over[i] : mainGroupGeometry(sn, nb);
      }
    }
    // Conjugation: amide/aniline-type N becomes planar; ester/phenol-type O is flagged for planar torsions.
    for (let i = 0; i < n; i++) {
      const A = atoms[i];
      if (A.metal || !(Z[i] === 7 || Z[i] === 8) || A.sn !== 4) continue;
      if (A.list.length !== nbr[i].length || A.list.some(j => metal[j])) continue;
      if (!ord[i].every(o => o === 1)) continue;
      if (!A.list.some(j => atoms[j].pi)) continue;
      if (Z[i] === 7 && A.lp === 1 && A.nb === 3) { A.sn = 3; A.lp = 0; A.hyb = 2; A.theta0 = 120; A.planarN = true; }
      else if (Z[i] === 8 && A.lp === 2 && A.nb === 2) { A.conjO = true; A.theta0 = 112; }
    }
    // Polyhedral centres: choose the polyhedron and assign each neighbour to a site.
    for (let i = 0; i < n; i++) {
      const A = atoms[i];
      if (!(Z[i] > 0) || A.list.length < 2) continue;
      const poly = usesPolyhedron(Z[i]);
      if (!poly && !A.geom) continue;
      let dirs = null;
      if (haveXYZ) {
        dirs = A.list.map(j => dirTo(i, j));
        if (dirs.some(d => d[0] === 0 && d[1] === 0 && d[2] === 0)) dirs = null;
      }
      if (poly) A.geom = metalGeometry(Z[i], A.list.length, dirs, over[i]);
      if (!A.geom) continue;
      const sites = GEOMETRIES[A.geom].vectors;
      if (sites.length < A.list.length) { A.geom = null; continue; }
      const perm = dirs ? assignSites(dirs, sites).perm : A.list.map((_, t) => t);
      A.sites = perm.map(p => sites[p]);
    }
    return { atoms, nbr, ord, bonds: adj.bonds, metal };
  }

  // ------------------------------------------------------------------------------------------- setup
  /** Build the force-field context (topology, parameters, neighbour list) for `mol`.
   * opts: {frozen: number[], geometries: {atomIndex: geometryId}} — geometries forces a coordination polyhedron. */
  function setup(mol, opts) {
    opts = opts || {};
    const n = mol.Z.length;
    const Z = Int32Array.from(mol.Z, z => (z > 0 ? z | 0 : 0));
    const x0 = Float64Array.from(mol.xyz);
    if (x0.length !== 3 * n) throw new Error('ff.setup: mol.xyz must hold 3 coordinates per atom.');
    const adj = buildAdjacency(mol);
    const P = perceive(mol, { xyz: x0, adj, geometries: opts.geometries });
    const { nbr, ord } = adj, atoms = P.atoms, metal = P.metal;
    const bonded = (a, b) => nbr[a].indexOf(b) >= 0;
    const raw = { bonds: [], angles: [], tors: [], oop: [] };

    for (const [i, j, o] of adj.bonds) {
      raw.bonds.push([i, j, (metal[i] || metal[j]) ? K_BOND_METAL : bondK(o), idealBondLength(Z[i], Z[j], o)]);
    }

    // Angles.  3-membered rings: no angle term (the bonds fix the triangle); 4-membered rings: θ0 = 90°.
    const small4 = (a, j, b) => {
      if (metal[a] || metal[j] || metal[b]) return false;
      for (const m of nbr[a]) if (m !== j && !metal[m] && nbr[b].indexOf(m) >= 0) return true;
      return false;
    };
    for (let j = 0; j < n; j++) {
      const A = atoms[j], L = A.list;
      if (Z[j] <= 1 || L.length < 2) continue;
      if (A.sites) {
        for (let p = 0; p < L.length; p++) for (let q = p + 1; q < L.length; q++) {
          if (bonded(L[p], L[q])) continue;
          const c = Math.max(-1, Math.min(1, dot3(A.sites[p], A.sites[q])));
          raw.angles.push([L[p], j, L[q], K_ANGLE_POLY, Math.acos(c)]);
        }
      } else if (A.theta0 > 0) {
        for (let p = 0; p < L.length; p++) for (let q = p + 1; q < L.length; q++) {
          const a = L[p], b = L[q];
          if ((metal[a] && metal[b]) || bonded(a, b)) continue;
          const th = (A.sn >= 3 && small4(a, j, b)) ? 90 : A.theta0;
          const k = (metal[a] || metal[b]) ? K_ANGLE_DONOR : th > 179 ? K_ANGLE_LINEAR : K_ANGLE;
          raw.angles.push([a, j, b, k, th * RAD]);
        }
      }
    }

    // Torsions about bonds between sp2/sp3 non-metal atoms; barrier V is shared by all dihedrals of the bond.
    const OS = new Set([8, 16]);
    const doubleTo = (c, els) => {
      for (let t = 0; t < nbr[c].length; t++) if (ord[c][t] >= 1.75 && els.has(Z[nbr[c][t]])) return true;
      return false;
    };
    for (const [j, k, o] of adj.bonds) {
      if (metal[j] || metal[k]) continue;
      const aj = atoms[j], ak = atoms[k];
      if (aj.hyb < 2 || ak.hyb < 2 || nbr[j].length < 2 || nbr[k].length < 2) continue;
      let ring3 = false;
      for (const m of nbr[j]) if (m !== k && nbr[k].indexOf(m) >= 0) { ring3 = true; break; }
      if (ring3) continue;
      const s2j = aj.hyb === 2, s2k = ak.hyb === 2;
      let V, nf, d;
      if (o > 1.25 && s2j && s2k) { V = o >= 1.75 ? V_DOUBLE : V_AROMATIC; nf = 2; d = 1; }
      else if (s2j && s2k) {
        const amide = (aj.planarN && doubleTo(k, OS)) || (ak.planarN && doubleTo(j, OS));
        V = amide ? V_AMIDE : V_CONJ; nf = 2; d = 1;
      } else if ((s2j && ak.conjO) || (s2k && aj.conjO)) {
        V = doubleTo(s2j ? j : k, OS) ? V_ESTER : V_ARYL_ETHER; nf = 2; d = 1;
      } else if (aj.hyb === 3 && ak.hyb === 3) {
        if (CHALCOGEN.has(Z[j]) && CHALCOGEN.has(Z[k])) { V = V_SP3; nf = 2; d = -1; }
        else {
          const cc = Z[j] === 6 && Z[k] === 6, cn = (Z[j] === 6 && Z[k] === 7) || (Z[j] === 7 && Z[k] === 6);
          V = V_SP3 * (cc ? 1 : cn ? 0.8 : 0.5); nf = 3; d = -1;
        }
      } else { V = V_SP2_SP3; nf = 6; d = 1; }
      const list = [];
      for (const i of nbr[j]) if (i !== k) for (const l of nbr[k]) if (l !== j && l !== i) list.push([i, l]);
      for (const [i, l] of list) raw.tors.push([i, j, k, l, V / list.length, nf, d]);
    }

    // Out-of-plane term for trigonal (sp2 / aromatic / amide N) centres.
    for (let c = 0; c < n; c++) {
      const A = atoms[c];
      if (A.metal || A.hyb !== 2 || A.list.length !== 3) continue;
      raw.oop.push([c, A.list[0], A.list[1], A.list[2], K_OOP]);
    }

    // Nonbonded exclusions (1-2, 1-3) and scaled 1-4 pairs.  Metals without angle terms (CN ≥ 8, η-ligands)
    // do not exclude their ligands from one another, so the ligands keep apart by repulsion alone.
    const ex = new Set(), ex14 = new Set();
    const key = (a, b) => (a < b ? a * n + b : b * n + a);
    const open = j => metal[j] && !atoms[j].sites;
    for (const [i, j] of adj.bonds) ex.add(key(i, j));
    for (let j = 0; j < n; j++) {
      if (open(j)) continue;
      const L = nbr[j];
      for (let p = 0; p < L.length; p++) for (let q = p + 1; q < L.length; q++) ex.add(key(L[p], L[q]));
    }
    for (const [j, k] of adj.bonds) {
      if (open(j) || open(k)) continue;
      for (const i of nbr[j]) {
        if (i === k) continue;
        for (const l of nbr[k]) {
          if (l === j || l === i) continue;
          const kk = key(i, l);
          if (!ex.has(kk)) ex14.add(kk);
        }
      }
    }
    const R = new Float64Array(n), hbH = new Uint8Array(n), hbA = new Uint8Array(n);
    let rmax = 1;
    for (let i = 0; i < n; i++) {
      if (!(Z[i] > 0)) continue;
      R[i] = el.vdw[Z[i]] || 2.0;
      if (R[i] > rmax) rmax = R[i];
      if (Z[i] === 7 || Z[i] === 8 || Z[i] === 9) hbA[i] = 1;
      if (Z[i] === 1 && nbr[i].some(j => Z[j] === 7 || Z[j] === 8 || Z[j] === 9)) hbH[i] = 1;
    }
    const ctx = {
      n, Z, x0, adj, perception: P, raw, ex, ex14, R, hbH, hbA, geometries: opts.geometries || null,
      cutList: 2 * VDW_SCALE * rmax + SKIN, nl: { ref: null, nAct: 0, nIn: 0, P: null, W: null },
      free: null, frozen: [], nFree: 0, B: null, A: null, T: null, O: null
    };
    setFrozenCtx(ctx, opts.frozen || []);
    return ctx;
  }

  /** Terms touching at least one free atom are packed first, so optimizers can skip constant terms. */
  function packRows(rows, nIdx, nPar, free) {
    const act = [], ina = [];
    for (const r of rows) {
      let a = false;
      for (let t = 0; t < nIdx; t++) if (free[r[t]]) { a = true; break; }
      (a ? act : ina).push(r);
    }
    const all = act.concat(ina), n = all.length;
    const idx = new Int32Array(n * nIdx), par = new Float64Array(n * nPar);
    for (let t = 0; t < n; t++) {
      const r = all[t];
      for (let c = 0; c < nIdx; c++) idx[t * nIdx + c] = r[c];
      for (let c = 0; c < nPar; c++) par[t * nPar + c] = r[nIdx + c];
    }
    return { n, nAct: act.length, idx, par };
  }
  function setFrozenCtx(ctx, list) {
    const n = ctx.n, free = new Uint8Array(n);
    for (let i = 0; i < n; i++) free[i] = ctx.Z[i] > 0 ? 1 : 0;
    const frozen = new Set();
    for (const f of (list || [])) { const i = +f; if (Number.isInteger(i) && i >= 0 && i < n) { free[i] = 0; frozen.add(i); } }
    let nFree = 0;
    for (let i = 0; i < n; i++) nFree += free[i];
    ctx.free = free; ctx.frozen = [...frozen].sort((a, b) => a - b); ctx.nFree = nFree;
    ctx.B = packRows(ctx.raw.bonds, 2, 2, free);
    ctx.A = packRows(ctx.raw.angles, 3, 2, free);
    ctx.T = packRows(ctx.raw.tors, 4, 3, free);
    ctx.O = packRows(ctx.raw.oop, 4, 1, free);
    ctx.nl.ref = null;
  }

  // ------------------------------------------------------------------------------------------- neighbour list
  function buildPairs(ctx, x) {
    const n = ctx.n, Z = ctx.Z, R = ctx.R, free = ctx.free, ex = ctx.ex, ex14 = ctx.ex14, hbH = ctx.hbH, hbA = ctx.hbA;
    const cut = ctx.cutList, cut2 = cut * cut;
    let mnx = Infinity, mny = Infinity, mnz = Infinity;
    for (let i = 0; i < n; i++) {
      if (!(Z[i] > 0)) continue;
      if (x[3 * i] < mnx) mnx = x[3 * i];
      if (x[3 * i + 1] < mny) mny = x[3 * i + 1];
      if (x[3 * i + 2] < mnz) mnz = x[3 * i + 2];
    }
    const cell = new Int32Array(3 * n), next = new Int32Array(n).fill(-1), head = new Map();
    const K = (a, b, c) => (a * 4096 + b) * 4096 + c;
    for (let i = 0; i < n; i++) {
      if (!(Z[i] > 0)) continue;
      const a = Math.floor((x[3 * i] - mnx) / cut), b = Math.floor((x[3 * i + 1] - mny) / cut), c = Math.floor((x[3 * i + 2] - mnz) / cut);
      cell[3 * i] = a; cell[3 * i + 1] = b; cell[3 * i + 2] = c;
      const k = K(a, b, c), h = head.get(k);
      next[i] = h === undefined ? -1 : h;
      head.set(k, i);
    }
    const act = [], ina = [];
    for (let i = 0; i < n; i++) {
      if (!(Z[i] > 0)) continue;
      const xi = x[3 * i], yi = x[3 * i + 1], zi = x[3 * i + 2];
      for (let da = -1; da <= 1; da++) for (let db = -1; db <= 1; db++) for (let dc = -1; dc <= 1; dc++) {
        let j = head.get(K(cell[3 * i] + da, cell[3 * i + 1] + db, cell[3 * i + 2] + dc));
        if (j === undefined) continue;
        for (; j >= 0; j = next[j]) {
          if (j <= i) continue;
          const dx = x[3 * j] - xi, dy = x[3 * j + 1] - yi, dz = x[3 * j + 2] - zi;
          if (dx * dx + dy * dy + dz * dz > cut2) continue;
          const kk = i * n + j;
          if (ex.has(kk)) continue;
          let d0 = VDW_SCALE * (R[i] + R[j]);
          if ((hbH[i] && hbA[j]) || (hbH[j] && hbA[i])) d0 = Math.min(d0, HBOND_D0);
          (free[i] || free[j] ? act : ina).push(i, j, ex14.has(kk) ? SCALE_14 : 1, d0 * d0);
        }
      }
    }
    const nAct = act.length / 4, nIn = ina.length / 4, all = act.concat(ina), m = nAct + nIn;
    const P = new Int32Array(2 * m), W = new Float64Array(2 * m);
    for (let t = 0; t < m; t++) { P[2 * t] = all[4 * t]; P[2 * t + 1] = all[4 * t + 1]; W[2 * t] = all[4 * t + 2]; W[2 * t + 1] = all[4 * t + 3]; }
    ctx.nl = { ref: Float64Array.from(x), nAct, nIn, P, W };
  }
  function needRebuild(ctx, x) {
    const ref = ctx.nl.ref;
    if (!ref || ref.length !== x.length) return true;
    const lim = 0.25 * SKIN * SKIN, Z = ctx.Z;
    for (let i = 0; i < ctx.n; i++) {
      if (!(Z[i] > 0)) continue;
      const dx = x[3 * i] - ref[3 * i], dy = x[3 * i + 1] - ref[3 * i + 1], dz = x[3 * i + 2] - ref[3 * i + 2];
      if (dx * dx + dy * dy + dz * dz > lim) return true;
    }
    return false;
  }

  // ------------------------------------------------------------------------------------------- energy + gradient
  /** Energy and analytic gradient (into g).  activeOnly skips terms whose atoms are all fixed (constant energy). */
  function evalCore(ctx, x, g, activeOnly) {
    g.fill(0);
    if (needRebuild(ctx, x)) buildPairs(ctx, x);
    const out = { energy: 0, bond: 0, angle: 0, torsion: 0, oop: 0, vdw: 0 };
    let e, m, I, Q;

    // bonds: ½ k (r − r0)²
    ({ idx: I, par: Q } = ctx.B); m = activeOnly ? ctx.B.nAct : ctx.B.n; e = 0;
    for (let t = 0; t < m; t++) {
      const i = 3 * I[2 * t], j = 3 * I[2 * t + 1], k = Q[2 * t], r0 = Q[2 * t + 1];
      let dx = x[j] - x[i], dy = x[j + 1] - x[i + 1], dz = x[j + 2] - x[i + 2];
      const r = Math.sqrt(dx * dx + dy * dy + dz * dz), dr = r - r0;
      e += 0.5 * k * dr * dr;
      let f;
      if (r > 1e-10) f = k * dr / r;
      else { dx = 1; dy = 0; dz = 0; f = k * dr; } // coincident atoms: push apart along x
      g[j] += f * dx; g[j + 1] += f * dy; g[j + 2] += f * dz;
      g[i] -= f * dx; g[i + 1] -= f * dy; g[i + 2] -= f * dz;
    }
    out.bond = e;

    // angles: ½ kθ (θ − θ0)²
    ({ idx: I, par: Q } = ctx.A); m = activeOnly ? ctx.A.nAct : ctx.A.n; e = 0;
    for (let t = 0; t < m; t++) {
      const i = 3 * I[3 * t], j = 3 * I[3 * t + 1], k = 3 * I[3 * t + 2], kt = Q[2 * t], t0 = Q[2 * t + 1];
      const ux = x[i] - x[j], uy = x[i + 1] - x[j + 1], uz = x[i + 2] - x[j + 2];
      const vx = x[k] - x[j], vy = x[k + 1] - x[j + 1], vz = x[k + 2] - x[j + 2];
      const ru2 = ux * ux + uy * uy + uz * uz, rv2 = vx * vx + vy * vy + vz * vz;
      if (ru2 < 1e-20 || rv2 < 1e-20) continue;
      const ru = Math.sqrt(ru2), rv = Math.sqrt(rv2);
      const cx = uy * vz - uz * vy, cy = uz * vx - ux * vz, cz = ux * vy - uy * vx;
      const sn = Math.sqrt(cx * cx + cy * cy + cz * cz) / (ru * rv), cs = (ux * vx + uy * vy + uz * vz) / (ru * rv);
      const dth = Math.atan2(sn, cs) - t0;
      e += 0.5 * kt * dth * dth;
      const de = kt * dth;
      let gix, giy, giz, gkx, gky, gkz;
      if (sn > 1e-6) {
        const f = -de / sn, a1 = f / (ru * rv), ai = f * cs / ru2, ak = f * cs / rv2;
        gix = a1 * vx - ai * ux; giy = a1 * vy - ai * uy; giz = a1 * vz - ai * uz;
        gkx = a1 * ux - ak * vx; gky = a1 * uy - ak * vy; gkz = a1 * uz - ak * vz;
      } else {
        // (anti)parallel arms: bend along a fixed perpendicular so the term can still act
        let px, py, pz;
        if (Math.abs(ux) < 0.9 * ru) { px = 0; py = uz; pz = -uy; } else { px = -uz; py = 0; pz = ux; }
        const pl = Math.sqrt(px * px + py * py + pz * pz); px /= pl; py /= pl; pz /= pl;
        const si = (cs < 0 ? -de : de) / ru, sk = -de / rv;
        gix = si * px; giy = si * py; giz = si * pz; gkx = sk * px; gky = sk * py; gkz = sk * pz;
      }
      g[i] += gix; g[i + 1] += giy; g[i + 2] += giz;
      g[k] += gkx; g[k + 1] += gky; g[k + 2] += gkz;
      g[j] -= gix + gkx; g[j + 1] -= giy + gky; g[j + 2] -= giz + gkz;
    }
    out.angle = e;

    // torsions: ½ V (1 − d cos nφ)   (gradient after Blondel & Karplus 1996)
    ({ idx: I, par: Q } = ctx.T); m = activeOnly ? ctx.T.nAct : ctx.T.n; e = 0;
    for (let t = 0; t < m; t++) {
      const i = 3 * I[4 * t], j = 3 * I[4 * t + 1], k = 3 * I[4 * t + 2], l = 3 * I[4 * t + 3];
      const V = Q[3 * t], nf = Q[3 * t + 1], d = Q[3 * t + 2];
      const Fx = x[i] - x[j], Fy = x[i + 1] - x[j + 1], Fz = x[i + 2] - x[j + 2];
      const Gx = x[j] - x[k], Gy = x[j + 1] - x[k + 1], Gz = x[j + 2] - x[k + 2];
      const Hx = x[l] - x[k], Hy = x[l + 1] - x[k + 1], Hz = x[l + 2] - x[k + 2];
      const Ax = Fy * Gz - Fz * Gy, Ay = Fz * Gx - Fx * Gz, Az = Fx * Gy - Fy * Gx;
      const Bx = Hy * Gz - Hz * Gy, By = Hz * Gx - Hx * Gz, Bz = Hx * Gy - Hy * Gx;
      const A2 = Ax * Ax + Ay * Ay + Az * Az, B2 = Bx * Bx + By * By + Bz * Bz, G2 = Gx * Gx + Gy * Gy + Gz * Gz;
      if (A2 < 1e-12 || B2 < 1e-12 || G2 < 1e-12) continue;
      const Gn = Math.sqrt(G2), rAB = Math.sqrt(A2 * B2);
      const cosp = (Ax * Bx + Ay * By + Az * Bz) / rAB;
      const sinp = ((By * Az - Bz * Ay) * Gx + (Bz * Ax - Bx * Az) * Gy + (Bx * Ay - By * Ax) * Gz) / (rAB * Gn);
      const phi = Math.atan2(sinp, cosp);
      e += 0.5 * V * (1 - d * Math.cos(nf * phi));
      const dE = 0.5 * V * d * nf * Math.sin(nf * phi);
      const FG = Fx * Gx + Fy * Gy + Fz * Gz, HG = Hx * Gx + Hy * Gy + Hz * Gz;
      const ci = -Gn / A2 * dE, cl = Gn / B2 * dE, cfa = FG / (A2 * Gn) * dE, chb = HG / (B2 * Gn) * dE;
      g[i] += ci * Ax; g[i + 1] += ci * Ay; g[i + 2] += ci * Az;
      g[l] += cl * Bx; g[l + 1] += cl * By; g[l + 2] += cl * Bz;
      g[j] += (cfa - ci) * Ax - chb * Bx; g[j + 1] += (cfa - ci) * Ay - chb * By; g[j + 2] += (cfa - ci) * Az - chb * Bz;
      g[k] += (chb - cl) * Bx - cfa * Ax; g[k + 1] += (chb - cl) * By - cfa * Ay; g[k + 2] += (chb - cl) * Bz - cfa * Az;
    }
    out.torsion = e;

    // out-of-plane: K h², h = height of centre c over the plane (a, b, d)
    ({ idx: I, par: Q } = ctx.O); m = activeOnly ? ctx.O.nAct : ctx.O.n; e = 0;
    for (let t = 0; t < m; t++) {
      const c = 3 * I[4 * t], a = 3 * I[4 * t + 1], b = 3 * I[4 * t + 2], d = 3 * I[4 * t + 3], K = Q[t];
      const e1x = x[b] - x[a], e1y = x[b + 1] - x[a + 1], e1z = x[b + 2] - x[a + 2];
      const e2x = x[d] - x[a], e2y = x[d + 1] - x[a + 1], e2z = x[d + 2] - x[a + 2];
      const nx = e1y * e2z - e1z * e2y, ny = e1z * e2x - e1x * e2z, nz = e1x * e2y - e1y * e2x;
      const nn = Math.sqrt(nx * nx + ny * ny + nz * nz);
      if (nn < 1e-10) continue;
      const hx = nx / nn, hy = ny / nn, hz = nz / nn;
      const dvx = x[c] - x[a], dvy = x[c + 1] - x[a + 1], dvz = x[c + 2] - x[a + 2];
      const h = dvx * hx + dvy * hy + dvz * hz;
      e += K * h * h;
      const dE = 2 * K * h;
      const gnx = (dvx - h * hx) / nn, gny = (dvy - h * hy) / nn, gnz = (dvz - h * hz) / nn;
      const bx = e2y * gnz - e2z * gny, by = e2z * gnx - e2x * gnz, bz = e2x * gny - e2y * gnx;
      const qx = gny * e1z - gnz * e1y, qy = gnz * e1x - gnx * e1z, qz = gnx * e1y - gny * e1x;
      g[c] += dE * hx; g[c + 1] += dE * hy; g[c + 2] += dE * hz;
      g[b] += dE * bx; g[b + 1] += dE * by; g[b + 2] += dE * bz;
      g[d] += dE * qx; g[d + 1] += dE * qy; g[d + 2] += dE * qz;
      g[a] -= dE * (hx + bx + qx); g[a + 1] -= dE * (hy + by + qy); g[a + 2] -= dE * (hz + bz + qz);
    }
    out.oop = e;

    // soft repulsion: w ε (1 − r²/d0²)² for r < d0
    const nl = ctx.nl, PP = nl.P, W = nl.W;
    m = activeOnly ? nl.nAct : nl.nAct + nl.nIn; e = 0;
    for (let t = 0; t < m; t++) {
      const d0sq = W[2 * t + 1], i = 3 * PP[2 * t], j = 3 * PP[2 * t + 1];
      const dx = x[j] - x[i], dy = x[j + 1] - x[i + 1], dz = x[j + 2] - x[i + 2];
      const r2 = dx * dx + dy * dy + dz * dz;
      if (r2 >= d0sq) continue;
      const w = W[2 * t] * EPS_VDW, q = 1 - r2 / d0sq;
      e += w * q * q;
      const f = -4 * w * q / d0sq;
      g[j] += f * dx; g[j + 1] += f * dy; g[j + 2] += f * dz;
      g[i] -= f * dx; g[i + 1] -= f * dy; g[i + 2] -= f * dz;
    }
    out.vdw = e;
    out.energy = out.bond + out.angle + out.torsion + out.oop + out.vdw;
    return out;
  }
  function maskGrad(ctx, g) {
    const f = ctx.free;
    for (let i = 0; i < ctx.n; i++) if (!f[i]) { g[3 * i] = 0; g[3 * i + 1] = 0; g[3 * i + 2] = 0; }
  }
  function rmsGrad(ctx, g) {
    if (!ctx.nFree) return 0;
    let s = 0;
    for (let k = 0; k < g.length; k++) s += g[k] * g[k];
    return Math.sqrt(s / (3 * ctx.nFree));
  }

  /** Energy (all terms) and gradient at xyz; gradient components of frozen and dummy atoms are zero. */
  function evaluate(ctx, xyz) {
    const x = xyz || ctx.x0;
    if (x.length !== 3 * ctx.n) throw new Error('ff.evaluate: xyz must hold 3 coordinates per atom.');
    const g = new Float64Array(3 * ctx.n);
    const r = evalCore(ctx, x, g, false);
    maskGrad(ctx, g);
    return { energy: r.energy, grad: g, terms: { bond: r.bond, angle: r.angle, torsion: r.torsion, oop: r.oop, vdw: r.vdw } };
  }

  // ------------------------------------------------------------------------------------------- L-BFGS
  const HISTORY = 8;
  /** Nudge free atoms that sit on top of another atom (< 0.01 Å): bond directions and angles are undefined
   * there, so no gradient method could separate them.  Deterministic directions; modifies x in place. */
  function separateCoincident(ctx, x) {
    const n = ctx.n, Z = ctx.Z, free = ctx.free, cell = 0.05, lim2 = 1e-4;
    const cellKey = (a, b, c) => (a & 1023) * 1048576 + (b & 1023) * 1024 + (c & 1023); // wraps; distances are checked
    const nudge = i => {
      const z = 1 - 2 * (((i + 1) * 0.6180339887) % 1), r = Math.sqrt(Math.max(0, 1 - z * z)), ph = (i + 1) * 2.399963229728653;
      x[3 * i] += 0.3 * r * Math.cos(ph); x[3 * i + 1] += 0.3 * r * Math.sin(ph); x[3 * i + 2] += 0.3 * z;
    };
    for (let pass = 0; pass < 4; pass++) {
      const grid = new Map();
      let moved = false;
      for (let i = 0; i < n; i++) {
        if (!(Z[i] > 0)) continue;
        const a = Math.floor(x[3 * i] / cell), b = Math.floor(x[3 * i + 1] / cell), c = Math.floor(x[3 * i + 2] / cell);
        let hit = -1;
        for (let da = -1; da <= 1 && hit < 0; da++) for (let db = -1; db <= 1 && hit < 0; db++) for (let dc = -1; dc <= 1 && hit < 0; dc++) {
          const list = grid.get(cellKey(a + da, b + db, c + dc));
          if (!list) continue;
          for (const j of list) {
            const dx = x[3 * i] - x[3 * j], dy = x[3 * i + 1] - x[3 * j + 1], dz = x[3 * i + 2] - x[3 * j + 2];
            if (dx * dx + dy * dy + dz * dz < lim2) { hit = j; break; }
          }
        }
        if (hit >= 0 && (free[i] || free[hit])) { nudge(free[i] ? i : hit); moved = true; }
        const k = cellKey(Math.floor(x[3 * i] / cell), Math.floor(x[3 * i + 1] / cell), Math.floor(x[3 * i + 2] / cell));
        const list = grid.get(k);
        if (list) list.push(i); else grid.set(k, [i]);
      }
      if (!moved) break;
    }
  }
  function makeState(ctx, x, opts) {
    const n3 = 3 * ctx.n;
    x = Float64Array.from(x);
    separateCoincident(ctx, x);
    const st = {
      ctx, x, g: new Float64Array(n3), E: 0, eConst: 0, rms: 0,
      d: new Float64Array(n3), xn: new Float64Array(n3), gn: new Float64Array(n3), alpha: new Float64Array(HISTORY),
      S: [], Y: [], RHO: [], pool: [],
      maxStep: opts.maxStep > 0 ? opts.maxStep : 0.3, tol: opts.tol > 0 ? opts.tol : 0.02,
      steps: 0, converged: false, stalled: false
    };
    refresh(st);
    return st;
  }
  function resetHistory(st) {
    while (st.S.length) st.pool.push(st.S.pop(), st.Y.pop());
    st.RHO.length = 0;
  }
  /** Re-evaluate after positions or the frozen set changed (drops the curvature history). */
  function refresh(st) {
    resetHistory(st);
    const ctx = st.ctx;
    const full = evalCore(ctx, st.x, st.gn, false).energy;
    st.E = evalCore(ctx, st.x, st.g, true).energy;
    maskGrad(ctx, st.g);
    st.eConst = full - st.E;               // terms among fixed atoms only: constant while they stay put
    st.rms = rmsGrad(ctx, st.g);
    st.converged = ctx.nFree === 0 || st.rms < st.tol;
    st.stalled = !Number.isFinite(st.E);
  }
  function direction(st) {
    const g = st.g, d = st.d, S = st.S, Y = st.Y, RHO = st.RHO, alpha = st.alpha, n3 = g.length, k = S.length;
    for (let t = 0; t < n3; t++) d[t] = -g[t];
    if (!k) { for (let t = 0; t < n3; t++) d[t] /= 300; return; }   // first step: scaled steepest descent
    for (let h = k - 1; h >= 0; h--) {
      const s = S[h], y = Y[h];
      let a = 0;
      for (let t = 0; t < n3; t++) a += s[t] * d[t];
      a *= RHO[h]; alpha[h] = a;
      for (let t = 0; t < n3; t++) d[t] -= a * y[t];
    }
    const sl = S[k - 1], yl = Y[k - 1];
    let sy = 0, yy = 0;
    for (let t = 0; t < n3; t++) { sy += sl[t] * yl[t]; yy += yl[t] * yl[t]; }
    const gamma = sy / yy;
    for (let t = 0; t < n3; t++) d[t] *= gamma;
    for (let h = 0; h < k; h++) {
      const s = S[h], y = Y[h];
      let b = 0;
      for (let t = 0; t < n3; t++) b += y[t] * d[t];
      const c = alpha[h] - RHO[h] * b;
      for (let t = 0; t < n3; t++) d[t] += c * s[t];
    }
  }
  /** One L-BFGS iteration with a per-atom step cap and a backtracking (Armijo) line search. */
  function iterate(st) {
    const ctx = st.ctx, x = st.x, g = st.g, d = st.d, xn = st.xn, gn = st.gn, n3 = x.length, free = ctx.free;
    if (ctx.nFree === 0) { st.converged = true; return false; }
    direction(st); maskGrad(ctx, d);
    let gd = 0;
    for (let t = 0; t < n3; t++) gd += g[t] * d[t];
    if (!(gd < 0)) {
      resetHistory(st); direction(st); maskGrad(ctx, d);
      gd = 0;
      for (let t = 0; t < n3; t++) gd += g[t] * d[t];
      if (!(gd < 0)) { st.stalled = true; return false; }
    }
    let mx = 0;
    for (let i = 0; i < ctx.n; i++) {
      if (!free[i]) continue;
      const s2 = d[3 * i] * d[3 * i] + d[3 * i + 1] * d[3 * i + 1] + d[3 * i + 2] * d[3 * i + 2];
      if (s2 > mx) mx = s2;
    }
    mx = Math.sqrt(mx);
    if (mx > st.maxStep) { const s = st.maxStep / mx; for (let t = 0; t < n3; t++) d[t] *= s; gd *= s; }
    let a = 1, En = NaN, ok = false;
    for (let ls = 0; ls < 12; ls++) {
      for (let t = 0; t < n3; t++) xn[t] = x[t] + a * d[t];
      En = evalCore(ctx, xn, gn, true).energy;
      if (En <= st.E + 1e-4 * a * gd) { ok = true; break; }
      const den = 2 * (En - st.E - a * gd);
      const an = Number.isFinite(den) && den > 0 ? -gd * a * a / den : 0.5 * a;
      a = Math.min(0.5 * a, Math.max(0.1 * a, an));
    }
    if (!ok) {
      if (st.S.length) { resetHistory(st); return true; }   // retry from steepest descent
      st.stalled = true; return false;
    }
    maskGrad(ctx, gn);
    const sv = st.pool.pop() || new Float64Array(n3), yv = st.pool.pop() || new Float64Array(n3);
    let sy = 0;
    for (let t = 0; t < n3; t++) { sv[t] = xn[t] - x[t]; yv[t] = gn[t] - g[t]; sy += sv[t] * yv[t]; }
    if (sy > 1e-10) {
      if (st.S.length === HISTORY) { st.pool.push(st.S.shift(), st.Y.shift()); st.RHO.shift(); }
      st.S.push(sv); st.Y.push(yv); st.RHO.push(1 / sy);
    } else st.pool.push(sv, yv);
    x.set(xn); g.set(gn); st.E = En; st.steps++;
    st.rms = rmsGrad(ctx, g);
    st.converged = st.rms < st.tol;
    return true;
  }

  /** Minimise mol's geometry.  opts: {frozen, maxSteps = 2000, tol = 0.02 (rms gradient, kcal/mol/Å),
   * maxStep = 0.3 (Å per atom per step), geometries}.  Returns {xyz, energy, steps, converged, rmsGrad}. */
  function optimize(mol, opts) {
    opts = opts || {};
    const ctx = setup(mol, opts);
    const st = makeState(ctx, ctx.x0, opts);
    const maxSteps = opts.maxSteps >= 0 ? Math.floor(opts.maxSteps) : 2000;
    let tries = 0;
    while (!st.converged && !st.stalled && st.steps < maxSteps && tries++ < 2 * maxSteps + 10) iterate(st);
    return { xyz: Float64Array.from(st.x), energy: st.E + st.eConst, steps: st.steps, converged: st.converged, rmsGrad: st.rms };
  }

  /** Incremental optimizer for live clean-up while atoms are dragged.  Topology changes only via rebuild(). */
  function createOptimizer(mol, opts) {
    opts = Object.assign({}, opts || {});
    let ctx = setup(mol, opts), st = makeState(ctx, ctx.x0, opts);
    const snapshot = () => ({ xyz: Float64Array.from(st.x), energy: st.E + st.eConst, rmsGrad: st.rms, converged: st.converged });
    return {
      step(nSteps) {
        const N = nSteps === undefined ? 1 : Math.max(0, Math.floor(nSteps) || 0);
        st.stalled = !Number.isFinite(st.E);
        let done = 0, tries = 0;
        while (done < N && !st.converged && !st.stalled && tries++ < 2 * N + 4) {
          const s0 = st.steps;
          iterate(st);
          if (st.steps > s0) done++;
        }
        return snapshot();
      },
      setFrozen(list) { opts.frozen = Array.from(list || []); setFrozenCtx(ctx, opts.frozen); refresh(st); },
      setPositions(xyz) {
        if (!xyz || xyz.length !== 3 * ctx.n) throw new Error('ff optimizer: setPositions expects ' + 3 * ctx.n + ' coordinates.');
        st.x = Float64Array.from(xyz);
        separateCoincident(ctx, st.x);
        refresh(st);
      },
      rebuild(newMol, o) {
        o = o || {};
        if (o.frozen) opts.frozen = Array.from(o.frozen);
        if (o.geometries !== undefined) opts.geometries = o.geometries;
        const n = newMol.Z.length;
        opts.frozen = (opts.frozen || []).filter(i => i >= 0 && i < n);
        ctx = setup(newMol, opts);
        st = makeState(ctx, ctx.x0, opts);
      },
      get xyz() { return Float64Array.from(st.x); },
      get context() { return ctx; }
    };
  }

  BM.ff = {
    setup, evaluate, optimize, createOptimizer,
    perceive, idealBondLength, isMetal, usesPolyhedron, valence, period,
    assignSites, metalGeometry, mainGroupGeometry, GEOMETRIES, GEOMETRY_IDS: GEOMETRY_LIST.map(g => g[0])
  };
})(typeof window !== 'undefined' ? (window.BM = window.BM || {}) : (globalThis.BM = globalThis.BM || {}));
