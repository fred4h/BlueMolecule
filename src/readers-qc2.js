/* BlueMolecule — more quantum-chemistry readers (2.3): Molpro and Dalton outputs, Turbomole program outputs and
 * GAMESS punch files.
 *
 *   molpro-out      Molpro output (.out) and the optimisation log next to it (.log)
 *   dalton-out      Dalton output (.out)
 *   turbomole-out   Turbomole program outputs: dscf / ridft (SCF), escf (excited states), aoforce (vibrations),
 *                   ricc2 (CC2, ADC(2) states) and the job.last / jobex logs of optimisations
 *   gamess-punch    GAMESS punch file (.dat): $DATA, the steps of an optimisation, energies and orbital energies
 *
 * Datasets follow parsers.js (Angstrom, Hartree, Debye, cm-1), with the helpers of readers-qc.js
 * (BM.qcReaders._h). Written from the layouts of public example outputs; no code from the programs themselves.
 */
(function (BM) {
  'use strict';
  const H = BM.qcReaders && BM.qcReaders._h;
  if (!H) return;
  const { BOHR, AU2DEBYE, HARTREE2EV, num, int, isNum, isInt, tok, lines, last, elZ, addJT, newThermo, Frames, setOpt,
    state, buildVib, finish, reg, isLinear, gamessOps, gamessExpand } = H;

  /** An XYZ block in Angstrom: the atom count, a title line, then "El x y z" rows. */
  function xyzBlock(L, j) {
    const n = int(tok(L[j] || '')[0]);
    if (!(n > 0) || tok(L[j]).length !== 1) return null;
    const Z = [], xyz = [];
    for (let k = j + 2; k < j + 2 + n && k < L.length; k++) {
      const r = tok(L[k]);
      if (r.length < 4 || !r.slice(1, 4).every(isNum)) return null;
      Z.push(elZ(r[0])); xyz.push(num(r[1]), num(r[2]), num(r[3]));
    }
    return Z.length === n && j + 2 + n < L.length ? { Z, xyz, title: (L[j + 1] || '').trim(), end: j + 2 + n } : null;
  }
  /** Sorted orbital energies with their occupations: {energies, nOcc}. */
  function orbitalSet(list) {
    if (!list || !list.length) return null;
    const s = list.slice().sort((a, b) => a.e - b.e);
    return { energies: s.map(o => o.e), nOcc: s.filter(o => o.occ > 0.5).length };
  }

  // =================================================================== Molpro
  const MOLPRO_SIG = /\*{3}\s+PROGRAM SYSTEM MOLPRO\s+\*{3}/;
  // the optimisation log next to the output has no banner: its program sections tell
  const MOLPRO_LOG = /^1?\s?PROGRAM \* (?:ALASKA|SEWARD|OPTG?|RHF-SCF|RKS-SCF|UHF-SCF|UKS-SCF)\b/m;
  const MOLPRO_SCF = /^(R|U|RO)?(HF|KS)(-SCF)?$/i;
  /** "ATOMIC COORDINATES": a header (NR ATOM CHARGE X Y Z) and rows in Bohr. */
  function molproAtoms(L, i) {
    let j = i + 1;
    while (j < L.length && j < i + 5 && !/^\s*NR\s+ATOM\s+CHARGE/.test(L[j])) j++;
    if (!/^\s*NR\s+ATOM\s+CHARGE/.test(L[j] || '')) return null;
    for (j++; j < L.length && !L[j].trim(); j++);
    const Z = [], xyz = [];
    for (; j < L.length; j++) {
      const r = tok(L[j]);
      if (r.length < 6 || !isInt(r[0]) || !isNum(r[2]) || !r.slice(3, 6).every(isNum)) break;
      let z = elZ(r[1]);
      if (!z && !/^(Q|X|DUM)/i.test(r[1])) z = Math.max(0, Math.round(num(r[2])));
      Z.push(z); xyz.push(num(r[3]) * BOHR, num(r[4]) * BOHR, num(r[5]) * BOHR);
    }
    // a table cut while the program was writing it (a job still running) makes no frame
    if (Z.length && (j >= L.length || (j === L.length - 1 && L[j] === ''))) return { Z: [], xyz: [], end: j };
    return { Z, xyz, end: j };
  }
  function parseMolpro(text, name) {
    const ds = BM.emptyDataset('molpro-out', name), meta = ds.meta, F = new Frames(ds), L = lines(text), n = L.length;
    ds.program = 'Molpro';
    const cmds = [];
    let inEcho = false, nucCharge = NaN, nA = NaN, nB = NaN, basis = null, libBasis = null, summary = null, normal = null, ks = false;
    let optHead = null, optRows = null, optDone = false, noConv = false, thr = {}, zpe = null;
    const freqs = [], modes = [], ir = [], sym = [];
    let orbs = null, orbA = null, orbB = null;
    const pending = {};
    let topKind = null;
    let barrier = -1;   // frames before the optimisation table: energies printed after it belong to a later geometry
    const energy = (kind, v, label) => {
      const f = F.cur();
      if (!f) { pending[kind] = v; return; }
      if (ds.frames.length <= barrier) return;   // (Molpro 2006 prints no final geometry: its energy is in the table)
      if (label && label !== '1.1' && f.energies && f.energies[kind] !== undefined) return;   // other states of a multi-state job
      f.energies = f.energies || {};
      f.energies[kind] = v;
      f.energy = v;   // each program prints its result after the one it builds on: the last is the job's level
      topKind = kind;
    };
    const addFrame = (Z, xyz, label) => {
      const f = F.add(Z, xyz, label);
      for (const k of Object.keys(pending)) { energy(k, pending[k]); delete pending[k]; }
      return f;
    };
    for (let i = 0; i < n; i++) {
      const s = L[i], t = s.trim();
      if (!t) continue;
      // the input, echoed after the header (Molpro 2006 and 2012) or after "Variables initialized" (2015 on), up to
      // "Commands initialized"
      if (/^default implementation of scratch files|^Variables initialized/.test(t)) { inEcho = true; continue; }
      if (inEcho) {
        if (/^Commands\s+initialized/.test(t)) inEcho = false;
        else if (t.replace(/!.*$/, '').trim()) cmds.push(t.replace(/!.*$/, '').trim());   // (comments start with "!")
        continue;
      }
      if (/^Optimization point \d+/.test(t)) { addJT(meta, 'opt'); continue; }
      let m;
      if ((m = /^Version\s+(\d{4}\.\d+)/.exec(t))) { ds.program = 'Molpro ' + m[1]; continue; }
      if (t === 'ATOMIC COORDINATES') {
        const g = molproAtoms(L, i);
        if (g) { if (g.Z.length) addFrame(g.Z, g.xyz); i = g.end - 1; }
        continue;
      }
      if (/^Current geometry \(xyz format, in Angstrom\)/.test(t)) {
        let j = i + 1;
        while (j < n && !L[j].trim()) j++;
        const g = xyzBlock(L, j);
        if (g) { addFrame(g.Z, g.xyz); i = g.end - 1; }
        continue;
      }
      if (t[0] === '!') {
        if ((m = /^!\s*([A-Z][\w()+\-[\]]*?)\s+(?:STATE\s+(\S+)\s+)?(?:total\s+)?energy\s+(-?\d+\.\d+)/i.exec(t))) {
          const meth = m[1].toUpperCase();
          if (/\[T\]|-T$/.test(meth)) continue;   // Molpro 2006 also prints CCSD[T] and CCSD-T before CCSD(T)
          if (MOLPRO_SCF.test(meth) && /KS/.test(meth)) ks = true;
          energy(MOLPRO_SCF.test(meth) ? 'SCF' : meth, num(m[3]), m[2]);
          continue;
        }
        if ((m = /^!\s*\S+\s+STATE\s+(\S+)\s+Dipole moment\s+(\S+)\s+(\S+)\s+(\S+)/i.exec(t)) && (m[1] === '1.1' || !ds.dipole)) {
          const x = num(m[2]) * AU2DEBYE, y = num(m[3]) * AU2DEBYE, z = num(m[4]) * AU2DEBYE;
          if ([x, y, z].every(Number.isFinite)) ds.dipole = { x, y, z, total: Math.hypot(x, y, z) };
        }
        continue;
      }
      if ((m = /^SETTING BASIS\s*=\s*(\S+)/.exec(t))) { basis = m[1]; continue; }
      if ((m = /^Library entry\s+\S+\s+\S+\s+(\S+)\s+selected/.exec(t))) { if (!libBasis) libBasis = m[1]; continue; }   // (the default basis)
      if ((m = /^NUCLEAR CHARGE:\s*(\d+)/.exec(t))) { nucCharge = int(m[1]); continue; }
      if ((m = /^NUMBER OF ELECTRONS:\s*(\d+)\+\s*(\d+)-/.exec(t))) { nA = int(m[1]); nB = int(m[2]); continue; }
      if (t === 'THRESHOLDS:') {
        for (let j = i + 1; j < n && j < i + 14; j++) {
          const re = /\b(OPTGRAD|OPTENERG|OPTSTEP)\s*=\s*(\S+)/g;
          let q;
          while ((q = re.exec(L[j]))) thr[q[1]] = num(q[2]);
        }
        continue;
      }
      if (/^ITER\.\s+ENERGY\(OLD\)\s+ENERGY\(NEW\)/.test(t)) {
        optHead = tok(t);
        optRows = [];
        barrier = ds.frames.length;
        let j = i + 1;
        for (; j < n; j++) {
          const r = tok(L[j]);
          if (!r.length) break;
          if (!isInt(r[0])) { if (/^[A-Za-z]/.test(r[0])) continue; break; }   // "Freezing grid" between the rows
          if (r.length < optHead.length - 1 || !r.slice(1, 10).every(isNum)) break;
          optRows.push(r.map(num));
        }
        i = j - 1;
        continue;
      }
      if (/^END OF GEOMETRY OPTIMIZATION/.test(t)) { optDone = true; continue; }
      if (/No convergence/i.test(t) && optHead) { noConv = true; continue; }
      if (/^Normal Modes/.test(t)) {
        // five modes per block: mode numbers (and symmetries), wavenumbers, intensities, then 3N displacement rows
        // "of low/zero frequencies": translations and rotations (0.00), and in Molpro 2006 the lowest vibrations too
        const low = /low\/zero/i.test(t), imag = /imaginary/i.test(t);
        const fq = low ? [] : freqs, md = low ? [] : modes, ix = low ? [] : ir, sy = low ? [] : sym;
        let cols = null, j = i + 1;
        for (; j < n; j++) {
          const u = L[j], v = u.trim();
          if (!v) { cols = null; if (j + 1 < n && !L[j + 1].trim()) break; continue; }
          let q;
          if ((q = /^Wavenumbers \[cm-1\]\s+(.*)$/.exec(v))) {
            const f = tok(q[1]).map(x => (imag ? -Math.abs(num(x)) : num(x)));
            cols = f.map((x, k) => fq.length + k);
            f.forEach(x => { fq.push(x); md.push([]); ix.push(NaN); sy.push(''); });
            const hd = tok(L[j - 1] || '');   // "1 Au 2 Bg ..." (symmetries) above
            if (hd.length === 2 * f.length) f.forEach((x, k) => { sy[cols[k]] = hd[2 * k + 1]; });
            continue;
          }
          if (!cols) { if (/^(\d+\s*\S*\s*)+$/.test(v)) continue; break; }
          if ((q = /^Intensities \[km\/mol\]\s+(.*)$/.exec(v))) { tok(q[1]).forEach((x, k) => { if (k < cols.length) ix[cols[k]] = num(x); }); continue; }
          if (/^Intensities/.test(v)) continue;
          if ((q = /^([A-Za-z]+?)([XYZ])(\d+)\s+(.*)$/.exec(v))) {
            const a = int(q[3]) - 1, ax = 'XYZ'.indexOf(q[2]);
            tok(q[4]).forEach((x, k) => { if (k < cols.length) md[cols[k]][3 * a + ax] = num(x); });
            continue;
          }
          break;
        }
        if (low) {
          const keep = fq.map((x, k) => k).filter(k => Math.abs(fq[k]) > 1);
          freqs.unshift(...keep.map(k => fq[k])); modes.unshift(...keep.map(k => md[k]));
          ir.unshift(...keep.map(k => ix[k])); sym.unshift(...keep.map(k => sy[k]));
        }
        i = j - 1;
        continue;
      }
      if ((m = /^Zero point energy:\s*(-?\d+\.\d+)\s*\[H\]/.exec(t))) { zpe = num(m[1]); continue; }
      if (/^Unique atom\s+s\s+p\b.*\bCharge$/.test(t)) {   // Mulliken populations (PROGRAM * POP): the charge is the last column
        const q = [];
        let j = i + 1;
        for (; j < n; j++) {
          const r = tok(L[j]);
          if (r.length < 4 || !isInt(r[0])) break;
          const c = /^[-+]$/.test(r[r.length - 2]) ? num(r[r.length - 2] + r[r.length - 1]) : num(r[r.length - 1]);
          if (!Number.isFinite(c)) break;
          q.push(c);
        }
        const fl = F.cur();
        if (fl && q.length === fl.Z.filter(z => z > 0).length) ds.charges = { Mulliken: q.length === fl.Z.length ? q : fl.Z.map(() => 0) };
        i = j - 1;
        continue;
      }
      if (/^ELECTRON ORBITALS/.test(t)) {
        orbs = [];
        if (/BETA|NEGATIVE/i.test(t)) orbB = orbs; else orbA = orbs;
        continue;
      }
      if (orbs) {
        // "1.1   2.00000   -10.01" (or "2", or "+" / "-" for one electron in the UHF tables of Molpro 2006 and 2012)
        if ((m = /^([1-9]\d*)\.([1-8])\s+([+-]|\d+(?:\.\d+)?)\s+(-?\d+\.\d+)(\s|$)/.exec(t))) { orbs.push({ occ: /^[+-]$/.test(m[3]) ? 1 : num(m[3]), e: num(m[4]) }); continue; }
        if (/^(HOMO|LUMO|\*{5}|PROGRAM|NATURAL\b|[A-Z][A-Z ]*ORBITALS\b|Orbital\s+Occupation)/.test(t)) { if (!/^Orbital/.test(t)) orbs = null; continue; }
        if (/^\d+\s+\d*[a-z]/.test(t) || /^=+$/.test(t)) continue;
      }
      if ((m = /^([\w()+\-*]+(?:-[\w()]+)*)\/(\S+)\s+energy=\s*(-?\d+\.\d+)/.exec(t))) { summary = { method: m[1], basis: m[2] }; continue; }
      if (/^Molpro calculation terminated|^Variable memory released/.test(t)) { normal = true; continue; }
      if (/^\? Error|^GLOBAL ERROR/.test(t)) normal = false;
    }
    // the input: title, basis, methods, job types, charge and spin
    const words = [];
    for (const c of cmds) for (const w of c.replace(/[{}]/g, ';').split(';')) { const x = w.trim(); if (x) words.push(x); }
    let method = null, fnl = null;
    for (const w of words) {
      const lw = w.toLowerCase();
      let q;
      if ((q = /^\*{3}\s*,\s*(.*)$/.exec(w))) { if (q[1].trim()) ds.title = q[1].trim(); continue; }
      if ((q = /^basis\s*[=,]\s*([^,;\s{]+)/i.exec(w)) && !basis) basis = q[1];
      if (/^optg\b/.test(lw)) addJT(meta, 'opt');
      if (/^(freq|frequencies)\b/.test(lw)) addJT(meta, 'freq');
      if ((q = /^set\s*,\s*charge\s*=\s*(-?\d+)/.exec(lw))) meta.charge = int(q[1]);
      const head = lw.split(/[,\s]/)[0];
      if (/^(r|u)?ks$/.test(head)) { const fn = w.split(',')[1]; fnl = fn && fn.trim() ? fn.trim().toUpperCase() : 'DFT'; }
      else if (/^(r|u|ro)?hf$/.test(head)) method = 'HF';
    }
    // the method: the level of the last energy printed ("!MP2 total energy"...), or the functional of a KS job
    if (topKind && topKind !== 'SCF') method = topKind;
    else if (fnl) method = fnl;
    else if (topKind === 'SCF' && !method) method = ks ? 'DFT' : 'HF';
    if (summary) basis = summary.basis;   // (the summary names the basis as written, but the last program run)
    if (!basis) basis = libBasis;
    meta.method = method;
    meta.basis = basis;
    if (Number.isFinite(nA) && Number.isFinite(nB)) {
      meta.multiplicity = Math.abs(nA - nB) + 1;
      if (Number.isFinite(nucCharge)) meta.charge = nucCharge - nA - nB;
    }
    meta.normalTermination = normal;
    if (!ds.frames.length) throw new Error('No ATOMIC COORDINATES found in this Molpro output.');
    // the log of a numerical frequency job: the reference geometry first, then one geometry per displacement
    if (/^\s*Cartesian Coordinates are in \[bohr\]/.test(L[0] || '') && meta.jobTypes.indexOf('opt') < 0 && ds.frames.length > 2) {
      ds.warnings.push('This is the log of a numerical frequency calculation: only its reference geometry is shown, not the ' +
        (ds.frames.length - 1) + ' displaced ones. The vibrations are in the output (.out).');
      ds.frames.length = 1;
      addJT(meta, 'freq');
    }
    // the optimisation: Molpro prints the steps as a table (energies and convergence criteria) and the geometries
    // only in its log, so the output has the first and the last geometry, and the table as cycles
    if (optRows && optRows.length) {
      addJT(meta, 'opt');
      const c = k => optHead.indexOf(k), g = (r, k) => (c(k) >= 0 ? r[c(k)] : NaN);
      const tg = Number.isFinite(thr.OPTGRAD) ? thr.OPTGRAD : 3e-4, te = Number.isFinite(thr.OPTENERG) ? thr.OPTENERG : 1e-6;
      const ts = Number.isFinite(thr.OPTSTEP) ? thr.OPTSTEP : 3e-4;
      const steps = optRows.map(r => {
        const o = {};
        setOpt(o, 'maxForce', g(r, 'GRADMAX'), tg);
        setOpt(o, 'gradNorm', g(r, 'GRADNORM'));
        setOpt(o, 'rmsForce', g(r, 'GRADRMS'));
        setOpt(o, 'maxDisp', g(r, 'STEPMAX'), ts);
        setOpt(o, 'rmsDisp', g(r, 'STEPRMS'));
        setOpt(o, 'stepNorm', g(r, 'STEPLEN'));
        setOpt(o, 'energyChange', g(r, 'DE'), te);
        return { step: r[0], energy: g(r, 'ENERGY(OLD)'), opt: o };
      });
      const lr = last(optRows);
      steps.push({ step: lr[0] + 1, energy: g(lr, 'ENERGY(NEW)'), opt: null, converged: optDone && !noConv });
      ds.optSteps = steps;
      const fl = last(ds.frames);
      if (ds.frames.length > 1) { fl.converged = optDone && !noConv; fl.opt = Object.assign({}, last(steps.slice(0, -1)).opt); }
      if (!/\.log$/i.test(name || '') && ds.frames.length === 1) {
        ds.warnings.push('This Molpro version prints the optimised geometry only as Z-matrix variables: the structure shown is the first one. Open the .log file of the job for every step.');
      }
    } else if (meta.jobTypes.indexOf('opt') >= 0 && ds.frames.length > 1) last(ds.frames).converged = optDone && !noConv;
    if (freqs.length) {
      const N = last(ds.frames).Z.length;
      const mv = modes.map(v => { const a = new Array(3 * N).fill(0); v.forEach((x, k) => { if (k < 3 * N && Number.isFinite(x)) a[k] = x; }); return a; });
      ds.vibrations = buildVib(ds, ds.frames.length - 1, freqs, mv, { ir: ir.some(Number.isFinite) ? ir : null, sym: sym.some(Boolean) ? sym : null }, false);
      addJT(meta, 'freq');
    }
    if (zpe !== null) { ds.thermo = newThermo(); ds.thermo.zpe = zpe; }
    if (orbA) ds.orbitals = { alpha: orbitalSet(orbA), beta: orbB ? orbitalSet(orbB) : null };
    if (ds.orbitals && !ds.orbitals.alpha) ds.orbitals = null;
    return finish(ds);
  }

  // =================================================================== Dalton
  const DALTON_SIG = /Dalton - An Electronic Structure Program|This is output from DALTON/i;
  /** A geometry block of Dalton ("Molecular geometry (au)", "Next geometry (au)", "Final geometry (au)"): a title,
   *  dashes, a blank line, then rows "C   _1   x y z" (or "C1  x y z") in Bohr. */
  function daltonGeom(L, i) {
    let j = i + 1;
    if (/^\s*-+\s*$/.test(L[j] || '')) j++;
    while (j < L.length && !L[j].trim()) j++;
    const Z = [], xyz = [];
    for (; j < L.length; j++) {
      const r = tok(L[j]);
      if (r.length < 4 || r.length > 5 || !r.slice(-3).every(isNum) || !/^[A-Za-z]/.test(r[0])) break;
      Z.push(elZ(r[0])); xyz.push(num(r[r.length - 3]) * BOHR, num(r[r.length - 2]) * BOHR, num(r[r.length - 1]) * BOHR);
    }
    if (Z.length && (j >= L.length || (j === L.length - 1 && L[j] === ''))) return { Z: [], xyz: [], end: j };
    return { Z, xyz, end: j };
  }
  /** "Cartesian Coordinates (a.u.)": rows "C   / 1 :   1  x  1.3498   2  y  2.3494   3  z  0.0000" (Bohr). */
  function daltonCart(L, i) {
    const Z = [], xyz = [];
    let j = i + 1;
    for (; j < L.length && j < i + 6 && !/Total number of coordinates/.test(L[j]); j++);
    for (j++; j < L.length; j++) {
      const m = /^\s*([A-Za-z][^:]*?)\s*:\s*\d+\s+x\s+(\S+)\s+\d+\s+y\s+(\S+)\s+\d+\s+z\s+(\S+)/.exec(L[j]);
      if (!m) { if (!L[j].trim() && !Z.length) continue; break; }
      Z.push(elZ(m[1].split(/[\s/]/)[0])); xyz.push(num(m[2]) * BOHR, num(m[3]) * BOHR, num(m[4]) * BOHR);
    }
    if (Z.length && (j >= L.length || (j === L.length - 1 && L[j] === ''))) return { Z: [], xyz: [], end: j };
    return { Z, xyz, end: j };
  }
  /** The molecule from the echo of the .mol file (when the output prints no coordinates): Angstrom or Bohr rows
   *  "C  x y z" under "Charge=6.0 Atoms=6" lines. The echo is cut at 50 columns, so the last column may be rounded. */
  function daltonMolEcho(echo) {
    const k = echo.findIndex(t => /^(Atomtypes|AtomTypes)\s*=/i.test(t) || /^ATOMBASIS|^BASIS|^INTGRL/i.test(t));
    if (k < 0) return null;
    const head = echo.find(t => /^Atomtypes\s*=/i.test(t)) || '';
    if (/Generators?\s*=/i.test(head) && !/Generators?\s*=\s*0/i.test(head)) return null;   // symmetry-unique atoms only
    const ang = /\bAngstrom\b/i.test(head);
    const Z = [], xyz = [];
    let cut = false, want = 0;
    for (const t of echo.slice(k)) {
      const r = tok(t), q = /\bAtoms\s*=\s*(\d+)/i.exec(t);
      if (q) want += int(q[1]);
      if (r.length === 4 && /^[A-Za-z]/.test(r[0]) && r.slice(1).every(isNum)) {
        const f = ang ? 1 : BOHR;
        Z.push(elZ(r[0])); xyz.push(num(r[1]) * f, num(r[2]) * f, num(r[3]) * f);
        if (t.length >= 49) cut = true;
      } else if (/^\*\*/.test(t)) break;
    }
    return Z.length && Z.length === want ? { Z, xyz, cut } : null;   // (every atom of the "Atoms=" counts)
  }
  function parseDalton(text, name) {
    const ds = BM.emptyDataset('dalton-out', name), meta = ds.meta, F = new Frames(ds), L = lines(text), n = L.length;
    ds.program = 'Dalton';
    const echo = [];
    let inEcho = false, echoDone = false, basis = null, method = null, normal = null, iterRows = null, conv = null;
    const pending = {};
    const energy = (kind, v) => {
      const f = F.cur();
      if (!f) { pending[kind] = v; return; }
      F.energy(kind, v, f);
    };
    const add = g => {
      if (!g || !g.Z.length) return;
      F.add(g.Z, g.xyz);
      for (const k of Object.keys(pending)) { energy(k, pending[k]); delete pending[k]; }
    };
    let freqRows = null, raman = null, nc = null, states = [], orbE = null, dip = null, nElec = NaN, zpe = null;
    for (let i = 0; i < n; i++) {
      const s = L[i], t = s.trim();
      if (!t) continue;
      let m;
      if (/^Content of the \.(dal|mol) (input )?file/.test(t)) { inEcho = true; continue; }
      if (inEcho) {
        if (/^\*{10,}/.test(t) || /Output from DALTON general input processing/.test(t)) { inEcho = false; echoDone = true; continue; }
        if (!/^-+$/.test(t)) echo.push(t);
        continue;
      }
      if ((m = /This is output from DALTON\s*\(?(?:Release\s+)?(?:Dalton)?\s*([\w.]+)/i.exec(t))) { ds.program = 'Dalton ' + m[1].replace(/^DALTON/i, ''); continue; }
      if (/^Cartesian Coordinates \(a\.u\.\)/.test(t)) { const g = daltonCart(L, i); add(g); i = g.end - 1; continue; }
      if (/^(Molecular|Next|Final) geometry \(au\)$/.test(t)) { const g = daltonGeom(L, i); add(g); i = g.end - 1; continue; }
      if ((m = /^@?\s*Final (HF|DFT|MCSCF|SCF) energy:\s*(-?\d+\.\d+)/.exec(t))) { energy(m[1] === 'MCSCF' ? 'MCSCF' : 'SCF', num(m[2])); if (!method && m[1] === 'HF') method = 'HF'; continue; }
      if ((m = /^@?\s*=\s*MP2 second order energy\s*:\s*(-?\d+\.\d+)/.exec(t))) { energy('MP2', num(m[1])); continue; }
      if ((m = /^Total (MP2|CCSD|CC2|CCS|CC3|CCSDT)\s+energy:\s*(-?\d+\.\d+)/.exec(t))) { energy(m[1], num(m[2])); continue; }
      if ((m = /^Total energy (CCSD\(T\)|CCSD\[T\]):\s*(-?\d+\.\d+)/.exec(t))) { if (m[1] === 'CCSD(T)') energy('CCSD(T)', num(m[2])); continue; }
      if ((m = /^@\s*Total charge of the molecule\s+(-?\d+)/.exec(t))) { meta.charge = int(m[1]); continue; }
      if ((m = /^@\s*Spin multiplicity(?: and 2 M_S)?:?\s+(\d+)/.exec(t))) { meta.multiplicity = int(m[1]); continue; }
      if ((m = /^Basis set used is "([^"]+)" from the basis set library/.exec(t))) { basis = m[1]; continue; }
      if ((m = /^Number of electrons :\s*(\d+)/.exec(t))) { nElec = int(m[1]); continue; }
      if (/^Sym\s+(Hartree-Fock|Kohn-Sham)\s+orbital energies/.test(t)) {
        // per symmetry: "1 Ag  e e e e e" and continuation lines; the occupied ones are the lowest of each symmetry
        orbE = [];
        let j = i + 1, cur = null;
        for (; j < n; j++) {
          const u = L[j].trim();
          if (!u) { if (cur && j + 1 < n && !L[j + 1].trim() && j + 2 < n && !L[j + 2].trim()) break; continue; }
          const r = tok(u);
          if (isInt(r[0]) && r.length >= 3 && /^[A-Za-z]/.test(r[1]) && r.slice(2).every(isNum)) { cur = []; orbE.push(cur); cur.push(...r.slice(2).map(num)); continue; }
          if (cur && r.every(isNum)) { cur.push(...r.map(num)); continue; }
          break;
        }
        i = j - 1;
        continue;
      }
      if (/^Dipole moment$/.test(t)) {
        for (let j = i + 1; j < n && j < i + 6; j++) {
          const r = tok(L[j]);
          if (r.length === 3 && r.every(isNum)) { dip = { total: num(r[0]) }; break; }   // au, Debye, C m
        }
        continue;
      }
      if (/^Dipole moment components$/.test(t) && dip) {
        const c = {};
        for (let j = i + 1; j < n && j < i + 10; j++) {
          const r = tok(L[j]);
          if (/zero by symmetry/.test(L[j])) { c.x = c.y = c.z = 0; break; }
          if (r.length >= 3 && /^[xyz]$/.test(r[0]) && isNum(r[1])) c[r[0]] = num(r[1]);   // (au)
        }
        dip.x = c.x || 0; dip.y = c.y || 0; dip.z = c.z || 0;
        continue;
      }
      if (/^Total Molecular Energy$/.test(t)) {   // electronic, vibrational (the zero-point energy) and total, in Hartrees
        for (let j = i + 1; j < n && j < i + 8; j++) {
          const r = tok(L[j]);
          if (r.length === 4 && r[3] === 'Hartrees' && r.slice(0, 3).every(isNum)) { zpe = { zpe: num(r[1]), eZPE: num(r[2]) }; break; }
        }
        continue;
      }
      if (/^@ Iter\s+Energy\s+Change\s+GradNorm/.test(t)) {
        iterRows = [];
        let j = i + 1;
        for (; j < n; j++) {
          const r = tok(L[j].replace(/^@/, '').replace(/\*/g, ''));
          if (!/^@/.test(L[j].trim())) break;
          if (r.length >= 6 && isInt(r[0]) && r.slice(1, 6).every(isNum)) iterRows.push(r.map(num));
          else if (iterRows.length && !/^-+$/.test(r.join(''))) break;
        }
        i = j - 1;
        continue;
      }
      if (/^@ Geometry converged in/.test(t)) { conv = true; continue; }
      if (/Geometry optimization did not converge|Maximum number of optimization iterations reached/i.test(t)) { conv = false; continue; }
      if (/^Vibrational Frequencies and IR Intensities$/.test(t)) {
        freqRows = [];
        let j = i + 1;
        while (j < n && j < i + 8 && !/\bcm-1\b/.test(L[j])) j++;   // the units line, then dashes, then the rows
        for (j++; j < n && /^\s*-+\s*$/.test(L[j]); j++);
        for (; j < n; j++) {
          const r = tok(L[j]);
          if (r.length < 5 || !isInt(r[0])) break;
          const f = /i$/.test(r[2]) ? -num(r[2].slice(0, -1)) : num(r[2]);
          freqRows.push({ f, sym: r[1], ir: num(r[4]) });
        }
        i = j - 1;
        continue;
      }
      if (/^Normal Coordinates \(bohrs\*amu\*\*\(1\/2\)\):/.test(t) && freqRows) {
        // blocks of five modes: per atom three rows "C   x  v v v v v", atoms in order
        nc = freqRows.map(() => []);
        let j = i + 1, cols = null, atomRow = 0;
        for (; j < n; j++) {
          const u = L[j].trim();
          if (!u || /^-+$/.test(u)) continue;
          const r = tok(u);
          if (r.length >= 2 && r.length % 2 === 0 && r.every(isNum) && isInt(r[0])) { cols = []; for (let k = 0; k < r.length; k += 2) cols.push(int(r[k]) - 1); atomRow = 0; continue; }
          if (cols && r.length === cols.length + 2 && /^[xyz]$/.test(r[1]) && r.slice(2).every(isNum)) {
            const ax = 'xyz'.indexOf(r[1]), a = Math.floor(atomRow / 3);
            r.slice(2).forEach((v, k) => { if (nc[cols[k]]) nc[cols[k]][3 * a + ax] = num(v); });
            atomRow++;
            continue;
          }
          break;
        }
        i = j - 1;
        continue;
      }
      if (/^Raman related properties for freq\.\s+0?\.0+\s+au/.test(t) && !raman) {
        raman = new Map();
        let j = i + 1;
        for (; j < n && j < i + 6 && !/^Mode\s+Freq\./.test(L[j].trim()); j++);
        for (j++; j < n; j++) {
          const r = tok(L[j]);
          if (!r.length) { if (raman.size) break; continue; }
          if (r.length < 7 || !isInt(r[0])) break;
          raman.set(int(r[0]), { act: num(r[4]) + num(r[5]), dep: num(r[6]) });
        }
        i = j - 1;
        continue;
      }
      if ((m = /^@ Excited state no:\s*(\d+)\s+in symmetry\s+\d+\s*\(\s*(\S+)\s*\)/.exec(t))) { states.push({ n: int(m[1]), sym: m[2], ev: NaN, f: 0 }); continue; }
      if ((m = /^@ Excitation energy\s*:\s*(-?\d+\.\d+)\s+au/.exec(t)) && states.length) { last(states).ev = num(m[1]) * HARTREE2EV; continue; }
      // (one line per dipole component, XDIPLEN, YDIPLEN, ZDIPLEN: the oscillator strength is their sum)
      if ((m = /^@ Oscillator strength \(LENGTH\)\s*:\s*(\S+)/.exec(t)) && states.length) { const v = num(m[1]); if (Number.isFinite(v)) last(states).f += v; continue; }
      if (/^(>{3,}\s*)?Total (CPU|wall)\s+time used in DALTON/.test(t)) { normal = true; continue; }
      if (/^--- SEVERE ERROR|^\*\*\* ERROR/.test(t)) normal = false;
    }
    // the input: basis, wave function, job type
    let dft = null;
    for (let k = 0; k < echo.length; k++) {
      const e = echo[k].toUpperCase();
      if (e === 'BASIS' && echo[k + 1] && !basis) basis = echo[k + 1].trim();
      if (e === '.DFT' && echo[k + 1]) dft = echo[k + 1].trim().toUpperCase();
      if (/^\.OPTIMIZE|^\.OPTIMI/.test(e)) addJT(meta, 'opt');
      if (/^\.VIBANA|^\.RAMAN/.test(e)) addJT(meta, 'freq');
      if (/^\.EXCITA|^\*\*RESPONSE|^\.RUN RESPONSE/.test(e)) addJT(meta, 'td');
      if (/^\.HF$/.test(e) && !method) method = 'HF';
      if (/^\.MP2$/.test(e)) method = 'MP2';
      if (/^\.CC$/.test(e)) method = method || 'CC';
      if (/^\.(CCSD\(T\)|CCSD|CC2|CC3)$/.test(e)) method = e.slice(1);
    }
    if (!ds.frames.length && echoDone) {
      const g = daltonMolEcho(echo);
      if (g) {
        add(g);
        if (g.cut) ds.warnings.push('This Dalton output prints the coordinates only in its echo of the .mol file, cut at 50 columns: the last coordinate may be rounded. Open the .mol file for the exact geometry.');
      }
    }
    if (!ds.frames.length) throw new Error('No geometry found in this Dalton output.');
    const top = last(ds.frames).energies || {};
    if (top['CCSD(T)'] !== undefined) method = 'CCSD(T)';
    else if (top.CCSD !== undefined && method !== 'CCSD(T)') method = method && /^CC/.test(method) ? method : 'CCSD';
    else if (top.MP2 !== undefined && !method) method = 'MP2';
    meta.method = method && method !== 'HF' ? method : dft || method || null;
    meta.basis = basis;
    meta.normalTermination = normal;
    // an optimisation keeps every geometry; other jobs (numerical Raman, response) show their reference geometry
    if (meta.jobTypes.indexOf('opt') < 0 && ds.frames.length > 1) ds.frames.length = 1;
    if (iterRows && iterRows.length) {
      addJT(meta, 'opt');
      // "@ Iter Energy Change GradNorm Index StepLen TrustRad #Rej": one row per geometry from the first
      const F2 = ds.frames;
      iterRows.forEach((r, k) => {
        const f = F2[k];
        if (!f) return;
        const o = f.opt = {};
        setOpt(o, 'energyChange', r[2]);
        setOpt(o, 'gradNorm', r[3]);
        setOpt(o, 'stepNorm', r[5]);
      });
      if (conv !== null) last(ds.frames).converged = conv;
    }
    if (freqRows && freqRows.length) {
      // Dalton lists the highest frequency first
      const N = last(ds.frames).Z.length;
      const order = freqRows.map((r, k) => k).sort((a, b) => freqRows[a].f - freqRows[b].f);
      const mv = order.map(k => { const a = new Array(3 * N).fill(0); if (nc && nc[k]) nc[k].forEach((v, q) => { if (q < 3 * N && Number.isFinite(v)) a[q] = v; }); return a; });
      const ram = raman ? order.map(k => (raman.has(k + 1) ? raman.get(k + 1).act : NaN)) : null;
      const dep = raman ? order.map(k => (raman.has(k + 1) ? raman.get(k + 1).dep : NaN)) : null;
      ds.vibrations = buildVib(ds, ds.frames.length - 1, order.map(k => freqRows[k].f), mv,
        { ir: order.map(k => freqRows[k].ir), sym: order.map(k => freqRows[k].sym) }, false);
      if (ds.vibrations && ram) { ds.vibrations.raman = ram; ds.vibrations.depol = dep; }
      addJT(meta, 'freq');
    }
    if (states.length) {
      const st = states.filter(x => Number.isFinite(x.ev)).sort((a, b) => a.ev - b.ev);
      ds.excitedStates = st.map((x, k) => state(k + 1, x.ev, x.f, x.sym + ', state ' + x.n));
      addJT(meta, 'td');
    }
    if (zpe) { ds.thermo = newThermo(); ds.thermo.zpe = zpe.zpe; ds.thermo.eZPE = zpe.eZPE; }
    if (dip && Number.isFinite(dip.total)) {
      ds.dipole = { x: (dip.x || 0) * AU2DEBYE, y: (dip.y || 0) * AU2DEBYE, z: (dip.z || 0) * AU2DEBYE, total: dip.total * AU2DEBYE };
    }
    if (orbE && orbE.length) {
      const all = [].concat(...orbE).sort((a, b) => a - b);
      const q = meta.charge || 0, mult = meta.multiplicity || 1;
      const Ztot = last(ds.frames).Z.reduce((a, z) => a + z, 0);
      const ne = Number.isFinite(nElec) && mult === 1 ? nElec : Ztot - q;
      const na = Math.round((ne + mult - 1) / 2), nb = Math.round((ne - mult + 1) / 2);
      ds.orbitals = { alpha: { energies: all, nOcc: na }, beta: mult > 1 ? { energies: all, nOcc: nb } : null };
    }
    return finish(ds);
  }

  // =================================================================== Turbomole program outputs
  // every program starts with "<name> (<host>) : TURBOMOLE ..." (older versions: "<name>(<host>) : TURBOMOLE V5-9-0")
  const TM_SIG = /^\s*([a-z][\w]*)\s*\([^)]*\)\s*:\s*TURBOMOLE\b/m;
  const TM_PROGRAMS = { dscf: 'SCF', ridft: 'SCF (RI)', escf: 'excited states', egrad: 'excited-state gradient', grad: 'gradient',
    rdgrad: 'gradient (RI)', force: 'vibrations', aoforce: 'vibrations', ricc2: 'CC2 / ADC(2)', ccsdf12: 'coupled cluster',
    rimp2: 'MP2 (RI)', mpgrad: 'MP2 gradient', statpt: 'optimisation step', relax: 'optimisation step', mpshift: 'NMR shieldings' };
  function parseTurbomoleOut(text, name) {
    const ds = BM.emptyDataset('turbomole-out', name), meta = ds.meta, F = new Frames(ds), L = lines(text), n = L.length;
    const progs = [];
    let version = null, basis = null, fnl = null, normal = null, uhf = false, unpaired = null, charge = null;
    const pending = {};
    // the frame energy is the highest level computed (an MP3 or CC2 run prints its MP2 energy first)
    const LEVELS = ['CCSD(T)', 'CCSD', 'CC2', 'ADC(2)', 'MP4', 'MP3', 'MP2', 'CIS(D)', 'CCS', 'SCF'];
    const energy = (kind, v) => {
      const f = F.cur();
      if (!f) { pending[kind] = v; return; }
      f.energies = f.energies || {};
      f.energies[kind] = v;
      const top = LEVELS.find(k => f.energies[k] !== undefined);
      f.energy = top ? f.energies[top] : v;
    };
    const freqs = [], ir = [], sym = [], modes = [];
    let es = null, esCur = null, ccRows = null, ccModel = null, ccOpen = false, ccF = new Map(), ccKey = null, zpe = null, dipAu = null;
    for (let i = 0; i < n; i++) {
      const s = L[i], t = s.trim();
      if (!t) continue;
      let m;
      if ((m = TM_SIG.exec(s)) && /TURBOMOLE/.test(s)) {
        progs.push(m[1]);
        const v = /TURBOMOLE\s+(?:rev\.\s*)?V?(\d+[.-]\d+(?:[.-]\d+)?)/.exec(s);
        if (v) version = v[1].replace(/-/g, '.').replace(/\.0$/, '');
        continue;
      }
      if (/^atomic coordinates\s+atom\b/.test(t)) {   // Bohr: x y z element [shells] charge [pseudo] isotop
        const Z = [], xyz = [];
        let j = i + 1;
        for (; j < n; j++) {
          const r = tok(L[j]);
          if (r.length < 5 || !r.slice(0, 3).every(isNum) || !/^[a-z]{1,2}$/i.test(r[3])) break;
          Z.push(elZ(r[3])); xyz.push(num(r[0]) * BOHR, num(r[1]) * BOHR, num(r[2]) * BOHR);
        }
        if (Z.length && j < n - 1) {
          F.add(Z, xyz);
          for (const k of Object.keys(pending)) { energy(k, pending[k]); delete pending[k]; }
        }
        i = j - 1;
        continue;
      }
      if ((m = /^\|?\s*total energy\s*=\s*(-?\d+\.\d+)/.exec(t))) { energy('SCF', num(m[1])); continue; }
      if ((m = /^\*?\s*Final (MP2|CC2|CCSD\(T\)|CCSD|ADC\(2\)|CIS\(D\)|CIS|MP3|MP4|CCS) energy\s*:\s*(-?\d+\.\d+)/.exec(t))) { energy(m[1], num(m[2])); continue; }
      if ((m = /^EMP2\s*:\s*(-?\d+\.\d+)/.exec(t))) { energy('MP2', num(m[1])); continue; }   // (rimp2)
      if (/^UHF mode switched on/.test(t)) { uhf = true; continue; }
      if (/^type\s+atoms\s+prim\s+cont\s+basis/.test(t)) {
        for (let j = i + 2; j < n; j++) {
          const r = tok(L[j]);
          if (r.length < 5 || !isInt(r[1])) break;
          if (!basis) basis = r[4];
        }
        continue;
      }
      if (/^density functional$/.test(t)) {   // the next line names it: "B3-LYP hybrid functional", "...: B3-LYP"
        for (let j = i + 1; j < n && j < i + 4; j++) {
          const u = L[j].trim();
          if (!u || /^-+$/.test(u)) continue;
          const q = /:\s*(\S+)\s*$/.exec(u) || /^(\S+)\s+(?:hybrid\s+|meta-GGA\s+|GGA\s+|LDA\s+)?functional/i.exec(u);
          if (q) fnl = q[1].replace(/-/g, '').toUpperCase().replace(/^BP86$|^BP$/, 'BP86');
          break;
        }
        continue;
      }
      if (/^charge$/.test(t)) {   // electrostatic moments: nuc, elec -> total (the molecular charge); "a-b" (unpaired electrons)
        for (let j = i + 1; j < n && j < i + 5; j++) {
          const r = tok(L[j]);
          if (r.length === 3 && r.every(isNum)) charge = Math.round(num(r[2])) || 0;
          else if (r[0] === 'a-b' && isNum(r[1])) unpaired = Math.round(Math.abs(num(r[1])));
        }
        continue;
      }
      if (/^dipole moment$/.test(t)) {
        const c = {};
        for (let j = i + 1; j < n && j < i + 6; j++) {
          const r = tok(L[j]);
          if (r.length === 4 && /^[xyz]$/.test(r[0]) && isNum(r[3])) c[r[0]] = num(r[3]);
        }
        if (c.x !== undefined && c.y !== undefined && c.z !== undefined) dipAu = c;
        continue;
      }
      // ---- escf: one block per state
      if ((m = /^(\d+)\s+(?:(singlet|doublet|triplet|quartet|quintet)\s+)?(\S+)\s+excitation$/i.exec(t))) {   // (no multiplicity: unrestricted)
        es = es || [];
        esCur = { n: int(m[1]), mult: m[2] ? m[2][0].toUpperCase() + m[2].slice(1).toLowerCase() : '', sym: m[3], ev: NaN, f: null, R: null, tr: [], part: null };
        es.push(esCur);
        continue;
      }
      if (esCur) {
        if ((m = /^Excitation energy \/ eV:\s*(-?\d+\.\d+(?:[EeDd][-+]?\d+)?)/.exec(t))) { esCur.ev = num(m[1]); continue; }
        if (/^Oscillator strength:/.test(t)) { esCur.part = 'f'; continue; }
        if (/^Rotatory strength:/.test(t)) { esCur.part = 'R'; continue; }
        if ((m = /^length representation:\s*(\S+)/.exec(t)) && esCur.part === 'f' && esCur.f === null) { esCur.f = num(m[1]); continue; }
        if ((m = /^length rep\. \/ 10\^\(-40\)erg\*cm\^3:\s*(\S+)/.exec(t)) && esCur.part === 'R') { esCur.R = num(m[1]); continue; }
        if (/^Dominant contributions:/.test(t)) { esCur.part = 'tr'; continue; }
        if (esCur.part === 'tr' && (m = /^(\d+)\s+(\S+)\s+-?\d+\.\d+\s+(\d+)\s+(\S+)\s+-?\d+\.\d+\s+(\d+\.\d+)$/.exec(t))) {
          esCur.tr.push({ from: m[1] + m[2], to: m[3] + m[4], coeff: Math.sqrt(num(m[5]) / 100) });
          continue;
        }
        if (/^Change of electron number|^Electric transition dipole/.test(t)) { esCur.part = null; continue; }
      }
      // ---- ricc2: the excitation-energy tables (the last model printed is the final one) and the transition strengths
      if ((m = /^\|\s*sym\s*\|\s*multi\s*\|\s*state\s*\|\s*(\S+)\s+excitation energies/.exec(t))) { ccModel = m[1]; ccRows = new Map(); ccOpen = true; continue; }
      if (ccOpen) {   // rows "| a | 1 | 1 | E (Hartree) | E (eV) | ..." until the table ends
        if ((m = /^\|\s*(\S+)\s*\|\s*(\d+|-)\s*\|\s*(\d+)\s*\|\s*(-?\d+\.\d+)\s*\|\s*(-?\d+\.\d+)\s*\|/.exec(t))) {   // (multiplicity "-": unrestricted)
          ccRows.set(m[1] + ' ' + m[2] + ' ' + m[3], { sym: m[1], mult: int(m[2]), n: int(m[3]), ev: num(m[5]) });
          continue;
        }
        if (!/^[|+]/.test(t) || (/^\+=+/.test(t) && ccRows.size)) ccOpen = false;   // (the border after the rows ends it)
      }
      if ((m = /number, symmetry, multiplicity:\s*(\d+)\s+(\S+)\s+(\d+|-)/.exec(t))) { ccKey = m[2] + ' ' + m[3] + ' ' + m[1]; continue; }
      if ((m = /^oscillator strength \(length gauge\)\s*:\s*(\S+)/.exec(t)) && ccKey) { ccF.set(ccKey, Math.max(0, num(m[1]))); ccKey = null; continue; }
      // ---- aoforce: blocks of six modes
      if (/^mode\s+\d+(\s+\d+)*$/.test(t) && progs.some(p => /force/.test(p))) {
        let cols = null, j = i + 1, atom = -1;
        for (; j < n; j++) {
          const u = L[j], v = u.trim();
          if (!v) continue;
          let q;
          if ((q = /^frequency\s+(.*)$/.exec(v))) {
            const fv = tok(q[1]).map(x => (/^i/.test(x) ? -num(x.slice(1)) : num(x)));
            cols = fv.map((x, k) => freqs.length + k);
            fv.forEach(x => { freqs.push(x); ir.push(NaN); sym.push(''); modes.push([]); });
            continue;
          }
          if (!cols) break;
          if ((q = /^symmetry\s+(.*)$/.exec(v))) { tok(q[1]).forEach((x, k) => { if (k < cols.length) sym[cols[k]] = x; }); continue; }
          if ((q = /^intensity \(km\/mol\)\s+(.*)$/.exec(v))) { tok(q[1]).forEach((x, k) => { if (k < cols.length) ir[cols[k]] = num(x); }); continue; }
          if (/^(IR|RAMAN|\|dDIP\/dQ\||intensity \(\s*%\s*\))/.test(v)) continue;
          if (/^symmetry$/.test(v)) continue;
          if ((q = /^(\d+)\s+[a-z]{1,2}\s+([xyz])\s+(.*)$/i.exec(v))) { atom = int(q[1]) - 1; const ax = 'xyz'.indexOf(q[2]); tok(q[3]).forEach((x, k) => { if (k < cols.length) modes[cols[k]][3 * atom + ax] = num(x); }); continue; }
          if ((q = /^([xyz])\s+(.*)$/.exec(v)) && atom >= 0) { const ax = 'xyz'.indexOf(q[1]); tok(q[2]).forEach((x, k) => { if (k < cols.length) modes[cols[k]][3 * atom + ax] = num(x); }); continue; }
          if (/^reduced mass/.test(v)) { j++; break; }
          break;
        }
        i = j - 1;
        continue;
      }
      if ((m = /zero point VIBRATIONAL energy\s*:\s*(-?\d+\.\d+)\s*Hartree/i.exec(t))) { zpe = num(m[1]); continue; }
      if (/ended normally\s*$/.test(t) || /^(\*{4}\s+)?\w+\s*:\s*all done(\s+\*{4})?$/.test(t)) { normal = true; continue; }
      if (/ended abnormally|^ERROR:|error termination/i.test(t)) normal = false;
    }
    const uniq = [...new Set(progs)];
    ds.program = 'Turbomole' + (version ? ' ' + version : '') + (uniq.length === 1 ? ' (' + uniq[0] + ')' : '');
    if (!ds.frames.length) {
      throw new Error(uniq.length ? 'This Turbomole output (' + uniq.join(', ') + ') prints no atomic coordinates: open the coord file of the job next to it.'
        : 'No atomic coordinates found in this Turbomole output.');
    }
    meta.charge = charge;
    meta.multiplicity = unpaired !== null ? unpaired + 1 : (charge !== null ? 1 : null);
    meta.basis = basis;
    meta.normalTermination = normal;
    if (uhf && unpaired === null) meta.multiplicity = null;
    if (freqs.length) {
      const keep = freqs.map((f, k) => k).filter(k => freqs[k] !== 0);   // (0.00: translations and rotations)
      const N = last(ds.frames).Z.length;
      ds.vibrations = buildVib(ds, ds.frames.length - 1, keep.map(k => freqs[k]),
        keep.map(k => { const a = new Array(3 * N).fill(0); modes[k].forEach((v, q) => { if (q < 3 * N && Number.isFinite(v)) a[q] = v; }); return a; }),
        { ir: keep.map(k => ir[k]), sym: keep.some(k => sym[k]) ? keep.map(k => sym[k]) : null }, false);
      addJT(meta, 'freq');
    }
    if (zpe !== null) { ds.thermo = newThermo(); ds.thermo.zpe = zpe; }
    if (es && es.length) {
      const list = es.filter(x => Number.isFinite(x.ev)).sort((a, b) => a.ev - b.ev);
      ds.excitedStates = list.map((x, k) => { const st = state(k + 1, x.ev, x.f, (x.mult + ' ' + x.sym).trim()); st.transitions = x.tr; if (Number.isFinite(x.R)) st.R = x.R; return st; });
      addJT(meta, 'td');
    } else if (ccRows && ccRows.size) {
      const MULT = { 1: 'Singlet', 2: 'Doublet', 3: 'Triplet' };
      const list = [...ccRows.entries()].map(([k, x]) => Object.assign({ key: k }, x)).sort((a, b) => a.ev - b.ev);
      ds.excitedStates = list.map((x, k) => state(k + 1, x.ev, ccF.has(x.key) ? ccF.get(x.key) : null, ((MULT[x.mult] || '') + ' ' + x.sym).trim() + ' (' + ccModel + ')'));
      addJT(meta, 'td');
    }
    if (dipAu) ds.dipole = { x: dipAu.x * AU2DEBYE, y: dipAu.y * AU2DEBYE, z: dipAu.z * AU2DEBYE, total: Math.hypot(dipAu.x, dipAu.y, dipAu.z) * AU2DEBYE };
    const top = last(ds.frames).energies || {};
    meta.method = ['CCSD(T)', 'CCSD', 'CC2', 'ADC(2)', 'MP4', 'MP3', 'MP2', 'CIS(D)'].find(k => top[k] !== undefined) || fnl || (top.SCF !== undefined ? 'HF' : null);
    if (progs.some(p => p === 'statpt' || p === 'relax') || ds.frames.length > 1) addJT(meta, 'opt');
    return finish(ds);
  }

  // =================================================================== GAMESS punch files (.dat)
  // GAMESS writes the symmetry-unique atoms ($DATA, and each step of an optimisation); all atoms only in the AIMPAC
  // section (Bohr), when it was asked for. The .dat extension is shared with Psi4 inputs and LAMMPS data files.
  const PUNCH_SIG = /^\s*-{3,}\s*(CLOSED|OPEN) SHELL ORBITALS\s*-{3,}|^-{10,}\s*DATA FROM NSERCH=|^-{5}\s*RESULTS FROM SUCCESSFUL|^E\([A-Z0-9-]+\)=\s*-?\d+\.\d+,\s*E\(NUC\)=/m;
  function punchUnique(L, j) {   // rows "C  6.0  x y z" under "COORDINATES OF SYMMETRY UNIQUE ATOMS (ANGS)"
    while (j < L.length && !/^\s*-{20,}/.test(L[j])) j++;
    const Z = [], xyz = [];
    for (j++; j < L.length; j++) {
      const r = tok(L[j]);
      if (r.length !== 5 || !r.slice(1).every(isNum)) break;
      Z.push(Math.round(num(r[1])) || elZ(r[0])); xyz.push(num(r[2]), num(r[3]), num(r[4]));
    }
    if (Z.length && (j >= L.length || (j === L.length - 1 && L[j] === ''))) return { Z: [], xyz: [], end: j };   // (cut while written)
    return { Z, xyz, end: j };
  }
  function parseGamessPunch(text, name) {
    const ds = BM.emptyDataset('gamess-punch', name), meta = ds.meta, L = lines(text), n = L.length;
    ds.program = 'GAMESS (punch file)';
    let group = 'C1', naxis = 1, data = null, title = null, aim = null, aimE = null, orb = [], molplt = null, masses = null;
    const steps = [];   // {Z, xyz (unique, Angstrom), energy, label}
    let curStep = null, method = null, finalSeen = false, eLast = null, inModes = false, dipole = null;
    for (let i = 0; i < n; i++) {
      const s = L[i], t = s.trim();
      if (!t) continue;
      let m;
      if (/^\$DATA\b/i.test(t) && !data) {   // title, group (and NAXIS), a blank line unless C1, atoms with their basis lines
        title = (L[i + 1] || '').trim();
        const g = tok(L[i + 2] || '');
        group = (g[0] || 'C1').toUpperCase(); naxis = int(g[1]) || 1;
        data = { Z: [], xyz: [] };
        let j = i + 3;
        for (; j < n && !/^\s*\$END\b/i.test(L[j]); j++) {
          const r = tok(L[j]);
          if (r.length === 5 && r.slice(1).every(isNum) && /^[A-Za-z]/.test(r[0])) {
            data.Z.push(Math.round(num(r[1])) || elZ(r[0])); data.xyz.push(num(r[2]), num(r[3]), num(r[4]));
          }
        }
        if (j >= n) data = { Z: [], xyz: [] };   // no $END: the group was cut
        i = j;
        continue;
      }
      if ((m = /^-+\s*DATA FROM NSERCH=\s*(\d+)/.exec(t))) { curStep = { Z: null, xyz: null, energy: null, label: 'Step ' + (int(m[1]) + 1) }; steps.push(curStep); continue; }
      if (/^-+\s*RESULTS FROM SUCCESSFUL/.test(t)) { curStep = { Z: null, xyz: null, energy: null, label: 'Optimised geometry', final: true }; steps.push(curStep); finalSeen = true; continue; }
      if (/^COORDINATES OF SYMMETRY UNIQUE ATOMS/.test(t) && curStep) {
        const u = punchUnique(L, i + 1);
        if (u.Z.length) { curStep.Z = u.Z; curStep.xyz = u.xyz; }
        i = u.end - 1;
        continue;
      }
      if ((m = /^E\(([^)]+)\)=\s*(-?\d+\.\d+)/.exec(t))) {   // E(R-B3LYP)= ..., E(NUC)= ..., 17 ITERS
        const v = num(m[2]);
        if (!method) method = m[1].replace(/^(R|U|RO)-?/, '').replace(/^HF$/, 'HF') || null;
        if (v !== 0) { eLast = v; if (curStep && curStep.energy === null) curStep.energy = v; }   // (UHF natural orbitals: E = 0)
        continue;
      }
      if ((m = /^E=\s*(-?\d+\.\d+)/.exec(t))) { const v = num(m[1]); if (v !== 0) { eLast = v; if (curStep && curStep.energy === null) curStep.energy = v; } continue; }
      if ((m = /^DIPOLE\s+(\S+)\s+(\S+)\s+(\S+)\s*$/.exec(t))) {   // (Debye, as in the output)
        const x = num(m[1]), y = num(m[2]), z = num(m[3]);
        if ([x, y, z].every(Number.isFinite)) dipole = { x, y, z, total: Math.hypot(x, y, z) };
        continue;
      }
      if (/^-+\s*START OF NORMAL MODES FOR -MOLPLT- PROGRAM/.test(t)) { molplt = []; inModes = true; continue; }
      if (inModes) {
        if (/^-+\s*END OF NORMAL MODES/.test(t)) { inModes = false; continue; }
        if (/^ATOMIC MASSES/.test(t)) { masses = []; continue; }
        if ((m = /^MODE\s+(\d+)\s+FREQUENCY=\s*(-?\d+\.\d+)/.exec(t))) { molplt.push({ f: num(m[2]), v: [] }); continue; }
        const r = tok(t);
        if (r.every(isNum)) { if (molplt.length) last(molplt).v.push(...r.map(num)); else if (masses) masses.push(...r.map(num)); }
        continue;
      }
      if (/^-+\s*TOP OF INPUT FILE FOR BADER'S AIMPAC PROGRAM/.test(t)) { aim = { Z: [], xyz: [], n: NaN }; continue; }
      if (aim && (m = /^GAUSSIAN\s+\d+\s+MOL ORBITALS\s+\d+\s+PRIMITIVES\s+(\d+)\s+NUCLEI/.exec(t))) { aim.n = int(m[1]); continue; }
      if (aim && (m = /^\s*([A-Za-z]{1,2})\s+\d+\s+\(CENTRE\s+\d+\)\s+(\S+)\s+(\S+)\s+(\S+)\s+CHARGE\s*=\s*(\S+)/.exec(s))) {
        aim.Z.push(Math.round(num(m[5])) || elZ(m[1])); aim.xyz.push(num(m[2]) * BOHR, num(m[3]) * BOHR, num(m[4]) * BOHR);
        continue;
      }
      if (aim && (m = /^MO\s+\d+\s+OCC NO =\s*(\S+)\s+ORB\. ENERGY =\s*(\S+)/.exec(t))) { orb.push({ occ: num(m[1]), e: num(m[2]) }); continue; }
      if (aim && (m = /^\S+\s+ENERGY =\s*(-?\d+\.\d+)/.exec(t))) { aimE = num(m[1]); continue; }
    }
    if (aim && Number.isFinite(aim.n) && aim.Z.length !== aim.n) aim.Z = [];   // (cut before its last atom)
    if (aim && !aim.Z.length) aim = null;
    if (!data && !aim) throw new Error('No $DATA group found in this GAMESS punch file.');
    ds.title = title || ds.title;
    // symmetry: expand the unique atoms with the group's operations (C1 and "CN 1" have none)
    const c1 = group === 'C1' || (group === 'CN' && naxis === 1);
    const ops = c1 ? null : gamessOps(group, naxis);
    const expand = g => (ops ? gamessExpand(g.Z, g.xyz, ops) : { Z: g.Z, xyz: g.xyz });
    let geoms = steps.filter(st => st.Z).map(st => Object.assign(expand(st), { energy: st.energy, label: st.label, final: !!st.final }));
    if (!geoms.length && data && data.Z.length && (c1 || ops)) geoms = [Object.assign(expand(data), { energy: eLast, label: null })];
    if (!c1 && !ops && !aim) throw new Error('The point group ' + group + ' of this punch file is not supported, and the file has no AIMPAC section with all atoms.');
    // GAMESS's own atom order, from the AIMPAC section (all atoms, final geometry): reorder the expanded geometries to it
    let ordered = !ops;
    if (aim && aim.Z.length) {
      const fin = geoms.length ? last(geoms) : null;
      let perm = null;
      if (fin && fin.Z.length === aim.Z.length) {
        perm = aim.Z.map((z, a) => {
          let best = -1, bd = 0.02;
          for (let b = 0; b < fin.Z.length; b++) {
            if (fin.Z[b] !== z) continue;
            const d = Math.hypot(fin.xyz[3 * b] - aim.xyz[3 * a], fin.xyz[3 * b + 1] - aim.xyz[3 * a + 1], fin.xyz[3 * b + 2] - aim.xyz[3 * a + 2]);
            if (d < bd) { bd = d; best = b; }
          }
          return best;
        });
        if (perm.some(b => b < 0) || new Set(perm).size !== perm.length) perm = null;
      }
      if (perm) { geoms = geoms.map(g => Object.assign({}, g, { Z: perm.map(b => g.Z[b]), xyz: [].concat(...perm.map(b => [g.xyz[3 * b], g.xyz[3 * b + 1], g.xyz[3 * b + 2]])) })); ordered = true; }
      else if (!geoms.length || (fin && fin.Z.length !== aim.Z.length)) { geoms = [{ Z: aim.Z, xyz: aim.xyz, energy: aimE !== null ? aimE : eLast, label: null }]; ordered = true; }
    }
    if (!geoms.length) throw new Error('No atoms found in this GAMESS punch file.');
    const F = new Frames(ds);
    for (const g of geoms) {
      const f = F.add(g.Z, g.xyz, g.label);
      if (Number.isFinite(g.energy)) F.energy('SCF', g.energy, f);
      if (g.final) f.converged = true;
    }
    if (steps.length) addJT(meta, 'opt');
    if (finalSeen) last(ds.frames).converged = true;
    meta.method = method;
    // vibrations (a Hessian job): every mode with its frequency, translations and rotations included. GAMESS orders
    // them by eigenvalue, so the modes before the six (five) of translation and rotation are imaginary.
    if (molplt && molplt.length && !ordered) {
      ds.warnings.push('The normal modes are not shown: they follow GAMESS\'s atom order, which this punch file does not give for the atoms generated by symmetry. Open the output for the vibrations.');
    } else if (molplt && molplt.length) {
      const fr = last(ds.frames), N = fr.Z.length, ntr = isLinear(fr.xyz) ? 5 : 6;
      const ok = molplt.filter(md => md.v.length === 3 * N);
      if (ok.length >= ntr + 1) {
        const tr = ok.map((md, k) => k).sort((a, b) => Math.abs(ok[a].f) - Math.abs(ok[b].f)).slice(0, ntr);
        const first = Math.min(...tr), drop = new Set(tr);
        const keep = ok.map((md, k) => k).filter(k => !drop.has(k));
        const fq = keep.map(k => (k < first ? -Math.abs(ok[k].f) : ok[k].f));
        const order = keep.map((k, q) => q).sort((a, b) => fq[a] - fq[b]);
        ds.vibrations = buildVib(ds, ds.frames.length - 1, order.map(q => fq[q]), order.map(q => ok[keep[q]].v), {}, false);
        addJT(meta, 'freq');
      }
    }
    if (dipole) ds.dipole = dipole;
    if (orb.length && orb.some(o => o.e !== 0)) {
      const e = orb.map(o => o.e).sort((a, b) => a - b);
      ds.orbitals = { alpha: { energies: e, nOcc: orb.filter(o => o.occ > 0.5).length }, beta: null };
    }
    if (!ordered) ds.warnings.push('The atoms of this ' + group + (/N/.test(group) ? ' (n = ' + naxis + ')' : '') + ' punch file were generated from its symmetry-unique ones: they may be numbered differently from the GAMESS output.');
    return finish(ds);
  }

  reg('molpro-out', 'Molpro output', ['out', 'log', 'pro'], 'output',
    h => (MOLPRO_SIG.test(h) ? 0.95 : (MOLPRO_LOG.test(h) || /^\s*Cartesian Coordinates are in \[bohr\]\s*$/m.test(h.slice(0, 200))) && /ATOMIC COORDINATES|Optimization point \d/.test(h) ? 0.9 : 0), parseMolpro);
  reg('dalton-out', 'Dalton output', ['out', 'log'], 'output', h => (DALTON_SIG.test(h) ? 0.95 : 0), parseDalton);
  reg('turbomole-out', 'Turbomole program output', ['out', 'log', 'last'], 'output', h => (TM_SIG.test(h) && /TURBOMOLE/.test(h) ? 0.95 : 0), parseTurbomoleOut);
  // (a punch file has no $CONTRL group: that is a GAMESS input)
  reg('gamess-punch', 'GAMESS punch file (.dat)', ['dat', 'pun'], 'output',
    h => (/^\s*\$DATA\b/im.test(h) && !/^\s*\$CONTRL\b/im.test(h) && PUNCH_SIG.test(h) ? 0.9 : 0), parseGamessPunch);

  BM.qcReaders.molpro = parseMolpro;
  BM.qcReaders.dalton = parseDalton;
  BM.qcReaders.turbomole = parseTurbomoleOut;
  BM.qcReaders.gamessPunch = parseGamessPunch;
})(typeof window !== 'undefined' ? (window.BM = window.BM || {}) : (globalThis.BM = globalThis.BM || {}));
