/* BlueMolecule — exporters.
 *
 * BM.writers.xyz({Z, xyz, comment})                         -> string
 * BM.writers.xyzMulti([{Z, xyz, comment}, ...])             -> string (multi-frame XYZ)
 * BM.writers.gaussianInput({Z, xyz, route, title, charge, multiplicity, link0, cell})
 * BM.writers.orcaInput({Z, xyz, keywords, charge, multiplicity, nprocs, maxcore, extraBlocks})
 * BM.writers.mol({Z, xyz, bonds, title})                    -> MDL V2000 molfile
 * BM.writers.pdb({Z, xyz, atomNames, resNames, resIds, cell})
 *
 * Z: atomic numbers (0 = dummy/ghost), xyz: flat array (Float64Array or Array) in Angstrom.
 * All writers end the text with a newline.
 *
 * Quantum-chemistry job inputs, all taking
 *   opts = {Z, xyz, charge = 0, multiplicity = 1, method, basis, task: 'sp'|'opt'|'freq'|'optfreq',
 *           nprocs, memoryMB (total), title, solvent (implicit-solvent name), extra (raw text), baseName}
 * BM.writers.psi4(opts) / qchem / nwchem / gamess / mopac / molpro -> input text
 * BM.writers.turbomole(opts)  -> `coord` file ($coord in bohr ... $end)
 * BM.writers.xtb(opts)        -> {files: [{name, text}], command}
 * Coordinates are written with 10 decimals (Å) or 14 decimals (bohr). Dummy atoms (Z = 0) are left
 * out by these writers. `extra` goes where each program expects additional input: Psi4 before the
 * energy/optimize call, Q-Chem $rem lines (or whole sections if it starts with "$"), NWChem before the
 * task lines, GAMESS groups before $DATA, MOPAC keywords, xtb command-line flags (or an xcontrol file if
 * it starts with "$"), Molpro before the method, Turbomole data groups before $end.
 *
 * BM.inputPrograms: [{id, name, ext, defaultMethod, defaultBasis, methods, bases, tasks,
 *                     write(opts) -> {text, files?, command?, note?}}]
 *   ids: gaussian, orca, psi4, qchem, nwchem, gamess, mopac, xtb, molpro, turbomole.
 *   gaussian/orca are built on gaussianInput/orcaInput (task keywords, solvent, %nprocshared/%mem,
 *   %pal/%maxcore per core; their `extra` goes after the geometry / as % blocks).
 */
(function (BM) {
  'use strict';

  function sym(z, dummy) { return z > 0 ? BM.el.symbol(z) : (dummy || 'X'); }
  function fx(x, w, d) {
    const s = Number.isFinite(x) ? (Math.abs(x) < 0.5 * Math.pow(10, -d) ? 0 : x).toFixed(d) : 'NaN';
    return s.length >= w ? s : ' '.repeat(w - s.length) + s;
  }
  function ix(v, w) { const s = String(v); return s.length >= w ? s : ' '.repeat(w - s.length) + s; }
  function oneLine(s) { return String(s === undefined || s === null ? '' : s).replace(/[\r\n]+/g, ' ').trim(); }
  function check(m) {
    if (!m || !m.Z || !m.xyz || m.xyz.length < 3 * m.Z.length) throw new Error('Writer needs Z and xyz arrays of matching size.');
  }
  /** fx for free-format (space-separated) files: a number wider than its field still gets a separating space. */
  function fs(x, w, d) { const s = fx(x, w, d); return s[0] === ' ' ? s : ' ' + s; }
  function atomLines(Z, xyz, dummy, w, d) {
    const out = [];
    for (let i = 0; i < Z.length; i++) {
      out.push(sym(Z[i], dummy).padEnd(3) + fs(xyz[3 * i], w, d) + fs(xyz[3 * i + 1], w, d) + fs(xyz[3 * i + 2], w, d));
    }
    return out;
  }

  function xyz(m) {
    check(m);
    return [String(m.Z.length), oneLine(m.comment)].concat(atomLines(m.Z, m.xyz, 'X', 15, 8)).join('\n') + '\n';
  }

  function xyzMulti(frames) {
    if (!Array.isArray(frames) || !frames.length) throw new Error('xyzMulti needs at least one frame.');
    return frames.map(xyz).join('');
  }

  function gaussianInput(m) {
    check(m);
    const lines = [];
    for (let l of m.link0 || []) { l = oneLine(l); if (l) lines.push(l[0] === '%' ? l : '%' + l); }
    let route = String(m.route || '#p sp hf/sto-3g').trim();
    const rl = route.split(/\r?\n/).map(s => s.trim()).filter(Boolean);
    if (rl[0][0] !== '#') rl[0] = '# ' + rl[0];
    lines.push(...rl, '', oneLine(m.title) || 'Title', '');
    lines.push((m.charge === undefined || m.charge === null ? 0 : m.charge) + ' ' +
      (m.multiplicity === undefined || m.multiplicity === null ? 1 : m.multiplicity));
    lines.push(...atomLines(m.Z, m.xyz, 'Bq', 15, 8));
    if (m.cell) for (const v of m.cell) if (Math.hypot(v[0], v[1], v[2]) > 1e-8) lines.push('Tv ' + fs(v[0], 15, 8) + fs(v[1], 15, 8) + fs(v[2], 15, 8));
    return lines.join('\n') + '\n\n\n';   // blank line terminates the molecule specification
  }

  function orcaInput(m) {
    check(m);
    let kw = Array.isArray(m.keywords) ? m.keywords.join(' ') : String(m.keywords || 'HF def2-SVP');
    const L = kw.split(/\r?\n/).map(s => s.trim()).filter(Boolean).map(s => (s[0] === '!' ? s : '! ' + s));
    if (m.nprocs && m.nprocs > 1) L.push('%pal nprocs ' + Math.round(m.nprocs) + ' end');
    if (m.maxcore) L.push('%maxcore ' + Math.round(m.maxcore));
    const extra = Array.isArray(m.extraBlocks) ? m.extraBlocks : (m.extraBlocks ? [m.extraBlocks] : []);
    for (const b of extra) { const s = String(b).replace(/\s+$/, ''); if (s) L.push(s); }
    L.push('', '* xyz ' + (m.charge === undefined || m.charge === null ? 0 : m.charge) + ' ' +
      (m.multiplicity === undefined || m.multiplicity === null ? 1 : m.multiplicity));
    for (const s of atomLines(m.Z, m.xyz, 'DA', 15, 8)) L.push('  ' + s);
    L.push('*');
    return L.join('\n') + '\n';
  }

  function guessBonds(Z, X) {
    const out = [];
    const n = Z.length;
    for (let i = 0; i < n; i++) {
      for (let j = i + 1; j < n; j++) {
        if (!Z[i] || !Z[j]) continue;
        const d = Math.hypot(X[3 * i] - X[3 * j], X[3 * i + 1] - X[3 * j + 1], X[3 * i + 2] - X[3 * j + 2]);
        if (d > 0.1 && d < 1.15 * (BM.el.covalent[Z[i]] + BM.el.covalent[Z[j]])) out.push([i, j, 1]);
      }
    }
    return out;
  }

  function mol(m) {
    check(m);
    const n = m.Z.length;
    if (n > 999) throw new Error('V2000 molfiles hold at most 999 atoms.');
    const bonds = m.bonds || guessBonds(m.Z, m.xyz);
    if (bonds.length > 999) throw new Error('V2000 molfiles hold at most 999 bonds.');
    outOfRange(m, -9999.9999, 99999.9999, 'Molfile');
    const d = new Date();
    const p2 = v => String(v).padStart(2, '0');
    const stamp = p2(d.getMonth() + 1) + p2(d.getDate()) + p2(d.getFullYear() % 100) + p2(d.getHours()) + p2(d.getMinutes());
    const L = [oneLine(m.title).slice(0, 80), '  BlueMol ' + stamp + '3D', ''];
    L.push(ix(n, 3) + ix(bonds.length, 3) + '  0  0  0  0  0  0  0  0999 V2000');
    for (let i = 0; i < n; i++) {
      L.push(fx(m.xyz[3 * i], 10, 4) + fx(m.xyz[3 * i + 1], 10, 4) + fx(m.xyz[3 * i + 2], 10, 4) + ' ' +
        (m.Z[i] > 0 ? BM.el.symbol(m.Z[i]) : 'Q').padEnd(3) + ' 0  0  0  0  0  0  0  0  0  0  0  0');
    }
    for (const b of bonds) {
      const o = b[2] === 1.5 ? 4 : Math.max(1, Math.min(3, Math.round(b[2] || 1)));
      L.push(ix(b[0] + 1, 3) + ix(b[1] + 1, 3) + ix(o, 3) + '  0');
    }
    L.push('M  END');
    return L.join('\n') + '\n';
  }

  const PROTEIN = new Set('ALA ARG ASN ASP CYS GLN GLU GLY HIS ILE LEU LYS MET PHE PRO SER THR TRP TYR VAL'.split(' '));

  /** Fixed-column formats cannot hold every coordinate: refuse rather than write fields that run together. */
  function outOfRange(m, lo, hi, what) {
    for (let k = 0; k < 3 * m.Z.length; k++) {
      const v = m.xyz[k];
      if (Number.isFinite(v) && (v < lo || v > hi)) {
        throw new Error(what + ' coordinates must lie between ' + lo + ' and ' + hi + ' Å (atom ' + (Math.floor(k / 3) + 1) + ' is at ' + v.toFixed(3) + '); use XYZ or another free-format file.');
      }
    }
  }

  function pdb(m) {
    check(m);
    outOfRange(m, -999.999, 9999.999, 'PDB');
    const L = ['REMARK   1 WRITTEN BY BLUEMOLECULE'];
    if (m.cell) {
      const len = v => Math.hypot(v[0], v[1], v[2]);
      const ang = (u, v) => Math.acos(Math.max(-1, Math.min(1, (u[0] * v[0] + u[1] * v[1] + u[2] * v[2]) / (len(u) * len(v) || 1)))) * 180 / Math.PI;
      const [a, b, c] = m.cell;
      if (len(a) > 0 && len(b) > 0 && len(c) > 0) {
        L.push('CRYST1' + fx(len(a), 9, 3) + fx(len(b), 9, 3) + fx(len(c), 9, 3) + fx(ang(b, c), 7, 2) +
          fx(ang(a, c), 7, 2) + fx(ang(a, b), 7, 2) + ' P 1           1');
      }
    }
    const counts = {};
    for (let i = 0; i < m.Z.length; i++) {
      const el = m.Z[i] > 0 ? BM.el.symbol(m.Z[i]) : 'X';
      counts[el] = (counts[el] || 0) + 1;
      let name = (m.atomNames && m.atomNames[i]) || (el + counts[el]);
      name = String(name).slice(0, 4);
      const nm4 = name.length < 4 && el.length === 1 ? (' ' + name).padEnd(4) : name.padEnd(4);
      const res = String((m.resNames && m.resNames[i]) || 'MOL').slice(0, 4);
      const rid = (m.resIds && Number.isFinite(m.resIds[i]) ? m.resIds[i] : 1) % 10000;
      const rec = PROTEIN.has(res.toUpperCase()) ? 'ATOM  ' : 'HETATM';
      L.push(rec + ix((i + 1) % 100000, 5) + ' ' + nm4 + ' ' + res.padStart(3).padEnd(4) + 'A' + ix(rid, 4) + '    ' +
        fx(m.xyz[3 * i], 8, 3) + fx(m.xyz[3 * i + 1], 8, 3) + fx(m.xyz[3 * i + 2], 8, 3) + '  1.00  0.00          ' +
        el.toUpperCase().padStart(2) + '  ');
    }
    L.push('END');
    return L.join('\n') + '\n';
  }

  // ================================================================== quantum-chemistry job inputs
  const BOHR = 0.529177210903;             // Å per bohr (CODATA 2018, as in the parsers)
  const TASKS = ['sp', 'opt', 'freq', 'optfreq'];
  const TASK_ALIAS = { sp: 'sp', energy: 'sp', singlepoint: 'sp', opt: 'opt', optimize: 'opt', optimise: 'opt',
    optimization: 'opt', optimisation: 'opt', freq: 'freq', frequency: 'freq', frequencies: 'freq', hessian: 'freq',
    optfreq: 'optfreq', optfrequency: 'optfreq', optfrequencies: 'optfreq' };
  const MULT_NAMES = ['', 'SINGLET', 'DOUBLET', 'TRIPLET', 'QUARTET', 'QUINTET', 'SEXTET', 'SEPTET', 'OCTET', 'NONET'];
  const SEMIEMPIRICAL = /^(pm[367](-d3h4x?|-d3|-org|-ts)?|am1|rm1|pm3|mndod?|pddg|zindo(\/s)?|om[123])$/i;
  const XTB_METHOD = /^(gfn[012](-?xtb)?|gfn-?ff|xtb[012]?)$/i;
  const HF_RE = /^(r|u|ro)?(hf|scf)$/i;
  const MP2_RE = /^(r|u)?mp2$/i;
  const CC_RE = /^ccsd(\(t\))?$/i;
  const WAVEFUNCTION = /^(r|u|ro)?(hf|scf|mp[234]|mp4\(sdq\)|ccsd|ccsd\(t\)|qcisd|qcisd\(t\)|cisd|cis|casscf|mrci|caspt2|nevpt2|ri-mp2|scs-mp2|lmp2|dlpno-mp2|dlpno-ccsd\(t\)|eom-ccsd)$/i;
  const methodHead = s => String(s || '').trim().split(/\s+/)[0];
  /** Composite (-3c), semi-empirical and tight-binding methods carry their own basis. */
  const needsBasis = s => { const h = methodHead(s); return !(SEMIEMPIRICAL.test(h) || XTB_METHOD.test(h) || /-3c$/i.test(h)); };
  const solventKey = s => String(s).toLowerCase().replace(/[\s,_()-]/g, '');
  const F = v => fx(v, 16, 10);
  const jobAtomLine = a => a.s.padEnd(3) + F(a.p[0]) + F(a.p[1]) + F(a.p[2]);

  function intOpt(v, def, what) {
    if (v === undefined || v === null || v === '') return def;
    const n = Number(v);
    if (!Number.isInteger(n)) throw new Error(what + ' must be an integer (got "' + v + '").');
    return n;
  }
  function fileBase(s, def) {
    const t = String(s || '').replace(/\.[A-Za-z0-9]{1,4}$/, '').replace(/[^A-Za-z0-9_.-]+/g, '_').replace(/^[._-]+|[._-]+$/g, '');
    return (t || def).slice(0, 40);
  }

  /** Validated, normalised job options shared by all generators. */
  function job(m, defMethod, defBasis) {
    check(m);
    const raw = m.task === undefined || m.task === null || m.task === '' ? 'sp' : String(m.task);
    const task = TASK_ALIAS[raw.toLowerCase().replace(/[^a-z]/g, '')];
    if (!task) throw new Error('Unknown task "' + raw + '": use sp, opt, freq or optfreq.');
    const charge = intOpt(m.charge, 0, 'The charge');
    const mult = intOpt(m.multiplicity, 1, 'The multiplicity');
    if (mult < 1) throw new Error('The multiplicity must be 1 or larger.');
    const atoms = [];
    let dummies = 0, electrons = -charge;
    for (let i = 0; i < m.Z.length; i++) {
      const z = m.Z[i];
      if (!(z > 0)) { dummies++; continue; }
      const p = [+m.xyz[3 * i], +m.xyz[3 * i + 1], +m.xyz[3 * i + 2]];
      if (!p.every(Number.isFinite)) throw new Error('Atom ' + (i + 1) + ' has a non-numeric coordinate.');
      electrons += z;
      atoms.push({ Z: z, s: BM.el.symbol(z), p });
    }
    if (!atoms.length) throw new Error('There are no real atoms to write.');
    const unpaired = mult - 1;
    return {
      atoms, dummies, charge, mult, task, electrons,
      opt: task === 'opt' || task === 'optfreq', freq: task === 'freq' || task === 'optfreq',
      method: oneLine(m.method) || defMethod, basis: oneLine(m.basis) || defBasis,
      nprocs: Math.max(1, Math.round(Number(m.nprocs) || 1)),
      memoryMB: Number(m.memoryMB) > 0 ? Math.round(Number(m.memoryMB)) : 0,
      title: oneLine(m.title), solvent: oneLine(m.solvent),
      extra: m.extra === undefined || m.extra === null ? '' : String(m.extra).replace(/\r\n?/g, '\n').replace(/\s+$/, ''),
      parityOK: electrons >= unpaired && (((electrons - unpaired) % 2) + 2) % 2 === 0
    };
  }

  // ------------------------------------------------------------------ Psi4
  const PCM_SOLVENTS = { water: 'Water', h2o: 'Water', methanol: 'Methanol', ethanol: 'Ethanol', acetonitrile: 'Acetonitrile',
    mecn: 'Acetonitrile', ch3cn: 'Acetonitrile', dmso: 'Dimethylsulfoxide', dimethylsulfoxide: 'Dimethylsulfoxide',
    thf: 'Tetrahydrofurane', tetrahydrofuran: 'Tetrahydrofurane', dichloromethane: 'Methylenechloride', dcm: 'Methylenechloride',
    ch2cl2: 'Methylenechloride', chloroform: 'Chloroform', chcl3: 'Chloroform', benzene: 'Benzene', toluene: 'Toluene',
    acetone: 'Acetone', cyclohexane: 'Cyclohexane', heptane: 'N-heptane', nheptane: 'N-heptane', aniline: 'Aniline',
    nitromethane: 'Nitromethane', chlorobenzene: 'Chlorobenzene', ccl4: 'Carbon Tetrachloride',
    carbontetrachloride: 'Carbon Tetrachloride', '12dichloroethane': '1,2-Dichloroethane', dichloroethane: '1,2-Dichloroethane',
    dioxane: '1,4-Dioxane', '14dioxane': '1,4-Dioxane', propylenecarbonate: 'Propylene Carbonate' };

  function psi4(m) {
    const j = job(m, 'B3LYP', 'def2-SVP');
    const name = methodHead(j.method).toLowerCase();
    const L = ['# ' + (j.title || 'BlueMolecule job')];
    if (j.memoryMB) L.push('memory ' + j.memoryMB + ' MB');
    if (j.nprocs > 1) L.push('set_num_threads(' + j.nprocs + ')');
    L.push('', 'molecule {', j.charge + ' ' + j.mult);
    for (const a of j.atoms) L.push(jobAtomLine(a));
    L.push('units angstrom', '}', '');
    if (needsBasis(name)) L.push('set basis ' + j.basis);
    if (j.mult > 1) L.push('set reference ' + (WAVEFUNCTION.test(name) ? 'uhf' : 'uks'));
    if (j.solvent) {
      L.push('set pcm true', 'set pcm_scf_type total', '', 'pcm = {', '  Units = Angstrom', '  Medium {', '    SolverType = IEFPCM',
        '    Solvent = ' + (PCM_SOLVENTS[solventKey(j.solvent)] || j.solvent), '  }', '  Cavity {', '    RadiiSet = UFF',
        '    Type = GePol', '    Scaling = False', '    Area = 0.3', '    Mode = Implicit', '  }', '}');
    }
    if (j.extra) L.push('', j.extra);
    L.push('');
    const call = "('" + name.replace(/'/g, '') + "')";
    if (j.task === 'sp') L.push('energy' + call);
    if (j.opt) L.push('optimize' + call);
    if (j.freq) L.push('frequency' + call);
    return L.join('\n') + '\n';
  }

  // ------------------------------------------------------------------ Q-Chem
  function qchem(m) {
    const j = job(m, 'B3LYP', 'def2-SVP');
    const sections = /^\s*\$/.test(j.extra) ? j.extra : '';
    const remExtra = j.extra && !sections ? j.extra.split('\n').map(s => s.trim()).filter(Boolean) : [];
    const rem = (jobtype, readGuess) => {
      const R = [['JOBTYPE', jobtype], ['METHOD', methodHead(j.method).toLowerCase()]];
      if (needsBasis(j.method)) R.push(['BASIS', j.basis]);
      if (j.memoryMB) R.push(['MEM_TOTAL', j.memoryMB]);
      if (j.solvent) R.push(['SOLVENT_METHOD', 'PCM']);
      if (readGuess) R.push(['SCF_GUESS', 'read']);
      return ['$rem'].concat(R.map(([k, v]) => '   ' + k.padEnd(16) + v), remExtra.map(s => '   ' + s), ['$end']);
    };
    const tail = [];
    if (j.solvent) tail.push('', '$solvent', '   SolventName      ' + j.solvent, '$end');
    if (sections) tail.push('', sections);
    const L = ['$comment', j.title || 'BlueMolecule job', '$end', '', '$molecule', j.charge + ' ' + j.mult];
    for (const a of j.atoms) L.push('   ' + jobAtomLine(a));
    L.push('$end', '', ...rem(j.task === 'optfreq' ? 'opt' : j.task, false), ...tail);
    if (j.task === 'optfreq') L.push('', '@@@', '', '$molecule', '   read', '$end', '', ...rem('freq', true), ...tail);
    return L.join('\n') + '\n';
  }

  // ------------------------------------------------------------------ NWChem
  const NW_XC = { b3lyp: 'b3lyp', pbe0: 'pbe0', pbe1pbe: 'pbe0', pbe: 'xpbe96 cpbe96', blyp: 'becke88 lyp',
    bp86: 'becke88 perdew86', tpss: 'xtpss03 ctpss03', 'm06-2x': 'm06-2x', m062x: 'm06-2x', m06: 'm06', 'm06-l': 'm06-l',
    m06l: 'm06-l', svwn: 'slater vwn_5', lda: 'slater vwn_5',
    'cam-b3lyp': 'xcamb88 1.00 lyp 0.81 vwn_5 0.19 hfexch 1.00\n  cam 0.33 cam_alpha 0.19 cam_beta 0.46' };
  /** NWChem basis-library name: lower case, Pople polarisation in star notation. */
  function nwBasis(b) {
    let s = b.toLowerCase();
    const m = /^(\d-\d+\+{0,2}g)\((d|d,p)\)$/.exec(s);
    if (m) s = m[1] + (m[2] === 'd' ? '*' : '**');
    return s;
  }

  function nwchem(m) {
    const j = job(m, 'B3LYP', 'def2-SVP');
    const lo = methodHead(j.method).toLowerCase();
    const L = ['echo', 'start ' + fileBase(m.baseName || j.title, 'job'),
      'title "' + (j.title || 'BlueMolecule job').replace(/"/g, "'") + '"'];
    if (j.memoryMB) L.push('memory total ' + Math.max(1, Math.floor(j.memoryMB / j.nprocs)) + ' mb');
    L.push('charge ' + j.charge, '', 'geometry units angstrom');
    for (const a of j.atoms) L.push('  ' + jobAtomLine(a));
    L.push('end', '', 'basis', '  * library ' + nwBasis(j.basis), 'end', '');
    const scf = () => ['scf'].concat(j.mult > 1 ? ['  uhf', '  nopen ' + (j.mult - 1)] : ['  rhf'], ['end', '']);
    let theory;
    if (HF_RE.test(lo)) { theory = 'scf'; L.push(...scf()); }
    else if (MP2_RE.test(lo)) { theory = 'mp2'; L.push(...scf()); }
    else if (CC_RE.test(lo)) { theory = lo; L.push(...scf()); }
    else {
      theory = 'dft';
      let xc = lo, disp = 0, d;
      if ((d = /^(.+?)-?d3\(?bj\)?$/.exec(xc))) { xc = d[1]; disp = 4; }
      else if ((d = /^(.+?)-d3$/.exec(xc))) { xc = d[1]; disp = 3; }
      else if ((d = /^(.+?)-d2$/.exec(xc))) { xc = d[1]; disp = 2; }
      L.push('dft', '  xc ' + (NW_XC[xc] || xc));
      if (disp) L.push('  disp vdw ' + disp);
      L.push('  mult ' + j.mult);
      if (j.mult > 1) L.push('  odft');
      L.push('end', '');
    }
    if (j.solvent) L.push('cosmo', '  solvent ' + j.solvent.toLowerCase().replace(/\s+/g, ''), 'end', '');
    if (j.extra) L.push(j.extra, '');
    if (j.task === 'sp') L.push('task ' + theory + ' energy');
    if (j.opt) L.push('task ' + theory + ' optimize');
    if (j.freq) L.push('task ' + theory + ' freq');
    return L.join('\n') + '\n';
  }

  // ------------------------------------------------------------------ GAMESS (US)
  const GAMESS_SOLVENTS = { water: 'WATER', h2o: 'WATER', methanol: 'CH3OH', ethanol: 'C2H5OH', chloroform: 'CLFORM',
    chcl3: 'CLFORM', dichloromethane: 'METHYCL', dcm: 'METHYCL', ch2cl2: 'METHYCL', carbontetrachloride: 'CTCL', ccl4: 'CTCL',
    '12dichloroethane': '12DCLET', dichloroethane: '12DCLET', benzene: 'BENZENE', toluene: 'TOLUENE', chlorobenzene: 'CLBENZ',
    nitromethane: 'NITMET', heptane: 'NEPTANE', nheptane: 'NEPTANE', cyclohexane: 'CYCHEX', aniline: 'ANILINE',
    acetone: 'ACETONE', thf: 'THF', tetrahydrofuran: 'THF', dmso: 'DMSO', dimethylsulfoxide: 'DMSO', acetonitrile: 'CH3CN',
    ch3cn: 'CH3CN' };
  const GAMESS_DFT = { 'CAM-B3LYP': 'CAMB3LYP', PBE1PBE: 'PBE0', 'WB97X-D': 'wB97X-D', 'WB97X': 'wB97X' };

  /** $BASIS keywords for a basis-set name: {tokens, note}. */
  function gamessBasis(b) {
    const s = b.replace(/\s+/g, '').toUpperCase();
    let m;
    if ((m = /^STO-?(\d)G$/.exec(s))) return { tokens: ['GBASIS=STO', 'NGAUSS=' + m[1]] };
    if ((m = /^(3-21|6-31|6-311)(\+{0,2})G(\*{0,2}|\(([^)]*)\))$/.exec(s))) {
      const t = ['GBASIS=' + { '3-21': 'N21', '6-31': 'N31', '6-311': 'N311' }[m[1]], 'NGAUSS=' + (m[1] === '3-21' ? 3 : 6)];
      let heavy = '', light = '';
      if (m[3] === '*') heavy = 'D';
      else if (m[3] === '**') { heavy = 'D'; light = 'P'; }
      else if (m[4]) { const parts = m[4].split(','); heavy = parts[0] || ''; light = parts[1] || ''; }
      const count = (str, letter) => { const q = new RegExp('(\\d*)' + letter).exec(str); return q ? (q[1] ? +q[1] : 1) : 0; };
      const nd = count(heavy, 'D'), nf = count(heavy, 'F'), np = count(light, 'P');
      if (nd) t.push('NDFUNC=' + nd);
      if (nf) t.push('NFFUNC=' + nf);
      if (np) t.push('NPFUNC=' + np);
      if (m[2].length >= 1) t.push('DIFFSP=.TRUE.');
      if (m[2].length === 2) t.push('DIFFS=.TRUE.');
      return { tokens: t, note: count(light, 'D') ? 'GAMESS $BASIS has no keyword for d functions on hydrogen; they were left out.' : null };
    }
    if ((m = /^(AUG-)?CC-P(C?)V([DTQ56])Z$/.exec(s))) return { tokens: ['GBASIS=' + (m[1] ? 'A' : '') + 'CC' + m[3] + (m[2] ? 'C' : '')] };
    if ((m = /^(AUG-)?PCSEG-([0-4])$/.exec(s))) return { tokens: ['GBASIS=' + (m[1] ? 'APCSEG-' : 'PCSEG-') + m[2]] };
    return { tokens: ['GBASIS=' + s], note: 'Check that GBASIS=' + s + ' is a basis set your GAMESS version knows.' };
  }
  /** ' $NAME A=1 B=2 $END', continued on further lines so that no line passes column 79. */
  function gamessGroup(name, items) {
    const out = [];
    let line = ' $' + name;
    for (const it of items) {
      if (line.length + it.length + 1 > 74) { out.push(line); line = '  '; }
      line += ' ' + it;
    }
    if (line.length + 5 > 79) { out.push(line); line = ' '; }
    out.push(line + ' $END');
    return out;
  }
  function gamessParts(m) {
    const j = job(m, 'B3LYP', '6-31G(d)');
    const lo = methodHead(j.method).toLowerCase();
    const semi = /^(am1|pm3|mndo|rm1)$/.test(lo);
    const ctrl = ['SCFTYP=' + (j.mult > 1 ? 'UHF' : 'RHF'), 'RUNTYP=' + (j.opt ? 'OPTIMIZE' : j.freq ? 'HESSIAN' : 'ENERGY'),
      'ICHARG=' + j.charge, 'MULT=' + j.mult];
    let basis, note = null;
    if (semi) basis = ['GBASIS=' + lo.toUpperCase()];
    else {
      const b = gamessBasis(j.basis);
      basis = b.tokens; note = b.note;
      if (MP2_RE.test(lo)) ctrl.push('MPLEVL=2');
      else if (CC_RE.test(lo)) ctrl.push('CCTYP=' + lo.toUpperCase());
      else if (!HF_RE.test(lo)) { const up = methodHead(j.method).toUpperCase(); ctrl.push('DFTTYP=' + (GAMESS_DFT[up] || up)); }
    }
    ctrl.push('MAXIT=100');
    const L = gamessGroup('CONTRL', ctrl);
    if (j.memoryMB) L.push(...gamessGroup('SYSTEM', ['MWORDS=' + Math.max(1, Math.floor(j.memoryMB * 1048576 / 8e6 / j.nprocs))]));
    L.push(...gamessGroup('BASIS', basis));
    if (j.opt) L.push(...gamessGroup('STATPT', ['NSTEP=100'].concat(j.freq ? ['HSSEND=.TRUE.'] : [])));
    if (j.freq && !semi && !(HF_RE.test(lo) && j.mult === 1)) L.push(...gamessGroup('FORCE', ['METHOD=SEMINUM']));
    if (j.solvent) L.push(...gamessGroup('PCM', ['SOLVNT=' + (GAMESS_SOLVENTS[solventKey(j.solvent)] || j.solvent.toUpperCase())]));
    if (j.extra) for (const s of j.extra.split('\n')) L.push(/^\s/.test(s) ? s : ' ' + s);   // groups start in column 2
    L.push(' $DATA', (j.title || 'BlueMolecule job').slice(0, 80), 'C1');
    for (const a of j.atoms) L.push(a.s.padEnd(3) + a.Z.toFixed(1).padStart(6) + F(a.p[0]) + F(a.p[1]) + F(a.p[2]));
    L.push(' $END');
    return { text: L.join('\n') + '\n', note, j };
  }
  const gamess = m => gamessParts(m).text;

  // ------------------------------------------------------------------ MOPAC
  const DIELECTRIC = { water: 78.36, h2o: 78.36, methanol: 32.61, ethanol: 24.85, acetonitrile: 35.69, dmso: 46.83,
    dimethylsulfoxide: 46.83, dmf: 37.22, dimethylformamide: 37.22, nndimethylformamide: 37.22, acetone: 20.49, thf: 7.43,
    tetrahydrofuran: 7.43, dichloromethane: 8.93, dcm: 8.93, ch2cl2: 8.93, chloroform: 4.71, chcl3: 4.71, toluene: 2.37,
    benzene: 2.27, hexane: 1.88, nhexane: 1.88, heptane: 1.91, nheptane: 1.91, cyclohexane: 2.02, diethylether: 4.24,
    ether: 4.24, ccl4: 2.23, carbontetrachloride: 2.23, nitromethane: 36.56, ethylacetate: 5.99, pyridine: 12.98,
    chlorobenzene: 5.70, dichloroethane: 10.13, '12dichloroethane': 10.13, dioxane: 2.21, '14dioxane': 2.21 };

  function mopac(m) {
    const j = job(m, 'PM7', '');
    if (j.mult > 9) throw new Error('MOPAC multiplicity keywords go up to NONET (9).');
    const kw = [j.method.toUpperCase(), 'CHARGE=' + j.charge, MULT_NAMES[j.mult]];
    if (j.task === 'sp') kw.push('1SCF');
    if (j.task === 'freq') kw.push('FORCE');
    if (j.task === 'optfreq') kw.push('OPT', 'FORCE');
    if (j.nprocs > 1) kw.push('THREADS=' + j.nprocs);
    if (j.solvent) {
      const eps = /^\d+(\.\d*)?$/.test(j.solvent) ? +j.solvent : DIELECTRIC[solventKey(j.solvent)];
      if (!eps) throw new Error('MOPAC needs the dielectric constant of "' + j.solvent + '" (e.g. solvent "78.4").');
      kw.push('EPS=' + eps);
    }
    if (j.extra) kw.push(j.extra.replace(/\s+/g, ' ').trim());
    // One keyword line, or two joined by "+" (MOPAC then still reads two title lines).
    const L = [];
    let line = '';
    for (const k of kw) {
      if (line && line.length + k.length + 1 > 78 && L.length === 0) { L.push(line + ' +'); line = ''; }
      line += (line ? ' ' : '') + k;
    }
    L.push(line, j.title || 'BlueMolecule job', 'BlueMolecule: ' + { sp: 'single point', opt: 'geometry optimisation',
      freq: 'vibrational frequencies', optfreq: 'optimisation + frequencies' }[j.task]);
    const flag = j.task === 'sp' ? ' 0' : ' 1';
    for (const a of j.atoms) L.push(a.s.padEnd(3) + F(a.p[0]) + flag + F(a.p[1]) + flag + F(a.p[2]) + flag);
    return L.join('\n') + '\n';
  }

  // ------------------------------------------------------------------ xtb
  const ALPB = { water: 'water', h2o: 'water', acetone: 'acetone', acetonitrile: 'acetonitrile', mecn: 'acetonitrile',
    aniline: 'aniline', benzaldehyde: 'benzaldehyde', benzene: 'benzene', dichloromethane: 'ch2cl2', dcm: 'ch2cl2',
    ch2cl2: 'ch2cl2', chloroform: 'chcl3', chcl3: 'chcl3', carbondisulfide: 'cs2', cs2: 'cs2', dioxane: 'dioxane',
    '14dioxane': 'dioxane', dmf: 'dmf', dimethylformamide: 'dmf', dmso: 'dmso', dimethylsulfoxide: 'dmso', ether: 'ether',
    diethylether: 'ether', ethylacetate: 'ethylacetate', furan: 'furane', furane: 'furane', hexadecane: 'hexandecane',
    hexandecane: 'hexandecane', hexane: 'hexane', methanol: 'methanol', nitromethane: 'nitromethane', octanol: 'octanol',
    woctanol: 'woctanol', phenol: 'phenol', toluene: 'toluene', thf: 'thf', tetrahydrofuran: 'thf' };

  function xtb(m) {
    const j = job(m, 'GFN2-xTB', '');
    const lo = methodHead(j.method).toLowerCase().replace(/[^a-z0-9]/g, '');
    let gfn;
    if (lo === 'gfnff') gfn = ['--gfnff'];
    else {
      const g = /^gfn([012])(xtb)?$/.exec(lo) || /^xtb([012])?$/.exec(lo);
      if (!g) throw new Error('xtb runs GFN2-xTB, GFN1-xTB, GFN0-xTB or GFN-FF, not "' + j.method + '".');
      gfn = ['--gfn', g[1] || '2'];
    }
    const base = fileBase(m.baseName, 'mol');
    const X = [String(j.atoms.length), j.title || 'BlueMolecule'];
    for (const a of j.atoms) X.push(jobAtomLine(a));
    const files = [{ name: base + '.xyz', text: X.join('\n') + '\n' }];
    const cmd = ['xtb', base + '.xyz'];
    if (j.task !== 'sp') cmd.push({ opt: '--opt', freq: '--hess', optfreq: '--ohess' }[j.task]);
    cmd.push('--chrg', String(j.charge), '--uhf', String(j.mult - 1), ...gfn);
    if (j.solvent) cmd.push('--alpb', ALPB[solventKey(j.solvent)] || j.solvent.toLowerCase().replace(/\s+/g, ''));
    if (j.nprocs > 1) cmd.push('-P', String(j.nprocs));
    if (j.extra) {
      if (/^\s*\$/.test(j.extra)) {
        files.push({ name: 'xcontrol', text: j.extra + (/\$end\s*$/i.test(j.extra) ? '' : '\n$end') + '\n' });
        cmd.push('--input', 'xcontrol');
      } else cmd.push(j.extra.replace(/\s+/g, ' ').trim());
    }
    return { files, command: cmd.join(' ') };
  }

  // ------------------------------------------------------------------ Molpro
  function molpro(m) {
    const j = job(m, 'B3LYP', 'cc-pVDZ');
    const lo = methodHead(j.method).toLowerCase();
    const safe = s => String(s).replace(/[,;!{}]/g, ' ').replace(/\s+/g, ' ').trim();
    const title = safe(j.title) || 'BlueMolecule job';
    const L = ['***,' + title];
    if (j.memoryMB) L.push('memory,' + Math.max(1, Math.floor(j.memoryMB / j.nprocs / 8)) + ',m   ! megawords per process');
    L.push('geomtyp=xyz', 'geometry={', String(j.atoms.length), title);
    for (const a of j.atoms) L.push(jobAtomLine(a));
    L.push('}');
    if (needsBasis(j.method)) L.push('basis=' + j.basis);
    L.push('set,charge=' + j.charge, 'set,spin=' + (j.mult - 1));
    if (j.solvent) L.push('! implicit solvent "' + safe(j.solvent) + '" requested: add the COSMO/PCM directive for your Molpro version');
    if (j.extra) L.push(j.extra);
    const open = j.mult > 1;
    if (HF_RE.test(lo)) L.push(open ? '{uhf}' : '{hf}');
    else if (MP2_RE.test(lo)) L.push(open ? '{rhf}' : '{hf}', open ? '{rmp2}' : '{mp2}');
    else if (CC_RE.test(lo)) L.push(open ? '{rhf}' : '{hf}', '{' + (open ? 'u' : '') + lo + '}');
    else L.push('{' + (open ? 'uks' : 'ks') + ',' + lo + '}');
    if (j.opt) L.push('optg');
    if (j.freq) L.push('{frequencies}');
    L.push('---');
    return L.join('\n') + '\n';
  }

  // ------------------------------------------------------------------ Turbomole
  function turbomole(m) {
    const j = job(m, '', '');
    const L = ['$coord'];
    for (const a of j.atoms) {
      L.push(fx(a.p[0] / BOHR, 22, 14) + fx(a.p[1] / BOHR, 22, 14) + fx(a.p[2] / BOHR, 22, 14) + '      ' + a.s.toLowerCase());
    }
    if (j.extra) L.push(j.extra);
    L.push('$end');
    return L.join('\n') + '\n';
  }

  // ------------------------------------------------------------------ program registry
  function gaussBasis(b) { return /^def2-/i.test(b) ? 'Def2' + b.slice(5) : b; }
  function notes(j, list) {
    const N = [];
    if (!j.parityOK) N.push('Check charge and multiplicity: ' + j.electrons + ' electrons cannot have multiplicity ' + j.mult + '.');
    for (const s of list) if (s) N.push(s);
    return N.length ? N.join(' ') : undefined;
  }
  const dummyNote = j => (j.dummies ? j.dummies + ' dummy/ghost atom(s) were left out.' : null);
  const out = (text, j, command, list) => {
    const r = { text, command };
    const n = notes(j, list || []);
    if (n) r.note = n;
    return r;
  };

  function writeGaussian(o) {
    const j = job(o, 'B3LYP', '6-31G(d)');
    const parts = j.method.split(/\s+/), head = parts[0];
    const route = ['#p', head.indexOf('/') < 0 && needsBasis(head) ? head + '/' + gaussBasis(j.basis) : head].concat(parts.slice(1));
    if (j.opt) route.push('opt');
    if (j.freq) route.push('freq');
    if (j.solvent) route.push('scrf=(smd,solvent=' + j.solvent.replace(/\s+/g, '') + ')');
    const link0 = [];
    if (j.nprocs > 1) link0.push('%nprocshared=' + j.nprocs);
    if (j.memoryMB) link0.push('%mem=' + j.memoryMB + 'MB');
    let text = gaussianInput({ Z: o.Z, xyz: o.xyz, route: route.join(' '), title: j.title || 'BlueMolecule job',
      charge: j.charge, multiplicity: j.mult, link0 });
    if (j.extra) text = text.replace(/\n\n\n$/, '\n\n' + j.extra + '\n\n');
    return out(text, j, 'g16 ' + fileBase(o.baseName, 'job') + '.gjf',
      [j.dummies ? j.dummies + ' dummy atom(s) written as Bq ghost atoms.' : null]);
  }

  function writeOrca(o) {
    const j = job(o, 'B3LYP', 'def2-SVP');
    const kw = [j.method];
    if (needsBasis(j.method)) kw.push(j.basis);
    if (j.opt) kw.push('Opt');
    if (j.freq) kw.push('Freq');
    if (j.solvent) kw.push('CPCM(' + j.solvent.replace(/\s+/g, '') + ')');
    // The title goes in a comment after the "!" lines: a leading "#" line would look like a Gaussian route.
    const text = orcaInput({ Z: o.Z, xyz: o.xyz, keywords: kw.join(' '), charge: j.charge, multiplicity: j.mult,
      nprocs: j.nprocs, maxcore: j.memoryMB ? Math.max(1, Math.floor(j.memoryMB / j.nprocs)) : 0,
      extraBlocks: j.extra ? [j.extra] : [] })
      .replace(/^((?:![^\n]*\n)+)/, '$1# ' + (j.title || 'BlueMolecule job').replace(/#/g, '') + '\n');
    const base = fileBase(o.baseName, 'job');
    return out(text, j, (j.nprocs > 1 ? '$(which orca) ' : 'orca ') + base + '.inp > ' + base + '.out',
      [j.dummies ? j.dummies + ' dummy atom(s) written as DA.' : null]);
  }

  function program(id, name, ext, defaultMethod, defaultBasis, methods, bases, write) {
    return { id, name, ext, defaultMethod, defaultBasis, methods, bases, tasks: TASKS.slice(), write };
  }
  const withDefaults = (o, p) => Object.assign({}, o, {
    method: oneLine(o && o.method) || p.defaultMethod || undefined,
    basis: oneLine(o && o.basis) || p.defaultBasis || undefined
  });

  const PROGRAMS = [
    program('gaussian', 'Gaussian', 'gjf', 'B3LYP', '6-31G(d)',
      ['B3LYP', 'PBE1PBE', 'M062X', 'wB97XD', 'CAM-B3LYP', 'B2PLYPD3', 'HF', 'MP2', 'CCSD(T)', 'PM7', 'AM1'],
      ['6-31G(d)', '6-31+G(d,p)', '6-311+G(d,p)', 'def2SVP', 'def2TZVP', 'cc-pVDZ', 'cc-pVTZ', 'aug-cc-pVTZ', 'LANL2DZ', 'SDD', 'Gen', 'GenECP'],
      writeGaussian),
    program('orca', 'ORCA', 'inp', 'B3LYP', 'def2-SVP',
      ['B3LYP', 'PBE0', 'r2SCAN-3c', 'B97-3c', 'wB97X-D3', 'M06-2X', 'TPSSh', 'HF', 'RI-MP2', 'DLPNO-CCSD(T)', 'XTB2'],
      ['def2-SVP', 'def2-TZVP', 'def2-TZVPP', 'def2-QZVP', 'ma-def2-TZVP', 'cc-pVDZ', 'cc-pVTZ', 'aug-cc-pVTZ', '6-31G(d)'],
      writeOrca),
    program('psi4', 'Psi4', 'in', 'B3LYP', 'def2-SVP',
      ['B3LYP', 'B3LYP-D3BJ', 'PBE0', 'wB97X-D', 'M06-2X', 'HF', 'MP2', 'CCSD(T)'],
      ['def2-SVP', 'def2-TZVP', 'cc-pVDZ', 'cc-pVTZ', 'aug-cc-pVDZ', 'aug-cc-pVTZ', 'jun-cc-pVDZ', '6-31G*'],
      o => {
        const j = job(o, 'B3LYP', 'def2-SVP'), b = fileBase(o.baseName, 'job');
        return out(psi4(o), j, 'psi4 ' + (j.nprocs > 1 ? '-n ' + j.nprocs + ' ' : '') + b + '.in ' + b + '.out',
          [dummyNote(j), j.solvent ? 'PCM needs Psi4 built with PCMSolver.' : null]);
      }),
    program('qchem', 'Q-Chem', 'in', 'B3LYP', 'def2-SVP',
      ['B3LYP', 'wB97X-D', 'wB97M-V', 'PBE0', 'M06-2X', 'HF', 'MP2', 'CCSD(T)'],
      ['def2-SVP', 'def2-TZVP', '6-31G*', '6-311+G**', 'cc-pVDZ', 'cc-pVTZ', 'aug-cc-pVTZ'],
      o => {
        const j = job(o, 'B3LYP', 'def2-SVP'), b = fileBase(o.baseName, 'job');
        return out(qchem(o), j, 'qchem ' + (j.nprocs > 1 ? '-nt ' + j.nprocs + ' ' : '') + b + '.in ' + b + '.out', [dummyNote(j)]);
      }),
    program('nwchem', 'NWChem', 'nw', 'B3LYP', 'def2-SVP',
      ['B3LYP', 'PBE0', 'PBE', 'M06-2X', 'CAM-B3LYP', 'TPSS', 'B3LYP-D3BJ', 'HF', 'MP2', 'CCSD(T)'],
      ['def2-SVP', 'def2-TZVP', '6-31G*', '6-311+G**', 'cc-pVDZ', 'cc-pVTZ', 'aug-cc-pVDZ'],
      o => {
        const j = job(o, 'B3LYP', 'def2-SVP'), b = fileBase(o.baseName, 'job');
        return out(nwchem(o), j, (j.nprocs > 1 ? 'mpirun -np ' + j.nprocs + ' ' : '') + 'nwchem ' + b + '.nw > ' + b + '.out',
          [dummyNote(j)]);
      }),
    program('gamess', 'GAMESS (US)', 'inp', 'B3LYP', '6-31G(d)',
      ['B3LYP', 'PBE0', 'M06-2X', 'wB97X-D', 'HF', 'MP2', 'CCSD(T)', 'AM1', 'PM3'],
      ['6-31G(d)', '6-31G(d,p)', '6-311G(d,p)', '6-31+G(d)', 'STO-3G', 'cc-pVDZ', 'cc-pVTZ', 'aug-cc-pVDZ', 'PCseg-1'],
      o => {
        const g = gamessParts(o), b = fileBase(o.baseName, 'job');
        return out(g.text, g.j, 'rungms ' + b + '.inp 00 ' + g.j.nprocs + ' > ' + b + '.log', [dummyNote(g.j), g.note]);
      }),
    program('mopac', 'MOPAC', 'mop', 'PM7', null,
      ['PM7', 'PM6', 'PM6-D3H4', 'PM6-ORG', 'RM1', 'AM1', 'PM3', 'MNDO'], [],
      o => {
        const j = job(o, 'PM7', ''), b = fileBase(o.baseName, 'job');
        return out(mopac(o), j, 'mopac ' + b + '.mop', [dummyNote(j),
          j.task === 'optfreq' ? 'FORCE results are only meaningful if the optimisation converged (small gradient norm).' : null]);
      }),
    program('xtb', 'xtb', 'xyz', 'GFN2-xTB', null, ['GFN2-xTB', 'GFN1-xTB', 'GFN0-xTB', 'GFN-FF'], [],
      o => {
        const j = job(o, 'GFN2-xTB', ''), r = xtb(o);
        const res = out(r.files[0].text, j, r.command, [dummyNote(j)]);
        res.files = r.files;
        return res;
      }),
    program('molpro', 'Molpro', 'inp', 'B3LYP', 'cc-pVDZ',
      ['B3LYP', 'PBE0', 'PBE', 'HF', 'MP2', 'CCSD(T)'],
      ['cc-pVDZ', 'cc-pVTZ', 'aug-cc-pVDZ', 'aug-cc-pVTZ', 'def2-SVP', 'def2-TZVP'],
      o => {
        const j = job(o, 'B3LYP', 'cc-pVDZ'), b = fileBase(o.baseName, 'job');
        return out(molpro(o), j, 'molpro ' + (j.nprocs > 1 ? '-n ' + j.nprocs + ' ' : '') + b + '.inp', [dummyNote(j),
          j.solvent ? 'Implicit solvent is not written for Molpro (see the comment in the input).' : null]);
      }),
    program('turbomole', 'Turbomole (coord)', '', null, null, [], [],
      o => {
        const j = job(o, '', ''), text = turbomole(o);
        const run = { sp: 'ridft (or dscf)', opt: 'jobex -ri', freq: 'aoforce', optfreq: 'jobex -ri, then aoforce' }[j.task];
        const res = out(text, j, undefined, [dummyNote(j), 'Only the coord file (bohr) is written: set method, basis, ' +
          'charge and multiplicity with define' + (j.solvent ? ' and the solvent with cosmoprep' : '') + ', then run ' + run + '.']);
        res.files = [{ name: 'coord', text }];
        return res;
      })
  ];
  // Methods and bases left empty take each program's defaults.
  for (const p of PROGRAMS) { const w = p.write; p.write = o => w(withDefaults(o || {}, p)); }

  BM.writers = { xyz, xyzMulti, gaussianInput, orcaInput, mol, pdb, psi4, qchem, nwchem, gamess, mopac, xtb, molpro, turbomole };
  BM.inputPrograms = PROGRAMS;

  // ------------------------------------------------------------------ "Save as" registry
  // BM.registerSaveFormat({ id, name, ext, kind: 'structure'|'trajectory'|'volume'|'vibration'|'orbitals',
  //   needs?: ['cell'|'volume'|'vibrations'|'frames'|'mos'], binary?: false, write(ctx) -> string | {files:[{name,text}]} })
  // ctx = { Z, xyz, bonds: [[i,j,order]] | null, formalCharges, cell, title, atomNames, resNames, resIds, chains,
  //         frames: [{Z, xyz, energy, label}], frameIndex, volume, vibrations, basis, mos, charges, meta, baseName }
  const SAVE = [];
  function registerSaveFormat(f) {
    if (!f || !f.id || !f.name || !f.ext || typeof f.write !== 'function') throw new Error('registerSaveFormat: id, name, ext and write() are required');
    const e = Object.assign({ kind: 'structure', needs: [] }, f);
    const k = SAVE.findIndex(x => x.id === e.id);
    if (k >= 0) SAVE[k] = e; else SAVE.push(e);
    return e;
  }
  const pairsToBonds = c => (c.bonds && c.bonds.length ? c.bonds : null);
  registerSaveFormat({ id: 'xyz', name: 'XYZ', ext: 'xyz', write: c => xyz({ Z: c.Z, xyz: c.xyz, comment: c.title || '' }) });
  registerSaveFormat({ id: 'xyz-traj', name: 'XYZ trajectory (all geometries)', ext: 'xyz', kind: 'trajectory', needs: ['frames'],
    write: c => xyzMulti((c.frames || []).map((f, i) => ({ Z: f.Z, xyz: f.xyz, comment: (c.title || '') + ' geometry ' + (i + 1) + (f.energy !== null && f.energy !== undefined && isFinite(f.energy) ? ' E = ' + f.energy.toFixed(10) : '') }))) });
  registerSaveFormat({ id: 'mol', name: 'MDL Molfile', ext: 'mol', write: c => mol({ Z: c.Z, xyz: c.xyz, bonds: pairsToBonds(c), title: c.title || '' }) });
  registerSaveFormat({ id: 'pdb', name: 'Protein Data Bank', ext: 'pdb', write: c => pdb({ Z: c.Z, xyz: c.xyz, atomNames: c.atomNames, resNames: c.resNames, resIds: c.resIds, cell: c.cell }) });
  registerSaveFormat({ id: 'turbomole-coord', name: 'Turbomole coord (Bohr)', ext: 'coord', write: c => turbomole({ Z: c.Z, xyz: c.xyz }) });
  BM.saveFormats = SAVE;
  BM.registerSaveFormat = registerSaveFormat;
})(typeof window !== 'undefined' ? (window.BM = window.BM || {}) : (globalThis.BM = globalThis.BM || {}));

/* "Save as" formats for quantum-chemistry exchange: Molden (atoms + vibrations), Gaussian cube,
 * Gaussian input and ORCA input built from the current structure and dataset metadata. */
(function (BM) {
  'use strict';
  if (!BM.registerSaveFormat || !BM.writers) return;
  const BOHR = 0.529177210903;
  const sym = z => (z > 0 ? BM.el.symbol(z) : 'X');
  const line1 = s => String(s === undefined || s === null ? '' : s).replace(/[\r\n]+/g, ' ').trim();
  const fx = (x, w, d) => (Number.isFinite(x) ? (Math.abs(x) < 0.5 * Math.pow(10, -d) ? 0 : x) : 0).toFixed(d).padStart(w);
  const ix = (v, w) => String(v).padStart(w);
  /** Fortran E13.5 ("  1.23456E-03"). */
  function e13(x) {
    if (!Number.isFinite(x)) x = 0;
    const s = x.toExponential(5);
    const k = s.indexOf('e');
    let e = s.slice(k + 1);
    const sign = e[0] === '-' ? '-' : '+';
    e = e.replace(/^[-+]/, '');
    return (s.slice(0, k) + 'E' + sign + (e.length < 2 ? '0' + e : e)).padStart(13);
  }

  function molden(c) {
    if (!c || !c.Z || !c.xyz) throw new Error('Molden export needs a structure.');
    const L = ['[Molden Format]', '[Title]', ' ' + (line1(c.title) || 'BlueMolecule'), '[Atoms] Angs'];
    c.Z.forEach((z, i) => L.push(sym(z).padEnd(4) + ix(i + 1, 6) + ix(Math.max(0, z), 5) +
      fx(c.xyz[3 * i], 18, 10) + fx(c.xyz[3 * i + 1], 18, 10) + fx(c.xyz[3 * i + 2], 18, 10)));
    const v = c.vibrations;
    if (v && v.freqs && v.freqs.length && v.modes && v.modes.length) {
      let g = c.frames && c.frames[v.frameIndex];
      if (!g || !g.Z || 3 * g.Z.length !== v.modes[0].length) g = { Z: c.Z, xyz: c.xyz };
      const N = g.Z.length, nm = Math.min(v.freqs.length, v.modes.length);
      if (3 * N !== v.modes[0].length) throw new Error('The normal modes do not match the structure (' + N + ' atoms).');
      L.push('[FREQ]');
      for (let k = 0; k < nm; k++) L.push(fx(v.freqs[k], 14, 4));
      L.push('[FR-COORD]');
      for (let i = 0; i < N; i++) L.push(sym(g.Z[i]).padEnd(4) + fx(g.xyz[3 * i] / BOHR, 18, 10) + fx(g.xyz[3 * i + 1] / BOHR, 18, 10) + fx(g.xyz[3 * i + 2] / BOHR, 18, 10));
      L.push('[FR-NORM-COORD]');
      for (let k = 0; k < nm; k++) {
        L.push(' vibration' + ix(k + 1, 6));
        const m = v.modes[k];
        for (let i = 0; i < N; i++) L.push(fx(m[3 * i], 14, 8) + fx(m[3 * i + 1], 14, 8) + fx(m[3 * i + 2], 14, 8));
      }
      const fin = a => a && a.length >= nm && a.slice(0, nm).some(Number.isFinite);
      if (fin(v.ir)) {
        L.push('[INT]');
        for (let k = 0; k < nm; k++) L.push(fx(v.ir[k], 14, 4) + (fin(v.raman) ? fx(v.raman[k], 14, 4) : ''));
      }
    }
    return L.join('\n') + '\n';
  }

  function cube(c) {
    const v = c && c.volume;
    if (!v || !v.data || !v.dims || !v.axes || !v.origin) throw new Error('There is no volume (cube) data to save.');
    const [nx, ny, nz] = v.dims;
    if (v.data.length < nx * ny * nz) throw new Error('Volume data are shorter than the grid.');
    const Z = c.Z || [], X = c.xyz || [];
    const mo = v.mo !== undefined && v.mo !== null && Number.isFinite(+v.mo);
    const f12 = x => fx(x, 12, 6);
    const L = [line1(v.name || c.title) || 'BlueMolecule volume', line1((v.kind || 'volume') + ' grid written by BlueMolecule, lengths in Bohr')];
    L.push(ix(mo ? -Z.length : Z.length, 5) + f12(v.origin[0] / BOHR) + f12(v.origin[1] / BOHR) + f12(v.origin[2] / BOHR));
    for (let k = 0; k < 3; k++) L.push(ix(v.dims[k], 5) + f12(v.axes[k][0] / BOHR) + f12(v.axes[k][1] / BOHR) + f12(v.axes[k][2] / BOHR));
    Z.forEach((z, i) => L.push(ix(z, 5) + f12(z) + f12(X[3 * i] / BOHR) + f12(X[3 * i + 1] / BOHR) + f12(X[3 * i + 2] / BOHR)));
    if (mo) L.push(ix(1, 5) + ix(+v.mo, 5));
    const d = v.data;
    for (let r = 0; r < nx * ny; r++) {
      const base = r * nz;
      for (let z0 = 0; z0 < nz; z0 += 6) {
        let s = '';
        for (let z = z0; z < Math.min(nz, z0 + 6); z++) s += e13(d[base + z]);
        L.push(s);
      }
    }
    return L.join('\n') + '\n';
  }

  function gaussianRoute(meta) {
    const r = meta && typeof meta.route === 'string' && /^\s*#/.test(meta.route) ? meta.route : '#p';
    // geometry/guess from a checkpoint or a connectivity section would not match the exported file
    return r.replace(/\b(geom|guess)\s*=\s*(\([^)]*\)|\S+)/gi, ' ').replace(/\s+/g, ' ').trim() || '#p';
  }
  const cm = (m, k, d) => (m && Number.isFinite(m[k]) ? m[k] : d);

  const R = BM.registerSaveFormat;
  R({ id: 'molden', name: 'Molden (atoms + vibrations)', ext: 'molden', write: molden });
  R({ id: 'cube', name: 'Gaussian cube (Bohr)', ext: 'cube', kind: 'volume', needs: ['volume'], write: cube });
  R({ id: 'gaussian-input', name: 'Gaussian input', ext: 'gjf', write: c => BM.writers.gaussianInput({ Z: c.Z, xyz: c.xyz,
    route: gaussianRoute(c.meta), title: line1(c.title) || 'BlueMolecule', charge: cm(c.meta, 'charge', 0),
    multiplicity: cm(c.meta, 'multiplicity', 1), cell: c.cell }) });
  R({ id: 'orca-input', name: 'ORCA input', ext: 'inp', write: c => {
    const m = c.meta || {};
    const kw = typeof m.route === 'string' && /^\s*!/.test(m.route) ? m.route.replace(/^\s*!\s*/, '')
      : ([m.method, m.basis].filter(Boolean).join(' ') || 'HF def2-SVP');
    return BM.writers.orcaInput({ Z: c.Z, xyz: c.xyz, keywords: kw, charge: cm(m, 'charge', 0), multiplicity: cm(m, 'multiplicity', 1) });
  } });
  BM.qcWriters = { molden, cube, gaussianRoute };
})(typeof window !== 'undefined' ? (window.BM = window.BM || {}) : (globalThis.BM = globalThis.BM || {}));
