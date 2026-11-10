/* BlueMolecule — atom selection language, residue classification and the representation model.
 *
 *   BM.reps.context({frame, bonds, layers, selection, frozen, charges, ssRecords, lightBg}) -> ctx
 *   BM.reps.parse(text) -> AST (throws BM.reps.SelectionError {message, pos, end})
 *   BM.reps.evaluate(text | AST, ctx) -> Uint8Array mask (1 = selected);  count(mask)
 *   BM.reps.classify(ctx) -> residues + per-atom class (protein, nucleic, water, ions, lipid, carbohydrate, ligand)
 *   BM.reps.sugars(ctx) -> monosaccharide rings: SNFG type, label (GlcNS6S, IdoA2S…), shape and colour
 *   BM.reps.resolve(rows, ctx) -> per-atom style codes, colours, radius scales, opacities, overlays and markers
 *   BM.reps.presetRows(id, ctx, view), autoPreset(ctx), complete(text, caret, ctx), contents(ctx), describe(ctx, i)
 *
 * Representation rows: {name, selection, style, color, opacity, scale, visible}. Later rows override earlier ones
 * for per-atom styles (ballstick, sticks, wireframe, lines, spacefill, ghost, hidden, lipid); overlays (cartoon,
 * trace, rings, snfg, markers) add on top. Atoms that no row styles stay hidden. `color` is a scheme id or '#rrggbb'.
 * Pure data code: runs in the browser and in Node (tests). */
(function (BM) {
  'use strict';
  const el = BM.el;
  const W = s => s.trim().split(/\s+/);

  // ------------------------------------------------------------------ residue tables
  const AA20 = W('ALA ARG ASN ASP CYS GLN GLU GLY HIS ILE LEU LYS MET PHE PRO SER THR TRP TYR VAL');
  const AAVAR = W('HID HIE HIP HSD HSE HSP HSC HISD HISE HISH HISA HISB HIS1 HIS2 CYX CYM CYS2 CYSH ASH ASPH GLH GLUH LYN LYSH LSN ARN TYM');
  const AMINO = new Set([...AA20, ...AAVAR, ...W('SEC PYL MSE HYP SEP TPO PTR CSO CSD CME MLY KCX ALY M3L NLE ORN ABA DAL DLE DVA ACE NME NHE NH2 NMA FOR'),
    ...[...AA20, ...AAVAR].flatMap(r => ['N' + r, 'C' + r])]);
  const NUCLEIC = new Set(W('A C G U T I DA DC DG DT DU DI RA RC RG RU ADE CYT GUA THY URA URI DA5 DA3 DC5 DC3 DG5 DG3 DT5 DT3 RA5 RA3 RC5 RC3 RG5 RG3 RU5 RU3 ' +
    'A5 A3 C5 C3 G5 G3 U5 U3 DAN DCN DGN DTN RAN RCN RGN RUN PSU 5MC OMC OMG 1MA 2MG M2G 7MG H2U 5MU'));
  const WATER = new Set(W('HOH WAT SOL TIP3 TIP3P T3P TP3 TIP4 TIP4P T4P TP4 TIP5 T5P TP5 SPC SPCE H2O DOD OPC OPC3 TIP WTR SWM4 HO4 TIP2'));
  const IONS = new Set(W('NA CL K MG CA ZN FE FE2 FE3 CU CU1 CU2 MN MN3 CO NI CD HG LI RB CS BA SR F BR I IOD AL CR AG AU PB PT SOD CLA POT CAL CES LIT RUB ' +
    'BAR ZN2 MG2 CA2 NA+ CL- K+ LI+ RB+ CS+ F- BR- I- MG2+ CA2+ ZN2+ IB NIO CIO KIO NAI CLI'));
  const STEROL = new Set(W('CHL1 CHOL CHL ERG SITO STIG LANO'));
  const LIPID = new Set([...STEROL, ...W('POPC POPE POPS POPG POPA POPI DPPC DOPC DMPC DLPC DSPC DAPC DEPC DIPC DYPC DPPE DOPE DMPE DLPE DSPE DYPE DPPS DOPS DMPS ' +
    'DLPS DSPS DPPG DOPG DMPG DLPG DSPG DYPG DPPA DOPA DMPA DLPA DSPA PLPC PLPE PLPS PLPG PLPA SOPC SOPE SOPS SOPG SOPA SAPC SAPE SAPS SDPC SDPE SDPS PAPC PAPE ' +
    'PAPS OSPC OSPE YOPC YOPE YOPG YOPS PMPE PMPG PVPE PVPG QMPE PPPE PYPE PYPG PSPC PSM SSM NSM LSM OSM BSM ASM CER CER2 CER160 CER180 TOCL TOCL1 TOCL2 TLCL ' +
    'TLCL1 TLCL2 POCL1 POCL2 CDL0 CDL1 CDL2 PIP PIP2 PIP3 POP2 POP3 SAPI SAPI24 SAPI25 SAP2 DPC SDS DDM LMT BOG LDAO LDA C12E8 DHPC OLEO PALM PALA STEA MYR LAU ' +
    // Amber Lipid14/21 residues (heads and tails are separate residues)
    'PC PE PS PGR PGS PH- PI SPM PA OL ST MY LA LAL LEO LEN AR DHA')]);
  // PDB chemical-component ids of monosaccharides -> stereo family (kind and modifications come from the topology)
  const CCD = {
    GLC: 'Glc', BGC: 'Glc', AGC: 'Glc', NAG: 'Glc', NDG: 'Glc', GCS: 'Glc', PA1: 'Glc', SGN: 'Glc', NGS: 'Glc', BDP: 'Glc', GCU: 'Glc', GC4: 'Glc', G6S: 'Glc', QUI: 'Glc',
    MAN: 'Man', BMA: 'Man', BM3: 'Man', MAV: 'Man', BEM: 'Man', RAM: 'Man', RM4: 'Man', XXR: 'Man',
    GAL: 'Gal', GLA: 'Gal', NGA: 'Gal', A2G: 'Gal', ASG: 'Gal', NG6: 'Gal', GZL: 'Gal', GAD: 'Gal', ADA: 'Gal', GTR: 'Gal', FUC: 'Gal', FUL: 'Gal', FCA: 'Gal', FCB: 'Gal',
    XYS: 'Xyl', XYP: 'Xyl', XYZ: 'Xyl', LXC: 'Xyl', ARA: 'Ara', ARB: 'Ara', AHR: 'Ara', FUB: 'Ara', BXY: 'Ara', RIB: 'Rib', LDY: 'Lyx',
    IDR: 'Ido', IDS: 'Ido', IDO: 'Ido', SIA: 'Neu5Ac', SLB: 'Neu5Ac', NGC: 'Neu5Gc', NGE: 'Neu5Gc', KDN: 'Kdn', KDO: 'Kdo', FRU: 'Fru', BDF: 'Fru', SUC: null,
    MAL: null, LAT: null, TRE: null, CBI: null, UAP: 'dUA', GCD: 'dUA', ALL: 'All', AFD: 'All', SHD: 'Alt', TAL: 'Tal', GUP: 'Gul'
  };
  // CHARMM / CHARMM-GUI glycan residue names (AGLC, BGLCNA, AIDOA, BGLCA…; truncated to 4-5 characters in PDB/GRO files)
  const CHARMM_SUGAR = /^[AB](GLC|GAL|MAN|FUC|XYL|IDO|ALT|ALL|GUL|TAL|RHM|RHA|ARB|ARA|LYX|RIB|FRU|NE5|SIA|NEU|KDO|KDN|QUI)/;
  const CHARMM_FAM = { GLC: 'Glc', GAL: 'Gal', MAN: 'Man', FUC: 'Gal', XYL: 'Xyl', IDO: 'Ido', ALT: 'Alt', ALL: 'All', GUL: 'Gul', TAL: 'Tal', RHM: 'Man', RHA: 'Man',
    ARB: 'Ara', ARA: 'Ara', LYX: 'Lyx', RIB: 'Rib', FRU: 'Fru', NE5: 'Neu5Ac', SIA: 'Neu5Ac', NEU: 'Neu5Ac', KDO: 'Kdo', KDN: 'Kdn', QUI: 'Glc' };

  const C_OTHER = 0, C_PROTEIN = 1, C_NUCLEIC = 2, C_WATER = 3, C_ION = 4, C_LIPID = 5, C_SUGAR = 6;
  const CLASS_IDS = ['ligand', 'protein', 'nucleic', 'water', 'ions', 'lipid', 'carbohydrate'];
  const NAME_CLASS = new Map();
  for (const [set, c] of [[LIPID, C_LIPID], [NUCLEIC, C_NUCLEIC], [AMINO, C_PROTEIN], [IONS, C_ION], [WATER, C_WATER]]) for (const r of set) NAME_CLASS.set(r, c);
  for (const r of Object.keys(CCD)) NAME_CLASS.set(r, C_SUGAR);
  const PROT_BB = new Set(W('N CA C O OXT OT1 OT2 OC1 OC2 O1 O2 H HN H1 H2 H3 HT1 HT2 HT3 HA HA1 HA2 HA3'));
  const NUC_BB = new Set(W("P OP1 OP2 OP3 O1P O2P O3P O5' C5' C4' C3' O3' O5* C5* C4* C3* O3* H5' H5'' H4' H3' H5T H3T"));
  // nonmetals and metalloids; every other element is a metal
  const NONMETAL = new Set([0, 1, 2, 5, 6, 7, 8, 9, 10, 14, 15, 16, 17, 18, 32, 33, 34, 35, 36, 51, 52, 53, 54, 85, 86]);
  const isMetal = z => z > 0 && !NONMETAL.has(z);
  const HALOGEN = new Set([9, 17, 35, 53]);

  // Renderer style codes (identical to BM.Renderer.STYLE)
  const CODE = { hidden: 0, ballstick: 1, sticks: 2, spacefill: 3, wireframe: 4, lines: 5, ghost: 6 };
  const STYLES = [
    { id: 'ballstick', label: 'Ball & stick', kind: 'atom' }, { id: 'sticks', label: 'Sticks', kind: 'atom' }, { id: 'lines', label: 'Lines', kind: 'atom' },
    { id: 'wireframe', label: 'Wire (thin tubes)', kind: 'atom' }, { id: 'spacefill', label: 'Spacefill', kind: 'atom' }, { id: 'ghost', label: 'Ghost (translucent)', kind: 'atom' },
    { id: 'hidden', label: 'Hidden', kind: 'atom' }, { id: 'lipid', label: 'Lipid (tails as lines, P/N spheres)', kind: 'atom' },
    { id: 'cartoon', label: 'Cartoon', kind: 'overlay' }, { id: 'trace', label: 'Trace (CA / P tube)', kind: 'overlay' }, { id: 'rings', label: 'Filled rings', kind: 'overlay' },
    { id: 'snfg', label: 'SNFG symbols', kind: 'overlay' }, { id: 'markers', label: 'Markers (rings)', kind: 'overlay' }
  ];
  const STYLE_BY_ID = Object.fromEntries(STYLES.map(s => [s.id, s]));
  const COLORS = [
    { id: 'element', label: 'Element' }, { id: 'layer', label: 'QM/MM layer' }, { id: 'chain', label: 'Chain' }, { id: 'residue', label: 'Residue name' },
    { id: 'restype', label: 'Residue type' }, { id: 'ss', label: 'Secondary structure' }, { id: 'rainbow', label: 'Rainbow (N→C)' }, { id: 'sugar', label: 'Monosaccharide (SNFG)' },
    { id: 'charge', label: 'Charge' }, { id: 'bfactor', label: 'B-factor' }, { id: 'occupancy', label: 'Occupancy' }, { id: 'uniform', label: 'Uniform colour' }
  ];

  // ------------------------------------------------------------------ colours
  const hex = h => [parseInt(h.slice(1, 3), 16) / 255, parseInt(h.slice(3, 5), 16) / 255, parseInt(h.slice(5, 7), 16) / 255];
  const toHex = c => '#' + c.map(v => Math.round(Math.max(0, Math.min(1, v)) * 255).toString(16).padStart(2, '0')).join('');
  const SN = { white: [1, 1, 1], blue: hex('#0090bc'), green: hex('#00a651'), yellow: hex('#ffd400'), lightblue: hex('#8fcce9'), pink: hex('#f69ea1'),
    purple: hex('#a54399'), brown: hex('#a17a4d'), orange: hex('#f47920'), red: hex('#ed1c24') };
  const FAM_COL = { Glc: SN.blue, Man: SN.green, Gal: SN.yellow, Gul: SN.orange, Alt: SN.pink, All: SN.purple, Tal: SN.lightblue, Ido: SN.brown };
  const PEN_COL = { Ara: SN.green, Lyx: SN.yellow, Xyl: SN.orange, Rib: SN.pink };
  const DEOXY = { Glc: ['Qui', SN.blue], Man: ['Rha', SN.green], Gal: ['Fuc', SN.red], Alt: ['6dAlt', SN.pink], Tal: ['6dTal', SN.lightblue], Gul: ['6dGul', SN.orange] };
  const LAYER_COL = { H: hex('#ffb347'), M: hex('#4fd1c5'), L: hex('#8c9bb5') };
  const CHAIN_COL = ['#4f8cff', '#f2994a', '#3cc27a', '#eb5757', '#a66cf0', '#4fd0f0', '#f2c94c', '#e56bb5', '#8bd36b', '#ff9e7a', '#6d8cf5', '#c7a35c'].map(hex);
  const SS_COL = { helix: hex('#e2507a'), sheet: hex('#f4c542'), coil: hex('#d2d9e3'), coilLight: hex('#8a94a6') };
  const RT_COL = { nonpolar: hex('#c3c9d1'), polar: hex('#5ec26a'), acidic: hex('#e5534b'), basic: hex('#4a86e8'), nucleic: hex('#f39c4a'), lipid: hex('#c8a872'),
    carbohydrate: hex('#4fb3d9'), water: hex('#7fd3f0'), ions: hex('#b48ef0'), ligand: hex('#7ed957') };
  const RES_TYPE = {};
  for (const r of W('ALA VAL LEU ILE MET PHE TRP PRO GLY MSE NLE')) RES_TYPE[r] = 'nonpolar';
  for (const r of W('SER THR ASN GLN TYR CYS CYX CYM SEC HYP')) RES_TYPE[r] = 'polar';
  for (const r of W('ASP GLU ASH GLH ASPH GLUH SEP TPO PTR')) RES_TYPE[r] = 'acidic';
  for (const r of W('LYS ARG HIS HID HIE HIP HSD HSE HSP HISD HISE HISH LYN ARN')) RES_TYPE[r] = 'basic';
  const AA_COL = { ALA: '#c8c8c8', ARG: '#3d5df0', ASN: '#36c9c0', ASP: '#e63c3c', CYS: '#e8d43a', GLN: '#2fb8b0', GLU: '#d23a3a', GLY: '#ebebeb', HIS: '#8282d2', ILE: '#1f9e3a',
    LEU: '#2ab04a', LYS: '#2a5ad8', MET: '#e6c23a', PHE: '#4646a8', PRO: '#dc9682', SER: '#fa9600', THR: '#f0a030', TRP: '#b45ab4', TYR: '#3232aa', VAL: '#36a84e' };
  const NUC_COL = { A: hex('#e04f4f'), G: hex('#3fae5a'), C: hex('#f2c94c'), T: hex('#4f7fe0'), U: hex('#49c1d6') };
  const elColor = z => (z === 0 ? [1.0, 0.42, 0.78] : el.color[z] || [1, 0.08, 0.6]);
  function hashCol(s) { let h = 2166136261; for (let k = 0; k < s.length; k++) h = Math.imul(h ^ s.charCodeAt(k), 16777619); return CHAIN_COL[(h >>> 0) % CHAIN_COL.length]; }
  function hsv(h, s, v) { const f = (n, k = (n + h / 60) % 6) => v - v * s * Math.max(Math.min(k, 4 - k, 1), 0); return [f(5), f(3), f(1)]; }
  const ramp = t => (t < 0.5 ? [0.23 + 1.54 * t, 0.44 + 1.12 * t, 0.85 + 0.3 * t] : [1, 1 - 1.46 * (t - 0.5), 1 - 1.54 * (t - 0.5)]); // blue - white - red

  // ------------------------------------------------------------------ context
  function context(o) {
    const f = o.frame;
    return { frame: f, n: f.Z.length, Z: f.Z, xyz: o.xyz || f.xyz, bonds: o.bonds || new Uint32Array(0), layers: o.layers || null, ssRecords: o.ssRecords || null,
      charges: o.charges || null, selection: o.selection || [], frozen: o.frozen || [], lightBg: !!o.lightBg, cache: new Map() };
  }
  /** Bond graph as CSR {off, nbr} (cached on the context). */
  function adjacency(ctx) {
    if (ctx._adj) return ctx._adj;
    const n = ctx.n, b = ctx.bonds, off = new Int32Array(n + 1);
    for (let k = 0; k + 1 < b.length; k += 2) { const i = b[k], j = b[k + 1]; if (i < n && j < n && i !== j) { off[i + 1]++; off[j + 1]++; } }
    for (let i = 0; i < n; i++) off[i + 1] += off[i];
    const fill = off.slice(0, n), nbr = new Int32Array(off[n]);
    for (let k = 0; k + 1 < b.length; k += 2) { const i = b[k], j = b[k + 1]; if (i < n && j < n && i !== j) { nbr[fill[i]++] = j; nbr[fill[j]++] = i; } }
    ctx._adj = { off, nbr };
    return ctx._adj;
  }

  // ------------------------------------------------------------------ classification
  const topoCache = new WeakMap(); // frame.Z array -> {n, names, res, ids, ch, nb, cls}
  function classify(ctx) {
    if (ctx._cls) return ctx._cls;
    const f = ctx.frame, Z = f.Z;
    const e = topoCache.get(Z);
    if (e && e.n === ctx.n && e.names === (f.atomNames || null) && e.res === (f.resNames || null) && e.ids === (f.resIds || null) && e.ch === (f.chains || null) &&
      (f.resNames || e.nb === ctx.bonds.length)) { ctx._cls = e.cls; return e.cls; }
    const cls = buildClass(ctx);
    topoCache.set(Z, { n: ctx.n, names: f.atomNames || null, res: f.resNames || null, ids: f.resIds || null, ch: f.chains || null, nb: ctx.bonds.length, cls });
    ctx._cls = cls;
    return cls;
  }
  function buildClass(ctx) {
    const f = ctx.frame, Z = f.Z, n = ctx.n;
    const names = f.atomNames && f.atomNames.length === n ? f.atomNames : null;
    const rn = f.resNames && f.resNames.length === n ? f.resNames : null;
    const ri = f.resIds && f.resIds.length === n ? f.resIds : null, chs = f.chains && f.chains.length === n ? f.chains : null;
    const resOf = new Int32Array(n);
    let nRes = 0;
    if (rn) {
      let pn, pi, pc;
      for (let i = 0; i < n; i++) {
        const a = rn[i], b0 = ri ? ri[i] : 0, b = b0 === b0 ? b0 : null, c = chs ? chs[i] : null;
        if (i === 0 || a !== pn || b !== pi || c !== pc) { nRes++; pn = a; pi = b; pc = c; }
        resOf[i] = nRes - 1;
      }
    } else if (n) {
      // no residue information (QC outputs, XYZ…): residues are the molecules (connected components)
      const par = new Int32Array(n);
      for (let i = 0; i < n; i++) par[i] = i;
      const find = x => { while (par[x] !== x) { par[x] = par[par[x]]; x = par[x]; } return x; };
      const b = ctx.bonds;
      for (let k = 0; k + 1 < b.length; k += 2) {
        if (b[k] >= n || b[k + 1] >= n) continue;
        const x = find(b[k]), y = find(b[k + 1]);
        if (x !== y) { if (x < y) par[y] = x; else par[x] = y; }
      }
      const id = new Int32Array(n).fill(-1);
      for (let i = 0; i < n; i++) { const r = find(i); if (id[r] < 0) id[r] = nRes++; resOf[i] = id[r]; }
    }
    const resOff = new Int32Array(nRes + 1);
    for (let i = 0; i < n; i++) resOff[resOf[i] + 1]++;
    for (let r = 0; r < nRes; r++) resOff[r + 1] += resOff[r];
    const fill = resOff.slice(0, nRes), resAtoms = new Int32Array(n);
    for (let i = 0; i < n; i++) resAtoms[fill[resOf[i]]++] = i;
    const resName = new Array(nRes), resId = new Array(nRes), chain = new Array(nRes), resClass = new Uint8Array(nRes);
    const cls = { n, nRes, resOf, resOff, resAtoms, resName, resId, chain, resClass, atomClass: new Uint8Array(n), backbone: new Uint8Array(n), names, rings: new Map(), masks: new Map() };
    const unknown = [];
    for (let r = 0; r < nRes; r++) {
      const a = resAtoms[resOff[r]];
      const nm = rn ? String(rn[a] === undefined || rn[a] === null ? '' : rn[a]).trim().toUpperCase() : '';
      resName[r] = nm; resId[r] = ri ? ri[a] : r + 1; chain[r] = chs && chs[a] !== null && chs[a] !== undefined ? String(chs[a]).trim() : '';
      let c = nm ? NAME_CLASS.get(nm) : undefined;
      if (c === undefined && nm && CHARMM_SUGAR.test(nm)) c = C_SUGAR;
      if (c === undefined) unknown.push(r); else resClass[r] = c;
    }
    // topology fallbacks for residues whose names say nothing
    const lipidCand = [];
    for (const r of unknown) {
      const c = fallbackClass(ctx, cls, r);
      if (c === -1) lipidCand.push(r); else resClass[r] = c;
    }
    if (lipidCand.length) {
      const cnt = new Map();
      for (const r of lipidCand) cnt.set(resName[r], (cnt.get(resName[r]) || 0) + 1);
      for (const r of lipidCand) resClass[r] = cnt.get(resName[r]) >= 3 ? C_LIPID : C_OTHER;
    }
    const counts = {};
    for (const id of CLASS_IDS) counts[id] = { res: 0, atoms: 0 };
    for (let r = 0; r < nRes; r++) { const k = CLASS_IDS[resClass[r]]; counts[k].res++; counts[k].atoms += resOff[r + 1] - resOff[r]; }
    cls.counts = counts;
    for (let i = 0; i < n; i++) {
      const c = resClass[resOf[i]];
      cls.atomClass[i] = c;
      if (names && (c === C_PROTEIN || c === C_NUCLEIC)) {
        const nm = String(names[i] || '').trim().toUpperCase();
        if (c === C_PROTEIN ? PROT_BB.has(nm) : NUC_BB.has(nm)) cls.backbone[i] = 1;
      }
    }
    return cls;
  }
  function atomNameSet(cls, r) {
    const s = new Set();
    if (!cls.names) return s;
    for (let t = cls.resOff[r]; t < cls.resOff[r + 1]; t++) s.add(String(cls.names[cls.resAtoms[t]] || '').trim().toUpperCase());
    return s;
  }
  /** Class of a residue with an unknown name; -1 = lipid candidate (a lipid only when its name occurs ≥ 3 times). */
  function fallbackClass(ctx, cls, r) {
    const Z = ctx.Z, a0 = cls.resOff[r], a1 = cls.resOff[r + 1];
    let heavy = 0, nO = 0, nH = 0, nC = 0, lastHeavy = -1;
    for (let t = a0; t < a1; t++) { const z = Z[cls.resAtoms[t]]; if (z === 1) nH++; else if (z > 1) { heavy++; lastHeavy = cls.resAtoms[t]; if (z === 8) nO++; if (z === 6) nC++; } }
    if (heavy === 1 && nO === 1 && nH <= 2) return C_WATER;
    if (heavy === 1 && nH === 0 && a1 - a0 === 1) { const z = Z[lastHeavy]; if (isMetal(z) || HALOGEN.has(z)) { const A = adjacency(ctx); if (A.off[lastHeavy + 1] === A.off[lastHeavy]) return C_ION; } }
    if (cls.names) {
      const s = atomNameSet(cls, r);
      if (s.has('CA') && s.has('N') && s.has('C')) return C_PROTEIN;
      if ((s.has("C4'") || s.has('C4*')) && (s.has("C1'") || s.has('C1*')) && (s.has('P') || s.has("O4'") || s.has('O4*'))) return C_NUCLEIC;
    }
    if (nC >= 4 && nO >= 3 && residueRings(ctx, cls, r).length) return C_SUGAR;
    if (heavy >= 14 && nC >= 10 && chainCarbons(ctx, cls, r) >= 8) return -1;
    return C_OTHER;
  }
  /** Carbons whose heavy neighbours are exactly two carbons (aliphatic chain / ring CH2). */
  function chainCarbons(ctx, cls, r) {
    const Z = ctx.Z, A = adjacency(ctx);
    let k = 0;
    for (let t = cls.resOff[r]; t < cls.resOff[r + 1]; t++) {
      const a = cls.resAtoms[t];
      if (Z[a] !== 6) continue;
      let hc = 0, other = 0;
      for (let q = A.off[a]; q < A.off[a + 1]; q++) { const z = Z[A.nbr[q]]; if (z === 6) hc++; else if (z > 1) other++; }
      if (hc === 2 && !other) k++;
    }
    return k;
  }
  /** Sugar rings of residue r: [ringO, C, C, …] for 5- and 6-membered rings of carbons with one ring oxygen and at least
   *  two ring carbons carrying an exocyclic O or N (cached on the classification). */
  function residueRings(ctx, cls, r) {
    if (cls.rings.has(r)) return cls.rings.get(r);
    const Z = ctx.Z, A = adjacency(ctx), resOf = cls.resOf, out = [];
    for (let t = cls.resOff[r]; t < cls.resOff[r + 1]; t++) {
      const o = cls.resAtoms[t];
      if (Z[o] !== 8) continue;
      const cs = [];
      for (let q = A.off[o]; q < A.off[o + 1]; q++) { const c = A.nbr[q]; if (Z[c] === 6 && resOf[c] === r) cs.push(c); }
      if (cs.length !== 2) continue;
      const path = carbonPath(Z, A, resOf, r, cs[0], cs[1], 5);
      if (!path || path.length < 4) continue;
      const inR = new Set(path);
      let exo = 0;
      for (const c of path) {
        for (let q = A.off[c]; q < A.off[c + 1]; q++) { const x = A.nbr[q]; if ((Z[x] === 8 || Z[x] === 7) && x !== o && !inR.has(x)) { exo++; break; } }
      }
      if (exo >= 2) out.push([o].concat(path));
    }
    cls.rings.set(r, out);
    return out;
  }
  function carbonPath(Z, A, resOf, r, s, t, maxLen) {
    const prev = new Map([[s, -1]]);
    let front = [s];
    for (let depth = 1; depth < maxLen && front.length; depth++) {
      const next = [];
      for (const u of front) {
        for (let q = A.off[u]; q < A.off[u + 1]; q++) {
          const v = A.nbr[q];
          if (Z[v] !== 6 || resOf[v] !== r || prev.has(v)) continue;
          prev.set(v, u);
          if (v === t) { const p = []; for (let x = v; x !== -1; x = prev.get(x)) p.push(x); return p.reverse(); }
          next.push(v);
        }
      }
      front = next;
    }
    return null;
  }

  // ------------------------------------------------------------------ monosaccharides (SNFG)
  const HEX_FAM = { dud: 'Glc', duu: 'Gal', uud: 'Man', ddd: 'All', udd: 'Alt', ddu: 'Gul', udu: 'Ido', uuu: 'Tal' };
  const PEN_FAM = { dud: 'Xyl', udu: 'Xyl', udd: 'Ara', duu: 'Ara', ddd: 'Rib', uuu: 'Rib', uud: 'Lyx', ddu: 'Lyx' };
  function nameFamily(nm) {
    if (Object.prototype.hasOwnProperty.call(CCD, nm)) return CCD[nm];
    const m = CHARMM_SUGAR.exec(nm);
    return m ? CHARMM_FAM[m[1]] : null;
  }
  /** Every sugar ring of the structure: {res, ring: [O, C1…], centroid?, kind, family, dl, key, label, shape, color, color2, inv, mods}. */
  function sugars(ctx) {
    const cls = classify(ctx);
    if (cls.sugars && cls.sugarsN === ctx.bonds.length) return cls.sugars;
    const out = [];
    for (let r = 0; r < cls.nRes; r++) {
      if (cls.resClass[r] !== C_SUGAR) continue;
      const rings = residueRings(ctx, cls, r);
      for (const ring of rings) out.push(identifyRing(ctx, cls, r, ring, rings.length === 1 ? nameFamily(cls.resName[r]) : null));
    }
    cls.sugars = out; cls.sugarsN = ctx.bonds.length;
    cls.sugarOfRes = new Map();
    for (const s of out) if (!cls.sugarOfRes.has(s.res)) cls.sugarOfRes.set(s.res, s);
    return out;
  }
  function identifyRing(ctx, cls, r, ring, nameFam) {
    const Z = ctx.Z, xyz = ctx.xyz, A = adjacency(ctx), resOf = cls.resOf;
    const inRing = new Set(ring), O5 = ring[0];
    const nb = a => { const s = []; for (let q = A.off[a]; q < A.off[a + 1]; q++) s.push(A.nbr[q]); return s; };
    const exo = a => nb(a).filter(x => !inRing.has(x) && Z[x] > 1);
    const isON = x => Z[x] === 7 || Z[x] === 8;
    const countO = c => nb(c).filter(x => Z[x] === 8).length;
    const bondedTo = (x, z) => nb(x).some(y => Z[y] === z);
    const score = c => (exo(c).some(isON) ? 2 : 0) - (exo(c).some(x => Z[x] === 6) ? 1 : 0);
    const C = ring.slice(1);
    if (score(C[C.length - 1]) > score(C[0])) C.reverse();
    const size = ring.length;
    const info = { res: r, ring: [O5].concat(C), size, anomeric: C[0], kind: 'hex', family: null, dl: '', mods: [], amine: false, acyl: null, ns: false, furanose: size === 5 };
    const acylOf = N => {
      for (const c of nb(N)) {
        if (Z[c] !== 6 || inRing.has(c) || !nb(c).some(y => Z[y] === 8)) continue;
        const m = nb(c).find(y => Z[y] === 6 && y !== N);
        if (m !== undefined) return nb(m).some(y => Z[y] === 8) ? 'Gc' : 'Ac';
      }
      return null;
    };
    const anomC = exo(C[0]).filter(x => Z[x] === 6);
    let first = 1; // sugar numbering of C[0]
    if (anomC.length && exo(C[0]).some(isON)) {
      first = 2;
      const c1 = anomC[0];
      if (countO(c1) >= 2 && size === 6) {
        const n5 = exo(C[3]).find(x => Z[x] === 7);
        if (n5 !== undefined) { info.kind = 'sia'; info.acyl = acylOf(n5); }
        else {
          // side-chain carbons beyond the last ring carbon: 3 for Kdn (C7-C9), 2 for Kdo (C7-C8)
          let len = 0, prev = C[C.length - 1], cur = exo(prev).find(x => Z[x] === 6);
          while (cur !== undefined && len < 6) { len++; const nx = exo(cur).find(x => Z[x] === 6 && x !== prev); prev = cur; cur = nx; }
          info.kind = len >= 3 ? 'kdn' : 'kdo';
        }
      } else info.kind = 'ket';
    } else {
      const last = C[C.length - 1], exoL = exo(last).filter(x => Z[x] === 6);
      if (size === 6) {
        if (!exoL.length) info.kind = 'pen';
        else { const nO = countO(exoL[0]); info.kind = nO >= 2 ? 'uronic' : nO === 0 ? 'deoxy' : 'hex'; info.c6 = exoL[0]; }
      } else if (exoL.length) {
        const c5 = exoL[0], beyond = nb(c5).filter(x => Z[x] === 6 && x !== last);
        info.kind = beyond.length ? 'hex' : 'pen';
        if (beyond.length) info.c6 = beyond[0];
      } else info.kind = 'pen';
      const n2 = exo(C[1]).find(x => Z[x] === 7);
      if (n2 !== undefined) { info.amine = true; info.acyl = acylOf(n2); info.ns = bondedTo(n2, 16); }
      if (size === 6) {
        const side = k => {
          const c = C[k], prev = k === 0 ? O5 : C[k - 1], next = k === C.length - 1 ? O5 : C[k + 1];
          const ex = exo(c);
          let sub = k === C.length - 1 ? ex.find(x => Z[x] === 6) : ex.find(isON);
          if (sub === undefined) sub = ex[0];
          if (sub === undefined) return 0;
          const P = i => [xyz[3 * i], xyz[3 * i + 1], xyz[3 * i + 2]], pc = P(c), pp = P(prev), pn = P(next), ps = P(sub);
          const a = [pc[0] - pp[0], pc[1] - pp[1], pc[2] - pp[2]], b = [pn[0] - pc[0], pn[1] - pc[1], pn[2] - pc[2]];
          const nl = [a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]];
          const d = (ps[0] - pc[0]) * nl[0] + (ps[1] - pc[1]) * nl[1] + (ps[2] - pc[2]) * nl[2];
          return Math.abs(d) < 1e-6 ? 0 : d < 0 ? 1 : -1; // +1 = up in the Haworth projection
        };
        const s2 = side(1), s3 = side(2), s4 = side(3), s5 = info.kind === 'pen' ? 0 : side(4);
        const code = v => v.map(x => (x > 0 ? 'u' : 'd')).join('');
        if (info.kind === 'pen') { if (s2 && s3 && s4) { info.family = PEN_FAM[code([s2, s3, s4])] || null; info.dl = s4 < 0 ? 'D' : 'L'; } }
        else if (s2 && s3 && s4 && s5) { info.family = HEX_FAM[code(s5 > 0 ? [s2, s3, s4] : [-s2, -s3, -s4])] || null; info.dl = s5 > 0 ? 'D' : 'L'; }
      }
    }
    if (nameFam === 'Neu5Ac' || nameFam === 'Neu5Gc') { if (info.kind !== 'sia') { info.kind = 'sia'; info.acyl = nameFam === 'Neu5Gc' ? 'Gc' : 'Ac'; } }
    else if (nameFam === 'Kdn' || nameFam === 'Kdo') info.kind = nameFam.toLowerCase();
    else if (nameFam === 'Fru') { info.kind = 'ket'; info.family = 'Fru'; }
    else if (nameFam === 'dUA') { info.kind = 'uronic'; info.family = null; info.special = 'ΔUA'; }
    else if (nameFam) info.family = nameFam;
    // O-sulfates / phosphates on ring carbons and on the exocyclic C6 (C5 of furanoses)
    const mods = new Set();
    C.forEach((c, k) => { for (const x of exo(c)) if (Z[x] === 8 && resOf[x] === r) { if (bondedTo(x, 16)) mods.add((k + first) + 'S'); else if (bondedTo(x, 15)) mods.add((k + first) + 'P'); } });
    if (info.c6 !== undefined) for (const x of nb(info.c6)) if (Z[x] === 8 && bondedTo(x, 16)) mods.add('6S');
    info.mods = [...mods].sort((a, b) => parseInt(a, 10) - parseInt(b, 10));
    Object.assign(info, snfgType(info));
    return info;
  }
  function snfgType(info) {
    const f = info.family;
    let key, shape, color, inv = false;
    switch (info.kind) {
      case 'sia': key = info.acyl === 'Gc' ? 'Neu5Gc' : info.acyl === 'Ac' ? 'Neu5Ac' : 'Neu'; shape = 'diamond'; color = key === 'Neu5Gc' ? SN.lightblue : key === 'Neu5Ac' ? SN.purple : SN.brown; break;
      case 'kdn': key = 'Kdn'; shape = 'diamond'; color = SN.green; break;
      case 'kdo': key = 'Kdo'; shape = 'hexagon'; color = SN.yellow; break;
      case 'ket': key = f === 'Fru' ? 'Fru' : 'Ket'; shape = 'pentagon'; color = f === 'Fru' ? SN.green : SN.white; break;
      case 'pen': key = f || 'Pen'; shape = 'star'; color = PEN_COL[f] || SN.white; break;
      case 'deoxy': { const d = DEOXY[f]; key = d ? d[0] : 'dHex'; shape = 'cone'; color = d ? d[1] : SN.white; break; }
      case 'uronic': key = info.special || (f || 'Hex') + 'A'; shape = 'ddiamond'; color = (!info.special && FAM_COL[f]) || SN.white; inv = info.dl === 'L' || f === 'Ido'; break;
      default:
        if (info.amine) { const nac = info.acyl === 'Ac'; key = (f || 'Hex') + (nac ? 'NAc' : info.acyl === 'Gc' ? 'NGc' : 'N'); shape = nac ? 'cube' : 'dcube'; }
        else { key = f || 'Hex'; shape = 'sphere'; }
        color = FAM_COL[f] || SN.white;
    }
    let label = key;
    if (info.furanose && f && info.kind !== 'ket') label = label.replace(f, f + 'f');
    if (info.amine && info.ns && /N$/.test(label)) label += 'S';
    label += info.mods.join('');
    return { key, label, shape, color, color2: SN.white, inv, known: color !== SN.white };
  }

  // ------------------------------------------------------------------ secondary structure (via BM.cartoon)
  /** Per-residue secondary structure (0 coil, 1 helix, 2 sheet) or null when BM.cartoon is not loaded. */
  function secondary(ctx) {
    const cls = classify(ctx);
    if (cls.ss !== undefined) return cls.ss;
    cls.ss = BM.cartoon && BM.cartoon.secondary ? BM.cartoon.secondary(ctx) : null;
    return cls.ss;
  }

  // ------------------------------------------------------------------ selection language
  class SelectionError extends Error {
    constructor(msg, pos, end) { super(msg); this.name = 'SelectionError'; this.pos = pos || 0; this.end = end === undefined ? this.pos : end; }
  }
  const FLAGS = { all: 'all', everything: 'all', none: 'none', protein: 'protein', backbone: 'backbone', sidechain: 'sidechain', nucleic: 'nucleic', dna: 'nucleic', rna: 'nucleic',
    water: 'water', waters: 'water', solvent: 'water', ions: 'ions', ion: 'ions', lipid: 'lipid', lipids: 'lipid', carbohydrate: 'carbohydrate', carbohydrates: 'carbohydrate',
    sugar: 'carbohydrate', sugars: 'carbohydrate', glycan: 'carbohydrate', glycans: 'carbohydrate', ligand: 'ligand', ligands: 'ligand', hydrogen: 'hydrogen', hydrogens: 'hydrogen',
    heavy: 'heavy', metal: 'metal', metals: 'metal', frozen: 'frozen', active: 'active', link: 'link', selected: 'selected', selection: 'selected',
    helix: 'helix', sheet: 'sheet', strand: 'sheet', coil: 'coil', qm: 'layerH', high: 'layerH', medium: 'layerM', qm2: 'layerM', low: 'layerL', mm: 'layerL' };
  const ARGS = { chain: 'chain', segid: 'chain', segname: 'chain', resname: 'resname', resn: 'resname', resid: 'resid', resi: 'resid', resnum: 'resid', name: 'name',
    element: 'element', elem: 'element', index: 'index', serial: 'serial', layer: 'layer' };
  const STOP = new Set(['and', 'or', 'not', 'within', 'same']);
  const KEYWORDS = ['all', 'none', 'protein', 'backbone', 'sidechain', 'nucleic', 'water', 'ions', 'lipid', 'carbohydrate', 'sugar', 'glycan', 'ligand', 'hydrogen', 'heavy', 'metal',
    'helix', 'sheet', 'coil', 'layer', 'qm', 'mm', 'high', 'medium', 'low', 'frozen', 'active', 'link', 'selected', 'chain', 'resname', 'resid', 'name', 'element', 'index', 'serial',
    'within', 'same residue as', 'and', 'or', 'not'];
  const KW_HINT = { protein: 'amino-acid residues', backbone: 'N CA C O + nucleic P/sugar backbone', nucleic: 'DNA / RNA', water: 'HOH WAT SOL TIP3…', ions: 'NA CL SOD CLA POT CAL…',
    lipid: 'POPC CHL1 PSM, Lipid21 PA/PC/OL…', carbohydrate: 'sugars (alias sugar, glycan)', ligand: 'other hetero groups', layer: 'layer H M L (QM/MM)', link: 'link atoms',
    chain: 'chain A B', resname: 'resname POPC CHL1', resid: 'resid 10-20 35', name: 'name CA N C O (wildcards * ?)', element: 'element C N', index: 'index 0-99 (0-based)',
    serial: 'serial 1-100 (1-based)', within: 'within 5 of ligand', 'same residue as': 'same residue as within 5 of ligand', selected: 'atoms selected in the view',
    frozen: 'frozen atoms', active: 'active (not frozen) atoms', qm: 'layer H', mm: 'layer L', metal: 'metal atoms', helix: 'protein in helices', sheet: 'protein in sheets' };
  const LAYER_ALIAS = { h: 'H', high: 'H', qm: 'H', m: 'M', medium: 'M', qm2: 'M', mid: 'M', l: 'L', low: 'L', mm: 'L' };

  function tokenize(s) {
    const t = [];
    let i = 0;
    while (i < s.length) {
      const c = s[i];
      if (/\s/.test(c)) { i++; continue; }
      if (c === '(' || c === ')') { t.push({ v: c, p: i, e: i + 1, kind: c }); i++; continue; }
      if (c === '"') {
        const j = s.indexOf('"', i + 1);
        if (j < 0) throw new SelectionError('Missing closing quote "', i, s.length);
        t.push({ v: s.slice(i + 1, j), p: i, e: j + 1, kind: 'str' }); i = j + 1; continue;
      }
      let j = i;
      while (j < s.length && !/[\s()"]/.test(s[j])) j++;
      t.push({ v: s.slice(i, j), p: i, e: j, kind: 'word' }); i = j;
    }
    return t;
  }
  function lev(a, b) {
    const d = Array.from({ length: a.length + 1 }, (_, i) => [i]);
    for (let j = 1; j <= b.length; j++) d[0][j] = j;
    for (let i = 1; i <= a.length; i++) for (let j = 1; j <= b.length; j++) d[i][j] = Math.min(d[i - 1][j] + 1, d[i][j - 1] + 1, d[i - 1][j - 1] + (a[i - 1] === b[j - 1] ? 0 : 1));
    return d[a.length][b.length];
  }
  function suggest(w) {
    let best = null, bd = 3;
    for (const k of Object.keys(FLAGS).concat(Object.keys(ARGS), ['within', 'same'])) { const d = lev(w, k); if (d < bd) { bd = d; best = k; } }
    return best && bd <= (w.length > 4 ? 2 : 1) ? best : null;
  }
  const parseCache = new Map();
  /** Selection text -> AST. Precedence: or < and < not / within D of / same residue as < terms and parentheses. */
  function parse(text) {
    text = String(text === undefined || text === null ? '' : text);
    const hit = parseCache.get(text);
    if (hit) return hit;
    const toks = tokenize(text);
    if (!toks.length) throw new SelectionError('Empty selection: type e.g. protein, resname POPC or layer H', 0, 0);
    let k = 0;
    const peek = () => toks[k];
    const low = t => (t && t.kind === 'word' ? t.v.toLowerCase() : null);
    const fail = (msg, t) => { throw new SelectionError(msg, t ? t.p : text.length, t ? t.e : text.length); };
    const flags = { coords: false, sel: false, frozen: false };
    function orE() { let a = andE(); while (low(peek()) === 'or') { k++; if (!peek()) fail('"or" needs a selection after it'); a = { op: 'or', a, b: andE() }; } return a; }
    function andE() { let a = un(); while (low(peek()) === 'and') { k++; if (!peek()) fail('"and" needs a selection after it'); a = { op: 'and', a, b: un() }; } return a; }
    function un() {
      const t = peek();
      if (!t) fail('The selection ends too early');
      const w = low(t);
      if (w === 'not') { k++; if (!peek()) fail('"not" needs a selection after it', t); return { op: 'not', a: un() }; }
      if (w === 'within') {
        k++;
        const d = peek(), v = d && d.kind === 'word' ? Number(d.v) : NaN;
        if (!d || !isFinite(v) || v < 0) fail('"within" needs a distance in Å, e.g. within 5 of ligand', d || t);
        k++;
        if (low(peek()) !== 'of') fail('Expected "of" after "within ' + d.v + '"', peek() || d);
        k++;
        if (!peek()) fail('"within ' + d.v + ' of" needs a selection after it', d);
        flags.coords = true;
        return { op: 'within', d: v, a: un() };
      }
      if (w === 'same') {
        k++;
        if (low(peek()) !== 'residue') fail('Expected "same residue as …"', peek() || t);
        k++;
        if (low(peek()) !== 'as') fail('Expected "as" after "same residue"', peek() || t);
        k++;
        if (!peek()) fail('"same residue as" needs a selection after it', t);
        return { op: 'same', a: un() };
      }
      return prim();
    }
    function prim() {
      const t = peek();
      if (t.kind === '(') {
        k++;
        if (peek() && peek().kind === ')') fail('Empty parentheses', peek());
        const e = orE();
        if (!peek() || peek().kind !== ')') fail('Missing ")" for the "(" at position ' + (t.p + 1), t);
        k++;
        return e;
      }
      if (t.kind === ')') fail('Unexpected ")"', t);
      if (t.kind === 'str') fail('Expected a keyword before "' + t.v + '" (e.g. name "' + t.v + '")', t);
      const w = t.v.toLowerCase();
      if (FLAGS[w]) {
        k++;
        const f = FLAGS[w];
        if (f === 'selected') flags.sel = true;
        if (f === 'frozen' || f === 'active') flags.frozen = true;
        return f.startsWith('layer') ? { op: 'layer', vals: [f.slice(5)] } : { op: 'flag', f };
      }
      if (ARGS[w]) {
        k++;
        const kw = ARGS[w], vals = [];
        while (peek() && (peek().kind === 'str' || (peek().kind === 'word' && !STOP.has(low(peek()))))) { vals.push(peek()); k++; }
        if (!vals.length) fail('"' + t.v + '" needs ' + ({ chain: 'chain ids, e.g. chain A B', resname: 'residue names, e.g. resname POPC CHL1', resid: 'numbers or ranges, e.g. resid 10-20 35',
          name: 'atom names, e.g. name CA N C O', element: 'element symbols, e.g. element C N', index: 'numbers or ranges, e.g. index 0-99', serial: 'numbers or ranges, e.g. serial 1-100',
          layer: 'H, M or L' })[kw], t);
        return argNode(kw, vals, fail);
      }
      if (w === 'and' || w === 'or') fail('"' + t.v + '" needs a selection before it', t);
      if (w === 'of' || w === 'as' || w === 'to' || w === 'residue') fail('Unexpected "' + t.v + '"', t);
      if (/^-?\d/.test(w)) fail('Unexpected number "' + t.v + '": put it after a keyword such as resid, index or serial', t);
      const s = suggest(w);
      fail('Unknown keyword "' + t.v + '"' + (s ? ' — did you mean "' + s + '"?' : ''), t);
    }
    const ast = orE();
    if (k < toks.length) {
      const t = toks[k];
      fail(t.kind === ')' ? 'Unexpected ")"' : 'Unexpected "' + t.v + '": join terms with and / or', t);
    }
    ast.flags = flags;
    if (parseCache.size > 500) parseCache.clear();
    parseCache.set(text, ast);
    return ast;
  }
  function globRe(v) { return new RegExp('^' + v.replace(/[.+^${}()|[\]\\]/g, '\\$&').replace(/\*/g, '.*').replace(/\?/g, '.') + '$', 'i'); }
  function argNode(kw, toks, fail) {
    if (kw === 'resid' || kw === 'index' || kw === 'serial') {
      const ranges = [];
      for (let q = 0; q < toks.length; q++) {
        const v = toks[q].v;
        let m;
        if (/^-?\d+$/.test(v) && toks[q + 1] && toks[q + 1].v.toLowerCase() === 'to' && toks[q + 2] && /^-?\d+$/.test(toks[q + 2].v)) { ranges.push([+v, +toks[q + 2].v]); q += 2; }
        else if ((m = /^(-?\d+)(?:-|:|\.\.)(-?\d+)$/.exec(v))) ranges.push([+m[1], +m[2]]);
        else if (/^-?\d+$/.test(v)) ranges.push([+v, +v]);
        else fail('"' + kw + '" expects numbers or ranges like 10-20, not "' + v + '"', toks[q]);
      }
      return { op: kw, ranges: ranges.map(([a, b]) => (a <= b ? [a, b] : [b, a])) };
    }
    if (kw === 'layer') {
      const vals = toks.map(t => { const L = LAYER_ALIAS[t.v.toLowerCase()]; if (!L) fail('Layer must be H, M or L (or qm, mm, high, medium, low), not "' + t.v + '"', t); return L; });
      return { op: 'layer', vals };
    }
    if (kw === 'element') {
      const zs = toks.map(t => {
        const v = /^[A-Za-z]{1,3}$/.test(t.v) ? t.v[0].toUpperCase() + t.v.slice(1).toLowerCase() : '', z = v ? el.toZ(v) : -1;
        if (!(z >= 0) || el.symbol(z) !== v) fail('Unknown element "' + t.v + '"', t);
        return z;
      });
      return { op: 'element', zs: new Set(zs) };
    }
    const exact = new Set(), res = [];
    for (const t of toks) { if (/[*?]/.test(t.v)) res.push(globRe(t.v)); else exact.add(kw === 'chain' ? t.v : t.v.toUpperCase()); }
    return { op: kw, exact, res };
  }

  // ---- evaluation
  /** Mask of the atoms a selection (text or AST) picks; static selections are cached per topology. */
  function evaluate(sel, ctx) {
    const ast = typeof sel === 'string' ? parse(sel) : sel;
    const cls = classify(ctx);
    const dyn = ast.flags && (ast.flags.coords || ast.flags.sel || ast.flags.frozen);
    const key = typeof sel === 'string' ? sel : null;
    if (key !== null) {
      const store = dyn ? ctx.cache : cls.masks;
      const hit = store.get(key);
      if (hit) return hit;
      const m = evalNode(ast, ctx, cls);
      if (store.size > 200) store.clear();
      store.set(key, m);
      return m;
    }
    return evalNode(ast, ctx, cls);
  }
  function evalNode(node, ctx, cls) {
    const n = ctx.n, out = new Uint8Array(n), Z = ctx.Z, f = ctx.frame;
    switch (node.op) {
      case 'and': { const a = evalNode(node.a, ctx, cls), b = evalNode(node.b, ctx, cls); for (let i = 0; i < n; i++) out[i] = a[i] & b[i]; return out; }
      case 'or': { const a = evalNode(node.a, ctx, cls), b = evalNode(node.b, ctx, cls); for (let i = 0; i < n; i++) out[i] = a[i] | b[i]; return out; }
      case 'not': { const a = evalNode(node.a, ctx, cls); for (let i = 0; i < n; i++) out[i] = a[i] ? 0 : 1; return out; }
      case 'within': return withinMask(ctx, node.d, evalNode(node.a, ctx, cls));
      case 'same': {
        const a = evalNode(node.a, ctx, cls), hit = new Uint8Array(cls.nRes);
        for (let i = 0; i < n; i++) if (a[i]) hit[cls.resOf[i]] = 1;
        for (let i = 0; i < n; i++) out[i] = hit[cls.resOf[i]];
        return out;
      }
      case 'flag': return flagMask(node.f, ctx, cls);
      case 'layer': {
        const L = ctx.layers && ctx.layers.atomLayer;
        if (!L) return out;
        const want = new Set(node.vals);
        for (let i = 0; i < n && i < L.length; i++) if (want.has(L[i])) out[i] = 1;
        return out;
      }
      case 'element': for (let i = 0; i < n; i++) if (node.zs.has(Z[i])) out[i] = 1; return out;
      case 'index': case 'serial': {
        const off = node.op === 'serial' ? 1 : 0;
        for (const [a, b] of node.ranges) for (let i = Math.max(0, a - off); i <= b - off && i < n; i++) out[i] = 1;
        return out;
      }
      case 'resid': {
        const ids = f.resIds && f.resIds.length === n ? f.resIds : null;
        for (let i = 0; i < n; i++) { const v = ids ? ids[i] : cls.resId[cls.resOf[i]]; for (const [a, b] of node.ranges) if (v >= a && v <= b) { out[i] = 1; break; } }
        return out;
      }
      case 'name': case 'resname': case 'chain': {
        const src = node.op === 'name' ? f.atomNames : node.op === 'resname' ? f.resNames : f.chains;
        if (!src || src.length !== n) return out;
        const memo = new Map();
        for (let i = 0; i < n; i++) {
          const raw = src[i];
          let v = memo.get(raw);
          if (v === undefined) {
            const s = String(raw === null || raw === undefined ? '' : raw).trim(), u = node.op === 'chain' ? s : s.toUpperCase();
            v = node.exact.has(u) || node.res.some(re => re.test(s)) ? 1 : 0;
            memo.set(raw, v);
          }
          out[i] = v;
        }
        return out;
      }
    }
    return out;
  }
  function flagMask(fl, ctx, cls) {
    const n = ctx.n, out = new Uint8Array(n), Z = ctx.Z, ac = cls.atomClass, L = ctx.layers;
    const cl = { protein: C_PROTEIN, nucleic: C_NUCLEIC, water: C_WATER, ions: C_ION, lipid: C_LIPID, carbohydrate: C_SUGAR, ligand: C_OTHER }[fl];
    if (cl !== undefined) { for (let i = 0; i < n; i++) if (ac[i] === cl) out[i] = 1; return out; }
    switch (fl) {
      case 'all': out.fill(1); break;
      case 'none': break;
      case 'backbone': out.set(cls.backbone); break;
      case 'sidechain': for (let i = 0; i < n; i++) if (ac[i] === C_PROTEIN && !cls.backbone[i]) out[i] = 1; break;
      case 'hydrogen': for (let i = 0; i < n; i++) if (Z[i] === 1) out[i] = 1; break;
      case 'heavy': for (let i = 0; i < n; i++) if (Z[i] > 1) out[i] = 1; break;
      case 'metal': for (let i = 0; i < n; i++) if (isMetal(Z[i])) out[i] = 1; break;
      case 'frozen': case 'active': {
        const fr = new Uint8Array(n);
        if (L && L.frozen) for (let i = 0; i < n && i < L.frozen.length; i++) if (L.frozen[i]) fr[i] = 1;
        for (const i of ctx.frozen || []) if (i >= 0 && i < n) fr[i] = 1;
        if (fl === 'frozen') return fr;
        if (L && L.active) { for (let i = 0; i < n && i < L.active.length; i++) if (L.active[i] && !fr[i]) out[i] = 1; }
        else if (L || (ctx.frozen && ctx.frozen.length)) for (let i = 0; i < n; i++) out[i] = fr[i] ? 0 : 1;
        break;
      }
      case 'link': if (L && L.linkAtoms) for (const a of L.linkAtoms) if (a && a.atom >= 0 && a.atom < n) out[a.atom] = 1; break;
      case 'selected': for (const i of ctx.selection || []) if (i >= 0 && i < n) out[i] = 1; break;
      case 'helix': case 'sheet': case 'coil': {
        const ss = secondary(ctx), want = fl === 'helix' ? 1 : fl === 'sheet' ? 2 : 0;
        for (let i = 0; i < n; i++) if (ac[i] === C_PROTEIN && (ss ? ss[cls.resOf[i]] : 0) === want) out[i] = 1;
        break;
      }
    }
    return out;
  }
  /** Atoms within d Å of any atom of src (uniform grid over the source atoms). */
  function withinMask(ctx, d, src) {
    const n = ctx.n, xyz = ctx.xyz, out = new Uint8Array(n);
    let mnx = Infinity, mny = Infinity, mnz = Infinity, mxx = -Infinity, mxy = -Infinity, mxz = -Infinity, ns = 0;
    for (let i = 0; i < n; i++) {
      if (!src[i]) continue;
      ns++;
      const x = xyz[3 * i], y = xyz[3 * i + 1], z = xyz[3 * i + 2];
      if (x < mnx) mnx = x; if (x > mxx) mxx = x; if (y < mny) mny = y; if (y > mxy) mxy = y; if (z < mnz) mnz = z; if (z > mxz) mxz = z;
    }
    if (!ns) return out;
    if (!(d > 0)) { out.set(src); return out; }
    const cell = Math.max(d, 0.5), nx = Math.floor((mxx - mnx) / cell) + 1, ny = Math.floor((mxy - mny) / cell) + 1, nz = Math.floor((mxz - mnz) / cell) + 1;
    const dense = nx * ny * nz <= 4e6;
    const head = dense ? new Int32Array(nx * ny * nz).fill(-1) : null, map = dense ? null : new Map(), next = new Int32Array(n).fill(-1);
    const key = (ix, iy, iz) => (iz * ny + iy) * nx + ix;
    for (let i = 0; i < n; i++) {
      if (!src[i]) continue;
      const k = key(Math.floor((xyz[3 * i] - mnx) / cell), Math.floor((xyz[3 * i + 1] - mny) / cell), Math.floor((xyz[3 * i + 2] - mnz) / cell));
      if (dense) { next[i] = head[k]; head[k] = i; } else { next[i] = map.has(k) ? map.get(k) : -1; map.set(k, i); }
    }
    const d2 = d * d;
    for (let i = 0; i < n; i++) {
      if (src[i]) { out[i] = 1; continue; }
      const x = xyz[3 * i], y = xyz[3 * i + 1], z = xyz[3 * i + 2];
      if (x < mnx - d || x > mxx + d || y < mny - d || y > mxy + d || z < mnz - d || z > mxz + d) continue;
      const cx = Math.floor((x - mnx) / cell), cy = Math.floor((y - mny) / cell), cz = Math.floor((z - mnz) / cell);
      let found = false;
      for (let iz = cz - 1; iz <= cz + 1 && !found; iz++) {
        if (iz < 0 || iz >= nz) continue;
        for (let iy = cy - 1; iy <= cy + 1 && !found; iy++) {
          if (iy < 0 || iy >= ny) continue;
          for (let ix = cx - 1; ix <= cx + 1 && !found; ix++) {
            if (ix < 0 || ix >= nx) continue;
            const k = key(ix, iy, iz);
            for (let j = dense ? head[k] : (map.has(k) ? map.get(k) : -1); j >= 0; j = next[j]) {
              const dx = xyz[3 * j] - x, dy = xyz[3 * j + 1] - y, dz = xyz[3 * j + 2] - z;
              if (dx * dx + dy * dy + dz * dz <= d2) { found = true; break; }
            }
          }
        }
      }
      if (found) out[i] = 1;
    }
    return out;
  }
  function count(mask) { let c = 0; for (let i = 0; i < mask.length; i++) c += mask[i]; return c; }
  function indices(mask) { const out = []; for (let i = 0; i < mask.length; i++) if (mask[i]) out.push(i); return out; }
  /** {ok: true, count} or {ok: false, error, pos, end} for the UI. */
  function check(text, ctx) {
    try { const m = evaluate(String(text), ctx); return { ok: true, count: count(m), mask: m }; }
    catch (e) { if (e instanceof SelectionError) return { ok: false, error: e.message, pos: e.pos, end: e.end }; throw e; }
  }

  // ---- autocomplete
  /** Suggestions for the word at the caret: {from, to, items: [{text, hint}]}. After chain/resname/name/element/layer
   *  the values present in the structure are offered; elsewhere keywords and operators. */
  function complete(text, caret, ctx) {
    text = String(text || '');
    if (caret === undefined || caret === null) caret = text.length;
    let from = caret;
    while (from > 0 && !/[\s()"]/.test(text[from - 1])) from--;
    let to = caret;
    while (to < text.length && !/[\s()"]/.test(text[to])) to++;
    const word = text.slice(from, caret).toLowerCase();
    let toks = [];
    try { toks = tokenize(text.slice(0, from)); } catch (e) { toks = []; }
    let argKw = null;
    for (let q = toks.length - 1; q >= 0; q--) {
      const w = toks[q].kind === 'word' ? toks[q].v.toLowerCase() : null;
      if (toks[q].kind !== 'word' && toks[q].kind !== 'str') break;
      if (ARGS[w]) { argKw = ARGS[w]; break; }
      if (STOP.has(w) || FLAGS[w] || w === 'of' || w === 'as' || w === 'residue') break;
    }
    let items = [];
    const last = toks.length ? toks[toks.length - 1] : null, lw = last && last.kind === 'word' ? last.v.toLowerCase() : null;
    if (lw === 'within') items = ['3', '4', '5', '6', '8'].map(t => ({ text: t, hint: 'Å' }));
    else if (lw === 'same') items = [{ text: 'residue', hint: '' }];
    else if (lw === 'residue' && toks.length > 1 && toks[toks.length - 2].v.toLowerCase() === 'same') items = [{ text: 'as', hint: '' }];
    else if (argKw && ctx) items = valuesFor(argKw, ctx);
    else if (argKw === 'layer') items = ['H', 'M', 'L'].map(t => ({ text: t, hint: '' }));
    if (!items.length && !argKw && lw !== 'within') items = KEYWORDS.map(k => ({ text: k, hint: KW_HINT[k] || '' }));
    const out = items.filter(it => !word || it.text.toLowerCase().startsWith(word));
    return { from, to, items: out.slice(0, 40) };
  }
  function valuesFor(kw, ctx) {
    const f = ctx.frame, n = ctx.n;
    if (kw === 'layer') {
      const L = ctx.layers && ctx.layers.atomLayer;
      const lab = ctx.layers && ctx.layers.labels || {};
      return ['H', 'M', 'L'].filter(x => !L || L.includes(x)).map(x => ({ text: x, hint: lab[x] || '' }));
    }
    if (kw === 'element') {
      const c = new Map();
      for (let i = 0; i < n; i++) c.set(ctx.Z[i], (c.get(ctx.Z[i]) || 0) + 1);
      return [...c.entries()].sort((a, b) => b[1] - a[1]).map(([z, k]) => ({ text: el.symbol(z), hint: k + ' atoms' }));
    }
    const src = kw === 'name' ? f.atomNames : kw === 'resname' ? f.resNames : kw === 'chain' ? f.chains : null;
    if (!src) return [];
    const c = new Map();
    const step = n > 200000 ? 7 : 1;
    for (let i = 0; i < n; i += step) { const v = String(src[i] === null || src[i] === undefined ? '' : src[i]).trim(); if (v) c.set(v, (c.get(v) || 0) + 1); }
    return [...c.entries()].sort((a, b) => b[1] - a[1]).slice(0, 80).map(([v, k]) => ({ text: /\s/.test(v) ? '"' + v + '"' : v, hint: (step > 1 ? '~' + k * step : k) + (kw === 'resname' ? ' atoms' : '') }));
  }

  // ------------------------------------------------------------------ colour schemes
  /** Colours (Float32Array 3n) of every atom in a scheme ('element', 'layer', … or '#rrggbb'); `mask` limits the work. */
  function colorize(scheme, ctx, mask) {
    const n = ctx.n, Z = ctx.Z, out = new Float32Array(3 * n);
    const put = (i, c) => { out[3 * i] = c[0]; out[3 * i + 1] = c[1]; out[3 * i + 2] = c[2]; };
    const want = i => !mask || mask[i];
    for (let i = 0; i < n; i++) if (want(i)) put(i, elColor(Z[i]));
    if (!scheme || scheme === 'element') return out;
    if (/^#[0-9a-f]{6}$/i.test(scheme)) { const c = hex(scheme); for (let i = 0; i < n; i++) if (want(i)) put(i, c); return out; }
    const cls = classify(ctx);
    switch (scheme) {
      case 'layer': {
        const L = ctx.layers && ctx.layers.atomLayer;
        if (L) for (let i = 0; i < n; i++) if (want(i) && LAYER_COL[L[i]]) put(i, LAYER_COL[L[i]]);
        break;
      }
      case 'chain': {
        const ids = new Map();
        for (let r = 0; r < cls.nRes; r++) { const c = cls.chain[r] || ''; if (!ids.has(c)) ids.set(c, ids.size); }
        if (ids.size <= 1 && !ctx.frame.chains) break;
        for (let i = 0; i < n; i++) if (want(i)) put(i, CHAIN_COL[ids.get(cls.chain[cls.resOf[i]] || '') % CHAIN_COL.length]);
        break;
      }
      case 'residue': {
        const memo = new Map();
        for (let r = 0; r < cls.nRes; r++) {
          const nm = cls.resName[r];
          let c = memo.get(nm);
          if (!c) {
            const rc = cls.resClass[r];
            if (rc === C_PROTEIN) c = AA_COL[nm] ? hex(AA_COL[nm]) : AA_COL[nm.slice(-3)] ? hex(AA_COL[nm.slice(-3)]) : hashCol(nm);
            else if (rc === C_NUCLEIC) { const b = nm.replace(/^[DR]/, '').replace(/[35N]$/, '')[0]; c = NUC_COL[b] || RT_COL.nucleic; }
            else c = nm ? hashCol(nm) : null;
            memo.set(nm, c || false);
          }
          if (!c) continue;
          for (let t = cls.resOff[r]; t < cls.resOff[r + 1]; t++) { const i = cls.resAtoms[t]; if (want(i)) put(i, c); }
        }
        break;
      }
      case 'restype': {
        for (let i = 0; i < n; i++) {
          if (!want(i)) continue;
          const r = cls.resOf[i], rc = cls.resClass[r];
          const c = rc === C_PROTEIN ? RT_COL[RES_TYPE[cls.resName[r]] || RES_TYPE[cls.resName[r].slice(-3)] || 'polar'] : RT_COL[CLASS_IDS[rc]];
          if (c) put(i, c);
        }
        break;
      }
      case 'ss': {
        const ss = secondary(ctx), coil = ctx.lightBg ? SS_COL.coilLight : SS_COL.coil;
        for (let i = 0; i < n; i++) {
          if (!want(i)) continue;
          const r = cls.resOf[i];
          if (cls.resClass[r] !== C_PROTEIN) continue;
          const s = ss ? ss[r] : 0;
          put(i, s === 1 ? SS_COL.helix : s === 2 ? SS_COL.sheet : coil);
        }
        break;
      }
      case 'rainbow': {
        // per chain of polymer residues: blue (first) -> red (last)
        let r = 0;
        while (r < cls.nRes) {
          const rc = cls.resClass[r];
          if (rc !== C_PROTEIN && rc !== C_NUCLEIC) { r++; continue; }
          let e = r;
          while (e + 1 < cls.nRes && cls.resClass[e + 1] === rc && cls.chain[e + 1] === cls.chain[r]) e++;
          for (let q = r; q <= e; q++) {
            const c = hsv(240 * (1 - (e > r ? (q - r) / (e - r) : 0)), 0.75, 0.95);
            for (let t = cls.resOff[q]; t < cls.resOff[q + 1]; t++) { const i = cls.resAtoms[t]; if (want(i)) put(i, c); }
          }
          r = e + 1;
        }
        break;
      }
      case 'sugar': {
        sugars(ctx);
        for (const [r, s] of cls.sugarOfRes) for (let t = cls.resOff[r]; t < cls.resOff[r + 1]; t++) { const i = cls.resAtoms[t]; if (want(i)) put(i, s.color); }
        break;
      }
      case 'charge': {
        const q = ctx.charges;
        if (!q || q.length !== n) break;
        let m = 0.2;
        for (let i = 0; i < n; i++) if (Math.abs(q[i]) > m) m = Math.abs(q[i]);
        for (let i = 0; i < n; i++) {
          if (!want(i)) continue;
          const t = Math.max(-1, Math.min(1, q[i] / m));
          put(i, t < 0 ? [1 + 0.75 * t, 1 + 0.55 * t, 1] : [1, 1 - 0.7 * t, 1 - 0.75 * t]);
        }
        break;
      }
      case 'bfactor': case 'occupancy': {
        const v = scheme === 'bfactor' ? ctx.frame.bfactors : ctx.frame.occupancy;
        if (!v || v.length !== n) break;
        const s = Array.from(v).filter(Number.isFinite).sort((a, b) => a - b);
        if (!s.length) break;
        const lo = s[Math.floor(s.length * 0.02)], hi = s[Math.min(s.length - 1, Math.floor(s.length * 0.98))];
        for (let i = 0; i < n; i++) if (want(i) && Number.isFinite(v[i])) put(i, ramp(hi > lo ? Math.max(0, Math.min(1, (v[i] - lo) / (hi - lo))) : 0.5));
        break;
      }
    }
    return out;
  }

  // ------------------------------------------------------------------ representation rows and presets
  function row(name, selection, style, color, extra) {
    return Object.assign({ name, selection, style, color: color || 'element', opacity: 1, scale: 1, visible: true }, extra || {});
  }
  const PRESETS = [
    { id: 'all', label: 'All atoms', title: 'Every atom in the toolbar style (default for small molecules)' },
    { id: 'qmmm', label: 'QM/MM layers', title: 'High layer ball-and-stick, medium sticks, low lines; frozen atoms dimmed; link atoms marked; legend with the methods' },
    { id: 'protein', label: 'Protein–ligand', title: 'Cartoon by secondary structure, ligand ball-and-stick, residues within 5 Å of it as sticks, water hidden, ions as small spheres' },
    { id: 'membrane', label: 'Membrane', title: 'Protein cartoon, lipid tails as lines with P/N spheres, water hidden, ions as small spheres' },
    { id: 'glycan', label: 'Glycan / proteoglycan', title: 'Protein cartoon, sugars as sticks with filled rings and SNFG symbols' },
    { id: 'nucleic', label: 'Nucleic acid', title: 'Backbone tube with base slabs coloured by base, protein cartoon' }
  ];
  function presetRows(id, ctx, view) {
    view = view || {};
    const cls = ctx ? classify(ctx) : null, has = k => !!(cls && cls.counts[k].atoms);
    const tail = () => {
      const r = [];
      if (has('ions')) r.push(row('Ions', 'ions', 'spacefill', 'element', { scale: 0.5 }));
      if (has('water')) r.push(row('Water', 'water', 'hidden'));
      return r;
    };
    const ligand = () => (has('ligand') ? [row('Ligand', 'ligand', 'ballstick')] : []);
    switch (id) {
      case 'qmmm': {
        const L = ctx && ctx.layers, al = L && L.atomLayer, out = [];
        const any = x => !!(al && al.includes(x));
        const fr = L && L.frozen && L.frozen.some(Boolean);
        if (!al) return [row('All atoms', 'all', view.style || 'ballstick')];
        if (any('H')) out.push(row('QM (high)', 'layer H', 'ballstick'));
        if (any('M')) out.push(row('QM2 (medium)', 'layer M', 'sticks'));
        if (any('L')) out.push(row('MM (low)', 'layer L', 'lines'));
        if (fr && any('L')) out.push(row('Frozen MM', 'frozen and layer L', 'lines', 'element', { opacity: 0.35 }));
        if (fr && (any('H') || any('M')) && al.some((x, i) => x !== 'L' && L.frozen[i])) out.push(row('Frozen QM', 'frozen and not layer L', 'ghost', 'element', { opacity: 0.4 }));
        if (L.linkAtoms && L.linkAtoms.length) out.push(row('Link atoms', 'link', 'markers', '#ff7ad9'));
        return out;
      }
      case 'protein': return [row('Protein', 'protein', 'cartoon', 'ss')]
        .concat(has('nucleic') ? [row('Nucleic acid', 'nucleic', 'cartoon', 'residue')] : [])
        .concat(has('ligand') && has('protein') ? [row('Binding site', 'same residue as (protein and within 5 of ligand)', 'sticks')] : [])
        .concat(ligand(), has('carbohydrate') ? [row('Glycans', 'carbohydrate', 'sticks')] : [], has('lipid') ? [row('Lipids', 'lipid', 'lipid')] : [], tail());
      case 'membrane': return (has('protein') ? [row('Protein', 'protein', 'cartoon', 'ss')] : [])
        .concat([row('Lipids', 'lipid', 'lipid')], has('carbohydrate') ? [row('Glycans', 'carbohydrate', 'sticks')] : [], ligand(), tail());
      case 'glycan': return (has('protein') ? [row('Protein', 'protein', 'cartoon', 'ss')] : [])
        .concat(has('nucleic') ? [row('Nucleic acid', 'nucleic', 'cartoon', 'residue')] : [])
        .concat([row('Glycan', 'carbohydrate', 'sticks', 'element', { scale: 0.75 }), row('Sugar rings', 'carbohydrate', 'rings', 'sugar', { opacity: 0.85 }), row('SNFG', 'carbohydrate', 'snfg')],
          has('lipid') ? [row('Lipids', 'lipid', 'lipid')] : [], ligand(), tail());
      case 'nucleic': return [row('Nucleic acid', 'nucleic', 'cartoon', 'residue')]
        .concat(has('protein') ? [row('Protein', 'protein', 'cartoon', 'ss')] : [], ligand(), has('carbohydrate') ? [row('Glycans', 'carbohydrate', 'sticks')] : [], tail());
      default: return [row('All atoms', 'all', view.style || 'ballstick', view.color === 'charge' ? 'charge' : 'element', { base: true })];
    }
  }
  /** Preset for a newly opened structure: QM/MM layers, membrane, glycan, nucleic acid, protein or all atoms. */
  function autoPreset(ctx) {
    if (ctx.layers && ctx.layers.atomLayer && ctx.layers.atomLayer.length === ctx.n) return 'qmmm';
    const c = classify(ctx).counts;
    if (c.lipid.res >= 5) return 'membrane';
    if (c.carbohydrate.res >= 1 && (c.protein.res >= 1 || c.carbohydrate.res >= 3)) return 'glycan';
    if (c.nucleic.res >= 4) return 'nucleic';
    if (c.protein.res >= 20) return 'protein';
    return 'all';
  }
  const TOOLBAR_STYLES = new Set(['ballstick', 'sticks', 'spacefill', 'wireframe']);
  /** True for the plain one-row "All atoms" list that the renderer draws without per-atom data (toolbar style). */
  function isDefault(rows) {
    if (!rows || rows.length !== 1) return false;
    const r = rows[0];
    return !!r.base && r.visible !== false && String(r.selection).trim().toLowerCase() === 'all' && TOOLBAR_STYLES.has(r.style) &&
      !(r.opacity < 0.999) && (r.scale === undefined || r.scale === 1) && (r.color === 'element' || r.color === 'charge');
  }
  /** Rows -> {style: Uint8Array, scale, alpha: Float32Array, colors: Float32Array(3n), extent: Uint8Array,
   *  overlays: [{kind, mask, color, opacity, row, index}], markers: [{indices, color}], counts, errors}. */
  function resolve(rows, ctx) {
    const n = ctx.n, Z = ctx.Z;
    const style = new Uint8Array(n), scale = new Float32Array(n).fill(1), alpha = new Float32Array(n).fill(1), colors = new Float32Array(3 * n), extent = new Uint8Array(n);
    const ovExt = new Uint8Array(n); // atoms under an overlay count for fitting the view even when hidden
    const overlays = [], markers = [], counts = [], errors = [];
    let cls = null;
    rows.forEach((r, idx) => {
      counts.push(null); errors.push(null);
      let mask;
      try { mask = evaluate(String(r.selection), ctx); }
      catch (e) { if (e instanceof SelectionError) { errors[idx] = e.message; return; } throw e; }
      counts[idx] = count(mask);
      if (r.visible === false || !counts[idx]) return;
      const st = STYLE_BY_ID[r.style] || STYLE_BY_ID.ballstick;
      const op = r.opacity === undefined ? 1 : Math.max(0.05, Math.min(1, +r.opacity || 0));
      const sc = r.scale > 0 ? +r.scale : 1;
      if (st.kind === 'overlay') {
        if (r.style === 'markers') {
          const c = /^#[0-9a-f]{6}$/i.test(r.color) ? hex(r.color) : [1.0, 0.48, 0.85];
          markers.push({ indices: indices(mask), color: c });
        } else overlays.push({ kind: r.style, mask, color: r.color, opacity: op, row: r, index: idx });
        for (let i = 0; i < n; i++) if (mask[i]) ovExt[i] = 1;
        return;
      }
      const cols = colorize(r.color, ctx, mask);
      if (r.style === 'lipid') {
        cls = cls || classify(ctx);
        for (let i = 0; i < n; i++) {
          if (!mask[i]) continue;
          const z = Z[i], rn = cls.resName[cls.resOf[i]];
          let code = CODE.lines, s = sc;
          if (z === 1) code = CODE.hidden;
          else if (z === 15 || (z === 7 && cls.atomClass[i] === C_LIPID) || (z === 8 && STEROL.has(rn))) { code = CODE.spacefill; s = sc * (z === 15 ? 0.55 : 0.45); }
          style[i] = code; scale[i] = s; alpha[i] = op;
          colors[3 * i] = cols[3 * i]; colors[3 * i + 1] = cols[3 * i + 1]; colors[3 * i + 2] = cols[3 * i + 2];
        }
        return;
      }
      const code = CODE[r.style] === undefined ? CODE.ballstick : CODE[r.style];
      for (let i = 0; i < n; i++) {
        if (!mask[i]) continue;
        style[i] = code; scale[i] = sc; alpha[i] = op;
        colors[3 * i] = cols[3 * i]; colors[3 * i + 1] = cols[3 * i + 1]; colors[3 * i + 2] = cols[3 * i + 2];
      }
    });
    for (let i = 0; i < n; i++) extent[i] = style[i] || ovExt[i] ? 1 : 0;
    // atoms no row coloured keep element colours (a later overlay or the selection halo may still show them)
    for (let i = 0; i < n; i++) if (!style[i] && colors[3 * i] === 0 && colors[3 * i + 1] === 0 && colors[3 * i + 2] === 0) { const c = elColor(Z[i]); colors[3 * i] = c[0]; colors[3 * i + 1] = c[1]; colors[3 * i + 2] = c[2]; }
    return { style, scale, alpha, colors, extent, overlays, markers, counts, errors };
  }

  // ------------------------------------------------------------------ summaries
  /** Residue-class counts {protein: {res, atoms}, …}, chain count and QM/MM layer counts for the Info tab. */
  function contents(ctx) {
    const cls = classify(ctx), chains = new Set();
    for (let r = 0; r < cls.nRes; r++) if (cls.resClass[r] === C_PROTEIN || cls.resClass[r] === C_NUCLEIC) chains.add(cls.chain[r]);
    const layers = {};
    const L = ctx.layers && ctx.layers.atomLayer;
    if (L) for (const x of L) layers[x] = (layers[x] || 0) + 1;
    return { classes: cls.counts, chains: chains.size, layers, nRes: cls.nRes };
  }
  /** Short description of atom i's residue for the status bar: "GLY 45 A · helix", "BGLC 12 · GlcNS6S". */
  function describe(ctx, i) {
    if (!(i >= 0 && i < ctx.n)) return '';
    const cls = classify(ctx), r = cls.resOf[i], rc = cls.resClass[r];
    const parts = [];
    if (ctx.frame.resNames) parts.push((cls.resName[r] + ' ' + (cls.resId[r] === undefined ? '' : cls.resId[r]) + (cls.chain[r] ? ' ' + cls.chain[r] : '')).trim());
    if (rc === C_SUGAR) { sugars(ctx); const s = cls.sugarOfRes.get(r); parts.push(s ? s.label : 'carbohydrate'); }
    else if (rc === C_PROTEIN && cls.ss) parts.push(['coil', 'helix', 'sheet'][cls.ss[r]] || 'coil');
    else if (rc !== C_PROTEIN && (ctx.frame.resNames || rc !== C_OTHER)) parts.push(CLASS_IDS[rc]);
    if (ctx.layers && ctx.layers.atomLayer && ctx.layers.atomLayer[i]) parts.push('layer ' + ctx.layers.atomLayer[i]);
    return parts.filter(Boolean).join(' · ');
  }

  BM.reps = {
    CODE, STYLES, COLORS, PRESETS, KEYWORDS, KW_HINT, SN, LAYER_COL, SS_COL, CLASS_IDS, CHAIN_COL,
    SelectionError, context, adjacency, classify, sugars, residueRings, secondary, parse, evaluate, check, count, indices, complete,
    colorize, row, presetRows, autoPreset, isDefault, resolve, contents, describe, isMetal, toHex, hex,
    tables: { AMINO, NUCLEIC, WATER, IONS, LIPID, STEROL, CCD }
  };
})(typeof window !== 'undefined' ? (window.BM = window.BM || {}) : (globalThis.BM = globalThis.BM || {}));
