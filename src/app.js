/* BlueMolecule — application logic. */
(function (BM) {
  'use strict';
  const $ = (s, r) => (r || document).querySelector(s);
  const $$ = (s, r) => Array.from((r || document).querySelectorAll(s));
  const el = BM.el, geom = BM.geom;
  const HARTREE = { kcal: 627.509474, kj: 2625.49964, ev: 27.211386 };
  const UNIT_LABEL = { kcal: 'kcal/mol', kj: 'kJ/mol', ev: 'eV' };
  const BG = { navy: [0.063, 0.090, 0.133], black: [0, 0, 0], white: [1, 1, 1], gray: [0.91, 0.925, 0.945] };
  const BG2 = { navy: [0.16, 0.22, 0.32], black: [0.16, 0.16, 0.19], white: [0.92, 0.94, 0.97], gray: [1, 1, 1] };   // the top of a background gradient
  // scene presets (issue #15): background, material and effects together
  const SCENES = {
    default: { bg: 'navy', material: 'standard', ao: 0, outline: 0, dof: 0, shadows: false, floor: false, gradient: false, fog: true },
    publication: { bg: 'white', material: 'matte', ao: 0.7, outline: 1, dof: 0, shadows: false, floor: false, gradient: false, fog: false },
    studio: { bg: 'gray', material: 'glossy', ao: 0.6, outline: 0, dof: 0, shadows: true, floor: true, gradient: true, fog: false },
    dark: { bg: 'black', material: 'glossy', ao: 0.6, outline: 0, dof: 0, shadows: false, floor: false, gradient: true, fog: true },
    cartoon: { bg: 'white', material: 'flat', ao: 0, outline: 2, dof: 0, shadows: false, floor: false, gradient: false, fog: false }
  };
  const FX_KEYS = ['material', 'ao', 'outline', 'dof', 'shadows', 'floor', 'gradient'];
  const VERSION = '2.6.2';
  const esc = s => String(s === undefined || s === null ? '' : s).replace(/[&<>"]/g, c => ({ '&': '&amp;', '<': '&lt;', '>': '&gt;', '"': '&quot;' }[c]));
  const isNum = v => typeof v === 'number' && isFinite(v);
  const fx = (v, d) => (isNum(v) ? v.toFixed(d) : '—');
  const kit = BM.buildUI.kit;
  BM.version = VERSION;
  const PERCEIVE_MAX = kit.LIMIT.perceive;

  const state = {
    docs: [], active: null, nextId: 1,
    view: { style: 'ballstick', labels: 'none', color: 'element', showH: true, bondLengths: false, perspective: true, fog: true, axes: true, bg: 'navy', atomScale: 0.25, bondRadius: 0.11, vdwScale: 1, tol: 0.4, dummies: true, theme: 'dark', bondOrders: true, hbonds: false, dipole: false, cell: true,
      // scene (issue #15): a preset, or 'custom' once an effect is changed by hand
      scene: 'default', material: 'standard', ao: 0, outline: 0, dof: 0, shadows: false, floor: false, gradient: false, fxLarge: false },
    anim: { kind: null, raf: 0, last: 0, phase: 0, acc: 0 },
    // amplitude: largest atom shift of the animation and the two ends (Å); arrowLen: length of the longest arrow (Å);
    // ends: show both ends of the selected vibration as translucent structures
    vib: { amplitude: 0.35, speed: 1.0, arrows: false, arrowLen: 0.8, ends: false, playing: true, disp: 0.15 },
    fps: 6,
    plotCollapsed: false,
    build: {
      tool: 'navigate', element: 6, order: 1, adjustH: true, autoOpt: false, symTol: 0.1,
      tpl: { mode: 'molecules', molecule: 'benzene', group: 'methyl', ligand: 'NH3', metalZ: 78, geometry: 'square-planar', metalLigand: 'Cl' },
      fold: { align: true, clean: true, sym: true, coords: true, conf: true, crystal: true }, inputs: { el: '', smiles: '', libq: '', expand: '3' },
      conf: { max: 20, window: 10, rmsd: 0.5 },
      crys: { n: [2, 2, 2], hkl: [1, 0, 0], layers: 4, vacuum: 15, shift: 0, margin: 10, scale: true }
    }
  };
  let R = null, plot = null, ui = null;

  // ---------------- preferences ----------------
  function loadPrefs() {
    try {
      const p = JSON.parse(localStorage.getItem('bluemolecule.prefs') || '{}');
      if (p.view) Object.assign(state.view, p.view);
      if (p.vib) for (const k of ['amplitude', 'speed', 'arrows', 'arrowLen']) if (p.vib[k] !== undefined) state.vib[k] = p.vib[k];
      if (p.build) {
        const b = state.build;
        for (const k of ['element', 'order', 'adjustH', 'autoOpt', 'symTol']) if (p.build[k] !== undefined) b[k] = p.build[k];
        if (p.build.tpl) Object.assign(b.tpl, p.build.tpl);
        if (p.build.fold) Object.assign(b.fold, p.build.fold);
        // (numbers only: a damaged preference must not reach the tools)
        const nums = (dst, src) => { if (src && typeof src === 'object') for (const k of Object.keys(dst)) { const v = src[k]; if (Array.isArray(dst[k]) ? Array.isArray(v) && v.length === dst[k].length && v.every(Number.isFinite) : typeof dst[k] === 'boolean' ? typeof v === 'boolean' : Number.isFinite(v)) dst[k] = Array.isArray(v) ? v.slice() : v; } };
        nums(b.conf, p.build.conf); nums(b.crys, p.build.crys);
      }
      if (p.sidebar) document.documentElement.style.setProperty('--sidebar-w', p.sidebar);
      if (p.plotH) document.documentElement.style.setProperty('--plot-h', p.plotH);
    } catch (e) { /* storage unavailable */ }
  }
  function savePrefs() {
    try {
      const cs = document.documentElement.style, b = state.build;
      localStorage.setItem('bluemolecule.prefs', JSON.stringify({ view: state.view, vib: { amplitude: state.vib.amplitude, speed: state.vib.speed, arrows: state.vib.arrows, arrowLen: state.vib.arrowLen },
        build: { element: b.element, order: b.order, adjustH: b.adjustH, autoOpt: b.autoOpt, symTol: b.symTol, tpl: b.tpl, fold: b.fold, conf: b.conf, crys: b.crys },
        sidebar: cs.getPropertyValue('--sidebar-w'), plotH: cs.getPropertyValue('--plot-h') }));
    } catch (e) { /* ignore */ }
  }

  // ---------------- messages ----------------
  let toastTimer = 0, hintTimer = 0, hint = '';
  function toast(msg, isError) {
    const t = $('#toast');
    t.textContent = msg; t.classList.toggle('error', !!isError); t.classList.remove('hidden');
    clearTimeout(toastTimer); toastTimer = setTimeout(() => t.classList.add('hidden'), isError ? 6000 : 2600);
    status(msg);
  }
  /** Transient message in the status bar; the active tool's hint comes back after a few seconds. */
  /** Copy to the clipboard (async API, then the selection fallback); the text stays visible either way. */
  function copyText(text, okMsg) {
    const fallback = () => {
      try {
        const ta = document.createElement('textarea');
        ta.value = text; ta.setAttribute('readonly', ''); ta.style.position = 'fixed'; ta.style.left = '-9999px';
        document.body.appendChild(ta); ta.select();
        const ok = document.execCommand('copy');
        ta.remove();
        if (ok) { toast(okMsg || 'Copied to clipboard'); return; }
      } catch (e) { /* ignore */ }
      const blocked = () => toast('Copying is blocked here; select the text and press Ctrl/⌘+C.');
      const nat = BM.native;
      if (nat && typeof nat.copyText === 'function') nat.copyText(text).then(ok => (ok ? toast(okMsg || 'Copied to clipboard') : blocked()), blocked);
      else blocked();
    };
    try {
      if (navigator.clipboard && navigator.clipboard.writeText) navigator.clipboard.writeText(text).then(() => toast(okMsg || 'Copied to clipboard'), fallback);
      else fallback();
    } catch (e) { fallback(); }
  }
  function status(msg) {
    $('#status-msg').textContent = msg;
    clearTimeout(hintTimer);
    hintTimer = setTimeout(() => { hintTimer = 0; showHint(); }, 7000);
  }
  /** The active tool's hint; `now` shows it at once (tool or option changed), otherwise after any message. */
  function setHint(h, now) { hint = h || ''; if (now) { clearTimeout(hintTimer); hintTimer = 0; } if (!hintTimer) showHint(); }
  function showHint() { $('#status-msg').textContent = hint || 'Ready'; }

  // ---------------- documents ----------------
  function firstChargeScheme(data) { return data.charges ? Object.keys(data.charges).find(k => Array.isArray(data.charges[k])) || null : null; }

  function makeSurface(v, i) {
    const st = BM.iso.stats(v);
    if (v.kind === 'molsurf') {   // a molecular surface: the zero level of its field (BM.molsurf)
      return { vol: v, stats: st, iso: 0, visible: i === 0, opacity: 1, both: false, colorPos: '#b9c9de', colorNeg: '#ff8c42', mesh: null, meshIso: null, esp: { src: 'none', range: null } };
    }
    let iso;
    if (v.kind === 'density') iso = st.absMax > 0.02 ? 0.002 : st.absMax * 0.2;
    else if (v.kind === 'spin') iso = st.absMax > 0.02 ? 0.002 : st.absMax * 0.2;
    else if (v.kind === 'orbital') iso = st.absMax > 0.06 ? 0.02 : st.absMax * 0.3;
    else iso = st.absMax * 0.1;
    iso = Number(iso.toPrecision(2)) || 0.001;
    // both signs: orbitals and other signed volumes; the spin density whenever it has a negative part (spin polarisation)
    const both = v.kind === 'spin' ? st.min < 0 : st.min < -iso * 0.5 && v.kind !== 'density';
    return { vol: v, stats: st, iso, visible: i === 0, opacity: 1, both, colorPos: '#3b7bff', colorNeg: '#ff8c42', mesh: null, meshIso: null, esp: { src: 'none', range: null } };
  }

  function makeDoc(name, data) {
    const d = {
      id: state.nextId++, name, data,
      frameIndex: Math.max(0, data.frames.length - 1),
      selection: [], vibIndex: -1, stateIndex: -1, bonds: new Map(), undo: [], redo: [], frozen: [],
      tab: 'info', plotType: null, convergedOnly: false, built: data.format === 'built', job: -1,
      reps: null, repPreset: null, repSel: 0, // representation rows (Display tab), set on first display
      pending: !!(data.meta && data.meta.topologyOnly), moSel: null,
      chargeScheme: firstChargeScheme(data),
      spec: {
        ir: { shape: 'lorentzian', fwhm: 12, scale: 1, yMode: 'relative', invert: false, peaks: true, sticks: true, xmin: null, xmax: null },
        raman: { shape: 'lorentzian', fwhm: 12, scale: 1, yMode: 'activity', laser: 532, temperature: 298.15, peaks: true, sticks: true, xmin: null, xmax: null },
        uv: { unit: 'nm', sigma: 0.4, peaks: true, sticks: true, xmin: null, xmax: null },
        ecd: { unit: 'nm', sigma: 0.3, peaks: false, sticks: true, xmin: null, xmax: null },
        nmr: { z: null, ref: 0, merge: true, peaks: true },
        mo: { window: 10 },
        opt: { log: true },
        energy: { unit: 'kcal', ref: 'min' }
      },
      surfaces: (data.volumes || []).map(makeSurface)
    };
    if (data.vibrations && data.vibrations.frameIndex < data.frames.length) d.frameIndex = data.vibrations.frameIndex;
    const v = data.vibrations;
    if (v && v.ir && v.freqs.some(f => f > 0)) d.plotType = 'ir';
    else if (v && v.raman) d.plotType = 'raman';
    else if (data.excitedStates && data.excitedStates.length) d.plotType = 'uv';
    else if (data.nmr && data.nmr.iso) d.plotType = 'nmr';
    if (d.pending) d.tab = 'info';
    else if (data.vibrations) d.tab = 'vib';
    else if (data.excitedStates && data.excitedStates.length) d.tab = 'states';
    else if (d.surfaces.length) d.tab = 'surf';
    else if (data.frames.length > 1) d.tab = 'steps';
    else if (d.built) d.tab = 'build';
    return d;
  }

  /** An empty document for building (one edited frame with explicit topology). */
  function addBuiltDoc(name) {
    if (!name) {
      const taken = new Set(state.docs.map(d => d.name));
      name = 'New molecule';
      for (let k = 2; taken.has(name); k++) name = 'New molecule ' + k;
    }
    const data = {
      format: 'built', program: 'BlueMolecule', title: name,
      meta: { route: null, keywords: null, method: null, basis: null, charge: null, multiplicity: null, jobTypes: [], normalTermination: null, solvent: null },
      frames: [{ Z: [], xyz: new Float64Array(0), energy: null, energies: null, label: 'Built', converged: false, opt: null, atomNames: null, resNames: null, resIds: null, edited: true, bonds: [], formalCharges: [] }],
      vibrations: null, excitedStates: null, charges: null, spins: null, dipole: null, orbitals: null, thermo: null, volumes: null, bonds: null, cell: null, warnings: []
    };
    const d = makeDoc(name, data);
    state.docs.push(d);
    activate(d);
    R.setOrientation([0, 0, 0, 1]);
    ui.fitForBuilding();
    status('New molecule: click in the view to place atoms, or insert SMILES or a template from the Build tab');
    return d;
  }

  function addDoc(name, data, note) {
    if (!data.frames || !data.frames.length) throw new Error('No atoms found in ' + name + '.');
    const d = makeDoc(name, data);
    state.docs.push(d);
    activate(d);
    const w = (data.warnings || []).length;
    const nf = data.frames.length;
    if (d.pending) { toast(`Opened ${name}: ${data.frames[0].Z.length} atoms, coordinates pending — open a coordinate or trajectory file next`); return d; }
    // frequencies with imaginary modes (or of a transition-state search): say how many, and show the largest one
    // unless it is only numerical noise (issue #14)
    let vibNote = '';
    const v = data.vibrations;
    if (v && v.frameIndex === d.frameIndex) {
      const st = BM.vibtools.stationary(v), ts = tsJob(d);
      if (st.imag.length || ts) vibNote = ' — ' + BM.vibtools.describe(v, st, ts).short;
      if (st.imag.length && (ts || st.large.length) && modeOf(d, st.imag[0])) { selectVib(d, st.imag[0]); renderPanel(); vibNote += ', selected'; }
    }
    toast(`Opened ${name}${note || ''}: ${data.frames[d.frameIndex].Z.length} atoms${nf > 1 ? `, ${nf} geometries` : ''}${v ? `, ${v.freqs.length} vibrations` : ''}${w ? ` (${w} warning${w > 1 ? 's' : ''})` : ''}${vibNote}`);
    return d;
  }

  function activate(d) {
    stopAnim(true);
    state.active = d;
    if (d && d.reps && d.reps.length === 1 && d.reps[0].base) { d.reps[0].style = state.view.style; d.reps[0].color = state.view.color === 'charge' ? 'charge' : 'element'; }
    $('#empty').classList.toggle('hidden', !!d);
    document.title = d ? `${d.name} — BlueMolecule` : 'BlueMolecule';
    if (ui) ui.onActivate(d);
    refreshScene(false);
    renderDocs(); renderTabs(); renderPanel(); renderPlot(); renderPlaybar(); updateStatus();
  }

  function closeDoc(d) {
    const i = state.docs.indexOf(d);
    if (i < 0) return;
    state.docs.splice(i, 1);
    if (state.active === d) activate(state.docs[Math.min(i, state.docs.length - 1)] || null);
    else renderDocs();
    if (!state.docs.length) { $('#empty').classList.remove('hidden'); R.setScene([], new Float64Array(0), new Uint32Array(0)); R.setSurfaces([]); refreshExtras(); }
  }

  const curFrame = d => d.data.frames[d.frameIndex];
  /** The unit cell of a frame: an edited frame's own (null once removed), otherwise the frame's or the file's. */
  const cellOf = (d, f) => { f = f || curFrame(d); return f && f.edited && f.cell !== undefined ? f.cell : ((f && f.cell) || d.data.cell || null); };

  /** Topology of frame fi, cached per frame: one stable flat pairs array plus orders.
   *  Edited frames use their explicit bonds; frame 0 of MOL/SDF files the file's bonds; others distance bonds
   *  whose orders are perceived lazily (up to PERCEIVE_MAX atoms). Invalidated by invalidateTopo(). */
  function topoFor(d, fi) {
    const f = d.data.frames[fi];
    const fileBonds = !f.bonds && d.data.bonds && d.data.bonds.length && !f.edited && (fi === 0 || d.topoAll) && Math.abs(state.view.tol - 0.4) < 1e-9;
    const key = fi + '|' + (f.bonds ? 'x' : fileBonds ? 'f' : state.view.tol);
    let t = d.bonds.get(key);
    if (t) return t;
    const list = f.bonds || (fileBonds ? d.data.bonds : null);
    if (list) {
      const pairs = new Uint32Array(list.length * 2), orders = new Float32Array(list.length);
      list.forEach((b, k) => { pairs[2 * k] = b[0]; pairs[2 * k + 1] = b[1]; orders[k] = b[2] === undefined ? 1 : +b[2] || 1; });
      t = { pairs, orders, perceived: false, charges: f.bonds ? f.formalCharges || null : null };
    } else t = { pairs: geom.computeBonds(f.Z, f.xyz, { tolerance: state.view.tol }), orders: undefined, perceived: true, charges: null };
    d.bonds.set(key, t);
    return t;
  }
  function bondsFor(d, fi) { return topoFor(d, fi).pairs; }
  /** Bond orders parallel to bondsFor (perceived on demand), or null when unknown (all drawn single). */
  function ordersFor(d, fi) {
    const t = topoFor(d, fi);
    if (t.orders === undefined) {
      t.orders = null;
      const f = d.data.frames[fi];
      if (f.Z.length <= PERCEIVE_MAX && t.pairs.length) {
        try {
          const p = kit.perceive(f.Z, f.xyz, t.pairs, null);
          t.orders = Float32Array.from(p.bonds, b => b[2]); t.charges = p.charges;
        } catch (e) { console.warn('[BlueMolecule] bond-order perception failed:', e); }
      }
    }
    return t.orders;
  }
  function invalidateTopo(d, fi) { for (const k of [...d.bonds.keys()]) if (k.startsWith(fi + '|')) d.bonds.delete(k); }
  /** {Z, xyz, bonds: [[i, j, order]], formalCharges} of the current frame (explicit, or perceived for file frames). */
  function topoFrame(d) {
    const fi = d.frameIndex, f = d.data.frames[fi], t = topoFor(d, fi), orders = ordersFor(d, fi);
    return { Z: f.Z, xyz: f.xyz, bonds: f.bonds || kit.pairsToList(t.pairs, orders), formalCharges: f.formalCharges || t.charges || new Array(f.Z.length).fill(0) };
  }

  // ---------------- file opening ----------------
  const FORMAT_SUMMARY = 'QC outputs and inputs (Gaussian, ORCA, GAMESS, NWChem, Q-Chem, Psi4, Molpro, Dalton, xtb, MOPAC, CP2K, Turbomole), fchk, Molden, cube; ' +
    'XYZ, MOL/SDF, PDB, mmCIF, CIF, Mol2, CML, CJSON, VASP, Quantum ESPRESSO, CASTEP, XSF, GRO, LAMMPS; prmtop/PSF with inpcrd, Amber NetCDF, DCD, TRR, XTC; .gz';
  const conciseError = msg => String(msg).replace(/Supported: [\s\S]*$/, 'Supported: ' + FORMAT_SUMMARY + '.');
  // Files of BIG bytes or more (and .gz files and streamed launcher files) are read with BM.bigfile.readLines and
  // parsed from lines (BM.parseLines): V8 cannot hold a string above ~512 MiB. Smaller files take the fast path.
  const BIGF = BM.bigfile || null;
  const BIG0 = BIGF ? BIGF.BIG : 64 * 1048576;
  const bigSize = () => (window.__bmBigSize > 0 ? window.__bmBigSize : BIG0);   // (tests lower it to read small files the big way)
  const MBs = b => (b >= 1e9 ? (b / 1073741824).toFixed(2) + ' GB' : (b / 1048576).toFixed(b < 10 * 1048576 ? 1 : 0) + ' MB');
  const nextPaint = () => new Promise(r => (typeof requestAnimationFrame === 'function' ? requestAnimationFrame(() => setTimeout(r, 0)) : setTimeout(r, 0)));
  let progEl = null;
  /** Slim progress card with Cancel (top centre); begin(label, size) -> {signal, onProgress, set, parsing, done}. */
  function progressCard() {
    if (progEl) return progEl;
    const el = document.createElement('div');
    el.id = 'bm-progress';
    el.setAttribute('role', 'progressbar'); el.setAttribute('aria-valuemin', '0'); el.setAttribute('aria-valuemax', '100');
    el.style.cssText = 'position:fixed;left:50%;top:54px;transform:translateX(-50%);z-index:130;width:min(460px,86vw);box-sizing:border-box;' +
      'background:var(--panel);color:var(--text);border:1px solid var(--line-2);border-radius:10px;padding:8px 12px 10px;' +
      'box-shadow:0 12px 30px rgba(0,0,0,.35);font-size:12px;display:none';
    el.innerHTML = '<div style="display:flex;align-items:center;gap:10px"><span class="bmp-t mono" style="flex:1;min-width:0;overflow:hidden;' +
      'text-overflow:ellipsis;white-space:nowrap"></span><button type="button" class="btn small bmp-x">Cancel</button></div>' +
      '<div style="height:4px;margin-top:7px;border-radius:2px;background:var(--line);overflow:hidden"><div class="bmp-f" ' +
      'style="height:100%;width:0;background:var(--accent);transition:width .15s"></div></div>';
    document.body.appendChild(el);
    progEl = el;
    return el;
  }
  function beginProgress(label, size, delay) {
    const el = progressCard(), ctl = new AbortController();
    const t = $('.bmp-t', el), f = $('.bmp-f', el), x = $('.bmp-x', el);
    let timer = 0, closed = false;
    x.disabled = false;
    x.onclick = () => { ctl.abort(); x.disabled = true; t.textContent = 'Cancelling…'; };
    const show = () => { if (!closed) el.style.display = ''; };
    if (delay) timer = setTimeout(show, delay); else show();
    const set = (frac, text) => {
      if (frac !== null && frac !== undefined && isFinite(frac)) { const v = Math.max(0, Math.min(1, frac)) * 100; f.style.width = v.toFixed(1) + '%'; el.setAttribute('aria-valuenow', v.toFixed(0)); }
      if (text) { t.textContent = text; status(text); }
    };
    set(0, 'Reading ' + label + '…');
    return {
      signal: ctl.signal, set,
      onProgress: p => { const tot = p.total || size; set(tot ? p.bytes / tot : null, `Reading ${label} — ${MBs(p.bytes)}${tot ? ' of ' + MBs(tot) : ''}`); },
      parsing: async text => { set(1, text); x.disabled = true; show(); await nextPaint(); },   // parsing is synchronous: no Cancel
      done: () => { closed = true; clearTimeout(timer); el.style.display = 'none'; }
    };
  }
  /** Whatever the shells and tests hand over (File/Blob, {name, size, readChunk}, {name, bytes}, stream) as a source. */
  function asSource(x, name) {
    if (BIGF) return BIGF.toSource(x, name);
    return { name: name || x.name || 'untitled', size: x.size, kind: 'blob', data: x };
  }
  async function readBytes(src, gz, prog) {
    if (gz) {
      if (typeof DecompressionStream === 'undefined') throw new Error('This browser cannot decompress .gz files; please unzip it first.');
      const st = BIGF ? BIGF.gunzip(src.data) : src.data.stream().pipeThrough(new DecompressionStream('gzip'));
      return new Uint8Array(await new Response(st).arrayBuffer());
    }
    if (src.kind === 'bytes') return src.data;
    if (src.kind === 'text') return new TextEncoder().encode(src.data);
    if (src.kind === 'blob' && !prog) return new Uint8Array(await src.data.arrayBuffer());
    return BIGF.readBytes(src.data, { onProgress: prog && prog.onProgress, signal: prog && prog.signal, total: src.size });
  }
  /** A memory estimate that cannot fit the page's JavaScript heap stops before reading (Chromium: jsHeapSizeLimit). */
  function heapLimit() {
    try { const m = performance.memory; return m && m.jsHeapSizeLimit > 0 ? m.jsHeapSizeLimit : 0; } catch (e) { return 0; }
  }
  /** Memory the computer can give to the number arrays of a big file (bytes): what the Linux app reports as free, else
   *  half of the browser's hint (navigator.deviceMemory: rounded, at most 8 GB), else unknown (no limit). */
  function ramLimit() {
    const q = window.BM_QT || window.BM_TAURI;   // the Linux Qt window or the desktop app
    if (q && q.memAvailable > 0) return 0.8 * q.memAvailable;
    const dm = typeof navigator !== 'undefined' ? navigator.deviceMemory : 0;
    return dm > 0 ? 0.5 * dm * 1073741824 : undefined;
  }
  /** Large text file: streamed line by line. Huge program outputs and XYZ trajectories go through a condenser that keeps
   *  only what their reader needs (condense.js, issue #6), with the choices of the loading dialog (load); other formats
   *  are collected into lines. */
  async function readLinesAndParse(src, name, gz, prog, load) {
    const lim = heapLimit();
    const opts = { signal: prog.signal, total: src.size, memoryLimit: lim ? 0.85 * lim : undefined, ramLimit: ramLimit() };
    let input = src.data;
    if (gz) {
      if (typeof DecompressionStream === 'undefined') throw new Error('This browser cannot decompress .gz files; please unzip it first.');
      opts.memoryLimit = undefined;            // the uncompressed size is not known in advance
      input = BIGF.gunzip(src.data, n => prog.onProgress({ bytes: n, total: src.size }));   // progress on the compressed bytes
    } else opts.onProgress = prog.onProgress;
    const pick = head => {
      const id = BM.detectFormat(name, head);
      return BM.condense && BM.condense.formats.indexOf(id) >= 0 ? BM.condense.create(id, { load }) : null;
    };
    let lines = BIGF.readCondensed ? await BIGF.readCondensed(input, pick, opts) : await BIGF.readLines(input, opts);
    input = null;
    const nLines = lines.length, condensed = !!lines.bmSide;
    await prog.parsing(`Parsing ${src.name} — ${nLines.toLocaleString()} ${condensed ? 'condensed ' : ''}lines…`);
    const t0 = performance.now();
    const ds = BM.parseLines(name, lines);
    console.info(`[BlueMolecule] parsed ${name} from ${nLines} ${condensed ? 'condensed ' : ''}lines in ${(performance.now() - t0).toFixed(0)} ms`);
    lines.length = 0; lines = null;
    return ds;
  }

  // ---------------- loading part of a large file (issue #7) ----------------
  /** From this size the loading dialog asks what to load (tests lower it with window.__bmAskSize). */
  const askSize = () => (window.__bmAskSize > 0 ? window.__bmAskSize : 1024 * 1048576);
  const ordinal = n => n + (n % 100 >= 11 && n % 100 <= 13 ? 'th' : ({ 1: 'st', 2: 'nd', 3: 'rd' }[n % 10] || 'th'));
  /** What the start and the end of a big output say (a few MB each): atoms, jobs, the number of steps, frequencies. */
  async function probeText(src, name) {
    const R = BIGF && BIGF.byteReader(src.data);
    if (!R || !R.size) return null;
    const dec = u8 => new TextDecoder('utf-8').decode(u8), MB4 = 4 * 1048576;
    const head = dec(await R.read(0, Math.min(R.size, MB4)));
    const tail = R.size > 2 * MB4 ? dec(await R.read(R.size - MB4, MB4)) : '';
    const format = BM.detectFormat(name, head.slice(0, 65536));
    if (!BM.condense || BM.condense.formats.indexOf(format) < 0) return null;
    const info = { format, size: R.size, atoms: null, jobs: [], steps: null, freq: false, qmmm: false, frames: null, bytesPerValue: format === 'gaussian-log' ? 1 : 4 };
    const lastOf = (re, text) => { let m, last = null; re.lastIndex = 0; while ((m = re.exec(text))) last = m; return last; };
    /** rows of the first table after a header line (Gaussian: after the second dashes line; ORCA: after the dashes) */
    const tableRows = (text, header, skipDashes, rowRe) => {
      const k = text.indexOf(header);
      if (k < 0) return null;
      const L = text.slice(k, k + 8 * 1048576).split(/\r?\n/);
      let j = 1, d = 0;
      while (j < L.length && d < skipDashes) { if (/^\s*-{10,}\s*$/.test(L[j])) d++; j++; }
      let n = 0;
      for (; j < L.length && rowRe.test(L[j]); j++) n++;
      return j < L.length ? n : null;   // (a table cut by the end of the sample says nothing)
    };
    if (format === 'gaussian-log') {
      const route = (/^\s*#.*$/m.exec(head) || [''])[0];
      for (const [re, j] of [[/\bopt\b/i, 'opt'], [/\bfreq\b/i, 'freq'], [/\bscan\b|modredundant/i, 'scan'], [/\birc\b/i, 'irc']]) if (re.test(route)) info.jobs.push(j);
      info.freq = info.jobs.indexOf('freq') >= 0 || /Harmonic frequencies/.test(tail);
      info.qmmm = /oniom\s*\(/i.test(route);
      const na = /NAtoms=\s*(\d+)/.exec(head);
      info.atoms = na ? +na[1] : tableRows(head, ' orientation:', 2, /^\s*\d+\s+-?\d+\s+/);
      const st = lastOf(/Step number\s+(\d+)\s+out of a maximum of/g, tail) || lastOf(/Step number\s+(\d+)\s+out of a maximum of/g, head);
      if (st) info.steps = +st[1];
    } else if (format === 'orca-out') {
      const kw = [];
      for (const m of head.matchAll(/\|\s*\d+>\s*!(.*)/g)) kw.push(m[1]);
      const k = kw.join(' ');
      for (const [re, j] of [[/\b(opt|copt|zopt|optts|tightopt|verytightopt|looseopt)\b/i, 'opt'], [/\b(freq|numfreq|anfreq)\b/i, 'freq'], [/\birc\b/i, 'irc']]) if (re.test(k)) info.jobs.push(j);
      info.freq = info.jobs.indexOf('freq') >= 0 || /NORMAL MODES/.test(tail);
      info.qmmm = /\b(qmmm|qm\/xtb|qm\/qm2|qm\/mm)\b/i.test(k) || /%qmmm/i.test(head);
      info.atoms = tableRows(head, 'CARTESIAN COORDINATES (ANGSTROEM)', 1, /^\s*[A-Za-z][\w:>*]*\s+-?\d/);
      const cy = lastOf(/GEOMETRY OPTIMIZATION CYCLE\s+(\d+)/g, tail) || lastOf(/GEOMETRY OPTIMIZATION CYCLE\s+(\d+)/g, head);
      if (cy) info.steps = +cy[1];
    } else {   // XYZ trajectory: frames from the size of the first one
      const L = head.split(/\r?\n/);
      let i = 0;
      while (i < L.length && !L[i].trim()) i++;
      const n = parseInt(L[i], 10);
      if (n > 0 && L.length > i + n + 2) {
        info.atoms = n;
        const bytes = new TextEncoder().encode(L.slice(i, i + n + 2).join('\n') + '\n').length;
        info.frames = Math.max(1, Math.round(R.size / bytes));
      }
    }
    return info;
  }
  const MBs1 = b => (b >= 1073741824 ? (b / 1073741824).toFixed(1) + ' GB' : Math.max(1, Math.round(b / 1048576)) + ' MB');
  /**
   * The loading dialog: what to load of a big output (info.kind 'text': geometries and vibrations) or trajectory
   * (info.kind 'traj': frames). Resolves to the choice (text: {steps, modes} for condense.js; traj: {stride, first,
   * last, tail} for BM.parseSource) or null when cancelled. prev: the choice made before (Load other parts).
   */
  function askLoad(fileName, info, prev) {
    return new Promise(resolve => {
      const dlg = $('#dlg');
      let done = false;
      const finish = v => { if (done) return; done = true; dlg.removeEventListener('close', onClose); resolve(v); };
      const onClose = () => finish(null);
      const traj = info.kind === 'traj', xyz = info.format === 'xyz', frames = traj || xyz;
      const N = info.atoms || 0, perStep = 24 * N * (frames ? 1 : 2);
      const total = frames ? info.frames : info.steps;
      const unit = frames ? 'frame' : 'step';
      const budget = Math.min(2 * 1073741824, (ramLimit() || 4 * 1073741824) * 0.5);
      const ps = prev && prev.steps ? prev.steps : null, pm = prev && prev.modes ? prev.modes : null;
      // defaults: everything that fits
      let defMode = 'all', every = 10;
      if (total && perStep && total * perStep > budget) { defMode = 'every'; every = Math.max(2, Math.ceil(total * perStep / budget)); }
      const sMode = ps ? ps.mode : defMode;
      const nModes = info.freq && N ? 3 * N - 6 : 0, modeBytes = nModes * 3 * N * (info.bytesPerValue || 4);
      let mMode = pm ? pm.mode : 'all', firstK = pm && pm.n ? pm.n : 100;
      if (!pm && modeBytes > budget) { mMode = 'first'; firstK = Math.max(20, Math.floor(budget / (3 * N * (info.bytesPerValue || 4)))); }
      const facts = [MBs1(info.size), traj ? (info.format || 'trajectory') : formatName({ format: info.format }) + (info.qmmm ? ' (QM/MM)' : ''),
        N ? N.toLocaleString() + ' atoms' : '', info.jobs && info.jobs.length ? info.jobs.join(', ') : '',
        total ? (info.exact === false || xyz ? 'about ' : '') + total.toLocaleString() + ' ' + unit + 's' : ''].filter(Boolean).join(' · ');
      const num = (k, v, min) => `<input type="number" data-ld="${k}" value="${v}" min="${min || 1}" step="1" style="width:84px">`;
      const radio = (group, value, checked, html) => `<label class="ld-opt"><input type="radio" name="${group}" value="${value}" ${checked ? 'checked' : ''}> <span>${html}</span></label>`;
      const lastN = ps && ps.mode === 'last' ? ps.n : 20, from = ps && ps.mode === 'range' ? ps.from : 1, to = ps && ps.mode === 'range' && ps.to ? ps.to : (total || 100);
      const body = `<p class="note" style="margin-top:0"><b>${esc(fileName)}</b>: ${esc(facts)}.</p>
        <div class="ld-h">${frames ? 'Frames' : 'Geometries'}</div>
        ${radio('ld-steps', 'all', sMode === 'all', 'all ' + unit + 's' + (total && perStep ? ` <span class="dim">(about ${MBs1(total * perStep)})</span>` : ''))}
        ${radio('ld-steps', 'every', sMode === 'every', 'every ' + num('every', ps && ps.mode === 'every' ? ps.n : every, 2) + ` ${unit}`)}
        ${radio('ld-steps', 'last', sMode === 'last', 'only the last ' + num('last', lastN) + ` ${unit}s`)}
        ${radio('ld-steps', 'range', sMode === 'range', `${unit}s ` + num('from', from) + ' to ' + num('to', to))}
        ${frames ? '' : '<p class="note">The energies and convergence data of every step are read in any case, and the first and last geometry of each job are always loaded.</p>'}
        ${nModes ? `<div class="ld-h">Vibrations <span class="dim">(${nModes.toLocaleString()})</span></div>
          ${radio('ld-modes', 'all', mMode === 'all', 'all, with their displacements' + (modeBytes ? ` <span class="dim">(about ${MBs1(modeBytes)})</span>` : ''))}
          ${radio('ld-modes', 'first', mMode === 'first', 'displacements of the first ' + num('first', firstK) + ' only')}
          ${radio('ld-modes', 'none', mMode === 'none', 'frequencies and intensities only (no animation)')}` : ''}
        <p class="note" data-ld="est"></p>`;
      const read = form => {
        const q = sel => form.querySelector(sel);
        const val = k => Math.max(1, Math.floor(+q(`[data-ld="${k}"]`).value || 1));
        const sm = (q('input[name="ld-steps"]:checked') || { value: 'all' }).value;
        const steps = sm === 'every' ? { mode: 'every', n: Math.max(2, val('every')) } : sm === 'last' ? { mode: 'last', n: val('last') }
          : sm === 'range' ? { mode: 'range', from: Math.min(val('from'), val('to')), to: Math.max(val('from'), val('to')) } : { mode: 'all' };
        if (traj) {   // BM.parseSource options (0-based frames)
          if (steps.mode === 'every') return { stride: steps.n, first: 0, last: null, tail: null, choice: steps };
          if (steps.mode === 'last') return { stride: 1, tail: steps.n, choice: steps };
          if (steps.mode === 'range') return { stride: 1, first: steps.from - 1, last: steps.to - 1, tail: null, choice: steps };
          return { stride: 1, first: 0, last: null, tail: null, choice: steps };
        }
        const mm = q('input[name="ld-modes"]:checked');
        const modes = !mm || mm.value === 'all' ? { mode: 'all' } : mm.value === 'none' ? { mode: 'none' } : { mode: 'first', n: val('first') };
        return { steps, modes };
      };
      const estimate = form => {
        const c = read(form), st = traj ? c.choice : c.steps;
        let n = total || 0;
        if (n) n = st.mode === 'every' ? Math.ceil(n / st.n) : st.mode === 'last' ? Math.min(n, st.n) : st.mode === 'range' ? Math.max(0, Math.min(n, st.to) - st.from + 1) : n;
        let b = n * perStep;
        if (!traj && c.modes && nModes) b += c.modes.mode === 'all' ? modeBytes : c.modes.mode === 'first' ? Math.min(nModes, c.modes.n) * 3 * N * (info.bytesPerValue || 4) : 0;
        const e = form.querySelector('[data-ld="est"]');
        if (e) e.textContent = b ? 'Memory for what is kept: about ' + MBs1(b) + (total ? '' : ' plus the geometries') + '.' : '';
      };
      dialog('Load part of a large file', body, [
        { label: 'Cancel', onClick: () => { finish(null); } },
        { label: 'Load', primary: true, onClick: form => { finish(read(form)); } }
      ], form => {
        const bodyEl = form.querySelector('.dlg-body');
        bodyEl.addEventListener('input', e => {   // typing a number chooses its option
          const lab = e.target.closest('.ld-opt');
          if (lab && e.target.type === 'number') { const r = lab.querySelector('input[type=radio]'); if (r) r.checked = true; }
          estimate(form);
        });
        bodyEl.addEventListener('change', () => estimate(form));
        estimate(form);
        const pb = form.querySelector('.dlg-foot .btn.primary'); if (pb) pb.focus();
      });
      dlg.addEventListener('close', onClose);
    });
  }
  /** The first bytes of a file whose extension may hide a binary format (an Amber NetCDF restart named .rst7). */
  async function looksBinary(src, name) {
    if (!BM.maybeBinary || !BM.maybeBinary(name)) return false;
    try {
      let head = null;
      if (src.kind === 'bytes') head = src.data.subarray(0, 16);
      else if ((src.kind === 'blob' || src.kind === 'chunks') && BIGF) { const R = BIGF.byteReader(src.data); head = R && R.size ? await R.read(0, Math.min(16, R.size)) : null; }
      return !!head && BM.sniffBinary(name, head);
    } catch (e) { return false; }
  }
  /** Reads a file into a dataset: the fast path for small files, slices and lines for big ones (with load choices). */
  async function readSource(src, name, gz, size, choice, onProg) {
    let data, note = '', prog = null;
    try {
      if (BM.isBinaryFormat(name) || (!gz && await looksBinary(src, name))) {
        if (!gz && size !== null && size >= bigSize() && BM.hasSourceReader && BM.hasSourceReader(name) && (src.kind === 'blob' || src.kind === 'chunks')) {
          prog = onProg(src.name, size);   // DCD / TRR / XTC / NetCDF: frames read from slices
          data = await BM.parseSource(name, src.data, Object.assign({ onProgress: prog.onProgress, signal: prog.signal }, choice || {}));
          data.meta.sourceSize = size;
          if (choice && choice.choice) data.meta.choice = choice.choice;
        } else {
          if (src.kind === 'stream' || (size !== null && size >= bigSize())) prog = onProg(src.name, size);
          data = BM.parseFileBinary(name, await readBytes(src, gz, prog));
        }
      } else if (BIGF && (gz || src.kind === 'stream' || size === null || size >= bigSize())) {
        prog = onProg(src.name, size, size !== null && size < bigSize() ? 300 : 0);
        try { data = await readLinesAndParse(src, name, gz, prog, choice); }
        catch (e) {
          // unrecognised small files: Open Babel needs the bytes (re-read; streams cannot be re-read)
          if (!/not recognised|binary file/i.test(e.message || '') || src.kind === 'stream' || (size !== null && size >= bigSize())) throw e;
          const conv = await viaOpenBabel(name, await readBytes(src, gz), e);
          if (!conv) throw e;
          data = conv; note = ' (converted with Open Babel)';
        }
      } else {
        let bytes = await readBytes(src, gz);
        let text = null;
        try { text = new TextDecoder('utf-8').decode(bytes); } catch (e) { text = null; }   // too large for one string
        if (text === null) {
          if (!BIGF) throw new Error('This file is too large to open as text.');
          bytes = null;
          prog = onProg(src.name, size);
          data = await readLinesAndParse(src, name, false, prog, choice);
        } else {
          try { data = BM.parseFile(name, text); }
          catch (e) {
            text = null;
            const conv = await viaOpenBabel(name, bytes, e);
            if (!conv) throw e;
            data = conv; note = ' (converted with Open Babel)';
          }
        }
        text = null; bytes = null;
      }
    } finally {
      if (prog) prog.done();
    }
    return { data, note };
  }
  /** Asks what to load when a big output or trajectory is large enough; null: cancelled; {}: no question asked. */
  async function chooseLoad(src, name, gz, size, prev) {
    const seekable = !gz && size !== null && (src.kind === 'blob' || src.kind === 'chunks');
    if (!seekable || (!prev && size < askSize())) return {};
    if (BM.isBinaryFormat(name)) {
      if (!(size >= bigSize() && BM.hasSourceReader && BM.hasSourceReader(name) && BM.probeSource)) return {};
      let pr = null;
      try { pr = await BM.probeSource(name, src.data); } catch (e) { pr = null; }
      if (!pr) return {};
      const fmt = (BM.formats.find(f => f.extensions.indexOf(name.split('.').pop().toLowerCase()) >= 0) || {}).name;
      const choice = await askLoad(src.name, Object.assign({ kind: 'traj', size, format: fmt }, pr), prev);
      return choice ? { traj: choice } : null;
    }
    if (size < bigSize()) return {};
    let info = null;
    try { info = await probeText(src, name); } catch (e) { info = null; }
    if (!info) return {};
    const choice = await askLoad(src.name, Object.assign({ kind: 'text' }, info), prev);
    return choice ? { text: choice } : null;
  }
  async function openFile(file, nameHint) {
    let src;
    try { src = asSource(file, nameHint); } catch (e) { toast(e.message || String(e), true); return null; }
    const gz = /\.gz$/i.test(src.name), name = src.name.replace(/\.gz$/i, '');
    const size = Number.isFinite(src.size) ? src.size : null;
    status('Reading ' + src.name + '…');
    try {
      const ask = await chooseLoad(src, name, gz, size, null);
      if (!ask) { layersRequest = null; toast('Opening ' + src.name + ' was cancelled.'); return null; }
      const choice = ask.traj || ask.text || null;
      const t0 = performance.now();
      const { data, note } = await readSource(src, name, gz, size, choice, beginProgress);
      console.info(`[BlueMolecule] parsed ${name} (${data.format}) in ${(performance.now() - t0).toFixed(0)} ms`);
      const d = openData(name, data, note);
      if (layersRequest) layersOpened(d);
      if (d && d.data === data && data.meta && data.meta.loaded && !gz && (src.kind === 'blob' || src.kind === 'chunks')) {
        d.loadSource = { src, name: src.name, size, choice: ask.text || { steps: { mode: 'all' }, modes: { mode: 'all' } } };
      }
      if (d && data.meta && data.meta.totalFrames > data.frames.length && !ask.traj) {
        toast(`${name}: ${data.frames.length} of ${data.meta.totalFrames} frames loaded (every ${ordinal(data.meta.stride)})`);
      }
      return d;
    } catch (e) {
      layersRequest = null;
      if (BIGF && BIGF.isAbort(e)) { toast('Opening ' + src.name + ' was cancelled.'); return null; }
      console.error(e);
      toast(conciseError(BIGF ? BIGF.describeError(e, src.name, size) : (e.message || String(e))), true);
      return null;
    }
  }
  /** A document's data replaced by a new reading of the same file (Load other parts). */
  function replaceData(d, data) {
    stopAnim(true);
    d.data = data;
    d.frameIndex = Math.max(0, data.frames.length - 1);
    if (data.vibrations && data.vibrations.frameIndex < data.frames.length) d.frameIndex = data.vibrations.frameIndex;
    d.selection = []; d.vibIndex = -1; d.stateIndex = -1; d.job = -1; d.convergedOnly = false;
    d.bonds.clear(); d._conv = null; d._prof = null; d._vibVec = null; d._animXyz = null; d.undo = []; d.redo = [];
    d.chargeScheme = firstChargeScheme(data);
    if (state.active === d) { refreshScene(true); renderDocs(); renderTabs(); renderPanel(); renderPlot(); renderPlaybar(); updateStatus(); }
  }
  /** Info tab "Load other parts…": the loading dialog again, then the same file read with the new choice. */
  async function reloadPart(d) {
    if (!d) return;
    const L = d.loadSource, T = d.trajLoad;
    try {
      if (L) {
        const name = L.name.replace(/\.gz$/i, '');
        const ask = await chooseLoad(L.src, name, false, L.size, L.choice);
        if (!ask || !ask.text) return;
        const { data } = await readSource(L.src, name, false, L.size, ask.text, beginProgress);
        L.choice = ask.text;
        replaceData(d, data);
        toast(`${L.name}: ${data.frames.length} geometr${data.frames.length === 1 ? 'y' : 'ies'} loaded`);
      } else if (T && T.reload) {
        const choice = await askLoad(T.name, { kind: 'traj', size: T.size, format: T.format, frames: T.total, atoms: T.atoms, perFrame: 24 * T.atoms, exact: T.exact }, { steps: T.choice });
        if (!choice) return;
        const prog = beginProgress(T.name, T.size);
        let td;
        try { td = await T.reload(Object.assign({ onProgress: prog.onProgress, signal: prog.signal }, choice)); } finally { prog.done(); }
        const t = d.data.frames[0];
        const frames = td.frames.map((f, k) => ({ Z: t.Z, xyz: f.xyz, energy: isNum(f.energy) ? f.energy : null, energies: null, label: f.label || `${T.name} ${k + 1}`,
          converged: false, opt: null, atomNames: t.atomNames, resNames: t.resNames, resIds: t.resIds, chains: t.chains || null, cell: f.cell || null, time: f.time }));
        d.data.frames.splice(T.start, T.count, ...frames);
        Object.assign(T, { count: frames.length, total: td.meta.totalFrames, stride: td.meta.stride, first: td.meta.first, last: td.meta.last, choice: choice.choice, reload: td.meta.reload || T.reload });
        stopAnim(true);
        d.frameIndex = T.start; d.bonds.clear(); d._conv = null; d._prof = null; d.job = -1;
        if (state.active === d) { refreshScene(true); renderTabs(); renderPanel(); renderPlot(); renderPlaybar(); updateStatus(); }
        toast(`${T.name}: ${frames.length} of ${T.total} frames loaded`);
      }
    } catch (e) {
      if (BIGF && BIGF.isAbort(e)) { toast('Loading was cancelled.'); return; }
      console.error(e);
      toast(conciseError(BIGF ? BIGF.describeError(e, (L || T).name, (L || T).size) : (e.message || String(e))), true);
    }
  }
  /** Formats BlueMolecule does not read: convert with Open Babel (to SDF, else CML) when the host provides it. */
  async function viaOpenBabel(name, bytes, err) {
    const ob = BM.native && BM.native.obabel, m = /\.([A-Za-z0-9_+-]+)$/.exec(name);
    if (!ob || !m || !/not recognised|binary file/i.test((err && err.message) || '')) return null;
    let ok = false;
    try { ok = await ob.available(); } catch (e) { ok = false; }
    if (!ok) return null;
    status('Converting ' + name + ' with Open Babel…');
    let last = null;
    for (const fmt of ['sdf', 'cml']) {
      try {
        const out = await ob.convert(bytes, m[1].toLowerCase(), fmt);
        if (!out || !out.trim()) continue;
        const ds = BM.parseFile(name.replace(/\.[^.]+$/, '') + '.' + fmt, out);
        ds.title = name; ds.convertedBy = 'Open Babel';
        ds.warnings = (ds.warnings || []).concat(['Converted with Open Babel (' + m[1] + ' → ' + fmt.toUpperCase() + ').']);
        return ds;
      } catch (e) { last = e; }
    }
    throw new Error('Open Babel could not convert "' + name + '"' + (last ? ': ' + last.message : '.'));
  }
  /** Several files at once: topologies and structures before coordinate-only files (DCD, TRR, XTC, NetCDF, inpcrd…). */
  const coordOnly = f => /\.(dcd|trr|xtc|nc|ncdf|netcdf|ncrst|inpcrd|rst7|restrt|rst|mdcrd|crd)(\.gz)?$/i.test(f.name) ? 1 : 0;
  async function openFiles(list) {
    for (const f of Array.from(list).sort((a, b) => coordOnly(a) - coordOnly(b))) await openFile(f);
  }
  /** A document pushed by a native shell: cb(name, Uint8Array) for small files, cb(name, {size, readChunk}) or
   *  cb({name, size, readChunk}) for big ones (read in slices, never copied whole). */
  function hostItem(name, data) {
    if (name && typeof name === 'object') return name;
    if (data && typeof data.readChunk === 'function') {
      return { name: String(name || data.name || 'untitled'), size: data.size, readChunk: (o, n) => data.readChunk(o, n) };
    }
    if (data && (data instanceof Uint8Array || ArrayBuffer.isView(data) || data instanceof ArrayBuffer)) return new File([data], String(name || 'untitled'));
    return BIGF ? BIGF.toSource(data, name) : new File([data], String(name || 'untitled'));
  }
  /** Open button, Ctrl/Cmd+O: the shell's own dialog when it has one (big files then arrive as chunk objects). */
  async function openDialog() {
    const nat = BM.native;
    if (!nat || typeof nat.openDialog !== 'function') { $('#file-input').click(); return 'input'; }
    let r;
    try { r = await nat.openDialog({ multiple: true }); }
    catch (e) { console.error(e); toast('The open dialog failed: ' + (e && e.message ? e.message : e), true); return; }
    if (r === null) { $('#file-input').click(); return 'input'; }          // not handled by this host
    if (r === undefined || typeof r === 'boolean') return;          // cancelled, or the files come through onOpenFile
    const list = Array.isArray(r) ? r : (typeof r === 'object' && typeof r.length === 'number' && !r.readChunk) ? Array.from(r) : [r];
    const items = [];
    for (const x of list) { try { items.push(hostItem(x)); } catch (e) { toast(e.message || String(e), true); } }
    if (items.length) await openFiles(items);
  }
  /** Documents from the macOS app (Finder, Dock, open -a) arrive one message each, in no useful order: collect the
   *  ones that come together and open them like a multi-file drop (topology before trajectory, one at a time). */
  let hostFiles = [], hostTimer = 0;
  function openHostFile(file) {
    hostFiles.push(file);
    clearTimeout(hostTimer);
    hostTimer = setTimeout(() => { const list = hostFiles; hostFiles = []; openFiles(list); }, 80);
  }
  function openText(name, text) {
    try { return openData(name, BM.parseFile(name, text)); } catch (e) { toast(conciseError(e.message), true); return null; }
  }
  /** A parsed dataset: coordinates without atoms go to a matching document, the rest opens as a new document. */
  function openData(name, data, note) {
    if (data.meta && data.meta.needsTopology) return attachCoordinates(name, data);
    const a = state.active, n0 = data.frames[0] && data.frames[0].Z ? data.frames[0].Z.length : -1;
    if (a && a.pending && !(data.meta && data.meta.topologyOnly) && n0 === a.data.frames[0].Z.length) return attachFrames(a, name, data);
    return addDoc(name, data, note);
  }
  function attachCoordinates(name, data) {
    const n = data.meta.natoms || (data.frames[0] ? data.frames[0].xyz.length / 3 : 0);
    const a = state.active;
    const target = a && a.data.frames[0].Z.length === n && curFrame(a).Z.length === n ? a : state.docs.find(d => d.pending && d.data.frames[0].Z.length === n);
    if (!target) {
      toast(`${name} has coordinates for ${n} atoms but no atom types: open the topology/structure first (prmtop, psf, pdb, gro…) with ${n} atoms.`, true);
      return null;
    }
    return attachFrames(target, name, data);
  }
  /** Append (or, for a topology waiting for coordinates, set) frames from another file, with the document's atoms. */
  function attachFrames(d, name, data) {
    const ds = d.data, t = ds.frames[0];
    const frames = data.frames.map((f, k) => ({
      Z: t.Z, xyz: f.xyz, energy: isNum(f.energy) ? f.energy : null, energies: null, label: f.label || (data.frames.length > 1 ? `${name} ${k + 1}` : name),
      converged: false, opt: null, atomNames: t.atomNames, resNames: t.resNames, resIds: t.resIds, chains: t.chains || null, cell: f.cell || null, time: f.time
    }));
    const wasPending = d.pending;
    if (wasPending) { ds.frames = frames; d.pending = false; ds.meta.topologyOnly = false; d.frameIndex = 0; }
    else { ds.frames = ds.frames.concat(frames); d.frameIndex = ds.frames.length - frames.length; }
    const m = data.meta || {};
    if (typeof m.reload === 'function' && m.totalFrames) {   // a big trajectory read in slices: "Load other parts" can read it again
      d.trajLoad = { name, reload: m.reload, total: m.totalFrames, stride: m.stride, first: m.first, last: m.last, start: d.frameIndex, count: frames.length,
        atoms: m.natoms || t.Z.length, size: m.sourceSize || null, format: formatName(data), exact: true, choice: m.choice || null };
    }
    if (data.cell && !ds.cell) ds.cell = data.cell;
    if (data.cellOrigin && !ds.cellOrigin) ds.cellOrigin = data.cellOrigin;
    if (ds.bonds && ds.bonds.length) d.topoAll = true;
    d.bonds.clear(); d._conv = null; d._prof = null; d.job = -1; d.convergedOnly = false;
    ds.sources = (ds.sources || []).concat([name]);
    if (frames.length > 1 && d.tab === 'info') d.tab = 'steps';
    if (state.active !== d) activate(d);
    else { stopAnim(true); refreshScene(!wasPending); renderDocs(); renderTabs(); renderPanel(); renderPlot(); renderPlaybar(); updateStatus(); }
    toast(`${name}: ${frames.length} geometr${frames.length > 1 ? 'ies' : 'y'} ${wasPending ? 'loaded into' : 'added to'} ${d.name}`);
    return d;
  }

  // ---------------- saving ----------------
  /** Every save goes through here: the test hook, the native app's Save dialog, the browser's file picker, a download. */
  async function saveFile(name, data, mime) {
    const blob = data instanceof Blob ? data : new Blob([data], { type: mime || 'text/plain' });
    if (window.__bmSaveHook) { await window.__bmSaveHook(name, blob); status('Saved ' + name); return true; }
    const nat = BM.native;
    if (nat && typeof nat.save === 'function') {
      let r = null;
      try { r = await nat.save(name, blob); } catch (e) { toast('Could not save ' + name + ': ' + e.message, true); return false; }
      if (r === true) { const p = nat.lastSave && nat.lastSave.path; toast('Saved to ' + (p || name)); return true; }
      if (r === false) { const er = nat.lastSave && nat.lastSave.error; if (er) toast('Could not save ' + name + ': ' + er, true); else status('Save cancelled'); return false; }
    }
    if (window.showSaveFilePicker && window.isSecureContext !== false) {
      try {
        const h = await window.showSaveFilePicker({ suggestedName: name });
        const w = await h.createWritable(); await w.write(blob); await w.close();
        toast('Saved to ' + h.name); return true;
      } catch (e) { if (e && e.name === 'AbortError') { status('Save cancelled'); return false; } }
    }
    const a = document.createElement('a');
    a.href = URL.createObjectURL(blob); a.download = name;
    document.body.appendChild(a); a.click();
    setTimeout(() => { URL.revokeObjectURL(a.href); a.remove(); }, 1500);
    toast(name + ' downloaded (see your browser’s downloads)');
    return true;
  }
  const baseName = d => d.name.replace(/\.[^.]+$/, '');

  // ---------------- scene ----------------
  function renderOptions() {
    const v = state.view;
    return { style: v.style, showH: v.showH, perspective: v.perspective, fog: v.fog, axes: v.axes, bg: BG[v.bg] || BG.navy, atomScale: v.atomScale, bondRadius: v.bondRadius, spacefillScale: v.vdwScale, showDummies: v.dummies,
      material: v.material, ao: +v.ao || 0, outline: +v.outline || 0, dof: +v.dof || 0, shadows: !!v.shadows, floor: !!v.floor, gradient: !!v.gradient, bg2: BG2[v.bg] || null, effectsLarge: !!v.fxLarge };
  }
  function refreshScene(keepView) {
    const d = state.active;
    R.setOptions(renderOptions());
    syncFxNote();
    $('#viewport').style.background = 'rgb(' + (BG[state.view.bg] || BG.navy).map(c => Math.round(c * 255)).join(',') + ')';
    $('#pending').classList.toggle('hidden', !(d && d.pending));
    if (d && d.pending) $('#pending-msg').textContent = `${d.name}: ${d.data.frames[0].Z.length} atoms. Open a coordinate or trajectory file with the same atom count (inpcrd, NetCDF, DCD, TRR, XTC, PDB, GRO…).`;
    if (!d || d.pending) { R.setRepresentation(null, true); R.setScene([], new Float64Array(0), new Uint32Array(0)); pushReps(null); R.setSurfaces([]); R.setArrows(null); refreshEnds(); refreshExtras(); updateMeasure(); return; }
    const f = curFrame(d);
    if (!keepView) R.setOrientation(principalQuat(f.Z, f.xyz));
    const fi = d.frameIndex;
    const rr = computeReps(d);
    R.setRepresentation(rr.rep, true);
    R.setScene(f.Z, d._animXyz && state.anim.kind === 'vib' ? d._animXyz : f.xyz, bondsFor(d, fi), keepView, state.view.bondOrders ? ordersFor(d, fi) : null);
    pushReps(rr);
    R.setSelection(d.selection);
    refreshArrows();
    refreshEnds();
    refreshSurfaces();
    refreshExtras();
    updateMeasure();
  }
  /** Hydrogen bonds, dipole arrow, unit cell and frozen-atom markers for the active document (they persist in the
   *  renderer across setScene, so they are always reset here). */
  function refreshExtras() {
    const d = state.active;
    if (!d || d.pending) { R.setDashedBonds([]); R.setExtraArrows([]); R.setCell(null); R.setAtomMarkers([]); return; }
    const f = curFrame(d), n = f.Z.length, v = state.view, xyz = R.xyz && R.xyz.length === 3 * n ? R.xyz : f.xyz;
    R.setDashedBonds(v.hbonds && n ? kit.hydrogenBonds(f.Z, xyz, bondsFor(d, d.frameIndex)) : [], [0.42, 0.85, 1.0]);
    const dp = d.data.dipole;
    if (v.dipole && dp && !f.edited && n && [dp.x, dp.y, dp.z].every(isNum)) {
      let cx = 0, cy = 0, cz = 0;
      for (let i = 0; i < n; i++) { cx += f.xyz[3 * i]; cy += f.xyz[3 * i + 1]; cz += f.xyz[3 * i + 2]; }
      const c = [cx / n, cy / n, cz / n], s = 0.6;
      R.setExtraArrows([{ from: c, to: [c[0] + dp.x * s, c[1] + dp.y * s, c[2] + dp.z * s], color: [1.0, 0.58, 0.18], radius: 0.06 }]);
    } else R.setExtraArrows([]);
    const cell = cellOf(d, f);
    R.setCell(v.cell && cell ? cell : null, d.data.cellOrigin || [0, 0, 0]);
    R.setAtomMarkers(d.frozen.filter(i => i < n));
  }
  /** Orientation that shows the molecule's widest face (principal axes), with a slight tilt for depth. */
  function principalQuat(Z, xyz) {
    const n = Z.length;
    if (n < 3) return BM.quat.norm(BM.quat.mul(BM.quat.axis([0, 1, 0], 0.35), BM.quat.axis([1, 0, 0], -0.3)));
    let cx = 0, cy = 0, cz = 0;
    for (let i = 0; i < n; i++) { cx += xyz[3 * i]; cy += xyz[3 * i + 1]; cz += xyz[3 * i + 2]; }
    cx /= n; cy /= n; cz /= n;
    const A = [[0, 0, 0], [0, 0, 0], [0, 0, 0]];
    for (let i = 0; i < n; i++) {
      const p = [xyz[3 * i] - cx, xyz[3 * i + 1] - cy, xyz[3 * i + 2] - cz];
      for (let a = 0; a < 3; a++) for (let b = 0; b < 3; b++) A[a][b] += p[a] * p[b];
    }
    const V = [[1, 0, 0], [0, 1, 0], [0, 0, 1]];
    for (let sweep = 0; sweep < 30; sweep++) {
      let off = Math.abs(A[0][1]) + Math.abs(A[0][2]) + Math.abs(A[1][2]);
      if (off < 1e-12) break;
      for (const [p, q] of [[0, 1], [0, 2], [1, 2]]) {
        if (Math.abs(A[p][q]) < 1e-15) continue;
        const th = 0.5 * Math.atan2(2 * A[p][q], A[q][q] - A[p][p]), c = Math.cos(th), s = Math.sin(th);
        for (let k = 0; k < 3; k++) { const akp = A[k][p], akq = A[k][q]; A[k][p] = c * akp - s * akq; A[k][q] = s * akp + c * akq; }
        for (let k = 0; k < 3; k++) { const apk = A[p][k], aqk = A[q][k]; A[p][k] = c * apk - s * aqk; A[q][k] = s * apk + c * aqk; }
        for (let k = 0; k < 3; k++) { const vkp = V[k][p], vkq = V[k][q]; V[k][p] = c * vkp - s * vkq; V[k][q] = s * vkp + c * vkq; }
      }
    }
    const order = [0, 1, 2].sort((a, b) => A[b][b] - A[a][a]);
    const rows = order.map(k => [V[0][k], V[1][k], V[2][k]]);
    const det = rows[0][0] * (rows[1][1] * rows[2][2] - rows[1][2] * rows[2][1]) - rows[0][1] * (rows[1][0] * rows[2][2] - rows[1][2] * rows[2][0]) + rows[0][2] * (rows[1][0] * rows[2][1] - rows[1][1] * rows[2][0]);
    if (det < 0) rows[2] = rows[2].map(v => -v);
    const m = rows, tr = m[0][0] + m[1][1] + m[2][2];
    let q;
    if (tr > 0) { const S = Math.sqrt(tr + 1) * 2; q = [(m[2][1] - m[1][2]) / S, (m[0][2] - m[2][0]) / S, (m[1][0] - m[0][1]) / S, 0.25 * S]; }
    else if (m[0][0] > m[1][1] && m[0][0] > m[2][2]) { const S = Math.sqrt(1 + m[0][0] - m[1][1] - m[2][2]) * 2; q = [0.25 * S, (m[0][1] + m[1][0]) / S, (m[0][2] + m[2][0]) / S, (m[2][1] - m[1][2]) / S]; }
    else if (m[1][1] > m[2][2]) { const S = Math.sqrt(1 + m[1][1] - m[0][0] - m[2][2]) * 2; q = [(m[0][1] + m[1][0]) / S, 0.25 * S, (m[1][2] + m[2][1]) / S, (m[0][2] - m[2][0]) / S]; }
    else { const S = Math.sqrt(1 + m[2][2] - m[0][0] - m[1][1]) * 2; q = [(m[0][2] + m[2][0]) / S, (m[1][2] + m[2][1]) / S, 0.25 * S, (m[1][0] - m[0][1]) / S]; }
    const tilt = BM.quat.mul(BM.quat.axis([0, 1, 0], 0.22), BM.quat.axis([1, 0, 0], -0.32));
    return BM.quat.norm(BM.quat.mul(tilt, BM.quat.norm(q)));
  }
  function resetView() { const d = state.active; if (!d) return R.resetView(); const f = curFrame(d); R.setOrientation(principalQuat(f.Z, f.xyz)); R.fit(); }
  function chargesFor(d) {
    if (!d || !d.data.charges || !d.chargeScheme) return null;
    const q = d.data.charges[d.chargeScheme];
    const f = curFrame(d);
    return q && q.length === f.Z.length ? q : null;
  }
  function applyColors() {
    const d = state.active;
    if (!d) return;
    if (state.view.color === 'charge') {
      const q = chargesFor(d);
      if (q) {
        const m = Math.max(0.2, ...q.map(Math.abs));
        const c = new Float32Array(q.length * 3);
        q.forEach((v, i) => {
          const t = Math.max(-1, Math.min(1, v / m));
          const col = t < 0 ? [1 + 0.75 * t, 1 + 0.55 * t, 1] : [1, 1 - 0.7 * t, 1 - 0.75 * t];
          c[3 * i] = col[0]; c[3 * i + 1] = col[1]; c[3 * i + 2] = col[2];
        });
        R.setAtomColors(c); return;
      }
    }
    R.setAtomColors(null);
  }
  /** Displacement arrows of the selected vibration: the longest is state.vib.arrowLen Å, whatever the amplitude. */
  function refreshArrows() {
    const d = state.active;
    const a = d && d.vibIndex >= 0 && state.vib.arrows && d._vibVec ? vibVec(d, d.vibIndex, state.vib.arrowLen) : null;
    R.setArrows(a);
  }
  /** Translucent structures at both ends of the selected vibration (the animation's amplitude), when shown. */
  function refreshEnds() {
    const d = state.active;
    const on = !!(d && !d.pending && state.vib.ends && d.vibIndex >= 0 && d._vibVec);
    const v = on ? state.view : null;
    const key = on ? [d.id, d.vibIndex, state.vib.amplitude, v.style, v.atomScale, v.bondRadius, v.showH].join(':') : '';
    if (key === refreshEnds.key) return;
    refreshEnds.key = key;
    if (!on) { R.setVibEnds([]); return; }
    const eq = d.data.frames[d.data.vibrations.frameIndex], small = v.style === 'sticks' || v.style === 'wireframe';
    R.setVibEnds(BM.vibtools.endMeshes(eq.Z, eq.xyz, d._vibVec, {
      radius: i => (eq.Z[i] === 1 && !v.showH ? 0 : small ? 0.16 : Math.max(0.12, 0.85 * R.atomRadius(i))), bondRadius: small ? 0.1 : 0.065
    }));
  }
  refreshEnds.key = '';
  function surfaceMesh(s, sign) {
    const key = sign + ':' + s.iso;
    s.mesh = s.mesh || {};
    if (!s.mesh[key]) s.mesh = { [key]: BM.iso.extract(s.vol, sign * s.iso), ...(s.mesh[(-sign) + ':' + s.iso] ? { [(-sign) + ':' + s.iso]: s.mesh[(-sign) + ':' + s.iso] } : {}) };
    return s.mesh[key];
  }
  const hex2rgb = h => [parseInt(h.slice(1, 3), 16) / 255, parseInt(h.slice(3, 5), 16) / 255, parseInt(h.slice(5, 7), 16) / 255];
  function refreshSurfaces() {
    const d = state.active;
    if (!d) { R.setSurfaces([]); return; }
    const list = [];
    for (const s of d.surfaces) {
      if (!s.visible) continue;
      const mesh = surfaceMesh(s, 1);
      let colors = null;
      try { colors = espColors(d, s, mesh); } catch (e) { console.warn('[BlueMolecule] surface colours failed:', e); }
      list.push({ mesh, color: hex2rgb(s.colorPos), opacity: s.opacity, colors });
      if (s.both) list.push({ mesh: surfaceMesh(s, -1), color: hex2rgb(s.colorNeg), opacity: s.opacity });
    }
    R.setSurfaces(list);
  }

  // ---------------- representations (Display tab) ----------------
  const REPS = BM.reps;
  const STYLE_SHORT = { ballstick: 'Ball & stick', sticks: 'Sticks', lines: 'Lines', wireframe: 'Wire', spacefill: 'Spacefill', ghost: 'Ghost', hidden: 'Hidden', lipid: 'Lipid',
    cartoon: 'Cartoon', trace: 'Trace', rings: 'Rings', snfg: 'SNFG', markers: 'Markers' };
  const TOOLBAR_STYLES = ['ballstick', 'sticks', 'spacefill', 'wireframe'];
  const LAYER_SEL = /\b(layer|qm2?|mm|high|medium|low)\b/i;
  const baseRow = () => REPS.presetRows('all', null, state.view)[0];
  /** Layers of the document when they describe the current frame's atoms. */
  function layersFor(d, f) { const L = d.data.layers; return L && L.atomLayer && L.atomLayer.length === f.Z.length ? L : null; }
  /** Selection context of the current frame (cached per frame, bonds and positions; dynamic masks are cleared on use). */
  function repCtx(d) {
    const f = curFrame(d), bonds = bondsFor(d, d.frameIndex);
    const xyz = d._animXyz && state.anim.kind === 'vib' && d._animXyz.length === f.xyz.length ? d._animXyz : f.xyz;
    let c = d._rc;
    if (!c || c.frame !== f || c.bonds !== bonds || c.xyz !== xyz || c.n !== f.Z.length) {
      c = REPS.context({ frame: f, bonds, xyz, layers: layersFor(d, f), ssRecords: d.data.secondaryStructure || null });
      d._rc = c;
    }
    c.selection = d.selection; c.frozen = d.frozen; c.charges = chargesFor(d); c.lightBg = state.view.bg === 'white' || state.view.bg === 'gray';
    return c;
  }
  /** The first display of a document picks its preset: QM/MM layers, membrane, glycan, nucleic acid, protein or all atoms. */
  function ensureReps(d) {
    if (d.reps) return;
    const c = repCtx(d);
    d.repPreset = REPS.autoPreset(c);
    d.reps = REPS.presetRows(d.repPreset, c, state.view);
    d.repSel = 0;
  }
  /** Resolve the rows for the current frame: renderer representation (null = toolbar style), overlay meshes, pick
   *  proxies, markers, legend, per-row atom counts and selection errors. */
  function computeReps(d) {
    ensureReps(d);
    const c = repCtx(d), n = c.n;
    c.cache.clear();
    let rows = d.reps;
    if (!c.layers && rows.some(r => r.visible !== false && LAYER_SEL.test(r.selection)) && !rows.some(r => !LAYER_SEL.test(r.selection) && r.visible !== false && r.style !== 'hidden' && r.style !== 'markers'))
      rows = [baseRow()]; // layer rows without layer data for this frame (e.g. after deleting atoms): show every atom
    if (REPS.isDefault(rows)) { d._repRes = { rep: null, meshes: [], proxies: null, proxyR: 0, markers: [], legend: null, counts: rows === d.reps ? [n] : d.reps.map(() => null), errors: d.reps.map(() => null) }; return d._repRes; }
    const t0 = performance.now(), res = REPS.resolve(rows, c), meshes = [], proxies = [];
    let pr = 1.2;
    for (const ov of res.overlays) {
      const colors = REPS.colorize(ov.color, c, ov.mask);
      const g = BM.cartoon.build(ov.kind, c, ov.mask, colors, { cylinders: !!ov.row.cylinders, scheme: ov.color, scale: ov.row.scale > 0 ? +ov.row.scale : 1 });
      for (const m of g.meshes) { m.opacity = ov.opacity; meshes.push(m); }
      for (const i of g.proxies) proxies.push(i);
      if (g.proxies.length) pr = Math.max(pr, g.proxyRadius);
    }
    d._repRes = { rep: { style: res.style, scale: res.scale, alpha: res.alpha, colors: res.colors, extent: res.extent }, meshes, proxies, proxyR: pr, markers: res.markers,
      legend: legendFor(d, c, rows), counts: rows === d.reps ? res.counts : d.reps.map(() => null), errors: rows === d.reps ? res.errors : d.reps.map(() => null), ms: performance.now() - t0 };
    return d._repRes;
  }
  /** Overlays, proxies, markers, legend and colours of a computeReps result (the representation itself goes with setScene). */
  function pushReps(rr) {
    if (rr && rr.rep) { if (R.atomColors) R.setAtomColors(null); } else applyColors();
    R.setOverlays(rr ? rr.meshes : []);
    R.setPickProxies(rr ? rr.proxies : null, rr ? rr.proxyR : 0);
    R.setRepMarkers(rr ? rr.markers : []);
    R.setLegend(rr ? rr.legend : null);
  }
  /** Re-apply the rows of the active document without resetting the scene. */
  function refreshReps(d) {
    if (!d || d !== state.active || d.pending) return;
    const rr = computeReps(d);
    R.setRepresentation(rr.rep);
    pushReps(rr);
    updateRadius();
    updateDisplayCounts(d);
  }
  let repTimer = 0;
  function scheduleReps(d, ms) { clearTimeout(repTimer); repTimer = setTimeout(() => refreshReps(d), ms === undefined ? 60 : ms); }
  /** Overlay meshes follow animated positions (vibrations, drags) for small structures. */
  function refreshOverlays(d) {
    if (!d || !d._repRes || !d._repRes.meshes.length || curFrame(d).Z.length > 20000) return;
    const c = repCtx(d), res = REPS.resolve(d.reps, c), meshes = [];
    for (const ov of res.overlays) {
      const g = BM.cartoon.build(ov.kind, c, ov.mask, REPS.colorize(ov.color, c, ov.mask), { cylinders: !!ov.row.cylinders, scheme: ov.color, scale: ov.row.scale > 0 ? +ov.row.scale : 1 });
      for (const m of g.meshes) { m.opacity = ov.opacity; meshes.push(m); }
    }
    d._repRes.meshes = meshes;
    R.setOverlays(meshes);
  }
  const layerOfSel = sel => { const m = /^\s*(?:layer\s+([hml])|(qm|high)|(qm2|medium)|(mm|low))\s*$/i.exec(sel); return !m ? null : m[1] ? m[1].toUpperCase() : m[2] ? 'H' : m[3] ? 'M' : 'L'; };
  /** Legend of the QM/MM layers (style glyph, method, atom count) when rows select whole layers. */
  function legendFor(d, c, rows) {
    const L = c.layers;
    if (!L) return null;
    const vis = rows.filter(r => r.visible !== false);
    if (!vis.some(r => layerOfSel(r.selection))) return null;
    const cnt = {};
    for (const x of L.atomLayer) cnt[x] = (cnt[x] || 0) + 1;
    const oniom = L.scheme === 'oniom', name = { H: oniom ? 'High' : 'QM', M: oniom ? 'Medium' : 'QM2', L: oniom ? 'Low' : 'MM' };
    const items = [];
    for (const k of ['H', 'M', 'L']) {
      if (!cnt[k]) continue;
      const row = vis.find(r => layerOfSel(r.selection) === k), st = row ? row.style : 'hidden';
      const col = row && /^#[0-9a-f]{6}$/i.test(row.color) ? row.color : row && row.color === 'layer' ? REPS.toHex(REPS.LAYER_COL[k]) : null;
      items.push({ label: name[k], sub: [L.labels && L.labels[k], cnt[k] + ' atom' + (cnt[k] > 1 ? 's' : '')].filter(Boolean).join(' · '), glyph: st === 'wireframe' ? 'lines' : st, color: col });
    }
    const fz = vis.find(r => /^\s*frozen\b/i.test(r.selection));
    const nf = L.frozen ? L.frozen.filter(Boolean).length : 0;
    if (fz && nf) items.push({ label: 'Frozen', sub: nf + ' atoms', glyph: fz.style === 'lines' ? 'dim' : fz.style });
    const lk = vis.find(r => /^\s*link\s*$/i.test(r.selection) && r.style === 'markers');
    if (lk && L.linkAtoms && L.linkAtoms.length) items.push({ label: 'Link atoms', sub: String(L.linkAtoms.length), glyph: 'markers', color: /^#[0-9a-f]{6}$/i.test(lk.color) ? lk.color : '#ff7ad9' });
    return { title: oniom ? 'ONIOM layers' : L.scheme === 'qm-qm2' ? 'QM/QM2 layers' : 'QM/MM layers', items };
  }
  function applyPreset(d, id) {
    if (!d || d.pending) return;
    const c = repCtx(d);
    d.reps = REPS.presetRows(id, c, state.view); d.repPreset = id; d.repSel = 0;
    refreshReps(d);
    if (d.tab === 'display' || d.tab === 'info') renderPanel();
    const p = REPS.PRESETS.find(x => x.id === id);
    status('Display preset: ' + (p ? p.label : id));
  }
  /** Select the atoms of a selection expression in the view. */
  function selectBy(d, text) {
    const r = REPS.check(text, repCtx(d));
    if (!r.ok) { toast(r.error, true); return; }
    d.selection = REPS.indices(r.mask);
    R.setSelection(d.selection); updateMeasure();
    status(d.selection.length + ' atoms selected (' + text + ')');
  }

  // ---- Display tab
  const SWATCH = {
    element: 'conic-gradient(#909090 0 40%, #3050f8 40% 60%, #ff0d0d 60% 80%, #ffffff 80%)', layer: 'linear-gradient(90deg,#ffb347 0 33%,#4fd1c5 33% 66%,#8c9bb5 66%)',
    chain: 'linear-gradient(90deg,#4f8cff 0 33%,#f2994a 33% 66%,#3cc27a 66%)', residue: 'linear-gradient(90deg,#e63c3c 0 25%,#2ab04a 25% 50%,#3d5df0 50% 75%,#e8d43a 75%)',
    restype: 'linear-gradient(90deg,#c3c9d1 0 25%,#5ec26a 25% 50%,#e5534b 50% 75%,#4a86e8 75%)', ss: 'linear-gradient(90deg,#e2507a 0 40%,#f4c542 40% 70%,#d2d9e3 70%)',
    rainbow: 'linear-gradient(90deg,#3d5df0,#3cc2c2,#8bd36b,#f2c94c,#eb5757)', sugar: 'linear-gradient(90deg,#0090bc 0 25%,#00a651 25% 50%,#ffd400 50% 75%,#a17a4d 75%)',
    charge: 'linear-gradient(90deg,#4064ff,#ffffff,#ff4d3b)', bfactor: 'linear-gradient(90deg,#3b70d9,#ffffff,#d9453b)', occupancy: 'linear-gradient(90deg,#3b70d9,#ffffff,#d9453b)'
  };
  const swatchOf = r => (/^#[0-9a-f]{6}$/i.test(r.color) ? r.color : SWATCH[r.color] || SWATCH.element);
  const fmtN = v => (v === null || v === undefined ? '' : v >= 10000 ? (v / 1000).toFixed(v >= 100000 ? 0 : 1) + 'k' : String(v));
  const dbtn = (act, html, title, enabled, cls) => `<button class="btn small ${cls || ''}" data-act="${act}" title="${esc(title)}" ${enabled === false ? 'disabled' : ''}>${html}</button>`;
  function repRowsHtml(d) {
    const rr = d._repRes || {}, counts = rr.counts || [], errors = rr.errors || [];
    return d.reps.map((r, k) => {
      const err = errors[k], n = counts[k];
      return `<div class="rp ${k === d.repSel ? 'cur' : ''} ${r.visible === false ? 'off' : ''}" data-rp-row="${k}">
        <button class="rp-eye" data-rp-vis="${k}" title="${r.visible === false ? 'Show' : 'Hide'}" aria-label="Show or hide"><svg class="ic sm"><use href="#${r.visible === false ? 'i-eye-off' : 'i-eye'}"/></svg></button>
        <span class="rp-sw" style="background:${swatchOf(r)}" title="${esc(REPS.COLORS.find(c => c.id === r.color) ? REPS.COLORS.find(c => c.id === r.color).label : r.color)}"></span>
        <span class="rp-main"><span class="rp-name">${esc(r.name || 'Representation')}</span><span class="rp-q mono">${esc(r.selection)}</span></span>
        <span class="rp-style">${esc(STYLE_SHORT[r.style] || r.style)}</span>
        <span class="rp-n mono ${err ? 'err' : ''}" title="${err ? esc(err) : (n === null || n === undefined ? '' : n + ' atoms')}">${err ? '!' : fmtN(n)}</span>
      </div>`;
    }).join('');
  }
  function panelDisplay(d) {
    ensureReps(d);
    const r = d.reps[d.repSel], k = d.repSel, nr = d.reps.length;
    const presets = REPS.PRESETS.map(p => `<button data-act="rp-preset" data-preset="${p.id}" title="${esc(p.title)}">${esc(p.label)}${d.repPreset === p.id ? '<kbd>✓</kbd>' : ''}</button>`).join('');
    let editor = '<p class="note">No representations: everything is hidden. Add one or pick a preset.</p>';
    if (r) {
      const st = REPS.STYLES, uni = /^#[0-9a-f]{6}$/i.test(r.color);
      const opt = (list, cur) => list.map(x => `<option value="${x.id}" ${x.id === cur ? 'selected' : ''}>${esc(x.label)}</option>`).join('');
      const colorSel = REPS.COLORS.map(x => `<option value="${x.id}" ${(uni ? 'uniform' : r.color) === x.id ? 'selected' : ''}>${esc(x.label)}</option>`).join('');
      const overlay = (REPS.STYLES.find(x => x.id === r.style) || {}).kind === 'overlay';
      editor = `<div class="sec rp-edit">
        <div class="rp-qbox"><input type="text" id="rp-q" class="mono" spellcheck="false" autocomplete="off" value="${esc(r.selection)}" placeholder="e.g. protein and within 5 of ligand" title="Selection: keywords (protein, water, lipid, carbohydrate, ligand, layer H…), chain/resname/resid/name/element/index lists, within D of …, same residue as …, combined with and / or / not and parentheses"><span class="rp-qn mono" id="rp-qn"></span></div>
        <div class="rp-ac hidden" id="rp-ac" role="listbox"></div>
        <div class="rp-err" id="rp-err"></div>
        <div class="brow"><span class="blbl">Name</span><span class="bctl"><input type="text" id="rp-name" value="${esc(r.name || '')}" spellcheck="false"></span></div>
        <div class="brow"><span class="blbl">Style</span><span class="bctl"><select id="rp-style"><optgroup label="Atoms">${opt(st.filter(x => x.kind === 'atom'), r.style)}</optgroup><optgroup label="Overlays">${opt(st.filter(x => x.kind === 'overlay'), r.style)}</optgroup></select></span></div>
        <div class="brow"><span class="blbl">Colour</span><span class="bctl"><select id="rp-color">${colorSel}</select>${uni ? `<input type="color" id="rp-uni" value="${r.color}">` : ''}</span></div>
        <div class="brow"><span class="blbl">Opacity</span><span class="bctl rp-rng"><input type="range" id="rp-op" min="0.1" max="1" step="0.05" value="${r.opacity === undefined ? 1 : r.opacity}"><output>${Math.round((r.opacity === undefined ? 1 : r.opacity) * 100)}%</output></span></div>
        ${r.style === 'rings' || r.style === 'markers' || r.style === 'hidden' ? '' : `<div class="brow"><span class="blbl">Size</span><span class="bctl rp-rng"><input type="range" id="rp-sc" min="0.25" max="2" step="0.05" value="${r.scale || 1}"><output>${(r.scale || 1).toFixed(2)}×</output></span></div>`}
        ${r.style === 'cartoon' ? `<div class="brow"><span class="blbl"></span><span class="bctl"><label class="chk"><input type="checkbox" id="rp-cyl" ${r.cylinders ? 'checked' : ''}> Helices as cylinders</label></span></div>` : ''}
      </div>`;
    }
    const kw = ['protein', 'backbone', 'sidechain', 'nucleic', 'water', 'ions', 'lipid', 'carbohydrate', 'ligand', 'hydrogen', 'metal', 'helix', 'sheet', 'layer H', 'layer M', 'layer L', 'frozen', 'link', 'selected',
      'chain A', 'resname', 'resid 1-10', 'name CA', 'element C', 'within 5 of', 'same residue as', 'and', 'or', 'not'];
    return `<div class="sec dsp-head">
        <div class="dsp-bar">
          <div class="dropdown" id="dd-presets"><button class="btn small" data-dd title="Representation presets"><span>Presets</span><svg class="ic sm"><use href="#i-chev"/></svg></button><div class="menu">${presets}</div></div>
          <span class="bgrp">${dbtn('rp-add', '+ Add', 'Add a representation (of the selected atoms when there is a selection)')}${dbtn('rp-dup', '<svg class="ic sm"><use href="#i-dup"/></svg>', 'Duplicate this representation', !!r, 'icon-only')}${dbtn('rp-del', '<svg class="ic sm"><use href="#i-close"/></svg>', 'Delete this representation', !!r, 'icon-only')}</span>
          <span class="bgrp">${dbtn('rp-up', '<svg class="ic sm"><use href="#i-collapse"/></svg>', 'Move up (drawn earlier: later rows override atom styles)', k > 0, 'icon-only')}${dbtn('rp-down', '<svg class="ic sm"><use href="#i-chev"/></svg>', 'Move down', k < nr - 1, 'icon-only')}</span>
        </div>
      </div>
      <div class="replist" id="rp-list">${repRowsHtml(d)}</div>
      ${editor}
      <div class="sec bsec ${state.build.fold.kw ? 'folded' : ''}" data-sec="kw"><h3 class="fold" data-act="b-fold" data-sec="kw"><span>Selection keywords</span><svg class="ic sm chev"><use href="#i-chev"/></svg></h3>
        <div class="bbody"><div class="kwcloud">${kw.map(x => `<button class="kw" data-kw="${esc(x)}" title="${esc(REPS.KW_HINT[x.split(' ')[0]] || '')}">${esc(x)}</button>`).join('')}</div>
        <p class="note">Later rows override earlier ones for atom styles; overlays (cartoon, rings, SNFG, markers) add on top. Atoms no row shows are hidden.</p></div></div>`;
  }
  /** Counts and errors in the Display tab after a change, without re-rendering the editor. */
  function updateDisplayCounts(d) {
    if (!d || d.tab !== 'display') return;
    const list = $('#rp-list');
    if (list) list.innerHTML = repRowsHtml(d);
    updateQueryInfo(d);
  }
  function updateQueryInfo(d, text) {
    const q = $('#rp-q'), qn = $('#rp-qn'), er = $('#rp-err');
    if (!q || !qn || !er) return;
    const res = REPS.check(text === undefined ? q.value : text, repCtx(d));
    qn.textContent = res.ok ? res.count.toLocaleString() + ' atom' + (res.count === 1 ? '' : 's') : '';
    er.textContent = res.ok ? '' : res.error;
    q.classList.toggle('bad', !res.ok);
    return res;
  }
  let acState = null, qTimer = 0;
  function showAc(d) {
    const q = $('#rp-q'), box = $('#rp-ac');
    if (!q || !box) return;
    const c = REPS.complete(q.value, q.selectionStart, repCtx(d));
    const word = q.value.slice(c.from, q.selectionStart);
    const items = c.items.filter(it => it.text.toLowerCase() !== word.toLowerCase()).slice(0, 8);
    if (!items.length || document.activeElement !== q) { hideAc(); return; }
    acState = { from: c.from, to: c.to, items, k: word ? 0 : -1 };
    box.innerHTML = items.map((it, k) => `<button class="ac-it ${k === acState.k ? 'on' : ''}" data-ac="${k}" tabindex="-1"><span class="mono">${esc(it.text)}</span><span class="dim">${esc(it.hint || '')}</span></button>`).join('');
    box.classList.remove('hidden');
  }
  function hideAc() { acState = null; const box = $('#rp-ac'); if (box) { box.classList.add('hidden'); box.innerHTML = ''; } }
  function acceptAc(d, k) {
    const q = $('#rp-q');
    if (!q || !acState || !acState.items[k]) return;
    const it = acState.items[k], v = q.value, ins = it.text + (it.text === 'same residue as' || !/[()]$/.test(it.text) ? ' ' : '');
    q.value = v.slice(0, acState.from) + ins + v.slice(acState.to).replace(/^\s+/, '');
    const pos = acState.from + ins.length;
    q.setSelectionRange(pos, pos);
    hideAc();
    onQueryInput(d);
    showAc(d);
  }
  function onQueryInput(d) {
    const q = $('#rp-q'), r = d.reps[d.repSel];
    if (!q || !r) return;
    clearTimeout(qTimer);
    qTimer = setTimeout(() => {
      const res = updateQueryInfo(d);
      if (res && res.ok && q.value !== r.selection) { r.selection = q.value; refreshReps(d); }
    }, 140);
  }
  /** Data-act handlers of the Display tab ('rp-…'). */
  function displayAction(act, el) {
    const d = state.active;
    if (!d || d.pending) return;
    ensureReps(d);
    const k = d.repSel, rows = d.reps;
    switch (act) {
      case 'rp-preset': applyPreset(d, el.dataset.preset); return;
      case 'rp-add': {
        const sel = d.selection.length > 0;
        rows.splice(k + 1, 0, REPS.row(sel ? 'Selection' : 'New', sel ? 'selected' : 'none', 'ballstick'));
        d.repSel = Math.min(rows.length - 1, k + 1); d.repPreset = null;
        break;
      }
      case 'rp-dup': if (rows[k]) { const c = Object.assign({}, rows[k], { name: (rows[k].name || 'Copy') + ' (copy)', base: false }); rows.splice(k + 1, 0, c); d.repSel = k + 1; d.repPreset = null; } break;
      case 'rp-del': if (rows[k]) { rows.splice(k, 1); d.repSel = Math.max(0, Math.min(k, rows.length - 1)); d.repPreset = null; } break;
      case 'rp-up': if (k > 0) { [rows[k - 1], rows[k]] = [rows[k], rows[k - 1]]; d.repSel = k - 1; d.repPreset = null; } break;
      case 'rp-down': if (k < rows.length - 1) { [rows[k + 1], rows[k]] = [rows[k], rows[k + 1]]; d.repSel = k + 1; d.repPreset = null; } break;
      default: return;
    }
    refreshReps(d); renderPanel();
    if (act === 'rp-add') { const q = $('#rp-q'); if (q) { q.focus(); q.select(); } }
  }
  /** Clicks inside the Display tab (rows, eye buttons, autocomplete, keyword chips); true when handled. */
  function displayClick(e) {
    const d = state.active;
    if (!d || d.tab !== 'display') return false;
    const ac = e.target.closest('[data-ac]');
    if (ac) { e.preventDefault(); acceptAc(d, +ac.dataset.ac); return true; }
    const vis = e.target.closest('[data-rp-vis]');
    if (vis) { const r = d.reps[+vis.dataset.rpVis]; if (r) { r.visible = r.visible === false; refreshReps(d); } return true; }
    const row = e.target.closest('[data-rp-row]');
    if (row) { const k = +row.dataset.rpRow; if (k !== d.repSel) { d.repSel = k; renderPanel(); } return true; }
    const kw = e.target.closest('[data-kw]');
    if (kw) {
      const q = $('#rp-q');
      if (!q) return true;
      const v = q.value.trim(), w = kw.dataset.kw, op = /^(and|or|not)$/.test(w);
      q.value = !v || v === 'none' ? w : op || /\b(and|or|not|of|as)$/.test(v) ? v + ' ' + w : v + ' and ' + w;
      q.focus(); q.setSelectionRange(q.value.length, q.value.length);
      onQueryInput(d);
      return true;
    }
    return false;
  }
  function displayInput(e) {
    const d = state.active, t = e.target;
    if (!d || d.tab !== 'display' || !t.id || !t.id.startsWith('rp-')) return false;
    const r = d.reps[d.repSel];
    if (!r) return true;
    const syncBase = () => {
      if (r.base && d.reps.length === 1 && TOOLBAR_STYLES.includes(r.style)) { state.view.style = r.style; if (r.color === 'element' || r.color === 'charge') state.view.color = r.color; syncToolbar(); savePrefs(); }
    };
    switch (t.id) {
      case 'rp-q': onQueryInput(d); showAc(d); return true;
      case 'rp-name': r.name = t.value; { const nm = $(`[data-rp-row="${d.repSel}"] .rp-name`); if (nm) nm.textContent = t.value; } return true;
      case 'rp-style': r.style = t.value; syncBase(); refreshReps(d); renderPanel(); return true;
      case 'rp-color': r.color = t.value === 'uniform' ? (r.uniform || '#4f8cff') : t.value; syncBase(); refreshReps(d); renderPanel(); return true;
      case 'rp-uni': r.color = r.uniform = t.value; { const sw = $(`[data-rp-row="${d.repSel}"] .rp-sw`); if (sw) sw.style.background = t.value; } scheduleReps(d, 40); return true;
      case 'rp-op': r.opacity = parseFloat(t.value); t.nextElementSibling.textContent = Math.round(r.opacity * 100) + '%'; scheduleReps(d, 40); return true;
      case 'rp-sc': r.scale = parseFloat(t.value); t.nextElementSibling.textContent = r.scale.toFixed(2) + '×'; scheduleReps(d, 40); return true;
      case 'rp-cyl': r.cylinders = t.checked; refreshReps(d); return true;
    }
    return true;
  }
  function displayKey(e) {
    const d = state.active, t = e.target;
    if (!d || t.id !== 'rp-q') return false;
    if (acState && (e.key === 'ArrowDown' || e.key === 'ArrowUp')) {
      e.preventDefault();
      const n = acState.items.length;
      acState.k = ((acState.k < 0 ? (e.key === 'ArrowDown' ? -1 : 0) : acState.k) + (e.key === 'ArrowDown' ? 1 : -1) + n) % n;
      $$('#rp-ac .ac-it').forEach((b, k) => b.classList.toggle('on', k === acState.k));
      return true;
    }
    if (acState && acState.k >= 0 && (e.key === 'Enter' || e.key === 'Tab')) { e.preventDefault(); acceptAc(d, acState.k); return true; }
    if (e.key === 'Escape') { if (acState) { e.preventDefault(); hideAc(); } return true; }
    if (e.key === 'Enter') { e.preventDefault(); clearTimeout(qTimer); const r = d.reps[d.repSel], res = updateQueryInfo(d); if (res && res.ok && r) { r.selection = t.value; refreshReps(d); } hideAc(); return true; }
    return false;
  }

  function labelProvider() {
    const d = state.active;
    if (!d) return [];
    const f = curFrame(d), Z = f.Z, out = [], mode = state.view.labels;
    const n = Z.length;
    if (mode !== 'none' && n <= 5000) {
      const q = mode === 'charge' ? chargesFor(d) : null;
      for (let i = 0; i < n; i++) {
        let t;
        if (mode === 'symbol') t = el.symbol(Z[i]);
        else if (mode === 'number') t = String(i + 1);
        else if (mode === 'symnum') t = el.symbol(Z[i]) + (i + 1);
        else if (mode === 'charge') { if (!q) break; t = (q[i] >= 0 ? '+' : '') + q[i].toFixed(3); }
        out.push({ atom: i, text: t, size: 11.5 });
      }
    }
    if (state.view.bondLengths) {
      const b = bondsFor(d, d.frameIndex), xyz = R.xyz;
      for (let k = 0; k < b.length && k < 6000; k += 2) {
        const i = b[k], j = b[k + 1];
        if (!R.visible(i) || !R.visible(j)) continue;
        out.push({ pos: [(xyz[3 * i] + xyz[3 * j]) / 2, (xyz[3 * i + 1] + xyz[3 * j + 1]) / 2, (xyz[3 * i + 2] + xyz[3 * j + 2]) / 2], text: geom.dist(xyz, i, j).toFixed(3), size: 10.5, color: state.view.bg === 'white' || state.view.bg === 'gray' ? '#9a3412' : '#ffd38a' });
      }
    }
    const fc = f.formalCharges;
    if (fc && n <= 5000) {
      for (let i = 0; i < n; i++) {
        const q = fc[i];
        if (q) out.push({ atom: i, text: (Math.abs(q) > 1 ? Math.abs(q) : '') + (q > 0 ? '+' : '−'), size: 12, bold: true, dx: 10, dy: -10, color: q > 0 ? '#8fc3ff' : '#ff9d8f' });
      }
    }
    const sel = d.selection;
    if (sel.length && sel.length <= 4 && mode === 'none') sel.forEach(i => out.push({ atom: i, text: el.symbol(Z[i]) + (i + 1), size: 12, bold: true, dy: -16 }));
    const m = measurement(d);
    if (m && sel.length >= 2) {
      const xyz = R.xyz, P = i => [xyz[3 * i], xyz[3 * i + 1], xyz[3 * i + 2]];
      let p;
      if (sel.length === 2) p = P(sel[0]).map((v, k) => (v + P(sel[1])[k]) / 2);
      else if (sel.length === 3) p = P(sel[1]);
      else p = P(sel[1]).map((v, k) => (v + P(sel[2])[k]) / 2);
      out.push({ pos: p, text: m.text, size: 13, bold: true, color: '#ffd45e', dy: sel.length === 3 ? 18 : 0 });
    }
    return out;
  }

  // ---------------- selection & measurement ----------------
  /** Distance, angle or dihedral for 2–4 selected atoms (null otherwise). */
  function measurement(d) {
    const s = d.selection, xyz = R.xyz, n = xyz.length / 3;
    if (s.length < 2 || s.length > 4 || s.some(i => i >= n)) return null;
    if (s.length === 2) { const v = geom.dist(xyz, s[0], s[1]); return { kind: 'distance', value: v, text: v.toFixed(4) + ' Å', unit: 'Å' }; }
    if (s.length === 3) { const v = geom.angle(xyz, s[0], s[1], s[2]); return { kind: 'angle', value: v, text: v.toFixed(2) + '°', unit: '°' }; }
    const v = geom.dihedral(xyz, s[0], s[1], s[2], s[3]); return { kind: 'dihedral', value: v, text: v.toFixed(2) + '°', unit: '°' };
  }
  /** Measurement chip and lines; `light` (animation frames) skips the Build-panel refresh. */
  function updateMeasure(light) {
    const d = state.active, chip = $('#measure'), over = ui && d ? ui.chipHtml() : null;
    const s = d ? d.selection : [];
    if (over) { chip.innerHTML = over; chip.classList.remove('hidden'); }
    else if (!d || !s.length) chip.classList.add('hidden');
    else {
      const Z = curFrame(d).Z, xyz = R.xyz, m = measurement(d);
      let html;
      if (s.length > 4) html = `<span class="val first">${s.length}</span> <span>atoms selected</span>`;
      else {
        html = `<span class="atoms">${esc(s.map(i => el.symbol(Z[i]) + (i + 1)).join(' – '))}</span>`;
        if (m) html += ` <span>${m.kind}</span><span class="val">${m.text}</span>`;
        else html += `<span class="val">${fx(xyz[3 * s[0]], 4)}, ${fx(xyz[3 * s[0] + 1], 4)}, ${fx(xyz[3 * s[0] + 2], 4)}</span>`;
      }
      chip.innerHTML = html; chip.classList.remove('hidden');
    }
    const lines = [];
    if (s.length >= 2 && s.length <= 4) for (let k = 0; k + 1 < s.length; k++) lines.push([s[k], s[k + 1]]);
    R.setMeasureLines(lines);
    if (!light && d && ui) ui.refreshSelection(d);
    if (!light && d && d.reps && d.reps.some(r => r.visible !== false && /\bselect(ed|ion)\b/i.test(r.selection))) scheduleReps(d);
  }
  function onAtomClick(i, e) {
    const d = state.active; if (!d) return;
    if (ui && ui.atomClick(i, e)) return;
    const s = d.selection, k = s.indexOf(i);
    if (k >= 0) s.splice(k, 1);
    else { if (s.length >= 4) s.length = 0; s.push(i); }
    R.setSelection(s); updateMeasure();
  }
  function clearSelection() { const d = state.active; if (!d) return; d.selection = []; R.setSelection([]); updateMeasure(); }
  function selectAll() { const d = state.active; if (!d) return; d.selection = curFrame(d).Z.map((_, i) => i); R.setSelection(d.selection); updateMeasure(); status(d.selection.length + ' atoms selected'); }

  // ---------------- frames & animation ----------------
  /** Frames shown in the Steps list, playbar and profiles: the selected job's range, optionally converged points only. */
  function visibleFrames(d) {
    const fr = d.data.frames, [a, b] = jobRange(d);
    const all = [];
    for (let i = a; i <= b; i++) all.push(i);
    if (!d.convergedOnly) return all;
    const idx = all.filter(i => fr[i].converged);
    return idx.length ? idx : all;
  }
  function jobRange(d) { return BM.analysis.jobRange(d.data, d.job); }
  function setJob(d, j) {
    d.job = j;
    const [a, b] = jobRange(d);
    d._conv = null; d._prof = null;
    // renderPlot first: it replaces a plot type the new job does not have (e.g. an optimisation criterion),
    // which setFrame would otherwise try to draw
    renderTabs(); renderPanel(); renderPlot();
    if (d.frameIndex < a || d.frameIndex > b) setFrame(d, b);
    renderPlaybar();
  }
  function setFrame(d, fi, fromAnim) {
    if (fi < 0 || fi >= d.data.frames.length) return;
    if (!fromAnim) stopAnim(false);
    if (d.vibIndex >= 0 && fi !== d.data.vibrations.frameIndex) { d.vibIndex = -1; d._vibVec = null; d._animXyz = null; }
    const prevN = curFrame(d).Z.length;
    d.frameIndex = fi;
    if (curFrame(d).Z.length !== prevN) d.selection = [];
    refreshScene(true);
    highlightRow('#steps-table', fi);
    updatePlaybarFrame();
    if (framePlot(d.plotType)) renderPlotDraw();
    updateStatus();
    if (d.tab === 'build') { clearTimeout(buildTimer); if (fromAnim) buildTimer = setTimeout(renderPanel, 250); else renderPanel(); }
  }
  let buildTimer = 0;
  function selectState(d, i) {
    d.stateIndex = i;
    renderPanel();
    if (d.plotType !== 'uv' && d.plotType !== 'ecd' && plotTypes(d).some(p => p[0] === 'uv')) { d.plotType = 'uv'; renderPlot(); } else renderPlotDraw();
  }
  /** Up/Down (and Page Up/Down, Home/End) in the Steps, Vibrations and States lists move to another row, as Left/Right
   *  do in the plot (issue #4). False when the active tab has no such list, so the key keeps its usual meaning. */
  function listKey(d, key) {
    const tab = d.tab, updown = key === 'ArrowUp' || key === 'ArrowDown', back = key === 'ArrowUp' || key === 'PageUp' || key === 'Home';
    if (tab === 'steps' && d.data.frames.length > 1) {
      if (updown) { stepFrame(back ? -1 : 1); return true; }
      const vis = visibleFrames(d), k = Math.max(0, vis.indexOf(d.frameIndex));
      const to = key === 'PageUp' ? Math.max(0, k - 10) : key === 'PageDown' ? Math.min(vis.length - 1, k + 10) : key === 'Home' ? 0 : vis.length - 1;
      setFrame(d, vis[to]);
      return true;
    }
    const n = tab === 'vib' && d.data.vibrations ? d.data.vibrations.freqs.length
      : tab === 'states' && d.data.excitedStates ? d.data.excitedStates.length : 0;
    if (!n) return false;
    const cur = tab === 'vib' ? d.vibIndex : d.stateIndex;
    let k;
    if (updown) k = cur < 0 ? (back ? n - 1 : 0) : (cur + (back ? -1 : 1) + n) % n;
    else if (key === 'PageUp' || key === 'PageDown') k = Math.max(0, Math.min(n - 1, Math.max(0, cur) + (back ? -10 : 10)));
    else k = back ? 0 : n - 1;
    if (tab === 'vib') { selectVib(d, k); renderPanel(); } else selectState(d, k);
    return true;
  }
  function stepFrame(delta) {
    const d = state.active; if (!d) return;
    const vis = visibleFrames(d);
    let k = vis.indexOf(d.frameIndex);
    if (k < 0) k = 0;
    k = (k + delta + vis.length) % vis.length;
    setFrame(d, vis[k], state.anim.kind === 'frames');
  }
  function startAnim(kind) {
    state.anim.kind = kind; state.anim.last = 0; state.anim.acc = 0;
    if (!state.anim.raf) state.anim.raf = requestAnimationFrame(tick);
    renderPlaybar();
  }
  function stopAnim(resetVib) {
    const d = state.active;
    const wasVib = state.anim.kind === 'vib';
    state.anim.kind = null;
    if (state.anim.raf) { cancelAnimationFrame(state.anim.raf); state.anim.raf = 0; }
    if (d && resetVib && d.vibIndex >= 0) { d.vibIndex = -1; d._vibVec = null; d._animXyz = null; }
    if (d && wasVib) { R.updatePositions(curFrame(d).xyz); refreshArrows(); }
  }
  function tick(ts) {
    state.anim.raf = 0;
    const d = state.active;
    if (!d || !state.anim.kind) return;
    const dt = state.anim.last ? Math.min(0.1, (ts - state.anim.last) / 1000) : 0;
    state.anim.last = ts;
    if (state.anim.kind === 'vib' && d.vibIndex >= 0 && d._vibVec) {
      state.anim.phase += dt * 2 * Math.PI * state.vib.speed * 0.8;
      const eq = d.data.frames[d.data.vibrations.frameIndex].xyz, u = d._vibVec, x = d._animXyz, s = Math.sin(state.anim.phase);
      for (let i = 0; i < x.length; i++) x[i] = eq[i] + s * u[i];
      R.updatePositions(x);
      refreshOverlays(d);
      if (d.selection.length) updateMeasure(true);
    } else if (state.anim.kind === 'frames') {
      state.anim.acc += dt * state.fps;
      if (state.anim.acc >= 1) { state.anim.acc = 0; stepFrame(1); }
    }
    if (state.anim.kind) state.anim.raf = requestAnimationFrame(tick);
  }
  /** The displacement vector of vibration i, or null when it was not loaded (part of a huge output, issue #7). */
  const modeOf = (d, i) => { const v = d.data.vibrations; return v && i >= 0 && v.modes ? v.modes[i] || null : null; };
  const vibVec = (d, i, amp) => { const m = modeOf(d, i); return m ? geom.normalizedMode(m, amp) : null; };
  const NOT_LOADED = 'The displacements of this vibration were not loaded (Info tab: Load other parts of the file).';
  function selectVib(d, i, play) {
    const v = d.data.vibrations;
    if (!v || i < 0 || i >= v.freqs.length) return;
    stopAnim(false);
    if (d.frameIndex !== v.frameIndex) { d.frameIndex = v.frameIndex; }
    d.vibIndex = i;
    d._vibVec = vibVec(d, i, state.vib.amplitude);
    d._animXyz = Float64Array.from(d.data.frames[v.frameIndex].xyz);
    const still = d._vibVec && !d._vibVec.some(x => x !== 0);   // the file gave no displacements for it
    if (still) d._vibVec = null;
    if (!d._vibVec) {
      status(still ? 'The file has no atom displacements for this vibration, so it cannot be animated.' : NOT_LOADED);
      refreshScene(true); highlightRow('#vib-table', i); renderPlaybar();
      if (d.plotType === 'ir' || d.plotType === 'raman') renderPlotDraw();
      return;
    }
    state.anim.phase = 0;
    refreshScene(true);
    highlightRow('#vib-table', i);
    if (play !== false && state.vib.playing) startAnim('vib');
    else renderPlaybar();
    if (d.plotType === 'ir' || d.plotType === 'raman') renderPlotDraw();
    if (d.tab === 'vib') { const r = $('#vib-table tr.on'); if (r) scrollRowIntoView(r); }
  }
  /** A new amplitude for the animation, the two ends and the distances in the bond table. */
  function setAmplitude(d, amp) {
    if (!d || d.vibIndex < 0 || !d._vibVec) return;
    d._vibVec = vibVec(d, d.vibIndex, amp);
    refreshEnds();
    refreshVibBonds(d);
  }
  /** Pause the animation of the selected vibration at its centre (the geometry of the frequency job). */
  function holdVib(d) {
    state.vib.playing = false;   // also for the next vibration selected
    if (!d || d.vibIndex < 0) return;
    if (state.anim.kind === 'vib') stopAnim(false); else R.updatePositions(curFrame(d).xyz);
    state.anim.phase = 0;
    refreshOverlays(d);
    if (d.selection.length) updateMeasure(true);
    renderPlaybar();
  }
  function toggleVibPlay() {
    const d = state.active; if (!d || d.vibIndex < 0) return;
    if (state.anim.kind === 'vib') { state.vib.playing = false; state.anim.kind = null; renderPlaybar(); }
    else { state.vib.playing = true; startAnim('vib'); }
  }
  function stopVib() {
    const d = state.active; if (!d) return;
    stopAnim(true); refreshScene(true); renderPlaybar(); highlightRow('#vib-table', -1);
    if (d.plotType === 'ir' || d.plotType === 'raman') renderPlotDraw();
  }
  function togglePlay() {
    const d = state.active; if (!d) return;
    if (d.vibIndex >= 0) return toggleVibPlay();
    if (state.anim.kind === 'frames') { stopAnim(false); renderPlaybar(); }
    else if (d.data.frames.length > 1) startAnim('frames');
  }

  // ---------------- playbar ----------------
  function renderPlaybar() {
    const d = state.active, pb = $('#playbar');
    if (!d) { pb.classList.add('hidden'); return; }
    const v = d.data.vibrations;
    if (d.vibIndex >= 0 && v) {
      const f = v.freqs[d.vibIndex];
      const playing = state.anim.kind === 'vib';
      pb.innerHTML = `
        <button data-pb="play" title="Play / pause (Space)"><svg class="ic"><use href="#i-${playing ? 'pause' : 'play'}"/></svg></button>
        <span class="lbl">Mode ${d.vibIndex + 1} · <b class="mono">${f < 0 ? Math.abs(f).toFixed(1) + 'i' : f.toFixed(1)}</b> cm⁻¹${v.sym && v.sym[d.vibIndex] ? ' · ' + esc(v.sym[d.vibIndex]) : ''}</span>
        <span class="grp" title="Amplitude (largest atomic displacement)"><span class="lbl">Amp</span><input type="range" data-pb="amp" min="0.05" max="1.2" step="0.05" value="${state.vib.amplitude}"></span>
        <span class="grp" title="Speed"><span class="lbl">Speed</span><input type="range" data-pb="speed" min="0.2" max="3" step="0.1" value="${state.vib.speed}"></span>
        <button data-pb="arrows" class="${state.vib.arrows ? 'on' : ''}" title="Displacement arrows">↗</button>
        <button data-pb="stop" title="Stop (Esc)"><svg class="ic"><use href="#i-close"/></svg></button>`;
      pb.classList.remove('hidden');
      return;
    }
    const n = d.data.frames.length;
    if (n > 1) {
      const vis = visibleFrames(d);
      const k = Math.max(0, vis.indexOf(d.frameIndex));
      const playing = state.anim.kind === 'frames';
      pb.innerHTML = `
        <button data-pb="prev" title="Previous geometry (←)"><svg class="ic"><use href="#i-prev"/></svg></button>
        <button data-pb="play" title="Play / pause (Space)"><svg class="ic"><use href="#i-${playing ? 'pause' : 'play'}"/></svg></button>
        <button data-pb="next" title="Next geometry (→)"><svg class="ic"><use href="#i-next"/></svg></button>
        <input type="range" data-pb="frame" min="0" max="${vis.length - 1}" step="1" value="${k}">
        <span class="lbl" id="pb-label"></span>
        <span class="grp" title="Playback speed (frames per second)"><span class="lbl">fps</span><input type="range" data-pb="fps" min="1" max="30" step="1" value="${state.fps}"></span>`;
      pb.classList.remove('hidden');
      updatePlaybarFrame();
      return;
    }
    pb.classList.add('hidden');
  }
  function updatePlaybarFrame() {
    const d = state.active; if (!d) return;
    const lab = $('#pb-label'), sl = $('#playbar input[data-pb=frame]');
    const vis = visibleFrames(d), k = vis.indexOf(d.frameIndex);
    if (lab) lab.innerHTML = `Geometry <b class="mono">${d.frameIndex + 1}</b> / ${d.data.frames.length}${isNum(curFrame(d).energy) ? ' · <span class="mono">' + curFrame(d).energy.toFixed(6) + '</span> Eh' : inProgress(d, d.frameIndex) ? ' · in progress' : ''}`;
    if (sl && k >= 0) sl.value = k;
  }

  // ---------------- sidebar ----------------
  function renderDocs() {
    const box = $('#doc-list');
    box.innerHTML = state.docs.map(d => `
      <div class="doc ${d === state.active ? 'active' : ''}" data-doc="${d.id}" title="${esc(d.name)}">
        <span class="name">${esc(d.name)}</span><span class="badge">${esc(formatName(d.data))}</span>
        <button class="x" data-close="${d.id}" title="Close"><svg class="ic"><use href="#i-close"/></svg></button>
      </div>`).join('');
    const a = box.querySelector('.doc.active');
    if (a) { const br = box.getBoundingClientRect(), ar = a.getBoundingClientRect(); if (ar.top < br.top) box.scrollTop -= br.top - ar.top + 6; else if (ar.bottom > br.bottom) box.scrollTop += ar.bottom - br.bottom + 6; }
  }
  const FORMAT_SHORT = { 'gaussian-log': 'Gaussian', 'gaussian-input': 'Gaussian input', 'orca-out': 'ORCA', 'orca-input': 'ORCA input', xyz: 'XYZ', cube: 'Cube', molden: 'Molden', fchk: 'fchk', mol: 'MOL', pdb: 'PDB', gro: 'GRO', built: 'Built',
    'gamess-out': 'GAMESS', 'nwchem-out': 'NWChem', 'qchem-out': 'Q-Chem', 'psi4-out': 'Psi4', 'xtb-out': 'xtb', 'mopac-out': 'MOPAC', 'cp2k-out': 'CP2K', 'turbomole-coord': 'Turbomole', 'turbomole-gradient': 'Turbomole', 'turbomole-out': 'Turbomole',
    'molpro-out': 'Molpro', 'dalton-out': 'Dalton', 'gamess-punch': 'GAMESS punch',
    cif: 'CIF', mmcif: 'mmCIF', mol2: 'Mol2', cml: 'CML', cjson: 'CJSON', poscar: 'VASP', xdatcar: 'XDATCAR', outcar: 'OUTCAR', 'qe-in': 'QE input', 'qe-out': 'QE', 'castep-cell': 'CASTEP', xsf: 'XSF',
    prmtop: 'prmtop', inpcrd: 'inpcrd', psf: 'PSF', dcd: 'DCD', 'lammps-data': 'LAMMPS', 'lammps-dump': 'LAMMPS', trr: 'TRR', xtc: 'XTC' };
  function formatName(data) {
    const id = String(data.format || ''), m = FORMAT_SHORT[id];
    if (m) return m + (data.convertedBy ? ' (obabel)' : '');
    const f = (BM.formats || []).find(x => x.id === id);
    return f ? f.name.split(/[\s(]/)[0] : id.replace(/-(out|input|in)$/, '');
  }
  function tabsFor(d) {
    const t = [['info', 'Info']];
    if (!d.pending) t.push(['display', 'Display']);
    if (d.data.frames.length > 1) t.push(['steps', 'Steps', visibleFrames(d).length]);
    if (d.data.vibrations) t.push(['vib', 'Vibrations', d.data.vibrations.freqs.length]);
    if (d.data.excitedStates && d.data.excitedStates.length) t.push(['states', 'States', d.data.excitedStates.length]);
    if (d.surfaces.length || (!d.pending && curFrame(d).Z.length)) t.push(['surf', 'Surfaces', d.surfaces.length || 0]);
    if (!d.pending) t.push(['build', 'Build']);
    return t;
  }
  function renderTabs() {
    const d = state.active, nav = $('#tabs');
    if (!d) { nav.innerHTML = ''; $('#panel').innerHTML = ''; return; }
    const t = tabsFor(d);
    if (!t.some(x => x[0] === d.tab)) d.tab = 'info';
    nav.innerHTML = t.map(([id, name, c]) => `<button data-tab="${id}" class="${d.tab === id ? 'on' : ''}">${name}${c ? `<span class="count">${c}</span>` : ''}</button>`).join('');
  }
  function renderPanel() {
    const d = state.active, p = $('#panel');
    if (!d) { p.innerHTML = ''; return; }
    const fn = { info: panelInfo, display: panelDisplay, steps: panelSteps, vib: panelVib, states: panelStates, surf: panelSurf, build: ui.panelBuild }[d.tab] || panelInfo;
    const keep = p.scrollTop;
    p.innerHTML = fn(d);
    p.scrollTop = keep;
    if (d.tab === 'steps') highlightRow('#steps-table', d.frameIndex, true);
    if (d.tab === 'vib') highlightRow('#vib-table', d.vibIndex, true);
    if (d.tab === 'states') highlightRow('#states-table', d.stateIndex, true);
    if (d.tab === 'display') updateQueryInfo(d);
  }
  function highlightRow(sel, idx, scroll) {
    const t = $(sel); if (!t) return;
    $$('tr.row.on', t).forEach(r => r.classList.remove('on'));
    const r = t.querySelector(`tr.row[data-i="${idx}"]`);
    if (r) { r.classList.add('on'); if (scroll !== false) scrollRowIntoView(r); }
  }

  function scrollRowIntoView(r) {
    const p = $('#panel'), pr = p.getBoundingClientRect(), rr = r.getBoundingClientRect();
    if (rr.top < pr.top + 28) p.scrollTop -= (pr.top + 28 - rr.top);
    else if (rr.bottom > pr.bottom) p.scrollTop += rr.bottom - pr.bottom + 4;
  }
  function kv(rows) { const ok = rows.filter(r => r && r[1] !== null && r[1] !== undefined && r[1] !== ''); return ok.length ? '<div class="kv">' + ok.map(r => `<div class="k">${r[0]}</div><div class="v ${r[2] || ''}">${r[1]}</div>`).join('') + '</div>' : ''; }
  const sec = (title, body) => (body ? `<div class="sec"><h3>${title}</h3>${body}</div>` : '');
  function transWeight(d, t) {
    const unrestricted = /[AB]$/i.test(String(t.from)) || /[AB]$/i.test(String(t.to));
    return (d.data.format === 'gaussian-log' && !unrestricted ? 2 : 1) * t.coeff * t.coeff;
  }

  // ---------------- QM/MM layers taken from another document ----------------
  // Outputs that do not list their layers (Gaussian Geom=Check / AllCheck, ORCA without QMAtoms, xtb with the inner
  // region in a file, CP2K when the QM atoms cannot be matched) take them from an open document with the same atoms.
  const sameAtoms = (a, b) => a.length === b.length && a.every((z, i) => z === b[i]);
  /** Open documents whose layers fit d: the same number of atoms with the same elements in the same order. */
  function layerDonors(d) {
    const f = curFrame(d), N = f.Z.length;
    return state.docs.filter(o => o !== d && !o.pending && o.data.layers && o.data.layers.atomLayer && o.data.layers.atomLayer.length === N &&
      o.data.frames.some(fr => fr.Z.length === N && sameAtoms(fr.Z, f.Z)));
  }
  function layersSourceHtml(d) {
    const data = d.data, L = data.layers, m = data.meta || {};
    if (!m.qmmm || !L || d.pending) return '';
    const N = curFrame(d).Z.length, known = L.atomLayer && L.atomLayer.length === N;
    if (known && ['rebuilt', 'population', 'document'].indexOf(L.source) < 0) return '';
    const title = L.scheme === 'oniom' ? 'ONIOM layers' : L.scheme === 'qm-qm2' ? 'QM/QM2 layers' : 'QM/MM layers';
    const what = L.source === 'document' ? `Taken from <b>${esc(L.from)}</b>.`
      : L.source === 'population' ? 'Read from the population analysis of the high-level calculation (see the reader notes).'
      : L.source === 'rebuilt' ? 'Rebuilt from the link atoms (see the reader notes).' : 'Not listed in this file, so the layers are not shown.';
    const donors = layerDonors(d);
    const from = data.format === 'gaussian-log' ? 'the input (.gjf, .com) or output (.log) that started the series' + (m.oldchk ? ` (this job continued from ${esc(m.oldchk)})` : '') : 'the input or another output of the same system that lists them';
    const pick = donors.length
      ? `<div class="field"><span>From</span><select data-take-layers>${donors.map(o => `<option value="${o.id}">${esc(o.name)}</option>`).join('')}</select><button class="btn" data-take-layers-go title="Copy the layers (and residue names) of that document to this one">Use</button></div>`
      : '';
    const hint = !known || L.source === 'rebuilt' ? `<p class="note">${donors.length ? 'Or open' : 'Open'} ${from}, with the same ${N} atoms in the same order: its layers are then copied to this document.</p>` : '';
    const open = '<div class="btnrow" style="margin-top:6px"><button class="btn small" data-act="layers-open" title="Open a file with the same atoms: its layers are copied to this document">Open file with the layers…</button></div>';
    return `<div class="sec" id="sec-layers"><h3>${title}</h3><p class="note">${what}${donors.length && L.source !== 'document' ? ' Take them from an open document with the same atoms:' : ''}</p>${pick}${hint}${open}</div>`;
  }
  // "Open file with the layers…": the next document opened from a file (within ten minutes) gives d its layers
  let layersRequest = null;
  async function openLayersFile(d) {
    if (!d) return;
    const req = layersRequest = { doc: d, until: Date.now() + 600e3 };
    // the shell's own dialog returns once its files are open (or it was cancelled); the browser's file input later
    if ((await openDialog()) !== 'input' && layersRequest === req) layersRequest = null;
  }
  /** A document was opened from a file: when it was opened for another document's layers, copy them there. */
  function layersOpened(o) {
    const r = layersRequest;
    layersRequest = null;
    if (!r || !o || o === r.doc || Date.now() > r.until || state.docs.indexOf(r.doc) < 0) return;
    const d = r.doc;
    if (layerDonors(d).indexOf(o) < 0) {
      toast(o.data.layers && o.data.layers.atomLayer
        ? `${o.name} does not have the same ${curFrame(d).Z.length} atoms in the same order as ${d.name}: its layers were not copied.`
        : `${o.name} does not list QM/MM layers either: nothing was copied to ${d.name}.`, true);
      return;
    }
    activate(d);
    takeLayers(d, o);
  }
  /** Copy the layers of document o to d (and residue names when d has none); d keeps its own method labels, freeze
   *  codes, MM types and charges when it has them. */
  function takeLayers(d, o) {
    const f = curFrame(d), N = f.Z.length, src = o && o.data.layers;
    const sf = o && o.data.frames.find(fr => fr.Z.length === N && sameAtoms(fr.Z, f.Z));
    if (!src || !src.atomLayer || src.atomLayer.length !== N || !sf) { toast('The layers of ' + (o ? o.name : 'that document') + ' do not fit these atoms.', true); return; }
    const old = d.data.layers || {}, ol = old.labels || {}, sl = src.labels || {};
    const pick = (a, b) => (a !== null && a !== undefined ? a : b);
    d.data.layers = Object.assign({}, src, {
      scheme: old.scheme || src.scheme,
      labels: { H: pick(ol.H, sl.H), M: pick(ol.M, sl.M), L: pick(ol.L, sl.L) },
      nLayers: old.nLayers || src.nLayers,
      embedding: old.embedding || src.embedding,
      chargeMult: old.chargeMult && old.chargeMult.length ? old.chargeMult : src.chargeMult,
      frozen: old.frozen || src.frozen, active: old.frozen ? old.active : src.active, freezeCodes: old.freezeCodes || src.freezeCodes || null,
      mmTypes: old.mmTypes || src.mmTypes, mmCharges: old.mmCharges || src.mmCharges,
      linkAtoms: src.linkAtoms && src.linkAtoms.length ? src.linkAtoms : (old.linkAtoms || []),
      source: 'document', from: o.name, rebuilt: null, population: null
    });
    // the reader's note about the missing (or rebuilt) layers no longer applies
    d.data.warnings = (d.data.warnings || []).filter(w => !/take its layers in the Info tab/.test(w));
    if (sf.resNames && !f.resNames) {
      for (const fr of d.data.frames) if (fr.Z.length === N) { fr.atomNames = sf.atomNames || null; fr.resNames = sf.resNames; fr.resIds = sf.resIds || null; if (sf.chains) fr.chains = sf.chains; }
    }
    d.bonds.clear(); d._rc = null; d.reps = null;
    refreshScene(true); renderTabs(); renderPanel(); updateStatus();
    const cnt = {};
    for (const x of d.data.layers.atomLayer) cnt[x] = (cnt[x] || 0) + 1;
    const nm = d.data.layers.scheme === 'oniom' ? { H: 'high', M: 'medium', L: 'low' } : { H: 'QM', M: 'QM2', L: 'MM' };
    toast(`Layers taken from ${o.name}: ` + ['H', 'M', 'L'].filter(k => cnt[k]).map(k => `${cnt[k]} ${nm[k]}`).join(', ') +
      (d.data.layers.linkAtoms.length ? `, ${d.data.layers.linkAtoms.length} link atom${d.data.layers.linkAtoms.length === 1 ? '' : 's'}` : ''));
  }

  function panelInfo(d) {
    const data = d.data, m = data.meta || {}, f = curFrame(d);
    const nAtoms = f.Z.filter(z => z > 0).length;
    const lastE = [...data.frames].reverse().find(x => isNum(x.energy));
    let term = null;
    if (m.normalTermination === true) term = '<span class="pill ok">Normal termination</span>';
    else if (m.normalTermination === false) term = '<span class="pill bad">Error / incomplete</span>';
    else if (m.running) term = '<span class="pill warn" title="No termination message yet: the job is still running, or it was stopped">Not finished: still running?</span>';
    const vib = data.vibrations;
    let vibTxt = null;
    if (vib) {
      const ni = vib.freqs.filter(x => x < 0).length;
      vibTxt = `${vib.freqs.length} modes · ` + (ni ? `<span class="pill warn">${ni} imaginary</span>` : '<span class="pill ok">no imaginary</span>');
    }
    const conv = data.frames.some(x => x.converged) ? '<span class="pill ok">converged</span>' : (data.frames.length > 1 && (m.jobTypes || []).some(j => /opt/i.test(j)) ? '<span class="pill warn">not converged</span>' : null);
    let html = contentsHtml(d);
    if (d.pending) html += `<div class="sec"><div class="route pending-note">Coordinates pending: open a coordinate or trajectory file with ${f.Z.length} atoms (inpcrd, NetCDF, DCD, TRR, XTC, PDB, GRO…).</div></div>`;
    html += sec('File', kv([
      ['Name', esc(d.name)], ['Format', esc(formatName(data))], ['Program', esc(data.program)], ['Title', esc(data.title)],
      ['Downloaded', data.downloaded ? `<span title="${esc(data.downloaded.url)}">${esc(data.downloaded.from)} · ${esc(data.downloaded.date)}</span>` : null],
      ['Status', term], ['Jobs', esc((m.jobTypes || []).join(', '))]]));
    html += loadedHtml(d);
    const route = m.route || m.keywords;
    if (route) html += `<div class="sec"><h3>${m.route ? 'Route' : 'Keywords'}</h3><div class="route">${esc(route)}</div></div>`;
    html += `<div class="sec"><h3>Molecule</h3>${kv([
      ['Formula', esc(geom.formula(f.Z))], ['Atoms', String(nAtoms) + (nAtoms !== f.Z.length ? ` (+${f.Z.length - nAtoms} dummy)` : '')],
      ['Mol. weight', geom.molWeight(f.Z).toFixed(3) + ' g/mol'], ['Method', esc(m.method)], ['Basis', esc(m.basis)], ['Solvent', esc(m.solvent)],
      ['Charge', isNum(m.charge) ? String(m.charge) : null], ['Multiplicity', isNum(m.multiplicity) ? String(m.multiplicity) : null]])}</div>`;
    // QM/MM and ONIOM: the components of this geometry's total energy (issue #10)
    const parts = m.qmmm && f.energies && isNum(f.energy) ? BM.analysis.energyComponents(data, [f], 1) : [];
    // the other levels the job computed for this geometry (SCF, MP2, CCSD..., an excited state), and which one is the
    // geometry's energy (the highest level: docs/analysing-outputs.md)
    const lev = !m.qmmm && f.energies && isNum(f.energy) ? Object.keys(f.energies).filter(k => isNum(f.energies[k])) : [];
    const same = k => Math.abs(f.energies[k] - f.energy) < 1e-10;
    const other = lev.filter(k => !same(k)), used = lev.filter(same).sort((a, b) => (a === 'Final') - (b === 'Final'))[0];
    const levName = k => (k === 'TD' ? 'Excited state' : BM.analysis.componentLabel(k));   // (Gaussian and fchk call it TD)
    html += sec('Energy', kv([
      ['This geometry', isNum(f.energy) ? `<span class="mono">${f.energy.toFixed(8)}</span> Eh${other.length && used ? ` <span class="dim">(${esc(levName(used))})</span>` : ''}` : null],
      ...other.map(k => ['<span class="dim">' + esc(levName(k)) + '</span>', `<span class="mono">${f.energies[k].toFixed(8)}</span> Eh`]),
      ...parts.map(cp => ['<span class="dim">' + esc(cp.label) + '</span>', `<span class="mono">${cp.values[0].toFixed(8)}</span> Eh`]),
      ['Final', lastE ? `<span class="mono">${lastE.energy.toFixed(8)}</span> Eh` : null],
      ['Geometries', data.frames.length > 1 ? `${data.frames.length} ${conv || ''}` : null],
      ['Vibrations', vibTxt], ['Excited states', data.excitedStates && data.excitedStates.length ? String(data.excitedStates.length) : null]]));
    const orb = data.orbitals;
    if (orb && orb.alpha && orb.alpha.energies && orb.alpha.nOcc > 0) {
      const e = orb.alpha.energies, h = e[orb.alpha.nOcc - 1], l = e[orb.alpha.nOcc];
      const rows = [['HOMO', `<span class="mono">${fx(h, 5)}</span> Eh (${fx(h * HARTREE.ev, 3)} eV)`]];
      if (isNum(l)) rows.push(['LUMO', `<span class="mono">${fx(l, 5)}</span> Eh (${fx(l * HARTREE.ev, 3)} eV)`], ['Gap', `<b>${fx((l - h) * HARTREE.ev, 3)}</b> eV`]);
      if (orb.beta && orb.beta.energies && orb.beta.nOcc > 0) {
        const hb = orb.beta.energies[orb.beta.nOcc - 1], lb = orb.beta.energies[orb.beta.nOcc];
        rows.push(['β HOMO / LUMO', `${fx(hb * HARTREE.ev, 3)} / ${fx(lb * HARTREE.ev, 3)} eV`]);
      }
      html += `<div class="sec" id="sec-orb"><h3>Orbitals</h3>${kv(rows)}</div>`;
    }
    if (data.dipole) html += `<div class="sec" id="sec-dipole"><h3>Dipole moment</h3>${kv([['Total', `<b class="mono">${fx(data.dipole.total, 4)}</b> D`], ['x, y, z', `<span class="mono">${fx(data.dipole.x, 4)}, ${fx(data.dipole.y, 4)}, ${fx(data.dipole.z, 4)}</span>`]])}</div>`;
    const t = data.thermo;
    if (t) html += `<div class="sec" id="sec-thermo"><h3>Thermochemistry${isNum(t.temperature) ? ` (${t.temperature.toFixed(2)} K)` : ''}</h3>${kv([
      ['ZPE', isNum(t.zpe) ? `<span class="mono">${t.zpe.toFixed(6)}</span> Eh (${(t.zpe * HARTREE.kcal).toFixed(2)} kcal/mol)` : null],
      ['E + ZPE', isNum(t.eZPE) ? `<span class="mono">${t.eZPE.toFixed(6)}</span> Eh` : null],
      ['Enthalpy H', isNum(t.H) ? `<span class="mono">${t.H.toFixed(6)}</span> Eh` : null],
      ['Gibbs G', isNum(t.G) ? `<b class="mono">${t.G.toFixed(6)}</b> Eh` : null]])}</div>`;
    if (data.charges) {
      const schemes = Object.keys(data.charges).filter(k => Array.isArray(data.charges[k]));
      if (schemes.length) {
        const q = data.charges[d.chargeScheme] || data.charges[schemes[0]];
        const Z = data.frames[data.frames.length - 1].Z;
        html += `<div class="sec" id="sec-charges"><h3>Atomic charges</h3><div class="field"><span>Scheme</span><select data-charge-scheme>${schemes.map(s => `<option ${s === d.chargeScheme ? 'selected' : ''}>${s}</option>`).join('')}</select><span></span></div>
          <table class="list" style="margin-top:8px"><tr><th>Atom</th><th>Charge</th>${data.spins && data.spins.Mulliken ? '<th>Spin</th>' : ''}</tr>
          ${q.slice(0, 600).map((v, i) => `<tr class="row" data-atom="${i}"><td>${el.symbol(Z[i] || 0)}${i + 1}</td><td>${fx(v, 4)}</td>${data.spins && data.spins.Mulliken ? `<td>${fx(data.spins.Mulliken[i], 4)}</td>` : ''}</tr>`).join('')}</table>
          </div>`;
      }
    }
    if (data.cell) html += cellHtml(data.cell);
    html += layersSourceHtml(d);
    if (data.warnings && data.warnings.length) html += `<div class="sec" id="sec-warn"><h3>Reader notes</h3><ul class="warnings">${data.warnings.map(w => `<li>${esc(w)}</li>`).join('')}</ul></div>`;
    return html;
  }

  /** What was loaded of a big file (issue #7), with "Load other parts…" when the file can be read again. */
  function loadedHtml(d) {
    const L = d.data.meta && d.data.meta.loaded, T = d.trajLoad, rows = [];
    const how = s => (s.mode === 'every' ? 'every ' + ordinal(s.n) : s.mode === 'last' ? 'the last ' + s.n : s.mode === 'range' ? s.from + '–' + (s.to || '…') : '');
    if (L && L.steps) {
      const s = L.steps;
      rows.push([s.frames ? 'Frames' : 'Geometries', `${s.loaded.toLocaleString()} of ${s.total.toLocaleString()}` + (how(s) ? ' (' + how(s) + ')' : '') +
        (!s.frames && s.loaded < s.total ? '; the plots show the energies of every step' : '')]);
    }
    if (L && L.modes) rows.push(['Displacements', L.modes.loaded >= L.modes.total ? `all ${L.modes.total.toLocaleString()} vibrations` : `${L.modes.loaded.toLocaleString()} of ${L.modes.total.toLocaleString()} vibrations`]);
    if (L && L.bytes) rows.push(['Memory', 'about ' + MBs1(L.bytes) + ' in number arrays']);
    if (T) rows.push(['Frames of ' + esc(T.name), `${T.count.toLocaleString()} of ${T.total.toLocaleString()}` + (T.stride > 1 ? ` (every ${ordinal(T.stride)})` : T.count < T.total ? ` (${T.first + 1}–${T.last + 1})` : '')]);
    if (!rows.length) return '';
    const can = !!(d.loadSource || (T && T.reload));
    return `<div class="sec" id="sec-loaded"><h3>Loaded part of the file</h3>${kv(rows)}${can ? '<div class="btnrow" style="margin-top:6px"><button class="btn small" data-act="reload-part" title="Choose again what to load and read the file again">Load other parts…</button></div>' : ''}</div>`;
  }

  /** Everything the file holds, with counts; each entry jumps to its tab or plot. */
  function contentsItems(d) {
    const data = d.data, it = [], types = plotTypes(d), has = id => types.some(t => t[0] === id);
    const n = data.frames.length, c = convFor(d), prof = profFor(d);
    if (n > 1 || c.source === 'optSteps') {
      const sc = prof.find(p => p.kind === 'scan'), ir = prof.find(p => p.kind === 'irc');
      const opt = c.isOpt && c.steps.length > 1;
      const sub = [sc ? sc.points.length + ' scan points' : '', ir ? ir.points.length + ' IRC points' : '', opt ? c.steps.length + (c.source === 'optSteps' ? ' opt cycles' : ' opt steps') : ''].filter(Boolean).join(' · ');
      it.push({ name: 'Geometries', count: n, sub, tab: n > 1 ? 'steps' : '', plot: sc ? 'scan' : ir ? 'irc' : opt ? 'opt:dE' : has('energy') ? 'energy' : '' });
    }
    if (data.jobs && data.jobs.length > 1) it.push({ name: 'Jobs', count: data.jobs.length, sub: data.jobs.map(j => (j.kinds || []).join('+') || 'job').join(', '), anchor: 'ct-job' });
    const v = data.vibrations;
    if (v) { const ni = v.freqs.filter(x => x < 0).length; it.push({ name: 'Vibrations', count: v.freqs.length, sub: ni ? ni + ' imaginary' : 'no imaginary', warn: ni > 0, tab: 'vib', plot: has('ir') ? 'ir' : has('raman') ? 'raman' : '' }); }
    const S = data.excitedStates;
    if (S && S.length) it.push({ name: 'Excited states', count: S.length, sub: has('ecd') ? 'ECD available' : '', tab: 'states', plot: has('uv') ? 'uv' : has('ecd') ? 'ecd' : '' });
    if (data.nmr && data.nmr.iso) { const nc = (data.nmr.couplings || []).length; it.push({ name: 'NMR shieldings', count: data.nmr.iso.filter(isNum).length, sub: nc ? nc + ' couplings' : data.nmr.method || '', plot: has('nmr') ? 'nmr' : '' }); }
    const sch = data.charges ? Object.keys(data.charges).filter(k => Array.isArray(data.charges[k])) : [];
    if (sch.length) it.push({ name: 'Charges', count: sch.length, sub: sch.join(', '), plot: has('charges') ? 'charges' : '', anchor: 'sec-charges' });
    const mo = BM.analysis.moLevels(data);
    if (mo.length) {
      const a = mo[0], gap = a.homo >= 0 && a.lumo >= 0 ? (a.levels[a.lumo].e - a.levels[a.homo].e) * HARTREE.ev : null;
      it.push({ name: 'Orbital energies', count: mo.reduce((x, c2) => x + c2.levels.length, 0), sub: (mo.length > 1 ? 'α + β' : '') + (isNum(gap) ? (mo.length > 1 ? ' · ' : '') + 'gap ' + gap.toFixed(2) + ' eV' : ''), plot: 'mo', anchor: 'sec-orb' });
    }
    const bm = data.basis && data.mos && BM.orbitals ? BM.orbitals.model(data) : null;
    if (bm) it.push({ name: 'Basis set and orbitals', count: bm.nbf, sub: 'functions · orbitals and densities can be drawn', tab: 'surf' });
    if (d.surfaces.length) it.push({ name: 'Volumes', count: d.surfaces.length, sub: [...new Set(d.surfaces.map(x => KIND_NAMES[x.vol.kind] || x.vol.kind))].join(', '), tab: 'surf' });
    const T = data.thermo;
    if (T) it.push({ name: 'Thermochemistry', count: '', sub: isNum(T.G) ? 'G = ' + T.G.toFixed(5) : isNum(T.zpe) ? 'ZPE' : '', anchor: 'sec-thermo' });
    if (data.dipole) it.push({ name: 'Dipole', count: '', sub: isNum(data.dipole.total) ? data.dipole.total.toFixed(3) + ' D' : '', anchor: 'sec-dipole' });
    if (data.cell) it.push({ name: 'Unit cell', count: '', sub: data.cell.map(c2 => Math.hypot(...c2).toFixed(2)).join(' × ') + ' Å', anchor: 'sec-cell' });
    if (data.warnings && data.warnings.length) it.push({ name: 'Warnings', count: data.warnings.length, warn: true, anchor: 'sec-warn' });
    // QM/MM layers and residue classes: a click applies the matching preset or selects the atoms
    const f = curFrame(d);
    if (!d.pending && f.Z.length) {
      const c = repCtx(d), L = c.layers, DL = data.layers;
      if (L) {
        const cnt = {};
        for (const x of L.atomLayer) cnt[x] = (cnt[x] || 0) + 1;
        const nm = L.scheme === 'oniom' ? { H: 'High layer', M: 'Medium layer', L: 'Low layer' } : { H: 'QM region', M: 'QM2 region', L: 'MM region' };
        for (const k of ['H', 'M', 'L']) if (cnt[k]) it.push({ name: nm[k], count: cnt[k], sub: (L.labels && L.labels[k]) || 'layer ' + k, preset: 'qmmm', select: 'layer ' + k });
      } else if (data.meta && data.meta.qmmm && DL) {
        // a QM/MM file that does not list its layers: the way to take them from another file
        it.push({ name: DL.scheme === 'oniom' ? 'ONIOM layers' : 'QM/MM layers', count: '', sub: 'not listed in this file: take them from another', warn: true, anchor: 'sec-layers' });
      }
      if (f.resNames) {
        const cc = REPS.contents(c), cl = cc.classes;
        const add = (key, name, preset, select, unit) => { const v = cl[key]; if (v.res) it.push({ name, count: v.res, sub: (unit ? v.res + ' ' + (v.res === 1 ? unit.replace(/s$/, '') : unit) + ' · ' : '') + v.atoms.toLocaleString() + ' atom' + (v.atoms === 1 ? '' : 's'), preset, select }); };
        add('protein', 'Protein', 'protein', null, cc.chains > 1 ? cc.chains + ' chains · residues' : 'residues');
        add('nucleic', 'Nucleic acid', 'nucleic', null, 'residues');
        add('lipid', 'Lipids', 'membrane', null, 'residues');
        add('carbohydrate', 'Glycans', 'glycan', null, 'sugar residues');
        add('ligand', 'Ligands', null, 'ligand', 'groups');
        add('water', 'Water', null, 'water', 'molecules');
        add('ions', 'Ions', null, 'ions', '');
      }
    }
    return it;
  }
  function contentsHtml(d) {
    const it = contentsItems(d), jobs = d.data.jobs;
    if (!it.length) return '';
    const chips = it.map(x => `<button class="ct ${x.warn ? 'warn' : ''}" data-go-tab="${x.tab || ''}" data-go-plot="${x.plot || ''}" data-go-anchor="${x.anchor || ''}" data-go-preset="${x.preset || ''}" data-go-select="${esc(x.select || '')}" title="${esc(x.name + (x.sub ? ': ' + x.sub : '') + (x.preset ? ' — click for the ' + (REPS.PRESETS.find(p => p.id === x.preset) || {}).label + ' preset' : x.select ? ' — click to select' : ''))}"><span class="ct-n">${esc(x.name)}</span><span class="ct-c mono">${x.count === '' ? '' : x.count}</span><span class="ct-s">${esc(x.sub || '')}</span></button>`).join('');
    return `<div class="sec"><h3>Contents</h3><div class="contents">${chips}</div>${jobs && jobs.length > 1 ? jobSelect(d, 'ct-job') : ''}</div>`;
  }
  function jobSelect(d, id) {
    const jobs = d.data.jobs;
    return `<label class="field ct-jobsel" ${id ? `id="${id}"` : ''}><span>Job</span><select data-job title="Restrict steps and plots to one job"><option value="-1">All jobs (${d.data.frames.length} geometries)</option>${jobs.map((j, k) => `<option value="${k}" ${d.job === k ? 'selected' : ''}>${k + 1} · ${esc((j.kinds || []).join('+') || 'job')} · ${j.frameEnd - j.frameStart + 1} geom.${j.normalTermination === false ? ' · error' : ''}</option>`).join('')}</select><span></span></label>`;
  }
  function cellHtml(c) {
    const L = c.map(v => Math.hypot(v[0], v[1], v[2]));
    const ang = (u, v, a, b) => Math.acos(Math.max(-1, Math.min(1, (u[0] * v[0] + u[1] * v[1] + u[2] * v[2]) / (a * b)))) * 180 / Math.PI;
    const vol = Math.abs(c[0][0] * (c[1][1] * c[2][2] - c[1][2] * c[2][1]) - c[0][1] * (c[1][0] * c[2][2] - c[1][2] * c[2][0]) + c[0][2] * (c[1][0] * c[2][1] - c[1][1] * c[2][0]));
    return `<div class="sec" id="sec-cell"><h3>Unit cell</h3>${kv([['a, b, c', `<span class="mono">${L.map(x => x.toFixed(4)).join(', ')}</span> Å`],
      ['α, β, γ', `<span class="mono">${[ang(c[1], c[2], L[1], L[2]), ang(c[0], c[2], L[0], L[2]), ang(c[0], c[1], L[0], L[1])].map(x => x.toFixed(2)).join(', ')}</span>°`],
      ['Volume', `<span class="mono">${vol.toFixed(3)}</span> Å³`]])}</div>`;
  }
  function goTo(d, tab, plotId, anchor) {
    if (plotId && plotTypes(d).some(t => t[0] === plotId)) {
      d.plotType = plotId;
      if (state.plotCollapsed) { state.plotCollapsed = false; $('#center').classList.remove('plot-collapsed'); }
      renderPlot();
    }
    if (tab && tabsFor(d).some(t => t[0] === tab)) { d.tab = tab; renderTabs(); renderPanel(); }
    else if (anchor) { const e = document.getElementById(anchor); if (e) e.scrollIntoView({ block: 'start' }); }
  }

  function energyIn(d, e, ref) { return (e - ref) * HARTREE[d.spec.energy.unit]; }
  function energyRef(d) {
    const E = d.data.frames.map(f => f.energy).filter(isNum);
    if (!E.length) return 0;
    if (d.spec.energy.ref === 'first') return E[0];
    if (d.spec.energy.ref === 'last') return E[E.length - 1];
    return Math.min(...E);
  }
  /** The last geometry of a job that is still running, whose energy is not printed yet. */
  function inProgress(d, i) {
    const fr = d.data.frames;
    return !!(d.data.meta && d.data.meta.running) && i === fr.length - 1 && fr.length > 1 && !isNum(fr[i].energy);
  }
  function panelSteps(d) {
    const fr = d.data.frames, ref = energyRef(d), u = UNIT_LABEL[d.spec.energy.unit];
    const vis = visibleFrames(d);
    const anyConv = fr.some(f => f.converged);
    const rows = vis.slice(0, 5000).map(i => {
      const f = fr[i];
      const busy = inProgress(d, i);
      return `<tr class="row" data-i="${i}"><td>${i + 1}${f.label ? ` <span class="lab" title="${esc(f.label)}">${esc(f.label)}</span>` : ''}${busy ? ' <span class="lab busy" title="The job is still running: this step has no energy yet">in progress</span>' : ''}</td><td>${fx(f.energy, 6)}</td><td>${isNum(f.energy) ? energyIn(d, f.energy, ref).toFixed(2) : '—'}</td><td>${f.converged ? '<span class="ok-mark">●</span>' : ''}</td></tr>`;
    }).join('');
    return `<div class="sec"><div class="controls">
        ${d.data.jobs && d.data.jobs.length > 1 ? jobSelect(d) : ''}
        <div class="field"><span>ΔE unit</span><select data-energy="unit">${Object.keys(UNIT_LABEL).map(k => `<option value="${k}" ${k === d.spec.energy.unit ? 'selected' : ''}>${UNIT_LABEL[k]}</option>`).join('')}</select><span></span></div>
        <div class="field"><span>Relative to</span><select data-energy="ref"><option value="min" ${d.spec.energy.ref === 'min' ? 'selected' : ''}>lowest energy</option><option value="first" ${d.spec.energy.ref === 'first' ? 'selected' : ''}>first geometry</option><option value="last" ${d.spec.energy.ref === 'last' ? 'selected' : ''}>last geometry</option></select><span></span></div>
        ${anyConv ? `<label class="field"><span>Filter</span><span><input type="checkbox" data-conv-only ${d.convergedOnly ? 'checked' : ''}> converged points only (scans, IRC)</span><span></span></label>` : ''}
      </div>
      <table class="list" id="steps-table"><tr><th>#</th><th>E (Eh)</th><th>ΔE (${u})</th><th>${anyConv ? 'conv.' : ''}</th></tr>${rows}</table>
      ${vis.length > 5000 ? `<p class="note">Showing the first 5000 of ${vis.length} geometries; use the slider to reach the rest.</p>` : ''}
      ${opDetails(d)}</div>`;
  }
  const CRIT_NAMES = { maxForce: 'Max force', rmsForce: 'RMS force', maxDisp: 'Max displacement', rmsDisp: 'RMS displacement', stepNorm: 'Step length', energyChange: 'Energy change', gradNorm: 'Gradient norm' };
  function opDetails(d) {
    const f = curFrame(d);
    if (!f.opt) return '';
    const o = f.opt, rows = [];
    for (const k of Object.keys(CRIT_NAMES)) if (isNum(o[k])) {
      const thr = o[k + 'Thr'], cv = o[k + 'Conv'];
      rows.push([CRIT_NAMES[k], `<span class="mono">${o[k].toExponential(2)}</span>${isNum(thr) ? ` <span class="dim mono">/ ${thr.toExponential(1)}</span>` : ''}${cv === true ? ' <span class="ok-mark">✓</span>' : cv === false ? ' <span class="no-mark">✗</span>' : ''}`]);
    }
    if (isNum(o.negEig)) rows.push(['Negative Hessian eigenvalues', `<span class="mono">${o.negEig}</span>${tsJob(d) ? (o.negEig === 1 ? ' <span class="ok-mark">✓</span>' : ' <span class="no-mark">✗</span> <span class="dim">(1 for a transition-state search)</span>') : ''}`]);
    return rows.length ? `<div class="sec" style="margin-top:14px"><h3>Convergence at geometry ${d.frameIndex + 1}</h3>${kv(rows)}<p class="note"><a href="#" data-go-plot="opt:all" data-go-tab="" data-go-anchor="">Plot all criteria / threshold</a></p></div>` : '';
  }

  function panelVib(d) {
    const v = d.data.vibrations;
    const hasIR = !!v.ir, hasR = !!v.raman, hasSym = v.sym && v.sym.some(Boolean);
    // huge outputs (QM/MM frequencies): a window of rows around the selected mode
    const nV = v.freqs.length, W = 2000, lo = nV <= W ? 0 : Math.max(0, Math.min(nV - W, (d.vibIndex < 0 ? 0 : d.vibIndex) - W / 2)), hi = Math.min(nV, lo + W);
    const partial = v.modesLoaded !== undefined && v.modesLoaded < nV;
    const hasMode = i => (v.modes && v.modes.has ? v.modes.has(i) : !!v.modes[i]);
    const rows = v.freqs.slice(lo, hi).map((f, i0) => { const i = lo + i0; return `<tr class="row${partial && !hasMode(i) ? ' nodisp' : ''}" data-i="${i}"${partial && !hasMode(i) ? ' title="Displacements not loaded"' : ''}><td>${i + 1}</td><td class="${f < 0 ? 'imag' : ''}">${f < 0 ? Math.abs(f).toFixed(2) + 'i' : f.toFixed(2)}</td>${hasSym ? `<td class="dim">${esc(v.sym[i] || '')}</td>` : ''}${hasIR ? `<td>${fx(v.ir[i], 2)}</td>` : ''}${hasR ? `<td>${fx(v.raman[i], 2)}</td>` : ''}</tr>`; }).join('');
    const st = BM.vibtools.stationary(v), desc = BM.vibtools.describe(v, st, tsJob(d)), top = st.imag.length ? st.imag[0] : -1;
    const show = top >= 0 && d.vibIndex !== top ? ` <a href="#" data-act="vib-imag">Show ${st.imag.length > 1 ? 'the largest' : 'it'}</a>` : '';
    const off = d.vibIndex < 0 ? 'disabled' : '';
    return `<div class="sec">
      <p class="vib-summary ${desc.ok ? 'ok' : 'warn'}" id="vib-summary">${esc(desc.text)}${show}</p>
      <div class="controls">
        <label class="field" title="Largest atom shift in the animation and at the two ends"><span>Amplitude</span><input type="range" data-vib="amplitude" min="0.05" max="1.2" step="0.05" value="${state.vib.amplitude}"><output>${state.vib.amplitude.toFixed(2)} Å</output></label>
        <label class="field"><span>Speed</span><input type="range" data-vib="speed" min="0.2" max="3" step="0.1" value="${state.vib.speed}"><output>${state.vib.speed.toFixed(1)}×</output></label>
        <label class="field"><span>Arrows</span><span><input type="checkbox" data-vib="arrows" ${state.vib.arrows ? 'checked' : ''}> show displacement vectors</span><span></span></label>
        <label class="field" title="Length of the longest arrow"><span>Arrow length</span><input type="range" data-vib="arrowLen" min="0.2" max="3" step="0.1" value="${state.vib.arrowLen}"><output>${state.vib.arrowLen.toFixed(1)} Å</output></label>
        <label class="field" title="Both ends of the motion as translucent structures, at the amplitude above"><span>Both ends</span><span><input type="checkbox" data-vib="ends" ${state.vib.ends ? 'checked' : ''}> <span class="end-key minus">−</span> and <span class="end-key plus">+</span> ends, still</span><span></span></label>
      </div>
      <table class="list" id="vib-table"><tr><th>#</th><th>ν (cm⁻¹)</th>${hasSym ? '<th>Sym.</th>' : ''}${hasIR ? '<th>IR (' + esc(v.irUnit || 'km/mol') + ')</th>' : ''}${hasR ? '<th>Raman (Å⁴/amu)</th>' : ''}</tr>${rows}</table>
      <p class="note">Click a mode to animate it, or click a peak in the spectrum below. Imaginary modes are shown in orange.</p>
      </div>
      <div class="sec" id="vib-bonds">${vibBondsHtml(d)}</div>
      <div class="sec"><h3>Displace along the selected mode</h3>
        <div class="inline"><input type="number" id="disp-amt" value="${state.vib.disp || 0.15}" step="0.05" min="0.01" style="width:90px"> <span class="note">Å (largest atom shift)</span></div>
        <div class="btnrow" style="margin-top:8px">
          <button class="btn small" data-act="displace-minus" ${off}>− direction</button>
          <button class="btn small" data-act="displace-plus" ${off}>+ direction</button>
          <button class="btn small" data-act="export-vib" ${off}>Save animation (.xyz)</button>
          <button class="btn small" data-act="vib-inputs" ${off} title="Two calculation inputs, displaced by − and + this amount, to optimise towards both minima">Forward and reverse inputs…</button>
        </div>
        <p class="note">− / + direction creates a new “Edited” geometry, e.g. to push a structure off a saddle point. Forward and reverse inputs writes one optimisation input for each direction (<i>name</i>_minus and <i>name</i>_plus).</p>
      </div>`;
  }
  /** Is the document's job (or one of its jobs) a transition-state search? */
  const tsJob = d => ((d.data.meta && d.data.meta.jobTypes) || []).indexOf('ts') >= 0 || (d.data.jobs || []).some(j => (j.kinds || []).indexOf('ts') >= 0);
  /** Atom pairs whose distance changes most along the selected vibration, with the distances at both ends. */
  function vibBondsHtml(d) {
    const v = d.data.vibrations;
    if (!v || d.vibIndex < 0 || !d._vibVec) return '';
    const eq = d.data.frames[v.frameIndex], f = v.freqs[d.vibIndex], imag = f < 0;
    const rows = BM.vibtools.bondChanges(eq.Z, eq.xyz, d._vibVec, { imaginary: imag });
    if (!rows.length) return '';
    const lab = i => el.symbol(eq.Z[i]) + (i + 1), sel = d.selection.length === 2 ? d.selection.slice().sort((a, b) => a - b).join(',') : '';
    const body = rows.map(r => `<tr class="row${sel === r.i + ',' + r.j ? ' on' : ''}" data-pair="${r.i},${r.j}"><td>${lab(r.i)}–${lab(r.j)}</td><td class="mono">${r.dMinus.toFixed(3)}</td><td class="mono">${r.d0.toFixed(3)}</td><td class="mono">${r.dPlus.toFixed(3)}</td><td class="chg ${r.change}">${r.change}</td></tr>`).join('');
    return `<h3>Bonds along mode ${d.vibIndex + 1}</h3>
      <table class="list" id="vib-bonds-table"><tr><th>Atoms</th><th>− end</th><th>Here</th><th>+ end</th><th>Toward +</th></tr>${body}</table>
      <p class="note">Distances in Å at both ends of the motion (±${state.vib.amplitude.toFixed(2)} Å, the amplitude) and at this geometry, largest change first.${imag ? ' Along the imaginary mode of a transition state, the bonds that form and break are those of the reaction.' : ''} Click a row to measure it.</p>`;
  }
  function refreshVibBonds(d) { const box = $('#vib-bonds'); if (box && d && d.tab === 'vib') box.innerHTML = vibBondsHtml(d); }

  function panelStates(d) {
    const S = d.data.excitedStates;
    const rows = S.map((s, i) => `<tr class="row" data-i="${i}"><td>${s.index}</td><td>${fx(s.energyEV, 4)}</td><td>${fx(s.wavelengthNM, 2)}</td><td>${s.f === null ? '—' : fx(s.f, 4)}</td><td class="dim">${esc(s.label || '')}</td></tr>`).join('');
    const cur = d.stateIndex >= 0 ? S[d.stateIndex] : null;
    let det = '';
    if (cur) {
      const tr = (cur.transitions || []).slice().sort((a, b) => Math.abs(b.coeff) - Math.abs(a.coeff));
      det = `<div class="sec" style="margin-top:12px"><h3>State ${cur.index}: ${fx(cur.energyEV, 3)} eV · ${fx(cur.wavelengthNM, 1)} nm</h3>
        ${tr.length ? tr.slice(0, 20).map(t => `<div class="trans">${esc(t.from)} → ${esc(t.to)}  <b style="color:var(--text)">${fx(t.coeff, 5)}</b>  (${(transWeight(d, t) * 100).toFixed(1)}%)</div>`).join('') : '<p class="note">No orbital transitions listed for this state.</p>'}
        ${isNum(cur.s2) ? `<p class="note">⟨S²⟩ = ${cur.s2.toFixed(3)}</p>` : ''}</div>`;
    }
    return `<div class="sec"><table class="list" id="states-table"><tr><th>#</th><th>E (eV)</th><th>λ (nm)</th><th>f</th><th>Label</th></tr>${rows}</table>
      <p class="note">Contributions are 2c² for closed-shell Gaussian output and c² otherwise.</p></div>${det}`;
  }

  // ---------------- Surfaces tab: computed orbitals and densities, molecular surfaces, ESP colouring ----------------
  const ORB_GRIDS = { coarse: { h: 0.45, max: 1.5e6, name: 'Coarse' }, medium: { h: 0.3, max: 4e6, name: 'Medium' }, fine: { h: 0.2, max: 1e7, name: 'Fine' } };
  const KIND_NAMES = { orbital: 'orbital', density: 'density', spin: 'spin density', molsurf: 'surface', potential: 'potential', other: 'volume' };
  /** The orbital section's choices for a document (created on first use, HOMO selected). */
  function orbState(d, M) {
    if (!d.orb) {
      const f = BM.orbitals.frontier(M, 'alpha');
      d.orb = { kind: 'mo', spin: 'alpha', index: Math.max(0, f.homo), grid: 'medium' };
    }
    if (!M.sets[d.orb.spin]) d.orb.spin = 'alpha';
    return d.orb;
  }
  function orbOptions(M, spin, sel) {
    const items = BM.orbitals.list(M).filter(x => x.spin === spin), EV = HARTREE.ev;
    return items.map(x => `<option value="${x.index}" ${x.index === sel ? 'selected' : ''}>${x.index + 1}${x.tag ? ' · ' + x.tag : ''}${isNum(x.energy) ? ' · ' + (x.energy * EV).toFixed(2) + ' eV' : ''}${x.occ ? ' · occ ' + (+x.occ.toFixed(3)) : ''}${x.sym ? ' · ' + esc(x.sym) : ''}</option>`).join('');
  }
  /** Progress text: percent, and the time left once it can be told. */
  function jobText(job) {
    if (job.progress === null) return 'working…';
    const f = job.progress || 0, pc = Math.round(100 * f), el = (Date.now() - job.t0) / 1000;
    if (f < 0.05 || el < 2) return pc + '%';
    const left = el * (1 - f) / f;
    return pc + '% · ' + (left < 90 ? Math.ceil(left) + ' s' : Math.ceil(left / 60) + ' min') + ' left';
  }
  function jobBar(job, act) {
    const pc = job.progress === null ? 100 : Math.round(100 * (job.progress || 0));
    return `<div class="sf-job"><span class="sf-bar${job.progress === null ? ' busy' : ''}"><i style="width:${pc}%"></i></span><span class="note mono" data-sf-pc>${esc(jobText(job))}</span><button class="btn small" data-act="${act}">Cancel</button></div>`;
  }
  function panelOrbitals(d, M) {
    const o = orbState(d, M), det = M.detected, uhf = !!M.sets.beta, spin = BM.orbitals.hasSpin(M);
    const kinds = [['mo', 'Orbital'], ['density', 'Electron density']].concat(spin ? [['spin', 'Spin density']] : []);
    const busy = d.orbJob;
    const conv = !det ? 'Checking the orbitals…' : det.skipped ? `The orbitals were not checked (${M.nbf} basis functions)` : det.error <= 1e-3 ? `MOs orthonormal to ${det.error.toExponential(0)} (${esc(BM.orbitals.describe(det))})` : `<span class="warn">The MOs are not orthonormal in any convention BlueMolecule knows (deviation ${det.error.toExponential(1)}); the result may be wrong.</span>`;
    return `<div class="sec" id="sf-orb"><h3>Orbitals and densities</h3>
      <div class="controls">
        <label class="field"><span>Show</span><select data-orb="kind">${kinds.map(([k, n]) => `<option value="${k}" ${o.kind === k ? 'selected' : ''}>${n}</option>`).join('')}</select><span></span></label>
        ${o.kind === 'mo' ? `${uhf ? `<label class="field"><span>Spin</span><select data-orb="spin"><option value="alpha" ${o.spin === 'alpha' ? 'selected' : ''}>α</option><option value="beta" ${o.spin === 'beta' ? 'selected' : ''}>β</option></select><span></span></label>` : ''}
        <label class="field"><span>Orbital</span><select data-orb="index">${orbOptions(M, o.spin, o.index)}</select><span></span></label>` : ''}
        <label class="field" title="Grid spacing: coarse 0.45, medium 0.3, fine 0.2 Bohr"><span>Grid</span><select data-orb="grid">${Object.keys(ORB_GRIDS).map(k => `<option value="${k}" ${o.grid === k ? 'selected' : ''}>${ORB_GRIDS[k].name}</option>`).join('')}</select><span></span></label>
      </div>
      ${busy ? jobBar(busy, 'sf-orb-cancel') : `<div class="btnrow"><button class="btn small primary" data-act="sf-orb">Compute</button>${o.kind === 'mo' ? `<button class="btn small" data-act="sf-homo" title="Compute the highest occupied orbital">HOMO</button><button class="btn small" data-act="sf-lumo" title="Compute the lowest unoccupied orbital">LUMO</button>` : ''}</div>`}
      <p class="note">Computed from the basis set and orbitals in the file, as cubegen and orca_plot do.${conv ? ' ' + conv + '.' : ''}</p>
    </div>`;
  }
  function panelMolSurf(d) {
    const m = d.msurf || (d.msurf = { type: 'ses', probe: 1.4 });
    const busy = d.msJob;
    return `<div class="sec" id="sf-ms"><h3>Molecular surface</h3>
      <div class="controls">
        <label class="field"><span>Surface</span><select data-ms="type">${[['vdw', 'van der Waals'], ['sas', 'Solvent accessible'], ['ses', 'Solvent excluded']].map(([k, n]) => `<option value="${k}" ${m.type === k ? 'selected' : ''}>${n}</option>`).join('')}</select><span></span></label>
        ${m.type !== 'vdw' ? `<label class="field" title="1.4 Å is water"><span>Probe radius</span><input type="number" data-ms="probe" value="${m.probe}" min="0.1" max="5" step="0.1"><span class="note">Å</span></label>` : ''}
      </div>
      ${busy ? jobBar(busy, 'sf-ms-cancel') : `<div class="btnrow"><button class="btn small primary" data-act="sf-ms">Compute</button></div>`}
    </div>`;
  }
  /** Sources of a surface's colour: atomic charges of the current geometry, and potential (or other) volumes of the
   *  open documents on the same atoms. */
  function espSources(d, s) {
    const out = [['none', 'Single colour']];
    for (const k of chargeSchemes(d)) out.push(['q:' + k, 'ESP of the ' + k + ' charges']);
    const f = curFrame(d);
    for (const doc of state.docs) {
      const g = doc.data.frames[doc.frameIndex];
      const same = doc === d || (g && g.xyz && f.xyz && g.Z.length === f.Z.length && g.Z.every((z, i) => z === f.Z[i]) && Array.from(g.xyz).every((x, i) => Math.abs(x - f.xyz[i]) < 2e-3));
      if (!same) continue;
      doc.surfaces.forEach((x, i) => {
        if (x === s || x.vol.kind === 'molsurf' || x.vol.kind === 'orbital' || x.vol.kind === 'density') return;
        out.push(['v:' + doc.id + ':' + i, (x.vol.kind === 'potential' ? 'Potential: ' : 'Volume: ') + (x.vol.name || 'volume') + (doc === d ? '' : ' (' + doc.name + ')')]);
      });
    }
    return out;
  }
  function surfCard(d, s, i) {
    const k = s.vol.kind, mx = s.stats.absMax || 1, lo = Math.max(mx * 1e-4, 1e-6), pos = Math.log(s.iso / lo) / Math.log(mx / lo);
    const canColour = k === 'density' || k === 'molsurf';
    const srcs = canColour ? espSources(d, s) : null;
    const e = s.esp || (s.esp = { src: 'none', range: null });
    if (srcs && !srcs.some(x => x[0] === e.src)) e.src = 'none';
    const r = e.src !== 'none' ? (e.range > 0 ? e.range : e.auto) : null;
    const info = k === 'molsurf' ? `${esc(s.vol.name)} · grid ${s.vol.surface.spacing.toFixed(2)} Å` + (s.frameLabel ? ' · ' + esc(s.frameLabel) : '')
      : `grid ${s.vol.dims.join(' × ')} · min ${s.stats.min.toExponential(2)} · max ${s.stats.max.toExponential(2)}` + (s.vol.computed ? ` · spacing ${(s.vol.computed.spacing * 0.529177).toFixed(2)} Å` : '');
    return `<div class="surf" data-surf="${i}">
        <div class="surf-head"><input type="checkbox" data-s="visible" ${s.visible ? 'checked' : ''}><span class="t" title="${esc(s.vol.name)}">${esc(s.vol.name || 'Volume ' + (i + 1))}</span><span class="badge pill">${esc(KIND_NAMES[k] || k)}</span>${s.computed ? `<button class="btn small icon-only" data-act="sf-del" title="Remove">×</button>` : ''}</div>
        <div class="controls" style="margin:4px 0 0">
          ${k === 'molsurf' ? '' : `<label class="field"><span>Isovalue</span><input type="range" data-s="isoSlider" min="0" max="1" step="0.001" value="${Math.max(0, Math.min(1, pos)).toFixed(3)}"><input type="number" data-s="iso" value="${s.iso}" step="any" style="width:84px"></label>`}
          <label class="field"><span>Opacity</span><input type="range" data-s="opacity" min="0.15" max="1" step="0.05" value="${s.opacity}"><output>${Math.round(s.opacity * 100)}%</output></label>
          ${canColour ? `<label class="field"><span>Colour</span><select data-s="esp">${srcs.map(([v, n]) => `<option value="${esc(v)}" ${v === e.src ? 'selected' : ''}>${esc(n)}</option>`).join('')}</select><span></span></label>` : ''}
          ${e.src !== 'none' ? `<label class="field" title="Potential at the ends of the colour scale (atomic units; empty: from the surface)"><span>Range ±</span><input type="number" data-s="espRange" value="${e.range > 0 ? e.range : ''}" placeholder="${isNum(e.auto) ? e.auto.toPrecision(2) : 'auto'}" step="any" min="0"><span class="note">a.u.</span></label>
            <div class="esp-scale"><span>${isNum(r) ? '−' + r.toPrecision(2) : ''}</span><i></i><span>${isNum(r) ? '+' + r.toPrecision(2) : ''}</span></div>
            <p class="note">${isNum(r) ? `Red: negative potential (electron rich), blue: positive. ±${r.toPrecision(2)} a.u. = ±${(r * HARTREE.kcal).toFixed(1)} kcal/mol.` : ''}${e.src.startsWith('q:') ? ' From point charges at the atoms: an approximation of the real potential.' : ''}</p>` : `<div class="field"><span>Colours</span><span class="inline"><input type="color" data-s="colorPos" value="${s.colorPos}" title="Positive">${s.stats.min < 0 && k !== 'molsurf' ? `<input type="color" data-s="colorNeg" value="${s.colorNeg}" title="Negative"><label class="inline"><input type="checkbox" data-s="both" ${s.both ? 'checked' : ''}> ± lobes</label>` : ''}</span><span></span></div>`}
        </div>
        <p class="note mono">${info}</p>
      </div>`;
  }
  function panelSurf(d) {
    const M = BM.orbitals ? BM.orbitals.model(d.data) : null;
    // the orthonormality check (an overlap matrix: up to a few seconds for thousands of functions), after the panel
    if (M && !M.detected && !M._checking) {
      M._checking = true;
      setTimeout(() => {
        try { BM.orbitals.detect(M); } catch (e) { console.warn('[BlueMolecule] orbital check failed:', e); }
        M._checking = false;
        if (state.active === d && d.tab === 'surf' && !d.orbJob) renderPanel();
      }, 30);
    }
    const cards = d.surfaces.map((s, i) => surfCard(d, s, i)).join('');
    const tip = M ? '' : d.surfaces.some(s => s.vol.kind === 'orbital') ? '' : '<p class="note">To view orbitals, open the calculation’s fchk or Molden file, or a cube from cubegen (Gaussian) or orca_plot (ORCA).</p>';
    return (M ? panelOrbitals(d, M) : '') + panelMolSurf(d) + (cards ? `<div class="sec"><h3>Surfaces</h3>${cards}</div>` : '') + tip;
  }
  /** Index of a frame with the basis set's geometry (the orbitals belong to it), or -1. */
  function basisFrame(d, M) {
    const fr = d.data.frames, B = 0.529177210903;
    const cen = new Map();
    for (const sh of M.shells) if (!cen.has(sh.atom)) cen.set(sh.atom, sh.center);
    const fits = f => f && f.xyz && Array.from(cen).every(([a, c]) => a < f.Z.length && Math.abs(f.xyz[3 * a] - c[0] * B) < 2e-3 && Math.abs(f.xyz[3 * a + 1] - c[1] * B) < 2e-3 && Math.abs(f.xyz[3 * a + 2] - c[2] * B) < 2e-3);
    if (fits(fr[d.frameIndex])) return d.frameIndex;
    for (let k = fr.length - 1; k >= 0; k--) if (fits(fr[k])) return k;
    return -1;
  }
  function addComputed(d, vol, key) {
    const s = makeSurface(vol, 0);
    s.visible = true; s.computed = true; s.key = key;
    // the new surface alone (the others stay in the list, to be shown again with their check boxes)
    for (const x of d.surfaces) x.visible = false;
    d.surfaces.push(s);
    // keep memory in check: at most 12 computed volumes, the oldest hidden ones go first
    let comp = d.surfaces.filter(x => x.computed);
    while (comp.length > 12) { const old = comp.find(x => !x.visible) || comp[0]; d.surfaces.splice(d.surfaces.indexOf(old), 1); comp = d.surfaces.filter(x => x.computed); }
    return s;
  }
  function computeOrbital(d) {
    const M = BM.orbitals.model(d.data);
    if (!M || d.orbJob) return;
    const o = orbState(d, M);
    const what = o.kind === 'mo' ? { kind: 'mo', spin: o.spin, index: o.index } : { kind: o.kind };
    const key = (what.kind === 'mo' ? 'mo:' + what.spin + ':' + what.index : what.kind) + '@' + o.grid;
    const have = d.surfaces.find(x => x.key === key);
    const fi = basisFrame(d, M);
    if (fi >= 0 && fi !== d.frameIndex) { setFrame(d, fi); status('Showing the geometry of the orbitals (frame ' + (fi + 1) + ').'); }
    if (have) {
      for (const x of d.surfaces) x.visible = x === have;
      renderTabs(); renderPanel(); refreshSurfaces(); return;
    }
    let job;
    try {
      const det = BM.orbitals.detect(M);
      if (det.error > 1e-3) toast('The orbitals of this file are not orthonormal in any convention BlueMolecule knows; the surface may be wrong.', true);
      const g = ORB_GRIDS[o.grid] || ORB_GRIDS.medium, grid = BM.orbitals.gridFor(M, { spacing: g.h, maxPoints: g.max });
      job = BM.orbitals.compute(M, what, grid, { onProgress: f => jobProgress(d, d.orbJob, f) });
    } catch (e) { toast('Could not compute: ' + e.message, true); return; }
    const t0 = Date.now();
    d.orbJob = { cancel: job.cancel, progress: 0, t0: Date.now() };
    renderPanel();
    job.promise.then(vol => {
      d.orbJob = null;
      const s = addComputed(d, vol, key);
      if (state.active === d) { renderTabs(); renderPanel(); refreshSurfaces(); includeSurface(s); }
      const th = vol.computed && vol.computed.workers;
      status(vol.name + ': ' + vol.dims.join(' × ') + ' points in ' + ((Date.now() - t0) / 1000).toFixed(1) + ' s' + (th ? ' (' + th + (th > 1 ? ' threads)' : ' thread)') : ''));
    }).catch(e => {
      d.orbJob = null;
      if (state.active === d) renderPanel();
      if (e && e.message !== 'cancelled') toast('Could not compute the orbital: ' + e.message, true);
    });
  }
  function computeMolSurface(d) {
    if (d.msJob) return;
    const m = d.msurf || (d.msurf = { type: 'ses', probe: 1.4 }), g = geomFor(d);
    let job;
    try { job = BM.molsurf.fieldAsync(g.Z, g.xyz, { type: m.type, probe: m.probe }); } catch (e) { toast('Could not compute the surface: ' + e.message, true); return; }
    const t0 = Date.now(), frame = d.frameIndex;
    d.msJob = { cancel: job.cancel, progress: null, t0: Date.now() };
    renderPanel();
    job.promise.then(vol => {
      d.msJob = null;
      const s = addComputed(d, vol, null);
      if (d.data.frames.length > 1) s.frameLabel = 'frame ' + (frame + 1);
      if (state.active === d) { renderTabs(); renderPanel(); refreshSurfaces(); includeSurface(s); }
      status(vol.name + ' in ' + ((Date.now() - t0) / 1000).toFixed(1) + ' s');
    }).catch(e => {
      d.msJob = null;
      if (state.active === d) renderPanel();
      if (e && e.message !== 'cancelled') toast('Could not compute the surface: ' + e.message, true);
    });
  }
  /** Zoom out, if needed, so that a new surface is in view. */
  function includeSurface(s) {
    if (!s || !s.mesh) return;
    for (const k of Object.keys(s.mesh)) if (s.mesh[k] && s.mesh[k].positions) R.include(s.mesh[k].positions);
  }
  function jobProgress(d, job, f) {
    if (!job) return;
    job.progress = f;
    if (state.active !== d || d.tab !== 'surf') return;
    const b = $('#sf-orb .sf-bar i'), t = $('#sf-orb [data-sf-pc]');
    if (b) b.style.width = Math.round(100 * f) + '%';
    if (t) t.textContent = jobText(job);
  }
  function surfAction(act, a) {
    const d = state.active;
    if (!d) return;
    const M = BM.orbitals.model(d.data);
    if (act === 'sf-orb') computeOrbital(d);
    else if (act === 'sf-orb-cancel') { if (d.orbJob) d.orbJob.cancel(); }
    else if (act === 'sf-ms') computeMolSurface(d);
    else if (act === 'sf-ms-cancel') { if (d.msJob) d.msJob.cancel(); }
    else if ((act === 'sf-homo' || act === 'sf-lumo') && M) {
      const o = orbState(d, M), f = BM.orbitals.frontier(M, o.spin), k = act === 'sf-homo' ? f.homo : f.lumo;
      if (k >= 0) { o.index = k; renderPanel(); computeOrbital(d); }
    } else if (act === 'sf-del') {
      const card = a.closest('[data-surf]'), i = card ? +card.dataset.surf : -1;
      if (i >= 0 && d.surfaces[i] && d.surfaces[i].computed) { d.surfaces.splice(i, 1); renderTabs(); renderPanel(); refreshSurfaces(); }
    }
  }
  /** Show an orbital picked in the MO level plot (id: column * 100000 + 1-based index). */
  function showPlotOrbital(d, id) {
    const M = BM.orbitals.model(d.data);
    if (!M || id === null || id === undefined) return;
    const cols = BM.analysis.moLevels(d.data), c = cols[Math.floor(id / 100000)];
    if (!c) return;
    const o = orbState(d, M);
    o.kind = 'mo'; o.spin = c.spin === 'beta' && M.sets.beta ? 'beta' : 'alpha'; o.index = id % 100000 - 1;
    if (d.tab !== 'surf') { d.tab = 'surf'; renderTabs(); }
    renderPanel();
    computeOrbital(d);
  }
  /** Per-vertex colours of a surface from its potential source, or null. */
  function espColors(d, s, mesh) {
    const e = s.esp;
    if (!e || e.src === 'none' || !mesh || !mesh.positions.length) return null;
    const key = e.src + '|' + e.range + '|' + s.iso + '|' + mesh.positions.length;
    if (s._espKey === key && s._espColors) return s._espColors;
    let vals = null;
    if (e.src.startsWith('q:')) {
      const q = d.data.charges && d.data.charges[e.src.slice(2)], g = geomFor(d);
      if (q && q.length === g.Z.length) vals = BM.molsurf.espCharges(g.Z, g.xyz, q, mesh.positions);
    } else if (e.src.startsWith('v:')) {
      const [, id, i] = e.src.split(':'), doc = state.docs.find(x => x.id === +id), src = doc && doc.surfaces[+i];
      if (src) vals = BM.molsurf.sample(src.vol, mesh.positions);
    }
    if (!vals) return null;
    // the scale: the 98th percentile of |V| on the surface, unless set
    const abs = Array.from(vals).filter(Number.isFinite).map(Math.abs).sort((a, b) => a - b);
    e.auto = abs.length ? Number((abs[Math.min(abs.length - 1, Math.floor(abs.length * 0.98))] || 1e-3).toPrecision(2)) : null;
    const r = e.range > 0 ? e.range : e.auto || 0.05;
    s._espColors = BM.molsurf.colors(vals, -r, r);
    s._espKey = key;
    return s._espColors;
  }

  // ---------------- editing ----------------
  /** Explicit bonds (with orders) and formal charges for an edited frame, from its current topology. */
  function ensureTopology(d, fi) {
    const f = d.data.frames[fi], n = f.Z.length;
    if (f.bonds && f.formalCharges && f.formalCharges.length === n) return;
    if (!f.bonds) {
      const t = topoFor(d, fi);
      let bonds = null, fc = null;
      if (!t.perceived) bonds = kit.pairsToList(t.pairs, t.orders);
      else if (n <= PERCEIVE_MAX && t.pairs.length) {
        try { const p = kit.perceive(f.Z, f.xyz, t.pairs, null); bonds = p.bonds; fc = p.charges; } catch (e) { console.warn('[BlueMolecule] bond-order perception failed:', e); }
      }
      f.bonds = bonds || kit.pairsToList(t.pairs, null);
      if (!f.formalCharges) f.formalCharges = fc;
    }
    if (!f.formalCharges || f.formalCharges.length !== n) f.formalCharges = new Array(n).fill(0);
    invalidateTopo(d, fi);
  }
  function ensureEditFrame(d) {
    const f = curFrame(d);
    if (f.edited) { ensureTopology(d, d.frameIndex); return f; }
    const t = topoFor(d, d.frameIndex);
    let bonds = null, fc = null;
    if (!t.perceived) bonds = kit.pairsToList(t.pairs, t.orders);
    else if (f.Z.length <= PERCEIVE_MAX && t.pairs.length) {
      if (t.orders) { bonds = kit.pairsToList(t.pairs, t.orders); fc = t.charges; }
      else { try { const p = kit.perceive(f.Z, f.xyz, t.pairs, null); bonds = p.bonds; fc = p.charges; } catch (e) { console.warn('[BlueMolecule] bond-order perception failed:', e); } }
    }
    const nf = {
      Z: f.Z.slice(), xyz: Float64Array.from(d._animXyz && state.anim.kind === 'vib' ? d._animXyz : f.xyz), energy: null, energies: null, label: 'Edited', converged: false, opt: null,
      atomNames: f.atomNames ? f.atomNames.slice() : null, resNames: f.resNames ? f.resNames.slice() : null, resIds: f.resIds ? f.resIds.slice() : null, edited: true,
      cell: (c => (c ? c.map(v => v.slice()) : null))(cellOf(d, f)),
      bonds: bonds || kit.pairsToList(t.pairs, null), formalCharges: fc ? fc.slice() : new Array(f.Z.length).fill(0)
    };
    d.data.frames.push(nf);
    d.frameIndex = d.data.frames.length - 1;
    return nf;
  }
  const copyArr = a => (a ? Array.from(a) : null);
  function snapshotFrame(d) {
    const f = curFrame(d);
    return { fi: d.frameIndex, Z: f.Z.slice(), xyz: Float64Array.from(f.xyz), bonds: f.bonds ? f.bonds.map(b => b.slice()) : null, formalCharges: copyArr(f.formalCharges),
      atomNames: copyArr(f.atomNames), resNames: copyArr(f.resNames), resIds: copyArr(f.resIds), selection: d.selection.slice(),
      cell: f.cell === undefined ? undefined : f.cell ? f.cell.map(v => v.slice()) : null };
  }
  function restoreFrame(d, s) {
    const f = d.data.frames[s.fi];
    f.Z = s.Z; f.xyz = s.xyz; f.bonds = s.bonds; f.formalCharges = s.formalCharges; f.atomNames = s.atomNames; f.resNames = s.resNames; f.resIds = s.resIds;
    if (s.cell !== undefined) f.cell = s.cell;
    d.frameIndex = s.fi;
    const n = f.Z.length;
    d.selection = (s.selection || d.selection).filter(i => i < n);
    d.frozen = d.frozen.filter(i => i < n);
  }
  /** Engine-ready copy of a frame: {Z, xyz, bonds, formalCharges} (engines never mutate their input). */
  function frameMol(f) {
    const n = f.Z.length;
    return { Z: Array.from(f.Z), xyz: Float64Array.from(f.xyz), bonds: (f.bonds || []).map(b => [b[0], b[1], b[2] === undefined ? 1 : +b[2]]),
      formalCharges: f.formalCharges && f.formalCharges.length === n ? Array.from(f.formalCharges, q => +q || 0) : new Array(n).fill(0) };
  }
  /** Store an engine result in frame f; map (old → new index, −1 deleted; null = unchanged) remaps names, selection and frozen atoms. */
  function commitMol(f, d, mol, map) {
    const n0 = f.Z.length, n1 = mol.Z.length, idx = i => (map ? map[i] : i);
    const remap = a => { if (!a) return null; const out = new Array(n1).fill(null); for (let i = 0; i < n0 && i < a.length; i++) { const j = idx(i); if (j >= 0 && j < n1) out[j] = a[i]; } return out; };
    const remapList = l => [...new Set(l.map(idx).filter(j => j >= 0 && j < n1))];
    f.atomNames = remap(f.atomNames); f.resNames = remap(f.resNames); f.resIds = remap(f.resIds);
    f.Z = Array.from(mol.Z); f.xyz = Float64Array.from(mol.xyz);
    f.bonds = (mol.bonds || []).map(b => [b[0], b[1], b[2] === undefined ? 1 : +b[2]]);
    f.formalCharges = mol.formalCharges && mol.formalCharges.length === n1 ? Array.from(mol.formalCharges, q => +q || 0) : new Array(n1).fill(0);
    d.selection = remapList(d.selection); d.frozen = remapList(d.frozen);
  }
  /** One undoable edit of the current geometry (an "Edited" copy is made first for file geometries).
   *  fn(f, d) may return false to cancel without an undo step. opts: {fit, auto} (see ui.afterEdit). */
  function doEdit(fn, msg, opts) {
    const d = state.active; if (!d) return false;
    stopAnim(true);
    if (ui) ui.relaxStop(true);
    const before = d.data.frames.length, fi0 = d.frameIndex;
    const f = ensureEditFrame(d);
    const snap = snapshotFrame(d), frozen0 = d.frozen.slice();
    const cancel = () => {
      restoreFrame(d, snap); d.frozen = frozen0;
      if (d.data.frames.length > before) { d.data.frames.length = before; d.frameIndex = fi0; invalidateTopo(d, before); }
    };
    let r;
    try { r = fn(f, d); } catch (e) { console.warn('[BlueMolecule] edit failed:', e); cancel(); toast(e.message || String(e), true); return false; }
    if (r === false) { cancel(); return false; }
    d.undo.push(snap); d.redo = [];
    if (d.undo.length > 200) d.undo.shift();
    invalidateTopo(d, d.frameIndex);
    refreshScene(true); updateRadius(); renderTabs(); renderPanel(); renderPlaybar(); updateStatus();
    ui.afterEdit(d, opts || {});
    const m = typeof msg === 'function' ? msg() : msg;
    if (m) status(m);
    return true;
  }
  /** Start a drag-style edit (manipulate, bond rotation, clean-up run): one undo step, positions change in place. */
  function beginContinuousEdit(d) {
    stopAnim(true);
    if (ui) ui.relaxStop(true);
    const before = d.data.frames.length;
    const f = ensureEditFrame(d);
    d.undo.push(snapshotFrame(d)); d.redo = [];
    if (d.undo.length > 200) d.undo.shift();
    if (d.data.frames.length !== before) { refreshScene(true); renderTabs(); renderPlaybar(); }
    return f;
  }
  function endContinuousEdit(d, msg) {
    if (d !== state.active) return;
    d._sym = null; d._ff = null; d._smiles = null;
    afterPositions(d);
    if (msg) status(msg);
  }
  /** Positions changed outside doEdit (drags, clean-up): refresh everything that depends on them. */
  function afterPositions(d) {
    if (!d || d !== state.active) return;
    R.updatePositions(curFrame(d).xyz);
    refreshOverlays(d);
    refreshExtras(); updateRadius(); updateMeasure(true); renderPanel(); updateStatus();
  }
  /** Keep the renderer's scene radius (clipping and depth cueing) in step with a growing structure. */
  function updateRadius() { if (R && typeof R._updateRadius === 'function') { R._updateRadius(); R.requestRender(); } }
  function undo(redo) {
    const d = state.active; if (!d) return;
    if (ui) { ui.relaxStop(true); ui.cancelDrag(); }
    const from = redo ? d.redo : d.undo, to = redo ? d.undo : d.redo;
    if (!from.length) { status(redo ? 'Nothing to redo' : 'Nothing to undo'); return; }
    stopAnim(true);
    const s = from.pop();
    d.frameIndex = s.fi;
    to.push(snapshotFrame(d));
    restoreFrame(d, s);
    invalidateTopo(d, s.fi);
    d._sym = null; d._ff = null; d._smiles = null; d._coordDirty = false;
    refreshScene(true); updateRadius(); renderTabs(); renderPanel(); renderPlaybar(); updateStatus();
    status(redo ? 'Redone' : 'Undone');
  }
  function editSet() {
    const d = state.active, val = parseFloat($('#edit-val').value);
    if (!isNum(val)) return toast('Enter a number first.', true);
    const s = d.selection.slice();
    doEdit(f => {
      const b = bondsFor(d, d.frameIndex);
      let r;
      if (s.length === 2) r = geom.setDistance(f.xyz, b, s[0], s[1], val);
      else if (s.length === 3) r = geom.setAngle(f.xyz, b, s[0], s[1], s[2], val);
      else if (s.length === 4) r = geom.setDihedral(f.xyz, b, s[0], s[1], s[2], s[3], val);
      else throw new Error('Select 2, 3 or 4 atoms.');
      f.xyz = r.xyz;
      if (r.ring) toast('Atoms are in a ring, so only the last selected atom was moved.');
    }, 'Geometry edited');
  }
  function editElement() {
    const d = state.active, inp = $('#edit-el'), z = el.toZ(inp ? inp.value : '');
    if (!d || !d.selection.length) return;
    if (z <= 0) return toast('Unknown element symbol.', true);
    const s = d.selection.slice();
    doEdit(f => { for (const i of s) f.Z[i] = z; }, `Changed ${s.length} atom${s.length > 1 ? 's' : ''} to ${el.symbol(z)}`);
  }
  function editDelete() {
    const d = state.active; if (!d) return;
    const list = [...new Set(d.selection)];
    if (!list.length) return;
    if (list.length === curFrame(d).Z.length && !d.built) return toast('Cannot delete every atom of a file geometry (start a New molecule to build from scratch).', true);
    doEdit(f => {
      const r = BM.chem.deleteAtoms(frameMol(f), list);
      commitMol(f, d, r.mol, r.map);
      d.selection = [];
    }, `Deleted ${list.length} atom${list.length > 1 ? 's' : ''}`);
    R.setSelection(d.selection); updateMeasure();
  }
  function editCenter() {
    doEdit(f => {
      const n = f.Z.length; let cx = 0, cy = 0, cz = 0;
      if (!n) return false;
      for (let i = 0; i < n; i++) { cx += f.xyz[3 * i]; cy += f.xyz[3 * i + 1]; cz += f.xyz[3 * i + 2]; }
      cx /= n; cy /= n; cz /= n;
      for (let i = 0; i < n; i++) { f.xyz[3 * i] -= cx; f.xyz[3 * i + 1] -= cy; f.xyz[3 * i + 2] -= cz; }
    }, 'Centred at origin');
    R.fit();
  }
  function displace(sign) {
    const d = state.active, v = d && d.data.vibrations;
    if (!v || d.vibIndex < 0) return;
    const amt = parseFloat($('#disp-amt') ? $('#disp-amt').value : state.vib.disp) || state.vib.disp || 0.15;
    const mode = d.vibIndex, u = vibVec(d, mode, amt * sign);
    if (!u) return toast(NOT_LOADED, true);
    const eq = d.data.frames[v.frameIndex];
    stopAnim(true);
    const nf = { Z: eq.Z.slice(), xyz: Float64Array.from(eq.xyz, (x, i) => x + u[i]), energy: null, energies: null, label: `ν${mode + 1} ${sign > 0 ? '+' : '−'}${amt}`, converged: false, opt: null, atomNames: null, resNames: null, resIds: null, edited: true };
    d.data.frames.push(nf);
    d.frameIndex = d.data.frames.length - 1;
    refreshScene(true); renderTabs(); renderPanel(); renderPlaybar(); updateStatus();
    toast(`New geometry displaced ${sign > 0 ? '+' : '−'}${amt} Å along mode ${mode + 1}. Export it from the Export menu.`);
  }

  // ---------------- plot ----------------
  const framePlot = t => t === 'energy' || t === 'scan' || t === 'irc' || (!!t && t.startsWith('opt:'));
  function convFor(d) {
    const n = d.data.frames.length;
    if (!d._conv || d._conv.job !== d.job || d._conv.n !== n) d._conv = { job: d.job, n, c: BM.analysis.convergenceSeries(d.data, { job: d.job }) };
    return d._conv.c;
  }
  function profFor(d) {
    const n = d.data.frames.length;
    if (!d._prof || d._prof.job !== d.job || d._prof.n !== n) d._prof = { job: d.job, n, p: BM.analysis.profiles(d.data, { job: d.job }) };
    return d._prof.p;
  }
  /** Atoms the NMR shieldings refer to: the current geometry if it matches, else the last one. */
  function nmrZ(d) {
    const n = d.data.nmr && d.data.nmr.iso ? d.data.nmr.iso.length : 0, f = curFrame(d);
    return f.Z.length === n ? f.Z : d.data.frames[d.data.frames.length - 1].Z;
  }
  function chargeSchemes(d) {
    const c = d.data.charges, n = curFrame(d).Z.length;
    return c ? Object.keys(c).filter(k => Array.isArray(c[k]) && c[k].length === n) : [];
  }
  /** The plot picker: [[group, [[id, label], ...]], ...] for what the document holds. */
  function plotGroups(d) {
    const g = [], data = d.data, v = data.vibrations, c = convFor(d), prof = profFor(d);
    const opt = c.isOpt && c.steps.length >= 2;
    if (opt) {
      const items = c.series.map(x => ['opt:' + x.key, x.label + (x.unit ? ' (' + x.unit + ')' : '')]);
      if (c.components.length) items.splice(Math.min(2, items.length), 0, ['opt:parts', 'Energy components (ΔE vs final)']);
      if (c.ratios.length) items.push(['opt:all', 'All criteria / threshold']);
      g.push(['Optimisation', items]);
    }
    const geo = prof.map(p => [p.kind, p.kind === 'scan' ? 'Scan: energy vs ' + p.xLabel.replace(/ \(.*\)$/, '') : 'IRC: energy vs ' + p.xLabel.toLowerCase()]);
    if (!opt && visibleFrames(d).filter(i => isNum(data.frames[i].energy)).length > 1) geo.push(['energy', 'Energy profile']);
    if (geo.length) g.push([prof.length ? 'Scan / IRC' : 'Geometries', geo]);
    const sp = [], S = data.excitedStates;
    if (v && v.ir && v.freqs.some(f => f > 0)) sp.push(['ir', 'IR']);
    if (v && v.raman && v.freqs.some(f => f > 0)) sp.push(['raman', 'Raman']);
    if (S && S.some(s => isNum(s.f))) sp.push(['uv', 'UV-Vis']);
    if (S && S.some(s => isNum(s.R) && s.R !== 0)) sp.push(['ecd', 'ECD']);
    if (data.nmr && data.nmr.iso && BM.analysis.nmrNuclei(data, nmrZ(d)).length) sp.push(['nmr', 'NMR']);
    if (sp.length) g.push(['Spectra', sp]);
    if (BM.analysis.moLevels(data).length) g.push(['Orbitals', [['mo', 'MO energy levels']]]);
    if (chargeSchemes(d).length && curFrame(d).Z.length <= 5000) g.push(['Charges', [['charges', 'Atomic charges']]]);
    return g;
  }
  function plotTypes(d) { return d && !d.pending ? plotGroups(d).reduce((a, x) => a.concat(x[1]), []) : []; }
  function defaultPlot(types) {
    const has = id => types.some(t => t[0] === id);
    for (const id of ['ir', 'raman', 'uv', 'nmr', 'scan', 'irc', 'opt:dE', 'opt:energy', 'energy', 'ecd', 'mo', 'charges']) if (has(id)) return id;
    return types[0][0];
  }
  function renderPlot() {
    const d = state.active, center = $('#center');
    const groups = d && !d.pending ? plotGroups(d) : [];
    const types = groups.reduce((a, x) => a.concat(x[1]), []);
    center.classList.toggle('no-plot', !types.length);
    if (!types.length) { $('#plot-pick').innerHTML = ''; $('#plot-controls').innerHTML = ''; plot.set(null); return; }
    if (!types.some(t => t[0] === d.plotType)) d.plotType = defaultPlot(types);
    $('#plot-pick').innerHTML = groups.map(([g, items]) => `<optgroup label="${g}">${items.map(([id, n]) => `<option value="${id}" ${id === d.plotType ? 'selected' : ''}>${esc(n)}</option>`).join('')}</optgroup>`).join('');
    $('#plot-count').textContent = types.length > 1 ? types.length : '';
    $('#plot-pick-wrap').title = types.length + ' plot' + (types.length > 1 ? 's' : '') + ' for this file: ' + groups.map(x => x[0]).join(', ');
    $('#plot-controls').innerHTML = plotControls(d, groups);
    renderPlotDraw();
  }
  const numIn = (k, v, step, w) => `<input type="number" data-pc="${k}" value="${v === null || v === undefined ? '' : v}" step="${step || 'any'}" ${w ? `style="width:${w}px"` : ''}>`;
  const selIn = (k, v, opts) => `<select data-pc="${k}">${opts.map(([a, b]) => `<option value="${a}" ${String(a) === String(v) ? 'selected' : ''}>${b}</option>`).join('')}</select>`;
  const chk = (k, v, label, title) => `<label ${title ? `title="${title}"` : ''}><input type="checkbox" data-pc="${k}" ${v ? 'checked' : ''}>${label}</label>`;
  function specKey(t) { return t === 'energy' || t === 'scan' || t === 'irc' ? 'energy' : t.startsWith('opt:') ? 'opt' : t; }
  function plotControls(d, groups) {
    const t = d.plotType, s = d.spec[specKey(t)];
    const optGroup = groups.find(g => g[0] === 'Optimisation');
    const go = optGroup && !t.startsWith('opt:') ? `<button class="mini go" data-plot-go="${optGroup[1].some(x => x[0] === 'opt:all') ? 'opt:all' : 'opt:dE'}" title="Optimisation plots: energy, criteria vs thresholds, RMSD">Convergence ›</button>` : '';
    if (t.startsWith('opt:')) {
      const c = convFor(d), ser = c.series.find(x => 'opt:' + x.key === t);
      const st = c.converged >= 0 ? `<span class="pill ok">converged at ${c.source === 'optSteps' ? 'cycle' : 'step'} ${c.steps[c.converged].n}</span>` : '<span class="pill warn">not converged</span>';
      return (ser && ser.log ? chk('log', s.log, 'log y') : '') + (ser && ser.threshold ? `<span class="dim">threshold <span class="mono">${ser.threshold.toExponential(2)}</span></span>` : '') + st +
        (c.source === 'optSteps' ? '<span class="dim" title="This output reports optimisation cycles without their geometries">cycles from the output</span>' : '');
    }
    if (t === 'energy' || t === 'scan' || t === 'irc') return go + `<label>Unit ${selIn('unit', s.unit, Object.entries(UNIT_LABEL))}</label><label>Relative to ${selIn('ref', s.ref, [['min', 'lowest'], ['first', 'first'], ['last', 'last']])}</label>`;
    if (t === 'uv' || t === 'ecd') return go + `<label>Axis ${selIn('unit', s.unit, [['nm', 'λ (nm)'], ['eV', 'E (eV)']])}</label><label title="Gaussian half-width at 1/e">σ (eV) ${numIn('sigma', s.sigma, 0.05)}</label><label>Range ${numIn('xmin', s.xmin, 10)}–${numIn('xmax', s.xmax, 10)}</label>${chk('sticks', s.sticks, t === 'uv' ? 'f sticks' : 'R sticks')}${chk('peaks', s.peaks, 'labels')}`;
    if (t === 'nmr') {
      const nuc = BM.analysis.nmrNuclei(d.data, nmrZ(d));
      if (!nuc.some(n => n.z === s.z)) s.z = nuc[0].z;
      return go + `<label>Nucleus ${selIn('z', s.z, nuc.map(n => [n.z, n.label + ' (' + n.count + ')']))}</label><label title="Reference shielding σref (ppm): δ = σref − σ. Leave 0 to plot shieldings.">σref ${numIn('ref', s.ref, 0.1, 70)}</label>${chk('merge', s.merge, 'merge', 'Merge near-equal shieldings (within 0.05 ppm) into one stick')}${chk('peaks', s.peaks, 'labels')}`;
    }
    if (t === 'mo') return go + `<label title="Energy window below the HOMO and above the LUMO">Window ± ${numIn('window', s.window, 1, 52)} eV</label><span class="pc-read" id="mo-read">${d.moSel !== null ? esc(moLabel(d, d.moSel)) : '<span class="dim">click a level</span>'}</span>` +
      (BM.orbitals && BM.orbitals.model(d.data) ? `<button class="btn small" data-show-mo ${d.moSel === null ? 'disabled' : ''} title="Compute the selected orbital and show it (Surfaces tab)">Show orbital</button>` : '');
    if (t === 'charges') return go + `<label>Scheme ${selIn('scheme', d.chargeScheme, chargeSchemes(d).map(x => [x, x]))}</label><span class="dim">click a bar to select the atom</span>`;
    let h = go + `<label>${selIn('shape', s.shape, [['lorentzian', 'Lorentzian'], ['gaussian', 'Gaussian']])}</label><label>FWHM ${numIn('fwhm', s.fwhm, 1, 52)}</label><label title="Frequency scaling factor">Scale ${numIn('scale', s.scale, 0.001, 62)}</label>`;
    if (t === 'ir') h += `${d.data.vibrations.irUnit ? '' : `<label>${selIn('yMode', s.yMode, [['relative', 'Relative'], ['epsilon', 'ε']])}</label>`}${chk('invert', s.invert, 'transmittance')}`;
    else h += `<label>${selIn('yMode', s.yMode, [['activity', 'Activity'], ['intensity', 'Intensity']])}</label>${s.yMode === 'intensity' ? `<label>λ₀ ${numIn('laser', s.laser, 1, 54)} nm</label><label>T ${numIn('temperature', s.temperature, 1, 58)} K</label>` : ''}`;
    h += `<label>Range ${numIn('xmin', s.xmin, 50, 58)}–${numIn('xmax', s.xmax, 50, 58)}</label>${chk('sticks', s.sticks, 'sticks')}${chk('peaks', s.peaks, 'labels')}`;
    return h;
  }
  function uvRange(s, S) {
    let xmin, xmax;
    if (s.unit === 'eV') { const e = S.map(x => x.energyEV); xmin = isNum(s.xmin) ? s.xmin : Math.max(0.5, Math.floor(Math.min(...e) - 1)); xmax = isNum(s.xmax) ? s.xmax : Math.ceil(Math.max(...e) + 1); }
    else { const l = S.map(x => 1239.84198 / x.energyEV); xmin = isNum(s.xmin) ? s.xmin : Math.max(100, Math.floor((Math.min(...l) - 40) / 10) * 10); xmax = isNum(s.xmax) ? s.xmax : Math.min(1500, Math.ceil((Math.max(...l) + 80) / 10) * 10); }
    return [xmin, xmax];
  }
  function spectrumData(d) {
    const t = d.plotType, s = d.spec[t], v = d.data.vibrations;
    if (t === 'ir' || t === 'raman') {
      const fr = v.freqs.filter(f => f > 0).map(f => f * s.scale);
      const xmin = isNum(s.xmin) ? s.xmin : Math.max(0, Math.min(t === 'ir' ? 400 : 0, Math.floor((Math.min(...fr) - 150) / 100) * 100));
      const xmax = isNum(s.xmax) ? s.xmax : Math.max(t === 'ir' ? 4000 : 3600, Math.ceil((Math.max(...fr) + 150) / 100) * 100);
      const o = Object.assign({}, s, { xmin, xmax, n: 2400 });
      const r = t === 'ir' ? BM.spectrum.ir(v, o) : BM.spectrum.raman(v, o);
      return Object.assign(r, { xmin, xmax, xLabel: t === 'ir' ? 'Wavenumber (cm⁻¹)' : 'Raman shift (cm⁻¹)' });
    }
    if (t === 'uv' || t === 'ecd') {
      const S = d.data.excitedStates.filter(x => isNum(x.energyEV) && (t === 'uv' || isNum(x.R)));
      const [xmin, xmax] = uvRange(s, S);
      const r = t === 'uv' ? BM.spectrum.uv(S, { xmin, xmax, unit: s.unit, sigmaEV: s.sigma, n: 1400 }) : BM.spectrum.ecd(S, { xmin, xmax, unit: s.unit, sigmaEV: s.sigma, n: 1400 });
      return Object.assign(r, { xmin, xmax, xLabel: s.unit === 'eV' ? 'Excitation energy (eV)' : 'Wavelength (nm)' });
    }
    return null;
  }
  const csvNum = v => (isNum(v) ? String(v) : '');
  /** CSV rows follow the zoom of the plot (issue #9): inZ(z) tests x (or y, for the orbital levels) against what is shown. */
  function inZ(z, axis) {
    if (!z) return () => true;
    const a = axis === 'y' ? z.y0 : z.x0, b = axis === 'y' ? z.y1 : z.x1, tol = (b - a) * 1e-9;
    return v => isNum(v) && v >= a - tol && v <= b + tol;
  }
  /** Adds a line under the title of a zoomed CSV saying which rows it holds. */
  function zoomCSV(text, z, what, axis) {
    if (!z) return text;
    const a = axis === 'y' ? z.y0 : z.x0, b = axis === 'y' ? z.y1 : z.x1, i = text.indexOf('\n') + 1;
    return text.slice(0, i) + `# zoomed in: only ${what} from ${+a.toPrecision(7)} to ${+b.toPrecision(7)}\n` + text.slice(i);
  }
  const fmtCrit = (key, v) => (!isNum(v) ? '—' : key === 'energy' ? v.toFixed(6) + ' Eh' : key === 'dE' ? v.toFixed(3) + ' kcal/mol' : key === 'rmsd' || key === 'maxStep' ? v.toFixed(4) + ' Å' : key === 'negEig' ? String(v) : v.toExponential(3));
  /** {spec, csv()} for the document's current plot. */
  function buildPlot(d) {
    const t = d.plotType;
    if (t === 'energy') {
      const ref = energyRef(d), vis = visibleFrames(d), fr = d.data.frames, u = UNIT_LABEL[d.spec.energy.unit];
      const pts = vis.filter(i => isNum(fr[i].energy)).map(i => ({ x: i + 1, y: energyIn(d, fr[i].energy, ref), id: i, mark: fr[i].converged, label: `Geometry ${i + 1} · E = ${fr[i].energy.toFixed(6)} Eh · ΔE = ${energyIn(d, fr[i].energy, ref).toFixed(2)} ${u}` }));
      const xs = pts.map(p => p.x);
      return {
        spec: { xRange: [Math.min(...xs) - 0.5, Math.max(...xs) + 0.5], integerX: true, xLabel: 'Geometry', yLabel: `ΔE (${u})`, points: pts, selected: d.frameIndex, xFormat: x => 'Geometry ' + Math.round(x) },
        csv: z => { const ok = inZ(z); return zoomCSV(`# BlueMolecule energy profile: ${d.name}\ngeometry,energy_hartree,delta_${d.spec.energy.unit}\n` + vis.filter(i => isNum(fr[i].energy) && ok(i + 1)).map(i => `${i + 1},${fr[i].energy.toFixed(10)},${energyIn(d, fr[i].energy, ref).toFixed(6)}`).join('\n') + '\n', z, 'geometries'); }
      };
    }
    if (t.startsWith('opt:')) return optPlot(d, t.slice(4));
    if (t === 'scan' || t === 'irc') return profilePlot(d, t);
    if (t === 'nmr') return nmrPlot(d);
    if (t === 'mo') return moPlot(d);
    if (t === 'charges') return chargePlot(d);
    // spectra
    const s = d.spec[t], r = spectrumData(d), vib = d.data.vibrations, uvLike = t === 'uv' || t === 'ecd';
    const sticks = r.sticks.map(p => Object.assign(p, {
      label: uvLike ? `State ${p.id} · ${(s.unit === 'eV' ? p.x.toFixed(3) + ' eV' : p.x.toFixed(1) + ' nm')} · ${t === 'uv' ? 'f = ' + p.y.toFixed(4) : 'R = ' + p.y.toFixed(3) + ' (10⁻⁴⁰ cgs)'}`
        : `Mode ${p.id + 1} · ${p.x.toFixed(1)} cm⁻¹ · ${t === 'ir' ? fx(vib.ir[p.id], 2) + ' ' + (vib.irUnit || 'km/mol') : fx(vib.raman[p.id], 2) + ' Å⁴/amu'}`
    }));
    let yLabel = r.yLabel, invertY = false;
    if (t === 'ir' && s.invert) { yLabel = yLabel.replace('Relative absorbance', 'Transmittance (rel.)'); invertY = true; }
    const selected = uvLike ? (d.stateIndex >= 0 ? d.data.excitedStates[d.stateIndex].index : null) : (d.vibIndex >= 0 ? d.vibIndex : null);
    d._lastSpectrum = r;
    return {
      spec: {
        xRange: [r.xmin, r.xmax], reverseX: t === 'ir', xLabel: r.xLabel, yLabel, series: [{ x: r.curve.x, y: r.curve.y, fill: 'rgba(79,140,255,0.14)' }], sticks, showSticks: s.sticks, peakLabels: s.peaks, invertY,
        selected, peakFormat: uvLike ? (x => (s.unit === 'eV' ? x.toFixed(2) : Math.round(x).toString())) : null,
        xFormat: uvLike ? (x => (s.unit === 'eV' ? x.toFixed(3) + ' eV' : x.toFixed(1) + ' nm')) : (x => x.toFixed(1) + ' cm⁻¹')
      },
      csv: z => {
        const ok = inZ(z);
        const head = t === 'uv' ? `# BlueMolecule UV-Vis spectrum: ${d.name}; Gaussian sigma = ${s.sigma} eV\n${s.unit === 'eV' ? 'energy_eV' : 'wavelength_nm'},epsilon_L_per_mol_cm`
          : t === 'ecd' ? `# BlueMolecule ECD spectrum: ${d.name}; Gaussian sigma = ${s.sigma} eV\n${s.unit === 'eV' ? 'energy_eV' : 'wavelength_nm'},delta_epsilon_L_per_mol_cm`
            : `# BlueMolecule ${t === 'ir' ? 'IR' : 'Raman'} spectrum: ${d.name}; ${s.shape} FWHM = ${s.fwhm} cm-1; scale = ${s.scale}\nwavenumber_cm-1,${t === 'ir' ? (/^ε/.test(r.yLabel) ? 'epsilon_L_per_mol_cm' : 'relative_absorbance') : 'relative_' + s.yMode}`;
        const rows = [];
        for (let i = 0; i < r.curve.x.length; i++) if (ok(r.curve.x[i])) rows.push(`${r.curve.x[i].toFixed(4)},${r.curve.y[i].toExponential(6)}`);
        return zoomCSV(head + '\n' + rows.join('\n') +
          '\n\n# sticks\n' + (t === 'uv' ? 'state,x,oscillator_strength\n' : t === 'ecd' ? 'state,x,rotatory_strength_1e-40_cgs\n' : `mode,wavenumber_cm-1,${t === 'ir' ? (vib.irUnit ? 'intensity_D2' : 'intensity_km_per_mol') : 'raman_' + s.yMode}\n`) +
          r.sticks.filter(p => ok(p.x)).map(p => `${uvLike ? p.id : p.id + 1},${p.x.toFixed(4)},${p.y.toExponential(6)}`).join('\n') + '\n', z, uvLike ? (s.unit === 'eV' ? 'energies (eV)' : 'wavelengths (nm)') : 'wavenumbers (cm-1)');
      }
    };
  }
  function optPlot(d, key) {
    const c = convFor(d), steps = c.steps, xs = steps.map(s => s.n);
    const unitWord = c.source === 'optSteps' ? 'cycle' : 'step';
    const base = { xRange: [Math.min(...xs) - 0.5, Math.max(...xs) + 0.5], integerX: true, xLabel: 'Optimisation ' + unitWord, selected: d.frameIndex, xFormat: x => unitWord[0].toUpperCase() + unitWord.slice(1) + ' ' + Math.round(x) };
    const head = `# BlueMolecule optimisation convergence: ${d.name}${d.job >= 0 ? ' (job ' + (d.job + 1) + ')' : ''}\n`;
    if (key === 'all') {
      const lines = c.ratios.map(r => ({ name: r.label, pts: steps.map((s, k) => (isNum(r.values[k]) && r.values[k] > 0 ? { x: s.n, y: r.values[k], id: s.frame, mark: k === c.converged, label: `${s.label} · ${r.label} = ${r.values[k].toPrecision(3)} × threshold` } : null)).filter(Boolean) }));
      return {
        spec: Object.assign(base, { yLabel: 'Criterion / threshold', logY: true, lines, hlines: [{ y: 1, label: 'converged below 1' }] }),
        csv: z => { const ok = inZ(z); return zoomCSV(head + [unitWord, 'frame'].concat(c.ratios.map(r => r.key + '_over_threshold')).join(',') + '\n' + steps.map((s, k) => (ok(s.n) ? [s.n, s.frame === null ? '' : s.frame + 1].concat(c.ratios.map(r => csvNum(r.values[k]))).join(',') : null)).filter(x => x !== null).join('\n') + '\n', z, unitWord + 's'); }
      };
    }
    if (key === 'parts') {   // the total and the QM/MM or ONIOM components, each relative to its own last value
      const tot = c.series.find(x => x.key === 'energy');
      const all = [{ key: 'total', label: 'Total', values: tot ? tot.values : [] }].concat(c.components);
      const lines = all.map(cp => {
        const last = [...cp.values].reverse().find(isNum);
        return { name: cp.label, pts: steps.map((s, k) => (isNum(cp.values[k]) ? { x: s.n, y: (cp.values[k] - last) * HARTREE.kcal, id: s.frame, mark: k === c.converged,
          label: `${s.label} · ${cp.label}: ${cp.values[k].toFixed(8)} Eh (${((cp.values[k] - last) * HARTREE.kcal).toFixed(3)} kcal/mol vs final)` } : null)).filter(Boolean) };
      });
      return {
        spec: Object.assign(base, { yLabel: 'ΔE vs final (kcal/mol)', lines }),
        csv: z => { const ok = inZ(z); return zoomCSV(head + [unitWord, 'frame'].concat(all.map(cp => cp.key.replace(/[^A-Za-z0-9]+/g, '_') + '_hartree')).join(',') + '\n' +
          steps.map((s, k) => (ok(s.n) ? [s.n, s.frame === null ? '' : s.frame + 1].concat(all.map(cp => csvNum(cp.values[k]))).join(',') : null)).filter(x => x !== null).join('\n') + '\n', z, unitWord + 's'); }
      };
    }
    const ser = c.series.find(x => x.key === key);
    const log = ser.log && d.spec.opt.log;
    const pts = steps.map((s, k) => ({ x: s.n, y: ser.values[k], id: s.frame, mark: k === c.converged, label: `${s.label}${s.frame === null ? ' (geometry not loaded)' : ''} · ${ser.label} = ${fmtCrit(key, ser.values[k])}${ser.threshold ? ' (threshold ' + ser.threshold.toExponential(1) + ')' : ''}` })).filter(p => isNum(p.y) && (!log || p.y > 0));
    return {
      spec: Object.assign(base, { yLabel: ser.label + (ser.unit ? ' (' + ser.unit + ')' : ''), logY: log, points: pts, hlines: ser.threshold ? [{ y: ser.threshold, label: 'threshold ' + ser.threshold.toExponential(1) }] : [] },
        ser.integer ? { integerY: true, yRange: [-0.5, Math.max(1, ...pts.map(p => p.y)) + 0.5] } : {}),
      csv: z => { const ok = inZ(z); return zoomCSV(head + `${unitWord},frame,${key}${ser.unit ? '_' + ser.unit.replace(/[^A-Za-z]+/g, '_') : ''}${ser.threshold ? ',threshold' : ''}\n` + steps.map((s, k) => (ok(s.n) ? [s.n, s.frame === null ? '' : s.frame + 1, csvNum(ser.values[k])].concat(ser.threshold ? [ser.threshold] : []).join(',') : null)).filter(x => x !== null).join('\n') + '\n', z, unitWord + 's'); }
    };
  }
  function profilePlot(d, kind) {
    const p = profFor(d).find(x => x.kind === kind), u = d.spec.energy.unit, E = p.points.map(q => q.energy);
    const ref = d.spec.energy.ref === 'first' ? E[0] : d.spec.energy.ref === 'last' ? E[E.length - 1] : Math.min(...E);
    const pts = p.points.map(q => ({ x: q.x, y: (q.energy - ref) * HARTREE[u], id: q.frame, label: `${q.label} · E = ${q.energy.toFixed(6)} Eh · ΔE = ${((q.energy - ref) * HARTREE[u]).toFixed(2)} ${UNIT_LABEL[u]}` }));
    const xs = pts.map(q => q.x), lo = Math.min(...xs), hi = Math.max(...xs), pad = (hi - lo) * 0.05 || 0.5;
    const integer = /point|step/i.test(p.xLabel);
    return {
      spec: { xRange: integer ? [lo - 0.5, hi + 0.5] : [lo - pad, hi + pad], integerX: integer, xLabel: p.xLabel, yLabel: `ΔE (${UNIT_LABEL[u]})`, points: pts, selected: d.frameIndex, xFormat: x => p.xLabel.replace(/ \(.*\)$/, '') + ' ' + (integer ? Math.round(x) : x.toFixed(3)) },
      csv: z => { const ok = inZ(z); return zoomCSV(`# BlueMolecule ${kind === 'scan' ? 'scan' : 'IRC'} profile: ${d.name}\n${p.xLabel.replace(/[^A-Za-z0-9]+/g, '_').replace(/^_|_$/g, '')},frame,energy_hartree,delta_${u}\n` + p.points.filter(q => ok(q.x)).map(q => `${q.x},${q.frame === null ? '' : q.frame + 1},${q.energy.toFixed(10)},${((q.energy - ref) * HARTREE[u]).toFixed(6)}`).join('\n') + '\n', z, 'points with ' + p.xLabel); }
    };
  }
  function nmrPlot(d) {
    const s = d.spec.nmr, Z = nmrZ(d), nuc = BM.analysis.nmrNuclei(d.data, Z);
    if (!nuc.some(n => n.z === s.z)) s.z = nuc[0].z;
    const lab = BM.analysis.nucleusLabel(s.z), shift = isNum(s.ref) && s.ref !== 0;
    const st = BM.analysis.nmrSticks(d.data, Z, s.z, { ref: shift ? s.ref : 0, merge: s.merge });
    const names = a => a.map(i => el.symbol(Z[i]) + (i + 1)).join(', ');
    const sticks = st.map(p => ({ x: p.x, y: p.count, id: p.atoms[0], atoms: p.atoms, label: `${names(p.atoms)} · σ = ${p.sigma.toFixed(3)} ppm${shift ? ' · δ = ' + p.x.toFixed(3) + ' ppm' : ''}${p.count > 1 ? ' · ×' + p.count : ''}` }));
    const xs = st.map(p => p.x), lo = Math.min(...xs), hi = Math.max(...xs), pad = Math.max((hi - lo) * 0.08, s.z === 1 ? 0.5 : 5);
    const sel = d.selection.length ? d.selection[d.selection.length - 1] : null, selStick = sticks.find(p => p.atoms.indexOf(sel) >= 0);
    return {
      spec: { xRange: [lo - pad, hi + pad], reverseX: shift, xLabel: shift ? `δ(${lab}) (ppm, σref = ${s.ref})` : `Isotropic shielding σ(${lab}) (ppm)`, yLabel: 'Nuclei', sticks, showSticks: true, peakLabels: s.peaks, peakThreshold: 0, maxPeakLabels: 24, peakFormat: x => x.toFixed(2), selected: selStick ? selStick.id : null, xFormat: x => x.toFixed(2) + ' ppm' },
      csv: z => { const ok = inZ(z); return zoomCSV(`# BlueMolecule NMR sticks: ${d.name}; nucleus ${lab}${shift ? '; reference shielding ' + s.ref + ' ppm' : ''}\natoms,shielding_ppm,${shift ? 'shift_ppm,' : ''}count\n` + st.filter(p => ok(p.x)).map(p => `"${names(p.atoms)}",${p.sigma.toFixed(4)},${shift ? p.x.toFixed(4) + ',' : ''}${p.count}`).join('\n') + '\n', z, shift ? 'shifts (ppm)' : 'shieldings (ppm)'); }
    };
  }
  function moLabel(d, id) {
    const cols = BM.analysis.moLevels(d.data), ci = Math.floor(id / 100000), c = cols[ci], l = c && c.levels[id % 100000 - 1];
    if (!l) return '';
    const k = l.index - 1, tag = k === c.homo ? ' (HOMO)' : k === c.lumo ? ' (LUMO)' : '';
    return `${c.label === 'MO' ? 'MO' : c.label} ${l.index}${tag}: ${l.e.toFixed(5)} Eh = ${(l.e * HARTREE.ev).toFixed(3)} eV · ${l.occ ? 'occupied' : 'virtual'}`;
  }
  function moPlot(d) {
    const cols = BM.analysis.moLevels(d.data), EV = HARTREE.ev, w = Math.max(1, d.spec.mo.window || 10);
    let lo = Infinity, hi = -Infinity;
    for (const c of cols) {
      const h = c.homo >= 0 ? c.levels[c.homo].e * EV : null, l = c.lumo >= 0 ? c.levels[c.lumo].e * EV : null;
      const a = isNum(h) ? h : l, b = isNum(l) ? l : h;
      if (isNum(a)) lo = Math.min(lo, a - w);
      if (isNum(b)) hi = Math.max(hi, b + w);
    }
    if (!isFinite(lo)) { lo = -20; hi = 5; }
    const columns = cols.map((c, ci) => ({
      name: c.label === 'MO' ? 'Orbitals' : c.label + ' spin',
      levels: c.levels.filter(l => l.e * EV >= lo - 1 && l.e * EV <= hi + 1).map(l => {
        const k = l.index - 1, tag = k === c.homo ? 'HOMO ' + (l.e * EV).toFixed(2) : k === c.lumo ? 'LUMO ' + (l.e * EV).toFixed(2) : '';
        return { y: l.e * EV, occ: l.occ, id: ci * 100000 + l.index, tag, label: moLabel(d, ci * 100000 + l.index) };
      })
    }));
    const gaps = cols.map((c, ci) => (c.homo >= 0 && c.lumo >= 0 ? { col: ci, from: c.levels[c.homo].e * EV, to: c.levels[c.lumo].e * EV, text: ((c.levels[c.lumo].e - c.levels[c.homo].e) * EV).toFixed(2) + ' eV' } : null)).filter(Boolean);
    return {
      spec: { type: 'levels', columns, gaps, yRange: [lo, hi], yLabel: 'Orbital energy (eV)', selected: d.moSel },
      csv: z => { const ok = inZ(z, 'y'); return zoomCSV(`# BlueMolecule orbital energies: ${d.name}\nspin,index,energy_hartree,energy_eV,occupied\n` + cols.map(c => c.levels.filter(l => ok(l.e * EV)).map(l => `${c.spin},${l.index},${l.e},${(l.e * EV).toFixed(6)},${l.occ ? 1 : 0}`).join('\n')).filter(Boolean).join('\n') + '\n', z, 'orbital energies (eV)', 'y'); }
    };
  }
  function chargePlot(d) {
    const sch = chargeSchemes(d);
    if (sch.indexOf(d.chargeScheme) < 0) d.chargeScheme = sch[0];
    const q = d.data.charges[d.chargeScheme], Z = curFrame(d).Z, lab = i => (Z[i] > 0 ? el.symbol(Z[i]) : 'X') + (i + 1);
    const bars = q.map((v, i) => ({ x: i + 1, y: v, id: i, label: `${lab(i)} · ${d.chargeScheme} q = ${fx(v, 4)}`, color: v >= 0 ? '#e8684a' : '#4f8cff' }));
    return {
      spec: { xRange: [0.4, q.length + 0.6], bars, barWidth: 0.72, xTickLabels: q.map((_, i) => ({ x: i + 1, label: lab(i) })), xLabel: 'Atom', yLabel: `${d.chargeScheme} charge (e)`, selected: d.selection.length ? d.selection[d.selection.length - 1] : null },
      csv: z => { const ok = inZ(z); return zoomCSV(`# BlueMolecule atomic charges: ${d.name}; scheme ${d.chargeScheme}\natom,element,charge\n` + q.map((v, i) => (ok(i + 1) ? `${i + 1},${Z[i] > 0 ? el.symbol(Z[i]) : 'X'},${v}` : null)).filter(x => x !== null).join('\n') + '\n', z, 'atoms'); }
    };
  }
  /** The plot settings that change what an axis means (a unit, a reference, log or linear) or set its range: a zoom
   *  starts again when one of them changes, and is kept through the others (a line width, a scale factor, labels). */
  const ZOOM_SETTINGS = { ir: ['yMode', 'xmin', 'xmax'], raman: ['yMode', 'xmin', 'xmax'], uv: ['unit', 'xmin', 'xmax'], ecd: ['unit', 'xmin', 'xmax'],
    energy: ['unit', 'ref'], opt: ['log'], nmr: ['z', 'ref'], mo: ['window'] };
  /** Names a plot for its zoom (issue #9): the zoom is kept while the same plot is redrawn (a selection, a running
   *  job's new steps, another line width) and starts again for another file, plot, job or axis setting. */
  function zoomKey(d) {
    const t = d.plotType, k = specKey(t), s = d.spec[k];
    return [d.id, t, d.job].concat((ZOOM_SETTINGS[k] || []).map(n => (s ? String(s[n]) : ''))).join('|');
  }
  function renderPlotDraw() {
    const d = state.active;
    if (!d || !d.plotType || d.pending) return;
    let r;
    try { r = buildPlot(d); } catch (e) { console.error(e); r = { spec: { empty: 'This plot could not be drawn: ' + e.message }, csv: null }; }
    d._plotCsv = r.csv || null;
    r.spec.zoomKey = zoomKey(d);
    plot.set(r.spec);
  }
  function plotCSV() {
    const d = state.active; if (!d || !d.plotType) return;
    let r;
    try { r = buildPlot(d); } catch (e) { toast(e.message, true); return; }
    if (r.csv) saveFile(`${baseName(d)}_${d.plotType.replace(':', '_')}.csv`, r.csv(plot.zoomed()), 'text/csv');
  }
  /** The zoom buttons show whether the plot is zoomed: the exports then hold only what is shown. */
  function onPlotZoom(zoomed) {
    const fit = $('#btn-plot-fit'), out = $('#btn-plot-zoom-out');
    if (fit) { fit.classList.toggle('on', zoomed); fit.disabled = !zoomed; }
    if (out) out.disabled = !zoomed;
  }

  // ---------------- dialogs ----------------
  let dlgSeq = 0;
  function dialog(title, body, buttons, onOpen) {
    const dlg = $('#dlg'), form = $('#dlg-form'), seq = ++dlgSeq;
    form.innerHTML = `<div class="dlg-head">${title}</div><div class="dlg-body">${body}</div><div class="dlg-foot">${buttons.map((b, i) => `<button type="button" class="btn ${b.primary ? 'primary' : ''}" data-dbtn="${i}">${b.label}</button>`).join('')}</div>`;
    $$('[data-dbtn]', form).forEach(btn => btn.addEventListener('click', async () => {
      const b = buttons[+btn.dataset.dbtn];
      const keep = b.onClick ? await b.onClick(form) : false;
      if (!keep && seq === dlgSeq) dlg.close(); // a newer dialog may have replaced this one meanwhile
    }));
    if (!dlg.open) dlg.showModal();
    if (!form.querySelector('.dlg-body input, .dlg-body select, .dlg-body textarea')) { const pb = form.querySelector('.dlg-foot .btn.primary'); if (pb) pb.focus(); }
    if (onOpen) onOpen(form);
  }
  function geomFor(d) { const f = curFrame(d); const xyz = d._animXyz && state.anim.kind === 'vib' ? d._animXyz : f.xyz; return { Z: f.Z, xyz, f }; }
  const TASK_NAMES = { sp: 'Single point', opt: 'Optimization', freq: 'Frequencies', optfreq: 'Optimization + frequencies' };
  function loadCalcPrefs() {
    const def = { program: 'gaussian', task: 'optfreq', nprocs: 8, memoryMB: 16000, methods: {}, bases: {} };
    try { const p = JSON.parse(localStorage.getItem('bluemolecule.calc') || '{}'); return Object.assign(def, p, { methods: Object.assign({}, p.methods), bases: Object.assign({}, p.bases) }); } catch (e) { return def; }
  }
  function saveCalcPrefs(p) { try { localStorage.setItem('bluemolecule.calc', JSON.stringify(p)); } catch (e) { /* ignore */ } }
  /** Charge and multiplicity defaults: the file's values; for built molecules the suggested charge (formal charges,
   *  plus ligand charges and metal oxidation states for complexes) and multiplicity 1. */
  function chargeDefaults(d) {
    const m = d.data.meta || {};
    if (!d.built && isNum(m.charge)) return { charge: m.charge, mult: isNum(m.multiplicity) ? m.multiplicity : 1, note: '' };
    const t = topoFrame(d), s = kit.suggestCharge(t), ox = kit.oxText(s.metals);
    if (!s.determinable) return { charge: t.formalCharges.reduce((a, q) => a + (q || 0), 0), mult: 1, note: 'Check the charge: ' + ox + '.' };
    return { charge: s.total, mult: 1, note: s.metals.length ? 'Charge suggested from formal and ligand charges' + (ox ? ' (assuming ' + ox + ')' : '') + '; please check it.' : '' };
  }
  function parityNote(Z, charge, mult) {
    let e = -charge;
    for (const z of Z) if (z > 0) e += z;
    return e < mult - 1 || (((e - (mult - 1)) % 2) + 2) % 2 ? 'Check charge and multiplicity: ' + e + ' electrons cannot have multiplicity ' + mult + '.' : '';
  }
  /** The calculation-input dialog (input files for BM.inputPrograms; programId preselects a program). pair = {Z, minus, plus, amount, mode, freq} writes two inputs instead, for the
   *  geometries displaced along a vibration (forward and reverse optimisations from a transition state, issue #14):
   *  name_minus and name_plus, the Task set to Optimization and a TS route or keyword line turned into an optimisation. */
  function exportCalc(programId, pair) {
    const d = state.active; if (!d) return toast('Open a file first.', true);
    const g = pair ? { Z: pair.Z, xyz: pair.plus, f: null } : geomFor(d);
    if (!g.Z.some(z => z > 0)) return toast('There are no atoms to write yet.', true);
    const progs = BM.inputPrograms, prefs = loadCalcPrefs(), m = d.data.meta || {};
    let pid = programId || (d.data.format === 'orca-out' ? 'orca' : /^gaussian/.test(d.data.format) ? 'gaussian' : prefs.program);
    if (!progs.some(p => p.id === pid)) pid = 'gaussian';
    // the route / keywords of the file; geom= and guess= would read the old geometry from a checkpoint file
    const fileRoute = /^gaussian/.test(d.data.format) && m.route
      ? (pair ? BM.vibtools.minimumRoute(m.route) : String(m.route).replace(/\b(geom|geometry|guess)\s*=\s*(\([^)]*\)|\S+)/gi, ' ').replace(/\s+/g, ' ').trim()) : null;
    const fileKw = d.data.format === 'orca-out' && m.keywords
      ? (pair ? BM.vibtools.minimumKeywords(m.keywords) : Array.isArray(m.keywords) ? m.keywords.join(' ') : String(m.keywords)) : null;
    const cd = chargeDefaults(d);
    // the Job name is the title of the input and, made safe, the name of its files (issue #3)
    const jobName = d.built ? d.name : baseName(d);
    const body = `
      <div class="two"><label class="field"><span>Program</span><select data-c="program">${progs.map(p => `<option value="${p.id}" ${p.id === pid ? 'selected' : ''}>${esc(p.name)}</option>`).join('')}</select></label>
        <label class="field"><span>Task</span><select data-c="task"></select></label></div>
      ${fileRoute ? `<div data-c-row="gaussian" class="crow"><label class="field cw"><span></span><span class="inline"><input type="checkbox" data-c="useRoute" checked> Use the route from this file${pair ? ', as an optimisation' : ''}</span></label>
        <label class="field cw" data-c-row="route"><span>Route</span><input type="text" data-c="route" value="${esc(fileRoute)}" spellcheck="false"></label></div>` : ''}
      ${fileKw ? `<div data-c-row="orca" class="crow"><label class="field cw"><span></span><span class="inline"><input type="checkbox" data-c="useKw" checked> Use the keywords from this file${pair ? ', as an optimisation' : ''}</span></label>
        <label class="field cw" data-c-row="kw"><span>Keywords</span><input type="text" data-c="kw" value="${esc(fileKw)}" spellcheck="false"></label></div>` : ''}
      <div class="two"><label class="field"><span>Method</span><input type="text" data-c="method" list="calc-methods" spellcheck="false"></label>
        <label class="field"><span>Basis set</span><input type="text" data-c="basis" list="calc-bases" spellcheck="false"></label></div>
      <datalist id="calc-methods"></datalist><datalist id="calc-bases"></datalist>
      <div class="two"><label class="field"><span>Charge</span><input type="number" data-c="charge" value="${cd.charge}" step="1"></label>
        <label class="field"><span>Multiplicity</span><input type="number" data-c="mult" value="${cd.mult}" min="1" step="1"></label></div>
      <div class="two"><label class="field"><span>Cores</span><input type="number" data-c="nprocs" value="${prefs.nprocs}" min="1" step="1"></label>
        <label class="field"><span>Memory (MB)</span><input type="number" data-c="mem" value="${prefs.memoryMB}" min="0" step="500"></label></div>
      <div class="two"><label class="field"><span>Solvent</span><input type="text" data-c="solvent" placeholder="none (gas phase)" spellcheck="false"></label>
        <label class="field" title="Used as the file name (and the checkpoint name) and as the title of the input"><span>Job name</span><input type="text" data-c="title" value="${esc(jobName)}" spellcheck="false"></label></div>
      <label class="field cw"><span>Extra keywords</span><input type="text" data-c="extra" placeholder="optional; written where the program expects extra input" spellcheck="false"></label>
      <textarea data-c="preview" spellcheck="false"${pair ? ' readonly' : ''}></textarea>
      <p class="note" data-c="note"></p>`;
    let result = null;
    const V = (form, k) => $(`[data-c="${k}"]`, form);
    const prog = form => progs.find(p => p.id === V(form, 'program').value) || progs[0];
    const usesFile = (form, p) => (p.id === 'gaussian' && fileRoute && V(form, 'useRoute').checked) || (p.id === 'orca' && fileKw && V(form, 'useKw').checked);
    const baseOf = form => BM.writers.fileBase(V(form, 'title').value, 'job');
    /** One input: {text, files?, command?, note?}. */
    const make = (form, xyz, base, title) => {
      const p = prog(form), n = s => V(form, s).value;
      const o = { Z: g.Z, xyz, charge: parseInt(n('charge'), 10) || 0, multiplicity: parseInt(n('mult'), 10) || 1, method: n('method').trim(), basis: n('basis').trim(), task: n('task'),
        nprocs: Math.max(1, parseInt(n('nprocs'), 10) || 1), memoryMB: Math.max(0, parseInt(n('mem'), 10) || 0), title, solvent: n('solvent').trim() || undefined, extra: n('extra').trim() || undefined, baseName: base };
      if (p.id === 'gaussian' && usesFile(form, p)) {
        const link0 = ['%chk=' + base + '.chk'];
        if (o.nprocs > 1) link0.push('%nprocshared=' + o.nprocs);
        if (o.memoryMB) link0.push('%mem=' + o.memoryMB + 'MB');
        let text = BM.writers.gaussianInput({ Z: g.Z, xyz, route: n('route'), title, charge: o.charge, multiplicity: o.multiplicity, link0, cell: cellOf(d) });
        if (o.extra) text = text.replace(/\n\n\n$/, '\n\n' + o.extra + '\n\n');
        return { text, command: 'g16 ' + base + '.gjf', note: parityNote(g.Z, o.charge, o.multiplicity) };
      }
      if (p.id === 'orca' && usesFile(form, p)) {
        const text = BM.writers.orcaInput({ Z: g.Z, xyz, keywords: n('kw'), charge: o.charge, multiplicity: o.multiplicity, nprocs: o.nprocs,
          maxcore: o.memoryMB ? Math.max(1, Math.floor(o.memoryMB / o.nprocs)) : null, extraBlocks: o.extra ? [o.extra] : [] });
        return { text, command: 'orca ' + base + '.inp > ' + base + '.out', note: parityNote(g.Z, o.charge, o.multiplicity) };
      }
      return p.write(o);
    };
    /** The files of both inputs of a pair: [{name, text}]; files both share (an xcontrol) once, other names that
     *  clash (Turbomole's coord) get _minus / _plus. */
    const pairFiles = (p, r, base) => {
      const list = (x, b) => (x.files && x.files.length ? x.files.map(q => ({ name: q.name, text: q.text })) : [{ name: b + '.' + (p.ext || 'inp'), text: x.text }]);
      const a = list(r.minus, base + '_minus'), b = list(r.plus, base + '_plus'), out = a.slice();
      for (const q of b) {
        const k = out.findIndex(x => x.name === q.name);
        if (k < 0) out.push(q);
        else if (out[k].text !== q.text) { out[k] = { name: q.name + '_minus', text: out[k].text }; out.push({ name: q.name + '_plus', text: q.text }); }
      }
      return out;
    };
    const build = form => {
      const p = prog(form), base = baseOf(form), title = V(form, 'title').value, notes = [cd.note];
      try {
        if (pair) {
          const nu = pair.freq < 0 ? Math.abs(pair.freq).toFixed(1) + 'i' : pair.freq.toFixed(1), what = `mode ${pair.mode + 1} (${nu} cm⁻¹)`;
          // input titles in plain ASCII (programs read them as such)
          const t = s => `${title}, mode ${pair.mode + 1} (${nu} cm-1) displaced by ${s}${pair.amount} A`;
          result = { minus: make(form, pair.minus, base + '_minus', t('-')), plus: make(form, pair.plus, base + '_plus', t('+')) };
          result.files = pairFiles(p, result, base);
          V(form, 'preview').value = result.files.map(x => `----- ${x.name} -----\n${x.text}`).join('\n');
          notes.push(result.plus.note);
          notes.push(`Saves ${result.files.map(x => x.name).join(', ')}: the geometry displaced by −${pair.amount} and +${pair.amount} Å (largest atom shift) along ${what}.`);
          if (result.plus.command) notes.push('Run: <code>' + esc(result.minus.command) + '</code> and <code>' + esc(result.plus.command) + '</code>');
        } else {
          result = make(form, g.xyz, base, title);
          V(form, 'preview').value = result.text;
          notes.push(result.note);
          if (result.files && result.files.length > 1) notes.push('Saves ' + result.files.map(x => x.name).join(' and ') + '.');
          else if (result.files && result.files.length) notes.push('Saves ' + result.files[0].name + '.');
          else notes.push('Saves ' + base + '.' + (p.ext || 'inp') + '.');
          if (result.command) notes.push('Run: <code>' + esc(result.command) + '</code>');
        }
        V(form, 'note').innerHTML = notes.filter(Boolean).map(x => (/^Run: /.test(x) ? x : esc(x))).join(' ') || `The preview is editable. Geometry ${d.frameIndex + 1} of ${d.data.frames.length} is used.`;
        V(form, 'note').classList.remove('err');
      } catch (e) {
        result = null;
        V(form, 'preview').value = '';
        V(form, 'note').textContent = e.message || String(e);
        V(form, 'note').classList.add('err');
      }
    };
    const setProgram = form => {
      const p = prog(form), task = V(form, 'task'), keep = task.value || (pair ? 'opt' : prefs.task), mi = V(form, 'method'), bi = V(form, 'basis');
      if (form._pid !== p.id) {
        form._pid = p.id;
        task.innerHTML = p.tasks.map(t => `<option value="${t}" ${t === (p.tasks.indexOf(keep) >= 0 ? keep : 'optfreq') ? 'selected' : ''}>${TASK_NAMES[t] || t}</option>`).join('');
        mi.value = prefs.methods[p.id] !== undefined ? prefs.methods[p.id] : p.defaultMethod || '';
        bi.value = prefs.bases[p.id] !== undefined ? prefs.bases[p.id] : p.defaultBasis || '';
      }
      $('#calc-methods', form).innerHTML = (p.methods || []).map(x => `<option value="${esc(x)}">`).join('');
      $('#calc-bases', form).innerHTML = (p.bases || []).map(x => `<option value="${esc(x)}">`).join('');
      const noMethod = !p.defaultMethod && !(p.methods || []).length, noBasis = !p.defaultBasis && !(p.bases || []).length;
      mi.placeholder = noMethod ? 'set later in the program' : ''; bi.placeholder = noBasis ? 'not used' : '';
      const fileMode = usesFile(form, p);
      mi.disabled = noMethod || fileMode; bi.disabled = noBasis || fileMode; task.disabled = fileMode; V(form, 'solvent').disabled = fileMode;
      for (const id of ['gaussian', 'orca']) { const row = $(`[data-c-row="${id}"]`, form); if (row) row.classList.toggle('hidden', p.id !== id); }
      for (const [id, chk] of [['route', 'useRoute'], ['kw', 'useKw']]) { const row = $(`[data-c-row="${id}"]`, form), c = V(form, chk); if (row && c) row.classList.toggle('hidden', !c.checked); }
      const btn = $('[data-dbtn="2"]', form);
      if (btn) btn.textContent = pair ? 'Save both' : p.ext ? 'Save .' + p.ext : (p.id === 'turbomole' ? 'Save coord' : 'Save');
      build(form);
    };
    const remember = form => {
      const p = prog(form);
      prefs.program = p.id;
      if (!pair) prefs.task = V(form, 'task').value || prefs.task;   // the pair's Optimization is not a new default
      prefs.nprocs = parseInt(V(form, 'nprocs').value, 10) || prefs.nprocs; prefs.memoryMB = parseInt(V(form, 'mem').value, 10) || 0;
      if (!V(form, 'method').disabled) prefs.methods[p.id] = V(form, 'method').value.trim();
      if (!V(form, 'basis').disabled) prefs.bases[p.id] = V(form, 'basis').value.trim();
      saveCalcPrefs(prefs);
    };
    dialog(pair ? 'Forward and reverse inputs' : 'Calculation input', body, [
      { label: 'Copy', onClick: form => { copyText(V(form, 'preview').value, pair ? 'Both inputs copied to the clipboard' : 'Input copied to the clipboard'); return true; } },
      { label: 'Cancel' },
      {
        label: 'Save', primary: true, onClick: form => {
          if (!result) { toast('Fix the problem shown below the preview first.', true); return true; }
          const p = prog(form), text = V(form, 'preview').value, files = result.files, command = pair ? null : result.command;
          remember(form);
          if (files && files.length) {
            (async () => {
              for (let k = 0; k < files.length; k++) await saveFile(files[k].name, k === 0 && !pair ? text : files[k].text);
              toast('Saved ' + files.map(x => x.name).join(', ') + (command ? ' · run: ' + command : ''));
            })();
          } else saveFile(baseOf(form) + '.' + (p.ext || 'inp'), text);
        }
      }
    ], form => {
      setProgram(form);
      // listen on the dialog body, which dialog() rebuilds each time: listeners on the shared form would pile up
      // and keep running for later dialogs with this document's (stale) geometry and route/keyword fields
      const dbody = form.querySelector('.dlg-body') || form;
      dbody.addEventListener('input', e => { const k = e.target.dataset.c; if (!k || k === 'preview') return; if (k === 'program' || k === 'useRoute' || k === 'useKw') setProgram(form); else build(form); });
      dbody.addEventListener('change', e => { const k = e.target.dataset.c; if (k === 'program' || k === 'useRoute' || k === 'useKw') setProgram(form); });
    });
  }
  /** Forward and reverse inputs: the geometry of the frequency job displaced by -/+ the Displace amount along the
   *  selected vibration, each as an optimisation input. */
  function vibInputs() {
    const d = state.active, v = d && d.data.vibrations;
    if (!v || d.vibIndex < 0) return toast('Select a vibration first (Vibrations tab).', true);
    const box = $('#disp-amt'), amount = (box && parseFloat(box.value) > 0 ? parseFloat(box.value) : 0) || state.vib.disp || 0.15;
    const mode = modeOf(d, d.vibIndex);
    if (!mode) return toast(NOT_LOADED, true);
    const eq = d.data.frames[v.frameIndex], e = BM.vibtools.ends(eq.xyz, mode, amount);
    exportCalc(null, { Z: eq.Z, minus: e.minus, plus: e.plus, amount, mode: d.vibIndex, freq: v.freqs[d.vibIndex] });
  }
  function exportGaussian() { exportCalc('gaussian'); }
  function exportOrca() { exportCalc('orca'); }
  // ---------------- image export (issue #16) ----------------
  const IMG_UNITS = { px: 1, cm: 2.54, in: 1 };
  function loadImagePrefs() {
    const def = { format: 'png', unit: 'px', width: 0, dpi: 300, transparent: false, quality: 92, level: 6, profile: 'srgb' };
    try { return Object.assign(def, JSON.parse(localStorage.getItem('bluemolecule.image') || '{}')); } catch (e) { return def; }
  }
  function saveImagePrefs(p) { try { localStorage.setItem('bluemolecule.image', JSON.stringify(p)); } catch (e) { /* ignore */ } }
  const MAX_IMAGE_PIXELS = 250e6;
  /** Save image: PNG, JPEG, WebP, TIFF or PDF of the current view (labels included) at a size in pixels or a print
   *  size, with BlueMolecule's own encoders (BM.imagex); large images are drawn in tiles. */
  function exportImage() {
    const d = state.active; if (!d) return toast('Open a file first.', true);
    const cw = Math.max(1, $('#gl').clientWidth), ch = Math.max(1, $('#gl').clientHeight), aspect = ch / cw;
    const p = loadImagePrefs(), F = BM.imagex.FORMATS;
    if (!F.some(f => f.id === p.format)) p.format = 'png';
    if (!IMG_UNITS[p.unit]) p.unit = 'px';
    const pxOf = (v, unit, dpi) => (unit === 'px' ? v : unit === 'cm' ? v / 2.54 * dpi : v * dpi);
    const ofPx = (px, unit, dpi) => (unit === 'px' ? px : unit === 'cm' ? px / dpi * 2.54 : px / dpi);
    const w0 = p.width > 0 ? p.width : ofPx(2 * cw, p.unit, p.dpi);   // the first time: twice the view
    const num = (v, unit) => (unit === 'px' ? String(Math.round(v)) : String(Math.round(v * 100) / 100));
    const body = `
      <div class="two"><label class="field"><span>Format</span><select data-i="format">${F.map(f => `<option value="${f.id}" ${f.id === p.format ? 'selected' : ''}>${esc(f.name)}</option>`).join('')}</select></label>
        <label class="field"><span>Colour profile</span><select data-i="profile"><option value="srgb" ${p.profile === 'srgb' ? 'selected' : ''}>sRGB (recommended)</option><option value="p3" ${p.profile === 'p3' ? 'selected' : ''}>Display P3</option><option value="" ${!p.profile ? 'selected' : ''}>None</option></select></label></div>
      <div class="two"><label class="field"><span>Width</span><span class="inline"><input type="number" data-i="width" min="0.01" step="any" value="${num(w0, p.unit)}" style="width:96px">
          <select data-i="unit"><option value="px" ${p.unit === 'px' ? 'selected' : ''}>pixels</option><option value="cm" ${p.unit === 'cm' ? 'selected' : ''}>cm</option><option value="in" ${p.unit === 'in' ? 'selected' : ''}>inches</option></select></span></label>
        <label class="field"><span>Resolution</span><span class="inline"><input type="number" data-i="dpi" min="10" max="4800" step="1" value="${p.dpi}" style="width:80px"> DPI</span></label></div>
      <div class="two"><label class="field"><span>Background</span><span class="inline"><input type="checkbox" data-i="transparent" ${p.transparent ? 'checked' : ''}> transparent</span></label>
        <label class="field" data-i-row="quality"><span>Quality</span><span class="inline"><input type="range" data-i="quality" min="1" max="100" step="1" value="${p.quality}"><output>${p.quality}</output></span></label>
        <label class="field" data-i-row="level"><span>Compression</span><select data-i="level"><option value="0" ${p.level === 0 ? 'selected' : ''}>none (fastest)</option><option value="1" ${p.level === 1 ? 'selected' : ''}>fast</option><option value="6" ${p.level === 6 ? 'selected' : ''}>normal</option><option value="9" ${p.level === 9 ? 'selected' : ''}>best (slowest)</option></select></label></div>
      <p class="note" data-i="size"></p>
      <p class="note">The image is the current view with its labels and legend, in the proportions of the view (resize the window to change them). Choose a white background in the toolbar for publications.</p>`;
    let plan = null;
    const V = (form, k) => $(`[data-i="${k}"]`, form);
    const update = form => {
      const f = F.find(x => x.id === V(form, 'format').value) || F[0], unit = V(form, 'unit').value, dpi = Math.max(10, Math.min(4800, parseFloat(V(form, 'dpi').value) || 300));
      const W = Math.round(pxOf(parseFloat(V(form, 'width').value) || 0, unit, dpi)), H = Math.round(W * aspect);
      const tr = V(form, 'transparent');
      tr.disabled = !f.alpha;
      $('[data-i-row="quality"]', form).classList.toggle('hidden', f.id !== 'jpeg');
      $('[data-i-row="level"]', form).classList.toggle('hidden', f.id === 'jpeg' || f.id === 'webp');
      const cmW = W / dpi * 2.54, cmH = H / dpi * 2.54, mp = W * H / 1e6;
      let err = '';
      if (!(W >= 16 && H >= 16)) err = 'The image must be at least 16 pixels wide and high.';
      else if (Math.max(W, H) > f.max) err = `${f.name} images are at most ${f.max.toLocaleString('en-US')} pixels wide and high.`;
      else if (W * H > MAX_IMAGE_PIXELS) err = `At most ${MAX_IMAGE_PIXELS / 1e6} million pixels (this is ${Math.round(mp)} million).`;
      const tiles = Math.ceil(W / 2048) * Math.ceil(H / 2048);
      const el2 = V(form, 'size');
      el2.textContent = err || `${W.toLocaleString('en-US')} × ${H.toLocaleString('en-US')} pixels (${mp < 10 ? mp.toFixed(1) : Math.round(mp)} megapixels) · ${cmW.toFixed(1)} × ${cmH.toFixed(1)} cm (${(cmW / 2.54).toFixed(2)} × ${(cmH / 2.54).toFixed(2)} in) at ${Math.round(dpi)} DPI${tiles > 1 ? ` · drawn in ${tiles} tiles` : ''}`;
      el2.classList.toggle('err', !!err);
      const btn = $('[data-dbtn="1"]', form);
      if (btn) btn.textContent = 'Save .' + f.ext;
      plan = err ? null : { f, W, H, dpi, unit, transparent: f.alpha && tr.checked, quality: parseInt(V(form, 'quality').value, 10) || 92, level: parseInt(V(form, 'level').value, 10), profile: V(form, 'profile').value || null };
    };
    dialog('Save image', body, [{ label: 'Cancel' }, {
      label: 'Save', primary: true, onClick: form => {
        update(form);
        if (!plan) { toast(V(form, 'size').textContent, true); return true; }
        const q = plan;
        saveImagePrefs({ format: q.f.id, unit: q.unit, width: parseFloat(V(form, 'width').value), dpi: q.dpi, transparent: V(form, 'transparent').checked, quality: q.quality, level: q.level, profile: q.profile || '' });
        renderAndSave(d, q);
        return false;
      }
    }], form => {
      update(form);
      const body2 = form.querySelector('.dlg-body') || form;
      body2.addEventListener('input', e => {
        const k = e.target.dataset.i;
        if (k === 'quality' && e.target.nextElementSibling) e.target.nextElementSibling.textContent = e.target.value;
        update(form);
      });
      let lastUnit = V(form, 'unit').value;
      body2.addEventListener('change', e => {
        const k = e.target.dataset.i;
        if (k === 'unit') {   // keep the same number of pixels in the new unit
          const dpi = parseFloat(V(form, 'dpi').value) || 300, px = pxOf(parseFloat(V(form, 'width').value) || 0, lastUnit, dpi), u = V(form, 'unit').value;
          V(form, 'width').value = num(ofPx(px, u, dpi), u);
          lastUnit = u;
        }
        update(form);
      });
    });
  }
  /** Draw the view at the planned size and save it in the planned format. */
  async function renderAndSave(d, q) {
    const tick = () => new Promise(r => setTimeout(r, 0));
    try {
      status(`Drawing the image (${q.W} × ${q.H} pixels)…`);
      await tick();
      const img = await R.renderImage({ width: q.W, height: q.H, transparent: q.transparent, onProgress: (k, n) => { if (n > 1) status(`Drawing the image: tile ${k} of ${n}…`); return tick(); } });
      status(`Writing the ${q.f.name} file…`);
      await tick();
      const bytes = BM.imagex.encode(q.f.id, img, { dpi: q.dpi, profile: q.profile, quality: q.quality, level: q.level, alpha: q.transparent, title: d.name });
      await saveFile(baseName(d) + '.' + q.f.ext, new Blob([bytes], { type: q.f.mime }), q.f.mime);
    } catch (e) {
      toast('Could not make the image: ' + (e && e.message ? e.message : e), true);
    }
  }
  function exportSimple(kind) {
    const d = state.active; if (!d) return toast('Open a file first.', true);
    const g = geomFor(d), f = g.f;
    try {
      if (kind === 'xyz') saveFile(baseName(d) + (f.edited ? '_edited' : '') + '.xyz', BM.writers.xyz({ Z: g.Z, xyz: g.xyz, comment: `${d.name} geometry ${d.frameIndex + 1}${isNum(f.energy) ? ' E = ' + f.energy.toFixed(10) : ''}` }));
      else if (kind === 'xyz-all') saveFile(baseName(d) + '_all.xyz', BM.writers.xyzMulti(d.data.frames.map((x, i) => ({ Z: x.Z, xyz: x.xyz, comment: `${d.name} geometry ${i + 1}${isNum(x.energy) ? ' E = ' + x.energy.toFixed(10) : ''}` }))));
      else if (kind === 'mol') saveFile(baseName(d) + '.mol', BM.writers.mol({ Z: g.Z, xyz: g.xyz, bonds: bondPairs(d), title: baseName(d) }));
      else if (kind === 'pdb') saveFile(baseName(d) + '.pdb', BM.writers.pdb({ Z: g.Z, xyz: g.xyz, atomNames: f.atomNames, resNames: f.resNames, resIds: f.resIds, cell: cellOf(d) }));
    } catch (e) { toast(e.message, true); }
  }
  /** Bonds of the current frame with orders (explicit, from the file, or perceived) for MOL export. */
  function bondPairs(d) {
    return topoFrame(d).bonds.map(b => [b[0], b[1], b[2] === undefined ? 1 : b[2]]);
  }
  function exportVibration() {
    const d = state.active, v = d && d.data.vibrations;
    if (!v || d.vibIndex < 0) return toast('Select a vibration first (Vibrations tab).', true);
    const eq = d.data.frames[v.frameIndex], u = vibVec(d, d.vibIndex, state.vib.amplitude), N = 20, frames = [];
    if (!u) return toast(NOT_LOADED, true);
    const f = v.freqs[d.vibIndex];
    for (let k = 0; k < N; k++) {
      const s = Math.sin(2 * Math.PI * k / N);
      frames.push({ Z: eq.Z, xyz: Float64Array.from(eq.xyz, (x, i) => x + s * u[i]), comment: `mode ${d.vibIndex + 1} ${f.toFixed(2)} cm-1 frame ${k + 1}/${N}` });
    }
    saveFile(`${baseName(d)}_mode${d.vibIndex + 1}.xyz`, BM.writers.xyzMulti(frames));
  }
  // ---------------- Save As ----------------
  const SA_GROUPS = [
    ['Molecule', ['xyz', 'extxyz', 'mol', 'sdf', 'mol2', 'cml', 'cjson', 'pdb', 'smiles']],
    ['Crystal / MD', ['cif', 'poscar', 'qe-in', 'cp2k-in', 'xsf', 'gro', 'lammps-data']],
    ['Quantum chemistry', ['molden', 'cube', 'gaussian-input', 'orca-input', 'turbomole-coord']]
  ];
  const SA_PAIRS = { xyz: 'xyz-traj', sdf: 'sdf-frames', extxyz: 'extxyz-frames' };
  const SA_NEEDS = { cell: 'needs a unit cell', volume: 'needs volume data', vibrations: 'needs vibrations', frames: 'needs several geometries', mos: 'needs orbitals' };
  function saveNeeds(d, f) {
    for (const n of f.needs || []) {
      if (n === 'cell' && !cellOf(d)) return SA_NEEDS.cell;
      if (n === 'volume' && !d.surfaces.length) return SA_NEEDS.volume;
      if (n === 'vibrations' && !d.data.vibrations) return SA_NEEDS.vibrations;
      if (n === 'frames' && visibleFrames(d).length < 2) return SA_NEEDS.frames;
      if (n === 'mos' && !d.data.mos) return SA_NEEDS.mos;
    }
    return null;
  }
  /** Writer context (see writers.js) for the current geometry, or for all shown geometries. */
  function saveContext(d, all) {
    const g = geomFor(d), f = g.f, topo = topoFrame(d), fr = d.data.frames, vis = visibleFrames(d);
    const list = all ? vis.map(i => fr[i]) : fr;
    const surf = d.surfaces.find(x => x.visible) || d.surfaces[0];
    return {
      Z: g.Z, xyz: g.xyz, bonds: topo.bonds.map(b => [b[0], b[1], b[2] === undefined ? 1 : b[2]]), formalCharges: topo.formalCharges,
      cell: cellOf(d, f), title: d.built ? d.name : (d.data.title || baseName(d)), atomNames: f.atomNames, resNames: f.resNames, resIds: f.resIds, chains: f.chains || null,
      frames: list.map(x => ({ Z: x.Z, xyz: x.xyz, energy: x.energy, label: x.label, cell: x.cell || null })), frameIndex: all ? Math.max(0, vis.indexOf(d.frameIndex)) : d.frameIndex,
      volume: surf ? surf.vol : null, vibrations: d.data.vibrations, charges: d.data.charges, meta: d.data.meta, baseName: baseName(d), mos: d.data.mos || null,
      orbitalModel: BM.orbitals ? BM.orbitals.model(d.data) : null
    };
  }
  function runWriter(d, id, all) {
    const f = (BM.saveFormats || []).find(x => x.id === id);
    if (!f) throw new Error('Unknown format ' + id);
    return f.write(saveContext(d, all));
  }
  function saveAs() {
    const d = state.active; if (!d) return toast('Open a file first.', true);
    if (d.pending) return toast('This topology has no coordinates yet.', true);
    const fmts = BM.saveFormats || [], byId = id => fmts.find(x => x.id === id);
    const trajIds = new Set(Object.values(SA_PAIRS)), listed = new Set(), groups = [];
    for (const [g, ids] of SA_GROUPS) { const items = ids.map(byId).filter(Boolean); items.forEach(x => listed.add(x.id)); if (items.length) groups.push([g, items]); }
    const rest = fmts.filter(x => !listed.has(x.id) && !trajIds.has(x.id));
    if (rest.length) groups.push(['Other', rest]);
    let prefs = {};
    try { prefs = JSON.parse(localStorage.getItem('bluemolecule.saveas') || '{}'); } catch (e) { /* ignore */ }
    const nfr = visibleFrames(d).length;
    let fmt = prefs.fmt && (byId(prefs.fmt) || /^ob:/.test(prefs.fmt)) && !(byId(prefs.fmt) && saveNeeds(d, byId(prefs.fmt))) ? prefs.fmt : 'xyz';
    if (/^ob:/.test(fmt)) fmt = 'xyz';   // Open Babel formats are listed once obabel answers
    const opt = x => { const why = saveNeeds(d, x); return `<option value="${x.id}" ${why ? 'disabled' : ''} ${x.id === fmt ? 'selected' : ''}>${esc(x.name)} · .${x.ext}${why ? ' — ' + why : ''}</option>`; };
    const body = `<div class="sa">
        <label class="field"><span>Format</span><select data-sa="fmt">${groups.map(([g, items]) => `<optgroup label="${g}">${items.map(opt).join('')}</optgroup>`).join('')}</select></label>
        <div class="field" data-sa="frrow"><span>Geometries</span><span class="sa-seg"><label><input type="radio" name="sa-fr" value="cur" ${prefs.all ? '' : 'checked'}> current (${d.frameIndex + 1})</label><label><input type="radio" name="sa-fr" value="all" ${prefs.all ? 'checked' : ''}> all ${nfr}</label></span></div>
        <label class="field"><span>File name</span><input type="text" data-sa="name" spellcheck="false"></label>
        <div class="field"><span>Open Babel</span><span class="dim" data-sa="ob">checking…</span></div>
      </div>
      <textarea data-sa="preview" class="sa-prev" readonly spellcheck="false" wrap="off"></textarea>
      <p class="note" data-sa="note"></p>`;
    let nameEdited = false;
    const cur = form => {
      const id = $('[data-sa=fmt]', form).value, ob = /^ob:/.test(id);
      const all = nfr > 1 && $('input[name=sa-fr]:checked', form).value === 'all';
      const base = ob ? 'sdf' : id, pairOk = !!SA_PAIRS[base];
      return { id, ob, code: ob ? id.slice(3) : null, all: all && (pairOk || ob), writer: ob ? (all ? 'sdf-frames' : 'sdf') : (all && pairOk ? SA_PAIRS[base] : id), pairOk };
    };
    const update = form => {
      const c = cur(form), f = byId(c.writer) || byId(c.id);
      $('[data-sa=frrow]', form).classList.toggle('hidden', !(nfr > 1 && (c.pairOk || c.ob)));
      const ext = c.ob ? c.code : (byId(c.id) || f).ext;
      const nm = $('[data-sa=name]', form);
      if (!nameEdited) nm.value = baseName(d).replace(/[^A-Za-z0-9_.+-]+/g, '_') + (c.all ? '_all' : '') + '.' + ext;
      const note = $('[data-sa=note]', form), pv = $('[data-sa=preview]', form);
      try {
        const out = runWriter(d, c.writer, c.all);
        const text = typeof out === 'string' ? out : out && out.files ? out.files.map(x => '# ' + x.name + '\n' + x.text).join('\n') : '';
        const lines = text.split('\n');
        pv.value = lines.slice(0, 40).join('\n') + (lines.length > 41 ? `\n… ${lines.length - 40} more lines` : '');
        note.textContent = c.ob ? `Open Babel converts this SD file to “${c.code}” when saving.` : (f && f.kind === 'volume' ? 'Writes the visible volume.' : '');
        note.classList.remove('err');
      } catch (e) { pv.value = ''; note.textContent = e.message; note.classList.add('err'); }
    };
    dialog('Save as', body, [{ label: 'Cancel' }, {
      label: 'Save', primary: true, onClick: async form => {
        const c = cur(form), name = $('[data-sa=name]', form).value.trim();
        if (!name) { toast('Enter a file name.', true); return true; }
        try {
          let out = runWriter(d, c.writer, c.all);
          if (c.ob) { status('Converting with Open Babel…'); out = await BM.native.obabel.convert(out, 'sdf', c.code); }
          if (typeof out === 'string') await saveFile(name, out);
          else if (out && out.files && out.files.length) for (let k = 0; k < out.files.length; k++) await saveFile(k === 0 ? name : out.files[k].name, out.files[k].text);
          else throw new Error('The writer produced no output.');
          try { localStorage.setItem('bluemolecule.saveas', JSON.stringify({ fmt: c.id, all: c.all })); } catch (e) { /* ignore */ }
        } catch (e) { toast(e.message || String(e), true); return true; }
        return false;
      }
    }], form => {
      $('[data-sa=fmt]', form).addEventListener('change', () => update(form));
      $$('input[name=sa-fr]', form).forEach(r => r.addEventListener('change', () => update(form)));
      $('[data-sa=name]', form).addEventListener('input', () => { nameEdited = true; });
      update(form);
      const ob = BM.native && BM.native.obabel, obInfo = $('[data-sa=ob]', form);
      const none = () => { obInfo.textContent = 'not available here (needs the app or launcher, and obabel)'; };
      if (!ob || typeof ob.available !== 'function') { none(); return; }
      Promise.resolve(ob.available()).then(okay => {
        if (!okay) return none();
        return ob.formats().then(list => {
          const w = list.filter(x => x.write && /^[A-Za-z0-9_-]{1,16}$/.test(x.code));
          const sel = $('[data-sa=fmt]', form);
          if (!sel || !document.body.contains(sel)) return;
          const og = document.createElement('optgroup');
          og.label = 'Open Babel (' + w.length + ')';
          og.innerHTML = w.map(x => `<option value="ob:${esc(x.code)}">${esc(x.code)} — ${esc(x.description)}</option>`).join('');
          sel.appendChild(og);
          obInfo.textContent = w.length + ' more formats at the end of the list';
          if (/^ob:/.test(prefs.fmt || '') && w.some(x => 'ob:' + x.code === prefs.fmt)) { sel.value = prefs.fmt; update(form); }
        });
      }).catch(none);
    });
  }

  function showHelp() {
    const k = (a, b) => `<div><kbd>${a}</kbd></div><div>${b}</div>`;
    const h = t => `<div class="keys-h">${t}</div>`;
    dialog('BlueMolecule ' + VERSION, `
      <p class="note" style="margin-top:0"><b>A vibe-coded project:</b> BlueMolecule's code was written by an AI assistant (Claude, by Anthropic) from its author's instructions and tested mostly on synthetic files and public example outputs. Check numbers that matter against the original output files. Documentation, source and bug reports: <a href="https://github.com/fred4h/BlueMolecule" target="_blank" rel="noopener">github.com/fred4h/BlueMolecule</a> · MIT License.</p>
      ${h('View')}
      <div class="keys cols2">
        ${k('Left drag', 'Rotate')}${k('Right / Shift+drag', 'Move (pan)')}${k('Alt/Option+drag', 'Rotate in the screen plane')}${k('Wheel / pinch', 'Zoom')}
        ${k('Click atom', 'Select; 2–4 atoms measure')}${k('Double-click atom', 'Centre view on atom')}${k('Click empty space', 'Clear selection')}${k('Ctrl/⌘ + O', 'Open files')}
        ${k('Ctrl/⌘ + Shift + S', 'Save as… (any format)')}${k('?', 'This help')}
        ${k('Space', 'Play / pause')}${k('← / →', 'Previous / next geometry')}${k('↑ / ↓', 'Previous / next row of the Steps, Vibrations or States list')}${k('Home / End', 'First / last geometry (or row)')}${k('1 2 3 4', 'Ball & stick, sticks, spacefill, wire')}
        ${k('H', 'Toggle hydrogens')}${k('L', 'Cycle labels')}${k('B', 'Toggle bond lengths')}${k('P', 'Perspective / orthographic')}
        ${k('R', 'Reset view')}${k('Esc', 'Stop animation or clean-up, clear selection')}
      </div>
      ${h('Build tools (palette at the left of the 3D view)')}
      <div class="keys">
        ${k('N', 'Navigate: rotate, pan, zoom; click atoms to measure')}
        ${k('S', 'Select: drag a box (Shift adds, Ctrl/⌘ removes), click toggles, double-click selects a fragment')}
        ${k('D', 'Draw: click to add an atom, drag from an atom for a bonded atom, click a bond to cycle 1→2→3, right-click or Alt/Option-click deletes an atom or bond')}
        ${k('T', 'Template: insert molecules and metal centres, attach groups, replace ligands')}
        ${k('K', 'Bond-centric: drag a bond to rotate its smaller side, Shift-drag to stretch it')}
        ${k('M', 'Manipulate: drag an atom to move its fragment or the selection, Shift-drag to rotate')}
        ${k('Right-click atom', 'Delete it (Draw, Template, Bond-centric and Manipulate; right-drag still pans)')}
        ${k('Ctrl/⌘ + A', 'Select all atoms')}${k('Delete / Backspace', 'Delete the selected atoms')}
        ${k('Ctrl/⌘ + Z', 'Undo (Ctrl/⌘ + Shift + Z or Ctrl + Y: redo)')}
      </div>
      ${h('Contents and plots')}
      <div class="keys">
        ${k('Info → Contents', 'Everything in the file with counts; click an entry to open its tab or plot. Several jobs: pick one to restrict steps and plots')}
        ${k('Plot picker', 'Optimisation: energy, ΔE vs final, QM/MM or ONIOM energy components, each criterion with its threshold (log scale), RMSD to the final geometry, largest atomic step, all criteria / threshold')}
        ${k('', 'Scan / IRC energy profiles · IR, Raman, UV-Vis, ECD (Δε from rotatory strengths), NMR (shieldings, or shifts δ = σref − σ) · MO level diagram · atomic charges')}
        ${k('Click the plot', 'Select a geometry, mode, state, level, atom or NMR signal; CSV, PNG and SVG export the current plot')}
        ${k('Drag on the plot', 'Zoom into the box (a flat box: x only); the magnifiers zoom in, out and fit; double-click shows the whole plot. Exports follow the zoom')}
      </div>
      ${h('Files')}
      <div class="keys">
        ${k('Open', FORMAT_SUMMARY)}
        ${k('MD files', 'Open a topology (prmtop, PSF) or structure first, then coordinates or trajectories (inpcrd, Amber NetCDF, DCD, TRR, XTC) with the same atom count')}
        ${k('Files of 1 GB+', 'A dialog asks what to load: every Nth step or frame, a range or the last N, and which vibrations; Info → Load other parts… reads it again')}
        ${k('Other formats', 'Converted with Open Babel when the BlueMolecule app or launcher finds obabel')}
        ${k('Download', 'A molecule from PubChem by name or CID, or an entry or biological assembly from the PDB by its ID; off until you allow it in the dialog')}
        ${k('Export', 'Save as… (structure, crystal, MD and QC formats, plus Open Babel), calculation inputs for ten programs, images, vibration animations')}
      </div>
      <p class="note">The Build tab adds SMILES, a fragment library, hydrogens, bond orders, alignment, a clean-up force field (not a QM method), symmetry and a coordinate editor. BlueMolecule is an independent viewer for quantum-chemistry work; everything runs locally and files never leave your computer.</p>`,
    (BM.native && BM.native.updates ? [{ label: 'Check for updates…', onClick: () => { checkUpdates(true); return true; } }] : []).concat([{ label: 'Close', primary: true }]));
  }

  // ---------------- updates (the desktop app, issue #12) ----------------
  const UPDATES_KEY = 'bluemolecule.updates.auto';
  function updatesAuto() { try { return localStorage.getItem(UPDATES_KEY) === '1'; } catch (e) { return false; } }
  /** Release notes are the release's Markdown: shown as plain text, without its emphasis, code and link marks. */
  const plainNotes = s => String(s).replace(/\r\n?/g, '\n').replace(/^#{1,6}\s+/gm, '').replace(/\*\*([^*\n]+)\*\*/g, '$1')
    .replace(/`([^`\n]+)`/g, '$1').replace(/\[([^\]\n]+)\]\((?:https?:[^)\s]+)\)/g, '$1').replace(/^(\s*)[-*]\s+/gm, '$1• ').trim();
  function setUpdatesAuto(on) { try { localStorage.setItem(UPDATES_KEY, on ? '1' : '0'); } catch (e) { /* ignore */ } }
  /** Help → Check for updates (manual), or the opt-in check at start-up (quiet unless there is a new version). The
   *  check asks GitHub for the latest release; the update is installed only after its signature is checked. */
  async function checkUpdates(manual) {
    const up = BM.native && BM.native.updates;
    if (!up) return;
    if (manual) status('Checking for updates…');
    let r;
    try { r = await up.check(); } catch (e) { r = { supported: up.supported, available: false, error: String((e && e.message) || e) }; }
    r = r || {};
    if (!manual) {
      if (r.available) toast(`BlueMolecule ${r.version} is available: Help → Check for updates.`);
      return;
    }
    status(r.available ? 'An update is available' : 'Ready');
    const box = r.supported === false ? '' : `<label class="upd-auto"><input type="checkbox" data-upd="auto"${updatesAuto() ? ' checked' : ''}> Check for updates when BlueMolecule starts</label>
      <p class="note">A check asks GitHub for the latest release and sends nothing else. It never happens unless you ask or tick this box. (BlueMolecule's only other use of the network is Download, from PubChem and the PDB, which is off until you allow it.)</p>`;
    const keep = form => { const c = form.querySelector('[data-upd="auto"]'); if (c) setUpdatesAuto(c.checked); };
    if (r.available) {
      dialog(`BlueMolecule ${esc(r.version)} is available`, `<p>You have version ${esc(r.current)}${r.date ? '; this one was released on ' + esc(String(r.date).slice(0, 10)) : ''}.</p>
        ${r.notes ? `<div class="upd-notes">${esc(plainNotes(r.notes))}</div>` : ''}${up.selfInstall === false ? '<p class="note">This is a .deb installation: install the new .deb package from the Releases page to update.</p>' : ''}${box}`, [
        { label: 'Later', onClick: f => { keep(f); } },
        up.selfInstall === false
          ? { label: 'Open the Releases page', primary: true, onClick: f => { keep(f); up.openReleases(); } }
          : { label: 'Install and restart', primary: true, onClick: async f => {
            keep(f);
            status('Downloading and installing the update…');
            try { await up.install(); } catch (e) { toast('The update was not installed: ' + ((e && e.message) || e), true); }
          } }]);
    } else {
      // (a build without the update key says so as a plain note: nothing went wrong)
      const msg = r.supported === false && r.error ? esc(r.error)
        : r.error ? 'The check did not succeed: ' + esc(r.error) : `BlueMolecule ${esc(r.current)} is the latest version.`;
      dialog('Updates', `<p>${msg}</p>${box}`, [{ label: 'Close', primary: true, onClick: f => { keep(f); } }]);
    }
  }

  // ---------------- downloads from PubChem and the PDB (issue #22) ----------------
  // Off until the user ticks "Allow downloads…" in the dialog, and even then a request goes out only when they press
  // Download. What is sent is what they typed, as part of the address (BM.fetchdb).
  const DOWNLOADS_KEY = 'bluemolecule.downloads.allow', DOWNLOADS_PREFS = 'bluemolecule.downloads';
  function downloadsAllowed() { try { return localStorage.getItem(DOWNLOADS_KEY) === '1'; } catch (e) { return false; } }
  function setDownloadsAllowed(on) { try { localStorage.setItem(DOWNLOADS_KEY, on ? '1' : '0'); } catch (e) { /* ignore */ } }
  function downloadPrefs() {
    let p = {};
    try { p = JSON.parse(localStorage.getItem(DOWNLOADS_PREFS) || '{}') || {}; } catch (e) { p = {}; }
    return { source: p.source === 'pdb' ? 'pdb' : 'pubchem', assembly: !!p.assembly };
  }
  const sizeText = n => n >= 1048576 ? (n / 1048576).toFixed(1) + ' MB' : Math.max(1, Math.round(n / 1024)) + ' kB';
  let download = null;   // the running download: {ctl: AbortController}
  function stopDownload() { if (download) { const j = download; download = null; j.ctl.abort(); } }
  /** The Download dialog: PubChem by name or CID, or the PDB by ID (optionally the first biological assembly). */
  function downloadDialog(source) {
    const p = downloadPrefs();
    if (source === 'pdb' || source === 'pubchem') p.source = source;
    let busy = false;
    const body = `<div class="fdb">
        <div class="field"><span>Database</span><span class="sa-seg">
          <label><input type="radio" name="fdb-src" value="pubchem" ${p.source === 'pubchem' ? 'checked' : ''}> PubChem (molecules)</label>
          <label><input type="radio" name="fdb-src" value="pdb" ${p.source === 'pdb' ? 'checked' : ''}> Protein Data Bank</label></span></div>
        <label class="field"><span data-fdb="ql">Name or CID</span><input type="text" data-fdb="q" spellcheck="false" autocomplete="off"></label>
        <label class="field" data-fdb="asmrow"><span></span><span class="inline"><input type="checkbox" data-fdb="asm" ${p.assembly ? 'checked' : ''}> biological assembly 1, instead of the deposited coordinates</span></label>
        <label class="fdb-allow"><input type="checkbox" data-fdb="allow" ${downloadsAllowed() ? 'checked' : ''}> Allow downloads from PubChem and the Protein Data Bank</label>
        <p class="note">When you press Download, BlueMolecule sends what you typed above, and nothing else, to pubchem.ncbi.nlm.nih.gov (PubChem) or files.rcsb.org (RCSB PDB). Your files never leave your computer. Untick the box to turn downloads off again.</p>
        <p class="note" data-fdb="msg" aria-live="polite"></p>
      </div>`;
    const src = form => { const r = $('input[name=fdb-src]:checked', form); return r ? r.value : 'pubchem'; };
    const say = (form, text, err) => { const m = $('[data-fdb=msg]', form); if (m) { m.textContent = text || ''; m.classList.toggle('err', !!err); } };
    const sync = form => {
      const s = src(form), q = $('[data-fdb=q]', form);
      $('[data-fdb=ql]', form).textContent = s === 'pdb' ? 'PDB ID' : 'Name or CID';
      q.placeholder = s === 'pdb' ? 'e.g. 1CRN or 4HHB' : 'e.g. caffeine, ibuprofen or 2244';
      $('[data-fdb=asmrow]', form).classList.toggle('hidden', s !== 'pdb');
      $$('input', form).forEach(i => { if (i.dataset.fdb !== 'allow') i.disabled = busy; });
      const go = $('[data-dbtn="1"]', form);
      if (go) go.disabled = busy || !$('[data-fdb=allow]', form).checked;
    };
    const unticked = 'Tick “Allow downloads…” first.';
    dialog('Download a structure', body, [
      { label: 'Cancel', onClick: () => { stopDownload(); } },
      { label: 'Download', primary: true, onClick: async form => {
        if (busy) return true;
        const s = src(form), q = $('[data-fdb=q]', form).value, asm = $('[data-fdb=asm]', form).checked;
        if (!$('[data-fdb=allow]', form).checked) { say(form, unticked, true); return true; }
        const chk = BM.fetchdb.parse(s, q);
        if (chk.error) { say(form, chk.error, true); $('[data-fdb=q]', form).focus(); return true; }
        try { localStorage.setItem(DOWNLOADS_PREFS, JSON.stringify({ source: s, assembly: asm })); } catch (e) { /* ignore */ }
        const where = s === 'pdb' ? 'the Protein Data Bank' : 'PubChem';
        stopDownload();
        const job = { ctl: new AbortController() };
        download = job; busy = true; sync(form);
        say(form, 'Downloading from ' + where + '…'); status('Downloading from ' + where + '…');
        const onProgress = n => { if (download === job) say(form, 'Downloading from ' + where + '… ' + sizeText(n)); };
        try {
          const opts = { signal: job.ctl.signal, onProgress, assembly: asm };
          const r = s === 'pdb' ? await BM.fetchdb.pdb(q, opts) : await BM.fetchdb.pubchem(q, opts);
          if (download !== job) return true;             // cancelled meanwhile (Cancel and Esc close the dialog themselves)
          download = null; busy = false;
          status('Ready');
          openDownload(r);
          return false;
        } catch (e) {
          if (download !== job) return true;
          download = null; busy = false; sync(form);
          status('Ready');
          say(form, e.message || String(e), true);
          return true;
        }
      } }], form => {
      sync(form);
      $$('input[name=fdb-src]', form).forEach(r => r.addEventListener('change', () => { say(form, ''); sync(form); $('[data-fdb=q]', form).focus(); }));
      $('[data-fdb=allow]', form).addEventListener('change', e => {
        setDownloadsAllowed(e.target.checked);
        if (!e.target.checked && busy) { stopDownload(); busy = false; status('Ready'); say(form, 'Download stopped.'); }
        else say(form, '');
        sync(form);
      });
      const qi = $('[data-fdb=q]', form);
      qi.addEventListener('keydown', e => {
        if (e.key !== 'Enter') return;
        e.preventDefault();                                // (Enter would otherwise submit, and close, the dialog)
        const go = $('[data-dbtn="1"]', form);
        if (go && !go.disabled) go.click(); else if (!busy) say(form, unticked, true);
      });
      qi.focus();
      $('#dlg').addEventListener('close', stopDownload, { once: true });
    });
  }
  /** A downloaded structure opens like a file of its format; a PubChem compound without a 3D conformer is built from
   *  its SMILES with BlueMolecule's own builder, and says so. */
  function openDownload(r) {
    const from = r.source === 'pdb' ? 'RCSB Protein Data Bank' + (r.assembly ? ', biological assembly 1' : '') : 'PubChem' + (r.cid ? ' CID ' + r.cid : '');
    let ds;
    if (r.text) {
      try { ds = BM.parseFile(r.fileName, r.text); } catch (e) { toast(conciseError(e.message), true); return null; }
      if (r.source === 'pubchem') ds.title = r.name + (r.cid && r.name !== 'CID ' + r.cid ? ' · PubChem CID ' + r.cid : '');
    } else {
      let b;
      try { b = BM.chem.smilesTo3D(r.smiles, { seed: 1 }); }
      catch (e) { toast('PubChem has no 3D structure of ' + r.name + ', and its SMILES could not be built: ' + e.message, true); return null; }
      const mol = b.mol, arom = new Set(b.aromaticBonds || []);
      ds = BM.emptyDataset('built', r.fileName);
      ds.program = 'BlueMolecule, from the PubChem SMILES';
      ds.title = r.name + (r.cid && r.name !== 'CID ' + r.cid ? ' · PubChem CID ' + r.cid : '');
      ds.bonds = mol.bonds.map((bd, k) => [bd[0], bd[1], arom.has(k) ? 1.5 : bd[2]]);
      const fr = { Z: Array.from(mol.Z), xyz: Float64Array.from(mol.xyz), energy: null, energies: null, label: 'Built from the PubChem SMILES', converged: false, opt: null, atomNames: null, resNames: null, resIds: null };
      if (mol.formalCharges && mol.formalCharges.some(q => q)) fr.formalCharges = Array.from(mol.formalCharges, q => q || 0);
      ds.frames = [fr];
      ds.warnings.push('PubChem has no 3D conformer of this compound, so BlueMolecule built one from its SMILES, ' + r.smiles +
        ', with its own builder and clean-up force field (not a PubChem conformer).' + (b.warnings.length ? ' ' + b.warnings.join(' ') : ''));
    }
    ds.downloaded = { from, url: r.url, date: new Date().toISOString().slice(0, 10) };
    const d = openData(r.fileName, ds, ' from ' + (r.source === 'pdb' ? 'the PDB' : 'PubChem'));
    if (d && !r.text && state.active === d) { d.tab = 'info'; renderTabs(); renderPanel(); }
    return d;
  }

  // ---------------- status ----------------
  function updateStatus() {
    const d = state.active;
    if (!d) { $('#status-info').textContent = ''; return; }
    const f = curFrame(d);
    $('#status-info').textContent = `${f.Z.length} atoms · ${geom.formula(f.Z)}${isNum(f.energy) ? ' · E = ' + f.energy.toFixed(6) + ' Eh' : ''}`;
  }
  function onHover(i) {
    const d = state.active;
    if (!d || i < 0) { $('#status-hover').textContent = ''; return; }
    const f = curFrame(d), x = R.xyz, q = chargesFor(d);
    let desc = '';
    if (f.resNames || layersFor(d, f)) { try { desc = REPS.describe(repCtx(d), i); } catch (e) { desc = ''; } }
    $('#status-hover').textContent = `${el.symbol(f.Z[i])}${i + 1}${f.atomNames && f.atomNames[i] ? ' ' + f.atomNames[i] : ''}${desc ? ' · ' + desc : ''}  ${x[3 * i].toFixed(4)} ${x[3 * i + 1].toFixed(4)} ${x[3 * i + 2].toFixed(4)}${q ? '  q = ' + q[i].toFixed(4) : ''}`;
  }

  // ---------------- toolbar sync ----------------
  function syncToolbar() {
    const v = state.view;
    $$('#style-seg button').forEach(b => b.classList.toggle('on', b.dataset.style === v.style));
    $('#sel-labels').value = v.labels; $('#sel-color').value = v.color; $('#sel-bg').value = v.bg;
    $('#tg-h').classList.toggle('on', v.showH); $('#tg-bl').classList.toggle('on', v.bondLengths);
    $('#tg-persp').classList.toggle('on', v.perspective); $('#tg-fog').classList.toggle('on', v.fog); $('#tg-axes').classList.toggle('on', v.axes);
    document.documentElement.dataset.theme = v.theme;
    const setR = (id, val) => { const i = $(id); i.value = val; i.nextElementSibling.textContent = Number(val).toFixed(2); };
    setR('#set-atom', v.atomScale); setR('#set-bond', v.bondRadius); setR('#set-vdw', v.vdwScale); setR('#set-tol', v.tol);
    $('#set-dummy').checked = v.dummies;
    $('#set-orders').checked = v.bondOrders; $('#set-hbonds').checked = v.hbonds; $('#set-dipole').checked = v.dipole; $('#set-cell').checked = v.cell;
    $('#set-scene').value = SCENES[v.scene] ? v.scene : 'custom'; $('#set-material').value = v.material;
    const setF = (id, val, d) => { const i = $(id); i.value = val; i.nextElementSibling.textContent = +val ? Number(val).toFixed(d) : 'off'; };
    setF('#set-ao', v.ao, 2); setF('#set-outline', v.outline, 1); setF('#set-dof', v.dof, 2);
    $('#set-shadows').checked = v.shadows; $('#set-floor').checked = v.floor; $('#set-gradient').checked = v.gradient; $('#set-fxlarge').checked = v.fxLarge;
    syncFxNote();
  }
  /** Says when shadows, ambient occlusion or depth of field wait because the structure is very large. */
  function syncFxNote() {
    const v = state.view, d = state.active, n = d && !d.pending ? curFrame(d).Z.length : 0, note = $('#set-fxnote');
    const wait = n > BM.Renderer.LARGE_SYSTEM && !v.fxLarge && (v.shadows || v.ao > 0 || v.dof > 0);
    note.textContent = wait ? `Shadows, ambient occlusion and depth of field are paused for this structure (${n.toLocaleString('en-US')} atoms, over ${BM.Renderer.LARGE_SYSTEM.toLocaleString('en-US')}).` : '';
    note.classList.toggle('hidden', !wait);
  }
  /** A scene preset: its background, material and effects. */
  function applyScene(name) {
    const s = SCENES[name];
    if (!s) return;
    setView(Object.assign({ scene: name }, s));
  }
  function setView(patch, rebuildBonds) {
    Object.assign(state.view, patch);
    syncToolbar();
    const d = state.active;
    if (d && d.reps && !d.pending && ('style' in patch || 'color' in patch)) {
      if (d.reps.length === 1 && d.reps[0].base) { d.reps[0].style = state.view.style; d.reps[0].color = state.view.color === 'charge' ? 'charge' : 'element'; }
      else { d.reps = [baseRow()]; d.repPreset = 'all'; d.repSel = 0; status('Display: all atoms in the toolbar style (Display ▸ Presets brings the preset back)'); }
      if (!rebuildBonds) { R.setOptions(renderOptions()); refreshReps(d); if (d.tab === 'display') renderPanel(); savePrefs(); return; }
    }
    if (rebuildBonds && state.active) refreshScene(true);
    else {
      R.setOptions(renderOptions()); applyColors(); $('#viewport').style.background = 'rgb(' + (BG[state.view.bg] || BG.navy).map(c => Math.round(c * 255)).join(',') + ')';
      if ('hbonds' in patch || 'dipole' in patch || 'cell' in patch || 'showH' in patch) refreshExtras();
    }
    savePrefs();
  }

  // ---------------- events ----------------
  function bindEvents() {
    $('#btn-open').onclick = () => openDialog();
    $('#btn-open-2').onclick = () => openDialog();
    $('#btn-open-3').onclick = () => openDialog();
    $('#btn-fetch').onclick = () => downloadDialog();
    $('#btn-fetch-2').onclick = () => downloadDialog();
    $('#file-input').onchange = e => { openFiles(e.target.files); e.target.value = ''; };
    $('#file-input').addEventListener('cancel', () => { layersRequest = null; });
    $('#btn-image').onclick = exportImage;
    $('#btn-help').onclick = showHelp;
    $('#btn-theme').onclick = () => { setView({ theme: state.view.theme === 'dark' ? 'light' : 'dark' }); renderPlotDraw(); };
    $$('#style-seg button').forEach(b => b.onclick = () => setView({ style: b.dataset.style }));
    $('#sel-labels').onchange = e => { setView({ labels: e.target.value }); if (e.target.value === 'charge' && !chargesFor(state.active)) toast('This file has no atomic charges for the current geometry.'); };
    $('#sel-color').onchange = e => { setView({ color: e.target.value }); if (e.target.value === 'charge' && !chargesFor(state.active)) toast('This file has no atomic charges for the current geometry.'); };
    // the background and depth cueing are part of a scene: a change that leaves the preset's value keeps the scene
    const sceneIf = (key, val) => SCENES[state.view.scene] && SCENES[state.view.scene][key] === val ? state.view.scene : 'custom';
    $('#sel-bg').onchange = e => setView({ bg: e.target.value, scene: sceneIf('bg', e.target.value) });
    $('#tg-h').onclick = () => setView({ showH: !state.view.showH });
    $('#tg-bl').onclick = () => setView({ bondLengths: !state.view.bondLengths });
    $('#tg-persp').onclick = () => setView({ perspective: !state.view.perspective });
    $('#tg-fog').onclick = () => setView({ fog: !state.view.fog, scene: sceneIf('fog', !state.view.fog) });
    $('#tg-axes').onclick = () => setView({ axes: !state.view.axes });
    const rng = (id, key, rebuild) => { $(id).oninput = e => { e.target.nextElementSibling.textContent = Number(e.target.value).toFixed(2); setView({ [key]: parseFloat(e.target.value) }, rebuild); }; };
    rng('#set-atom', 'atomScale'); rng('#set-bond', 'bondRadius'); rng('#set-vdw', 'vdwScale'); rng('#set-tol', 'tol', true);
    $('#set-dummy').onchange = e => setView({ dummies: e.target.checked });
    $('#set-orders').onchange = e => setView({ bondOrders: e.target.checked }, true);
    $('#set-hbonds').onchange = e => setView({ hbonds: e.target.checked });
    $('#set-dipole').onchange = e => { setView({ dipole: e.target.checked }); const d = state.active; if (e.target.checked && d && !d.data.dipole) toast('This file has no dipole moment.'); else if (e.target.checked && d && curFrame(d).edited) toast('The dipole arrow is shown for file geometries only.'); };
    $('#set-cell').onchange = e => { setView({ cell: e.target.checked }); const d = state.active; if (e.target.checked && d && !cellOf(d)) toast('This structure has no unit cell.'); };
    // scene: a preset, or each effect by hand (which makes the scene "custom")
    $('#set-scene').onchange = e => { if (SCENES[e.target.value]) applyScene(e.target.value); else setView({ scene: 'custom' }); };
    $('#set-material').onchange = e => setView({ material: e.target.value, scene: 'custom' });
    for (const [id, key] of [['#set-ao', 'ao'], ['#set-outline', 'outline'], ['#set-dof', 'dof']]) $(id).oninput = e => setView({ [key]: parseFloat(e.target.value), scene: 'custom' });
    for (const [id, key] of [['#set-shadows', 'shadows'], ['#set-floor', 'floor'], ['#set-gradient', 'gradient']]) $(id).onchange = e => setView({ [key]: e.target.checked, scene: 'custom' });
    $('#set-fxlarge').onchange = e => setView({ fxLarge: e.target.checked });
    $('#btn-new').onclick = () => ui.newMolecule();
    $('#btn-new-2').onclick = () => ui.newMolecule();
    const emptySmiles = () => { const t = $('#empty-smiles').value.trim(); if (!t) { toast('Type a SMILES string, e.g. CC(=O)Oc1ccccc1C(=O)O', true); return; } if (ui.insertSmiles(t)) $('#empty-smiles').value = ''; };
    $('#btn-smiles-2').onclick = emptySmiles;
    $('#empty-smiles').addEventListener('keydown', e => { if (e.key === 'Enter') { e.preventDefault(); emptySmiles(); } });

    // dropdowns
    document.addEventListener('click', e => {
      const t = e.target.closest('[data-dd]');
      $$('.dropdown.open').forEach(dd => { if (!t || dd !== t.parentElement) { if (!e.target.closest('.dropdown.open .menu.settings')) dd.classList.remove('open'); } });
      if (t) t.parentElement.classList.toggle('open');
    });
    // actions
    document.addEventListener('click', e => {
      const a = e.target.closest('[data-act]');
      if (!a || a.disabled) return;
      const act = a.dataset.act;
      if (a.tagName === 'A') e.preventDefault();
      const dd = a.closest('.dropdown'); if (dd) dd.classList.remove('open');
      const needDoc = () => { if (!state.active) { toast('Open a file first.', true); return false; } return true; };
      switch (act) {
        case 'save-as': if (needDoc()) saveAs(); break;
        case 'export-calc': if (needDoc()) exportCalc(); break;
        case 'export-gjf': if (needDoc()) exportGaussian(); break;
        case 'export-orca': if (needDoc()) exportOrca(); break;
        case 'export-xyz': exportSimple('xyz'); break;
        case 'export-xyz-all': exportSimple('xyz-all'); break;
        case 'export-mol': exportSimple('mol'); break;
        case 'export-pdb': exportSimple('pdb'); break;
        case 'export-vib': exportVibration(); break;
        case 'export-image': exportImage(); break;
        case 'view-reset': resetView(); break;
        case 'view-x': R.viewAlong('x'); break;
        case 'view-y': R.viewAlong('y'); break;
        case 'view-z': R.viewAlong('z'); break;
        case 'view-center-sel': { const d = state.active; if (d && d.selection.length) R.centerOn(d.selection[d.selection.length - 1]); break; }
        case 'plot-csv': plotCSV(); break;
        case 'plot-zoom-in': plot.zoomIn(); break;
        case 'plot-zoom-out': plot.zoomOut(); break;
        case 'plot-zoom-fit': plot.fit(); break;
        case 'plot-png': if (state.active && state.active.plotType) { const nm = `${baseName(state.active)}_${state.active.plotType.replace(':', '_')}.png`; plot.toPNG(2).then(b => saveFile(nm, b, 'image/png')); } break;
        case 'plot-svg': if (state.active && state.active.plotType) saveFile(`${baseName(state.active)}_${state.active.plotType.replace(':', '_')}.svg`, plot.toSVG(960, 440), 'image/svg+xml'); break;
        case 'displace-plus': displace(1); break;
        case 'reload-part': reloadPart(state.active); break;
        case 'layers-open': openLayersFile(state.active); break;
        case 'displace-minus': displace(-1); break;
        case 'vib-inputs': vibInputs(); break;
        case 'vib-imag': { const d = state.active, v = d && d.data.vibrations, st = v ? BM.vibtools.stationary(v) : null; if (st && st.imag.length) { selectVib(d, st.imag[0]); renderPanel(); } break; }
        case 'edit-set': editSet(); break;
        case 'edit-element': editElement(); break;
        case 'edit-delete': editDelete(); break;
        case 'edit-center': editCenter(); break;
        case 'undo': undo(false); break;
        case 'redo': undo(true); break;
        case 'example': { const ex = (BM.examples || []).find(x => x.id === a.dataset.ex); if (ex) openText(ex.name, ex.text); break; }
        default: if (act.startsWith('rp-')) displayAction(act, a); else if (act.startsWith('b-')) ui.action(act, a, e); else if (act.startsWith('sf-')) surfAction(act, a);
      }
    });
    // documents
    $('#doc-list').addEventListener('click', e => {
      const x = e.target.closest('[data-close]');
      if (x) { e.stopPropagation(); closeDoc(state.docs.find(d => d.id === +x.dataset.close)); return; }
      const it = e.target.closest('[data-doc]');
      if (it) { const d = state.docs.find(q => q.id === +it.dataset.doc); if (d && d !== state.active) activate(d); }
    });
    $('#tabs').addEventListener('click', e => {
      const b = e.target.closest('[data-tab]'); if (!b || !state.active) return;
      state.active.tab = b.dataset.tab; renderTabs(); renderPanel();
    });
    // panel
    const panel = $('#panel');
    panel.addEventListener('click', e => {
      const d = state.active; if (!d) return;
      if (e.target.closest('[data-take-layers-go]')) {
        const sel = $('[data-take-layers]', panel);
        takeLayers(d, sel ? state.docs.find(o => o.id === +sel.value) : null);
        return;
      }
      if (displayClick(e)) return;
      const go = e.target.closest('[data-go-plot]');
      if (go && (go.dataset.goPreset || go.dataset.goSelect)) {
        e.preventDefault();
        if (go.dataset.goPreset && d.repPreset !== go.dataset.goPreset) applyPreset(d, go.dataset.goPreset);
        if (go.dataset.goSelect) selectBy(d, go.dataset.goSelect);
        return;
      }
      if (go) { e.preventDefault(); goTo(d, go.dataset.goTab, go.dataset.goPlot, go.dataset.goAnchor); return; }
      const row = e.target.closest('tr.row');
      if (!row) return;
      const t = row.closest('table');
      if (t.id === 'steps-table') setFrame(d, +row.dataset.i);
      else if (t.id === 'vib-table') { const i = +row.dataset.i; if (d.vibIndex === i && state.anim.kind === 'vib') toggleVibPlay(); else { selectVib(d, i); renderPanel(); } }
      else if (t.id === 'states-table') selectState(d, +row.dataset.i);
      else if (t.id === 'vib-bonds-table') {   // measure the pair (its distance follows the animation)
        d.selection = row.dataset.pair.split(',').map(Number);
        R.setSelection(d.selection); updateMeasure();
        for (const r of $$('tr.row', t)) r.classList.toggle('on', r === row);
      }
      else if (row.dataset.atom !== undefined) { onAtomClick(+row.dataset.atom); }
    });
    panel.addEventListener('input', e => {
      const d = state.active; if (!d) return;
      const t = e.target;
      if (ui.panelInput(e)) return;
      if (displayInput(e)) return;
      if (t.dataset.vib) {
        const k = t.dataset.vib;
        if (k === 'arrows') { state.vib.arrows = t.checked; refreshArrows(); renderPlaybar(); }
        else if (k === 'ends') { state.vib.ends = t.checked; if (t.checked) holdVib(d); refreshEnds(); }
        else {
          state.vib[k] = parseFloat(t.value);
          if (t.nextElementSibling) t.nextElementSibling.textContent = k === 'amplitude' ? state.vib.amplitude.toFixed(2) + ' Å' : k === 'arrowLen' ? state.vib.arrowLen.toFixed(1) + ' Å' : state.vib.speed.toFixed(1) + '×';
          if (k === 'amplitude') setAmplitude(d, state.vib.amplitude);
          if (k === 'arrowLen') refreshArrows();
          renderPlaybar();
        }
        savePrefs();
        return;
      }
      if (t.id === 'disp-amt') { const x = parseFloat(t.value); if (x > 0) state.vib.disp = x; return; }
      if (t.dataset.orb) {
        const o = d.orb, k = t.dataset.orb;
        if (!o) return;
        if (k === 'index') o.index = +t.value;
        else o[k] = t.value;
        if (k === 'kind' || k === 'spin') renderPanel();
        return;
      }
      if (t.dataset.ms) {
        const m = d.msurf, k = t.dataset.ms;
        if (!m) return;
        if (k === 'probe') { const v = parseFloat(t.value); if (v > 0 && v <= 5) m.probe = v; }
        else { m[k] = t.value; renderPanel(); }
        return;
      }
      if (t.dataset.s) {
        const card = t.closest('[data-surf]'), s = d.surfaces[+card.dataset.surf], k = t.dataset.s;
        if (k === 'esp') { s.esp.src = t.value; s._espKey = null; refreshSurfaces(); renderPanel(); return; }
        if (k === 'espRange') { const v = parseFloat(t.value); s.esp.range = v > 0 ? v : null; s._espKey = null; scheduleSurfaces(); return; }
        if (k === 'visible') s.visible = t.checked;
        else if (k === 'both') s.both = t.checked;
        else if (k === 'opacity') { s.opacity = parseFloat(t.value); t.nextElementSibling.textContent = Math.round(s.opacity * 100) + '%'; }
        else if (k === 'colorPos' || k === 'colorNeg') s[k] = t.value;
        else if (k === 'iso') { const v = Math.abs(parseFloat(t.value)); if (v > 0) s.iso = v; else return; }
        else if (k === 'isoSlider') {
          const mx = s.stats.absMax || 1, lo = Math.max(mx * 1e-4, 1e-6);
          s.iso = Number((lo * Math.pow(mx / lo, parseFloat(t.value))).toPrecision(3));
          $('[data-s=iso]', card).value = s.iso;
        }
        scheduleSurfaces();
        return;
      }
      if (t.dataset.energy) { d.spec.energy[t.dataset.energy] = t.value; renderPanel(); if (specKey(d.plotType || '') === 'energy') renderPlot(); return; }
      if (t.dataset.job !== undefined) { setJob(d, +t.value); return; }
      if (t.dataset.convOnly !== undefined) { d.convergedOnly = t.checked; renderPanel(); renderPlaybar(); renderPlotDraw(); return; }
      if (t.dataset.chargeScheme !== undefined) { d.chargeScheme = t.value; applyColors(); R.requestRender(); if (d.plotType === 'charges') renderPlot(); return; }
    });
    // plot
    $('#plot-pick').addEventListener('change', e => { const d = state.active; if (!d) return; d.plotType = e.target.value; renderPlot(); e.target.blur(); });
    $('#plot-controls').addEventListener('click', e => {
      const b = e.target.closest('[data-plot-go]'), d = state.active;
      if (b && d) { d.plotType = b.dataset.plotGo; renderPlot(); return; }
      if (e.target.closest('[data-show-mo]') && d && d.moSel !== null) showPlotOrbital(d, d.moSel);
    });
    $('#plot-controls').addEventListener('change', e => {
      const d = state.active, t = e.target; if (!d || !t.dataset.pc) return;
      const key = specKey(d.plotType), spec = d.spec[key], k = t.dataset.pc;
      if (k === 'scheme') { d.chargeScheme = t.value; applyColors(); R.requestRender(); if (d.tab === 'info') renderPanel(); renderPlot(); return; }
      if (t.type === 'checkbox') spec[k] = t.checked;
      else if (t.type === 'number') { const v = parseFloat(t.value); spec[k] = isNum(v) ? v : null; if ((k === 'fwhm' || k === 'sigma') && !(spec[k] > 0)) spec[k] = k === 'fwhm' ? 12 : 0.4; if (k === 'scale' && !(spec[k] > 0)) spec[k] = 1; if (k === 'ref' && !isNum(spec[k])) spec[k] = 0; if (k === 'window' && !(spec[k] > 0)) spec[k] = 10; }
      else spec[k] = k === 'z' ? +t.value : t.value;
      if (k === 'unit' && (d.plotType === 'uv' || d.plotType === 'ecd')) { spec.xmin = null; spec.xmax = null; }
      if (key === 'energy' && d.tab === 'steps') renderPanel();
      renderPlot();
    });
    plot.onPick = (id, item) => {
      const d = state.active; if (!d) return;
      const t = d.plotType;
      if (framePlot(t)) { if (isNum(id) && id >= 0 && id < d.data.frames.length) setFrame(d, id); }
      else if (t === 'uv' || t === 'ecd') { d.stateIndex = d.data.excitedStates.findIndex(s => s.index === id); if (d.tab !== 'states') d.tab = 'states'; renderTabs(); renderPanel(); renderPlotDraw(); }
      else if (t === 'nmr') { d.selection = (item && item.atoms ? item.atoms : [id]).slice(); R.setSelection(d.selection); updateMeasure(); renderPlotDraw(); status(item && item.label ? item.label : ''); }
      else if (t === 'charges') { onAtomClick(id); renderPlotDraw(); }
      else if (t === 'mo') { d.moSel = id; const r = $('#mo-read'); if (r) r.textContent = moLabel(d, id); const sb = $('#plot-controls [data-show-mo]'); if (sb) sb.disabled = false; status(moLabel(d, id)); renderPlotDraw(); }
      else { if (d.tab !== 'vib') { d.tab = 'vib'; renderTabs(); } selectVib(d, id); renderPanel(); }
    };
    $('#btn-plot-collapse').onclick = () => { state.plotCollapsed = !state.plotCollapsed; $('#center').classList.toggle('plot-collapsed', state.plotCollapsed); };
    // playbar
    const pb = $('#playbar');
    pb.addEventListener('click', e => {
      const b = e.target.closest('button[data-pb]'); if (!b) return;
      const k = b.dataset.pb;
      if (k === 'play') togglePlay();
      else if (k === 'prev') stepFrame(-1);
      else if (k === 'next') stepFrame(1);
      else if (k === 'stop') stopVib();
      else if (k === 'arrows') { state.vib.arrows = !state.vib.arrows; refreshArrows(); renderPlaybar(); if (state.active.tab === 'vib') renderPanel(); savePrefs(); }
    });
    pb.addEventListener('input', e => {
      const t = e.target, k = t.dataset.pb, d = state.active; if (!d) return;
      if (k === 'frame') { const vis = visibleFrames(d); setFrame(d, vis[+t.value]); }
      else if (k === 'fps') state.fps = +t.value;
      else if (k === 'amp') {
        state.vib.amplitude = +t.value;
        setAmplitude(d, state.vib.amplitude);
        const s = $('#panel [data-vib="amplitude"]');   // keep the Vibrations tab's slider in step
        if (s) { s.value = state.vib.amplitude; if (s.nextElementSibling) s.nextElementSibling.textContent = state.vib.amplitude.toFixed(2) + ' Å'; }
        savePrefs();
      }
      else if (k === 'speed') { state.vib.speed = +t.value; savePrefs(); }
    });
    panel.addEventListener('keydown', e => { if (displayKey(e)) return; ui.panelKey(e); });
    panel.addEventListener('focusout', e => { if (e.target.id === 'rp-q') setTimeout(() => { if (document.activeElement !== $('#rp-q')) hideAc(); }, 120); });
    panel.addEventListener('focusin', e => { if (e.target.id === 'rp-q' && state.active) showAc(state.active); });
    panel.addEventListener('mousedown', e => { if (e.target.closest('#rp-ac')) e.preventDefault(); });
    // renderer callbacks (the active build tool gets the first say)
    R.onAtomClick = (i, e) => onAtomClick(i, e);
    R.onBackgroundClick = e => { if (!ui.backgroundClick(e)) clearSelection(); };
    R.onAtomDblClick = (i, e) => { if (!ui.atomDblClick(i, e) && i >= 0) R.centerOn(i); };
    R.overlayHook = (ctx, project) => ui.overlay(ctx, project);
    R.onHover = onHover;
    R.labelProvider = labelProvider;

    // drag & drop
    let dragDepth = 0;
    const hasFiles = e => e.dataTransfer && Array.from(e.dataTransfer.types || []).includes('Files');
    window.addEventListener('dragenter', e => { if (!hasFiles(e)) return; e.preventDefault(); dragDepth++; $('#drop').classList.remove('hidden'); });
    window.addEventListener('dragover', e => { if (hasFiles(e)) e.preventDefault(); });
    window.addEventListener('dragleave', e => { if (!hasFiles(e)) return; dragDepth = Math.max(0, dragDepth - 1); if (!dragDepth) $('#drop').classList.add('hidden'); });
    window.addEventListener('drop', e => { if (!hasFiles(e)) return; e.preventDefault(); dragDepth = 0; $('#drop').classList.add('hidden'); openFiles(e.dataTransfer.files); });

    // keyboard
    window.addEventListener('keydown', e => {
      const tag = (e.target.tagName || '').toLowerCase();
      const mod = e.ctrlKey || e.metaKey;
      if (mod && e.key.toLowerCase() === 'o') { e.preventDefault(); openDialog(); return; }
      if (mod && e.shiftKey && e.key.toLowerCase() === 's') { e.preventDefault(); if (state.active && !$('#dlg').open) saveAs(); return; }
      if (tag === 'input' || tag === 'textarea' || tag === 'select' || $('#dlg').open) return;
      if (ui.keydown(e)) { e.preventDefault(); return; }
      const d = state.active;
      if (mod && e.key.toLowerCase() === 'z') { e.preventDefault(); undo(e.shiftKey); return; }
      if (mod && e.key.toLowerCase() === 'y') { e.preventDefault(); undo(true); return; }
      if (mod && e.key.toLowerCase() === 'a') { if (d) { e.preventDefault(); selectAll(); } return; }
      if (mod || e.altKey) return;
      switch (e.key) {
        case ' ': e.preventDefault(); togglePlay(); break;
        case 'ArrowLeft': if (d) { e.preventDefault(); stepFrame(-1); } break;
        case 'ArrowRight': if (d) { e.preventDefault(); stepFrame(1); } break;
        case 'ArrowUp': case 'ArrowDown': case 'PageUp': case 'PageDown': if (d && listKey(d, e.key)) e.preventDefault(); break;
        case 'Home': if (d && !listKey(d, 'Home')) setFrame(d, visibleFrames(d)[0]); break;
        case 'End': if (d && !listKey(d, 'End')) { const v = visibleFrames(d); setFrame(d, v[v.length - 1]); } break;
        case '1': setView({ style: 'ballstick' }); break;
        case '2': setView({ style: 'sticks' }); break;
        case '3': setView({ style: 'spacefill' }); break;
        case '4': setView({ style: 'wireframe' }); break;
        case 'h': case 'H': setView({ showH: !state.view.showH }); break;
        case 'b': case 'B': setView({ bondLengths: !state.view.bondLengths }); break;
        case 'p': case 'P': setView({ perspective: !state.view.perspective }); break;
        case 'l': case 'L': { const o = ['none', 'symbol', 'number', 'symnum', 'charge']; setView({ labels: o[(o.indexOf(state.view.labels) + 1) % o.length] }); break; }
        case 'r': case 'R': resetView(); break;
        case 'n': case 'N': ui.setTool('navigate'); break;
        case 's': case 'S': ui.setTool('select'); break;
        case 'd': case 'D': ui.setTool('draw'); break;
        case 't': case 'T': ui.setTool('template'); break;
        case 'k': case 'K': ui.setTool('bond'); break;
        case 'm': case 'M': ui.setTool('manipulate'); break;
        case 'Escape': if (d && d.vibIndex >= 0) stopVib(); else if (state.anim.kind) { stopAnim(false); renderPlaybar(); } else clearSelection(); break;
        case 'Delete': case 'Backspace': if (d && d.selection.length) { e.preventDefault(); editDelete(); } break;
        case '?': showHelp(); break;
      }
    });

    // resizers
    const drag = (rz, onMove) => {
      rz.addEventListener('pointerdown', e => {
        rz.setPointerCapture(e.pointerId); rz.classList.add('drag');
        const mv = ev => onMove(ev);
        const upf = () => { rz.classList.remove('drag'); rz.removeEventListener('pointermove', mv); rz.removeEventListener('pointerup', upf); savePrefs(); };
        rz.addEventListener('pointermove', mv); rz.addEventListener('pointerup', upf);
      });
    };
    drag($('#rz-side'), e => { const w = Math.max(240, Math.min(620, e.clientX)); document.documentElement.style.setProperty('--sidebar-w', w + 'px'); });
    drag($('#rz-plot'), e => { const r = $('#center').getBoundingClientRect(); const h = Math.max(120, Math.min(r.height - 160, r.bottom - e.clientY)); document.documentElement.style.setProperty('--plot-h', h + 'px'); });
  }
  let surfTimer = 0;
  function scheduleSurfaces() { clearTimeout(surfTimer); surfTimer = setTimeout(() => refreshSurfaces(), 40); }

  function renderExamples() {
    const ex = BM.examples || [];
    $('#examples').innerHTML = ex.length ? '<div class="lbl">Try an example</div>' + ex.map(x => `<button data-act="example" data-ex="${x.id}" title="${esc(x.note || '')}">${esc(x.label)}</button>`).join('') : '';
  }

  async function loadLaunchFiles() {
    const L = window.BM_LAUNCH;
    if (!L) return;
    if (L.ping) { const ping = () => fetch(L.ping, { method: 'POST', keepalive: true }).catch(() => {}); ping(); setInterval(ping, 4000); }
    for (const f of (L.files || []).slice().sort((a, b) => coordOnly(a) - coordOnly(b))) {   // topology before trajectory
      try {
        const r = await fetch(f.url);
        if (!r.ok) throw new Error('HTTP ' + r.status);
        // text files stream straight into lines (no whole-file Blob or string); binary formats need random access
        const len = +(r.headers.get('Content-Length') || 0) || null;
        if (!r.body || !BIGF || BM.isBinaryFormat(f.name.replace(/\.gz$/i, ''))) await openFile(new File([await r.blob()], f.name));
        else await openFile({ name: f.name, size: len, stream: r.body });
      } catch (e) { toast(`Could not load ${f.name}: ${e.message}`, true); }
    }
  }

  function init() {
    loadPrefs();
    try { R = new BM.Renderer($('#gl'), $('#overlay')); }
    catch (e) { $('#empty .empty-card p').textContent = e.message; console.error(e); return; }
    plot = new BM.Plot($('#plot'));
    plot.onZoom = onPlotZoom;
    ui = BM.buildUI.create({
      state, R: () => R, esc, fx, toast, status, setHint, curFrame, cellOf, topoFor, bondsFor, ordersFor, topoFrame, measurement,
      doEdit, beginContinuousEdit, endContinuousEdit, afterPositions, frameMol, commitMol, addBuiltDoc, openData,
      renderPanel, updateMeasure, savePrefs, editSet, editElement, copyText
    });
    syncToolbar();
    bindEvents();
    renderExamples();
    R.setOptions(renderOptions());
    $('#viewport').style.background = 'rgb(' + (BG[state.view.bg] || BG.navy).map(c => Math.round(c * 255)).join(',') + ')';
    renderPlot();
    ui.setTool('navigate', true);
    $('#formats').innerHTML = '<b>QC</b> Gaussian, ORCA, GAMESS, NWChem, Q-Chem, Psi4, Molpro, Dalton, xtb, MOPAC, CP2K, Turbomole, fchk, Molden, cube<br>' +
      '<b>Structures</b> XYZ, MOL/SDF, PDB, mmCIF, CIF, Mol2, CML, CJSON, VASP, QE, CASTEP, XSF, GRO, LAMMPS<br>' +
      '<b>MD</b> prmtop/PSF + inpcrd, DCD, TRR, XTC, Amber NetCDF · .gz · others via Open Babel';
    const mac = /Mac|iPhone|iPad/.test(navigator.platform || navigator.userAgent || '');
    $$('[data-act="save-as"] kbd').forEach(k => { k.textContent = mac ? '⇧⌘S' : 'Ctrl+Shift+S'; });
    if (BM.native && BM.native.updates && BM.native.updates.supported && updatesAuto()) setTimeout(() => checkUpdates(false), 8000);
    if (BM.native && typeof BM.native.onOpenFile === 'function') {
      BM.native.onOpenFile((name, data) => { try { openHostFile(hostItem(name, data)); } catch (e) { toast(e.message || String(e), true); } });
    }
    loadLaunchFiles();
    BM.app = {
      state, openText, openFile, openFiles, openDialog, loadLaunchFiles, openData, downloadDialog, openDownload, takeLayers, layerDonors, activate, selectVib, setFrame, setView, R: () => R, plot: () => plot, exportSimple, exportGaussian, exportOrca, exportCalc, plotCSV, displace, onAtomClick, editSet, undo,
      saveAs, saveFile, saveContext, plotTypes, plotGroups, renderPlot, goTo, setJob, convFor, version: VERSION,
      curFrame, cellOf, topoFor, bondsFor, ordersFor, topoFrame, selectAll, editDelete, newMolecule: () => ui.newMolecule(), setTool: id => ui.setTool(id), ui: () => ui, kit,
      applyPreset, refreshReps, repCtx, selectBy, renderPanel, renderTabs
    };
  }
  if (document.readyState === 'loading') document.addEventListener('DOMContentLoaded', init); else init();
})(window.BM = window.BM || {});
