/* BlueMolecule — condensed reading of huge program outputs (BM.condense, issue #6).
 *
 * The outputs of big QM/MM jobs reach several GB. Most of such a file is geometry tables, normal-mode rows and
 * tables that no reader uses (forces, MO coefficients, density matrices...). A JavaScript engine holds at most
 * about 4 GB of text, so these files cannot be kept as lines. A condenser receives the lines one at a time:
 *   - geometry tables go to typed arrays (Float64Array coordinates, shared element lists);
 *   - normal-mode rows go to one integer array per mode (1 byte per value for Gaussian's 2-decimal modes);
 *   - a run of rows starting with a digit that the reader would skip anyway becomes one short line ("0");
 *   - every other line is kept.
 * The condensed lines carry a .bmSide store and are read by BM.parseLines like the full line array: lines that
 * start with U+0001 stand for the compacted tables (see readGaussOrientation and the ORCA / XYZ readers in
 * parsers.js). Each condenser classifies lines exactly as its reader does; tests/test_condense.js checks that
 * every example output reads the same either way.
 *
 * What to load (opts.load, the "Load part of a large file" dialog, issue #7):
 *   steps: {mode: 'all' | 'every' | 'last' | 'range', n, from, to}   (1-based steps, counted over the whole file)
 *     A step that is not loaded keeps its energies and convergence data (ds.steps; the plots show every step),
 *     only its geometry is not stored. The first and last geometry of every job are always loaded. For XYZ
 *     trajectories, frames that are not loaded are left out entirely.
 *   modes: {mode: 'all' | 'first' | 'none', n}   displacements of the first n vibrations (frequencies and
 *     intensities of all vibrations are always read).
 *
 *   BM.condense.create(formatId, opts) -> condenser | null      formats: BM.condense.formats
 *   condenser.push(line) ... condenser.finish(lastLine) -> lines (with .bmSide, .bmPartial)
 *   condenser.bytes() -> bytes in typed arrays;  condenser.textChars(), condenser.lineCount() -> the kept text
 *   BM.condense.own(s) -> a copy of s that does not keep the large text it was cut from alive
 */
(function (BM) {
  'use strict';
  const P = BM._parse;
  if (!P) throw new Error('condense.js needs parsers.js');
  const { tok, isIntTok, isNumTok, num, int, floatsIn, modRedCoord, DEC10, FREQ_LABEL, STARS } = P;
  const PH = '\u0001';
  const DROPPED = '0';   // stands for a run of dropped rows starting with a digit (the readers skip such rows)

  /** V8 keeps a substring as a slice of the 32 MB piece of text it was split from: a kept line must be a copy. */
  const own = s => (s.length < 13 ? s : (' ' + s).slice(1));

  // ------------------------------------------------------------------ fast, exact number scanning
  // A plain decimal with d digits after the point is the integer M / 10^d, which the division rounds exactly as
  // Number() does (both are the double nearest to the same rational): so these values equal the readers' num().
  const SP = new Int32Array(64);   // token spans [start, end) of the line being scanned
  /** Whitespace-separated tokens of s into SP; the count, or -1 (more than max, or characters beyond ASCII). */
  function spans(s, max) {
    let n = 0, i = 0;
    const L = s.length;
    for (;;) {
      let c = s.charCodeAt(i);
      while (i < L && (c === 32 || (c >= 9 && c <= 13))) c = s.charCodeAt(++i);
      if (i >= L) return n;
      if (n >= max) return -1;
      const a = i;
      while (i < L && !(c === 32 || (c >= 9 && c <= 13))) { if (c > 127) return -1; c = s.charCodeAt(++i); }
      SP[2 * n] = a; SP[2 * n + 1] = i; n++;
    }
  }
  let DM = 0, DD = 0;   // mantissa and decimals of the last plain decimal
  /** s[a, b) as a plain decimal ([-+]digits[.digits], at most 15 digits): sets DM, DD; false otherwise. */
  function dec(s, a, b) {
    let i = a, c = s.charCodeAt(i), neg = false;
    if (c === 45 || c === 43) { neg = c === 45; i++; }
    let m = 0, d = 0, digits = 0, dot = false;
    for (; i < b; i++) {
      c = s.charCodeAt(i);
      if (c >= 48 && c <= 57) { m = m * 10 + (c - 48); digits++; if (dot) d++; } else if (c === 46 && !dot) dot = true; else return false;
    }
    if (!digits || digits > 15) return false;
    DM = neg ? -m : m; DD = d;
    return true;
  }
  /** s[a, b) as a plain integer token (isIntTok): its value, or NaN. */
  function intAt(s, a, b) {
    let i = a, c = s.charCodeAt(i), neg = false;
    if (c === 45 || c === 43) { neg = c === 45; i++; }
    if (i >= b) return NaN;
    let m = 0;
    for (; i < b; i++) { c = s.charCodeAt(i); if (c < 48 || c > 57) return NaN; m = m * 10 + (c - 48); }
    return neg ? -m : m;
  }
  const decValue = () => (DD ? DM / DEC10[DD] : DM);
  const ROWV = new Int32Array(64);   // values of the last plain row (rowFast)
  /**
   * The hot path of huge frequency sections: "atom AN v v v ..." with nv plain values of exactly d decimals (Gaussian
   * prints F7.2), parsed in one pass into ROWV. Returns the first integer (the atom), or -1 when the line is anything
   * else (the general path then reads it exactly as the reader does).
   */
  function rowFast(s, nv, d) {
    const L = s.length;
    let i = 0, c = 0, atom = 0, k = 0;
    while (i < L && (c = s.charCodeAt(i)) === 32) i++;
    if (i === L || c < 48 || c > 57) return -1;
    while (i < L && (c = s.charCodeAt(i)) >= 48 && c <= 57) { atom = atom * 10 + c - 48; i++; }
    if (i === L || c !== 32) return -1;
    while (i < L && (c = s.charCodeAt(i)) === 32) i++;
    if (i === L || c < 48 || c > 57) return -1;
    while (i < L && (c = s.charCodeAt(i)) >= 48 && c <= 57) i++;
    for (; k < nv; k++) {
      if (i === L || s.charCodeAt(i) !== 32) return -1;
      while (i < L && (c = s.charCodeAt(i)) === 32) i++;
      let neg = false, m = 0, q = 0;
      if (c === 45) { neg = true; c = s.charCodeAt(++i); }
      if (c < 48 || c > 57) return -1;
      while (i < L && (c = s.charCodeAt(i)) >= 48 && c <= 57) { m = m * 10 + c - 48; i++; }
      if (c !== 46) return -1;
      i++;
      while (i < L && (c = s.charCodeAt(i)) >= 48 && c <= 57) { m = m * 10 + c - 48; i++; q++; }
      if (q !== d) return -1;
      ROWV[k] = neg ? -m : m;
    }
    while (i < L && s.charCodeAt(i) === 32) i++;
    return i === L ? atom : -1;
  }
  const DASHES10 = /^\s*-{10,}\s*$/;

  // ------------------------------------------------------------------ growable typed storage
  /** A vector of printed values kept as integers M with d decimals (Int8 -> Int16 -> Int32 as needed); a value
   *  that does not fit, or another number of decimals, turns it into a Float64Array of the values (d = -1). */
  function Vec(cap) { this.a = new Int8Array(Math.max(16, cap)); this.d = -2; this.len = 0; }
  Vec.prototype.room = function (i) {
    if (i < this.a.length) return;
    let n = this.a.length;
    while (n <= i) n *= 2;
    const b = new this.a.constructor(n); b.set(this.a); this.a = b;
  };
  Vec.prototype.toFloat = function () {
    const b = new Float64Array(this.a.length), s = this.d >= 0 ? DEC10[this.d] : 1;
    for (let k = 0; k < this.len; k++) b[k] = this.a[k] / s;
    this.a = b; this.d = -1;
  };
  /** value M / 10^d at index i */
  Vec.prototype.set = function (i, m, d) {
    if (this.d === -2) this.d = d;
    if (this.d === -1 || d !== this.d) { if (this.d !== -1) this.toFloat(); this.room(i); this.a[i] = d ? m / DEC10[d] : m; }
    else {
      const A = this.a;
      if (A instanceof Int8Array ? (m > 127 || m < -128) : A instanceof Int16Array ? (m > 32767 || m < -32768) : (m > 2147483647 || m < -2147483648)) {
        const T = A instanceof Int8Array && m <= 32767 && m >= -32768 ? Int16Array : A instanceof Int32Array || m > 2147483647 || m < -2147483648 ? null : Int32Array;
        if (T) { const b = new T(A.length); b.set(A); this.a = b; } else { this.toFloat(); this.a[i] = m / DEC10[d]; if (i >= this.len) this.len = i + 1; return; }
      }
      this.room(i);
      this.a[i] = m;
    }
    if (i >= this.len) this.len = i + 1;
  };
  /** a value that is not a plain decimal (NaN, as the readers' vibrations use, becomes 0) */
  Vec.prototype.setRaw = function (i, v) {
    if (this.d !== -1) { if (this.d === -2) this.d = 0; this.toFloat(); }
    this.room(i);
    this.a[i] = Number.isNaN(v) ? 0 : v;
    if (i >= this.len) this.len = i + 1;
  };
  Vec.prototype.done = function () { return { vals: this.a.length === this.len ? this.a : this.a.slice(0, this.len), d: this.d < 0 ? -1 : this.d }; };

  /** Geometry rows: element numbers (or label codes) and coordinates, growable. */
  function Rows(cap) { this.z = new Int32Array(Math.max(16, cap)); this.x = new Float64Array(3 * Math.max(16, cap)); this.n = 0; }
  Rows.prototype.add = function (z, x, y, w) {
    if (this.n === this.z.length) {
      const z2 = new Int32Array(this.z.length * 2); z2.set(this.z); this.z = z2;
      const x2 = new Float64Array(this.x.length * 2); x2.set(this.x); this.x = x2;
    }
    const k = this.n++;
    this.z[k] = z; this.x[3 * k] = x; this.x[3 * k + 1] = y; this.x[3 * k + 2] = w;
  };

  // ------------------------------------------------------------------ options, steps, the shared base
  function normOpts(opts) {
    const L = (opts && opts.load) || {};
    const s = L.steps || {}, m = L.modes || {};
    const pos = (v, d) => (Number.isFinite(+v) && +v >= 1 ? Math.floor(+v) : d);
    const steps = { mode: ['every', 'last', 'range'].indexOf(s.mode) >= 0 ? s.mode : 'all', n: pos(s.n, 1), from: pos(s.from, 1), to: pos(s.to, Infinity) };
    if (steps.mode === 'every' && steps.n <= 1) steps.mode = 'all';
    const modes = { mode: ['first', 'none'].indexOf(m.mode) >= 0 ? m.mode : 'all', n: Number.isFinite(+m.n) && +m.n >= 0 ? Math.floor(+m.n) : Infinity };
    return { steps, modes };
  }

  function Base(format, opts) {
    this.format = format;
    this.out = [];
    this.side = { format, geoms: [], modes: [], loaded: null };
    this.opts = normOpts(opts);
    this.ext = 0;          // bytes in typed arrays
    this.chars = 0;        // characters of the kept lines
    this.dropping = false;
    this.partial = false;
    // steps: a new step starts at a geometry after a step-closing line; its geometry is settled (kept or not)
    // when the next step starts, so that the last step of a job is known and always kept
    this.step = null; this.closed = true; this.jobFirst = true; this.nSteps = 0; this.steps = []; this.lastQ = [];
    this.zCache = [];
    this.modesSeen = 0; this.modesStored = 0;
  }
  Base.prototype.keep = function (line) { this.out.push(own(line)); this.chars += line.length; this.dropping = false; };
  Base.prototype.drop = function () { if (!this.dropping) { this.out.push(DROPPED); this.dropping = true; } };
  Base.prototype.bytes = function () { return this.ext; };
  Base.prototype.textChars = function () { return this.chars; };
  Base.prototype.lineCount = function () { return this.out.length; };
  /** One element list per distinct sequence (frames share it, as FrameSink does). */
  Base.prototype.internZ = function (z, n) {
    for (const c of this.zCache) {
      if (c.length !== n) continue;
      let same = true;
      for (let k = 0; k < n; k++) if (c[k] !== z[k]) { same = false; break; }
      if (same) return c;
    }
    const a = Array.from(z.subarray(0, n));
    this.zCache.unshift(a);
    if (this.zCache.length > 4) this.zCache.pop();
    return a;
  };
  /** Stores a geometry and puts its placeholder line in place of the rows. */
  Base.prototype.geometry = function (g) {
    const k = this.side.geoms.length;
    this.side.geoms.push(g);
    if (g.xyz) this.ext += g.xyz.byteLength;
    this.out.push(PH + 'G' + k);
    this.dropping = false;
    if (this.closed || !this.step) this.newStep();
    this.step.tables.push(k);
    return k;
  };
  Base.prototype.newStep = function () {
    if (this.step) this.settle(this.step, false);
    this.step = { n: this.nSteps++, tables: [], first: this.jobFirst, unloaded: false };
    this.steps.push(this.step);
    this.jobFirst = false; this.closed = false;
  };
  Base.prototype.closeStep = function () { this.closed = true; };
  Base.prototype.newJob = function () {
    if (this.step) this.settle(this.step, true);
    this.step = null; this.closed = true; this.jobFirst = true;
  };
  Base.prototype.settle = function (st, last) {
    const o = this.opts.steps;
    if (o.mode === 'all') return;
    if (o.mode === 'last') {
      // a window of the last N steps of the file, the last step of a job included; a job's last step stays loaded
      // after the window moves on. A job's first step is always loaded and takes no place in the window (it is
      // often the previous job's final geometry printed again, as in opt + freq).
      if (st.first) return;
      st.pinned = last;
      this.lastQ.push(st);
      while (this.lastQ.length > o.n) { const old = this.lastQ.shift(); if (!old.pinned) this.unload(old); }
      return;
    }
    if (last || st.first) return;
    const keep = o.mode === 'every' ? st.n % o.n === 0 : st.n + 1 >= o.from && st.n + 1 <= o.to;
    if (!keep) this.unload(st);
  };
  Base.prototype.unload = function (st) {
    for (const k of st.tables) { const g = this.side.geoms[k]; if (g.xyz) { this.ext -= g.xyz.byteLength; g.xyz = null; } }
    st.unloaded = true;
  };
  /** Mode k of the file (0-based, in the order printed) is stored? */
  Base.prototype.wantMode = function (rank) {
    const o = this.opts.modes;
    return o.mode === 'all' || (o.mode === 'first' && rank < o.n);
  };
  /** Ends the input: lastLine is the text after the last line break ('' for a complete file). */
  Base.prototype.finish = function (lastLine) {
    if (lastLine) this.partial = true;   // a job still running: its incomplete last line is not read
    this.close();
    if (this.step) this.settle(this.step, true);
    this.step = null;
    this.out.push('');
    const lines = this.out;
    const st = this.steps.length, kept = this.steps.filter(s => !s.unloaded).length;
    this.side.loaded = {
      steps: Object.assign({ total: st, loaded: kept }, this.opts.steps, { to: Number.isFinite(this.opts.steps.to) ? this.opts.steps.to : null }),
      modes: Object.assign({ total: this.modesSeen, loaded: this.modesStored }, this.opts.modes, { n: Number.isFinite(this.opts.modes.n) ? this.opts.modes.n : null }),
      bytes: this.ext
    };
    Object.defineProperty(lines, 'bmSide', { value: this.side, configurable: true });
    Object.defineProperty(lines, 'bmPartial', { value: this.partial, configurable: true });
    return lines;
  };
  Base.prototype.close = function () {};

  // ================================================================== Gaussian
  const M_PASS = 0, M_ROUTE = 1, M_OHEAD = 2, M_OROWS = 3, M_FREQ = 4, M_SHIELD = 5, M_MODRED = 6, M_TRI = 7, M_BERNY = 8;
  const G_EXC = /^Excited State\s+(\d+):\s+(\S+)\s+(\S+?)\s*eV\s*(\S+?)\s*nm\s+f=\s*(\S+)(?:\s+<S\*\*2>=\s*(\S+))?/;
  const G_TRANS = /^\d+[AB]?\s*(->|<-)\s*\d+[AB]?\s+\S+/;
  const G_SHIELD_ROW = /^\s*(\d+)\s+\S+\s+Isotropic\s*=\s*(\S+)\s+Anisotropy\s*=\s*(\S+)/, G_SHIELD_SUB = /^\s*(XX=|XY=|XZ=|Eigenvalues:)/;
  const G_BERNY_HEAD = /^Variable\s+Old X\s+-DE\/DX/, G_BERNY_ROW = /^\s*[A-Z][A-Za-z]?\d+\s+-?\d+\.\d+(?:\s+-?\d+\.\d+){4,6}\s*$/;
  const G_FREQ_TEXT = /^(activities|incident light|and normal coordinates|depolarization|reduced masses)/i;
  const G_LABEL = FREQ_LABEL;   // the reader's own rule (parsers.js)

  function Gaussian(opts) {
    Base.call(this, 'gaussian-log', opts);
    this.mode = M_PASS;
    this.prevDash = false; this.afterDash = -1;
    this.inState = false; this.statesSeen = false;
    this.ctx = null;       // a section whose digit rows a reader consumes: {digit: seen a digit row, idle: lines}
    this.n = 0;            // lines in the current special mode
    this.tri = false;
    this.rows = null; this.dashes = 0; this.lastN = 64;
    this.fr = null;        // frequency section state
    this.hpAt = -1;        // geometry count when the last high-precision mode section ended
    this.modeN = 0;        // atoms of the last mode block
  }
  Gaussian.prototype = Object.create(Base.prototype);

  Gaussian.prototype.push = function (line) {
    switch (this.mode) {
      case M_PASS: return this.pass(line);
      case M_ROUTE:
        this.keep(line);
        if (/^\s*-{3,}\s*$/.test(line)) { this.mode = M_PASS; this.prevDash = true; this.afterDash = 0; }
        return;
      case M_OHEAD:
        if (DASHES10.test(line) && ++this.dashes === 2) { this.mode = M_OROWS; this.rows = new Rows(this.lastN); this.tv = []; }
        return;   // the reader skips everything up to the second dashes line
      case M_OROWS: if (this.orientRow(line)) return; this.endTable(DASHES10.test(line)); return;   // the end line is consumed too
      case M_FREQ: return this.freq(line);
      case M_SHIELD:
        if (G_SHIELD_ROW.test(line) || G_SHIELD_SUB.test(line)) return this.keep(line);
        this.mode = M_PASS; return this.pass(line);
      case M_MODRED: {
        const tk = tok(line);
        if (this.n >= 199 || !tk.length || modRedCoord(tk) === undefined) { this.mode = M_PASS; return this.pass(line); }
        this.n++; return this.keep(line);
      }
      case M_TRI: {   // readLowerTriangle: header rows of column numbers, then "row value value ..."
        const tk = tok(line);
        if (tk.length && (tk.every(isIntTok) ? (this.tri = true) : this.tri && isIntTok(tk[0]) && tk.slice(1).every(isNumTok))) return this.keep(line);
        this.mode = M_PASS; return this.pass(line);
      }
      case M_BERNY:
        if (G_BERNY_ROW.test(line)) return this.drop();
        if (/^\s*\(Linear\)/.test(line)) return this.keep(line);
        this.mode = M_PASS; return this.pass(line);
    }
  };

  Gaussian.prototype.pass = function (line) {
    let p = 0, c0 = line.charCodeAt(0);
    while (c0 === 32) c0 = line.charCodeAt(++p);
    if (c0 >= 48 && c0 <= 57) {   // the reader skips such rows unless a section reader or an excited state takes them
      this.prevDash = false;
      if (this.ctx) { this.ctx.digit = true; return this.keep(line); }
      if (this.inState) {
        if (this.statesSeen && G_TRANS.test(line.slice(p))) return this.keep(line);
        this.inState = false;
      }
      if (this.afterDash >= 0 && this.afterDash < 5) { this.afterDash++; return this.keep(line); }   // title lines
      this.afterDash = -1;
      return this.drop();
    }
    const t = line.slice(p).trimStart();
    this.inState = false;
    if (this.ctx && (this.ctx.digit || ++this.ctx.idle > 20)) this.ctx = null;
    if (!t) { this.prevDash = false; this.afterDash = this.afterDash >= 0 && this.afterDash < 5 ? this.afterDash + 1 : -1; return this.keep(line); }
    const c = t[0];
    if (c === '-' && /^-{3,}\s*$/.test(t)) { this.keep(line); this.prevDash = true; this.afterDash = 0; return; }
    const wasDash = this.prevDash;
    this.prevDash = false;
    this.afterDash = this.afterDash >= 0 && this.afterDash < 5 ? this.afterDash + 1 : -1;
    if (c === '#' && wasDash) { this.newJob(); this.hpAt = -1; this.mode = M_ROUTE; return this.keep(line); }
    if ((c === 'S' || c === 'I' || c === 'Z') && (t.startsWith('Standard orientation:') || t.startsWith('Input orientation:') || t.startsWith('Z-Matrix orientation:'))) {
      this.keep(line);
      this.mode = M_OHEAD; this.dashes = 0;
      return;
    }
    if (c === 'H' && t.startsWith('Harmonic frequencies')) {
      this.keep(line);
      this.mode = M_FREQ;
      this.fr = { block: null, expectSym: false, freqs: false, hp: false, run: null, skip: this.hpAt === this.side.geoms.length, rank: 0, rowsIn: false, odd: 0, lastA: -1 };
      return;
    }
    if (c === 'E' && t.startsWith('Excited State')) {
      this.keep(line);
      if (G_EXC.test(t)) { this.inState = true; this.statesSeen = true; }
      return;
    }
    if ((c === 'S' || c === 'C' || c === 'I' || c === 'M' || c === 'G') && t.indexOf('Magnetic shielding tensor (ppm)') >= 0) { this.keep(line); this.mode = M_SHIELD; return; }
    if (c === 'T' && t.startsWith('Total nuclear spin-spin coupling J (Hz)')) { this.keep(line); this.mode = M_TRI; this.tri = false; return; }
    if (c === 'T' && t.startsWith('The following ModRedundant input section has been read')) { this.keep(line); this.mode = M_MODRED; this.n = 0; return; }
    if ((c === 'M' && t.startsWith('Mulliken')) || (c === 'E' && t.startsWith('ESP charges:')) || (c === 'C' && t.startsWith('Charges from ESP fit')) ||
        (c === 'A' && t.startsWith('APT ')) || (c === 'H' && (t.startsWith('Hirshfeld charges') || t.startsWith('Hirshfeld spin densities'))) ||
        (c === 'R' && t.startsWith('Rotatory Strengths (R) in cgs'))) {
      this.keep(line); this.ctx = { digit: false, idle: 0 }; return;
    }
    if (c === 'V' && G_BERNY_HEAD.test(t)) { this.keep(line); this.mode = M_BERNY; return; }
    if ((c === 'S' && t.startsWith('Step number')) || (c === 'P' && t.startsWith('Point Number:')) ||
        (c === 'I' && t.startsWith('Item') && /Value\s+Threshold\s+Converged/.test(t))) this.closeStep();
    // the structure printed again after convergence belongs to the converged step (the reader merges it there)
    if (c === 'O' && t.startsWith('Optimization completed') && this.step) this.closed = false;
    this.keep(line);
  };

  /** One row of an orientation table (readGaussOrientation), or false when the table has ended. */
  Gaussian.prototype.orientRow = function (line) {
    const n = spans(line, 8);
    let an, x, y, z;
    if ((n === 5 || n === 6) && !Number.isNaN(intAt(line, SP[0], SP[1])) && !Number.isNaN(an = intAt(line, SP[2], SP[3]))) {
      let ok = true;
      for (let k = 2; k < n - 3 && ok; k++) ok = dec(line, SP[2 * k], SP[2 * k + 1]);
      if (ok && dec(line, SP[2 * n - 6], SP[2 * n - 5])) { x = decValue(); if (dec(line, SP[2 * n - 4], SP[2 * n - 3])) { y = decValue(); if (dec(line, SP[2 * n - 2], SP[2 * n - 1])) z = decValue(); } }
      if (z === undefined) an = undefined;
    } else an = undefined;
    if (an === undefined) {   // anything unusual (glued fields, exponents): exactly as the reader does it
      const t = tok(line);
      if (t.length < 3 || !isIntTok(t[0]) || !isIntTok(t[1])) return false;
      const f = floatsIn(line);
      if (f.length < 5) return false;
      an = int(t[1]); x = f[f.length - 3]; y = f[f.length - 2]; z = f[f.length - 1];
    }
    if (an === -2) this.tv.push([x, y, z]);
    else if (an >= 0) this.rows.add(an, x, y, z);
    return true;
  };
  Gaussian.prototype.endTable = function (complete) {
    const r = this.rows || new Rows(16);
    this.rows = null;
    this.mode = M_PASS;
    this.lastN = Math.max(16, r.n);
    const Z = this.internZ(r.z, r.n);
    this.geometry({ Z, xyz: r.x.slice(0, 3 * r.n), tv: this.tv || [], complete });
    this.tv = null;
  };

  /** The frequency section (parseGaussFreqSection): the atom rows of each block go to per-mode vectors. */
  Gaussian.prototype.freq = function (line) {
    const F = this.fr;
    let p = 0, c0 = line.charCodeAt(0);
    while (c0 === 32 || c0 === 9) c0 = line.charCodeAt(++p);
    if (F.block && !F.hp && c0 >= 48 && c0 <= 57) {   // the hot path: a standard row of 2-decimal values
      const bl = F.block.length, a1 = rowFast(line, 3 * bl, 2);
      if (a1 > 0) {
        const R = this.run(), V = R.vecs, a = a1 - 1, base = 3 * a;
        if (a1 > R.nAtoms) R.nAtoms = a1;
        F.rowsIn = true; F.lastA = a;
        for (let k = 0; k < bl; k++) {
          const v = V[k];
          if (!v) continue;
          if (v.d === 2 && v.a instanceof Int8Array) {
            if (base + 2 >= v.a.length) v.room(base + 2);
            const A = v.a;
            for (let c = 0; c < 3; c++) {
              const m = ROWV[3 * k + c];
              if (m > 127 || m < -128) v.set(base + c, m, 2); else A[base + c] = m;
            }
            if (base + 3 > v.len) v.len = base + 3;
          } else for (let c = 0; c < 3; c++) v.set(base + c, ROWV[3 * k + c], 2);
        }
        return;
      }
    }
    if (F.block && c0 >= 48 && c0 <= 57) {   // fast path for atom rows (no other kind of line in a block looks like one)
      const bl = F.block.length, n = spans(line, 2 + 3 * bl + 1);
      const point = q => { const k = line.indexOf('.', SP[2 * q]); return k >= 0 && k < SP[2 * q + 1]; };   // (all integers: a header)
      if (n === 2 + 3 * bl) {
        const a = intAt(line, SP[0], SP[1]) - 1;
        if (!Number.isNaN(a) && !Number.isNaN(intAt(line, SP[2], SP[3])) && point(2)) return this.stdRow(line, a, bl);
      } else if (F.hp && n === 3 + bl) {
        const cc = intAt(line, SP[0], SP[1]) - 1, a = intAt(line, SP[2], SP[3]) - 1;
        if (!Number.isNaN(intAt(line, SP[4], SP[5])) && cc >= 0 && cc < 3 && a >= 0 && point(3)) return this.hpRow(line, cc, a, bl);
      }
    }
    const t = line.trim();
    if (!t) { if (F.freqs) return this.endFreq(line); return this.keepFreq(line); }
    if (!F.block && G_FREQ_TEXT.test(t)) return this.keepFreq(line);
    const tk = t.split(/\s+/);
    if (tk.every(isIntTok)) {
      this.flushRun();
      F.block = tk.map(x => int(x) - 1); F.expectSym = true; F.rowsIn = false; F.odd = 0; F.lastA = -1;
      return this.keepFreq(line);
    }
    const lm = G_LABEL.exec(t);
    if (lm && F.block) {
      F.expectSym = false;
      if (lm[2] === '---') F.hp = true;
      if (lm[1].toLowerCase().startsWith('frequencies')) F.freqs = true;
      return this.keepFreq(line);
    }
    if (F.expectSym && F.block && tk.length === F.block.length) { F.expectSym = false; return this.keepFreq(line); }
    if (/^Atom\s+AN\b/.test(t) || /^Coord\s+Atom\s+Element/.test(t)) return this.keepFreq(line);
    if (F.block && tk.length >= 3 && (isIntTok(tk[0]) || (STARS.test(tk[0]) && F.lastA >= 0)) && isIntTok(tk[1])) {
      if (F.hp && isIntTok(tk[0]) && isIntTok(tk[2]) && tk.length === 3 + F.block.length) {
        const cc = int(tk[0]) - 1, a = int(tk[1]) - 1;
        if (cc >= 0 && cc < 3 && a >= 0) return this.hpRow(tk, cc, a, F.block.length);
      }
      if (tk.length === 2 + 3 * F.block.length) return this.stdRow(tk, isIntTok(tk[0]) ? int(tk[0]) - 1 : F.lastA + 1, F.block.length);
    }
    if (/^Harmonic frequencies/.test(t)) return this.endFreq(line);
    if (F.block && !F.rowsIn && ++F.odd <= 8) return this.keepFreq(line);   // as the reader: an unknown line among the labels
    return this.endFreq(line);
  };
  Gaussian.prototype.keepFreq = function (line) { this.flushRun(); this.keep(line); };
  Gaussian.prototype.run = function () {
    const F = this.fr;
    if (!F.run) {
      F.run = { vecs: F.block.map(m => (!F.skip && this.wantMode(m) ? new Vec(3 * (this.modeN || 64)) : null)), nAtoms: 0 };
      this.modesSeen += F.block.length;
    }
    return F.run;
  };
  /** value token q of the row (the line with its spans in SP, or a token array) into vector v at index i */
  const putValue = (v, i, src, q) => {
    if (typeof src === 'string') { const s0 = SP[2 * q], s1 = SP[2 * q + 1]; if (dec(src, s0, s1)) v.set(i, DM, DD); else v.setRaw(i, num(src.slice(s0, s1))); }
    else { const x = src[q]; if (dec(x, 0, x.length)) v.set(i, DM, DD); else v.setRaw(i, num(x)); }
  };
  Gaussian.prototype.stdRow = function (src, a, bl) {
    const R = this.run(), V = R.vecs;
    this.fr.rowsIn = true; this.fr.lastA = a;
    if (a + 1 > R.nAtoms) R.nAtoms = a + 1;
    for (let k = 0; k < bl; k++) { const v = V[k]; if (v) for (let c = 0; c < 3; c++) putValue(v, 3 * a + c, src, 2 + 3 * k + c); }
  };
  Gaussian.prototype.hpRow = function (src, cc, a, bl) {
    const R = this.run(), V = R.vecs;
    this.fr.rowsIn = true;
    if (a + 1 > R.nAtoms) R.nAtoms = a + 1;
    for (let k = 0; k < bl; k++) { const v = V[k]; if (v) putValue(v, 3 * a + cc, src, 3 + k); }
  };
  Gaussian.prototype.flushRun = function () {
    const F = this.fr;
    if (!F || !F.run) return;
    const R = F.run;
    F.run = null;
    const vecs = R.vecs.map(v => (v ? v.done() : null));
    for (const v of vecs) if (v) { this.ext += v.vals.byteLength; this.modesStored++; }
    if (R.nAtoms) this.modeN = R.nAtoms;
    const k = this.side.modes.length;
    this.side.modes.push({ nAtoms: R.nAtoms, vecs });
    this.out.push(PH + 'M' + k);
    this.dropping = false;
  };
  Gaussian.prototype.endFreq = function (line) {
    this.flushRun();
    if (this.fr.hp) this.hpAt = this.side.geoms.length;
    this.fr = null;
    this.mode = M_PASS;
    this.closeStep();
    // a line that starts with a digit would be dropped as "0", which the reader would take for the mode numbers of a
    // next block: the line that ends the section is kept as it is, so the section ends there for the reader too
    const n = this.out.length;
    this.pass(line);
    if (this.out.length === n + 1 && this.out[n] === DROPPED) { this.out[n] = own(line); this.chars += line.length; this.dropping = false; }
  };
  Gaussian.prototype.close = function () {
    if (this.mode === M_OHEAD || this.mode === M_OROWS) { if (!this.rows) { this.rows = new Rows(16); this.tv = []; } this.endTable(false); }
    if (this.mode === M_FREQ) { this.flushRun(); this.fr = null; }
    this.mode = M_PASS;
  };

  // ================================================================== ORCA
  const O_PASS = 0, O_COORD = 1, O_MODES = 2, O_INTERNAL = 3;
  const O_CTX = new Set(['CD SPECTRUM', 'CHEMICAL SHIELDING SUMMARY (ppm)', 'VIBRATIONAL FREQUENCIES', 'IR SPECTRUM', 'RAMAN SPECTRUM',
    'ABSORPTION SPECTRUM VIA TRANSITION ELECTRIC DIPOLE MOMENTS', 'HIRSHFELD ANALYSIS', 'ORBITAL ENERGIES']);
  function Orca(opts) {
    Base.call(this, 'orca-out', opts);
    this.mode = O_PASS;
    this.ctx = null;
    this.rows = null; this.first = false; this.labels = []; this.labelCode = new Map(); this.codeCache = [];
    this.lastN = 64;
    this.md = null;        // normal-mode section state
    this.freqRank = null;  // mode index -> rank among the vibrations (non-zero frequencies), from VIBRATIONAL FREQUENCIES
    this.inFreqs = false;
  }
  Orca.prototype = Object.create(Base.prototype);

  Orca.prototype.push = function (line) {
    switch (this.mode) {
      case O_PASS: return this.pass(line);
      case O_COORD: return this.coord(line);
      case O_MODES: return this.modes(line);
      case O_INTERNAL:
        if (line.trim() && !/^\s*-+\s*$/.test(line)) return this.drop();
        this.mode = O_PASS; return this.pass(line);
    }
  };
  Orca.prototype.pass = function (line) {
    let p = 0, c0 = line.charCodeAt(0);
    while (c0 === 32 || c0 === 9) c0 = line.charCodeAt(++p);
    if (c0 >= 48 && c0 <= 57) {   // the reader skips such rows; section readers take theirs
      if (this.ctx) {
        this.ctx.digit = true;
        if (this.inFreqs) { const m = /^\s*(\d+):\s+(-?\d+\.\d*)\s+cm\*\*-1/.exec(line); if (m) this.freqRow(int(m[1]), num(m[2])); }
        return this.keep(line);
      }
      return this.drop();
    }
    const t = line.trim();
    if (this.ctx && (this.ctx.digit || ++this.ctx.idle > 30)) { this.ctx = null; this.inFreqs = false; }
    if (!t) return this.keep(line);
    const c = t[0];
    if (c === 'C' && t === 'CARTESIAN COORDINATES (ANGSTROEM)') { this.keep(line); this.mode = O_COORD; this.first = true; this.rows = null; return; }
    if (c === 'N' && t === 'NORMAL MODES') { this.keep(line); this.mode = O_MODES; this.md = { head: null, run: null, blank: false }; return; }
    if (c === 'I' && (t === 'INTERNAL COORDINATES (ANGSTROEM)' || t === 'INTERNAL COORDINATES (A.U.)')) { this.keep(line); this.mode = O_INTERNAL; return; }
    if (O_CTX.has(t) || (c === 'C' && t.startsWith('CD SPECTRUM VIA TRANSITION')) || (c === 'S' && (t.startsWith('STATE') || /SPIN (UP|DOWN) ORBITALS/.test(t))) ||
        (c === 'M' && t.startsWith('MULLIKEN ATOMIC CHARGES')) || (c === 'L' && t.startsWith('LOEWDIN ATOMIC CHARGES')) || (c === '*' && t.indexOf('IRC PATH SUMMARY') >= 0)) {
      this.keep(line);
      this.ctx = { digit: false, idle: 0 };
      if (t === 'VIBRATIONAL FREQUENCIES') { this.inFreqs = true; this.freqRank = []; this.nonZero = 0; }
      return;
    }
    if (c === '*' && (t.indexOf('GEOMETRY OPTIMIZATION CYCLE') >= 0 || t.indexOf('RELAXED SURFACE SCAN STEP') >= 0)) this.closeStep();
    // the final evaluation repeats the converged structure: it belongs to the last cycle (the reader merges it there)
    if (c === '*' && t.indexOf('FINAL ENERGY EVALUATION AT THE STATIONARY POINT') >= 0 && this.step) this.closed = false;
    if (c === '$' && /JOB NUMBER\s+(\d+)/.test(t)) this.newJob();
    this.keep(line);
  };
  Orca.prototype.freqRow = function (k, f) {
    if (!this.freqRank) return;
    this.freqRank[k] = Math.abs(f) > 1e-3 ? this.nonZero++ : -1;
  };
  /** CARTESIAN COORDINATES (ANGSTROEM): rows "label x y z" (label codes interned), then the reader's end line. */
  Orca.prototype.coord = function (line) {
    if (this.first) {
      this.first = false;
      if (/^-+$/.test(line.trim())) return this.keep(line);   // the reader skips one dashes line
    }
    const n = spans(line, 6);
    let label = null, x, y, z;
    if (n === 4 && dec(line, SP[2], SP[3])) {
      const c = line.charCodeAt(SP[0]);
      if (((c >= 65 && c <= 90) || (c >= 97 && c <= 122)) && line.charCodeAt(SP[2]) !== 58 && line.charCodeAt(SP[2]) !== 62) {
        x = decValue();
        if (dec(line, SP[4], SP[5])) { y = decValue(); if (dec(line, SP[6], SP[7])) { z = decValue(); label = line.slice(SP[0], SP[1]); } }
      }
    }
    if (label === null) {   // anything else exactly as the reader: ghost "H :", embedding "Zr>", point charges, exponents
      let tk = tok(line);
      if (tk.length > 1 && (tk[1] === ':' || tk[1] === '>')) tk = [tk[0] + tk[1]].concat(tk.slice(2));
      const k = tk.length;
      if (k < 4 || k > 5 || !/^[A-Za-z]/.test(tk[0]) || !tk.slice(1).every(isNumTok)) { this.endCoord(); return this.pass(line); }
      label = tk[0]; x = num(tk[k - 3]); y = num(tk[k - 2]); z = num(tk[k - 1]);
    }
    let code = this.labelCode.get(label);
    if (code === undefined) { code = this.labels.length; this.labels.push(own(label)); this.labelCode.set(this.labels[code], code); }
    if (!this.rows) this.rows = new Rows(this.lastN);
    this.rows.add(code, x, y, z);
  };
  Orca.prototype.endCoord = function () {
    const r = this.rows;
    this.rows = null;
    this.mode = O_PASS;
    if (!r || !r.n) return;
    this.lastN = Math.max(16, r.n);
    let codes = null;
    for (const c of this.codeCache) if (c.length === r.n && c.every((v, k) => v === r.z[k])) { codes = c; break; }
    if (!codes) { codes = Uint16Array.from(r.z.subarray(0, r.n)); this.codeCache.unshift(codes); if (this.codeCache.length > 4) this.codeCache.pop(); }
    this.geometry({ n: r.n, codes, labels: this.labels, xyz: r.x.slice(0, 3 * r.n) });
  };
  /** NORMAL MODES: blocks of column numbers followed by rows "coordinate value value ..." (6 columns). */
  Orca.prototype.modes = function (line) {
    const D = this.md;
    if (D.blank) {   // a blank line inside the section: the reader looks at the next line to know whether it ends
      D.blank = false;
      const nx = line.trim();
      this.keep('');
      if (!nx || !/^\d/.test(nx)) { this.mode = O_PASS; this.md = null; return this.pass(line); }
    }
    let p = 0, c0 = line.charCodeAt(0);
    while (c0 === 32 || c0 === 9) c0 = line.charCodeAt(++p);
    if (D.head && c0 >= 48 && c0 <= 57) {   // fast path for rows (plain decimals, at least one with a point)
      const bl = D.head.length, n = spans(line, bl + 2);
      if (n === bl + 1) {
        const r = intAt(line, SP[0], SP[1]);
        if (!Number.isNaN(r)) {
          let ok = true, point = false;
          for (let q = 1; q <= bl && ok; q++) { ok = dec(line, SP[2 * q], SP[2 * q + 1]); if (DD) point = true; }
          if (ok && point) return this.modeRow(line, r, bl);
        }
      }
    }
    const s = line.trim();
    if (!s) {
      if (D.head) { this.flushModes(); D.blank = true; return; }   // (the rows' placeholder must come before the blank)
      return this.keep(line);
    }
    const tk = s.split(/\s+/);
    if (tk.every(isIntTok)) {
      this.flushModes();
      D.head = tk.map(int);
      return this.keep(line);
    }
    if (D.head && isIntTok(tk[0]) && tk.length === D.head.length + 1 && tk.slice(1).every(isNumTok)) return this.modeRow(tk, int(tk[0]), D.head.length);
    if (D.head) { this.flushModes(); this.mode = O_PASS; this.md = null; return this.pass(line); }
    this.keep(line);
  };
  /** rank of mode k among the vibrations (-1: a translation or rotation, which the reader drops) */
  Orca.prototype.rank = function (k) { return this.freqRank && this.freqRank[k] !== undefined ? this.freqRank[k] : k; };
  Orca.prototype.modeRow = function (src, r, bl) {
    const D = this.md;
    if (!D.run) {
      D.run = { vecs: D.head.map(k => { const q = this.rank(k); return q >= 0 && this.wantMode(q) ? new Vec(3 * (this.modeN || 64)) : null; }), rows: 0 };
      this.modesSeen += D.head.filter(k => this.rank(k) >= 0).length;
    }
    const R = D.run;
    if (r + 1 > R.rows) R.rows = r + 1;
    for (let q = 0; q < bl; q++) { const v = R.vecs[q]; if (v) putValue(v, r, src, q + 1); }
  };
  Orca.prototype.flushModes = function () {
    const D = this.md;
    if (!D || !D.run) return;
    const R = D.run;
    D.run = null;
    const vecs = R.vecs.map(v => (v ? v.done() : null));
    for (const v of vecs) if (v) { this.ext += v.vals.byteLength; this.modesStored++; }
    this.modeN = Math.ceil(R.rows / 3);
    const k = this.side.modes.length;
    this.side.modes.push({ nAtoms: this.modeN, vecs });
    this.out.push(PH + 'M' + k);
    this.dropping = false;
  };
  Orca.prototype.close = function () {
    if (this.mode === O_COORD) this.endCoord();
    if (this.mode === O_MODES) { this.flushModes(); if (this.md && this.md.blank) this.keep(''); this.md = null; }
    this.mode = O_PASS;
  };

  // ================================================================== multi-frame XYZ
  // parseXYZ: a count line, a comment line and that many atom rows, repeated (blank lines between frames are
  // skipped). The atom rows of each frame are compacted; frames that are not selected are left out entirely.
  const X_START = 0, X_COUNT = 1, X_COMMENT = 2, X_ROWS = 3, X_COPY = 4, X_DONE = 5;
  function Xyz(opts) {
    Base.call(this, 'xyz', opts);
    this.state = X_START;
    this.head = [];        // lines before the first frame
    this.frames = [];      // {count, comment, g} of the frames kept (assembled in finish: 'last' needs the end)
    this.tail = [];        // what the reader reads after the last frame (a bad or incomplete frame), as it is
    this.cur = null;
    this.nFrames = 0;
    this.labels = []; this.lab = new Map(); this.codeCache = [];
  }
  Xyz.prototype = Object.create(Base.prototype);
  Xyz.prototype.want = function (k) {
    const o = this.opts.steps;
    return o.mode === 'all' || o.mode === 'last' || (o.mode === 'every' ? k % o.n === 0 : k + 1 >= o.from && k + 1 <= o.to);
  };
  Xyz.prototype.push = function (line) {
    switch (this.state) {
      case X_DONE: return;   // the reader stops reading before this
      case X_COPY: this.tail.push(own(line)); this.chars += line.length; return;
      case X_START: case X_COUNT: {
        if (!line.trim()) { if (this.state === X_START) this.head.push(own(line)); return; }
        const h = tok(line);
        if (this.state === X_START && !(isIntTok(h[0]) && (h.length === 1 || !isNumTok(h[1])))) { this.state = X_COPY; return this.push(line); }   // one structure without a count line
        if (!isIntTok(h[0])) { this.tail.push(own(line)); this.state = X_DONE; return; }   // the reader stops here
        this.cur = { count: own(line), na: int(h[0]), comment: null, rows: null, k: this.nFrames++ };
        this.state = X_COMMENT;
        return;
      }
      case X_COMMENT:
        this.cur.comment = own(line);
        if (this.cur.na > 0) this.state = X_ROWS; else this.endFrame();
        return;
    }
    // atom rows "El x y z ..." (the reader converts the element, with its warnings for unknown symbols)
    const F = this.cur, n = spans(line, 64);
    let x, y, z, label = null;
    if (n >= 4 && dec(line, SP[2], SP[3])) {
      x = decValue();
      if (dec(line, SP[4], SP[5])) { y = decValue(); if (dec(line, SP[6], SP[7])) { z = decValue(); label = line.slice(SP[0], SP[1]); } }
    }
    if (label === null) {
      const t = tok(line);
      if (t.length < 4 || !isNumTok(t[1]) || !isNumTok(t[2]) || !isNumTok(t[3])) {   // the reader skips this frame and stops
        this.tail.push(F.count, F.comment, own(line));
        this.cur = null; this.state = X_DONE;
        return;
      }
      label = t[0]; x = num(t[1]); y = num(t[2]); z = num(t[3]);
    }
    let code = this.lab.get(label);
    if (code === undefined) { code = this.labels.length; this.labels.push(own(label)); this.lab.set(this.labels[code], code); }
    if (!F.rows) F.rows = new Rows(Math.min(F.na, 1 << 20));
    F.rows.add(code, x, y, z);
    if (F.rows.n === F.na) this.endFrame();
  };
  Xyz.prototype.endFrame = function () {
    const F = this.cur;
    this.cur = null;
    this.state = X_COUNT;
    if (!this.want(F.k)) return;
    const r = F.rows || new Rows(1);
    let codes = null;
    for (const c of this.codeCache) if (c.length === r.n && c.every((v, k) => v === r.z[k])) { codes = c; break; }
    if (!codes) { codes = Uint16Array.from(r.z.subarray(0, r.n)); this.codeCache.unshift(codes); if (this.codeCache.length > 4) this.codeCache.pop(); }
    const g = { n: r.n, codes, labels: this.labels, xyz: r.x.slice(0, 3 * r.n) };
    this.frames.push({ count: F.count, comment: F.comment, g });
    this.ext += g.xyz.byteLength;
    if (this.opts.steps.mode === 'last' && this.frames.length > this.opts.steps.n) this.ext -= this.frames.shift().g.xyz.byteLength;
  };
  Xyz.prototype.finish = function (lastLine) {
    if (lastLine) this.push(lastLine);   // (not an output: a last line without a line break is a complete line)
    if (this.cur) {   // the file ends inside a frame (an MD still running): the reader reports it as incomplete
      this.tail.push(this.cur.count);
      if (this.cur.comment !== null) this.tail.push(this.cur.comment);
      this.cur = null;
    }
    const out = this.head.slice();
    for (const f of this.frames) {
      out.push(f.count, f.comment, PH + 'G' + this.side.geoms.length);
      this.side.geoms.push(f.g);
    }
    for (const t of this.tail) out.push(t);
    out.push('');
    this.out = out;
    const o = this.opts.steps;
    this.side.loaded = { steps: Object.assign({}, o, { total: this.nFrames, loaded: this.frames.length, frames: true, to: Number.isFinite(o.to) ? o.to : null }), modes: null, bytes: this.ext };
    Object.defineProperty(out, 'bmSide', { value: this.side, configurable: true });
    Object.defineProperty(out, 'bmPartial', { value: this.partial, configurable: true });
    return out;
  };
  Xyz.prototype.lineCount = function () { return 3 * this.frames.length + this.head.length + this.tail.length; };

  const MAKERS = { 'gaussian-log': o => new Gaussian(o), 'orca-out': o => new Orca(o), xyz: o => new Xyz(o) };
  BM.condense = {
    formats: Object.keys(MAKERS),
    create: (format, opts) => (MAKERS[format] ? MAKERS[format](opts || {}) : null),
    own,
    /** Condenses a whole line array (tests, and files that are already in memory). */
    lines(format, lines, opts) {
      const c = this.create(format, opts);
      if (!c) return null;
      const n = lines.length;
      for (let k = 0; k < n - 1; k++) c.push(lines[k]);
      return c.finish(n ? lines[n - 1] : '');
    }
  };
})(typeof window !== 'undefined' ? (window.BM = window.BM || {}) : (globalThis.BM = globalThis.BM || {}));
