/* BlueMolecule — spectrum simulation (IR, Raman, UV-Vis) and a small canvas/SVG plotting component. */
(function (BM) {
  'use strict';
  const LN10 = Math.log(10);
  const EV_TO_CM = 8065.544;
  const HC_K = 1.438777; // cm K

  function lineshape(shape, dx, fwhm) {
    if (shape === 'gaussian') {
      const sig = fwhm / (2 * Math.sqrt(2 * Math.log(2)));
      return Math.exp(-dx * dx / (2 * sig * sig)) / (sig * Math.sqrt(2 * Math.PI));
    }
    const g = fwhm / 2;
    return g / (Math.PI * (dx * dx + g * g));
  }

  function linspace(a, b, n) { const out = new Float64Array(n); for (let i = 0; i < n; i++) out[i] = a + (b - a) * i / (n - 1); return out; }

  function broaden(sx, sy, opt) {
    const n = opt.n || 1600, xs = linspace(opt.xmin, opt.xmax, n), ys = new Float64Array(n);
    const cut = opt.shape === 'gaussian' ? opt.fwhm * 4 : opt.fwhm * 60;
    for (let t = 0; t < sx.length; t++) {
      const x0 = sx[t], h = sy[t];
      if (!isFinite(x0) || !isFinite(h) || h === 0) continue;
      for (let i = 0; i < n; i++) { const dx = xs[i] - x0; if (dx > cut || dx < -cut) continue; ys[i] += h * lineshape(opt.shape, dx, opt.fwhm); }
    }
    return { x: xs, y: ys };
  }
  const maxOf = a => { let m = 0; for (const v of a) if (isFinite(v) && v > m) m = v; return m; };

  /** IR: returns {curve:{x,y}, sticks:[{x,y,id}], yLabel} */
  function ir(vib, o) {
    const sx = [], sy = [], ids = [];
    vib.freqs.forEach((f, i) => { if (f > 0 && vib.ir && isFinite(vib.ir[i])) { sx.push(f * o.scale); sy.push(vib.ir[i]); ids.push(i); } });
    const curve = broaden(sx, sy, o);
    let yLabel;
    if (o.yMode === 'epsilon') { const k = 100 / LN10; for (let i = 0; i < curve.y.length; i++) curve.y[i] *= k; yLabel = 'ε (L mol⁻¹ cm⁻¹)'; }
    else { const m = maxOf(curve.y) || 1; for (let i = 0; i < curve.y.length; i++) curve.y[i] /= m; yLabel = 'Relative absorbance'; }
    return { curve, sticks: sx.map((x, t) => ({ x, y: sy[t], id: ids[t] })), yLabel };
  }

  function raman(vib, o) {
    const sx = [], sy = [], ids = [];
    const nu0 = 1e7 / (o.laser || 532);
    const T = o.temperature || 298.15;
    vib.freqs.forEach((f, i) => {
      if (!(f > 0) || !vib.raman || !isFinite(vib.raman[i])) return;
      const nu = f * o.scale;
      let y = vib.raman[i];
      if (o.yMode === 'intensity') y = Math.pow(nu0 - nu, 4) * y / (nu * (1 - Math.exp(-HC_K * nu / T)));
      sx.push(nu); sy.push(y); ids.push(i);
    });
    const curve = broaden(sx, sy, o);
    const m = maxOf(curve.y) || 1;
    for (let i = 0; i < curve.y.length; i++) curve.y[i] /= m;
    return { curve, sticks: sx.map((x, t) => ({ x, y: sy[t], id: ids[t] })), yLabel: o.yMode === 'intensity' ? 'Relative intensity' : 'Relative activity' };
  }

  /** UV-Vis from oscillator strengths. o.unit = 'nm'|'eV', o.sigmaEV = half-width at 1/e (eV). */
  function uv(states, o) {
    const n = o.n || 1200;
    const xs = linspace(o.xmin, o.xmax, n), ys = new Float64Array(n);
    const sig = (o.sigmaEV || 0.4) * EV_TO_CM;
    const st = states.filter(s => isFinite(s.energyEV) && s.f !== null && isFinite(s.f));
    for (let i = 0; i < n; i++) {
      const nu = o.unit === 'eV' ? xs[i] * EV_TO_CM : 1e7 / xs[i];
      let e = 0;
      for (const s of st) { const d = (nu - s.energyEV * EV_TO_CM) / sig; if (d * d < 40) e += 1.3062974e8 * s.f / sig * Math.exp(-d * d); }
      ys[i] = e;
    }
    const sticks = st.map(s => ({ x: o.unit === 'eV' ? s.energyEV : 1239.84198 / s.energyEV, y: s.f, id: s.index }));
    return { curve: { x: xs, y: ys }, sticks, yLabel: 'ε (L mol⁻¹ cm⁻¹)' };
  }

  /** ECD from rotatory strengths R (10⁻⁴⁰ cgs): Δε(E) = Σ E_k R_k exp(−((E − E_k)/σ)²) / (22.96 √π σ), E and σ in eV.
   *  o.unit = 'nm'|'eV', o.sigmaEV = half-width at 1/e (eV). Sticks carry the signed R. */
  function ecd(states, o) {
    const n = o.n || 1200;
    const xs = linspace(o.xmin, o.xmax, n), ys = new Float64Array(n);
    const sig = o.sigmaEV || 0.3, pref = 1 / (22.96 * Math.sqrt(Math.PI) * sig);
    const st = states.filter(s => isFinite(s.energyEV) && s.R !== null && s.R !== undefined && isFinite(s.R));
    for (let i = 0; i < n; i++) {
      const E = o.unit === 'eV' ? xs[i] : 1239.84198 / xs[i];
      let y = 0;
      for (const s of st) { const d = (E - s.energyEV) / sig; if (d * d < 40) y += pref * s.energyEV * s.R * Math.exp(-d * d); }
      ys[i] = y;
    }
    const sticks = st.map(s => ({ x: o.unit === 'eV' ? s.energyEV : 1239.84198 / s.energyEV, y: s.R, id: s.index }));
    return { curve: { x: xs, y: ys }, sticks, yLabel: 'Δε (L mol⁻¹ cm⁻¹)' };
  }

  // ---------- plotting ----------
  function niceStep(range, target) {
    const raw = range / Math.max(1, target), mag = Math.pow(10, Math.floor(Math.log10(raw))), f = raw / mag;
    return (f < 1.5 ? 1 : f < 3 ? 2 : f < 7 ? 5 : 10) * mag;
  }
  function ticks(a, b, target) {
    const lo = Math.min(a, b), hi = Math.max(a, b);
    if (!(hi > lo)) return [lo];
    const st = niceStep(hi - lo, target), out = [];
    for (let v = Math.ceil(lo / st) * st; v <= hi + st * 1e-9; v += st) out.push(Math.abs(v) < st * 1e-9 ? 0 : v);
    return out;
  }
  function fmt(v, step) {
    const d = step >= 1 ? 0 : Math.min(6, Math.ceil(-Math.log10(step) - 1e-9));
    return v.toFixed(d);
  }
  const SUP = { '-': '⁻', 0: '⁰', 1: '¹', 2: '²', 3: '³', 4: '⁴', 5: '⁵', 6: '⁶', 7: '⁷', 8: '⁸', 9: '⁹' };
  const decade = k => (k === 0 ? '1' : k === 1 ? '10' : '10' + String(k).split('').map(c => SUP[c]).join(''));
  const esc = s => String(s).replace(/&/g, '&amp;').replace(/</g, '&lt;').replace(/>/g, '&gt;');
  const textW = (t, size) => String(t).length * size * 0.56;
  const PALETTE = ['#4f8cff', '#ff8a4c', '#34c98a', '#d66bd9', '#e2c044', '#4fc6d9', '#ff6b8a', '#9aa5b8'];
  const PRINT = { bg: '#ffffff', fg: '#1f2937', grid: '#d5dbe5', accent: '#1d4ed8', stick: '#6b7280', sel: '#d97706', ok: '#15803d', legendBg: 'rgba(255,255,255,0.92)',
    palette: ['#1d4ed8', '#ea580c', '#15803d', '#a21caf', '#b45309', '#0e7490', '#be123c', '#475569'] };

  /** Canvas plot of spectra, profiles and bar/level charts. spec (all optional):
   *   xRange [a,b], reverseX, integerX, xLabel, yLabel, yRange, logY, invertY, xFormat(x), empty (message)
   *   series [{x, y, fill, color, width}]         curves (spectra)
   *   sticks [{x, y, id, label}], showSticks, peakLabels, peakFormat(x), selected
   *   points [{x, y, id, mark, label}]             one profile with markers (energy vs geometry)
   *   lines  [{name, color, pts: [...points]}]     several profiles; named lines get a legend
   *   hlines [{y, label, color}]                   dashed horizontal lines (thresholds)
   *   bars   [{x, y, id, label, color}], barWidth, xTickLabels [{x, label}]
   *   type: 'levels', columns [{name, levels: [{y, occ, id, label, tag}]}], gaps [{col, from, to, text}]
   * onPick(id, item) is called when a point, stick, bar or level is clicked. */
  class Plot {
    constructor(canvas) {
      this.canvas = canvas; this.ctx = canvas.getContext('2d');
      this.spec = null; this.hover = null; this.onPick = null;
      this.m = { l: 64, r: 18, t: 16, b: 40 };
      canvas.addEventListener('mousemove', e => this._move(e));
      canvas.addEventListener('mouseleave', () => { this.hover = null; this.draw(); });
      canvas.addEventListener('click', e => this._click(e));
      if (typeof ResizeObserver !== 'undefined') new ResizeObserver(() => this.draw()).observe(canvas.parentElement);
    }
    set(spec) { this.spec = spec; this.hover = null; this.draw(); }
    _lines() {
      const s = this.spec;
      if (s.lines) return s.lines;
      return s.points && s.points.length ? [{ pts: s.points, main: true }] : [];
    }
    _layout(w, h) {
      const s = this.spec, m = this.m;
      let [x0, x1] = s.xRange || [0, 1];
      if (s.reverseX) [x0, x1] = [x1, x0];
      let y0, y1;
      if (s.yRange) [y0, y1] = s.yRange;
      else {
        let mn = Infinity, mx = -Infinity;
        const add = v => { if (isFinite(v) && (!s.logY || v > 0)) { if (v < mn) mn = v; if (v > mx) mx = v; } };
        for (const ser of s.series || []) for (let i = 0; i < ser.y.length; i++) add(ser.y[i]);
        for (const ln of this._lines()) for (const p of ln.pts) add(p.y);
        for (const b of s.bars || []) add(b.y);
        for (const hl of s.hlines || []) add(hl.y);
        if (s.type === 'levels') for (const c of s.columns || []) for (const l of c.levels) add(l.y);
        if (!(s.series && s.series.length) && s.sticks) for (const p of s.sticks) add(p.y);
        const zero = !s.logY && ((s.series && s.series.length) || (s.sticks && s.sticks.length) || (s.bars && s.bars.length));
        if (zero) { add(0); }
        if (!isFinite(mn)) { mn = s.logY ? 1e-6 : 0; mx = s.logY ? 1 : 1; }
        if (s.logY) {
          let a = Math.log10(mn), b = Math.log10(mx);
          if (b - a < 1) { const c = (a + b) / 2; a = c - 0.5; b = c + 0.5; }
          y0 = Math.pow(10, a - (b - a) * 0.06); y1 = Math.pow(10, b + (b - a) * 0.08);
        } else {
          if (mx === mn) { const d = Math.abs(mx) * 0.01 || 1; mn -= d; mx += d; }
          const pad = (mx - mn) * 0.08;
          y0 = mn - (zero && mn >= 0 ? 0 : pad); y1 = mx + pad * (s.peakLabels ? 2.2 : 1);
          if (s.legend !== false && this._lines().filter(l => l.name).length > 1) y1 += (mx - mn) * 0.12;
        }
      }
      const pw = Math.max(10, w - m.l - m.r), ph = Math.max(10, h - m.t - m.b);
      const X = x => m.l + (x - x0) / (x1 - x0) * pw;
      let Y;
      if (s.logY) { const la = Math.log10(y0), lb = Math.log10(y1); Y = y => (y > 0 ? m.t + ph - (Math.log10(y) - la) / (lb - la) * ph : NaN); }
      else Y = s.invertY ? (y => m.t + (y - y0) / (y1 - y0) * ph) : (y => m.t + ph - (y - y0) / (y1 - y0) * ph);
      return { x0, x1, y0, y1, pw, ph, X, Y };
    }
    _stickScale(L) {
      const s = this.spec;
      if (!s.sticks || !s.sticks.length) return 1;
      const ms = Math.max(...s.sticks.map(p => Math.abs(p.y))) || 1;
      if (s.series && s.series.length) {
        const top = s.sticks.some(p => p.y < 0) ? Math.max(L.y1, -L.y0) * 0.85 : Math.max(L.y1 - (L.y1 - L.y0) * 0.12, 0);
        return top / ms * (s.stickFraction || 0.9);
      }
      return 1;
    }
    _theme() {
      const cs = getComputedStyle(this.canvas);
      const v = (n, d) => (cs.getPropertyValue(n) || '').trim() || d;
      const accent = v('--accent', '#4f8cff');
      return { bg: v('--plot-bg', '#0d1420'), fg: v('--plot-fg', '#c9d4e3'), grid: v('--plot-grid', 'rgba(150,170,200,0.14)'), accent, stick: v('--plot-stick', '#8aa4c8'), sel: v('--plot-sel', '#ffcc40'), ok: v('--ok', '#33d17a'), legendBg: v('--plot-bg', '#0d1420'), palette: [accent].concat(PALETTE.slice(1)) };
    }
    /** Items under the cursor: the nearest point, stick, bar or level (pixel distance). */
    _hit(px, py) {
      const s = this.spec; if (!s || s.empty) return null;
      const w = this.canvas.clientWidth, h = this.canvas.clientHeight, L = this._layout(w, h);
      let best = null, bd = Infinity;
      if (s.type === 'levels') {
        for (const g of this._levelGeom(L)) { const dx = px < g.x0 ? g.x0 - px : px > g.x1 ? px - g.x1 : 0, d = Math.hypot(dx * 2, g.y - py); if (d < bd && d < 12) { bd = d; best = g.item; } }
        return best;
      }
      const lines = this._lines();
      if (lines.length) {
        const multi = lines.length > 1;
        for (const ln of lines) for (const p of ln.pts) {
          const X = L.X(p.x), Y = L.Y(p.y); if (!isFinite(Y)) continue;
          const d = multi ? Math.hypot(X - px, Y - py) : Math.abs(X - px);
          if (d < bd && d < (multi ? 14 : 9)) { bd = d; best = p; }
        }
        return best;
      }
      if (s.bars && s.bars.length) {
        const bw = Math.abs(L.X((s.barWidth || 0.7) / 2) - L.X(0)) + 2;
        for (const b of s.bars) { const d = Math.abs(L.X(b.x) - px); if (d < bd && d <= Math.max(bw, 4)) { bd = d; best = b; } }
        return best;
      }
      for (const p of s.sticks || []) { const d = Math.abs(L.X(p.x) - px); if (d < bd && d < 9) { bd = d; best = p; } }
      return best;
    }
    _move(e) {
      if (!this.spec) return;
      const r = this.canvas.getBoundingClientRect();
      const px = e.clientX - r.left, py = e.clientY - r.top;
      const hit = this._hit(px, py);
      const w = this.canvas.clientWidth, h = this.canvas.clientHeight, L = this._layout(w, h);
      const inside = px >= this.m.l && px <= w - this.m.r && py >= this.m.t && py <= h - this.m.b;
      const xVal = L.x0 + (px - this.m.l) / L.pw * (L.x1 - L.x0);
      this.hover = inside ? { px, py, hit, xVal } : null;
      this.canvas.style.cursor = hit ? 'pointer' : 'crosshair';
      this.draw();
    }
    _click(e) {
      const r = this.canvas.getBoundingClientRect();
      const hit = this._hit(e.clientX - r.left, e.clientY - r.top);
      if (hit && this.onPick && hit.id !== undefined && hit.id !== null) this.onPick(hit.id, hit);
    }
    draw() {
      const cv = this.canvas, dpr = Math.min(window.devicePixelRatio || 1, 3);
      const w = cv.clientWidth, h = cv.clientHeight;
      if (!w || !h) return;
      if (cv.width !== Math.round(w * dpr) || cv.height !== Math.round(h * dpr)) { cv.width = Math.round(w * dpr); cv.height = Math.round(h * dpr); }
      const ctx = this.ctx; ctx.setTransform(dpr, 0, 0, dpr, 0, 0);
      const th = this._theme();
      ctx.clearRect(0, 0, w, h);
      ctx.fillStyle = th.bg; ctx.fillRect(0, 0, w, h);
      const s = this.spec;
      if (!s) return;
      if (s.empty) { ctx.fillStyle = th.fg; ctx.globalAlpha = 0.6; ctx.font = '13px ui-sans-serif, -apple-system, sans-serif'; ctx.textAlign = 'center'; ctx.fillText(s.empty, w / 2, h / 2); ctx.globalAlpha = 1; return; }
      this._render(new CanvasPainter(ctx), w, h, th);
      if (this.hover) {
        const hv = this.hover;
        ctx.strokeStyle = th.grid; ctx.lineWidth = 1; ctx.setLineDash([3, 3]);
        if (s.type !== 'levels' && !(s.bars && s.bars.length)) { ctx.beginPath(); ctx.moveTo(hv.px, this.m.t); ctx.lineTo(hv.px, h - this.m.b); ctx.stroke(); }
        ctx.setLineDash([]);
        const text = hv.hit && hv.hit.label ? hv.hit.label : (s.type === 'levels' || (s.bars && s.bars.length) ? '' : s.xFormat ? s.xFormat(hv.xVal) : hv.xVal.toFixed(1));
        if (!text) return;
        ctx.font = '12px ui-sans-serif, -apple-system, sans-serif';
        const tw = ctx.measureText(text).width + 14;
        let tx = hv.px + 10; if (tx + tw > w - 4) tx = hv.px - tw - 10;
        if (tx < 2) tx = 2;
        ctx.fillStyle = 'rgba(12,18,28,0.92)';
        ctx.beginPath(); if (ctx.roundRect) ctx.roundRect(tx, this.m.t + 4, tw, 22, 5); else ctx.rect(tx, this.m.t + 4, tw, 22); ctx.fill();
        ctx.fillStyle = '#e8eef7'; ctx.textAlign = 'left'; ctx.textBaseline = 'middle'; ctx.fillText(text, tx + 7, this.m.t + 15);
        ctx.textBaseline = 'alphabetic';
      }
    }
    _axes(P, L, th, w, h) {
      const s = this.spec, m = this.m;
      const yt = [];
      if (s.logY) {
        const a = Math.ceil(Math.log10(L.y0) - 1e-9), b = Math.floor(Math.log10(L.y1) + 1e-9), n = b - a + 1, step = Math.max(1, Math.ceil(n / Math.max(2, Math.floor(L.ph / 26))));
        for (let k = a; k <= b; k += step) yt.push([Math.pow(10, k), decade(k)]);
      } else {
        const t = ticks(L.y0, L.y1, Math.max(2, Math.floor(L.ph / 45))), st = t.length > 1 ? Math.abs(t[1] - t[0]) : 1;
        for (const v of t) yt.push([v, fmt(v, st)]);
      }
      for (const [v, lab] of yt) { const y = L.Y(v); if (!isFinite(y)) continue; P.line(m.l, y, m.l + L.pw, y, th.grid, 1); P.text(lab, m.l - 8, y + 4, th.fg, 11, 'right'); }
      if (s.type === 'levels') {
        const cols = s.columns || [];
        cols.forEach((c, i) => P.text(c.name, m.l + L.pw * (i + 0.5) / cols.length, m.t + L.ph + 16, th.fg, 12, 'center'));
      } else if (s.xTickLabels) {
        const labs = s.xTickLabels, slot = labs.length > 1 ? Math.abs(L.X(labs[1].x) - L.X(labs[0].x)) : L.pw;
        const widest = Math.max(...labs.map(t => textW(t.label, 10.5))) + 6, stride = Math.max(1, Math.ceil(widest / Math.max(1, slot)));
        labs.forEach((t, k) => { if (k % stride === 0) P.text(t.label, L.X(t.x), m.t + L.ph + 15, th.fg, 10.5, 'center'); });
      } else {
        let xt = ticks(L.x0, L.x1, Math.max(3, Math.floor(L.pw / 80)));
        if (s.integerX) { xt = xt.filter(v => Math.abs(v - Math.round(v)) < 1e-9); if (xt.length < 2) xt = ticks(Math.ceil(Math.min(L.x0, L.x1)), Math.floor(Math.max(L.x0, L.x1)), 2).map(Math.round); }
        const xs = xt.length > 1 ? Math.abs(xt[1] - xt[0]) : 1;
        for (const v of xt) { const x = L.X(v); P.line(x, m.t, x, m.t + L.ph, th.grid, 1); P.text(fmt(v, xs), x, m.t + L.ph + 15, th.fg, 11, 'center'); }
      }
      P.rect(m.l, m.t, L.pw, L.ph, th.grid);
      if (s.xLabel) P.text(s.xLabel, m.l + L.pw / 2, h - 6, th.fg, 12, 'center');
      if (s.yLabel) P.textV(s.yLabel, 14, m.t + L.ph / 2, th.fg, 12);
    }
    /** Segment geometry of the level diagram: [{x0, x1, y, item, col}] (degenerate levels side by side). */
    _levelGeom(L) {
      const s = this.spec, cols = s.columns || [], out = [], m = this.m;
      const span = Math.abs(L.y1 - L.y0) || 1, tol = span * 0.004;
      cols.forEach((c, ci) => {
        const cx = m.l + L.pw * (ci + 0.5) / cols.length, hw = Math.min(L.pw / cols.length * 0.22, 60);
        const lv = c.levels.slice().sort((a, b) => a.y - b.y);
        for (let i = 0; i < lv.length;) {
          let j = i + 1;
          while (j < lv.length && lv[j].y - lv[i].y <= tol) j++;
          const k = j - i, gap = 4, sw = (2 * hw - gap * (k - 1)) / k;
          for (let t = 0; t < k; t++) { const x0 = cx - hw + t * (sw + gap); out.push({ x0, x1: x0 + sw, y: L.Y(lv[i + t].y), item: lv[i + t], col: ci, cx, hw }); }
          i = j;
        }
      });
      return out;
    }
    _legend(P, L, th, lines, colors) {
      const named = lines.map((l, k) => ({ name: l.name, color: colors[k] })).filter(l => l.name);
      if (named.length < 2 || this.spec.legend === false) return;
      const m = this.m, size = 10.5, rowH = 15, wMax = Math.max(...named.map(l => textW(l.name, size))) + 30;
      const cols = Math.max(1, Math.min(named.length, Math.floor((L.pw - 20) / (wMax + 8)))), rows = Math.ceil(named.length / cols);
      const bw = cols * (wMax + 8) + 6, bh = rows * rowH + 8, x0 = m.l + L.pw - bw - 6, y0 = m.t + 6;
      P.fillRect(x0, y0, bw, bh, th.legendBg, 0.85);
      named.forEach((l, k) => {
        const cx = x0 + 6 + (k % cols) * (wMax + 8), cy = y0 + 4 + Math.floor(k / cols) * rowH + rowH / 2;
        P.line(cx, cy, cx + 18, cy, l.color, 2);
        P.circle(cx + 9, cy, 2.6, l.color);
        P.text(l.name, cx + 23, cy + 4, th.fg, size, 'left');
      });
    }
    _render(P, w, h, th) {
      const s = this.spec, L = this._layout(w, h), m = this.m;
      this._axes(P, L, th, w, h);
      P.clip(m.l, m.t, L.pw, L.ph);
      const baseY = s.logY ? m.t + L.ph : L.Y(Math.max(L.y0, Math.min(0, L.y1)));
      // level diagram
      if (s.type === 'levels') {
        const geo = this._levelGeom(L);
        for (const g of geo) {
          const it = g.item, sel = s.selected !== undefined && s.selected !== null && s.selected === it.id;
          P.line(g.x0, g.y, g.x1, g.y, sel ? th.sel : it.occ ? th.accent : th.stick, sel ? 3 : it.occ ? 2.2 : 1.6);
        }
        for (const g of geo) if (g.item.tag) P.text(g.item.tag, g.cx + g.hw + 6, g.y + 4, th.fg, 10.5, 'left');
        for (const gp of s.gaps || []) {
          const g0 = geo.find(g => g.col === gp.col), cx = g0 ? g0.cx - g0.hw - 12 : m.l + 20, ya = L.Y(gp.from), yb = L.Y(gp.to);
          P.line(cx, ya, cx, yb, th.fg, 1, [3, 3]);
          P.line(cx - 4, ya, cx + 4, ya, th.fg, 1); P.line(cx - 4, yb, cx + 4, yb, th.fg, 1);
          P.text(gp.text, cx - 6, (ya + yb) / 2 + 4, th.fg, 10.5, 'right');
        }
        P.unclip();
        return;
      }
      // bars
      if (s.bars && s.bars.length) {
        const bw = Math.max(1, Math.abs(L.X((s.barWidth || 0.7)) - L.X(0)));
        const y0 = L.Y(0);
        for (const b of s.bars) {
          const x = L.X(b.x), y = L.Y(b.y), sel = s.selected !== undefined && s.selected !== null && s.selected === b.id;
          P.fillRect(x - bw / 2, Math.min(y, y0), bw, Math.abs(y - y0), sel ? th.sel : b.color || th.accent, 1);
        }
        P.line(m.l, y0, m.l + L.pw, y0, th.fg, 1);
      }
      // sticks
      const k = this._stickScale(L);
      if (s.sticks && s.showSticks !== false) for (const p of s.sticks) {
        const x = L.X(p.x), sel = s.selected !== undefined && s.selected === p.id;
        P.line(x, baseY, x, L.Y(p.y * k), sel ? th.sel : th.stick, sel ? 2.2 : 1.2);
      }
      // curves
      for (const ser of s.series || []) {
        const pts = [];
        for (let i = 0; i < ser.x.length; i++) pts.push([L.X(ser.x[i]), L.Y(ser.y[i])]);
        if (ser.fill) P.fillPath(pts, baseY, ser.fill);
        P.path(pts, ser.color || th.accent, ser.width || 1.8);
      }
      if (s.series && s.series.length && s.sticks && s.sticks.some(p => p.y < 0)) P.line(m.l, L.Y(0), m.l + L.pw, L.Y(0), th.grid, 1);
      // threshold lines
      for (const hl of s.hlines || []) {
        const y = L.Y(hl.y); if (!isFinite(y)) continue;
        P.line(m.l, y, m.l + L.pw, y, hl.color || th.sel, 1.3, [6, 4]);
        if (hl.label) P.text(hl.label, m.l + L.pw - 6, y - 4, hl.color || th.sel, 10.5, 'right');
      }
      // profiles (one or several lines with markers)
      const lines = this._lines();
      const colors = lines.map((ln, li) => ln.color || (lines.length > 1 ? th.palette[li % th.palette.length] : th.accent));
      lines.forEach((ln, li) => {
        const col = colors[li];
        let seg = [];
        const flush = () => { if (seg.length > 1) P.path(seg, col, 1.6); seg = []; };
        for (const p of ln.pts) { const X = L.X(p.x), Y = L.Y(p.y); if (isFinite(X) && isFinite(Y)) seg.push([X, Y]); else flush(); }
        flush();
        for (const p of ln.pts) {
          const X = L.X(p.x), Y = L.Y(p.y); if (!isFinite(Y)) continue;
          const sel = s.selected !== undefined && s.selected !== null && s.selected === p.id;
          P.circle(X, Y, sel ? 5 : p.mark ? 3.6 : 2.6, sel ? th.sel : p.mark ? th.ok : col);
        }
      });
      // selected marker for spectra
      if (s.sticks && s.selected !== undefined && s.selected !== null) {
        const p = s.sticks.find(q => q.id === s.selected);
        if (p) P.line(L.X(p.x), m.t, L.X(p.x), m.t + L.ph, th.sel, 1, [4, 4]);
      }
      // peak labels
      if (s.peakLabels && s.sticks && s.sticks.length) {
        const ms = Math.max(...s.sticks.map(p => Math.abs(p.y))) || 1;
        const strong = s.sticks.filter(p => Math.abs(p.y) > ms * (s.peakThreshold || 0.08)).sort((a, b) => Math.abs(b.y) - Math.abs(a.y)).slice(0, s.maxPeakLabels || 14);
        const placed = [];
        for (const p of strong) {
          const x = L.X(p.x);
          if (placed.some(q => Math.abs(q - x) < 30)) continue;
          placed.push(x);
          const neg = p.y < 0;
          let yv = p.y * k;
          if (s.series && s.series[0]) { // put the label on the curve near this x
            const ser = s.series[0]; let best = 0;
            for (let i = 0; i < ser.x.length; i++) if (Math.abs(ser.x[i] - p.x) < Math.abs((L.x1 - L.x0) / 300)) best = neg ? Math.min(best, ser.y[i]) : Math.max(best, ser.y[i]);
            yv = neg ? Math.min(yv, best) : Math.max(yv, best);
          }
          const y = L.Y(yv) + (s.invertY !== neg ? 14 : -7);
          P.text(s.peakFormat ? s.peakFormat(p.x) : Math.round(p.x).toString(), x, y, th.fg, 10, 'center');
        }
      }
      P.unclip();
      this._legend(P, L, th, lines, colors);
    }

    toSVG(w, h) {
      w = w || 900; h = h || 420;
      const P = new SvgPainter(w, h);
      P.bg(PRINT.bg);
      if (this.spec && !this.spec.empty) this._render(P, w, h, PRINT);
      return P.done();
    }
    toPNG(scale) {
      scale = scale || 2;
      const w = this.canvas.clientWidth || 900, h = this.canvas.clientHeight || 420;
      const c = document.createElement('canvas'); c.width = w * scale; c.height = h * scale;
      const ctx = c.getContext('2d'); ctx.setTransform(scale, 0, 0, scale, 0, 0);
      ctx.fillStyle = PRINT.bg; ctx.fillRect(0, 0, w, h);
      if (this.spec && !this.spec.empty) this._render(new CanvasPainter(ctx), w, h, PRINT);
      return new Promise(res => c.toBlob(res, 'image/png'));
    }
  }

  class CanvasPainter {
    constructor(ctx) { this.c = ctx; }
    line(x1, y1, x2, y2, col, w, dash) { const c = this.c; c.strokeStyle = col; c.lineWidth = w; c.setLineDash(dash || []); c.beginPath(); c.moveTo(x1, y1); c.lineTo(x2, y2); c.stroke(); c.setLineDash([]); }
    rect(x, y, w, h, col) { const c = this.c; c.strokeStyle = col; c.lineWidth = 1; c.strokeRect(x + 0.5, y + 0.5, w, h); }
    fillRect(x, y, w, h, col, alpha) { const c = this.c; c.globalAlpha = alpha === undefined ? 1 : alpha; c.fillStyle = col; c.fillRect(x, y, w, h); c.globalAlpha = 1; }
    text(t, x, y, col, size, align) { const c = this.c; c.fillStyle = col; c.font = size + 'px ui-sans-serif, -apple-system, "Segoe UI", sans-serif'; c.textAlign = align || 'left'; c.textBaseline = 'alphabetic'; c.fillText(t, x, y); }
    textV(t, x, y, col, size) { const c = this.c; c.save(); c.translate(x, y); c.rotate(-Math.PI / 2); this.text(t, 0, 0, col, size, 'center'); c.restore(); }
    path(pts, col, w) { if (!pts.length) return; const c = this.c; c.strokeStyle = col; c.lineWidth = w; c.lineJoin = 'round'; c.beginPath(); c.moveTo(pts[0][0], pts[0][1]); for (let i = 1; i < pts.length; i++) c.lineTo(pts[i][0], pts[i][1]); c.stroke(); }
    fillPath(pts, base, col) { if (!pts.length) return; const c = this.c; c.fillStyle = col; c.beginPath(); c.moveTo(pts[0][0], base); for (const p of pts) c.lineTo(p[0], p[1]); c.lineTo(pts[pts.length - 1][0], base); c.closePath(); c.fill(); }
    circle(x, y, r, col) { const c = this.c; c.fillStyle = col; c.beginPath(); c.arc(x, y, r, 0, Math.PI * 2); c.fill(); }
    clip(x, y, w, h) { const c = this.c; c.save(); c.beginPath(); c.rect(x, y, w, h); c.clip(); }
    unclip() { this.c.restore(); }
  }
  class SvgPainter {
    constructor(w, h) { this.w = w; this.h = h; this.out = []; this.clipId = 0; this.open = 0; }
    bg(col) { this.out.push(`<rect x="0" y="0" width="${this.w}" height="${this.h}" fill="${col}"/>`); }
    line(x1, y1, x2, y2, col, w, dash) { this.out.push(`<line x1="${x1.toFixed(2)}" y1="${y1.toFixed(2)}" x2="${x2.toFixed(2)}" y2="${y2.toFixed(2)}" stroke="${col}" stroke-width="${w}"${dash ? ` stroke-dasharray="${dash.join(' ')}"` : ''}/>`); }
    rect(x, y, w, h, col) { this.out.push(`<rect x="${x + 0.5}" y="${y + 0.5}" width="${w}" height="${h}" fill="none" stroke="${col}"/>`); }
    fillRect(x, y, w, h, col, alpha) { this.out.push(`<rect x="${x.toFixed(2)}" y="${y.toFixed(2)}" width="${w.toFixed(2)}" height="${h.toFixed(2)}" fill="${col}"${alpha !== undefined && alpha < 1 ? ` fill-opacity="${alpha}"` : ''}/>`); }
    text(t, x, y, col, size, align) { const a = align === 'center' ? 'middle' : align === 'right' ? 'end' : 'start'; this.out.push(`<text x="${x.toFixed(2)}" y="${y.toFixed(2)}" fill="${col}" font-size="${size}" font-family="Helvetica, Arial, sans-serif" text-anchor="${a}">${esc(t)}</text>`); }
    textV(t, x, y, col, size) { this.out.push(`<text transform="translate(${x},${y}) rotate(-90)" fill="${col}" font-size="${size}" font-family="Helvetica, Arial, sans-serif" text-anchor="middle">${esc(t)}</text>`); }
    path(pts, col, w) { if (!pts.length) return; this.out.push(`<polyline fill="none" stroke="${col}" stroke-width="${w}" stroke-linejoin="round" points="${pts.map(p => p[0].toFixed(2) + ',' + p[1].toFixed(2)).join(' ')}"/>`); }
    fillPath(pts, base, col) { if (!pts.length) return; const d = [`M${pts[0][0].toFixed(2)},${base.toFixed(2)}`].concat(pts.map(p => `L${p[0].toFixed(2)},${p[1].toFixed(2)}`), [`L${pts[pts.length - 1][0].toFixed(2)},${base.toFixed(2)}Z`]).join(''); this.out.push(`<path d="${d}" fill="${col}"/>`); }
    circle(x, y, r, col) { this.out.push(`<circle cx="${x.toFixed(2)}" cy="${y.toFixed(2)}" r="${r}" fill="${col}"/>`); }
    clip(x, y, w, h) { const id = 'c' + (this.clipId++); this.out.push(`<clipPath id="${id}"><rect x="${x}" y="${y}" width="${w}" height="${h}"/></clipPath><g clip-path="url(#${id})">`); this.open++; }
    unclip() { if (this.open) { this.out.push('</g>'); this.open--; } }
    done() { while (this.open) this.unclip(); return `<?xml version="1.0" encoding="UTF-8"?>\n<svg xmlns="http://www.w3.org/2000/svg" width="${this.w}" height="${this.h}" viewBox="0 0 ${this.w} ${this.h}">\n${this.out.join('\n')}\n</svg>\n`; }
  }

  BM.spectrum = { ir, raman, uv, ecd, broaden, lineshape, EV_TO_CM };
  BM.Plot = Plot;
})(typeof window !== 'undefined' ? (window.BM = window.BM || {}) : (globalThis.BM = globalThis.BM || {}));
