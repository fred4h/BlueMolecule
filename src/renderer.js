/* BlueMolecule — WebGL2 molecule renderer (instanced spheres and cylinders, surfaces, picking, overlay labels,
 * bond orders, tool input hooks, dashed bonds, unit cell, extra arrows and atom markers).
 * Per-atom representations: setRepresentation({style: Uint8Array (Renderer.STYLE codes), scale, alpha, colors, extent})
 * draws each atom (and each half of every bond) in its own style — ball-and-stick, sticks, wire, spacefill, screen-space
 * lines (cheap enough for 300k atoms), translucent ghosts (separate depth-pre-pass blend) or hidden; null restores the
 * single toolbar style. setOverlays([mesh]) adds vertex-coloured meshes (cartoons, rings, SNFG symbols),
 * setPickProxies(indices, r) makes hidden atoms under a cartoon pickable, setRepMarkers and setLegend add markers and a
 * legend box. Sphere and cylinder detail drop with the number of drawn instances. */
(function (BM) {
  'use strict';
  const el = BM.el;

  // ---------- math ----------
  function quatMul(a, b) {
    return [a[3] * b[0] + a[0] * b[3] + a[1] * b[2] - a[2] * b[1],
      a[3] * b[1] - a[0] * b[2] + a[1] * b[3] + a[2] * b[0],
      a[3] * b[2] + a[0] * b[1] - a[1] * b[0] + a[2] * b[3],
      a[3] * b[3] - a[0] * b[0] - a[1] * b[1] - a[2] * b[2]];
  }
  function quatAxis(axis, ang) {
    const l = Math.hypot(axis[0], axis[1], axis[2]) || 1, s = Math.sin(ang / 2) / l;
    return [axis[0] * s, axis[1] * s, axis[2] * s, Math.cos(ang / 2)];
  }
  function quatNorm(q) { const l = Math.hypot(q[0], q[1], q[2], q[3]) || 1; return [q[0] / l, q[1] / l, q[2] / l, q[3] / l]; }
  /** 3x3 rotation, row-major r[row*3+col] */
  function quatMat(q) {
    const [x, y, z, w] = q;
    return [1 - 2 * (y * y + z * z), 2 * (x * y - z * w), 2 * (x * z + y * w),
      2 * (x * y + z * w), 1 - 2 * (x * x + z * z), 2 * (y * z - x * w),
      2 * (x * z - y * w), 2 * (y * z + x * w), 1 - 2 * (x * x + y * y)];
  }
  const cross3 = (a, b) => [a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]];
  const lerp3 = (a, b, t) => [a[0] + (b[0] - a[0]) * t, a[1] + (b[1] - a[1]) * t, a[2] + (b[2] - a[2]) * t];

  /** Entry distance of the ray o + t·d (|d| = 1) into a sphere, or Infinity on a miss. */
  function raySphere(o, d, cx, cy, cz, r) {
    const px = o[0] - cx, py = o[1] - cy, pz = o[2] - cz;
    const b = px * d[0] + py * d[1] + pz * d[2], h = b * b - (px * px + py * py + pz * pz - r * r);
    return h < 0 ? Infinity : -b - Math.sqrt(h);
  }
  /** Entry distance of the ray into the capsule of radius r around segment a–b (ray–segment distance ≤ r), or Infinity. */
  function rayCapsule(o, d, ax, ay, az, bx, by, bz, r) {
    const bax = bx - ax, bay = by - ay, baz = bz - az, oax = o[0] - ax, oay = o[1] - ay, oaz = o[2] - az;
    const baba = bax * bax + bay * bay + baz * baz, bard = bax * d[0] + bay * d[1] + baz * d[2];
    const baoa = bax * oax + bay * oay + baz * oaz, rdoa = d[0] * oax + d[1] * oay + d[2] * oaz;
    const A = baba - bard * bard;
    if (A > 1e-9 * baba) {
      const B = baba * rdoa - baoa * bard, C = baba * (oax * oax + oay * oay + oaz * oaz) - baoa * baoa - r * r * baba;
      const h = B * B - A * C;
      if (h < 0) return Infinity;
      const t = (-B - Math.sqrt(h)) / A, y = baoa + t * bard;
      if (y > 0 && y < baba) return t;
    }
    return Math.min(raySphere(o, d, ax, ay, az, r), raySphere(o, d, bx, by, bz, r));
  }
  /** Nearer of an atom hit and a bond hit along the ray; atoms win ties. */
  function nearerHit(a, b) {
    if (a.index >= 0 && (b.index < 0 || a.t <= b.t + 1e-6)) return { kind: 'atom', index: a.index };
    if (b.index >= 0) return { kind: 'bond', index: b.index };
    return { kind: 'none', index: -1 };
  }
  /** Shortest path i → … → j avoiding the direct i–j edge (the smallest ring through that bond, ≤ maxRing atoms), or null. */
  function ringPath(adj, i, j, maxRing) {
    const prev = new Map([[i, -1]]);
    let front = [i];
    for (let depth = 1; depth < maxRing && front.length; depth++) {
      const next = [];
      for (const u of front) {
        for (const v of adj[u]) {
          if ((u === i && v === j) || prev.has(v)) continue;
          prev.set(v, u);
          if (v === j) { const path = [j]; for (let x = u; x !== -1; x = prev.get(x)) path.push(x); return path.reverse(); }
          next.push(v);
        }
      }
      front = next;
    }
    return null;
  }

  // ---------- cylinder instance helpers (start xyz, end xyz, radius, rgb) ----------
  function pushCyl(out, a, b, r, c) { out.push(a[0], a[1], a[2], b[0], b[1], b[2], r, c[0], c[1], c[2]); }
  /** Cylinder a→b split at its midpoint: colour ci on a's half, cj on b's half. */
  function pushSplit(out, a, b, r, ci, cj) { const m = lerp3(a, b, 0.5); pushCyl(out, a, m, r, ci); pushCyl(out, m, b, r, cj); }
  /** Evenly spaced dashes between fractions t0..t1 of a→b (duty = dash share of each period); a dash
   *  crossing the midpoint is split so colours change there (cj null = a single colour). */
  function pushDashed(out, a, b, r, ci, cj, t0, t1, period, duty) {
    const len = Math.hypot(b[0] - a[0], b[1] - a[1], b[2] - a[2]) * (t1 - t0);
    if (!(len > 1e-6)) return;
    const nd = Math.max(1, Math.round(len / period)), step = (t1 - t0) / nd, half = step * duty / 2;
    for (let q = 0; q < nd; q++) {
      const c = t0 + step * (q + 0.5), s0 = c - half, s1 = c + half;
      if (cj && s0 < 0.5 && s1 > 0.5) { pushCyl(out, lerp3(a, b, s0), lerp3(a, b, 0.5), r, ci); pushCyl(out, lerp3(a, b, 0.5), lerp3(a, b, s1), r, cj); }
      else pushCyl(out, lerp3(a, b, s0), lerp3(a, b, s1), r, cj && c > 0.5 ? cj : ci);
    }
  }

  // ---------- meshes ----------
  function icosphere(sub) {
    const t = (1 + Math.sqrt(5)) / 2;
    let v = [[-1, t, 0], [1, t, 0], [-1, -t, 0], [1, -t, 0], [0, -1, t], [0, 1, t], [0, -1, -t], [0, 1, -t], [t, 0, -1], [t, 0, 1], [-t, 0, -1], [-t, 0, 1]]
      .map(p => { const l = Math.hypot(...p); return [p[0] / l, p[1] / l, p[2] / l]; });
    let f = [[0, 11, 5], [0, 5, 1], [0, 1, 7], [0, 7, 10], [0, 10, 11], [1, 5, 9], [5, 11, 4], [11, 10, 2], [10, 7, 6], [7, 1, 8],
      [3, 9, 4], [3, 4, 2], [3, 2, 6], [3, 6, 8], [3, 8, 9], [4, 9, 5], [2, 4, 11], [6, 2, 10], [8, 6, 7], [9, 8, 1]];
    for (let s = 0; s < sub; s++) {
      const cache = new Map(), nf = [];
      const mid = (a, b) => {
        const k = a < b ? a * 100000 + b : b * 100000 + a;
        if (cache.has(k)) return cache.get(k);
        const p = [(v[a][0] + v[b][0]) / 2, (v[a][1] + v[b][1]) / 2, (v[a][2] + v[b][2]) / 2];
        const l = Math.hypot(...p); v.push([p[0] / l, p[1] / l, p[2] / l]);
        cache.set(k, v.length - 1); return v.length - 1;
      };
      for (const [a, b, c] of f) { const ab = mid(a, b), bc = mid(b, c), ca = mid(c, a); nf.push([a, ab, ca], [b, bc, ab], [c, ca, bc], [ab, bc, ca]); }
      f = nf;
    }
    return { pos: new Float32Array(v.flat()), idx: new Uint16Array(f.flat()) };
  }
  /** Unit-radius tube along +y from y=0 to y=1; `capped` closes both ends with flat discs. */
  function cylinderMesh(seg, capped) {
    const pos = [], nrm = [], idx = [];
    for (let i = 0; i <= seg; i++) {
      const a = (i / seg) * Math.PI * 2, c = Math.cos(a), s = Math.sin(a);
      pos.push(c, 0, s, c, 1, s); nrm.push(c, 0, s, c, 0, s);
    }
    for (let i = 0; i < seg; i++) { const a = 2 * i, b = a + 1, c = a + 2, d = a + 3; idx.push(a, b, c, b, d, c); }
    if (capped) {
      for (const y of [0, 1]) {
        const base = pos.length / 3, ny = y ? 1 : -1;
        pos.push(0, y, 0); nrm.push(0, ny, 0);
        for (let i = 0; i <= seg; i++) { const a = (i / seg) * Math.PI * 2; pos.push(Math.cos(a), y, Math.sin(a)); nrm.push(0, ny, 0); }
        for (let i = 0; i < seg; i++) { if (y) idx.push(base, base + 2 + i, base + 1 + i); else idx.push(base, base + 1 + i, base + 2 + i); }
      }
    }
    return { pos: new Float32Array(pos), nrm: new Float32Array(nrm), idx: new Uint16Array(idx) };
  }
  function coneMesh(seg) {
    const pos = [], nrm = [], idx = [];
    const k = 1 / Math.SQRT2;
    for (let i = 0; i <= seg; i++) {
      const a = (i / seg) * Math.PI * 2, c = Math.cos(a), s = Math.sin(a);
      pos.push(c, 0, s, 0, 1, 0); nrm.push(c * k, k, s * k, c * k, k, s * k);
    }
    for (let i = 0; i < seg; i++) { const a = 2 * i; idx.push(a, a + 1, a + 2); }
    const base = pos.length / 3;
    pos.push(0, 0, 0); nrm.push(0, -1, 0);
    for (let i = 0; i <= seg; i++) { const a = (i / seg) * Math.PI * 2; pos.push(Math.cos(a), 0, Math.sin(a)); nrm.push(0, -1, 0); }
    for (let i = 0; i < seg; i++) idx.push(base, base + 1 + i, base + 2 + i);
    return { pos: new Float32Array(pos), nrm: new Float32Array(nrm), idx: new Uint16Array(idx) };
  }

  // ---------- shaders ----------
  const FRAG = `#version 300 es
precision highp float;
in vec3 vNormal; in vec3 vEye; in vec3 vColor;
uniform vec3 uLight; uniform float uAmbient; uniform float uSpec; uniform float uShine;
uniform vec3 uFogColor; uniform vec2 uFogRange; uniform float uFogOn;
uniform float uAlpha; uniform float uOrtho; uniform float uTwoSided;
in float vAlpha;
out vec4 outColor;
void main(){
  vec3 N = normalize(vNormal);
  vec3 V = uOrtho > 0.5 ? vec3(0.0,0.0,1.0) : normalize(-vEye);
  if (uTwoSided > 0.5 && dot(N, V) < 0.0) N = -N;
  float diff = max(dot(N, uLight), 0.0);
  float fill = max(dot(N, normalize(vec3(-0.6, -0.2, 0.7))), 0.0) * 0.22;
  float rim = pow(1.0 - max(dot(N, V), 0.0), 3.0) * 0.12;
  vec3 H = normalize(uLight + V);
  float sp = pow(max(dot(N, H), 0.0), uShine) * uSpec;
  vec3 c = vColor * (uAmbient + (1.0 - uAmbient) * diff + fill) + vec3(sp) + vColor * rim;
  float depth = -vEye.z;
  float f = clamp((depth - uFogRange.x) / max(uFogRange.y - uFogRange.x, 1e-3), 0.0, 1.0) * uFogOn;
  c = mix(c, uFogColor, f * 0.7);
  float A = uAlpha * vAlpha;
  outColor = vec4(c * A, A);
}`;
  // sphere / cylinder instances: rgb colour, or rgba (translucent ghost atoms) when `alpha`
  const vsSphere = alpha => `#version 300 es
layout(location=0) in vec3 aPos;
layout(location=1) in vec4 aCenterRadius;
layout(location=2) in ${alpha ? 'vec4' : 'vec3'} aColor;
uniform mat4 uView; uniform mat4 uProj; uniform float uGrow; uniform float uGrowAbs;
out vec3 vNormal; out vec3 vEye; out vec3 vColor; out float vAlpha;
void main(){
  float r = aCenterRadius.w * uGrow + uGrowAbs;
  vec4 eye = uView * vec4(aCenterRadius.xyz + aPos * r, 1.0);
  vEye = eye.xyz; vNormal = mat3(uView) * aPos; vColor = aColor.rgb; vAlpha = ${alpha ? 'aColor.a' : '1.0'};
  gl_Position = uProj * eye;
}`;
  const VS_SPHERE = vsSphere(false);
  const vsCyl = alpha => `#version 300 es
layout(location=0) in vec3 aPos;
layout(location=1) in vec3 aNrm;
layout(location=2) in vec3 aStart;
layout(location=3) in vec4 aEndRadius;
layout(location=4) in ${alpha ? 'vec4' : 'vec3'} aColor;
uniform mat4 uView; uniform mat4 uProj;
out vec3 vNormal; out vec3 vEye; out vec3 vColor; out float vAlpha;
void main(){
  vec3 axis = aEndRadius.xyz - aStart; float len = length(axis);
  vec3 Y = axis / max(len, 1e-6);
  vec3 t = abs(Y.y) < 0.99 ? vec3(0.0, 1.0, 0.0) : vec3(1.0, 0.0, 0.0);
  vec3 X = normalize(cross(t, Y)); vec3 Zv = cross(X, Y);
  float r = aEndRadius.w;
  vec3 world = aStart + X * (aPos.x * r) + Y * (aPos.y * len) + Zv * (aPos.z * r);
  vec4 eye = uView * vec4(world, 1.0);
  vEye = eye.xyz; vNormal = mat3(uView) * normalize(X * aNrm.x + Y * aNrm.y + Zv * aNrm.z); vColor = aColor.rgb; vAlpha = ${alpha ? 'aColor.a' : '1.0'};
  gl_Position = uProj * eye;
}`;
  const VS_CYL = vsCyl(false);
  const VS_MESH = `#version 300 es
layout(location=0) in vec3 aPos;
layout(location=1) in vec3 aNrm;
uniform mat4 uView; uniform mat4 uProj; uniform vec3 uColor;
out vec3 vNormal; out vec3 vEye; out vec3 vColor; out float vAlpha;
void main(){
  vec4 eye = uView * vec4(aPos, 1.0);
  vEye = eye.xyz; vNormal = mat3(uView) * aNrm; vColor = uColor; vAlpha = 1.0;
  gl_Position = uProj * eye;
}`;
  // overlay meshes (cartoons, rings, SNFG symbols): per-vertex colours
  const VS_CMESH = `#version 300 es
layout(location=0) in vec3 aPos;
layout(location=1) in vec3 aNrm;
layout(location=2) in vec3 aColor;
uniform mat4 uView; uniform mat4 uProj;
out vec3 vNormal; out vec3 vEye; out vec3 vColor; out float vAlpha;
void main(){
  vec4 eye = uView * vec4(aPos, 1.0);
  vEye = eye.xyz; vNormal = mat3(uView) * aNrm; vColor = aColor; vAlpha = 1.0;
  gl_Position = uProj * eye;
}`;
  // lines: one instanced screen-space quad per segment (constant pixel width, depth-tested), colour split at the middle
  const VS_LINE = `#version 300 es
layout(location=0) in vec2 aCorner;
layout(location=1) in vec3 aA;
layout(location=2) in vec3 aB;
layout(location=3) in vec3 aCA;
layout(location=4) in vec3 aCB;
uniform mat4 uView; uniform mat4 uProj; uniform vec2 uViewport; uniform float uWidth; uniform float uNear;
out float vT; out float vDepth; flat out vec3 vCA; flat out vec3 vCB;
void main(){
  vec4 ea = uView * vec4(aA, 1.0), eb = uView * vec4(aB, 1.0);
  float lim = -uNear;
  vCA = aCA; vCB = aCB;
  if (ea.z > lim && eb.z > lim) { gl_Position = vec4(2.0, 2.0, 2.0, 1.0); vT = 0.0; vDepth = 0.0; return; }
  if (ea.z > lim) ea = mix(eb, ea, (lim - eb.z) / (ea.z - eb.z));
  if (eb.z > lim) eb = mix(ea, eb, (lim - ea.z) / (eb.z - ea.z));
  vec4 ca = uProj * ea, cb = uProj * eb;
  vec2 sa = ca.xy / ca.w * uViewport, sb = cb.xy / cb.w * uViewport, d = sb - sa;
  float L = length(d);
  vec2 dir = L > 1e-5 ? d / L : vec2(1.0, 0.0), nrm = vec2(-dir.y, dir.x);
  bool atEnd = aCorner.x > 0.5;
  vec4 c = atEnd ? cb : ca;
  vec2 off = nrm * (aCorner.y * uWidth * 0.5) + dir * ((atEnd ? 0.5 : -0.5) * uWidth * 0.5);
  c.xy += off / uViewport * c.w;
  gl_Position = c;
  vT = atEnd ? 1.0 : 0.0; vDepth = -(atEnd ? eb.z : ea.z);
}`;
  const FS_LINE = `#version 300 es
precision highp float;
in float vT; in float vDepth; flat in vec3 vCA; flat in vec3 vCB;
uniform vec3 uFogColor; uniform vec2 uFogRange; uniform float uFogOn;
out vec4 outColor;
void main(){
  vec3 c = vT < 0.5 ? vCA : vCB;
  float f = clamp((vDepth - uFogRange.x) / max(uFogRange.y - uFogRange.x, 1e-3), 0.0, 1.0) * uFogOn;
  outColor = vec4(mix(c, uFogColor, f * 0.7), 1.0);
}`;
  // Atom markers: camera-facing dashed rings around the atom (instanced quads, same instance layout as spheres).
  const VS_RING = `#version 300 es
layout(location=0) in vec2 aCorner;
layout(location=1) in vec4 aCenterRadius;
layout(location=2) in vec3 aColor;
uniform mat4 uView; uniform mat4 uProj;
out vec2 vUV; out vec3 vColor; out float vDepth;
void main(){
  vec4 eye = uView * vec4(aCenterRadius.xyz, 1.0);
  eye.xy += aCorner * aCenterRadius.w;
  vUV = aCorner; vColor = aColor; vDepth = -eye.z;
  gl_Position = uProj * eye;
}`;
  const FS_RING = `#version 300 es
precision highp float;
in vec2 vUV; in vec3 vColor; in float vDepth;
uniform float uInner; uniform float uDashes; uniform float uAlpha;
uniform vec3 uFogColor; uniform vec2 uFogRange; uniform float uFogOn;
out vec4 outColor;
void main(){
  float r = length(vUV);
  float aa = max(fwidth(r), 1e-4);
  float a = smoothstep(uInner - aa, uInner + aa, r) * (1.0 - smoothstep(1.0 - 2.0 * aa, 1.0, r));
  if (uDashes > 0.5) {
    float s = cos(atan(vUV.y, vUV.x) * uDashes);
    float sa = max(fwidth(s), 1e-3);
    a *= smoothstep(-0.3 - sa, -0.3 + sa, s);
  }
  a *= uAlpha;
  if (a < 0.003) discard;
  float band = clamp((r - uInner) / (1.0 - uInner), 0.0, 1.0);
  float rim = 1.0 - smoothstep(0.0, 0.32, min(band, 1.0 - band));
  vec3 c = vColor * (1.0 - 0.45 * rim);
  float f = clamp((vDepth - uFogRange.x) / max(uFogRange.y - uFogRange.x, 1e-3), 0.0, 1.0) * uFogOn;
  c = mix(c, uFogColor, f * 0.7);
  outColor = vec4(c * a, a);
}`;

  function compile(gl, vs, fs) {
    const mk = (type, src) => {
      const s = gl.createShader(type); gl.shaderSource(s, src); gl.compileShader(s);
      if (!gl.getShaderParameter(s, gl.COMPILE_STATUS)) throw new Error('Shader error: ' + gl.getShaderInfoLog(s));
      return s;
    };
    const p = gl.createProgram();
    gl.attachShader(p, mk(gl.VERTEX_SHADER, vs)); gl.attachShader(p, mk(gl.FRAGMENT_SHADER, fs));
    gl.linkProgram(p);
    if (!gl.getProgramParameter(p, gl.LINK_STATUS)) throw new Error('Program error: ' + gl.getProgramInfoLog(p));
    const u = {};
    const n = gl.getProgramParameter(p, gl.ACTIVE_UNIFORMS);
    for (let i = 0; i < n; i++) { const info = gl.getActiveUniform(p, i); u[info.name] = gl.getUniformLocation(p, info.name); }
    return { p, u };
  }

  const DEFAULTS = {
    style: 'ballstick', atomScale: 0.25, bondRadius: 0.11, spacefillScale: 1.0, showH: true, showDummies: true,
    bg: [0.063, 0.090, 0.133], perspective: true, fog: true, axes: true, fov: 30,
    selectionColor: [1.0, 0.80, 0.25], arrowColor: [0.15, 0.85, 0.45], measureColor: '#ffd04a'
  };
  const DEFAULT_DASH = [0.45, 0.75, 1.0], DEFAULT_MARKER = [0.25, 0.80, 1.0], DEFAULT_EXTRA_ARROW = [1.0, 0.58, 0.18];
  // per-atom style codes (BM.reps.CODE uses the same numbers)
  const STYLE = { hidden: 0, ballstick: 1, sticks: 2, spacefill: 3, wireframe: 4, lines: 5, ghost: 6 };
  const S_HIDDEN = 0, S_BALL = 1, S_STICK = 2, S_SPACE = 3, S_WIRE = 4, S_LINES = 5, S_GHOST = 6;
  const GHOST_ALPHA = 0.3, WIRE_AS_LINES = 60000; // wire tubes become screen-space lines above this many bonds
  /** Growable Float32Array of fixed-size records. */
  class Rec {
    constructor(stride, cap) { this.s = stride; this.a = new Float32Array(Math.max(1, cap) * stride); this.n = 0; }
    next() { if ((this.n + 1) * this.s > this.a.length) { const b = new Float32Array(this.a.length * 2); b.set(this.a); this.a = b; } return this.n++ * this.s; }
    data() { return this.a.subarray(0, this.n * this.s); }
  }

  class Renderer {
    constructor(canvas, overlay) {
      this.canvas = canvas; this.overlay = overlay;
      const gl = canvas.getContext('webgl2', { antialias: true, alpha: true, premultipliedAlpha: true, preserveDrawingBuffer: false });
      if (!gl) throw new Error('WebGL 2 is not available in this browser. Please use a recent Chrome, Edge, Firefox or Safari.');
      this.gl = gl;
      this.ctx2d = overlay.getContext('2d');
      this.opt = Object.assign({}, DEFAULTS);
      this.q = [0, 0, 0, 1]; this.center = [0, 0, 0]; this.dist = 20; this.radius = 5;
      this.Z = []; this.xyz = new Float64Array(0); this.bonds = new Uint32Array(0);
      this.atomColors = null; this.selection = []; this.arrows = null; this.surfaces = []; this.measureLines = [];
      this.labelProvider = null; this.hoverAtom = -1;
      this.onAtomClick = null; this.onAtomDblClick = null; this.onBackgroundClick = null; this.onHover = null; this.onViewChange = null;
      // bond orders, extra geometry and tool hooks
      this.bondOrders = null; this._bondsLen = 0; this._topo = null;
      this.dashedBonds = []; this.dashColor = DEFAULT_DASH; this.extraArrows = []; this.cell = null; this.cellOrigin = [0, 0, 0];
      this.markers = []; this.markerColor = DEFAULT_MARKER; this.highlightedBond = -1; this.highlightColor = null;
      this.overlayHook = null; this._tool = null; this._cursor = '';
      // per-atom representations, overlay meshes, pick proxies, representation markers, legend
      this.rep = null; this.overlays = []; this._proxy = null; this._proxySet = null; this.repMarkers = []; this.legend = null;
      this._dirty = false; this._raf = 0; this._lod = -1; this._cylLod = -1;
      this._initGL();
      this._initInput();
      this._resize();
      if (typeof ResizeObserver !== 'undefined') { this._ro = new ResizeObserver(() => { this._resize(); this.requestRender(); }); this._ro.observe(canvas.parentElement); }
    }

    _initGL() {
      const gl = this.gl;
      this.progSphere = compile(gl, VS_SPHERE, FRAG);
      this.progCyl = compile(gl, VS_CYL, FRAG);
      this.progMesh = compile(gl, VS_MESH, FRAG);
      this.progRing = compile(gl, VS_RING, FS_RING);
      this.progSphereA = compile(gl, vsSphere(true), FRAG);
      this.progCylA = compile(gl, vsCyl(true), FRAG);
      this.progCMesh = compile(gl, VS_CMESH, FRAG);
      this.progLine = compile(gl, VS_LINE, FS_LINE);
      this.sphereInst = gl.createBuffer(); this.cylInst = gl.createBuffer(); this.haloInst = gl.createBuffer();
      this.arrowCylInst = gl.createBuffer(); this.arrowConeInst = gl.createBuffer();
      this.multiInst = gl.createBuffer(); this.decoInst = gl.createBuffer(); this.xConeInst = gl.createBuffer();
      this.bondHaloInst = gl.createBuffer(); this.markerInst = gl.createBuffer();
      this.ghostSphInst = gl.createBuffer(); this.ghostCylInst = gl.createBuffer(); this.lineInst = gl.createBuffer();
      const cone = coneMesh(20), capped = cylinderMesh(16, true);
      this.lineVAO = this._makeLineVAO(this.lineInst);
      this._setCylLOD(20);
      this.arrowCylVAO = this._makeCylVAO(cylinderMesh(12), this.arrowCylInst);
      this.coneVAO = this._makeCylVAO(cone, this.arrowConeInst);
      this.multiVAO = this._makeCylVAO(capped, this.multiInst);      // thin lines of double/triple/aromatic bonds
      this.decoVAO = this._makeCylVAO(capped, this.decoInst);        // dashed bonds, cell edges, extra arrow shafts
      this.xConeVAO = this._makeCylVAO(coneMesh(24), this.xConeInst); // extra arrow heads
      this.bondHaloVAO = this._makeCylVAO(cylinderMesh(24, true), this.bondHaloInst);
      this.markerVAO = this._makeRingVAO(this.markerInst);
      this.nSphere = 0; this.nCyl = 0; this.nHalo = 0; this.nArrowCyl = 0; this.nCone = 0;
      this.nMainCyl = 0; this.nMultiCyl = 0; this.nDeco = 0; this.nXCone = 0; this.nBondHalo = 0; this.nMarker = 0;
      this.nLine = 0; this.nGhostSph = 0; this.nGhostCyl = 0;
      this._setLOD(3);
    }
    /** Sphere detail (icosphere subdivisions) for atoms, halos and ghosts. */
    _setLOD(sub) {
      if (this._lod === sub) return;
      const gl = this.gl;
      this._lod = sub;
      const s = icosphere(sub);
      for (const k of ['sphereVAO', 'haloVAO', 'ghostSphVAO']) if (this[k]) gl.deleteVertexArray(this[k]);
      this.sphereMesh = s;
      this.sphereVAO = this._makeSphereVAO(s, this.sphereInst);
      this.haloVAO = this._makeSphereVAO(s, this.haloInst);
      this.ghostSphVAO = this._makeSphereVAO(s, this.ghostSphInst, true);
    }
    /** Bond cylinder detail (sides) for the main and ghost bonds. */
    _setCylLOD(seg) {
      if (this._cylLod === seg) return;
      const gl = this.gl, m = cylinderMesh(seg);
      this._cylLod = seg;
      for (const k of ['cylVAO', 'ghostCylVAO']) if (this[k]) gl.deleteVertexArray(this[k]);
      this.cylVAO = this._makeCylVAO(m, this.cylInst);
      this.ghostCylVAO = this._makeCylVAO(m, this.ghostCylInst, true);
    }
    _makeSphereVAO(mesh, inst, alpha) {
      const gl = this.gl, vao = gl.createVertexArray(), st = alpha ? 32 : 28;
      gl.bindVertexArray(vao);
      const vb = gl.createBuffer(); gl.bindBuffer(gl.ARRAY_BUFFER, vb); gl.bufferData(gl.ARRAY_BUFFER, mesh.pos, gl.STATIC_DRAW);
      gl.enableVertexAttribArray(0); gl.vertexAttribPointer(0, 3, gl.FLOAT, false, 0, 0);
      const ib = gl.createBuffer(); gl.bindBuffer(gl.ELEMENT_ARRAY_BUFFER, ib); gl.bufferData(gl.ELEMENT_ARRAY_BUFFER, mesh.idx, gl.STATIC_DRAW);
      gl.bindBuffer(gl.ARRAY_BUFFER, inst);
      gl.enableVertexAttribArray(1); gl.vertexAttribPointer(1, 4, gl.FLOAT, false, st, 0); gl.vertexAttribDivisor(1, 1);
      gl.enableVertexAttribArray(2); gl.vertexAttribPointer(2, alpha ? 4 : 3, gl.FLOAT, false, st, 16); gl.vertexAttribDivisor(2, 1);
      gl.bindVertexArray(null);
      vao.count = mesh.idx.length;
      return vao;
    }
    _makeCylVAO(mesh, inst, alpha) {
      const gl = this.gl, vao = gl.createVertexArray(), st = alpha ? 44 : 40;
      gl.bindVertexArray(vao);
      const vb = gl.createBuffer(); gl.bindBuffer(gl.ARRAY_BUFFER, vb); gl.bufferData(gl.ARRAY_BUFFER, mesh.pos, gl.STATIC_DRAW);
      gl.enableVertexAttribArray(0); gl.vertexAttribPointer(0, 3, gl.FLOAT, false, 0, 0);
      const nb = gl.createBuffer(); gl.bindBuffer(gl.ARRAY_BUFFER, nb); gl.bufferData(gl.ARRAY_BUFFER, mesh.nrm, gl.STATIC_DRAW);
      gl.enableVertexAttribArray(1); gl.vertexAttribPointer(1, 3, gl.FLOAT, false, 0, 0);
      const ib = gl.createBuffer(); gl.bindBuffer(gl.ELEMENT_ARRAY_BUFFER, ib); gl.bufferData(gl.ELEMENT_ARRAY_BUFFER, mesh.idx, gl.STATIC_DRAW);
      gl.bindBuffer(gl.ARRAY_BUFFER, inst);
      gl.enableVertexAttribArray(2); gl.vertexAttribPointer(2, 3, gl.FLOAT, false, st, 0); gl.vertexAttribDivisor(2, 1);
      gl.enableVertexAttribArray(3); gl.vertexAttribPointer(3, 4, gl.FLOAT, false, st, 12); gl.vertexAttribDivisor(3, 1);
      gl.enableVertexAttribArray(4); gl.vertexAttribPointer(4, alpha ? 4 : 3, gl.FLOAT, false, st, 28); gl.vertexAttribDivisor(4, 1);
      gl.bindVertexArray(null);
      vao.count = mesh.idx.length;
      return vao;
    }
    _makeLineVAO(inst) {
      const gl = this.gl, vao = gl.createVertexArray();
      gl.bindVertexArray(vao);
      const vb = gl.createBuffer(); gl.bindBuffer(gl.ARRAY_BUFFER, vb); gl.bufferData(gl.ARRAY_BUFFER, new Float32Array([0, -1, 1, -1, 1, 1, 0, 1]), gl.STATIC_DRAW);
      gl.enableVertexAttribArray(0); gl.vertexAttribPointer(0, 2, gl.FLOAT, false, 0, 0);
      const ib = gl.createBuffer(); gl.bindBuffer(gl.ELEMENT_ARRAY_BUFFER, ib); gl.bufferData(gl.ELEMENT_ARRAY_BUFFER, new Uint16Array([0, 1, 2, 0, 2, 3]), gl.STATIC_DRAW);
      gl.bindBuffer(gl.ARRAY_BUFFER, inst);
      for (let k = 0; k < 4; k++) { gl.enableVertexAttribArray(1 + k); gl.vertexAttribPointer(1 + k, 3, gl.FLOAT, false, 48, 12 * k); gl.vertexAttribDivisor(1 + k, 1); }
      gl.bindVertexArray(null);
      vao.count = 6;
      return vao;
    }
    _makeRingVAO(inst) {
      const gl = this.gl, vao = gl.createVertexArray();
      gl.bindVertexArray(vao);
      const vb = gl.createBuffer(); gl.bindBuffer(gl.ARRAY_BUFFER, vb); gl.bufferData(gl.ARRAY_BUFFER, new Float32Array([-1, -1, 1, -1, 1, 1, -1, 1]), gl.STATIC_DRAW);
      gl.enableVertexAttribArray(0); gl.vertexAttribPointer(0, 2, gl.FLOAT, false, 0, 0);
      const ib = gl.createBuffer(); gl.bindBuffer(gl.ELEMENT_ARRAY_BUFFER, ib); gl.bufferData(gl.ELEMENT_ARRAY_BUFFER, new Uint16Array([0, 1, 2, 0, 2, 3]), gl.STATIC_DRAW);
      gl.bindBuffer(gl.ARRAY_BUFFER, inst);
      gl.enableVertexAttribArray(1); gl.vertexAttribPointer(1, 4, gl.FLOAT, false, 28, 0); gl.vertexAttribDivisor(1, 1);
      gl.enableVertexAttribArray(2); gl.vertexAttribPointer(2, 3, gl.FLOAT, false, 28, 16); gl.vertexAttribDivisor(2, 1);
      gl.bindVertexArray(null);
      vao.count = 6;
      return vao;
    }

    // ---------- scene ----------
    /** `orders` (optional) runs parallel to the bond pairs; when omitted, orders are kept only if `bonds` is the same array as before. */
    setScene(Z, xyz, bonds, keepView, orders) {
      this.Z = Z; this.xyz = xyz; this._assignBonds(bonds, orders);
      const n = Z.length;
      this.selection = this.selection.filter(i => i < n);
      this._rebuild();
      if (!keepView) this.fit();
      this.requestRender();
    }
    updatePositions(xyz) { this.xyz = xyz; this._rebuild(); this.requestRender(); }
    setBonds(bonds, orders) { this._assignBonds(bonds, orders); this._rebuild(); this.requestRender(); }
    _assignBonds(bonds, orders) {
      const b = bonds || new Uint32Array(0);
      const changed = b !== this.bonds || b.length !== this._bondsLen;
      this.bonds = b; this._bondsLen = b.length;
      if (orders !== undefined) this.bondOrders = orders || null;
      else if (changed) this.bondOrders = null;
      if (changed) this.highlightedBond = -1;
      if (changed || orders !== undefined) this._topo = null;
    }
    /** Bond orders parallel to the bond pairs: 1, 2, 3 or 1.5 (aromatic); null = all single. */
    setBondOrders(orders) { this.bondOrders = orders || null; this._topo = null; this._rebuild(); this.requestRender(); }
    setOptions(o) {
      const styleChanged = o.style !== undefined && o.style !== this.opt.style;
      Object.assign(this.opt, o);
      this._rebuild();
      if (styleChanged && this.opt.style === 'spacefill') this._updateRadius();
      this.requestRender();
    }
    setAtomColors(c) { this.atomColors = c; this._rebuild(); this.requestRender(); }
    setSelection(sel) { this.selection = sel.slice(); this._rebuildHalo(); this.requestRender(); }
    setArrows(vectors) { this.arrows = vectors; this._rebuildArrows(); this.requestRender(); }
    setMeasureLines(lines) { this.measureLines = lines || []; this.requestRender(); }
    /** Dashed thin cylinders between atom pairs ([[i, j], …] or flat [i0, j0, i1, j1, …]); they follow the atoms. */
    setDashedBonds(pairs, color) {
      const list = [];
      if (pairs && pairs.length) {
        if (typeof pairs[0] === 'number') { for (let k = 0; k + 1 < pairs.length; k += 2) list.push([pairs[k], pairs[k + 1]]); }
        else for (const p of pairs) if (p && p.length >= 2) list.push([p[0], p[1]]);
      }
      this.dashedBonds = list; this.dashColor = color || DEFAULT_DASH;
      this._rebuildDeco(); this.requestRender();
    }
    /** World-space arrows: [{from: [x,y,z], to: [x,y,z], color?: [r,g,b], radius?: shaft radius in Å}]. */
    setExtraArrows(list) { this.extraArrows = list ? Array.from(list).filter(Boolean) : []; this._rebuildDeco(); this.requestRender(); }
    /** Unit cell edges: vectors [[ax,ay,az],[bx,by,bz],[cx,cy,cz]] (or 9 numbers) from `origin`; null hides it. */
    setCell(vectors, origin) {
      let c = null;
      if (vectors && vectors.length) {
        const v = typeof vectors[0] === 'number' ? [vectors.slice(0, 3), vectors.slice(3, 6), vectors.slice(6, 9)] : vectors;
        c = [0, 1, 2].map(k => { const e = v[k] || [0, 0, 0]; return [+e[0] || 0, +e[1] || 0, +e[2] || 0]; });
        if (!c.some(e => Math.hypot(e[0], e[1], e[2]) > 1e-6)) c = null;
      }
      this.cell = c;
      this.cellOrigin = origin ? [+origin[0] || 0, +origin[1] || 0, +origin[2] || 0] : [0, 0, 0];
      this._rebuildDeco(); this.requestRender();
    }
    /** Dashed rings around atoms (e.g. frozen atoms), distinct from the selection halo. */
    setAtomMarkers(indices, color) {
      this.markers = indices ? Array.from(new Set(Array.from(indices))) : [];
      this.markerColor = color || DEFAULT_MARKER;
      this._rebuildMarkers(); this.requestRender();
    }
    /** Translucent halo around one bond (pair index), or -1 for none. Optional colour defaults to the selection colour. */
    highlightBond(index, color) {
      this.highlightedBond = index >= 0 ? Math.floor(index) : -1;
      this.highlightColor = color || null;
      this._rebuildBondHalo(); this.requestRender();
    }
    /** Per-atom representation: {style: Uint8Array (Renderer.STYLE codes), scale?: Float32Array (radius factors),
     *  alpha?: Float32Array (opacity < 1 = translucent), colors?: Float32Array(3n), extent?: Uint8Array (atoms that
     *  count for fitting the view, e.g. hidden atoms under a cartoon)}; null = every atom in opt.style.
     *  `defer` stores it for the next setScene without rebuilding now. */
    setRepresentation(rep, defer) { this.rep = rep || null; if (defer) return; this._rebuild(); this.requestRender(); }
    /** Vertex-coloured triangle meshes drawn with the molecule: [{positions, normals, colors, indices (Uint32Array), opacity?, twoSided?}]. */
    setOverlays(list) {
      const gl = this.gl;
      for (const s of this.overlays) { gl.deleteVertexArray(s._vao); for (const b of s._bufs) gl.deleteBuffer(b); }
      this.overlays = [];
      for (const m of list || []) {
        if (!m || !m.indices || !m.indices.length) continue;
        const vao = gl.createVertexArray(), bufs = [];
        gl.bindVertexArray(vao);
        [[m.positions, 0], [m.normals, 1], [m.colors, 2]].forEach(([data, loc]) => {
          const b = gl.createBuffer(); gl.bindBuffer(gl.ARRAY_BUFFER, b); gl.bufferData(gl.ARRAY_BUFFER, data, gl.STATIC_DRAW);
          gl.enableVertexAttribArray(loc); gl.vertexAttribPointer(loc, 3, gl.FLOAT, false, 0, 0); bufs.push(b);
        });
        const ib = gl.createBuffer(); gl.bindBuffer(gl.ELEMENT_ARRAY_BUFFER, ib); gl.bufferData(gl.ELEMENT_ARRAY_BUFFER, m.indices instanceof Uint32Array ? m.indices : Uint32Array.from(m.indices), gl.STATIC_DRAW); bufs.push(ib);
        gl.bindVertexArray(null);
        this.overlays.push({ _vao: vao, _bufs: bufs, _count: m.indices.length, opacity: m.opacity === undefined ? 1 : +m.opacity, twoSided: !!m.twoSided, vertices: m.positions.length / 3 });
      }
      this.requestRender();
    }
    /** Hidden atoms that stay pickable (and show selection halos) within `radius` Å, e.g. CA atoms under a cartoon. */
    setPickProxies(indices, radius) {
      this._proxy = indices && indices.length ? { list: Int32Array.from(indices), r: radius > 0 ? radius : 1.2 } : null;
      this._proxySet = this._proxy ? new Set(this._proxy.list) : null;
      this._rebuildPick(); this._rebuildHalo(); this.requestRender();
    }
    /** Extra marker rings from representations: [{indices, color: [r,g,b]}] (independent of setAtomMarkers). */
    setRepMarkers(list) { this.repMarkers = Array.isArray(list) ? list.filter(m => m && m.indices) : []; this._rebuildMarkers(); this.requestRender(); }
    /** Legend box in the top-right corner: {title, items: [{label, sub, glyph: style id, color: '#hex'}]} or null. */
    setLegend(legend) { this.legend = legend && legend.items && legend.items.length ? legend : null; this.requestRender(); }
    setSurfaces(list) {
      const gl = this.gl;
      for (const s of this.surfaces) { if (s._vao) { gl.deleteVertexArray(s._vao); gl.deleteBuffer(s._pb); gl.deleteBuffer(s._nb); gl.deleteBuffer(s._ib); } }
      this.surfaces = [];
      for (const s of list || []) {
        if (!s.mesh || !s.mesh.indices.length) { this.surfaces.push(Object.assign({}, s)); continue; }
        const vao = gl.createVertexArray(); gl.bindVertexArray(vao);
        const pb = gl.createBuffer(); gl.bindBuffer(gl.ARRAY_BUFFER, pb); gl.bufferData(gl.ARRAY_BUFFER, s.mesh.positions, gl.STATIC_DRAW);
        gl.enableVertexAttribArray(0); gl.vertexAttribPointer(0, 3, gl.FLOAT, false, 0, 0);
        const nb = gl.createBuffer(); gl.bindBuffer(gl.ARRAY_BUFFER, nb); gl.bufferData(gl.ARRAY_BUFFER, s.mesh.normals, gl.STATIC_DRAW);
        gl.enableVertexAttribArray(1); gl.vertexAttribPointer(1, 3, gl.FLOAT, false, 0, 0);
        const ib = gl.createBuffer(); gl.bindBuffer(gl.ELEMENT_ARRAY_BUFFER, ib); gl.bufferData(gl.ELEMENT_ARRAY_BUFFER, s.mesh.indices, gl.STATIC_DRAW);
        gl.bindVertexArray(null);
        this.surfaces.push(Object.assign({}, s, { _vao: vao, _pb: pb, _nb: nb, _ib: ib, _count: s.mesh.indices.length }));
      }
      this.requestRender();
    }

    /** Drawn at all (not hidden by showH/showDummies or a hidden representation style). */
    visible(i) {
      if (this._st && this._st.length === this.Z.length) return this._st[i] !== S_HIDDEN;
      const z = this.Z[i];
      if (z === 1 && !this.opt.showH) return false;
      if (z === 0 && !this.opt.showDummies) return false;
      return true;
    }
    /** Drawn sphere radius of atom i (0 for hidden atoms and for line-style atoms). */
    atomRadius(i) {
      if (this._rad && this._st && this._rad.length === this.Z.length) return this._rad[i];
      return this._styleRadius(STYLE[this.opt.style] === undefined ? S_BALL : STYLE[this.opt.style], i);
    }
    _styleRadius(s, i) {
      const z = this.Z[i], o = this.opt;
      if (s === S_HIDDEN) return 0;
      if (s === S_LINES) return 0;
      if (z === 0) return 0.18;
      switch (s) {
        case S_SPACE: return el.vdw[z] * o.spacefillScale;
        case S_STICK: return o.bondRadius * 1.45;
        case S_WIRE: return this._bonded && this._bonded[i] ? 0.045 : 0.16;
        default: return el.vdw[z] * o.atomScale; // ball-and-stick, ghost
      }
    }
    atomColor(i) {
      const r = this.rep, n = this.Z.length;
      if (r && r.colors && r.colors.length === 3 * n) return [r.colors[3 * i], r.colors[3 * i + 1], r.colors[3 * i + 2]];
      if (this.atomColors) return [this.atomColors[3 * i], this.atomColors[3 * i + 1], this.atomColors[3 * i + 2]];
      const z = this.Z[i];
      return z === 0 ? [1.0, 0.42, 0.78] : el.color[z] || [1, 0.08, 0.6];
    }
    /** Colours of every atom (representation colours, setAtomColors or element colours) as Float32Array(3n). */
    _colorArray() {
      const n = this.Z.length, r = this.rep;
      if (r && r.colors && r.colors.length === 3 * n) return r.colors;
      if (this.atomColors && this.atomColors.length >= 3 * n) return this.atomColors;
      if (this._elemCol && this._elemColZ === this.Z && this._elemCol.length === 3 * n) return this._elemCol;
      const c = new Float32Array(3 * n);
      for (let i = 0; i < n; i++) { const z = this.Z[i], e = z === 0 ? [1.0, 0.42, 0.78] : el.color[z] || [1, 0.08, 0.6]; c[3 * i] = e[0]; c[3 * i + 1] = e[1]; c[3 * i + 2] = e[2]; }
      this._elemCol = c; this._elemColZ = this.Z;
      return c;
    }

    // ---------- bond orders ----------
    /** 1 single, 2 double, 3 triple, 1.5 aromatic (values are rounded to the nearest of these). */
    _orderClass(k) {
      const o = this.bondOrders;
      if (!o || k >= o.length) return 1;
      const v = +o[k];
      return v >= 2.5 ? 3 : v >= 1.75 ? 2 : v >= 1.25 ? 1.5 : 1;
    }
    /** Radii (r*) and axis offsets (d*) of the thin lines; null for styles that always draw single cylinders. */
    _multiGeom(style) {
      const o = this.opt;
      style = style || o.style;
      if (style === 'ballstick') { const r = o.bondRadius; return { r, r2: 0.6 * r, d2: 0.95 * r, r3: 0.46 * r, d3: 1.4 * r, ra: 0.45 * r, da: 2.0 * r, dash: Math.max(0.12, 1.8 * r) }; }
      if (style === 'wireframe') { const r = 0.045; return { r, r2: 0.75 * r, d2: 1.6 * r, r3: 0.62 * r, d3: 2.3 * r, ra: 0.7 * r, da: 2.6 * r, dash: 0.16 }; }
      return null;
    }
    /** Drawn outer radius of bond k in the current style (0 in spacefill). */
    _bondExtent(k, g) {
      const o = this.opt;
      if (this.rep && this._st && this._st.length === this.Z.length) {
        const i = this.bonds[2 * k], j = this.bonds[2 * k + 1], st = this._st, sc = this.rep.scale && this.rep.scale.length === this.Z.length ? this.rep.scale : null;
        const half = a => { const s = st[a], f = sc ? sc[a] : 1; return s === S_STICK ? o.bondRadius * 1.45 * f : s === S_BALL || s === S_GHOST ? o.bondRadius * f : s === S_WIRE ? 0.045 * f : 0; };
        let e = Math.max(half(i), half(j));
        if (this.bondOrders && st[i] === st[j] && (st[i] === S_BALL || st[i] === S_WIRE)) {
          const gg = this._multiGeom(st[i] === S_BALL ? 'ballstick' : 'wireframe'), cls = this._orderClass(k);
          if (cls === 2) e = gg.d2 + gg.r2; else if (cls === 3) e = gg.d3 + gg.r3; else if (cls === 1.5) e = gg.da + gg.ra;
        }
        return e;
      }
      if (o.style === 'spacefill') return 0;
      if (o.style === 'sticks') return o.bondRadius * 1.45;
      g = g || this._multiGeom();
      const cls = this.bondOrders ? this._orderClass(k) : 1;
      return cls === 2 ? g.d2 + g.r2 : cls === 3 ? g.d3 + g.r3 : cls === 1.5 ? g.da + g.ra : g.r;
    }
    /** Per multiple bond: candidate neighbours [atom, neighbour, …] that define its plane, ring neighbours first. Cached per topology. */
    _topology() {
      const bonds = this.bonds, orders = this.bondOrders, Z = this.Z, n = Z.length, c = this._topo;
      if (c && c.bonds === bonds && c.orders === orders && c.Z === Z && c.n === n && c.len === bonds.length) return c;
      const adj = Array.from({ length: n }, () => []);
      for (let b = 0; b + 1 < bonds.length; b += 2) {
        const i = bonds[b], j = bonds[b + 1];
        if (i < n && j < n && i !== j) { adj[i].push(j); adj[j].push(i); }
      }
      const heavyFirst = list => list.slice().sort((p, q) => (Z[p] > 1 ? 0 : 1) - (Z[q] > 1 ? 0 : 1));
      const cand = new Map();
      for (let k = 0, np = bonds.length >> 1; k < np; k++) {
        if (this._orderClass(k) === 1) continue;
        const i = bonds[2 * k], j = bonds[2 * k + 1];
        if (i >= n || j >= n || i === j) continue;
        const list = [];
        const ring = ringPath(adj, i, j, 8);
        if (ring) list.push(i, ring[1], j, ring[ring.length - 2]);
        for (const nb of heavyFirst(adj[i])) if (nb !== j) list.push(i, nb);
        for (const nb of heavyFirst(adj[j])) if (nb !== i) list.push(j, nb);
        cand.set(k, list);
      }
      this._topo = { bonds, orders, Z, n, len: bonds.length, cand };
      return this._topo;
    }
    /** Unit offset ⟂ bond axis u, in the plane of a neighbouring atom (pointing towards it, i.e. into a ring);
     *  without a usable neighbour, ⟂ to the current view direction. V = view rotation (row-major). */
    _offsetDir(topo, k, u, V) {
      const xyz = this.xyz, cand = topo.cand.get(k);
      if (cand) {
        for (let c = 0; c < cand.length; c += 2) {
          const a = cand[c], nb = cand[c + 1];
          const vx = xyz[3 * nb] - xyz[3 * a], vy = xyz[3 * nb + 1] - xyz[3 * a + 1], vz = xyz[3 * nb + 2] - xyz[3 * a + 2];
          const d = vx * u[0] + vy * u[1] + vz * u[2];
          const px = vx - d * u[0], py = vy - d * u[1], pz = vz - d * u[2], pl = Math.hypot(px, py, pz);
          if (pl > 1e-4 && pl > 0.17 * Math.hypot(vx, vy, vz)) return [px / pl, py / pl, pz / pl];
        }
      }
      let w = cross3(u, [V[6], V[7], V[8]]);
      if (Math.hypot(w[0], w[1], w[2]) < 1e-3) w = cross3(u, [V[3], V[4], V[5]]);
      const l = Math.hypot(w[0], w[1], w[2]) || 1;
      return [w[0] / l, w[1] / l, w[2] / l];
    }
    /** Thin parallel lines of a double/triple bond, or the dashed inner line of an aromatic bond. */
    _pushMultiBond(out, k, i, j, cls, g, topo, V, ci, cj) {
      const xyz = this.xyz;
      const a = [xyz[3 * i], xyz[3 * i + 1], xyz[3 * i + 2]], b = [xyz[3 * j], xyz[3 * j + 1], xyz[3 * j + 2]];
      const L = Math.hypot(b[0] - a[0], b[1] - a[1], b[2] - a[2]);
      if (!(L > 1e-6)) return;
      const u = [(b[0] - a[0]) / L, (b[1] - a[1]) / L, (b[2] - a[2]) / L];
      const w = this._offsetDir(topo, k, u, V);
      const shift = (p, s) => [p[0] + w[0] * s, p[1] + w[1] * s, p[2] + w[2] * s];
      if (cls === 2) { for (const s of [-g.d2, g.d2]) pushSplit(out, shift(a, s), shift(b, s), g.r2, ci, cj); }
      else if (cls === 3) { for (const s of [-g.d3, 0, g.d3]) pushSplit(out, shift(a, s), shift(b, s), g.r3, ci, cj); }
      else {
        // aromatic: dashes beside the normal cylinder on the ring side, spread over the gap between the atom spheres
        const hidden = R => Math.sqrt(Math.max(0, R * R - g.da * g.da)) / L;
        const t0 = Math.min(0.45, Math.max(0.15, hidden(this.atomRadius(i)))), t1 = Math.max(0.55, Math.min(0.85, 1 - hidden(this.atomRadius(j))));
        pushDashed(out, shift(a, g.da), shift(b, g.da), g.ra, ci, cj, t0, t1, g.dash, 0.58);
      }
    }

    _rebuild() {
      const gl = this.gl, n = this.Z.length, xyz = this.xyz, o = this.opt, bonds = this.bonds, Z = this.Z;
      const rep = this.rep && this.rep.style && this.rep.style.length === n ? this.rep : null;
      const gcode = STYLE[o.style] === undefined ? S_BALL : STYLE[o.style];
      // per-atom style: hidden hydrogens / dummies first, then the representation (or the toolbar style)
      const st = this._st && this._st.length === n ? this._st : (this._st = new Uint8Array(n));
      for (let i = 0; i < n; i++) { const z = Z[i]; st[i] = (z === 1 && !o.showH) || (z === 0 && !o.showDummies) ? S_HIDDEN : rep ? rep.style[i] : gcode; }
      const scale = rep && rep.scale && rep.scale.length === n ? rep.scale : null;
      const alpha = rep && rep.alpha && rep.alpha.length === n ? rep.alpha : null;
      const bonded = new Uint8Array(n), visB = [];
      for (let b = 0; b + 1 < bonds.length; b += 2) {
        const i = bonds[b], j = bonds[b + 1];
        if (i >= n || j >= n || !st[i] || !st[j]) continue;
        bonded[i] = 1; bonded[j] = 1; visB.push(b >> 1);
      }
      this._bonded = bonded;
      // huge systems: thin wire tubes become screen-space lines
      this._wireLines = visB.length > WIRE_AS_LINES;
      if (this._wireLines) for (let i = 0; i < n; i++) if (st[i] === S_WIRE) st[i] = S_LINES;
      const C = this._colorArray();
      const rad = this._rad && this._rad.length === n ? this._rad : (this._rad = new Float32Array(n));
      for (let i = 0; i < n; i++) rad[i] = this._styleRadius(st[i], i) * (scale ? scale[i] : 1);
      const alphaOf = i => (st[i] === S_GHOST ? (alpha && alpha[i] < 0.999 ? alpha[i] : GHOST_ALPHA) : alpha && alpha[i] < 0.999 ? alpha[i] : 1);
      const bg = o.bg; // translucent lines fade towards the background
      let nL = 0;
      for (let i = 0; i < n; i++) if (st[i] === S_LINES) nL++;
      // spheres (opaque and translucent) and lines
      const sph = new Rec(7, n), gsp = new Rec(8, 16), lin = new Rec(12, nL ? Math.min(2 * visB.length, 2 * nL + 64) + 3 * nL : 16);
      for (let i = 0; i < n; i++) {
        const s = st[i];
        if (!s) continue;
        const x = xyz[3 * i], y = xyz[3 * i + 1], z = xyz[3 * i + 2];
        if (s === S_LINES) {
          if (bonded[i]) continue;
          const a = alphaOf(i), r = 0.28 * (scale ? scale[i] : 1), c0 = C[3 * i], c1 = C[3 * i + 1], c2 = C[3 * i + 2];
          for (let ax = 0; ax < 3; ax++) { // a small 3D cross marks an atom without bonds
            const k = lin.next(), A = lin.a;
            A[k] = x - (ax === 0 ? r : 0); A[k + 1] = y - (ax === 1 ? r : 0); A[k + 2] = z - (ax === 2 ? r : 0);
            A[k + 3] = x + (ax === 0 ? r : 0); A[k + 4] = y + (ax === 1 ? r : 0); A[k + 5] = z + (ax === 2 ? r : 0);
            A[k + 6] = A[k + 9] = a < 1 ? bg[0] + (c0 - bg[0]) * a : c0; A[k + 7] = A[k + 10] = a < 1 ? bg[1] + (c1 - bg[1]) * a : c1; A[k + 8] = A[k + 11] = a < 1 ? bg[2] + (c2 - bg[2]) * a : c2;
          }
          continue;
        }
        if (!(rad[i] > 0)) continue;
        const a = alphaOf(i);
        if (a < 1) { const k = gsp.next(), A = gsp.a; A[k] = x; A[k + 1] = y; A[k + 2] = z; A[k + 3] = rad[i]; A[k + 4] = C[3 * i]; A[k + 5] = C[3 * i + 1]; A[k + 6] = C[3 * i + 2]; A[k + 7] = a; }
        else { const k = sph.next(), A = sph.a; A[k] = x; A[k + 1] = y; A[k + 2] = z; A[k + 3] = rad[i]; A[k + 4] = C[3 * i]; A[k + 5] = C[3 * i + 1]; A[k + 6] = C[3 * i + 2]; }
      }
      // bonds: each half in the style of the atom on that side; double/triple/aromatic lines where both halves are
      // opaque, unscaled ball-and-stick (or wire)
      const cyl = new Rec(10, 2 * visB.length), gcy = new Rec(11, 16), mul = [];
      const gB = this.bondOrders ? this._multiGeom('ballstick') : null, gW = this.bondOrders ? this._multiGeom('wireframe') : null;
      let topo = null, V = null;
      const segCyl = (ax, ay, az, bx, by, bz, r, i, a) => {
        if (a < 1) { const k = gcy.next(), A = gcy.a; A[k] = ax; A[k + 1] = ay; A[k + 2] = az; A[k + 3] = bx; A[k + 4] = by; A[k + 5] = bz; A[k + 6] = r; A[k + 7] = C[3 * i]; A[k + 8] = C[3 * i + 1]; A[k + 9] = C[3 * i + 2]; A[k + 10] = a; }
        else { const k = cyl.next(), A = cyl.a; A[k] = ax; A[k + 1] = ay; A[k + 2] = az; A[k + 3] = bx; A[k + 4] = by; A[k + 5] = bz; A[k + 6] = r; A[k + 7] = C[3 * i]; A[k + 8] = C[3 * i + 1]; A[k + 9] = C[3 * i + 2]; }
      };
      const segLine = (ax, ay, az, bx, by, bz, i, ai, j, aj) => {
        const k = lin.next(), A = lin.a;
        A[k] = ax; A[k + 1] = ay; A[k + 2] = az; A[k + 3] = bx; A[k + 4] = by; A[k + 5] = bz;
        A[k + 6] = ai < 1 ? bg[0] + (C[3 * i] - bg[0]) * ai : C[3 * i]; A[k + 7] = ai < 1 ? bg[1] + (C[3 * i + 1] - bg[1]) * ai : C[3 * i + 1]; A[k + 8] = ai < 1 ? bg[2] + (C[3 * i + 2] - bg[2]) * ai : C[3 * i + 2];
        A[k + 9] = aj < 1 ? bg[0] + (C[3 * j] - bg[0]) * aj : C[3 * j]; A[k + 10] = aj < 1 ? bg[1] + (C[3 * j + 1] - bg[1]) * aj : C[3 * j + 1]; A[k + 11] = aj < 1 ? bg[2] + (C[3 * j + 2] - bg[2]) * aj : C[3 * j + 2];
      };
      const halfR = (s, f) => (s === S_STICK ? o.bondRadius * 1.45 : s === S_WIRE ? 0.045 : s === S_BALL || s === S_GHOST ? o.bondRadius : 0) * f;
      for (let t = 0; t < visB.length; t++) {
        const kb = visB[t], i = bonds[2 * kb], j = bonds[2 * kb + 1], si = st[i], sj = st[j];
        const ix = xyz[3 * i], iy = xyz[3 * i + 1], iz = xyz[3 * i + 2], jx = xyz[3 * j], jy = xyz[3 * j + 1], jz = xyz[3 * j + 2];
        const mx = (ix + jx) / 2, my = (iy + jy) / 2, mz = (iz + jz) / 2;
        const ai = alphaOf(i), aj = alphaOf(j), fi = scale ? scale[i] : 1, fj = scale ? scale[j] : 1;
        if (si === S_LINES && sj === S_LINES) { segLine(ix, iy, iz, jx, jy, jz, i, ai, j, aj); continue; }
        if (si === sj && (si === S_BALL || si === S_WIRE) && ai === 1 && aj === 1 && fi === 1 && fj === 1) {
          const g = si === S_BALL ? gB : gW, cls = g ? this._orderClass(kb) : 1;
          if (cls === 1 || cls === 1.5) { const r = halfR(si, 1); segCyl(ix, iy, iz, mx, my, mz, r, i, 1); segCyl(mx, my, mz, jx, jy, jz, r, j, 1); }
          if (cls !== 1) {
            if (!topo) { topo = this._topology(); V = quatMat(this.q); }
            this._pushMultiBond(mul, kb, i, j, cls, g, topo, V, [C[3 * i], C[3 * i + 1], C[3 * i + 2]], [C[3 * j], C[3 * j + 1], C[3 * j + 2]]);
          }
          continue;
        }
        if (si === S_LINES) segLine(ix, iy, iz, mx, my, mz, i, ai, i, ai);
        else { const r = halfR(si, fi); if (r > 0) segCyl(ix, iy, iz, mx, my, mz, r, i, ai); }
        if (sj === S_LINES) segLine(mx, my, mz, jx, jy, jz, j, aj, j, aj);
        else { const r = halfR(sj, fj); if (r > 0) segCyl(mx, my, mz, jx, jy, jz, r, j, aj); }
      }
      const up = (buf, data) => { gl.bindBuffer(gl.ARRAY_BUFFER, buf); gl.bufferData(gl.ARRAY_BUFFER, data, gl.DYNAMIC_DRAW); };
      up(this.sphereInst, sph.data()); up(this.ghostSphInst, gsp.data()); up(this.cylInst, cyl.data()); up(this.ghostCylInst, gcy.data());
      up(this.multiInst, new Float32Array(mul)); up(this.lineInst, lin.data());
      this.nSphere = sph.n; this.nGhostSph = gsp.n; this.nMainCyl = cyl.n; this.nGhostCyl = gcy.n; this.nMultiCyl = mul.length / 10; this.nLine = lin.n;
      this.nCyl = this.nMainCyl + this.nMultiCyl; // all opaque bond cylinder instances
      // detail follows the number of drawn instances
      const nsAll = sph.n + gsp.n, ncAll = cyl.n + gcy.n;
      this._setLOD(nsAll > 150000 ? 0 : nsAll > 20000 ? 1 : nsAll > 3000 ? 2 : 3);
      this._setCylLOD(ncAll > 200000 ? 6 : ncAll > 40000 ? 10 : 20);
      this._ext = rep && rep.extent && rep.extent.length === n ? rep.extent : null;
      this._rebuildPick();
      this._rebuildHalo();
      this._rebuildArrows();
      this._rebuildDeco();
      this._rebuildMarkers();
      this._rebuildBondHalo();
    }
    /** Pick radius per atom: drawn atoms max(radius, 0.28 Å), proxies their proxy radius, others 0. */
    _rebuildPick() {
      const n = this.Z.length, st = this._st, rad = this._rad;
      const pr = this._pickR && this._pickR.length === n ? this._pickR : (this._pickR = new Float32Array(n));
      if (!st || st.length !== n) { pr.fill(0); return; }
      for (let i = 0; i < n; i++) pr[i] = st[i] ? Math.max(rad[i], 0.28) : 0;
      if (this._proxy) for (const i of this._proxy.list) if (i >= 0 && i < n && !st[i]) pr[i] = this._proxy.r;
    }
    _rebuildHalo() {
      const gl = this.gl, xyz = this.xyz, px = this._proxySet;
      const sel = this.selection.filter(i => i < this.Z.length && (this.visible(i) || (px && px.has(i))));
      const c = this.opt.selectionColor;
      const buf = new Float32Array(sel.length * 7);
      sel.forEach((i, t) => buf.set([xyz[3 * i], xyz[3 * i + 1], xyz[3 * i + 2], this.visible(i) ? this.atomRadius(i) : this._proxy.r * 0.45, c[0], c[1], c[2]], t * 7));
      gl.bindBuffer(gl.ARRAY_BUFFER, this.haloInst); gl.bufferData(gl.ARRAY_BUFFER, buf, gl.DYNAMIC_DRAW);
      this.nHalo = sel.length;
    }
    _rebuildArrows() {
      const gl = this.gl, xyz = this.xyz, v = this.arrows, n = this.Z.length;
      if (!v || v.length !== 3 * n) { this.nArrowCyl = 0; this.nCone = 0; return; }
      const c = this.opt.arrowColor;
      const cyl = [], cone = [];
      for (let i = 0; i < n; i++) {
        if (!this.visible(i)) continue;
        const dx = v[3 * i], dy = v[3 * i + 1], dz = v[3 * i + 2], L = Math.hypot(dx, dy, dz);
        if (L < 0.03) continue;
        const ux = dx / L, uy = dy / L, uz = dz / L;
        const r0 = this.atomRadius(i) * 0.9;
        const head = Math.min(0.32, L * 0.45), shaftR = 0.045, headR = 0.11;
        const sx = xyz[3 * i] + ux * r0, sy = xyz[3 * i + 1] + uy * r0, sz = xyz[3 * i + 2] + uz * r0;
        const tx = sx + dx, ty = sy + dy, tz = sz + dz;
        const hx = tx - ux * head, hy = ty - uy * head, hz = tz - uz * head;
        cyl.push(sx, sy, sz, hx, hy, hz, shaftR, c[0], c[1], c[2]);
        cone.push(hx, hy, hz, tx, ty, tz, headR, c[0], c[1], c[2]);
      }
      gl.bindBuffer(gl.ARRAY_BUFFER, this.arrowCylInst); gl.bufferData(gl.ARRAY_BUFFER, new Float32Array(cyl), gl.DYNAMIC_DRAW);
      gl.bindBuffer(gl.ARRAY_BUFFER, this.arrowConeInst); gl.bufferData(gl.ARRAY_BUFFER, new Float32Array(cone), gl.DYNAMIC_DRAW);
      this.nArrowCyl = cyl.length / 10; this.nCone = cone.length / 10;
    }
    /** Dashed bonds, unit cell edges and extra arrows (capped cylinders + cones). */
    _rebuildDeco() {
      const gl = this.gl, n = this.Z.length, xyz = this.xyz, o = this.opt;
      const out = [], cones = [];
      const P = i => [xyz[3 * i], xyz[3 * i + 1], xyz[3 * i + 2]];
      if (this.dashedBonds.length) {
        const r = o.style === 'wireframe' ? 0.022 : Math.min(0.05, Math.max(0.025, o.bondRadius * 0.32));
        for (const [i, j] of this.dashedBonds) {
          if (!(i >= 0 && j >= 0 && i < n && j < n) || i === j || !this.visible(i) || !this.visible(j)) continue;
          const a = P(i), b = P(j), L = Math.hypot(b[0] - a[0], b[1] - a[1], b[2] - a[2]);
          if (!(L > 1e-6)) continue;
          const t0 = this.atomRadius(i) / L, t1 = 1 - this.atomRadius(j) / L; // from sphere surface to sphere surface
          if ((t1 - t0) * L < 0.08) continue;
          pushDashed(out, a, b, r, this.dashColor, null, t0, t1, 0.25, 0.55);
        }
      }
      if (this.cell) {
        const vec = this.cell, O = this.cellOrigin, lens = vec.map(v => Math.hypot(v[0], v[1], v[2]));
        const r = Math.min(0.1, Math.max(0.02, Math.max(...lens) * 0.003));
        const grey = this._isLightBg() ? [0.36, 0.39, 0.45] : [0.70, 0.74, 0.80];
        const axisCol = [[0.92, 0.32, 0.32], [0.32, 0.80, 0.38], [0.36, 0.52, 1.0]]; // a, b, c edges at the origin
        for (let ax = 0; ax < 3; ax++) {
          if (lens[ax] < 1e-6) continue;
          const e = vec[ax], x = [e[0] / lens[ax] * r, e[1] / lens[ax] * r, e[2] / lens[ax] * r]; // overlap at corners
          const others = [0, 1, 2].filter(q => q !== ax), p1 = vec[others[0]], p2 = vec[others[1]];
          for (let f = 0; f < 4; f++) {
            const fa = f & 1, fb = f >> 1;
            if ((fa && lens[others[0]] < 1e-6) || (fb && lens[others[1]] < 1e-6)) continue;
            const s = [O[0] + fa * p1[0] + fb * p2[0], O[1] + fa * p1[1] + fb * p2[1], O[2] + fa * p1[2] + fb * p2[2]];
            pushCyl(out, [s[0] - x[0], s[1] - x[1], s[2] - x[2]], [s[0] + e[0] + x[0], s[1] + e[1] + x[1], s[2] + e[2] + x[2]], r, f === 0 ? axisCol[ax] : grey);
          }
        }
      }
      for (const A of this.extraArrows) {
        if (!A.from || !A.to) continue;
        const f = A.from, t = A.to, dx = t[0] - f[0], dy = t[1] - f[1], dz = t[2] - f[2], L = Math.hypot(dx, dy, dz);
        if (!(L > 1e-6)) continue;
        const sr = A.radius > 0 ? A.radius : 0.06, hr = Math.max(sr * 2.3, sr + 0.05), hl = Math.min(L * 0.45, hr * 2.6);
        const c = A.color || DEFAULT_EXTRA_ARROW;
        const hb = [t[0] - dx / L * hl, t[1] - dy / L * hl, t[2] - dz / L * hl];
        pushCyl(out, [f[0], f[1], f[2]], hb, sr, c); pushCyl(cones, hb, [t[0], t[1], t[2]], hr, c);
      }
      gl.bindBuffer(gl.ARRAY_BUFFER, this.decoInst); gl.bufferData(gl.ARRAY_BUFFER, new Float32Array(out), gl.DYNAMIC_DRAW);
      gl.bindBuffer(gl.ARRAY_BUFFER, this.xConeInst); gl.bufferData(gl.ARRAY_BUFFER, new Float32Array(cones), gl.DYNAMIC_DRAW);
      this.nDeco = out.length / 10; this.nXCone = cones.length / 10;
    }
    _rebuildMarkers() {
      const gl = this.gl, xyz = this.xyz, n = this.Z.length, c = this.markerColor, out = [];
      for (const i of this.markers) if (i >= 0 && i < n && this.visible(i)) out.push(xyz[3 * i], xyz[3 * i + 1], xyz[3 * i + 2], Math.max(0.3, this.atomRadius(i) * 1.22 + 0.1), c[0], c[1], c[2]);
      for (const m of this.repMarkers || []) {
        const mc = m.color || DEFAULT_MARKER;
        for (const i of m.indices) if (i >= 0 && i < n) out.push(xyz[3 * i], xyz[3 * i + 1], xyz[3 * i + 2], Math.max(0.5, this.atomRadius(i) * 1.3 + 0.15), mc[0], mc[1], mc[2]);
      }
      gl.bindBuffer(gl.ARRAY_BUFFER, this.markerInst); gl.bufferData(gl.ARRAY_BUFFER, new Float32Array(out), gl.DYNAMIC_DRAW);
      this.nMarker = out.length / 7;
    }
    _rebuildBondHalo() {
      const gl = this.gl, k = this.highlightedBond, b = this.bonds, n = this.Z.length, xyz = this.xyz;
      this.nBondHalo = 0;
      if (!(k >= 0) || 2 * k + 1 >= b.length) return;
      const i = b[2 * k], j = b[2 * k + 1];
      if (i >= n || j >= n || !this.visible(i) || !this.visible(j)) return;
      const r = this._bondExtent(k) * 1.35 + 0.07, c = this.highlightColor || this.opt.selectionColor;
      gl.bindBuffer(gl.ARRAY_BUFFER, this.bondHaloInst);
      gl.bufferData(gl.ARRAY_BUFFER, new Float32Array([xyz[3 * i], xyz[3 * i + 1], xyz[3 * i + 2], xyz[3 * j], xyz[3 * j + 1], xyz[3 * j + 2], r, c[0], c[1], c[2]]), gl.DYNAMIC_DRAW);
      this.nBondHalo = 1;
    }

    // ---------- camera ----------
    _updateRadius() {
      const ext = this._ext && this._ext.length === this.Z.length ? this._ext : null;
      const b = BM.geom.bounds(this.xyz, this.Z, ext ? (i => ext[i] === 1 || this.visible(i)) : (i => this.visible(i)));
      this.radius = b.radius;
      return b;
    }
    fit() {
      const b = this._updateRadius();
      this.center = b.center;
      const half = (this.opt.fov / 2) * Math.PI / 180;
      const aspect = Math.max(0.2, this.canvas.clientWidth / Math.max(1, this.canvas.clientHeight));
      const fitR = this.radius * 1.05;
      this.dist = Math.max(fitR / Math.sin(half), fitR / Math.sin(Math.atan(Math.tan(half) * aspect)));
      this.requestRender();
    }
    resetView() { this.q = [0, 0, 0, 1]; this.fit(); }
    setOrientation(q) { this.q = quatNorm(q); this.requestRender(); }
    viewAlong(axis) {
      // look down +axis: rotate so that axis maps to +z (towards viewer)
      const map = { x: quatAxis([0, 1, 0], -Math.PI / 2), y: quatAxis([1, 0, 0], Math.PI / 2), z: [0, 0, 0, 1] };
      this.q = map[axis] || [0, 0, 0, 1]; this.requestRender();
    }
    rotateBy(axisView, ang) { this.q = quatNorm(quatMul(quatAxis(axisView, ang), this.q)); this.requestRender(); }
    centerOn(i) { this.center = [this.xyz[3 * i], this.xyz[3 * i + 1], this.xyz[3 * i + 2]]; this.requestRender(); }

    _matrices(w, h) {
      const R = quatMat(this.q), c = this.center, d = this.dist;
      // column-major view matrix
      const tx = -(R[0] * c[0] + R[1] * c[1] + R[2] * c[2]);
      const ty = -(R[3] * c[0] + R[4] * c[1] + R[5] * c[2]);
      const tz = -(R[6] * c[0] + R[7] * c[1] + R[8] * c[2]) - d;
      const view = new Float32Array([R[0], R[3], R[6], 0, R[1], R[4], R[7], 0, R[2], R[5], R[8], 0, tx, ty, tz, 1]);
      const aspect = w / Math.max(1, h);
      const half = (this.opt.fov / 2) * Math.PI / 180;
      const near = Math.max(0.05, d - this.radius * 3), far = d + this.radius * 3 + 10;
      let proj;
      if (this.opt.perspective) {
        const f = 1 / Math.tan(half);
        proj = new Float32Array([f / aspect, 0, 0, 0, 0, f, 0, 0, 0, 0, (far + near) / (near - far), -1, 0, 0, (2 * far * near) / (near - far), 0]);
      } else {
        const hh = d * Math.tan(half), hw = hh * aspect, n2 = -this.radius * 4 - 10, f2 = d + this.radius * 4 + 10;
        proj = new Float32Array([1 / hw, 0, 0, 0, 0, 1 / hh, 0, 0, 0, 0, -2 / (f2 - n2), 0, 0, 0, -(f2 + n2) / (f2 - n2), 1]);
      }
      return { view, proj, R, aspect, half, near, far };
    }
    /** world -> CSS pixel coords {x, y, depth} */
    project(p, m, w, h) {
      const R = m.R, c = this.center;
      const vx = R[0] * (p[0] - c[0]) + R[1] * (p[1] - c[1]) + R[2] * (p[2] - c[2]);
      const vy = R[3] * (p[0] - c[0]) + R[4] * (p[1] - c[1]) + R[5] * (p[2] - c[2]);
      const vz = R[6] * (p[0] - c[0]) + R[7] * (p[1] - c[1]) + R[8] * (p[2] - c[2]) - this.dist;
      let nx, ny;
      if (this.opt.perspective) { const s = -vz * Math.tan(m.half); nx = vx / (s * m.aspect); ny = vy / s; }
      else { const hh = this.dist * Math.tan(m.half); nx = vx / (hh * m.aspect); ny = vy / hh; }
      return { x: (nx + 1) / 2 * w, y: (1 - ny) / 2 * h, z: vz };
    }
    /** World point -> CSS pixels relative to the canvas: {x, y, z}, z = eye-space z (negative in front of the camera). */
    worldToScreen(p) {
      const w = this.canvas.clientWidth, h = this.canvas.clientHeight;
      return this.project(p, this._matrices(w, h), w, h);
    }
    /** CSS pixel -> world point on the plane parallel to the screen through `throughPoint` (default: the view centre). */
    screenToWorld(px, py, throughPoint) {
      const { o, d } = this.ray(px, py), R = quatMat(this.q), c = throughPoint || this.center;
      const p = [c[0], c[1], c[2]], nx = R[6], ny = R[7], nz = R[8];
      const den = nx * d[0] + ny * d[1] + nz * d[2];
      if (Math.abs(den) < 1e-12) return p;
      const t = (nx * (p[0] - o[0]) + ny * (p[1] - o[1]) + nz * (p[2] - o[2])) / den;
      return [o[0] + t * d[0], o[1] + t * d[1], o[2] + t * d[2]];
    }
    /** CSS pixel -> world ray */
    ray(px, py) {
      const w = this.canvas.clientWidth, h = this.canvas.clientHeight;
      const m = this._matrices(w, h), R = m.R;
      const nx = (2 * px / w) - 1, ny = 1 - (2 * py / h);
      let o, d;
      if (this.opt.perspective) { const t = Math.tan(m.half); o = [0, 0, 0]; d = [nx * t * m.aspect, ny * t, -1]; }
      else { const hh = this.dist * Math.tan(m.half); o = [nx * hh * m.aspect, ny * hh, 0]; d = [0, 0, -1]; }
      const toWorld = v => [R[0] * v[0] + R[3] * v[1] + R[6] * v[2], R[1] * v[0] + R[4] * v[1] + R[7] * v[2], R[2] * v[0] + R[5] * v[1] + R[8] * v[2]];
      const ow = toWorld([o[0], o[1], o[2] + this.dist]);
      const dw = toWorld(d), dl = Math.hypot(...dw);
      return { o: [ow[0] + this.center[0], ow[1] + this.center[1], ow[2] + this.center[2]], d: [dw[0] / dl, dw[1] / dl, dw[2] / dl] };
    }
    pick(px, py) { const { o, d } = this.ray(px, py); return this._pickAtomRay(o, d).index; }
    /** Bond pair index under the pixel (ray within the drawn radius, at least 0.15 Å, of the bond segment), or -1. */
    pickBond(px, py) { const { o, d } = this.ray(px, py); return this._pickBondRay(o, d).index; }
    /** Nearest atom or bond along the ray: {kind: 'atom'|'bond'|'none', index}; atoms win ties. */
    pickAny(px, py) { const { o, d } = this.ray(px, py); return nearerHit(this._pickAtomRay(o, d), this._pickBondRay(o, d)); }
    _pickAtomRay(o, d) {
      const xyz = this.xyz, persp = this.opt.perspective, pr = this._pickR && this._pickR.length === this.Z.length ? this._pickR : null;
      let best = -1, bestT = Infinity;
      for (let i = 0; i < this.Z.length; i++) {
        const r = pr ? pr[i] : this.visible(i) ? Math.max(this.atomRadius(i), 0.28) : 0;
        if (!(r > 0)) continue;
        const lx = xyz[3 * i] - o[0], ly = xyz[3 * i + 1] - o[1], lz = xyz[3 * i + 2] - o[2];
        const tca = lx * d[0] + ly * d[1] + lz * d[2];
        const d2 = lx * lx + ly * ly + lz * lz - tca * tca;
        if (d2 > r * r) continue;
        const s = Math.sqrt(r * r - d2), t = tca - s;
        if (persp && tca + s < 0) continue; // entirely behind the camera
        if (t < bestT) { bestT = t; best = i; }
      }
      return { index: best, t: bestT };
    }
    _pickBondRay(o, d) {
      const xyz = this.xyz, b = this.bonds, n = this.Z.length, persp = this.opt.perspective, g = this._multiGeom();
      let best = -1, bestT = Infinity;
      for (let k = 0, np = b.length >> 1; k < np; k++) {
        const i = b[2 * k], j = b[2 * k + 1];
        if (i >= n || j >= n || !this.visible(i) || !this.visible(j)) continue;
        const r = Math.max(this._bondExtent(k, g), 0.15);
        const t = rayCapsule(o, d, xyz[3 * i], xyz[3 * i + 1], xyz[3 * i + 2], xyz[3 * j], xyz[3 * j + 1], xyz[3 * j + 2], r);
        if (t === Infinity || (persp && t < 0)) continue;
        if (t < bestT) { bestT = t; best = k; }
      }
      return { index: best, t: bestT };
    }
    /** Visible atoms whose projected centres lie inside the CSS-pixel rectangle (corners in any order). */
    atomsInRect(x0, y0, x1, y1) {
      const w = this.canvas.clientWidth, h = this.canvas.clientHeight, m = this._matrices(w, h), xyz = this.xyz, out = [];
      const xa = Math.min(x0, x1), xb = Math.max(x0, x1), ya = Math.min(y0, y1), yb = Math.max(y0, y1), p = [0, 0, 0];
      for (let i = 0; i < this.Z.length; i++) {
        if (!this.visible(i)) continue;
        p[0] = xyz[3 * i]; p[1] = xyz[3 * i + 1]; p[2] = xyz[3 * i + 2];
        const s = this.project(p, m, w, h);
        if (this.opt.perspective && s.z > -1e-6) continue; // behind the camera
        if (s.x >= xa && s.x <= xb && s.y >= ya && s.y <= yb) out.push(i);
      }
      return out;
    }
    /** {kind, index, px, py, world}: world lies on the screen-parallel plane through the hit atom, else through the view centre. */
    _hitAt(px, py) {
      const { o, d } = this.ray(px, py);
      return this._makeHit(px, py, nearerHit(this._pickAtomRay(o, d), this._pickBondRay(o, d)));
    }
    _makeHit(px, py, r) {
      const xyz = this.xyz, through = r.kind === 'atom' ? [xyz[3 * r.index], xyz[3 * r.index + 1], xyz[3 * r.index + 2]] : null;
      return { kind: r.kind, index: r.index, px, py, world: this.screenToWorld(px, py, through) };
    }

    // ---------- input ----------
    /** Active tool: null or {onDown(ev, hit) → boolean, onMove(ev, hit), onUp(ev, hit), onHover?(ev, hit), cursor?}. */
    get tool() { return this._tool; }
    set tool(t) { this._tool = t || null; this._updateCursor(null); }
    _updateCursor(hit) {
      const t = this._tool;
      let c = '';
      if (t && t.cursor) c = typeof t.cursor === 'function' ? (hit ? String(t.cursor(hit) || '') : '') : String(t.cursor);
      if (c !== this._cursor) { this._cursor = c; this.canvas.style.cursor = c; }
    }
    _initInput() {
      const cv = this.canvas;
      let down = null, moved = false, lastClick = 0;
      const local = e => { const rect = cv.getBoundingClientRect(); return [e.clientX - rect.left, e.clientY - rect.top]; };
      // a drag whose onDown returned true sends the tool onMove/onUp until release (dropped if the tool is replaced meanwhile)
      const toolEvent = (name, e) => {
        const t = down && down.tool;
        if (!t || t !== this._tool || typeof t[name] !== 'function') return;
        const [px, py] = local(e), hit = this._hitAt(px, py);
        t[name](e, hit);
        this._updateCursor(hit);
        this.requestRender();
      };
      cv.addEventListener('contextmenu', e => e.preventDefault());
      cv.addEventListener('pointerdown', e => {
        try { cv.setPointerCapture(e.pointerId); } catch (err) { /* synthetic event without an active pointer */ }
        down = { x: e.clientX, y: e.clientY, button: e.button, shift: e.shiftKey, alt: e.altKey, meta: e.metaKey || e.ctrlKey, tool: null };
        moved = false;
        const tool = this._tool;
        if (tool && e.button === 0 && typeof tool.onDown === 'function') { // right/middle drag always pans
          const [px, py] = local(e), hit = this._hitAt(px, py);
          if (tool.onDown(e, hit)) { down.tool = tool; this.requestRender(); }
          this._updateCursor(hit);
        }
      });
      cv.addEventListener('pointermove', e => {
        if (down && down.tool) { toolEvent('onMove', e); return; }
        if (!down) {
          const [px, py] = local(e), tool = this._tool;
          const { o, d } = this.ray(px, py), a = this._pickAtomRay(o, d), i = a.index;
          if (i !== this.hoverAtom) { this.hoverAtom = i; if (this.onHover) this.onHover(i); this.requestRender(); }
          if (tool && (typeof tool.onHover === 'function' || typeof tool.cursor === 'function')) {
            const hit = this._makeHit(px, py, nearerHit(a, this._pickBondRay(o, d)));
            if (typeof tool.onHover === 'function') tool.onHover(e, hit);
            this._updateCursor(hit);
          }
          return;
        }
        const dx = e.clientX - down.x, dy = e.clientY - down.y;
        if (!moved && Math.hypot(dx, dy) < 3) return;
        moved = true;
        down.x = e.clientX; down.y = e.clientY;
        const h = cv.clientHeight || 1;
        if (down.button === 2 || down.button === 1 || down.shift) {
          const m = this._matrices(cv.clientWidth, h);
          const wpp = 2 * this.dist * Math.tan(m.half) / h; // world units per pixel at the target depth
          const R = m.R;
          const vx = -dx * wpp, vy = dy * wpp;
          this.center = [this.center[0] + R[0] * vx + R[3] * vy, this.center[1] + R[1] * vx + R[4] * vy, this.center[2] + R[2] * vx + R[5] * vy];
        } else if (down.alt) {
          this.rotateBy([0, 0, 1], -dx * 0.01);
        } else {
          const ang = Math.hypot(dx, dy) * 0.0085;
          this.rotateBy([dy, dx, 0], ang);
        }
        if (this.onViewChange) this.onViewChange();
        this.requestRender();
      });
      const up = e => {
        if (!down) return;
        if (down.tool) { try { toolEvent('onUp', e); } finally { down = null; } return; }
        const wasMoved = moved, btn = down.button;
        down = null;
        if (wasMoved || btn !== 0) return;
        const rect = cv.getBoundingClientRect();
        const i = this.pick(e.clientX - rect.left, e.clientY - rect.top);
        const now = performance.now();
        const dbl = now - lastClick < 320; lastClick = now;
        if (dbl && this.onAtomDblClick) { this.onAtomDblClick(i, e); return; }
        if (i >= 0) { if (this.onAtomClick) this.onAtomClick(i, e); }
        else if (this.onBackgroundClick) this.onBackgroundClick(e);
      };
      cv.addEventListener('pointerup', up);
      cv.addEventListener('pointercancel', e => { if (down && down.tool) { try { toolEvent('onUp', e); } finally { down = null; } return; } down = null; });
      cv.addEventListener('pointerleave', e => {
        if (down) return;
        if (this.hoverAtom !== -1) { this.hoverAtom = -1; if (this.onHover) this.onHover(-1); this.requestRender(); }
        const tool = this._tool;
        if (tool && typeof tool.onHover === 'function') { const [px, py] = local(e); tool.onHover(e, { kind: 'none', index: -1, px, py, world: this.screenToWorld(px, py) }); }
      });
      cv.addEventListener('wheel', e => {
        e.preventDefault();
        const k = e.ctrlKey ? 0.01 : 0.0015;
        this.dist *= Math.exp(e.deltaY * k);
        this.dist = Math.min(Math.max(this.dist, 0.5), 5000);
        if (this.onViewChange) this.onViewChange();
        this.requestRender();
      }, { passive: false });
    }

    _resize() {
      const dpr = Math.min(window.devicePixelRatio || 1, 3);
      const w = Math.max(1, Math.round(this.canvas.clientWidth * dpr)), h = Math.max(1, Math.round(this.canvas.clientHeight * dpr));
      if (this.canvas.width !== w || this.canvas.height !== h) { this.canvas.width = w; this.canvas.height = h; }
      if (this.overlay.width !== w || this.overlay.height !== h) { this.overlay.width = w; this.overlay.height = h; }
      this.dpr = dpr;
    }

    requestRender() {
      if (this._raf) return;
      this._raf = requestAnimationFrame(() => { this._raf = 0; this.render(); });
    }

    // ---------- drawing ----------
    _drawScene(wPx, hPx, transparent) {
      const gl = this.gl, o = this.opt;
      const cssW = this.canvas.clientWidth || wPx, cssH = this.canvas.clientHeight || hPx;
      const m = this._matrices(cssW, cssH);
      gl.viewport(0, 0, wPx, hPx);
      const bg = o.bg;
      gl.clearColor(transparent ? 0 : bg[0], transparent ? 0 : bg[1], transparent ? 0 : bg[2], transparent ? 0 : 1);
      gl.clear(gl.COLOR_BUFFER_BIT | gl.DEPTH_BUFFER_BIT);
      gl.enable(gl.DEPTH_TEST); gl.depthFunc(gl.LEQUAL); gl.depthMask(true);
      gl.disable(gl.BLEND);
      gl.enable(gl.CULL_FACE); gl.cullFace(gl.BACK); gl.frontFace(gl.CCW);
      const light = (() => { const v = [0.35, 0.55, 0.75], l = Math.hypot(...v); return [v[0] / l, v[1] / l, v[2] / l]; })();
      const fogNear = this.dist - this.radius * 0.3, fogFar = this.dist + this.radius * 1.1;
      const common = (prog, alpha, twoSided) => {
        gl.useProgram(prog.p);
        gl.uniformMatrix4fv(prog.u.uView, false, m.view); gl.uniformMatrix4fv(prog.u.uProj, false, m.proj);
        gl.uniform3fv(prog.u.uLight, light); gl.uniform1f(prog.u.uAmbient, 0.28); gl.uniform1f(prog.u.uSpec, 0.38); gl.uniform1f(prog.u.uShine, 48);
        gl.uniform3fv(prog.u.uFogColor, bg); gl.uniform2f(prog.u.uFogRange, fogNear, fogFar); gl.uniform1f(prog.u.uFogOn, o.fog ? 1 : 0);
        gl.uniform1f(prog.u.uAlpha, alpha); gl.uniform1f(prog.u.uOrtho, o.perspective ? 0 : 1); gl.uniform1f(prog.u.uTwoSided, twoSided ? 1 : 0);
      };
      const cylinders = (vao, count) => {
        if (!count) return;
        common(this.progCyl, 1, false);
        gl.bindVertexArray(vao);
        gl.drawElementsInstanced(gl.TRIANGLES, vao.count, gl.UNSIGNED_SHORT, 0, count);
      };
      // atoms
      if (this.nSphere) {
        common(this.progSphere, 1, false);
        gl.uniform1f(this.progSphere.u.uGrow, 1); gl.uniform1f(this.progSphere.u.uGrowAbs, 0);
        gl.bindVertexArray(this.sphereVAO);
        gl.drawElementsInstanced(gl.TRIANGLES, this.sphereVAO.count, gl.UNSIGNED_SHORT, 0, this.nSphere);
      }
      cylinders(this.cylVAO, this.nMainCyl);
      if (this.nArrowCyl) {
        common(this.progCyl, 1, false);
        gl.bindVertexArray(this.arrowCylVAO);
        gl.drawElementsInstanced(gl.TRIANGLES, this.arrowCylVAO.count, gl.UNSIGNED_SHORT, 0, this.nArrowCyl);
        gl.bindVertexArray(this.coneVAO);
        gl.drawElementsInstanced(gl.TRIANGLES, this.coneVAO.count, gl.UNSIGNED_SHORT, 0, this.nCone);
      }
      cylinders(this.multiVAO, this.nMultiCyl);
      cylinders(this.decoVAO, this.nDeco);
      cylinders(this.xConeVAO, this.nXCone);
      // lines (screen-space quads, constant pixel width)
      if (this.nLine) {
        const p = this.progLine;
        gl.useProgram(p.p);
        gl.uniformMatrix4fv(p.u.uView, false, m.view); gl.uniformMatrix4fv(p.u.uProj, false, m.proj);
        gl.uniform2f(p.u.uViewport, wPx / 2, hPx / 2);
        gl.uniform1f(p.u.uWidth, (o.lineWidth || 1.6) * (wPx / Math.max(1, cssW)));
        gl.uniform1f(p.u.uNear, o.perspective ? m.near * 1.001 : -1e9);
        gl.uniform3fv(p.u.uFogColor, bg); gl.uniform2f(p.u.uFogRange, fogNear, fogFar); gl.uniform1f(p.u.uFogOn, o.fog ? 1 : 0);
        gl.disable(gl.CULL_FACE);
        gl.bindVertexArray(this.lineVAO);
        gl.drawElementsInstanced(gl.TRIANGLES, 6, gl.UNSIGNED_SHORT, 0, this.nLine);
        gl.enable(gl.CULL_FACE);
      }
      // overlay meshes (cartoons, rings, symbols)
      const meshDraw = (ov, alpha) => {
        common(this.progCMesh, alpha, ov.twoSided);
        if (ov.twoSided) gl.disable(gl.CULL_FACE); else gl.enable(gl.CULL_FACE);
        gl.bindVertexArray(ov._vao); gl.drawElements(gl.TRIANGLES, ov._count, gl.UNSIGNED_INT, 0);
      };
      for (const ov of this.overlays) if (ov.opacity >= 0.999) meshDraw(ov, 1);
      gl.enable(gl.CULL_FACE);
      // surfaces
      const surfs = this.surfaces.filter(s => s._vao && s.visible !== false);
      gl.disable(gl.CULL_FACE);
      for (const s of surfs.filter(s => (s.opacity === undefined ? 1 : s.opacity) >= 0.999)) {
        common(this.progMesh, 1, true);
        gl.uniform3fv(this.progMesh.u.uColor, s.color);
        gl.bindVertexArray(s._vao); gl.drawElements(gl.TRIANGLES, s._count, gl.UNSIGNED_INT, 0);
      }
      gl.enable(gl.CULL_FACE);
      // translucent atoms and bonds (ghosts): depth pre-pass, then only the front-most layer is blended
      if (this.nGhostSph || this.nGhostCyl) {
        const ghosts = () => {
          if (this.nGhostSph) {
            common(this.progSphereA, 1, false);
            gl.uniform1f(this.progSphereA.u.uGrow, 1); gl.uniform1f(this.progSphereA.u.uGrowAbs, 0);
            gl.bindVertexArray(this.ghostSphVAO); gl.drawElementsInstanced(gl.TRIANGLES, this.ghostSphVAO.count, gl.UNSIGNED_SHORT, 0, this.nGhostSph);
          }
          if (this.nGhostCyl) { common(this.progCylA, 1, false); gl.bindVertexArray(this.ghostCylVAO); gl.drawElementsInstanced(gl.TRIANGLES, this.ghostCylVAO.count, gl.UNSIGNED_SHORT, 0, this.nGhostCyl); }
        };
        gl.colorMask(false, false, false, false); ghosts(); gl.colorMask(true, true, true, true);
        gl.enable(gl.BLEND); gl.blendFunc(gl.ONE, gl.ONE_MINUS_SRC_ALPHA); gl.depthMask(false);
        ghosts();
        gl.depthMask(true); gl.disable(gl.BLEND);
      }
      // translucent overlay meshes, same scheme
      const tov = this.overlays.filter(ov => ov.opacity < 0.999);
      if (tov.length) {
        gl.colorMask(false, false, false, false); for (const ov of tov) meshDraw(ov, 1); gl.colorMask(true, true, true, true);
        gl.enable(gl.BLEND); gl.blendFunc(gl.ONE, gl.ONE_MINUS_SRC_ALPHA); gl.depthMask(false);
        for (const ov of tov) meshDraw(ov, ov.opacity);
        gl.depthMask(true); gl.disable(gl.BLEND);
      }
      gl.disable(gl.CULL_FACE);
      const trans = surfs.filter(s => (s.opacity === undefined ? 1 : s.opacity) < 0.999);
      if (trans.length) {
        // depth pre-pass so only the front-most layer is blended
        gl.colorMask(false, false, false, false);
        for (const s of trans) { common(this.progMesh, 1, true); gl.uniform3fv(this.progMesh.u.uColor, s.color); gl.bindVertexArray(s._vao); gl.drawElements(gl.TRIANGLES, s._count, gl.UNSIGNED_INT, 0); }
        gl.colorMask(true, true, true, true);
        gl.enable(gl.BLEND); gl.blendFunc(gl.ONE, gl.ONE_MINUS_SRC_ALPHA); gl.depthMask(false);
        for (const s of trans) { common(this.progMesh, s.opacity, true); gl.uniform3fv(this.progMesh.u.uColor, s.color); gl.bindVertexArray(s._vao); gl.drawElements(gl.TRIANGLES, s._count, gl.UNSIGNED_INT, 0); }
        gl.depthMask(true); gl.disable(gl.BLEND);
      }
      gl.enable(gl.CULL_FACE);
      // selection halos
      if (this.nHalo) {
        gl.enable(gl.BLEND); gl.blendFunc(gl.ONE, gl.ONE_MINUS_SRC_ALPHA); gl.depthMask(false);
        common(this.progSphere, 0.42, false);
        gl.uniform1f(this.progSphere.u.uGrow, 1.18); gl.uniform1f(this.progSphere.u.uGrowAbs, 0.07);
        gl.bindVertexArray(this.haloVAO);
        gl.drawElementsInstanced(gl.TRIANGLES, this.haloVAO.count, gl.UNSIGNED_SHORT, 0, this.nHalo);
        gl.depthMask(true); gl.disable(gl.BLEND);
      }
      // highlighted bond
      if (this.nBondHalo) {
        gl.enable(gl.BLEND); gl.blendFunc(gl.ONE, gl.ONE_MINUS_SRC_ALPHA); gl.depthMask(false);
        common(this.progCyl, 0.42, false);
        gl.bindVertexArray(this.bondHaloVAO);
        gl.drawElementsInstanced(gl.TRIANGLES, this.bondHaloVAO.count, gl.UNSIGNED_SHORT, 0, this.nBondHalo);
        gl.depthMask(true); gl.disable(gl.BLEND);
      }
      // atom markers (camera-facing dashed rings)
      if (this.nMarker) {
        const p = this.progRing;
        gl.enable(gl.BLEND); gl.blendFunc(gl.ONE, gl.ONE_MINUS_SRC_ALPHA); gl.depthMask(false); gl.disable(gl.CULL_FACE);
        gl.useProgram(p.p);
        gl.uniformMatrix4fv(p.u.uView, false, m.view); gl.uniformMatrix4fv(p.u.uProj, false, m.proj);
        gl.uniform1f(p.u.uInner, 0.8); gl.uniform1f(p.u.uDashes, 8); gl.uniform1f(p.u.uAlpha, 0.95);
        gl.uniform3fv(p.u.uFogColor, bg); gl.uniform2f(p.u.uFogRange, fogNear, fogFar); gl.uniform1f(p.u.uFogOn, o.fog ? 1 : 0);
        gl.bindVertexArray(this.markerVAO);
        gl.drawElementsInstanced(gl.TRIANGLES, this.markerVAO.count, gl.UNSIGNED_SHORT, 0, this.nMarker);
        gl.enable(gl.CULL_FACE); gl.depthMask(true); gl.disable(gl.BLEND);
      }
      gl.bindVertexArray(null);
      return m;
    }

    _drawOverlay(ctx, m, w, h, scale, noHook) {
      ctx.save();
      ctx.setTransform(scale, 0, 0, scale, 0, 0);
      ctx.clearRect(0, 0, w, h);
      const light = this._isLightBg();
      const txt = light ? '#111827' : '#f3f6fb', halo = light ? 'rgba(255,255,255,0.85)' : 'rgba(5,10,20,0.75)';
      // measurement lines
      if (this.measureLines.length) {
        ctx.setLineDash([5, 4]); ctx.lineWidth = 1.5; ctx.strokeStyle = this.opt.measureColor;
        for (const [i, j] of this.measureLines) {
          if (i >= this.Z.length || j >= this.Z.length) continue;
          const a = this.project([this.xyz[3 * i], this.xyz[3 * i + 1], this.xyz[3 * i + 2]], m, w, h);
          const b = this.project([this.xyz[3 * j], this.xyz[3 * j + 1], this.xyz[3 * j + 2]], m, w, h);
          ctx.beginPath(); ctx.moveTo(a.x, a.y); ctx.lineTo(b.x, b.y); ctx.stroke();
        }
        ctx.setLineDash([]);
      }
      // labels
      const labels = this.labelProvider ? this.labelProvider() : [];
      if (labels.length) {
        ctx.textAlign = 'center'; ctx.textBaseline = 'middle';
        ctx.lineJoin = 'round';
        for (const L of labels) {
          let p = L.pos;
          if (L.atom !== undefined) { const i = L.atom; if (!this.visible(i)) continue; p = [this.xyz[3 * i], this.xyz[3 * i + 1], this.xyz[3 * i + 2]]; }
          const s = this.project(p, m, w, h);
          if (s.z > -0.05 && this.opt.perspective) continue;
          ctx.font = (L.bold ? '600 ' : '500 ') + (L.size || 12) + 'px ui-sans-serif, -apple-system, "Segoe UI", Roboto, sans-serif';
          ctx.lineWidth = 3.2; ctx.strokeStyle = halo; ctx.strokeText(L.text, s.x + (L.dx || 0), s.y + (L.dy || 0));
          ctx.fillStyle = L.color || txt; ctx.fillText(L.text, s.x + (L.dx || 0), s.y + (L.dy || 0));
        }
      }
      if (this.legend) this._drawLegend(ctx, w, light);
      // axes gizmo
      if (this.opt.axes) {
        const R = m.R, ox = 44, oy = h - 44, len = 26;
        const axes = [[R[0], R[3], R[6], '#ef4444', 'x'], [R[1], R[4], R[7], '#22c55e', 'y'], [R[2], R[5], R[8], '#3b82f6', 'z']];
        axes.sort((a, b) => a[2] - b[2]);
        ctx.lineWidth = 2.2; ctx.font = '600 11px ui-sans-serif, -apple-system, sans-serif'; ctx.textAlign = 'center'; ctx.textBaseline = 'middle';
        for (const [x, y, , col, name] of axes) {
          ctx.strokeStyle = col; ctx.beginPath(); ctx.moveTo(ox, oy); ctx.lineTo(ox + x * len, oy - y * len); ctx.stroke();
          ctx.fillStyle = col; ctx.fillText(name, ox + x * (len + 9), oy - y * (len + 9));
        }
      }
      ctx.restore();
      // tool overlay (rubber bands, drag previews): interactive renders only, never in snapshots
      if (this.overlayHook && !noHook) {
        ctx.save();
        ctx.setTransform(scale, 0, 0, scale, 0, 0);
        try { this.overlayHook(ctx, p => this.project(p, m, w, h)); } finally { ctx.restore(); }
      }
    }
    _isLightBg() { const b = this.opt.bg; return 0.299 * b[0] + 0.587 * b[1] + 0.114 * b[2] > 0.55; }
    /** Legend box (top right): a style glyph, a bold label and a muted description per row. */
    _drawLegend(ctx, w, light) {
      const L = this.legend, pad = 8, rowH = 19, gw = 26, font = 'ui-sans-serif, -apple-system, "Segoe UI", Roboto, sans-serif';
      ctx.save();
      ctx.font = '600 11.5px ' + font;
      const labW = Math.max(...L.items.map(it => ctx.measureText(it.label).width));
      ctx.font = '11.5px ' + font;
      const subW = Math.max(0, ...L.items.map(it => (it.sub ? ctx.measureText(it.sub).width : 0)));
      ctx.font = '650 10.5px ' + font;
      const titleW = L.title ? ctx.measureText(L.title.toUpperCase()).width : 0;
      const bw = Math.min(w - 24, Math.max(titleW, gw + labW + (subW ? 8 + subW : 0)) + 2 * pad), bh = pad * 2 + L.items.length * rowH + (L.title ? 16 : 0);
      const x0 = w - bw - 12, y0 = 12;
      ctx.beginPath();
      if (ctx.roundRect) ctx.roundRect(x0, y0, bw, bh, 9); else ctx.rect(x0, y0, bw, bh);
      ctx.fillStyle = light ? 'rgba(255,255,255,0.9)' : 'rgba(10,16,26,0.82)'; ctx.fill();
      ctx.strokeStyle = light ? 'rgba(60,80,110,0.25)' : 'rgba(120,150,200,0.25)'; ctx.lineWidth = 1; ctx.stroke();
      const txt = light ? '#182233' : '#e6edf7', mut = light ? '#5d6b80' : '#9fb0c6';
      let y = y0 + pad;
      ctx.textBaseline = 'middle'; ctx.textAlign = 'left';
      if (L.title) { ctx.font = '650 10.5px ' + font; ctx.fillStyle = mut; ctx.fillText(L.title.toUpperCase(), x0 + pad, y + 6); y += 16; }
      for (const it of L.items) {
        const cy = y + rowH / 2, gx = x0 + pad, col = it.color || txt;
        ctx.strokeStyle = col; ctx.fillStyle = col; ctx.setLineDash([]);
        switch (it.glyph) {
          case 'ballstick': ctx.lineWidth = 2.2; ctx.beginPath(); ctx.moveTo(gx + 4, cy); ctx.lineTo(gx + 17, cy); ctx.stroke(); ctx.beginPath(); ctx.arc(gx + 4, cy, 3.6, 0, 7); ctx.arc(gx + 18, cy, 3.6, 0, 7); ctx.fill(); break;
          case 'sticks': ctx.lineWidth = 4.5; ctx.lineCap = 'round'; ctx.beginPath(); ctx.moveTo(gx + 3, cy); ctx.lineTo(gx + 19, cy); ctx.stroke(); ctx.lineCap = 'butt'; break;
          case 'spacefill': ctx.beginPath(); ctx.arc(gx + 11, cy, 6, 0, 7); ctx.fill(); break;
          case 'ghost': ctx.globalAlpha = 0.45; ctx.beginPath(); ctx.arc(gx + 6, cy, 4, 0, 7); ctx.arc(gx + 16, cy, 4, 0, 7); ctx.fill(); ctx.globalAlpha = 1; break;
          case 'markers': ctx.lineWidth = 1.6; ctx.setLineDash([2.5, 2]); ctx.beginPath(); ctx.arc(gx + 11, cy, 5.5, 0, 7); ctx.stroke(); ctx.setLineDash([]); break;
          case 'dim': ctx.globalAlpha = 0.45; ctx.lineWidth = 1.4; ctx.beginPath(); ctx.moveTo(gx + 2, cy + 3); ctx.lineTo(gx + 9, cy - 3); ctx.lineTo(gx + 15, cy + 3); ctx.lineTo(gx + 21, cy - 3); ctx.stroke(); ctx.globalAlpha = 1; break;
          case 'cartoon': ctx.lineWidth = 4; ctx.beginPath(); ctx.moveTo(gx + 2, cy + 2); ctx.bezierCurveTo(gx + 8, cy - 7, gx + 14, cy + 7, gx + 20, cy - 2); ctx.stroke(); break;
          default: ctx.lineWidth = 1.4; ctx.beginPath(); ctx.moveTo(gx + 2, cy + 3); ctx.lineTo(gx + 9, cy - 3); ctx.lineTo(gx + 15, cy + 3); ctx.lineTo(gx + 21, cy - 3); ctx.stroke();
        }
        ctx.font = '600 11.5px ' + font; ctx.fillStyle = txt; ctx.fillText(it.label, gx + gw, cy);
        if (it.sub) { ctx.font = '11.5px ' + font; ctx.fillStyle = mut; ctx.fillText(it.sub, gx + gw + labW + 8, cy, Math.max(20, x0 + bw - pad - (gx + gw + labW + 8))); }
        y += rowH;
      }
      ctx.restore();
    }

    render() {
      this._resize();
      const m = this._drawScene(this.canvas.width, this.canvas.height, false);
      this._drawOverlay(this.ctx2d, m, this.canvas.clientWidth, this.canvas.clientHeight, this.dpr);
    }

    /** Render at higher resolution and return a PNG Blob. */
    snapshot(scale, transparent) {
      scale = scale || 2;
      const cw = this.canvas.clientWidth, ch = this.canvas.clientHeight;
      const maxDim = this.gl.getParameter(this.gl.MAX_RENDERBUFFER_SIZE) || 8192;
      const s = Math.min(scale, maxDim / Math.max(cw, ch), 8);
      const W = Math.round(cw * s), H = Math.round(ch * s);
      const oldW = this.canvas.width, oldH = this.canvas.height;
      this.canvas.width = W; this.canvas.height = H;
      const m = this._drawScene(W, H, !!transparent);
      const out = document.createElement('canvas'); out.width = W; out.height = H;
      const ctx = out.getContext('2d');
      ctx.drawImage(this.canvas, 0, 0);
      const ov = document.createElement('canvas'); ov.width = W; ov.height = H;
      this._drawOverlay(ov.getContext('2d'), m, cw, ch, s, true);
      ctx.drawImage(ov, 0, 0);
      this.canvas.width = oldW; this.canvas.height = oldH;
      this.render();
      return new Promise(res => out.toBlob(b => res(b), 'image/png'));
    }
  }

  Renderer.STYLE = STYLE;
  BM.Renderer = Renderer;
  BM.quat = { mul: quatMul, axis: quatAxis, mat: quatMat, norm: quatNorm };
})(typeof window !== 'undefined' ? (window.BM = window.BM || {}) : (globalThis.BM = globalThis.BM || {}));
