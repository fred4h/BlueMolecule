// Render tools/icon.svg into the PNG icons embedded in the <head> of src/index.html (favicons at
// 32, 192 and 512 px plus a 180 px apple-touch-icon) and replace the page's icon links with them.
// Usage: node native/make_favicons.js            (re-run whenever tools/icon.svg changes)
// Only the block between <!--BM:ICONS--> and <!--/BM:ICONS--> (or, the first time, the existing
// <link rel="icon"> line) is rewritten.
const { chromium } = require('playwright');
const fs = require('fs');
const path = require('path');
const zlib = require('zlib');

const ROOT = path.resolve(__dirname, '..');
const INDEX = path.join(ROOT, 'src', 'index.html');

/** Edge-preserving smoothing: where a 3x3 neighbourhood varies by at most `t` levels in every channel
 *  (smooth gradient plus the renderer's dither noise), use its mean. Edges are left untouched. */
function smooth(w, h, src, t) {
  const out = Buffer.from(src);
  for (let y = 1; y < h - 1; y++) {
    for (let x = 1; x < w - 1; x++) {
      let flat = true;
      const sum = [0, 0, 0, 0];
      for (let c = 0; c < 4 && flat; c++) {
        let lo = 255, hi = 0;
        for (let dy = -1; dy <= 1; dy++) for (let dx = -1; dx <= 1; dx++) {
          const v = src[((y + dy) * w + x + dx) * 4 + c];
          if (v < lo) lo = v; if (v > hi) hi = v; sum[c] += v;
        }
        if (hi - lo > t) flat = false;
      }
      if (flat) for (let c = 0; c < 4; c++) out[(y * w + x) * 4 + c] = Math.round(sum[c] / 9);
    }
  }
  return out;
}

/** Compact PNG (RGBA, per-row adaptive filter, deflate level 9). */
function encodePNG(w, h, rgba) {
  const stride = w * 4, raw = Buffer.alloc((stride + 1) * h), zero = Buffer.alloc(stride);
  for (let i = 0; i < rgba.length; i += 4) if (rgba[i + 3] === 0) rgba[i] = rgba[i + 1] = rgba[i + 2] = 0;
  const cand = [0, 1, 2, 3, 4].map(() => Buffer.alloc(stride));
  for (let y = 0; y < h; y++) {
    const row = rgba.subarray(y * stride, (y + 1) * stride), up = y ? rgba.subarray((y - 1) * stride, y * stride) : zero;
    let best = 0, bestSum = Infinity;
    for (let t = 0; t < 5; t++) {
      const out = cand[t];
      let sum = 0;
      for (let i = 0; i < stride; i++) {
        const a = i >= 4 ? row[i - 4] : 0, b = up[i], c = i >= 4 ? up[i - 4] : 0;
        let p;
        if (t === 0) p = 0; else if (t === 1) p = a; else if (t === 2) p = b; else if (t === 3) p = (a + b) >> 1;
        else { const q = a + b - c, pa = Math.abs(q - a), pb = Math.abs(q - b), pc = Math.abs(q - c); p = pa <= pb && pa <= pc ? a : pb <= pc ? b : c; }
        const v = (row[i] - p) & 255;
        out[i] = v; sum += v < 128 ? v : 256 - v;
      }
      if (sum < bestSum) { bestSum = sum; best = t; }
    }
    raw[y * (stride + 1)] = best;
    cand[best].copy(raw, y * (stride + 1) + 1);
  }
  const chunk = (type, data) => {
    const len = Buffer.alloc(4); len.writeUInt32BE(data.length);
    const td = Buffer.concat([Buffer.from(type, 'ascii'), data]);
    const crc = Buffer.alloc(4); crc.writeUInt32BE(zlib.crc32(td) >>> 0);
    return Buffer.concat([len, td, crc]);
  };
  const ihdr = Buffer.alloc(13);
  ihdr.writeUInt32BE(w, 0); ihdr.writeUInt32BE(h, 4); ihdr[8] = 8; ihdr[9] = 6;
  return Buffer.concat([Buffer.from([137, 80, 78, 71, 13, 10, 26, 10]), chunk('IHDR', ihdr),
    chunk('IDAT', zlib.deflateSync(raw, { level: 9, memLevel: 9 })), chunk('IEND', Buffer.alloc(0))]);
}

(async () => {
  const svg = fs.readFileSync(path.join(ROOT, 'tools', 'icon.svg'), 'utf8');
  const browser = await chromium.launch();
  const page = await browser.newPage();
  await page.setContent('<html><body style="margin:0;background:transparent"></body></html>');
  // favicons: cropped to the rounded square (the .icns keeps the macOS margin); touch icon: full bleed, opaque
  const jobs = [
    { rel: 'icon', size: 32, viewBox: '96 96 832 832' },
    { rel: 'icon', size: 192, viewBox: '96 96 832 832' },
    { rel: 'icon', size: 512, viewBox: '96 96 832 832' },
    { rel: 'apple-touch-icon', size: 180, viewBox: '100 100 824 824', fill: '#0f2547' },
  ];
  const links = [];
  for (const j of jobs) {
    const src = svg.replace(/viewBox="[^"]*"/, `viewBox="${j.viewBox}"`);
    const rgba = await page.evaluate(async ({ src, size, fill }) => {
      const img = new Image();
      img.src = 'data:image/svg+xml;base64,' + btoa(unescape(encodeURIComponent(src)));
      await img.decode();
      const c = document.createElement('canvas'); c.width = c.height = size;
      const ctx = c.getContext('2d'); ctx.imageSmoothingQuality = 'high';
      if (fill) { ctx.fillStyle = fill; ctx.fillRect(0, 0, size, size); }
      ctx.drawImage(img, 0, 0, size, size);
      return Array.from(ctx.getImageData(0, 0, size, size).data);
    }, { src, size: j.size, fill: j.fill || null });
    const png = encodePNG(j.size, j.size, j.size >= 128 ? smooth(j.size, j.size, Buffer.from(rgba), 6) : Buffer.from(rgba));
    console.log(`${j.rel} ${j.size}px: ${(png.length / 1024).toFixed(1)} KB`);
    const type = j.rel === 'icon' ? ' type="image/png"' : '';
    // the URL parser drops newlines, so the base64 is wrapped to keep index.html readable
    const b64 = png.toString('base64').replace(/.{1,160}/g, '$&\n').trimEnd();
    links.push(`<link rel="${j.rel}"${type} sizes="${j.size}x${j.size}" href="data:image/png;base64,\n${b64}">`);
  }
  await browser.close();

  const block = '<!--BM:ICONS-->\n' + links.join('\n') + '\n<!--/BM:ICONS-->';
  let html = fs.readFileSync(INDEX, 'utf8');
  if (/<!--BM:ICONS-->[\s\S]*?<!--\/BM:ICONS-->/.test(html)) html = html.replace(/<!--BM:ICONS-->[\s\S]*?<!--\/BM:ICONS-->/, () => block);
  else if (/<link rel="icon"[^>]*>/.test(html)) html = html.replace(/<link rel="icon"[^>]*>/, () => block);
  else throw new Error('no icon link found in the <head> of src/index.html');
  fs.writeFileSync(INDEX, html);
  console.log('icon links written to', INDEX);
})().catch(e => { console.error(e); process.exit(1); });
