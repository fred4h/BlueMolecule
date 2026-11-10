// Render tools/icon.svg to PNGs at the sizes needed for macOS (.icns) and Linux.
// Usage: node tools/render_icons.js OUTDIR
const { chromium } = require('playwright');
const fs = require('fs');
const path = require('path');

(async () => {
  const out = process.argv[2];
  fs.mkdirSync(out, { recursive: true });
  const svg = fs.readFileSync(path.join(__dirname, 'icon.svg'), 'utf8');
  const browser = await chromium.launch();
  const page = await browser.newPage();
  await page.setContent('<html><body style="margin:0;background:transparent"></body></html>');
  for (const size of [16, 32, 48, 64, 128, 256, 512, 1024]) {
    const b64 = await page.evaluate(async ({ svg, size }) => {
      const img = new Image();
      img.src = 'data:image/svg+xml;base64,' + btoa(unescape(encodeURIComponent(svg)));
      await img.decode();
      const c = document.createElement('canvas'); c.width = size; c.height = size;
      const ctx = c.getContext('2d'); ctx.imageSmoothingQuality = 'high';
      ctx.drawImage(img, 0, 0, size, size);
      return c.toDataURL('image/png').split(',')[1];
    }, { svg, size });
    fs.writeFileSync(path.join(out, `icon_${size}.png`), Buffer.from(b64, 'base64'));
  }
  await browser.close();
  console.log('icons written to', out);
})();
