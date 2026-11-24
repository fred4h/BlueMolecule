/* BlueMolecule desktop self-test (bluemolecule --self-test=report.json [files…]), run by the cloud builds on each
 * system. It checks what the desktop engine must do: the page and its WebGL 2 view start, the files named on the
 * command line open through the native bridge (one of them big), a save goes through the host, and a picture of the
 * view is taken. The host writes the report (and report.png) and exits: 0 when every check passed. */
(function () {
  'use strict';
  if (!window.BM_TAURI || !window.BM_TAURI.selfTest || window.top !== window) return;
  const rep = { started: new Date().toISOString(), checks: {}, errors: [], files: [] };
  const invoke = (cmd, args) => window.__TAURI__.core.invoke(cmd, args);
  window.addEventListener('error', e => rep.errors.push(String((e && e.message) || e)));
  window.addEventListener('unhandledrejection', e => rep.errors.push('unhandled rejection: ' + String((e && e.reason && e.reason.message) || (e && e.reason))));
  const wait = (test, ms, what) => new Promise((resolve, reject) => {
    const t0 = Date.now();
    (function poll() {
      let v = null;
      try { v = test(); } catch (e) { v = null; }
      if (v) return resolve(v);
      if (Date.now() - t0 > ms) return reject(new Error('timed out waiting for ' + what));
      setTimeout(poll, 100);
    })();
  });
  const toB64 = async blob => {
    const u8 = new Uint8Array(await blob.arrayBuffer());
    let s = '';
    for (let i = 0; i < u8.length; i += 0x8000) s += String.fromCharCode.apply(null, u8.subarray(i, i + 0x8000));
    return btoa(s);
  };
  async function run() {
    const t0 = performance.now();
    await wait(() => window.BM && window.BM.app && window.BM.native, 60000, 'the app');
    const BM = window.BM, app = BM.app;
    rep.version = app.version;
    rep.userAgent = navigator.userAgent;
    rep.nativeKind = BM.native.kind;
    rep.startSeconds = +((performance.now() - t0) / 1000).toFixed(2);
    rep.checks.bridge = BM.native.kind === 'tauri-app';
    const R = app.R();
    const gl = R && R.gl;
    rep.checks.webgl2 = !!gl;
    if (gl) {
      const dbg = gl.getExtension('WEBGL_debug_renderer_info');
      rep.gpu = String(gl.getParameter(dbg ? dbg.UNMASKED_RENDERER_WEBGL : gl.RENDERER));
      rep.maxTexture = gl.getParameter(gl.MAX_TEXTURE_SIZE);
    }
    rep.memory = { total: window.BM_TAURI.memTotal, available: window.BM_TAURI.memAvailable,
      jsHeapLimit: performance.memory ? performance.memory.jsHeapSizeLimit : null };
    // the files of the command line, opened by the host through bm_ready / bm_read
    const want = (window.BM_TAURI.selfTestFiles || []).slice();
    const t1 = performance.now();
    await wait(() => app.state.docs.length >= Math.max(1, want.length) && app.state.docs.every(d => !d.pending), 600000, 'the files');
    rep.openSeconds = +((performance.now() - t1) / 1000).toFixed(2);
    for (const d of app.state.docs) {
      const f = d.data.frames[d.data.frames.length - 1];
      rep.files.push({ name: d.name, frames: d.data.frames.length, atoms: f.Z ? f.Z.length : 0, vibrations: d.data.vibrations ? d.data.vibrations.freqs.length : 0,
        energy: f.energy, warnings: (d.data.warnings || []).slice(0, 5) });
    }
    rep.checks.files = rep.files.length > 0 && rep.files.every(f => f.atoms > 0);
    // a save through the host (written next to the report in self-test mode)
    rep.checks.save = await BM.native.save('self-test.xyz', new Blob(['3\nwater\nO 0 0 0.117\nH 0 0.757 -0.469\nH 0 -0.757 -0.469\n'], { type: 'text/plain' })) === true;
    // downloads from PubChem and the PDB in this system's web engine (issue #22). Reported, not required to pass:
    // the cloud runner's network and the two services are outside BlueMolecule.
    if (BM.fetchdb) {
      rep.downloads = {};
      const attempt = async (key, fn, sum) => {
        const s = performance.now();
        try { const r = await fn(); rep.downloads[key] = Object.assign({ ok: true, ms: Math.round(performance.now() - s) }, sum(r)); }
        catch (e) { rep.downloads[key] = { ok: false, ms: Math.round(performance.now() - s), error: String((e && e.message) || e) }; }
      };
      const read = r => BM.parseFile(r.fileName, r.text);
      await attempt('pubchem', () => BM.fetchdb.pubchem('aspirin', { timeout: 30000 }), r => { const ds = read(r); return { cid: r.cid, atoms: ds.frames[0].Z.length, formula: BM.geom.formula(ds.frames[0].Z) }; });
      await attempt('pubchemSmiles', () => BM.fetchdb.pubchem('sodium chloride', { timeout: 30000 }), r => ({ cid: r.cid, conformer: !!r.text, smiles: r.smiles || null }));
      await attempt('pdb', () => BM.fetchdb.pdb('1CRN', { timeout: 30000 }), r => { const ds = read(r); return { atoms: ds.frames[0].Z.length, title: ds.title, bytes: r.text.length }; });
      await attempt('pdbAssembly', () => BM.fetchdb.pdb('1CRN', { timeout: 30000, assembly: true }), r => ({ atoms: read(r).frames[0].Z.length }));
      await attempt('missing', () => BM.fetchdb.pdb('0XXX', { timeout: 30000 }), () => ({ note: 'found?' }));
    }
    // a picture of the smallest molecule
    let image = null;
    const small = app.state.docs.slice().sort((a, b) => a.data.frames[0].Z.length - b.data.frames[0].Z.length)[0];
    if (small && R && R.snapshot) {
      app.activate(small);
      await new Promise(r => setTimeout(r, 1500));
      const blob = await R.snapshot(1, false);
      if (blob && blob.size > 2000) { image = await toB64(blob); rep.imageBytes = blob.size; }
    }
    rep.checks.picture = !!image;
    rep.seconds = +((performance.now() - t0) / 1000).toFixed(2);
    rep.ok = Object.values(rep.checks).every(Boolean) && rep.errors.length === 0;
    await invoke('bm_self_test_done', { report: rep, image });
  }
  const start = () => run().catch(e => {
    rep.error = String((e && e.stack) || e);
    rep.ok = false;
    invoke('bm_self_test_done', { report: rep, image: null });
  });
  if (document.readyState === 'loading') document.addEventListener('DOMContentLoaded', start); else start();
})();
