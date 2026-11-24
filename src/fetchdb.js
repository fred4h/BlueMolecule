/* BlueMolecule — structures downloaded from PubChem and the RCSB Protein Data Bank (BM.fetchdb), when the user allows it.
 *
 *   BM.fetchdb.parse(source, text)   -> {kind: 'cid' | 'name' | 'pdb', value} or {error}   (source: 'pubchem' or 'pdb')
 *   BM.fetchdb.pubchem(query, opts)  -> Promise<{source, name, cid, url, fileName, text}>, the first 3D conformer as SDF,
 *                                       or <{source, name, cid, url, fileName, smiles}> when PubChem has no 3D conformer
 *   BM.fetchdb.pdb(id, opts)         -> Promise<{source, id, name, url, fileName, text, assembly}>, the entry as mmCIF;
 *                                       opts.assembly: the first biological assembly instead of the deposited coordinates
 *   opts: fetch (the fetch function; tests), signal (an AbortSignal; rejects with 'cancelled'), timeout (ms, 60 s
 *         without any data), onProgress(bytesReceived)
 *
 * Only what the user typed goes out, as part of the address: to pubchem.ncbi.nlm.nih.gov (PUG REST) or files.rcsb.org.
 * No cookies, referrer or other headers of BlueMolecule's own are sent. Both services answer requests from any page (CORS), so this works the same in
 * the desktop window and in the single HTML file. The app asks the user before the first download (app.js).
 */
(function (BM) {
  'use strict';
  const PUG = 'https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/';
  const RCSB = 'https://files.rcsb.org/download/';

  function parse(source, text) {
    const t = String(text == null ? '' : text).trim();
    if (source === 'pdb') {
      if (!t) return { error: 'Type a PDB ID, e.g. 1CRN.' };
      const id = t.replace(/\.(cif|pdb)$/i, '');
      if (/^[0-9][A-Za-z0-9]{3}$/.test(id)) return { kind: 'pdb', value: id.toUpperCase() };
      if (/^pdb_[0-9a-z]{8}$/i.test(id)) return { kind: 'pdb', value: id.toLowerCase() };
      return { error: 'A PDB ID is four characters, a digit and three letters or digits (e.g. 1CRN), or the extended form (e.g. pdb_00001crn).' };
    }
    if (!t) return { error: 'Type a compound name or a PubChem CID, e.g. caffeine or 2519.' };
    if (/^(cid[:\s]*)?\d{1,10}$/i.test(t)) {
      const v = t.replace(/\D/g, '').replace(/^0+(?=\d)/, '');
      return v === '0' ? { error: 'PubChem CIDs start at 1.' } : { kind: 'cid', value: v };
    }
    if (t.length > 300) return { error: 'That name is too long (300 characters at most).' };
    return { kind: 'name', value: t };
  }

  /** One GET request: {ok, status, text}. Network failures, time-outs and cancelling become plain errors. */
  async function get(url, opts) {
    opts = opts || {};
    const f = opts.fetch || (typeof fetch === 'function' ? fetch.bind(globalThis) : null);
    if (!f) throw new Error('Downloads are not possible here (no fetch).');
    if (opts.signal && opts.signal.aborted) throw new Error('cancelled');
    const ctl = typeof AbortController === 'function' ? new AbortController() : null;
    const ms = opts.timeout || 60000;
    let timedOut = false, timer = 0;
    const arm = () => { if (!ctl) return; clearTimeout(timer); timer = setTimeout(() => { timedOut = true; ctl.abort(); }, ms); };
    const onAbort = () => { if (ctl) ctl.abort(); };
    if (opts.signal) opts.signal.addEventListener('abort', onAbort);
    const fail = e => {
      if (opts.signal && opts.signal.aborted) return new Error('cancelled');
      if (timedOut) return new Error('The server did not answer in time.');
      return new Error('The server could not be reached: no network connection, or blocked by a firewall or proxy' + (e && e.message ? ' (' + e.message + ')' : '') + '.');
    };
    try {
      arm();
      let r;
      try { r = await f(url, { signal: ctl ? ctl.signal : undefined, credentials: 'omit', referrerPolicy: 'no-referrer', cache: 'no-store', redirect: 'follow' }); }
      catch (e) { throw fail(e); }
      arm();
      let text;
      try {
        if (opts.onProgress && r.body && typeof r.body.getReader === 'function' && typeof TextDecoder === 'function') {
          const rd = r.body.getReader(), dec = new TextDecoder(), parts = [];
          let got = 0;
          for (;;) {
            arm();                                     // the time-out counts from the last data
            const c = await rd.read();
            if (c.done) break;
            got += c.value.length;
            parts.push(dec.decode(c.value, { stream: true }));
            opts.onProgress(got);
          }
          parts.push(dec.decode());
          text = parts.join('');
        } else text = await r.text();
      } catch (e) { throw fail(e); }
      return { ok: !!r.ok, status: r.status, text };
    } finally {
      clearTimeout(timer);
      if (opts.signal) opts.signal.removeEventListener('abort', onAbort);
    }
  }
  const busy = s => s === 503 || s === 429;
  const safeName = s => String(s).replace(/[^\p{L}\p{N}.+-]+/gu, '_').replace(/^_+|_+$/g, '').slice(0, 80) || 'compound';

  async function pubchem(query, opts) {
    const q = parse('pubchem', query);
    if (q.error) throw new Error(q.error);
    const ref = q.kind === 'cid' ? 'cid/' + q.value : 'name/' + encodeURIComponent(q.value);
    const label = q.kind === 'cid' ? 'CID ' + q.value : '“' + q.value + '”';
    const named = cid => ({ name: q.kind === 'name' ? q.value : 'CID ' + cid, fileName: (q.kind === 'name' ? safeName(q.value) + '_' : '') + 'CID' + cid });
    const url = PUG + ref + '/SDF?record_type=3d';
    const r = await get(url, opts);
    if (r.ok) {
      const end = r.text.indexOf('$$$$');
      if (end < 0) throw new Error('PubChem sent something that is not an SD file.');
      // the first record only (a name can match several compounds); its first line is the CID
      const rec = r.text.slice(0, end + 4) + '\n';
      const first = (rec.split(/\r?\n/)[0] || '').trim();
      const cid = /^\d+$/.test(first) ? first : q.kind === 'cid' ? q.value : '';
      const n = named(cid || '?');
      return { source: 'pubchem', name: n.name, cid, url, fileName: n.fileName + '.sdf', text: rec };
    }
    if (busy(r.status)) throw new Error('PubChem is busy; try again in a moment.');
    if (r.status !== 404 && r.status !== 400) throw new Error('PubChem answered with an error (HTTP ' + r.status + ').');
    // no 3D conformer (salts and mixtures, large or very flexible molecules, some elements), or no such compound
    for (const prop of ['SMILES', 'IsomericSMILES']) {
      const purl = PUG + ref + '/property/' + prop + '/JSON';
      const s = await get(purl, opts);
      if (s.ok) {
        let p = null;
        try { p = JSON.parse(s.text).PropertyTable.Properties[0]; } catch (e) { p = null; }
        const smiles = p && (p.SMILES || p.IsomericSMILES || p.CanonicalSMILES || p.ConnectivitySMILES);
        if (smiles) {
          const cid = p.CID ? String(p.CID) : q.kind === 'cid' ? q.value : '';
          const n = named(cid || '?');
          return { source: 'pubchem', name: n.name, cid, url: purl, fileName: n.fileName, smiles: String(smiles) };
        }
      }
      if (busy(s.status)) throw new Error('PubChem is busy; try again in a moment.');
      if (s.status === 404) break;
    }
    throw new Error('PubChem has no compound ' + (q.kind === 'cid' ? 'with ' : 'named ') + label + '.');
  }

  async function pdb(id, opts) {
    const q = parse('pdb', id);
    if (q.error) throw new Error(q.error);
    const asm = !!(opts && opts.assembly);
    const file = q.value + (asm ? '-assembly1' : '') + '.cif';
    const url = RCSB + file;
    const r = await get(url, opts);
    if (r.ok) {
      if (!/^\s*data_/m.test(r.text)) throw new Error('The PDB sent something that is not an mmCIF file.');
      return { source: 'pdb', id: q.value, name: q.value + (asm ? ' assembly 1' : ''), url, fileName: file, text: r.text, assembly: asm };
    }
    if (r.status === 404) throw new Error(asm ? 'The PDB has no entry ' + q.value + ', or no biological assembly for it.' : 'The PDB has no entry ' + q.value + '.');
    if (busy(r.status)) throw new Error('The PDB is busy; try again in a moment.');
    throw new Error('The PDB answered with an error (HTTP ' + r.status + ').');
  }

  BM.fetchdb = { parse, pubchem, pdb, get, PUG, RCSB, hosts: ['pubchem.ncbi.nlm.nih.gov', 'files.rcsb.org'] };
})(typeof window !== 'undefined' ? (window.BM = window.BM || {}) : (globalThis.BM = globalThis.BM || {}));
