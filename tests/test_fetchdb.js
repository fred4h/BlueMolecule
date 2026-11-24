#!/usr/bin/env node
/* BlueMolecule tests for the downloads from PubChem and the PDB (src/fetchdb.js), with a stand-in for fetch: the
 * addresses asked for, what is (not) sent, the 3D-conformer and SMILES paths, errors, time-outs, cancelling and
 * progress. Nothing here touches the network.
 * Run:  node tests/test_fetchdb.js */
'use strict';
const fs = require('fs');
const path = require('path');
const ROOT = path.join(__dirname, '..');
for (const m of ['elements', 'bigfile', 'parsers', 'condense', 'readers-struct', 'writers', 'readers-qc', 'geom', 'chem', 'ff', 'fetchdb']) {
  require(path.join(ROOT, 'src', m + '.js'));
}
const BM = globalThis.BM, F = BM.fetchdb;

let pass = 0, fail = 0, current = '';
const failures = [];
function ok(c, msg) { if (c) pass++; else { fail++; failures.push('[' + current + '] ' + msg); } }
function eq(a, b, msg) { ok(a === b, msg + ': got ' + JSON.stringify(a) + ', expected ' + JSON.stringify(b)); }
function test(name, fn) { current = name; try { fn(); } catch (e) { fail++; failures.push('[' + name + '] threw: ' + (e && e.stack || e)); } }
async function testAsync(name, fn) { current = name; try { await fn(); } catch (e) { fail++; failures.push('[' + name + '] threw: ' + (e && e.stack || e)); } }
async function rejects(p, re, msg) {
  let m = null;
  try { await p; } catch (e) { m = e && e.message; }
  ok(m !== null && re.test(m), msg + ': ' + (m === null ? 'did not fail' : 'got “' + m + '”'));
}

const PUG = 'https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/', RCSB = 'https://files.rcsb.org/download/';

/** An SD record the way PubChem writes it: the CID as the title line, then the molfile and data items. */
function sdRecord(cid, smiles) {
  const m = BM.chem.smilesTo3D(smiles, { seed: 3 }).mol;
  const mol = BM.writers.mol({ Z: m.Z, xyz: m.xyz, bonds: m.bonds.map(b => [b[0], b[1], Math.round(b[2]) || 1]), title: String(cid) });
  return mol.replace(/\n?$/, '\n') + '> <PUBCHEM_COMPOUND_CID>\n' + cid + '\n\n> <PUBCHEM_CONFORMER_RMSD>\n0.6\n\n$$$$\n';
}
const ASPIRIN = sdRecord(2244, 'CC(=O)Oc1ccccc1C(=O)O');
const SALICYLIC = sdRecord(338, 'OC(=O)c1ccccc1O');
const MMCIF = fs.readFileSync(path.join(__dirname, 'fixtures', 'struct', 'small.mmcif'), 'utf8');

/** A stand-in for fetch: routes(url) -> {status, body} | Error | 'hang'; every call is recorded with its options. */
function mockFetch(routes) {
  const calls = [];
  const f = (url, init) => {
    calls.push({ url, init: init || {} });
    const r = routes(url);
    if (r instanceof Error) return Promise.reject(r);
    if (r === 'hang') {
      return new Promise((res, rej) => {
        const s = init && init.signal;
        if (s) s.addEventListener('abort', () => { const e = new Error('The operation was aborted.'); e.name = 'AbortError'; rej(e); });
      });
    }
    const status = r ? r.status : 404, body = r ? r.body : '{"Fault":{"Code":"PUGREST.NotFound"}}';
    return Promise.resolve({ ok: status >= 200 && status < 300, status, headers: { get: () => null }, body: r && r.stream ? r.stream : null, text: async () => body });
  };
  f.calls = calls;
  return f;
}
/** A response body read in pieces (bytes), like a real download. */
function streamOf(chunks) {
  let k = 0;
  return { getReader: () => ({ read: async () => (k < chunks.length ? { done: false, value: chunks[k++] } : { done: true, value: undefined }) }) };
}

(async () => {
  test('queries', () => {
    eq(F.parse('pdb', '1crn').value, '1CRN', 'PDB ID in capitals');
    eq(F.parse('pdb', ' 4hhb.cif ').value, '4HHB', 'a file name is accepted');
    eq(F.parse('pdb', 'PDB_00001CRN').value, 'pdb_00001crn', 'extended PDB ID');
    for (const bad of ['', 'crn1', '12345', '1cr', '1c-n', 'pdb_1crn', '1CRN; rm -rf']) ok(!!F.parse('pdb', bad).error, 'rejected PDB ID “' + bad + '”');
    eq(F.parse('pubchem', '2244').kind, 'cid', 'a number is a CID');
    eq(F.parse('pubchem', 'CID 2244').value, '2244', 'CID 2244');
    eq(F.parse('pubchem', 'cid:002244').value, '2244', 'cid:002244');
    ok(!!F.parse('pubchem', '0').error, 'CID 0');
    eq(F.parse('pubchem', '  caffeine ').value, 'caffeine', 'a name, trimmed');
    eq(F.parse('pubchem', '2,4-dinitrophenol').kind, 'name', 'a name with digits');
    ok(!!F.parse('pubchem', '').error && !!F.parse('pubchem', '   ').error, 'empty');
    ok(!!F.parse('pubchem', 'x'.repeat(301)).error, 'too long');
    ok(F.hosts.length === 2 && F.PUG.startsWith('https://' + F.hosts[0] + '/') && F.RCSB.startsWith('https://' + F.hosts[1] + '/'), 'two HTTPS hosts');
  });

  await testAsync('PubChem by name: the first 3D conformer', async () => {
    const f = mockFetch(u => (u === PUG + 'name/aspirin/SDF?record_type=3d' ? { status: 200, body: ASPIRIN + SALICYLIC } : null));
    const r = await F.pubchem('aspirin', { fetch: f });
    eq(f.calls.length, 1, 'one request');
    eq(f.calls[0].init.credentials, 'omit', 'no cookies');
    eq(f.calls[0].init.referrerPolicy, 'no-referrer', 'no referrer');
    ok(!f.calls[0].init.headers && !f.calls[0].init.body && !f.calls[0].init.method, 'a plain GET: no headers, no body');
    eq(r.source, 'pubchem', 'source'); eq(r.cid, '2244', 'CID from the title line'); eq(r.name, 'aspirin', 'name');
    eq(r.fileName, 'aspirin_CID2244.sdf', 'file name'); eq(r.url, PUG + 'name/aspirin/SDF?record_type=3d', 'url');
    ok(r.text.endsWith('$$$$\n') && r.text.indexOf('$$$$') === r.text.lastIndexOf('$$$$'), 'only the first record');
    const ds = BM.parseFile(r.fileName, r.text);
    eq(ds.frames.length, 1, 'one structure'); eq(ds.frames[0].Z.length, 21, 'aspirin has 21 atoms');
    eq(BM.geom.formula(ds.frames[0].Z), 'C9H8O4', 'formula');
  });

  await testAsync('PubChem by CID; names are encoded', async () => {
    const f = mockFetch(u => (/\/SDF\?record_type=3d$/.test(u) ? { status: 200, body: ASPIRIN } : null));
    const r = await F.pubchem('CID 2244', { fetch: f });
    eq(f.calls[0].url, PUG + 'cid/2244/SDF?record_type=3d', 'cid url'); eq(r.name, 'CID 2244', 'name'); eq(r.fileName, 'CID2244.sdf', 'file name');
    const cases = [['acetic acid', 'acetic%20acid', 'acetic_acid_CID2244.sdf'], ['2,4-D', '2%2C4-D', '2_4-D_CID2244.sdf'],
      ['β-carotene', '%CE%B2-carotene', 'β-carotene_CID2244.sdf'], ['a/b?c#d&e', 'a%2Fb%3Fc%23d%26e', 'a_b_c_d_e_CID2244.sdf']];
    for (const [q, enc, file] of cases) {
      const g = mockFetch(u => (/\/SDF\?record_type=3d$/.test(u) ? { status: 200, body: ASPIRIN } : null));
      const s = await F.pubchem(q, { fetch: g });
      eq(g.calls[0].url, PUG + 'name/' + enc + '/SDF?record_type=3d', 'encoded “' + q + '”');
      eq(s.fileName, file, 'file name for “' + q + '”');
    }
  });

  await testAsync('PubChem without a 3D conformer: the SMILES', async () => {
    const f = mockFetch(u => (u === PUG + 'name/sodium%20chloride/property/SMILES/JSON'
      ? { status: 200, body: '{"PropertyTable":{"Properties":[{"CID":5234,"SMILES":"[Na+].[Cl-]"}]}}' } : null));
    const r = await F.pubchem('sodium chloride', { fetch: f });
    eq(f.calls.length, 2, 'the conformer, then the SMILES');
    eq(f.calls[0].url, PUG + 'name/sodium%20chloride/SDF?record_type=3d', 'asked for the conformer first');
    eq(r.smiles, '[Na+].[Cl-]', 'smiles'); eq(r.cid, '5234', 'cid'); ok(!r.text, 'no file'); eq(r.fileName, 'sodium_chloride_CID5234', 'name');
    ok(BM.chem.smilesTo3D(r.smiles).mol.Z.length === 2, 'the SMILES builds');
    // the property's older name, when the newer one is refused
    const g = mockFetch(u => (/\/property\/SMILES\/JSON$/.test(u) ? { status: 400, body: '{"Fault":{"Code":"PUGREST.BadRequest"}}' }
      : /\/property\/IsomericSMILES\/JSON$/.test(u) ? { status: 200, body: '{"PropertyTable":{"Properties":[{"CID":702,"IsomericSMILES":"CCO"}]}}' } : null));
    const s = await F.pubchem('702', { fetch: g });
    eq(g.calls.length, 3, 'three requests'); eq(s.smiles, 'CCO', 'IsomericSMILES'); eq(s.cid, '702', 'cid'); eq(s.fileName, 'CID702', 'name');
  });

  await testAsync('PubChem errors', async () => {
    const none = mockFetch(() => null);
    await rejects(F.pubchem('xyzzy', { fetch: none }), /PubChem has no compound named “xyzzy”/, 'unknown name');
    eq(none.calls.length, 2, 'no third request after a 404');
    await rejects(F.pubchem('99999999', { fetch: mockFetch(() => null) }), /no compound with CID 99999999/, 'unknown CID');
    await rejects(F.pubchem('aspirin', { fetch: mockFetch(() => ({ status: 503, body: 'busy' })) }), /PubChem is busy/, '503');
    await rejects(F.pubchem('aspirin', { fetch: mockFetch(u => (/SDF/.test(u) ? null : { status: 503, body: 'busy' })) }), /PubChem is busy/, '503 on the SMILES');
    await rejects(F.pubchem('aspirin', { fetch: mockFetch(() => ({ status: 500, body: 'oops' })) }), /HTTP 500/, '500');
    await rejects(F.pubchem('aspirin', { fetch: mockFetch(() => ({ status: 200, body: '<html>Sign in to the network</html>' })) }), /not an SD file/, 'a captive portal page');
    await rejects(F.pubchem('aspirin', { fetch: mockFetch(() => new TypeError('Failed to fetch')) }), /could not be reached/, 'no network');
    await rejects(F.pubchem('', { fetch: mockFetch(() => null) }), /Type a compound name/, 'empty query');
    const nothing = mockFetch(() => null);
    await rejects(F.pdb('nope', { fetch: nothing }), /PDB ID is four characters/, 'bad PDB ID');
    eq(nothing.calls.length, 0, 'a bad query sends nothing');
  });

  await testAsync('PDB entries and assemblies', async () => {
    const f = mockFetch(u => (u === RCSB + '1ABC.cif' || u === RCSB + '1ABC-assembly1.cif' || u === RCSB + 'pdb_00001abc.cif' ? { status: 200, body: MMCIF } : null));
    const r = await F.pdb('1abc', { fetch: f });
    eq(f.calls[0].url, RCSB + '1ABC.cif', 'entry url'); eq(r.fileName, '1ABC.cif', 'file'); eq(r.id, '1ABC', 'id'); eq(r.assembly, false, 'not an assembly');
    const ds = BM.parseFile(r.fileName, r.text);
    eq(ds.format, 'mmcif', 'read as mmCIF'); eq(ds.title, 'Tiny test peptide', 'title');
    const a = await F.pdb('1ABC', { fetch: f, assembly: true });
    eq(f.calls[1].url, RCSB + '1ABC-assembly1.cif', 'assembly url'); eq(a.fileName, '1ABC-assembly1.cif', 'assembly file'); eq(a.name, '1ABC assembly 1', 'assembly name');
    const x = await F.pdb('pdb_00001ABC', { fetch: f });
    eq(f.calls[2].url, RCSB + 'pdb_00001abc.cif', 'extended ID url'); eq(x.fileName, 'pdb_00001abc.cif', 'extended ID file');
    await rejects(F.pdb('9ZZZ', { fetch: mockFetch(() => null) }), /no entry 9ZZZ/, '404');
    await rejects(F.pdb('9ZZZ', { fetch: mockFetch(() => null), assembly: true }), /no biological assembly/, '404 assembly');
    await rejects(F.pdb('1ABC', { fetch: mockFetch(() => ({ status: 503, body: '' })) }), /The PDB is busy/, '503');
    await rejects(F.pdb('1ABC', { fetch: mockFetch(() => ({ status: 200, body: '<html></html>' })) }), /not an mmCIF file/, 'not mmCIF');
  });

  await testAsync('time-out, cancelling, progress', async () => {
    const hang = mockFetch(() => 'hang');
    const t0 = Date.now();
    await rejects(F.pdb('1ABC', { fetch: hang, timeout: 60 }), /did not answer in time/, 'time-out');
    ok(Date.now() - t0 < 2000, 'the time-out is the one asked for');
    const ctl = new AbortController();
    const p = F.pubchem('aspirin', { fetch: mockFetch(() => 'hang'), signal: ctl.signal });
    setTimeout(() => ctl.abort(), 20);
    await rejects(p, /^cancelled$/, 'cancelled while waiting');
    const pre = new AbortController(); pre.abort();
    const idle = mockFetch(() => ({ status: 200, body: MMCIF }));
    await rejects(F.pdb('1ABC', { fetch: idle, signal: pre.signal }), /^cancelled$/, 'cancelled before');
    eq(idle.calls.length, 0, 'nothing sent once cancelled');
    // a body in pieces: the progress counts bytes, and a character split between pieces is decoded whole
    const bytes = Buffer.from(MMCIF.replace('Tiny test peptide', 'Tiny β-test peptide'), 'utf8');
    const cut = bytes.indexOf(Buffer.from('β', 'utf8')) + 1;
    const chunks = [bytes.subarray(0, cut), bytes.subarray(cut, cut + 100), bytes.subarray(cut + 100)].map(b => new Uint8Array(b));
    const seen = [];
    const st = mockFetch(() => ({ status: 200, body: 'unused', stream: streamOf(chunks) }));
    const r = await F.pdb('1ABC', { fetch: st, onProgress: n => seen.push(n) });
    eq(seen.join(','), [cut, cut + 100, bytes.length].join(','), 'progress in bytes');
    ok(r.text.indexOf('Tiny β-test peptide') > 0, 'UTF-8 across pieces');
    eq(BM.parseFile(r.fileName, r.text).title, 'Tiny β-test peptide', 'the streamed file reads');
  });

  console.log('BlueMolecule download tests: ' + pass + ' passed, ' + fail + ' failed');
  if (fail) { for (const f of failures.slice(0, 40)) console.log('  FAIL ' + f); process.exit(1); }
})();
