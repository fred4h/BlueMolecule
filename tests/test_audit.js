#!/usr/bin/env node
/* Regression tests for bugs found in the November 2026 audit.  Run:  node tests/test_audit.js */
'use strict';
const path = require('path');
const ROOT = path.join(__dirname, '..');
for (const m of ['elements', 'parsers', 'readers-struct', 'writers', 'readers-qc', 'geom', 'chem', 'writers-struct']) require(path.join(ROOT, 'src', m + '.js'));
const BM = globalThis.BM;

let pass = 0, fail = 0;
const failures = [];
function ok(cond, msg) { if (cond) pass++; else { fail++; failures.push(msg); } }
function throwsMsg(fn, re, msg) {
  try { fn(); ok(false, msg + ': no error thrown'); } catch (e) { ok(e instanceof Error && e.constructor === Error && re.test(e.message), msg + ': "' + e.message + '"'); }
}
const save = id => BM.saveFormats.find(f => f.id === id);
const ctx = (Z, xyz, extra) => Object.assign({ Z, xyz: Float64Array.from(xyz), bonds: [], cell: null, title: 't', meta: { charge: 0, multiplicity: 1 }, baseName: 't', frames: [], frameIndex: 0 }, extra || {});

// ---- readers: coordinates that overflow to Infinity are refused with a readable message (validate())
throwsMsg(() => BM.parseFile('big.xyz', '1\n\nH 1e999 0 0\n'), /not a finite number/, 'XYZ coordinate 1e999');
throwsMsg(() => BM.parseFile('big.xyz', '2\n\nH 0 0 0\nH 0 -1e400 0\n'), /atom 2/, 'non-finite coordinate names the atom');
ok(BM.parseFile('ok.xyz', '1\n\nH 1e300 0 0\n').frames[0].xyz[0] === 1e300, 'large but finite coordinates still read');
let ds = BM.parseFile('blowup.xyz', '1\na\nH 0 0 0\n1\nb\nH 0 0 1\n1\nc\nH 1e999 0 0\n1\nd\nH 0 0 3\n');
ok(ds.frames.length === 2 && ds.warnings.some(w => /Stopped at geometry 3/.test(w)), 'a trajectory keeps the geometries before the first non-finite one: ' + ds.frames.length + ' ' + ds.warnings.join('|'));

// ---- MOL/SDF: a blank line inside the atom block (truncated file) is an error, not a NaN atom
const sdf = 'acetic acid\n  BlueMol\n\n  2  1  0  0  0  0  0  0  0  0999 V2000\n   -1.3670    0.0570    0.0000 C   0  0\n';
throwsMsg(() => BM.parseFile('cut.sdf', sdf + '\n'), /truncated/, 'SDF atom block cut after a newline');
throwsMsg(() => BM.parseFile('cut.sdf', sdf), /truncated/, 'SDF atom block cut at the end of a line');

// ---- POSCAR: unknown species labels (the "X" BlueMolecule writes for dummy atoms) are opened as dummies
const pos = 'Si X\n 1.0\n 5.43 0 0\n 0 5.43 0\n 0 0 5.43\n Si X\n 1 1\nDirect\n 0 0 0\n 0.25 0.25 0.25\n';
ds = BM.parseFile('test.vasp', pos);
ok(ds.format === 'poscar' && JSON.stringify(ds.frames[0].Z) === '[14,0]', '.vasp with an "X" species opens (dummy atom): ' + ds.format + ' ' + JSON.stringify(ds.frames[0].Z));
ds = BM.parseFile('odd.vasp', pos.replace(/Si X\n 1 1/, 'Si Zz\n 1 1'));
ok(JSON.stringify(ds.frames[0].Z) === '[14,0]' && ds.warnings.some(w => /Unknown species/.test(w)), 'an unknown label ("Zz") opens as a dummy, with a warning');
const dummy = ctx([8, 1, 1, 0], [0, 0, 0.1173, 0, 0.7572, -0.4692, 0, -0.7572, -0.4692, 0, 0, 1.5]);
ds = BM.parseFile('dummy.vasp', save('poscar').write(dummy));
ok(ds.frames[0].Z.length === 4 && ds.frames[0].Z[3] === 0, 'POSCAR written with a dummy atom reads back');
ok(BM.detectFormat('notes.txt', 'just some words\n 1\n 2 3\n') !== 'poscar', 'loose text is not taken for a POSCAR');

// ---- TRR: a corrupt frame header (negative atom count) gives a readable error, not "Invalid typed array length"
function trr(nat, natField) {
  const i32 = v => { const b = Buffer.alloc(4); b.writeInt32BE(v); return b; };
  const f32 = v => { const b = Buffer.alloc(4); b.writeFloatBE(v); return b; };
  const parts = [i32(1993), i32(13), i32(12), Buffer.from('GMX_trn_file')];
  [0, 0, 36, 0, 0, 0, 0, 12 * nat, 0, 0, natField, 0, 0].forEach(v => parts.push(i32(v)));
  parts.push(f32(0), f32(0));
  for (let k = 0; k < 9; k++) parts.push(f32(k % 4 === 0 ? 3 : 0));
  for (let k = 0; k < 3 * nat; k++) parts.push(f32(0.1 * k));
  return new Uint8Array(Buffer.concat(parts));
}
ok(BM.parseFileBinary('ok.trr', trr(2, 2)).frames[0].xyz.length === 6, 'a small TRR still reads');
throwsMsg(() => BM.parseFileBinary('bad.trr', trr(2, -5)), /frame header is corrupt/, 'TRR with a negative atom count');

// ---- writers: free-format files keep a separator for very large coordinates; fixed-column formats refuse them
const far = ctx([6, 8], [-12345.678, 23456.789, -99999.5, 1, 2, 3]);
for (const id of ['xyz', 'gaussian-input', 'orca-input']) {
  const f = save(id), out = f.write(far), text = typeof out === 'string' ? out : out.files[0].text;
  let back = null;
  try { back = BM.parseFile('far.' + f.ext, text).frames[0]; } catch (e) { back = null; }
  ok(back && Math.abs(back.xyz[0] + 12345.678) < 1e-6 && Math.abs(back.xyz[2] + 99999.5) < 1e-6 && Math.abs(back.xyz[1] - 23456.789) < 1e-6,
    id + ' round-trips coordinates beyond -10000 A');
}
const xyzNormal = save('xyz').write(ctx([8], [1.5, -2.25, 0]));
ok(xyzNormal.split('\n')[2] === 'O       1.50000000    -2.25000000     0.00000000', 'ordinary XYZ columns unchanged: ' + JSON.stringify(xyzNormal.split('\n')[2]));
throwsMsg(() => save('pdb').write(ctx([6], [-1500, 0, 0])), /PDB coordinates must lie between/, 'PDB refuses x = -1500');
ok(/HETATM/.test(save('pdb').write(ctx([6], [-999.9, 9999.9, 0]))), 'PDB still writes the extreme representable values');
throwsMsg(() => save('mol').write(ctx([6], [0, 123456, 0])), /Molfile coordinates must lie between/, 'Molfile refuses y = 123456');

console.log('\nBlueMolecule audit regression tests: ' + pass + ' passed, ' + fail + ' failed');
if (fail) { failures.forEach(f => console.log('  FAIL: ' + f)); process.exit(1); }
