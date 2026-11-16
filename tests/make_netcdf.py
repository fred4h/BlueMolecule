#!/usr/bin/env python3
"""Writes the Amber NetCDF test files in tests/fixtures/netcdf with SciPy's NetCDF 3 writer (an implementation
independent of BlueMolecule's reader), following the Amber conventions (ambermd.org/netcdf/nctraj.xhtml).

    python3 tests/make_netcdf.py        (needs NumPy and SciPy; the files are committed, so the tests do not)

Files, all for the 4 atoms of tests/fixtures/struct/water.prmtop (O, H1, H2, Na+):
    traj_64bit.nc         12 frames, 64-bit offset format, time, box (changing), velocities between the variables
    traj_classic.nc       5 frames, classic format, no box, no time
    restart.ncrst         an AMBERRESTART file: double coordinates, velocities, box, time
    restart_netcdf.rst7   the same restart under the name pmemd gives it by default
    traj_cdf5.nc          6 frames in the 64-bit data format (CDF-5, which SciPy does not write: written here from the
                          format description with struct), with 64-bit and unsigned attributes
    expected.json         the coordinates, boxes and times as stored (float32 values written out exactly)
"""
import json
import shutil
import struct
import sys
from pathlib import Path

import numpy as np
from scipy.io import netcdf_file

OUT = Path(__file__).resolve().parent / 'fixtures' / 'netcdf'
N = 4
BASE = np.array([[0.000, 0.000, 0.117], [0.000, 0.757, -0.469], [0.000, -0.757, -0.469], [3.100, 0.400, 0.200]])


def frames(n, seed):
    """Coordinates that change from frame to frame (a vibrating water, a moving ion)."""
    out = []
    for k in range(n):
        t = 0.37 * k + seed
        d = np.array([[0.01 * np.sin(t), 0.02 * np.cos(t), 0.0],
                      [0.03 * np.sin(2 * t), -0.02 * np.sin(t), 0.01 * np.cos(t)],
                      [-0.03 * np.sin(2 * t), 0.02 * np.sin(t), 0.01 * np.cos(t)],
                      [0.25 * k, -0.1 * k, 0.05 * k]])
        out.append(BASE + d + np.array([10.0, 11.0, 12.0]))
    return np.array(out)


def amber_globals(f, conv):
    f.Conventions = conv
    f.ConventionVersion = '1.0'
    f.program = 'sander'
    f.programVersion = '18.0'
    f.title = 'BlueMolecule NetCDF test (written by SciPy)'
    f.application = 'AMBER'


def labels(f, box):
    f.createDimension('spatial', 3)
    v = f.createVariable('spatial', 'c', ('spatial',))
    v[:] = np.array(list('xyz'), dtype='S1')
    if box:
        f.createDimension('cell_spatial', 3)
        f.createDimension('cell_angular', 3)
        f.createDimension('label', 5)
        v = f.createVariable('cell_spatial', 'c', ('cell_spatial',))
        v[:] = np.array(list('abc'), dtype='S1')
        v = f.createVariable('cell_angular', 'c', ('cell_angular', 'label'))
        v[:] = np.array([list('alpha'), list('beta '), list('gamma')], dtype='S1')


def trajectory(path, version, n, box, time, vel):
    X = frames(n, 0.0 if version == 2 else 1.3)
    f = netcdf_file(str(path), 'w', version=version)
    amber_globals(f, 'AMBER')
    f.createDimension('frame', None)
    f.createDimension('atom', N)
    labels(f, box)
    exp = {'atoms': N, 'frames': n, 'xyz': [], 'cell': [], 'time': []}
    if time:
        t = f.createVariable('time', 'f', ('frame',))
        t.units = 'picosecond'
    if vel:   # a record variable between the time and the coordinates, which the reader must skip
        v = f.createVariable('velocities', 'f', ('frame', 'atom', 'spatial'))
        v.units = 'angstrom/picosecond'
        v.scale_factor = np.float32(20.455)
    c = f.createVariable('coordinates', 'f', ('frame', 'atom', 'spatial'))
    c.units = 'angstrom'
    if box:
        cl = f.createVariable('cell_lengths', 'd', ('frame', 'cell_spatial'))
        cl.units = 'angstrom'
        ca = f.createVariable('cell_angles', 'd', ('frame', 'cell_angular'))
        ca.units = 'degree'
    for k in range(n):
        c[k] = X[k].astype(np.float32)
        exp['xyz'].append([float(x) for x in X[k].astype(np.float32).ravel()])
        if time:
            t[k] = np.float32(2.0 * (k + 1))
            exp['time'].append(float(np.float32(2.0 * (k + 1))))
        if vel:
            v[k] = np.full((N, 3), 0.5 * k, dtype=np.float32)
        if box:
            L = [30.0 + 0.1 * k, 31.0 - 0.05 * k, 32.0]
            A = [90.0, 90.0, 90.0] if k % 2 == 0 else [90.0, 109.4712206, 90.0]
            cl[k] = L
            ca[k] = A
            exp['cell'].append(L + A)
    f.close()
    return exp


def restart(path):
    X = frames(1, 2.1)[0]
    f = netcdf_file(str(path), 'w', version=2)
    amber_globals(f, 'AMBERRESTART')
    f.createDimension('atom', N)
    labels(f, True)
    t = f.createVariable('time', 'd', ())
    t.units = 'picosecond'
    t[()] = 1234.5
    c = f.createVariable('coordinates', 'd', ('atom', 'spatial'))
    c.units = 'angstrom'
    c[:] = X
    v = f.createVariable('velocities', 'd', ('atom', 'spatial'))
    v.units = 'angstrom/picosecond'
    v.scale_factor = 20.455
    v[:] = np.full((N, 3), 0.25)
    cl = f.createVariable('cell_lengths', 'd', ('cell_spatial',))
    cl.units = 'angstrom'
    cl[:] = [40.0, 41.0, 42.0]
    ca = f.createVariable('cell_angles', 'd', ('cell_angular',))
    ca.units = 'degree'
    ca[:] = [90.0, 90.0, 90.0]
    f.close()
    return {'atoms': N, 'frames': 1, 'xyz': [[float(x) for x in X.ravel()]], 'cell': [[40.0, 41.0, 42.0, 90.0, 90.0, 90.0]], 'time': [1234.5]}


def cdf5(path):
    """A CDF-5 trajectory written byte by byte (Unidata/PnetCDF format description): 64-bit counts, sizes and offsets."""
    n = 6
    X = frames(n, 0.7).astype(np.float32)
    T = np.array([0.5 * (k + 1) for k in range(n)], dtype=np.float32)
    i32 = lambda v: struct.pack('>i', v)
    i64 = lambda v: struct.pack('>q', v)
    pad = lambda b: b + b'\0' * (-len(b) % 4)
    name = lambda t: i64(len(t)) + pad(t.encode())
    def att(nm, typ, vals):
        if typ == 2:
            body, cnt = vals.encode(), len(vals)
        else:
            fmt = {4: 'i', 5: 'f', 6: 'd', 9: 'I', 10: 'q', 11: 'Q'}[typ]
            body, cnt = struct.pack('>' + fmt * len(vals), *vals), len(vals)
        return name(nm) + i32(typ) + i64(cnt) + pad(body)
    def atts(lst):
        return i32(12) + i64(len(lst)) + b''.join(att(*a) for a in lst) if lst else i32(0) + i64(0)
    dims = [('frame', 0), ('spatial', 3), ('atom', N)]
    gatts = [('Conventions', 2, 'AMBER'), ('ConventionVersion', 2, '1.0'), ('program', 2, 'cdf5-writer'), ('big_number', 10, [2 ** 40 + 5]),
             ('unsigned', 9, [4000000000]), ('ratio', 6, [0.25])]
    # variables: spatial (non-record char), time (record float), coordinates (record float)
    vars_ = [('spatial', [1], [], 2, 3), ('time', [0], [('units', 2, 'picosecond')], 5, 4), ('coordinates', [0, 2, 1], [('units', 2, 'angstrom')], 5, N * 3 * 4)]
    def header(begins):
        h = b'CDF\x05' + i64(n) + i32(10) + i64(len(dims)) + b''.join(name(d) + i64(l) for d, l in dims) + atts(gatts)
        h += i32(11) + i64(len(vars_))
        for (nm, ids, va, typ, size), b in zip(vars_, begins):
            h += name(nm) + i64(len(ids)) + b''.join(i64(i) for i in ids) + atts(va) + i32(typ) + i64(size + (-size % 4)) + i64(b)
        return h
    hlen = len(header([0, 0, 0]))
    spatial_at = hlen
    rec0 = spatial_at + 4                           # the 3 chars, padded to 4
    begins = [spatial_at, rec0, rec0 + 4]           # time then coordinates in each record
    rec_size = 4 + N * 3 * 4
    data = bytearray(header(begins))
    data += pad(b'xyz')
    for k in range(n):
        data += struct.pack('>f', T[k]) + struct.pack('>' + 'f' * (3 * N), *X[k].ravel())
    assert len(data) == rec0 + n * rec_size
    Path(path).write_bytes(bytes(data))
    return {'atoms': N, 'frames': n, 'xyz': [[float(x) for x in X[k].ravel()] for k in range(n)], 'cell': [], 'time': [float(t) for t in T]}


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    exp = {
        'traj_64bit.nc': trajectory(OUT / 'traj_64bit.nc', 2, 12, True, True, True),
        'traj_classic.nc': trajectory(OUT / 'traj_classic.nc', 1, 5, False, False, False),
        'restart.ncrst': restart(OUT / 'restart.ncrst'),
    }
    shutil.copyfile(OUT / 'restart.ncrst', OUT / 'restart_netcdf.rst7')
    exp['traj_cdf5.nc'] = cdf5(OUT / 'traj_cdf5.nc')
    exp['restart_netcdf.rst7'] = exp['restart.ncrst']
    (OUT / 'expected.json').write_text(json.dumps(exp, indent=1) + '\n')
    for name in sorted(exp):
        print('%-22s %6d bytes' % (name, (OUT / name).stat().st_size))
    return 0


if __name__ == '__main__':
    sys.exit(main())
