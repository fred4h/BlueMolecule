#!/usr/bin/env python3
"""Synthetic test outputs for BlueMolecule's large-file and QM/MM readers.

  python3 tools/make_big_fixture.py gaussian OUT.log [--mb 600]   Gaussian 16 optimisation log of about N MB
  python3 tools/make_big_fixture.py orca OUT.out [--mb 600]       ORCA 5 optimisation output of about N MB
  python3 tools/make_big_fixture.py qmmm DIR                      small QM/MM / ONIOM outputs (tests/fixtures/qmmm)

The big logs repeat realistic optimisation cycles (orientation tables, SCF iterations, full MO coefficient prints,
population analysis, forces, convergence tables); geometry and energy change every cycle. OUT.json beside the log
records what a reader must find: number of geometries, atoms, the energy series and the convergence flags.
The layouts follow the public appearance of Gaussian / ORCA outputs; no program code is involved.
"""
import json
import math
import os
import sys

# ----------------------------------------------------------------------------- molecules
SYM = {1: 'H', 6: 'C', 7: 'N', 8: 'O'}
NBF = {1: 2, 6: 15, 7: 15, 8: 15}          # 6-31G(d), 6D


def chain_molecule(n_heavy):
    """A zig-zag carbon chain with an oxygen every 5th position, hydrogens on carbons: (Z list, xyz list)."""
    Z, X = [], []
    for k in range(n_heavy):
        z = 8 if k % 5 == 4 else 6
        x, y = 1.25 * k, 0.75 * (k % 2)
        Z.append(z); X.append([x, y, 0.0])
    for k in range(n_heavy):
        if Z[k] != 6:
            continue
        x, y, _ = X[k]
        dy = 0.63 if k % 2 else -0.63
        Z.append(1); X.append([x, y + dy, 0.89])
        Z.append(1); X.append([x, y + dy, -0.89])
    return Z, X


def perturbed(X, k, amp, tau):
    f = amp * math.exp(-k / tau)
    return [[c + f * math.sin(1.7 * i + 0.9 * d + 0.37 * k) for d, c in enumerate(p)] for i, p in enumerate(X)]


def energy_series(n, e_final, de, tau):
    return [e_final + de * math.exp(-k / tau) * (1 + 0.05 * math.cos(k)) for k in range(n - 1)] + [e_final]


# ----------------------------------------------------------------------------- Gaussian pieces
GLINE = ' ' + '-' * 69 + '\n'


def g_orient(kind, Z, X):
    out = ['                          %s orientation:                          \n' % kind, GLINE,
           ' Center     Atomic      Atomic             Coordinates (Angstroms)\n',
           ' Number     Number       Type             X           Y           Z\n', GLINE]
    for i, (z, p) in enumerate(zip(Z, X)):
        out.append('%7d%11d%12d    %12.6f%12.6f%12.6f\n' % (i + 1, z, 0, p[0], p[1], p[2]))
    out.append(GLINE)
    return ''.join(out)


def g_scf_cycles(n, e0, e1, method='RB3LYP'):
    out = [' Closed shell SCF:\n', ' Using DIIS extrapolation, IDIIS=  1040.\n',
           ' Integral symmetry usage will be decided dynamically.\n', ' Keep R1 ints in memory in canonical form, NReq=1234567.\n']
    for c in range(1, n + 1):
        e = e1 + (e0 - e1) * math.exp(-0.8 * c)
        err = 0.05 * math.exp(-0.9 * c)
        out.append(' Cycle %3d  Pass 1  IDiag  1:\n' % c)
        out.append(' E= %20.12f     Delta-E=       %15.12f Rises=F Damp=F\n' % (e, -err))
        out.append(' DIIS: error= %8.2E at cycle %3d NSaved= %3d.\n' % (err, c, min(c, 20)))
        out.append(' NSaved=%3d IEnMin=%3d EnMin= %20.12f     IErMin=%3d ErrMin= %8.2E\n' % (min(c, 20), c, e, c, err))
        out.append(' ErrMax= %8.2E  0.00D+00 EMaxC= 1.00D-01 BMatC= %8.2E BMatP= %8.2E\n' % (err, err * err, err * err))
        out.append(' IDIUse=1 WtCom= 1.00D+00 WtEn= 0.00D+00\n')
        out.append(' Coeff-Com: ' + ''.join(' %10.3E' % (0.1 * (j + 1)) for j in range(min(c, 6))) + '\n')
        out.append(' Coeff:     ' + ''.join(' %10.3E' % (0.1 * (j + 1)) for j in range(min(c, 6))) + '\n')
        out.append(' Gap=     0.234 Goal=   None    Shift=    0.000\n')
        out.append(' RMSDP=%8.2E MaxDP=%8.2E DE=%10.2E OVMax= %8.2E\n\n' % (err / 30, err / 3, -err, err))
    out.append(' SCF Done:  E(%s) =  %.12f     A.U. after %4d cycles\n' % (method, e1, n))
    out.append('            NFock=%3d  Conv=0.43D-08     -V/T= 2.0098\n' % n)
    return ''.join(out)


_MO_CACHE = {}


def g_mo_block(Z, nocc, variant):
    key = (len(Z), nocc, variant)
    if key in _MO_CACHE:
        return _MO_CACHE[key]
    labels = []
    for a, z in enumerate(Z):
        names = ['1S', '2S', '2PX', '2PY', '2PZ', '3S', '3PX', '3PY', '3PZ', '4XX', '4YY', '4ZZ', '4XY', '4XZ', '4YZ'] if z > 1 else ['1S', '2S']
        for k, nm in enumerate(names):
            labels.append((a + 1, SYM[z], nm, k == 0))
    nbf = len(labels)
    eig = [-20.0 + 20.3 * (i / nbf) ** 0.6 + 0.001 * variant for i in range(nbf)]
    out = ['     Molecular Orbital Coefficients:\n']
    for b0 in range(0, nbf, 5):
        cols = range(b0, min(nbf, b0 + 5))
        out.append('                       ' + ''.join('%10d' % (c + 1) for c in cols) + '\n')
        out.append('                       ' + ''.join('%10s' % ('O' if c < nocc else 'V') for c in cols) + '\n')
        out.append('     Eigenvalues -- ' + ''.join('%10.5f' % eig[c] for c in cols) + '\n')
        for r, (a, el, nm, first) in enumerate(labels):
            head = '%4d %-3d %-2s %-6s' % (r + 1, a, el, nm) if first else '%4d        %-6s' % (r + 1, nm)
            vals = ''.join('%10.5f' % (0.9 * math.sin(0.37 * r * (c + 1) + 0.11 * variant)) for c in cols)
            out.append(head.ljust(21) + vals + '\n')
    s = ''.join(out)
    _MO_CACHE[key] = s
    return s


def g_population(Z, nocc, variant, charges):
    nbf = sum(NBF[z] for z in Z)
    eig = [-20.0 + 20.3 * (i / nbf) ** 0.6 for i in range(nbf)]
    out = [' Copying SCF densities to generalized density rwf, IOpCl= 0 IROHF=0.\n\n',
           ' **********************************************************************\n\n',
           '            Population analysis using the SCF Density.\n\n',
           ' **********************************************************************\n\n',
           ' Orbital symmetries:\n       Occupied  ' + ' '.join(['(A)'] * min(nocc, 15)) + '\n',
           ' The electronic state is 1-A.\n']
    occ, virt = eig[:nocc], eig[nocc:]
    for k in range(0, len(occ), 5):
        out.append(' Alpha  occ. eigenvalues -- ' + ''.join('%10.5f' % v for v in occ[k:k + 5]) + '\n')
    for k in range(0, len(virt), 5):
        out.append(' Alpha virt. eigenvalues -- ' + ''.join('%10.5f' % v for v in virt[k:k + 5]) + '\n')
    out.append(g_mo_block(Z, nocc, variant))
    out.append('     Density Matrix:\n')
    n = min(nbf, 40)
    for b0 in range(0, n, 5):
        out.append('                       ' + ''.join('%10d' % (c + 1) for c in range(b0, min(n, b0 + 5))) + '\n')
        for r in range(b0, n):
            out.append('%4d %-15s' % (r + 1, '') + ''.join('%10.5f' % (0.01 * ((r * 7 + c) % 13)) for c in range(b0, min(r + 1, b0 + 5))) + '\n')
    out.append('     Condensed to atoms (all electrons):\n')
    N = len(Z)
    for b0 in range(0, N, 6):
        out.append('            ' + ''.join('%11d' % (c + 1) for c in range(b0, min(N, b0 + 6))) + '\n')
        for r in range(N):
            out.append('%6d  %-2s' % (r + 1, SYM[Z[r]]) + ''.join('%11.6f' % (0.001 * ((r + c) % 17)) for c in range(b0, min(N, b0 + 6))) + '\n')
    out.append(' Mulliken charges:\n               1\n')
    for i, z in enumerate(Z):
        out.append('%6d  %-2s%12.6f\n' % (i + 1, SYM[z], charges[i]))
    out.append(' Sum of Mulliken charges =   0.00000\n')
    out.append(' Dipole moment (field-independent basis, Debye):\n')
    out.append('    X=              1.2345    Y=             -0.1234    Z=              0.0000  Tot=              1.2407\n')
    return ''.join(out)


def g_forces(Z, k):
    out = [' ***** Axes restored to original set *****\n', ' ' + '-' * 67 + '\n',
           ' Center     Atomic                   Forces (Hartrees/Bohr)\n',
           ' Number     Number              X              Y              Z\n', ' ' + '-' * 67 + '\n']
    s = math.exp(-k / 6.0)
    for i, z in enumerate(Z):
        out.append('%7d%9d       %15.9f%15.9f%15.9f\n' % (i + 1, z, 0.01 * s * math.sin(i), -0.01 * s * math.cos(i), 0.002 * s))
    out.append(' ' + '-' * 67 + '\n')
    out.append(' Cartesian Forces:  Max     %.9f RMS     %.9f\n' % (0.01 * s, 0.004 * s))
    return ''.join(out)


def g_conv(step, maxstep, k, n_steps):
    s = math.exp(-k / 6.0)
    last = k == n_steps - 1
    vals = [(0.02 * s, 0.000450), (0.008 * s, 0.000300), (0.05 * s, 0.001800), (0.02 * s, 0.001200)]
    if last:
        vals = [(0.000012, 0.000450), (0.000004, 0.000300), (0.000150, 0.001800), (0.000060, 0.001200)]
    names = ['Maximum Force', 'RMS     Force', 'Maximum Displacement', 'RMS     Displacement']
    out = [' GradGradGradGradGradGradGradGradGradGradGradGradGradGradGradGradGradGrad\n', ' Berny optimization.\n',
           ' Using GEDIIS/GDIIS optimizer.\n', ' Internal  Forces:  Max     %.9f RMS     %.9f\n' % (vals[0][0], vals[1][0]),
           ' Search for a local minimum.\n', ' Step number %3d out of a maximum of %4d\n' % (step, maxstep),
           ' All quantities printed in internal units (Hartrees-Bohrs-Radians).\n',
           '         Item               Value     Threshold  Converged?\n']
    for nm, (v, t) in zip(names, vals):
        out.append(' %-20s %10.6f %12.6f     %s\n' % (nm, v, t, 'YES' if v < t else 'NO '))
    out.append(' Predicted change in Energy=-%.6fD-04\n' % (1.234 * s))
    if last:
        out.append(' Optimization completed.\n    -- Stationary point found.\n')
    out.append(' GradGradGradGradGradGradGradGradGradGradGradGradGradGradGradGradGradGrad\n')
    return ''.join(out)


def gaussian_big(path, mb):
    Z, X0 = chain_molecule(24)              # 24 heavy atoms + 40 H = 64 atoms
    nocc = sum(Z) // 2
    charges = [(-0.35 if z == 8 else -0.12 if z == 6 else 0.11) for z in Z]

    def block(k, n_steps, e):
        X = perturbed(X0, k, 0.08, n_steps / 5.0)
        return ''.join([' (Enter /opt/g16/l202.exe)\n', g_orient('Input', Z, X),
                        ' Stoichiometry    C20H40O4\n Framework group  C1[X(C20H40O4)]\n Deg. of freedom   186\n',
                        ' Full point group                 C1      NOp   1\n', g_orient('Standard', Z, X),
                        ' Rotational constants (GHZ):           0.4567890           0.0567890           0.0512345\n',
                        ' Leave Link  202 at Sun Nov  9 12:00:00 2026, MaxMem=  1073741824 cpu:               0.1 elap:               0.0\n',
                        ' (Enter /opt/g16/l301.exe)\n Standard basis: 6-31G(d) (6D, 7F)\n',
                        ' There are   %d symmetry adapted cartesian basis functions of A   symmetry.\n' % sum(NBF[z] for z in Z),
                        '    %d alpha electrons       %d beta electrons\n' % (nocc, nocc),
                        '       nuclear repulsion energy      %14.10f Hartrees.\n' % (1234.5 + k * 0.01),
                        ' (Enter /opt/g16/l502.exe)\n', g_scf_cycles(14, e + 0.6, e),
                        ' (Enter /opt/g16/l601.exe)\n', g_population(Z, nocc, k % 3, charges),
                        ' (Enter /opt/g16/l716.exe)\n', g_forces(Z, k), ' (Enter /opt/g16/l103.exe)\n\n',
                        g_conv(k + 1, n_steps + 10, k, n_steps)])

    n_steps = max(3, int(round(mb * 1e6 / len(block(1, 100, -925.0)))))
    E = energy_series(n_steps, -925.4183124, 0.035, n_steps / 5.0)
    route = '#P B3LYP/6-31G(d) Opt(MaxCycles=%d) Pop=Full' % (n_steps + 10)
    head = [' Entering Gaussian System, Link 0=g16\n', ' Input=big.gjf\n', ' Output=big.log\n',
            ' Copyright (c) 1988-2019, Gaussian, Inc.  All Rights Reserved.\n', ' This is part of the Gaussian(R) 16 program.\n',
            ' ******************************************\n', ' Gaussian 16:  ES64L-G16RevC.01  3-Jul-2019\n',
            '                 9-Nov-2026 \n', ' ******************************************\n',
            ' %chk=big.chk\n', ' %mem=8GB\n', ' ' + '-' * (len(route) + 1) + '\n', ' ' + route + '\n', ' ' + '-' * (len(route) + 1) + '\n',
            ' 1/14=-1,18=20,19=15,38=1/1,3;\n', ' 2/9=110,12=2,17=6,18=5,40=1/2;\n', ' -------------------\n',
            ' Big synthetic chain\n', ' -------------------\n', ' Symbolic Z-matrix:\n', ' Charge =  0 Multiplicity = 1\n']
    for z, p in zip(Z, X0):
        head.append(' %-24s%14.8f%14.8f%14.8f\n' % (SYM[z], p[0], p[1], p[2]))
    head.append(' \n GradGradGradGradGradGradGradGradGradGradGradGradGradGradGradGradGradGrad\n Berny optimization.\n Initialization pass.\n')
    with open(path, 'w') as f:
        f.write(''.join(head))
        for k in range(n_steps):
            f.write(block(k, n_steps, E[k]))
        f.write(' Normal termination of Gaussian 16 at Sun Nov  9 12:34:56 2026.\n')
    meta = {'program': 'gaussian', 'frames': n_steps, 'atoms': len(Z), 'energies': E,
            'converged_last': True, 'bytes': os.path.getsize(path)}
    with open(os.path.splitext(path)[0] + '.json', 'w') as f:
        json.dump(meta, f)
    return meta


# ----------------------------------------------------------------------------- ORCA pieces
def o_coords(Z, X, au=True, internal=True):
    out = ['---------------------------------\n', 'CARTESIAN COORDINATES (ANGSTROEM)\n', '---------------------------------\n']
    for z, p in zip(Z, X):
        out.append('  %-2s  %12.6f%12.6f%12.6f\n' % (SYM.get(z, 'H'), p[0], p[1], p[2]))
    out.append('\n')
    if au:
        out += ['----------------------------\n', 'CARTESIAN COORDINATES (A.U.)\n', '----------------------------\n',
                '  NO LB      ZA    FRAG     MASS         X           Y           Z\n']
        for i, (z, p) in enumerate(zip(Z, X)):
            out.append('%4d %-2s  %8.4f %4d %9.3f %11.6f %11.6f %11.6f\n' % (i, SYM.get(z, 'H'), z, 0, 2.0 * z, p[0] / 0.529177, p[1] / 0.529177, p[2] / 0.529177))
        out.append('\n')
    if internal:
        out += ['--------------------------------\n', 'INTERNAL COORDINATES (ANGSTROEM)\n', '--------------------------------\n']
        for i, z in enumerate(Z):
            out.append(' %-2s  %4d %3d %3d     %14.12f   %12.8f   %12.8f\n' % (SYM.get(z, 'H'), max(0, i), max(0, i - 1), max(0, i - 2), 1.5 if i else 0.0, 109.5 if i > 1 else 0.0, 180.0 if i > 2 else 0.0))
        out.append('\n')
    return ''.join(out)


def o_scf(n, e1):
    out = ['--------------\n', 'SCF ITERATIONS\n', '--------------\n',
           'ITER       Energy         Delta-E        Max-DP      RMS-DP      [F,P]     Damp\n']
    prev = e1 + 0.8
    for c in range(n):
        e = e1 + 0.8 * math.exp(-0.9 * c)
        out.append('%3d %17.10f %17.12f %10.8f %11.8f %10.7f %6.4f\n' % (c, e, e - prev, 0.04 * math.exp(-c), 0.002 * math.exp(-c), 0.1 * math.exp(-c), 0.7))
        if c == 0:
            out.append('                               ***Turning on DIIS***\n')
        if c == n // 2:
            out.append('               *** Restarting incremental Fock matrix formation ***\n')
        prev = e
    out += ['\n', '                  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n',
            '                  !        SCF CONVERGED AFTER %3d CYCLES        !\n' % n,
            '                  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n\n',
            '----------------\n', 'TOTAL SCF ENERGY\n', '----------------\n\n',
            'Total Energy       :     %17.8f Eh          %14.5f eV\n\n' % (e1, e1 * 27.211386),
            'Components:\n', 'Nuclear Repulsion  :       %17.8f Eh          %14.5f eV\n' % (1234.5, 1234.5 * 27.2114)]
    return ''.join(out)


_OMO = {}


def o_mos(Z, nocc, variant):
    key = (len(Z), nocc, variant)
    if key in _OMO:
        return _OMO[key]
    labels = []
    for a, z in enumerate(Z):
        names = ['s', 's', 'pz', 'px', 'py', 's', 'pz', 'px', 'py', 'dz2', 'dxz', 'dyz', 'dx2y2', 'dxy'] if z > 1 else ['s', 's', 'pz', 'px', 'py']
        for nm in names:
            labels.append(('%3d%-2s' % (a, SYM[z]), nm))
    nbf = len(labels)
    out = ['------------------\n', 'MOLECULAR ORBITALS\n', '------------------\n']
    for b0 in range(0, nbf, 6):
        cols = range(b0, min(nbf, b0 + 6))
        out.append('                  ' + ''.join('%10d' % c for c in cols) + '   \n')
        out.append('                  ' + ''.join('%10.5f' % (-20 + 20.3 * (c / nbf) ** 0.6) for c in cols) + '\n')
        out.append('                  ' + ''.join('%10.5f' % (2.0 if c < nocc else 0.0) for c in cols) + '\n')
        out.append('                  ' + '  --------' * len(cols) + '\n')
        for r, (at, nm) in enumerate(labels):
            out.append('%s %-6s' % (at, nm) + ''.join('%10.6f' % (0.9 * math.sin(0.29 * r * (c + 1) + 0.13 * variant)) for c in cols) + '\n')
    s = ''.join(out)
    _OMO[key] = s
    return s


def o_pop(Z, charges):
    out = ['\n', '                    ********************************\n', '                    * MULLIKEN POPULATION ANALYSIS *\n',
           '                    ********************************\n\n', '-----------------------\n', 'MULLIKEN ATOMIC CHARGES\n', '-----------------------\n']
    for i, z in enumerate(Z):
        out.append('%4d %-2s:%12.6f\n' % (i, SYM.get(z, 'H'), charges[i]))
    out.append('Sum of atomic charges:    0.0000000\n\n')
    out += ['--------------------------------\n', 'MULLIKEN REDUCED ORBITAL CHARGES\n', '--------------------------------\n']
    for i, z in enumerate(Z):
        out.append('%3d %-2s s       :     %8.6f  s :     %8.6f\n' % (i, SYM.get(z, 'H'), 3.1, 3.1))
        if z > 1:
            out.append('      pz      :     1.012345  p :     3.012345\n      px      :     1.000000\n      py      :     1.000000\n')
            out.append('      dz2     :     0.001234  d :     0.012345\n')
    out += ['\n', '                     *****************************\n', '                     * MAYER POPULATION ANALYSIS *\n',
            '                     *****************************\n\n']
    for i in range(0, len(Z) - 2, 3):
        out.append('B(%3d-%-2s,%3d-%-2s) :   0.9576 B(%3d-%-2s,%3d-%-2s) :   0.9576 B(%3d-%-2s,%3d-%-2s) :   0.9576\n'
                   % (i, SYM.get(Z[i], 'H'), i + 1, SYM.get(Z[i + 1], 'H'), i, SYM.get(Z[i], 'H'), i + 2, SYM.get(Z[i + 2], 'H'), i + 1, SYM.get(Z[i + 1], 'H'), i + 2, SYM.get(Z[i + 2], 'H')))
    return ''.join(out)


def o_grad(Z, k):
    s = math.exp(-k / 6.0)
    out = ['------------------\n', 'CARTESIAN GRADIENT\n', '------------------\n\n']
    for i, z in enumerate(Z):
        out.append('%4d   %-2s  :  %14.9f %14.9f %14.9f\n' % (i + 1, SYM.get(z, 'H'), 0.01 * s * math.sin(i), -0.01 * s * math.cos(i), 0.002 * s))
    out.append('\nDifference to translation invariance:\n           :   0.0000012   -0.0000003    0.0000000\n\n')
    return ''.join(out)


def o_conv(k, n, dE):
    s = math.exp(-k / 6.0)
    last = k == n - 1
    rows = [('RMS gradient', 0.004 * s, 0.0001), ('MAX gradient', 0.01 * s, 0.0003), ('RMS step', 0.01 * s, 0.002), ('MAX step', 0.04 * s, 0.004)]
    if last:
        rows = [(a, t * 0.3, t) for a, _, t in rows]
    out = ['                                .--------------------.\n',
           '          ----------------------|Geometry convergence|-------------------------\n',
           '          Item                value                   Tolerance       Converged\n',
           '          ---------------------------------------------------------------------\n']
    if k > 0:
        out.append('          %-18s%15.10f %18.10f      %s\n' % ('Energy change', dE, 0.000005, 'YES' if abs(dE) < 5e-6 else 'NO'))
    for nm, v, t in rows:
        out.append('          %-18s%15.10f %18.10f      %s\n' % (nm, v, t, 'YES' if v < t else 'NO'))
    out += ['          ........................................................\n',
            '          Max(Bonds)      0.0123      Max(Angles)    1.23\n',
            '          ---------------------------------------------------------------------\n']
    return ''.join(out)


def orca_header(inp_lines, version='5.0.4'):
    out = ['\n                                 *****************\n', '                                 * O   R   C   A *\n',
           '                                 *****************\n\n',
           '                           Program Version %s -  RELEASE  -\n\n' % version,
           '================================================================================\n',
           '                                       INPUT FILE\n',
           '================================================================================\n', 'NAME = job.inp\n']
    for i, s in enumerate(inp_lines):
        out.append('|%3d> %s\n' % (i + 1, s))
    out.append('|%3d>                          ****END OF INPUT****\n' % (len(inp_lines) + 1))
    out.append('================================================================================\n\n')
    return ''.join(out)


def o_cycle_head(c):
    return ('\n                         *************************************************************\n'
            '                         *                GEOMETRY OPTIMIZATION CYCLE %3d            *\n'
            '                         *************************************************************\n' % c)


def orca_big(path, mb):
    Z, X0 = chain_molecule(24)
    nocc = sum(Z) // 2
    charges = [(-0.35 if z == 8 else -0.12 if z == 6 else 0.11) for z in Z]
    nbf = sum(14 if z > 1 else 5 for z in Z)

    def block(k, n_steps, E):
        X = perturbed(X0, k, 0.08, n_steps / 5.0)
        blk = [o_cycle_head(k + 1), o_coords(Z, X), o_scf(16, E[k]), '\n----------------\nORBITAL ENERGIES\n----------------\n\n  NO   OCC          E(Eh)            E(eV) \n']
        for i in range(nbf):
            e = -20 + 20.3 * (i / nbf) ** 0.6
            blk.append('%4d   %6.4f   %14.6f   %14.4f \n' % (i, 2.0 if i < nocc else 0.0, e, e * 27.2114))
        blk += ['\n', o_mos(Z, nocc, k % 3), o_pop(Z, charges),
                '\n-------------------------   --------------------\nFINAL SINGLE POINT ENERGY      %18.12f\n-------------------------   --------------------\n\n' % E[k],
                o_grad(Z, k), o_conv(k, n_steps, E[k] - E[k - 1] if k else 0.0)]
        if k == n_steps - 1:
            blk.append('\n                    ***********************HURRAY********************\n'
                       '                    ***        THE OPTIMIZATION HAS CONVERGED     ***\n'
                       '                    *************************************************\n\n'
                       '                       *******************************************************\n'
                       '                       *** FINAL ENERGY EVALUATION AT THE STATIONARY POINT ***\n'
                       '                       *******************************************************\n\n')
            blk.append(o_coords(Z, X, au=False, internal=False))
            blk.append(o_scf(4, E[k]))
            blk.append('\n-------------------------   --------------------\nFINAL SINGLE POINT ENERGY      %18.12f\n-------------------------   --------------------\n\n' % E[k])
        return ''.join(blk)

    n_steps = max(3, int(round(mb * 1e6 / len(block(1, 100, [-925.0] * 3)))))
    E = energy_series(n_steps, -925.3812871, 0.035, n_steps / 5.0)
    inp = ['! B3LYP def2-SVP Opt TightSCF', '%output', '  Print[P_MOs] 1', 'end', '%%geom MaxIter %d end' % (n_steps + 20), '* xyz 0 1']
    inp += ['  %-2s %12.6f %12.6f %12.6f' % (SYM[z], p[0], p[1], p[2]) for z, p in zip(Z, X0)] + ['*', '']
    with open(path, 'w') as f:
        f.write(orca_header(inp))
        f.write('                       *****************************\n                       * Geometry Optimization Run *\n'
                '                       *****************************\n\n')
        for k in range(n_steps):
            f.write(block(k, n_steps, E))
        f.write('\n                             ****ORCA TERMINATED NORMALLY****\nTOTAL RUN TIME: 0 days 1 hours 2 minutes 3 seconds 456 msec\n')
    meta = {'program': 'orca', 'frames': n_steps, 'atoms': len(Z), 'energies': E, 'converged_last': True, 'bytes': os.path.getsize(path)}
    with open(os.path.splitext(path)[0] + '.json', 'w') as f:
        json.dump(meta, f)
    return meta


# ----------------------------------------------------------------------------- QM/MM and ONIOM fixtures
# Propanoic acid + two waters (17 atoms). QM / high layer: CH2-COOH (atoms 5-11, 1-based); the C5-C1 bond is cut
# (link H on C1); the second water (15-17) is frozen.
ACID = [  # el, MM type, charge, layer, freeze, x, y, z
    ('C', 'CT', -0.0932, 'L', 0, 0.000, 0.000, 0.000), ('H', 'HC', 0.0400, 'L', 0, -0.363, -1.027, 0.000),
    ('H', 'HC', 0.0400, 'L', 0, -0.363, 0.513, 0.890), ('H', 'HC', 0.0400, 'L', 0, -0.363, 0.513, -0.890),
    ('C', 'CT', -0.0500, 'H', 0, 1.530, 0.000, 0.000), ('H', 'HC', 0.0500, 'H', 0, 1.893, -0.513, 0.890),
    ('H', 'HC', 0.0500, 'H', 0, 1.893, -0.513, -0.890), ('C', 'C', 0.7000, 'H', 0, 2.040, 1.420, 0.000),
    ('O', 'O', -0.5800, 'H', 0, 1.330, 2.410, 0.000), ('O', 'OH', -0.6100, 'H', 0, 3.370, 1.520, 0.000),
    ('H', 'HO', 0.4400, 'H', 0, 3.640, 2.460, 0.000), ('O', 'OW', -0.8340, 'L', 0, 0.400, 4.700, 0.600),
    ('H', 'HW', 0.4170, 'L', 0, 0.700, 3.800, 0.400), ('H', 'HW', 0.4170, 'L', 0, -0.560, 4.640, 0.640),
    ('O', 'OW', -0.8340, 'L', -1, 5.100, 3.800, -0.200), ('H', 'HW', 0.4170, 'L', -1, 4.350, 3.180, -0.150),
    ('H', 'HW', 0.4170, 'L', -1, 5.800, 3.150, -0.300)]
ACID_LINK = (0, 4)          # 0-based: atom replaced by the link H, atom it is bonded to
GFAC = 0.723886
ZOF = {'H': 1, 'C': 6, 'N': 7, 'O': 8}


def qm_steps(X0, frozen, n, amp=0.035):
    out = []
    for k in range(n):
        f = amp * math.exp(-0.9 * k) if k < n - 1 else 0.0
        out.append([[c + (0.0 if frozen[i] else f * math.sin(1.3 * i + 0.7 * d + 0.5 * k)) for d, c in enumerate(p)] for i, p in enumerate(X0)])
    return out


def link_pos(X, host, qm):
    return [X[qm][d] + GFAC * (X[host][d] - X[qm][d]) for d in range(3)]


def g_header(route, title, link0=('%chk=job.chk', '%mem=4GB')):
    out = [' Entering Gaussian System, Link 0=g16\n', ' Input=job.gjf\n', ' Output=job.log\n',
           ' Copyright (c) 1988-2019, Gaussian, Inc.  All Rights Reserved.\n', ' This is part of the Gaussian(R) 16 program.\n',
           ' ******************************************\n', ' Gaussian 16:  ES64L-G16RevC.01  3-Jul-2019\n',
           '                 9-Nov-2026 \n', ' ******************************************\n']
    out += [' %s\n' % l for l in link0]
    out += [' ' + '-' * (len(route) + 1) + '\n', ' ' + route + '\n', ' ' + '-' * (len(route) + 1) + '\n',
            ' 1/14=-1,18=20,19=15,38=1/1,3;\n', ' 2/9=110,12=2,17=6,18=5,40=1/2;\n',
            ' ' + '-' * (len(title) + 1) + '\n', ' ' + title + '\n', ' ' + '-' * (len(title) + 1) + '\n']
    return ''.join(out)


def g_mm_energy(e):
    return ' Energy= %20.12f NIter=   0.\n Dipole moment=      0.123456      -0.234567       0.012345\n' % e


def g_scf_done(method, e, cyc=11):
    return (' Standard basis: 6-31G(d) (6D, 7F)\n Closed shell SCF:\n Requested convergence on RMS density matrix=1.00D-08 within 128 cycles.\n'
            ' SCF Done:  E(%s) =  %.12f     A.U. after %3d cycles\n            NFock=%3d  Conv=0.52D-08     -V/T= 2.0087\n' % (method, e, cyc, cyc))


def write_g16_oniom_amber(path, n=6):
    """ONIOM(B3LYP/6-31G(d):Amber)=EmbedCharge optimisation; every step prints three sub-calculations (MM model,
    QM model with its own orientation table of the model system incl. the link atom, MM real)."""
    Z = [ZOF[a[0]] for a in ACID]
    frozen = [a[4] == -1 for a in ACID]
    steps = qm_steps([list(a[5:8]) for a in ACID], frozen, n)
    route = '#p oniom(b3lyp/6-31g(d):amber)=embedcharge opt nosymm geom=connectivity'
    out = [g_header(route, 'propanoic acid + 2 H2O, QM/MM test'), ' Symbolic Z-matrix:\n',
           ' Charge =  0 Multiplicity = 1 for low   level calculation on real  system.\n',
           ' Charge =  0 Multiplicity = 1 for high  level calculation on model system.\n',
           ' Charge =  0 Multiplicity = 1 for low   level calculation on model system.\n']
    for i, a in enumerate(ACID):
        spec = '%s-%s-%s' % (a[0], a[1], ('%.4f' % a[2]))
        link = '    H-HC    %d' % (ACID_LINK[1] + 1) if i == ACID_LINK[0] else ''
        out.append(' %-28s%3d %13.6f %13.6f %13.6f %s%s\n' % (spec, a[4], a[5], a[6], a[7], a[3], link))
    out += [' \n', ' ONIOM: Cut between C /H  %d and C  %d factor= %.6f %.6f\n' % (ACID_LINK[0] + 1, ACID_LINK[1] + 1, GFAC, GFAC),
            ' NAtoms=     17 NQM=        7 NQMF=       0 NMMI=      0 NMMIF=      0\n', '                NMic=       10 NMicF=      0.\n',
            '                    Isotopes and Nuclear Properties:\n',
            ' GradGradGradGradGradGradGradGradGradGradGradGradGradGradGradGradGradGrad\n Berny optimization.\n Initialization pass.\n']
    qm = [i for i, a in enumerate(ACID) if a[3] == 'H']
    energies, layers = [], []
    for k, X in enumerate(steps):
        e_tot = -268.400000 + 0.012 * math.exp(-0.9 * k) if k < n - 1 else -268.400000
        lr = -0.0281 + 0.002 * math.exp(-k)
        lm = 0.0135 + 0.001 * math.exp(-k)
        hm = e_tot - lr + lm
        energies.append(e_tot)
        layers.append({'high-model': hm, 'low-model': lm, 'low-real': lr})
        model_Z = [Z[i] for i in qm] + [1]
        model_X = [X[i] for i in qm] + [link_pos(X, ACID_LINK[0], ACID_LINK[1])]
        out += [g_orient('Input', Z, X), ' ONIOM: generating point charges with Q=   0.000000 Q1=   0.000000 Q2=   0.000000 Q3=   0.000000\n',
                g_mm_energy(lm), g_orient('Input', model_Z, model_X), g_scf_done('RB3LYP', hm), g_mm_energy(lr),
                ' ONIOM: calculating energy.\n',
                ' ONIOM: gridpoint  1 method:  low   system:  model energy: %20.12f\n' % lm,
                ' ONIOM: gridpoint  2 method:  high  system:  model energy: %20.12f\n' % hm,
                ' ONIOM: gridpoint  3 method:  low   system:  real  energy: %20.12f\n' % lr,
                ' ONIOM: extrapolated energy = %20.12f\n' % e_tot,
                ' ONIOM: Dipole moment (Debye):\n    X=      1.2345 Y=     -0.4321 Z=      0.0123  Tot=      1.3080\n',
                g_forces(Z, k), g_conv(k + 1, 100, k, n)]
    out.append(' Normal termination of Gaussian 16 at Sun Nov  9 12:34:56 2026.\n')
    open(path, 'w').write(''.join(out))
    return {'frames': n, 'atoms': 17, 'energies': energies, 'layerEnergies': layers,
            'atomLayer': [a[3] for a in ACID], 'frozen': frozen, 'mmTypes': [a[1] for a in ACID], 'mmCharges': [a[2] for a in ACID],
            'linkAtoms': [{'atom': ACID_LINK[0], 'bondedTo': ACID_LINK[1], 'element': 'H', 'type': 'HC'}],
            'labels': {'H': 'B3LYP/6-31G(d)', 'M': None, 'L': 'AMBER'}, 'embedding': 'electronic'}


# 3-layer: CHO (high) - CH2 (medium) - CH3 + water (low)
THREE = [('C', 'C_R', 'H', 0.000, 0.000, 0.000, None), ('O', 'O_2', 'H', -0.600, 1.050, 0.000, None),
         ('H', 'H_', 'H', -0.540, -0.960, 0.000, None), ('C', 'C_3', 'M', 1.510, 0.000, 0.000, 1),
         ('H', 'H_', 'M', 1.870, -0.520, 0.890, None), ('H', 'H_', 'M', 1.870, -0.520, -0.890, None),
         ('C', 'C_3', 'L', 2.040, 1.440, 0.000, 4), ('H', 'H_', 'L', 3.130, 1.440, 0.000, None),
         ('H', 'H_', 'L', 1.680, 1.960, 0.890, None), ('H', 'H_', 'L', 1.680, 1.960, -0.890, None),
         ('O', 'O_3', 'L', -1.900, 3.300, 0.300, None), ('H', 'H_', 'L', -1.500, 2.450, 0.200, None),
         ('H', 'H_', 'L', -2.850, 3.150, 0.350, None)]


def write_g16_oniom3(path, n=5):
    """3-layer ONIOM(B3LYP/6-31G(d):PM6:UFF) optimisation with Input + Standard orientation of the real system, the
    model and mid systems' own orientation tables, a re-printed real geometry, five sub-calculations per step."""
    Z = [ZOF[a[0]] for a in THREE]
    steps = qm_steps([list(a[3:6]) for a in THREE], [False] * len(THREE), n, 0.03)
    route = '#p oniom(b3lyp/6-31g(d):pm6:uff) opt'
    out = [g_header(route, '3-layer ONIOM test'), ' Symbolic Z-matrix:\n',
           ' Charge =  0 Multiplicity = 1 for low   level calculation on real  system.\n',
           ' Charge =  0 Multiplicity = 1 for med   level calculation on mid   system.\n',
           ' Charge =  0 Multiplicity = 1 for low   level calculation on mid   system.\n',
           ' Charge =  0 Multiplicity = 1 for high  level calculation on model system.\n',
           ' Charge =  0 Multiplicity = 1 for med   level calculation on model system.\n']
    for a in THREE:
        link = '    H-H_    %d' % a[6] if a[6] else ''
        out.append(' %-28s%14.8f%14.8f%14.8f %s%s\n' % ('%s-%s' % (a[0], a[1]), a[3], a[4], a[5], a[2], link))
    out += [' \n', ' GradGradGradGradGradGradGradGradGradGradGradGradGradGradGradGradGradGrad\n Berny optimization.\n Initialization pass.\n']
    rot = lambda X: [[-p[1], p[0], p[2]] for p in X]
    hi = [0, 1, 2]
    mid = [0, 1, 2, 3, 4, 5]
    energies, layers = [], []
    for k, X in enumerate(steps):
        e_tot = -114.600000 + 0.009 * math.exp(-0.8 * k) if k < n - 1 else -114.600000
        lr, lmid, mm_, mmod = 0.0034 + 0.001 * math.exp(-k), 0.0123 + 0.0005 * math.exp(-k), -0.0456 - 0.0004 * math.exp(-k), -0.0234 - 0.0002 * math.exp(-k)
        hm = e_tot - (mm_ - mmod + lr - lmid)
        energies.append(e_tot)
        layers.append({'high-model': hm, 'med-mid': mm_, 'med-model': mmod, 'low-real': lr, 'low-mid': lmid})
        S = rot(X)
        modZ = [Z[i] for i in hi] + [1]
        modX = [S[i] for i in hi] + [link_pos(S, 3, 0)]
        midZ = [Z[i] for i in mid] + [1]
        midX = [S[i] for i in mid] + [link_pos(S, 6, 3)]
        out += [g_orient('Input', Z, X), g_orient('Standard', Z, S),
                ' ONIOM: generating point  1 -- low level on real system.\n', g_mm_energy(lr),
                ' ONIOM: generating point  2 -- med level on mid system.\n', g_orient('Standard', midZ, midX), g_scf_done('RPM6', mm_),
                ' ONIOM: generating point  3 -- low level on mid system.\n', g_mm_energy(lmid),
                g_orient('Input', Z, X),
                ' ONIOM: generating point  4 -- high level on model system.\n', g_orient('Standard', modZ, modX), g_scf_done('RB3LYP', hm),
                ' ONIOM: generating point  5 -- med level on model system.\n', g_orient('Standard', modZ, modX), g_scf_done('RPM6', mmod),
                ' ONIOM: calculating energy.\n',
                ' ONIOM: gridpoint  1 method:  low   system:  mid   energy: %20.12f\n' % lmid,
                ' ONIOM: gridpoint  2 method:  med   system:  mid   energy: %20.12f\n' % mm_,
                ' ONIOM: gridpoint  3 method:  med   system:  model energy: %20.12f\n' % mmod,
                ' ONIOM: gridpoint  4 method:  high  system:  model energy: %20.12f\n' % hm,
                ' ONIOM: gridpoint  5 method:  low   system:  real  energy: %20.12f\n' % lr,
                ' ONIOM: extrapolated energy = %20.12f\n' % e_tot, g_forces(Z, k), g_conv(k + 1, 100, k, n)]
    out.append(' Normal termination of Gaussian 16 at Sun Nov  9 12:34:56 2026.\n')
    open(path, 'w').write(''.join(out))
    return {'frames': n, 'atoms': len(THREE), 'energies': energies, 'layerEnergies': layers,
            'atomLayer': [a[2] for a in THREE], 'mmTypes': [a[1] for a in THREE],
            'linkAtoms': [{'atom': 3, 'bondedTo': 0, 'element': 'H', 'type': 'H_'}, {'atom': 6, 'bondedTo': 3, 'element': 'H', 'type': 'H_'}],
            'labels': {'H': 'B3LYP/6-31G(d)', 'M': 'PM6', 'L': 'UFF'}}


def o_cart(Z, X, title='CARTESIAN COORDINATES (ANGSTROEM)'):
    out = ['-' * len(title) + '\n', title + '\n', '-' * len(title) + '\n']
    for z, p in zip(Z, X):
        out.append('  %-2s  %12.6f%12.6f%12.6f\n' % ({1: 'H', 6: 'C', 7: 'N', 8: 'O'}[z], p[0], p[1], p[2]))
    return ''.join(out) + '\n'


def o_conv_qm(k, n, dE):
    return o_conv(k, n, dE)


def write_orca_qmmm(path, n=5):
    """ORCA 6 QM/MM optimisation (%qmmm with an ORCAFF force field, *pdbfile input): every cycle prints the full
    system, then the QM region with its link atom, the QM SCF energy and the QM/MM energy summary."""
    Z = [ZOF[a[0]] for a in ACID]
    frozen = [a[4] == -1 for a in ACID]
    steps = qm_steps([list(a[5:8]) for a in ACID], frozen, n)
    qm = [i for i, a in enumerate(ACID) if a[3] == 'H']
    inp = ['! QMMM B3LYP def2-SVP D3BJ TightSCF Opt', '%maxcore 2000', '%qmmm', '  ORCAFFFilename "acid_water.ORCAFF.prms"',
           '  QMAtoms {4:10} end', '  ActiveAtoms {0:13} end', '  Charge_Total 0', '  Mult_Total 1', 'end', '*pdbfile 0 1 acid_water.pdb', '']
    out = [orca_header(inp, '6.0.1'),
           '------------------------------------------------------------------------------\n',
           '                                QM/MM SETUP\n',
           '------------------------------------------------------------------------------\n',
           'Total system size                               ...       17\n',
           'Number of QM atoms                              ...        7\n',
           'Number of link atoms                            ...        1\n',
           'Number of active atoms                          ...       14\n',
           'Embedding type                                  ... electrostatic\n',
           'Link atom   1: QM atom    4 (C) - MM atom    0 (C), g = %.6f\n\n' % GFAC,
           '                       *****************************\n                       * Geometry Optimization Run *\n'
           '                       *****************************\n\n']
    energies, comps = [], []
    for k, X in enumerate(steps):
        e = -268.43 + 0.011 * math.exp(-0.9 * k) if k < n - 1 else -268.43
        eqm, emm = e + 0.0447 - 0.0012 * math.exp(-k), -0.0412 + 0.0008 * math.exp(-k)
        eqmmm = e - eqm - emm
        energies.append(e)
        comps.append({'QM': eqm, 'MM': emm, 'QM-MM': eqmmm})
        subZ = [Z[i] for i in qm] + [1]
        subX = [X[i] for i in qm] + [link_pos(X, ACID_LINK[0], ACID_LINK[1])]
        last = k == n - 1
        blocks = [(o_cycle_head(k + 1), True)] + ([(None, False)] if last else [])
        for head, real_cycle in blocks:
            if head:
                out.append(head)
            else:
                out.append('\n                    ***********************HURRAY********************\n'
                           '                    ***        THE OPTIMIZATION HAS CONVERGED     ***\n'
                           '                    *************************************************\n\n'
                           '                       *******************************************************\n'
                           '                       *** FINAL ENERGY EVALUATION AT THE STATIONARY POINT ***\n'
                           '                       *******************************************************\n\n')
            out += [o_cart(Z, X), 'QM/MM: QM region with link atoms (7 QM atoms, 1 link atoms)\n', o_cart(subZ, subX),
                    o_scf(9, eqm),
                    '\n------------------------------------------------------------------------------\n',
                    '                              QM/MM ENERGY\n',
                    '------------------------------------------------------------------------------\n',
                    'E(QM)                                  ...   %18.10f Eh\n' % eqm,
                    'E(MM)                                  ...   %18.10f Eh\n' % emm,
                    'E(QM-MM)                               ...   %18.10f Eh\n' % eqmmm,
                    'E(QM/MM)                               ...   %18.10f Eh\n' % e,
                    '------------------------------------------------------------------------------\n\n',
                    '-------------------------   --------------------\nFINAL SINGLE POINT ENERGY      %18.12f\n-------------------------   --------------------\n\n' % e]
            if real_cycle:
                out.append(o_conv(k, n, e - energies[k - 1] if k else 0.0))
    out.append('\n                             ****ORCA TERMINATED NORMALLY****\n')
    open(path, 'w').write(''.join(out))
    return {'frames': n, 'atoms': 17, 'energies': energies, 'components': comps,
            'atomLayer': ['H' if a[3] == 'H' else 'L' for a in ACID], 'active': [i <= 13 for i in range(17)],
            'linkAtoms': [{'atom': ACID_LINK[0], 'bondedTo': ACID_LINK[1], 'element': 'H'}],
            'labels': {'H': 'B3LYP/def2-SVP', 'M': None, 'L': 'MM (acid_water.ORCAFF.prms)'}}


ETOH = [('C', -0.000, 0.000, 0.000), ('H', -0.370, 1.030, 0.000), ('H', -0.370, -0.510, 0.890), ('H', -0.370, -0.510, -0.890),
        ('C', 1.520, 0.000, 0.000), ('H', 1.890, -0.520, 0.890), ('H', 1.890, -0.520, -0.890), ('O', 2.030, 1.330, 0.000),
        ('H', 2.990, 1.300, 0.000), ('O', 3.500, 3.600, 0.200), ('H', 2.900, 2.850, 0.150), ('H', 4.400, 3.300, 0.200)]


def write_orca_qmxtb(path, n=4):
    """ORCA 5 QM/XTB (QM/QM2) optimisation, inline *xyz input; cycles print the full system, the QM region with its
    link atom and the subtractive QM/QM2 energy contributions."""
    Z = [ZOF[a[0]] for a in ETOH]
    steps = qm_steps([list(a[1:4]) for a in ETOH], [False] * len(ETOH), n, 0.03)
    qm = [4, 5, 6, 7, 8]
    inp = ['! QM/XTB B3LYP def2-SVP Opt', '%qmmm', '  QMAtoms {4:8} end', '  Charge_Medium 0', 'end', '* xyz 0 1']
    inp += ['  %-2s %12.6f %12.6f %12.6f' % (a[0], a[1], a[2], a[3]) for a in ETOH] + ['*', '']
    out = [orca_header(inp, '5.0.4'),
           '                       *****************************\n                       * Geometry Optimization Run *\n'
           '                       *****************************\n\n']
    energies, comps = [], []
    for k, X in enumerate(steps):
        e = -175.2 + 0.008 * math.exp(-0.9 * k) if k < n - 1 else -175.2
        lr, lm = -20.1234 + 0.002 * math.exp(-k), -9.8765 + 0.001 * math.exp(-k)
        hm = e - lr + lm
        energies.append(e)
        comps.append({'low-real': lr, 'low-model': lm, 'high-model': hm})
        subZ = [Z[i] for i in qm] + [1]
        subX = [X[i] for i in qm] + [link_pos(X, 0, 4)]
        out += [o_cycle_head(k + 1), o_cart(Z, X), o_cart(subZ, subX), o_scf(8, hm),
                '\n---------------------------------------------\n', '        QM/QM2 ENERGY CONTRIBUTIONS\n', '---------------------------------------------\n',
                'E(QM2, full system)               ...   %16.10f Eh\n' % lr,
                'E(QM2, QM region)                 ...   %16.10f Eh\n' % lm,
                'E(QM, QM region)                  ...   %16.10f Eh\n' % hm,
                'E(QM/QM2)                         ...   %16.10f Eh\n' % e,
                '---------------------------------------------\n\n',
                '-------------------------   --------------------\nFINAL SINGLE POINT ENERGY      %18.12f\n-------------------------   --------------------\n\n' % e,
                o_conv(k, n, e - energies[k - 1] if k else 0.0)]
        if k == n - 1:
            out.append('\n                    ***        THE OPTIMIZATION HAS CONVERGED     ***\n\n')
    out.append('\n                             ****ORCA TERMINATED NORMALLY****\n')
    open(path, 'w').write(''.join(out))
    return {'frames': n, 'atoms': len(ETOH), 'energies': energies, 'components': comps,
            'atomLayer': ['H' if i in qm else 'M' for i in range(len(ETOH))],
            'linkAtoms': [{'atom': 0, 'bondedTo': 4, 'element': 'H'}], 'labels': {'H': 'B3LYP/def2-SVP', 'M': 'GFN2-xTB', 'L': None}}


def write_cp2k_qmmm(path, n=4):
    """CP2K QM/MM geometry optimisation: the MM (FIST) coordinates of the whole system, the QM (QUICKSTEP)
    coordinates of the QM region + link atom shifted into the QM box, QMMM energies per step."""
    Z = [ZOF[a[0]] for a in ACID]
    X = [list(a[5:8]) for a in ACID]
    qm = [i for i, a in enumerate(ACID) if a[3] == 'H']
    shift = [5.0, 4.0, 6.0]
    out = [' DBCSR| Multiplication driver                                               XSMM\n',
           ' **** **** ******  **  PROGRAM STARTED AT               2026-11-09 12:00:00.000\n',
           ' CP2K| version string:                                          CP2K version 2024.1\n',
           ' GLOBAL| Run type                                                          GEO_OPT\n\n',
           ' QMMM| Number of QM atoms                                                        %d\n' % len(qm),
           ' QMMM| Number of link atoms                                                      1\n\n',
           ' MODULE FIST:  ATOMIC COORDINATES IN angstrom\n\n',
           '  Atom  Kind  Element       X           Y           Z          Z(eff)       Mass\n\n']
    mass = {1: 1.0079, 6: 12.0107, 8: 15.9994}
    for i, (z, p) in enumerate(zip(Z, X)):
        out.append('  %6d %5d %-2s %4d %11.6f %11.6f %11.6f %10.4f %11.4f\n' % (i + 1, 1 + [1, 6, 8].index(z), {1: 'H', 6: 'C', 8: 'O'}[z], z, p[0], p[1], p[2], 0.0, mass[z]))
    out += ['\n MODULE QUICKSTEP:  ATOMIC COORDINATES IN angstrom\n\n',
            '  Atom  Kind  Element       X           Y           Z          Z(eff)       Mass\n\n']
    subZ = [Z[i] for i in qm] + [1]
    subX = [X[i] for i in qm] + [link_pos(X, ACID_LINK[0], ACID_LINK[1])]
    for i, (z, p) in enumerate(zip(subZ, subX)):
        out.append('  %6d %5d %-2s %4d %11.6f %11.6f %11.6f %10.4f %11.4f\n' % (i + 1, 1 + [1, 6, 8].index(z), {1: 'H', 6: 'C', 8: 'O'}[z], z,
                                                                           p[0] + shift[0], p[1] + shift[1], p[2] + shift[2], {1: 1.0, 6: 4.0, 8: 6.0}[z], mass[z]))
    energies = [-45.1234 - 0.004 * (1 - math.exp(-0.9 * k)) for k in range(n + 1)]
    out.append('\n  ENERGY| Total FORCE_EVAL ( QMMM ) energy [a.u.]:              %18.12f\n' % energies[0])
    for k in range(1, n + 1):
        out += [' \n --------------------------\n OPTIMIZATION STEP: %6d\n --------------------------\n\n' % k,
                '  ENERGY| Total FORCE_EVAL ( QMMM ) energy [a.u.]:              %18.12f\n\n' % energies[k],
                ' --------  Informations at step = %5d ------------\n' % k,
                '  Optimization Method        =                 BFGS\n',
                '  Total Energy               =       %16.10f\n' % energies[k],
                '  Real energy change         =       %16.10f\n' % (energies[k] - energies[k - 1]),
                '  Max. step size             =         %.10f\n' % (0.02 * math.exp(-k)),
                '  Conv. limit for step size  =         0.0030000000\n',
                '  Convergence in step size   =                   %s\n' % ('YES' if k == n else 'NO'),
                ' ---------------------------------------------------\n']
    out += [' \n *******************************************************************************\n',
            ' ***                    GEOMETRY OPTIMIZATION COMPLETED                      ***\n',
            ' *******************************************************************************\n',
            '\n **** **** ******  **  PROGRAM ENDED AT                 2026-11-09 12:01:00.000\n']
    open(path, 'w').write(''.join(out))
    return {'steps': n, 'atoms': 17, 'energies': energies[1:], 'atomLayer': ['H' if a[3] == 'H' else 'L' for a in ACID],
            'linkAtoms': [{'atom': ACID_LINK[0], 'bondedTo': ACID_LINK[1], 'element': 'H'}]}


def write_qmmm(d):
    os.makedirs(d, exist_ok=True)
    exp = {
        'g16_oniom_amber_opt.log': write_g16_oniom_amber(os.path.join(d, 'g16_oniom_amber_opt.log')),
        'g16_oniom3_opt.log': write_g16_oniom3(os.path.join(d, 'g16_oniom3_opt.log')),
        'orca6_qmmm_opt.out': write_orca_qmmm(os.path.join(d, 'orca6_qmmm_opt.out')),
        'orca5_qmxtb_opt.out': write_orca_qmxtb(os.path.join(d, 'orca5_qmxtb_opt.out')),
        'cp2k_qmmm_geoopt.out': write_cp2k_qmmm(os.path.join(d, 'cp2k_qmmm_geoopt.out')),
    }
    # xtb ONIOM: the xtb fixture of the QC tests, run as "--oniom gfn2:gfnff 1,2" (inner region O, H1)
    src = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'tests', 'fixtures', 'qc', 'xtb_water_ohess.out')
    if os.path.exists(src):
        t = open(src).read().replace('xtb water.xyz --ohess', 'xtb water.xyz --oniom gfn2:gfnff 1,2 --ohess', 1)
        open(os.path.join(d, 'xtb_oniom_water.out'), 'w').write(t)
        exp['xtb_oniom_water.out'] = {'frames': 1, 'atoms': 3, 'atomLayer': ['H', 'H', 'L'], 'labels': {'H': 'GFN2-xTB', 'M': None, 'L': 'GFN-FF'}}
    with open(os.path.join(d, 'expected.json'), 'w') as f:
        json.dump(exp, f, indent=1)
    return {k: os.path.getsize(os.path.join(d, k)) for k in exp}


def main(argv):
    if len(argv) < 3 or argv[1] not in ('gaussian', 'orca', 'qmmm'):
        print(__doc__)
        return 2
    mb = 600.0
    if '--mb' in argv:
        mb = float(argv[argv.index('--mb') + 1])
    if argv[1] == 'gaussian':
        print(json.dumps({k: v for k, v in gaussian_big(argv[2], mb).items() if k != 'energies'}))
    elif argv[1] == 'orca':
        print(json.dumps({k: v for k, v in orca_big(argv[2], mb).items() if k != 'energies'}))
    else:
        print(json.dumps(write_qmmm(argv[2])))
    return 0


if __name__ == '__main__':
    sys.exit(main(sys.argv))
