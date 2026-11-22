#!/usr/bin/env python3
"""Generate the test fixtures of the 2.3 readers (tests/fixtures/qc/):

  Molpro       molpro_water_optfreq.out
  Dalton       dalton_water_optfreq.out, dalton_water_tddft.out
  Turbomole    turbomole_water_dscf.out, turbomole_water_escf.out, turbomole_water_aoforce.out,
               turbomole_water_ricc2.out
  GAMESS       gamess_water_opt.dat (C2v, symmetry-unique atoms), gamess_water_hess.dat (C1, normal modes)
  fchk         g16_water_opt.fchk (optimisation history), g16_water_td.fchk (ETran arrays),
               qchem_water_cis.fchk (Q-Chem excitation arrays)

Run:  python3 tools/make_fixtures_qc2.py   (needs numpy; uses the water model of tools/make_fixtures.py)

The numbers are synthetic but self-consistent (the harmonic water surface of make_fixtures.py); the layouts follow
public example outputs of each program (the test files of the cclib project). tests/fixtures/qc/expected_qc2.json
receives the values as printed, for tests/test_readers_qc2.js.
"""
import json
import math
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from make_fixtures import (BOHR, CM_TO_EH, GW_E0, GW_F, GW_OCC, GW_VIRT, MASS, fchk_iarr, fchk_int, fchk_rarr,  # noqa: E402
                           fchk_real, harmonic_energy, normal_modes, water_geom)

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
QC = os.path.join(ROOT, 'tests', 'fixtures', 'qc')
EXPECT = {}
Z = [8, 1, 1]
SYMB = {8: 'O', 1: 'H'}
XE = water_geom(0.9686, 103.65)
PATH = [water_geom(0.9600, 104.50), water_geom(0.9672, 103.95), XE]   # three optimisation steps, the last converged
EN = [harmonic_energy(X, XE, GW_F, GW_E0) for X in PATH]
MODES = normal_modes(XE, GW_F)            # bend, symmetric stretch, antisymmetric stretch (ascending)
FREQ = [m['freq'] for m in MODES]
IRS = [m['ir'] for m in MODES]
SYMS = ['A1', 'A1', 'B2']
ZPE = 0.5 * sum(FREQ) * CM_TO_EH
# TD-DFT states of the fixtures: (eV, f, symmetry, multiplicity)
STATES = [(7.1032, 0.0128, 'B1', 1), (8.9054, 0.0000, 'A2', 1), (9.6510, 0.0754, 'A1', 1), (11.2063, 0.0332, 'B2', 1)]
EV = 27.211386245988


def write(name, text):
    with open(os.path.join(QC, name), 'w', newline='\n') as fh:
        fh.write(text)


def lc(m):
    """Cartesian displacement with sum(m x^2) = 1 (Dalton, Molpro)."""
    d = m['disp'].flatten()
    mass = np.repeat([MASS[z] for z in Z], 3)
    return d / math.sqrt(float(np.sum(mass * d * d)))


def unit(m):
    d = m['disp'].flatten()
    return d / np.linalg.norm(d)


def ex_frames(paths, energies):
    return [dict(xyz=[float(v) for v in X.flatten()], energy=e) for X, e in zip(paths, energies)]


# ------------------------------------------------------------------ Molpro
def make_molpro():
    E = [round(e, 12) for e in EN]
    L = [' Primary working directories    : /tmp/molpro.fixture', ' Secondary working directories  : /tmp/molpro.fixture',
         ' Wavefunction directory         : /home/user/wfu/', ' Main file repository           : /tmp/molpro.fixture/', '',
         ' id        : fixture', '', ' Nodes     nprocs', ' node1        1', '',
         ' default implementation of scratch files=df  ', '', '',
         ' Variables initialized (930), CPU time= 0.00 sec',
         ' ***,water optimisation and frequencies', ' basis=cc-pvdz', ' geometry={angstrom;', '  O', '  H1 O 0.96', '  H2 O 0.96 H1 104.5', ' }',
         ' ks,b3lyp', ' optg', ' {freq;print,hessian}                                                     !harmonic frequencies',
         ' Commands  initialized (760), CPU time= 0.00 sec, 617 directives.', ' Default parameters read. Elapsed time= 0.03 sec',
         '', '                                         ***  PROGRAM SYSTEM MOLPRO  ***',
         '                                       Copyright, TTI GmbH Stuttgart, 2015',
         '                                    Version 2018.2 linked Nov 10 2018 21:55:55', '',
         ' THRESHOLDS:', '',
         ' ZERO    =  1.00D-12  ONEINT  =  1.00D-12  TWOINT  =  1.00D-11  PREFAC  =  1.00D-14  XXXXXX  =  1.00D-09  EORDER  =  1.00D-04',
         ' ORBITAL =  1.00D-05  CIVEC   =  1.00D-05  COEFF   =  1.00D-04  PRINTCI =  5.00D-02  PUNCHCI =  9.90D+01  OPTGRAD =  3.00D-04',
         ' OPTENERG=  1.00D-06  OPTSTEP =  3.00D-04  THRGRAD =  1.00D-03  COMPRESS=  1.00D-11  VARMIN  =  1.00D-07  VARMAX  =  1.00D-03',
         '', ' SETTING BASIS          =    CC-PVDZ', '', ' Point group  C2v ', '', '', ' ATOMIC COORDINATES', '',
         ' NR  ATOM    CHARGE       X              Y              Z', '']
    for k, (z, lab) in enumerate(zip(Z, ['O', 'H1', 'H2'])):
        x, y, zz = PATH[0][k] / BOHR
        L.append('%4d  %-3s %8.2f %14.9f %14.9f %14.9f' % (k + 1, lab, z, x, y, zz))
    L += ['', ' NUCLEAR CHARGE:                   10', ' NUMBER OF CONTRACTIONS:           24   (   11A1  +    4B1  +    7B2  +    2A2  )', '',
          ' PROGRAM * RKS-SCF (Restricted closed shell Kohn-Sham)     Authors: W. Meyer, H.-J. Werner, P.J. Knowles, 1993', '',
          ' NUMBER OF ELECTRONS:       5+    5-    SPACE SYMMETRY=1    SPIN SYMMETRY: Singlet ', '',
          ' !RKS STATE  1.1 Energy              %.12f' % E[0], ' Nuclear energy                         9.18758113',
          ' !RKS STATE  1.1 Dipole moment          0.00000000     0.00000000    -0.82560000',
          ' Dipole moment /Debye                   0.00000000     0.00000000    -2.09841000', '', '',
          ' ELECTRON ORBITALS', ' =================', '', '  Orbital  Occupation    Energy     Coefficients']
    orbs = [(1, 1, 2.0, GW_OCC[0]), (2, 1, 2.0, GW_OCC[1]), (1, 3, 2.0, GW_OCC[2]), (3, 1, 2.0, GW_OCC[3]), (1, 2, 2.0, GW_OCC[4]),
            (4, 1, 0.0, GW_VIRT[0]), (2, 3, 0.0, GW_VIRT[1])]
    for n, sy, occ, e in orbs:
        L.append('   %3d.%d     %7.5f  %10.5f     1 1s    0.70000     2 1s   -0.10000' % (n, sy, occ, e))
        L.append('                                    3 1s    0.10000')
    L += ['', '', ' HOMO      1.2    %.6f =     %.4feV' % (GW_OCC[4], GW_OCC[4] * EV), ' LUMO      4.1     %.6f =       %.4feV' % (GW_VIRT[0], GW_VIRT[0] * EV), '',
          '1PROGRAM * POP (Mulliken population analysis)', '', '', ' Density matrix read from record         2100.2  Type=RKS/CHARGE (state 1.1)', '',
          ' Population analysis by basis function type', '', ' Unique atom        s        p        d        f        g    Total    Charge',
          '   1  O       3.78100  4.65500  0.01300  0.00000  0.00000  8.44900  - 0.44900',
          '   2  H       0.77300  0.00250  0.00000  0.00000  0.00000  0.77550  + 0.22450',
          '   3  H       0.77300  0.00250  0.00000  0.00000  0.00000  0.77550  + 0.22450', '', '',
          ' Entering Rational Function Geometry Optimization', '',
          ' ITER.   ENERGY(OLD)    ENERGY(NEW)      DE          GRADMAX     GRADNORM    GRADRMS     STEPMAX     STEPLEN     STEPRMS   CPU-time']
    opt = [(E[0], E[1], 1.2e-2, 1.6e-2, 5.3e-3, 1.4e-2, 1.9e-2, 6.3e-3), (E[1], E[2], 2.4e-4, 3.1e-4, 1.0e-4, 2.6e-4, 3.4e-4, 1.1e-4)]
    for k, (eo, en, gm, gn, gr, sm, sl, sr) in enumerate(opt):
        L.append('%4d  %13.8f  %13.8f  %13.8f  %10.8f  %10.8f  %10.8f  %10.8f  %10.8f  %10.8f     %5.2f' % (k + 1, eo, en, en - eo, gm, gn, gr, sm, sl, sr, 1.2 * (k + 1)))
        if k == 0:
            L.append(' Freezing grid')
    L += ['', ' END OF GEOMETRY OPTIMIZATION.', '', ' Current geometry (xyz format, in Angstrom)', '', '    3', ' KS-SCF000/CC-PVDZ  ENERGY=%.8f' % E[2]]
    for k, z in enumerate(Z):
        L.append(' %-2s %19.10f %19.10f %19.10f' % (SYMB[z], *PATH[2][k]))
    L += ['', ' PROGRAM * RKS-SCF (Restricted closed shell Kohn-Sham)     Authors: W. Meyer, H.-J. Werner, P.J. Knowles, 1993', '',
          ' NUMBER OF ELECTRONS:       5+    5-    SPACE SYMMETRY=1    SPIN SYMMETRY: Singlet ', '',
          ' !RKS STATE  1.1 Energy              %.12f' % E[2], ' !RKS STATE  1.1 Dipole moment          0.00000000     0.00000000    -0.82560000', '',
          ' Normal Modes', '', '                              1 A1        2 A1        3 B2',
          ' Wavenumbers [cm-1]   ' + ''.join('%12.2f' % f for f in FREQ),
          ' Intensities [km/mol] ' + ''.join('%12.2f' % i for i in IRS),
          ' Intensities [relative]' + ''.join('%11.2f' % (100 * i / max(IRS)) for i in IRS)]
    disp = [lc(m) for m in MODES]
    for a in range(3):
        for ax, c in enumerate('XYZ'):
            L.append('           %-12s' % (SYMB[Z[a]] + c + str(a + 1)) + ''.join('%12.5f' % d[3 * a + ax] for d in disp))
    L += ['', ' Normal Modes of low/zero frequencies', '', '                                  1           2           3           4           5           6',
          ' Wavenumbers [cm-1]            0.00        0.00        0.00        0.00        0.00        0.00',
          ' Intensities [km/mol]          0.00        0.00        0.00        0.00        0.00        0.00']
    for a in range(3):
        for c in 'XYZ':
            L.append('           %-12s' % (SYMB[Z[a]] + c + str(a + 1)) + '     0.10000' * 6)
    L += ['', '', ' Zero point energy:  %.8f [H]    %.2f [1/CM]' % (ZPE, ZPE / CM_TO_EH), '',
          ' B3LYP/cc-pVDZ energy=    %.12f' % E[2], '', ' Molpro calculation terminated']
    write('molpro_water_optfreq.out', '\n'.join(L) + '\n')
    EXPECT['molpro'] = dict(first=[float(v) for v in PATH[0].flatten()], last=[float(v) for v in PATH[2].flatten()], E=E,
                            freqs=[round(f, 2) for f in FREQ], ir=[round(i, 2) for i in IRS], zpe=round(ZPE, 8),
                            cycles=[E[0], E[1], E[2]], mulliken=[-0.449, 0.2245, 0.2245], occ=GW_OCC + GW_VIRT[:2])


# ------------------------------------------------------------------ Dalton
def dalton_head(dal, molrows):
    L = ['', '     *************** Dalton - An Electronic Structure Program ***************', '',
         '    This is output from DALTON 2015.0', '', '', '   Content of the .dal input file', ' ----------------------------------', '']
    L += molrows + dal
    L += ['', '', '       *******************************************************************',
          '       *********** Output from DALTON general input processing ***********',
          '       *******************************************************************', '']
    return L


def dalton_geom(title, X):
    L = ['', '                             %s' % title, '                             ' + '-' * len(title), '']
    for k, z in enumerate(Z):
        L.append(' %-2s %s %19.10f %19.10f %19.10f' % (SYMB[z], '   ' if z == 8 else '_%d ' % k, *(X[k] / BOHR)))
    return L + ['', '']


def make_dalton():
    E = [round(e, 12) for e in EN]
    mol = ['BASIS', 'cc-pVDZ', 'water', 'synthetic fixture', 'AtomTypes=2 Angstrom NoSymmetry', 'Charge=8.0 Atoms=1',
           'O   %15.10f %15.10f %15.10f' % tuple(PATH[0][0]), 'Charge=1.0 Atoms=2',
           'H   %15.10f %15.10f %15.10f' % tuple(PATH[0][1]), 'H   %15.10f %15.10f %15.10f' % tuple(PATH[0][2])]
    mol = [m[:50] for m in mol]
    L = dalton_head(['**DALTON', '.OPTIMIZE', '**WAVE FUNCTIONS', '.DFT', 'B3LYP', '**PROPERTIES', '.VIBANA', '**END OF'], mol)
    L += ['  Basis set used is "cc-pVDZ" from the basis set library.', '', '  Cartesian Coordinates (a.u.)', '  ----------------------------', '',
          '  Total number of coordinates:    9']
    for k, z in enumerate(Z):
        x, y, zz = PATH[0][k] / BOHR
        lab = 'O' if z == 8 else 'H%d' % k
        L.append('  %-6s:  %3d  x %14.10f  %3d  y %14.10f  %3d  z %14.10f' % (lab, 3 * k + 1, x, 3 * k + 2, y, 3 * k + 3, zz))
    L += ['', '@    Total charge of the molecule               0', '@    Spin multiplicity and 2 M_S                1         0', '',
          '@    Final DFT energy:          %18.12f' % E[0], '', ' Number of electrons :   10', '', ' Sym       Kohn-Sham orbital energies', '',
          '1 A     ' + ''.join('%15.8f' % e for e in (GW_OCC + GW_VIRT[:3])[:5]),
          '        ' + ''.join('%15.8f' % e for e in (GW_OCC + GW_VIRT[:3])[5:]), '', '',
          '    E(LUMO) :     %.8f au (symmetry 1)' % GW_VIRT[0], '', '                                  Dipole moment',
          '                                  -------------', '', '                 au               Debye          C m (/(10**-30)',
          '              0.825600           2.098512           7.000000', '', '',
          '                             Dipole moment components', '                             ------------------------', '',
          '                 au               Debye          C m (/(10**-30)', '',
          '      z     -0.82560000        -2.09851200        -7.00000000', '', '']
    L += ['   ***************************************************************************',
          '   ************************ FINAL RESULTS from ABACUS ************************',
          '   ***************************************************************************']
    L += dalton_geom('Molecular geometry (au)', PATH[0])
    L += ['           ::::::::::::::: Optimization Control Center :::::::::::::::']
    L += dalton_geom('Next geometry (au)', PATH[1])
    L += ['@    Final DFT energy:          %18.12f' % E[1], '']
    L += dalton_geom('Next geometry (au)', PATH[2])
    L += ['@    Final DFT energy:          %18.12f' % E[2], '']
    L += dalton_geom('Final geometry (au)', PATH[2])
    L += dalton_geom('Molecular geometry (au)', PATH[2])
    L += ['@ Iter     Energy        Change       GradNorm  Index   StepLen    TrustRad #Rej',
          '@ ------------------------------------------------------------------------------']
    rows = [(0, E[0], 0.0, 0.024512, 0.031820), (1, E[1], E[1] - E[0], 0.000458, 0.000610), (2, E[2], E[2] - E[1], 0.000012, 0.000015)]
    for it, e, de, g, st in rows:
        L.append('@ %3d    %12.6f %12.6f %12.6f    0   %8.6f*   0.500000   0 ' % (it, e, de, g, st))
    L += ['@  *) Newton step taken.', '', '@ Geometry converged in  3 iterations!', '',
          '  Vibrational Frequencies and IR Intensities', '  ------------------------------------------', '',
          '  mode   irrep        frequency             IR intensity', '  ============================================================',
          '                  cm-1       hartrees     km/mol   (D/A)**2/amu', '  ------------------------------------------------------------']
    order = [2, 1, 0]   # Dalton lists the highest first
    for k, q in enumerate(order):
        L.append('  %2d      %-4s %9.2f    %8.6f   %8.3f   %6.4f' % (k + 1, SYMS[q], FREQ[q], FREQ[q] * CM_TO_EH, IRS[q], IRS[q] / 42.2561))
    L += ['', '', '', '  Normal Coordinates (bohrs*amu**(1/2)):', '  --------------------------------------', '', '',
          '              ' + ''.join('%3d %5d   ' % (k + 1, round(FREQ[q])) for k, q in enumerate(order)),
          '  ----------------------------------------------------------------------', '']
    disp = [lc(MODES[q]) for q in order]
    for a in range(3):
        for ax, c in enumerate('xyz'):
            L.append('    %-6s %s ' % (SYMB[Z[a]], c) + ''.join('%12.6f' % d[3 * a + ax] for d in disp))
        L.append('')
    L += ['', '                              Total Molecular Energy', '                              ----------------------', '',
          '             electronic     vibrational           total    energy ', '',
          '            %11.6f       %9.6f     %11.6f    Hartrees' % (E[2], ZPE, E[2] + ZPE), '', '',
          ' >>>> Total CPU  time used in DALTON:   2.40 seconds']
    write('dalton_water_optfreq.out', '\n'.join(L) + '\n')
    EXPECT['dalton'] = dict(frames=ex_frames(PATH, E), freqs=[round(f, 2) for f in FREQ], ir=[round(i, 3) for i in IRS], zpe=round(ZPE, 6),
                            dipole=round(0.8256 * 2.541746, 6), occ=GW_OCC + GW_VIRT[:3])
    # TD-DFT: coordinates only in the echo of the .mol file (cut at 50 columns), states in three symmetries
    molt = [m[:50] for m in ['BASIS', 'cc-pVDZ', 'water', 'synthetic fixture', 'AtomTypes=2 Angstrom NoSymmetry', 'Charge=8.0 Atoms=1',
                             'O   %15.10f %15.10f %15.10f' % tuple(XE[0]), 'Charge=1.0 Atoms=2',
                             'H   %15.10f %15.10f %15.10f' % tuple(XE[1]), 'H   %15.10f %15.10f %15.10f' % tuple(XE[2])]]
    L = dalton_head(['**DALTON', '.RUN RESPONSE', '**WAVE FUNCTIONS', '.DFT', 'B3LYP', '**RESPONSE', '*LINEAR', '.SINGLE RESIDUE', '.ROOTP', ' 2', '**END OF'], molt)
    L += ['  Basis set used is "cc-pVDZ" from the basis set library.', '', '@    Final DFT energy:          %18.12f' % EN[2], '']
    bysym = {}
    for ev, f, sy, mult in STATES:
        bysym.setdefault(sy, []).append((ev, f))
    for s_i, sy in enumerate(['A1', 'B1', 'B2', 'A2']):
        if sy not in bysym:
            continue
        L += ['', ' *** @ Excit. operator sym %d & ref. state sym 1 => excited state symmetry %d  ( %s ) ***' % (s_i + 1, s_i + 1, sy), '']
        for k, (ev, f) in enumerate(bysym[sy]):
            L += [' @ Excited state no:    %d in symmetry %d  ( %s )' % (k + 1, s_i + 1, sy), ' ' + '-' * 50, '',
                  '@ Excitation energy :  %.8f     au' % (ev / EV), '@                      %.7f     eV;   %9.3f     cm-1' % (ev, ev * 8065.544), '',
                  '@ Total energy :      %.5f     au' % (EN[2] + ev / EV), '']
            if f > 0:   # the components of a dipole-allowed state: two lines for B states in this fixture
                parts = [('ZDIPLEN', f)] if sy == 'A1' else [('XDIPLEN', f * 0.25), ('YDIPLEN', f * 0.75)]
                for op, fp in parts:
                    L += ['@ Operator type:    %s ' % op, '@ Oscillator strength (LENGTH)   :  %.8E  (Transition moment :  0.10000000     )' % fp, '']
    L += ['', ' >>>> Total CPU  time used in DALTON:   1.10 seconds']
    write('dalton_water_tddft.out', '\n'.join(L) + '\n')
    EXPECT['dalton_td'] = dict(xyz=[float(v) for v in XE.flatten()],
                               states=sorted([[ev, round(f * 0.25, 8) + round(f * 0.75, 8) if sy != 'A1' else f] for ev, f, sy, m in STATES]))


# ------------------------------------------------------------------ Turbomole
def tm_head(prog, title):
    return ['  ========================================================================', '',
            ' %s (fixture-host) : TURBOMOLE rev. V7.4.1 (2bfdd732) compiled 4 Nov 2019 at 15:12:03' % prog,
            ' Copyright (C) 2019 TURBOMOLE GmbH, Karlsruhe', '', '', '                             %s' % title, '']


def tm_coords(X):
    L = ['              +--------------------------------------------------+', '              | Atomic coordinate, charge and isotop information |',
         '              +--------------------------------------------------+', '', '                    atomic coordinates            atom    charge  isotop']
    for k, z in enumerate(Z):
        L.append('    %14.8f %14.8f %14.8f    %-2s   %6.3f     0' % (*(X[k] / BOHR), SYMB[z].lower(), z))
    return L + [' ', '       center of nuclear mass  :    0.00000000    0.00000000    0.00000000', '',
                '   type   atoms  prim   cont   basis', '   ---------------------------------------------------------------------------',
                '    o        1     27     14   def2-SVP   [3s2p1d|7s4p1d]', '    h        2      7      5   def2-SVP   [2s1p|4s1p]',
                '   ---------------------------------------------------------------------------', '']


def make_turbomole():
    E = EN[2]
    L = tm_head('dscf', 'd s c f - program') + tm_coords(XE)
    L += ['           ------------------', '           density functional', '           ------------------', ' B3-LYP hybrid functional',
          ' exchange:    LDA + Becke + 0.2 HF', '', '    mo occupation :', '   irrep   mo\'s   occupied', '    a       24        5', '',
          ' convergence criteria satisfied after    10 iterations', '', '                  ------------------------------------------ ',
          '                 |  total energy      =    %.11f  |' % E, '                  ------------------------------------------ ', '',
          ' ==============================================================================', '                           electrostatic moments',
          ' ==============================================================================', '', '              nuc           elec       ->  total',
          ' ------------------------------------------------------------------------------', '                          charge      ',
          ' ------------------------------------------------------------------------------', '          10.000000     -10.000000      -0.000000', ' ',
          ' ------------------------------------------------------------------------------', '                       dipole moment  ',
          ' ------------------------------------------------------------------------------', '   x       0.000000       0.000000       0.000000',
          '   y       0.000000       0.000000       0.000000', '   z       0.120000      -0.945600      -0.825600', ' ',
          '   | dipole moment | =     0.8256 a.u. =     2.0985 debye ', '', '    dscf ended normally']
    write('turbomole_water_dscf.out', '\n'.join(L) + '\n')
    # escf: one block per state, with the rotatory strength (ECD) and the dominant contributions
    L = tm_head('escf', 'e s c f') + tm_coords(XE)
    for k, (ev, f, sy, mult) in enumerate(sorted(STATES)):
        L += ['', '', '                         %d singlet %s excitation' % (k + 1, sy.lower()), '', '',
              ' Total energy:                           %.16f    ' % (E + ev / EV), '',
              ' Excitation energy:                      %.16f    ' % (ev / EV), '', ' Excitation energy / eV:                  %.15f    ' % ev, '',
              ' Excitation energy / nm:                  %.13f    ' % (1239.84198 / ev), '', '', ' Oscillator strength:', '',
              '    velocity representation:             %.16f    ' % (f * 1.05), '', '    length representation:               %.16E' % f, '',
              '    mixed representation:                %.16f    ' % (f * 1.02), '', '', ' Rotatory strength:', '',
              '    velocity representation:             0.1282566570726007E-10', '',
              '    velocity rep. / 10^(-40)erg*cm^3:    %.16E' % (0.5 * (k + 1)), '', '    length representation:               0.1242930160103103E-10', '',
              '    length rep. / 10^(-40)erg*cm^3:      %.16E' % (0.4 * (k + 1) * (-1) ** k), '', '', ' Dominant contributions:', '',
              '      occ. orbital   energy / eV   virt. orbital     energy / eV   |coeff.|^2*100',
              '        5 a1            -7.95           6 a1               1.75       %4.1f' % (96.5 - k), ' ', '', ' Change of electron number:', '', '             0.000000', '']
    L += ['', '    escf ended normally']
    write('turbomole_water_escf.out', '\n'.join(L) + '\n')
    EXPECT['turbomole_escf'] = dict(states=[[ev, f] for ev, f, sy, m in sorted(STATES)], R=[0.4 * (k + 1) * (-1) ** k for k in range(4)],
                                     weight=[(96.5 - k) / 100 for k in range(4)])
    # aoforce: six translations and rotations (0.00), then the three vibrations, six columns per block
    L = tm_head('force', 'a o f o r c e - program') + tm_coords(XE)
    L += ['          ---------------------------------------------------', '          NORMAL MODES and VIBRATIONAL FREQUENCIES (cm**(-1))',
          '          ---------------------------------------------------', '']
    fr = [0.0] * 6 + FREQ
    irs = [0.0] * 6 + IRS
    sy = [''] * 6 + [s.lower() for s in SYMS]
    disp = [np.full(9, 0.1)] * 6 + [unit(m) for m in MODES]
    for b in range(0, 9, 6):
        cols = range(b, min(b + 6, 9))
        L += ['', '       mode     ' + ''.join('%9d' % (c + 1) for c in cols), '', '     frequency  ' + ''.join('%9.2f' % fr[c] for c in cols), '',
              '     symmetry   ' + ''.join('%9s' % sy[c] for c in cols), '', '        IR      ' + ''.join('%9s' % ('YES' if fr[c] else '-') for c in cols),
              '|dDIP/dQ|   (a.u.)' + ''.join('%9.4f' % 0.01 for c in cols), 'intensity (km/mol)' + ''.join('%9.2f' % irs[c] for c in cols),
              'intensity (  %   )' + ''.join('%9.2f' % 1.0 for c in cols), '', '       RAMAN    ' + ''.join('%9s' % '-' for c in cols), '']
        for a in range(3):
            for ax, c in enumerate('xyz'):
                head = '%3d   %-2s        %s ' % (a + 1, SYMB[Z[a]].lower(), c) if ax == 0 else '                  %s ' % c
                L.append(head + ''.join('%9.5f' % disp[q][3 * a + ax] for q in cols))
        L += ['', 'reduced mass(g/mol)' + ''.join('%9.3f' % 1.0 for c in cols)]
    L += ['', '         *****************************************************', '         *  zero point VIBRATIONAL energy  :      %.7f  Hartree  *' % ZPE,
          '         *****************************************************', '', '    force ended normally']
    write('turbomole_water_aoforce.out', '\n'.join(L) + '\n')
    EXPECT['turbomole_force'] = dict(freqs=[round(f, 2) for f in FREQ], ir=[round(i, 2) for i in IRS], zpe=round(ZPE, 7))
    # ricc2: MP2 and CC2 ground state, CC2 states (a table per model; the last is final) and their strengths
    L = tm_head('ricc2', 'R I C C 2 - PROGRAM') + tm_coords(XE)
    L += ['              *   Final MP2 energy      :   %.10f       *' % (E - 0.2101), '     *   Final CC2 energy                        :   %.10f      *' % (E - 0.2153), '',
          ' +================================================================================+',
          ' | sym | multi | state |          CCS excitation energies           |',
          ' +================================================================================+']
    for k, (ev, f, sy, mult) in enumerate(sorted(STATES)):
        L.append(' | a   |   1   |  %2d   |  %11.7f |  %11.5f |  %10.3f |' % (k + 1, (ev + 0.3) / EV, ev + 0.3, (ev + 0.3) * 8065.5))
    L += [' +================================================================================+', '',
          ' | sym | multi | state |          CC2 excitation energies       |  %t1   |  %t2   |',
          ' |     |       |       +----------------------------------------+--------+--------+',
          ' |     |       |       |   Hartree    |    eV      |    cm-1    |    %   |    %   |',
          ' +================================================================================+']
    for k, (ev, f, sy, mult) in enumerate(sorted(STATES)):
        L.append(' | a   |   1   |  %2d   |  %11.7f | %10.5f | %10.3f |  94.89 |   5.11 |' % (k + 1, ev / EV, ev, ev * 8065.5))
    L += [' +================================================================================+', '',
          '     | sym | multi | state |    E / eV  | % singles | % doubles |']
    for k, (ev, f, sy, mult) in enumerate(sorted(STATES)):
        L.append('     | a   |   1   |  %2d   |  %9.5f |   95.51   |    4.49   |' % (k + 1, ev + 9))
    L += ['', ' **************************************************************', ' *<<<<<<<<<<<<  ONE-PHOTON ABSORPTION STRENGTHS  >>>>>>>>>>>>>*',
          ' **************************************************************', '']
    for k, (ev, f, sy, mult) in enumerate(sorted(STATES)):
        L += ['  +=========================================================================+', '  |  Transition                                             model: CC2      |',
              '  |    number, symmetry, multiplicity:  %3d a    1                          |' % (k + 1),
              '  |    frequency :   %.10f a.u.     %.5f e.V.     %.1f rcm     |' % (ev / EV, ev, ev * 8065.5),
              '  +=========================================================================+', '', '       oscillator strength (length gauge)   :      %.8f' % f, '']
    L += ['', '    ricc2 : all done']
    write('turbomole_water_ricc2.out', '\n'.join(L) + '\n')
    EXPECT['turbomole_ricc2'] = dict(cc2=round(E - 0.2153, 10), mp2=round(E - 0.2101, 10), states=[[ev, f] for ev, f, sy, m in sorted(STATES)])
    EXPECT['turbomole'] = dict(xyz=[float(v) for v in XE.flatten()], energy=round(E, 11), dipole=round(0.8256 * 2.541746, 5))


# ------------------------------------------------------------------ GAMESS punch
def punch_atoms(X, idx):
    return [' %-2s          %3.1f %17.10f %17.10f %17.10f' % (SYMB[Z[k]], Z[k], *X[k]) for k in idx]


def make_gamess():
    E = [round(e, 10) for e in EN]
    # C2v in GAMESS's frame (C2 = z, a mirror plane xz): water lies in the yz plane, the unique atoms are O and one H
    L = [' $DATA  ', 'water C2v optimisation'.ljust(80), 'CNV      2', '']
    for k in (0, 1):
        L += ['%-2s          %3.1f %17.10f %17.10f %17.10f' % (SYMB[Z[k]], Z[k], *PATH[0][k]), '   N31     6', '        ']
    L += [' $END      ']
    for step in range(3):
        L += ['-------------------- DATA FROM NSERCH=%4d --------------------' % step, ' COORDINATES OF SYMMETRY UNIQUE ATOMS (ANGS)',
              '   ATOM   CHARGE       X              Y              Z', ' ------------------------------------------------------------']
        L += punch_atoms(PATH[step], (0, 1))
        L += ['--- CLOSED SHELL ORBITALS --- GENERATED AT Mon Aug  5 12:54:36 2019', 'water C2v optimisation'.ljust(80),
              'E(R-B3LYP)=%20.10f, E(NUC)=    9.1875811300,   12 ITERS' % E[step], ' $VEC   ', ' 1  1 1.00000000E+00', ' $END   ']
    L += ['----- RESULTS FROM SUCCESSFUL RHF      GEOMETRY SEARCH -----', '----- COORDS, ORBS, GRADIENT, AND APPROX. HESSIAN -----',
          ' COORDINATES OF SYMMETRY UNIQUE ATOMS (ANGS)', '   ATOM   CHARGE       X              Y              Z',
          ' ------------------------------------------------------------']
    L += punch_atoms(PATH[2], (0, 1))
    L += ['--- OPTIMIZED RHF      MO-S --- GENERATED AT Mon Aug  5 12:55:16 2019', 'E=%20.10f, E(NUC)=    9.1875811300' % E[2], ' $VEC   ', ' $END   ']
    write('gamess_water_opt.dat', '\n'.join(L) + '\n')
    EXPECT['gamess_opt'] = dict(frames=ex_frames(PATH, E))
    # C1 Hessian job: every atom, the MOLPLT normal modes (translations and rotations included) and the AIMPAC section
    L = [' $DATA  ', 'water hessian'.ljust(80), 'C1       0']
    for k in range(3):
        L += ['%-2s          %3.1f %17.10f %17.10f %17.10f' % (SYMB[Z[k]], Z[k], *XE[k]), '   N31     6', '        ']
    L += [' $END      ', '--- CLOSED SHELL ORBITALS --- GENERATED AT Mon Aug  5 13:05:02 2019', 'water hessian'.ljust(80),
          'E(R-B3LYP)=%20.10f, E(NUC)=    9.1875811300,   12 ITERS' % E[2], ' $VEC   ', ' $END   ',
          ' MOMENTS AT POINT    1 X,Y,Z=  0.000000  0.000000  0.000000', ' DIPOLE       0.000000  0.000000 -2.098512',
          '----- START OF NORMAL MODES FOR -MOLPLT- PROGRAM -----', 'ATOMIC MASSES', '    15.99491     1.00783     1.00783']
    low = [3.21, 1.07, 0.02, 0.51, 24.6, 35.2]
    allm = [(low[k], np.roll(np.full(9, 0.05), k)) for k in range(6)] + [(f, unit(m)) for f, m in zip(FREQ, MODES)]
    for k, (f, v) in enumerate(allm):
        L.append('MODE %5d   FREQUENCY=%10.5f (CM**-1)' % (k + 1, f))
        for a in range(3):
            L.append(' %16.9E %16.9E %16.9E' % tuple(v[3 * a:3 * a + 3]))
    L += ['----- END OF NORMAL MODES FOR -MOLPLT- PROGRAM -----', "----- TOP OF INPUT FILE FOR BADER'S AIMPAC PROGRAM -----", 'water hessian',
          'GAUSSIAN              5 MOL ORBITALS     24 PRIMITIVES        3 NUCLEI']
    for k in range(3):
        L.append('  %-2s %4d    (CENTRE %2d)  %11.8f %11.8f %11.8f  CHARGE = %4.1f' % (SYMB[Z[k]], k + 1, k + 1, *(XE[k] / BOHR), Z[k]))
    L += ['CENTRE ASSIGNMENTS    1  1  1', 'EXPONENTS']
    for k, e in enumerate(GW_OCC):
        L += ['MO %2d                     OCC NO =   2.00000000 ORB. ENERGY = %12.8f' % (k + 1, e), '']
    L += ['END DATA', 'RHF      ENERGY =   %18.10f   VIRIAL(-V/T)  =   2.00000000' % E[2], "----- END OF INPUT FILE FOR BADER'S AIMPAC PROGRAM -----"]
    write('gamess_water_hess.dat', '\n'.join(L) + '\n')
    EXPECT['gamess_hess'] = dict(xyz=[float(v) for v in XE.flatten()], freqs=[round(f, 5) for f in FREQ], energy=E[2], occ=GW_OCC)


# ------------------------------------------------------------------ fchk
def make_fchk_files():
    # an optimisation: "Opt point 1" with every geometry and energy; the last geometry is the current one
    head = ['Water optimisation'.ljust(72), '%-10s%-60s%-20s' % ('FOpt', 'RB3LYP', '6-31G(d)')]
    L = head + fchk_int('Number of atoms', 3) + fchk_int('Charge', 0) + fchk_int('Multiplicity', 1)
    L += ['Route                                      C   N=           2', '#p b3lyp/6-31g(d) opt                                       ']
    L += fchk_int('Number of alpha electrons', 5) + fchk_int('Number of beta electrons', 5) + fchk_iarr('Atomic numbers', Z)
    L += fchk_rarr('Current cartesian coordinates', list((PATH[2] / BOHR).flatten()))
    L += fchk_real('SCF Energy', EN[2]) + fchk_real('Total Energy', EN[2])
    L += fchk_int('Optimization Num results per geometry', 2) + fchk_int('Optimization Num geometry variables', 9)
    L += fchk_rarr('Opt point       1 Results for each geome', [v for e in EN for v in (e, 0.0)])
    L += fchk_rarr('Opt point       1 Geometries', [v for X in PATH for v in (X / BOHR).flatten()])
    L += fchk_iarr('Optimization Number of geometries', [3])
    write('g16_water_opt.fchk', '\n'.join(L) + '\n')
    EXPECT['fchk_opt'] = dict(E=[float('%.8E' % e) for e in EN[:2]] + [EN[2]], last=[float('%.8E' % v) * BOHR for v in (PATH[2] / BOHR).flatten()])
    # a TD job: ETran arrays (16 values per state: energy, transition dipoles...), singlets (0) and triplets (2)
    states = [(5.9031, 0.0, 2), (7.1032, 0.0128, 0), (7.6640, 0.0, 2), (9.6510, 0.0754, 0)]
    vals = []
    for ev, f, sp in states:
        de = ev / EV
        mu = math.sqrt(1.5 * f / de) if f > 0 else 0.0
        vals += [EN[2] + de, 0.0, 0.0, mu] + [0.0] * 12
    L = ['Water TD-DFT'.ljust(72), '%-10s%-60s%-20s' % ('SP', 'RB3LYP TD-FC', '6-31G(d)')]
    L += fchk_int('Number of atoms', 3) + fchk_int('Charge', 0) + fchk_int('Multiplicity', 1) + fchk_iarr('Atomic numbers', Z)
    L += fchk_rarr('Current cartesian coordinates', list((XE / BOHR).flatten()))
    L += fchk_real('SCF Energy', EN[2]) + fchk_real('CIS Energy', EN[2] + states[1][0] / EV) + fchk_real('Total Energy', EN[2] + states[1][0] / EV)
    L += fchk_iarr('ETran scalars', [4, 16, 1, 10, 2] + [0] * 15) + fchk_iarr('ETran spin', [s[2] for s in states]) + fchk_iarr('ETran sym', [0] * 4)
    L += fchk_rarr('ETran state values', vals)
    write('g16_water_td.fchk', '\n'.join(L) + '\n')
    EXPECT['fchk_td'] = dict(states=[[ev, f, 'Singlet' if sp == 0 else 'Triplet'] for ev, f, sp in states])
    # Q-Chem: excitation energies (Eh) and oscillator strengths, triplets first as Q-Chem writes them
    qs = [(0.2169, 0.0), (0.2816, 0.0), (0.2610, 0.0128), (0.3547, 0.0754)]
    L = ['Jobname.Temp', '%-10s%-30s%-30s' % ('SP', 'R', '6-31G*')]
    L += fchk_int('Number of atoms', 3) + fchk_int('Charge', 0) + fchk_int('Multiplicity', 1) + fchk_iarr('Atomic numbers', Z)
    L += fchk_rarr('Current cartesian coordinates', list((XE / BOHR).flatten())) + fchk_int('Number of Excited States', 4)
    L += fchk_rarr('Excitation Energies', [q[0] for q in qs]) + fchk_rarr('Oscillator Strengths', [q[1] for q in qs])
    write('qchem_water_cis.fchk', '\n'.join(L) + '\n')
    EXPECT['fchk_qchem'] = dict(states=[[float('%.8E' % q[0]) * EV, q[1]] for q in qs])


def main():
    make_molpro()
    make_dalton()
    make_turbomole()
    make_gamess()
    make_fchk_files()
    with open(os.path.join(QC, 'expected_qc2.json'), 'w') as fh:
        json.dump(EXPECT, fh, indent=1)
    print('written to', QC)


if __name__ == '__main__':
    main()
