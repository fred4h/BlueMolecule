#!/usr/bin/env python3
"""Generate the BlueMolecule parser test fixtures (tests/fixtures/*).

Run:  python3 tools/make_fixtures.py

All numbers are synthetic but self-consistent:
  * water vibrations come from a valence force field (stretch, stretch-stretch,
    bend, stretch-bend) turned into a Cartesian Hessian H = B^T F B and
    diagonalised with masses, so frequencies, reduced masses, force constants
    and displacement vectors agree with each other;
  * IR intensities come from fixed point charges (dmu/dQ);
  * optimisation steps use the same harmonic surface for energies and forces;
  * cube data are analytic functions sampled on the grid.
Layouts follow the public text formats of Gaussian 09/16, ORCA 4/5/6, Molden,
Gaussian fchk, MDL V2000, PDB and GROMACS .gro files.
tests/fixtures/expected.json receives reference values (rounded exactly as
printed in the fixtures) for tests/test_parsers.js.
"""
import json
import math
import os

import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
FIX = os.path.join(ROOT, 'tests', 'fixtures')
BOHR = 0.529177210903
MDYNA_TO_EH = 1e-18 / 4.3597447222071e-18      # 1 mdyn*Angstrom in Hartree
MDYN_TO_EH_BOHR = 1e-8 / 8.2387234983e-8        # 1 mdyn in Hartree/Bohr
CM_TO_EH = 1.0 / 219474.6313705
NU_CONST = math.sqrt(100.0 / 1.66053906660e-27) / (2 * math.pi * 2.99792458e10)
IR_CONST = 974.8801                              # (e/amu^0.5)^2 -> km/mol
EXPECT = {}

SYM = {0: 'X', 1: 'H', 6: 'C', 7: 'N', 8: 'O', 11: 'Na', 17: 'Cl', 30: 'Zn', 78: 'Pt'}
MASS = {1: 1.00783, 6: 12.0, 7: 14.00307, 8: 15.99491, 17: 34.96885, 78: 194.96479}


def write(name, text):
    with open(os.path.join(FIX, name), 'w', newline='\n') as fh:
        fh.write(text)


def r4(x, n):
    return float(round(float(x), n)) + 0.0


# ---------------------------------------------------------------- water model
def water_geom(r, theta_deg):
    """C2v water in the Gaussian standard orientation (YZ plane, C2 = z,
    centre of nuclear charge at the origin)."""
    half = math.radians(theta_deg) / 2
    y, dz = r * math.sin(half), r * math.cos(half)
    zo = 2 * dz / 10.0
    return np.array([[0.0, 0.0, zo], [0.0, y, zo - dz], [0.0, -y, zo - dz]])


def internals(X):
    v1, v2 = X[1] - X[0], X[2] - X[0]
    r1, r2 = np.linalg.norm(v1), np.linalg.norm(v2)
    return np.array([r1, r2, math.acos(np.dot(v1, v2) / (r1 * r2))])


def bmatrix(X, h=1e-6):
    B = np.zeros((3, 9))
    x = X.flatten()
    for k in range(9):
        xp, xm = x.copy(), x.copy()
        xp[k] += h
        xm[k] -= h
        B[:, k] = (internals(xp.reshape(3, 3)) - internals(xm.reshape(3, 3))) / (2 * h)
    return B


def force_field(frr, frr2, ftt, frt):
    return np.array([[frr, frr2, frt], [frr2, frr, frt], [frt, frt, ftt]])


def normal_modes(X, F, Z=(8, 1, 1), charges=(-0.66, 0.33, 0.33)):
    B = bmatrix(X)
    H = B.T @ F @ B                       # mdyn/Angstrom
    m = np.repeat([MASS[z] for z in Z], 3)
    w, V = np.linalg.eigh(H / np.sqrt(np.outer(m, m)))
    modes = []
    for i in np.argsort(w)[-3:]:
        lam, L = w[i], V[:, i]
        lc = L / np.sqrt(m)               # Angstrom per amu^0.5
        mu = 1.0 / np.sum(lc ** 2)
        d = lc * math.sqrt(mu)            # normalised Cartesian displacement
        big = np.argmax(np.abs(d))
        if d[big] < 0:
            d, lc = -d, -lc
        dmu = np.zeros(3)
        for a in range(3):
            dmu += charges[a] * lc[3 * a:3 * a + 3]
        modes.append(dict(freq=NU_CONST * math.sqrt(lam), mu=mu, k=lam * mu, dmu=dmu,
                          ir=IR_CONST * float(dmu @ dmu), disp=d.reshape(3, 3)))
    return modes


def harmonic_energy(X, Xe, F, e0):
    dq = internals(X) - internals(Xe)
    return e0 + 0.5 * dq @ F @ dq * MDYNA_TO_EH


def cart_gradient(X, Xe, F):
    dq = internals(X) - internals(Xe)
    return (bmatrix(X).T @ (F @ dq)) * MDYN_TO_EH_BOHR    # Hartree/Bohr


def thermo_from(freqs, e_el, T=298.15):
    zpe = 0.5 * sum(freqs) * CM_TO_EH
    rt = 3.166811563e-6 * T
    ecorr = zpe + 3 * rt + 0.000003
    hcorr = ecorr + rt
    gcorr = hcorr - T * 45.10 / 627509.474          # T*S with S = 45.10 cal/(mol K)
    return dict(zpe=zpe, eCorr=ecorr, hCorr=hcorr, gCorr=gcorr, eZPE=e_el + zpe,
                eTotal=e_el + ecorr, H=e_el + hcorr, G=e_el + gcorr, temperature=T)


def rotation(ax, ay, az):
    cx, sx, cy, sy, cz, sz = (math.cos(ax), math.sin(ax), math.cos(ay), math.sin(ay),
                              math.cos(az), math.sin(az))
    Rx = np.array([[1, 0, 0], [0, cx, -sx], [0, sx, cx]])
    Ry = np.array([[cy, 0, sy], [0, 1, 0], [-sy, 0, cy]])
    Rz = np.array([[cz, -sz, 0], [sz, cz, 0], [0, 0, 1]])
    return Rz @ Ry @ Rx


# ------------------------------------------------------------ Gaussian pieces
def fd(x, n=6):
    """Fortran D-exponent string, e.g. -1.234567D-08."""
    return ('%.*E' % (n, x)).replace('E', 'D')


def clean(v, eps=5e-6):
    return 0.0 if abs(v) < eps else v


def g_banner(ver, rev):
    year = {'16': '2019', '09': '2013'}[ver]
    build = {'16': 'ES64L-G16Rev%s  3-Jul-2019', '09': 'EM64L-G09Rev%s 24-Apr-2013'}[ver] % rev
    pub = {'16': '2016', '09': '2013'}[ver]
    return [' Entering Gaussian System, Link 0=g%s' % ver,
            ' Input=job.gjf', ' Output=job.log', ' Initial command:',
            ' /opt/g%s/l1.exe "/scratch/Gau-12345.inp" -scrdir="/scratch/"' % ver,
            ' Entering Link 1 = /opt/g%s/l1.exe PID=     12345.' % ver, '  ',
            ' Copyright (c) 1988-%s, Gaussian, Inc.  All Rights Reserved.' % year, '  ',
            ' This is part of the Gaussian(R) %s program.  It is based on' % ver,
            ' the Gaussian(R) 09 system (copyright 2009, Gaussian, Inc.),', ' ... (legal notice abridged) ...',
            '  ', ' Cite this work as:', ' Gaussian %s, Revision %s,' % (ver, rev),
            ' M. J. Frisch, G. W. Trucks, H. B. Schlegel, G. E. Scuseria, et al.',
            ' Gaussian, Inc., Wallingford CT, %s.' % pub, ' ',
            ' ******************************************',
            ' Gaussian %s:  %s' % (ver, build),
            '                 4-Nov-2026 ',
            ' ******************************************']


def g_box(lines):
    d = ' ' + '-' * max(len(s) for s in lines)
    return [d] + [' ' + s for s in lines] + [d]


def g_job_head(link0, route_lines, title, charge, mult, Z, X, iops=None):
    L = []
    for l0 in link0:
        L.append(' ' + l0)
        if l0.lower().startswith('%nproc'):
            L.append(' Will use up to    4 processors via shared memory.')
    L += g_box(route_lines)
    L += iops or [' 1/18=20,19=15,26=3,38=1/1,3;', ' 2/9=110,12=2,17=6,18=5,40=1/2;',
                  ' 3/5=1,6=6,7=1,11=2,25=1,30=1,71=1,74=-5/1,2,3;', ' 4//1;',
                  ' 5/5=2,38=5/2;', ' 6/7=2,8=2,9=2,10=2,28=1/1;', ' 99/5=1,9=1/99;']
    L += g_box([title])
    L += [' Symbolic Z-matrix:', ' Charge = %2d Multiplicity = %d' % (charge, mult)]
    for z, x in zip(Z, X):
        L.append(' %-2s %20.5f%10.5f%10.5f ' % (SYM[z], x[0], x[1], x[2]))
    L.append(' ')
    return L


def g_orient(kind, Z, X):
    head = {'Input': '                          Input orientation:                          ',
            'Standard': '                         Standard orientation:                         ',
            'Z-Matrix': '                         Z-Matrix orientation:                         '}[kind]
    dash = ' ' + '-' * 69
    L = [head, dash, ' Center     Atomic      Atomic             Coordinates (Angstroms)',
         ' Number     Number       Type             X           Y           Z', dash]
    for i, (z, x) in enumerate(zip(Z, X)):
        L.append('%7d%11d%12d    %12.6f%12.6f%12.6f' % (i + 1, z, 0, clean(x[0], 5e-7),
                                                      clean(x[1], 5e-7), clean(x[2], 5e-7)))
    L.append(dash)
    return L


def g_scf(method, e, cycles=10):
    return [' SCF Done:  E(%s) = %16.10f     A.U. after %4d cycles' % (method, e, cycles),
            '            NFock= %2d  Conv=0.12D-08     -V/T= 2.0101' % cycles]


def eig_lines(tag, vals):
    return [' %s eigenvalues --%s' % (tag, ''.join('%10.5f' % v for v in vals[i:i + 5]))
            for i in range(0, len(vals), 5)]


def g_atom_table(header, Z, cols, fmt='%11.6f', sum_label=None):
    L = [' ' + header, '               1' + ''.join('%11d' % (c + 2) for c in range(len(cols) - 1))]
    for a, z in enumerate(Z):
        L.append('%6d  %-2s' % (a + 1, SYM[z]) + ''.join(fmt % c[a] for c in cols))
    if sum_label:
        L.append(' Sum of %s = %9.5f' % (sum_label, sum(cols[0])))
    return L


def g_pop(Z, occ, virt, mull, dip, state='1-A1', dip_old=False, mull_header='Mulliken charges:',
          hirsh=None, extra=None):
    L = [' **********************************************************************', '',
         '            Population analysis using the SCF Density.', '',
         ' **********************************************************************', '',
         ' Orbital symmetries:', '       Occupied  (A1) (A1) (B2) (A1) (B1)',
         '       Virtual   (A1) (B2) (B2) (A1) (B1) (A1) (B2) (A1) (A2) (A1)',
         ' The electronic state is %s.' % state]
    L += eig_lines('Alpha  occ.', occ) + eig_lines('Alpha virt.', virt)
    L += ['          Condensed to atoms (all electrons):']
    L += g_atom_table(mull_header, Z, [mull], sum_label='Mulliken charges' if 'atomic' not in mull_header
                      else 'Mulliken atomic charges')
    L += [' Mulliken charges with hydrogens summed into heavy atoms:', '               1']
    L.append('%6d  %-2s%11.6f' % (1, SYM[Z[0]], sum(mull)))
    L += [' Electronic spatial extent (au):  <R**2>=             19.0123',
          ' Charge=              0.0000 electrons',
          ' Dipole moment (field-independent basis, Debye):']
    w = 11 if dip_old else 20
    L.append('    X=%*.4f    Y=%*.4f    Z=%*.4f  Tot=%*.4f' % (w, dip[0], w, dip[1], w, dip[2], w,
                                                           math.sqrt(sum(d * d for d in dip))))
    L += [' Quadrupole moment (field-independent basis, Debye-Ang):',
          '   XX=             -7.1234   YY=             -4.0123   ZZ=             -6.0123']
    if extra:
        L += extra
    if hirsh:
        L += [' Hirshfeld charges, spin densities, dipoles, and CM5 charges using IRadAn=      4:',
              '              Q-H        S-H        Dx         Dy         Dz        Q-CM5']
        for a, z in enumerate(Z):
            L.append('%6d  %-2s' % (a + 1, SYM[z]) + ''.join('%11.6f' % v for v in hirsh[a]))
        L.append('       Tot ' + ''.join('%11.6f' % sum(h[c] for h in hirsh) for c in range(6)))
    return L


def g_conv(vals, pred):
    thr = [0.000450, 0.000300, 0.001800, 0.001200]
    names = ['Maximum Force', 'RMS     Force', 'Maximum Displacement', 'RMS     Displacement']
    L = ['         Item               Value     Threshold  Converged?']
    for n, v, t in zip(names, vals, thr):
        L.append(' %-20s %12.6f %12.6f     %s' % (n, v, t, 'YES' if v <= t else 'NO '))
    L.append(' Predicted change in Energy=%s' % fd(pred))
    return L


def g_berny(step, maxstep=20, scan=None):
    s = ' Step number %3d out of a maximum of %4d' % (step, maxstep)
    if scan:
        s += ' on scan point %5d out of %5d' % scan
    return [' GradGradGradGradGradGradGradGradGradGradGradGradGradGradGradGradGradGradGradGrad',
            ' Berny optimization.', ' Search for a local minimum.', s,
            ' All quantities printed in internal units (Hartrees-Bohrs-Radians).']


def g_completed(params):
    L = [' Optimization completed.', '    -- Stationary point found.',
         '                           ----------------------------',
         '                           !   Optimized Parameters   !',
         '                           ! (Angstroms and Degrees)  !',
         ' --------------------------                            --------------------------',
         ' ! Name  Definition              Value          Derivative Info.                !',
         ' --------------------------------------------------------------------------------']
    for name, defin, val in params:
        L.append(' ! %-5s %-16s %13.4f         -DE/DX =    0.0                 !' % (name, defin, val))
    L += [' --------------------------------------------------------------------------------',
          ' GradGradGradGradGradGradGradGradGradGradGradGradGradGradGradGradGradGradGradGrad']
    return L


FREQ_HEAD = [' Harmonic frequencies (cm**-1), IR intensities (KM/Mole), Raman scattering',
             ' activities (A**4/AMU), depolarization ratios for plane and unpolarized',
             ' incident light, reduced masses (AMU), force constants (mDyne/A),',
             ' and normal coordinates:']


def g_freq_regular(modes, syms, Z, raman=None, depol=None):
    L = list(FREQ_HEAD)
    for b in range(0, len(modes), 3):
        blk = list(range(b, min(b + 3, len(modes))))

        def row(label, vals):
            return ' %-14s' % label + (' ' * 11).join('%12.4f' % v for v in vals)
        L.append(''.join('%23d' % (i + 1) for i in blk))
        L.append(''.join('%23s' % syms[i] for i in blk))
        L.append(row('Frequencies --', [modes[i]['freq'] for i in blk]))
        L.append(row('Red. masses --', [modes[i]['mu'] for i in blk]))
        L.append(row('Frc consts  --', [modes[i]['k'] for i in blk]))
        L.append(row('IR Inten    --', [modes[i]['ir'] for i in blk]))
        if raman:
            L.append(row('Raman Activ --', [raman[i] for i in blk]))
            L.append(row('Depolar (P) --', [depol[i] for i in blk]))
            L.append(row('Depolar (U) --', [2 * depol[i] / (1 + depol[i]) for i in blk]))
        L.append('  Atom  AN      X      Y      Z' + '        X      Y      Z' * (len(blk) - 1))
        for a, z in enumerate(Z):
            s = '%6d%4d' % (a + 1, z)
            for i in blk:
                s += '  ' + ''.join('%7.2f' % clean(v, 0.005) for v in modes[i]['disp'][a])
            L.append(s)
    return L


def g_freq_hp(modes, syms, Z, raman=None, depol=None):
    L = list(FREQ_HEAD)
    for b in range(0, len(modes), 5):
        blk = list(range(b, min(b + 5, len(modes))))

        def row(label, vals):
            return '%22s ' % label + ''.join('%10.4f' % v for v in vals)
        L.append(' ' * 23 + ''.join('%10d' % (i + 1) for i in blk))
        L.append(' ' * 23 + ''.join('%10s' % syms[i] for i in blk))
        L.append(row('Frequencies ---', [modes[i]['freq'] for i in blk]))
        L.append(row('Reduced masses ---', [modes[i]['mu'] for i in blk]))
        L.append(row('Force constants ---', [modes[i]['k'] for i in blk]))
        L.append(row('IR Intensities ---', [modes[i]['ir'] for i in blk]))
        if raman:
            L.append(row('Raman Activities ---', [raman[i] for i in blk]))
            L.append(row('Depolar (P) ---', [depol[i] for i in blk]))
            L.append(row('Depolar (U) ---', [2 * depol[i] / (1 + depol[i]) for i in blk]))
        L.append(' Coord Atom Element:')
        for a, z in enumerate(Z):
            for c in range(3):
                L.append('%4d%6d%6d       ' % (c + 1, a + 1, z) +
                         ''.join('%10.5f' % clean(modes[i]['disp'][a][c]) for i in blk))
    return L


def g_thermo(Z, th):
    L = [' -------------------', ' - Thermochemistry -', ' -------------------',
         ' Temperature %9.3f Kelvin.  Pressure   1.00000 Atm.' % th['temperature']]
    for a, z in enumerate(Z):
        L.append(' Atom %5d has atomic number %2d and mass %9.5f' % (a + 1, z, MASS[z]))
    L += [' Molecular mass:    18.01056 amu.',
          ' Zero-point vibrational energy      55424.3 (Joules/Mol)',
          '                                   13.24671 (Kcal/Mol)',
          ' Zero-point correction=                   %16.6f (Hartree/Particle)' % th['zpe'],
          ' Thermal correction to Energy=            %16.6f' % th['eCorr'],
          ' Thermal correction to Enthalpy=          %16.6f' % th['hCorr'],
          ' Thermal correction to Gibbs Free Energy= %16.6f' % th['gCorr'],
          ' Sum of electronic and zero-point Energies=        %16.6f' % th['eZPE'],
          ' Sum of electronic and thermal Energies=           %16.6f' % th['eTotal'],
          ' Sum of electronic and thermal Enthalpies=         %16.6f' % th['H'],
          ' Sum of electronic and thermal Free Energies=      %16.6f' % th['G']]
    return L


def g_end(ver='16', ok=True, link1=False):
    L = [' 1\\1\\GINC-NODE01\\SP\\RB3LYP\\6-31G(d)\\H2O1\\USER\\04-Nov-2026\\0\\\\#p\\\\title\\\\0,1\\O,0.,0.,0.1\\@',
         '', ' THE ONLY WAY TO HAVE A FRIEND IS TO BE ONE.  -- EMERSON',
         ' Job cpu time:       0 days  0 hours  0 minutes 12.3 seconds.',
         ' Elapsed time:       0 days  0 hours  0 minutes  3.4 seconds.',
         ' File lengths (MBytes):  RWF=      6 Int=      0 D2E=      0 Chk=      1 Scr=      1']
    if ok:
        L.append(' Normal termination of Gaussian %s at Wed Nov  4 21:00:00 2026.' % ver)
    else:
        L += [' Error termination request processed by link 9999.',
              ' Error termination via Lnk1e in /opt/g%s/l9999.exe at Wed Nov  4 21:00:00 2026.' % ver]
    if link1:
        L.append(' Link1:  Proceeding to internal job step number  2.')
    return L


# ------------------------------------------------------- Gaussian fixtures
GW_F = force_field(8.45, -0.10, 0.70, 0.22)
GW_E0 = -76.4089533224
GW_Z = [8, 1, 1]
GW_RAMAN = [4.8971, 105.3392, 42.1957]
GW_DEPOL = [0.6934, 0.1521, 0.7500]
GW_OCC = [-19.13798, -0.99766, -0.51936, -0.37203, -0.29215]
GW_VIRT = [0.06424, 0.14857, 0.79948, 0.84744, 0.88786, 0.89907, 1.03163, 1.12895, 1.71745,
           1.74055, 1.80276, 2.40306, 2.41802, 3.31468]
GW_MULL = [-0.868866, 0.434433, 0.434433]
GW_HIRSH = [[-0.330523, 0.0, 0.0, 0.0, -0.188839, -0.648178],
            [0.165262, 0.0, 0.0, 0.105386, 0.038936, 0.324089],
            [0.165261, 0.0, 0.0, -0.105386, 0.038936, 0.324089]]
GW_DIP = [0.0, 0.0, -2.0985]
ROT_IN = rotation(0.3, -0.2, 0.5)


def to_input(X):
    return X @ ROT_IN.T + np.array([0.1, -0.2, 0.3])


def mode_syms(modes):
    return ['B2' if abs(m['disp'][0][1]) > 1e-4 else 'A1' for m in modes]


def mode_expect(modes, nd):
    return dict(freqs=[r4(m['freq'], 4) for m in modes], redMass=[r4(m['mu'], 4) for m in modes],
                forceConst=[r4(m['k'], 4) for m in modes], ir=[r4(m['ir'], 4) for m in modes],
                modes=[[r4(clean(v, 0.5 * 10 ** -nd), nd) for v in m['disp'].flatten()] for m in modes])


def low_freq_lines(modes):
    return [' Full mass-weighted force constant matrix:',
            ' Low frequencies ---   -0.0012    0.0009    0.0012   10.5236   15.0987   19.0000',
            ' Low frequencies --- ' + ' '.join('%9.4f' % m['freq'] for m in modes)]


def make_gaussian_water():
    Xe = water_geom(0.9686, 103.65)
    steps = [(0.9500, 105.50), (0.9712, 103.10), (0.9688, 103.62), (0.96862, 103.648)]
    geoms = [water_geom(r, t) for r, t in steps]
    energies = [harmonic_energy(X, Xe, GW_F, GW_E0) for X in geoms]
    route = '#p opt freq=raman b3lyp/6-31g(d) pop=hirshfeld scrf=(pcm,solvent=water) int=ultrafine'
    title = 'Water optimization and frequencies'
    L = g_banner('16', 'C.01')
    L += g_job_head(['%nprocshared=4', '%mem=4GB', '%chk=water.chk'], [route[:70], route[70:]], title,
                    0, 1, GW_Z, to_input(geoms[0]))
    L += [' GradGradGradGradGradGradGradGradGradGradGradGradGradGradGradGradGradGradGradGrad',
          ' Berny optimization.', ' Initialization pass.']
    opts = []
    for k, X in enumerate(geoms):
        L += g_orient('Input', GW_Z, to_input(X))
        L += [' Stoichiometry    H2O', ' Framework group  C2V[C2(O),SGV(H2)]', ' Deg. of freedom     2',
              ' Full point group                 C2V     NOp   4']
        L += g_orient('Standard', GW_Z, X)
        L += [' Rotational constants (GHZ):    823.8149880    437.8516440    285.8733330']
        if k == 0:
            L += [' Standard basis: 6-31G(d) (6D, 7F)',
                  '    19 basis functions,    36 primitive gaussians,    19 cartesian basis functions',
                  '     5 alpha electrons        5 beta electrons',
                  '       nuclear repulsion energy         9.1681207385 Hartrees.',
                  ' ------------------------------------------------------------------------------',
                  ' Polarizable Continuum Model (PCM)',
                  ' Solvent              : Water, Eps=  78.355300 Eps(inf)=   1.777849',
                  ' ------------------------------------------------------------------------------']
        L += g_scf('RB3LYP', energies[k], 12 if k == 0 else 8)
        if k == 0:
            L += g_pop(GW_Z, [v - 0.002 for v in GW_OCC], GW_VIRT, [-0.872345, 0.436172, 0.436173],
                       [0.0, 0.0, -2.1234])
        g = cart_gradient(X, Xe, GW_F)
        if k + 1 < len(geoms):
            dx = (geoms[k + 1] - X).flatten() / BOHR
        else:
            dx = np.array([0.000052, -0.000035, 0.000021])
        vals = [float(np.max(np.abs(g))), float(np.sqrt(np.mean(g ** 2))),
                float(np.max(np.abs(dx))), float(np.sqrt(np.mean(dx ** 2)))]
        if k == len(geoms) - 1:
            vals = [max(v, 1e-6) for v in vals]
        opts.append([r4(v, 6) for v in vals])
        L += g_berny(k + 1)
        L += g_conv(vals, -(energies[k] - energies[min(k + 1, len(geoms) - 1)]) - 1e-9)
        if k < len(geoms) - 1:
            L.append(' GradGradGradGradGradGradGradGradGradGradGradGradGradGradGradGradGradGradGradGrad')
    final = geoms[-1]
    q = internals(final)
    L += g_completed([('R1', 'R(1,2)', q[0]), ('R2', 'R(1,3)', q[1]), ('A1', 'A(2,1,3)', math.degrees(q[2]))])
    L += [' Largest change from initial coordinates is atom    2       0.031 Angstoms.']
    L += g_orient('Input', GW_Z, to_input(final)) + g_orient('Standard', GW_Z, final)
    L += g_pop(GW_Z, GW_OCC, GW_VIRT, GW_MULL, GW_DIP, hirsh=GW_HIRSH)
    L += g_end(link1=True)
    # ---- job 2: frequencies at the optimised geometry
    L += g_box(['#P Geom=AllCheck Guess=TCheck SCRF=Check GenChk RB3LYP/6-31G(d) Freq=Raman'])
    L += [' 1/10=4,29=7,30=1,38=1,40=1/1,3;', ' 2/12=2,40=1/2;', ' 99//99;',
          ' Structure from the checkpoint file:  "water.chk"']
    L += g_box([title])
    L += [' Charge =  0 Multiplicity = 1', ' Redundant internal coordinates found in file.  (old form).',
          ' O,0,0.,0.,%.10f' % final[0][2], ' H,0,0.,%.10f,%.10f' % (final[1][1], final[1][2]),
          ' H,0,0.,%.10f,%.10f' % (final[2][1], final[2][2]), ' Recover connectivity data from disk.']
    L += g_orient('Input', GW_Z, to_input(final)) + g_orient('Standard', GW_Z, final)
    L += g_scf('RB3LYP', energies[-1], 1)
    L += g_pop(GW_Z, GW_OCC, GW_VIRT, GW_MULL, GW_DIP, hirsh=GW_HIRSH)
    modes = normal_modes(final, GW_F)
    L += low_freq_lines(modes)
    L += g_freq_regular(modes, mode_syms(modes), GW_Z, GW_RAMAN, GW_DEPOL)
    th = thermo_from([m['freq'] for m in modes], energies[-1])
    L += g_thermo(GW_Z, th)
    L += g_end()
    write('g16_water_optfreq.log', '\n'.join(L) + '\n')
    e = mode_expect(modes, 2)
    e.update(energies=[r4(x, 10) for x in energies], opt=opts, sym=mode_syms(modes),
             final=[r4(clean(v, 5e-7), 6) for v in final.flatten()],
             thermo={k: r4(v, 6) for k, v in th.items()}, route=route, title=title)
    EXPECT['g16_water_optfreq'] = e
    return modes, final, energies[-1]


def make_gaussian_hpmodes():
    X = water_geom(0.9686, 103.65)
    modes = normal_modes(X, GW_F)
    L = g_banner('16', 'A.03')
    L += g_job_head(['%chk=water_hp.chk'], ['#p freq=(raman,hpmodes) b3lyp/6-31g(d)'], 'Water hpmodes',
                    0, 1, GW_Z, X)
    L += g_orient('Input', GW_Z, X) + g_orient('Standard', GW_Z, X)
    L += [' Standard basis: 6-31G(d) (6D, 7F)']
    L += g_scf('RB3LYP', GW_E0, 11)
    L += g_pop(GW_Z, GW_OCC, GW_VIRT, GW_MULL, GW_DIP)
    L += low_freq_lines(modes)
    L += g_freq_hp(modes, mode_syms(modes), GW_Z, GW_RAMAN, GW_DEPOL)
    L += [' ']
    L += g_freq_regular(modes, mode_syms(modes), GW_Z, GW_RAMAN, GW_DEPOL)
    L += g_thermo(GW_Z, thermo_from([m['freq'] for m in modes], GW_E0))
    L += g_end()
    write('g16_water_hpmodes.log', '\n'.join(L) + '\n')
    EXPECT['g16_water_hpmodes'] = mode_expect(modes, 5)


TD_STATES = [  # index, label, eV, f, S2, [(from, arrow, to, coeff)]
    (1, 'Triplet-B1', 6.5432, 0.0000, 2.000, [(5, '->', 6, 0.70123)]),
    (2, 'Singlet-B1', 7.1178, 0.0156, 0.000, [(5, '->', 6, 0.70488), (5, '<-', 6, -0.10234)]),
    (3, 'Triplet-A1', 8.1234, 0.0000, 2.000, [(4, '->', 6, 0.12345), (5, '->', 7, 0.69012)]),
    (4, 'Singlet-A2', 8.9012, 0.0000, 0.000, [(5, '->', 7, 0.70655)])]


def make_gaussian_td():
    X = water_geom(0.9686, 103.65)
    e0 = -76.4229876543
    L = g_banner('16', 'C.01')
    L += g_job_head(['%chk=water_td.chk'], ['#p td=(nstates=2,50-50) b3lyp/6-31+g(d) pop=(nbo,mk)'],
                    'Water TD-DFT', 0, 1, GW_Z, X)
    L += g_orient('Input', GW_Z, X) + g_orient('Standard', GW_Z, X)
    L += [' Standard basis: 6-31+G(d) (6D, 7F)']
    L += g_scf('RB3LYP', e0, 13)
    L += [' **********************************************************************', '',
          '            Excited states from <AA,BB:AA,BB> singles matrix:', '',
          ' **********************************************************************', '',
          ' Excitation energies and oscillator strengths:', '']
    for idx, lab, ev, f, s2, tr in TD_STATES:
        L.append(' Excited State %3d:      %-10s%8.4f eV %7.2f nm  f=%6.4f  <S**2>=%5.3f'
                 % (idx, lab, ev, 1239.84193 / ev, f, s2))
        for a, arrow, b, c in tr:
            L.append('%8d %s %3d %15.5f' % (a, arrow, b, c))
        if idx == 1:
            L += [' This state for optimization and/or second-order correction.',
                  ' Total Energy, E(TD-HF/TD-DFT) =  %.10f' % (e0 + ev / 27.211386245988),
                  ' Copying the excited state density for this state as the 1-particle RhoCI density.']
        L.append(' ')
    L += [' SavETr:  write IOETrn=   770 NScale= 10 NData=  16 NLR=1 NState=    4 LETran=      82.']
    L += g_pop(GW_Z, GW_OCC, GW_VIRT, [-0.912345, 0.456172, 0.456173], [0.0, 0.0, -2.2345])
    esp = [-0.789123, 0.394561, 0.394562]
    L += [' Fitting point charges to electrostatic potential',
          ' Charges from ESP fit, RMS=   0.00210 RRMS=   0.17227:',
          ' Charge=   0.00000 Dipole=     0.0000     0.0000    -0.8456 Tot=     0.8456']
    L += g_atom_table('ESP charges:', GW_Z, [esp], sum_label='ESP charges')
    npa = [-0.91935, 0.45968, 0.45967]
    L += ['', ' Summary of Natural Population Analysis:', '',
          '                                       Natural Population',
          '                Natural  -----------------------------------------------',
          '    Atom  No    Charge         Core      Valence    Rydberg      Total',
          ' -----------------------------------------------------------------------']
    for a, z in enumerate(GW_Z):
        L.append('%7s%5d%11.5f%13.5f%12.5f%11.5f%12.5f' % (SYM[z], a + 1, npa[a], 1.99982 if z > 1 else 0,
                                                       6.90869 if z > 1 else 0.53818, 0.01084, z - npa[a]))
    L += [' =======================================================================',
          '   * Total *    0.00000      1.99982     7.98505    0.01513    10.00000']
    L += g_end()
    write('g16_water_tddft.log', '\n'.join(L) + '\n')
    EXPECT['g16_water_tddft'] = dict(
        states=[dict(index=s[0], label=s[1], eV=s[2], nm=r4(1239.84193 / s[2], 2), f=s[3], s2=s[4],
                     transitions=[[str(t[0]), str(t[2]), t[3]] for t in s[5]]) for s in TD_STATES],
        esp=esp, npa=npa, energy=e0)


def cisplatin(pt_cl=2.33, pt_n=2.05):
    """cis-[PtCl2(NH3)2], square planar in the xy plane."""
    def u(deg):
        return np.array([math.cos(math.radians(deg)), math.sin(math.radians(deg)), 0.0])
    Z, X = [78, 17, 17], [np.zeros(3), pt_cl * u(135), pt_cl * u(225)]
    for ang in (45, -45):
        n = pt_n * u(ang)
        a = n / np.linalg.norm(n)
        p1 = np.array([0.0, 0.0, 1.0])
        p2 = np.cross(a, p1)
        Z.append(7)
        X.append(n)
        for k in range(3):
            phi = math.radians(120 * k + 30)
            t = math.radians(180 - 109.5)
            X.append(n + 1.02 * (math.cos(t) * a + math.sin(t) * (math.cos(phi) * p1 + math.sin(phi) * p2)))
            Z.append(1)
    return Z, np.array(X)


CP_MULL = [0.351234, -0.452101, -0.452099, -0.851234, 0.375905, 0.375905, 0.375901, -0.851230,
           0.375906, 0.375907, 0.375906]


def make_gaussian_g09_cisplatin():
    Z, X1 = cisplatin()
    _, X2 = cisplatin(2.345, 2.062)
    energies = [-1153.4123456789, -1153.4145678901]
    occ = ([-101.56789, -101.56788, -14.41234, -14.41230, -9.48765, -9.48761] +
           [-7.24321, -7.24320, -7.24318, -7.24317, -7.24310, -7.24309] +
           [-3.41234, -2.05123, -2.05001, -2.04876] +
           [round(-1.02 + 0.04 * i, 5) for i in range(20)])
    virt = [round(-0.07 + 0.07 * i, 5) for i in range(10)]
    esp = [0.612345, -0.523456, -0.523457, -0.912345, 0.391234, 0.391233, 0.391235, -0.912344,
           0.391234, 0.391236, 0.391235]
    L = g_banner('09', 'D.01')
    L += g_job_head(['%nprocshared=4', '%chk=cisplatin.chk'], ['#p opt=(maxcycles=2) b3lyp/genecp nosymm pop=mk'],
                    'cisplatin opt (nosymm)', 0, 1, Z, X1)
    for k, X in enumerate((X1, X2)):
        L += g_orient('Input', Z, X)
        L += [' Stoichiometry    Cl2H6N2Pt', ' Framework group  C1[X(Cl2H6N2Pt)]', ' Deg. of freedom    27',
              ' Full point group                 C1      NOp   1',
              ' Rotational constants (GHZ):      1.4567890      1.2345678      0.6789012']
        L += g_scf('RB3LYP', energies[k], 15)
        if k == 0:
            L += g_pop(Z, occ, virt, CP_MULL, [-3.1234, 3.1234, 0.0], state='1-A', dip_old=True,
                       mull_header='Mulliken atomic charges:')
            L += [' Charges from ESP fit, RMS=   0.00321 RRMS=   0.21234:',
                  ' Charge=   0.00000 Dipole=    -8.1234     5.4321     0.0000 Tot=     9.7724',
                  '              1']
            for a, z in enumerate(Z):
                L.append('%6d  %-2s%11.6f' % (a + 1, SYM[z], esp[a]))
            L.append(' -----------------------------------------------------------------')
        L += g_berny(k + 1, 2)
        L += g_conv([0.023456, 0.008765, 0.087654, 0.023456] if k == 0 else
                    [0.004567, 0.001234, 0.012345, 0.004567], -0.002345 if k == 0 else -0.000123)
    L += [' Optimization stopped.', '    -- Number of steps exceeded,  NStep=   2',
          '    -- Flag reset to prevent archiving.']
    L += g_end('09', ok=False)
    write('g09_cisplatin_input_orientation.log', '\n'.join(L) + '\n')
    EXPECT['g09_cisplatin'] = dict(Z=Z, x1=[r4(clean(v, 5e-7), 6) for v in X1.flatten()],
                                   x2=[r4(clean(v, 5e-7), 6) for v in X2.flatten()], energies=energies,
                                   occ=occ, virt=virt, mull=CP_MULL, esp=esp)
    return Z, X1


def make_gaussian_scan():
    Xe = water_geom(0.9686, 103.65)
    pts = [[(0.9700, 100.0), (0.9690, 100.0)], [(0.9690, 105.0), (0.9684, 105.0)]]
    L = g_banner('16', 'C.01')
    L += g_job_head(['%chk=scan.chk'], ['#p opt=modredundant b3lyp/6-31g(d)'], 'Water bend scan', 0, 1,
                    GW_Z, water_geom(*pts[0][0]))
    energies, labels = [], []
    for p, steps in enumerate(pts):
        for s, (r, t) in enumerate(steps):
            X = water_geom(r, t)
            e = harmonic_energy(X, Xe, GW_F, GW_E0)
            energies.append(e)
            labels.append('Scan point %d, step %d' % (p + 1, s + 1))
            L += g_orient('Input', GW_Z, X) + g_orient('Standard', GW_Z, X)
            L += g_scf('RB3LYP', e, 9)
            L += g_berny(s + 1, 20, (p + 1, 2))
            conv = s == len(steps) - 1
            L += g_conv([0.000012, 0.000008, 0.000045, 0.000030] if conv else
                        [0.004512, 0.003008, 0.012345, 0.008030], -1e-7 if conv else -2.3e-5)
            if conv:
                L += g_completed([('R1', 'R(1,2)', r), ('R2', 'R(1,3)', r), ('A1', 'A(2,1,3)', t)])
    L += [' Summary of Optimized Potential Surface Scan (add -76.0 to energies):',
          '                           1         2',
          '   Eigenvalues --  %9.5f %9.5f' % (energies[1] + 76, energies[3] + 76),
          '           R1         0.96900   0.96840']
    L += g_end()
    write('g16_water_scan.log', '\n'.join(L) + '\n')
    EXPECT['g16_scan'] = dict(energies=[r4(e, 10) for e in energies], labels=labels)


# ----------------------------------------------------------------- ORCA pieces
def o_banner(ver):
    return ['', '                                 *****************',
            '                                 * O   R   C   A *',
            '                                 *****************', '',
            '           --- An Ab Initio, DFT and Semiempirical electronic structure package ---', '',
            '                  #######################################################',
            '                  #                        -***-                        #',
            '                  #          Department of theory and spectroscopy      #',
            '                  #    Directorship and core code : Frank Neese         #',
            '                  #        Max Planck Institute fuer Kohlenforschung    #',
            '                  #                Kaiser Wilhelm Platz 1               #',
            '                  #                 D-45470 Muelheim/Ruhr               #',
            '                  #                      Germany                        #',
            '                  #                                                     #',
            '                  #                  All rights reserved                #',
            '                  #                        -***-                        #',
            '                  #######################################################', '', '',
            '                         Program Version %s -  RELEASE  -' % ver, '']


def o_input(name, lines):
    L = ['=' * 80, '                                       INPUT FILE', '=' * 80, 'NAME = %s' % name]
    for i, s in enumerate(lines + ['', '                         ****END OF INPUT****']):
        L.append('|%3d> %s' % (i + 1, s))
    return L + ['=' * 80, '']


def o_coords(Z, X):
    L = ['', '---------------------------------', 'CARTESIAN COORDINATES (ANGSTROEM)',
         '---------------------------------']
    for z, x in zip(Z, X):
        L.append('  %-2s %13.6f%12.6f%12.6f' % (SYM[z], clean(x[0], 5e-7), clean(x[1], 5e-7), clean(x[2], 5e-7)))
    L += ['', '----------------------------', 'CARTESIAN COORDINATES (A.U.)', '----------------------------',
          '  NO LB      ZA    FRAG     MASS         X           Y           Z']
    for i, (z, x) in enumerate(zip(Z, X)):
        b = x / BOHR
        L.append('%4d %-2s %8.4f %4d %9.3f %11.6f %11.6f %11.6f' % (i, SYM[z], z, 0, MASS[z], b[0], b[1], b[2]))
    return L + ['']


def o_settings(charge, mult, nel):
    return ['-' * 20, ' General Settings:', '-' * 20,
            ' Integral files         IntName         .... water',
            ' Hartree-Fock type      HFTyp           .... RHF',
            ' Total Charge           Charge          .... %4d' % charge,
            ' Multiplicity           Mult            .... %4d' % mult,
            ' Number of Electrons    NEL             .... %4d' % nel,
            ' Basis Dimension        Dim             ....   24', '']


def o_scf(e, cycles=9):
    return ['', '               *****************************************************',
            '               *                     SUCCESS                       *',
            '               *           SCF CONVERGED AFTER %3d CYCLES          *' % cycles,
            '               *****************************************************', '',
            '----------------', 'TOTAL SCF ENERGY', '----------------', '',
            'Total Energy       :  %20.8f Eh         %16.5f eV' % (e, e * 27.211386245988), '',
            'Components:', 'Nuclear Repulsion  :            9.16812074 Eh             249.47652 eV', '']


def o_orbitals(occ, virt, nocc_occ=2.0):
    L = ['----------------', 'ORBITAL ENERGIES', '----------------', '',
         '  NO   OCC          E(Eh)            E(eV) ']
    for i, v in enumerate(occ + virt):
        L.append('%4d   %6.4f %16.6f %16.4f ' % (i, nocc_occ if i < len(occ) else 0.0, v, v * 27.211386245988))
    return L + ['']


def o_charges(title, Z, q):
    L = ['-' * len(title), title, '-' * len(title)]
    for i, z in enumerate(Z):
        L.append('%4d %-2s:%12.6f' % (i, SYM[z], q[i]))
    return L + ['Sum of atomic charges:%12.7f' % sum(q), '']


def o_dipole(d):
    mag = math.sqrt(sum(v * v for v in d))
    return ['', '-------------', 'DIPOLE MOMENT', '-------------',
            '                                X             Y             Z',
            'Electronic contribution:  %12.5f  %12.5f  %12.5f' % (0.0, -d[1] * 0.4, d[2] - 0.94282),
            'Nuclear contribution   :  %12.5f  %12.5f  %12.5f' % (0.0, d[1] * 1.4, 0.94282),
            '                        -----------------------------------------',
            'Total Dipole Moment    :  %12.5f  %12.5f  %12.5f' % tuple(d),
            '                        -----------------------------------------',
            'Magnitude (a.u.)       :  %12.5f' % mag, 'Magnitude (Debye)      :  %12.5f' % (mag * 2.541746), '']


def o_fsp(e):
    return ['', '-------------------------   --------------------',
            'FINAL SINGLE POINT ENERGY    %20.12f' % e, '-------------------------   --------------------', '']


def o_cycle(n):
    return ['', '                         ' + '*' * 60,
            '                         *                GEOMETRY OPTIMIZATION CYCLE %3d            *' % n,
            '                         ' + '*' * 60]


def o_conv(vals):
    names = ['Energy change', 'RMS gradient', 'MAX gradient', 'RMS step', 'MAX step']
    tols = [0.0000050000, 0.0001000000, 0.0003000000, 0.0020000000, 0.0040000000]
    L = ['', '                                .--------------------.',
         '          ----------------------|Geometry convergence|-------------------------',
         '          Item                value                   Tolerance       Converged',
         '          ---------------------------------------------------------------------']
    for n, v, t in zip(names, vals, tols):
        if v is not None:
            L.append('          %-15s %16.10f %20.10f      %s' % (n, v, t, 'YES' if abs(v) <= t else 'NO'))
    return L + ['          ........................................................',
                '          Max(Bonds)      0.0046      Max(Angles)    0.35',
                '          Max(Dihed)        0.00      Max(Improp)    0.00',
                '          ---------------------------------------------------------------------', '']


O_HURRAY = ['', '                    ***********************HURRAY********************',
            '                    ***        THE OPTIMIZATION HAS CONVERGED     ***',
            '                    *************************************************', '']
O_FINAL_EVAL = ['', '                       ******************************************************',
                '                       *          FINAL ENERGY EVALUATION AT THE STATIONARY POINT          *',
                '                       ******************************************************']
O_END = ['', '                             ****ORCA TERMINATED NORMALLY****',
         'TOTAL RUN TIME: 0 days 0 hours 0 minutes 12 seconds 345 msec']


def o_freqs(modes, natoms, ir_layout='5', raman=None, depol=None):
    n3 = 3 * natoms
    nz = n3 - len(modes)
    allf = [0.0] * nz + [m['freq'] for m in modes]
    L = ['', '-----------------------', 'VIBRATIONAL FREQUENCIES', '-----------------------', '',
         'Scaling factor for frequencies =  1.000000000  (already applied!)', '']
    for i, f in enumerate(allf):
        L.append('%4d: %12.2f cm**-1' % (i, f) + ('  ***imaginary mode***' if f < 0 else ''))
    L += ['', '', '------------', 'NORMAL MODES', '------------', '',
          'These modes are the Cartesian displacements weighted by the diagonal matrix',
          'M(i,i)=1/sqrt(m[i]) where m[i] is the mass of the displaced atom',
          'Thus, these vectors are normalized but *not* orthogonal', '']
    M = np.zeros((n3, n3))
    for j, m in enumerate(modes):
        M[:, nz + j] = m['disp'].flatten()
    for b in range(0, n3, 6):
        cols = list(range(b, min(b + 6, n3)))
        L.append('          ' + ''.join('%11d' % c for c in cols) + '    ')
        for r in range(n3):
            L.append('%6d     ' % r + ''.join('%11.6f' % clean(M[r, c], 5e-7) for c in cols))
    L += ['', '', '-----------', 'IR SPECTRUM', '-----------', '']
    if ir_layout == '4':
        L += [' Mode    freq (cm**-1)   T**2         TX         TY         TZ',
              '-------------------------------------------------------------------']
        for j, m in enumerate(modes):
            t = m['dmu'] / max(np.linalg.norm(m['dmu']), 1e-30) * math.sqrt(m['ir'])
            L.append('%4d: %12.2f %11.6f  (%11.6f %11.6f %11.6f)' % (nz + j, m['freq'], m['ir'], t[0], t[1], t[2]))
    else:
        L += [' Mode   freq       eps      Int      T**2         TX        TY        TZ',
              '       cm**-1   L/(mol*cm) km/mol    a.u.',
              '----------------------------------------------------------------------------']
        for j, m in enumerate(modes):
            t2 = m['ir'] / 26320.0
            t = m['dmu'] / max(np.linalg.norm(m['dmu']), 1e-30) * math.sqrt(t2)
            L.append('%3d: %9.2f %10.6f %7.2f %9.6f  (%9.6f %9.6f %9.6f)'
                     % (nz + j, m['freq'], m['ir'] / 5053.0, m['ir'], t2, t[0], t[1], t[2]))
        L += ['', '* The epsilon (eps) is given for a Dirac delta lineshape.',
              '** The dipole moment derivative (T) already includes vibrational overlap.']
    L += ['', 'The first frequency considered to be a vibration is %d' % nz,
          'The total number of vibrations considered is %d' % len(modes), '']
    if raman:
        L += ['', '--------------', 'RAMAN SPECTRUM', '--------------', '',
              ' Mode    freq (cm**-1)   Activity   Depolarization',
              '-------------------------------------------------------------------']
        for j, m in enumerate(modes):
            L.append('%4d: %12.2f %13.6f %13.6f' % (nz + j, m['freq'], raman[j], depol[j]))
        L += ['', 'The first frequency considered to be a vibration is %d' % nz, '']
    return L


def o_thermo(e_el, th):
    rt = 3.166811563e-6 * th['temperature']
    ts = th['hCorr'] - th['gCorr']
    thermal = th['eCorr'] - th['zpe']
    k = 627.509474
    return ['', '--------------------------', 'THERMOCHEMISTRY AT 298.15K', '--------------------------', '',
            'Temperature         ...   298.15 K', 'Pressure            ...     1.00 atm',
            'Total Mass          ...    18.02 AMU', '', 'Throughout the following assumptions are being made:',
            '  (1) The electronic state is orbitally nondegenerate', '', '------------', 'INNER ENERGY',
            '------------', '', 'The inner energy is: U= E(el) + E(ZPE) + E(vib) + E(rot) + E(trans)',
            '     E(el)   - is the total energy from the electronic structure calculation', '',
            'Summary of contributions to the inner energy U:',
            'Electronic energy                ... %16.8f Eh' % e_el,
            'Zero point energy                ... %16.8f Eh %10.2f kcal/mol' % (th['zpe'], th['zpe'] * k),
            'Thermal vibrational correction   ... %16.8f Eh %10.2f kcal/mol' % (thermal - 2 * rt, 0.0),
            'Thermal rotational correction    ... %16.8f Eh %10.2f kcal/mol' % (rt, rt * k),
            'Thermal translational correction ... %16.8f Eh %10.2f kcal/mol' % (rt, rt * k),
            '-----------------------------------------------------------------------',
            'Total thermal energy                 %16.8f Eh' % th['eTotal'], '',
            'Summary of corrections to the electronic energy:', '(perhaps to be used in another calculation)',
            'Total thermal correction           %16.8f Eh %10.2f kcal/mol' % (thermal, thermal * k),
            'Non-thermal (ZPE) correction       %16.8f Eh %10.2f kcal/mol' % (th['zpe'], th['zpe'] * k),
            '-----------------------------------------------------------------------',
            'Total correction                   %16.8f Eh %10.2f kcal/mol' % (th['eCorr'], th['eCorr'] * k),
            '', '--------', 'ENTHALPY', '--------', '', 'The enthalpy is H = U + kB*T',
            '                kB*T = %12.8f Eh %10.2f kcal/mol' % (rt, rt * k),
            'Total Enthalpy                    ... %16.8f Eh' % th['H'], '', '-------', 'ENTROPY',
            '-------', '', 'Final entropy term                ... %16.8f Eh %10.2f kcal/mol' % (ts, ts * k), '',
            '-------------------', 'GIBBS FREE ENERGY', '-------------------', '',
            'The Gibbs free energy is G = H - T*S', '',
            'Total enthalpy                    ... %16.8f Eh' % th['H'],
            'Total entropy correction          ... %16.8f Eh %10.2f kcal/mol' % (-ts, -ts * k),
            '-----------------------------------------------------------------------',
            'Final Gibbs free energy         ... %16.8f Eh' % th['G'], '',
            'For completeness - the Gibbs free energy minus the electronic energy',
            'G-E(el)                           ... %16.8f Eh %10.2f kcal/mol' % (th['gCorr'], th['gCorr'] * k), '']


# ----------------------------------------------------------- ORCA fixtures
OW_F = force_field(8.20, -0.06, 0.68, 0.21)
OW_E0 = -76.3179423435
OW_EQ = (0.9643, 102.7)
OW_OCC = [-18.817393, -0.998734, -0.522553, -0.368372, -0.292829]
OW_VIRT = [0.058823, 0.140112, 0.812345, 0.854321, 1.123456]
OW_RAMAN = [2.184578, 87.149025, 21.880203]
OW_DEPOL = [0.723104, 0.172085, 0.750000]
OW_MULL = [-0.371632, 0.185816, 0.185816]
OW_LOEW = [-0.213456, 0.106728, 0.106728]
OW_HIRSH = [-0.312345, 0.156172, 0.156173]
OW_DIP = [0.0, 0.0, -0.82399]


def o_input_lines(keys, X, extra=()):
    return [keys] + list(extra) + ['* xyz 0 1'] + ['  %-2s %12.6f %12.6f %12.6f' % (SYM[z], x[0], x[1], x[2])
                                                    for z, x in zip(GW_Z, X)] + ['*']


def o_cycle_body(X, e, first, mull, dip, conv=None):
    L = o_coords(GW_Z, X)
    if first:
        L += o_settings(0, 1, 10) + ['Your calculation utilizes the basis: def2-SVP', '']
    L += o_scf(e) + o_orbitals(OW_OCC, OW_VIRT)
    L += o_charges('MULLIKEN ATOMIC CHARGES', GW_Z, mull) + o_charges('LOEWDIN ATOMIC CHARGES', GW_Z, OW_LOEW)
    L += o_fsp(e) + o_dipole(dip)
    if conv:
        L += o_conv(conv)
    return L


def make_orca5_water():
    Xe = water_geom(*OW_EQ)
    steps = [(0.9500, 105.0), (0.9672, 102.30), (0.96432, 102.705)]
    geoms = [water_geom(r, t) for r, t in steps]
    energies = [harmonic_energy(X, Xe, OW_F, OW_E0) for X in geoms]
    keys = '! B3LYP def2-SVP Opt Freq'
    L = o_banner('5.0.4') + o_input('water.inp', o_input_lines(keys, geoms[0], ['%pal nprocs 4 end',
                                                                                 '%elprop Polar 1 end']))
    L += ['', '                       *****************************',
          '                       * Geometry Optimization Run *',
          '                       *****************************', '']
    opts = []
    for k, X in enumerate(geoms):
        g = cart_gradient(X, Xe, OW_F)
        dx = ((geoms[k + 1] - X).flatten() / BOHR) if k + 1 < len(geoms) else np.array([0.00011, -0.00007, 0.0])
        conv = [None if k == 0 else energies[k] - energies[k - 1], float(np.sqrt(np.mean(g ** 2))),
                float(np.max(np.abs(g))), float(np.sqrt(np.mean(dx ** 2))), float(np.max(np.abs(dx)))]
        opts.append([None if v is None else r4(v, 10) for v in conv])
        mull = [OW_MULL[0] - 0.01 * (2 - k), OW_MULL[1] + 0.005 * (2 - k), OW_MULL[2] + 0.005 * (2 - k)]
        L += o_cycle(k + 1) + o_cycle_body(X, energies[k], k == 0, mull, OW_DIP, conv)
    L += O_HURRAY + O_FINAL_EVAL + o_cycle_body(geoms[-1], energies[-1], False, OW_MULL, OW_DIP)
    L += ['------------------', 'HIRSHFELD ANALYSIS', '------------------', '',
          'Total integrated alpha density =      4.999999', 'Total integrated beta density  =      4.999999', '',
          '  ATOM     CHARGE      SPIN    ']
    for i, z in enumerate(GW_Z):
        L.append('%4d %-2s %11.6f %11.6f' % (i, SYM[z], OW_HIRSH[i], 0.0))
    L += ['', '  TOTAL  -0.000000    0.000000', '']
    modes = normal_modes(geoms[-1], OW_F)
    L += o_freqs(modes, 3, '5', OW_RAMAN, OW_DEPOL)
    th = thermo_from([m['freq'] for m in modes], energies[-1])
    L += o_thermo(energies[-1], th) + O_END
    write('orca5_water_optfreq.out', '\n'.join(L) + '\n')
    EXPECT['orca5_water'] = dict(
        energies=[r4(e, 12) for e in energies], opt=opts, freqs=[r4(m['freq'], 2) for m in modes],
        ir=[r4(m['ir'], 2) for m in modes], modes=[[r4(clean(v, 5e-7), 6) for v in m['disp'].flatten()]
                                                  for m in modes],
        final=[r4(clean(v, 5e-7), 6) for v in geoms[-1].flatten()],
        thermo={k: r4(v, 8) for k, v in th.items()}, raman=OW_RAMAN, depol=OW_DEPOL)
    return modes, geoms[-1], energies[-1]


def make_orca4_ir():
    X = water_geom(*OW_EQ)
    modes = normal_modes(X, OW_F)
    L = o_banner('4.2.1') + o_input('water_freq.inp', o_input_lines('! B3LYP def2-SVP Freq', X))
    L += o_cycle_body(X, OW_E0, True, OW_MULL, OW_DIP)
    L += o_freqs(modes, 3, '4')
    L += o_thermo(OW_E0, thermo_from([m['freq'] for m in modes], OW_E0)) + O_END
    write('orca4_water_freq.out', '\n'.join(L) + '\n')
    EXPECT['orca4_ir'] = dict(freqs=[r4(m['freq'], 2) for m in modes], ir=[r4(m['ir'], 6) for m in modes])


O_TD_S = [(1, 61944.4, 0.016487919, [('4a', '5a', 0.997345, -0.99867164)]),
          (2, 74847.4, 0.000000000, [('3a', '5a', 0.999048, 0.99952387)]),
          (3, 79943.2, 0.071234567, [('4a', '6a', 0.981234, -0.99057253), ('2a', '5a', 0.012345, 0.11110)])]
O_TD_T = [(1, 56355.1, [('4a', '5a', 0.996123, 0.99805962)]),
          (2, 70693.6, [('3a', '5a', 0.998765, -0.99938231)])]


def o_states(title, states, s2, mult=None):
    L = ['', '-' * len(title), title, '-' * len(title), '',
         'the weight of the individual excitations are printed if larger than 1.0e-02', '']
    for st in states:
        n, cm, tr = st[0], st[1], st[-1]
        line = 'STATE %2d:  E= %10.6f au %10.3f eV %10.1f cm**-1 <S**2> = %11.6f' % (
            n, cm / 219474.6313705, cm / 8065.543937, cm, s2)
        if mult:
            line += ' Mult %d' % mult
        L.append(line)
        for a, b, w, c in tr:
            L.append('%6s -> %4s  : %12.6f (c= %11.8f)' % (a, b, w, c))
        L.append('')
    return L


def make_orca_td(ver):
    X = water_geom(*OW_EQ)
    six = ver.startswith('6')
    extra = ['%tddft nroots 3' + ('' if six else ' triplets true') + ' end']
    L = o_banner(ver) + o_input('water_td.inp', o_input_lines('! B3LYP def2-SVP', X, extra))
    L += o_cycle_body(X, OW_E0, True, OW_MULL, OW_DIP)
    L += ['', '-' * 79, '                          ORCA TD-DFT/TDA CALCULATION', '-' * 79]
    L += o_states('TD-DFT/TDA EXCITED STATES (SINGLETS)', O_TD_S, 0.0, 1 if six else None)
    if not six:
        L += o_states('TD-DFT/TDA EXCITED STATES (TRIPLETS)', O_TD_T, 2.0)
    if six:
        dash = '-' * 100
        L += ['', dash, '                     ABSORPTION SPECTRUM VIA TRANSITION ELECTRIC DIPOLE MOMENTS', dash,
              '     Transition      Energy     Energy  Wavelength fosc(D2)      D2        DX        DY        DZ   ',
              '                      (eV)      (cm-1)    (nm)                 (au**2)    (au)      (au)      (au)  ',
              dash]
        for n, cm, f, _ in O_TD_S:
            L.append('  0-1A  -> %2d-1A %11.6f %9.1f %7.1f %13.9f %9.5f %9.5f %9.5f %9.5f'
                     % (n, cm / 8065.543937, cm, 1e7 / cm, f, f * 5.3, 0.0, -math.sqrt(f * 5.3), 0.0))
        L += ['']
    else:
        for kind, fac in (('ELECTRIC', 1.0), ('VELOCITY', 0.75)):
            dash = '-' * 77
            L += ['', dash, '         ABSORPTION SPECTRUM VIA TRANSITION %s DIPOLE MOMENTS' % kind, dash,
                  'State   Energy    Wavelength  fosc         T2        TX        TY        TZ  ',
                  '        (cm-1)      (nm)                 (au**2)    (au)      (au)      (au) ', dash]
            for n, cm, f, _ in O_TD_S:
                L.append('%4d %9.1f %8.1f %13.9f %9.5f %9.5f %9.5f %9.5f'
                         % (n, cm, 1e7 / cm, f * fac, f * 5.3, 0.0, -math.sqrt(f * 5.3), 0.0))
            L += ['']
    L += O_END
    write('orca%s_water_tddft.out' % ver[0], '\n'.join(L) + '\n')
    EXPECT['orca%s_td' % ver[0]] = dict(
        singlets=[dict(index=s[0], eV=s[1] / 8065.543937, nm=1e7 / s[1], f=s[2],
                       transitions=[[t[0], t[1], t[3]] for t in s[3]]) for s in O_TD_S],
        triplets=[] if six else [dict(index=s[0], eV=s[1] / 8065.543937) for s in O_TD_T])


def make_orca_scan():
    Xe = water_geom(*OW_EQ)
    pts = [[(0.9700, 100.0), (0.9650, 100.0)], [(0.9650, 105.0), (0.9638, 105.0)]]
    X0 = water_geom(*pts[0][0])
    L = o_banner('5.0.3') + o_input('scan.inp', o_input_lines('! B3LYP def2-SVP Opt', X0,
                                                              ['%geom scan A 1 0 2 = 100.0, 105.0, 2 end end']))
    energies, cyc = [], 0
    for p, steps in enumerate(pts):
        L += ['', '         *************************************************************',
              '         *               RELAXED SURFACE SCAN STEP %3d               *' % (p + 1),
              '         *                                                           *',
              '         *   Angle (  1,   0,   2)  : %14.8f                  *' % steps[0][1],
              '         *************************************************************']
        for s, (r, t) in enumerate(steps):
            X = water_geom(r, t)
            e = harmonic_energy(X, Xe, OW_F, OW_E0)
            energies.append(e)
            last = s == len(steps) - 1
            conv = [None if s == 0 else -2.1e-6, 0.00003 if last else 0.0021, 0.00005 if last else 0.0034,
                    0.0004 if last else 0.0051, 0.0007 if last else 0.0089]
            L += o_cycle(s + 1) + o_cycle_body(X, e, p == 0 and s == 0, OW_MULL, OW_DIP, conv)
            if last:
                L += O_HURRAY
    L += ['', '                             *** RELAXED SURFACE SCAN RESULTS ***', '',
          "The Calculated Surface using the 'Actual Energy'",
          '  100.00000000 %16.8f' % energies[1], '  105.00000000 %16.8f' % energies[3], ''] + O_END
    write('orca5_water_scan.out', '\n'.join(L) + '\n')
    EXPECT['orca_scan'] = dict(energies=[r4(e, 12) for e in energies])


# --------------------------------------------------------- other formats
def write_cube(name, titles, origin, axes, dims, atoms, vols, mo_ids=None, angstrom=False):
    L = list(titles)
    na = len(atoms)
    L.append('%5d%12.6f%12.6f%12.6f' % (-na if mo_ids else na, origin[0], origin[1], origin[2]))
    for d, ax in zip(dims, axes):
        L.append('%5d%12.6f%12.6f%12.6f' % (-d if angstrom else d, ax[0], ax[1], ax[2]))
    for z, x in atoms:
        L.append('%5d%12.6f%12.6f%12.6f%12.6f' % (z, float(z), x[0], x[1], x[2]))
    if mo_ids:
        L.append('%5d' % len(mo_ids) + ''.join('%5d' % m for m in mo_ids))
    n1, n2, n3 = dims
    for i in range(n1):
        for j in range(n2):
            row = [v[i, j, k] for k in range(n3) for v in vols]
            for s in range(0, len(row), 6):
                L.append(''.join('%13.5E' % x for x in row[s:s + 6]))
    write(name, '\n'.join(L) + '\n')


def grid(origin, axes, dims):
    i, j, k = np.meshgrid(*[np.arange(d) for d in dims], indexing='ij')
    o, a = np.array(origin), np.array(axes)
    return [o[c] + i * a[0][c] + j * a[1][c] + k * a[2][c] for c in range(3)]


def make_cubes():
    Xb = water_geom(0.9686, 103.65) / BOHR
    atoms = list(zip(GW_Z, Xb))
    exp = {}
    # 1: single MO, Bohr, non-cubic voxels
    o, ax, dims = [-4.0, -4.2, -3.9], [[0.40, 0, 0], [0, 0.45, 0], [0, 0, 0.50]], [20, 18, 16]
    x, y, z = grid(o, ax, dims)
    psi = (z - 0.1) * np.exp(-0.6 * (x ** 2 + y ** 2 + z ** 2))
    write_cube('water_homo_bohr.cube', [' water HOMO (synthetic p-type orbital)', ' MO coefficients'], o, ax, dims,
               atoms, [psi], mo_ids=[5])
    exp['single'] = dict(dims=dims, origin=[v * BOHR for v in o], axes=[[v * BOHR for v in a] for a in ax],
                         probe=[[3, 5, 7], float('%.5E' % psi[3, 5, 7])], last=float('%.5E' % psi[-1, -1, -1]),
                         first=float('%.5E' % psi[0, 0, 0]))
    # 2: two MOs interleaved per point
    o, ax, dims = [-3.0, -3.0, -3.0], [[0.55, 0, 0], [0, 0.6, 0], [0, 0, 0.65]], [12, 11, 10]
    x, y, z = grid(o, ax, dims)
    r2 = x ** 2 + y ** 2 + z ** 2
    m5, m6 = np.exp(-r2 / 2), x * np.exp(-r2 / 2)
    write_cube('water_mo5_mo6.cube', [' water MOs 5 and 6', ' MO coefficients'], o, ax, dims, atoms, [m5, m6],
               mo_ids=[5, 6])
    exp['two'] = dict(dims=dims, probe5=[[2, 3, 4], float('%.5E' % m5[2, 3, 4])],
                      probe6=[[7, 1, 9], float('%.5E' % m6[7, 1, 9])])
    # 3: Angstrom units (negative voxel counts), density
    Xa = water_geom(0.9686, 103.65)
    o, ax, dims = [-1.5, -1.4, -1.2], [[0.3, 0, 0], [0, 0.32, 0], [0.05, 0, 0.3]], [10, 9, 8]
    x, y, z = grid(o, ax, dims)
    rho = np.exp(-(x ** 2 + y ** 2 + z ** 2))
    write_cube('water_density_angstrom.cube', [' water electron density (synthetic)', ' SCF Total Density'], o, ax,
               dims, list(zip(GW_Z, Xa)), [rho], angstrom=True)
    exp['ang'] = dict(dims=dims, origin=o, axes=ax, probe=[[4, 4, 4], float('%.5E' % rho[4, 4, 4])],
                      atom1=[r4(v, 6) for v in Xa[1]])
    EXPECT['cubes'] = exp


def make_molden(modes, X, energies):
    L = ['[Molden Format]', '[Title]', ' water opt+freq (synthetic, from ORCA data)', '[Atoms] Angs']
    for i, (z, x) in enumerate(zip(GW_Z, X)):
        L.append('%-2s %5d %4d %14.8f %14.8f %14.8f' % (SYM[z], i + 1, z, x[0], x[1], x[2]))
    L += ['[GTO]', '  1 0', ' s    3 1.00', '   0.1307093214D+03  0.1543289673D+00',
          '   0.2380886605D+02  0.5353281423D+00', '   0.6443608313D+01  0.4446345422D+00', '',
          '[MO]', ' Sym=     1a', ' Ene= -18.817393', ' Spin= Alpha', ' Occup= 2.000000', '   1   0.994',
          '[GEOCONV]', 'energy'] + ['%18.10f' % e for e in energies]
    L += ['max-force', '  0.0123', '  0.0021', '  0.0001', 'rms-force', '  0.0081', '  0.0012', '  0.00005']
    L += ['[GEOMETRIES] XYZ']
    for k, e in enumerate(energies):
        Xk = X if k == len(energies) - 1 else X * (1 + 0.01 * (len(energies) - 1 - k))
        L += ['%6d' % 3, ' step %d' % (k + 1)]
        for z, x in zip(GW_Z, Xk):
            L.append(' %-2s %14.8f %14.8f %14.8f' % (SYM[z], x[0], x[1], x[2]))
    L += ['[FREQ]'] + ['%12.4f' % m['freq'] for m in modes] + ['[FR-COORD]']
    for z, x in zip(GW_Z, X / BOHR):
        L.append('%-2s %14.8f %14.8f %14.8f' % (SYM[z], x[0], x[1], x[2]))
    L += ['[FR-NORM-COORD]']
    for n, m in enumerate(modes):
        L.append(' vibration %5d' % (n + 1))
        for d in m['disp']:
            L.append('%14.8f %14.8f %14.8f' % (clean(d[0], 5e-9), clean(d[1], 5e-9), clean(d[2], 5e-9)))
    L += ['[INT]'] + ['%12.4f %12.4f' % (m['ir'], OW_RAMAN[n]) for n, m in enumerate(modes)]
    write('water_freq.molden', '\n'.join(L) + '\n')
    EXPECT['molden'] = dict(freqs=[r4(m['freq'], 4) for m in modes], ir=[r4(m['ir'], 4) for m in modes],
                            energies=[r4(e, 10) for e in energies],
                            mode0=[r4(clean(v, 5e-9), 8) for v in modes[0]['disp'].flatten()],
                            final=[r4(v, 8) for v in X.flatten()])


def fchk_int(name, v):
    return ['%-40s   I     %12d' % (name, v)]


def fchk_real(name, v):
    return ['%-40s   R     %22.15E' % (name, v)]


def fchk_iarr(name, vals):
    return ['%-40s   I   N=%12d' % (name, len(vals))] + [''.join('%12d' % v for v in vals[i:i + 6])
                                                         for i in range(0, len(vals), 6)]


def fchk_rarr(name, vals):
    return ['%-40s   R   N=%12d' % (name, len(vals))] + [''.join('%16.8E' % v for v in vals[i:i + 5])
                                                         for i in range(0, len(vals), 5)]


def make_fchk(modes, X, e):
    nm = len(modes)
    e2 = ([m['freq'] for m in modes] + [m['mu'] for m in modes] + [m['k'] for m in modes] +
          [m['ir'] for m in modes] + GW_RAMAN + GW_DEPOL + [0.0] * (8 * nm))
    L = ['Water optimization and frequencies'.ljust(72), '%-10s%-30s%-30s' % ('Freq', 'RB3LYP', '6-31G(d)')]
    L += fchk_int('Number of atoms', 3) + fchk_iarr('Info1-9', [21, 15, 0, 0, 0, 110, 1, 18, -10])
    L += fchk_int('Charge', 0) + fchk_int('Multiplicity', 1) + fchk_int('Number of electrons', 10)
    L += fchk_int('Number of alpha electrons', 5) + fchk_int('Number of beta electrons', 5)
    L += fchk_int('Number of basis functions', 19) + fchk_iarr('Atomic numbers', GW_Z)
    L += fchk_rarr('Nuclear charges', [8.0, 1.0, 1.0])
    L += fchk_rarr('Current cartesian coordinates', list((X / BOHR).flatten()))
    L += fchk_rarr('Real atomic weights', [MASS[z] for z in GW_Z])
    L += fchk_real('SCF Energy', e) + fchk_real('Total Energy', e)
    L += fchk_rarr('Alpha Orbital Energies', GW_OCC + GW_VIRT)
    L += fchk_rarr('Mulliken Charges', GW_MULL)
    L += fchk_rarr('Dipole Moment', [d / 2.541746 for d in GW_DIP])
    L += fchk_int('Vib-NDim', 14) + fchk_int('Vib-NDim0', 14) + fchk_iarr('Vib-AtMap', [1, 2, 3])
    L += fchk_rarr('Vib-E2', e2) + fchk_rarr('Vib-Modes', [v for m in modes for v in m['disp'].flatten()])
    write('water_freq.fchk', '\n'.join(L) + '\n')
    EXPECT['fchk'] = dict(freqs=[float('%.8E' % m['freq']) for m in modes],
                          ir=[float('%.8E' % m['ir']) for m in modes],
                          redMass=[float('%.8E' % m['mu']) for m in modes],
                          xyz=[float('%.8E' % v) * BOHR for v in (X / BOHR).flatten()], energy=e,
                          dipole=[float('%.8E' % (d / 2.541746)) * 2.541746 for d in GW_DIP])


def make_sdf():
    atoms = [('C', -1.3670, 0.0570, 0.0000), ('C', 0.1370, -0.0330, 0.0000), ('O', 0.7840, -1.0610, 0.0000),
             ('O', 0.7200, 1.1830, 0.0000), ('H', -1.7520, -0.9620, 0.0000), ('H', -1.7110, 0.5900, 0.8900),
             ('H', -1.7110, 0.5900, -0.8900), ('H', 1.6800, 1.0200, 0.0000)]
    bonds = [(1, 2, 1), (2, 3, 2), (2, 4, 1), (1, 5, 1), (1, 6, 1), (1, 7, 1), (4, 8, 1)]
    L = ['acetic acid', '  BlueMol11042621003D', 'synthetic test record',
         '%3d%3d  0  0  0  0  0  0  0  0999 V2000' % (len(atoms), len(bonds))]
    for s, x, y, z in atoms:
        L.append('%10.4f%10.4f%10.4f %-3s 0  0  0  0  0  0  0  0  0  0  0  0' % (x, y, z, s))
    for a, b, o in bonds:
        L.append('%3d%3d%3d  0' % (a, b, o))
    L += ['M  END', '> <PUBCHEM_COMPOUND_CID>', '176', '', '$$$$', 'water', '  BlueMol', '',
          '  3  2  0  0  0  0  0  0  0  0999 V2000',
          '    0.0000    0.0000    0.1197 O   0  0  0  0  0  0  0  0  0  0  0  0',
          '    0.0000    0.7614   -0.4789 H   0  0  0  0  0  0  0  0  0  0  0  0',
          '    0.0000   -0.7614   -0.4789 H   0  0  0  0  0  0  0  0  0  0  0  0',
          '  1  2  1  0', '  1  3  1  0', 'M  END', '$$$$']
    write('acetic_acid.sdf', '\n'.join(L) + '\n')
    EXPECT['sdf'] = dict(n=8, bonds=[[a - 1, b - 1, o] for a, b, o in bonds], first=list(atoms[0][1:]))


def pdb_atoms():
    Zc, Xc = cisplatin()
    ala = [('N', 'N', (-0.677, -1.230, -0.491)), ('CA', 'C', (-0.001, 0.064, -0.491)),
           ('C', 'C', (1.499, -0.110, -0.491)), ('O', 'O', (2.030, -1.227, -0.502)),
           ('CB', 'C', (-0.509, 0.856, 0.727)), ('H', 'H', (-0.270, -1.993, -0.990))]
    rows = [('ATOM', n, 'ALA', 1, e, x) for n, e, x in ala]
    rows += [('HETATM', 'O', 'HOH', 101, 'O', (5.0, 5.0, 5.0)), ('HETATM', 'H1', 'HOH', 101, 'H', (5.757, 5.586, 5.0)),
             ('HETATM', 'H2', 'HOH', 101, 'H', (4.243, 5.586, 5.0)), ('HETATM', 'NA', 'NA', 201, 'NA', (-5.0, 0, 0)),
             ('HETATM', 'CL', 'CL', 202, 'CL', (-8.0, 0, 0)), ('HETATM', 'ZN', 'ZN', 203, 'ZN', (0, -6.0, 0))]
    names = ['PT1', 'CL1', 'CL2', 'N1', 'H11', 'H12', 'H13', 'N2', 'H21', 'H22', 'H23']
    for nm, z, x in zip(names, Zc, Xc):
        rows.append(('HETATM', nm, 'CPT', 301, SYM[z].upper(), tuple(x + np.array([10.0, 0, 0]))))
    return rows


def pdb_line(rec, serial, nm4, res, resi, x, el):
    return '%-6s%5d %4s %3s A%4d    %8.3f%8.3f%8.3f%6.2f%6.2f          %2s  ' % (
        rec, serial, nm4, res, resi, x[0], x[1], x[2], 1.0, 0.0, el)


def make_pdbs():
    rows = pdb_atoms()
    L = ['HEADER    BLUEMOLECULE PARSER TEST                04-NOV-26   XXXX              ',
         'CRYST1   30.000   30.000   30.000  90.00  90.00  90.00 P 1           1          ']
    for i, (rec, nm, res, resi, el, x) in enumerate(rows):
        nm4 = nm.ljust(4) if (len(el) == 2 or len(nm) == 4) else (' ' + nm).ljust(4)
        L.append(pdb_line(rec, i + 1, nm4, res, resi, x, el.rjust(2)))
    L += ['TER      %2d      CPT A 301' % (len(rows) + 1), 'END']
    write('peptide_ions_elements.pdb', '\n'.join(L) + '\n')
    # GROMACS-style names, no element columns, two models, plus a Ca2+ ion and an old-style 1HB name
    rows2 = rows[:5] + [('ATOM', '1HB', 'ALA', 1, 'H', (-1.2, 1.5, 0.4))] + rows[6:12] + \
        [('HETATM', 'CA', 'CA', 204, 'CA', (0, 6.0, 0))] + rows[12:]
    L = ['REMARK    GENERATED BY A GROMACS-LIKE WRITER (NO ELEMENT COLUMNS)']
    for m in (1, 2):
        L.append('MODEL     %4d' % m)
        for i, (rec, nm, res, resi, el, x) in enumerate(rows2):
            if nm == 'CA' and res == 'ALA':
                nm4 = 'CA  '                      # left-justified (non-standard but common)
            elif nm == 'CA':
                nm4 = 'CA  '
            else:
                nm4 = nm.ljust(4) if len(nm) == 4 else (' ' + nm).ljust(4)
            xx = tuple(v + 0.1 * (m - 1) for v in x)
            L.append(pdb_line(rec, i + 1, nm4, res, resi, xx, '  ').rstrip())
        L += ['TER', 'ENDMDL']
    L.append('END')
    write('peptide_ions_noelements.pdb', '\n'.join(L) + '\n')
    sym2z = {'H': 1, 'C': 6, 'N': 7, 'O': 8, 'NA': 11, 'CL': 17, 'ZN': 30, 'CA': 20, 'PT': 78}
    EXPECT['pdb'] = dict(Z1=[sym2z[r[4]] for r in rows], Z2=[sym2z[r[4]] for r in rows2],
                         names1=[r[1] for r in rows], res1=[r[2] for r in rows], ids1=[r[3] for r in rows],
                         x_last=list(rows[-1][5]))


def make_gro():
    Zc, Xc = cisplatin()
    rows = [(1, 'SOL', 'OW', (0.126, 1.624, 1.679)), (1, 'SOL', 'HW1', (0.190, 1.661, 1.747)),
            (1, 'SOL', 'HW2', (0.177, 1.568, 1.613)), (2, 'NA', 'NA', (1.000, 1.000, 1.000)),
            (3, 'CL', 'CL', (1.500, 1.000, 1.000))]
    for nm, x in zip(['PT', 'CL1', 'CL2', 'N1', 'H11', 'H12', 'H13', 'N2', 'H21', 'H22', 'H23'], Xc):
        rows.append((4, 'CPT', nm, tuple(x / 10.0 + 1.5)))
    L = []
    for t in (0, 1):
        L += ['BlueMolecule GRO test t= %8.5f step= %d' % (10.0 * t, 5000 * t), '%5d' % len(rows)]
        for i, (ri, rn, an, x) in enumerate(rows):
            v = (0.1227, -0.0580, 0.0434) if i % 2 == 0 else (-0.9045, -2.6469, 1.3180)
            L.append('%5d%-5s%5s%5d%8.3f%8.3f%8.3f%8.4f%8.4f%8.4f' % (ri, rn, an, i + 1, x[0] + 0.01 * t,
                                                                  x[1], x[2], v[0], v[1], v[2]))
        L.append('%10.5f%10.5f%10.5f' % (3.0, 3.1, 3.2))
    write('water_ions_cisplatin.gro', '\n'.join(L) + '\n')
    EXPECT['gro'] = dict(Z=[8, 1, 1, 11, 17] + Zc, first=[r4(v * 10, 2) for v in rows[0][3]],
                         pt=[round(v, 3) * 10 for v in rows[5][3]])


def make_inputs():
    L = ['%chk=methanol.chk', '%mem=2GB', '%nprocshared=4', '#p opt freq b3lyp/6-31g(d)', ' geom=connectivity',
         '', 'Methanol, mixed Cartesian label styles', '', '0 1',
         'C-CA--0.25        -0.046520    0.663200    0.000000',
         'O(Fragment=1)     -0.046520   -0.754800    0.000000',
         'H1                 1.004700    0.980600    0.000000',
         '1         0       -0.549500    1.042200    0.894700',
         'H        -1       -0.549500    1.042200   -0.894700',
         'H4                 0.869500   -1.071900    0.000000', '',
         ' 1 2 1.0 3 1.0 4 1.0 5 1.0', ' 2 6 1.0', ' 3', ' 4', ' 5', ' 6', '', '']
    write('methanol_cartesian.gjf', '\n'.join(L))
    L = ['%chk=nh3.chk', '#p hf/sto-3g opt', '', 'Ammonia Z-matrix with dummy atom and variables', '', '0 1',
         'X', 'N1  1  1.0', 'H1  N1 nh  1  hnx', 'H2  N1 nh  1  hnx  H1  d1', 'H3  N1 nh  1  hnx  H1  -d1',
         'Variables:', 'nh   1.0120', 'hnx  112.15', 'Constants:', 'd1   120.0', '', 'B 2 3 F', '', '']
    write('nh3_zmatrix.com', '\n'.join(L))
    Zc, Xc = cisplatin()
    L = ['%nprocshared=8', '%mem=8GB', '%chk=cisplatin.chk', '#p opt freq b3lyp/genecp scrf=(smd,solvent=water)', '',
         'cisplatin, LANL2DZ on Pt', '', '0 1']
    L += ['%-2s %14.8f %14.8f %14.8f' % (SYM[z], x[0], x[1], x[2]) for z, x in zip(Zc, Xc)]
    L += ['', 'C H N Cl 0', '6-31G(d)', '****', 'Pt 0', 'LANL2DZ', '****', '', 'Pt 0', 'LANL2DZ', '', '']
    write('cisplatin.gjf', '\n'.join(L))
    L = ['%d' % len(Zc), 'cisplatin, cis-[PtCl2(NH3)2]']
    L += ['%-2s %14.8f %14.8f %14.8f' % (SYM[z], x[0], x[1], x[2]) for z, x in zip(Zc, Xc)]
    write('cisplatin.xyz', '\n'.join(L) + '\n')
    EXPECT['cisplatin'] = dict(Z=Zc, xyz=[r4(v, 8) for v in Xc.flatten()])


def make_xyz_traj():
    comments = ['water opt step 1 E = -76.40712345', 'energy: -76.40845678',
                'Coordinates from ORCA-job water E -76.408953', '-76.4089533224']
    L = []
    for k, c in enumerate(comments):
        X = water_geom(0.95 + 0.006 * k, 105.0 - 0.45 * k)
        L += ['3', c]
        for z, x in zip(GW_Z, X):
            lab = str(z) if k == 3 else SYM[z]
            extra = '   0.123  -0.456' if k == 1 else ''
            sep = '\t' if k == 2 else '   '
            L.append(sep.join([lab, '%.6f' % x[0], '%.6f' % x[1], '%.6f' % x[2]]) + extra)
    write('water_traj.xyz', '\r\n'.join(L) + '\r\n')
    EXPECT['xyz'] = dict(energies=[-76.40712345, -76.40845678, -76.408953, -76.4089533224],
                         last=[r4(v, 6) for v in water_geom(0.968, 103.65).flatten()])


def main():
    os.makedirs(FIX, exist_ok=True)
    make_gaussian_water()
    make_gaussian_hpmodes()
    make_gaussian_td()
    make_gaussian_g09_cisplatin()
    make_gaussian_scan()
    omodes, ofinal, oenergy = make_orca5_water()
    make_orca4_ir()
    make_orca_td('5.0.4')
    make_orca_td('6.0.1')
    make_orca_scan()
    make_cubes()
    make_molden(omodes, ofinal, [oenergy + 0.0021, oenergy + 0.0003, oenergy])
    gmodes = normal_modes(water_geom(0.9686, 103.65), GW_F)
    make_fchk(gmodes, water_geom(0.9686, 103.65), GW_E0)
    make_sdf()
    make_pdbs()
    make_gro()
    make_inputs()
    make_xyz_traj()
    with open(os.path.join(FIX, 'expected.json'), 'w') as fh:
        json.dump(EXPECT, fh, indent=1)
    print('fixtures written to', FIX)
    for n in sorted(os.listdir(FIX)):
        print('  %-40s %8d bytes' % (n, os.path.getsize(os.path.join(FIX, n))))


if __name__ == '__main__':
    main()
