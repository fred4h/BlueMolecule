#!/usr/bin/env python3
"""Minimal Hartree-Fock (RHF and UHF) used to make orbital test fixtures for BlueMolecule.

    python3 tools/mini_scf.py            # writes tests/fixtures/orbitals/
    python3 tools/mini_scf.py --check    # also runs self-checks of the integrals

Everything here is written from the textbook formulas:
  * Cartesian Gaussian integrals by the McMurchie-Davidson scheme (Hermite expansion
    coefficients E, Hermite Coulomb integrals R, Boys function by series / asymptotics);
  * real solid harmonics (pure d, f, g) from their explicit polynomial shapes, normalised
    numerically against the Cartesian monomials;
  * closed-shell RHF and open-shell UHF with DIIS.

Only numpy is needed.  Output: Gaussian-style .fchk and Molden files plus, for every file,
a <file>.json holding the basis/MO data a reader should produce from it ("ds") and the
reference values the tests compare against ("expect").
"""
import json
import math
import os
import sys

import numpy as np

BOHR = 0.529177210903          # Angstrom per bohr
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.join(ROOT, 'tests', 'fixtures', 'orbitals')


def dfact(n):
    """n!! with (-1)!! = (0)!! = 1."""
    r = 1
    while n > 1:
        r *= n
        n -= 2
    return r


# --------------------------------------------------------------------------- Cartesian orders
def cart_order(l, fmt):
    """Cartesian components (a, b, c) of a shell in the order the file format lists them."""
    if l == 0:
        return [(0, 0, 0)]
    if l == 1:
        return [(1, 0, 0), (0, 1, 0), (0, 0, 1)]
    if l == 2:
        return [(2, 0, 0), (0, 2, 0), (0, 0, 2), (1, 1, 0), (1, 0, 1), (0, 1, 1)]
    if l == 3:
        return [(3, 0, 0), (0, 3, 0), (0, 0, 3), (1, 2, 0), (2, 1, 0), (2, 0, 1), (1, 0, 2), (0, 1, 2), (0, 2, 1), (1, 1, 1)]
    if l == 4 and fmt == 'molden':
        names = ['xxxx', 'yyyy', 'zzzz', 'xxxy', 'xxxz', 'yyyx', 'yyyz', 'zzzx', 'zzzy', 'xxyy', 'xxzz', 'yyzz', 'xxyz', 'yyxz', 'zzxy']
        return [(s.count('x'), s.count('y'), s.count('z')) for s in names]
    # Gaussian order for g and higher: x power rising, then y power rising (zzzz, yzzz, ..., xxxx)
    return [(a, b, l - a - b) for a in range(l + 1) for b in range(l - a + 1)]


# --------------------------------------------------------------------------- solid harmonics
def _p(*terms):
    d = {}
    for c, k in terms:
        d[k] = d.get(k, 0.0) + c
    return d


def _mul(p, q):
    r = {}
    for k1, v1 in p.items():
        for k2, v2 in q.items():
            k = (k1[0] + k2[0], k1[1] + k2[1], k1[2] + k2[2])
            r[k] = r.get(k, 0.0) + v1 * v2
    return r


def _add(*pairs):
    r = {}
    for c, p in pairs:
        for k, v in p.items():
            r[k] = r.get(k, 0.0) + c * v
    return r


X, Y, Z = _p((1, (1, 0, 0))), _p((1, (0, 1, 0))), _p((1, (0, 0, 1)))
X2, Y2, Z2 = _mul(X, X), _mul(Y, Y), _mul(Z, Z)
R2 = _add((1, X2), (1, Y2), (1, Z2))


def pure_shapes(l):
    """Real solid harmonics in the order m = 0, +1, -1, +2, -2, ... (unnormalised shapes)."""
    if l == 2:
        return [_add((2, Z2), (-1, X2), (-1, Y2)), _mul(X, Z), _mul(Y, Z), _add((1, X2), (-1, Y2)), _mul(X, Y)]
    if l == 3:
        a = _add((4, Z2), (-1, X2), (-1, Y2))
        return [_mul(Z, _add((2, Z2), (-3, X2), (-3, Y2))), _mul(X, a), _mul(Y, a),
                _mul(Z, _add((1, X2), (-1, Y2))), _mul(_mul(X, Y), Z),
                _mul(X, _add((1, X2), (-3, Y2))), _mul(Y, _add((3, X2), (-1, Y2)))]
    if l == 4:
        z4, r4 = _mul(Z2, Z2), _mul(R2, R2)
        b = _add((7, Z2), (-3, R2))
        c = _add((7, Z2), (-1, R2))
        return [_add((35, z4), (-30, _mul(Z2, R2)), (3, r4)),
                _mul(_mul(X, Z), b), _mul(_mul(Y, Z), b),
                _mul(_add((1, X2), (-1, Y2)), c), _mul(_mul(X, Y), c),
                _mul(_mul(X, Z), _add((1, X2), (-3, Y2))), _mul(_mul(Y, Z), _add((3, X2), (-1, Y2))),
                _add((1, _mul(X2, X2)), (-6, _mul(X2, Y2)), (1, _mul(Y2, Y2))),
                _mul(_mul(X, Y), _add((1, X2), (-1, Y2)))]
    raise ValueError('pure shells implemented for l = 2..4')


def _gram(k1, k2):
    a, b, c = k1[0] + k2[0], k1[1] + k2[1], k1[2] + k2[2]
    if a % 2 or b % 2 or c % 2:
        return 0.0
    return float(dfact(a - 1) * dfact(b - 1) * dfact(c - 1))


def pure_matrix(l, comps):
    """T[m, k]: pure function m as a combination of normalised Cartesian functions comps[k]."""
    shapes = pure_shapes(l)
    T = np.zeros((2 * l + 1, len(comps)))
    for m, P in enumerate(shapes):
        q = sum(v1 * v2 * _gram(k1, k2) for k1, v1 in P.items() for k2, v2 in P.items())
        for k, comp in enumerate(comps):
            if comp in P:
                T[m, k] = P[comp] * math.sqrt(_gram(comp, comp) / q)
    return T


# --------------------------------------------------------------------------- basis sets
STO3G_1S = [0.15432897, 0.53532814, 0.44463454]
STO3G_2S = [-0.09996723, 0.39951283, 0.70011547]
STO3G_2P = [0.15591627, 0.60768372, 0.39195739]
BASES = {
    'STO-3G': {
        1: [('s', [3.42525091, 0.62391373, 0.16885540], STO3G_1S)],
        7: [('s', [99.1061690, 18.0523120, 4.8856602], STO3G_1S),
            ('sp', [3.7804559, 0.8784966, 0.2857144], STO3G_2S, STO3G_2P)],
        8: [('s', [130.7093200, 23.8088610, 6.4436083], STO3G_1S),
            ('sp', [5.0331513, 1.1695961, 0.3803890], STO3G_2S, STO3G_2P)],
    },
    '6-31G': {
        1: [('s', [18.7311370, 2.8253937, 0.6401217], [0.03349460, 0.23472695, 0.81375733]),
            ('s', [0.1612778], [1.0])],
        7: [('s', [4173.5110, 627.45790, 142.90210, 40.234330, 12.820210, 4.3904370],
             [0.0018348, 0.0139950, 0.0685870, 0.2322410, 0.4690700, 0.3604550]),
            ('sp', [11.626358, 2.7162800, 0.7722180], [-0.1149610, -0.1691180, 1.1458520], [0.0675800, 0.3239070, 0.7408950]),
            ('sp', [0.2120313], [1.0], [1.0])],
        8: [('s', [5484.6717, 825.23495, 188.04696, 52.964500, 16.897570, 5.7996353],
             [0.0018311, 0.0139501, 0.0684451, 0.2327143, 0.4701930, 0.3585209]),
            ('sp', [15.539616, 3.5999336, 1.0137618], [-0.1107775, -0.1480263, 1.1307670], [0.0708743, 0.3397528, 0.7271586]),
            ('sp', [0.2700058], [1.0], [1.0])],
    },
}
D_POL = {7: 0.8, 8: 0.8}
SYMB = {1: 'H', 6: 'C', 7: 'N', 8: 'O', 26: 'Fe'}
LCHAR = 'spdfghi'


class Shell:
    def __init__(self, atom, center, l, exps, coefs, pure=False, sp=0):
        self.atom, self.center, self.l = atom, np.asarray(center, float), l
        self.exps, self.coefs = np.asarray(exps, float), np.asarray(coefs, float)
        self.pure = bool(pure) and l >= 2
        self.sp = sp          # 1: s half of an SP shell, 2: p half, 0: plain shell
        self.comps = None     # Cartesian components, set per file format

    @property
    def ncart(self):
        return (self.l + 1) * (self.l + 2) // 2

    @property
    def nfunc(self):
        return 2 * self.l + 1 if self.pure else self.ncart


def build_shells(Zs, xyz_bohr, basis, d_pol=None, pure_d=False):
    shells = []
    for ia, z in enumerate(Zs):
        for ent in BASES[basis][z]:
            if ent[0] == 'sp':
                shells.append(Shell(ia, xyz_bohr[ia], 0, ent[1], ent[2], sp=1))
                shells.append(Shell(ia, xyz_bohr[ia], 1, ent[1], ent[3], sp=2))
            else:
                shells.append(Shell(ia, xyz_bohr[ia], 'spdfg'.index(ent[0]), ent[1], ent[2]))
        if d_pol and z in D_POL:
            shells.append(Shell(ia, xyz_bohr[ia], 2, [D_POL[z]], [1.0], pure=pure_d))
    for s in shells:
        s.comps = cart_order(s.l, 'fchk')
    return shells


# --------------------------------------------------------------------------- Boys function
def boys(nmax, T):
    """F_n(T) for n = 0..nmax; returns array of shape (nmax+1,) + T.shape."""
    T = np.asarray(T, float)
    F = np.zeros((nmax + 1,) + T.shape)
    small = T < 33.0
    eT = np.exp(-T)
    if np.any(small):
        Ts = T[small]
        term = np.ones_like(Ts) / (2 * nmax + 1)
        acc = term.copy()
        k = 0
        while True:
            k += 1
            term = term * 2 * Ts / (2 * nmax + 2 * k + 1)
            acc += term
            if np.all(term < 1e-17 * acc) or k > 400:
                break
        Fs = np.zeros((nmax + 1,) + Ts.shape)
        Fs[nmax] = acc * np.exp(-Ts)
        for n in range(nmax, 0, -1):
            Fs[n - 1] = (2 * Ts * Fs[n] + np.exp(-Ts)) / (2 * n - 1)
        F[(slice(None),) + (small,)] = Fs
    big = ~small
    if np.any(big):
        Tb = T[big]
        Fb = np.zeros((nmax + 1,) + Tb.shape)
        Fb[0] = 0.5 * np.sqrt(np.pi / Tb)          # erf(sqrt(T)) == 1 to double precision here
        e = np.exp(-Tb)
        for n in range(nmax):
            Fb[n + 1] = ((2 * n + 1) * Fb[n] - e) / (2 * Tb)
        F[(slice(None),) + (big,)] = Fb
    del eT
    return F


# --------------------------------------------------------------------------- McMurchie-Davidson
def hermite_E(la, lb, a, b, Ax, Bx):
    """E[i, j, t] (Hermite expansion of x_A^i x_B^j exp(..)), vectorised over primitive pairs."""
    p = a + b
    XPA = b * (Bx - Ax) / p
    XPB = a * (Ax - Bx) / p
    E = np.zeros((la + 1, lb + 1, la + lb + 1) + p.shape)
    E[0, 0, 0] = np.exp(-a * b / p * (Ax - Bx) ** 2)
    h = 0.5 / p
    for i in range(la + 1):
        for j in range(lb + 1):
            if i == 0 and j == 0:
                continue
            src, Xp = (E[i - 1, j], XPA) if i > 0 else (E[i, j - 1], XPB)
            top = i + j - 1
            for t in range(i + j + 1):
                v = Xp * src[t] if t <= top else 0.0
                if t > 0:
                    v = v + h * src[t - 1]
                if t + 1 <= top:
                    v = v + (t + 1) * src[t + 1]
                E[i, j, t] = v
    return E


def hermite_R(L, alpha, Xc, Yc, Zc):
    """Hermite Coulomb integrals R_tuv (t+u+v <= L) for exponent alpha and vector P-C."""
    T = alpha * (Xc * Xc + Yc * Yc + Zc * Zc)
    F = boys(L, T)
    Rn = [dict() for _ in range(L + 1)]
    fac = np.ones_like(T)
    for n in range(L + 1):
        Rn[n][(0, 0, 0)] = fac * F[n]
        fac = fac * (-2 * alpha)
    for s in range(1, L + 1):
        for n in range(0, L - s + 1):
            nx = Rn[n + 1]
            for t in range(s + 1):
                for u in range(s - t + 1):
                    v = s - t - u
                    if t > 0:
                        val = Xc * nx[(t - 1, u, v)]
                        if t > 1:
                            val = val + (t - 1) * nx[(t - 2, u, v)]
                    elif u > 0:
                        val = Yc * nx[(t, u - 1, v)]
                        if u > 1:
                            val = val + (u - 1) * nx[(t, u - 2, v)]
                    else:
                        val = Zc * nx[(t, u, v - 1)]
                        if v > 1:
                            val = val + (v - 1) * nx[(t, u, v - 2)]
                    Rn[n][(t, u, v)] = val
    return Rn[0]


def herm_index(L):
    return [(t, u, v) for t in range(L + 1) for u in range(L + 1 - t) for v in range(L + 1 - t - u)]


def prim_weights(sh):
    """Contraction weights for x^a y^b z^c exp(-a r^2) primitives (component factor left out)."""
    a = sh.exps
    return sh.coefs * (2 * a / np.pi) ** 0.75 * (4 * a) ** (sh.l / 2)


class Pair:
    """Bra or ket pair data for two shells (Hermite expansion per Cartesian component pair)."""

    def __init__(self, A, B, extra=0):
        a = A.exps[:, None]
        b = B.exps[None, :]
        self.A, self.B = A, B
        self.p = a + b
        self.P = [(a * A.center[k] + b * B.center[k]) / self.p for k in range(3)]
        self.w = prim_weights(A)[:, None] * prim_weights(B)[None, :]
        self.E = [hermite_E(A.l, B.l + extra, a, b, A.center[k], B.center[k]) for k in range(3)]
        self.L = A.l + B.l
        self.herm = herm_index(self.L)

    def EAB(self, sign=False):
        A, B = self.A, self.B
        out = np.zeros((len(A.comps) * len(B.comps), len(self.herm)) + self.p.shape)
        Ex, Ey, Ez = self.E
        for ia, ca in enumerate(A.comps):
            for ib, cb in enumerate(B.comps):
                row = ia * len(B.comps) + ib
                for h, (t, u, v) in enumerate(self.herm):
                    if t > ca[0] + cb[0] or u > ca[1] + cb[1] or v > ca[2] + cb[2]:
                        continue
                    val = Ex[ca[0], cb[0], t] * Ey[ca[1], cb[1], u] * Ez[ca[2], cb[2], v]
                    out[row, h] = -val if (sign and (t + u + v) % 2) else val
        return out


def one_electron(shells, Zs, xyz):
    """Overlap, kinetic and nuclear-attraction matrices for the (unnormalised-component) Cartesian basis."""
    offs = np.cumsum([0] + [len(s.comps) for s in shells])
    n = offs[-1]
    S, T, V = np.zeros((n, n)), np.zeros((n, n)), np.zeros((n, n))
    for i, A in enumerate(shells):
        for j, B in enumerate(shells[:i + 1]):
            pr = Pair(A, B, extra=2)
            Ex, Ey, Ez = pr.E
            b = B.exps[None, :]
            sq = np.sqrt(np.pi / pr.p)

            def s1(E, i1, j1):
                return E[i1, j1, 0] * sq if j1 >= 0 else 0.0

            def t1(E, i1, j1):
                v = -2 * b * b * s1(E, i1, j1 + 2) + b * (2 * j1 + 1) * s1(E, i1, j1)
                if j1 >= 2:
                    v = v - 0.5 * j1 * (j1 - 1) * s1(E, i1, j1 - 2)
                return v
            ER = pr.EAB()
            Vb = np.zeros(ER.shape[0])
            for ic, z in enumerate(Zs):
                C = xyz[ic]
                R = hermite_R(pr.L, pr.p, pr.P[0] - C[0], pr.P[1] - C[1], pr.P[2] - C[2])
                Rm = np.array([R[h] for h in pr.herm])
                Vb += -z * np.einsum('xhab,hab,ab->x', ER, Rm, 2 * np.pi / pr.p * pr.w)
            for ia, ca in enumerate(A.comps):
                for ib, cb in enumerate(B.comps):
                    sx, sy, sz = s1(Ex, ca[0], cb[0]), s1(Ey, ca[1], cb[1]), s1(Ez, ca[2], cb[2])
                    tx, ty, tz = t1(Ex, ca[0], cb[0]), t1(Ey, ca[1], cb[1]), t1(Ez, ca[2], cb[2])
                    mu, nu = offs[i] + ia, offs[j] + ib
                    S[mu, nu] = S[nu, mu] = np.sum(pr.w * sx * sy * sz)
                    T[mu, nu] = T[nu, mu] = np.sum(pr.w * (tx * sy * sz + sx * ty * sz + sx * sy * tz))
                    V[mu, nu] = V[nu, mu] = Vb[ia * len(B.comps) + ib]
    return S, T, V


def two_electron(shells):
    offs = np.cumsum([0] + [len(s.comps) for s in shells])
    n = offs[-1]
    eri = np.zeros((n, n, n, n))
    ns = len(shells)
    pairs = {}
    for i in range(ns):
        for j in range(i + 1):
            pr = Pair(shells[i], shells[j])
            pairs[(i, j)] = (pr, pr.EAB(), pr.EAB(sign=True))
    keys = sorted(pairs)
    for x, (i, j) in enumerate(keys):
        bra, EB, _ = pairs[(i, j)]
        for (k, l) in keys[:x + 1]:
            ket, _, EK = pairs[(k, l)]
            p = bra.p[:, :, None, None]
            q = ket.p[None, None, :, :]
            alpha = p * q / (p + q)
            PQ = [bra.P[d][:, :, None, None] - ket.P[d][None, None, :, :] for d in range(3)]
            L = bra.L + ket.L
            R = hermite_R(L, alpha, PQ[0], PQ[1], PQ[2])
            Rm = np.empty((len(bra.herm), len(ket.herm)) + alpha.shape)
            for h1, (t, u, v) in enumerate(bra.herm):
                for h2, (t2, u2, v2) in enumerate(ket.herm):
                    Rm[h1, h2] = R[(t + t2, u + u2, v + v2)]
            pref = 2 * np.pi ** 2.5 / (p * q * np.sqrt(p + q)) * bra.w[:, :, None, None] * ket.w[None, None, :, :]
            blk = np.einsum('xTab,TUabcd,abcd,yUcd->xy', EB, Rm, pref, EK, optimize=True)
            A, B, C, D = shells[i], shells[j], shells[k], shells[l]
            blk = blk.reshape(len(A.comps), len(B.comps), len(C.comps), len(D.comps))
            sa, sb, sc, sd = (slice(offs[m], offs[m] + len(shells[m].comps)) for m in (i, j, k, l))
            eri[sa, sb, sc, sd] = blk
            eri[sb, sa, sc, sd] = blk.transpose(1, 0, 2, 3)
            eri[sa, sb, sd, sc] = blk.transpose(0, 1, 3, 2)
            eri[sb, sa, sd, sc] = blk.transpose(1, 0, 3, 2)
            eri[sc, sd, sa, sb] = blk.transpose(2, 3, 0, 1)
            eri[sd, sc, sa, sb] = blk.transpose(3, 2, 0, 1)
            eri[sc, sd, sb, sa] = blk.transpose(2, 3, 1, 0)
            eri[sd, sc, sb, sa] = blk.transpose(3, 2, 1, 0)
    return eri


def cart_to_ao(shells):
    """Block-diagonal T (nao x ncart): normalised Cartesian -> AO (pure for pure shells)."""
    nc = sum(len(s.comps) for s in shells)
    na = sum(s.nfunc for s in shells)
    T = np.zeros((na, nc))
    r = c = 0
    for s in shells:
        k = len(s.comps)
        if s.pure:
            T[r:r + 2 * s.l + 1, c:c + k] = pure_matrix(s.l, s.comps)
            r += 2 * s.l + 1
        else:
            T[r:r + k, c:c + k] = np.eye(k)
            r += k
        c += k
    return T


def integrals(shells, Zs, xyz, need_eri=True):
    S, Tk, V = one_electron(shells, Zs, xyz)
    nrm = 1 / np.sqrt(np.diag(S))       # every Cartesian component normalised on its own
    Tao = cart_to_ao(shells) * nrm[None, :]
    out = {'S': Tao @ S @ Tao.T, 'T': Tao @ Tk @ Tao.T, 'V': Tao @ V @ Tao.T, 'Tao': Tao}
    if need_eri:
        eri = two_electron(shells)
        out['eri'] = np.einsum('pi,qj,rk,sl,ijkl->pqrs', Tao, Tao, Tao, Tao, eri, optimize=True)
    return out


def nuclear_repulsion(Zs, xyz):
    e = 0.0
    for i in range(len(Zs)):
        for j in range(i):
            e += Zs[i] * Zs[j] / np.linalg.norm(xyz[i] - xyz[j])
    return e


# --------------------------------------------------------------------------- SCF
class DIIS:
    def __init__(self, n=8):
        self.F, self.E, self.n = [], [], n

    def push(self, F, err):
        self.F.append(F)
        self.E.append(err)
        if len(self.F) > self.n:
            self.F.pop(0)
            self.E.pop(0)
        m = len(self.F)
        if m < 2:
            return F
        B = -np.ones((m + 1, m + 1))
        B[m, m] = 0
        for i in range(m):
            for j in range(m):
                B[i, j] = np.vdot(self.E[i], self.E[j])
        rhs = np.zeros(m + 1)
        rhs[m] = -1
        try:
            c = np.linalg.solve(B, rhs)[:m]
        except np.linalg.LinAlgError:
            return F
        return sum(ci * Fi for ci, Fi in zip(c, self.F))


def scf(ints, nalpha, nbeta, enuc, unrestricted):
    S, H, eri = ints['S'], ints['T'] + ints['V'], ints['eri']
    w, U = np.linalg.eigh(S)
    Xm = U @ np.diag(w ** -0.5) @ U.T

    def solve(F):
        e, Cp = np.linalg.eigh(Xm.T @ F @ Xm)
        return e, Xm @ Cp
    e, C = solve(H)
    Ca = Cb = C
    diis_a, diis_b = DIIS(), DIIS()
    Eold = 0.0
    for it in range(300):
        Pa = Ca[:, :nalpha] @ Ca[:, :nalpha].T
        Pb = Cb[:, :nbeta] @ Cb[:, :nbeta].T
        J = np.einsum('pqrs,rs->pq', eri, Pa + Pb)
        Fa = H + J - np.einsum('prqs,rs->pq', eri, Pa)
        Fb = H + J - np.einsum('prqs,rs->pq', eri, Pb)
        E = 0.5 * (np.sum((Pa + Pb) * H) + np.sum(Pa * Fa) + np.sum(Pb * Fb)) + enuc
        ea = Fa @ Pa @ S - S @ Pa @ Fa
        eb = Fb @ Pb @ S - S @ Pb @ Fb
        err = max(np.abs(ea).max(), np.abs(eb).max())
        if abs(E - Eold) < 1e-11 and err < 1e-8:
            break
        Eold = E
        if unrestricted:
            ea_, Ca = solve(diis_a.push(Fa, ea))
            eb_, Cb = solve(diis_b.push(Fb, eb))
        else:
            ea_, Ca = solve(diis_a.push(Fa, ea))
            Cb = Ca
    else:
        raise RuntimeError('SCF did not converge')
    ea, Ca = solve(Fa)
    eb, Cb = solve(Fb)
    if not unrestricted:
        eb, Cb = ea, Ca
    for C in (Ca, Cb):                       # deterministic phases: largest coefficient positive
        for m in range(C.shape[1]):
            k = np.argmax(np.abs(C[:, m]))
            if C[k, m] < 0:
                C[:, m] *= -1
    return {'E': E, 'iters': it, 'ea': ea, 'Ca': Ca, 'eb': eb, 'Cb': Cb, 'Pa': Pa, 'Pb': Pb}


# --------------------------------------------------------------------------- evaluation at points
def ao_values(shells, pts):
    """Values of the normalised AOs (rows) at points (bohr), independent of the JS code."""
    cols = []
    for s in shells:
        d = pts - s.center[None, :]
        r2 = np.sum(d * d, axis=1)
        rad = np.zeros(len(pts))
        for a, w in zip(s.exps, prim_weights(s)):
            rad += w * np.exp(-a * r2)
        for c in s.comps:
            cols.append(rad * d[:, 0] ** c[0] * d[:, 1] ** c[1] * d[:, 2] ** c[2])
    return np.array(cols)


def eval_mos(shells, ints, C, pts):
    """C: nao x nmo (columns = MOs) for normalised AOs; returns nmo x npts."""
    nrmT = ints['Tao']                   # includes 1/sqrt(S_cart) normalisation
    return C.T @ (nrmT @ ao_values(shells, pts))


# --------------------------------------------------------------------------- writers
def fmt16(x):
    return '%16.8E' % x


def printed(x, f):
    return float(f % x)


def fchk_text(title, method, basis, Zs, xyz, charge, mult, shells, res, unrestricted, pure_d, pure_f):
    na, nb = res['na'], res['nb']
    nbf = sum(s.nfunc for s in shells)
    lines = ['%-72s' % title[:72], '%-10s%-30s%-30s' % ('SP', method, basis)]

    def I(name, v):
        lines.append('%-40s   I     %12d' % (name, v))

    def R(name, v):
        lines.append('%-40s   R     %22.15E' % (name, v))

    def IA(name, arr):
        lines.append('%-40s   I   N=%12d' % (name, len(arr)))
        for k in range(0, len(arr), 6):
            lines.append(''.join('%12d' % v for v in arr[k:k + 6]))

    def RA(name, arr):
        lines.append('%-40s   R   N=%12d' % (name, len(arr)))
        for k in range(0, len(arr), 5):
            lines.append(''.join(fmt16(v) for v in arr[k:k + 5]))
    # file shells: SP pairs merged
    fsh = []
    k = 0
    while k < len(shells):
        s = shells[k]
        if s.sp == 1:
            p = shells[k + 1]
            fsh.append({'type': -1, 'atom': s.atom, 'exps': s.exps, 'c': s.coefs, 'pc': p.coefs, 'center': s.center})
            k += 2
        else:
            t = s.l if (not s.pure or s.l < 2) else -s.l
            fsh.append({'type': t, 'atom': s.atom, 'exps': s.exps, 'c': s.coefs, 'pc': None, 'center': s.center})
            k += 1
    has_sp = any(f['type'] == -1 for f in fsh)
    I('Number of atoms', len(Zs))
    I('Charge', charge)
    I('Multiplicity', mult)
    I('Number of electrons', na + nb)
    I('Number of alpha electrons', na)
    I('Number of beta electrons', nb)
    I('Number of basis functions', nbf)
    I('Number of independent functions', nbf)
    IA('Atomic numbers', Zs)
    RA('Nuclear charges', [float(z) for z in Zs])
    RA('Current cartesian coordinates', list(np.asarray(xyz).ravel()))
    I('Number of contracted shells', len(fsh))
    I('Number of primitive shells', sum(len(f['exps']) for f in fsh))
    I('Pure/Cartesian d shells', 0 if pure_d else 1)
    I('Pure/Cartesian f shells', 0 if pure_f else 1)
    I('Highest angular momentum', max(abs(f['type']) for f in fsh))
    I('Largest degree of contraction', max(len(f['exps']) for f in fsh))
    IA('Shell types', [f['type'] for f in fsh])
    IA('Number of primitives per shell', [len(f['exps']) for f in fsh])
    IA('Shell to atom map', [f['atom'] + 1 for f in fsh])
    RA('Primitive exponents', [e for f in fsh for e in f['exps']])
    RA('Contraction coefficients', [c for f in fsh for c in f['c']])
    if has_sp:
        RA('P(S=P) Contraction coefficients', [(f['pc'][i] if f['pc'] is not None else 0.0)
                                               for f in fsh for i in range(len(f['exps']))])
    RA('Coordinates of each shell', [x for f in fsh for x in f['center']])
    R('Total Energy', res['E'])
    R('SCF Energy', res['E'])
    RA('Alpha Orbital Energies', list(res['ea']))
    if unrestricted:
        RA('Beta Orbital Energies', list(res['eb']))
    RA('Alpha MO coefficients', list(res['Ca'].T.ravel()))
    if unrestricted:
        RA('Beta MO coefficients', list(res['Cb'].T.ravel()))
    P = res.get('Ptot')
    if P is not None:
        RA('Total SCF Density', [P[i, j] for i in range(nbf) for j in range(i + 1)])
        if unrestricted:
            Ps = res['Pspin']
            RA('Spin SCF Density', [Ps[i, j] for i in range(nbf) for j in range(i + 1)])
    if res.get('mulliken') is not None:
        RA('Mulliken Charges', list(res['mulliken']))
    return '\n'.join(lines) + '\n'


def fchk_ds(shells, res, unrestricted, nbf):
    """What a reader should return for the fchk file: printed values, SP split into s + p."""
    sh = []
    for s in shells:
        sh.append({'atom': s.atom, 'center': [printed(x, '%16.8E') for x in s.center], 'l': s.l, 'pure': s.pure,
                   'exps': [printed(x, '%16.8E') for x in s.exps], 'coefs': [printed(x, '%16.8E') for x in s.coefs]})
    return {'basis': {'source': 'fchk', 'shells': sh, 'nbf': nbf},
            'mos': mos_ds(res, unrestricted, '%16.8E', None)}


def occupations(res, unrestricted):
    n = res['Ca'].shape[1]
    na, nb = res['na'], res['nb']
    if unrestricted:
        return [1.0 if i < na else 0.0 for i in range(n)], [1.0 if i < nb else 0.0 for i in range(n)]
    return [2.0 if i < nb else (1.0 if i < na else 0.0) for i in range(n)], None


def mos_ds(res, unrestricted, cfmt, efmt, syms=None):
    oa, ob = occupations(res, unrestricted)
    efmt = efmt or '%16.8E'

    def spin(e, C, occ):
        return {'energies': [printed(x, efmt) for x in e], 'occupations': occ,
                'symmetries': syms, 'coeffs': [printed(x, cfmt) for x in C.T.ravel()]}
    out = {'restricted': not unrestricted, 'alpha': spin(res['ea'], res['Ca'], oa), 'beta': None}
    if unrestricted:
        out['beta'] = spin(res['eb'], res['Cb'], ob)
    return out


def molden_text(title, Zs, xyz, shells, res, unrestricted, style='generic', cfmt=None):
    """style 'generic': E-format numbers, 'sp' shells kept; 'orca': ORCA-like layout (separate s/p)."""
    L = ['[Molden Format]', '[Title]', ' ' + title, '[Atoms] AU']
    for i, z in enumerate(Zs):
        L.append('%-2s %5d %3d %20.10f %20.10f %20.10f' % (SYMB[z], i + 1, z, xyz[i][0], xyz[i][1], xyz[i][2]))
    L.append('[GTO]')
    ef = '%20.10E' if style == 'generic' else '%20.10f'
    cf = '%20.10E' if style == 'generic' else '%20.10f'
    for ia in range(len(Zs)):
        L.append('%4d 0' % (ia + 1))
        mine = [s for s in shells if s.atom == ia]
        k = 0
        while k < len(mine):
            s = mine[k]
            if s.sp == 1 and style == 'generic':
                p = mine[k + 1]
                L.append(' sp %4d 1.00' % len(s.exps))
                for e, c1, c2 in zip(s.exps, s.wcoefs, p.wcoefs):
                    L.append((ef % e) + (cf % c1) + (cf % c2))
                k += 2
                continue
            L.append(' %s %4d 1.00' % (LCHAR[s.l], len(s.exps)))
            for e, c in zip(s.exps, s.wcoefs):
                L.append((ef % e) + (cf % c))
            k += 1
        L.append('')
    pd = any(s.pure and s.l == 2 for s in shells)
    pf = any(s.pure and s.l == 3 for s in shells)
    pg = any(s.pure and s.l == 4 for s in shells)
    hasf = any(s.l == 3 for s in shells)
    if pd and (pf or not hasf):
        L.append('[5D7F]' if style == 'orca' else '[5D]')
    elif pd:
        L.append('[5D10F]')
    elif pf:
        L.append('[7F]')
    if pg:
        L.append('[9G]')
    L.append('[MO]')
    oa, ob = occupations(res, unrestricted)
    cfmt = cfmt or '%20.12E'
    for spin, e, C, occ in (('Alpha', res['ea'], res['Cmolden_a'], oa), ('Beta', res['eb'], res['Cmolden_b'], ob)):
        if spin == 'Beta' and not unrestricted:
            break
        for m in range(C.shape[1]):
            L.append(' Sym= %8s' % ('%da' % (m + 1)))
            L.append(' Ene= %20.10E' % e[m])
            L.append(' Spin= %s' % spin)
            L.append(' Occup= %10.6f' % occ[m])
            for mu in range(C.shape[0]):
                L.append('%5d ' % (mu + 1) + (cfmt % C[mu, m]).strip().rjust(20))
    return '\n'.join(L) + '\n'


def molden_ds(shells, res, unrestricted, nbf, style, cfmt):
    ef = '%20.10E' if style == 'generic' else '%20.10f'
    cf = '%20.10E' if style == 'generic' else '%20.10f'
    sh = []
    for s in shells:
        sh.append({'atom': s.atom, 'center': [printed(x, '%20.10f') for x in s.center], 'l': s.l, 'pure': s.pure,
                   'exps': [printed(x, ef) for x in s.exps], 'coefs': [printed(x, cf) for x in s.wcoefs]})
    oa, ob = occupations(res, unrestricted)

    def spin(e, C, occ):
        return {'energies': [printed(x, '%20.10E') for x in e], 'occupations': [printed(x, '%10.6f') for x in occ],
                'symmetries': ['%da' % (m + 1) for m in range(C.shape[1])],
                'coeffs': [printed(x, cfmt) for x in C.T.ravel()]}
    mos = {'restricted': not unrestricted, 'alpha': spin(res['ea'], res['Cmolden_a'], oa), 'beta': None}
    if unrestricted:
        mos['beta'] = spin(res['eb'], res['Cmolden_b'], ob)
    return {'basis': {'source': 'molden', 'shells': sh, 'nbf': nbf}, 'mos': mos}


# --------------------------------------------------------------------------- convention variants
def axis_factors(shells):
    """Per-AO factor f with C_axis = C_component * f (Cartesian l >= 2 shells)."""
    f = []
    for s in shells:
        if s.pure:
            f += [1.0] * (2 * s.l + 1)
        else:
            for c in s.comps:
                f.append(math.sqrt(dfact(2 * s.l - 1) / (dfact(2 * c[0] - 1) * dfact(2 * c[1] - 1) * dfact(2 * c[2] - 1))))
    return np.array(f)


def orca_signs(shells):
    f = []
    for s in shells:
        if s.pure:
            ms = [0] + [m for k in range(1, s.l + 1) for m in (k, -k)]
            f += [-1.0 if (s.l >= 3 and abs(m) in (3, 4)) else 1.0 for m in ms]
        else:
            f += [1.0] * len(s.comps)
    return np.array(f)


def raw_coefs(s):
    """Contraction coefficients that include the primitive normalisation (x^l component) and the
    contraction normalisation, as some programs print them."""
    a = s.exps
    n = (2 * a / np.pi) ** 0.75 * (4 * a) ** (s.l / 2) / math.sqrt(dfact(2 * s.l - 1))
    c = s.coefs * n
    # contraction normalisation for the x^l component
    ov = 0.0
    for i in range(len(a)):
        for j in range(len(a)):
            ov += s.coefs[i] * s.coefs[j] * (2 * math.sqrt(a[i] * a[j]) / (a[i] + a[j])) ** (s.l + 1.5)
    return c / math.sqrt(ov)


def renorm_contraction(s):
    a = s.exps
    ov = sum(s.coefs[i] * s.coefs[j] * (2 * math.sqrt(a[i] * a[j]) / (a[i] + a[j])) ** (s.l + 1.5)
             for i in range(len(a)) for j in range(len(a)))
    return math.sqrt(ov)


def permute_to_molden(shells, C):
    """Reorder Cartesian g (and higher) rows of C (nao x nmo) from fchk to Molden order; returns C copy."""
    rows = []
    off = 0
    for s in shells:
        n = s.nfunc
        if not s.pure and s.l >= 4:
            src = cart_order(s.l, 'fchk')
            dst = cart_order(s.l, 'molden')
            rows += [off + src.index(c) for c in dst]
        else:
            rows += list(range(off, off + n))
        off += n
    return C[rows, :]


# --------------------------------------------------------------------------- jobs
def sample_points(xyz, n=24, seed=7):
    rng = np.random.default_rng(seed)
    pts = []
    for k in range(n):
        a = xyz[k % len(xyz)]
        pts.append(a + rng.normal(scale=1.2, size=3))
    return np.array(pts)


def save(name, text=None, ds=None, expect=None):
    if text is not None:
        with open(os.path.join(OUT, name), 'w') as f:
            f.write(text)
    with open(os.path.join(OUT, name + '.json'), 'w') as f:
        json.dump({'file': name, 'ds': ds, 'expect': expect}, f, separators=(',', ':'))


def mulliken(res, S, shells, Zs):
    P = res['Ptot']
    PS = P @ S
    q = np.array(Zs, float)
    off = 0
    for s in shells:
        for k in range(s.nfunc):
            q[s.atom] -= PS[off + k, off + k]
        off += s.nfunc
    return q


def run_molecule(tag, Zs, xyz_ang, basis, d_pol, pure_d, charge, mult, variants, log):
    xyz = np.array(xyz_ang, float) / BOHR
    shells = build_shells(Zs, xyz, basis, d_pol=d_pol, pure_d=pure_d)
    nel = sum(Zs) - charge
    na, nb = (nel + mult - 1) // 2, (nel - mult + 1) // 2
    unres = mult != 1
    ints = integrals(shells, Zs, xyz)
    enuc = nuclear_repulsion(Zs, xyz)
    res = scf(ints, na, nb, enuc, unres)
    res['na'], res['nb'] = na, nb
    nbf = ints['S'].shape[0]
    Pa = res['Ca'][:, :na] @ res['Ca'][:, :na].T
    Pb = res['Cb'][:, :nb] @ res['Cb'][:, :nb].T
    res['Ptot'], res['Pspin'] = Pa + Pb, Pa - Pb
    res['mulliken'] = mulliken(res, ints['S'], shells, Zs)
    dev = max(np.abs(res['Ca'].T @ ints['S'] @ res['Ca'] - np.eye(nbf)).max(),
              np.abs(res['Cb'].T @ ints['S'] @ res['Cb'] - np.eye(nbf)).max())
    log.append('%-22s %-7s E = %.8f  (%d iterations, nbf %d, |CSC-1| %.1e)' % (tag, basis + ('*' if d_pol else ''), res['E'], res['iters'], nbf, dev))
    pts = sample_points(xyz)
    vals = {'alpha': eval_mos(shells, ints, res['Ca'], pts)}
    if unres:
        vals['beta'] = eval_mos(shells, ints, res['Cb'], pts)
    homo_a = na - 1
    expect = {'energy': res['E'], 'nelec': nel, 'nalpha': na, 'nbeta': nb, 'homo': homo_a, 'nbf': nbf,
              'points_bohr': pts.tolist(),
              'values': {sp: {str(m): vals[sp][m].tolist() for m in sorted({0, 1, homo_a, homo_a + 1, nbf - 1})}
                         for sp in vals}}
    for v in variants:
        name, kind = v['file'], v['kind']
        exp = dict(expect)
        exp['convention'] = v['convention']
        if kind == 'fchk':
            text = fchk_text(tag, ('U' if unres else 'R') + 'HF', basis + ('*' if d_pol else ''), Zs, xyz, charge, mult,
                             shells, res, unres, pure_d, pure_d)
            save(name, text, fchk_ds(shells, res, unres, nbf), exp)
        else:
            style = v.get('style', 'generic')
            fac = np.ones(nbf)
            if v.get('axis'):
                fac = fac * axis_factors(shells)
            if v.get('orca_signs'):
                fac = fac * orca_signs(shells)
            res['Cmolden_a'] = res['Ca'] * fac[:, None]
            res['Cmolden_b'] = res['Cb'] * fac[:, None]
            for s in shells:
                s.wcoefs = raw_coefs(s) if v.get('raw') else s.coefs
            cfmt = v.get('cfmt', '%20.12E')
            text = molden_text(tag, Zs, xyz, shells, res, unres, style=style, cfmt=cfmt)
            save(name, text, molden_ds(shells, res, unres, nbf, style, cfmt), exp)
    return res, ints, shells


def random_fg(log):
    """Basis with pure / Cartesian d, f and g shells and Loewdin-orthonormalised random MOs."""
    Zs = [26, 8, 1]
    xyz = np.array([[0.0, 0.0, 0.0], [0.4, -0.3, 3.1], [2.2, 1.1, -1.0]])
    spec = {
        0: [(0, [9.5, 1.3], [0.4, 0.7]), (1, [6.0, 0.9], [0.5, 0.6]), (2, [4.1, 1.2, 0.35], [0.2, 0.5, 0.5]),
            (2, [0.6], [1.0]), (3, [2.4, 0.7], [0.45, 0.65]), (4, [1.5], [1.0])],
        1: [(0, [12.0, 2.5, 0.5], [0.15, 0.55, 0.5]), (1, [3.3, 0.6], [0.4, 0.7]), (2, [1.1], [1.0]), (3, [0.9], [1.0])],
        2: [(0, [5.0, 0.8], [0.3, 0.8]), (1, [0.75], [1.0])],
    }
    rng = np.random.default_rng(20261105)
    out = {}
    for pure in (True, False):
        shells = []
        for ia in range(3):
            for l, e, c in spec[ia]:
                s = Shell(ia, xyz[ia], l, e, c, pure=pure)
                s.comps = cart_order(l, 'fchk')
                shells.append(s)
        ints = integrals(shells, Zs, xyz, need_eri=False)
        S = ints['S']
        n = S.shape[0]
        C0 = rng.normal(size=(n, n))
        M = C0.T @ S @ C0
        w, U = np.linalg.eigh(M)
        C = C0 @ (U @ np.diag(w ** -0.5) @ U.T)       # columns orthonormal: C^T S C = 1
        dev = np.abs(C.T @ S @ C - np.eye(n)).max()
        nocc = 12
        res = {'E': -1234.5, 'ea': np.linspace(-3.0, 2.0, n), 'eb': np.linspace(-3.0, 2.0, n), 'Ca': C, 'Cb': C,
               'na': nocc, 'nb': nocc, 'Ptot': None}
        pts = sample_points(xyz, seed=11)
        vals = eval_mos(shells, ints, C, pts)
        expect = {'energy': None, 'nelec': 2 * nocc, 'nalpha': nocc, 'nbeta': nocc, 'homo': nocc - 1, 'nbf': n,
                  'points_bohr': pts.tolist(),
                  'values': {'alpha': {str(m): vals[m].tolist() for m in (0, 3, nocc - 1, nocc, n - 1)}}}
        log.append('random f/g basis %-9s nbf %d, |CSC-1| %.1e' % ('pure' if pure else 'Cartesian', n, dev))
        tag = 'fg_pure' if pure else 'fg_cart'
        exp = dict(expect)
        exp['convention'] = {'prim': 'normalized', 'renorm': True, 'cart': 'component', 'sign': 'standard'}
        save(tag + '.fchk', fchk_text(tag + ' random orthonormal MOs', 'RHF', 'Gen', Zs, xyz, 0, 1, shells, res, False, pure, pure),
             fchk_ds(shells, res, False, n), exp)
        for s in shells:
            s.comps = cart_order(s.l, 'molden')
        Cm = permute_to_molden(shells, C)
        variants = [(tag + '.molden', {}, exp['convention'])]
        if pure:
            variants.append((tag + '_orca.molden', {'raw': True, 'orca_signs': True, 'style': 'orca', 'cfmt': '%12.6f'},
                             {'prim': 'raw', 'renorm': True, 'sign': 'orca'}))
        else:
            variants.append((tag + '_axis.molden', {'axis': True}, {'prim': 'normalized', 'renorm': True, 'cart': 'axis'}))
        for fname, opt, conv in variants:
            fac = np.ones(n)
            if opt.get('axis'):
                fac = fac * axis_factors(shells)
            if opt.get('orca_signs'):
                fac = fac * orca_signs(shells)
            res['Cmolden_a'] = res['Cmolden_b'] = Cm * fac[:, None]
            for s in shells:
                s.wcoefs = raw_coefs(s) if opt.get('raw') else s.coefs
            e2 = dict(expect)
            e2['convention'] = conv
            style = opt.get('style', 'generic')
            cfmt = opt.get('cfmt', '%20.12E')
            save(fname, molden_text(tag, Zs, xyz, shells, res, False, style=style, cfmt=cfmt),
                 molden_ds(shells, res, False, n, style, cfmt), e2)
        out[tag] = (shells, ints, C)
    return out


# --------------------------------------------------------------------------- self checks
def self_check(log):
    """Integrals against brute-force quadrature and the pure transform against orthonormality."""
    xyz = np.array([[0.0, 0.0, 0.0], [0.3, -0.4, 1.1]])
    sh = [Shell(0, xyz[0], 2, [0.9, 0.3], [0.6, 0.5]), Shell(1, xyz[1], 3, [0.7], [1.0]), Shell(1, xyz[1], 1, [1.4], [1.0])]
    for s in sh:
        s.comps = cart_order(s.l, 'fchk')
    ints = integrals(sh, [3, 2], xyz, need_eri=False)
    g = np.linspace(-9, 9, 81)
    h = g[1] - g[0]
    P = np.stack(np.meshgrid(g, g, g, indexing='ij'), axis=-1).reshape(-1, 3)
    A = ints['Tao'] @ ao_values(sh, P)
    Sq = (A @ A.T) * h ** 3
    err = np.abs(Sq - ints['S']).max()
    log.append('self-check: overlap vs quadrature max error %.1e' % err)
    assert err < 1e-6, err
    for l in (2, 3, 4):
        comps = cart_order(l, 'fchk')
        s1 = Shell(0, [0, 0, 0], l, [1.0], [1.0], pure=True)
        s1.comps = comps
        S1 = integrals([s1], [1], np.zeros((1, 3)), need_eri=False)['S']
        e = np.abs(S1 - np.eye(2 * l + 1)).max()
        log.append('self-check: pure l=%d orthonormality %.1e' % (l, e))
        assert e < 1e-12


def main():
    os.makedirs(OUT, exist_ok=True)
    log = []
    if '--check' in sys.argv:
        self_check(log)
    h2o = [8, 1, 1]
    h2o_xyz = [[0.0, 0.0, 0.1173], [0.0, 0.7572, -0.4692], [0.0, -0.7572, -0.4692]]
    nh3 = [7, 1, 1, 1]
    nh3_xyz = [[0.0, 0.0, 0.1162], [0.0, 0.9397, -0.2711], [0.8138, -0.4699, -0.2711], [-0.8138, -0.4699, -0.2711]]
    nh2 = [7, 1, 1]
    nh2_xyz = [[0.0, 0.0, 0.1430], [0.0, 0.8015, -0.5006], [0.0, -0.8015, -0.5006]]
    std = {'prim': 'normalized', 'renorm': True}
    run_molecule('h2o_sto3g', h2o, h2o_xyz, 'STO-3G', False, False, 0, 1, [
        {'file': 'h2o_sto3g.fchk', 'kind': 'fchk', 'convention': std},
        {'file': 'h2o_sto3g.molden', 'kind': 'molden', 'convention': std}], log)
    run_molecule('nh3_631g', nh3, nh3_xyz, '6-31G', False, False, 0, 1, [
        {'file': 'nh3_631g.fchk', 'kind': 'fchk', 'convention': std},
        {'file': 'nh3_631g.molden', 'kind': 'molden', 'convention': std}], log)
    run_molecule('h2o_631gd_cart', h2o, h2o_xyz, '6-31G', True, False, 0, 1, [
        {'file': 'h2o_631gd_cart.fchk', 'kind': 'fchk', 'convention': dict(std, cart='component')},
        {'file': 'h2o_631gd_cart_axis.molden', 'kind': 'molden', 'axis': True, 'convention': dict(std, cart='axis')}], log)
    run_molecule('h2o_631gd_pure', h2o, h2o_xyz, '6-31G', True, True, 0, 1, [
        {'file': 'h2o_631gd_pure.fchk', 'kind': 'fchk', 'convention': std},
        {'file': 'h2o_631gd_pure.molden', 'kind': 'molden', 'convention': std},
        {'file': 'h2o_631gd_pure_orca.molden', 'kind': 'molden', 'raw': True, 'orca_signs': True, 'style': 'orca',
         'cfmt': '%12.6f', 'convention': {'prim': 'raw', 'renorm': True}}], log)
    run_molecule('nh2_uhf_631g', nh2, nh2_xyz, '6-31G', False, False, 0, 2, [
        {'file': 'nh2_uhf_631g.fchk', 'kind': 'fchk', 'convention': std},
        {'file': 'nh2_uhf_631g.molden', 'kind': 'molden', 'convention': std}], log)
    random_fg(log)
    print('\n'.join(log))
    print('wrote fixtures to', OUT)


if __name__ == '__main__':
    main()
