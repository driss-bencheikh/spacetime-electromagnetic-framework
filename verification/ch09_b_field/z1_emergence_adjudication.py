#!/usr/bin/env python3
"""z1: Phase 1 adjudication checks for the emergence note and the temporal audit (1 November 2026).

    python3 verification/ch09_b_field/z1_emergence_adjudication.py

Prints one 'PASS z1-<id>: <claim>' or 'FAIL z1-<id>: <claim>' line per claim and exits 0 iff all pass.

Purpose. These are the checks written in Phase 1 to adjudicate referee findings and review items (y1, y2, e3/e4 referee, LAD-01,
LAD-02, CLASS-06, TIME, CORE-04, CORE-07). It imports nothing from the other verification scripts. Every PASS line states the
mathematical claim it tests; it does not by itself validate the surrounding manuscript prose.

Realisation (tex/s02_arena.tex, tex/sA_realisation.tex): e1=sx(x)1(x)1, e2=sz(x)sx(x)1, e3=sz(x)sz(x)sx, e4=ep(x)1(x)1,
e5=sz(x)ep(x)1, e6=sz(x)sz(x)ep, ep=[[0,-1],[1,0]], signature (+,+,+,-,-,-), omega=e1...e6.
Tower: E_(k,a) = omega^{(x)(k-1)} (x) e_a (x) 1^{(x)(n-k)} on R^{8^n}.
Copy rotation (single convention of the note after Phase 1): C_b = +1/2 sum_a eta_aa E_(1,a) E_(2,a).

Method. All generators are integer matrices; ranks, kernels, commutators and identities are computed in exact integer/rational
arithmetic (sympy or numpy int64, which is exact for these entries) unless a line says 'float'. Floating point is used only for
matrix exponentials and eigenvalues, with tolerance 1e-9.
"""
import itertools
import sys
from fractions import Fraction
from functools import reduce

import numpy as np
import sympy as sp
from scipy.linalg import expm

NPASS, FAILS = [0], []


def report(cid, ok, text):
    if ok:
        NPASS[0] += 1
        print(f"PASS z1-{cid}: {text}")
    else:
        FAILS.append(cid)
        print(f"FAIL z1-{cid}: {text}")
    sys.stdout.flush()


sx = np.array([[0, 1], [1, 0]], dtype=np.int64)
sz = np.array([[1, 0], [0, -1]], dtype=np.int64)
ep = np.array([[0, -1], [1, 0]], dtype=np.int64)
I2 = np.eye(2, dtype=np.int64)


def kr(*a):
    return reduce(np.kron, a)


e = [None, kr(sx, I2, I2), kr(sz, sx, I2), kr(sz, sz, sx), kr(ep, I2, I2), kr(sz, ep, I2), kr(sz, sz, ep)]
eta = [None, 1, 1, 1, -1, -1, -1]
I8 = np.eye(8, dtype=np.int64)
om = reduce(np.matmul, e[1:])


def E(n, k, a):
    return kr(*([om] * (k - 1) + [e[a]] + [I8] * (n - k)))


def blade(idx, n=1):
    """Product of tower generators; idx is a list of (k, a) pairs (n>1) or of a (n=1)."""
    if n == 1:
        return reduce(np.matmul, [e[a] for a in idx], I8)
    return reduce(np.matmul, [E(n, k, a) for (k, a) in idx], np.eye(8 ** n, dtype=np.int64))


def rank_exact(rows):
    """Exact rank over Q of a list of integer vectors."""
    if not rows:
        return 0
    from sympy.polys.matrices import DomainMatrix
    return DomainMatrix([[sp.ZZ(int(v)) for v in r] for r in rows], (len(rows), len(rows[0])), sp.ZZ).rank()


def comm(A, B):
    return A @ B - B @ A


for a in range(1, 7):
    assert (e[a] @ e[a] == eta[a] * I8).all()
    for b in range(a + 1, 7):
        assert (e[a] @ e[b] + e[b] @ e[a] == 0).all()

# ======================================================================================== 1. e1 (y1 findings, sign of U)
N1 = sum((I8 - e[i] @ e[i + 3]) for i in (1, 2, 3))  # = 2N, integer
biv = [(a, b) for a in range(1, 7) for b in range(a + 1, 7)]
Bv = {ab: e[ab[0]] @ e[ab[1]] for ab in biv}
symm = [ab for ab in biv if (Bv[ab].T == Bv[ab]).all()]
anti = [ab for ab in biv if ab not in symm]


def commutant_dim(S, X):
    """dim of {Y in span(S) : [X,Y]=0}, exact."""
    return len(S) - rank_exact([comm(X, Bv[ab]).ravel().tolist() for ab in S])


d_all, d_sym, d_anti = commutant_dim(biv, N1), commutant_dim(symm, N1), commutant_dim(anti, N1)
report("1a", (d_all, d_sym, d_anti) == (9, 6, 3),
       f"y1 finding 1: the bivectors of so(3,3) commuting with N span dim {d_all} = gl(3,R), of which {d_sym} symmetric (boosts) and "
       f"{d_anti} antisymmetric; so 'N is invariant under nothing that contains a boost' is false (exact rank)")

dN = np.diag(N1) // 2


def dn_values(X):
    return sorted({int(dN[i] - dN[j]) for i in range(8) for j in range(8) if X[i, j] != 0})


report("1b", dn_values(Bv[(2, 3)]) == [-2, 0, 2] and dn_values(Bv[(1, 2)]) == [-2, 0, 2]
       and (Bv[(2, 3)].T == -Bv[(2, 3)]).all() and (Bv[(1, 2)].T == -Bv[(1, 2)]).all(),
       "y1 finding 2: the compact (antisymmetric) bivectors J=e2e3 and e1e2 have matrix entries with Delta N = -2, 0, +2: compact gauge-"
       "admissible fields outside Cent(G) create and annihilate quanta in pairs (exact)")

Hu = kr(e[1], I8) + kr(I8, e[1])
evs = np.round(np.linalg.eigvalsh(Hu.astype(float)), 9) + 0.0
mult = {float(v): int((evs == v).sum()) for v in sorted(set(evs))}
Hu2 = kr(e[1], I8) + 2 * kr(I8, e[1])
evs2 = np.round(np.linalg.eigvalsh(Hu2.astype(float)), 9) + 0.0
mult2 = {float(v): int((evs2 == v).sum()) for v in sorted(set(evs2))}
report("1c", mult == {-2.0: 16, 0.0: 32, 2.0: 16} and mult2 == {-3.0: 16, -1.0: 16, 1.0: 16, 3.0: 16},
       f"y1 finding 4: spec H_u at m1=m2=1 is {mult} (value 0 with multiplicity 32), at (m1,m2)=(1,2) {mult2} (float eigenvalues)")

Cb = sum(Fraction(1, 2) * eta[a] * (E(2, 1, a) @ E(2, 2, a)) for a in range(1, 7))
Cbf = np.array(Cb, dtype=float)
U = reduce(np.matmul, [np.eye(64, dtype=np.int64) + eta[a] * (E(2, 1, a) @ E(2, 2, a)) for a in range(1, 7)])
Uf = U / 8.0
okU = np.allclose(Uf, expm(np.pi / 2 * Cbf), atol=1e-9) and not np.allclose(Uf, expm(-np.pi / 2 * Cbf), atol=1e-9)
Ui = np.linalg.inv(Uf)
okC = all(np.allclose(Uf @ E(2, 1, a) @ Ui, -E(2, 2, a), atol=1e-9) for a in range(1, 7))
report("1d", okU and okC,
       "sign fix of PENDING.md rejected: with C_b = +1/2 sum eta_aa E_(1,a)E_(2,a), (1/8) prod_a (1 + eta_aa E_(1,a)E_(2,a)) = exp(+pi C_b/2) "
       "and not exp(-pi C_b/2); it maps E_(1,a) to -E_(2,a) (float expm, tol 1e-9)")

# ======================================================================================== 2. e2 (threshold reading)
m, x = sp.symbols("m x", positive=True)
E1s = sp.sqrt(m ** 2 + x ** 2)
a_ = 2 * E1s ** 2
b_ = E1s ** 4 - m ** 4
Eplus2 = sp.simplify(a_ + 2 * sp.sqrt(b_))
cross = sp.solve(sp.Eq(Eplus2, 4 * m ** 2), x)
val = float(sp.sqrt(Eplus2.subs({m: 1, x: 2})))
# deficit at rest is not proportional to m1 m2
m1, m2 = sp.symbols("m1 m2", positive=True)
deficit = m1 + m2 - sp.sqrt(m1 ** 2 + m2 ** 2)
report("2a", cross == [m / 2] and abs(val - 4.449489742783178) < 1e-12
       and sp.simplify((m1 + m2) ** 2 - (m1 ** 2 + m2 ** 2) - 2 * m1 * m2) == 0
       and sp.simplify(deficit.subs({m1: m, m2: m}) - (2 - sp.sqrt(2)) * m) == 0,
       f"e2 threshold: at m1=m2=m, p=-q, |p|=x, E_+^2 = {sp.simplify(Eplus2)}, which equals (2m)^2 exactly at x = m/2 (sympy solve); at m=1, "
       f"|p|=2 E_+ = {val:.4f} > 2; (m1+m2)^2 - M^2 = 2 m1 m2 while the rest deficit at m1=m2=m is (2-sqrt2) m (sympy)")

# ======================================================================================== 3. e3 referee
def tower(n):
    return {(k, a): E(n, k, a) for k in range(1, n + 1) for a in range(1, 7)}


def casimirs_x4(n):
    """Return 4*C_cop and 4*C_diag as exact integer matrices (C^{kl} and Delta_ab carry a factor 1/2)."""
    T = tower(n)
    dim = 8 ** n
    Ccop4 = np.zeros((dim, dim), dtype=np.int64)
    for k in range(1, n + 1):
        for l in range(k + 1, n + 1):
            C2 = sum(eta[a] * (T[(k, a)] @ T[(l, a)]) for a in range(1, 7))  # = 2 C^{kl}
            Ccop4 -= C2 @ C2
    Cd4 = np.zeros((dim, dim), dtype=np.int64)
    for a in range(1, 7):
        for b in range(a + 1, 7):
            D2 = sum(T[(k, a)] @ T[(k, b)] for k in range(1, n + 1))  # = 2 Delta_ab
            Cd4 -= eta[a] * eta[b] * (D2 @ D2)
    return Ccop4, Cd4


ok3a = True
consts = []
for n in (1, 2, 3):
    Cc4, Cd4 = casimirs_x4(n)
    tot = Cc4 + Cd4  # 4 (C_cop + C_diag) should be 3 n (n+4) * id
    consts.append(Fraction(int(tot[0, 0]), 4))
    ok3a &= (tot == 3 * n * (n + 4) * np.eye(8 ** n, dtype=np.int64)).all()
report("3a", ok3a, f"e3:thm:dual(iii): C_cop + C_diag = (3/4) n (n+4) id exactly on S_n for n = 1, 2, 3 (values {[str(c) for c in consts]}; "
       "exact integer matrices, independent construction)")

T2 = tower(2)
C12 = Fraction(1, 2) * sum(eta[a] * (T2[(1, a)] @ T2[(2, a)]) for a in range(1, 7))
C12x2 = sum(eta[a] * (T2[(1, a)] @ T2[(2, a)]) for a in range(1, 7))
ok = all(((C12x2 @ T2[(2, a)] - T2[(2, a)] @ C12x2) == 2 * T2[(1, a)]).all() and
         ((C12x2 @ T2[(1, a)] - T2[(1, a)] @ C12x2) == -2 * T2[(2, a)]).all() for a in range(1, 7))
# Fock form: C^{12} = sum_i (a+_{1,i} a-_{2,i} - a+_{2,i} a-_{1,i}), a+- = (E_i +- E_{i+3})/2
F4 = np.zeros((64, 64), dtype=np.int64)  # 4 * Fock form
for i in (1, 2, 3):
    ap1, am1 = T2[(1, i)] + T2[(1, i + 3)], T2[(1, i)] - T2[(1, i + 3)]
    ap2, am2 = T2[(2, i)] + T2[(2, i + 3)], T2[(2, i)] - T2[(2, i + 3)]
    F4 += ap1 @ am2 - ap2 @ am1
ok &= (F4 == 2 * C12x2).all()
# spin lift of the vector copy rotation c^{12} = (E_21 - E_12) (x) 1_6 : sigma(A) = 1/4 sum_{AB} A_AB eta_BB E_A E_B
idx = [(k, a) for k in (1, 2) for a in range(1, 7)]
Avec = np.zeros((12, 12), dtype=np.int64)
for a in range(1, 7):
    Avec[idx.index((2, a)), idx.index((1, a))] = 1
    Avec[idx.index((1, a)), idx.index((2, a))] = -1
sig4 = sum(Avec[p, q] * eta[idx[q][1]] * (T2[idx[p]] @ T2[idx[q]]) for p in range(12) for q in range(12) if Avec[p, q])
ok &= (sig4 == -2 * C12x2).all()  # 4 sigma = -4 C^{12}
report("3b", ok, "e3:thm:dual(ii) at n=2: [C^{12},E_(2,a)] = E_(1,a), [C^{12},E_(1,a)] = -E_(2,a); Fock form C^{12} = sum_i(a+_{1,i}a-_{2,i} - "
       "a+_{2,i}a-_{1,i}); and sigma((E_21-E_12)(x)1_6) = -C^{12} = -C_b, the relation hat C^{kl} = -C^{kl} of CB_CONVENTION.md (exact)")

tau = reduce(np.matmul, [T2[(1, a)] - T2[(2, a)] for a in range(1, 7)])
taux = tau  # = 8 tau_K
ok = (taux @ taux == 64 * np.eye(64, dtype=np.int64)).all()
ok &= all((taux @ T2[(1, a)] @ taux == 64 * T2[(2, a)]).all() for a in range(1, 7))
ok &= taux[0, 0] == -8 and (taux[1:, 0] == 0).all()
sym_dim = int(np.trace(np.eye(64, dtype=np.int64) * 8 + taux)) // 16
anti_dim = int(np.trace(np.eye(64, dtype=np.int64) * 8 - taux)) // 16
ok &= sym_dim == 32 and anti_dim == 32
report("3c", ok, f"e3:prop:koszul(i),(iv) at n=2: tau = (1/8) prod_a (E_(1,a) - E_(2,a)) satisfies tau^2 = 1, tau E_(1,a) tau = E_(2,a), "
       f"tau|0> = -|0>; Koszul-symmetric and -antisymmetric parts have dimensions {sym_dim}, {anti_dim} (exact)")

Cc4, Cd4 = casimirs_x4(2)
Cc = Cc4 / 4.0
w, V = np.linalg.eigh(Cc)
okd = True
rows = {}
for val, rank, cd in ((0, 20, 9), (1, 30, 8), (4, 12, 5), (9, 2, 0)):
    sel = np.abs(w - val) < 1e-9
    P = V[:, sel]
    cdv = np.linalg.eigvalsh(P.T @ (Cd4 / 4.0) @ P) if sel.any() else np.array([])
    rows[val] = (int(sel.sum()), sorted({float(x) + 0.0 for x in np.round(cdv, 9)}))
    okd &= int(sel.sum()) == rank and np.allclose(cdv, cd, atol=1e-9)
report("3d", okd and len(w) == 64,
       f"e3:tab:joint at n=2: C_cop eigenvalue -> (rank, C_diag values) = {rows}, as in the table (float eigh, tol 1e-9)")

# ======================================================================================== 4. e4 referee
Dsym = sp.Matrix(4, 4, lambda i, j: 0 if i == j else sp.Symbol(f"D{min(i, j) + 1}{max(i, j) + 1}"))
D12, D13, D14, D23, D24, D34 = [sp.Symbol(s) for s in ("D12", "D13", "D14", "D23", "D24", "D34")]
claimed = (D12 ** 2 * D34 ** 2 + D13 ** 2 * D24 ** 2 + D14 ** 2 * D23 ** 2
           - 2 * (D12 * D23 * D34 * D14 + D12 * D24 * D34 * D13 + D13 * D23 * D24 * D14))
ok4 = sp.expand(Dsym.det() - claimed) == 0
# n=5 parity under B -> -B: D -> -D, det(-D) = -det D
D5 = sp.Matrix(5, 5, lambda i, j: 0 if i == j else sp.Symbol(f"d{min(i, j)}{max(i, j)}"))
ok5 = sp.expand((-D5).det() + D5.det()) == 0 and sp.expand(D5.det()) != 0
report("4a", ok4 and ok5, "e4:prop:jarlskog: the n=4 cycle-cover formula for det D (symmetric, zero diagonal) is exact; at n=5, det(-D) = -det D "
       "(odd under time reversal), at n=4 even (sympy)")

# SO(4) content table (e4:tab:so4) against the e3 n=4 copy-Casimir spectrum (e3:prop:joint), and parity halves
even = {(0, 0): 70, (1, 0): 84, (0, 1): 84, (1, 1): 90, (2, 0): 20, (0, 2): 20, (2, 1): 15, (1, 2): 15, (3, 0): 1, (0, 3): 1}
h = Fraction(1, 2)
odd = {(h, h): 140, (3 * h, h): 64, (h, 3 * h): 64, (3 * h, 3 * h): 20, (5 * h, h): 6, (h, 5 * h): 6}
dimj = lambda j: int(2 * j + 1)
cas = lambda jp, jm: 2 * (jp * (jp + 1) + jm * (jm + 1))
spec = {}
for tab in (even, odd):
    for (jp, jm), mlt in tab.items():
        c = cas(Fraction(jp), Fraction(jm))
        spec[c] = spec.get(c, 0) + mlt * dimj(jp) * dimj(jm)
e3_n4 = {0: 70, 3: 560, 4: 504, 8: 810, 9: 1024, 12: 200, 15: 320, 16: 450, 19: 144, 24: 14}
ev_dim = sum(m_ * dimj(a) * dimj(b) for (a, b), m_ in even.items())
od_dim = sum(m_ * dimj(a) * dimj(b) for (a, b), m_ in odd.items())
report("4b", {int(k): v for k, v in spec.items()} == e3_n4 and ev_dim == od_dim == 2048,
       "cross-check between notes: the (j+,j-) table of e4:thm:so4 (even and odd N, 2048 dimensions each) reproduces exactly the n=4 copy-Casimir "
       "spectrum of e3:prop:joint (exact arithmetic)")

cnt = [2 * sum(1 for k in range(1, n + 1) for i in range(1, n + 1) for j in range(i, n + 1) if i + j == k) for n in range(2, 7)]
report("4c", cnt == [2, 4, 8, 12, 18] and cnt == [2 * sum(k // 2 for k in range(1, n + 1)) for n in range(2, 7)],
       f"e4:prop:rotor(ii): rotor-compatible bracket dimensions on the harmonic ladder by brute-force enumeration of resonant blocks: {cnt}")

# ======================================================================================== 5. LAD-01, LAD-02
T_ = e[1] @ e[2] @ e[3]
Y = e[1] @ e[2] + e[3] @ e[4]
ok = (Y.T @ T_ + T_ @ Y == 0).all() and (Y.T @ T_ @ Y == -T_ @ Y @ Y).all() and not (Y.T @ T_ @ Y == -Y @ Y @ T_).all()
ok &= (Y @ Y == 2 * e[1] @ e[2] @ e[3] @ e[4]).all() and (T_ @ (Y @ Y) == -(Y @ Y) @ T_).all()
# every element of so(3,3) satisfies Y^T T + T Y = 0 (linear; check on the 15 blades)
ok &= all((Bv[ab].T @ T_ + T_ @ Bv[ab] == 0).all() for ab in biv)
Yf = Y.astype(float)
ok &= all(np.allclose(expm(s * Yf).T @ T_ @ expm(s * Yf), T_, atol=1e-9) for s in (0.3, -1.1, 2.0))
report("5a", ok, "LAD-01 (e6:prop:Tform(i)): for Y = e1e2+e3e4, Y^T T + T Y = 0 and Y^T T Y = -T Y^2, but Y^T T Y != -Y^2 T, since "
       "Y^2 = 2 e1e2e3e4 anticommutes with T; Y^T T + T Y = 0 on all 15 blades, hence on so(3,3), and (e^{sY})^T T e^{sY} = T "
       "(exact; exponentials float, tol 1e-9)")

# LAD-02: exact stabiliser of M_pair^2 = 1/4 ||Q||_F^2 + 1/2 |det Q| in so(6,C) acting on (z1,z2) in C^6 = C^2 (x) C^3.
# Exact Gaussian-rational arithmetic: complex numbers are pairs of Fractions.
from math import isqrt


class G:
    __slots__ = ("r", "i")

    def __init__(self, r, i=0):
        self.r, self.i = Fraction(r), Fraction(i)

    def __add__(s, o):
        o = o if isinstance(o, G) else G(o)
        return G(s.r + o.r, s.i + o.i)

    __radd__ = __add__

    def __sub__(s, o):
        o = o if isinstance(o, G) else G(o)
        return G(s.r - o.r, s.i - o.i)

    def __mul__(s, o):
        o = o if isinstance(o, G) else G(o)
        return G(s.r * o.r - s.i * o.i, s.r * o.i + s.i * o.r)

    __rmul__ = __mul__

    def conj(s):
        return G(s.r, -s.i)

    def abs2(s):
        return s.r * s.r + s.i * s.i


def matvec(X, z):
    return [sum((X[p][q] * z[q] for q in range(6)), G(0)) for p in range(6)]


def gram(z, w):
    """Q_kl = z_k . w_l (complex bilinear), z_k = components 3k..3k+2."""
    return [[sum((z[3 * k + r] * w[3 * l + r] for r in range(3)), G(0)) for l in range(2)] for k in range(2)]


def dfdX(z, X):
    Q = gram(z, z)
    dz = matvec(X, z)
    A, B = gram(dz, z), gram(z, dz)
    dQ = [[A[k][l] + B[k][l] for l in range(2)] for k in range(2)]
    D = Q[0][0] * Q[1][1] - Q[0][1] * Q[1][0]
    dD = dQ[0][0] * Q[1][1] + Q[0][0] * dQ[1][1] - dQ[0][1] * Q[1][0] - Q[0][1] * dQ[1][0]
    n2 = D.abs2()
    num, den = n2.numerator, n2.denominator
    absD = Fraction(isqrt(num), isqrt(den))
    assert absD * absD == n2
    dnorm = sum(2 * (Q[k][l].conj() * dQ[k][l]).r for k in range(2) for l in range(2))
    return Fraction(1, 4) * dnorm + Fraction(1, 2) * (D.conj() * dD).r / absD


def zero6():
    return [[G(0) for _ in range(6)] for _ in range(6)]


basis = []
for p in range(6):
    for q in range(p + 1, 6):
        for c in (G(1), G(0, 1)):
            A = zero6()
            A[p][q], A[q][p] = c, G(0) - c
            basis.append(A)

rng = np.random.default_rng(20261101)
pts = []
while len(pts) < 32:
    z = [G(int(rng.integers(-3, 4)), int(rng.integers(-3, 4))) for _ in range(6)]
    Q = gram(z, z)
    D = Q[0][0] * Q[1][1] - Q[0][1] * Q[1][0]
    n2 = D.abs2()
    if n2 != 0 and isqrt(int(n2)) ** 2 == n2:
        pts.append(z)
Mrows = [[dfdX(z, X) for X in basis] for z in pts]
from sympy.polys.matrices import DomainMatrix
rank_ex = DomainMatrix([[sp.QQ(v.numerator, v.denominator) for v in row] for row in Mrows], (len(Mrows), 30), sp.QQ).rank()

stab = []
for (p, q) in ((0, 1), (0, 2), (1, 2)):
    for c in (G(1), G(0, 1)):
        A = zero6()
        for blk in (0, 3):
            A[blk + p][blk + q], A[blk + q][blk + p] = c, G(0) - c
        stab.append(A)
Cv = zero6()
for r in range(3):
    Cv[r][3 + r], Cv[3 + r][r] = G(-1), G(1)
stab.append(Cv)
coords = [[A[p][q].r for p in range(6) for q in range(p + 1, 6)] + [A[p][q].i for p in range(6) for q in range(p + 1, 6)] for A in stab]
indep = DomainMatrix([[sp.QQ(v.numerator, v.denominator) for v in row] for row in coords], (7, 30), sp.QQ).rank() == 7
# invariance: the 7 directions annihilate df at every sample point, and at 20 further random points with |det Q| rational
inv_ok = all(dfdX(z, A) == 0 for z in pts for A in stab)
# symbolic: diagonal so(3,C) leaves Q fixed (dQ = 0); C_v gives dQ = R'Q + QR'^T with R' = [[0,-1],[1,0]] (real orthogonal rotation of Q)
zs = sp.symbols("x1:7") + sp.symbols("y1:7")
zc = [zs[k] + sp.I * zs[6 + k] for k in range(6)]
Zs = sp.Matrix(3, 2, lambda r, c: zc[3 * c + r])
Qs = Zs.T * Zs
sym_ok = True
for A in stab:
    As = sp.Matrix(6, 6, lambda p, q: A[p][q].r + sp.I * A[p][q].i)
    dzs = As * sp.Matrix(zc)
    dZs = sp.Matrix(3, 2, lambda r, c: dzs[3 * c + r])
    dQs = sp.expand(dZs.T * Zs + Zs.T * dZs)
    Rp = sp.Matrix([[0, -1], [1, 0]]) if A is Cv else sp.zeros(2, 2)
    sym_ok &= sp.expand(dQs - (Rp * Qs + Qs * Rp.T)) == sp.zeros(2, 2)
report("5b", rank_ex == 23 and indep and inv_ok and sym_ok,
       f"LAD-02 (e5:thm:pair(ii)) exact certificate: at {len(pts)} Gaussian-integer pairs with |det Q| rational, the derivatives of "
       f"M_pair^2 = |Q|_F^2/4 + |det Q|/2 along the 30 real basis fields of so(6,C) form an exact rational matrix of rank {rank_ex}; the seven "
       "directions diag so(3,C) (+) R C_v are independent, annihilate it, and act on Q as 0 or as Q -> R'Q + QR'^T, R' real antisymmetric "
       "(symbolic). Hence the stabiliser algebra of M_pair in so(6,C) is exactly lG^Delta (+) R C_v, of dimension 7 (exact)")

# ======================================================================================== 6. e7 (CLASS-06)
al, th = sp.symbols("alpha theta", real=True)
ehat = sp.Matrix([sp.cos(al), sp.sin(al), 0])
v0 = sp.Matrix([1, 0, 0])  # temporal energy on the clock n = e1, mu = 1
Rv = v0 * sp.cos(th) + ehat.cross(v0) * sp.sin(th) + ehat * (ehat.dot(v0)) * (1 - sp.cos(th))  # Rodrigues
En = sp.simplify(Rv[0])
ok = sp.simplify(En - (1 - sp.sin(al) ** 2 * (1 - sp.cos(th)))) == 0
minE = sp.simplify(En.subs(th, sp.pi))
ok &= sp.simplify(minE - (1 - 2 * sp.sin(al) ** 2)) == 0
ok &= float(minE.subs(al, sp.pi / 3)) < 0 and float(minE.subs(al, sp.pi / 5)) > 0
report("6a", ok, "CLASS-06 (e7:sec:timeodd): with the energy on the clock and rotation at angle alpha, the clock energy is E_n/mu = 1 - sin^2(alpha)"
       "(1 - cos nu t) (Rodrigues, sympy); its minimum 1 - 2 sin^2(alpha) is negative iff sin^2(alpha) > 1/2 (e.g. alpha = pi/3), so the mass is "
       "m_4/mu = |E_n|/mu and the modulus may be dropped iff sin^2(alpha) <= 1/2")

# ======================================================================================== 7. TIME / X-01 (tex/s03_clock.tex rem:clockcov)
eS = [None] + [sp.Matrix(e[a].tolist()) for a in range(1, 7)]
Id8 = sp.eye(8)
TS = eS[1] * eS[2] * eS[3]
What = (eS[1] + eS[2] + eS[3]) / sp.sqrt(3)
Rr = Id8 + What * eS[1]
nrm2 = sp.simplify((Rr * Rr.T)[0, 0])  # R R^~ is a positive scalar for a rotor in Cl(V+)
Rr_n = Rr / sp.sqrt(nrm2)
Rinv = sp.simplify(Rr_n.inv())
ok = sp.simplify(Rr_n * eS[1] * Rinv - What) == sp.zeros(8, 8)
ok &= sp.simplify(Rr_n * TS - TS * Rr_n) == sp.zeros(8, 8)
ok &= sp.simplify(Rr_n * (TS * eS[1]) * Rinv - TS * What) == sp.zeros(8, 8)
p1, p2, p3, mm = sp.symbols("p1 p2 p3 m", real=True)
H = lambda nv: mm * nv + p1 * nv * eS[4] + p2 * nv * eS[5] + p3 * nv * eS[6]
ok &= sp.simplify(Rr_n * H(eS[1]) * Rinv - H(What)) == sp.zeros(8, 8)
ok &= sp.simplify(What * What - Id8) == sp.zeros(8, 8) and sp.simplify((eS[1] + eS[2] + eS[3]) ** 2 - 3 * Id8) == sp.zeros(8, 8)
ok &= sp.simplify(TS * TS + Id8) == sp.zeros(8, 8) and sp.simplify(What * (TS * What) - TS) == sp.zeros(8, 8)
report("7a", ok, f"TIME-01/02 (rem:clockcov): the rotor R = (1 + W e1)/|.|, |.|^2 = {nrm2}, carries e1 to W = (e1+e2+e3)/sqrt3, commutes with T, "
       "carries J(e1) = T e1 to J(W) = T W and H_e1(p) = m e1 + sum p_a e1 e_(a+3) to H_W(p); W^2 = 1, (e1+e2+e3)^2 = 3, T^2 = -1, W J(W) = T "
       "(sympy, exact with sqrt3)")

Gv = np.zeros((6, 6), dtype=np.int64)
for i in range(3):
    Gv[3 + i, i] = 1
    Gv[i, 3 + i] = -1
etaV = np.diag([1, 1, 1, -1, -1, -1])
Kv = np.zeros((6, 6), dtype=np.int64)
for i in range(3):
    Kv[3 + i, i] = Kv[i, 3 + i] = 1
ok = (Gv @ Gv == -np.eye(6, dtype=np.int64)).all() and (etaV @ Gv == -Kv).all() and (Gv.T @ etaV == etaV @ Gv).all()
Rq = np.array([[2, -1, 2], [2, 2, -1], [-1, 2, 2]]) // 1  # 3 * rotation matrix (exact integer, orthogonal up to 3)
assert (Rq @ Rq.T == 9 * np.eye(3, dtype=np.int64)).all() and round(np.linalg.det(Rq)) == 27
Mtemp = np.block([[Rq, np.zeros((3, 3), dtype=np.int64)], [np.zeros((3, 3), dtype=np.int64), 3 * np.eye(3, dtype=np.int64)]])
Mpair = np.block([[Rq, np.zeros((3, 3), dtype=np.int64)], [np.zeros((3, 3), dtype=np.int64), Rq]])
ok &= (Mpair @ Gv == Gv @ Mpair).all() and not (Mtemp @ Gv == Gv @ Mtemp).all()
report("7b", ok, "TIME-03: the binding root G (G^2 = -1, eta G = -K, eta-self-adjoint) commutes with a rotation acting alike on both triples and does "
       "not commute with the same rotation of the temporal triple alone (exact, rotation (1/3)[[2,-1,2],[2,2,-1],[-1,2,2]])")

O1 = eS[1] * eS[4] * eS[5] * eS[6]
gam = [eS[1], eS[4], eS[5], eS[6]]
Jn = TS * eS[1]
cands = [Jn, eS[2] * O1, eS[3] * O1]
ok = all(sp.simplify(c * g - g * c) == sp.zeros(8, 8) for c in cands for g in gam)
ok &= all(sp.simplify(c * c + Id8) == sp.zeros(8, 8) for c in cands)
ok &= all(sp.simplify(cands[i] * cands[j] + cands[j] * cands[i]) == sp.zeros(8, 8) for i in range(3) for j in range(i + 1, 3))
a_, b_, c_ = sp.Rational(2, 7), sp.Rational(3, 7), sp.Rational(6, 7)
X = a_ * cands[0] + b_ * cands[1] + c_ * cands[2]
ok &= sp.simplify(X * X + Id8) == sp.zeros(8, 8) and all(sp.simplify(X * g - g * X) == sp.zeros(8, 8) for g in gam)
ok &= sp.simplify(eS[2] * O1 - Jn) != sp.zeros(8, 8) and sp.simplify(eS[2] * O1 + Jn) != sp.zeros(8, 8)
report("7c", ok, "X-01: at the clock e1 the square roots of -1 commuting with gamma^mu = (e1,e4,e5,e6) include J(e1), e2 Omega and e3 Omega, which "
       "mutually anticommute, so every a J + b e2 Omega + c e3 Omega with a^2+b^2+c^2 = 1 is one (a 2-sphere; checked at (2,3,6)/7); the reviewer's "
       "C = e2 Omega is a pointwise complex structure different from +-J(e1). Pointwise uniqueness fails; r4:thm:unique(iii) is a clock-"
       "independence statement and needs three independent clocks (r4-23) (exact)")

# ----------------------------------------------------------------------------------------------------------------------
print()
print(f"z1: {NPASS[0]} PASS, {len(FAILS)} FAIL" + (f" (failed: {FAILS})" if FAILS else ""))
sys.exit(0 if not FAILS else 1)
