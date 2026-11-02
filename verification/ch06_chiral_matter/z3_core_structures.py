#!/usr/bin/env python3
"""z3: Phase 3 adjudication checks for the global manuscript (tex/), review items CORE-*, and one check for the literature audit
(z3-9, the author's literature record), 1 November 2026.

    python3 verification/core/z3_core.py

Prints one 'PASS z3-<id>: <claim>' or 'FAIL z3-<id>: <claim>' line per claim; exits 0 iff all pass.

Realisation (tex/s02_arena.tex, tex/sA_realisation.tex): e1=sx(x)1(x)1, e2=sz(x)sx(x)1, e3=sz(x)sz(x)sx, e4=ep(x)1(x)1, e5=sz(x)ep(x)1,
e6=sz(x)sz(x)ep, ep=[[0,-1],[1,0]], signature (+,+,+,-,-,-). Exact integer or sympy arithmetic unless a line says 'float'.
Each PASS line states the mathematical assertion tested; it does not by itself validate the surrounding prose.
"""
import itertools
import sys
from functools import reduce

import numpy as np
import sympy as sp

NPASS, FAILS = [0], []


def report(cid, ok, text):
    if ok:
        NPASS[0] += 1
        print(f"PASS z3-{cid}: {text}")
    else:
        FAILS.append(cid)
        print(f"FAIL z3-{cid}: {text}")
    sys.stdout.flush()


sx = np.array([[0, 1], [1, 0]], dtype=np.int64)
sz = np.array([[1, 0], [0, -1]], dtype=np.int64)
ep = np.array([[0, -1], [1, 0]], dtype=np.int64)
I2 = np.eye(2, dtype=np.int64)
kr = lambda *a: reduce(np.kron, a)
e = [None, kr(sx, I2, I2), kr(sz, sx, I2), kr(sz, sz, sx), kr(ep, I2, I2), kr(sz, ep, I2), kr(sz, sz, ep)]
I8 = np.eye(8, dtype=np.int64)
bl = lambda S: reduce(np.matmul, [e[i] for i in S], I8)
eta = [None, 1, 1, 1, -1, -1, -1]
for a in range(1, 7):
    assert (e[a] @ e[a] == eta[a] * I8).all()

# ======================================================================================== CORE-05 (Example 2.2)
B23, B14, B12, B56 = bl((2, 3)), bl((1, 4)), bl((1, 2)), bl((5, 6))
ok = (B23 @ B14 == B14 @ B23).all() and (B23 @ B56 == B56 @ B23).all() and (B14 @ B12 == -B12 @ B14).all()
ok &= (B23 @ B23 == -I8).all() and (B14 @ B14 == I8).all() and (B23.T == -B23).all() and (B14.T == B14).all()
report("1", ok, "CORE-05 (ex:blades): e2e3 (rotation, square -1) and e1e4 (boost, square +1) share no index and commute; e2e3 and e5e6 commute; "
       "e1e4 and e1e2 share one index and anticommute (exact)")

# ======================================================================================== CORE-02 (lem:tracepot, thm:twophases)
rots = [bl((4, 5)), bl((4, 6)), bl((5, 6))]
boosts = [bl((1, 4)), bl((1, 5)), bl((1, 6))]
Kbl = []
for r in range(7):
    for S in itertools.combinations(range(1, 7), r):
        B = bl(S)
        if (B.T == -B).all() and all((B @ R == R @ B).all() for R in rots) and all((B @ X == -X @ B).all() for X in boosts):
            Kbl.append(S)
ok = Kbl == [(1, 2), (1, 3), (1, 2, 3), (2, 4, 5, 6), (3, 4, 5, 6), (2, 3, 4, 5, 6)]
eS = [None] + [sp.Matrix(m.tolist()) for m in e[1:]]
blS = lambda S: reduce(lambda a, b: a * b, [eS[i] for i in S], sp.eye(8))
c = sp.symbols("c1:7", real=True)
M = sum((ci * blS(S) for ci, S in zip(c, Kbl)), sp.zeros(8))
A = sp.expand(M.T * M)
N = sp.expand(A.trace())
D = sp.expand(A - N / 8 * sp.eye(8))
D2 = sp.expand(D * D)
u = sp.Matrix([c[0], c[1], c[5]])
v = sp.Matrix([c[3], c[4], c[2]])
ok &= sp.expand(D.trace()) == 0 and sp.expand(D2 - D2[0, 0] * sp.eye(8)) == sp.zeros(8)
ok &= sp.expand(D2[0, 0] - 4 * (u.cross(v)).dot(u.cross(v))) == 0
ok &= sp.expand(N - 8 * (u.dot(u) + v.dot(v))) == 0
Q = sp.expand((A * A).trace())
ok &= sp.expand(Q - (N ** 2 / 8 + 8 * D2[0, 0])) == 0
report("2", ok, "CORE-02 (lem:tracepot): the kernel space K of the interior so(1,3) is spanned by e1e2, e1e3, e1e2e3, e2e4e5e6, e3e4e5e6, "
       "e2e3e4e5e6; M^T M = (N/8) 1 + D with D traceless, D^2 = d^2 1, d^2 = 4|u x v|^2, N = 8(|u|^2+|v|^2) for u = (c1,c2,c6), "
       "v = (c4,c5,c3); Q = N^2/8 + 8 d^2 (sympy, exact). Hence q = Q/N^2 in [1/8, 1/4], q = 1/8 iff u || v, q = 1/4 iff u perp v, |u| = |v|")

t = sp.symbols("t", real=True)
Mt = t * blS((1, 2, 3))
Nt = sp.expand((Mt.T * Mt).trace())
Qt = sp.expand(((Mt.T * Mt) ** 2).trace())
Vt = sp.expand(-1 * Nt + (-1) * Nt ** 2 + 1 * Qt)
ok = sp.expand(Vt - (-8 * t ** 2 - 56 * t ** 4)) == 0
# radial minimum formula
mu2, al, be, qq, NN = sp.symbols("mu2 alpha beta q N", positive=True)
Vrad = -mu2 * NN + (al + be * qq) * NN ** 2
Nst = sp.solve(sp.diff(Vrad, NN), NN)[0]
ok &= sp.simplify(Nst - mu2 / (2 * (al + be * qq))) == 0 and sp.simplify(Vrad.subs(NN, Nst) + mu2 ** 2 / (4 * (al + be * qq))) == 0
report("3", ok, "CORE-02: at mu^2 = 1, alpha = -1, beta = 1, on the ray M = t gamma^0 J = t e1e2e3, V = -8t^2 - 56t^4 (unbounded below); for "
       "a(q) = alpha + beta q > 0 the radial minimum is N = mu^2/(2a(q)) with value -mu^4/(4a(q)) (sympy)")

# orbit of g_I = <J, Omega, e2 Omega, e3 Omega> through the Gamma=+1 Majorana line and through gamma^0 J: tangent dimension 3,
# equal to the dimension of the minimum set at fixed N
J = bl((2, 3))
Om = bl((1, 4, 5, 6))
gI = [J, Om, e[2] @ Om, e[3] @ Om]
Gam = bl((1, 2, 3, 4, 5, 6))
Kmats = [bl(S).astype(float) for S in Kbl]


def coords(Mf):
    return np.array([np.trace(B.T @ Mf) / 8.0 for B in Kmats])


def orbit_dim(Mf):
    tang = [coords(Mf @ X - X @ Mf) for X in gI]
    return np.linalg.matrix_rank(np.array(tang), tol=1e-9)


T = bl((1, 2, 3)).astype(float)
# a kernel supported on Gamma = +1: M = (1+Gamma)/2 * K * (1+Gamma)/2 restricted to K; find one null z with M^T M prop. (1+Gamma)/2
# u perp v, |u| = |v|: u = (1,0,0), v = (0,1,0) -> c1 = 1, c4 = 1
cand = []
for (uu, vv) in (((1, 0, 0), (0, 1, 0)), ((1, 0, 0), (0, -1, 0)), ((0, 1, 0), (1, 0, 0)), ((0, 1, 0), (-1, 0, 0))):
    cc = [uu[0], uu[1], vv[2], vv[0], vv[1], uu[2]]
    Mf = sum(ci * B for ci, B in zip(cc, Kmats))
    AA = Mf.T @ Mf
    cand.append((cc, np.allclose(AA, 2 * (np.eye(8) + Gam), atol=1e-12), np.allclose(AA, 2 * (np.eye(8) - Gam), atol=1e-12)))
plus = [cc for cc, p, m in cand if p]
ok = len(plus) >= 1
Mplus = sum(ci * B for ci, B in zip(plus[0], Kmats))
ok &= orbit_dim(T) == 3 and orbit_dim(Mplus) == 3
# the orbit leaves the Gamma=+1 support: the tangent direction [M, e2 Omega] is not supported on Gamma=+1
Xt = Mplus @ (e[2] @ Om) - (e[2] @ Om) @ Mplus
ok &= not np.allclose(Xt, 0.25 * (np.eye(8) + Gam) @ Xt @ (np.eye(8) + Gam), atol=1e-9)
report("4", ok, "CORE-02 (thm:twophases): the g_I-orbit through gamma^0 J and through a Gamma=+1 Majorana kernel (u = e_1, v = e_2 in the "
       "coordinates of z3-2) has tangent dimension 3, equal to the dimension of the minimum set at fixed N; the direction [M, e2 Omega] "
       "leaves the Gamma=+1 support (float, exact entries)")

# ======================================================================================== CORE-01 (vector-mass normalisation)
# n = 3 realisation and the conventions of tex/s08_pstar.tex, s19_kernels.tex (as used in the photon-selection checks)
I512 = np.eye(512, dtype=np.int64)
E3 = [None]
om = reduce(np.matmul, e[1:])
for kk in range(3):
    for a in range(1, 7):
        E3.append(kr(*([om] * kk + [e[a]] + [I8] * (2 - kk))))
bl3 = lambda S: reduce(np.matmul, [E3[i] for i in S], I512)
interior = [1, 4, 5, 6]
internal = [i for i in range(1, 19) if i not in interior]
Widx, Cidx = [7, 13, 14, 15], [10, 11, 12, 16, 17, 18]
J3 = bl3(internal)
Gam3 = J3 @ bl3(interior)
chiW, chiC, qA = bl3(Widx), bl3(Cidx), E3[2] @ E3[3]
X_, Z_, U_ = Gam3 @ chiW, -(J3 @ chiC), -(J3 @ qA)
Pst = ((I512 + X_) @ (I512 + Z_) @ (I512 + U_)).astype(float) / 8
w = [E3[i] for i in Widx]
tau = [-(w[1] @ w[2]) / 2, -(w[2] @ w[0]) / 2, -(w[0] @ w[1]) / 2]
tL = [((I512 - chiW) @ t) / 2 for t in tau]
tR = [((I512 + chiW) @ t) / 2 for t in tau]
HC = -(E3[10] @ E3[16] + E3[11] @ E3[17] + E3[12] @ E3[18])
colour = [(E3[i] @ E3[j]).astype(float) for i, j in itertools.combinations(Cidx, 2)]
PS = [t.astype(float) for t in tL + tR] + colour
Jf, Gf = J3.astype(float), Gam3.astype(float)
TR3h = -Jf @ tR[2]
BmL = -Jf @ HC.astype(float) / 3
Id = np.eye(512)
# P_nu = block x (Gamma=+1) x (hat T^R_3 = +1/2) x (B-L = -1); spectral projectors on the block
Pnu = Pst @ ((Id + Gf) / 2) @ (2 * TR3h @ (TR3h + Id / 2)) @ ((Id / 3 - BmL) * 0.75)
Pnu = Pst @ Pnu @ Pst
ok = np.allclose(Pnu @ Pnu, Pnu, atol=1e-12) and abs(np.trace(Pnu) - 4) < 1e-9
Phi0 = Pnu @ Pst @ (E3[1] @ bl3([3, 9, 13, 14, 16, 17, 18])).astype(float) @ Pst @ Pnu
ok &= np.allclose(Phi0.T @ Phi0, Pnu, atol=1e-12) and np.allclose(Phi0.T, -Phi0, atol=1e-12) and np.allclose(Jf @ Phi0, -Phi0 @ Jf, atol=1e-12)
G = np.array([[np.trace(a.T @ b) / 512 for b in PS] for a in PS])
Mm = np.array([[np.trace((a @ Phi0 - Phi0 @ a).T @ (b @ Phi0 - Phi0 @ b)) for b in PS] for a in PS]) / np.trace(Phi0.T @ Phi0)
evs = np.sort(np.linalg.eigvals(np.linalg.solve(G, Mm)).real)
vals = sorted({round(float(x), 9) + 0.0 for x in evs})
mult = {vv: int(np.sum(np.abs(evs - vv) < 1e-7)) for vv in vals}
ok &= mult == {0.0: 12, 4.0: 8, 20.0: 1}
ok &= np.allclose(np.diag(G)[:6], 1 / 8) and np.allclose(np.diag(G)[6:], 1)
report("5", ok, f"CORE-01 (thm:compat(b)): with Phi_0 of eq:Phi0 (Phi_0^T Phi_0 = P_nu, rank 4, J Phi_0 = -Phi_0 J) and the 21 Pati-Salam generators "
       f"(module metric G = 1/8 on t^L, t^R and 1 on the colour bivectors), M_ab = <[t_a,Phi0],[t_b,Phi0]>/|Phi0|^2 relative to G has spectrum "
       f"{mult} (exact entries; float eigenvalues)")
# canonical masses from eq:Saction: for module-orthonormal generators and canonical fields with coupling g_A, the mass term is
# (1/2) g_A^2 <[A,Phi],[A,Phi]> = (1/2) g_A^2 |Phi|^2 M, |Phi|^2 = 2 v^2  ->  m^2 = 2 v^2 g_A^2 * eig(G^-1 M) = 8 g_A^2 v^2 and 40 g_A^2 v^2
m2_over_gA2v2 = sorted({round(2 * x, 9) + 0.0 for x in evs})
ok = m2_over_gA2v2 == [0.0, 8.0, 40.0]
# Example 25.3 at (v_R, kappa) -> (20, 0): v^2 = 2 v_R^2 = 800; with g_A = 1 the X and Z' masses^2 are 6400 and 32000
ok &= abs(2 * 800 * 4 - 6400) < 1e-9 and abs(2 * 800 * 20 - 32000) < 1e-9
# conventional normalisation: a colour generator with charge eigenvalues +-1/2 on two colours (lambda_3/2 of su(4)) has module norm 1/8, as t^R_i
Xc = (E3[10] @ E3[16] - E3[11] @ E3[17]).astype(float) / 4
chg = np.round(np.linalg.eigvalsh(Pst @ (-(Jf @ Xc)) @ Pst), 9)
ok &= abs(np.trace(Xc.T @ Xc) / 512 - 1 / 8) < 1e-12 and set(np.unique(chg)) <= {-0.5, 0.0, 0.5} and (chg == 0.5).sum() > 0
ok &= abs(np.trace(tR[2].T @ tR[2]) / 512 - 1 / 8) < 1e-12
report("6", ok, f"CORE-01: generators of conventional normalisation (Tr_fund T T = 1/2: t^R_i, and the colour generator (E10E16-E11E17)/4 with "
       f"charges +-1/2) have module norm 1/8; eq:Saction gives squared masses {m2_over_gA2v2} x g_A^2 v^2 for module-orthonormal canonical gauge fields of "
       "coupling g_A; thm:compat(b) states g^2 v^2 and 5 g^2 v^2, so its g is g = sqrt(8) g_A; ex:EWspectrum's 6400 and 32032.1 at v_R = 20 "
       "(v^2 = 800) are the canonical values at g_A = 1, i.e. g^2 = 8")

# ======================================================================================== CORE-03 (closure uniqueness, n >= 2)
etaV = np.diag([1, 1, 1, -1, -1, -1])
Gv = np.zeros((6, 6), dtype=np.int64)
for i in range(3):
    Gv[3 + i, i], Gv[i, 3 + i] = 1, -1
so33 = []
for a_ in range(6):
    for b_ in range(a_ + 1, 6):
        Xm = np.zeros((6, 6), dtype=np.int64)
        Xm[a_, b_] = 1
        Xm[b_, a_] = -etaV[a_, a_] * etaV[b_, b_]  # X^T eta + eta X = 0
        so33.append(Xm)
assert all((Xm.T @ etaV + etaV @ Xm == 0).all() for Xm in so33)
# l_G = Cent(G) in so(3,3): exact null space
from sympy import Matrix as SM
Acomm = SM(np.array([(Gv @ Xm - Xm @ Gv).ravel() for Xm in so33], dtype=object).T.tolist())
nsG = Acomm.nullspace()
lG = [sum((int(cvec[i]) * so33[i] for i in range(15)), np.zeros((6, 6), dtype=np.int64)) for cvec in [v_ * sp.ilcm(*[x.q for x in v_]) for v_ in nsG]]
ok = len(lG) == 6
# commutant of the diagonal l_G in End(R^2 (x) R^6), exact
n_ = 2
diag = [np.kron(np.eye(n_, dtype=np.int64), L) for L in lG]
dimE = (6 * n_) ** 2
rows = []
for L in diag:
    # [L, Y] = 0 as a linear map on vec(Y): (I (x) L - L^T (x) I) vec(Y)
    rows.append(np.kron(np.eye(6 * n_, dtype=np.int64), L) - np.kron(L.T, np.eye(6 * n_, dtype=np.int64)))
Mbig = np.vstack(rows)
from sympy.polys.matrices import DomainMatrix
rk = DomainMatrix([[sp.ZZ(int(x)) for x in r] for r in Mbig], Mbig.shape, sp.ZZ).rank()
dim_comm = dimE - rk
cross = np.kron(np.array([[0, 1], [0, 0]], dtype=np.int64), np.eye(6, dtype=np.int64))  # maps constituent 2 into constituent 1
ok &= dim_comm == 2 * n_ ** 2 and all((cross @ L == L @ cross).all() for L in diag)
ok &= all((np.kron(np.eye(n_, dtype=np.int64), Gv) @ L == L @ np.kron(np.eye(n_, dtype=np.int64), Gv)).all() for L in diag)
report("7", ok, f"CORE-03: l_G = Cent(G) in so(3,3) has dimension {len(lG)}; under the diagonal action on V_2 = R^2 (x) R^6 its commutant in End(V_2) has "
       f"dimension {dim_comm} = 2n^2, i.e. {{A (x) 1 + B (x) G}} = M_2(C), and contains the cross-constituent map E_12 (x) 1_6: equivariant maps "
       "between distinct constituents do not vanish under a common Lorentz action (exact ranks)")

# ======================================================================================== CORE-06 (complex structure on the kernel space)
def idx_par(Sa, Sb):  # blades e_Sa, e_Sb commute iff |Sa||Sb| - |Sa n Sb| is even
    return (len(Sa) * len(Sb) - len(set(Sa) & set(Sb))) % 2 == 0


all18 = set(range(1, 19))
Xset = sorted(all18 - set(Widx))                     # X = Gamma chi_W
Zset = sorted(set(internal) - set(Cidx))             # Z = -J chi_C
Uset = sorted(set(internal) - {2, 3})                # U = -J q_A
Jset = internal
rngk = np.random.default_rng(7)
found = {"dirac": None, "majorana": None}
tries = 0
while (found["dirac"] is None or found["majorana"] is None) and tries < 20000:
    tries += 1
    I_ = sorted(rngk.choice(internal, size=int(rngk.integers(1, 14)), replace=False).tolist())
    S_ = [1] + I_
    if not all(idx_par(S_, T_) for T_ in (Xset, Zset, Uset)):
        continue
    B = bl3(S_)
    if not (B.T == -B).all():
        continue
    kind = "dirac" if idx_par(S_, Jset) else "majorana"
    if found[kind] is None:
        found[kind] = S_
ok = found["dirac"] is not None and found["majorana"] is not None
Mdir = bl3(found["dirac"]).astype(float) @ Pst
Mmaj = bl3(found["majorana"]).astype(float) @ Pst
rot = (E3[4] @ E3[5]).astype(float)
boo = (E3[1] @ E3[4]).astype(float)


def is_kernel(Mk):
    return (np.allclose(Mk.T, -Mk, atol=1e-12) and np.allclose(Mk @ rot, rot @ Mk, atol=1e-12) and np.allclose(Mk @ boo, -boo @ Mk, atol=1e-12)
            and np.allclose(Pst @ Mk @ Pst, Mk, atol=1e-12) and np.linalg.norm(Mk) > 1e-6)


ok &= is_kernel(Mdir) and is_kernel(Mmaj)
ok &= np.allclose(Jf @ Mdir, Mdir @ Jf) and np.allclose(Gf @ Mdir, -Mdir @ Gf) and np.allclose(Jf @ Mmaj, -Mmaj @ Jf) and np.allclose(Gf @ Mmaj, Mmaj @ Gf)
JMd = Jf @ Mdir
ok &= np.allclose(JMd.T, JMd) and not is_kernel(JMd)          # J times a Dirac kernel is symmetric: leaves K
Cd = Jf @ Mdir @ Gf
Cm = Jf @ Mmaj
ok &= is_kernel(Cd) and is_kernel(Cm)
ok &= np.allclose(Jf @ Cd @ Gf, -Mdir) and np.allclose(Jf @ Cm, -Mmaj)  # C^2 = -1
tt = PS[3]  # t^R_1: C commutes with the Pati-Salam action M -> [t, M]
ok &= np.allclose(Jf @ (tt @ Mdir - Mdir @ tt) @ Gf, tt @ Cd - Cd @ tt) and np.allclose(Jf @ (tt @ Mmaj - Mmaj @ tt), tt @ Cm - Cm @ tt)
# realified dimensions of the complex decomposition, and the 10 + 126 bookkeeping
cdims = {"(10,3,1)": 30, "(10,1,3)": 30, "(15,2,2)": 60, "2(6,1,1)": 12, "(1,2,2)": 4}
ok &= sum(cdims.values()) == 136 and 2 * sum(cdims.values()) == 272 and 6 + 30 + 30 + 60 == 126 and 6 + 4 == 10
report("8", ok, f"CORE-06 (thm:kernelspace): explicit kernels gamma^0 e_I P_* of Dirac type (I = {found['dirac'][1:]}) and Majorana type "
       f"(I = {found['majorana'][1:]}); J M is symmetric for a Dirac kernel, so M -> JM does not preserve K; C(M) = JM on Majorana and JM Gamma "
       "on Dirac kernels preserves K, squares to -1 and commutes with the Pati-Salam action; complex dimensions 30+30+60+12+4 = 136 = 272/2, "
       "= 10 + 126 (float, exact entries)")

# ======================================================================================== literature audit (Krasnov 2018, the author's literature record)
# Krasnov's SO(7,7) model has the commutant SO(4,6) of the Lorentz group, with compact part SO(4) x SO(6). In the n = 3 realisation the
# directions W (4) and C (6) carry opposite signs; the question is which bivectors of so(W + C) preserve the chiral block P_*.
E3f = {i: E3[i].astype(float) for i in Widx + Cidx}           # signed permutation matrices: float products are exact
sq = {i: int(round((E3f[i] @ E3f[i])[0, 0])) for i in Widx + Cidx}
ok = all(np.array_equal(E3f[i] @ E3f[i], sq[i] * Id) for i in Widx + Cidx)
sig_W, sig_C = {sq[i] for i in Widx}, {sq[i] for i in Cidx}
ok &= len(sig_W) == 1 and len(sig_C) == 1 and sig_W != sig_C
same = [(i, j) for S in (Widx, Cidx) for i, j in itertools.combinations(S, 2)]
mixed = [(i, j) for i in Widx for j in Cidx]
ok &= len(same) == 21 and len(mixed) == 24
ok &= all(idx_par([i, j], T_) for i, j in same for T_ in (Xset, Zset, Uset))               # commute with X, Z, U
ok &= all(not idx_par([i, j], Xset) and not idx_par([i, j], Zset) and idx_par([i, j], Uset) for i, j in mixed)
biv = {(i, j): E3f[i] @ E3f[j] for i, j in same + mixed}
ok &= all(np.array_equal(biv[k] @ biv[k], -Id) for k in same) and all(np.array_equal(biv[k] @ biv[k], Id) for k in mixed)
leak = max(np.abs(Pst @ biv[k] @ Pst).max() for k in mixed)
keep = max(np.abs(biv[k] @ Pst - Pst @ biv[k]).max() for k in same)
ok &= leak < 1e-12 and keep < 1e-12
report("9", ok, f"literature (Krasnov 2018): the generators of W = {Widx} square to {sig_W.pop():+d} and those of C = {Cidx} to {sig_C.pop():+d}, "
       "so so(W + C) = so(4,6); its 21 same-block bivectors (so(4) + so(6), each squaring to -1, compact) commute with X, Z, U and with P_*; "
       "its 24 mixed W-C bivectors (squaring to +1, the non-compact part) anticommute with X and Z, and P_* B P_* = 0 for each: P_* is "
       "preserved by the compact Pati-Salam part only (blade rule exact; projections float, tol 1e-12)")

print()
print(f"z3: {NPASS[0]} PASS, {len(FAILS)} FAIL" + (f" (failed: {FAILS})" if FAILS else ""))
sys.exit(0 if not FAILS else 1)
