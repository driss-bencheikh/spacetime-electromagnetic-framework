#!/usr/bin/env python3
"""c1: the selector, the light cone, type selection and kinetic signatures of Chapter 12 (Phase 9C).

    python3 verification/ch12_causality/c1_selector_types.py

thm:selector, the pairing omega_S, thm:cone, thm:selector3, prop:types, thm:degree, thm:flat, prop:kineticB, prop:positivealt. Exact arithmetic
(sympy) in the 8x8 realisation of Appendix A; the three-constituent statements use the 512x512 signed-permutation realisation with integer
entries (exact), and two statements are tested at random spinors (marked numerical).
"""
import itertools
import os
import sys
from functools import reduce

import numpy as np
import sympy as sp

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "common"))
from clifford8 import I8, bl, comm, e, eta, vec, vec_coords  # noqa: E402

NPASS, FAILS = [0], []


def report(cid, ok, text):
    if ok:
        NPASS[0] += 1
        print(f"PASS c1-{cid}: {text}")
    else:
        FAILS.append(cid)
        print(f"FAIL c1-{cid}: {text}")
    sys.stdout.flush()


Z8 = sp.zeros(8, 8)
BIVS = [(a, b) for a, b in itertools.combinations(range(1, 7), 2)]
biv = [e[a] * e[b] for a, b in BIVS]


def wedge(a, b):
    A, Bm = vec(a), vec(b)
    return (A * Bm - Bm * A) / 2


def gram_det(a, b):
    ip = lambda p, q: sum(eta[i + 1] * p[i] * q[i] for i in range(6))
    return ip(a, a) * ip(b, b) - ip(a, b) ** 2


def centraliser_dim(C):
    M = sp.Matrix.hstack(*[sp.Matrix(list(comm(C, X))) for X in biv])
    return 15 - M.rank()


# ---------------------------------------------------------------- thm:selector
av = sp.symbols("a1:7"); bv = sp.symbols("b1:7")
C = wedge(av, bv)
ok = (C * C + gram_det(av, bv) * I8).applyfunc(sp.expand) == Z8
cases = {"positive <e2,e3>": ([0, 1, 0, 0, 0, 0], [0, 0, 1, 0, 0, 0]), "negative <e5,e6>": ([0, 0, 0, 0, 1, 0], [0, 0, 0, 0, 0, 1]),
         "indefinite <e1,e4>": ([1, 0, 0, 0, 0, 0], [0, 0, 0, 1, 0, 0]), "tilted positive": ([0, 2, 0, 1, 0, 0], [0, 0, 1, 0, 0, 0]),
         "tilted indefinite": ([1, 1, 0, 0, 0, 0], [0, 0, 0, 1, 0, 1]), "degenerate <e1+e4,e2>": ([1, 0, 0, 1, 0, 0], [0, 1, 0, 0, 0, 0])}
res = {}
for name, (a, b) in cases.items():
    Cn = wedge(a, b)
    gd = gram_det(a, b)
    Cu = Cn / sp.sqrt(abs(gd)) if gd != 0 else Cn
    sq = (Cu * Cu)
    typ = "-1" if sq == -I8 else ("+1" if sq == I8 else ("0" if sq == Z8 else "?"))
    cd = centraliser_dim(Cn)
    antisym = (Cu.T == -Cu)
    res[name] = (typ, cd, antisym)
ok &= res["positive <e2,e3>"] == ("-1", 7, True) and res["negative <e5,e6>"] == ("-1", 7, True) and res["indefinite <e1,e4>"][:2] == ("+1", 7)
ok &= res["tilted positive"][:2] == ("-1", 7) and res["tilted positive"][2] is False and res["tilted indefinite"][:2] == ("+1", 7)
ok &= res["degenerate <e1+e4,e2>"][0] == "0"
ok &= e[2] * e[3] * e[2] * e[3] == -I8
report("1", ok, f"thm:selector: (a^b)^2 = -det Gram(a,b) for general a, b (exact identity); normalised planes: {res} (type of C^2, centraliser dimension "
       "in so(3,3), C antisymmetric on S): C^2 = -1 exactly for definite planes (Lorentzian complement), +1 for indefinite ones ((2,2) complement), 0 "
       "for a degenerate plane; centralisers of nondegenerate planes have dimension 7, so each type of choice is 15 - 7 = 8 dimensional; coordinate "
       "planes with C^2 = -1 give antisymmetric (orthogonal) complex structures, the tilted positive plane does not; the surplus plane gives J = e2e3")

# ---------------------------------------------------------------- the pairing omega_S and thm:cone
wS = e[4] * e[5] * e[6]
ok = all(wS * e[a] * wS.inv() == -e[a].T for a in range(1, 7))
ok &= all((X.T * wS + wS * X) == Z8 for X in biv)
Vol = bl(1, 2, 3, 4, 5, 6)
ok &= Vol * Vol == I8
Splus = (I8 + Vol).columnspace(); Sminus = (I8 - Vol).columnspace()
Sp = sp.Matrix.hstack(*Splus); Sm = sp.Matrix.hstack(*Sminus)
ok &= Sp.T * wS * Sp == sp.zeros(4, 4)
mA = [wS * e[a] for a in range(1, 7)]
ok &= all(m.T == -m for m in mA)
ok &= all(X.T * mA[a - 1] + mA[a - 1] * X == -wS * comm(X, e[a]) for X in biv for a in range(1, 7))
restr = sp.Matrix.hstack(*[sp.Matrix(list(Sp.T * m * Sp)) for m in mA])
ok &= restr.rank() == 6
J = e[2] * e[3]
ok &= (Sp.T * J * J * Sp) == -(Sp.T * Sp) and all(comm(J, e[a] * e[b]) == Z8 for a, b in [(1, 4), (1, 5), (1, 6), (4, 5), (5, 6), (4, 6)])
report("2", ok, "pairing and thm:cone (i)-(ii): omega_S = e4e5e6 satisfies omega_S e_A omega_S^-1 = -e_A^T, so psi^T omega_S phi is Spin(3,3)-invariant "
       "(X^T omega_S + omega_S X = 0 for all 15 bivectors) and vanishes on S+ x S+; m_A = omega_S e_A are antisymmetric with X^T m_A + m_A X = -omega_S[X, e_A], "
       "and their restrictions to S+ span the 6-dimensional Lambda^2(S+)*; on S+, J^2 = -1 and the interior so(1,3) commutes with J (exact)")

psi = sp.Matrix(sp.symbols("p1:5"))
ps = Sp * psi


def vC(Cm, s=ps):
    return [sp.expand((s.T * wS * e[a] * Cm * s)[0]) for a in range(1, 7)]


def eta_sq(v):
    return sp.expand(sum(eta[i + 1] * v[i] ** 2 for i in range(6)))


ok = True
info = {}
for name, (a, b) in cases.items():
    Cn = wedge(a, b)
    gd = gram_det(a, b)
    v = vC(Cn)
    null = eta_sq(v) == 0
    Nmat = sp.Matrix([a, b])
    in_perp = all(sp.expand(sum(Nmat[k, i] * v[i] for i in range(6))) == 0 for k in range(2))   # v_A has a lower index: N^A v_A = 0
    jac = sp.Matrix([[sp.diff(vi, pj) for pj in psi] for vi in v]).subs({psi[0]: 1, psi[1]: 2, psi[2]: -1, psi[3]: 3})
    info[name] = (null, in_perp, jac.rank())
    ok &= null and in_perp and jac.rank() == 3
# (iv) definiteness: C = J (complement <e1,e4,e5,e6>, time e1) vs C = e1e4 (complement <e2,e3,e5,e6>)
qJ = sp.hessian(vC(J)[0], psi) / 2
evJ = qJ.eigenvals()
definite = all(sp.sign(k_) == sp.sign(list(evJ)[0]) and k_ != 0 for k_ in evJ)
Cb = e[1] * e[4]
vb = vC(Cb)
indef = all(len({sp.sign(k_) for k_ in sp.hessian(vb[i], psi).eigenvals() if k_ != 0}) == 2 for i in (1, 2, 4, 5))
ok &= definite and indef
# (v) v_{-C} = -v_C, other nappe on S-, v_J = -j with j_A = psi^T gamma^0 e_A psi
jA = [sp.expand((ps.T * e[1] * e[a] * ps)[0]) for a in range(1, 7)]
ok &= all(sp.expand(x + y) == 0 for x, y in zip(vC(J), jA))
psm = Sm * psi
vm = vC(J, psm)[0]
sgn_p = sp.sign(list(evJ)[0]); sgn_m = sp.sign(list((sp.hessian(vm, psi) / 2).eigenvals())[0])
ok &= sgn_p == -sgn_m and all(sp.expand(x + y) == 0 for x, y in zip(vC(-J), vC(J)))
report("3", ok, f"thm:cone (iii)-(v): for every plane tested, v_C(psi) = (psi^T omega_S e_A C psi)_A on S+ is in N-perp, eta(v_C,v_C) = 0 identically and the image "
       f"has dimension 3 (Jacobian rank at a generic point): {info}; for C = J the e1-component is a definite quadratic form on S+ (one nappe), for "
       "C = e1e4 every component of the (2,2) complement is indefinite; v_J = -j with j_A = psi^T gamma^0 e_A psi, the J-twisted square on S- lies in the "
       "other nappe and v_{-C} = -v_C (exact). That the image fills the whole nappe is the Klein correspondence (inherited), not recomputed")

# ---------------------------------------------------------------- thm:selector3 (n = 3)
sx = np.array([[0, 1], [1, 0]]); sz = np.array([[1, 0], [0, -1]]); ep = np.array([[0, -1], [1, 0]]); I2 = np.eye(2, dtype=int)
kr = lambda *a: reduce(np.kron, a)
e8 = [None, kr(sx, I2, I2), kr(sz, sx, I2), kr(sz, sz, sx), kr(ep, I2, I2), kr(sz, ep, I2), kr(sz, sz, ep)]
om = reduce(np.matmul, e8[1:])
E = [None] + [kr(*([om] * k + [e8[a]] + [np.eye(8, dtype=int)] * (2 - k))).astype(float) for k in range(3) for a in range(1, 7)]
Id = np.eye(512)
blN = lambda S: reduce(np.matmul, [E[i] for i in S], Id)
sq = [None] + [int(round(np.trace(E[i] @ E[i]) / 512)) for i in range(1, 19)]
interior = [1, 4, 5, 6]; W3 = [i for i in range(1, 19) if i not in interior]
Jn, Om = blN(W3), blN(interior)
Gam = Jn @ Om
ok = np.allclose(Jn @ Jn, -Id) and np.allclose(Om @ Om, -Id) and np.allclose(Jn @ Om, Om @ Jn) and np.allclose(Gam @ Gam, Id)
Pg = (Id + Gam) / 2
ok &= np.allclose(Om @ Pg, -Jn @ Pg)


def blade_sq(S):
    k = len(S)
    return (-1) ** (k * (k - 1) // 2) * int(np.prod([sq[i] for i in S]))


cnt_ok = 0
for Ep in itertools.combinations(range(1, 19), 4):
    comp = [i for i in range(1, 19) if i not in Ep]
    p = int(np.prod([sq[i] for i in Ep]))
    if blade_sq(list(Ep)) == p and blade_sq(comp) == p:
        cnt_ok += 1
neg_in = lambda Ep: sum(1 for i in Ep if sq[i] < 0)
lor = sum(1 for Ep in itertools.combinations(range(1, 19), 4) if neg_in(Ep) in (1, 3))
ok &= cnt_ok == 3060
nb = sum(1 for i, j in itertools.combinations(range(1, 19), 2) if ((i in W3) == (j in W3)))
ok &= nb == 97 and 153 - nb == 56
neg = [i for i in range(1, 19) if sq[i] < 0]; posg = [i for i in range(1, 19) if sq[i] > 0]
bodd, beven = blN(neg), blN(posg)
ok &= np.allclose(bodd.T, -bodd) and np.allclose(beven.T, beven)
inv_ok = all(np.allclose((E[i] @ E[j]).T @ b + b @ (E[i] @ E[j]), 0) for b in (bodd, beven) for i, j in itertools.combinations(range(1, 19), 2))
sym_tw = {nm: [a for a in range(1, 19) if np.allclose((b @ E[a] @ Jn).T, b @ E[a] @ Jn)] for nm, b in (("odd", bodd), ("even", beven))}
ok &= inv_ok and sym_tw["odd"] == W3 and sym_tw["even"] == W3
ok &= all(np.allclose((bodd @ E[a]).T, bodd @ E[a]) for a in range(1, 19))
rng = np.random.default_rng(9)
nonnull = 0
w_, V_ = np.linalg.eigh(Pg)
Sp3 = V_[:, w_ > 0.5]                       # an orthonormal basis of S+ = {Gamma = +1}, 256-dimensional
ok &= Sp3.shape[1] == 256
for _ in range(5):
    s = Sp3 @ rng.normal(size=256)          # random spinors of S+ (v1.0.1: the sample is restricted to S+)
    vA = np.array([s @ bodd @ E[a] @ s for a in range(1, 19)])
    nonnull += abs(sum(sq[a] * vA[a - 1] ** 2 for a in range(1, 19))) > 1e-6
ok &= nonnull == 5
report("4", ok, f"thm:selector3: at n = 3, J^2 = Omega^2 = -1, [J,Omega] = 0, Gamma = J Omega with Gamma^2 = 1 and Omega = -J on Gamma = +1; for all 3060 coordinate "
       f"four-planes J_E'^2 = Omega_E'^2 = prod_{{i in E'}} e_i^2 (so the complement volume is a complex structure exactly for the {lor} Lorentzian ones); the "
       "centraliser of J among the 153 bivectors has dimension 97 = 6 + 91 (56-dimensional choices); beta_odd (antisymmetric) and beta_even (symmetric) are "
       "Spin(9,9)-invariant, beta e_A J is symmetric exactly for the 14 directions of W_3, and beta_odd e_A is symmetric for all 18 (exact); at 5 random spinors "
       "of S+ (seed 9) the untwisted square is not null (numerical)")

# Dirac current at n = 3: future causal (numerical sample) and the matrix argument
g0 = E[1]; gs = [E[4], E[5], E[6]]
Ms = [g0 @ g for g in gs]
ok = all(np.allclose(M, M.T) for M in Ms) and all(np.allclose(Ms[i] @ Ms[j] + Ms[j] @ Ms[i], 2 * Id * (i == j)) for i in range(3) for j in range(3))
caus = 0
for _ in range(20):
    s = rng.normal(size=512)
    j0 = s @ s; jv = np.array([s @ M @ s for M in Ms])
    caus += j0 >= np.linalg.norm(jv) - 1e-9
ok &= caus == 20
# product states psi = phi1 (x) phi2 (x) phi3: null when the first factor is chiral, timelike for |000> + |100> in the first factor
om8 = om.astype(float)
chi1 = np.kron(np.kron(om8, np.eye(8)), np.eye(8))      # chirality of the first constituent (volume element of its Cl(3,3))
cur = lambda p: (p @ p, np.array([p @ M @ p for M in Ms]))  # noqa: E731
b = np.eye(8)
psi_c = np.kron(np.kron(b[0], b[0]), b[0])
psi_n = np.kron(np.kron(b[0] + b[4], b[0]), b[0])
j0c, jvc = cur(psi_c)
j0n, jvn = cur(psi_n)
chiral_first = abs(abs(psi_c @ chi1 @ psi_c) - psi_c @ psi_c) < 1e-12 and abs(psi_n @ chi1 @ psi_n) < 1e-12
ok &= chiral_first and abs(j0c ** 2 - jvc @ jvc) < 1e-12 and abs(j0n - 2) < 1e-12 and np.allclose(jvn, 0)
w8, V8 = np.linalg.eigh(om8)
nulls = 0
for _ in range(5):
    f1 = V8[:, w8 > 0] @ rng.normal(size=int(np.sum(w8 > 0)))
    p = np.kron(np.kron(f1, rng.normal(size=8)), rng.normal(size=8))
    j0, jv = cur(p)
    nulls += abs(j0 ** 2 - jv @ jv) < 1e-8 * j0 ** 2
ok &= nulls == 5
report("5", ok, "n = 3 Dirac current j_mu = psi^T gamma^0 gamma_mu psi: the three matrices gamma^0 gamma^k are symmetric and anticommute with square 1, so "
       "j^0 = |psi|^2 >= |j| for every psi (future causal; exact argument from the matrix identities, checked on 512x512 matrices) and 20 random spinors "
       "confirm it (numerical). It uses gamma^0 and the transpose, i.e. the positive inner product, as the text says. On product states the current is null "
       "when the first factor is chiral (|000>, and 5 random states with a chiral first factor), and (2,0,0,0), timelike, for |000> + |100> in the first factor")

# ---------------------------------------------------------------- prop:types, thm:degree on S+ = sl(4,R)
def restrict(X):
    return (Sp.T * Sp).inv() * Sp.T * X * Sp


ok = sp.Matrix.hstack(*[sp.Matrix(list(restrict(X))) for X in biv]).rank() == 15 and all(restrict(X).trace() == 0 for X in biv)
y = sp.symbols("y")


def invs(Xr):
    return tuple(sp.nsimplify((Xr ** kx).trace()) for kx in (2, 3, 4))


tp, tn, tb = invs(restrict(e[2] * e[3])), invs(restrict(e[5] * e[6])), invs(restrict(e[1] * e[4]))
ok &= tp == (-4, 0, 4) and tn == (-4, 0, 4) and tb == (4, 0, 4)
ok &= sp.factor(restrict(e[2] * e[3]).charpoly(y).as_expr()) == (y ** 2 + 1) ** 2 and sp.factor(restrict(e[1] * e[4]).charpoly(y).as_expr()) == (y - 1) ** 2 * (y + 1) ** 2


def orient(Xr):
    v1 = sp.Matrix([1, 0, 0, 0]); v2 = sp.Matrix([0, 1, 0, 0])
    if sp.Matrix.hstack(v1, Xr * v1, v2, Xr * v2).det() == 0:
        v2 = sp.Matrix([0, 0, 1, 0])
    return sp.sign(sp.Matrix.hstack(v1, Xr * v1, v2, Xr * v2).det())


ok &= orient(restrict(e[2] * e[3])) == -orient(restrict(e[5] * e[6]))
report("6", ok, f"prop:types: so(3,3) restricted to S+ is sl(4,R) (rank 15, traceless); rotations of a positive and of a negative plane both have "
       f"(t2,t3,t4) = {tp} = {tn} and characteristic polynomial (y^2+1)^2, boosts {tb} and (y-1)^2(y+1)^2; the two complex structures induce opposite "
       "orientations of S+ (sign of det(v, Cv, w, Cw)), so they lie on distinct SL(4,R)-orbits with identical invariants (exact)")

T = lambda X, kx: (X ** kx).trace()
A0 = sp.diag(1, -1, 0, 0); A0[2, 3] = -1; A0[3, 2] = 1
B0 = sp.diag(sp.Matrix([[1, -1], [1, 1]]), sp.Matrix([[-1, -1], [1, -1]]))
C0 = sp.diag(-1, -3, 0, 0); C0[2, 2] = 2; C0[2, 3] = -3; C0[3, 2] = 3; C0[3, 3] = 2
D0 = sp.diag(1, -1, 2, -2)
R12 = sp.diag(sp.Matrix([[0, -1], [1, 0]]), sp.Matrix([[0, -2], [2, 0]]))
Jc = restrict(e[2] * e[3]); Kp = restrict(e[1] * e[4])
vals = {"A0": [T(A0, 2), T(A0, 3), T(A0, 4)], "B0": [T(B0, 2), T(B0, 3), T(B0, 4)], "C0": [T(C0, 2), T(C0, 3)],
        "J": [T(Jc, 2), T(Jc, 4)], "K": [T(Kp, 2), T(Kp, 4)], "D0": [T(D0, 2), T(D0, 3), T(D0, 4)], "R12": [T(R12, 2), T(R12, 3), T(R12, 4)]}
ok = vals == {"A0": [0, 0, 4], "B0": [0, 0, -16], "C0": [0, -120], "J": [-4, 4], "K": [4, 4], "D0": [10, 0, 34], "R12": [-10, 0, 34]}
ok &= all(X.trace() == 0 for X in (A0, B0, C0, D0, R12))
eE, fF, gG = sp.symbols("e f g")
top = lambda t2, t4: eE * t2 ** 3 + gG * t2 * t4
sol = sp.solve([top(-4, 4) + top(4, 4), top(10, 34) + top(-10, 34)], [eE, gG], dict=True)
ok &= top(-4, 4) == -64 * eE - 16 * gG and top(4, 4) == 64 * eE + 16 * gG and top(10, 34) == 1000 * eE + 340 * gG
ok &= sp.solve([64 * eE + 16 * gG, 1000 * eE + 340 * gG], [eE, gG]) == {eE: 0, gG: 0}
Bw = sp.Rational(3, 5) * e[2] * e[3] + sp.Rational(4, 5) * e[5] * e[6]
tw = invs(restrict(Bw))
ok &= tw == (-4, 0, sp.Rational(4804, 625))
grade4 = (Bw * Bw - (Bw * Bw).trace() / 8 * I8)
ok &= grade4 != Z8
report("7", ok, f"thm:degree: the test elements of the proof are traceless and have the stated invariants {vals}; non-negativity of the degree-6 top part "
       "e t2^3 + f t3^2 + g t2 t4 at the complex and para-complex structures forces 64e + 16g = 0 and at D0 and the two-speed rotation forces 1000e + 340g = 0, "
       f"hence e = g = 0; B_w = (3/5)e2e3 + (4/5)e5e6 has (t2,t3,t4) = {tw}, so it shares t2 = -4, t3 = 0 with J but is not simple (B_w^2 has a grade-4 "
       "part) (exact). The bounded-below argument itself is the proof in the text; the script checks its numerical inputs")

bs = sp.symbols("q0:15")
Xg = sum((bs[i] * restrict(biv[i]) for i in range(15)), sp.zeros(4, 4))
f = (T(Jc + Xg, 2) + 4) ** 2 + T(Jc + Xg, 3) ** 2 + (T(Jc + Xg, 4) - 4) ** 2
Hf = sp.hessian(f, bs).subs({b: 0 for b in bs})
ok = Hf.rank() == 1 and f.subs({b: 0 for b in bs}) == 0
report("8", ok, "thm:degree (iv): F = (t2+4)^2 + t3^2 + (t4-4)^2 vanishes at J and its Hessian at J on the 15-dimensional so(3,3) has rank one (exact)")

# ---------------------------------------------------------------- thm:flat
nn = (e[1] + e[4]) * e[5]
ok = comm(J, nn) == Z8 and nn * nn == Z8
cd = centraliser_dim(J + nn)
interior_biv = [e[1] * e[4], e[1] * e[5], e[1] * e[6], e[4] * e[5], e[5] * e[6], e[4] * e[6]]
cent = sp.Matrix.hstack(*[sp.Matrix(list(comm(J + nn, X))) for X in biv]).nullspace()
cent_ops = [sum((v[i] * biv[i] for i in range(15)), Z8) for v in cent]
span = sp.Matrix.hstack(*[sp.Matrix(list(X)) for X in cent_ops])
contains = all(span.rank() == sp.Matrix.hstack(span, sp.Matrix(list(X))).rank() for X in interior_biv)
Kb = e[1] * e[4] / 2
ok &= cd == 3 and not contains and comm(Kb, nn) == -nn and comm(Kb, J) == Z8
tt = sp.symbols("t", real=True)
eK = sp.cosh(tt / 2) * I8 + sp.sinh(tt / 2) * e[1] * e[4]
eKi = sp.cosh(tt / 2) * I8 - sp.sinh(tt / 2) * e[1] * e[4]
ok &= (eK * (J + nn) * eKi - (J + sp.exp(-tt) * nn)).applyfunc(lambda q: sp.simplify(q.rewrite(sp.exp))) == Z8
# (ii) Stab(J)-invariant symmetric forms on so(1,3)
stab = [X for X in biv if comm(X, J) == Z8]
L = interior_biv


def coords_in(X, basis):
    Mb = sp.Matrix.hstack(*[sp.Matrix(list(Y)) for Y in basis])
    return (Mb.T * Mb).inv() * Mb.T * sp.Matrix(list(X))


adm = [sp.Matrix.hstack(*[coords_in(comm(S_, Y), L) for Y in L]) for S_ in stab]
qs = sp.symbols("w0:21")
Bf = sp.zeros(6, 6); kq = 0
for i in range(6):
    for j in range(i, 6):
        Bf[i, j] = Bf[j, i] = qs[kq]; kq += 1
eqs = []
for Am in adm:
    eqs += list(Am.T * Bf + Bf * Am)
sol = sp.linsolve([q for q in eqs if q != 0], qs)
solv = list(sol)[0]
nfree = len(set().union(*[sp.sympify(s_).free_symbols for s_ in solv]))
Bgen = Bf.subs(dict(zip(qs, solv)))
fs = sorted(set().union(*[sp.sympify(s_).free_symbols for s_ in solv]), key=str)
sigs = []
for vals_ in ([1, 0], [0, 1], [2, 3], [-1, 5]):
    Bn = Bgen.subs(dict(zip(fs, vals_)))
    evs = Bn.eigenvals()
    sigs.append((sum(m for v_, m in evs.items() if v_ > 0), sum(m for v_, m in evs.items() if v_ < 0)))
ok &= len(stab) == 7 and nfree == 2 and all(s_ == (3, 3) for s_ in sigs)
report("9", ok, f"thm:flat: n = (e1+e4)e5 commutes with J, n^2 = 0; the centraliser of J + n in so(3,3) has dimension {cd} and does not contain so(1,3); "
       "K = e1e4/2 has [K,n] = -n, [K,J] = 0 and e^{tK}(J+n)e^{-tK} = J + e^{-t} n, so J lies in the closure of the orbit of J + n; the Stab(J)-invariant "
       f"symmetric forms on so(1,3) form a {nfree}-dimensional space and nonzero ones have signature {sigs[0]} (four samples), so a positive semidefinite "
       "invariant Hessian vanishes there (exact)")

# ---------------------------------------------------------------- prop:kineticB and prop:positivealt
Kf = sp.Matrix(15, 15, lambda i, j: (biv[i] * biv[j]).trace())
evK = Kf.eigenvals()
sigK = (sum(m for v_, m in evK.items() if v_ > 0), sum(m for v_, m in evK.items() if v_ < 0))
gold = [X for X in biv if comm(X, J) != Z8]
rotations = lambda Xs: sum(1 for X in Xs if X.T == -X)
gt = [comm(X, J) for X in gold]
Kg = sp.Matrix(len(gt), len(gt), lambda i, j: (gt[i] * gt[j]).trace())
evg = Kg.eigenvals()
sigG = (sum(m for v_, m in evg.items() if v_ > 0), sum(m for v_, m in evg.items() if v_ < 0))
ok = sigK == (9, 6) and len(gold) == 8 and rotations(gold) == 2 and rotations(interior_biv) == 3 and J.T == -J and sigG == (6, 2)
report("10", ok, f"prop:kineticB: the trace form (a multiple of the Killing form) on the 15 bivectors has signature {sigK} (boosts one sign, rotations the "
       f"other); at J the radial direction J is a rotation, the 8 Goldstone directions are 2 rotations and 6 boosts, the so(1,3) directions 3 of each, and the "
       f"induced form on the tangent space of the 8-dimensional space of choices has signature {sigG} (exact). Whatever the overall sign, at least six "
       "directions have the wrong sign; ghost freedom of an interacting theory is not addressed by this count")

ip = sp.Matrix(15, 15, lambda i, j: (biv[i] * biv[j].T).trace() / 8)
ok = ip == sp.eye(15)
inv_alg = [v for v in sp.Matrix.hstack(*[sp.Matrix([ (ip * sp.Matrix([coords_in(comm(X, Y), biv)[k] for k in range(15)]))[i] for i in range(15)])
                                          for Y in biv]).T.nullspace()] if False else None
compact = [i for i in range(15) if biv[i].T == -biv[i]]
ad = lambda X: sp.Matrix.hstack(*[coords_in(comm(X, Y), biv) for Y in biv])
anti = [i for i in range(15) if ad(biv[i]).T == -ad(biv[i])]
ok &= anti == compact and len(compact) == 6
T3 = [biv[BIVS.index(p)] for p in [(1, 2), (1, 3), (2, 3)]]
S3 = [biv[BIVS.index(p)] for p in [(4, 5), (4, 6), (5, 6)]]
xs = sp.symbols("x0:15")
Bx = sum((xs[i] * biv[i] for i in range(15)), Z8)
cc = coords_in(Bx, biv)
iT = [BIVS.index(p) for p in [(1, 2), (1, 3), (2, 3)]]; iS = [BIVS.index(p) for p in [(4, 5), (4, 6), (5, 6)]]
iP = [i for i in range(15) if i not in iT + iS]
Vp = (sum(cc[i] ** 2 for i in iT) - 1) ** 2 + sum(cc[i] ** 2 for i in iS) + sum(cc[i] ** 2 for i in iP)
pt = {xs[i]: (1 if BIVS[i] == (2, 3) else 0) for i in range(15)}
Hv = sp.hessian(Vp, xs).subs(pt)
evH = Hv.eigenvals()
ok &= evH == {8: 1, 0: 2, 2: 12} and Vp.subs(pt) == 0
ok &= e[1] * (e[1] * e[2] * e[3]) == e[2] * e[3]
ok &= all(ad(biv[BIVS.index(p)]).T != -ad(biv[BIVS.index(p)]) for p in [(1, 4), (1, 5), (1, 6)])
report("11", ok, "prop:positivealt: <X,Y> = tr(XY^T)/8 is the identity on the 15 bivectors and its invariance algebra is the 6 compact rotations so3_T + so3_S; "
       "V_+ = (|B_T|^2 - 1)^2 + |B_S|^2 + |B_p|^2 vanishes at J = e2e3 with Hessian eigenvalues 8, 0 (x2), 2 (x12); J = gamma^0 omega_T with gamma^0 = e1; the "
       "boosts e1e4, e1e5, e1e6 are not in the invariance algebra (exact)")

# ---------------------------------------------------------------- thm:selector (iii): the degenerate stratum (v1.0.1)
a_null = vec([1, 0, 0, 1, 0, 0]); b_null = vec([0, 1, 0, 0, 1, 0])
C0 = a_null * b_null
ok = (C0 * C0 == Z8) and (a_null * b_null + b_null * a_null == Z8)
cs = sp.symbols("x0:15")
X = sum((c * bv for c, bv in zip(cs, biv)), Z8)
eqs = [x for x in comm(X, C0)]
sol = sp.linsolve(eqs, cs)
(vec_sol,) = list(sol)
free = set().union(*[sp.sympify(v).free_symbols for v in vec_sol])
cent_null = len(free)
Xj = sum((c * bv for c, bv in zip(cs, biv)), Z8)
(vj,) = list(sp.linsolve([x for x in comm(Xj, e[2] * e[3])], cs))
cent_J = len(set().union(*[sp.sympify(v).free_symbols for v in vj]))
ok &= cent_null == 9 and cent_J == 7
report("12", ok, f"thm:selector (iii), degenerate stratum: for the totally null plane <e1+e4, e2+e5>, C0 = (e1+e4)(e2+e5) has C0^2 = 0 and its centraliser "
       f"in so(3,3) has dimension {cent_null}, against {cent_J} for the nondegenerate C = J (exact)")

print(f"{NPASS[0]}/{NPASS[0] + len(FAILS)} checks passed")
sys.exit(0 if not FAILS else 1)
