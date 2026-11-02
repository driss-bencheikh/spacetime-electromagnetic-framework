"""Phase 9D, batch 5: adversarial audit of selected Chapter 7-9 statements that had neither a cited check nor a written proof.

Exact combinatorics of Clifford blades (commutation by index parity), exact integer 512x512 arithmetic, exact rational arithmetic for the
stated numerical examples, the 8x8 realisation for the one-constituent photon statements, and float64 linear algebra on the exact operators
of verification/vacuum/kernel_space.py for the string statement (labelled numerical).
Each check prints PASS/FAIL a5-<n>; the script exits 0 iff all pass.
"""
import itertools
import os
import sys
from functools import reduce

import numpy as np
import sympy as sp

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(HERE, "..", "common"))

NPASS, FAILS = [0], []


def report(cid, ok, text):
    if ok:
        NPASS[0] += 1
        print(f"PASS a5-{cid}: {text}")
    else:
        FAILS.append(cid)
        print(f"FAIL a5-{cid}: {text}")
    sys.stdout.flush()


# ---------------------------------------------------------------- a5-1: thm:fourclasses (exact blade combinatorics in Cl(9,9))
ALL = frozenset(range(1, 19))
interior = frozenset([1, 4, 5, 6])
internal = ALL - interior
W, C, qA = frozenset([7, 13, 14, 15]), frozenset([10, 11, 12, 16, 17, 18]), frozenset([2, 3])
sdiff = lambda a, b: a ^ b  # noqa: E731  (blade product up to sign)
Jb = internal
Gb = ALL  # Gamma = J e1e4e5e6 is the volume blade up to sign
Xb, Zb, Ub = sdiff(Gb, W), sdiff(Jb, C), sdiff(Jb, qA)
rots = [frozenset([4, 5]), frozenset([5, 6]), frozenset([4, 6])]
boosts = [frozenset([1, a]) for a in (4, 5, 6)]


def commutes(S, T):
    return (len(S) * len(T) - len(S & T)) % 2 == 0


count, classes = 0, {}
for r in range(19):
    for S in itertools.combinations(range(1, 19), r):
        S = frozenset(S)
        if all(commutes(S, T) for T in (Xb, Zb, Ub)) and all(commutes(S, R) for R in rots):
            count += 1
            key = (commutes(S, Jb), commutes(S, Gb), all(commutes(S, B) for B in boosts), all(not commutes(S, B) for B in boosts))
            classes[key] = classes.get(key, 0) + 1
named = {"internal": classes.get((True, True, True, False), 0), "P": classes.get((True, False, False, True), 0),
         "C": classes.get((False, False, True, False), 0), "CP": classes.get((False, True, False, True), 0)}
ok = count == 8192 and named == {"internal": 2048, "P": 2048, "C": 2048, "CP": 2048} and sum(classes.values()) == 8192 and len(classes) == 4
report("1", ok, f"thm:fourclasses: the blades of Cl(9,9) commuting with X, Z, U and the three rotations number {count}, in exactly four classes "
       f"{named} by (J-linearity, Gamma, boosts) (exact, by the parity rule for blade commutation over all 2^18 blades)")

# ---------------------------------------------------------------- the 512 module (integer)
sx = np.array([[0, 1], [1, 0]]); sz = np.array([[1, 0], [0, -1]]); ep = np.array([[0, -1], [1, 0]]); I2 = np.eye(2, dtype=np.int64)
kr = lambda *a: reduce(np.kron, a)  # noqa: E731
e8 = [None, kr(sx, I2, I2), kr(sz, sx, I2), kr(sz, sz, sx), kr(ep, I2, I2), kr(sz, ep, I2), kr(sz, sz, ep)]
I8 = np.eye(8, dtype=np.int64)
om = reduce(np.matmul, e8[1:])
E = [None] + [kr(*([om] * k + [e8[a]] + [I8] * (2 - k))).astype(np.int64) for k in range(3) for a in range(1, 7)]
I512 = np.eye(512, dtype=np.int64)
Z512 = np.zeros((512, 512), dtype=np.int64)
blq = lambda S: reduce(np.matmul, [E[i] for i in sorted(S)], I512)  # noqa: E731
Jn = blq(internal)
Gam = Jn @ blq([1, 4, 5, 6])
Xs, Zs, Us = Gam @ blq(W), -(Jn @ blq(C)), -(Jn @ (E[2] @ E[3]))
Ps8 = (I512 + Xs) @ (I512 + Zs) @ (I512 + Us)
g0 = E[1]
IB = list(itertools.combinations(sorted(internal), 2))
ok = all(np.array_equal(g0 @ (E[a] @ E[b]), (E[a] @ E[b]) @ g0) for a, b in IB) and np.array_equal(g0 @ Jn, Jn @ g0)
ok &= all(np.array_equal(g0 @ (E[a] @ E[b]), (E[a] @ E[b]) @ g0) for a, b in [(4, 5), (5, 6), (4, 6)])
ok &= np.array_equal(g0 @ Gam, -Gam @ g0) and all(np.array_equal(g0 @ (E[1] @ E[a]), -(E[1] @ E[a]) @ g0) for a in (4, 5, 6))
Pm8 = g0 @ Ps8 @ g0.T
ok &= np.array_equal(Pm8 @ Ps8, Z512) and np.trace(Pm8) == 512 and np.array_equal(g0.T @ g0, I512)
report("2", ok, "prop:mirror: gamma^0 commutes with the 91 internal bivectors, J and the rotations, anticommutes with Gamma and the boosts, and maps V_* "
       "isometrically onto gamma^0 V_*, of rank 64 and orthogonal to V_* (exact integer arithmetic)")

PS = [E[a] @ E[b] for a, b in itertools.combinations(sorted(W), 2)] + [E[a] @ E[b] for a, b in itertools.combinations(sorted(C), 2)]
Kmu = [Ps8 @ g0 @ g @ Ps8 for g in (E[1], E[4], E[5], E[6])]
JG = Jn @ Gam
ok = len(PS) == 21 and np.array_equal(Jn.T, -Jn) and np.array_equal(JG.T, -JG)
ok &= all(np.array_equal(Jn @ k, k @ Jn) and np.array_equal(JG @ k, k @ JG) for k in Kmu)
ok &= all(np.array_equal(Jn @ t, t @ Jn) and np.array_equal(JG @ t, t @ JG) for t in PS)
ok &= np.array_equal(Jn @ Ps8, Ps8 @ Jn) and np.array_equal(Gam @ Ps8, Ps8 @ Gam)
report("3", ok, "prop:fermionnumber, prop:global (algebraic part): J and J Gamma are antisymmetric and commute with the kinetic kernels P_* gamma^0 gamma^mu "
       "P_* and the 21 Pati-Salam generators (exact integer arithmetic); the arithmetic F = 3B + L, 4B = F + (B - L) is immediate")

# ---------------------------------------------------------------- a5-4: ex:sumrule, thm:tracepot, ex:EWspectrum ratio (exact rationals)
lam = [sp.Rational(x) for x in ("1", "1.7", "2.3", "3.1", "2.9")]
D = [lam[k] - lam[0] for k in range(1, 5)]
ok = D == [sp.Rational(7, 10), sp.Rational(13, 10), sp.Rational(21, 10), sp.Rational(19, 10)]
mee = sp.Rational(4, 3) * D[0]
mue = D[3]
mde = sp.Rational(2, 3) * D[0] + D[3]
muu = sp.Rational(4, 3) * D[1]
mud = sp.Rational(2, 3) * D[1] + D[3]
ok &= mee == sp.Rational(14, 15) and mue == sp.Rational(19, 10) and mde == sp.Rational(71, 30) and mde == mue + mee / 2
ok &= muu == sp.Rational(26, 15) and mud == sp.Rational(83, 30) and mud == mue + muu / 2
al, be = sp.symbols("alpha beta")
a = [1, -sp.Rational(1, 2), -sp.Rational(1, 2), sp.Rational(1, 4), 0]
lk = [16 * al + 4 * be * x for x in a]
ok &= [sp.expand(lk[k] - lk[0]) for k in range(1, 5)] == [sp.expand(-4 * be * x) for x in (sp.Rational(3, 2), sp.Rational(3, 2), sp.Rational(3, 4), 1)]
ok &= abs(79.360 / 127.872 - 0.62062) < 5e-6
report("4", ok, "ex:sumrule, thm:tracepot (coefficient map), ex:EWspectrum (ratio): Delta = (0.7,1.3,2.1,1.9); m_ee^2 = 14/15, m_de^2 = 71/30 = 1.9 + "
       "m_ee^2/2, m_uu^2 = 26/15, m_ud^2 = 83/30 = 1.9 + m_uu^2/2; lambda_k = 16 alpha + 4 beta a_k gives Delta = -4 beta (3/2,3/2,3/4,1) (exact); "
       "79.360/127.872 = 0.62062 (the spectrum itself is not recomputed here)")

# ---------------------------------------------------------------- a5-5: thm:strings (numerical, on the operators of kernel_space)
import kernel_space as ks  # noqa: E402
from scipy.linalg import expm  # noqa: E402

Pst = ks.Pst
Hs = 3 * ks.Jf @ ks.BmL
Phi0 = ks.Phi0_512
w_, V_ = np.linalg.eigh(Pst)
Bk = V_[:, w_ > 0.5]
BmLb = Bk.T @ ks.BmL @ Bk
evs = np.round(np.linalg.eigvalsh((BmLb + BmLb.T) / 2), 9)
ok = set(np.unique(evs)) <= {1 / 3, -1.0} or set(np.round(np.unique(evs), 6)) == {round(1 / 3, 6), -1.0}
quark = Bk @ np.linalg.eigh((BmLb + BmLb.T) / 2)[1][:, evs > 0]
lepton = Bk @ np.linalg.eigh((BmLb + BmLb.T) / 2)[1][:, evs < 0]
th = 0.37
U = expm(th * Hs)
ok &= np.allclose(U @ quark, (np.cos(th) * np.eye(512) + np.sin(th) * ks.Jf) @ quark, atol=1e-10)
ok &= np.allclose(U @ lepton, (np.cos(3 * th) * np.eye(512) - np.sin(3 * th) * ks.Jf) @ lepton, atol=1e-10)
gs = expm(np.pi * Hs / 3)
ok &= np.allclose(gs @ Phi0 @ np.linalg.inv(gs), Phi0, atol=1e-10) and np.allclose(gs @ lepton, -lepton, atol=1e-10)
path = U @ Phi0 @ np.linalg.inv(U)
ok &= np.allclose(path, np.cos(6 * th) * Phi0 - np.sin(6 * th) * ks.Jf @ Phi0, atol=1e-10)
ok &= np.linalg.norm(Phi0) > 1e-6
report("5", ok, "thm:strings (numerical, float64 on the exact operators, 1e-10): on the block H = 3J(B-L) acts on quarks as cos t + J sin t and on leptons "
       "as cos 3t - J sin 3t; g_* = e^{pi H/3} fixes Phi_0 and is -1 on every lepton; e^{tH} Phi_0 e^{-tH} = cos 6t Phi_0 - sin 6t J Phi_0")

# ---------------------------------------------------------------- a5-6: Chapter 8 at one constituent (exact, 8x8)
from clifford8 import I8 as I8s, bl, comm, e  # noqa: E402

Z8 = sp.zeros(8, 8)
JN = (e[2] * e[3] + e[3] * e[1] + e[1] * e[2]) / sp.sqrt(3)
Q = JN / 2
Wh = (e[1] + e[2] + e[3]) / sp.sqrt(3)
fp = [Wh, e[4], e[5], e[6]]  # the four-plane of the clock What, with which J_N commutes (lem:JN)
lint = [fp[a] * fp[b] for a, b in itertools.combinations(range(4), 2)]
xs = sp.symbols("c0:64")
X = sp.Matrix(8, 8, xs)
A, _ = sp.linear_eq_to_matrix([x for L in lint + [Q] for x in (X * L - L * X)], xs)
sols = [sp.Matrix(8, 8, list(v)) for v in A.nullspace()]
Vv = bl(1, 2, 3, 4, 5, 6)
sym = [s for s in sols if s.T == s]
symspace = sp.Matrix([list((s + s.T) / 2) for s in sols])
ok = len(sols) == 4 and symspace.rank() == 2 and sp.Matrix([list(I8s), list(Vv)] + [list((s + s.T) / 2) for s in sols]).rank() == 2
ok &= sp.simplify(JN * JN + I8s) == Z8
# thm:massless1: generic rank-4 projector commuting with Q: span of complex 2-planes of C^4 (2Q complex structure)
rng = np.random.default_rng(3)
Qn = np.array(Q.evalf(), dtype=float)
Jc = 2 * Qn
# complex structure: build a basis of R^8 as C^4 with i = Jc; choose a random complex 2-plane
v = rng.standard_normal((8, 2))
basis = np.hstack([v, Jc @ v])
Qr, _ = np.linalg.qr(basis)
Pi0 = Qr @ Qr.T
ok &= np.allclose(Pi0 @ Qn, Qn @ Pi0) and abs(np.trace(Pi0) - 4) < 1e-12
so33 = [np.array((e[a] * e[b]).evalf(), dtype=float) / 2 for a, b in itertools.combinations(range(1, 7), 2)]
tang = np.array([(B @ Pi0 - Pi0 @ B).reshape(-1) for B in so33])
r14 = np.linalg.matrix_rank(tang, tol=1e-9)
ch = [B @ Pi0 - Pi0 @ B for B in so33]
# subspace of the tangent span with [Q, pi] != 0: dimension of the span minus the dimension of its Q-commuting part
Tm = tang.T
coef_null = np.linalg.svd(np.array([(Qn @ t.reshape(8, 8) - t.reshape(8, 8) @ Qn).reshape(-1) for t in tang]).T, compute_uv=False)
r_charged = np.linalg.matrix_rank(np.array([(Qn @ t.reshape(8, 8) - t.reshape(8, 8) @ Qn).reshape(-1) for t in tang]), tol=1e-9)
compact_idx = [i for i, (a, b) in enumerate(itertools.combinations(range(1, 7), 2)) if (a <= 3) == (b <= 3)]
tc = tang[compact_idx]
r_c = np.linalg.matrix_rank(tc, tol=1e-9)
r_cc = np.linalg.matrix_rank(np.array([(Qn @ t.reshape(8, 8) - t.reshape(8, 8) @ Qn).reshape(-1) for t in tc]), tol=1e-9)
ok &= r14 == 14 and r_charged == 8 and r_c == 5 and r_cc == 2
report("6", ok, "prop:rigid (exact; interior action of <What,e4,e5,e6>, the four-plane of J_N) and thm:massless1 (numerical, a random member of Gr_C(2,4)): the commutant of {l_int, Q = J_N/2} has dimension 4 with "
       "symmetric part span{1, V}; J_N^2 = -1; at a generic rank-4 projector commuting with Q the span of [B, Pi_0] over so(3,3) has dimension "
       f"{r14} and the image of its charged part under ad Q has dimension {r_charged}; compact-generated directions {r_c}, charged among them {r_cc}")

# ---------------------------------------------------------------- a5-7: bf:qym(i) (exact identity for su(4))
Tg = []
for i in range(4):
    for j in range(i + 1, 4):
        Mx = sp.zeros(4); Mx[i, j] = Mx[j, i] = sp.Rational(1, 2); Tg.append(Mx)
        Mx = sp.zeros(4); Mx[i, j] = -sp.I / 2; Mx[j, i] = sp.I / 2; Tg.append(Mx)
for k in range(1, 4):
    d = [1] * k + [-k] + [0] * (3 - k)
    Tg.append(sp.diag(*d) / sp.sqrt(2 * k * (k + 1)))
ok = len(Tg) == 15 and all(sp.simplify((Tg[a] * Tg[b]).trace() - (sp.Rational(1, 2) if a == b else 0)) == 0 for a in range(15) for b in range(15))
Zm = sp.Matrix(4, 4, sp.symbols("z0:16"))
lhs = sp.simplify(sum((T * Zm * T for T in Tg), sp.zeros(4)))
ok &= sp.simplify(lhs - (Zm.trace() * sp.eye(4) - Zm / 4) / 2) == sp.zeros(4)
ok &= sp.simplify(sum((T * T for T in Tg), sp.zeros(4)) - sp.Rational(15, 8) * sp.eye(4)) == sp.zeros(4)
report("7", ok, "bf:qym(i): for su(4) with tr(T^A T^B) = delta/2, sum_A T^A Z T^A = (tr Z - Z/4)/2 for every 4x4 Z, and C_F = 15/8 (exact)")

print(f"{NPASS[0]}/{NPASS[0] + len(FAILS)} checks passed")
sys.exit(0 if not FAILS else 1)
