"""Release audit (v1.0.0): exact certificates for statements of Chapter 7 that had no recomputation.

    python3 verification/ch07_mass_vacuum/a6_kernels_cp_certificates.py

a6-1  thm:subkernels      kernel spaces M_W, internal algebras g_W, decomposition of M_W under g_W and invariant kernels, for every Lorentz
                          algebra so(W) of the theorem (one constituent, 8x8, exact rational arithmetic on blades).
a6-2  thm:vacCP           CP Phi_0 CP^T = Phi_0; C Phi_0 C^T is the nu_L nu_L kernel and P Phi_0 P^T the e_L e_L kernel, both in (10,3,1);
                          exact integer arithmetic on the 512-dimensional module with Phi_0 rebuilt exactly (no floating point).
a6-3  prop:coupled3       the algebraic identities behind H(p)^2 = |p|^2 + y^2 r^2 P_nu on the block (exact integer arithmetic).
a6-4  thm:Bcounts         Pati-Salam-invariant polynomials of degrees 2 and 4 on the kernel space graded by fermion number F, by the exact
                          Molien-Weyl computation of z4 read at every F charge (exact rational arithmetic).
a6-5  thm:Bcounts         the same counts graded by (F, J Gamma charge): 8 + 228 invariants of U(1)_{1-Gamma}, 6 + 162 with F and J Gamma
                          charge both imposed, 31 + 31 quartics of F = 0 and J Gamma charge +-8 (exact; charges of the components from a8-9).
Each check prints PASS/FAIL a6-<n>; exit status 0 iff all pass.
"""
import contextlib
from collections import Counter
import io
import itertools
import os
import runpy
import sys
from fractions import Fraction
from functools import reduce

import numpy as np
import sympy as sp

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "common"))
from clifford8 import I8, bl, comm, e  # noqa: E402

NPASS, FAILS = [0], []


def report(cid, ok, text):
    if ok:
        NPASS[0] += 1
        print(f"PASS a6-{cid}: {text}")
    else:
        FAILS.append(cid)
        print(f"FAIL a6-{cid}: {text}")
    sys.stdout.flush()


# ---------------------------------------------------------------- a6-1: thm:subkernels (n = 1, exact)
ETA = {a: (1 if a <= 3 else -1) for a in range(1, 7)}
BLADES = [S for r in range(7) for S in itertools.combinations(range(1, 7), r)]
BM = {S: (bl(*S) if S else I8) for S in BLADES}
Z8 = sp.zeros(8, 8)


def antisym(S):
    return BM[S].T == -BM[S]


def anticomm_parity(S, T):  # e_S e_T = (-1)^p e_T e_S
    return (len(S) * len(T) - len(set(S) & set(T))) % 2


def kernel_and_algebra(Widx):
    gens = [(a, b) for a, b in itertools.combinations(Widx, 2)]
    rot = [g for g in gens if ETA[g[0]] == ETA[g[1]]]
    boo = [g for g in gens if ETA[g[0]] != ETA[g[1]]]
    M = [S for S in BLADES if antisym(S) and all(anticomm_parity(S, g) == 0 for g in rot) and all(anticomm_parity(S, g) == 1 for g in boo)]
    g = [S for S in BLADES if antisym(S) and all(anticomm_parity(S, x) == 0 for x in gens)]
    return M, g


def coords(X, basis):
    return [sp.Rational((BM[S].inv() * X).trace(), 8) for S in basis]


def action_matrices(M, g):
    """matrices of ad(X) restricted to span M, X in g (M and g spanned by blades; ad maps blades to multiples of blades)."""
    mats = []
    for X in g:
        cols = [coords(comm(BM[X], BM[S]), M) for S in M]
        mats.append(sp.Matrix(cols).T)
    return mats


def commutant(mats, n):
    xs = sp.symbols(f"c0:{n * n}")
    Y = sp.Matrix(n, n, xs)
    eqs = [x for A in mats for x in (A * Y - Y * A)]
    if not eqs:
        return [sp.Matrix(n, n, lambda i, j: 1 if i == j else 0)]
    A_, _ = sp.linear_eq_to_matrix(eqs, xs)
    return [sp.Matrix(n, n, list(v)) for v in A_.nullspace()]


def invariant_vectors(mats, n):
    if not mats:
        return [sp.eye(n)[:, i] for i in range(n)]
    return sp.Matrix.vstack(*mats).nullspace()


T3, S3 = (1, 2, 3), (4, 5, 6)
mixed = [W for W in itertools.combinations(range(1, 7), 3) if W not in (T3, S3)]
results = {}
ok = True
for name, W in [("so T", T3), ("so S", S3)] + [(f"so{W}", W) for W in mixed] + [("so(3,3)", tuple(range(1, 7))), ("interior so(1,3)", (1, 4, 5, 6))]:
    M, g = kernel_and_algebra(W)
    mats = action_matrices(M, g) if M and g else []
    if M and g:
        # g preserves M: every ad(X)(blade of M) lies in span M (the coordinates reproduce the commutator)
        for X in g:
            for S in M:
                Y = comm(BM[X], BM[S])
                ok &= Y == sum((c * BM[T] for c, T in zip(coords(Y, M), M)), Z8)
    com = commutant(mats, len(M)) if M else []
    inv = invariant_vectors(mats, len(M)) if M else []
    results[name] = (len(M), len(g), len(com), len(inv), M, g, mats, com, inv)
# so T and so S: dim 10, 10; M absolutely irreducible (commutant R); no invariant vector; g simple (perfect, centre 0) and M ~ adjoint
for name in ("so T", "so S"):
    dM, dg, dc, di, M, g, mats, com, inv = results[name]
    ok &= (dM, dg, dc, di) == (10, 10, 1, 0)
    gm = [BM[S] for S in g]
    span = sp.Matrix([list(x) for x in gm])
    der = sp.Matrix([list(comm(a, b)) for a, b in itertools.combinations(gm, 2)])
    ok &= der.rank() == 10 and span.rank() == 10 and sp.Matrix.vstack(span, der).rank() == 10
    # Hom_g(g, M) != 0: g acts on itself by ad; solve L ad_g(X) = ad_M(X) L
    adg = [sp.Matrix([coords(comm(BM[X], BM[Y]), g) for Y in g]).T for X in g]
    xs = sp.symbols("h0:100")
    L = sp.Matrix(10, 10, xs)
    A_, _ = sp.linear_eq_to_matrix([x for a, b in zip(adg, mats) for x in (L * a - b * L)], xs)
    ok &= len(A_.nullspace()) == 1
    # the compact form: the Killing form of g is negative definite (Tr(ad X ad X) < 0) -> sp(2) = so(5) compact real form of dimension 10
    kill = sp.Matrix(10, 10, lambda i, j: (adg[i] * adg[j]).trace())
    ok &= all(ev < 0 for ev in kill.eigenvals())
mixed_ok = True
for W in mixed:
    dM, dg, dc, di, M, g, mats, com, inv = results[f"so{W}"]
    mixed_ok &= (dM, dg, di) == (10, 6, 1)
    # the invariant kernel is +- gamma^0_W omega_W, gamma^0_W the generator of minority signature
    signs = [ETA[a] for a in W]
    minority = [a for a in W if signs.count(ETA[a]) == 1][0]
    target = BM[(minority,)] * BM[W]
    v = sum((c * BM[S] for c, S in zip(inv[0], M)), Z8)
    lam = [sp.nsimplify((target.T * v).trace() / (target.T * target).trace())][0]
    mixed_ok &= lam != 0 and v == lam * target
    # 9 + 1: the commutant is 2-dimensional and contains the projector onto the invariant line; the complement (9) is irreducible
    mixed_ok &= dc == 2
    Pinv = None
    for C in com:
        evs = C.eigenvals()
        if len(evs) == 2:
            Pinv = C
    mixed_ok &= Pinv is not None
ok &= mixed_ok
dM, dg, dc, di, M, g, mats, com, inv = results["so(3,3)"]
ok &= dM == 1 and dg == 0 and BM[M[0]] == e[1] * e[2] * e[3]
dM, dg, dc, di, M, g, mats, com, inv = results["interior so(1,3)"]
ok &= dM == 6 and dg == 4 and di == 0
# irreducible: the commutant is 2-dimensional and isomorphic to C (a traceless element squares to a negative scalar), so it contains no
# nontrivial idempotent and M has no proper invariant subspace
Cn = [C - (C.trace() / 6) * sp.eye(6) for C in com]
Cn = [C for C in Cn if C != sp.zeros(6)]
ok &= dc == 2 and len(Cn) >= 1 and (Cn[0] * Cn[0]) == (Cn[0] * Cn[0])[0, 0] * sp.eye(6) and (Cn[0] * Cn[0])[0, 0] < 0
gI = sorted(g)
J1, Om1 = e[2] * e[3], e[1] * e[4] * e[5] * e[6]
target = [J1, Om1, e[2] * Om1, e[3] * Om1]
ok &= sp.Matrix([list(BM[S]) for S in gI]).rank() == sp.Matrix([list(BM[S]) for S in gI] + [list(t) for t in target]).rank() == 4
report("1", ok, "thm:subkernels: so T and so S: dim M_W = dim g_W = 10, g_W perfect with negative-definite Killing form (sp(2)), M_W absolutely "
       "irreducible and isomorphic to the adjoint, no invariant kernel; each of the 18 mixed coordinate three-planes: dim M_W = 10, dim g_W = 6, "
       "M_W = 9 + 1 with the invariant line spanned by gamma^0_W omega_W (gamma^0_W the generator of minority signature); so(3,3): M = R e1e2e3, "
       "g = 0; interior so(1,3): dim M = 6, g_I = <J, Omega, e2 Omega, e3 Omega>, M irreducible (commutant C), no invariant kernel "
       "(exact rational arithmetic on blades)")

# ---------------------------------------------------------------- the 512 module, exact (integers; float64 BLAS only where every value is
# an integer below 2^52, checked)
sx = np.array([[0, 1], [1, 0]]); sz = np.array([[1, 0], [0, -1]]); ep = np.array([[0, -1], [1, 0]]); I2 = np.eye(2, dtype=np.int64)
kr = lambda *a: reduce(np.kron, a)  # noqa: E731
e8 = [None, kr(sx, I2, I2), kr(sz, sx, I2), kr(sz, sz, sx), kr(ep, I2, I2), kr(sz, ep, I2), kr(sz, sz, ep)]
I8n = np.eye(8, dtype=np.int64)
om = reduce(np.matmul, e8[1:])
E = [None] + [kr(*([om] * k + [e8[a]] + [I8n] * (2 - k))).astype(np.int64) for k in range(3) for a in range(1, 7)]
I512 = np.eye(512, dtype=np.int64)
Z512 = np.zeros((512, 512), dtype=np.int64)
blq = lambda S: reduce(np.matmul, [E[i] for i in S], I512)  # noqa: E731


def mm(*ms):
    """exact product of integer matrices through float64 BLAS (asserts that every partial sum is an exact integer below 2^52)"""
    out = ms[0].astype(np.float64)
    for m in ms[1:]:
        out = out @ m.astype(np.float64)
        assert np.abs(out).max() < 2 ** 52
    r = np.rint(out)
    assert np.array_equal(r, out)
    return r.astype(np.int64)


interior = [1, 4, 5, 6]
internal = [i for i in range(1, 19) if i not in interior]
Widx, Cidx = [7, 13, 14, 15], [10, 11, 12, 16, 17, 18]
J = blq(internal)
Gam = mm(J, blq(interior))
chiW, chiC, qA = blq(Widx), blq(Cidx), mm(E[2], E[3])
A8 = mm(I512 + mm(Gam, chiW), I512 - mm(J, chiC), I512 - mm(J, qA))          # 8 P_*
w4 = [E[i] for i in Widx]
tau2 = [-mm(w4[1], w4[2]), -mm(w4[2], w4[0]), -mm(w4[0], w4[1])]               # 2 tau_k
tL4 = [mm(I512 - chiW, t) for t in tau2]                                        # 4 t^L_k
tR4 = [mm(I512 + chiW, t) for t in tau2]                                        # 4 t^R_k
HC = -(mm(E[10], E[16]) + mm(E[11], E[17]) + mm(E[12], E[18]))
B3 = -mm(J, HC)                                                                 # 3 (B-L)
T4 = -mm(J, tR4[2])                                                             # 4 T^R_3 (charge operator)
# P_nu = P_* ((1+Gamma)/2) (2 T3R (T3R + 1/2)) ((1 - 3(B-L))/4) P_*  (right-handed, T3R = +1/2, lepton)
Nnu = mm(A8, A8, I512 + Gam, T4, T4 + 2 * I512, I512 - B3, A8)                  # = 2^15 P_nu * 64 ... scale below
SC_NU = 8 * 8 * 2 * 8 * 4 * 8                                                   # A8 A8 (I+G) T4(T4+2) (I-B3) A8 -> P_nu * SC_NU
Pnu_check = mm(Nnu, Nnu)
ok_nu = np.array_equal(Pnu_check, SC_NU * Nnu) and np.trace(Nnu) == SC_NU * 4
X = mm(E[1], blq([3, 9, 13, 14, 16, 17, 18]))
PHI = mm(Nnu, A8, X, A8, Nnu)                                                    # = Phi_0 * SC_PHI
SC_PHI = SC_NU * 8 * 8 * SC_NU
G0 = int(np.gcd.reduce(np.abs(PHI[PHI != 0])))
PHI = PHI // G0
SC_PHI = Fraction(SC_PHI, G0)                                                    # PHI = SC_PHI * Phi_0 with a small integer PHI


def kernel_from(PHIs):
    return PHIs


C_op = blq([2, 7, 10, 11, 12])
P_op = blq([1, 7, 8])
CP_op = blq([1, 2, 7, 8, 10, 11, 12, 14])
conj = lambda O, M: mm(O, M, O.T)  # noqa: E731
ok = ok_nu and PHI.T.tolist() == (-PHI).tolist() and np.any(PHI)
ok &= np.array_equal(conj(CP_op, PHI), PHI)
CPHI, PPHI = conj(C_op, PHI), conj(P_op, PHI)
# colour su(3): the derived algebra of the centraliser of HC in colour so(6), computed exactly (no hand-written generators)
colour15 = [mm(E[a], E[b]) for a, b in itertools.combinations(Cidx, 2)]
flat = np.array([(mm(c, HC) - mm(HC, c)).reshape(-1) for c in colour15], dtype=np.float64)
gram = flat @ flat.T
assert np.array_equal(gram, np.rint(gram)) and np.abs(gram).max() < 2 ** 52
nsp = sp.Matrix(np.rint(gram).astype(np.int64)).nullspace()
cent = []
for v in nsp:
    den = sp.ilcm(*[sp.fraction(x)[1] for x in v])
    cent.append(sum(int(x * den) * c for x, c in zip(v, colour15)))
derived = [mm(a, b) - mm(b, a) for a, b in itertools.combinations(cent, 2)]
dflat = np.array([d.reshape(-1) for d in derived], dtype=np.float64)
dgram = dflat @ dflat.T
assert np.array_equal(dgram, np.rint(dgram)) and np.abs(dgram).max() < 2 ** 52
su3 = derived                                                                    # spans su(3); commuting with all of them = commuting with su(3)
ok &= len(cent) == 9 and sp.Matrix(np.rint(dgram).astype(np.int64)).rank() == 8
Qgen12 = 3 * tL4[2] + 3 * tR4[2] + 2 * HC                                       # 12 (t^L_3 + t^R_3 + (B-L)gen/2), (B-L)gen = HC/3
commutes = lambda M, g: np.array_equal(mm(M, g), mm(g, M))  # noqa: E731
# sanity: Phi_0 itself commutes with su3, su2_L and Q, not with su2_R
ok &= all(commutes(PHI, s) for s in su3) and all(commutes(PHI, t) for t in tL4) and commutes(PHI, Qgen12)
ok &= not all(commutes(PHI, t) for t in tR4)
ok &= all(commutes(CPHI, s) for s in su3) and all(commutes(CPHI, t) for t in tR4) and commutes(CPHI, Qgen12)
ok &= not all(commutes(CPHI, t) for t in tL4)
# supports: nu_L nu_L and e_L e_L (lepton, Gamma = -1, T3L = +1/2 resp. -1/2)
TL4 = -mm(J, tL4[2])                                                            # 4 T^L_3
lep4 = I512 - B3                                                                # 4 P_lepton on the block
left2 = I512 - Gam                                                              # 2 P_left
up8 = mm(TL4, TL4 + 2 * I512)                                                   # 8 P(T3L = +1/2) on doublets
dn8 = mm(TL4, TL4 - 2 * I512)                                                   # 8 P(T3L = -1/2)
P_nuL = mm(A8, lep4, left2, up8, A8)                                            # 8*4*2*8*8 P_nuL
P_eL = mm(A8, lep4, left2, dn8, A8)
SCP = 8 * 4 * 2 * 8 * 8
ok &= np.array_equal(mm(P_nuL, P_nuL), SCP * P_nuL) and np.trace(P_nuL) == SCP * 4 and np.array_equal(mm(P_eL, P_eL), SCP * P_eL)
ok &= np.array_equal(mm(P_nuL, CPHI, P_nuL), SCP * SCP * CPHI) and np.any(CPHI)
ok &= np.array_equal(mm(P_eL, PPHI, P_eL), SCP * SCP * PPHI) and np.any(PPHI)
# both in (10,3,1): su(2)_R singlets, su(2)_L triplets (Casimir 2), colour Casimir equal to that of Phi_0 (the 10)
def cas(gens, M, scale):
    acc = sum(mm(g, mm(g, M) - mm(M, g)) - mm(mm(g, M) - mm(M, g), g) for g in gens)
    return -acc, scale
for Mx in (CPHI, PPHI):
    ok &= all(commutes(Mx, t) for t in tR4)
    cL, _ = cas(tL4, Mx, 16)
    ok &= np.array_equal(cL, 16 * 2 * Mx)
cc_phi, _ = cas(colour15, PHI, 1)
cc_c, _ = cas(colour15, CPHI, 1)
lam_phi = Fraction(int(np.sum(cc_phi * PHI)), int(np.sum(PHI * PHI)))
lam_c = Fraction(int(np.sum(cc_c * CPHI)), int(np.sum(CPHI * CPHI)))
ok &= np.array_equal(cc_phi * lam_phi.denominator, lam_phi.numerator * PHI) and lam_phi == lam_c
ok &= np.array_equal(cc_c * lam_c.denominator, lam_c.numerator * CPHI)
# nu_L nu_L and e_L e_L on one orbit: R = (1 - P_L') + 2 t^L_1 (the su(2)_L rotation by pi about the first axis) maps one to the other
PLp2 = I512 - chiW                                                              # 2 P_L'
R2 = 2 * I512 - PLp2 + tL4[0]                                                   # 2 R
ok &= np.array_equal(mm(R2, R2.T), 4 * I512)
img = mm(R2, CPHI, R2.T)
ok &= np.array_equal(img, 4 * PPHI) or np.array_equal(img, -4 * PPHI)
report("2", ok, "thm:vacCP: Phi_0 rebuilt exactly (P_nu idempotent of rank 4); CP Phi_0 CP^T = Phi_0 for CP = E1E2E7E8E10E11E12E14; C Phi_0 C^T "
       "(C = E2E7E10E11E12) is supported on nu_L and commutes with su(3), su(2)_R and Q but not with su(2)_L; P Phi_0 P^T (P = E1E7E8) is supported "
       f"on e_L; both are su(2)_R singlets and su(2)_L triplets with the colour Casimir of Phi_0 ({lam_phi}), i.e. in (10,3,1), and an su(2)_L "
       "rotation by pi maps one to the other (exact integer arithmetic)")

# ---------------------------------------------------------------- a6-3: prop:coupled3 (algebraic identities on the block, exact)
K = [mm(A8, E[1], E[a], A8) for a in (4, 5, 6)]                                 # 64 K^i
ok = all(np.array_equal(mm(K[i], K[j]) + mm(K[j], K[i]), (2 * 64 * 8 if i == j else 0) * A8) for i in range(3) for j in range(3))
ok &= all(np.array_equal(mm(Ki, PHI) + mm(PHI, Ki), Z512) for Ki in K)
lhs = -mm(PHI, PHI) * SC_NU                                                       # SC_NU * SC_PHI^2 * (-Phi_0^2)
rhs = SC_PHI * SC_PHI                                                             # times Nnu = SC_NU P_nu
ok &= rhs.denominator == 1 and np.array_equal(lhs, int(rhs) * Nnu) if rhs.denominator == 1 else \
    np.array_equal(lhs * rhs.denominator, rhs.numerator * Nnu)
report("3", ok, "prop:coupled3 (algebraic core): on the block {K^i, K^j} = 2 delta^ij P_*, {K^i, Phi_0} = 0 and -Phi_0^2 = P_nu, so "
       "H(p)^2 = |p|^2 + y^2 r^2 P_nu for H(p) = p_i K^i + i y r Phi_0 (exact integer arithmetic; the Euler-Lagrange equations and the "
       "inherited one-loop Coleman-Weinberg formula are not part of this check)")

# ---------------------------------------------------------------- a6-4: thm:Bcounts, F-graded (exact Molien-Weyl, reusing z4)
src = open(os.path.join(HERE, "z4_invariant_counts.py")).read()
src = src[:src.index("c0 = counts(False)")]            # definitions only (the z4 checks themselves run in z4)
z4 = {"__file__": os.path.join(HERE, "z4_invariant_counts.py"), "__name__": "z4_definitions"}
with contextlib.redirect_stdout(io.StringIO()):
    exec(compile(src, "z4_invariant_counts.py", "exec"), z4)
KC, scale, mult, WF, W_ORDER = z4["KC"], z4["scale"], z4["mult"], z4["WF"], z4["W_ORDER"]


def coeff_F(poly, f):
    """sum over Weyl-factor weights w of WF(w) * poly(-w_PS, F = f), divided by |W|: the number of PS invariants of F-charge f"""
    s = 0
    for w, c in WF.items():
        key = tuple(-x for x in w[:5]) + (f,)
        s += c * poly.get(key, 0)
    return Fraction(s, W_ORDER)


ch = dict(KC)
p = {k: scale(ch, k) for k in (1, 2, 3, 4)}
p11 = mult(p[1], p[1])
h2 = {}
for poly, cf in ((p11, Fraction(1, 2)), (p[2], Fraction(1, 2))):
    for k, v in poly.items():
        h2[k] = h2.get(k, 0) + cf * v
p1111, p112, p22, p13 = mult(p11, p11), mult(p11, p[2]), mult(p[2], p[2]), mult(p[1], p[3])
h4 = {}
for poly, cf in ((p1111, Fraction(1, 24)), (p112, Fraction(1, 4)), (p22, Fraction(1, 8)), (p13, Fraction(1, 3)), (p[4], Fraction(1, 4))):
    for k, v in poly.items():
        h4[k] = h4.get(k, 0) + cf * v
Fvals = sorted({k[5] for k in h4})
d2 = {f_: coeff_F(h2, f_) for f_ in sorted({k[5] for k in h2})}
d4 = {f_: coeff_F(h4, f_) for f_ in Fvals}
d2 = {k: v for k, v in d2.items() if v}
d4 = {k: v for k, v in d4.items() if v}
unit = min(abs(k) for k in d2 if k) if any(d2) else 1
ok = sum(d2.values()) == 18 and sum(d4.values()) == 686 and d2.get(0) == 12 and d4.get(0) == 402
nz2 = sorted((k, int(v)) for k, v in d2.items() if k)
nz4 = sorted((k, int(v)) for k, v in d4.items() if k)
ok &= [v for _, v in nz2] == [3, 3] and sorted(v for _, v in nz4) == [15, 15, 127, 127]
ok &= all(v.denominator == 1 for v in list(d2.values()) + list(d4.values()))
report("4", ok, f"thm:Bcounts (F grading): degree 2 by F charge {sorted((k, int(v)) for k, v in d2.items())}, degree 4 {sorted((k, int(v)) for k, v in d4.items())} "
       "(F in the units of z4): 12 + 3 + 3 = 18 and 402 + 127 + 127 + 15 + 15 = 686 as stated (exact Molien-Weyl integrals)")

# ---------------------------------------------------------------- a6-5: thm:Bcounts, U(1)_{1-Gamma} and (F, J Gamma) gradings (exact Molien-Weyl)
# Charges of the components under ad J (fermion number F, z4 units) and ad J Gamma (q), relative to the complex structure of V.  They are
# read from the realisation by a8-9 (verification/ch07_mass_vacuum/a8_electroweak_numerical.py): on (10,1,3) ad(J Gamma) = +ad J, on
# (10,3,1) ad(J Gamma) = -ad J, on 2(6,1,1) +ad J on one copy and -ad J on the other, and on the Dirac-type (15,2,2), (1,2,2) ad J = 0 and
# ad(J Gamma) a complex structure (charge +-2; these representations are self-conjugate, so the sign is immaterial).
# U(1)_{1-Gamma} is generated by J - J Gamma (Chapter 7), acting on kernels as ad J - ad J Gamma: its invariants are those with F = q.
w10, w6, w15, w1 = z4["w10"], z4["w6"], z4["w15"], z4["w1"]
sp1, sp12, sp0 = [-2, 0, 2], [-1, 1], [0]


def irrep7(ws4, wl, wr, f, q):
    return [a + (l_, r_, f, q) for a in ws4 for l_ in wl for r_ in wr]


V7 = (irrep7(w10, sp1, sp0, 2, -2) + irrep7(w10, sp0, sp1, 2, 2) + irrep7(w15, sp12, sp12, 0, 2) + irrep7(w6, sp0, sp0, 2, 2)
      + irrep7(w6, sp0, sp0, 2, -2) + irrep7(w1, sp12, sp12, 0, 2))
KC7 = Counter(V7 + [tuple(-x for x in w) for w in V7])
ok = len(V7) == 136 and sum(KC7.values()) == 272
ok &= Counter(w[:6] for w in KC7.elements()) == Counter(KC.elements())               # forgetting q gives the character of z4
WF7 = {w + (0,): c for w, c in WF.items()}


def coeff_Fq(poly, f, q):
    s = 0
    for w, c in WF7.items():
        s += c * poly.get(tuple(-x for x in w[:5]) + (f, q), 0)
    return Fraction(s, W_ORDER)


ch7 = dict(KC7)
p7 = {k: scale(ch7, k) for k in (1, 2, 3, 4)}
q11 = mult(p7[1], p7[1])
h2q, h4q = {}, {}
for poly, cf in ((q11, Fraction(1, 2)), (p7[2], Fraction(1, 2))):
    for k, v in poly.items():
        h2q[k] = h2q.get(k, 0) + cf * v
for poly, cf in ((mult(q11, q11), Fraction(1, 24)), (mult(q11, p7[2]), Fraction(1, 4)), (mult(p7[2], p7[2]), Fraction(1, 8)),
                 (mult(p7[1], p7[3]), Fraction(1, 3)), (p7[4], Fraction(1, 4))):
    for k, v in poly.items():
        h4q[k] = h4q.get(k, 0) + cf * v
g2 = {fq: coeff_Fq(h2q, *fq) for fq in sorted({k[5:] for k in h2q})}
g4 = {fq: coeff_Fq(h4q, *fq) for fq in sorted({k[5:] for k in h4q})}
g2 = {k: v for k, v in g2.items() if v}
g4 = {k: v for k, v in g4.items() if v}
ok &= all(v.denominator == 1 for v in list(g2.values()) + list(g4.values()))
ok &= sum(g2.values()) == 18 and sum(g4.values()) == 686
ok &= sum(v for (f_, _), v in g2.items() if f_ == 0) == 12 and sum(v for (f_, _), v in g4.items() if f_ == 0) == 402
n2_1g = sum(v for (f_, q_), v in g2.items() if f_ == q_)
n4_1g = sum(v for (f_, q_), v in g4.items() if f_ == q_)
n2_fq = g2.get((0, 0), 0)
n4_fq = g4.get((0, 0), 0)
n4_q8 = (g4.get((0, 8), 0), g4.get((0, -8), 0))
ok &= (n2_1g, n4_1g, n2_fq, n4_fq) == (8, 228, 6, 162) and n4_q8 == (31, 31)
report("5", ok, f"thm:Bcounts (U(1)_{{1-Gamma}} and J Gamma gradings): with charges (F, q) under (ad J, ad J Gamma) the invariants of degree 2 "
       f"and 4 with F = q (U(1)_{{1-Gamma}}) number {int(n2_1g)} + {int(n4_1g)}, and those with F = q = 0 number {int(n2_fq)} + {int(n4_fq)}, as "
       f"stated; the F = 0 quartics of J Gamma charge +-8 number {int(n4_q8[0])} + {int(n4_q8[1])} (thm:mA discussion) (exact Molien-Weyl "
       "integrals; the component charges are read from the realisation in a8-9)")

print(f"{NPASS[0]}/{NPASS[0] + len(FAILS)} checks passed")
sys.exit(0 if not FAILS else 1)
