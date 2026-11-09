#!/usr/bin/env python3
"""Exact certificates for the tensor model of Chapter 7 (errata items 59 and 60).

    python3 verification/ch07_mass_vacuum/a13_tensor_certificates.py

a13-1  thm:potential (ii)  The holomorphic quartic D depends on S only through the (20',1) component w of S (x) S, and
                           |D(S)| <= (3/8) I_4(S) with equality at 1 (x) 1, the constant computed exactly (item 60).
a13-2  prop:tensor         An explicit Pati-Salam-equivariant isomorphism iota from the (10,1,3) kernel component to the tensor space T,
                           with Tr Phi^T Phi = 4 |iota(Phi)|^2, iota(Phi_0) = S_0 and Phi^T Phi = (S^dagger S) (x) 1_4 in an adapted basis (item 59).

Conventions are those of v1_tensor_model.py: T = Sym^2(C^4) (x) Sym^2(C^2) inside the complex symmetric 8 x 8 matrices, index (a,s) -> 2a+s,
|S|^2 = Tr S^dagger S, P_{ab,ef}(S) = eps_su eps_tv S_{as,et} S_{bu,fv}, D(S) = (1/24) eps^{abcd} eps^{efgh} P_{ab,ef} P_{cd,gh}.

a13-1, the method. Sym^2 T = (35,5) + (35,1) + (20',5) + (20',1) + (45,3) (inherited plethysm), so the (20',1) component W is the tensor product
of the su(2) singlet s_{st,uv} = eps_su eps_tv + eps_sv eps_tu of Sym^2(Sym^2 C^2) with the 20' of Sym^2(Sym^2 C^4), the tensors T_{ab,cd} symmetric
in (ab), in (cd) and under the pair exchange whose total symmetrisation vanishes. W has an integer basis B (20 vectors of the 4096-dimensional
space (C^8)^{(x)4}) and integer Gram matrix G. D is the restriction to S (x) S of the bilinear form D~(X,Y) = P(X)^T E P(Y) on T (x) T. Exact
integer identities certify: (a) D~(X,Y) = D~(pi X, pi Y) for X, Y in Sym^2 T, pi the orthogonal projector onto W, on a basis of the 465-dimensional
Sym^2 T; (b) Q G^{-1} Q = (3/8)^2 G for the matrix Q of D~ on B, so every Takagi value of D~ on W, measured with the norm of W, is 3/8 and
|D~(w,w)| <= (3/8)|w|^2 for every complex w in W; (c) D(1 (x) 1) = 4 and I_4(1 (x) 1) = |pi(1 (x) 1)|^2 = 32/3. All arithmetic is integer or rational.

a13-2, the method. V_+ = im P_*(1+Gamma) has an integer basis R (32 columns of 16 P_* P_+). An intertwiner U: V_+ -> C^2 (x) C^4 (x) C^2 with
U M_g = rho(g) U for the six Lorentz bivectors, the fifteen colour bivectors, the three su(2)_R generators and J -> i is found as a numerical
null space (dimension one, for the chiralities realised on V_+), rationalised, and then verified in exact rational arithmetic, together with
Re U^dagger U = H/16 (H = R^T R). For each S of a real basis of T the real antilinear map v -> (eps (x) S) conj(v), pulled back by U, is verified
to be an antisymmetric Lorentz-invariant kernel; Tr Phi_i^T Phi_j = 4 Re Tr S_i^dagger S_j on the basis; Phi_0 (exact, up to a positive factor)
lies in their span with a rank-one image of B-L = -2, T^R_3 = 1; and Phi_0^T Phi_0 equals the realification of 1 (x) S S^dagger. Since the
kernel space contains (10,1,3) once (inherited count, v2-1), the span is that component, and iota(Phi) = S.
"""
import itertools
import os
import sys
from fractions import Fraction
from functools import reduce

import numpy as np
import sympy as sp
from sympy.polys.domains import QQ
from sympy.polys.matrices import DomainMatrix

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

FAILS = []


def report(cid, ok, text):
    print(f"{'PASS' if ok else 'FAIL'} {cid}: {text}")
    if not ok:
        FAILS.append(cid)


def objmat(a):
    return np.array(a, dtype=object)


def iadj_det(G):
    """integer adjugate and determinant of an integer matrix (exact, sympy)."""
    M = sp.Matrix(G.tolist())
    d = M.det()
    return objmat((M.adjugate()).tolist()), int(d)


# ------------------------------------------------------------------------------------------------ tensors
eps2 = np.array([[0, 1], [-1, 0]], dtype=np.int64)
eps4 = np.zeros((4, 4, 4, 4), dtype=np.int64)
for p in itertools.permutations(range(4)):
    eps4[p] = round(np.linalg.det(np.eye(4)[list(p)]))

# s_{st,uv} = eps_su eps_tv + eps_sv eps_tu
s2 = np.einsum("su,tv->stuv", eps2, eps2) + np.einsum("sv,tu->stuv", eps2, eps2)

# the 20' inside Sym^2(Sym^2 C^4): integer basis by exact null space
sym_pairs = [(a, b) for a in range(4) for b in range(a, 4)]                  # basis of Sym^2 C^4 (10)
sym_sq = [(i, j) for i in range(10) for j in range(i, 10)]                    # basis of Sym^2(Sym^2 C^4) (55)


def unit_sym(a, b):
    M = np.zeros((4, 4), dtype=np.int64)
    M[a, b] += 1
    M[b, a] += 1 if a != b else 0
    return M


def T4(i, j):
    A, Bm = unit_sym(*sym_pairs[i]), unit_sym(*sym_pairs[j])
    return np.einsum("ab,cd->abcd", A, Bm) + (np.einsum("ab,cd->abcd", Bm, A) if i != j else 0)


cand = [T4(i, j) for i, j in sym_sq]
perms = list(itertools.permutations(range(4)))
tot = lambda X: sum(np.transpose(X, p) for p in perms)                          # noqa: E731  total symmetrisation (times 24)
Mtot = sp.Matrix([[int(x) for x in tot(c).ravel()] for c in cand]).T            # 256 x 55
ns = Mtot.nullspace()
ok_dim = len(ns) == 20
W4 = []
for v in ns:
    den = sp.ilcm(*[sp.Rational(x).q for x in v])
    v = [int(x * den) for x in v]
    W4.append(sum(c * X for c, X in zip(v, cand)))                              # integer tensors T_{ab,cd} in the 20'

# W in (C^8)^{(x)4}, index order [a,s,b,t,c,u,d,v] for X_{(as),(bt),(cu),(dv)}
Bw = [np.einsum("abcd,stuv->asbtcudv", X, s2).reshape(-1).astype(np.int64) for X in W4]
Bm_ = np.array(Bw).T                                                            # 4096 x 20
G = objmat((Bm_.T @ Bm_).tolist())
adjG, detG = iadj_det(G)

# Sym^2 T: integer basis of T (30) and its symmetrised products (465)
Tbasis = []
for (a, b) in sym_pairs:
    for (s, t) in [(0, 0), (0, 1), (1, 1)]:
        A, Bs = unit_sym(a, b), np.zeros((2, 2), dtype=np.int64)
        Bs[s, t] += 1
        Bs[t, s] += 1 if s != t else 0
        Tbasis.append(np.einsum("ab,st->asbt", A, Bs))                          # S_{as,bt}
Ys = []
for i in range(30):
    for j in range(i, 30):
        Y = np.einsum("asbt,cudv->asbtcudv", Tbasis[i], Tbasis[j]) + np.einsum("asbt,cudv->asbtcudv", Tbasis[j], Tbasis[i])
        Ys.append(Y.reshape(-1))
Ym = np.array(Ys).T                                                              # 4096 x 465


def Pmap(Xflat):
    """P(X)_{ab,ef} = eps_su eps_tv X_{(as),(et),(bu),(fv)}, X with index order [a,s,e,t,b,u,f,v]."""
    X = Xflat.reshape(4, 2, 4, 2, 4, 2, 4, 2)
    return np.einsum("su,tv,asetbufv->abef", eps2, eps2, X).reshape(-1)


E24 = np.einsum("abcd,efgh->abefcdgh", eps4, eps4).reshape(256, 256)          # 24 E: D~(X,Y) = P(X)^T E24 P(Y) / 24

LY = np.array([Pmap(Ym[:, k]) for k in range(Ym.shape[1])]).T                 # 256 x 465
LB = np.array([Pmap(Bm_[:, k]) for k in range(20)]).T                           # 256 x 20
M1 = objmat((LY.T @ E24 @ LY).tolist())                                         # 24 D~ on Sym^2 T
Q24 = objmat((LB.T @ E24 @ LB).tolist())                                        # 24 D~ on W
C = objmat((Bm_.T @ Ym).tolist())                                               # coordinates <B_i, Y_k>
# (a) 24 D~(Y_k, Y_l) = 24 D~(pi Y_k, pi Y_l), pi Y = B G^{-1} B^T Y:   det^2 * M1 == C^T adj Q24 adj C
lhs = M1 * (detG * detG)
rhs = C.T.dot(adjG).dot(Q24).dot(adjG).dot(C)
ok_a = bool(np.all(lhs == rhs))
# (b) Q G^{-1} Q = (9/64) G  <=>  Q24 adj Q24 = 576 (9/64) det G * G = 81 det G * G
ok_b = bool(np.all(Q24.dot(adjG).dot(Q24) == G * (81 * detG)))
# (c) the point 1 (x) 1
one = np.einsum("ab,st->asbt", np.eye(4, dtype=np.int64), np.eye(2, dtype=np.int64)).reshape(8, 8)
X1 = np.einsum("ij,kl->ijkl", one, one).reshape(-1)
p1 = Pmap(X1)
D1 = sp.Rational(int(p1 @ E24 @ p1), 24)
c1 = objmat((Bm_.T @ X1).tolist())
I4 = sp.Rational(int(c1.dot(adjG).dot(c1)), detG)
ok_c = D1 == 4 and I4 == sp.Rational(32, 3) and abs(D1) == sp.Rational(3, 8) * I4
# sanity: the 1 (x) 1 value of D agrees with the direct contraction of the manuscript
S4 = one.reshape(4, 2, 4, 2)
Pt = np.einsum("su,tv,aset,bufv->abef", eps2, eps2, S4, S4)
ok_c &= sp.Rational(int(np.einsum("abcd,efgh,abef,cdgh->", eps4, eps4, Pt, Pt)), 24) == D1
report("a13-1", ok_dim and ok_a and ok_b and ok_c,
       f"thm:potential (ii), exact: the 20' of Sym^2(Sym^2 C^4) has dimension {len(ns)} and with the su(2) singlet spans the (20',1) component W "
       f"(integer basis, Gram determinant {detG}); on the 465-dimensional Sym^2 T the quartic form satisfies D~(X,Y) = D~(pi X, pi Y), pi the projector "
       f"onto W ({ok_a}); its matrix Q on W satisfies Q G^-1 Q = (3/8)^2 G ({ok_b}), so |D(S)| <= (3/8) I_4(S) for every S; D(1(x)1) = {D1} and "
       f"I_4(1(x)1) = {I4}, so the constant 3/8 is attained (integer and rational arithmetic)")



# ================================================================================================ a13-2  prop:tensor (item 59)
# V_+ = im P_* (1+Gamma), real dimension 32, with the complex structure J.  An intertwiner U: V_+ -> C^2 (x) C^4 (x) C^2 (Lorentz Weyl
# (x) colour (x) su(2)_R) is found from the linear equations, rationalised, and then verified exactly; every kernel of the (10,1,3)
# component is then v -> (eps (x) S) conj(v) in the model, and iota(Phi) = S.
import kernel_space as ks  # noqa: E402

Ei = ks.E
I512 = np.eye(512, dtype=np.int64)
P8 = (I512 + ks.Xi) @ (I512 + ks.Zi) @ (I512 + ks.Ui)                          # 8 P_*
Pp = P8 @ (I512 + ks.Gam)                                                      # 16 P_* P_+ (integer)
cols = []
Pf = Pp.astype(float)
for j in range(512):
    if np.linalg.matrix_rank(Pf[:, cols + [j]]) > len(cols):
        cols.append(j)
    if len(cols) == 32:
        break
R = Pp[:, cols]                                                                # integer basis of V_+
dm = lambda A: DomainMatrix([[QQ(int(x)) for x in row] for row in np.asarray(A).tolist()], np.asarray(A).shape, QQ)  # noqa: E731
H = dm(R.T @ R)
Hinv = H.inv()
ok_basis = dm(R).rank() == 32 and int(np.linalg.matrix_rank(Pf)) == 32


def rep(g):
    """matrix of the integer operator g on V_+ in the basis R (exact; asserts invariance)."""
    M = Hinv * dm(R.T @ g @ R)
    assert dm(g @ R) == dm(R) * M
    return M


def gammas(m, signs):
    s1 = np.array([[0, 1], [1, 0]]); s2 = np.array([[0, -1j], [1j, 0]]); s3 = np.diag([1, -1]); i2 = np.eye(2)
    kr = lambda *a: reduce(np.kron, a)  # noqa: E731
    gs = []
    for k in range(m):
        gs += [kr(*([s3] * k + [s1] + [i2] * (m - k - 1))), kr(*([s3] * k + [s2] + [i2] * (m - k - 1)))]
    return [g * (1 if s > 0 else 1j) for g, s in zip(gs, signs)]


# model Clifford generators with the squares of e_1,e_4,e_5,e_6 (+,-,-,-), of the colour directions (all -1) and of the weak ones (all +1)
gL, gC, gW = gammas(2, [1, -1, -1, -1]), gammas(3, [-1] * 6), gammas(2, [1] * 4)


def half(gs, val):
    d = np.diag(reduce(np.matmul, gs))
    return np.where(np.isclose(d, val))[0]


hL, hC, hW = half(gL, -1j), half(gC, 1j), half(gW, 1)                       # the chiralities realised on V_+ (the others admit no intertwiner)
kr3 = lambda XL, XC, XW: np.kron(np.kron(XL, XC), XW)                          # noqa: E731
I2c, I4c = np.eye(2), np.eye(4)
mL = lambda X: kr3(X[np.ix_(hL, hL)], I4c, I2c)                                # noqa: E731
mC = lambda X: kr3(I2c, X[np.ix_(hC, hC)], I2c)                                # noqa: E731
mW = lambda X: kr3(I2c, I4c, X[np.ix_(hW, hW)])                                # noqa: E731
Lint = [1, 4, 5, 6]
gens = []                                                                      # (block operator, model operator)
for a, b in itertools.combinations(range(4), 2):
    gens.append((Ei[Lint[a]] @ Ei[Lint[b]], mL(gL[a] @ gL[b])))
for a, b in itertools.combinations(range(6), 2):
    gens.append((Ei[ks.Cidx[a]] @ Ei[ks.Cidx[b]], mC(gC[a] @ gC[b])))
w4 = [Ei[i] for i in ks.Widx]
for (j, k) in ((1, 2), (2, 0), (0, 1)):                                        # 2 tau_i = -w_j w_k (su(2)_R on V_+)
    gens.append((-(w4[j] @ w4[k]), mW(-(gW[j] @ gW[k]))))
gens.append((ks.J, 1j * np.eye(16)))
Ms = [rep(g) for g, _ in gens]
# numerical null space of U M_g = rho(g) U, then rational reconstruction
Mf = [np.array(M.to_Matrix().tolist(), dtype=float) for M in Ms]
A = np.vstack([np.kron(np.eye(16), M.T) - np.kron(r, np.eye(32)) for M, (_, r) in zip(Mf, gens)])
_, sv, Vh = np.linalg.svd(A)
null_dim = int(np.sum(sv < 1e-8))
u = Vh[-1].conj().reshape(16, 32)
u = u / u.flat[np.argmax(np.abs(u))]
rat = lambda x: Fraction(x).limit_denominator(64)                               # noqa: E731
Ure = [[rat(x.real) for x in row] for row in u]
Uim = [[rat(x.imag) for x in row] for row in u]
dmq = lambda rows: DomainMatrix([[QQ(x.numerator, x.denominator) for x in r] for r in rows], (len(rows), len(rows[0])), QQ)  # noqa: E731
Ur, Ui_ = dmq(Ure), dmq(Uim)


def cplx(X):
    """real and imaginary parts of a model matrix with entries in Z[i]/2 (exact)."""
    re = [[Fraction(round(2 * z.real), 2) for z in row] for row in X]
    im = [[Fraction(round(2 * z.imag), 2) for z in row] for row in X]
    assert np.allclose(np.array(re, float) + 1j * np.array(im, float), X)
    return dmq(re), dmq(im)


ok_U = null_dim == 1
for M, (_, r) in zip(Ms, gens):                                                # U M = r U, real and imaginary parts
    rr, ri = cplx(r)
    ok_U &= (Ur * M == rr * Ur - ri * Ui_) and (Ui_ * M == rr * Ui_ + ri * Ur)
# U is a scaled isometry: Re(U^dagger U) = lam H
G_U = Ur.transpose() * Ur + Ui_.transpose() * Ui_
lam = G_U.to_Matrix()[0, 0] / H.to_Matrix()[0, 0]
ok_iso = G_U == H * QQ(int(sp.Rational(lam).p), int(sp.Rational(lam).q))
Ureal = DomainMatrix.vstack(Ur, Ui_)                                           # R^32 (V_+ coordinates) -> R^32 = C^16 as (Re; Im)
Ureal_inv = Ureal.inv()
Cj = DomainMatrix.diag([QQ(1)] * 16 + [QQ(-1)] * 16, QQ)


def realify(Are, Aim):
    return DomainMatrix.vstack(DomainMatrix.hstack(Are, -Aim), DomainMatrix.hstack(Aim, Are))


eps = np.array([[0, 1], [-1, 0]])


def F_of(Sre, Sim):
    """the kernel on V_+ (coordinates of R) whose model form is v -> (eps (x) S) conj(v)."""
    return Ureal_inv * realify(dm(np.kron(eps, Sre)), dm(np.kron(eps, Sim))) * Cj * Ureal


Tb = [np.einsum("ab,st->asbt", Mx, Bx).reshape(8, 8) for Mx in [unit_sym(a, b) for (a, b) in sym_pairs]
      for Bx in [np.array([[2, 0], [0, 0]]), np.array([[0, 1], [1, 0]]), np.array([[0, 0], [0, 2]])]]
Z8 = np.zeros((8, 8), dtype=np.int64)
basis60 = [(S, Z8) for S in Tb] + [(Z8, S) for S in Tb]                        # real basis of T: S and i S
Fs = [F_of(Sre, Sim) for Sre, Sim in basis60]
Lor = Ms[:6]
ok_ker = True
for F in Fs:
    K = H * F                                                                  # matrix of the bilinear form on V_+
    ok_ker &= K.transpose() == -K                                              # antisymmetric kernel
    ok_ker &= all((Mx.transpose() * K + K * Mx).is_zero_matrix for Mx in Lor)   # Lorentz invariant
# Tr Phi_i^T Phi_j = c Re Tr S_i^dagger S_j on the basis (polarised norm identity)
tr = lambda D: sum(D.to_Matrix()[i, i] for i in range(D.shape[0]))             # noqa: E731
Fstar = [Hinv * F.transpose() * H for F in Fs]
GF = sp.Matrix(60, 60, lambda i, j: tr(Fstar[i] * Fs[j]))
GS = sp.Matrix(60, 60, lambda i, j: int(np.trace(basis60[i][0].T @ basis60[j][0] + basis60[i][1].T @ basis60[j][1])))
cst = GF[0, 0] / GS[0, 0]
ok_norm = GF == cst * GS and GS.rank() == 60
# Phi_0 (exact, up to a positive factor) and its image
J, Gam = ks.J, ks.Gam
tR3x2 = ((I512 + (w4[0] @ w4[1] @ w4[2] @ w4[3])) @ (-(w4[0] @ w4[1])))        # 4 t^R_3 (integer); chi_W = w_7 w_13 w_14 w_15
HCi = -(Ei[10] @ Ei[16] + Ei[11] @ Ei[17] + Ei[12] @ Ei[18])
TR3x4 = -J @ tR3x2                                                              # 4 T^R_3
BmLx3 = -J @ HCi                                                                # 3 (B-L)
Pnu = P8 @ (I512 + Gam) @ (TR3x4 @ (TR3x4 + 2 * I512)) @ (I512 - BmLx3) @ P8     # positive multiple of P_nu (integer)
Pnu = P8 @ Pnu @ P8
assert np.abs(Pnu).max() < 2 ** 40
Mphi = Ei[1] @ reduce(np.matmul, [Ei[i] for i in [3, 9, 13, 14, 16, 17, 18]])
Phi0 = Pnu @ P8 @ Mphi @ P8 @ Pnu                                              # positive multiple of Phi_0 (float-free integer arithmetic)
assert np.abs(Phi0).max() < 2 ** 62
ok_phi0 = not np.any(Phi0 @ (I512 - Gam)) and not np.any((I512 - Gam) @ Phi0)    # supported on Gamma = +1
F0 = Hinv * dm(R.T @ Phi0 @ R)
ok_phi0 &= dm(Phi0 @ R) == dm(R) * F0
# solve F0 = sum c_j F_j exactly from the trace pairing, then verify
rhs = sp.Matrix(60, 1, lambda i, _: tr(Fstar[i] * F0))
coef = GF.LUsolve(rhs)
Fsum = DomainMatrix.zeros((32, 32), QQ)
for c, F in zip(coef, Fs):
    c = sp.Rational(c)
    Fsum = Fsum + F * QQ(int(c.p), int(c.q))
ok_phi0 &= Fsum == F0
S0re = sum((sp.Rational(c) * sp.Matrix(b[0]) for c, b in zip(coef, basis60)), sp.zeros(8, 8))
S0im = sum((sp.Rational(c) * sp.Matrix(b[1]) for c, b in zip(coef, basis60)), sp.zeros(8, 8))
S0m = S0re + sp.I * S0im
ok_phi0 &= S0m.rank() == 1 and not S0m.is_zero_matrix
# its charges: B-L = -i HC/3 on the colour factor, T^R_3 = -i tau_3 on the weak factor; on S: h S + S h^T
HCm = -(gC[0] @ gC[3] + gC[1] @ gC[4] + gC[2] @ gC[5])[np.ix_(hC, hC)]
BLm = sp.Matrix(np.round(-1j * HCm / 3 * 3).astype(complex).tolist()) / 3
T3m = sp.Matrix(np.round(-1j * (-(gW[0] @ gW[1]) / 2)[np.ix_(hW, hW)] * 2).astype(complex).tolist()) / 2
BL8 = sp.kronecker_product(sp.nsimplify(BLm), sp.eye(2))
T38 = sp.kronecker_product(sp.eye(4), sp.nsimplify(T3m))
qBL = sp.simplify(BL8 * S0m + S0m * BL8.T)
qT3 = sp.simplify(T38 * S0m + S0m * T38.T)
ok_phi0 &= sp.simplify(qBL + 2 * S0m).is_zero_matrix and sp.simplify(qT3 - S0m).is_zero_matrix
# the spectrum: F* F is the realification of conj((eps (x) S)^dagger (eps (x) S)) = 1 (x) S S^dagger, checked at Phi_0
SSd = S0m * S0m.H
SSd_re = np.array(sp.re(SSd).tolist(), dtype=object); SSd_im = np.array(sp.im(SSd).tolist(), dtype=object)
dmS = lambda X: DomainMatrix([[QQ(int(sp.Rational(x).p), int(sp.Rational(x).q)) for x in row] for row in X.tolist()], X.shape, QQ)  # noqa: E731
lhs0 = (Hinv * F0.transpose() * H) * F0
rhs0 = Ureal_inv * realify(dmS(np.kron(np.eye(2, dtype=object), SSd_re)), dmS(np.kron(np.eye(2, dtype=object), SSd_im))) * Ureal
ok_spec = lhs0 == rhs0
print('   a13-2 parts:', ok_basis, ok_U, ok_iso, ok_ker, ok_norm, ok_phi0, ok_spec, flush=True)
report("a13-2", ok_basis and ok_U and ok_iso and ok_ker and ok_norm and ok_phi0 and ok_spec,
       f"prop:tensor, exact: an intertwiner U from V_+ = im P_*(1+Gamma) (J -> i) to C^2 (x) C^4 (x) C^2 for the six Lorentz, fifteen colour and "
       f"three su(2)_R generators (null space of dimension {null_dim}, rationalised and verified: {ok_U}), a scaled isometry, Re U^dagger U = {lam} H "
       f"({ok_iso}); for every S in T, v -> (eps (x) S) conj(v) is an antisymmetric Lorentz-invariant kernel ({ok_ker}), so iota(Phi) = S is a "
       f"Pati-Salam-equivariant isomorphism onto the (10,1,3) component, with Tr Phi^T Phi = {cst} |S|^2 on a basis ({ok_norm}); Phi_0 lies in it with "
       f"iota(Phi_0) of rank one, B-L = -2 and T^R_3 = 1, the tensor S_0 up to scale ({ok_phi0}); Phi^T Phi is the realification of 1 (x) S S^dagger, "
       f"so its spectrum is that of S^dagger S, each eigenvalue four times (checked at Phi_0: {ok_spec})")

sys.exit(0 if not FAILS else 1)
