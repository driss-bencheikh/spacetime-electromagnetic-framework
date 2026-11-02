"""Phase 9D, batch 4: adversarial audit of the Chapter 5-6 statements that had neither a cited check nor a written proof.

Integer realisations of Appendix A at n = 2 (64x64) and n = 3 (512x512). Ranks are computed modulo a large prime (a lower bound for the
rational rank) and, where a dimension is asserted, matched with explicit exact solutions; spectra of exact integer matrices are read by
float64 diagonalisation where stated (numerical). Polynomial identities are exact (sympy).
Each check prints PASS/FAIL a4-<n>; the script exits 0 iff all pass.
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
        print(f"PASS a4-{cid}: {text}")
    else:
        FAILS.append(cid)
        print(f"FAIL a4-{cid}: {text}")
    sys.stdout.flush()


PRIME = 2_147_483_629


def rank_mod_p(Mx, pr_=PRIME):
    Mx = np.array(Mx, dtype=np.int64) % pr_
    rows, cols = Mx.shape
    r = 0
    for c in range(cols):
        nz = np.nonzero(Mx[r:, c])[0]
        if nz.size == 0:
            continue
        piv = r + nz[0]
        Mx[[r, piv]] = Mx[[piv, r]]
        inv = pow(int(Mx[r, c]), pr_ - 2, pr_)
        Mx[r] = (Mx[r] * inv) % pr_
        others = np.nonzero(Mx[:, c])[0]
        others = others[others != r]
        if others.size:
            f = Mx[others, c].reshape(-1, 1)
            Mx[others] = (Mx[others] - (f * Mx[r]) % pr_) % pr_
        r += 1
        if r == rows:
            break
    return r


sx = np.array([[0, 1], [1, 0]]); sz = np.array([[1, 0], [0, -1]]); ep = np.array([[0, -1], [1, 0]]); I2 = np.eye(2, dtype=np.int64)
kr = lambda *a: reduce(np.kron, a)  # noqa: E731
e8 = [None, kr(sx, I2, I2), kr(sz, sx, I2), kr(sz, sz, sx), kr(ep, I2, I2), kr(sz, ep, I2), kr(sz, sz, ep)]
I8 = np.eye(8, dtype=np.int64)
om = reduce(np.matmul, e8[1:])
ETA = {a: (1 if a <= 3 else -1) for a in range(1, 7)}

# ---------------------------------------------------------------- a4-1: Chapter 5 tables (exact arithmetic)
n_ = sp.symbols("n", positive=True, integer=True)
ok = [sp.binomial(6 * n - 4, 2) for n in range(1, 6)] == [1, 28, 91, 190, 325]
census = {1: (2, 0), 2: (5, 3), 3: (8, 6), 4: (11, 9)}
ok &= all(census[n] == (3 * n - 1, 3 * n - 3) for n in census)
dims = [sp.binomial(p, 2) + sp.binomial(q, 2) for p, q in census.values()]
ranks = [p // 2 + q // 2 for p, q in census.values()]
ok &= dims == [1, 13, 43, 91] and ranks == [1, 3, 7, 9]
ok &= 2 + 1 + 1 == 4 and 8 + 3 + 1 == 12  # rank and dimension of su3+su2+u1
ok &= [sp.binomial(6 * n - 4, 2) - 12 * (n - 1) for n in range(1, 6)] == [1, 16, 67, 154, 277]
report("1", ok, "prop:frame, thm:census, lem:rank, cor:tradeoff: dim so(W_n) = C(6n-4,2) = 1,28,91,190,325; g_int = so(3n-1,3n-3); maximal compact "
       "dimensions 1,13,43,91 and ranks 1,3,7,9 against rank 4 and dimension 12 of su3+su2+u1; centraliser of e2e3 = C(6n-4,2) - 12(n-1) = "
       "1,16,67,154,277 (exact)")

# ---------------------------------------------------------------- a4-2: thm:mixeddirac, thm:coupled at n = 2 (64x64, exact rationals)
E2 = [None] + [np.kron(e8[a], I8) for a in range(1, 7)] + [np.kron(om, e8[a]) for a in range(1, 7)]
I64 = np.eye(64, dtype=np.int64)
S_ = lambda B: sp.Matrix(B)  # noqa: E731
cyc = [(1, 2, 3), (2, 3, 1), (3, 1, 2)]
g1 = [E2[1], E2[4], E2[5], E2[6]]
Sk = [g1[j] @ g1[k] for i, j, k in cyc]          # 2 S_k (rotations of the interior four-plane)
Tk = [g1[0] @ g1[a] for a in (1, 2, 3)]           # 2 T_k
Lk = [E2[6 + j + 3] @ E2[6 + k + 3] - E2[6 + j] @ E2[6 + k] for i, j, k in cyc]   # 2 L_k on constituent 2
Kk = [E2[6 + k] @ E2[6 + j + 3] - E2[6 + j] @ E2[6 + k + 3] for i, j, k in cyc]   # 2 K_k on constituent 2
Smix = [(Sk[i] + Lk[i]) for i in range(3)]
Tmix = [(Tk[i] + Kk[i]) for i in range(3)]
ok = True
# brackets of so(1,3) (in units of 2x generators: [2A,2B] = 2(2[A,B]))
cm = lambda A, B: A @ B - B @ A  # noqa: E731
for i, j, k in cyc:
    ok &= np.array_equal(cm(Smix[i - 1], Smix[j - 1]), 2 * Smix[k - 1]) or np.array_equal(cm(Smix[i - 1], Smix[j - 1]), -2 * Smix[k - 1])
C4 = sum(x @ x for x in Tmix) - sum(x @ x for x in Smix)  # = 4 C_mix
ok &= np.array_equal((C4 - 6 * I64) @ (C4 - 22 * I64), np.zeros((64, 64), dtype=np.int64))
P16 = (22 * I64 - C4)  # 16 P
ok &= np.array_equal(P16 @ P16, 16 * P16) and np.trace(P16) == 16 * 16
ok &= all(np.array_equal(P16 @ x, x @ P16) for x in Smix + Tmix)
J1 = E2[2] @ E2[3]
etaU = [1, -1, -1, -1]
gl = [etaU[m] * g1[m] for m in range(4)]
Pq = P16 / 16.0  # dyadic rationals: float64 arithmetic below is exact
gh = [2 * Pq @ g @ Pq for g in gl]
for mu in range(4):
    for nu in range(4):
        ok &= np.array_equal(gh[mu] @ gh[nu] + gh[nu] @ gh[mu], 2.0 * etaU[mu] * (1 if mu == nu else 0) * Pq)
ok &= all(np.array_equal(x @ J1, J1 @ x) for x in gh) and np.array_equal(J1 @ J1, -I64)
for (i, j, k), Sm, Tm, a in zip(cyc, Smix, Tmix, (1, 2, 3)):
    lhsS, rhsS = Pq @ (Sm / 2.0) @ Pq, gh[j] @ gh[k] / 2.0
    lhsT, rhsT = Pq @ (Tm / 2.0) @ Pq, gh[0] @ gh[a] / 2.0
    ok &= (np.array_equal(lhsS, rhsS) or np.array_equal(lhsS, -rhsS)) and (np.array_equal(lhsT, rhsT) or np.array_equal(lhsT, -rhsT))
# flavour: commutant of {ghat_mu, J1} on the block, its antisymmetric part
Pf = P16 / 16.0
w_, V_ = np.linalg.eigh(Pf)
Bk = V_[:, w_ > 0.5]
gens = [Bk.T @ g @ Bk for g in gh] + [Bk.T @ J1 @ Bk]
Acom = np.vstack([np.kron(np.eye(16), x) - np.kron(x.T, np.eye(16)) for x in gens])
s = np.linalg.svd(Acom, compute_uv=False)
dcom = int(np.sum(s < 1e-9))
ok &= dcom == 8
report("2", ok, "thm:mixeddirac (i)-(ii) and the flavour count of (iii) at n = 2: l_mix has the so(1,3) brackets; (C_mix - 3/2)(C_mix - 11/2) = 0; "
       "P = (11 - 2C_mix)/8 is the rank-16 projector commuting with the mixed generators; {ghat_mu, ghat_nu} = 2 eta_mu_nu P, [ghat, J1] = 0, "
       "J1^2 = -1; P S_k P and P T_k P are +-(1/2) ghat ghat (exact: integer and dyadic arithmetic); the commutant of {ghat, J1} on the block has real "
       f"dimension {dcom} = dim M_2(C), whose antisymmetric part is u(2) (numerical SVD)")

# thm:coupled: algebra of the real Lagrangian in flavour coordinates R^8 (x) R^2 (exact, sympy)
e8s = [None] + [sp.Matrix(x) for x in e8[1:]]
gam = [e8s[1], e8s[4], e8s[5], e8s[6]]
Jm = e8s[2] * e8s[3]
I2s = sp.eye(2)
gh_ = [-sp.kronecker_product(etaU[m] * gam[m], I2s) for m in range(4)]  # ghat_mu = -gamma_mu (x) 1 (lower index)
ghu = [etaU[m] * gh_[m] for m in range(4)]  # ghat^mu
J1s = sp.kronecker_product(Jm, I2s)
U = [ghu[0] * ghu[m] for m in range(4)]
ok = U[0] == sp.eye(16) and all(u.T == u for u in U) and all((J1s * u).T == -(J1s * u) for u in U) and all(J1s * u == u * J1s for u in U)
ok &= (ghu[0]).T == ghu[0]
psi = sp.Matrix(sp.symbols("y0:16"))
dpsi = [sp.Matrix(sp.symbols(f"d{m}_0:16")) for m in range(4)]
A = sp.symbols("A0:4")
q_, m_ = sp.symbols("q m")
Dpsi = [dpsi[mu] + q_ * A[mu] * J1s * psi for mu in range(4)]
Lf = sum(((psi.T * J1s * U[mu] * Dpsi[mu])[0] for mu in range(4))) - m_ * (psi.T * ghu[0] * psi)[0]
# Euler-Lagrange for psi with constant A (algebraic part): dL/dpsi - d_mu dL/d(d_mu psi), the second term = -sum_mu d_mu (U^mu^T J1^T psi) -> use symbolic derivative fields
dLdpsi = sp.Matrix([sp.diff(Lf, x) for x in psi])
dLdd = [sp.Matrix([sp.diff(Lf, x) for x in dpsi[mu]]) for mu in range(4)]
# d_mu of dL/d(d_mu psi) = (J1 U^mu)^T d_mu psi
EL = dLdpsi - sum(((J1s * U[mu]).T * dpsi[mu] for mu in range(4)), sp.zeros(16, 1))
target = 2 * (sum((J1s * U[mu] * Dpsi[mu] for mu in range(4)), sp.zeros(16, 1)) - m_ * ghu[0] * psi)
ok &= sp.expand(EL - target) == sp.zeros(16, 1)
# J1 U^mu D_mu psi = m ghat^0 psi  <=>  J1 ghat^mu D_mu psi = m psi (multiply by ghat^0)
ok &= all(sp.expand(ghu[0] * J1s * U[mu] - J1s * ghu[mu]) == sp.zeros(16) for mu in range(4))
# current: dL/dA_nu = q psi^T J1 U^nu J1 psi = -q psi^T U^nu psi ; j^0 = |psi|^2
ok &= all(sp.expand(sp.diff(Lf, A[nu]) + q_ * (psi.T * U[nu] * psi)[0]) == 0 for nu in range(4))
ok &= sp.expand((psi.T * U[0] * psi)[0] - (psi.T * psi)[0]) == 0
# exact invariance under the global rotation exp(theta J1): J1^T J1 U J1 = J1 U J1^T J1 ...
th = sp.symbols("theta", real=True)
R = sp.cos(th) * sp.eye(16) + sp.sin(th) * J1s
ok &= all(sp.simplify(R.T * J1s * U[mu] * R - J1s * U[mu]) == sp.zeros(16) for mu in range(4)) and sp.simplify(R.T * ghu[0] * R - ghu[0]) == sp.zeros(16)
report("3", ok, "thm:coupled (algebraic part, exact): U^0 = 1, U^mu symmetric, J1 U^mu antisymmetric and commuting; the Euler-Lagrange expression of eq:coupledL "
       "is 2(J1 U^mu D_mu psi - m ghat^0 psi), equivalent to J1 ghat^mu D_mu psi = m psi; dL/dA_nu = -q psi^T U^nu psi with j^0 = |psi|^2; the "
       "Lagrangian is invariant under exp(theta J1). The sign of q j in Maxwell's equation follows the convention dL/dA = -q j")

# ---------------------------------------------------------------- the n = 3 module
E = [None] + [kr(*([om] * k + [e8[a]] + [I8] * (2 - k))).astype(np.int64) for k in range(3) for a in range(1, 7)]
I512 = np.eye(512, dtype=np.int64)
Z512 = np.zeros((512, 512), dtype=np.int64)
blq = lambda S: reduce(np.matmul, [E[i] for i in S], I512)  # noqa: E731
eta18 = {A: (1 if (A - 1) % 6 < 3 else -1) for A in range(1, 19)}
interior = [1, 4, 5, 6]
internal = [i for i in range(1, 19) if i not in interior]
Jn = blq(internal)
Gam = Jn @ blq(interior)
J = {k: E[6 * k + 2] @ E[6 * k + 3] for k in range(3)}
omk = {k: blq(range(6 * k + 1, 6 * k + 7)) for k in range(3)}
X1 = J[0] @ J[1]
X2 = omk[2]
X3 = omk[0] @ (E[13] @ E[16])
Pr = {"Pfam": (I512 + X1) @ (I512 + X2) @ (I512 + X3), "P23": 2 * (I512 + X2) @ (I512 + X3), "P13": 2 * (I512 + X1) @ (I512 + X3),
      "P12": 2 * (I512 + X1) @ (I512 + X2), "P1": 4 * (I512 + X1)}  # all = 8 P
rankP = {k: np.trace(v) // 8 for k, v in Pr.items()}
ok = rankP == {"Pfam": 64, "P23": 128, "P13": 128, "P12": 128, "P1": 256} and all(np.array_equal(v @ v, 8 * v) for v in Pr.values())
ok &= all(np.array_equal(x @ y, y @ x) for x, y in itertools.combinations([X1, X2, X3], 2))
IB = list(itertools.combinations(internal, 2))
Bv = {ab: E[ab[0]] @ E[ab[1]] for ab in IB}
compact = {ab: eta18[ab[0]] == eta18[ab[1]] for ab in IB}
rng = np.random.default_rng(11)
probe = rng.integers(-3, 4, size=(512, 5))
gl3 = [E[1], E[4], E[5], E[6]]
table = {}


def gram(mats):
    """exact integer Gram matrix Tr(A^T B) of a list of integer matrices (float64 BLAS; every partial sum is an integer below 2^53)."""
    flat = np.array([m.reshape(-1) for m in mats], dtype=np.float64)
    g = flat @ flat.T
    assert np.abs(g).max() < 2 ** 52 and np.array_equal(g, np.round(g))
    return sp.Matrix(np.round(g).astype(np.int64))


Bf = {ab: Bv[ab].astype(np.float64) for ab in IB}


for name, P8i in Pr.items():
    P8 = P8i.astype(np.float64)
    res = []
    for part in ([ab for ab in IB if compact[ab]], [ab for ab in IB if not compact[ab]]):
        Gc = gram([Bf[ab] @ P8 - P8 @ Bf[ab] for ab in part])          # rank of X -> [X,P] on the span of this part
        Nsp = Gc.nullspace()                                          # the commuting subspace (exact rational coefficients)
        if not Nsp:
            res.append((0, 0))
            continue
        Nm = sp.Matrix.hstack(*Nsp)
        H = gram([Bf[ab] @ P8 for ab in part])                        # Tr(P B_k^T B_l P) up to the factor 64
        img = (Nm.T * H * Nm).rank()                                  # dimension of the image on the block
        res.append((Nm.shape[1], img))
    (ck, ik), (cp, ip) = res
    axial = all(np.array_equal(P8i @ g @ P8i, Z512) for g in gl3)
    table[name] = (rankP[name], ck + cp, ck, cp, (ck - ik) + (cp - ip), axial, ik)
expected = {"Pfam": (64, 17, 9, 8, 3, True), "P23": (128, 23, 9, 14, 0, True), "P13": (128, 33, 17, 16, 3, True), "P12": (128, 27, 15, 12, 3, False),
            "P1": (256, 51, 27, 24, 3, False)}
keff_dim = {"Pfam": 3 + 3, "P23": 3 + 3 + 3, "P13": 1 + 3 + 10, "P12": 4 * 3, "P1": 15 + 3 * 3}  # u1^3+su2, u1^3+su2+su2, u1+su2+so5, su2^4, so6+su2^3
ok &= all(table[k][:6] == expected[k] and table[k][6] == keff_dim[k] for k in expected)
report("4", ok, "thm:table: for the five projectors (rank, dim g_phys, dim k, dim p, dim ker, axial) = "
       f"{ {k: v[:6] for k, v in table.items()} }; the commuting space is the exact nullspace of X -> [X,P] on the span of the 91 internal bivectors "
       "(Frobenius Gram matrices over the integers, sympy ranks), split into its compact and non-compact parts; the columns g_phys, k, p count "
       "the commuting algebra before the quotient by its kernel ker, and the image of k on the block has dimension k - ker = dim k_eff = "
       f"{ {k: v[6] for k, v in table.items()} } as listed (exact)")

ok = np.array_equal((J[0] + J[1]) @ (I512 + X1), Z512) and not np.array_equal((J[0] + J[1]) @ Pr["P23"], Z512)
ok &= np.array_equal(J[0] @ J[0], -I512) and np.array_equal(J[0] @ J[1], J[1] @ J[0])
f = [10, 11, 12, 16, 17, 18]
fb = [E[a] @ E[b] for a, b in itertools.combinations(f, 2)]
P1 = Pr["P1"]
ok &= all(np.array_equal(x @ x, -I512) and np.array_equal(x.T, -x) and np.array_equal(x @ X1, X1 @ x) for x in fb)
ok &= rank_mod_p(np.array([((P1 @ x @ P1) @ probe).reshape(-1) for x in fb]).T) == 15
Om6 = E[10] @ E[16] + E[11] @ E[17] + E[12] @ E[18]
cols = [(x @ Om6 - Om6 @ x).reshape(-1)[::3] for x in fb]
cdim = 15 - rank_mod_p(np.array(cols).T)
cent = [x for x in fb if np.array_equal(x @ Om6, Om6 @ x)]
# u(3): the centraliser contains the 3 Cartan-like f_k f_{k+3}'s and combinations; build it exactly as the nullspace via sympy on the 15x15 adjoint
adj = sp.Matrix([[sp.Rational(int(np.sum((fb[j] @ Om6 - Om6 @ fb[j]) * fb[i])), -512) for j in range(15)] for i in range(15)])
nsp = adj.nullspace()
ok &= cdim == 9 and len(nsp) == 9
centM = [sum(int(v[j]) * fb[j] for j in range(15)) if all(x.is_integer for x in v) else None for v in [vv * sp.lcm([sp.fraction(x)[1] for x in vv]) for vv in nsp]]
dercols = [((a @ b - b @ a) @ probe).reshape(-1) for a, b in itertools.combinations(centM, 2)]
der = rank_mod_p(np.array(dercols).T)
ok &= der == 8
ok &= rank_mod_p(np.array([((P1 @ (a @ b - b @ a) @ P1) @ probe).reshape(-1) for a, b in itertools.combinations(centM, 2)]).T) == 8
ok &= all(np.array_equal(X1 @ g, g @ X1) for g in gl3)
report("5", ok, "ex:annihilator, thm:su3, thm:p1: (J1+J2)(1+X1) = 0 but (J1+J2)P23 != 0; the 15 bivectors f_a f_b square to -1, are antisymmetric, commute "
       "with X1 and act faithfully on the P1 block (rank 15); the centraliser of Omega_6 in so(6) has dimension 9 and derived dimension 8, faithful on "
       "the block; X1 commutes with every interior gamma (exact integer arithmetic; ranks mod p)")

# thm:PS, cor:PSfamily, prop:PSledger, prop:PSanom
chiW, chiC, qA = blq([7, 13, 14, 15]), blq([10, 11, 12, 16, 17, 18]), E[2] @ E[3]
Xs, Zs, Us = Gam @ chiW, -(Jn @ chiC), -(Jn @ qA)
ok = all(np.array_equal(x @ x, I512) and np.array_equal(x, x.T) and np.trace(x) == 0 for x in (Xs, Zs, Us))
ok &= all(np.array_equal(x @ y, y @ x) for x, y in itertools.combinations([Xs, Zs, Us], 2))
ok &= all(np.trace(reduce(np.matmul, c)) == 0 for r in (2, 3) for c in itertools.combinations([Xs, Zs, Us], r))
Ps8 = (I512 + Xs) @ (I512 + Zs) @ (I512 + Us)
ok &= np.array_equal(Ps8 @ Ps8, 8 * Ps8) and np.trace(Ps8) == 512
Sig = [gl3[m] @ gl3[n] for m, n in itertools.combinations(range(4), 2)]
ok &= all(np.array_equal(Ps8 @ s_, s_ @ Ps8) for s_ in Sig) and np.array_equal(Ps8 @ Jn, Jn @ Ps8)
ok &= all(np.array_equal(Ps8 @ g @ Ps8, Z512) for g in gl3)
ok &= np.trace(Ps8 @ (I512 + Gam)) == 2 * 8 * 32 and np.trace(Ps8 @ (I512 - Gam)) == 2 * 8 * 32
ok &= np.array_equal(chiW @ Ps8, Gam @ Ps8)
ok &= np.array_equal(chiC @ chiC, -I512)
GC = Jn @ chiC
ok &= np.array_equal(GC @ GC, I512) and np.array_equal(GC @ Ps8, -Ps8) and np.array_equal(chiC @ Ps8, Jn @ Ps8)
Rm = Ps8 @ E[1] @ E[7] @ E[8]  # 8 R
ok &= np.array_equal(Rm @ Rm, -8 * Ps8)  # (8R)^2 = -64 P_* = -8 (8 P_*)
colour = [E[a] @ E[b] for a, b in itertools.combinations([10, 11, 12, 16, 17, 18], 2)]
ok &= all(np.array_equal(Rm @ c, c @ Rm) for c in colour) and np.array_equal(Rm @ Gam, -Gam @ Rm)
Hc = -(E[10] @ E[16] + E[11] @ E[17] + E[12] @ E[18])
Qc = -(Jn @ Hc)
ok &= np.array_equal((Qc - I512) @ (Qc + 3 * I512) @ Ps8, Z512)
mult1 = rank_mod_p(((Qc + 3 * I512) @ Ps8) @ rng.integers(-3, 4, size=(512, 64)))
mult3 = rank_mod_p(((Qc - I512) @ Ps8) @ rng.integers(-3, 4, size=(512, 64)))
ok &= mult1 == 48 and mult3 == 16
report("6", ok, "thm:PS, cor:PSfamily, prop:PSledger: X, Z, U commuting symmetric traceless involutions with traceless products; P_* projector of trace 64 "
       "commuting with Sigma^mu nu and J, P_* gamma P_* = 0, Weyl blocks of rank 32, chi_W P_* = Gamma P_*; chi_C^2 = -1, Gamma_C^2 = 1, Gamma_C P_* = "
       "-P_*; R = P_* gamma^0 e^(2)_1 e^(2)_2 has R^2 = -P_*, commutes with colour and exchanges the Weyl blocks; (Q_c - 1)(Q_c + 3)P_* = 0 with "
       f"multiplicities {mult1} (+1) and {mult3} (-3) (exact integer arithmetic)")

# content per chirality and anomalies (numerical spectra of exact integer matrices; anomaly algebra exact)
W = [7, 13, 14, 15]
Rw = {(a, b): E[W[a]] @ E[W[b]] / 2 for a, b in itertools.combinations(range(4), 2)}
Rf = lambda a, b: Rw[(a, b)] if a < b else -Rw[(b, a)]  # noqa: E731
Xp = [(Rf(1, 2) + Rf(0, 3)) / 2, (Rf(2, 0) + Rf(1, 3)) / 2, (Rf(0, 1) + Rf(2, 3)) / 2]
Xm = [(Rf(1, 2) - Rf(0, 3)) / 2, (Rf(2, 0) - Rf(1, 3)) / 2, (Rf(0, 1) - Rf(2, 3)) / 2]
CL, CR = -sum(x @ x for x in Xp), -sum(x @ x for x in Xm)
Ccol = -sum((c / 2) @ (c / 2) for c in colour)
w_, V_ = np.linalg.eigh(Ps8 / 8.0)
Bk = V_[:, w_ > 0.5]
Mix = Bk.T @ (CL + np.e * CR + np.pi * Gam + 0.37 * Ccol) @ Bk
_, Vv = np.linalg.eigh(Mix)
basis = Bk @ Vv
vals = [np.einsum("ij,ij->j", basis, O @ basis) for O in (CL, CR, Gam, Ccol)]
content = {}
ok = True
for a, b, g, c in zip(*vals):
    key = (round(4 * a) / 4, round(4 * b) / 4, int(round(g)), round(4 * c) / 4)
    ok &= abs(a - key[0]) < 1e-9 and abs(b - key[1]) < 1e-9 and abs(c - key[3]) < 1e-9
    content[key] = content.get(key, 0) + 1
ok &= content == {(0.75, 0.0, -1, 3.75): 32, (0.0, 0.75, 1, 3.75): 32} or content == {(0.75, 0.0, 1, 3.75): 32, (0.0, 0.75, -1, 3.75): 32}
Lsign = [k[2] for k in content if k[0] == 0.75][0]
a_, b_ = sp.symbols("a b")
qL, qR = a_ + b_ * Lsign, a_ - b_ * Lsign  # charge a + b Gamma on each half
qRc = -qR  # right-handed fields written left-handed
nW = 8  # Weyl fermions per half: (4,2) and (4,2)
anom = {"SU2L^2": sp.expand(4 * sp.Rational(1, 2) * qL), "SU2R^2": sp.expand(4 * sp.Rational(1, 2) * qRc), "SU4^2": sp.expand(2 * sp.Rational(1, 2) * (qL + qRc)),
        "cubic": sp.expand(nW * (qL ** 3 + qRc ** 3)), "grav": sp.expand(nW * (qL + qRc))}
exp_anom = {"SU2L^2": 2 * (a_ - b_), "SU2R^2": -2 * (a_ + b_), "SU4^2": -2 * b_, "cubic": sp.expand(-16 * (3 * a_ ** 2 * b_ + b_ ** 3)), "grav": -16 * b_}
ok &= Lsign == -1 and all(sp.expand(anom[k] - exp_anom[k]) == 0 for k in anom)
ok &= sp.solve([anom["grav"], anom["SU2L^2"]], [a_, b_], dict=True) == [{a_: 0, b_: 0}]
report("7", ok, "cor:PSfamily, prop:PSanom: on the P_* block the Gamma = -1 half carries su2_L Casimir 3/4, su2_R 0, colour Casimir 15/4 (the 4 of su4) and the "
       "Gamma = +1 half the reverse: (4,2,1)_L + (4,1,2)_R (numerical spectra of exact integer matrices, 1e-9); for the central charge a + b Gamma the "
       "anomaly coefficients follow exactly: 2(a-b), -2(a+b), -2b, -16(3a^2 b + b^3), -16b, vanishing only at a = b = 0")

# prop:orbits: stabilisers of P_fam and P_* in the compact connection algebra so8 + so6 (43 compact internal bivectors)
comp = [ab for ab in IB if compact[ab]]
ok = len(comp) == 43
stab = {}
for name, P8i in (("Pfam", Pr["Pfam"]), ("Pst", Ps8)):
    P8 = P8i.astype(np.float64)
    Gc = gram([Bf[ab] @ P8 - P8 @ Bf[ab] for ab in comp])
    stab[name] = len(comp) - Gc.rank()
orbit = {k: 43 - v for k, v in stab.items()}
report("8", stab == {"Pfam": 9, "Pst": 23}, f"prop:orbits: the stabilisers of P_fam and P_* in so8+so6 (the exact commuting compact subspaces, annihilators "
       f"included) have dimensions {stab['Pfam']} and {stab['Pst']}, so the orbits have dimensions {orbit['Pfam']} and {orbit['Pst']} (exact)")
ORBIT = orbit

# ---------------------------------------------------------------- a4-9: thm:variety, thm:nonu, prop:integrality, ex:multiset (exact)
Nc, q, l, yu, yd, yn, ye, h = sp.symbols("N_c q l y_u y_d y_nu y_e h")
# left-handed: Q_L (Nc colours, doublet) hypercharge q, L_L doublet l; right-handed u, d (Nc each), nu, e
lin1 = 2 * Nc * q + 2 * l - Nc * (yu + yd) - (yn + ye)  # gravitational
lin2 = Nc * q + l  # SU(2)^2 Y   (doublets)
lin3 = 2 * q - yu - yd  # SU(3)^2 Y
solL = sp.solve([lin2, lin3, lin1], [l, yd, ye], dict=True)[0]
ok = solL == {l: -Nc * q, yd: 2 * q - yu, ye: 2 * (-Nc * q) - yn}
cub = 2 * Nc * q ** 3 + 2 * l ** 3 - Nc * (yu ** 3 + yd ** 3) - (yn ** 3 + ye ** 3)
cubs = sp.expand(cub.subs(solL))
hq, hl = sp.symbols("h_q h_l")
cubs2 = sp.expand(cubs.subs({yu: q + hq, yn: -Nc * q + hl}))
ok &= sp.expand(cubs2 - 6 * Nc * q * (hl ** 2 - hq ** 2)) == 0
ok &= sp.expand(2 * q ** 3 - yu ** 3 - (2 * q - yu) ** 3 + 6 * q * (yu - q) ** 2) == 0
# without nu_R
cub_nonu = sp.expand((2 * Nc * q ** 3 + 2 * l ** 3 - Nc * (yu ** 3 + yd ** 3) - ye ** 3).subs({l: -Nc * q, yd: 2 * q - yu, ye: -2 * Nc * q}).subs(yu, q + hq))
ok &= sp.expand(cub_nonu - 6 * Nc * q * (Nc ** 2 * q ** 2 - hq ** 2)) == 0
sm = {q: sp.Rational(1, 6), Nc: 3}
ok &= sp.solve(cub_nonu.subs(sm), hq) == [-sp.Rational(1, 2), sp.Rational(1, 2)]
lam = sp.symbols("lambda")
Qnu, Qe, Qu = -lam, -1 - lam, sp.Rational(2, 3) + lam / 3
Qd = Qu - 1
ok &= sp.expand(3 * Qd - (lam - 1)) == 0 and sp.expand(3 * Qu - (lam + 2)) == 0
ok &= sp.solve(Qnu * Qe, lam) == [-1, 0]
ms = [-3] * 5 + [1] * 10 + [5]
ok &= sum(ms) == 0 and sum(x ** 3 for x in ms) == 0 and sorted(ms) != sorted(-x for x in ms)
signs = [((1, 1, 1), 1), ((-1, -1, 1), 10), ((1, -1, -1), 5)]
tr = [sum(s[i] * mlt for s, mlt in signs) for i in range(3)]
Y = [s[0] + s[1] + 3 * s[2] for s, _ in signs]
ok &= tr == [-4, -14, 6] and sorted(sum(([y] * mlt for y, (_, mlt) in zip(Y, signs)), [])) == sorted(ms)
report("9", ok, "thm:variety, thm:nonu, prop:integrality, ex:multiset (exact): the linear conditions give l = -N_c q, y_d = 2q - y_u, y_e = 2l - y_nu; the "
       "cubic is 6 N_c q (h_l^2 - h_q^2); without nu_R it is 6 N_c q (N_c^2 q^2 - h_q^2), at N_c = 3, q = 1/6 h_q = +-1/2 (the SM values); "
       "Q(ddd) = lambda - 1, Q(uuu) = lambda + 2; a neutral elementary lepton iff lambda in {0,-1}; the multiset and its three involutions "
       "with traces (-4,-14,6) (exact)")

print(f"{NPASS[0]}/{NPASS[0] + len(FAILS)} checks passed")
sys.exit(0 if not FAILS else 1)
