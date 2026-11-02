#!/usr/bin/env python3
"""v1: the tensor model T = Sym^2(C^4) (x) Sym^2(C^2) of the (10,1,3) component (Phase 9B, 2 November 2026).

    python3 verification/vacuum/v1_tensor_model.py

Rebuilds, from the definitions of Section 7.3 of the manuscript (subsections "The tensor model" to "The fluctuations decide"), the
representation of su(4)+su(2)_R on T, the projectors P_1..P_5 of Sym^2 T, the invariants I_k, the holomorphic quartic D, the tree potential
(eq:Vgeneral) and the one-loop vector and fermion traces, and tests: thm:potential(i)-(ii), the recoupling identity, lem:hw, thm:select0 (orbit
and stabiliser), the benchmark statements after thm:selectk, thm:secondvar and eq:sumrules, thm:oneloop1013 and the masses quoted after it.

Method: floating-point linear algebra (numpy) on exactly specified matrices; every comparison prints its tolerance. Random inputs come from
numpy.random.default_rng with the seeds printed. One PASS/FAIL line per check; exit 0 iff all pass. A PASS shows that the coded assertion
holds; it does not by itself prove the theorem around it (see the author's verification baseline).
"""
import itertools
import sys

import numpy as np

NPASS, FAILS = [0], []


def report(cid, ok, text):
    if ok:
        NPASS[0] += 1
        print(f"PASS v1-{cid}: {text}")
    else:
        FAILS.append(cid)
        print(f"FAIL v1-{cid}: {text}")
    sys.stdout.flush()


rng = np.random.default_rng(20261102)
TOL = 1e-10

# ------------------------------------------------------------------ generators, Tr(T_a T_b) = 1/2 delta_ab in the fundamental
def su(n):
    gens = []
    for a, b in itertools.combinations(range(n), 2):
        X = np.zeros((n, n), complex); X[a, b] = X[b, a] = 0.5; gens.append(X)
        Y = np.zeros((n, n), complex); Y[a, b] = -0.5j; Y[b, a] = 0.5j; gens.append(Y)
    for k in range(1, n):
        H = np.zeros((n, n), complex)
        H[:k, :k] = np.eye(k); H[k, k] = -k
        gens.append(H / np.sqrt(2 * k * (k + 1)))
    return gens


g4, g2 = su(4), su(2)
I4c, I2c = np.eye(4), np.eye(2)
G8 = [np.kron(X, I2c) for X in g4] + [np.kron(I4c, Y) for Y in g2]       # 15 + 3 Hermitian generators on C^8 = C^4 (x) C^2
assert all(abs(np.trace(X @ X) - 0.5) < 1e-14 for X in g4 + g2)

# ------------------------------------------------------------------ orthonormal basis of T (complex dim 30), Tr S^dagger S' = delta
def sym_basis(n):
    B = []
    for a in range(n):
        M = np.zeros((n, n)); M[a, a] = 1; B.append(M)
    for a, b in itertools.combinations(range(n), 2):
        M = np.zeros((n, n)); M[a, b] = M[b, a] = 1 / np.sqrt(2); B.append(M)
    return B


FB = [np.kron(A, B) for A in sym_basis(4) for B in sym_basis(2)]          # 30 matrices 8x8, index (a,s) -> 2a+s, s=0 is e_+
NT = len(FB)
FBm = np.array(FB)                                                          # 30 x 8 x 8


def to_mat(s):
    return np.einsum("j,jab->ab", s, FBm)


def to_vec(S):
    return np.einsum("jab,ab->j", FBm.conj(), S)


def rho(X):  # representation matrix on T for a Hermitian generator X of su4+su2: S -> X S + S X^T
    return np.array([[np.trace(F1.conj().T @ (X @ F2 + F2 @ X.T)) for F2 in FB] for F1 in FB])


R = [rho(X) for X in G8]
R4, R2 = R[:15], R[15:]
ok = NT == 30 and all(np.allclose(r, r.conj().T, atol=1e-14) for r in R)
ok &= all(np.allclose(to_vec(F), np.eye(NT)[j]) for j, F in enumerate(FB))

# ------------------------------------------------------------------ Sym^2 T: the five components by their Casimirs
I30 = np.eye(NT)
def tens(r):
    return np.kron(r, I30) + np.kron(I30, r)


C4 = sum(tens(r) @ tens(r) for r in R4)
C2 = sum(tens(r) @ tens(r) for r in R2)
SW = np.zeros((NT * NT, NT * NT))
for i in range(NT):
    for j in range(NT):
        SW[i * NT + j, j * NT + i] = 1
Psym = (np.eye(NT * NT) + SW) / 2
Hc = C4 + np.pi * C2
ev, U = np.linalg.eigh(Psym @ Hc @ Psym + 1000 * (np.eye(NT * NT) - Psym))
groups = []
for val in np.unique(np.round(ev, 7)):
    if val > 999:
        continue
    cols = U[:, np.abs(ev - val) < 1e-6]
    groups.append((cols.shape[1], float(np.real(np.trace(cols.conj().T @ C4 @ cols)) / cols.shape[1]),
                   float(np.real(np.trace(cols.conj().T @ C2 @ cols)) / cols.shape[1]), cols))
dims = sorted(g[0] for g in groups)
order = {175: 0, 35: 1, 100: 2, 20: 3, 135: 4}                           # (35,5), (35,1), (20',5), (20',1), (45,3)
ok &= dims == sorted(order)
P = [None] * 5
labels = [None] * 5
for d, c4, c2, cols in groups:
    P[order[d]] = cols @ cols.conj().T
    labels[order[d]] = (round(c4, 6), round(c2, 6))
ok &= abs(sum(np.trace(p).real for p in P) - 465) < 1e-8
# su(2) Casimir of spin j with Tr = 1/2 normalisation is j(j+1): 5 -> 6, 1 -> 0, 3 -> 2; su(4) labels must separate 35, 20', 45
ok &= [l[1] for l in labels] == [6.0, 0.0, 6.0, 0.0, 2.0] and labels[0][0] == labels[1][0] and labels[2][0] == labels[3][0]
report("1", ok, f"T = Sym^2(C^4)(x)Sym^2(C^2) has complex dimension 30 and Sym^2 T = (35,5)+(35,1)+(20',5)+(20',1)+(45,3) with dimensions "
       f"175+35+100+20+135 = 465, separated by the joint Casimir values (su4, su2) = {labels} (float, tol 1e-6)")


def Ik(s):
    w = np.kron(s, s)
    return np.array([np.real(w.conj() @ (p @ w)) for p in P])


def rand_s():
    return rng.normal(size=NT) + 1j * rng.normal(size=NT)


S0 = np.zeros((8, 8), complex); S0[6, 6] = 1.0                            # u_0 = e_4 (x) e_+ has index 2*3+0 = 6
s0 = to_vec(S0)
one = to_vec(np.eye(8).astype(complex))
ok = np.allclose(to_mat(s0), S0) and np.allclose(to_mat(one), np.eye(8))

# ------------------------------------------------------------------ recoupling identity and lem:hw
a_rec = np.array([1, -0.5, -0.5, 0.25, 0])
err = 0.0
for _ in range(20):
    s = rand_s(); S = to_mat(s)
    lhs = np.real(np.trace((S.conj().T @ S) @ (S.conj().T @ S)))
    err = max(err, abs(lhs - a_rec @ Ik(s)) / np.linalg.norm(s) ** 4, abs(Ik(s).sum() - np.linalg.norm(s) ** 4) / np.linalg.norm(s) ** 4)
ok &= err < TOL
report("2", ok, f"recoupling identity Tr(S^dag S)^2 = I1 - I2/2 - I3/2 + I4/4 and sum_k I_k = |S|^4 at 20 random S (seed 20261102), "
       f"max relative error {err:.1e} (tol {TOL:g})")

f0 = Ik(s0)
fr = []
for _ in range(20):
    a = rng.normal(size=4) + 1j * rng.normal(size=4); b = rng.normal(size=2) + 1j * rng.normal(size=2)
    u = np.kron(a, b); fr.append(Ik(to_vec(np.outer(u, u))) / np.linalg.norm(np.outer(u, u)) ** 4)
gen = [Ik(s) / np.linalg.norm(s) ** 4 for s in (rand_s() for _ in range(20))]
ok = np.allclose(f0, [1, 0, 0, 0, 0], atol=TOL) and all(np.allclose(f, [1, 0, 0, 0, 0], atol=1e-9) for f in fr)
ok &= max(f[0] for f in gen) < 1 - 1e-3
report("3", ok, f"lem:hw: f(S_0) = {np.round(f0, 12).tolist()} and f = (1,0,0,0,0) on 20 random rank-one S = (a(x)b)(a(x)b)^T; at 20 random S "
       f"f_1 <= {max(f[0] for f in gen):.3f} < 1 (float)")

# ------------------------------------------------------------------ thm:select0: orbit dimension and stabiliser of S_0
def real_cols(vs):  # complex vectors -> real matrix (Re; Im)
    return np.vstack([np.real(vs), np.imag(vs)])


tang = real_cols(np.array([r @ s0 for r in R]).T)                        # 60 x 18
rk = np.linalg.matrix_rank(tang, tol=1e-9)
null = np.linalg.svd(tang)[2][rk:]                                         # stabiliser coefficients (real combinations of the 18)
stab = [sum(c * X for c, X in zip(v, G8)) for v in null]
# su(3) on e_1..e_3 and Y = T3R + (B-L)/2 must lie in the stabiliser
BmL = np.kron(np.diag([1 / 3, 1 / 3, 1 / 3, -1]), I2c); T3R = np.kron(I4c, np.diag([0.5, -0.5]))
Yop = T3R + BmL / 2
cands = [np.kron(np.pad(X, ((0, 1), (0, 1))), I2c) for X in su(3)] + [Yop]
span = np.array([np.concatenate([np.real(X).ravel(), np.imag(X).ravel()]) for X in stab]).T
inside = all(np.linalg.norm(span @ np.linalg.lstsq(span, np.concatenate([np.real(X).ravel(), np.imag(X).ravel()]), rcond=None)[0]
                            - np.concatenate([np.real(X).ravel(), np.imag(X).ravel()])) < 1e-9 for X in cands)
ok = rk == 9 and len(stab) == 9 and inside
report("4", ok, "thm:select0: the orbit of S_0 under su(4)+su(2)_R has real dimension 9 (rank of the 18 tangent vectors) and the stabiliser "
       "has dimension 9 and contains su(3) and Y = T^R_3 + (B-L)/2, so it is su(3)+u(1)_Y (float, tol 1e-9)")

# ------------------------------------------------------------------ the holomorphic quartic D
eps2 = np.array([[0, 1], [-1, 0]])
eps4 = np.zeros((4,) * 4)
for p in itertools.permutations(range(4)):
    eps4[p] = np.linalg.det(np.eye(4)[list(p)])


def Dq(S):
    S4 = S.reshape(4, 2, 4, 2)
    Pt = np.einsum("su,tv,aset,bufv->abef", eps2, eps2, S4, S4)
    return np.einsum("abcd,efgh,abef,cdgh->", eps4, eps4, Pt, Pt) / 24


errs = []
for _ in range(5):
    s = rand_s(); S = to_mat(s)
    th = rng.normal()
    errs.append(abs(Dq(np.exp(1j * th) * S) - np.exp(4j * th) * Dq(S)) / np.linalg.norm(s) ** 4)
    Xh = sum(rng.normal() * X for X in G8)
    Ug = np.linalg.eigh(Xh)[1] @ np.diag(np.exp(1j * np.linalg.eigvalsh(Xh))) @ np.linalg.eigh(Xh)[1].conj().T
    errs.append(abs(Dq(Ug @ S @ Ug.T) - Dq(S)) / np.linalg.norm(s) ** 4)
ok = max(errs) < TOL
D1 = Dq(np.eye(8).astype(complex))
I1n = Ik(one)
ok &= abs(abs(D1) - 4) < TOL and abs(I1n[3] - 32 / 3) < 1e-9
report("5", ok, f"thm:potential(i): D(e^(i th) S) = e^(4 i th) D(S) and D(U S U^T) = D(S) for U in SU(4)xSU(2)_R (max rel. error {max(errs):.1e}); "
       f"|D(1(x)1)| = {abs(D1):.12f} = (3/8) I_4(1(x)1) with I_4(1(x)1) = {I1n[3]:.12f} = 32/3 (float, tol {TOL:g})")

# D depends on S only through the (20',1) component, and |D| <= (3/8) I_4 by Schur's lemma
Wb = np.linalg.eigh(P[3])[1][:, -20:]                                       # orthonormal basis of the (20',1) component in T (x) T
RW = [Wb.conj().T @ tens(r) @ Wb for r in R]                                 # 20 x 20 action
nb = 20
idx = [(i, j) for i in range(nb) for j in range(i, nb)]
rows = []
for r in RW:   # invariance of a complex-bilinear symmetric beta: r^T beta + beta r = 0
    blk = np.zeros((nb * nb, len(idx)), complex)
    for k, (i, j) in enumerate(idx):
        Bm = np.zeros((nb, nb)); Bm[i, j] = Bm[j, i] = 1
        blk[:, k] = (r.T @ Bm + Bm @ r).ravel()
    rows.append(blk)
A_inv = np.vstack(rows)
sv = np.linalg.svd(A_inv, compute_uv=False)
nul = int(np.sum(sv < 1e-9))
vb = np.linalg.svd(A_inv)[2][-1].conj()
beta = np.zeros((nb, nb), complex)
for k, (i, j) in enumerate(idx):
    beta[i, j] = beta[j, i] = vb[k]
ratios = []
for _ in range(8):
    s = rand_s(); om = Wb.conj().T @ np.kron(s, s)
    ratios.append(Dq(to_mat(s)) / (om @ beta @ om))
cst = np.mean(ratios)
sig = np.linalg.svd(beta, compute_uv=False)
bound = abs(cst) * sig.max()
ok = nul == 1 and np.allclose(ratios, cst, rtol=1e-9) and np.allclose(sig, sig[0], rtol=1e-9) and abs(bound - 3 / 8) < 1e-9
report("6", ok, f"thm:potential(ii): the (20',1) component carries a one-dimensional space of invariant symmetric bilinear forms beta (null space "
       f"of the invariance equations: {nul}); D(S) = c*beta(w,w) with w = P_4(S(x)S) at 8 random S (spread {np.ptp(np.abs(ratios)) / abs(cst):.1e}); "
       f"all 20 singular values of beta are equal (Schur), so |D(S)| <= |c| sigma I_4(S) with |c| sigma = {bound:.12f} = 3/8 (float, tol 1e-9). "
       "With v1-5 the bound is attained at 1(x)1. This proves the inequality given the numerically identified constant")

# D, its first and second derivatives vanish at S_0
co = []
for _ in range(5):
    dl = to_mat(rand_s())
    ts = np.array([-2, -1, 0.5, 1, 2])
    vals = np.array([Dq(S0 + t * dl) for t in ts])
    co.append(np.abs(np.linalg.solve(np.vander(ts, 5, increasing=True), vals)[:3]).max() / np.linalg.norm(dl) ** 4)
ok = max(co) < 1e-9
report("7", ok, f"D and its first and second derivatives vanish at S_0: the coefficients of t^0, t^1, t^2 in D(S_0 + t delta) are at most "
       f"{max(co):.1e} for 5 random delta (float, tol 1e-9); hence the second variation at S_0 does not depend on kappa (thm:selectk(d))")

# ------------------------------------------------------------------ benchmark after thm:selectk
f1 = I1n / np.linalg.norm(one) ** 4
d1 = abs(D1) / np.linalg.norm(one) ** 4
lam = np.array([1, 1.7, 2.3, 3.1, 2.9])
Dl = lam[1:] - lam[0]
gap = Dl @ f1[1:]
ok = np.allclose(f1, [1 / 3, 1 / 6, 1 / 3, 1 / 6, 0], atol=TOL) and abs(d1 - 1 / 16) < TOL and abs(gap - 0.9) < TOL
ok &= abs(4 / 3 * Dl[2] - 2.8) < 1e-12 and abs(4 / 3 * (Dl[0] + 2 * Dl[1] + Dl[2]) - 7.2) < 1e-12 and abs(4 / 3 * lam[3] - 4.1333333333333) < 1e-9
ok &= abs(4 / 3 * (2 * lam[0] + lam[1] + 2 * lam[2] + lam[3]) - 15.2) < 1e-12
report("8", ok, f"benchmark lambda = (1,1.7,2.3,3.1,2.9): f(1(x)1) = (1/3,1/6,1/3,1/6,0), d(1(x)1) = 1/16, so Lambda_kappa(1(x)1) - Lambda_kappa(S_0) = "
       f"{gap:.12f} - |kappa|/8, zero at |kappa| = 7.2; the quoted bounds 2.8, 7.2, 4.13, 15.2 follow from the stated formulas (float, tol {TOL:g})")

# ------------------------------------------------------------------ Hessians on the 60 real coordinates
Dmat = np.hstack([np.eye(NT), 1j * np.eye(NT)])                            # real coordinates (x, y) -> complex s = x + i y


def hess_I(s):
    """Hessians of I_1..I_5 at s in the 60 real coordinates."""
    W1 = np.kron(s[:, None], Dmat) + np.kron(Dmat, s[:, None])            # 900 x 60: s(x)delta + delta(x)s
    out = []
    w = np.kron(s, s)
    for p in P:
        h1 = 2 * np.real(W1.conj().T @ p @ W1)
        Cm = np.conj((p @ w).reshape(NT, NT))
        h2 = 2 * np.real(Dmat.T @ (Cm + Cm.T) @ Dmat)
        out.append(h1 + h2)
    return out


HI0 = hess_I(s0)


def sphere_hess(HF, Fval):  # Riemannian Hessian at s0 (|s0| = 1) of a degree-4 homogeneous F, on the 59-dim tangent space
    x0 = np.concatenate([np.real(s0), np.imag(s0)])
    Pt = np.eye(2 * NT) - np.outer(x0, x0)
    Ht = Pt @ (HF - 4 * Fval * np.eye(2 * NT)) @ Pt
    e = np.sort(np.linalg.eigvalsh(Ht))
    # drop one zero belonging to the radial direction (removed by the projector)
    k = np.argmin(np.abs(e))
    return np.delete(e, k)


def table_secondvar(Dl, l1, dd_coef):
    D2, D3, D4, D5 = Dl
    rows = [(4 * l1, 1), (0, 9), (4 / 3 * D2, 2), (D5, 6), (2 / 3 * D2 + D5, 6), (4 / 3 * D3, 12), (2 / 3 * D3 + D5, 12),
            (dd_coef[0] * D2 + dd_coef[1] * D3 + dd_coef[2] * D4 + D5, 12)]
    return np.sort(np.concatenate([[m] * n for m, n in rows]))


def tree_masses(lam, mu2):
    v2 = mu2 / (2 * lam[0])
    s = np.sqrt(v2) * s0
    H = -2 * mu2 * np.eye(2 * NT) + sum(l * h for l, h in zip(lam, hess_I(s)))
    return np.sort(np.linalg.eigvalsh(H)) / 2 / v2                          # m^2 / v^2


printed = (4 / 9, 2 / 9, 4 / 9)
corrected = (2 / 9, 2 / 9, 4 / 9)
dev_p, dev_c = 0.0, 0.0
for trial in range(6):
    lam_r = np.array([1.0, *(1 + 2 * rng.random(4))]) if trial else lam
    m = tree_masses(lam_r, 2.0)
    Dl_r = lam_r[1:] - lam_r[0]
    dev_p = max(dev_p, np.abs(m - table_secondvar(Dl_r, lam_r[0], printed)).max())
    dev_c = max(dev_c, np.abs(m - table_secondvar(Dl_r, lam_r[0], corrected)).max())
ok = dev_c < 1e-9
report("9", ok, f"thm:secondvar: the full 60x60 mass matrix (m^2 = Hessian/2 in coordinates orthonormal for Re Tr S^dag S', divided by v^2, mu^2 = 2) "
       f"at the benchmark and at 5 random couplings (seed 20261102) agrees with the table to {dev_c:.1e} when the d_Rd_R entry is "
       f"(2/9)D2 + (2/9)D3 + (4/9)D4 + D5; with the printed entry (4/9)D2 + ... the deviation is {dev_p:.2f}: THE PRINTED d_Rd_R ENTRY IS WRONG "
       "(the local-minimum condition 2D2+2D3+4D4+9D5 > 0 of the same theorem is 9 times the corrected entry) (float, tol 1e-9)")
ok = dev_p > 1e-3
report("10", ok, f"the printed d_Rd_R entry (4/9)D2+(2/9)D3+(4/9)D4+D5 disagrees with the computed spectrum by up to {dev_p:.2f} v^2 (this line records "
       "the error found in 9B; it passes when the disagreement is reproduced)")
m = tree_masses(lam, 2.0)
D2, D3, D4, D5 = Dl
ee, ue, de, uu, ud = 4 / 3 * D2, D5, 2 / 3 * D2 + D5, 4 / 3 * D3, 2 / 3 * D3 + D5
ok = abs(de - (ue + ee / 2)) < 1e-12 and abs(ud - (ue + uu / 2)) < 1e-12 and abs(m[-1] - max(table_secondvar(Dl, 1, corrected))) < 1e-9
ok &= abs(ee - 0.9333333333333) < 1e-9 and abs(de - 2.3666666666667) < 1e-9 and abs(uu - 1.7333333333333) < 1e-9 and abs(ud - 2.7666666666667) < 1e-9
report("11", ok, "eq:sumrules m_de^2 = m_ue^2 + m_ee^2/2 and m_ud^2 = m_ue^2 + m_uu^2/2 hold for the table entries, and ex:sumrule's values "
       "0.933, 1.9, 2.367, 1.733, 2.767 at the benchmark are reproduced by the 60x60 matrix (float)")

# ------------------------------------------------------------------ one loop on T: vector and fermion traces
def MV2(S):
    Ts = [X @ S + S @ X.T for X in G8]
    return np.array([[2 * np.real(np.trace(A.conj().T @ B)) for B in Ts] for A in Ts])


rows, rhs, rhsF = [], [], []
for _ in range(15):
    s = rand_s(); S = to_mat(s)
    M2 = MV2(S)
    rows.append(Ik(s)); rhs.append(np.trace(M2 @ M2)); rhsF.append(np.real(np.trace(np.linalg.matrix_power(S.conj().T @ S, 2))))
rows = np.array(rows)
cfit, res, *_ = np.linalg.lstsq(rows, np.array(rhs), rcond=None)
afit = np.linalg.lstsq(rows, np.array(rhsF), rcond=None)[0]
resid = np.abs(rows @ cfit - rhs).max() / np.abs(rhs).max()
ok = np.allclose(cfit, [33, 3, -3, 39, 5], atol=1e-8) and resid < 1e-10 and np.allclose(afit, a_rec, atol=1e-9)
c_vec = np.array([33, 3, -3, 39, 5])
Gk = -3 * (c_vec[1:] - c_vec[0]); Fk = -2 * (a_rec[1:] - a_rec[0])
ok &= np.allclose(Gk, [90, 108, -18, 84]) and np.allclose(Fk, [3, 3, 1.5, 2])
report("12", ok, f"thm:oneloop1013: with (M_V^2)_ab = 2 Re Tr(T_aS+ST_a^T)^dag(T_bS+ST_b^T) (g = 1, one coupling for su4 and su2_R), Tr M_V^4 = "
       f"sum c_k I_k with c = {np.round(cfit, 9).tolist()} (least squares on 15 random S, residual {resid:.1e}) and Tr(S^dag S)^2 = sum a_k I_k; "
       "hence G = -3(c_k - c_1) = (90,108,-18,84) and F = -2(a_k - a_1) = (3,3,3/2,2) (float)")


def oneloop_masses(x):
    ell = 2 * x * a_rec - 3 * c_vec                                         # coefficients of I_k in delta V / (L g^4 / 64 pi^2)
    HF = sum(l * h for l, h in zip(ell, HI0))
    return sphere_hess(HF, ell @ f0) / 2                                    # m^2 in units of L g^4 v^2 / (64 pi^2)


def oneloop_table(x):
    rows = [(0, 9), (120 - 4 * x, 2), (84 - 2 * x, 6), (144 - 4 * x, 6), (144 - 4 * x, 12), (156 - 4 * x, 12), (120 - 4 * x, 12)]
    return np.sort(np.concatenate([[m] * n for m, n in rows]))


devs = [np.abs(oneloop_masses(x) - oneloop_table(x)).max() for x in (0, 5, 17.5, 29)]
m0 = oneloop_masses(0)
ok = max(devs) < 1e-8 and np.sum(np.abs(m0 - 120) < 1e-8) == 14 and np.sum(np.abs(m0 - 140) < 1e-6) == 0
report("13", ok, f"one-loop angular masses on T at S_0 (units L g^4 v^2/64pi^2), x = 0, 5, 17.5, 29: nine Goldstones, e_Re_R 120-4x (2), u_Re_R 84-2x (6), "
       f"d_Re_R 144-4x (6), u_Ru_R 144-4x (12), u_Rd_R 156-4x (12), d_Rd_R 120-4x (12), max deviation {max(devs):.1e} (tol 1e-8). The quoted "
       "x = 0 masses 120, 84, 144, 144, 156, 120 are reproduced; d_Rd_R is 120 (not 140, the value the printed entry of thm:secondvar would give)")

xs = np.linspace(0, 32, 321)
first_neg = next((x for x in xs if oneloop_masses(x).min() < -1e-9), None)
exc = lambda x: (2 * x * a_rec - 3 * c_vec) @ (f1 - f0)
ok = abs(first_neg - 30.1) < 1e-9 and abs(exc(0) - 48) < 1e-9 and abs(exc(192 / 7)) < 1e-9 and abs(exc(1) - (48 - 7 / 4)) < 1e-9
report("14", ok, f"thm:oneloop1013: on T the Standard-Model orbit is a non-degenerate local minimum for x < 30 (first negative mass on the grid "
       f"of step 0.1 at x = {first_neg:.1f}), and 1(x)1 has excess 48 - 7x/4 over it, zero at x = 192/7 (float)")

print(f"{NPASS[0]}/{NPASS[0] + len(FAILS)} checks passed")
sys.exit(0 if not FAILS else 1)
