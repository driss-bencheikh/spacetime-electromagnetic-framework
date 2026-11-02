"""Release audit (v1.0.0): recomputation of the electroweak-stage statements of Chapter 7 that had no check (numerical).

    python3 verification/ch07_mass_vacuum/a8_electroweak_numerical.py

Float64 linear algebra on the exactly specified integer operators of kernel_space.py (the 64-dimensional block V_* = im P_*, the
272-dimensional kernel space, the 21 Pati-Salam generators, J, Gamma, Phi_0).  Every dimension below is a numerical rank: the script
asserts a spectral gap and prints it (direct SVDs: largest discarded singular value < 1e-9; commutants, whose singular values are square
roots of Gram eigenvalues: < 1e-6; smallest kept > 1e-3 in both cases).  These are numerical recomputations, not exact certificates.

a8-1  prop:tensor     the (10,1,3) component is J-invariant of real dimension 60 with PS-commutant span{1, J} (complex irreducible); for
                      Phi in it Phi^T Phi vanishes on the Gamma = -1 half and its spectrum on the Gamma = +1 half is fourfold; Phi_0 lies
                      in it with Phi_0^T Phi_0 of a single eigenvalue 1, four times.
a8-2  prop:EW (i)     kernels commuting with su(3) and Q = T^L_3 + T^R_3 + (B-L)/2: dimension 12 = 2 + 2 + 4 + 4 over (10,1,3), (10,3,1),
                      (15,2,2), (1,2,2); up-type and down-type spaces of dimension two in each bidoublet; charge sectors 24, 24, 8, 8.
a8-3  prop:EW (ii)    (1,2,2): m_u = m_nu, m_d = m_e; (15,2,2) aligned with (1,2,2) on the quarks gives -3 on the leptons; the mass formulas.
a8-4  prop:EW (iii)   seesaw singular values (sqrt(v_R^2 + 4 kappa^2) -+ v_R)/2, each four times, on the neutral sector.
a8-5  prop:EW (iv)    stabiliser of v_R Phi_0 + v_u K_1u + v_d K_1d among the 21 generators: dimension 9, Q annihilates Phi.
a8-6  thm:mA          a, e_u, e_d kernels, orthogonal, norms 4, 32, 32; t = [J Gamma, Phi] = 2 v_R a + 2 v_u e_u + 2 v_d e_d; the gauge
                      tangent space meets span{a, e_u, e_d} in span{a, z} and splits orthogonally; |P_perp t|^2 = 16 n_K v_u^2 v_d^2/(v_u^2+v_d^2).
a8-7  prop:rankfive   on a spin part (32 dimensions) the commutant of J, Gamma and the 12 Standard-Model generators has dimension 12, its
                      symmetric part dimension 6, commutative, commuting with T^R_3 and B-L, with multiplicities 12, 4, 6, 6, 2, 2.
a8-8  prop:bidoublet  three PS- and J-invariant symmetric forms on (1,2,2), one J Gamma-invariant; CP eigenvalues (+1, +1, -1), the
                      invariant form CP-even; along e^{theta ad J Gamma} the breaking forms are pure fourth harmonics of amplitude
                      proportional to v_u v_d and B_inv is constant.
a8-9  (for a6-5)      charges of the components under ad J (fermion number) and ad J Gamma, used by the graded count of thm:Bcounts.
Each check prints PASS/FAIL a8-<n>; exit status 0 iff all pass.
"""
import itertools
import os
import sys
from functools import reduce

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import kernel_space as ks  # noqa: E402

NPASS, FAILS = [0], []
GAPS = []


def report(cid, ok, text):
    if ok:
        NPASS[0] += 1
        print(f"PASS a8-{cid}: {text}")
    else:
        FAILS.append(cid)
        print(f"FAIL a8-{cid}: {text}")
    sys.stdout.flush()


def nullsp(M, lo=1e-9, hi=1e-3):
    """Orthonormal null space of M with an asserted spectral gap (numerical rank)."""
    _, s, Vt = np.linalg.svd(M, full_matrices=M.shape[0] < M.shape[1])
    s_full = np.concatenate([s, np.zeros(Vt.shape[0] - s.size)])
    kept = s_full > np.sqrt(lo * hi)
    small, big = (s_full[~kept].max() if (~kept).any() else 0.0), (s_full[kept].min() if kept.any() else np.inf)
    GAPS.append((small, big))
    assert small < lo and big > hi, (small, big)
    return Vt[~kept].T


def commutant(mats, n):
    """Null space of X -> [X, g] over all g, from the accumulated Gram matrix sum L_g^T L_g (eigenvalues = squared singular values)."""
    G = np.zeros((n * n, n * n))
    for g in mats:
        L = np.kron(g, np.eye(n)) - np.kron(np.eye(n), g.T)
        G += L.T @ L
    w, V = np.linalg.eigh(G)
    sv_ = np.sqrt(np.clip(w, 0, None))
    kept = sv_ > 1e-6
    small, big = (sv_[~kept].max() if (~kept).any() else 0.0), (sv_[kept].min() if kept.any() else np.inf)
    GAPS.append((small, big))
    assert small < 1e-6 and big > 1e-3, (small, big)
    return V[:, ~kept]


def eigsp(M, val, basis, tol=1e-8):
    Mr = basis.T @ M @ basis
    w, V = np.linalg.eigh((Mr + Mr.T) / 2)
    return basis @ V[:, np.abs(w - val) < tol]


A, d, B = ks.A, ks.dimK, ks.B
ip = lambda X, Y: float(np.sum(X * Y))  # noqa: E731
CL = sum(A[i] @ A[i] for i in range(3))
CR = sum(A[i] @ A[i] for i in range(3, 6))
C4 = sum(A[i] @ A[i] for i in range(6, 21))
comp = {}
for name, (l_, r_, c_) in {"1013": (0, -2, -36), "1031": (-2, 0, -36), "1522": (-.75, -.75, -32), "122": (-.75, -.75, 0),
                           "611": (0, 0, -20)}.items():
    b = np.eye(d)
    b = eigsp(CL, l_, b)
    b = eigsp(CR, r_, b)
    comp[name] = eigsp(C4, c_, b)
dims = {k: v.shape[1] for k, v in comp.items()}

Cidx = [10, 11, 12, 16, 17, 18]
cpairs = list(itertools.combinations(Cidx, 2))
HCc = ks.cp(ks.HC)
col = ks.tPS[6:]
cen = nullsp(np.array([(c @ HCc - HCc @ c).ravel() for c in col]).T).T                 # centraliser of HC in su(4): su(3) + u(1)
cenM = [sum(x * c for x, c in zip(v, col)) for v in cen]
_, s, Vt = np.linalg.svd(np.array([(a @ b - b @ a).ravel() for a, b in itertools.combinations(cenM, 2)]), full_matrices=False)
su3 = [Vt[k].reshape(64, 64) for k in range(int(np.sum(s > 1e-6)))]
adK = lambda t: np.array([ks.coords(t @ B[i] - B[i] @ t) for i in range(d)]).T  # noqa: E731
Yg = ks.tPS[5] + HCc / 6                                                              # T^R_3 + (B-L)/2 (generator form; B-L = -J HC/3)
Qg = ks.tPS[2] + Yg
SM = su3 + ks.tPS[0:3] + [Yg]
JG = ks.Jv @ ks.Gv

# ---------------------------------------------------------------- a8-1 prop:tensor
c13 = comp["1013"]
Jl = c13.T @ np.array([ks.coords(ks.Jv @ ks.mat(c13[:, j])) for j in range(60)]).T
ok = dims == {"1013": 60, "1031": 60, "1522": 120, "122": 8, "611": 24} and len(cen) == 9 and len(su3) == 8
ok &= max(np.linalg.norm(ks.Jv @ ks.mat(c13[:, j]) - ks.mat(c13 @ Jl[:, j])) for j in range(60)) < 1e-10        # J-invariant
gens13 = [c13.T @ A[a] @ c13 for a in range(21)]
comm13 = commutant(gens13, 60)
ok &= comm13.shape[1] == 2                                                          # span{1, J}: complex irreducible
rng = np.random.default_rng(2026)
spec_ok = True
for _ in range(3):
    Phi = ks.mat(c13 @ rng.normal(size=60))
    M = Phi.T @ Phi
    spec_ok &= np.linalg.norm((np.eye(64) - ks.Gv) / 2 @ M) < 1e-10
    ev = np.sort(np.linalg.eigvalsh(M))[32:]
    spec_ok &= np.allclose(ev.reshape(8, 4), ev.reshape(8, 4)[:, :1], atol=1e-9) and np.min(np.diff(ev.reshape(8, 4)[:, 0])) > 1e-4
ok &= spec_ok
ok &= np.linalg.norm(c13 @ (c13.T @ ks.c0) - ks.c0) < 1e-10
ev0 = np.sort(np.linalg.eigvalsh(ks.Phi0.T @ ks.Phi0))
ok &= np.allclose(ev0[-4:], 1) and np.allclose(ev0[:-4], 0, atol=1e-10)
report("1", ok, f"prop:tensor: components {dims}; (10,1,3) is J-invariant, its Pati-Salam commutant is span{{1, J}} (complex irreducible, "
       "so the equivariant isomorphism to T is unique up to a complex scale); Phi^T Phi vanishes on Gamma = -1 and has fourfold spectrum on "
       "Gamma = +1 (three random Phi); Phi_0 is in (10,1,3) with Phi_0^T Phi_0 = 1 four times (numerical)")

# ---------------------------------------------------------------- a8-2 prop:EW (i)
AQ = adK(Qg)
nul = nullsp(np.vstack([adK(t) for t in su3] + [AQ]))
inter = {k: int(np.linalg.matrix_rank(v.T @ nul, 1e-8)) for k, v in comp.items()}
q2, Vq = np.linalg.eigh(-(Qg @ Qg))
Psec = {k: Vq[:, np.abs(q2 - v) < 1e-6] for k, v in dict(u=4 / 9, d=1 / 9, nu=0.0, e=1.0).items()}
Pm = {k: v @ v.T for k, v in Psec.items()}


def typed(name, kill):
    U, s_, _ = np.linalg.svd(nul @ (nul.T @ comp[name]), full_matrices=False)
    sp_ = U[:, s_ > 0.5]
    cond = np.vstack([np.array([(Pm[k] @ ks.mat(sp_[:, j])).ravel() for j in range(sp_.shape[1])]).T for k in kill])
    return sp_ @ nullsp(cond)


up1, dn1, up15, dn15 = typed("122", ("d", "e")), typed("122", ("u", "nu")), typed("1522", ("d", "e")), typed("1522", ("u", "nu"))
ok = nul.shape[1] == 12 and inter == {"1013": 2, "1031": 2, "1522": 4, "122": 4, "611": 0}
ok &= {k: v.shape[1] for k, v in Psec.items()} == {"u": 24, "d": 24, "nu": 8, "e": 8}
ok &= [x.shape[1] for x in (up1, dn1, up15, dn15)] == [2, 2, 2, 2]
ok &= np.linalg.norm(AQ @ ks.c0) < 1e-10
report("2", ok, f"prop:EW (i): the kernels commuting with su(3) and Q form a {nul.shape[1]}-dimensional space, by component {inter}; the "
       "charge sectors u, d, nu, e have real dimensions 24, 24, 8, 8; in each bidoublet the up-type (no d, e entries) and down-type (no u, nu "
       "entries) kernels form two-dimensional spaces; Q annihilates Phi_0 (numerical)")

# ---------------------------------------------------------------- a8-3 prop:EW (ii)
sv = lambda M, k: np.linalg.svd(Psec[k].T @ M @ Psec[k], compute_uv=False)  # noqa: E731
unit32 = lambda M: M * np.sqrt(32) / np.linalg.norm(M)  # noqa: E731
K1u, K1d = unit32(ks.mat(up1[:, 0])), unit32(ks.mat(dn1[:, 0]))
ok = all(np.allclose(sv(K1u, k), 1) for k in ("u", "nu")) and all(np.allclose(sv(K1u, k), 0) for k in ("d", "e"))
ok &= all(np.allclose(sv(K1d, k), 1) for k in ("d", "e")) and all(np.allclose(sv(K1d, k), 0) for k in ("u", "nu"))


def align(sp15, K1, q, l_):
    Vq_, Vl = Psec[q], Psec[l_]
    Mx = np.array([(Vq_.T @ ks.mat(sp15[:, j]) @ Vq_).ravel() for j in range(sp15.shape[1])]).T
    x, *_ = np.linalg.lstsq(Mx, (Vq_.T @ K1 @ Vq_).ravel(), rcond=None)
    K15 = ks.mat(sp15 @ x)
    return K15, np.linalg.norm(Mx @ x - (Vq_.T @ K1 @ Vq_).ravel()), np.linalg.norm(Vl.T @ K15 @ Vl + 3 * Vl.T @ K1 @ Vl)


K15u, r1, r2 = align(up15, K1u, "u", "nu")
K15d, r3, r4 = align(dn15, K1d, "d", "e")
ok &= max(r1, r2, r3, r4) < 1e-10
for k1, k15 in [(1.0, 0.3), (0.7, -1.1), (-0.4, 2.5)]:
    ok &= np.allclose(sv(k1 * K1u + k15 * K15u, "u"), abs(k1 + k15)) and np.allclose(sv(k1 * K1u + k15 * K15u, "nu"), abs(k1 - 3 * k15))
    ok &= np.allclose(sv(k1 * K1d + k15 * K15d, "d"), abs(k1 + k15)) and np.allclose(sv(k1 * K1d + k15 * K15d, "e"), abs(k1 - 3 * k15))
report("3", ok, "prop:EW (ii): with |K|^2 = n_K = 32 the (1,2,2) kernels have all singular values 1 on u and nu (up type) or d and e (down "
       "type), so m_u = m_nu^D and m_d = m_e; the (15,2,2) kernel equal to K_1 on the quarks is -3 K_1 on the leptons, and "
       "m_quark = |k_1 + k_15|, m_lepton = |k_1 - 3 k_15| at three sample points (numerical)")

# ---------------------------------------------------------------- a8-4 prop:EW (iii)
ok = True
for vR, ka in [(20.0, 3.0), (5.0, 1.7), (1.0, 4.0)]:
    s_ = np.sort(sv(vR * ks.Phi0 + ka * K1u, "nu"))
    r = np.sqrt(vR ** 2 + 4 * ka ** 2)
    ok &= np.allclose(s_[:4], (r - vR) / 2) and np.allclose(s_[4:], (r + vR) / 2)
report("4", ok, "prop:EW (iii): on the neutral sector v_R Phi_0 + kappa K_1u has singular values (sqrt(v_R^2 + 4 kappa^2) -+ v_R)/2, "
       "each four times, at (v_R, kappa) = (20, 3), (5, 1.7), (1, 4); at (20, 3): 0.440307 and 20.440307 as in ex:EWspectrum (numerical)")

# ---------------------------------------------------------------- a8-5 prop:EW (iv) and a8-6 thm:mA
ok5, ok6 = True, True
details = []
for (vR, vu, vd), (cu, cd) in [((20.0, 1.0, 0.5), (0.0, 0.0)), ((3.0, 2.0, 1.3), (0.4, 1.9)), ((1.0, 0.3, 2.2), (2.5, -0.7))]:
    Ku = unit32(ks.mat(up1 @ np.array([np.cos(cu), np.sin(cu)])))
    Kd = unit32(ks.mat(dn1 @ np.array([np.cos(cd), np.sin(cd)])))
    Phi = vR * ks.Phi0 + vu * Ku + vd * Kd
    c = ks.coords(Phi)
    Z = np.array([A[a] @ c for a in range(21)])
    stab = nullsp(Z.T)
    ok5 &= stab.shape[1] == 9 and np.linalg.norm(Qg @ Phi - Phi @ Qg) < 1e-10
    a_, eu, ed = ks.Jv @ ks.Phi0, JG @ Ku, JG @ Kd
    ok6 &= max(np.linalg.norm(X - ks.mat(ks.coords(X))) for X in (a_, eu, ed)) < 1e-10
    ok6 &= np.allclose([ip(a_, a_), ip(eu, eu), ip(ed, ed)], [4, 32, 32]) and max(abs(ip(a_, eu)), abs(ip(a_, ed)), abs(ip(eu, ed))) < 1e-10
    t = JG @ Phi - Phi @ JG
    ok6 &= np.linalg.norm(t - 2 * vR * a_ - 2 * vu * eu - 2 * vd * ed) < 1e-10
    U, sg, _ = np.linalg.svd(Z.T, full_matrices=False)
    Ug = U[:, sg > 1e-6]
    S = np.array([ks.coords(X) for X in (a_, eu, ed)]).T
    PG, PS_ = Ug @ Ug.T, S @ np.linalg.pinv(S)
    ok6 &= nullsp(np.hstack([S, -Ug])).shape[1] == 2 and np.linalg.norm(PG @ PS_ @ Ug - PS_ @ Ug) < 1e-10
    zc = ks.coords(vu * eu - vd * ed)
    ok6 &= np.linalg.norm(PG @ ks.coords(a_) - ks.coords(a_)) < 1e-10 and np.linalg.norm(PG @ zc - zc) < 1e-10
    tc = ks.coords(t)
    perp = tc - PG @ tc
    val, target = float(perp @ perp), 16 * 32 * vu ** 2 * vd ** 2 / (vu ** 2 + vd ** 2)
    ok6 &= abs(val - target) < 1e-8 * max(1.0, target)
    details.append(f"{val:.6f}")
report("5", ok5, "prop:EW (iv): for v_R Phi_0 + v_u K_1u + v_d K_1d (three points, independent phases of K_1u and K_1d in their planes) the "
       "stabiliser among the 21 generators has dimension 9 (su(3) + u(1)_Q) and [Q, Phi] = 0 (numerical)")
report("6", ok6, "thm:mA: a = J Phi_0, e_u = J Gamma K_1u, e_d = J Gamma K_1d are kernels, mutually orthogonal, |a|^2 = 4, |e_u|^2 = "
       "|e_d|^2 = 32; t = [J Gamma, Phi] = 2 v_R a + 2 v_u e_u + 2 v_d e_d; the gauge tangent space meets span{a, e_u, e_d} in span{a, z} and "
       f"splits orthogonally; |P_perp t|^2 = 512 v_u^2 v_d^2/(v_u^2 + v_d^2) at three points ({', '.join(details)}) (numerical)")

# ---------------------------------------------------------------- a8-7 prop:rankfive
S3 = ks.Jv @ ks.rot[0]                                                                 # symmetric, eigenvalues +-1: spin up/down
w, V = np.linalg.eigh((S3 + S3.T) / 2)
ok = np.abs(S3 - S3.T).max() < 1e-12 and np.allclose(np.sort(np.abs(w)), 1)
mults_all = []
for sgn in (1, -1):
    Up = V[:, np.abs(w - sgn) < 1e-6]
    ops = [ks.Jv, ks.Gv] + SM
    ok &= Up.shape[1] == 32 and max(np.abs(o @ S3 - S3 @ o).max() for o in ops) < 1e-12
    R = [Up.T @ o @ Up for o in ops]
    co = commutant(R, 32)
    Ms = [co[:, k].reshape(32, 32) for k in range(co.shape[1])]
    sym = nullsp(np.array([(M - M.T).ravel() for M in Ms]).T)
    Sy = [sum(x * M for x, M in zip(v, Ms)) for v in sym.T]
    tr3, hc = Up.T @ ks.tPS[5] @ Up, Up.T @ HCc @ Up
    ok &= len(Ms) == 12 and len(Sy) == 6
    ok &= max(np.abs(a @ b - b @ a).max() for a in Sy for b in Sy) < 1e-10
    ok &= max(np.abs(a @ tr3 - tr3 @ a).max() + np.abs(a @ hc - hc @ a).max() for a in Sy) < 1e-10
    g = sum(x * M for x, M in zip(rng.normal(size=len(Sy)), Sy))
    ev = np.sort(np.linalg.eigvalsh(g))
    groups = []
    for x in ev:
        if groups and abs(x - groups[-1][0]) < 1e-6:
            groups[-1][1] += 1
        else:
            groups.append([x, 1])
    mults = sorted(m for _, m in groups)
    mults_all.append(mults)
    ok &= mults == [2, 2, 4, 6, 6, 12]
report("7", ok, "prop:rankfive: on each spin half of the block (eigenspaces of J e4e5, 32 dimensions) the operators commuting with J, Gamma "
       "and the 12 Standard-Model generators form a 12-dimensional algebra with 6-dimensional symmetric part, commutative and commuting with "
       f"T^R_3 and B-L; a generic symmetric element has six eigenspaces of dimensions {mults_all[0]} (Q, L, u, d, e, nu: pairwise "
       "inequivalent) (numerical)")

# ---------------------------------------------------------------- a8-8 prop:bidoublet
bas = comp["122"]


def onc(f):
    return bas.T @ np.array([ks.coords(f(ks.mat(bas[:, j]))) for j in range(8)]).T


gens = [bas.T @ A[a] @ bas for a in range(21)]
Jc = onc(lambda M: ks.Jv.T @ M @ ks.Jv)
JGc = onc(lambda M: JG @ M)                                                            # on (1,2,2): ad(J Gamma) = 2 J Gamma
CP512 = reduce(np.matmul, [ks.E[i] for i in (1, 2, 7, 8, 10, 11, 12, 14)], ks.I512).astype(float)
CPv = ks.W.T @ CP512 @ ks.W
ok = np.linalg.norm(ks.W @ CPv @ ks.W.T - ks.Pst @ CP512 @ ks.Pst) < 1e-10 and np.linalg.norm(CPv @ CPv.T - np.eye(64)) < 1e-10
CPc = onc(lambda M: CPv @ M @ CPv.T)
ok &= np.linalg.norm(ks.mat(bas @ CPc[:, 0]) - CPv @ ks.mat(bas[:, 0]) @ CPv.T) < 1e-10
ok &= np.linalg.norm(ks.mat(bas @ (2 * JGc[:, 0])) - (JG @ ks.mat(bas[:, 0]) - ks.mat(bas[:, 0]) @ JG)) < 1e-10
symb = []
for i in range(8):
    for j in range(i, 8):
        Sb = np.zeros((8, 8))
        Sb[i, j] = Sb[j, i] = 1
        symb.append(Sb)
lin = lambda f: np.array([f(Sb).ravel() for Sb in symb]).T  # noqa: E731
inv = nullsp(np.vstack([lin(lambda Sb, X=X: X.T @ Sb + Sb @ X) for X in gens] + [lin(lambda Sb: Jc.T @ Sb @ Jc - Sb)]))
Bs = [sum(x * Sb for x, Sb in zip(v, symb)) for v in inv.T]
ok &= len(Bs) == 3
Bm = np.array([b.ravel() for b in Bs]).T
coef = lambda M: np.linalg.lstsq(Bm, M.ravel(), rcond=None)[0]  # noqa: E731
Lad = np.array([coef(2 * (JGc.T @ b + b @ JGc)) for b in Bs]).T                       # Lie derivative of forms along ad J Gamma
Cmat = np.array([coef(CPc.T @ b @ CPc) for b in Bs]).T
evL = np.sort_complex(np.round(np.linalg.eigvals(Lad), 9))
ok &= np.allclose(evL, [-4j, 0, 4j])
B_inv_c = nullsp(Lad)[:, 0]
ok &= np.allclose(Cmat @ B_inv_c, B_inv_c)                                             # the invariant form is CP-even
ok &= np.allclose(np.sort(np.linalg.eigvals(Cmat).real), [-1, 1, 1])
ok &= np.allclose(Cmat @ Lad, -Lad @ Cmat)                                             # CP reverses the J Gamma charge
# breaking space: the Lad-image (charges +-4); CP eigenvectors there
brk = np.linalg.svd(Lad)[0][:, :2]
cw_, cv_ = np.linalg.eig(brk.T @ Cmat @ brk)
Be_c = brk @ cv_[:, np.argmax(cw_.real)].real
Bo_c = -Lad @ Be_c / 4                                                                 # so that d/dtheta B_e(Phi_theta) = -4 B_o(Phi_theta)
Be = sum(x * b for x, b in zip(Be_c, Bs))
Bo = sum(x * b for x, b in zip(Bo_c, Bs))
Binv = sum(x * b for x, b in zip(B_inv_c, Bs))
ok &= np.allclose(Cmat @ Be_c, Be_c) and np.allclose(Cmat @ Bo_c, -Bo_c)
xu, xd = bas.T @ up1[:, 0], bas.T @ dn1[:, 0]
ok &= all(abs(x @ F @ x) < 1e-10 for x in (xu, xd) for F in (Be, Bo))                 # no up-up or down-down term: amplitude ~ v_u v_d
th = np.linspace(0, np.pi / 2, 9)
Rot = lambda t_: np.cos(2 * t_) * np.eye(8) + np.sin(2 * t_) * JGc  # noqa: E731      # e^{theta ad J Gamma} on (1,2,2)
ok &= np.allclose(JGc @ JGc, -np.eye(8))
for vu, vd in [(1.0, 0.5), (2.0, -1.3)]:
    x0 = vu * xu + vd * xd
    be = np.array([(Rot(t_) @ x0) @ Be @ (Rot(t_) @ x0) for t_ in th])
    bo = np.array([(Rot(t_) @ x0) @ Bo @ (Rot(t_) @ x0) for t_ in th])
    bi = np.array([(Rot(t_) @ x0) @ Binv @ (Rot(t_) @ x0) for t_ in th])
    amp, ph = np.hypot(be[0], bo[0]), np.arctan2(bo[0], be[0])
    ok &= np.allclose(be, amp * np.cos(4 * th + ph)) and np.allclose(bo, amp * np.sin(4 * th + ph)) and np.allclose(bi, bi[0])
    ok &= amp > 1e-6
    ref = amp / abs(vu * vd)
ok &= abs(np.hypot(*[(xu + xd) @ F @ (xu + xd) for F in (Be, Bo)]) - ref) < 1e-9 * ref     # amplitude proportional to |v_u v_d|
report("8", ok, "prop:bidoublet: on (1,2,2) the Pati-Salam- and J-invariant symmetric forms are 3-dimensional; under ad J Gamma they carry "
       "charges 0, +-4, the invariant one is CP-even, CP has eigenvalues (+1, +1, -1) and reverses the charge, so the breaking space splits "
       "into a CP-even B_e and a CP-odd B_o; along e^{theta ad J Gamma}(v_u K_u + v_d K_d) they are A cos(4 theta + phi) and A sin(4 theta + phi) "
       "with A proportional to |v_u v_d| (no K_u K_u or K_d K_d term), and B_inv is constant (numerical; the normalisation and the phase "
       "phi = pi of the text fix the scale of B_e and the relative phase of K_1u, K_1d)")

# ---------------------------------------------------------------- a8-9 charges of the components under ad J and ad J Gamma (for a6-5)
def ad_on(t, bas):
    return bas.T @ np.array([ks.coords(t @ ks.mat(bas[:, j]) - ks.mat(bas[:, j]) @ t) for j in range(bas.shape[1])]).T


ok = True
signs = {}
for k_, bas in comp.items():
    aJ, aG = ad_on(ks.Jv, bas), ad_on(JG, bas)
    ok &= np.abs(aJ @ aG - aG @ aJ).max() < 1e-10 and np.allclose(aG @ aG, -4 * np.eye(bas.shape[1]))
    if k_ in ("1522", "122"):
        ok &= np.abs(aJ).max() < 1e-10                                                 # Dirac type: F = 0
        signs[k_] = "F=0, q=+-2"
    else:
        ok &= np.allclose(aJ @ aJ, -4 * np.eye(bas.shape[1]))                           # Majorana type: F = 2
        ev = np.linalg.eigvalsh((aG @ aJ.T + aJ @ aG.T) / 8)                            # = eps on the eps-eigenspaces of ad(J Gamma) ad(J)^-1
        signs[k_] = (int(np.sum(np.isclose(ev, 1))), int(np.sum(np.isclose(ev, -1))))
        ok &= signs[k_][0] + signs[k_][1] == bas.shape[1]
ok &= signs["1013"] == (60, 0) and signs["1031"] == (0, 60) and signs["611"] == (12, 12)
report("9", ok, f"charges for thm:Bcounts (a6-5): ad J and ad J Gamma commute and square to -4 or 0; ad(J Gamma) = eps ad J on the "
       f"Majorana-type components with (dim eps=+1, dim eps=-1) = {signs} (numerical)")

print(f"spectral gaps (largest discarded, smallest kept singular value) over all rank decisions: "
      f"max discarded {max(g[0] for g in GAPS):.2e}, min kept {min(g[1] for g in GAPS):.2e}")
print(f"{NPASS[0]}/{NPASS[0] + len(FAILS)} checks passed")
sys.exit(0 if not FAILS else 1)
