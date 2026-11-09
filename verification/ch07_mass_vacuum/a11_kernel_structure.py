#!/usr/bin/env python3
"""Corrections audit (v1.0.1): the complex structure of the kernel space (Theorem 7.17) and the v_R = 0 strata of Theorem 7.42.

    python3 verification/ch07_mass_vacuum/a11_kernel_structure.py

Uses the kernel space of kernel_space.py (the 272 Lorentz-invariant mass kernels of the chiral block V_* = im P_*, an orthonormal float64
basis on the 64-dimensional block). Numerical: float64 linear algebra on exactly specified integer matrices, tolerances as printed.

a11-1  thm:kernelspace   C~(M) = J M Gamma maps every kernel to a kernel and squares to -1 on the kernel space; it commutes with the action
                         M -> [t, M] of the 21 Pati-Salam generators, so (M, C~) is a complex representation; for the colour Cartan element H
                         the complex cubic trace (1/2) tr_R((-C~ L_H)^3) vanishes, as on Sym^2(16), while the structure of version 1.0.0
                         (J M on the Majorana-type kernels, J M Gamma on the Dirac-type ones) gives +-1152.
a11-2  thm:axion         rows v_R = 0 of the table: for Phi = v_u K_1u + v_d K_1d with |K|^2 = 32 (K_1d swept over its two-dimensional
                         space of rotations), the stabilisers in the 21 Pati-Salam generators and in the 23 generators (with J and J Gamma),
                         hence the massive vectors, the orbit dimension and the Goldstone count, are (16,17,5,6,1) for v_u v_d != 0 and
                         |v_u| != |v_d|, (18,19,3,4,1) for |v_u| = |v_d| != 0, (16,18,5,5,0) when one value is zero, and (21,23,0,0,0) at Phi = 0.
Each check prints PASS/FAIL a11-<n>; exit status 0 iff all pass.
"""
import os
import sys
from collections import Counter

import numpy as np
from scipy.linalg import null_space, orth

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import kernel_space as ks  # noqa: E402

NPASS, FAILS = [0], []


def report(cid, ok, text):
    if ok:
        NPASS[0] += 1
        print(f"PASS a11-{cid}: {text}")
    else:
        FAILS.append(cid)
        print(f"FAIL a11-{cid}: {text}")
    sys.stdout.flush()


def nullsp(A, tol=1e-9):
    """null space of A; for a tall A through the Gram matrix A^T A (same null space), relative tolerance on its eigenvalues"""
    A = np.asarray(A)
    if A.shape[0] > A.shape[1]:
        w, V = np.linalg.eigh(A.T @ A)
        return V[:, w < tol * max(1.0, float(w.max()))]
    return null_space(A, rcond=tol)


B, Bv, dimK, n = ks.B, ks.Bv, ks.dimK, ks.n
assert dimK == 272
Jc, Gc = ks.Jv, ks.Gv                     # J and Gamma compressed to the block (64 x 64)
coord = lambda M: Bv @ M.ravel()          # noqa: E731
mats = [B[k] for k in range(dimK)]


def opmat(f):
    return np.array([coord(f(M)) for M in mats]).T


def resid_in_space(f):
    worst = 0.0
    for M in mats:
        X = f(M)
        worst = max(worst, np.linalg.norm(Bv.T @ coord(X) - X.ravel()) / max(1.0, np.linalg.norm(X)))
    return worst


# ------------------------------------------------------------------ a11-1
Ct = opmat(lambda M: Jc @ M @ Gc)
r_in = resid_in_space(lambda M: Jc @ M @ Gc)
r_sq = np.linalg.norm(Ct @ Ct + np.eye(dimK))
act = [opmat(lambda M, t=t: t @ M - M @ t) for t in ks.tPS]
r_comm = max(np.linalg.norm(Ct @ L - L @ Ct) for L in act)
H = ks.cp(ks.HC)
LH = opmat(lambda M: H @ M - M @ H)
hh = -Ct @ LH
hh = (hh + hh.T) / 2
cub_new = float(np.sum(np.linalg.eigvalsh(hh) ** 3) / 2)
# the structure of version 1.0.0: J M on Majorana-type kernels (anticommuting with J), J M Gamma on Dirac-type kernels (commuting with J)
maj = null_space(np.array([(M @ Jc + Jc @ M).ravel() for M in mats]).T, rcond=1e-10)
dir_ = null_space(np.array([(M @ Jc - Jc @ M).ravel() for M in mats]).T, rcond=1e-10)
Pm, Pd = maj @ maj.T, dir_ @ dir_.T
Cold = opmat(lambda M: Jc @ M) @ Pm + Ct @ Pd
ho = -Cold @ LH
ho = (ho + ho.T) / 2
cub_old = float(np.sum(np.linalg.eigvalsh(ho) ** 3) / 2)
ok1 = (r_in < 1e-9 and r_sq < 1e-8 and r_comm < 1e-8 and abs(cub_new) < 1e-6 and abs(abs(cub_old) - 1152) < 1e-6
       and maj.shape[1] + dir_.shape[1] == dimK)
report(1, ok1, f"C~(M) = J M Gamma preserves the kernel space (residual {r_in:.1e}), C~^2 = -1 ({r_sq:.1e}), commutes with the 21 Pati-Salam "
               f"generators ({r_comm:.1e}); complex cubic trace of the colour Cartan element H over (M, C~) = {cub_new:.2e} (zero, as on Sym^2(16)); "
               f"the version 1.0.0 structure (Majorana {maj.shape[1]}, Dirac {dir_.shape[1]}) gives {cub_old:.3f} (numerical, float64)")

# ------------------------------------------------------------------ a11-2
I512 = np.eye(512)
colour = ks.colour
HC = ks.HC
cen = [sum(v[i] * colour[i] for i in range(15)) for v in nullsp(np.array([(c @ HC - HC @ c).ravel() for c in colour]).T).T]
su3 = [v.reshape(512, 512) for v in orth(np.array([(a @ b - b @ a).ravel() for a in cen for b in cen]).T).T]
tY = ks.tR[2] + HC / 6
tQ = ks.tL[2] + tY
cp = ks.cp
sm_gens = [cp(t) for t in su3] + [cp(tQ)]
sm_cond = lambda M: np.concatenate([(t @ M - M @ t).ravel() for t in sm_gens])  # noqa: E731
# the (1,2,2) component: su(4) Casimir zero
cas4 = np.array([coord(sum(cp(t) @ (cp(t) @ M - M @ cp(t)) - (cp(t) @ M - M @ cp(t)) @ cp(t) for t in colour)) for M in mats]).T
w4, V4 = np.linalg.eigh((cas4 + cas4.T) / 2)
one22 = V4[:, np.abs(w4) < 1e-6]
Akn = np.array([sm_cond((Bv.T @ one22[:, k]).reshape(n, n)) for k in range(one22.shape[1])]).T
ns = one22 @ nullsp(Akn)
Qh = cp(-ks.Jf @ tQ)
Qh = (Qh + Qh.T) / 2
wq, Vq = np.linalg.eigh(Qh)
sect = lambda q: Vq[:, np.abs(wq - q) < 1e-6]  # noqa: E731
Pde = sum(s @ s.T for s in (sect(-1 / 3), sect(-1.0)))
Pun = sum(s @ s.T for s in (sect(2 / 3), sect(0.0)))


def typ(cond):
    Bm = np.array([cond((Bv.T @ ns[:, k]).reshape(n, n)).ravel() for k in range(ns.shape[1])]).T
    return ns @ nullsp(Bm)


up = typ(lambda M: np.concatenate([(Pde @ M).ravel(), (M @ Pde).ravel()]))
dn = typ(lambda M: np.concatenate([(Pun @ M).ravel(), (M @ Pun).ravel()]))
PS21 = [cp(t) for t in ks.PS512]
G23 = PS21 + [ks.Jv, ks.Jv @ ks.Gv]


def stab(Phi, gens):
    A = np.array([(t @ Phi - Phi @ t).ravel() for t in gens])
    return len(gens) - np.linalg.matrix_rank(A, tol=1e-7)


def norm32(K):
    return K * np.sqrt(32) / np.linalg.norm(K)


Ku = norm32((Bv.T @ up[:, 0]).reshape(n, n))
rows = {}
for phi in np.linspace(0, 2 * np.pi, 25)[:-1]:
    Kd = norm32((Bv.T @ (np.cos(phi) * dn[:, 0] + np.sin(phi) * dn[:, 1])).reshape(n, n))
    for lab, (vu, vd) in {"generic": (0.7, 0.4), "equal": (0.5, 0.5), "vd=0": (0.7, 0.0), "vu=0": (0.0, 0.4), "zero": (0.0, 0.0)}.items():
        Phi = vu * Ku + vd * Kd
        s21, s23 = stab(Phi, PS21), stab(Phi, G23)
        rows.setdefault(lab, set()).add(tuple(int(x) for x in (s21, s23, 21 - s21, 23 - s23, (23 - s23) - (21 - s21))))
expected = {"generic": {(16, 17, 5, 6, 1)}, "equal": {(18, 19, 3, 4, 1)}, "vd=0": {(16, 18, 5, 5, 0)}, "vu=0": {(16, 18, 5, 5, 0)},
            "zero": {(21, 23, 0, 0, 0)}}
ok2 = (one22.shape[1] == 8 and ns.shape[1] == 4 and up.shape[1] == 2 and dn.shape[1] == 2 and rows == expected)
report(2, ok2, f"v_R = 0 strata of thm:axion over 24 rotations of K_1d (|K|^2 = 32): (stab21, stab23, massive, orbit, Goldstone) = "
               f"{ {k: sorted(v) for k, v in rows.items()} }; the (1,2,2) component has 8 real dimensions, 4 Standard-Model neutral, 2 up- and "
               f"2 down-type (numerical ranks, tolerance 1e-7)")

# ---------------------------------------------------------------- a11-3: the stratum v_R != 0, v_u = v_d = 0 (high-scale condensate)
P0 = ks.Phi0
s21, s23 = stab(P0, PS21), stab(P0, G23)
ok3 = (s21, s23) == (12, 14)
report(3, ok3, f"thm:axion, stratum v_R != 0, v_u = v_d = 0 (Phi_0 alone): stabiliser {s21} of the 21 Pati-Salam generators, {s23} of the 23; "
               f"orbits {21 - s21} and {23 - s23}: nine massive vectors and no uneaten direction (numerical ranks, tolerance 1e-7)")

sys.exit(0 if not FAILS and NPASS[0] > 0 else 1)
