#!/usr/bin/env python3
"""v2: the kernel space and the one-loop second variation on all of it (Phase 9B, 2 November 2026).

    python3 verification/vacuum/v2_kernel_hessian.py

Tests, on the 272-dimensional kernel space built in kernel_space.py: thm:kernelspace (dimension, components, types, Phi_0), the
normalisation E_x(u_0) = 2x - 99 of the leading-logarithm energy E_x = 8xQ - (3/16)T of sec:fullK, thm:fullsecond (i)-(iii) with every entry
of tab:271, the thresholds 27, 28.5, 315/11, 30, the agreement of the (10,1,3) block with the tensor model of v1, prop:LR and the great-circle
formula. Floating-point linear algebra; tolerances printed. One PASS/FAIL line per check; exit 0 iff all pass.
"""
import sys

import numpy as np

import kernel_space as k

NPASS, FAILS = [0], []


def report(cid, ok, text):
    if ok:
        NPASS[0] += 1
        print(f"PASS v2-{cid}: {text}")
    else:
        FAILS.append(cid)
        print(f"FAIL v2-{cid}: {text}")
    sys.stdout.flush()


dK = k.dimK
ok = dK == 272 and k.gap_null > 1 and k.A_resid < 1e-10
report("1", ok, f"def:kernel: the antisymmetric 64x64 matrices on V_* that commute with the three interior rotations and anticommute with the three "
       f"boosts form a space of dimension {dK} (null space of the Gram matrix of the 6 conditions; next eigenvalue {k.gap_null:.3f}); the 21 "
       f"Pati-Salam generators preserve it (residual {k.A_resid:.1e}) (float)")

# components by Casimirs
A = k.A
CL = sum(A[a] @ A[a] for a in range(0, 3)); CR = sum(A[a] @ A[a] for a in range(3, 6)); C4 = sum(A[a] @ A[a] for a in range(6, 21))
Hc = C4 + np.pi * CL + np.e * CR
ev, U = np.linalg.eigh(Hc)
comps = {}
for val in np.unique(np.round(ev, 6)):
    cols = U[:, np.abs(ev - val) < 1e-5]
    lab = tuple(round(float(np.trace(cols.T @ C @ cols)) / cols.shape[1], 6) for C in (C4, CL, CR))
    comps.setdefault(lab, []).append(cols)
comps = {lab: np.hstack(v) for lab, v in comps.items()}
dims = {lab: v.shape[1] for lab, v in comps.items()}
want = {(-36.0, -2.0, 0.0): 60, (-36.0, 0.0, -2.0): 60, (-32.0, -0.75, -0.75): 120, (-20.0, 0.0, 0.0): 24, (0.0, -0.75, -0.75): 8}
ok = dims == want
Jv, Gv = k.Jv, k.Gv
typ = {lab: (np.allclose(Jv @ k.mat(v[:, 0]), -k.mat(v[:, 0]) @ Jv, atol=1e-9), np.allclose(Gv @ k.mat(v[:, 0]), k.mat(v[:, 0]) @ Gv, atol=1e-9))
       for lab, v in comps.items()}
ok &= all(typ[l] == (True, True) for l in [(-36.0, -2.0, 0.0), (-36.0, 0.0, -2.0), (-20.0, 0.0, 0.0)])
ok &= all(typ[l] == (False, False) for l in [(-32.0, -0.75, -0.75), (0.0, -0.75, -0.75)])
report("2", ok, f"thm:kernelspace(ii),(iii): Casimir labels (su4, su2_L, su2_R) and real dimensions {dims}, as stated; the (10,.,.) and 2(6,1,1) components "
       "are of Majorana type (anticommute with J, commute with Gamma), (15,2,2) and (1,2,2) of Dirac type (float, tol 1e-5)")
P1013 = comps[(-36.0, 0.0, -2.0)]; P1031 = comps[(-36.0, -2.0, 0.0)]; P611 = comps[(-20.0, 0.0, 0.0)]

c0 = k.c0
ok = np.linalg.norm(k.mat(c0) - k.u0) < 1e-12 and np.linalg.norm(P1013 @ (P1013.T @ c0) - c0) < 1e-12
ok &= np.allclose(k.Phi0.T @ k.Phi0, k.Pnu, atol=1e-12) and abs(np.trace(k.Pnu) - 4) < 1e-12
Qv, Tv = k.Q(c0), k.Tq(c0)
ok &= abs(Qv - 0.25) < 1e-12 and abs(Tv - 528) < 1e-9
report("3", ok, f"Phi_0 of eq:Phi0 lies in the kernel space and in (10,1,3), Phi_0^T Phi_0 = P_nu of rank 4; at u_0 = Phi_0/|Phi_0|, Q = {Qv:.12f} = 1/4 and "
       f"T = {Tv:.9f} = 528, so E_x(u_0) = 8x/4 - (3/16)528 = 2x - 99 (float, tol 1e-9)")

# second variation on the tangent space
HQ, HT = k.hessQ(c0), k.hessT(c0)
x0grad = k.grad(c0, 0.0)
gradQ = (k.grad(c0, 1.0) - x0grad) / 8
radial = np.linalg.norm(k.grad(c0, 10.0) - (c0 @ k.grad(c0, 10.0)) * c0)


def H_of(x):
    Hfull = 8 * x * HQ - 3 / 16 * HT
    gfull = 8 * x * gradQ + x0grad
    return k.sphere_hessian(Hfull, gfull, c0)


H0, Tb = H_of(0.0)
H1, _ = H_of(1.0)
HF = H1 - H0
HG = H0
# Pi_nu: kernels pairing nu_R with a Weyl fermion: (1 - P_nu) M (1 - P_nu) = 0
Qn = np.eye(64) - k.Pnu
img = np.einsum("ab,ibc,cd->iad", Qn, k.B, Qn).reshape(dK, -1)
gg_ev, gg_U = np.linalg.eigh(img @ img.T)                                  # kernels M with (1 - P_nu) M (1 - P_nu) = 0
nulnu = gg_U[:, gg_ev < 1e-9].T
Pi_nu = nulnu.T @ nulnu
Pi_nu_t = Tb.T @ Pi_nu @ Tb
ok = radial < 1e-9 and nulnu.shape[0] == 32 and np.allclose(HF, -8 * (np.eye(dK - 1) - Pi_nu_t), atol=1e-9)
report("4", ok, f"thm:fullsecond(i): u_0 is a critical point on the sphere (tangential gradient {radial:.1e}); H(x) = H_G + x H_F with H_F = -8(1 - Pi_nu) on the "
       f"271-dimensional tangent space, Pi_nu the projector onto the kernels pairing nu_R with a fermion of the block, of rank {nulnu.shape[0]} on the kernel "
       "space (31 on the tangent space) (float, tol 1e-9)")
# H_G preserves the components
blocks = [Tb.T @ v for v in comps.values()]
offd = max(np.linalg.norm(b1.T @ HG @ b2) for i, b1 in enumerate(blocks) for j, b2 in enumerate(blocks) if i < j)
ok = offd < 1e-8
report("5", ok, f"thm:fullsecond(i): H_G preserves each Pati-Salam component of the kernel space (largest off-diagonal block norm {offd:.1e}) (float)")


def table271(x):
    r = [(0, 9), (120 - 4 * x, 2), (144 - 4 * x, 6), (144 - 4 * x, 12), (156 - 4 * x, 12), (120 - 4 * x, 12),
         (132 - 2 * x - 2 * np.sqrt(x * x + 576), 6), (132 - 2 * x + 2 * np.sqrt(x * x + 576), 6), (180 - 4 * x, 6), (180 - 4 * x, 12),
         (108 - 4 * x, 6), (156 - 4 * x, 18), (180 - 4 * x, 36), (174 - 4 * x, 64), (162 - 4 * x, 24), (114 - 4 * x, 12), (114, 12),
         (138 - 4 * x, 4), (186 - 4 * x, 4), (162 - 2 * x - 2 * np.sqrt(x * x - 12 * x + 144), 4), (162 - 2 * x + 2 * np.sqrt(x * x - 12 * x + 144), 4)]
    assert sum(m for _, m in r) == 271
    return np.sort(np.concatenate([[v] * m for v, m in r]))


devs = {}
for x in (0.0, 10.0, 24.0, 27.0, 29.0):
    Hx, _ = H_of(x)
    m2 = np.sort(np.linalg.eigvalsh(Hx)) / 2
    devs[x] = float(np.abs(m2 - table271(x)).max())
ok = max(devs.values()) < 1e-8
report("6", ok, f"thm:fullsecond(ii), tab:271: m^2 = (1/2) eigenvalues of H(x) agree with every entry of tab:271 (271 modes, with the mixing matrices M_3 and "
       f"M_2) at x = 0, 10, 24, 27, 29; max deviation {max(devs.values()):.1e} (float, tol 1e-8)")


def low(x):
    Hx, _ = H_of(x)
    return np.sort(np.linalg.eigvalsh(Hx))


l26, l27, l275 = low(26.999), low(27.0), low(27.001)
ok = np.sum(np.abs(l26) < 1e-7) == 9 and l26.min() > -1e-9 and np.sum(np.abs(l27) < 1e-7) == 15 and np.sum(l275 < -1e-6) == 6
ok &= abs(low(26.0)[9] - 2 * 4.0) < 1e-7
report("7", ok, "thm:fullsecond(iii): H(x) >= 0 with kernel exactly the nine Goldstone directions for x < 27 (at x = 26.999: 9 zero modes, none negative); "
       "at x = 27 six more modes become massless and for x > 27 they are negative (6 at x = 27.001), the l_L l_L triplet 108 - 4x (float, tol 1e-7)")

x3 = 315 / 11
m3low = 132 - 2 * x3 - 2 * np.sqrt(x3 * x3 + 576)
ok = abs(m3low) < 1e-12 and abs(114 - 4 * 28.5) < 1e-12 and abs(120 - 4 * 30) < 1e-12 and abs(108 - 4 * 27) < 1e-12
ok &= abs(84 / 2 - 42) < 1e-12
report("8", ok, "the thresholds of tab:271 are 27 (l_L l_L), 28.5 (l_L-d_R), 315/11 = 28.636... (lower eigenvalue of M_3, which moves the u_Re_R threshold "
       "of the tensor model from 42 to 315/11) and 30 (e_Re_R, d_Rd_R) (exact arithmetic on the table entries, which v2-6 ties to H(x))")

# the (10,1,3) block agrees with the tensor model of v1
_u, _s, _ = np.linalg.svd(Tb.T @ P1013, full_matrices=False)
Tb1013 = _u[:, _s > 0.5]                                                     # tangent part of (10,1,3): 59 dims (u_0 itself lies in (10,1,3))
devT = 0.0
for x in (0.0, 5.0, 17.5, 29.0):
    Hx, _ = H_of(x)
    m2 = np.sort(np.linalg.eigvalsh(Tb1013.T @ Hx @ Tb1013)) / 2
    tab = np.sort(np.concatenate([[0] * 9, [120 - 4 * x] * 2, [84 - 2 * x] * 6, [144 - 4 * x] * 6, [144 - 4 * x] * 12, [156 - 4 * x] * 12,
                                  [120 - 4 * x] * 12]))
    devT = max(devT, np.abs(m2 - tab).max())
ok = Tb1013.shape[1] == 59 and devT < 1e-8
report("9", ok, f"the second variation restricted to (10,1,3) (59 tangent directions) reproduces the tensor-model spectrum of thm:oneloop1013 "
       f"(v1-13: d_Rd_R = 120 - 4x, u_Re_R = 84 - 2x) at x = 0, 5, 17.5, 29, max deviation {devT:.1e}; the full-space normalisation (module metric, "
       "T, Q) and the tensor-model normalisation (Tr T_aT_b = 1/2, M_V, Tr(S^dag S)^2) agree (float)")

# prop:LR and the great circle
O = k.cp((k.E[1] @ k.E[7] @ k.E[8]).astype(float))
Phi0p = O @ k.Phi0 @ O.T
cp0 = k.coords(Phi0p / np.linalg.norm(Phi0p))
ok = np.allclose(O @ O.T, np.eye(64), atol=1e-12) and np.linalg.norm(k.mat(cp0) - Phi0p / np.linalg.norm(Phi0p)) < 1e-12
ok &= abs(cp0 @ c0) < 1e-12 and np.linalg.norm(P1031 @ (P1031.T @ cp0) - cp0) < 1e-12
devc = 0.0
for th in np.linspace(0, np.pi / 2, 7):
    c = np.cos(th) * c0 + np.sin(th) * cp0
    s2 = np.sin(2 * th) ** 2
    devc = max(devc, abs(k.Q(c) - (0.25 - s2 / 8)), max(abs(k.Ex(c, x) - (2 * x - 99 + (27 - x) * s2)) for x in (0.0, 13.0, 27.0, 31.0)))
ok &= devc < 1e-9
report("10", ok, f"prop:LR: O = E1E7E8 restricted to V_* is orthogonal and maps Phi_0 to a kernel in (10,3,1) orthogonal to Phi_0; along the great circle "
       f"Q = 1/4 - sin^2(2 theta)/8 and E_x = 2x - 99 + (27 - x) sin^2(2 theta) (7 angles, x = 0, 13, 27, 31; max deviation {devc:.1e}) (float)")

print(f"{NPASS[0]}/{NPASS[0] + len(FAILS)} checks passed")
sys.exit(0 if not FAILS else 1)
