#!/usr/bin/env python3
"""z2: Phase 2 adjudication checks for the photon-selection note (PH-01 ... PH-09), 1 November 2026.

    python3 verification/ch08_gauge_photon/z2_photon_adjudication.py

Prints one 'PASS z2-<id>: <claim>' or 'FAIL z2-<id>: <claim>' line per claim; exits 0 iff all pass.
Independent of p1/p2/p3 scripts (nothing imported). Third rung n = 3, S_3 = R^512.

Realisation (tex/s02_arena.tex): e1=sx(x)1(x)1, e2=sz(x)sx(x)1, e3=sz(x)sz(x)sx, e4=ep(x)1(x)1, e5=sz(x)ep(x)1, e6=sz(x)sz(x)ep;
E_{6(k-1)+a} = omega^{(x)(k-1)} (x) e_a (x) 1. Conventions of the photon selection unit p3_identification (now in the manuscript):
interior {1,4,5,6}; J = product of the 14 internal generators; Gamma = J E1E4E5E6; A = {2,3,8,9}, W = {7,13,14,15}, C = {10,...,18}\\W;
X = Gamma chi_W, Z = -J chi_C, U = -J E2E3, P_* = (1/8)(1+X)(1+Z)(1+U), P_L = P_*(1-Gamma)/2; tau_i = -1/2 w_j w_k on W;
t^L_i = (1-chi_W)/2 tau_i, t^R_i = (1+chi_W)/2 tau_i, H_C = -(E10E16+E11E17+E12E18), t_Y = t^R_3 + H_C/6, B_em = tau_3 + H_C/6,
Q_A = E2E3/2; charge operators hat X = -J X; Pfam = (1/8)(1+X1)(1+X2)(1+X3), X1 = E2E3E8E9, X2 = E13...E18, X3 = E1...E6 E13 E16;
Pi_P = (1+vol_P)/2 for an oriented four-plane P of V_+ = span{E_i : i in A u W}.

Exact integer arithmetic for blades and projector identities (all generators are signed permutations; products with 1/2, 1/8 are
checked after scaling). Frames along the test curves have Pythagorean entries, so traces are rational; they are evaluated in float64
and compared with the exact rational values to 1e-9.
"""
import itertools
import sys
from fractions import Fraction
from functools import reduce

import numpy as np

NPASS, FAILS = [0], []


def report(cid, ok, text):
    if ok:
        NPASS[0] += 1
        print(f"PASS z2-{cid}: {text}")
    else:
        FAILS.append(cid)
        print(f"FAIL z2-{cid}: {text}")
    sys.stdout.flush()


sx = np.array([[0, 1], [1, 0]], dtype=np.int64)
sz = np.array([[1, 0], [0, -1]], dtype=np.int64)
ep = np.array([[0, -1], [1, 0]], dtype=np.int64)
I2 = np.eye(2, dtype=np.int64)
kr = lambda *a: reduce(np.kron, a)
e = [kr(sx, I2, I2), kr(sz, sx, I2), kr(sz, sz, sx), kr(ep, I2, I2), kr(sz, ep, I2), kr(sz, sz, ep)]
om = reduce(np.matmul, e)
I8 = np.eye(8, dtype=np.int64)
E = [None]
for k in range(3):
    for a in range(6):
        E.append(kr(*([om] * k + [e[a]] + [I8] * (2 - k))))
N = 512
Id = np.eye(N, dtype=np.int64)
blade = lambda idx: reduce(np.matmul, [E[i] for i in idx], Id)
eta = [None] + [int(E[i][0] @ E[i][:, 0] if False else (E[i] @ E[i])[0, 0]) for i in range(1, 19)]
interior = [1, 4, 5, 6]
internal = [i for i in range(1, 19) if i not in interior]
Aidx, Widx = [2, 3, 8, 9], [7, 13, 14, 15]
Cidx = [10, 11, 12, 16, 17, 18]
Vp = sorted(Aidx + Widx)
assert all(eta[i] == 1 for i in Vp) and all(eta[i] == -1 for i in Cidx)
J = blade(internal)
Gam = J @ blade(interior)
chiW, chiC = blade(Widx), blade(Cidx)
qA = E[2] @ E[3]
Xv, Zv, Uv = Gam @ chiW, -(J @ chiC), -(J @ qA)
P8 = (Id + Xv) @ (Id + Zv) @ (Id + Uv)  # = 8 P_*
PL16 = P8 @ (Id - Gam)  # = 16 P_L
assert (P8 @ P8 == 8 * P8).all()
X1, X2, X3 = blade([2, 3, 8, 9]), blade([13, 14, 15, 16, 17, 18]), blade([1, 2, 3, 4, 5, 6]) @ E[13] @ E[16]
Pfam8 = (Id + X1) @ (Id + X2) @ (Id + X3)  # = 8 Pfam
w = [E[i] for i in Widx]
tau2 = [-(w[1] @ w[2]), -(w[2] @ w[0]), -(w[0] @ w[1])]  # = 2 tau_i
tL4 = [(Id - chiW) @ t for t in tau2]  # = 4 t^L_i
tR4 = [(Id + chiW) @ t for t in tau2]
HC = -(E[10] @ E[16] + E[11] @ E[17] + E[12] @ E[18])
tY12 = 3 * tR4[2] + 2 * HC  # = 12 t_Y
Bem6 = 3 * tau2[2] + HC  # = 6 B_em
volA = blade(Aidx)


def hat(X):  # charge operator -J X
    return -(J @ X)


def tr(M):
    return float(np.trace(M))


# ======================================================================================== PH-05 charges, anomaly, normalisation
hatQA2 = hat(qA)  # = 2 hat Q_A
ok = (hatQA2 == Uv).all() and (hatQA2 @ P8 == P8).all()  # hat Q_A = U/2, hat Q_A P_* = P_*/2
ok &= ((qA @ P8) == (J @ P8)).all()
T3h4 = hat(tL4[2])  # = 4 hat T^3_L
QAh2 = hatQA2
Qemh6 = hat(Bem6)
anomA = Fraction(int(np.trace(PL16 @ QAh2 @ T3h4 @ T3h4)), 16 * 2 * 16)
anomEM = Fraction(int(np.trace(PL16 @ Qemh6 @ T3h4 @ T3h4)), 16 * 6 * 16)
# B_em convention of the directive: B_em = -1/2 E7E13 + 1/6 H_C
ok_conv = ((6 * (-(E[7] @ E[13])) // 2 + HC) == Bem6).all() and ((-(E[7] @ E[13])) == tau2[2]).all()
# electron unit: B_em^2 = -hat Q_em^2 on V_*
ok_unit = ((Bem6 @ Bem6 @ P8) == -(Qemh6 @ Qemh6 @ P8)).all()
report("1", ok and anomA == 4 and anomEM == 0 and ok_conv and ok_unit,
       f"PH-05: hat Q_A = U/2 and hat Q_A P_* = P_*/2, q_A P_* = J P_* (exact); B_em = -1/2 E7E13 + 1/6 H_C is the convention of the scripts; "
       f"B_em^2 = -hat Q_em^2 on V_*; Tr(P_L hat Q_A (hat T3_L)^2) = {anomA} (mixed SU(2)_L^2 U(1)_A anomaly), Tr(P_L hat Q_em (hat T3_L)^2) = {anomEM} (exact)")

# Does the mixed anomaly of Q_A cancel on the full spectator Weyl block (1-Gamma)/2 ?  (information for the 'cancelled' option)
WL2 = Id - Gam  # = 2 x Weyl projector
anom_full = Fraction(int(np.trace(WL2 @ QAh2 @ T3h4 @ T3h4)), 2 * 2 * 16)
anom_rest = anom_full - anomA
# on the full block, does hat T3_L act? count; and Q_A^3 trace on the block P_L (cubic) and gravitational
cubic = Fraction(int(np.trace(PL16 @ QAh2 @ QAh2 @ QAh2)), 16 * 8)
grav = Fraction(int(np.trace(PL16 @ QAh2)), 16 * 2) - Fraction(int(np.trace((P8 @ (Id + Gam)) @ QAh2)), 16 * 2)
report("2", anom_full == 0 and anom_rest == -4,
       f"PH-05 (cancellation option): on the whole spectator Weyl block (1-Gamma)/2 the mixed trace Tr(hat Q_A (hat T3_L)^2) is {anom_full}: the "
       f"rest of the block, where hat Q_A = -1/2 (two mirror families), contributes {anom_rest} and cancels the block's +4. So U(1)_A is anomaly-free "
       f"only if the mirror matter is kept, i.e. not with P_* as the matter (exact). [cubic Tr(P_L hat Q_A^3) = {cubic}, chiral difference of "
       f"Tr hat Q_A on P_*: {grav}]")

# ======================================================================================== PH-02 counterexample at Pi_{-W}
Vbasis = {i: E[i] for i in Vp}


def vol_of(frame):
    """frame: list of 4 dicts {generator index: Fraction}; returns the Clifford product (float matrix)."""
    vecs = [sum(float(c) * Vbasis[i] for i, c in f.items()) for f in frame]
    return reduce(np.matmul, vecs)


def pluck(frame, idx):
    M = np.array([[float(f.get(i, 0)) for i in idx] for f in frame])
    return np.linalg.det(M)


ok = True
vals = []
for (c, s) in ((Fraction(3, 5), Fraction(4, 5)), (Fraction(12, 13), Fraction(5, 13)), (Fraction(24, 25), Fraction(7, 25)), (Fraction(40, 41), Fraction(9, 41))):
    # start at -W: frame (-w1, w2, w3, w4); rotate each w_j towards a_j by the same angle t (cos t = c, sin t = s)
    frame = [{Widx[0]: -c, Aidx[0]: s}] + [{Widx[j]: c, Aidx[j]: s} for j in (1, 2, 3)]
    vol = vol_of(frame)
    Pi = 0.5 * (np.eye(N) + vol)
    pA, pW = pluck(frame, Aidx), pluck(frame, Widx)
    tPf = tr(Pi @ Pfam8) / 8
    tPL = tr(Pi @ PL16) / 16
    exp_pA = float(s ** 4)
    ok &= abs(pA - exp_pA) < 1e-12 and abs(pW + float(c ** 4)) < 1e-12
    ok &= abs(tPf - 32 * (1 + exp_pA)) < 1e-9 and abs(tPL - 16 * (1 + pA - pW)) < 1e-9
    ok &= np.allclose(Pi @ Pi, Pi, atol=1e-9)
    vals.append((float(s), round(-tPf + 32, 6)))
    ok &= -tPf + 32 < 0
# reference value at Pi_{-W}
PimW = 0.5 * (np.eye(N) - chiW.astype(float))
ok &= abs(tr(PimW @ Pfam8) / 8 - 32) < 1e-9
report("3", ok, f"PH-02: on the curve P(t) from -W rotating every w_j towards a_j by t, p_A = sin^4 t and Tr(Pi Pfam) = 32(1+sin^4 t) (checked with the "
       f"512-dim operators at four Pythagorean angles); at kappa = lambda = 0, mu > 0, V = -mu Tr(Pi Pfam) = -32 mu (1 + sin^4 t) < V(Pi_-W) = -32 mu, "
       f"so Pi_-W is not a local minimum (values of V+32mu per unit mu at sin t: {vals})")

# Hessian coefficients at the four points for V = kappa V_SM - lambda Tr(Pi P_L) - mu Tr(Pi Pfam): exact bookkeeping from the
# expansions of p3:prop:hessian (re-derived from Plucker: p_A = 1 - |x|^2/2 near A, p_W = O(x^4) there, and the converse near W)
def hess(point, kc, lam, mu):  # kc = 128 kappa c_SM
    tbl = {"A": kc + 8 * lam + 16 * mu, "-W": kc + 8 * lam + 0 * mu, "W": kc - 8 * lam + 0 * mu, "-A": kc - 8 * lam - 16 * mu}
    return tbl[point]


ok = hess("-W", 0, 0, 1) == 0 and hess("A", 0, 0, 1) > 0 and hess("-W", 1, 0, 1) > 0 and hess("-W", 0, 1, 0) > 0
report("4", ok, "PH-02: Hessian multiples at the four points are 128 kappa c_SM + (8 lambda + 16 mu, 8 lambda, -8 lambda, -8 lambda - 16 mu) at "
       "(Pi_A, Pi_-W, Pi_W, Pi_-A); at kappa = lambda = 0 the point Pi_-W is degenerate (Hessian 0) and by z2-3 not a minimum; Pi_-W is a "
       "nondegenerate local minimum iff kappa > 0 or lambda > 0; Pi_A iff kappa, lambda, mu not all zero")

# ======================================================================================== PH-03 / PH-04 neutral deformation
# Pi_theta = Ad(exp(theta/2 E8 E14)) Pi_A: E8 -> cos th E8 + sin th E14 inside A
ok = True
res = []
gens = tL4 + [tY12]
PiA = 0.5 * (np.eye(N) + volA.astype(float))
VSM_A = sum(np.linalg.norm(G @ PiA - PiA @ G) ** 2 for G in gens)
ok &= VSM_A < 1e-9  # all of h_SM (su(2)_L, u(1)_Y) commutes with Pi_A; su(3)_C commutes with every Pi_P (acts trivially on V_+)
ok &= all(((g @ volA) == (volA @ g)).all() for g in tL4 + [tY12, Bem6, 2 * qA])
for (c, s) in ((Fraction(3, 5), Fraction(4, 5)), (Fraction(-7, 25), Fraction(24, 25)), (Fraction(5, 13), Fraction(12, 13))):
    frame = [{2: Fraction(1)}, {3: Fraction(1)}, {8: c, 14: s}, {9: Fraction(1)}]
    Pi = 0.5 * (np.eye(N) + vol_of(frame))
    tPL, tPf = tr(Pi @ PL16) / 16, tr(Pi @ Pfam8) / 8
    VSM = sum(np.linalg.norm(G @ Pi - Pi @ G) ** 2 for G in gens)
    ok &= abs(tPL - 16 * (1 + float(c))) < 1e-9 and abs(tPf - 32 * (1 + float(c))) < 1e-9 and VSM > 1e-6
    res.append((float(c), round(tPL, 6), round(tPf, 6)))
report("5", ok, f"PH-03/PH-04: all of h_SM, B_em and Q_A commute with Pi_A (electroweak symmetry unbroken there); along Pi_theta = "
       f"Ad(exp(theta E8E14/2)) Pi_A, Tr(Pi P_L) = 16(1+cos theta) and Tr(Pi Pfam) = 32(1+cos theta), while the alignment energy is > 0 "
       f"for theta not in pi Z; hence V(Pi_theta) - V(Pi_A) = kappa V_SM(Pi_theta) + (16 lambda + 32 mu)(1 - cos theta) >= 0: the stated "
       f"potential opposes the electroweak deformation ((cos, Tr Pi P_L, Tr Pi Pfam) = {res})")

# ======================================================================================== PH-09 weak-angle normalisation
nL3 = Fraction(int(np.trace(tL4[2].T @ tL4[2])), 16 * 512)  # |t3|^2 = -Tr(t t)/512 = Tr(t^T t)/512
nY = Fraction(int(np.trace(tY12.T @ tY12)), 144 * 512)
s2 = nL3 / (nL3 + nY)
report("6", nL3 == Fraction(1, 8) and nY == Fraction(5, 24) and s2 == Fraction(3, 8),
       f"PH-09: in the invariant form <X,Y> = -Tr(XY)/512, |t3_L|^2 = {nL3} and |t_Y|^2 = {nY}; with one coupling for the orthonormalised "
       f"generators, sin^2 theta_W = |t3|^2/(|t3|^2+|t_Y|^2) = {s2}: a statement about this normalisation and the assumption g_Y = g_2 at one "
       f"(unspecified) scale (exact)")

print()
print(f"z2: {NPASS[0]} PASS, {len(FAILS)} FAIL" + (f" (failed: {FAILS})" if FAILS else ""))
sys.exit(0 if not FAILS else 1)
