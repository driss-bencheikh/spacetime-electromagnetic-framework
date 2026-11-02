"""Shared construction for the Phase 9B vacuum checks: the chiral block V_* = im P_* of Cl(9,9) on R^512, the 272-dimensional space of
Lorentz-invariant mass kernels (def:kernel), the 21 Pati-Salam generators, the module metric, Phi_0, P_nu, and the radiative energy
E_x = 8x Q - (3/16) T of Section 7.3.6 (sec:fullK) with its gradient and Hessian.

The realisation is that of Appendix A (and of verification/core/z3_core.py, which this module re-implements without importing it):
e1=sx(x)1(x)1, e2=sz(x)sx(x)1, e3=sz(x)sz(x)sx, e4=ep(x)1(x)1, e5=sz(x)ep(x)1, e6=sz(x)sz(x)ep, ep=[[0,-1],[1,0]]; the 18 generators of the
third rung are om^k (x) e_a (x) 1 with om = e1...e6. Everything below is float64 linear algebra on exactly specified integer matrices.
"""
import itertools
from functools import reduce

import numpy as np

sx = np.array([[0, 1], [1, 0]]); sz = np.array([[1, 0], [0, -1]]); ep = np.array([[0, -1], [1, 0]]); I2 = np.eye(2, dtype=int)
kr = lambda *a: reduce(np.kron, a)
e = [None, kr(sx, I2, I2), kr(sz, sx, I2), kr(sz, sz, sx), kr(ep, I2, I2), kr(sz, ep, I2), kr(sz, sz, ep)]
I8 = np.eye(8, dtype=int)
om = reduce(np.matmul, e[1:])
E = [None] + [kr(*([om] * k + [e[a]] + [I8] * (2 - k))) for k in range(3) for a in range(1, 7)]
I512 = np.eye(512, dtype=int)
bl = lambda S: reduce(np.matmul, [E[i] for i in S], I512)
interior = [1, 4, 5, 6]
internal = [i for i in range(1, 19) if i not in interior]
Widx, Cidx = [7, 13, 14, 15], [10, 11, 12, 16, 17, 18]
J = bl(internal); Gam = J @ bl(interior)
chiW, chiC, qA = bl(Widx), bl(Cidx), E[2] @ E[3]
Xi, Zi, Ui = Gam @ chiW, -(J @ chiC), -(J @ qA)
Pst = ((I512 + Xi) @ (I512 + Zi) @ (I512 + Ui)).astype(float) / 8
w4 = [E[i] for i in Widx]
tau = [-(w4[1] @ w4[2]) / 2, -(w4[2] @ w4[0]) / 2, -(w4[0] @ w4[1]) / 2]
tL = [((I512 - chiW) @ t) / 2 for t in tau]
tR = [((I512 + chiW) @ t) / 2 for t in tau]
colour = [(E[i] @ E[j]).astype(float) for i, j in itertools.combinations(Cidx, 2)]
PS512 = [t.astype(float) for t in tL + tR] + colour                       # 3 + 3 + 15
HC = -(E[10] @ E[16] + E[11] @ E[17] + E[12] @ E[18]).astype(float)
Jf, Gf = J.astype(float), Gam.astype(float)
Id = np.eye(512)
TR3h = -Jf @ tR[2]
BmL = -Jf @ HC / 3
Pnu512 = Pst @ ((Id + Gf) / 2) @ (2 * TR3h @ (TR3h + Id / 2)) @ ((Id / 3 - BmL) * 0.75)
Pnu512 = Pst @ Pnu512 @ Pst
Phi0_512 = Pnu512 @ Pst @ (E[1] @ bl([3, 9, 13, 14, 16, 17, 18])).astype(float) @ Pst @ Pnu512
Gmod = np.array([[np.trace(a.T @ b) / 512 for b in PS512] for a in PS512])

# orthonormal basis of V_* and compressions
_ev, _U = np.linalg.eigh(Pst)
W = _U[:, _ev > 0.5]                                                        # 512 x 64
cp = lambda X: W.T @ X @ W
rot = [cp((E[4] @ E[5]).astype(float)), cp((E[5] @ E[6]).astype(float)), cp((E[6] @ E[4]).astype(float))]
boo = [cp((E[1] @ E[a]).astype(float)) for a in (4, 5, 6)]
tPS = [cp(t) for t in PS512]
Jv, Gv, Pnu, Phi0 = cp(Jf), cp(Gf), cp(Pnu512), cp(Phi0_512)
K0 = cp(E[1].astype(float) @ E[1].astype(float))

# ------------------------------------------------------------------ the kernel space: null space of the defining conditions
n = 64
pairs = list(itertools.combinations(range(n), 2))
Eanti = np.zeros((len(pairs), n, n))
for k, (a, b) in enumerate(pairs):
    Eanti[k, a, b] = 1 / np.sqrt(2); Eanti[k, b, a] = -1 / np.sqrt(2)
Gram = np.zeros((len(pairs), len(pairs)))
for Rm in rot:
    Cm = (Eanti @ Rm - Rm @ Eanti).reshape(len(pairs), -1)
    Gram += Cm @ Cm.T
for Bm in boo:
    Cm = (Eanti @ Bm + Bm @ Eanti).reshape(len(pairs), -1)
    Gram += Cm @ Cm.T
gev, gU = np.linalg.eigh(Gram)
NULL_TOL = 1e-9
dimK = int(np.sum(gev < NULL_TOL))
B = np.einsum("ki,kab->iab", gU[:, gev < NULL_TOL], Eanti)                 # orthonormal basis of the kernel space (dimK x 64 x 64)
Bv = B.reshape(dimK, -1)
gap_null = float(gev[dimK]) if dimK < len(gev) else float("inf")


def coords(M):
    return Bv @ M.ravel()


def mat(c):
    return np.einsum("i,iab->ab", c, B)


# adjoint action of the Pati-Salam generators on the kernel space
A = np.array([Bv @ (t @ B - B @ t).reshape(dimK, -1).T for t in tPS])      # 21 x dimK x dimK, A[a][j,i] = <B_j,[t_a,B_i]>
A_resid = max(np.linalg.norm(np.einsum("ji,jab->iab", A[a], B) - (tPS[a] @ B - B @ tPS[a])) for a in range(len(tPS)))
ginv = np.linalg.inv(Gmod)
u0 = Phi0 / np.linalg.norm(Phi0)
c0 = coords(u0)


def Q(c):
    u = mat(c); U = u.T @ u
    return np.sum(U * U)


def Tq(c):
    Z = np.einsum("aji,i->aj", A, c); K = Z @ Z.T
    return np.trace(ginv @ K @ ginv @ K)


def Ex(c, x):
    return 8 * x * Q(c) - 3 / 16 * Tq(c)


def grad(c, x):
    u = mat(c); U = u.T @ u
    gQ = 4 * Bv @ (u @ U).ravel()
    Z = np.einsum("aji,i->aj", A, c); K = Z @ Z.T; N = ginv @ K @ ginv
    gT = 4 * np.einsum("aji,aj->i", A, N @ Z)
    return 8 * x * gQ - 3 / 16 * gT


def hessQ(c):
    u = mat(c); U = u.T @ u
    X = np.einsum("ba,ibc->iac", u, B) + np.einsum("iba,bc->iac", B, u)     # u^T B_i + B_i^T u
    Xv = X.reshape(dimK, -1)
    Mq = Bv @ (B @ U).reshape(dimK, -1).T
    return 2 * Xv @ Xv.T + 2 * (Mq + Mq.T)


def hessT(c):
    Z = np.einsum("aji,i->aj", A, c); K = Z @ Z.T; N = ginv @ K @ ginv
    dK = np.einsum("ak,bki->iab", Z, A); dK = dK + dK.transpose(0, 2, 1)
    gdKg = np.einsum("ab,ibc,cd->iad", ginv, dK, ginv)
    H1 = 2 * gdKg.reshape(dimK, -1) @ dK.reshape(dimK, -1).T
    H2 = 4 * np.einsum("ab,aki,bkj->ij", N, A, A)
    return H1 + H2


def sphere_hessian(Hfull, gfull, c):
    """Riemannian Hessian on the unit sphere at unit c, in an orthonormal basis of the tangent space (dimK-1)."""
    Tb = np.linalg.svd(np.eye(dimK) - np.outer(c, c))[0][:, : dimK - 1]
    return Tb.T @ (Hfull - (c @ gfull) * np.eye(dimK)) @ Tb, Tb
