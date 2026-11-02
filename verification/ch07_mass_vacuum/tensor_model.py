"""Shared construction of the tensor model T = Sym^2(C^4)(x)Sym^2(C^2) for the Phase 9B minimisations (v3). It re-implements the
construction of v1_tensor_model.py as importable functions, with analytic gradients on the unit sphere of T (60 real coordinates).
"""
import itertools

import numpy as np


def su(n):
    gens = []
    for a, b in itertools.combinations(range(n), 2):
        X = np.zeros((n, n), complex); X[a, b] = X[b, a] = 0.5; gens.append(X)
        Y = np.zeros((n, n), complex); Y[a, b] = -0.5j; Y[b, a] = 0.5j; gens.append(Y)
    for k in range(1, n):
        H = np.zeros((n, n), complex); H[:k, :k] = np.eye(k); H[k, k] = -k
        gens.append(H / np.sqrt(2 * k * (k + 1)))
    return gens


def sym_basis(n):
    B = []
    for a in range(n):
        M = np.zeros((n, n)); M[a, a] = 1; B.append(M)
    for a, b in itertools.combinations(range(n), 2):
        M = np.zeros((n, n)); M[a, b] = M[b, a] = 1 / np.sqrt(2); B.append(M)
    return B


G8 = [np.kron(X, np.eye(2)) for X in su(4)] + [np.kron(np.eye(4), Y) for Y in su(2)]
FB = np.array([np.kron(A, Bm) for A in sym_basis(4) for Bm in sym_basis(2)])
NT = 30
to_mat = lambda s: np.einsum("j,jab->ab", s, FB)
to_vec = lambda S: np.einsum("jab,ab->j", FB.conj(), S)
R = [np.array([[np.trace(F1.conj().T @ (X @ F2 + F2 @ X.T)) for F2 in FB] for F1 in FB]) for X in G8]
I30 = np.eye(NT)
tens = lambda r: np.kron(r, I30) + np.kron(I30, r)
C4 = sum(tens(r) @ tens(r) for r in R[:15]); C2 = sum(tens(r) @ tens(r) for r in R[15:])
SW = np.zeros((NT * NT, NT * NT))
for i in range(NT):
    for j in range(NT):
        SW[i * NT + j, j * NT + i] = 1
Psym = (np.eye(NT * NT) + SW) / 2
_ev, _U = np.linalg.eigh(Psym @ (C4 + np.pi * C2) @ Psym + 1000 * (np.eye(NT * NT) - Psym))
_order = {175: 0, 35: 1, 100: 2, 20: 3, 135: 4}
P = [None] * 5
for val in np.unique(np.round(_ev, 7)):
    if val < 999:
        cols = _U[:, np.abs(_ev - val) < 1e-6]
        P[_order[cols.shape[1]]] = cols @ cols.conj().T
a_rec = np.array([1, -0.5, -0.5, 0.25, 0])
c_vec = np.array([33, 3, -3, 39, 5])
S0 = np.zeros((8, 8), complex); S0[6, 6] = 1.0
s0 = to_vec(S0)
one = to_vec(np.eye(8).astype(complex))
eps2 = np.array([[0, 1], [-1, 0]])
eps4 = np.zeros((4,) * 4)
for p in itertools.permutations(range(4)):
    eps4[p] = np.linalg.det(np.eye(4)[list(p)])


def cplx(xr):
    return xr[:NT] + 1j * xr[NT:]


def Ik(s):
    w = np.kron(s, s)
    return np.array([np.real(w.conj() @ (p @ w)) for p in P])


def grad_Ik(s):
    """Gradients of I_1..I_5 in the 60 real coordinates (Re s, Im s)."""
    w = np.kron(s, s)
    out = []
    for p in P:
        Pm = (p @ w).reshape(NT, NT)
        h = np.conj(Pm).T @ s
        out.append(4 * np.concatenate([np.real(h), -np.imag(h)]))
    return out


def Dq(S):
    S4 = S.reshape(4, 2, 4, 2)
    Pt = np.einsum("su,tv,aset,bufv->abef", eps2, eps2, S4, S4)
    return np.einsum("abcd,efgh,abef,cdgh->", eps4, eps4, Pt, Pt) / 24


def dD(S):
    """Holomorphic derivative dD/dS_{as,et} as a 4x2x4x2 array (D(S + dS) = D(S) + sum G dS + O(dS^2))."""
    S4 = S.reshape(4, 2, 4, 2)
    Pt = np.einsum("su,tv,aset,bufv->abef", eps2, eps2, S4, S4)
    return np.einsum("abcd,efgh,su,tv,bufv,cdgh->aset", eps4, eps4, eps2, eps2, S4, Pt) / 6


def grad_reD(s, kappa):
    """Gradient of 2 Re(kappa D) in the 60 real coordinates."""
    G = dD(to_mat(s)).reshape(8, 8)
    h = np.einsum("jab,ab->j", FB, G)                                          # dD/ds_j
    return 2 * np.concatenate([np.real(kappa * h), -np.imag(kappa * h)])


def energy(xr, ell, kappa=0.0):
    """F(s/|s|) for F = sum ell_k I_k + 2 Re(kappa D) (degree 4), with its gradient in the 60 real coordinates."""
    nrm = np.linalg.norm(xr); u = xr / nrm; s = cplx(u)
    I = Ik(s); F = ell @ I + (2 * np.real(kappa * Dq(to_mat(s))) if kappa else 0.0)
    g = sum(l * gi for l, gi in zip(ell, grad_Ik(s))) + (grad_reD(s, kappa) if kappa else 0.0)
    g = g - (u @ g) * u
    return F, g / nrm
