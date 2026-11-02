"""The 8x8 realisation of Cl(3,3) of Appendix A, exact (sympy), with helpers for Chapters 11 and 12.
e1=sx(x)1(x)1, e2=sz(x)sx(x)1, e3=sz(x)sz(x)sx, e4=ep(x)1(x)1, e5=sz(x)ep(x)1, e6=sz(x)sz(x)ep, ep=[[0,-1],[1,0]]; eta = diag(1,1,1,-1,-1,-1).
"""
import itertools
from functools import reduce

import sympy as sp

sx = sp.Matrix([[0, 1], [1, 0]]); sz = sp.Matrix([[1, 0], [0, -1]]); ep = sp.Matrix([[0, -1], [1, 0]]); I2 = sp.eye(2)


def kron(*ms):
    return reduce(lambda a, b: sp.kronecker_product(a, b), ms)


e = [None, kron(sx, I2, I2), kron(sz, sx, I2), kron(sz, sz, sx), kron(ep, I2, I2), kron(sz, ep, I2), kron(sz, sz, ep)]
I8 = sp.eye(8)
eta = [None, 1, 1, 1, -1, -1, -1]


def bl(*idx):
    return reduce(lambda a, b: a * b, [e[i] for i in idx], I8)


def vec(c):
    """Clifford image of the vector sum c_a e_a (c indexed 1..6 via a list of length 6)."""
    return sum((c[a - 1] * e[a] for a in range(1, 7)), sp.zeros(8, 8))


def vec_coords(X):
    """coefficients of X in the vectors e_1..e_6 (trace inner product; e_a^{-1} = eta_a e_a)."""
    return [sp.nsimplify(sp.simplify((eta[a] * e[a] * X).trace() / 8)) for a in range(1, 7)]


BIV = [(a, b) for a, b in itertools.combinations(range(1, 7), 2)]


def biv(a, b):
    return e[a] * e[b]


def comm(A, B):
    return A * B - B * A


def blade_coeffs(X):
    """coefficients of X in the 64 blades e_S (S increasing), via Tr(e_S^{-1} X)/8."""
    out = {}
    for k in range(0, 7):
        for S in itertools.combinations(range(1, 7), k):
            b = bl(*S) if S else I8
            c = sp.nsimplify(sp.simplify((b.inv() * X).trace() / 8))
            if c != 0:
                out[S] = c
    return out
