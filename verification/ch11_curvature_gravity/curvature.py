"""Exact curvature of a metric given in coordinates (sympy), with the conventions of Chapter 11:
Gamma^a_bc = (1/2) g^ad (d_b g_dc + d_c g_db - d_d g_bc); R^a_bcd = d_c Gamma^a_db - d_d Gamma^a_cb + Gamma^a_ce Gamma^e_db - Gamma^a_de Gamma^e_cb;
R_bd = R^a_bad; R = g^bd R_bd; Einstein G_ab = R_ab - (1/2) R g_ab. Signature as given by the metric.
"""
import itertools

import sympy as sp


class Geometry:
    def __init__(self, g, coords, simplify=sp.cancel, ginv=None):
        self.g, self.x, self.n, self.simp = sp.Matrix(g), list(coords), len(coords), simplify
        self.ginv = (sp.Matrix(ginv) if ginv is not None else self.g.inv(method="LU")).applyfunc(simplify)
        n, x, gi = self.n, self.x, self.ginv
        dg = [[[sp.diff(self.g[i, j], x[k]) for k in range(n)] for j in range(n)] for i in range(n)]
        self.Gam = [[[simplify(sum(gi[a, d] * (dg[d][c][b] + dg[d][b][c] - dg[b][c][d]) for d in range(n)) / 2)
                      for c in range(n)] for b in range(n)] for a in range(n)]
        self._riem = {}

    def riemann(self, a, b, c, d):
        key = (a, b, c, d)
        if key not in self._riem:
            G, x, n = self.Gam, self.x, self.n
            v = (sp.diff(G[a][d][b], x[c]) - sp.diff(G[a][c][b], x[d])
                 + sum(G[a][c][e] * G[e][d][b] - G[a][d][e] * G[e][c][b] for e in range(n)))
            self._riem[key] = self.simp(v)
        return self._riem[key]

    def ricci(self):
        if not hasattr(self, "_ric"):
            n = self.n
            self._ric = sp.Matrix(n, n, lambda b, d: self.simp(sum(self.riemann(a, b, a, d) for a in range(n))))
        return self._ric

    def scalar(self):
        Ric = self.ricci()
        return self.simp(sum(self.ginv[b, d] * Ric[b, d] for b in range(self.n) for d in range(self.n)))

    def einstein(self):
        return (self.ricci() - self.scalar() / 2 * self.g).applyfunc(self.simp)


def maxwell(h, A, coords):
    """F_mn = d_m A_n - d_n A_m, F^2 = F_mn F^mn, W_mn = F_mr F_n^r - (1/4) h_mn F^2 for a 4-metric h and potential components A."""
    n = len(coords); hi = sp.Matrix(h).inv()
    F = sp.Matrix(n, n, lambda m, k: sp.diff(A[k], coords[m]) - sp.diff(A[m], coords[k]))
    Fu = hi * F * hi
    F2 = sp.cancel(sum(F[m, k] * Fu[m, k] for m in range(n) for k in range(n)))
    W = sp.Matrix(n, n, lambda m, k: sp.cancel(sum(F[m, r] * (F * hi)[k, r] for r in range(n)) - sp.Rational(1, 4) * h[m, k] * F2))
    return F, F2, W


def frame_metric(h, fibres):
    """Kaluza-Klein form g = h + sum_k c_k (theta^k)^2 with theta^k = dy^k + a^k_mu dx^mu (a^k may depend on all coordinates).
    fibres: list of (c_k, [a^k_0, ..., a^k_3]). Returns (g, ginv) with ginv = M^-1 D^-1 M^-T exactly, M the unit-triangular coframe matrix."""
    h = sp.Matrix(h); m = len(fibres); n = 4 + m
    M = sp.eye(n); D = sp.zeros(n, n); D[:4, :4] = h
    for k, (c, a) in enumerate(fibres):
        for mu in range(4):
            M[4 + k, mu] = a[mu]
        D[4 + k, 4 + k] = c
    Minv = sp.eye(n)
    Minv[4:, :4] = -M[4:, :4]
    Dinv = sp.zeros(n, n); Dinv[:4, :4] = h.inv()
    for k, (c, a) in enumerate(fibres):
        Dinv[4 + k, 4 + k] = 1 / sp.sympify(c)
    g = (M.T * D * M).applyfunc(sp.expand)
    ginv = (Minv * Dinv * Minv.T).applyfunc(sp.expand)
    assert (g * ginv).applyfunc(lambda q: sp.simplify(sp.expand(q))) == sp.eye(n)
    return g, ginv
