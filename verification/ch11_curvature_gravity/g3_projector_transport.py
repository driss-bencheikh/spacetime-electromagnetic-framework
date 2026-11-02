#!/usr/bin/env python3
"""g3: the Clifford and projector identities of Chapter 11 (Phase 9C).

    python3 verification/ch11_curvature/g3_projector_transport.py

prop:parahermitian with the grading root of prop:twoK, lem:JN, the holonomy on the bound orbits, prop:obstruct, prop:adapted, the adapted
generator, the expansion of Pi J_N Pi, ex:stueckelberg and the trace form on the spin orbit (ex:ghostPi). Exact arithmetic in the 8x8
realisation (clifford8.py) and on R^{3,3}; sqrt(3) is kept exact.
"""
import sys

import sympy as sp

import os  # noqa: E402
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "common"))
from clifford8 import I8, bl, blade_coeffs, comm, e, eta, vec, vec_coords

NPASS, FAILS = [0], []


def report(cid, ok, text):
    if ok:
        NPASS[0] += 1
        print(f"PASS g3-{cid}: {text}")
    else:
        FAILS.append(cid)
        print(f"FAIL g3-{cid}: {text}")
    sys.stdout.flush()


Z8 = sp.zeros(8, 8)
s3 = sp.sqrt(3)
etaM = sp.diag(1, 1, 1, -1, -1, -1)
varpi = e[1] * e[4] + e[2] * e[5] + e[3] * e[6]
Vol = bl(1, 2, 3, 4, 5, 6)

# grading root K(v) = (1/2)[varpi, v] on vectors
Kmat = sp.zeros(6, 6)
for a in range(1, 7):
    col = vec_coords(comm(varpi, e[a]) / 2)
    for b in range(6):
        Kmat[b, a - 1] = col[b]
Kexp = sp.Matrix([[0, 0, 0, -1, 0, 0], [0, 0, 0, 0, -1, 0], [0, 0, 0, 0, 0, -1], [-1, 0, 0, 0, 0, 0], [0, -1, 0, 0, 0, 0], [0, 0, -1, 0, 0, 0]])
ok = Kmat == Kexp and Kmat * Kmat == sp.eye(6) and Kmat.T * etaM * Kmat == -etaM
plus = [sp.Matrix([1 if j == i else (-1 if j == i + 3 else 0) for j in range(6)]) for i in range(3)]
minus = [sp.Matrix([1 if j in (i, i + 3) else 0 for j in range(6)]) for i in range(3)]
ok &= all(Kmat * p == p for p in plus) and all(Kmat * m == -m for m in minus)
ok &= all((p.T * etaM * q)[0] == 0 for p in plus for q in plus) and all((p.T * etaM * q)[0] == 0 for p in minus for q in minus)
Om = etaM * Kmat
ok &= Om.T == -Om
report("1", ok, "prop:parahermitian with prop:twoK: K(v) = (1/2)[varpi, v] is the matrix K e_i = -e_{i+3}, K e_{i+3} = -e_i, K^2 = 1, eta(K.,K.) = -eta; "
       "its +1 and -1 eigenspaces are spanned by e_i - e_{i+3} and e_i + e_{i+3} and are totally null; Omega(X,Y) = eta(KX,Y) is antisymmetric "
       "(exact). The para-Kaehler criterion (nabla K = 0) is inherited (Cruceanu-Fortuny-Gadea)")

# lem:JN
u1 = (e[1] - e[2]) / sp.sqrt(2); u2 = (e[1] + e[2] - 2 * e[3]) / sp.sqrt(6)
JN = u1 * u2
JN_claim = (e[2] * e[3] + e[3] * e[1] + e[1] * e[2]) / s3
W = (e[1] + e[2] + e[3]) / s3
ok = (JN - JN_claim).applyfunc(sp.simplify) == Z8 and (JN * JN + I8).applyfunc(sp.simplify) == Z8
ok &= all(comm(JN, X).applyfunc(sp.simplify) == Z8 for X in (W, e[4], e[5], e[6]))
report("2", ok, "lem:JN: J_N = u1 u2 = (e2e3 + e3e1 + e1e2)/sqrt3, J_N^2 = -1, and J_N commutes with the symmetric clock W = (e1+e2+e3)/sqrt3 and with "
       "e4, e5, e6 (exact)")

# holonomy on the bound orbits (R^{3,3}, 6x6)
G = sp.Matrix([[0, 0, 0, -1, 0, 0], [0, 0, 0, 0, -1, 0], [0, 0, 0, 0, 0, -1], [1, 0, 0, 0, 0, 0], [0, 1, 0, 0, 0, 0], [0, 0, 1, 0, 0, 0]])
U1 = sp.Matrix([1, -1, 0, 0, 0, 0]) / sp.sqrt(2); U2 = sp.Matrix([1, 1, -2, 0, 0, 0]) / sp.sqrt(6)
RN = U2 * U1.T - U1 * U2.T                                          # R_N u1 = u2, R_N u2 = -u1, zero on N-perp (all in V^+, eta = +1 there)
Phi = sp.symbols("Phi", real=True)
hol = sp.eye(6) + sp.sin(Phi) * RN + (1 - sp.cos(Phi)) * RN * RN
ok = (RN ** 3 + RN).applyfunc(sp.simplify) == sp.zeros(6, 6) and (RN.T * etaM + etaM * RN) == sp.zeros(6, 6)
ok &= (hol - (RN * Phi).exp()).applyfunc(lambda x: sp.simplify(x.rewrite(sp.cos))) == sp.zeros(6, 6)
D = (hol * G * hol.inv() - G).applyfunc(sp.simplify)
A1 = sp.Matrix(6, 6, lambda i, j: sp.diff(D[i, j], sp.sin(Phi)) if False else 0)
# write D = sin(Phi) A1 + (1 - cos(Phi)) A2 by evaluating at two angles
A1 = D.subs(Phi, sp.pi / 2) - D.subs(Phi, sp.pi) / 2
A2 = D.subs(Phi, sp.pi) / 2
ok &= (D - (sp.sin(Phi) * A1 + (1 - sp.cos(Phi)) * A2)).applyfunc(sp.simplify) == sp.zeros(6, 6)
ok &= sp.Matrix.hstack(sp.Matrix(list(A1)), sp.Matrix(list(A2))).rank() == 2 and comm(RN, G) != sp.zeros(6, 6)
kap, tau = sp.symbols("kappa tau", real=True)
X0 = sp.Matrix(sp.symbols("x1:7"))
eG = sp.cos(kap * tau) * sp.eye(6) - sp.sin(kap * tau) * G                # e^{-kappa G tau}, G^2 = -1
Xt = eG * X0
Y = hol * Xt
alpha_plus = sp.diff(Y, tau) + kap * G * Y
ok &= (alpha_plus - kap * (G - hol * G * hol.inv()) * hol * Xt).applyfunc(sp.simplify) == sp.zeros(6, 1)
report("3", ok, "holonomy: R_N (rotation of the normal plane N = <u1,u2>) is eta-antisymmetric, R_N^3 = -R_N, exp(Phi R_N) = 1 + sin(Phi) R_N + (1-cos(Phi)) R_N^2, "
       "and hGh^-1 - G = sin(Phi) A1 + (1-cos(Phi)) A2 with A1, A2 linearly independent, so hGh^-1 = G iff Phi in 2 pi Z; R_N does not commute with G; for the "
       "zero-momentum orbit X = e^{-kappa G tau} X0, alpha_+(hX) = kappa (G - hGh^-1) hX (exact)")

# prop:obstruct and the projector
Pi = I8 / 2 - (varpi + Vol) / 4
ok = Pi * Pi == Pi and Pi.T == Pi and Pi.rank() == 4
CJ = comm(JN, Pi).applyfunc(sp.simplify)
ok &= CJ.rank() == 4
ok &= comm(e[4] * e[5], Pi) == (e[2] * e[4] - e[1] * e[5]) / 2
ok &= all(comm(Vol, e[a] * e[b]) == Z8 for a in range(1, 7) for b in range(a + 1, 7))
report("4", ok, "prop:obstruct: Pi = 1/2 - (varpi + V)/4 is an orthogonal projector of rank 4; rank [J_N, Pi] = 4; [e4e5, Pi] = (e2e4 - e1e5)/2; V = e1...e6 "
       "commutes with every bivector (exact). With prop:kk (R_{uv mu nu} = F_{mu nu} on Sigma) the spin curvature on Sigma is (1/2) F J_N, so D Pi = 0 is "
       "impossible where F != 0; the tangent-level statement (J_N fixes E pointwise) is g3-2")

# prop:adapted: [[Pi, delta], Pi] = -delta for off-diagonal delta; the correction
d = sp.Matrix(8, 8, lambda i, j: sp.Rational((3 * i + 5 * j) % 7 - 3, 1 + (i * j) % 4))
offd = Pi * d * (I8 - Pi) + (I8 - Pi) * d * Pi
ok = comm(comm(Pi, offd), Pi) == -offd
Cadj = comm(Pi, offd)
ok &= comm(Cadj, Pi) == -offd                                       # D-hat Pi = delta + [C, Pi] = 0
Kv = Kmat
Xr = sp.Matrix(6, 6, lambda i, j: sp.Rational((2 * i + 3 * j) % 5 - 2, 1 + (i + j) % 3))
Xr = Xr - etaM * Xr.T * etaM                                         # eta-antisymmetric
Dk = (Xr - Kv * Xr * Kv) / 2                                         # anticommutes with K: a tangent vector to the involutions
ok &= Kv * Dk + Dk * Kv == sp.zeros(6, 6) and Kv * Dk * Kv == -Dk
corr = Kv * Dk / 2
ok &= Dk + comm(corr, Kv) == sp.zeros(6, 6) and (etaM * corr).T == -(etaM * corr)
report("5", ok, "prop:adapted: for an off-diagonal delta (Pi delta Pi = 0 = (1-Pi) delta (1-Pi)), [[Pi,delta],Pi] = -delta, so C = [Pi, D Pi] gives D-hat Pi = 0; "
       "on tangent vectors, for K^2 = 1 and D K + K D = 0, K D K = -D and the correction (1/2) K D is eta-antisymmetric and satisfies D + [(1/2)KD, K] = 0 "
       "(exact, on explicit rational test matrices; the identities are linear in delta, D)")

# the adapted generator and Pi J_N Pi
Jm = (e[5] * e[6] + e[6] * e[4] + e[4] * e[5]) / s3
ok = comm(JN - Jm, Pi).applyfunc(sp.simplify) == Z8
ok &= (comm((JN - Jm) / 2, e[4]) - (e[6] - e[5]) / s3).applyfunc(sp.simplify) == Z8
Wc = (e[1] + e[2] + e[3]) / s3
biv_model = [Wc * e[4], Wc * e[5], Wc * e[6], e[4] * e[5], e[5] * e[6], e[6] * e[4]]
ok &= all(comm(X, Pi).applyfunc(sp.simplify) != Z8 for X in biv_model)
axial = [(a, b) for a, b in [(1, 4), (1, 5), (1, 6), (4, 5), (5, 6), (4, 6)] if comm(e[a] * e[b], Pi) == Z8]
ok &= axial == [(1, 4)]
coef = blade_coeffs((4 * s3 * Pi * JN * Pi).applyfunc(sp.simplify))
want = {(1, 2): 1, (1, 3): -1, (2, 3): 1, (4, 5): -1, (4, 6): 1, (5, 6): -1, (1, 2, 3, 4): -1, (1, 2, 3, 5): -1, (1, 2, 3, 6): -1,
        (1, 4, 5, 6): 1, (2, 4, 5, 6): 1, (3, 4, 5, 6): 1}
ok &= coef == want
report("6", ok, "adapted generator: [J_+ - J_-, Pi] = 0 with J_- = (e5e6 + e6e4 + e4e5)/sqrt3, but [(J_+ - J_-)/2, e4] = (e6 - e5)/sqrt3; all six bivectors of the "
       "model's four-plane <W, e4, e5, e6> fail to commute with Pi, whereas for the axial clock e1 the boost e1e4, a pole operator, commutes with Pi (the other "
       "five do not); 4 sqrt3 Pi J_N Pi = e1e2 - e1e3 + e2e3 - e4e5 + e4e6 - e5e6 - (e1e2e3e4 + e1e2e3e5 + e1e2e3e6) + (e1e4e5e6 + e2e4e5e6 + "
       "e3e4e5e6), with grade-four terms of the same size as the bivector terms (exact blade expansion)")

# ex:stueckelberg and ex:ghostPi
B = JN / 2
T = comm(B, Pi).applyfunc(sp.simplify)
ok = T.T == T and sp.simplify((T * T).trace()) == 1
ok &= comm(Pi, comm(B, comm(B, Pi))).applyfunc(sp.simplify) == Z8
tang = [comm(e[a] * e[b], Pi) for a in range(1, 7) for b in range(a + 1, 7)]
Gram = sp.Matrix(15, 15, lambda i, j: (tang[i] * tang[j]).trace())
ev = Gram.eigenvals()
pos = sum(m for v, m in ev.items() if v > 0); neg = sum(m for v, m in ev.items() if v < 0)
ok &= (pos, neg) == (3, 3)
report("7", ok, f"ex:stueckelberg: T = [J_N/2, Pi] is symmetric with Tr(T^2) = 1, and [Pi,[B,[B,Pi]]] = 0 for B = J_N/2, so the one-angle ansatz is consistent; "
       f"ex:ghostPi: the trace form Tr([X,Pi][Y,Pi]) on the spin orbit of Pi has signature ({pos},{neg}), the compact directions positive and the boosts negative "
       "(exact). A negative direction of this quadratic form is a ghost of the projector kinetic term only once that term is adopted")

print(f"{NPASS[0]}/{NPASS[0] + len(FAILS)} checks passed")
sys.exit(0 if not FAILS else 1)
