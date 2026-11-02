#!/usr/bin/env python3
"""rb2: the identification of the closure with a rigid body, tested (Chapter 4, Section "The identification tested"; Phase 7).

    python3 verification/rigidbody/rb2_es_map.py

Prints one 'PASS rb2-<id>: <claim>' or 'FAIL rb2-<id>: <claim>' line per claim; exits 0 iff all pass. Checks rb2-1 to rb2-14 are exact
(sympy); rb2-15 is numerical (scipy DOP853, tolerances stated in the line).

Conventions (Chapter 4, Section "Closing the force law"): V = R^{3,3}, eta = diag(1,1,1,-1,-1,-1), binding root G e_i = e_{i+3},
G e_{i+3} = -e_i (G^T = -G, G^2 = -1, eta G = -G eta); closure s = R ghat x|_G V with [ghat, X] = G X, [X, Y] = 0; coordinates (h, p) on s*;
bracket {F1,F2} = p^T (F2_h G grad F1 - F1_h G grad F2) (eq. restLPgen); Hamiltonian H = h^2/(2 I0) + p^T N p/2; reconstruction
thetadot = h/I0, xdot = e^{theta G} N p (Theorem "The frame under the closure"). Lorentz-class inertia N = (nu1 eta^{-1} + nu2 G eta^{-1})/m,
w(p) = C - i D with C = p^T eta^{-1} p, D = p^T G eta^{-1} p (eq. wdef). Units: [X] = L, [tau] = T, [m] = M, [I0] = M L^2.
"""
import sys

import numpy as np
import sympy as sp
from scipy.integrate import solve_ivp

NPASS, FAILS = [0], []


def report(cid, ok, text):
    if ok:
        NPASS[0] += 1
        print(f"PASS rb2-{cid}: {text}")
    else:
        FAILS.append(cid)
        print(f"FAIL rb2-{cid}: {text}")
    sys.stdout.flush()


def is_zero(expr):
    return sp.simplify(sp.expand(expr)) == 0


eta = sp.diag(1, 1, 1, -1, -1, -1)
G = sp.zeros(6, 6)
for i in range(3):
    G[i + 3, i] = 1
    G[i, i + 3] = -1
assert G.T == -G and G * G == -sp.eye(6) and eta * G == -G * eta

h = sp.symbols("h", real=True)
p = sp.Matrix(sp.symbols("p1:7", real=True))
I0, m = sp.symbols("I0 m", positive=True)
coords = [h] + list(p)


def grad_p(F):
    return sp.Matrix([sp.diff(F, x) for x in p])


def br(F1, F2):
    return sp.expand((p.T * (sp.diff(F2, h) * G * grad_p(F1) - sp.diff(F1, h) * G * grad_p(F2)))[0])


# 1. Jacobi identity on s* (coordinate functions and three generic functions quadratic in (h, p))
jac = all(is_zero(br(a, br(b, c)) + br(b, br(c, a)) + br(c, br(a, b))) for a in coords for b in coords for c in coords)
rng = np.random.default_rng(7)
quad = []
for k in range(3):
    cs = rng.integers(-3, 4, size=(7, 7))
    quad.append(sp.expand(sum(int(cs[i, j]) * coords[i] * coords[j] for i in range(7) for j in range(i, 7))))
A, B, C = quad
jac = jac and is_zero(br(A, br(B, C)) + br(B, br(C, A)) + br(C, br(A, B)))
ok1 = jac and all(br(h, p[j]) == (G * p)[j] for j in range(6)) and all(br(p[i], p[j]) == 0 for i in range(6) for j in range(6))
report(1, ok1, "the bracket of eq. (restLPgen) on s* has {h,p_j} = (Gp)_j, {p_i,p_j} = 0 and satisfies the Jacobi identity (coordinate "
       "functions and three integer quadratic functions of (h,p); exact)")

# 2. the Lie-Poisson equations of H = h^2/(2 I0) + p^T N p/2, N generic symmetric (21 symbols), are eq. (LProtor)
nsym = sp.symbols("n0:21", real=True)
N = sp.zeros(6, 6)
k = 0
for i in range(6):
    for j in range(i, 6):
        N[i, j] = N[j, i] = nsym[k]
        k += 1
H = h**2 / (2 * I0) + sp.Rational(1, 2) * (p.T * N * p)[0]
pdot = sp.Matrix([br(p[i], H) for i in range(6)])
hdot = br(h, H)
ok2 = all(is_zero(pdot[i] - ((h / I0) * G.T * p)[i]) for i in range(6)) and is_zero(hdot + (p.T * G * N * p)[0])
report(2, ok2, "for every symmetric N the Lie-Poisson equations of H = h^2/(2 I0) + p^T N p/2 are pdot = (h/I0) G^T p and hdot = -p^T G N p, "
       "eq. (LProtor) (exact)")

# 3. ranks of the Poisson tensors at a generic point: s* rank 2 (5 Casimirs), so(3)* rank 2 (1 Casimir), se(2)* rank 2 (1 Casimir)
P7 = sp.Matrix(7, 7, lambda i, j: br(coords[i], coords[j]))
pt = {h: sp.Rational(3, 2), **{p[i]: sp.Rational(v) for i, v in enumerate([1, -2, 3, 5, -1, 4])}}
r_s = P7.subs(pt).rank()
M3 = sp.Matrix(sp.symbols("M1:4", real=True))
P3 = sp.Matrix([[0, -M3[2], M3[1]], [M3[2], 0, -M3[0]], [-M3[1], M3[0], 0]])     # {M_i,M_j} = -eps_ijk M_k
r_so3 = P3.subs({M3[0]: 1, M3[1]: 2, M3[2]: 3}).rank()
q1, q2 = sp.symbols("q1 q2", real=True)
P2 = sp.Matrix([[0, -q2, q1], [q2, 0, 0], [-q1, 0, 0]])                          # (h, p1, p2) of se(2)*, rb1-6
r_se2 = P2.subs({q1: 1, q2: 2}).rank()
# five functionally independent Casimirs: the real and imaginary parts of z_j conj(z_k), z = z(pi), pi = eta^{-1} p
pi = eta.inv() * p
z = [pi[j] + sp.I * pi[j + 3] for j in range(3)]
cas = []
for j in range(3):
    for kk in range(j, 3):
        e = sp.expand(z[j] * sp.conjugate(z[kk]))
        cas += [sp.re(e), sp.im(e)] if j != kk else [sp.re(e)]
casJ = sp.Matrix([[sp.diff(c, x) for x in p] for c in cas]).subs(pt)
ok3 = r_s == 2 and r_so3 == 2 and r_se2 == 2 and casJ.rank() == 5 and all(is_zero(br(c, h)) for c in cas)
report(3, ok3, f"at a generic rational point the Poisson tensor of s* has rank {r_s} (leaves of dimension 2, 7-2 = 5 Casimirs), those of so(3)* "
       f"and se(2)* rank {r_so3} and {r_se2} (one Casimir each); the nine components of z conj(z)^T, z = z(pi), Poisson-commute with h and "
       f"have Jacobian rank {casJ.rank()} (exact)")

# 4. restricted to a pairing plane (h, p1, p4) the bracket is that of se(2)* with J = G restricted to the plane
G14 = sp.Matrix([[G[0, 0], G[0, 3]], [G[3, 0], G[3, 3]]])
ok4 = G14 == sp.Matrix([[0, -1], [1, 0]]) and br(h, p[0]) == -p[3] and br(h, p[3]) == p[0] and br(p[0], p[3]) == 0
report(4, ok4, "on the pairing plane (e1,e4) G acts as J = [[0,-1],[1,0]] and {h,p1} = -p4, {h,p4} = p1, {p1,p4} = 0: the bracket of se(2)* "
       "of rb1-6 (exact)")

# 5. spatial momenta: p_s = e^{theta G} p and J_s = h + p_s^T G x are conserved for every symmetric N
th = sp.symbols("theta", real=True)
E = sp.cos(th) * sp.eye(6) + sp.sin(th) * G                                       # e^{theta G}
x = sp.Matrix(sp.symbols("x1:7", real=True))
thdot = h / I0
d_ps = E * pdot + thdot * E * G * p                                               # d/dtau (e^{theta G} p)
xdot = E * N * p
ps = E * p
dJs = hdot + (ps.T * G * xdot)[0] + (d_ps.T * G * x)[0]
ok5 = all(is_zero(sp.expand_trig(e)) for e in d_ps) and is_zero(sp.expand(dJs))
report(5, ok5, "with thetadot = h/I0 and xdot = e^{theta G} N p the spatial momentum p_s = e^{theta G} p and J_s = h + p_s^T G x are constants "
       "of motion for every symmetric N, the analogue of linear and total angular momentum (exact)")

# 6. duality class [N,G] = 0 gives hdot = 0 identically; Lorentz class gives hdot = Im(nu w(p))/m and p^T N p/2 = Re(nu w(p))/(2m)
eqs = list(N * G - G * N)
solc = sp.solve(eqs, nsym, dict=True)[0]
Nc = N.subs(solc)
free_c = sorted(Nc.free_symbols, key=str)
ok6a = is_zero((p.T * G * Nc * p)[0]) and len(free_c) == 9
nu1, nu2 = sp.symbols("nu1 nu2", real=True)
NL = (nu1 * eta.inv() + nu2 * G * eta.inv()) / m
Cw = (p.T * eta.inv() * p)[0]
Dw = (p.T * G * eta.inv() * p)[0]
nu = nu1 + sp.I * nu2
w = Cw - sp.I * Dw
ok6b = is_zero(-(p.T * G * NL * p)[0] - sp.im(sp.expand(nu * w)) / m) and is_zero(sp.Rational(1, 2) * (p.T * NL * p)[0] - sp.re(sp.expand(nu * w)) / (2 * m))
ok6b = ok6b and NL * G + G * NL == sp.zeros(6, 6) and NL.T == NL
report(6, ok6a and ok6b, f"duality class: the symmetric N commuting with G form a {len(free_c)}-dimensional space and give hdot = -p^T G N p = 0 for "
       "every p (rotation and frame decouple); Lorentz class N = (nu1 eta^-1 + nu2 G eta^-1)/m anticommutes with G and gives "
       "hdot = Im(nu w(p))/m, p^T N p/2 = Re(nu w(p))/(2m) (exact)")

# 7. Lorentz-class inertia is indefinite: characteristic polynomial (lam^2 - |nu|^2/m^2)^3
lam = sp.symbols("lam")
cpoly = sp.factor(NL.charpoly(lam).as_expr())
ok7 = is_zero(NL.charpoly(lam).as_expr() - (lam**2 - (nu1**2 + nu2**2) / m**2)**3)
report(7, ok7, f"the Lorentz-class inertia has characteristic polynomial {cpoly}: eigenvalues +-|nu|/m, three each, so it is indefinite of "
       "signature (3,3) for nu != 0 and no nonzero Lorentz-class inertia is positive definite (exact)")

# 8. H is unbounded below on s*, bounded below on every leaf: min over the leaf circle of p^T N p/2 is -|nu w(p0)|/(2m)
v_neg = None
for vec in NL.subs({nu1: 1, nu2: 0, m: 1}).eigenvects():
    if vec[0] == -1:
        v_neg = vec[2][0]
t = sp.symbols("t", positive=True)
Hline = (t * v_neg).T * NL.subs({nu1: 1, nu2: 0, m: 1}) * (t * v_neg) / 2
unbounded = sp.limit(Hline[0], t, sp.oo) == -sp.oo
# on the leaf p = e^{-theta G} p0: w(p) = e^{2 i theta} w(p0) (Theorem "The momentum field is frozen"), so Re(nu w) >= -|nu w0|
p0 = sp.Matrix([2, -1, 1, 1, 3, -2])
pl = (sp.cos(th) * sp.eye(6) - sp.sin(th) * G) * p0
wl = sp.expand((pl.T * eta.inv() * pl)[0] - sp.I * (pl.T * G * eta.inv() * pl)[0])
w0 = (p0.T * eta.inv() * p0)[0] - sp.I * (p0.T * G * eta.inv() * p0)[0]
ok8 = unbounded and is_zero(sp.expand(sp.expand_complex(wl - sp.exp(2 * sp.I * th) * w0).rewrite(sp.cos)))
nuv = 2 - sp.I
U = sp.re(sp.expand(nuv * sp.exp(2 * sp.I * th) * w0)) / 2
Umin = -sp.Abs(nuv * w0) / 2
crit = sp.solve(sp.diff(sp.expand_complex(U), th), th)
vals = [sp.nsimplify(sp.simplify(U.subs(th, c))) for c in crit]
ok8 = ok8 and any(sp.simplify(vv - Umin) == 0 for vv in vals) and all(sp.N(vv - Umin) >= -1e-30 for vv in vals)
report(8, ok8, f"Lorentz class: H = h^2/(2 I0) + p^T N p/2 is unbounded below on s* (along an eigenvector of eigenvalue -|nu|/m), while on a leaf "
       f"w(p) = e^(2 i theta) w(p0) and the frame energy is bounded below by -|nu w(p0)|/(2m) (example: p0 = {list(p0)}, nu = 2 - i, minimum "
       f"{Umin} attained; exact)")

# 9. the rotor group is the centre U(1) of U(3) on C^3 = R^6: it commutes with every element of gl(3,C) = {A : [A,G] = 0}, fixes no nonzero
#    vector for theta not in 2 pi Z (det(e^{theta G} - 1) = (2 - 2 cos theta)^3), while every rotation of R^3 fixes an axis
gsyms = sp.symbols("g0:36", real=True)
Ag = sp.Matrix(6, 6, gsyms)
solg = sp.solve(list(Ag * G - G * Ag), gsyms, dict=True)[0]
Agc = Ag.subs(solg)
dim_glc = len(Agc.free_symbols)
comm = all(is_zero(e) for e in (E * Agc - Agc * E))
detE = sp.factor(sp.simplify((E - sp.eye(6)).det()))
ok_det = is_zero(detE - (2 - 2 * sp.cos(th))**3)
Sk = sp.Matrix([[0, -3, 2], [3, 0, -1], [-2, 1, 0]])
Rcay = (sp.eye(3) - Sk).inv() * (sp.eye(3) + Sk)                                   # a rational rotation (Cayley transform)
ok9 = dim_glc == 18 and comm and ok_det and Rcay.T * Rcay == sp.eye(3) and Rcay.det() == 1 and (Rcay - sp.eye(3)).det() == 0
report(9, ok9, f"the rotor group e^(theta G) commutes with the {dim_glc}-dimensional algebra gl(3,C) = {{A : [A,G] = 0}} (it is the centre U(1) of "
       f"U(3) on C^3 = R^6), and det(e^(theta G) - 1) = {detE}, so it fixes no nonzero vector unless theta is in 2 pi Z; a rotation of R^3 "
       "(here a rational Cayley rotation) always fixes its axis, det(R - 1) = 0 (exact)")

# 10. e^{theta G} is an isometry of eta only for theta in pi Z: (e^{theta G})^T eta e^{theta G} = cos 2theta eta + sin 2theta eta G
lhs = sp.simplify(E.T * eta * E)
ok10 = all(is_zero(sp.expand_trig(e)) for e in (lhs - (sp.cos(2 * th) * eta + sp.sin(2 * th) * eta * G)))
report(10, ok10, "(e^(theta G))^T eta e^(theta G) = cos(2 theta) eta + sin(2 theta) eta G: the rotor rotates the three pairing planes of the "
       "Euclidean R^6 but is an isometry of eta only at theta in pi Z (restates Proposition cl:norel(ii); exact)")

# 11. the clamped rotor theta = kappa tau solves the closure iff hdot = 0 along it: always in the duality class, and in the Lorentz class iff
#     nu w(p0) = 0 (hdot = Im(nu e^{2 i theta} w(p0))/m vanishes for all theta iff nu w(p0) = 0)
A1, B1 = sp.symbols("A1 B1", real=True)                                            # nu w(p0) = A1 + i B1
hd = sp.im(sp.expand_complex((A1 + sp.I * B1) * (sp.cos(2 * th) + sp.I * sp.sin(2 * th)))) / m
sol11 = sp.solve([hd.subs(th, 0), hd.subs(th, sp.pi / 4)], [A1, B1], dict=True)
ok11 = sol11 == [{A1: 0, B1: 0}] and ok6a
report(11, ok11, "the clamped rotor (theta = kappa tau, h constant) is a solution of the closure iff hdot vanishes along it: always for duality-class "
       "inertia, and for Lorentz-class inertia iff nu w(p0) = 0, i.e. nu = 0 or a null momentum field (exact)")

# 12. units: with [X] = L, [tau] = T, [m] = M, [I0] = M L^2 the dictionary is dimensionally consistent
dim = {"M": sp.Matrix([1, 0, 0]), "L": sp.Matrix([0, 1, 0]), "T": sp.Matrix([0, 0, 1])}
d_p = dim["M"] + dim["L"] - dim["T"]                                               # p = m eta(alpha_+ X)
d_h = dim["M"] + 2 * dim["L"] - dim["T"]                                           # h = I0 thetadot
d_I0 = dim["M"] + 2 * dim["L"]
d_N = -dim["M"]
d_H = dim["M"] + 2 * dim["L"] - 2 * dim["T"]
ok12 = (d_p - dim["T"] == (d_h - d_I0) + d_p) and (d_h - dim["T"] == 2 * d_p + d_N) and (2 * d_h - d_I0 == d_H) and (2 * d_p + d_N == d_H)
ok12 = ok12 and (d_h == dim["M"] + 2 * dim["L"] - dim["T"]) and (d_h - d_I0 == -dim["T"])
report(12, ok12, "units: [p] = M L T^-1, [h] = M L^2 T^-1 (an angular momentum), [I0] = M L^2 (a moment of inertia), [N] = M^-1, [H] = M L^2 T^-2; "
       "both equations of eq. (LProtor) and both terms of H are dimensionally consistent, and h/I0 is a rate (exact)")

# 13. dimension count of the three systems: (configuration group, dual of the algebra, generic leaf, Casimirs)
counts = {"SO(3)": (3, 3, r_so3, 3 - r_so3), "SE(2)": (3, 3, r_se2, 3 - r_se2), "S": (7, 7, r_s, 7 - r_s)}
ok13 = counts == {"SO(3)": (3, 3, 2, 1), "SE(2)": (3, 3, 2, 1), "S": (7, 7, 2, 5)}
report(13, ok13, f"dimension of the configuration group, of the dual Lie algebra, of a generic leaf, and number of Casimirs: {counts}; in all "
       "three the reduced motion on a leaf has one degree of freedom (exact)")

# 14. the closure algebra s is solvable ([s,s] = V, [V,V] = 0) while so(3) is perfect and has no two-dimensional subalgebra; hence every Lie
#     algebra homomorphism so(3) -> s is zero and every one s -> so(3) has image of dimension at most one
def bracket_s(u, v):                     # u = (a, X) = a ghat + X: [u, v] = (0, a G Y - b G X)
    return (0, u[0] * (G * v[1]) - v[0] * (G * u[1]))
basis_s = [(1, sp.zeros(6, 1))] + [(0, sp.eye(6)[:, j]) for j in range(6)]
der1 = sp.Matrix.hstack(*[sp.Matrix([b_[0]] + list(b_[1])) for b_ in (bracket_s(u, v) for u in basis_s for v in basis_s)])
rank1 = der1.rank()
ghat_comp = all(e == 0 for e in der1.row(0))                         # the derived algebra lies in V
derived = [(0, sp.Matrix(der1[1:, j])) for j in range(der1.shape[1])]
der2 = [bracket_s(u, v) for u in derived for v in derived]
second_zero = all(b_[1] == sp.zeros(6, 1) for b_ in der2)
ex = [sp.Matrix([1, 0, 0]), sp.Matrix([0, 1, 0]), sp.Matrix([0, 0, 1])]
der_so3 = sp.Matrix.hstack(*[a_.cross(b_) for a_ in ex for b_ in ex]).rank()
# no 2-dim subalgebra of so(3): for independent a, b the product a x b is orthogonal to both and nonzero, so not in span(a, b)
a_, b_ = sp.Matrix(sp.symbols("a1:4", real=True)), sp.Matrix(sp.symbols("b1:4", real=True))
c_ = a_.cross(b_)
orth = is_zero(c_.dot(a_)) and is_zero(c_.dot(b_)) and is_zero(c_.dot(c_) - (a_.dot(a_) * b_.dot(b_) - a_.dot(b_)**2))
report(14, rank1 == 6 and ghat_comp and second_zero and der_so3 == 3 and orth,
       f"[s,s] = V (rank {rank1}, no ghat component) and [V,V] = 0, so s is solvable; [so(3),so(3)] has dimension {der_so3}, and for independent "
       "a, b the product a x b is nonzero (|a x b|^2 = |a|^2|b|^2 - (a.b)^2) and orthogonal to span(a,b), so so(3) has no two-dimensional "
       "subalgebra: every homomorphism so(3) -> s is zero and every homomorphism s -> so(3) has image of dimension at most one (exact)")

# 15. numerical: the full closure (theta, h, p, x) with Lorentz-class inertia keeps H, the Casimirs, p_s and J_s while h and the frame
#     exchange angular momentum
Gn = np.array(G.tolist(), dtype=float)
etan = np.array(eta.tolist(), dtype=float)
NLn = (1.0 * np.linalg.inv(etan) + 0.5 * Gn @ np.linalg.inv(etan)) / 1.0          # nu = 1 + 0.5 i, m = 1
I0n = 1.0


def rhs(t, y):
    th_, h_ = y[0], y[1]
    p_ = y[2:8]
    Eth = np.cos(th_) * np.eye(6) + np.sin(th_) * Gn
    return np.concatenate([[h_ / I0n, -p_ @ Gn @ NLn @ p_], (h_ / I0n) * Gn.T @ p_, Eth @ NLn @ p_])


p0n = np.array([1.0, 0.3, -0.2, 0.1, 0.4, -0.5])
y0 = np.concatenate([[0.0, 1.1], p0n, np.zeros(6)])
sol = solve_ivp(rhs, (0.0, 200.0), y0, method="DOP853", rtol=1e-12, atol=1e-12, dense_output=False, t_eval=np.linspace(0, 200, 2001))


def invariants(y):
    th_, h_ = y[0], y[1]
    p_, x_ = y[2:8], y[8:14]
    Eth = np.cos(th_) * np.eye(6) + np.sin(th_) * Gn
    zz = (np.linalg.inv(etan) @ p_)[:3] + 1j * (np.linalg.inv(etan) @ p_)[3:]
    ps_ = Eth @ p_
    return dict(H=h_**2 / (2 * I0n) + 0.5 * p_ @ NLn @ p_, p2=p_ @ p_, E0=abs(zz @ zz), ps=ps_, Js=h_ + ps_ @ Gn @ x_)


inv0 = invariants(sol.y[:, 0])
drift = {k: 0.0 for k in inv0}
for j in range(sol.y.shape[1]):
    iv = invariants(sol.y[:, j])
    for k in inv0:
        drift[k] = max(drift[k], float(np.max(np.abs(np.atleast_1d(iv[k] - inv0[k])))))
h_range = float(np.ptp(sol.y[1]))
ok14 = sol.success and max(drift.values()) < 1e-9 and h_range > 0.1
report(15, ok14, "numerical (scipy DOP853, rtol = atol = 1e-12, tau in [0,200], 2001 samples; nu = 1 + i/2, m = I0 = 1, theta(0) = 0, h(0) = 1.1, "
       "p(0) = (1, 0.3, -0.2, 0.1, 0.4, -0.5), x(0) = 0): max |drift| of H, |p|^2, |w(p)| = 2 E0(pi), p_s and J_s = "
       + ", ".join(f"{drift[k]:.1e}" for k in ("H", "p2", "E0", "ps", "Js")) + f", while h varies over a range {h_range:.3f} (< 1e-9 required)")

print(f"\n{NPASS[0]} PASS, {len(FAILS)} FAIL")
sys.exit(1 if FAILS else 0)
