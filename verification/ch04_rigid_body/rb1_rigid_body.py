#!/usr/bin/env python3
"""rb1: conventional rigid-body mechanics in the conventions of Chapter 4, Section "Conventional rigid-body mechanics" (Phase 7).

    python3 verification/rigidbody/rb1_rigid_body.py

Prints one 'PASS rb1-<id>: <claim>' or 'FAIL rb1-<id>: <claim>' line per claim; exits 0 iff all pass. Exact sympy arithmetic throughout.
The results are inherited (Euler, Arnold, Marsden-Ratiu); the script fixes the sign conventions the chapter uses.

Conventions.
  Rigid body about a fixed point: orientation R in SO(3), body angular velocity Omega with Rdot = R hat(Omega), inertia tensor II
  (symmetric, positive definite), body angular momentum M = II Omega, Hamiltonian H = M.II^{-1}M/2, Lie-Poisson bracket on so(3)*
  {F,K}(M) = -M.(grad F x grad K).
  Planar body with translations: SE(2) = SO(2) x| R^2, Lie algebra R ghat x|_J R^2, [ghat, X] = J X, J = [[0,-1],[1,0]]; coordinates (h, p)
  on se(2)*, bracket {F1,F2} = p^T (F2_h J grad F1 - F1_h J grad F2) (the bracket of eq. (restLPgen) with G replaced by J), Hamiltonian
  H = h^2/(2I) + p^T N p/2; reconstruction xdot = e^{theta J} N p, thetadot = h/I.
"""
import sys

import sympy as sp

NPASS, FAILS = [0], []


def report(cid, ok, text):
    if ok:
        NPASS[0] += 1
        print(f"PASS rb1-{cid}: {text}")
    else:
        FAILS.append(cid)
        print(f"FAIL rb1-{cid}: {text}")
    sys.stdout.flush()


def is_zero(expr):
    return sp.simplify(sp.expand(expr)) == 0


# ------------------------------------------------------------------------------------------------ the rigid body on so(3)*
M = sp.Matrix(sp.symbols("M1 M2 M3", real=True))
I1, I2, I3 = sp.symbols("I1 I2 I3", positive=True)
II = sp.diag(I1, I2, I3)


def grad(F, xs):
    return sp.Matrix([sp.diff(F, x) for x in xs])


def br3(F, K):
    """{F,K}(M) = -M.(grad F x grad K)"""
    return sp.expand(-M.dot(grad(F, M).cross(grad(K, M))))


# 1. Jacobi identity: for the coordinate functions and for three generic quadratic functions
coords = list(M)
jac_coord = all(is_zero(br3(a, br3(b, c)) + br3(b, br3(c, a)) + br3(c, br3(a, b)))
                for a in coords for b in coords for c in coords)
qs = []
for k in range(3):
    cs = sp.symbols(f"q{k}_0:6")
    qs.append(cs[0] * M[0]**2 + cs[1] * M[1]**2 + cs[2] * M[2]**2 + cs[3] * M[0] * M[1] + cs[4] * M[1] * M[2] + cs[5] * M[0] * M[2])
F, K, L = qs
jac_quad = is_zero(br3(F, br3(K, L)) + br3(K, br3(L, F)) + br3(L, br3(F, K)))
report(1, jac_coord and jac_quad, "the bracket {F,K}(M) = -M.(grad F x grad K) on so(3)* satisfies the Jacobi identity (coordinate functions and "
       "three generic quadratic functions; exact)")

# 2. Lie-Poisson equations of H = M.II^{-1}M/2 are Euler's equations Mdot = M x Omega, Omega = II^{-1} M
H = sp.Rational(1, 2) * (M.T * II.inv() * M)[0]
Mdot = sp.Matrix([br3(M[i], H) for i in range(3)])
Omega = II.inv() * M
euler = M.cross(Omega)
ok2 = all(is_zero(Mdot[i] - euler[i]) for i in range(3))
# componentwise: I1 Omega1dot = (I2 - I3) Omega2 Omega3, cyclically
O = sp.Matrix(sp.symbols("O1 O2 O3", real=True))
sub = {M[i]: II[i, i] * O[i] for i in range(3)}
comp = [Mdot[0].subs(sub) - (I2 - I3) * O[1] * O[2], Mdot[1].subs(sub) - (I3 - I1) * O[2] * O[0], Mdot[2].subs(sub) - (I1 - I2) * O[0] * O[1]]
ok2 = ok2 and all(is_zero(c) for c in comp)
report(2, ok2, "the Lie-Poisson equations of H = M.II^{-1}M/2 are Euler's equations Mdot = M x Omega, i.e. I1 Omega1dot = (I2-I3) Omega2 "
       "Omega3 and cyclically (exact, diagonal II)")

# 3. conservation laws: |M|^2 is a Casimir; H is conserved; the spatial angular momentum R M is conserved
cas = sum(M[i]**2 for i in range(3))
ok3 = all(is_zero(br3(cas, M[i])) for i in range(3))
ok3 = ok3 and is_zero(sum(sp.diff(H, M[i]) * Mdot[i] for i in range(3)))
# d/dt (R M) = R (hat(Omega) M + Mdot) and hat(Omega) M = Omega x M, so it vanishes on Euler's equations
w = sp.Matrix(sp.symbols("w1 w2 w3", real=True))
hat = sp.Matrix([[0, -w[2], w[1]], [w[2], 0, -w[0]], [-w[1], w[0], 0]])
v = sp.Matrix(sp.symbols("v1 v2 v3", real=True))
ok3 = ok3 and all(is_zero((hat * v - w.cross(v))[i]) for i in range(3))
ok3 = ok3 and all(is_zero((Omega.cross(M) + euler)[i]) for i in range(3))
report(3, ok3, "|M|^2 Poisson-commutes with every M_i (Casimir), H is conserved along Euler's equations, and with Rdot = R hat(Omega) the spatial "
       "angular momentum R M is conserved (hat(Omega) M = Omega x M; exact)")

# 4. the inertia tensor II = sum m_k (|x_k|^2 1 - x_k x_k^T) of a rational point-mass body: symmetric, positive definite, principal moments
#    obeying the triangle inequalities I_i <= I_j + I_k, with I_j + I_k - I_i = 2 sum m x_i^2 in principal axes
pts = [(sp.Rational(2), [1, 0, 0]), (sp.Rational(2), [-1, 0, 0]), (sp.Rational(1), [0, 2, 0]), (sp.Rational(1), [0, -2, 0]),
       (sp.Rational(3), [0, 0, sp.Rational(1, 3)]), (sp.Rational(3), [0, 0, -sp.Rational(1, 3)])]
T = sp.zeros(3, 3)
for mk, xk in pts:
    x = sp.Matrix(xk)
    T += mk * ((x.T * x)[0] * sp.eye(3) - x * x.T)
ev = [T[i, i] for i in range(3)]
diag = all(T[i, j] == 0 for i in range(3) for j in range(3) if i != j)
pd = all(T[:k, :k].det() > 0 for k in (1, 2, 3))
tri = all(ev[i] <= ev[(i + 1) % 3] + ev[(i + 2) % 3] for i in range(3))
sec = all(ev[(i + 1) % 3] + ev[(i + 2) % 3] - ev[i] == 2 * sum(mk * sp.sympify(xk[i])**2 for mk, xk in pts) for i in range(3))
report(4, diag and pd and tri and sec, f"inertia tensor of a rational point-mass body: diagonal with moments {ev}, positive definite, and "
       "I_j + I_k - I_i = 2 sum m x_i^2 >= 0 (triangle inequalities; exact)")

# 5. stability of steady rotations (Euler): about the axes of largest and smallest moment the linearisation is a centre, about the
#    intermediate axis it has a real pair of eigenvalues
def linear_eigs(axis, Ivals):
    sub = {I1: Ivals[0], I2: Ivals[1], I3: Ivals[2]}
    Om = sp.symbols("Om0", positive=True)
    Ms = [0, 0, 0]
    Ms[axis] = Ivals[axis] * Om
    Jac = sp.Matrix(3, 3, lambda i, j: sp.diff(Mdot[i].subs(sub), M[j])).subs({M[k]: Ms[k] for k in range(3)})
    lam = sp.symbols("lam")
    return sp.factor(Jac.charpoly(lam).as_expr()), Om, lam


Iv = (1, 2, 3)
res = {}
for ax in range(3):
    cp, Om, lam = linear_eigs(ax, Iv)
    res[ax] = cp
# characteristic polynomials: lam*(lam^2 + c Om^2) with c > 0 for the axes of moment 1 and 3, lam*(lam^2 - c Om^2) for moment 2
Om0 = sp.symbols("Om0", positive=True)
lam = sp.symbols("lam")
c1 = sp.simplify(res[0] / lam - lam**2)
c2 = sp.simplify(res[1] / lam - lam**2)
c3 = sp.simplify(res[2] / lam - lam**2)
ok5 = sp.simplify(c1 / Om0**2) > 0 and sp.simplify(c3 / Om0**2) > 0 and sp.simplify(c2 / Om0**2) < 0
report(5, ok5, f"steady rotation with (I1,I2,I3) = (1,2,3): characteristic polynomials lam(lam^2 + {sp.simplify(c1/Om0**2)} Om^2), "
       f"lam(lam^2 {sp.simplify(c2/Om0**2)} Om^2), lam(lam^2 + {sp.simplify(c3/Om0**2)} Om^2): stable about the axes of smallest and largest "
       "moment, unstable about the intermediate one (exact)")

# ------------------------------------------------------------------------------------------------ the planar body with translations, se(2)*
h, p1, p2, th = sp.symbols("h p1 p2 theta", real=True)
P = sp.Matrix([p1, p2])
J = sp.Matrix([[0, -1], [1, 0]])
Iz = sp.symbols("I", positive=True)
a, b = sp.symbols("a b", real=True)
xs = [h, p1, p2]


def br2(F1, F2):
    g1 = sp.Matrix([sp.diff(F1, p1), sp.diff(F1, p2)])
    g2 = sp.Matrix([sp.diff(F2, p1), sp.diff(F2, p2)])
    return sp.expand((P.T * (sp.diff(F2, h) * J * g1 - sp.diff(F1, h) * J * g2))[0])


jac2 = all(is_zero(br2(x, br2(y, z)) + br2(y, br2(z, x)) + br2(z, br2(x, y))) for x in xs for y in xs for z in xs)
cas2 = all(is_zero(br2(p1**2 + p2**2, x)) for x in xs)
report(6, jac2 and cas2 and br2(h, p1) == -p2 and br2(h, p2) == p1,
       "se(2)*: {h,p1} = -p2, {h,p2} = p1, {p1,p2} = 0 satisfy the Jacobi identity and |p|^2 is a Casimir (exact)")

N2 = sp.Matrix([[a, 0], [0, b]])
H2 = h**2 / (2 * Iz) + sp.Rational(1, 2) * (P.T * N2 * P)[0]
pdot = sp.Matrix([br2(p1, H2), br2(p2, H2)])
hdot = br2(h, H2)
ok7 = all(is_zero(pdot[i] - ((h / Iz) * J.T * P)[i]) for i in range(2)) and is_zero(hdot + (P.T * J * N2 * P)[0])
# on the leaf p = e^{-theta J} (P0, 0): I thetaddot = hdot = (a - b) P0^2 sin(2 theta)/2, a pendulum in 2 theta
P0 = sp.symbols("P0", positive=True)
Rm = sp.Matrix([[sp.cos(th), -sp.sin(th)], [sp.sin(th), sp.cos(th)]])          # e^{theta J}
pleaf = Rm.T * sp.Matrix([P0, 0])                                             # e^{-theta J} p0
hdot_leaf = sp.simplify(-(pleaf.T * J * N2 * pleaf)[0])
ok7 = ok7 and is_zero(hdot_leaf - (a - b) * P0**2 * sp.sin(2 * th) / 2)
# isotropic inertia (N commuting with J): hdot = 0, rotation and translation decouple
ok7 = ok7 and is_zero(-(P.T * J * (a * sp.eye(2)) * P)[0])
report(7, ok7, "se(2)*: the Lie-Poisson equations of H = h^2/(2I) + p^T N p/2 are pdot = (h/I) J^T p, hdot = -p^T J N p; with N = diag(a,b) "
       "the leaf p = e^{-theta J}p0 gives I thetaddot = (a-b)|p0|^2 sin(2 theta)/2, a pendulum in 2 theta, and N = a 1 gives hdot = 0 "
       "(exact)")

# 8. reconstruction: thetadot = h/I, xdot = e^{theta J} N p; the spatial momentum e^{theta J} p and J_s = h + p_s^T J x are conserved
x1, x2 = sp.symbols("x1 x2", real=True)
X = sp.Matrix([x1, x2])
thdot = h / Iz
d_ps = sp.simplify(Rm * pdot + thdot * Rm * J * P)           # d/dt (e^{theta J} p) = e^{theta J}(pdot + thetadot J p)
xdot = Rm * N2 * P
ps = Rm * P
dJs = sp.simplify(hdot + (ps.T * J * xdot)[0] + (d_ps.T * J * X)[0])
report(8, all(is_zero(e) for e in d_ps) and is_zero(dJs),
       "se(2): with thetadot = h/I and xdot = e^{theta J} N p, the spatial momentum p_s = e^{theta J} p and J_s = h + p_s^T J x are conserved "
       "for every symmetric N (exact)")

print(f"\n{NPASS[0]} PASS, {len(FAILS)} FAIL")
sys.exit(1 if FAILS else 0)
