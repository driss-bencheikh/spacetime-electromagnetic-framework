#!/usr/bin/env python3
"""rb3: from rigid-body rotation to planetary orbital dynamics (Chapter 4, Section "From Rigid-Body Rotation to Planetary Orbital
Dynamics"; Phase 7).

    python3 verification/rigidbody/rb3_orbital.py [--write-data]

Prints one 'PASS rb3-<id>: <claim>' or 'FAIL rb3-<id>: <claim>' line per claim; exits 0 iff all pass. rb3-1 to rb3-8 are exact (sympy);
rb3-9 is a quadrature (scipy quad, tolerance stated); rb3-10 to rb3-13 are numerical integrations (scipy DOP853, tolerances stated).
With --write-data the script also writes the data of the spin-rate figure to manuscript/figures/ch04/spinrate_{so,es}.dat.

Benchmarks (assumed, not derived from the arena): the Newtonian two-body problem with potential -G M1 M2 / r; relative motion
rddot = -mu r/|r|^3, mu = G(M1+M2); the planar spin-orbit problem of a rigid triaxial body (principal moments A < B < C, spin axis normal to
the orbital plane) in the test-body limit, C thetaddot = -(3/2)(G M/r^3)(B - A) sin 2(theta - f), f the true anomaly.
Arena side (Chapter 4): in a pairing plane (e1,e4) the split-form interior with the binding root, X = x e1 + y e4, z = x + i y,
L = (m/2)[(xdot^2 - ydot^2) - kappa^2 (x^2 - y^2)] - m kappa d(xy)/dtau, solutions z = c_- e^{-i kappa tau} + c_+ e^{i kappa tau};
the complex invariant q(z) = z.z (= z^2 in one plane); the duality pendulum I0 thetaddot = (R0/m) sin(2 theta + phi0) on a leaf of the
closure. Units of the numerical runs: dimensionless, G M = a = 1, so that the mean motion is n = 1 and the orbital period 2 pi.
"""
import pathlib
import sys

import numpy as np
import sympy as sp
from scipy.integrate import quad, solve_ivp

ROOT = pathlib.Path(__file__).resolve().parents[2]
WRITE = "--write-data" in sys.argv
NPASS, FAILS = [0], []


def report(cid, ok, text):
    if ok:
        NPASS[0] += 1
        print(f"PASS rb3-{cid}: {text}")
    else:
        FAILS.append(cid)
        print(f"FAIL rb3-{cid}: {text}")
    sys.stdout.flush()


def is_zero(expr):
    return sp.simplify(sp.expand(expr)) == 0


tau, kap, m, s = sp.symbols("tau kappa m s", positive=True)
x, y = sp.Function("x")(tau), sp.Function("y")(tau)

# 1. the pairing-plane interior is the isotropic Hooke oscillator; xydot - yxdot is conserved as the Noether charge of the boost of the mixed
#    plane, while the Euclidean rotation, which maps solutions to solutions, is not a symmetry of the Lagrangian
L = m / 2 * ((x.diff(tau)**2 - y.diff(tau)**2) - kap**2 * (x**2 - y**2)) - m * kap * (x * y).diff(tau)
el = [sp.simplify(sp.diff(L, q) - sp.diff(sp.diff(L, q.diff(tau)), tau)) for q in (x, y)]
hooke = is_zero(el[0] + m * (x.diff(tau, 2) + kap**2 * x)) and is_zero(el[1] - m * (y.diff(tau, 2) + kap**2 * y))
# variation of L under a vector field (dx, dy) depending on (x, y): delta L = sum dL/dq dq + dL/dqdot dqdot
def delta(Lag, dx, dy):
    return sp.expand(sp.diff(Lag, x) * dx + sp.diff(Lag, y) * dy + sp.diff(Lag, x.diff(tau)) * dx.diff(tau) + sp.diff(Lag, y.diff(tau)) * dy.diff(tau))


dL_boost = delta(L, y, x)
K_boost = -m * kap * (x**2 + y**2)                     # delta L = dK/dtau
boost_sym = is_zero(dL_boost - K_boost.diff(tau))
Q_boost = sp.expand(sp.diff(L, x.diff(tau)) * y + sp.diff(L, y.diff(tau)) * x - K_boost)
noether = is_zero(Q_boost - m * (x.diff(tau) * y - y.diff(tau) * x))
dL_rot = delta(L, -y, x)
# a total derivative has identically vanishing Euler-Lagrange expression; the rotation's delta L does not
el_rot = sp.simplify(sp.diff(dL_rot, x) - sp.diff(sp.diff(dL_rot, x.diff(tau)), tau))
report(1, hooke and boost_sym and noether and not is_zero(el_rot),
       "in a pairing plane the split-form interior gives xddot = -kappa^2 x, yddot = -kappa^2 y (isotropic Hooke law); the boost "
       "(dx,dy) = (y,x) of the mixed plane changes L by a total derivative and has Noether charge m(xdot y - ydot x); the Euclidean rotation "
       "(dx,dy) = (-y,x) does not (its delta L has Euler-Lagrange expression " + str(el_rot) + "; exact)")

# 2. isochrony against Kepler's third law
c1, c2, a1, b1, a2, b2 = sp.symbols("c1 c2 a1 b1 a2 b2", real=True)
cm, cp = a1 + sp.I * b1, a2 + sp.I * b2
z = cm * sp.exp(-sp.I * kap * tau) + cp * sp.exp(sp.I * kap * tau)
period = is_zero(sp.expand(z.subs(tau, tau + 2 * sp.pi / kap) - z))
half = is_zero(sp.expand(z.subs(tau, tau + sp.pi / kap) + z))
amu, aa, lamb = sp.symbols("mu a lambda", positive=True)
TK = 2 * sp.pi * sp.sqrt(aa**3 / amu)
report(2, period and half and sp.simplify(TK.subs(aa, lamb * aa) / TK) == lamb**sp.Rational(3, 2),
       "every interior orbit has period 2 pi/kappa whatever its size and satisfies z(tau + pi/kappa) = -z(tau), so the ellipse is centred on the "
       "force centre; a Kepler orbit has period 2 pi sqrt(a^3/mu), multiplied by lambda^(3/2) when the orbit is scaled by lambda (exact)")

# 3. Bohlin: q = z^2 has foci 0 and 4 c_- c_+ and major axis 2(|c_-|^2 + |c_+|^2): |q| + |q - 4 c_- c_+| = 2(|c_-|^2 + |c_+|^2)
ph = sp.symbols("phi", real=True)
zz = cm * sp.exp(-sp.I * ph) + cp * sp.exp(sp.I * ph)
ident1 = is_zero(sp.expand(zz**2 - 4 * cm * cp - (cm * sp.exp(-sp.I * ph) - cp * sp.exp(sp.I * ph))**2))
mod2 = lambda e: sp.expand(e * sp.conjugate(e))
sumdist = sp.simplify(sp.expand_complex(mod2(zz) + mod2(cm * sp.exp(-sp.I * ph) - cp * sp.exp(sp.I * ph))))
ident2 = is_zero(sumdist - 2 * (a1**2 + b1**2 + a2**2 + b2**2))
report(3, ident1 and ident2,
       "Bohlin's map: along z = c_- e^(-i phi) + c_+ e^(i phi), q = z^2 satisfies q - 4 c_- c_+ = (c_- e^(-i phi) - c_+ e^(i phi))^2, hence "
       "|q| + |q - 4 c_- c_+| = |z|^2 + |c_- e^(-i phi) - c_+ e^(i phi)|^2 = 2(|c_-|^2 + |c_+|^2): an ellipse with one focus at the origin, "
       "semi-major axis |c_-|^2 + |c_+|^2 and eccentricity 2|c_-||c_+|/(|c_-|^2 + |c_+|^2) (exact)")

# 4. Levi-Civita: with dt/dtau = lambda |z|^2, Q(t) = z(tau)^2 solves d^2Q/dt^2 = -mu Q/|Q|^3 with mu = (4 kappa^2/lambda^2)(|c_-|^2+|c_+|^2),
#    Kepler energy -2 kappa^2/lambda^2 and semi-major axis |c_-|^2 + |c_+|^2
zt = sp.expand_complex(z)
zr, zi = sp.re(zt), sp.im(zt)
n2 = sp.expand(zr**2 + zi**2)
Qr, Qi = sp.expand(zr**2 - zi**2), sp.expand(2 * zr * zi)
dt = lamb * n2
dQr, dQi = sp.diff(Qr, tau) / dt, sp.diff(Qi, tau) / dt
ddQr, ddQi = sp.diff(dQr, tau) / dt, sp.diff(dQi, tau) / dt
S2 = a1**2 + b1**2 + a2**2 + b2**2
muK = 4 * kap**2 * S2 / lamb**2
# |Q|^3 = n2^3
res_r = sp.simplify(sp.expand_trig(sp.simplify(ddQr * n2**3 + muK * Qr)))
res_i = sp.simplify(sp.expand_trig(sp.simplify(ddQi * n2**3 + muK * Qi)))
EK = sp.simplify((dQr**2 + dQi**2) / 2 - muK / n2)
report(4, res_r == 0 and res_i == 0 and is_zero(sp.expand_trig(EK + 2 * kap**2 / lamb**2)),
       "Levi-Civita: with dt/dtau = lambda|z|^2 the image Q(t) = z(tau)^2 solves d^2Q/dt^2 = -mu Q/|Q|^3 with mu = 4 kappa^2 (|c_-|^2 + |c_+|^2)"
       "/lambda^2, and its Kepler energy is -2 kappa^2/lambda^2 for every interior orbit (exact, general complex c_+-)")

# 5. no common inverse-square law: at fixed lambda, mu = 2|E_K| a_K grows with the size of the interior orbit
mu_of = lambda cmv, cpv, lam=1, k=1: 4 * k**2 * (abs(cmv)**2 + abs(cpv)**2) / lam**2
mus = [sp.nsimplify(mu_of(1, sp.Rational(2, 5))), sp.nsimplify(mu_of(1, sp.Rational(4, 5)))]
report(5, mus[0] != mus[1] and mus == [sp.Rational(116, 25), sp.Rational(164, 25)],
       f"at kappa = lambda = 1 the interior orbits (c_-,c_+) = (1, 2/5) and (1, 4/5) are carried to Kepler orbits of the inverse-square "
       f"couplings mu = {mus[0]} and {mus[1]}: the image of the family is a family of Kepler problems, one coupling per orbit, all of energy "
       f"-2 kappa^2/lambda^2 (exact)")

# 6. beyond one plane: for z in C^3 and q(z) = z.z = A e^(-2i phi) + C + B e^(2i phi) (A = c_-.c_-, B = c_+.c_+, C = 2 c_-.c_+), the foci are
#    C +- 2 sqrt(AB); the origin is a focus iff (c_-.c_+)^2 = (c_-.c_-)(c_+.c_+). Counterexample c_- = e1, c_+ = e2
sA, tB = sp.symbols("sA tB")                   # sqrt(A) e^(-i phi), sqrt(B) e^(i phi); their product sqrt(AB) is constant
ident6 = is_zero(sp.expand((sA**2 + tB**2 - 2 * sA * tB) - (sA - tB)**2)) and is_zero(sp.expand((sA**2 + tB**2 + 2 * sA * tB) - (sA + tB)**2))
e1, e2 = sp.Matrix([1, 0, 0]), sp.Matrix([0, 1, 0])
zv = e1 * sp.exp(-sp.I * ph) + e2 * sp.exp(sp.I * ph)
qv = sp.simplify(sp.expand_complex((zv.T * zv)[0]))
hv = sp.simplify(sp.expand_complex((zv.H * zv)[0]))
gram = sp.Matrix([[(e1.T * e1)[0], (e1.T * e2)[0]], [(e1.T * e2)[0], (e2.T * e2)[0]]]).det()
report(6, ident6 and is_zero(qv - 2 * sp.cos(2 * ph)) and hv == 2 and gram != 0,
       f"for z in C^3 the image q(z) = z.z of an interior orbit is the ellipse with foci C +- 2 sqrt(AB) (|q - F1| + |q - F2| = "
       f"2(|A| + |B|)), so the origin is a focus iff the Gram determinant (c_-.c_-)(c_+.c_+) - (c_-.c_+)^2 vanishes, which holds in one plane; "
       f"for c_- = e1, c_+ = e2 the Gram determinant is {gram}, q = {qv} passes through the origin while h(z,z) = {hv}: no Levi-Civita "
       "image of a Kepler orbit (exact)")

# 7. spin map: the duality pendulum I0 thetaddot = (R0/m) sin(2 theta + phi0) is the circular-orbit spin-orbit equation
#    C gammaddot = -(3/2)(GM/a^3)(B - A) sin 2 gamma under theta = gamma + (pi - phi0)/2, I0 = C, R0/m = (3/2)(GM/a^3)(B - A)
t = sp.symbols("t", real=True)
gam = sp.Function("gamma")(t)
phi0, I0, Cm, BA, GMa3 = sp.symbols("phi0 I0 C BmA GMa3", positive=True)
R0m = sp.Rational(3, 2) * GMa3 * BA
theta_of = gam + (sp.pi - phi0) / 2
es_eq = I0 * sp.diff(theta_of, t, 2) - R0m * sp.sin(2 * theta_of + phi0)
so_eq = Cm * sp.diff(gam, t, 2) + sp.Rational(3, 2) * GMa3 * BA * sp.sin(2 * gam)
same = is_zero(sp.expand_trig(es_eq.subs(I0, Cm) - so_eq))
Om2 = 2 * R0m / (Cm)
w02 = 3 * GMa3 * BA / Cm
report(7, same and is_zero(Om2 - w02),
       "the duality pendulum I0 thetaddot = (R0/m) sin(2 theta + phi0) becomes the circular-orbit spin-orbit equation "
       "C gammaddot = -(3/2)(GM/a^3)(B - A) sin 2 gamma under theta = gamma + (pi - phi0)/2, I0 = C, R0/m = (3/2)(GM/a^3)(B - A); the stable "
       "point psi = pi is gamma = 0 (long axis towards the primary: synchronous rotation) and the frequencies 2R0/(m I0) and 3(GM/a^3)(B-A)/C "
       "agree (exact)")

# 8. what survives as spin-orbit separation: with B = A the torque vanishes (the spin decouples), as hdot = 0 for duality-class inertia (rb2-6)
report(8, is_zero(so_eq.subs(BA, 0) - Cm * sp.diff(gam, t, 2)) and is_zero(R0m.subs(BA, 0)),
       "with B = A (a body axisymmetric about its spin axis) the spin-orbit torque vanishes and the spin decouples from the orbit, the analogue "
       "of R0 = 0 or duality-class inertia in the closure (exact)")


# ------------------------------------------------------------------------------------------------ numerical part
def kepler_E(M_, e):
    E_ = M_ if e < 0.8 else np.pi
    for _ in range(60):
        E_ = E_ - (E_ - e * np.sin(E_) - M_) / (1 - e * np.cos(E_))
    return E_


def a_over_r_cube_cos(M_, p_, e):
    E_ = kepler_E(M_, e)
    r = 1 - e * np.cos(E_)
    f = 2 * np.arctan2(np.sqrt(1 + e) * np.sin(E_ / 2), np.sqrt(1 - e) * np.cos(E_ / 2))
    return r**-3 * np.cos(2 * f - 2 * p_ * M_)


# 9. resonance strengths of the eccentric spin-orbit problem, H(p,e) = (1/2pi) int (a/r)^3 cos(2f - 2pM) dM, against the series of
#    Goldreich-Peale / Murray-Dermott at small e: H(1/2) = -e/2 + e^3/16, H(1) = 1 - 5e^2/2 + 13e^4/16, H(3/2) = 7e/2 - 123e^3/16
def Hpe(p_, e):
    val, err = quad(lambda M_: a_over_r_cube_cos(M_, p_, e), 0, 2 * np.pi, limit=400, epsabs=1e-13, epsrel=1e-13)
    return val / (2 * np.pi)


ser = {0.5: lambda e: -e / 2 + e**3 / 16, 1.0: lambda e: 1 - 5 * e**2 / 2 + 13 * e**4 / 16, 1.5: lambda e: 7 * e / 2 - 123 * e**3 / 16}
e_small = 0.05
dev = max(abs(Hpe(p_, e_small) - ser[p_](e_small)) for p_ in ser)
H02 = {p_: Hpe(p_, 0.2) for p_ in (0.5, 1.0, 1.5, 2.0)}
H00 = {p_: Hpe(p_, 0.0) for p_ in (0.5, 1.0, 1.5)}
ok9 = dev < 1e-5 and abs(H00[1.0] - 1) < 1e-12 and abs(H00[0.5]) < 1e-12 and abs(H00[1.5]) < 1e-12 and H02[1.5] > 0.6
report(9, ok9, f"quadrature (scipy quad, tolerance 1e-13): H(p,e) agrees with the series to {dev:.1e} at e = 0.05 (truncation O(e^5)); at e = 0 only "
       f"p = 1 survives; at e = 0.2, H(1/2) = {H02[0.5]:.4f}, H(1) = {H02[1.0]:.4f}, H(3/2) = {H02[1.5]:.4f}, H(2) = {H02[2.0]:.4f}: the "
       "eccentric problem has resonances at half-integer spin rates, the autonomous duality pendulum only the one equilibrium family")

TOL = dict(method="DOP853", rtol=1e-12, atol=1e-12)


# 10. the Kepler benchmark: e = 0.5, mu = a = 1, 100 orbits; drift of energy, angular momentum and the Laplace-Runge-Lenz vector
def kep_rhs(t_, u, mu_=1.0):
    r3 = (u[0]**2 + u[1]**2)**1.5
    return [u[2], u[3], -mu_ * u[0] / r3, -mu_ * u[1] / r3]


e_b = 0.5
u0 = [1 - e_b, 0.0, 0.0, np.sqrt((1 + e_b) / (1 - e_b))]           # pericentre, mu = a = 1
Tk = 2 * np.pi
solk = solve_ivp(kep_rhs, (0, 100 * Tk), u0, dense_output=True, **TOL)


def kep_inv(u, mu_=1.0):
    r = np.hypot(u[0], u[1])
    E_ = 0.5 * (u[2]**2 + u[3]**2) - mu_ / r
    Lz = u[0] * u[3] - u[1] * u[2]
    A_ = np.array([u[3] * Lz - mu_ * u[0] / r, -u[2] * Lz - mu_ * u[1] / r])
    return E_, Lz, A_


E0k, L0k, A0k = kep_inv(u0)
tt = np.linspace(0, 100 * Tk, 20001)
dE = dL = dA = 0.0
for tk in tt:
    Ek, Lk, Ak = kep_inv(solk.sol(tk))
    dE, dL, dA = max(dE, abs(Ek / E0k - 1)), max(dL, abs(Lk / L0k - 1)), max(dA, float(np.linalg.norm(Ak - A0k)))
back = solk.sol(100 * Tk)
ret = float(np.hypot(back[0] - u0[0], back[1] - u0[1]))
report(10, solk.success and dE < 1e-9 and dL < 1e-9 and dA < 1e-8 and ret < 1e-6,
       f"Kepler benchmark (DOP853, rtol = atol = 1e-12; mu = a = 1, e = 0.5, start at pericentre, 100 periods 2 pi a^(3/2)/sqrt(mu)): max relative "
       f"drift of the energy {dE:.1e} and of the angular momentum {dL:.1e}, max drift of the Laplace-Runge-Lenz vector {dA:.1e} (no apsidal "
       f"precession beyond integration error), return to the pericentre after 100 periods within {ret:.1e} (accumulated phase error; "
       f"limits 1e-9, 1e-9, 1e-8, 1e-6)")

# 11. the interior orbit and its Levi-Civita image, integrated independently: (c_-, c_+) = (1, 2/5), kappa = lambda = 1, mu = 116/25
kap_n, lam_n, cmv, cpv = 1.0, 1.0, 1.0, 0.4


def hooke_rhs(tau_, u):                       # u = (x, y, xdot, ydot, t)
    return [u[2], u[3], -kap_n**2 * u[0], -kap_n**2 * u[1], lam_n * (u[0]**2 + u[1]**2)]


zx0, zy0 = cmv + cpv, 0.0                      # z(0) = c_- + c_+
vx0, vy0 = 0.0, kap_n * (cpv - cmv)            # zdot(0) = i kappa (c_+ - c_-)
solh = solve_ivp(hooke_rhs, (0, 10 * np.pi / kap_n), [zx0, zy0, vx0, vy0, 0.0], dense_output=True, **TOL)
mu_n = 4 * kap_n**2 * (cmv**2 + cpv**2) / lam_n**2
Q0 = complex(zx0, zy0)**2
dQ0 = 2 * complex(zx0, zy0) * complex(vx0, vy0) / (lam_n * (zx0**2 + zy0**2))
T_end = solh.y[4, -1]
solq = solve_ivp(lambda t_, u: kep_rhs(t_, u, mu_n), (0, T_end), [Q0.real, Q0.imag, dQ0.real, dQ0.imag], dense_output=True, **TOL)
taus = np.linspace(0, 10 * np.pi / kap_n, 4001)
maxdev = 0.0
for tv in taus:
    u = solh.sol(tv)
    Qes = complex(u[0], u[1])**2
    qk = solq.sol(u[4])
    maxdev = max(maxdev, abs(Qes - complex(qk[0], qk[1])))
aK = cmv**2 + cpv**2
T_kep = 2 * np.pi * np.sqrt(aK**3 / mu_n)
report(11, solh.success and solq.success and maxdev < 1e-8 and abs(T_end - 10 * T_kep) < 1e-8,
       f"numerical (DOP853, rtol = atol = 1e-12): the interior orbit (c_-,c_+) = (1, 2/5), kappa = lambda = 1, integrated in tau together with "
       f"dt/dtau = lambda|z|^2, and the Kepler orbit with mu = {mu_n:.2f} from the image initial data, integrated in t, agree to {maxdev:.1e} "
       f"over ten Kepler periods; tau = 10 pi/kappa corresponds to t = {T_end:.10f} = 10 x 2 pi sqrt(a^3/mu)")


# 12. the circular-orbit spin-orbit problem and the duality pendulum of the closure, integrated independently
eps_ = 1e-2                                    # (B - A)/C
k_so = 1.5 * eps_                              # (3/2)(GM/a^3)(B - A)/C with GM = a = 1


def so_rhs(t_, u, asph):                       # u = (X, Y, VX, VY, theta, thetadot): orbit (GM = 1) and spin; asph = (B - A)/C
    r = np.hypot(u[0], u[1])
    f = np.arctan2(u[1], u[0])
    return [u[2], u[3], -u[0] / r**3, -u[1] / r**3, u[5], -1.5 * asph / r**3 * np.sin(2 * (u[4] - f))]


def so_run(e, theta0, thetadot0, T, asph):     # start at pericentre (a = 1), f(0) = 0
    u0_ = [1 - e, 0.0, 0.0, np.sqrt((1 + e) / (1 - e)), theta0, thetadot0]
    return solve_ivp(lambda t_, u: so_rhs(t_, u, asph), (0, T), u0_, dense_output=True, method="DOP853", rtol=1e-11, atol=1e-11)


# closure data: p0 = e1 (w = 1), nu real, m = I0 = 1: R0 = nu, phi0 = 0, so theta_es = gamma + pi/2
Gm = np.zeros((6, 6))
for i in range(3):
    Gm[i + 3, i], Gm[i, i + 3] = 1.0, -1.0
etam = np.diag([1.0, 1.0, 1.0, -1.0, -1.0, -1.0])


def closure_rhs(t_, u, NL):                    # u = (theta, h, p1..p6), I0 = 1
    th_, h_, p_ = u[0], u[1], u[2:8]
    return np.concatenate([[h_, -p_ @ Gm @ NL @ p_], h_ * Gm.T @ p_])


def closure_run(theta0, h0, T, k_es, t_eval=None):
    NL = k_es * np.linalg.inv(etam)            # nu = k_es, m = 1: R0/m = k_es
    p0 = np.array([1.0, 0, 0, 0, 0, 0])
    Eth = np.cos(theta0) * np.eye(6) - np.sin(theta0) * Gm     # e^{-theta0 G} p0: the leaf point at angle theta0
    u0_ = np.concatenate([[theta0, h0], Eth @ p0])
    return solve_ivp(lambda t_, u: closure_rhs(t_, u, NL), (0, T), u0_, dense_output=True, method="DOP853", rtol=1e-11, atol=1e-11,
                     t_eval=t_eval)


gamma0 = 0.5
T12 = 50 * 2 * np.pi
so0 = so_run(0.0, gamma0, 1.0, T12, eps_)      # spin angle theta = gamma + f, f = t on the circular orbit, thetadot = n + gammadot
cl0 = closure_run(gamma0 + np.pi / 2, 0.0, T12, k_so)
dmax = 0.0
for tv in np.linspace(0, T12, 5001):
    us, uc = so0.sol(tv), cl0.sol(tv)
    f_ = np.arctan2(us[1], us[0])
    g_so = np.angle(np.exp(1j * (us[4] - f_)))
    g_es = np.angle(np.exp(1j * (uc[0] - np.pi / 2)))
    dmax = max(dmax, abs(np.angle(np.exp(1j * (g_so - g_es)))))
Hcl = lambda u: 0.5 * u[1]**2 + 0.5 * k_so * (u[2:8] @ np.linalg.inv(etam) @ u[2:8])
dH = max(abs(Hcl(cl0.sol(tv)) - Hcl(cl0.sol(0))) for tv in np.linspace(0, T12, 501))
report(12, so0.success and cl0.success and dmax < 1e-7 and dH < 1e-9,
       f"numerical (DOP853, rtol = atol = 1e-11): on a circular orbit (e = 0, GM = a = 1, (B-A)/C = 0.01) the spin-orbit angle gamma = theta - f "
       f"from gamma(0) = {gamma0}, gammadot(0) = 0, and the rotor of the closure (p0 = e1, Lorentz-class N = 0.015 eta^-1, m = I0 = 1, "
       f"theta(0) = gamma(0) + pi/2, h(0) = 0), integrated independently, agree to {dmax:.1e} rad over 50 orbits (libration: synchronous "
       f"rotation); drift of the closure Hamiltonian {dH:.1e}")

# 13. an eccentric orbit breaks the equivalence: e = 0.2, (B-A)/C = 1e-3, start at pericentre with the long axis towards the primary and spin
#     rate 1.52 n: the spin-orbit problem is captured in the 3:2 resonance, the duality pendulum (the circular-orbit problem) is not
eps13 = 1e-3
T13 = 1000 * 2 * np.pi
so13 = so_run(0.2, 0.0, 1.52, T13, eps13)
cl13 = closure_run(np.pi / 2, 0.52, T13, 1.5 * eps13)        # gamma(0) = theta(0) - f(0) = 0, gammadot(0) = 1.52 - 1
t1, t2 = 200 * 2 * np.pi, T13
mean_so = (so13.sol(t2)[4] - so13.sol(t1)[4]) / (t2 - t1)
mean_es = 1.0 + (cl13.sol(t2)[0] - cl13.sol(t1)[0]) / (t2 - t1)     # spin rate = n + gammadot
sig = [so13.sol(tv)[4] - 1.5 * tv for tv in np.linspace(0, T13, 4001)]   # resonant angle theta - (3/2) M, M = n t
sig_amp = float(np.max(np.abs(np.angle(np.exp(2j * np.array(sig))) / 2)))
orbE = [0.5 * (u[2]**2 + u[3]**2) - 1 / np.hypot(u[0], u[1]) for u in (so13.sol(tv) for tv in np.linspace(0, T13, 2001))]
dEo = max(abs(np.array(orbE) / orbE[0] - 1))
# the mean rate of the rotating pendulum gammaddot = -k sin 2 gamma, k = (3/2)(B-A)/C, in closed form: pi / int_0^pi dgamma/gammadot
k13 = 1.5 * eps13
E13 = 0.5 * 0.52**2 - 0.5 * k13
T_pi, _ = quad(lambda g: 1.0 / np.sqrt(2 * E13 + k13 * np.cos(2 * g)), 0, np.pi, epsabs=1e-13, epsrel=1e-13)
mean_es_exact = 1.0 + np.pi / T_pi
ok13 = (so13.success and cl13.success and abs(mean_so - 1.5) < 2e-3 and sig_amp < 1.0 and abs(mean_es - mean_es_exact) < 1e-5
        and abs(mean_es - 1.5) > 10 * abs(mean_so - 1.5) and dEo < 1e-6)
report(13, ok13,
       f"numerical (DOP853, rtol = atol = 1e-11, 1000 orbits; GM = a = 1, e = 0.2, (B-A)/C = 1e-3, theta(0) = f(0) = 0, thetadot(0) = 1.52 n): "
       f"the eccentric spin-orbit problem is captured in the 3:2 resonance (mean spin rate over orbits 200-1000: {mean_so:.5f} n; the resonant "
       f"angle theta - 3M/2 librates with amplitude {sig_amp:.2f} rad), while the duality pendulum with the same initial data keeps the mean "
       f"rate {mean_es:.5f} n, which is the closed-form mean rate {mean_es_exact:.5f} n of the circular-orbit pendulum; relative drift of the "
       f"orbital energy {dEo:.1e} (limit 1e-6)")

if WRITE:
    out = ROOT / "manuscript" / "figures" / "ch04"
    out.mkdir(parents=True, exist_ok=True)
    # orbit-averaged spin rates (theta(t + 2 pi) - theta(t))/(2 pi): the average over one orbit removes the short-period terms
    grid = np.linspace(0, 200 * 2 * np.pi, 401)                                  # the first 200 orbits, two samples per orbit
    with open(out / "spinrate_so.dat", "w") as fh:
        fh.write("% generated by: python3 verification/rigidbody/rb3_orbital.py --write-data (check rb3-13)\n"
                 "% t/(2pi)  orbit-averaged spin rate / n of the spin-orbit problem, e = 0.2\n")
        for tv in grid:
            fh.write(f"{tv / (2 * np.pi):.4f} {(so13.sol(tv + 2 * np.pi)[4] - so13.sol(tv)[4]) / (2 * np.pi):.6f}\n")
    with open(out / "spinrate_es.dat", "w") as fh:
        fh.write("% generated by: python3 verification/rigidbody/rb3_orbital.py --write-data (check rb3-13)\n"
                 "% t/(2pi)  orbit-averaged spin rate / n of the duality pendulum, 1 + (gamma(t + 2 pi) - gamma(t))/(2 pi)\n")
        for tv in grid:
            fh.write(f"{tv / (2 * np.pi):.4f} {1.0 + (cl13.sol(tv + 2 * np.pi)[0] - cl13.sol(tv)[0]) / (2 * np.pi):.6f}\n")
    print(f"data written to {out}")

print(f"\n{NPASS[0]} PASS, {len(FAILS)} FAIL")
sys.exit(1 if FAILS else 0)
