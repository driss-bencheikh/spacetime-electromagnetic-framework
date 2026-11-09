#!/usr/bin/env python3
"""c2: dispersion, tachyon counts, the surplus constraints and the worldline brackets of Chapters 11-12 (Phase 9C).

    python3 verification/ch12_causality/c2_dynamics_brackets.py

prop:toy, prop:orient, ex:tachyons, the mode content of p_u = p_v = 0 (CORE-07), prop:constraints, thm:wlgauge (brackets, closure, Jacobi,
gauge invariance of the extended action), prop:wltau, prop:isometryonly. Exact symbolic arithmetic (sympy) with general functions h^{mu nu}(x)
and A_mu(x) where the text claims generality; exact integer counting for the lattice statements.
"""
import itertools
import sys

import sympy as sp

NPASS, FAILS = [0], []


def report(cid, ok, text):
    if ok:
        NPASS[0] += 1
        print(f"PASS c2-{cid}: {text}")
    else:
        FAILS.append(cid)
        print(f"FAIL c2-{cid}: {text}")
    sys.stdout.flush()


# ------------------------------------------------------------------ prop:toy
t, u, v, x, y, z = sp.symbols("t u v x y z", real=True)
w, q1, q2, k1, k2, k3, m, a, q = sp.symbols("omega q1 q2 k1 k2 k3 m a q", real=True)
phi = sp.Function("phi")(t, u, v, x, y, z)
etaA = [1, 1, 1, -1, -1, -1]; XA = [t, u, v, x, y, z]
L = sp.Rational(1, 2) * (sum(etaA[i] * sp.diff(phi, XA[i]) ** 2 for i in range(6)) - m ** 2 * phi ** 2)
EL = sp.euler_equations(L, [phi], XA)[0].lhs
wave = sp.exp(sp.I * (-w * t + q1 * u + q2 * v + k1 * x + k2 * y + k3 * z))
disp = sp.simplify(EL.subs(phi, wave).doit() / wave)
kk2 = k1 ** 2 + k2 ** 2 + k3 ** 2
ok = sp.simplify(sp.solve(disp, w ** 2)[0] - (kk2 + m ** 2 - q1 ** 2 - q2 ** 2)) == 0 if disp.has(w) else False
pi_ = sp.diff(phi, t)
Hden = pi_ * sp.diff(phi, t) - L
Hclaim = sp.Rational(1, 2) * (pi_ ** 2 - sp.diff(phi, u) ** 2 - sp.diff(phi, v) ** 2 + sum(sp.diff(phi, c_) ** 2 for c_ in (x, y, z)) + m ** 2 * phi ** 2)
ok &= sp.simplify(Hden - Hclaim) == 0
qp = sp.symbols("q_p", positive=True)
cfg = a * sp.sin(qp * u)
avg = sp.simplify(sp.integrate(Hclaim.subs(phi, cfg).doit(), (u, 0, 2 * sp.pi / qp)) * qp / (2 * sp.pi))
ok &= sp.simplify(avg - a ** 2 * (m ** 2 - qp ** 2) / 4) == 0
s_ = sp.symbols("s", positive=True)
gam = sp.sqrt(q ** 2 - kk2 - m ** 2)
ok &= sp.limit((1 + q ** 2) ** (-s_ / 2) * sp.exp(sp.sqrt(q ** 2 - 1)), q, sp.oo) == sp.oo
report("1", ok, "prop:toy: plane waves of the free scalar on R^{3,3} (t,u,v temporal) obey omega^2 = |k|^2 + m^2 - q1^2 - q2^2; the canonical energy density is "
       "(1/2)(pi^2 - phi_u^2 - phi_v^2 + |grad phi|^2 + m^2 phi^2) and averages to a^2(m^2 - q^2)/4 on phi = a sin(qu) (unbounded below as q grows); for "
       "q^2 > |k|^2 + m^2 the growth e^gamma beats every Sobolev weight (1+q^2)^(-s/2) (exact). Hadamard ill-posedness is the inherited consequence")

# ------------------------------------------------------------------ prop:orient
sv = sp.symbols("s", real=True)
Xs = [sp.cos(sv), sp.sin(sv)] + [sp.S(0)] * 4
ok = sp.simplify(sum(etaA[i] * Xs[i] ** 2 for i in range(6)) - 1) == 0
ok &= [c_.subs(sv, sp.pi) for c_ in Xs] == [-c_.subs(sv, 0) for c_ in Xs]
report("2", ok, "prop:orient: the path X(s) = (cos s, sin s, 0; 0,0,0) has eta(X,X) = 1 throughout and joins X(0) to -X(0), so the timelike set of R^{3,3} is "
       "connected and has no past/future splitting (exact); for one time the two components are the textbook statement")

# ------------------------------------------------------------------ ex:tachyons and the mode content of the constraints (CORE-07, 9C.8)
def count(N, mR=3):
    pts = [(n1, n2) for n1 in range(-N, N + 1) for n2 in range(-N, N + 1)]
    tach = sum(1 for n1, n2 in pts if n1 * n1 + n2 * n2 > mR * mR)
    stable_nonzero = sum(1 for n1, n2 in pts if 0 < n1 * n1 + n2 * n2 <= mR * mR)
    return tach, stable_nonzero, len(pts)


c5, c10, c20 = count(5), count(10), count(20)
ok = (c5[0], c10[0], c20[0]) == (92, 412, 1652)
mm, R, n1, n2 = sp.symbols("m R n1 n2", positive=True)
ok &= c5[1] == 28 and c10[1] == 28
# imaginary frequency depends on the spatial momentum: a tachyonic mass m_n^2 < 0 is unstable only for |k|^2 < -m_n^2
mn2 = 1 - sp.Rational(25, 9)                                              # m = 1, R = 3, n = (5,0): m_n^2 = -16/9
ok &= mn2 < 0
unstable_k1 = (1 + mn2) < 0; stable_k2 = (2 + mn2) > 0                     # omega^2 = |k|^2 + m_n^2
ok &= unstable_k1 and stable_k2
report("3", ok, f"ex:tachyons and 9C.8: with compact surplus times (m R = 3) the masses m_n^2 = m^2 - |n|^2/R^2 are negative for {c5[0]}, {c10[0]}, {c20[0]} lattice "
       f"points with |n_i| <= 5, 10, 20, as stated, and infinitely many in all; (1) negative m_n^2 is a mass statement; (2) a mode is unstable (imaginary "
       f"omega) only for |k|^2 < -m_n^2: for n = (5,0), m = 1, R = 3 it is unstable at |k|^2 = 1 and stable at |k|^2 = 2; (3) the constraints p_u = p_v = 0 "
       f"remove every n != 0, i.e. the {c5[0]} tachyonic AND the {c5[1]} non-tachyonic nonzero modes in |n_i| <= 5; (4) only n = 0 is retained (exact counting)")

k0 = sp.Rational(0); qs = (sp.Rational(1, 2), 0); m1 = 1
om2 = k0 + m1 ** 2 - (qs[0] ** 2 + qs[1] ** 2)
ok = om2 == sp.Rational(3, 4) and om2 > 0
report("4", ok, "CORE-07 (9C.8): in the non-compact case the constraints keep only q = 0; the mode k = 0, m = 1, q = (1/2, 0) has omega^2 = 3/4 > 0, a stable "
       "oscillating mode, and it is removed together with the unstable band q^2 > |k|^2 + m^2: the constraints are sufficient for stability of the free "
       "spectrum, not minimal, and the minimal restriction q^2 <= |k|^2 + m^2 depends on the four-momentum (exact). Stability of an interacting theory is "
       "not addressed")

# ------------------------------------------------------------------ prop:constraints, thm:wlgauge: Poisson brackets
X4 = sp.symbols("x0:4", real=True)
P4 = sp.symbols("p0:4", real=True)
pu, pv = sp.symbols("p_u p_v", real=True)
Q = list(X4) + [u, v]; Pm = list(P4) + [pu, pv]


def PB(F, Gf):
    return sp.expand(sum(sp.diff(F, Q[i]) * sp.diff(Gf, Pm[i]) - sp.diff(F, Pm[i]) * sp.diff(Gf, Q[i]) for i in range(6)))


A = [sp.Function(f"A{i}")(*X4) for i in range(4)]
hu = [[sp.Function(f"h{min(i, j)}{max(i, j)}")(*X4) for j in range(4)] for i in range(4)]   # general inverse metric h^{mu nu}(x)
ell = v * pu - u * pv
kv = [P4[i] - A[i] * ell for i in range(4)]
H = sp.Rational(1, 2) * (sum(hu[i][j] * kv[i] * kv[j] for i in range(4) for j in range(4)) + pu ** 2 + pv ** 2)
C = sum(A[i] * hu[i][j] * kv[j] for i in range(4) for j in range(4))
ok = PB(pu, pv) == 0 and sp.expand(PB(pu, H) + C * pv) == 0 and sp.expand(PB(pv, H) - C * pu) == 0 and PB(ell, H) == 0
Hred = sp.Rational(1, 2) * sum(hu[i][j] * P4[i] * P4[j] for i in range(4) for j in range(4))
ok &= sp.expand(H.subs({pu: 0, pv: 0}) - Hred) == 0
report("5", ok, "prop:constraints and thm:wlgauge: for a general inverse metric h^{mu nu}(x) and general A_mu(x), {p_u,p_v} = 0, {p_u,H} = -C p_v, {p_v,H} = C p_u "
       "with C = A_mu h^{mu nu} k_nu, and {l,H} = 0 for l = v p_u - u p_v; so p_u, p_v, H form a first-class set closing with structure functions; on "
       "p_u = p_v = 0, H = (1/2) h^{mu nu} p_mu p_nu (exact). The reduced phase space has dimension 12 - 2x2 = 8 (two first-class constraints, two gauge "
       "orbits); the theorem's restriction to flat h is not needed for the brackets")

ok = True
F1, F2_, F3 = pu, H, u * P4[0] + A[1] * pv
jac = sp.expand(PB(F1, PB(F2_, F3)) + PB(F2_, PB(F3, F1)) + PB(F3, PB(F1, F2_)))
ok &= jac == 0 and sp.expand(PB(F1, F2_) + PB(F2_, F1)) == 0
G1 = u * X4[1]
ok &= sp.expand(PB(F1, F2_ * G1) - (PB(F1, F2_) * G1 + F2_ * PB(F1, G1))) == 0
Nf = sp.Function("N"); l1, l2 = sp.Function("lambda1"), sp.Function("lambda2")
tau = sp.symbols("tau", real=True)
mm = sp.symbols("m_mass", positive=True)
weak = sp.expand(PB(pu, Nf(tau) * (H - mm ** 2 / 2) + l1(tau) * pu + l2(tau) * pv))
ok &= sp.expand(weak + Nf(tau) * C * pv) == 0
report("6", ok, "9C.9: the canonical bracket is antisymmetric, satisfies the Leibniz rule and the Jacobi identity (checked on the triple p_u, H, u p_0 + A_1 p_v "
       "with general h, A); the total Hamiltonian N (H - m^2/2) + lambda1 p_u + lambda2 p_v preserves the constraints weakly, d p_u/d tau = -N C p_v ~ 0 (exact). "
       "This is constraint preservation for one test worldline in a background; it says nothing about interacting fields")

# gauge invariance of the extended action (thm:wlgauge), flat or general h
e1, e2 = sp.Function("epsilon1")(tau), sp.Function("epsilon2")(tau)
qf = [sp.Function(f"q{i}")(tau) for i in range(6)]; pf = [sp.Function(f"P{i}")(tau) for i in range(6)]
sub = {Q[i]: qf[i] for i in range(6)}; sub.update({Pm[i]: pf[i] for i in range(6)})
Hf = H.subs(sub, simultaneous=True); Cf = C.subs(sub, simultaneous=True)
N_, L1, L2 = Nf(tau), l1(tau), l2(tau)
Lag = sum(pf[i] * sp.diff(qf[i], tau) for i in range(6)) - N_ * (Hf - mm ** 2 / 2) - L1 * pf[4] - L2 * pf[5]
dl1 = sp.diff(e1, tau) + N_ * Cf * e2; dl2 = sp.diff(e2, tau) - N_ * Cf * e1
var = (sum(pf[i] * sp.diff(dq, tau) for i, dq in ((4, e1), (5, e2)))
       - N_ * (sp.diff(Hf, qf[4]) * e1 + sp.diff(Hf, qf[5]) * e2) - dl1 * pf[4] - dl2 * pf[5])
ok = sp.simplify(sp.expand(var)) == 0
report("7", ok, "thm:wlgauge: under delta u = eps1(tau), delta v = eps2(tau), delta lambda1 = eps1' + N C eps2, delta lambda2 = eps2' - N C eps1 (all other variables "
       "fixed) the first variation of the extended Lagrangian p.x' - N (H - m^2/2) - lambda1 p_u - lambda2 p_v vanishes identically, not only up to a boundary term, "
       "for general h(x) and A(x) (exact). The gauge principle is chosen: the multipliers are added by hand")

# ------------------------------------------------------------------ prop:wltau
f = sp.Function("f")
dx = [sp.diff(PB(X4[i], Nf(tau) * H), 1) if False else PB(X4[i], Nf(tau) * H) for i in range(4)]
ok = all(sp.expand(dx[i] - Nf(tau) * sum(hu[i][j] * kv[j] for j in range(4))) == 0 for i in range(4))
hl = sp.Matrix(4, 4, lambda i, j: sp.Symbol(f"g{min(i, j)}{max(i, j)}"))
hupper = hl.inv()
mm2 = sp.symbols("m2", positive=True)
xdot = [Nf(tau) * sum(hupper[i, j] * P4[j] for j in range(4)) for i in range(4)]
norm = sp.simplify(sum(hl[i, j] * xdot[i] * xdot[j] for i in range(4) for j in range(4)) - Nf(tau) ** 2 * sum(hupper[i, j] * P4[i] * P4[j] for i in range(4) for j in range(4)))
ok &= norm == 0
tp = sp.symbols("tau_p")
ok &= sp.simplify((Nf(tau) / sp.diff(f(tau), tau)) * sp.diff(f(tau), tau) - Nf(tau)) == 0
# v1.0.1: the mass shell is derived by varying the lapse, not inserted
dN_new = sp.expand(sp.diff(Lag, N_))                                   # Euler-Lagrange for N (no N' in the Lagrangian)
ok &= sp.expand(dN_new + (Hf - mm ** 2 / 2)) == 0
Lag_old = Lag + N_ * (Hf - mm ** 2 / 2) - N_ * Hf
ok &= sp.expand(sp.diff(Lag_old, N_) + Hf) == 0                          # with -N H alone the variation gives H = 0, the massless shell
Hsym = sp.Symbol("H_val")
ok &= sp.solve(-(Hsym - mm ** 2 / 2), Hsym) == [mm ** 2 / 2] and sp.solve(-Hsym, Hsym) == [0]
# on the reduced shell (l = 0, k = p) with 2H = h^{mu nu} p p = m^2: h_{mu nu} x' x' = N^2 m^2 and d tau_h = |N| m d tau
Nr = sp.symbols("N_r", real=True)
ok &= sp.simplify(sp.sqrt(Nr ** 2 * mm ** 2) - sp.Abs(Nr) * mm) == 0
report("8", ok, "prop:wltau: varying N in p.x' - N (H - m^2/2) - lambda1 p_u - lambda2 p_v gives the mass shell 2H = m^2, while -N H alone would give "
       "H = 0 (the shell is derived, not inserted); Hamilton's equation for x^mu is x'^mu = N h^{mu nu} k_nu (general h, A); on the reduced shell "
       "(l = 0, k = p) h_{mu nu} x'^mu x'^nu = N^2 h^{mu nu} p_mu p_nu = N^2 m^2 for a general metric, so d tau_h = |N| m d tau; N d tau is invariant under "
       "tau -> f(tau), N -> N/f' (exact)")

# ------------------------------------------------------------------ prop:isometryonly
r, ph = sp.symbols("r phi", positive=True)
pr, pph = sp.symbols("p_r p_phi", real=True)
pu_pol = sp.cos(ph) * pr - sp.sin(ph) / r * pph
pv_pol = sp.sin(ph) * pr + sp.cos(ph) / r * pph
uu, vv = r * sp.cos(ph), r * sp.sin(ph)
ok = sp.simplify(vv * pu_pol - uu * pv_pol + pph) == 0 and sp.simplify(pu_pol ** 2 + pv_pol ** 2 - pr ** 2 - pph ** 2 / r ** 2) == 0
J0 = sp.besselj(0, q * sp.sqrt(u ** 2 + v ** 2))
lap = sp.simplify(sp.diff(J0, u, 2) + sp.diff(J0, v, 2) + q ** 2 * J0)
ok &= sp.simplify(lap.rewrite(sp.besselj)) == 0 or sp.N(lap.subs({q: 1.3, u: 0.7, v: -0.4})) == 0
report("9", ok, "prop:isometryonly: in polar normal coordinates l = -p_phi and p_u^2 + p_v^2 = p_r^2 + p_phi^2/r^2, so gauging only the rotation leaves p_r; "
       "J0(q r) satisfies (d_u^2 + d_v^2) J0(qr) = -q^2 J0(qr), so the rotation-invariant modes keep omega^2 = |k|^2 + m^2 - q^2 (exact)")

print(f"{NPASS[0]}/{NPASS[0] + len(FAILS)} checks passed")
sys.exit(0 if not FAILS else 1)
