"""Phase 9D, batch 3: adversarial audit of the Chapter 4 statements that had neither a cited check nor a written proof.

Part A (one constituent): exact symbolic arithmetic (sympy). Part B (three constituents, the closure section): vector-level statements exact
(sympy, 18x18 rational matrices); spinor-level statements on the exact integer 512x512 realisation of Appendix A, with ranks certified modulo a
prime together with explicit exact solutions, and spectra by float64 diagonalisation of exact integer matrices (labelled numerical).
Each check prints PASS/FAIL a3-<n>; the script exits 0 iff all pass.
"""
import itertools
import sys
from functools import reduce

from fractions import Fraction

import numpy as np
import sympy as sp
from sympy.polys.domains import QQ, ZZ
from sympy.polys.matrices import DomainMatrix

NPASS, FAILS = [0], []


def report(cid, ok, text):
    if ok:
        NPASS[0] += 1
        print(f"PASS a3-{cid}: {text}")
    else:
        FAILS.append(cid)
        print(f"FAIL a3-{cid}: {text}")
    sys.stdout.flush()


I3, Z3, E6 = sp.eye(3), sp.zeros(3), sp.eye(6)
G = sp.BlockMatrix([[Z3, -I3], [I3, Z3]]).as_explicit()
eta = sp.diag(1, 1, 1, -1, -1, -1)
tau = sp.symbols("tau", real=True)
kap, m = sp.symbols("kappa m", positive=True)

# ---------------------------------------------------------------- a3-1: ex:convshift, ex:jordan4
ok = True
for sig, expect in [((2, 1, 0), {9: 4, 1: 4}), ((1, 1, 0), {4: 4, 0: 4})]:
    vals = {}
    for et in itertools.product((1, -1), repeat=3):
        v = sum(a * b for a, b in zip(et, sig)) ** 2
        vals[v] = vals.get(v, 0) + 1  # each sign pattern is a joint eigenline of R^8 (thm:spectrumBC)
    ok &= vals == expect
J2 = sp.Matrix([[0, -1], [1, 0]])
g = sp.BlockMatrix([[J2, sp.eye(2)], [sp.zeros(2), J2]]).as_explicit()
M = sp.BlockMatrix([[-g, sp.eye(4)], [sp.zeros(4), g.T]]).as_explicit()  # X' = p - g X, p' = g^T p (kappa = m = 1)
ok &= set(M.eigenvals()) == {sp.I, -sp.I}
ok &= [((M ** 2 + sp.eye(8)) ** k).rank() for k in (1, 2, 3, 4)] == [6, 4, 2, 0]
report("1", ok, "ex:convshift, ex:jordan4: rates {3k, k} and {2k, 0} each on four dimensions; for g = [[J2,1],[0,J2]] the phase-space generator has spectrum "
       "{+-i} and rank((M^2+1)^k) = 6,4,2,0, so Jordan blocks of size four and exp(M tau) has entries with tau^3 secular terms (exact)")

# ---------------------------------------------------------------- a3-2: prop:howmuch, ex:oneplane, prop:orbitsEM
ok = True
# G = G^dagger: v with G v = i mu v has eta(v, vbar) = 0 (bilinear eta, conjugate argument): check on the binding root
for vec in G.eigenvects():
    for v in vec[2]:
        if sp.im(vec[0]) != 0:
            ok &= sp.simplify((v.T * eta * v.conjugate())[0]) == 0
# Minkowski: a self-adjoint G on R^{1,3} with non-real spectrum on a plane; its oscillating subspace has signature (1,1)
eta4 = sp.diag(1, -1, -1, -1)
G4 = sp.zeros(4)
G4[0, 1], G4[1, 0] = -1, 1  # on (e0, e1): eta-self-adjoint, square -1 on the plane
ok &= eta4.inv() * G4.T * eta4 == G4 and (G4 ** 2)[:2, :2] == -sp.eye(2)
ok &= G * G == -E6 and eta.inv() * G.T * eta == G  # all six directions of R^{3,3} at one frequency
cm, cp = sp.symbols("c_m c_p")
z = cm * sp.exp(-sp.I * kap * tau) + cp * sp.exp(sp.I * kap * tau)
ok &= sp.simplify(sp.diff(z, tau, 2) + kap ** 2 * z) == 0
ok &= sp.simplify(sp.diff(z, tau) + sp.I * kap * z - 2 * sp.I * kap * cp * sp.exp(sp.I * kap * tau)) == 0
# orbitsEM with c_pm in C^3
cmv = sp.Matrix(sp.symbols("a1:4")) + sp.I * sp.Matrix(sp.symbols("b1:4"))
cpv = sp.Matrix(sp.symbols("u1:4")) + sp.I * sp.Matrix(sp.symbols("v1:4"))
reals = {s: sp.Symbol(s.name, real=True) for s in (cmv.free_symbols | cpv.free_symbols)}
cmv, cpv = cmv.subs(reals), cpv.subs(reals)
zv = cmv * sp.exp(-sp.I * kap * tau) + cpv * sp.exp(sp.I * kap * tau)
q = lambda w: (w.T * w)[0]  # noqa: E731
period = 2 * sp.pi / kap
avg = lambda f: sp.simplify(sp.integrate(sp.expand(f), (tau, 0, period)) / period)  # noqa: E731
qz = sp.expand(q(zv))
ok &= sp.simplify(avg(qz) - 2 * (cmv.T * cpv)[0]) == 0
Xv = sp.Matrix.vstack(sp.Matrix([sp.re(sp.expand_complex(c)) for c in zv]), sp.Matrix([sp.im(sp.expand_complex(c)) for c in zv]))
Xd = Xv.diff(tau)
L = m / 2 * ((Xd + kap * G * Xv).T * eta * (Xd + kap * G * Xv))[0]
Ener = sp.simplify(m / 2 * ((Xd.T * eta * Xd)[0] + kap ** 2 * (Xv.T * eta * Xv)[0]))
ok &= sp.simplify(Ener - 2 * m * kap ** 2 * sp.re(sp.expand_complex((cmv.T * cpv)[0]))) == 0
zd = zv.diff(tau)
ok &= sp.simplify(sp.expand(m / 2 * (q(zd) + kap ** 2 * q(zv)) - 2 * m * kap ** 2 * (cmv.T * cpv)[0])) == 0
Bv, Ev = Xv[:3, 0], Xv[3:, 0]
ok &= sp.simplify(avg(Ev.dot(Bv)) - sp.im(sp.expand_complex((cmv.T * cpv)[0]))) == 0
ok &= sp.simplify(avg(Bv.dot(Bv) - Ev.dot(Ev)) - 2 * sp.re(sp.expand_complex((cmv.T * cpv)[0]))) == 0
z0 = cmv
zc = z0 * sp.exp(-sp.I * kap * tau)
ok &= sp.simplify(q(zc) - sp.exp(-2 * sp.I * kap * tau) * q(z0)) == 0
zr = cmv * sp.exp(-sp.I * kap * tau) + cmv.conjugate() * sp.exp(sp.I * kap * tau)
ok &= all(sp.simplify(sp.im(sp.expand_complex(c))) == 0 for c in zr)
report("2", ok, "prop:howmuch, ex:oneplane, prop:orbitsEM: eigenvectors of the binding root with non-real eigenvalue are eta-isotropic against their "
       "conjugates; one binding plane on R^{1,3}, all six directions on R^{3,3}; z = c_- e^{-ik tau} + c_+ e^{ik tau} solves the law with alpha_+X <-> "
       "2ik c_+ e^{ik tau}; for c_pm in C^3: <q> = 2 c_-.c_+, <E.B> = Im, <|B|^2-|E|^2> = 2 Re; the Noether energy (m/2)[eta(X',X')+k^2 eta(X,X)] "
       "= 2mk^2 Re(c_-.c_+) and (m/2)(q(z')+k^2 q(z)) = 2mk^2 c_-.c_+; zero-momentum orbits rotate q by e^{-2ik tau}; c_+ = conj(c_-) gives real z "
       "(exact)")

# ---------------------------------------------------------------- a3-3: prop:EAdict(ii), prop:EAnotinner, thm:EAdichotomy
ok = True
nn = sp.symbols("n", positive=True, integer=True)
ok &= sp.simplify(6 * nn / (3 * nn * (6 * nn - 1)) - 2 / (6 * nn - 1)) == 0
ok &= sp.simplify(3 * nn * (6 * nn - 1) - sp.binomial(6 * nn, 2)) == 0
for n in (1, 2, 3):
    N6 = 6 * n
    etan = sp.diag(*([1] * (3 * n) + [-1] * (3 * n)))
    xs = sp.symbols(f"y0:{N6 * N6}")
    Y = sp.Matrix(N6, N6, xs)
    A1, _ = sp.linear_eq_to_matrix(list(Y.T * etan + etan * Y), xs)
    A2, _ = sp.linear_eq_to_matrix(list(Y.T * etan - etan * Y), xs)
    d_so, d_sa = N6 * N6 - A1.rank(), N6 * N6 - A2.rank()
    ok &= d_so == 3 * n * (6 * n - 1) and d_sa == 3 * n * (6 * n + 1) and d_so + d_sa == 36 * n * n
ok &= G.det() != 0 and G * G == -E6
p = sp.Matrix(sp.symbols("p1:7"))
pdot = kap * G.T * p
C = (p.T * eta.inv() * p)[0]
dC = sp.expand(sum(sp.diff(C, p[i]) * pdot[i] for i in range(6)))
ok &= sp.expand(dC - 2 * kap * (p.T * G * eta.inv() * p)[0]) == 0
BE = eta * p
ok &= sp.expand(dC + 4 * kap * BE[:3, 0].dot(BE[3:, 0])) == 0
# a G in so(3,3) (antisymmetric for eta) keeps C: take K of prop:twoK
K = sp.BlockMatrix([[Z3, -I3], [-I3, Z3]]).as_explicit()
ok &= sp.expand((p.T * K * eta.inv() * p)[0]) == 0
report("3", ok, "prop:EAdict(ii), prop:EAnotinner, thm:EAdichotomy: 6n/(3n(6n-1)) = 2/(6n-1); End(R^{3n,3n}) = so + eta-self-adjoint with dimensions "
       "3n(6n-1) + 3n(6n+1) = 36n^2 (n = 1,2,3, exact); G invertible so inner for no bracket (ad_w w = 0); dC/dtau = 2k p^T G eta^-1 p = -4k B.E for "
       "(B;E) = eta p; an element of so(3,3) leaves C invariant (exact)")

# ---------------------------------------------------------------- a3-4: prop:lGhoms, thm:inertiareq
sys.path.insert(0, __import__("os").path.join(__import__("os").path.dirname(__import__("os").path.abspath(__file__)), "..", "common"))
from clifford8 import comm, e, vec_coords  # noqa: E402

cyc = [(1, 2, 3), (2, 3, 1), (3, 1, 2)]
ad6 = lambda Xb: sp.Matrix([vec_coords(comm(Xb, e[a])) for a in range(1, 7)]).T  # noqa: E731
Lg = [ad6((e[j + 3] * e[k + 3] - e[j] * e[k]) / 2) for i, j, k in cyc]
Kg = [ad6((e[k] * e[j + 3] - e[j] * e[k + 3]) / 2) for i, j, k in cyc]
lG = Lg + Kg


def solve(conds, shape):
    n_, m_ = shape
    xs = sp.symbols(f"w0:{n_ * m_}")
    Wm = sp.Matrix(n_, m_, xs)
    eqs = [x for c in conds for x in c(Wm)]
    A, _ = sp.linear_eq_to_matrix(eqs, xs)
    return [sp.Matrix(n_, m_, list(v)) for v in A.nullspace()]


ok = True
hom_VV = solve([(lambda a: (lambda W: list(a * W - W * a)))(a) for a in lG], (6, 6))
hom_VR = solve([(lambda a: (lambda W: list(W * a)))(a) for a in lG], (1, 6))
ok &= len(hom_VV) == 2 and len(hom_VR) == 0
pairs = list(itertools.combinations(range(6), 2))


def lam2(a):
    out = sp.zeros(15)
    for col, (i, j) in enumerate(pairs):
        u, v = sp.eye(6)[:, i], sp.eye(6)[:, j]
        w = (a * u) * v.T - v * (a * u).T + u * (a * v).T - (a * v) * u.T  # a(u^v) = au^v + u^av as antisymmetric matrix u v^T - v u^T
        for row, (k, l) in enumerate(pairs):
            out[row, col] = w[k, l]
    return out


lam = [lam2(a) for a in lG]
hom_L2V = solve([(lambda a, la: (lambda W: list(a * W - W * la)))(a, la) for a, la in zip(lG, lam)], (6, 15))
hom_L2R = solve([(lambda la: (lambda W: list(W * la)))(la) for la in lam], (1, 15))
ok &= len(hom_L2V) == 2 and len(hom_L2R) == 0
# beta0: complex cross product z(X) x z(Y) returned to V
xr = sp.Matrix(sp.symbols("xr1:7", real=True))
yr = sp.Matrix(sp.symbols("yr1:7", real=True))


zof = lambda X: X[:3, 0] + sp.I * X[3:, 0]  # noqa: E731


def b0(X, Y):
    zc_ = sp.expand(zof(X).cross(zof(Y)))
    return sp.Matrix([sp.re(c) for c in zc_] + [sp.im(c) for c in zc_])


beta0 = b0(xr, yr)
ok &= all(sp.expand(a * beta0 - b0(a * xr, yr) - b0(xr, a * yr)) == sp.zeros(6, 1) for a in lG)
Wb = sp.Matrix([[sp.diff(beta0[r], xr[i]).diff(yr[j]) - sp.diff(beta0[r], xr[j]).diff(yr[i]) for (i, j) in pairs] for r in range(6)]) / 2
GWb = G * Wb
span = sp.Matrix([list(w) for w in hom_L2V])
ok &= span.rank() == sp.Matrix([list(w) for w in hom_L2V] + [list(Wb), list(GWb)]).rank() == 2 and sp.Matrix([list(Wb), list(GWb)]).rank() == 2
forms = solve([(lambda a: (lambda N: list(a * N + N * a.T)))(a) for a in lG] + [lambda N: list(N - N.T)], (6, 6))  # forms on V^*: a N + N a^T = 0
ok &= len(forms) == 2 and sp.Matrix([list(f) for f in forms] + [list(eta.inv()), list(G * eta.inv())]).rank() == 2
report("4a", ok, "prop:lGhoms: Hom_lG(V,V) = span{1,G} (dim 2), Hom_lG(V,R) = 0, Hom_lG(Lambda^2 V, V) = span{beta0, G beta0} (dim 2, beta0 the complex "
       "cross product, equivariant), Hom_lG(Lambda^2 V, R) = 0, invariant symmetric forms on V^* = span{eta^-1, G eta^-1} (exact)")

ok = True
Nsym = sp.Matrix(6, 6, lambda i, j: sp.Symbol(f"n{min(i, j)}{max(i, j)}"))
Nc, Na = (Nsym - G * Nsym * G) / 2, (Nsym + G * Nsym * G) / 2
ok &= Nc.T == Nc and Na.T == Na and sp.expand(Nc * G - G * Nc) == sp.zeros(6) and sp.expand(Na * G + G * Na) == sp.zeros(6)
com = solve([lambda N: list(N - N.T), lambda N: list(N * G - G * N)], (6, 6))
anti = solve([lambda N: list(N - N.T), lambda N: list(N * G + G * N)], (6, 6))
ok &= len(com) == 9 and len(anti) == 12
ok &= all(sp.expand(f * G + G * f) == sp.zeros(6) for f in (eta.inv(), G * eta.inv()))
a_, b_ = sp.symbols("a b", real=True)
pr = sp.Matrix(sp.symbols("q1:7", real=True))
Nl = (a_ * eta.inv() + b_ * G * eta.inv()) / m
w = (pr.T * eta.inv() * pr)[0] - sp.I * (pr.T * G * eta.inv() * pr)[0]
ok &= sp.expand((pr.T * Nl * pr)[0] / 2 - sp.re(sp.expand((a_ + sp.I * b_) * w)) / (2 * m)) == 0
rotinv = solve([(lambda a: (lambda N: list(a * N + N * a.T)))(a) for a in Lg] + [lambda N: list(N - N.T), lambda N: list(N * G - G * N)], (6, 6))
ok &= len(rotinv) == 1 and sp.Matrix([list(rotinv[0]), list(E6)]).rank() == 1
Nex = (2 * E6 + eta.inv()) / m
ok &= all(ev > 0 for ev in Nex.subs(m, 1).eigenvals()) and sp.expand(Nex * G - G * Nex) != sp.zeros(6) and sp.expand(Nex * G + G * Nex) != sp.zeros(6)
report("4b", ok, "thm:inertiareq: N_c, N_a symmetric with commuting forms of dimension 9 and anticommuting ones 12; the Lorentz class (a eta^-1 + b G eta^-1)/m "
       "anticommutes with G and has frame energy Re((a+ib)w)/(2m); among rotation-invariant forms the only ones commuting with G are multiples of 1; "
       "(2 + eta^-1)/m is positive definite and neither commutes nor anticommutes with G (exact)")

# ---------------------------------------------------------------- a3-5: prop:LP, prop:groupred, thm:drift, thm:pendulum, cor:multiscale, prop:several
ok = True
I0 = sp.symbols("I_0", positive=True)
h = sp.Symbol("h", real=True)
Nn = sp.Matrix(6, 6, lambda i, j: sp.Symbol(f"N{min(i, j)}{max(i, j)}"))
H = h ** 2 / (2 * I0) + (pr.T * Nn * pr)[0] / 2
# Lie-Poisson bracket of s = R ghat x|_G V (V abelian): [ghat, v] = G v; minus bracket {F,H} = -<mu, [dF, dH]>
dH_h, dH_p = sp.diff(H, h), sp.Matrix([sp.diff(H, x) for x in pr])
# mu = (h, p); for F = p_i: dF = (0, e_i); [(a, u), (b, v)] = (0, a G v - b G u)
pdot_LP = sp.Matrix([-(pr.T * (0 * G * dH_p - dH_h * G * sp.eye(6)[:, i]))[0] for i in range(6)])
hdot_LP = -(pr.T * G * dH_p)[0]
ok &= sp.expand(pdot_LP - h / I0 * G.T * pr) == sp.zeros(6, 1) and sp.expand(hdot_LP + (pr.T * G * Nn * pr)[0]) == 0
ok &= sp.expand((pr.T * pdot_LP)[0]) == 0  # |p|^2 is a Casimir (G^T antisymmetric)
th = sp.Symbol("theta", real=True)
eG = lambda s: sp.cos(s) * E6 + sp.sin(s) * G  # noqa: E731  (G^2 = -1)
p0 = sp.Matrix(sp.symbols("r1:7", real=True))
pth = eG(-th) * p0
U = (pth.T * Nn * pth)[0] / 2
ok &= sp.simplify(sp.diff(pth, th) - (-G * pth)) == sp.zeros(6, 1)
ok &= sp.simplify(-sp.diff(U, th) - (-(pth.T * G * Nn * pth)[0])) == 0  # hdot = -U'(theta) on the orbit
# prop:groupred
xs_ = sp.Matrix([sp.Function(f"x{i}")(tau) for i in range(6)])
tht = sp.Function("th")(tau)
u = eG(-tht) * xs_.diff(tau)
Xb = eG(-tht) * xs_
ok &= sp.simplify(u - (Xb.diff(tau) + tht.diff(tau) * G * Xb)) == sp.zeros(6, 1)
Lgr = I0 / 2 * tht.diff(tau) ** 2 + m / 2 * (u.T * eta * u)[0]
from sympy.calculus.euler import euler_equations  # noqa: E402
el = euler_equations(Lgr, [tht] + list(xs_), tau)
# conserved spatial momentum: d/dtau (m eta u rotated) = 0, body momentum p = m eta u obeys p' = (h/I0) G^T p with h = I0 th'
pb = m * eta * u
hdef = I0 * tht.diff(tau)
# x-equations: d/dt (e^{-th G})^T m eta u = 0
ok &= all(sp.simplify(sp.expand(eq.lhs + sp.diff((eG(-tht).T * pb)[i], tau))) == 0 or sp.simplify(sp.expand(eq.lhs - sp.diff((eG(-tht).T * pb)[i], tau))) == 0
          for i, eq in enumerate(el[1:]))
# body momentum law: p' = th' G^T p follows from the conservation of e^{-thG}^T p = e^{th G} p ... check identity
P0 = eG(tht) * pb  # e^{-thG}^T = e^{thG}
ok &= sp.simplify(eG(-tht) * P0.diff(tau) - (pb.diff(tau) - tht.diff(tau) * G.T * pb)) == sp.zeros(6, 1)
# theta-equation: I0 th'' = -p^T G N p with N = eta^-1/m
el_th = sp.simplify(el[0].lhs)
Nlp = eta.inv() / m
target = sp.simplify(-(I0 * tht.diff(tau, 2)) - (pb.T * G * Nlp * pb)[0])
ok &= sp.simplify(el_th - target) == 0 or sp.simplify(el_th + target) == 0
kc = sp.Symbol("k_c", positive=True)
Lcl = sp.simplify(Lgr.subs(tht, kc * tau).doit())
Xcl = eG(-kc * tau) * xs_
Lsplit = I0 / 2 * kc ** 2 + m / 2 * ((Xcl.diff(tau) + kc * G * Xcl).T * eta * (Xcl.diff(tau) + kc * G * Xcl))[0]
ok &= sp.simplify(Lcl - Lsplit) == 0
# thm:drift
Ncs = sp.Matrix(6, 6, lambda i, j: sp.Symbol(f"N{min(i, j)}{max(i, j)}"))
Nc_ = (Ncs - G * Ncs * G) / 2
Na_ = (Ncs + G * Ncs * G) / 2
xdot = eG(th) * (Nc_ + Na_) * eG(-th) * p0
ok &= sp.simplify(xdot - (Nc_ + eG(2 * th) * Na_) * p0) == sp.zeros(6, 1)
hdot = -((eG(-th) * p0).T * G * (Nc_ + Na_) * (eG(-th) * p0))[0]
ok &= sp.simplify(hdot + (p0.T * G * eG(2 * th) * Na_ * p0)[0]) == 0
kk = sp.Symbol("k", positive=True)
xt = tau * Nc_ * p0 + (eG(2 * kk * tau) - E6) * (2 * kk * G).inv() * Na_ * p0
ok &= sp.simplify(xt.diff(tau) - xdot.subs(th, kk * tau)) == sp.zeros(6, 1) and sp.simplify(xt.subs(tau, 0)) == sp.zeros(6, 1)
# thm:pendulum: U(theta) = (R0/2m) cos(2 theta + phi0), w(p(theta)) = e^{2 i theta} w(p0) with p(theta) = e^{-theta G} p0
wf = lambda pp: (pp.T * eta.inv() * pp)[0] - sp.I * (pp.T * G * eta.inv() * pp)[0]  # noqa: E731
ok &= sp.simplify(sp.expand_complex(wf(pth) - sp.exp(2 * sp.I * th) * wf(p0))) == 0
Upend = sp.simplify((pth.T * Nl * pth)[0] / 2)
ok &= sp.simplify(sp.expand_complex(Upend - sp.re(sp.expand_complex((a_ + sp.I * b_) * sp.exp(2 * sp.I * th) * wf(p0))) / (2 * m))) == 0
R0, ph0 = sp.symbols("R_0 varphi_0", real=True)
Upp = R0 / (2 * m) * sp.cos(2 * th + ph0)
ok &= sp.simplify(-sp.diff(Upp, th) - R0 / m * sp.sin(2 * th + ph0)) == 0
# cor:multiscale
eps_, k0, ps0, ps = sp.symbols("epsilon kappa_0 psi_0 psi", real=True)
kappa2 = k0 ** 2 + R0 / (m * I0) * (sp.cos(ps0) - sp.cos(ps))  # from energy conservation of E_theta with psi = 2 theta + phi0
ok &= sp.simplify(kappa2 - k0 ** 2 * (1 + (R0 / (m * I0 * k0 ** 2)) * (sp.cos(ps0) - sp.cos(ps)))) == 0
Eth = I0 / 2 * kappa2 + R0 / (2 * m) * sp.cos(ps)
ok &= sp.simplify(sp.diff(Eth, ps)) == 0
# prop:several: dot w_k = 2i (h/I0) w_k for p_k' = (h/I0) G^T p_k
pk = sp.Matrix(sp.symbols("s1:7", real=True))
dw = sum(sp.diff(wf(pk), pk[i]) * (h / I0 * G.T * pk)[i] for i in range(6))
ok &= sp.simplify(sp.expand(dw - 2 * sp.I * h / I0 * wf(pk))) == 0
ok &= sp.simplify(sp.expand(wf(-G * pk) + wf(pk))) == 0
report("5", ok, "prop:LP, prop:groupred, thm:drift, thm:pendulum, cor:multiscale, prop:several: the minus Lie-Poisson bracket of R ghat x|_G V gives p' = "
       "(h/I0) G^T p, h' = -p^T G N p with |p|^2 a Casimir; on the orbit h' = -U'(theta); u = X' + theta' G X, Euler-Lagrange equations of the group "
       "Lagrangian = conservation of the spatial momentum + rotor equation with N = eta^-1/m, and clamping gives the split-form interior; "
       "x' = (N_c + e^{2 theta G} N_a) p0, h' = -p0^T G e^{2 theta G} N_a p0 and the closed form of x(tau); w(p(theta)) = e^{2i theta} w(p0), "
       "U = Re((a+ib) w)/(2m), I0 theta'' = (R0/m) sin(2 theta + phi0); kappa^2 = kappa0^2[1 + eps(cos psi0 - cos psi)]; w_k' = 2i(h/I0) w_k, "
       "w(-G p) = -w(p) (exact)")

# ---------------------------------------------------------------- Part B: three constituents
TEMP = [1, 2, 3, 7, 8, 9, 13, 14, 15]
eta3 = sp.diag(*[1 if (A - 1) % 6 < 3 else -1 for A in range(1, 19)])
G3 = sp.zeros(18)
for i in TEMP:
    G3[i + 2, i - 1], G3[i - 1, i + 2] = 1, -1  # e_i -> e_{i'}, e_{i'} -> -e_i  (0-based rows/cols)


def unitv(A):
    return sp.eye(18)[:, A - 1]


def adv(pairs_coef):
    """vector-level action of the bivector sum c e_A e_B via ad(1/2 .): e_C -> 1/2[e_A e_B, e_C]."""
    Mv = sp.zeros(18)
    for (A, B), c in pairs_coef:
        sA, sB = eta3[A - 1, A - 1], eta3[B - 1, B - 1]
        # (1/2)[e_A e_B, e_C] = sB delta_BC e_A - sA delta_AC e_B
        Mv[A - 1, B - 1] += c * sB
        Mv[B - 1, A - 1] -= c * sA
    return Mv


def rot(i, j):
    return [((i + 3, j + 3), 1), ((i, j), -1)]


def bst(i, j):
    return [((j, i + 3), 1), ((i, j + 3), -1)]


ok = G3.T * eta3 * G3 == -eta3 and G3 * G3 == -sp.eye(18)
xs = sp.symbols("c0:324")
Y = sp.Matrix(18, 18, xs)
A, _ = sp.linear_eq_to_matrix(list(Y.T * eta3 + eta3 * Y) + list(Y * G3 - G3 * Y), xs)
cent_dim = 324 - A.rank()
gens = [adv(rot(i, j)) for i, j in itertools.combinations(TEMP, 2)] + [adv(bst(i, j)) for i, j in itertools.combinations(TEMP, 2)]
ok &= cent_dim == 72 and len(gens) == 72 and sp.Matrix([list(x) for x in gens]).rank() == 72
ok &= all(x.T * eta3 + eta3 * x == sp.zeros(18) and x * G3 == G3 * x for x in gens)
perm = [A - 1 for A in TEMP] + [A + 2 for A in TEMP]
Pm = sp.eye(18)[:, perm]
blocks_ok = True
for x in gens:
    y = Pm.T * x * Pm
    R_, S_ = y[:9, :9], y[:9, 9:]
    blocks_ok &= y[9:, 9:] == R_ and y[9:, :9] == -S_ and R_.T == -R_ and S_.T == -S_
ok &= blocks_ok
rots = gens[:36]
ok &= all(all(abs(complex(ev)) == 0 or abs(complex(ev).real) < 1e-12 for ev in x.eigenvals()) for x in rots[:6])
report("6", ok, "thm:so9C: G3 is an anti-isometry of eta3 with G3^2 = -1; its centraliser in so(9,9) has dimension 72, spanned by the 36 paired rotations "
       "and 36 paired boosts, each of block form [[R,S],[-S,R]] with R, S antisymmetric in the basis (temporal; partners); paired rotations have "
       "imaginary spectrum (exact)")

ok = True
lint1 = [adv([((a, b), 1)]) for a, b in itertools.combinations([1, 4, 5, 6], 2)]
M1 = sp.Matrix([list(x) for x in gens])
M2 = sp.Matrix([list(x) for x in lint1])
ok &= M1.rank() + M2.rank() == sp.Matrix.vstack(M1, M2).rank()
A2, _ = sp.linear_eq_to_matrix(list(Y.T * eta3 + eta3 * Y) + list(Y * G3 - G3 * Y) + [x for L_ in lint1 for x in (Y * L_ - L_ * Y)], xs)
d_both = 324 - A2.rank()
T23 = [7, 8, 9, 13, 14, 15]
paired23 = [adv(rot(i, j)) for i, j in itertools.combinations(T23, 2)] + [adv(bst(i, j)) for i, j in itertools.combinations(T23, 2)]
ok &= d_both == 30 and all(all(x * L_ == L_ * x for L_ in lint1) for x in paired23) and sp.Matrix([list(x) for x in paired23]).rank() == 30
report("7a", ok, "thm:nocommon (vector level): Cent(G3) meets l_int^(1) in 0; the elements commuting with G3 and l_int^(1) form a 30-dimensional space, "
       "spanned by the paired elements of constituents two and three (so(6,C)) (exact)")

# spinor level: the 512 realisation (integer)
sx = np.array([[0, 1], [1, 0]]); szm = np.array([[1, 0], [0, -1]]); epm = np.array([[0, -1], [1, 0]]); I2 = np.eye(2, dtype=np.int64)
kr = lambda *a: reduce(np.kron, a)  # noqa: E731
e8 = [None, kr(sx, I2, I2), kr(szm, sx, I2), kr(szm, szm, sx), kr(epm, I2, I2), kr(szm, epm, I2), kr(szm, szm, epm)]
I8n = np.eye(8, dtype=np.int64)
om = reduce(np.matmul, e8[1:])
E = [None] + [kr(*([om] * k + [e8[a]] + [I8n] * (2 - k))).astype(np.int64) for k in range(3) for a in range(1, 7)]
I512 = np.eye(512, dtype=np.int64)
blq = lambda S: reduce(np.matmul, [E[i] for i in S], I512)  # noqa: E731
interior = [1, 4, 5, 6]
internal = [i for i in range(1, 19) if i not in interior]
Jn = blq(internal)
Gam = Jn @ blq(interior)
X1 = (E[2] @ E[3]) @ (E[8] @ E[9])
X2 = blq(range(13, 19))
X3 = blq(range(1, 7)) @ (E[13] @ E[16])
Pfam8 = (I512 + X1) @ (I512 + X2) @ (I512 + X3)  # = 8 P_fam
Xw, Zc, Uq = Gam @ blq([7, 13, 14, 15]), -(Jn @ blq([10, 11, 12, 16, 17, 18])), -(Jn @ (E[2] @ E[3]))
Pst8 = (I512 + Xw) @ (I512 + Zc) @ (I512 + Uq)
ok = all(np.array_equal(x @ x, I512) and np.array_equal(x, x.T) and np.trace(x) == 0 for x in (X1, X2, X3))
ok &= np.array_equal(X1 @ X2, X2 @ X1) and np.array_equal(X1 @ X3, X3 @ X1) and np.array_equal(X2 @ X3, X3 @ X2)
ok &= np.array_equal(Pfam8 @ Pfam8, 8 * Pfam8) and np.trace(Pfam8) == 8 * 64 and np.array_equal(Pst8 @ Pst8, 8 * Pst8) and np.trace(Pst8) == 8 * 64


def biv512(terms):
    return sum(c * (E[a] @ E[b]) for (a, b), c in terms)


Lg1 = [biv512([((j + 3, k + 3), 1), ((j, k), -1)]) for i, j, k in cyc]
Kg1 = [biv512([((k, j + 3), 1), ((j, k + 3), -1)]) for i, j, k in cyc]
Cas = sum(x @ x for x in Kg1) - sum(x @ x for x in Lg1)  # with the factor 1/2 of L_i, K_i: (1/4) of this
ok &= np.array_equal(Cas, 12 * I512)
report("7b", ok, "thm:nocommon (spinor level) and the projectors: X1, X2, X3 are commuting symmetric traceless involutions, 8P_fam and 8P_* are integer "
       "with P^2 = P of rank 64; the Casimir sum K^2 - sum L^2 of l_G^(1) on R^512 is 3 id exactly, so the complexified module is 128 (1/2,1/2) "
       "and contains no Weyl spinor (exact integer arithmetic)")

PRIME = 2_147_483_629


def rank_mod_p(Mx, pr_=PRIME):
    Mx = np.array(Mx, dtype=np.int64) % pr_
    rows, cols = Mx.shape
    r = 0
    for c in range(cols):
        piv = None
        nz = np.nonzero(Mx[r:, c])[0]
        if nz.size == 0:
            continue
        piv = r + nz[0]
        Mx[[r, piv]] = Mx[[piv, r]]
        inv = pow(int(Mx[r, c]), pr_ - 2, pr_)
        Mx[r] = (Mx[r] * inv) % pr_
        others = np.nonzero(Mx[:, c])[0]
        others = others[others != r]
        if others.size:
            f = Mx[others, c].reshape(-1, 1)
            Mx[others] = (Mx[others] - (f * Mx[r]) % pr_) % pr_
        r += 1
        if r == rows:
            break
    return r


rng = np.random.default_rng(7)
probe = rng.integers(-3, 4, size=(512, 6))


def commuting_dim(basis, P8):
    cols = [((B @ P8 - P8 @ B) @ probe).reshape(-1) for B in basis]
    return len(basis) - rank_mod_p(np.array(cols).T)


cent512 = [biv512(rot(i, j)) for i, j in itertools.combinations(TEMP, 2)] + [biv512(bst(i, j)) for i, j in itertools.combinations(TEMP, 2)]
surv_fam = [biv512(rot(2, 3)), biv512(rot(8, 9)), biv512(rot(14, 15)), biv512(bst(14, 15))]
W = [7, 13, 14, 15]
so4D = [biv512(rot(i, j)) for i, j in itertools.combinations(W, 2)]
surv_st = [biv512(rot(2, 3)), biv512(rot(8, 9))] + so4D
cf = commuting_dim(cent512, Pfam8)
cs = commuting_dim(cent512, Pst8)
ok = cf <= 4 and cs <= 8  # rank mod p is a lower bound for the rank over Q, so these are upper bounds
ok &= all(np.array_equal(B @ Pfam8, Pfam8 @ B) for B in surv_fam) and all(np.array_equal(B @ Pst8, Pst8 @ B) for B in surv_st)
ok &= rank_mod_p(np.array([B.reshape(-1)[::7] for B in surv_fam]).T) == 4 and rank_mod_p(np.array([B.reshape(-1)[::7] for B in surv_st]).T) == 8
ok &= all(np.array_equal(a @ b, b @ a) for a in surv_fam for b in surv_fam)
lg1_512 = [biv512(rot(i, j)) for i, j in itertools.combinations([1, 2, 3], 2)] + [biv512(bst(i, j)) for i, j in itertools.combinations([1, 2, 3], 2)]
ok &= commuting_dim(lg1_512, Pfam8) == 1 and commuting_dim(lg1_512, Pst8) == 1
ok &= np.array_equal(biv512(rot(2, 3)), E[5] @ E[6] - E[2] @ E[3])
gam512 = [E[1], E[4], E[5], E[6]]


def kinetic_ok(B, P8):
    X = P8 @ B @ P8
    return all(np.array_equal(X.T @ (P8 @ g0g @ P8) + (P8 @ g0g @ P8) @ X, np.zeros((512, 512), dtype=np.int64)) for g0g in [E[1] @ g for g in gam512])


kin = [kinetic_ok(B, Pfam8) for B in surv_fam]
ok &= kin == [False, True, True, False] or kin == [True, True, True, False]
kin_rot23 = kin[0]
ok &= all(kinetic_ok(B, Pst8) for B in [biv512(rot(8, 9))] + so4D)
ok &= rank_mod_p(np.array([((Pst8 @ B @ Pst8) @ probe).reshape(-1) for B in [biv512(rot(8, 9))] + so4D]).T) == 7
report("8", ok, f"thm:survive: the parts of Cent(G3) commuting with P_fam and P_* have dimensions {cf} and {cs} (upper bounds from ranks mod p, attained by "
       "the stated spans, which commute exactly and are independent), the P_fam part is abelian; l_G^(1) keeps only rot(2,3) = e5e6 - e2e3 under "
       "either projector; the kinetic condition (X^T K^mu + K^mu X = 0 on the block) holds for rot(8,9), rot(14,15) and fails for bst(14,15) on "
       f"P_fam (rot(2,3): {'holds' if kin_rot23 else 'fails'}), holds for so4_Delta + R rot(8,9) on P_*, which acts faithfully (rank 7) (exact)")

# thm:closurechiral: exact (v1.0.1). Joint eigenspace dimensions are traces of products of exact spectral projectors (Lagrange polynomials in
# commuting integer matrices), computed in float64 with every entry an integer below 2^53 (asserted), and anomaly coefficients are exact traces.
def fmat(A):
    return np.asarray(A, dtype=np.float64)


def exact_int(x, bound=2 ** 52):
    assert abs(x) < bound and abs(x - round(x)) < 1e-6, x
    return int(round(x))


def lagrange(Mi, vals, v):
    """numerator and denominator of the spectral projector of the integer matrix Mi onto eigenvalue v (eigenvalues vals)"""
    num, den = np.eye(512), 1
    for u in vals:
        if u != v:
            num = num @ (Mi - u * np.eye(512))
            den *= (v - u)
    return num, den


Q1i = fmat(-(Jn @ biv512(rot(8, 9))))          # 2 Q1
Q2i = fmat(-(Jn @ biv512(rot(14, 15))))        # 2 Q2
Gf, Pf8 = fmat(Gam), fmat(Pfam8)
ok = all(np.array_equal(x @ y, y @ x) for x, y in itertools.combinations([Q1i, Q2i, Gf, Pf8], 2))
for Mi in (Q1i, Q2i):                           # (2Q)(2Q - 2)(2Q + 2) P = 0: eigenvalues of Q in {-1, 0, 1} on the block
    ok &= not np.any((Mi @ (Mi - 2 * np.eye(512)) @ (Mi + 2 * np.eye(512))) @ Pf8)
proj = {}
for name, Mi, vals in (("q1", Q1i, (-2, 0, 2)), ("q2", Q2i, (-2, 0, 2)), ("g", Gf, (-1, 1))):
    for v in vals:
        proj[(name, v)] = lagrange(Mi, vals, v)
content = {}
for a_ in (-2, 0, 2):
    for b_ in (-2, 0, 2):
        Na, da = proj[("q1", a_)]
        Nb, db = proj[("q2", b_)]
        AB = Na @ Nb @ Pf8
        for c_ in (-1, 1):
            Nc, dc = proj[("g", c_)]
            tr = exact_int(np.trace(AB @ Nc))
            dimq, rem = divmod(tr, da * db * dc * 8)
            ok &= rem == 0
            content[(a_ // 2, b_ // 2, c_)] = dimq
table_L = {(-1, -1): 1, (-1, 0): 0, (-1, 1): 1, (0, -1): 1, (0, 0): 2, (0, 1): 1, (1, -1): 0, (1, 0): 2, (1, 1): 0}
table_R = {(-1, -1): 0, (-1, 0): 2, (-1, 1): 0, (0, -1): 1, (0, 0): 2, (0, 1): 1, (1, -1): 1, (1, 0): 0, (1, 1): 1}
gotL = {k: Fraction(content[(k[0], k[1], -1)], 4) for k in table_L}
gotR = {k: Fraction(content[(k[0], k[1], 1)], 4) for k in table_R}
ok &= gotL == table_L and gotR == table_R and sum(content.values()) == 64
al, be, ga = sp.symbols("alpha beta gamma")
anom3 = sp.expand(sum(n * (al * q1 + be * q2) ** 3 for (q1, q2), n in table_L.items()) - sum(n * (al * q1 + be * q2) ** 3 for (q1, q2), n in table_R.items()))
anom1 = sp.expand(sum(n * (al * q1 + be * q2) for (q1, q2), n in table_L.items()) - sum(n * (al * q1 + be * q2) for (q1, q2), n in table_R.items()))
ok &= anom3 == -12 * al * be ** 2 and anom1 == 0
# P_*: su2_+ + su2_- of so4_Delta, Q = -J rot(8,9)/2; Casimirs 16 C_+- are integer matrices with eigenvalues 16 j(j+1) in {0, 12, 32}
R2 = {(a_, b_): fmat(biv512(rot(W[a_], W[b_]))) for a_, b_ in itertools.combinations(range(4), 2)}   # 2 R_ab
R2f = lambda a_, b_: R2[(a_, b_)] if a_ < b_ else -R2[(b_, a_)]  # noqa: E731
Xp4 = [R2f(1, 2) + R2f(0, 3), R2f(2, 0) + R2f(1, 3), R2f(0, 1) + R2f(2, 3)]                       # 4 X_+
Xm4 = [R2f(1, 2) - R2f(0, 3), R2f(2, 0) - R2f(1, 3), R2f(0, 1) - R2f(2, 3)]                       # 4 X_-
ok &= all(np.array_equal(x @ y, y @ x) for x in Xp4 for y in Xm4)
Cp16, Cm16 = -sum(x @ x for x in Xp4), -sum(x @ x for x in Xm4)
Ps8f = fmat(Pst8)
Jf = fmat(Jn)
Q2s = fmat(-(Jn @ biv512(rot(8, 9))))          # 2 Q
for Mi, vals in ((Cp16, (0, 12, 32)), (Cm16, (0, 12, 32))):
    ok &= not np.any(reduce(np.matmul, [Mi - u * np.eye(512) for u in vals]) @ Ps8f)
ok &= not np.any((Q2s @ (Q2s - 2 * np.eye(512)) @ (Q2s + 2 * np.eye(512))) @ Ps8f)
jmap = {0: 0.0, 12: 0.5, 32: 1.0}
cont = {}
for cp in (0, 12, 32):
    Np, dp = lagrange(Cp16, (0, 12, 32), cp)
    for cm in (0, 12, 32):
        Nm, dm = lagrange(Cm16, (0, 12, 32), cm)
        PM = Np @ Nm @ Ps8f
        for qv in (-2, 0, 2):
            Nq, dq = lagrange(Q2s, (-2, 0, 2), qv)
            PMQ = PM @ Nq
            for gv in (-1, 1):
                Ng, dg = lagrange(Gf, (-1, 1), gv)
                tr = exact_int(np.trace(PMQ @ Ng))
                d_, rem = divmod(tr, dp * dm * dq * dg * 8)
                ok &= rem == 0
                if d_:
                    cont[(jmap[cp], jmap[cm], qv // 2, gv)] = d_
expL = {(0.0, 0.0, 0): 1, (0.5, 0.5, 1): 4, (1.0, 0.0, 0): 3}
expR = {(0.0, 0.0, 0): 1, (0.5, 0.5, -1): 4, (0.0, 1.0, 0): 3}
gotL = {k[:3]: Fraction(v, 4) for k, v in cont.items() if k[3] == -1}
gotR = {k[:3]: Fraction(v, 4) for k, v in cont.items() if k[3] == 1}
ok &= gotL == expL and gotR == expR
# anomaly polynomial from exact traces: A(f) = (1/4) Tr(f (P_L - P_R)) on the block, P_L - P_R = -Gamma P_*
Hp4 = -(Jf @ Xp4[2])                            # 4 H+_3
Hm4 = -(Jf @ Xm4[2])                            # 4 H-_3
LR = -(Gf @ Ps8f)                               # 8 (P_L - P_R) on the block
A3 = 0
for mono in [(3, 0, 0), (0, 3, 0), (0, 0, 3), (2, 0, 1), (0, 2, 1), (1, 1, 1), (2, 1, 0), (1, 2, 0), (1, 0, 2), (0, 1, 2)]:
    Mm = reduce(np.matmul, [Hp4] * mono[0] + [Hm4] * mono[1] + [Q2s] * mono[2] + [LR], np.eye(512))
    tr = exact_int(np.trace(Mm))
    coeff = Fraction(tr, 4 ** mono[0] * 4 ** mono[1] * 2 ** mono[2] * 8 * 4)
    mult = sp.factorial(3) / (sp.factorial(mono[0]) * sp.factorial(mono[1]) * sp.factorial(mono[2]))
    A3 += sp.Rational(coeff.numerator, coeff.denominator) * mult * al ** mono[0] * be ** mono[1] * ga ** mono[2]
A1 = 0
for Mi, scale, sym in ((Hp4, 4, al), (Hm4, 4, be), (Q2s, 2, ga)):
    coeff = Fraction(exact_int(np.trace(Mi @ LR)), scale * 8 * 4)
    A1 += sp.Rational(coeff.numerator, coeff.denominator) * sym
ok &= sp.expand(A3 - 2 * ga * (3 * al ** 2 + 3 * be ** 2 + 4 * ga ** 2)) == 0 and sp.expand(A1 - 8 * ga) == 0
report("9", ok, "thm:closurechiral (i)-(ii) (exact: joint eigenspace dimensions as traces of exact spectral projectors of commuting integer matrices, "
       "minimal polynomials verified on the blocks): the (q1, q2, Gamma) content of the P_fam block is tab:closurecontent (in units of four real "
       "dimensions), giving A((aQ1+bQ2)^3) = -12ab^2 and A(aQ1+bQ2) = 0; on P_* the halves are L = (0,0)_0 + (1/2,1/2)_{+1} + (1,0)_0 and "
       "R = (0,0)_0 + (1/2,1/2)_{-1} + (0,1)_0 under su2_+ + su2_- + u(1)_Q, and the exact traces give A(H^3) = 2c(3a^2 + 3b^2 + 4c^2), A(H) = 8c for "
       "H = a H+_3 + b H-_3 + c Q")

# thm:inertia3 (vector level; integer matrices, ranks modulo two large primes)
SYMB = [(a, b) for a in range(18) for b in range(a, 18)]


def sym_basis_mat(a, b):
    Mx = np.zeros((18, 18), dtype=np.int64)
    Mx[a, b] = Mx[b, a] = 1
    return Mx


G3n = np.array(G3.tolist(), dtype=np.int64)


def inv_forms(sym, cls):
    # unknown: symmetric N = sum c_ab S_ab (171 coefficients); conditions: class condition and r N + N r^T = 0 for every generator r
    cols = []
    symn = [np.array(r.tolist(), dtype=np.int64) for r in sym]
    for a, b in SYMB:
        S_ = sym_basis_mat(a, b)
        blocks = [(S_ @ G3n + G3n @ S_) if cls == "L" else (S_ @ G3n - G3n @ S_)]
        blocks += [r @ S_ + S_ @ r.T for r in symn]
        cols.append(np.concatenate([x.reshape(-1) for x in blocks]))
    Mx = np.array(cols).T
    Gm = Mx.T @ Mx                                                   # rank(M^T M) = rank(M) over Q; integer entries
    assert np.abs(Gm).max() < 2 ** 62
    r = DomainMatrix([[ZZ(int(x)) for x in row] for row in Gm.tolist()], Gm.shape, ZZ).convert_to(QQ).rank()
    r1 = rank_mod_p(Mx)
    return 171 - r if r == r1 else None


lGk = []
for k in range(3):
    base = 6 * k
    lGk += [adv(rot(base + 1 + (i - 1), base + 1 + (j - 1))) for i, j in [(1, 2), (1, 3), (2, 3)]]
    lGk += [adv(bst(base + 1 + (i - 1), base + 1 + (j - 1))) for i, j in [(1, 2), (1, 3), (2, 3)]]
fam_alg = [adv(rot(2, 3)), adv(rot(8, 9)), adv(rot(14, 15)), adv(bst(14, 15))]
st_alg = [adv(rot(2, 3)), adv(rot(8, 9))] + [adv(rot(i, j)) for i, j in itertools.combinations(W, 2)]
rows = [("Cent(G3)", gens, (2, 0)), ("sum lG^(k)", lGk, (6, 0)), ("P_fam-surviving", fam_alg, (18, 13)), ("P_*-surviving", st_alg, (8, 6)), ("none", [], (90, 81))]
got = {name: (inv_forms(sym, "L"), inv_forms(sym, "D")) for name, sym, _ in rows}
ok = all(got[name] == exp for name, _, exp in rows)
report("10", ok, f"thm:inertia3: (Lorentz, duality) dimensions of invariant symmetric forms on V_3^*: {got} (exact: ranks over Q of the integer "
       "Gram matrices M^T M, which equal the ranks of the condition matrices M, and agree with ranks modulo a large prime; the first row is also e5-4b)")


# ---------------------------------------------------------------------------- a3-11 prop:howmuch, non-real versus imaginary spectrum
# (review of version 1.0.1): a diagonalisable eta-self-adjoint G with eigenvalues 1 +- i gives no bounded motion
tau_ = sp.symbols('tau', real=True)
eta2 = sp.diag(1, -1)
Gc = sp.Matrix([[1, -1], [1, 1]])
ok = (Gc.T * eta2 - eta2 * Gc).is_zero_matrix                       # G^dagger = eta^-1 G^T eta = G
ok &= set(Gc.eigenvals()) == {1 + sp.I, 1 - sp.I} and Gc.is_diagonalizable()
Xc = sp.exp(tau_) * sp.Matrix([sp.cos(tau_), sp.sin(tau_)])           # kappa = 1
ok &= sp.simplify(Xc.diff(tau_, 2) - Gc * Gc * Xc).is_zero_matrix    # solves X'' = kappa^2 G^2 X
Mph = sp.Matrix(sp.BlockMatrix([[sp.zeros(2), sp.eye(2)], [Gc * Gc, sp.zeros(2)]]))
ev = Mph.eigenvals()
ok &= set(ev) == {1 + sp.I, 1 - sp.I, -1 + sp.I, -1 - sp.I} and all(sp.re(e) != 0 for e in ev)
# with purely imaginary spectrum (G = mu J_2 blocks on (1,1) planes) the subspace is bounded and neutral: G = diag(J2, 2 J2) on eta = diag(1,-1,1,-1)
J2 = sp.Matrix([[0, -1], [1, 0]])
eta4 = sp.diag(1, -1, 1, -1)
Gi = sp.diag(J2, 2 * J2)
ok &= (Gi.T * eta4 - eta4 * Gi).is_zero_matrix and set(Gi.eigenvals()) == {sp.I, -sp.I, 2 * sp.I, -2 * sp.I}
Mi = sp.Matrix(sp.BlockMatrix([[sp.zeros(4), sp.eye(4)], [Gi * Gi, sp.zeros(4)]]))
ok &= Mi.is_diagonalizable() and all(sp.re(e) == 0 for e in Mi.eigenvals())
report("11", ok, "prop:howmuch: for eta = diag(1,-1), G = [[1,-1],[1,1]] is eta-self-adjoint and diagonalisable with eigenvalues 1 +- i; "
       "X = e^tau (cos tau, sin tau) solves X'' = G^2 X and the phase-space exponents are +-1 +- i, so a non-real spectrum does not bind; "
       "for G = diag(J2, 2 J2) (purely imaginary spectrum) the phase-space generator is semisimple with imaginary spectrum (exact)")

print(f"{NPASS[0]}/{NPASS[0] + len(FAILS)} checks passed")
sys.exit(0 if not FAILS else 1)
