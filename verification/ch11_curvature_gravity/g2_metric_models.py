#!/usr/bin/env python3
"""g2: the explicit metric models of Chapter 11 (Phase 9C).

    python3 verification/ch11_curvature/g2_metric_models.py

Independent computation (Christoffel symbols -> Riemann -> Ricci -> scalar -> Einstein, curvature.py) of the statements of Section 11.1:
prop:kk, thm:ein6, cor:ein6, ex:nullwave, the two-modulus formula, the off-slice section metric and ex:landau, and of thm:ghostvec and
rem:signconv. Exact symbolic arithmetic (sympy). These are geometric calculations for given metrics; none of them solves or derives a
gravitational field equation.
"""
import sys

import sympy as sp

from curvature import Geometry, frame_metric, maxwell

NPASS, FAILS = [0], []


def report(cid, ok, text):
    if ok:
        NPASS[0] += 1
        print(f"PASS g2-{cid}: {text}")
    else:
        FAILS.append(cid)
        print(f"FAIL g2-{cid}: {text}")
    sys.stdout.flush()


def zero(M):
    return all(sp.simplify(e) == 0 for e in (M if isinstance(M, (list, tuple)) else list(M)))


eta = sp.diag(1, -1, -1, -1)

# ---------------------------------------------------------------------- prop:kk
t, x, y, z, u, v = sp.symbols("t x y z u v", real=True)
B = sp.symbols("B", real=True)
X4 = [t, x, y, z]


def kk_metric(h, A):
    """ds^2 = h + (du + v A)^2 + (dv - u A)^2 in coordinates (t,x,y,z,u,v), A = A_mu dx^mu; returns (g, g^-1)."""
    return frame_metric(h, [(1, [v * a for a in A]), (1, [-u * a for a in A])])


g6, g6i = kk_metric(eta, [0, 0, 0, B * x])
Gk = Geometry(g6, X4 + [u, v], ginv=g6i)
on = {u: 0, v: 0}
sig = sp.Matrix(g6.subs(on)).eigenvals()
ok = sorted(sp.sign(e) for e, mlt in sig.items() for _ in range(mlt)) == [-1, -1, -1, 1, 1, 1]
ok &= sp.simplify(Gk.riemann(4, 5, 1, 3).subs(on) - B) == 0
ok &= all(sp.simplify(Gk.riemann(a, b, c, d).subs(on)) == 0 for a in range(4) for b in range(4) for c in range(4) for d in range(4))
# general A_mu(x) on a flat base: R^u_{v mu nu} on Sigma is (dA)_{mu nu}
Af = [sp.Function(f"A{i}")(*X4) for i in range(4)]
gA, gAi = kk_metric(eta, Af)
Gg = Geometry(gA, X4 + [u, v], simplify=sp.expand, ginv=gAi)
F_A = [[sp.diff(Af[n], X4[m]) - sp.diff(Af[m], X4[n]) for n in range(4)] for m in range(4)]
okA = all(sp.simplify(Gg.riemann(4, 5, m, n).subs(on) - F_A[m][n]) == 0 for m in range(4) for n in range(4) if m < n)
report("1", ok and okA, "prop:kk: ds^2 = eta + (du + vA)^2 + (dv - uA)^2 has signature (3,3); with the convention R^a_bcd = d_c Gamma^a_db - ..., "
       "R^u_{v mu nu} on Sigma = {u = v = 0} equals (dA)_{mu nu} for arbitrary A_mu(x) on a flat base, and for A = B x dz, R^u_{vxz} = B while every "
       "R^mu_{nu rho sigma} vanishes on Sigma (exact)")

# ---------------------------------------------------------------------- thm:ein6 on the Friedmann base of the proof
r, ph = sp.symbols("r phi", positive=True)
hF = sp.diag(1, -t ** 2, -t ** 2, -t ** 2)
AF = [z, 0, x ** 2 / 2, 0]
g6F, g6Fi = frame_metric(hF, [(1, [0, 0, 0, 0]), (r ** 2, [-a for a in AF])])   # coordinates (t,x,y,z,r,phi)
G6 = Geometry(g6F, X4 + [r, ph], ginv=g6Fi)
Ein6 = G6.einstein(); R6 = G6.scalar()
G4 = Geometry(hF, X4)
Ein4 = G4.einstein(); R4 = G4.scalar()
F, F2, W = maxwell(hF, AF, X4)
hi = hF.inv()
sq = sp.sqrt(-hF.det())
Fuu = hi * F * hi                                                  # nabla_rho F^{rho nu} = d_rho F^{rho nu} + Gamma^rho_{rho s} F^{s nu}
divu = [sum(sp.diff(Fuu[rr, n], X4[rr]) for rr in range(4)) + sum(G4.Gam[rr][rr][s_] * Fuu[s_, n] for rr in range(4) for s_ in range(4))
        for n in range(4)]
divF = [sp.simplify(sum(hF[m, n] * divu[n] for n in range(4))) for m in range(4)]
Ev = [sp.Matrix([1 if j == m else 0 for j in range(4)] + [0, AF[m]]) for m in range(4)]   # E_mu = d_mu + A_mu d_phi
er = sp.Matrix([0, 0, 0, 0, 1, 0]); ep = sp.Matrix([0, 0, 0, 0, 0, 1 / r])
Gfr = lambda U, V: sp.simplify((U.T * Ein6 * V)[0])
ok = all(sp.simplify(Gfr(Ev[m], Ev[n]) - (Ein4[m, n] - r ** 2 / 2 * W[m, n])) == 0 for m in range(4) for n in range(4))
ok &= sp.simplify(Gfr(er, er) - (-R4 / 2 + r ** 2 * F2 / 8)) == 0
ok &= sp.simplify(Gfr(ep, ep) - (-R4 / 2 + 3 * r ** 2 * F2 / 8)) == 0
sgn = [sp.simplify(Gfr(Ev[m], ep) / (r / 2 * divF[m])) if divF[m] != 0 else None for m in range(4)]
ok_mixed = all(sp.simplify(Gfr(Ev[m], ep) - r / 2 * divF[m]) == 0 for m in range(4))
ok &= all(sp.simplify(Gfr(Ev[m], er)) == 0 for m in range(4)) and sp.simplify(Gfr(er, ep)) == 0
ok &= sp.simplify(R6 - (R4 - r ** 2 * F2 / 4)) == 0
report("2", ok and ok_mixed, f"thm:ein6 on the Friedmann base h = diag(1,-t^2,-t^2,-t^2), A = z dt + x^2 dy/2 (electric and magnetic, non-constant F): in the frame "
       "E_mu = d_mu + A_mu d_phi, d_r, r^-1 d_phi the six-dimensional Einstein tensor is G_mn[h] - (r^2/2) W_mn, G_rr = -R/2 + r^2F^2/8, "
       "G_phiphi = -R/2 + 3r^2F^2/8, G_{mu phi} = (r/2) nabla^rho F_{rho mu}, G_{mu r} = G_{r phi} = 0, R^(6) = R[h] - r^2F^2/4, for all r > 0 "
       f"(exact; mixed-component ratios {sgn})")
F2s = sp.symbols("F2"); Rh = sp.symbols("R_h"); rc = sp.symbols("r_c", positive=True)
ok = sp.simplify(sp.integrate(sp.integrate((Rh - r ** 2 * F2s / 4) * r, (r, 0, rc)), (ph, 0, 2 * sp.pi))
                 - (sp.pi * rc ** 2 * Rh - sp.pi * rc ** 4 / 8 * F2s)) == 0
# -W_00 = (|E|^2 + |B|^2)/2 for a general constant field in (+,-,-,-); variation of sqrt|g| F^2 at eta gives 2 W
f01, f02, f03, f12, f13, f23 = sp.symbols("f01 f02 f03 f12 f13 f23", real=True)
Fc = sp.Matrix([[0, f01, f02, f03], [-f01, 0, f12, f13], [-f02, -f12, 0, f23], [-f03, -f13, -f23, 0]])
Fcu = eta.inv() * Fc * eta.inv()
F2c = sum(Fc[m, n] * Fcu[m, n] for m in range(4) for n in range(4))
Wc = sp.Matrix(4, 4, lambda m, n: sum(Fc[m, rr] * (Fc * eta.inv())[n, rr] for rr in range(4)) - eta[m, n] * F2c / 4)
ok &= sp.expand(-Wc[0, 0] - (f01 ** 2 + f02 ** 2 + f03 ** 2 + f12 ** 2 + f13 ** 2 + f23 ** 2) / 2) == 0
dl = sp.symbols("d0:10"); eps_ = sp.symbols("epsilon")
Dm = sp.zeros(4, 4); n_ = 0
for i in range(4):
    for j in range(i, 4):
        Dm[i, j] = Dm[j, i] = dl[n_]; n_ += 1
ginv_e = eta.inv() + eps_ * Dm
g_e1 = eta - eps_ * eta * Dm * eta                                  # (eta^-1 + eps D)^-1 to first order in eps (exact at that order)
F2e = sum(Fc[m, n] * (ginv_e * Fc * ginv_e)[m, n] for m in range(4) for n in range(4))
mdet = sp.expand(-g_e1.det())
dsq = sp.diff(mdet, eps_).subs(eps_, 0) / (2 * sp.sqrt(mdet.subs(eps_, 0)))   # d sqrt|g| / d eps at 0
var = sp.expand(dsq * F2e.subs(eps_, 0) + sp.diff(F2e, eps_).subs(eps_, 0))
ok &= sp.expand(var - 2 * sum(Wc[m, n] * Dm[m, n] for m in range(4) for n in range(4))) == 0
report("3", ok, "cor:ein6: (ii) the integral of r R^(6) over a normal disc of radius r_c is pi r_c^2 R - pi r_c^4 F^2/8, and the variation of "
       "sqrt|h| F^2 with respect to h^{mn} is 2 sqrt|h| W_mn, so the Euler-Lagrange equation of sqrt|h|(a1 R + a2 F^2) is a1 G + 2 a2 W = 0 and the "
       "reduced ratio a2/a1 = -r_c^2/8 has the opposite sign to Einstein-Maxwell (+1); (i) -W_00 = (|E|^2 + |B|^2)/2 in (+,-,-,-) for a general "
       "constant field (exact). Rigidity (iii) follows from g2-2: G_rr - G_phiphi = -r^2 F^2/4")

# ---------------------------------------------------------------------- ex:nullwave
U, V, Xn, Yn, psi, chi = sp.symbols("U V X Y psi chi", real=True)
ell, mm = sp.symbols("ell m", positive=True)
Hf = sp.Function("H")(U, Xn, Yn); af = sp.Function("a")(U)
hN = sp.Matrix([[-Hf, 1, 0, 0], [1, 0, 0, 0], [0, 0, -1, 0], [0, 0, 0, -1]])
g6N, g6Ni = frame_metric(hN, [(ell ** 2, [-af * Xn, 0, 0, 0]), (mm ** 2, [0, 0, 0, 0])])
GN = Geometry(g6N, [U, V, Xn, Yn, psi, chi], simplify=sp.simplify, ginv=g6Ni)
EN = GN.einstein()
target = -sp.Rational(1, 2) * (sp.diff(Hf, Xn, 2) + sp.diff(Hf, Yn, 2) - ell ** 2 * af ** 2)
ok = sp.simplify(EN[0, 0] - target) == 0 and all(sp.simplify(EN[i, j]) == 0 for i in range(6) for j in range(6) if (i, j) != (0, 0))
Hsol = ell ** 2 * af ** 2 * (Xn ** 2 + Yn ** 2) / 4
ok &= sp.simplify(target.subs(Hf, Hsol).doit()) == 0
G4N = Geometry(hN.subs(Hf, Hsol), [U, V, Xn, Yn], simplify=sp.simplify).einstein()
_, _, WN = maxwell(hN, [af * Xn, 0, 0, 0], [U, V, Xn, Yn])
ok &= sp.simplify(G4N[0, 0] + ell ** 2 * af ** 2 / 2) == 0 and sp.simplify(-2 * WN[0, 0] - 2 * af ** 2) == 0
sig = sp.Matrix(g6N.subs({Hf: 0, af: 0})).eigenvals()
ok &= sum(m_ for e_, m_ in sig.items() if e_.is_positive) == 3
report("4", ok, "ex:nullwave: for h = 2dUdV - H dU^2 - dX^2 - dY^2, A = a(U) X dU and the fibre l^2(dpsi - A)^2 + m^2 dchi^2 (signature (3,3)) the only "
       "nonzero component of the six-dimensional Einstein tensor is G_UU = -(1/2)(H_XX + H_YY - l^2 a^2); H = l^2 a^2 (X^2+Y^2)/4 solves the "
       "six-dimensional vacuum with F != 0; there G_UU[h] = -l^2 a^2/2 < 0 while a Maxwell null wave has -2 W_UU = +2a^2 > 0 (exact, general H, a)")

# ---------------------------------------------------------------------- two moduli: G_psipsi - G_chichi
al = sp.Function("alpha")(t, x); be = sp.Function("beta")(t, x); Az = sp.Function("Az")(t, x)
gM, gMi = frame_metric(eta, [(sp.exp(2 * al), [0, 0, 0, -Az]), (sp.exp(2 * be), [0, 0, 0, 0])])
GM = Geometry(gM, X4 + [psi, chi], simplify=sp.simplify, ginv=gMi)
EM = GM.einstein()
lhs = sp.simplify(sp.exp(-2 * al) * EM[4, 4] - sp.exp(-2 * be) * EM[5, 5])
_, F2M, _ = maxwell(eta, [0, 0, 0, Az], X4)
dd = al - be; ss = al + be
box = sum(eta.inv()[m, m] * sp.diff(dd, X4[m], 2) for m in range(4))
grad = sum(eta.inv()[m, m] * sp.diff(dd, X4[m]) * sp.diff(ss, X4[m]) for m in range(4))
ok = sp.simplify(lhs - (sp.exp(2 * al) * F2M / 4 - box - grad)) == 0
report("5", ok, "two moduli: for ds^2 = eta + e^{2 alpha}(dpsi - A)^2 + e^{2 beta} dchi^2 with alpha, beta, A_z depending on (t,x), "
       "G_psipsi - G_chichi (orthonormal) = (1/4) e^{2 alpha} F^2 - Box(alpha - beta) - grad(alpha - beta).grad(alpha + beta) (exact)")

# ---------------------------------------------------------------------- off the slice, ex:landau
hs = sp.symbols("h0:10"); bs = sp.symbols("b0:4"); r0 = sp.symbols("r_0", positive=True)
Hm = sp.zeros(4, 4); n_ = 0
for i in range(4):
    for j in range(i, 4):
        Hm[i, j] = Hm[j, i] = hs[n_]; n_ += 1
Bv = sp.Matrix(bs)
ok = sp.simplify((Hm + r0 ** 2 * Bv * Bv.T).det() - Hm.det() * (1 + r0 ** 2 * (Bv.T * Hm.inv() * Bv)[0])) == 0
kq, bq, Rc = sp.symbols("k b R_c", positive=True)
hL = sp.diag(1, -1, -1, -(1 - kq ** 2 * x ** 2))
RL = Geometry(hL, X4, simplify=sp.simplify).scalar()
ok &= sp.simplify(RL + 2 * kq ** 2 / (1 - kq ** 2 * x ** 2) ** 2) == 0
th_ = sp.symbols("vartheta")
xc, zc = Rc * sp.cos(th_), Rc * sp.sin(th_)
circ = sp.integrate((-bq * xc) * sp.diff(zc, th_), (th_, 0, 2 * sp.pi))  # loop integral of B = dchi - A, A = b x dz
ok &= sp.simplify(circ + sp.pi * bq * Rc ** 2) == 0
report("6", ok, "off the slice: det(h + r_0^2 B B) = det h (1 + r_0^2 h^{mn} B_m B_n) for a general symmetric h (exact); ex:landau: the section metric "
       "dt^2 - dx^2 - dy^2 - (1 - k^2x^2) dz^2 has R = -2k^2/(1 - k^2x^2)^2, and the circulation of B = dchi - A around a circle of radius R_c in an "
       "x-z plane is -pi b R_c^2 whatever chi (exact); the bound sum B_i^2 >= b^2 R_c^2/4 somewhere is then Cauchy-Schwarz on the loop")

# ---------------------------------------------------------------------- thm:ghostvec
sgm = sp.symbols("sigma")
b1, b2 = sp.symbols("b1 b2", real=True)
x0, x1, x2, x3 = X4
gG, gGi = frame_metric(eta, [(sgm, [0, 0, b1 * x1, 0]), (sgm, [0, 0, 0, b2 * x0])])
RG = Geometry(gG, X4 + [u, v], simplify=sp.simplify, ginv=gGi).scalar()
FF = 2 * b1 ** 2 * eta[1, 1] * eta[2, 2] + 2 * b2 ** 2 * eta[0, 0] * eta[3, 3]
ok = sp.simplify(RG + sgm / 4 * FF) == 0
a1, a2 = sp.symbols("a1 a2", real=True)
Arot = [a1 * x1, 0, 0, a2 * x1]
gR, gRi = frame_metric(eta, [(sgm, [v * a for a in Arot]), (sgm, [-u * a for a in Arot])])
RR = Geometry(gR, X4 + [u, v], simplify=sp.simplify, ginv=gRi).scalar()
_, F2r, _ = maxwell(eta, Arot, X4)
ok &= sp.simplify((RR + sgm * (u ** 2 + v ** 2) * F2r / 4).subs(sgm ** 2, 1)) == 0
report("7", ok, "thm:ghostvec: for ds^2 = eta + sigma[(du + b1 x^1 dx^2)^2 + (dv + b2 x^0 dx^3)^2], R^(6) = -(sigma/4) sum F_mn F^mn, and for the rotation "
       "ansatz eta + sigma[(du + vA)^2 + (dv - uA)^2], R^(6) = -sigma r^2 F^2/4 (sigma = +-1), odd in sigma (exact). With thm:selector(ii) this "
       "is the sign statement; it concerns the scalar curvature of given metrics, not a derived reduced action")

# ---------------------------------------------------------------------- rem:signconv
Gp = Geometry(hF, X4); Gm = Geometry(-hF, X4)
_, _, Wp = maxwell(hF, AF, X4); _, F2m, Wm = maxwell(-hF, AF, X4)
ok = sp.simplify(Gm.scalar() + Gp.scalar()) == 0 and zero(Gm.einstein() - Gp.einstein()) and zero(Wm + Wp) and sp.simplify(F2m - F2) == 0
report("8", ok, "rem:signconv: under g -> -g (on the Friedmann example) R and W_mn change sign while the Einstein tensor and F^2 do not; so the healthy "
       "ratio a2/a1 of a1 R + a2 F^2 changes sign with the convention, and so does the reduced ratio: the mismatch is convention independent (exact)")

print(f"{NPASS[0]}/{NPASS[0] + len(FAILS)} checks passed")
sys.exit(0 if not FAILS else 1)
