"""Phase 9D, batch 2: adversarial audit of the Chapter 3 statements that had neither a cited check nor a written proof.

Exact symbolic arithmetic (sympy) on the 8x8 realisation of Appendix A and on explicit field configurations.
Each check prints PASS/FAIL a2-<n>; the script exits 0 iff all pass.
"""
import itertools
import os
import sys
from functools import reduce

import sympy as sp

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "common"))
from clifford8 import I8, bl, comm, e, vec_coords  # noqa: E402

NPASS, FAILS = [0], []


def report(cid, ok, text):
    if ok:
        NPASS[0] += 1
        print(f"PASS a2-{cid}: {text}")
    else:
        FAILS.append(cid)
        print(f"FAIL a2-{cid}: {text}")
    sys.stdout.flush()


Z8 = sp.zeros(8, 8)
eps = sp.LeviCivita
cyc = [(1, 2, 3), (2, 3, 1), (3, 1, 2)]


def ad6(Xb):
    """6x6 matrix of v -> [Xb, v] on the vectors e1..e6."""
    return sp.Matrix([vec_coords(comm(Xb, e[a])) for a in range(1, 7)]).T


def nullspace_dim(conds, n):
    """dimension of the solution space of linear conditions on n x n matrices (list of functions X -> matrix)."""
    xs = sp.symbols(f"q0:{n * n}")
    X = sp.Matrix(n, n, xs)
    eqs = [x for c in conds for x in c(X)]
    A, _ = sp.linear_eq_to_matrix(eqs, xs)
    return n * n - A.rank(), A, xs


def solve_mats(conds, n):
    d, A, xs = nullspace_dim(conds, n)
    return [sp.Matrix(n, n, list(v)) for v in A.nullspace()]


Lb = {i: (e[j + 3] * e[k + 3] - e[j] * e[k]) / 2 for i, j, k in cyc}
Kb = {i: (e[k] * e[j + 3] - e[j] * e[k + 3]) / 2 for i, j, k in cyc}
L6 = {i: ad6(Lb[i]) for i in Lb}
K6 = {i: ad6(Kb[i]) for i in Kb}
I3, Z3 = sp.eye(3), sp.zeros(3)
G6 = sp.BlockMatrix([[Z3, -I3], [I3, Z3]]).as_explicit()
eta6 = sp.diag(1, 1, 1, -1, -1, -1)
lG6 = list(L6.values()) + list(K6.values())
gam = [e[1], e[4], e[5], e[6]]
LINT8 = [gam[m] * gam[n] for m, n in itertools.combinations(range(4), 2)]
lint6 = [ad6(x) for x in LINT8]

# ---------------------------------------------------------------- a2-1: thm:emfield, ex:emboost, prop:notconj, cor:Gselected, prop:emexclusive
ok = True
for i in (1, 2, 3):
    for j in (1, 2, 3):
        sL = sum((eps(i, j, k) * Lb[k] for k in (1, 2, 3)), Z8)
        sK = sum((eps(i, j, k) * Kb[k] for k in (1, 2, 3)), Z8)
        ok &= comm(Lb[i], Lb[j]) == sL and comm(Lb[i], Kb[j]) == sK and comm(Kb[i], Kb[j]) == -sL
ok &= all(comm(G6, x) == sp.zeros(6) and x.T * eta6 + eta6 * x == sp.zeros(6) for x in lG6)
Bs, Es, vs = sp.Matrix(sp.symbols("B1:4")), sp.Matrix(sp.symbols("E1:4")), sp.Matrix(sp.symbols("v1:4"))
X6 = sp.Matrix.vstack(Bs, Es)
rk = {k: sp.Matrix(3, 3, lambda a, c: -eps(k, a + 1, c + 1)) for k in (1, 2, 3)}  # standard rotation generator, r_k u = e_k x u
ok &= all(L6[k] * X6 == sp.Matrix.vstack(rk[k] * Bs, rk[k] * Es) for k in (1, 2, 3))
dX = sum((vs[k - 1] * K6[k] for k in (1, 2, 3)), sp.zeros(6)) * X6
ok &= sp.expand(dX[3:, 0] - vs.cross(Bs)) == sp.zeros(3, 1) and sp.expand(dX[:3, 0] + vs.cross(Es)) == sp.zeros(3, 1)
# G is the Hodge star: F_X with F_0i = E_i, F_jk = eps_ijk B_i; (*F)_ab = 1/2 eps_abcd F^cd, eps_0123 = +1


def Fmat(B, E):
    F = sp.zeros(4)
    for i in range(3):
        F[0, i + 1], F[i + 1, 0] = E[i], -E[i]
    for i in range(3):
        for j in range(3):
            for k in range(3):
                F[j + 1, k + 1] += eps(i + 1, j + 1, k + 1) * B[i]
    return F


def hodge(F, eta4):
    Fu = eta4 * F * eta4
    return sp.Matrix(4, 4, lambda a, b: sum(eps(a, b, c, d) * Fu[c, d] for c in range(4) for d in range(4)) / 2)


GX = G6 * X6
hodge_ok = {}
for name, eta4 in [("(+,-,-,-)", sp.diag(1, -1, -1, -1)), ("(-,+,+,+)", sp.diag(-1, 1, 1, 1))]:
    hodge_ok[name] = sp.expand(Fmat(GX[:3, 0], GX[3:, 0]) - hodge(Fmat(Bs, Es), eta4)) == sp.zeros(4)
ok &= all(hodge_ok.values())
ok &= sp.expand((X6.T * eta6 * X6)[0] - (Bs.dot(Bs) - Es.dot(Es))) == 0 and sp.expand((X6.T * eta6 * GX)[0] + 2 * Es.dot(Bs)) == 0
quad = solve_mats([(lambda x: (lambda Q: list(x.T * Q + Q * x)))(x) for x in lG6] + [lambda Q: list(Q - Q.T)], 6)
ok &= len(quad) == 2 and sp.Matrix([list(q) for q in quad]).rank() == sp.Matrix([list(q) for q in quad] + [list(eta6), list(eta6 * G6 + (eta6 * G6).T)]).rank()
comG = solve_mats([(lambda x: (lambda Q: list(x * Q - Q * x)))(x) for x in lG6], 6)
comI = solve_mats([(lambda x: (lambda Q: list(x * Q - Q * x)))(x) for x in lint6], 6)
ok &= len(comG) == 2 and len(comI) == 5
fixG = sp.Matrix.vstack(*lG6).nullspace()
fixI = sp.Matrix.vstack(*lint6).nullspace()
ok &= len(fixG) == 0 and len(fixI) == 2 and all(v[0] == v[3] == v[4] == v[5] == 0 for v in fixI)
ok &= comm(Kb[3], e[1]) == e[5] and comm(Kb[3], e[5]) == e[1] and Kb[3] == (e[2] * e[4] - e[1] * e[5]) / 2
inter = sp.Matrix([list(x) for x in LINT8 + list(Lb.values()) + list(Kb.values())]).rank()
ok &= inter == 12
report("1", ok, "thm:emfield, ex:emboost, prop:notconj, cor:Gselected, prop:emexclusive: the brackets of L_i, K_i; B, E rotate as vectors; "
       "delta E = v x B, delta B = -v x E; F_{GX} = *F_X for both metric signatures; eta(X,X) = |B|^2-|E|^2 and eta(X,GX) = -2 E.B span the "
       "invariant symmetric forms (dimension 2); commutants of l_G and l_int on R^6 of dimensions 2 and 5; common fixed spaces of dimensions 0 and "
       "2 (the surplus plane); [K_3,e1] = e5, [K_3,e5] = e1; l_int and l_G meet in zero (exact)")

# ---------------------------------------------------------------- a2-2: prop:bivectorF
Om = gam[0] * gam[1] * gam[2] * gam[3]
L2E = LINT8
span = sp.Matrix([list(x) for x in L2E])
ok = all(span.rank() == sp.Matrix([list(x) for x in L2E] + [list(comm(a, b))]).rank() for a in LINT8 for b in L2E)
tf = sp.Matrix(6, 6, lambda i, j: (L2E[i] * L2E[j]).trace() / 8)
evs = tf.eigenvals()
ok &= sum(m for v, m in evs.items() if v > 0) == 3 and sum(m for v, m in evs.items() if v < 0) == 3
ok &= Om * Om == -I8 and all(comm(Om, x) == Z8 for x in LINT8)
ok &= all(span.rank() == sp.Matrix([list(x) for x in L2E] + [list(Om * b)]).rank() for b in L2E)
# Omega gamma^0 gamma^a is a spatial bivector
ok &= all(any(Om * gam[0] * gam[a] == s * gam[b] * gam[c] for s in (1, -1) for b, c in itertools.combinations(range(1, 4), 2)) for a in (1, 2, 3))
eta4 = sp.diag(1, -1, -1, -1)
fs = sp.symbols("f0:6")
pairs = list(itertools.combinations(range(4), 2))
Fl = sp.zeros(4)
for (a, b), f in zip(pairs, fs):
    Fl[a, b], Fl[b, a] = f, -f


def biv_of(F):  # 1/2 F_ab gamma^a gamma^b
    return sum((F[a, b] * gam[a] * gam[b] for a in range(4) for b in range(4)), Z8) / 2


hs = hodge(Fl, eta4)
sign_star = [s for s in (1, -1) if sp.expand(Om * biv_of(Fl) - s * biv_of(hs)) == Z8]
ok &= len(sign_star) == 1
report("2", ok, "prop:bivectorF: Lambda^2 E is ad(l_int)-closed; (1/8)Tr(XY) has signature (3,3); Omega^2 = -1 commutes with l_int, preserves Lambda^2 E by "
       f"left multiplication, sends gamma^0 gamma^a to spatial bivectors, and acts as (corrected text: minus the Hodge star) {'+' if sign_star == [1] else '-'}* (eta = (+,-,-,-), eps_0123 = +1) (exact)")

# ---------------------------------------------------------------- a2-3: thm:twoops
Lor4 = []
for a, b in pairs:  # generators of so(1,3) on V: (M_ab)^c_d = delta^c_a eta_bd - delta^c_b eta_ad
    M = sp.zeros(4)
    for c in range(4):
        for d in range(4):
            M[c, d] = (1 if c == a else 0) * eta4[b, d] - (1 if c == b else 0) * eta4[a, d]
    Lor4.append(M)
# Lambda^2 (lower indices) basis pairs; generator action on lower-index tensors is -M^T
def lam2(M):
    out = sp.zeros(6)
    for col, (a, b) in enumerate(pairs):
        T = sp.zeros(4)
        T[a, b], T[b, a] = 1, -1
        Tn = -M.T * T - T * M
        for row, (c, d) in enumerate(pairs):
            out[row, col] = Tn[c, d]
    return out


def kron(A, B):
    return sp.kronecker_product(A, B)


ops24 = [kron(lam2(M), sp.eye(4)) + kron(sp.eye(6), -M.T) for M in Lor4]  # on Lambda^2 (x) V^* (the derivative index is lower)
xs = sp.symbols("h0:96")
Hm = sp.Matrix(4, 24, xs)
eqs = []
for M, O in zip(Lor4, ops24):
    eqs += list((-M.T) * Hm - Hm * O)  # target: lower-index vector (nu)
A_, _ = sp.linear_eq_to_matrix(eqs, xs)
homdim = 96 - A_.rank()
cas = sum(((-1 if (a == 0) != (b == 0) else 1) * O * O for O, (a, b) in zip(ops24, pairs)), sp.zeros(24))
cev = {sp.nsimplify(k): v for k, v in cas.eigenvals().items()}
ok = homdim == 2 and sorted(cev.values()) == [8, 16]
t, x1, x2, x3 = sp.symbols("t x1 x2 x3", real=True)
co = (t, x1, x2, x3)
Ff = sp.zeros(4)
for (a, b) in pairs:
    f = sp.Function(f"F{a}{b}")(*co)
    Ff[a, b], Ff[b, a] = f, -f
Fu = eta4 * Ff * eta4
delta = lambda F: [sp.expand(sum(sp.diff((eta4 * F * eta4)[m, n], co[m]) for m in range(4))) for n in range(4)]  # noqa: E731
dF = lambda F: {(m, n, r): sp.diff(F[n, r], co[m]) + sp.diff(F[r, m], co[n]) + sp.diff(F[m, n], co[r]) for m, n, r in itertools.combinations(range(4), 3)}  # noqa: E731
dstar = delta(hodge(Ff, eta4))
dd = dF(Ff)
ratios = set()
for n in range(4):
    comp = [mrs for mrs in dd if n not in mrs][0]
    for s_ in (1, -1):
        if sp.expand(dstar[n] - s_ * dd[comp]) == 0:
            ratios.add((n, s_))
ok &= len(ratios) == 4 and len({n for n, _ in ratios}) == 4
report("3", ok, "thm:twoops: Hom_{so(1,3)}(Lambda^2 (x) V, V) has dimension 2 (equivariance solved over all 96 coefficients); the Casimir on Lambda^2 (x) V "
       "has two eigenspaces of dimensions 16 and 8 ((3/2,1/2)+(1/2,3/2) and 2(1/2,1/2)); delta(*F) = +-*dF componentwise, so d and delta realise the two "
       "copies (exact)")

# ---------------------------------------------------------------- a2-4: prop:RS, eq:RSsource, continuity; cor:maxwellform duality
Bf = sp.Matrix([sp.Function(f"B{i}")(*co) for i in (1, 2, 3)])
Ef = sp.Matrix([sp.Function(f"E{i}")(*co) for i in (1, 2, 3)])
X = (x1, x2, x3)
curl = lambda v: sp.Matrix([sp.diff(v[2], X[1]) - sp.diff(v[1], X[2]), sp.diff(v[0], X[2]) - sp.diff(v[2], X[0]), sp.diff(v[1], X[0]) - sp.diff(v[0], X[1])])  # noqa: E731
div = lambda v: sum(sp.diff(v[i], X[i]) for i in range(3))  # noqa: E731
FX = Fmat(Bf, Ef)
dFX = dF(FX)
dSX = dF(hodge(FX, eta4))
# expected: dF = 0 <=> div B = 0 and dt B = curl E;  d*F = 0 <=> div E = 0 and dt E = -curl B
faraday = sp.Matrix([sp.diff(Bf[i], t) for i in range(3)]) - curl(Ef)
ampere = sp.Matrix([sp.diff(Ef[i], t) for i in range(3)]) + curl(Bf)
vals_dF = [sp.expand(v) for v in dFX.values()]
vals_dS = [sp.expand(v) for v in dSX.values()]
def same_span(vals, target):
    M1 = sp.Matrix([[sp.expand(v).coeff(sym) for sym in basis_syms] for v in vals])
    M2 = sp.Matrix([[sp.expand(v).coeff(sym) for sym in basis_syms] for v in target])
    return M1.rank() == M2.rank() == sp.Matrix.vstack(M1, M2).rank()


derivs = sorted({d for v in vals_dF + vals_dS + list(faraday) + list(ampere) for d in sp.expand(v).atoms(sp.Derivative)}, key=str)
basis_syms = derivs
ok = same_span(vals_dF, list(faraday) + [div(Bf)]) and same_span(vals_dS, list(ampere) + [div(Ef)])
rho_, j_ = sp.Function("rho")(*co), sp.Matrix([sp.Function(f"j{i}")(*co) for i in (1, 2, 3)])
src = delta(FX)  # d_mu F^{mu nu}
ok &= sp.expand(src[0] - div(Ef)) == 0 and all(sp.expand(src[i + 1] - (-sp.diff(Ef[i], t) - curl(Bf)[i])) == 0 for i in range(3))
z = Bf + sp.I * Ef
# with d_mu F^{mu nu} = (rho, j): div E = rho, dt E = -curl B - j; then dt z = -i curl z - i j, div z = i rho
sub = {sp.Derivative(Ef[i], t): -curl(Bf)[i] - j_[i] for i in range(3)}
sub.update({sp.Derivative(Bf[i], t): curl(Ef)[i] for i in range(3)})
lhs = sp.Matrix([sp.diff(z[i], t) for i in range(3)]).subs(sub)
ok &= sp.expand(lhs - (-sp.I * curl(z) - sp.I * j_)) == sp.zeros(3, 1)
# continuity: d_t rho = d_t div E = div(-curl B - j) = -div j
ok &= sp.expand(div(sp.Matrix([-curl(Bf)[i] - j_[i] for i in range(3)])) + div(j_)) == 0
ok &= all(hodge_ok.values())
report("4", ok, "prop:RS, eq:RSsource, cor:maxwellform: for F_X, dF = 0 is {dt B = curl E, div B = 0} and d*F = 0 is {dt E = -curl B, div E = 0} (equal "
       "spans of the component equations); d_mu F^{mu nu} = (div E, -dt E - curl B), so with sources dt z = -i curl z - i j, div z = i rho, and "
       "the divergence of the first gives continuity; G = * exchanges the two operators (exact)")

# ---------------------------------------------------------------- a2-5: ex:circular, ex:threewaves, prop:orthframe, ex:counter, prop:EBlocal
k = sp.symbols("k", positive=True)
ph = k * (x3 - t)


def BE_of(zz):
    return sp.Matrix([sp.re(sp.expand_complex(c)) for c in zz]), sp.Matrix([sp.im(sp.expand_complex(c)) for c in zz])


def solves_RS(zz):
    return sp.simplify(sp.Matrix([sp.diff(c, t) for c in zz]) + sp.I * curl(zz)) == sp.zeros(3, 1) and sp.simplify(div(zz)) == 0


zc = sp.Matrix([1, sp.I, 0]) * sp.exp(sp.I * ph)
B_, E_ = BE_of(zc)
ok = solves_RS(zc) and sp.simplify(curl(zc) - k * zc) == sp.zeros(3, 1) and sp.simplify(zc.dot(zc)) == 0
pv = sp.simplify(sp.I / 2 * sp.Matrix([sp.conjugate(c) for c in zc]).cross(zc))
ok &= sp.simplify(pv - sp.Matrix([0, 0, -1])) == sp.zeros(3, 1) and sp.simplify(E_.cross(B_) - pv) == sp.zeros(3, 1)
ok &= sp.simplify(E_.cross(-B_) - sp.Matrix([0, 0, 1])) == sp.zeros(3, 1)
zs = [sp.Matrix([1, sp.I, 0]) * sp.exp(sp.I * ph), sp.Matrix([1, sp.I, 0]) * sp.exp(-sp.I * ph), sp.Matrix([1, sp.I, 0]) * sp.cos(ph)]
for zz, hel in zip(zs, (k, -k, 0)):
    B_, E_ = BE_of(zz)
    ok &= solves_RS(zz) and sp.simplify(B_.dot(E_)) == 0 and sp.simplify(B_.dot(B_) - E_.dot(E_)) == 0
    bxe = sp.simplify(B_.cross(E_))
    ok &= bxe[0] == 0 and bxe[1] == 0 and sp.simplify(bxe[2] - B_.dot(B_)) == 0
    Edot = sp.Matrix([sp.diff(c, t) for c in E_])
    ok &= sp.simplify(E_.cross(Edot)[2].subs(x3, 0) - hel * E_.dot(E_).subs(x3, 0)) == 0
zmix = sp.Matrix([1, sp.I, 0]) * sp.cos(ph) + sp.Matrix([0, 1, sp.I]) * sp.cos(k * (x1 - t))
B_, E_ = BE_of(zmix)
ok &= solves_RS(zmix) and sp.simplify(B_.dot(B_) - E_.dot(E_)) == 0 and sp.simplify(B_.dot(E_) - sp.cos(k * (x1 - t)) * sp.cos(ph)) == 0
zg = Bs + sp.I * Es
ok &= sp.expand(zg.dot(zg) - ((X6.T * eta6 * X6)[0] - sp.I * (X6.T * eta6 * GX)[0])) == 0
ok &= sp.expand((GX.T * eta6 * GX)[0] + (X6.T * eta6 * X6)[0]) == 0
# ex:counter (textbook variables)
z_ = x3
c1, c2 = sp.cos(k * z_ - k * t), sp.cos(k * z_ + k * t)
Ec, Bc = sp.Matrix([c1, c2, 0]), sp.Matrix([c2, c1, 0])
Ep, Bp = sp.Matrix([c1 + c2, 0, 0]), sp.Matrix([0, c1 - c2, 0])


def maxwell_tb(E, B):
    return (sp.simplify(div(E)) == 0 and sp.simplify(div(B)) == 0 and sp.simplify(sp.Matrix([sp.diff(c, t) for c in B]) + curl(E)) == sp.zeros(3, 1)
            and sp.simplify(sp.Matrix([sp.diff(c, t) for c in E]) - curl(B)) == sp.zeros(3, 1))


ok &= maxwell_tb(Ec, Bc) and maxwell_tb(Ep, Bp)
ok &= sp.simplify(Ec.dot(Bc) - (sp.cos(2 * k * z_) + sp.cos(2 * k * t))) == 0 and sp.simplify(Bc.dot(Bc) - Ec.dot(Ec)) == 0
ok &= sp.simplify(Ep.dot(Bp)) == 0 and sp.simplify(Bp.dot(Bp) - Ep.dot(Ep) + 2 * (sp.cos(2 * k * z_) + sp.cos(2 * k * t))) == 0
ok &= sp.simplify(Ec.dot(Ec) + Bc.dot(Bc) - Ep.dot(Ep) - Bp.dot(Bp)) == 0 and sp.simplify(Ec.cross(Bc) - Ep.cross(Bp)) == sp.zeros(3, 1)
Rz = sp.Matrix([[1, -1, 0], [1, 1, 0], [0, 0, sp.sqrt(2)]]) / sp.sqrt(2)
al = -sp.pi / 4
Er, Br = Rz * Ep, Rz * Bp
ok &= sp.simplify(sp.cos(al) * Er + sp.sin(al) * Br - Ec) == sp.zeros(3, 1) and sp.simplify(sp.cos(al) * Br - sp.sin(al) * Er - Bc) == sp.zeros(3, 1)
# prop:EBlocal, frame variables and potentials
jj = j_
sub2 = {sp.Derivative(Bf[i], t): curl(Ef)[i] for i in range(3)}
sub2.update({sp.Derivative(Ef[i], t): -curl(Bf)[i] - jj[i] for i in range(3)})
q = z.dot(z)
ok &= sp.expand(sp.diff(q, t).subs(sub2) - (-2 * sp.I * z.dot(curl(z)) - 2 * sp.I * z.dot(jj))) == 0
ok &= sp.expand(sp.diff(Ef.dot(Bf), t).subs(sub2) - (Ef.dot(curl(Ef)) - Bf.dot(curl(Bf)) - jj.dot(Bf))) == 0
Bt = -Bf
ok &= sp.expand(sp.diff(Ef.dot(Bt), t).subs(sub2) - (Bt.dot(curl(Bt)) - Ef.dot(curl(Ef)) - jj.dot(Bt))) == 0
Av = sp.Matrix([sp.Function(f"A{i}")(*co) for i in (1, 2, 3)])
phi = sp.Function("phi")(*co)
Cv = sp.Matrix([sp.Function(f"C{i}")(*co) for i in (1, 2, 3)])
chi = sp.Function("chi")(*co)
grad = lambda f: sp.Matrix([sp.diff(f, xx) for xx in X])  # noqa: E731
BtA, EA = curl(Av), -sp.Matrix([sp.diff(c, t) for c in Av]) - grad(phi)
ok &= sp.expand(sp.diff(Av.dot(BtA), t) + div(phi * BtA + EA.cross(Av)) + 2 * EA.dot(BtA)) == 0
EC, BtC = -curl(Cv), -sp.Matrix([sp.diff(c, t) for c in Cv]) + grad(chi)
ok &= sp.expand(sp.diff(Cv.dot(EC), t) - div(BtC.cross(Cv) + chi * EC) + 2 * EC.dot(BtC)) == 0
report("5", ok, "ex:circular, ex:threewaves, prop:orthframe, ex:counter, prop:EBlocal: the circular wave solves the frame system with curl z = k z, q = 0 and "
       "Poynting (i/2) zbar x z = E x B = (0,0,-1), E x Btext = (0,0,1); z_+, z_-, z_lin solve, are null with B x E along +e3 and (E x Edot).e3 = "
       "+k|E|^2, -k|E|^2, 0 at x3 = 0; the mixed field solves with |B| = |E|, B.E = cos cos; z.z = eta(X,X) - i eta(X,GX); both counter-propagating "
       "configurations solve Maxwell with the stated invariants, equal energy density and Poynting vector, and crossed = duality(-pi/4) o R_z(pi/4) "
       "(parallel) in the convention E' = cos a E + sin a B; the local laws of q and E.B, and the two potential identities hold identically (exact)")

# ---------------------------------------------------------------- a2-6: prop:planes, prop:modes, prop:twomodes


def pf(F):
    return F[0, 1] * F[2, 3] - F[0, 2] * F[1, 3] + F[0, 3] * F[1, 2]


Fs = Fmat(Bs, Es)
ok = sp.expand(pf(Fs) - Es.dot(Bs)) == 0 and sp.expand(Fs.det() - Es.dot(Bs) ** 2) == 0
a_, b_ = sp.symbols("a b", positive=True)
Xab = sp.Matrix([b_, 0, 0, a_, 0, 0])
th = sp.atan(a_ / b_)
Xr = sp.simplify((sp.cos(th) * sp.eye(6) - sp.sin(th) * G6) * Xab)
ok &= sp.simplify(Xr - sp.Matrix([sp.sqrt(a_ ** 2 + b_ ** 2), 0, 0, 0, 0, 0])) == sp.zeros(6, 1)
ok &= Fmat(Xab[:3, 0], Xab[3:, 0]) == sp.Matrix([[0, a_, 0, 0], [-a_, 0, 0, 0], [0, 0, 0, b_], [0, 0, -b_, 0]])
wedge = lambda u, v: u * v.T - v * u.T  # noqa: E731
U = [sp.Matrix(sp.symbols(f"{n}0:4")) for n in "kbuw"]
ok &= sp.expand(pf(wedge(U[0], U[1]))) == 0
ok &= sp.expand(pf(wedge(U[0], U[1]) + wedge(U[2], U[3])) - sp.Matrix.hstack(*U).det()) == 0
kk, bb = U[0], U[1]
dot = lambda u, v: (u.T * eta4 * v)[0]  # noqa: E731
Fk = wedge(kk, bb)  # lower-index two-form
inv = sp.expand(sum(Fk[m, n] * (eta4 * Fk * eta4)[m, n] for m in range(4) for n in range(4)) / 2)
ok &= sp.expand(inv - (dot(kk, kk) * dot(bb, bb) - dot(kk, bb) ** 2)) == 0
# the frame invariant Re q = |B|^2 - |E|^2 equals F_mn F^mn / 2 for F_X
FXs = Fmat(Bs, Es)
ok &= sp.expand(sum(FXs[m, n] * (eta4 * FXs * eta4)[m, n] for m in range(4) for n in range(4)) / 2 - (Bs.dot(Bs) - Es.dot(Es))) == 0
# Proca: k timelike, b orthogonal to k is spacelike, so (k.k)(b.b) <= 0
mm = sp.symbols("m", positive=True)
bsp = sp.Matrix([0] + list(sp.symbols("y1:4")))
ok &= sp.expand(mm ** 2 * dot(bsp, bsp) + mm ** 2 * sum(y ** 2 for y in bsp[1:])) == 0
lv = sp.Matrix([1, 0, 0, 1])
om_ = sp.symbols("w1:3")
b1, b2 = sp.Matrix(sp.symbols("p0:4")), sp.Matrix(sp.symbols("s0:4"))
Ftot = om_[0] * wedge(lv, b1) + om_[1] * wedge(lv, b2)
ok &= sp.expand(Ftot - wedge(lv, om_[0] * b1 + om_[1] * b2)) == sp.zeros(4) and sp.expand(pf(Ftot)) == 0
report("6", ok, "prop:planes, prop:modes, prop:twomodes: Pf F_X = E.B and det F_X = (E.B)^2; e^{-theta G} with tan theta = a/b carries a th01 + b th23 to "
       "sqrt(a^2+b^2) th23; Pf(k^b) = 0 and Pf(k1^b1 + k2^b2) = det(k1,b1,k2,b2); for F = k^b, F.F/2 = (k.k)(b.b) - (k.b)^2, equal to |B|^2-|E|^2; "
       "a Proca mode in the rest frame has q = m^2 b.b <= 0; co-propagating waves give one null plane (exact)")

# ---------------------------------------------------------------- a2-7: prop:gradeone, prop:samerotors, prop:censusint, prop:spinmatrices, prop:halfspinors
Sk = [gam[2] * gam[3] / 2, gam[3] * gam[1] / 2, gam[1] * gam[2] / 2]
Tk = [gam[0] * gam[a] / 2 for a in (1, 2, 3)]
J = e[2] * e[3]
ok = all(comm(Sk[i], Sk[j]) == sum((eps(i + 1, j + 1, kq + 1) * Sk[kq] for kq in range(3)), Z8) for i in range(3) for j in range(3))
adS = [ad6(s) for s in Sk]
adT = [ad6(s) for s in Tk]
Cas1 = sum((x * x for x in adT), sp.zeros(6)) - sum((x * x for x in adS), sp.zeros(6))
ok &= Cas1 == sp.diag(3, 0, 0, 3, 3, 3)
adJ = ad6(J / 2)
ok &= adJ[1:3, 1:3] == sp.Matrix([[0, 1], [-1, 0]]) and all(adJ[i, j] == 0 for i in range(6) for j in range(6) if not (i in (1, 2) and j in (1, 2)))
ok &= -sum((x * x for x in adS), sp.zeros(6)) == sp.diag(0, 0, 0, 2, 2, 2)
ok &= all(adS[kq][3:, 3:] == rk[kq + 1] for kq in range(3)) and all(adS[kq][:3, :] == sp.zeros(3, 6) for kq in range(3))
ok &= adS[2].eigenvals() == {sp.I: 1, -sp.I: 1, 0: 4}
ok &= -sum((s * s for s in Sk), Z8) == sp.Rational(3, 4) * I8 and Sk[2].eigenvals() == {sp.I / 2: 4, -sp.I / 2: 4}
ok &= sp.simplify((2 * sp.pi * adS[2]).exp() - sp.eye(6)) == sp.zeros(6) and sp.simplify((2 * sp.pi * Sk[2]).exp() + I8) == Z8
ok &= sum((x * x for x in Tk), Z8) - sum((x * x for x in Sk), Z8) == sp.Rational(3, 2) * I8
# census on all grades, by the adjoint action on the blades of grade r
census = {}
for r in range(7):
    blades = [bl(*S) for S in itertools.combinations(range(1, 7), r)] if r else [I8]

    def coords(Y, blades=blades):
        return [sp.nsimplify((b.inv() * Y).trace() / 8) for b in blades]

    def adr(Xb, blades=blades, coords=coords):
        return sp.Matrix([coords(comm(Xb, b)) for b in blades]).T

    C = sum((adr(x) ** 2 for x in Tk), sp.zeros(len(blades))) - sum((adr(x) ** 2 for x in Sk), sp.zeros(len(blades)))
    spec = {sp.nsimplify(kq): v for kq, v in C.eigenvals().items()}
    four = (C - 4 * sp.eye(len(blades))).nullspace()
    charges = set()
    if four:
        Pm = sp.Matrix.hstack(*four)
        AJ = adr(J / 2)
        red = (Pm.T * Pm).inv() * Pm.T * AJ * Pm
        charges = {sp.nsimplify(-sp.I * ev) for ev in red.eigenvals()}
    census[r] = (spec, len(four) // 6, charges)
ok &= census[0][0] == {0: 1} and census[6][0] == {0: 1} and census[1][0] == {3: 4, 0: 2} and census[5][0] == {3: 4, 0: 2}
ok &= census[2][0] == {4: 6, 3: 8, 0: 1} and census[4][0] == {4: 6, 3: 8, 0: 1} and census[3][0] == {4: 12, 3: 8}
ok &= census[2][1] == 1 and census[4][1] == 1 and census[3][1] == 2 and census[2][2] == {0} and census[3][2] == {1, -1}
# prop:spinmatrices on V
ok &= all(K6[kq] == G6 * L6[kq] for kq in (1, 2, 3)) and sum((K6[kq] ** 2 for kq in (1, 2, 3)), sp.zeros(6)) - sum((L6[kq] ** 2 for kq in (1, 2, 3)), sp.zeros(6)) == 4 * sp.eye(6)
s_ = {kq: sp.I * rk[kq] for kq in (1, 2, 3)}
ok &= all(s_[kq] == sp.Matrix(3, 3, lambda a, c: -sp.I * eps(kq, a + 1, c + 1)) and s_[kq].H == s_[kq] for kq in (1, 2, 3))
ok &= all(comm(s_[i], s_[j]) == sum((sp.I * eps(i, j, kq) * s_[kq] for kq in (1, 2, 3)), sp.zeros(3)) for i in (1, 2, 3) for j in (1, 2, 3))
ok &= sum((s_[kq] ** 2 for kq in (1, 2, 3)), sp.zeros(3)) == 2 * sp.eye(3) and s_[3].eigenvals() == {-1: 1, 0: 1, 1: 1}
zf = sp.Matrix([sp.Function(f"z{i}")(*co) for i in (1, 2, 3)])
sp_z = sum((s_[kq] * (-sp.I) * sp.Matrix([sp.diff(c, X[kq - 1]) for c in zf]) for kq in (1, 2, 3)), sp.zeros(3, 1))
ok &= sp.expand(sp_z - curl(zf)) == sp.zeros(3, 1) and s_[3].nullspace() == [sp.Matrix([0, 0, 1])]
# action of L_k, K_k on z = B + iE is r_k and i r_k
ok &= all(sp.expand((L6[kq] * X6)[:3, 0] + sp.I * (L6[kq] * X6)[3:, 0] - rk[kq] * zg) == sp.zeros(3, 1) for kq in (1, 2, 3))
ok &= all(sp.expand((K6[kq] * X6)[:3, 0] + sp.I * (K6[kq] * X6)[3:, 0] - s_[kq] * zg) == sp.zeros(3, 1) for kq in (1, 2, 3))
# half-spinors: S^pm = +-1 eigenspaces of the volume element e1...e6 (square +1)
Vol = bl(1, 2, 3, 4, 5, 6)
ok &= Vol * Vol == I8
Sp_ = sp.Matrix.hstack(*(Vol - I8).nullspace())
Sm_ = sp.Matrix.hstack(*(Vol + I8).nullspace())


def restrict(Xb, P):
    return (P.T * P).inv() * P.T * Xb * P


lint_p = [restrict(x, Sp_) for x in LINT8]
lG_p = [restrict(x, Sp_) for x in list(Lb.values()) + list(Kb.values())]
CasI = [sum((restrict(x, P) ** 2 for x in Tk), sp.zeros(4)) - sum((restrict(x, P) ** 2 for x in Sk), sp.zeros(4)) for P in (Sp_, Sm_)]
ok &= all(c == sp.Rational(3, 2) * sp.eye(4) for c in CasI)
CasG = sum((Kb[kq] ** 2 for kq in (1, 2, 3)), Z8) - sum((Lb[kq] ** 2 for kq in (1, 2, 3)), Z8)
ok &= CasG == 3 * I8
# invariant forms on S^+ with respect to the Euclidean basis of S^+ obtained from the nullspace (forms: X^T Q + Q X = 0)
Pp = sp.Matrix.hstack(*sp.GramSchmidt([Sp_[:, i] for i in range(4)], True))
lint_p = [Pp.T * x * Pp for x in LINT8]
lG_p = [Pp.T * x * Pp for x in list(Lb.values()) + list(Kb.values())]
anti_I = solve_mats([(lambda x: (lambda Q: list(x.T * Q + Q * x)))(x) for x in lint_p] + [lambda Q: list(Q + Q.T)], 4)
sym_I = solve_mats([(lambda x: (lambda Q: list(x.T * Q + Q * x)))(x) for x in lint_p] + [lambda Q: list(Q - Q.T)], 4)
anti_G = solve_mats([(lambda x: (lambda Q: list(x.T * Q + Q * x)))(x) for x in lG_p] + [lambda Q: list(Q + Q.T)], 4)
sym_G = solve_mats([(lambda x: (lambda Q: list(x.T * Q + Q * x)))(x) for x in lG_p] + [lambda Q: list(Q - Q.T)], 4)
ok &= len(anti_I) == 2 and len(sym_I) == 0 and len(anti_G) == 0 and len(sym_G) == 1
sg = sym_G[0].eigenvals()
npos = sum(m for v, m in sg.items() if v > 0)
ok &= sorted([npos, 4 - npos]) == [1, 3]
ok &= sp.simplify((2 * sp.pi * Lb[3]).exp() - I8) == Z8
report("7", ok, "prop:gradeone, prop:samerotors, prop:censusint, prop:spinmatrices, prop:halfspinors: Casimir sum ad(T)^2 - sum ad(S)^2 = diag(3,0,0,3,3,3) "
       "on grade one and 3/2 on S; ad(J/2) as stated; [S_i,S_j] = eps S_k, -sum (ad S)^2 = diag(0,0,0,2,2,2), spec ad S3 = {+-i,0^4}, S3 = +-i/2, "
       "exp(2 pi ad S3) = 1, exp(2 pi S3) = -1; the full grade census (Casimir spectra, field-strength modules 1,2 with charges 0, +-1); K_k = G L_k, "
       "sum K^2 - sum L^2 = 4 on V, s_k = i r_k Hermitian with [s_i,s_j] = i eps s_k, sum s^2 = 2, spec s3 = {-1,0,1}, (s.p)z = curl z, ker s3 = e3; "
       "on S^+ (the +1 eigenspace of e1...e6): l_int has Casimir 3/2, 2 invariant antisymmetric and 0 symmetric forms; l_G has Casimir 3 on S, "
       "1 invariant symmetric form of signature (1,3) and no antisymmetric one; exp(2 pi L3) = +1 on S (exact)")

# ---------------------------------------------------------------- a2-8: thm:nablaF, ex:univectorwave
etaU = [1, -1, -1, -1]
gl = [etaU[m] * gam[m] for m in range(4)]  # gamma_mu
mons = {S: reduce(lambda a, b: a * b, [gam[i] for i in S], I8) for r in range(5) for S in itertools.combinations(range(4), r)}


def grade(Y, r):
    return sum(((mons[S].inv() * Y).trace() / 8 * mons[S] for S in mons if len(S) == r), Z8)


nabla = lambda Y: sum((gam[m] * Y.applyfunc(lambda f: sp.diff(f, co[m])) for m in range(4)), Z8)  # noqa: E731
Au = [sp.Function(f"a{m}")(*co) for m in range(4)]
Avec = sum((Au[m] * gl[m] for m in range(4)), Z8)
NA = nabla(Avec)
divA = sum(sp.diff(Au[m], co[m]) for m in range(4))
ok = sp.expand(grade(NA, 0) - divA * I8) == Z8 and sp.expand(NA - grade(NA, 0) - grade(NA, 2)) == Z8
Fa = grade(NA, 2)
NF = nabla(Fa)
box = lambda f: sp.diff(f, t, 2) - sum(sp.diff(f, xx, 2) for xx in X)  # noqa: E731
up = lambda m, f: etaU[m] * f  # noqa: E731
g1 = sum(((box(Au[n]) - up(n, sp.diff(divA, co[n]))) * gl[n] for n in range(4)), Z8)
ok &= sp.expand(grade(NF, 1) - g1) == Z8 and sp.expand(grade(NF, 3)) == Z8
sig = [gl[kq] * gl[0] for kq in (1, 2, 3)]
Ips = gl[0] * gl[1] * gl[2] * gl[3]
ok &= Ips == -Om
Ek = [-sp.diff(Au[0], X[kq]) - sp.diff(Au[kq + 1], t) for kq in range(3)]
Bk = list(curl(sp.Matrix(Au[1:])))
ok &= sp.expand(Fa - sum((Ek[kq] * sig[kq] + Bk[kq] * Ips * sig[kq] for kq in range(3)), Z8)) == Z8
Fgen = sum((Ef[kq] * sig[kq] + Bf[kq] * Ips * sig[kq] for kq in range(3)), Z8)
NFg = nabla(Fgen)
cB, cE = curl(Bf), curl(Ef)
rhs = (div(Ef) * gl[0] + sum(((cB[kq] - sp.diff(Ef[kq], t)) * gl[kq + 1] for kq in range(3)), Z8)
       + Ips * (-div(Bf) * gl[0] + sum(((sp.diff(Bf[kq], t) + cE[kq]) * gl[kq + 1] for kq in range(3)), Z8)))
ok &= sp.expand(NFg - rhs) == Z8
kv = gl[0] + gl[3]
ok &= kv * kv == Z8
eps_ = sp.symbols("c0:4")
ev = eps_[0] * gl[0] + eps_[1] * gl[1] + eps_[2] * gl[2] + eps_[0] * gl[3]
Fj = [kv * gl[1], kv * gl[2]]
ok &= sp.expand(kv * ev - eps_[1] * Fj[0] - eps_[2] * Fj[1]) == Z8
ok &= all(f * f == Z8 for f in Fj) and Fj[0] * Fj[1] + Fj[1] * Fj[0] == Z8
S3u = gam[1] * gam[2] / 2


def on_span(Op):
    cols = []
    for f in Fj:
        y = Op(f)
        a1, a2 = sp.symbols("u v")
        sol = sp.solve(list(y - a1 * Fj[0] - a2 * Fj[1]), [a1, a2], dict=True)
        if not sol:
            return None
        cols.append([sol[0][a1], sol[0][a2]])
    return sp.Matrix(cols).T


adS3 = on_span(lambda f: comm(S3u, f))
OmF = on_span(lambda f: Ips * f)
ok &= adS3 == sp.Matrix([[0, -1], [1, 0]]) and OmF == sp.Matrix([[0, 1], [-1, 0]])
# F^2 = E^2 - B^2 + 2 Omega E.B in the space-time algebra
Fsq = sp.expand(Fgen * Fgen)
ok &= sp.expand(Fsq - ((Ef.dot(Ef) - Bf.dot(Bf)) * I8 + 2 * Ips * Ef.dot(Bf))) == Z8
report("8", ok, "thm:nablaF (i)-(iii) and ex:univectorwave (as corrected: pseudoscalar I_4 = gamma_0gamma_1gamma_2gamma_3 = -Omega): nabla A = d.A + nabla^A; grade one of nabla F is (box A^nu - d^nu d.A) gamma_nu and grade three "
       "vanishes; F = E + I B with E_k = -d_k A^0 - d_t A^k, B = curl A; the formula for nabla F with no other grades (general E, B); k^2 = 0, "
       "k eps = eps^1 F_1 + eps^2 F_2, F_j^2 = 0, F_1F_2 + F_2F_1 = 0; ad S3 and I on span{F1,F2} are the stated matrices (Omega = +ad S3 there); "
       "F^2 = E^2 - B^2 + 2 I E.B (exact)")

# ---------------------------------------------------------------- a2-9: thm:SS, ex:Gcharge
XA = sp.Matrix([sp.Function(f"X{i}")(*co) for i in range(6)])
curl6 = lambda v: sp.Matrix.vstack(curl(v[:3, 0]), curl(v[3:, 0]))  # noqa: E731
XAd = sp.Matrix([sp.diff(c, t) for c in XA])
Lag = (curl6(XA).T * G6 * XAd)[0] / 2 - (curl6(XA).T * curl6(XA))[0] / 2
from sympy.calculus.euler import euler_equations  # noqa: E402
EL = euler_equations(Lag, list(XA), list(co))
ELv = sp.Matrix([eq.lhs for eq in EL])
target = G6 * curl6(XAd + G6 * curl6(XA))
sgnEL = [s for s in (1, -1) if sp.expand(ELv - s * target) == sp.zeros(6, 1)]
ok = len(sgnEL) == 1
ok &= sp.expand(curl6(-G6 * curl6(XA)) + G6 * curl6(curl6(XA))) == sp.zeros(6, 1)  # curl commutes with G, so X_F obeys dt X_F = -G curl X_F
l1, l2 = sp.Function("l1")(*co), sp.Function("l2")(*co)
gvec = sp.Matrix.vstack(grad(l1), grad(l2))
ok &= sp.expand(curl6(gvec)) == sp.zeros(6, 1)
dL = sp.expand(Lag.subs({XA[i]: XA[i] + gvec[i] for i in range(6)}).doit() - Lag)
# dL is a divergence: it equals d_i( (1/2) (curl X)_a (G d_t grad l)_a ... ); check by the Euler operator, which annihilates divergences
ELd = euler_equations(dL, [l1, l2], list(co)) if dL != 0 else []
ok &= all(sp.expand(eq.lhs) == 0 for eq in ELd)
pmom = sp.Matrix([sp.diff(Lag, XAd[i]) for i in range(6)])
Ham = sp.expand((pmom.T * XAd)[0] - Lag)
ok &= sp.expand(Ham - (curl6(XA).T * curl6(XA))[0] / 2) == 0
thv = sp.symbols("theta", real=True)
R = sp.cos(thv) * sp.eye(6) + sp.sin(thv) * G6
ok &= sp.simplify((R * G6 - G6 * R)) == sp.zeros(6) and sp.simplify(R.T * R) == sp.eye(6)
LR = (curl6(R * XA).T * G6 * (R * XAd))[0] / 2 - (curl6(R * XA).T * curl6(R * XA))[0] / 2
ok &= sp.simplify(sp.expand(LR - Lag)) == 0
j0 = sp.expand((pmom.T * G6 * XA)[0])
ok &= sp.expand(j0 + (XA.T * curl6(XA))[0] / 2) == 0
kz = sp.symbols("k", positive=True)
kx = lambda a: sp.Matrix.vstack(sp.Matrix([0, 0, kz]).cross(a[:3, 0]), sp.Matrix([0, 0, kz]).cross(a[3:, 0]))  # noqa: E731
Mw = sp.Matrix.hstack(*[G6 * kx(sp.eye(6)[:, i]) for i in range(6)])
ok &= Mw.eigenvals() == {kz: 2, -kz: 2, 0: 2}
v1, v2 = sp.Matrix([1, sp.I, 0]), sp.Matrix([1, -sp.I, 0])
a1 = sp.Matrix.vstack(v1, -sp.I * v1)
a2 = sp.Matrix.vstack(v2, sp.I * v2)
s3_6 = sp.diag(s_[3], s_[3])
ok &= Mw * a1 == kz * a1 and Mw * a2 == kz * a2 and G6 * a1 == sp.I * a1 and G6 * a2 == -sp.I * a2
ok &= s3_6 * a1 == a1 and s3_6 * a2 == -a2
out = []
for a in (a1, a2, a1 + a2):
    Xr_ = sp.Matrix([sp.re(sp.expand_complex(c * sp.exp(sp.I * kz * (x3 - t)))) for c in a])
    out.append((sp.simplify((Xr_.T * curl6(Xr_))[0] / 2), sp.simplify((curl6(Xr_).T * curl6(Xr_))[0] / 2)))
ok &= sp.simplify(out[0][0] - kz) == 0 and sp.simplify(out[0][1] - kz ** 2) == 0
ok &= sp.simplify(out[1][0] + kz) == 0 and sp.simplify(out[1][1] - kz ** 2) == 0
ok &= sp.simplify(out[2][0]) == 0 and sp.simplify(out[2][1] - 4 * kz ** 2 * sp.sin(kz * (x3 - t)) ** 2) == 0
report("9", ok, f"thm:SS, ex:Gcharge: the Euler-Lagrange expression of eq:SSL is {'+' if sgnEL == [1] else '-'}G curl(dX_A + G curl X_A); curl commutes with G, so "
       "X_F obeys the frame system; gauge shifts leave X_F unchanged and change L by a term annihilated by the Euler operator (a divergence); "
       "H = |curl X_A|^2/2; L(e^{theta G} X_A) = L(X_A) and the Noether density of GX_A is -X_A.curl X_A/2; omega a = G(k x a) has eigenvalues "
       "+-k (twice) and 0 (twice); a1, a2 as stated with G = +-i, s3 = +-1; helicity and energy densities (k,k^2), (-k,k^2), (0, 4k^2 sin^2) (exact)")

# ---------------------------------------------------------------- a2-10: thm:klein
BIV = [e[a] * e[b] for a, b in itertools.combinations(range(1, 7), 2)]
forms = solve_mats([(lambda x: (lambda Q: list(x.T * Q + Q * x)))(x) for x in BIV], 8)
ok = len(forms) == 2
xv = sp.symbols("y1:7")
Xv = sum((xv[a - 1] * e[a] for a in range(1, 7)), Z8)
etaX = sum(xv[a] ** 2 for a in range(3)) - sum(xv[a] ** 2 for a in range(3, 6))
found = []
for beta in forms:
    W = Pp.T * beta * Xv * Pp
    if sp.expand(W + W.T) == sp.zeros(4):
        p_ = sp.expand(pf(W))
        ratio = sp.simplify(p_ / etaX)
        if ratio.is_number and ratio != 0:
            found.append(ratio)
ok &= len(found) >= 1
# equivariance of W: W(Ad_g X) = rho(g)^T-twisted; infinitesimally, for every bivector b: beta(psi,[b,X]chi) = -(beta(b psi, X chi) + beta(psi, X b chi)) on S^+
beta = [f for f in forms if sp.expand((Pp.T * f * Xv * Pp) + (Pp.T * f * Xv * Pp).T) == sp.zeros(4)][0]
ok &= all(sp.expand(Pp.T * (beta * comm(b, Xv)) * Pp + Pp.T * (b.T * beta * Xv + beta * Xv * b) * Pp) == sp.zeros(4) for b in BIV)
report("10", ok, f"thm:klein: the Spin(3,3)-invariant bilinear forms on S form a 2-dimensional space; for one of them W(X)(psi,chi) = beta(psi, X chi) is "
       f"antisymmetric on S^+ and equivariant, with Pf W(X) = c eta(X,X), c = {found[0] if found else '?'} in an orthonormal basis of S^+ (c > 0 can be "
       "normalised to 1 by scaling beta; the sign of Pf depends on the orientation of the basis); so W is injective off the cone and Pf W = 0 exactly "
       "on the null cone (exact)")

# ---------------------------------------------------------------- a2-11: en:SI dimensions and the split identity
zz = sp.Matrix(sp.symbols("r1:4", real=True)) + sp.I * sp.Matrix(sp.symbols("i1:4", real=True))
zb = sp.Matrix([sp.conjugate(c) for c in zz])
n2 = sp.expand(zz.dot(zb))
cr = zb.cross(zz)
lhs = sp.expand(n2 ** 2)
rhs = sp.expand(sp.expand(cr.dot(sp.Matrix([sp.conjugate(c) for c in cr]))) + sp.expand(zz.dot(zz) * sp.conjugate(zz.dot(zz))))
ok = sp.simplify(lhs - rhs) == 0
from sympy.physics.units import farad, meter, volt, joule, kilogram, second, tesla, henry, convert_to  # noqa: E402
eps0, mu0 = farad / meter, henry / meter
Eu, Bu, cu = volt / meter, tesla, meter / second
dims = [convert_to(eps0 * Eu ** 2, [joule, meter, kilogram, second]), convert_to(eps0 * (cu * Bu) ** 2, [joule, meter, kilogram, second]),
        convert_to(cu * (Eu * Bu / mu0) / cu ** 2, [joule, meter, kilogram, second]), convert_to(eps0 * Eu ** 2 / cu ** 2, [kilogram, meter])]
dims = [convert_to(x, [joule, meter, second]) for x in (eps0 * Eu ** 2, eps0 * (cu * Bu) ** 2, cu * (Eu * Bu / mu0) / cu ** 2)] + \
    [convert_to(eps0 * Eu ** 2 / cu ** 2, [kilogram, meter, second])]
ok &= all(d == joule / meter ** 3 for d in dims[:3]) and dims[3] == kilogram / meter ** 3
report("11", ok, "en:SI: |z|^4 = |zbar x z|^2 + |z.z|^2 for complex 3-vectors (the identity behind u^2 = c^2|g|^2 + u0^2); [eps0 E^2] = [eps0 c^2 B^2] = "
       "[c g] = J m^-3 and [u0/c^2] = kg m^-3 (exact, sympy units)")

print(f"{NPASS[0]}/{NPASS[0] + len(FAILS)} checks passed")
sys.exit(0 if not FAILS else 1)
