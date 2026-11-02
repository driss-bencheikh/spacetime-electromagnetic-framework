"""Phase 9D, batch 1: adversarial audit of the Chapter 1-2 statements that had neither a cited check nor a written proof.

Exact arithmetic throughout (sympy rationals/integers on the 8x8 realisation of Appendix A; integer 16x16 matrices for the module).
Each check prints PASS/FAIL a1-<n>; the script exits 0 iff all pass.
"""
import itertools
import os
import sys
from functools import reduce

import sympy as sp

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "common"))
from clifford8 import I8, bl, comm, e  # noqa: E402

NPASS, FAILS = [0], []


def report(cid, ok, text):
    if ok:
        NPASS[0] += 1
        print(f"PASS a1-{cid}: {text}")
    else:
        FAILS.append(cid)
        print(f"FAIL a1-{cid}: {text}")
    sys.stdout.flush()


Z8 = sp.zeros(8, 8)
BLADES = [S for k in range(7) for S in itertools.combinations(range(1, 7), k)]
BASIS = [bl(*S) if S else I8 for S in BLADES]
g0, ga = e[1], [e[4], e[5], e[6]]
gam = [g0] + ga
J = e[2] * e[3]
T = e[1] * e[2] * e[3]
Om = e[1] * e[4] * e[5] * e[6]
Gam = J * Om
LOR = [gam[m] * gam[n] for m, n in itertools.combinations(range(4), 2)]


def solve_space(conds):
    """Real solution space of linear conditions L(X)=0 on M8(R), as a list of matrices (exact nullspace in the blade basis)."""
    rows = []
    for L in conds:
        cols = [list(L(B)) for B in BASIS]
        rows.append(sp.Matrix(cols).T)
    M = sp.Matrix.vstack(*rows)
    return [sum((v[i] * BASIS[i] for i in range(64)), Z8) for v in M.nullspace()]


def span_eq(A, B):
    """span(A) == span(B) for lists of 8x8 matrices."""
    fa = sp.Matrix([list(x) for x in A]) if A else sp.zeros(0, 64)
    fb = sp.Matrix([list(x) for x in B]) if B else sp.zeros(0, 64)
    ra, rb = fa.rank(), fb.rank()
    return ra == rb == sp.Matrix.vstack(fa, fb).rank()


# ---------------------------------------------------------------- Chapter 1
biv = {(a, b): e[a] * e[b] for a, b in itertools.combinations(range(1, 7), 2)}
eta = [None, 1, 1, 1, -1, -1, -1]
ok = all(B * B == -eta[a] * eta[b] * I8 for (a, b), B in biv.items())
same = [p for p in biv if eta[p[0]] == eta[p[1]]]
mixed = [p for p in biv if eta[p[0]] != eta[p[1]]]
ok &= len(same) == 6 and len(mixed) == 9
ok &= all((biv[p].T == -biv[p]) == (p in same) and (biv[p].T == biv[p]) == (p in mixed) for p in biv)
tt, ss = [p for p in same if p[0] <= 3], [p for p in same if p[0] >= 4]
spanT, spanS = [biv[p] for p in tt], [biv[p] for p in ss]
ok &= all(span_eq(spanT, spanT + [comm(biv[p], biv[q])]) for p in tt for q in tt)
ok &= all(span_eq(spanS, spanS + [comm(biv[p], biv[q])]) for p in ss for q in ss)
ok &= all(comm(biv[p], biv[q]) == Z8 for p in tt for q in ss)
th = sp.symbols("theta", real=True)
ok &= all(sp.simplify((th * biv[p]).exp() - (sp.cosh(th) * I8 + sp.sinh(th) * biv[p])) == Z8 for p in [(1, 4), (3, 6)])
ok &= sp.simplify((th * biv[(2, 3)]).exp() - (sp.cos(th) * I8 + sp.sin(th) * biv[(2, 3)])) == Z8
ok &= comm(biv[(2, 3)], biv[(1, 4)]) == Z8 and comm(biv[(2, 3)], biv[(5, 6)]) == Z8
ok &= biv[(1, 4)] * biv[(1, 2)] == -biv[(1, 2)] * biv[(1, 4)]
report("1", ok, "thm:sig, cor:noboost, ex:blades: (e_a e_b)^2 = -e_a^2 e_b^2 for all 15 bivectors; the 6 same-sign ones are the antisymmetric ones and "
       "close into two commuting so(3) triples (time-time, space-space); the 9 mixed ones are symmetric with exp(tB) = cosh t + sinh t B; "
       "e2e3 rotation of period 2pi; e2e3 commutes with e1e4 and e5e6; e1e4 anticommutes with e1e2 (exact)")

om = [e[1] * e[4], e[2] * e[5], e[3] * e[6]]
vp = om[0] + om[1] + om[2]
V = bl(1, 2, 3, 4, 5, 6)
ok = all(o * o == I8 for o in om) and all(comm(a, b) == Z8 for a in om for b in om)
ok &= V == -om[0] * om[1] * om[2]
ev = vp.eigenvals()
ok &= ev == {3: 1, 1: 3, -1: 3, -3: 1}
sgn = (vp + V) / 2
ok &= sgn == (13 * vp - vp ** 3) / 12
Pi = I8 / 2 - (vp + V) / 4
ok &= Pi * Pi == Pi and Pi.T == Pi and Pi.rank() == 4
# the range of Pi is the sum of the strata varpi = -1, -3
ok &= all((vp * v).T * v != sp.zeros(1, 1) and ((vp * v).T * v)[0] < 0 for v in Pi.columnspace())
report("2", ok, "ex:pairedmode: pole operators commute and square to 1; V = -w1w2w3; spec varpi = {3,1,-1,-3} with multiplicities 1,3,3,1; "
       "(varpi+V)/2 = (13 varpi - varpi^3)/12 = sign(varpi); Pi = 1/2 - (varpi+V)/4 is the orthogonal projector of rank 4 onto varpi < 0 (exact)")

# ---------------------------------------------------------------- Chapter 2, clock section
c1, c2, c3 = sp.symbols("c1 c2 c3", positive=True)
C = sp.diag(c1, c2, c3)
KC = sp.BlockMatrix([[sp.zeros(3), C.inv()], [C, sp.zeros(3)]]).as_explicit()
ok = sp.simplify(KC * KC - sp.eye(6)) == sp.zeros(6)
L = [sp.Matrix([[0, 0, 0], [0, 0, -1], [0, 1, 0]]), sp.Matrix([[0, 0, 1], [0, 0, 0], [-1, 0, 0]]), sp.Matrix([[0, -1, 0], [1, 0, 0], [0, 0, 0]])]
cs = sp.symbols("x0:9")
X3 = sp.Matrix(3, 3, cs)
sol = sp.solve([x for Lk in L for x in (X3 * Lk - Lk * X3)], cs, dict=True)[0]
X3s = X3.subs(sol)
ok &= len(X3s.free_symbols) == 1 and X3s == X3s[0, 0] * sp.eye(3)
report("3", ok, "prop:KC, thm:onec: K_C^2 = 1 for symbolic diag(c1,c2,c3) with distinct entries; the commutant of the vector representation of so(3) "
       "in M3(R) is R*id (the Schur step of thm:onec, exact)")

I3, Z3 = sp.eye(3), sp.zeros(3)
eta6 = sp.diag(1, 1, 1, -1, -1, -1)
K = sp.BlockMatrix([[Z3, -I3], [-I3, Z3]]).as_explicit()
G = sp.BlockMatrix([[Z3, -I3], [I3, Z3]]).as_explicit()
dag = lambda X: eta6.inv() * X.T * eta6  # noqa: E731
E6 = sp.eye(6)
ok = K * K == E6 and G * G == -E6 and dag(K) == -K and dag(G) == G
ok &= K.T * eta6 * K == -eta6 and G.T * eta6 * G == -eta6
ok &= all(K * E6[:, i] == -E6[:, i + 3] and K * E6[:, i + 3] == -E6[:, i] and G * E6[:, i] == E6[:, i + 3] and G * E6[:, i + 3] == -E6[:, i] for i in range(3))
for s in (1, -1):
    vs = [E6[:, i] - s * E6[:, i + 3] for i in range(3)]
    ok &= all(K * v == s * v for v in vs) and all((u.T * eta6 * v)[0] == 0 for u in vs for v in vs)
report("4", ok, "prop:twoK: K^2 = 1, G^2 = -1, K^dag = -K (K in so(3,3)), G^dag = G (outside so(3,3)), both anti-isometries (X^T eta X = -eta), the "
       "stated actions on e_i, e_{i+3}; the eigenspaces of K are the totally null three-planes e_i -+ e_{i+3} (exact)")

w1, w2 = sp.symbols("w1 w2", positive=True)
w3 = 1 - w1 - w2
S = -(w1 * sp.log(w1) + w2 * sp.log(w2) + w3 * sp.log(w3))
grad = [sp.diff(S, v) for v in (w1, w2)]
Hs = sp.hessian(S, (w1, w2)).subs({w1: sp.Rational(1, 3), w2: sp.Rational(1, 3)})
ok = all(sp.simplify(g.subs({w1: sp.Rational(1, 3), w2: sp.Rational(1, 3)})) == 0 for g in grad)
ok &= sp.simplify(S.subs({w1: sp.Rational(1, 3), w2: sp.Rational(1, 3)}) - sp.log(3)) == 0
ok &= Hs.eigenvals() == {-9: 1, -3: 1}
Hgen = sp.hessian(S, (w1, w2))
# H = -(diag(1/w1, 1/w2) + (1/w3) 11^T): minus (positive definite + positive semidefinite), so S is strictly concave on the open simplex
ok &= sp.simplify(Hgen + sp.diag(1 / w1, 1 / w2) + sp.ones(2, 2) / w3) == sp.zeros(2)
Wh = (e[1] + e[2] + e[3]) / sp.sqrt(3)
ok &= sp.simplify(Wh * Wh - I8) == Z8
t = sp.symbols("t1:4")
ok &= sp.simplify(sum(t) / sp.sqrt(3) - sp.sqrt(3) * sum(t) / 3) == 0
perms = list(itertools.permutations(range(3)))
pt = [sp.Rational(1, 3)] * 3
ok &= all([pt[p[i]] for i in range(3)] == pt for p in perms)
wa, wb = sp.symbols("wa wb", positive=True)
fixed = sp.solve([wa - wb, wb - (1 - wa - wb)], [wa, wb], dict=True)
ok &= fixed == [{wa: sp.Rational(1, 3), wb: sp.Rational(1, 3)}]
report("5", ok, "prop:clock: the barycentre is the unique permutation-fixed point and a critical point of the entropy with value log 3 and Hessian "
       "eigenvalues -9, -3 in the chart (w1,w2); the entropy is strictly concave, so the maximum is unique; What^2 = 1; lambda = sqrt3 * mean t (exact)")


def cayley(a, b, c):
    A = sp.Matrix([[0, -c, b], [c, 0, -a], [-b, a, 0]])
    return (sp.eye(3) - A) * (sp.eye(3) + A).inv()


def BC(Cm):
    return sum((Cm[i, a] * e[i + 1] * e[a + 4] for i in range(3) for a in range(3)), Z8)


ok = True
for sig in [(5, 3, 1), (3, 2, 1), (4, 0, 0), (2, 2, 0)]:
    Uo, Vo = cayley(sp.Rational(1, 2), sp.Rational(1, 3), -1), cayley(2, sp.Rational(-1, 5), sp.Rational(1, 7))
    Cm = Uo * sp.diag(*sig) * Vo.T
    B = BC(Cm)
    B2 = B * B
    Bk = [sum((Uo[i, k] * e[i + 1] for i in range(3)), Z8) * sum((Vo[a, k] * e[a + 4] for a in range(3)), Z8) for k in range(3)]
    ok &= all(sp.simplify(b * b - I8) == Z8 for b in Bk) and all(sp.simplify(comm(Bk[k], Bk[l])) == Z8 for k in range(3) for l in range(3))
    ok &= sp.simplify(B - sum((sig[k] * Bk[k] for k in range(3)), Z8)) == Z8
    frob = sum(s ** 2 for s in sig)
    ok &= sp.simplify(B2 - frob * I8 - 2 * sum((sig[k] * sig[l] * Bk[k] * Bk[l] for k, l in [(0, 1), (0, 2), (1, 2)]), Z8)) == Z8
    ok &= sp.simplify(Bk[0] * Bk[1] * Bk[0] * Bk[2] - Bk[1] * Bk[2]) == Z8
    expected = {}
    for et in itertools.product((1, -1), repeat=3):
        val = (et[0] * sig[0] + et[1] * sig[1] + et[2] * sig[2]) ** 2
        expected[val] = expected.get(val, 0) + 1
    got = {sp.nsimplify(k): v for k, v in sp.Matrix(sp.simplify(B2)).eigenvals().items()}
    ok &= got == expected
    is_scalar = sp.simplify(B2 - B2[0, 0] * I8) == Z8
    ok &= is_scalar == (sum(1 for s in sig if s != 0) <= 1)
    has_zero = 0 in got
    s_sorted = sorted(sig, reverse=True)
    ok &= has_zero == (s_sorted[0] == s_sorted[1] + s_sorted[2])
for n in [sp.Matrix([1, 2, 2]) / 3, sp.Matrix([1, 1, 1]) / sp.sqrt(3), sp.Matrix([1, 0, 0])]:
    for a in range(3):
        Ga = sum((n[i] * e[i + 1] for i in range(3)), Z8) * e[a + 4]
        Cm = n * sp.Matrix([[1 if b == a else 0 for b in range(3)]])
        ok &= sp.simplify(Ga * Ga - I8) == Z8 and Cm.rank() == 1 and sorted(sp.simplify(x) for x in (Cm.T * Cm).eigenvals(multiple=True)) == [0, 0, 1]
report("6", ok, "thm:spectrumBC, cor:rankclock, ex:blindclock: for C = U diag(sigma) V^T (rational orthogonal U, V; sigma = (5,3,1), (3,2,1), (4,0,0), "
       "(2,2,0)) the B_k are commuting involutions, B_C = sum sigma_k B_k, B_C^2 = |C|_F^2 + 2 sum sigma_k sigma_l Omega_kl, Omega_12 Omega_13 = Omega_23, "
       "spec B_C^2 = {(eta.sigma)^2} with the stated multiplicities; B_C^2 is scalar iff rank C <= 1 and singular iff sigma_max = sum of the others; "
       "G_a(n)^2 = 1 with C of rank one and singular values (1,0,0) for three clocks including What (exact)")

# ---------------------------------------------------------------- Chapter 2, Dirac section
mono = [reduce(lambda a, b: a * b, [gam[i] for i in S], I8) for k in range(5) for S in itertools.combinations(range(4), k)]
comJ = solve_space([lambda X: X * J - J * X])
ok = len(comJ) == 32 and all(comm(g, J) == Z8 for g in gam)
ok &= span_eq(comJ, mono + [J * m for m in mono]) and sp.Matrix([list(x) for x in mono + [J * m for m in mono]]).rank() == 32
report("7", ok, "thm:realification (supports cor:nopredict): J commutes with the four gamma^mu; the 16 gamma monomials and their J-multiples are 32 linearly "
       "independent matrices spanning the commutant of J in M8(R), of dimension 32 (exact)")

Kmu = [(g0 * J) * (J * g) for g in gam]
ok = (g0 * J).T == -(g0 * J) and g0 * J == T
ok &= all(Kmu[m] == -g0 * gam[m] for m in range(4)) and Kmu[0] == -I8 and all(k.T == k for k in Kmu)
ok &= all(((g0 * J) * g).T == -((g0 * J) * g) for g in gam)
report("8", ok, "prop:kernels1: the mass form gamma^0 J = T is antisymmetric; K^mu = (gamma^0 J)(J gamma^mu) = -gamma^0 gamma^mu are symmetric with K^0 = -1; "
       "(gamma^0 J) gamma^mu are antisymmetric (exact)")

mass = solve_space([lambda X: X + X.T, lambda X: X * J - J * X] + [(lambda Lg: (lambda X: Lg.T * X + X * Lg))(Lg) for Lg in LOR])
ok = span_eq(mass, [T, T * Om])
ok &= all(((T * Om) * J * g).T == -((T * Om) * J * g) for g in gam)
ok &= Om * Om == -I8 and all(Om * g == -g * Om for g in gam)
com_mass = solve_space([lambda X: X + X.T, lambda X: X * J - J * X] + [(lambda Lg: (lambda X: Lg * X - X * Lg))(Lg) for Lg in LOR])
report("9", ok, "prop:massangle: the antisymmetric J-linear kernels A with Lorentz-invariant form (L^T A + A L = 0 for the six interior bivectors) are "
       f"exactly span{{T, T Omega}}; T Omega J gamma^mu is antisymmetric for every mu; Omega^2 = -1 anticommuting with every gamma^mu (exact). "
       f"Note: 'invariant' means invariance of the form; the antisymmetric J-linear elements that commute with the Lorentz algebra form a different "
       f"space, of dimension {len(com_mass)}")

p3, mm = (0, 0, 3), 4
Hex = sum((p3[a] * g0 * ga[a] for a in range(3)), Z8) + mm * g0
Sig3 = J * e[4] * e[5] / 2
ok = Hex * Hex == 25 * I8 and Hex.eigenvals() == {5: 4, -5: 4} and Hex.T == Hex and comm(Hex, J) == Z8
ok &= comm(Sig3, Hex) == Z8 and Sig3 * Sig3 == I8 / 4
for lam in (5, -5):
    Pl = (I8 + Hex / lam) / 2
    for s in (sp.Rational(1, 2), -sp.Rational(1, 2)):
        ok &= (Pl * (Sig3 + s * I8)).rank() == 2  # each sign of Sigma_3 picks a real 2-plane inside the 4-dim eigenspace
report("10", ok, "ex:Hp, eq:Hp: H = 3 e1e6 + 4 e1 is symmetric, commutes with J, H^2 = 25, eigenvalues +-5 of multiplicity 4; Sigma_3 = J e4e5/2 commutes with "
       "H, squares to 1/4 and splits each eigenspace into two real planes (exact)")

# ---------------------------------------------------------------- Chapter 2, discrete symmetries
Cs = solve_space([lambda X: X * J + J * X] + [(lambda g: (lambda X: X * J * g - J * g * X))(g) for g in gam])
ok = span_eq(Cs, [e[2], e[3]])
pv = sp.symbols("p1:4", real=True)
m_ = sp.symbols("m", positive=True)


def Hp(p, m):
    return sum((p[a] * g0 * ga[a] for a in range(3)), Z8) + m * g0


Hs_, Hm_ = Hp(pv, m_), Hp([-x for x in pv], m_)
Es = sp.sqrt(sum(x ** 2 for x in pv) + m_ ** 2)
Sa = [ga[1] * ga[2] / 2, ga[2] * ga[0] / 2, ga[0] * ga[1] / 2]
al = sp.symbols("alpha", real=True)
for Cc in [e[2], e[3], sp.cos(al) * e[2] + sp.sin(al) * e[3]]:
    Ci = sp.simplify(Cc.inv())
    ok &= sp.simplify(Cc * Hs_ * Ci + Hm_) == Z8 and sp.simplify(Cc * J * Ci + J) == Z8 and sp.simplify(Cc * Gam * Ci + Gam) == Z8
    ok &= sp.simplify(Cc * (I8 + Hs_ / Es) / 2 * Ci - (I8 - Hm_ / Es) / 2) == Z8
    ok &= all(sp.simplify(comm(Cc, Lg)) == Z8 for Lg in LOR)
    ok &= sp.simplify(Cc.T * T * Cc - T) == Z8 and all(sp.simplify(Cc.T * k * Cc - k) == Z8 for k in Kmu)
    ok &= sp.simplify(-sum(((Cc * s * Ci) ** 2 for s in Sa), Z8) - sp.Rational(3, 4) * I8) == Z8
    ok &= all(sp.simplify(Cc.T * (k * J + (k * J).T) * Cc + (k * J + (k * J).T)) == Z8 for k in Kmu)
    ok &= all(sp.simplify(Cc.T * (k * J) * Cc + k * J) == Z8 for k in Kmu)
ok &= -sum((s * s for s in Sa), Z8) == sp.Rational(3, 4) * I8
report("11", ok, "thm:C and its corollary: the antilinear intertwiners of J gamma^mu d_mu - m are exactly span{e2, e3}; every unit C in U satisfies C H(p) C^-1 = "
       "-H(-p), C J C^-1 = -J, C Gamma C^-1 = -Gamma, C Lambda_+(p) C^-1 = Lambda_-(-p) (symbolic p, m), commutes with the interior Lorentz algebra, "
       "preserves T and every K^mu, preserves the Casimir -sum S_a^2 = 3/4; C^T (K^mu J) C = -K^mu J, so the current Psi^T K^mu J Psi vanishes when "
       "C Psi = Psi (exact)")

Bv = ga[0] * ga[1] * ga[2]
par = solve_space([lambda X: X * J - J * X, lambda X: X * g0 - g0 * X] + [(lambda g: (lambda X: X * g0 * g + g0 * g * X))(g) for g in ga])
trv = solve_space([lambda X: X * J + J * X, lambda X: X * g0 - g0 * X] + [(lambda g: (lambda X: X * g0 * g + g0 * g * X))(g) for g in ga])
ok = span_eq(par, [g0, g0 * J]) and span_eq(trv, [e[2] * Bv, e[3] * Bv])
a_, b_ = sp.symbols("a b", real=True)
P2 = sp.expand((a_ * g0 + b_ * g0 * J) ** 2)
ok &= sp.simplify(P2 - ((a_ ** 2 - b_ ** 2) * I8 + 2 * a_ * b_ * J)) == Z8
ok &= (g0 * J) ** 2 == -I8 and Bv * Bv == I8
Cg = sp.cos(al) * e[2] + sp.sin(al) * e[3]
ok &= sp.simplify(Cg * Bv + Bv * Cg) == Z8 and sp.simplify((Cg * Bv) ** 2 + I8) == Z8
ok &= all(sp.simplify(x * Hs_ - Hm_ * x) == Z8 for x in [g0, g0 * J, e[2] * Bv, e[3] * Bv])
report("12", ok, "thm:PT: the linear X with X H(p) = H(-p) X for all p, m form span{gamma^0, gamma^0 J}; (a gamma^0 + b gamma^0 J)^2 = (a^2-b^2) + 2ab J, "
       "so P^2 = +1 only at +-gamma^0 and P = +-gamma^0 J gives -1; the antilinear ones form U.gamma^1gamma^2gamma^3; B^2 = 1, B anticommutes with U, "
       "T^2 = -|C|^2 (exact)")

ok = True
for Cc in [e[2], e[3], Cg]:
    Tt = Cc * Bv
    ok &= sp.simplify(Cc * g0 * Tt + Om) == Z8 and sp.simplify(Cc * g0 * J * Tt - Gam) == Z8
ok &= Om * Om == -I8 and Gam * Gam == I8 and Gam == bl(1, 2, 3, 4, 5, 6)
ok &= sp.simplify((sp.cos(th) * I8 + sp.sin(th) * J) * Om - (sp.cos(th) * Om + sp.sin(th) * J * Om)) == Z8
report("13", ok, "thm:CPT: with C any unit vector of U and T = C gamma^1gamma^2gamma^3, P = gamma^0 gives CPT = -Omega with (CPT)^2 = -1 and P = gamma^0 J gives "
       "CPT = +Gamma = e1...e6 with (CPT)^2 = +1; e^{theta J} Omega = cos theta Omega + sin theta J Omega (exact)")

gs = solve_space([(lambda k: (lambda X: X.T * k + k * X))(k) for k in Kmu])
gens = [Om, J, e[2] * Om, e[3] * Om]
ok = span_eq(gs, gens)
ok &= all(sp.Matrix([(x * y).trace() for y in gens]) == sp.Matrix([-8 if x is y else 0 for y in gens]) for x in gens)
su2 = gens[1:]
ok &= all(span_eq(su2, su2 + [comm(x, y)]) for x in su2 for y in su2) and all(comm(Om, x) == Z8 for x in su2)
ok &= any(comm(x, y) != Z8 for x in su2 for y in su2)
gs_flip = solve_space([(lambda k: (lambda X: X.T * (-k) + (-k) * X))(k) for k in Kmu])
ok &= span_eq(gs, gs_flip)
gm = solve_space([(lambda k: (lambda X: X.T * k + k * X))(k) for k in Kmu] + [lambda X: X.T * T + T * X])
ok &= span_eq(gm, [J])
poles = [e[1] * e[4], e[2] * e[5], e[3] * e[6]]
ok &= sp.Matrix([list(x) for x in gens + poles]).rank() == 7
ok &= all(comm(x, Lg) == Z8 for x in gens for Lg in LOR) and all(x.T == -x for x in gens)
report("14", ok, "prop:gauge1: the solutions of X^T K^mu + K^mu X = 0 over all of M8(R) are span{Omega, J, e2 Omega, e3 Omega}, antisymmetric and commuting with "
       "the interior Lorentz algebra; trace form -8 id on this basis; {J, e2 Omega, e3 Omega} closes into su(2) with Omega central; the space is "
       "unchanged under K^mu -> -K^mu (the P* convention of Chapters 6-7, K^0 = +1); adding mass-form invariance leaves R J; it meets the pole span in 0 (exact)")


# ---------------------------------------------------------------- the module (thm:module, thm:spectra, prop:Cstates)
P = {"I": sp.eye(2), "X": sp.Matrix([[0, 1], [1, 0]]), "Z": sp.Matrix([[1, 0], [0, -1]]), "E": sp.Matrix([[0, -1], [1, 0]])}
words = [w for w in itertools.product("IXZE", repeat=4) if w.count("E") % 2 == 1]


def anti(u, v):
    return sum(1 for a, b in zip(u, v) if a != "I" and b != "I" and a != b) % 2 == 1


def clique(cur, start):
    if len(cur) == 8:
        return cur
    for i in range(start, len(words)):
        if all(anti(words[i], w) for w in cur):
            r = clique(cur + [words[i]], i + 1)
            if r:
                return r
    return None


W8 = clique([], 0)
I16 = sp.eye(16)
th8 = [reduce(sp.kronecker_product, [P[c] for c in w]) for w in W8]
ok = all(x.T == -x and x * x == -I16 for x in th8)
ok &= all(th8[i] * th8[j] + th8[j] * th8[i] == (-2 * I16 if i == j else sp.zeros(16)) for i in range(8) for j in range(8))
mons = [reduce(lambda a, b: a * b, [th8[i] for i in S], I16) for k in range(9) for S in itertools.combinations(range(8), k)]
gram_ok = all((mons[i].T * mons[j]).trace() == (16 if i == j else 0) for i in range(256) for j in range(i, 256))
ok &= len(mons) == 256 and gram_ok
ok &= all((th8[i] * th8[j]).T == -(th8[i] * th8[j]) for i in range(8) for j in range(i + 1, 8))


def rho(A):
    return -sum((A[i, j] * th8[i] * th8[j] for i in range(8) for j in range(8)), sp.zeros(16)) / 4


ok &= all(sp.simplify(rho(comm(A, B)) - comm(rho(A), rho(B))) == sp.zeros(16) for A, B in [(J, Om), (J, e[2] * Om), (e[2] * Om, e[3] * Om), (J, J * Hex)])
ok &= all(sp.simplify(comm(rho(A), th8[k]) - sum((A[i, k] * th8[i] for i in range(8)), sp.zeros(16))) == sp.zeros(16) for A in [J, Om] for k in range(8))
report("15", ok, f"thm:module: eight integer 16x16 matrices (Kronecker words {' '.join(''.join(w) for w in W8)}) are antisymmetric with theta^2 = -1 and pairwise "
       "anticommuting; the 256 monomials are trace-orthogonal (Gram = 16 id), so they span M16(R) and the commutant is R id; theta_i theta_j (i<j) are "
       "antisymmetric; rho(A) = -1/4 sum A_ij theta_i theta_j is a Lie homomorphism with [rho(A), theta_k] = sum_i A_ik theta_i (exact)")

x = sp.symbols("x")
RJ = rho(J)
ok = sp.factor((2 * RJ).charpoly(x).as_expr()) == sp.factor(x ** 6 * (x ** 2 + 4) ** 4 * (x ** 2 + 16))
for p_, m0 in [((0, 0, 3), 4), ((1, 4, 8), 4), ((2, 3, 6), 0)]:
    E2 = sum(q ** 2 for q in p_) + m0 ** 2
    RH = rho(J * Hp(p_, m0))
    ok &= sp.expand((2 * RH).charpoly(x).as_expr() - x ** 6 * (x ** 2 + 4 * E2) ** 4 * (x ** 2 + 16 * E2)) == 0
    ok &= comm(RJ, RH) == sp.zeros(16)
RH = rho(J * Hex)
table = {}
for h in (-10, -5, 0, 5, 10):
    for q in (-2, -1, 0, 1, 2):
        Mjoint = sp.Matrix.vstack(RH - sp.I * h * I16, RJ - sp.I * q * I16)
        d = 16 - Mjoint.rank()
        if d:
            table[(h // 5, q)] = d
expected = {(-2, 0): 1, (-1, 1): 2, (-1, -1): 2, (0, 2): 1, (0, -2): 1, (0, 0): 4, (1, 1): 2, (1, -1): 2, (2, 0): 1}
ok &= table == expected
report("16", ok, "thm:spectra: chi_{2rho(J)} = x^6 (x^2+4)^4 (x^2+16); chi_{2rho(JH)} = x^6 (x^2+4E^2)^4 (x^2+16E^2) for three momenta (E^2 = 25, 97, 49); "
       "rho(J) and rho(JH) commute; at E = 5 the joint (H/E, Q) eigenspace dimensions on F are exactly the table "
       f"{sorted(table.items())} (exact)")

# lift of C = e2 to Spin(8): product of four reflections in an orthonormal basis of its -1 eigenspace
Cm8 = e[2]
ok = Cm8.T * Cm8 == I8 and Cm8.det() == 1 and len((Cm8 + I8).nullspace()) == 4
vs = sp.Matrix.hstack(*(Cm8 + I8).nullspace())
vs = sp.Matrix.hstack(*sp.GramSchmidt([vs[:, k] for k in range(4)], True))
refl = [I8 - 2 * vs[:, k] * vs[:, k].T for k in range(4)]
ok &= all(all(int(v) == v for v in (sp.eye(8) - (vs[:, k] * vs[:, k].T) * 2)) for k in range(4))
ok &= sp.simplify(reduce(lambda a, b: a * b, refl) - Cm8) == Z8
Chat = sp.simplify(reduce(lambda a, b: a * b, [sum((vs[i, k] * th8[i] for i in range(8)), sp.zeros(16)) for k in range(4)]))
ok &= sp.simplify(Chat.T * Chat - I16) == sp.zeros(16) and sp.simplify(Chat * Chat - I16) == sp.zeros(16)
Ci = Chat.T
ok &= all(sp.simplify(Chat * th8[k] * Ci - sum((Cm8[i, k] * th8[i] for i in range(8)), sp.zeros(16))) == sp.zeros(16) for k in range(8))
ok &= sp.simplify(Chat * RJ * Ci + RJ) == sp.zeros(16)
RHm = rho(J * Hp((0, 0, -3), 4))
ok &= sp.simplify(Chat * RH * Ci - RHm) == sp.zeros(16)
plus = (Chat - I16).nullspace()
minus = (Chat + I16).nullspace()
kp = 16 - sp.Matrix.vstack(Chat - I16, RJ).rank()
km = 16 - sp.Matrix.vstack(Chat + I16, RJ).rank()
ok &= len(plus) == 8 and len(minus) == 8 and kp == 3 and km == 3
report("17", ok, "prop:Cstates: C = e2 is in SO(8) with a 4-dim -1 eigenspace, the product of four integer reflection matrices; its lift Chat (product of the "
       "four unit Clifford vectors) is orthogonal with Chat^2 = 1, implements C on the theta_i, anticommutes with Q and carries rho(JH(p)) to rho(JH(-p)); "
       f"the +1 eigenspace has real dimension {len(plus)} and Q annihilates a {kp}-dimensional part of it (the -1 eigenspace likewise: {len(minus)}, {km}) (exact)")

RO = rho(Om)
ok = comm(RJ, RO) == sp.zeros(16) and (RJ * RO).rank() == 8
joint = {}
for qj in range(-2, 3):
    for qo in range(-2, 3):
        d = 16 - sp.Matrix.vstack(RJ - sp.I * qj * I16, RO - sp.I * qo * I16).rank()
        if d:
            joint[(qj, qo)] = d
ok &= sum(joint.values()) == 16
ok &= all((qj + qo) % 2 == 0 for (qj, qo) in joint)
report("18", ok, "remark (many-body chirality): rho(J) rho(Omega) has rank 8; rho(J), rho(Omega) commute and every joint eigenvector has q_J + q_Omega even, "
       f"so Gamma_F = exp[pi/2 (rho(J)+rho(Omega))] with eigenvalue i^(q_J+q_Omega) is an involution; joint multiplicities {sorted(joint.items())} (exact)")

# ---------------------------------------------------------------- the K^mu sign convention (Chapter 2 versus Chapters 6-7)
Hfl = -RH
ok = sp.factor(Hfl.charpoly(x).as_expr()) == sp.factor(RH.charpoly(x).as_expr())
ok &= span_eq(gs, gs_flip)
psi = sp.Matrix(sp.symbols("u0:8"))
dpsi = [sp.Matrix(sp.symbols(f"d{m}_0:8")) for m in range(4)]
mval = sp.Symbol("mu_m")
eq_ch2 = sum((Kmu[m] * dpsi[m] for m in range(4)), sp.zeros(8, 1)) + mval * T * psi
eq_ch6 = sum((-(Kmu[m]) * dpsi[m] for m in range(4)), sp.zeros(8, 1)) - mval * T * psi
ok &= sp.expand(eq_ch2 + eq_ch6) == sp.zeros(8, 1)
report("19", ok, "K^mu sign convention: Chapter 2 uses K^mu = -gamma^0 gamma^mu (K^0 = -1), Chapters 6-7 the compressed P* gamma^0 gamma^mu P* (K^0 = +1). An overall "
       "sign of the fermion Lagrangian leaves the field equation unchanged, leaves the massless gauge algebra unchanged (a1-14) and leaves the "
       "one-momentum spectrum unchanged (H -> -H has the same characteristic polynomial, the spectrum being symmetric); a convention, reconciled in "
       "the text after r1:prop:dirac (exact)")

print(f"{NPASS[0]}/{NPASS[0] + len(FAILS)} checks passed")
sys.exit(0 if not FAILS else 1)
