#!/usr/bin/env python3
"""g1: the curvature operator on the frame, the phase family, the Helmholtz condition and the Palatini-Holst reading (Phase 9C).

    python3 verification/ch11_curvature/g1_curvature_operator.py

Chapter 11, Section "Curvature in the field-strength reading" (prop:curvop, thm:phase, ex:dS, cor:phasesel, thm:helmholtz, prop:holst).
Conventions (from the chapter and from thm:emfield): signature (+,-,-,-); the frame Lambda^2 R^{1,3} in the basis
(e^23, e^31, e^12, e^01, e^02, e^03), which is the basis (e_1,...,e_6) of R^{3,3} with B on the first and E on the last three; induced
metric <w,w'> = (1/2) w_ab w'^ab = diag(1,1,1,-1,-1,-1); Hodge star (*w)_ab = (1/2) eps_abcd w^cd with eps_0123 = +1, which is G
(G e_i = e_{i+3}); (R w)_ab = (1/2) R_ab^cd w_cd; Ric(R)_bd = R^a_bad; R = g^bd Ric_bd.
Exact symbolic arithmetic (sympy) on a general operator with 21 (or 36) free entries: every PASS is an identity, not a sample.
"""
import itertools
import sys

import sympy as sp

NPASS, FAILS = [0], []


def report(cid, ok, text):
    if ok:
        NPASS[0] += 1
        print(f"PASS g1-{cid}: {text}")
    else:
        FAILS.append(cid)
        print(f"FAIL g1-{cid}: {text}")
    sys.stdout.flush()


g = sp.diag(1, -1, -1, -1)
ginv = g.inv()
pairs = [(2, 3), (3, 1), (1, 2), (0, 1), (0, 2), (0, 3)]
eta6 = sp.diag(1, 1, 1, -1, -1, -1)
eps = lambda *ix: sp.LeviCivita(*ix)                       # eps_{0123} = +1 (lower indices)


def basis2(I):
    """lower components w_ab of the I-th basis 2-form."""
    a, b = pairs[I]
    w = sp.zeros(4, 4); w[a, b] = 1; w[b, a] = -1
    return w


def raise2(w):
    return ginv * w * ginv


def coords2(w):
    """coordinates of a 2-form (lower components) in the basis."""
    return sp.Matrix([w[a, b] for a, b in pairs])


# induced metric and Hodge star in the basis
ok = all(sp.Rational(1, 2) * sum(basis2(I)[a, b] * raise2(basis2(J))[a, b] for a in range(4) for b in range(4)) == eta6[I, J]
         for I in range(6) for J in range(6))


def hodge(w):
    wu = raise2(w)
    return sp.Matrix(4, 4, lambda a, b: sp.Rational(1, 2) * sum(eps(a, b, c, d) * wu[c, d] for c in range(4) for d in range(4)))


G = sp.Matrix.hstack(*[coords2(hodge(basis2(J))) for J in range(6)])
ok &= G == sp.Matrix([[0, 0, 0, -1, 0, 0], [0, 0, 0, 0, -1, 0], [0, 0, 0, 0, 0, -1], [1, 0, 0, 0, 0, 0], [0, 1, 0, 0, 0, 0], [0, 0, 1, 0, 0, 0]])
ok &= G * G == -sp.eye(6) and (eta6 * G).T == eta6 * G
report("1", ok, "conventions: in the basis (e^23,e^31,e^12,e^01,e^02,e^03) the induced metric (1/2)w_ab w'^ab is diag(1,1,1,-1,-1,-1) = eta, and the Hodge "
       "star with eps_0123 = +1 is G (G e_i = e_{i+3}, G e_{i+3} = -e_i), G^2 = -1, eta-self-adjoint (exact)")


def tensor(M):
    """R_abcd (all lower) of the operator M: R_ab^{cd} = sum_I M_IJ (E_I)_ab for the pair J = (c,d)."""
    Rud = [[[[sp.S(0)] * 4 for _ in range(4)] for _ in range(4)] for _ in range(4)]
    for J, (c, d) in enumerate(pairs):
        for a in range(4):
            for b in range(4):
                val = sum(M[I, J] * basis2(I)[a, b] for I in range(6))
                Rud[a][b][c][d] = val; Rud[a][b][d][c] = -val
    return [[[[sp.expand(sum(g[c, cc] * g[d, dd] * Rud[a][b][cc][dd] for cc in range(4) for dd in range(4)))
               for d in range(4)] for c in range(4)] for b in range(4)] for a in range(4)]


def ric(R):
    return sp.Matrix(4, 4, lambda b, d: sp.expand(sum(ginv[a, c] * R[c][b][a][d] for a in range(4) for c in range(4))))


def scal(Ric):
    return sp.expand(sum(ginv[b, d] * Ric[b, d] for b in range(4) for d in range(4)))


m = sp.symbols("m0:36")
Mgen = sp.Matrix(6, 6, m)
# pair symmetry <=> eta-self-adjoint
Rg = tensor(Mgen)
pair_eqs = [sp.expand(Rg[a][b][c][d] - Rg[c][d][a][b]) for a, b, c, d in itertools.product(range(4), repeat=4)]
sa_eqs = list(eta6 * Mgen - (eta6 * Mgen).T)
sol_pair = sp.linsolve([e for e in pair_eqs if e != 0], m)
sol_sa = sp.linsolve([e for e in sa_eqs if e != 0], m)
ok = sol_pair == sol_sa
free = len(set().union(*[sp.sympify(x).free_symbols for x in list(sol_pair)[0]]))
ok &= free == 21
# a general element of S (21 parameters): M = eta^{-1} Ssym with Ssym symmetric
s = sp.symbols("s0:21")
Ssym = sp.zeros(6, 6)
for k_, (i, j) in enumerate([(i, j) for i in range(6) for j in range(i, 6)]):
    Ssym[i, j] = Ssym[j, i] = s[k_]
M = eta6 * Ssym                                             # eta^{-1} = eta
R = tensor(M); Ric = ric(R); Rs = scal(Ric)
trM = sp.expand(M.trace()); trGM = sp.expand((G * M).trace())
ok &= sp.expand(trM - Rs / 2) == 0
cyc_ok = all(sp.expand(R[a][b][c][d] + R[a][c][d][b] + R[a][d][b][c] + sp.Rational(1, 2) * trGM * eps(a, b, c, d)) == 0
             for a, b, c, d in itertools.product(range(4), repeat=4))
ok &= cyc_ok
GG = (G * G).trace()
ok &= GG == -6 and (eta6 * G).T == eta6 * G
report("2", ok, f"prop:curvop: R_abcd = R_cdab iff the operator is eta-self-adjoint (the two linear systems on the 36 entries have the same solution "
       f"space, of dimension {free} = 21); on that space Tr R = R/2 and R_abcd + R_acdb + R_adbc = -(1/2) Tr(G R) eps_abcd identically, so the "
       "first Bianchi identity is the single condition Tr(G R) = 0 (dimension 20); G is eta-self-adjoint with Tr(G G) = -6, so G is in S but not "
       "in V (exact, general 21-parameter operator)")

# thm:phase
th = sp.symbols("theta", real=True)
c_, s_ = sp.cos(th), sp.sin(th)
E = c_ * sp.eye(6) + s_ * G
Mth = E * M * E
ok = (eta6 * Mth).T == sp.expand(eta6 * Mth) or sp.simplify((eta6 * Mth).T - eta6 * Mth) == sp.zeros(6, 6)
Ric_th = ric(tensor(Mth))
target = Ric - sp.Rational(1, 2) * s_ ** 2 * Rs * g + sp.Rational(1, 2) * sp.sin(2 * th) * trGM * g
diff = (Ric_th - target).applyfunc(lambda e: sp.simplify(sp.expand_trig(sp.expand(e))))
ok &= diff == sp.zeros(4, 4)
lem1 = ric(tensor(G * M * G)) - (Ric - Rs / 2 * g)
lem2 = ric(tensor(G * M)) + ric(tensor(M * G)) - trGM * g
ok &= lem1.applyfunc(sp.expand) == sp.zeros(4, 4) and lem2.applyfunc(sp.expand) == sp.zeros(4, 4)
report("3", ok, "thm:phase: for every eta-self-adjoint operator R (21 free entries, on or off the Bianchi shell) and every theta, e^{theta G} R e^{theta G} is "
       "eta-self-adjoint and Ric(R_theta)_ab = R_ab - (1/2) sin^2(theta) R g_ab + (1/2) sin(2 theta) Tr(G R) g_ab; the two lemmas Ric(GRG) = Ric(R) - (R/2) g "
       "and Ric(GR) + Ric(RG) = Tr(GR) g hold (exact symbolic identities)")

# ex:dS
k = sp.symbols("k")
RdS = [[[[k * (g[a, c] * g[b, d] - g[a, d] * g[b, c]) for d in range(4)] for c in range(4)] for b in range(4)] for a in range(4)]
MdS = sp.zeros(6, 6)
for I, (a, b) in enumerate(pairs):
    for J, (c, d) in enumerate(pairs):
        MdS[I, J] = sum(ginv[c, cc] * ginv[d, dd] * RdS[a][b][cc][dd] for cc in range(4) for dd in range(4))
ok = MdS == k * sp.eye(6) and tensor(MdS) == [[[[sp.expand(x) for x in y] for y in z] for z in w] for w in RdS]
RicdS = ric(RdS); RdSs = scal(RicdS)
ok &= RdSs == 12 * k and sp.expand((G * MdS).trace()) == 0 and G * MdS * G == -k * sp.eye(6)
ok &= ric(tensor(G * MdS * G)) == (RicdS - RdSs / 2 * g) and RicdS - RdSs / 2 * g == -3 * k * g
report("4", ok, "ex:dS: R_abcd = k(g_ac g_bd - g_ad g_bc) has curvature operator +k id (with these conventions), R = 12k, Tr(G R) = 0, G R G = -k id, "
       "and Ric(GRG) = R_ab - (R/2) g_ab = -3k g_ab (exact)")

# cor:phasesel (ii), (iii)
trGMth = sp.simplify(sp.expand_trig(sp.expand((G * Mth).trace())))
ok = sp.simplify(trGMth - (sp.cos(2 * th) * trGM - sp.Rational(1, 2) * sp.sin(2 * th) * Rs)) == 0
Qv = sp.expand(sp.Rational(1, 2) * (M * G * M * G).trace())
Ru = [[[[sum(ginv[a, aa] * ginv[b, bb] * ginv[c, cc] * ginv[d, dd] * R[aa][bb][cc][dd] for aa in range(4) for bb in range(4) for cc in range(4)
              for dd in range(4)) for d in range(4)] for c in range(4)] for b in range(4)] for a in range(4)]
Riem2 = sp.expand(sum(R[a][b][c][d] * Ru[a][b][c][d] for a, b, c, d in itertools.product(range(4), repeat=4)))
Ricu = ginv * Ric * ginv
Ric2 = sp.expand(sum(Ric[a, b] * Ricu[a, b] for a in range(4) for b in range(4)))
ok &= sp.expand(Qv + sp.Rational(1, 8) * (Rs ** 2 - 4 * Ric2 + Riem2)) == 0
Mbi = M.subs(s[3], sp.solve(trGM, s[3])[0]) if trGM.has(s[3]) else None
ok &= Mbi is not None and sp.expand((G * (G * Mbi * G)).trace()) == 0
report("5", ok, "cor:phasesel (ii) Tr(G R_theta) = cos(2 theta) Tr(G R) - (1/2) sin(2 theta) R; (iii) Q(R) = (1/2) Tr(R G R G) = -(1/8)(R^2 - 4 R_ab R^ab + "
       "R_abcd R^abcd) for every pair-symmetric R, and G R G stays on the Bianchi shell when R does (exact)")

# thm:helmholtz: antisymmetric part of the linearisation at flat space on plane waves
kk = sp.symbols("k0:4")
hs = sp.symbols("h0:10"); qs = sp.symbols("q0:10")
def symm(v):
    Mx = sp.zeros(4, 4); n = 0
    for i in range(4):
        for j in range(i, 4):
            Mx[i, j] = Mx[j, i] = v[n]; n += 1
    return Mx
h = symm(hs); q = symm(qs)
kl = sp.Matrix(kk)                                          # lower components k_a
ku = ginv * kl


def dH(hm):
    """symbol of delta H_ab(theta) at g = eta on h_ab e^{ik.x} (d -> ik, overall sign of i^2 = -1 included)."""
    hu = ginv * hm                                          # h^c_b
    trh = (ginv * hm).trace()
    k2 = (kl.T * ku)[0]
    dR = sp.Matrix(4, 4, lambda a, b: sp.Rational(1, 2) * (-sum(kl[c] * kl[a] * hu[c, b] for c in range(4))
                                                          - sum(kl[c] * kl[b] * hu[c, a] for c in range(4)) + k2 * hm[a, b] + kl[a] * kl[b] * trh))
    dRs = sum(ginv[a, b] * dR[a, b] for a in range(4) for b in range(4))
    return dR - sp.Rational(1, 2) * s_ ** 2 * dRs * g


qu = ginv * q * ginv; hu_ = ginv * h * ginv
asym = sp.expand(sum(qu[a, b] * dH(h)[a, b] - hu_[a, b] * dH(q)[a, b] for a in range(4) for b in range(4)))
kqk = (ku.T * q * ku)[0]; khk = (ku.T * h * ku)[0]
target = sp.Rational(1, 2) * c_ ** 2 * (kqk * (ginv * h).trace() - khk * (ginv * q).trace())
ok = sp.simplify(sp.expand(asym - target).subs(s_ ** 2, 1 - c_ ** 2)) == 0
ok &= sp.expand(asym.subs(th, sp.pi / 2)) == 0
report("6", ok, "thm:helmholtz: the linearisation of H(theta) = Ric - (1/2) sin^2(theta) R g at g = eta on plane waves has antisymmetric part "
       "q^ab dH_ab[h] - h^ab dH_ab[q] = (1/2) cos^2(theta) [(k.q.k) tr h - (k.h.k) tr q], which vanishes identically only at cos(theta) = 0 (exact, general "
       "symmetric h, q and covector k). With the necessity half of the Helmholtz conditions (inherited, Vainberg, Olver) no Lagrangian exists for "
       "cos(theta) != 0; at cos(theta) = 0 the Lagrangian R sqrt|g| is the classical one. This is a statement about the linearised operator at flat space")

# prop:holst: forms at a point, orthonormal coframe e^a = dx^a
def four(c, d, e_, f):
    return eps(c, d, e_, f)                                 # e^c e^d e^e e^f = [cdef] vol


def wedge_pair(X, Y):
    """<X ^ Y> = (1/2) X_ab ^ Y^ab for Lambda^2-valued 2-forms with components X[a][b][c][d] (internal ab upper, form cd), coefficient of vol."""
    tot = 0
    for a, b, c, d, e_, f in itertools.product(range(4), repeat=6):
        if a == b or c == d or e_ == f:
            continue
        Xl = g[a, a] * g[b, b] * X[a][b][c][d]
        tot += sp.Rational(1, 2) * sp.Rational(1, 4) * Xl * Y[a][b][e_][f] * four(c, d, e_, f)
    return sp.expand(tot)


def Sigma():
    return [[[[(1 if (a, b) == (c, d) else 0) - (1 if (a, b) == (d, c) else 0) for d in range(4)] for c in range(4)] for b in range(4)] for a in range(4)]


def Fform(Mx):
    Rt = tensor(Mx)
    return [[[[sp.expand(ginv[a, a] * ginv[b, b] * Rt[a][b][c][d]) for d in range(4)] for c in range(4)] for b in range(4)] for a in range(4)]


def Gint(X):
    """Hodge star on the internal (upper) indices: (GX)^ab = (1/2) eps^ab_cd X^cd."""
    return [[[[sp.expand(sum(sp.Rational(1, 2) * ginv[a, a] * ginv[b, b] * eps(a, b, p, r) * X[p][r][c][d] for p in range(4) for r in range(4)))
               for d in range(4)] for c in range(4)] for b in range(4)] for a in range(4)]


Sg, Fg = Sigma(), Fform(M)
pal = sp.expand(sum(sp.Rational(1, 4) * eps(a, b, c, d) * sp.Rational(1, 2) * Fg[c][d][e_][f] * four(a, b, e_, f)
                    for a, b, c, d, e_, f in itertools.product(range(4), repeat=6) if len({a, b, c, d}) == 4 and e_ != f))
GS = Gint(Sg)
ok = sp.expand(wedge_pair(GS, Fg) - pal) == 0 and sp.expand(pal - trM) == 0
hol = sp.expand(sum(g[a, a] * g[b, b] * sp.Rational(1, 2) * Fg[a][b][e_][f] * four(a, b, e_, f) for a, b, e_, f in itertools.product(range(4), repeat=4)
                    if a != b and e_ != f))
ok &= sp.expand(hol - 2 * wedge_pair(Sg, Fg)) == 0 and sp.expand(hol + 2 * trGM) == 0
Eth = lambda X: [[[[sp.expand(c_ * X[a][b][c][d] + s_ * Gint(X)[a][b][c][d]) for d in range(4)] for c in range(4)] for b in range(4)] for a in range(4)]
lhs = wedge_pair(Eth(Sg), Eth(Fg))
ok &= sp.simplify(sp.expand_trig(lhs - (sp.cos(2 * th) * wedge_pair(Sg, Fg) + sp.sin(2 * th) * wedge_pair(GS, Fg)))) == 0
ok &= sp.simplify(sp.expand_trig(lhs + trGMth)) == 0
report("7", ok, "prop:holst (i)-(iii): (1/4) eps_abcd e^a e^b F^cd = <G Sigma ^ F> = Tr(R) vol; e_a e_b F^ab = 2 <Sigma ^ F> = -2 Tr(G R) vol; "
       "<e^{theta G} Sigma ^ e^{theta G} F> = cos(2 theta) <Sigma ^ F> + sin(2 theta) <G Sigma ^ F> = -Tr(G R_theta) vol, for every eta-self-adjoint R "
       "(exact; forms at a point with an orthonormal coframe, F^ab = (1/2) R^ab_cd e^c e^d)")

# gamma = -tan(2 theta): L_theta = sin2t * P + cos2t * (Ho/2) = (sin2t/2) [ (1/2) eps e e F - gamma^{-1} e e F ]  with  -1/gamma = cot(2 theta)
gam = sp.symbols("gamma")
L_theta = sp.sin(2 * th) * pal + sp.cos(2 * th) * hol / 2
L_holst = sp.sin(2 * th) / 2 * (2 * pal - hol / gam)
ok = sp.simplify(L_theta - L_holst.subs(gam, -sp.tan(2 * th))) == 0
# Holst term inert on the Bianchi shell: e_b ^ F^ab = 0 iff Tr(G R) = 0 (3-form components)
three = []
for a in range(4):
    for p, q_, r in itertools.combinations(range(4), 3):
        val = 0
        for b in range(4):
            for c, d in itertools.permutations(range(4), 2):
                perm = [b, c, d]
                if sorted(perm) == [p, q_, r]:
                    val += g[b, b] * sp.Rational(1, 2) * Fg[a][b][c][d] * sp.combinatorics.Permutation([perm.index(x) for x in [p, q_, r]]).signature()
        three.append(sp.expand(val))
ratio_ok = all(sp.simplify(t.subs(s[3], sp.solve(trGM, s[3])[0])) == 0 for t in three)
ok &= ratio_ok and any(t != 0 for t in three)
report("8", ok, "prop:holst (iii)-(iv): the first-order Lagrangian <e^{theta G} Sigma ^ e^{theta G} F> equals (sin 2theta / 2)[(1/2) eps e e F - gamma^{-1} e e F] "
       "with gamma = -tan(2 theta) (Holst's normalisation); the Holst 3-forms e_b ^ F^ab vanish exactly on the Bianchi shell Tr(G R) = 0, so there the "
       "Holst term does not contribute to the coframe equation, which is sin(2 theta) times the Palatini (Einstein) equation (exact). The connection "
       "equation and the variation itself are not recomputed here")

print(f"{NPASS[0]}/{NPASS[0] + len(FAILS)} checks passed")
sys.exit(0 if not FAILS else 1)
