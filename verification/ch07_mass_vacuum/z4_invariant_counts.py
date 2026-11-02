#!/usr/bin/env python3
"""z4: independent count of the Pati-Salam-invariant polynomials on the kernel space K (CORE-09), 1 November 2026.

    python3 verification/core/z4_invariants.py

Claim tested (tex/s20_selection.tex, subsection 'The full kernel space'): the numbers of linearly independent Pati-Salam-invariant real
polynomials on K of degrees 1, 2, 3, 4 are 0, 18, 0, 686, and with fermion number imposed 0, 12, 0, 402.

Method (independent of any matrix realisation): K (x) C = V + conj(V) with V the complex representation of thm:kernelspace(ii)
(complex structure C of Phase 3), V = (10,3,1) + (10,1,3) + (15,2,2) + 2(6,1,1) + (1,2,2) of su(4)+su(2)_L+su(2)_R, and fermion number
acting with charge 2 on the Majorana-type components (10,3,1), (10,1,3), 2(6,1,1) and 0 on the Dirac-type ones (prop:fermionnumber).
dim Sym^d(K)^G = (1/|W|) CT[ h_d(chi) * prod_{roots}(1 - e^alpha) ] (Molien-Weyl), h_d from the power sums p_k(x) = chi(x^k) by Newton's
identities. All Laurent polynomials are exact integer dictionaries; only the coefficients needed by the Weyl factor are formed.
Also computed: dim Sym^d(V (x) C)^G for the complex V alone is NOT what is claimed; the claim is about real polynomials on K.
Prints PASS/FAIL lines and exits 0 iff all pass.
"""
import itertools
import sys
from collections import Counter, defaultdict
from fractions import Fraction

NPASS, FAILS = [0], []


def report(cid, ok, text):
    if ok:
        NPASS[0] += 1
        print(f"PASS z4-{cid}: {text}")
    else:
        FAILS.append(cid)
        print(f"FAIL z4-{cid}: {text}")
    sys.stdout.flush()


# torus coordinates: su(4) in the epsilon basis modulo (1,1,1,1), stored as (w1-w4, w2-w4, w3-w4); su(2)_L, su(2)_R as doubled weights;
# u(1)_F charge.  A weight is a 6-tuple of integers.
def su4(w):
    return (w[0] - w[3], w[1] - w[3], w[2] - w[3])


eps = [tuple(1 if i == j else 0 for i in range(4)) for j in range(4)]
add = lambda a, b: tuple(x + y for x, y in zip(a, b))
sub = lambda a, b: tuple(x - y for x, y in zip(a, b))
w10 = [su4(add(eps[i], eps[j])) for i in range(4) for j in range(i, 4)]          # Sym^2 4
w6 = [su4(add(eps[i], eps[j])) for i in range(4) for j in range(i + 1, 4)]       # Lambda^2 4
w15 = [su4(sub(eps[i], eps[j])) for i in range(4) for j in range(4) if i != j] + [(0, 0, 0)] * 3
w1 = [(0, 0, 0)]
sp1, sp12, sp0 = [-2, 0, 2], [-1, 1], [0]


def irrep(ws4, wl, wr, f):
    return [a + (l, r, f) for a in ws4 for l in wl for r in wr]


V = (irrep(w10, sp1, sp0, 2) + irrep(w10, sp0, sp1, 2) + irrep(w15, sp12, sp12, 0) + 2 * irrep(w6, sp0, sp0, 2)
     + irrep(w1, sp12, sp12, 0))
assert len(V) == 136
Vbar = [tuple(-x for x in w) for w in V]
KC = Counter(V + Vbar)  # character of K (x) C, 272 weights
assert sum(KC.values()) == 272

# roots
roots4 = [su4(sub(eps[i], eps[j])) for i in range(4) for j in range(4) if i != j]
roots = [a + (0, 0, 0) for a in roots4] + [(0, 0, 0, 2, 0, 0), (0, 0, 0, -2, 0, 0), (0, 0, 0, 0, 2, 0), (0, 0, 0, 0, -2, 0)]
W_ORDER = 24 * 2 * 2


def weyl_factor(rs):
    poly = {(0,) * 6: 1}
    for a in rs:
        new = defaultdict(int)
        for k, c in poly.items():
            new[k] += c
            new[add(k, a)] -= c
        poly = {k: c for k, c in new.items() if c}
    return poly


WF = weyl_factor(roots)


def mult(p, q):
    out = defaultdict(int)
    for a, x in p.items():
        for b, y in q.items():
            out[add(a, b)] += x * y
    return dict(out)


def scale(ch, k):
    out = defaultdict(int)
    for w, c in ch.items():
        out[tuple(k * x for x in w)] += c
    return dict(out)


def project_F(ch, keepF):
    """drop the u(1)_F coordinate (set it to 0) unless fermion number is imposed"""
    if keepF:
        return ch
    out = defaultdict(int)
    for w, c in ch.items():
        out[w[:5] + (0,)] += c
    return dict(out)


def ct_with_weyl(poly):
    """constant term of poly * WF, divided by |W|"""
    s = 0
    for w, c in WF.items():
        s += c * poly.get(tuple(-x for x in w), 0)
    return Fraction(s, W_ORDER)


def coeff_product_at(p, q, u):
    s = 0
    for a, x in p.items():
        y = q.get(sub(u, a))
        if y:
            s += x * y
    return s


def ct_product(p, q):
    """CT[(p*q) * WF]/|W| without forming p*q"""
    s = 0
    for w, c in WF.items():
        s += c * coeff_product_at(p, q, tuple(-x for x in w))
    return Fraction(s, W_ORDER)


def counts(keepF):
    ch = project_F(dict(KC), keepF)
    p = {k: scale(ch, k) for k in (1, 2, 3, 4)}
    p11 = mult(p[1], p[1])
    res = {}
    res[1] = ct_with_weyl(p[1])
    # h2 = (p1^2 + p2)/2
    res[2] = (ct_with_weyl(p11) + ct_with_weyl(p[2])) / 2
    # h3 = (p1^3 + 3 p1 p2 + 2 p3)/6
    res[3] = (ct_product(p11, p[1]) + 3 * ct_product(p[1], p[2]) + 2 * ct_with_weyl(p[3])) / 6
    # h4 = (p1^4 + 6 p1^2 p2 + 3 p2^2 + 8 p1 p3 + 6 p4)/24
    res[4] = (ct_product(p11, p11) + 6 * ct_product(p11, p[2]) + 3 * ct_product(p[2], p[2]) + 8 * ct_product(p[1], p[3])
              + 6 * ct_with_weyl(p[4])) / 24
    return res


c0 = counts(False)
report("1", c0[1] == 0 and c0[2] == 18 and c0[3] == 0,
       f"CORE-09: dimensions of the Pati-Salam invariants in Sym^d(K), d = 1, 2, 3: {[str(c0[d]) for d in (1, 2, 3)]} (claimed 0, 18, 0)")
report("2", c0[4] == 686, f"CORE-09: dimension of the Pati-Salam invariants in Sym^4(K): {c0[4]} (claimed 686)")
cF = counts(True)
report("3", cF[1] == 0 and cF[2] == 12 and cF[3] == 0 and cF[4] == 402,
       f"CORE-09: with fermion number imposed, d = 1..4: {[str(cF[d]) for d in (1, 2, 3, 4)]} (claimed 0, 12, 0, 402)")

# tex/s20_selection.tex after thm:potential: the trivial representation occurs in T(x)Tbar, Sym^2 T, Sym^3 T, Sym^2 T (x) Tbar, Sym^4 T,
# Sym^3 T (x) Tbar, Sym^2 T (x) Sym^2 Tbar with multiplicities 1,0,0,0,1,0,5, for T = (10,1,3) (complex, 30-dimensional).
T = Counter(irrep(w10, sp0, sp1, 0))
Tb = Counter(tuple(-x for x in w) for w in T.elements())
assert sum(T.values()) == 30


def h_polys(ch):
    p = {k: scale(dict(ch), k) for k in (1, 2, 3, 4)}
    one = {(0,) * 6: 1}
    h = {0: one, 1: p[1]}
    lin = lambda terms: {k: v for k, v in _sum(terms).items() if v}
    def _sum(terms):
        out = defaultdict(Fraction)
        for coef, poly in terms:
            for k, v in poly.items():
                out[k] += coef * v
        return out
    p11 = mult(p[1], p[1])
    h[2] = lin([(Fraction(1, 2), p11), (Fraction(1, 2), p[2])])
    p111 = mult(p11, p[1])
    h[3] = lin([(Fraction(1, 6), p111), (Fraction(1, 2), mult(p[1], p[2])), (Fraction(1, 3), p[3])])
    h[4] = lin([(Fraction(1, 24), mult(p11, p11)), (Fraction(1, 4), mult(p11, p[2])), (Fraction(1, 8), mult(p[2], p[2])),
                (Fraction(1, 3), mult(p[1], p[3])), (Fraction(1, 4), p[4])])
    return h


hT, hTb = h_polys(T), h_polys(Tb)
pairs = [(1, 1), (2, 0), (3, 0), (2, 1), (4, 0), (3, 1), (2, 2)]
mults = [ct_product(hT[a], hTb[b]) for a, b in pairs]
report("4", [int(m) for m in mults] == [1, 0, 0, 0, 1, 0, 5] and all(m.denominator == 1 for m in mults),
       f"thm:potential: multiplicities of the trivial representation in Sym^a T (x) Sym^b Tbar for (a,b) = {pairs}: "
       f"{[str(m) for m in mults]} (claimed 1,0,0,0,1,0,5)")

print()
print(f"z4: {NPASS[0]} PASS, {len(FAILS)} FAIL" + (f" (failed: {FAILS})" if FAILS else ""))
sys.exit(0 if not FAILS else 1)
