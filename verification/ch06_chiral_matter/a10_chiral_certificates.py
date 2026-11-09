#!/usr/bin/env python3
"""Corrections audit (v1.0.1): exact certificates for Theorem 6.5 (thm:table), Theorem 6.11 (thm:PSalg) and Corollary 6.12 (cor:PSfamily).

    python3 verification/ch06_chiral_matter/a10_chiral_certificates.py

Integer realisation of Appendix A at n = 3 (512 x 512): e1..e6 as in verification/common/clifford8.py, E_(k,a) = om^{(x)(k-1)} (x) e_a (x) 1,
indexed 1..18; interior directions 1, 4, 5, 6; the 14 internal directions the others; weak directions W = {7,13,14,15}, colour directions
C = {10,11,12,16,17,18}. Every matrix below is an integer matrix; Gram matrices Tr(A^T B) are formed in float64 only where every partial sum is an
integer below 2^53 (asserted), and all ranks and null spaces are computed exactly (sympy, rational arithmetic).

a10-1  thm:PSalg     every one of the 91 internal bivector blades commutes or anticommutes with each of the involutions X, Z, U of P_*; exactly 23
                     commute with all three, all 23 compact; the exact null space of B -> [B, P_*] on the 91-dimensional span has dimension 23, so
                     these 23 blades are a basis of it (completeness); they act faithfully on the block (Gram rank 23), close under the bracket,
                     and their algebra has a derived algebra of dimension 21 and a centre of dimension 2.
a10-2  thm:table     for P_fam, P_23, P_13, P_12, P_1: the compact commuting algebra k (exact null space), its kernel on the block, and the centre and
                     derived algebra of k_eff = k / ker: centres 3, 3, 1, 0, 0 and derived algebras 3, 6, 13, 12, 24.
a10-3  cor:PSfamily  the quadratic Casimir of the colour so(6) is the constant 15/4 on the whole module, so it cannot tell 4 from 4bar; the colour
                     volume chi_C commutes with the colour algebra and with J, the cubic invariant sum over the 15 pairings of the six colour
                     directions is a nonzero multiple of chi_C, and chi_C P_* = J P_* on both Weyl blocks: the sign of the colour volume relative to J
                     is the invariant that fixes which of the two four-dimensional representations the block carries.
Each check prints PASS/FAIL a10-<n>; exit status 0 iff all pass.
"""
import itertools
import os
import sys
from functools import reduce

import numpy as np
import sympy as sp

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "common"))
from clifford8 import e as E8  # noqa: E402

NPASS, FAILS = [0], []


def report(cid, ok, text):
    if ok:
        NPASS[0] += 1
        print(f"PASS a10-{cid}: {text}")
    else:
        FAILS.append(cid)
        print(f"FAIL a10-{cid}: {text}")
    sys.stdout.flush()


e8 = [None] + [np.array(E8[a].tolist(), dtype=np.int64) for a in range(1, 7)]
I8 = np.eye(8, dtype=np.int64)
om = reduce(np.matmul, e8[1:])
kr = lambda *a: reduce(np.kron, a)  # noqa: E731
E = [None] + [kr(*([om] * k + [e8[a]] + [I8] * (2 - k))).astype(np.int64) for k in range(3) for a in range(1, 7)]
I512 = np.eye(512, dtype=np.int64)
Z512 = np.zeros((512, 512), dtype=np.int64)
blq = lambda S: reduce(np.matmul, [E[i] for i in S], I512)  # noqa: E731
eta18 = {A: (1 if (A - 1) % 6 < 3 else -1) for A in range(1, 19)}
interior = [1, 4, 5, 6]
internal = [i for i in range(1, 19) if i not in interior]
Jn = blq(internal)
Gam = Jn @ blq(interior)
IB = list(itertools.combinations(internal, 2))
Bv = {ab: E[ab[0]] @ E[ab[1]] for ab in IB}
compact = {ab: eta18[ab[0]] == eta18[ab[1]] for ab in IB}


def comm(A, B):
    return A @ B - B @ A


def gram(mats):
    """exact integer Gram matrix Tr(A^T B) of integer matrices (float64 BLAS; every partial sum is an integer below 2^53)."""
    flat = np.array([m.reshape(-1) for m in mats], dtype=np.float64)
    g = flat @ flat.T
    assert np.abs(g).max() < 2 ** 52 and np.array_equal(g, np.round(g))
    return sp.Matrix(np.round(g).astype(np.int64))


# ------------------------------------------------------------------ a10-1: thm:PSalg
chiW, chiC, qA = blq([7, 13, 14, 15]), blq([10, 11, 12, 16, 17, 18]), E[2] @ E[3]
X, Z, U = Gam @ chiW, -(Jn @ chiC), -(Jn @ qA)
P8 = (I512 + X) @ (I512 + Z) @ (I512 + U)                                   # 8 P_*
assert np.array_equal(P8 @ P8, 8 * P8) and np.trace(P8) == 8 * 64
pattern, ok_cac = {}, True
for ab in IB:
    sig = []
    for Y in (X, Z, U):
        c, a = comm(Bv[ab], Y), Bv[ab] @ Y + Y @ Bv[ab]
        if not c.any():
            sig.append(0)
        elif not a.any():
            sig.append(1)
        else:
            ok_cac = False
            sig.append(None)
    pattern[ab] = tuple(sig)
basis = [ab for ab in IB if pattern[ab] == (0, 0, 0)]
G_all = gram([comm(Bv[ab], P8) for ab in IB])
null_dim = len(IB) - G_all.rank()
ok1 = ok_cac and len(basis) == 23 and null_dim == 23 and all(compact[ab] for ab in basis)
ok1 &= all(np.array_equal(Bv[ab].T, -Bv[ab]) for ab in basis)
faithful = gram([Bv[ab] @ P8 for ab in basis]).rank()
ok1 &= faithful == 23
# structure constants: [B_i, B_j] is 0 or +-2 B_k with B_k in the basis (closure), computed exactly
idx = {ab: i for i, ab in enumerate(basis)}
cols, closed = [], True
for ab, cd in itertools.combinations(basis, 2):
    C = comm(Bv[ab], Bv[cd])
    if not C.any():
        continue
    hit = None
    for ef in basis:
        for s in (2, -2):
            if np.array_equal(C, s * Bv[ef]):
                hit = (ef, s)
    if hit is None:
        closed = False
        continue
    v = [0] * 23
    v[idx[hit[0]]] = hit[1]
    cols.append(v)
derived = sp.Matrix(cols).rank() if cols else 0
# centre: x = sum c_i B_i with [x, B_j] = 0 for all j; ad matrices from the structure constants
ad_rows = []
for j, cd in enumerate(basis):
    M = sp.zeros(23, 23)
    for i, ab in enumerate(basis):
        C = comm(Bv[ab], Bv[cd])
        if C.any():
            for ef in basis:
                for s in (2, -2):
                    if np.array_equal(C, s * Bv[ef]):
                        M[idx[ef], i] = s
    ad_rows.append(M)
centre = 23 - sp.Matrix.vstack(*ad_rows).rank()
ok1 &= closed and derived == 21 and centre == 2
report(1, ok1, f"every internal bivector blade commutes or anticommutes with X, Z, U ({ok_cac}); {len(basis)} commute with all three, all compact "
               f"and antisymmetric; the exact null space of B -> [B,P_*] on the 91 blades has dimension {null_dim}, so they span it; they act "
               f"faithfully on the block (rank {faithful}), close under the bracket ({closed}), with derived algebra of dimension {derived} and centre "
               f"of dimension {centre} (exact)")

# ------------------------------------------------------------------ a10-2: thm:table
J = {k: E[6 * k + 2] @ E[6 * k + 3] for k in range(3)}
omk = {k: blq(range(6 * k + 1, 6 * k + 7)) for k in range(3)}
X1, X2, X3 = J[0] @ J[1], omk[2], omk[0] @ (E[13] @ E[16])
Pr = {"Pfam": (I512 + X1) @ (I512 + X2) @ (I512 + X3), "P23": 2 * (I512 + X2) @ (I512 + X3), "P13": 2 * (I512 + X1) @ (I512 + X3),
      "P12": 2 * (I512 + X1) @ (I512 + X2), "P1": 4 * (I512 + X1)}         # all = 8 P
comp = [ab for ab in IB if compact[ab]]                                   # the 43 compact internal bivector blades
Bf = {ab: Bv[ab].astype(np.float64) for ab in comp}
cidx = {ab: i for i, ab in enumerate(comp)}
# exact structure constants of the compact bivector algebra on the blade basis: [B_a, B_b] = 0 or +-2 B_c
SC = {}
for ab, cd in itertools.combinations(comp, 2):
    C = comm(Bv[ab], Bv[cd])
    if not C.any():
        continue
    hit = [(ef, s_) for ef in comp for s_ in (2, -2) if np.array_equal(C, s_ * Bv[ef])]
    assert len(hit) == 1
    SC[(ab, cd)] = hit[0]


def bracket(u, v):
    """bracket of two coefficient vectors (sympy column vectors of length 43) in the compact blade basis"""
    out = sp.zeros(len(comp), 1)
    for (ab, cd), (ef, s_) in SC.items():
        i, j, k_ = cidx[ab], cidx[cd], cidx[ef]
        c = u[i] * v[j] - u[j] * v[i]
        if c != 0:
            out[k_] += s_ * c
    return out


result = {}
for name, Pi in Pr.items():
    Pf = Pi.astype(np.float64)
    Gc = gram([Bf[ab] @ Pf - Pf @ Bf[ab] for ab in comp])
    kb = Gc.nullspace()                                                    # basis of k in coefficient space (exact)
    k = len(kb)
    mats_k = []
    for v in kb:
        den = sp.ilcm(*[sp.fraction(x)[1] for x in v])
        mats_k.append(sum(int(x * den) * Bv[ab] for x, ab in zip(v, comp) if x != 0))
    # kernel on the block: coefficient vectors of k annihilating the block
    Gk = gram([m @ Pi for m in mats_k])
    kerc = Gk.nullspace()
    kerdim = len(kerc)
    Kmat = sp.Matrix.hstack(*kb)                                           # 43 x k
    ker_vecs = [Kmat * c for c in kerc]
    # derived algebra of k_eff: dim([k,k] + ker) - dim ker
    brs = [bracket(kb[i], kb[j]) for i, j in itertools.combinations(range(k), 2)]
    span = brs + ker_vecs
    derived = (sp.Matrix.hstack(*span).rank() if span else 0) - kerdim
    # centre of k_eff: {x in k : [x, k] within ker} / ker
    if ker_vecs:
        Nperp = sp.Matrix.hstack(*ker_vecs).T.nullspace()                  # basis of the orthogonal complement of ker in R^43
        Nt = sp.Matrix.hstack(*Nperp).T
    else:
        Nt = sp.eye(len(comp))
    rows = []
    for j in range(k):
        cols_ = [Nt * bracket(kb[i], kb[j]) for i in range(k)]
        rows.append(sp.Matrix.hstack(*cols_))
    Mcen = sp.Matrix.vstack(*rows)
    centre = (k - Mcen.rank()) - kerdim
    result[name] = (k, kerdim, centre, derived)
expected = {"Pfam": (9, 3, 3, 3), "P23": (9, 0, 3, 6), "P13": (17, 3, 1, 13), "P12": (15, 3, 0, 12), "P1": (27, 3, 0, 24)}
ok2 = result == expected
report(2, ok2, "dim k, dim ker, centre and derived algebra of k_eff for (P_fam, P_23, P_13, P_12, P_1) = "
               f"{[result[k] for k in expected]}: centres 3,3,1,0,0 and derived algebras 3,6,13,12,24 (exact null spaces and ranks)")

# ------------------------------------------------------------------ a10-3: cor:PSfamily, colour volume vs quadratic Casimir
Cdir = [10, 11, 12, 16, 17, 18]
col = {ab: E[ab[0]] @ E[ab[1]] for ab in itertools.combinations(Cdir, 2)}
quad = sum(c @ c for c in col.values())                                    # sum of B_ab^2 = -15 * 1, so -(1/4) sum = 15/4
ok3 = np.array_equal(quad, -15 * I512)
ok3 &= all(not comm(chiC, c).any() for c in col.values()) and not comm(chiC, Jn).any() and np.array_equal(chiC @ chiC, -I512)
Cf = {ab: c.astype(np.float64) for ab, c in col.items()}


def matchings(s):
    if not s:
        yield []
        return
    a = s[0]
    for i in range(1, len(s)):
        b = s[i]
        rest = s[1:i] + s[i + 1:]
        for m in matchings(rest):
            yield [(a, b)] + m


def perm_sign(seq):
    seq, sgn = list(seq), 1
    for i in range(len(seq)):
        for j in range(i + 1, len(seq)):
            if seq[i] > seq[j]:
                sgn = -sgn
    return sgn


cubic = np.zeros((512, 512))
nterms = 0
for m in matchings(Cdir):
    for order in itertools.permutations(m):
        flat = [x for pair in order for x in pair]
        cubic += perm_sign(flat) * (Cf[order[0]] @ Cf[order[1]] @ Cf[order[2]])
        nterms += 1
ratio = np.trace(cubic.T @ chiC.astype(np.float64)) / np.trace(chiC.T @ chiC)
ok3 &= nterms == 90 and np.array_equal(cubic, np.round(ratio) * chiC.astype(np.float64)) and round(ratio) != 0
PL, PR = P8 @ (I512 - Gam), P8 @ (I512 + Gam)                             # 16 P_* (1 -+ Gamma)/2
ok3 &= np.array_equal(chiC @ P8, Jn @ P8) and np.array_equal(chiC @ PL, Jn @ PL) and np.array_equal(chiC @ PR, Jn @ PR)
ok3 &= np.trace(PL) == 16 * 32 and np.trace(PR) == 16 * 32
report(3, ok3, f"sum of the squared colour bivectors = -15, so the colour Casimir is 15/4 on the whole module and cannot tell 4 from 4bar; chi_C^2 = -1 "
               f"commutes with the 15 colour generators and with J; the cubic invariant over the {nterms} ordered pairings equals {int(round(ratio))} chi_C; "
               "chi_C P_* = J P_* on both Weyl blocks (rank 32 each): on the block the cubic invariant is a fixed multiple of J, and its sign, the sign of "
               "the colour volume relative to J, fixes the four-dimensional representation (exact integer arithmetic)")

sys.exit(0 if not FAILS and NPASS[0] > 0 else 1)
