#!/usr/bin/env python3
"""Corrections audit (v1.0.1): exact certificates for the two restricted identities in the proof of Theorem 5.13 (thm:blockgamma).

    python3 verification/ch05_clifford_tower/a9_blockgamma_certificates.py

Conventions (Appendix A): e1..e6 the 8x8 realisation of verification/common/clifford8.py, omega = e1...e6; at n = 3 constituents
E_(k,a) = omega^{(x)(k-1)} (x) e_a (x) 1^{(x)(3-k)}; interior gammas gamma^(k)_mu = E_(k,a) with a = 1,4,5,6 for mu = 0,1,2,3,
eta_mu nu = diag(1,-1,-1,-1); the diagonal Lorentz generators Delta_ab = (1/2) sum_k gamma^(k)_a gamma^(k)_b (a != b); the Lorentz Casimir
C = sum_{a<b} eta_aa eta_bb Delta_ab^2 (any normalisation: only its eigenspaces are used); the Dirac block D is its eigenspace of
dimension 160; omega_S = E_(1,2)E_(1,3)E_(2,2)E_(2,3)E_(3,2)E_(3,3), the product of the six surplus generators;
Q_mu = sum_nu eta_nu nu (gamma^(1)_mu gamma^(2)_nu gamma^(3)_nu + gamma^(1)_nu gamma^(2)_mu gamma^(3)_nu + gamma^(1)_nu gamma^(2)_nu gamma^(3)_mu).
All matrices are integer matrices (4C and 2 Delta included); the projector onto D is P = prod_{l != l_D}(C - l)/(l_D - l), and the
script works with the integer multiple P' = prod_{l != l_D}(4C - 4l), so every identity below is an exact integer identity.

a9-1  thm:blockgamma   C has exactly three eigenvalues, with multiplicities 64, 288 and 160; Q_mu commutes with C, so it preserves D;
                       and ({Q_mu, Q_nu} + 72 eta_mu nu) P' = 0 for all mu, nu: {Q_mu,Q_nu}|_D = -72 eta_mu nu, exactly.
a9-2  thm:blockgamma   omega_S^2 = -1, [omega_S, Q_mu] = 0 and [omega_S, C] = 0; with G_mu = (1/6) omega_S Q_mu|_D this gives
                       [G_a, G_b] = -(1/36)[Q_a, Q_b] on D, and ([Q_a, Q_b] + 144 Delta_ab) P' = 0 for all a < b:
                       (1/4)[G_a, G_b] = Delta_ab|_D, exactly.
Each check prints PASS/FAIL a9-<n>; exit status 0 iff all pass.
"""
import os
import sys
from functools import reduce
from itertools import combinations

import numpy as np
import sympy as sp

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "common"))
from clifford8 import e as E8  # noqa: E402

NPASS, FAILS = [0], []


def report(cid, ok, text):
    if ok:
        NPASS[0] += 1
        print(f"PASS a9-{cid}: {text}")
    else:
        FAILS.append(cid)
        print(f"FAIL a9-{cid}: {text}")
    sys.stdout.flush()


def to_int(M):
    return np.array(M.tolist(), dtype=np.int64)


e1 = [None] + [to_int(E8[a]) for a in range(1, 7)]
I8 = np.eye(8, dtype=np.int64)
omega = reduce(np.matmul, [e1[a] for a in range(1, 7)])
assert np.array_equal(omega @ omega, I8)


def kron(*ms):
    return reduce(np.kron, ms)


def E(k, a, n=3):
    return kron(*([omega] * (k - 1) + [e1[a]] + [I8] * (n - k)))


N = 512
ID = np.eye(N, dtype=np.int64)
MU = [1, 4, 5, 6]
ETA = [1, -1, -1, -1]
gam = {(k, m): E(k, MU[m]) for k in (1, 2, 3) for m in range(4)}
for k in (1, 2, 3):
    for m in range(4):
        for n in range(4):
            ac = gam[k, m] @ gam[k, n] + gam[k, n] @ gam[k, m]
            assert np.array_equal(ac, 2 * ETA[m] * ID if m == n else 0 * ID)

# 2 Delta_ab (integer) and 4C (integer)
D2 = {(a, b): sum(gam[k, a] @ gam[k, b] for k in (1, 2, 3)) for a, b in combinations(range(4), 2)}
C4 = sum(ETA[a] * ETA[b] * (D2[a, b] @ D2[a, b]) for a, b in combinations(range(4), 2))


def comm(A, B):
    return A @ B - B @ A


def maxabs(A):
    return int(np.abs(A).max())


# eigenvalues of 4C: exact characteristic data via sympy on the small set of distinct values found from the trace powers
vals = np.linalg.eigvalsh(C4.astype(float))
cand = sorted(set(int(round(v)) for v in vals))
mins = {}
for lam in cand:
    mins[lam] = int(np.sum(np.abs(vals - lam) < 1e-6))
ok_spec = (len(cand) == 3 and sorted(mins.values()) == [64, 160, 288])
# exact confirmation: the product of (4C - lambda) over the three candidates vanishes, and the eigenspace dimensions are the traces of the
# exact projectors
prod_all = reduce(np.matmul, [C4 - l * ID for l in cand])
ok_spec &= not prod_all.any()
lamD = [l for l in cand if mins[l] == 160][0]
others = [l for l in cand if l != lamD]
Pp = reduce(np.matmul, [C4 - l * ID for l in others])  # = d * P_D with d below
d = reduce(lambda x, y: x * y, [lamD - l for l in others])
ok_spec &= (int(np.trace(Pp)) == 160 * d) and np.array_equal(Pp @ Pp, d * Pp)
assert maxabs(Pp) < 10 ** 9

# Q_mu
Q = []
for m in range(4):
    acc = np.zeros((N, N), dtype=np.int64)
    for n in range(4):
        acc += ETA[n] * (gam[1, m] @ gam[2, n] @ gam[3, n] + gam[1, n] @ gam[2, m] @ gam[3, n] + gam[1, n] @ gam[2, n] @ gam[3, m])
    Q.append(acc)

ok1 = ok_spec and all(not comm(Q[m], C4).any() for m in range(4))
bad = []
for m in range(4):
    for n in range(m, 4):
        lhs = Q[m] @ Q[n] + Q[n] @ Q[m] + (72 * ETA[m] * ID if m == n else 0 * ID)
        assert maxabs(lhs) * maxabs(Pp) * N < 2 ** 62
        if (lhs @ Pp).any():
            bad.append((m, n))
ok1 &= not bad
report(1, ok1, f"4C has eigenvalues {cand} with multiplicities {[mins[l] for l in cand]}; the Dirac block D is the eigenspace "
              f"4C = {lamD} of dimension 160 (trace of the exact projector); [Q_mu, C] = 0 and ({{Q_mu,Q_nu}} + 72 eta_mu nu) P_D = 0 for all "
              f"mu, nu: {{Q_mu,Q_nu}}|_D = -72 eta_mu nu (exact integer arithmetic){'' if not bad else ' -- failing pairs ' + str(bad)}")

omS = E(1, 2) @ E(1, 3) @ E(2, 2) @ E(2, 3) @ E(3, 2) @ E(3, 3)
ok2 = np.array_equal(omS @ omS, -ID) and not comm(omS, C4).any() and all(not comm(omS, Q[m]).any() for m in range(4))
bad2 = []
for a, b in combinations(range(4), 2):
    lhs = comm(Q[a], Q[b]) + 72 * D2[a, b]          # 144 Delta_ab = 72 * (2 Delta_ab)
    assert maxabs(lhs) * maxabs(Pp) * N < 2 ** 62
    if (lhs @ Pp).any():
        bad2.append((a, b))
ok2 &= not bad2
report(2, ok2, "omega_S^2 = -1, [omega_S, Q_mu] = 0 and [omega_S, C] = 0, so G_mu = (1/6) omega_S Q_mu|_D satisfies [G_a, G_b] = -(1/36)[Q_a, Q_b] "
              "on D; ([Q_a, Q_b] + 144 Delta_ab) P_D = 0 for all a < b: (1/4)[G_a, G_b] = Delta_ab|_D (exact integer arithmetic)"
              + ('' if not bad2 else ' -- failing pairs ' + str(bad2)))

sys.exit(0 if not FAILS and NPASS[0] > 0 else 1)
