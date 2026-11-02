#!/usr/bin/env python3
"""v3: random-start minimisation of the vacuum energies and the crossing x_* (Phase 9B, 2 November 2026).

    python3 verification/vacuum/v3_minimisation.py

Reproduces, with archived seeds, start distribution, optimisers and tolerances, the numerical statements of Section 7.3 of the manuscript:
the gradient flows after thm:selectk (tree potential with the holomorphic quartic at the benchmark), the projected flows after thm:oneloop1013
(one loop on the tensor model), and prop:randomstart (one loop on the 272-dimensional kernel space), including the crossing value x_*.

Method. Each run minimises F(c/|c|) over c in R^d from a standard Gaussian start (numpy default_rng with the printed seed) by scipy's L-BFGS-B
(maxcor 30, gtol 1e-12, ftol 1e-16), and at x = 25 also by scipy's nonlinear conjugate gradient; the residual is the norm of the Riemannian
gradient at the end point. Endpoints are grouped by energy (rounded to 1e-6). These are numerical experiments: they support, and do not prove,
statements about global minima. A local Hessian test (v2) does not establish global minimality either.
"""
import collections
import sys
import time

import numpy as np
from scipy.optimize import brentq, minimize

import kernel_space as k
import tensor_model as tm

NPASS, FAILS = [0], []
T0 = time.time()


def report(cid, ok, text):
    if ok:
        NPASS[0] += 1
        print(f"PASS v3-{cid}: {text}")
    else:
        FAILS.append(cid)
        print(f"FAIL v3-{cid}: {text}")
    sys.stdout.flush()


LB = dict(maxiter=50000, gtol=1e-12, ftol=1e-16, maxcor=30)

# ======================================================================== tensor model, tree potential with kappa (after thm:selectk)
lam = np.array([1, 1.7, 2.3, 3.1, 2.9])
out = {}
for kap in (6.0, 7.0, 7.5, 8.0):
    ends = []
    for r in range(16):
        x0 = np.random.default_rng(100 + r).normal(size=60)
        res = minimize(tm.energy, x0, args=(lam, kap), jac=True, method="L-BFGS-B", options=LB)
        ends.append(res.fun)
    out[kap] = ends
okk = True
desc = []
for kap, ends in out.items():
    eS0, e11 = 1.0, 1.9 - kap / 8
    lo = min(eS0, e11)
    found = sorted({round(e, 6) for e in ends})
    okk &= min(ends) > lo - 1e-8 and abs(min(ends) - lo) < 1e-6
    desc.append(f"|kappa| = {kap}: end energies {found} (S_0 orbit {eS0}, 1(x)1 orbit {e11:.4f})")
okk &= (1.9 - 7.0 / 8 > 1.0) and (1.9 - 7.5 / 8 < 1.0)
report("1", okk, "tree potential at the benchmark lambda = (1,1.7,2.3,3.1,2.9) on the unit sphere of T, 16 starts per |kappa| (seeds 100-115): "
       + "; ".join(desc) + ". No end point lies below min(Lambda_kappa(S_0), Lambda_kappa(1(x)1)); the lower of the two changes at |kappa| = 7.2 "
       "(numerical; the claim 'one of two starts ends on 1(x)1 for 6 <= |kappa| <= 7.5' depends on the start distribution and is not tested)")

# ======================================================================== tensor model, one loop (after thm:oneloop1013)
desc, okk = [], True
for x in (0.0, 20.0, 27.0, 28.0, 29.0):
    ell = 2 * x * tm.a_rec - 3 * tm.c_vec
    ends = []
    for r in range(16):
        x0 = np.random.default_rng(200 + r).normal(size=60)
        ends.append(minimize(tm.energy, x0, args=(ell,), jac=True, method="L-BFGS-B", options=LB).fun)
    sm, e11 = 2 * x - 99, 2 * x - 99 + 48 - 7 * x / 4
    found = sorted({round(e, 6) for e in ends})
    if x < 192 / 7:
        okk &= min(ends) > sm - 1e-8
    else:
        okk &= abs(min(ends) - e11) < 1e-6 and e11 < sm
    desc.append(f"x = {x:g}: {found} (SM {sm:g}, 1(x)1 {e11:.4f})")
report("2", okk, "one-loop energy sum_k(2x a_k - 3c_k) f_k on the unit sphere of T, 16 starts per x (seeds 200-215): " + "; ".join(desc)
       + ". For x = 0, 20, 27 nothing lies below the Standard-Model orbit; for x = 28, 29 (> 192/7) the lowest end point is the 1(x)1 orbit "
       "(numerical evidence for the conjecture on T; not a proof)")

# ======================================================================== the kernel space
def fg(c, x):
    nrm = np.linalg.norm(c); u = c / nrm
    E = k.Ex(u, x); g = k.grad(u, x); g = g - (u @ g) * u
    return E, g / nrm


def run(x, seed, method="L-BFGS-B"):
    c = np.random.default_rng(seed).normal(size=k.dimK)
    opts = LB if method == "L-BFGS-B" else dict(maxiter=200000, gtol=1e-10)
    r = minimize(fg, c, args=(x,), jac=True, method=method, options=opts)
    u = r.x / np.linalg.norm(r.x)
    return r.fun, u, float(np.linalg.norm(fg(u, x)[1]))


# components (as in v2) for the weights of end points
A = k.A
CL = sum(A[a] @ A[a] for a in range(0, 3)); CR = sum(A[a] @ A[a] for a in range(3, 6)); C4 = sum(A[a] @ A[a] for a in range(6, 21))
ev, U = np.linalg.eigh(C4 + np.pi * CL + np.e * CR)
comp = {}
for val in np.unique(np.round(ev, 6)):
    cols = U[:, np.abs(ev - val) < 1e-5]
    lab = tuple(int(round(float(np.trace(cols.T @ C @ cols)) / cols.shape[1] * 4)) for C in (C4, CL, CR))
    comp.setdefault(lab, []).append(cols)
names = {(-144, -8, 0): "(10,3,1)", (-144, 0, -8): "(10,1,3)", (-128, -3, -3): "(15,2,2)", (-80, 0, 0): "2(6,1,1)", (0, -3, -3): "(1,2,2)"}
comp = {names[l]: np.hstack(v) for l, v in comp.items()}


def weights(u):
    return {n: float(np.sum((b.T @ u) ** 2)) for n, b in comp.items()}


def classify(E, u, x):
    if abs(E - (2 * x - 99)) < 1e-6:
        wt = weights(u)
        return "SM" if wt["(10,1,3)"] > 1 - 1e-6 else ("mirror" if wt["(10,3,1)"] > 1 - 1e-6 else "SM-energy other")
    return f"E={E:.6f}"


tally, resid_max, keep = {}, 0.0, {}
mirror_forms = collections.defaultdict(collections.Counter)
plan = {0.0: 40, 20.0: 40, 22.0: 40, 24.5: 40, 25.0: 100, 26.0: 40}
seed0 = {0.0: 1000, 20.0: 2000, 22.0: 3000, 24.5: 4000, 25.0: 5000, 26.0: 6000}
for x, nrun in plan.items():
    cnt = collections.Counter()
    for r in range(nrun):
        E, u, res = run(x, seed0[x] + r)
        resid_max = max(resid_max, res)
        lab = classify(E, u, x)
        cnt[lab] += 1
        if lab.startswith("E=") and E < 2 * x - 99 - 1e-6:
            wt_ = weights(u)
            form = "L" if wt_["(10,1,3)"] >= wt_["(10,3,1)"] else "R"
            mirror_forms[x][form] += 1
            if form == "L" and (x not in keep or E < keep[x][0] - 1e-9):
                keep[x] = (E, u)
    tally[x] = cnt
summ = "; ".join(f"x = {x:g} ({plan[x]} runs, seeds {seed0[x]}-{seed0[x] + plan[x] - 1}): {dict(c)}" for x, c in tally.items())
# Corrected in Phase 10: the robust statement is that no run ends below the Standard-Model energy for x = 0, 20, 22. Which end points
# above it are reached depends on the floating-point path (one run of 40 at x = 22 stops at a higher critical value when BLAS runs on
# one thread); the tally is reported, not asserted.
def below_sm(lab, x):
    return lab.startswith("E=") and float(lab[2:]) < 2 * x - 99 - 1e-6


ok = all(not any(below_sm(k, x) for k in tally[x]) for x in (0.0, 20.0, 22.0)) and resid_max < 1e-5
ok &= all(tally[x]["SM"] + tally[x]["mirror"] >= 0.9 * plan[x] for x in (0.0, 20.0, 22.0))
report("3", ok, f"prop:randomstart, x = 0, 20, 22: no run ends below the energy 2x - 99 of the Standard-Model orbit, and at least 90% end at that "
       f"orbit or its mirror (end points classified by energy and by their weight on (10,1,3) or (10,3,1)); largest final Riemannian gradient "
       f"{resid_max:.1e}. All runs: {summ} (numerical)")

E25, u25 = keep[25.0]
frac = sum(v for kk, v in tally[25.0].items() if kk.startswith("E=") and float(kk[2:]) < 2 * 25 - 99 - 1e-6) / plan[25.0]
wt = weights(u25)
Z = np.einsum("aji,i->aj", A, u25)
sv = np.linalg.svd(Z, compute_uv=False)
rank = int(np.sum(sv > 1e-6 * sv.max()))
stab = np.linalg.svd(Z.T)[2][rank:] if rank < 21 else np.zeros((0, 21))
stab = np.linalg.svd(Z.T, full_matrices=True)[2][rank:]
contains_L = np.linalg.norm(Z[0:3]) < 1e-6
colour_rank = int(np.linalg.matrix_rank(stab[:, 6:21], tol=1e-6)) if len(stab) else 0
um = k.mat(u25)
s_u = np.linalg.svd(um, compute_uv=False)
nz = s_u[s_u > 1e-6]
levels = sorted(collections.Counter(np.round(nz / nz.min(), 4)).items())
ok = abs(E25 - (-49.13987)) < 5e-6 and abs(wt["(10,1,3)"] - 0.968) < 5e-4 and abs(wt["2(6,1,1)"] - 0.032) < 5e-4
ok &= len(stab) == 9 and contains_L and colour_rank == 6 and 0.15 < frac < 0.5
report("4", ok, f"prop:randomstart, x = 25: {frac:.0%} of 100 runs end below the Standard-Model energy -49, at E = {E25:.7f}, in two forms exchanged by the involution of prop:LR ({mirror_forms[25.0]['L']} with su(2)_L unbroken, {mirror_forms[25.0]['R']} mirror images); the form with su(2)_L unbroken has weight "
       f"{wt['(10,1,3)']:.4f} on (10,1,3) and {wt['2(6,1,1)']:.4f} on 2(6,1,1), a stabiliser of dimension {len(stab)} among the 21 Pati-Salam generators "
       f"that contains su(2)_L and projects onto colour with rank {colour_rank} (< 8, so colour is broken); singular values of the kernel, relative to "
       f"the smallest, {levels} (multiplicity: real components). The manuscript's 'a third of the runs' is {frac:.0%} with this start distribution (numerical)")

# a second optimiser
cg = [run(25.0, 5000 + r, method="CG") for r in range(40)]
cg_low = min(e for e, _, _ in cg)
ok = abs(cg_low - E25) < 1e-7
report("5", ok, f"a second optimiser (nonlinear conjugate gradient, 40 of the same starts) reaches the same lowest end point at x = 25, E = {cg_low:.9f} "
       f"(difference {abs(cg_low - E25):.1e}) (numerical)")

# local minimality of the competitor and the crossing x_*
def refine(x, u):
    r = minimize(fg, u, args=(x,), jac=True, method="L-BFGS-B", options=LB)
    v = r.x / np.linalg.norm(r.x)
    return r.fun, v


state = {"u": u25}


def delta(x):
    E, v = refine(x, state["u"]); state["u"] = v
    return E - (2 * x - 99)


d25, d248 = delta(25.0), delta(24.80)
state["u"] = u25
xs = brentq(delta, 24.80, 25.0, xtol=1e-10, rtol=1e-12)
E_star, u_star = refine(xs, state["u"])
Hf = 8 * xs * k.hessQ(u_star) - 3 / 16 * k.hessT(u_star)
gf = k.grad(u_star, xs)
Hs, _ = k.sphere_hessian(Hf, gf, u_star)
eh = np.sort(np.linalg.eigvalsh(Hs))
nzero = int(np.sum(np.abs(eh) < 1e-6))
ok = abs(xs - 24.88405) < 5e-6 and eh[0] > -1e-6 and d248 > 0 > d25
report("6", ok, f"crossing: following the competitor in x (warm-started L-BFGS-B) and solving E_comp(x) = 2x - 99 by Brent's method (xtol 1e-10) gives "
       f"x_* = {xs:.7f}, against the quoted 24.88405 (difference {abs(xs - 24.88405):.1e}); E_comp - E_SM = {d248:+.6f} at x = 24.80 and {d25:+.6f} at "
       f"x = 25; at x_* the competitor is a local minimum on the sphere (lowest Riemannian Hessian eigenvalue {eh[0]:.1e}, {nzero} zero modes) (numerical)")

nb = sum(v for kk, v in tally[24.5].items() if kk.startswith("E=") and float(kk[2:]) < 2 * 24.5 - 99 - 1e-6)
ok = nb == 0 and any(kk.startswith("E=") for kk in tally[26.0])
report("7", ok, f"either side of x_*: at x = 24.5 no run of 40 ends below the Standard-Model energy ({dict(tally[24.5])}); at x = 26 the competitor branch is found "
       f"({dict(tally[26.0])}). This is evidence for, not a proof of, the conjecture that the Standard-Model orbit and its mirror are the global minima "
       "exactly for x < x_* (numerical)")

print(f"{NPASS[0]}/{NPASS[0] + len(FAILS)} checks passed ({time.time() - T0:.0f} s)")
sys.exit(0 if not FAILS else 1)
