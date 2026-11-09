#!/usr/bin/env python3
"""Boundary cases and counterexamples from the review of version 1.0.1 (Chapters 4 and 9), as exact checks.

    python3 verification/ch09_b_field/a12_review_boundary_cases.py

a12-1  cl:frozen, e5:prop:carry (i): E_0 = |z.z|/2 is conserved but not smooth at nonzero null fields; E_0^2 is a polynomial.
a12-2  e6:rem:compare (b): elements of so(3) of different norms are not conjugate; elements of one nonzero norm are.
a12-3  e4:prop:rotor and the current-algebra remark after it: at rates (1,2,3) the brackets with a != 0, b = 0 give a two-step
       nilpotent frame (lower central series 18, 6, 0), unlike the case a b != 0 (18, 12, 6, 0), the truncated current algebra.

The companion checks of the same review are a3-11 (prop:howmuch), d6-8 (M_pair = |sin phi|) and d7-4c (d7:prop:boost at zero rates).
All arithmetic is exact (sympy). Prints one PASS/FAIL line per check and exits 0 iff all pass.
"""
import sys

import sympy as sp

FAILS = []


def report(cid, ok, text):
    print(f"{'PASS' if ok else 'FAIL'} {cid}: {text}")
    if not ok:
        FAILS.append(cid)


# ---------------------------------------------------------------------------------------------------------- a12-1
t = sp.symbols('t', real=True)
z = sp.Matrix([1, sp.I + t, 0])
q = sp.expand((z.T * z)[0])                                        # z.z = t^2 + 2 i t
E0 = sp.sqrt(sp.re(q) ** 2 + sp.im(q) ** 2) / 2
ok = sp.simplify(q - (t ** 2 + 2 * sp.I * t)) == 0
ok &= sp.simplify(E0 - sp.Abs(t) * sp.sqrt(t ** 2 + 4) / 2) == 0
right = sp.limit(E0 / t, t, 0, '+')
left = sp.limit(E0 / t, t, 0, '-')
ok &= right == 1 and left == -1
E0sq = sp.expand(E0 ** 2)
ok &= sp.Poly(E0sq, t).as_expr() == sp.expand((t ** 4 + 4 * t ** 2) / 4)
report("a12-1", ok, "along z(t) = (1, i+t, 0), z.z = t^2 + 2it and E_0 = |t| sqrt(t^2+4)/2: the one-sided derivatives at the null field "
       "t = 0 are -1 and +1, so E_0 is not smooth there, while E_0^2 = (t^4 + 4t^2)/4 is a polynomial (exact)")

# ---------------------------------------------------------------------------------------------------------- a12-2
def rot(i, j):
    A = sp.zeros(3, 3)
    A[i, j], A[j, i] = -1, 1
    return A


J1, J2 = rot(1, 2), rot(2, 0)                                       # rotations about e_1 and e_2
ok = (J1 * J1).trace() == -2 and ((2 * J1) * (2 * J1)).trace() == -8  # trace form is conjugation invariant
Rm = sp.Matrix([[0, -1, 0], [1, 0, 0], [0, 0, 1]])                   # rotation taking e_1 to e_2
ok &= Rm.T * Rm == sp.eye(3) and Rm.det() == 1 and Rm * J1 * Rm.T == J2
report("a12-2", ok, "so(3): tr(J^2) = -2 and tr((2J)^2) = -8, so J and 2J are not conjugate; a rotation taking e_1 to e_2 conjugates "
       "J(e_1) to J(e_2): the Cartan directions, not all elements, are conjugate (exact)")

# ---------------------------------------------------------------------------------------------------------- a12-3
def cross(x, y):
    return sp.Matrix([x[1] * y[2] - x[2] * y[1], x[2] * y[0] - x[0] * y[2], x[0] * y[1] - x[1] * y[0]])


def bracket(X, Y, a, b):                                           # X = (x1, x2, x3) in (C^3)^3, rates (1,2,3)
    x1, x2 = X[0:3, 0], X[3:6, 0]
    y1, y2 = Y[0:3, 0], Y[3:6, 0]
    return sp.Matrix.vstack(sp.zeros(3, 1), a * cross(x1, y1), b * (cross(x1, y2) + cross(x2, y1)))


def realify(vs):
    return sp.Matrix.hstack(*[sp.Matrix.vstack(v.applyfunc(sp.re), v.applyfunc(sp.im)) for v in vs]) if vs else sp.zeros(18, 0)


def real_basis(cols):                                               # an R-basis of the real span of complex 9-vectors
    if not cols:
        return []
    M = realify(cols)
    piv = M.rref()[1]
    return [cols[k] for k in piv]


BASIS = []
for k in range(9):
    e = sp.zeros(9, 1)
    e[k] = 1
    BASIS += [e, sp.I * e]


def lcs(a, b):
    dims = [18]
    cur = BASIS
    while cur:
        new = real_basis([bracket(X, W, a, b) for X in BASIS for W in cur])
        dims.append(len(new))
        if len(new) == len(cur) and dims[-2] == len(new):
            break
        cur = new
    return dims


d10, d11 = lcs(1, 0), lcs(1, 1)
ok = d10 == [18, 6, 0] and d11 == [18, 12, 6, 0]
report("a12-3", ok, f"rates (1,2,3): the lower central series of the frame has real dimensions {d10} for (a, b) = (1, 0), with V_3 central and "
       f"decoupled, against {d11} for (a, b) = (1, 1), the truncated current algebra; not every resonant coefficient choice is the "
       "standard truncation (exact)")

sys.exit(0 if not FAILS else 1)
