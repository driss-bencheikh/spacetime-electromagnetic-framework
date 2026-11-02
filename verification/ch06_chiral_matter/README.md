# Checks for the global manuscript (`tex/`)

These scripts were written in Phase 3 (1 November 2026). They adjudicate the CORE review items and check one point of the literature
audit. They do not import any other verification script. Each prints one `PASS` or `FAIL` line per check and exits 0 iff all pass.

```bash
python3 verification/core/z3_core.py
```

```bash
python3 verification/core/z4_invariants.py
```

| script | checks | what it tests | review items |
|---|---|---|---|
| `z3_core.py` | z3-1 … z3-9 (about 1 min) | blade commutation (`ex:blades`); trace potential on the kernel space (`lem:tracepot`, `thm:twophases`); orbit dimensions; vector-mass normalisation (`thm:compat`, eq. (45)); cross-copy intertwiners (`thm:closureunique`); complex structure on the kernel space (`thm:kernelspace`); which bivectors of so(4,6) preserve $P_*$ (Krasnov comparison) | CORE-01, 02, 03, 05, 06; literature |
| `z4_invariants.py` | z4-1 … z4-4 (about 20 s) | Molien–Weyl counts of Pati–Salam-invariant polynomials on the kernel space (degrees 1–4, with and without fermion number) and the multiplicities of the tensor model, in exact integer arithmetic | CORE-09 |

**Requirements.** Python 3 with numpy and sympy.

**Not covered.** The 271-dimensional angular Hessian and the minimisation runs of §24 (CORE-09, CORE-10) are not reproduced here. They are
planned for the Phase 9 verification package.

The outputs of the Phase 3 runs are in `reproducibility/phase3_runs/`.
