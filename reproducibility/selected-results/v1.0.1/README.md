# Regression of version 1.0.1 (8 November 2026)

**Command** (from the repository root; the scripts run are byte-identical to the recorded revision, as `environment.txt` records no
modification of `verification/`):

```bash
python3 verification/run_all.py --out reproducibility/selected-results/v1.0.1 -j 3
```

**Environment** (`environment.txt`): Python 3.13.18, numpy 2.5.3, scipy 1.18.1, sympy 1.14.0, mpmath 1.3.0; Linux x86_64 (glibc 2.39);
3 workers; one BLAS thread per script. The longest script, `ch06_chiral_matter/a10_chiral_certificates.py`, takes 23 min. In this run one
process of `ch03_readings_of_x/x1_unitary_energy_referee.py` (normally 17 s) hung after 19 of its 28 checks and was stopped after 48 min; the
script was then rerun on its own through the same `run_all.run` routine and environment, passed all 28 checks in 16.8 s, and its record and
summary row are that rerun. The scripts are unchanged.

**Result** (`summary.tsv`): **66 scripts, 1405 PASS, 0 FAIL**; every exit status is zero and consistent with the PASS/FAIL lines. One script
(`ch08_gauge_photon/p3_identification.py`) writes its timing lines to standard error. Against the run record of version 1.0.0
(`../v1.0.0/`), every check identifier of that record is reproduced with the same outcome; sixteen identifiers are new (a9-1, a9-2, a10-1
to a10-3, a11-1 to a11-3, a12-1 to a12-3, a13-1, a13-2, a3-11, c1-12, d7-4c), and the claim texts differ only where version 1.0.1 changed a check (a3-9,
a3-10, c1-4, c1-5, c1-8, c2-6 to c2-8, d6-8, e7-26) and, for r4-15, in an elapsed-time figure. The random-start minimisations of
`ch07_mass_vacuum/v3_minimisation.py` give the same end points as in version 1.0.0: one run at x = 22 stops at the higher critical value
E = −51.562388, and 23 of 100 runs at x = 25 end below the Standard-Model energy, as stated in Chapter 7.

**New or changed in version 1.0.1** (errata items in brackets): the exact certificates `ch05_clifford_tower/a9_blockgamma_certificates.py`
(a9-1, a9-2; item 37) and `ch06_chiral_matter/a10_chiral_certificates.py` (a10-1 to a10-3; item 52); the numerical
`ch07_mass_vacuum/a11_kernel_structure.py` (a11-1 to a11-3; items 53, 56, 185); the exact
`ch07_mass_vacuum/a13_tensor_certificates.py` (a13-1, a13-2; items 60, 59); a3-9 and a3-10 made exact (item 32); c1-12 added, c1-4 sampled on
S⁺ and c1-5 extended to product states (items 120, 121); c2-8 derives the mass shell from the lapse variation (item 117); `p2_alignment.py`
exits on the aggregate result (item 84); e7-26 compares mean with mean only (item 105); the boundary cases and counterexamples of the review
of the prepared version, all exact: `ch09_b_field/a12_review_boundary_cases.py` (a12-1 to a12-3; items 174, 175 and the current-algebra
remark after Proposition 9.36), a3-11 (Proposition 4.9, item 22), d7-4c (Proposition 9.22(ii)) and d6-8 with M_pair = |sin φ| (item 173).

**Generated from this run**: `verification/MANIFEST.tsv` (1405 checks: 1234 exact, 171 numerical; 662 cited in the manuscript by
identifier), `reproducibility/manifest.json`, `verification/COVERAGE.tsv` (`coverage_summary.txt` here), the theorem-to-code table
`verification/THEOREM_CODE.tsv` (464 labelled statements, 362 with a check) and the script list `verification/SCRIPTS.tsv`, the indexes
`../EXACT_CHECKS.tsv` and `../NUMERICAL_CHECKS.tsv`, and the claim tables of `verification/ch11_curvature_gravity/` and
`verification/ch12_causality/`.
