# Regression of version 1.0.0 (2 November 2026)

**Command** (from the repository root; the scripts run are byte-identical to the recorded revision, as `environment.txt` records no
modification of `verification/`):

```bash
python3 verification/run_all.py --out reproducibility/selected-results/v1.0.0 -j 3
```

**Environment** (`environment.txt`): Python 3.13.18, numpy 2.5.3, scipy 1.18.1, sympy 1.14.0, mpmath 1.3.0; Linux x86_64 (glibc 2.39);
3 workers; one BLAS thread per script. Wall time 17 min 57 s.

**Result** (`summary.tsv`): **61 scripts, 1389 PASS, 0 FAIL**; every exit status is zero and consistent with the PASS/FAIL lines. One script
(`ch08_gauge_photon/p3_identification.py`) writes its timing lines to standard error. An earlier run of the same version gave the same
check identifiers, outcomes and claim texts (one claim differs only in an elapsed-time figure). The random-start minimisations of
`ch07_mass_vacuum/v3_minimisation.py` include one run at x = 22 that stops at the higher critical value E = −51.562388 and 23 of 100 runs
below the Standard-Model energy at x = 25, as stated in Chapter 7.

**Generated from this run**: `verification/MANIFEST.tsv` (1389 checks: 1220 exact, 169 numerical; 643 cited in the manuscript by
identifier), `reproducibility/manifest.json`, `verification/COVERAGE.tsv` (`coverage_summary.txt` here), `verification/THEOREM_CODE.tsv`
and Appendix D tables D.1–D.2 (464 labelled statements, 360 with a check), the indexes `../EXACT_CHECKS.tsv` and `../NUMERICAL_CHECKS.tsv`,
and the claim tables of `verification/ch11_curvature_gravity/` and `verification/ch12_causality/`. `audit_initial_status.tsv` gives the
current evidence class of 124 statements that had neither a check nor a written proof at an earlier stage of the work.
