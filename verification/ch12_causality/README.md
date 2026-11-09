# Chapter 12 checks — Three Temporal Directions and Causal Reduction (Phase 9C, 2 November 2026)

```bash
python3 verification/ch12_causality/c1_selector_types.py
```

```bash
python3 verification/ch12_causality/c2_dynamics_brackets.py
```

| script | checks | runtime | content |
|---|---|---|---|
| `c1_selector_types.py` | c1-1 … c1-12 | ~40 s | thm:selector; the pairing ω_S; thm:cone; thm:selector3 (512×512, exact; random spinors of S⁺, product-state currents); prop:types; thm:degree; thm:flat; prop:kineticB; prop:positivealt; the null centraliser of thm:selector (iii) (c1-12, v1.0.1) |
| `c2_dynamics_brackets.py` | c2-1 … c2-9 | ~10 s | prop:toy; prop:orient; ex:tachyons and the mode content of the constraints; prop:constraints and thm:wlgauge (brackets, Jacobi, weak preservation, gauge invariance); prop:wltau, with the mass shell derived by varying the lapse (v1.0.1); prop:isometryonly |

Uses `verification/common/clifford8.py`. Exact symbolic or integer arithmetic, except the random-spinor samples of c1-4 and c1-5
(numerical, seed printed). Claim-level manifest: `CLAIMS.tsv`. The checks establish algebraic and kinematic statements and the bracket
algebra of a test worldline; they do not establish a causal, stable or consistently interacting theory.
