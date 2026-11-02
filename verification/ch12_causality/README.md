# Chapter 12 checks — Three Temporal Directions and Causal Reduction (Phase 9C, 2 November 2026)

```bash
python3 verification/ch12_causal/c1_selector_types.py
```

```bash
python3 verification/ch12_causal/c2_dynamics_brackets.py
```

| script | checks | runtime | content |
|---|---|---|---|
| `c1_selector_types.py` | c1-1 … c1-11 | ~40 s | thm:selector; the pairing ω_S; thm:cone; thm:selector3 (512×512, exact; two random-spinor samples); prop:types; thm:degree; thm:flat; prop:kineticB; prop:positivealt |
| `c2_dynamics_brackets.py` | c2-1 … c2-9 | ~10 s | prop:toy; prop:orient; ex:tachyons and the mode content of the constraints (CORE-07); prop:constraints and thm:wlgauge (brackets, Jacobi, weak preservation, gauge invariance); prop:wltau; prop:isometryonly |

Uses `verification/ch11_curvature/clifford8.py`. Exact symbolic or integer arithmetic, except the random-spinor samples of c1-4 and c1-5
(numerical, seed printed). Claim-level manifest: `CLAIMS.tsv`. The checks establish algebraic and kinematic statements and the bracket
algebra of a test worldline; they do not establish a causal, stable or consistently interacting theory.
