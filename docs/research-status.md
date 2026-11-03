# Research status (version 1.0.0)

**Status: preprint, not peer reviewed.** Version 1.0.0, prepared 2 November 2026. DOI: [10.5281/zenodo.23900465](https://doi.org/10.5281/zenodo.23900465) (Zenodo record). This is a mathematical research
manuscript; it is **not** a completed unified physical theory, and it claims no experimental confirmation.

## What the manuscript establishes (mathematically)

- At one constituent, the joint commutant of the interior Lorentz algebra and the internal bivectors fixes the imaginary unit up to sign,
  and real Dirac theory on Cl(3,3) is ordinary complex Dirac theory; charge conjugation is the surplus plane (Chapter 2).
- One vector of R^{3,3} has two readings, a displacement and an electromagnetic field strength with the binding root G as Hodge star;
  they share the metric and the Clifford action (Chapter 3).
- At three constituents an axial projector of rank 64 carries one chiral Pati–Salam family with a right-handed neutrino; the electric
  charge is identified among its compact generators (Chapters 6–8).
- Curvature identities of the split frame reproduce the Einstein tensor of a given metric; the selector, light-cone and constraint
  algebras are computed exactly (Chapters 11–12).
- A list of proved obstructions: see [`known-obstructions.md`](known-obstructions.md).

## Evidence classes

Every labelled statement of the manuscript is mapped to its evidence in Appendix D (`verification/THEOREM_CODE.tsv`).

| class | meaning | count (464 labelled statements) |
|---|---|---|
| exact | at least one exact check (integer, rational or symbolic arithmetic, exact modular ranks, blade parity) | 300 |
| exact + numerical | exact and floating-point checks | 31 |
| numerical | floating-point checks only (stated tolerances, spectral gaps, sampled inputs, minimisations) | 29 |
| proof in text | a written proof, no check | 86 |
| argument in text | an argument given where the statement is made, no check | 14 |
| inherited | a result of the literature, used with its hypotheses | 3 |
| conjecture | stated as a conjecture | 1 |

The verification suite has 61 scripts and 1389 checks (1220 exact, 169 numerical); all pass in `reproducibility/selected-results/v1.0.0/`.
A passing check certifies the coded assertion, and a statement only to the extent that the assertion is the statement. Numerical checks
do not prove what they test.

Partial evidence that remains: the electroweak-stage statements of Chapter 7 (Propositions 7.22, 7.35, 7.38, 7.43 and the structural part
of Theorem 7.44) are recomputed numerically, not certified exactly; in Proposition 7.22 the identification of the spectrum of ΦᵀΦ with that
of S†S rests on an explicit isomorphism that no check reproduces; the one-loop formula of Proposition 7.21 is inherited; Example 7.15 has
no check.

## Conjectures

- Global minimality of the selected vacuum below the crossing value x_* (Proposition 7.33 records the numerical searches that support it).
- A transmuted binding length (Conjecture 9.17 [bf:transmute]).

## Unresolved physical questions

- No gravitational dynamics, matter coupling or low-energy limit is derived; the gravitational action is adopted (Chapter 11).
- The photon is identified only under an adopted potential; the W/Z/photon mass splitting is not obtained (Chapter 8).
- The chiral projector P_* and the clock reduction are chosen, not derived.
- Causality is algebraic; constraint closure is shown for a test worldline only; no interacting multi-time theory and no local clock
  symmetry are constructed (Chapter 12).
- Quantisation of the field-strength reading is not addressed.
- The open problems of Section 16.3 and the hypotheses listed in [`mathematical-assumptions.md`](mathematical-assumptions.md).

## AI assistance

The manuscript and its verification package were prepared with substantial assistance from AI systems, under the author's direction.
Chapter 15 describes how: which steps were assisted, how claims were classified, how the executable checks were written and run, and
what a passing check does and does not show. All results are the author's responsibility.
