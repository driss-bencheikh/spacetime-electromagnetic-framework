# Mathematical assumptions

Everything the manuscript derives is derived inside the real Clifford algebra of the split space R^{3,3}. Everything it does not derive is
stated as an assumption, an adopted input, an inherited result or a hypothesis. This page lists them. Numbers refer to
`release/manuscript-v1.0.1.pdf` (the numbering of version 1.0.0 is unchanged); labels in brackets are the LaTeX labels in `manuscript/`.

## The arena (definitions, not assumptions)

- V = R^{3,3} with an orthonormal basis e₁, …, e₆, e₁–e₃ of square +1 and e₄–e₆ of square −1, and the pairing eᵢ ↦ eᵢ₊₃.
- Cl(3,3) ≅ M₈(R) and the graded tower Cl(3n,3n) ≅ M_{2^{3n}}(R), with the explicit integer realisation of Appendix A.
- Blade commutation by index parity (eq. bladecomm); the conventions of the Introduction ("Conventions", "Notation").

## Explicit assumptions of the text

| assumption | where | what it supplies |
|---|---|---|
| Clock reduction | Assumption 2.16 [assum:clockred] | the choice of one clock (a unit temporal direction) through which the surplus times are reduced |
| A quantised leaf | Assumption 4.51 [cl:quant] | the quantisation used for the closure of Chapter 4 |
| Newtonian benchmark | Assumption 4.55 [assum:newton] | the classical benchmark for the orbital extension of Chapter 4 |
| One field, one Yukawa coupling | Assumption 7.19 [ass:field] | a single scalar field on the kernel space with one Yukawa coupling and one kinetic normalisation, although symmetry allows one per irreducible component |

## Selections adopted rather than derived

- **The chiral projector P_*** (Section 6.2 [sec:pstar]): chosen to cut one chiral Pati–Salam family out of the third rung; no dynamical
  selection of P_* is derived.
- **The entropy selection of the symmetric clock Ŵ** (Chapter 2): rests on Jaynes' maximum-entropy principle, used as a hypothesis.
- **The vacuum potential and the Maxwell term** (Chapters 7–8): the selection of the Standard-Model orbit (Section 7.5 [sec:selection]) and
  the photon identification hold for an adopted potential and an adopted Maxwell term. The single gauge coupling of the adopted action
  (equation 7.4 [eq:Saction]) is a choice: it fixes equal normalised couplings at an unspecified scale, the relation that the spectral action
  of Chamseddine, Connes and Marcolli derives as a boundary condition at its cut-off.
- **The electroweak stage** (Proposition 7.35 [prop:EW]): tree level, for a given bidoublet value, with adopted couplings.
- **The gauge principle that removes the surplus times** (Chapter 12): weak emergence for scalar matter holds under the chosen gauge
  principle (Corollary 12.14 [cor:weak]), whose constraint preservation is shown for a test worldline only.
- **The gravitational action** (Proposition 11.14 [prop:holst]): the Palatini–Holst action is adopted, not derived.

## Inherited results (used with their hypotheses)

Results of the literature are used as such and marked *inherited* in the text, for example the Helmholtz conditions behind the
no-Lagrangian statement (Theorem 11.13 [thm:helmholtz]), the Coleman–Weinberg one-loop potential in Proposition 7.21
[prop:coupled3], the Klein–Gordon theory behind Corollary 12.14, Wigner's classification (Chapter 3), the reductive dual pair and its Casimir
relation (Howe 1989; Neergård 2019), and standard representation theory; Section 14.5 [sec:rel:pointers] lists the comparisons with the
literature chapter by chapter. Appendix D lists, for every labelled statement, whether its evidence is a check, a written proof, an inherited
result or a conjecture.

## Hypotheses of the author (not results)

- **The reading of E = mc²** (Introduction, Section I.2.1): c² as the symmetric component and m as the symmetry-breaking component,
  hypothesised to arise through the b-particle associated with the Yang–Mills b-field.
- **The energy cascade** (Section I.2.2): the b-particle, induced by the b-field through the interaction of two or more manifolds, proposed
  as the mechanism enabling an energy cascade from Cl(3,3) to Cl(6,6) and beyond.
- **A universal evolution parameter** (Section 12.6 [sec:universaltau]): prior art in the parametrised relativistic dynamics of Horwitz and
  Piron (1973) and Horwitz (2020); what may be specific here is its connection with the clock reduction. It is stated with its requirements
  and is not a result.

The manuscript tests what the algebra fixes, admits or obstructs under these hypotheses; it does not establish them.
