# Errata for version 1.0.0

Known issues in `release/manuscript-v1.0.0.pdf` (DOI [10.5281/zenodo.23900465](https://doi.org/10.5281/zenodo.23900465)), found by
independent AI-assisted reviews of the abstract, the Introduction and Chapters 1 to 14 and 16, and by an AI-assisted chapter-by-chapter
derivation check, the author's review of the comparison with the literature and reviews of the prepared version 1.0.1 (Chapters 4
and 9, items 168–179; Chapters 1–3, 5–16 and the Introduction, items 180–191), and checked by the author. Page numbers are the printed
page numbers of version 1.0.0.

**Status (9 November 2026):** every item is corrected in version 1.0.1 (`release/manuscript-v1.0.1.pdf`); items 59 and 60 by the exact
checks a13-2 and a13-1 of `verification/ch07_mass_vacuum/a13_tensor_certificates.py`. Theorem, proposition, table and
figure numbers are those of version 1.0.0 (Figure 13.2 is new); page numbers differ.

**None of these issues invalidates a main theorem.** Each correction is local. Where a later chapter repeats an affected statement or cites
it as justification, the item says so (items 53, 65, 68 and 77); in every such case the later results hold once the correction is made.
Otherwise a later chapter that uses an affected result uses it in a case where the result holds, or already states the missing hypothesis.

Items are numbered in the order in which they were recorded and listed by chapter. Items are marked **Error** (a statement false as written), **Hypothesis** (a missing condition or domain restriction), **Wording**
(imprecise or misleading phrasing), **Attribution** (a missing or incomplete reference) or **Evidence** (a proof or check weaker
than the statement it supports).

## Abstract and Introduction

153. **Wording. Abstract.** Several statements need their scope: the imaginary unit is fixed up to sign under stated commutation requirements
     with the chosen gamma matrices; "no equivariant map joins" the two readings means no nonzero Lorentz-equivariant linear intertwiner between
     the specified modules; binding "only in neutral signature" and the unique rotor closure hold for the factorised shift model and within the
     closure ansatz; the electric charge is identified at selected vacua and the anomaly is computed on the chosen matter block; selection results
     hold for specified potentials and sectors while the global radiative claims are conjectural; "without a gravitational action" means that no
     gravitational dynamics is derived (a Palatini–Holst action is adopted in Chapter 11); and the gauge principle removing the surplus times is
     formulated for a test worldline. The keyword "computer-verified mathematics" overstates the evidence (Chapter 15). The abstract will be
     rewritten in version 1.0.1.
154. **Error. §I.1.2 (p. 11).** "The finite-dimensional statements are exact: every blade … is a signed permutation matrix, so products,
     commutators, traces and ranks are integer computations" is false as a universal statement: linear combinations such as
     (e₁+e₂+e₃)/√3, parameter-dependent spectra, minimisations and numerical ranks are not integer computations, and 169 of the 1389 checks are
     numerical (Chapter 15).
155. **Error. §I.2.4 (p. 16).** "Each section has a verification script … and a second, adversarial check script … exit with status zero only if
     all pass" contradicts Chapter 15: two sections were refereed by reading only, no referee script exists for p1 and p2, and `p2_alignment.py`
     has no exit statement (items 84, 144).
156. **Wording. §I.1.3 (p. 11).** The gauging/grading dichotomy needs the restrictions of items 65 and 66. "At one constituent this centraliser is
     the single bivector J, so the photon is not one symmetry of the interior structure but the only one", and "the electromagnetic u(1)", identify
     the interior U(1) with the photon, which is not proved (item 67): a compact U(1) candidate still needs a connection, a kinetic term, a matter
     coupling and physical degrees of freedom. "The kinetic forms alone fix the massless gauge algebra" is subject to item 8.
157. **Wording. §I.1.3 (p. 12).** E₀ is the infimum of E₆ over boosts, attained for non-null fields only (item 130), and a field-point identity is not
     a particle mass; the extension of the force law is unique within its ansatz; and "the quantum state space must be formed by complexifying
     from outside" concerns the real CAR representation and observables constructed in Chapter 2 (item 9).
158. **Wording. §I.1.3 (pp. 12–13).** The kernel decomposition is subject to item 53; the spectrum of the 271 angular modes holds within the adopted
     leading-logarithm model (item 62); "an electroweak stage" (a given bidoublet value in the kernel model) and the electroweak-symmetric vacuum of
     the projector model (Chapter 8) belong to different models, whose common realisation is not established; and numerics support, not locate,
     the transition between the sufficient and necessary conditions.
159. **Wording. §I.1.3 (p. 13).** Ill-posedness concerns unrestricted six-dimensional evolution (constrained ultrahyperbolic problems can be well
     posed, Craig and Weinstein 2009); "no bounded-below invariant potential" means invariant polynomial, and "every invariant kinetic term"
     the unconstrained adjoint-scalar ansatz (item 125); and the surplus-time constraints follow from a worldline gauge symmetry that is chosen,
     not derived, and is shown for a test worldline only (Chapter 12).
160. **Wording. §I.1.2 and §I.2.4.** The manuscript uses three classification schemes: four classes with three provenance labels (§I.1.2), six
     classes (§I.2.4) and five link statuses (Chapter 13); they will be reconciled. The principle that a statement "carries information only if it
     could have come out otherwise" is a requirement on physical claims (that they add a constraint beyond a change of notation or a freely chosen
     embedding), not a definition of mathematical content.
161. **Wording. §I.2 (pp. 14–16).** The second section restarts under the heading "The question" and gives the historical rest-energy hypothesis
     prominence without saying at once what survived the tests (Chapter 13, §13.5); its dictionary entries (c² as a symmetric component, b as the
     Yang–Mills field) are conjectural readings. It will be retitled as historical motivation and shortened.

## Chapter 1 — Split-Signature Clifford Geometry

65. **Error. Corollary 1.4 and the paragraph after it (p. 18).** That exp(θB) = cosh θ·1 + sinh θ·B generates ℝ and not U(1) is correct. The
    paragraph's inference, that a positive-definite kinetic metric on the gauge algebra excludes a gauge field along a boost, fails for an
    abelian algebra: on ℝB every scalar product is invariant, and A_μ = a_μB has the ordinary healthy abelian kinetic term. The valid
    obstruction concerns the action on S: no positive-definite H satisfies BᵀH + HB = 0 (if Bv = ±v then v·Hv = 0), so a boost is not
    unitary for a fixed positive form on the module, and the gauge form −Tr has the wrong sign on B. "Gauge-admissible" should be defined as
    compact, as a condition of the model or through unitarity on the matter module, and the corollary retitled "A mixed bivector generates
    no compact U(1)". Theorem 5.39(iii) (p. 154), the gauge-admissibility paragraph of Chapter 9 (p. 267) and its later citations, and
    Chapter 13 use the corollary in exactly this sense; their conclusions stand under that definition.
66. **Wording. p. 18.** "Square −1 gauges and square +1 grades, and no bivector does both" is too strong: conjugation by J = e₂e₃ (J² = −1)
    is an involutive automorphism of Cl(3,3) with ±1 eigenspaces of dimension 32 each. The distinction is between operators on S (J is a
    complex structure, each ωᵢ has real ±1 eigenspaces) and conjugation on the algebra (either kind induces an involution). Proposition 1.5
    is unaffected.
67. **Wording. Figure 1.1, caption.** "J … will become the electromagnetic charge" anticipates more than is proved. Proposition 1.7 shows that
    J spans the centraliser of the interior so(1,3) in so(3,3) (in End(S) the centraliser is eight-dimensional, as the text notes). J is the
    compact internal rotation used as complex structure in Chapter 2; the electromagnetic generator is constructed in Chapter 8. Later chapters
    use "the charge J" for the U(1) of the Dirac reading and do not identify it with electromagnetism.
68. **Wording. Chapter opening (p. 17) and Theorem 1.13 (p. 23).** "This decides the sign of every kinetic term later on" holds only for
    actions built on this trace metric; the theorem does not fix an arbitrary action, spacetime contraction, constraint or propagator residue.
    Theorems 1.3 and 1.13 classify bivectors relative to the fixed Cartan decomposition (compact = antisymmetric for the chosen positive form);
    a mixed generator gives a tangent to the orthogonal projectors exactly when its symmetric part commutes with P, so generators should be
    counted modulo the stabiliser. The same qualification applies to the sentence "exactly when B is compact" in Theorem 8.54 (p. 259), whose
    dimension counts are unaffected. The opening's "makes no physical reading" should also note that Proposition 1.10 uses the Fisher–Rao clock
    of Chapter 2 and equation (1.6) uses the E/B reading of Chapter 3.
69. **Wording. Example 1.6 (p. 19).** The paired-mode projector should carry its limitation where it is introduced: among the interior Lorentz
    generators (e₁e₄, e₁e₅, e₁e₆, e₄e₅, e₄e₆, e₅e₆) only e₁e₄ commutes with Π (commutator ranks 0, 4, 4, 4, 4, 4), and [Π, J] has rank four. Its
    stabiliser in so(3,3) is nine-dimensional, which counts linear combinations, not only commuting basis bivectors.
70. **Wording. p. 19.** The four-sector grading of Proposition 1.5(iii) is the coordinate parity Q(D₃)/2ℤ³, of order four, not the quotient of
    the root lattice by twice itself, Q(D₃)/2Q(D₃), of order eight.
71. **Wording. Proposition 1.10 (p. 22).** If k coordinates of x vanish, the top eigenspace of B_x has dimension 2^k (8 at x = 0), and the map
    from the clock sphere to the simplex has 2^{3−k} preimages over a face: a quotient with boundary stabilisers, not an eight-sheeted covering.
    The radius-two coordinates xᵢ = 2√wᵢ give B_x = 2ϖ/√3 at the symmetric point, the unit clock ϖ/√3; the top eigenspaces agree.
72. **Wording. Proposition 1.12 (p. 22).** The spin lift acts as the derived action dΛ^w(A) − ½ tr A, not the exterior power of the matrix A;
    for A = 1 the eigenvalue on degree r is r − 3/2. On GL⁺(L) the integrated action is ρ(g) = (det g)^{−1/2}Λ*(g).
73. **Wording. Notation in Chapter 1.** The identity ω = 𝒱 = Γ should be stated once, with a different symbol for the vector space V. The
    paragraph after Proposition 1.8 writes I = −ω₁ω₂ω₃ for what the proposition calls Γ. The letter w is 3 − N in Proposition 1.8(v) and N in
    Proposition 1.11; T is e₁e₂e₃ in the proof of Proposition 1.9 (defined only in Chapter 2), the commutator [B, P] in Theorem 1.13 and the
    boosts T_a. In Proposition 1.8(v), the antecedent of "which is the chirality γ₅" should be made explicit.

## Chapter 2 — Lorentz Structures, Complex Units and Observational Spacetime

1. **Error. Theorem 2.7 (p. 26).** "Four values, each of multiplicity two" holds only when the four values are distinct. When values
   coincide their multiplicities add: for C = I₃ the spectrum of B_C² is 9 (twice) and 1 (six times). Correct statement: the characteristic
   polynomial is ∏_{ε₂,ε₃=±1} [t − (σ₁ + ε₂σ₂ + ε₃σ₃)²]². Corollary 2.8 is unaffected.
2. **Wording. Example 2.9 (p. 27).** The squared eigenvalues are non-negative, so the factorised equation has growing or affine modes,
   not bounded oscillations; read "up to four rates", with X ∈ S ≅ R⁸.
3. **Hypothesis. §2.1.3 (p. 25).** The Fisher metric is singled out by invariance under *congruent* Markov morphisms (Čencov), not under
   arbitrary stochastic maps; it is a smooth metric on the open simplex, and the closed octant is its metric completion.
4. **Wording. Proposition 2.13(iv) (p. 29).** Reversing the clock, n → −n, at fixed six-momentum is a clock-orientation reversal, distinct
   from the physical time reversal of Theorem 2.27.
5. **Hypothesis. Theorem 2.31 (p. 36) and the projectors Λ± = (I ± H/E)/2.** Require E > 0; at E = 0 the projectors are undefined and all
   sixteen energies coincide.
6. **Wording. Theorems 2.25 and 2.27 (p. 34).** The normalisation of C, P and T should be part of the statements.
7. **Wording. Theorem 2.28 and Figure 2.2 (p. 35).** The circle e^{θJ}Ω is an orbit under left multiplication (a choice of phase), not under
   conjugation, since [J, Ω] = 0; "four real points" should read "four distinguished points".
8. **Wording. Proposition 2.29 (p. 35).** The algebra computed is the internal symmetry algebra of the free kinetic action; it becomes a gauge
   algebra only once a connection is introduced and quantum consistency is checked.
9. **Wording. Theorem 2.30 (p. 36).** The absence of a commuting complex structure is a statement about the chosen real representation;
   the CAR convention is clarified by χᵢ = iθᵢ, which satisfy {χᵢ, χⱼ} = 2δᵢⱼ.
10. **Error. Proposition 2.32 (p. 37).** "A product of four integral reflections" holds for the representative C = e₂, not for every
    C = cos θ e₂ + sin θ e₃; in general four orthogonal reflections suffice. The stated dimensions 8 and 3 are real dimensions on R¹⁶; on the
    complexified space F they are complex dimensions.

Added from the derivation check (5 November 2026):

74. **Wording. Theorem 2.27 (p. 35).** "Anticommutes with every vector C of the clock plane" should read "of the surplus plane U".
75. **Wording. Theorems 2.27 and 2.28.** T denotes both the temporal trivector e₁e₂e₃ and the time-reversal operator; in P = γ⁰J the two
    meanings meet in one formula.

## Chapter 3 — Classical and Electromagnetic Readings of X

11. **Error. Lemma 3.43 (p. 61).** A single Lᵢ does not annihilate exactly the two pole states: its kernel is four-dimensional,
    span{vac, full, fᵢ, gᵢ}. The joint kernel of L₁, L₂, L₃ is span{vac, full}. The representation decomposition that follows is unaffected.
12. **Wording. Corollary 3.9 (p. 42).** The source-free equation δF = 0 alone also needs no external data. The correct selection statement:
    among nonzero real subspaces of the two operator equations, only the whole space is invariant under the continuous duality rotation.
13. **Wording. Corollary 3.50 (p. 64).** The generators Qᵢ and −T are odd; finite transformations generated by them are in general of mixed
    parity (for example e^{−θT} = cos θ I − sin θ T).
14. **Hypothesis. Theorem 3.95(ii) (p. 94).** For signed rapidity read |S| = E₀ |sinh 2χ| and v_E = |tanh 2χ|, or restrict to χ ≥ 0.
15. **Hypothesis. Theorem 3.97(i) (p. 95).** "Every stabiliser element satisfies N³ = 0 ≠ N²" applies to the nonzero elements.
16. **Hypothesis. §3.8.4 (p. 89).** E·B ≠ 0 implies an electric-type invariant, but vacuum decay into pairs is established for constant
    (or controlled locally constant) backgrounds; for general backgrounds, spatial extent, time dependence and available energy matter.
17. **Error. §3.8.4 (p. 89).** The Lagrangian L = S + c₁S² + c₂P² is not duality invariant in general: its constitutive duality defect is
    4(c₁ − c₂)SP + 4c₁(c₁ − 2c₂)S²P − 4c₂²P³, which vanishes to first order only if c₁ = c₂. The phase-velocity formula
    1 − cᵢΘ^{μν}k_μk_ν needs a normalised propagation direction.
18. **Hypothesis. Theorem 3.21 (p. 48).** L₀ is non-differentiable at *nonzero* null fields (at F = 0 its derivative vanishes).
19. **Wording. Proposition 3.23 (p. 49).** The ellipticity statement applies to a gauge-fixed (or transverse) system.
20. **Wording. Proposition 3.91 (p. 91).** Distinguish the group U(3) from its Lie algebra u(3), and write the effective action as
    O(3,C)·U(1).
21. **Attribution. §3.9.** The local invariant u₀ = ½|z·z| with an inertia and rest-energy interpretation, energy-flow velocity and
    standing-wave examples appears in G. Kaiser, *J. Phys. A* 44 (2011) 345206, which should be cited and compared directly.

## Chapter 4 — Classical Rigid-Body Realization

22. **Error. Proposition 4.9 (p. 104).** False as stated: an eigenvector pair of an η-self-adjoint G need not span a plane of signature
    (1,1) (it can be totally null), and for non-diagonalisable G the bounded oscillating subspace need not have signature (k,k). Counterexample:
    η = [[0,K],[K,0]] with K = diag(1,−1) and G = [[J,I],[0,J]], whose bounded motions lie in a plane on which η vanishes. Correct statement:
    for diagonalisable G the oscillating subspace is nondegenerate of signature (k,k) with k ≤ min(P,Q), so the conclusions for G² = −I,
    which are the only ones used later, stand. *Amended after a review of the prepared version 1.0.1:* the oscillating subspace is the
    real span of the eigenvectors with nonzero purely imaginary eigenvalues; a non-real spectrum is not enough (η = diag(1,−1),
    G = [[1,−1],[1,1]], eigenvalues 1 ± i, gives X = e^τ(cos τ, sin τ) for κ = 1).
23. **Error. Canonical transformation after Theorem 4.12 (p. 106).** The shift that removes the boundary term is p ↦ p − mκηGX, not
    p ↦ p + 2mκηGX (the latter relates the +κ and −κ gauges). Theorem 4.12 and the absence of a ground state are unaffected.
24. **Hypothesis. Proposition 4.33 (p. 113).** "Every frame moves on a bounded orbit" on the balance surface requires h(0) ≠ 0; with
    h(0) = 0 the frames drift linearly. (Theorem 5.52 and the corresponding statements of Chapter 9 already include h(0) ≠ 0.)
25. **Error. §4.4.7 (p. 115).** "A quadratic Lie algebra whose invariant form is η must have a maximal compact subalgebra of half its
    dimension" holds for *simple* algebras only (so(4) with opposite-sign forms is a counterexample). The classification that follows is
    correct: on R^{3,3} exactly so(1,3), so(4) and so(2,2) carry an invariant form of signature (3,3), only so(1,3) with η proportional to the
    Killing form, and no exceptional real form is balanced.
26. **Hypothesis. Theorem 4.21 (p. 109).** Exclude the identically zero rate; the statement concerns a nontrivial force law on an open set.
27. **Hypothesis. Theorem 4.40 (p. 116).** E₀ = |w|/2 is not smooth at w = 0; E₀² is the smooth Casimir, and the independence of the four
    functions holds generically.
28. **Error. Definition 4.54 and Proposition 4.57 (pp. 125–127).** The spin angle is measured from the line to the primary, so the spin
    angular momentum is h + Cn on a circular orbit (h + C ḟ on an eccentric one), not h. The pendulum equations still agree.
29. **Hypothesis. Proposition 4.56(iii) (p. 126).** The oscillator–Kepler map holds away from collisions (z ≠ 0); the equal-amplitude orbit is
    a regularised collision.
30. **Wording. §4.7.4, Table 4.3, Figure 4.3.** "Captured" overstates the result: the model has no tidal dissipation, and the run shows
    libration in the 3:2 resonance for the chosen initial data.
31. **Attribution. Theorems 4.8 and 4.34.** The linear algebra (a complex structure anti-isometric for η forces neutral signature; the
    complex-orthogonal centraliser) is that of Norden (anti-Hermitian) geometry, e.g. Borowiec, Francaviglia and Volovich,
    *Differential Geom. Appl.* 12 (2000) 281; it should be presented as inherited.
32. **Evidence. Theorems 4.38 and 4.39.** The check of Theorem 4.39 uses ranks modulo two primes without an exact kernel certificate, and
    the representation content in Theorem 4.38 is read by floating-point diagonalisation. The values of Theorem 4.39 have since been confirmed
    by independent exact rational elimination; exact certificates will be added.
168. **Hypothesis. Theorem 4.21 (p. 109).** The state-dependent rate must not vanish identically on the open set being compared; a rate
     that vanishes there and not elsewhere would let the zero Hamiltonian coincide with the force law locally. The Casimir argument proves the
     local statement (continuity gives an open subset on which the rate is nonzero); the computation at twelve null vectors is a check, not the
     proof.
169. **Hypothesis. Theorem 4.53(vi) (p. 125).** "H is unbounded below" needs ν ≠ 0: at ν = 0, N = 0 and H = h²/(2I₀) is bounded below.
170. **Wording. §4.2.3 (p. 103).** The stress tensor is the usual one when the fibre metric and the bundle connection are held fixed in the
     metric variation (a metric-dependent connection needs further terms), and its trace "vanishes only in D = 2" should read "vanishes
     identically for arbitrary fields only in D = 2".
171. **Wording. Opening of §4.2 (p. 101).** "One piece of interior dynamics without further input" and "a single choice" overstate what the
     algebra gives: for the adopted shifted quadratic action, the comparison concerns the choice of a positive or a split quadratic form.
172. **Attribution. After Theorem 4.53 (p. 124).** The conserved total momentum J_s is the momentum map of semidirect-product mechanics (Holm,
     Marsden and Ratiu 1998); "new" should read "computed here for this closure".

## Chapter 5 — The Clifford Tower and Constituent Duplication

33. **Attribution. Theorem 5.11 (p. 134) and Theorem 5.68(iii) (p. 180).** The constants n(n+2)/2 and 3n(n+4)/4, and the general form
    nd(d+n−2)/8, are special cases of the relation C_o(d) + C_o(2k) = ¼kd(d+2k−2) between the Casimirs of the o(d)–o(2k) dual pair on a
    fermionic Fock space, K. Neergård, *J. Math. Phys.* 60 (2019) 081705, eq. (73). The statement that 3n(n+4)/4 is new is withdrawn, here
    and in the claims table of Chapter 9 (p. 319). The Clifford-algebra proof stands as an independent derivation in the real realisation.
34. **Error. Theorem 5.77 (p. 186).** The identity ∂_μj^μ + ψᵀJ₁𝓔_ψ = 0 holds with 𝓔_ψ the Euler–Lagrange derivative of the action,
    𝓔_ψ = 2γ̂⁰(J₁γ̂^μD_μψ − mψ), not with the Dirac expression J₁γ̂^μD_μψ − mψ itself. The field equations are unaffected.
35. **Error. Proposition 5.44(ii) (p. 160).** "Each on 16 dimensions when b > 0" holds when the four values ±E₊, ±E₋ are distinct; when
    they coincide the multiplicities add. E₋ = 0 exactly when m₁m₂ = 0 and E₁ = E₂: for m₁ = m₂ = 0 and p = q = (1,0,0) the spectrum is
    −2 (16), 0 (32), 2 (16). The eigenvalue formulas are correct.
36. **Error. Theorem 5.57 and Example 5.58 (pp. 171–172).** M_pair, the infimum over observers of the momentum density of the summed
    stress tensor at a point, is a local rest-energy density, not the invariant mass of two photons. For two null fields it equals
    √(u₁u₂)(1 − cos ϑ), whereas two photons have invariant mass √(2E₁E₂(1 − cos ϑ)); for z₁ = (1,i,0) and z₂ = (0,1,i) the values are 1
    and √2, and the counter-propagating example agrees only by coincidence. The sentence "the invariant mass of the pair is M_pair" is
    withdrawn, the Breit–Wheeler passages are kinematic analogies, not a production process, and Chapter 9 uses the same name in the same
    sense. Theorem 5.57 itself is unaffected.
37. **Evidence. Theorem 5.13 (p. 135).** The proof asserts {Q_μ,Q_ν}|_D = −72η_μν and ¼[G_a,G_b] = Δ_ab|_D without deriving them. Both are
    correct (confirmed by an independent construction); a derivation or an exact certificate will be added.
38. **Wording. Theorem 5.39(ii) (p. 154).** "Conserves no part of it" means that a nonzero paired boost X has no component commuting with N,
    so [N,X] ≠ 0. It may still conserve partial occupation numbers: (e₃e₅ − e₂e₆)/2 conserves N₁ and N₂ − N₃.
39. **Wording. p. 148.** "Copying it is impossible" refers to the cloning of an unknown quantum state (Wootters–Zurek); the linear map
    z ↦ (z,z) of classical field data is well defined.
40. **Hypothesis. Proposition 5.63 (p. 175).** The cancellation of zero-point contributions concerns the nonzero-frequency lines of a
    Hamiltonian that is unbounded below; at the zero-frequency resonance there is no oscillator vacuum. No vacuum-energy conclusion follows
    (compare Theorem 4.12).
41. **Wording. Theorem 5.15 (p. 136).** "Compact" means contained in the chosen compact orthogonal subgroup, the one compatible with a
    positive kinetic form, in the stated embedding; it is not a statement about abstract Lie algebras. A compact candidate centraliser is
    not by itself a gauge theory.
42. **Wording. Particle language (e.g. p. 156, p. 176).** "The rung carries one particle" and "the composite is a fermion" are
    interpretations or spin-parity statements under the stated assumptions, as the chapter's closing synthesis says. "No interacting
    theory is constructed" means no derived quantum theory of tower-changing processes; Theorem 5.77 constructs a classical interacting model.

## Chapter 6 — Chiral Matter and Internal Symmetries

43. **Error. Figure 6.1, caption (p. 194).** "Every rotation inside each group survives the projection" is false for the auxiliary group:
    it enters P_* through the two-plane volume q_A = e₂⁽¹⁾e₃⁽¹⁾, not through the volume of its four directions, and only the plane rotations
    E₂E₃ and E₈E₉ survive (E₂E₈ anticommutes with −Jq_A and P_*E₂E₈P_* = 0). The surviving bivectors are the 15 colour rotations, the 6
    weak rotations, E₂E₃ and E₈E₉, which is the 23-dimensional algebra of Theorem 6.11; the theorem is unaffected.
44. **Error. Proposition 6.22 (p. 198).** The values λ = 0 and λ = −1 give (Q_u, Q_d, Q_ν, Q_e) = (2/3, −1/3, 0, −1) and
    (1/3, −2/3, 1, 0). Exchanging ν and e alone does not map one to the other; charge conjugation together with u ↔ d and ν ↔ e does. The
    integrality and neutrality statements are correct.
45. **Hypothesis. Several generations (p. 198).** "With N_g generations the standard branch has 2N_g parameters" needs the family to be
    specified. If anomalies cancel only after summing over generations, the anomaly variety has local dimension 6N_g − 4 at the
    Standard-Model point (14 for three generations). A 2N_g-dimensional family is obtained with one common Yukawa difference h
    (y_u = q_i + h, y_d = q_i − h, y_ν = l_i + h, y_e = l_i − h, and Σ_i(N_c q_i + l_i) = 0), which contains L_e − L_μ; imposing the one-family
    branch generation by generation also gives 2N_g but excludes L_e − L_μ.
46. **Wording. Electric charge of the dressed benchmark (p. 198).** From Q_em = t³ + B₀ only the projectors P_L, P_R on the weak axis
    cancel; B₀ = (ε/6)Γ(q₁ + q₂ + q₃) still contains Γ. "The chirality dressing cancels identically" should be narrowed accordingly.
47. **Wording. Theorem 6.23 (p. 199).** The twelve solutions are the six differences N_i − N_j (c = 0) and the six signed sums
    ±(N_i + N_j − 1) (c = ∓1); the printed "Q = N_i + N_j − 1 with c = −1" names only one sign of the second orbit. Hodge conjugation maps
    exterior degree k to 3 − k and so exchanges the even and odd halves: the relation CQC⁻¹ = −Q, equivalent to Σq_i + 2c = 0, holds on the
    whole exterior algebra, not as an operation inside the even four-state sector, and the charge skeleton alone does not identify neutrinos
    and electrons.
48. **Hypothesis. Proposition 6.14 (p. 196).** "The largest anomaly-free gauge algebra containing Pati–Salam" is proved inside the
    23-dimensional effective bivector algebra of Theorem 6.11 and for this matter content, as Chapter 7 states it; it is not a classification
    of all algebras acting on the block. The restriction belongs in the statement.
49. **Wording. Proposition 6.14 and the following paragraph (p. 196).** The fermion-number current has a triangle anomaly with SU(2)_L and
    SU(2)_R already in perturbation theory (Adler 1969); what holds to all orders of perturbation theory is the corresponding charge selection
    rule about the trivial topological sector, violated nonperturbatively by instantons ('t Hooft 1976). The two central directions are
    classical global symmetries with mixed anomalies, not exact quantum symmetries. The coefficients are correct.
50. **Hypothesis. Theorems 6.20 and 6.21 (p. 198).** Theorem 6.20 assumes N_c ≥ 2 (for N_c = 1 there is no colour mixed condition), and the
    title of Theorem 6.21 should carry the qualification "on the generic branch" already present in its statement.
51. **Wording. §6.1 (p. 191).** "To obtain chiral matter one must project" holds for the undressed internal bivector action of
    Propositions 6.1 and 6.2; Theorem 6.9 already obtains a chiral interaction with dressing and spectator states.
52. **Evidence. Theorems 6.5 and 6.11 and Corollary 6.12.** The decisive step is the statement that a nullspace has the stated dimension.
    Exact certificates will be added: the 23-generator basis of item 43 and the completeness argument (every blade commutes or anticommutes
    with each of the three involutions), the centre and derived-algebra dimensions for each row of Theorem 6.5, and the colour-volume sign
    rather than the quadratic Casimir to tell 4 from 4̄. The statements themselves have been confirmed by an independent exact computation.

## Chapter 7 — Mass Kernels and Conditional Vacuum Selection

53. **Error. Theorem 7.17(ii) and the paragraph after it (pp. 212–213).** With the complex structure 𝒞(M) = JM on Majorana-type and JMΓ on
    Dirac-type kernels, the two colour tens carry the same conjugation, so the complex representation is not the Pati–Salam branching of
    10 ⊕ 126, nor of its conjugate. A label-independent test: for the SU(4) charge h = diag(1,1,1,−3), tr h³ = −192 on 10, and the cubic trace
    over the kernel space is ±1152 for 𝒞, while it vanishes on Sym²(16). The repair is 𝒞̃(M) = JMΓ on every kernel. It preserves the kernel
    space, squares to −1, commutes with the Pati–Salam action and conjugates the left ten, giving
    (10̄,3,1) ⊕ (10,1,3) ⊕ (15,2,2) ⊕ 2(6,1,1) ⊕ (1,2,2), the conjugate of the all-left branching; with 𝒞̃ the cubic trace is 0. The labels
    will be updated in Theorem 7.17, its figure, the Introduction (p. 12), Chapter 14 (p. 368) and Table C.2 (p. 384), whose caption already
    writes (10̄,1,3) inside the 126. The real kernel space, the potentials, the Hessians and Φ₀ ∈ (10,1,3) are unaffected.
54. **Error. Lemma 7.11(iv) (pp. 208–209).** "Otherwise V is unbounded below, or flat along a ray, and has no isolated minimum" fails on the
    boundary min{a(1/8), a(1/4)} = 0 with μ² < 0: for μ² = −1, α = β = 0, V = N has the isolated minimum M = 0. On that boundary μ² < 0 gives
    the unique minimum M = 0, μ² = 0 flat rays and μ² > 0 unboundedness. The phase results with both a(q) > 0 are unaffected.
55. **Hypothesis. Theorem 7.14 (p. 210).** The unique Fock ground state needs E > 0: A = 0, E = 0 satisfies every printed hypothesis and all
    sixteen states have zero energy. The proof divides by E.
56. **Error. Theorem 7.42, table (p. 226).** The row v_R = 0 describes the generic point of that stratum. When |v_u| = |v_d| the stabiliser
    among the 21 generators enlarges from 16 to 18 (a diagonal SU(2)), and at the zero vacuum it is 21. The strata will be listed.
57. **Hypothesis. Theorem 7.44 (p. 227).** The formula |P_⊥t|² = 16n_Kv_u²v_d²/(v_u² + v_d²) needs v_u² + v_d² > 0 (its continuous extension
    is 0). At v_R = 0 the gauge tangent space has no Majorana component, so its part in span{a, e_u, e_d} is span{z}, not span{a, z}; the
    formula still holds there because ⟨t, a⟩ = 8v_R = 0.
58. **Wording. Theorem 7.32(i) (p. 220).** Π_ν has rank 32 on the 272-dimensional kernel space and rank 31 on the 271-dimensional tangent
    space of the sphere.
59. **Evidence. Proposition 7.22 (p. 215).** As the text states, the identification of the spectrum of ΦᵀΦ with that of S†S uses an explicit
    isomorphism of the original computation that the version 1.0.0 suite does not reproduce. The map and its adjoint and product identities
    will be published with an exact check. *Done:* a13-2 gives the map exactly (v ↦ (ε⊗S)v̄ in a Weyl basis of the Γ = +1 half).
60. **Evidence. Theorem 7.23(ii) (p. 216).** The reduction of |𝒟| ≤ (3/8)I₄ to one constant is exact (Schur's lemma), but the constant 3/8 is
    obtained numerically; an exact contraction will replace it. *Done:* a13-1 (Q G⁻¹ Q = (3/8)² G on an integer basis of the (20′,1) component).
61. **Wording. Theorem 7.12(b) (p. 209) and Proposition 7.20 (p. 214).** J − JΓ is called "the left-handed fermion number" in Theorem 7.12,
    while its charge 1 − Γ is "twice the left-handed fermion number" in Proposition 7.20; the latter is correct.
62. **Wording. Model assumptions (§§7.5–7.6).** The thresholds x < 27 and x_* ≈ 24.884, and the ratios M_N/M_{W_R} < 27^{1/4} ≈ 2.280
    (local) and ≈ 2.233 (conjectured global) for one dominant family, hold within the stated model: an O(272)-symmetric tree potential at the
    matching scale, the trace kinetic normalisation, a common gauge coupling, a universal Yukawa coupling, fermion number, and the
    leading-logarithm approximation. The invariant counts (18 quadratic and 686 quartic invariants, 12 and 402 with fermion number) show that
    symmetry does not force this small potential; its special coefficients are matching conditions, and their stability under running and
    finite corrections is not analysed.
63. **Wording. p. 216 and p. 214.** The value 4.7 × 10⁸ s is the free-neutron oscillation time inferred from the Super-Kamiokande search
    with neutrons bound in ¹⁶O, not a lifetime of bound neutrons. By the Nielsen identities the value of the effective potential at a
    stationary point is independent of the gauge fixing; the background coordinates and the raw Hessian are not.
64. **Wording. Notation in Chapter 7.** Several symbols carry more than one meaning: κ (the frame rate c/ℓ, the holomorphic-quartic
    coupling, the bidoublet values κ₁, κ₁₅, κ_{u,d}); Q (the quartic invariant Tr(MᵀM)², the electric charge, the charge of Theorem 7.14);
    H (a Hamiltonian and the colour generator, both in Theorem 7.34); T (the mass form e₁e₂e₃ and T(û) in §7.6); K (kinetic kernels,
    bidoublet kernels K_{1u}, K_{1d}, frame boosts).

## Chapter 8 — Gauge Geometry and Photon Identification

76. **Error. Proposition 8.25 (p. 245).** V_𝔥 = Σ_a g_a²‖[t_a, Π]‖² is invariant under the group generated by 𝔥 only when the couplings form
    an Ad(H)-invariant quadratic form: one coupling on each simple ideal (the abelian centre may carry its own). For arbitrary positive g_a the
    claim fails: with 𝔥 = so⟨E₂,E₃,E₈⟩, t = ½(E₂E₃, E₂E₈, E₃E₈) and g = (1, 2, 3), the H-related planes ⟨E₂,E₇,E₁₃,E₁₄⟩ and ⟨E₃,E₇,E₁₃,E₁₄⟩ give
    V_𝔥 = 640 and 1280. The positivity of M² and the kernel statement are correct, and the Standard-Model and Pati–Salam calculations use one
    coupling per ideal.
77. **Error. pp. 246 and 262.** "The only other critical set of V_𝔥SM" is false. The plane P_e = ⟨E₂,E₃,E₉,E₁₄⟩, the θ = π/2 member of the
    neutral deformation of Proposition 8.49, is a critical point of the SU(2)_L term and of the hypercharge term separately, hence of V_𝔥SM for
    all positive couplings. It leaves su(3)_C ⊕ u(1)_em unbroken; with unit raw couplings and f = 1 its electroweak mass matrix has eigenvalues
    0, 32, 32, 64 (massive W and Z, massless photon). It is a saddle: relative to Tr(dΠ²) its Hessian has eigenvalues −2 (once), 0 (nine
    times), 1 (six times), and in the normalisation of §8.4 V_𝔥SM(Π_θ) = 128c_SM sin²θ is maximal at θ = π/2. The conclusion that alignment
    produces no stable electroweak-breaking vacuum stands; the claim of a complete list of critical sets is withdrawn (random-start
    minimisation cannot certify the absence of saddles).
78. **Error. Proposition 8.28 (p. 246) and the conventions of §§8.3–8.4.** The printed eigenvalues ½(g_Y² − g₂²) (×8) and 2g_Y² (×2) hold for
    the generators t^L_i and t_Y, whose norms in −Tr(XY)/128 are ½ and 5/6, not for the orthonormal basis that §8.3 declares. In that basis
    (√2 t^L_i, √(6/5) t_Y) the eigenvalues are 0 (×6), (3/5)g_Y² − g₂² (×8) and (12/5)g_Y² (×2), and the metastability threshold is
    g_Y/g₂ > √(5/3); in the form −Tr(XY)/512 of §8.4 the entries gain a factor four and the threshold is unchanged. Equal normalised couplings
    remain outside the metastable region. One convention for generators, couplings, the kinetic scale f (M²_phys = f²M²) and Hessians will be
    used throughout, with a conversion table; the two invariant forms of §§8.3 and 8.4, which differ by a factor four, will be unified.
79. **Wording. Proposition 8.47 (p. 255).** The quoted "Hessian coefficient" 128ξ₁c_SM + h is the coefficient of |x|² in V; the coordinate
    Hessian is twice it, and canonical scalar masses are the Hessian divided by the kinetic metric 128f². Signs and conclusions are unchanged.
80. **Wording. After Theorem 8.10 (p. 239) and Theorem 8.54 (p. 259).** "[B, Π] stays among the orthogonal projectors exactly when B is
    compact" holds modulo the stabiliser: the condition is [B + Bᵀ, Π] = 0. At Π_A the boost E₁₃E₁₀/2 commutes with the projector, so
    E₂E₇/2 + E₁₃E₁₀/2 is non-compact and gives the same tangent as E₂E₇/2 (see item 68). The dimension counts are unaffected.
81. **Wording. Remark 8.34 (p. 249).** A Stückelberg mass does not by itself cancel the anomalous Ward identity of the normal U(1); decoupling
    requires a specified Green–Schwarz or Wess–Zumino completion with its transformations and coefficients. Cancellation by the mirror block
    is one available completion, not a statement excluding other anomaly-cancelling matter.
82. **Hypothesis. Proposition 8.1 (p. 231).** 𝒫_r is non-empty only for ranks r divisible by four.
83. **Wording. Theorem 8.48 (p. 256).** "Global" is relative to the sixteen-dimensional family F₃, not to the manifold of all
    Lorentz-compatible rank-256 projectors, Gr_ℂ(64,128), of real dimension 8192. At ξ₁ = 0 the equality wording of the proof should be
    tightened; the conclusion holds because the overlap maxima already fix the points.
84. **Evidence. `verification/ch08_gauge_photon/p2_alignment.py`.** The header promises "exits 0 iff all pass", but the script has no final
    exit statement, so a failing check would appear only as a FAIL line; `run_all.py` counts FAIL lines and reports such a mismatch, and the
    version 1.0.0 record (20 PASS, no FAIL) is unaffected. An exit on the aggregate result will be added. Check p2-3 samples couplings that
    are constant on each simple ideal, and the minimisations of p2-16 cannot establish a complete list of stationary points.
85. **Wording. Notation in Chapter 8.** Ω denotes the interior volume and the connection Ω = Ω_Lorentz + A^aT_a; Q the second-rung charge ½ab,
    the electric charge and Q_A, Q_c, Q_v; H the second-rung vacuum blade and H_C, H(p); 𝒞, C and 𝖢 the Lorentz commutant and the current 𝒞_μ,
    the colour group and an element of the Klein group; D the duality map and the covariant derivative; T the tangents T_i, T_aw and the
    trivector.

## Chapter 9 — Inter-Constituent Dynamics and the b-Field

86. **Error. Proposition 9.22 (pp. 279–280).** With the nonnegative rates allowed in §9.2, item (i) fails at (κ₁, κ₂, β′) = (0, 1, 0), which
    is below the threshold although the first copy is free and M is not semisimple, and item (ii) fails at (1, 1, 0), which is at the
    threshold although the two copies are uncoupled oscillators and M is semisimple. Item (i) needs κ₁κ₂ + β′² > 0, item (ii) needs β′ ≠ 0,
    and the statement that the preserved positive form is not E₆ also needs β′ ≠ 0. The eigenvalue formula is correct.
87. **Error. Remark 9.57(a) (p. 310).** "The two factors coincide exactly in the isotropic case" is false: CCᵀ = CᵀC holds exactly when C is
    normal (after identifying the two triples), for instance C = diag(2, 3, 4). Isotropy is sufficient, not necessary. The orbit and energy
    statements do not use the converse.
88. **Error. Proposition 9.60(ii) (p. 312).** "Rotating (m₁, m₂) changes the spectrum unless sin 2θ = 0" holds for the initial pair (M, 0)
    only. For (m₁, m₂) = (2, 1) and cos θ = 4/5 the rotated pair is (1, 2), with the same spectrum, although sin 2θ = 24/25. In general the
    spectrum is unchanged exactly when the unordered pair (|m₁|, |m₂|) is. The conclusion that no continuous family at fixed spectrum changes
    (|m₁|, |m₂|) stands.
89. **Wording. Theorem 9.6(iv) (p. 268).** For the commuting boost the exponents are ±β ± iκ and also ±iκ, from the factor λ² + κ² of (ii).
90. **Wording. Theorem 9.9 (p. 270).** The necessary condition a·d = 0 has codimension one in the six-dimensional space of compact
    connections, so bounded motion is not robust under compact perturbations: for κ = 1, a = (0, 0, 2), d = (1, 0, ε),
    tr(η(S†)²) = −8a·d = −16ε, and M is not semisimple at ε = 1/10 although S is real antisymmetric. A physical model must explain what
    restricts the connection to this set.
91. **Evidence. Proposition 9.10 (p. 270).** The ten-dimensional bounded set and its singular locus rest on finite-difference Jacobian ranks,
    and rank loss of a defining map does not by itself make its zero set singular (x² = 0). These statements are numerical observations until an
    exact rank and local-geometry argument is given.
92. **Hypothesis. Theorem 9.49 (p. 300).** In (i), "dF is onto" needs nonzero momentum (at p = 0, dF = 0); (ii) and (iii) divide by κ₀² and
    need h(0) ≠ 0; in (iv) a frame with N_k p^{(k)}(0) = 0 is stationary, not unbounded. The secular-instability result stands.
93. **Wording. Theorem 9.27 (p. 284) and the copy field.** Nonnegative rates do not make H positive: for the bounded seesaw example, rates
    (0, 3, 5) and couplings (1, ½, ⅓), H has the eigenvalue −0.326115, so the empty state is not the lowest state of dΓ(H), which is bounded
    below with its negative modes filled. This exterior-algebra model is distinct from the canonical bosonic system of Theorem 9.12 and is not
    its quantisation. The removability of b refers to removing it while keeping the prescribed constant binding tensor.
94. **Wording. Proposition 9.36(ii) (p. 290).** "Dimension 8 for (1, 2, 3, 4)" counts the space of resonant coefficients; the Jacobi identity
    (ad = bc, item (v)) leaves a bracket variety of real dimension 6.
95. **Wording. Notation in Chapter 9.** b denotes an index and the connection; S_n the spinor module and, in Theorem 9.62, the spatial volume;
    N the number operator, the inertias N_k and a nilpotent part; H the Noether energy, the copy Hamiltonian and the operator 𝖧; C the
    conversion, the copy rotations C_±, C^{kl}, Ĉ, C_v, C_b and a complex symmetric matrix; M the phase-space generator and masses.
173. **Error. Theorem 9.40(iii) (p. 294) and Example 5.59 (p. 172).** M_pair, half a sum of singular values, is |sin φ|, not sin φ: at
     φ = −π/2 the singular values are 1 and 1. The chain E₆ ≥ M_pair ≥ ΣE₀ of Example 5.59 holds with |sin φ|.
174. **Wording. Proposition 9.51(i) (p. 301).** E₀ and |Q_kl| are conserved invariant functions, smooth only where they do not vanish (along
     z(t) = (1, i+t, 0), E₀ = |t|√(t²+4)/2); the globally smooth Casimirs are E₀² and |Q_kl|², as in item 27 for Theorem 4.40.
175. **Error. Remark 9.63(b) (p. 315).** "All elements of so(3)₊ are conjugate" is false (J and 2J are not); all Cartan directions, or all
     elements of a fixed nonzero norm, are. The conclusion that the algebra selects no Cartan direction stands.
176. **Wording. Table 9.3 (p. 292).** The rotor column gives the real dimension of the space of resonant coefficients before the Jacobi
     identity; the entries at n = 5, 6 are not dimensions of Lie-bracket varieties.
177. **Wording. Chapter 9 opening, §9.4 summary and Table 9.3 (pp. 265, 291–292); also the subsection "The copy rotation moves one quantum" of Chapter 5 (p. 167).** "Physical exactly between unequal
     rates" means non-removable while the prescribed binding tensor is kept fixed; a spacetime gauge action, observables and physical states
     are further structure.
178. **Wording. §9.3.5 (p. 281).** The plaquette is the holonomy of the copy connection around a triangle in copy space, the copy-space
     analogue of a Wilson loop, not the Wilson loop of a spacetime connection, as §9.4 already says.

## Chapter 10 — Classical and Astrophysical Consistency

96. **Error. "The analogy with the tower" (p. 325).** The paragraph says that every element outside the centraliser of K_n, "every boost
    included", creates quanta in pairs, and that "one algebraic dichotomy, compact against non-compact, decides both questions". The dichotomy
    of Lemma 5.46 is commuting against anticommuting with K_n, not compact against non-compact: the compact charge J = e₂e₃ creates a pair from
    the Fock vacuum, while the pole operators ωᵢ are boosts that commute with N. The paragraph will state the commutator criterion; the warning
    that tower quanta are not electrons stands.
97. **Error. Atomic clocks (p. 329).** With the bound Aω < 3.7×10⁻¹⁶ yr⁻¹ as printed, the amplitudes are A < 5.1×10⁻⁶ for ω equal to the
    inverse age and A < 3.7×10⁻¹⁶ for ω = 1 yr⁻¹. The printed values 2.6×10⁻⁶ and 1.8×10⁻¹⁶ are those for the libration frequency Ω = ω/2.
    Note also that ω = 1 yr⁻¹ is a period of about 6.3 years, which is not long compared with the years of comparison.
98. **Hypothesis. Clock and quasar constraints (p. 329).** A bound on a drift constrains Aω only under an assumption on the phase: for
    R = 1 − A + A cos(ωt + φ) the slope vanishes at an extremum. For two epochs, R(t₂) − R(t₁) = −2A sin(ω(t₁+t₂)/2 + φ) sin(ω(t₂−t₁)/2),
    which vanishes whenever the separation is a whole number of periods. The figures are order-of-magnitude sensitivities; a constraint on an
    oscillation needs a frequency-dependent analysis of a time series.
99. **Hypothesis. Birefringence rate (p. 328).** Dividing the angle by the age assumes the unwrapped, small-angle branch: a rotation by β and
    by β + nπ leaves linear polarisation unchanged, so for β = 0.35° the first alias is κ = 7.2×10⁻¹⁸ s⁻¹, about 290 times the quoted bound.
    The bound refers to that branch; excluding the others needs the visibility function, reionisation and the observation window. A
    time-dependent duality rotation is moreover not a symmetry of the free Maxwell equations, so the coupling hypothesis needs a stated
    propagation law.
100. **Wording. Proton stability (p. 327).** "A weak decay with MeV energy release has a lifetime of minutes unless it is forbidden" is not a
     general statement: lifetimes depend on matrix elements, selection rules and couplings. The argument is kinematic and motivates, but does
     not derive, the one-clock selection.
101. **Wording. Rainville et al. (p. 326).** The value (−1.4 ± 4.4)×10⁻⁷ combines silicon and sulfur; silicon alone gives (−9.7 ± 8.0)×10⁻⁷.
     The angle derived from it is a sensitivity estimate, not a silicon-specific bound.
102. **Wording. Smeared fields (p. 323).** Components smeared with different test functions "need not commute" (concentric radial test
     functions can commute); the smeared quantum distribution includes vacuum fluctuations and is not by itself a statement about a renormalised
     local rest-energy observable.
103. **Wording. p. 324.** Pair-production thresholds are not "the edges of the gamma-ray horizon": the horizon is an optical-depth condition
     that also needs the background spectrum, the cross section, the distance and cosmological evolution.
104. **Wording. Closing paragraph (p. 331).** A length and a coupling would not by themselves make the angles and the rate predicted; the
     dynamics and the states or boundary conditions that fix amplitudes, phases and clock orientations are needed as well.
105. **Evidence. Check e7-26.** The archived output of e7-26 still compares the local galactic magnetic energy density with the cosmological
     mean, the comparison that the text replaced by the local-with-local comparison of z5-1. The output will be corrected or marked as
     superseded.
106. **Attribution. Proposition 10.3 (p. 323).** The law f(t) = 2t is the eigenvalue law of a real Wishart matrix (N = 2, K = 3): Życzkowski and
     Sommers, *J. Phys. A* 34 (2001) 7111, eq. (4.5). It should be presented as inherited; the electromagnetic reading is the arena's.

## Chapter 11 — Curvature, Geometry and Gravitational Obstructions

107. **Error. p. 339.** "In (11.1) every mixed Christoffel symbol vanishes on Σ" is false: on Σ, Γ^u_{zv} = A_z and Γ^v_{zu} = −A_z (for A = bx dz,
     Γ^u_{zv} = bx); these are the normal connection coefficients. What vanishes are the symbols that mix E with N, so the conclusion that the
     splitting TM = E ⊕ N is Levi-Civita parallel along Σ stands.
108. **Wording. Proposition 11.14 and Figure 11.1 (pp. 338–339).** At the pure-Holst points the vacuum equations are not "empty": no equation for the
     metric or coframe remains after the connection is eliminated, but the connection variation still gives DΣ = 0 and, for an invertible tetrad,
     vanishing torsion.
109. **Hypothesis. p. 339 and Table 11.1.** That the phase enters an effective four-fermion interaction assumes minimal coupling of the fermions
     (Perez and Rovelli 2006); with a non-minimal coupling the Holst term can combine into the Nieh–Yan term and γ can drop out (Mercuri 2006).
110. **Wording. Theorem 11.13 (p. 337).** The density normalisation and the variables varied should be stated. The theorem concerns the displayed
     operator H(θ) = Ric − ½ sin²θ R g: since tr H = cos 2θ R, for cos 2θ ≠ 0 the equations H = 0 and Ric = 0 have the same solutions, and at
     sin²θ = ½ the equation is the traceless (unimodular) one. It does not exclude other actions, constrained variations or equivalent
     formulations of the vacuum equations.
111. **Wording. §11.1 (p. 332).** "The para-complex structure need not be posited": 𝒦 is identified with the grading root once a pairing of
     positive and negative directions is chosen; the pairing is not canonical given the neutral metric alone, and a global construction needs the
     bundle reductions and spin assumptions stated.
112. **Wording. Theorem 11.10 and Corollary 11.12.** ℛ ↦ e^{θG}ℛe^{θG} is a congruence, not a similarity, and is not a map between metrics or
     solutions; the gradient of the Gauss–Bonnet density is taken with respect to the pairing Tr(AB).
113. **Wording. Theorem 11.5, Corollary 11.6 and §11.4.** The wrong sign is relative to the adopted Einstein–Hilbert action and the stated
     connection-metric ansatz, and integrating over a finite normal disc is not a consistent compactification (boundary terms and the
     discarded equations are not treated). These hypotheses should appear in every summary of the result.
114. **Wording. Example 11.8 (p. 335).** The radii should be written with |b| and |k| (or b > 0 fixed by orientation).
115. **Wording. Example 11.17 (pp. 340–341).** A Maxwell term added by hand does not cure the wrong-sign term obtained from the ambient
     Einstein–Hilbert action, and the positive norm of the angle mode does not make the projector sector stable (its trace form has signature
     (3,3)). The two Stückelberg combinations of the chapter use opposite conventions, B = dχ − A with χ ↦ χ + λ (p. 335) and b = dθ + A with
     θ ↦ θ − λ (p. 340); one convention, or the relation θ = −χ, should be stated.
116. **Wording. Notation in Chapter 11.** B denotes the magnetic constant, the Stückelberg combination, a generic bivector and ½J_N; G the Hodge
     star, the shift of Theorem 4.2 and G_N; Ω the para-Hermitian two-form and the spin connection; 𝒱 the volume element and the space of
     curvature operators; T a tangent vector and the stress tensor; Q the Gauss–Bonnet functional.

## Chapter 12 — Three Temporal Directions and Causal Reduction

117. **Error. Theorem 12.12 (p. 350) and Proposition 12.16 (p. 352).** With the action ∫dτ (p·ẋ − NH − λ₁p_u − λ₂p_v), varying N gives H = 0, not
     the mass shell 2H = m². The action should read ∫dτ (p·ẋ − N(H − m²/2) − λ₁p_u − λ₂p_v); the brackets and the gauge symmetry are unchanged.
     Then dτ_h = |N|m dτ, the gauge N = 1/m needs m > 0, and reparametrisations are generated on shell by the total first-class Hamiltonian,
     including the multipliers. The check c2-8 inserts the mass shell instead of deriving it; a test of the N-variation will be added.
118. **Attribution. §12.6 (p. 353).** The statement that no proposal for a universal evolution parameter was found and that the concept is
     introduced here as an original proposal is withdrawn. An invariant evolution parameter is the subject of Horwitz and Piron (1973), cited
     in the same paragraph, and a universal parameter for many-body systems is stated in Horwitz, *Entropy* 22 (2020) 804. What may be specific
     here is its relation to the clock reduction of this framework.
119. **Error. §12.6, item (3).** If the action is not reparametrisation invariant, τ_U need not be a preferred time that breaks the Lorentz
     algebra or SO(3,3): ∫dτ_U (p·dX/dτ_U − η^{AB}p_Ap_B/2M) is SO(3,3)-invariant with τ_U a scalar. A self-adjoint generator gives unitary
     evolution in τ_U even when it is unbounded below, so unitarity and the existence of a ground state remain separate requirements.
120. **Hypothesis. Theorem 12.4(iii) (p. 346).** The centraliser of dimension seven requires N nondegenerate: for the totally null plane
     ⟨e₁+e₄, e₂+e₅⟩, C₀ = (e₁+e₄)(e₂+e₅) has C₀² = 0 and a nine-dimensional centraliser. The degenerate strata should be treated separately.
121. **Hypothesis. After Theorem 12.6 (p. 348).** The Dirac current at three constituents is null on product states whose first factor is
     chiral; for |000⟩ + |100⟩ in the first factor it is (2, 0, 0, 0), timelike. The check of the untwisted square in Theorem 12.6(iii) samples
     the whole module and will be restricted to S⁺.
122. **Wording. Closing paragraph (p. 355).** That a gauge mechanism "would make the would-be ghosts gauge … and remove the six flat ones" is a
     research direction, not a result: the linearised gauge variation [ε, J] has rank eight, and the six so(1,3) directions are not among the
     gauge directions. Likewise, the suggestion that positivity of the Hilbert-space norm fixes causal order is a conjecture.
123. **Wording. After Proposition 12.2 (p. 345).** "Causal order … is absent" should read that the split metric supplies no canonical future cone
     invariant under the connected group; a cone or order can be introduced by additional data.
124. **Wording. Proposition 12.1 (p. 344).** Plane waves and a sin(qu) are not of finite norm on the noncompact Cauchy surface; the ill-posedness
     argument should use wave packets concentrated at large surplus momentum, or a periodic box.
125. **Wording. Proposition 12.7, Theorem 12.9 and Proposition 12.10 (pp. 348–350).** "No invariant potential can prefer (1,3) to (3,1)" holds for
     invariant polynomials; smooth invariant functions can separate the two orientations, and the stronger obstruction is the orbit closure of
     Theorem 12.9. The six modes of Theorem 12.9 are non-Goldstone flat directions (their vanishing masses are enforced by invariance and
     stability). Proposition 12.10 concerns an unconstrained adjoint scalar with a factorised kinetic term, which should be stated.
126. **Hypothesis. Corollary 12.14 and Theorem 12.15 (p. 351).** The reduced Hilbert space and the choice of positive frequencies should be
     specified (constant functions of u, v are not normalisable). The vanishing of the mixed kernel follows from translation invariance in the
     flat product; the rotating normal metric with general A is not translation invariant and needs the reduced operator itself. Theorem 12.15
     concerns normal vectors obtained by the Einstein–Hilbert and Kaluza–Klein route, and its statement about spinless matter the
     zero-normal-momentum sector.
127. **Wording. Identifiers in the printed text (Chapters 3, 10, 12–16).** Codes such as CORE-07, X-03, TIME-08 and CLASS-08 (PH-07, PH-09 and
     CORE-13 in Chapters 14 and 16), and LAD-01,
     LAD-04, PH-02, z1-1d in Chapter 15, refer to an internal review register that is not published; they will be removed from the printed text
     or listed with one-line descriptions.
128. **Evidence. Verification READMEs.** The READMEs of `ch04_rigid_body`, `ch06_chiral_matter`, `ch11_curvature_gravity` and `ch12_causality`
     give paths that no longer exist (for example `verification/ch12_causal/…` for `verification/ch12_causality/…`, and helpers now in
     `verification/common/` and `verification/tools/`). The scripts and the version 1.0.0 record are unaffected. Tests will be added for the
     N-variation of item 117, the null centraliser of item 120 and a chiral sample at three constituents.
129. **Wording. Notation in Chapter 12.** In Theorem 12.8, f is both the degree-eight function and a coefficient of the proof, beside e, g, m;
     B denotes the order parameter, the Kaluza–Klein potentials and the magnetic field; C the selector and the constraint function; N the normal
     plane and the lapse; n the null rotation and a lattice vector; K = ½e₁e₄ beside 𝒦 and the kernels.

167. **Wording. Table 12.1 (p. 354), row "universal τ_U".** "No proposal" repeats the claim withdrawn in item 118. It should read: prior art
     (Horwitz–Piron); no mechanism in this arena; requirements listed.

## Chapter 13 — Compatibility of the Mathematical Constructions

130. **Error. §13.5 (pp. 363, 365–366) and Table 13.4.** "The minimum of E₆ over boosts, attained exactly in the rest frame" fails for nonzero null
     fields, which have no rest frame: for E = (1,0,0), B = (0,1,0), E₀ = 0 while a boost along the propagation direction gives E₆ = e^{−2r} > 0 for
     every finite rapidity. The correct statement is E₀ = inf over boosts of E₆, attained exactly at vanishing Poynting vector for non-null fields
     and not attained for nonzero null fields (Theorem 3.95 states this correctly). Likewise E₀ classifies the Lorentz-and-duality orbits of
     nonzero fields; the zero field is a separate orbit.
131. **Error. §13.5 (p. 364).** "A field point with E₀ ≠ 0 keeps only {±1} of [the duality group]": under duality alone a nonzero field has trivial
     stabiliser, since e^{πG}F = −F. The ℤ₂ is the stabiliser of the invariant q, or of the Lorentz orbit it labels, as the preceding paragraph says.
132. **Hypothesis. §13.5 (pp. 363–365).** U(3) is the group of transformations that preserve E₆ and commute with G (Theorem 3.91); the full linear
     symmetry group of E₆ is O(6). Within U(3) the real structure is preserved by O(3) and E₀ by U(1)·O(3).
133. **Hypothesis. §13.6 and Table 13.5 (pp. 366–367).** The connection between manifolds of unequal rates is not removable while the rate matrix
     is kept fixed; on an open interval of the worldline the gauge e^{−τb} removes it at the price of a time-dependent rate matrix, with
     d/dτ K′ = [b, K] at τ = 0. That the rates are constant follows from the fixed-rate model, not from gauge invariance; closed worldlines need
     a holonomy statement.
134. **Error. §13.3, item (v) (p. 362).** The transported projector field Π(x) of Chapter 11 acts on the spinor module ℝ⁸, the observational
     projection pr_n̂ on the vector space ℝ⁶; Π(x) is not its curved counterpart, and Chapter 11 shows that the two behave differently under
     transport. They will be listed separately.
135. **Wording. Figure 13.1 and Table 13.2 (p. 359).** The link "field-strength reading → charge and gauge geometry" marked established combines
     a representation theorem of Chapter 3 with the charge identification of Chapter 8 without a map or common action between them; it should be
     marked conditional, in line with the figure's stated rule of drawing the weakest status.
136. **Wording. §13.5 (pp. 364–367).** Schur's lemma fixes the isotropic form C = c·1 of the conversion, not the value of c; for the explicit
     Dirac mass, the three broken directions are directions of the mass orbit, not Goldstone excitations; and a linear copy z ↦ (z, z) of classical
     field data is not forbidden by the no-cloning theorem, which concerns unknown quantum states (item 39).
137. **Wording. Summaries in §§13.1–13.4.** The corrections of items 8, 22, 62, 65, 66 and 125 apply to their summaries here. In addition: §13.1.3
     says the repair at the second rung costs "the whole internal algebra" (p. 357), while §13.4 (p. 362) correctly keeps the rotation of one
     surplus pair; "the electroweak stage" means the analysis with a given bidoublet value (Proposition 7.35) on p. 357 and the missing breaking
     mechanism on pp. 357, 359 and 362; and the gap between the sufficient and the necessary conditions of Chapter 7 is supported, not resolved,
     by numerics.
138. **Wording. Tables 13.2 and 13.3.** "The gauge action (45)" refers to equation (7.4). In the field-octant row of Table 13.3 the complex
     structure entry "none" should say whether it means no additional choice, since the duality lift of Chapter 3 exists.

## Chapter 14 — Comparison with Related Geometric Approaches

139. **Error. §14.2 (p. 369).** An intertwiner between V_* and Krasnov's even forms does not require identifying a block of Cl(9,9) with Cl(7,7):
     Cl(7,7) ≅ M₁₂₈(ℝ) has no nonzero module of dimension below 128, and the 64-dimensional half-spin space carries only the even part
     Cl⁰(7,7) ≅ M₆₄(ℝ) ⊕ M₆₄(ℝ). What is needed is an invertible map intertwining the explicitly identified common generators, then a comparison
     of the kinetic forms. The absence of the non-compact so(4,6) generators on V_* obstructs that ambient extension only.
140. **Error. §14.2 (p. 369).** Dropping a positive four-plane from (9,9) leaves signature (5,9), not (7,7); an orthogonal complement of signature
     (7,7) requires removing a (2,2)-plane.
141. **Wording. §14.2 (p. 369).** SO(4) × SO(6) is the Pati–Salam group only as Lie-algebra shorthand: Spin(4) ≅ SU(2) × SU(2), Spin(6) ≅ SU(4), and
     on the family (4,2,1) ⊕ (4̄,1,2) the faithful group is (SU(4) × SU(2)_L × SU(2)_R)/ℤ₂.
142. **Attribution. §14.3 (pp. 369–370).** The selection of 𝒜_F in Chamseddine and Connes also uses the condition [D, Z(𝒜)] ≠ 0 and maximality
     among the admissible subalgebras of 𝒜^ev (Why the Standard Model, §4); the order-one condition alone is vacuous at D = 0.
143. **Wording. §14.4.** The novelty labels should separate the graded tensor-product identity, standard centraliser computations, the count
     2 dim_ℂ Sym²(ℂ¹⁶) = 272 and the energy–inertia identities (Kaiser 2011, item 21) from the explicit realisations, projectors, obstructions and
     coefficient-dependent selections; a universal evolution parameter keeps its Horwitz–Piron comparison (item 118).

162. **Wording. Table 14.1 (p. 370), row "one generation", and Remark 6.24 (p. 199).** "Same count; equivalence as modules not constructed"
     understates what is known. Under so(1,3) ⊕ so(4) ⊕ so(6), the block V_* (Corollary 6.12) and Krasnov's Λ^even ℝ⁷ (his §6) are both the real
     form of (2,1) ⊗ [(4,2,1) ⊕ (4̄,1,2)]. By complete reducibility they are equivalent representations of that algebra, once its two embeddings
     are identified compatibly with the labels. What is open is the explicit identification of the generators and conventions, the map itself,
     and whether it carries the kinetic structures. Remark 6.24 should say the same. The non-compact generators of so(4,6) obstruct only an
     extension to the ambient action (item 139).
163. **Wording. Table 14.1 (p. 370), row "kinetic term".** "Different" should read "kinetic equivalence not established". The opposite signs of
     Krasnov's eq. (219), which he states coincide with the Standard-Model kinetic term, and K⁰ = id on V_* can be compared only after common
     choices of complex structure, conjugation, Grassmann ordering and field variables. The displayed signs establish neither a ghost on one
     side nor its absence on the other.
164. **Attribution. §14.3 (p. 369).** The classification of Chamseddine and Connes makes the dimension per generation a square, k². Sixteen
     fermions (k = 4) do not follow from the symplectic–unitary hypothesis alone: they follow from the quaternionic hypothesis together with the
     choice of its simplest non-trivially graded case, of quaternionic dimension two (Why the Standard Model, §3). The finite algebra is, up to
     an automorphism of 𝒜^ev, the unique involutive subalgebra of maximal dimension that admits off-diagonal Dirac operators satisfying the
     order-one condition (ibid., Theorem 4.1; see item 142).
165. **Wording. §14.3 (p. 370), last sentence; also the paragraph after Theorem 7.34 (p. 223).** "A normalisation and not a prediction" conflates
     three steps:
     - normalising the hypercharge generator, T_Y = √(3/5) Y, which is a convention and gives g₁ = √(5/3) g_Y;
     - the equality g₁ = g₂ at a scale, which is a physical boundary condition;
     - the spectrum, thresholds and running.

     In the spectral action the boundary condition follows from the action's gauge kinetic coefficients, so it is a conditional prediction at
     the cut-off (Chamseddine, Connes and Marcolli 2007). In this manuscript it comes from the single gauge coupling of the adopted action
     (7.4), which is a choice. Chapter 7's phrase "this is a normalisation statement" should name that choice. Its conclusion, that no low-energy
     angle is predicted, stands.
166. **Attribution. §14.3.** The spectral Pati–Salam model of Chamseddine, Connes and van Suijlekom should be compared:
     - *Beyond the spectral Standard Model: emergence of Pati–Salam unification*, JHEP 11 (2013) 132;
     - the inner fluctuations without the first-order condition, arXiv:1304.7583;
     - the one-loop running of its couplings, JHEP 11 (2015) 011.

     It has the same family, (2,1,4) ⊕ (1,2,4), so the Pati–Salam content is shared prior art there as well. The distinguishing questions are
     realisation, selection and bosonic dynamics.

## Chapter 15 — AI-Assisted Research and Verification Methodology

144. **Wording. §15.2 (p. 373).** "Exits with status zero only if every check passes" is not true of `p2_alignment.py` until the fix of item 84.
145. **Wording. §15.2 (p. 374).** "Fifty referee findings … none rejected, and no theorem refuted" describes the internal review step of the rest-energy
     scripts; the external reviews later found false statements in the same chapters (for example items 22 and 86–88).
146. **Wording. §15.2.3 (p. 375).** "The review adjudications of Chapter 13": Chapter 13 contains no review adjudications; the reference should be
     replaced (see item 127).

## Chapter 16 — Conclusions and Open Mathematical Problems

147. **Error. §16.4 (p. 379).** sin²θ_W = 3/8 does not follow from the charge content: with the same charges, g_Y = g₂ gives 1/2. The value 3/8
     needs g_Y² = (3/5)g₂² (equal normalised couplings g₁ = g₂) at a stated scale; the traces Tr Y² = 10/3 and Tr (T³_L)² = 2 per family give the
     ratio only under a common kinetic normalisation. Chapter 14 states this correctly.
148. **Error. §16.3, item (5), and the next steps (p. 378).** An index is a signed count of zero modes, not a spectrum: on flat four-tori the index
     is zero while the eigenvalues depend on the radii. A Lorentzian Dirac operator is not elliptic, so the base, signature, bundles, domain and
     boundary conditions must be specified before an index theorem applies (for example Bär and Strohmaier 2019); the spectrum is a separate
     problem.
149. **Wording. §16.3, item (4), and the next steps (p. 378).** No finiteness of the relevant embeddings has been proved, so the "finite search"
     is not yet a decision procedure; and an algebraic answer would not by itself settle a shared physical theory, which also needs compatible
     dynamics.
150. **Wording. §16.2, items (5)–(7) (p. 377).** The four objections to the normal U(1) concern different models and should be attributed each
     to its own; the scalar no-go concerns an unconstrained adjoint scalar with an invariant bilinear kinetic term (item 125); ill-posedness concerns
     unrestricted ultrahyperbolic data (Craig and Weinstein treat constrained evolution); and the failure of the Ricci–DeTurck method is not the
     non-existence of every split-signature flow.
151. **Wording. §16.4 (p. 378).** "The arena is derived" should read that the arena is chosen and its consequences derived. A theorem need not
     remove a choice without adding another to be new mathematics; reducing the independent choices is the programme's aim, not a definition.
152. **Attribution. §16.1, items (2) and (6) (pp. 376–377), and §13.1.1 (p. 357).** The reductive dual pair and the identity
     C_cop + C_Δ = ½n(n+2) are inherited (Howe 1989; Neergård 2019, item 33) and should be marked so in these summaries; the explicit complex
     structure of the kernel space is the corrected 𝒞̃ of item 53.
179. **Wording. §13.1 and §16.1 (pp. 357, 376); Appendix A (p. 380).** "The force law closes in one way" holds within the stated equivariant
     one-generator ansatz. The opening of Appendix A should distinguish the integer realisation of the generators from parameter-dependent
     spectral, variational and numerical statements, as §I.1.2 does (item 154).

## Corrections from the reviews of the prepared version 1.0.1

These were found in independent reviews of the prepared version 1.0.1 (8–9 November 2026) and corrected before its release; where the same
text stood in version 1.0.0 it is affected as well. Statement numbers are those of both versions.

180. **Chapter 1.** After Corollary 1.4: an indefinite matter form removes the positive-metric obstruction, and whether ghosts occur depends on
     the action and its constraint reduction; gauge-admissibility is defined relative to the matter metric H (XᵀH + HX = 0). Theorem 1.3 is
     retitled "The signature rule for coordinate bivectors"; three paired planes (eᵢ, eᵢ₊₃), not six; the grading root is defined in
     Proposition 1.12; the blade identities concern coordinate blades.
181. **Chapter 2.** Probability weights alone do not fix the Fisher–Rao metric: compatibility with congruent Markov embeddings is a stated
     requirement (Ay, Jost, Lê and Schwachhöfer 2015). Theorem 2.19 identifies the algebras; the Maxwell–Dirac and QED statements of
     Corollary 2.20 transport the standard action and quantisation (Hestenes 1967). Moretti and Oppio's result has its own hypotheses
     (Remark 2.18). The involution S induces the Cartan involution.
182. **Chapter 3.** Table 3.1: the interior oscillator action is O(3,3)-invariant up to a boundary term. Figure 3.2: I₄ = −Ω = ⋆, not Ω = ⋆.
     Theorem 3.95(iii): |v_rest| = |β|, and a null field has no timelike rest observer. Proposition 3.31(i) needs E ≠ 0; Theorem 3.70(i) needs
     two linearly independent clocks. Summaries: generators of boosts and duality are odd, finite transformations of mixed parity. The
     higher-derivative ghost is the Podolsky-type extension (Podolsky 1942). Proposition 3.100(ii): an infinite plane wave has no finite
     integrated four-momentum. The classical/quantum status of the Maxwell–Dirac sector and of the axial rotation in Corollary 3.59 is stated.
183. **Chapter 5.** The raw duplication z ↦ (z,z) doubles E₆; z ↦ (z,z)/√2 divides it. Table 5.9: constant *nonzero* rate. An invariant number
     operator exists exactly at even rungs *under the common binding action* (Theorem 5.50(iii) gives another action). Interior boosts also
     mix sectors differing by two quanta (Proposition 5.47). Neergård's parameters are d_N = n, k_N = d/2. "No joint mass shell" means no
     simultaneous nonzero-mass graded Dirac constraints; zero masses and zero lines are treated separately (Propositions 5.34, 5.63);
     Corollary 5.2 is scoped to the grade-two connection algebra; centralisers on vectors versus the O(n) duality on spinors; T(z₁)+T(z₂) is
     not T(z₁+z₂) (Example 5.58); Theorem 5.57(ii),(iv): priority unestablished.
184. **Chapter 6.** The charge table of §6.3.2 uses ε = +1; Remark 6.24: underlying real representation (realification); Theorem 6.20 concerns
     the perturbative hypercharge-anomaly variety, of dimension 6N_g − 4 before an overall normalisation.
185. **Chapter 7.** Theorem 7.42 lacked the stratum v_R ≠ 0, v_u = v_d = 0 (stabiliser 12 in 21, 14 in 23, nine massive vectors; a11-3), and
     Proposition 7.35(iv) needs v_u² + v_d² > 0. C = e₂ fixes Π_W only when e₂ ∈ W (§7.3.2), and Table 7.3 is a selection of representatives.
     Proposition 7.40: an anomaly allows a pseudo-Goldstone mass but does not compute it. The first loss of positivity is the left-handed
     triplet at x = 27, not "the lowest line"; zero-quartic boundary in the summary of Lemma 7.11; v = √2|r|; Fock-vacuum summaries need E > 0.
186. **Chapter 8.** Proposition 8.27 and the summaries need a nonempty fixed set of the gauged subalgebra. Remark 8.29: the cutoff term is
     3Λ²f²V_𝔥/(32π²). Theorem 8.54: the eight charged directions are the image of ad_Q on the orbit tangent space. Anomaly cancellation in
     the larger spectator block is one completion, not the only one; the normal-connection field couples to half the fermion number; the
     extra critical set leaves a massless neutral vector; the coefficient 128c_SM sin²θ is in the conventions of §8.4.
187. **Chapter 10.** The birefringence bound holds on the unwrapped small-angle branch of the stated propagation hypothesis; the aliases
     β + nπ (first alias 7.2 × 10⁻¹⁸ s⁻¹) are not excluded by dividing the angle by the age, so the exclusion of particle-linked rates is
     restricted to that branch (synopsis, §10.1.5, Table 10.2, answer). Uphill kinematics motivates the imposed common clock but bounds no
     proton lifetime. The time-odd invariant is scoped to the three-copy interior model. u₀ is the infimum over frames, attained for non-null
     fields. Table 10.2 labels the pair row by number conservation, [X,N] = 0.
188. **Chapter 11.** The vacuum count of general relativity needs sin 2θ ≠ 0 (an invertible tetrad and a nonzero Palatini coefficient); at the
     pure-Holst points the coframe equation is empty (§11.4(ii)). Spinless matter decouples from the normal connection only in the
     zero-normal-momentum reduced sector (Table 11.1, §11.1). The dichotomy after §11.3 is left open, not resolved by the torsion discussion.
189. **Chapter 12.** Corollary 12.14 states its background and quantisation: flat h with a time orientation and m² ≥ 0 for the Fock statements;
     for globally hyperbolic curved h only the microcausality of the algebraic construction follows (Fewster and Verch 2012). The proper-time
     gauge of Proposition 12.16 is N = ±1/m on a segment with nonvanishing lapse. The no-go summary of §12.5.1 is scoped to an unconstrained
     adjoint order parameter with an invariant potential and the factorised kinetic term. Theorem 12.4(ii): the degenerate case follows from (i).
190. **Introduction, Chapters 13, 15 and 16.** Chapter 13 repeated the duplication error of item 183 (z ↦ (z,z) doubles the energy; z ↦ (az,bz)
     with |a|² + |b|² = 1 preserves it) and said "the algebra fixes c" (only the isotropic form of the conversion map is fixed). Figure 13.2
     put the faithful 23-dimensional algebra under "anomaly-free gauge action" (the anomaly-free gauged algebra is Pati–Salam, dimension 21).
     The restricted tensor-potential selection, the full 272-mode local result and the conjectural global minimality are separated. The
     qualifications of items 186, 188 and 189 (reduced-sector decoupling, flat h, sin 2θ ≠ 0, nonzero linear intertwiner, closure ansatz,
     fixed-rate model, Ricci–DeTurck method) are carried into the figure and the summaries of Chapters 13 and 16; a Lorentzian four-plane
     corresponds to a normalised simple bivector with an orientation; the classification "fixed by the algebra" reads "proved from explicitly
     stated hypotheses"; a decision procedure needs an effective enumeration, not finiteness alone; the run command of Chapter 15 is typeset
     with literal hyphens.
191. **Chapter 14.** §14.6(iv): the colour obstruction concerns the effective bivector algebras of the axial members of the five-projector
     family (Corollary 6.6); P₁ carries colour (Theorem 6.7). An intertwiner for the common Lorentz–Pati–Salam algebra exists after compatible
     identification; what is open is an explicit map and the comparison of kinetic forms (§14.3, §14.6(v), Table 14.1); "realification" for
     "real form". §14.4: the fermion count is k², the finite Hilbert space has complex dimension 2k² (32 for k = 4). The birefringence entry
     points to §10.1.5; the vacuum-selection status separates the conjecture of Chapter 7 from the conditional global result of Theorem 8.48;
     the worldline gauge removal is an input of the selected model, not "rather than by fiat".

## Related work to be added (all chapters)

The comparison of Chapter 14 will be extended to: two-time physics (I. Bars and collaborators); an SO(3,3) BF gauge theory on a
six-dimensional spacetime of signature (3,3) (arXiv:2602.19151, 2026); three-dimensional time (quant-ph/0510010); real Clifford algebras
and spinors in all signatures (arXiv:1908.02235); and Atiyah, Bott and Shapiro, *Clifford modules* (1964). Hestenes, *Real spinor fields* (1967), already cited, will be
compared directly with Chapter 2. For Chapter 1: Gualtieri, *Ann. of Math.* 174 (2011) 75, for split spinors as forms and the half-density factor, and Krasnov (2018, §4) for SO(3,3) with its Lorentz embedding and extra SO(2). For Chapter 5: Neergård (2019, 2020) and Howe (1989) on orthogonal dual pairs and their Casimir
relations; the Jordan–Wigner construction for the chirality strings; Gaillard and Zumino (1981) on duality rotations of several Maxwell
fields; and Kaiser (2011) and Thompson (1979, singular values and diagonal entries of complex symmetric matrices) for the pair invariant
of Theorem 5.57, whose priority is not yet established. For Chapter 6: the constructions of one chiral family by Maraner (2004), Nesti and
Percacci (2010) and Krasnov (2018, §§6.4–6.7), compared in Chapter 14, will be cited where the Pati–Salam block is constructed
(Theorems 6.10–6.13), with a statement of what this realisation adds; Remark 6.24 distinguishes the ambient actions only. For Chapter 7: Bertolini, Di Luzio and Malinsky (2010, 2012) on the quantum vacuum of
nonsupersymmetric SO(10), as the comparison for the radiative vacuum analysis; Konaka et al. (1986) and Peccei, Wu and Yanagida (1986) for the
exclusion of an electroweak-scale axion. For Chapter 8: Coleman, Wess and Zumino, and Callan, Coleman, Wess and Zumino (1969), for
sigma models on homogeneous spaces and their gauging; Preskill (1991) on anomalies and massive anomalous gauge fields. For Chapter 9: Flanders and Wimmer (1977) for an invariant proof of Roth's theorem; Bravyi (2005) on
quadratic fermion dynamics; Zusmanovich (2010) on current Lie algebras, which contain the truncated graded brackets of Proposition 9.36; and
Alekseevsky and Cortés (1997) on invariant bilinear forms of spinor modules (Theorem 9.62). For Chapter 10: Fedderke, Graham and Rajendran (2019) on oscillating
birefringence and its washout; Van Tilburg et al. (2015) on searches for oscillations with atomic clocks; and the ACT DR6 analysis of
Diego-Palazuelos and Komatsu (2026), 0.215° ± 0.074°, with its authors' caution on systematics. For Chapter 11: Mercuri (2006) on non-minimal fermion coupling in the Holst
action; Henneaux and Teitelboim (2005) on duality-invariant linearised gravity and (1989) on the unimodular formulation. For Chapter 12: Penrose
(1960) on spinor coordinates of the null cone; Richardson and Slodowy (1990) and Jablonski (2012) for real orbits and invariants; Bars, Deliduman
and Andreev (1998) on gauge constraints in two-time physics; Horwitz (2020) on a universal evolution parameter; and Dartora and Cabrera (2010),
a split-signature six-dimensional Dirac model with local selection of a time direction, as prior art for the general architecture. For Chapters 13–16: Chamseddine, Connes and Marcolli (2007) for the coupling relation at the
cut-off of the spectral action, derived there from the action; Chamseddine, Connes and van Suijlekom (2013, 2015) for the spectral Pati–Salam model; Bär and Strohmaier (2019) for an index theorem on Lorentzian manifolds; Kaiser (2011), item 21. For Chapter 12 and the Introduction: Feynman, *Rev. Mod. Phys.* 20 (1948) 367, for the
sum-over-histories motivation of the causal programme (motivation only), and Krasnov (2018) retained in the Introduction.
