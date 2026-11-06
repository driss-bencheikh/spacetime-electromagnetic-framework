# Errata for version 1.0.0

Known issues in `release/manuscript-v1.0.0.pdf` (DOI [10.5281/zenodo.23900465](https://doi.org/10.5281/zenodo.23900465)), found by
independent AI-assisted reviews of Chapters 1 to 9, and by an AI-assisted chapter-by-chapter
derivation check, and checked by the author. They will be corrected in version 1.0.1. Page numbers are
the printed page numbers.

**None of these issues invalidates a main theorem.** Each correction is local. Where a later chapter repeats an affected statement or cites
it as justification, the item says so (items 53, 65, 68 and 77); in every such case the later results hold once the correction is made.
Otherwise a later chapter that uses an affected result uses it in a case where the result holds, or already states the missing hypothesis.

Items are numbered in the order in which they were recorded and listed by chapter. Items are marked **Error** (a statement false as written), **Hypothesis** (a missing condition or domain restriction), **Wording**
(imprecise or misleading phrasing), **Attribution** (a missing or incomplete reference) or **Evidence** (a proof or check weaker
than the statement it supports).

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
    which are the only ones used later, stand.
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
    will be published with an exact check.
60. **Evidence. Theorem 7.23(ii) (p. 216).** The reduction of |𝒟| ≤ (3/8)I₄ to one constant is exact (Schur's lemma), but the constant 3/8 is
    obtained numerically; an exact contraction will replace it.
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
Alekseevsky and Cortés (1997) on invariant bilinear forms of spinor modules (Theorem 9.62).
