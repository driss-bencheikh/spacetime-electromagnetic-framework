# Referee-style review of version 1.0.1

**Manuscript:** *The Spacetime–Electromagnetic Framework: Mathematical Structures and Obstructions in Split-Signature Clifford
Geometry*, version 1.0.1 (8 November 2026), `release/manuscript-v1.0.1.pdf`, 428 pages, built from `manuscript/main.tex` at commit
`75af386`. The SHA-256 of the PDF matches the value in `release/README.md`.

**Status of this document.** This is an AI-assisted review written in the form of a referee report. It is not formal peer review and
certifies nothing; in the vocabulary of Chapter 15 it is evidence class (vi), a source of objections to be adjudicated. Every factual
statement below about the manuscript or the repository can be checked against the cited file, line or command.

**What was examined.**

- Read in full: abstract, Introduction, Chapters 1, 2, 6, 15 and 16, Appendix D, §§7.5.1–7.5.2 (the tensor model and the selection theorems),
  §14.6, and the release, status and errata documents.
- Read at the level of every formal statement (theorem, proposition, lemma, corollary, definition, assumption), with proofs read
  selectively: Chapters 3, 4, 8, 10, 11, 12, Chapter 7 from §7.4, and about two thirds of Chapter 5.
- Read at the level of statement titles only: Chapter 9 and the last third of Chapter 5. Chapters 13 and 14 were read in part.
- Recomputed independently, in code written for this review from equations (1.1)–(1.2) only, without importing anything from
  `verification/`: the items of §4 below.
- Re-ran the whole verification suite (`verification/run_all.py -j 3`): result in §5.
- Not done: a line-by-line check of all proofs; a rebuild of the PDF (no TeX installation was available); a systematic literature
  search for priority.

---

## 1. Summary of the submission

The manuscript studies the real Clifford algebra Cl(3,3) ≅ M₈(ℝ) of ℝ^{3,3} and its graded powers Cl(3n,3n), and asks which ingredients
of relativistic field theory this algebra fixes, admits or obstructs. Its main threads are:

1. a chosen Lorentzian four-plane in ℝ^{3,3}, whose orthogonal positive two-plane supplies a complex structure J, so that the real
   spinor module with J is ordinary complex Dirac theory (Chapters 1–2);
2. a second, inequivalent Lorentz algebra, the centraliser of a complex structure G, under which ℝ^{3,3} is Λ²ℝ^{1,3} with G the Hodge
   star, and a large body of Maxwell-theory identities written in these variables (Chapter 3);
3. an "interior dynamics", a six-dimensional oscillator with an η-indefinite kinetic form, and a Lie–Poisson closure of pendulum type
   (Chapter 4);
4. the tower Cl(3n,3n): centralisers, dual pairs, Fock gradings, and many statements about two and more "constituents" (Chapters 5, 9);
5. at n = 3, a rank-64 projector P_* built from three commuting involutions whose block carries one Pati–Salam family, the
   272-dimensional space of its Lorentz-invariant mass kernels, and a conditional vacuum analysis (Chapters 6–7);
6. a projector-valued field on a Grassmannian and the identification of an unbroken U(1) at aligned points (Chapter 8);
7. curvature identities, a Kaluza–Klein ansatz with a wrong-sign vector, and obstructions to selecting the Lorentzian four-plane by an
   adjoint scalar (Chapters 11–12);
8. a synthesis, a comparison with the literature, a disclosure of the AI-assisted workflow and a verification package of 66 scripts
   (Chapters 13–16, Appendix D).

The manuscript states plainly that it is not a unified theory, derives no action, makes no prediction, and that every selection beyond
the algebra is an assumption.

## 2. Overall assessment

**Recommendation: not suitable for journal submission in its present form; major restructuring recommended.** The reason is not
mathematical error. In everything I recomputed or checked by hand the stated results were correct, including the headline counts
(§4). The reasons are (i) the absence of an identifiable central new result at the scale of a 428-page work, (ii) a proof standard that
replaces written proofs by scripts for a large fraction of the formal statements, (iii) reliability claims about the verification
package that the package does not support as worded, and (iv) length, redundancy and residue of the drafting process.

**Strengths.**

- *Candour.* The text separates what is proved, constructed, assumed and open more carefully than most work in this area. Chapter 16
  and §14.6 disclaim novelty and prediction explicitly, and the "mechanisms that fail" are reported as results.
- *Correctness at the level tested.* Every finite-dimensional claim I rebuilt independently agreed with the text (§4), and about
  twenty hand checks of formulas across Chapters 2, 5, 6, 7, 10 and 12 found no error of substance.
- *Reproducibility infrastructure.* The realisation is explicit, most scripts rebuild it from scratch instead of sharing code (only 9
  of 66 import `common/clifford8.py`), the run records are archived, and the errata of version 1.0.0 are public and itemised.
- *Disclosure.* The description of the AI-assisted workflow and of its limits in Chapter 15 is exemplary and should be kept.

**Weaknesses, in order of importance,** are developed in §3.

## 3. Major concerns

### M1. What is the contribution?

The manuscript's own deflationary principle (§I.1.2, Corollary 2.20) says that a statement carries information only if it could have
come out otherwise. Applied to the manuscript as a whole, it leaves little:

- §14.6 classifies as *inherited* the Clifford classification, the real form of Dirac theory, the Maxwell identities, the dual pair
  and its Casimir relation, the Pati–Salam content and its anomaly freedom, the count 272 = 2 dim_ℂ Sym²(ℂ¹⁶), Kaluza–Klein reduction
  and the Holst action; and it says of the remainder that "*constructed* means exhibited in this realisation, not new", and that any
  novelty "rests on a targeted search, not on a survey of the literature".
- The results listed as *fixed* are either elementary once stated in standard language or equivalent to an input. Theorem 5.1 is the
  fact that the centraliser of so(E) in so(V) is so(E^⊥). Theorem 5.10 and Corollary 5.12 are the Clebsch–Gordan decomposition of the
  cube of a Dirac spinor, organised by an O(3) × O(1,3) dual pair. Theorem 3.3 is the classical identification of Λ²ℝ^{1,3} with ℂ³
  under SO(3,ℂ), with the Riemann–Silberstein vector. Theorem 4.8 is the remark that an η-self-adjoint complex structure is an
  anti-isometry, hence requires neutral signature.
- Theorem 2.17 ("selection of the imaginary unit") is correct, but the selection is equivalent to the choice of the four-plane. At
  n = 1 the internal bivectors consist of J alone, so the hypothesis "commutes with every internal bivector" is "commutes with J", and
  the theorem says that among the complex structures in the commutant of the gamma matrices (a two-sphere in ℍ, Remark 2.18) the one
  commuting with J is ±J. In Clifford terms: Cl(3,3) = Cl(E) ⊗̂ Cl(U), and the unit bivector of the positive two-plane U commutes with
  Cl(E). The clock and the four-plane are assumptions (Assumption 2.16), so nothing is selected that was not chosen. The sentence
  "the imaginary unit is selected rather than adjoined" (Introduction, Chapters 13 and 16) should be weakened accordingly.
- Chapters 6–8 are, once P_* is fixed, one-family Pati–Salam model building with a scalar in the conjugate of Sym²(16) = 10 ⊕ 126 and
  with adopted potentials. Nothing in §§7.5–7.7 depends on Cl(3,3) beyond the choice of the block; Remark 6.16 ("exhibiting a known
  representation content is a consistency result, not evidence") already says so. The projector itself is chosen by hand, discards
  seven eighths of the module, and puts chirality in by construction (Proposition 7.37).
- The obstruction results are the most promising part (Theorems 12.8 and 12.9, Theorem 6.17 with Example 6.18, the exclusion of the
  two readings in Theorem 3.15), but each holds inside a narrow ansatz that the text itself delimits.

The title promises a "Spacetime–Electromagnetic Framework". The manuscript proves that the two readings of ℝ^{3,3} the title joins
admit no nonzero equivariant intertwiner and cannot share a Lorentz action at one constituent (Proposition 3.6, Theorem 3.15), and lists
their joint realisation as open problem 4 of §16.3. The subtitle describes the content; the title names something the content
obstructs.

**Request.** State, in one page at the front, which results the author claims as new, each with the closest prior result and the
difference. If the honest list is short, the work should be reorganised around it (see §7).

### M2. Proof standard

`verification/THEOREM_CODE.tsv`, the repository's own table, gives for the 464 labelled statements: 194 with no written proof; by
evidence class, 302 "exact" computation, 84 "proof in text", 33 "exact+numerical", 27 "numerical", 14 "argument in text", 3 "inherited",
1 conjecture. `verification/COVERAGE.tsv` shows 102 of the 437 theorems, propositions, lemmas and corollaries with neither a written
proof nor a check cited in the text. (Chapter 15 says 438; the table has 437.)

- Statements called *Theorem* whose only support is a script are computational certificates, not theorems in the usual sense. This
  includes central ones: Theorem 6.10 (P_*), Theorem 3.3 (the field-strength reading), Theorems 2.27, 2.28, 2.30, 2.31, 3.84, 4.28, 5.4.
  Most have short proofs. For Theorem 6.10 a blade-parity argument of a few lines suffices; for the dimension 272 in Theorem 7.17 the
  manuscript itself gives the reason in §14.6 (Sym² of sixteen Weyl fermions) but not at the theorem.
- Twenty-seven statements rest on floating-point computation only, among them Theorem 7.17 (per the table; Chapter 15 nevertheless
  lists 272 among "exact ranks"), Theorems 7.25, 7.26, 7.34 and 8.48, and Lemma 7.24. Theorem 7.25 has a two-line proof that defers its
  two substantive steps to Lemma 7.24 (cited, checked numerically) and to check v1-4 (numerical).
- Proposition 7.33 ("Random-start minimisation") is a report of numerical experiments whose outcome, by its own statement, depends on
  the floating-point path. It is not a proposition. Proposition 9.16 ("Inherited facts") is a list of textbook formulas.
- Many universal statements are tested on samples ("20 random pairs", "at random points of F₃", "for six rational clocks"). An exact
  computation on a sample establishes the sample.

**Request.** Use *Theorem/Proposition* only where a complete proof is printed or cited; introduce a separate environment
(*Computation*, *Numerical observation*) for script-supported statements, with the check identifier in the environment header; write
out the short proofs of the central structural results; and move check identifiers out of running prose (M5).

### M3. The verification package: what its numbers mean

The package is a real asset, but three of the claims made about it are not supported as worded.

**(a) The split "1234 exact, 171 numerical" is produced by a keyword match on the printed claim text, not by the arithmetic used.**
`verification/tools/make_manifest.py` (line 61) labels a check *numerical* if its claim text matches
`float|tol\b|tolerance|numerical|Monte Carlo|quadrature|integration`, and *exact* otherwise. Consequences:

- `ch08_gauge_photon/p1_families.py` has all 47 checks labelled exact. The script imports only NumPy and SciPy, sets `TOL = 1e-9`, and
  its docstring says floating point is used "at random points of the families (orthonormal frames, matrix exponentials), for ranks of
  real spans and for null spaces". Checks p1-34 and p1-40 say "at random points" in their own text.
- d6-2c is labelled exact and prints "R_eff=1.1e-15, |h-h0| <= 1.4e-14"; d4-25 is labelled exact and prints "(numeric)".
- A keyword scan of the archived run records finds at least 140 of the 1234 checks labelled exact whose own claim text describes a
  random sample, a floating-point residual or a numerical procedure. This is a lower bound from wording alone.

Appendix D defines *numerically supported* as "floating-point arithmetic, a random sample or a numerical integration". By that
definition the figure "about one check in eight is numerical" (§D.5), repeated in the README, the release notes and Chapter 15, is too
low, and the per-statement evidence classes in `THEOREM_CODE.tsv` inherit the same error. The classification should be declared in the
code of each check, at the point where the arithmetic is chosen, and sample-based exact checks should be a third class.

**(b) A passing suite has not, historically, meant correct statements.** `docs/errata.md` lists 191 issues in version 1.0.0, 39 of
them marked *Error* (a statement false as written), found by reading while all 1389 checks of that version passed; items 168–191 record
further counterexamples found in the prepared version 1.0.1. Chapter 15 says this candidly. The implication should be drawn in the
abstract-level claims: the residual error rate of statements that go beyond what a check literally asserts is unknown and, on this
record, not small. The keyword "computer-assisted algebra" is appropriate; any wording that suggests the 428 pages are verified is
not.

**(c) Reproducibility is good but not unconditional.** The archived run of version 1.0.1 records that
`x1_unitary_energy_referee.py` hung after 19 of 28 checks, was stopped after 48 minutes and rerun by hand
(`reproducibility/selected-results/v1.0.1/README.md`). The cause is not identified. The runner has no per-script timeout, so an
unattended run can stall indefinitely.

Smaller points: one PASS line often bundles a dozen assertions, so a failure is hard to localise and the count 1405 measures neither
claims nor assertions; script docstrings and claim texts still refer to superseded names ("Phase 9D", "verify/r4_bridge.py",
"RE Thm 8.1", "CORE-06"); Appendix D, first paragraph, still says the scripts are those "at the commit of version 1.0.0".

### M4. Physical content

The manuscript is explicit that it makes no physical claim. A referee must still ask whether the physical readings justify the
framing. As presented they rest on a stack of independent assumptions, each acknowledged individually:

- three time directions, with the restriction to one clock imposed (Assumption 2.16) and the unrestricted six-dimensional field
  ill-posed and without ground state (Proposition 12.1);
- two readings of ℝ^{3,3} that exclude one another;
- an interior dynamics whose energy has signature (6,6) and no ground state (Corollary 4.15), and whose scale ℓ nothing fixes
  (Proposition 4.49); the relation ℏκ = mc² is added by hand (Proposition 4.52(iii)), and Conjecture 9.17 is speculation;
- a matter block chosen by hand, one family, with hypercharge "adopted, not derived" (§6.3.2);
- two different scalar sectors, the kernel field of Chapter 7 and the projector field of Chapter 8, each with an adopted potential
  (the general Pati–Salam-invariant potential on the kernel space has 12 quadratic and 402 quartic fermion-number-conserving
  couplings, Theorem 7.47), which are not combined; the selected vacuum of Chapter 8 is electroweak symmetric; the sign hypothesis (H3)
  of Theorem 8.48 is assumed;
- no gravitational dynamics, and a normal-plane vector of the wrong sign (Theorem 12.15).

Given this, Chapter 10 (astronomy and cosmology of "field rest energy"), §4.7 (planetary orbits and spin–orbit coupling as an image of
the interior dynamics) and the pair-creation readings of Chapter 5 do not earn their place: they attach observational language to
objects whose physical status the manuscript itself leaves open. The conclusion of Chapter 10 is mostly negative ("field rest energy
is not dark matter"), which is correct and could be one remark.

The one-loop analysis of §7.5 (thresholds x < 30, x < 27, x_* ≈ 24.884) is a computation in a specific model: O(272)-symmetric tree
potential, one gauge coupling, one Yukawa coupling, leading logarithm, fixed radius. The text says so (errata item 62), but the
numbers are then quoted in Chapters 13 and 16 with less qualification than they carry in Chapter 7.

### M5. Length, organisation and residue of the drafting process

- **Redundancy.** The chapters are visibly assembled from separate documents (label prefixes `r1:`–`r4:`, `d5:`–`d7:`, `e1:`–`e7:`,
  `p1:`–`p3:`; several sections with their own "Introduction", "Answer" and "Assessment"). Maxwell's equations are derived three times
  (Proposition 3.10, Theorem 3.53, Theorem 3.74); the Rainich identity is stated at least three times (Propositions 3.20(iii), 3.35(i),
  3.94) and again as Proposition 10.1(i); E₆² = E₀² + |S|² at least four times; the statement that surplus-time momentum lowers the
  mass three times (Propositions 3.27, 3.31, 3.101). Chapter 3 has 87 formal statements over 62 pages, most of them inherited.
- **Process residue in the printed text.** The release notes say internal review identifiers were removed. The PDF still contains
  about 350 parenthetical check identifiers in running text and inside theorem statements ("(p1-21)", "(r3-1, r3-2, r3-3)"), about ten
  paragraphs that open with the bare status word "Fixed.", status tags inside theorem titles ("…; fixed", "…; restated", "…; fixed, and restated"), lines
  "Executable check: a4-4 and a10-2", and remarks about "an earlier and superseded version of the programme" (Remark 2.15, §3.5).
  Chapter 15 refers to "the review checklist that governs the present revision" and "the CORE review items", which the reader has no
  access to.
- **Notation.** Table I.1 promises one meaning per symbol. T is the temporal trivector, but the stress tensor appears as T, T_M, T₀,
  T^{(γ)}, T(z), T_tot beside Θ (see §6); ω is both the volume element and, indexed, the pole operators; E denotes the four-plane, an
  energy, field components, matrix units E_ij, generators E_A and E₆; the sign triple (−,−,−) has three meanings (Remark 2.15).
- **Classification overload.** Four classes, three provenance labels and five link statuses are applied sentence by sentence. The
  effect is that the reader cannot tell a main result from a remark. Sorting belongs in one table per chapter.
- **Decorative structure.** The Fisher–Rao clock simplex and the entropy maximum (§2.1.3, Proposition 2.5) are followed by the
  statement that every unit clock is equivalent and none is selected (Remark 2.6). Theorem 2.2 is Schur's lemma. These can go.

My estimate is that the mathematics a reader needs occupies well under a quarter of the 428 pages.

### M6. Related work

Chapter 14 is careful about Krasnov, two-time physics and the spectral Pati–Salam model. Two further works found in a brief search
appear, from their abstracts, close to the method of Chapters 5–6 and are not cited; I have not read them in full:

- arXiv:2108.08274, a unification scheme in Cl(7,7) in which seven commuting elements with binary eigenvalues label the fermion states,
  which is the same device as building P_* and the charges from commuting blade involutions (and as Theorem 6.17);
- arXiv:2308.12293, which indexes Spin(N) spinors by bit codes and studies Spin(10) ⊂ Spin(11,1) with the Lorentz group, comparable
  to the Fock/octant description of §1.2 and to the embedding question of §14.3.

For §4.7 the map z ↦ z² between Hooke and Kepler orbits is cited (Bohlin, Levi-Civita); for Chapter 3 the Pauli-algebra form of
Maxwell's equation in Theorem 3.53(iv) deserves a reference to the algebra-of-physical-space literature alongside Hestenes.

---

## 4. Independent verification carried out for this review

**Recomputed in separate code** (NumPy/SciPy, matrices built from equations (1.1)–(1.2); nothing imported from the repository). All
agree with the manuscript.

| statement | result recomputed |
|---|---|
| Remark 2.18 | commutant of the four gamma matrices in M₈(ℝ) has dimension 4 |
| Theorem 2.17 at n = 1 | joint commutant of the interior Lorentz algebra and J has dimension 4 |
| Proposition 2.29 | solutions of XᵀK^μ + K^μX = 0 have dimension 4; with the mass form, dimension 1 |
| Theorem 2.25 | antilinear intertwiners of the Dirac operator form a space of dimension 2 |
| Example 1.6 | the paired-mode projector has rank 4 and a 9-dimensional stabiliser in so(3,3) |
| Theorem 3.1 | the centraliser of G in so(3,3) has dimension 6 |
| Theorem 6.10 | P_* is an orthogonal projector of rank 64, commutes with the Lorentz generators and J, and P_*γ^μP_* = 0 |
| Theorem 6.11 | the internal bivectors commuting with P_* form a 23-dimensional space (full linear solve over the 91), all compact, acting faithfully on the block |
| Theorem 7.17(i) | the Lorentz-invariant antisymmetric kernels on the block form a space of dimension 272 |
| Theorem 5.10 | diagonal Casimir spectrum {3/2: 160, 11/2: 288, 15/2: 64} on ℝ⁵¹² |
| Theorem 5.9 | the centraliser of the diagonal Lorentz algebra in so(9,9) has dimension 18 |

**Checked by hand.** No error found.

- Theorem 2.17: the parity argument, and (ω_I)² = −1 at every rung.
- Proposition 2.5: the Hessian eigenvalues −9, −3.
- Theorem 3.53: the expansion of (∂_t − D)Φ_oct term by term; Proposition 3.34(ii): the dimensions 0, 2, 3, 1 by representation theory.
- Theorem 5.10 and Corollary 5.12: the multiplicities 64, 288, 160 and the Casimir values from the decomposition of D ⊗ D ⊗ D.
- Theorem 5.31(v): boundedness at κ₁ = κ₂ = |β|, where a zero frequency occurs but the kernel of S† lies in the image of S.
- Proposition 5.17(ii): the charge tables at n = 3, 5, 7.
- Proposition 6.14: all five anomaly coefficients. Theorems 6.20 and 6.21: the factorisations 6N_c q(h_ℓ² − h_q²) and
  6N_c q(N_c²q² − h_q²). Proposition 6.22 and Example 6.18: all numbers.
- §7.5: equality in |𝒟| ≤ (3/8)I₄ at the identity; the recoupling identity evaluated at the identity; the benchmark difference
  0.9 − |κ|/8; the necessary and sufficient conditions of Theorem 7.26(a)–(c); the trace-potential dictionary; the sum rules and the
  mode count 60 in Theorem 7.27; the thresholds x < 30 and 48 − (7/4)x in Theorem 7.30; the vector-mass ratio 5 in Theorem 7.34(b),
  which agrees with the left–right literature for g₄ = g_R.
- Proposition 10.3: the density 2t from the Wishart eigenvalue law. Proposition 10.6: the numerical example.
- Theorem 11.13: the value ½ at the linearised level. Theorem 12.8(i): t₃ and t₄ are excluded by boundedness below; (iii):
  t₄ = 4804/625.

**Bibliography.** Three less familiar entries were checked against arXiv and exist as cited (0901.4230, 2509.13654 including the value
0.215° ± 0.074°, 2602.19151). No sign of fabricated references was found in the sample.

## 5. Re-run of the verification suite

Run on commit `75af386` with `python3 verification/run_all.py --out <scratch> -j 3`, Python 3.13.18, NumPy 2.5.3, SciPy 1.18.1,
SymPy 1.14.0, mpmath 1.3.0 (the versions of `requirements.txt`), Linux x86_64, one BLAS thread per script.

- **Result: 66 scripts, 1405 PASS, 0 FAIL**; every exit status is zero and consistent with the PASS/FAIL lines.
- Compared line by line with the archived records in `reproducibility/selected-results/v1.0.1/`: the same 1405 check identifiers, the
  same status for each, and identical claim text for all but r4-15, which prints an elapsed time. This includes the printed
  floating-point residuals and the end points of the random-start minimisations (23 of 100 runs below the Standard-Model energy at
  x = 25, x_* = 24.8840501).
- No script hung in this run. The longest was `a10_chiral_certificates.py` at 28 minutes; the summed script time was 87 minutes.
- One script (`p3_identification.py`) writes to standard error, as the archived README notes.

So the release claim "66 scripts, 1405 checks, all PASS" is reproduced exactly. What the 1405 PASS lines establish is the subject of
M2 and M3.

## 6. Specific corrections

Line numbers refer to the files in `manuscript/chapters/`.

1. **Corollaries 5.2 and 5.6** (`05_clifford_tower.tex`, lines 39–60). The title "The photon at one constituent", the sentence "The
   electromagnetic u(1) is not a symmetry of the interior structure but the symmetry", and "From n = 2 onwards the photon is one
   direction in a space that grows without bound; its uniqueness is lost…" identify the bivector e₂e₃ with the photon. Errata items 67
   and 156 withdraw exactly this identification elsewhere, the caption of Figure 1.1 now denies it, and Chapter 8 shows that this
   charge is half the fermion number and anomalous on the block. The two corollaries were missed. The column "u(1) alone" of
   Theorem 5.4 has the same problem in milder form.
2. **Clock symbol without hat** (search-and-replace remnants; the letter n is glued to the next symbol). `02_lorentz_structures.tex`
   lines 346 ("nu₁u₂ = I_T"), 351 and 367 ("Ω_n̂ := ns₁s₂s₃"), 374 ("c ns_ap_a"); `03_readings_of_x.tex` lines 1780, 1781, 1854
   ("R = nu₂"), 1866, 1925; also "pr_n̂ n = n̂" in the proof of Proposition 2.11. About seventeen occurrences.
3. **Stress tensor written T instead of Θ.** `03_readings_of_x.tex` lines 587–588 ("Θ^μ_λ T^λ_ν = E₀²δ", "E₀ = ½√(T_{μν}Θ^{μν})"),
   620–652 (T_M, T₀, T^M_{μν}, T^{(γ)}_{μν}), 730–731 (T_M); `05_clifford_tower.tex` lines 657–658, 2117–2230 (T(z), T_tot). T is the
   temporal trivector everywhere else.
4. **Chapter 15**, "Of the 438 theorems, propositions, lemmas and corollaries": `COVERAGE.tsv` has 437 (143 + 240 + 20 + 34), and
   `THEOREM_CODE.tsv` has 33 corollaries where `COVERAGE.tsv` has 34 (one corollary of Chapter 2, after Theorem 2.25, has no label).
5. **Chapter 15**, class (ii): "Kernel dimensions such as 272 … are exact ranks". The checks that compute 272 (v2-1, z3-8) are
   floating-point, and `THEOREM_CODE.tsv` classes Theorem 7.17 as numerical. Either add the exact blade count that Theorem 7.17(i)
   describes as a check, or give the two-line representation-theoretic proof.
6. **Appendix D**, first paragraph: "at the commit of version 1.0.0" should read 1.0.1.
7. **Chapter 2**, line 362: the sentence "Checks: r2-3a and r2-3c … r2-3d)." is a fragment with an unmatched parenthesis.
8. **Introduction**, summary paragraph on the arena: "the entropy maximum of the clock simplex under the Fisher–Rao metric". The
   entropy maximum does not involve a metric; the two properties are separate.
9. **Chapter 15**: "developed over roughly three years" and "the repository history of the programme begins in June 2026" should be
   reconciled in one sentence, since the second invites the question of where the earlier record is.
10. **Proposition 7.33** should not be a proposition (M2). **Proposition 9.16** should be a paragraph of citations.
11. **Theorem titles** carrying status tags ("; fixed", "; restated", "; inherited, with proof") in Chapter 8: move the tag out of the
    title.

## 7. A possible path to publication

The material would be better served by a small number of focused papers, each with complete written proofs, a stated novelty claim
checked against the literature, and the scripts as supplementary material:

1. **Obstructions to selecting a Lorentzian four-plane in ℝ^{3,3}** (Theorems 12.4, 12.8, 12.9, Propositions 12.10–12.11, the
   exclusion of the two readings). Short, self-contained, and the part most likely to be new and of interest to readers of the
   multi-time literature.
2. **The (10,1,3) scalar potential of Pati–Salam models** (Theorems 7.23–7.30): the complete renormalisable potential with the
   holomorphic ΔB = 2 quartic, the necessary and the sufficient conditions for the Standard-Model orbit, the second variation and its
   sum rules. This is ordinary model building, independent of Cl(3,3), and should be written and compared as such with the left–right
   and Pati–Salam literature on scalar potentials. Its value depends on a priority check I have not made.
3. **A chiral Pati–Salam block in Cl(9,9) from commuting involutions** (Chapter 6 with Theorem 7.17), as a short note that includes
   the explicit comparison with Krasnov's module left open in Remark 6.24.
4. **An expository account** of Chapters 1–3, if wanted, clearly labelled as such.
5. **A methods paper** on the AI-assisted workflow and its measured error rates: 191 errata, 39 false statements under a fully passing
   suite, what kind of review found them. This is unusual data and of independent interest.

The interior dynamics (Chapter 4), the multi-constituent chapters (most of 5 and 9), Chapter 10 and §4.7 are, in my assessment, not at
present publishable results, and removing them would strengthen what remains.

## 8. Questions for the author

1. Which three results do you regard as the main new theorems, and what is the nearest published statement to each?
2. Is there any statement in Chapters 6–8 whose truth depends on the arena being Cl(9,9) and not on the Pati–Salam content alone?
3. Proposition 7.37 shows chirality is imposed by the choice of P_*. Is there a candidate principle, even conjectural, that
   distinguishes P_* from γ⁰P_*γ⁰, or is this to be recorded as a ninth row of Table 13.1?
4. For the 194 statements without a written proof: which would you be unable to prove by hand in under a page? Those are the ones a
   reader most needs flagged.
5. What caused the hang of `x1_unitary_energy_referee.py` in the archived run?
