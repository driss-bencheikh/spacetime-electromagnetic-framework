# Vacuum selection checks (Phase 9B, 2 November 2026)

Reproducible checks of Section 7.3 of the manuscript (`manuscript/chapters/07_mass_vacuum.tex`, from "The tensor model" to "The full kernel
space"): the renormalizable potential on (10,1,3), its second variation, the one-loop selection on the tensor model, the 271-dimensional
angular Hessian on the full kernel space, and the random-start minimisations with the crossing x\*. They close the evidence gaps of review
items CORE-09 and CORE-10.

```bash
python3 verification/ch07_mass_vacuum/v1_tensor_model.py
```

```bash
python3 verification/ch07_mass_vacuum/v2_kernel_hessian.py
```

```bash
python3 verification/ch07_mass_vacuum/v3_minimisation.py
```

| script | checks | runtime | what it tests |
|---|---|---|---|
| `v1_tensor_model.py` | v1-1 … v1-14 | ~30 s | Sym²T components, recoupling identity, lem:hw, thm:select0 (orbit, stabiliser), the holomorphic quartic (invariance, degree, the bound \|D\| ≤ (3/8) I₄ by Schur's lemma, vanishing to second order at S₀), the benchmark numbers, thm:secondvar (60×60 mass matrix; **the d_Rd_R entry of the 26 October version is wrong**), eq:sumrules, thm:oneloop1013 (c = (33,3,−3,39,5), G, F, masses, x < 30, 192/7) |
| `v2_kernel_hessian.py` | v2-1 … v2-10 | ~1 min | the 272-dimensional kernel space from its defining conditions, its components and types, Φ₀, E_x(u₀) = 2x − 99, thm:fullsecond (H = H_G + xH_F, H_F = −8(1 − Π_ν), every entry of tab:271, x < 27), the thresholds, agreement with the tensor model, prop:LR and the great circle |
| `v3_minimisation.py` | v3-1 … v3-7 | ~10 min | random-start minimisations (archived seeds, Gaussian starts, L-BFGS-B and conjugate gradient, tolerances, residuals) on T (tree potential with κ; one loop) and on the kernel space (one loop), the competitor at x = 25 and the crossing x\* |

Shared modules: `kernel_space.py` (the block, the kernel space, E_x with analytic gradient and Hessian, checked against finite differences)
and `tensor_model.py` (the tensor model with analytic gradients).

**What the checks establish and what they do not.** All three scripts use floating-point linear algebra on exactly specified matrices;
tolerances are printed. v1 and v2 test identities, spectra and Hessians: they establish the stated local statements within the stated
approximation (leading logarithm, curvature masses, common couplings). v3 reports numerical searches; a random search, however large, does
not establish a global minimum, and the conjectures of the manuscript remain conjectures. The invariant-theory multiplicities
(1,0,0,0,1,0,5) quoted after thm:potential, and the invariant counts on the whole kernel space, are certified exactly by
`verification/ch07_mass_vacuum/z4_invariant_counts.py` (z4-1 to z4-4), not here.

## Release audit (version 1.0.0)

Two scripts recompute the Chapter 7 statements that the Phase 9D audit had left without a recomputation:

```bash
python3 verification/ch07_mass_vacuum/a6_kernels_cp_certificates.py
```

```bash
python3 verification/ch07_mass_vacuum/a8_electroweak_numerical.py
```

| script | checks | runtime | what it tests |
|---|---|---|---|
| `a6_kernels_cp_certificates.py` | a6-1 … a6-5 | ~2 min | **exact**: thm:subkernels (blades of Cl(3,3)); thm:vacCP (Φ₀ rebuilt in integer arithmetic on R^512); the algebraic core of prop:coupled3; thm:Bcounts graded by F (a6-4) and by U(1)_{1−Γ} and the JΓ charge (a6-5, Molien–Weyl) |
| `a8_electroweak_numerical.py` | a8-1 … a8-9 | ~2 min | **numerical** on `kernel_space.py`: prop:tensor, prop:EW (i)–(iv), thm:mA (structural facts), prop:rankfive, prop:bidoublet, and the component charges used by a6-5; every rank is read with an asserted spectral gap (discarded < 1e-9 for direct SVDs, < 1e-6 for commutants from Gram matrices; kept > 1e-3) |

The Chapter 8 counterpart is `verification/ch08_gauge_photon/a7_projector_geometry_certificates.py` (exact: prop:config, thm:holonomy,
prop:chiralpreserved, thm:towervac, thm:S4, thm:towercurved, prop:towermatter, and prop:spectator at n = 3). What remains without a
check is listed in `docs/research-status.md` ("Partial evidence that remains").

## Corrections audit (v1.0.1)

```bash
python3 verification/ch07_mass_vacuum/a11_kernel_structure.py
```

a11-1 and a11-2, numerical (float64 on exactly specified integer matrices, about one minute): the complex structure C~(M) = J M Gamma of
thm:kernelspace and its vanishing cubic trace (the structure of version 1.0.0 gives +-1152), and the rows with v_R = 0 of the strata table of
thm:axion.
a11-3 (numerical): the stratum v_R != 0, v_u = v_d = 0 of the same table (stabilisers 12 and 14).

```bash
python3 verification/ch07_mass_vacuum/a13_tensor_certificates.py
```

a13-1 and a13-2, exact (rational arithmetic, about three minutes): the constant 3/8 of thm:potential (ii) on an integer basis of the (20',1)
component, and the explicit isomorphism iota of prop:tensor with its norm and spectral identities (errata items 59 and 60).
