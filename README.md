# The Spacetime–Electromagnetic Framework

**Mathematical Structures and Obstructions in Split-Signature Clifford Geometry** — Driss Bencheikh, Independent Researcher.

**Version 1.0.0 — Preprint, not peer reviewed** (2 November 2026).

**DOI:** [10.5281/zenodo.23900465](https://doi.org/10.5281/zenodo.23900465) (Zenodo record).

This repository contains one research manuscript, its LaTeX source, and the computational verification package that checks its
statements.

- **PDF:** [`release/manuscript-v1.0.0.pdf`](release/manuscript-v1.0.0.pdf)
- **Source:** [`manuscript/main.tex`](manuscript/main.tex) (16 chapters and Appendices A–E; bibliography `references.bib`)
- **Verification:** [`verification/`](verification/) (61 scripts, 1389 checks, one directory per chapter)
- **Results:** [`reproducibility/manifest.json`](reproducibility/manifest.json) and [`reproducibility/selected-results/`](reproducibility/selected-results/)
- **Documentation:** [`docs/mathematical-assumptions.md`](docs/mathematical-assumptions.md),
  [`docs/known-obstructions.md`](docs/known-obstructions.md), [`docs/research-status.md`](docs/research-status.md)
- **Errata:** [`docs/errata.md`](docs/errata.md) — known issues in version 1.0.0, to be corrected in version 1.0.1

## What the manuscript investigates

The real Clifford algebra Cl(3,3) ≅ M₈(R) of the split space R^{3,3} — three temporal directions paired with three spatial ones — and its
graded tower Cl(3n,3n): which ingredients of relativistic field theory the algebra *fixes*, which it merely *admits*, and which it
*obstructs*.

**It is not a completed unified physical theory.** Every selection beyond the algebra (clock, chiral projector, vacuum potential,
couplings, gauge principle) is an explicit assumption; no gravitational dynamics, matter coupling or low-energy limit is derived; global
minimality of the vacuum is a conjecture; the photon is identified only under an adopted potential. No experimental confirmation is
claimed. See [`docs/research-status.md`](docs/research-status.md).

## Repository layout

| path | content |
|---|---|
| `manuscript/` | `main.tex`, `preamble.tex`, `chapters/` (chapters and appendices), `figures/`, `references.bib`, `esframework.bst` |
| `release/` | the PDF of version 1.0.0 and its checksum |
| `verification/` | `run_all.py`, the scripts by chapter (`chNN_<topic>/`), shared code (`common/`), generators (`tools/`), `MANIFEST.tsv`, `THEOREM_CODE.tsv`, `COVERAGE.tsv`, `requirements.txt` |
| `reproducibility/` | `manifest.json` (one record per check) and `selected-results/` (the run records of version 1.0.0 and the indexes of exact and numerical checks) |
| `docs/` | assumptions, obstructions, research status |

## Build the PDF

TeX Live with the packages listed in [`reproducibility/selected-results/latex-environment.md`](reproducibility/selected-results/latex-environment.md):

```bash
cd manuscript && pdflatex main && bibtex main && pdflatex main && pdflatex main
```

## Run the verification

Python 3 with numpy, scipy, sympy and mpmath:

```bash
pip install -r verification/requirements.txt
```

```bash
python3 verification/run_all.py --out runs/local -j 3
```

The whole suite takes about 20 minutes on three cores; a single script runs on its own, for example
`python3 verification/ch03_readings_of_x/r1_lagrangians.py`. Each check prints `PASS` or `FAIL`, and each script exits with status zero
only if all its checks pass. `python3 verification/tools/make_manifest.py --runs runs/local` rebuilds the manifest from a run.

## AI assistance

The work was carried out with substantial assistance from AI systems under the author's direction; Chapter 15 of the manuscript describes
the workflow, the author's role and the limits of AI review. All results are the author's responsibility.

## How to cite

Driss Bencheikh, *The Spacetime–Electromagnetic Framework: Mathematical Structures and Obstructions in Split-Signature Clifford
Geometry*, version 1.0.0, preprint (not peer reviewed), 2026, doi:10.5281/zenodo.23900465.
Zenodo record: https://doi.org/10.5281/zenodo.23900465. Source repository: https://github.com/driss-bencheikh/spacetime-electromagnetic-framework.
Machine-readable metadata: [`CITATION.cff`](CITATION.cff).

## Licence

- Manuscript, LaTeX source, PDF, documentation and result records: **CC BY 4.0** ([`LICENSE-CC-BY-4.0.txt`](LICENSE-CC-BY-4.0.txt)).
- Verification software (every Python file): **MIT** ([`LICENSE-MIT.txt`](LICENSE-MIT.txt)).

[`LICENSE.md`](LICENSE.md) states which files each licence covers.
