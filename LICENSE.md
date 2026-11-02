# Licences

Copyright © 2026 Driss Bencheikh. This repository is released under two licences, chosen by the kind of material.

| material | licence | files |
|---|---|---|
| **Manuscript, text, results and documentation** | [Creative Commons Attribution 4.0 International (CC BY 4.0)](LICENSE-CC-BY-4.0.txt), SPDX `CC-BY-4.0` | `manuscript/` (LaTeX source, bibliography, bibliography style, figure data), `release/` (the PDF and its README), `reproducibility/` (manifest and run records), `docs/`, and every `*.md`, `*.tsv`, `*.json`, `*.txt` and `*.cff` file not listed below |
| **Verification software** | [MIT](LICENSE-MIT.txt), SPDX `MIT` | every `*.py` file (`verification/`, including `verification/common/` and `verification/tools/`) |

Notes.

- The data files read by two provenance checks under `verification/common/provenance/` (label lists and one LaTeX excerpt of an earlier
  version of the manuscript by the same author) are text, under CC BY 4.0.
- No third-party code or text is bundled in this repository; the cited literature is referenced, not reproduced. The Python packages the
  scripts import (NumPy, SciPy, SymPy, mpmath) and the TeX Live packages used to build the PDF are installed separately under their own
  licences.
- Attribution for CC BY 4.0: *Driss Bencheikh, “The Spacetime–Electromagnetic Framework: Mathematical Structures and Obstructions in
  Split-Signature Clifford Geometry”, version 1.0.0 (preprint, not peer reviewed), 2026,
  https://github.com/driss-bencheikh/spacetime-electromagnetic-framework*, with an indication of any changes.
