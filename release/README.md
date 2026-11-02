# Version 1.0.0 — Preprint, not peer reviewed

`manuscript-v1.0.0.pdf` — *The Spacetime–Electromagnetic Framework: Mathematical Structures and Obstructions in Split-Signature
Clifford Geometry*, Driss Bencheikh (Independent Researcher). 412 pages: title, abstract and keywords, contents, introduction with
notation and dependency map, Chapters 1–16, Appendices A–E (Appendix D: verification and reproducibility), bibliography (281 entries).
Prepared 2 November 2026.

- **Status:** preprint, not peer reviewed; no DOI.
- **Source:** `manuscript/main.tex` in the same commit as this file. Build from `manuscript/` with
  `pdflatex main && bibtex main && pdflatex main && pdflatex main` (TeX Live 2023, pdfTeX 1.40.25).
- **Build checks:** no undefined references or citations; no overfull boxes. Embedded PDF metadata: title, author, subject
  ("Version 1.0.0. Preprint, not peer reviewed"), keywords.
- **Verification:** 61 scripts, 1389 checks (1220 exact, 169 numerical), all PASS in `reproducibility/selected-results/v1.0.0/`.
- **SHA-256:** `4fd6ffa39b38911bfd4c550e630f564aec404a930d5173e2f3e9dbd928af9bcb` (3 243 693 bytes). A rebuild changes the embedded
  creation date and therefore the checksum.
- **Licence:** CC BY 4.0 (`LICENSE-CC-BY-4.0.txt`); the verification software is under MIT. See `LICENSE.md`.
