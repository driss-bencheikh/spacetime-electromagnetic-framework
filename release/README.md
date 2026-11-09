# Version 1.0.1 — Preprint, not peer reviewed

`manuscript-v1.0.1.pdf` — *The Spacetime–Electromagnetic Framework: Mathematical Structures and Obstructions in Split-Signature
Clifford Geometry*, Driss Bencheikh (Independent Researcher). 428 pages: title, abstract and keywords, contents, introduction with
notation and dependency map, Chapters 1–16, Appendices A–E (Appendix D: verification and reproducibility), bibliography (326 entries).
Prepared 8 November 2026. `manuscript-v1.0.0.pdf`, the first version, is kept unchanged.

- **Status:** preprint, not peer reviewed.
- **What changed:** the issues listed in `docs/errata.md` are corrected (all items; items 59 and 60 by new exact checks). In particular:
  - errors, missing hypotheses, imprecise wording and missing attributions are corrected in every chapter;
  - a new abstract and a new opening of the Introduction; the rest-energy hypothesis is presented as historical motivation, with what
    survived the tests stated first;
  - Chapter 12 opens with its motivation and scope, uses the corrected worldline action, and treats the universal evolution parameter as
    prior art (Horwitz and Piron);
  - Chapter 13: one link of Figure 13.1 is marked conditional, and a new full-page Figure 13.2 shows the five research chains in detail;
  - Chapter 14 is rewritten as a comparison with the literature, including two-time physics, SO(3,3) gauge theories, Krasnov's
    construction and the spectral Pati–Salam model;
  - internal review identifiers are removed from the printed text; 41 references are added;
  - the verification package gains exact certificates (a9, a10; a3-9 and a3-10 made exact; c1-12), a numerical check of the corrected
    kernel structure (a11), a derived mass shell in c2-8, a corrected record e7-26 and an exit status for `p2_alignment.py`; the boundary
    cases and counterexamples of the reviews of the prepared version (items 168–191) are exact or numerical checks (a12, a3-11, d7-4c, d6-8,
    a11-3); the isomorphism of Proposition 7.22 and the constant 3/8 of Theorem 7.23(ii) are certified exactly (a13, items 59, 60);
  - the theorem-to-code table and the list of scripts (Section D.7 and Tables D.1–D.2 of version 1.0.0) are no longer printed: they are
    `verification/THEOREM_CODE.tsv` and `verification/SCRIPTS.tsv`; the title page carries the title, author, affiliation and date only.
- **Numbering:** theorem, proposition, lemma, corollary, definition, example, remark, table and figure numbers are those of version 1.0.0;
  Figure 13.2 is new, Section D.7 and Tables D.1–D.2 are removed; Sections 14.2–14.4 of version 1.0.0 are now 14.3, 14.4 and 14.6; equations (7.6)–(7.7) are now (7.7)–(7.8). Page numbers
  differ.
- **DOI:** version 1.0.0: [10.5281/zenodo.23900465](https://doi.org/10.5281/zenodo.23900465) (Zenodo record). Version 1.0.1 receives its
  own DOI when it is deposited.
- **Cite as:** Driss Bencheikh, *The Spacetime–Electromagnetic Framework: Mathematical Structures and Obstructions in
  Split-Signature Clifford Geometry*, version 1.0.1, preprint (not peer reviewed), 2026.
- **Source:** `manuscript/main.tex` in the same commit as this file. Build from `manuscript/` with
  `pdflatex main && bibtex main && pdflatex main && pdflatex main` (TeX Live 2023, pdfTeX 1.40.25).
- **Build checks:** no undefined references or citations; no overfull boxes. Embedded PDF metadata: title, author, subject
  ("Version 1.0.1. Preprint, not peer reviewed"), keywords.
- **Verification:** 66 scripts, 1405 checks (1234 exact, 171 numerical), all PASS in `reproducibility/selected-results/v1.0.1/`.
- **SHA-256:** `c0e86ccf9a4aa180e823d5d17fb0bddf857ec16648e0e5e0e78239a15ca3d7a2` (3 413 923 bytes). A rebuild changes the embedded creation date and therefore the checksum.
- **Licence:** CC BY 4.0 (`LICENSE-CC-BY-4.0.txt`); the verification software is under MIT. See `LICENSE.md`.
