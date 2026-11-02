# LaTeX environment for the global manuscript

Set up in the cloud container on 1 November 2026 (Ubuntu 24.04, TeX Live 2023/Debian, pdfTeX 3.141592653-2.6-1.40.25):

```bash
apt-get install -y --no-install-recommends texlive-latex-base texlive-latex-recommended texlive-latex-extra \
  texlive-pictures texlive-science texlive-fonts-recommended texlive-plain-generic lmodern
```

`lmodern` is required (the first build failed with "File `lmodern.sty' not found" without it). `poppler-utils` (pdftoppm) was added
for visual inspection of rendered pages.

Build, from `manuscript/`:

```bash
pdflatex main && bibtex main && pdflatex main && pdflatex main
```

The bibliography is `references.bib`, formatted by `esframework.bst` (in the same directory).
