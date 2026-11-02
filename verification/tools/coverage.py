#!/usr/bin/env python3
"""Map the claim statements of the global manuscript to the evidence cited for them (Phase 9A, 2 November 2026).

    python3 verification/coverage.py

For every theorem, proposition, lemma, corollary, conjecture and example of manuscript/chapters/ (chapters and appendices) it takes
the span from the start of the environment to the start of the next such environment or the next sectioning command, and records:
- the check identifiers of verification/MANIFEST.tsv cited in that span, split into exact and numerical;
- whether a written proof (a proof environment) occurs in the span;
- whether the statement itself cites the literature or carries the words inherited, restated, interpretive or conjecture;
- its provenance: the source unit of the nearest preceding '% --- source:' marker left by the migration (tex/ is the 26 October
  manuscript; the others are the investigation notes), or 'hand' when the file has no such marker before it.
It writes verification/COVERAGE.tsv and prints a summary per chapter, together with the check identifiers cited in the manuscript
that are not in the manifest.

What this measures is citation, not proof. A statement with a cited check is one for which some coded assertion exists; whether the
assertion is the statement is a matter for the audit (Phase 9D). A statement with no cited check may still be proved in the text, be
inherited from the literature, or be checked by a script that the text does not cite by identifier. Identifier ranges written in
the text (such as d1-4a--c) count only through the identifiers spelled out in full.
"""
import pathlib
import re

ROOT = pathlib.Path(__file__).resolve().parents[2]
ENVS = ("theorem", "proposition", "lemma", "corollary", "conjecture", "example")
BEGIN = re.compile(r"\\begin\{(" + "|".join(ENVS) + r")\}")
SECT = re.compile(r"\\(?:chapter|section|subsection)\*?\{")
GENERIC_ID = re.compile(r"(?<![\w-])((?:rb|[a-z])\d+-\d+[a-z]?)(?![\w-])")


def manifest():
    rows = (ROOT / "verification" / "MANIFEST.tsv").read_text().splitlines()[1:]
    return {r.split("\t")[0]: r.split("\t")[4] for r in rows if r.strip()}


def files():
    ms = ROOT / "manuscript"
    return sorted(f for f in (ms / "chapters").glob("*.tex") if f.stem not in ("appendix_d_table", "appendix_d_scripts"))


def main():
    method = manifest()
    ids = sorted(method, key=len, reverse=True)
    idre = re.compile(r"(?<![\w-])(" + "|".join(map(re.escape, ids)) + r")(?![\w-])")
    rows, unknown, cited_all = [], {}, set()
    for f in files():
        text = f.read_text()
        cuts = sorted([m.start() for m in BEGIN.finditer(text)] + [m.start() for m in SECT.finditer(text)] + [len(text)])
        for m in BEGIN.finditer(text):
            env, start = m.group(1), m.start()
            end = next(c for c in cuts if c > start)
            span = text[start:end]
            stmt = span[: span.find(f"\\end{{{env}}}")] if f"\\end{{{env}}}" in span else span
            lab = re.search(r"\\label\{([^}]*)\}", stmt)
            label = lab.group(1) if lab else f"(unlabelled, line {text.count(chr(10), 0, start) + 1})"
            found = sorted(set(idre.findall(span)))
            cited_all.update(found)
            exact = [c for c in found if method[c] == "exact"]
            numer = [c for c in found if method[c] != "exact"]
            proof = "\\begin{proof}" in span
            lit = "\\cite" in stmt
            src = [s for s in re.finditer(r"% --- source: ([^#\s]+)", text[:start])]
            prov = src[-1].group(1) if src else "hand"
            marks = ",".join(w for w in ("inherited", "restated", "interpretive", "conjecture") if w in stmt.lower())
            rows.append((f.stem, env, label, prov, len(exact), len(numer), "yes" if proof else "no", "yes" if lit else "no",
                         marks, " ".join(found)))
        for g in GENERIC_ID.findall(text):
            if g not in method:
                unknown.setdefault(g, set()).add(f.stem)
    head = "file\tenvironment\tlabel\tprovenance\texact_checks\tnumerical_checks\twritten_proof\tcites_literature\tmarkers\tcheck_ids\n"
    (ROOT / "verification" / "COVERAGE.tsv").write_text(head + "".join("\t".join(map(str, r)) + "\n" for r in rows))
    print("file\tclaims\twith_check\tnumerical_only\tproof_no_check\tneither")
    tot = [0, 0, 0, 0, 0]
    for stem in dict.fromkeys(r[0] for r in rows):
        rs = [r for r in rows if r[0] == stem and r[1] != "example"]
        c = [len(rs), sum(r[4] + r[5] > 0 for r in rs), sum(r[4] == 0 and r[5] > 0 for r in rs),
             sum(r[4] + r[5] == 0 and r[6] == "yes" for r in rs), sum(r[4] + r[5] == 0 and r[6] == "no" for r in rs)]
        tot = [a + b for a, b in zip(tot, c)]
        print(stem + "\t" + "\t".join(map(str, c)))
    print("TOTAL\t" + "\t".join(map(str, tot)))
    print("provenance\tclaims\twith_check\tproof_no_check\tneither")
    for src in ("tex", "rest_energy", "two_readings", "emergence", "photon_selection", "hand"):
        rs = [r for r in rows if r[1] != "example" and r[3].split("/")[0] == src]
        print(f"{src}\t{len(rs)}\t{sum(r[4] + r[5] > 0 for r in rs)}\t{sum(r[4] + r[5] == 0 and r[6] == 'yes' for r in rs)}"
              f"\t{sum(r[4] + r[5] == 0 and r[6] == 'no' for r in rs)}")
    print(f"manifest checks cited inside claim spans: {len(cited_all)} of {len(method)}")
    print("identifiers cited in the manuscript but absent from the manifest:",
          "; ".join(f"{k} ({','.join(sorted(v))})" for k, v in sorted(unknown.items())) or "none")


if __name__ == "__main__":
    main()
