#!/usr/bin/env python3
"""Generate verification/MANIFEST.tsv, the check-level manifest of the verification package.

    python3 verification/tools/make_manifest.py [--runs reproducibility/selected-results/v1.0.1]

For every PASS/FAIL line in the run record <runs>/<chapter_dir>__<stem>.out of each script verification/<chapter_dir>/<stem>.py it
writes one row: check id, status, script, command, exact or numerical, the manuscript chapters and labels that cite the check, the run
record, and the claim text. A check is attributed to a manuscript label when its identifier occurs in a chapter file and a \\label
precedes it within the same environment or paragraph (the nearest \\label within 4000 characters). The generated tables of Appendix D
are not counted as citations. Reads only; exit status 0.
"""
import pathlib
import re

ROOT = pathlib.Path(__file__).resolve().parents[2]
import sys
RUNS = ROOT / (sys.argv[sys.argv.index("--runs") + 1] if "--runs" in sys.argv else "reproducibility/selected-results/v1.0.1")
LINE = re.compile(r"^(PASS|FAIL) ([A-Za-z0-9]+-[0-9A-Za-z.,+-]*|[A-Za-z0-9]+):\s*(.*)$")


def latest_outputs():
    return {out.stem: out for out in sorted(RUNS.glob("*__*.out"))}


def manuscript_index():
    files = sorted((ROOT / "manuscript").glob("chapters/*.tex"))  # chapters and appendices
    files = [f for f in files if f.stem not in ("appendix_d_table", "appendix_d_scripts")]  # generated tables of version 1.0.0 (no longer printed): not citations
    return [(f, f.read_text()) for f in files]


def cited(check, index):
    hits = []
    pat = re.compile(r"(?<![\w-])" + re.escape(check) + r"(?![\w])")
    for f, text in index:
        for m in pat.finditer(text):
            labs = list(re.finditer(r"\\label\{([^}]*)\}", text[max(0, m.start() - 4000):m.start()]))
            hits.append((f.stem, labs[-1].group(1) if labs else ""))
    return hits


def main():
    index = manuscript_index()
    rows = []
    for stem, out in sorted(latest_outputs().items()):
        d, _, name = stem.partition("__")
        script = f"verification/{d}/{name}.py"
        for line in out.read_text().splitlines():
            m = LINE.match(line)
            if not m:
                continue
            status, cid, claim = m.groups()
            numeric = bool(re.search(r"float|tol\b|tolerance|numerical|Monte Carlo|quadrature|integration", claim, re.I))
            hits = cited(cid, index)
            chapters = sorted({h[0] for h in hits})
            labels = sorted({h[1] for h in hits if h[1]})
            rows.append([cid, status, script, f"python3 {script}", "numerical" if numeric else "exact", ";".join(chapters),
                         ";".join(labels), str(out.relative_to(ROOT)), claim[:300].replace("\t", " ")])
    head = ["check", "status", "script", "command", "method", "manuscript_files", "manuscript_labels", "run_record", "claim"]
    (ROOT / "verification" / "MANIFEST.tsv").write_text("\t".join(head) + "\n" + "\n".join("\t".join(r) for r in rows) + "\n")
    n = len(rows)
    print(f"{n} checks; {sum(r[1] == 'PASS' for r in rows)} PASS; cited in the manuscript: {sum(bool(r[5]) for r in rows)}; "
          f"with a label: {sum(bool(r[6]) for r in rows)}; numerical: {sum(r[4] == 'numerical' for r in rows)}")


if __name__ == "__main__":
    main()
