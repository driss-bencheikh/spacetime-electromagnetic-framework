#!/usr/bin/env python3
"""Generate the theorem-to-code table described in Appendix D (Phase 9F, 2 November 2026; not printed since version 1.0.1).

    python3 verification/make_appendix_d.py

Inputs: verification/MANIFEST.tsv (every check, its script, method and the manuscript labels whose statement cites it),
verification/COVERAGE.tsv (every labelled statement of the manuscript with its provenance and whether a written proof is present),
the audit checks (verification/chNN_*/a<k>_*.py: Phase 9D, the v1.0.0 release audit and the v1.0.1 corrections audit; their claim text names the labels they
test), and
reproducibility/selected-results/uncited_related_checks.tsv (checks that test a statement without being cited by it).
Outputs: verification/THEOREM_CODE.tsv (label, chapter, environment, written proof, cited, audit and related checks, evidence, scripts),
verification/SCRIPTS.tsv (script, exact and numerical checks) and a printed coverage summary.

The table records evidence, not proof: a check listed against a label tests the coded assertion stated in its claim text.
"""
import csv
import pathlib
import re
from collections import defaultdict

ROOT = pathlib.Path(__file__).resolve().parents[2]
csv.field_size_limit(10 ** 7)


def read_tsv(p):
    with open(p, newline="") as f:
        return list(csv.DictReader(f, delimiter="\t"))


manifest = read_tsv(ROOT / "verification" / "MANIFEST.tsv")
coverage = read_tsv(ROOT / "verification" / "COVERAGE.tsv")
related = read_tsv(ROOT / "reproducibility" / "selected-results" / "uncited_related_checks.tsv")
labels = [r for r in coverage if r["label"] and not r["label"].startswith("(")]
known = {r["label"] for r in labels}

method = {r["check"]: r["method"] for r in manifest}
script = {r["check"]: r["script"] for r in manifest}
status = {r["check"]: r["status"] for r in manifest}
ev = defaultdict(lambda: {"cited": set(), "audit": set(), "related": set()})
for r in manifest:
    for lab in filter(None, r["manuscript_labels"].split(";")):
        ev[lab]["cited"].add(r["check"])
    if r["check"].startswith("a") and re.match(r"a\d+-", r["check"]):
        head = r["claim"]
        for lab in known:
            if re.search(r"(?<![\w:])" + re.escape(lab) + r"(?!\w|:\w)", head):
                ev[lab]["audit"].add(r["check"])
for r in related:
    ids = (r.get("checks naming the label", "") + " " + r.get("checks citing its 26 October theorem number", "")).split()
    for c in ids:
        if c in method:
            ev[r["label"]]["related"].add(c)

rows = []
for r in labels:
    lab = r["label"]
    e = ev.get(lab, {"cited": set(), "audit": set(), "related": set()})
    checks = sorted(e["cited"] | e["audit"] | e["related"])
    meth = sorted({method[c] for c in checks})
    if checks:
        kind = "exact" if meth == ["exact"] else ("numerical" if meth == ["numerical"] else "exact+numerical")
    elif r["written_proof"] == "yes":
        kind = "proof in text"
    elif r["cites_literature"] == "yes" or "inherited" in r["markers"]:
        kind = "inherited"
    elif "conjecture" in r["markers"] or r["environment"] == "conjecture":
        kind = "conjecture"
    else:
        kind = "argument in text"
    rows.append({"label": lab, "chapter": r["file"], "environment": r["environment"], "written_proof": r["written_proof"],
                 "cited_checks": " ".join(sorted(e["cited"])), "audit_checks": " ".join(sorted(e["audit"])),
                 "related_checks": " ".join(sorted(e["related"])), "evidence": kind,
                 "scripts": " ".join(sorted({script[c].rsplit("/", 1)[-1] for c in checks}))})

with open(ROOT / "verification" / "THEOREM_CODE.tsv", "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=list(rows[0]), delimiter="\t")
    w.writeheader()
    w.writerows(rows)


# scripts table
by_script = defaultdict(lambda: [0, 0])
for r in manifest:
    by_script[r["script"]][0 if r["method"] == "exact" else 1] += 1
with open(ROOT / "verification" / "SCRIPTS.tsv", "w", newline="") as f:
    f.write("script\texact\tnumerical\n")
    for s_, (a, b) in sorted(by_script.items()):
        f.write(f"{s_}\t{a}\t{b}\n")

n = len(rows)
cnt = defaultdict(int)
for r in rows:
    cnt[r["evidence"]] += 1
with_check = sum(1 for r in rows if r["cited_checks"] or r["audit_checks"] or r["related_checks"])
print(f"{n} labelled statements; with a check: {with_check}; by evidence: {dict(cnt)}")
print(f"checks: {len(manifest)} ({sum(1 for r in manifest if r['status'] == 'PASS')} PASS); scripts: {len(by_script)}")
