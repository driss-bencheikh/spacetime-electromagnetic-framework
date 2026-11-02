#!/usr/bin/env python3
"""Run the complete verification suite of the global manuscript and record the results.

    python3 verification/run_all.py --out reproducibility/selected-results/<name> [-j 4] [--only ch07_mass_vacuum]

Runs every verification script of the chapter directories verification/chNN_*/ (helper modules excluded), each in its own process
from the repository root, with the interpreter that runs this file. For each script it writes <out>/<chapter_dir>__<stem>.out (stdout
and stderr) and one line of <out>/summary.tsv: script, exit status, PASS and FAIL counts, lines that are neither, wall time, and whether
the exit status agrees with the counts (exit 0 exactly when there is at least one PASS and no FAIL). It also writes
<out>/environment.txt. The runner reports and does not judge what a PASS establishes (Appendix D of the manuscript).
"""
import argparse
import concurrent.futures as cf
import os
import pathlib
import platform
import re
import subprocess
import sys
import time

ROOT = pathlib.Path(__file__).resolve().parents[1]
HELPERS = {"kernel_space", "tensor_model", "curvature"}  # imported modules, not scripts
LINE = re.compile(r"^(PASS|FAIL) ")
THREADS1 = True  # one BLAS thread per script, so that -j workers do not oversubscribe the cores


def scripts(only=None):
    for d in sorted((ROOT / "verification").glob("ch[0-9][0-9]_*")):
        if only and d.name not in only:
            continue
        for f in sorted(d.glob("*.py")):
            if f.stem not in HELPERS:
                yield d.name + "__", f


def run(prefix, f, out):
    t0 = time.time()
    env = dict(os.environ, OMP_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1", MKL_NUM_THREADS="1") if THREADS1 else None
    p = subprocess.run([sys.executable, str(f)], cwd=ROOT, capture_output=True, text=True, env=env)
    dt = time.time() - t0
    text = p.stdout + (("\n" + p.stderr) if p.stderr else "")
    (out / f"{prefix}{f.stem}.out").write_text(text)
    lines = [l for l in p.stdout.splitlines() if l.strip()]
    n_pass = sum(l.startswith("PASS ") for l in lines)
    n_fail = sum(l.startswith("FAIL ") for l in lines)
    other = sum(not LINE.match(l) for l in lines)
    consistent = (p.returncode == 0) == (n_pass > 0 and n_fail == 0)
    return (str(f.relative_to(ROOT)), p.returncode, n_pass, n_fail, other, round(dt, 1), consistent, bool(p.stderr.strip()))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", required=True)
    ap.add_argument("-j", type=int, default=1)
    ap.add_argument("--only", nargs="*")
    a = ap.parse_args()
    out = ROOT / a.out
    out.mkdir(parents=True, exist_ok=True)
    jobs = list(scripts(a.only))
    with cf.ThreadPoolExecutor(max_workers=a.j) as ex:
        res = list(ex.map(lambda pf: run(pf[0], pf[1], out), jobs))
    head = "script\texit\tPASS\tFAIL\tother_lines\tseconds\tconsistent\tstderr\n"
    (out / "summary.tsv").write_text(head + "".join("\t".join(map(str, r)) + "\n" for r in res))
    mods = []
    for m in ("numpy", "scipy", "sympy", "mpmath"):
        try:
            mods.append(f"{m} {__import__(m).__version__}")
        except Exception as e:  # pragma: no cover
            mods.append(f"{m} unavailable ({e})")
    rev = subprocess.run(["git", "rev-parse", "HEAD"], cwd=ROOT, capture_output=True, text=True).stdout.strip()
    dirty = subprocess.run(["git", "status", "--porcelain", "--", "verification"], cwd=ROOT,
                           capture_output=True, text=True).stdout.strip()
    (out / "environment.txt").write_text(
        f"revision {rev}{' (verification files modified)' if dirty else ''}\npython {platform.python_version()}\n"
        + "\n".join(mods) + f"\nplatform {platform.platform()}\nworkers {a.j}\nBLAS threads per script: {1 if THREADS1 else 'default'}\n")
    tp, tf = sum(r[2] for r in res), sum(r[3] for r in res)
    bad = [r[0] for r in res if r[1] != 0 or r[3] or not r[6]]
    print(f"scripts {len(res)}; PASS {tp}; FAIL {tf}; nonzero exit or inconsistent: {len(bad)}")
    for b in bad:
        print("  ", b)
    return 0 if not bad else 1


if __name__ == "__main__":
    sys.exit(main())
