#!/usr/bin/env python3
"""
Run the whole data pipeline (combined PDB + AFDB dataset), step by step.

Steps: split -> check -> download -> chains (PDB IDs only) -> fasta -> check_fasta

Usage:
  python scripts/run_pipeline.py                         # everything
  python scripts/run_pipeline.py --steps split check     # only the split + audit
  python scripts/run_pipeline.py --steps download --splits test
"""

import argparse
import subprocess
import sys

from common import SPLITS

STEPS = ["split", "check", "download", "chains", "fasta", "check_fasta"]


def run(*args):
    cmd = [sys.executable, *args]
    print(f"\n>>> {' '.join(cmd)}", flush=True)
    subprocess.run(cmd, check=True)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--steps", nargs="+", choices=STEPS, default=STEPS)
    ap.add_argument("--splits", nargs="+", choices=SPLITS, default=SPLITS)
    ap.add_argument("--workers", default="8")
    args = ap.parse_args()

    for step in STEPS:
        if step not in args.steps:
            continue
        if step == "split":
            run("scripts/prepare_splits.py")
        elif step == "check":
            run("scripts/check_leakage.py")
        elif step == "download":
            for s in args.splits:
                run("scripts/download_structures.py", s, "--workers", args.workers)
        elif step == "chains":
            for s in args.splits:
                run("scripts/extract_chain.py", "-", f"data/structures/{s}/pdb", f"data/chains/{s}",
                    "--ids-file", f"splits/{s}.txt", "--workers", args.workers)
        elif step == "fasta":
            for s in args.splits:
                run("scripts/make_fasta.py", s)
        elif step == "check_fasta":
            run("scripts/check_leakage.py", "--fasta-dir", "data/fasta")


if __name__ == "__main__":
    main()
