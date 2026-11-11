#!/usr/bin/env python3
"""
Download mmCIF structures for a split (Python/Windows-friendly version of
download_structures.sh; no wget required). PDB and AFDB IDs can be mixed;
files go to data/structures/<split>/pdb/ and data/structures/<split>/af/.

Usage:
  python scripts/download_structures.py train
  python scripts/download_structures.py test

ID preprocessing:
  PDB:  1ABC_2              -> https://files.rcsb.org/download/1ABC.cif
  AFDB: AF_AFP12345F1_1     -> https://alphafold.ebi.ac.uk/files/AF-P12345-F1-model_v6.cif
"""

import argparse
import sys
import urllib.error
import urllib.request
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

from common import SPLITS, id_to_url, read_split, structure_path


def fetch(url: str, dest: Path) -> tuple[str, str]:
    if dest.exists() and dest.stat().st_size > 0:
        return url, "cached"
    try:
        with urllib.request.urlopen(url, timeout=60) as resp:
            dest.write_bytes(resp.read())
        return url, "ok"
    except (urllib.error.URLError, TimeoutError) as e:
        return url, f"FAILED: {e}"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("split", choices=SPLITS)
    ap.add_argument("--splits-dir", default="splits")
    ap.add_argument("--out-root", default="data/structures")
    ap.add_argument("--workers", type=int, default=8)
    args = ap.parse_args()

    ids = read_split(args.splits_dir, args.split)
    jobs = {id_to_url(i): structure_path(i, args.out_root, args.split) for i in ids if id_to_url(i)}
    for dest in set(jobs.values()):
        dest.parent.mkdir(parents=True, exist_ok=True)
    print(f"[i] {len(ids)} IDs -> {len(jobs)} unique files -> {Path(args.out_root, args.split)}")

    failed = []
    with ThreadPoolExecutor(args.workers) as pool:
        for n, (url, status) in enumerate(pool.map(lambda u: fetch(u, jobs[u]), sorted(jobs)), 1):
            if status.startswith("FAILED"):
                failed.append(f"{url}\t{status}")
            if n % 100 == 0 or n == len(jobs):
                print(f"  {n}/{len(jobs)} (failed: {len(failed)})", flush=True)

    if failed:
        log = Path(args.out_root, args.split, "failed.tsv")
        log.write_text("\n".join(failed) + "\n")
        print(f"[!] {len(failed)} downloads failed, see {log}", file=sys.stderr)


if __name__ == "__main__":
    main()
