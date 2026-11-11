#!/usr/bin/env python3
"""
Audit train/eval/test splits for redundancy and leakage.

Checks (within each split and between every pair of splits):
  * duplicate IDs
  * shared source files (same PDB entry / same UniProt accession in AFDB)
  * shared 50% clusters (the provided clusters-by-entity-50.txt)
  * (optional, after make_fasta.py) identical sequences, using real sequences
  * (optional, --levels 100 95 30) shared RCSB clusters at other identities
    (those files are downloaded on demand)
Also prints the PDB/AFDB composition of each split.

Writes a report to splits/leakage_report.txt; exit code 1 if any
between-split leak is found.

Usage:
  python scripts/check_leakage.py
  python scripts/check_leakage.py --fasta-dir data/fasta
  python scripts/check_leakage.py --levels 100 95 30
"""

import argparse
import sys
from collections import Counter
from itertools import combinations
from pathlib import Path

from common import SPLITS, cluster_labels, download_clusters, entry_of, read_split, source_of


def read_fasta(path: Path) -> dict[str, str]:
    seqs, name = {}, None
    for line in path.read_text().splitlines():
        if line.startswith(">"):
            name = line[1:].split()[0]
            seqs[name] = ""
        elif name:
            seqs[name] += line.strip()
    return seqs


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--splits-dir", default="splits")
    ap.add_argument("--clusters", default="data/raw/clusters-by-entity-50.txt")
    ap.add_argument("--levels", type=int, nargs="*", default=[], choices=[100, 95, 90, 70, 40, 30],
                    help="extra RCSB identity levels to check")
    ap.add_argument("--data-dir", default="data/raw", help="where extra cluster files are downloaded")
    ap.add_argument("--fasta-dir", help="also compare real sequences from <split>.fasta")
    args = ap.parse_args()

    splits = {s: read_split(args.splits_dir, s) for s in SPLITS}
    keys = {"source file": entry_of}
    paths = {50: Path(args.clusters)}
    paths |= {lvl: download_clusters(lvl, Path(args.data_dir))[0] for lvl in args.levels}
    for ident in sorted(paths, reverse=True):
        labels = cluster_labels(paths[ident])
        keys[f"{ident}% cluster"] = lambda i, lab=labels: lab.get(i, i)
    if args.fasta_dir:
        seqs = {}
        for s in SPLITS:
            seqs.update(read_fasta(Path(args.fasta_dir, f"{s}.fasta")))
        keys["identical sequence"] = lambda i: seqs.get(i, i)

    lines = ["Leakage report", "", "Composition:"]
    for s, ids in splits.items():
        comp = Counter(source_of(i) for i in ids)
        lines.append(f"  {s:5s} " + "  ".join(f"{k}={v}" for k, v in sorted(comp.items())))

    lines += ["", "Within split (redundant members = n - number of distinct keys):"]
    for s, ids in splits.items():
        dup_ids = sum(v - 1 for v in Counter(ids).values())
        cols = [f"ids={len(ids)}", f"dup_ids={dup_ids}"]
        cols += [f"{k}={len(ids) - len({f(i) for i in ids})}" for k, f in keys.items()]
        lines.append(f"  {s:5s} " + "  ".join(cols))

    lines += ["", "Between splits (number of shared keys; must be 0):"]
    leaks = 0
    for a, b in combinations(SPLITS, 2):
        n_ids = len(set(splits[a]) & set(splits[b]))
        leaks += n_ids
        cols = [f"ids={n_ids}"]
        for k, f in keys.items():
            n = len({f(i) for i in splits[a]} & {f(i) for i in splits[b]})
            leaks += n
            cols.append(f"{k}={n}")
        lines.append(f"  {a}-{b:5s} " + "  ".join(cols))
    lines += ["", "RESULT: " + ("OK, no leakage between splits" if leaks == 0 else f"LEAKAGE FOUND ({leaks})")]

    report = "\n".join(lines)
    print(report)
    Path(args.splits_dir, "leakage_report.txt").write_text(report + "\n")
    sys.exit(1 if leaks else 0)


if __name__ == "__main__":
    main()
