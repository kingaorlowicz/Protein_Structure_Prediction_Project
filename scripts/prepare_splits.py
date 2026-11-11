#!/usr/bin/env python3
"""
Train/eval/test split for Protein Structure Prediction (PSP).
One combined dataset: train/eval mix PDB + AFDB, test is PDB only.

Input: only data/raw/clusters-by-entity-50.txt (--clusters), RCSB clusters at
50% sequence identity, one cluster per line.

Procedure (wg "Trainevaltest split preparation.txt", adapted):
  1. Read the 50% clusters; drop MA_ (ModelArchive) members.
  2. Shuffle clusters (fixed seed).
  3. Fill splits in order test -> eval -> train with whole clusters:
       test        : members from --test-source  (default: pdb)
       eval, train : members from --train-source (default: all = pdb + af);
                     --eval-source pdb makes eval look like test
     A cluster goes to the first split that still needs clusters and in which
     it has >= 10 usable members; 10 random members are picked.
  4. Usable member = its source file is not used anywhere else in the dataset
     (PDB entry 1ABC for 1ABC_1/1ABC_2, UniProt accession for AFDB F1/F2
     fragments). So no identical .cif file appears twice in a split or in two
     splits, and no 50% cluster is shared between splits.
  5. 400 / 50 / 50 clusters -> 4000 / 500 / 500 proteins.
  6. Save splits/train.txt, eval.txt, test.txt (one ID per line) + metadata.json.

Usage:
  python scripts/prepare_splits.py
  python scripts/prepare_splits.py --train-source all --test-source pdb --seed 42
"""

import argparse
import json
import random
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path

from common import SPLITS, entry_of, file_info, keep_member, read_clusters, source_of

SEED = 42
MEMBERS_PER_CLUSTER = 10
SPLIT_SIZES = {"test": 500, "eval": 500, "train": 4000}  # in proteins; filled in this order


def pick_members(members, rng, used_entries, k):
    """Pick k members from distinct, not yet used source files.
    Returns None (and reserves nothing) if not possible."""
    members = members[:]
    rng.shuffle(members)
    picked, entries = [], set()
    for m in members:
        e = entry_of(m)
        if e in entries or e in used_entries:
            continue
        picked.append(m)
        entries.add(e)
        if len(picked) == k:
            used_entries |= entries
            return picked
    return None


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--clusters", default="data/raw/clusters-by-entity-50.txt")
    ap.add_argument("--train-source", choices=["all", "pdb", "af"], default="all")
    ap.add_argument("--eval-source", choices=["all", "pdb", "af"], help="default: same as --train-source")
    ap.add_argument("--test-source", choices=["all", "pdb", "af"], default="pdb")
    ap.add_argument("--seed", type=int, default=SEED)
    ap.add_argument("--out-dir", default="splits")
    args = ap.parse_args()

    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    source = {"test": args.test_source, "eval": args.eval_source or args.train_source, "train": args.train_source}
    targets = {s: n // MEMBERS_PER_CLUSTER for s, n in SPLIT_SIZES.items()}

    rng = random.Random(args.seed)
    clusters = read_clusters(Path(args.clusters))
    rng.shuffle(clusters)

    used_entries = set()
    cluster_log = {s: [] for s in SPLIT_SIZES}
    for c in clusters:
        for split in SPLIT_SIZES:
            if len(cluster_log[split]) == targets[split]:
                continue
            members = [m for m in c if keep_member(m, source[split])]
            if len(members) < MEMBERS_PER_CLUSTER:
                continue
            picked = pick_members(members, rng, used_entries, MEMBERS_PER_CLUSTER)
            if picked:
                cluster_log[split].append({"cluster_size": len(c), "representative": c[0], "members": picked})
                break
        if all(len(cluster_log[s]) == targets[s] for s in SPLIT_SIZES):
            break
    else:
        missing = {s: targets[s] - len(cluster_log[s]) for s in SPLIT_SIZES}
        raise SystemExit(f"Not enough clusters, still missing: {missing}")

    splits = {s: [m for c in cluster_log[s] for m in c["members"]] for s in SPLITS}

    # Sanity checks (full audit: scripts/check_leakage.py)
    all_ids = [i for ids in splits.values() for i in ids]
    assert len(all_ids) == len(set(all_ids)), "duplicate IDs"
    assert len({entry_of(i) for i in all_ids}) == len(all_ids), "shared source files"
    for s, ids in splits.items():
        assert all(keep_member(i, source[s]) for i in ids), f"wrong source in {s}"

    composition = {}
    for s, ids in splits.items():
        composition[s] = dict(Counter(source_of(i) for i in ids))
        (out_dir / f"{s}.txt").write_text("\n".join(ids) + "\n")
        parts = ", ".join(f"{k}={v}" for k, v in sorted(composition[s].items()))
        print(f"[+] {out_dir / f'{s}.txt'}: {len(ids)} proteins, {len(cluster_log[s])} clusters ({parts})")

    meta = {
        "seed": args.seed,
        "sources": source,
        "members_per_cluster": MEMBERS_PER_CLUSTER,
        "split_sizes": SPLIT_SIZES,
        "composition": composition,
        "clusters_file": file_info(Path(args.clusters)),
        "created_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "clusters": {s: cluster_log[s] for s in SPLITS},
    }
    (out_dir / "metadata.json").write_text(json.dumps(meta, indent=2))
    print(f"[+] {out_dir / 'metadata.json'}")


if __name__ == "__main__":
    main()
