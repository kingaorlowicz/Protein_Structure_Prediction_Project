#!/usr/bin/env python3
"""
Build the FASTA file (model input) for a split (PDB and AFDB IDs mixed).

The sequence is the full canonical entity sequence from the original mmCIF
(_entity_poly.pdbx_seq_one_letter_code_can), not only the residues that were
resolved in the crystal: that is what the model gets as input at test time.

Usage:
  python scripts/make_fasta.py train      # -> data/fasta/train.fasta
  python scripts/make_fasta.py test
"""

import argparse
from pathlib import Path

from Bio.PDB.MMCIF2Dict import MMCIF2Dict

from common import SPLITS, read_split, source_of, structure_path


def as_list(x):
    return x if isinstance(x, list) else [x]


def entity_sequence(cif: Path, member_id: str) -> str | None:
    d = MMCIF2Dict(str(cif))
    ids = as_list(d["_entity_poly.entity_id"])
    seqs = as_list(d["_entity_poly.pdbx_seq_one_letter_code_can"])
    entity = member_id.rsplit("_", 1)[1] if source_of(member_id) == "pdb" else ids[0]
    for e, s in zip(ids, seqs):
        if e == entity:
            return "".join(s.split())
    return None


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("split", choices=SPLITS)
    ap.add_argument("--splits-dir", default="splits")
    ap.add_argument("--structures-root", default="data/structures")
    ap.add_argument("--out-dir", default="data/fasta")
    args = ap.parse_args()

    ids = read_split(args.splits_dir, args.split)
    out_path = Path(args.out_dir, f"{args.split}.fasta")
    out_path.parent.mkdir(parents=True, exist_ok=True)

    missing = []
    with open(out_path, "w") as out:
        for member_id in ids:
            cif = structure_path(member_id, args.structures_root, args.split)
            seq = entity_sequence(cif, member_id) if cif.exists() else None
            if seq:
                out.write(f">{member_id} source={source_of(member_id)} len={len(seq)}\n{seq}\n")
            else:
                missing.append(member_id)
    print(f"[+] {len(ids) - len(missing)}/{len(ids)} sequences -> {out_path}")
    if missing:
        print(f"[!] missing/failed: {len(missing)} (e.g. {missing[:5]})")


if __name__ == "__main__":
    main()
