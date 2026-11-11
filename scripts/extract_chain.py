#!/usr/bin/env python3
"""
Team 1 (PDB): cut each entry (PDBID_ENTITY) out of the full mmCIF file.

By default only ONE chain of the entity (first in the file) and only the
first model are kept: a homo-oligomer would otherwise give several identical
copies of the same chain in one training example.

Usage (single entry, as in the instructions):
  python scripts/extract_chain.py 1ABC_2 data/structures/train/pdb data/chains/train

Usage (whole split, parallel; AFDB IDs in the list are skipped):
  python scripts/extract_chain.py - data/structures/train/pdb data/chains/train --workers 8 < splits/train.txt
  python scripts/extract_chain.py - data/structures/train/pdb data/chains/train --ids-file splits/train.txt

Changes vs. the original snippet:
  * entity -> chain mapping uses _atom_site (label_entity_id -> auth_asym_id),
    because Biopython's chain.id is the *author* chain ID, while
    _struct_asym.id is the *label* ID (they often differ).
  * one chain + first model only (--all-chains keeps every copy).
  * accepts a list of IDs on stdin and processes them in parallel.
"""

import argparse
import sys
from concurrent.futures import ProcessPoolExecutor
from functools import partial
from pathlib import Path

from Bio.PDB import MMCIFIO, MMCIFParser, Select
from Bio.PDB.MMCIF2Dict import MMCIF2Dict


class Keep(Select):
    def __init__(self, chains):
        self.chains = chains

    def accept_model(self, model):
        return model.id == 0

    def accept_chain(self, chain):
        return chain.id in self.chains


def extract(entry_id: str, in_dir: str, out_dir: str, all_chains: bool = False) -> str:
    out_path = Path(out_dir, f"{entry_id}.cif")
    if out_path.exists():
        return f"{entry_id}\tcached"
    pdb_id, entity_id = entry_id.rsplit("_", 1)
    cif = Path(in_dir, f"{pdb_id}.cif")
    if not cif.exists():
        return f"{entry_id}\tMISSING {cif}"
    try:
        d = MMCIF2Dict(str(cif))
        # dict.fromkeys keeps file order
        chains = list(dict.fromkeys(a for a, e in zip(d["_atom_site.auth_asym_id"], d["_atom_site.label_entity_id"])
                                    if e == entity_id))
        if not chains:
            return f"{entry_id}\tNO CHAINS for entity {entity_id}"
        keep = set(chains if all_chains else chains[:1])
        st = MMCIFParser(QUIET=True).get_structure(pdb_id, str(cif))
        io = MMCIFIO()
        io.set_structure(st)
        io.save(str(out_path), Keep(keep))
        return f"{entry_id}\tok ({','.join(sorted(keep))})"
    except Exception as e:  # noqa: BLE001 - log and continue with the rest
        return f"{entry_id}\tERROR {e}"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("entry_id", help="PDBID_ENTITY, or '-' to read IDs from stdin")
    ap.add_argument("in_dir")
    ap.add_argument("out_dir")
    ap.add_argument("--workers", type=int, default=4)
    ap.add_argument("--ids-file", help="file with one ID per line (instead of stdin; works in PowerShell)")
    ap.add_argument("--all-chains", action="store_true", help="keep every chain of the entity")
    args = ap.parse_args()

    Path(args.out_dir).mkdir(parents=True, exist_ok=True)
    if args.entry_id != "-":
        ids = [args.entry_id]
    else:
        ids = (Path(args.ids_file).read_text() if args.ids_file else sys.stdin.read()).split()
    ids = [i for i in ids if not i.startswith(("AF_", "MA_"))]  # AFDB models are single-chain already

    job = partial(extract, in_dir=args.in_dir, out_dir=args.out_dir, all_chains=args.all_chains)
    with ProcessPoolExecutor(args.workers) as pool:
        results = list(pool.map(job, ids))

    problems = [r for r in results if not (r.endswith("cached") or "\tok" in r)]
    for r in problems:
        print(r, file=sys.stderr)
    print(f"[i] {len(ids) - len(problems)}/{len(ids)} entries extracted to {args.out_dir}")
    if problems:
        Path(args.out_dir, "problems.tsv").write_text("\n".join(problems) + "\n")


if __name__ == "__main__":
    main()
