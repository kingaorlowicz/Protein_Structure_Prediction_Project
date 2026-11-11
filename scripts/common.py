"""Shared helpers: cluster files (with provenance), ID parsing.

IDs in the RCSB cluster file:
  PDB:  1ABC_2            (PDB entry + entity number)
  AFDB: AF_AFP12345F1_1   (AlphaFold model: UniProt P12345, fragment F1, entity 1)
  MA_*                    (ModelArchive) -- not used
"""

import hashlib
import json
import re
import urllib.request
from datetime import datetime, timezone
from pathlib import Path

CLUSTERS_URL = "https://cdn.rcsb.org/resources/sequence/clusters/clusters-by-entity-{}.txt"
AF_VERSION = "v6"
AF_RE = re.compile(r"AF_AF(.+?)(F\d+)_\d+$")
SPLITS = ["train", "eval", "test"]


def download_clusters(identity: int, data_dir: Path) -> tuple[Path, dict]:
    """Download clusters-by-entity-<identity>.txt once; return (path, provenance)."""
    url = CLUSTERS_URL.format(identity)
    dest = data_dir / f"clusters-by-entity-{identity}.txt"
    info_path = dest.with_suffix(".download.json")
    data_dir.mkdir(parents=True, exist_ok=True)
    if dest.exists() and info_path.exists():
        return dest, json.loads(info_path.read_text())

    print(f"[i] Downloading {url}")
    with urllib.request.urlopen(url) as resp:
        data = resp.read()
        headers = resp.headers
    dest.write_bytes(data)
    info = {
        "url": url,
        "access_date_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "last_modified": headers.get("Last-Modified"),
        "etag": headers.get("ETag"),
        "md5": hashlib.md5(data).hexdigest(),
        "n_bytes": len(data),
    }
    info_path.write_text(json.dumps(info, indent=2))
    return dest, info


def file_info(path: Path) -> dict:
    """Provenance of a local cluster file (the one provided by the team)."""
    path = Path(path)
    data = path.read_bytes()
    info = {"path": str(path), "md5": hashlib.md5(data).hexdigest(), "n_bytes": len(data),
            "mtime_utc": datetime.fromtimestamp(path.stat().st_mtime, timezone.utc).isoformat(timespec="seconds")}
    sidecar = path.with_suffix(".download.json")
    if sidecar.exists():
        info["download"] = json.loads(sidecar.read_text())
    return info


def read_clusters(path: Path) -> list[list[str]]:
    with open(path) as fh:
        return [line.split() for line in fh]


def cluster_labels(path: Path) -> dict[str, int]:
    """member ID -> cluster index (line number)."""
    labels = {}
    with open(path) as fh:
        for i, line in enumerate(fh):
            for member in line.split():
                labels[member] = i
    return labels


def source_of(member_id: str) -> str:
    """'af' (AlphaFold DB), 'ma' (ModelArchive) or 'pdb' (experimental)."""
    if member_id.startswith("AF_"):
        return "af"
    if member_id.startswith("MA_"):
        return "ma"
    return "pdb"


def keep_member(member_id: str, source: str) -> bool:
    """source: 'pdb', 'af' or 'all' (= pdb + af, never MA_)."""
    s = source_of(member_id)
    return s != "ma" and (source == "all" or s == source)


def entry_of(member_id: str) -> str:
    """Physical source file of a member: PDB entry (1ABC_2 -> pdb:1ABC) or
    UniProt accession for AFDB (AF_AFP12345F1_1 -> af:P12345, all fragments)."""
    if source_of(member_id) == "af":
        m = AF_RE.match(member_id)
        return f"af:{m.group(1) if m else member_id}"
    return f"pdb:{member_id.rsplit('_', 1)[0]}"


def id_to_url(member_id: str) -> str | None:
    if source_of(member_id) == "pdb":
        return f"https://files.rcsb.org/download/{member_id.rsplit('_', 1)[0]}.cif"
    m = AF_RE.match(member_id)
    if m:
        return f"https://alphafold.ebi.ac.uk/files/AF-{m.group(1)}-{m.group(2)}-model_{AF_VERSION}.cif"
    return None


def structure_path(member_id: str, root: Path, split: str) -> Path:
    """Where download_structures.py stores the file: <root>/<split>/<pdb|af>/<name>.cif"""
    return Path(root, split, source_of(member_id), id_to_url(member_id).rsplit("/", 1)[1])


def read_split(splits_dir, split: str) -> list[str]:
    return Path(splits_dir, f"{split}.txt").read_text().split()
