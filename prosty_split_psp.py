#!/usr/bin/env python3
"""
Prosta wersja zadania z Meeting 2 (PSP) - jeden plik.

Kroki:
  1. wczytanie klastrów 50% (clusters-by-entity-50.txt; pobiera, jeśli brak)
  2. wybór źródła: SOURCE = "pdb" (Team 1), "af" (Team 2) albo "mixed"
     (train/eval: PDB + AFDB, test: tylko PDB)
  3. klastry z >= 10 członkami, tasowanie (SEED), 10 losowych członków
     + filtr redundancji (DEDUP_100): ten sam plik .cif ani ta sama sekwencja
       (klastry 100%, np. model AF i struktura PDB tego samego białka) nie
       występuje w zbiorze dwa razy
  4. całe klastry do splitów: 4000 / 500 / 500 białek
  5. zapis ID: train_<src>.txt, eval_<src>.txt, test_<src>.txt
     + metadata_<src>.json (seed, wersja AFDB, md5 i data dostępu pliku klastrów)
  6. pobranie plików .cif (structures/<split>/<pdb|af>/)
  7. wycięcie łańcucha dla ID z PDB (chains/<split>/)

Uruchomienie (z folderu, w którym jest ten plik):
  python prosty_split_psp.py               # tylko podział
  python prosty_split_psp.py --download    # + pobranie struktur
  python prosty_split_psp.py --download --chains
Wymaga: biopython (tylko do --chains)
"""

import argparse
import hashlib
import json
import random
import re
import urllib.error
import urllib.request
from concurrent.futures import ThreadPoolExecutor
import time
from datetime import datetime, timezone
from email.utils import parsedate_to_datetime
from pathlib import Path

SOURCE = "mixed"          # "pdb" | "af" | "mixed"
SEED = 42
K = 10                    # członków na klaster
SIZES = {"test": 500, "eval": 500, "train": 4000}   # kolejność wypełniania
CLUSTERS_URL = "https://cdn.rcsb.org/resources/sequence/clusters/clusters-by-entity-50.txt"
CLUSTERS_FILE = Path("clusters-by-entity-50.txt")
AF_RE = re.compile(r"AF_AF(.+?)(F\d+)_\d+$")
AF_VERSION = "v6"

# Filtr redundancji: z jednej grupy identycznych sekwencji (klastry 100%) bierzemy
# najwyżej jedno białko w całym zbiorze -> brak par PDB+AF / PDB+PDB tego samego białka.
# False = podział bez filtra (identyczny z 2wersja).
DEDUP_100 = True
CLUSTERS100_URL = "https://cdn.rcsb.org/resources/sequence/clusters/clusters-by-entity-100.txt"
CLUSTERS100_FILE = Path("clusters-by-entity-100.txt")

# Pliki klastrów dołączone do projektu (żeby dało się odtworzyć dokładnie ten podział)
KNOWN_CLUSTERS = {
    "058f24f6704227377e7d0a5e6ba4305c": {      # clusters-by-entity-50.txt (wydanie 2026-11-09)
        "access_date_utc": "2026-11-11T22:46:47+00:00",
        "last_modified": "Mon, 09 Nov 2026 18:21:29 GMT",
    },
    "e61dccb7a61e0f43d62590b79702c978": {      # stary clusters-by-entity-50.txt (2026-10-24), kopia w stary_split_2026-10-24/
        "access_date_utc": "2026-10-29T23:42:05+00:00",
        "last_modified": "Sat, 24 Oct 2026 23:19:49 GMT",
    },
    "9a41c9048d60900c512a69552fe57cd1": {      # clusters-by-entity-100.txt
        "access_date_utc": "2026-11-11T14:32:51+00:00",
        "last_modified": "Mon, 09 Nov 2026 18:23:37 GMT",
    },
}


def src(i):                       # skąd pochodzi ID
    return "af" if i.startswith("AF_") else "ma" if i.startswith("MA_") else "pdb"


def allowed(i, split):            # czy ID pasuje do źródła danego splitu
    s = src(i)
    if SOURCE == "pdb":
        return s == "pdb"
    if SOURCE == "af":
        return s == "af"
    return s == "pdb" if split == "test" else s in ("pdb", "af")   # mixed


def file_key(i):                  # fizyczny plik .cif, z którego pochodzi ID
    if src(i) == "af":
        m = AF_RE.match(i)
        return "af:" + (m.group(1) if m else i)
    return "pdb:" + i.rsplit("_", 1)[0]       # 1ABC_1 i 1ABC_2 -> ten sam plik 1ABC.cif


def url_of(i):
    if src(i) == "pdb":
        return f"https://files.rcsb.org/download/{i.rsplit('_', 1)[0]}.cif"
    m = AF_RE.match(i)
    return f"https://alphafold.ebi.ac.uk/files/AF-{m.group(1)}-{m.group(2)}-model_{AF_VERSION}.cif" if m else None


# ---------------------------------------------------------------- 1-5: podział
def load_clusters(url, path):
    """Wczytuje (i w razie braku pobiera) plik klastrów; zwraca (klastry, metadane)."""
    provenance = {"url": url}
    if not path.exists():
        print("Pobieram", url)
        with urllib.request.urlopen(url) as resp:
            path.write_bytes(resp.read())
            provenance["last_modified"] = resp.headers.get("Last-Modified")
        provenance["access_date_utc"] = datetime.now(timezone.utc).isoformat(timespec="seconds")
    data = path.read_bytes()
    md5 = hashlib.md5(data).hexdigest()
    provenance.update(KNOWN_CLUSTERS.get(md5, {}))
    provenance.update({"md5": md5, "n_bytes": len(data)})
    if "access_date_utc" not in provenance:
        print(f"[!] Nieznany plik {path} (inny md5) - data dostępu nieznana, podział będzie inny")
    return [line.split() for line in data.decode().splitlines()], provenance


def make_splits():
    # Pliki 50% i 100% muszą pochodzić z tego samego wydania RCSB -> jeśli brakuje
    # jednego, pobieramy oba na nowo.
    if DEDUP_100 and CLUSTERS_FILE.exists() != CLUSTERS100_FILE.exists():
        print("[!] Brak jednego z plików klastrów - pobieram oba z tego samego wydania")
        CLUSTERS_FILE.unlink(missing_ok=True)
        CLUSTERS100_FILE.unlink(missing_ok=True)

    clusters, provenance = load_clusters(CLUSTERS_URL, CLUSTERS_FILE)

    seq_group = {}                                    # ID -> numer grupy identycznych sekwencji
    provenance100 = None
    if DEDUP_100:
        clusters100, provenance100 = load_clusters(CLUSTERS100_URL, CLUSTERS100_FILE)
        lm50, lm100 = provenance.get("last_modified"), provenance100.get("last_modified")
        if lm50 and lm100:
            gap = abs((parsedate_to_datetime(lm50) - parsedate_to_datetime(lm100)).total_seconds())
            if gap > 24 * 3600:
                print(f"[!] Pliki klastrów z różnych wydań RCSB: 50% = {lm50}, 100% = {lm100}")
        for n, line in enumerate(clusters100):
            for m in line:
                seq_group[m] = n
        del clusters100

    def keys(m):                                      # wszystko, co nie może się powtórzyć
        k = {file_key(m)}
        if DEDUP_100:
            k.add(f"seq:{seq_group.get(m, m)}")       # ID spoza pliku 100% = własna grupa
        return k

    rng = random.Random(SEED)
    rng.shuffle(clusters)

    need = {s: n // K for s, n in SIZES.items()}      # 50 / 50 / 400 klastrów
    splits = {s: [] for s in SIZES}
    used = set()            # użyte pliki .cif (i grupy identycznych sekwencji, jeśli DEDUP_100)

    for c in clusters:
        for split in SIZES:                           # test -> eval -> train
            if need[split] == 0:
                continue
            members = [m for m in c if allowed(m, split)]
            if len(members) < K:                      # klaster z >= 10 członkami
                continue
            rng.shuffle(members)
            picked, taken = [], set()
            for m in members:
                k = keys(m)
                if not (k & used) and not (k & taken):
                    picked.append(m)
                    taken |= k
                if len(picked) == K:
                    break
            if len(picked) == K:                      # cały klaster -> jeden split
                splits[split] += picked
                used |= taken
                need[split] -= 1
                break
        if not any(need.values()):
            break

    for s, ids in splits.items():
        Path(f"{s}_{SOURCE}.txt").write_text("\n".join(ids) + "\n")
        n_af = sum(src(i) == "af" for i in ids)
        print(f"{s}_{SOURCE}.txt: {len(ids)} białek (PDB={len(ids) - n_af}, AF={n_af})")

    meta = {
        "run_date_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "seed": SEED,
        "source": SOURCE,
        "members_per_cluster": K,
        "split_sizes": SIZES,
        "alphafold_db_version": AF_VERSION,
        "clusters_file": provenance,
        "dedup_identical_sequences": DEDUP_100,
        "clusters100_file": provenance100,
        "composition": {s: {"pdb": sum(src(i) == "pdb" for i in ids), "af": sum(src(i) == "af" for i in ids)}
                        for s, ids in splits.items()},
    }
    Path(f"metadata_{SOURCE}.json").write_text(json.dumps(meta, indent=2))
    print(f"metadata_{SOURCE}.json: seed, wersje i data dostępu do danych")
    return splits


# ---------------------------------------------------------------- 6: pobranie
def download(splits):
    def fetch(job):
        url, dest = job
        if dest.exists():
            return True
        for attempt in range(6):                      # ponawianie po 429 / timeoucie
            try:
                data = urllib.request.urlopen(url, timeout=60).read()
                tmp = dest.with_suffix(".part")       # niepełny plik nie udaje gotowego
                tmp.write_bytes(data)
                tmp.replace(dest)
                return True
            except urllib.error.HTTPError as e:
                if e.code == 404:
                    print("BŁĄD 404", url)
                    return False
                err = e
            except Exception as e:
                err = e
            time.sleep(5 * 2 ** attempt)              # 5, 10, 20, 40, 80, 160 s
        print("BŁĄD", url, err)
        return False

    for s, ids in splits.items():
        jobs = {(url_of(i), Path("structures", s, src(i), url_of(i).rsplit("/", 1)[1])) for i in ids if url_of(i)}
        for _, d in jobs:
            d.parent.mkdir(parents=True, exist_ok=True)
        with ThreadPoolExecutor(3) as pool:           # mniej wątków = mniej 429 z AlphaFold DB
            ok = sum(pool.map(fetch, sorted(jobs)))
        print(f"{s}: pobrano {ok}/{len(jobs)} plików")


# ---------------------------------------------------------------- 7: łańcuchy
def extract_chains(splits):
    from Bio.PDB import MMCIFIO, MMCIFParser, Select
    from Bio.PDB.MMCIF2Dict import MMCIF2Dict

    for s, ids in splits.items():
        out = Path("chains", s)
        out.mkdir(parents=True, exist_ok=True)
        for i in ids:
            if src(i) != "pdb":                       # model AF = już jeden łańcuch
                continue
            pdb_id, entity = i.rsplit("_", 1)
            cif = Path("structures", s, "pdb", f"{pdb_id}.cif")
            if not cif.exists():
                continue
            d = MMCIF2Dict(str(cif))
            # Biopython używa auth_asym_id jako chain.id -> mapujemy encję przez _atom_site
            chains = [a for a, e in zip(d["_atom_site.auth_asym_id"], d["_atom_site.label_entity_id"]) if e == entity]
            if not chains:
                continue
            keep = chains[0]                          # jedna kopia łańcucha

            class Keep(Select):
                def accept_model(self, model): return model.id == 0
                def accept_chain(self, chain): return chain.id == keep

            io = MMCIFIO()
            io.set_structure(MMCIFParser(QUIET=True).get_structure(pdb_id, str(cif)))
            io.save(str(out / f"{i}.cif"), Keep())
        print(f"{s}: łańcuchy zapisane w {out}")


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--download", action="store_true")
    ap.add_argument("--chains", action="store_true")
    args = ap.parse_args()
    sp = make_splits()
    if args.download:
        download(sp)
    if args.chains:
        extract_chains(sp)
