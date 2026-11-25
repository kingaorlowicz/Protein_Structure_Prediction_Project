import json
import urllib.request
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path


#train_mixed.txt, eval_mixed.txt, test_mixed.txt
INPUT_DIR = Path(".")
CACHE_DIR = Path("sequence_check/cache")
OUTPUT_DIR = Path("fasta")

WORKERS = 4


def get_sequence(identifier):
    """Odczytuje sekwencję z cache lub pobiera ją z RCSB."""
    cache_file = CACHE_DIR / f"{identifier}.json"

    if cache_file.exists():
        record = json.loads(cache_file.read_text(encoding="utf-8"))

        if record["id"] != identifier:
            raise ValueError(f"Nieprawidłowe ID w cache: {identifier}")

        sequence = record["sequence"]

    else:
        # ID ma identyfikator wpisu oraz numer encji
        entry, entity = identifier.rsplit("_", 1)
        url = (
            "https://data.rcsb.org/rest/v1/core/polymer_entity/"
            f"{entry}/{entity}"
        )

        with urllib.request.urlopen(url, timeout=30) as response:
            data = json.load(response)

        if data["rcsb_id"].upper() != identifier.upper():
            raise ValueError(f"Nieprawidłowe ID w odpowiedzi: {identifier}")

        sequence = data["entity_poly"]["pdbx_seq_one_letter_code_can"]

    sequence = "".join(sequence.split()).upper()
    if not sequence or any(letter not in "ABCDEFGHIJKLMNOPQRSTUVWXYZ" for letter in sequence):
        raise ValueError(f"Nieprawidłowa sekwencja: {identifier}")

    if not cache_file.exists():
        record = {"id": identifier, "sequence": sequence}
        temporary = cache_file.with_suffix(".part")
        temporary.write_text(json.dumps(record), encoding="utf-8")
        temporary.replace(cache_file)

    return sequence


def main():
    CACHE_DIR.mkdir(parents=True, exist_ok=True)
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    prepared = {}

    for split in ("train", "eval", "test"):
        input_file = INPUT_DIR / f"{split}_mixed.txt"
        ids = input_file.read_text(encoding="utf-8-sig").split()

        if not ids:
            raise ValueError(f"Pusta lista ID: {input_file}")

        print(f"{split}: odczytuję sekwencje dla {len(ids)} białek...", flush=True)

        #dla przyspieszenia - równolegle pobieramy sekwencje z RCSB
        with ThreadPoolExecutor(max_workers=WORKERS) as pool:
            sequences = list(pool.map(get_sequence, ids))

        prepared[split] = list(zip(ids, sequences))

    for split, records in prepared.items():
        output_file = OUTPUT_DIR / f"{split}.fasta"

        with output_file.open("w", encoding="utf-8") as output:
            for identifier, sequence in records:
                output.write(f">{identifier}\n{sequence}\n")

        print(f"Zapisano {len(records)} sekwencji do {output_file}")


if __name__ == "__main__":
    main()