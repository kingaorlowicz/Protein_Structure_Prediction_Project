from pathlib import Path
import skbio
import torch
from torch import nn
from torch.nn import functional as F

FASTA_FILE = Path("fasta/test.fasta")
OUTPUT_DIR = Path("dummy_test")

BATCH_SIZE = 16
MAX_LENGTH = 1024
LIMIT = 5               # None = wszystkie białka
SEED = 42


# Słownik: 20 aminokwasów + X + PAD
STANDARD = "ACDEFGHIKLMNPQRSTVWY"
VOCAB = list(STANDARD) + ["X", "PAD"]
AA_TO_INDEX = {aa: i for i, aa in enumerate(VOCAB)}

# Nazwy aminokwasów używane w pliku PDB
AA3 = dict(zip(STANDARD, (
    "ALA CYS ASP GLU PHE GLY HIS ILE LYS LEU MET ASN "
    "PRO GLN ARG SER THR VAL TRP TYR"
).split()))
AA3["X"] = "UNK"


def read_fasta(path):
    """Wczytuje ID i sekwencje białek."""
    records = []

    for record in skbio.io.read(str(path), format="fasta"):
        protein_id = record.metadata["id"]
        sequence = str(record).upper()

        if not sequence:
            raise ValueError(f"Pusta sekwencja: {protein_id}")

        # Niestandardowe aminokwasy zamieniamy na X.
        sequence = "".join(
            aa if aa in STANDARD else "X"
            for aa in sequence
        )

        records.append((protein_id, sequence))

    if not records:
        raise ValueError("Plik FASTA jest pusty.")

    return records


def encode_batch(sequences):
    """Zamiana sekwencji na wejście modelu: Tworzy wejście (B, D, 22) oraz maskę (B, D)."""
    # Początkowo wszystkie pozycje to PAD.
    tokens = torch.full(            #tablcia do wypełnienia o wymairach liczba białek x MAX_LENGTH
        (len(sequences), MAX_LENGTH),
        AA_TO_INDEX["PAD"],
        dtype=torch.long
    )

    for row, sequence in enumerate(sequences):
        sequence = sequence[:MAX_LENGTH]
        indices = [AA_TO_INDEX[aa] for aa in sequence]

        tokens[row, :len(indices)] = torch.tensor(indices)

    # True = aminokwas, False = dopełnienie.
    mask = tokens != AA_TO_INDEX["PAD"]

    encoded = F.one_hot(tokens, num_classes=len(VOCAB)).float()

    return encoded, mask


class DummyPredictor(nn.Module):
    """Generuje losowe współrzędne zamiast prawdziwej predykcji."""

    def forward(self, encoded, mask):
        batch_size, max_length, vocabulary_size = encoded.shape

        coordinates = torch.randn(
            batch_size, max_length, 3,
            device=encoded.device
        ) * 10.0

        # Każde białko otrzymuje macierz (L, 3), bez paddingu.
        predictions = []

        for i in range(batch_size):
            predictions.append(coordinates[i][mask[i]])

        return predictions


def write_pdb(path, sequence, coordinates):
    """Zapisuje jeden atom CA na każdy aminokwas."""
    coordinates = coordinates.detach().cpu().tolist()

    if len(sequence) != len(coordinates):
        raise ValueError("Liczba aminokwasów i współrzędnych jest różna.")

    with path.open("w", encoding="ascii") as output:
        #output.write("REMARK 900 DUMMY RANDOM COORDINATES\n")

        for i, (aa, xyz) in enumerate(zip(sequence, coordinates), start=1):
            x, y, z = xyz

            #odstępy bo PDB ma ustalone kolumny
            line = (
                f"ATOM  {i:5d}  CA  {AA3[aa]:3s} A{i:4d}    "
                f"{x:8.3f}{y:8.3f}{z:8.3f}"
                f"{1.0:6.2f}{0.0:6.2f}          "
                f"{'C':>2s}  \n"
            )
            output.write(line)

        output.write(
            f"TER   {len(sequence) + 1:5d}      "
            f"{AA3[sequence[-1]]:3s} A{len(sequence):4d}\n"
        )
        output.write("END\n")


def main():
    if not 1 <= MAX_LENGTH <= 9999 or BATCH_SIZE < 1:
        raise ValueError("MAX_LENGTH: 1–9999, BATCH_SIZE: co najmniej 1.")

    if LIMIT is not None and LIMIT < 1:
        raise ValueError("LIMIT musi być dodatni albo równy None.")

    torch.manual_seed(SEED)

    records = read_fasta(FASTA_FILE)

    if LIMIT is not None:
        records = records[:LIMIT]

    #nowy folder dla każdego uruchomienia
    OUTPUT_DIR.mkdir(parents=True, exist_ok=False)

    model = DummyPredictor().eval()

    truncated = sum(len(sequence) > MAX_LENGTH for _, sequence in records)
    print(f"Białka: {len(records)}, do obcięcia: {truncated}")

    with torch.inference_mode():        #wyłączamy obliczanie gradientów, bo nie potrzebujemy ich do predykcji
        for start in range(0, len(records), BATCH_SIZE):
            batch = records[start:start + BATCH_SIZE]
            sequences = [sequence for _, sequence in batch]

            encoded, mask = encode_batch(sequences)
            predictions = model(encoded, mask)

            for number, ((protein_id, sequence), coordinates) in enumerate(
                zip(batch, predictions), start=start + 1
            ):
                sequence = sequence[:MAX_LENGTH]        #obcięcie do zapisu

                #nazwa pliku
                safe_id = "".join(
                    character if character.isalnum() else "_"
                    for character in protein_id
                )[:80]

                path = OUTPUT_DIR / f"{number:05d}_{safe_id}.pdb"
                write_pdb(path, sequence, coordinates)

            print(f"Zapisano {start + len(batch)}/{len(records)}")

    print(f"Gotowe: {OUTPUT_DIR.resolve()}")


if __name__ == "__main__":
    main()