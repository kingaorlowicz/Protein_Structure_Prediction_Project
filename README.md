# PSP – Protein Structure Prediction

**Cel projektu:** framework end-to-end, który przewiduje strukturę 3D białka (plik PDB/mmCIF) na podstawie sekwencji aminokwasów.

**Meeting 2:** przygotowanie zbioru treningowego, walidacyjnego i testowego, **bez powtórzonych plików i bez wycieku danych** między splitami.

**Założenie:** jeden wspólny zbiór. **Train i eval mieszają struktury eksperymentalne (PDB) i modele AlphaFold (AFDB), a test zawiera tylko PDB.** Model ocenia się więc wyłącznie na prawdziwych, eksperymentalnych strukturach.

## Struktura projektu

```
2wersja/
├── README.md
├── requirements.txt
├── .gitignore
├── .vscode/                      # zadania (Ctrl+Shift+B) i konfiguracje uruchamiania (F5)
├── docs/                         # temat, instrukcja podziału, opis projektów (PDF)
├── scripts/
│   ├── common.py                 # wczytywanie klastrów (md5/metadane), parsowanie ID PDB/AFDB
│   ├── prepare_splits.py         # klastry 50% -> train/eval/test
│   ├── check_leakage.py          # audyt: skład, duplikaty w splitach, wycieki między nimi
│   ├── download_structures.py    # pobieranie .cif (Python, działa na Windows)
│   ├── download_structures.sh    # to samo w bash + wget (Linux/WSL)
│   ├── extract_chain.py          # ID z PDB: wycięcie jednego łańcucha encji
│   ├── make_fasta.py             # FASTA (wejście modelu) z pełną sekwencją encji
│   └── run_pipeline.py           # uruchamia powyższe kroki po kolei
├── data/
│   ├── raw/                      # clusters-by-entity-50.txt (dostarczony plik) + metadane pobrania
│   ├── structures/<split>/pdb/   # pobrane wpisy PDB (.cif)
│   ├── structures/<split>/af/    # pobrane modele AlphaFold (.cif)
│   ├── chains/<split>/           # pojedyncze łańcuchy wycięte z wpisów PDB
│   └── fasta/<split>.fasta       # sekwencje wejściowe
├── splits/
│   ├── train.txt / eval.txt / test.txt   # ID białek, jedno na linię
│   ├── metadata.json             # seed, źródła, skład, md5 pliku, klastry per split
│   └── leakage_report.txt        # wynik audytu
├── src/psp/                      # kod modelu (kolejne etapy)
├── notebooks/                    # eksploracja / wykresy
├── models/                       # checkpointy (poza gitem)
└── results/                      # predykcje .pdb, metryki (poza gitem)
```

## Uruchomienie

```bash
pip install -r requirements.txt
python scripts/run_pipeline.py                        # wszystko po kolei
python scripts/run_pipeline.py --steps split check    # tylko podział + audyt
```

W VS Code: **Ctrl+Shift+B** uruchamia podział z audytem, **Terminal → Run Task…** pokazuje pozostałe kroki. Wszystkie polecenia uruchamiamy z katalogu `2wersja/`.

## 1. Podział train/eval/test (`prepare_splits.py`)

**Wejście:** tylko dostarczony plik `data/raw/clusters-by-entity-50.txt`, czyli klastry RCSB przy 50% identyczności sekwencji, jeden klaster na linię (md5 `e61dccb7a61e0f43d62590b79702c978`, wersja z 24.10.2026).

Procedura:
1. Wczytanie klastrów i odrzucenie członków `MA_` (ModelArchive).
2. Przetasowanie klastrów (`seed = 42`).
3. Wypełnianie splitów **całymi klastrami** w kolejności test → eval → train:
   - **test:** tylko członkowie z PDB (`--test-source pdb`);
   - **eval, train:** członkowie z PDB i AFDB razem (`--train-source all`; `--eval-source` domyślnie taki sam).

   Klaster trafia do pierwszego splitu, który jeszcze potrzebuje klastrów i w którym klaster ma ≥ 10 użytecznych członków. Z niego losujemy 10.
4. Użyteczny członek to taki, którego **plik źródłowy nie jest użyty nigdzie indziej w zbiorze**. Dla PDB plikiem jest wpis (`1ABC` dla `1ABC_1`/`1ABC_2`), a dla AFDB accession UniProt (fragmenty `F1`/`F2` tego samego białka się nakładają).
5. 400 / 50 / 50 klastrów, czyli **4000 / 500 / 500 białek**.
6. Zapis `splits/train.txt`, `eval.txt`, `test.txt` i `metadata.json`.

Gwarancje:
- ten sam plik `.cif` nie występuje dwa razy w splicie ani w dwóch splitach;
- żaden klaster 50% nie jest w dwóch splitach, więc między splitami nie ma pary białek o ≥ 50% identyczności. Dotyczy to także pary model AlphaFold w train i struktura PDB tego samego białka w teście.

### Dlaczego pilnujemy plików

Pierwsza wersja (`../1wersja`), która losowała po prostu 10 członków z klastra, miała w PDB **47 tych samych plików `.cif` jednocześnie w train i test** oraz 28 w train i eval. Encje `1ABC_1` i `1ABC_2` pochodzą z jednego pliku, a mogą należeć do różnych klastrów, np. łańcuch ciężki i lekki przeciwciała. Model widziałby wtedy w treningu fragment struktury testowej, co zawyża wynik i sprzyja przeuczeniu.

### Wynik (2026-10-30, seed 42)

| split | białka | klastry | PDB | AFDB |
|-------|-------:|--------:|----:|-----:|
| train | 4000   | 400     | 1394 | 2606 |
| eval  | 500    | 50      | 53   | 447  |
| test  | 500    | 50      | 500  | 0    |

Audyt (`splits/leakage_report.txt`): między każdą parą splitów wspólnych ID = 0, wspólnych plików = 0, wspólnych klastrów 50% = 0. Wewnątrz splitów powtórzonych plików = 0.

### Na co uważać

- **Eval ma inny rozkład niż test** (89% AFDB, a test to 100% PDB). Jeśli eval ma służyć do wyboru modelu pod test, lepiej użyć `python scripts/prepare_splits.py --eval-source pdb`, co daje eval 500 PDB i train 1099 PDB + 2901 AFDB.
- **Test wybierany jest z klastrów, które mają ≥ 10 różnych wpisów PDB**, czyli z dobrze zbadanych rodzin białek.
- **Identyczne sekwencje wewnątrz splitu.** Sam plik 50% nie zawiera sekwencji, więc nie odróżni różnych wpisów PDB tego samego białka ani modelu AF i struktury PDB tego samego białka. Kontrolny pomiar z klastrami RCSB 100% (poza pipeline’em) pokazał takie powtórzenia: train 1124 / 4000, eval 107 / 500, test 339 / 500. Między splitami jest ich 0, więc nie ma wycieku, ale część białek ma większą wagę w treningu. Sprawdzenie na prawdziwych sekwencjach: krok 4.
- Modele typu AlphaFold/ESMFold widziały w treningu większość PDB. Jeśli użyjecie ich jako baseline’u albo do fine-tuningu, wynik na teście będzie optymistyczny.

## 2. Pobranie struktur

```bash
python scripts/download_structures.py train      # oraz eval, test
```

Przekształcanie ID (wzory z `download_structures.sh` w instrukcji, `docs/split_instructions.txt`):
- PDB: `1ABC_2` → `https://files.rcsb.org/download/1ABC.cif` → `data/structures/<split>/pdb/`
- AFDB: `AF_AFP12345F1_1` → `https://alphafold.ebi.ac.uk/files/AF-P12345-F1-model_v6.cif` → `data/structures/<split>/af/`

Pliki już pobrane są pomijane. Nieudane pobrania zapisują się do `data/structures/<split>/failed.tsv`.

## 3. Wycięcie łańcucha z wpisów PDB

```bash
python scripts/extract_chain.py - data/structures/train/pdb data/chains/train --ids-file splits/train.txt --workers 8
```

Analogicznie dla `eval` i `test`. ID z AFDB są pomijane, bo model AlphaFold to już jeden łańcuch.

Poprawki względem fragmentu z instrukcji:
- Mapowanie encja → łańcuch idzie przez `_atom_site.label_entity_id` → `_atom_site.auth_asym_id`, bo Biopython używa **auth** chain ID, a `_struct_asym.id` to **label** ID.
- Zapisywany jest **jeden łańcuch** encji (homo-oligomer dałby kilka identycznych kopii) i tylko pierwszy model (NMR). Opcja `--all-chains` zachowuje wszystkie kopie.
- Lista ID jest wczytywana z pliku (`--ids-file`, działa w PowerShellu) albo ze stdin i przetwarzana równolegle.

## 4. Sekwencje FASTA i końcowa kontrola

```bash
python scripts/make_fasta.py train               # -> data/fasta/train.fasta (oraz eval, test)
python scripts/check_leakage.py --fasta-dir data/fasta
```

Sekwencja to pełna kanoniczna sekwencja encji (`_entity_poly.pdbx_seq_one_letter_code_can`), a nie tylko reszty widoczne w krysztale. Nagłówek ma postać `>ID source=pdb|af len=N`. Audyt z `--fasta-dir` porównuje rzeczywiste sekwencje, w splitach i między nimi.

## Reprodukowalność

- Stały seed (`--seed`, domyślnie 42). Ten sam seed i ten sam plik klastrów dają identyczny podział.
- Dla pliku klastrów zapisane są md5, rozmiar, URL, data dostępu, `Last-Modified` i ETag (`data/raw/clusters-by-entity-50.download.json` oraz `metadata.json`).
- Wersja modeli AFDB ustalona na `v6`.
- `metadata.json` zawiera parametry, skład i pełną listę wylosowanych członków każdego klastra w każdym splicie.
- Plik klastrów RCSB aktualizuje się co tydzień, więc do odtworzenia podziału potrzebny jest dokładnie dostarczony plik (sprawdźcie md5).
