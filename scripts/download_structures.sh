#!/usr/bin/env bash
# Script from the instructions, adapted to the combined split (Linux/WSL; requires wget).
# On Windows use: python scripts/download_structures.py <split>
#
# Usage (run from project root): ./scripts/download_structures.sh train
#                                ./scripts/download_structures.sh test

SPLIT=$1   # train / eval / test
IDS="splits/${SPLIT}.txt"
mkdir -p "data/structures/${SPLIT}/pdb" "data/structures/${SPLIT}/af"

# PDB IDs (everything that is not AF_ / MA_)
grep -v -E '^(AF_|MA_)' "$IDS" \
    | awk -F'_' '{print "https://files.rcsb.org/download/"$1".cif"}' \
    | sort -u | wget -nc -i - -P "data/structures/${SPLIT}/pdb"

# AFDB IDs
python3 -c "
import re
for line in open('${IDS}'):
    m = re.match(r'AF_AF(.+?)(F\d+)_\d+$', line.strip())
    if m:
        print(f'https://alphafold.ebi.ac.uk/files/AF-{m.group(1)}-{m.group(2)}-model_v6.cif')
" | sort -u | wget -nc -i - -P "data/structures/${SPLIT}/af"
