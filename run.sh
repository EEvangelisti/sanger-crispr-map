#!/usr/bin/env bash
# Example run of analyze_edits.py on .ab1 chromatograms, with gel bands
# sequenced separately. Creates a virtual environment on first use.
set -euo pipefail

VENV=crispr

if [ ! -d "$VENV" ]; then
    python3 -m venv "$VENV"
fi
# shellcheck disable=SC1091
source "$VENV/bin/activate"
python -c "import Bio, matplotlib" 2>/dev/null || python -m pip install --upgrade pip biopython matplotlib

python analyze_edits.py \
    --template template.fas \
    --reads chromatograms \
    --out results \
    --ab1 \
    --split-bands \
    --pam-motif TTTV \
    --gap-open -10 \
    --gap-extend -0.01 \
    --min-flank 15 \
    --min-aligned 10 \
    --min-indel 5 \
    --font-size 15 \
    --row-height 0.38 \
    --samples 2_75 2_04 2_23 2_24 2_37 2_42 2_73 3_16 3_36 \
    --samples-1 2_75 3_16 \
    --samples-2 2_04 2_23 2_24 2_37 2_42 2_73 3_16 3_36 \
    --verbose > alignment_debug.txt 2>&1

deactivate
