#!/usr/bin/env bash
# Run the analysis (folders 1-5, reading data/input/) in dependency order, from the repo root, with the
# active Python env. Given folders as arguments, runs those instead (run_input_prep.sh passes 0_prep).
# Usage: conda activate v1-multiome-continua && bash run_all.sh
# One log per script in results/logs/. Stops at the first failing script.
set -euo pipefail

cd "$(dirname "$0")"
export PYTHONNOUSERSITE=1          # ignore ~/.local packages; use only the env's (numpy 1.26.4)
export PYTHONPATH=src

# Every folder reads only from itself, lower-numbered folders and data/input/.
if [ $# -gt 0 ]; then FOLDERS=("$@"); else FOLDERS=(1_mouse_human 2_development 3_regulon 4_dr 5_easifish); fi

python -c "import numpy; assert numpy.__version__ == '1.26.4', f'numpy {numpy.__version__}, expected 1.26.4'"
mkdir -p results/logs
for f in "${FOLDERS[@]}"; do
    for s in scripts/"$f"/[0-9][0-9].*.py; do
        log="results/logs/$f.$(basename "$s" .py | cut -d. -f1).log"
        echo "$(date +%T) $s"
        python -u "$s" > "$log" 2>&1 || { echo "FAILED: $s (see $log)"; tail -20 "$log"; exit 1; }
    done
done
# README.md is the figure report: every PDF in figures/ must be linked from it (skipped for 0_prep alone)
if [ "${FOLDERS[*]}" != "0_prep" ]; then python check_readme_figures.py; fi
echo "$(date +%T) done"
