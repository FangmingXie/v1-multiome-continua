#!/usr/bin/env bash
# Run 0_prep: build data/input/'s generated files from data/raw/ (see README "Inputs").
# Usage: conda activate v1-multiome-continua && bash run_input_prep.sh
# One log per script in results/logs/. Stops at the first failing script.
set -euo pipefail
exec bash "$(dirname "$0")/run_all.sh" 0_prep
