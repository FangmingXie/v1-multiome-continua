# CLAUDE.md

Repository layout, inputs and the script index are in README.md.

## Environment

- Use this conda env to run this project: `v1-multiome-continua`, that means to run any python script with `conda run -n v1-multiome-continua`
- For heavy-lifting scripts, prefer running unbuffered at both levels so output (progress, logs) streams in real time: use `conda run --no-capture-output -n v1-multiome-continua python -u <script>.py`
- Always prefix with `PYTHONNOUSERSITE=1` (as `run_all.sh` does): without it, numpy 2 from `~/.local` shadows the env's numpy 1.26.4.
- The shared code `src/l23/` is installed in the env in editable mode (`pip install --no-deps -e .`), so scripts import `l23` from anywhere and edits to `src/l23/` take effect without reinstalling.

## coding styles
- Define all file paths (input and output files) in the beginning of each script as much as possible. Capitalize the variables that store these file paths.

**Simplify Relentlessly**: Remove complexity aggressively - the simplest design that works is usually best

#### Fail-Fast, No Fallbacks
- **No Silent Fallbacks**: Code must fail immediately when expected conditions aren't met. Silent fallback behavior masks bugs and creates unpredictable systems.
- **Explicit Error Messages**: When something goes wrong, stop execution with clear error messages explaining what failed and what was expected.
- **Example**: `raise ValueError(f"Required model {model_name} not found")` instead of falling back to first available model.

## Script conventions
- Scripts are `scripts/<folder>/NN.name.py`, numbered in run order within a folder; `run_all.sh` runs every `[0-9][0-9].*.py` in each folder.
- Each output file name starts with the number of the script that writes it: tables to `results/<folder>/`, figures to `figures/<folder>/`, the log to `results/logs/<folder>.NN.log`.
- Computing scripts come first and write tables; figure scripts come after them and only read tables.
- Each script starts with a docstring stating what it does, what it reads and what it writes.

## Reproducibility
- Do not change package versions in `environment.yml`; numerical notes on the pinned env (numpy, OpenBLAS, DESeq2) are in REPRODUCIBILITY.md.
- Seed every random step (shuffles, subsampling, random initializations) so reruns give the same numbers.

## Imported development scripts
- Scripts copied from scomics (e.g. `scripts/2_development/07.*`) are kept as is, with only the changes needed to run here, each marked `# CHANGED` and listed at the top of the file, until they are cleaned up.
