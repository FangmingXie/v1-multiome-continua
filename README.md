# v1-multiome-continua

Code to reproduce the L2/3 archetype analyses of mouse V1:

1. Mouse vs human L2/3
2. Development of mouse L2/3 under normal rearing (NR, P6-P21)
3. Regulon enrichment of the L2/3 archetype genes
4. Effect of dark rearing (DR, P12-P21)
5. EASI-FISH projection labeling along L2/3 depth

Work in progress; see `plan/00.release_plan.md`.

## Layout
- `data/`: small inputs (tracked). `data/raw/`: large inputs (not tracked; download links TBD).
- `src/l23/`: shared code.
- `scripts/<k>_<result>/`: analysis scripts, numbered in run order.
- `results/`, `figures/`: outputs (not tracked).

## Environment
```
conda env create -f environment.yml
conda activate v1-multiome-continua
pip install -e .
```
