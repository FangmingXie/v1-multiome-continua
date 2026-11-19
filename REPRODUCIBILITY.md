# Reproducibility notes

Numerical details behind the pinned environment (`environment.yml`) and how far results move
outside it.

- **numpy:** the pipeline runs on numpy 1.26.4 (`environment.yml`). The archetype fit (`py_pcha`)
  uses `np.mat`, which numpy 2 removed. `run_all.sh` sets `PYTHONNOUSERSITE=1`, so that packages
  in `~/.local` can't shadow the env's.
- **OpenBLAS:** `environment.yml` pins OpenBLAS 0.3.32. With 0.3.34, the PCA/varimax coordinates
  shift by about 1e-4, a few tied top-cell ranks swap, and DESeq2 p-values change. With the pinned
  env, two runs of `run_all.sh` give identical tables.
- **DESeq2 p-values:** the DESeq2 p-values of `2_development/03` depend slightly on the numpy version.
  Under numpy 2.2.6, padj differs by up to 0.13, and about 0.5% of significance calls flip. Fold changes
  agree to 4 decimals, and the genes and stars shown in `4_dr/02` don't change.
