# v1-multiome-continua

Code to reproduce the L2/3 archetype analyses of mouse primary visual cortex (V1). L2/3 neurons
vary along continua between archetypes A', B' and C'. The archetype genes are rebuilt from raw
counts (`0_prep/02`, identical to Table S4). Five results build on them:

| # | Result | Folder |
|---|---|---|
| 1 | Mouse vs human L2/3 | `scripts/1_mouse_human` |
| 2 | Development of mouse L2/3 under normal rearing (NR, P6-P21) | `scripts/2_development` |
| 3 | Regulon enrichment of the L2/3 archetype genes | `scripts/3_regulon` |
| 4 | Effect of dark rearing (DR, P12-P21) | `scripts/4_dr` |
| 5 | EASI-FISH projection labeling along L2/3 depth | `scripts/5_easifish` |

## Quick start
```
conda env create -f environment.yml
conda activate v1-multiome-continua
# put the inputs in data/raw/ (see "Inputs")
bash run_all.sh
```
`run_all.sh` runs every script in dependency order: 0_prep, 3_regulon, 1_mouse_human,
2_development, 4_dr, 5_easifish. It writes tables to `results/<folder>/`, figures to
`figures/<folder>/` and one log per script to `results/logs/`. The whole run takes about 10 minutes.
Each output file name starts with the number of the script that writes it.

## Scripts and figures

### 0. Shared prep
| Script | What it does |
|---|---|
| `0_prep/01.extract_yoo25_l23_p21.py` | Yoo25 P21 L2/3 neurons, NR and DR, raw counts |
| `0_prep/02.archetype_genes.py` | Per-cell A/B/C scores, archetype cells and archetype genes per condition (= Table S4) |

### 1. Mouse vs human L2/3
| Script | What it does | Figures |
|---|---|---|
| `01`, `02` | Mouse (Cheng22 P28 NR) and human (Jorstad23) L2/3: PCA + varimax | |
| `03` | Mouse L2/3 archetypes, markers and scores | |
| `04` | Mouse archetype scores on human cells | |
| `05` | Cross-species CCA of the varimax axes | |
| `06` | A'/B'/C' scores along CCA1/CCA2, human vs mouse | `06.*normalized*.pdf` |
| `07`-`10` | The same varimax and CCA at top-4000 HVGs, loadings extended to all genes | |
| `11` | Conserved A vs C genes in situ: mouse P28 MERFISH vs human MERFISH (Fang22) | `11.*.pdf` |
| `12`, `13` | Human (Wang25) AP-1/EGR regulon targets over development | `13.*.pdf` |
| `14` | B' (AP-1/EGR) regulon targets over development, mouse vs human | `14.*.pdf` |

### 2. Development (NR)
| Script | What it does | Figures |
|---|---|---|
| `01` | Number of archetype genes per age, NR vs DR; checks the exact match to Table S4 | `01.*.pdf` |
| `02` | Archetype cell groups and pseudobulk counts, P6-P21 | |
| `03` | Pseudobulk DESeq2: DR effect and its change with age (P12-P21) | |
| `04`, `05` | Archetype gene modules M1-M5 over age, NR vs DR | `05.*.pdf` |
| `06`, `07` | Archetype chromatin (ATAC) region modules over age, NR vs DR | `07.*.pdf` |

### 3. Regulon enrichment
| Script | What it does | Figures |
|---|---|---|
| `01` | Yoo25 L2/3 SCENIC+ regulon -> target table | |
| `02`, `03` | Adult L2/3 archetypes (Cheng22 + Yoo25, Harmony) and markers; sets the gene universe | |
| `04` | Expression-stratified regulon x archetype enrichment | |
| `05`, `06` | Dot plot of enriched regulons; archetype-gene overlap bars | `05.*.pdf`, `06.*.pdf` |
| `07` | TF -> target network | `07.*.pdf` |

### 4. Dark rearing
Results 2's `03`, `05` and `07` also show NR vs DR.

| Script | What it does | Figures |
|---|---|---|
| `01` | Selected genes over age, NR vs DR, per archetype cell group | `01.*.pdf` |
| `02` | Target activity of archetype-enriched regulons, NR vs DR | `02.*.pdf` |

### 5. EASI-FISH
| Script | What it does | Figures |
|---|---|---|
| `01` | LM vs RL retrograde-label depth profiles (3 samples) | `01.*.pdf` |
| `02` | Per-cell A'-C' axis, Yoo25 P21 NR | |
| `03` | Archetype scores along the A'-C' axis (snRNA) and depth (P21 MERFISH) next to the EASI-FISH profiles | `03.*.pdf` |

## Inputs
All inputs go in `data/raw/`, which is not tracked:

| File in `data/raw/` | Content | Used by | Source |
|---|---|---|---|
| `superdupermegaRNA_hasraw.h5ad` | Cheng22 + Yoo25 mouse V1 snRNA-seq, raw counts | 0_prep, 2/02, 1/14, 4/02 | TBD |
| `superdupermegaRNA_cheng22_IT_P28NR.h5ad` | Cheng22 P28 NR IT neurons | 1/01, 1/03, 1/08, 1/09, 3/02-04 | TBD |
| `superdupermegaRNA_yoo25_IT_P21.h5ad` | Yoo25 P21 IT neurons | 3/02-04 | TBD |
| `jorstad23_human_WithinArea_L23IT.h5ad` | Human L2/3 IT (Jorstad et al. 2023) | 1/02, 1/04, 1/07, 1/09 | TBD |
| `wang25_human_multiome_gex_v1_l23.h5ad` | Human developing V1 L2/3, multiome GEX (Wang et al. 2025) | 1/13 | TBD |
| `regulons/yoo25_L2_3/` | Yoo25 L2/3 SCENIC+ eRegulon metadata | 3/01 | TBD |
| `regulons/human_regulons_SuppTable13.xlsx` | Wang25 human eRegulons (Supp. Table 13) | 1/12 | TBD |
| `atac/pmat_l23concensus_{age}.h5ad` | Yoo25 L2/3 ATAC peak matrices, 11 conditions | 2/06 | TBD |
| `P28NRDR_v1glut_rna_merfish_240729.h5ad` | Mouse P28 NR/DR V1 MERFISH | 1/11 | TBD |
| `P28NR{a,b}_{ant,pos}_l2_v1_240723.h5ad` | Mouse P28 NR MERFISH sections | 1/11 | TBD |
| `fang22_human_merfish/*.csv` | Human cortex MERFISH, 9 sections (Fang et al. 2022) | 1/11 | TBD |
| `P8P14P21NR_v1l23glut_rna_merfish_250411.h5ad` | Mouse P8/P14/P21 NR V1 L2/3 MERFISH | 5/03 | TBD |
| `easifish/easifish_sample{1,2,3}_res.csv` | EASI-FISH retrograde labels (LM, RL), 3 samples | 5/01, 5/03 | TBD |
| `Cheng22_Xie25_L23_ABC_genes_n286.csv` | Adult L2/3 archetype genes, the seed of the archetype genes | 0_prep/02 | TBD |
| `gene_sets_l23{a,b,c}_250423.csv` | Module (M1-M5) of each archetype gene | 2/04 | TBD |
| `atac_modules/L23{A,B,C}_regions_M{1-5}.bed` | Archetype chromatin region modules | 2/06 | TBD |
| `human_mouse_orthologs.tsv` | Human-mouse ortholog pairs | 1/04, 1/05, 1/09-11 | TBD |
| `TableS4_L23_archetypal_genes.tsv` | Published Table S4, used only for the exact-match check | 2/01 | TBD |

## Environment and numerical notes
- **numpy:** the pipeline runs on numpy 1.26.4 (`environment.yml`). The archetype fit (`py_pcha`)
  uses `np.mat`, which numpy 2 removed. `run_all.sh` sets `PYTHONNOUSERSITE=1`, so that packages
  in `~/.local` can't shadow the env's.
- **OpenBLAS:** `environment.yml` pins OpenBLAS 0.3.32. With 0.3.34, the PCA/varimax coordinates
  shift by about 1e-4, a few tied top-cell ranks swap, and DESeq2 p-values change. With the pinned
  env, two runs of `run_all.sh` give identical tables. The one exception is the float32 marker
  statistics of `3_regulon/03`, which vary by about 1e-6.
- **DESeq2 p-values:** the DESeq2 p-values of `2_development/03` depend slightly on the numpy version.
  Under numpy 2.2.6, padj differs by up to 0.13, and about 0.5% of significance calls flip. Fold changes
  agree to 4 decimals, and the genes and stars shown in `4_dr/01` don't change.

## Layout
- `data/raw/`: inputs (not tracked).
- `src/l23/`: shared code:
  - archetype analysis (`archetype.py`, `utils.py`, `hvg.py`);
  - cross-species varimax/CCA (`crossspecies.py`);
  - regulon enrichment (`regulon.py`, `regulon_plot.py`).
- `scripts/<k>_<result>/NN.*.py`: analysis scripts, numbered in run order.
- `plan/`: the release plan.
- `results/`, `figures/`: outputs (not tracked).
