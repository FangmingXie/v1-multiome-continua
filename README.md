# v1-multiome-continua

Code to reproduce the L2/3 archetype analyses of mouse primary visual cortex (V1). L2/3 neurons
vary along continua between archetypes A', B' and C'. The archetype genes are rebuilt from raw
counts (`2_development/01`, identical to Table S4) and are also used by results 3-5. Five results:

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
pip install --no-deps -e .     # shared code src/l23/, importable by every script
# put the inputs in data/raw/ and data/input/ (see "Inputs")
bash run_input_prep.sh         # 0_prep: data/raw/ -> data/input/; skip if data/input/ is complete
bash run_all.sh                # the analysis, folders 1-5
```
`run_input_prep.sh` runs 0_prep, which writes its outputs to `data/input/`. `run_all.sh` runs the
analysis scripts in number order, folder by folder (1_mouse_human to 5_easifish); each folder reads
only from itself, lower-numbered folders and `data/input/`. They write tables to `results/<folder>/`,
figures to `figures/<folder>/`, and both runners write one log per script to `results/logs/`.
`run_all.sh` takes about 4 minutes; `run_input_prep.sh` well under one once `data/raw/` is cached.
Each output file name starts with the number of the script that writes it. Within a folder, the
computing scripts come first and write tables; the figure scripts (`*figure*.py`, and those listed
with figures below) come after them and only read tables. Numerical notes on the pinned
environment (numpy, OpenBLAS, DESeq2) are in [REPRODUCIBILITY.md](REPRODUCIBILITY.md).

## Scripts and figures

### 0. Input prep
| Script | What it does |
|---|---|
| `0_prep/01.extract_l23.py` | L2/3 neurons with raw counts: Yoo25 of every age (NR P6-P21, DR P12-P21) and Cheng22 P28 NR; every later step reads these extracts |
| `0_prep/02.extract_merfish_l23.py` | L2/3 MERFISH cells, P28 NR (for 1/04) and P21 NR (for 5/01); later steps read these extracts |
| `0_prep/03.yoo25_l23_regulon_targets.py` | Yoo25 L2/3 SCENIC+ regulon -> target table |
| `0_prep/04.human_wang25_regulon_targets.py` | Human (Wang25) regulon -> target table |

### 1. Mouse vs human L2/3
| Script | What it does | Figures |
|---|---|---|
| `01` | Mouse (Cheng22 P28 NR) and human (Jorstad23) L2/3: PCA + varimax | |
| `02` | Mouse L2/3 archetypes, markers and scores; the mouse archetype scores on human cells | |
| `03` | Cross-species CCA of the varimax axes; A'/B'/C' scores along CCA1/CCA2, human vs mouse: binned profiles and human-mouse r; per-cell CCA coordinates | |
| `04` | Conserved A vs C genes in situ, mouse P28 MERFISH vs human MERFISH (Fang22): relative L2/3 depth, profiles | |
| `05` | Figures of 03: scores along CCA1/CCA2, mouse vs human; mouse archetype scores (B', C' − A') of human and mouse cells on the CCA1 x CCA2 embedding, with the mouse archetype triangle | `05.*.pdf` |
| `06` | Figures of 04 | `06.*.pdf` |

### 2. Development (NR)
| Script | What it does | Figures |
|---|---|---|
| `01` | Per-cell A/B/C scores, archetype cells and archetype genes per age (NR and DR); checks the exact match to Table S4 | |
| `02` | Cells in P8 / P21 PC space of the global HVGs; variance per age along the P21 archetype A-C and B axes | |
| `03` | Archetype cell groups and pseudobulk counts, P6-P21 | |
| `04` | Pseudobulk DESeq2: DR effect and its change with age (P12-P21) | |
| `05` | Archetype gene modules M1-M5 over age, NR vs DR | |
| `06` | Archetype chromatin (ATAC) region modules over age, NR vs DR | |
| `07` | Figures of 01 and 02: cells in P21 PC space by archetype score, variance along the A-C and B axes, number of archetype genes per age (NR vs DR) | `07.*.pdf` |
| `08` | Figures of 05 (gene modules) and 06 (region modules) | `08.*.pdf` |

### 3. Regulon enrichment
| Script | What it does | Figures |
|---|---|---|
| `01` | Expression-stratified regulon x archetype enrichment over the adult L2/3 gene universe (Cheng22 + Yoo25); enriched TFs | |
| `02` | TF -> target network: edges and layout | |
| `03` | AP-1/EGR regulon target activity over development, human (Wang25) and mouse (B' TFs) | |
| `04`, `05` | Dot plot of enriched regulons; archetype-gene overlap bars | `04.*.pdf`, `05.*.pdf` |
| `06` | Figure of 02 | `06.regulon_network_free.pdf` |
| `07` | Figure of 03, mouse vs human | `07.*.pdf` |

### 4. Dark rearing
Results 2's `04`, `05` and `06` also show NR vs DR.

| Script | What it does | Figures |
|---|---|---|
| `01` | Target activity of archetype-enriched regulons, NR vs DR | |
| `02` | Selected genes over age, NR vs DR, per archetype cell group | `02.*.pdf` |
| `03` | Figures of 01 | `03.*.pdf` |

### 5. EASI-FISH
| Script | What it does | Figures |
|---|---|---|
| `01` | Archetype scores along the Yoo25 P21 NR A'-C' axis (snRNA) and depth (P21 MERFISH); EASI-FISH LM vs RL, 3 samples: depth profiles, x-y projection bias, per x-y bin mean-depth difference vs a label-shuffle null | |
| `02` | Figure of 01: score panels next to the EASI-FISH depth profiles | `02.*.pdf` |
| `03` | Figures of 01: EASI-FISH cell maps, depth profiles, x-y maps, null histograms | `03.*.pdf` |

## Inputs
Inputs are split by the stage that reads them; neither folder is tracked.

### `data/raw/`: sources of 0_prep
Needed only to rerun 0_prep.

| File in `data/raw/` | Content | Used by | Source |
|---|---|---|---|
| `superdupermegaRNA_hasraw.h5ad` | Cheng22 + Yoo25 mouse V1 snRNA-seq, raw counts | 0_prep/01 | TBD |
| `P28NRDR_v1glut_rna_merfish_240729.h5ad` | Mouse P28 NR/DR V1 MERFISH; only the NR sections' L2/3 labels are used | 0_prep/02 | TBD |
| `P28NR{a,b}_{ant,pos}_l2_v1_240723.h5ad` | Mouse P28 NR MERFISH sections (500 genes) | 0_prep/02 | TBD |
| `P8P14P21NR_v1l23glut_rna_merfish_250411.h5ad` | Mouse P8/P14/P21 NR V1 L2/3 MERFISH; only the P21 NR sections are used | 0_prep/02 | TBD |
| `regulons/yoo25_L2_3/L2_3_{direct,extended}_e_regulon_metadata.csv` | Yoo25 L2/3 SCENIC+ eRegulon metadata | 0_prep/03 | TBD |
| `regulons/human_regulons_SuppTable13.xlsx` | Wang25 human eRegulons (Supp. Table 13) | 0_prep/04 | TBD |

### `data/input/`: everything scripts 1-5 read
0_prep writes its outputs here (Source = the 0_prep script); the other files are external. A file
with a number prefix was written by 0_prep, one without is external. With `data/input/` complete
(0_prep's outputs included), every figure can be reproduced without `data/raw/`: run the folders
after 0_prep.

Datasets:

| File in `data/input/` | Content | Used by | Source |
|---|---|---|---|
| `01.yoo25_l23.h5ad` | Yoo25 mouse V1 L2/3 snRNA-seq, every age (NR P6-P21, DR P12-P21), raw counts | 2/01, 2/02, 2/03, 3/01, 3/03, 4/01, 5/01 | 0_prep/01 |
| `01.cheng22_p28nr_l23.h5ad` | Cheng22 mouse V1 L2/3 snRNA-seq, P28 NR, raw counts | 1/01, 1/02, 3/01 | 0_prep/01 |
| `02.merfish_p28nr_l23.h5ad` | Mouse P28 NR V1 MERFISH, L2/3 cells (500 genes) | 1/04 | 0_prep/02 |
| `02.merfish_p21nr_l23.h5ad` | Mouse P21 NR V1 MERFISH, L2/3 cells (361 genes) | 5/01 | 0_prep/02 |
| `jorstad23_human_WithinArea_L23IT.h5ad` | Human L2/3 IT (Jorstad et al. 2023) | 1/01, 1/02 | TBD |
| `fang22_human_merfish/*.csv` | Human cortex MERFISH, 6 sections (Fang et al. 2022) | 1/04 | TBD |
| `atac/pmat_l23concensus_{age}.h5ad` | Yoo25 L2/3 ATAC peak matrices, 11 conditions | 2/06 | TBD |
| `wang25_human_multiome_gex_v1_l23.h5ad` | Human developing V1 L2/3, multiome GEX (Wang et al. 2025) | 3/03 | TBD |
| `easifish/easifish_sample{1,2,3}_res.csv` | EASI-FISH retrograde labels (LM, RL), 3 samples | 5/01, 5/03 | TBD |

Small tables and annotations:

| File in `data/input/` | Content | Used by | Source |
|---|---|---|---|
| `03.yoo25_l23_regulon_targets.tsv` | Yoo25 L2/3 SCENIC+ regulon -> target table | 3/01, 3/02, 3/03, 4/01 | 0_prep/03 |
| `04.human_wang25_regulon_targets.tsv` | Wang25 human regulon -> target table | 3/03 | 0_prep/04 |
| `human_mouse_orthologs.tsv` | Human-mouse ortholog pairs | 1/02, 1/03, 1/04 | TBD |
| `Cheng22_Xie25_L23_ABC_genes_n286.csv` | Adult L2/3 archetype genes, the seed of the archetype genes | 2/01 | TBD |
| `TableS4_L23_archetypal_genes.tsv` | Published Table S4; exact-match check only | 2/01 | TBD |
| `gene_sets_l23{a,b,c}_250423.csv` | Module (M1-M5) of each archetype gene | 2/05 | TBD |
| `atac_modules/L23{A,B,C}_regions_M{1-5}.bed` | Archetype chromatin region modules | 2/06 | TBD |

## Layout
```
v1-multiome-continua/
├── README.md
├── REPRODUCIBILITY.md      numerical notes on the pinned environment
├── run_input_prep.sh       runs 0_prep (data/raw/ -> data/input/)
├── run_all.sh              runs the analysis, folders 1-5, in dependency order
├── environment.yml         pinned conda env (see REPRODUCIBILITY.md)
├── data/                   not tracked; see "Inputs"
│   ├── raw/                sources of 0_prep
│   └── input/              inputs of scripts 1-5, incl. 0_prep's outputs
├── src/l23/                shared code, imported by the scripts
│   ├── constants.py        Yoo25 ages, NR vs DR test days, colours
│   ├── expr.py             log2(CP10k+1), percentile-scaled gene-set scores
│   ├── timecourse.py       NR vs DR OLS test and time-course panel
│   ├── archetype.py        archetype analysis (PCA + PCHA)
│   ├── utils.py            helpers of archetype.py: normalization, projection, PCHA
│   ├── crossspecies.py     cross-species varimax / CCA
│   ├── regulon.py          regulon enrichment, regulon tables, target activity
│   └── regulon_plot.py     shared encoding of the enrichment figures
├── scripts/
│   ├── 0_prep/             NN.*.py, numbered in run order
│   ├── 1_mouse_human/
│   ├── 2_development/
│   ├── 3_regulon/
│   ├── 4_dr/
│   └── 5_easifish/
├── results/<folder>/       tables (not tracked); results/logs/ = one log per script
├── figures/<folder>/       figures (tracked)
└── plan/                   development plans (not tracked)
```
