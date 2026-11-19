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
`figures/<folder>/` and one log per script to `results/logs/`. The whole run takes about 3 minutes.
Each output file name starts with the number of the script that writes it. Within a folder, the
computing scripts come first and write tables; the figure scripts (`*figure*.py`, and those listed
with figures below) come after them and only read tables. Numerical notes on the pinned
environment (numpy, OpenBLAS, DESeq2) are in [REPRODUCIBILITY.md](REPRODUCIBILITY.md).

## Scripts and figures

### 0. Shared prep
| Script | What it does |
|---|---|
| `0_prep/01.extract_yoo25_l23.py` | Yoo25 L2/3 neurons of every age (NR P6-P21, DR P12-P21), raw counts; every later step reads this extract |
| `0_prep/02.archetype_genes.py` | Per-cell A/B/C scores, archetype cells and archetype genes per condition; checks the exact match to Table S4 |

### 1. Mouse vs human L2/3
| Script | What it does | Figures |
|---|---|---|
| `01` | Mouse (Cheng22 P28 NR) and human (Jorstad23) L2/3: PCA + varimax | |
| `02` | Mouse L2/3 archetypes, markers and scores | |
| `03` | Mouse archetype scores on human cells | |
| `04` | Cross-species CCA of the varimax axes | |
| `05` | A'/B'/C' scores along CCA1/CCA2, human vs mouse: binned profiles and human-mouse r | |
| `06` | Conserved A vs C genes in situ, mouse P28 MERFISH vs human MERFISH (Fang22): relative L2/3 depth, profiles | |
| `07` | Human (Wang25) regulon -> target table | |
| `08` | AP-1/EGR regulon target activity over development, human (Wang25) and mouse (B' TFs) | |
| `09` | Figure of 05 | `09.l23_abc_normalized_compact.pdf` |
| `10` | Figures of 06 | `10.*.pdf` |
| `11` | Figure of 08, mouse vs human | `11.*.pdf` |

### 2. Development (NR)
| Script | What it does | Figures |
|---|---|---|
| `01` | Archetype cell groups and pseudobulk counts, P6-P21 | |
| `02` | Pseudobulk DESeq2: DR effect and its change with age (P12-P21) | |
| `03` | Archetype gene modules M1-M5 over age, NR vs DR | |
| `04` | Archetype chromatin (ATAC) region modules over age, NR vs DR | |
| `05` | Number of archetype genes per age, NR vs DR | `05.*.pdf` |
| `06` | Figures of 03 (gene modules) and 04 (region modules) | `06.*.pdf` |

### 3. Regulon enrichment
| Script | What it does | Figures |
|---|---|---|
| `01` | Yoo25 L2/3 SCENIC+ regulon -> target table | |
| `02` | Gene universe: genes expressed in adult L2/3 (Cheng22 + Yoo25) and their expression stats | |
| `03` | Expression-stratified regulon x archetype enrichment; enriched TFs | |
| `04` | TF -> target network: edges and layout | |
| `05`, `06` | Dot plot of enriched regulons; archetype-gene overlap bars | `05.*.pdf`, `06.*.pdf` |
| `07` | Figure of 04 | `07.regulon_network_free.pdf` |

### 4. Dark rearing
Results 2's `02`, `03` and `04` also show NR vs DR.

| Script | What it does | Figures |
|---|---|---|
| `01` | Target activity of archetype-enriched regulons, NR vs DR | |
| `02` | Selected genes over age, NR vs DR, per archetype cell group | `02.*.pdf` |
| `03` | Figures of 01 | `03.*.pdf` |

### 5. EASI-FISH
| Script | What it does | Figures |
|---|---|---|
| `01` | Per-cell A'-C' axis, Yoo25 P21 NR | |
| `02` | Archetype scores along the A'-C' axis (snRNA) and depth (P21 MERFISH); EASI-FISH LM vs RL depth profiles (3 samples) | |
| `03` | Figure of 02 | `03.*.pdf` |

## Inputs
All inputs go in `data/raw/`, which is not tracked:

| File in `data/raw/` | Content | Used by | Source |
|---|---|---|---|
| `superdupermegaRNA_hasraw.h5ad` | Cheng22 + Yoo25 mouse V1 snRNA-seq, raw counts | 0_prep/01 | TBD |
| `superdupermegaRNA_cheng22_IT_P28NR.h5ad` | Cheng22 P28 NR IT neurons | 1/01, 1/02, 3/02 | TBD |
| `superdupermegaRNA_yoo25_IT_P21.h5ad` | Yoo25 P21 IT neurons | 3/02 | TBD |
| `jorstad23_human_WithinArea_L23IT.h5ad` | Human L2/3 IT (Jorstad et al. 2023) | 1/01, 1/03 | TBD |
| `wang25_human_multiome_gex_v1_l23.h5ad` | Human developing V1 L2/3, multiome GEX (Wang et al. 2025) | 1/08 | TBD |
| `regulons/yoo25_L2_3/` | Yoo25 L2/3 SCENIC+ eRegulon metadata | 3/01 | TBD |
| `regulons/human_regulons_SuppTable13.xlsx` | Wang25 human eRegulons (Supp. Table 13) | 1/07 | TBD |
| `atac/pmat_l23concensus_{age}.h5ad` | Yoo25 L2/3 ATAC peak matrices, 11 conditions | 2/04 | TBD |
| `P28NRDR_v1glut_rna_merfish_240729.h5ad` | Mouse P28 NR/DR V1 MERFISH | 1/06 | TBD |
| `P28NR{a,b}_{ant,pos}_l2_v1_240723.h5ad` | Mouse P28 NR MERFISH sections | 1/06 | TBD |
| `fang22_human_merfish/*.csv` | Human cortex MERFISH, 9 sections (Fang et al. 2022) | 1/06 | TBD |
| `P8P14P21NR_v1l23glut_rna_merfish_250411.h5ad` | Mouse P8/P14/P21 NR V1 L2/3 MERFISH | 5/02 | TBD |
| `easifish/easifish_sample{1,2,3}_res.csv` | EASI-FISH retrograde labels (LM, RL), 3 samples | 5/02 | TBD |
| `Cheng22_Xie25_L23_ABC_genes_n286.csv` | Adult L2/3 archetype genes, the seed of the archetype genes | 0_prep/02 | TBD |
| `gene_sets_l23{a,b,c}_250423.csv` | Module (M1-M5) of each archetype gene | 2/03 | TBD |
| `atac_modules/L23{A,B,C}_regions_M{1-5}.bed` | Archetype chromatin region modules | 2/04 | TBD |
| `human_mouse_orthologs.tsv` | Human-mouse ortholog pairs | 1/03, 1/04, 1/06 | TBD |
| `TableS4_L23_archetypal_genes.tsv` | Published Table S4, used only for the exact-match check | 0_prep/02 | TBD |

## Layout
```
v1-multiome-continua/
├── README.md
├── REPRODUCIBILITY.md      numerical notes on the pinned environment
├── run_all.sh              runs every script in dependency order
├── environment.yml         pinned conda env (see REPRODUCIBILITY.md)
├── data/raw/               inputs (not tracked; see "Inputs")
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
├── figures/<folder>/       figures (not tracked)
└── plan/                   development plans
```
