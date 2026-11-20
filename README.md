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
bash download_inputs.sh        # data/input/ from the input bundles (see "Data")
bash run_all.sh                # the analysis, folders 1-5
```
`download_inputs.sh` downloads the input bundles from Google Drive with `curl` (about 1.9 GB), checks their
sha256 and unpacks them into `data/input/`. `run_all.sh` runs the analysis scripts in number order,
folder by folder (1_mouse_human to 5_easifish); each folder reads only from itself, lower-numbered
folders and `data/input/`. They write tables to `results/<folder>/`,
figures to `figures/<folder>/`, and one log per script to `results/logs/`. It takes about 4 minutes.
Each output file name starts with the number of the script that writes it. Within a folder, the
computing scripts come first and write tables; the figure scripts (`*figure*.py`, and those listed
with figures below) come after them and only read tables. Numerical notes on the pinned
environment (numpy, OpenBLAS, DESeq2) are in [REPRODUCIBILITY.md](REPRODUCIBILITY.md).

## Data
`download_inputs.sh` installs `data/input/` from two archives (version 2026-11-20) on Google Drive:

| Bundle | Content | Size (packed) | Needed by |
|---|---|---|---|
| [core](https://drive.google.com/file/d/10kRY471qCVIJj3j6GRAzAmrYdjXbCWNS/view?usp=sharing) | every input except the three external human datasets | 1.7 GB (741 MB) | folders 2-5 |
| [external](https://drive.google.com/file/d/1yhkWYA0xS0GC7xOS-CQS1K3ZN45ALzog/view?usp=sharing) | Jorstad23, Wang25 and Fang22 human datasets | 4 GB (1.1 GB) | folder 1, 3_regulon/03 and its figure 07 |

The script checks each archive against the sha256 recorded in `download_inputs.sh`, and that every file
listed in the archive's `MANIFEST_<bundle>.tsv` (path, size and sha256 of each file) is unpacked with its
size. `--core-only` skips the external bundle; `--from <dir>` installs archives already downloaded (e.g.
from the links above in a browser).

Each file's content, the scripts that read it, its source and licence: [DATA_ATTRIBUTION.md](DATA_ATTRIBUTION.md).
Rebuilding `data/input/` from `data/raw/` (0_prep) and publishing new bundles (maintainers):
[INPUT_PREP.md](INPUT_PREP.md).

## Scripts and figures

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

## Layout
```
v1-multiome-continua/
├── README.md
├── REPRODUCIBILITY.md      numerical notes on the pinned environment
├── DATA_ATTRIBUTION.md     every input file: content, scripts that read it, source, licence
├── INPUT_PREP.md           (maintainers) 0_prep, data/raw/, publishing bundles
├── download_inputs.sh      installs data/input/ from the input bundles
├── make_input_bundle.sh    (maintainers) packs data/input/ into the input bundles
├── run_input_prep.sh       (maintainers) runs 0_prep (data/raw/ -> data/input/)
├── run_all.sh              runs the analysis, folders 1-5, in dependency order
├── environment.yml         pinned conda env (see REPRODUCIBILITY.md)
├── data/                   not tracked
│   ├── input/              inputs of scripts 1-5 (DATA_ATTRIBUTION.md)
│   └── raw/                (maintainers) sources of 0_prep (INPUT_PREP.md)
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
│   ├── 0_prep/             (maintainers) data/raw/ -> data/input/ (INPUT_PREP.md)
│   ├── 1_mouse_human/
│   ├── 2_development/
│   ├── 3_regulon/
│   ├── 4_dr/
│   └── 5_easifish/
├── results/<folder>/       tables (not tracked); results/logs/ = one log per script
├── figures/<folder>/       figures (tracked)
└── plan/                   development plans (not tracked)
```
