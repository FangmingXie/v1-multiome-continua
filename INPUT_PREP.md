# Input prep (maintainers)

Users get `data/input/` from the input bundles (`download_inputs.sh`, see [README.md](README.md)) and
never need this file. It covers how `data/input/` is built from `data/raw/` (not public), and how new
bundles are published.

## 0_prep: `data/raw/` -> `data/input/`
```
conda activate v1-multiome-continua
bash run_input_prep.sh         # runs every scripts/0_prep/NN.*.py; one log per script in results/logs/
```
0_prep reads only `data/raw/` and writes its outputs to `data/input/`. The name of each output starts
with the number of the script that writes it. `run_input_prep.sh` takes well under a minute once
`data/raw/` is cached. The other files in `data/input/` are copied in as they are; sources are in
[DATA_ATTRIBUTION.md](DATA_ATTRIBUTION.md).

| Script | What it does |
|---|---|
| `0_prep/01.extract_l23.py` | L2/3 neurons with raw counts: Yoo25 of every age (NR P6-P21, DR P12-P21) and Cheng22 P28 NR; every later step reads these extracts |
| `0_prep/02.extract_merfish_l23.py` | L2/3 MERFISH cells, P28 NR (for 1/04) and P21 NR (for 5/01); later steps read these extracts |
| `0_prep/03.yoo25_l23_regulon_targets.py` | Yoo25 L2/3 SCENIC+ regulon -> target table |
| `0_prep/04.human_wang25_regulon_targets.py` | Human (Wang25) regulon -> target table |

### `data/raw/`: sources of 0_prep

| File in `data/raw/` | Content | Used by | Source |
|---|---|---|---|
| `superdupermegaRNA_hasraw.h5ad` | Cheng22 + Yoo25 mouse V1 snRNA-seq, raw counts | 0_prep/01 | TBD |
| `P28NRDR_v1glut_rna_merfish_240729.h5ad` | Mouse P28 NR/DR V1 MERFISH; only the NR sections' L2/3 labels are used | 0_prep/02 | TBD |
| `P28NR{a,b}_{ant,pos}_l2_v1_240723.h5ad` | Mouse P28 NR MERFISH sections (500 genes) | 0_prep/02 | TBD |
| `P8P14P21NR_v1l23glut_rna_merfish_250411.h5ad` | Mouse P8/P14/P21 NR V1 L2/3 MERFISH; only the P21 NR sections are used | 0_prep/02 | TBD |
| `regulons/yoo25_L2_3/L2_3_{direct,extended}_e_regulon_metadata.csv` | Yoo25 L2/3 SCENIC+ eRegulon metadata | 0_prep/03 | TBD |
| `regulons/human_regulons_SuppTable13.xlsx` | Wang25 human eRegulons (Supp. Table 13) | 0_prep/04 | TBD |

## Publishing new input bundles
1. `bash make_input_bundle.sh data/bundle`. It packs `data/input/` into the core and external bundles,
   resolving symlinks. Each archive holds a `MANIFEST_<bundle>.tsv` (path, size, sha256 of each file) and
   `ATTRIBUTION.md` (a copy of DATA_ATTRIBUTION.md), and gets a `.sha256` file next to it. The script's
   CORE and EXTERNAL lists must cover every entry of `data/input/`, or it stops.
2. Upload both archives to Google Drive, shared as "anyone with the link".
3. In `download_inputs.sh`, update `VERSION`, the two Drive IDs and the two sha256 values. Update the
   Drive links in the README's Data table.
