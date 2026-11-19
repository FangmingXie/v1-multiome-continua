"""Clean regulon -> target-gene table from Wang25 human SCENIC+ eRegulons.

Human analog of 3_regulon/01 (mouse Yoo25 L2/3); the merge is l23.regulon.merge_direct_extended.
The input here is a single Excel
supplementary table (Wang25, human cortex multiome) rather than the SCENIC+ CSVs used
for the mouse data, so the parsing differs:

  - Source: `Supplementary Table 13b` (eRegulon target GENES) inside
    human_regulons_SuppTable13.xlsx. Each row is one eRegulon with a single
    comma-separated `Target genes` cell.
  - The eRegulon name is encoded as `<TF>[_extended]_<TF2Gsign>_<R2Gsign>`
    (e.g. `ETV1_+_+` = direct, `RORB_extended_+_+` = extended, `NKX6-2_-_+`).
    There is NO explicit "direct" tag: a regulon is direct unless it carries the
    `_extended_` tag. Signs are joined with `/` to match the mouse output
    (`+_+` -> `+/+`). TF names may contain hyphens (e.g. NKX6-2) but not
    underscores, so the last two underscore tokens are always the signs.

The Excel file is read with the Python standard library (zipfile + ElementTree)
to avoid an openpyxl dependency.

Same regulon model and dedup rule as the mouse script:
  - A regulon is identified by TF + sign pattern (direct/extended normalized away),
    e.g. `RORB_+/+`.
  - If a regulon exists in both direct and extended, keep its DIRECT targets;
    otherwise keep its EXTENDED targets.
  - The comma-separated gene list is exploded to unique (regulon, Gene) pairs.

Reads:
  data/raw/regulons/human_regulons_SuppTable13.xlsx
Outputs:
  results/1_mouse_human/07.human_wang25_regulon_targets.tsv
"""

import os
import zipfile
from xml.etree import ElementTree as ET

import pandas as pd

from l23.regulon import REGULON_COLS, merge_direct_extended

REPO_ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

XLSX = os.path.join(REPO_ROOT, "data", "raw", "regulons", "human_regulons_SuppTable13.xlsx")
OUT_TSV = os.path.join(REPO_ROOT, "results", "1_mouse_human", "07.human_wang25_regulon_targets.tsv")

GENE_SHEET = "xl/worksheets/sheet2.xml"  # "Supplementary Table 13b" — target genes

NS = "{http://schemas.openxmlformats.org/spreadsheetml/2006/main}"

# keep only regulons with a positive R2G (second) sign; drop +/- and -/-
KEEP_DIRECTIONS = {"+/+", "-/+"}


def read_sheet_rows(xlsx_path, sheet_xml):
    """Return the sheet as a list of row-cell-string lists (stdlib xlsx reader)."""
    with zipfile.ZipFile(xlsx_path) as z:
        shared = [
            "".join(t.text or "" for t in si.iter(NS + "t"))
            for si in ET.fromstring(z.read("xl/sharedStrings.xml"))
        ]
        sheet = ET.fromstring(z.read(sheet_xml))

    def cell_text(cell):
        v = cell.find(NS + "v")
        if v is None:
            return ""
        return shared[int(v.text)] if cell.get("t") == "s" else v.text

    rows = []
    for row in sheet.find(NS + "sheetData"):
        rows.append([cell_text(c) for c in row])
    return rows


def parse_regulon(name):
    """Parse '<TF>[_extended]_<s1>_<s2>' -> (regulon_key, TF, direction, source)."""
    source = "extended" if "_extended_" in name else "direct"
    toks = name.replace("_extended", "").split("_")
    s1, s2 = toks[-2], toks[-1]
    if not {s1, s2} <= {"+", "-"}:
        raise ValueError(f"Unexpected sign tokens in eRegulon name: {name!r}")
    tf = "_".join(toks[:-2])
    direction = f"{s1}/{s2}"
    return f"{tf}_{direction}", tf, direction, source


def main():
    rows = read_sheet_rows(XLSX, GENE_SHEET)
    # row 0 = title, row 1 = header; data starts at row 2.
    # cols: [eRegulon, eRegulon_(number of target genes), Target genes]
    records = []
    for cells in rows[2:]:
        if not cells or not cells[0]:
            continue
        name, _count, genes = cells[0], cells[1], cells[2]
        regulon, tf, direction, source = parse_regulon(name.strip())
        # keep only regulons whose R2G (second) sign is positive
        if direction not in KEEP_DIRECTIONS:
            continue
        for gene in genes.split(","):
            gene = gene.strip()
            if gene:
                records.append((regulon, tf, direction, source, gene))

    table = merge_direct_extended(pd.DataFrame.from_records(records, columns=REGULON_COLS))

    # sanity checks
    n_direct = table.loc[table["source"] == "direct", "regulon"].nunique()
    n_extended = table.loc[table["source"] == "extended", "regulon"].nunique()
    n_pos = table.loc[table["regulation_direction"] == "+/+", "regulon"].nunique()
    n_neg = table.loc[table["regulation_direction"] == "-/+", "regulon"].nunique()

    os.makedirs(os.path.dirname(OUT_TSV), exist_ok=True)
    table.to_csv(OUT_TSV, sep="\t", index=False)

    print(f"[human] regulons total         : {table['regulon'].nunique()}")
    print(f"[human]   from direct          : {n_direct}")
    print(f"[human]   from extended (only) : {n_extended}")
    print(f"[human]   +/+ (activating)     : {n_pos}")
    print(f"[human]   -/+ (repressing)     : {n_neg}")
    print(f"[human] unique target genes    : {table['Gene'].nunique()}")
    print(f"[human] target pairs (rows)    : {len(table)}")
    print(f"[human] wrote -> {OUT_TSV}")


if __name__ == "__main__":
    main()
