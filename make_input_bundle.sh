#!/usr/bin/env bash
# Maintainers: pack data/input/ into the two input bundles that download_inputs.sh fetches.
#   core      every file the analysis reads except the external datasets (incl. 0_prep's outputs)
#   external  the external human datasets (Jorstad23, Wang25, Fang22)
# Each bundle holds its files (paths relative to data/input/, symlinks resolved), MANIFEST_<bundle>.tsv
# (path, bytes, sha256 of every file) and ATTRIBUTION.md (= DATA_ATTRIBUTION.md); next to it goes a
# .sha256 file. Fails if a listed file is missing or if data/input/ holds an entry in neither list.
# Usage: bash make_input_bundle.sh <out_dir> [version]      (version defaults to today, YYYY-MM-DD)
set -euo pipefail

REPO_ROOT="$(cd "$(dirname "$0")" && pwd)"
INPUT_DIR="$REPO_ROOT/data/input"
ATTRIBUTION="$REPO_ROOT/DATA_ATTRIBUTION.md"
OUT_DIR="${1:?usage: bash make_input_bundle.sh <out_dir> [version]}"
VERSION="${2:-$(date +%F)}"
PREFIX="v1-multiome-continua_input"

# entries of data/input/ (files or directories) per bundle
CORE=(
    01.yoo25_l23.h5ad
    01.cheng22_p28nr_l23.h5ad
    02.merfish_p28nr_l23.h5ad
    02.merfish_p21nr_l23.h5ad
    03.yoo25_l23_regulon_targets.tsv
    04.human_wang25_regulon_targets.tsv
    atac
    atac_modules
    easifish
    Cheng22_Xie25_L23_ABC_genes_n286.csv
    TableS4_L23_archetypal_genes.tsv
    gene_sets_l23a_250423.csv
    gene_sets_l23b_250423.csv
    gene_sets_l23c_250423.csv
    human_mouse_orthologs.tsv
)
EXTERNAL=(
    jorstad23_human_WithinArea_L23IT.h5ad
    wang25_human_multiome_gex_v1_l23.h5ad
    fang22_human_merfish
)

# every listed entry exists, and every entry of data/input/ is listed
for e in "${CORE[@]}" "${EXTERNAL[@]}"; do
    [ -e "$INPUT_DIR/$e" ] || { echo "ERROR: $INPUT_DIR/$e not found" >&2; exit 1; }
done
for p in "$INPUT_DIR"/*; do
    e="$(basename "$p")"
    printf '%s\n' "${CORE[@]}" "${EXTERNAL[@]}" | grep -qxF "$e" \
        || { echo "ERROR: $INPUT_DIR/$e is in neither bundle list; add it to CORE or EXTERNAL" >&2; exit 1; }
done
[ -f "$ATTRIBUTION" ] || { echo "ERROR: $ATTRIBUTION not found" >&2; exit 1; }

mkdir -p "$OUT_DIR"
OUT_DIR="$(cd "$OUT_DIR" && pwd)"
STAGE="$(mktemp -d "$OUT_DIR/stage.XXXXXX")"
trap 'rm -rf "$STAGE"' EXIT
cp "$ATTRIBUTION" "$STAGE/ATTRIBUTION.md"

pack() {   # pack <bundle name> <entries...>
    local name="$1"; shift
    local archive="${PREFIX}_${name}_${VERSION}.tar.gz"
    local manifest="MANIFEST_${name}.tsv"
    echo "$(date +%T) $archive: hashing files"
    (cd "$INPUT_DIR" && printf 'path\tbytes\tsha256\n' \
        && find -L "$@" -type f | LC_ALL=C sort | while read -r f; do
               printf '%s\t%s\t%s\n' "$f" "$(stat -L -c %s "$f")" "$(sha256sum < "$f" | cut -d' ' -f1)"
           done) > "$STAGE/$manifest"
    echo "$(date +%T) $archive: $(($(wc -l < "$STAGE/$manifest") - 1)) files; packing"
    tar -czhf "$OUT_DIR/$archive" -C "$STAGE" "$manifest" ATTRIBUTION.md -C "$INPUT_DIR" "$@"
    (cd "$OUT_DIR" && sha256sum "$archive" > "$archive.sha256")
    echo "$(date +%T) wrote $OUT_DIR/$archive ($(du -h "$OUT_DIR/$archive" | cut -f1)), sha256 $(cut -d' ' -f1 "$OUT_DIR/$archive.sha256")"
}

pack core "${CORE[@]}"
pack external "${EXTERNAL[@]}"
