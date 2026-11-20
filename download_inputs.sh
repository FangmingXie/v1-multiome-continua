#!/usr/bin/env bash
# Download the input bundles, check them, and unpack them into data/input/, so run_all.sh can run without
# data/raw/ and 0_prep. Bundles (made by make_input_bundle.sh; sources and licences in DATA_ATTRIBUTION.md):
#   core      every input except the external human datasets            -> folders 2, 3 (but 03), 4, 5
#   external  Jorstad23, Wang25, Fang22 human datasets (CC BY 4.0 / CC0)   -> folder 1, 3_regulon/03
# Downloads use gdown (pip install gdown), which handles Google Drive's large-file confirmation.
# Usage: bash download_inputs.sh [--core-only] [--from <dir>]
#   --core-only  skip the external bundle
#   --from DIR   use the archives already in DIR (e.g. downloaded in a browser) instead of downloading
set -euo pipefail

REPO_ROOT="$(cd "$(dirname "$0")" && pwd)"
INPUT_DIR="$REPO_ROOT/data/input"
DOWNLOAD_DIR="$REPO_ROOT/data/downloads"

VERSION="2026-11-20"
CORE_ARCHIVE="v1-multiome-continua_input_core_${VERSION}.tar.gz"
CORE_DRIVE_ID=""          # Google Drive file ID, set after upload
CORE_SHA256="840c9af2b632a876f44c56b9ee2da55d62ae809130150d9826d820f801c4c693"
EXTERNAL_ARCHIVE="v1-multiome-continua_input_external_${VERSION}.tar.gz"
EXTERNAL_DRIVE_ID=""      # Google Drive file ID, set after upload
EXTERNAL_SHA256="6dae15d06fd0143ede8fb34729e8d219322c299f749f32af99715da392bacf41"

CORE_ONLY=0
FROM=""
while [ $# -gt 0 ]; do
    case "$1" in
        --core-only) CORE_ONLY=1; shift ;;
        --from) FROM="$(cd "${2:?--from needs a directory}" && pwd)"; shift 2 ;;
        *) echo "ERROR: unknown argument $1 (usage: bash download_inputs.sh [--core-only] [--from <dir>])" >&2; exit 1 ;;
    esac
done

fetch() {   # fetch <archive> <drive id>: print the local path of the archive
    local archive="$1" id="$2"
    if [ -n "$FROM" ]; then
        [ -f "$FROM/$archive" ] || { echo "ERROR: $FROM/$archive not found" >&2; exit 1; }
        echo "$FROM/$archive"
        return
    fi
    [ -n "$id" ] || { echo "ERROR: no Google Drive ID set for $archive in $0" >&2; exit 1; }
    command -v gdown > /dev/null || { echo "ERROR: gdown not found; run: pip install gdown" >&2; exit 1; }
    mkdir -p "$DOWNLOAD_DIR"
    gdown "$id" -O "$DOWNLOAD_DIR/$archive" >&2
    echo "$DOWNLOAD_DIR/$archive"
}

install_bundle() {   # install_bundle <name> <archive> <drive id> <sha256>
    local name="$1" archive="$2" id="$3" sha="$4" path manifest n
    echo "$(date +%T) $name bundle: $archive"
    path="$(fetch "$archive" "$id")"
    echo "$sha  $path" | sha256sum -c --quiet - \
        || { echo "ERROR: $path does not match the expected sha256 $sha" >&2; exit 1; }
    mkdir -p "$INPUT_DIR"
    tar -xzf "$path" -C "$INPUT_DIR"
    manifest="$INPUT_DIR/MANIFEST_${name}.tsv"
    n=0
    while IFS=$'\t' read -r f bytes _; do
        [ -f "$INPUT_DIR/$f" ] && [ "$(stat -c %s "$INPUT_DIR/$f")" = "$bytes" ] \
            || { echo "ERROR: $INPUT_DIR/$f is missing or has the wrong size (see $manifest)" >&2; exit 1; }
        n=$((n + 1))
    done < <(tail -n +2 "$manifest")
    echo "$(date +%T) $name bundle: $n files in $INPUT_DIR"
    if [ -z "$FROM" ]; then rm -f "$path"; fi
}

install_bundle core "$CORE_ARCHIVE" "$CORE_DRIVE_ID" "$CORE_SHA256"
if [ "$CORE_ONLY" = 1 ]; then
    echo "Skipped the external bundle: folder 1 (1_mouse_human) and 3_regulon/03 (and its figure 07) need it."
else
    install_bundle external "$EXTERNAL_ARCHIVE" "$EXTERNAL_DRIVE_ID" "$EXTERNAL_SHA256"
fi
echo "data/input/ is ready; sources and licences: DATA_ATTRIBUTION.md"
