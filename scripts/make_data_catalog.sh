#!/usr/bin/env bash
set -euo pipefail

# Write docs/data_catalog.md from the local data/ folder.
#
# data/ is gitignored, so a fresh clone has no way to know what the S3 bucket
# holds. This catalogue is committed to git. Read it, choose the folders you need,
# then pull only those folders. That avoids a 100 GB transfer.
#
# Run this script again after you add or remove a dataset.
#
# Usage:
#   ./scripts/make_data_catalog.sh

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
DATA_DIR="$SCRIPT_DIR/data"
OUT="$SCRIPT_DIR/docs/data_catalog.md"

if [[ ! -d "$DATA_DIR" ]]; then
    echo "Error: there is no data/ folder. Run ./sync.sh pull first."
    exit 1
fi

TODAY="$(date +%Y-%m-%d)"

# Print one Markdown table row for every child folder of the given directory.
write_rows() {
    local parent="$1"
    local prefix="$2"
    local child size files fcs

    for child in "$parent"/*/; do
        [[ -d "$child" ]] || continue
        name="$(basename "$child")"
        size="$(du -sh "$child" 2>/dev/null | cut -f1)"
        files="$(find "$child" -type f 2>/dev/null | wc -l)"
        fcs="$(find "$child" -iname '*.fcs' -type f 2>/dev/null | wc -l)"
        printf '| `%s%s` | %s | %s | %s |\n' "$prefix" "$name" "$size" "$files" "$fcs"
    done
}

{
    echo "# Data catalogue"
    echo ""
    echo "The \`data/\` folder is gitignored and lives in S3. This page is committed, so"
    echo "you can choose what to pull before you transfer anything."
    echo ""
    echo "Generated on $TODAY by \`scripts/make_data_catalog.sh\`."
    echo ""
    echo "## How to pull one folder"
    echo ""
    echo '```bash'
    echo "./sync.sh catalog                                  # the same table, read from S3"
    echo "./sync.sh pull datasets/flowrepository/FR-FCM-ZZZU  # one accession"
    echo "./sync.sh pull literature repositories             # two folders at once"
    echo "./sync.sh pull                                     # everything, about 100 GB"
    echo '```'
    echo ""
    echo "A full pull transfers about 100 GB. Name the folders you need instead."
    echo ""
    echo "## Top level"
    echo ""
    echo "| Folder | Size | Files | FCS files |"
    echo "|---|---|---|---|"
    write_rows "$DATA_DIR" ""

    if [[ -d "$DATA_DIR/datasets" ]]; then
        echo ""
        echo "## datasets/"
        echo ""
        echo "| Folder | Size | Files | FCS files |"
        echo "|---|---|---|---|"
        write_rows "$DATA_DIR/datasets" "datasets/"
    fi

    if [[ -d "$DATA_DIR/datasets/flowrepository" ]]; then
        echo ""
        echo "## datasets/flowrepository/"
        echo ""
        echo "These are FlowRepository downloads. The source site stopped accepting new"
        echo "experiments in 2025 and its TLS certificate expired on 18 March 2023, so"
        echo "several of these accessions are hard to download again. Treat this copy as"
        echo "the working copy."
        echo ""
        echo "| Folder | Size | Files | FCS files |"
        echo "|---|---|---|---|"
        write_rows "$DATA_DIR/datasets/flowrepository" "datasets/flowrepository/"
    fi

    echo ""
    echo "## Total"
    echo ""
    printf '| Measure | Value |\n|---|---|\n'
    printf '| Size of `data/` | %s |\n' "$(du -sh "$DATA_DIR" 2>/dev/null | cut -f1)"
    printf '| Files | %s |\n' "$(find "$DATA_DIR" -type f 2>/dev/null | wc -l)"
    printf '| FCS files | %s |\n' "$(find "$DATA_DIR" -iname '*.fcs' -type f 2>/dev/null | wc -l)"
} > "$OUT"

echo "Wrote $OUT"
