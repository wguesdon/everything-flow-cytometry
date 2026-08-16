#!/usr/bin/env bash
set -euo pipefail

# Copy the flow cytometry archive from the WD1 external drive into data/.
#
# The source is the collection that was organised in 2025. It holds FlowRepository
# downloads, tool example data, cloned analysis repositories and one paper. Several
# of the FlowRepository accessions can no longer be downloaded from the source site,
# because FlowRepository stopped accepting new experiments in 2025 and its TLS
# certificate expired on 18 March 2023. This copy is therefore the working copy.
#
# The script copies. It does not delete the source. Run it again to resume an
# interrupted transfer, because rsync skips files that are already complete.
#
# Usage:
#   ./scripts/import_from_wd1.sh            # copy everything
#   ./scripts/import_from_wd1.sh --dry-run  # list what would move

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
SOURCE="/media/will/WD1/Organized_2025/Work/flow_cytometry"
DEST="$SCRIPT_DIR/data"

RSYNC_FLAGS=(-a --info=progress2 --human-readable)
for arg in "$@"; do
    case "$arg" in
        --dry-run) RSYNC_FLAGS+=(--dry-run) ;;
        *) echo "Unknown option: $arg"; exit 1 ;;
    esac
done

if [[ ! -d "$SOURCE" ]]; then
    echo "Error: the source is not present at $SOURCE"
    echo "Mount the WD1 drive first, then run this script again."
    exit 1
fi

mkdir -p "$DEST"

# Each pair is a source folder and the name it takes under data/.
copy_folder() {
    local from="$1"
    local to="$2"
    if [[ ! -d "$SOURCE/$from" ]]; then
        echo "Skipping $from, it is not in the source."
        return
    fi
    echo ""
    echo "Copying $from -> data/$to"
    mkdir -p "$DEST/$to"
    rsync "${RSYNC_FLAGS[@]}" "$SOURCE/$from/" "$DEST/$to/"
}

copy_folder "data" "datasets"
copy_folder "repositories" "repositories"
copy_folder "literature" "literature"

if [[ -f "$SOURCE/README.md" ]]; then
    echo ""
    echo "Copying README.md -> data/SOURCE_README.md"
    rsync "${RSYNC_FLAGS[@]}" "$SOURCE/README.md" "$DEST/SOURCE_README.md"
fi

echo ""
echo "Import complete. The source on WD1 was not changed."
du -sh "$DEST" 2>/dev/null || true
