#!/usr/bin/env bash
set -euo pipefail

# Move data/ to and from S3 so the archive can travel between machines.
#
# data/ is gitignored. It holds about 103 GB of FCS files, cloned analysis
# repositories and one paper. Git holds the code and the documentation. S3 holds
# the data.
#
# Several FlowRepository accessions in data/datasets/flowrepository/ can no longer
# be downloaded from the source site. FlowRepository stopped accepting new
# experiments in 2025 and its TLS certificate expired on 18 March 2023. Treat this
# archive as the working copy and keep the S3 copy current.
#
# The bucket is read from config.sh. To override it on one machine, export
# FLOWCYTO_S3_URI before you run this script.
#
# Usage:
#   ./sync.sh push                      # upload all of data/
#   ./sync.sh pull                      # download all of data/, about 100 GB
#   ./sync.sh pull datasets/flowrepository/FR-FCM-ZZZU
#                                       # download one accession only
#   ./sync.sh pull datasets/flowrepository/FR-FCM-ZZZU datasets/flowrepository/FR-FCM-ZZZV
#                                       # download two accessions
#   ./sync.sh catalog                   # size of every dataset folder in the bucket
#   ./sync.sh push --dry-run            # list the changes, transfer nothing
#   ./sync.sh list                      # show the folders at one level
#   ./sync.sh size                      # object count and total size
#
# A full pull transfers about 100 GB and costs about 9 US dollars in egress. Name
# the folders you need instead. Run 'catalog' first to see what each one costs, or
# read docs/data_catalog.md, which is committed to git and needs no transfer.
#
# --delete is push only and asks for confirmation first. It is refused on pull,
# because it would remove local files that are missing from the bucket and data/
# is gitignored, so this working copy can be the only copy.

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
LOCAL_ROOT="$SCRIPT_DIR/data"

usage() {
    echo "Usage: $0 {push|pull|list|size|catalog} [subpath ...] [--dry-run] [--delete]"
    echo ""
    echo "  push     upload data/ to \$FLOWCYTO_S3_URI"
    echo "  pull     download data/ from \$FLOWCYTO_S3_URI"
    echo "  list     show the folders at one level"
    echo "  size     show the object count and the total size"
    echo "  catalog  show the size of every dataset folder, so you can choose"
    echo ""
    echo "  subpath    limit the transfer to one or more folders under data/,"
    echo "             for example datasets/flowrepository/FR-FCM-ZZZU"
    echo "  --dry-run  show the changes without a transfer"
    echo "  --delete   remove remote files that are missing locally. Push only."
    echo ""
    echo "A full pull transfers about 100 GB. Name the folders you need instead."
    exit 1
}

if [[ $# -lt 1 ]]; then
    usage
fi

ACTION="$1"
shift

case "$ACTION" in
    push|pull|list|size|catalog) ;;
    *) usage ;;
esac

# Fall back to the committed configuration when the variable is not already set.
if [[ -z "${FLOWCYTO_S3_URI:-}" && -f "$SCRIPT_DIR/config.sh" ]]; then
    # shellcheck source=/dev/null
    source "$SCRIPT_DIR/config.sh"
fi

if [[ -z "${FLOWCYTO_S3_URI:-}" ]]; then
    echo "Error: FLOWCYTO_S3_URI is not set."
    echo "Set it in config.sh or export it, for example:"
    echo "  export FLOWCYTO_S3_URI=\"s3://my-bucket/everything-flow-cytometry\""
    exit 1
fi

if ! command -v aws >/dev/null 2>&1; then
    echo "Error: the aws CLI is not installed or is not on PATH."
    exit 1
fi

S3_BASE="${FLOWCYTO_S3_URI%/}"

SUBPATHS=()
SYNC_FLAGS=()
DRY_RUN=0
DELETE=0
for arg in "$@"; do
    case "$arg" in
        --dry-run) SYNC_FLAGS+=(--dryrun); DRY_RUN=1 ;;
        --delete)  DELETE=1 ;;
        --*) echo "Unknown option: $arg"; usage ;;
        *)
            clean="${arg#/}"
            clean="${clean%/}"
            # Accept a path copied from a shell prompt with the data/ prefix on it.
            clean="${clean#data/}"
            if [[ "$clean" == *".."* ]]; then
                echo "Error: '..' is not allowed in a subpath. Got '$arg'."
                exit 1
            fi
            if [[ -z "$clean" ]]; then
                echo "Error: '$arg' resolves to an empty subpath."
                exit 1
            fi
            SUBPATHS+=("$clean")
            ;;
    esac
done

# With no subpath the whole of data/ moves. The empty string makes the loop below
# run once over the root.
if [[ ${#SUBPATHS[@]} -eq 0 ]]; then
    SUBPATHS=("")
fi

if [[ ${#SUBPATHS[@]} -gt 1 && "$DELETE" -eq 1 ]]; then
    echo "Error: --delete takes one subpath at a time, so the list you confirm is"
    echo "the list that is removed. Run the command once per folder."
    exit 1
fi

SUBPATH="${SUBPATHS[0]}"

# aws s3 ls exits 1 when a prefix matches no object. That is an empty prefix and
# not a failure, so report it and exit 0. A real error, such as a missing bucket or
# a bad credential, exits with a higher code and is passed through.
run_ls() {
    local output status
    set +e
    output=$(aws s3 ls "$@" 2>&1)
    status=$?
    set -e

    if [[ $status -eq 0 ]]; then
        printf '%s\n' "$output"
        return 0
    fi

    if [[ $status -eq 1 ]]; then
        if [[ -n "$output" ]]; then
            printf '%s\n' "$output"
        fi
        echo "The prefix holds no object yet. Run '$0 push' to upload data/."
        return 0
    fi

    printf '%s\n' "$output" >&2
    return $status
}

if [[ "$ACTION" == "list" ]]; then
    run_ls "$S3_BASE/${SUBPATH:+$SUBPATH/}"
    exit $?
fi

if [[ "$ACTION" == "size" ]]; then
    run_ls --recursive --summarize "$S3_BASE/${SUBPATH:+$SUBPATH/}" \
        | grep -E "Total (Objects|Size)|holds no object"
    exit ${PIPESTATUS[0]}
fi

# catalog walks one level of folders and reports the size of each, so you can
# choose what to pull before any large transfer starts.
if [[ "$ACTION" == "catalog" ]]; then
    root="$S3_BASE/${SUBPATH:+$SUBPATH/}"
    echo "Folders under ${SUBPATH:-data}/ in $S3_BASE"
    echo ""
    printf '%-52s %10s %8s\n' "FOLDER" "SIZE" "OBJECTS"

    folders=$(aws s3 ls "$root" 2>/dev/null | awk '$1 == "PRE" {print $2}' || true)
    if [[ -z "$folders" ]]; then
        echo "No folder found. Run '$0 push' first, or give a subpath."
        exit 0
    fi

    for folder in $folders; do
        summary=$(aws s3 ls --recursive --summarize "$root$folder" 2>/dev/null || true)
        objects=$(printf '%s\n' "$summary" | awk '/Total Objects:/ {print $3}')
        bytes=$(printf '%s\n' "$summary" | awk '/Total Size:/ {print $3}')
        human=$(numfmt --to=iec --suffix=B "${bytes:-0}" 2>/dev/null || echo "${bytes:-0}")
        printf '%-52s %10s %8s\n' "${folder%/}" "$human" "${objects:-0}"
    done

    echo ""
    echo "Pull one with:  $0 pull ${SUBPATH:+$SUBPATH/}<folder>"
    exit 0
fi

LOCAL_PATH="$LOCAL_ROOT${SUBPATH:+/$SUBPATH}"
REMOTE_PATH="$S3_BASE${SUBPATH:+/$SUBPATH}"

# A pull with --delete removes local files that are absent from the bucket. There
# is no undo for that, so refuse it.
if [[ "$DELETE" -eq 1 && "$ACTION" == "pull" ]]; then
    echo "Error: --delete is not allowed on pull."
    echo ""
    echo "It would delete local files that are missing from S3, and data/ is"
    echo "gitignored, so this machine can hold the only copy."
    echo "To see what a pull would bring down, use:  $0 pull --dry-run"
    echo "To remove stale remote files instead, use: $0 push --delete"
    exit 1
fi

# Count what --delete would remove, show the list, and ask before the transfer.
if [[ "$DELETE" -eq 1 && "$DRY_RUN" -eq 0 ]]; then
    echo "Checking what --delete would remove from $REMOTE_PATH ..."
    pending=$(aws s3 sync "$LOCAL_PATH" "$REMOTE_PATH" --dryrun --delete 2>/dev/null \
        | grep '^(dryrun) delete:' || true)

    if [[ -z "$pending" ]]; then
        echo "Nothing to delete. The push continues as a normal push."
        DELETE=0
    else
        count=$(printf '%s\n' "$pending" | wc -l)
        echo ""
        printf '%s\n' "$pending"
        echo ""
        echo "This deletes $count remote object(s) from the list above."
        echo "The bucket has versioning enabled, so a deletion can be restored."
        echo "Read the list before you continue."

        if [[ ! -t 0 ]]; then
            echo "Error: --delete needs an interactive terminal. The script stops here."
            exit 1
        fi

        read -r -p "Type 'delete' to continue: " reply
        if [[ "$reply" != "delete" ]]; then
            echo "Stopped. Nothing was deleted and nothing was uploaded."
            exit 1
        fi
    fi
fi

if [[ "$DELETE" -eq 1 ]]; then
    SYNC_FLAGS+=(--delete)
fi

if [[ "$ACTION" == "push" ]]; then
    SYNC_FLAGS+=(--storage-class "${FLOWCYTO_STORAGE_CLASS:-STANDARD}")
fi

for sub in "${SUBPATHS[@]}"; do
    local_path="$LOCAL_ROOT${sub:+/$sub}"
    remote_path="$S3_BASE${sub:+/$sub}"

    if [[ "$ACTION" == "push" ]]; then
        if [[ ! -d "$local_path" ]]; then
            echo "Error: there is no local folder at $local_path"
            exit 1
        fi
        echo ""
        echo "Pushing ${sub:-data} -> $remote_path"
        aws s3 sync "$local_path" "$remote_path" ${SYNC_FLAGS[@]+"${SYNC_FLAGS[@]}"}
    else
        echo ""
        echo "Pulling $remote_path -> data/${sub}"
        mkdir -p "$local_path"
        aws s3 sync "$remote_path" "$local_path" ${SYNC_FLAGS[@]+"${SYNC_FLAGS[@]}"}
    fi
done

echo ""
echo "Sync ($ACTION) complete for ${#SUBPATHS[@]} path(s)."
