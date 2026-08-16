#!/usr/bin/env bash
set -euo pipefail

# Start JupyterLab in the analysis container.
#
# The port is published on 127.0.0.1 only. The notebook is therefore not reachable
# from the network, and you connect from another computer through an SSH tunnel.
# On that other computer:
#
#   ssh -N -L 8888:localhost:8888 will@<this-host>
#
# Then open the URL that this script prints, on the other computer, at
# http://localhost:8888. The token is in the URL.
#
# Usage:
#   ./containers/run.sh                 # start in the foreground
#   ./containers/run.sh --detach        # start in the background
#   ./containers/run.sh --port 9999     # use another host port
#   ./containers/run.sh --shell         # open a shell instead of the notebook
#   ./containers/run.sh --stop          # stop and remove the container
#
# Two folders are mounted. examples/ and notebooks/ are read and write, so your
# work is saved on the host. data/ is mounted read only, because the archive is
# the input and nothing in an example should change it.

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"

IMAGE="everything-flow-cytometry:latest"
NAME="flow-cytometry-lab"
HOST_PORT="8888"
DETACH=0
SHELL_MODE=0

while [[ $# -gt 0 ]]; do
    case "$1" in
        --detach|-d) DETACH=1; shift ;;
        --shell)     SHELL_MODE=1; shift ;;
        --port)      HOST_PORT="${2:?--port needs a number}"; shift 2 ;;
        --stop)
            podman rm -f "$NAME" 2>/dev/null && echo "Stopped and removed $NAME" \
                || echo "No container named $NAME is running."
            exit 0
            ;;
        *) echo "Unknown option: $1"; exit 1 ;;
    esac
done

if ! command -v podman >/dev/null 2>&1; then
    echo "Error: podman is not installed or is not on PATH."
    exit 1
fi

if ! podman image exists "$IMAGE"; then
    echo "Error: the image $IMAGE does not exist."
    echo "Build it first with:  ./containers/build.sh"
    exit 1
fi

if podman container exists "$NAME"; then
    echo "A container named $NAME already exists."
    echo "Stop it with:  $0 --stop"
    exit 1
fi

mkdir -p "$REPO_ROOT/notebooks"

MOUNTS=(
    --volume "$REPO_ROOT/examples:/home/analysis/examples:z"
    --volume "$REPO_ROOT/notebooks:/home/analysis/notebooks:z"
    --volume "$REPO_ROOT/docs:/home/analysis/docs:ro,z"
)

# data/ is only mounted when it is present, so the container still starts on a
# machine that has not pulled the archive.
if [[ -d "$REPO_ROOT/data" ]]; then
    MOUNTS+=(--volume "$REPO_ROOT/data:/home/analysis/data:ro,z")
else
    echo "Note: data/ is not present. Pull what you need with ./sync.sh pull <folder>."
fi

RUN_FLAGS=(
    --name "$NAME"
    --publish "127.0.0.1:${HOST_PORT}:8888"
    --userns keep-id
    "${MOUNTS[@]}"
)

if [[ "$SHELL_MODE" -eq 1 ]]; then
    echo "Opening a shell in $IMAGE"
    podman run --rm -it "${RUN_FLAGS[@]}" "$IMAGE" /bin/bash
    exit 0
fi

echo "Starting JupyterLab"
echo "  image      $IMAGE"
echo "  bound to   127.0.0.1:${HOST_PORT}, not to the network"
echo "  examples   read and write"
echo "  data       read only"
echo ""
echo "From another computer, open the tunnel first:"
echo "  ssh -N -L ${HOST_PORT}:localhost:${HOST_PORT} ${USER}@$(hostname)"
echo "then open the printed URL there, with localhost in place of the host name."
echo ""

if [[ "$DETACH" -eq 1 ]]; then
    podman run --detach "${RUN_FLAGS[@]}" "$IMAGE"
    echo "Started in the background. Read the token with:"
    echo "  podman logs $NAME 2>&1 | grep token="
    echo "Stop it with:  $0 --stop"
else
    podman run --rm -it "${RUN_FLAGS[@]}" "$IMAGE"
fi
