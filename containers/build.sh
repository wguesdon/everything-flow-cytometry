#!/usr/bin/env bash
set -euo pipefail

# Build the analysis image with Podman.
#
# The build installs the Bioconductor cytometry stack from source, so it takes a
# long time on the first run. Later builds reuse the layer cache.
#
# Usage:
#   ./containers/build.sh              # build with the pinned tag
#   ./containers/build.sh --no-cache   # rebuild every layer
#
# The image tag carries the Bioconductor release and the CRAN snapshot date, so
# two builds with different pins do not overwrite one another.

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"

BIOC_VERSION="3.23"
CRAN_SNAPSHOT="2026-08-15"
PYTHON_VERSION="3.12"

IMAGE="everything-flow-cytometry"
TAG="bioc${BIOC_VERSION}-cran${CRAN_SNAPSHOT}"

BUILD_FLAGS=()
for arg in "$@"; do
    case "$arg" in
        --no-cache) BUILD_FLAGS+=(--no-cache) ;;
        *) echo "Unknown option: $arg"; exit 1 ;;
    esac
done

if ! command -v podman >/dev/null 2>&1; then
    echo "Error: podman is not installed or is not on PATH."
    exit 1
fi

echo "Building ${IMAGE}:${TAG}"
echo "  Bioconductor  ${BIOC_VERSION}"
echo "  CRAN snapshot ${CRAN_SNAPSHOT}"
echo "  Python        ${PYTHON_VERSION}"
echo ""

podman build \
    "${BUILD_FLAGS[@]+"${BUILD_FLAGS[@]}"}" \
    --file "$REPO_ROOT/containers/Containerfile" \
    --build-arg "BIOC_VERSION=${BIOC_VERSION}" \
    --build-arg "CRAN_SNAPSHOT=${CRAN_SNAPSHOT}" \
    --build-arg "PYTHON_VERSION=${PYTHON_VERSION}" \
    --tag "${IMAGE}:${TAG}" \
    --tag "${IMAGE}:latest" \
    "$REPO_ROOT"

echo ""
echo "Built ${IMAGE}:${TAG}"
podman images "${IMAGE}" --format "  {{.Repository}}:{{.Tag}}  {{.Size}}"
echo ""
echo "Start it with:  ./containers/run.sh"
