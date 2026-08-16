#!/usr/bin/env bash
set -uo pipefail

# Retrieve the paper that describes each package in docs/packages.md.
#
# CLAUDE.md forbids a citation written from memory. This script queries Europe
# PMC and prints the record, so every entry in the "Papers behind the packages"
# section of docs/literature.md traces to a command that you can run again.
#
# The script does not edit any document. It writes a JSON file and a table. You
# compare the table against docs/literature.md by hand.
#
# A query that returns hitCount=0 means that Europe PMC holds no article with
# that title. Record the package under "Packages with no paper".
#
# Usage:
#   ./scripts/verify_package_papers.sh                     # table to the terminal
#   ./scripts/verify_package_papers.sh out/papers.json     # table plus a JSON file

API="https://www.ebi.ac.uk/europepmc/webservices/rest/search"
OUT="${1:-}"

if ! command -v jq > /dev/null 2>&1; then
    echo "Error: this script needs jq. Install it and run the script again." >&2
    exit 1
fi

# package<TAB>Europe PMC query
QUERIES=$(cat <<'EOF'
flowCore	TITLE:"flowCore"
flowWorkspace	TITLE:"flowWorkspace"
CytoML	TITLE:"CytoML"
ggcyto	TITLE:"ggCyto"
PeacoQC	TITLE:"PeacoQC"
flowAI	TITLE:"flowAI"
CytoNorm	TITLE:"CytoNorm"
openCyto	TITLE:"OpenCyto"
flowDensity	TITLE:"flowDensity"
flowClust	TITLE:"flowClust"
flowStats	TITLE:"flowStats"
FlowSOM	TITLE:"FlowSOM: Using self-organizing maps"
FlowSOM_python	TITLE:"Efficient cytometry analysis with FlowSOM in Python"
diffcyt	TITLE:"diffcyt"
CATALYST_spillover	TITLE:"Compensation of Signal Spillover in Suspension and Imaging Mass Cytometry"
CATALYST_pipeline	TITLE:"An R-based reproducible and user-friendly preprocessing pipeline for CyTOF data"
cydar	TITLE:"Testing for differential abundance in mass cytometry data"
HDCytoData	TITLE:"HDCytoData"
FlowRepositoryR	TITLE:"FlowRepositoryR"
FlowKit	TITLE:"FlowKit"
readfcs	TITLE:"readfcs"
Pytometry	TITLE:"Pytometry"
cytoflow	TITLE:"Cytoflow"
scanpy	TITLE:"SCANPY: large-scale single-cell gene expression data analysis"
EOF
)

TMP=$(mktemp)
trap 'rm -f "$TMP"' EXIT

echo "[" > "$TMP"
FIRST=1
printf 'package\thits\tpmid\tpmcid\tyear\tjournal\tvolume(issue):pages\tdoi\ttitle\n'

while IFS=$'\t' read -r PKG QUERY; do
    [[ -z "$PKG" ]] && continue
    ENCODED=$(printf '%s' "$QUERY" | jq -sRr @uri)
    BODY=$(curl -s --max-time 30 "${API}?query=${ENCODED}&format=json&pageSize=5&resultType=core")
    if [[ -z "$BODY" ]]; then
        echo "Error: no answer from Europe PMC for $PKG. Check the network." >&2
        continue
    fi
    HITS=$(printf '%s' "$BODY" | jq -r '.hitCount // 0')

    printf '%s' "$BODY" | jq -r --arg pkg "$PKG" --arg hits "$HITS" '
        if (.resultList.result | length) == 0 then
            "\($pkg)\t\($hits)\t-\t-\t-\t-\t-\t-\tNO RECORD"
        else
            .resultList.result[] |
            "\($pkg)\t\($hits)\t\(.pmid // .id)\t\(.pmcid // "-")\t\(.pubYear)\t" +
            "\(.journalInfo.journal.title // "preprint")\t" +
            "\(.journalInfo.volume // "-")(\(.journalInfo.issue // "-")):\(.pageInfo // "-")\t" +
            "\(.doi // "-")\t\(.title)"
        end'

    [[ "$FIRST" -eq 0 ]] && echo "," >> "$TMP"
    FIRST=0
    printf '%s' "$BODY" | jq --arg pkg "$PKG" --arg q "$QUERY" \
        '{package: $pkg, query: $q, hitCount: (.hitCount // 0), results: .resultList.result}' >> "$TMP"
    sleep 1
done <<< "$QUERIES"

echo "]" >> "$TMP"

if [[ -n "$OUT" ]]; then
    mkdir -p "$(dirname "$OUT")"
    cp "$TMP" "$OUT"
    echo "Wrote the full records to $OUT" >&2
fi
