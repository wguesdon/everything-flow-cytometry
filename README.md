# everything-flow-cytometry

This repository collects the literature, the software packages and the example code for computational
flow cytometry analysis in R and Python.

## Why this repository exists

Panels grew faster than the analysis method. A three colour assay gives 3 possible bi-axial plots. A
ten colour panel gives 45. A 45 marker spectral panel gives 990. Manual gating does not scale to that
number, and the operator who draws the gates is a large source of variation between laboratories.

The software to replace manual gating exists. Most of it is command line R, and most clinical
laboratories do not use it. This repository puts the reading list, the package list and runnable
example code in one place.

## Status

The repository is at an early stage. The table below states what is ready and what is not.

| Item | State |
|---|---|
| `docs/literature.md` | Ready. Every citation was verified against a retrieved record on 2026-08-13. |
| `docs/packages.md` | Partial. The package list is complete for the tools in daily use. The version and maintenance status is checked for some entries only. |
| `docs/datasets.md` | Ready. Every repository count was checked on 2026-08-13. |
| `examples/` | Empty. It needs a public dataset and a container. |
| Agent skills | Not started. They come after the example code runs. |

## Contents

| Path | Content |
|---|---|
| `docs/literature.md` | Key reviews, the gating variability studies and the data sharing papers |
| `docs/packages.md` | R and Python packages, grouped by the step they perform |
| `docs/datasets.md` | Public repositories that hold FCS files, with their current state |

## Planned work

1. Add one worked example on a public dataset. The example reads FCS files, applies compensation and
   transformation, runs quality control, gates the data and finds populations.
2. Put the example in a container so the result is the same on every machine. Use Podman.
3. Add agent skills for Claude Code, Codex and OpenCode. Each skill runs the example code and reports
   the result.

## Provenance

The citations and the repository counts come from a literature search that was run on 2026-08-13.
Each claim was checked against Europe PMC, PubMed, Crossref, Semantic Scholar or the PMC full text. A
claim that failed verification is not in this repository.
