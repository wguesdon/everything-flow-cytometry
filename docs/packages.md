# Packages

This page lists the R and Python packages for flow cytometry analysis. The packages are grouped by the
step that they perform.

Read the status column before you start a project. A package that has no recent release can still work,
but you must pin the version and test it.

## Status of this page

The list of packages in daily use is complete. The status column is checked for some entries only. The
date of the check is in the column. An entry with no date was not checked.

## R packages

### Read and prepare the data

| Package | Source | Step | Status |
|---|---|---|---|
| flowCore | Bioconductor | Read FCS files. Hold the data as a flowFrame or a flowSet. | |
| flowWorkspace | Bioconductor | Hold a gating hierarchy as a GatingSet. | |
| CytoML | Bioconductor | Import a FlowJo, Diva or Cytobank workspace into a GatingSet. | |
| ggcyto | Bioconductor | Plot a flowSet or a GatingSet with ggplot2 grammar. | |

### Quality control and normalisation

| Package | Source | Step | Status |
|---|---|---|---|
| PeacoQC | Bioconductor | Remove unstable events before analysis. | |
| flowAI | Bioconductor | Remove anomalies in flow rate, signal and dynamic range. | |
| CytoNorm | GitHub only, saeyslab/CytoNorm | Normalise batches with control samples. | Checked 2026-08-13. It is not on CRAN or Bioconductor. |

### Gating

| Package | Source | Step | Status |
|---|---|---|---|
| openCyto | Bioconductor | Run sequential automated gating from a template. | |
| flowDensity | Bioconductor | Set a gate from the density of one or two channels. | |
| flowClust | Bioconductor | Model based clustering for gating. | |

### Find populations and test differences

| Package | Source | Step | Status |
|---|---|---|---|
| FlowSOM | Bioconductor | Cluster events with a self-organizing map. Build a minimum spanning tree. | |
| diffcyt | Bioconductor | Test for differential abundance and differential state between groups. | |
| CATALYST | Bioconductor | Preprocess and analyse mass cytometry data. | |
| cydar | Bioconductor | Test differential abundance in hyperspheres. | |

### Data access

| Package | Source | Step | Status |
|---|---|---|---|
| HDCytoData | Bioconductor | Load ten standard benchmark datasets from ExperimentHub. | Checked 2026-08-13. Version 1.32.1. The maintainer is Lukas Weber. |
| FlowRepositoryR | Removed | Download from FlowRepository. | Checked 2026-08-13. Bioconductor removed it at release 3.14 in October 2021. No supported client exists. |

## Python packages

| Package | Source | Step | Status |
|---|---|---|---|
| FlowKit | PyPI | Read FCS files. Apply compensation and transformation. Import a FlowJo workspace. | |
| readfcs | PyPI | Read an FCS file into an AnnData object. | |
| Pytometry | PyPI | Analyse flow and mass cytometry data in the scanpy ecosystem. | Checked 2026-08-13. The last PyPI release is from September 2024. |
| cytoflow | PyPI and GitHub | Build a reproducible analysis workflow with a graphical interface. | Checked 2026-08-13. The PyPI record is four years behind the GitHub release. Install from GitHub. |
| scanpy | PyPI | Downstream analysis on an AnnData object. | |

## Packages that no longer work

Do not start a new project with these. The check was run on 2026-08-13.

| Package | Language | State |
|---|---|---|
| CytoPy | Python | Dead or removed. |
| FlowCytometryTools | Python | Dead or removed. |
| Rphenograph | R | Dead or removed. |
| tidytof | R | Dead or removed. |
| cytofkit | R | Dead or removed. |
| FlowRepositoryR | R | Removed from Bioconductor in October 2021. |

## What the field actually uses

Cheung and colleagues surveyed 51 tools and 49 clinical laboratory respondents in 2021. 59% of the
tools are written in R. 41% of the tools have a graphical interface. 26 of the 49 respondents never
use automated software.

That gap is the reason this repository exists. The packages above are almost all command line R, and
most laboratories do not run them.

## To do

1. Re-run the package inventory against live Bioconductor, CRAN, PyPI and GitHub records.
2. Record the version and the last release date for every package in the table.
3. Write the inventory as a script so the check is repeatable.
