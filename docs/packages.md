# Packages

This page lists the R and Python packages for flow cytometry analysis. The packages are grouped by the
step that they perform.

Read the status column before you start a project. A package that has no recent release can still work,
but you must pin the version and test it.

The paper column gives the first author and the year of the paper that describes the package. The
section "Papers behind the packages" in [`literature.md`](literature.md#papers-behind-the-packages)
holds the full citation, the PMID and the open access state of each one. An empty paper cell means
that a search of Europe PMC on 2026-08-16 found no article for that package.

## Status of this page

The list of packages in daily use is complete. The status column is checked for some entries only. The
date of the check is in the column. An entry with no date was not checked.

The paper column is complete. Every reference in it was retrieved from Europe PMC on 2026-08-16.

## R packages

### Read and prepare the data

| Package | Source | Step | Paper | Status |
|---|---|---|---|---|
| flowCore | Bioconductor | Read FCS files. Hold the data as a flowFrame or a flowSet. | Hahne 2009 | |
| flowWorkspace | Bioconductor | Hold a gating hierarchy as a GatingSet. | None. Cite Finak 2018 or Finak 2014. | |
| CytoML | Bioconductor | Import a FlowJo, Diva or Cytobank workspace into a GatingSet. | Finak 2018 | |
| ggcyto | Bioconductor | Plot a flowSet or a GatingSet with ggplot2 grammar. | Van 2018 | |

### Quality control and normalisation

| Package | Source | Step | Paper | Status |
|---|---|---|---|---|
| PeacoQC | Bioconductor | Remove unstable events before analysis. | Emmaneel 2022 | |
| flowAI | Bioconductor | Remove anomalies in flow rate, signal and dynamic range. | Monaco 2016 | |
| CytoNorm | GitHub only, saeyslab/CytoNorm | Normalise batches with control samples. | Van Gassen 2020. Quintelier 2025 covers batches with no control sample. | Checked 2026-08-13. It is not on CRAN or Bioconductor. |

### Gating

| Package | Source | Step | Paper | Status |
|---|---|---|---|---|
| openCyto | Bioconductor | Run sequential automated gating from a template. | Finak 2014 | |
| flowDensity | Bioconductor | Set a gate from the density of one or two channels. | Malek 2015 | |
| flowClust | Bioconductor | Model based clustering for gating. | Lo 2009 | |

### Find populations and test differences

| Package | Source | Step | Paper | Status |
|---|---|---|---|---|
| FlowSOM | Bioconductor | Cluster events with a self-organizing map. Build a minimum spanning tree. | Van Gassen 2015 | |
| diffcyt | Bioconductor | Test for differential abundance and differential state between groups. | Weber 2019, diffcyt | |
| CATALYST | Bioconductor | Preprocess and analyse mass cytometry data. | Chevrier 2018 for spillover. Crowell 2020 for the pipeline. | |
| cydar | Bioconductor | Test differential abundance in hyperspheres. | Lun 2017 | |

### Additional methods

| Package | Source | Step | Paper | Status |
|---|---|---|---|---|
| flowStats | Bioconductor | Analyse flow data beyond the basic infrastructure that flowCore provides. | None | Checked 2026-08-16. Version 4.24.0 at Bioconductor 3.23. The maintainers are Greg Finak and Mike Jiang. |

### Data access

| Package | Source | Step | Paper | Status |
|---|---|---|---|---|
| HDCytoData | Bioconductor | Load ten standard benchmark datasets from ExperimentHub. | Weber 2019, HDCytoData | Checked 2026-08-13. Version 1.32.1. The maintainer is Lukas Weber. |
| FlowRepositoryR | Removed | Download from FlowRepository. | None. Cite Spidlen 2012 for the repository. | Checked 2026-08-13. Bioconductor removed it at release 3.14 in October 2021. No supported client exists. |

## Python packages

| Package | Source | Step | Paper | Status |
|---|---|---|---|---|
| FlowKit | PyPI | Read FCS files. Apply compensation and transformation. Import a FlowJo workspace. | White 2021 | |
| readfcs | PyPI | Read an FCS file into an AnnData object. | None | |
| flowsom | PyPI | Run the FlowSOM algorithm in Python. The PyPI summary calls it the complete FlowSOM package known from R. | Couckuyt 2024 | Checked 2026-08-16. PyPI version 0.2.2. |
| Pytometry | PyPI | Analyse flow and mass cytometry data in the scanpy ecosystem. | Büttner 2022, a preprint | Installed at 0.1.6. Used by `reports/omip58_pytometry.qmd`. Checked 2026-08-18. |
| cytoflow | PyPI and GitHub | Build a reproducible analysis workflow with a graphical interface. | Teague 2026 | Checked 2026-08-13. The PyPI record is four years behind the GitHub release. Install from GitHub. |
| scanpy | PyPI | Downstream analysis on an AnnData object. | Wolf 2018 | |

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
   `scripts/verify_package_papers.sh` does this for the paper column already.
