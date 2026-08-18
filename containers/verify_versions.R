#!/usr/bin/env Rscript

# Assert the version of every pinned R package and stop the build on a mismatch.
#
# The Containerfile pins CRAN to a dated snapshot and Bioconductor to one release.
# Those two pins should give the versions below every time. This script turns that
# expectation into a check, so an upstream change breaks the build instead of
# changing a result quietly.
#
# The versions were read on 2026-08-16 from two sources:
#   CRAN         https://packagemanager.posit.co/cran/2026-08-15
# The four figure and manifest packages were read on 2026-08-18 from the image
# built against that same snapshot, because they arrived transitively before
# they were named here.
#   Bioconductor https://bioconductor.org/packages/3.23/bioc/VIEWS
#
# To move a pin, change the number here and change the snapshot date or the
# Bioconductor release in the Containerfile at the same time.

expected <- c(
  # CRAN, from the 2026-08-15 snapshot
  IRkernel    = "1.3.2",
  IRdisplay   = "1.1",
  repr        = "1.1.7",
  ggplot2     = "4.0.3",
  dplyr       = "1.2.1",
  tidyr       = "1.3.2",
  readr       = "2.2.0",
  here        = "1.0.2",
  patchwork   = "1.3.2",
  matrixStats = "1.5.0",
  uwot        = "0.2.4",
  Rtsne       = "0.17",
  pheatmap    = "1.0.13",
  ggpubr      = "1.0.0",
  rstatix     = "1.1.0",
  # R/figures.R writes every raster figure through ragg and picks the font
  # through systemfonts, and textshaping draws the glyphs. A change in any of
  # the three changes a figure, so all three are pinned. R/cytokit.R writes a
  # bundle manifest with jsonlite.
  ragg        = "1.5.2",
  systemfonts = "1.3.2",
  textshaping = "1.0.5",
  jsonlite    = "2.0.0",
  testthat    = "3.3.2",
  roxygen2    = "8.1.0",
  knitr       = "1.51",
  rmarkdown   = "2.31",
  withr       = "3.0.3",

  # Bioconductor 3.23
  flowCore             = "2.24.0",
  flowWorkspace        = "4.24.0",
  openCyto             = "2.24.0",
  CytoML               = "2.24.0",
  ggcyto               = "1.40.0",
  flowAI               = "1.42.0",
  PeacoQC              = "1.22.0",
  flowClust            = "3.50.0",
  flowDensity          = "1.46.0",
  flowStats            = "4.24.0",
  flowViz              = "1.76.0",
  FlowSOM              = "2.20.0",
  diffcyt              = "1.32.1",
  CATALYST             = "1.36.0",
  SingleCellExperiment = "1.34.0",
  BiocParallel         = "1.46.0",
  ExperimentHub        = "3.2.0",
  HDCytoData           = "1.32.1"
)

installed <- installed.packages()[, "Version"]

missing <- setdiff(names(expected), names(installed))
if (length(missing) > 0) {
  cat("MISSING PACKAGES:\n")
  cat(paste0("  ", missing, collapse = "\n"), "\n")
  quit(status = 1)
}

found <- installed[names(expected)]
mismatched <- names(expected)[found != expected]

cat(sprintf("%-24s %-12s %-12s %s\n", "PACKAGE", "EXPECTED", "INSTALLED", "STATUS"))
for (pkg in names(expected)) {
  status <- if (found[[pkg]] == expected[[pkg]]) "ok" else "MISMATCH"
  cat(sprintf("%-24s %-12s %-12s %s\n", pkg, expected[[pkg]], found[[pkg]], status))
}

cat("\nR version:            ", as.character(getRversion()), "\n")
cat("Bioconductor version: ", as.character(BiocManager::version()), "\n")

if (length(mismatched) > 0) {
  cat("\nBuild stopped.", length(mismatched), "package(s) do not match the pin.\n")
  cat("Either the CRAN snapshot date moved or the Bioconductor release changed.\n")
  cat("Update containers/verify_versions.R and the Containerfile together.\n")
  quit(status = 1)
}

cat("\nAll", length(expected), "pinned packages match.\n")
