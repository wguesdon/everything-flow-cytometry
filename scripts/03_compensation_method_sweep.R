#!/usr/bin/env Rscript

# Why the spillover computation uses median and not the flowStats default.
#
# flowStats::spillover_ng() defaults to method = "mode". That suits a bead
# control, where nearly every event is positive. OMIP-39 uses cell controls, and
# four of its eleven markers sit on a minority of cells, so the mode lands on
# the
# negative population and the ratio it produces is meaningless.
#
# This script measures that instead of asserting it. It computes the matrix
# under
# each combination of method, pregate and normalisation filter, and scores every
# result against the matrix the instrument stored for the same experiment. The
# stored matrix is a fair reference: it was built from these same controls, by
# the
# acquisition software, with an operator checking it.
#
# The output feeds the compensation section of
# reports/omip39_automated_gating.qmd,
# so the recommendation in that report rests on a table this script produced.
#
# Run it in the container:
#   podman run --rm -v "$PWD:/work:z" -w /work everything-flow-cytometry:latest
# \
#     Rscript scripts/03_compensation_method_sweep.R

suppressPackageStartupMessages({
  library(flowCore)
  library(flowStats)
})

for (file in list.files("R", pattern = "\\.R$", full.names = TRUE)) {
  source(file)
}

kDataDir <- file.path(
  "data", "datasets", "flowrepository", "OMIP-39",
  "FlowRepository_FR-FCM-ZYY6_files"
)
kOutputDir <- file.path("output", "omip39")

dir.create(kOutputDir, recursive = TRUE, showWarnings = FALSE)

Log <- function(...) {
  cat(format(Sys.time(), "%H:%M:%S"), "|", ..., "\n")
}

if (!dir.exists(kDataDir)) {
  stop("OMIP-39 is not present. Pull it with:\n",
       "  ./sync.sh pull datasets/flowrepository/OMIP-39")
}

Log("Reading the controls and the sample")
control_set <- read.flowSet(
  list.files(kDataDir, pattern = "^Single stainings.*\\.fcs$", full.names = TRUE),
  truncate_max_range = FALSE
)
sample_set <- read.flowSet(
  list.files(kDataDir, pattern = "^Samples_.*\\.fcs$", full.names = TRUE),
  truncate_max_range = FALSE
)

match_table <- MatchControlsToChannels(control_set, unstained_pattern = "unstained")
match_file <- file.path(kOutputDir, "spillover_match.csv")
WriteMatchFile(match_table, match_file)

stored <- ExtractSpillover(sample_set[[1]])
Log("Reference is the stored matrix,", nrow(stored), "by", ncol(stored))

#' Score a computed matrix against the stored one
#'
#' @param computed A computed spillover matrix.
#' @param reference The stored spillover matrix.
#' @return A one row `data.frame` of summary statistics over the detectors the
#'   two matrices share, with every difference in percentage points.
ScoreAgainstStored <- function(computed, reference) {
  shared <- intersect(colnames(computed), colnames(reference))
  a <- computed[shared, shared]
  b <- reference[shared, shared]
  off_diagonal <- row(a) != col(a)

  data.frame(
    channels = length(shared),
    max_difference_points = 100 * max(abs(a[off_diagonal] - b[off_diagonal])),
    median_difference_points = 100 * stats::median(abs(a[off_diagonal] - b[off_diagonal])),
    values_above_100_percent = sum(a[off_diagonal] > 1),
    correlation = stats::cor(a[off_diagonal], b[off_diagonal]),
    stringsAsFactors = FALSE
  )
}

settings <- expand.grid(
  method = c("mode", "median"),
  pregate = c(TRUE, FALSE),
  use_norm_filter = c(TRUE, FALSE),
  stringsAsFactors = FALSE
)

Log("Running", nrow(settings), "configurations")

rows <- lapply(seq_len(nrow(settings)), function(i) {
  setting <- settings[i, ]
  label <- sprintf("method=%s pregate=%s normFilter=%s",
                   setting$method, setting$pregate, setting$use_norm_filter)

  outcome <- tryCatch({
    computed <- flowStats::spillover_ng(
      control_set,
      fsc = "FSC-A",
      ssc = "SSC-A",
      matchfile = match_file,
      method = setting$method,
      pregate = setting$pregate,
      useNormFilt = setting$use_norm_filter,
      plot = FALSE
    )
    score <- ScoreAgainstStored(computed, stored)
    score$status <- "ok"
    score
  }, error = function(e) {
    data.frame(
      channels = NA_integer_,
      max_difference_points = NA_real_,
      median_difference_points = NA_real_,
      values_above_100_percent = NA_integer_,
      correlation = NA_real_,
      status = paste("failed:", conditionMessage(e)),
      stringsAsFactors = FALSE
    )
  })

  # A configuration can return without an error and still produce a matrix that
  # holds no usable value, which shows up as a missing correlation.
  if (is.na(outcome$correlation) && outcome$status == "ok") {
    outcome$status <- "produced no usable estimate"
  }

  Log(" ", label, "->", outcome$status)
  cbind(setting, outcome)
})

sweep <- do.call(rbind, rows)
sweep <- sweep[order(sweep$correlation, decreasing = TRUE, na.last = TRUE), ]
rownames(sweep) <- NULL

write.csv(sweep, file.path(kOutputDir, "method_sweep.csv"), row.names = FALSE)

cat("\n")
print(sweep, digits = 3)

usable <- sweep[!is.na(sweep$correlation), ]
if (nrow(usable) > 0) {
  best <- usable[1, ]
  Log("Best configuration: method =", best$method,
      "pregate =", best$pregate,
      "useNormFilt =", best$use_norm_filter,
      sprintf("(correlation %.3f, %d value(s) above 100 percent)",
              best$correlation, best$values_above_100_percent))
}

Log("Wrote", file.path(kOutputDir, "method_sweep.csv"))
