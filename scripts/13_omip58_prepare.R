# The R half of the OMIP-058 workflow.
#
# FR-FCM-ZYRN holds two PBMC files on a 28 colour panel and 59 compensation
# controls, 55 of them single stains. This script computes the spillover matrix
# from those
# controls, applies it, removes the out of range events, the doublets, the
# debris and the dead cells, splits the T cells from the rest, and writes the
# result to disk as FCS.
#
# The Python half reads what this writes. Nothing below CD3 is gated here, and
# scripts/14_omip58_pytometry.py explains why that boundary sits where it does.
#
# Run it in the container:
#   podman run --rm -v "$PWD:/work:z" -w /work everything-flow-cytometry:latest
# \
#     Rscript scripts/13_omip58_prepare.R

suppressPackageStartupMessages({
  library(flowCore)
  library(flowStats)
  library(ggplot2)
  library(robustbase)
})

source(file.path("R", "figures.R"))
source(file.path("R", "io.R"))
source(file.path("R", "compensation.R"))
source(file.path("R", "spillover_compute.R"))
source(file.path("R", "naive_memory.R"))
source(file.path("R", "panels.R"))
source(file.path("R", "omip58.R"))

# Several steps below draw random subsets. robustbase::covMcd does, and so does
# flowStats::norm2Filter, which spillover_ng calls while it gates each control.
# Without a seed the matrix moves a little between runs, the compensated values
# move with it, and a count in this report cannot be reproduced. Two unseeded
# runs differed by 6 live lymphocytes of 598,880 and swung the fitted Vg9 cut
# from 60.97 to 79.76 percent of T cells.
kSeed <- 42
set.seed(kSeed)

kDepositDir <- file.path("data", "datasets", "flowrepository", "FR-FCM-ZYRN")
kOutputDir <- file.path("output", "omip58")
kHandoffDir <- file.path(kOutputDir, "handoff")
kUnstainedCells <- "Comp_Cells_unstained_D8_D08_057.fcs"
kUnstainedBeads <- "Comp_Beads_unstained mouse_D4_D04_028.fcs"

dir.create(kHandoffDir, recursive = TRUE, showWarnings = FALSE)
Write <- function(x, name) {
  utils::write.csv(x, file.path(kOutputDir, name), row.names = FALSE)
  invisible(x)
}
Say <- function(...) cat(..., "\n", sep = "")

# ---------------------------------------------------------------------------
# Part 1. What the deposit says about its own compensation.
# ---------------------------------------------------------------------------

Say("Part 1: the compensation state of the deposit")

sample_files <- list.files(kDepositDir, pattern = "^PBMC", full.names = TRUE)
if (length(sample_files) != 2) {
  stop("The deposit should hold two PBMC files. Found ", length(sample_files),
       ". Pull the accession with ./sync.sh first.")
}
donors <- sub("^PBMC_NK and T cell_Donor_([0-9]+)_.*$", "\\1",
              basename(sample_files))

states <- do.call(rbind, lapply(seq_along(sample_files), function(index) {
  frame <- read.FCS(sample_files[index], truncate_max_range = FALSE,
                    which.lines = 100)
  cbind(donor = donors[index], ReadCompensationState(frame))
}))
Write(states, "compensation_state.csv")
print(states, row.names = FALSE)

# ---------------------------------------------------------------------------
# Part 2. The matrix, computed from the single stains.
# ---------------------------------------------------------------------------

Say("\nPart 2: the spillover matrix from the single stains")

cell_result <- ComputeOmip58Spillover(
  kDepositDir, "^Comp_Cells", file.path(kOutputDir, "match_cells.csv"),
  keep_unstained = kUnstainedCells
)
bead_result <- ComputeOmip58Spillover(
  kDepositDir, "^Comp_Beads", file.path(kOutputDir, "match_beads.csv"),
  keep_unstained = kUnstainedBeads
)

cell_result$gateable$control_set <- "cells"
bead_result$gateable$control_set <- "beads"
gateable <- rbind(cell_result$gateable, bead_result$gateable)
Write(gateable, "control_gateability.csv")
Say("  cell controls gateable: ", sum(cell_result$gateable$gateable), " of ",
    nrow(cell_result$gateable))
Say("  bead controls gateable: ", sum(bead_result$gateable$gateable), " of ",
    nrow(bead_result$gateable))
dropped <- gateable[!gateable$gateable, c("control_set", "stain", "channel")]
if (nrow(dropped) > 0) {
  Say("  dropped:")
  print(dropped, row.names = FALSE)
}

spillover <- cell_result$spillover
Write(cbind(detector = rownames(spillover), as.data.frame(spillover)),
      "spillover_cells.csv")
Write(cbind(detector = rownames(bead_result$spillover),
            as.data.frame(bead_result$spillover)), "spillover_beads.csv")

shared <- intersect(rownames(spillover), rownames(bead_result$spillover))
agreement <- data.frame(
  shared_detectors = length(shared),
  correlation = stats::cor(as.vector(spillover[shared, shared]),
                           as.vector(bead_result$spillover[shared, shared])),
  cells_only = paste(setdiff(rownames(spillover),
                             rownames(bead_result$spillover)), collapse = " "),
  beads_only = paste(setdiff(rownames(bead_result$spillover),
                             rownames(spillover)), collapse = " "),
  stringsAsFactors = FALSE
)
Write(agreement, "spillover_agreement.csv")
print(agreement, row.names = FALSE)

heatmap <- PlotSpilloverHeatmap(
  spillover, title = "OMIP-058 spillover, computed from the cell controls"
)
SaveFigure(heatmap, file.path(kOutputDir, "spillover_heatmap.svg"), width = 9,
  height = 8)

# ---------------------------------------------------------------------------
# Part 3. What compensation does to the marker correlations.
# ---------------------------------------------------------------------------

Say("\nPart 3: marker correlation before and after compensation")

# Spillover makes two detectors move together whatever the biology says. The
# median correlation over every marker pair therefore measures how much
# spillover is left, and it needs no gate and no threshold to compute.
CorrelationSummary <- function(frame, spillover_matrix, label, donor) {
  if (!is.null(spillover_matrix)) {
    frame <- compensate(frame, spillover_matrix)
  }
  channels <- ResolveOmip58Channels(frame)
  scatter <- PanelScatterChannels(frame)
  values <- ArcsinhTransform(exprs(frame), channels$channel)
  in_range <- InScatterRange(values, scatter)
  singlets <- in_range
  singlets[in_range] <- RatioSingletMask(values[in_range, , drop = FALSE],
                                          scatter)
  lymphocytes <- singlets
  lymphocytes[singlets] <- LymphocyteMask(values[singlets, , drop = FALSE],
                                          scatter)
  events <- values[lymphocytes, channels$channel, drop = FALSE]
  colnames(events) <- channels$name
  result <- MarkerCorrelation(events, seed = 42, threshold = 0.5)
  correlation <- result$matrix
  data.frame(
    donor = donor, compensation = label,
    median_absolute_r = result$summary$median_absolute_r,
    pairs_above_half = result$summary$positive_pairs_above_threshold,
    pairs = result$summary$pairs,
    cd4_against_cd8 = correlation["CD4", "CD8"],
    cd3_against_cd8 = correlation["CD3", "CD8"],
    cd16_against_viability = correlation["CD16", "viability"],
    stringsAsFactors = FALSE
  )
}

correlations <- list()
for (index in seq_along(sample_files)) {
  frame <- read.FCS(sample_files[index], truncate_max_range = FALSE)
  correlations[[length(correlations) + 1]] <-
    CorrelationSummary(frame, NULL, "as deposited", donors[index])
  correlations[[length(correlations) + 1]] <-
    CorrelationSummary(frame, bead_result$spillover, "bead matrix",
                       donors[index])
  correlations[[length(correlations) + 1]] <-
    CorrelationSummary(frame, spillover, "cell matrix", donors[index])
  rm(frame)
  invisible(gc())
}
correlations <- do.call(rbind, correlations)
Write(correlations, "marker_correlation.csv")
print(correlations, row.names = FALSE, digits = 3)

# ---------------------------------------------------------------------------
# Part 4. The gate, and the handoff to Python.
# ---------------------------------------------------------------------------

Say("\nPart 4: gating and the handoff")

all_counts <- list()
all_cuts <- list()
all_rare <- list()
manifest <- list()

for (index in seq_along(sample_files)) {
  donor <- donors[index]
  Say("  donor ", donor)
  gated <- GateOmip58File(sample_files[index], spillover)

  all_counts[[index]] <- cbind(donor = donor, gated$counts)
  all_cuts[[index]] <- cbind(donor = donor, gated$cuts)
  all_rare[[index]] <- cbind(donor = donor, gated$rare)
  print(gated$counts[, c("population", "events", "percent_of_parent")],
        row.names = FALSE, digits = 4)

  for (population in kOmip58Handoff) {
    name <- paste0("donor_", donor, "_", population, ".fcs")
    written <- WriteGatedPopulation(gated$frame, gated$masks[[population]],
                                    file.path(kHandoffDir, name))
    manifest[[length(manifest) + 1]] <- data.frame(
      file = name, donor = donor, population = population, events = written,
      source_file = basename(sample_files[index]),
      compensation = "cell single stains, applied by scripts/13",
      transform = "none, the Python side transforms",
      stringsAsFactors = FALSE
    )
    Say("    wrote ", name, ", ", written, " events")
  }
  rm(gated)
  invisible(gc())
}

Write(do.call(rbind, all_counts), "gate_counts.csv")
Write(do.call(rbind, all_cuts), "gate_cuts.csv")
Write(do.call(rbind, all_rare), "one_dimensional_cuts.csv")

manifest <- do.call(rbind, manifest)
utils::write.csv(manifest, file.path(kHandoffDir, "manifest.csv"),
                 row.names = FALSE)
print(manifest[, c("file", "donor", "population", "events")], row.names = FALSE)

Say("\nWhat a one dimensional cut does to the rare markers:")
print(do.call(rbind, all_rare), row.names = FALSE, digits = 4)

Say("\nDone. Run scripts/14_omip58_pytometry.py next.")
