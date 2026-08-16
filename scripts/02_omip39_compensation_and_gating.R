#!/usr/bin/env Rscript

# OMIP-39: computed compensation, automated gating, and a comparison against the
# published manual gates.
#
# This is the positive control for the repository. OMIP-39 is a published panel,
# and its FlowRepository deposit ships the FlowJo workspace that the authors used.
# The manual result therefore does not come from this project. The script gates
# the same file from a template and puts the two side by side.
#
# Three parts:
#   1. Compute a spillover matrix from the 12 single stained controls, and compare
#      it against the matrix the instrument stored in the sample file.
#   2. Gate the sample with an openCyto template that mirrors the published
#      hierarchy.
#   3. Compare every mapped population against the manual gates.
#
# Run it in the container:
#   podman run --rm -v "$PWD:/work:z" -w /work everything-flow-cytometry:latest \
#     Rscript scripts/02_omip39_compensation_and_gating.R
#
# Reference: Hammer Q, Romagnani C. OMIP-039: Detection and analysis of human
# adaptive NKG2C+ natural killer cells. Cytometry A 2017;997-1000.
# PMID 28715616. doi:10.1002/cyto.a.23168.

suppressPackageStartupMessages({
  library(flowCore)
  library(flowStats)
  library(flowWorkspace)
  library(openCyto)
  library(CytoML)
  library(ggcyto)
  library(ggplot2)
})

for (file in list.files("R", pattern = "\\.R$", full.names = TRUE)) {
  source(file)
}

kDataDir <- file.path(
  "data", "datasets", "flowrepository", "OMIP-39",
  "FlowRepository_FR-FCM-ZYY6_files"
)
kWorkspace <- file.path(kDataDir, "attachments", "OMIP_Hammer.wsp")
kTemplatePath <- file.path("gating", "omip39_gating_template.csv")
kPopulationMap <- file.path("gating", "omip39_population_map.csv")
kOutputDir <- file.path("output", "omip39")
kSamplePattern <- "^Samples_.*\\.fcs$"
kControlPattern <- "^Single stainings.*\\.fcs$"
kTolerancePoints <- 5

dir.create(kOutputDir, recursive = TRUE, showWarnings = FALSE)

Log <- function(...) {
  cat(format(Sys.time(), "%H:%M:%S"), "|", ..., "\n")
}

if (!dir.exists(kDataDir)) {
  stop(
    "OMIP-39 is not present. Pull it with:\n",
    "  ./sync.sh pull datasets/flowrepository/OMIP-39"
  )
}

# ---------------------------------------------------------------------------
# 1. Compensation computed from the single stained controls
# ---------------------------------------------------------------------------

Log("Reading the single stained controls")
control_files <- list.files(kDataDir, pattern = kControlPattern, full.names = TRUE)
Log("Found", length(control_files), "control files")

control_set <- read.flowSet(control_files, truncate_max_range = FALSE)

Log("Matching every control to the channel it stains")
match_table <- MatchControlsToChannels(control_set, unstained_pattern = "unstained")
print(match_table[, c("stain", "channel", "marker", "matched_by")], right = FALSE)

write.csv(match_table, file.path(kOutputDir, "control_match_table.csv"),
          row.names = FALSE)

matched <- sum(match_table$matched_by != "none")
Log("Matched", matched, "of", nrow(match_table), "controls")
Log("Match passes used:",
    paste(names(table(match_table$matched_by)), table(match_table$matched_by),
          sep = "=", collapse = ", "))

match_file <- file.path(kOutputDir, "spillover_match.csv")
WriteMatchFile(match_table, match_file)

Log("Checking whether each control can support a spillover estimate")
control_quality <- CheckControlQuality(control_set, match_table)
write.csv(control_quality, file.path(kOutputDir, "control_quality.csv"),
          row.names = FALSE)
print(control_quality[, c("stain", "positive_percent", "primary_is_brightest",
                          "verdict")], right = FALSE)
Log("Controls judged weak:",
    sum(control_quality$verdict == "weak"), "of", nrow(control_quality))

Log("Computing the spillover matrix from the controls")
# method = "median", not the flowStats default of "mode". These are cell
# controls, and several markers sit on a minority of cells, so the mode lands on
# the negative population. Measured against the stored matrix, mode gives a
# correlation of 0.619 and two spillover values above 100 percent, while median
# gives 0.975 and none. See ?ComputeSpilloverFromControls.
computed_spillover <- ComputeSpilloverFromControls(
  control_set,
  match_file = match_file,
  fsc = "FSC-A",
  ssc = "SSC-A",
  method = "median",
  pregate = TRUE
)
write.csv(computed_spillover,
          file.path(kOutputDir, "spillover_computed.csv"))
Log("Computed a", nrow(computed_spillover), "by", ncol(computed_spillover),
    "matrix")

computed_top <- SummariseSpillover(computed_spillover, top = 10)
write.csv(computed_top, file.path(kOutputDir, "spillover_computed_top10.csv"),
          row.names = FALSE)
Log("Largest computed spillover:", computed_top$from[1], "into",
    computed_top$to[1], "at", sprintf("%.1f%%", 100 * computed_top$spill[1]))

heatmap_plot <- PlotSpilloverHeatmap(
  computed_spillover,
  title = "OMIP-39 spillover, computed from 12 single stained controls"
)
ggsave(file.path(kOutputDir, "spillover_heatmap.png"), heatmap_plot,
       width = 10, height = 8, dpi = 150)
Log("Wrote spillover_heatmap.png")

# ---------------------------------------------------------------------------
# 2. Read the sample and compare the two matrices
# ---------------------------------------------------------------------------

Log("Reading the sample file")
sample_files <- list.files(kDataDir, pattern = kSamplePattern, full.names = TRUE)
sample_set <- read.flowSet(sample_files, truncate_max_range = FALSE)
Log("Read", length(sample_set), "sample(s):", sampleNames(sample_set))

stored_spillover <- tryCatch(
  ExtractSpillover(sample_set[[1]]),
  error = function(e) {
    Log("The sample carries no stored matrix:", conditionMessage(e))
    NULL
  }
)

if (!is.null(stored_spillover)) {
  write.csv(stored_spillover, file.path(kOutputDir, "spillover_stored.csv"))
  comparison_matrices <- CompareSpilloverMatrices(
    computed_spillover, stored_spillover, top = 20
  )
  write.csv(comparison_matrices,
            file.path(kOutputDir, "spillover_computed_vs_stored.csv"),
            row.names = FALSE)
  Log("Largest disagreement between the computed and the stored matrix:",
      sprintf("%.2f percentage points", max(abs(comparison_matrices$difference))))
} else {
  Log("Skipping the matrix comparison, because the sample stores no matrix.")
}

# The computed matrix is the one the analysis uses, because it was derived from
# the controls that were acquired with this sample.
Log("Applying the computed matrix to the sample")
compensated_set <- compensate(sample_set, computed_spillover)

# ---------------------------------------------------------------------------
# 3. Transform and gate from the template
# ---------------------------------------------------------------------------

Log("Estimating the logicle transform")
transform_result <- ApplyLogicleTransform(compensated_set)
transformed_set <- transform_result$data
Log("Transformed", length(transform_result$channels), "fluorescence channels")
writeLines(transform_result$channels,
           file.path(kOutputDir, "transformed_channels.txt"))

Log("Running the automated gating template")
template <- ReadGatingTemplate(kTemplatePath)
automated_gs <- RunAutomatedGating(transformed_set, template)

automated_paths <- gs_get_pop_paths(automated_gs, path = "auto")
Log("Automated hierarchy:", paste(automated_paths, collapse = ", "))
writeLines(automated_paths, file.path(kOutputDir, "automated_paths.txt"))

automated_stats <- CollectPopulationStats(automated_gs)
write.csv(automated_stats, file.path(kOutputDir, "automated_stats.csv"),
          row.names = FALSE)

png(file.path(kOutputDir, "automated_tree.png"), width = 900, height = 600)
plot(automated_gs)
dev.off()
Log("Wrote automated_tree.png")

# ---------------------------------------------------------------------------
# 4. The published manual gates
# ---------------------------------------------------------------------------

Log("Importing the published FlowJo workspace")
manual_gs <- ImportFlowJoGates(kWorkspace, kDataDir, group = "Sample")
manual_paths <- gs_get_pop_paths(manual_gs, path = "full")
Log("Manual hierarchy holds", length(manual_paths), "populations")
writeLines(manual_paths, file.path(kOutputDir, "manual_paths.txt"))

manual_stats <- CollectPopulationStats(manual_gs)
write.csv(manual_stats, file.path(kOutputDir, "manual_stats.csv"),
          row.names = FALSE)

png(file.path(kOutputDir, "manual_tree.png"), width = 1400, height = 900)
plot(manual_gs, boolean = FALSE)
dev.off()
Log("Wrote manual_tree.png")

# ---------------------------------------------------------------------------
# 5. The comparison
# ---------------------------------------------------------------------------

Log("Comparing the two results")
population_map <- ReadPopulationMap(kPopulationMap)
comparison <- CompareGatingResults(manual_stats, automated_stats, population_map)
comparison <- JudgeAgreement(comparison, tolerance_points = kTolerancePoints)

write.csv(comparison, file.path(kOutputDir, "comparison.csv"), row.names = FALSE)

print(comparison[, c("manual", "automated", "manual_percent",
                     "automated_percent", "difference", "verdict")])

Log("Agreement within", kTolerancePoints, "percentage points:",
    sum(comparison$verdict == "agree", na.rm = TRUE), "of", nrow(comparison))

comparison_plot <- ggplot(
  comparison[!is.na(comparison$difference), ],
  aes(x = manual_percent, y = automated_percent, colour = verdict)
) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey50") +
  geom_point(size = 4, alpha = 0.85) +
  geom_text(aes(label = automated), vjust = -1.1, size = 3, show.legend = FALSE) +
  scale_colour_manual(values = c(agree = "#2166ac", differ = "#b2182b")) +
  labs(
    title = "Automated template against published manual gates, OMIP-39",
    subtitle = paste0(
      "One donor file. The dashed line is exact agreement. Tolerance is ",
      kTolerancePoints, " percentage points."
    ),
    x = "Manual, percent of parent",
    y = "Automated, percent of parent",
    colour = "Verdict"
  ) +
  theme_bw()

ggsave(file.path(kOutputDir, "comparison_scatter.png"), comparison_plot,
       width = 9, height = 7, dpi = 150)
Log("Wrote comparison_scatter.png")

Log("Done. Output is in", kOutputDir)
