#!/usr/bin/env Rscript

# Compensation and automated gating of the FlowJo Basic Tutorial PBMC dataset.
#
# The script runs the whole analysis and writes every intermediate result to
# output/. The Quarto report reads those files rather than repeating the work, so
# rendering the report is fast and the numbers in the report are the numbers this
# script produced.
#
# Run it in the container:
#   podman run --rm -v "$PWD:/work:z" -w /work everything-flow-cytometry:latest \
#     Rscript scripts/01_automated_gating_pbmc.R
#
# The dataset is 8 FCS files of about 250,000 events each, from two donors under
# four stimulation conditions, acquired on an LSRII on 18 January 2013.

suppressPackageStartupMessages({
  library(flowCore)
  library(flowWorkspace)
  library(openCyto)
  library(ggcyto)
  library(ggplot2)
})

for (file in list.files("R", pattern = "\\.R$", full.names = TRUE)) {
  source(file)
}

kDataDir <- file.path(
  "data", "datasets", "flowjo",
  "Basic_Tutorial+Data_20180102T225815Z-001",
  "TQC Basic Tutorial Data"
)
kTemplatePath <- file.path("gating", "pbmc_gating_template.csv")
kOutputDir <- file.path("output", "automated_gating_pbmc")
kBlankMarker <- "Blank"

dir.create(kOutputDir, recursive = TRUE, showWarnings = FALSE)

Log <- function(...) {
  cat(format(Sys.time(), "%H:%M:%S"), "|", ..., "\n")
}

# ---------------------------------------------------------------------------
# 1. Read the files and derive the sample sheet
# ---------------------------------------------------------------------------

Log("Reading FCS files from", kDataDir)
if (!dir.exists(kDataDir)) {
  stop(
    "The dataset is not present. Pull it with:\n",
    "  ./sync.sh pull datasets/flowjo"
  )
}

raw_set <- ReadTutorialFlowSet(kDataDir)
sample_sheet <- pData(raw_set)

Log("Read", length(raw_set), "samples")
print(sample_sheet[, c("donor", "condition", "well", "stimulated")])

write.csv(
  sample_sheet,
  file.path(kOutputDir, "sample_sheet.csv"),
  row.names = FALSE
)

channels <- DescribeChannels(raw_set[[1]])
write.csv(channels, file.path(kOutputDir, "panel.csv"), row.names = FALSE)
Log("Panel holds", sum(channels$is_marker), "named markers of",
    nrow(channels), "parameters")

event_counts <- data.frame(
  sample = sampleNames(raw_set),
  events = as.integer(fsApply(raw_set, nrow)),
  stringsAsFactors = FALSE
)
write.csv(event_counts, file.path(kOutputDir, "event_counts.csv"),
          row.names = FALSE)

# ---------------------------------------------------------------------------
# 2. Compensation
# ---------------------------------------------------------------------------

Log("Applying the spillover matrix stored in the files")
spillover <- ExtractSpillover(raw_set[[1]])
Log("Spillover matrix is", nrow(spillover), "by", ncol(spillover))

write.csv(spillover, file.path(kOutputDir, "spillover_matrix.csv"))

spillover_summary <- SummariseSpillover(spillover, top = 10)
write.csv(spillover_summary, file.path(kOutputDir, "spillover_top10.csv"),
          row.names = FALSE)
Log("Largest spillover:", spillover_summary$from[1], "into",
    spillover_summary$to[1], "at",
    sprintf("%.1f%%", 100 * spillover_summary$spill[1]))

compensated_set <- ApplyCompensation(raw_set)

# ---------------------------------------------------------------------------
# 3. Transformation
# ---------------------------------------------------------------------------

Log("Estimating one logicle transform on the first sample")
blank_channel <- tryCatch(
  ChannelForMarker(raw_set[[1]], kBlankMarker),
  error = function(e) character()
)

transform_result <- ApplyLogicleTransform(
  compensated_set,
  channels = FluorescenceChannels(raw_set[[1]], exclude = blank_channel)
)
transformed_set <- transform_result$data
Log("Transformed", length(transform_result$channels), "fluorescence channels")

writeLines(
  transform_result$channels,
  file.path(kOutputDir, "transformed_channels.txt")
)

# ---------------------------------------------------------------------------
# 4. Automated gating from the template
# ---------------------------------------------------------------------------

Log("Reading the gating template", kTemplatePath)
template <- ReadGatingTemplate(kTemplatePath)

Log("Running automated gating. Every gate is fitted per sample.")
gating_set <- RunAutomatedGating(transformed_set, template)

population_paths <- gs_get_pop_paths(gating_set, path = "auto")
Log("Gated", length(population_paths), "populations:",
    paste(population_paths, collapse = ", "))

writeLines(population_paths, file.path(kOutputDir, "population_paths.txt"))

save_gs(gating_set, file.path(kOutputDir, "gating_set"))
Log("Saved the GatingSet so the report can reload it without re-gating")

# ---------------------------------------------------------------------------
# 5. Statistics
# ---------------------------------------------------------------------------

Log("Collecting counts and frequencies")
population_stats <- CollectPopulationStats(gating_set, sample_sheet)
write.csv(population_stats, file.path(kOutputDir, "population_stats.csv"),
          row.names = FALSE)

spread_all <- SummarisePopulationSpread(population_stats)
write.csv(spread_all, file.path(kOutputDir, "population_spread.csv"),
          row.names = FALSE)

spread_by_condition <- SummarisePopulationSpread(
  population_stats,
  group_by = "condition"
)
write.csv(
  spread_by_condition,
  file.path(kOutputDir, "population_spread_by_condition.csv"),
  row.names = FALSE
)

Log("Population spread across all eight samples:")
print(spread_all)

# ---------------------------------------------------------------------------
# 6. Figures
# ---------------------------------------------------------------------------

Log("Writing figures")

SavePlot <- function(plot_object, name, width = 9, height = 6) {
  path <- file.path(kOutputDir, paste0(name, ".png"))
  ggsave(path, plot = plot_object, width = width, height = height, dpi = 150)
  Log("  wrote", basename(path))
}

# The gate hierarchy as a tree.
png(file.path(kOutputDir, "gating_tree.png"), width = 900, height = 700)
plot(gating_set)
dev.off()
Log("  wrote gating_tree.png")

# Every gate of the hierarchy on the first sample.
for (population in setdiff(population_paths, "root")) {
  gate_plot <- tryCatch(
    autoplot(gating_set[[1]], population, bins = 128),
    error = function(e) NULL
  )
  if (!is.null(gate_plot)) {
    path <- file.path(
      kOutputDir,
      paste0("gate_", gsub("[^A-Za-z0-9]+", "_", population), ".png")
    )
    png(path, width = 700, height = 600)
    print(gate_plot)
    dev.off()
    Log("  wrote", basename(path))
  }
}

# The frequency of every population, by condition.
frequency_plot <- ggplot(
  population_stats,
  aes(x = condition, y = percent_of_parent, colour = donor)
) +
  geom_point(size = 3, alpha = 0.8) +
  facet_wrap(~population, scales = "free_y") +
  labs(
    title = "Population frequency, one automated template across eight samples",
    x = "Stimulation condition",
    y = "Percent of parent population",
    colour = "Donor"
  ) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
SavePlot(frequency_plot, "population_frequencies", width = 10, height = 7)

Log("Done. Output is in", kOutputDir)
