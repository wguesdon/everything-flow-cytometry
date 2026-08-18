#!/usr/bin/env Rscript

# A third route to the OMIP-039 populations: cluster the events instead of
# gating
# them, and draw a UMAP.
#
# Debris, dead cells and doublets are removed by gating first, and the
# clustering
# runs on the lymphocyte gate only. A cluster of debris is not a finding, and it
# stretches the embedding so that everything else is squeezed into a corner.
#
# The clustering never sees the gate hierarchy below that point. It is given the
# lineage channels and a number of metaclusters, and each metacluster is then
# labelled by scoring its median marker expression against
# gating/omip39_cell_type_definitions.csv. No cluster is named by eye.
#
# The output is a three way comparison of the frequency of each population:
# manual gating by the authors, automated gating by the template, and
# clustering.
#
# Run it in the container:
#   podman run --rm -v "$PWD:/work:z" -w /work everything-flow-cytometry:latest
# \
#     Rscript scripts/05_omip39_clustering_umap.R

suppressPackageStartupMessages({
  library(flowCore)
  library(flowStats)
  library(flowWorkspace)
  library(openCyto)
  library(CytoML)
  library(FlowSOM)
  library(uwot)
  library(ggplot2)
  library(withr)
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
kDefinitionsPath <- file.path("gating", "omip39_cell_type_definitions.csv")
kOutputDir <- file.path("output", "omip39")

# Cluster on lineage markers only. NKp30, CD7 and the activation related markers
# are left out, so a cluster is defined by identity rather than by state.
kLineageChannels <- c(
  "CD3 PECy5",
  "CD56 PEDazzle594",
  "NKG2C PE",
  "NKG2A PE-Vio770",
  "CD57 pure + aIgM BV605",
  "ILT2 APC",
  "Siglec-7 APCVio770",
  "CD2 PerCPCy55"
)

kPreGate <- "Lymphocytes"
kSubsampleSize <- 50000
kMetaclusters <- 12
kSeed <- 42

dir.create(kOutputDir, recursive = TRUE, showWarnings = FALSE)

Log <- function(...) {
  cat(format(Sys.time(), "%H:%M:%S"), "|", ..., "\n")
}

if (!dir.exists(kDataDir)) {
  stop("OMIP-39 is not present. Pull it with:\n",
       "  ./sync.sh pull datasets/flowrepository/OMIP-39")
}

# ---------------------------------------------------------------------------
# 1. Pre-gate, so the clustering sees lymphocytes and not debris
# ---------------------------------------------------------------------------

Log("Computing compensation from the single stained controls")
control_set <- read.flowSet(
  list.files(kDataDir, pattern = "^Single stainings.*\\.fcs$",
             full.names = TRUE),
  truncate_max_range = FALSE
)
match_table <- MatchControlsToChannels(control_set,
                                       unstained_pattern = "unstained")
match_file <- file.path(kOutputDir, "spillover_match.csv")
WriteMatchFile(match_table, match_file)
computed_spillover <- ComputeSpilloverFromControls(
  control_set, match_file = match_file, method = "median", pregate = TRUE
)

sample_set <- read.flowSet(
  list.files(kDataDir, pattern = "^Samples_.*\\.fcs$", full.names = TRUE),
  truncate_max_range = FALSE
)
compensated_set <- compensate(sample_set, computed_spillover)
transformed_set <- ApplyLogicleTransform(compensated_set)$data

Log("Gating down to", kPreGate)
template <- ReadGatingTemplate(kTemplatePath)
gating_set <- RunAutomatedGating(transformed_set, template)

events <- ExtractGatedEvents(gating_set, kPreGate)
Log("Lymphocyte gate holds", nrow(events), "events")

# The event matrix carries detector names. Every file downstream, including the
# cell type definitions, names markers. Translate once, here.
events <- RenameChannelsToMarkers(events, sample_set[[1]])

subsample <- SubsampleEvents(events, n = kSubsampleSize, seed = kSeed)
Log("Clustering on", nrow(subsample), "events, sampled from",
    attr(subsample, "sampled_from"))

# ---------------------------------------------------------------------------
# 2. Cluster and embed
# ---------------------------------------------------------------------------

Log("Running FlowSOM with", kMetaclusters, "metaclusters")
clustering <- RunFlowSomClustering(
  subsample,
  channels = kLineageChannels,
  grid_size = 10,
  n_metaclusters = kMetaclusters,
  seed = kSeed
)

Log("Running UMAP")
embedding <- RunUmapEmbedding(subsample, channels = kLineageChannels,
                              seed = kSeed)

Log("Summarising median expression per cluster")
median_expression <- ClusterMedianExpression(
  subsample, clustering$metacluster, kLineageChannels
)
write.csv(median_expression,
          file.path(kOutputDir, "cluster_median_expression.csv"),
          row.names = FALSE)

# ---------------------------------------------------------------------------
# 3. Label the clusters from the definitions file
# ---------------------------------------------------------------------------

Log("Annotating clusters against", kDefinitionsPath)
definitions <- ReadCellTypeDefinitions(kDefinitionsPath)
annotation <- AnnotateClusters(median_expression, definitions)
write.csv(annotation, file.path(kOutputDir, "cluster_annotation.csv"),
          row.names = FALSE)

cat("\n=== Cluster annotation ===\n")
print(annotation[, c("cluster", "events", "percent_of_total", "cell_type",
                     "runner_up", "margin")], digits = 3)

# The score is a weighted mean on a 0 to 1 scale, so a margin of 0.1 is already
# a
# clear separation between the best and the second best label.
kCloseCallMargin <- 0.1
close_calls <- annotation[annotation$margin < kCloseCallMargin, ]
if (nrow(close_calls) > 0) {
  Log("Close calls, margin under", kCloseCallMargin, ":",
      paste(close_calls$cluster, collapse = ", "))
}

cell_types <- SummariseCellTypes(annotation)
write.csv(cell_types, file.path(kOutputDir, "clustering_cell_types.csv"),
          row.names = FALSE)

cat("\n=== Cell type frequencies from clustering ===\n")
print(cell_types, digits = 3)

# ---------------------------------------------------------------------------
# 4. Figures
# ---------------------------------------------------------------------------

plot_data <- cbind(
  embedding,
  cluster = factor(clustering$metacluster),
  cell_type = annotation$cell_type[
    match(clustering$metacluster, annotation$cluster)
  ]
)
# The CSV carries a draw of the embedding for the report to plot, and without a
# seed that draw changes on every run.
coordinate_rows <- withr::with_seed(
  kSeed, sample.int(nrow(plot_data), min(20000, nrow(plot_data)))
)
write.csv(
  plot_data[coordinate_rows, ],
  file.path(kOutputDir, "umap_coordinates.csv"),
  row.names = FALSE
)

Log("Writing figures")

umap_by_type <- ggplot(plot_data, aes(x = umap_1, y = umap_2,
                       colour = cell_type)) +
  geom_point(size = 0.35, alpha = 0.6, shape = 16) +
  ScaleColourPublication(name = "Cell type") +
  LegendPoints() +
  labs(
    title = "UMAP of the lymphocyte gate, coloured by the annotated cell type",
    subtitle = paste0(
      CountLabels(nrow(subsample)), " events, FlowSOM with ", kMetaclusters,
      " metaclusters, labelled from the definitions file"
    ),
    x = "UMAP 1", y = "UMAP 2"
  ) +
  ThemeEmbedding()
SaveFigure(umap_by_type, file.path(kOutputDir, "umap_cell_types.png"),
  width = 10, height = 7)

umap_by_cluster <- ggplot(plot_data, aes(x = umap_1, y = umap_2,
                          colour = cluster)) +
  geom_point(size = 0.35, alpha = 0.6, shape = 16) +
  ScaleColourPublication(name = "Metacluster") +
  LegendPoints() +
  labs(
    title = "The same embedding, coloured by metacluster",
    subtitle = "Several metaclusters can carry one cell type label",
    x = "UMAP 1", y = "UMAP 2"
  ) +
  ThemeEmbedding()
SaveFigure(umap_by_cluster, file.path(kOutputDir, "umap_metaclusters.png"),
  width = 10, height = 7)

# One panel per lineage marker, so the annotation can be checked against the
# data.
marker_long <- do.call(rbind, lapply(kLineageChannels, function(channel) {
  data.frame(
    umap_1 = embedding$umap_1,
    umap_2 = embedding$umap_2,
    marker = channel,
    value = subsample[, channel],
    stringsAsFactors = FALSE
  )
}))

umap_by_marker <- ggplot(marker_long, aes(x = umap_1, y = umap_2,
                         colour = value)) +
  geom_point(size = 0.25, alpha = 0.6, shape = 16) +
  facet_wrap(~marker, ncol = 4) +
  scale_colour_viridis_c(option = "viridis", name = "Intensity",
                         guide = ColourbarGuide()) +
  labs(
    title = "The same embedding, coloured by each lineage marker",
    subtitle = paste(
      "Read this against the annotation table. A label that does not match",
      "the marker panels is wrong."
    ),
    x = "UMAP 1", y = "UMAP 2"
  ) +
  ThemeEmbedding()
SaveFigure(umap_by_marker, file.path(kOutputDir, "umap_markers.png"),
  width = 14, height = 7)

expression_long <- do.call(rbind, lapply(kLineageChannels, function(channel) {
  data.frame(
    cluster = factor(median_expression$cluster),
    marker = channel,
    median = median_expression[[channel]],
    stringsAsFactors = FALSE
  )
}))
expression_long$scaled <- stats::ave(
  expression_long$median, expression_long$marker,
  FUN = function(x) {
    span <- max(x) - min(x)
    if (span == 0) rep(0.5, length(x)) else (x - min(x)) / span
  }
)

heatmap_plot <- ggplot(expression_long, aes(x = marker, y = cluster,
                       fill = scaled)) +
  geom_tile(colour = "white", linewidth = 0.3) +
  scale_fill_viridis_c(option = "viridis", name = "Scaled\nmedian",
                       guide = ColourbarGuide()) +
  labs(
    title = "Median marker expression per metacluster",
    subtitle = paste(
      "Scaled to run from 0 to 1 within each marker.",
      "This table drives the annotation."
    ),
    x = NULL, y = "Metacluster"
  ) +
  ThemePublication() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
SaveFigure(heatmap_plot, file.path(kOutputDir, "cluster_heatmap.svg"),
  width = 10, height = 6)

Log("Wrote four figures")

# ---------------------------------------------------------------------------
# 5. The three way comparison
# ---------------------------------------------------------------------------

Log("Collecting the gated frequencies for comparison")

automated_stats <- CollectPopulationStats(gating_set)
manual_gs <- ImportFlowJoGates(kWorkspace, kDataDir, group = "Sample")
manual_stats <- CollectPopulationStats(manual_gs)

# Every route is expressed as a percentage of the lymphocyte gate, so the three
# numbers are comparable. A percentage of a different parent is not.
PercentOfLymphocytes <- function(stats, ends_with, lymphocyte_path) {
  paths <- as.character(stats$population)
  lymphocyte_hit <- which(endsWith(paths, lymphocyte_path))
  if (length(lymphocyte_hit) == 0) return(NA_real_)
  lymphocyte_count <- stats$count[lymphocyte_hit[1]]

  hit <- which(endsWith(paths, ends_with))
  if (length(hit) == 0) return(NA_real_)

  100 * stats$count[hit[1]] / lymphocyte_count
}

# CD56bright is included even though the manual workspace has no gate for it.
# The
# clustering assigns events to it, so leaving the row out would hide where those
# events went and make the clustering look as though it lost them.
comparison <- data.frame(
  cell_type = c("T cells", "NKT-like cells", "CD56bright NK",
                "CD56dim NKG2C+ NK", "CD56dim NKG2C- NK", "Other lymphocytes"),
  manual_path = c("/CD56- CD3+", "/CD56+CD3+", NA,
                  "/CD56dim NKG2C+", "/CD56dim NKG2C-", "/CD56- CD3-"),
  automated_path = c("/CD56negCD3pos", "/CD56posCD3pos", "/CD56bright",
                     "/NKG2Cpos", "/NKG2Cneg", "/CD56negCD3neg"),
  stringsAsFactors = FALSE
)

SafePercent <- function(path, stats, lymphocyte_path) {
  if (is.na(path)) return(NA_real_)
  PercentOfLymphocytes(stats, path, lymphocyte_path)
}

comparison$manual_percent <- vapply(
  comparison$manual_path,
  SafePercent,
  numeric(1),
  stats = manual_stats,
  lymphocyte_path = "/Lymphocytes"
)
comparison$automated_percent <- vapply(
  comparison$automated_path,
  SafePercent,
  numeric(1),
  stats = automated_stats,
  lymphocyte_path = "/Lymphocytes"
)
comparison$clustering_percent <- cell_types$percent_of_total[
  match(comparison$cell_type, cell_types$cell_type)
]

# The three routes are compared again on the NK compartment as a whole. How the
# NK cells are split between bright and dim is a separate question from how many
# NK cells there are, and the two disagree very differently.
nk_rows <- comparison$cell_type %in%
  c("CD56bright NK", "CD56dim NKG2C+ NK", "CD56dim NKG2C- NK")
nk_total <- data.frame(
  route = c("manual", "automated", "clustering"),
  nk_percent = c(
    PercentOfLymphocytes(manual_stats, "/CD56+ CD3-", "/Lymphocytes"),
    PercentOfLymphocytes(automated_stats, "/CD56posCD3neg", "/Lymphocytes"),
    sum(comparison$clustering_percent[nk_rows], na.rm = TRUE)
  ),
  stringsAsFactors = FALSE
)
write.csv(nk_total, file.path(kOutputDir, "nk_compartment_total.csv"),
          row.names = FALSE)

cat("\n=== The NK compartment as a whole, percent of lymphocytes ===\n")
print(nk_total, digits = 3)

write.csv(comparison, file.path(kOutputDir, "three_way_comparison.csv"),
          row.names = FALSE)

cat("\n=== Three routes to the same populations, percent of lymphocytes ===\n")
print(comparison[, c("cell_type", "manual_percent", "automated_percent",
                     "clustering_percent")], digits = 3)

long_comparison <- do.call(rbind, lapply(
  c("manual_percent", "automated_percent", "clustering_percent"),
  function(column) {
    data.frame(
      cell_type = comparison$cell_type,
      route = sub("_percent", "", column),
      percent = comparison[[column]],
      stringsAsFactors = FALSE
    )
  }
))

route_plot <- ggplot(
  long_comparison[!is.na(long_comparison$percent), ],
  aes(x = cell_type, y = percent, fill = route)
) +
  geom_col(position = position_dodge(width = 0.8), width = 0.75) +
  coord_flip() +
  ScaleFillPublication(name = "Route") +
  labs(
    title = "Manual gating, automated gating and clustering on one file",
    subtitle = "Percent of the lymphocyte gate. All three read the same events.",
    x = NULL, y = "Percent of lymphocytes"
  ) +
  ThemePublication()
SaveFigure(route_plot, file.path(kOutputDir, "three_way_comparison.svg"),
  width = 10, height = 6)

Log("Done. Output is in", kOutputDir)
