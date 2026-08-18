# OMIP-051, a 28 colour B cell and myeloid panel.
#
# FR-FCM-ZYN4 holds one PBMC file and 59 compensation controls. This script
# computes the spillover matrix from the single stains, measures what applying
# it does to the marker correlations, gates the hierarchy of Figure 1, and
# judges the paper's claims against the result.
#
# The deposit comes from the same laboratory and the same instrument as
# FR-FCM-ZYRN, eight days earlier, and it stores the same identity spillover
# matrix. Whether the values carry their spillover is measured here rather than
# read off the keyword.
#
# Run it in the container:
#   podman run --rm -v "$PWD:/work:z" -w /work everything-flow-cytometry:latest \
#     Rscript scripts/15_omip51_bcell_dc.R

suppressPackageStartupMessages({
  library(flowCore)
  library(flowStats)
  library(ggplot2)
  library(robustbase)
  library(withr)
  library(patchwork)
  library(FlowSOM)
  library(uwot)
})

for (module in c("figures", "io", "compensation", "spillover_compute",
                 "naive_memory", "panels", "clustering", "omip51")) {
  source(file.path("R", paste0(module, ".R")))
}

# covMcd and the scatter filter inside spillover_ng both draw random subsets, so
# without a seed the matrix moves between runs and a count in the report cannot
# be reproduced.
kSeed <- 42
set.seed(kSeed)

kDepositDir <- file.path("data", "datasets", "flowrepository", "OMIP-51",
                         "FlowRepository_FR-FCM-ZYN4_files")
kOutputDir <- file.path("output", "omip51")
kGatingDir <- "gating"
kUnstainedCells <- "Comp_Cells_unstained_F8_F08_061.fcs"
kUnstainedBeads <- "Comp_Beads_unstained mouse_F4_F04_030.fcs"

dir.create(kOutputDir, recursive = TRUE, showWarnings = FALSE)
Write <- function(x, name) {
  utils::write.csv(x, file.path(kOutputDir, name), row.names = FALSE)
  invisible(x)
}
Say <- function(...) cat(..., "\n", sep = "")

claims <- utils::read.csv(file.path(kGatingDir, "omip51_paper_claims.csv"),
                          stringsAsFactors = FALSE)

# ---------------------------------------------------------------------------
# Part 1. The panel, and the two markers a name match cannot resolve.
# ---------------------------------------------------------------------------

Say("Part 1: the panel")

sample_files <- list.files(kDepositDir, pattern = "^PBMC", full.names = TRUE)
if (length(sample_files) != 1) {
  stop("The deposit should hold one PBMC file. Found ", length(sample_files),
       ". Pull the accession with ./sync.sh first.")
}

header <- read.FCS(sample_files[1], truncate_max_range = FALSE,
                   which.lines = 100)
channels <- ResolveOmip51Channels(header)
Write(channels, "channel_resolution.csv")
print(channels, row.names = FALSE)

ambiguous <- DescribeChannels(header)
ambiguous <- ambiguous[grepl(" or ", ambiguous$marker), ]
Write(ambiguous[, c("channel", "marker")], "ambiguous_markers.csv")
Say("\n  detectors whose marker name lists two antibodies: ", nrow(ambiguous))
print(ambiguous[, c("channel", "marker")], row.names = FALSE)

state <- cbind(file = basename(sample_files[1]), ReadCompensationState(header))
Write(state, "compensation_state.csv")
Say("")
print(state, row.names = FALSE)

# ---------------------------------------------------------------------------
# Part 2. The matrix, computed from the single stains.
# ---------------------------------------------------------------------------

Say("\nPart 2: the spillover matrix from the single stains")

cell_result <- ComputeOmip51Spillover(
  kDepositDir, "^Comp_Cells|^Comp_Live", file.path(kOutputDir, "match_cells.csv"),
  keep_unstained = kUnstainedCells
)
bead_result <- ComputeOmip51Spillover(
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
dropped <- gateable[!gateable$gateable, c("control_set", "stain", "channel",
                                          "gated_events")]
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

spillover_figure <- PlotSpilloverHeatmap(
  spillover, title = "OMIP-051 spillover, from the cell controls"
)
SaveFigure(spillover_figure, file.path(kOutputDir, "spillover_heatmap.svg"),
           width = 9, height = 8)

# ---------------------------------------------------------------------------
# Part 3. What compensation does to the marker correlations.
# ---------------------------------------------------------------------------

Say("\nPart 3: marker correlation before and after compensation")

# Two detectors that carry spillover move together whatever the biology says, so
# the median correlation over every marker pair measures how much spillover is
# left without needing a gate or a threshold.
CorrelationSummary <- function(frame, matrix_used, label) {
  if (!is.null(matrix_used)) {
    frame <- compensate(frame, matrix_used)
  }
  panel <- ResolveOmip51Channels(frame)
  scatter <- PanelScatterChannels(frame)
  values <- ArcsinhTransform(exprs(frame), panel$channel)
  in_range <- InScatterRange(values, scatter)
  singlets <- in_range
  singlets[in_range] <- RatioSingletMask(values[in_range, , drop = FALSE],
                                         scatter)
  events <- values[singlets, panel$channel, drop = FALSE]
  colnames(events) <- panel$name
  rows <- withr::with_seed(kSeed, sample(nrow(events),
                                         min(100000, nrow(events))))
  correlation <- stats::cor(events[rows, ])
  pairs <- data.frame(
    a = rep(rownames(correlation), times = ncol(correlation)),
    b = rep(colnames(correlation), each = nrow(correlation)),
    r = as.vector(correlation), stringsAsFactors = FALSE
  )
  pairs <- pairs[pairs$a < pairs$b, ]
  data.frame(
    compensation = label,
    median_absolute_r = stats::median(abs(pairs$r)),
    pairs_above_half = sum(pairs$r > 0.5), pairs = nrow(pairs),
    cd19_against_cd14 = correlation["CD19", "CD14"],
    cd20_against_cd19 = correlation["CD20", "CD19"],
    igd_against_cd27 = correlation["IgD", "CD27"],
    stringsAsFactors = FALSE
  )
}

full_frame <- read.FCS(sample_files[1], truncate_max_range = FALSE)
correlations <- rbind(
  CorrelationSummary(full_frame, NULL, "as deposited"),
  CorrelationSummary(full_frame, bead_result$spillover, "bead matrix"),
  CorrelationSummary(full_frame, spillover, "cell matrix")
)
Write(correlations, "marker_correlation.csv")
print(correlations, row.names = FALSE, digits = 3)
rm(full_frame)
invisible(gc())

# ---------------------------------------------------------------------------
# Part 4. The gate hierarchy of Figure 1.
# ---------------------------------------------------------------------------

Say("\nPart 4: gating")

# Every lineage threshold comes from the deposit's unstained cell control at its
# 99.9th percentile. A fitted density cut is recorded beside each one and is not
# used, because it fails on CD14 and on the viability channel of this file.
kUnstainedControl <- file.path(kDepositDir, kUnstainedCells)
thresholds <- UnstainedThresholds(kUnstainedControl, channels, spillover)
Write(data.frame(marker = names(thresholds), threshold = unname(thresholds)),
      "unstained_thresholds.csv")
Say("  thresholds from ", basename(kUnstainedControl), ":")
print(data.frame(marker = names(thresholds),
                 threshold = round(unname(thresholds), 3)), row.names = FALSE)

gated <- GateOmip51File(sample_files[1], spillover, thresholds = thresholds)
Write(gated$counts, "gate_counts.csv")
Write(gated$cuts, "gate_cuts.csv")
Write(gated$rare, "one_dimensional_cuts.csv")
print(gated$counts[, c("population", "parent", "events", "percent_of_parent")],
      row.names = FALSE, digits = 4)
Say("")
print(gated$cuts[, c("label", "marker", "cut", "density_cut", "mixture_cut",
                     "parent_events")], row.names = FALSE, digits = 4)

counts <- gated$counts
Events <- function(name) counts$events[counts$population == name]
Percent <- function(name) counts$percent_of_parent[counts$population == name]

scatter_removed <- 100 * (Events("b_cells") - Events("b_lymphocytes")) /
  max(Events("b_cells"), 1)
Say("\n  high scatter events inside the B cell gate: ",
    round(scatter_removed, 2), " percent")

SaveFigure(
  PlotGateTree(counts, title = "OMIP-051 gate hierarchy, one PBMC file"),
  file.path(kOutputDir, "gate_tree.svg"), width = 11, height = 8)

# The gating flow as Figure 1A draws it, one panel per boundary, with the
# threshold that was applied.
transformed <- ArcsinhTransform(exprs(gated$frame), channels$channel,
                                gated$cofactor)
scatter <- PanelScatterChannels(gated$frame)
Ch <- function(name) channels$channel[channels$name == name]
Threshold <- function(name) {
  row <- gated$cuts[gated$cuts$marker == name, ]
  if (nrow(row) == 0) NA_real_ else row$cut[1]
}
masks <- gated$masks

panels <- list(
  PlotGatePair(transformed, masks[["in_scatter_range"]],
               scatter[["forward_area"]], scatter[["forward_height"]],
               "FSC-A", "FSC-H", title = "Singlets"),
  PlotGatePair(transformed, masks[["singlets"]], scatter[["forward_area"]],
               scatter[["side_area"]], "FSC-A", "SSC-A", title = "Cells"),
  PlotGatePair(transformed, masks[["cells"]], scatter[["forward_area"]],
               Ch("viability"), "FSC-A", "Live Dead UV Blue",
               y_threshold = Threshold("viability"), title = "Viable"),
  PlotGatePair(transformed, masks[["viable"]], Ch("CD14"),
               scatter[["side_area"]], "CD14", "SSC-A",
               x_threshold = Threshold("CD14"), title = "Monocytes"),
  PlotGatePair(transformed, masks[["cd14_negative"]], Ch("HLADR"), Ch("CD20"),
               "HLA-DR", "CD20", x_threshold = Threshold("HLADR"),
               y_threshold = Threshold("CD20"),
               title = "B and myeloid cells"),
  PlotGatePair(transformed, masks[["b_and_myeloid"]], Ch("CD19"), Ch("CD20"),
               "CD19", "CD20", x_threshold = Threshold("CD19"),
               y_threshold = Threshold("CD20"), title = "B cells"),
  PlotGatePair(transformed, masks[["b_cells"]], scatter[["side_area"]],
               scatter[["forward_area"]], "SSC-A", "FSC-A",
               title = "Lymphocyte gate, after the B cell gate"),
  PlotGatePair(transformed, masks[["dendritic_parent"]], Ch("CD11c"),
               Ch("CD123"), "CD11c", "CD123",
               x_threshold = Threshold("CD11c"),
               y_threshold = Threshold("CD123"), title = "Dendritic cells")
)
gating_flow <- patchwork::wrap_plots(panels, ncol = 2) +
  patchwork::plot_annotation(tag_levels = "A",
                             theme = ThemePublication())
SaveFigure(gating_flow, file.path(kOutputDir, "gating_flow.png"),
           width = 9.5, height = 16)
Say("  wrote gate_tree.png and gating_flow.png")

# ---------------------------------------------------------------------------
# Part 5. Clustering, a second route to the same subsets.
# ---------------------------------------------------------------------------

Say("\nPart 5: clustering")

# The gate reads one marker at a time and needs a threshold for each. Clustering
# reads every marker at once and needs none, so it is a second route to the
# subsets whose thresholds the unstained control had to supply.
ClusterRoute <- function(mask, markers, definitions_path, label, title,
                         n_metaclusters = 12, subsample = 50000) {
  detectors <- channels$channel[channels$name %in% markers]
  events <- transformed[mask, detectors, drop = FALSE]
  colnames(events) <- channels$name[match(detectors, channels$channel)]
  sampled <- SubsampleEvents(events, n = subsample, seed = kSeed)

  clustering <- RunFlowSomClustering(sampled, colnames(sampled),
                                     n_metaclusters = n_metaclusters,
                                     scale_channels = TRUE, seed = kSeed)
  medians <- ClusterMedianExpression(sampled, clustering$metacluster,
                                     colnames(sampled))
  definitions <- ReadCellTypeDefinitions(definitions_path)
  annotation <- AnnotateClusters(medians, definitions)

  embedding <- RunUmapEmbedding(sampled, colnames(sampled), seed = kSeed)
  embedding$cell_type <- annotation$cell_type[
    match(clustering$metacluster, annotation$cluster)
  ]

  figure <- ggplot(embedding, aes(x = umap_1, y = umap_2, colour = cell_type)) +
    geom_point(size = 0.35, alpha = 0.6, shape = 16) +
    ScaleColourPublication(name = NULL) +
    LegendPoints() +
    labs(
      title = title,
      subtitle = paste0(
        CountLabels(nrow(sampled)), " events, ", n_metaclusters,
        " FlowSOM metaclusters, each labelled by its median expression"
      ),
      x = "UMAP 1", y = "UMAP 2"
    ) +
    ThemeEmbedding(base_size = 11)
  SaveFigure(figure, file.path(kOutputDir, paste0("umap_", label, ".png")),
    width = 9, height = 6.5)

  list(medians = medians, annotation = annotation,
       summary = SummariseCellTypes(annotation))
}

b_markers <- c("CD10", "CD20", "CD27", "IgD", "IgM", "IgG", "IgA", "CD21",
               "CD85j")
b_route <- ClusterRoute(gated$masks[["b_lymphocytes"]], b_markers,
                        file.path(kGatingDir, "omip51_bcell_definitions.csv"),
                        "b_cells", "B cell subsets inside the CD19 CD20 gate")
Write(cbind(cluster = rownames(b_route$medians),
            as.data.frame(b_route$medians)), "bcell_cluster_medians.csv")
Write(b_route$annotation, "bcell_cluster_annotation.csv")
Write(b_route$summary, "bcell_cluster_summary.csv")
Say("  B cell clusters:")
print(b_route$annotation[, c("cluster", "events", "percent_of_total",
                             "cell_type", "margin")],
      row.names = FALSE, digits = 3)

myeloid_markers <- c("CD19", "CD20", "HLADR", "CD123", "CD11c", "CD1c",
                     "CD141", "CD16")
myeloid_route <- ClusterRoute(
  gated$masks[["dendritic_parent"]], myeloid_markers,
  file.path(kGatingDir, "omip51_myeloid_definitions.csv"), "myeloid",
  "Myeloid subsets inside the HLA-DR positive CD19 CD20 negative gate",
  n_metaclusters = 10
)
Write(cbind(cluster = rownames(myeloid_route$medians),
            as.data.frame(myeloid_route$medians)), "myeloid_cluster_medians.csv")
Write(myeloid_route$annotation, "myeloid_cluster_annotation.csv")
Write(myeloid_route$summary, "myeloid_cluster_summary.csv")
Say("\n  myeloid clusters:")
print(myeloid_route$annotation[, c("cluster", "events", "percent_of_total",
                                   "cell_type", "margin")],
      row.names = FALSE, digits = 3)

# The two routes side by side, where a population carries a name in both.
gate_percent <- function(name) {
  value <- counts$percent_of_parent[counts$population == name]
  if (length(value) == 0) NA_real_ else value
}
route_map <- data.frame(
  cell_type = c("Transitional B cells", "Naive B cells",
                "Marginal zone B cells", "IgD only memory B cells",
                "Plasmablasts"),
  gate_population = c("transitional_b", "naive_b", "marginal_zone_b",
                      "igd_only_memory_b", "plasmablasts"),
  stringsAsFactors = FALSE
)
route_map$clustering_percent <- b_route$summary$percent_of_total[
  match(route_map$cell_type, b_route$summary$cell_type)
]
route_map$gate_percent <- vapply(route_map$gate_population, gate_percent,
                                 numeric(1))
route_map$clustering_percent[is.na(route_map$clustering_percent)] <- 0
Write(route_map, "route_comparison.csv")
Say("\n  the two routes, as a percentage of their own parent:")
print(route_map, row.names = FALSE, digits = 3)

# ---------------------------------------------------------------------------
# Part 6. The claims.
# ---------------------------------------------------------------------------

Say("\nPart 6: the claims")

Percent <- function(name) counts$percent_of_parent[counts$population == name]
Events <- function(name) counts$events[counts$population == name]
PercentOfBCells <- function(name) {
  100 * Events(name) / max(Events("b_lymphocytes"), 1)
}
Verdict <- function(passed) {
  if (is.na(passed)) "unresolved" else if (passed) "reproduced" else
    "not reproduced"
}

rows <- list()
Add <- function(id, observed, passed) {
  rows[[length(rows) + 1]] <<- data.frame(
    claim_id = id, observed = observed, verdict = Verdict(passed),
    stringsAsFactors = FALSE
  )
}

mz <- PercentOfBCells("marginal_zone_b")
Add(1, paste0("marginal zone B cells are ", round(mz, 2),
              " percent of B cells"), mz <= 20)

plasmablasts <- Events("plasmablasts")
clustered_plasmablasts <- b_route$summary$events[
  b_route$summary$cell_type == "Plasmablasts"
]
if (length(clustered_plasmablasts) == 0) clustered_plasmablasts <- 0
Add(2, paste0("the threshold gate holds ", plasmablasts,
              " events and the clustering labels one cluster of ",
              clustered_plasmablasts), clustered_plasmablasts > 0)

# The threshold comes from an unstained control and not from a fluorescence
# minus one control, so it carries none of the spillover the other 27 dyes put
# into a channel. Where the frequency it produces is far from what the paper
# describes, the claim is unresolved rather than refuted.
transitional <- Percent("transitional_b")
clustered_transitional <- b_route$summary$percent_of_total[
  b_route$summary$cell_type == "Transitional B cells"
]
if (length(clustered_transitional) == 0) clustered_transitional <- 0
Add(3, paste0("the CD10 threshold selects ", round(transitional, 2),
              " percent of B cells and the clustering labels ",
              round(clustered_transitional, 2), " percent"),
    clustered_transitional > 0)

Add(4, paste0("naive cells are ", round(Percent("naive_b"),
              2), " percent, IgD positive CD27 positive cells ",
              round(Percent("igd_cd27_b"), 2),
              " percent and IgD negative cells ",
              round(Percent("igd_negative_b"), 2),
              " percent of mature B cells"),
    Percent("naive_b") > 0 && Percent("igd_cd27_b") > 0 &&
      Percent("igd_negative_b") > 0)

Add(5, paste0("monocytes are ", round(Percent("monocytes"), 2),
              " percent of viable events"), Percent("monocytes") > 0)
Add(6, paste0("CD14 negative cells carrying HLA-DR or CD20 are ",
              round(Percent("b_and_myeloid"), 2),
              " percent of CD14 negative events"),
    Percent("b_and_myeloid") > 0)
Add(7, paste0("B cells are ", round(Percent("b_cells"), 2),
              " percent of the preselected population"),
    Percent("b_cells") > 0)
Add(8, paste0(round(100 - Percent("b_lymphocytes"), 2),
              " percent of the B cell gate is removed by the lymphocyte gate"),
    (100 - Percent("b_lymphocytes")) > 0)
Add(9, paste0("plasmacytoid cells are ",
              round(Percent("plasmacytoid_dendritic_cells"), 2),
              " percent and myeloid cells ",
              round(Percent("myeloid_dendritic_cells"), 2),
              " percent of the dendritic cell parent"),
    Percent("plasmacytoid_dendritic_cells") > 0 &&
      Percent("myeloid_dendritic_cells") > 0)

# CD1c and CD141 are two of the markers whose fitted cut selects most of the
# parent, so the three subsets cannot be separated on this file.
mdc_rows <- gated$rare[gated$rare$parent == "myeloid_dendritic_cells", ]
myeloid_found <- myeloid_route$summary$cell_type
subsets_found <- sum(c("CD1c positive myeloid dendritic cells",
                       "CD141 positive myeloid dendritic cells",
                       "Double negative myeloid dendritic cells") %in%
                       myeloid_found)
Add(10, paste0("a fitted cut selects ",
               round(mdc_rows$percent_selected[mdc_rows$marker == "CD1c"], 1),
               " percent of myeloid dendritic cells on CD1c and ",
               round(mdc_rows$percent_selected[mdc_rows$marker == "CD141"], 1),
               " percent on CD141, and the clustering labels ", subsets_found,
               " of the three subsets"), subsets_found == 3)

Add(11, paste0(nrow(channels), " of the 28 markers were resolved to a detector"),
    NA)

verdicts <- merge(claims, do.call(rbind, rows), by = "claim_id")
verdicts <- verdicts[order(verdicts$claim_id), ]
Write(verdicts[, c("claim_id", "short_name", "expected", "observed",
                   "verdict")], "claim_verdicts.csv")
print(verdicts[, c("claim_id", "short_name", "verdict")], row.names = FALSE)

Write(as.data.frame(table(verdict = verdicts$verdict)), "verdict_counts.csv")
Say("")
print(table(verdicts$verdict))

Say("\nDone. Tables and figures are in ", kOutputDir, ".")
