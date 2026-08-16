#!/usr/bin/env Rscript

# OMIP-043: antibody secreting cells across four tissues.
#
# The target population is defined as CD38 high and CD27 positive. That word
# "high" is the difficulty: CD38 high is the bright tail of one continuous
# distribution rather than a separated peak, and a density based gating method has
# no valley to find.
#
# Five parts:
#   1. The manual gates, which are the reference throughout.
#   2. The paper's counting claims. It states a 5 percent CV target and 400 to
#      2,000 events in the ASC gate, and both are arithmetic that can be checked.
#   3. Automated gating from a template that mirrors the manual hierarchy. It
#      fails, and the script records the failure rather than tuning until it
#      agrees.
#   4. Clustering, which recovers the population that gating missed.
#   5. The paper's phenotype claims about scatter, CD20, CD19 and HLA-DR.
#
# Every tissue group holds seven files that share one $SRC and one acquisition
# date, so they are replicate acquisitions of one preparation rather than seven
# donors. Donor variation is therefore absent by design, and the spread across
# them measures acquisition plus analysis.
#
# Run it in the container:
#   podman run --rm -v "$PWD:/work:z" -w /work everything-flow-cytometry:latest \
#     Rscript scripts/06_omip43_asc_analysis.R
#
# Reference: Carrell J, Groves CJ. OMIP-043: Identification of human antibody
# secreting cell subsets. Cytometry A 2018;93A:190-193. PMID 29286577.
# doi:10.1002/cyto.a.23305.

suppressPackageStartupMessages({
  library(flowCore)
  library(flowWorkspace)
  library(openCyto)
  library(CytoML)
  library(FlowSOM)
  library(ggplot2)
  library(withr)
})

for (file in list.files("R", pattern = "\\.R$", full.names = TRUE)) {
  source(file)
}

kDataDir <- file.path(
  "data", "datasets", "flowrepository", "OMIP-43",
  "FlowRepository_FR-FCM-ZYBP_files"
)
kWorkspace <- file.path(kDataDir, "attachments", "04-Aug-2017_OMIP.wsp")
kTemplatePath <- file.path("gating", "omip43_gating_template.csv")
kClaimsPath <- file.path("gating", "omip43_paper_claims.csv")
kDefinitionsPath <- file.path("gating", "omip43_cell_type_definitions.csv")
kOutputDir <- file.path("output", "omip43")

kTissues <- c("PBMC", "Bone Marrow", "Spleen", "Tonsil")
kParentPath <- "/time/live/scatter/sing1/sing2/dumped"
kAscPath <- paste0(kParentPath, "/PC")

# The paper's own targets.
kTargetCvPercent <- 5
kTargetEventsLow <- 400
kTargetEventsHigh <- 2000

# Channel names carry a Comp- prefix after CytoML applies the workspace matrix.
kChannels <- c(
  CD19 = "Comp-APC-A",
  CD20 = "Comp-BV421-A",
  IgD = "Comp-PerCP-Cy5-5-A",
  CD27 = "Comp-BV786-A",
  CD38 = "Comp-BUV395-A",
  Ig = "Comp-FITC-A",
  `HLA-DR` = "Comp-PE-A"
)
kClusterChannels <- unname(kChannels)
kMetaclusters <- 25
# How many of the CD38 highest clusters make up the ASC population. Chosen by the
# sweep in script 07, which ranks configurations on their worst tissue.
kAscClusters <- 2
kSubsampleSize <- 50000
kSeed <- 42

dir.create(kOutputDir, recursive = TRUE, showWarnings = FALSE)

Log <- function(...) {
  cat(format(Sys.time(), "%H:%M:%S"), "|", ..., "\n")
}

if (!dir.exists(kDataDir)) {
  stop("OMIP-43 is not present. Pull it with:\n",
       "  ./sync.sh pull datasets/flowrepository/OMIP-43")
}

claims <- ReadPaperClaims2 <- utils::read.csv(kClaimsPath, stringsAsFactors = FALSE,
                                              check.names = FALSE)
Log("Read", nrow(claims), "claims from the paper")

# ---------------------------------------------------------------------------
# 1. The manual gates
# ---------------------------------------------------------------------------

Log("Importing the deposited FlowJo workspace, one group per tissue")
workspace <- CytoML::open_flowjo_xml(kWorkspace)

gating_sets <- list()
for (tissue in kTissues) {
  gating_sets[[tissue]] <- CytoML::flowjo_to_gatingset(
    workspace, name = tissue, path = kDataDir
  )
  Log(" ", tissue, "imported,", length(gating_sets[[tissue]]), "replicates")
}

manual_rows <- list()
for (tissue in kTissues) {
  counts <- as.data.frame(
    flowWorkspace::gs_pop_get_stats(gating_sets[[tissue]], type = "count")
  )
  for (sample_name in unique(counts$sample)) {
    parent <- counts$count[counts$sample == sample_name &
                             counts$pop == kParentPath]
    asc <- counts$count[counts$sample == sample_name & counts$pop == kAscPath]
    root <- counts$count[counts$sample == sample_name & counts$pop == "root"]
    if (length(parent) == 0 || length(asc) == 0) next

    manual_rows[[length(manual_rows) + 1]] <- data.frame(
      tissue = tissue,
      sample = sample_name,
      total_events = root,
      parent_events = parent,
      asc_events = asc,
      asc_percent = 100 * asc / parent,
      poisson_cv_percent = PoissonCv(asc),
      stringsAsFactors = FALSE
    )
  }
}
manual <- do.call(rbind, manual_rows)
write.csv(manual, file.path(kOutputDir, "manual_asc_counts.csv"), row.names = FALSE)

cat("\n=== ASC per replicate, manual gates ===\n")
print(manual[, c("tissue", "asc_events", "parent_events", "asc_percent",
                 "poisson_cv_percent")], digits = 3, row.names = FALSE)

# ---------------------------------------------------------------------------
# 2. The paper's counting claims
# ---------------------------------------------------------------------------

Log("Testing the counting claims")

spread_rows <- lapply(kTissues, function(tissue) {
  piece <- manual[manual$tissue == tissue, ]
  out <- CompareSpreadToPoisson(piece$asc_events, piece$asc_percent)
  cbind(data.frame(tissue = tissue, stringsAsFactors = FALSE), out)
})
spread <- do.call(rbind, spread_rows)
write.csv(spread, file.path(kOutputDir, "asc_spread.csv"), row.names = FALSE)

cat("\n=== Spread across replicates, against the counting floor ===\n")
print(spread, digits = 3, row.names = FALSE)

counting_claims <- data.frame(
  tissue = kTissues,
  mean_events = spread$mean_count,
  in_target_range = vapply(kTissues, function(tissue) {
    piece <- manual[manual$tissue == tissue, ]
    all(piece$asc_events >= kTargetEventsLow &
          piece$asc_events <= kTargetEventsHigh)
  }, logical(1)),
  samples_in_range = vapply(kTissues, function(tissue) {
    piece <- manual[manual$tissue == tissue, ]
    sum(piece$asc_events >= kTargetEventsLow &
          piece$asc_events <= kTargetEventsHigh)
  }, integer(1)),
  poisson_cv_percent = spread$poisson_cv_percent,
  meets_cv_target = spread$poisson_cv_percent <= kTargetCvPercent,
  stringsAsFactors = FALSE
)
write.csv(counting_claims, file.path(kOutputDir, "counting_claims.csv"),
          row.names = FALSE)

cat("\n=== The paper's targets: 400 to 2000 events, 5 percent CV ===\n")
print(counting_claims, digits = 3, row.names = FALSE)

Log("Replicates inside the stated 400 to 2000 range:",
    sum(counting_claims$samples_in_range), "of", nrow(manual))

# ---------------------------------------------------------------------------
# 3. Automated gating, and its failure
# ---------------------------------------------------------------------------

Log("Running the automated template on one replicate per tissue")

automated_rows <- list()
for (tissue in kTissues) {
  files <- list.files(
    kDataDir,
    pattern = paste0("^", switch(tissue, PBMC = "PBMC", `Bone Marrow` = "BM",
                                 Spleen = "Spleen", Tonsil = "Tonsil"), "_.*\\.fcs$"),
    full.names = TRUE
  )
  if (length(files) == 0) {
    Log("  no raw file found for", tissue, ", skipping")
    next
  }

  flow_set <- read.flowSet(files[1], truncate_max_range = FALSE)
  spillover <- tryCatch(ExtractSpillover(flow_set[[1]]), error = function(e) NULL)
  if (!is.null(spillover)) flow_set <- compensate(flow_set, spillover)
  transformed <- ApplyLogicleTransform(flow_set)$data

  template <- ReadGatingTemplate(kTemplatePath)
  gating_set <- tryCatch(
    RunAutomatedGating(transformed, template),
    error = function(e) { Log("  gating failed for", tissue, ":",
                              conditionMessage(e)); NULL }
  )
  if (is.null(gating_set)) next

  counts <- as.data.frame(
    flowWorkspace::gs_pop_get_stats(gating_set, type = "count")
  )
  parent <- counts$count[grepl("/dumped$", counts$pop)]
  asc <- counts$count[grepl("/ASC$", counts$pop)]
  if (length(parent) == 0 || length(asc) == 0) next

  automated_rows[[length(automated_rows) + 1]] <- data.frame(
    tissue = tissue,
    file = basename(files[1]),
    parent_events = parent,
    asc_events = asc,
    asc_percent = 100 * asc / parent,
    stringsAsFactors = FALSE
  )
  Log("  ", tissue, sprintf("automated ASC = %.2f%%", 100 * asc / parent))
}

automated <- do.call(rbind, automated_rows)
manual_means <- aggregate(asc_percent ~ tissue, data = manual, FUN = mean)
automated <- merge(automated, manual_means, by = "tissue",
                   suffixes = c("_automated", "_manual"))
automated$fold_error <- automated$asc_percent_automated /
  automated$asc_percent_manual
write.csv(automated, file.path(kOutputDir, "automated_vs_manual.csv"),
          row.names = FALSE)

cat("\n=== Automated template against the manual mean ===\n")
print(automated[, c("tissue", "asc_percent_manual", "asc_percent_automated",
                    "fold_error")], digits = 3, row.names = FALSE)

# ---------------------------------------------------------------------------
# 4. Clustering, scored against the manual gate
# ---------------------------------------------------------------------------

Log("Clustering, with the manual gate as the reference")

cluster_rows <- list()
umap_saved <- FALSE

for (tissue in kTissues) {
  gating_set <- gating_sets[[tissue]]
  gating_hierarchy <- gating_set[[1]]

  in_parent <- flowWorkspace::gh_pop_get_indices(gating_hierarchy, kParentPath)
  in_asc <- flowWorkspace::gh_pop_get_indices(gating_hierarchy, kAscPath)
  truth_all <- in_asc[in_parent]

  events_all <- flowCore::exprs(
    flowWorkspace::gh_pop_get_data(gating_hierarchy, kParentPath)
  )

  withr::with_seed(kSeed, {
    keep <- sample.int(nrow(events_all), min(kSubsampleSize, nrow(events_all)))
  })
  events <- events_all[keep, , drop = FALSE]
  truth <- truth_all[keep]

  clustering <- RunFlowSomClustering(
    events, channels = kClusterChannels,
    grid_size = 10, n_metaclusters = kMetaclusters, seed = kSeed
  )

  # Purity and recall per metacluster, against the manual gate.
  per_cluster <- lapply(sort(unique(clustering$metacluster)), function(cluster) {
    in_cluster <- clustering$metacluster == cluster
    data.frame(
      tissue = tissue,
      cluster = cluster,
      events = sum(in_cluster),
      asc_events = sum(in_cluster & truth),
      purity_percent = 100 * sum(in_cluster & truth) / sum(in_cluster),
      recall_percent = 100 * sum(in_cluster & truth) / sum(truth),
      stringsAsFactors = FALSE
    )
  })
  per_cluster <- do.call(rbind, per_cluster)
  per_cluster <- per_cluster[order(-per_cluster$purity_percent), ]

  best <- per_cluster[1, ]
  # An ASC cluster is one where the manual gate calls most of its events ASC.
  asc_clusters <- per_cluster[per_cluster$purity_percent >= 50, ]

  cluster_rows[[length(cluster_rows) + 1]] <- data.frame(
    tissue = tissue,
    manual_asc = sum(truth),
    manual_percent = 100 * mean(truth),
    best_cluster_purity = best$purity_percent,
    best_cluster_recall = best$recall_percent,
    asc_clusters = nrow(asc_clusters),
    combined_recall = sum(asc_clusters$recall_percent),
    combined_percent = 100 * sum(asc_clusters$events) / length(truth),
    stringsAsFactors = FALSE
  )
  Log("  ", tissue, sprintf("best cluster %.1f%% pure, %.1f%% recall",
                            best$purity_percent, best$recall_percent))

  write.csv(per_cluster,
            file.path(kOutputDir, paste0("clusters_",
                                         gsub(" ", "_", tolower(tissue)), ".csv")),
            row.names = FALSE)

  # One UMAP, from the tissue with the most ASC, so the figure shows them clearly.
  if (tissue == "Spleen" && !umap_saved) {
    embedding <- RunUmapEmbedding(events, channels = kClusterChannels, seed = kSeed)
    plot_data <- cbind(
      embedding,
      manual_label = ifelse(truth, "ASC by the manual gate", "other"),
      cluster = factor(clustering$metacluster)
    )
    write.csv(plot_data[sample.int(nrow(plot_data), min(20000, nrow(plot_data))), ],
              file.path(kOutputDir, "umap_spleen.csv"), row.names = FALSE)

    umap_manual <- ggplot(plot_data[order(plot_data$manual_label == "other",
                                          decreasing = TRUE), ],
                          aes(x = umap_1, y = umap_2, colour = manual_label)) +
      geom_point(size = 0.3, alpha = 0.6) +
      scale_colour_manual(values = c("ASC by the manual gate" = "#b2182b",
                                     other = "grey75")) +
      guides(colour = guide_legend(override.aes = list(size = 3, alpha = 1))) +
      labs(
        title = "Spleen, UMAP of the dumped population",
        subtitle = "Red is the antibody secreting cell gate the authors drew. The embedding never saw that gate.",
        x = "UMAP 1", y = "UMAP 2", colour = NULL
      ) +
      theme_bw()
    ggsave(file.path(kOutputDir, "umap_spleen_manual.png"), umap_manual,
           width = 9, height = 7, dpi = 150)

    umap_cluster <- ggplot(plot_data, aes(x = umap_1, y = umap_2, colour = cluster)) +
      geom_point(size = 0.3, alpha = 0.6) +
      guides(colour = guide_legend(override.aes = list(size = 3, alpha = 1))) +
      labs(
        title = "The same embedding, coloured by metacluster",
        x = "UMAP 1", y = "UMAP 2", colour = "Metacluster"
      ) +
      theme_bw()
    ggsave(file.path(kOutputDir, "umap_spleen_clusters.png"), umap_cluster,
           width = 9, height = 7, dpi = 150)
    umap_saved <- TRUE
    Log("  wrote the spleen UMAP figures")
  }
}

clustering_summary <- do.call(rbind, cluster_rows)
write.csv(clustering_summary, file.path(kOutputDir, "clustering_summary.csv"),
          row.names = FALSE)

cat("\n=== Clustering against the manual gate ===\n")
print(clustering_summary, digits = 3, row.names = FALSE)

# ---------------------------------------------------------------------------
# 4b. Why the density gate cannot work, shown rather than asserted
# ---------------------------------------------------------------------------

Log("Drawing the CD38 distribution with both gate decisions on it")

evidence_rows <- list()
manual_gate <- NULL
for (tissue in kTissues) {
  gating_hierarchy <- gating_sets[[tissue]][[1]]
  parent_data <- flowWorkspace::gh_pop_get_data(gating_hierarchy, kParentPath)
  events <- flowCore::exprs(parent_data)

  in_parent <- flowWorkspace::gh_pop_get_indices(gating_hierarchy, kParentPath)
  in_asc <- flowWorkspace::gh_pop_get_indices(gating_hierarchy, kAscPath)
  truth <- in_asc[in_parent]

  gate <- flowWorkspace::gh_pop_get_gate(gating_hierarchy, kAscPath)
  if (tissue == "Spleen") manual_gate <- gate
  manual_low <- gate@min[[kChannels[["CD38"]]]]

  # Where mindensity would cut, on exactly the data the analyst looked at.
  found <- tryCatch(
    openCyto::gate_mindensity(parent_data, channel = kChannels[["CD38"]]),
    error = function(e) NULL
  )
  automatic_cut <- if (is.null(found)) NA_real_ else found@min[[1]]

  evidence_rows[[length(evidence_rows) + 1]] <- data.frame(
    tissue = tissue,
    manual_cd38_threshold = manual_low,
    mindensity_cd38_threshold = automatic_cut,
    percent_above_manual = 100 * mean(events[, kChannels[["CD38"]]] >= manual_low),
    percent_above_mindensity = if (is.na(automatic_cut)) NA_real_ else
      100 * mean(events[, kChannels[["CD38"]]] >= automatic_cut),
    manual_asc_percent = 100 * mean(truth),
    stringsAsFactors = FALSE
  )
}
cd38_evidence <- do.call(rbind, evidence_rows)
write.csv(cd38_evidence, file.path(kOutputDir, "cd38_thresholds.csv"),
          row.names = FALSE)

cat("\n=== Where each method puts the CD38 cut, on the same data ===\n")
print(cd38_evidence, digits = 4, row.names = FALSE)

spleen_hierarchy <- gating_sets[["Spleen"]][[1]]
spleen_events <- flowCore::exprs(
  flowWorkspace::gh_pop_get_data(spleen_hierarchy, kParentPath)
)
spleen_truth <- flowWorkspace::gh_pop_get_indices(spleen_hierarchy, kAscPath)[
  flowWorkspace::gh_pop_get_indices(spleen_hierarchy, kParentPath)
]
spleen_row <- cd38_evidence[cd38_evidence$tissue == "Spleen", ]

density_plot <- ggplot(data.frame(cd38 = spleen_events[, kChannels[["CD38"]]]),
                       aes(x = cd38)) +
  geom_density(fill = "grey85", colour = "grey40", linewidth = 0.4) +
  geom_vline(xintercept = spleen_row$manual_cd38_threshold,
             colour = "#b2182b", linewidth = 0.9) +
  geom_vline(xintercept = spleen_row$mindensity_cd38_threshold,
             colour = "#2166ac", linewidth = 0.9, linetype = "dashed") +
  annotate("text", x = spleen_row$manual_cd38_threshold, y = Inf,
           label = "  the analyst cut here", colour = "#b2182b",
           hjust = 0, vjust = 2, size = 3.6) +
  annotate("text", x = spleen_row$mindensity_cd38_threshold, y = Inf,
           label = "mindensity cuts here  ", colour = "#2166ac",
           hjust = 1, vjust = 4, size = 3.6) +
  labs(
    title = "Spleen: the CD38 distribution has no valley where the analyst cut",
    subtitle = "mindensity finds the dip between negative and positive. The CD38 high population sits far above it, on a smooth shoulder.",
    x = "CD38 intensity, on the workspace scale", y = "Density"
  ) +
  theme_bw()
ggsave(file.path(kOutputDir, "cd38_density.png"), density_plot,
       width = 10, height = 6, dpi = 150)

withr::with_seed(kSeed, {
  show <- sample.int(nrow(spleen_events), min(40000, nrow(spleen_events)))
})
biaxial_plot <- ggplot(
  data.frame(
    cd38 = spleen_events[show, kChannels[["CD38"]]],
    cd27 = spleen_events[show, kChannels[["CD27"]]],
    asc = ifelse(spleen_truth[show], "ASC", "other"),
    stringsAsFactors = FALSE
  ),
  aes(x = cd38, y = cd27)
) +
  geom_point(aes(colour = asc), size = 0.25, alpha = 0.35) +
  annotate("rect",
           xmin = manual_gate@min[[kChannels[["CD38"]]]],
           xmax = manual_gate@max[[kChannels[["CD38"]]]],
           ymin = manual_gate@min[[kChannels[["CD27"]]]],
           ymax = manual_gate@max[[kChannels[["CD27"]]]],
           fill = NA, colour = "#b2182b", linewidth = 0.8) +
  geom_vline(xintercept = spleen_row$mindensity_cd38_threshold,
             colour = "#2166ac", linewidth = 0.8, linetype = "dashed") +
  scale_colour_manual(values = c(ASC = "#b2182b", other = "grey70")) +
  guides(colour = guide_legend(override.aes = list(size = 3, alpha = 1))) +
  labs(
    title = "Spleen: CD38 against CD27, with both decisions drawn",
    subtitle = "The red box is the rectangleGate the analyst placed. The dashed line is where mindensity cuts CD38.",
    x = "CD38", y = "CD27", colour = NULL
  ) +
  theme_bw()
ggsave(file.path(kOutputDir, "cd38_cd27_biaxial.png"), biaxial_plot,
       width = 10, height = 7, dpi = 150)

Log("Wrote cd38_density.png and cd38_cd27_biaxial.png")

# ---------------------------------------------------------------------------
# 4c. What the metaclusters are
# ---------------------------------------------------------------------------

Log("Annotating the spleen metaclusters")

definitions <- ReadCellTypeDefinitions(kDefinitionsPath)

withr::with_seed(kSeed, {
  keep <- sample.int(nrow(spleen_events), min(kSubsampleSize, nrow(spleen_events)))
})
spleen_sub <- spleen_events[keep, , drop = FALSE]
spleen_truth_sub <- spleen_truth[keep]

spleen_clustering <- RunFlowSomClustering(
  spleen_sub, channels = kClusterChannels,
  grid_size = 10, n_metaclusters = kMetaclusters, seed = kSeed
)

median_expression <- ClusterMedianExpression(
  spleen_sub, spleen_clustering$metacluster, kClusterChannels
)
write.csv(median_expression,
          file.path(kOutputDir, "spleen_cluster_median_expression.csv"),
          row.names = FALSE)

# The score reads a scaled version of the median table, not the raw medians, so
# the scaled table is written out too. Without it a reader cannot check a label.
ScaleColumn <- function(x) {
  span <- max(x) - min(x)
  if (span == 0) rep(0.5, length(x)) else (x - min(x)) / span
}
scaled_expression <- median_expression[, c("cluster", "events", "percent_of_total")]
for (i in seq_along(kClusterChannels)) {
  scaled_expression[[names(kChannels)[i]]] <-
    round(ScaleColumn(median_expression[[kClusterChannels[i]]]), 3)
}
write.csv(scaled_expression,
          file.path(kOutputDir, "spleen_cluster_scaled_expression.csv"),
          row.names = FALSE)

annotation <- AnnotateClusters(median_expression, definitions)
annotation$asc_percent_by_manual_gate <- vapply(
  annotation$cluster,
  function(cluster) {
    in_cluster <- spleen_clustering$metacluster == cluster
    100 * sum(in_cluster & spleen_truth_sub) / sum(in_cluster)
  },
  numeric(1)
)
write.csv(annotation, file.path(kOutputDir, "spleen_cluster_annotation.csv"),
          row.names = FALSE)

cat("\n=== What the metaclusters are ===\n")
print(annotation[, c("cluster", "events", "percent_of_total", "cell_type",
                     "runner_up", "margin", "asc_percent_by_manual_gate")],
      digits = 3, row.names = FALSE)

marker_names <- names(kChannels)
expression_long <- do.call(rbind, lapply(seq_along(kClusterChannels), function(i) {
  data.frame(
    cluster = factor(median_expression$cluster),
    marker = marker_names[i],
    median = median_expression[[kClusterChannels[i]]],
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

heatmap_plot <- ggplot(expression_long,
                       aes(x = marker, y = cluster, fill = scaled)) +
  geom_tile(colour = "white") +
  scale_fill_viridis_c(option = "magma") +
  labs(
    title = "Spleen: median marker expression per metacluster",
    subtitle = "Scaled within each marker. This table drives the labels in the annotation table.",
    x = NULL, y = "Metacluster", fill = "Scaled\nmedian"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(file.path(kOutputDir, "spleen_cluster_heatmap.png"), heatmap_plot,
       width = 9, height = 6, dpi = 150)

spleen_embedding <- RunUmapEmbedding(spleen_sub, channels = kClusterChannels,
                                     seed = kSeed)
marker_long <- do.call(rbind, lapply(seq_along(kClusterChannels), function(i) {
  data.frame(
    umap_1 = spleen_embedding$umap_1,
    umap_2 = spleen_embedding$umap_2,
    marker = marker_names[i],
    value = spleen_sub[, kClusterChannels[i]],
    stringsAsFactors = FALSE
  )
}))

marker_plot <- ggplot(marker_long, aes(x = umap_1, y = umap_2, colour = value)) +
  geom_point(size = 0.15, alpha = 0.5) +
  facet_wrap(~marker, ncol = 4) +
  scale_colour_viridis_c(option = "magma") +
  labs(
    title = "Spleen: the same embedding, one panel per marker",
    subtitle = "Read the annotation table against these panels. A label that does not match them is wrong.",
    x = "UMAP 1", y = "UMAP 2", colour = "Intensity"
  ) +
  theme_bw()
ggsave(file.path(kOutputDir, "spleen_umap_markers.png"), marker_plot,
       width = 13, height = 7, dpi = 150)

type_plot <- ggplot(
  cbind(spleen_embedding,
        cell_type = annotation$cell_type[
          match(spleen_clustering$metacluster, annotation$cluster)
        ]),
  aes(x = umap_1, y = umap_2, colour = cell_type)
) +
  geom_point(size = 0.3, alpha = 0.6) +
  guides(colour = guide_legend(override.aes = list(size = 3, alpha = 1))) +
  labs(
    title = "Spleen: metaclusters labelled from the definitions file",
    subtitle = "No cluster was named by eye. Each was scored against gating/omip43_cell_type_definitions.csv.",
    x = "UMAP 1", y = "UMAP 2", colour = "Cell type"
  ) +
  theme_bw()
ggsave(file.path(kOutputDir, "spleen_umap_cell_types.png"), type_plot,
       width = 10, height = 7, dpi = 150)

Log("Wrote the annotation table and three more figures")

# ---------------------------------------------------------------------------
# 5. The paper's phenotype claims
# ---------------------------------------------------------------------------

Log("Testing the phenotype claims on the manual gates")

phenotype_rows <- list()
for (tissue in kTissues) {
  gating_hierarchy <- gating_sets[[tissue]][[1]]
  in_parent <- flowWorkspace::gh_pop_get_indices(gating_hierarchy, kParentPath)
  in_asc <- flowWorkspace::gh_pop_get_indices(gating_hierarchy, kAscPath)
  truth <- in_asc[in_parent]

  events <- flowCore::exprs(
    flowWorkspace::gh_pop_get_data(gating_hierarchy, kParentPath)
  )

  MedianIn <- function(channel, subset) stats::median(events[subset, channel])

  phenotype_rows[[length(phenotype_rows) + 1]] <- data.frame(
    tissue = tissue,
    asc_events = sum(truth),
    ssc_asc = MedianIn("SSC-A", truth),
    ssc_other = MedianIn("SSC-A", !truth),
    cd20_asc = MedianIn(kChannels[["CD20"]], truth),
    cd20_other = MedianIn(kChannels[["CD20"]], !truth),
    cd19_asc = MedianIn(kChannels[["CD19"]], truth),
    hladr_asc = MedianIn(kChannels[["HLA-DR"]], truth),
    cd38_asc = MedianIn(kChannels[["CD38"]], truth),
    cd38_other = MedianIn(kChannels[["CD38"]], !truth),
    stringsAsFactors = FALSE
  )
}
phenotype <- do.call(rbind, phenotype_rows)
phenotype$ssc_ratio <- phenotype$ssc_asc / phenotype$ssc_other
phenotype$cd20_difference <- phenotype$cd20_asc - phenotype$cd20_other
write.csv(phenotype, file.path(kOutputDir, "phenotype.csv"), row.names = FALSE)

cat("\n=== Phenotype of ASC against the rest of the dumped population ===\n")
print(phenotype[, c("tissue", "ssc_asc", "ssc_other", "ssc_ratio",
                    "cd20_asc", "cd20_other", "cd20_difference")],
      digits = 3, row.names = FALSE)

cat("\n=== Tissue variation among ASC, median intensity ===\n")
print(phenotype[, c("tissue", "cd19_asc", "hladr_asc", "cd38_asc")],
      digits = 3, row.names = FALSE)

tissue_variation <- data.frame(
  marker = c("CD19", "HLA-DR"),
  cv_across_tissues = c(MeasuredCv(phenotype$cd19_asc),
                        MeasuredCv(phenotype$hladr_asc)),
  stringsAsFactors = FALSE
)
write.csv(tissue_variation, file.path(kOutputDir, "tissue_variation.csv"),
          row.names = FALSE)

# ---------------------------------------------------------------------------
# 5b. Every claim in the file, with a verdict
# ---------------------------------------------------------------------------
#
# The tables above show the measurements. This section turns each one into a
# verdict against the claim it was meant to test, so no claim is quietly left
# without an answer.

Log("Scoring every recorded claim")

Verdict <- function(claim_id, observed, verdict, note) {
  data.frame(claim_id = claim_id, observed = observed, verdict = verdict,
             note = note, stringsAsFactors = FALSE)
}

pbmc_percent <- mean(manual$asc_percent[manual$tissue == "PBMC"])
other_percent <- vapply(setdiff(kTissues, "PBMC"), function(tissue) {
  mean(manual$asc_percent[manual$tissue == tissue])
}, numeric(1))

verdict_rows <- list(
  Verdict(
    1,
    sprintf("PBMC %.2f%%, others %.2f to %.2f%%", pbmc_percent,
            min(other_percent), max(other_percent)),
    if (all(pbmc_percent < other_percent)) "reproduced" else "not reproduced",
    "Blood is the lowest of the four tissues."
  ),
  Verdict(
    2,
    sprintf("%d of %d replicates inside 400 to 2000",
            sum(counting_claims$samples_in_range), nrow(manual)),
    if (sum(counting_claims$samples_in_range) == nrow(manual)) "reproduced" else
      "not reproduced",
    "Blood and bone marrow fall below the range, spleen and tonsil above it."
  ),
  Verdict(
    3,
    paste0(sum(counting_claims$meets_cv_target), " of ",
           nrow(counting_claims), " tissues at or below 5 percent"),
    if (all(counting_claims$meets_cv_target)) "reproduced" else "partly reproduced",
    "Met in spleen and tonsil, missed in blood and bone marrow."
  ),
  Verdict(
    4,
    sprintf("side scatter ratio %.2f to %.2f",
            min(phenotype$ssc_ratio), max(phenotype$ssc_ratio)),
    if (all(phenotype$ssc_ratio > 1)) "reproduced" else "not reproduced",
    "ASC scatter more light than the cells around them in every tissue."
  ),
  Verdict(
    5,
    sprintf("CD20 difference %.1f to %.1f",
            min(phenotype$cd20_difference), max(phenotype$cd20_difference)),
    if (all(phenotype$cd20_difference < 0)) "reproduced" else "not reproduced",
    "CD20 is lower on ASC than on the rest of the dumped population everywhere."
  ),
  Verdict(
    6,
    sprintf("CD19 median %.0f to %.0f across tissues, CV %.1f%%",
            min(phenotype$cd19_asc), max(phenotype$cd19_asc),
            MeasuredCv(phenotype$cd19_asc)),
    if (MeasuredCv(phenotype$cd19_asc) > 5) "reproduced" else "too small to call",
    "Lowest in bone marrow and highest in tonsil."
  ),
  Verdict(
    7,
    sprintf("HLA-DR median %.0f to %.0f across tissues, CV %.1f%%",
            min(phenotype$hladr_asc), max(phenotype$hladr_asc),
            MeasuredCv(phenotype$hladr_asc)),
    if (MeasuredCv(phenotype$hladr_asc) > 5) "reproduced" else "too small to call",
    "The paper also names Ki67 here. The PE channel in these files carries HLA-DR, so Ki67 is not measurable."
  )
)

verdicts <- merge(claims, do.call(rbind, verdict_rows), by = "claim_id")
verdicts <- verdicts[order(verdicts$claim_id), ]
write.csv(verdicts, file.path(kOutputDir, "claims_verdicts.csv"), row.names = FALSE)

cat("\n=== Every claim, with a verdict ===\n")
print(verdicts[, c("claim_id", "short_name", "expected", "observed", "verdict")],
      row.names = FALSE)

Log("Claims reproduced:", sum(verdicts$verdict == "reproduced"), "of",
    nrow(verdicts))

# ---------------------------------------------------------------------------
# 5c. The same claims again, from a population the manual gate never touched
# ---------------------------------------------------------------------------
#
# Everything above defines ASC by the gate the authors drew. That tests whether
# the paper's conclusions follow from its own analysis, which is worth knowing and
# is not the same as reproducing them independently.
#
# Here the ASC population is defined by clustering instead. Each tissue is
# clustered, every metacluster is labelled from the definitions file, and the
# clusters labelled as antibody secreting cells become the population. The manual
# gate is not consulted at any point, so the claims below are tested against cells
# selected by a route that shares nothing with the original analysis.

Log("Rebuilding the ASC population from clustering, per tissue")

cluster_phenotype_rows <- list()
for (tissue in kTissues) {
  gating_hierarchy <- gating_sets[[tissue]][[1]]
  events_all <- flowCore::exprs(
    flowWorkspace::gh_pop_get_data(gating_hierarchy, kParentPath)
  )
  truth_all <- flowWorkspace::gh_pop_get_indices(gating_hierarchy, kAscPath)[
    flowWorkspace::gh_pop_get_indices(gating_hierarchy, kParentPath)
  ]

  withr::with_seed(kSeed, {
    keep <- sample.int(nrow(events_all), min(kSubsampleSize, nrow(events_all)))
  })
  events <- events_all[keep, , drop = FALSE]
  truth <- truth_all[keep]

  clusters <- RunFlowSomClustering(
    events, channels = kClusterChannels,
    grid_size = 14, n_metaclusters = kMetaclusters, seed = kSeed
  )
  medians <- ClusterMedianExpression(events, clusters$metacluster, kClusterChannels)
  labels <- AnnotateClusters(medians, definitions)

  # The paper identifies ASC by CD38 being highest, not by the shape of the whole
  # profile: "very high CD38 expression is in fact considered adequate for basic
  # identification of ASC". Ranking clusters on CD38 encodes that sentence.
  # Scoring the profile instead selects a second CD38 positive population that
  # bone marrow carries and the other tissues do not, and the F1 there falls from
  # 90.9 to 12.6. Script 07 holds that measurement.
  asc_clusters <- SelectByHighestMarker(
    medians, kChannels[["CD38"]], n_clusters = kAscClusters
  )
  is_asc <- clusters$metacluster %in% asc_clusters

  if (sum(is_asc) < 50) {
    Log("  ", tissue, "no usable ASC cluster, skipping")
    next
  }

  MedianIn <- function(channel, subset) stats::median(events[subset, channel])

  cluster_phenotype_rows[[length(cluster_phenotype_rows) + 1]] <- data.frame(
    tissue = tissue,
    asc_clusters = length(asc_clusters),
    asc_events = sum(is_asc),
    asc_percent = 100 * mean(is_asc),
    manual_percent = 100 * mean(truth),
    agreement_percent = 100 * sum(is_asc & truth) / sum(is_asc),
    ssc_asc = MedianIn("SSC-A", is_asc),
    ssc_other = MedianIn("SSC-A", !is_asc),
    cd20_asc = MedianIn(kChannels[["CD20"]], is_asc),
    cd20_other = MedianIn(kChannels[["CD20"]], !is_asc),
    cd19_asc = MedianIn(kChannels[["CD19"]], is_asc),
    hladr_asc = MedianIn(kChannels[["HLA-DR"]], is_asc),
    stringsAsFactors = FALSE
  )
  Log("  ", tissue, sprintf("%d cluster(s), %.2f%% of the parent against %.2f%% manual",
                            length(asc_clusters), 100 * mean(is_asc), 100 * mean(truth)))
}

cluster_phenotype <- do.call(rbind, cluster_phenotype_rows)
cluster_phenotype$ssc_ratio <- cluster_phenotype$ssc_asc / cluster_phenotype$ssc_other
cluster_phenotype$cd20_difference <- cluster_phenotype$cd20_asc -
  cluster_phenotype$cd20_other
write.csv(cluster_phenotype, file.path(kOutputDir, "cluster_route_phenotype.csv"),
          row.names = FALSE)

cat("\n=== ASC defined by clustering alone ===\n")
print(cluster_phenotype[, c("tissue", "asc_percent", "manual_percent",
                            "agreement_percent", "ssc_ratio", "cd20_difference")],
      digits = 3, row.names = FALSE)

# The same claims, scored against this population.
cluster_pbmc <- cluster_phenotype$asc_percent[cluster_phenotype$tissue == "PBMC"]
cluster_others <- cluster_phenotype$asc_percent[cluster_phenotype$tissue != "PBMC"]

route_rows <- list(
  Verdict(1,
          sprintf("PBMC %.2f%%, others %.2f to %.2f%%", cluster_pbmc,
                  min(cluster_others), max(cluster_others)),
          if (all(cluster_pbmc < cluster_others)) "reproduced" else "not reproduced",
          "Tested on clusters, not on the manual gate."),
  Verdict(4,
          sprintf("side scatter ratio %.2f to %.2f",
                  min(cluster_phenotype$ssc_ratio), max(cluster_phenotype$ssc_ratio)),
          if (all(cluster_phenotype$ssc_ratio > 1)) "reproduced" else "not reproduced",
          "Tested on clusters, not on the manual gate."),
  Verdict(5,
          sprintf("CD20 difference %.1f to %.1f",
                  min(cluster_phenotype$cd20_difference),
                  max(cluster_phenotype$cd20_difference)),
          if (all(cluster_phenotype$cd20_difference < 0)) "reproduced" else
            "not reproduced",
          "Tested on clusters, not on the manual gate."),
  Verdict(6,
          sprintf("CD19 median %.0f to %.0f, CV %.1f%%",
                  min(cluster_phenotype$cd19_asc), max(cluster_phenotype$cd19_asc),
                  MeasuredCv(cluster_phenotype$cd19_asc)),
          if (MeasuredCv(cluster_phenotype$cd19_asc) > 5) "reproduced" else
            "too small to call",
          "Tested on clusters, not on the manual gate."),
  Verdict(7,
          sprintf("HLA-DR median %.0f to %.0f, CV %.1f%%",
                  min(cluster_phenotype$hladr_asc), max(cluster_phenotype$hladr_asc),
                  MeasuredCv(cluster_phenotype$hladr_asc)),
          if (MeasuredCv(cluster_phenotype$hladr_asc) > 5) "reproduced" else
            "too small to call",
          "Tested on clusters, not on the manual gate.")
)

route_two <- merge(claims, do.call(rbind, route_rows), by = "claim_id")
route_two <- route_two[order(route_two$claim_id), ]

both_routes <- rbind(
  cbind(verdicts[, c("claim_id", "short_name", "expected", "observed", "verdict")],
        route = "the authors' gates"),
  cbind(route_two[, c("claim_id", "short_name", "expected", "observed", "verdict")],
        route = "clustering")
)
both_routes <- both_routes[order(both_routes$claim_id), ]
write.csv(both_routes, file.path(kOutputDir, "claims_both_routes.csv"),
          row.names = FALSE)

cat("\n=== The same claims by two independent routes ===\n")
print(both_routes[, c("claim_id", "short_name", "route", "observed", "verdict")],
      row.names = FALSE)

# ---------------------------------------------------------------------------
# 6. Figures
# ---------------------------------------------------------------------------

Log("Writing figures")

count_plot <- ggplot(manual, aes(x = tissue, y = asc_events)) +
  annotate("rect", xmin = -Inf, xmax = Inf,
           ymin = kTargetEventsLow, ymax = kTargetEventsHigh,
           alpha = 0.15, fill = "#2166ac") +
  geom_point(size = 3, alpha = 0.8) +
  scale_y_log10() +
  labs(
    title = "Events collected in the ASC gate, by tissue",
    subtitle = paste0(
      "The band is the paper's stated range of ", kTargetEventsLow, " to ",
      kTargetEventsHigh, " events. Seven replicates per tissue."
    ),
    x = NULL, y = "Events in the ASC gate, log scale"
  ) +
  theme_bw()
ggsave(file.path(kOutputDir, "asc_event_counts.png"), count_plot,
       width = 9, height = 6, dpi = 150)

cv_long <- rbind(
  data.frame(tissue = spread$tissue, source = "Poisson floor",
             cv = spread$poisson_cv_percent, stringsAsFactors = FALSE),
  data.frame(tissue = spread$tissue, source = "Measured across replicates",
             cv = spread$measured_cv_percent, stringsAsFactors = FALSE)
)
cv_plot <- ggplot(cv_long, aes(x = tissue, y = cv, fill = source)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.7) +
  geom_hline(yintercept = kTargetCvPercent, linetype = "dashed") +
  annotate("text", x = 0.7, y = kTargetCvPercent + 0.4,
           label = "the paper's 5 percent target", hjust = 0, size = 3) +
  labs(
    title = "Counting noise against measured spread",
    subtitle = "The Poisson floor is set by the event count alone. The measured CV is what seven replicates did.",
    x = NULL, y = "Coefficient of variation, percent", fill = NULL
  ) +
  theme_bw()
ggsave(file.path(kOutputDir, "cv_comparison.png"), cv_plot,
       width = 9, height = 6, dpi = 150)

Log("Done. Output is in", kOutputDir)
