#!/usr/bin/env Rscript

# OMIP-043: antibody secreting cells across four tissues.
#
# This dataset asks a different question from OMIP-39. There the population was
# separated on one marker and an automated template found it. Here the population
# is the bright tail of a continuum, and it is not separable that way.
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
kMetaclusters <- 15
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
