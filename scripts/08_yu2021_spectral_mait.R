#!/usr/bin/env Rscript

# Reproduce the flow cytometry findings of Yu 2021 from the data the authors
# deposited, on a 35 marker spectral panel.
#
# The dataset is FlowRepository FR-FCM-Z3WR, which is 83 Cytek Aurora files from
# Yu C, et al. Mucosal-associated invariant T cell responses differ by sex in
# COVID-19. Med 2021;2:755-772.e5. PMID 33870241.
#
# Four things make this different from the OMIP-039 report.
#
# The instrument is spectral, so the file that left it was already unmixed by
# SpectroFlo. There is no spillover matrix to compute and no compensation step.
# APPLY COMPENSATION reads FALSE in every file.
#
# The claim is about a cohort rather than a file. It compares 83 samples from 45
# subjects across four severity ranks and two sexes, so the clinical grouping has
# to come from somewhere. It comes from Table S1 of the paper, joined to the files
# in gating/yu2021_sample_metadata.csv, and all 83 files match a metadata row.
#
# Every sample is gated at the same cut, listed in gating/yu2021_gate_cuts.csv. A
# cut fitted per sample moves the boundary for reasons that are not biology, and
# every claim here is a comparison between samples. The cuts were fitted on a
# pooled subsample of all 83 files and six of the seven sit on a real density
# minimum. The seventh, CD45RA, does not, and the sweep in step 6 shows what that
# costs.
#
# The paper reads MAIT cells as CD8+ CD161hi T cells, because the flow panel
# carries no Va7.2. It validates that surrogate with scRNA-seq, where the CD161hi
# cluster co-expresses KLRB1, CD3D, CD8A and TRAV1-2. This script tests the flow
# claims only, so it inherits the surrogate and does not check it.
#
# Run it in the container:
#   podman run --rm -v "$PWD:/work:z" -w /work everything-flow-cytometry:latest \
#     Rscript scripts/08_yu2021_spectral_mait.R

suppressPackageStartupMessages({
  library(flowCore)
  library(FlowSOM)
  library(uwot)
  library(ggplot2)
  library(withr)
})

for (file in list.files("R", pattern = "\\.R$", full.names = TRUE)) {
  source(file)
}

kDataDir <- file.path(
  "data", "datasets", "flowrepository", "FlowRepository_FR-FCM-Z3WR_files"
)
kSheetPath <- file.path("gating", "yu2021_sample_metadata.csv")
kCutsPath <- file.path("gating", "yu2021_gate_cuts.csv")
kClaimsPath <- file.path("gating", "yu2021_paper_claims.csv")
kOutputDir <- file.path("output", "yu2021")

# The paper pooled every participant and drew 3000 live CD45+ cells per sample for
# its clustering figure. The same draw is used here so the two are comparable.
kEventsPerSample <- 3000

# Cluster on lineage markers only, so a cluster is an identity rather than a
# state. These are the markers the paper names in Figures 1C and 1E.
kLineageMarkers <- c(
  "CD3", "CD4", "CD8", "CD56", "CD19", "CD20", "CD14", "CD16",
  "CD11c", "CD123", "CD161", "TCR gd", "CD45ra", "CCR7", "HLA-DR"
)

# The CD45RA cut has no density minimum to sit on, so the memory claims are run
# again at every value in this sweep.
kCd45raSweep <- seq(2.4, 4.2, by = 0.2)

kMetaclusters <- 15
kSeed <- 42

Log <- function(...) {
  cat(format(Sys.time(), "%H:%M:%S"), "|", ..., "\n")
}

dir.create(kOutputDir, recursive = TRUE, showWarnings = FALSE)

# ---------------------------------------------------------------------------
# Step 1: the cohort
# ---------------------------------------------------------------------------

Log("Reading the sample sheet")
sheet <- ReadYuSampleSheet(kSheetPath)

files <- file.path(kDataDir, sheet$file_name)
absent <- sheet$file_name[!file.exists(files)]
if (length(absent) > 0) {
  stop(
    "These files are named in the sample sheet and are not on disk: ",
    paste(utils::head(absent, 5), collapse = ", "),
    ". Pull the dataset with ./sync.sh pull ",
    "datasets/flowrepository/FlowRepository_FR-FCM-Z3WR_files"
  )
}
Log("The cohort is", nrow(sheet), "files from",
    length(unique(sheet$alias_subject_id)), "subjects")

design <- as.data.frame(table(sheet$severity_rank, sheet$sex))
colnames(design) <- c("severity_rank", "sex", "samples")
subjects <- as.data.frame(table(
  unique(sheet[, c("alias_subject_id", "severity_rank", "sex")])$severity_rank,
  unique(sheet[, c("alias_subject_id", "severity_rank", "sex")])$sex
))
colnames(subjects) <- c("severity_rank", "sex", "subjects")
design <- merge(design, subjects, by = c("severity_rank", "sex"))
write.csv(design, file.path(kOutputDir, "cohort_design.csv"), row.names = FALSE)
cat("\n=== The cohort ===\n")
print(design, row.names = FALSE)

# ---------------------------------------------------------------------------
# Step 2: the panel, the transform and the cut points
# ---------------------------------------------------------------------------

Log("Reading the panel from", basename(files[1]))
reference_frame <- read.FCS(files[1], truncate_max_range = FALSE)
panel <- DescribeChannels(reference_frame)
write.csv(panel, file.path(kOutputDir, "panel.csv"), row.names = FALSE)
Log("The panel holds", nrow(panel), "parameters and",
    sum(panel$is_marker), "named markers")

# The viability channel carries no marker name, so FluorescenceChannels() would
# drop it. SpectralChannels() keeps it, and the first gate needs it.
transform_channels <- SpectralChannels(reference_frame)
Log("Estimating one logicle transform for", length(transform_channels),
    "channels")
transform_result <- EstimateSpectralTransform(files[1], transform_channels)

cuts <- ReadGateCuts(kCutsPath)
write.csv(cuts, file.path(kOutputDir, "gate_cuts.csv"), row.names = FALSE)
cat("\n=== The gate cut points, shared by every sample ===\n")
print(cuts[, c("marker", "parent", "side", "cut", "source")], row.names = FALSE)

# ---------------------------------------------------------------------------
# Step 3: gate every file
# ---------------------------------------------------------------------------

all_counts <- list()
all_events <- list()
all_medians <- list()
failures <- list()

for (i in seq_len(nrow(sheet))) {
  file_name <- sheet$file_name[i]
  if (i %% 10 == 1) {
    Log(sprintf("Gating %d of %d", i, nrow(sheet)))
  }

  # The chosen CD45RA cut comes first, so counts[1, ] is the primary analysis and
  # the rest of the rows answer the sweep from the same read.
  gated <- GateSpectralFile(
    path = files[i],
    transform_list = transform_result$transform,
    cuts = cuts,
    panel = panel,
    cd45ra_cut = c(cuts$cut[cuts$marker == "CD45ra"], kCd45raSweep),
    subsample_n = kEventsPerSample,
    seed = kSeed + i
  )

  if (!is.na(gated$error)) {
    Log("  FAILED", file_name, substr(gated$error, 1, 120))
    failures[[file_name]] <- data.frame(
      file_name = file_name, message = gated$error, stringsAsFactors = FALSE
    )
    next
  }

  all_counts[[file_name]] <- gated$counts
  all_medians[[file_name]] <- gated$medians
  if (!is.null(gated$events)) {
    all_events[[file_name]] <- gated$events
  }
}

if (length(all_counts) == 0) {
  stop("Every file failed, so there is nothing to report.")
}

all_rows <- do.call(rbind, all_counts)
rownames(all_rows) <- NULL

primary_cut <- cuts$cut[cuts$marker == "CD45ra"]
counts <- all_rows[all_rows$cd45ra_cut == primary_cut, , drop = FALSE]
sweep_counts <- all_rows[all_rows$cd45ra_cut %in% kCd45raSweep, , drop = FALSE]

failure_table <- if (length(failures) > 0) {
  do.call(rbind, failures)
} else {
  data.frame(file_name = character(), message = character(),
             stringsAsFactors = FALSE)
}
write.csv(failure_table, file.path(kOutputDir, "gating_failures.csv"),
          row.names = FALSE)
Log(nrow(counts), "files gated,", nrow(failure_table), "failed")

# ---------------------------------------------------------------------------
# Step 4: per sample frequencies, joined to the clinical grouping
# ---------------------------------------------------------------------------

Log("Computing the frequencies")
counts <- AddCd8Frequencies(counts)
medians <- do.call(rbind, all_medians)
counts <- merge(counts, medians, by = "sample", all.x = TRUE)

measures <- merge(counts, sheet, by.x = "sample", by.y = "file_name",
                  all.x = TRUE)
measures <- measures[order(measures$sample), ]
write.csv(measures, file.path(kOutputDir, "sample_measures.csv"),
          row.names = FALSE)

cat("\n=== The gating hierarchy, median across the 83 samples ===\n")
hierarchy <- data.frame(
  step = c("Events in the file", "Singlets", "Live", "CD45 positive",
           "CD3 positive", "CD8 positive of CD3", "CD161hi of CD8",
           "Memory of CD8", "Naive of CD8", "NK of CD45 positive"),
  median_value = c(
    stats::median(measures$total_events),
    stats::median(measures$singlet_events),
    stats::median(measures$live_events),
    stats::median(measures$cd45pos_events),
    stats::median(measures$cd3pos_events),
    stats::median(measures$cd8_percent_of_cd3),
    stats::median(measures$cd161hi_percent_of_cd8),
    stats::median(measures$memory_percent_of_cd8),
    stats::median(measures$naive_percent_of_cd8),
    stats::median(measures$nk_percent_of_cd45)
  ),
  unit = c(rep("events", 5), rep("percent", 5)),
  stringsAsFactors = FALSE
)
print(hierarchy, digits = 4, row.names = FALSE)
write.csv(hierarchy, file.path(kOutputDir, "gating_hierarchy.csv"),
          row.names = FALSE)

group_medians <- do.call(rbind, lapply(
  split(measures, list(measures$severity_rank, measures$sex), drop = TRUE),
  function(piece) {
    data.frame(
      severity_rank = as.character(piece$severity_rank[1]),
      sex = piece$sex[1],
      samples = nrow(piece),
      median_cd161hi = stats::median(piece$cd161hi_percent_of_cd8, na.rm = TRUE),
      median_memory = stats::median(piece$memory_percent_of_cd8, na.rm = TRUE),
      median_naive = stats::median(piece$naive_percent_of_cd8, na.rm = TRUE),
      stringsAsFactors = FALSE
    )
  }
))
rownames(group_medians) <- NULL
write.csv(group_medians, file.path(kOutputDir, "group_medians.csv"),
          row.names = FALSE)
cat("\n=== CD161hi and memory by severity rank and sex ===\n")
print(group_medians, digits = 3, row.names = FALSE)

# ---------------------------------------------------------------------------
# Step 5: test the paper's claims
# ---------------------------------------------------------------------------

Log("Testing the claims")
claims <- ReadYuClaims(kClaimsPath)
verdicts <- TestYuClaims(claims, measures)
write.csv(verdicts, file.path(kOutputDir, "claim_verdicts.csv"),
          row.names = FALSE)

cat("\n=== The paper's claims against its own data ===\n")
print(verdicts[, c("claim_id", "short_name", "observed", "p_value", "verdict")],
      digits = 3, row.names = FALSE)
cat("\n=== Verdict counts ===\n")
print(table(verdicts$verdict))

# ---------------------------------------------------------------------------
# Step 6: what the CD45RA choice changes
# ---------------------------------------------------------------------------

Log("Sweeping the CD45RA cut across", length(kCd45raSweep), "values")

sweep_rows <- lapply(kCd45raSweep, function(ra_cut) {
  swept <- sweep_counts[sweep_counts$cd45ra_cut == ra_cut, , drop = FALSE]
  swept <- AddCd8Frequencies(swept)
  swept <- merge(swept, sheet, by.x = "sample", by.y = "file_name")

  memory_claims <- claims[claims$measure == "memory_percent_of_cd8", ]
  result <- TestYuClaims(memory_claims, swept)
  data.frame(
    cd45ra_cut = ra_cut,
    memory_percent_median = stats::median(swept$memory_percent_of_cd8,
                                          na.rm = TRUE),
    naive_percent_median = stats::median(swept$naive_percent_of_cd8,
                                         na.rm = TRUE),
    claim_4 = result$verdict[result$claim_id == 4],
    claim_6 = result$verdict[result$claim_id == 6],
    claim_8 = result$verdict[result$claim_id == 8],
    stringsAsFactors = FALSE
  )
})
sweep <- do.call(rbind, sweep_rows)
write.csv(sweep, file.path(kOutputDir, "cd45ra_sweep.csv"), row.names = FALSE)
cat("\n=== The memory claims at every CD45RA cut ===\n")
print(sweep, digits = 3, row.names = FALSE)

# ---------------------------------------------------------------------------
# Step 7: the figures the claims rest on
# ---------------------------------------------------------------------------

Log("Drawing the figures")

severity_plot <- ggplot(
  measures[!is.na(measures$cd161hi_percent_of_cd8), ],
  aes(x = severity_rank, y = cd161hi_percent_of_cd8, fill = sex)
) +
  geom_boxplot(outlier.shape = NA, alpha = 0.6) +
  geom_point(position = position_jitterdodge(jitter.width = 0.15), size = 1.4,
             alpha = 0.8) +
  labs(
    title = "CD8+ CD161hi T cells against COVID-19 severity",
    subtitle = paste(
      "FR-FCM-Z3WR, 35 marker Cytek Aurora panel, 83 samples.",
      "One cut point gated every sample."
    ),
    x = "Severity rank", y = "CD161hi, percent of CD8 T cells", fill = "Sex"
  ) +
  theme_bw()
ggsave(file.path(kOutputDir, "cd161hi_by_severity.png"), severity_plot,
       width = 9, height = 6, dpi = 150)

timepoint_levels <- c("control", "early", "middle", "late")
timepoint_data <- measures[measures$timepoint %in% timepoint_levels, ]
timepoint_data$timepoint <- factor(timepoint_data$timepoint,
                                   levels = timepoint_levels)

timepoint_plot <- ggplot(
  timepoint_data, aes(x = timepoint, y = cd161hi_percent_of_cd8, fill = sex)
) +
  geom_boxplot(outlier.shape = NA, alpha = 0.6) +
  geom_point(position = position_jitterdodge(jitter.width = 0.15), size = 1.4,
             alpha = 0.8) +
  labs(
    title = "CD8+ CD161hi T cells by time after symptom onset",
    subtitle = "Early is 14 days or less, middle is 15 to 21 days, late is over 21 days",
    x = "Time point", y = "CD161hi, percent of CD8 T cells", fill = "Sex"
  ) +
  theme_bw()
ggsave(file.path(kOutputDir, "cd161hi_by_timepoint.png"), timepoint_plot,
       width = 9, height = 6, dpi = 150)

memory_plot <- ggplot(
  timepoint_data, aes(x = timepoint, y = memory_percent_of_cd8, fill = sex)
) +
  geom_boxplot(outlier.shape = NA, alpha = 0.6) +
  geom_point(position = position_jitterdodge(jitter.width = 0.15), size = 1.4,
             alpha = 0.8) +
  labs(
    title = "CD8 memory T cells by time after symptom onset",
    subtitle = "Memory is the sum of the central memory, effector memory and EMRA gates",
    x = "Time point", y = "Memory, percent of CD8 T cells", fill = "Sex"
  ) +
  theme_bw()
ggsave(file.path(kOutputDir, "memory_by_timepoint.png"), memory_plot,
       width = 9, height = 6, dpi = 150)

slope_plot <- ggplot(
  measures[!is.na(measures$cd161hi_percent_of_cd8), ],
  aes(x = severity_index, y = cd161hi_percent_of_cd8, colour = sex)
) +
  geom_point(size = 1.8, alpha = 0.8) +
  geom_smooth(method = "lm", formula = y ~ x, se = TRUE, alpha = 0.15) +
  scale_x_continuous(
    breaks = 1:4,
    labels = c("normal", "exposed", "infected", "hospitalized")
  ) +
  labs(
    title = "The regression the paper draws in Figure 2C",
    subtitle = "A steeper fall in females is the claim under test",
    x = "Severity rank", y = "CD161hi, percent of CD8 T cells", colour = "Sex"
  ) +
  theme_bw()
ggsave(file.path(kOutputDir, "cd161hi_slope_by_sex.png"), slope_plot,
       width = 9, height = 6, dpi = 150)

sweep_long <- do.call(rbind, lapply(c("claim_4", "claim_6", "claim_8"),
  function(column) {
    data.frame(cd45ra_cut = sweep$cd45ra_cut, claim = column,
               verdict = sweep[[column]], stringsAsFactors = FALSE)
  }
))
sweep_plot <- ggplot(sweep_long, aes(x = cd45ra_cut, y = claim,
                                     fill = verdict)) +
  geom_tile(colour = "white", linewidth = 0.6) +
  labs(
    title = "What the CD45RA cut changes",
    subtitle = paste(
      "CD45RA has no density minimum inside the CD8 gate, so this cut is a",
      "choice and not a fit."
    ),
    x = "CD45RA cut, logicle scale", y = NULL, fill = "Verdict"
  ) +
  theme_bw()
ggsave(file.path(kOutputDir, "cd45ra_sweep.png"), sweep_plot,
       width = 9, height = 4, dpi = 150)

# ---------------------------------------------------------------------------
# Step 8: FlowSOM and UMAP on the pooled draw, as the paper did
# ---------------------------------------------------------------------------

if (length(all_events) < 2) {
  Log("Fewer than two samples produced events, so the clustering is skipped")
} else {
  Log("Pooling", length(all_events), "samples for the clustering")

  pooled <- do.call(rbind, all_events)
  origin <- rep(names(all_events), vapply(all_events, nrow, integer(1)))
  pooled <- RenameChannelsToMarkers(pooled, reference_frame)

  present <- intersect(kLineageMarkers, colnames(pooled))
  absent_markers <- setdiff(kLineageMarkers, colnames(pooled))
  if (length(absent_markers) > 0) {
    Log("These lineage markers are not in the panel and are skipped:",
        paste(absent_markers, collapse = ", "))
  }
  Log("Clustering", nrow(pooled), "events on", length(present), "markers")

  clustering <- RunFlowSomClustering(
    pooled, channels = present, grid_size = 10,
    n_metaclusters = kMetaclusters, seed = kSeed
  )

  embedding <- RunUmapEmbedding(pooled, channels = present, seed = kSeed)
  embedding$metacluster <- factor(clustering$metacluster)
  embedding$sample <- origin
  embedding <- merge(
    embedding, sheet[, c("file_name", "sex", "severity_rank", "timepoint")],
    by.x = "sample", by.y = "file_name", all.x = TRUE
  )

  median_expression <- ClusterMedianExpression(
    pooled, clustering$metacluster, present
  )
  write.csv(median_expression,
            file.path(kOutputDir, "cluster_median_expression.csv"),
            row.names = FALSE)

  # Clustering all leucocytes into 15 metaclusters does not isolate MAIT cells.
  # They are about 9 percent of CD8 T cells, which is under 2 percent of the
  # CD45+ parent, and a rare population is absorbed by a neighbouring
  # metacluster. The second route therefore clusters inside the CD8 gate, where
  # the subset is a tenth of the events rather than a fiftieth.
  cd8_mask <- pooled[, "CD3"] > cuts$cut[cuts$marker == "CD3"] &
    pooled[, "CD8"] > cuts$cut[cuts$marker == "CD8"]
  Log("Clustering the CD8 gate again:", sum(cd8_mask), "pooled events")

  cd8_pooled <- pooled[cd8_mask, , drop = FALSE]
  cd8_origin <- origin[cd8_mask]
  cd8_markers <- intersect(
    c("CD161", "CD45ra", "CCR7", "CD56", "CD16", "TCR gd", "HLA-DR", "CD4"),
    colnames(cd8_pooled)
  )
  cd8_clustering <- RunFlowSomClustering(
    cd8_pooled, channels = cd8_markers, grid_size = 8, n_metaclusters = 10,
    seed = kSeed
  )
  cd8_expression <- ClusterMedianExpression(
    cd8_pooled, cd8_clustering$metacluster, cd8_markers
  )
  write.csv(cd8_expression,
            file.path(kOutputDir, "cd8_cluster_median_expression.csv"),
            row.names = FALSE)
  cat("\n=== Metaclusters inside the CD8 gate ===\n")
  print(cd8_expression, digits = 3, row.names = FALSE)

  # The MAIT metacluster is the one with the highest median CD161. Nothing is
  # named by eye, and the runner up is reported so a close call is visible.
  ordered_by_cd161 <- order(cd8_expression$CD161, decreasing = TRUE)
  mait_cluster <- cd8_expression$cluster[ordered_by_cd161[1]]
  runner_up <- cd8_expression$cluster[ordered_by_cd161[2]]
  Log(sprintf(
    "The MAIT metacluster is %s with CD161 %.2f, against %.2f in the runner up %s",
    mait_cluster, cd8_expression$CD161[ordered_by_cd161[1]],
    cd8_expression$CD161[ordered_by_cd161[2]], runner_up
  ))

  cluster_frequency <- do.call(rbind, lapply(
    split(data.frame(sample = cd8_origin,
                     metacluster = cd8_clustering$metacluster,
                     stringsAsFactors = FALSE), cd8_origin),
    function(piece) {
      data.frame(
        sample = piece$sample[1],
        cd8_t_events = nrow(piece),
        mait_events = sum(piece$metacluster == mait_cluster),
        stringsAsFactors = FALSE
      )
    }
  ))
  cluster_frequency$cd161hi_percent_of_cd8 <- ifelse(
    cluster_frequency$cd8_t_events > 0,
    100 * cluster_frequency$mait_events / cluster_frequency$cd8_t_events,
    NA_real_
  )
  cluster_measures <- merge(cluster_frequency, sheet, by.x = "sample",
                            by.y = "file_name")
  write.csv(cluster_measures,
            file.path(kOutputDir, "clustering_measures.csv"), row.names = FALSE)

  cluster_claims <- claims[claims$measure == "cd161hi_percent_of_cd8", ]
  cluster_verdicts <- TestYuClaims(cluster_claims, cluster_measures)
  write.csv(cluster_verdicts,
            file.path(kOutputDir, "claim_verdicts_clustering.csv"),
            row.names = FALSE)
  cat("\n=== The same CD161hi claims, read from the clustering instead ===\n")
  print(cluster_verdicts[, c("claim_id", "short_name", "detail", "observed",
                             "p_value", "verdict")], digits = 3,
        row.names = FALSE)

  route_comparison <- merge(
    measures[, c("sample", "cd161hi_percent_of_cd8")],
    cluster_measures[, c("sample", "cd161hi_percent_of_cd8")],
    by = "sample", suffixes = c("_gating", "_clustering")
  )
  route_correlation <- stats::cor(
    route_comparison$cd161hi_percent_of_cd8_gating,
    route_comparison$cd161hi_percent_of_cd8_clustering,
    method = "spearman", use = "complete.obs"
  )
  Log(sprintf("The two routes correlate at Spearman rho = %.3f",
              route_correlation))
  write.csv(route_comparison,
            file.path(kOutputDir, "route_comparison.csv"), row.names = FALSE)

  route_plot <- ggplot(route_comparison,
                       aes(x = cd161hi_percent_of_cd8_gating,
                           y = cd161hi_percent_of_cd8_clustering)) +
    geom_point(size = 1.8, alpha = 0.8) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed",
                colour = "grey40") +
    labs(
      title = "Two routes to the same population, one sample per point",
      subtitle = sprintf(
        "Sequential gating against a FlowSOM metacluster. Spearman rho = %.3f.",
        route_correlation
      ),
      x = "CD161hi by gating, percent of CD8 T cells",
      y = "CD161hi by clustering, percent of CD8 T cells"
    ) +
    theme_bw()
  ggsave(file.path(kOutputDir, "route_comparison.png"), route_plot,
         width = 8, height = 6, dpi = 150)

  umap_plot <- ggplot(embedding, aes(x = umap_1, y = umap_2,
                                     colour = metacluster)) +
    geom_point(size = 0.15, alpha = 0.4) +
    guides(colour = guide_legend(override.aes = list(size = 3, alpha = 1))) +
    labs(
      title = paste(kMetaclusters, "FlowSOM metaclusters on the pooled draw"),
      subtitle = paste(
        format(nrow(pooled), big.mark = ","),
        "live CD45+ events,", kEventsPerSample, "per sample, as in Figure 1B"
      ),
      x = "UMAP 1", y = "UMAP 2", colour = "Metacluster"
    ) +
    theme_bw()
  ggsave(file.path(kOutputDir, "umap_metaclusters.png"), umap_plot,
         width = 9, height = 7, dpi = 150)

  umap_by_severity <- ggplot(
    embedding[!is.na(embedding$severity_rank), ], aes(x = umap_1, y = umap_2)
  ) +
    geom_point(size = 0.1, alpha = 0.3, colour = "grey30") +
    facet_grid(sex ~ severity_rank) +
    labs(
      title = "The same embedding, split by sex and severity rank",
      subtitle = "This is the layout of Figure 2A",
      x = "UMAP 1", y = "UMAP 2"
    ) +
    theme_bw()
  ggsave(file.path(kOutputDir, "umap_by_sex_and_severity.png"),
         umap_by_severity, width = 11, height = 6, dpi = 150)

  marker_long <- do.call(rbind, lapply(
    intersect(c("CD3", "CD8", "CD4", "CD56", "CD161", "CD19", "CD14", "CD45ra"),
              present),
    function(marker) {
      data.frame(
        umap_1 = embedding$umap_1, umap_2 = embedding$umap_2,
        marker = marker, value = pooled[, marker], stringsAsFactors = FALSE
      )
    }
  ))
  marker_plot <- ggplot(marker_long, aes(x = umap_1, y = umap_2,
                                         colour = value)) +
    geom_point(size = 0.1, alpha = 0.4) +
    facet_wrap(~ marker, ncol = 4) +
    scale_colour_viridis_c() +
    labs(
      title = "Marker expression on the embedding",
      subtitle = "Read the metacluster labels against these panels",
      x = "UMAP 1", y = "UMAP 2", colour = "Logicle"
    ) +
    theme_bw()
  ggsave(file.path(kOutputDir, "umap_markers.png"), marker_plot,
         width = 12, height = 7, dpi = 150)
}

Log("Done. Output is in", kOutputDir)
