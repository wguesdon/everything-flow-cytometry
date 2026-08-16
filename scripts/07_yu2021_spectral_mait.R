#!/usr/bin/env Rscript

# Reproduce the flow cytometry findings of Yu 2021 from the data the authors
# deposited, on a 35 marker spectral panel.
#
# The dataset is FlowRepository FR-FCM-Z3WR, which is 83 Cytek Aurora files from
# Yu C, et al. Mucosal-associated invariant T cell responses differ by sex in
# COVID-19. Med 2021;2:755-772.e5. PMID 33870241.
#
# Three things make this different from the OMIP-039 report.
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
# The paper reads MAIT cells as CD8+ CD161hi T cells, because the flow panel
# carries no Va7.2. It validates that surrogate with scRNA-seq, where the CD161hi
# cluster co-expresses KLRB1, CD3D, CD8A and TRAV1-2. This script tests the flow
# claims only, so it inherits the surrogate and does not check it.
#
# The cohort is gated one file at a time. 23 million events at 43 parameters cost
# about 8 GB if read together, and a failure on one sample would stop the run.
#
# Run it in the container:
#   podman run --rm -v "$PWD:/work:z" -w /work everything-flow-cytometry:latest \
#     Rscript scripts/07_yu2021_spectral_mait.R

suppressPackageStartupMessages({
  library(flowCore)
  library(flowWorkspace)
  library(openCyto)
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
kTemplatePath <- file.path("gating", "yu2021_gating_template.csv")
kClaimsPath <- file.path("gating", "yu2021_paper_claims.csv")
kOutputDir <- file.path("output", "yu2021")

# The paper pooled every participant and drew 3000 live CD45+ cells per sample for
# its clustering figure. The same draw is used here so the two are comparable.
kEventsPerSample <- 3000

# Cluster on lineage markers only, so a cluster is an identity rather than a
# state. These are the markers the paper names in Figure 1C and 1E.
kLineageMarkers <- c(
  "CD3", "CD4", "CD8", "CD56", "CD19", "CD20", "CD14", "CD16",
  "CD11c", "CD123", "CD161", "TCR gd", "CD45ra", "CCR7", "HLA-DR"
)

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
missing <- sheet$file_name[!file.exists(files)]
if (length(missing) > 0) {
  stop(
    "These files are named in the sample sheet and are not on disk: ",
    paste(utils::head(missing, 5), collapse = ", "),
    ". Pull the dataset with ./sync.sh pull ",
    "datasets/flowrepository/FlowRepository_FR-FCM-Z3WR_files"
  )
}
Log("The cohort is", nrow(sheet), "files from",
    length(unique(sheet$alias_subject_id)), "subjects")

design <- as.data.frame(table(sheet$severity_rank, sheet$sex))
colnames(design) <- c("severity_rank", "sex", "samples")
write.csv(design, file.path(kOutputDir, "cohort_design.csv"), row.names = FALSE)
cat("\n=== The cohort ===\n")
print(design)

# ---------------------------------------------------------------------------
# Step 2: the panel and the transform
# ---------------------------------------------------------------------------

Log("Reading the panel from the first file")
reference_frame <- read.FCS(files[1], truncate_max_range = FALSE)
panel <- DescribeChannels(reference_frame)
write.csv(panel, file.path(kOutputDir, "panel.csv"), row.names = FALSE)

Log("The panel holds", nrow(panel), "parameters and",
    sum(panel$is_marker), "named markers")

# The viability channel carries no marker name, so FluorescenceChannels() would
# drop it. SpectralChannels() keeps it, and the first gate needs it.
transform_channels <- SpectralChannels(reference_frame)
Log("Estimating one logicle transform on", basename(files[1]), "for",
    length(transform_channels), "channels")
transform_result <- EstimateSpectralTransform(files[1], transform_channels)

template <- ReadGatingTemplate(kTemplatePath)

# ---------------------------------------------------------------------------
# Step 3: gate every file
# ---------------------------------------------------------------------------

cd56_channel <- ChannelForMarker(reference_frame, "CD56")

all_stats <- list()
all_events <- list()
all_medians <- list()
failures <- list()

for (i in seq_len(nrow(sheet))) {
  file_name <- sheet$file_name[i]
  Log(sprintf("Gating %d of %d: %s", i, nrow(sheet), file_name))

  gated <- GateOneSpectralFile(
    path = files[i],
    template = template,
    transform_list = transform_result$transform,
    subsample_population = "CD45pos",
    subsample_n = kEventsPerSample,
    median_populations = c("NK", "CD8pos"),
    median_channels = cd56_channel,
    seed = kSeed + i
  )

  if (!is.na(gated$error)) {
    Log("  FAILED:", substr(gated$error, 1, 160))
    failures[[file_name]] <- data.frame(
      file_name = file_name, message = gated$error, stringsAsFactors = FALSE
    )
    next
  }

  gated$stats$sample <- file_name
  all_stats[[file_name]] <- gated$stats

  if (!is.null(gated$events)) {
    all_events[[file_name]] <- gated$events
  }
  if (!is.null(gated$medians)) {
    gated$medians$sample <- file_name
    all_medians[[file_name]] <- gated$medians
  }
}

if (length(all_stats) == 0) {
  stop("Every file failed to gate, so there is nothing to report.")
}

stats <- do.call(rbind, all_stats)
rownames(stats) <- NULL
write.csv(stats, file.path(kOutputDir, "population_stats.csv"), row.names = FALSE)

failure_table <- if (length(failures) > 0) {
  do.call(rbind, failures)
} else {
  data.frame(file_name = character(), message = character(),
             stringsAsFactors = FALSE)
}
write.csv(failure_table, file.path(kOutputDir, "gating_failures.csv"),
          row.names = FALSE)
Log(length(all_stats), "files gated,", nrow(failure_table), "failed")

# ---------------------------------------------------------------------------
# Step 4: per sample frequencies, joined to the clinical grouping
# ---------------------------------------------------------------------------

Log("Summarising the CD8 compartment")
compartment <- SummariseCd8Compartment(stats)

if (length(all_medians) > 0) {
  medians <- do.call(rbind, all_medians)
  wide_medians <- stats::reshape(
    medians[, c("sample", "population", "median_intensity")],
    idvar = "sample", timevar = "population", direction = "wide"
  )
  colnames(wide_medians) <- sub("median_intensity.NK", "nk_cd56_median",
                               colnames(wide_medians))
  colnames(wide_medians) <- sub("median_intensity.CD8pos", "cd8_cd56_median",
                               colnames(wide_medians))
  compartment <- merge(compartment, wide_medians, by = "sample", all.x = TRUE)
}

measures <- merge(compartment, sheet, by.x = "sample", by.y = "file_name",
                  all.x = TRUE)
measures <- measures[order(measures$sample), ]
write.csv(measures, file.path(kOutputDir, "sample_measures.csv"),
          row.names = FALSE)

cat("\n=== CD161hi as a percentage of CD8 T cells, by severity rank ===\n")
print(
  SummariseByGroup <- do.call(rbind, lapply(
    split(measures, list(measures$severity_rank, measures$sex), drop = TRUE),
    function(piece) {
      data.frame(
        severity_rank = as.character(piece$severity_rank[1]),
        sex = piece$sex[1],
        n = nrow(piece),
        median_cd161hi = stats::median(piece$cd161hi_percent_of_cd8, na.rm = TRUE),
        median_memory = stats::median(piece$memory_percent_of_cd8, na.rm = TRUE),
        stringsAsFactors = FALSE
      )
    }
  )),
  digits = 3, row.names = FALSE
)
write.csv(SummariseByGroup, file.path(kOutputDir, "group_medians.csv"),
          row.names = FALSE)

# ---------------------------------------------------------------------------
# Step 5: test the paper's claims
# ---------------------------------------------------------------------------

Log("Testing the claims")
claims <- ReadYuClaims(kClaimsPath)
verdicts <- TestYuClaims(claims, measures)
write.csv(verdicts, file.path(kOutputDir, "claim_verdicts.csv"), row.names = FALSE)

cat("\n=== The paper's claims against its own data ===\n")
print(verdicts[, c("claim_id", "short_name", "expected", "observed", "p_value",
                   "verdict")], digits = 3, row.names = FALSE)

cat("\n=== Verdict counts ===\n")
print(table(verdicts$verdict))

# ---------------------------------------------------------------------------
# Step 6: the figures the claims rest on
# ---------------------------------------------------------------------------

Log("Drawing the severity and sex figures")

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
      "FR-FCM-Z3WR, 35 marker Cytek Aurora panel.",
      "Each point is one sample, gated by one template."
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
    title = "The slope the paper reports in Figure 2C",
    subtitle = "A steeper fall in females is the claim under test",
    x = "Severity rank", y = "CD161hi, percent of CD8 T cells", colour = "Sex"
  ) +
  theme_bw()
ggsave(file.path(kOutputDir, "cd161hi_slope_by_sex.png"), slope_plot,
       width = 9, height = 6, dpi = 150)

# ---------------------------------------------------------------------------
# Step 7: FlowSOM and UMAP on the pooled draw, as the paper did
# ---------------------------------------------------------------------------

if (length(all_events) < 2) {
  Log("Fewer than two samples produced events, so the clustering is skipped")
} else {
  Log("Pooling", length(all_events), "samples for the clustering")

  pooled <- do.call(rbind, all_events)
  origin <- rep(names(all_events), vapply(all_events, nrow, integer(1)))
  pooled <- RenameChannelsToMarkers(pooled, reference_frame)

  present <- intersect(kLineageMarkers, colnames(pooled))
  absent <- setdiff(kLineageMarkers, colnames(pooled))
  if (length(absent) > 0) {
    Log("These lineage markers are not in the panel and are skipped:",
        paste(absent, collapse = ", "))
  }
  Log("Clustering", nrow(pooled), "events on", length(present), "markers")

  clustering <- RunFlowSomClustering(
    pooled, channels = present, grid_size = 10,
    n_metaclusters = kMetaclusters, seed = kSeed
  )

  embedding <- RunUmapEmbedding(pooled, channels = present, seed = kSeed)
  embedding$metacluster <- factor(clustering$metacluster)
  embedding$sample <- origin

  sample_lookup <- sheet[, c("file_name", "sex", "severity_rank", "timepoint")]
  embedding <- merge(embedding, sample_lookup, by.x = "sample",
                     by.y = "file_name", all.x = TRUE)

  write.csv(
    embedding[sample.int(nrow(embedding), min(nrow(embedding), 50000)), ],
    file.path(kOutputDir, "umap_embedding_sample.csv"), row.names = FALSE
  )

  median_expression <- ClusterMedianExpression(
    pooled, clustering$metacluster, present
  )
  write.csv(median_expression,
            file.path(kOutputDir, "cluster_median_expression.csv"),
            row.names = FALSE)

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
    embedding[!is.na(embedding$severity_rank), ],
    aes(x = umap_1, y = umap_2)
  ) +
    geom_point(size = 0.1, alpha = 0.3, colour = "grey30") +
    facet_grid(sex ~ severity_rank) +
    labs(
      title = "The same embedding, split by sex and severity rank",
      subtitle = "This is the layout of Figure 2A",
      x = "UMAP 1", y = "UMAP 2"
    ) +
    theme_bw()
  ggsave(file.path(kOutputDir, "umap_by_sex_and_severity.png"), umap_by_severity,
         width = 11, height = 6, dpi = 150)

  marker_long <- do.call(rbind, lapply(
    intersect(c("CD3", "CD8", "CD4", "CD56", "CD161", "CD19", "CD14"), present),
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
