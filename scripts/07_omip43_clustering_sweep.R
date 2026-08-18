#!/usr/bin/env Rscript

# Improving the clustering route on OMIP-043.
#
# The first attempt recovered antibody secreting cells well in tonsil and spleen,
# loosely in blood, and not at all in bone marrow. Bone marrow holds the fewest
# ASC of the four tissues, so the working explanation is that a rare population
# cannot claim a metacluster of its own when there are only fifteen to go round.
#
# Three changes are worth testing, and this script measures all of them rather
# than picking one.
#
#   Side scatter. The paper states that ASC carry more cytoplasm and sit above a
#   normal lymphocyte scatter gate, and Part 5 of the report confirms a ratio of
#   1.8 to 2.7. That is a strong discriminating feature and the first attempt did
#   not give it to the clustering at all.
#
#   More metaclusters. With fifteen groups over 44,000 bone marrow events, the
#   average group holds 2,900 events and the target population holds 195. It
#   cannot dominate anything.
#
#   Scaling. FlowSOM weights each channel by its spread. A marker with a wide
#   range then pulls harder than one with a narrow range, whether or not it
#   separates anything.
#
# The manual gate is used here to score configurations, which makes this a
# calibration step and not an independent test. The report must say so. Once a
# configuration is chosen, script 06 runs it and the independent comparison there
# stands or falls on its own.
#
# Run it in the container:
#   podman run --rm -v "$PWD:/work:z" -w /work everything-flow-cytometry:latest \
#     Rscript scripts/07_omip43_clustering_sweep.R

suppressPackageStartupMessages({
  library(flowCore)
  library(flowWorkspace)
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
kDefinitionsPath <- file.path("gating", "omip43_cell_type_definitions.csv")
kOutputDir <- file.path("output", "omip43")

kTissues <- c("PBMC", "Bone Marrow", "Spleen", "Tonsil")
kParentPath <- "/time/live/scatter/sing1/sing2/dumped"
kAscPath <- paste0(kParentPath, "/PC")

kMarkerChannels <- c(
  "Comp-APC-A", "Comp-BV421-A", "Comp-PerCP-Cy5-5-A",
  "Comp-BV786-A", "Comp-BUV395-A", "Comp-FITC-A", "Comp-PE-A"
)
kScatterChannels <- c("SSC-A", "FSC-A")
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

#' Score a selected population against the manual gate
#'
#' @param selected A logical vector, `TRUE` where the method calls an event ASC.
#' @param truth A logical vector, `TRUE` where the manual gate does.
#' @return A one row `data.frame` with `precision`, `recall` and `f1`, all as
#'   percentages. `f1` is the harmonic mean, which punishes a method that wins on
#'   one at the cost of the other.
ScoreSelection <- function(selected, truth) {
  true_positive <- sum(selected & truth)
  precision <- if (sum(selected) == 0) 0 else 100 * true_positive / sum(selected)
  recall <- if (sum(truth) == 0) NA_real_ else 100 * true_positive / sum(truth)
  f1 <- if (is.na(recall) || precision + recall == 0) {
    0
  } else {
    2 * precision * recall / (precision + recall)
  }

  data.frame(precision = precision, recall = recall, f1 = f1,
             selected_events = sum(selected), stringsAsFactors = FALSE)
}

Log("Importing the workspace")
workspace <- CytoML::open_flowjo_xml(kWorkspace)
definitions <- ReadCellTypeDefinitions(kDefinitionsPath)

# Read every tissue once, so the sweep does not repeat the slow part.
tissue_data <- list()
for (tissue in kTissues) {
  gating_set <- CytoML::flowjo_to_gatingset(workspace, name = tissue, path = kDataDir)
  gating_hierarchy <- gating_set[[1]]

  events_all <- flowCore::exprs(
    flowWorkspace::gh_pop_get_data(gating_hierarchy, kParentPath)
  )
  truth_all <- flowWorkspace::gh_pop_get_indices(gating_hierarchy, kAscPath)[
    flowWorkspace::gh_pop_get_indices(gating_hierarchy, kParentPath)
  ]

  withr::with_seed(kSeed, {
    keep <- sample.int(nrow(events_all), min(kSubsampleSize, nrow(events_all)))
  })

  tissue_data[[tissue]] <- list(
    events = events_all[keep, , drop = FALSE],
    truth = truth_all[keep]
  )
  Log(" ", tissue, nrow(tissue_data[[tissue]]$events), "events,",
      sum(tissue_data[[tissue]]$truth), "ASC",
      sprintf("(%.2f%%)", 100 * mean(tissue_data[[tissue]]$truth)))
}

# The selection rule is the fourth dimension and it turned out to matter more than
# the other three together. "profile" scores each cluster against the definitions
# file and takes every cluster labelled ASC. "cd38 top N" ranks clusters on median
# CD38 and takes the highest, which is what the paper's own identification
# criterion says: "very high CD38 expression is in fact considered adequate for
# basic identification of ASC".
settings <- expand.grid(
  metaclusters = c(15, 25, 40),
  use_scatter = c(FALSE, TRUE),
  scale_channels = c(FALSE, TRUE),
  selection = c("profile", "cd38 top 1", "cd38 top 2", "cd38 top 3"),
  stringsAsFactors = FALSE
)

Log("Running", nrow(settings), "configurations across", length(kTissues), "tissues")

sweep_rows <- list()
for (i in seq_len(nrow(settings))) {
  setting <- settings[i, ]
  channels <- if (setting$use_scatter) {
    c(kMarkerChannels, kScatterChannels)
  } else {
    kMarkerChannels
  }

  for (tissue in kTissues) {
    events <- tissue_data[[tissue]]$events
    truth <- tissue_data[[tissue]]$truth

    outcome <- tryCatch({
      frame <- flowCore::flowFrame(events)
      withr::with_seed(kSeed, {
        fsom <- FlowSOM::ReadInput(frame, transform = FALSE,
                                   scale = setting$scale_channels)
        fsom <- FlowSOM::BuildSOM(fsom, colsToUse = channels,
                                  xdim = 10, ydim = 10, silent = TRUE)
        fsom <- FlowSOM::BuildMST(fsom, silent = TRUE)
        consensus <- FlowSOM::metaClustering_consensus(
          fsom$map$codes, k = setting$metaclusters, seed = kSeed
        )
      })
      metacluster <- as.integer(consensus[FlowSOM::GetClusters(fsom)])

      medians <- ClusterMedianExpression(events, metacluster, channels)

      asc_clusters <- if (setting$selection == "profile") {
        labels <- AnnotateClusters(medians, definitions)
        labels$cluster[labels$cell_type == "Antibody secreting cells"]
      } else {
        top_n <- as.integer(sub("cd38 top ", "", setting$selection))
        SelectByHighestMarker(medians, "Comp-BUV395-A", n_clusters = top_n)
      }
      selected <- metacluster %in% asc_clusters

      ScoreSelection(selected, truth)
    }, error = function(e) {
      data.frame(precision = NA_real_, recall = NA_real_, f1 = NA_real_,
                 selected_events = NA_integer_, stringsAsFactors = FALSE)
    })

    sweep_rows[[length(sweep_rows) + 1]] <- cbind(
      data.frame(tissue = tissue, stringsAsFactors = FALSE),
      setting,
      outcome
    )
  }
  Log(sprintf("  metaclusters=%d scatter=%s scaled=%s selection=%s done",
              setting$metaclusters, setting$use_scatter, setting$scale_channels,
              setting$selection))
}

sweep <- do.call(rbind, sweep_rows)
write.csv(sweep, file.path(kOutputDir, "clustering_sweep.csv"), row.names = FALSE)

# A configuration is only useful if it works on every tissue, so rank on the
# worst tissue rather than the average. An average hides a total failure.
summary_rows <- lapply(seq_len(nrow(settings)), function(i) {
  setting <- settings[i, ]
  piece <- sweep[sweep$metaclusters == setting$metaclusters &
                   sweep$use_scatter == setting$use_scatter &
                   sweep$scale_channels == setting$scale_channels &
                   sweep$selection == setting$selection, ]
  cbind(setting, data.frame(
    mean_f1 = mean(piece$f1, na.rm = TRUE),
    worst_f1 = min(piece$f1, na.rm = TRUE),
    worst_tissue = piece$tissue[which.min(piece$f1)],
    stringsAsFactors = FALSE
  ))
})
sweep_summary <- do.call(rbind, summary_rows)
sweep_summary <- sweep_summary[order(-sweep_summary$worst_f1), ]
write.csv(sweep_summary, file.path(kOutputDir, "clustering_sweep_summary.csv"),
          row.names = FALSE)

cat("\n=== Configurations, ranked by their worst tissue ===\n")
print(sweep_summary, digits = 3, row.names = FALSE)

best <- sweep_summary[1, ]
Log(sprintf("Best: metaclusters=%d scatter=%s scaled=%s selection=%s, worst tissue F1 %.1f",
            best$metaclusters, best$use_scatter, best$scale_channels,
            best$selection, best$worst_f1))

cat("\n=== Per tissue, first configuration against the best ===\n")
first <- sweep[sweep$metaclusters == 15 & !sweep$use_scatter &
                 !sweep$scale_channels & sweep$selection == "profile", ]
chosen <- sweep[sweep$metaclusters == best$metaclusters &
                  sweep$use_scatter == best$use_scatter &
                  sweep$scale_channels == best$scale_channels &
                  sweep$selection == best$selection, ]
comparison <- data.frame(
  tissue = first$tissue,
  first_precision = first$precision,
  first_recall = first$recall,
  first_f1 = first$f1,
  best_precision = chosen$precision[match(first$tissue, chosen$tissue)],
  best_recall = chosen$recall[match(first$tissue, chosen$tissue)],
  best_f1 = chosen$f1[match(first$tissue, chosen$tissue)],
  stringsAsFactors = FALSE
)
write.csv(comparison, file.path(kOutputDir, "clustering_improvement.csv"),
          row.names = FALSE)
print(comparison, digits = 3, row.names = FALSE)

sweep$configuration <- sprintf("%s, %d clusters%s%s", sweep$selection,
                               sweep$metaclusters,
                               ifelse(sweep$use_scatter, ", scatter", ""),
                               ifelse(sweep$scale_channels, ", scaled", ""))

sweep_plot <- ggplot(sweep, aes(x = reorder(configuration, f1), y = f1,
                                fill = tissue)) +
  geom_col(position = position_dodge(width = 0.85), width = 0.8) +
  coord_flip() +
  ScaleFillPublication(name = "Tissue") +
  labs(
    title = "Recovering antibody secreting cells by clustering",
    subtitle = paste(
      "F1 against the manual gate. A configuration is only usable if it works",
      "on every tissue."
    ),
    x = NULL, y = "F1, percent"
  ) +
  ThemePublication()
SaveFigure(sweep_plot, file.path(kOutputDir, "clustering_sweep.svg"),
  width = 10, height = 8)

Log("Done. Output is in", kOutputDir)
