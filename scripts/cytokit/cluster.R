#!/usr/bin/env Rscript

# cytokit cluster: group the events with FlowSOM, and draw them with UMAP.
#
# Clustering finds populations that a gate has to be drawn for. It does not name
# them. The naming is a separate step, cytokit annotate, because a name is a
# claim about an antibody and the scientist owns that claim.
#
# The events come either from an FCS folder or from one gate of a hierarchy that
# cytokit gate saved. Clustering inside a gate is the usual case, because a
# clustering that runs over debris and doublets spends its clusters on them.
#
# Called through cli/cytokit, never directly.

suppressPackageStartupMessages({
  library(flowCore)
  library(flowWorkspace)
  library(FlowSOM)
  library(ggplot2)
})

for (module in c("figures", "io", "panels", "compensation", "transform",
                 "clustering", "cytokit")) {
  source(file.path("R", paste0(module, ".R")))
}

arguments <- ParseCytokitArguments(
  commandArgs(trailingOnly = TRUE),
  allowed = c("data", "gates", "parent", "out", "label", "markers",
              "metaclusters", "grid", "events", "seed", "sample"),
  required = character(0),
  flags = c("recursive", "no-compensate", "no-transform", "no-umap")
)

if (is.null(arguments$data) && is.null(arguments$gates)) {
  stop("Give --data with a folder of FCS files, or --gates with a bundle that ",
       "cytokit gate wrote.\nClustering inside a gate is the usual case, ",
       "because a clustering over debris spends its clusters on debris.")
}

seed <- SetCytokitSeed(arguments)
n_metaclusters <- if (is.null(arguments$metaclusters)) 12L else
  as.integer(arguments$metaclusters)
grid_size <- if (is.null(arguments$grid)) 10L else as.integer(arguments$grid)
max_events <- if (is.null(arguments$events)) 50000L else
  as.integer(arguments$events)
sample_index <- if (is.null(arguments$sample)) 1L else
  as.integer(arguments$sample)

Say <- function(...) cat(..., "\n", sep = "")

# The events, from a saved hierarchy or from the files.
if (!is.null(arguments$gates)) {
  gates_path <- file.path(arguments$gates, "gating_set")
  if (!dir.exists(gates_path)) {
    stop("No saved hierarchy is in ", DisplayPath(arguments$gates), ".\n",
         "Point --gates at the bundle that cytokit gate wrote. It holds a ",
         "gating_set folder.")
  }
  loaded <- CollectNotes(flowWorkspace::load_gs(gates_path))
  gating_set <- loaded$value
  available <- flowWorkspace::gs_get_pop_paths(gating_set, path = "auto")
  parent <- if (is.null(arguments$parent)) {
    utils::tail(available, 1)
  } else {
    arguments$parent
  }
  if (!parent %in% available) {
    stop("The population '", parent, "' is not in the hierarchy.\n",
         "It holds: ", paste(available, collapse = ", "))
  }
  label <- if (is.null(arguments$label)) ShortLabel(arguments$gates) else
    arguments$label
  events <- ExtractGatedEvents(gating_set, parent, sample = sample_index)
  frame <- flowWorkspace::gh_pop_get_data(gating_set[[sample_index]], parent)
  source_note <- paste0("the '", parent, "' gate of ", basename(gates_path))
  inputs <- gates_path
} else {
  files <- FcsFilesIn(arguments$data, recursive = isTRUE(arguments$recursive))
  label <- if (is.null(arguments$label)) basename(arguments$data) else
    arguments$label
  read_set <- CollectNotes(read.flowSet(files[sample_index],
                                        truncate_max_range = FALSE,
                                        transformation = FALSE))
  flow_set <- read_set$value
  if (!isTRUE(arguments$`no-compensate`)) {
    state <- ReadCompensationState(flow_set[[1]])
    if (identical(state$state, "matrix to apply")) {
      applied <- CollectNotes(tryCatch(ApplyCompensation(flow_set),
                                       error = function(e) e))
      if (!inherits(applied$value, "error")) {
        flow_set <- applied$value
      }
    }
  }
  if (!isTRUE(arguments$`no-transform`)) {
    result <- CollectNotes(tryCatch(ApplyLogicleTransform(flow_set),
                                    error = function(e) e))
    if (inherits(result$value, "error")) {
      stop("The logicle transform failed: ", conditionMessage(result$value),
           "\nAdd --no-transform when the values already sit on a scale.")
    }
    flow_set <- result$value$data
  }
  frame <- flow_set[[1]]
  events <- exprs(frame)
  parent <- "every event"
  source_note <- paste0("every event of ", basename(files[sample_index]))
  inputs <- files[sample_index]
}

out_root <- if (is.null(arguments$out)) kCytokitOutputRoot else arguments$out
bundle <- OpenCytokitBundle("cluster", label, out_root)

Say("cytokit cluster")
Say("  events from ", source_note)
Say("  events      ", nrow(events))
Say("  seed        ", seed)
Say("  bundle      ", DisplayPath(bundle), "\n")

# A cluster is only as meaningful as the channel names it was built from, so the
# marker name replaces the detector name wherever the file supplies one.
events <- RenameChannelsToMarkers(events, frame)
panel <- DescribeChannels(frame)
stain_channels <- colnames(events)[!colnames(events) %in%
                                     panel$channel[!panel$is_marker]]

markers <- if (is.null(arguments$markers)) {
  stain_channels
} else {
  asked <- trimws(strsplit(arguments$markers, ",")[[1]])
  unknown <- setdiff(asked, colnames(events))
  if (length(unknown) > 0) {
    stop("These markers are not in the data: ", paste(unknown, collapse = ", "),
         "\nThe data carries: ", paste(colnames(events), collapse = ", "))
  }
  asked
}
if (length(markers) < 2) {
  stop("Clustering needs two markers. This data offers ", length(markers), ".")
}
Say("Clustering on ", length(markers), " marker(s): ",
    paste(markers, collapse = ", "))

subsampled <- SubsampleEvents(events, n = max_events, seed = seed)
drawn_from <- attr(subsampled, "sampled_from")
if (!is.null(drawn_from) && drawn_from > nrow(subsampled)) {
  Say("Subsampled ", nrow(subsampled), " of ", drawn_from, " events.")
}

clustering <- CollectNotes(RunFlowSomClustering(
  subsampled, channels = markers, grid_size = grid_size,
  n_metaclusters = n_metaclusters, seed = seed))
clusters <- clustering$value$metacluster

medians <- ClusterMedianExpression(subsampled, clusters, markers)
WriteBundleTable(bundle, medians, "cluster_medians.csv")

Say("\n", nrow(medians), " cluster(s) from a ", grid_size, " by ", grid_size,
    " grid")
print(medians[, c("cluster", "events", "percent_of_total")], row.names = FALSE,
      digits = 4)

# A cluster that holds a handful of events cannot carry a cell type label, so it
# is named here rather than left to be found in the annotation.
small <- medians[medians$events < 20, , drop = FALSE]
if (nrow(small) > 0) {
  Say("\nWARNING: ", nrow(small), " cluster(s) hold fewer than 20 events.")
  Say("  ", paste(small$cluster, collapse = ", "))
  Say("  A median over that few events is noise. Lower --metaclusters, or")
  Say("  raise --events, before you read a label from one of them.")
}

# The heat map is the picture that says what separates the clusters.
long <- do.call(rbind, lapply(markers, function(marker) {
  data.frame(cluster = factor(medians$cluster, levels = medians$cluster),
             marker = marker, value = medians[[marker]],
             stringsAsFactors = FALSE)
}))
long$scaled <- stats::ave(long$value, long$marker, FUN = function(x) {
  span <- max(x) - min(x)
  if (span == 0) rep(0.5, length(x)) else (x - min(x)) / span
})
heatmap <- ggplot2::ggplot(long,
                           ggplot2::aes(x = marker, y = cluster,
                                        fill = .data$scaled)) +
  ggplot2::geom_tile() +
  ggplot2::scale_fill_viridis_c(name = "Scaled median",
                                guide = ColourbarGuide()) +
  ggplot2::labs(title = "Median expression per cluster", x = NULL,
                y = "Cluster") +
  ThemePublication() +
  ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
SaveFigure(heatmap, file.path(bundle, "cluster_medians.svg"),
           width = 2 + 0.4 * length(markers), height = 2 + 0.3 * nrow(medians))

if (!isTRUE(arguments$`no-umap`)) {
  Say("\nRunning UMAP over ", nrow(subsampled), " events")
  embedding <- CollectNotes(RunUmapEmbedding(subsampled, channels = markers,
                                             seed = seed))
  coordinates <- as.data.frame(embedding$value)
  colnames(coordinates)[1:2] <- c("umap_1", "umap_2")
  coordinates$cluster <- factor(clusters)
  WriteBundleTable(bundle, coordinates, "umap.csv")

  # A raster, because a vector records every one of tens of thousands of points.
  drawing <- ggplot2::ggplot(coordinates,
                             ggplot2::aes(x = .data$umap_1, y = .data$umap_2,
                                          colour = .data$cluster)) +
    ggplot2::geom_point(size = 0.3, alpha = 0.6) +
    ScaleColourPublication(name = "Cluster") +
    ggplot2::labs(title = paste("UMAP of", parent), x = "UMAP 1",
                  y = "UMAP 2") +
    ThemeEmbedding() +
    LegendPoints()
  SaveFigure(drawing, file.path(bundle, "umap.png"), width = 8, height = 7)
} else {
  embedding <- list(notes = NULL)
  Say("\nUMAP was skipped, because --no-umap was given.")
}

Say("\nName these clusters with:")
Say("  cytokit annotate --clusters ", DisplayPath(bundle),
    " --definitions <file>")

notes <- unique(do.call(rbind, Filter(Negate(is.null), list(
  if (exists("loaded")) loaded$notes,
  if (exists("read_set")) read_set$notes,
  clustering$notes,
  embedding$notes
))))
ReportNotes(notes, bundle)

arguments$seed <- seed
arguments$parent <- parent
CloseCytokitBundle(
  bundle, "cluster", arguments, inputs = inputs,
  command = paste("cytokit cluster",
                  if (is.null(arguments$gates)) {
                    paste("--data", DisplayPath(arguments$data))
                  } else {
                    paste("--gates", DisplayPath(arguments$gates),
                          "--parent", parent)
                  })
)

Say("\nWrote cluster_medians.csv and the figures to ", DisplayPath(bundle))
