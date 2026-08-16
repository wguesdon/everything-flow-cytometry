# Unsupervised clustering and a UMAP, as a third route to the same populations.
#
# Gating asks "how many cells fall inside this boundary". Clustering asks "how
# many groups are in this data and what does each one express". The second
# question needs no gate hierarchy, so a population that nobody thought to gate
# can still appear.
#
# The two routes are compared on the same pre-gated events. Debris, dead cells
# and doublets are removed by gating first, because a cluster of debris is not an
# interesting finding and it distorts the embedding.
#
# The annotation rule is in a file, not in this code. A cluster is labelled by
# scoring its median marker expression against gating/omip39_cell_type_definitions.csv,
# so no cluster is named by eye.

#' Pull the expression matrix of one gated population out of a GatingSet
#'
#' @param gating_set A gated `GatingSet`.
#' @param population The population path to extract, for example `"Lymphocytes"`.
#' @param sample The sample name or index. Defaults to the first.
#' @return A numeric matrix of events by channels.
#' @export
ExtractGatedEvents <- function(gating_set, population, sample = 1) {
  available <- flowWorkspace::gs_get_pop_paths(gating_set, path = "auto")
  if (!population %in% available) {
    stop(
      "The population '", population, "' is not in the GatingSet. ",
      "Available: ", paste(available, collapse = ", "), "."
    )
  }

  cytoframe <- flowWorkspace::gh_pop_get_data(
    gating_set[[sample]], y = population
  )
  flowCore::exprs(cytoframe)
}

#' Take a reproducible subsample of events
#'
#' A UMAP on a million events is slow and the picture is no clearer than one on a
#' hundred thousand. The seed is an argument rather than a global, so a caller
#' cannot forget to set it.
#'
#' @param events A numeric matrix of events by channels.
#' @param n The number of events to keep. When the matrix holds fewer, every
#'   event is kept and no warning is raised.
#' @param seed The random seed.
#' @return A matrix with at most `n` rows, carrying an attribute `sampled_from`
#'   holding the original row count.
#' @export
SubsampleEvents <- function(events, n = 50000, seed = 42) {
  total <- nrow(events)
  if (total <= n) {
    attr(events, "sampled_from") <- total
    return(events)
  }

  withr::with_seed(seed, {
    keep <- sample.int(total, n)
  })

  out <- events[keep, , drop = FALSE]
  attr(out, "sampled_from") <- total
  out
}

#' Cluster events with FlowSOM and reduce to metaclusters
#'
#' FlowSOM fits a self organising map, which gives many small nodes, and then a
#' consensus step merges those nodes into a smaller number of metaclusters. The
#' metacluster is the unit that gets a cell type label.
#'
#' @param events A numeric matrix of events by channels.
#' @param channels The channel names to cluster on. Use the lineage markers, not
#'   every channel, so a cluster is defined by identity and not by activation.
#' @param grid_size The side of the square SOM grid. 10 gives 100 nodes.
#' @param n_metaclusters The number of metaclusters to merge down to.
#' @param seed The random seed.
#' @return A list with `metacluster`, an integer per event, `node`, the SOM node
#'   per event, and `fsom`, the fitted FlowSOM object.
#' @export
RunFlowSomClustering <- function(events,
                                 channels,
                                 grid_size = 10,
                                 n_metaclusters = 12,
                                 seed = 42) {
  unknown <- setdiff(channels, colnames(events))
  if (length(unknown) > 0) {
    stop("These channels are not in the data: ",
         paste(unknown, collapse = ", "), ".")
  }
  if (n_metaclusters < 2) {
    stop("n_metaclusters must be 2 or more, not ", n_metaclusters, ".")
  }

  frame <- flowCore::flowFrame(events)

  withr::with_seed(seed, {
    fsom <- FlowSOM::ReadInput(frame, transform = FALSE, scale = FALSE)
    fsom <- FlowSOM::BuildSOM(
      fsom, colsToUse = channels, xdim = grid_size, ydim = grid_size, silent = TRUE
    )
    fsom <- FlowSOM::BuildMST(fsom, silent = TRUE)
    consensus <- FlowSOM::metaClustering_consensus(
      fsom$map$codes, k = n_metaclusters, seed = seed
    )
  })

  node <- FlowSOM::GetClusters(fsom)

  list(
    metacluster = as.integer(consensus[node]),
    node = node,
    fsom = fsom
  )
}

#' Embed events in two dimensions with UMAP
#'
#' @param events A numeric matrix of events by channels.
#' @param channels The channel names to embed on. Use the same set as the
#'   clustering, so the picture and the clusters describe the same space.
#' @param n_neighbors The UMAP neighbourhood size.
#' @param min_dist The UMAP minimum distance.
#' @param seed The random seed.
#' @return A `data.frame` with the columns `umap_1` and `umap_2`.
#' @export
RunUmapEmbedding <- function(events,
                             channels,
                             n_neighbors = 15,
                             min_dist = 0.2,
                             seed = 42) {
  unknown <- setdiff(channels, colnames(events))
  if (length(unknown) > 0) {
    stop("These channels are not in the data: ",
         paste(unknown, collapse = ", "), ".")
  }

  withr::with_seed(seed, {
    embedding <- uwot::umap(
      events[, channels, drop = FALSE],
      n_neighbors = n_neighbors,
      min_dist = min_dist,
      n_components = 2,
      verbose = FALSE
    )
  })

  data.frame(umap_1 = embedding[, 1], umap_2 = embedding[, 2])
}

#' Median expression of every marker in every cluster
#'
#' @param events A numeric matrix of events by channels.
#' @param clusters An integer cluster label per event.
#' @param channels The channels to summarise.
#' @return A `data.frame` with one row per cluster and one column per channel,
#'   plus `cluster`, `events` and `percent_of_total`.
#' @export
ClusterMedianExpression <- function(events, clusters, channels) {
  if (length(clusters) != nrow(events)) {
    stop("clusters has ", length(clusters), " entries and events has ",
         nrow(events), " rows. They must match.")
  }

  levels_present <- sort(unique(clusters))
  total <- length(clusters)

  rows <- lapply(levels_present, function(cluster) {
    in_cluster <- clusters == cluster
    medians <- vapply(
      channels,
      function(channel) stats::median(events[in_cluster, channel]),
      numeric(1)
    )

    out <- as.data.frame(as.list(medians), check.names = FALSE)
    cbind(
      data.frame(
        cluster = cluster,
        events = sum(in_cluster),
        percent_of_total = 100 * sum(in_cluster) / total,
        stringsAsFactors = FALSE
      ),
      out
    )
  })

  do.call(rbind, rows)
}

#' Read the marker signature that defines each cell type
#'
#' @param path Path to a CSV whose first column is `cell_type`, whose last column
#'   is `note`, and whose remaining columns are channel names holding `pos`,
#'   `neg`, `high` or an empty string. An empty string means the marker does not
#'   take part in defining that cell type.
#' @return A `data.frame` of definitions.
#' @export
ReadCellTypeDefinitions <- function(path) {
  if (!file.exists(path)) {
    stop("The cell type definitions file does not exist: ", path)
  }

  definitions <- utils::read.csv(path, stringsAsFactors = FALSE,
                                 check.names = FALSE)
  if (!"cell_type" %in% colnames(definitions)) {
    stop("The definitions file has no 'cell_type' column.")
  }

  marker_columns <- setdiff(colnames(definitions), c("cell_type", "note"))
  values <- unlist(definitions[, marker_columns])
  allowed <- c("pos", "neg", "high", "", NA)
  bad <- setdiff(values, allowed)
  if (length(bad) > 0) {
    stop("A definition value must be pos, neg, high or empty. Found: ",
         paste(unique(bad), collapse = ", "), ".")
  }

  definitions
}

#' Label each cluster by scoring it against the cell type definitions
#'
#' Every marker's median expression is scaled across clusters to run from 0 to 1,
#' so a bright marker and a dim one contribute equally. A definition then scores
#' each cluster by adding the scaled value where it expects `pos`, subtracting it
#' where it expects `neg`, and weighting `high` twice, which separates
#' CD56bright from CD56dim. The best scoring definition wins.
#'
#' The margin between the best and the second best score is returned, because a
#' small margin means the label is a close call and a report should say so.
#'
#' @param median_expression The output of [ClusterMedianExpression()].
#' @param definitions The output of [ReadCellTypeDefinitions()].
#' @return A `data.frame` with `cluster`, `events`, `percent_of_total`,
#'   `cell_type`, `score`, `runner_up` and `margin`.
#' @export
AnnotateClusters <- function(median_expression, definitions) {
  marker_columns <- intersect(
    setdiff(colnames(definitions), c("cell_type", "note")),
    colnames(median_expression)
  )
  if (length(marker_columns) == 0) {
    stop(
      "The definitions and the expression table share no marker column. ",
      "Definitions name: ",
      paste(setdiff(colnames(definitions), c("cell_type", "note")),
            collapse = ", "), "."
    )
  }

  ScaleColumn <- function(x) {
    span <- max(x) - min(x)
    if (span == 0) rep(0.5, length(x)) else (x - min(x)) / span
  }
  scaled <- as.data.frame(lapply(median_expression[, marker_columns, drop = FALSE],
                                 ScaleColumn))
  colnames(scaled) <- marker_columns

  # The score is a weighted mean over the markers a definition specifies, not a
  # sum. A sum rewards a definition simply for naming more markers: the adaptive
  # NK definition names eight and the T cell definition names two, so under a sum
  # the adaptive label wins every cluster, including the one that is 90 percent of
  # the file and obviously T cells. Dividing by the total weight removes that.
  score_matrix <- vapply(seq_len(nrow(definitions)), function(d) {
    total <- rep(0, nrow(scaled))
    weight <- 0

    for (marker in marker_columns) {
      expectation <- definitions[d, marker]
      if (is.na(expectation) || expectation == "") next

      contribution <- switch(
        expectation,
        pos = scaled[[marker]],
        high = 2 * scaled[[marker]],
        neg = 1 - scaled[[marker]],
        NULL
      )
      if (is.null(contribution)) next

      total <- total + contribution
      weight <- weight + if (expectation == "high") 2 else 1
    }

    if (weight == 0) rep(0, nrow(scaled)) else total / weight
  }, numeric(nrow(scaled)))

  colnames(score_matrix) <- definitions$cell_type

  best <- apply(score_matrix, 1, function(row) {
    ordered <- sort(row, decreasing = TRUE)
    c(best = ordered[1], second = ordered[2])
  })

  ranking <- apply(score_matrix, 1, function(row) {
    names(sort(row, decreasing = TRUE))
  })

  data.frame(
    cluster = median_expression$cluster,
    events = median_expression$events,
    percent_of_total = median_expression$percent_of_total,
    cell_type = ranking[1, ],
    score = best[1, ],
    runner_up = ranking[2, ],
    margin = best[1, ] - best[2, ],
    stringsAsFactors = FALSE
  )
}

#' Add up the cluster frequencies for each cell type
#'
#' Several clusters can carry the same label, which is normal. The comparison
#' against a gated population needs the total.
#'
#' @param annotation The output of [AnnotateClusters()].
#' @return A `data.frame` with `cell_type`, `clusters`, `events` and
#'   `percent_of_total`, sorted by frequency.
#' @export
SummariseCellTypes <- function(annotation) {
  pieces <- split(annotation, annotation$cell_type)

  rows <- lapply(names(pieces), function(cell_type) {
    piece <- pieces[[cell_type]]
    data.frame(
      cell_type = cell_type,
      clusters = nrow(piece),
      events = sum(piece$events),
      percent_of_total = sum(piece$percent_of_total),
      stringsAsFactors = FALSE
    )
  })

  result <- do.call(rbind, rows)
  result[order(-result$percent_of_total), ]
}

#' Rename the columns of an event matrix from detectors to markers
#'
#' An FCS file names its columns after detectors, such as `yg-LP635 670/30-B-A`.
#' Every downstream file here names markers, such as `CD3 PECy5`, because a
#' detector name means nothing without the panel beside it. This translates once,
#' at the point where events leave the flow objects.
#'
#' A channel with no marker, which is scatter or time, keeps its detector name.
#' So does a channel whose marker is a placeholder. OMIP-39 labels its two unused
#' detectors `Available`, and that word is not a marker.
#'
#' @param events A numeric matrix of events by channels.
#' @param frame A `flowFrame` from the same file, used for the panel.
#' @param placeholders Marker labels that mean "no marker here". Matching ignores
#'   case.
#' @return The matrix with renamed columns.
#' @export
RenameChannelsToMarkers <- function(events, frame,
                                    placeholders = c("Available", "Empty",
                                                     "Unused", "Blank")) {
  panel <- DescribeChannels(frame)
  is_placeholder <- tolower(trimws(panel$marker)) %in% tolower(placeholders)
  panel$marker[is_placeholder] <- ""

  lookup <- stats::setNames(panel$marker, panel$channel)

  current <- colnames(events)
  replacement <- unname(lookup[current])
  keep_detector <- is.na(replacement) | !nzchar(replacement)
  replacement[keep_detector] <- current[keep_detector]

  duplicated_names <- replacement[duplicated(replacement)]
  if (length(duplicated_names) > 0) {
    stop(
      "Two channels carry the same marker name, so renaming would lose one: ",
      paste(unique(duplicated_names), collapse = ", "), "."
    )
  }

  colnames(events) <- replacement
  events
}

#' Select the clusters with the highest median expression of one marker
#'
#' A definitions file scores a cluster on the shape of its whole marker profile.
#' That is the right rule when a population is defined by a combination, and the
#' wrong one when a paper defines a population as the extreme of a single marker.
#'
#' OMIP-043 is the second case. It states that antibody secreting cells "express
#' very high levels of the ectoenzyme CD38" and that "very high CD38 expression is
#' in fact considered adequate for basic identification of ASC". Ranking clusters
#' on CD38 encodes that sentence directly.
#'
#' The difference is large where a tissue holds a second population that is
#' positive but not highest. Scoring the profile selects both; ranking selects
#' only the top. On OMIP-043 bone marrow the profile rule reaches an F1 of 12.6
#' against the manual gate and this rule reaches 90.9.
#'
#' @param median_expression The output of [ClusterMedianExpression()].
#' @param channel The channel to rank on.
#' @param n_clusters How many of the top clusters to take. Two is the default,
#'   because one alone loses events in a tissue where the population spans two
#'   clusters, and three starts admitting a positive but not extreme population.
#' @return An integer vector of cluster identifiers, ordered from highest median
#'   downwards.
#' @examples
#' \dontrun{
#' SelectByHighestMarker(median_expression, "Comp-BUV395-A", n_clusters = 2)
#' }
#' @export
SelectByHighestMarker <- function(median_expression, channel, n_clusters = 2) {
  if (!channel %in% colnames(median_expression)) {
    stop(
      "The channel '", channel, "' is not in the expression table. ",
      "It holds: ", paste(colnames(median_expression), collapse = ", "), "."
    )
  }
  if (n_clusters < 1) {
    stop("n_clusters must be 1 or more, not ", n_clusters, ".")
  }
  if (n_clusters > nrow(median_expression)) {
    stop(
      "n_clusters is ", n_clusters, " but there are only ",
      nrow(median_expression), " clusters."
    )
  }

  ordered <- order(median_expression[[channel]], decreasing = TRUE)
  median_expression$cluster[ordered[seq_len(n_clusters)]]
}
