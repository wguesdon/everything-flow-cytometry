# The FlowCAP II sample classification challenges.
#
# Every other report in this repository tests a claim about a cell frequency.
# This one tests a claim about an algorithm. FlowCAP II asked whether a pipeline
# can read a class label out of a set of FCS files, it gave every entrant half of
# the labels for training, and it scored the entrants on the half it held back.
#
# Two of the three challenge datasets are in this archive.
#
#   FR-FCM-ZZZV  Challenge 3, HVTN. Gag stimulated against Env stimulated.
#   FR-FCM-ZZZU  Challenge 1, HEUvsUE. HIV exposed uninfected infants against
#                unexposed infants.
#
# A cytokine gate cannot be fitted the way every other gate in this repository is
# fitted. An unstimulated T cell is negative for all four cytokines, so the
# density carries one mode and a two component mixture splits the negative
# population down the middle. Both challenge deposits ship their own unstimulated
# controls, and the cut is therefore placed on the control and applied to the
# stimulated file of the same subject. That is how an intracellular staining
# panel is read at the bench.
#
# The classifier is deliberately the simplest one that can enter the benchmark.
# It scores every feature on the training half alone, takes the single best one,
# places a threshold on it, and applies that threshold to the held out half. A
# one feature rule cannot hide behind its own complexity, and it names the
# population it used, which is the part the paper's own discussion turns on.

# Transform every channel with one arcsinh.
#
# This is the only report in the repository that does not estimate a logicle per
# file, and the reason is the cytokine gate. A threshold taken from a subject's
# unstimulated control is applied to that subject's stimulated file, so the two
# files have to sit on one scale. A logicle estimated per file does not give
# that, and a handful of files in the HVTN deposit are refused by the logicle
# estimator at every decade count, which would leave those subjects on a
# different scale from their own controls.
#
# The arcsinh needs no fit, it never fails, and it is the standard transform for
# this kind of comparison. The cofactor is fixed and the same for every file.
TransformedEvents. <- function(frame, channels, cofactor = 150) {
  events <- flowCore::exprs(frame)
  for (channel in channels) {
    events[, channel] <- asinh(events[, channel] / cofactor)
  }
  list(events = events, transform = paste0("arcsinh/", cofactor))
}

# Match a marker name to a channel.
#
# The HEUvsUE deposit spells the same panel three ways across its 308 files. It
# writes `IL-6` and `IL6`, `IL-12` and `IL12`, `TNF-a` and `TNFa`, and `MHC II`
# and `MHCII`. Every separator is therefore removed from both sides before the
# comparison, which makes the three spellings one name.
ChannelFor. <- function(frame, marker) {
  Normalise <- function(x) toupper(gsub("[^A-Za-z0-9]", "", x))
  parameters <- flowCore::pData(flowCore::parameters(frame))
  hit <- which(Normalise(parameters$desc) == Normalise(marker))
  if (length(hit) == 0) NA_character_ else parameters$name[hit[1]]
}

#' Gate one file of the HVTN challenge down to its T cell subsets
#'
#' The hierarchy is the live cells, then CD3, then the CD4 and CD8 single
#' positive subsets. The cytokine values of each subset are returned rather than
#' cut, because the cut belongs to the subject's own unstimulated control.
#'
#' @param path Path to the FCS file.
#' @param cytokines The cytokine markers to return.
#' @return A list with `counts`, a one row `data.frame`, `cytokines`, a list of
#'   one matrix per subset, and `error`.
#' @examples
#' \dontrun{
#' GateHvtnSubsets(path)$counts
#' }
#' @export
GateHvtnSubsets <- function(path, cytokines = c("IL2", "IL4", "IFNg", "TNFa")) {
  Fail <- function(message) {
    list(counts = NULL, cytokines = NULL, error = message)
  }

  frame <- try(
    flowCore::read.FCS(path, truncate_max_range = FALSE, emptyValue = FALSE),
    silent = TRUE
  )
  if (methods::is(frame, "try-error")) {
    return(Fail(trimws(as.character(frame))))
  }

  markers <- c("ViViD", "CD3", "CD4", "CD8", cytokines)
  channels <- vapply(markers, function(m) ChannelFor.(frame, m), character(1))
  if (anyNA(channels)) {
    return(Fail(paste0("No channel carries: ",
                       paste(markers[is.na(channels)], collapse = ", "))))
  }

  transformed <- TransformedEvents.(frame, unname(channels))
  events <- transformed$events
  total <- nrow(events)

  # The viability dye does not resolve into two modes on every file of this
  # deposit. Losing a subject costs more than carrying a few dead cells into the
  # CD3 gate, so a file whose dye carries no cut keeps every event and says so.
  viability <- ResolveCut(events[, channels[["ViViD"]]])
  live <- if (viability$rule == "none") {
    events
  } else {
    events[events[, channels[["ViViD"]]] < viability$cut, , drop = FALSE]
  }
  if (nrow(live) < 1000) {
    return(Fail("Fewer than 1000 live events"))
  }

  cd3 <- ResolveCut(live[, channels[["CD3"]]])
  if (cd3$rule == "none") {
    return(Fail("No cut could be fitted on CD3"))
  }
  t_cells <- live[live[, channels[["CD3"]]] > cd3$cut, , drop = FALSE]
  if (nrow(t_cells) < 1000) {
    return(Fail("Fewer than 1000 T cells"))
  }

  cd4 <- ResolveCut(t_cells[, channels[["CD4"]]])
  cd8 <- ResolveCut(t_cells[, channels[["CD8"]]])
  if (cd4$rule == "none" || cd8$rule == "none") {
    return(Fail("No cut could be fitted on CD4 or on CD8"))
  }
  is_cd4 <- t_cells[, channels[["CD4"]]] > cd4$cut
  is_cd8 <- t_cells[, channels[["CD8"]]] > cd8$cut

  Keep <- function(mask) {
    piece <- t_cells[mask, channels[cytokines], drop = FALSE]
    colnames(piece) <- cytokines
    piece
  }
  subsets <- list(CD4 = Keep(is_cd4 & !is_cd8), CD8 = Keep(is_cd8 & !is_cd4))

  counts <- data.frame(
    file_name = basename(path),
    total_events = total,
    live_events = nrow(live),
    cd3_events = nrow(t_cells),
    cd4_events = nrow(subsets$CD4),
    cd8_events = nrow(subsets$CD8),
    CD3_percent = 100 * nrow(t_cells) / nrow(live),
    CD4_percent = 100 * nrow(subsets$CD4) / nrow(t_cells),
    CD8_percent = 100 * nrow(subsets$CD8) / nrow(t_cells),
    transform = transformed$transform,
    viability_rule = viability$rule,
    viability_cut = viability$cut,
    cd3_cut = cd3$cut,
    cd4_cut = cd4$cut,
    cd8_cut = cd8$cut,
    stringsAsFactors = FALSE
  )

  list(counts = counts, cytokines = subsets, error = NA_character_)
}

#' Gate one file of the HEUvsUE challenge down to its innate populations
#'
#' The panel carries no viability dye and no T cell marker, so the hierarchy runs
#' from a scatter debris cut to the monocytes, the myeloid dendritic cells and
#' the plasmacytoid dendritic cells.
#'
#' @param path Path to the FCS file.
#' @param cytokines The cytokine markers to return.
#' @return A list with `counts`, `cytokines` and `error`.
#' @examples
#' \dontrun{
#' GateHeuSubsets(path)$counts
#' }
#' @export
GateHeuSubsets <- function(path, cytokines = c("IFNa", "IL6", "IL12", "TNFa")) {
  Fail <- function(message) {
    list(counts = NULL, cytokines = NULL, error = message)
  }

  frame <- try(
    flowCore::read.FCS(path, truncate_max_range = FALSE, emptyValue = FALSE),
    silent = TRUE
  )
  if (methods::is(frame, "try-error")) {
    return(Fail(trimws(as.character(frame))))
  }

  markers <- c("CD14", "CD11c", "CD123", cytokines)
  channels <- vapply(markers, function(m) ChannelFor.(frame, m), character(1))
  if (anyNA(channels)) {
    return(Fail(paste0("No channel carries: ",
                       paste(markers[is.na(channels)], collapse = ", "))))
  }

  transformed <- TransformedEvents.(frame, unname(channels))
  events <- transformed$events
  total <- nrow(events)

  columns <- flowCore::colnames(frame)
  forward <- grep("^FSC-A$", columns, value = TRUE)
  if (length(forward) == 0) {
    return(Fail("The file carries no forward scatter area channel"))
  }
  debris <- DensityCut(events[, forward[1]])
  cells <- if (is.na(debris)) {
    events
  } else {
    events[events[, forward[1]] > debris, , drop = FALSE]
  }
  if (nrow(cells) < 1000) {
    return(Fail("Fewer than 1000 events above the debris cut"))
  }

  counts <- data.frame(
    file_name = basename(path),
    total_events = total,
    cell_events = nrow(cells),
    transform = transformed$transform,
    debris_cut = if (is.na(debris)) NA_real_ else debris,
    stringsAsFactors = FALSE
  )

  # CD14 and CD123 do not resolve into two modes on every file of this deposit,
  # so the whole cell population is always returned. A cytokine response can be
  # read from it even when a lineage gate cannot be placed.
  all_cells <- cells[, channels[cytokines], drop = FALSE]
  colnames(all_cells) <- cytokines
  subsets <- list(all_cells = all_cells)
  for (item in list(c("monocytes", "CD14"), c("myeloid_dc", "CD11c"),
                    c("plasmacytoid_dc", "CD123"))) {
    name <- item[1]
    marker <- item[2]
    result <- ResolveCut(cells[, channels[[marker]]])
    if (result$rule == "none") {
      counts[[paste0(name, "_percent")]] <- NA_real_
      subsets[[name]] <- NULL
      next
    }
    mask <- cells[, channels[[marker]]] > result$cut
    piece <- cells[mask, channels[cytokines], drop = FALSE]
    colnames(piece) <- cytokines
    subsets[[name]] <- piece
    counts[[paste0(name, "_percent")]] <- 100 * sum(mask) / nrow(cells)
    counts[[paste0(name, "_cut")]] <- result$cut
  }

  list(counts = counts, cytokines = subsets, error = NA_character_)
}

#' Place a cytokine cut on a control and apply it to a stimulated sample
#'
#' The threshold is a high quantile of the control distribution. Anything above
#' it in the stimulated sample is called positive, which is the rule an
#' intracellular staining panel is read by.
#'
#' @param control A numeric matrix of control events, one column per cytokine.
#' @param stimulated A numeric matrix of stimulated events, with the same
#'   columns.
#' @param quantile_level The quantile of the control that sets the threshold.
#' @param minimum The smallest population that a frequency is reported for.
#' @return A named numeric vector of percentages, one per column.
#' @examples
#' \dontrun{
#' ControlGatedFrequencies(control_matrix, stimulated_matrix)
#' }
#' @export
ControlGatedFrequencies <- function(control, stimulated,
                                    quantile_level = 0.999, minimum = 100) {
  if (is.null(control) || is.null(stimulated) ||
      nrow(control) < minimum || nrow(stimulated) < minimum) {
    columns <- if (is.null(stimulated)) character() else colnames(stimulated)
    return(stats::setNames(rep(NA_real_, length(columns)), columns))
  }
  thresholds <- apply(control, 2, stats::quantile, probs = quantile_level,
                      names = FALSE, na.rm = TRUE)
  vapply(colnames(stimulated), function(column) {
    100 * mean(stimulated[, column] > thresholds[[column]])
  }, numeric(1))
}

#' The area under the receiver operating characteristic curve
#'
#' Computed from the rank sum, which needs no package and handles ties.
#'
#' @param values A numeric vector.
#' @param labels A logical vector, `TRUE` for the positive class.
#' @return The area under the curve, or `NA` when a class is empty.
#' @examples
#' AreaUnderCurve(c(1, 2, 3, 4), c(FALSE, FALSE, TRUE, TRUE))
#' @export
AreaUnderCurve <- function(values, labels) {
  usable <- is.finite(values) & !is.na(labels)
  values <- values[usable]
  labels <- labels[usable]
  positive <- sum(labels)
  negative <- sum(!labels)
  if (positive == 0 || negative == 0) {
    return(NA_real_)
  }
  ranks <- rank(values)
  (sum(ranks[labels]) - positive * (positive + 1) / 2) / (positive * negative)
}

#' Score a classification against the truth
#'
#' @param predicted A logical vector of predictions.
#' @param actual A logical vector of the true labels.
#' @return A one row `data.frame` with `recall`, `precision`, `accuracy` and
#'   `f_measure`, which are the four measures the FlowCAP table reports.
#' @examples
#' ClassificationScore(c(TRUE, FALSE), c(TRUE, TRUE))
#' @export
ClassificationScore <- function(predicted, actual) {
  usable <- !is.na(predicted) & !is.na(actual)
  predicted <- predicted[usable]
  actual <- actual[usable]
  if (length(predicted) == 0) {
    return(data.frame(recall = NA_real_, precision = NA_real_,
                      accuracy = NA_real_, f_measure = NA_real_, samples = 0L))
  }
  true_positive <- sum(predicted & actual)
  false_positive <- sum(predicted & !actual)
  false_negative <- sum(!predicted & actual)
  recall <- if (true_positive + false_negative == 0) {
    NA_real_
  } else {
    true_positive / (true_positive + false_negative)
  }
  precision <- if (true_positive + false_positive == 0) {
    0
  } else {
    true_positive / (true_positive + false_positive)
  }
  f_measure <- if (is.na(recall) || precision + recall == 0) {
    0
  } else {
    2 * precision * recall / (precision + recall)
  }
  data.frame(
    recall = recall, precision = precision,
    accuracy = mean(predicted == actual), f_measure = f_measure,
    samples = length(predicted), stringsAsFactors = FALSE
  )
}

#' Choose one feature on the training half and score it on the held out half
#'
#' Every feature is scored by the area under the curve on the training samples,
#' the best one is taken, and a threshold is placed on it at the value that
#' maximises training accuracy. Nothing about the held out samples takes part in
#' either step.
#'
#' @param features A `data.frame` of one column per feature.
#' @param labels A logical vector, `TRUE` for the positive class.
#' @param training A logical vector, `TRUE` for a training sample.
#' @return A list with `selected`, `direction`, `threshold`, `training_score`,
#'   `testing_score` and `ranking`.
#' @examples
#' \dontrun{
#' SelectAndScore(features, labels, training)$selected
#' }
#' @export
SelectAndScore <- function(features, labels, training) {
  if (sum(training) < 4 || sum(!training) < 4) {
    stop("Both halves need at least four samples")
  }

  ranking <- do.call(rbind, lapply(names(features), function(name) {
    area <- AreaUnderCurve(features[[name]][training], labels[training])
    if (is.na(area)) {
      return(NULL)
    }
    data.frame(feature = name, training_auc = area,
               separation = abs(area - 0.5), stringsAsFactors = FALSE)
  }))
  if (is.null(ranking) || nrow(ranking) == 0) {
    stop("No feature could be scored on the training half")
  }
  ranking <- ranking[order(-ranking$separation), ]

  selected <- ranking$feature[1]
  values <- features[[selected]]
  high_is_positive <- ranking$training_auc[1] >= 0.5
  oriented <- if (high_is_positive) values else -values

  candidates <- sort(unique(oriented[training & is.finite(oriented)]))
  if (length(candidates) < 2) {
    stop("The selected feature takes fewer than two values on the training half")
  }
  midpoints <- (candidates[-1] + candidates[-length(candidates)]) / 2
  accuracy <- vapply(midpoints, function(threshold) {
    mean((oriented[training] > threshold) == labels[training], na.rm = TRUE)
  }, numeric(1))
  threshold <- midpoints[which.max(accuracy)]

  predicted <- oriented > threshold
  list(
    selected = selected,
    direction = if (high_is_positive) "high" else "low",
    threshold = if (high_is_positive) threshold else -threshold,
    training_score = ClassificationScore(predicted[training], labels[training]),
    testing_score = ClassificationScore(predicted[!training], labels[!training]),
    ranking = ranking
  )
}
