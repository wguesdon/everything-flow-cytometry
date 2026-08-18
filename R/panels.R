# Helpers that any deposit in this repository can use.
#
# These were written for FR-FCM-ZYRN and moved here when FR-FCM-ZYN4 needed the
# same steps. Nothing in this file knows which panel it is reading.
#
# The names avoid the ones already taken elsewhere in R/. `ScatterChannels` is
# defined in R/harmonisation.R, so the one here is `PanelScatterChannels`.
# tests/testthat.R sources every file in this folder into one environment, and
# the second definition of a name silently replaces the first.

#' Name the scatter channels of a frame
#'
#' @param frame A `flowFrame`.
#' @return A named character vector with the elements `forward_area`,
#'   `forward_height` and `side_area`.
#' @examples
#' \dontrun{
#' PanelScatterChannels(frame)
#' }
#' @export
PanelScatterChannels <- function(frame) {
  available <- flowCore::colnames(frame)
  wanted <- c(forward_area = "FSC-A", forward_height = "FSC-H",
              side_area = "SSC-A")
  missing <- wanted[!wanted %in% available]
  if (length(missing) > 0) {
    stop("The frame is missing the scatter channel(s): ",
         paste(missing, collapse = ", "), ".")
  }
  wanted
}

#' Transform fluorescence values with an inverse hyperbolic sine
#'
#' One fixed cofactor is used for every channel and every file. A per file
#' estimate is not used here, because the two donors are compared to each other
#' and to a cluster assignment, and a threshold only means the same thing on
#' both files when both sit on the same scale.
#'
#' @param values A numeric matrix of events by channels.
#' @param channels The columns to transform. The rest are copied through.
#' @param cofactor The arcsinh cofactor. 150 suits fluorescence from a BD
#'   instrument that digitises to 262144.
#' @return A matrix of the same shape.
#' @examples
#' ArcsinhTransform(matrix(c(0, 150, 1500), ncol = 1,
#'                        dimnames = list(NULL, "V510-A")), "V510-A")
#' @export
ArcsinhTransform <- function(values, channels, cofactor = 150) {
  if (!is.matrix(values)) {
    stop("values must be a matrix, not a ", class(values)[1], ".")
  }
  missing <- setdiff(channels, colnames(values))
  if (length(missing) > 0) {
    stop("The matrix is missing the channel(s): ",
         paste(missing, collapse = ", "), ".")
  }
  if (cofactor <= 0) {
    stop("cofactor must be positive, not ", cofactor, ".")
  }
  values[, channels] <- asinh(values[, channels, drop = FALSE] / cofactor)
  values
}

#' Remove the doublets on the ratio of forward scatter height to area
#'
#' A doublet passes the laser for longer than a singlet, so its pulse area grows
#' while its pulse height does not. The ratio of height to area is therefore
#' lower for a doublet, and the singlets form the upper mode. The cut is a
#' robust one sided distance below the median, because the population above the
#' median is the one to keep and a symmetric rule would trim the largest
#' singlets for no reason.
#'
#' [SingletMask()] in R/spectral.R is the symmetric rule that the other reports
#' use. The names differ because `tests/testthat.R` sources every file in `R/`
#' into one environment.
#'
#' @param values A numeric matrix of events by channels.
#' @param scatter The output of [PanelScatterChannels()].
#' @param deviations How many median absolute deviations below the median the
#'   cut sits.
#' @return A logical vector, one element per event.
#' @examples
#' \dontrun{
#' RatioSingletMask(values, scatter)
#' }
#' @export
RatioSingletMask <- function(values, scatter, deviations = 3) {
  area <- values[, scatter[["forward_area"]]]
  height <- values[, scatter[["forward_height"]]]
  ratio <- height / area
  centre <- stats::median(ratio)
  spread <- stats::mad(ratio)
  if (!is.finite(spread) || spread == 0) {
    return(rep(TRUE, nrow(values)))
  }
  ratio >= centre - deviations * spread
}

#' Split a marker name into lower case tokens
#'
#' @param x A character vector of marker names.
#' @return A list of character vectors, one per input element.
#' @examples
#' MarkerTokens("TCR Va7_2 BV711")
#' @export
MarkerTokens <- function(x) {
  strsplit(NormaliseMarkerName(x), " ", fixed = TRUE)
}

#' Read the compensation state that a file records
#'
#' A deposit can carry a real spillover matrix, an identity matrix or no matrix
#' at all, and the three cases call for three different actions. This reports
#' which one a file is in, so the choice is made from the file and not from a
#' convention.
#'
#' @param frame A `flowFrame`.
#' @return A one row `data.frame` with the columns `apply_keyword`,
#'   `matrix_size`, `largest_off_diagonal` and `state`. `state` is one of
#'   `"identity, no matrix supplied"`, `"matrix to apply"` or `"no matrix"`.
#'
#' An identity matrix is reported as no matrix supplied and not as compensated
#' data. Both deposits from this laboratory store one beside `APPLY
#' COMPENSATION` set to TRUE, and on FR-FCM-ZYRN the marker correlations show
#' the values still carry their spillover.
#' @examples
#' \dontrun{
#' ReadCompensationState(frame)
#' }
#' @export
ReadCompensationState <- function(frame) {
  keywords <- flowCore::keyword(frame)
  apply_keyword <- keywords[["APPLY COMPENSATION"]]
  spillover <- keywords[["SPILL"]]
  if (is.null(spillover)) {
    spillover <- keywords[["$SPILLOVER"]]
  }

  if (is.null(spillover) || !is.matrix(spillover)) {
    return(data.frame(
      apply_keyword = if (is.null(apply_keyword)) NA_character_ else
        as.character(apply_keyword),
      matrix_size = NA_integer_, largest_off_diagonal = NA_real_,
      state = "no matrix", stringsAsFactors = FALSE
    ))
  }

  off_diagonal <- spillover
  diag(off_diagonal) <- 0
  largest <- max(abs(off_diagonal))

  data.frame(
    apply_keyword = if (is.null(apply_keyword)) NA_character_ else
      as.character(apply_keyword),
    matrix_size = nrow(spillover),
    largest_off_diagonal = largest,
    state = if (largest == 0) "identity, no matrix supplied" else
      "matrix to apply",
    stringsAsFactors = FALSE
  )
}

#' Keep the events whose scatter lies inside the digitiser range
#'
#' A file read with `truncate_max_range = FALSE` keeps the events that the
#' instrument recorded above its own range. On this deposit forward scatter runs
#' to 115422144 against a stated range of 262144. Those events are not cells and
#' they move every robust centre that is fitted afterwards.
#'
#' @param values A numeric matrix of events by channels.
#' @param scatter The output of [PanelScatterChannels()].
#' @param limit The top of the digitiser range.
#' @return A logical vector, one element per event.
#' @examples
#' InScatterRange(matrix(c(1, 3e8), ncol = 1,
#'                       dimnames = list(NULL, "FSC-A")),
#'                c(forward_area = "FSC-A"))
#' @export
InScatterRange <- function(values, scatter, limit = 262144) {
  keep <- rep(TRUE, nrow(values))
  for (channel in scatter) {
    if (!channel %in% colnames(values)) {
      next
    }
    keep <- keep & values[, channel] > 0 & values[, channel] <= limit
  }
  keep
}

#' Gate the lymphocytes on forward and side scatter
#'
#' Side scatter is cut first at the density minimum above the lowest mode, which
#' is the boundary between the lymphocytes and the monocytes. A robust bivariate
#' centre is then fitted inside that population to remove the debris that
#' remains.
#'
#' @param values A numeric matrix of events by channels.
#' @param scatter The output of [PanelScatterChannels()].
#' `robustbase::covMcd` draws random subsets, so the fit moves between runs
#' unless the generator is set. Two runs of this hierarchy differed by 1,714
#' live lymphocytes of 600,000 before the seed was added, which is small but it
#' means a number in a report cannot be reproduced exactly.
#'
#' @param quantile_limit The chi squared quantile that bounds the Mahalanobis
#'   distance from the robust centre.
#' @param seed The seed for the robust fit.
#' @return A logical vector, one element per event.
#' @examples
#' \dontrun{
#' LymphocyteMask(values, scatter)
#' }
#' @export
LymphocyteMask <- function(values, scatter, quantile_limit = 0.95,
                           seed = 42) {
  side <- values[, scatter[["side_area"]]]
  side_cut <- DensityCut(side)
  low_side <- if (is.na(side_cut)) rep(TRUE, nrow(values)) else side < side_cut
  if (sum(low_side) < 100) {
    return(low_side)
  }

  columns <- c(scatter[["forward_area"]], scatter[["side_area"]])
  pair <- values[low_side, columns, drop = FALSE]
  centre <- withr::with_seed(seed, {
    try(robustbase::covMcd(pair, alpha = 0.75), silent = TRUE)
  })
  if (methods::is(centre, "try-error")) {
    return(low_side)
  }
  distance <- try(
    stats::mahalanobis(pair, centre$center, centre$cov), silent = TRUE
  )
  if (methods::is(distance, "try-error") || anyNA(distance)) {
    return(low_side)
  }

  keep <- low_side
  keep[low_side] <- distance <= stats::qchisq(quantile_limit, df = 2)
  if (sum(keep) < 100) {
    return(low_side)
  }
  keep
}

#' Split a channel at a fitted cut and report which rule placed it
#'
#' @param values A numeric matrix of events by channels.
#' @param channel The column to cut.
#' @param keep Either `"above"` or `"below"`.
#' @return A list with `mask`, a logical vector, `cut`, the cut point, and
#'   `rule`, one of `"density"`, `"mixture"` or `"none"`.
#' @examples
#' SplitOnChannel(matrix(c(rnorm(400), rnorm(400, 6)), ncol = 1,
#'                       dimnames = list(NULL, "x")), "x")$rule
#' @export
SplitOnChannel <- function(values, channel, keep = c("above", "below")) {
  keep <- match.arg(keep)
  if (!channel %in% colnames(values)) {
    stop("The matrix has no channel called '", channel, "'.")
  }
  column <- values[, channel]
  resolved <- ResolveCut(column)
  if (is.na(resolved$cut)) {
    return(list(mask = rep(NA, length(column)), cut = NA_real_, rule = "none"))
  }
  mask <- if (keep == "above") column > resolved$cut else column <= resolved$cut
  list(mask = mask, cut = resolved$cut, rule = resolved$rule)
}

#' Measure what a one dimensional cut does to each marker of a parent
#'
#' A density cut and a two component mixture both look for a second mode. A
#' marker carried by a few percent of the parent has no second mode, so the rule
#' either returns nothing or places the cut inside the negative population and
#' calls a large fraction of the parent positive. The second failure is the
#' dangerous one, because it returns a number.
#'
#' This runs the cut on every marker named and reports the fraction it would
#' select, so the failure is a measurement in the report rather than an
#' assertion.
#'
#' @param values A numeric matrix of transformed events by channels.
#' @param channels The output of [ResolveOmip58Channels()].
#' @param parent A logical vector selecting the parent population.
#' @param markers The marker names to assess.
#' @return A `data.frame` with the columns `marker`, `channel`, `cut`, `rule`
#'   and `percent_selected`.
#' @examples
#' \dontrun{
#' AssessOneDimensionalCuts(values, channels, cd3, c("Vd1", "Vd2"))
#' }
#' @export
AssessOneDimensionalCuts <- function(values, channels, parent, markers) {
  rows <- lapply(markers, function(marker) {
    channel <- channels$channel[channels$name == marker]
    split <- SplitOnChannel(values[parent, , drop = FALSE], channel, "above")
    data.frame(
      marker = marker, channel = channel, cut = split$cut, rule = split$rule,
      percent_selected = if (split$rule == "none") NA_real_ else
        100 * mean(split$mask),
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}

#' Write one gated population to an FCS file for the Python side
#'
#' The values written are the ones the instrument recorded, which on this
#' deposit are already compensated. No transform is applied, so the Python side
#' chooses its own and the file stays readable by any cytometry tool.
#'
#' @param frame The `flowFrame` that was gated.
#' @param mask A logical vector, one element per event of `frame`.
#' @param path Where to write the file.
#' @return The number of events written, invisibly.
#' @examples
#' \dontrun{
#' WriteGatedPopulation(frame, mask, "output/omip58/handoff/donor1_cd3.fcs")
#' }
#' @export
WriteGatedPopulation <- function(frame, mask, path) {
  if (length(mask) != nrow(flowCore::exprs(frame))) {
    stop("The mask holds ", length(mask), " elements and the frame holds ",
         nrow(flowCore::exprs(frame)), " events. They must agree.")
  }
  if (sum(mask) == 0) {
    stop("The mask selects no event, so there is nothing to write to ", path,
         ".")
  }
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  flowCore::write.FCS(frame[mask, ], path)
  invisible(sum(mask))
}


#' Gate the cell cloud on forward and side scatter
#'
#' A robust bivariate centre is fitted to forward and side scatter and the
#' events far from it are dropped. Unlike [LymphocyteMask()] this does not cut
#' side scatter first, so the monocytes stay in. It is the gate a panel needs
#' when the next step separates monocytes from everything else.
#'
#' @param values A numeric matrix of events by channels.
#' @param scatter The output of [PanelScatterChannels()].
#' @param quantile_limit The chi squared quantile that bounds the Mahalanobis
#'   distance from the robust centre.
#' @param seed The seed for the robust fit, which draws random subsets.
#' @return A logical vector, one element per event.
#' @examples
#' \dontrun{
#' ScatterCloudMask(values, scatter)
#' }
#' @export
ScatterCloudMask <- function(values, scatter, quantile_limit = 0.99,
                             seed = 42) {
  columns <- c(scatter[["forward_area"]], scatter[["side_area"]])
  pair <- values[, columns, drop = FALSE]
  if (nrow(pair) < 100) {
    return(rep(TRUE, nrow(pair)))
  }
  centre <- withr::with_seed(seed, {
    try(robustbase::covMcd(pair, alpha = 0.75), silent = TRUE)
  })
  if (methods::is(centre, "try-error")) {
    return(rep(TRUE, nrow(pair)))
  }
  distance <- try(stats::mahalanobis(pair, centre$center, centre$cov),
                  silent = TRUE)
  if (methods::is(distance, "try-error") || anyNA(distance)) {
    return(rep(TRUE, nrow(pair)))
  }
  keep <- distance <= stats::qchisq(quantile_limit, df = 2)
  if (sum(keep) < 100) {
    return(rep(TRUE, nrow(pair)))
  }
  keep
}

#' Place a threshold at a percentile of an unstained control
#'
#' A marker whose positive population is a small tail has no density minimum to
#' find, and a fitted cut then lands inside the negative population. Where the
#' deposit supplies an unstained control for that channel, the control is the
#' negative population, and a high percentile of it is a threshold that owes
#' nothing to the file being gated.
#'
#' @param path Path to the unstained control.
#' @param channel The detector to read.
#' @param spillover The spillover matrix to apply first, or `NULL`.
#' @param cofactor The arcsinh cofactor.
#' @param percentile The percentile of the control to use.
#' @return The threshold on the transformed scale.
#' @examples
#' \dontrun{
#' UnstainedThreshold(path, "U450-A", spillover)
#' }
#' @export
UnstainedThreshold <- function(path, channel, spillover = NULL, cofactor = 150,
                               percentile = 0.999) {
  frame <- flowCore::read.FCS(path, truncate_max_range = FALSE)
  if (!is.null(spillover)) {
    frame <- flowCore::compensate(frame, spillover)
  }
  if (!channel %in% flowCore::colnames(frame)) {
    stop("The control has no channel called '", channel, "'.")
  }
  values <- ArcsinhTransform(flowCore::exprs(frame), channel, cofactor)
  unname(stats::quantile(values[, channel], percentile))
}

#' Place a threshold on every channel from one unstained control
#'
#' @param path Path to the unstained control.
#' @param channels A `data.frame` with the columns `name` and `channel`.
#' @param spillover The spillover matrix to apply first, or `NULL`.
#' @param cofactor The arcsinh cofactor.
#' @param percentile The percentile of the control to use.
#' @return A named numeric vector of thresholds, named by `channels$name`.
#' @examples
#' \dontrun{
#' UnstainedThresholds(path, channels, spillover)
#' }
#' @export
UnstainedThresholds <- function(path, channels, spillover = NULL,
                                cofactor = 150, percentile = 0.999) {
  frame <- flowCore::read.FCS(path, truncate_max_range = FALSE)
  if (!is.null(spillover)) {
    frame <- flowCore::compensate(frame, spillover)
  }
  missing <- setdiff(channels$channel, flowCore::colnames(frame))
  if (length(missing) > 0) {
    stop("The control is missing the channel(s): ",
         paste(missing, collapse = ", "), ".")
  }
  values <- ArcsinhTransform(flowCore::exprs(frame), channels$channel, cofactor)
  stats::setNames(
    vapply(channels$channel,
           function(channel) unname(stats::quantile(values[, channel],
                                                    percentile)),
           numeric(1)),
    channels$name
  )
}

#' Draw a gate hierarchy from a table of counts
#'
#' A `GatingSet` draws its own tree, but a hierarchy built from event masks has
#' no object to ask. This lays the tree out from the `population` and `parent`
#' columns, so any analysis that records those two can show its flow.
#'
#' @param counts A `data.frame` with the columns `population`, `parent`,
#'   `events` and `percent_of_parent`. The root carries `NA` as its parent.
#' @param title The plot title.
#' @return A `ggplot` object.
#' @examples
#' \dontrun{
#' PlotGateTree(counts)
#' }
#' @export
PlotGateTree <- function(counts, title = "Gate hierarchy") {
  required <- c("population", "parent", "events", "percent_of_parent")
  missing <- setdiff(required, colnames(counts))
  if (length(missing) > 0) {
    stop("counts is missing the column(s): ", paste(missing, collapse = ", "),
         ".")
  }

  depth <- rep(NA_integer_, nrow(counts))
  depth[is.na(counts$parent)] <- 0L
  for (pass in seq_len(nrow(counts))) {
    for (index in which(is.na(depth))) {
      parent_at <- match(counts$parent[index], counts$population)
      if (!is.na(parent_at) && !is.na(depth[parent_at])) {
        depth[index] <- depth[parent_at] + 1L
      }
    }
    if (!anyNA(depth)) {
      break
    }
  }
  if (anyNA(depth)) {
    stop("These populations name a parent that is not in the table: ",
         paste(counts$population[is.na(depth)], collapse = ", "), ".")
  }

  # A row is placed by a depth first walk from the root rather than by depth,
  # so that a child sits under its own parent. Ordering by depth alone
  # interleaves the branches, and an edge then runs past an unrelated node that
  # a reader takes for the parent.
  VisitOrder. <- function(population) {
    children <- counts$population[!is.na(counts$parent) &
                                    counts$parent == population]
    c(population, unlist(lapply(children, VisitOrder.), use.names = FALSE))
  }
  roots <- counts$population[is.na(counts$parent)]
  walk <- unlist(lapply(roots, VisitOrder.), use.names = FALSE)
  position <- stats::setNames(seq_along(walk), walk)
  nodes <- data.frame(
    population = counts$population,
    depth = depth,
    y = -unname(position[counts$population]),
    label = sprintf("%s\n%s (%.1f%%)", counts$population,
                    format(counts$events, big.mark = ","),
                    counts$percent_of_parent),
    stringsAsFactors = FALSE
  )

  edges <- counts[!is.na(counts$parent), c("population", "parent")]
  edges$x <- nodes$depth[match(edges$parent, nodes$population)]
  edges$y <- nodes$y[match(edges$parent, nodes$population)]
  edges$xend <- nodes$depth[match(edges$population, nodes$population)]
  edges$yend <- nodes$y[match(edges$population, nodes$population)]

  ggplot2::ggplot() +
    ggplot2::geom_segment(
      data = edges,
      ggplot2::aes(x = .data$x, y = .data$y, xend = .data$xend,
                   yend = .data$yend),
      colour = "grey65"
    ) +
    ggplot2::geom_label(
      data = nodes,
      ggplot2::aes(x = .data$depth, y = .data$y, label = .data$label),
      hjust = 0, size = 2.9, linewidth = 0.25, fill = "white",
      colour = "grey15", family = FigureFont(), lineheight = 1.05
    ) +
    ggplot2::scale_x_continuous(expand = ggplot2::expansion(c(0.02, 0.25))) +
    ggplot2::labs(title = title, x = NULL, y = NULL) +
    ggplot2::theme_void(base_size = 11, base_family = FigureFont()) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(size = 12, hjust = 0,
                                         margin = ggplot2::margin(b = 8)),
      plot.margin = ggplot2::margin(8, 8, 8, 8)
    )
}

#' Draw one gate as a two dimensional plot with its thresholds
#'
#' @param values A numeric matrix of transformed events by channels.
#' @param parent A logical vector selecting the parent population.
#' @param x The channel on the horizontal axis.
#' @param y The channel on the vertical axis.
#' @param x_label The label for the horizontal axis.
#' @param y_label The label for the vertical axis.
#' @param x_threshold A vertical line to draw, or `NA`.
#' @param y_threshold A horizontal line to draw, or `NA`.
#' @param title The panel title.
#' @param max_points How many events to draw. The rest are dropped, because a
#'   panel of a million points is a black rectangle and a slow render.
#' @param seed The seed for the subsample.
#' @return A `ggplot` object.
#' @examples
#' \dontrun{
#' PlotGatePair(values, parent, "FSC-A", "SSC-A")
#' }
#' @export
PlotGatePair <- function(values, parent, x, y, x_label = x, y_label = y,
                         x_threshold = NA_real_, y_threshold = NA_real_,
                         title = NULL, max_points = 40000, seed = 42) {
  for (channel in c(x, y)) {
    if (!channel %in% colnames(values)) {
      stop("The matrix has no channel called '", channel, "'.")
    }
  }
  rows <- which(parent)
  if (length(rows) > max_points) {
    rows <- withr::with_seed(seed, sample(rows, max_points))
  }
  frame <- data.frame(x = values[rows, x], y = values[rows, y])

  plot <- PlotDensityScatter(frame, "x", "y", x_label = x_label,
                             y_label = y_label, title = title,
                             point_size = 0.2) +
    ggplot2::scale_x_continuous(labels = AxisLabels) +
    ggplot2::scale_y_continuous(labels = AxisLabels)

  # The line is drawn over the events so that a reader can see which side of it
  # a dense population sits on.
  if (!is.na(x_threshold)) {
    plot <- plot + ggplot2::geom_vline(xintercept = x_threshold,
                                       colour = "#D55E00", linewidth = 0.6)
  }
  if (!is.na(y_threshold)) {
    plot <- plot + ggplot2::geom_hline(yintercept = y_threshold,
                                       colour = "#D55E00", linewidth = 0.6)
  }
  plot
}
