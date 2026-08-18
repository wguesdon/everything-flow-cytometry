# The naive and memory T cell hierarchy of the FR-FCM-Z282 study.
#
# The study SOP gates in this order. It takes the window of stable acquisition,
# then singlets, then live cells for frozen PBMC or CD45 positive leukocytes for
# whole blood, then lymphocytes on scatter, then CD3 positive cells, then the
# CD4
# and CD8 single positive cells inside CD3, and finally the four naive and
# memory
# quadrants on CD45RA against CCR7 inside CD3, CD4 and CD8.
#
# That is fifteen reported frequencies. CD3 is a percentage of lymphocytes, CD4
# and CD8 are percentages of CD3, and the twelve quadrant values are percentages
# of their own parent.
#
# Two decisions are recorded here rather than buried in the code.
#
# The first is the cut rule. `DensityCut()` finds the deepest minimum between
# the
# two tallest modes of the kernel density. When a marker has no such minimum the
# function returns `NA` and the caller records a failure. It does not fall back
# on a quantile, because a silent fallback is what made the OMIP-043 template
# report a number for a gate that had kept everything.
#
# The second is where the CD45RA and CCR7 cuts are fitted. The paper describes
# both markers as continuous and poorly resolved, and the OMIP-039 and Yu 2021
# reports met the same problem. Both cuts are therefore fitted once per file
# inside the CD3 gate, which holds the most events, and the same two cuts gate
# the CD4 and the CD8 compartments of that file. Within one file the boundary
# between positive and negative for one reagent does not move between two T cell
# subsets.

#' Find the density minimum that separates two modes
#'
#' The rule is the one that `openCyto::gate_mindensity` implements, written out
#' here so that the failure is visible. The function estimates a kernel density,
#' finds every local maximum and every local minimum, takes the two tallest
#' maxima, and returns the lowest minimum that lies between them.
#'
#' One bandwidth does not fit this deposit. CCR7 and CD45RA are continuous and
#' poorly resolved, and the width that resolves two modes on one instrument
#' smooths them into one on another. The function therefore walks a fixed ladder
#' of bandwidths from the smoothest to the roughest and returns the first cut
#' that separates two modes by at least `depth`. The ladder is fixed and it is
#' the same for every file, so the rule is still one rule.
#'
#' @param values A numeric vector, already transformed.
#' @param adjust The ladder of bandwidth multipliers passed to
#'   [stats::density()], from the smoothest to the roughest.
#' @param trim The quantile trimmed from each tail before the density is
#'   estimated, so a handful of extreme events cannot create a mode.
#' @param depth How far the valley must sit below the shorter of the two modes,
#'   as a fraction of that mode. A ripple on a shoulder does not qualify.
#' @return The cut point, or `NA_real_` when no bandwidth in the ladder resolves
#'   two modes.
#' @examples
#' DensityCut(c(rnorm(500, 0), rnorm(500, 6)))
#' @export
DensityCut <- function(values,
                       adjust = c(2, 1.5, 1, 0.75, 0.5),
                       trim = 0.001,
                       depth = 0.05) {
  values <- values[is.finite(values)]
  if (length(values) < 100) {
    return(NA_real_)
  }
  limits <- stats::quantile(values, c(trim, 1 - trim), names = FALSE)
  if (!all(is.finite(limits)) || diff(limits) <= 0) {
    return(NA_real_)
  }
  values <- values[values >= limits[1] & values <= limits[2]]
  if (length(values) < 100) {
    return(NA_real_)
  }

  for (width in adjust) {
    cut <- DensityCutAt.(values, width, depth)
    if (!is.na(cut)) {
      return(cut)
    }
  }
  NA_real_
}

# One bandwidth of the ladder. Returns NA when this width resolves fewer than
# two
# modes, or when the valley between the two tallest is not deep enough.
DensityCutAt. <- function(values, adjust, depth) {
  estimate <- stats::density(values, adjust = adjust, n = 512)
  y <- estimate$y
  inner <- seq.int(2, length(y) - 1)
  peaks <- inner[y[inner] > y[inner - 1] & y[inner] > y[inner + 1]]
  valleys <- inner[y[inner] < y[inner - 1] & y[inner] < y[inner + 1]]
  if (length(peaks) < 2 || length(valleys) == 0) {
    return(NA_real_)
  }

  tallest <- peaks[order(y[peaks], decreasing = TRUE)][1:2]
  left <- min(tallest)
  right <- max(tallest)
  between <- valleys[valleys > left & valleys < right]
  if (length(between) == 0) {
    return(NA_real_)
  }
  lowest <- between[which.min(y[between])]
  if (y[lowest] > (1 - depth) * min(y[tallest])) {
    return(NA_real_)
  }
  estimate$x[lowest]
}

#' Split a marker with a two component mixture
#'
#' CCR7 and CD45RA are continuous. The paper says so in its own discussion, and
#' [DensityCut()] returns `NA` on many of the files for that reason. A mixture
#' does not need a valley. It fits a negative component and a positive component
#' and cuts where the two weighted densities cross, which is the boundary that
#' minimises the events assigned to the wrong component.
#'
#' The fit is an expectation maximisation loop written out here rather than
#' taken
#' from a package. `mclust::Mclust` needs its own package attached to the search
#' path before it runs, which makes a function that calls it fail in a test
#' session, and the model is twenty lines of arithmetic.
#'
#' @param values A numeric vector, already transformed.
#' @param max_events The largest number of events to fit. A larger sample is
#'   drawn down to this size, because the fit converges long before it.
#' @param seed The random seed for that draw.
#' @param iterations The largest number of expectation maximisation steps.
#' @param tolerance The change in the log likelihood that ends the loop.
#' @return The cut point, or `NA_real_` when the fit fails.
#' @examples
#' MixtureCut(c(rnorm(2000, 0), rnorm(2000, 3)))
#' @export
MixtureCut <- function(values, max_events = 20000, seed = 42,
                       iterations = 300, tolerance = 1e-8) {
  values <- values[is.finite(values)]
  if (length(values) < 100 || stats::var(values) == 0) {
    return(NA_real_)
  }
  if (length(values) > max_events) {
    withr::with_seed(seed, {
      values <- values[sample.int(length(values), max_events)]
    })
  }

  fit <- FitTwoComponents.(values, iterations, tolerance)
  if (is.null(fit)) {
    return(NA_real_)
  }

  grid <- seq(fit$means[1], fit$means[2], length.out = 1024)
  low <- fit$weights[1] * stats::dnorm(grid, fit$means[1],
                     sqrt(fit$variances[1]))
  high <- fit$weights[2] * stats::dnorm(grid, fit$means[2],
                      sqrt(fit$variances[2]))
  crossing <- which(diff(sign(low - high)) != 0)
  if (length(crossing) == 0) {
    return(NA_real_)
  }
  grid[crossing[1]]
}

# Fit a two component Gaussian mixture by expectation maximisation.
#
# The start is deterministic. The median splits the values, and each half gives
# the mean and the variance of one component. There is no random start, so the
# same input always gives the same cut.
#
# Returns NULL when the two components do not separate, which happens when one
# component collapses onto a single value.
FitTwoComponents. <- function(values, iterations, tolerance) {
  middle <- stats::median(values)
  lower <- values[values <= middle]
  upper <- values[values > middle]
  if (length(lower) < 10 || length(upper) < 10) {
    return(NULL)
  }

  means <- c(mean(lower), mean(upper))
  variances <- c(stats::var(lower), stats::var(upper))
  weights <- c(length(lower), length(upper)) / length(values)
  floor_variance <- stats::var(values) * 1e-6
  variances <- pmax(variances, floor_variance)

  previous <- -Inf
  for (step in seq_len(iterations)) {
    low <- weights[1] * stats::dnorm(values, means[1], sqrt(variances[1]))
    high <- weights[2] * stats::dnorm(values, means[2], sqrt(variances[2]))
    total <- low + high
    if (any(!is.finite(total)) || all(total == 0)) {
      return(NULL)
    }
    total[total == 0] <- .Machine$double.xmin

    responsibility <- low / total
    counts <- c(sum(responsibility), sum(1 - responsibility))
    if (min(counts) < 10) {
      return(NULL)
    }

    weights <- counts / length(values)
    means <- c(
      sum(responsibility * values) / counts[1],
      sum((1 - responsibility) * values) / counts[2]
    )
    variances <- pmax(c(
      sum(responsibility * (values - means[1])^2) / counts[1],
      sum((1 - responsibility) * (values - means[2])^2) / counts[2]
    ), floor_variance)

    likelihood <- sum(log(total))
    if (is.finite(likelihood) && abs(likelihood - previous) < tolerance) {
      break
    }
    previous <- likelihood
  }

  order_index <- order(means)
  means <- means[order_index]
  variances <- variances[order_index]
  weights <- weights[order_index]
  if (!all(is.finite(c(means, variances, weights))) || diff(means) <= 0) {
    return(NULL)
  }
  list(means = means, variances = variances, weights = weights)
}

#' Cut a marker by the density minimum, and fall back on a mixture
#'
#' @param values A numeric vector, already transformed.
#' @return A list with `cut`, the cut point, and `rule`, either `"density"`,
#'   `"mixture"` or `"none"`.
#' @examples
#' ResolveCut(c(rnorm(500, 0), rnorm(500, 6)))$rule
#' @export
ResolveCut <- function(values) {
  cut <- DensityCut(values)
  if (!is.na(cut)) {
    return(list(cut = cut, rule = "density"))
  }
  cut <- MixtureCut(values)
  if (!is.na(cut)) {
    return(list(cut = cut, rule = "mixture"))
  }
  list(cut = NA_real_, rule = "none")
}

# Return the fraction of a logical vector that is TRUE, or NA when the parent is
# empty. A zero event parent must not report zero percent, because that reads as
# a measurement rather than as a missing one.
Fraction. <- function(mask) {
  if (length(mask) == 0) {
    return(NA_real_)
  }
  100 * sum(mask) / length(mask)
}

#' Gate one file of the FR-FCM-Z282 deposit
#'
#' Every cut is fitted on the file itself, so the function is one automated
#' analyst applied to all 234 files. It never reads the operator's own result
#' and
#' it never reads the reference operator's result.
#'
#' @param path Path to the FCS file.
#' @param panel The table returned by [ReadZ282Panel()].
#' @param material Either `"PBMC"` or `"WB"`.
#' @param time_bins The number of time bins for [StableTimeWindow()].
#' @param fixed_cuts An optional named numeric vector of cut points on the
#'   logicle scale, with the names `NIR` or `CD45`, `CD3`, `CD4`, `CD8`,
#'   `CD45RA` and `CCR7`. A marker named here is not fitted. Supplying one cut
#'   per operator turns the function from thirteen automated analysts into one,
#'   which is what the reference operator does by hand.
#' @return A list with `counts`, a one row `data.frame` holding the fifteen
#'   frequencies and the event count at every step, `cuts`, a `data.frame` of
#' the
#'   cut point of each marker, `channels`, the output of
#'   [ResolveMarkerChannels()], and `error`, a message or `NA`.
#' @examples
#' \dontrun{
#' GateNaiveMemoryFile(path, panel, "WB")$counts
#' }
#' @export
GateNaiveMemoryFile <- function(path, panel, material, time_bins = 100,
                                fixed_cuts = NULL) {
  Fail <- function(message) {
    list(counts = NULL, cuts = NULL, channels = NULL, error = message)
  }

  frame <- try(
    flowCore::read.FCS(path, truncate_max_range = FALSE, emptyValue = FALSE),
    silent = TRUE
  )
  if (methods::is(frame, "try-error")) {
    return(Fail(trimws(as.character(frame))))
  }

  channels <- try(ResolveMarkerChannels(frame, panel, material), silent = TRUE)
  if (methods::is(channels, "try-error")) {
    return(Fail(trimws(as.character(channels))))
  }
  if (anyNA(channels$channel)) {
    return(Fail(paste0(
      "No channel resolves for: ",
      paste(channels$marker[is.na(channels$channel)], collapse = ", ")
    )))
  }

  scatter <- try(ScatterChannels(frame), silent = TRUE)
  if (methods::is(scatter, "try-error")) {
    return(Fail(trimws(as.character(scatter))))
  }

  settled <- SettleCompensation(frame)
  frame <- settled$frame

  # `estimateLogicle` solves for a linearisation width from the negative values
  # of each channel. When a channel holds a long negative tail the solution
  # needs
  # a wider decade count, and the function stops with "w is negative". The
  # ladder
  # widens `m` until the fit succeeds, and the width that worked is recorded.
  fluorescence <- channels$channel
  transform_list <- NULL
  logicle_m <- NA_real_
  for (width in c(4.5, 5.5, 6.5, 8, 10)) {
    attempt <- try(
      suppressWarnings(
        flowCore::estimateLogicle(frame, channels = fluorescence, m = width)
      ),
      silent = TRUE
    )
    if (!methods::is(attempt, "try-error")) {
      transform_list <- attempt
      logicle_m <- width
      break
    }
  }
  if (is.null(transform_list)) {
    return(Fail("The logicle transform failed at every decade count"))
  }

  windowed <- StableTimeWindow(frame, scatter[["time"]], bins = time_bins)
  transformed <- flowCore::transform(windowed$frame, transform_list)
  events <- flowCore::exprs(transformed)
  total <- nrow(flowCore::exprs(frame))

  Channel <- function(marker) {
    channels$channel[channels$marker == marker]
  }

  singlets <- events[
    SingletMask(events, area = scatter[["forward_area"]],
                height = scatter[["forward_height"]]), ,
    drop = FALSE
  ]

  # Live cells for frozen PBMC, CD45 positive leukocytes for whole blood. The
  # sixth colour carries a different reagent in each material.
  # A marker named in `fixed_cuts` is taken, not fitted.
  Cut <- function(marker, values) {
    if (!is.null(fixed_cuts) && marker %in% names(fixed_cuts) &&
        is.finite(fixed_cuts[[marker]])) {
      return(list(cut = fixed_cuts[[marker]], rule = "fixed"))
    }
    ResolveCut(values)
  }

  entry_marker <- if (material == "PBMC") "NIR" else "CD45"
  entry <- Cut(entry_marker, singlets[, Channel(entry_marker)])
  if (entry$rule == "none") {
    return(Fail(paste0("No cut could be fitted on ", entry_marker)))
  }
  entry_events <- if (material == "PBMC") {
    singlets[singlets[, Channel("NIR")] < entry$cut, , drop = FALSE]
  } else {
    singlets[singlets[, Channel("CD45")] > entry$cut, , drop = FALSE]
  }

  lymphocytes <- LymphocyteGate.(entry_events, scatter)
  if (is.null(lymphocytes)) {
    return(Fail("The lymphocyte gate could not be fitted"))
  }

  cd3 <- Cut("CD3", lymphocytes[, Channel("CD3")])
  if (cd3$rule == "none") {
    return(Fail("No cut could be fitted on CD3"))
  }
  cd3_events <- lymphocytes[lymphocytes[, Channel("CD3")] > cd3$cut, ,
                            drop = FALSE]

  cd4 <- Cut("CD4", cd3_events[, Channel("CD4")])
  cd8 <- Cut("CD8", cd3_events[, Channel("CD8")])
  if (cd4$rule == "none" || cd8$rule == "none") {
    return(Fail("No cut could be fitted on CD4 or on CD8"))
  }
  is_cd4 <- cd3_events[, Channel("CD4")] > cd4$cut
  is_cd8 <- cd3_events[, Channel("CD8")] > cd8$cut
  cd4_events <- cd3_events[is_cd4 & !is_cd8, , drop = FALSE]
  cd8_events <- cd3_events[is_cd8 & !is_cd4, , drop = FALSE]

  # Both quadrant cuts are fitted once, inside CD3, and reused in CD4 and CD8.
  ra <- Cut("CD45RA", cd3_events[, Channel("CD45RA")])
  ccr7 <- Cut("CCR7", cd3_events[, Channel("CCR7")])
  if (ra$rule == "none" || ccr7$rule == "none") {
    return(Fail("No cut could be fitted on CD45RA or on CCR7"))
  }
  ra_cut <- ra$cut
  ccr7_cut <- ccr7$cut

  Quadrants <- function(subset, prefix) {
    if (nrow(subset) == 0) {
      values <- rep(NA_real_, 4)
    } else {
      is_ra <- subset[, Channel("CD45RA")] > ra_cut
      is_ccr7 <- subset[, Channel("CCR7")] > ccr7_cut
      values <- c(
        Fraction.(is_ra & is_ccr7),
        Fraction.(!is_ra & is_ccr7),
        Fraction.(!is_ra & !is_ccr7),
        Fraction.(is_ra & !is_ccr7)
      )
    }
    stats::setNames(values, paste0(prefix, c("_N", "_CM", "_EM", "_TD")))
  }

  frequencies <- c(
    CD3 = Fraction.(lymphocytes[, Channel("CD3")] > cd3$cut),
    CD4 = Fraction.(is_cd4 & !is_cd8),
    CD8 = Fraction.(is_cd8 & !is_cd4),
    Quadrants(cd3_events, "CD3"),
    Quadrants(cd4_events, "CD4"),
    Quadrants(cd8_events, "CD8")
  )

  # The acquisition stamp travels with the counts. Two files of this deposit
  # that
  # carry the same date, the same start time and the same GUID are one
  # acquisition exported twice, and a report that calls them replicates would
  # understate the spread.
  keywords <- flowCore::keyword(frame)
  Keyword <- function(name) {
    value <- keywords[[name]]
    if (is.null(value)) NA_character_ else trimws(as.character(value))
  }

  counts <- data.frame(
    file_name = basename(path),
    material = material,
    compensation = settled$state,
    logicle_m = logicle_m,
    acquisition_date = Keyword("$DATE"),
    acquisition_time = Keyword("$BTIM"),
    acquisition_guid = Keyword("GUID"),
    total_events = total,
    time_kept = windowed$kept,
    time_window_applied = windowed$applied,
    singlet_events = nrow(singlets),
    entry_events = nrow(entry_events),
    lymphocyte_events = nrow(lymphocytes),
    cd3_events = nrow(cd3_events),
    cd4_events = nrow(cd4_events),
    cd8_events = nrow(cd8_events),
    stringsAsFactors = FALSE
  )
  counts <- cbind(counts, as.list(frequencies))

  cuts <- data.frame(
    file_name = basename(path),
    marker = c(entry_marker, "CD3", "CD4", "CD8", "CD45RA", "CCR7"),
    parent = c("singlets", "lymphocytes", "CD3", "CD3", "CD3", "CD3"),
    cut = c(entry$cut, cd3$cut, cd4$cut, cd8$cut, ra$cut, ccr7$cut),
    rule = c(entry$rule, cd3$rule, cd4$rule, cd8$rule, ra$rule, ccr7$rule),
    stringsAsFactors = FALSE
  )

  channels$file_name <- basename(path)
  list(counts = counts, cuts = cuts, channels = channels, error = NA_character_)
}

# Fit the lymphocyte gate on the scatter pair.
#
# A bivariate normal on its own does not work here. In whole blood the parent is
# every CD45 positive leukocyte, so the cloud holds granulocytes and monocytes
# as
# well, and a robust centre fitted to all three keeps nearly all of them. The
# gate therefore runs in two steps. It first cuts side scatter at the density
# minimum above the lowest mode, which is the boundary that separates
# lymphocytes from monocytes and granulocytes. It then fits a bivariate normal
# inside that population to remove the remaining debris.
#
# This is the one step of the hierarchy that only the CD3 percentage of
# lymphocytes depends on. The twelve naive and memory frequencies are counted
# inside CD3, CD4 and CD8, so they do not move with it.
LymphocyteGate. <- function(events, scatter, quantile_limit = 0.95) {
  if (nrow(events) < 100) {
    return(NULL)
  }
  side <- events[, scatter[["side_area"]]]
  side_cut <- DensityCut(side)
  low_side <- if (is.na(side_cut)) {
    events
  } else {
    events[side < side_cut, , drop = FALSE]
  }
  if (nrow(low_side) < 100) {
    return(NULL)
  }

  columns <- c(scatter[["forward_area"]], scatter[["side_area"]])
  pair <- low_side[, columns, drop = FALSE]
  centre <- try(robustbase::covMcd(pair, alpha = 0.75), silent = TRUE)
  if (methods::is(centre, "try-error")) {
    return(low_side)
  }
  distance <- try(
    stats::mahalanobis(pair, centre$center, centre$cov), silent = TRUE
  )
  if (methods::is(distance, "try-error") || anyNA(distance)) {
    return(low_side)
  }
  keep <- distance <= stats::qchisq(quantile_limit, df = 2)
  if (sum(keep) < 100) {
    return(low_side)
  }
  low_side[keep, , drop = FALSE]
}
