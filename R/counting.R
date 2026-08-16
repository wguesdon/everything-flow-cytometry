# Rare event counting.
#
# OMIP-043 is a rare cell panel, and it states its design target in counting
# terms: a 5 percent coefficient of variation, reached by collecting 400 to 2,000
# events in the antibody secreting cell gate. Both halves of that sentence are
# arithmetic, so both can be checked against the deposited data.
#
# Two coefficients of variation appear below and they answer different questions.
# The Poisson CV is the floor set by counting statistics alone, and no amount of
# care in the laboratory goes below it. The measured CV is what the replicates
# actually did. The gap between them is everything that is not counting noise.

#' Coefficient of variation from counting statistics alone
#'
#' A count of rare events follows a Poisson distribution, whose standard deviation
#' is the square root of the count. The relative spread is therefore
#' `sqrt(n) / n`, which is `1 / sqrt(n)`. Collecting 400 events gives 5 percent
#' and collecting 2,500 gives 2 percent.
#'
#' This is a floor, not an estimate. A measured CV below it means the replicates
#' are not independent.
#'
#' @param count The number of events in the gate. May be a vector.
#' @return The Poisson coefficient of variation as a percentage. A count of zero
#'   returns `NA`, because a relative spread is undefined there.
#' @examples
#' PoissonCv(c(400, 2000))
#' @export
PoissonCv <- function(count) {
  if (any(count < 0, na.rm = TRUE)) {
    stop("A count cannot be negative.")
  }
  ifelse(is.na(count) | count == 0, NA_real_, 100 / sqrt(count))
}

#' The number of events needed to reach a target coefficient of variation
#'
#' The inverse of [PoissonCv()]. Rearranging `cv = 100 / sqrt(n)` gives
#' `n = (100 / cv)^2`.
#'
#' @param target_cv_percent The wanted coefficient of variation, as a percentage.
#' @return The number of events, rounded up.
#' @examples
#' EventsForCv(5)
#' @export
EventsForCv <- function(target_cv_percent) {
  if (any(target_cv_percent <= 0)) {
    stop("target_cv_percent must be above zero.")
  }
  ceiling((100 / target_cv_percent)^2)
}

#' Coefficient of variation of a set of measurements
#'
#' @param values A numeric vector, for example one frequency per replicate.
#' @return The coefficient of variation as a percentage, or `NA` when the mean is
#'   zero or fewer than two values are given.
#' @export
MeasuredCv <- function(values) {
  values <- values[!is.na(values)]
  if (length(values) < 2) {
    return(NA_real_)
  }
  mean_value <- mean(values)
  if (mean_value == 0) {
    return(NA_real_)
  }
  100 * stats::sd(values) / mean_value
}

#' Compare the measured spread against the counting floor
#'
#' The measured CV holds counting noise plus everything else: pipetting, staining,
#' acquisition and the gate. Subtracting in quadrature separates the two, because
#' independent sources of variance add.
#'
#' @param counts The event count per replicate.
#' @param frequencies The frequency per replicate, as a percentage.
#' @return A one row `data.frame` with `replicates`, `mean_count`,
#'   `mean_percent`, `poisson_cv_percent`, `measured_cv_percent` and
#'   `excess_cv_percent`. The excess is `NA` when the measured CV sits below the
#'   Poisson floor, which is a result worth seeing rather than a number to force.
#' @export
CompareSpreadToPoisson <- function(counts, frequencies) {
  if (length(counts) != length(frequencies)) {
    stop("counts and frequencies must be the same length.")
  }

  mean_count <- mean(counts, na.rm = TRUE)
  poisson_cv <- PoissonCv(mean_count)
  measured_cv <- MeasuredCv(frequencies)

  excess <- if (is.na(measured_cv) || is.na(poisson_cv) ||
                measured_cv <= poisson_cv) {
    NA_real_
  } else {
    sqrt(measured_cv^2 - poisson_cv^2)
  }

  data.frame(
    replicates = length(counts),
    mean_count = mean_count,
    mean_percent = mean(frequencies, na.rm = TRUE),
    poisson_cv_percent = poisson_cv,
    measured_cv_percent = measured_cv,
    excess_cv_percent = excess,
    stringsAsFactors = FALSE
  )
}
