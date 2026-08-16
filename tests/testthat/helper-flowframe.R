# Test fixtures.
#
# A test builds its own small flowFrame rather than reading an FCS file. A test
# that reads a 17 MB file is slow, and it fails on a machine that has not pulled
# the archive from S3.

#' Build a small flowFrame for a test
#'
#' The frame carries two scatter channels with no marker and two fluorescence
#' channels with a marker, which is the smallest shape that exercises the channel
#' and marker code.
#'
#' @param n_events The number of events. Defaults to 200.
#' @param seed The seed for the random draw, so a test is deterministic.
#' @return A `flowFrame` with the channels `FSC-A`, `FSC-H`, `Ax700-A` and
#'   `PE-TxRed-A`, and the markers `CD3` and `CD4` on the last two.
MakeTestFlowFrame <- function(n_events = 200, seed = 42) {
  set.seed(seed)

  values <- cbind(
    `FSC-A` = stats::rnorm(n_events, mean = 100000, sd = 15000),
    `FSC-H` = stats::rnorm(n_events, mean = 95000, sd = 14000),
    `Ax700-A` = stats::rnorm(n_events, mean = 2000, sd = 500),
    `PE-TxRed-A` = stats::rnorm(n_events, mean = 1500, sd = 400)
  )

  frame <- flowCore::flowFrame(values)

  parameters <- flowCore::parameters(frame)
  flowCore::pData(parameters)$desc <- c(NA, NA, "CD3", "CD4")
  flowCore::parameters(frame) <- parameters

  frame
}

#' Build a small square spillover matrix for a test
#'
#' @param channels The detector names for the rows and the columns.
#' @param spill The single off diagonal value, as a fraction.
#' @return A numeric matrix with 1 on the diagonal and `spill` off it.
MakeTestSpillover <- function(channels = c("Ax700-A", "PE-TxRed-A"),
                              spill = 0.1) {
  n <- length(channels)
  matrix_out <- matrix(spill, nrow = n, ncol = n)
  diag(matrix_out) <- 1
  dimnames(matrix_out) <- list(channels, channels)
  matrix_out
}

#' Build a flowFrame that carries a spillover matrix in its keywords
#'
#' @param spill The single off diagonal value, as a fraction.
#' @return A `flowFrame` whose `SPILL` keyword holds a two by two matrix.
MakeTestFlowFrameWithSpillover <- function(spill = 0.1) {
  frame <- MakeTestFlowFrame()
  spillover <- MakeTestSpillover(spill = spill)
  flowCore::keyword(frame)[["SPILL"]] <- spillover
  frame
}
