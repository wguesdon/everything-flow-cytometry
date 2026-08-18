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

#' Build a flowFrame whose time channel carries a known acquisition rate
#'
#' `StableTimeWindow()` reads the time channel alone, so the frame needs nothing
#' else. With `gap = FALSE` the events arrive at an even rate. With `gap = TRUE`
#' the middle third of the run arrives four times as fast, which is the shape
#' that makes an automated window discard a large part of a file.
#'
#' @param n_events The number of events. Defaults to 4000.
#' @param gap Whether to make the middle third arrive faster.
#' @return A `flowFrame` with the channels `Time` and `FSC-A`.
MakeTimeFlowFrame <- function(n_events = 4000, gap = FALSE) {
  if (gap) {
    slow <- seq(0, 1000, length.out = round(n_events / 2))
    fast <- seq(1000, 1100, length.out = n_events - length(slow))
    time <- sort(c(slow, fast))
  } else {
    time <- seq(0, 1000, length.out = n_events)
  }

  values <- cbind(
    Time = time,
    `FSC-A` = rep(100000, length(time))
  )
  flowCore::flowFrame(values)
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

#' Build a small flowSet of compensation controls for a test
#'
#' The set holds one control per fluorescence channel plus one unstained control.
#' The control names are deliberately spelled in three different ways, so the
#' exact, collapsed and unstained matching passes are each exercised.
#'
#' @return A `flowSet` of three samples.
MakeTestControlSet <- function() {
  MakeNamedFrame <- function(seed) {
    frame <- MakeTestFlowFrame(n_events = 100, seed = seed)
    parameters <- flowCore::parameters(frame)
    # The panel spells the second fluorochrome without a separator, while the
    # control file below spells it with one.
    flowCore::pData(parameters)$desc <- c(NA, NA, "CD3 Ax700", "CD4 PETxRed")
    flowCore::parameters(frame) <- parameters
    frame
  }

  frames <- list(
    "Comp_Beads_CD3 Ax700_A1_A01_001.fcs" = MakeNamedFrame(1),
    "Comp_Beads_CD4 PE-TxRed_B1_B01_002.fcs" = MakeNamedFrame(2),
    "Comp_Beads_unstained_C1_C01_003.fcs" = MakeNamedFrame(3)
  )

  flowCore::flowSet(frames)
}

#' Build a small gated GatingSet for a test
#'
#' `CollectGateTree` reads a real hierarchy, so a test needs one. The gates are
#' plain rectangles added by hand rather than an openCyto template, because the
#' function under test reads the hierarchy and not the way it was built.
#'
#' The two gates are nested, and each one drops a known share of its parent, so
#' a test can assert the arithmetic as well as the shape.
#'
#' @param n_events The number of events. Defaults to 400.
#' @return A `GatingSet` with `nonDebris` under root and `singlets` under
#'   `nonDebris`.
MakeGatedSet <- function(n_events = 400) {
  frame <- MakeTestFlowFrame(n_events = n_events)
  flow_set <- flowCore::flowSet(sample = frame)
  gating_set <- flowWorkspace::GatingSet(flow_set)

  wide <- flowCore::rectangleGate(
    filterId = "nonDebris",
    list(`FSC-A` = c(80000, Inf), `FSC-H` = c(-Inf, Inf)))
  flowWorkspace::gs_pop_add(gating_set, wide, parent = "root")

  narrow <- flowCore::rectangleGate(
    filterId = "singlets",
    list(`FSC-A` = c(-Inf, Inf), `FSC-H` = c(85000, Inf)))
  flowWorkspace::gs_pop_add(gating_set, narrow, parent = "nonDebris")

  flowWorkspace::recompute(gating_set)
  gating_set
}

#' Write a small flowFrame to a temporary FCS file
#'
#' `DescribeFcsFile` and its siblings read a path rather than a frame, so a test
#' needs a file. The file is written into a temporary folder that the caller
#' owns, so nothing is left behind.
#'
#' @param directory The folder to write into.
#' @param name The file name. Defaults to `"sample.fcs"`.
#' @param keywords Extra FCS keywords to store, as a named list.
#' @param n_events The number of events.
#' @return The path that was written.
WriteTestFcs <- function(directory, name = "sample.fcs", keywords = list(),
                         n_events = 200) {
  frame <- MakeTestFlowFrame(n_events = n_events)
  for (key in names(keywords)) {
    flowCore::keyword(frame)[[key]] <- keywords[[key]]
  }
  path <- file.path(directory, name)
  suppressWarnings(flowCore::write.FCS(frame, path))
  path
}
