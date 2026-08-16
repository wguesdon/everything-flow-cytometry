# Spectral flow cytometry.
#
# A spectral analyser records the emission of every fluorochrome across all of its
# detectors and then solves for the contribution of each one. The file that leaves
# the instrument is therefore already unmixed, and it carries one column per
# fluorochrome rather than one column per detector.
#
# Two consequences shape this file.
#
# There is no spillover matrix to apply, because unmixing replaced compensation
# before export. `ApplyCompensation()` has nothing to do here and the pipeline
# starts at the transform.
#
# The viability channel of the Yu 2021 panel carries no marker name, so
# `FluorescenceChannels()` drops it. That function selects a channel only when
# `$PnS` is set, which is correct for a conventional panel and wrong here.
# `SpectralChannels()` keeps every channel that is not scatter and not time.

#' Select every fluorescence channel of a spectral flowFrame
#'
#' Unlike [FluorescenceChannels()], a channel with no marker name is kept. A
#' spectral panel routinely leaves the viability dye unnamed, and that channel
#' still needs a transform because the first gate uses it.
#'
#' @param frame A `flowFrame`.
#' @param exclude Extra channel names to drop. Matching ignores case.
#' @return A character vector of detector names.
#' @examples
#' \dontrun{
#' SpectralChannels(frame)
#' }
#' @export
SpectralChannels <- function(frame, exclude = character()) {
  channels <- DescribeChannels(frame)

  is_scatter <- grepl("^(FSC|SSC)", channels$channel, ignore.case = TRUE)
  is_time <- grepl("^time$", channels$channel, ignore.case = TRUE)
  is_excluded <- tolower(channels$channel) %in% tolower(exclude)

  channels$channel[!is_scatter & !is_time & !is_excluded]
}

#' Read the sample sheet that describes the Yu 2021 deposit
#'
#' The sheet is derived from Table S1 of the paper and it is committed, so the
#' clinical grouping travels with the repository rather than with a download.
#'
#' @param path Path to the CSV.
#' @return A `data.frame` with one row per FCS file. `severity_rank` is returned
#'   as an ordered factor running from `normal` to `hospitalized`, and
#'   `severity_index` gives the same order as an integer from 1 to 4.
#' @export
ReadYuSampleSheet <- function(path) {
  if (!file.exists(path)) {
    stop("The sample sheet does not exist: ", path)
  }

  sheet <- utils::read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)

  required <- c("file_name", "alias_subject_id", "sex", "severity_rank",
                "timepoint", "igg_result")
  missing <- setdiff(required, colnames(sheet))
  if (length(missing) > 0) {
    stop("The sample sheet is missing the column(s): ",
         paste(missing, collapse = ", "), ".")
  }

  levels_in_order <- c("normal", "exposed", "infected", "hospitalized")
  unknown <- setdiff(unique(sheet$severity_rank), levels_in_order)
  if (length(unknown) > 0) {
    stop("The sample sheet holds an unknown severity_rank: ",
         paste(unknown, collapse = ", "), ".")
  }

  sheet$severity_rank <- factor(sheet$severity_rank, levels = levels_in_order,
                                ordered = TRUE)
  sheet$severity_index <- as.integer(sheet$severity_rank)
  sheet
}

#' Estimate one logicle transform from a reference file
#'
#' Every file in a cohort must share one transform, otherwise a gate fitted on one
#' sample means something different on the next one. This reads a single file and
#' returns the transform, so the caller can apply it to all of them without
#' holding the cohort in memory.
#'
#' @param path Path to the FCS file used as the reference.
#' @param channels The detector names to transform. Defaults to
#'   [SpectralChannels()] of the reference frame.
#' @return A list with `transform`, a `transformList`, and `channels`.
#' @export
EstimateSpectralTransform <- function(path, channels = NULL) {
  if (!file.exists(path)) {
    stop("The reference file does not exist: ", path)
  }

  frame <- flowCore::read.FCS(path, truncate_max_range = FALSE)
  if (is.null(channels)) {
    channels <- SpectralChannels(frame)
  }

  unknown <- setdiff(channels, flowCore::colnames(frame))
  if (length(unknown) > 0) {
    stop("These channels are not in the reference file: ",
         paste(unknown, collapse = ", "), ".")
  }

  list(
    transform = flowCore::estimateLogicle(frame, channels = channels),
    channels = channels
  )
}

