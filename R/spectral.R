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


#' Read the gate cut points that the cohort shares
#'
#' Every sample is gated at the same cut. The alternative, a cut fitted per
#' sample, moves the boundary for reasons that have nothing to do with the
#' biology, and every claim in this analysis is a comparison between samples. One
#' cut for all 83 files means a difference between two samples is a difference in
#' the data.
#'
#' The price is that a sample whose staining drifted is gated at the wrong place
#' and nothing corrects for it. [CutPointDiagnostics()] reports that.
#'
#' @param path Path to a CSV with the columns `marker`, `parent`, `side`, `cut`,
#'   `source` and `note`. `side` is `"above"` or `"below"`.
#' @return A `data.frame` of cut points.
#' @export
ReadGateCuts <- function(path) {
  if (!file.exists(path)) {
    stop("The gate cut file does not exist: ", path)
  }

  cuts <- utils::read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
  required <- c("marker", "parent", "side", "cut", "source", "note")
  missing <- setdiff(required, colnames(cuts))
  if (length(missing) > 0) {
    stop("The gate cut file is missing the column(s): ",
         paste(missing, collapse = ", "), ".")
  }

  bad_side <- setdiff(cuts$side, c("above", "below"))
  if (length(bad_side) > 0) {
    stop("A gate side must be 'above' or 'below'. Found: ",
         paste(unique(bad_side), collapse = ", "), ".")
  }

  cuts
}

#' Keep the singlets of a transformed event matrix
#'
#' A doublet gives about twice the area for the same height, so singlets lie on a
#' line through the origin in FSC-A against FSC-H. The line is fitted per sample
#' and events further than `tolerance` robust deviations from it are dropped. The
#' fit is per sample because the scatter gain, not the biology, sets its slope.
#'
#' @param events A numeric matrix with the columns `FSC-A` and `FSC-H`.
#' @param tolerance The number of median absolute deviations to keep.
#' @return A logical vector, `TRUE` for an event to keep.
#' @export
SingletMask <- function(events, tolerance = 3) {
  needed <- c("FSC-A", "FSC-H")
  missing <- setdiff(needed, colnames(events))
  if (length(missing) > 0) {
    stop("The event matrix has no ", paste(missing, collapse = " and "),
         " column, so singlets cannot be gated.")
  }

  area <- events[, "FSC-A"]
  height <- events[, "FSC-H"]

  fit <- stats::lm(height ~ area)
  residual <- stats::residuals(fit)
  spread <- stats::mad(residual)
  if (!is.finite(spread) || spread == 0) {
    return(rep(TRUE, nrow(events)))
  }

  abs(residual) <= tolerance * spread
}

#' Gate one file against the shared cut points
#'
#' The hierarchy is singlets, live, CD45 positive, then CD3 positive and CD3
#' negative, then CD8 positive under CD3 positive and NK under CD3 negative. The
#' CD8 gate carries CD161hi and the CD45RA against CCR7 quadrant.
#'
#' @param path Path to the FCS file.
#' @param transform_list A `transformList`, from [EstimateSpectralTransform()].
#' @param cuts The output of [ReadGateCuts()].
#' @param panel The output of [DescribeChannels()] on the reference file, used to
#'   resolve a marker name to a detector name.
#' @param cd45ra_cut Overrides the CD45RA cut in `cuts`. A vector returns one row
#'   of counts per value, which is how the sweep in
#'   `scripts/07_yu2021_spectral_mait.R` shows what the choice changes without
#'   reading every file again.
#' @param subsample_n Live CD45 positive events to draw for the pooled
#'   clustering. Pass 0 to skip the draw.
#' @param seed The random seed for the draw.
#' @return A list with `counts`, a one row `data.frame`, `events`, the drawn
#'   matrix or `NULL`, `medians`, the CD56 median in the NK and CD8 gates, and
#'   `error`, a message when the file could not be read and `NA` when it could.
#' @export
GateSpectralFile <- function(path,
                             transform_list,
                             cuts,
                             panel,
                             cd45ra_cut = NULL,
                             subsample_n = 3000,
                             seed = 42) {
  Detector <- function(marker) {
    if (marker %in% panel$channel) {
      return(marker)
    }
    hit <- panel$channel[tolower(trimws(panel$marker)) == tolower(trimws(marker))]
    if (length(hit) == 0) {
      stop("No channel carries the marker '", marker, "'.")
    }
    hit[1]
  }

  Cut <- function(marker) {
    value <- cuts$cut[cuts$marker == marker]
    if (length(value) == 0) {
      stop("No cut point is recorded for '", marker, "'.")
    }
    value[1]
  }

  read_result <- try(
    flowCore::read.FCS(path, truncate_max_range = FALSE), silent = TRUE
  )
  if (methods::is(read_result, "try-error")) {
    return(list(counts = NULL, events = NULL, medians = NULL,
                error = trimws(as.character(read_result))))
  }

  events <- flowCore::exprs(flowCore::transform(read_result, transform_list))

  # Three files carry an extra SampleID parameter from a concatenation step.
  shared <- intersect(panel$channel, colnames(events))
  events <- events[, shared, drop = FALSE]

  ra_cuts <- if (is.null(cd45ra_cut)) Cut("CD45ra") else cd45ra_cut

  singlets <- events[SingletMask(events), , drop = FALSE]
  live <- singlets[singlets[, "LIVE DEAD Blue-A"] < Cut("LIVE DEAD Blue-A"), ,
                   drop = FALSE]
  cd45pos <- live[live[, Detector("CD45")] > Cut("CD45"), , drop = FALSE]

  is_cd3 <- cd45pos[, Detector("CD3")] > Cut("CD3")
  cd3pos <- cd45pos[is_cd3, , drop = FALSE]
  cd3neg <- cd45pos[!is_cd3, , drop = FALSE]

  nk <- cd3neg[cd3neg[, Detector("CD56")] > Cut("CD56"), , drop = FALSE]
  cd8pos <- cd3pos[cd3pos[, Detector("CD8")] > Cut("CD8"), , drop = FALSE]

  ra_value <- cd8pos[, Detector("CD45ra")]
  is_ccr7 <- cd8pos[, Detector("CCR7")] > Cut("CCR7")
  is_cd161 <- cd8pos[, Detector("CD161")] > Cut("CD161")

  # One row per CD45RA cut. The file is read once and every cut in the sweep is
  # answered from the same events, so the sweep costs one pass and not one pass
  # per cut.
  counts <- do.call(rbind, lapply(ra_cuts, function(ra_cut) {
    is_ra <- ra_value > ra_cut
    data.frame(
      sample = basename(path),
      total_events = nrow(events),
      singlet_events = nrow(singlets),
      live_events = nrow(live),
      cd45pos_events = nrow(cd45pos),
      cd3pos_events = nrow(cd3pos),
      cd3neg_events = nrow(cd3neg),
      nk_events = nrow(nk),
      cd8pos_events = nrow(cd8pos),
      cd161hi_events = sum(is_cd161),
      naive_events = sum(is_ra & is_ccr7),
      cm_events = sum(!is_ra & is_ccr7),
      em_events = sum(!is_ra & !is_ccr7),
      emra_events = sum(is_ra & !is_ccr7),
      cd45ra_cut = ra_cut,
      stringsAsFactors = FALSE
    )
  }))

  drawn <- NULL
  if (subsample_n > 0 && nrow(cd45pos) > 0) {
    withr::with_seed(seed, {
      keep <- if (nrow(cd45pos) <= subsample_n) {
        seq_len(nrow(cd45pos))
      } else {
        sample.int(nrow(cd45pos), subsample_n)
      }
    })
    drawn <- cd45pos[keep, , drop = FALSE]
  }

  medians <- data.frame(
    sample = basename(path),
    nk_cd56_median = if (nrow(nk) > 0) {
      stats::median(nk[, Detector("CD56")])
    } else {
      NA_real_
    },
    cd8_cd56_median = if (nrow(cd8pos) > 0) {
      stats::median(cd8pos[, Detector("CD56")])
    } else {
      NA_real_
    },
    stringsAsFactors = FALSE
  )

  list(counts = counts, events = drawn, medians = medians,
       error = NA_character_)
}

#' Turn the event counts into the frequencies the paper's claims are about
#'
#' @param counts The stacked output of [GateSpectralFile()].
#' @return The input with the percentage columns added. Memory is the sum of the
#'   central memory, effector memory and EMRA gates, which is the definition the
#'   paper states.
#' @export
AddCd8Frequencies <- function(counts) {
  required <- c("cd8pos_events", "cd3pos_events", "cd45pos_events",
                "cd161hi_events", "naive_events", "cm_events", "em_events",
                "emra_events")
  missing <- setdiff(required, colnames(counts))
  if (length(missing) > 0) {
    stop("counts is missing the column(s): ", paste(missing, collapse = ", "),
         ".")
  }

  Percent <- function(numerator, denominator) {
    ifelse(denominator > 0, 100 * numerator / denominator, NA_real_)
  }

  counts$memory_events <- counts$cm_events + counts$em_events +
    counts$emra_events
  counts$cd3_percent_of_cd45 <- Percent(counts$cd3pos_events,
                                        counts$cd45pos_events)
  counts$cd8_percent_of_cd3 <- Percent(counts$cd8pos_events,
                                       counts$cd3pos_events)
  counts$cd161hi_percent_of_cd8 <- Percent(counts$cd161hi_events,
                                           counts$cd8pos_events)
  counts$memory_percent_of_cd8 <- Percent(counts$memory_events,
                                          counts$cd8pos_events)
  counts$naive_percent_of_cd8 <- Percent(counts$naive_events,
                                         counts$cd8pos_events)
  counts$nk_percent_of_cd45 <- Percent(counts$nk_events, counts$cd45pos_events)
  counts
}
