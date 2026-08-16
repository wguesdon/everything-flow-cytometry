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

#' Gate one spectral file and return its population statistics
#'
#' The cohort is gated one file at a time rather than as a single `flowSet`. The
#' deposit holds 23 million events, and reading all of them at once costs about
#' 8 GB before any gate is fitted. Gating file by file keeps the peak memory at
#' one file and lets a failure on one sample be recorded instead of stopping the
#' run.
#'
#' @param path Path to the FCS file.
#' @param template A `gatingTemplate`, from [ReadGatingTemplate()].
#' @param transform_list A `transformList`, from [EstimateSpectralTransform()].
#' @param subsample_population The population to draw pooled events from, for the
#'   clustering step. Pass `NULL` to skip the draw.
#' @param subsample_n The number of events to draw. The paper drew 3000 live CD45+
#'   cells per sample for its clustering figure.
#' @param median_populations Populations to summarise by marker median. Pass
#'   `NULL` to skip.
#' @param median_channels The channels to take the median of, in the transformed
#'   scale.
#' @param seed The random seed for the draw.
#' @return A list with `stats`, a `data.frame` of population counts, `events`, the
#'   drawn matrix or `NULL`, `medians`, a `data.frame` of marker medians or
#'   `NULL`, and `error`, a message when the gating failed and `NA` when it did
#'   not.
#' @export
GateOneSpectralFile <- function(path,
                                template,
                                transform_list,
                                subsample_population = "CD45pos",
                                subsample_n = 3000,
                                median_populations = NULL,
                                median_channels = NULL,
                                seed = 42) {
  if (!file.exists(path)) {
    stop("The FCS file does not exist: ", path)
  }

  empty <- list(stats = NULL, events = NULL, medians = NULL,
                error = NA_character_)

  result <- try({
    frame <- flowCore::read.FCS(path, truncate_max_range = FALSE)
    transformed <- flowCore::transform(frame, transform_list)

    flow_set <- methods::as(list(transformed), "flowSet")
    flowCore::sampleNames(flow_set) <- basename(path)

    gating_set <- flowWorkspace::GatingSet(flow_set)
    openCyto::gt_gating(template, gating_set, mc.cores = 1,
                        parallel_type = "none")
    flowWorkspace::recompute(gating_set)
    gating_set
  }, silent = TRUE)

  if (methods::is(result, "try-error")) {
    empty$error <- trimws(as.character(result))
    return(empty)
  }

  stats <- CollectPopulationStats(result, sample_sheet = NULL)

  events <- NULL
  if (!is.null(subsample_population)) {
    drawn <- try(
      flowWorkspace::gh_pop_get_data(result[[1]], subsample_population),
      silent = TRUE
    )
    if (!methods::is(drawn, "try-error")) {
      matrix_data <- flowCore::exprs(drawn)
      withr::with_seed(seed, {
        keep <- if (nrow(matrix_data) <= subsample_n) {
          seq_len(nrow(matrix_data))
        } else {
          sample.int(nrow(matrix_data), subsample_n)
        }
      })
      events <- matrix_data[keep, , drop = FALSE]
    }
  }

  medians <- NULL
  if (!is.null(median_populations) && !is.null(median_channels)) {
    rows <- lapply(median_populations, function(population) {
      drawn <- try(
        flowWorkspace::gh_pop_get_data(result[[1]], population), silent = TRUE
      )
      if (methods::is(drawn, "try-error")) {
        return(NULL)
      }
      matrix_data <- flowCore::exprs(drawn)
      present <- intersect(median_channels, colnames(matrix_data))
      if (length(present) == 0 || nrow(matrix_data) == 0) {
        return(NULL)
      }
      data.frame(
        sample = basename(path),
        population = population,
        channel = present,
        events = nrow(matrix_data),
        median_intensity = vapply(
          present, function(ch) stats::median(matrix_data[, ch]), numeric(1)
        ),
        stringsAsFactors = FALSE
      )
    })
    medians <- do.call(rbind, rows)
    if (!is.null(medians)) {
      rownames(medians) <- NULL
    }
  }

  list(stats = stats, events = events, medians = medians,
       error = NA_character_)
}

#' Turn the long population table into one row per sample
#'
#' Every claim in the paper is a frequency within the CD8 compartment, so the
#' denominators are fixed here once rather than in each test.
#'
#' @param stats The stacked output of [GateOneSpectralFile()], carrying the
#'   columns `sample`, `population` and `count`.
#' @return A `data.frame` with one row per sample and the columns
#'   `cd161hi_percent_of_cd8`, `memory_percent_of_cd8`, `naive_percent_of_cd8`,
#'   `cd8_percent_of_cd3`, `cd3_percent_of_lymphocytes` and the event counts the
#'   percentages came from.
#' @export
SummariseCd8Compartment <- function(stats) {
  required <- c("sample", "population", "count")
  missing <- setdiff(required, colnames(stats))
  if (length(missing) > 0) {
    stop("stats is missing the column(s): ", paste(missing, collapse = ", "), ".")
  }

  stats$leaf <- basename(as.character(stats$population))

  CountFor <- function(piece, name) {
    hit <- piece$count[piece$leaf == name]
    if (length(hit) == 0) NA_real_ else as.numeric(hit[1])
  }

  pieces <- split(stats, stats$sample)
  rows <- lapply(names(pieces), function(sample_name) {
    piece <- pieces[[sample_name]]

    lymphocytes <- CountFor(piece, "Lymphocytes")
    cd3 <- CountFor(piece, "CD3pos")
    cd8 <- CountFor(piece, "CD8pos")
    cd161hi <- CountFor(piece, "CD161hi")
    naive <- CountFor(piece, "Naive")
    memory <- sum(
      CountFor(piece, "CM"), CountFor(piece, "EM"), CountFor(piece, "EMRA"),
      na.rm = TRUE
    )

    Percent <- function(numerator, denominator) {
      if (is.na(numerator) || is.na(denominator) || denominator == 0) {
        return(NA_real_)
      }
      100 * numerator / denominator
    }

    data.frame(
      sample = sample_name,
      lymphocyte_events = lymphocytes,
      cd3_events = cd3,
      cd8_events = cd8,
      cd161hi_events = cd161hi,
      memory_events = memory,
      naive_events = naive,
      cd3_percent_of_lymphocytes = Percent(cd3, lymphocytes),
      cd8_percent_of_cd3 = Percent(cd8, cd3),
      cd161hi_percent_of_cd8 = Percent(cd161hi, cd8),
      memory_percent_of_cd8 = Percent(memory, cd8),
      naive_percent_of_cd8 = Percent(naive, cd8),
      stringsAsFactors = FALSE
    )
  })

  result <- do.call(rbind, rows)
  rownames(result) <- NULL
  result
}
