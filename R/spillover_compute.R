# Compute a spillover matrix from single stained controls.
#
# The example in scripts/01 applies the matrix that the instrument wrote into the
# file. That covers the common case, where compensation was set on the cytometer
# in Diva. It does not cover the case where the files were acquired uncompensated,
# or where the stored matrix is wrong and has to be recomputed.
#
# FlowJo computes the matrix in its interface and shows an N by N grid so you can
# check every pair by eye. R can do both. flowStats::spillover_ng() computes the
# matrix, and it needs a match file that maps each control file to the channel it
# stains. FlowJo builds that mapping through its interface. The functions here
# build it from the FCS metadata instead, so the mapping is derived and not typed.

#' Strip the prefix and the well suffix from a control sample name
#'
#' A compensation control is usually named for the stain it carries, wrapped in a
#' prefix that names the control type and a suffix that names the plate well. The
#' stain in the middle is the part that has to match a marker in the panel.
#'
#' @param sample_names A character vector of sample names or file names.
#' @param prefixes Prefixes to remove, for example `"Comp_Beads_"`. Matching
#'   ignores case.
#' @return A character vector of stain names.
#' @examples
#' StripControlName("Comp_Beads_TCR Vd1 FITC_A1_A01_001.fcs")
#' @export
StripControlName <- function(sample_names,
                             prefixes = c("Comp_Beads_", "Comp_Cells_",
                                          "Compensation Controls_",
                                          "Single stainings_", "Single Stains_",
                                          "Controls_", "Comp_")) {
  stains <- basename(sample_names)
  stains <- sub("\\.fcs$", "", stains, ignore.case = TRUE)

  for (prefix in prefixes) {
    stains <- sub(paste0("^", prefix), "", stains, ignore.case = TRUE)
  }

  # Remove a trailing plate position and acquisition number, for example
  # "_A1_A01_001" or "_005".
  stains <- sub("_[A-H][0-9]{1,2}_[A-H][0-9]{2}_[0-9]+$", "", stains)
  stains <- sub("_[0-9]{3,}$", "", stains)

  trimws(stains)
}

#' Normalise a marker or stain name so two spellings of it compare equal
#'
#' Control file names and panel markers rarely agree character for character. The
#' same stain appears as `LD UV Blue` in one place and `Live Dead UV Blue` in
#' another, or as `TCR Va7_2` against `TCR Va7.2`. Normalising both sides catches
#' the difference that is only punctuation or case.
#'
#' @param x A character vector.
#' @return A character vector in lower case with every run of non alphanumeric
#'   characters replaced by a single space.
#' @export
NormaliseMarkerName <- function(x) {
  out <- tolower(x)
  out <- gsub("[^a-z0-9]+", " ", out)
  trimws(gsub("\\s+", " ", out))
}

#' Collapse a marker name to letters and digits only
#'
#' A stricter form of [NormaliseMarkerName()]. It removes every separator instead
#' of replacing it with a space, so `PerCP-Cy55` and `PerCPCy55` compare equal.
#' That pair is real: OMIP-39 writes the first form on the control file and the
#' second in the panel.
#'
#' @param x A character vector.
#' @return A character vector of lower case letters and digits.
#' @export
CollapseMarkerName <- function(x) {
  gsub("[^a-z0-9]+", "", tolower(x))
}

#' Take the antibody name off the front of a stain label
#'
#' A control is labelled with the antibody and then the fluorochrome, as in
#' `CD2 PerCP-Cy55`. When the fluorochrome is spelled differently on the two sides,
#' the antibody alone still identifies the channel. This returns that first token.
#'
#' @param x A character vector of stain or marker labels.
#' @return A character vector holding the first whitespace delimited token, in
#'   lower case.
#' @examples
#' AntibodyToken(c("CD2 PerCP-Cy55", "NKp30 eFluor450"))
#' @export
AntibodyToken <- function(x) {
  tolower(trimws(sub("\\s.*$", "", trimws(x))))
}

#' Match every compensation control to the channel it stains
#'
#' This is the step that FlowJo performs in its interface. Each control file
#' carries one stain, and that stain names one detector in the panel. The match is
#' made on the marker names that the FCS files already carry, so no mapping is
#' typed by hand.
#'
#' Matching runs in two passes. The first pass compares the stain to the marker
#' character for character. The second pass compares the normalised forms, which
#' catches a difference of case or punctuation. A control that neither pass
#' resolves is returned with `channel` set to `NA`, so it is visible rather than
#' silently dropped.
#'
#' @param flow_set A `flowSet` of compensation controls.
#' @param unstained_pattern A regular expression that identifies the unstained
#'   control. Matching ignores case.
#' @param reference The sample whose panel supplies the marker names. Defaults to
#'   the first sample.
#' @return A `data.frame` with the columns `filename`, `stain`, `channel`,
#'   `marker` and `matched_by`. `matched_by` is one of `"exact"`, `"normalised"`,
#'   `"unstained"` or `"none"`.
#' @export
MatchControlsToChannels <- function(flow_set,
                                    unstained_pattern = "unstain",
                                    reference = 1) {
  if (!methods::is(flow_set, "flowSet")) {
    stop("flow_set must be a flowSet, not a ", class(flow_set)[1], ".")
  }

  panel <- DescribeChannels(flow_set[[reference]])
  panel <- panel[panel$is_marker, ]

  filenames <- flowCore::sampleNames(flow_set)
  stains <- StripControlName(filenames)

  channel <- rep(NA_character_, length(stains))
  marker <- rep(NA_character_, length(stains))
  matched_by <- rep("none", length(stains))

  is_unstained <- grepl(unstained_pattern, stains, ignore.case = TRUE)
  channel[is_unstained] <- "unstained"
  matched_by[is_unstained] <- "unstained"

  exact_hit <- match(stains, panel$marker)
  use_exact <- !is_unstained & !is.na(exact_hit)
  channel[use_exact] <- panel$channel[exact_hit[use_exact]]
  marker[use_exact] <- panel$marker[exact_hit[use_exact]]
  matched_by[use_exact] <- "exact"

  # Pass 2, normalised: case and punctuation are ignored.
  normal_hit <- match(
    NormaliseMarkerName(stains),
    NormaliseMarkerName(panel$marker)
  )
  use_normal <- is.na(channel) & !is_unstained & !is.na(normal_hit)
  channel[use_normal] <- panel$channel[normal_hit[use_normal]]
  marker[use_normal] <- panel$marker[normal_hit[use_normal]]
  matched_by[use_normal] <- "normalised"

  # Pass 3, collapsed: every separator is removed, so PerCP-Cy55 equals PerCPCy55.
  collapsed_hit <- match(
    CollapseMarkerName(stains),
    CollapseMarkerName(panel$marker)
  )
  use_collapsed <- is.na(channel) & !is_unstained & !is.na(collapsed_hit)
  channel[use_collapsed] <- panel$channel[collapsed_hit[use_collapsed]]
  marker[use_collapsed] <- panel$marker[collapsed_hit[use_collapsed]]
  matched_by[use_collapsed] <- "collapsed"

  # Pass 4, antibody only: the fluorochrome is spelled differently on the two
  # sides, so match on the antibody name alone. This pass is used only when the
  # antibody name identifies exactly one channel, because a token that appears in
  # two markers cannot resolve a channel.
  panel_token <- AntibodyToken(panel$marker)
  unique_token <- panel_token[!duplicated(panel_token) &
                                !(panel_token %in% panel_token[duplicated(panel_token)])]
  token_hit <- match(AntibodyToken(stains), panel_token)
  token_is_unique <- AntibodyToken(stains) %in% unique_token
  use_token <- is.na(channel) & !is_unstained & !is.na(token_hit) & token_is_unique
  channel[use_token] <- panel$channel[token_hit[use_token]]
  marker[use_token] <- panel$marker[token_hit[use_token]]
  matched_by[use_token] <- "antibody"

  data.frame(
    filename = filenames,
    stain = stains,
    channel = channel,
    marker = marker,
    matched_by = matched_by,
    stringsAsFactors = FALSE
  )
}

#' Write the match table in the format that flowStats expects
#'
#' [flowStats::spillover_ng()] reads a CSV with the columns `filename` and
#' `channel`, and it needs exactly one row whose channel is `unstained`. This
#' function drops the unmatched controls, keeps one unstained control and writes
#' the file.
#'
#' @param match_table The output of [MatchControlsToChannels()].
#' @param path Where to write the CSV.
#' @param keep_unstained Which unstained control to keep when the set holds more
#'   than one. Give the file name, or leave it `NULL` to keep the first.
#' @return The path, invisibly. A warning names every control that was dropped.
#' @export
WriteMatchFile <- function(match_table, path, keep_unstained = NULL) {
  required <- c("filename", "channel", "matched_by")
  missing <- setdiff(required, colnames(match_table))
  if (length(missing) > 0) {
    stop("match_table is missing the column(s): ",
         paste(missing, collapse = ", "), ".")
  }

  unmatched <- match_table$filename[match_table$matched_by == "none"]
  if (length(unmatched) > 0) {
    warning(
      length(unmatched), " control(s) matched no channel and were dropped: ",
      paste(unmatched, collapse = ", "),
      ". Check the marker spelling in the panel against the control name."
    )
  }

  usable <- match_table[match_table$matched_by != "none", ]

  unstained_rows <- which(usable$channel == "unstained")
  if (length(unstained_rows) == 0) {
    stop(
      "No unstained control was found. flowStats::spillover_ng() needs one row ",
      "whose channel is 'unstained'."
    )
  }
  if (length(unstained_rows) > 1) {
    chosen <- if (is.null(keep_unstained)) {
      unstained_rows[1]
    } else {
      hit <- which(usable$filename == keep_unstained)
      if (length(hit) == 0) {
        stop("keep_unstained names a file that is not in the match table: ",
             keep_unstained)
      }
      hit[1]
    }
    dropped <- setdiff(unstained_rows, chosen)
    warning(
      "The set holds ", length(unstained_rows), " unstained controls. Keeping '",
      usable$filename[chosen], "' and dropping: ",
      paste(usable$filename[dropped], collapse = ", "), "."
    )
    usable <- usable[-dropped, ]
  }

  duplicated_channels <- usable$channel[
    duplicated(usable$channel) & usable$channel != "unstained"
  ]
  if (length(duplicated_channels) > 0) {
    stop(
      "More than one control matched the same channel: ",
      paste(unique(duplicated_channels), collapse = ", "),
      ". A spillover matrix needs one control per channel."
    )
  }

  utils::write.csv(
    usable[, c("filename", "channel")],
    path,
    row.names = FALSE,
    quote = TRUE
  )

  invisible(path)
}

#' Compute a spillover matrix from a set of single stained controls
#'
#' A thin wrapper over [flowStats::spillover_ng()] that keeps the match file and
#' the arguments together, so a report can state exactly how the matrix was made.
#'
#' @param flow_set A `flowSet` of compensation controls.
#' @param match_file Path to the CSV written by [WriteMatchFile()].
#' @param fsc The forward scatter channel used for the pregate.
#' @param ssc The side scatter channel used for the pregate.
#' @param method How the stained population is summarised. The default here is
#'   `"median"`, which is not the flowStats default of `"mode"`. See the note
#'   below, because the choice changes the result by more than 100 percentage
#'   points on a real dataset.
#' @param pregate Gate the control population before the matrix is computed. Keep
#'   this `TRUE`. With it off, the OMIP-39 controls produce 12 spillover values
#'   above 100 percent instead of none.
#' @return A numeric spillover matrix.
#'
#' @section Why the default is median and not mode:
#' `flowStats::spillover_ng()` defaults to `method = "mode"`. That works on bead
#' controls, where nearly every event is positive. It fails on cell controls for a
#' marker that only a minority of cells carry, because the mode then lands on the
#' negative population and the ratio it produces is meaningless.
#'
#' Measured on the OMIP-39 single stains, which are cells, against the matrix the
#' instrument stored for the same experiment:
#'
#' | method | correlation with stored | median difference | values above 100 percent |
#' | ------ | ----------------------- | ----------------- | ------------------------ |
#' | mode   | 0.619                   | 0.97 points       | 2                        |
#' | median | 0.975                   | 0.54 points       | 0                        |
#'
#' The two impossible values under `"mode"` came from NKG2C, positive on 7.0
#' percent of events, and NKG2A, positive on 17.2 percent. Run
#' [CheckControlQuality()] before you compute a matrix, so a control of that kind
#' is visible first.
#' @export
ComputeSpilloverFromControls <- function(flow_set,
                                         match_file,
                                         fsc = "FSC-A",
                                         ssc = "SSC-A",
                                         method = "median",
                                         pregate = TRUE) {
  if (!file.exists(match_file)) {
    stop("The match file does not exist: ", match_file)
  }

  available <- flowCore::colnames(flow_set[[1]])
  for (channel in c(fsc, ssc)) {
    if (!channel %in% available) {
      stop(
        "The scatter channel '", channel, "' is not in the data. ",
        "The file carries: ", paste(available, collapse = ", "), "."
      )
    }
  }

  # spillover_ng renames every sample of the set to the channel its match file
  # names, and it does not drop a sample that the match file leaves out. A set
  # that holds a second unstained control therefore reaches the naming step one
  # sample longer than the channel list, and the call fails with a length
  # mismatch. Subsetting to the match file first is what keeps the two in step.
  listed <- utils::read.csv(match_file, stringsAsFactors = FALSE)
  present <- flowCore::sampleNames(flow_set)
  missing <- setdiff(listed$filename, present)
  if (length(missing) > 0) {
    stop("The match file names ", length(missing),
         " file(s) that the set does not hold: ",
         paste(missing, collapse = ", "), ".")
  }
  flow_set <- flow_set[listed$filename]

  flowStats::spillover_ng(
    flow_set,
    fsc = fsc,
    ssc = ssc,
    matchfile = match_file,
    method = method,
    pregate = pregate,
    plot = FALSE
  )
}

#' Draw a spillover matrix as a heatmap
#'
#' FlowJo shows the matrix as a grid of numbers. A heatmap answers a different
#' question faster, which is where the large corrections sit. The diagonal is
#' dropped, because it is 1 everywhere and it would set the colour scale.
#'
#' @param spillover A spillover matrix.
#' @param title The plot title.
#' @param max_percent The top of the colour scale, as a percentage. Values above
#'   it are drawn at the top colour.
#' @return A `ggplot` object.
#' @export
PlotSpilloverHeatmap <- function(spillover,
                                 title = "Spillover matrix",
                                 max_percent = 25) {
  if (!is.matrix(spillover)) {
    stop("spillover must be a matrix, not a ", class(spillover)[1], ".")
  }

  long <- SummariseSpillover(spillover, top = nrow(spillover) * ncol(spillover))
  long$percent <- 100 * long$spill
  long$label <- ifelse(long$percent >= 1, sprintf("%.0f", long$percent), "")

  ggplot2::ggplot(
    long,
    ggplot2::aes(x = to, y = from, fill = pmin(percent, max_percent))
  ) +
    ggplot2::geom_tile(colour = "white", linewidth = 0.3) +
    ggplot2::geom_text(ggplot2::aes(label = label), size = 2.6,
                       colour = "grey15", family = FigureFont()) +
    ggplot2::scale_fill_gradient(
      low = "#f7f7f7", high = "#b2182b",
      name = "Spillover,\npercent",
      limits = c(0, max_percent),
      guide = ColourbarGuide()
    ) +
    ggplot2::labs(
      title = title,
      subtitle = paste(
        "The diagonal is dropped. A value is the percentage of the row",
        "detector\nthat the column detector also reads."
      ),
      x = "Into detector",
      y = "From detector"
    ) +
    ThemePublication() +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 90, hjust = 1, vjust = 0.5),
      panel.grid = ggplot2::element_blank()
    )
}

#' Compare a computed matrix against the one the instrument stored
#'
#' Use this when the files were compensated on the cytometer and you want to know
#' whether recomputing the matrix would change the result. A large difference on
#' one pair points at a control that was weak, wrongly matched, or acquired with a
#' different voltage.
#'
#' @param computed The matrix from [ComputeSpilloverFromControls()].
#' @param stored The matrix from [ExtractSpillover()].
#' @param top The number of rows to return, sorted by the size of the difference.
#' @return A `data.frame` with the columns `from`, `to`, `computed`, `stored` and
#'   `difference`, all as percentages, over the detectors the two matrices share.
#' @export
CompareSpilloverMatrices <- function(computed, stored, top = 20) {
  shared <- intersect(colnames(computed), colnames(stored))
  if (length(shared) == 0) {
    stop(
      "The two matrices share no detector name. Computed: ",
      paste(colnames(computed), collapse = ", "), ". Stored: ",
      paste(colnames(stored), collapse = ", "), "."
    )
  }

  computed_long <- SummariseSpillover(
    computed[shared, shared, drop = FALSE],
    top = length(shared)^2
  )
  stored_long <- SummariseSpillover(
    stored[shared, shared, drop = FALSE],
    top = length(shared)^2
  )

  merged <- merge(
    computed_long, stored_long,
    by = c("from", "to"),
    suffixes = c("_computed", "_stored")
  )

  result <- data.frame(
    from = merged$from,
    to = merged$to,
    computed = 100 * merged$spill_computed,
    stored = 100 * merged$spill_stored,
    stringsAsFactors = FALSE
  )
  result$difference <- result$computed - result$stored
  result <- result[order(-abs(result$difference)), ]
  rownames(result) <- NULL

  utils::head(result, top)
}

#' Check whether each single stained control can support a spillover estimate
#'
#' A spillover value is a ratio between the signal a fluorochrome puts into its
#' own detector and the signal it puts into another one. Both terms need a
#' positive population to measure. A bead control is almost entirely positive and
#' always supports the estimate. A cell control does not, because a marker that
#' sits on a minority of cells leaves most events negative.
#'
#' This is the check that explains the failure recorded in
#' [ComputeSpilloverFromControls()]. NKG2C is positive on 7.0 percent of the
#' OMIP-39 control events and NKG2A on 17.2 percent, and those two controls are
#' exactly the pair that produced spillover above 100 percent under the flowStats
#' default. Run this first and the problem is visible before the matrix is built.
#'
#' @param flow_set A `flowSet` of compensation controls.
#' @param match_table The output of [MatchControlsToChannels()].
#' @param quantile_cut The quantile of the unstained control that sets the
#'   positive threshold. The default of 0.999 is deliberately strict, so a dim
#'   spread in the negative population is not counted as positive.
#' @param min_positive_percent Below this percentage of positive events, the
#'   control is reported as weak.
#' @return A `data.frame` with one row per stained control, holding `filename`,
#'   `stain`, `channel`, `positive_percent`, `positive_events`,
#'   `brightest_channel`, `primary_is_brightest` and `verdict`. `verdict` is
#'   `"ok"`, `"weak"` when the positive fraction is below the threshold, or
#'   `"wrong channel"` when another detector reads the positive events more
#'   brightly than the assigned one.
#' @export
CheckControlQuality <- function(flow_set,
                                match_table,
                                quantile_cut = 0.999,
                                min_positive_percent = 20) {
  unstained_files <- match_table$filename[match_table$channel == "unstained" &
                                            !is.na(match_table$channel)]
  if (length(unstained_files) == 0) {
    stop("The match table names no unstained control, so no threshold can be set.")
  }

  unstained <- flowCore::exprs(flow_set[[unstained_files[1]]])
  fluorescence <- FluorescenceChannels(flow_set[[1]])

  stained <- match_table[!is.na(match_table$channel) &
                           match_table$channel != "unstained", ]
  if (nrow(stained) == 0) {
    stop("The match table holds no stained control.")
  }

  rows <- lapply(seq_len(nrow(stained)), function(i) {
    values <- flowCore::exprs(flow_set[[stained$filename[i]]])
    primary <- stained$channel[i]

    threshold <- stats::quantile(unstained[, primary], quantile_cut)
    is_positive <- values[, primary] > threshold
    positive_events <- sum(is_positive)
    positive_percent <- 100 * mean(is_positive)

    if (positive_events > 50) {
      medians <- vapply(
        fluorescence,
        function(channel) stats::median(values[is_positive, channel]),
        numeric(1)
      )
      brightest <- names(which.max(medians))
    } else {
      brightest <- NA_character_
    }

    primary_is_brightest <- !is.na(brightest) && brightest == primary

    verdict <- if (is.na(brightest)) {
      "weak"
    } else if (!primary_is_brightest) {
      "wrong channel"
    } else if (positive_percent < min_positive_percent) {
      "weak"
    } else {
      "ok"
    }

    data.frame(
      filename = stained$filename[i],
      stain = stained$stain[i],
      channel = primary,
      positive_percent = positive_percent,
      positive_events = positive_events,
      brightest_channel = brightest,
      primary_is_brightest = primary_is_brightest,
      verdict = verdict,
      stringsAsFactors = FALSE
    )
  })

  result <- do.call(rbind, rows)
  result[order(result$positive_percent), ]
}

#' Find the controls that the spillover gate cannot resolve
#'
#' [flowStats::spillover_ng()] gates each stained control on its own detector
#' with [flowStats::rangeGate()], asking for two peaks. A marker that only a
#' minority of the control cells carry does not present two peaks, the gate
#' fails, and the whole matrix fails with it. One failure out of 28 is enough to
#' return nothing at all.
#'
#' This runs the same two attempts that `spillover_ng` runs, so a control that
#' passes here passes there. Testing first turns a fatal error into a named list
#' of controls to leave out.
#'
#' @param flow_set A `flowSet` of compensation controls.
#' @param match_table The output of [MatchControlsToChannels()], or of a deposit
#'   specific wrapper around it.
#' A gate that succeeds can still leave a population too small to summarise.
#' `spillover()` fits a robust bivariate centre on whatever the gate returned,
#' and that fit needs more than three events. `min_events` rejects those cases
#' as well, because they fail later and with a message that names no control.
#'
#' @param scale_factor The width of the scatter pregate.
#' @param min_events The smallest gated population that is accepted.
#' @return A `data.frame` with the columns `filename`, `stain`, `channel`,
#'   `gated_events` and `gateable`.
#' @examples
#' \dontrun{
#' GateableControls(controls, matches)
#' }
#' @export
GateableControls <- function(flow_set, match_table, scale_factor = 2,
                             min_events = 500) {
  usable <- match_table[match_table$matched_by != "none" &
                          match_table$channel != "unstained", ]
  if (nrow(usable) == 0) {
    stop("No stained control was matched to a channel, so none can be tested.")
  }

  # The order below is the order `spillover_ng` uses. It gates on the stained
  # detector first and applies the scatter filter afterwards, so a test that
  # filters first passes controls that the real call then rejects.
  gated_events <- vapply(seq_len(nrow(usable)), function(index) {
    frame <- flow_set[[usable$filename[index]]]
    channel <- usable$channel[index]
    gate <- try(flowStats::rangeGate(frame, stain = channel, inBetween = TRUE,
                                     borderQuant = 0, absolute = FALSE,
                                     peakNr = 2), silent = TRUE)
    if (methods::is(gate, "try-error")) {
      gate <- try(flowStats::rangeGate(frame, stain = channel,
                                       inBetween = TRUE, borderQuant = 1,
                                       absolute = FALSE, peakNr = 2),
                  silent = TRUE)
    }
    if (methods::is(gate, "try-error")) {
      return(0L)
    }
    gated <- try(flowCore::Subset(frame, gate), silent = TRUE)
    if (methods::is(gated, "try-error")) {
      return(0L)
    }
    filtered <- try(flowCore::Subset(
      gated, flowStats::norm2Filter("FSC-A", "SSC-A",
                                    scale.factor = scale_factor)
    ), silent = TRUE)
    if (methods::is(filtered, "try-error")) {
      return(0L)
    }
    nrow(flowCore::exprs(filtered))
  }, integer(1))

  data.frame(
    filename = usable$filename, stain = usable$stain, channel = usable$channel,
    gated_events = gated_events, gateable = gated_events >= min_events,
    stringsAsFactors = FALSE
  )
}
