# Multicentre harmonisation.
#
# FlowRepository FR-FCM-Z282 holds 234 sample files from 13 operators at 5
# centres, acquired on 7 instruments of 3 models. Every earlier report in this
# repository gated one panel on one instrument, so a channel could be named once
# in a template CSV and reused. That is not possible here.
#
# The deposit carries 31 distinct channel profiles. Three makers name the same
# detector three ways, one operator configured the analyser by laser and filter,
# 52 files carry no marker name on at least one channel, and the sixth colour
# means a viability dye in the frozen PBMC panel and CD45 in the whole blood
# panel. A template that names `PE-Cy7-A` gates some of the deposit and fails
# silently on the rest.
#
# `ResolveMarkerChannels()` therefore resolves a channel per file, in two passes.
# The first pass reads the marker name that the operator typed. The second pass
# falls back on the fluorochrome, through the table in `gating/z282_panel.csv`.
# Every resolution records which pass produced it, so a weaker match is visible.
#
# The statistics in this file are the four indicators that the paper defines.
# They are the coefficient of variation, the bias against a reference value, the
# intraclass correlation, and the Z-score. Each one is a separate function, so a
# test asserts a number rather than the absence of an error.

# The order decides which marker a label resolves to when two tokens both match.
# `CD45RA` must come before `CD45`, and `CD45` before `CD4`, because a prefix
# match on `CD4` would otherwise claim the string `CD45RA APC`.
kMarkerOrder <- c("CD45RA", "CCR7", "CD45", "CD3", "CD8", "CD4", "NIR")

#' Read the panel table of the FR-FCM-Z282 deposit
#'
#' The table names the marker, the fluorochrome, the material the marker belongs
#' to, and every detector name that has been seen for that fluorochrome across
#' the three instrument makers.
#'
#' @param path Path to `gating/z282_panel.csv`.
#' @return A `data.frame` with the columns `marker`, `fluorochrome`, `material`,
#'   `detector_aliases` and `note`.
#' @examples
#' \dontrun{
#' ReadZ282Panel("gating/z282_panel.csv")
#' }
#' @export
ReadZ282Panel <- function(path) {
  if (!file.exists(path)) {
    stop("The panel table does not exist: ", path)
  }
  panel <- utils::read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
  required <- c("marker", "fluorochrome", "material", "detector_aliases")
  missing <- setdiff(required, names(panel))
  if (length(missing) > 0) {
    stop("The panel table lacks these columns: ", paste(missing, collapse = ", "))
  }
  panel
}

# Split a label into lower case words, with every separator treated as a space.
# `CD8 PE-Cy5.5` becomes c("cd8", "pe", "cy5", "5") and `CD3PC7` becomes
# c("cd3pc7"), which is why the caller also accepts a prefix match.
SplitLabel. <- function(label) {
  if (is.na(label) || !nzchar(label)) {
    return(character())
  }
  words <- strsplit(tolower(gsub("[^A-Za-z0-9]+", " ", label)), " +")[[1]]
  words[nzchar(words)]
}

# Return the marker that a `$PnS` label names, or NA when it names none.
MarkerFromLabel. <- function(label, order = kMarkerOrder) {
  words <- SplitLabel.(label)
  if (length(words) == 0) {
    return(NA_character_)
  }
  for (marker in order) {
    token <- tolower(marker)
    if (any(words == token | startsWith(words, token))) {
      return(marker)
    }
  }
  NA_character_
}

# Strip the `-A`, `-H` or `-W` suffix and every separator, so `PerCP-Cy5-5-A`
# and `PerCP-Cy5.5-A` compare equal.
NormaliseDetector. <- function(detector) {
  tolower(gsub("[^A-Za-z0-9]", "", sub("-[AHW]$", "", detector, ignore.case = TRUE)))
}

#' Match every marker of the panel to a channel of one file
#'
#' The first pass reads the marker name in `$PnS`. The second pass falls back on
#' the detector name, which is what rescues a file whose operator left `$PnS`
#' empty. Only area channels take part, because a height channel and a width
#' channel carry the same marker and would double the match.
#'
#' @param frame A `flowFrame`.
#' @param panel The table returned by [ReadZ282Panel()].
#' @param material Either `"PBMC"` or `"WB"`. It decides whether the sixth colour
#'   is read as the viability dye or as CD45.
#' @return A `data.frame` with one row per expected marker and the columns
#'   `marker`, `channel` and `resolved_by`. `channel` is `NA` when neither pass
#'   found it.
#' @examples
#' \dontrun{
#' ResolveMarkerChannels(frame, panel, "WB")
#' }
#' @export
ResolveMarkerChannels <- function(frame, panel, material) {
  if (!material %in% c("PBMC", "WB")) {
    stop("material must be PBMC or WB, not: ", material)
  }
  wanted <- panel[panel$material == "both" | panel$material == material, ]
  if (nrow(wanted) == 0) {
    stop("The panel table names no marker for the material: ", material)
  }

  parameters <- flowCore::pData(flowCore::parameters(frame))
  is_area <- grepl("-A$", parameters$name, ignore.case = TRUE)
  detectors <- parameters$name[is_area]
  labels <- parameters$desc[is_area]

  result <- data.frame(
    marker = wanted$marker,
    channel = NA_character_,
    resolved_by = NA_character_,
    stringsAsFactors = FALSE
  )

  # Pass 1. Read the label that the operator typed.
  for (index in seq_along(detectors)) {
    marker <- MarkerFromLabel.(labels[index])
    if (is.na(marker)) {
      next
    }
    row <- match(marker, result$marker)
    if (is.na(row) || !is.na(result$channel[row])) {
      next
    }
    result$channel[row] <- detectors[index]
    result$resolved_by[row] <- "marker"
  }

  # Pass 2. Fall back on the fluorochrome that the detector name carries.
  normalised <- NormaliseDetector.(detectors)
  for (row in seq_len(nrow(result))) {
    if (!is.na(result$channel[row])) {
      next
    }
    aliases <- strsplit(wanted$detector_aliases[row], "|", fixed = TRUE)[[1]]
    aliases <- tolower(gsub("[^A-Za-z0-9]", "", aliases))
    hit <- which(normalised %in% aliases & !(detectors %in% result$channel))
    if (length(hit) > 0) {
      result$channel[row] <- detectors[hit[1]]
      result$resolved_by[row] <- "detector"
    }
  }

  result
}

#' Name the scatter and time channels of one file
#'
#' Beckman Coulter writes `FS-A` and `SS-A` where Becton Dickinson writes `FSC-A`
#' and `SSC-A`, so a fixed name fails on a third of the deposit.
#'
#' @param frame A `flowFrame`.
#' @return A named character vector with the elements `forward_area`,
#'   `forward_height`, `side_area` and `time`.
#' @examples
#' \dontrun{
#' ScatterChannels(frame)
#' }
#' @export
ScatterChannels <- function(frame) {
  names_present <- flowCore::colnames(frame)
  Pick <- function(pattern, label) {
    hit <- grep(pattern, names_present, ignore.case = TRUE, value = TRUE)
    if (length(hit) == 0) {
      stop("This file carries no ", label, " channel")
    }
    hit[1]
  }
  c(
    forward_area = Pick("^FSC?-A$", "forward scatter area"),
    forward_height = Pick("^FSC?-H$", "forward scatter height"),
    side_area = Pick("^SSC?-A$", "side scatter area"),
    time = Pick("^time$", "time")
  )
}

#' Report and settle the compensation state of one file
#'
#' The deposit holds three states. A Becton Dickinson export carries `SPILL` with
#' `APPLY COMPENSATION` set to `TRUE`, which means the values are compensated
#' already. A Kaluza export carries `$SPILLOVER` and no such keyword, so the
#' matrix still has to be applied. Eighteen files carry no matrix at all.
#'
#' @param frame A `flowFrame`.
#' @return A list with `frame`, which is compensated when a matrix had to be
#'   applied, and `state`, which is one of `"already compensated"`,
#'   `"matrix applied"` or `"no matrix"`.
#' @examples
#' \dontrun{
#' SettleCompensation(frame)$state
#' }
#' @export
SettleCompensation <- function(frame) {
  keywords <- flowCore::keyword(frame)
  applied <- keywords[["APPLY COMPENSATION"]]
  already <- !is.null(applied) && toupper(trimws(applied)) == "TRUE"

  matrix_keyword <- intersect(c("$SPILLOVER", "SPILL", "$SPILL"), names(keywords))
  if (length(matrix_keyword) == 0) {
    return(list(frame = frame, state = "no matrix"))
  }
  if (already) {
    return(list(frame = frame, state = "already compensated"))
  }

  spillover <- keywords[[matrix_keyword[1]]]
  if (!is.matrix(spillover) || nrow(spillover) == 0) {
    return(list(frame = frame, state = "no matrix"))
  }
  list(
    frame = flowCore::compensate(frame, spillover),
    state = "matrix applied"
  )
}

#' Keep the stable part of the acquisition
#'
#' The study SOP asks the operator to gate the window of stable acquisition
#' before anything else. This function does the same by a rule rather than by
#' eye. It divides the run into equal time bins, takes the median number of
#' events per bin, and keeps the longest run of consecutive bins whose count
#' stays inside `tolerance`.
#'
#' Two guards keep the rule from doing harm.
#'
#' The first is the bin count. One operator exported a time channel with a range
#' of 83 units, so 100 bins leave 16 of them empty and the longest stable run
#' collapses. The bin count halves until fewer than a tenth of the bins are
#' empty.
#'
#' The second is `min_kept`. A smooth change of acquisition rate is not a fluidic
#' disturbance, and the rule cannot tell the two apart. When the window would
#' discard more than `1 - min_kept` of the events the function keeps the whole
#' file and reports `applied = FALSE`, so the report can count how often that
#' happened rather than lose a quarter of a file in silence.
#'
#' @param frame A `flowFrame`.
#' @param time_channel The name of the time channel.
#' @param bins The starting number of time bins. Defaults to 100.
#' @param tolerance A length two vector of multipliers on the median bin count. A
#'   bin passes when its count lies between them.
#' @param min_kept The smallest fraction of events the window may keep before it
#'   is abandoned.
#' @return A list with `frame`, the filtered `flowFrame`, `kept`, the fraction of
#'   events that survived, and `applied`, whether the window was used.
#' @examples
#' \dontrun{
#' StableTimeWindow(frame, "Time")$kept
#' }
#' @export
StableTimeWindow <- function(frame, time_channel, bins = 100,
                             tolerance = c(0.5, 2), min_kept = 0.75) {
  values <- flowCore::exprs(frame)[, time_channel]
  total <- length(values)
  if (total == 0 || diff(range(values)) == 0) {
    return(list(frame = frame, kept = 1, applied = FALSE))
  }

  repeat {
    breaks <- seq(min(values), max(values), length.out = bins + 1)
    bin_of <- cut(values, breaks = breaks, include.lowest = TRUE, labels = FALSE)
    counts <- tabulate(bin_of, nbins = bins)
    if (sum(counts == 0) <= 0.1 * bins || bins <= 20) {
      break
    }
    bins <- max(20, floor(bins / 2))
  }

  reference <- stats::median(counts[counts > 0])
  passes <- counts >= reference * tolerance[1] & counts <= reference * tolerance[2]

  runs <- rle(passes)
  if (!any(runs$values)) {
    return(list(frame = frame, kept = 1, applied = FALSE))
  }
  best <- which(runs$values)[which.max(runs$lengths[runs$values])]
  first <- if (best == 1) 1 else sum(runs$lengths[seq_len(best - 1)]) + 1
  last <- first + runs$lengths[best] - 1

  keep <- !is.na(bin_of) & bin_of >= first & bin_of <= last
  fraction <- sum(keep) / total
  if (fraction < min_kept) {
    return(list(frame = frame, kept = 1, applied = FALSE))
  }
  list(frame = frame[keep, ], kept = fraction, applied = TRUE)
}

#' Coefficient of variation
#'
#' @param x A numeric vector.
#' @return The standard deviation divided by the mean. `NA` when the mean is zero
#'   or when fewer than two finite values are present.
#' @examples
#' CoefficientOfVariation(c(10, 12, 11))
#' @export
CoefficientOfVariation <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) < 2) {
    return(NA_real_)
  }
  average <- mean(x)
  if (average == 0) {
    return(NA_real_)
  }
  stats::sd(x) / average
}

#' Relative bias against a reference value
#'
#' The paper defines bias as the absolute difference between an operator's mean
#' and the reference value, divided by the reference value.
#'
#' @param x A numeric vector of one operator's values.
#' @param reference The reference value.
#' @return The absolute relative bias. `NA` when the reference is zero.
#' @examples
#' RelativeBias(c(9, 11), reference = 10)
#' @export
RelativeBias <- function(x, reference) {
  x <- x[is.finite(x)]
  if (length(x) == 0 || !is.finite(reference) || reference == 0) {
    return(NA_real_)
  }
  abs((mean(x) - reference) / reference)
}

#' Z-score against a reference distribution
#'
#' @param x A numeric vector of observations to score.
#' @param mu The mean of the reference distribution.
#' @param sigma The standard deviation of the reference distribution.
#' @return A numeric vector of the same length as `x`. `NA` when `sigma` is zero
#'   or not finite.
#' @examples
#' ZScore(c(12, 8), mu = 10, sigma = 2)
#' @export
ZScore <- function(x, mu, sigma) {
  if (!is.finite(sigma) || sigma == 0) {
    return(rep(NA_real_, length(x)))
  }
  (x - mu) / sigma
}

#' Intraclass correlation from a two-way random effects model
#'
#' The paper defines agreement as the ratio of the biological variance, which is
#' the variance between donors, to the total variance. It fits a two-way design
#' with the analyst and the donor as random variables. This function fits
#' `value ~ (1 | donor) + (1 | analyst)` and returns that ratio.
#'
#' @param value A numeric vector, one value per analyst and donor.
#' @param analyst A vector that identifies the analyst.
#' @param donor A vector that identifies the donor.
#' @return The intraclass correlation. `NA` when the model cannot be fitted or
#'   when the total variance is zero.
#' @examples
#' \dontrun{
#' IntraclassCorrelation(results$percent, results$operator, results$vial)
#' }
#' @export
IntraclassCorrelation <- function(value, analyst, donor) {
  keep <- is.finite(value) & !is.na(analyst) & !is.na(donor)
  value <- value[keep]
  analyst <- as.factor(as.character(analyst[keep]))
  donor <- as.factor(as.character(donor[keep]))
  if (length(value) < 4 ||
      nlevels(analyst) < 2 ||
      nlevels(donor) < 2 ||
      stats::var(value) == 0) {
    return(NA_real_)
  }

  frame <- data.frame(value = value, analyst = analyst, donor = donor)
  model <- try(
    suppressMessages(suppressWarnings(
      lme4::lmer(value ~ (1 | donor) + (1 | analyst), data = frame)
    )),
    silent = TRUE
  )
  if (inherits(model, "try-error")) {
    return(NA_real_)
  }

  variances <- as.data.frame(lme4::VarCorr(model))
  between_donor <- variances$vcov[variances$grp == "donor"]
  total <- sum(variances$vcov)
  if (length(between_donor) == 0 || !is.finite(total) || total == 0) {
    return(NA_real_)
  }
  between_donor / total
}
