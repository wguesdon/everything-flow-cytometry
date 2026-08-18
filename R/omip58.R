# Gate the OMIP-058 deposit and hand the result to Python.
#
# FR-FCM-ZYRN holds two PBMC files acquired on a 28 colour panel, together with
# 59 compensation controls, 55 of them single stains. The panel is described in
# Liechti and Roederer, Cytometry A 2019;95(9):946-951.
#
# Every file of the deposit stores a spillover matrix that is exactly the
# identity, alongside the keyword `APPLY COMPENSATION` set to TRUE. The keyword
# is misleading. Measured inside the lymphocyte gate of the first donor, CD4 and
# CD8 correlate at 0.707 where they are mutually exclusive on a T cell, CD3 and
# CD8 at 0.938 on two adjacent violet detectors, and CD16 and the viability dye
# at 0.865 on two adjacent ultraviolet detectors where no biology joins them.
# The median absolute correlation over all 153 marker pairs is 0.582. The
# deposited values carry their spillover, and the identity matrix records that
# no matrix was supplied rather than that none is needed.
#
# The matrix is therefore computed from the deposited single stains, which is
# why the deposit ships 55 of them.
#
# The singlet gate here is one sided. [SingletMask()] in R/spectral.R fits the
# line of forward scatter height on area and drops the events far from it on
# either side, which is the rule the other reports use. This deposit keeps the
# events above the line, because a doublet lowers the ratio of height to area
# and nothing raises it, so trimming the upper tail removes the largest
# singlets for no reason. The two functions carry different names because
# tests/testthat.R sources every file in R/ into one environment, and the
# second definition of a name silently replaces the first.
#
# The gate hierarchy follows Figure 1A to 1D of the paper. It stops at the CD3
# split, because the populations below that split are what the Python side is
# asked to find without being told where they are.

#' The markers that the gate hierarchy needs, and the token that finds each one
#'
#' The deposit writes a marker as an antibody followed by a fluorochrome, for
#' example `TCR Vd1 FITC` or `CD161 PE-Cy5`. A token match on the antibody is
#' safer than a substring match, because `cd16` is a substring of `cd161` and a
#' substring rule silently resolves CD16 to the CD161 detector.
kOmip58Tokens <- c(
  viability = "live",
  CD3 = "cd3",
  CD4 = "cd4",
  CD8 = "cd8",
  CD16 = "cd16",
  CD56 = "cd56",
  CD161 = "cd161",
  HLADR = "hla",
  CCR7 = "ccr7",
  CD45RA = "cd45ra",
  CD27 = "cd27",
  CD28 = "cd28",
  CD95 = "cd95",
  tetramer = "tetramer",
  Va72 = "va7",
  Vd1 = "vd1",
  Vd2 = "vd2",
  Vg9 = "vg9"
)

#' The two populations that are written to disk for the Python side
kOmip58Handoff <- c("live_lymphocytes", "cd3_t_cells")


#' Resolve every marker of the OMIP-058 panel to its detector
#'
#' @param frame A `flowFrame` from the deposit.
#' @param tokens A named character vector of antibody tokens. Defaults to
#'   [kOmip58Tokens].
#' @return A `data.frame` with the columns `name`, `token`, `channel` and
#'   `marker`, in the order of `tokens`.
#' @examples
#' \dontrun{
#' ResolveOmip58Channels(frame)
#' }
#' @export
ResolveOmip58Channels <- function(frame, tokens = kOmip58Tokens) {
  panel <- DescribeChannels(frame)
  panel <- panel[panel$is_marker, ]
  if (nrow(panel) == 0) {
    stop("The frame carries no named marker, so no channel can be resolved.")
  }

  marker_tokens <- MarkerTokens(panel$marker)
  rows <- lapply(seq_along(tokens), function(index) {
    token <- tokens[[index]]
    hits <- which(vapply(marker_tokens, function(set) token %in% set,
                         logical(1)))
    if (length(hits) == 0) {
      stop("No marker of the panel carries the token '", token, "', which is ",
           "how ", names(tokens)[index], " is found. The panel holds: ",
           paste(panel$marker, collapse = ", "), ".")
    }
    if (length(hits) > 1) {
      stop("The token '", token, "' matches more than one marker: ",
           paste(panel$marker[hits], collapse = ", "),
           ". A token must resolve one detector.")
    }
    data.frame(
      name = names(tokens)[index], token = token,
      channel = panel$channel[hits], marker = panel$marker[hits],
      stringsAsFactors = FALSE
    )
  })

  do.call(rbind, rows)
}



#' Match the OMIP-058 compensation controls to their detectors
#'
#' [MatchControlsToChannels()] resolves 26 of the 27 antibody controls of this
#' deposit on its own. Two names defeat it, and both are named here rather than
#' left to a looser matching rule that would resolve them by accident.
#'
#' The viability control is filed as `LD UV Blue_stained` against a panel entry
#' of `Live Dead UV Blue`, so neither the spelling nor the antibody token
#' agrees. Its unstained partner is filed as `LD UV Blue_unstained`, which the
#' unstained rule already catches.
#'
#' @param flow_set A `flowSet` of compensation controls.
#' @param reference The sample whose panel supplies the marker names.
#' @return The `data.frame` of [MatchControlsToChannels()], with the viability
#'   control resolved and `matched_by` set to `"deposit"` on the rows this
#'   function fixed.
#' @examples
#' \dontrun{
#' MatchOmip58Controls(controls)
#' }
#' @export
MatchOmip58Controls <- function(flow_set, reference = 1) {
  matched <- MatchControlsToChannels(flow_set, reference = reference)
  panel <- DescribeChannels(flow_set[[reference]])
  panel <- panel[panel$is_marker, ]

  viability_channel <- panel$channel[
    vapply(MarkerTokens(panel$marker), function(set) "live" %in% set,
           logical(1))
  ]
  if (length(viability_channel) != 1) {
    stop("The panel holds ", length(viability_channel),
         " markers whose name carries the token 'live'. It must hold one.")
  }

  needs_fix <- matched$matched_by == "none" &
    grepl("^ld uv blue", NormaliseMarkerName(matched$stain)) &
    !grepl("unstain", matched$stain, ignore.case = TRUE)
  matched$channel[needs_fix] <- viability_channel
  matched$marker[needs_fix] <- panel$marker[panel$channel == viability_channel]
  matched$matched_by[needs_fix] <- "deposit"

  matched
}

#' Compute the spillover matrix of one OMIP-058 control set
#'
#' The cell controls are used rather than the bead controls, because the deposit
#' provides a viability control only as cells and a matrix that leaves the
#' viability detector out cannot correct the spillover between it and CD16.
#'
#' The summary method is the median and not the mode that
#' [flowStats::spillover_ng()] defaults to. A cell control for a marker that a
#' minority of cells carry puts the mode on the negative population, and the
#' ratio that follows is meaningless.
#'
#' @param control_dir The folder that holds the deposit.
#' @param pattern A regular expression that selects the control files.
#' @param match_path Where the match file is written.
#' @param keep_unstained The unstained control to keep when the set holds more
#'   than one.
#' @param subsample Read at most this many events per control. `NULL` reads all.
#' @return A list with `spillover`, the matrix, `matches`, the table that
#'   [MatchOmip58Controls()] produced, and `gateable`, the table that
#'   [GateableControls()] produced. A control that cannot be gated is marked
#'   unmatched, so it is left out of the matrix and named in `gateable` rather
#'   than stopping the computation.
#' @examples
#' \dontrun{
#' ComputeOmip58Spillover(dir, "^Comp_Cells", "output/omip58/match.csv")
#' }
#' @export
ComputeOmip58Spillover <- function(control_dir,
                                   pattern,
                                   match_path,
                                   keep_unstained = NULL,
                                   subsample = NULL) {
  files <- list.files(control_dir, pattern = pattern, full.names = TRUE)
  if (length(files) == 0) {
    stop("No control file in ", control_dir, " matches '", pattern, "'.")
  }
  controls <- if (is.null(subsample)) {
    flowCore::read.flowSet(files, truncate_max_range = FALSE,
                           alter.names = FALSE)
  } else {
    flowCore::read.flowSet(files, truncate_max_range = FALSE,
                           alter.names = FALSE, which.lines = subsample)
  }

  matches <- MatchOmip58Controls(controls)
  gateable <- GateableControls(controls, matches)
  dropped <- gateable$filename[!gateable$gateable]
  matches$matched_by[matches$filename %in% dropped] <- "none"

  dir.create(dirname(match_path), recursive = TRUE, showWarnings = FALSE)
  suppressWarnings(
    WriteMatchFile(matches, match_path, keep_unstained = keep_unstained)
  )

  spillover <- ComputeSpilloverFromControls(
    controls, match_path, fsc = "FSC-A", ssc = "SSC-A", method = "median"
  )
  list(spillover = spillover, matches = matches, gateable = gateable)
}







#' Gate one OMIP-058 PBMC file to the populations a cut can defend
#'
#' The hierarchy follows Figure 1A of the paper and stops at the CD3 split.
#' Scatter removes the out of range events, the doublets and the debris, the
#' viability channel removes the dead cells, and CD3 separates the T cells from
#' the natural killer cells. Each of those four boundaries is a separated pair of
#' modes that a fitted cut resolves.
#'
#' The hierarchy stops there on purpose. Inside CD3 positive cells the deposit's
#' rare markers are unimodal, the deposit carries no fluorescence minus one
#' control to place a threshold against, and a fitted cut lands inside the
#' negative population. [AssessOneDimensionalCuts()] records what each of those
#' cuts would have done. Finding those subsets is left to the clustering on the
#' Python side, which reads every marker at once instead of one at a time.
#'
#' @param path Path to the FCS file.
#' @param spillover The spillover matrix to apply, from
#'   [ComputeOmip58Spillover()]. The deposit stores no usable matrix of its own,
#'   so this argument is required.
#' @param cofactor The arcsinh cofactor used for every fitted cut.
#' @param min_events The smallest parent that a cut is attempted on.
#' @return A list with `counts`, a `data.frame` of one row per population,
#'   `cuts`, a `data.frame` of one row per fitted cut, `rare`, the output of
#'   [AssessOneDimensionalCuts()] inside the CD3 positive population, `masks`, a
#'   named list of logical vectors over the events of the file, `channels`,
#'   `state`, `frame` and `cofactor`.
#' @examples
#' \dontrun{
#' GateOmip58File(path, spillover)$counts
#' }
#' @export
GateOmip58File <- function(path, spillover, cofactor = 150,
                           min_events = 200) {
  frame <- flowCore::read.FCS(path, truncate_max_range = FALSE)

  # The argument is checked before any work is done on the file, so a matrix
  # that does not belong to this panel is reported as such rather than as a
  # missing marker.
  absent <- setdiff(rownames(spillover), flowCore::colnames(frame))
  if (length(absent) > 0) {
    stop("The spillover matrix names ", length(absent),
         " detector(s) that the file does not carry: ",
         paste(absent, collapse = ", "), ".")
  }

  state <- ReadCompensationState(frame)
  channels <- ResolveOmip58Channels(frame)
  scatter <- PanelScatterChannels(frame)
  frame <- flowCore::compensate(frame, spillover)

  values <- ArcsinhTransform(flowCore::exprs(frame), channels$channel, cofactor)
  total <- nrow(values)

  cuts <- list()
  Channel <- function(name) channels$channel[channels$name == name]

  # Fit one cut on one parent and lift the result back to the full event index,
  # so that every mask in this function indexes the same rows.
  Cut <- function(parent, name, keep, label) {
    empty <- rep(FALSE, total)
    if (sum(parent) < min_events) {
      cuts[[length(cuts) + 1]] <<- data.frame(
        label = label, marker = name, channel = Channel(name), cut = NA_real_,
        rule = "too few events", parent_events = sum(parent),
        stringsAsFactors = FALSE
      )
      return(empty)
    }
    split <- SplitOnChannel(values[parent, , drop = FALSE], Channel(name), keep)
    cuts[[length(cuts) + 1]] <<- data.frame(
      label = label, marker = name, channel = Channel(name), cut = split$cut,
      rule = split$rule, parent_events = sum(parent), stringsAsFactors = FALSE
    )
    if (split$rule == "none") {
      return(empty)
    }
    empty[which(parent)] <- split$mask
    empty
  }

  in_range <- InScatterRange(values, scatter)
  singlets <- in_range
  singlets[in_range] <- RatioSingletMask(values[in_range, , drop = FALSE],
                                          scatter)
  lymphocytes <- singlets
  lymphocytes[singlets] <- LymphocyteMask(values[singlets, , drop = FALSE],
                                          scatter)

  viable <- Cut(lymphocytes, "viability", "below", "viable")
  cd3 <- Cut(viable, "CD3", "above", "CD3 positive")
  not_cd3 <- viable & !cd3

  rare <- AssessOneDimensionalCuts(
    values, channels, cd3,
    c("tetramer", "Vd1", "Vd2", "Vg9", "Va72", "CD161")
  )
  rare_negative <- AssessOneDimensionalCuts(
    values, channels, not_cd3, c("CD16", "CD56", "HLADR")
  )
  rare$parent <- "cd3_t_cells"
  rare_negative$parent <- "cd3_negative"
  rare <- rbind(rare, rare_negative)

  hierarchy <- list(
    list("all_events", rep(TRUE, total), NA_character_),
    list("in_scatter_range", in_range, "all_events"),
    list("singlets", singlets, "in_scatter_range"),
    list("lymphocytes", lymphocytes, "singlets"),
    list("live_lymphocytes", viable, "lymphocytes"),
    list("cd3_t_cells", cd3, "live_lymphocytes"),
    list("cd3_negative", not_cd3, "live_lymphocytes")
  )

  masks <- stats::setNames(lapply(hierarchy, `[[`, 2),
                           vapply(hierarchy, `[[`, character(1), 1))
  counts <- do.call(rbind, lapply(hierarchy, function(row) {
    parent <- row[[3]]
    parent_events <- if (is.na(parent)) total else sum(masks[[parent]])
    data.frame(
      population = row[[1]], parent = parent, events = sum(row[[2]]),
      parent_events = parent_events,
      percent_of_parent = if (parent_events == 0) NA_real_ else
        100 * sum(row[[2]]) / parent_events,
      percent_of_lymphocytes = if (sum(masks[["live_lymphocytes"]]) == 0)
        NA_real_ else
        100 * sum(row[[2]]) / sum(masks[["live_lymphocytes"]]),
      stringsAsFactors = FALSE
    )
  }))

  list(
    counts = counts, cuts = do.call(rbind, cuts), rare = rare, masks = masks,
    channels = channels, state = state, frame = frame, cofactor = cofactor
  )
}

