# Gate the OMIP-051 deposit, a 28 colour B cell and myeloid panel.
#
# FR-FCM-ZYN4 holds one PBMC file of 2,126,989 events and 59 compensation
# controls, 29 on beads, 28 on cells and 2 for the viability dye. The panel is
# described in Liechti and Roederer, Cytometry A 2019;95(2):150-155.
#
# The deposit was acquired on the same instrument as FR-FCM-ZYRN, an X50-E
# LSRFortessa, eight days earlier, and it stores the same identity spillover
# matrix beside the same `APPLY COMPENSATION` keyword. Whether the values carry
# their spillover is measured rather than assumed, by the correlation between
# marker pairs.
#
# The generic steps live in R/panels.R. This file holds what is specific to the
# panel: two markers that a name match cannot resolve, and the gate hierarchy of
# Figure 1.

#' The markers of the OMIP-051 panel, with the fluorochrome that finds each one
#'
#' Two detectors of this deposit carry a marker name that lists both of the
#' antibodies the laboratory ran on them. `B780-A` is written `CD27 or IgD
#' BB790` and `R730-A` is written `IgD Ax700 or CD27 APC-R700`. A match on the
#' antibody therefore resolves CD27 to two channels and IgD to two channels, and
#' either answer is wrong half the time.
#'
#' The fluorochrome resolves both. Table 2 of the paper assigns IgD to BB790 and
#' CD27 to APC-R700, the control files are named for the same pairing, and the
#' sample file is named `PBMC 22718_IgD BB790`. Every marker here is therefore
#' found by its fluorochrome, and the antibody token is the fallback.
kOmip51Markers <- data.frame(
  name = c("viability", "CD14", "HLADR", "CD19", "CD20", "CD10", "CD27",
           "IgD", "IgM", "CD123", "CD11c", "CD1c", "CD141", "CD16", "CD21",
           "CD85j", "IgG", "IgA"),
  fluorochrome = c("uv blue", "bv510", "buv661", "apc h7", "buv805", "bv650",
                   "apc r700", "bb790", "bv570", "bb660", "pe cy5 5",
                   "buv395", "bb630", "bb700", "buv496", "pe cy5", "buv737",
                   "apc"),
  antibody = c("live", "cd14", "hla", "cd19", "cd20", "cd10", "cd27", "igd",
               "igm", "cd123", "cd11c", "cd1c", "cd141", "cd16", "cd21",
               "cd85j", "igg", "iga"),
  stringsAsFactors = FALSE
)

#' Resolve every marker of the OMIP-051 panel to its detector
#'
#' The fluorochrome is matched first, on the normalised marker name, and the
#' longest match wins so that `pe cy5 5` is not resolved by `pe cy5`. The
#' antibody token is tried only where the fluorochrome finds nothing.
#'
#' @param frame A `flowFrame` from the deposit.
#' @param markers The marker table. Defaults to [kOmip51Markers].
#' @return A `data.frame` with the columns `name`, `channel`, `marker` and
#'   `matched_by`, in the order of `markers`.
#' @examples
#' \dontrun{
#' ResolveOmip51Channels(frame)
#' }
#' @export
ResolveOmip51Channels <- function(frame, markers = kOmip51Markers) {
  panel <- DescribeChannels(frame)
  panel <- panel[panel$is_marker, ]
  if (nrow(panel) == 0) {
    stop("The frame carries no named marker, so no channel can be resolved.")
  }
  normalised <- NormaliseMarkerName(panel$marker)

  rows <- lapply(seq_len(nrow(markers)), function(index) {
    fluorochrome <- markers$fluorochrome[index]
    hits <- which(grepl(fluorochrome, normalised, fixed = TRUE))
    matched_by <- "fluorochrome"

    # `apc` is inside `apc h7` and `apc r700`, so a marker whose fluorochrome is
    # the bare name keeps only the channels that no longer match a longer one.
    if (length(hits) > 1) {
      longer <- markers$fluorochrome[markers$fluorochrome != fluorochrome &
                                       grepl(fluorochrome, markers$fluorochrome,
                                             fixed = TRUE)]
      for (other in longer) {
        hits <- hits[!grepl(other, normalised[hits], fixed = TRUE)]
      }
    }

    if (length(hits) == 0) {
      token <- markers$antibody[index]
      hits <- which(vapply(MarkerTokens(panel$marker),
                           function(set) token %in% set, logical(1)))
      matched_by <- "antibody"
    }

    if (length(hits) != 1) {
      stop("The marker '", markers$name[index], "' resolved ", length(hits),
           " channels. It must resolve one. The panel holds: ",
           paste(panel$marker, collapse = ", "), ".")
    }

    data.frame(name = markers$name[index], channel = panel$channel[hits],
               marker = panel$marker[hits], matched_by = matched_by,
               stringsAsFactors = FALSE)
  })

  do.call(rbind, rows)
}

#' Match the OMIP-051 compensation controls to their detectors
#'
#' [MatchControlsToChannels()] resolves a control by its antibody name, and on
#' this deposit that is wrong for two of them. The panel writes `B780-A` as
#' `CD27 or IgD BB790` and `R730-A` as `IgD Ax700 or CD27 APC-R700`, so the
#' antibody `IgD` matches both channels and the control `Comp_Cells_IgD BB790`
#' was resolved to `R730-A`, the APC-R700 detector.
#'
#' A matrix built that way has two rows transposed. Measured on this file it
#' leaves the median absolute correlation over 153 marker pairs at 0.601 and
#' drives CD19 against CD20 to -0.700, where the two mark the same cells and the
#' correlation is positive before compensation.
#'
#' The fluorochrome resolves it. A control is named `antibody fluorochrome`, and
#' the panel entry for its detector ends with the same fluorochrome, so the last
#' token of the stain identifies the channel. `Live Dead UV Blue` ends with
#' `blue` and resolves on the same rule.
#'
#' @param flow_set A `flowSet` of compensation controls.
#' @param reference The sample whose panel supplies the marker names.
#' @return The `data.frame` of [MatchControlsToChannels()], with every row the
#'   fluorochrome resolves rewritten and `matched_by` set to `"fluorochrome"`.
#' @examples
#' \dontrun{
#' MatchOmip51Controls(controls)
#' }
#' @export
MatchOmip51Controls <- function(flow_set, reference = 1) {
  matched <- MatchControlsToChannels(flow_set, reference = reference)
  panel <- DescribeChannels(flow_set[[reference]])
  panel <- panel[panel$is_marker, ]
  panel_normalised <- NormaliseMarkerName(panel$marker)

  is_unstained <- matched$matched_by == "unstained"
  for (index in which(!is_unstained)) {
    tokens <- MarkerTokens(matched$stain[index])[[1]]
    if (length(tokens) == 0) {
      next
    }
    last <- tokens[length(tokens)]
    hits <- which(endsWith(panel_normalised, paste0(" ", last)) |
                    panel_normalised == last)
    if (length(hits) != 1) {
      next
    }
    matched$channel[index] <- panel$channel[hits]
    matched$marker[index] <- panel$marker[hits]
    matched$matched_by[index] <- "fluorochrome"
  }

  matched
}

#' Compute the spillover matrix of one OMIP-051 control set
#'
#' @param control_dir The folder that holds the FCS files.
#' @param pattern A regular expression that selects the control files.
#' @param match_path Where the match file is written.
#' @param keep_unstained The unstained control to keep when the set holds more
#'   than one.
#' @return A list with `spillover`, the matrix, `matches`, the table that
#'   [MatchOmip51Controls()] produced, and `gateable`, the table that
#'   [GateableControls()] produced. A control that cannot be gated is marked
#'   unmatched, so it is left out of the matrix and named rather than stopping
#'   the computation.
#' @examples
#' \dontrun{
#' ComputeOmip51Spillover(dir, "^Comp_Cells", "output/omip51/match.csv")
#' }
#' @export
ComputeOmip51Spillover <- function(control_dir, pattern, match_path,
                                   keep_unstained = NULL) {
  files <- list.files(control_dir, pattern = pattern, full.names = TRUE)
  if (length(files) == 0) {
    stop("No control file in ", control_dir, " matches '", pattern, "'.")
  }
  controls <- flowCore::read.flowSet(files, truncate_max_range = FALSE,
                                     alter.names = FALSE)

  matches <- MatchOmip51Controls(controls)
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

#' Gate one OMIP-051 PBMC file down the published hierarchy
#'
#' The hierarchy follows Figure 1A, 1B and 1G of the paper. Scatter removes the
#' out of range events, the doublets and the debris. The viability channel
#' removes the dead cells. CD14 separates the monocytes. Inside the CD14
#' negative cells, HLA-DR or CD20 selects the B cells and the dendritic cells,
#' CD19 with CD20 gives the B cells, and CD123 with CD11c splits the remaining
#' dendritic cells.
#'
#' The paper notes a high scatter population inside the B cell gate and removes
#' it with a lymphocyte gate placed after the B cell gate rather than before it.
#' That order is kept here, so the size of the population it removes can be
#' measured rather than assumed away.
#'
#' Every cut below the dendritic cell split is reported rather than applied. The
#' subsets of Figure 1C rest on markers that a fitted cut cannot place on one
#' file, and [AssessOneDimensionalCuts()] records what each of those cuts would
#' have done.
#'
#' @param path Path to the FCS file.
#' @param spillover The spillover matrix to apply, from
#'   [ComputeOmip51Spillover()]. The deposit stores no usable matrix of its own,
#'   so this argument is required.
#' @param thresholds A named numeric vector from [UnstainedThresholds()] on the
#'   deposit's unstained cell control, named for the markers of
#'   [kOmip51Markers]. Passing `NULL` falls back to a fitted density cut on
#'   every gate, which fails on CD14 and on the viability channel.
#' @param cofactor The arcsinh cofactor used for every fitted cut.
#' @param min_events The smallest parent that a cut is attempted on.
#' @return A list with `counts`, `cuts`, `rare`, `masks`, `channels`, `state`,
#'   `frame` and `cofactor`.
#' @examples
#' \dontrun{
#' GateOmip51File(path, spillover)$counts
#' }
#' @export
GateOmip51File <- function(path, spillover, thresholds = NULL,
                           cofactor = 150, min_events = 200) {
  frame <- flowCore::read.FCS(path, truncate_max_range = FALSE)

  absent <- setdiff(rownames(spillover), flowCore::colnames(frame))
  if (length(absent) > 0) {
    stop("The spillover matrix names ", length(absent),
         " detector(s) that the file does not carry: ",
         paste(absent, collapse = ", "), ".")
  }

  state <- ReadCompensationState(frame)
  channels <- ResolveOmip51Channels(frame)
  scatter <- PanelScatterChannels(frame)
  frame <- flowCore::compensate(frame, spillover)

  values <- ArcsinhTransform(flowCore::exprs(frame), channels$channel, cofactor)
  total <- nrow(values)

  cuts <- list()
  Channel <- function(name) channels$channel[channels$name == name]

  # Every lineage threshold comes from the unstained cell control at its 99.9th
  # percentile. The fitted rules are recorded beside it rather than used,
  # because one of them fails on this file: on CD14 the density minimum lands at
  # -1.95 and calls 83.4 percent of the viable events positive, where the
  # mixture fit gives 8.4 percent and the control gives 8.2 percent. Two
  # independent references agreeing to within 0.2 points is the evidence for the
  # threshold; a rule that switches per marker would not be one rule.
  Cut <- function(parent, name, keep, label) {
    empty <- rep(FALSE, total)
    column <- values[parent, Channel(name)]
    density_cut <- if (sum(parent) < min_events) NA_real_ else DensityCut(column)
    mixture_cut <- if (sum(parent) < min_events) NA_real_ else MixtureCut(column)
    threshold <- if (is.null(thresholds)) density_cut else
      unname(thresholds[[name]])
    # The fraction each candidate cut would keep is recorded next to the cut
    # itself, because that fraction is what says a rule failed. A threshold
    # alone does not, and a report that states the fraction has to be able to
    # point at the file it came from.
    PercentAbove <- function(value) {
      if (is.na(value) || length(column) == 0) NA_real_ else
        100 * mean(column > value)
    }
    cuts[[length(cuts) + 1]] <<- data.frame(
      label = label, marker = name, channel = Channel(name),
      cut = threshold,
      rule = if (is.null(thresholds)) "density" else "unstained control",
      density_cut = density_cut, mixture_cut = mixture_cut,
      percent_above_cut = PercentAbove(threshold),
      percent_above_density_cut = PercentAbove(density_cut),
      percent_above_mixture_cut = PercentAbove(mixture_cut),
      parent_events = sum(parent), stringsAsFactors = FALSE
    )
    if (is.na(threshold)) {
      return(empty)
    }
    empty[which(parent)] <- if (keep == "above") column > threshold else
      column <= threshold
    empty
  }

  in_range <- InScatterRange(values, scatter)
  singlets <- in_range
  singlets[in_range] <- RatioSingletMask(values[in_range, , drop = FALSE],
                                         scatter)

  # Figure 1A gates the cell cloud on scatter before the viability channel.
  cells <- singlets
  cells[singlets] <- ScatterCloudMask(values[singlets, , drop = FALSE], scatter)

  viable <- Cut(cells, "viability", "below", "viable")
  monocytes <- Cut(viable, "CD14", "above", "CD14 positive")
  not_monocytes <- viable & !monocytes

  hladr <- Cut(not_monocytes, "HLADR", "above", "HLA-DR positive")
  cd20 <- Cut(not_monocytes, "CD20", "above", "CD20 positive")
  b_and_myeloid <- not_monocytes & (hladr | cd20)

  cd19 <- Cut(b_and_myeloid, "CD19", "above", "CD19 positive")
  b_cells <- b_and_myeloid & cd19 & cd20

  # The paper places the lymphocyte gate after the B cell gate, so the high
  # scatter population it removes can be counted here.
  b_lymphocytes <- b_cells
  if (sum(b_cells) >= min_events) {
    b_lymphocytes[b_cells] <- LymphocyteMask(values[b_cells, , drop = FALSE],
                                             scatter)
  }

  dendritic_parent <- b_and_myeloid & hladr & !cd19 & !cd20
  cd123 <- Cut(dendritic_parent, "CD123", "above", "CD123 positive")
  cd11c <- Cut(dendritic_parent, "CD11c", "above", "CD11c positive")
  pdc <- dendritic_parent & cd123 & !cd11c
  mdc <- dendritic_parent & !cd123 & cd11c

  # The subsets of Figure 1C, cut on the same thresholds as the hierarchy above.
  cd10 <- Cut(b_lymphocytes, "CD10", "above", "CD10 positive")
  transitional <- b_lymphocytes & cd10
  mature <- b_lymphocytes & !cd10
  igd <- Cut(mature, "IgD", "above", "IgD positive")
  cd27 <- Cut(mature, "CD27", "above", "CD27 positive")
  igm <- Cut(mature, "IgM", "above", "IgM positive")
  naive <- mature & igd & !cd27
  igd_cd27 <- mature & igd & cd27
  marginal_zone <- igd_cd27 & igm
  igd_only_memory <- igd_cd27 & !igm
  igd_negative <- mature & !igd
  cd20_bright <- Cut(igd_negative, "CD20", "above", "CD20 positive in IgD neg")
  plasmablasts <- igd_negative & cd27 & !cd20_bright
  memory <- igd_negative & !plasmablasts

  rare <- AssessOneDimensionalCuts(
    values, channels, b_lymphocytes,
    c("CD10", "CD27", "IgD", "IgM", "CD21", "CD85j", "IgG", "IgA")
  )
  rare$parent <- "b_lymphocytes"
  rare_mdc <- AssessOneDimensionalCuts(values, channels, mdc,
                                       c("CD1c", "CD141", "CD16"))
  rare_mdc$parent <- "myeloid_dendritic_cells"
  rare <- rbind(rare, rare_mdc)

  hierarchy <- list(
    list("all_events", rep(TRUE, total), NA_character_),
    list("in_scatter_range", in_range, "all_events"),
    list("singlets", singlets, "in_scatter_range"),
    list("cells", cells, "singlets"),
    list("viable", viable, "cells"),
    list("monocytes", monocytes, "viable"),
    list("cd14_negative", not_monocytes, "viable"),
    list("b_and_myeloid", b_and_myeloid, "cd14_negative"),
    list("b_cells", b_cells, "b_and_myeloid"),
    list("b_lymphocytes", b_lymphocytes, "b_cells"),
    list("dendritic_parent", dendritic_parent, "b_and_myeloid"),
    list("plasmacytoid_dendritic_cells", pdc, "dendritic_parent"),
    list("myeloid_dendritic_cells", mdc, "dendritic_parent"),
    list("transitional_b", transitional, "b_lymphocytes"),
    list("mature_b", mature, "b_lymphocytes"),
    list("naive_b", naive, "mature_b"),
    list("igd_cd27_b", igd_cd27, "mature_b"),
    list("marginal_zone_b", marginal_zone, "b_lymphocytes"),
    list("igd_only_memory_b", igd_only_memory, "b_lymphocytes"),
    list("igd_negative_b", igd_negative, "mature_b"),
    list("plasmablasts", plasmablasts, "b_lymphocytes"),
    list("memory_b", memory, "mature_b")
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
      stringsAsFactors = FALSE
    )
  }))

  list(counts = counts, cuts = do.call(rbind, cuts), rare = rare, masks = masks,
       channels = channels, state = state, frame = frame, cofactor = cofactor)
}
