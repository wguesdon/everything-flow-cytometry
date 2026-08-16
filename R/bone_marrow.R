# The bone marrow deposits of Oetjen 2018.
#
# Two accessions describe the same 20 donors. FR-FCM-ZYQ9 holds 13 colour flow
# cytometry on five panels, and FR-FCM-ZYQB holds a 49 marker Helios mass
# cytometry panel on eight of those donors.
#
# Four properties of the flow deposit shape the code.
#
# Every file came off one LSRFortessa, so the detector names are the same in all
# 132 files. A detector therefore identifies a marker once the panel is known,
# which is what `gating/oetjen2018_panels.csv` records. That matters because 46
# of the files carry no marker name at all.
#
# The T cell panel carries no CD45. Four panels put CD45 in the BV785 detector
# and the T panel puts CD3 there instead, so a CD45 percentage cannot come from
# the T panel.
#
# The deposit ships the depositors' own scatter gate and viability gate, in a
# cytoflow workflow file. Both are read from `gating/oetjen2018_manual_gates.csv`
# rather than refitted, so the entry into the analysis is the authors' own.
#
# The viability channel runs the opposite way to the usual convention. The
# depositors' gate is named `dead_lymphocytes` and it covers the LOW end of the
# V545 detector, so a live cell is bright and not dim. A pipeline that assumes
# the usual direction keeps the debris and discards the cells, and it does so
# without an error. `LiveMask()` reads the direction from the gate file.

#' Read the panel table of the Oetjen 2018 flow deposit
#'
#' @param path Path to `gating/oetjen2018_panels.csv`.
#' @return A `data.frame` with the columns `panel`, `detector`, `fluorochrome`
#'   and `marker`.
#' @examples
#' \dontrun{
#' ReadBoneMarrowPanels("gating/oetjen2018_panels.csv")
#' }
#' @export
ReadBoneMarrowPanels <- function(path) {
  if (!file.exists(path)) {
    stop("The panel table does not exist: ", path)
  }
  panels <- utils::read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
  required <- c("panel", "detector", "fluorochrome", "marker")
  missing <- setdiff(required, names(panels))
  if (length(missing) > 0) {
    stop("The panel table lacks these columns: ", paste(missing, collapse = ", "))
  }
  panels
}

#' Read the sample sheet of the Oetjen 2018 cohort
#'
#' The sheet is derived from Table 1 of the paper and it is committed, so the age
#' and the sex of each donor travel with the repository.
#'
#' @param path Path to `gating/oetjen2018_donor_metadata.csv`.
#' @return A `data.frame` with one row per sample.
#' @examples
#' \dontrun{
#' ReadBoneMarrowDonors("gating/oetjen2018_donor_metadata.csv")
#' }
#' @export
ReadBoneMarrowDonors <- function(path) {
  if (!file.exists(path)) {
    stop("The donor table does not exist: ", path)
  }
  donors <- utils::read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
  required <- c("sample", "donor", "sex", "age")
  missing <- setdiff(required, names(donors))
  if (length(missing) > 0) {
    stop("The donor table lacks these columns: ", paste(missing, collapse = ", "))
  }
  donors$age <- as.numeric(donors$age)
  donors
}

#' Read the gates that the depositors shipped with the experiment
#'
#' FR-FCM-ZYQ9 carries a cytoflow workflow file. It holds a scatter polygon named
#' `lymphocytes` and a second polygon named `dead_lymphocytes`. The vertices were
#' extracted once into a CSV, so the analysis does not parse a workflow format.
#'
#' @param path Path to `gating/oetjen2018_manual_gates.csv`.
#' @return A `data.frame` with one row per vertex.
#' @examples
#' \dontrun{
#' ReadBoneMarrowGates("gating/oetjen2018_manual_gates.csv")
#' }
#' @export
ReadBoneMarrowGates <- function(path) {
  if (!file.exists(path)) {
    stop("The manual gate table does not exist: ", path)
  }
  gates <- utils::read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
  required <- c("gate", "x_channel", "y_channel", "vertex", "x", "y")
  missing <- setdiff(required, names(gates))
  if (length(missing) > 0) {
    stop("The gate table lacks these columns: ", paste(missing, collapse = ", "))
  }
  gates
}

#' Test which points fall inside a polygon
#'
#' The ray casting rule. A point is inside when a ray drawn from it crosses the
#' boundary an odd number of times.
#'
#' @param x A numeric vector of x coordinates.
#' @param y A numeric vector of y coordinates.
#' @param vertex_x The x coordinates of the polygon, in order.
#' @param vertex_y The y coordinates of the polygon, in order.
#' @return A logical vector, `TRUE` for a point inside the polygon.
#' @examples
#' InPolygon(c(0.5, 5), c(0.5, 5), c(0, 1, 1, 0), c(0, 0, 1, 1))
#' @export
InPolygon <- function(x, y, vertex_x, vertex_y) {
  if (length(vertex_x) < 3) {
    stop("A polygon needs at least three vertices")
  }
  inside <- rep(FALSE, length(x))
  n <- length(vertex_x)
  previous <- n
  for (index in seq_len(n)) {
    crosses <- (vertex_y[index] > y) != (vertex_y[previous] > y)
    boundary <- (vertex_x[previous] - vertex_x[index]) *
      (y - vertex_y[index]) /
      (vertex_y[previous] - vertex_y[index]) + vertex_x[index]
    inside <- xor(inside, crosses & (x < boundary))
    previous <- index
  }
  inside
}

#' Mark the live cells of one file
#'
#' The direction comes from the deposited gate and not from a convention. The
#' polygon is named `dead_lymphocytes`, so an event is live when its viability
#' value sits above the top of that polygon.
#'
#' @param values The raw viability values.
#' @param gates The table returned by [ReadBoneMarrowGates()].
#' @return A logical vector, `TRUE` for a live event.
#' @examples
#' \dontrun{
#' LiveMask(flowCore::exprs(frame)[, "V545-A"], gates)
#' }
#' @export
LiveMask <- function(values, gates) {
  dead <- gates[gates$gate == "dead_lymphocytes", ]
  if (nrow(dead) == 0) {
    stop("The gate table names no dead cell gate")
  }
  values > max(dead$y)
}

#' Parse the panel and the sample out of a file name
#'
#' A stained file is named like `2-13-17 B cell Panel_B_E_C05_004.fcs`, where the
#' letters after the panel word are the panel code and the next field is the
#' sample. An unstained control is named `Unstained_Unstained_E_A05_005.fcs`. One
#' monocyte file carries commas in the sample field.
#'
#' @param file_names A character vector of file names.
#' @return A `data.frame` with the columns `file_name`, `panel` and `sample`.
#'   Both are `NA` when the name matches neither pattern.
#' @examples
#' ParseBoneMarrowFileNames("2-13-17 B cell Panel_B_E_C05_004.fcs")
#' @export
ParseBoneMarrowFileNames <- function(file_names) {
  stained <- "^(?:[0-9-]+ )?.+?[ _][Pp]anel_([A-Za-z]+)_([A-Za-z]+)[,0-9a-z]*_"
  unstained <- "^Unstained_Unstained_([A-Za-z]+)[,0-9a-z]*_"

  panel <- rep(NA_character_, length(file_names))
  sample <- rep(NA_character_, length(file_names))

  hit <- regmatches(file_names, regexec(stained, file_names))
  for (index in seq_along(hit)) {
    if (length(hit[[index]]) == 3) {
      panel[index] <- hit[[index]][2]
      sample[index] <- hit[[index]][3]
    }
  }

  hit <- regmatches(file_names, regexec(unstained, file_names))
  for (index in seq_along(hit)) {
    if (length(hit[[index]]) == 2) {
      panel[index] <- "Unstained"
      sample[index] <- hit[[index]][2]
    }
  }

  data.frame(
    file_name = file_names,
    panel = panel,
    sample = sample,
    stringsAsFactors = FALSE
  )
}

#' Match every marker of one panel to a channel of one file
#'
#' The first pass reads the marker name that the operator typed. The second pass
#' takes the detector name from the published panel design, which is what rescues
#' the 46 files whose marker fields are empty.
#'
#' @param frame A `flowFrame`.
#' @param panels The table returned by [ReadBoneMarrowPanels()].
#' @param panel The panel code, one of `"T"`, `"B"`, `"NK"`, `"Mono"` or `"DC"`.
#' @return A `data.frame` with the columns `marker`, `channel` and `resolved_by`.
#' @examples
#' \dontrun{
#' ResolveBoneMarrowChannels(frame, panels, "T")
#' }
#' @export
ResolveBoneMarrowChannels <- function(frame, panels, panel) {
  wanted <- panels[panels$panel == panel, ]
  if (nrow(wanted) == 0) {
    stop("The panel table names no marker for the panel: ", panel)
  }

  parameters <- flowCore::pData(flowCore::parameters(frame))
  result <- data.frame(
    marker = wanted$marker,
    channel = NA_character_,
    resolved_by = NA_character_,
    stringsAsFactors = FALSE
  )
  Clean <- function(x) toupper(trimws(x))

  for (row in seq_len(nrow(result))) {
    typed <- which(Clean(parameters$desc) == Clean(result$marker[row]))
    if (length(typed) > 0) {
      result$channel[row] <- parameters$name[typed[1]]
      result$resolved_by[row] <- "marker"
      next
    }
    detector <- which(Clean(parameters$name) == Clean(wanted$detector[row]))
    if (length(detector) > 0) {
      result$channel[row] <- parameters$name[detector[1]]
      result$resolved_by[row] <- "detector"
    }
  }
  result
}

#' Gate one flow cytometry file of the bone marrow deposit
#'
#' The entry is the depositors' own scatter gate and their own viability gate.
#' Everything below it is fitted by [ResolveCut()], which is the automated rule
#' this repository uses everywhere else.
#'
#' Every marker of the panel gets a cut and a percent positive, whether or not
#' the panel design intends it as a lineage marker. That is deliberate. The
#' spread of one marker across 22 samples of the same tissue is the measurement
#' that says whether an automated cut can be trusted on this deposit.
#'
#' @param path Path to the FCS file.
#' @param panels The table returned by [ReadBoneMarrowPanels()].
#' @param panel The panel code.
#' @param gates The table returned by [ReadBoneMarrowGates()].
#' @return A list with `counts`, a one row `data.frame` of the event counts,
#'   `markers`, one row per marker with its cut and its percent positive,
#'   `channels`, the output of [ResolveBoneMarrowChannels()], and `error`.
#' @examples
#' \dontrun{
#' GateBoneMarrowFile(path, panels, "T", gates)$counts
#' }
#' @export
GateBoneMarrowFile <- function(path, panels, panel, gates) {
  Fail <- function(message) {
    list(counts = NULL, markers = NULL, channels = NULL, error = message)
  }

  frame <- try(
    flowCore::read.FCS(path, truncate_max_range = FALSE, emptyValue = FALSE),
    silent = TRUE
  )
  if (methods::is(frame, "try-error")) {
    return(Fail(trimws(as.character(frame))))
  }

  channels <- try(ResolveBoneMarrowChannels(frame, panels, panel), silent = TRUE)
  if (methods::is(channels, "try-error")) {
    return(Fail(trimws(as.character(channels))))
  }
  if (anyNA(channels$channel)) {
    return(Fail(paste0(
      "No channel resolves for: ",
      paste(channels$marker[is.na(channels$channel)], collapse = ", ")
    )))
  }

  scatter <- try(ScatterChannels(frame), silent = TRUE)
  if (methods::is(scatter, "try-error")) {
    return(Fail(trimws(as.character(scatter))))
  }

  transform_list <- NULL
  for (width in c(4.5, 5.5, 6.5, 8, 10)) {
    attempt <- try(
      suppressWarnings(flowCore::estimateLogicle(
        frame, channels = channels$channel, m = width
      )),
      silent = TRUE
    )
    if (!methods::is(attempt, "try-error")) {
      transform_list <- attempt
      break
    }
  }
  if (is.null(transform_list)) {
    return(Fail("The logicle transform failed at every decade count"))
  }

  raw <- flowCore::exprs(frame)
  events <- flowCore::exprs(flowCore::transform(frame, transform_list))
  total <- nrow(raw)

  Channel <- function(marker) {
    channels$channel[channels$marker == marker]
  }

  singlets <- SingletMask(raw, area = scatter[["forward_area"]],
                          height = scatter[["forward_height"]])

  scatter_gate <- gates[gates$gate == "lymphocytes", ]
  scatter_gate <- scatter_gate[order(scatter_gate$vertex), ]
  in_scatter <- InPolygon(
    raw[, scatter_gate$x_channel[1]], raw[, scatter_gate$y_channel[1]],
    scatter_gate$x, scatter_gate$y
  )

  live <- LiveMask(raw[, Channel("LIVE/DEAD")], gates)
  keep <- singlets & in_scatter & live
  if (sum(keep) < 500) {
    return(Fail("Fewer than 500 events survive the deposited gates"))
  }
  entry <- events[keep, , drop = FALSE]

  # The paper counts the CD45 positive events it collected, before any scatter
  # gate. The deposited polygon is a lymphocyte gate, so a CD45 count taken after
  # it answers a different question. Both are recorded.
  wide <- events[singlets & live, , drop = FALSE]

  # Four panels carry CD45. The T panel does not, so its denominator is the
  # population the deposited gates produced.
  wide_cd45 <- NA_integer_
  if ("CD45" %in% channels$marker) {
    cd45 <- ResolveCut(entry[, Channel("CD45")])
    if (cd45$rule == "none") {
      return(Fail("No cut could be fitted on CD45"))
    }
    parent <- entry[entry[, Channel("CD45")] > cd45$cut, , drop = FALSE]
    cd45_cut <- cd45$cut
    cd45_rule <- cd45$rule
    outside <- ResolveCut(wide[, Channel("CD45")])
    if (outside$rule != "none") {
      wide_cd45 <- sum(wide[, Channel("CD45")] > outside$cut)
    }
  } else {
    parent <- entry
    cd45_cut <- NA_real_
    cd45_rule <- "no CD45 in this panel"
  }
  if (nrow(parent) < 500) {
    return(Fail("Fewer than 500 events in the parent gate"))
  }

  markers <- do.call(rbind, lapply(channels$marker, function(marker) {
    if (marker == "LIVE/DEAD") {
      return(NULL)
    }
    values <- parent[, Channel(marker)]
    result <- ResolveCut(values)
    data.frame(
      file_name = basename(path),
      panel = panel,
      marker = marker,
      cut = result$cut,
      rule = result$rule,
      percent_positive = if (result$rule == "none") {
        NA_real_
      } else {
        100 * mean(values > result$cut)
      },
      stringsAsFactors = FALSE
    )
  }))

  counts <- data.frame(
    file_name = basename(path),
    panel = panel,
    total_events = total,
    singlet_events = sum(singlets),
    scatter_events = sum(singlets & in_scatter),
    live_events = sum(keep),
    parent_events = nrow(parent),
    cd45_before_scatter = wide_cd45,
    cd45_cut = cd45_cut,
    cd45_rule = cd45_rule,
    stringsAsFactors = FALSE
  )

  list(
    counts = counts,
    markers = markers,
    channels = cbind(file_name = basename(path), channels),
    error = NA_character_
  )
}

#' Gate the T cell subsets of one file
#'
#' The T panel carries CD3, CD4, CD8, CD45RA and CCR7, so it holds the naive and
#' memory hierarchy that the mass cytometry panel also holds. This function
#' reports that hierarchy alone, which is what the flow against mass comparison
#' needs.
#'
#' @param path Path to the FCS file.
#' @param panels The table returned by [ReadBoneMarrowPanels()].
#' @param gates The table returned by [ReadBoneMarrowGates()].
#' @return A list with `counts` and `error`.
#' @examples
#' \dontrun{
#' GateTcellSubsets(path, panels, gates)$counts
#' }
#' @export
GateTcellSubsets <- function(path, panels, gates) {
  base <- GateBoneMarrowFile(path, panels, "T", gates)
  if (!is.na(base$error)) {
    return(list(counts = NULL, error = base$error))
  }

  frame <- flowCore::read.FCS(path, truncate_max_range = FALSE,
                              emptyValue = FALSE)
  channels <- ResolveBoneMarrowChannels(frame, panels, "T")
  scatter <- ScatterChannels(frame)
  transform_list <- suppressWarnings(
    flowCore::estimateLogicle(frame, channels = channels$channel)
  )
  raw <- flowCore::exprs(frame)
  events <- flowCore::exprs(flowCore::transform(frame, transform_list))
  Channel <- function(marker) channels$channel[channels$marker == marker]

  scatter_gate <- gates[gates$gate == "lymphocytes", ]
  scatter_gate <- scatter_gate[order(scatter_gate$vertex), ]
  keep <- SingletMask(raw, area = scatter[["forward_area"]],
                      height = scatter[["forward_height"]]) &
    InPolygon(raw[, scatter_gate$x_channel[1]], raw[, scatter_gate$y_channel[1]],
              scatter_gate$x, scatter_gate$y) &
    LiveMask(raw[, Channel("LIVE/DEAD")], gates)
  entry <- events[keep, , drop = FALSE]

  dump <- ResolveCut(entry[, Channel("CD14/19")])
  if (dump$rule == "none") {
    return(list(counts = NULL, error = "No cut on the CD14 and CD19 dump"))
  }
  undumped <- entry[entry[, Channel("CD14/19")] < dump$cut, , drop = FALSE]

  cd3 <- ResolveCut(undumped[, Channel("CD3")])
  if (cd3$rule == "none") {
    return(list(counts = NULL, error = "No cut on CD3"))
  }
  t_cells <- undumped[undumped[, Channel("CD3")] > cd3$cut, , drop = FALSE]
  if (nrow(t_cells) < 500) {
    return(list(counts = NULL, error = "Fewer than 500 T cells"))
  }

  cd4 <- ResolveCut(t_cells[, Channel("CD4")])
  cd8 <- ResolveCut(t_cells[, Channel("CD8")])
  ra <- ResolveCut(t_cells[, Channel("CD45RA")])
  ccr7 <- ResolveCut(t_cells[, Channel("CCR7")])
  if (any(c(cd4$rule, cd8$rule, ra$rule, ccr7$rule) == "none")) {
    return(list(counts = NULL, error = "No cut on CD4, CD8, CD45RA or CCR7"))
  }

  is_cd4 <- t_cells[, Channel("CD4")] > cd4$cut
  is_cd8 <- t_cells[, Channel("CD8")] > cd8$cut
  subsets <- list(
    CD4 = t_cells[is_cd4 & !is_cd8, , drop = FALSE],
    CD8 = t_cells[is_cd8 & !is_cd4, , drop = FALSE]
  )

  frequencies <- list()
  for (name in names(subsets)) {
    quadrant <- QuadrantFrequencies.(
      subsets[[name]], Channel("CD45RA"), Channel("CCR7"), ra$cut, ccr7$cut
    )
    names(quadrant) <- paste0(name, "_", names(quadrant))
    frequencies <- c(frequencies, as.list(quadrant))
  }

  counts <- data.frame(
    file_name = basename(path),
    technology = "flow",
    entry_events = nrow(entry),
    cd3_events = nrow(t_cells),
    cd4_events = nrow(subsets$CD4),
    cd8_events = nrow(subsets$CD8),
    CD3_percent_of_entry = 100 * nrow(t_cells) / nrow(entry),
    CD4_percent_of_CD3 = 100 * nrow(subsets$CD4) / nrow(t_cells),
    CD8_percent_of_CD3 = 100 * nrow(subsets$CD8) / nrow(t_cells),
    stringsAsFactors = FALSE
  )
  list(counts = cbind(counts, as.list(frequencies)), error = NA_character_)
}

# The four quadrants of a CD45RA against CCR7 plot, as percentages of the subset.
QuadrantFrequencies. <- function(events, ra_channel, ccr7_channel, ra_cut,
                                 ccr7_cut) {
  if (nrow(events) == 0) {
    return(stats::setNames(rep(NA_real_, 4), c("N", "CM", "EM", "TEMRA")))
  }
  is_ra <- events[, ra_channel] > ra_cut
  is_ccr7 <- events[, ccr7_channel] > ccr7_cut
  stats::setNames(
    100 * c(
      sum(is_ra & is_ccr7), sum(!is_ra & is_ccr7),
      sum(!is_ra & !is_ccr7), sum(is_ra & !is_ccr7)
    ) / nrow(events),
    c("N", "CM", "EM", "TEMRA")
  )
}

#' Gate one mass cytometry file of the bone marrow deposit
#'
#' A Helios file has no scatter and no viability dye of the usual kind. Intact
#' single cells carry two iridium DNA channels and a dead cell takes up the
#' platinum stain, so those channels replace the scatter gate and the viability
#' gate. The arcsinh transform with a cofactor of 5 is the convention for mass
#' cytometry.
#'
#' @param path Path to the FCS file.
#' @param cofactor The arcsinh cofactor. Defaults to 5.
#' @return A list with `counts` and `error`.
#' @examples
#' \dontrun{
#' GateMassCytometryFile(path)$counts
#' }
#' @export
GateMassCytometryFile <- function(path, cofactor = 5) {
  Fail <- function(message) list(counts = NULL, error = message)

  frame <- try(
    flowCore::read.FCS(path, truncate_max_range = FALSE, emptyValue = FALSE),
    silent = TRUE
  )
  if (methods::is(frame, "try-error")) {
    return(Fail(trimws(as.character(frame))))
  }

  parameters <- flowCore::pData(flowCore::parameters(frame))
  Find <- function(pattern) {
    hit <- grep(pattern, parameters$desc, ignore.case = TRUE)
    if (length(hit) == 0) NA_character_ else parameters$name[hit[1]]
  }
  needed <- c(
    dna1 = "191Ir", viability = "195Pt", CD45 = "_CD45$", CD3 = "_CD3$",
    CD4 = "_CD4$", CD8 = "_CD8a$", CD45RA = "_CD45RA$", CCR7 = "_CD197$"
  )
  channels <- vapply(needed, Find, character(1))
  if (anyNA(channels)) {
    return(Fail(paste0("No channel matches: ",
                       paste(names(channels)[is.na(channels)],
                             collapse = ", "))))
  }

  events <- asinh(flowCore::exprs(frame)[, channels, drop = FALSE] / cofactor)
  colnames(events) <- names(channels)
  total <- nrow(events)

  Step <- function(values, above = TRUE) {
    result <- ResolveCut(values)
    if (result$rule == "none") {
      return(NULL)
    }
    if (above) values > result$cut else values < result$cut
  }

  intact_mask <- Step(events[, "dna1"])
  if (is.null(intact_mask)) {
    return(Fail("No cut on the DNA channel"))
  }
  intact <- events[intact_mask, , drop = FALSE]

  live_mask <- Step(intact[, "viability"], above = FALSE)
  if (is.null(live_mask)) {
    return(Fail("No cut on the viability channel"))
  }
  live <- intact[live_mask, , drop = FALSE]

  cd45_mask <- Step(live[, "CD45"])
  if (is.null(cd45_mask)) {
    return(Fail("No cut on CD45"))
  }
  leukocytes <- live[cd45_mask, , drop = FALSE]

  cd3_mask <- Step(leukocytes[, "CD3"])
  if (is.null(cd3_mask)) {
    return(Fail("No cut on CD3"))
  }
  t_cells <- leukocytes[cd3_mask, , drop = FALSE]
  if (nrow(t_cells) < 500) {
    return(Fail("Fewer than 500 T cells"))
  }

  cd4 <- ResolveCut(t_cells[, "CD4"])
  cd8 <- ResolveCut(t_cells[, "CD8"])
  ra <- ResolveCut(t_cells[, "CD45RA"])
  ccr7 <- ResolveCut(t_cells[, "CCR7"])
  if (any(c(cd4$rule, cd8$rule, ra$rule, ccr7$rule) == "none")) {
    return(Fail("No cut on CD4, CD8, CD45RA or CCR7"))
  }

  is_cd4 <- t_cells[, "CD4"] > cd4$cut
  is_cd8 <- t_cells[, "CD8"] > cd8$cut
  subsets <- list(
    CD4 = t_cells[is_cd4 & !is_cd8, , drop = FALSE],
    CD8 = t_cells[is_cd8 & !is_cd4, , drop = FALSE]
  )

  frequencies <- list()
  for (name in names(subsets)) {
    quadrant <- QuadrantFrequencies.(subsets[[name]], "CD45RA", "CCR7",
                                     ra$cut, ccr7$cut)
    names(quadrant) <- paste0(name, "_", names(quadrant))
    frequencies <- c(frequencies, as.list(quadrant))
  }

  counts <- data.frame(
    file_name = basename(path),
    technology = "mass",
    entry_events = nrow(leukocytes),
    cd3_events = nrow(t_cells),
    cd4_events = nrow(subsets$CD4),
    cd8_events = nrow(subsets$CD8),
    CD3_percent_of_entry = 100 * nrow(t_cells) / nrow(leukocytes),
    CD4_percent_of_CD3 = 100 * nrow(subsets$CD4) / nrow(t_cells),
    CD8_percent_of_CD3 = 100 * nrow(subsets$CD8) / nrow(t_cells),
    total_events = total,
    stringsAsFactors = FALSE
  )
  list(counts = cbind(counts, as.list(frequencies)), error = NA_character_)
}
