# The COVID-19 intracellular cytokine deposit of Vanderbeke 2021.
#
# FR-FCM-Z2KP holds 49 files from a 34 colour intracellular cytokine panel on a
# BD Symphony, and it ships the FlowJo workspace the authors used. That
# workspace
# is the first one in this repository that gates every sample rather than one,
# so
# the manual arm covers the whole cohort instead of a single donor.
#
# Three properties of the deposit shape the code.
#
# The files are exported after the live cell gate, so the analysis starts inside
# a population that somebody already drew. The scatter, viability and debris
# steps of the earlier reports have nothing to act on.
#
# The values are compensated. Every channel is named `FJComp-`, which is what
# FlowJo writes when it exports compensated data, so the pipeline applies no
# matrix.
#
# The manual hierarchy is shallow and wide. CD3 carries CD4 and CD8 beneath it,
# while CD14 monocytes and CD19 sit beside CD3 rather than under it. The
# automated route reproduces that shape rather than a deeper one, so the two
# arms
# count the same denominators.

#' Read the manual gates of the deposited FlowJo workspace
#'
#' `CytoML` applies the workspace to the deposited files and returns the event
#' count of every population. Nothing here refits a gate.
#'
#' @param workspace_path Path to the `.wsp` file.
#' @param fcs_path The folder that holds the FCS files.
#' @param group The sample group index inside the workspace.
#' @return A `data.frame` with one row per sample and population, holding
#'   `sample`, `population` and `count`.
#' @examples
#' \dontrun{
#' ReadCovidManualGates(workspace, "data/datasets/flowrepository/FR-FCM-Z2KP")
#' }
#' @export
ReadCovidManualGates <- function(workspace_path, fcs_path, group = 1) {
  if (!file.exists(workspace_path)) {
    stop("The workspace does not exist: ", workspace_path)
  }
  workspace <- CytoML::open_flowjo_xml(workspace_path)
  gating_set <- CytoML::flowjo_to_gatingset(
    workspace, name = group, path = fcs_path, execute = TRUE
  )
  counts <- flowWorkspace::gs_pop_get_stats(gating_set, type = "count")
  counts <- as.data.frame(counts)
  names(counts)[names(counts) == "sample"] <- "sample"
  data.frame(
    sample = as.character(counts$sample),
    population = as.character(counts$pop),
    count = as.numeric(counts$count),
    stringsAsFactors = FALSE
  )
}

#' Read the clinical group out of a file name
#'
#' The deposit encodes the group in the file name. `HC` is a healthy control,
#' `W` is a ward patient and `ICU` is an intensive care patient.
#'
#' @param file_names A character vector of file names.
#' @return A `data.frame` with the columns `file_name`, `group` and `patient`.
#' @examples
#' ParseCovidFileNames("export_COVID19 samples 21_04_20_ST3_COVID19_ICU_005_A ST3 210420_080_Live_cells.fcs")
#' @export
ParseCovidFileNames <- function(file_names) {
  group <- rep(NA_character_, length(file_names))
  group[grepl("_HC_", file_names, fixed = TRUE)] <- "healthy"
  group[grepl("_W_", file_names, fixed = TRUE)] <- "ward"
  group[grepl("_ICU_", file_names, fixed = TRUE)] <- "intensive care"

  patient <- sub(".*_(HC|W|ICU)_([0-9]+).*", "\\1\\2", file_names)
  patient[is.na(group)] <- NA_character_

  data.frame(
    file_name = file_names,
    group = group,
    patient = patient,
    stringsAsFactors = FALSE
  )
}

#' Gate one file of the COVID-19 deposit
#'
#' The hierarchy matches the deposited workspace. CD3 sits at the top, CD4 and
#' CD8 sit inside it, and CD14 and CD19 sit beside it. The naive and memory
#' subsets and the cytokine gates are added below, because the panel carries
#' them
#' and the paper's claims need them.
#'
#' @param path Path to the FCS file.
#' @param minimum_events The smallest population that a subset frequency is
#'   reported for. A frequency counted in fewer events than this is returned as
#'   `NA` rather than as a number.
#' @return A list with `counts`, a one row `data.frame`, `cuts`, one row per
#'   marker, and `error`.
#' @examples
#' \dontrun{
#' GateCovidFile(path)$counts
#' }
#' @export
GateCovidFile <- function(path, minimum_events = 100) {
  Fail <- function(message) {
    list(counts = NULL, cuts = NULL, error = message)
  }

  frame <- try(
    flowCore::read.FCS(path, truncate_max_range = FALSE, emptyValue = FALSE),
    silent = TRUE
  )
  if (methods::is(frame, "try-error")) {
    return(Fail(trimws(as.character(frame))))
  }

  parameters <- flowCore::pData(flowCore::parameters(frame))
  Channel <- function(marker) {
    hit <- which(toupper(trimws(parameters$desc)) == toupper(marker))
    if (length(hit) == 0) NA_character_ else parameters$name[hit[1]]
  }

  needed <- c("CD3", "CD4", "CD8", "CD14", "CD19", "CD45RA", "CCR7", "PD-1",
              "IFNg", "livedead")
  channels <- vapply(needed, Channel, character(1))
  if (anyNA(channels)) {
    return(Fail(paste0("No channel carries: ",
                       paste(needed[is.na(channels)], collapse = ", "))))
  }

  transform_list <- NULL
  for (width in c(4.5, 5.5, 6.5, 8, 10)) {
    attempt <- try(
      suppressWarnings(flowCore::estimateLogicle(
        frame, channels = unname(channels), m = width
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

  events <- flowCore::exprs(flowCore::transform(frame, transform_list))
  total <- nrow(events)
  if (total < 500) {
    return(Fail(paste0("The file holds only ", total, " events")))
  }

  cuts <- list()
  Cut <- function(marker, values, parent) {
    result <- ResolveCut(values)
    cuts[[length(cuts) + 1]] <<- data.frame(
      file_name = basename(path), marker = marker, parent = parent,
      cut = result$cut, rule = result$rule, stringsAsFactors = FALSE
    )
    if (result$rule == "none") NULL else result$cut
  }

  Percent <- function(mask, parent_size) {
    if (parent_size < minimum_events) NA_real_ else 100 * sum(mask) / parent_size
  }

  cd3_cut <- Cut("CD3", events[, channels[["CD3"]]], "root")
  if (is.null(cd3_cut)) {
    return(Fail("No cut could be fitted on CD3"))
  }
  t_cells <- events[events[, channels[["CD3"]]] > cd3_cut, , drop = FALSE]

  cd14_cut <- Cut("CD14", events[, channels[["CD14"]]], "root")
  cd19_cut <- Cut("CD19", events[, channels[["CD19"]]], "root")

  counts <- data.frame(
    file_name = basename(path),
    total_events = total,
    cd3_events = nrow(t_cells),
    CD3_percent = 100 * nrow(t_cells) / total,
    CD14_percent = if (is.null(cd14_cut)) {
      NA_real_
    } else {
      100 * mean(events[, channels[["CD14"]]] > cd14_cut)
    },
    CD19_percent = if (is.null(cd19_cut)) {
      NA_real_
    } else {
      100 * mean(events[, channels[["CD19"]]] > cd19_cut)
    },
    stringsAsFactors = FALSE
  )

  if (nrow(t_cells) < minimum_events) {
    counts$CD4_percent <- NA_real_
    counts$CD8_percent <- NA_real_
    counts$cd4_events <- NA_integer_
    counts$cd8_events <- NA_integer_
    return(list(counts = counts, cuts = do.call(rbind, cuts),
                error = NA_character_))
  }

  cd4_cut <- Cut("CD4", t_cells[, channels[["CD4"]]], "CD3")
  cd8_cut <- Cut("CD8", t_cells[, channels[["CD8"]]], "CD3")
  if (is.null(cd4_cut) || is.null(cd8_cut)) {
    return(Fail("No cut could be fitted on CD4 or on CD8"))
  }
  is_cd4 <- t_cells[, channels[["CD4"]]] > cd4_cut
  is_cd8 <- t_cells[, channels[["CD8"]]] > cd8_cut
  cd4 <- t_cells[is_cd4 & !is_cd8, , drop = FALSE]
  cd8 <- t_cells[is_cd8 & !is_cd4, , drop = FALSE]

  counts$cd4_events <- nrow(cd4)
  counts$cd8_events <- nrow(cd8)
  counts$CD4_percent <- 100 * nrow(cd4) / nrow(t_cells)
  counts$CD8_percent <- 100 * nrow(cd8) / nrow(t_cells)

  # The naive and memory quadrant, and the two functional markers the paper's
  # claims need. Both cuts are fitted inside the CD3 gate and reused below it,
  # which is the rule the FR-FCM-Z282 report established.
  ra_cut <- Cut("CD45RA", t_cells[, channels[["CD45RA"]]], "CD3")
  ccr7_cut <- Cut("CCR7", t_cells[, channels[["CCR7"]]], "CD3")
  pd1_cut <- Cut("PD-1", t_cells[, channels[["PD-1"]]], "CD3")
  ifng_cut <- Cut("IFNg", t_cells[, channels[["IFNg"]]], "CD3")

  if (!is.null(ifng_cut)) {
    counts$Th1_percent_of_CD4 <- Percent(
      cd4[, channels[["IFNg"]]] > ifng_cut, nrow(cd4)
    )
    counts$Tc1_percent_of_CD8 <- Percent(
      cd8[, channels[["IFNg"]]] > ifng_cut, nrow(cd8)
    )
  } else {
    counts$Th1_percent_of_CD4 <- NA_real_
    counts$Tc1_percent_of_CD8 <- NA_real_
  }

  if (!is.null(ra_cut) && !is.null(ccr7_cut) && !is.null(pd1_cut)) {
    for (lineage in c("CD4", "CD8")) {
      subset <- if (lineage == "CD4") cd4 else cd8
      is_ra <- subset[, channels[["CD45RA"]]] > ra_cut
      is_ccr7 <- subset[, channels[["CCR7"]]] > ccr7_cut
      quadrants <- list(
        N = is_ra & is_ccr7, CM = !is_ra & is_ccr7,
        EM = !is_ra & !is_ccr7, TEMRA = is_ra & !is_ccr7
      )
      for (name in names(quadrants)) {
        piece <- subset[quadrants[[name]], , drop = FALSE]
        counts[[paste0(lineage, "_", name, "_percent")]] <-
          Percent(quadrants[[name]], nrow(subset))
        counts[[paste0(lineage, "_", name, "_PD1")]] <- Percent(
          piece[, channels[["PD-1"]]] > pd1_cut, nrow(piece)
        )
      }
    }
  }

  list(counts = counts, cuts = do.call(rbind, cuts), error = NA_character_)
}

#' Compare a measure across the three clinical groups
#'
#' @param values A numeric vector.
#' @param group A vector that names the group of each value.
#' @return A one row `data.frame` with the median of each group, the
#'   Kruskal-Wallis p across the three groups, and the Mann-Whitney p for every
#'   COVID-19 sample against the healthy controls.
#' @examples
#' \dontrun{
#' CompareByGroup(measures$CD4_percent, measures$group)
#' }
#' @export
CompareByGroup <- function(values, group) {
  usable <- is.finite(values) & !is.na(group)
  values <- values[usable]
  group <- group[usable]
  if (length(unique(group)) < 2) {
    return(NULL)
  }

  Median <- function(name) {
    piece <- values[group == name]
    if (length(piece) == 0) NA_real_ else stats::median(piece)
  }

  covid <- values[group != "healthy"]
  healthy <- values[group == "healthy"]
  severity <- values[group != "healthy"]
  severity_group <- group[group != "healthy"]

  covid_test <- if (length(healthy) >= 3 && length(covid) >= 3) {
    suppressWarnings(stats::wilcox.test(covid, healthy))$p.value
  } else {
    NA_real_
  }
  severity_test <- if (length(unique(severity_group)) >= 2) {
    suppressWarnings(
      stats::kruskal.test(severity, as.factor(severity_group))
    )$p.value
  } else {
    NA_real_
  }

  data.frame(
    healthy_median = Median("healthy"),
    ward_median = Median("ward"),
    intensive_care_median = Median("intensive care"),
    covid_against_healthy_p = covid_test,
    between_severity_p = severity_test,
    stringsAsFactors = FALSE
  )
}
