# Read FCS files and build the sample sheet that describes them.
#
# The example dataset is the FlowJo Basic Tutorial data. Its file names carry
# the
# experimental design, so the sample sheet is derived rather than typed. A typed
# sheet drifts from the files; a derived sheet cannot.

#' Parse the design out of a FlowJo tutorial file name
#'
#' The file names in the FlowJo Basic Tutorial dataset follow the pattern
#' `<donor>_<condition>_<well>_exp.fcs`, for example `LD1_NS+PI_C01_exp.fcs`.
#' The
#' condition holds two stimulation codes joined by a plus sign, where `NS` means
#' not stimulated and `PI` means stimulated with PMA and ionomycin.
#'
#' @param file_names A character vector of file names, with or without a path.
#' @return A `data.frame` with one row per file and the columns `file_name`,
#'   `donor`, `condition`, `stim_1`, `stim_2`, `well` and `stimulated`. The
#'   `stimulated` column is `TRUE` when either stimulation code is `PI`.
#' @examples
#' ParseTutorialFileNames(c("LD1_NS+NS_A01_exp.fcs", "LD2_PI+PI_D02_exp.fcs"))
#' @export
ParseTutorialFileNames <- function(file_names) {
  if (length(file_names) == 0) {
    stop("file_names is empty. Give at least one file name.")
  }

  base_names <- basename(file_names)
  pattern <- "^([^_]+)_([^_]+)\\+([^_]+)_([^_]+)_exp\\.fcs$"
  matched <- regmatches(base_names, regexec(pattern, base_names))

  unmatched <- base_names[lengths(matched) == 0]
  if (length(unmatched) > 0) {
    stop(
      "These file names do not follow the tutorial pattern ",
      "<donor>_<stim>+<stim>_<well>_exp.fcs: ",
      paste(unmatched, collapse = ", ")
    )
  }

  parts <- do.call(rbind, lapply(matched, function(x) x[2:5]))

  data.frame(
    file_name = base_names,
    donor = parts[, 1],
    condition = paste0(parts[, 2], "+", parts[, 3]),
    stim_1 = parts[, 2],
    stim_2 = parts[, 3],
    well = parts[, 4],
    stimulated = parts[, 2] == "PI" | parts[, 3] == "PI",
    stringsAsFactors = FALSE
  )
}

#' Read a folder of FCS files into a flowSet with its sample sheet attached
#'
#' The sample sheet is derived from the file names by
#' [ParseTutorialFileNames()] and written into the `phenoData` slot, so every
#' later step can group by donor or by condition without a second lookup.
#'
#' @param path Path to a folder that holds the FCS files.
#' @param pattern A regular expression that selects the files. Defaults to any
#'   name that ends in `.fcs`, in either case.
#' @param truncate_max_range Passed to [flowCore::read.flowSet()]. `FALSE` keeps
#'   events that sit above the stated range of a channel, which matters because
#'   the default silently drops them.
#' @return A `flowSet` whose `phenoData` holds the parsed design.
#' @examples
#' \dontrun{
#' ReadTutorialFlowSet("data/datasets/flowjo/.../TQC Basic Tutorial Data")
#' }
#' @export
ReadTutorialFlowSet <- function(path,
                                pattern = "\\.fcs$",
                                truncate_max_range = FALSE) {
  if (!dir.exists(path)) {
    stop("The folder does not exist: ", path)
  }

  files <- list.files(path, pattern = pattern, ignore.case = TRUE,
                      full.names = TRUE)
  if (length(files) == 0) {
    stop("No file in ", path, " matches the pattern ", pattern)
  }

  flow_set <- flowCore::read.flowSet(
    files = files,
    truncate_max_range = truncate_max_range
  )

  sample_sheet <- ParseTutorialFileNames(flowCore::sampleNames(flow_set))
  rownames(sample_sheet) <- flowCore::sampleNames(flow_set)
  flowCore::pData(flow_set) <- sample_sheet

  flow_set
}

#' Build a table of the channels and the markers in a flowFrame
#'
#' The FCS standard stores the detector in `$PnN` and the marker in `$PnS`. A
#' channel with no marker is a scatter or a time channel. Reading both keeps the
#' analysis code free of hard coded detector names.
#'
#' @param frame A `flowFrame`.
#' @return A `data.frame` with the columns `channel`, `marker` and `is_marker`.
#'   `is_marker` is `FALSE` for a channel with no `$PnS` value.
#' @export
DescribeChannels <- function(frame) {
  parameters <- flowCore::pData(flowCore::parameters(frame))

  markers <- as.character(parameters$desc)
  markers[is.na(markers)] <- ""

  data.frame(
    channel = as.character(parameters$name),
    marker = markers,
    is_marker = nzchar(markers),
    stringsAsFactors = FALSE
  )
}

#' Find the channel that carries a named marker
#'
#' @param frame A `flowFrame`.
#' @param marker The marker name to look for, for example `"CD3"`. The match
#'   ignores case and ignores surrounding white space.
#' @return The detector name as a single string.
#' @examples
#' \dontrun{
#' ChannelForMarker(frame, "CD3")
#' }
#' @export
ChannelForMarker <- function(frame, marker) {
  channels <- DescribeChannels(frame)
  wanted <- tolower(trimws(marker))
  hits <- which(tolower(trimws(channels$marker)) == wanted)

  if (length(hits) == 0) {
    stop(
      "No channel carries the marker '", marker, "'. Available markers: ",
      paste(channels$marker[channels$is_marker], collapse = ", ")
    )
  }
  if (length(hits) > 1) {
    stop(
      "More than one channel carries the marker '", marker, "': ",
      paste(channels$channel[hits], collapse = ", ")
    )
  }

  channels$channel[hits]
}
