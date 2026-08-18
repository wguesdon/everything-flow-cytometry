# Transformation.
#
# Fluorescence data spans several orders of magnitude and holds negative values
# after compensation, so a plain log transform fails on the negative part. The
# logicle transform is linear near zero and logarithmic above it, which keeps
# the
# negative population visible and keeps the positive population separated.
#
# Estimate the transform on one file and apply the same transform to every file.
# A per file estimate would move a gate between samples for a reason that has
# nothing to do with the biology.

#' Select the fluorescence channels of a flowFrame
#'
#' Scatter, time and any channel with no marker are excluded. Only a
#' fluorescence
#' channel needs a logicle transform.
#'
#' @param frame A `flowFrame`.
#' @param exclude Extra channel names to drop, for example a channel whose
#' marker
#'   is a blank. Matching ignores case.
#' @return A character vector of detector names.
#' @export
FluorescenceChannels <- function(frame, exclude = character()) {
  channels <- DescribeChannels(frame)

  is_scatter <- grepl("^(FSC|SSC)", channels$channel, ignore.case = TRUE)
  is_time <- grepl("^time$", channels$channel, ignore.case = TRUE)
  is_excluded <- tolower(channels$channel) %in% tolower(exclude)

  keep <- channels$is_marker & !is_scatter & !is_time & !is_excluded
  channels$channel[keep]
}

#' Estimate one logicle transform and apply it to every sample
#'
#' The transform is estimated on a single frame and reused, so a gate drawn on
#' one
#' sample means the same thing on the next one.
#'
#' @param x A `flowFrame` or a compensated `flowSet`.
#' @param channels The detector names to transform. Defaults to every
#'   fluorescence channel of the reference frame.
#' @param reference For a `flowSet`, the name or the index of the frame used to
#'   estimate the transform. Defaults to the first frame.
#' @return A list with two elements. `data` is the transformed object and
#'   `transform` is the `transformList` that was applied, kept so a report can
#'   state exactly what was done.
#' @examples
#' \dontrun{
#' result <- ApplyLogicleTransform(compensated_set)
#' transformed <- result$data
#' }
#' @export
ApplyLogicleTransform <- function(x, channels = NULL, reference = 1) {
  is_frame <- methods::is(x, "flowFrame")
  is_set <- methods::is(x, "flowSet")

  if (!is_frame && !is_set) {
    stop("x must be a flowFrame or a flowSet, not a ", class(x)[1], ".")
  }

  reference_frame <- if (is_frame) x else x[[reference]]

  if (is.null(channels)) {
    channels <- FluorescenceChannels(reference_frame)
  }
  if (length(channels) == 0) {
    stop("No fluorescence channel was selected, so there is nothing to transform.")
  }

  unknown <- setdiff(channels, flowCore::colnames(reference_frame))
  if (length(unknown) > 0) {
    stop(
      "These channels are not in the data: ", paste(unknown,
                                                    collapse = ", "), "."
    )
  }

  transform_list <- flowCore::estimateLogicle(reference_frame,
                                              channels = channels)

  list(
    data = flowCore::transform(x, transform_list),
    transform = transform_list,
    channels = channels
  )
}
