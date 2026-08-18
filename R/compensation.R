# Compensation.
#
# Spillover is the light from one fluorochrome that a second detector also
# reads.
# The correction is a linear one: measured = true %*% spillover, so the
# correction
# multiplies the measured values by the inverse of the spillover matrix.
#
# The FlowJo tutorial files carry the matrix in the FCS keyword, so the matrix
# is
# read and applied rather than computed. When single stained controls are the
# only
# source, flowStats::spillover_ng() computes one instead.

#' Read the spillover matrix out of an FCS file
#'
#' A cytometer writes the matrix into one of three keywords, and which one it
#' uses
#' depends on the acquisition software. [flowCore::spillover()] returns all
#' three
#' slots, of which at most one is filled, and it raises an error when the file
#' carries no matrix at all. Both outcomes are turned into one message here, so
#' a
#' caller does not have to know which shape flowCore chose.
#'
#' @param frame A `flowFrame`.
#' @return A numeric matrix whose row and column names are detector names.
#' @examples
#' \dontrun{
#' ExtractSpillover(frame)
#' }
#' @export
ExtractSpillover <- function(frame) {
  slots <- tryCatch(
    flowCore::spillover(frame),
    error = function(e) NULL
  )
  filled <- if (is.null(slots)) list() else Filter(Negate(is.null), slots)

  if (length(filled) == 0) {
    stop(
      "This flowFrame carries no spillover matrix. ",
      "Compute one from single stained controls with ",
      "ComputeSpilloverFromControls() ",
      "or with flowStats::spillover_ng()."
    )
  }

  matrix_out <- filled[[1]]
  if (!is.matrix(matrix_out)) {
    stop("The spillover keyword holds a ", class(matrix_out)[1],
         ", not a matrix.")
  }
  if (nrow(matrix_out) != ncol(matrix_out)) {
    stop(
      "The spillover matrix is not square. It has ", nrow(matrix_out),
      " rows and ", ncol(matrix_out), " columns."
    )
  }

  # A stored matrix often carries column names and no row names, because the
  # acquisition software writes the detector list once. The matrix is square
  # over
  # the same detectors, so the row names are the column names. Filling them here
  # means no caller downstream has to test for a NULL dimname.
  if (is.null(rownames(matrix_out)) && !is.null(colnames(matrix_out))) {
    rownames(matrix_out) <- colnames(matrix_out)
  }
  if (is.null(colnames(matrix_out)) && !is.null(rownames(matrix_out))) {
    colnames(matrix_out) <- rownames(matrix_out)
  }
  if (is.null(rownames(matrix_out))) {
    stop("The spillover matrix carries no detector names on either dimension.")
  }

  matrix_out
}

#' Check that a spillover matrix can be applied to a flowFrame
#'
#' Compensation fails when a matrix names a detector that the file does not
#' carry.
#' The usual cause is a name that the acquisition software wrote differently,
#' for
#' example `FITC-A` in the file against `FITC.A` in the matrix. This check
#' reports
#' the difference before [flowCore::compensate()] raises a less readable error.
#'
#' @param frame A `flowFrame`.
#' @param spillover A spillover matrix, as returned by [ExtractSpillover()].
#' @return `TRUE`, invisibly, when every matrix column is a channel of the
#' frame.
#' @export
CheckSpilloverChannels <- function(frame, spillover) {
  available <- flowCore::colnames(frame)
  missing <- setdiff(colnames(spillover), available)

  if (length(missing) > 0) {
    stop(
      "The spillover matrix names ", length(missing),
      " detector(s) that the file does not carry: ",
      paste(missing, collapse = ", "), ". ",
      "The file carries: ", paste(available, collapse = ", "), ". ",
      "Rename the matrix columns to match the file before you compensate."
    )
  }

  invisible(TRUE)
}

#' Apply the spillover matrix that the file carries
#'
#' @param x A `flowFrame` or a `flowSet`. For a `flowSet` the matrix is read
#' from
#'   the frame named by `reference`, because every file in one acquisition
#' shares
#'   one matrix.
#' @param reference For a `flowSet`, the name or the index of the frame whose
#'   matrix is used. Defaults to the first frame.
#' @return A compensated object of the same class as `x`.
#' @examples
#' \dontrun{
#' compensated <- ApplyCompensation(flow_set)
#' }
#' @export
ApplyCompensation <- function(x, reference = 1) {
  if (methods::is(x, "flowFrame")) {
    spillover <- ExtractSpillover(x)
    CheckSpilloverChannels(x, spillover)
    return(flowCore::compensate(x, spillover))
  }

  if (methods::is(x, "flowSet")) {
    spillover <- ExtractSpillover(x[[reference]])
    CheckSpilloverChannels(x[[reference]], spillover)
    return(flowCore::compensate(x, spillover))
  }

  stop("x must be a flowFrame or a flowSet, not a ", class(x)[1], ".")
}

#' Report how much spillover a matrix corrects
#'
#' The diagonal of a spillover matrix is 1 by construction. Every value off the
#' diagonal is the fraction of one fluorochrome that a second detector reads.
#' The
#' largest of those values tells you which pair of detectors the compensation
#' actually moves, which is the pair worth checking on a plot.
#'
#' @param spillover A spillover matrix.
#' @param top The number of rows to return. Defaults to 10.
#' @return A `data.frame` with the columns `from`, `to` and `spill`, sorted by
#'   `spill` in descending order. `spill` is a fraction, so 0.23 is 23 percent.
#' @export
SummariseSpillover <- function(spillover, top = 10) {
  if (!is.matrix(spillover)) {
    stop("spillover must be a matrix, not a ", class(spillover)[1], ".")
  }

  row_names <- rownames(spillover)
  col_names <- colnames(spillover)
  if (is.null(row_names)) row_names <- col_names
  if (is.null(col_names)) col_names <- row_names
  if (is.null(row_names)) {
    stop("spillover carries no detector names, so a pair cannot be named.")
  }

  long <- expand.grid(
    from = row_names,
    to = col_names,
    stringsAsFactors = FALSE
  )
  long$spill <- as.vector(spillover)
  off_diagonal <- long[long$from != long$to, ]
  ordered <- off_diagonal[order(-off_diagonal$spill), ]

  head(ordered, top)
}
