# Compare an automated gating result against a manual one.
#
# This is the positive control for the whole repository. OMIP-39 ships a FlowJo
# workspace with the gating that the authors published, so the manual result is
# not a result this project produced. CytoML reads those gates, openCyto gates the
# same file from a template, and the two frequencies sit side by side.
#
# The comparison is honest about what it can show. One donor file gives agreement
# per gate. It does not give a coefficient of variation across samples, and it
# does not show that either method is right. It shows whether an automated
# template lands where an expert landed.

#' Read the manual gates out of a FlowJo workspace
#'
#' @param workspace_path Path to the `.wsp` file.
#' @param fcs_path Folder that holds the FCS files the workspace refers to.
#' @param group The sample group inside the workspace, for example `"Sample"`.
#' @return A `GatingSet` carrying the gates the analyst drew.
#' @examples
#' \dontrun{
#' ImportFlowJoGates("attachments/OMIP_Hammer.wsp", "data/omip39", "Sample")
#' }
#' @export
ImportFlowJoGates <- function(workspace_path, fcs_path, group = "Sample") {
  if (!file.exists(workspace_path)) {
    stop("The FlowJo workspace does not exist: ", workspace_path)
  }
  if (!dir.exists(fcs_path)) {
    stop("The FCS folder does not exist: ", fcs_path)
  }

  workspace <- CytoML::open_flowjo_xml(workspace_path)
  CytoML::flowjo_to_gatingset(workspace, name = group, path = fcs_path)
}

#' Read a mapping between manual population names and template aliases
#'
#' The two gating strategies name the same population differently. The analyst
#' wrote `CD56+ CD3-` in FlowJo and the template writes `CD56posCD3neg`. That
#' mapping is a judgement, so it lives in a file that a reviewer can check rather
#' than inside a function.
#'
#' @param path Path to a CSV with the columns `manual` and `automated`. A row may
#'   also carry a `note` column, which is ignored here and printed by a report.
#' @return A `data.frame` with the columns `manual`, `automated` and `note`.
#' @export
ReadPopulationMap <- function(path) {
  if (!file.exists(path)) {
    stop("The population map does not exist: ", path)
  }

  map <- utils::read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
  required <- c("manual", "automated")
  missing <- setdiff(required, colnames(map))
  if (length(missing) > 0) {
    stop("The population map is missing the column(s): ",
         paste(missing, collapse = ", "), ".")
  }

  if (!"note" %in% colnames(map)) {
    map$note <- ""
  }

  duplicated_manual <- map$manual[duplicated(map$manual)]
  if (length(duplicated_manual) > 0) {
    stop("The population map names a manual population twice: ",
         paste(unique(duplicated_manual), collapse = ", "), ".")
  }

  map[, c("manual", "automated", "note")]
}

#' Put the manual and the automated frequency of each population side by side
#'
#' @param manual_stats The output of [CollectPopulationStats()] on the imported
#'   FlowJo `GatingSet`.
#' @param automated_stats The output of [CollectPopulationStats()] on the openCyto
#'   `GatingSet`.
#' @param population_map The output of [ReadPopulationMap()].
#' @return A `data.frame` with one row per mapped population, holding `manual`,
#'   `automated`, `manual_percent`, `automated_percent`, `difference` and
#'   `ratio`. `difference` is in percentage points. A population that either side
#'   did not produce keeps `NA` on that side, so a gap is visible.
#' @export
CompareGatingResults <- function(manual_stats, automated_stats, population_map) {
  LeafName <- function(paths) {
    basename(as.character(paths))
  }

  manual_lookup <- stats::setNames(
    manual_stats$percent_of_parent,
    LeafName(manual_stats$population)
  )
  automated_lookup <- stats::setNames(
    automated_stats$percent_of_parent,
    LeafName(automated_stats$population)
  )

  result <- data.frame(
    manual = population_map$manual,
    automated = population_map$automated,
    note = population_map$note,
    manual_percent = unname(manual_lookup[population_map$manual]),
    automated_percent = unname(automated_lookup[population_map$automated]),
    stringsAsFactors = FALSE
  )

  result$difference <- result$automated_percent - result$manual_percent
  result$ratio <- result$automated_percent / result$manual_percent

  result
}

#' Judge each population against an agreement threshold
#'
#' A comparison needs a stated rule, decided before the numbers are read. The
#' default of 5 percentage points is the tolerance used here, and it is a
#' judgement rather than a published standard, so a report must say so.
#'
#' @param comparison The output of [CompareGatingResults()].
#' @param tolerance_points The absolute difference, in percentage points, at or
#'   below which the two methods are called in agreement.
#' @return The input with two columns added. `agrees` is `TRUE`, `FALSE` or `NA`,
#'   and `verdict` is `"agree"`, `"differ"` or `"missing"`.
#' @export
JudgeAgreement <- function(comparison, tolerance_points = 5) {
  if (!"difference" %in% colnames(comparison)) {
    stop("comparison has no 'difference' column. Run CompareGatingResults first.")
  }
  if (tolerance_points <= 0) {
    stop("tolerance_points must be above zero, not ", tolerance_points, ".")
  }

  comparison$agrees <- abs(comparison$difference) <= tolerance_points
  comparison$verdict <- ifelse(
    is.na(comparison$difference), "missing",
    ifelse(comparison$agrees, "agree", "differ")
  )

  comparison
}
