# Compare an automated gating result against a manual one.
#
# This is the positive control for the whole repository. OMIP-39 ships a FlowJo
# workspace with the gating that the authors published, so the manual result is
# not a result this project produced. CytoML reads those gates, openCyto gates
# the
# same file from a template, and the two frequencies sit side by side.
#
# The comparison is honest about what it can show. One donor file gives
# agreement
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
#' mapping is a judgement, so it lives in a file that a reviewer can check
#' rather
#' than inside a function.
#'
#' @param path Path to a CSV with the columns `manual` and `automated`. A row
#' may
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
#' @param automated_stats The output of [CollectPopulationStats()] on the
#' openCyto
#'   `GatingSet`.
#' @param population_map The output of [ReadPopulationMap()].
#' @return A `data.frame` with one row per mapped population, holding `manual`,
#'   `automated`, `manual_percent`, `automated_percent`, `difference` and
#'   `ratio`. `difference` is in percentage points. A population that either
#' side
#'   did not produce keeps `NA` on that side, so a gap is visible.
#' @export
CompareGatingResults <- function(manual_stats, automated_stats,
                                 population_map) {
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
#' @return The input with two columns added. `agrees` is `TRUE`, `FALSE` or
#' `NA`,
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

#' Join a population count table to a sample metadata table
#'
#' A count is not a result. The result is the share of a parent that a
#' population holds in each sample, read beside what was done to that sample.
#' An FCS file carries no treatment, so the metadata comes from the scientist.
#'
#' A sample that the metadata does not name is kept and marked, because dropping
#' it silently is how a group loses a replicate without anybody noticing.
#'
#' @param stats A count table with the columns `sample`, `population` and
#'   `percent_of_parent`, as [CollectPopulationStats()] returns.
#' @param metadata A `data.frame` whose `sample_column` names the file.
#' @param sample_column The column of `metadata` that names the file. Defaults
#'   to `"sample"`.
#' @return A `data.frame` with one row per sample and population, carrying every
#'   metadata column and a `metadata_found` flag.
#' @examples
#' stats <- data.frame(sample = c("a.fcs", "b.fcs"), population = "T cells",
#'                     count = c(10, 20), percent_of_parent = c(30, 40))
#' metadata <- data.frame(sample = c("a.fcs", "b.fcs"),
#'                        treatment = c("control", "drug"))
#' PopulationProportions(stats, metadata)
#' @export
PopulationProportions <- function(stats, metadata, sample_column = "sample") {
  required <- c("sample", "population", "percent_of_parent")
  missing <- setdiff(required, colnames(stats))
  if (length(missing) > 0) {
    stop("stats is missing the column(s): ", paste(missing, collapse = ", "),
         ".")
  }
  if (!sample_column %in% colnames(metadata)) {
    stop("The metadata has no column called '", sample_column, "'.\n",
         "It carries: ", paste(colnames(metadata), collapse = ", "))
  }
  if (anyDuplicated(metadata[[sample_column]]) > 0) {
    duplicated_names <- unique(metadata[[sample_column]][
      duplicated(metadata[[sample_column]])])
    stop("The metadata names a sample twice: ",
         paste(duplicated_names, collapse = ", "),
         ".\nOne row of metadata describes one file.")
  }

  # A file name can carry a folder or an extension on one side and not the
  # other, so both are compared on the base name.
  key <- basename(as.character(stats$sample))
  metadata_key <- basename(as.character(metadata[[sample_column]]))
  position <- match(key, metadata_key)

  extra <- setdiff(colnames(metadata), sample_column)
  joined <- stats
  for (column in extra) {
    joined[[column]] <- metadata[[column]][position]
  }
  joined$metadata_found <- !is.na(position)
  rownames(joined) <- NULL
  joined
}

#' Compare one population's frequency between groups, and draw it
#'
#' The test is chosen from the design and not from the result. Two groups get a
#' Wilcoxon rank sum test and more than two get a Kruskal Wallis test, because a
#' cell frequency is bounded, skewed, and measured on a handful of samples.
#'
#' A group with one sample carries no spread, so no test runs on it. The
#' function
#' says so and still draws the points, because a picture of five samples is
#' useful and a p value over one sample per group is not.
#'
#' @param proportions The output of [PopulationProportions()].
#' @param population The population to compare.
#' @param group The metadata column that splits the samples.
#' @param value The column holding the frequency. Defaults to
#'   `"percent_of_parent"`.
#' @return A list with `data`, the rows compared, `summary`, one row per group,
#'   `test`, a one row `data.frame`, and `plot`, a `ggplot`.
#' @examples
#' \dontrun{
#' CompareProportions(proportions, "T cells", "treatment")
#' }
#' @export
CompareProportions <- function(proportions, population, group,
                               value = "percent_of_parent") {
  for (column in c("population", group, value)) {
    if (!column %in% colnames(proportions)) {
      stop("The table has no column called '", column, "'.\n",
           "It carries: ", paste(colnames(proportions), collapse = ", "))
    }
  }
  rows <- proportions[proportions$population == population, , drop = FALSE]
  if (nrow(rows) == 0) {
    stop("No row holds the population '", population, "'.\n",
         "The table holds: ",
         paste(utils::head(unique(proportions$population), 20),
               collapse = ", "))
  }
  rows <- rows[!is.na(rows[[group]]) & !is.na(rows[[value]]), , drop = FALSE]
  if (nrow(rows) == 0) {
    stop("Every row of '", population, "' has no value or no group.")
  }

  rows[[group]] <- factor(rows[[group]])
  counts <- table(rows[[group]])
  summary_table <- data.frame(
    group = names(counts),
    samples = as.integer(counts),
    median = as.numeric(tapply(rows[[value]], rows[[group]], stats::median)),
    mean = as.numeric(tapply(rows[[value]], rows[[group]], mean)),
    sd = as.numeric(tapply(rows[[value]], rows[[group]], stats::sd)),
    stringsAsFactors = FALSE
  )

  # The design decides whether a test is possible, and the design is the number
  # of samples per group.
  testable <- length(counts) >= 2 && all(counts >= 2)
  test <- if (!testable) {
    reason <- if (length(counts) < 2) {
      paste0("there is one group, '", names(counts)[1], "'")
    } else {
      paste0("these group(s) hold one sample: ",
             paste(names(counts)[counts < 2], collapse = ", "))
    }
    data.frame(test = "none", statistic = NA_real_, p_value = NA_real_,
               groups = length(counts), reason = reason,
               stringsAsFactors = FALSE)
  } else if (length(counts) == 2) {
    result <- stats::wilcox.test(rows[[value]] ~ rows[[group]], exact = FALSE)
    data.frame(test = "Wilcoxon rank sum",
               statistic = unname(result$statistic),
               p_value = result$p.value, groups = 2L,
               reason = NA_character_, stringsAsFactors = FALSE)
  } else {
    result <- stats::kruskal.test(rows[[value]] ~ rows[[group]])
    data.frame(test = "Kruskal Wallis", statistic = unname(result$statistic),
               p_value = result$p.value, groups = length(counts),
               reason = NA_character_, stringsAsFactors = FALSE)
  }

  subtitle <- if (testable) {
    sprintf("%s, p = %.4g", test$test, test$p_value)
  } else {
    paste0("No test ran, because ", test$reason)
  }

  # A box over three points draws a quartile that three points cannot support,
  # so the points are always drawn on top of it.
  drawing <- ggplot2::ggplot(
    rows, ggplot2::aes(x = .data[[group]], y = .data[[value]],
                       fill = .data[[group]])) +
    ggplot2::geom_boxplot(outlier.shape = NA, alpha = 0.55,
                          show.legend = FALSE) +
    ggplot2::geom_jitter(width = 0.12, height = 0, size = 2.2, alpha = 0.9,
                         show.legend = FALSE) +
    ScaleFillPublication() +
    ggplot2::labs(title = population, subtitle = subtitle, x = group,
                  y = "Percent of parent") +
    ThemePublication()

  list(data = rows, summary = summary_table, test = test, plot = drawing)
}
