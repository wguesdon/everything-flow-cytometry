# Automated gating.
#
# The gate hierarchy is not written in this file. It is written in
# gating/pbmc_gating_template.csv, and openCyto reads that file and fits every
# gate from the data. The separation is the point of the example. A gate drawn by
# hand lives in a workspace and cannot be read in a diff, while a template is a
# text file that any reviewer can read and any version control system can track.
#
# Every gate is still fitted per sample. The template fixes the method and the
# channels, not the cut point, so a sample with a different distribution still
# gets its own boundary.

#' Read and check an openCyto gating template
#'
#' @param path Path to the template CSV.
#' @return A `gatingTemplate` object.
#' @examples
#' \dontrun{
#' ReadGatingTemplate("gating/pbmc_gating_template.csv")
#' }
#' @export
ReadGatingTemplate <- function(path) {
  if (!file.exists(path)) {
    stop("The gating template does not exist: ", path)
  }

  required <- c("alias", "pop", "parent", "dims", "gating_method")
  header <- utils::read.csv(path, nrows = 1, check.names = FALSE)
  missing <- setdiff(required, colnames(header))
  if (length(missing) > 0) {
    stop(
      "The gating template is missing the column(s): ",
      paste(missing, collapse = ", "), "."
    )
  }

  openCyto::gatingTemplate(path)
}

#' Run an openCyto template against a transformed flowSet
#'
#' @param flow_set A compensated and transformed `flowSet`.
#' @param template A `gatingTemplate`, from [ReadGatingTemplate()].
#' @param n_cores The number of cores openCyto may use. Defaults to 1, because a
#'   higher value makes the error message from a failed gate harder to read.
#' @return A `GatingSet` with every population of the template applied.
#' @export
RunAutomatedGating <- function(flow_set, template, n_cores = 1) {
  if (!methods::is(flow_set, "flowSet")) {
    stop("flow_set must be a flowSet, not a ", class(flow_set)[1], ".")
  }

  gating_set <- flowWorkspace::GatingSet(flow_set)
  openCyto::gt_gating(template, gating_set, mc.cores = n_cores, parallel_type = "none")
  flowWorkspace::recompute(gating_set)

  gating_set
}

#' Collect the event count and the percentage for every population
#'
#' @param gating_set A gated `GatingSet`.
#' @param sample_sheet A `data.frame` with a `file_name` column, joined onto the
#'   result so the design travels with the numbers. Pass `NULL` to skip the join.
#' @return A `data.frame` with one row per sample and population, holding
#'   `sample`, `population`, `count`, `parent_count` and `percent_of_parent`.
#' @export
CollectPopulationStats <- function(gating_set, sample_sheet = NULL) {
  counts <- flowWorkspace::gs_pop_get_stats(gating_set, type = "count")
  percents <- flowWorkspace::gs_pop_get_stats(gating_set, type = "percent")

  counts <- as.data.frame(counts)
  percents <- as.data.frame(percents)

  merged <- merge(
    counts, percents,
    by = c("sample", "pop"),
    suffixes = c("_count", "_percent")
  )

  result <- data.frame(
    sample = merged$sample,
    population = merged$pop,
    count = merged$count,
    percent_of_parent = merged$percent * 100,
    stringsAsFactors = FALSE
  )
  result <- result[order(result$sample, result$population), ]

  if (!is.null(sample_sheet)) {
    result <- merge(
      result, sample_sheet,
      by.x = "sample", by.y = "file_name",
      all.x = TRUE
    )
  }

  rownames(result) <- NULL
  result
}

#' Report the spread of a population frequency across samples
#'
#' This is the number the whole example exists to show. When one automated
#' template gates every sample, the coefficient of variation of a population
#' frequency measures the biological and technical spread only. When each analyst
#' draws their own gates, the analyst is added to that spread.
#'
#' @param stats The output of [CollectPopulationStats()].
#' @param group_by A column name to split by before the CV is computed, for
#'   example `"condition"`. Pass `NULL` to pool every sample.
#' @return A `data.frame` with the columns `population`, the grouping column when
#'   one is given, `n`, `mean_percent`, `sd_percent` and `cv_percent`. The CV is
#'   reported as a percentage of the mean.
#' @export
SummarisePopulationSpread <- function(stats, group_by = NULL) {
  required <- c("population", "percent_of_parent")
  missing <- setdiff(required, colnames(stats))
  if (length(missing) > 0) {
    stop("stats is missing the column(s): ", paste(missing, collapse = ", "), ".")
  }

  keys <- c("population", group_by)
  unknown <- setdiff(keys, colnames(stats))
  if (length(unknown) > 0) {
    stop("stats has no column named ", paste(unknown, collapse = ", "), ".")
  }

  split_key <- interaction(stats[keys], drop = TRUE, sep = "\r")
  pieces <- split(stats, split_key)

  rows <- lapply(pieces, function(piece) {
    values <- piece$percent_of_parent
    mean_value <- mean(values)
    sd_value <- stats::sd(values)

    out <- piece[1, keys, drop = FALSE]
    out$n <- length(values)
    out$mean_percent <- mean_value
    out$sd_percent <- sd_value
    out$cv_percent <- if (mean_value == 0) NA_real_ else 100 * sd_value / mean_value
    out
  })

  result <- do.call(rbind, rows)
  rownames(result) <- NULL
  result[order(result$population), ]
}

#' Read the hierarchy of a gated set as a table with a parent column
#'
#' `CollectPopulationStats` reports one row per sample and population, which is
#' what a comparison needs. A drawing of the hierarchy needs the parent of each
#' population as well, and [PlotGateTree()] takes it in that shape.
#'
#' The percentage is of the parent, so a population that keeps every event of
#' its parent reads 100.
#'
#' @param gating_set A gated `GatingSet`.
#' @param sample The sample to read, by name or by index. Defaults to the first.
#' @return A `data.frame` with `population`, `parent`, `events`,
#'   `parent_events` and `percent_of_parent`, in the order the hierarchy walks.
#'   The first row is `all_events`, whose parent is `NA`, and it holds the count
#'   every gate starts from.
#' @examples
#' \dontrun{
#' CollectGateTree(gating_set)
#' }
#' @export
CollectGateTree <- function(gating_set, sample = 1) {
  paths <- flowWorkspace::gs_get_pop_paths(gating_set, path = "full")
  paths <- paths[paths != "root"]
  if (length(paths) == 0) {
    stop("The gating set holds no population below root.")
  }

  handle <- gating_set[[sample]]
  events <- vapply(paths, function(path) {
    as.numeric(flowWorkspace::gh_pop_get_count(handle, path))
  }, numeric(1))
  root_events <- as.numeric(
    flowWorkspace::gh_pop_get_count(handle, "root"))

  parent_path <- dirname(paths)
  # The root of the drawing carries no parent. PlotGateTree reads NA as the
  # root, and the root row gives a reader the count every gate starts from.
  parent <- ifelse(parent_path == "/", "all_events", basename(parent_path))
  parent_events <- ifelse(parent_path == "/", root_events,
                          events[match(parent_path, paths)])

  rbind(
    data.frame(population = "all_events", parent = NA_character_,
               events = root_events, parent_events = root_events,
               percent_of_parent = 100, stringsAsFactors = FALSE),
    data.frame(
      population = basename(paths),
      parent = parent,
      events = events,
      parent_events = parent_events,
      percent_of_parent = 100 * events / pmax(parent_events, 1),
      stringsAsFactors = FALSE
    )
  )
}
