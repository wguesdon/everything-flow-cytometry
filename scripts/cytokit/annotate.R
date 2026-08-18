#!/usr/bin/env Rscript

# cytokit annotate: put a cell type name on every cluster.
#
# The name comes from the scientist's definitions table and not from the tool. A
# definition says which markers a population is positive for and which it is
# negative for, and every cluster is scored against every definition.
#
# The margin between the best score and the second best is what decides whether
# a label is a fact or a close call. A cluster that beats its runner up by 0.02
# is not labelled, it is guessed, and the report says so.
#
# Called through cli/cytokit, never directly.

suppressPackageStartupMessages({
  library(ggplot2)
})

for (module in c("figures", "clustering", "cytokit")) {
  source(file.path("R", paste0(module, ".R")))
}

arguments <- ParseCytokitArguments(
  commandArgs(trailingOnly = TRUE),
  allowed = c("clusters", "definitions", "out", "label", "margin"),
  required = c("clusters", "definitions"),
  flags = character(0)
)

# A margin below this is reported as a close call. It is not a rule about the
# biology, it is a rule about how much of a lead counts as evidence.
close_call <- if (is.null(arguments$margin)) 0.1 else
  as.numeric(arguments$margin)

Say <- function(...) cat(..., "\n", sep = "")

medians_path <- file.path(arguments$clusters, "cluster_medians.csv")
if (!file.exists(medians_path)) {
  stop("No cluster_medians.csv is in ", DisplayPath(arguments$clusters), ".\n",
       "Point --clusters at the bundle that cytokit cluster wrote.")
}
medians <- utils::read.csv(medians_path, check.names = FALSE,
                           stringsAsFactors = FALSE)

definitions <- tryCatch(ReadCellTypeDefinitions(arguments$definitions),
                        error = function(e) e)
if (inherits(definitions, "error")) {
  stop("The definitions table does not parse: ", conditionMessage(definitions),
       "\nWrite a valid empty one with:  cytokit definitions --data <path> ",
       "--out <file>")
}

# A definition that names a marker the clustering never saw cannot score. When
# none of them overlaps, nothing can be scored at all, and the recipe stops
# before it opens a bundle. A folder that holds no result looks like one.
definition_markers <- setdiff(colnames(definitions), c("cell_type", "note"))
shared <- intersect(definition_markers, colnames(medians))
unmatched <- setdiff(definition_markers, colnames(medians))
if (length(shared) == 0) {
  stop("The definitions and the clustering share no marker.\n",
       "The definitions name: ", paste(definition_markers, collapse = ", "),
       "\nThe clustering carries: ",
       paste(setdiff(colnames(medians),
                     c("cluster", "events", "percent_of_total")),
             collapse = ", "),
       "\nWrite the definitions against the markers the clustering carries, ",
       "or cluster again with --markers.")
}

label <- if (is.null(arguments$label)) ShortLabel(arguments$clusters) else
  arguments$label
out_root <- if (is.null(arguments$out)) kCytokitOutputRoot else arguments$out
bundle <- OpenCytokitBundle("annotate", label, out_root)

Say("cytokit annotate")
Say("  clusters    ", nrow(medians))
Say("  definitions ", nrow(definitions), " cell type(s)")
Say("  bundle      ", DisplayPath(bundle), "\n")

Say("Markers in common: ", length(shared), " of ",
    length(definition_markers), " named by the definitions")
if (length(unmatched) > 0) {
  Say("  These are named by a definition and are not in the clustering:")
  Say("  ", paste(unmatched, collapse = ", "))
  Say("  A definition that rests on one of them scores on the rest, so its")
  Say("  label is weaker than it looks. Cluster on those markers, or take")
  Say("  them out of the definition.")
}

annotation <- AnnotateClusters(medians, definitions)
WriteBundleTable(bundle, annotation, "cluster_labels.csv")

Say("\nOne label per cluster")
print(annotation[, c("cluster", "events", "percent_of_total", "cell_type",
                     "score", "runner_up", "margin")],
      row.names = FALSE, digits = 3)

close <- annotation[annotation$margin < close_call, , drop = FALSE]
if (nrow(close) > 0) {
  WriteBundleTable(bundle, close, "close_calls.csv")
  Say("\nWARNING: ", nrow(close), " of ", nrow(annotation),
      " label(s) beat the runner up by less than ", close_call, ".")
  for (index in seq_len(nrow(close))) {
    Say("  cluster ", close$cluster[index], " reads ",
        close$cell_type[index], " over ", close$runner_up[index],
        " by ", sprintf("%.3f", close$margin[index]))
  }
  Say("  Report one of these as a candidate and not as a cell type.")
} else {
  Say("\nEvery label beats its runner up by at least ", close_call, ".")
}

# The cluster counts per sample carry a cluster number. Putting the label on
# them is what makes a cell type proportion per treatment possible, which is
# the question the whole chain exists to answer.
by_sample_path <- file.path(arguments$clusters, "cluster_counts_by_sample.csv")
if (file.exists(by_sample_path)) {
  by_sample <- utils::read.csv(by_sample_path, check.names = FALSE,
                               stringsAsFactors = FALSE)
  label_of <- stats::setNames(annotation$cell_type,
                              paste0("cluster_", annotation$cluster))
  by_sample$cell_type <- unname(label_of[by_sample$population])
  named <- by_sample[!is.na(by_sample$cell_type), , drop = FALSE]
  if (nrow(named) > 0) {
    totals <- stats::aggregate(count ~ sample + cell_type, data = named,
                               FUN = sum)
    per_sample_total <- stats::aggregate(count ~ sample, data = by_sample,
                                         FUN = sum)
    totals <- merge(totals, per_sample_total, by = "sample",
                    suffixes = c("", "_parent"))
    cell_type_counts <- data.frame(
      sample = totals$sample,
      population = totals$cell_type,
      count = totals$count,
      percent_of_parent = 100 * totals$count / totals$count_parent,
      stringsAsFactors = FALSE)
    cell_type_counts <- cell_type_counts[
      order(cell_type_counts$sample, cell_type_counts$population), ,
      drop = FALSE]
    WriteBundleTable(bundle, cell_type_counts,
                     "cell_type_counts_by_sample.csv")
    Say("\nOne row per sample and cell type is in ",
        "cell_type_counts_by_sample.csv,")
    Say("over ", length(unique(cell_type_counts$sample)), " sample(s). ",
        "Compare a cell type between treatments with:")
    Say("  cytokit proportions --counts ",
        DisplayPath(file.path(bundle, "cell_type_counts_by_sample.csv")),
        " --metadata <file>")
  }
} else {
  Say("\nThe cluster bundle holds no per sample count, so a cell type")
  Say("proportion per treatment is not reachable from it. Run cytokit cluster")
  Say("without --one-sample to get one.")
}

summary_table <- SummariseCellTypes(annotation)
WriteBundleTable(bundle, summary_table, "cell_type_summary.csv")
Say("\nOne row per cell type")
print(summary_table, row.names = FALSE, digits = 4)

# A definition that won nothing is either wrong for this panel or the population
# is absent. Both are worth knowing, and neither shows in the table above.
missing <- setdiff(definitions$cell_type, annotation$cell_type)
if (length(missing) > 0) {
  Say("\n", length(missing), " definition(s) won no cluster: ",
      paste(missing, collapse = ", "))
  Say("  Either the population is absent, or the definition does not suit")
  Say("  this panel. Raise --metaclusters on cytokit cluster to split a")
  Say("  cluster that holds two populations.")
}

drawing <- ggplot2::ggplot(
  summary_table,
  ggplot2::aes(x = stats::reorder(.data$cell_type, .data$percent_of_total),
               y = .data$percent_of_total, fill = .data$cell_type)) +
  ggplot2::geom_col(show.legend = FALSE) +
  ggplot2::coord_flip() +
  ScaleFillPublication() +
  ggplot2::labs(title = "Share of the clustered events",
                x = NULL, y = "Percent of events") +
  ThemePublication()
SaveFigure(drawing, file.path(bundle, "cell_type_summary.svg"),
           width = 8, height = 2 + 0.4 * nrow(summary_table))

CloseCytokitBundle(
  bundle, "annotate", arguments,
  inputs = c(medians_path, arguments$definitions),
  command = paste("cytokit annotate --clusters",
                  DisplayPath(arguments$clusters), "--definitions",
                  DisplayPath(arguments$definitions))
)

Say("\nWrote cluster_labels.csv and cell_type_summary.csv to ",
    DisplayPath(bundle))
