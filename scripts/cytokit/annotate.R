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

label <- if (is.null(arguments$label)) ShortLabel(arguments$clusters) else
  arguments$label
out_root <- if (is.null(arguments$out)) kCytokitOutputRoot else arguments$out
bundle <- OpenCytokitBundle("annotate", label, out_root)

Say("cytokit annotate")
Say("  clusters    ", nrow(medians))
Say("  definitions ", nrow(definitions), " cell type(s)")
Say("  bundle      ", DisplayPath(bundle), "\n")

# A definition that names a marker the clustering never saw cannot score, so the
# overlap is reported before the scores are read.
definition_markers <- setdiff(colnames(definitions), c("cell_type", "note"))
shared <- intersect(definition_markers, colnames(medians))
unmatched <- setdiff(definition_markers, colnames(medians))
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
