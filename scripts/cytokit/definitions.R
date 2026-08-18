#!/usr/bin/env Rscript

# cytokit definitions: write an empty cell type definitions table for a panel.
#
# A definitions table turns a cluster number into a cell type name. It has one
# column per marker and one row per population, and each cell says what that
# population does with that marker: pos, neg, high, or empty for "does not
# matter". AnnotateClusters scores every cluster against every row and reports
# the margin over the runner up.
#
# The columns come from the panel, so an agent fills a valid table rather than
# inventing a marker name that no detector carries.
#
# Called through cli/cytokit, never directly.

suppressPackageStartupMessages({
  library(flowCore)
})

for (module in c("io", "clustering", "panels", "cytokit")) {
  source(file.path("R", paste0(module, ".R")))
}

arguments <- ParseCytokitArguments(
  commandArgs(trailingOnly = TRUE),
  allowed = c("data", "out", "markers", "populations"),
  required = c("data", "out"),
  flags = c("force", "recursive")
)

if (file.exists(arguments$out) && !isTRUE(arguments$force)) {
  stop("The file already exists: ", arguments$out,
       "\nAdd --force to overwrite it. A definitions table holds decisions ",
       "that took work, so it is not overwritten by default.")
}

files <- FcsFilesIn(arguments$data,
                    recursive = isTRUE(arguments$recursive))
panel <- DescribeFcsPanel(files[1])
available <- PanelMarkers(panel)

Say <- function(...) cat(..., "\n", sep = "")

# A definitions table with 28 columns is unreadable, and a population is defined
# by a handful of markers. --markers narrows it, and on a panel that names no
# antibody the same argument carries the mapping instead.
source_of <- PanelMarkerSource(panel)
chosen <- ParseMarkerArgument(arguments$markers, available, panel$channel)
if (is.null(chosen)) {
  chosen <- data.frame(column = available, from = available,
                       stringsAsFactors = FALSE)
}
markers <- chosen$column

# A column that still carries a detector name is the fault this recipe has to
# state. A cell type is a claim about an antibody, and APC-A is not one.
chosen$source <- source_of$source[match(chosen$from, source_of$marker)]
chosen$source[chosen$from != chosen$column] <- "antibody"
chosen$source[is.na(chosen$source)] <- "antibody"
guessed <- chosen[chosen$source == "detector", , drop = FALSE]

populations <- if (is.null(arguments$populations)) {
  character(0)
} else {
  trimws(strsplit(arguments$populations, ",")[[1]])
}

header <- c("cell_type", markers, "note")
rows <- vapply(populations, function(population) {
  paste(c(population, rep("", length(markers)), ""), collapse = ",")
}, character(1), USE.NAMES = FALSE)

# The CSV carries no comment lines. ReadCellTypeDefinitions calls read.csv,
# whose comment character is empty by default, so a leading comment becomes the
# header and the table stops parsing. The guidance goes in a sibling file.
dir.create(dirname(arguments$out), recursive = TRUE, showWarnings = FALSE)
writeLines(c(paste(header, collapse = ","), rows), arguments$out)

notes_path <- sub("\\.csv$", "_notes.md", arguments$out)
if (identical(notes_path, arguments$out)) {
  notes_path <- paste0(arguments$out, "_notes.md")
}
writeLines(c(
  paste0("# Cell type definitions for ", basename(files[1])),
  "",
  paste0("Written by `cytokit definitions` beside `", basename(arguments$out),
         "`."),
  "",
  "One row is one population. A marker cell takes `pos`, `neg` or `high`, or",
  "stays empty when that marker does not decide the population.",
  "`AnnotateClusters` scores every cluster against every row and reports the",
  "margin over the runner up, so a small margin means the label is a close",
  "call and not a fact.",
  "",
  "Name only the markers that define a population. A row that names more",
  "markers is not scored higher for it, because the score is a weighted mean",
  "and not a sum.",
  "",
  if (nrow(guessed) > 0) "## Read this before you fill a cell" else NULL,
  if (nrow(guessed) > 0) "" else NULL,
  if (nrow(guessed) > 0) {
    paste0(nrow(guessed), " of the ", length(markers), " marker columns below ",
           "are detector names and not antibodies. The file leaves `$PnS`")
  } else NULL,
  if (nrow(guessed) > 0) {
    "empty for them, so nothing in the data says what was stained."
  } else NULL,
  if (nrow(guessed) > 0) "" else NULL,
  if (nrow(guessed) > 0) paste0("- ", guessed$column) else NULL,
  if (nrow(guessed) > 0) "" else NULL,
  if (nrow(guessed) > 0) {
    "A cell type is a claim about an antibody. Ask the scientist which antibody"
  } else NULL,
  if (nrow(guessed) > 0) {
    "sits on each detector, then write the file again with the mapping:"
  } else NULL,
  if (nrow(guessed) > 0) "" else NULL,
  if (nrow(guessed) > 0) "```" else NULL,
  if (nrow(guessed) > 0) {
    paste0("cytokit definitions --data <path> --out <file> --force \\\n",
           "  --markers \"", paste(paste0(guessed$column, "=<antibody>"),
                                  collapse = ","), "\"")
  } else NULL,
  if (nrow(guessed) > 0) "```" else NULL,
  if (nrow(guessed) > 0) "" else NULL,
  # A mapping the scientist supplied is provenance. Losing it leaves a column
  # called CD3 with nothing to say which detector it came from.
  if (any(chosen$from != chosen$column)) "## The mapping you supplied" else NULL,
  if (any(chosen$from != chosen$column)) "" else NULL,
  if (any(chosen$from != chosen$column)) {
    "The file names no antibody for these detectors. The names below came from"
  } else NULL,
  if (any(chosen$from != chosen$column)) {
    "the command line and not from the data."
  } else NULL,
  if (any(chosen$from != chosen$column)) "" else NULL,
  if (any(chosen$from != chosen$column)) "| Detector | Antibody |" else NULL,
  if (any(chosen$from != chosen$column)) "|---|---|" else NULL,
  if (any(chosen$from != chosen$column)) {
    mapped_rows <- chosen[chosen$from != chosen$column, , drop = FALSE]
    paste0("| ", mapped_rows$from, " | ", mapped_rows$column, " |")
  } else NULL,
  if (any(chosen$from != chosen$column)) "" else NULL,
  "## The markers in this file",
  "",
  paste0(length(markers), " of the ", length(available), " in the panel."),
  "",
  paste0("- ", markers),
  if (length(markers) < length(available)) "" else NULL,
  if (length(markers) < length(available)) {
    paste0("Not included: ",
           paste(setdiff(available, markers), collapse = ", "))
  } else NULL,
  "",
  "## Example row",
  "",
  "```",
  paste0("Naive B cells,", paste(rep("", length(markers)), collapse = ","),
         ",IgD positive and CD27 negative"),
  "```"
), notes_path)

Say("cytokit definitions")
Say("  panel from  ", basename(files[1]))
Say("  markers     ", length(markers), " of ", length(available),
    " in the panel")
Say("  populations ", length(populations))
if (nrow(guessed) > 0) {
  Say("")
  Say("  WARNING: ", nrow(guessed), " of the ", length(markers),
      " marker columns are detector names, not antibodies.")
  Say("  ", paste(guessed$column, collapse = ", "))
  Say("  A cell type is a claim about an antibody. Ask the scientist which")
  Say("  antibody sits on each detector, then write the file again with")
  Say("  --force and --markers \"", paste(paste0(guessed$column, "=<antibody>"),
                                          collapse = ","), "\".")
  Say("")
}
Say("  wrote       ", DisplayPath(arguments$out))
Say("  notes       ", DisplayPath(notes_path))

# The file is read back with the function that will consume it, so a scaffold
# that does not parse is caught here rather than three steps later.
check <- tryCatch(ReadCellTypeDefinitions(arguments$out),
                  error = function(e) conditionMessage(e))
if (is.character(check)) {
  stop("The file that was written does not parse: ", check)
}
markers_in_file <- setdiff(colnames(check), c("cell_type", "note"))
Say("  parses      yes, ", nrow(check), " population(s), ",
    length(markers_in_file), " marker column(s)")

Say("")
Say("Fill a cell with pos, neg or high, then annotate a clustering with:")
Say("  cytokit annotate --clusters <bundle> --definitions ",
    DisplayPath(arguments$out))
