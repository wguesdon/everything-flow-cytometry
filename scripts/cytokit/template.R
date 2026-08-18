#!/usr/bin/env Rscript

# cytokit template: write an empty openCyto gating template for this panel.
#
# A gating template is per panel and nobody arrives with one. This writes the
# correct header, lists the panel's detectors and markers as comments, and adds
# the two scatter steps that almost every panel starts with. An agent fills the
# rest with the scientist, then runs cytokit gate, which validates it.
#
# The file is a starting point and not a gate. Every threshold in it still has
# to be chosen and checked against the data.
#
# Called through cli/cytokit, never directly.

suppressPackageStartupMessages({
  library(flowCore)
})

for (module in c("io", "gating", "panels", "cytokit")) {
  source(file.path("R", paste0(module, ".R")))
}

arguments <- ParseCytokitArguments(
  commandArgs(trailingOnly = TRUE),
  allowed = c("data", "out"),
  required = c("data", "out"),
  flags = c("force", "recursive")
)

if (file.exists(arguments$out) && !isTRUE(arguments$force)) {
  stop("The file already exists: ", arguments$out,
       "\nAdd --force to overwrite it. A template holds decisions that took ",
       "work, so it is not overwritten by default.")
}

files <- FcsFilesIn(arguments$data,
                    recursive = isTRUE(arguments$recursive))
panel <- DescribeFcsPanel(files[1])
markers <- PanelMarkers(panel)
scatter <- panel$channel[panel$kind == "scatter"]

Say <- function(...) cat(..., "\n", sep = "")

# The header openCyto reads. ReadGatingTemplate checks these names, so they are
# written from one place rather than typed twice.
header <- c("alias", "pop", "parent", "dims", "gating_method", "gating_args",
            "collapseDataForGating", "groupBy", "preprocessing_method",
            "preprocessing_args")

Row <- function(alias, pop, parent, dims, method, args = "") {
  values <- c(alias, pop, parent, dims, method, args, "", "", "", "")
  paste(ifelse(grepl(",", values), paste0("\"", values, "\""), values),
        collapse = ",")
}

rows <- character(0)

# Two steps that nearly every panel begins with, so the file is runnable and a
# scientist can see the shape before adding their own populations.
forward_area <- grep("^FSC-A$", scatter, value = TRUE, ignore.case = TRUE)
forward_height <- grep("^FSC-H$", scatter, value = TRUE, ignore.case = TRUE)
side_area <- grep("^SSC-A$", scatter, value = TRUE, ignore.case = TRUE)

if (length(forward_area) == 1 && length(side_area) == 1) {
  rows <- c(rows, Row("nonDebris", "+", "root",
                      paste(forward_area, side_area, sep = ","),
                      "flowClust", "K=2"))
  if (length(forward_height) == 1) {
    rows <- c(rows, Row("singlets", "+", "nonDebris",
                        paste(forward_area, forward_height, sep = ","),
                        "singletGate"))
  }
} else {
  Say("Note: no FSC-A and SSC-A pair was found, so no scatter step was added.")
}

# The CSV carries no comment lines. ReadGatingTemplate calls read.csv, whose
# comment character is empty by default, so a leading comment becomes the header
# and the template stops parsing. The guidance goes in a sibling file instead.
dir.create(dirname(arguments$out), recursive = TRUE, showWarnings = FALSE)
writeLines(c(paste(header, collapse = ","), rows), arguments$out)

notes_path <- sub("\\.csv$", "_notes.md", arguments$out)
if (identical(notes_path, arguments$out)) {
  notes_path <- paste0(arguments$out, "_notes.md")
}
writeLines(c(
  paste0("# Gating template for ", basename(files[1])),
  "",
  paste0("Written by `cytokit template` beside `", basename(arguments$out),
         "`."),
  "",
  "One row is one population. `parent` names the row above it in the",
  "hierarchy, `dims` names the channel or the marker to cut on, and",
  "`gating_method` names the rule. `mindensity` finds a density minimum,",
  "`quantileGate` cuts at a quantile, and `flowClust` fits a cluster.",
  "",
  "A threshold that a rule fits is not a threshold that is right. Check each",
  "one against the data before you trust the count it produces.",
  "",
  "## The panel",
  "",
  paste0("Scatter detectors: ",
         if (length(scatter) > 0) paste(scatter, collapse = ", ") else "none"),
  "",
  paste0("Markers, for the `dims` column (", length(markers), "):"),
  "",
  paste0("- ", markers),
  "",
  "## Example rows",
  "",
  "```",
  "live,-,singlets,<viability marker>,mindensity,,,,,",
  "CD3,+,live,CD3,mindensity,,,,,",
  "```"
), notes_path)

# The file is read back with the function that will consume it, so a scaffold
# that does not parse is caught here rather than when a gate is run. openCyto
# rejects a template with no rows, which is why the scatter steps are written.
check <- tryCatch({
  ReadGatingTemplate(arguments$out)
  "yes"
}, error = function(e) conditionMessage(e))

Say("cytokit template")
Say("  panel from ", basename(files[1]))
Say("  markers    ", length(markers))
Say("  wrote      ", DisplayPath(arguments$out), " with ", length(rows),
    " starting row(s)")
Say("  parses     ", check)
Say("  notes      ", DisplayPath(notes_path))
Say("")
Say("Markers available for the dims column:")
Say("  ", paste(markers, collapse = ", "))
Say("")
Say("Add one row per population, then check it with:")
Say("  cytokit gate --data ", DisplayPath(arguments$data),
    " --template ", DisplayPath(arguments$out))
