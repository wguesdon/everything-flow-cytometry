#!/usr/bin/env Rscript

# cytokit inspect: read a panel and report what is in it.
#
# Every later recipe reads what this one prints. It reports the markers so that
# an agent can draft a gating template and a cell type definitions table, and it
# reports the compensation state by measuring the stored matrix rather than by
# trusting the keyword. On FR-FCM-ZYRN an identity matrix with
# APPLY COMPENSATION = TRUE meant that no matrix was supplied, and not that the
# values were compensated.
#
# Called through cli/cytokit, never directly.

suppressPackageStartupMessages({
  library(flowCore)
})

for (module in c("io", "panels", "cytokit")) {
  source(file.path("R", paste0(module, ".R")))
}

arguments <- ParseCytokitArguments(
  commandArgs(trailingOnly = TRUE),
  allowed = c("data", "out", "label"),
  required = "data",
  flags = "recursive"
)

files <- FcsFilesIn(arguments$data, recursive = isTRUE(arguments$recursive))
label <- if (is.null(arguments$label)) basename(arguments$data) else
  arguments$label
out_root <- if (is.null(arguments$out)) kCytokitOutputRoot else arguments$out
bundle <- OpenCytokitBundle("inspect", label, out_root)

Say <- function(...) cat(..., "\n", sep = "")

Say("cytokit inspect")
Say("  path   ", DisplayPath(arguments$data))
Say("  files  ", length(files))
Say("  bundle ", bundle, "\n")

# One row per file. A study whose files disagree on the parameter count or on
# the compensation state needs that difference resolved before a gate is drawn,
# so the table reports each file rather than a summary.
summary_table <- do.call(rbind, lapply(files, DescribeFcsFile))
WriteBundleTable(bundle, summary_table, "files.csv")

Say("Files")
print(summary_table, row.names = FALSE)

# The panel is read from the first file, and every other file is checked against
# it. A panel that differs between files cannot be gated by one template.
panel <- DescribeFcsPanel(files[1])
panel$file <- basename(files[1])
WriteBundleTable(bundle, panel, "panel.csv")

Say("\nPanel, from ", basename(files[1]))
print(panel[, c("channel", "marker", "kind", "range")], row.names = FALSE)

naming <- PanelNamingState(panel)
markers <- PanelMarkers(panel)
writeLines(markers, file.path(bundle, "markers.txt"))
WriteBundleTable(bundle, naming, "naming.csv")

Say("\nMarker names: ", naming$state,
    " (", naming$named, " named, ", naming$unnamed, " unnamed)")
Say("  ", paste(markers, collapse = ", "))
if (naming$unnamed > 0) {
  Say("")
  Say("  The names above fall back to the detector for ", naming$unnamed,
      " detector(s),")
  Say("  because the file leaves $PnS empty. A detector name is not an")
  Say("  antibody. Supply the mapping before you read a cell type label.")
}

if (length(files) > 1) {
  other_panels <- lapply(files[-1], function(path) {
    DescribeFcsPanel(path)$channel
  })
  differs <- vapply(other_panels,
                    function(channels) !identical(channels, panel$channel),
                    logical(1))
  if (any(differs)) {
    Say("\nWARNING: ", sum(differs), " of ", length(files),
        " files carry a different set of detectors than the first.")
    Say("  One gating template cannot cover them. Split the study by panel.")
    Say("  ", paste(basename(files[-1])[differs], collapse = ", "))
  } else {
    Say("\nEvery file carries the same ", nrow(panel), " detectors.")
  }
}

states <- unique(summary_table$compensation_state)
Say("\nCompensation state: ", paste(states, collapse = ", "))
if (any(states == "identity, no matrix supplied")) {
  Say("  An identity matrix means no matrix was supplied. It does not mean the")
  Say("  values are compensated. Run cytokit compensate to measure that.")
}

CloseCytokitBundle(
  bundle, "inspect", arguments, inputs = files,
  command = paste("cytokit inspect --data", DisplayPath(arguments$data))
)

Say("\nWrote files.csv, panel.csv and markers.txt to ", bundle)
