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
Say("  bundle ", DisplayPath(bundle), "\n")

# One row per file. A study whose files disagree on the parameter count or on
# the compensation state needs that difference resolved before a gate is drawn,
# so the table reports each file rather than a summary.
read_files <- CollectNotes(do.call(rbind, lapply(files, DescribeFcsFile)))
summary_table <- read_files$value
WriteBundleTable(bundle, summary_table, "files.csv")

Say("Files")
print(summary_table, row.names = FALSE)

# The panel is read from the first file, and every other file is checked against
# it. A panel that differs between files cannot be gated by one template.
read_panel <- CollectNotes(DescribeFcsPanel(files[1]))
panel <- read_panel$value
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
  read_others <- CollectNotes(lapply(files[-1], function(path) {
    DescribeFcsPanel(path)$channel
  }))
  other_panels <- read_others$value
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

# A mass cytometry file has no scatter and no spillover. An agent that reads
# "66 detectors" and reaches for a scatter gate gets nothing, so the kind is
# stated before the compensation state rather than left to be inferred.
acquisition <- AcquisitionKind(panel, summary_table$cytometer)
Say("\nAcquisition: ", acquisition$kind, ", because ", acquisition$reason)
if (identical(acquisition$kind, "mass")) {
  Say("  A mass cytometer counts an isotope, so there is no spillover to")
  Say("  compensate and no forward or side scatter to gate on. Skip compensate.")
  Say("  Gate on the DNA channel and on event length instead.")
} else if (!acquisition$has_scatter) {
  Say("  WARNING: this panel carries no scatter detector, so a scatter gate")
  Say("  cannot be drawn on it. Check that the export kept every parameter.")
}

states <- unique(summary_table$compensation_state)
Say("\nCompensation state: ", paste(states, collapse = ", "))
if (any(states == "identity, no matrix supplied")) {
  Say("  An identity matrix means no matrix was supplied. It does not mean the")
  Say("  values are compensated. Run cytokit compensate to measure that.")
}

# One malformed header raises the same warning once per file. The count is the
# fact a reader needs, so it is stated once rather than repeated.
notes <- unique(rbind(read_files$notes, read_panel$notes))
if (exists("read_others")) {
  notes <- unique(rbind(notes, read_others$notes))
}
if (nrow(notes) > 0) {
  notes <- stats::aggregate(times ~ note, data = notes, FUN = sum)
  Say("\nThe FCS reader raised ", sum(notes$times), " note(s) while reading:")
  for (index in seq_len(nrow(notes))) {
    Say("  ", notes$times[index], "x  ", notes$note[index])
  }
  WriteBundleTable(bundle, notes, "reader_notes.csv")
}

CloseCytokitBundle(
  bundle, "inspect", arguments, inputs = files,
  command = paste("cytokit inspect --data", DisplayPath(arguments$data))
)

Say("\nWrote files.csv, panel.csv and markers.txt to ", DisplayPath(bundle))
