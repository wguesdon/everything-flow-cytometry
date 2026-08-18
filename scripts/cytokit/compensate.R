#!/usr/bin/env Rscript

# cytokit compensate: read the spillover matrix, and measure whether it worked.
#
# The question a scientist arrives with is "is my data compensated". The
# keyword does not answer it. On FR-FCM-ZYRN, FR-FCM-ZZCA and one more deposit
# an identity matrix sits beside APPLY COMPENSATION = TRUE, which means that no
# matrix was supplied and not that the values are compensated.
#
# The measurement is the correlation between every pair of markers, before and
# after the stored matrix is applied. Spillover shows as a correlation between
# two detectors that measure two different antibodies, and compensation should
# drive it down. A pair that is driven far negative is over-compensated, which
# is why the count uses the absolute value.
#
# With --controls the recipe also computes a matrix from the single stain
# controls and compares it with the stored one.
#
# Called through cli/cytokit, never directly.

suppressPackageStartupMessages({
  library(flowCore)
  library(ggplot2)
})

for (module in c("figures", "io", "panels", "compensation", "spillover_compute",
                 "cytokit")) {
  source(file.path("R", paste0(module, ".R")))
}

arguments <- ParseCytokitArguments(
  commandArgs(trailingOnly = TRUE),
  allowed = c("data", "out", "label", "controls", "cofactor", "unstained",
              "seed"),
  required = "data",
  flags = "recursive"
)

seed <- SetCytokitSeed(arguments)

files <- FcsFilesIn(arguments$data, recursive = isTRUE(arguments$recursive))
label <- if (is.null(arguments$label)) basename(arguments$data) else
  arguments$label
out_root <- if (is.null(arguments$out)) kCytokitOutputRoot else arguments$out
cofactor <- if (is.null(arguments$cofactor)) 150 else
  as.numeric(arguments$cofactor)

Say <- function(...) cat(..., "\n", sep = "")

read_panel <- CollectNotes(DescribeFcsPanel(files[1]))
panel <- read_panel$value
read_files <- CollectNotes(do.call(rbind, lapply(files, DescribeFcsFile)))
summary_table <- read_files$value

# Both refusals come before the bundle is opened. A recipe that stops after it
# has created a folder leaves an empty bundle that looks like a result.
#
# A mass cytometer counts an isotope. There is no spillover to compute, and a
# matrix computed from single stains would be an answer to a question nobody
# asked.
acquisition <- AcquisitionKind(panel, summary_table$cytometer)
if (identical(acquisition$kind, "mass")) {
  stop("This is a mass cytometry run, because ", acquisition$reason, ".\n",
       "A mass cytometer counts an isotope, so there is no spillover to ",
       "compensate.\nSignal spillover on a mass cytometer is corrected at ",
       "acquisition, not by a matrix here.")
}

markers <- panel$channel[panel$kind %in% c("stain", "unnamed")]
if (length(markers) < 2) {
  stop("This panel carries ", length(markers), " marker channel(s). ",
       "A correlation needs two.")
}

bundle <- OpenCytokitBundle("compensate", label, out_root)
Say("cytokit compensate")
Say("  path   ", DisplayPath(arguments$data))
Say("  files  ", length(files))
Say("  seed   ", seed)
Say("  bundle ", DisplayPath(bundle), "\n")

# The correlation reads every event of one file. Saying which file, and that
# the others were not read, keeps the number from reading as a study wide one.
if (length(files) > 1) {
  Say("The correlation below is measured on ", basename(files[1]), " only.")
  Say("The other ", length(files) - 1, " file(s) are read for their header.")
  Say("A matrix that suits one file can still fail on another, so measure a")
  Say("second file when the run spans more than one day or one instrument.\n")
}

# The stored matrix, read from the first file and checked against the rest.
states <- unique(summary_table$compensation_state)
Say("Compensation state: ", paste(states, collapse = ", "))
if (length(states) > 1) {
  Say("  WARNING: the files disagree. One matrix cannot cover them, so split")
  Say("  the study by state before you go on.")
}

frame <- read.FCS(files[1], truncate_max_range = FALSE, transformation = FALSE)
stored <- tryCatch(ExtractSpillover(frame), error = function(e) NULL)
has_stored <- !is.null(stored) && is.matrix(stored) &&
  max(abs(stored - diag(nrow(stored)))) > 0

if (has_stored) {
  Say("  The file carries a ", nrow(stored), " by ", ncol(stored),
      " matrix with a largest off diagonal spill of ",
      sprintf("%.1f percent", 100 * max(abs(stored - diag(nrow(stored))))), ".")
  WriteBundleTable(bundle, as.data.frame(stored), "stored_matrix.csv")
  WriteBundleTable(bundle, SummariseSpillover(stored, top = 20),
                   "stored_matrix_top_spills.csv")
  SaveFigure(PlotSpilloverHeatmap(stored, title = "Stored spillover matrix"),
             file.path(bundle, "stored_matrix.svg"), width = 9, height = 8)
} else {
  Say("  No matrix was supplied. An identity matrix beside APPLY COMPENSATION")
  Say("  = TRUE means the same thing, and it does not mean the values are")
  Say("  compensated.")
}

# The measurement. The same events are correlated twice, so the number to read
# is the change and not the value on its own.
Correlate. <- function(matrix_values) {
  events <- ArcsinhTransform(matrix_values, markers, cofactor = cofactor)
  events <- events[, markers, drop = FALSE]
  named <- PanelMarkers(panel)
  if (length(named) == length(markers)) {
    colnames(events) <- named
  }
  MarkerCorrelation(events, seed = seed, threshold = 0.5)
}

before <- Correlate.(exprs(frame))
rows <- list(cbind(state = "as stored", before$summary))

if (has_stored) {
  compensated <- tryCatch(ApplyCompensation(frame), error = function(e) e)
  if (inherits(compensated, "error")) {
    Say("\nThe stored matrix could not be applied: ",
        conditionMessage(compensated))
  } else {
    after <- Correlate.(exprs(compensated))
    rows[[length(rows) + 1]] <- cbind(state = "stored matrix applied",
                                      after$summary)
    WriteBundleTable(bundle, after$pairs, "pairs_after.csv")
  }
}
WriteBundleTable(bundle, before$pairs, "pairs_before.csv")

# A matrix computed from the single stain controls, when they were given.
if (!is.null(arguments$controls)) {
  control_files <- FcsFilesIn(arguments$controls, recursive = TRUE)
  Say("\nControls: ", length(control_files), " file(s) in ",
      DisplayPath(arguments$controls))
  read_controls <- CollectNotes(
    read.flowSet(control_files, truncate_max_range = FALSE,
                 transformation = FALSE))
  controls <- read_controls$value
  pattern <- if (is.null(arguments$unstained)) "unstain" else
    arguments$unstained
  # The matcher reports one line per channel. On a 27 colour panel that is 27
  # lines that hide the result, so they are collected like any other note.
  read_match <- CollectNotes(
    MatchControlsToChannels(controls, unstained_pattern = pattern))
  match_table <- read_match$value
  WriteBundleTable(bundle, match_table, "control_match.csv")

  unmatched <- sum(match_table$matched_by == "none")
  Say("  matched ", sum(match_table$matched_by != "none"), " of ",
      nrow(match_table), " control file(s) to a channel")
  if (unmatched > 0) {
    Say("  ", unmatched, " file(s) matched nothing. Rename them after the")
    Say("  antibody, or the matrix will be computed from a smaller panel.")
  }

  match_file <- file.path(bundle, "control_match_used.csv")
  WriteMatchFile(match_table, match_file)
  read_computed <- CollectNotes(tryCatch(
    ComputeSpilloverFromControls(controls, match_file),
    error = function(e) e))
  computed <- read_computed$value
  if (inherits(computed, "error")) {
    Say("  The matrix could not be computed: ", conditionMessage(computed))
  } else {
    WriteBundleTable(bundle, as.data.frame(computed), "computed_matrix.csv")
    SaveFigure(PlotSpilloverHeatmap(computed,
                                    title = "Matrix from the controls"),
               file.path(bundle, "computed_matrix.svg"), width = 9, height = 8)
    Say("  computed a ", nrow(computed), " by ", ncol(computed), " matrix")

    if (has_stored) {
      comparison <- CompareSpilloverMatrices(computed, stored)
      WriteBundleTable(bundle, comparison, "matrix_comparison.csv")
      Say("  the two matrices are compared in matrix_comparison.csv")
    }

    applied <- tryCatch(compensate(frame, computed), error = function(e) e)
    if (!inherits(applied, "error")) {
      control_result <- Correlate.(exprs(applied))
      rows[[length(rows) + 1]] <- cbind(state = "control matrix applied",
                                        control_result$summary)
    }
  }
}

correlation_table <- do.call(rbind, rows)
WriteBundleTable(bundle, correlation_table, "marker_correlation.csv")

Say("\nMarker correlation, over ", length(markers), " channels and ",
    correlation_table$pairs[1], " pairs")
print(correlation_table[, c("state", "median_absolute_r",
                            "pairs_above_threshold", "max_absolute_r",
                            "most_correlated")], row.names = FALSE)
Say("")
Say("  The events are not gated, so read the change between the rows and not")
Say("  the value on its own. Compensation should lower the median. A pair")
Say("  driven far negative is over-compensated.")

if (nrow(correlation_table) > 1) {
  change <- correlation_table$median_absolute_r[1] -
    correlation_table$median_absolute_r[nrow(correlation_table)]
  verdict <- if (change > 0.05) {
    "the matrix lowers the correlation, so it is doing work"
  } else if (change < -0.05) {
    "the matrix raises the correlation, which is a fault in the matrix"
  } else {
    "the matrix changes almost nothing, so the values may already be compensated"
  }
  Say("  Verdict: ", verdict, ".")
} else {
  # One row means there was no matrix to apply, so the recipe measured the
  # state and cannot improve it. Saying what to do next is the useful part.
  worst <- before$pairs[1, ]
  Say("  There is one row because no matrix was supplied, so nothing could be")
  Say("  applied. The strongest pair is ", worst$a, " against ", worst$b,
      " at r = ", sprintf("%.2f", worst$r), ".")
  Say("  Two detectors that measure two different antibodies should not")
  Say("  correlate that far. Give the single stain controls with --controls,")
  Say("  or send the matrix that the acquisition software exported.")
}

notes <- unique(rbind(read_files$notes, read_panel$notes))
for (name in c("read_controls", "read_match", "read_computed")) {
  if (exists(name)) {
    notes <- unique(rbind(notes, get(name)$notes))
  }
}
ReportNotes(notes, bundle)

# The computed matrix depends on the control files as much as on the samples,
# so both are checksummed. A manifest that names only the samples cannot say
# which controls produced the matrix beside it.
inputs <- files
if (exists("control_files")) {
  inputs <- c(inputs, control_files)
}
command <- paste("cytokit compensate --data", DisplayPath(arguments$data))
if (!is.null(arguments$controls)) {
  command <- paste(command, "--controls", DisplayPath(arguments$controls))
}
arguments$seed <- seed
CloseCytokitBundle(bundle, "compensate", arguments, inputs = inputs,
                   command = command)

Say("\nWrote marker_correlation.csv and the matrices to ", DisplayPath(bundle))
