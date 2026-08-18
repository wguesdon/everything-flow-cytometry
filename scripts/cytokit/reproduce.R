#!/usr/bin/env Rscript

# cytokit reproduce: list the analyses in this repository, and run one.
#
# This is the secondary use of the tool. The first use is the scientist's own
# data. This one rebuilds an analysis that is already here, which is useful to
# check that the container still produces the same numbers and to read a worked
# example before writing your own template.
#
# Every analysis names its script, its deposit and its report. The script writes
# into output/ and the report reads from there.
#
# Called through cli/cytokit, never directly.

for (module in c("cytokit")) {
  source(file.path("R", paste0(module, ".R")))
}

arguments <- ParseCytokitArguments(
  commandArgs(trailingOnly = TRUE),
  allowed = c("analysis", "out", "label"),
  required = character(0),
  flags = "list"
)

Say <- function(...) cat(..., "\n", sep = "")

map_path <- file.path("docs", "analyses.csv")
if (!file.exists(map_path)) {
  stop("The analysis map is missing: ", map_path)
}
analyses <- utils::read.csv(map_path, check.names = FALSE,
                            stringsAsFactors = FALSE)

# Nothing to run means list, because a caller who does not know the names
# cannot ask for one.
if (is.null(arguments$analysis) || isTRUE(arguments$list)) {
  Say("The analyses in this repository\n")
  present <- file.exists(analyses$script)
  for (index in seq_len(nrow(analyses))) {
    Say(sprintf("%-32s %s", analyses$analysis[index],
                if (present[index]) analyses$script[index] else
                  paste(analyses$script[index], "(the script is missing)")))
    Say("  ", analyses$what_it_shows[index])
    Say("  deposit ", analyses$deposit[index], ", report ",
        analyses$report[index])
    Say("")
  }
  Say("Run one with:  cytokit reproduce --analysis <name>")
  Say("Each one reads from data/ and writes into output/.")
  quit(save = "no", status = 0)
}

wanted <- analyses[analyses$analysis == arguments$analysis, , drop = FALSE]
if (nrow(wanted) == 0) {
  stop("No analysis is called '", arguments$analysis, "'.\n",
       "The names are: ", paste(analyses$analysis, collapse = ", "))
}
script <- wanted$script[1]
if (!file.exists(script)) {
  stop("The script is missing: ", script)
}

label <- if (is.null(arguments$label)) wanted$analysis[1] else arguments$label
out_root <- if (is.null(arguments$out)) kCytokitOutputRoot else arguments$out
bundle <- OpenCytokitBundle("reproduce", label, out_root)

Say("cytokit reproduce")
Say("  analysis ", wanted$analysis[1])
Say("  script   ", script)
Say("  deposit  ", wanted$deposit[1])
Say("  report   ", wanted$report[1])
Say("  bundle   ", DisplayPath(bundle), "\n")
Say(wanted$what_it_shows[1], "\n")

# The script is sourced in its own environment rather than called through
# Rscript, so that one container run does the work and the log lands here.
log_path <- file.path(bundle, "run_log.txt")
Say("Running the script. The log goes to run_log.txt.\n")
# The sinks are unwound on the way out whatever happens, because a session left
# with an open sink prints nothing and looks hung.
connection <- file(log_path, open = "wt")
on.exit({
  while (sink.number(type = "message") > 2) sink(type = "message")
  while (sink.number() > 0) sink()
  close(connection)
}, add = TRUE)
sink(connection, split = TRUE)
sink(connection, type = "message")
status <- tryCatch({
  source(script, echo = FALSE, local = new.env())
  "finished"
}, error = function(e) paste("failed:", conditionMessage(e)))
sink(type = "message")
sink()

Say("\nThe script ", status, ".")

outcome <- data.frame(analysis = wanted$analysis[1], script = script,
                      deposit = wanted$deposit[1], report = wanted$report[1],
                      status = status, stringsAsFactors = FALSE)
WriteBundleTable(bundle, outcome, "outcome.csv")

if (identical(status, "finished")) {
  Say("Render its report with:")
  Say("  quarto render ", wanted$report[1])
} else {
  Say("The log is in ", DisplayPath(log_path), ".")
}

CloseCytokitBundle(
  bundle, "reproduce", arguments, inputs = c(map_path, script),
  command = paste("cytokit reproduce --analysis", wanted$analysis[1])
)

Say("\nWrote outcome.csv and run_log.txt to ", DisplayPath(bundle))
if (!identical(status, "finished")) {
  quit(save = "no", status = 1)
}
