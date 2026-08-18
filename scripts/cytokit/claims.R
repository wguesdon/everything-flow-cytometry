#!/usr/bin/env Rscript

# cytokit claims: put one verdict on every claim, and never invent one.
#
# A claims table is what a paper or a scientist says. A results table is what
# the data says. This recipe joins the two and writes supported, contradicted or
# unresolved beside each claim.
#
# "unresolved" is a real answer. A claim whose measure nobody computed has not
# been contradicted, and recording it as anything else overstates what was
# checked. A run where every claim reads supported and nothing reads unresolved
# is more often a sign that the claims were written from the results.
#
# Called through cli/cytokit, never directly.

suppressPackageStartupMessages({
  library(ggplot2)
})

for (module in c("figures", "reproduce", "cytokit")) {
  source(file.path("R", paste0(module, ".R")))
}

arguments <- ParseCytokitArguments(
  commandArgs(trailingOnly = TRUE),
  allowed = c("claims", "results", "out", "label", "tolerance"),
  required = c("claims", "results"),
  flags = character(0)
)

tolerance <- if (is.null(arguments$tolerance)) 0.05 else
  as.numeric(arguments$tolerance)

Say <- function(...) cat(..., "\n", sep = "")

claims <- utils::read.csv(arguments$claims, check.names = FALSE,
                          stringsAsFactors = FALSE)

results_path <- if (dir.exists(arguments$results)) {
  file.path(arguments$results, "results.csv")
} else {
  arguments$results
}
if (!file.exists(results_path)) {
  stop("No results table is at ", DisplayPath(results_path), ".\n",
       "It needs the columns measure and value, one row per number that a ",
       "claim can be judged against.")
}
results <- utils::read.csv(results_path, check.names = FALSE,
                           stringsAsFactors = FALSE)

label <- if (is.null(arguments$label)) {
  ShortLabel(tools::file_path_sans_ext(basename(arguments$claims)))
} else {
  arguments$label
}
out_root <- if (is.null(arguments$out)) kCytokitOutputRoot else arguments$out
bundle <- OpenCytokitBundle("claims", label, out_root)

Say("cytokit claims")
Say("  claims  ", nrow(claims), " from ", DisplayPath(arguments$claims))
Say("  results ", nrow(results), " measure(s) from ", DisplayPath(results_path))
Say("  bundle  ", DisplayPath(bundle), "\n")

judged <- JudgeClaims(claims, results, tolerance = tolerance)
WriteBundleTable(bundle, judged, "verdicts.csv")

print(judged[, c("claim_id", "measure", "test", "expected", "observed",
                 "verdict")], row.names = FALSE, digits = 4)

tally <- data.frame(
  verdict = kClaimVerdicts,
  claims = vapply(kClaimVerdicts,
                  function(name) sum(judged$verdict == name), integer(1)),
  stringsAsFactors = FALSE
)
WriteBundleTable(bundle, tally, "verdict_tally.csv")

Say("")
print(tally, row.names = FALSE)

unresolved <- judged[judged$verdict == "unresolved", , drop = FALSE]
if (nrow(unresolved) > 0) {
  Say("\n", nrow(unresolved), " claim(s) are unresolved:")
  for (index in seq_len(nrow(unresolved))) {
    Say("  ", unresolved$claim_id[index], ": ", unresolved$reason[index])
  }
  Say("  Unresolved is not contradicted. Compute the measure, or take the")
  Say("  claim out of the table.")
}

contradicted <- judged[judged$verdict == "contradicted", , drop = FALSE]
if (nrow(contradicted) > 0) {
  Say("\n", nrow(contradicted), " claim(s) are contradicted:")
  for (index in seq_len(nrow(contradicted))) {
    Say("  ", contradicted$claim_id[index], ": ",
        contradicted$reason[index])
  }
}
if (nrow(unresolved) == 0 && nrow(contradicted) == 0) {
  Say("\nEvery claim is supported. Check that the claims were written before")
  Say("the results and not from them.")
}

drawing <- ggplot2::ggplot(
  tally, ggplot2::aes(x = factor(.data$verdict, levels = kClaimVerdicts),
                      y = .data$claims, fill = .data$verdict)) +
  ggplot2::geom_col(show.legend = FALSE) +
  ScaleFillPublication() +
  ggplot2::labs(title = "One verdict per claim", x = NULL, y = "Claims") +
  ThemePublication()
SaveFigure(drawing, file.path(bundle, "verdicts.svg"), width = 6, height = 5)

CloseCytokitBundle(
  bundle, "claims", arguments, inputs = c(arguments$claims, results_path),
  command = paste("cytokit claims --claims", DisplayPath(arguments$claims),
                  "--results", DisplayPath(arguments$results))
)

Say("\nWrote verdicts.csv to ", DisplayPath(bundle))
