#!/usr/bin/env Rscript

# cytokit compare: draw one population's frequency by group, with the test.
#
# The test is chosen from the design and not from the result. Two groups get a
# Wilcoxon rank sum test and more than two get a Kruskal Wallis test, because a
# cell frequency is bounded, skewed, and measured on a handful of samples.
#
# A group that holds one sample carries no spread, so no test runs. The figure
# is still drawn, because a picture of five samples is useful and a p value over
# one sample per group is not.
#
# Called through cli/cytokit, never directly.

suppressPackageStartupMessages({
  library(ggplot2)
})

for (module in c("figures", "comparison", "cytokit")) {
  source(file.path("R", paste0(module, ".R")))
}

arguments <- ParseCytokitArguments(
  commandArgs(trailingOnly = TRUE),
  allowed = c("proportions", "population", "group", "out", "label", "value"),
  required = c("proportions", "group"),
  flags = "all-populations"
)

value_column <- if (is.null(arguments$value)) "percent_of_parent" else
  arguments$value

Say <- function(...) cat(..., "\n", sep = "")

proportions_path <- if (dir.exists(arguments$proportions)) {
  file.path(arguments$proportions, "proportions.csv")
} else {
  arguments$proportions
}
if (!file.exists(proportions_path)) {
  stop("No proportions table is at ", DisplayPath(proportions_path), ".\n",
       "Point --proportions at the bundle that cytokit proportions wrote.")
}
proportions <- utils::read.csv(proportions_path, check.names = FALSE,
                               stringsAsFactors = FALSE)

if (!arguments$group %in% colnames(proportions)) {
  stop("The table has no column called '", arguments$group, "'.\n",
       "It carries: ", paste(colnames(proportions), collapse = ", "))
}

wanted <- if (isTRUE(arguments$`all-populations`)) {
  unique(proportions$population)
} else if (is.null(arguments$population)) {
  stop("Give --population with the population to compare, or ",
       "--all-populations to compare every one.\nThe table holds: ",
       paste(utils::head(unique(proportions$population), 20),
             collapse = ", "))
} else {
  trimws(strsplit(arguments$population, ",")[[1]])
}

label <- if (!is.null(arguments$label)) {
  arguments$label
} else if (dir.exists(arguments$proportions)) {
  ShortLabel(arguments$proportions)
} else {
  ShortLabel(tools::file_path_sans_ext(basename(proportions_path)))
}
out_root <- if (is.null(arguments$out)) kCytokitOutputRoot else arguments$out
bundle <- OpenCytokitBundle("compare", label, out_root)

Say("cytokit compare")
Say("  proportions ", DisplayPath(proportions_path))
Say("  group       ", arguments$group)
Say("  populations ", length(wanted))
Say("  bundle      ", DisplayPath(bundle), "\n")

# A file name has to survive a population path such as /nonDebris/singlets.
SafeName. <- function(name) {
  cleaned <- gsub("[^A-Za-z0-9]+", "_", name)
  cleaned <- gsub("^_+|_+$", "", cleaned)
  if (nzchar(cleaned)) cleaned else "population"
}

results <- lapply(wanted, function(population) {
  outcome <- tryCatch(
    CompareProportions(proportions, population, arguments$group,
                       value = value_column),
    error = function(e) e)
  if (inherits(outcome, "error")) {
    Say("Skipped '", population, "': ", conditionMessage(outcome))
    return(NULL)
  }
  name <- SafeName.(population)
  SaveFigure(outcome$plot, file.path(bundle, paste0(name, ".svg")),
             width = 2 + 1.6 * nrow(outcome$summary), height = 6)
  WriteBundleTable(bundle, outcome$summary, paste0(name, "_by_group.csv"))
  cbind(population = population, outcome$test, stringsAsFactors = FALSE)
})
results <- Filter(Negate(is.null), results)
if (length(results) == 0) {
  stop("No population could be compared.")
}

tests <- do.call(rbind, results)

# Several populations from one table is several tests on one experiment, and a
# p value read one at a time then overstates the evidence.
tested <- !is.na(tests$p_value)
tests$p_value_adjusted <- NA_real_
if (sum(tested) > 1) {
  tests$p_value_adjusted[tested] <- stats::p.adjust(tests$p_value[tested],
                                                    method = "BH")
}
WriteBundleTable(bundle, tests, "tests.csv")

Say("One row per population")
print(tests[, c("population", "test", "groups", "p_value",
                "p_value_adjusted")],
      row.names = FALSE, digits = 4)

skipped <- tests[tests$test == "none", , drop = FALSE]
if (nrow(skipped) > 0) {
  Say("\n", nrow(skipped), " population(s) carry no test.")
  Say("  ", unique(skipped$reason))
  Say("  The figure is still drawn, and it shows one point per sample. Report")
  Say("  it as a picture of the samples and not as a difference.")
}
if (sum(tested) > 1) {
  Say("\n", sum(tested), " tests ran on one experiment, so p_value_adjusted")
  Say("  holds the Benjamini and Hochberg value. Read that column and not the")
  Say("  raw one.")
}

CloseCytokitBundle(
  bundle, "compare", arguments, inputs = proportions_path,
  command = paste("cytokit compare --proportions",
                  DisplayPath(arguments$proportions), "--group",
                  arguments$group)
)

Say("\nWrote tests.csv and ", length(results), " figure(s) to ",
    DisplayPath(bundle))
