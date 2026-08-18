#!/usr/bin/env Rscript

# cytokit proportions: join a count table to what was done to each sample.
#
# A count is not a result. The result is the share of a parent that a population
# holds in each sample, read beside the treatment of that sample. An FCS file
# carries no treatment, so the metadata table comes from the scientist.
#
# The recipe reads the counts that cytokit gate wrote. A sample that the
# metadata does not name is kept and reported, because dropping it silently is
# how a group loses a replicate.
#
# Called through cli/cytokit, never directly.

for (module in c("figures", "gating", "comparison", "cytokit")) {
  source(file.path("R", paste0(module, ".R")))
}

arguments <- ParseCytokitArguments(
  commandArgs(trailingOnly = TRUE),
  allowed = c("counts", "metadata", "out", "label", "sample-column"),
  required = c("counts", "metadata"),
  flags = character(0)
)

sample_column <- if (is.null(arguments$`sample-column`)) "sample" else
  arguments$`sample-column`

Say <- function(...) cat(..., "\n", sep = "")

counts_path <- if (dir.exists(arguments$counts)) {
  file.path(arguments$counts, "population_stats.csv")
} else {
  arguments$counts
}
if (!file.exists(counts_path)) {
  stop("No count table is at ", DisplayPath(counts_path), ".\n",
       "Point --counts at the bundle that cytokit gate wrote, or at a CSV ",
       "with the columns sample, population and percent_of_parent.")
}
stats <- utils::read.csv(counts_path, check.names = FALSE,
                         stringsAsFactors = FALSE)
metadata <- utils::read.csv(arguments$metadata, check.names = FALSE,
                            stringsAsFactors = FALSE)

label <- if (!is.null(arguments$label)) {
  arguments$label
} else if (dir.exists(arguments$counts)) {
  ShortLabel(arguments$counts)
} else {
  ShortLabel(tools::file_path_sans_ext(basename(counts_path)))
}
out_root <- if (is.null(arguments$out)) kCytokitOutputRoot else arguments$out
bundle <- OpenCytokitBundle("proportions", label, out_root)

Say("cytokit proportions")
Say("  counts     ", DisplayPath(counts_path))
Say("  metadata   ", DisplayPath(arguments$metadata))
Say("  bundle     ", DisplayPath(bundle), "\n")

joined <- PopulationProportions(stats, metadata, sample_column = sample_column)
WriteBundleTable(bundle, joined, "proportions.csv")

samples <- unique(joined$sample)
unmatched <- unique(joined$sample[!joined$metadata_found])
Say(length(samples), " sample(s) and ",
    length(unique(joined$population)), " population(s)")
if (length(unmatched) > 0) {
  Say("\nWARNING: the metadata names no row for ", length(unmatched),
      " sample(s):")
  Say("  ", paste(unmatched, collapse = ", "))
  Say("  They are kept, and every metadata column is empty for them. A")
  Say("  comparison drops them, so a group can lose a replicate here.")
  Say("  Check that the ", sample_column, " column matches the file names.")
}

described <- setdiff(colnames(metadata), sample_column)
Say("\nWhat the metadata describes")
for (column in described) {
  values <- unique(joined[[column]][joined$metadata_found])
  Say("  ", column, ": ", length(values), " value(s): ",
      paste(utils::head(values, 8), collapse = ", "))
}

# A group that holds one sample cannot carry a test, and that is a fact about
# the design that is better known now than after the box plot is drawn.
design <- do.call(rbind, lapply(described, function(column) {
  per_sample <- unique(joined[joined$metadata_found,
                              c("sample", column), drop = FALSE])
  counts_per_group <- table(per_sample[[column]])
  data.frame(column = column, groups = length(counts_per_group),
             smallest_group = if (length(counts_per_group) == 0) 0L else
               min(as.integer(counts_per_group)),
             testable = length(counts_per_group) >= 2 &&
               all(counts_per_group >= 2),
             stringsAsFactors = FALSE)
}))
if (!is.null(design) && nrow(design) > 0) {
  WriteBundleTable(bundle, design, "design.csv")
  Say("\nWhich columns can carry a test")
  print(design, row.names = FALSE)
  untestable <- design[!design$testable, , drop = FALSE]
  if (nrow(untestable) > 0) {
    Say("  A test needs two groups and two samples in each. These columns do")
    Say("  not have that: ", paste(untestable$column, collapse = ", "), ".")
    Say("  A box plot over one sample per group draws one point per box.")
  }
}

Say("\nCompare one population between groups with:")
Say("  cytokit compare --proportions ", DisplayPath(bundle),
    " --population <name> --group <column>")

CloseCytokitBundle(
  bundle, "proportions", arguments,
  inputs = c(counts_path, arguments$metadata),
  command = paste("cytokit proportions --counts", DisplayPath(arguments$counts),
                  "--metadata", DisplayPath(arguments$metadata))
)

Say("\nWrote proportions.csv to ", DisplayPath(bundle))
