#!/usr/bin/env Rscript

# cytokit gate: run an openCyto template and report what each gate kept.
#
# The template is the scientist's decision, written down. This recipe runs it,
# reports the count and the percentage of the parent for every population, and
# draws the hierarchy. A gate that keeps 100 percent of its parent or that keeps
# almost nothing is a gate that did not work, and the report says so rather than
# leaving it to be noticed three steps later.
#
# Compensation and the transform come first. A gate drawn on uncompensated
# values counts spillover as signal, so the recipe applies the stored matrix
# unless it is told not to.
#
# Called through cli/cytokit, never directly.

suppressPackageStartupMessages({
  library(flowCore)
  library(flowWorkspace)
  library(openCyto)
  library(ggplot2)
})

for (module in c("figures", "io", "panels", "compensation", "transform",
                 "gating", "cytokit")) {
  source(file.path("R", paste0(module, ".R")))
}

arguments <- ParseCytokitArguments(
  commandArgs(trailingOnly = TRUE),
  allowed = c("data", "out", "label", "template", "seed", "cores"),
  required = c("data", "template"),
  flags = c("recursive", "no-compensate", "no-transform", "no-save-gates")
)

seed <- SetCytokitSeed(arguments)
files <- FcsFilesIn(arguments$data, recursive = isTRUE(arguments$recursive))
label <- if (is.null(arguments$label)) basename(arguments$data) else
  arguments$label
out_root <- if (is.null(arguments$out)) kCytokitOutputRoot else arguments$out
cores <- if (is.null(arguments$cores)) 1 else as.integer(arguments$cores)

Say <- function(...) cat(..., "\n", sep = "")

# The template is read before anything expensive runs, because a template that
# does not parse is the most common way for this recipe to fail.
# openCyto prints one line per row while it parses, so the lines are collected
# and reported with the rest.
template_rows <- nrow(utils::read.csv(arguments$template, check.names = FALSE))
read_template <- CollectNotes(tryCatch(ReadGatingTemplate(arguments$template),
                                       error = function(e) e))
template <- read_template$value
if (inherits(template, "error")) {
  stop("The template does not parse: ", conditionMessage(template),
       "\nWrite a valid empty one with:  cytokit template --data ",
       DisplayPath(arguments$data), " --out <file>")
}

bundle <- OpenCytokitBundle("gate", label, out_root)
Say("cytokit gate")
Say("  path     ", DisplayPath(arguments$data))
Say("  files    ", length(files))
Say("  template ", DisplayPath(arguments$template))
Say("  seed     ", seed)
Say("  bundle   ", DisplayPath(bundle), "\n")

read_set <- CollectNotes(
  read.flowSet(files, truncate_max_range = FALSE, transformation = FALSE))
flow_set <- read_set$value

read_panel <- CollectNotes(DescribeFcsPanel(files[1]))
panel <- read_panel$value

# Compensation, then the transform. openCyto cuts on a transformed scale, and a
# cut on a linear scale lands in the wrong place.
compensated <- flow_set
if (!isTRUE(arguments$`no-compensate`)) {
  state <- ReadCompensationState(flow_set[[1]])
  if (identical(state$state, "matrix to apply")) {
    applied <- CollectNotes(tryCatch(ApplyCompensation(flow_set),
                                     error = function(e) e))
    if (inherits(applied$value, "error")) {
      Say("The stored matrix could not be applied: ",
          conditionMessage(applied$value))
      Say("The gates below run on uncompensated values.")
    } else {
      compensated <- applied$value
      Say("Applied the stored ", state$matrix_size, " by ", state$matrix_size,
          " matrix.")
    }
  } else {
    Say("No matrix was applied, because the state is: ", state$state, ".")
    Say("Run cytokit compensate to see whether that costs you a gate.")
  }
} else {
  Say("Compensation was skipped, because --no-compensate was given.")
}

transformed <- compensated
if (!isTRUE(arguments$`no-transform`)) {
  result <- CollectNotes(tryCatch(ApplyLogicleTransform(compensated),
                                  error = function(e) e))
  if (inherits(result$value, "error")) {
    stop("The logicle transform failed: ", conditionMessage(result$value),
         "\nAdd --no-transform when the values are already on a gating scale.")
  }
  transformed <- result$value$data
  Say("Applied the logicle transform.")
} else {
  Say("The transform was skipped, because --no-transform was given.")
}

Say("\nRunning ", template_rows, " template row(s) over ",
    length(files), " file(s)")
gated <- CollectNotes(tryCatch(
  RunAutomatedGating(transformed, template, n_cores = cores),
  error = function(e) e))
if (inherits(gated$value, "error")) {
  stop("The template did not run: ", conditionMessage(gated$value),
       "\nCheck that every dims value in the template names a channel that ",
       "the panel carries.")
}
gating_set <- gated$value

stats <- CollectPopulationStats(gating_set)
WriteBundleTable(bundle, stats, "population_stats.csv")

tree <- CollectGateTree(gating_set)
WriteBundleTable(bundle, tree, "gate_tree.csv")

Say("\nHierarchy, from ", basename(files[1]))
print(tree, row.names = FALSE, digits = 4)

# A gate that keeps everything did not cut, and a gate that keeps almost
# nothing cut in the wrong place. Both are reported, because both look like a
# result until somebody checks.
gates <- tree[!is.na(tree$parent), , drop = FALSE]
kept_all <- gates[gates$percent_of_parent > 99.5, , drop = FALSE]
kept_none <- gates[gates$percent_of_parent < 1, , drop = FALSE]
warnings_table <- rbind(
  if (nrow(kept_all) > 0) {
    data.frame(population = kept_all$population, fault = "keeps every event",
               percent_of_parent = kept_all$percent_of_parent,
               stringsAsFactors = FALSE)
  },
  if (nrow(kept_none) > 0) {
    data.frame(population = kept_none$population, fault = "keeps almost none",
               percent_of_parent = kept_none$percent_of_parent,
               stringsAsFactors = FALSE)
  }
)
if (!is.null(warnings_table) && nrow(warnings_table) > 0) {
  WriteBundleTable(bundle, warnings_table, "gates_to_check.csv")
  Say("\nWARNING: ", nrow(warnings_table), " gate(s) need a second look.")
  print(warnings_table, row.names = FALSE, digits = 4)
  Say("  A gate that keeps every event did not cut. A gate that keeps almost")
  Say("  none cut in the wrong place. Both look like a result in a table.")
} else {
  Say("\nEvery gate kept between 1 and 99.5 percent of its parent.")
}

# The spread across samples is what one template buys over one analyst per
# sample, so it is reported whenever there is more than one file.
if (length(files) > 1) {
  spread <- SummarisePopulationSpread(stats)
  WriteBundleTable(bundle, spread, "population_spread.csv")
  Say("\nSpread across ", length(files), " files")
  print(spread, row.names = FALSE, digits = 4)
}

# The hierarchy is saved so that cluster can work inside one of these gates.
# Without it a scientist has to gate twice to cluster once.
if (!isTRUE(arguments$`no-save-gates`)) {
  gates_path <- file.path(bundle, "gating_set")
  saved <- CollectNotes(tryCatch(
    flowWorkspace::save_gs(gating_set, gates_path), error = function(e) e))
  if (inherits(saved$value, "error")) {
    Say("\nThe hierarchy could not be saved: ", conditionMessage(saved$value))
  } else {
    size <- sum(file.info(list.files(gates_path, recursive = TRUE,
                                     full.names = TRUE))$size, na.rm = TRUE)
    Say("\nSaved the hierarchy to gating_set, ",
        format(structure(size, class = "object_size"), units = "auto"), ".")
    Say("  Cluster inside one of its gates with:")
    Say("    cytokit cluster --gates ", DisplayPath(bundle), " --parent <name>")
  }
}

SaveFigure(PlotGateTree(tree, title = paste("Gate hierarchy,",
                                            basename(files[1]))),
           file.path(bundle, "gate_tree.svg"),
           width = 10, height = 2 + 0.4 * nrow(tree))

notes <- unique(rbind(read_template$notes, read_set$notes,
                      read_panel$notes, gated$notes))
ReportNotes(notes, bundle)

arguments$seed <- seed
CloseCytokitBundle(
  bundle, "gate", arguments, inputs = c(files, arguments$template),
  command = paste("cytokit gate --data", DisplayPath(arguments$data),
                  "--template", DisplayPath(arguments$template))
)

Say("\nWrote population_stats.csv, gate_tree.csv and gate_tree.svg to ",
    DisplayPath(bundle))
