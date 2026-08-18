#!/usr/bin/env Rscript

# Reproduce the findings of OMIP-039 from the data the authors deposited.
#
# The paper states seven directional differences between CD56dim NKG2C+ and
# CD56dim NKG2C- NK cells, and one composite phenotype it calls prototypic for
# adaptive NK cells. Those claims are in gating/omip39_paper_claims.csv with the
# sentence each came from.
#
# The claims are tested twice.
#
#   Manual route. The deposited FlowJo workspace holds 64 boolean populations
#   under each NKG2C branch, which is Figure 1C. Every marker frequency in
#   Figure 1B is derived by summing those populations, so this route uses the
#   authors' own gates and adds no gating decision.
#
#   Automated route. The openCyto template gates the same file from scratch and
#   the same claims are tested against it.
#
# The first route asks whether the paper's numbers follow from its own gates.
# The
# second asks whether an automated analysis reaches the same conclusions.
#
# Run it in the container:
#   podman run --rm -v "$PWD:/work:z" -w /work everything-flow-cytometry:latest
# \
#     Rscript scripts/04_omip39_reproduce_paper.R
#
# Reference: Hammer Q, Romagnani C. OMIP-039: Detection and analysis of human
# adaptive NKG2C+ natural killer cells. Cytometry A 2017;91A:997-1000.
# PMID 28715616. doi:10.1002/cyto.a.23168.

suppressPackageStartupMessages({
  library(flowCore)
  library(flowStats)
  library(flowWorkspace)
  library(openCyto)
  library(CytoML)
  library(ggplot2)
})

for (file in list.files("R", pattern = "\\.R$", full.names = TRUE)) {
  source(file)
}

kDataDir <- file.path(
  "data", "datasets", "flowrepository", "OMIP-39",
  "FlowRepository_FR-FCM-ZYY6_files"
)
kWorkspace <- file.path(kDataDir, "attachments", "OMIP_Hammer.wsp")
kClaimsPath <- file.path("gating", "omip39_paper_claims.csv")
kTemplatePath <- file.path("gating", "omip39_gating_template.csv")
kOutputDir <- file.path("output", "omip39")

# The six markers the paper combines in Figure 1C. CD7 is not among them,
# because
# the paper states CD7 as an intensity difference and not as a frequency.
kBooleanMarkers <- c("CD2", "CD57", "ILT2", "NKG2A", "NKp30", "Siglec-7")

# The phenotype the paper calls prototypic for adaptive NK cells.
kAdaptivePhenotype <- c(
  CD2 = TRUE, CD57 = TRUE, ILT2 = TRUE,
  NKG2A = FALSE, NKp30 = FALSE, `Siglec-7` = FALSE
)

kManualNKG2Cpos <- "CD56dim NKG2C+"
kManualNKG2Cneg <- "CD56dim NKG2C-"

dir.create(kOutputDir, recursive = TRUE, showWarnings = FALSE)

Log <- function(...) {
  cat(format(Sys.time(), "%H:%M:%S"), "|", ..., "\n")
}

if (!dir.exists(kDataDir)) {
  stop("OMIP-39 is not present. Pull it with:\n",
       "  ./sync.sh pull datasets/flowrepository/OMIP-39")
}

claims <- ReadPaperClaims(kClaimsPath)
Log("Read", nrow(claims), "claims from the paper")

# ---------------------------------------------------------------------------
# Route 1: the authors' own gates
# ---------------------------------------------------------------------------

Log("Importing the deposited FlowJo workspace")
manual_gs <- ImportFlowJoGates(kWorkspace, kDataDir, group = "Sample")
manual_stats <- CollectPopulationStats(manual_gs)

# The boolean populations sit directly under each NKG2C branch and their names
# hold the marker signs. A single marker gate such as "CD2+" is excluded,
# because
# summing it alongside the combinations would count the same events twice.
IsBooleanLeaf <- function(population, parent) {
  leaf <- basename(as.character(population))
  in_branch <- grepl(parent, as.character(population), fixed = TRUE)
  # A combination names every one of the six markers; a single gate names one.
  names_all_markers <- vapply(leaf, function(x) {
    all(vapply(kBooleanMarkers, function(m) grepl(m, x, fixed = TRUE),
        logical(1)))
  }, logical(1))
  in_branch & names_all_markers
}

manual_pos <- manual_stats[IsBooleanLeaf(manual_stats$population,
                           kManualNKG2Cpos), ]
manual_neg <- manual_stats[IsBooleanLeaf(manual_stats$population,
                           kManualNKG2Cneg), ]

Log("Boolean populations found:", nrow(manual_pos), "under", kManualNKG2Cpos,
    "and", nrow(manual_neg), "under", kManualNKG2Cneg)

if (nrow(manual_pos) == 0 || nrow(manual_neg) == 0) {
  stop("No boolean population was found. Check the branch names against ",
       file.path(kOutputDir, "manual_paths.txt"), ".")
}

manual_pos_freq <- MarkerFrequencyFromBooleans(manual_pos, kBooleanMarkers)
manual_neg_freq <- MarkerFrequencyFromBooleans(manual_neg, kBooleanMarkers)

Log("Events in the NKG2C+ branch:", manual_pos_freq$total_events[1],
    "and in the NKG2C- branch:", manual_neg_freq$total_events[1])

manual_result <- TestPaperClaims(
  claims[claims$measure == "percent_positive", ],
  manual_pos_freq,
  manual_neg_freq
)
manual_result$route <- "manual, the authors' gates"

write.csv(manual_result, file.path(kOutputDir, "claims_manual.csv"),
          row.names = FALSE)

cat("\n=== Route 1: the authors' own gates ===\n")
print(manual_result[, c("marker", "direction", "test_percent",
                        "reference_percent", "difference_points", "verdict")],
      digits = 3)

# The composite phenotype of Figure 1C.
manual_phenotype_pos <- FindPhenotype(manual_pos, kAdaptivePhenotype)
manual_phenotype_neg <- FindPhenotype(manual_neg, kAdaptivePhenotype)

phenotype_table <- data.frame(
  branch = c(kManualNKG2Cpos, kManualNKG2Cneg),
  population = c(manual_phenotype_pos$population,
                 manual_phenotype_neg$population),
  count = c(manual_phenotype_pos$count, manual_phenotype_neg$count),
  parent_events = c(manual_phenotype_pos$total_events,
                    manual_phenotype_neg$total_events),
  percent_of_parent = c(manual_phenotype_pos$percent_of_parent,
                        manual_phenotype_neg$percent_of_parent),
  stringsAsFactors = FALSE
)
write.csv(phenotype_table, file.path(kOutputDir,
          "adaptive_phenotype_manual.csv"),
          row.names = FALSE)

cat("\n=== Figure 1C: the prototypic adaptive phenotype ===\n")
cat("NKG2A- NKp30- Siglec-7- ILT2+ CD57+ CD2+\n")
print(phenotype_table, digits = 3)

if (!any(is.na(phenotype_table$percent_of_parent))) {
  enrichment <- phenotype_table$percent_of_parent[1] /
    phenotype_table$percent_of_parent[2]
  Log(sprintf(
    "The phenotype is %.1f times as frequent in NKG2C+ as in NKG2C-",
    enrichment
  ))
}

# Every boolean combination, so a report can show the full distribution.
combination_table <- rbind(
  cbind(branch = "NKG2C+", manual_pos[, c("population", "count",
                                          "percent_of_parent")]),
  cbind(branch = "NKG2C-", manual_neg[, c("population", "count",
                                          "percent_of_parent")])
)
combination_table$population <- basename(as.character(combination_table$population))
write.csv(combination_table,
          file.path(kOutputDir, "boolean_combinations_manual.csv"),
          row.names = FALSE)
Log("Wrote all", nrow(combination_table), "boolean combinations")

# ---------------------------------------------------------------------------
# Route 2: the automated template
# ---------------------------------------------------------------------------

Log("Reading the controls and computing compensation")
control_set <- read.flowSet(
  list.files(kDataDir, pattern = "^Single stainings.*\\.fcs$",
             full.names = TRUE),
  truncate_max_range = FALSE
)
match_table <- MatchControlsToChannels(control_set,
                                       unstained_pattern = "unstained")
match_file <- file.path(kOutputDir, "spillover_match.csv")
WriteMatchFile(match_table, match_file)

computed_spillover <- ComputeSpilloverFromControls(
  control_set, match_file = match_file, method = "median", pregate = TRUE
)

sample_set <- read.flowSet(
  list.files(kDataDir, pattern = "^Samples_.*\\.fcs$", full.names = TRUE),
  truncate_max_range = FALSE
)
compensated_set <- compensate(sample_set, computed_spillover)
transformed_set <- ApplyLogicleTransform(compensated_set)$data

Log("Running the extended template, down to the marker gates")
template <- ReadGatingTemplate(kTemplatePath)
automated_gs <- RunAutomatedGating(transformed_set, template)

automated_paths <- gs_get_pop_paths(automated_gs, path = "auto")
Log("Automated hierarchy holds", length(automated_paths), "populations")
writeLines(automated_paths, file.path(kOutputDir,
           "automated_paths_extended.txt"))

automated_stats <- CollectPopulationStats(automated_gs)
write.csv(automated_stats, file.path(kOutputDir,
          "automated_stats_extended.csv"),
          row.names = FALSE)

# The template gates each marker directly inside each NKG2C branch, so the
# frequency is read rather than derived from combinations. gs_pop_get_stats
# returns the full path, so the lookup matches on the branch and the leaf
# instead
# of building a path by hand.
AutomatedMarkerFrequency <- function(stats, branch, markers) {
  # The template alias for Siglec-7 drops the hyphen, because openCyto reads a
  # hyphen inside an alias as part of a population path.
  alias <- gsub("-", "", markers)
  paths <- as.character(stats$population)

  rows <- lapply(seq_along(markers), function(i) {
    wanted <- paste0("/", branch, "/", alias[i])
    hit <- which(endsWith(paths, wanted))

    if (length(hit) == 0) {
      warning("No automated population ends with '", wanted, "'.")
    }

    data.frame(
      marker = markers[i],
      percent_positive = if (length(hit) == 0) {
        NA_real_
      } else {
        stats$percent_of_parent[hit[1]]
      },
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}

automated_pos_freq <- AutomatedMarkerFrequency(
  automated_stats, "NKG2Cpos", kBooleanMarkers
)
automated_neg_freq <- AutomatedMarkerFrequency(
  automated_stats, "NKG2Cneg", kBooleanMarkers
)

automated_result <- TestPaperClaims(
  claims[claims$measure == "percent_positive", ],
  automated_pos_freq,
  automated_neg_freq
)
automated_result$route <- "automated, openCyto template"

write.csv(automated_result, file.path(kOutputDir, "claims_automated.csv"),
          row.names = FALSE)

cat("\n=== Route 2: the automated template ===\n")
print(automated_result[, c("marker", "direction", "test_percent",
                           "reference_percent", "difference_points",
                           "verdict")],
      digits = 3)

# ---------------------------------------------------------------------------
# Both routes together
# ---------------------------------------------------------------------------

both <- rbind(manual_result, automated_result)
write.csv(both, file.path(kOutputDir, "claims_both_routes.csv"),
          row.names = FALSE)

summary_table <- as.data.frame(table(both$route, both$verdict))
colnames(summary_table) <- c("route", "verdict", "claims")
write.csv(summary_table, file.path(kOutputDir, "claims_summary.csv"),
          row.names = FALSE)

cat("\n=== Summary ===\n")
print(summary_table)

plot_data <- both[!is.na(both$difference_points), ]
plot_data$expected <- ifelse(plot_data$direction == "higher",
                             "higher in NKG2C+",
                             "lower in NKG2C+")

claim_plot <- ggplot(
  plot_data,
  aes(x = reorder(marker, difference_points), y = difference_points,
      fill = verdict)
) +
  geom_col() +
  geom_hline(yintercept = 0, colour = "grey30") +
  facet_wrap(~route) +
  coord_flip() +
  scale_fill_manual(values = c(reproduced = "#2166ac", opposite = "#b2182b",
                               `too small to call` = "grey70")) +
  labs(
    title = "OMIP-039 Figure 1B, tested against the deposited data",
    subtitle = "Difference in percent positive, CD56dim NKG2C+ minus CD56dim NKG2C-",
    x = NULL,
    y = "Difference, percentage points",
    fill = "Verdict"
  ) +
  ThemePublication()

SaveFigure(claim_plot, file.path(kOutputDir, "claims_plot.svg"), width = 11,
  height = 5)
Log("Wrote claims_plot.png")

Log("Done. Output is in", kOutputDir)
