# Multicentre harmonisation, FlowRepository FR-FCM-Z282.
#
# Three arms measure the same 234 files.
#
#   local      13 operators, each one analysing the files that operator acquired
#   central    one reference operator, analysing all 234 files in Kaluza
#   automated  this repository, analysing all 234 files with one rule
#
# The first two arms are read from the spreadsheet the authors attached to the
# deposit. The third is computed here and it exists nowhere else.
#
# Run it in the container:
#   podman run --rm -v "$PWD:/work:z" -w /work everything-flow-cytometry:latest \
#     Rscript scripts/09_z282_harmonisation.R

suppressPackageStartupMessages({
  library(flowCore)
  library(ggplot2)
  library(lme4)
  library(withr)
})

source(file.path("R", "harmonisation.R"))
source(file.path("R", "spectral.R"))
source(file.path("R", "naive_memory.R"))

kDeposit <- file.path(
  "data", "datasets", "flowrepository", "FR-FCM-Z282",
  "FlowRepository_FR-FCM-Z282_files"
)
kGatingDir <- "gating"
kOutputDir <- file.path("output", "z282")

# The fifteen frequencies the study reports, in the order the paper counts them.
kParameters <- c(
  "CD3", "CD4", "CD8",
  "CD3_N", "CD4_N", "CD8_N",
  "CD3_CM", "CD4_CM", "CD8_CM",
  "CD3_EM", "CD4_EM", "CD8_EM",
  "CD3_TD", "CD4_TD", "CD8_TD"
)

# The acceptability thresholds the paper sets for each indicator.
kCvThreshold <- 0.20
kBiasThreshold <- 0.20
kIccThreshold <- 0.75
kZUnsatisfactory <- 3

dir.create(kOutputDir, recursive = TRUE, showWarnings = FALSE)

Write <- function(x, name) {
  utils::write.csv(x, file.path(kOutputDir, name), row.names = FALSE)
  invisible(x)
}

Say <- function(...) {
  cat(..., "\n", sep = "")
}

# ---------------------------------------------------------------------------
# Part 1. The design, read from two independent sources.
# ---------------------------------------------------------------------------

Say("Part 1: the design")

samples <- utils::read.csv(
  file.path(kGatingDir, "z282_sample_metadata.csv"), stringsAsFactors = FALSE
)
results <- utils::read.csv(
  file.path(kGatingDir, "z282_operator_results.csv"), stringsAsFactors = FALSE
)
panel <- ReadZ282Panel(file.path(kGatingDir, "z282_panel.csv"))
claims <- utils::read.csv(
  file.path(kGatingDir, "z282_paper_claims.csv"), stringsAsFactors = FALSE
)

ObservationKey <- function(frame) {
  paste(frame$operator, frame$material, frame$vial, frame$round, sep = "|")
}

samples$key <- ObservationKey(samples)
results$key <- ObservationKey(results)

stopifnot(
  nrow(samples) == 234,
  length(setdiff(samples$key, results$key)) == 0,
  length(setdiff(results$key, samples$key)) == 0
)

design <- data.frame(
  Measure = c(
    "Sample files", "Observations in the spreadsheet", "Operators", "Centres",
    "Instruments", "Instrument models", "Vials per material", "Rounds",
    "Parameters per observation", "Excluded observations"
  ),
  Value = c(
    nrow(samples),
    length(unique(results$key)),
    length(unique(samples$operator)),
    length(unique(samples$centre)),
    length(unique(samples$instrument)),
    length(unique(samples$model)),
    length(unique(samples$vial)),
    length(unique(samples$round)),
    length(unique(results$parameter)),
    length(unique(results$key[results$excluded == "Yes"]))
  ),
  stringsAsFactors = FALSE
)
Write(design, "design.csv")
print(design)

operators <- unique(samples[, c("operator", "centre", "instrument", "model")])
operators <- operators[order(operators$operator), ]
Write(operators, "operators.csv")

# ---------------------------------------------------------------------------
# Part 2. The automated arm. One rule, all 234 files.
# ---------------------------------------------------------------------------

Say("\nPart 2: gating 234 files")

gate_results <- vector("list", nrow(samples))
started <- Sys.time()
for (index in seq_len(nrow(samples))) {
  row <- samples[index, ]
  gate_results[[index]] <- GateNaiveMemoryFile(
    file.path(kDeposit, row$file_name), panel, row$material
  )
  if (index %% 25 == 0 || index == nrow(samples)) {
    Say("  ", index, " of ", nrow(samples), " files, ",
        round(as.numeric(difftime(Sys.time(), started, units = "mins")), 1),
        " minutes")
  }
}

failed <- vapply(gate_results, function(x) !is.na(x$error), logical(1))
failures <- data.frame(
  file_name = samples$file_name[failed],
  material = samples$material[failed],
  operator = samples$operator[failed],
  reason = vapply(gate_results[failed], function(x) x$error, character(1)),
  stringsAsFactors = FALSE
)
Write(failures, "gating_failures.csv")
Say("  files gated: ", sum(!failed), " of ", nrow(samples),
    ". Files that failed: ", sum(failed))

automated <- do.call(rbind, lapply(gate_results[!failed], function(x) x$counts))
automated <- merge(
  automated,
  samples[, c("file_name", "operator", "centre", "instrument", "model",
              "vial", "round", "key")],
  by = "file_name"
)
Write(automated, "automated_measures.csv")

cuts <- do.call(rbind, lapply(gate_results[!failed], function(x) x$cuts))
Write(cuts, "automated_cuts.csv")

resolution <- do.call(rbind, lapply(gate_results[!failed], function(x) x$channels))
Write(resolution, "channel_resolution.csv")

resolution_summary <- as.data.frame(
  table(marker = resolution$marker, resolved_by = resolution$resolved_by)
)
resolution_summary <- resolution_summary[resolution_summary$Freq > 0, ]
Write(resolution_summary, "channel_resolution_summary.csv")

cut_rules <- as.data.frame(table(marker = cuts$marker, rule = cuts$rule))
cut_rules <- cut_rules[cut_rules$Freq > 0, ]
Write(cut_rules, "cut_rule_summary.csv")

pipeline_states <- data.frame(
  Measure = c(
    "Files whose data was already compensated",
    "Files where the stored matrix was applied",
    "Files that carry no matrix",
    "Files where the time window was applied",
    "Files that needed a wider logicle than the default",
    "Files that needed the fluorochrome fallback on a channel",
    "Channels resolved by the marker name",
    "Channels resolved by the fluorochrome",
    "Cuts fitted by the density minimum",
    "Cuts fitted by the two component mixture"
  ),
  Value = c(
    sum(automated$compensation == "already compensated"),
    sum(automated$compensation == "matrix applied"),
    sum(automated$compensation == "no matrix"),
    sum(automated$time_window_applied),
    sum(automated$logicle_m > 4.5),
    length(unique(resolution$file_name[resolution$resolved_by == "detector"])),
    sum(resolution$resolved_by == "marker"),
    sum(resolution$resolved_by == "detector"),
    sum(cuts$rule == "density"),
    sum(cuts$rule == "mixture")
  ),
  stringsAsFactors = FALSE
)
Write(pipeline_states, "pipeline_states.csv")
print(pipeline_states)

# Two files that share a date, a start time and a GUID are one acquisition
# exported twice. The design calls them separate replicates, so the point has to
# be counted before any spread is reported.
automated$acquisition <- paste(automated$acquisition_date,
                               automated$acquisition_time,
                               automated$acquisition_guid)
repeats <- automated[
  automated$acquisition %in%
    names(which(table(automated$acquisition) > 1)),
  c("file_name", "operator", "material", "vial", "round", "acquisition_date",
    "acquisition_time", "total_events")
]
repeats <- repeats[order(repeats$operator, repeats$material, repeats$vial,
                         repeats$round), ]
Write(repeats, "repeated_acquisitions.csv")
Say("  files that share an acquisition stamp with another file: ", nrow(repeats))

# ---------------------------------------------------------------------------
# Part 2b. The same rule, with one cut per operator instead of one per file.
#
# The first pass fits every cut on the file it gates. That is thirteen automated
# analysts, one per operator, and each of them redraws its gates on every file.
# The reference operator does not work that way. One person sets a boundary and
# holds it across a whole cohort.
#
# This pass takes the median of an operator's own fitted cuts, per marker and per
# material, and gates every file of that operator at that one value. Nothing
# outside the operator's own files is consulted, so the pass adds no information
# from the reference analysis.
# ---------------------------------------------------------------------------

Say("\nPart 2b: gating again with one cut per operator")

cuts_with_design <- merge(
  cuts, samples[, c("file_name", "operator", "material")], by = "file_name"
)
shared_cuts <- aggregate(
  cut ~ operator + material + marker, data = cuts_with_design,
  FUN = stats::median, na.rm = TRUE
)
Write(shared_cuts, "shared_cuts.csv")

SharedCutsFor <- function(operator, material) {
  piece <- shared_cuts[shared_cuts$operator == operator &
                         shared_cuts$material == material, ]
  stats::setNames(piece$cut, piece$marker)
}

shared_results <- vector("list", nrow(samples))
started <- Sys.time()
for (index in seq_len(nrow(samples))) {
  row <- samples[index, ]
  shared_results[[index]] <- GateNaiveMemoryFile(
    file.path(kDeposit, row$file_name), panel, row$material,
    fixed_cuts = SharedCutsFor(row$operator, row$material)
  )
  if (index %% 50 == 0 || index == nrow(samples)) {
    Say("  ", index, " of ", nrow(samples), " files, ",
        round(as.numeric(difftime(Sys.time(), started, units = "mins")), 1),
        " minutes")
  }
}

shared_failed <- vapply(shared_results, function(x) !is.na(x$error), logical(1))
Say("  files gated: ", sum(!shared_failed), " of ", nrow(samples))

shared <- do.call(rbind, lapply(shared_results[!shared_failed],
                                function(x) x$counts))
shared <- merge(
  shared,
  samples[, c("file_name", "operator", "centre", "instrument", "model",
              "vial", "round", "key")],
  by = "file_name"
)
Write(shared, "shared_cut_measures.csv")

# ---------------------------------------------------------------------------
# Part 3. One long table that holds all three arms.
# ---------------------------------------------------------------------------

Say("\nPart 3: the three arms in one table")

published <- results[, c("operator", "centre", "instrument", "material", "vial",
                         "round", "excluded", "parameter", "percent", "analysis",
                         "key")]
published$analysis <- ifelse(published$analysis == "centr", "central", "local")
published$percent <- suppressWarnings(as.numeric(published$percent))

LongForm <- function(frame, label) {
  do.call(rbind, lapply(kParameters, function(parameter) {
    data.frame(
      operator = frame$operator,
      centre = frame$centre,
      instrument = frame$instrument,
      material = frame$material,
      vial = frame$vial,
      round = frame$round,
      excluded = "No",
      parameter = parameter,
      percent = frame[[parameter]],
      analysis = label,
      key = frame$key,
      stringsAsFactors = FALSE
    )
  }))
}

automated_long <- rbind(
  LongForm(automated, "automated"),
  LongForm(shared, "automated shared")
)

# The paper drops six whole blood samples as outliers before every statistic.
# The same six are dropped from all three arms, so the arms stay comparable.
excluded_keys <- unique(published$key[published$excluded == "Yes"])
automated_long$excluded <- ifelse(
  automated_long$key %in% excluded_keys, "Yes", "No"
)

arms <- rbind(published, automated_long)
arms$parameter <- factor(arms$parameter, levels = kParameters)
arms <- arms[!is.na(arms$parameter), ]
Write(arms, "three_arms_long.csv")

kept <- arms[arms$excluded == "No" & is.finite(arms$percent), ]
Say("  rows kept after the paper's exclusions: ", nrow(kept))
Say("  rows per arm: ",
    paste(names(table(kept$analysis)), table(kept$analysis),
          sep = " = ", collapse = ", "))

# The donor mean. Every indicator below starts from the mean of an operator's
# replicates on one donor, which is what the paper does.
DonorMeans <- function(frame) {
  aggregate(
    percent ~ analysis + material + parameter + operator + centre + instrument +
      vial,
    data = frame, FUN = mean, na.rm = TRUE
  )
}

donor_means <- DonorMeans(kept)
Write(donor_means, "donor_means.csv")

# ---------------------------------------------------------------------------
# Part 4. The four indicators the paper defines.
# ---------------------------------------------------------------------------

Say("\nPart 4: the indicators")

combinations <- unique(donor_means[, c("analysis", "material", "parameter")])
combinations <- combinations[
  order(combinations$analysis, combinations$material, combinations$parameter),
]

# The reference value of a parameter and a donor is the mean over operators of
# the central analysis. It is the same reference for every arm, so a bias is
# always measured against the reference operator's reading of the whole cohort.
reference <- aggregate(
  percent ~ material + parameter + vial,
  data = donor_means[donor_means$analysis == "central", ],
  FUN = mean, na.rm = TRUE
)
names(reference)[names(reference) == "percent"] <- "reference"
Write(reference, "reference_values.csv")

# Interoperator CV. For each donor, the spread across operators. The reported
# value is the median of the three donor specific values.
interoperator <- do.call(rbind, lapply(seq_len(nrow(combinations)), function(i) {
  row <- combinations[i, ]
  piece <- donor_means[
    donor_means$analysis == row$analysis &
      donor_means$material == row$material &
      donor_means$parameter == row$parameter,
  ]
  per_vial <- vapply(split(piece$percent, piece$vial), CoefficientOfVariation,
                     numeric(1))
  with_reference <- merge(
    piece, reference[reference$material == row$material, ],
    by = c("material", "parameter", "vial")
  )
  bias <- abs((with_reference$percent - with_reference$reference) /
                with_reference$reference)
  data.frame(
    analysis = row$analysis,
    material = row$material,
    parameter = row$parameter,
    cv = stats::median(per_vial, na.rm = TRUE),
    bias = stats::median(bias[is.finite(bias)], na.rm = TRUE),
    icc = IntraclassCorrelation(piece$percent, piece$operator, piece$vial),
    operators = length(unique(piece$operator)),
    stringsAsFactors = FALSE
  )
}))
Write(interoperator, "interoperator.csv")

interoperator_counts <- do.call(rbind, lapply(
  split(interoperator, list(interoperator$analysis, interoperator$material),
        drop = TRUE),
  function(piece) {
    data.frame(
      analysis = piece$analysis[1],
      material = piece$material[1],
      parameters = nrow(piece),
      cv_acceptable = sum(piece$cv < kCvThreshold, na.rm = TRUE),
      bias_acceptable = sum(piece$bias < kBiasThreshold, na.rm = TRUE),
      icc_acceptable = sum(piece$icc > kIccThreshold, na.rm = TRUE),
      stringsAsFactors = FALSE
    )
  }
))
Write(interoperator_counts, "interoperator_counts.csv")
print(interoperator_counts)

# Intraoperator CV. For one operator and one parameter, the spread across the
# three replicates of a donor. The reported value is the median of the three.
intraoperator <- do.call(rbind, lapply(
  split(kept, list(kept$analysis, kept$material, kept$parameter, kept$operator),
        drop = TRUE),
  function(piece) {
    per_vial <- vapply(split(piece$percent, piece$vial), CoefficientOfVariation,
                       numeric(1))
    data.frame(
      analysis = piece$analysis[1],
      material = piece$material[1],
      parameter = as.character(piece$parameter[1]),
      operator = piece$operator[1],
      cv = stats::median(per_vial, na.rm = TRUE),
      stringsAsFactors = FALSE
    )
  }
))
Write(intraoperator, "intraoperator.csv")

intraoperator_counts <- do.call(rbind, lapply(
  split(intraoperator, list(intraoperator$analysis, intraoperator$material),
        drop = TRUE),
  function(piece) {
    data.frame(
      analysis = piece$analysis[1],
      material = piece$material[1],
      measurements = nrow(piece),
      above_threshold = sum(piece$cv > kCvThreshold, na.rm = TRUE),
      median_cv = stats::median(piece$cv, na.rm = TRUE),
      max_cv = max(piece$cv, na.rm = TRUE),
      stringsAsFactors = FALSE
    )
  }
))
Write(intraoperator_counts, "intraoperator_counts.csv")
print(intraoperator_counts)

# Z-score. The reference distribution of a parameter and a donor is the set of
# centrally analysed operator means. The score of an operator is the median of
# its three donor specific scores.
z_scores <- do.call(rbind, lapply(seq_len(nrow(combinations)), function(i) {
  row <- combinations[i, ]
  piece <- donor_means[
    donor_means$analysis == row$analysis &
      donor_means$material == row$material &
      donor_means$parameter == row$parameter,
  ]
  central <- donor_means[
    donor_means$analysis == "central" &
      donor_means$material == row$material &
      donor_means$parameter == row$parameter,
  ]
  per_vial <- lapply(split(piece, piece$vial), function(vial_piece) {
    pool <- central$percent[central$vial == vial_piece$vial[1]]
    data.frame(
      operator = vial_piece$operator,
      z = ZScore(vial_piece$percent, mean(pool, na.rm = TRUE),
                 stats::sd(pool, na.rm = TRUE)),
      stringsAsFactors = FALSE
    )
  })
  per_vial <- do.call(rbind, per_vial)
  scores <- vapply(split(per_vial$z, per_vial$operator), stats::median,
                   numeric(1), na.rm = TRUE)
  data.frame(
    analysis = row$analysis,
    material = row$material,
    parameter = row$parameter,
    operator = names(scores),
    z = as.numeric(scores),
    stringsAsFactors = FALSE
  )
}))
Write(z_scores, "z_scores.csv")

z_counts <- do.call(rbind, lapply(
  split(z_scores, list(z_scores$analysis, z_scores$material), drop = TRUE),
  function(piece) {
    data.frame(
      analysis = piece$analysis[1],
      material = piece$material[1],
      measurements = nrow(piece),
      unsatisfactory = sum(abs(piece$z) >= kZUnsatisfactory, na.rm = TRUE),
      questionable = sum(abs(piece$z) >= 2 & abs(piece$z) < 3, na.rm = TRUE),
      median_absolute_z = stats::median(abs(piece$z), na.rm = TRUE),
      stringsAsFactors = FALSE
    )
  }
))
Write(z_counts, "z_counts.csv")
print(z_counts)

operator_performance <- do.call(rbind, lapply(
  split(z_scores, list(z_scores$analysis, z_scores$material, z_scores$operator),
        drop = TRUE),
  function(piece) {
    data.frame(
      analysis = piece$analysis[1],
      material = piece$material[1],
      operator = piece$operator[1],
      median_absolute_z = stats::median(abs(piece$z), na.rm = TRUE),
      unsatisfactory = sum(abs(piece$z) >= kZUnsatisfactory, na.rm = TRUE),
      stringsAsFactors = FALSE
    )
  }
))
Write(operator_performance, "operator_performance.csv")

# ---------------------------------------------------------------------------
# Part 5. Agreement after the poor performers are removed.
# ---------------------------------------------------------------------------

Say("\nPart 5: agreement after an exclusion")

poor <- unique(z_scores$operator[
  z_scores$analysis == "central" & z_scores$material == "WB" &
    abs(z_scores$z) > 1.5
])
Say("  operators outside a Z of 1.5 in central whole blood: ",
    paste(sort(poor), collapse = ", "))

icc_after <- do.call(rbind, lapply(kParameters, function(parameter) {
  piece <- donor_means[
    donor_means$analysis == "central" & donor_means$material == "WB" &
      donor_means$parameter == parameter & !(donor_means$operator %in% poor),
  ]
  data.frame(
    parameter = parameter,
    icc_all = interoperator$icc[
      interoperator$analysis == "central" & interoperator$material == "WB" &
        interoperator$parameter == parameter
    ],
    icc_selected = IntraclassCorrelation(piece$percent, piece$operator,
                                         piece$vial),
    stringsAsFactors = FALSE
  )
}))
Write(icc_after, "icc_after_exclusion.csv")
print(icc_after)

# ---------------------------------------------------------------------------
# Part 6. Principal components on the whole blood central analysis.
# ---------------------------------------------------------------------------

Say("\nPart 6: principal components")

PcaVariance <- function(parameters) {
  piece <- donor_means[
    donor_means$analysis == "central" & donor_means$material == "WB" &
      donor_means$parameter %in% parameters,
  ]
  wide <- stats::reshape(
    piece[, c("operator", "vial", "parameter", "percent")],
    idvar = c("operator", "vial"), timevar = "parameter", direction = "wide"
  )
  matrix_values <- as.matrix(wide[, grep("^percent", names(wide))])
  colnames(matrix_values) <- sub("^percent\\.", "", colnames(matrix_values))
  matrix_values <- matrix_values[stats::complete.cases(matrix_values), ,
                                 drop = FALSE]
  fit <- stats::prcomp(matrix_values, center = TRUE, scale. = TRUE)
  variance <- fit$sdev^2 / sum(fit$sdev^2)
  list(fit = fit, first_two = 100 * sum(variance[1:2]), rows = nrow(matrix_values))
}

parent_pca <- PcaVariance(c("CD3", "CD4", "CD8"))
subset_pca <- PcaVariance(setdiff(kParameters, c("CD3", "CD4", "CD8")))

pca_summary <- data.frame(
  group = c("Parent populations", "Naive and memory subsets"),
  variables = c(3, 12),
  observations = c(parent_pca$rows, subset_pca$rows),
  first_two_percent = c(parent_pca$first_two, subset_pca$first_two),
  stringsAsFactors = FALSE
)
Write(pca_summary, "pca_summary.csv")
print(pca_summary)

# ---------------------------------------------------------------------------
# Part 7. The claims.
# ---------------------------------------------------------------------------

Say("\nPart 7: the claims")

Count <- function(frame, analysis, material, column) {
  frame[[column]][frame$analysis == analysis & frame$material == material]
}

Verdict <- function(observed, expected_local, expected_central,
                    observed_local, observed_central) {
  hits <- sum(observed_local == expected_local,
              observed_central == expected_central)
  if (hits == 2) {
    "reproduced"
  } else if (hits == 1) {
    "partly reproduced"
  } else {
    "not reproduced"
  }
}

CountClaim <- function(claim_id, material, column, expected_local,
                       expected_central, frame) {
  observed_local <- Count(frame, "local", material, column)
  observed_central <- Count(frame, "central", material, column)
  data.frame(
    claim_id = claim_id,
    observed = paste0(observed_local, " of 15 local and ",
                      observed_central, " of 15 central"),
    verdict = Verdict(NULL, expected_local, expected_central,
                      observed_local, observed_central),
    stringsAsFactors = FALSE
  )
}

verdicts <- rbind(
  CountClaim(1, "PBMC", "cv_acceptable", 8, 12, interoperator_counts),
  CountClaim(2, "WB", "cv_acceptable", 6, 13, interoperator_counts),
  CountClaim(3, "PBMC", "bias_acceptable", 9, 13, interoperator_counts),
  CountClaim(4, "WB", "bias_acceptable", 7, 14, interoperator_counts),
  CountClaim(5, "PBMC", "icc_acceptable", 5, 7, interoperator_counts),
  CountClaim(6, "WB", "icc_acceptable", 2, 11, interoperator_counts)
)

# Claim 7. The maximum and the median intraoperator CV in frozen PBMC.
pbmc_intra <- intraoperator_counts[intraoperator_counts$material == "PBMC", ]
verdicts <- rbind(verdicts, data.frame(
  claim_id = 7,
  observed = sprintf(
    "maximum %.2f local against %.2f central, median %.2f against %.2f",
    Count(pbmc_intra, "local", "PBMC", "max_cv"),
    Count(pbmc_intra, "central", "PBMC", "max_cv"),
    Count(pbmc_intra, "local", "PBMC", "median_cv"),
    Count(pbmc_intra, "central", "PBMC", "median_cv")
  ),
  verdict = if (Count(pbmc_intra, "local", "PBMC", "max_cv") >
                Count(pbmc_intra, "central", "PBMC", "max_cv")) {
    "reproduced"
  } else {
    "not reproduced"
  },
  stringsAsFactors = FALSE
))

# Claims 8 and 9. How many of the 195 intraoperator CVs sit above 0.20.
for (item in list(list(id = 8, material = "PBMC", local = 52, central = 56),
                  list(id = 9, material = "WB", local = 5, central = 2))) {
  piece <- intraoperator_counts[intraoperator_counts$material == item$material, ]
  observed_local <- Count(piece, "local", item$material, "above_threshold")
  observed_central <- Count(piece, "central", item$material, "above_threshold")
  verdicts <- rbind(verdicts, data.frame(
    claim_id = item$id,
    observed = paste0(observed_local, " of ",
                      Count(piece, "local", item$material, "measurements"),
                      " local and ", observed_central, " of ",
                      Count(piece, "central", item$material, "measurements"),
                      " central"),
    verdict = if ((observed_local > observed_central) ==
                  (item$local > item$central)) {
      "reproduced"
    } else {
      "not reproduced"
    },
    stringsAsFactors = FALSE
  ))
}

# Claims 10 and 11. Unsatisfactory Z-scores, central and local.
central_z <- z_counts[z_counts$analysis == "central", ]
local_z <- z_counts[z_counts$analysis == "local", ]
verdicts <- rbind(verdicts, data.frame(
  claim_id = 10,
  observed = paste0(
    sum(central_z$unsatisfactory), " unsatisfactory in central analysis, ",
    "PBMC ", central_z$unsatisfactory[central_z$material == "PBMC"],
    " and WB ", central_z$unsatisfactory[central_z$material == "WB"]
  ),
  verdict = if (sum(central_z$unsatisfactory) == 0) {
    "reproduced"
  } else {
    "not reproduced"
  },
  stringsAsFactors = FALSE
))
verdicts <- rbind(verdicts, data.frame(
  claim_id = 11,
  observed = paste0(
    "PBMC ", local_z$unsatisfactory[local_z$material == "PBMC"],
    " of ", local_z$measurements[local_z$material == "PBMC"],
    " and WB ", local_z$unsatisfactory[local_z$material == "WB"],
    " of ", local_z$measurements[local_z$material == "WB"]
  ),
  verdict = if (local_z$unsatisfactory[local_z$material == "WB"] >
                local_z$unsatisfactory[local_z$material == "PBMC"]) {
    "reproduced"
  } else {
    "not reproduced"
  },
  stringsAsFactors = FALSE
))

# Claim 12. Where the reference operator ranks.
rop_rank <- do.call(rbind, lapply(
  split(operator_performance[operator_performance$analysis == "local", ],
        operator_performance$material[operator_performance$analysis == "local"]),
  function(piece) {
    piece <- piece[order(piece$median_absolute_z), ]
    data.frame(
      material = piece$material[1],
      rop_rank = which(piece$operator == "M"),
      operators = nrow(piece),
      stringsAsFactors = FALSE
    )
  }
))
Write(rop_rank, "reference_operator_rank.csv")
verdicts <- rbind(verdicts, data.frame(
  claim_id = 12,
  observed = paste0(
    "operator M ranks ",
    paste(rop_rank$rop_rank, "of", rop_rank$operators, "in", rop_rank$material,
          collapse = ", ")
  ),
  verdict = if (all(rop_rank$rop_rank <= 5)) "reproduced" else "not reproduced",
  stringsAsFactors = FALSE
))

# Claim 13. CD8 naive is the least reliable parameter in frozen PBMC.
pbmc_central <- interoperator[
  interoperator$analysis == "central" & interoperator$material == "PBMC",
]
cd8n <- pbmc_central[pbmc_central$parameter == "CD8_N", ]
verdicts <- rbind(verdicts, data.frame(
  claim_id = 13,
  observed = sprintf("CD8_N cv %.2f bias %.2f icc %.2f",
                     cd8n$cv, cd8n$bias, cd8n$icc),
  verdict = if (nrow(cd8n) == 1 && cd8n$cv > kCvThreshold &&
                cd8n$bias > kBiasThreshold && cd8n$icc < kIccThreshold) {
    "reproduced"
  } else {
    "not reproduced"
  },
  stringsAsFactors = FALSE
))

# Claim 14. CD4 terminally differentiated carries the worst CV and bias.
worst <- do.call(rbind, lapply(
  split(interoperator, list(interoperator$analysis, interoperator$material),
        drop = TRUE),
  function(piece) {
    data.frame(
      analysis = piece$analysis[1],
      material = piece$material[1],
      worst_cv = piece$parameter[which.max(piece$cv)],
      worst_bias = piece$parameter[which.max(piece$bias)],
      stringsAsFactors = FALSE
    )
  }
))
Write(worst, "worst_parameters.csv")
verdicts <- rbind(verdicts, data.frame(
  claim_id = 14,
  observed = paste0(
    "the worst CV is ",
    paste(worst$analysis, worst$material, worst$worst_cv, collapse = ", ")
  ),
  verdict = if (sum(as.character(worst$worst_cv) == "CD4_TD") >= 2) {
    "reproduced"
  } else {
    "partly reproduced"
  },
  stringsAsFactors = FALSE
))

# Claim 15. Precision follows population size.
size_check <- do.call(rbind, lapply(
  split(interoperator, list(interoperator$analysis, interoperator$material),
        drop = TRUE),
  function(piece) {
    parents <- piece$cv[piece$parameter %in% c("CD3", "CD4", "CD8")]
    subsets <- piece$cv[!piece$parameter %in% c("CD3", "CD4", "CD8")]
    data.frame(
      analysis = piece$analysis[1],
      material = piece$material[1],
      parent_median_cv = stats::median(parents, na.rm = TRUE),
      subset_median_cv = stats::median(subsets, na.rm = TRUE),
      stringsAsFactors = FALSE
    )
  }
))
Write(size_check, "cv_by_population_size.csv")
verdicts <- rbind(verdicts, data.frame(
  claim_id = 15,
  observed = sprintf(
    "the parent median CV is %.3f against %.3f for the subsets, pooled",
    stats::median(size_check$parent_median_cv),
    stats::median(size_check$subset_median_cv)
  ),
  verdict = if (all(size_check$parent_median_cv < size_check$subset_median_cv)) {
    "reproduced"
  } else {
    "partly reproduced"
  },
  stringsAsFactors = FALSE
))

# Claim 16. Agreement after the poor performers are removed.
still_low <- icc_after$parameter[
  is.finite(icc_after$icc_selected) & icc_after$icc_selected < kIccThreshold
]
verdicts <- rbind(verdicts, data.frame(
  claim_id = 16,
  observed = paste0(length(still_low), " of 15 still below 0.75: ",
                    paste(still_low, collapse = ", ")),
  verdict = if (length(still_low) <=
                sum(icc_after$icc_all < kIccThreshold, na.rm = TRUE)) {
    "reproduced"
  } else {
    "not reproduced"
  },
  stringsAsFactors = FALSE
))

# Claims 17 and 18. The variance in the first two principal components.
for (item in list(list(id = 17, group = "Parent populations", expected = 99.4),
                  list(id = 18, group = "Naive and memory subsets",
                       expected = 63.9))) {
  observed <- pca_summary$first_two_percent[pca_summary$group == item$group]
  verdicts <- rbind(verdicts, data.frame(
    claim_id = item$id,
    observed = sprintf("%.1f percent", observed),
    verdict = if (abs(observed - item$expected) <= 10) {
      "reproduced"
    } else {
      "not reproduced"
    },
    stringsAsFactors = FALSE
  ))
}

# Claim 19. Dispersion is wider in locally analysed whole blood, and
# centralisation closes more of it there.
spread <- do.call(rbind, lapply(
  split(donor_means, list(donor_means$analysis, donor_means$material),
        drop = TRUE),
  function(piece) {
    per_parameter <- vapply(
      split(piece, list(piece$parameter, piece$vial), drop = TRUE),
      function(cell) CoefficientOfVariation(cell$percent),
      numeric(1)
    )
    data.frame(
      analysis = piece$analysis[1],
      material = piece$material[1],
      median_cv = stats::median(per_parameter, na.rm = TRUE),
      stringsAsFactors = FALSE
    )
  }
))
Write(spread, "dispersion_by_arm.csv")
pbmc_drop <- spread$median_cv[spread$analysis == "local" &
                                spread$material == "PBMC"] -
  spread$median_cv[spread$analysis == "central" & spread$material == "PBMC"]
wb_drop <- spread$median_cv[spread$analysis == "local" &
                              spread$material == "WB"] -
  spread$median_cv[spread$analysis == "central" & spread$material == "WB"]
verdicts <- rbind(verdicts, data.frame(
  claim_id = 19,
  observed = sprintf(
    "centralisation lowers the median CV by %.3f in whole blood and by %.3f in PBMC",
    wb_drop, pbmc_drop
  ),
  verdict = if (wb_drop > pbmc_drop) "reproduced" else "not reproduced",
  stringsAsFactors = FALSE
))

# Claim 20. The frozen PBMC preparation carries an unexpected phenotype.
phenotype <- do.call(rbind, lapply(
  split(donor_means[donor_means$analysis == "central", ],
        donor_means$material[donor_means$analysis == "central"]),
  function(piece) {
    Median <- function(parameter) {
      stats::median(piece$percent[piece$parameter == parameter], na.rm = TRUE)
    }
    data.frame(
      material = piece$material[1],
      cd4 = Median("CD4"),
      cd8 = Median("CD8"),
      cd4_over_cd8 = Median("CD4") / Median("CD8"),
      cd8_td = Median("CD8_TD"),
      stringsAsFactors = FALSE
    )
  }
))
Write(phenotype, "phenotype_by_material.csv")
verdicts <- rbind(verdicts, data.frame(
  claim_id = 20,
  observed = sprintf(
    "CD4 over CD8 is %.2f in PBMC against %.2f in whole blood, and CD8_TD is %.1f against %.1f percent",
    phenotype$cd4_over_cd8[phenotype$material == "PBMC"],
    phenotype$cd4_over_cd8[phenotype$material == "WB"],
    phenotype$cd8_td[phenotype$material == "PBMC"],
    phenotype$cd8_td[phenotype$material == "WB"]
  ),
  verdict = if (phenotype$cd4_over_cd8[phenotype$material == "PBMC"] <
                phenotype$cd4_over_cd8[phenotype$material == "WB"] &&
                phenotype$cd8_td[phenotype$material == "PBMC"] >
                phenotype$cd8_td[phenotype$material == "WB"]) {
    "reproduced"
  } else {
    "not reproduced"
  },
  stringsAsFactors = FALSE
))

claim_verdicts <- merge(claims, verdicts, by = "claim_id", all.x = TRUE)
claim_verdicts <- claim_verdicts[order(claim_verdicts$claim_id), ]
Write(claim_verdicts, "claim_verdicts.csv")
print(claim_verdicts[, c("claim_id", "short_name", "expected", "observed",
                         "verdict")])
Say("\n  verdicts: ",
    paste(names(table(claim_verdicts$verdict)), table(claim_verdicts$verdict),
          sep = " = ", collapse = ", "))

# ---------------------------------------------------------------------------
# Part 8. Figures.
# ---------------------------------------------------------------------------

Say("\nPart 8: figures")

kArmColours <- c(local = "#B2182B", central = "#2166AC",
                 automated = "#1B7837", `automated shared` = "#762A83")
kArmOrder <- c("local", "central", "automated", "automated shared")

Save <- function(plot, name, width = 9, height = 6) {
  ggplot2::ggsave(file.path(kOutputDir, name), plot, width = width,
                  height = height, dpi = 150)
}

interoperator$analysis <- factor(interoperator$analysis, levels = kArmOrder)
interoperator$parameter <- factor(interoperator$parameter, levels = kParameters)

Save(
  ggplot(interoperator, aes(parameter, cv, fill = analysis)) +
    geom_col(position = position_dodge(width = 0.8), width = 0.75) +
    geom_hline(yintercept = kCvThreshold, linetype = "dashed") +
    facet_wrap(~material, ncol = 1) +
    scale_fill_manual(values = kArmColours) +
    labs(x = NULL, y = "Interoperator coefficient of variation", fill = NULL) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)),
  "interoperator_cv.png", height = 7
)

Save(
  ggplot(interoperator, aes(parameter, icc, fill = analysis)) +
    geom_col(position = position_dodge(width = 0.8), width = 0.75) +
    geom_hline(yintercept = kIccThreshold, linetype = "dashed") +
    facet_wrap(~material, ncol = 1) +
    scale_fill_manual(values = kArmColours) +
    labs(x = NULL, y = "Intraclass correlation", fill = NULL) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)),
  "interoperator_icc.png", height = 7
)

Save(
  ggplot(interoperator, aes(parameter, bias, fill = analysis)) +
    geom_col(position = position_dodge(width = 0.8), width = 0.75) +
    geom_hline(yintercept = kBiasThreshold, linetype = "dashed") +
    facet_wrap(~material, ncol = 1) +
    scale_fill_manual(values = kArmColours) +
    labs(x = NULL, y = "Bias against the reference value", fill = NULL) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)),
  "interoperator_bias.png", height = 7
)

z_scores$analysis <- factor(z_scores$analysis, levels = kArmOrder)
z_scores$parameter <- factor(z_scores$parameter, levels = kParameters)
Save(
  ggplot(z_scores, aes(parameter, operator, fill = z)) +
    geom_tile() +
    facet_grid(material ~ analysis) +
    scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B",
                         midpoint = 0, limits = c(-4, 4), oob = scales::squish) +
    labs(x = NULL, y = "Operator", fill = "Z") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)),
  "z_score_heatmap.png", width = 11, height = 7
)

plot_frame <- donor_means
plot_frame$analysis <- factor(plot_frame$analysis, levels = kArmOrder)
plot_frame$parameter <- factor(plot_frame$parameter, levels = kParameters)
Save(
  ggplot(plot_frame, aes(parameter, percent, fill = analysis)) +
    geom_boxplot(outlier.size = 0.6, position = position_dodge(width = 0.8)) +
    facet_wrap(~material, ncol = 1, scales = "free_y") +
    scale_fill_manual(values = kArmColours) +
    labs(x = NULL, y = "Percent of parent", fill = NULL) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)),
  "distribution_by_arm.png", height = 8
)

intraoperator$analysis <- factor(intraoperator$analysis, levels = kArmOrder)
Save(
  ggplot(intraoperator, aes(operator, cv, fill = analysis)) +
    geom_boxplot(outlier.size = 0.5) +
    geom_hline(yintercept = kCvThreshold, linetype = "dashed") +
    facet_wrap(~material, ncol = 1, scales = "free_y") +
    scale_fill_manual(values = kArmColours) +
    labs(x = "Operator", y = "Intraoperator coefficient of variation",
         fill = NULL) +
    theme_bw(),
  "intraoperator_cv.png", height = 7
)

Say("\nDone. Every table is in ", kOutputDir)
