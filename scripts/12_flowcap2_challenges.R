# The FlowCAP II sample classification challenges.
#
# FR-FCM-ZZZV is Challenge 3, which asks a pipeline to tell a Gag stimulated
# sample from an Env stimulated one. FR-FCM-ZZZU is Challenge 1, which asks it to
# tell an HIV exposed uninfected infant from an unexposed one.
#
# Both are scored the way the benchmark scored them. A rule is fitted on the
# training half and it is applied, unchanged, to the half that was held back.
#
# Run it in the container:
#   podman run --rm -v "$PWD:/work:z" -w /work everything-flow-cytometry:latest \
#     Rscript scripts/12_flowcap2_challenges.R

suppressPackageStartupMessages({
  library(flowCore)
  library(ggplot2)
  library(withr)
})

source(file.path("R", "figures.R"))
source(file.path("R", "harmonisation.R"))
source(file.path("R", "spectral.R"))
source(file.path("R", "naive_memory.R"))
source(file.path("R", "flowcap.R"))

kHvtnDir <- file.path("data", "datasets", "flowrepository",
                      "FlowRepository_FR-FCM-ZZZV_files")
kHeuDir <- file.path("data", "datasets", "flowrepository", "FR-FCM-ZZZU")
kGatingDir <- "gating"
kOutputDir <- file.path("output", "flowcap2")

kCytokinesT <- c("IL2", "IL4", "IFNg", "TNFa")
kCytokinesInnate <- c("IFNa", "IL6", "IL12", "TNFa")

dir.create(kOutputDir, recursive = TRUE, showWarnings = FALSE)
Write <- function(x, name) {
  utils::write.csv(x, file.path(kOutputDir, name), row.names = FALSE)
  invisible(x)
}
Say <- function(...) cat(..., "\n", sep = "")

claims <- utils::read.csv(
  file.path(kGatingDir, "flowcap2_paper_claims.csv"), stringsAsFactors = FALSE
)

# A gating result carries a column only when the gate it belongs to resolved, so
# the rows are padded to a common set before they are stacked.
BindRows <- function(frames) {
  frames <- Filter(Negate(is.null), frames)
  if (length(frames) == 0) {
    return(NULL)
  }
  columns <- unique(unlist(lapply(frames, names)))
  do.call(rbind, lapply(frames, function(frame) {
    for (name in setdiff(columns, names(frame))) {
      frame[[name]] <- NA
    }
    frame[, columns, drop = FALSE]
  }))
}

# ---------------------------------------------------------------------------
# Part 1. Challenge 3, the HVTN deposit.
# ---------------------------------------------------------------------------

Say("Part 1: Challenge 3, the HVTN deposit")

hvtn_meta <- utils::read.csv(
  file.path(kHvtnDir, "attachments", "Challenge3Metadata.csv"),
  stringsAsFactors = FALSE, check.names = FALSE
)
names(hvtn_meta) <- c("treatment", "subject", "file_name", "split")
Write(hvtn_meta, "hvtn_metadata.csv")

hvtn_design <- as.data.frame(table(treatment = hvtn_meta$treatment,
                                   split = hvtn_meta$split))
Write(hvtn_design, "hvtn_design.csv")
print(hvtn_design)

present <- file.exists(file.path(kHvtnDir, hvtn_meta$file_name))
Say("  files named in the metadata that are present: ", sum(present), " of ",
    nrow(hvtn_meta))
Say("  subjects: ", length(unique(hvtn_meta$subject)),
    ", training subjects: ",
    length(unique(hvtn_meta$subject[hvtn_meta$split == "training"])),
    ", testing subjects: ",
    length(unique(hvtn_meta$subject[hvtn_meta$split == "testing"])))

Say("\n  gating ", nrow(hvtn_meta), " files, subject by subject")

subjects <- unique(hvtn_meta$subject)
started <- Sys.time()
hvtn_rows <- list()
hvtn_counts <- list()
hvtn_failures <- list()

for (index in seq_along(subjects)) {
  subject <- subjects[index]
  piece <- hvtn_meta[hvtn_meta$subject == subject, ]
  gated <- lapply(piece$file_name, function(name) {
    GateHvtnSubsets(file.path(kHvtnDir, name), cytokines = kCytokinesT)
  })
  names(gated) <- piece$file_name

  for (position in seq_along(gated)) {
    if (is.na(gated[[position]]$error)) {
      hvtn_counts[[length(hvtn_counts) + 1]] <- cbind(
        subject = subject, treatment = piece$treatment[position],
        split = piece$split[position], gated[[position]]$counts
      )
    } else {
      hvtn_failures[[length(hvtn_failures) + 1]] <- data.frame(
        file_name = piece$file_name[position], subject = subject,
        treatment = piece$treatment[position],
        reason = gated[[position]]$error, stringsAsFactors = FALSE
      )
    }
  }

  controls <- which(piece$treatment == "negctrl" &
                      vapply(gated, function(x) is.na(x$error), logical(1)))
  if (length(controls) < 2) {
    next
  }

  # One control sets the threshold and the other estimates the background at that
  # same threshold, so a stimulated frequency is reported above a background that
  # was measured rather than assumed.
  reference <- gated[[controls[1]]]$cytokines
  background_source <- gated[[controls[2]]]$cytokines

  for (position in seq_along(gated)) {
    treatment <- piece$treatment[position]
    if (treatment == "negctrl" || !is.na(gated[[position]]$error)) {
      next
    }
    row <- data.frame(
      subject = subject, treatment = treatment,
      split = piece$split[position], file_name = piece$file_name[position],
      stringsAsFactors = FALSE
    )
    for (subset in c("CD4", "CD8")) {
      background <- ControlGatedFrequencies(reference[[subset]],
                                            background_source[[subset]])
      stimulated <- ControlGatedFrequencies(
        reference[[subset]], gated[[position]]$cytokines[[subset]]
      )
      for (cytokine in kCytokinesT) {
        row[[paste0(subset, "_", cytokine)]] <-
          stimulated[[cytokine]] - background[[cytokine]]
      }
    }
    hvtn_rows[[length(hvtn_rows) + 1]] <- row
  }

  if (index %% 8 == 0 || index == length(subjects)) {
    Say("    ", index, " of ", length(subjects), " subjects, ",
        round(as.numeric(difftime(Sys.time(), started, units = "mins")), 1),
        " minutes")
  }
}

hvtn_features <- do.call(rbind, hvtn_rows)
Write(hvtn_features, "hvtn_features.csv")
Write(BindRows(hvtn_counts), "hvtn_counts.csv")
hvtn_failure_table <- if (length(hvtn_failures) > 0) {
  do.call(rbind, hvtn_failures)
} else {
  data.frame(file_name = character(), subject = character(),
             treatment = character(), reason = character(),
             stringsAsFactors = FALSE)
}
Write(hvtn_failure_table, "hvtn_failures.csv")
Say("  files that failed to gate: ", nrow(hvtn_failure_table))
Say("  stimulated samples with features: ", nrow(hvtn_features))

# ---------------------------------------------------------------------------
# Part 2. The population the paper names.
# ---------------------------------------------------------------------------

Say("\nPart 2: Env against Gag")

antigen <- hvtn_features[hvtn_features$treatment %in%
                           c("ENV-1-PTEG", "GAG-1-PTEG"), ]
antigen$is_env <- antigen$treatment == "ENV-1-PTEG"

feature_names <- grep("^(CD4|CD8)_", names(antigen), value = TRUE)
paired <- do.call(rbind, lapply(feature_names, function(name) {
  env <- antigen[antigen$is_env, c("subject", name)]
  gag <- antigen[!antigen$is_env, c("subject", name)]
  both <- merge(env, gag, by = "subject", suffixes = c("_env", "_gag"))
  usable <- is.finite(both[[paste0(name, "_env")]]) &
    is.finite(both[[paste0(name, "_gag")]])
  both <- both[usable, ]
  if (nrow(both) < 5) {
    return(NULL)
  }
  test <- suppressWarnings(stats::wilcox.test(
    both[[paste0(name, "_env")]], both[[paste0(name, "_gag")]], paired = TRUE
  ))
  data.frame(
    feature = name,
    subjects = nrow(both),
    env_median = stats::median(both[[paste0(name, "_env")]]),
    gag_median = stats::median(both[[paste0(name, "_gag")]]),
    subjects_env_higher = sum(both[[paste0(name, "_env")]] >
                                both[[paste0(name, "_gag")]]),
    p_value = test$p.value,
    stringsAsFactors = FALSE
  )
}))
paired <- paired[order(paired$p_value), ]
Write(paired, "env_against_gag.csv")
print(paired)

# ---------------------------------------------------------------------------
# Part 3. The classification, scored the way the benchmark scored it.
# ---------------------------------------------------------------------------

Say("\nPart 3: the Challenge 3 classification")

hvtn_matrix <- antigen[, feature_names]
hvtn_labels <- antigen$is_env
hvtn_training <- antigen$split == "training"
usable <- stats::complete.cases(hvtn_matrix)
Say("  samples with a complete feature vector: ", sum(usable), " of ",
    nrow(antigen))

hvtn_result <- SelectAndScore(hvtn_matrix[usable, ], hvtn_labels[usable],
                              hvtn_training[usable])
Write(hvtn_result$ranking, "hvtn_feature_ranking.csv")
Say("  selected feature: ", hvtn_result$selected,
    " (", hvtn_result$direction, " means Env)")
print(rbind(
  cbind(half = "training", hvtn_result$training_score),
  cbind(half = "testing", hvtn_result$testing_score)
))
Write(rbind(
  cbind(challenge = "HVTN", half = "training", hvtn_result$training_score),
  cbind(challenge = "HVTN", half = "testing", hvtn_result$testing_score)
), "hvtn_scores.csv")

# ---------------------------------------------------------------------------
# Part 4. Challenge 1, the HEUvsUE deposit.
# ---------------------------------------------------------------------------

Say("\nPart 4: Challenge 1, the HEUvsUE deposit")

heu_meta <- utils::read.csv(
  file.path(kHeuDir, "attachments", "HEUvsUE.csv"),
  stringsAsFactors = FALSE, check.names = FALSE
)
names(heu_meta) <- c("file_name", "condition", "dose", "individual")
Write(heu_meta, "heu_metadata.csv")

heu_design <- as.data.frame(table(condition = heu_meta$condition,
                                  dose = heu_meta$dose))
Write(heu_design, "heu_design.csv")
Say("  individuals: ", length(unique(heu_meta$individual)),
    ", conditions: ", paste(sort(unique(heu_meta$dose)), collapse = ", "))
Say("  files: ", nrow(heu_meta), ", HEU ", sum(heu_meta$condition == "HEU"),
    ", UE ", sum(heu_meta$condition == "UE"))

heu_present <- file.exists(file.path(kHeuDir, heu_meta$file_name))
Say("  files named in the metadata that are present: ", sum(heu_present),
    " of ", nrow(heu_meta))

Say("\n  gating ", nrow(heu_meta), " files, individual by individual")

individuals <- sort(unique(heu_meta$individual))
started <- Sys.time()
heu_rows <- list()
heu_counts <- list()
heu_failures <- list()

for (index in seq_along(individuals)) {
  person <- individuals[index]
  piece <- heu_meta[heu_meta$individual == person, ]
  gated <- lapply(piece$file_name, function(name) {
    GateHeuSubsets(file.path(kHeuDir, name), cytokines = kCytokinesInnate)
  })
  names(gated) <- piece$file_name

  for (position in seq_along(gated)) {
    if (is.na(gated[[position]]$error)) {
      heu_counts[[length(heu_counts) + 1]] <- cbind(
        individual = person, condition = piece$condition[position],
        dose = piece$dose[position], gated[[position]]$counts
      )
    } else {
      heu_failures[[length(heu_failures) + 1]] <- data.frame(
        file_name = piece$file_name[position], individual = person,
        dose = piece$dose[position], reason = gated[[position]]$error,
        stringsAsFactors = FALSE
      )
    }
  }

  control <- which(piece$dose == "unstim" &
                     vapply(gated, function(x) is.na(x$error), logical(1)))
  if (length(control) == 0) {
    next
  }
  reference <- gated[[control[1]]]$cytokines

  for (position in seq_along(gated)) {
    if (piece$dose[position] == "unstim" || !is.na(gated[[position]]$error)) {
      next
    }
    row <- data.frame(
      individual = person, condition = piece$condition[position],
      dose = piece$dose[position], file_name = piece$file_name[position],
      stringsAsFactors = FALSE
    )
    for (subset in names(reference)) {
      values <- ControlGatedFrequencies(
        reference[[subset]], gated[[position]]$cytokines[[subset]]
      )
      for (cytokine in kCytokinesInnate) {
        row[[paste0(subset, "_", cytokine)]] <- if (cytokine %in% names(values)) {
          values[[cytokine]]
        } else {
          NA_real_
        }
      }
    }
    heu_rows[[length(heu_rows) + 1]] <- row
  }

  if (index %% 8 == 0 || index == length(individuals)) {
    Say("    ", index, " of ", length(individuals), " individuals, ",
        round(as.numeric(difftime(Sys.time(), started, units = "mins")), 1),
        " minutes")
  }
}

heu_long <- BindRows(heu_rows)
Write(heu_long, "heu_features_long.csv")
Write(BindRows(heu_counts), "heu_counts.csv")
heu_failure_table <- if (length(heu_failures) > 0) {
  do.call(rbind, heu_failures)
} else {
  data.frame(file_name = character(), individual = character(),
             dose = character(), reason = character(), stringsAsFactors = FALSE)
}
Write(heu_failure_table, "heu_failures.csv")
Say("  files that failed to gate: ", nrow(heu_failure_table))

# One feature vector per individual, which is every cytokine of every population
# under every stimulation.
measure_names <- setdiff(names(heu_long),
                         c("individual", "condition", "dose", "file_name"))
heu_wide <- NULL
for (dose in sort(unique(heu_long$dose))) {
  piece <- heu_long[heu_long$dose == dose,
                    c("individual", "condition", measure_names)]
  names(piece)[-(1:2)] <- paste0(dose, "_", measure_names)
  heu_wide <- if (is.null(heu_wide)) {
    piece
  } else {
    merge(heu_wide, piece[, -2], by = "individual", all = TRUE)
  }
}
Write(heu_wide, "heu_features.csv")
Say("  individuals with a feature vector: ", nrow(heu_wide),
    ", features per individual: ", ncol(heu_wide) - 2)

# The deposit carries no training and testing split, so one is made here. It
# splits on the individual, it is deterministic, and it is stated in the report.
heu_wide <- heu_wide[order(heu_wide$individual), ]
heu_wide$split <- ifelse(seq_len(nrow(heu_wide)) %% 2 == 1, "training",
                         "testing")
heu_features_only <- heu_wide[, setdiff(names(heu_wide),
                                        c("individual", "condition", "split"))]
keep_feature <- vapply(heu_features_only, function(column) {
  sum(is.finite(column)) >= 0.8 * nrow(heu_wide)
}, logical(1))
heu_features_only <- heu_features_only[, keep_feature, drop = FALSE]
usable_heu <- stats::complete.cases(heu_features_only)
Say("  features kept after the completeness filter: ", ncol(heu_features_only))
Say("  individuals with a complete feature vector: ", sum(usable_heu))

heu_result <- SelectAndScore(
  heu_features_only[usable_heu, ],
  heu_wide$condition[usable_heu] == "HEU",
  heu_wide$split[usable_heu] == "training"
)
Write(heu_result$ranking, "heu_feature_ranking.csv")
Say("  selected feature: ", heu_result$selected,
    " (", heu_result$direction, " means HIV exposed)")
print(rbind(
  cbind(half = "training", heu_result$training_score),
  cbind(half = "testing", heu_result$testing_score)
))
Write(rbind(
  cbind(challenge = "HEUvsUE", half = "training", heu_result$training_score),
  cbind(challenge = "HEUvsUE", half = "testing", heu_result$testing_score)
), "heu_scores.csv")

scores <- rbind(
  cbind(challenge = "HVTN", half = "training", hvtn_result$training_score),
  cbind(challenge = "HVTN", half = "testing", hvtn_result$testing_score),
  cbind(challenge = "HEUvsUE", half = "training", heu_result$training_score),
  cbind(challenge = "HEUvsUE", half = "testing", heu_result$testing_score)
)
Write(scores, "scores.csv")
print(scores)

# ---------------------------------------------------------------------------
# Part 5. The claims.
# ---------------------------------------------------------------------------

Say("\nPart 5: the claims")

Verdict <- function(id, observed, verdict) {
  data.frame(claim_id = id, observed = observed, verdict = verdict,
             stringsAsFactors = FALSE)
}

hvtn_test <- hvtn_result$testing_score
heu_test <- heu_result$testing_score
heu_train <- heu_result$training_score

verdicts <- rbind(
  Verdict(1, sprintf(
    "one feature reaches an F-measure of %.2f and an accuracy of %.2f on the held out half",
    hvtn_test$f_measure, hvtn_test$accuracy
  ), if (hvtn_test$f_measure >= 0.8) {
    "reproduced"
  } else if (hvtn_test$f_measure >= 0.65) {
    "partly reproduced"
  } else {
    "not reproduced"
  }),
  Verdict(2, sprintf(
    "the same rule reaches an F-measure of %.2f and an accuracy of %.2f on the held out half",
    heu_test$f_measure, heu_test$accuracy
  ), if (heu_test$accuracy <= 0.7) "reproduced" else "not reproduced"),
  Verdict(6, sprintf(
    "training accuracy %.2f against testing accuracy %.2f",
    heu_train$accuracy, heu_test$accuracy
  ), if (heu_train$accuracy > heu_test$accuracy) {
    "reproduced"
  } else {
    "not reproduced"
  })
)

cd4_il2 <- paired[paired$feature == "CD4_IL2", ]
verdicts <- rbind(verdicts, Verdict(
  3, sprintf(
    "CD4 IL2 is %.3f in Env against %.3f in Gag, higher in %d of %d subjects, p %.4f",
    cd4_il2$env_median, cd4_il2$gag_median, cd4_il2$subjects_env_higher,
    cd4_il2$subjects, cd4_il2$p_value
  ),
  if (cd4_il2$env_median > cd4_il2$gag_median && cd4_il2$p_value < 0.05) {
    "reproduced"
  } else if (cd4_il2$env_median > cd4_il2$gag_median) {
    "partly reproduced"
  } else {
    "not reproduced"
  }
))

verdicts <- rbind(verdicts, Verdict(
  4, sprintf("the selected feature is %s", hvtn_result$selected),
  if (hvtn_result$selected == "CD4_IL2") {
    "reproduced"
  } else if (grepl("IL2", hvtn_result$selected)) {
    "partly reproduced"
  } else {
    "not reproduced"
  }
))

split_counts <- table(hvtn_meta$split)
verdicts <- rbind(verdicts, Verdict(
  5, sprintf("the HVTN deposit is %d training and %d testing files",
             split_counts[["training"]], split_counts[["testing"]]),
  if (split_counts[["training"]] == split_counts[["testing"]]) {
    "reproduced"
  } else {
    "not reproduced"
  }
))

claim_verdicts <- merge(claims, verdicts, by = "claim_id", all.x = TRUE)
claim_verdicts <- claim_verdicts[order(claim_verdicts$claim_id), ]
Write(claim_verdicts, "claim_verdicts.csv")
print(claim_verdicts[, c("claim_id", "short_name", "expected", "observed",
                         "verdict")])
Say("\n  verdicts: ",
    paste(names(table(claim_verdicts$verdict)),
          table(claim_verdicts$verdict), sep = " = ", collapse = ", "))

# ---------------------------------------------------------------------------
# Part 6. Figures.
# ---------------------------------------------------------------------------

Say("\nPart 6: figures")

Save <- function(plot, name, width = 9, height = 6) {
  SaveFigure(plot, file.path(kOutputDir, name), width = width, height = height)
}

long_antigen <- do.call(rbind, lapply(feature_names, function(name) {
  data.frame(feature = name, value = antigen[[name]],
             treatment = ifelse(antigen$is_env, "Env", "Gag"),
             stringsAsFactors = FALSE)
}))
Save(
  ggplot(long_antigen[is.finite(long_antigen$value), ],
         aes(treatment, value, fill = treatment)) +
    geom_boxplot(outlier.size = 0.6) +
    facet_wrap(~feature, scales = "free_y") +
    scale_fill_manual(values = c(Env = "#B2182B", Gag = "#2166AC")) +
    labs(x = NULL, y = "Percent above the control threshold, background removed",
         fill = NULL) +
    ThemePublication(),
  "env_against_gag.svg", width = 10, height = 7
)

selected <- hvtn_result$selected
Save(
  ggplot(antigen[is.finite(antigen[[selected]]), ],
         aes(split, .data[[selected]], colour = ifelse(is_env, "Env", "Gag"))) +
    geom_hline(yintercept = hvtn_result$threshold, linetype = "dashed") +
    geom_jitter(width = 0.15, height = 0, size = 2.5) +
    scale_colour_manual(values = c(Env = "#B2182B", Gag = "#2166AC")) +
    labs(x = NULL, y = paste(selected, "percent"), colour = NULL,
         title = "The selected feature and the threshold fitted on the training half") +
    ThemePublication(),
  "hvtn_selected_feature.svg", width = 8, height = 5
)

score_long <- do.call(rbind, lapply(c("accuracy", "f_measure"),
                                    function(measure) {
  data.frame(challenge = scores$challenge, half = scores$half,
             measure = measure, value = scores[[measure]],
             stringsAsFactors = FALSE)
}))
Save(
  ggplot(score_long, aes(half, value, fill = challenge)) +
    geom_col(position = position_dodge(width = 0.8), width = 0.7) +
    geom_hline(yintercept = 0.5, linetype = "dashed") +
    facet_wrap(~measure) +
    scale_fill_manual(values = c(HVTN = "#1B7837", HEUvsUE = "#762A83")) +
    ylim(0, 1) +
    labs(x = NULL, y = NULL, fill = NULL) +
    ThemePublication(),
  "scores.svg", width = 8, height = 5
)

Say("\nDone. Every table is in ", kOutputDir)
