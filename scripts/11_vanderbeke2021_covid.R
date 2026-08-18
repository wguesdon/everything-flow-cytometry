# COVID-19 intracellular cytokine staining, Vanderbeke 2021.
#
# FR-FCM-Z2KP holds 49 files and the FlowJo workspace the authors used. Two arms
# measure the same files.
#
#   manual     the deposited workspace, applied by CytoML
#   automated  one rule, fitted on each file
#
# This is the first deposit in the repository whose workspace gates every sample,
# so the manual arm covers the whole cohort rather than one donor.
#
# Run it in the container:
#   podman run --rm -v "$PWD:/work:z" -w /work everything-flow-cytometry:latest \
#     Rscript scripts/11_vanderbeke2021_covid.R

suppressPackageStartupMessages({
  library(flowCore)
  library(flowWorkspace)
  library(CytoML)
  library(ggplot2)
  library(withr)
})

source(file.path("R", "figures.R"))
source(file.path("R", "harmonisation.R"))
source(file.path("R", "spectral.R"))
source(file.path("R", "naive_memory.R"))
source(file.path("R", "covid_ics.R"))

kDeposit <- file.path("data", "datasets", "flowrepository", "FR-FCM-Z2KP")
kWorkspace <- file.path(kDeposit, "attachments",
                        "01-May-2020_Human_COVID_analysis_template.wsp")
kGatingDir <- "gating"
kOutputDir <- file.path("output", "vanderbeke2021")

kGroupOrder <- c("healthy", "ward", "intensive care")

dir.create(kOutputDir, recursive = TRUE, showWarnings = FALSE)
Write <- function(x, name) {
  utils::write.csv(x, file.path(kOutputDir, name), row.names = FALSE)
  invisible(x)
}
Say <- function(...) cat(..., "\n", sep = "")

# ---------------------------------------------------------------------------
# Part 1. The deposit.
# ---------------------------------------------------------------------------

Say("Part 1: the deposit")

claims <- utils::read.csv(
  file.path(kGatingDir, "vanderbeke2021_paper_claims.csv"),
  stringsAsFactors = FALSE
)

files <- ParseCovidFileNames(list.files(kDeposit, pattern = "\\.fcs$"))
files$group <- factor(files$group, levels = kGroupOrder)
Write(files, "file_index.csv")

design <- as.data.frame(table(group = files$group))
Write(design, "design.csv")
print(design)

headers <- do.call(rbind, lapply(files$file_name, function(name) {
  frame <- flowCore::read.FCS(file.path(kDeposit, name), which.lines = 1,
                              truncate_max_range = FALSE, emptyValue = FALSE)
  keywords <- flowCore::keyword(frame)
  parameters <- flowCore::pData(flowCore::parameters(frame))
  data.frame(
    file_name = name,
    total_events = as.numeric(keywords[["$TOT"]]),
    parameters = nrow(parameters),
    named_markers = sum(nzchar(trimws(parameters$desc)) &
                          trimws(parameters$desc) != "-", na.rm = TRUE),
    compensated_channels = sum(grepl("^FJComp", parameters$name)),
    cytometer = trimws(as.character(keywords[["$CYT"]])),
    stringsAsFactors = FALSE
  )
}))
headers <- merge(headers, files, by = "file_name")
Write(headers, "file_headers.csv")

panel <- flowCore::pData(flowCore::parameters(
  flowCore::read.FCS(file.path(kDeposit, files$file_name[1]), which.lines = 1,
                     truncate_max_range = FALSE, emptyValue = FALSE)
))
panel <- data.frame(channel = panel$name, marker = panel$desc,
                    stringsAsFactors = FALSE)
Write(panel, "panel.csv")

Say("  files: ", nrow(files), ", cytometer: ",
    paste(unique(headers$cytometer), collapse = ", "))
Say("  events per file: min ", min(headers$total_events), ", median ",
    stats::median(headers$total_events), ", max ", max(headers$total_events))
Say("  files above 500,000 events: ", sum(headers$total_events > 5e5))
Say("  files whose channels are all compensated: ",
    sum(headers$compensated_channels == headers$parameters -
          sum(grepl("^(FSC|SSC|Time)", panel$channel))))

# ---------------------------------------------------------------------------
# Part 2. The manual arm, from the deposited workspace.
# ---------------------------------------------------------------------------

Say("\nPart 2: the deposited workspace")

manual <- ReadCovidManualGates(kWorkspace, kDeposit)
Write(manual, "manual_counts.csv")
Say("  populations in the workspace: ",
    paste(sort(unique(manual$population)), collapse = ", "))
Say("  samples the workspace gated: ", length(unique(manual$sample)))

Wide <- function(frame) {
  roots <- frame[frame$population == "root", c("sample", "count")]
  names(roots)[2] <- "root_events"
  out <- roots
  for (population in setdiff(unique(frame$population), "root")) {
    piece <- frame[frame$population == population, c("sample", "count")]
    names(piece)[2] <- population
    out <- merge(out, piece, by = "sample", all.x = TRUE)
  }
  out
}

manual_wide <- Wide(manual)
manual_wide$CD3_percent <- 100 * manual_wide$`/CD3+` / manual_wide$root_events
manual_wide$CD4_percent <- 100 * manual_wide$`/CD3+/CD4` / manual_wide$`/CD3+`
manual_wide$CD8_percent <- 100 * manual_wide$`/CD3+/CD8` / manual_wide$`/CD3+`
manual_wide$CD14_percent <- 100 * manual_wide$`/CD14 monocytes` /
  manual_wide$root_events
manual_wide$CD19_percent <- 100 * manual_wide$`/CD19` / manual_wide$root_events
manual_wide$arm <- "manual"

# The workspace names a sample by the file name plus the root event count, so the
# trailing count is removed before the join.
#
# Six deposited files carry `ICU_changedW` where the workspace still carries `W`.
# Those patients moved from the ward to intensive care after the workspace was
# built, and the file name records the move while the workspace does not. The
# join key removes the marker, so the two arms still meet on those six files.
JoinKey <- function(x) sub("ICU_changedW", "W", x, fixed = TRUE)

manual_wide$file_name <- sub("_[0-9]+$", "", manual_wide$sample)
manual_wide$join_key <- JoinKey(manual_wide$file_name)
Write(manual_wide, "manual_frequencies.csv")

files$join_key <- JoinKey(files$file_name)
relabelled <- files[grepl("ICU_changedW", files$file_name, fixed = TRUE), ]
Write(relabelled, "relabelled_patients.csv")
Say("  patients relabelled from ward to intensive care: ", nrow(relabelled))

matched <- sum(manual_wide$join_key %in% files$join_key)
Say("  workspace samples that match a deposited file: ", matched, " of ",
    nrow(manual_wide))

# ---------------------------------------------------------------------------
# Part 3. The automated arm.
# ---------------------------------------------------------------------------

Say("\nPart 3: gating ", nrow(files), " files")

results <- lapply(seq_len(nrow(files)), function(index) {
  GateCovidFile(file.path(kDeposit, files$file_name[index]))
})
failed <- vapply(results, function(x) !is.na(x$error), logical(1))
Write(data.frame(
  file_name = files$file_name[failed],
  group = as.character(files$group[failed]),
  reason = vapply(results[failed], function(x) x$error, character(1)),
  stringsAsFactors = FALSE
), "gating_failures.csv")
Say("  files gated: ", sum(!failed), " of ", nrow(files))

# A file whose CD3 gate holds too few events returns fewer columns, so the rows
# are padded to a common set before they are stacked.
BindRows <- function(frames) {
  columns <- unique(unlist(lapply(frames, names)))
  do.call(rbind, lapply(frames, function(frame) {
    missing <- setdiff(columns, names(frame))
    for (name in missing) {
      frame[[name]] <- NA
    }
    frame[, columns, drop = FALSE]
  }))
}

automated <- BindRows(lapply(results[!failed], function(x) x$counts))
automated$arm <- "automated"
automated <- merge(automated, files, by = "file_name")
automated$join_key <- JoinKey(automated$file_name)
Write(automated, "automated_frequencies.csv")

cuts <- do.call(rbind, lapply(results[!failed], function(x) x$cuts))
Write(cuts, "automated_cuts.csv")
cut_rules <- as.data.frame(table(marker = cuts$marker, rule = cuts$rule))
Write(cut_rules[cut_rules$Freq > 0, ], "cut_rule_summary.csv")

# ---------------------------------------------------------------------------
# Part 4. Manual against automated.
# ---------------------------------------------------------------------------

Say("\nPart 4: the two arms on the same files")

kShared <- c("CD3_percent", "CD4_percent", "CD8_percent", "CD14_percent",
             "CD19_percent")

comparison <- merge(
  manual_wide[, c("join_key", kShared)],
  automated[, c("join_key", "file_name", kShared, "group")],
  by = "join_key", suffixes = c("_manual", "_automated")
)
Write(comparison, "arm_comparison.csv")
Say("  files measured by both arms: ", nrow(comparison))

agreement <- do.call(rbind, lapply(kShared, function(measure) {
  manual_values <- comparison[[paste0(measure, "_manual")]]
  automated_values <- comparison[[paste0(measure, "_automated")]]
  usable <- is.finite(manual_values) & is.finite(automated_values)
  if (sum(usable) < 4) {
    return(NULL)
  }
  test <- suppressWarnings(stats::cor.test(
    manual_values[usable], automated_values[usable], method = "spearman"
  ))
  data.frame(
    measure = measure,
    files = sum(usable),
    manual_median = stats::median(manual_values[usable]),
    automated_median = stats::median(automated_values[usable]),
    median_difference = stats::median(
      automated_values[usable] - manual_values[usable]
    ),
    spearman_rho = unname(test$estimate),
    p_value = test$p.value,
    stringsAsFactors = FALSE
  )
}))
Write(agreement, "arm_agreement.csv")
print(agreement)

# ---------------------------------------------------------------------------
# Part 5. The claims.
# ---------------------------------------------------------------------------

Say("\nPart 5: the claims")

manual_with_group <- merge(manual_wide, files, by = "file_name")
manual_with_group$CD4_over_CD8 <- manual_with_group$CD4_percent /
  manual_with_group$CD8_percent
automated$CD4_over_CD8 <- automated$CD4_percent / automated$CD8_percent

kTested <- c("CD3_percent", "CD4_percent", "CD8_percent", "CD14_percent",
             "CD19_percent", "CD4_over_CD8")

by_group <- do.call(rbind, lapply(kTested, function(measure) {
  rows <- lapply(
    list(manual = manual_with_group, automated = automated),
    function(frame) {
      if (!measure %in% names(frame)) {
        return(NULL)
      }
      result <- CompareByGroup(frame[[measure]], as.character(frame$group))
      if (is.null(result)) {
        return(NULL)
      }
      cbind(measure = measure, arm = if (identical(frame, automated)) {
        "automated"
      } else {
        "manual"
      }, result)
    }
  )
  do.call(rbind, rows)
}))
Write(by_group, "measures_by_group.csv")
print(by_group)

kFunctional <- c("Th1_percent_of_CD4", "Tc1_percent_of_CD8")
functional <- do.call(rbind, lapply(kFunctional, function(measure) {
  result <- CompareByGroup(automated[[measure]], as.character(automated$group))
  if (is.null(result)) {
    return(NULL)
  }
  cbind(measure = measure, arm = "automated", result)
}))
Write(functional, "functional_by_group.csv")
print(functional)

pd1_columns <- grep("_PD1$", names(automated), value = TRUE)
pd1 <- do.call(rbind, lapply(pd1_columns, function(column) {
  values <- automated[[column]]
  data.frame(
    subset = sub("_PD1$", "", column),
    samples = sum(is.finite(values)),
    median_pd1_percent = stats::median(values, na.rm = TRUE),
    stringsAsFactors = FALSE
  )
}))
pd1 <- pd1[order(-pd1$median_pd1_percent), ]
Write(pd1, "pd1_by_subset.csv")
print(pd1)

Verdict <- function(id, observed, verdict) {
  data.frame(claim_id = id, observed = observed, verdict = verdict,
             stringsAsFactors = FALSE)
}

Row <- function(frame, measure, arm) {
  frame[frame$measure == measure & frame$arm == arm, ]
}

verdicts <- NULL

cd4_manual <- Row(by_group, "CD4_percent", "manual")
cd8_manual <- Row(by_group, "CD8_percent", "manual")
verdicts <- rbind(verdicts, Verdict(
  1,
  sprintf("CD4 between severity p %.3f and CD8 p %.3f, by the deposited gates",
          cd4_manual$between_severity_p, cd8_manual$between_severity_p),
  if (cd4_manual$between_severity_p > 0.05 &&
      cd8_manual$between_severity_p > 0.05) {
    "reproduced"
  } else {
    "not reproduced"
  }
))

th1 <- functional[functional$measure == "Th1_percent_of_CD4", ]
verdicts <- rbind(verdicts, Verdict(
  2,
  sprintf("Th1 healthy median %.2f against COVID-19, p %.3f",
          th1$healthy_median, th1$covid_against_healthy_p),
  if (nrow(th1) == 1 && is.finite(th1$covid_against_healthy_p) &&
      th1$covid_against_healthy_p < 0.05 &&
      th1$healthy_median > mean(c(th1$ward_median,
                                  th1$intensive_care_median), na.rm = TRUE)) {
    "reproduced"
  } else if (nrow(th1) == 1 &&
             th1$healthy_median > mean(c(th1$ward_median,
                                         th1$intensive_care_median),
                                       na.rm = TRUE)) {
    "partly reproduced"
  } else {
    "not reproduced"
  }
))

verdicts <- rbind(verdicts, Verdict(
  3,
  sprintf("the highest PD-1 subset is %s at a median of %.1f percent",
          pd1$subset[1], pd1$median_pd1_percent[1]),
  if (pd1$subset[1] == "CD8_TEMRA") {
    "reproduced"
  } else if (grepl("^CD8", pd1$subset[1])) {
    "partly reproduced"
  } else {
    "not reproduced"
  }
))

ratio_manual <- Row(by_group, "CD4_over_CD8", "manual")
verdicts <- rbind(verdicts, Verdict(
  4,
  sprintf("the CD4 over CD8 ratio is %.2f healthy, %.2f ward and %.2f intensive care, p %.3f against healthy",
          ratio_manual$healthy_median, ratio_manual$ward_median,
          ratio_manual$intensive_care_median,
          ratio_manual$covid_against_healthy_p),
  if (ratio_manual$intensive_care_median > ratio_manual$healthy_median &&
      ratio_manual$ward_median > ratio_manual$healthy_median) {
    "reproduced"
  } else if (ratio_manual$intensive_care_median >
             ratio_manual$healthy_median) {
    "partly reproduced"
  } else {
    "not reproduced"
  }
))

verdicts <- rbind(verdicts, Verdict(
  5,
  sprintf("the largest file holds %s events and %d files exceed 500,000",
          format(max(headers$total_events), big.mark = ","),
          sum(headers$total_events > 5e5)),
  if (sum(headers$total_events > 5e5) == 0) "reproduced" else "not reproduced"
))

fluorescence <- sum(!grepl("^(FSC|SSC|Time)", panel$channel))
verdicts <- rbind(verdicts, Verdict(
  6,
  sprintf("%d of %d fluorescence channels are named FJComp in every file",
          min(headers$compensated_channels), fluorescence),
  if (all(headers$compensated_channels == fluorescence)) {
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

long <- do.call(rbind, lapply(kShared, function(measure) {
  data.frame(
    measure = measure,
    manual = comparison[[paste0(measure, "_manual")]],
    automated = comparison[[paste0(measure, "_automated")]],
    group = comparison$group,
    stringsAsFactors = FALSE
  )
}))
Save(
  ggplot(long, aes(manual, automated, colour = group)) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
    geom_point(size = 2) +
    facet_wrap(~measure, scales = "free") +
    labs(x = "The deposited gates, percent", y = "The automated rule, percent",
         colour = NULL) +
    ThemePublication(),
  "arm_comparison.svg", width = 10, height = 7
)

group_long <- do.call(rbind, lapply(kTested, function(measure) {
  rbind(
    data.frame(measure = measure, arm = "manual",
               value = manual_with_group[[measure]],
               group = as.character(manual_with_group$group),
               stringsAsFactors = FALSE),
    data.frame(measure = measure, arm = "automated",
               value = automated[[measure]],
               group = as.character(automated$group),
               stringsAsFactors = FALSE)
  )
}))
group_long$group <- factor(group_long$group, levels = kGroupOrder)
Save(
  ggplot(group_long[is.finite(group_long$value), ],
         aes(group, value, fill = arm)) +
    geom_boxplot(outlier.size = 0.6) +
    facet_wrap(~measure, scales = "free_y") +
    scale_fill_manual(values = c(manual = "#2166AC", automated = "#1B7837")) +
    labs(x = NULL, y = "Percent", fill = NULL) +
    ThemePublication() +
    theme(axis.text.x = element_text(angle = 20, hjust = 1)),
  "measures_by_group.svg", width = 10, height = 7
)

pd1_long <- do.call(rbind, lapply(pd1_columns, function(column) {
  data.frame(subset = sub("_PD1$", "", column), value = automated[[column]],
             stringsAsFactors = FALSE)
}))
Save(
  ggplot(pd1_long[is.finite(pd1_long$value), ],
         aes(stats::reorder(subset, value, stats::median), value)) +
    geom_boxplot(fill = "#B2182B") +
    coord_flip() +
    labs(x = NULL, y = "PD-1 positive, percent of the subset") +
    ThemePublication(),
  "pd1_by_subset.svg", width = 8, height = 5
)

Say("\nDone. Every table is in ", kOutputDir)
