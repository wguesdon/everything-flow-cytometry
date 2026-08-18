# Human bone marrow by flow cytometry and by mass cytometry, Oetjen 2018.
#
# FR-FCM-ZYQ9 holds 132 flow cytometry files, which is five 13 colour panels and
# one unstained control for each of 22 samples from 20 donors. FR-FCM-ZYQB holds
# a 49 marker Helios panel for eight of those donors.
#
# The entry into the analysis is the depositors' own scatter gate and their own
# viability gate, both read from the cytoflow workflow they shipped. Everything
# below it is the automated cut rule this repository uses elsewhere.
#
# Run it in the container:
#   podman run --rm -v "$PWD:/work:z" -w /work everything-flow-cytometry:latest
# \
#     Rscript scripts/10_oetjen2018_bone_marrow.R

suppressPackageStartupMessages({
  library(flowCore)
  library(ggplot2)
  library(withr)
})

source(file.path("R", "figures.R"))
source(file.path("R", "harmonisation.R"))
source(file.path("R", "spectral.R"))
source(file.path("R", "naive_memory.R"))
source(file.path("R", "bone_marrow.R"))

kFlowDir <- file.path("data", "datasets", "flowrepository", "FR-FCM-ZYQ9")
kMassDir <- file.path("data", "datasets", "flowrepository", "FR-FCM-ZYQB")
kGatingDir <- "gating"
kOutputDir <- file.path("output", "oetjen2018")

kStainedPanels <- c("T", "B", "NK", "Mono", "DC")

dir.create(kOutputDir, recursive = TRUE, showWarnings = FALSE)

Write <- function(x, name) {
  utils::write.csv(x, file.path(kOutputDir, name), row.names = FALSE)
  invisible(x)
}
Say <- function(...) cat(..., "\n", sep = "")

# ---------------------------------------------------------------------------
# Part 1. The deposit and the cohort.
# ---------------------------------------------------------------------------

Say("Part 1: the deposit")

donors <- ReadBoneMarrowDonors(file.path(kGatingDir,
                                         "oetjen2018_donor_metadata.csv"))
panels <- ReadBoneMarrowPanels(file.path(kGatingDir, "oetjen2018_panels.csv"))
gates <- ReadBoneMarrowGates(file.path(kGatingDir,
                                       "oetjen2018_manual_gates.csv"))
claims <- utils::read.csv(
  file.path(kGatingDir, "oetjen2018_paper_claims.csv"), stringsAsFactors = FALSE
)

files <- ParseBoneMarrowFileNames(list.files(kFlowDir, pattern = "\\.fcs$"))
Write(files, "file_index.csv")

coverage <- as.data.frame(table(panel = files$panel))
Write(coverage, "panel_coverage.csv")
print(coverage)

unique_donors <- donors[!duplicated(donors$donor), ]
cohort <- data.frame(
  Measure = c("Samples in Table 1", "Distinct donors", "Female donors",
              "Male donors", "Youngest age", "Oldest age",
              "Median age of the donors", "Median age of the samples",
              "Donors with mass cytometry"),
  Value = c(
    nrow(donors), nrow(unique_donors),
    sum(unique_donors$sex == "F"), sum(unique_donors$sex == "M"),
    min(donors$age), max(donors$age),
    stats::median(unique_donors$age), stats::median(donors$age),
    length(unique(donors$donor[donors$mass_cytometry == "yes"]))
  ),
  stringsAsFactors = FALSE
)
Write(cohort, "cohort.csv")
print(cohort)

Say("  samples in the deposit but not in Table 1: ",
    paste(setdiff(unique(files$sample), donors$sample), collapse = ", "))
Say("  samples in Table 1 but not in the deposit: ",
    paste(setdiff(donors$sample, unique(files$sample)), collapse = ", "))

# ---------------------------------------------------------------------------
# Part 2. Gate every stained file behind the depositors' own entry gates.
# ---------------------------------------------------------------------------

Say("\nPart 2: gating ", sum(files$panel %in% kStainedPanels), " stained files")

stained <- files[files$panel %in% kStainedPanels, ]
results <- vector("list", nrow(stained))
started <- Sys.time()
for (index in seq_len(nrow(stained))) {
  results[[index]] <- GateBoneMarrowFile(
    file.path(kFlowDir, stained$file_name[index]), panels,
    stained$panel[index], gates
  )
  if (index %% 20 == 0 || index == nrow(stained)) {
    Say("  ", index, " of ", nrow(stained), ", ",
        round(as.numeric(difftime(Sys.time(), started, units = "mins")), 1),
        " minutes")
  }
}

failed <- vapply(results, function(x) !is.na(x$error), logical(1))
failures <- data.frame(
  file_name = stained$file_name[failed],
  panel = stained$panel[failed],
  sample = stained$sample[failed],
  reason = vapply(results[failed], function(x) x$error, character(1)),
  stringsAsFactors = FALSE
)
Write(failures, "gating_failures.csv")
Say("  files gated: ", sum(!failed), " of ", nrow(stained))

counts <- do.call(rbind, lapply(results[!failed], function(x) x$counts))
counts <- merge(counts, stained[, c("file_name", "sample")], by = "file_name")
counts <- merge(counts, donors[, c("sample", "donor", "sex", "age")],
                by = "sample", all.x = TRUE)
Write(counts, "gating_counts.csv")

markers <- do.call(rbind, lapply(results[!failed], function(x) x$markers))
markers <- merge(markers, stained[, c("file_name", "sample")], by = "file_name")
Write(markers, "marker_cuts.csv")

resolution <- do.call(rbind, lapply(results[!failed], function(x) x$channels))
resolution_summary <- as.data.frame(
  table(panel = sub("_.*", "", resolution$file_name),
        resolved_by = resolution$resolved_by)
)
Write(as.data.frame(table(resolved_by = resolution$resolved_by)),
      "channel_resolution_summary.csv")
Say("  channels resolved by the typed marker name: ",
    sum(resolution$resolved_by == "marker"),
    ", by the published detector: ", sum(resolution$resolved_by == "detector"))
Say("  files that needed the detector fallback: ",
    length(unique(resolution$file_name[resolution$resolved_by == "detector"])))

# How many CD45 positive events each panel collected. The T panel carries no
# CD45, so it is reported separately.
with_cd45 <- counts[counts$panel != "T", ]
event_summary <- data.frame(
  Measure = c("Panels that carry CD45", "Files with a CD45 gate",
              "Median CD45 positive events before the scatter gate",
              "25th percentile", "75th percentile",
              "Median events in the T panel entry gate"),
  Value = c(
    length(unique(with_cd45$panel)), nrow(with_cd45),
    round(stats::median(with_cd45$cd45_before_scatter, na.rm = TRUE)),
    round(stats::quantile(with_cd45$cd45_before_scatter, 0.25, na.rm = TRUE)),
    round(stats::quantile(with_cd45$cd45_before_scatter, 0.75, na.rm = TRUE)),
    round(stats::median(counts$parent_events[counts$panel == "T"]))
  ),
  stringsAsFactors = FALSE
)
Write(event_summary, "event_summary.csv")
print(event_summary)

hierarchy <- do.call(rbind, lapply(split(counts, counts$panel), function(piece) {
  data.frame(
    panel = piece$panel[1],
    files = nrow(piece),
    median_total = round(stats::median(piece$total_events)),
    median_singlets = round(stats::median(piece$singlet_events)),
    median_after_scatter = round(stats::median(piece$scatter_events)),
    median_live = round(stats::median(piece$live_events)),
    median_parent = round(stats::median(piece$parent_events)),
    stringsAsFactors = FALSE
  )
}))
Write(hierarchy, "hierarchy.csv")
print(hierarchy)

# ---------------------------------------------------------------------------
# Part 3. Whether an automated cut can be trusted on this deposit.
# ---------------------------------------------------------------------------

Say("\nPart 3: the spread of every marker across the samples")

stability <- do.call(rbind, lapply(
  split(markers, list(markers$panel, markers$marker), drop = TRUE),
  function(piece) {
    values <- piece$percent_positive[is.finite(piece$percent_positive)]
    if (length(values) < 3) {
      return(NULL)
    }
    data.frame(
      panel = piece$panel[1],
      marker = piece$marker[1],
      samples = length(values),
      median_percent = stats::median(values),
      minimum_percent = min(values),
      maximum_percent = max(values),
      range_percent = max(values) - min(values),
      cv = CoefficientOfVariation(values),
      stringsAsFactors = FALSE
    )
  }
))
stability <- stability[order(-stability$range_percent), ]
Write(stability, "marker_stability.csv")
print(utils::head(stability, 12))

verdict <- data.frame(
  Measure = c(
    "Marker and panel pairs measured",
    "Pairs whose percent positive spans more than 50 points",
    "Pairs whose percent positive spans more than 80 points",
    "Pairs with a coefficient of variation above 1",
    "Median range across all pairs, in percentage points"
  ),
  Value = c(
    nrow(stability),
    sum(stability$range_percent > 50),
    sum(stability$range_percent > 80),
    sum(stability$cv > 1, na.rm = TRUE),
    round(stats::median(stability$range_percent), 1)
  ),
  stringsAsFactors = FALSE
)
Write(verdict, "stability_verdict.csv")
print(verdict)

# ---------------------------------------------------------------------------
# Part 4. T cell subsets, by flow cytometry and by mass cytometry.
# ---------------------------------------------------------------------------

Say("\nPart 4: T cell subsets by two technologies")

t_files <- files[files$panel == "T", ]
t_results <- lapply(seq_len(nrow(t_files)), function(index) {
  GateTcellSubsets(file.path(kFlowDir, t_files$file_name[index]), panels, gates)
})
t_failed <- vapply(t_results, function(x) !is.na(x$error), logical(1))
Write(data.frame(
  file_name = t_files$file_name[t_failed],
  sample = t_files$sample[t_failed],
  reason = vapply(t_results[t_failed], function(x) x$error, character(1)),
  stringsAsFactors = FALSE
), "t_panel_failures.csv")
flow_subsets <- do.call(rbind, lapply(which(!t_failed), function(index) {
  cbind(sample = t_files$sample[index], t_results[[index]]$counts)
}))
Say("  flow T panel files gated: ", nrow(flow_subsets), " of ", nrow(t_files))
if (any(t_failed)) {
  Say("  reasons: ", paste(unique(vapply(t_results[t_failed],
                                         function(x) x$error, character(1))),
                           collapse = " | "))
}

mass_files <- list.files(kMassDir, pattern = "\\.fcs$")
mass_subsets <- do.call(rbind, lapply(mass_files, function(name) {
  result <- GateMassCytometryFile(file.path(kMassDir, name))
  if (!is.na(result$error)) {
    Say("  mass file failed: ", name, " - ", result$error)
    return(NULL)
  }
  cbind(sample = sub("^[0-9]+_([A-Za-z]+)_.*$", "\\1", name), result$counts)
}))
Say("  mass cytometry files gated: ", nrow(mass_subsets), " of ",
    length(mass_files))

subsets <- rbind(
  flow_subsets[, intersect(names(flow_subsets), names(mass_subsets))],
  mass_subsets[, intersect(names(flow_subsets), names(mass_subsets))]
)
subsets <- merge(subsets, donors[, c("sample", "donor", "sex", "age")],
                 by = "sample", all.x = TRUE)
Write(subsets, "t_cell_subsets.csv")

kSubsetColumns <- c(
  "CD4_N", "CD4_CM", "CD4_EM", "CD4_TEMRA",
  "CD8_N", "CD8_CM", "CD8_EM", "CD8_TEMRA",
  "CD4_percent_of_CD3", "CD8_percent_of_CD3"
)

paired <- merge(
  subsets[subsets$technology == "flow", c("sample", kSubsetColumns)],
  subsets[subsets$technology == "mass", c("sample", kSubsetColumns)],
  by = "sample", suffixes = c("_flow", "_mass")
)
Write(paired, "flow_against_mass.csv")
Say("  donors measured by both technologies: ", nrow(paired))

technology_agreement <- do.call(rbind, lapply(kSubsetColumns, function(column) {
  flow <- paired[[paste0(column, "_flow")]]
  mass <- paired[[paste0(column, "_mass")]]
  usable <- is.finite(flow) & is.finite(mass)
  if (sum(usable) < 4) {
    return(NULL)
  }
  test <- suppressWarnings(
    stats::cor.test(flow[usable], mass[usable], method = "spearman")
  )
  data.frame(
    measure = column,
    donors = sum(usable),
    flow_median = stats::median(flow[usable]),
    mass_median = stats::median(mass[usable]),
    spearman_rho = unname(test$estimate),
    p_value = test$p.value,
    stringsAsFactors = FALSE
  )
}))
Write(technology_agreement, "technology_agreement.csv")
print(technology_agreement)

# ---------------------------------------------------------------------------
# Part 5. The age trend and the second aspiration.
# ---------------------------------------------------------------------------

Say("\nPart 5: age and the repeated aspiration")

flow_only <- subsets[subsets$technology == "flow" & !is.na(subsets$age), ]
flow_only$CD4_memory <- flow_only$CD4_CM + flow_only$CD4_EM
flow_only$CD8_memory <- flow_only$CD8_CM + flow_only$CD8_EM

age_trend <- do.call(rbind, lapply(
  c("CD4_N", "CD4_memory", "CD4_TEMRA", "CD8_N", "CD8_memory", "CD8_TEMRA"),
  function(column) {
    values <- flow_only[[column]]
    usable <- is.finite(values)
    if (sum(usable) < 6) {
      return(NULL)
    }
    test <- suppressWarnings(stats::cor.test(
      values[usable], flow_only$age[usable], method = "spearman"
    ))
    data.frame(
      measure = column,
      samples = sum(usable),
      spearman_rho = unname(test$estimate),
      p_value = test$p.value,
      stringsAsFactors = FALSE
    )
  }
))
Write(age_trend, "age_trend.csv")
print(age_trend)

repeated <- flow_only[flow_only$donor %in%
                        names(which(table(flow_only$donor) > 1)), ]
repeat_agreement <- do.call(rbind, lapply(
  split(repeated, repeated$donor),
  function(piece) {
    if (nrow(piece) < 2) {
      return(NULL)
    }
    do.call(rbind, lapply(kSubsetColumns, function(column) {
      values <- piece[[column]]
      data.frame(
        donor = piece$donor[1],
        measure = column,
        first = values[1],
        second = values[2],
        difference = abs(values[1] - values[2]),
        stringsAsFactors = FALSE
      )
    }))
  }
))
Write(repeat_agreement, "repeated_aspiration.csv")
if (!is.null(repeat_agreement)) {
  Say("  donors with two aspirations: ",
      length(unique(repeat_agreement$donor)),
      ", median absolute difference: ",
      round(stats::median(repeat_agreement$difference, na.rm = TRUE), 1),
      " percentage points")
}

# ---------------------------------------------------------------------------
# Part 6. The claims.
# ---------------------------------------------------------------------------

Say("\nPart 6: the claims")

Verdict <- function(id, observed, verdict) {
  data.frame(claim_id = id, observed = observed, verdict = verdict,
             stringsAsFactors = FALSE)
}

verdicts <- rbind(
  Verdict(1, sprintf(
    "%d male and %d female donors, ages %.0f to %.0f, median %.1f",
    sum(unique_donors$sex == "M"), sum(unique_donors$sex == "F"),
    min(donors$age), max(donors$age), stats::median(unique_donors$age)
  ), if (sum(unique_donors$sex == "M") == 10 &&
         sum(unique_donors$sex == "F") == 10 &&
         min(donors$age) == 24 && max(donors$age) == 84 &&
         abs(stats::median(unique_donors$age) - 57) <= 1) {
    "reproduced"
  } else {
    "partly reproduced"
  }),
  Verdict(2, sprintf(
    "median %s CD45 positive events before the scatter gate, quartiles %s and %s",
    format(round(stats::median(with_cd45$cd45_before_scatter, na.rm = TRUE)),
           big.mark = ","),
    format(round(stats::quantile(with_cd45$cd45_before_scatter, 0.25,
                                 na.rm = TRUE)), big.mark = ","),
    format(round(stats::quantile(with_cd45$cd45_before_scatter, 0.75,
                                 na.rm = TRUE)), big.mark = ",")
  ), if (abs(stats::median(with_cd45$cd45_before_scatter, na.rm = TRUE) -
             196000) < 100000) {
    "reproduced"
  } else {
    "not reproduced"
  }),
  Verdict(8, sprintf(
    "%d stained files and %d unstained controls for %d samples",
    sum(files$panel %in% kStainedPanels), sum(files$panel == "Unstained"),
    length(unique(files$sample))
  ), if (sum(files$panel %in% kStainedPanels) == 110 &&
         length(unique(files$sample)) == 22) {
    "reproduced"
  } else {
    "partly reproduced"
  })
)

# Claims 3 and 7 need a reliable lineage frequency from every panel. Part 3
# measures whether an automated cut can supply one.
unstable <- sum(stability$range_percent > 80)
verdicts <- rbind(verdicts, Verdict(
  3,
  sprintf("%d of %d marker and panel pairs span more than 80 percentage points across the samples",
          unstable, nrow(stability)),
  "not measurable"
))
verdicts <- rbind(verdicts, Verdict(
  7,
  "The lineage frequency of a panel is not stable enough to rank the populations",
  "not measurable"
))

if (!is.null(technology_agreement) && nrow(technology_agreement) > 0) {
  strong <- sum(technology_agreement$spearman_rho > 0.5, na.rm = TRUE)
  verdicts <- rbind(verdicts, Verdict(
    6,
    sprintf("%d of %d subset measures correlate above rho 0.5 across %d donors",
            strong, nrow(technology_agreement), max(technology_agreement$donors)),
    if (strong >= nrow(technology_agreement) / 2) {
      "reproduced"
    } else {
      "not reproduced"
    }
  ))
} else {
  verdicts <- rbind(verdicts, Verdict(
    6, "No donor was gated by both technologies", "not measurable"
  ))
}

if (!is.null(age_trend) && nrow(age_trend) > 0) {
  memory_rows <- age_trend[grepl("memory", age_trend$measure), ]
  verdicts <- rbind(verdicts, Verdict(
    4,
    paste(sprintf("%s rho %.2f p %.3f", memory_rows$measure,
                  memory_rows$spearman_rho, memory_rows$p_value),
          collapse = "; "),
    if (any(memory_rows$spearman_rho > 0 & memory_rows$p_value < 0.05)) {
      "reproduced"
    } else if (any(memory_rows$spearman_rho > 0)) {
      "partly reproduced"
    } else {
      "not reproduced"
    }
  ))
} else {
  verdicts <- rbind(verdicts, Verdict(4, "No usable T cell subsets",
                                      "not measurable"))
}

if (!is.null(repeat_agreement) && nrow(repeat_agreement) > 0) {
  verdicts <- rbind(verdicts, Verdict(
    5,
    sprintf("median absolute difference of %.1f percentage points across %d measures",
            stats::median(repeat_agreement$difference, na.rm = TRUE),
            nrow(repeat_agreement)),
    if (stats::median(repeat_agreement$difference, na.rm = TRUE) < 10) {
      "reproduced"
    } else {
      "not reproduced"
    }
  ))
} else {
  verdicts <- rbind(verdicts, Verdict(5, "No donor was gated twice",
                                      "not measurable"))
}

claim_verdicts <- merge(claims, verdicts, by = "claim_id", all.x = TRUE)
claim_verdicts <- claim_verdicts[order(claim_verdicts$claim_id), ]
Write(claim_verdicts, "claim_verdicts.csv")
print(claim_verdicts[, c("claim_id", "short_name", "expected", "observed",
                         "verdict")])
Say("\n  verdicts: ",
    paste(names(table(claim_verdicts$verdict)),
          table(claim_verdicts$verdict), sep = " = ", collapse = ", "))

# ---------------------------------------------------------------------------
# Part 7. Figures.
# ---------------------------------------------------------------------------

Say("\nPart 7: figures")

Save <- function(plot, name, width = 9, height = 6) {
  SaveFigure(plot, file.path(kOutputDir, name), width = width, height = height)
}

Save(
  ggplot(markers, aes(stats::reorder(marker, percent_positive, stats::median,
                                     na.rm = TRUE),
                      percent_positive)) +
    geom_boxplot(outlier.size = 0.6, fill = "#4393C3") +
    facet_wrap(~panel, scales = "free_y") +
    coord_flip() +
    labs(x = NULL, y = "Percent positive of the parent gate") +
    ThemePublication(),
  "marker_spread.svg", width = 11, height = 8
)

if (nrow(paired) > 0) {
  long <- do.call(rbind, lapply(kSubsetColumns, function(column) {
    data.frame(
      measure = column,
      flow = paired[[paste0(column, "_flow")]],
      mass = paired[[paste0(column, "_mass")]],
      stringsAsFactors = FALSE
    )
  }))
  Save(
    ggplot(long, aes(flow, mass)) +
      geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
      geom_point(size = 2, colour = "#B2182B") +
      facet_wrap(~measure, scales = "free") +
      labs(x = "Flow cytometry, percent", y = "Mass cytometry, percent") +
      ThemePublication(),
    "flow_against_mass.svg", width = 10, height = 7
  )
}

if (nrow(flow_only) > 0) {
  age_long <- do.call(rbind, lapply(
    c("CD4_N", "CD4_memory", "CD4_TEMRA", "CD8_N", "CD8_memory", "CD8_TEMRA"),
    function(column) {
      data.frame(measure = column, age = flow_only$age,
                 percent = flow_only[[column]], stringsAsFactors = FALSE)
    }
  ))
  Save(
    ggplot(age_long, aes(age, percent)) +
      geom_point(size = 2, colour = "#1B7837") +
      geom_smooth(method = "lm", formula = y ~ x, se = FALSE,
                  colour = "#762A83") +
      facet_wrap(~measure, scales = "free_y") +
      labs(x = "Donor age, years", y = "Percent of the T cell subset") +
      ThemePublication(),
    "age_trend.svg", width = 10, height = 6
  )
}

Say("\nDone. Every table is in ", kOutputDir)
