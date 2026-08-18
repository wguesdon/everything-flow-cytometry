# Test the claims of Yu 2021 against the data the authors deposited.
#
# The paper states its findings as directions between groups, not as effect
# sizes.
# "Healthy females had greater frequencies of CD161hi cells relative to males"
# is a
# direction. The claims live in gating/yu2021_paper_claims.csv with the sentence
# each one came from, so a reader can check the claim against the paper before
# checking the number against the data.
#
# Two columns are reported for every claim and they answer different questions.
# `verdict` compares the observed direction against the stated direction.
# `p_value`
# is the test the paper itself used, which is Mann-Whitney for a two group
# comparison and Spearman for a correlation. A direction can reproduce while the
# test is not significant, and reporting only one of the two hides that.

#' Read the claims a paper makes, with the sentence each one came from
#'
#' @param path Path to a CSV with the columns `claim_id`, `short_name`,
#' `measure`,
#'   `test`, `expected`, `quote` and `figure`.
#' @return A `data.frame` of claims.
#' @export
ReadYuClaims <- function(path) {
  if (!file.exists(path)) {
    stop("The claims file does not exist: ", path)
  }

  claims <- utils::read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
  required <- c("claim_id", "short_name", "measure", "test", "expected",
                "quote", "figure")
  missing <- setdiff(required, colnames(claims))
  if (length(missing) > 0) {
    stop("The claims file is missing the column(s): ",
         paste(missing, collapse = ", "), ".")
  }

  claims
}

#' Compare one measure between the sexes
#'
#' The paper used a Mann-Whitney test for every two group comparison of a
#' frequency, so this uses the same test.
#'
#' @param values A numeric vector of the measure, one entry per sample.
#' @param sex A character vector of `"Female"` and `"Male"`, the same length.
#' @return A one row `data.frame` with `n_female`, `n_male`, `female_median`,
#'   `male_median`, `difference`, `observed`, and `p_value`. `observed` is
#'   `"female higher"` or `"male higher"`, and it is `NA` when either group is
#'   empty.
#' @examples
#' CompareBySex(c(4, 5, 2, 1), c("Female", "Female", "Male", "Male"))
#' @export
CompareBySex <- function(values, sex) {
  usable <- !is.na(values) & !is.na(sex)
  values <- values[usable]
  sex <- sex[usable]

  female <- values[sex == "Female"]
  male <- values[sex == "Male"]

  if (length(female) == 0 || length(male) == 0) {
    return(data.frame(
      n_female = length(female), n_male = length(male),
      female_median = NA_real_, male_median = NA_real_,
      difference = NA_real_, observed = NA_character_, p_value = NA_real_,
      stringsAsFactors = FALSE
    ))
  }

  female_median <- stats::median(female)
  male_median <- stats::median(male)
  p_value <- NA_real_
  if (length(female) >= 2 && length(male) >= 2) {
    test <- try(stats::wilcox.test(female, male, exact = FALSE), silent = TRUE)
    if (!methods::is(test, "try-error")) {
      p_value <- test$p.value
    }
  }

  data.frame(
    n_female = length(female), n_male = length(male),
    female_median = female_median, male_median = male_median,
    difference = female_median - male_median,
    observed = if (female_median >= male_median) "female higher" else "male higher",
    p_value = p_value,
    stringsAsFactors = FALSE
  )
}

#' Correlate one measure against the severity rank
#'
#' @param values A numeric vector of the measure.
#' @param severity_index An integer severity rank, 1 for `normal` up to 4 for
#'   `hospitalized`.
#' @return A one row `data.frame` with `n`, `rho`, `p_value` and `observed`,
#' where
#'   `observed` is `"falls with severity"` or `"rises with severity"`.
#' @export
CorrelateWithSeverity <- function(values, severity_index) {
  usable <- !is.na(values) & !is.na(severity_index)
  values <- values[usable]
  severity_index <- severity_index[usable]

  if (length(values) < 3) {
    return(data.frame(n = length(values), rho = NA_real_, p_value = NA_real_,
                      observed = NA_character_, stringsAsFactors = FALSE))
  }

  test <- suppressWarnings(
    stats::cor.test(values, severity_index, method = "spearman")
  )

  data.frame(
    n = length(values),
    rho = unname(test$estimate),
    p_value = test$p.value,
    observed = if (test$estimate < 0) "falls with severity" else "rises with severity",
    stringsAsFactors = FALSE
  )
}

#' Fit one slope against the severity rank for each sex
#'
#' Figure 2C of the paper is a linear regression of the CD161hi frequency
#' against
#' the severity rank, drawn once per sex. The claim is about which slope is
#' steeper downwards.
#'
#' @param values A numeric vector of the measure.
#' @param severity_index An integer severity rank.
#' @param sex A character vector of `"Female"` and `"Male"`.
#' @return A one row `data.frame` with `female_slope`, `male_slope`,
#'   `slope_difference` and `observed`. `observed` is
#'   `"females fall faster"` or `"males fall faster"`.
#' @export
CompareSlopesBySex <- function(values, severity_index, sex) {
  SlopeFor <- function(which_sex) {
    keep <- sex == which_sex & !is.na(values) & !is.na(severity_index)
    if (sum(keep) < 3 || length(unique(severity_index[keep])) < 2) {
      return(NA_real_)
    }
    unname(stats::coef(stats::lm(values[keep] ~ severity_index[keep]))[2])
  }

  female_slope <- SlopeFor("Female")
  male_slope <- SlopeFor("Male")

  observed <- NA_character_
  if (!is.na(female_slope) && !is.na(male_slope)) {
    observed <- if (female_slope < male_slope) {
      "females fall faster"
    } else {
      "males fall faster"
    }
  }

  data.frame(
    female_slope = female_slope, male_slope = male_slope,
    slope_difference = female_slope - male_slope,
    observed = observed, stringsAsFactors = FALSE
  )
}

#' Test every claim of the paper against the measured frequencies
#'
#' @param claims The output of [ReadYuClaims()].
#' @param measures A `data.frame` with one row per sample, holding the columns
#' the
#'   claims name in their `measure` column, joined to the sample sheet.
#' @param alpha The significance level used to judge claim 9, which states that
#' a
#'   measure does not change.
#' @return A `data.frame` with one row per claim, holding the claim, the group
#'   medians, the test statistic, the p value and a `verdict` of `"reproduced"`,
#'   `"opposite"`, `"partly reproduced"` or `"not measured"`.
#' @export
TestYuClaims <- function(claims, measures, alpha = 0.05) {
  infected <- measures[measures$severity_rank != "normal", , drop = FALSE]

  Row <- function(claim, detail, observed, statistic, p_value, verdict) {
    data.frame(
      claim_id = claim$claim_id,
      short_name = claim$short_name,
      measure = claim$measure,
      figure = claim$figure,
      expected = claim$expected,
      detail = detail,
      observed = observed,
      statistic = statistic,
      p_value = p_value,
      verdict = verdict,
      stringsAsFactors = FALSE
    )
  }

  rows <- lapply(seq_len(nrow(claims)), function(i) {
    claim <- claims[i, ]
    values <- measures[[claim$measure]]

    if (is.null(values) && claim$test != "positive_control") {
      return(Row(claim, "the measure is not in the data", NA_character_,
                 NA_real_, NA_real_, "not measured"))
    }

    if (claim$test == "spearman_vs_severity") {
      result <- CorrelateWithSeverity(values, measures$severity_index)
      verdict <- if (is.na(result$observed)) {
        "not measured"
      } else if (result$observed == "falls with severity") {
        "reproduced"
      } else {
        "opposite"
      }
      return(Row(claim, sprintf("n = %d samples", result$n), result$observed,
                 result$rho, result$p_value, verdict))
    }

    if (claim$test == "sex_difference_in_normal") {
      keep <- measures$severity_rank == "normal"
      result <- CompareBySex(values[keep], measures$sex[keep])
      wanted <- if (grepl("^female", claim$expected)) {
        "female higher"
      } else {
        "male higher"
      }
      verdict <- if (is.na(result$observed)) {
        "not measured"
      } else if (result$observed == wanted) {
        "reproduced"
      } else {
        "opposite"
      }
      return(Row(
        claim,
        sprintf("female median %.2f against male median %.2f, n = %d and %d",
                result$female_median, result$male_median,
                result$n_female, result$n_male),
        result$observed, result$difference, result$p_value, verdict
      ))
    }

    if (claim$test == "regression_slope_by_sex") {
      result <- CompareSlopesBySex(values, measures$severity_index, measures$sex)
      verdict <- if (is.na(result$observed)) {
        "not measured"
      } else if (result$observed == "females fall faster") {
        "reproduced"
      } else {
        "opposite"
      }
      return(Row(
        claim,
        sprintf("female slope %.3f against male slope %.3f",
                result$female_slope, result$male_slope),
        result$observed, result$slope_difference, NA_real_, verdict
      ))
    }

    if (claim$test == "sex_difference_by_timepoint") {
      # Claims 5 and 6 are about whether a difference is detectable in each
      # window, not about which sex is numerically higher. The paper states the
      # criterion in the legend of Figure 2: significance was determined by a
      # Mann-Whitney test, and an asterisk means p < 0.05. A difference that is
      # "lost" is one that stops being significant, and it keeps a direction
      # while it is lost. Testing the direction instead answers a question the
      # paper never asked.
      per_timepoint <- lapply(c("early", "middle", "late"), function(tp) {
        keep <- measures$timepoint == tp
        result <- CompareBySex(values[keep], measures$sex[keep])
        result$timepoint <- tp
        result
      })
      combined <- do.call(rbind, per_timepoint)

      if (any(is.na(combined$p_value))) {
        return(Row(claim, "a window holds too few samples to test",
                   NA_character_, NA_real_, NA_real_, "not measured"))
      }

      combined$significant <- combined$p_value < alpha
      # Claim 5 expects no significant difference at early and middle and a
      # significant one at late. Claim 6 expects one in every window.
      expected_significant <- if (claim$claim_id == 5) {
        c(early = FALSE, middle = FALSE, late = TRUE)
      } else {
        c(early = TRUE, middle = TRUE, late = TRUE)
      }
      agree <- combined$significant == expected_significant[combined$timepoint]

      verdict <- if (all(agree)) {
        "reproduced"
      } else if (!any(agree)) {
        "opposite"
      } else {
        "partly reproduced"
      }

      detail <- paste(
        sprintf("%s p %.3f %s, %s higher", combined$timepoint,
                combined$p_value, ifelse(combined$significant, "*", "ns"),
                sub(" higher", "", combined$observed)),
        collapse = "; "
      )

      return(Row(
        claim, detail,
        paste(sprintf("%s %s", combined$timepoint,
                      ifelse(combined$significant, "significant",
                             "not significant")), collapse = "; "),
        sum(agree), max(combined$p_value), verdict
      ))
    }

    if (claim$test %in% c("sex_difference_before_seroconversion",
                          "sex_difference_after_seroconversion")) {
      wanted_result <- if (grepl("before", claim$test)) "NEGATIVE" else "POSITIVE"
      keep <- !is.na(infected$igg_result) & infected$igg_result == wanted_result
      subset_values <- infected[[claim$measure]][keep]
      result <- CompareBySex(subset_values, infected$sex[keep])
      wanted <- if (grepl("^female", claim$expected)) {
        "female higher"
      } else {
        "male higher"
      }
      verdict <- if (is.na(result$observed)) {
        "not measured"
      } else if (result$observed == wanted) {
        "reproduced"
      } else {
        "opposite"
      }
      return(Row(
        claim,
        sprintf("IgG %s, female median %.2f against male median %.2f, n = %d and %d",
                tolower(wanted_result), result$female_median, result$male_median,
                result$n_female, result$n_male),
        result$observed, result$difference, result$p_value, verdict
      ))
    }

    if (claim$test == "no_change_across_severity") {
      result <- CorrelateWithSeverity(values, measures$severity_index)
      verdict <- if (is.na(result$p_value)) {
        "not measured"
      } else if (result$p_value >= alpha) {
        "reproduced"
      } else {
        "opposite"
      }
      return(Row(claim, sprintf("Spearman against severity rank, n = %d", result$n),
                 sprintf("rho = %.3f", result$rho), result$rho, result$p_value,
                 verdict))
    }

    if (claim$test == "positive_control") {
      nk <- measures$nk_cd56_median
      cd8 <- measures$cd8_cd56_median
      if (is.null(nk) || is.null(cd8) || all(is.na(nk)) || all(is.na(cd8))) {
        return(Row(claim, "no CD56 median was collected", NA_character_,
                   NA_real_, NA_real_, "not measured"))
      }
      difference <- stats::median(nk, na.rm = TRUE) - stats::median(cd8, na.rm = TRUE)
      observed <- if (difference > 0) "NK higher" else "CD8 T higher"
      return(Row(
        claim,
        sprintf("NK median %.3f against CD8 T median %.3f",
                stats::median(nk, na.rm = TRUE), stats::median(cd8, na.rm = TRUE)),
        observed, difference, NA_real_,
        if (difference > 0) "reproduced" else "opposite"
      ))
    }

    Row(claim, paste("unknown test:", claim$test), NA_character_, NA_real_,
        NA_real_, "not measured")
  })

  result <- do.call(rbind, rows)
  rownames(result) <- NULL
  result
}
