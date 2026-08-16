# Test a paper's stated findings against numbers computed from its own data.
#
# OMIP-039 makes seven directional claims about how CD56dim NKG2C+ NK cells
# differ from CD56dim NKG2C- NK cells, and it states a prototypic adaptive
# phenotype as a boolean combination of six markers. The claims live in
# gating/omip39_paper_claims.csv with the sentence each one came from, so a
# reader can check the claim against the paper before checking the number
# against the data.
#
# The FlowJo workspace the authors deposited holds 64 boolean populations under
# each NKG2C branch. Every marker frequency in Figure 1B can therefore be derived
# by summing those populations, which means the comparison uses the authors' own
# gates and adds no gating decision of its own.

#' Read the claims a paper makes, with the sentence each one came from
#'
#' @param path Path to a CSV with the columns `claim_id`, `marker`, `measure`,
#'   `direction` and `quote`. `direction` is `"higher"` or `"lower"`, and it
#'   describes the test population relative to the reference population.
#' @return A `data.frame` of claims.
#' @export
ReadPaperClaims <- function(path) {
  if (!file.exists(path)) {
    stop("The claims file does not exist: ", path)
  }

  claims <- utils::read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
  required <- c("claim_id", "marker", "measure", "direction", "quote")
  missing <- setdiff(required, colnames(claims))
  if (length(missing) > 0) {
    stop("The claims file is missing the column(s): ",
         paste(missing, collapse = ", "), ".")
  }

  bad_direction <- setdiff(claims$direction, c("higher", "lower"))
  if (length(bad_direction) > 0) {
    stop("A claim direction must be 'higher' or 'lower'. Found: ",
         paste(unique(bad_direction), collapse = ", "), ".")
  }

  claims
}

#' Read the marker signs out of a boolean population name
#'
#' FlowJo writes a boolean population as a run of marker names each followed by
#' `+` or `-`, for example `CD2+CD57-ILT2+NKG2A-NKp30-Siglec-7-`. The marker
#' names themselves can hold a hyphen, as `Siglec-7` does, so the parse anchors on
#' the known marker names rather than splitting on punctuation.
#'
#' @param population_name A single boolean population name.
#' @param markers The marker names to look for.
#' @return A named logical vector, `TRUE` where the marker carries `+`. A marker
#'   that does not appear in the name is `NA`.
#' @examples
#' ParseBooleanPopulation(
#'   "CD2+CD57-ILT2+NKG2A-NKp30-Siglec-7-",
#'   c("CD2", "CD57", "ILT2", "NKG2A", "NKp30", "Siglec-7")
#' )
#' @export
ParseBooleanPopulation <- function(population_name, markers) {
  leaf <- basename(as.character(population_name))

  # Fixed string matching, not a regular expression. A marker name can hold a
  # hyphen, as `Siglec-7` does, and escaping that reliably for every marker is
  # more fragile than searching for the two literal strings.
  signs <- vapply(markers, function(marker) {
    has_positive <- grepl(paste0(marker, "+"), leaf, fixed = TRUE)
    has_negative <- grepl(paste0(marker, "-"), leaf, fixed = TRUE)

    if (has_positive && !has_negative) {
      return(TRUE)
    }
    if (has_negative && !has_positive) {
      return(FALSE)
    }
    # Neither sign, or both, which means the name does not resolve this marker.
    NA
  }, logical(1))

  stats::setNames(signs, markers)
}

#' Derive the frequency of every single marker from the boolean populations
#'
#' A boolean population set partitions the parent, so the frequency of one marker
#' is the sum of the counts of every population where that marker is positive,
#' divided by the sum over all of them. Deriving the frequency this way uses only
#' the gates the analyst already drew.
#'
#' @param stats The output of [CollectPopulationStats()], filtered to the boolean
#'   populations of one parent.
#' @param markers The marker names to derive.
#' @return A `data.frame` with the columns `marker`, `positive_events`,
#'   `total_events` and `percent_positive`.
#' @export
MarkerFrequencyFromBooleans <- function(stats, markers) {
  if (nrow(stats) == 0) {
    stop("No boolean population was given, so no frequency can be derived.")
  }
  if (!"count" %in% colnames(stats)) {
    stop("stats has no 'count' column.")
  }

  signs <- lapply(stats$population, ParseBooleanPopulation, markers = markers)
  sign_matrix <- do.call(rbind, signs)

  total <- sum(stats$count)

  rows <- lapply(markers, function(marker) {
    is_positive <- sign_matrix[, marker]
    usable <- !is.na(is_positive)
    positive_events <- sum(stats$count[usable & is_positive])

    data.frame(
      marker = marker,
      positive_events = positive_events,
      total_events = total,
      percent_positive = if (total == 0) NA_real_ else 100 * positive_events / total,
      stringsAsFactors = FALSE
    )
  })

  do.call(rbind, rows)
}

#' Test each claim against the two populations it compares
#'
#' @param claims The output of [ReadPaperClaims()].
#' @param test_frequencies Marker frequencies in the population the claim is about,
#'   from [MarkerFrequencyFromBooleans()]. For OMIP-039 this is CD56dim NKG2C+.
#' @param reference_frequencies Marker frequencies in the population it is compared
#'   against. For OMIP-039 this is CD56dim NKG2C-.
#' @param min_difference_points A difference smaller than this, in percentage
#'   points, is reported as `"too small to call"` rather than as a direction. The
#'   default of 1 point is a judgement and a report must say so.
#' @return A `data.frame` with one row per claim, holding the claim, both
#'   measured values, the difference, and a `verdict` of `"reproduced"`,
#'   `"opposite"`, `"too small to call"` or `"not measured"`.
#' @export
TestPaperClaims <- function(claims,
                            test_frequencies,
                            reference_frequencies,
                            min_difference_points = 1) {
  test_lookup <- stats::setNames(
    test_frequencies$percent_positive, test_frequencies$marker
  )
  reference_lookup <- stats::setNames(
    reference_frequencies$percent_positive, reference_frequencies$marker
  )

  result <- claims
  result$test_percent <- unname(test_lookup[claims$marker])
  result$reference_percent <- unname(reference_lookup[claims$marker])
  result$difference_points <- result$test_percent - result$reference_percent

  observed <- ifelse(
    is.na(result$difference_points), NA_character_,
    ifelse(result$difference_points > 0, "higher", "lower")
  )
  result$observed_direction <- observed

  result$verdict <- ifelse(
    is.na(result$difference_points), "not measured",
    ifelse(
      abs(result$difference_points) < min_difference_points, "too small to call",
      ifelse(observed == result$direction, "reproduced", "opposite")
    )
  )

  result
}

#' Find the prototypic adaptive phenotype among the boolean populations
#'
#' OMIP-039 states one combination as prototypic for adaptive NK cells:
#' NKG2A negative, NKp30 negative, Siglec-7 negative, ILT2 positive, CD57
#' positive and CD2 positive. This locates that population in each branch and
#' reports its frequency, which is the claim of Figure 1C.
#'
#' @param stats The output of [CollectPopulationStats()], filtered to the boolean
#'   populations of one parent.
#' @param phenotype A named logical vector giving the required sign of every
#'   marker.
#' @return A one row `data.frame` with `population`, `count`, `total_events` and
#'   `percent_of_parent`. `population` is `NA` when no population matches.
#' @export
FindPhenotype <- function(stats, phenotype) {
  markers <- names(phenotype)
  signs <- lapply(stats$population, ParseBooleanPopulation, markers = markers)
  sign_matrix <- do.call(rbind, signs)

  matches <- apply(sign_matrix, 1, function(row) {
    !any(is.na(row)) && all(row == phenotype)
  })

  total <- sum(stats$count)

  if (!any(matches)) {
    return(data.frame(
      population = NA_character_, count = NA_integer_,
      total_events = total, percent_of_parent = NA_real_,
      stringsAsFactors = FALSE
    ))
  }

  hit <- which(matches)[1]
  data.frame(
    population = basename(as.character(stats$population[hit])),
    count = stats$count[hit],
    total_events = total,
    percent_of_parent = if (total == 0) NA_real_ else 100 * stats$count[hit] / total,
    stringsAsFactors = FALSE
  )
}
