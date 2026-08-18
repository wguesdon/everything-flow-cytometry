
MakeClaims. <- function() {
  data.frame(
    claim_id = 1:6,
    claim = paste("claim", 1:6),
    measure = c("a", "a", "a", "a", "a", "absent"),
    test = c("at_least", "at_most", "between", "equals", "greater_than",
             "at_least"),
    expected = c("10", "10", "5,15", "20", "b", "1"),
    stringsAsFactors = FALSE)
}

test_that("JudgeClaims puts one verdict on every claim", {
  results <- data.frame(measure = c("a", "b"), value = c(12, 30))
  judged <- JudgeClaims(MakeClaims.(), results)
  expect_equal(nrow(judged), 6)
  expect_true(all(judged$verdict %in% kClaimVerdicts))
})

test_that("JudgeClaims reads at_least and at_most against the same value", {
  results <- data.frame(measure = "a", value = 12)
  judged <- JudgeClaims(MakeClaims.()[1:2, ], results)
  expect_equal(judged$verdict, c("supported", "contradicted"))
})

test_that("JudgeClaims reads a between claim with two numbers", {
  results <- data.frame(measure = "a", value = 12)
  expect_equal(JudgeClaims(MakeClaims.()[3, ], results)$verdict, "supported")
  results$value <- 20
  expect_equal(JudgeClaims(MakeClaims.()[3, ], results)$verdict,
               "contradicted")
})

test_that("JudgeClaims reads an equals claim inside the tolerance", {
  results <- data.frame(measure = "a", value = 20.5)
  expect_equal(JudgeClaims(MakeClaims.()[4, ], results,
                           tolerance = 0.05)$verdict, "supported")
  expect_equal(JudgeClaims(MakeClaims.()[4, ], results,
                           tolerance = 0.001)$verdict, "contradicted")
})

test_that("JudgeClaims compares one measure against another", {
  results <- data.frame(measure = c("a", "b"), value = c(40, 30))
  expect_equal(JudgeClaims(MakeClaims.()[5, ], results)$verdict, "supported")
  results$value <- c(20, 30)
  expect_equal(JudgeClaims(MakeClaims.()[5, ], results)$verdict,
               "contradicted")
})

test_that("JudgeClaims calls a claim unresolved when nothing measured it", {
  # Unresolved is a real answer. Calling it contradicted would say the data
  # disagreed, when the data was never asked.
  results <- data.frame(measure = "a", value = 12)
  judged <- JudgeClaims(MakeClaims.()[6, ], results)
  expect_equal(judged$verdict, "unresolved")
  expect_match(judged$reason, "no result names the measure 'absent'")
  expect_true(is.na(judged$observed))
})

test_that("JudgeClaims calls a claim unresolved when the other measure is absent", {
  results <- data.frame(measure = "a", value = 40)
  judged <- JudgeClaims(MakeClaims.()[5, ], results)
  expect_equal(judged$verdict, "unresolved")
  expect_match(judged$reason, "to compare against")
})

test_that("JudgeClaims calls a claim unresolved when the test is unknown", {
  claims <- MakeClaims.()[1, ]
  claims$test <- "vibes"
  judged <- JudgeClaims(claims, data.frame(measure = "a", value = 12))
  expect_equal(judged$verdict, "unresolved")
  expect_match(judged$reason, "is not a test this reads")
})

test_that("JudgeClaims calls a claim unresolved when the expected is not a number", {
  claims <- MakeClaims.()[1, ]
  claims$expected <- "quite high"
  judged <- JudgeClaims(claims, data.frame(measure = "a", value = 12))
  expect_equal(judged$verdict, "unresolved")
  expect_match(judged$reason, "is not a number")
})

test_that("JudgeClaims supports a present claim whenever the measure exists", {
  claims <- MakeClaims.()[1, ]
  claims$test <- "present"
  expect_equal(JudgeClaims(claims, data.frame(measure = "a",
                                              value = 12))$verdict,
               "supported")
})

test_that("JudgeClaims rejects a claims table with no test column", {
  claims <- MakeClaims.()
  claims$test <- NULL
  expect_error(JudgeClaims(claims, data.frame(measure = "a", value = 1)),
               "missing the column\\(s\\): test")
})

test_that("JudgeClaims rejects a results table with no value column", {
  expect_error(JudgeClaims(MakeClaims.(), data.frame(measure = "a")),
               "needs the columns measure and value")
})

test_that("ReadPaperClaims accepts higher and lower as a direction", {
  path <- withr::local_tempfile(fileext = ".csv")
  writeLines(c("claim_id,marker,measure,direction,quote",
               "1,CD57,percent_positive,higher,a sentence",
               "2,NKG2A,percent_positive,lower,another sentence"), path)
  claims <- ReadPaperClaims(path)
  expect_equal(nrow(claims), 2)
  expect_equal(claims$direction, c("higher", "lower"))
})

test_that("ReadPaperClaims rejects a direction it cannot read", {
  # "increased" reads to a person and not to the comparison, so it would pass
  # through and match nothing.
  path <- withr::local_tempfile(fileext = ".csv")
  writeLines(c("claim_id,marker,measure,direction,quote",
               "1,CD57,percent_positive,increased,a sentence"), path)
  expect_error(ReadPaperClaims(path), "must be 'higher' or 'lower'")
})

test_that("ReadPaperClaims rejects a file with no quote column", {
  path <- withr::local_tempfile(fileext = ".csv")
  writeLines(c("claim_id,marker,measure,direction", "1,CD57,x,higher"), path)
  expect_error(ReadPaperClaims(path), "missing the column\\(s\\): quote")
})

test_that("ParseBooleanPopulation reads the sign of every marker", {
  signs <- ParseBooleanPopulation("CD57+NKG2A-CD16+", c("CD57", "NKG2A",
                                                       "CD16"))
  expect_equal(unname(signs), c(TRUE, FALSE, TRUE))
  expect_equal(names(signs), c("CD57", "NKG2A", "CD16"))
})

test_that("ParseBooleanPopulation reads a marker whose name holds a hyphen", {
  # Siglec-7 is why the search is on a fixed string and not a pattern.
  signs <- ParseBooleanPopulation("Siglec-7+CD57-", c("Siglec-7", "CD57"))
  expect_equal(unname(signs), c(TRUE, FALSE))
})

test_that("ParseBooleanPopulation gives NA when the name resolves nothing", {
  signs <- ParseBooleanPopulation("CD57+", c("CD57", "NKG2A"))
  expect_true(is.na(signs[["NKG2A"]]))
})

test_that("ParseBooleanPopulation reads the leaf of a full path", {
  signs <- ParseBooleanPopulation("/Live/NK/CD57+NKG2A-", c("CD57", "NKG2A"))
  expect_equal(unname(signs), c(TRUE, FALSE))
})

MakeBooleanStats. <- function() {
  data.frame(
    population = c("CD57+NKG2A+", "CD57+NKG2A-", "CD57-NKG2A+",
                   "CD57-NKG2A-"),
    count = c(10, 20, 30, 40), stringsAsFactors = FALSE)
}

test_that("MarkerFrequencyFromBooleans adds up the positive populations", {
  frequencies <- MarkerFrequencyFromBooleans(MakeBooleanStats.(),
                                             c("CD57", "NKG2A"))
  expect_equal(frequencies$marker, c("CD57", "NKG2A"))
  expect_equal(frequencies$positive_events, c(30, 40))
  expect_equal(frequencies$total_events, c(100, 100))
  expect_equal(frequencies$percent_positive, c(30, 40))
})

test_that("MarkerFrequencyFromBooleans rejects an empty table", {
  expect_error(MarkerFrequencyFromBooleans(MakeBooleanStats.()[0, ], "CD57"),
               "No boolean population was given")
})

test_that("MarkerFrequencyFromBooleans rejects a table with no count", {
  stats <- MakeBooleanStats.()
  stats$count <- NULL
  expect_error(MarkerFrequencyFromBooleans(stats, "CD57"),
               "no 'count' column")
})

MakeDirectionalClaims. <- function() {
  data.frame(claim_id = 1:3, marker = c("CD57", "NKG2A", "NKp30"),
             measure = "percent_positive",
             direction = c("higher", "lower", "higher"),
             quote = "a sentence", stringsAsFactors = FALSE)
}

test_that("TestPaperClaims calls a claim reproduced when the direction matches",
          {
  test <- data.frame(marker = c("CD57", "NKG2A"),
                     percent_positive = c(80, 10), stringsAsFactors = FALSE)
  reference <- data.frame(marker = c("CD57", "NKG2A"),
                          percent_positive = c(40, 50),
                          stringsAsFactors = FALSE)
  judged <- TestPaperClaims(MakeDirectionalClaims.(), test, reference)
  expect_equal(judged$verdict[1:2], c("reproduced", "reproduced"))
})

test_that("TestPaperClaims calls a claim opposite when the direction inverts", {
  test <- data.frame(marker = "CD57", percent_positive = 10,
                     stringsAsFactors = FALSE)
  reference <- data.frame(marker = "CD57", percent_positive = 40,
                          stringsAsFactors = FALSE)
  judged <- TestPaperClaims(MakeDirectionalClaims.()[1, ], test, reference)
  expect_equal(judged$verdict, "opposite")
})

test_that("TestPaperClaims will not call a difference below the floor", {
  # Half a point in either direction is measurement and not a finding.
  test <- data.frame(marker = "CD57", percent_positive = 40.4,
                     stringsAsFactors = FALSE)
  reference <- data.frame(marker = "CD57", percent_positive = 40,
                          stringsAsFactors = FALSE)
  judged <- TestPaperClaims(MakeDirectionalClaims.()[1, ], test, reference,
                            min_difference_points = 1)
  expect_equal(judged$verdict, "too small to call")
})

test_that("TestPaperClaims calls a marker nobody measured not measured", {
  test <- data.frame(marker = "CD57", percent_positive = 80,
                     stringsAsFactors = FALSE)
  reference <- data.frame(marker = "CD57", percent_positive = 40,
                          stringsAsFactors = FALSE)
  judged <- TestPaperClaims(MakeDirectionalClaims.(), test, reference)
  expect_equal(judged$verdict[3], "not measured")
  expect_true(is.na(judged$difference_points[3]))
})

test_that("FindPhenotype finds the population that matches every sign", {
  found <- FindPhenotype(MakeBooleanStats.(),
                         c(CD57 = TRUE, NKG2A = FALSE))
  expect_equal(found$population, "CD57+NKG2A-")
  expect_equal(found$count, 20)
  expect_equal(found$percent_of_parent, 20)
})

test_that("FindPhenotype returns NA when no population matches", {
  found <- FindPhenotype(MakeBooleanStats.()[1, ],
                         c(CD57 = TRUE, NKG2A = FALSE))
  expect_true(is.na(found$population))
  expect_true(is.na(found$percent_of_parent))
  expect_equal(found$total_events, 10)
})
