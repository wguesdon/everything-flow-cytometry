
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
