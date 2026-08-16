# Tests for R/counting.R
#
# The paper states a 5 percent CV target and a 400 to 2000 event range. Those two
# numbers are consistent with each other only if the arithmetic below is right, so
# the arithmetic is tested against the paper's own pair.

test_that("PoissonCv reproduces the paper's own pairing of 400 events and 5 percent", {
  expect_equal(PoissonCv(400), 5)
})

test_that("PoissonCv falls as the count rises", {
  expect_equal(PoissonCv(2500), 2)
  expect_equal(PoissonCv(10000), 1)
  expect_lt(PoissonCv(2000), PoissonCv(400))
})

test_that("PoissonCv is undefined for a count of zero", {
  expect_true(is.na(PoissonCv(0)))
})

test_that("PoissonCv works on a vector", {
  result <- PoissonCv(c(400, 2500, 0))

  expect_equal(length(result), 3)
  expect_equal(result[1], 5)
  expect_equal(result[2], 2)
  expect_true(is.na(result[3]))
})

test_that("PoissonCv rejects a negative count", {
  expect_error(PoissonCv(-1), "cannot be negative")
})

test_that("EventsForCv inverts PoissonCv", {
  expect_equal(EventsForCv(5), 400)
  expect_equal(EventsForCv(2), 2500)
  expect_equal(PoissonCv(EventsForCv(5)), 5)
})

test_that("EventsForCv rejects a target of zero or less", {
  expect_error(EventsForCv(0), "above zero")
  expect_error(EventsForCv(-3), "above zero")
})

test_that("MeasuredCv computes the spread of a set of values", {
  values <- c(3.266, 3.217, 3.491, 3.391, 3.201, 3.434, 3.505)
  expected <- 100 * stats::sd(values) / mean(values)

  expect_equal(MeasuredCv(values), expected)
})

test_that("MeasuredCv returns NA when it cannot be computed", {
  expect_true(is.na(MeasuredCv(5)))
  expect_true(is.na(MeasuredCv(numeric())))
  expect_true(is.na(MeasuredCv(c(0, 0, 0))))
})

test_that("MeasuredCv ignores missing values", {
  expect_equal(MeasuredCv(c(10, 12, NA, 14)), MeasuredCv(c(10, 12, 14)))
})

test_that("CompareSpreadToPoisson separates counting noise from the rest", {
  # A measured CV of 5 with a Poisson floor of 3 leaves 4 in quadrature,
  # because 3 squared plus 4 squared is 5 squared.
  counts <- rep((100 / 3)^2, 4)
  frequencies <- c(95, 100, 105, 100)
  result <- CompareSpreadToPoisson(counts, frequencies)

  expect_equal(result$poisson_cv_percent, 3)
  expect_equal(
    result$excess_cv_percent,
    sqrt(result$measured_cv_percent^2 - 9)
  )
  expect_equal(result$replicates, 4)
})

test_that("CompareSpreadToPoisson returns NA when the spread is below the floor", {
  # Identical frequencies give a measured CV of zero, which is below any floor.
  result <- CompareSpreadToPoisson(rep(400, 4), rep(2.5, 4))

  expect_equal(result$measured_cv_percent, 0)
  expect_true(is.na(result$excess_cv_percent))
})

test_that("CompareSpreadToPoisson rejects mismatched lengths", {
  expect_error(
    CompareSpreadToPoisson(c(1, 2, 3), c(1, 2)),
    "same length"
  )
})
