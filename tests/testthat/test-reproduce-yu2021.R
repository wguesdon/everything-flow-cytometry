# Tests for R/reproduce_yu2021.R
#
# The claim readers and the two statistical helpers are pure, so they are tested
# on small vectors. The functions that gate an FCS file are exercised by
# scripts/08_yu2021_spectral_mait.R, because they need a real spectral file.

test_that("ReadYuClaims reads a claims file that carries every column", {
  path <- withr::local_tempfile(fileext = ".csv")
  writeLines(c("claim_id,short_name,measure,test,expected,quote,figure",
               "1,MAIT fall,mait_percent,lower,a fall,a sentence,Fig 1"), path)
  claims <- ReadYuClaims(path)
  expect_equal(nrow(claims), 1)
  expect_equal(claims$measure, "mait_percent")
})

test_that("ReadYuClaims rejects a claims file with no figure column", {
  path <- withr::local_tempfile(fileext = ".csv")
  writeLines(c("claim_id,short_name,measure,test,expected,quote",
               "1,MAIT fall,mait_percent,lower,a fall,a sentence"), path)
  expect_error(ReadYuClaims(path), "missing the column\\(s\\): figure")
})

test_that("ReadYuClaims rejects a path that does not exist", {
  expect_error(ReadYuClaims(file.path(tempdir(), "absent.csv")),
               "does not exist")
})

test_that("CompareBySex reports a median for each sex and a p value", {
  values <- c(10, 12, 11, 30, 32, 31)
  sex <- c(rep("Female", 3), rep("Male", 3))
  result <- CompareBySex(values, sex)
  expect_equal(result$n_female, 3L)
  expect_equal(result$n_male, 3L)
  expect_equal(result$female_median, 11)
  expect_equal(result$male_median, 31)
  expect_true(is.finite(result$p_value))
})

test_that("CompareBySex drops a sample with no value or no sex", {
  values <- c(10, 12, NA, 30, 32)
  sex <- c("Female", "Female", "Female", "Male", NA)
  result <- CompareBySex(values, sex)
  expect_equal(result$n_female, 2L)
  expect_equal(result$n_male, 1L)
})

test_that("CompareBySex runs no test when one sex is absent", {
  result <- CompareBySex(c(10, 12), c("Female", "Female"))
  expect_equal(result$n_male, 0L)
  expect_true(is.na(result$p_value))
})

test_that("CompareBySex runs no test when a sex holds one sample", {
  result <- CompareBySex(c(10, 12, 30), c("Female", "Female", "Male"))
  expect_equal(result$n_male, 1L)
  expect_true(is.na(result$p_value))
})

test_that("CorrelateWithSeverity reports the direction of the correlation", {
  rising <- CorrelateWithSeverity(c(1, 2, 3, 4, 5), c(1, 2, 3, 4, 5))
  expect_equal(rising$observed, "rises with severity")
  expect_gt(rising$rho, 0.9)

  falling <- CorrelateWithSeverity(c(5, 4, 3, 2, 1), c(1, 2, 3, 4, 5))
  expect_equal(falling$observed, "falls with severity")
  expect_lt(falling$rho, -0.9)
})

test_that("CorrelateWithSeverity runs no test on fewer than three samples", {
  result <- CorrelateWithSeverity(c(1, 2), c(1, 2))
  expect_equal(result$n, 2L)
  expect_true(is.na(result$rho))
})

test_that("CorrelateWithSeverity drops a sample with no value", {
  result <- CorrelateWithSeverity(c(1, NA, 3, 4), c(1, 2, 3, 4))
  expect_equal(result$n, 3L)
})

test_that("CompareSlopesBySex fits one slope for each sex", {
  values <- c(1, 2, 3, 4, 10, 8, 6, 4)
  severity <- rep(1:4, times = 2)
  sex <- rep(c("Female", "Male"), each = 4)
  result <- CompareSlopesBySex(values, severity, sex)
  expect_gt(result$female_slope, 0)
  expect_lt(result$male_slope, 0)
})

test_that("CompareSlopesBySex gives NA when a sex holds too few samples", {
  values <- c(1, 2, 3, 4, 10)
  severity <- c(1, 2, 3, 4, 1)
  sex <- c(rep("Female", 4), "Male")
  result <- CompareSlopesBySex(values, severity, sex)
  expect_true(is.na(result$male_slope))
  expect_true(is.na(result$observed))
})

test_that("CompareSlopesBySex gives NA when severity never varies", {
  values <- c(1, 2, 3, 4, 5, 6)
  severity <- rep(2, 6)
  sex <- rep(c("Female", "Male"), each = 3)
  result <- CompareSlopesBySex(values, severity, sex)
  expect_true(is.na(result$female_slope))
})
