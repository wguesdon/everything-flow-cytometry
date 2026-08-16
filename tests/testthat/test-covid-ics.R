# Tests for R/covid_ics.R

test_that("ParseCovidFileNames reads the clinical group", {
  result <- ParseCovidFileNames(c(
    "export_COVID19 samples 21_04_20_ST3_COVID19_ICU_005_A ST3 210420_080_Live_cells.fcs",
    "export_COVID19 samples 23_04_20_ST3_COVID19_HC_001 ST3 230420_017_Live_cells.fcs",
    "export_COVID19 samples 23_04_20_ST3_COVID19_W_022_O ST3 230420_006_Live_cells.fcs"
  ))
  expect_equal(result$group, c("intensive care", "healthy", "ward"))
  expect_equal(result$patient, c("ICU005", "HC001", "W022"))
})

test_that("ParseCovidFileNames reads a relabelled patient as intensive care", {
  result <- ParseCovidFileNames(
    "export_COVID19 samples 21_04_20_ST3_COVID19_ICU_changedW_002_O ST3 210420_043_Live_cells.fcs"
  )
  expect_equal(result$group, "intensive care")
})

test_that("ParseCovidFileNames returns NA for a name with no group", {
  result <- ParseCovidFileNames("something_else.fcs")
  expect_true(is.na(result$group))
  expect_true(is.na(result$patient))
})

test_that("CompareByGroup returns a median per group and two p values", {
  values <- c(rep(1, 6), rep(5, 20), rep(6, 23))
  group <- c(rep("healthy", 6), rep("ward", 20), rep("intensive care", 23))
  result <- CompareByGroup(values, group)
  expect_equal(result$healthy_median, 1)
  expect_equal(result$ward_median, 5)
  expect_equal(result$intensive_care_median, 6)
  expect_true(result$covid_against_healthy_p < 0.001)
  expect_true(result$between_severity_p < 0.001)
})

test_that("CompareByGroup refuses a single group", {
  expect_null(CompareByGroup(c(1, 2, 3), rep("healthy", 3)))
})

test_that("CompareByGroup ignores values it cannot use", {
  values <- c(1, NA, 3, Inf, 5, 7, 9, 11)
  group <- c(rep("healthy", 4), rep("ward", 4))
  result <- CompareByGroup(values, group)
  expect_equal(result$healthy_median, 2)
  expect_true(is.na(result$intensive_care_median))
})

test_that("ReadCovidManualGates refuses a workspace that is not there", {
  expect_error(
    ReadCovidManualGates(file.path(tempdir(), "absent.wsp"), tempdir()),
    "does not exist"
  )
})
