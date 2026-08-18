# Tests for R/gating.R
#
# The gating functions that need openCyto are exercised by the pipeline script,
# because a real GatingSet needs a real flowSet. The tests here cover the parts
# that can be checked without one: the template contract and the arithmetic that
# turns counts into a spread.

test_that("ReadGatingTemplate rejects a path that does not exist", {
  expect_error(
    ReadGatingTemplate("/no/such/template.csv"),
    "The gating template does not exist"
  )
})

test_that("ReadGatingTemplate names the columns that are missing", {
  path <- withr::local_tempfile(fileext = ".csv")
  utils::write.csv(
    data.frame(alias = "CD3", pop = "+", parent = "root"),
    path, row.names = FALSE
  )

  expect_error(ReadGatingTemplate(path), "missing the column")
  expect_error(ReadGatingTemplate(path), "dims")
  expect_error(ReadGatingTemplate(path), "gating_method")
})

test_that("the repository template holds every required column", {
  path <- testthat::test_path("..", "..", "gating", "pbmc_gating_template.csv")
  skip_if_not(file.exists(path),
              "the gating template is not at the expected path")

  template <- utils::read.csv(path, check.names = FALSE)
  required <- c("alias", "pop", "parent", "dims", "gating_method")

  expect_true(all(required %in% colnames(template)))
  expect_gt(nrow(template), 0)
})

test_that("the repository template names a parent that is defined before it", {
  path <- testthat::test_path("..", "..", "gating", "pbmc_gating_template.csv")
  skip_if_not(file.exists(path),
              "the gating template is not at the expected path")

  template <- utils::read.csv(path, check.names = FALSE)
  defined <- "root"

  for (i in seq_len(nrow(template))) {
    expect_true(
      template$parent[i] %in% defined,
      info = paste0(
        "Row ", i, " gates '", template$alias[i], "' from parent '",
        template$parent[i], "', which is not defined above it."
      )
    )
    defined <- c(defined, template$alias[i])
  }
})

test_that("RunAutomatedGating rejects an object that is not a flowSet", {
  expect_error(RunAutomatedGating(data.frame(a = 1), NULL), "must be a flowSet")
})

test_that("SummarisePopulationSpread computes the CV of a frequency", {
  stats <- data.frame(
    sample = paste0("s", 1:4),
    population = "CD3",
    percent_of_parent = c(40, 50, 60, 50),
    stringsAsFactors = FALSE
  )
  result <- SummarisePopulationSpread(stats)

  expect_equal(nrow(result), 1)
  expect_equal(result$n, 4)
  expect_equal(result$mean_percent, 50)
  # sd of 40, 50, 60, 50 is 8.164966, so the CV is 16.33 percent.
  expect_equal(result$sd_percent, stats::sd(c(40, 50, 60, 50)))
  expect_equal(result$cv_percent, 100 * stats::sd(c(40, 50, 60, 50)) / 50)
})

test_that("SummarisePopulationSpread splits by a grouping column", {
  stats <- data.frame(
    sample = paste0("s", 1:4),
    population = "CD3",
    condition = c("NS", "NS", "PI", "PI"),
    percent_of_parent = c(40, 60, 45, 55),
    stringsAsFactors = FALSE
  )
  result <- SummarisePopulationSpread(stats, group_by = "condition")

  expect_equal(nrow(result), 2)
  expect_equal(sort(result$condition), c("NS", "PI"))
  expect_equal(result$n, c(2, 2))
  expect_true(all(result$mean_percent == 50))
})

test_that("SummarisePopulationSpread returns NA for a zero mean", {
  stats <- data.frame(
    sample = c("s1", "s2"),
    population = "CD3",
    percent_of_parent = c(0, 0),
    stringsAsFactors = FALSE
  )
  result <- SummarisePopulationSpread(stats)

  expect_true(is.na(result$cv_percent))
})

test_that("SummarisePopulationSpread names a column it cannot find", {
  stats <- data.frame(
    sample = "s1", population = "CD3", percent_of_parent = 50,
    stringsAsFactors = FALSE
  )

  expect_error(
    SummarisePopulationSpread(stats, group_by = "treatment"),
    "no column named treatment"
  )
})

test_that("SummarisePopulationSpread rejects a table with no frequency", {
  expect_error(
    SummarisePopulationSpread(data.frame(population = "CD3")),
    "missing the column"
  )
})

test_that("CollectGateTree puts the root first and gives it no parent", {
  # PlotGateTree reads a missing parent as the root of the drawing, so the
  # first row has to carry NA and not the string root.
  gating_set <- MakeGatedSet()
  tree <- CollectGateTree(gating_set)
  expect_equal(tree$population[1], "all_events")
  expect_true(is.na(tree$parent[1]))
  expect_equal(tree$percent_of_parent[1], 100)
})

test_that("CollectGateTree names the parent of every population below root", {
  gating_set <- MakeGatedSet()
  tree <- CollectGateTree(gating_set)
  below <- tree[!is.na(tree$parent), , drop = FALSE]
  expect_true(all(below$parent %in% tree$population))
})

test_that("CollectGateTree reports the percentage of the parent", {
  gating_set <- MakeGatedSet()
  tree <- CollectGateTree(gating_set)
  below <- tree[!is.na(tree$parent), , drop = FALSE]
  expect_equal(below$percent_of_parent,
               100 * below$events / below$parent_events)
  expect_true(all(below$percent_of_parent <= 100))
})

test_that("CollectGateTree gives PlotGateTree a table it can draw", {
  gating_set <- MakeGatedSet()
  expect_s3_class(PlotGateTree(CollectGateTree(gating_set)), "ggplot")
})

test_that("CollectGateTree rejects a set with no population below root", {
  flow_set <- flowCore::flowSet(sample = MakeTestFlowFrame())
  empty <- flowWorkspace::GatingSet(flow_set)
  expect_error(CollectGateTree(empty), "holds no population below root")
})

test_that("CollectPopulationStats reports one row per sample and population", {
  gating_set <- MakeGatedSet()
  stats <- CollectPopulationStats(gating_set)
  expect_true(all(c("sample", "population", "count", "percent_of_parent") %in%
                    colnames(stats)))
  expect_true(all(c("/nonDebris", "/nonDebris/singlets") %in%
                    stats$population))
})

test_that("CollectPopulationStats reports the percentage as a percentage", {
  # gs_pop_get_stats returns a fraction, and a report that prints it as a
  # percentage without the multiplication is wrong by a factor of one hundred.
  gating_set <- MakeGatedSet()
  stats <- CollectPopulationStats(gating_set)
  expect_true(all(stats$percent_of_parent <= 100))
  expect_true(any(stats$percent_of_parent > 1))
})

test_that("CollectPopulationStats joins a sample sheet when it is given", {
  gating_set <- MakeGatedSet()
  sheet <- data.frame(file_name = "sample", condition = "control",
                      stringsAsFactors = FALSE)
  stats <- CollectPopulationStats(gating_set, sample_sheet = sheet)
  expect_true("condition" %in% colnames(stats))
  expect_true(all(stats$condition == "control"))
})
