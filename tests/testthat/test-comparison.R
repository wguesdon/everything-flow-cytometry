
test_that("PopulationProportions carries every metadata column across", {
  stats <- data.frame(sample = c("a.fcs", "b.fcs"), population = "T cells",
                      count = c(10, 20), percent_of_parent = c(30, 40),
                      stringsAsFactors = FALSE)
  metadata <- data.frame(sample = c("a.fcs", "b.fcs"),
                         treatment = c("control", "drug"),
                         donor = c("d1", "d2"), stringsAsFactors = FALSE)
  result <- PopulationProportions(stats, metadata)
  expect_equal(result$treatment, c("control", "drug"))
  expect_equal(result$donor, c("d1", "d2"))
  expect_true(all(result$metadata_found))
})

test_that("PopulationProportions keeps a sample the metadata never names", {
  # Dropping it silently is how a group loses a replicate.
  stats <- data.frame(sample = c("a.fcs", "orphan.fcs"), population = "T cells",
                      count = c(10, 20), percent_of_parent = c(30, 40),
                      stringsAsFactors = FALSE)
  metadata <- data.frame(sample = "a.fcs", treatment = "control",
                         stringsAsFactors = FALSE)
  result <- PopulationProportions(stats, metadata)
  expect_equal(nrow(result), 2)
  expect_equal(result$metadata_found, c(TRUE, FALSE))
  expect_true(is.na(result$treatment[2]))
})

test_that("PopulationProportions matches on the base name", {
  stats <- data.frame(sample = "/a/long/path/a.fcs", population = "T cells",
                      count = 10, percent_of_parent = 30,
                      stringsAsFactors = FALSE)
  metadata <- data.frame(sample = "a.fcs", treatment = "control",
                         stringsAsFactors = FALSE)
  expect_true(PopulationProportions(stats, metadata)$metadata_found)
})

test_that("PopulationProportions rejects a metadata row that repeats a sample", {
  stats <- data.frame(sample = "a.fcs", population = "T cells", count = 10,
                      percent_of_parent = 30, stringsAsFactors = FALSE)
  metadata <- data.frame(sample = c("a.fcs", "a.fcs"),
                         treatment = c("control", "drug"),
                         stringsAsFactors = FALSE)
  expect_error(PopulationProportions(stats, metadata), "names a sample twice")
})

test_that("PopulationProportions rejects a metadata with no sample column", {
  stats <- data.frame(sample = "a.fcs", population = "T cells", count = 10,
                      percent_of_parent = 30, stringsAsFactors = FALSE)
  metadata <- data.frame(file = "a.fcs", treatment = "control",
                         stringsAsFactors = FALSE)
  expect_error(PopulationProportions(stats, metadata),
               "has no column called 'sample'")
})

MakeProportions. <- function(per_group = 5, difference = 14) {
  withr::with_seed(11, {
    control <- 20 + stats::runif(per_group, -3, 3)
    drug <- 20 + difference + stats::runif(per_group, -3, 3)
  })
  data.frame(
    sample = paste0("s", seq_len(2 * per_group), ".fcs"),
    population = "T cells",
    percent_of_parent = c(control, drug),
    treatment = rep(c("control", "drug"), each = per_group),
    stringsAsFactors = FALSE)
}

test_that("CompareProportions runs a Wilcoxon test on two groups", {
  result <- CompareProportions(MakeProportions.(), "T cells", "treatment")
  expect_equal(result$test$test, "Wilcoxon rank sum")
  expect_equal(result$test$groups, 2L)
  expect_lt(result$test$p_value, 0.05)
  expect_equal(result$summary$samples, c(5L, 5L))
})

test_that("CompareProportions runs a Kruskal Wallis test on three groups", {
  proportions <- MakeProportions.()
  proportions$treatment <- rep(c("a", "b", "c"), length.out = nrow(proportions))
  result <- CompareProportions(proportions, "T cells", "treatment")
  expect_equal(result$test$test, "Kruskal Wallis")
  expect_equal(result$test$groups, 3L)
})

test_that("CompareProportions runs no test when a group holds one sample", {
  # This is the design every deposit of one file per condition has.
  result <- CompareProportions(MakeProportions.(per_group = 1), "T cells",
                               "treatment")
  expect_equal(result$test$test, "none")
  expect_true(is.na(result$test$p_value))
  expect_match(result$test$reason, "hold one sample")
})

test_that("CompareProportions still draws the figure with no test", {
  result <- CompareProportions(MakeProportions.(per_group = 1), "T cells",
                               "treatment")
  expect_s3_class(result$plot, "ggplot")
  expect_match(result$plot$labels$subtitle, "No test ran")
})

test_that("CompareProportions runs no test when there is one group", {
  proportions <- MakeProportions.()
  proportions$treatment <- "control"
  result <- CompareProportions(proportions, "T cells", "treatment")
  expect_equal(result$test$test, "none")
  expect_match(result$test$reason, "one group")
})

test_that("CompareProportions rejects a population that is not in the table", {
  expect_error(CompareProportions(MakeProportions.(), "B cells", "treatment"),
               "No row holds the population 'B cells'")
})

test_that("CompareProportions rejects a group column that is not there", {
  expect_error(CompareProportions(MakeProportions.(), "T cells", "arm"),
               "has no column called 'arm'")
})

test_that("CompareProportions drops a row with no group before it counts", {
  proportions <- MakeProportions.()
  proportions$treatment[1] <- NA
  result <- CompareProportions(proportions, "T cells", "treatment")
  expect_equal(sum(result$summary$samples), nrow(proportions) - 1)
})
