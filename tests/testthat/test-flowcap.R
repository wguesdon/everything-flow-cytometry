# Tests for R/flowcap.R

test_that("AreaUnderCurve separates two classes that do not overlap", {
  expect_equal(AreaUnderCurve(c(1, 2, 3, 4), c(FALSE, FALSE, TRUE, TRUE)), 1)
  expect_equal(AreaUnderCurve(c(4, 3, 2, 1), c(FALSE, FALSE, TRUE, TRUE)), 0)
})

test_that("AreaUnderCurve is a half when the two classes are identical", {
  expect_equal(AreaUnderCurve(c(1, 1, 1, 1), c(TRUE, TRUE, FALSE, FALSE)), 0.5)
})

test_that("AreaUnderCurve refuses a single class", {
  expect_true(is.na(AreaUnderCurve(c(1, 2, 3), c(TRUE, TRUE, TRUE))))
})

test_that("ClassificationScore reports the four benchmark measures", {
  perfect <- ClassificationScore(c(TRUE, TRUE, FALSE, FALSE),
                                 c(TRUE, TRUE, FALSE, FALSE))
  expect_equal(perfect$recall, 1)
  expect_equal(perfect$precision, 1)
  expect_equal(perfect$accuracy, 1)
  expect_equal(perfect$f_measure, 1)
  expect_equal(perfect$samples, 4L)
})

test_that("ClassificationScore handles a prediction that misses every positive",
          {
  result <- ClassificationScore(c(FALSE, FALSE), c(TRUE, TRUE))
  expect_equal(result$recall, 0)
  expect_equal(result$precision, 0)
  expect_equal(result$f_measure, 0)
  expect_equal(result$accuracy, 0)
})

test_that("ClassificationScore returns NA for an empty comparison", {
  result <- ClassificationScore(logical(), logical())
  expect_equal(result$samples, 0L)
  expect_true(is.na(result$accuracy))
})

test_that("ControlGatedFrequencies calls a shifted sample positive", {
  control <- cbind(IL2 = withr::with_seed(1, stats::rnorm(5000)))
  stimulated <- cbind(IL2 = withr::with_seed(2, stats::rnorm(5000) + 6))
  result <- ControlGatedFrequencies(control, stimulated)
  expect_true(result[["IL2"]] > 95)
})

test_that("ControlGatedFrequencies calls an unshifted sample negative", {
  control <- cbind(IL2 = withr::with_seed(3, stats::rnorm(5000)))
  stimulated <- cbind(IL2 = withr::with_seed(4, stats::rnorm(5000)))
  result <- ControlGatedFrequencies(control, stimulated)
  expect_true(result[["IL2"]] < 1)
})

test_that("ControlGatedFrequencies refuses a population that is too small", {
  control <- cbind(IL2 = stats::rnorm(10))
  stimulated <- cbind(IL2 = stats::rnorm(10))
  expect_true(is.na(ControlGatedFrequencies(control, stimulated)[["IL2"]]))
})

test_that("SelectAndScore picks the feature that separates the classes", {
  labels <- rep(c(TRUE, FALSE), each = 10)
  features <- data.frame(
    noise = withr::with_seed(5, stats::rnorm(20)),
    signal = rep(c(10, 0), each = 10) + withr::with_seed(6, stats::rnorm(20, 0,
                 0.1))
  )
  training <- rep(c(TRUE, FALSE), 10)
  result <- SelectAndScore(features, labels, training)
  expect_equal(result$selected, "signal")
  expect_equal(result$direction, "high")
  expect_equal(result$testing_score$accuracy, 1)
})

test_that("SelectAndScore reads a feature that runs the other way", {
  labels <- rep(c(TRUE, FALSE), each = 10)
  features <- data.frame(
    inverted = rep(c(0, 10), each = 10) +
      withr::with_seed(7, stats::rnorm(20, 0, 0.1))
  )
  training <- rep(c(TRUE, FALSE), 10)
  result <- SelectAndScore(features, labels, training)
  expect_equal(result$direction, "low")
  expect_equal(result$testing_score$accuracy, 1)
})

test_that("SelectAndScore refuses a split that is too small", {
  features <- data.frame(x = 1:6)
  labels <- rep(c(TRUE, FALSE), 3)
  expect_error(SelectAndScore(features, labels, c(rep(TRUE, 2), rep(FALSE, 4))),
               "at least four samples")
})

test_that("ChannelFor. finds the detector that carries a marker", {
  frame <- MakeTestFlowFrame()
  expect_equal(ChannelFor.(frame, "CD3"), "Ax700-A")
  expect_equal(ChannelFor.(frame, "CD4"), "PE-TxRed-A")
})

test_that("ChannelFor. ignores case and punctuation in a marker name", {
  # A panel writes CD8a, cd-8a and CD8A for one antibody.
  frame <- MakeTestFlowFrame()
  expect_equal(ChannelFor.(frame, "cd-3"), "Ax700-A")
  expect_equal(ChannelFor.(frame, "Cd3"), "Ax700-A")
})

test_that("ChannelFor. gives NA when no detector carries the marker", {
  frame <- MakeTestFlowFrame()
  expect_true(is.na(ChannelFor.(frame, "CD19")))
})

test_that("TransformedEvents. transforms only the channels it is given", {
  frame <- MakeTestFlowFrame()
  before <- flowCore::exprs(frame)
  result <- TransformedEvents.(frame, "Ax700-A", cofactor = 150)
  expect_equal(result$events[, "Ax700-A"],
               asinh(before[, "Ax700-A"] / 150))
  expect_equal(result$events[, "FSC-A"], before[, "FSC-A"])
})

test_that("TransformedEvents. records the cofactor it used", {
  frame <- MakeTestFlowFrame()
  expect_equal(TransformedEvents.(frame, "Ax700-A", cofactor = 500)$transform,
               "arcsinh/500")
})
