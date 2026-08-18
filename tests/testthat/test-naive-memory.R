# Tests for R/naive_memory.R

test_that("DensityCut lands between two separated modes", {
  values <- withr::with_seed(1, c(stats::rnorm(4000, 0), stats::rnorm(4000, 6)))
  cut <- DensityCut(values)
  expect_true(cut > 1 && cut < 5)
  expect_equal(round(mean(values > cut), 1), 0.5)
})

test_that("DensityCut refuses one mode rather than guessing", {
  values <- withr::with_seed(2, stats::rnorm(4000))
  expect_true(is.na(DensityCut(values)))
})

test_that("DensityCut refuses a sample that is too small", {
  expect_true(is.na(DensityCut(c(1, 2, 3))))
  expect_true(is.na(DensityCut(rep(1, 500))))
})

test_that("DensityCut walks the bandwidth ladder when one width over smooths", {
  values <- withr::with_seed(3,
    c(stats::rnorm(4000, 0, 1), stats::rnorm(4000, 3, 1))
  )
  expect_true(is.na(DensityCutAt.(values, adjust = 6, depth = 0.05)))
  expect_false(is.na(DensityCut(values, adjust = c(6, 2, 1))))
})

test_that("MixtureCut splits two components that share a shoulder", {
  values <- withr::with_seed(4,
    c(stats::rnorm(6000, 0, 1), stats::rnorm(4000, 2.5, 1))
  )
  cut <- MixtureCut(values)
  expect_true(cut > 0 && cut < 2.5)
  expect_true(abs(mean(values > cut) - 0.4) < 0.1)
})

test_that("MixtureCut refuses a constant vector", {
  expect_true(is.na(MixtureCut(rep(2, 1000))))
  expect_true(is.na(MixtureCut(c(1, 2, 3))))
})

test_that("ResolveCut names the rule that produced the cut", {
  separated <- withr::with_seed(5,
    c(stats::rnorm(4000, 0), stats::rnorm(4000, 6))
  )
  expect_equal(ResolveCut(separated)$rule, "density")

  continuous <- withr::with_seed(6, stats::rnorm(4000))
  expect_equal(ResolveCut(continuous)$rule, "mixture")

  expect_equal(ResolveCut(rep(1, 1000))$rule, "none")
})

test_that("Fraction. reports a percentage and refuses an empty parent", {
  expect_equal(Fraction.(c(TRUE, TRUE, FALSE, FALSE)), 50)
  expect_true(is.na(Fraction.(logical())))
})

test_that("StableTimeWindow keeps a run of even acquisition", {
  frame <- MakeTimeFlowFrame(n_events = 4000)
  result <- StableTimeWindow(frame, "Time")
  expect_true(result$applied)
  expect_equal(result$kept, 1)
  expect_equal(nrow(flowCore::exprs(result$frame)), 4000)
})

test_that("StableTimeWindow abandons a window that would cost too much", {
  frame <- MakeTimeFlowFrame(n_events = 4000, gap = TRUE)
  result <- StableTimeWindow(frame, "Time", min_kept = 0.9)
  expect_false(result$applied)
  expect_equal(result$kept, 1)
  expect_equal(nrow(flowCore::exprs(result$frame)), 4000)
})

test_that("FitTwoComponents. separates two modes", {
  values <- withr::with_seed(21, c(stats::rnorm(400, 1, 0.3),
                                   stats::rnorm(400, 6, 0.3)))
  fit <- FitTwoComponents.(values, iterations = 200, tolerance = 1e-6)
  expect_false(is.null(fit))
  expect_equal(length(fit$means), 2)
  expect_lt(fit$means[1], fit$means[2])
  expect_gt(fit$means[2] - fit$means[1], 3)
})

test_that("FitTwoComponents. gives NULL when a side holds too few values", {
  expect_null(FitTwoComponents.(c(1, 1, 2, 2), iterations = 50,
                                tolerance = 1e-6))
})

test_that("LymphocyteGate. gives NULL on too few events to fit", {
  events <- cbind(`FSC-A` = rnorm(50), `SSC-A` = rnorm(50))
  expect_null(LymphocyteGate.(events, c(forward_area = "FSC-A",
                                        side_area = "SSC-A")))
})

test_that("LymphocyteGate. keeps the low side scatter mode", {
  events <- withr::with_seed(22, rbind(
    cbind(`FSC-A` = stats::rnorm(700, 70000, 8000),
          `SSC-A` = stats::rnorm(700, 30000, 5000)),
    cbind(`FSC-A` = stats::rnorm(300, 90000, 9000),
          `SSC-A` = stats::rnorm(300, 120000, 12000))))
  kept <- LymphocyteGate.(events, c(forward_area = "FSC-A",
                                    side_area = "SSC-A"))
  expect_false(is.null(kept))
  expect_lt(nrow(kept), nrow(events))
  expect_lt(stats::median(kept[, "SSC-A"]),
            stats::median(events[, "SSC-A"]))
})
