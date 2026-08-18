# Tests for AnnotateClusters in R/clustering.R.
#
# SelectByHighestMarker, from the same file, is tested in
# test-select-highest.R.

MakeExpression <- function() {
  data.frame(
    cluster = 1:3,
    events = c(100L, 200L, 300L),
    percent_of_total = c(16.667, 33.333, 50),
    A = c(10, 0, 5),
    B = c(5, 0, 10),
    C = c(5, 0, 10),
    stringsAsFactors = FALSE
  )
}

MakeDefinitions <- function() {
  data.frame(
    cell_type = c("short", "long"),
    A = c("pos", "pos"),
    B = c("", "pos"),
    C = c("", "pos"),
    stringsAsFactors = FALSE
  )
}

test_that("AnnotateClusters divides by the total weight and does not add up", {
  annotation <- AnnotateClusters(MakeExpression(), MakeDefinitions())

  # Each column is scaled from 0 to 1 across the clusters, so cluster 1 scores
  # 1 on A and 0.5 on B and C. The short definition names A alone and scores 1.
  # The long definition names three markers and scores (1 + 0.5 + 0.5) / 3.
  expect_equal(annotation$cell_type[1], "short")
  expect_equal(annotation$score[1], 1)
  expect_equal(annotation$margin[1], 1 - 2 / 3)

  # Under a sum the long definition would score 2 against 1 and would take this
  # cluster for naming more markers rather than for fitting it better. The
  # division is what stops that.
  expect_gt(1, 2 / 3)
})

test_that("AnnotateClusters still lets a long definition win on merit", {
  annotation <- AnnotateClusters(MakeExpression(), MakeDefinitions())
  expect_equal(annotation$cell_type[3], "long")
  expect_equal(annotation$score[3], (0.5 + 1 + 1) / 3)
  expect_equal(annotation$runner_up[3], "short")
})

test_that("AnnotateClusters doubles the weight of a high expectation", {
  definitions <- data.frame(
    cell_type = c("plain", "emphatic"),
    A = c("pos", "high"),
    stringsAsFactors = FALSE
  )
  annotation <- AnnotateClusters(MakeExpression(), definitions)

  # 'high' contributes twice the scaled value and carries twice the weight, so
  # the score is the same as 'pos' and the label is decided by the order.
  expect_equal(annotation$score[1], 1)
  expect_equal(annotation$margin[1], 0)
})

test_that("AnnotateClusters carries the cluster identity through unchanged", {
  annotation <- AnnotateClusters(MakeExpression(), MakeDefinitions())
  expect_equal(annotation$cluster, 1:3)
  expect_equal(annotation$events, c(100L, 200L, 300L))
  expect_equal(nrow(annotation), 3)
})

test_that("AnnotateClusters names the markers when none of them match", {
  definitions <- data.frame(
    cell_type = "unrelated", CD99 = "pos", stringsAsFactors = FALSE
  )
  expect_error(AnnotateClusters(MakeExpression(), definitions),
               "share no marker column")
  expect_error(AnnotateClusters(MakeExpression(), definitions), "CD99")
})

test_that("AnnotateClusters gives every cluster the same score on a flat marker", {
  expression <- MakeExpression()
  expression$A <- c(4, 4, 4)
  definitions <- data.frame(cell_type = c("only", "other"), A = c("pos", "neg"),
                            stringsAsFactors = FALSE)
  annotation <- AnnotateClusters(expression, definitions)

  # A column with no spread scales to 0.5 everywhere, so 'pos' and 'neg' tie.
  expect_equal(unique(annotation$score), 0.5)
  expect_equal(unique(annotation$margin), 0)
})
