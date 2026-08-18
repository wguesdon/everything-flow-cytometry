# Tests for SelectByHighestMarker in R/clustering.R

MakeMedianTable <- function(cd38) {
  data.frame(
    cluster = seq_along(cd38),
    events = rep(100, length(cd38)),
    percent_of_total = rep(100 / length(cd38), length(cd38)),
    `Comp-BUV395-A` = cd38,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
}

test_that("SelectByHighestMarker returns the top cluster", {
  table <- MakeMedianTable(c(10, 90, 50))

  expect_equal(SelectByHighestMarker(table, "Comp-BUV395-A", 1), 2)
})

test_that("SelectByHighestMarker returns clusters in descending order", {
  table <- MakeMedianTable(c(10, 90, 50, 70))

  expect_equal(SelectByHighestMarker(table, "Comp-BUV395-A", 3), c(2, 4, 3))
})

test_that("SelectByHighestMarker defaults to two clusters", {
  table <- MakeMedianTable(c(10, 90, 50))

  expect_equal(length(SelectByHighestMarker(table, "Comp-BUV395-A")), 2)
})

test_that("SelectByHighestMarker rejects a channel that is absent", {
  table <- MakeMedianTable(c(10, 90))

  expect_error(SelectByHighestMarker(table, "Comp-PE-A"),
               "is not in the expression table")
})

test_that("SelectByHighestMarker rejects an impossible count", {
  table <- MakeMedianTable(c(10, 90))

  expect_error(SelectByHighestMarker(table, "Comp-BUV395-A", 0),
               "must be 1 or more")
  expect_error(SelectByHighestMarker(table, "Comp-BUV395-A", 5),
               "only 2 clusters")
})
