# Tests for R/compensation.R

test_that("ExtractSpillover returns the matrix from the SPILL keyword", {
  frame <- MakeTestFlowFrameWithSpillover(spill = 0.12)
  result <- ExtractSpillover(frame)

  expect_true(is.matrix(result))
  expect_equal(dim(result), c(2, 2))
  expect_equal(diag(result), c(`Ax700-A` = 1, `PE-TxRed-A` = 1))
  expect_equal(result["Ax700-A", "PE-TxRed-A"], 0.12)
})

test_that("ExtractSpillover fails with a usable message on a frame with none", {
  frame <- MakeTestFlowFrame()

  expect_error(ExtractSpillover(frame), "carries no spillover matrix")
  expect_error(ExtractSpillover(frame), "spillover_ng")
})

test_that("CheckSpilloverChannels passes when every detector is present", {
  frame <- MakeTestFlowFrame()
  spillover <- MakeTestSpillover()

  expect_true(CheckSpilloverChannels(frame, spillover))
})

test_that("CheckSpilloverChannels names the detectors that are missing", {
  frame <- MakeTestFlowFrame()
  spillover <- MakeTestSpillover(channels = c("FITC.A", "PE.A"))

  expect_error(CheckSpilloverChannels(frame, spillover), "FITC.A")
  expect_error(CheckSpilloverChannels(frame, spillover), "PE.A")
  expect_error(
    CheckSpilloverChannels(frame, spillover),
    "Rename the matrix columns"
  )
})

test_that("ApplyCompensation changes the values it is meant to change", {
  frame <- MakeTestFlowFrameWithSpillover(spill = 0.2)
  before <- flowCore::exprs(frame)
  compensated <- ApplyCompensation(frame)
  after <- flowCore::exprs(compensated)

  # Compensation multiplies by the inverse of the spillover matrix, so the two
  # fluorescence channels move and the two scatter channels do not.
  expect_false(isTRUE(all.equal(before[, "Ax700-A"], after[, "Ax700-A"])))
  expect_equal(before[, "FSC-A"], after[, "FSC-A"])
  expect_equal(dim(before), dim(after))
})

test_that("ApplyCompensation with an identity matrix leaves the data alone", {
  frame <- MakeTestFlowFrameWithSpillover(spill = 0)
  before <- flowCore::exprs(frame)
  after <- flowCore::exprs(ApplyCompensation(frame))

  expect_equal(before, after)
})

test_that("ApplyCompensation rejects an object of the wrong class", {
  expect_error(ApplyCompensation(data.frame(a = 1)), "must be a flowFrame")
})

test_that("SummariseSpillover drops the diagonal and sorts by size", {
  spillover <- matrix(
    c(1.0, 0.30, 0.05,
      0.10, 1.00, 0.02,
      0.01, 0.04, 1.00),
    nrow = 3, byrow = TRUE,
    dimnames = list(c("A", "B", "C"), c("A", "B", "C"))
  )
  result <- SummariseSpillover(spillover)

  expect_equal(nrow(result), 6)
  expect_false(any(result$from == result$to))
  expect_equal(result$spill[1], 0.30)
  expect_equal(result$from[1], "A")
  expect_equal(result$to[1], "B")
  expect_true(all(diff(result$spill) <= 0))
})

test_that("SummariseSpillover honours the top argument", {
  spillover <- MakeTestSpillover(channels = c("A", "B", "C"), spill = 0.1)
  result <- SummariseSpillover(spillover, top = 2)

  expect_equal(nrow(result), 2)
})

test_that("SummariseSpillover rejects a non matrix", {
  expect_error(SummariseSpillover(list(1, 2)), "must be a matrix")
})

test_that("ExtractSpillover fills in row names that the file left empty", {
  frame <- MakeTestFlowFrame()
  spillover <- MakeTestSpillover()
  rownames(spillover) <- NULL
  flowCore::keyword(frame)[["SPILL"]] <- spillover

  result <- ExtractSpillover(frame)

  expect_equal(rownames(result), colnames(result))
  expect_equal(rownames(result), c("Ax700-A", "PE-TxRed-A"))
})

test_that("SummariseSpillover works on a matrix with no row names", {
  spillover <- matrix(
    c(1.0, 0.30,
      0.10, 1.00),
    nrow = 2, byrow = TRUE,
    dimnames = list(NULL, c("A", "B"))
  )
  result <- SummariseSpillover(spillover)

  expect_equal(nrow(result), 2)
  expect_equal(result$spill[1], 0.30)
  expect_equal(result$from[1], "A")
  expect_equal(result$to[1], "B")
})
