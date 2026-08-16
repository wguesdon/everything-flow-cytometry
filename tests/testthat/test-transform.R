# Tests for R/transform.R

test_that("FluorescenceChannels keeps a marker channel and drops scatter", {
  frame <- MakeTestFlowFrame()
  result <- FluorescenceChannels(frame)

  expect_equal(result, c("Ax700-A", "PE-TxRed-A"))
  expect_false("FSC-A" %in% result)
  expect_false("FSC-H" %in% result)
})

test_that("FluorescenceChannels drops a channel named in exclude", {
  frame <- MakeTestFlowFrame()
  result <- FluorescenceChannels(frame, exclude = "PE-TxRed-A")

  expect_equal(result, "Ax700-A")
})

test_that("FluorescenceChannels ignores case in exclude", {
  frame <- MakeTestFlowFrame()
  result <- FluorescenceChannels(frame, exclude = "pe-txred-a")

  expect_equal(result, "Ax700-A")
})

test_that("ApplyLogicleTransform returns the data and the transform it used", {
  frame <- MakeTestFlowFrame()
  result <- ApplyLogicleTransform(frame)

  expect_named(result, c("data", "transform", "channels"))
  expect_s4_class(result$data, "flowFrame")
  expect_equal(result$channels, c("Ax700-A", "PE-TxRed-A"))
})

test_that("ApplyLogicleTransform compresses the fluorescence range", {
  frame <- MakeTestFlowFrame()
  before <- flowCore::exprs(frame)
  after <- flowCore::exprs(ApplyLogicleTransform(frame)$data)

  # A logicle transform maps a raw range of thousands onto a scale of a few
  # units, so the spread of the transformed channel is far smaller.
  expect_lt(diff(range(after[, "Ax700-A"])), diff(range(before[, "Ax700-A"])))
  # Scatter is untouched, because it is not a fluorescence channel.
  expect_equal(before[, "FSC-A"], after[, "FSC-A"])
})

test_that("ApplyLogicleTransform accepts an explicit channel list", {
  frame <- MakeTestFlowFrame()
  result <- ApplyLogicleTransform(frame, channels = "Ax700-A")

  expect_equal(result$channels, "Ax700-A")
  after <- flowCore::exprs(result$data)
  before <- flowCore::exprs(frame)
  expect_equal(before[, "PE-TxRed-A"], after[, "PE-TxRed-A"])
})

test_that("ApplyLogicleTransform rejects a channel that is not in the data", {
  frame <- MakeTestFlowFrame()

  expect_error(
    ApplyLogicleTransform(frame, channels = "APC-A"),
    "not in the data"
  )
})

test_that("ApplyLogicleTransform rejects an empty channel list", {
  frame <- MakeTestFlowFrame()

  expect_error(
    ApplyLogicleTransform(frame, channels = character()),
    "nothing to transform"
  )
})

test_that("ApplyLogicleTransform rejects an object of the wrong class", {
  expect_error(ApplyLogicleTransform(1:10), "must be a flowFrame or a flowSet")
})
