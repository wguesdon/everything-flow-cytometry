# Tests for R/harmonisation.R

test_that("CoefficientOfVariation returns the standard deviation over the mean", {
  expect_equal(CoefficientOfVariation(c(10, 10, 10)), 0)
  expect_equal(CoefficientOfVariation(c(8, 10, 12)), stats::sd(c(8, 10, 12)) / 10)
})

test_that("CoefficientOfVariation refuses a case it cannot answer", {
  expect_true(is.na(CoefficientOfVariation(5)))
  expect_true(is.na(CoefficientOfVariation(c(-1, 1))))
  expect_true(is.na(CoefficientOfVariation(c(NA, NA))))
})

test_that("RelativeBias is the absolute relative difference from the reference", {
  expect_equal(RelativeBias(c(9, 11), reference = 10), 0)
  expect_equal(RelativeBias(c(12, 12), reference = 10), 0.2)
  expect_equal(RelativeBias(c(8, 8), reference = 10), 0.2)
})

test_that("RelativeBias refuses a zero reference", {
  expect_true(is.na(RelativeBias(c(1, 2), reference = 0)))
  expect_true(is.na(RelativeBias(numeric(), reference = 10)))
})

test_that("ZScore counts standard deviations from the reference mean", {
  expect_equal(ZScore(c(12, 8, 10), mu = 10, sigma = 2), c(1, -1, 0))
})

test_that("ZScore returns NA when the reference has no spread", {
  expect_true(all(is.na(ZScore(c(1, 2), mu = 1, sigma = 0))))
})

test_that("IntraclassCorrelation is near one when donors differ and analysts agree", {
  donor <- rep(c("d1", "d2", "d3"), each = 6)
  analyst <- rep(paste0("a", 1:6), times = 3)
  value <- rep(c(10, 40, 70), each = 6) + rep(c(0.1, -0.1), times = 9)
  icc <- IntraclassCorrelation(value, analyst, donor)
  expect_true(icc > 0.99)
})

test_that("IntraclassCorrelation is near zero when analysts differ and donors agree", {
  donor <- rep(c("d1", "d2", "d3"), each = 6)
  analyst <- rep(paste0("a", 1:6), times = 3)
  value <- rep(c(10, 25, 40, 55, 70, 85), times = 3)
  icc <- IntraclassCorrelation(value, analyst, donor)
  expect_true(icc < 0.05)
})

test_that("IntraclassCorrelation refuses a design it cannot fit", {
  expect_true(is.na(IntraclassCorrelation(c(1, 2, 3), c("a", "a", "a"),
                                          c("d1", "d2", "d3"))))
  expect_true(is.na(IntraclassCorrelation(rep(5, 12),
                                          rep(c("a", "b"), 6),
                                          rep(c("d1", "d2"), each = 6))))
})

test_that("MarkerFromLabel. reads the marker out of a fluorochrome label", {
  expect_equal(MarkerFromLabel.("CD8 PerCP-Cy5-5-A"), "CD8")
  expect_equal(MarkerFromLabel.("CCR7 PE"), "CCR7")
  expect_equal(MarkerFromLabel.("FL8 NIR"), "NIR")
  expect_equal(MarkerFromLabel.("CD3PC7"), "CD3")
  expect_true(is.na(MarkerFromLabel.("")))
  expect_true(is.na(MarkerFromLabel.("Time")))
})

test_that("MarkerFromLabel. keeps CD45RA, CD45 and CD4 apart", {
  expect_equal(MarkerFromLabel.("CD45RA APC"), "CD45RA")
  expect_equal(MarkerFromLabel.("CD45 APC H7"), "CD45")
  expect_equal(MarkerFromLabel.("CD4 FITC"), "CD4")
})

test_that("NormaliseDetector. makes two spellings of one detector compare equal", {
  expect_equal(NormaliseDetector.("PerCP-Cy5-5-A"),
               NormaliseDetector.("PerCP-Cy5.5-A"))
  expect_equal(NormaliseDetector.("APC-Cy7-A"), "apccy7")
})

test_that("ReadZ282Panel rejects a table that lacks a column", {
  path <- tempfile(fileext = ".csv")
  on.exit(unlink(path), add = TRUE)
  utils::write.csv(data.frame(marker = "CD3"), path, row.names = FALSE)
  expect_error(ReadZ282Panel(path), "lacks these columns")
  expect_error(ReadZ282Panel(file.path(tempdir(), "absent.csv")),
               "does not exist")
})

test_that("SplitLabel. breaks a label into its lowercase words", {
  expect_equal(SplitLabel.("CD3 FITC"), c("cd3", "fitc"))
  expect_equal(SplitLabel.("CD8-BV421"), c("cd8", "bv421"))
})

test_that("SplitLabel. returns nothing for an empty or missing label", {
  expect_equal(SplitLabel.(NA_character_), character())
  expect_equal(SplitLabel.(""), character())
})

test_that("ResolveMarkerChannels rejects a material it does not know", {
  frame <- MakeTestFlowFrame()
  panel <- data.frame(marker = "CD3", material = "both",
                      stringsAsFactors = FALSE)
  expect_error(ResolveMarkerChannels(frame, panel, "saliva"),
               "must be PBMC or WB")
})

test_that("ResolveMarkerChannels rejects a panel with nothing for the material", {
  frame <- MakeTestFlowFrame()
  panel <- data.frame(marker = "CD3", material = "WB",
                      stringsAsFactors = FALSE)
  expect_error(ResolveMarkerChannels(frame, panel, "PBMC"),
               "names no marker for the material")
})

test_that("ResolveMarkerChannels finds the channel that carries a marker", {
  frame <- MakeTestFlowFrame()
  panel <- data.frame(marker = c("CD3", "CD4"), material = "both",
                      stringsAsFactors = FALSE)
  resolved <- ResolveMarkerChannels(frame, panel, "PBMC")
  expect_equal(resolved$marker, c("CD3", "CD4"))
  expect_equal(resolved$channel, c("Ax700-A", "PE-TxRed-A"))
})

test_that("SettleCompensation says no matrix when the file carries none", {
  frame <- MakeTestFlowFrame()
  settled <- SettleCompensation(frame)
  expect_equal(settled$state, "no matrix")
})

test_that("SettleCompensation applies a matrix that has not been applied", {
  frame <- MakeTestFlowFrameWithSpillover(spill = 0.2)
  settled <- SettleCompensation(frame)
  expect_equal(settled$state, "matrix applied")
  # Compensation subtracts the spill, so the second channel has to drop.
  before <- mean(flowCore::exprs(frame)[, "PE-TxRed-A"])
  after <- mean(flowCore::exprs(settled$frame)[, "PE-TxRed-A"])
  expect_lt(after, before)
})

test_that("SettleCompensation leaves a frame the instrument already compensated", {
  # APPLY COMPENSATION = TRUE beside a real matrix means the values are done.
  frame <- MakeTestFlowFrameWithSpillover(spill = 0.2)
  flowCore::keyword(frame)[["APPLY COMPENSATION"]] <- "TRUE"
  settled <- SettleCompensation(frame)
  expect_equal(settled$state, "already compensated")
  expect_equal(flowCore::exprs(settled$frame), flowCore::exprs(frame))
})
