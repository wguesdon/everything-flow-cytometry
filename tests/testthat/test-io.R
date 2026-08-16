# Tests for R/io.R

test_that("ParseTutorialFileNames reads the design out of a file name", {
  result <- ParseTutorialFileNames("LD1_NS+PI_C01_exp.fcs")

  expect_equal(nrow(result), 1)
  expect_equal(result$donor, "LD1")
  expect_equal(result$condition, "NS+PI")
  expect_equal(result$stim_1, "NS")
  expect_equal(result$stim_2, "PI")
  expect_equal(result$well, "C01")
  expect_true(result$stimulated)
})

test_that("ParseTutorialFileNames marks an unstimulated sample as FALSE", {
  result <- ParseTutorialFileNames("LD2_NS+NS_A02_exp.fcs")

  expect_false(result$stimulated)
  expect_equal(result$donor, "LD2")
})

test_that("ParseTutorialFileNames handles the full eight file design", {
  files <- c(
    "LD1_NS+NS_A01_exp.fcs", "LD1_NS+PI_C01_exp.fcs",
    "LD1_PI+NS_B01_exp.fcs", "LD1_PI+PI_D01_exp.fcs",
    "LD2_NS+NS_A02_exp.fcs", "LD2_NS+PI_C02_exp.fcs",
    "LD2_PI+NS_B02_exp.fcs", "LD2_PI+PI_D02_exp.fcs"
  )
  result <- ParseTutorialFileNames(files)

  expect_equal(nrow(result), 8)
  expect_equal(sort(unique(result$donor)), c("LD1", "LD2"))
  expect_equal(length(unique(result$condition)), 4)
  # Six of the eight files carry PI in one position or the other.
  expect_equal(sum(result$stimulated), 6)
})

test_that("ParseTutorialFileNames strips a leading path", {
  result <- ParseTutorialFileNames("/some/where/LD1_NS+NS_A01_exp.fcs")

  expect_equal(result$file_name, "LD1_NS+NS_A01_exp.fcs")
  expect_equal(result$donor, "LD1")
})

test_that("ParseTutorialFileNames rejects a name that does not match", {
  expect_error(
    ParseTutorialFileNames("Specimen_001_Tube_001.fcs"),
    "do not follow the tutorial pattern"
  )
})

test_that("ParseTutorialFileNames rejects an empty input", {
  expect_error(ParseTutorialFileNames(character()), "file_names is empty")
})

test_that("ReadTutorialFlowSet rejects a folder that does not exist", {
  expect_error(
    ReadTutorialFlowSet("/no/such/folder"),
    "The folder does not exist"
  )
})

test_that("ReadTutorialFlowSet rejects a folder with no matching file", {
  empty_dir <- withr::local_tempdir()

  expect_error(
    ReadTutorialFlowSet(empty_dir),
    "No file in .* matches the pattern"
  )
})

test_that("DescribeChannels separates a marker from a scatter channel", {
  frame <- MakeTestFlowFrame()
  result <- DescribeChannels(frame)

  expect_equal(nrow(result), 4)
  expect_equal(result$channel, c("FSC-A", "FSC-H", "Ax700-A", "PE-TxRed-A"))
  expect_equal(result$is_marker, c(FALSE, FALSE, TRUE, TRUE))
  expect_equal(result$marker[3], "CD3")
})

test_that("ChannelForMarker finds the detector for a marker", {
  frame <- MakeTestFlowFrame()

  expect_equal(ChannelForMarker(frame, "CD3"), "Ax700-A")
  expect_equal(ChannelForMarker(frame, "CD4"), "PE-TxRed-A")
})

test_that("ChannelForMarker ignores case and white space", {
  frame <- MakeTestFlowFrame()

  expect_equal(ChannelForMarker(frame, "cd3"), "Ax700-A")
  expect_equal(ChannelForMarker(frame, "  CD3  "), "Ax700-A")
})

test_that("ChannelForMarker names the available markers when it fails", {
  frame <- MakeTestFlowFrame()

  expect_error(ChannelForMarker(frame, "CD19"), "No channel carries the marker")
  expect_error(ChannelForMarker(frame, "CD19"), "CD3")
})
