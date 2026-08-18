

test_that("SpectralChannels drops scatter and time", {
  frame <- MakeTestFlowFrame()
  expect_equal(SpectralChannels(frame), c("Ax700-A", "PE-TxRed-A"))
})

test_that("SpectralChannels drops a detector the caller excludes", {
  frame <- MakeTestFlowFrame()
  expect_equal(SpectralChannels(frame, exclude = "Ax700-A"), "PE-TxRed-A")
})

test_that("SpectralChannels ignores the case of an excluded name", {
  frame <- MakeTestFlowFrame()
  expect_equal(SpectralChannels(frame, exclude = "ax700-a"), "PE-TxRed-A")
})

MakeCutFile. <- function(side = "above") {
  path <- withr::local_tempfile(fileext = ".csv", .local_envir = parent.frame())
  writeLines(c("marker,parent,side,cut,source,note",
               paste0("CD3,root,", side, ",1.5,the paper,a note")), path)
  path
}

test_that("ReadGateCuts accepts above and below as a side", {
  cuts <- ReadGateCuts(MakeCutFile.("below"))
  expect_equal(cuts$side, "below")
  expect_equal(cuts$cut, 1.5)
})

test_that("ReadGateCuts rejects a side it cannot read", {
  # "high" reads to a person and not to the comparison that applies the cut.
  expect_error(ReadGateCuts(MakeCutFile.("high")),
               "must be 'above' or 'below'")
})

test_that("ReadGateCuts rejects a file with no source column", {
  path <- withr::local_tempfile(fileext = ".csv")
  writeLines(c("marker,parent,side,cut,note", "CD3,root,above,1.5,a note"),
             path)
  expect_error(ReadGateCuts(path), "missing the column\\(s\\): source")
})

test_that("ReadGateCuts rejects a path that does not exist", {
  expect_error(ReadGateCuts(file.path(tempdir(), "absent.csv")),
               "does not exist")
})

MakeCounts. <- function() {
  data.frame(cd45pos_events = 1000, cd3pos_events = 600, cd8pos_events = 200,
             cd161hi_events = 20, naive_events = 50, cm_events = 40,
             em_events = 30, emra_events = 30, nk_events = 100,
             stringsAsFactors = FALSE)
}

test_that("AddCd8Frequencies divides each count by its own parent", {
  counts <- AddCd8Frequencies(MakeCounts.())
  expect_equal(counts$cd3_percent_of_cd45, 60)
  expect_equal(counts$cd8_percent_of_cd3, 100 * 200 / 600)
  expect_equal(counts$cd161hi_percent_of_cd8, 10)
  expect_equal(counts$memory_events, 100)
  expect_equal(counts$memory_percent_of_cd8, 50)
})

test_that("AddCd8Frequencies gives NA rather than a division by zero", {
  counts <- MakeCounts.()
  counts$cd8pos_events <- 0
  result <- AddCd8Frequencies(counts)
  expect_true(is.na(result$cd161hi_percent_of_cd8))
  expect_true(is.na(result$memory_percent_of_cd8))
})

test_that("AddCd8Frequencies rejects a table missing a count", {
  counts <- MakeCounts.()
  counts$naive_events <- NULL
  expect_error(AddCd8Frequencies(counts),
               "missing the column\\(s\\): naive_events")
})

test_that("ReadYuSampleSheet reads the severity as an ordered rank", {
  path <- withr::local_tempfile(fileext = ".csv")
  writeLines(c(paste("file_name,alias_subject_id,sex,severity_rank,timepoint,",
                     "igg_result", sep = ""),
               "a.fcs,s1,F,normal,T1,negative",
               "b.fcs,s2,M,hospitalized,T1,positive"), path)
  sheet <- ReadYuSampleSheet(path)
  expect_equal(nrow(sheet), 2)
  expect_true(is.ordered(sheet$severity_rank) ||
                is.numeric(sheet$severity_rank))
})

test_that("ReadYuSampleSheet rejects a severity it does not know", {
  path <- withr::local_tempfile(fileext = ".csv")
  writeLines(c(paste("file_name,alias_subject_id,sex,severity_rank,timepoint,",
                     "igg_result", sep = ""),
               "a.fcs,s1,F,very ill,T1,negative"), path)
  expect_error(ReadYuSampleSheet(path))
})

test_that("ReadYuSampleSheet rejects a path that does not exist", {
  expect_error(ReadYuSampleSheet(file.path(tempdir(), "absent.csv")),
               "does not exist")
})
