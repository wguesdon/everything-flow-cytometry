# Tests for R/bone_marrow.R

test_that("ParseBoneMarrowFileNames reads the panel and the sample", {
  result <- ParseBoneMarrowFileNames(c(
    "2-13-17 B cell Panel_B_E_C05_004.fcs",
    "3-2-17 NK panel_NK_A_B08_016.fcs",
    "Unstained_Unstained_Ck_A01_005.fcs",
    "3-17-17 Monocyte Panel_Mono_Sk,5c,_E03_003.fcs"
  ))
  expect_equal(result$panel, c("B", "NK", "Unstained", "Mono"))
  expect_equal(result$sample, c("E", "A", "Ck", "Sk"))
})

test_that("ParseBoneMarrowFileNames returns NA for a name it cannot read", {
  result <- ParseBoneMarrowFileNames("something_else.fcs")
  expect_true(is.na(result$panel))
  expect_true(is.na(result$sample))
})

test_that("InPolygon marks the points inside a square", {
  square_x <- c(0, 10, 10, 0)
  square_y <- c(0, 0, 10, 10)
  inside <- InPolygon(c(5, 15, 5, -1), c(5, 5, 15, 5), square_x, square_y)
  expect_equal(inside, c(TRUE, FALSE, FALSE, FALSE))
})

test_that("InPolygon refuses a shape that is not a polygon", {
  expect_error(InPolygon(1, 1, c(0, 1), c(0, 1)), "at least three vertices")
})

test_that("LiveMask reads the direction from the deposited dead cell gate", {
  gates <- data.frame(
    gate = rep("dead_lymphocytes", 4),
    x_channel = "SSC-H", y_channel = "V545-A",
    vertex = 1:4, x = c(1, 2, 3, 4), y = c(-200, 2500, 2500, -200),
    stringsAsFactors = FALSE
  )
  expect_equal(LiveMask(c(0, 1000, 3000), gates), c(FALSE, FALSE, TRUE))
})

test_that("LiveMask refuses a gate table with no dead cell gate", {
  gates <- data.frame(gate = "lymphocytes", y = 1, stringsAsFactors = FALSE)
  expect_error(LiveMask(1, gates), "names no dead cell gate")
})

test_that("ReadBoneMarrowPanels and ReadBoneMarrowDonors reject a bad table", {
  path <- tempfile(fileext = ".csv")
  on.exit(unlink(path), add = TRUE)
  utils::write.csv(data.frame(panel = "T"), path, row.names = FALSE)
  expect_error(ReadBoneMarrowPanels(path), "lacks these columns")
  expect_error(ReadBoneMarrowDonors(path), "lacks these columns")
  expect_error(ReadBoneMarrowPanels(file.path(tempdir(), "absent.csv")),
               "does not exist")
})

test_that("QuadrantFrequencies. splits a square into four equal quadrants", {
  events <- cbind(
    ra = c(2, 2, 0, 0),
    ccr7 = c(2, 0, 0, 2)
  )
  result <- QuadrantFrequencies.(events, "ra", "ccr7", 1, 1)
  expect_equal(unname(result), c(25, 25, 25, 25))
  expect_equal(names(result), c("N", "CM", "EM", "TEMRA"))
})

test_that("QuadrantFrequencies. returns NA for an empty subset", {
  events <- matrix(numeric(0), ncol = 2, dimnames = list(NULL, c("ra", "ccr7")))
  expect_true(all(is.na(QuadrantFrequencies.(events, "ra", "ccr7", 1, 1))))
})

test_that("ResolveBoneMarrowChannels falls back on the published detector", {
  frame <- MakeTestFlowFrame()
  panels <- data.frame(
    panel = "X", detector = c("Ax700-A", "PE-TxRed-A"),
    fluorochrome = c("Ax700", "PE-TxRed"), marker = c("CD3", "CD99"),
    stringsAsFactors = FALSE
  )
  result <- ResolveBoneMarrowChannels(frame, panels, "X")
  expect_equal(result$channel, c("Ax700-A", "PE-TxRed-A"))
  expect_equal(result$resolved_by, c("marker", "detector"))
})

test_that("ResolveBoneMarrowChannels refuses a panel it does not know", {
  frame <- MakeTestFlowFrame()
  panels <- data.frame(
    panel = "X", detector = "Ax700-A", fluorochrome = "Ax700", marker = "CD3",
    stringsAsFactors = FALSE
  )
  expect_error(ResolveBoneMarrowChannels(frame, panels, "Y"),
               "names no marker for the panel")
})
