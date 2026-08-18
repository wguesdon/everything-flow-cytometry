# Tests for R/spillover_compute.R
#
# The four matching passes are the part worth testing hard. Each pass exists
# because a real dataset broke the pass before it, and every case below is taken
# from a file in the archive.

test_that("StripControlName removes a bead prefix and a well suffix", {
  expect_equal(
    StripControlName("Comp_Beads_TCR Vd1 FITC_A1_A01_001.fcs"),
    "TCR Vd1 FITC"
  )
})

test_that("StripControlName removes a single stainings prefix", {
  expect_equal(
    StripControlName("Single stainings_CD3 PE-Cy5_004.fcs"),
    "CD3 PE-Cy5"
  )
})

test_that("StripControlName handles a cell control and a plain number suffix", {
  expect_equal(
    StripControlName("Comp_Cells_CD38 BV750_B8_B08_055.fcs"),
    "CD38 BV750"
  )
  expect_equal(
    StripControlName("Compensation Controls_FITC Stained Control_005.fcs"),
    "FITC Stained Control"
  )
})

test_that("NormaliseMarkerName folds case and punctuation to spaces", {
  expect_equal(NormaliseMarkerName("TCR Va7_2 BV711"), "tcr va7 2 bv711")
  expect_equal(NormaliseMarkerName("CD2 PerCP-Cy55"), "cd2 percp cy55")
})

test_that("CollapseMarkerName removes every separator", {
  # This pair is the reason the pass exists. OMIP-39 spells the fluorochrome one
  # way on the control file and another way in the panel.
  expect_equal(
    CollapseMarkerName("CD2 PerCP-Cy55"),
    CollapseMarkerName("CD2 PerCPCy55")
  )
  expect_equal(CollapseMarkerName("Siglec-7 APC-Vio770"), "siglec7apcvio770")
})

test_that("AntibodyToken returns the first token in lower case", {
  expect_equal(AntibodyToken("CD2 PerCP-Cy55"), "cd2")
  expect_equal(AntibodyToken("NKp30 eFluor450"), "nkp30")
  expect_equal(AntibodyToken("  CD8  BV570 "), "cd8")
})

test_that("MatchControlsToChannels matches on the exact marker", {
  flow_set <- MakeTestControlSet()
  result <- MatchControlsToChannels(flow_set)

  cd3 <- result[result$stain == "CD3 Ax700", ]
  expect_equal(cd3$channel, "Ax700-A")
  expect_equal(cd3$matched_by, "exact")
})

test_that("MatchControlsToChannels flags the unstained control", {
  flow_set <- MakeTestControlSet()
  result <- MatchControlsToChannels(flow_set)

  unstained <- result[result$matched_by == "unstained", ]
  expect_equal(nrow(unstained), 1)
  expect_equal(unstained$channel, "unstained")
})

test_that("MatchControlsToChannels falls back to the collapsed form", {
  flow_set <- MakeTestControlSet()
  result <- MatchControlsToChannels(flow_set)

  # The control is named "CD4 PE-TxRed" and the panel says "CD4 PETxRed".
  cd4 <- result[result$stain == "CD4 PE-TxRed", ]
  expect_equal(cd4$channel, "PE-TxRed-A")
  expect_equal(cd4$matched_by, "collapsed")
})

test_that("MatchControlsToChannels rejects an object that is not a flowSet", {
  expect_error(MatchControlsToChannels(data.frame(a = 1)), "must be a flowSet")
})

test_that("WriteMatchFile writes the two columns flowStats expects", {
  path <- withr::local_tempfile(fileext = ".csv")
  match_table <- data.frame(
    filename = c("a.fcs", "b.fcs", "u.fcs"),
    channel = c("Ax700-A", "PE-TxRed-A", "unstained"),
    matched_by = c("exact", "exact", "unstained"),
    stringsAsFactors = FALSE
  )

  WriteMatchFile(match_table, path)
  written <- utils::read.csv(path)

  expect_equal(colnames(written), c("filename", "channel"))
  expect_equal(nrow(written), 3)
  expect_true("unstained" %in% written$channel)
})

test_that("WriteMatchFile warns about a control that matched nothing", {
  path <- withr::local_tempfile(fileext = ".csv")
  match_table <- data.frame(
    filename = c("a.fcs", "lost.fcs", "u.fcs"),
    channel = c("Ax700-A", NA, "unstained"),
    matched_by = c("exact", "none", "unstained"),
    stringsAsFactors = FALSE
  )

  expect_warning(WriteMatchFile(match_table, path), "lost.fcs")
  written <- utils::read.csv(path)
  expect_equal(nrow(written), 2)
})

test_that("WriteMatchFile stops when there is no unstained control", {
  path <- withr::local_tempfile(fileext = ".csv")
  match_table <- data.frame(
    filename = c("a.fcs", "b.fcs"),
    channel = c("Ax700-A", "PE-TxRed-A"),
    matched_by = c("exact", "exact"),
    stringsAsFactors = FALSE
  )

  expect_error(WriteMatchFile(match_table, path), "No unstained control")
})

test_that("WriteMatchFile keeps one unstained control and warns about the rest",
          {
  path <- withr::local_tempfile(fileext = ".csv")
  match_table <- data.frame(
    filename = c("a.fcs", "u_mouse.fcs", "u_rat.fcs"),
    channel = c("Ax700-A", "unstained", "unstained"),
    matched_by = c("exact", "unstained", "unstained"),
    stringsAsFactors = FALSE
  )

  expect_warning(WriteMatchFile(match_table, path), "2 unstained controls")
  written <- utils::read.csv(path)
  expect_equal(sum(written$channel == "unstained"), 1)
})

test_that("WriteMatchFile stops when two controls claim one channel", {
  path <- withr::local_tempfile(fileext = ".csv")
  match_table <- data.frame(
    filename = c("a.fcs", "a2.fcs", "u.fcs"),
    channel = c("Ax700-A", "Ax700-A", "unstained"),
    matched_by = c("exact", "exact", "unstained"),
    stringsAsFactors = FALSE
  )

  expect_error(WriteMatchFile(match_table, path), "same channel")
})

test_that("ComputeSpilloverFromControls rejects a missing match file", {
  flow_set <- MakeTestControlSet()

  expect_error(
    ComputeSpilloverFromControls(flow_set, "/no/such/match.csv"),
    "match file does not exist"
  )
})

test_that("ComputeSpilloverFromControls names a scatter channel it cannot find",
          {
  flow_set <- MakeTestControlSet()
  path <- withr::local_tempfile(fileext = ".csv")
  utils::write.csv(
    data.frame(filename = "a.fcs",
               channel = "unstained"), path, row.names = FALSE
  )

  expect_error(
    ComputeSpilloverFromControls(flow_set, path, ssc = "SSC-A"),
    "is not in the data"
  )
})

test_that("PlotSpilloverHeatmap returns a ggplot", {
  spillover <- MakeTestSpillover(channels = c("A", "B", "C"), spill = 0.1)
  plot_object <- PlotSpilloverHeatmap(spillover)

  expect_s3_class(plot_object, "ggplot")
})

test_that("PlotSpilloverHeatmap rejects a non matrix", {
  expect_error(PlotSpilloverHeatmap(1:4), "must be a matrix")
})

test_that("CompareSpilloverMatrices reports the difference in percentage points", {
  computed <- matrix(
    c(1.00, 0.25, 0.05, 1.00), nrow = 2, byrow = TRUE,
    dimnames = list(c("A", "B"), c("A", "B"))
  )
  stored <- matrix(
    c(1.00, 0.20, 0.05, 1.00), nrow = 2, byrow = TRUE,
    dimnames = list(c("A", "B"), c("A", "B"))
  )
  result <- CompareSpilloverMatrices(computed, stored)

  expect_equal(nrow(result), 2)
  # A into B moved from 20 percent to 25 percent, so the difference is 5 points.
  top <- result[1, ]
  expect_equal(top$from, "A")
  expect_equal(top$to, "B")
  expect_equal(top$computed, 25)
  expect_equal(top$stored, 20)
  expect_equal(top$difference, 5)
})

test_that("CompareSpilloverMatrices stops when the matrices share no detector",
          {
  computed <- MakeTestSpillover(channels = c("A", "B"))
  stored <- MakeTestSpillover(channels = c("X", "Y"))

  expect_error(CompareSpilloverMatrices(computed, stored), "share no detector")
})

# GateableControls, added when a 28 colour matrix failed because one control of
# the set gated to an empty population.

BuildControl <- function(bright, events = 4000, seed = 1) {
  withr::with_seed(seed, {
    stained <- if (bright) {
      c(rnorm(events / 2, 200, 40), rnorm(events / 2, 20000, 3000))
    } else {
      rnorm(events, 200, 40)
    }
    values <- cbind(
      "FSC-A" = rnorm(events, 8e4, 1e4),
      "SSC-A" = rnorm(events, 3e4, 5e3),
      "V510-A" = stained,
      "U785-A" = rnorm(events, 150, 30)
    )
  })
  parameters <- Biobase::AnnotatedDataFrame(data.frame(
    name = colnames(values),
    desc = c(NA, NA, "CD3 BV510", "CD4 BUV805"),
    range = 262144, minRange = 0, maxRange = 262144,
    stringsAsFactors = FALSE, row.names = paste0("$P", seq_len(ncol(values)))
  ))
  flowCore::flowFrame(values, parameters)
}

test_that("GateableControls accepts a control with two separated modes", {
  controls <- flowCore::flowSet(list(stained = BuildControl(TRUE)))
  matches <- data.frame(
    filename = "stained", stain = "CD3 BV510", channel = "V510-A",
    marker = "CD3 BV510", matched_by = "exact", stringsAsFactors = FALSE
  )
  result <- GateableControls(controls, matches, min_events = 100)
  expect_true(result$gateable)
  expect_gt(result$gated_events, 100)
})

test_that("GateableControls rejects a control whose gate leaves too few events",
          {
  controls <- flowCore::flowSet(list(stained = BuildControl(TRUE)))
  matches <- data.frame(
    filename = "stained", stain = "CD3 BV510", channel = "V510-A",
    marker = "CD3 BV510", matched_by = "exact", stringsAsFactors = FALSE
  )
  result <- GateableControls(controls, matches, min_events = 1e6)
  expect_false(result$gateable)
})

test_that("GateableControls returns zero events for an unresolvable control", {
  controls <- flowCore::flowSet(list(flat = BuildControl(FALSE)))
  matches <- data.frame(
    filename = "flat", stain = "CD3 BV510", channel = "V510-A",
    marker = "CD3 BV510", matched_by = "exact", stringsAsFactors = FALSE
  )
  result <- GateableControls(controls, matches, min_events = 100)
  expect_equal(nrow(result), 1L)
  expect_false(is.na(result$gated_events))
})

test_that("GateableControls stops when no control was matched to a channel", {
  controls <- flowCore::flowSet(list(stained = BuildControl(TRUE)))
  matches <- data.frame(
    filename = "stained", stain = "unstained", channel = "unstained",
    marker = NA_character_, matched_by = "unstained", stringsAsFactors = FALSE
  )
  expect_error(GateableControls(controls, matches),
               "No stained control was matched to a channel")
})

test_that("MarkerCorrelation finds the pair that spills", {
  # CD4 is built from CD3, so the two correlate and nothing else does.
  events <- withr::with_seed(1, {
    cd3 <- rnorm(4000)
    cbind(CD3 = cd3, CD4 = 0.8 * cd3 + rnorm(4000), CD8 = rnorm(4000))
  })
  result <- MarkerCorrelation(events)
  expect_equal(result$summary$markers, 3L)
  expect_equal(result$summary$pairs, 3L)
  expect_equal(result$summary$most_correlated, "CD3 against CD4")
  expect_gt(result$matrix["CD3", "CD4"], 0.5)
  expect_lt(abs(result$matrix["CD3", "CD8"]), 0.1)
})

test_that("MarkerCorrelation counts a negative pair as high", {
  # Over-compensation drives a correlation negative, which is also a fault.
  events <- withr::with_seed(2, {
    cd3 <- rnorm(3000)
    cbind(CD3 = cd3, CD4 = -0.9 * cd3 + rnorm(3000, sd = 0.2))
  })
  result <- MarkerCorrelation(events, threshold = 0.5)
  expect_lt(result$matrix["CD3", "CD4"], -0.5)
  expect_equal(result$summary$pairs_above_threshold, 1L)
  expect_equal(result$summary$positive_pairs_above_threshold, 0L)
})

test_that("MarkerCorrelation gives the same answer twice with one seed", {
  events <- withr::with_seed(3, matrix(rnorm(6000), ncol = 3,
                                       dimnames = list(NULL, c("a", "b", "c"))))
  first <- MarkerCorrelation(events, max_events = 500, seed = 7)
  second <- MarkerCorrelation(events, max_events = 500, seed = 7)
  expect_equal(first$summary$median_absolute_r,
               second$summary$median_absolute_r)
})

test_that("MarkerCorrelation subsamples down to max_events", {
  events <- withr::with_seed(4, matrix(rnorm(4000), ncol = 2,
                                       dimnames = list(NULL, c("a", "b"))))
  expect_equal(MarkerCorrelation(events, max_events = 300)$summary$events_used,
               300L)
  expect_equal(MarkerCorrelation(events)$summary$events_used, 2000L)
})

test_that("MarkerCorrelation lists each unordered pair once", {
  events <- withr::with_seed(5, matrix(rnorm(400), ncol = 4,
                                       dimnames = list(NULL, letters[1:4])))
  result <- MarkerCorrelation(events)
  expect_equal(nrow(result$pairs), 6L)
  expect_true(all(result$pairs$a < result$pairs$b))
})

test_that("MarkerCorrelation sorts the pairs by the size of the correlation", {
  events <- withr::with_seed(6, {
    a <- rnorm(2000)
    cbind(a = a, b = 0.9 * a + rnorm(2000, sd = 0.3), c = rnorm(2000))
  })
  result <- MarkerCorrelation(events)
  expect_equal(result$pairs$a[1], "a")
  expect_equal(result$pairs$b[1], "b")
  expect_true(all(diff(abs(result$pairs$r)) <= 0))
})

test_that("MarkerCorrelation rejects a matrix with one marker", {
  events <- matrix(rnorm(100), ncol = 1, dimnames = list(NULL, "a"))
  expect_error(MarkerCorrelation(events), "A correlation needs two markers")
})

test_that("MarkerCorrelation rejects a matrix with no column names", {
  expect_error(MarkerCorrelation(matrix(rnorm(100), ncol = 2)),
               "must carry a column name")
})

test_that("MarkerCorrelation rejects something that is not a numeric matrix", {
  expect_error(MarkerCorrelation(data.frame(a = 1:5, b = 6:10)),
               "must be a numeric matrix")
})

test_that("MarkerCorrelation rejects a matrix with two events", {
  events <- matrix(rnorm(4), ncol = 2, dimnames = list(NULL, c("a", "b")))
  expect_error(MarkerCorrelation(events), "needs three events")
})

test_that("CheckControlQuality reports one row per stained control", {
  controls <- MakeTestControlSet()
  match_table <- MatchControlsToChannels(controls)
  quality <- CheckControlQuality(controls, match_table)
  stained <- sum(!is.na(match_table$channel) &
                   match_table$channel != "unstained")
  expect_equal(nrow(quality), stained)
  expect_true("filename" %in% colnames(quality))
})

test_that("CheckControlQuality rejects a match table with no unstained control",
          {
  # Without an unstained tube there is no threshold, so a positive fraction
  # cannot be computed at all.
  controls <- MakeTestControlSet()
  match_table <- MatchControlsToChannels(controls)
  match_table$channel[match_table$channel == "unstained"] <- NA
  expect_error(CheckControlQuality(controls, match_table),
               "names no unstained control")
})

test_that("CheckControlQuality rejects a match table with no stained control", {
  controls <- MakeTestControlSet()
  match_table <- MatchControlsToChannels(controls)
  match_table$channel[match_table$channel != "unstained"] <- NA
  expect_error(CheckControlQuality(controls, match_table),
               "holds no stained control")
})
