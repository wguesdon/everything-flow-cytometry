# Tests for R/panels.R, the helpers that any deposit can use.
#
# Every test builds its own small flowFrame. None of them reads the deposit,
# because a 138 MB file has no place in a test suite.

BuildFrame <- function(events = 500, seed = 1) {
  withr::with_seed(seed, {
    values <- cbind(
      "FSC-A" = c(runif(events - 2, 2e4, 2e5), 5e6, -3),
      "FSC-H" = c(runif(events - 2, 2e4, 2e5), 5e6, 10),
      "SSC-A" = runif(events, 1e3, 5e4),
      "V510-A" = c(rnorm(events / 2, 100, 30), rnorm(events / 2, 4000, 600)),
      "U785-A" = rnorm(events, 500, 100),
      "U450-A" = c(rnorm(events / 2, 80, 20), rnorm(events / 2, 3000, 400))
    )
  })
  parameters <- Biobase::AnnotatedDataFrame(data.frame(
    name = colnames(values),
    desc = c(NA, NA, NA, "CD3 BV510", "CD4 BUV805", "Live Dead UV Blue"),
    range = 262144, minRange = 0, maxRange = 262144,
    stringsAsFactors = FALSE,
    row.names = paste0("$P", seq_len(ncol(values)))
  ))
  flowCore::flowFrame(values, parameters)
}

test_that("MarkerTokens splits a marker name on every separator", {
  expect_equal(MarkerTokens("TCR Va7_2 BV711")[[1]],
               c("tcr", "va7", "2", "bv711"))
  expect_equal(MarkerTokens("HLA-DR PE-Cy55")[[1]],
               c("hla", "dr", "pe", "cy55"))
})

test_that("a token match keeps CD16 and CD161 apart", {
  tokens <- MarkerTokens(c("CD16 BUV496", "CD161 PE-Cy5"))
  expect_true("cd16" %in% tokens[[1]])
  expect_false("cd16" %in% tokens[[2]])
  expect_true("cd161" %in% tokens[[2]])
})

test_that("PanelScatterChannels names the three scatter channels", {
  expect_equal(
    unname(PanelScatterChannels(BuildFrame())),
    c("FSC-A", "FSC-H", "SSC-A")
  )
})

test_that("PanelScatterChannels stops when a scatter channel is absent", {
  frame <- BuildFrame()
  expect_error(
    PanelScatterChannels(frame[, c("FSC-A", "SSC-A", "V510-A")]),
    "missing the scatter channel"
  )
})

test_that("ReadCompensationState reports an identity matrix as none supplied", {
  frame <- BuildFrame()
  channels <- c("V510-A", "U785-A", "U450-A")
  identity_matrix <- diag(3)
  dimnames(identity_matrix) <- list(channels, channels)
  flowCore::keyword(frame)[["SPILL"]] <- identity_matrix
  state <- ReadCompensationState(frame)
  expect_equal(state$state, "identity, no matrix supplied")
  expect_equal(state$largest_off_diagonal, 0)
  expect_equal(state$matrix_size, 3L)
})

test_that("ReadCompensationState reports a real matrix as one to apply", {
  frame <- BuildFrame()
  channels <- c("V510-A", "U785-A", "U450-A")
  matrix_with_spill <- diag(3)
  matrix_with_spill[1, 2] <- 0.34
  dimnames(matrix_with_spill) <- list(channels, channels)
  flowCore::keyword(frame)[["SPILL"]] <- matrix_with_spill
  state <- ReadCompensationState(frame)
  expect_equal(state$state, "matrix to apply")
  expect_equal(state$largest_off_diagonal, 0.34)
})

test_that("ReadCompensationState reports a missing matrix", {
  frame <- BuildFrame()
  flowCore::keyword(frame)[["SPILL"]] <- NULL
  expect_equal(ReadCompensationState(frame)$state, "no matrix")
})

test_that("ArcsinhTransform applies arcsinh to the channels named and no others", {
  values <- matrix(c(0, 150, 1500, 10, 20, 30), ncol = 2,
                   dimnames = list(NULL, c("V510-A", "FSC-A")))
  out <- ArcsinhTransform(values, "V510-A", cofactor = 150)
  expect_equal(out[, "V510-A"], asinh(c(0, 1, 10)))
  expect_equal(out[, "FSC-A"], c(10, 20, 30))
})

test_that("ArcsinhTransform rejects a cofactor that is not positive", {
  values <- matrix(1, ncol = 1, dimnames = list(NULL, "V510-A"))
  expect_error(ArcsinhTransform(values, "V510-A", cofactor = 0),
               "cofactor must be positive")
})

test_that("ArcsinhTransform rejects a channel the matrix does not carry", {
  values <- matrix(1, ncol = 1, dimnames = list(NULL, "V510-A"))
  expect_error(ArcsinhTransform(values, "U450-A"), "missing the channel")
})

test_that("InScatterRange drops the events above the digitiser range", {
  values <- matrix(c(1e4, 3e5, -5, 2e5), ncol = 1,
                   dimnames = list(NULL, "FSC-A"))
  keep <- InScatterRange(values, c(forward_area = "FSC-A"), limit = 262144)
  expect_equal(keep, c(TRUE, FALSE, FALSE, TRUE))
})

test_that("RatioSingletMask keeps the singlets and drops the doublets", {
  # Singlets sit on one ratio with a little noise. A doublet doubles the pulse
  # area while the height stays where it was, so its ratio is halved.
  values <- withr::with_seed(4, {
    area <- c(runif(200, 90, 110), runif(5, 180, 220))
    height <- c(runif(200, 88, 92), runif(5, 88, 92))
    cbind("FSC-A" = area, "FSC-H" = height)
  })
  keep <- RatioSingletMask(values, c(forward_area = "FSC-A",
                                      forward_height = "FSC-H"))
  expect_equal(sum(keep[1:200]), 200)
  expect_equal(sum(keep[201:205]), 0)
})

test_that("RatioSingletMask keeps every event when the ratio does not vary", {
  values <- cbind("FSC-A" = rep(100, 50), "FSC-H" = rep(90, 50))
  expect_true(all(RatioSingletMask(values, c(forward_area = "FSC-A",
                                              forward_height = "FSC-H"))))
})

test_that("SplitOnChannel finds the boundary between two separated modes", {
  values <- withr::with_seed(2, {
    matrix(c(rnorm(500, 0, 0.4), rnorm(500, 8, 0.4)), ncol = 1,
           dimnames = list(NULL, "x"))
  })
  split <- SplitOnChannel(values, "x", "above")
  expect_equal(split$rule, "density")
  expect_gt(split$cut, 1)
  expect_lt(split$cut, 7)
  expect_equal(sum(split$mask), 500, tolerance = 20)
})

test_that("SplitOnChannel keeps the low side when asked to", {
  values <- withr::with_seed(2, {
    matrix(c(rnorm(500, 0, 0.4), rnorm(500, 8, 0.4)), ncol = 1,
           dimnames = list(NULL, "x"))
  })
  expect_equal(sum(SplitOnChannel(values, "x", "below")$mask), 500,
               tolerance = 20)
})

test_that("SplitOnChannel stops when the channel is absent", {
  values <- matrix(1, ncol = 1, dimnames = list(NULL, "x"))
  expect_error(SplitOnChannel(values, "y"), "no channel called 'y'")
})

test_that("AssessOneDimensionalCuts reports the share a cut would select", {
  values <- withr::with_seed(3, {
    cbind("V510-A" = c(rnorm(500, 0, 0.4), rnorm(500, 8, 0.4)),
          "U785-A" = rnorm(1000, 0, 1))
  })
  channels <- data.frame(name = c("CD3", "CD4"),
                         channel = c("V510-A", "U785-A"),
                         stringsAsFactors = FALSE)
  assessed <- AssessOneDimensionalCuts(values, channels, rep(TRUE, 1000),
                                       c("CD3", "CD4"))
  expect_equal(assessed$marker, c("CD3", "CD4"))
  expect_equal(assessed$rule[1], "density")
  expect_equal(assessed$percent_selected[1], 50, tolerance = 3)
})

test_that("WriteGatedPopulation stops when the mask is the wrong length", {
  frame <- BuildFrame(events = 100)
  expect_error(
    WriteGatedPopulation(frame, rep(TRUE, 99), tempfile(fileext = ".fcs")),
    "The mask holds 99 elements"
  )
})

test_that("WriteGatedPopulation stops when the mask selects nothing", {
  frame <- BuildFrame(events = 100)
  expect_error(
    WriteGatedPopulation(frame, rep(FALSE, 100), tempfile(fileext = ".fcs")),
    "selects no event"
  )
})

test_that("WriteGatedPopulation writes the events the mask selects", {
  frame <- BuildFrame(events = 100)
  path <- tempfile(fileext = ".fcs")
  mask <- rep(c(TRUE, FALSE), 50)
  expect_equal(WriteGatedPopulation(frame, mask, path), 50)
  written <- flowCore::read.FCS(path, truncate_max_range = FALSE)
  expect_equal(nrow(flowCore::exprs(written)), 50)
  expect_equal(flowCore::colnames(written), flowCore::colnames(frame))
})


test_that("ScatterCloudMask keeps the cloud and drops the far events", {
  values <- withr::with_seed(5, {
    cbind("FSC-A" = c(rnorm(2000, 8e4, 1e4), 2.5e5, 1e3),
          "SSC-A" = c(rnorm(2000, 3e4, 5e3), 1.2e5, 500))
  })
  keep <- ScatterCloudMask(values, c(forward_area = "FSC-A",
                                     side_area = "SSC-A"))
  expect_equal(sum(keep[1:2000]), 2000, tolerance = 60)
  expect_false(keep[2001])
  expect_false(keep[2002])
})

test_that("ScatterCloudMask keeps everything when there are too few events", {
  values <- cbind("FSC-A" = rnorm(20, 8e4, 1e4), "SSC-A" = rnorm(20, 3e4, 5e3))
  expect_true(all(ScatterCloudMask(values, c(forward_area = "FSC-A",
                                             side_area = "SSC-A"))))
})

test_that("UnstainedThreshold reads the percentile of the control", {
  values <- matrix(seq_len(1000) * 150, ncol = 1,
                   dimnames = list(NULL, "U450-A"))
  parameters <- Biobase::AnnotatedDataFrame(data.frame(
    name = "U450-A", desc = "Live Dead UV Blue", range = 262144, minRange = 0,
    maxRange = 262144, stringsAsFactors = FALSE, row.names = "$P1"
  ))
  path <- tempfile(fileext = ".fcs")
  flowCore::write.FCS(flowCore::flowFrame(values, parameters), path)
  threshold <- UnstainedThreshold(path, "U450-A", percentile = 0.5)
  expect_equal(threshold, asinh(median(values) / 150), tolerance = 1e-4)
})

test_that("UnstainedThreshold stops when the channel is absent", {
  values <- matrix(1, nrow = 10, ncol = 1, dimnames = list(NULL, "U450-A"))
  parameters <- Biobase::AnnotatedDataFrame(data.frame(
    name = "U450-A", desc = "Live Dead UV Blue", range = 262144, minRange = 0,
    maxRange = 262144, stringsAsFactors = FALSE, row.names = "$P1"
  ))
  path <- tempfile(fileext = ".fcs")
  flowCore::write.FCS(flowCore::flowFrame(values, parameters), path)
  expect_error(UnstainedThreshold(path, "V510-A"), "no channel called 'V510-A'")
})

test_that("PlotGateTree places a child under its own parent", {
  counts <- data.frame(
    population = c("root", "a", "a1", "b", "b1"),
    parent = c(NA, "root", "a", "root", "b"),
    events = c(100L, 60L, 30L, 40L, 20L),
    percent_of_parent = c(100, 60, 50, 40, 50),
    stringsAsFactors = FALSE
  )
  plot <- PlotGateTree(counts)
  nodes <- plot$layers[[2]]$data
  y <- stats::setNames(nodes$y, nodes$population)

  # The walk is root, a, a1, b, b1, so the whole of branch a sits above b.
  expect_true(y[["a"]] > y[["a1"]])
  expect_true(y[["a1"]] > y[["b"]])
  expect_true(y[["b"]] > y[["b1"]])
  expect_equal(unname(y[["root"]]), -1)
})

test_that("PlotGateTree names a population whose parent is absent", {
  counts <- data.frame(
    population = c("root", "orphan"), parent = c(NA, "missing"),
    events = c(10L, 5L), percent_of_parent = c(100, 50),
    stringsAsFactors = FALSE
  )
  expect_error(PlotGateTree(counts), "orphan")
})

MakeScatterFrame. <- function(n_events = 3000, seed = 9) {
  withr::with_seed(seed, {
    lymphocytes <- cbind(`FSC-A` = stats::rnorm(n_events * 0.7, 70000, 8000),
                         `SSC-A` = stats::rnorm(n_events * 0.7, 30000, 5000))
    granulocytes <- cbind(`FSC-A` = stats::rnorm(n_events * 0.3, 90000, 9000),
                          `SSC-A` = stats::rnorm(n_events * 0.3, 120000, 12000))
  })
  rbind(lymphocytes, granulocytes)
}

test_that("LymphocyteMask keeps the low side scatter mode", {
  values <- MakeScatterFrame.()
  mask <- LymphocyteMask(values, c(forward_area = "FSC-A",
                                   side_area = "SSC-A"))
  expect_equal(length(mask), nrow(values))
  expect_true(is.logical(mask))
  expect_gt(sum(mask), 0)
  # The kept events have to sit lower on side scatter than the ones dropped.
  expect_lt(stats::median(values[mask, "SSC-A"]),
            stats::median(values[!mask, "SSC-A"]))
})

test_that("LymphocyteMask gives the same mask twice with one seed", {
  values <- MakeScatterFrame.()
  scatter <- c(forward_area = "FSC-A", side_area = "SSC-A")
  expect_equal(LymphocyteMask(values, scatter, seed = 3),
               LymphocyteMask(values, scatter, seed = 3))
})

test_that("LymphocyteMask keeps every event when there are too few to fit", {
  values <- MakeScatterFrame.(n_events = 40)
  mask <- LymphocyteMask(values, c(forward_area = "FSC-A",
                                   side_area = "SSC-A"))
  expect_equal(length(mask), nrow(values))
})

test_that("UnstainedThresholds reads one threshold per channel", {
  directory <- withr::local_tempdir()
  path <- WriteTestFcs(directory, n_events = 500)
  channels <- data.frame(channel = c("Ax700-A", "PE-TxRed-A"),
                         name = c("CD3", "CD4"), stringsAsFactors = FALSE)
  thresholds <- UnstainedThresholds(path, channels)
  expect_equal(names(thresholds), c("CD3", "CD4"))
  expect_true(all(is.finite(thresholds)))
})

test_that("UnstainedThresholds rises with the percentile it is given", {
  directory <- withr::local_tempdir()
  path <- WriteTestFcs(directory, n_events = 500)
  channels <- data.frame(channel = "Ax700-A", name = "CD3",
                         stringsAsFactors = FALSE)
  expect_gt(UnstainedThresholds(path, channels, percentile = 0.999),
            UnstainedThresholds(path, channels, percentile = 0.5))
})

test_that("UnstainedThresholds rejects a control that lacks a channel", {
  directory <- withr::local_tempdir()
  path <- WriteTestFcs(directory)
  channels <- data.frame(channel = "APC-A", name = "CD8",
                         stringsAsFactors = FALSE)
  expect_error(UnstainedThresholds(path, channels),
               "missing the channel\\(s\\): APC-A")
})

test_that("PlotGatePair draws a biaxial plot of the parent events", {
  values <- MakeScatterFrame.()
  parent <- rep(TRUE, nrow(values))
  drawing <- PlotGatePair(values, parent, "FSC-A", "SSC-A")
  expect_s3_class(drawing, "ggplot")
})

test_that("PlotGatePair draws the threshold line it is given", {
  values <- MakeScatterFrame.()
  parent <- rep(TRUE, nrow(values))
  plain <- PlotGatePair(values, parent, "FSC-A", "SSC-A")
  with_line <- PlotGatePair(values, parent, "FSC-A", "SSC-A",
                            x_threshold = 80000)
  expect_gt(length(with_line$layers), length(plain$layers))
})

test_that("PlotGatePair rejects a channel the matrix does not carry", {
  values <- MakeScatterFrame.()
  parent <- rep(TRUE, nrow(values))
  expect_error(PlotGatePair(values, parent, "FSC-A", "APC-A"),
               "no channel called 'APC-A'")
})
