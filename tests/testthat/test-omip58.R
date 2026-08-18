# Tests for R/omip58.R.
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

test_that("ResolveOmip58Channels names the detector of every marker asked for", {
  frame <- BuildFrame()
  resolved <- ResolveOmip58Channels(
    frame, c(CD3 = "cd3", CD4 = "cd4", viability = "live")
  )
  expect_equal(resolved$channel, c("V510-A", "U785-A", "U450-A"))
  expect_equal(resolved$name, c("CD3", "CD4", "viability"))
})

test_that("ResolveOmip58Channels stops when a token matches nothing", {
  expect_error(
    ResolveOmip58Channels(BuildFrame(), c(CD8 = "cd8")),
    "No marker of the panel carries the token 'cd8'"
  )
})

test_that("ResolveOmip58Channels stops when a token matches twice", {
  values <- matrix(1, nrow = 10, ncol = 2,
                   dimnames = list(NULL, c("V510-A", "U390-A")))
  parameters <- Biobase::AnnotatedDataFrame(data.frame(
    name = colnames(values), desc = c("CD3 BV510", "CD3 BUV395"),
    range = 262144, minRange = 0, maxRange = 262144,
    stringsAsFactors = FALSE, row.names = c("$P1", "$P2")
  ))
  expect_error(
    ResolveOmip58Channels(flowCore::flowFrame(values, parameters),
                          c(CD3 = "cd3")),
    "matches more than one marker"
  )
})

test_that("GateOmip58File stops when the matrix names a missing detector", {
  path <- tempfile(fileext = ".fcs")
  flowCore::write.FCS(BuildFrame(), path)
  wrong <- diag(2)
  dimnames(wrong) <- list(c("V510-A", "B999-A"), c("V510-A", "B999-A"))
  expect_error(GateOmip58File(path, wrong), "detector\\(s\\) that the file does not carry")
})
