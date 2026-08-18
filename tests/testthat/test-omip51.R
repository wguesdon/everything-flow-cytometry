# Tests for R/omip51.R.
#
# Every test builds its own small flowFrame. None of them reads the deposit.

BuildOmip51Frame <- function(events = 200, seed = 1) {
  withr::with_seed(seed, {
    values <- matrix(rnorm(events * 5, 500, 100), nrow = events)
  })
  colnames(values) <- c("FSC-A", "FSC-H", "SSC-A", "B780-A", "R730-A")
  parameters <- Biobase::AnnotatedDataFrame(data.frame(
    name = colnames(values),
    desc = c(NA, NA, NA, "CD27 or IgD BB790", "IgD Ax700 or CD27 APC-R700"),
    range = 262144, minRange = 0, maxRange = 262144,
    stringsAsFactors = FALSE, row.names = paste0("$P", 1:5)
  ))
  flowCore::flowFrame(values, parameters)
}

kAmbiguous <- data.frame(
  name = c("CD27", "IgD"),
  fluorochrome = c("apc r700", "bb790"),
  antibody = c("cd27", "igd"),
  stringsAsFactors = FALSE
)

test_that("the fluorochrome resolves two markers that share an antibody name", {
  resolved <- ResolveOmip51Channels(BuildOmip51Frame(), kAmbiguous)
  expect_equal(resolved$channel[resolved$name == "CD27"], "R730-A")
  expect_equal(resolved$channel[resolved$name == "IgD"], "B780-A")
  expect_true(all(resolved$matched_by == "fluorochrome"))
})

test_that("an antibody match alone would resolve both markers to one channel", {
  # This is the failure the fluorochrome rule exists to prevent. Both panel
  # entries carry the token cd27, so a token match is ambiguous.
  tokens <- MarkerTokens(c("CD27 or IgD BB790", "IgD Ax700 or CD27 APC-R700"))
  expect_true(all(vapply(tokens, function(set) "cd27" %in% set, logical(1))))
  expect_true(all(vapply(tokens, function(set) "igd" %in% set, logical(1))))
})

test_that("a longer fluorochrome name is not resolved by a shorter one", {
  frame <- BuildOmip51Frame()
  markers <- data.frame(name = c("CD27", "IgD"),
                        fluorochrome = c("apc", "bb790"),
                        antibody = c("cd27", "igd"), stringsAsFactors = FALSE)
  resolved <- ResolveOmip51Channels(frame, markers)
  expect_equal(resolved$channel[resolved$name == "CD27"], "R730-A")
})

test_that("ResolveOmip51Channels stops when a marker resolves no channel", {
  markers <- data.frame(name = "CD4", fluorochrome = "buv805",
                        antibody = "cd4", stringsAsFactors = FALSE)
  expect_error(ResolveOmip51Channels(BuildOmip51Frame(), markers),
               "resolved 0 channels")
})

test_that("ResolveOmip51Channels stops on a frame with no named marker", {
  values <- matrix(1, nrow = 10, ncol = 1, dimnames = list(NULL, "FSC-A"))
  parameters <- Biobase::AnnotatedDataFrame(data.frame(
    name = "FSC-A", desc = NA, range = 262144, minRange = 0, maxRange = 262144,
    stringsAsFactors = FALSE, row.names = "$P1"
  ))
  expect_error(
    ResolveOmip51Channels(flowCore::flowFrame(values, parameters), kAmbiguous),
    "carries no named marker"
  )
})

test_that("MatchOmip51Controls sends each control to its own fluorochrome", {
  frames <- list(
    "Comp_Cells_IgD BB790_A1_A01_001.fcs" = BuildOmip51Frame(),
    "Comp_Cells_CD27 APC-R700_A2_A02_002.fcs" = BuildOmip51Frame()
  )
  matched <- MatchOmip51Controls(flowCore::flowSet(frames))
  expect_equal(matched$channel[grepl("IgD", matched$stain)], "B780-A")
  expect_equal(matched$channel[grepl("CD27", matched$stain)], "R730-A")
  expect_true(all(matched$matched_by == "fluorochrome"))
})

test_that("MatchOmip51Controls gives the two controls different channels", {
  frames <- list(
    "Comp_Cells_IgD BB790_A1_A01_001.fcs" = BuildOmip51Frame(),
    "Comp_Cells_CD27 APC-R700_A2_A02_002.fcs" = BuildOmip51Frame()
  )
  matched <- MatchOmip51Controls(flowCore::flowSet(frames))
  expect_equal(length(unique(matched$channel)), 2L)
})

test_that("ComputeOmip51Spillover stops when no file matches the pattern", {
  expect_error(
    ComputeOmip51Spillover(tempdir(), "^NoSuchPrefix",
                           tempfile(fileext = ".csv")),
    "No control file in"
  )
})

test_that("GateOmip51File stops when the matrix names a missing detector", {
  path <- tempfile(fileext = ".fcs")
  flowCore::write.FCS(BuildOmip51Frame(), path)
  wrong <- diag(2)
  dimnames(wrong) <- list(c("B780-A", "B999-A"), c("B780-A", "B999-A"))
  expect_error(GateOmip51File(path, wrong),
               "detector\\(s\\) that the file does not carry")
})
