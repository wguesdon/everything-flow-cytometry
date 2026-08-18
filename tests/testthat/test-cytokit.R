# Tests for R/cytokit.R.

test_that("ParseCytokitArguments reads a key and its value", {
  parsed <- ParseCytokitArguments(c("--data", "a.fcs", "--out", "here"),
                                  allowed = c("data", "out"))
  expect_equal(parsed$data, "a.fcs")
  expect_equal(parsed$out, "here")
  expect_length(parsed, 2)
})

test_that("ParseCytokitArguments turns a flag into TRUE", {
  parsed <- ParseCytokitArguments(c("--data", "a.fcs", "--force"),
                                  allowed = "data", flags = "force")
  expect_true(parsed$force)
  expect_equal(parsed$data, "a.fcs")
})

test_that("ParseCytokitArguments rejects an argument it does not know", {
  # A misspelled flag would otherwise change nothing and report success.
  expect_error(
    ParseCytokitArguments(c("--dta", "a.fcs"), allowed = "data"),
    "Unknown argument '--dta'"
  )
  expect_error(
    ParseCytokitArguments(c("--dta", "a.fcs"), allowed = c("data", "out")),
    "data, out"
  )
})

test_that("ParseCytokitArguments rejects a value with no name", {
  expect_error(ParseCytokitArguments("a.fcs", allowed = "data"),
               "beginning with '--'")
})

test_that("ParseCytokitArguments rejects a name with no value", {
  expect_error(ParseCytokitArguments("--data", allowed = "data"),
               "needs a value")
  expect_error(
    ParseCytokitArguments(c("--data", "--out", "here"),
                          allowed = c("data", "out")),
    "'--data' needs a value"
  )
})

test_that("ParseCytokitArguments names every argument that is missing", {
  expect_error(
    ParseCytokitArguments(c("--data", "a.fcs"), allowed = c("data", "out"),
                          required = c("data", "out")),
    "--out"
  )
})

test_that("CytokitTimestamp writes a sortable stamp", {
  stamp <- CytokitTimestamp(as.POSIXct("2026-08-18 16:30:45", tz = "UTC"))
  expect_match(stamp, "^\\d{4}_\\d{2}_\\d{2}_\\d{6}$")
  expect_equal(substr(stamp, 1, 10), "2026_08_18")
})

test_that("OpenCytokitBundle names a folder from the recipe and the input", {
  root <- withr::local_tempdir()
  bundle <- OpenCytokitBundle("inspect", "my study/", out_root = root,
                              timestamp = "2026_08_18_163045")
  expect_true(dir.exists(bundle))
  expect_equal(basename(bundle), "inspect_my_study_2026_08_18_163045")
})

test_that("OpenCytokitBundle survives a label with no usable characters", {
  root <- withr::local_tempdir()
  bundle <- OpenCytokitBundle("inspect", "///", out_root = root,
                              timestamp = "2026_08_18_163045")
  expect_equal(basename(bundle), "inspect_input_2026_08_18_163045")
})

test_that("CloseCytokitBundle records the arguments and the input checksums", {
  root <- withr::local_tempdir()
  bundle <- OpenCytokitBundle("inspect", "study", out_root = root)
  input <- file.path(root, "sample.txt")
  writeLines("one", input)

  CloseCytokitBundle(bundle, "inspect", list(data = "study"), inputs = input,
                     command = "cytokit inspect --data study")

  for (name in c("manifest.json", "session_info.txt", "REPRODUCE.md")) {
    expect_true(file.exists(file.path(bundle, name)), info = name)
  }
  manifest <- jsonlite::fromJSON(file.path(bundle, "manifest.json"))
  expect_equal(manifest$recipe, "inspect")
  expect_equal(manifest$arguments$data, "study")
  expect_equal(manifest$inputs$sample.txt,
               unname(tools::md5sum(input)))
  expect_match(paste(readLines(file.path(bundle, "REPRODUCE.md")),
                     collapse = " "),
               "cytokit inspect --data study")
})

test_that("CloseCytokitBundle leaves out an input that is not on disk", {
  root <- withr::local_tempdir()
  bundle <- OpenCytokitBundle("inspect", "study", out_root = root)
  CloseCytokitBundle(bundle, "inspect", list(), inputs = "absent.fcs")
  manifest <- jsonlite::fromJSON(file.path(bundle, "manifest.json"))
  expect_length(manifest$inputs, 0)
})

test_that("FcsFilesIn returns a single file unchanged", {
  root <- withr::local_tempdir()
  path <- file.path(root, "one.fcs")
  writeLines("not really an FCS file", path)
  expect_equal(FcsFilesIn(path), path)
})

test_that("FcsFilesIn sorts a folder and ignores anything else", {
  root <- withr::local_tempdir()
  for (name in c("b.fcs", "a.fcs", "notes.txt")) {
    writeLines("x", file.path(root, name))
  }
  found <- FcsFilesIn(root)
  expect_equal(basename(found), c("a.fcs", "b.fcs"))
})

test_that("FcsFilesIn says how to look deeper when a folder holds nothing", {
  root <- withr::local_tempdir()
  dir.create(file.path(root, "inner"))
  writeLines("x", file.path(root, "inner", "a.fcs"))
  expect_error(FcsFilesIn(root), "--recursive")
  expect_length(FcsFilesIn(root, recursive = TRUE), 1)
})

test_that("FcsFilesIn names a path that does not exist", {
  expect_error(FcsFilesIn("no/such/place"), "does not exist")
})

MakePanel <- function() {
  data.frame(
    channel = c("FSC-A", "SSC-A", "Time", "APC-A", "FITC-A", "PE-A"),
    marker = c("", "", "", "CD3", "", "CD19"),
    kind = c("scatter", "scatter", "time", "stain", "unnamed", "stain"),
    range = rep(262144, 6),
    stringsAsFactors = FALSE
  )
}

test_that("PanelMarkers falls back to the detector when a marker is missing", {
  expect_equal(PanelMarkers(MakePanel()), c("CD3", "CD19", "FITC-A"))
})

test_that("PanelMarkers can report only the detectors the file named", {
  expect_equal(PanelMarkers(MakePanel(), fallback = FALSE), c("CD3", "CD19"))
})

test_that("PanelMarkers leaves out the scatter and the time detectors", {
  expect_false(any(c("FSC-A", "SSC-A", "Time") %in% PanelMarkers(MakePanel())))
})

test_that("PanelNamingState counts what the file named", {
  state <- PanelNamingState(MakePanel())
  expect_equal(state$named, 2)
  expect_equal(state$unnamed, 1)
  expect_equal(state$state, "some detectors unnamed")
})

test_that("PanelNamingState reports a file that names nothing", {
  panel <- MakePanel()
  panel$kind[panel$kind == "stain"] <- "unnamed"
  panel$marker <- ""
  expect_equal(PanelNamingState(panel)$state, "no marker names")
})

test_that("PanelNamingState reports a file that names everything", {
  panel <- MakePanel()
  panel$kind[panel$kind == "unnamed"] <- "stain"
  panel$marker[panel$channel == "FITC-A"] <- "CD4"
  expect_equal(PanelNamingState(panel)$state, "every detector named")
})

test_that("DisplayPath maps the container path back to the one that was typed", {
  withr::local_envvar(CYTOKIT_DATA_HOST = "/home/will/study",
                      CYTOKIT_OUT_HOST = "/home/will/out/template.csv")
  expect_equal(DisplayPath("/indata"), "/home/will/study")
  expect_equal(DisplayPath("/indata/a.fcs"), "/home/will/study/a.fcs")
  expect_equal(DisplayPath("/outdata/template.csv"),
               "/home/will/out/template.csv")
})

test_that("DisplayPath leaves a path alone when the mapping is unknown", {
  withr::local_envvar(CYTOKIT_DATA_HOST = "", CYTOKIT_OUT_HOST = "")
  expect_equal(DisplayPath("/indata/a.fcs"), "/indata/a.fcs")
  expect_equal(DisplayPath("relative/path.csv"), "relative/path.csv")
})
