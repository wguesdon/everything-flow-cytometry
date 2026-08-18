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
  # Both variables carry the folder that the CLI mounted, never the file that
  # --data or --out named. One substitution then covers every case.
  withr::local_envvar(CYTOKIT_DATA_HOST = "/home/will/study",
                      CYTOKIT_OUT_HOST = "/home/will/out")
  expect_equal(DisplayPath("/indata"), "/home/will/study")
  expect_equal(DisplayPath("/indata/a.fcs"), "/home/will/study/a.fcs")
  expect_equal(DisplayPath("/outdata"), "/home/will/out")
  expect_equal(DisplayPath("/outdata/template.csv"),
               "/home/will/out/template.csv")
  expect_equal(DisplayPath("/outdata/inspect_study_2026_08_18_1200"),
               "/home/will/out/inspect_study_2026_08_18_1200")
})

test_that("DisplayPath keeps a root path that is only a slash", {
  withr::local_envvar(CYTOKIT_DATA_HOST = "/", CYTOKIT_OUT_HOST = "")
  expect_equal(DisplayPath("/indata/a.fcs"), "/a.fcs")
})

test_that("DisplayPath leaves a path alone when the mapping is unknown", {
  withr::local_envvar(CYTOKIT_DATA_HOST = "", CYTOKIT_OUT_HOST = "")
  expect_equal(DisplayPath("/indata/a.fcs"), "/indata/a.fcs")
  expect_equal(DisplayPath("relative/path.csv"), "relative/path.csv")
})

test_that("AcquisitionKind reads a mass cytometer from the instrument keyword", {
  panel <- data.frame(channel = c("FSC-A", "SSC-A", "APC-A"),
                      kind = c("scatter", "scatter", "stain"),
                      stringsAsFactors = FALSE)
  result <- AcquisitionKind(panel, "DVSSCIENCES-CYTOF-6.7.1014")
  expect_equal(result$kind, "mass")
  expect_true(grepl("CYTOF", result$reason))
})

test_that("AcquisitionKind reads a mass cytometer from the detector names", {
  # FR-FCM-Z244 names every detector after the isotope it counts.
  panel <- data.frame(
    channel = c("Y89Di", "Pd102Di", "Ir193Di", "Event_length"),
    kind = c("stain", "stain", "stain", "unnamed"),
    stringsAsFactors = FALSE)
  result <- AcquisitionKind(panel, NA_character_)
  expect_equal(result$kind, "mass")
  expect_equal(result$reason, "3 detectors are named after an isotope")
  expect_false(result$has_scatter)
})

test_that("AcquisitionKind calls a fluorescence panel fluorescence", {
  panel <- data.frame(channel = c("FSC-A", "SSC-A", "B515-A"),
                      kind = c("scatter", "scatter", "stain"),
                      stringsAsFactors = FALSE)
  result <- AcquisitionKind(panel, c("LSRFortessa", "LSRFortessa"))
  expect_equal(result$kind, "fluorescence")
  expect_true(result$has_scatter)
})

test_that("AcquisitionKind needs three isotope names, not one", {
  # A single channel called Cd45Di in a fluorescence panel is a coincidence.
  panel <- data.frame(channel = c("FSC-A", "SSC-A", "Cd45Di"),
                      kind = c("scatter", "scatter", "stain"),
                      stringsAsFactors = FALSE)
  expect_equal(AcquisitionKind(panel, NA_character_)$kind, "fluorescence")
})

test_that("AcquisitionKind reports a fluorescence panel with no scatter", {
  panel <- data.frame(channel = c("B515-A", "R780-A"),
                      kind = c("stain", "stain"), stringsAsFactors = FALSE)
  result <- AcquisitionKind(panel, "FACSAria")
  expect_equal(result$kind, "fluorescence")
  expect_false(result$has_scatter)
})

test_that("CollectNotes counts a repeated warning instead of repeating it", {
  result <- CollectNotes({
    for (index in 1:3) warning("uneven number of tokens: 695")
    42
  })
  expect_equal(result$value, 42)
  expect_equal(nrow(result$notes), 1)
  expect_equal(result$notes$times, 3L)
  expect_equal(result$notes$note, "uneven number of tokens: 695")
})

test_that("CollectNotes catches a line the reader printed rather than warned", {
  # flowCore writes some of these with cat, so a condition handler alone misses
  # them and the report fills with copies of one line.
  result <- CollectNotes({
    cat("The last keyword is dropped.\n")
    cat("The last keyword is dropped.\n")
    "panel"
  })
  expect_equal(result$value, "panel")
  expect_equal(result$notes$times, 2L)
})

test_that("CollectNotes returns an empty table when nothing was reported", {
  result <- CollectNotes(sum(1:4))
  expect_equal(result$value, 10L)
  expect_equal(nrow(result$notes), 0)
})

test_that("CollectNotes orders the notes by how often each was reported", {
  result <- CollectNotes({
    warning("once")
    for (index in 1:2) warning("twice")
    NULL
  })
  expect_equal(result$notes$note, c("twice", "once"))
  expect_equal(result$notes$times, c(2L, 1L))
})
