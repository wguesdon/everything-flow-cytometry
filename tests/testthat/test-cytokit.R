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
  withr::local_envvar(CYTOKIT_DATA_HOST = "/home/scientist/study",
                      CYTOKIT_OUT_HOST = "/home/scientist/out")
  expect_equal(DisplayPath("/indata"), "/home/scientist/study")
  expect_equal(DisplayPath("/indata/a.fcs"), "/home/scientist/study/a.fcs")
  expect_equal(DisplayPath("/outdata"), "/home/scientist/out")
  expect_equal(DisplayPath("/outdata/template.csv"),
               "/home/scientist/out/template.csv")
  expect_equal(DisplayPath("/outdata/inspect_study_2026_08_18_1200"),
               "/home/scientist/out/inspect_study_2026_08_18_1200")
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

test_that("IdentitySplits calls a keyword with two values a grouping", {
  identity <- data.frame(
    file = c("a.fcs", "b.fcs", "c.fcs", "d.fcs"),
    specimen = c("PBMC", "PBMC", "PBMC_001", "PBMC_001"),
    stringsAsFactors = FALSE)
  result <- IdentitySplits(identity)
  expect_equal(result$role, "grouping")
  expect_equal(result$distinct_values, 2L)
  expect_equal(result$values, "PBMC, PBMC_001")
})

test_that("IdentitySplits calls a value per file an identifier", {
  identity <- data.frame(file = c("a.fcs", "b.fcs", "c.fcs"),
                         tube = c("1", "2", "13"), stringsAsFactors = FALSE)
  expect_equal(IdentitySplits(identity)$role, "identifier")
})

test_that("IdentitySplits calls one value across every file a constant", {
  identity <- data.frame(file = c("a.fcs", "b.fcs"),
                         experiment = c("Expt 4", "Expt 4"),
                         stringsAsFactors = FALSE)
  expect_equal(IdentitySplits(identity)$role, "constant")
})

test_that("IdentitySplits never calls a timestamp a grouping", {
  # Two files acquired one second apart are not two treatment arms.
  identity <- data.frame(file = c("a.fcs", "b.fcs", "c.fcs"),
                         started = c("06:52:13", "06:52:13", "06:52:14"),
                         stringsAsFactors = FALSE)
  expect_equal(IdentitySplits(identity)$role, "timing")
})

test_that("IdentitySplits drops a keyword that no file carries", {
  identity <- data.frame(file = c("a.fcs", "b.fcs"),
                         project = c(NA_character_, NA_character_),
                         specimen = c("PBMC", "PBMC"),
                         stringsAsFactors = FALSE)
  result <- IdentitySplits(identity)
  expect_equal(result$keyword, "specimen")
})

test_that("IdentitySplits puts a grouping before an identifier", {
  identity <- data.frame(
    file = c("a.fcs", "b.fcs", "c.fcs"),
    tube = c("1", "2", "3"),
    specimen = c("PBMC", "PBMC", "PBMC_001"),
    stringsAsFactors = FALSE)
  expect_equal(IdentitySplits(identity)$keyword, c("specimen", "tube"))
})

test_that("IdentitySplits returns an empty table when nothing is recorded", {
  identity <- data.frame(file = "a.fcs", specimen = NA_character_,
                         stringsAsFactors = FALSE)
  expect_equal(nrow(IdentitySplits(identity)), 0)
})

test_that("PanelMarkerSource says which name came from the file", {
  panel <- data.frame(channel = c("B515-A", "APC-A", "FSC-A"),
                      marker = c("CD3", "", ""),
                      kind = c("stain", "unnamed", "scatter"),
                      stringsAsFactors = FALSE)
  result <- PanelMarkerSource(panel)
  expect_equal(result$marker, c("CD3", "APC-A"))
  expect_equal(result$source, c("antibody", "detector"))
})

test_that("ParseMarkerArgument selects from a panel that names its markers", {
  result <- ParseMarkerArgument("CD3,CD4", c("CD3", "CD4", "CD8"),
                                c("B515-A", "R780-A", "V800-A"))
  expect_equal(result$column, c("CD3", "CD4"))
  expect_equal(result$from, c("CD3", "CD4"))
})

test_that("ParseMarkerArgument reads a detector to antibody mapping", {
  result <- ParseMarkerArgument("APC-A=CD3,PE-A=CD4", c("APC-A", "PE-A"),
                                c("FSC-A", "APC-A", "PE-A"))
  expect_equal(result$column, c("CD3", "CD4"))
  expect_equal(result$from, c("APC-A", "PE-A"))
})

test_that("ParseMarkerArgument takes a selection and a mapping together", {
  result <- ParseMarkerArgument("CD3,APC-A=CD8", c("CD3", "APC-A"),
                                c("B515-A", "APC-A"))
  expect_equal(sort(result$column), c("CD3", "CD8"))
})

test_that("ParseMarkerArgument rejects a detector the panel does not carry", {
  expect_error(
    ParseMarkerArgument("XYZ-A=CD3", c("APC-A"), c("FSC-A", "APC-A")),
    "This detector is not in the panel: XYZ-A")
})

test_that("ParseMarkerArgument points at the mapping when a name is unknown", {
  # This is the message a scientist meets on a panel with no marker names, so
  # it has to say what to do and not only what went wrong.
  expect_error(
    ParseMarkerArgument("CD3", c("APC-A", "PE-A"), c("APC-A", "PE-A")),
    "give the mapping instead")
})

test_that("ParseMarkerArgument rejects a mapping with no antibody", {
  expect_error(ParseMarkerArgument("APC-A=", c("APC-A"), c("APC-A")),
               "A mapping reads detector=antibody")
})

test_that("ParseMarkerArgument rejects the same column twice", {
  expect_error(
    ParseMarkerArgument("APC-A=CD3,PE-A=CD3", c("APC-A", "PE-A"),
                        c("APC-A", "PE-A")),
    "A column is named twice: CD3")
})

test_that("ParseMarkerArgument returns NULL when no argument was given", {
  expect_null(ParseMarkerArgument(NULL, c("CD3"), c("B515-A")))
})

test_that("ParseMarkerArgument rejects an empty argument", {
  expect_error(ParseMarkerArgument(" ", c("CD3"), c("B515-A")),
               "--markers is empty")
})

test_that("SetCytokitSeed makes two runs draw the same subset", {
  SetCytokitSeed()
  first <- sample(100, 5)
  SetCytokitSeed()
  expect_equal(sample(100, 5), first)
})

test_that("SetCytokitSeed takes the seed the caller named", {
  expect_equal(SetCytokitSeed(list(seed = "7")), 7L)
  first <- sample(100, 5)
  SetCytokitSeed(list(seed = "7"))
  expect_equal(sample(100, 5), first)
})

test_that("SetCytokitSeed rejects a seed that is not a number", {
  expect_error(SetCytokitSeed(list(seed = "later")),
               "--seed must be a whole number")
})

test_that("ReportNotes adds up the times a line was reported", {
  notes <- data.frame(note = c("uneven tokens", "uneven tokens", "dropped"),
                      times = c(2L, 3L, 1L), stringsAsFactors = FALSE)
  result <- ReportNotes(notes, NULL)
  expect_equal(nrow(result), 2)
  expect_equal(result$times[result$note == "uneven tokens"], 5L)
})

test_that("ReportNotes puts the most reported line first", {
  notes <- data.frame(note = c("rare", "common"), times = c(1L, 9L),
                      stringsAsFactors = FALSE)
  expect_equal(ReportNotes(notes, NULL)$note, c("common", "rare"))
})

test_that("ReportNotes says how many lines it did not print", {
  notes <- data.frame(note = paste("line", 1:8), times = rep(1L, 8),
                      stringsAsFactors = FALSE)
  printed <- capture.output(ReportNotes(notes, NULL, limit = 3))
  expect_true(any(grepl("and 5 more", printed)))
  expect_true(any(grepl("8 note\\(s\\) over 8 distinct", printed)))
})

test_that("ReportNotes prints nothing when nothing was reported", {
  empty <- data.frame(note = character(0), times = integer(0))
  expect_equal(length(capture.output(ReportNotes(empty, NULL))), 0)
})

test_that("ShortLabel strips the recipe name and the timestamp", {
  # A chained recipe uses the bundle it read as its label. Without this the new
  # bundle carries two recipe names and two timestamps.
  expect_equal(ShortLabel("gate_my_study_2026_08_18_205932"), "my_study")
  expect_equal(ShortLabel("cluster_my_study_2026_08_18_205932"), "my_study")
  expect_equal(ShortLabel("/a/path/annotate_my_study_2026_08_18_205932"),
               "my_study")
})

test_that("ShortLabel leaves a name that carries neither", {
  expect_equal(ShortLabel("my_study"), "my_study")
})

test_that("ShortLabel never returns an empty label", {
  expect_equal(ShortLabel("gate_2026_08_18_205932"), "study")
})

test_that("CloseCytokitBundle checksums the files inside a folder input", {
  # A saved hierarchy is a folder, and md5sum fails on a folder.
  bundle <- withr::local_tempdir()
  folder <- withr::local_tempdir()
  writeLines("a", file.path(folder, "one.txt"))
  writeLines("b", file.path(folder, "two.txt"))
  CloseCytokitBundle(bundle, "gate", list(data = folder), inputs = folder)
  manifest <- jsonlite::fromJSON(file.path(bundle, "manifest.json"))
  expect_equal(sort(names(manifest$inputs)), c("one.txt", "two.txt"))
})

test_that("CloseCytokitBundle writes a manifest with no input at all", {
  bundle <- withr::local_tempdir()
  CloseCytokitBundle(bundle, "inspect", list(data = "x"))
  manifest <- jsonlite::fromJSON(file.path(bundle, "manifest.json"))
  expect_equal(manifest$recipe, "inspect")
  expect_true(file.exists(file.path(bundle, "REPRODUCE.md")))
})

test_that("WriteBundleTable writes a CSV with no row names", {
  bundle <- withr::local_tempdir()
  frame <- data.frame(a = 1:2, b = c("x", "y"), stringsAsFactors = FALSE)
  path <- WriteBundleTable(bundle, frame, "table.csv")
  expect_true(file.exists(path))
  read_back <- utils::read.csv(path, stringsAsFactors = FALSE)
  expect_equal(read_back, frame)
  expect_equal(colnames(read_back)[1], "a")
})

test_that("DescribeFcsFile reads the event count from the header", {
  directory <- withr::local_tempdir()
  path <- WriteTestFcs(directory, n_events = 137)
  described <- DescribeFcsFile(path)
  expect_equal(described$file, "sample.fcs")
  expect_equal(described$events, 137L)
  expect_equal(described$parameters, 4L)
})

test_that("DescribeFcsFile reports no matrix when the file carries none", {
  directory <- withr::local_tempdir()
  described <- DescribeFcsFile(WriteTestFcs(directory))
  expect_equal(described$compensation_state, "no matrix")
})

test_that("DescribeFcsPanel names the kind of every channel", {
  directory <- withr::local_tempdir()
  panel <- DescribeFcsPanel(WriteTestFcs(directory))
  expect_equal(panel$channel, c("FSC-A", "FSC-H", "Ax700-A", "PE-TxRed-A"))
  expect_equal(panel$kind, c("scatter", "scatter", "stain", "stain"))
  expect_equal(panel$marker[3:4], c("CD3", "CD4"))
})

test_that("DescribeFcsIdentity reads the keywords that name the specimen", {
  directory <- withr::local_tempdir()
  path <- WriteTestFcs(directory,
                       keywords = list(`$SRC` = "PBMC_001",
                                       `TUBE NAME` = "13",
                                       `EXPERIMENT NAME` = "Expt 4"))
  identity <- DescribeFcsIdentity(path)
  expect_equal(identity$specimen, "PBMC_001")
  expect_equal(identity$tube, "13")
  expect_equal(identity$experiment, "Expt 4")
})

test_that("DescribeFcsIdentity gives NA for a keyword the file leaves out", {
  directory <- withr::local_tempdir()
  identity <- DescribeFcsIdentity(WriteTestFcs(directory))
  expect_true(is.na(identity$project))
  expect_equal(identity$file, "sample.fcs")
})

test_that("DescribeFcsIdentity treats an empty keyword as absent", {
  # An empty $SRC is not a specimen name, and reading it as one splits a folder
  # into a group whose name is the empty string.
  directory <- withr::local_tempdir()
  path <- WriteTestFcs(directory, keywords = list(`$SRC` = ""))
  expect_true(is.na(DescribeFcsIdentity(path)$specimen))
})
