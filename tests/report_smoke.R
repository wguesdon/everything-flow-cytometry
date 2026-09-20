source("R/cytokit.R")

bundle <- file.path(tempdir(), "cytokit_report_smoke")
dir.create(bundle, recursive = TRUE, showWarnings = FALSE)
utils::write.csv(
  data.frame(population = "T cells", count = 12L),
  file.path(bundle, "counts.csv"),
  row.names = FALSE
)
Sys.setenv(
  CYTOKIT_GIT_COMMIT = "smoke-test-commit",
  CYTOKIT_GIT_DIRTY = "false",
  CYTOKIT_IMAGE = "cytokit:test",
  CYTOKIT_IMAGE_DIGEST = "sha256:test"
)

CloseCytokitBundle(
  bundle,
  "inspect",
  list(seed = "42"),
  command = "cytokit inspect --data study"
)

expected <- c(
  "manifest.json",
  "report.qmd",
  "report.html",
  "report_provenance.json"
)
stopifnot(all(file.exists(file.path(bundle, expected))))
manifest <- jsonlite::read_json(file.path(bundle, "manifest.json"), simplifyVector = TRUE)
stopifnot(identical(manifest$git_commit, "smoke-test-commit"))
stopifnot("report.html" %in% manifest$outputs$file)
