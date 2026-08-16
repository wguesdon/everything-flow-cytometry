# Entry point for the test suite.
#
# Run it in the container:
#   podman run --rm -v "$PWD:/work:z" -w /work everything-flow-cytometry:latest \
#     Rscript tests/testthat.R

library(testthat)
library(flowCore)
library(flowWorkspace)
library(openCyto)
library(withr)

# The functions live in R/ and are sourced rather than installed, because this is
# an example repository and not a package. A helper file in tests/testthat/ builds
# the fixtures.
for (file in list.files("R", pattern = "\\.R$", full.names = TRUE)) {
  source(file)
}

results <- test_dir("tests/testthat", reporter = "summary", stop_on_failure = TRUE)
