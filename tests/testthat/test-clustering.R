# Tests for AnnotateClusters in R/clustering.R.
#
# SelectByHighestMarker, from the same file, is tested in
# test-select-highest.R.

MakeExpression <- function() {
  data.frame(
    cluster = 1:3,
    events = c(100L, 200L, 300L),
    percent_of_total = c(16.667, 33.333, 50),
    A = c(10, 0, 5),
    B = c(5, 0, 10),
    C = c(5, 0, 10),
    stringsAsFactors = FALSE
  )
}

MakeDefinitions <- function() {
  data.frame(
    cell_type = c("short", "long"),
    A = c("pos", "pos"),
    B = c("", "pos"),
    C = c("", "pos"),
    stringsAsFactors = FALSE
  )
}

test_that("AnnotateClusters divides by the total weight and does not add up", {
  annotation <- AnnotateClusters(MakeExpression(), MakeDefinitions())

  # Each column is scaled from 0 to 1 across the clusters, so cluster 1 scores
  # 1 on A and 0.5 on B and C. The short definition names A alone and scores 1.
  # The long definition names three markers and scores (1 + 0.5 + 0.5) / 3.
  expect_equal(annotation$cell_type[1], "short")
  expect_equal(annotation$score[1], 1)
  expect_equal(annotation$margin[1], 1 - 2 / 3)

  # Under a sum the long definition would score 2 against 1 and would take this
  # cluster for naming more markers rather than for fitting it better. The
  # division is what stops that.
  expect_gt(1, 2 / 3)
})

test_that("AnnotateClusters still lets a long definition win on merit", {
  annotation <- AnnotateClusters(MakeExpression(), MakeDefinitions())
  expect_equal(annotation$cell_type[3], "long")
  expect_equal(annotation$score[3], (0.5 + 1 + 1) / 3)
  expect_equal(annotation$runner_up[3], "short")
})

test_that("AnnotateClusters doubles the weight of a high expectation", {
  definitions <- data.frame(
    cell_type = c("plain", "emphatic"),
    A = c("pos", "high"),
    stringsAsFactors = FALSE
  )
  annotation <- AnnotateClusters(MakeExpression(), definitions)

  # 'high' contributes twice the scaled value and carries twice the weight, so
  # the score is the same as 'pos' and the label is decided by the order.
  expect_equal(annotation$score[1], 1)
  expect_equal(annotation$margin[1], 0)
})

test_that("AnnotateClusters carries the cluster identity through unchanged", {
  annotation <- AnnotateClusters(MakeExpression(), MakeDefinitions())
  expect_equal(annotation$cluster, 1:3)
  expect_equal(annotation$events, c(100L, 200L, 300L))
  expect_equal(nrow(annotation), 3)
})

test_that("AnnotateClusters names the markers when none of them match", {
  definitions <- data.frame(
    cell_type = "unrelated", CD99 = "pos", stringsAsFactors = FALSE
  )
  expect_error(AnnotateClusters(MakeExpression(), definitions),
               "share no marker column")
  expect_error(AnnotateClusters(MakeExpression(), definitions), "CD99")
})

test_that("AnnotateClusters gives every cluster the same score on a flat marker", {
  expression <- MakeExpression()
  expression$A <- c(4, 4, 4)
  definitions <- data.frame(cell_type = c("only", "other"), A = c("pos", "neg"),
                            stringsAsFactors = FALSE)
  annotation <- AnnotateClusters(expression, definitions)

  # A column with no spread scales to 0.5 everywhere, so 'pos' and 'neg' tie.
  expect_equal(unique(annotation$score), 0.5)
  expect_equal(unique(annotation$margin), 0)
})

test_that("SubsampleEvents records how many events it drew from", {
  events <- matrix(rnorm(2000), ncol = 2,
                   dimnames = list(NULL, c("CD3", "CD4")))
  drawn <- SubsampleEvents(events, n = 100)
  expect_equal(nrow(drawn), 100)
  expect_equal(attr(drawn, "sampled_from"), 1000)
})

test_that("SubsampleEvents returns every event when there are few enough", {
  events <- matrix(rnorm(40), ncol = 2,
                   dimnames = list(NULL, c("CD3", "CD4")))
  drawn <- SubsampleEvents(events, n = 100)
  expect_equal(nrow(drawn), 20)
  expect_equal(attr(drawn, "sampled_from"), 20)
})

test_that("SubsampleEvents draws the same rows twice with one seed", {
  events <- matrix(rnorm(2000), ncol = 2,
                   dimnames = list(NULL, c("CD3", "CD4")))
  expect_equal(SubsampleEvents(events, n = 50, seed = 7),
               SubsampleEvents(events, n = 50, seed = 7))
})

test_that("ClusterMedianExpression reports the median of every cluster", {
  events <- cbind(CD3 = c(1, 1, 1, 9, 9, 9), CD4 = c(2, 2, 2, 8, 8, 8))
  medians <- ClusterMedianExpression(events, c(1, 1, 1, 2, 2, 2),
                                     c("CD3", "CD4"))
  expect_equal(medians$cluster, c(1, 2))
  expect_equal(medians$CD3, c(1, 9))
  expect_equal(medians$events, c(3L, 3L))
  expect_equal(medians$percent_of_total, c(50, 50))
})

test_that("ClusterMedianExpression rejects a cluster vector of the wrong length", {
  events <- cbind(CD3 = 1:4)
  expect_error(ClusterMedianExpression(events, c(1, 2), "CD3"),
               "They must match")
})

test_that("SummariseCellTypes adds up the clusters that share a label", {
  annotation <- data.frame(
    cluster = 1:3, events = c(10, 20, 30),
    percent_of_total = c(10, 20, 30),
    cell_type = c("T cells", "T cells", "B cells"),
    stringsAsFactors = FALSE)
  summary_table <- SummariseCellTypes(annotation)
  # The row order is not part of the contract when two labels tie, so the
  # values are read by name.
  clusters <- stats::setNames(summary_table$clusters, summary_table$cell_type)
  events <- stats::setNames(summary_table$events, summary_table$cell_type)
  expect_equal(sort(summary_table$cell_type), c("B cells", "T cells"))
  expect_equal(unname(clusters[["T cells"]]), 2L)
  expect_equal(unname(clusters[["B cells"]]), 1L)
  expect_equal(unname(events[["T cells"]]), 30)
})

test_that("SummariseCellTypes sorts the largest cell type first", {
  annotation <- data.frame(
    cluster = 1:2, events = c(10, 90), percent_of_total = c(10, 90),
    cell_type = c("rare", "common"), stringsAsFactors = FALSE)
  expect_equal(SummariseCellTypes(annotation)$cell_type, c("common", "rare"))
})

test_that("ReadCellTypeDefinitions accepts pos, neg, high and empty", {
  path <- withr::local_tempfile(fileext = ".csv")
  writeLines(c("cell_type,CD3,CD4,note", "T cells,pos,high,",
               "B cells,neg,,a note"), path)
  definitions <- ReadCellTypeDefinitions(path)
  expect_equal(definitions$cell_type, c("T cells", "B cells"))
  expect_equal(definitions$CD3, c("pos", "neg"))
})

test_that("ReadCellTypeDefinitions rejects a value it cannot read", {
  # A typo such as "positive" would otherwise score as no expectation at all.
  path <- withr::local_tempfile(fileext = ".csv")
  writeLines(c("cell_type,CD3,note", "T cells,positive,"), path)
  expect_error(ReadCellTypeDefinitions(path),
               "must be pos, neg, high or empty")
})

test_that("ReadCellTypeDefinitions rejects a file with no cell_type column", {
  path <- withr::local_tempfile(fileext = ".csv")
  writeLines(c("name,CD3", "T cells,pos"), path)
  expect_error(ReadCellTypeDefinitions(path), "no 'cell_type' column")
})

test_that("ReadCellTypeDefinitions rejects a path that does not exist", {
  expect_error(ReadCellTypeDefinitions(file.path(tempdir(), "absent.csv")),
               "does not exist")
})

test_that("RenameChannelsToMarkers puts the marker name on the column", {
  frame <- MakeTestFlowFrame()
  events <- flowCore::exprs(frame)
  renamed <- RenameChannelsToMarkers(events, frame)
  expect_equal(colnames(renamed), c("FSC-A", "FSC-H", "CD3", "CD4"))
})

test_that("RenameChannelsToMarkers keeps a detector whose marker is a placeholder", {
  # OMIP-039 labels two unused detectors "Available", and that is not a marker.
  frame <- MakeTestFlowFrame()
  parameters <- flowCore::parameters(frame)
  flowCore::pData(parameters)$desc <- c(NA, NA, "Available", "CD4")
  flowCore::parameters(frame) <- parameters
  renamed <- RenameChannelsToMarkers(flowCore::exprs(frame), frame)
  expect_equal(colnames(renamed), c("FSC-A", "FSC-H", "Ax700-A", "CD4"))
})

test_that("ExtractGatedEvents returns the events inside one gate", {
  gating_set <- MakeGatedSet()
  events <- ExtractGatedEvents(gating_set, "nonDebris")
  expect_true(is.matrix(events))
  expect_lt(nrow(events), 400)
  expect_gt(nrow(events), 0)
})

test_that("ExtractGatedEvents rejects a population that is not there", {
  gating_set <- MakeGatedSet()
  expect_error(ExtractGatedEvents(gating_set, "Monocytes"),
               "is not in the GatingSet")
})

test_that("RunFlowSomClustering gives one metacluster label per event", {
  events <- withr::with_seed(3, {
    a <- cbind(CD3 = rnorm(400, 1), CD4 = rnorm(400, 1))
    b <- cbind(CD3 = rnorm(400, 9), CD4 = rnorm(400, 9))
    rbind(a, b)
  })
  result <- RunFlowSomClustering(events, c("CD3", "CD4"), grid_size = 4,
                                 n_metaclusters = 2, seed = 5)
  expect_equal(length(result$metacluster), nrow(events))
  expect_equal(length(unique(result$metacluster)), 2)
})

test_that("RunFlowSomClustering rejects a channel the data does not carry", {
  events <- cbind(CD3 = rnorm(50), CD4 = rnorm(50))
  expect_error(RunFlowSomClustering(events, c("CD3", "CD8")),
               "not in the data")
})

test_that("RunFlowSomClustering rejects fewer than two metaclusters", {
  events <- cbind(CD3 = rnorm(50), CD4 = rnorm(50))
  expect_error(RunFlowSomClustering(events, "CD3", n_metaclusters = 1),
               "must be 2 or more")
})

test_that("RunUmapEmbedding returns two coordinates per event", {
  events <- withr::with_seed(4, cbind(CD3 = rnorm(300), CD4 = rnorm(300),
                                      CD8 = rnorm(300)))
  embedding <- RunUmapEmbedding(events, c("CD3", "CD4", "CD8"),
                                n_neighbors = 10, seed = 6)
  expect_equal(nrow(embedding), 300)
  expect_equal(colnames(embedding), c("umap_1", "umap_2"))
})

test_that("RunUmapEmbedding rejects a channel the data does not carry", {
  events <- cbind(CD3 = rnorm(50), CD4 = rnorm(50))
  expect_error(RunUmapEmbedding(events, c("CD3", "CD8")), "not in the data")
})

test_that("RunFlowSomClustering works at two metaclusters", {
  # Its own guard allows two, and FlowSOM's wrapper fails at exactly two,
  # because ConsensusClusterPlus drops a one row matrix to a vector.
  events <- withr::with_seed(31, rbind(
    cbind(CD3 = rnorm(300, 1), CD4 = rnorm(300, 1)),
    cbind(CD3 = rnorm(300, 9), CD4 = rnorm(300, 9))))
  result <- RunFlowSomClustering(events, c("CD3", "CD4"), grid_size = 4,
                                 n_metaclusters = 2, seed = 8)
  expect_equal(length(unique(result$metacluster)), 2)
  expect_equal(length(result$metacluster), nrow(events))
})
