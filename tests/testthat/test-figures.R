# Tests for R/figures.R.

# The PNG header stores the width and the height as two big endian four byte
# integers at offset 16. Reading them directly avoids a dependency on an image
# package for a test that only needs the size.
PngSize <- function(path) {
  connection <- file(path, "rb")
  on.exit(close(connection))
  header <- readBin(connection, "raw", n = 24)
  Integer32 <- function(bytes) {
    sum(as.integer(bytes) * c(256^3, 256^2, 256, 1))
  }
  c(width = Integer32(header[17:20]), height = Integer32(header[21:24]))
}

test_that("FigureFont returns one family name", {
  family <- FigureFont()
  expect_type(family, "character")
  expect_length(family, 1)
})

test_that("PublicationPalette returns the leading colours unchanged", {
  expect_equal(PublicationPalette(3), c("#0072B2", "#D55E00", "#009E73"))
  expect_equal(PublicationPalette(8)[7:8], c("#000000", "#F0E442"))
  expect_equal(PublicationPalette(1), "#0072B2")
  expect_equal(PublicationPalette(20), kPublicationPalette)
})

test_that("PublicationPalette interpolates past the named colours", {
  colours <- PublicationPalette(28)
  expect_length(colours, 28)
  expect_true(all(grepl("^#[0-9A-Fa-f]{6}$", colours)))
})

test_that("PublicationPalette rejects a size below one", {
  expect_error(PublicationPalette(0), "1 or more")
  expect_error(PublicationPalette(NA), "1 or more")
})

test_that("CountLabels writes a thousands separator and no exponent", {
  expect_equal(CountLabels(c(50000, 100000)), c("50,000", "100,000"))
  expect_equal(CountLabels(250), "250")
})

test_that("PointDensity ranks a dense cluster above an isolated point", {
  x <- c(rep(0, 200), 10)
  y <- c(rep(0, 200), 10)
  density <- PointDensity(x, y, bins = 16)
  expect_length(density, 201)
  expect_gt(density[1], density[201])
  expect_equal(density[201], 1)
})

test_that("PointDensity survives a constant axis and an empty input", {
  expect_equal(PointDensity(rep(1, 5), 1:5, bins = 4), c(3, 3, 4, 3, 2))
  expect_equal(PointDensity(numeric(0), numeric(0)), numeric(0))
})

test_that("PointDensity rejects two vectors of different lengths", {
  expect_error(PointDensity(1:3, 1:4), "same length")
})

test_that("PlotDensityScatter names the column that is missing", {
  frame <- data.frame(a = 1:5, b = 1:5)
  expect_error(PlotDensityScatter(frame, "a", "z"), "no column called 'z'")
})

test_that("PlotDensityScatter draws every event it is given", {
  frame <- data.frame(a = rnorm(50), b = rnorm(50))
  plot <- PlotDensityScatter(frame, "a", "b")
  expect_s3_class(plot, "ggplot")
  expect_equal(nrow(plot$data), 50)
  expect_true("shade" %in% colnames(plot$data))
})

test_that("PlotDensityScatter cuts the shade into the number of steps asked", {
  frame <- withr::with_seed(1, data.frame(a = rnorm(4000), b = rnorm(4000)))
  plot <- PlotDensityScatter(frame, "a", "b", shades = 16)
  expect_lte(length(unique(plot$data$shade)), 16)
  expect_equal(min(plot$data$shade), 0)
  expect_equal(max(plot$data$shade), 15)
})

test_that("PlotDensityScatter gives one shade when every event is stacked", {
  frame <- data.frame(a = rep(1, 20), b = rep(2, 20))
  plot <- PlotDensityScatter(frame, "a", "b")
  expect_equal(unique(plot$data$shade), 0)
})

test_that("ThemePublication drops the minor grid and keeps a panel border", {
  theme <- ThemePublication()
  expect_s3_class(theme$panel.grid.minor, "element_blank")
  expect_s3_class(theme$panel.border, "element_rect")
})

test_that("ThemeEmbedding drops the axis text that ThemePublication keeps", {
  expect_s3_class(ThemeEmbedding()$axis.text, "element_blank")
  expect_s3_class(ThemePublication()$axis.text, "element_text")
})

test_that("SaveFigure writes the pixel count that the dpi asks for", {
  path <- withr::local_tempfile(fileext = ".png")
  plot <- ggplot2::ggplot(data.frame(x = 1:3, y = 1:3),
                          ggplot2::aes(x = .data$x, y = .data$y)) +
    ggplot2::geom_point()
  SaveFigure(plot, path, width = 4, height = 3, dpi = 300)
  expect_true(file.exists(path))
  expect_equal(unname(PngSize(path)), c(1200, 900))
})

test_that("OpenFigureDevice sizes a base graphics plot in inches", {
  path <- withr::local_tempfile(fileext = ".png")
  OpenFigureDevice(path, width = 5, height = 4, res = 200)
  graphics::plot(1:10, 1:10)
  CloseFigureDevice()
  expect_equal(unname(PngSize(path)), c(1000, 800))
})

test_that("the default resolution is the one a journal accepts", {
  expect_equal(kFigureDpi, 600)
  expect_equal(kFigureRes, 600)
})

test_that("AxisLabels separates thousands and leaves a small scale alone", {
  expect_equal(AxisLabels(c(50000, 100000)), c("50,000", "100,000"))
  expect_equal(AxisLabels(c(-2.5, 0, 2.5)), c("-2.5", "0", "2.5"))
  expect_equal(AxisLabels(numeric(0)), character(0))
})

test_that("CountLabels and AxisLabels keep a missing break missing", {
  expect_true(is.na(CountLabels(c(1000, NA))[2]))
  expect_true(is.na(AxisLabels(c(1000, NA))[2]))
})

test_that("SaveFigure writes a vector file when the path ends in svg", {
  path <- withr::local_tempfile(fileext = ".svg")
  plot <- ggplot2::ggplot(data.frame(x = 1:3, y = 1:3),
                          ggplot2::aes(x = .data$x, y = .data$y)) +
    ggplot2::geom_point()
  SaveFigure(plot, path, width = 4, height = 3)
  expect_true(file.exists(path))

  first <- readLines(path, n = 4, warn = FALSE)
  expect_true(any(grepl("<svg", first)))
})

test_that("SaveFigure still writes a raster file for any other path", {
  path <- withr::local_tempfile(fileext = ".png")
  plot <- ggplot2::ggplot(data.frame(x = 1:3, y = 1:3),
                          ggplot2::aes(x = .data$x, y = .data$y)) +
    ggplot2::geom_point()
  SaveFigure(plot, path, width = 4, height = 3, dpi = 200)
  expect_equal(unname(PngSize(path)), c(800, 600))
})

test_that("OpenFigureDevice opens a vector device for an svg path", {
  path <- withr::local_tempfile(fileext = ".svg")
  OpenFigureDevice(path, width = 5, height = 4)
  graphics::plot(1:10, 1:10)
  CloseFigureDevice()
  expect_true(any(grepl("<svg", readLines(path, n = 4, warn = FALSE))))
})

test_that("ColourbarGuide asks for a bar drawn as rectangles", {
  guide <- ColourbarGuide()
  expect_equal(guide$params$display, "rectangles")
})

test_that("a continuous fill writes an svg with no raster fallback", {
  path <- withr::local_tempfile(fileext = ".svg")
  frame <- expand.grid(x = letters[1:4], y = LETTERS[1:4])
  frame$v <- seq_len(nrow(frame))
  plot <- ggplot2::ggplot(frame, ggplot2::aes(.data$x, .data$y,
                                              fill = .data$v)) +
    ggplot2::geom_tile() +
    ggplot2::scale_fill_viridis_c(guide = ColourbarGuide()) +
    ThemePublication()
  SaveFigure(plot, path, width = 5, height = 4)

  # cairo wraps the whole plot in an feImage filter when it cannot express an
  # element, and several renderers draw that filter wrongly.
  written <- readLines(path, warn = FALSE)
  expect_false(any(grepl("<filter", written)))
  expect_false(any(grepl("<image", written)))
})
