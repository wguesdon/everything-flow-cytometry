# Figure style shared by every script in this repository.
#
# A report figure has to hold up at the size a journal prints it, so the
# resolution, the type and the colour are set in one place and not in each
# script. Three faults motivated this file. A ggplot saved at 150 dpi is
# upscaled by a browser on a high density display. A base graphics device opened
# without `res` renders its text at 72 dpi against a canvas measured in pixels,
# which makes the type small and soft. The default discrete palette of ggplot2
# walks the colour wheel at one lightness, so a reader with a colour vision
# deficiency cannot separate the classes.
#
# tests/testthat.R sources every file in this folder into one environment, and
# the second definition of a name silently replaces the first. Check that a name
# is free before you add it.

# 600 dpi is what a journal asks for line art, and 300 is the floor it accepts
# for a photograph. Nothing here is tracked in git, so the file size does not
# constrain the choice and the higher number applies everywhere.
#
# A figure without a large point cloud does better as a vector than at any
# resolution, because it stays sharp at any zoom. Give SaveFigure a path ending
# in .svg and it writes one. A scatter of tens of thousands of events is the
# exception: cairo writes about 400 bytes for each point, so a UMAP of 50,000
# events is a 20 MB file that a browser struggles to draw. Those stay raster.
kFigureDpi <- 600

# The base graphics devices take a resolution rather than a dpi argument. The
# same number applies, because both mean pixels per inch.
kFigureRes <- 600

# The first eight are the Okabe and Ito (2008) set, which separates under the
# three common forms of colour vision deficiency. Yellow sits eighth rather than
# where Okabe and Ito put it, because a 0.35 point marker in #F0E442 on a white
# panel falls below the contrast a reader needs, and a scale of seven classes
# should not have to use it. The set is complete at eight either way, so the
# guarantee holds. Positions nine to twenty come from the Paul Tol qualitative
# sets and cover a metacluster plot. Past twenty no discrete palette stays
# readable and the scale interpolates instead.
kPublicationPalette <- c(
  "#0072B2", "#D55E00", "#009E73", "#CC79A7", "#E69F00",
  "#56B4E9", "#000000", "#F0E442", "#332288", "#117733",
  "#882255", "#44AA99", "#999933", "#AA4499", "#661100",
  "#6699CC", "#CC6677", "#88CCEE", "#DDCC77", "#AA7744"
)

# The container installs the URW and TeX Gyre families. Nimbus Sans is the
# Helvetica metric clone, which is the face that Nature and Cell set their
# figures in. The fallback matters because a contributor may run a script on a
# host that has neither.
kFigureFontOrder <- c("Nimbus Sans", "TeX Gyre Heros", "Helvetica", "Arial")

#' Name the sans serif face to set a figure in
#'
#' The face is resolved against the fonts that the system reports rather than
#' named as a constant, because a missing family makes a device warn on every
#' string that it draws.
#'
#' @return A font family name. `""` asks the device for its own default.
#' @examples
#' FigureFont()
#' @export
FigureFont <- function() {
  installed <- tryCatch(systemfonts::system_fonts()$family,
                        error = function(e) character(0))
  found <- kFigureFontOrder[kFigureFontOrder %in% installed]
  if (length(found) == 0) {
    return("")
  }
  found[1]
}

#' The theme that every figure in this repository uses
#'
#' The panel keeps a border and drops the minor grid, so a reader reads a value
#' off the axis rather than off a background rule. The text is darker than the
#' ggplot2 default, because grey50 on white falls below the contrast that a
#' printed figure needs.
#'
#' @param base_size The base font size in points.
#' @param base_family The font family. Defaults to the resolved sans serif face.
#' @return A `ggplot2` theme.
#' @examples
#' ThemePublication()
#' @export
ThemePublication <- function(base_size = 11, base_family = FigureFont()) {
  ggplot2::theme_bw(base_size = base_size, base_family = base_family) +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      panel.grid.major = ggplot2::element_line(colour = "grey92",
                                               linewidth = 0.3),
      panel.border = ggplot2::element_rect(colour = "grey20", fill = NA,
                                           linewidth = 0.5),
      axis.text = ggplot2::element_text(colour = "grey20",
                                        size = ggplot2::rel(0.85)),
      axis.title = ggplot2::element_text(colour = "black"),
      axis.ticks = ggplot2::element_line(colour = "grey20", linewidth = 0.3),
      strip.background = ggplot2::element_rect(fill = "grey94",
                                               colour = "grey20",
                                               linewidth = 0.5),
      strip.text = ggplot2::element_text(
        colour = "black", size = ggplot2::rel(0.9),
        margin = ggplot2::margin(3, 3, 3, 3)
      ),
      plot.title = ggplot2::element_text(face = "plain",
                                         size = ggplot2::rel(1.05), hjust = 0,
                                         margin = ggplot2::margin(b = 6)),
      plot.subtitle = ggplot2::element_text(colour = "grey30",
                                            size = ggplot2::rel(0.9),
                                            hjust = 0),
      plot.caption = ggplot2::element_text(colour = "grey30",
                                           size = ggplot2::rel(0.8), hjust = 0),
      plot.tag = ggplot2::element_text(face = "bold",
                                       size = ggplot2::rel(1.2)),
      legend.key = ggplot2::element_blank(),
      legend.background = ggplot2::element_blank(),
      legend.title = ggplot2::element_text(size = ggplot2::rel(0.9)),
      legend.text = ggplot2::element_text(size = ggplot2::rel(0.85)),
      plot.margin = ggplot2::margin(6, 8, 6, 6)
    )
}

#' The theme for a UMAP or another embedding
#'
#' The two axes of an embedding carry no unit, so a tick and a grid line invite
#' a reader to measure a distance that has no meaning. This drops both and keeps
#' the axis titles.
#'
#' @param base_size The base font size in points.
#' @param base_family The font family.
#' @return A `ggplot2` theme.
#' @examples
#' ThemeEmbedding()
#' @export
ThemeEmbedding <- function(base_size = 11, base_family = FigureFont()) {
  ThemePublication(base_size = base_size, base_family = base_family) +
    ggplot2::theme(
      panel.grid.major = ggplot2::element_blank(),
      axis.text = ggplot2::element_blank(),
      axis.ticks = ggplot2::element_blank()
    )
}

#' Return the discrete colours for a number of classes
#'
#' @param n How many colours the scale needs.
#' @return A character vector of `n` hex colours.
#' @examples
#' PublicationPalette(3)
#' @export
PublicationPalette <- function(n) {
  if (!is.numeric(n) || length(n) != 1 || is.na(n) || n < 1) {
    stop("n must be a single number of 1 or more.")
  }
  n <- as.integer(n)
  if (n <= length(kPublicationPalette)) {
    return(kPublicationPalette[seq_len(n)])
  }
  grDevices::colorRampPalette(kPublicationPalette)(n)
}

#' The discrete colour scale for a figure
#'
#' @param ... Passed to `ggplot2::discrete_scale`, for example `name`.
#' @return A `ggplot2` scale.
#' @examples
#' ScaleColourPublication()
#' @export
ScaleColourPublication <- function(...) {
  ggplot2::discrete_scale("colour", palette = PublicationPalette, ...)
}

#' The discrete fill scale for a figure
#'
#' @param ... Passed to `ggplot2::discrete_scale`, for example `name`.
#' @return A `ggplot2` scale.
#' @examples
#' ScaleFillPublication()
#' @export
ScaleFillPublication <- function(...) {
  ggplot2::discrete_scale("fill", palette = PublicationPalette, ...)
}

#' The colour bar for a continuous scale
#'
#' ggplot2 draws a continuous colour bar as an embedded raster image by default.
#' The cairo SVG device cannot express that inside a vector file, so it wraps the
#' whole plot in an `feImage` filter and falls back to rasterising it. Several
#' renderers draw that filter wrongly, and one of them drew a heat map as a
#' black rectangle. `display = "rectangles"` draws the bar as a stack of
#' rectangles instead, which avoids the fallback and stays sharp at any zoom.
#'
#' @param ... Passed to `ggplot2::guide_colourbar`.
#' @return A `guide` to give to the `guide` argument of a continuous scale.
#' @examples
#' ColourbarGuide()
#' @export
ColourbarGuide <- function(...) {
  ggplot2::guide_colourbar(display = "rectangles", ...)
}

#' Enlarge the key of a colour legend
#'
#' A scatter drawn at `size = 0.2` puts a key in the legend that a reader cannot
#' see. The point in the legend is set apart from the point in the panel.
#'
#' @param size The point size in the legend.
#' @return A `guides` object to add to a plot.
#' @examples
#' LegendPoints()
#' @export
LegendPoints <- function(size = 2.6) {
  ggplot2::guides(colour = ggplot2::guide_legend(
    override.aes = list(size = size, alpha = 1)
  ))
}

#' Format an axis of counts with a thousands separator
#'
#' The ggplot2 default prints a forward scatter axis as `1e+05`, which is a
#' notation for a log scale and not for a linear one.
#'
#' @param x A numeric vector, supplied by the scale.
#' @return A character vector of labels.
#' @examples
#' CountLabels(c(50000, 100000))
#' @export
CountLabels <- function(x) {
  scales::label_number(big.mark = ",", accuracy = 1)(x)
}

#' Format an axis that may hold counts or may hold a transformed intensity
#'
#' A scatter channel runs to a quarter of a million and reads best with a
#' thousands separator. An arcsinh transformed channel runs from about -3 to
#' about 6 and reads worse with one. The break values decide which form applies,
#' because a single figure carries both kinds of axis.
#'
#' @param x A numeric vector, supplied by the scale.
#' @return A character vector of labels.
#' @examples
#' AxisLabels(c(50000, 100000))
#' AxisLabels(c(-2.5, 0, 2.5))
#' @export
AxisLabels <- function(x) {
  finite <- x[is.finite(x)]
  if (length(finite) > 0 && max(abs(finite)) >= 1000) {
    return(CountLabels(x))
  }
  scales::label_number(drop0trailing = TRUE)(x)
}

#' Estimate the local density of every point in a scatter
#'
#' A flow cytometry scatter of a million events saturates to a black rectangle,
#' and the structure that a reader needs is the density inside it. The estimate
#' is a two dimensional histogram smoothed with a box kernel, which needs no
#' package outside base R and costs one pass over the data.
#'
#' @param x A numeric vector.
#' @param y A numeric vector of the same length as `x`.
#' @param bins The number of bins on each axis.
#' @return A numeric vector of smoothed counts, one for each point.
#' @examples
#' PointDensity(rnorm(100), rnorm(100), bins = 8)
#' @export
PointDensity <- function(x, y, bins = 128) {
  if (length(x) != length(y)) {
    stop("x and y must have the same length.")
  }
  if (length(x) == 0) {
    return(numeric(0))
  }
  BinIndex <- function(values) {
    span <- range(values, finite = TRUE)
    if (!all(is.finite(span)) || span[1] == span[2]) {
      return(rep(1L, length(values)))
    }
    edges <- seq(span[1], span[2], length.out = bins + 1)
    index <- .bincode(values, edges, include.lowest = TRUE)
    index[is.na(index)] <- 1L
    index
  }
  x_bin <- BinIndex(x)
  y_bin <- BinIndex(y)
  counts <- matrix(0, nrow = bins, ncol = bins)
  tallied <- tabulate((y_bin - 1L) * bins + x_bin, nbins = bins * bins)
  counts[] <- tallied

  # A raw histogram at this many bins is noisy, so a three by three box kernel
  # smooths it. A shift that leaves the grid contributes zero rather than the
  # repeated edge cell, because repeating the edge counts a corner four times
  # and reports a lone event out there as a small cluster.
  smoothed <- matrix(0, nrow = bins, ncol = bins)
  for (row_shift in -1:1) {
    for (column_shift in -1:1) {
      rows_from <- seq_len(bins) + row_shift
      columns_from <- seq_len(bins) + column_shift
      rows_kept <- rows_from >= 1 & rows_from <= bins
      columns_kept <- columns_from >= 1 & columns_from <= bins
      smoothed[rows_kept, columns_kept] <-
        smoothed[rows_kept, columns_kept] +
        counts[rows_from[rows_kept], columns_from[columns_kept], drop = FALSE]
    }
  }
  smoothed[cbind(x_bin, y_bin)]
}

#' Draw a two dimensional scatter shaded by local density
#'
#' @param frame A `data.frame` with the columns named by `x` and `y`.
#' @param x The column on the horizontal axis.
#' @param y The column on the vertical axis.
#' @param x_label The label for the horizontal axis.
#' @param y_label The label for the vertical axis.
#' @param title The panel title.
#' @param point_size The point size.
#' @param bins The number of density bins on each axis.
#' @param shades How many steps the shading uses. The scale reads as continuous
#'   well below the number of events, so the steps cost nothing and they keep the
#'   colour count of the PNG in a range that compresses.
#' @return A `ggplot` object.
#' @examples
#' \dontrun{
#' PlotDensityScatter(frame, "FSC-A", "SSC-A")
#' }
#' @export
PlotDensityScatter <- function(frame, x, y, x_label = x, y_label = y,
                               title = NULL, point_size = 0.25, bins = 128,
                               shades = 192) {
  for (column in c(x, y)) {
    if (!column %in% colnames(frame)) {
      stop("The frame has no column called '", column, "'.")
    }
  }
  drawn <- data.frame(x = frame[[x]], y = frame[[y]])
  density <- PointDensity(drawn$x, drawn$y, bins = bins)

  # A count of events per cell is heavily skewed, so the square root is taken
  # before the scale is cut into steps. Without it almost every event lands in
  # the first step and the panel reads as one colour.
  shade <- sqrt(density)
  span <- range(shade)
  drawn$shade <- if (span[2] > span[1]) {
    round((shade - span[1]) / (span[2] - span[1]) * (shades - 1))
  } else {
    rep(0, length(shade))
  }

  # The sparse events are drawn first so that a rare population is not buried
  # under the dense core it sits beside.
  drawn <- drawn[order(drawn$shade), , drop = FALSE]

  ggplot2::ggplot(drawn, ggplot2::aes(x = .data$x, y = .data$y,
                                      colour = .data$shade)) +
    ggplot2::geom_point(size = point_size, shape = 16) +
    ggplot2::scale_colour_viridis_c(option = "viridis", guide = "none") +
    ggplot2::labs(title = title, x = x_label, y = y_label) +
    ThemePublication()
}

#' Save a figure, as a vector when the path asks for one
#'
#' A path ending in `.svg` writes a vector figure through the cairo device,
#' which stays sharp at any zoom. Any other path writes a raster figure through
#' `ragg`, which renders text through the system font stack and gives a cleaner
#' glyph than the default `png` device on Linux.
#'
#' Choose the vector form for a figure whose elements can be counted: a tree, a
#' bar chart, a box plot, a heat map. Choose the raster form for a scatter of
#' more than a few thousand events, because a vector records every one of them.
#'
#' @param plot A `ggplot` object.
#' @param path The file to write. The extension selects the device.
#' @param width The width in inches.
#' @param height The height in inches.
#' @param dpi The resolution of a raster figure. Defaults to `kFigureDpi`.
#' @return The path, invisibly.
#' @examples
#' \dontrun{
#' SaveFigure(plot, "output/scatter.png", width = 9, height = 6)
#' SaveFigure(plot, "output/tree.svg", width = 9, height = 6)
#' }
#' @export
SaveFigure <- function(plot, path, width, height, dpi = kFigureDpi) {
  if (identical(tolower(tools::file_ext(path)), "svg")) {
    grDevices::svg(path, width = width, height = height, bg = "white")
    on.exit(grDevices::dev.off(), add = TRUE)
    print(plot)
    return(invisible(path))
  }
  device <- if (requireNamespace("ragg", quietly = TRUE)) {
    ragg::agg_png
  } else {
    grDevices::png
  }
  ggplot2::ggsave(filename = path, plot = plot, width = width,
                  height = height, units = "in", dpi = dpi,
                  bg = "white", device = device)
  invisible(path)
}

#' Open a device for a plot that base graphics draws
#'
#' A `GatingSet` draws through base graphics and returns no object, so it needs
#' a device rather than a save call. A path ending in `.svg` opens a vector
#' device, which is what a gate tree wants: it is a few dozen labels and a few
#' dozen edges, and its text has to stay readable when a reader zooms in.
#'
#' The size is given in inches so that the type scales with the figure. A size
#' in pixels fixes the canvas and leaves the type at one size, which is the
#' fault that makes a tree unreadable.
#'
#' @param path The file to write. The extension selects the device.
#' @param width The width in inches.
#' @param height The height in inches.
#' @param res The resolution of a raster device. Defaults to `kFigureRes`.
#' @param pointsize The base type size in points.
#' @return The path, invisibly.
#' @examples
#' \dontrun{
#' OpenFigureDevice("output/tree.svg", width = 9, height = 6)
#' plot(gating_set)
#' CloseFigureDevice()
#' }
#' @export
OpenFigureDevice <- function(path, width, height, res = kFigureRes,
                             pointsize = 12) {
  if (identical(tolower(tools::file_ext(path)), "svg")) {
    grDevices::svg(path, width = width, height = height, bg = "white",
                   pointsize = pointsize)
    return(invisible(path))
  }
  if (requireNamespace("ragg", quietly = TRUE)) {
    ragg::agg_png(path, width = width, height = height, units = "in",
                  res = res, background = "white", pointsize = pointsize)
  } else {
    grDevices::png(path, width = width, height = height, units = "in",
                   res = res, bg = "white", pointsize = pointsize,
                   type = "cairo")
  }
  invisible(path)
}

#' Close the device that `OpenFigureDevice` opened
#'
#' @return `NULL`, invisibly.
#' @examples
#' \dontrun{
#' CloseFigureDevice()
#' }
#' @export
CloseFigureDevice <- function() {
  grDevices::dev.off()
  invisible(NULL)
}
