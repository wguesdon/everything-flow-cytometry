"""Figure style for the Python arm of the analysis.

The R arm sets its style in ``R/figures.R``. A report that carries figures from
both languages has to look like one report, so the resolution, the type and the
discrete palette are the same values here. Keep the two files in step: a change
to one is a change to both.
"""

from __future__ import annotations

from pathlib import Path

import matplotlib
import matplotlib.pyplot as plt
from matplotlib.figure import Figure

# 600 dpi matches kFigureDpi in R/figures.R. Nothing here is tracked in git, so
# the file size does not constrain the choice.
#
# A figure without a large point cloud does better as a vector. matplotlib picks
# its writer from the file extension, so save_figure needs no argument for it:
# pass a path ending in .svg and the figure stays sharp at any zoom.
FIGURE_DPI = 600

# The same twenty colours as kPublicationPalette in R/figures.R. The first eight
# are the Okabe and Ito set, with yellow held back to eighth because a small
# marker in #F0E442 on a white panel is hard to see.
PUBLICATION_PALETTE = [
    "#0072B2", "#D55E00", "#009E73", "#CC79A7", "#E69F00",
    "#56B4E9", "#000000", "#F0E442", "#332288", "#117733",
    "#882255", "#44AA99", "#999933", "#AA4499", "#661100",
    "#6699CC", "#CC6677", "#88CCEE", "#DDCC77", "#AA7744",
]

# Nimbus Sans is the Helvetica metric clone that the container installs. DejaVu
# Sans is the matplotlib default and is always present, so it closes the list.
FONT_ORDER = ["Nimbus Sans", "TeX Gyre Heros", "Helvetica", "Arial",
              "DejaVu Sans"]


def available_font() -> str:
    """Name the first font of FONT_ORDER that matplotlib can see.

    Returns:
        A font family name. A missing family makes matplotlib warn once for
        every string it draws, so the name is resolved rather than assumed.
    """
    installed = {font.name for font in matplotlib.font_manager.fontManager.ttflist}
    for family in FONT_ORDER:
        if family in installed:
            return family
    return "sans-serif"


def apply_publication_style() -> None:
    """Set the matplotlib defaults that every figure in this repository uses.

    The call changes global state, so make it once at the start of a script and
    not inside a plotting function.
    """
    plt.rcParams.update({
        "figure.dpi": 110,
        "savefig.dpi": FIGURE_DPI,
        "savefig.bbox": "tight",
        "savefig.facecolor": "white",
        "figure.facecolor": "white",
        "axes.facecolor": "white",
        "font.family": available_font(),
        "font.size": 10,
        "axes.titlesize": 11.5,
        "axes.labelsize": 10.5,
        "axes.labelcolor": "black",
        "axes.edgecolor": "#333333",
        "axes.linewidth": 0.8,
        "axes.titlelocation": "left",
        "axes.grid": False,
        "xtick.color": "#333333",
        "ytick.color": "#333333",
        "xtick.labelsize": 9,
        "ytick.labelsize": 9,
        "xtick.major.width": 0.6,
        "ytick.major.width": 0.6,
        "legend.frameon": False,
        "legend.fontsize": 9,
        "legend.title_fontsize": 9.5,
        "axes.prop_cycle": matplotlib.cycler(color=PUBLICATION_PALETTE),
    })


def save_figure(figure: Figure, path: Path, dpi: int = FIGURE_DPI) -> Path:
    """Write a figure, as a vector when the path asks for one.

    matplotlib selects its writer from the file extension, so a path ending in
    .svg produces a vector figure and the dpi is then ignored.

    Args:
        figure: The matplotlib figure to write.
        path: The destination file.
        dpi: The resolution in dots per inch, for a raster path only.

    Returns:
        The path that was written.
    """
    figure.savefig(path, dpi=dpi, bbox_inches="tight", facecolor="white")
    plt.close(figure)
    return path
