"""Read the gated OMIP-058 populations and find their subsets by clustering.

The R half of this workflow, `scripts/13_omip58_prepare.R`, computes the
spillover matrix from the deposited single stains, applies it, removes the
debris and the dead cells, splits the T cells from the rest, and writes each
population to disk as FCS. This module reads those files with pytometry and
looks for the subsets that a one dimensional cut could not resolve.

The functions here hold the analysis. `scripts/14_omip58_pytometry.py` runs it.
"""

from __future__ import annotations

import re
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import pytometry as pm
import scanpy as sc

# The token that identifies each marker of the panel, matching kOmip58Tokens in
# R/omip58.R. A token match beats a substring match, because "cd16" is a
# substring of "cd161" and a substring rule resolves CD16 to the CD161 detector.
MARKER_TOKENS: dict[str, str] = {
    "viability": "live",
    "CD3": "cd3",
    "CD4": "cd4",
    "CD8": "cd8",
    "CD16": "cd16",
    "CD56": "cd56",
    "CD161": "cd161",
    "HLADR": "hla",
    "CCR7": "ccr7",
    "CD45RA": "cd45ra",
    "CD27": "cd27",
    "CD28": "cd28",
    "CD95": "cd95",
    "tetramer": "tetramer",
    "Va72": "va7",
    "Vd1": "vd1",
    "Vd2": "vd2",
    "Vg9": "vg9",
}

SCATTER_PATTERN = re.compile(r"^(FSC|SSC|Time)", re.IGNORECASE)


def marker_tokens(name: str) -> list[str]:
    """Split a marker name into lower case tokens.

    Args:
        name: A marker name as the FCS file writes it, such as `TCR Va7_2 BV711`.

    Returns:
        The lower case tokens of the name, with every run of non alphanumeric
        characters treated as a separator.
    """
    return [token for token in re.split(r"[^a-z0-9]+", name.lower()) if token]


def short_marker_names(var: pd.DataFrame) -> dict[str, str]:
    """Map each panel marker to the short name the definitions table uses.

    Args:
        var: The `.var` table of an AnnData object read from an FCS file. It
            must carry a `marker` column.

    Returns:
        A dictionary from the marker name in the file to the short name.

    Raises:
        ValueError: If a token matches no marker or more than one.
    """
    resolved: dict[str, str] = {}
    for short, token in MARKER_TOKENS.items():
        hits = [name for name in var.index if token in marker_tokens(str(name))]
        if len(hits) != 1:
            raise ValueError(
                f"The token '{token}', which finds {short}, matches {len(hits)} markers: "
                f"{hits}. It must match one."
            )
        resolved[hits[0]] = short
    return resolved


def read_handoff(
    handoff_dir: Path,
    cofactor: float = 150.0,
    populations: list[str] | None = None,
) -> ad.AnnData:
    """Read every population the R half wrote and stack them into one object.

    The values on disk are compensated and untransformed, so the transform is
    applied here. The cofactor matches the one the R half used for its own cuts,
    which is what lets a threshold mean the same thing on both sides.

    Args:
        handoff_dir: The folder holding the FCS files and `manifest.csv`.
        cofactor: The arcsinh cofactor.
        populations: Read only these populations. `None` reads every one. The
            populations nest, so reading two of them counts an event twice.

    Returns:
        One AnnData object carrying every file, with `donor`, `population` and
        `source_file` in `.obs` and the short marker name in `.var["short"]`.

    Raises:
        FileNotFoundError: If the manifest or a file it names is missing.
    """
    manifest_path = handoff_dir / "manifest.csv"
    if not manifest_path.exists():
        raise FileNotFoundError(
            f"{manifest_path} is missing. Run scripts/13_omip58_prepare.R first."
        )
    manifest = pd.read_csv(manifest_path)
    if populations is not None:
        manifest = manifest[manifest["population"].isin(populations)]
        if manifest.empty:
            raise FileNotFoundError(
                f"The manifest names no population in {populations}."
            )

    pieces = []
    for row in manifest.itertuples():
        path = handoff_dir / row.file
        if not path.exists():
            raise FileNotFoundError(f"The manifest names a file that is missing: {path}")
        adata = pm.io.read_fcs(str(path))
        adata.obs["donor"] = str(row.donor)
        adata.obs["population"] = row.population
        adata.obs["source_file"] = row.source_file
        pieces.append(adata)

    merged = ad.concat(pieces, index_unique="-")
    merged.var = pieces[0].var.copy()

    keep = [name for name in merged.var_names if not SCATTER_PATTERN.match(str(name))]
    merged = merged[:, keep].copy()

    short = short_marker_names(merged.var)
    merged.var["short"] = [short.get(name, str(name)) for name in merged.var_names]

    pm.tl.normalize_arcsinh(merged, cofactor=cofactor)
    return merged


def subsample(adata: ad.AnnData, per_donor: int, seed: int = 42) -> ad.AnnData:
    """Take a fixed number of events from each donor.

    A neighbour graph over a million events costs more memory than it returns in
    resolution. Sampling per donor rather than over the pooled object keeps the
    two donors equally represented whatever their event counts.

    Args:
        adata: The object to sample.
        per_donor: How many events to keep for each donor.
        seed: The seed for the generator.

    Returns:
        A copy holding the sampled events.
    """
    generator = np.random.default_rng(seed)
    chosen: list[np.ndarray] = []
    positions = np.arange(adata.n_obs)
    for donor in sorted(adata.obs["donor"].unique()):
        available = positions[(adata.obs["donor"] == donor).to_numpy()]
        take = min(per_donor, available.size)
        chosen.append(generator.choice(available, size=take, replace=False))
    return adata[np.sort(np.concatenate(chosen))].copy()


def cluster(
    adata: ad.AnnData,
    resolution: float = 1.0,
    n_neighbors: int = 30,
    seed: int = 42,
) -> ad.AnnData:
    """Build the neighbour graph, find the communities and embed them.

    Args:
        adata: The transformed object.
        resolution: The leiden resolution. A higher value gives more clusters.
        n_neighbors: The size of the local neighbourhood.
        seed: The seed for the graph and the embedding.

    Returns:
        The same object, with the transformed values kept in the `arcsinh`
        layer, `leiden` in `.obs` and `X_umap` in `.obsm`. `.X` holds the scaled
        values that the graph was built from, so any threshold a caller wants to
        read on the transformed scale must come from the layer.
    """
    adata.layers["arcsinh"] = adata.X.copy()
    sc.pp.scale(adata, max_value=10)
    sc.pp.pca(adata, n_comps=min(20, adata.n_vars - 1), random_state=seed)
    sc.pp.neighbors(adata, n_neighbors=n_neighbors, random_state=seed)
    sc.tl.leiden(adata, resolution=resolution, key_added="leiden",
                 flavor="igraph", n_iterations=2, directed=False,
                 random_state=seed)
    sc.tl.umap(adata, random_state=seed)
    return adata


def cluster_medians(
    adata: ad.AnnData, key: str = "leiden", layer: str | None = "arcsinh"
) -> pd.DataFrame:
    """Return the median expression of every marker in every cluster.

    Args:
        adata: A clustered object.
        key: The `.obs` column holding the cluster label.
        layer: The layer to read. `arcsinh` holds the transformed values, which
            are the ones a reader can interpret. `None` reads `.X`, which after
            clustering holds scaled values.

    Returns:
        A table of clusters by short marker names.
    """
    matrix = adata.X if layer is None else adata.layers[layer]
    values = pd.DataFrame(
        np.asarray(matrix), columns=list(adata.var["short"]), index=adata.obs_names
    )
    values[key] = adata.obs[key].to_numpy()
    return values.groupby(key, observed=True).median()


def two_level_threshold(values: np.ndarray) -> float:
    """Split one marker's cluster medians into a low and a high group.

    The split maximises the variance between the two groups, which is Otsu's
    rule. It is used rather than a two component Gaussian mixture because a
    mixture assigns by posterior probability, and where the two components have
    different variances the wider one claims values on the far side of the
    narrower one. The groups then overlap, and a gap between them comes out
    negative. Measured on this panel that put CD16 at a separation of -0.399,
    which is not a quantity that exists.

    A threshold keeps the two groups contiguous by construction, and it needs no
    seed.

    Args:
        values: The cluster medians of one marker.

    Returns:
        The threshold. Values above it are the high group.
    """
    ordered = np.sort(np.asarray(values, dtype=float))
    if ordered.size < 2 or ordered[0] == ordered[-1]:
        return float(ordered[-1])
    best_split, best_variance = ordered[0], -1.0
    for index in range(1, ordered.size):
        low, high = ordered[:index], ordered[index:]
        variance = (
            low.size * high.size * (low.mean() - high.mean()) ** 2
        ) / ordered.size ** 2
        if variance > best_variance:
            best_variance, best_split = variance, ordered[index - 1]
    return float(best_split)


def marker_separation(
    medians: pd.DataFrame, seed: int = 42, min_range: float = 0.5
) -> pd.DataFrame:
    """Measure how cleanly each marker splits the clusters into two levels.

    A marker that takes one level across every cluster carries no information
    about which cluster is which, and letting it vote turns a label into a coin
    toss. The measure is the gap between the dimmest cluster of the high group
    and the brightest of the low group, divided by the full range of the marker,
    so it does not depend on how bright the fluorochrome is.

    The ratio alone is not enough. It has no scale, so a marker whose cluster
    medians all sit within a hundredth of a unit still scores well if those
    points happen to fall in two tight groups. `min_range` rejects that case,
    because a marker that moves by less than half an arcsinh unit across every
    cluster is measuring noise.

    Args:
        medians: The output of `cluster_medians`.
        seed: Accepted for a stable signature. The split is deterministic.
        min_range: The smallest spread across clusters, on the transformed
            scale, that a marker must cover to score above zero.

    Returns:
        A table with `marker`, `positive_clusters`, `gap` and `separation`,
        sorted from the best separated marker to the worst. `gap` is never
        negative, because the two groups are separated by a threshold.
    """
    rows = []
    for column in medians.columns:
        values = medians[column].to_numpy(dtype=float)
        spread = values.max() - values.min()
        if np.unique(values).size < 2 or spread < min_range:
            rows.append((column, 0, 0.0, 0.0))
            continue
        threshold = two_level_threshold(values)
        positive = values[values > threshold]
        negative = values[values <= threshold]
        gap = float(positive.min() - negative.max())
        rows.append((column, int(positive.size), gap, gap / spread))
    return (
        pd.DataFrame(rows,
                     columns=["marker", "positive_clusters", "gap", "separation"])
        .sort_values("separation", ascending=False)
        .reset_index(drop=True)
    )


def binarise_medians(medians: pd.DataFrame, seed: int = 42) -> pd.DataFrame:
    """Call each marker positive or negative in each cluster.

    A cluster median is an average over thousands of events, so the two levels a
    marker takes across clusters separate even where the single event
    distribution does not. Splitting the cluster medians of one marker into two
    groups is therefore a far easier problem than splitting its events, and it
    is the problem a reader solves by eye when looking at a heatmap.

    Args:
        medians: The output of `cluster_medians`.
        seed: Accepted for a stable signature. The split is deterministic.

    Returns:
        A table of the same shape holding 0 for negative, 1 for positive and 2
        for the brighter half of the positive clusters, which is what a `high`
        expectation asks for.
    """
    called = pd.DataFrame(0, index=medians.index, columns=medians.columns,
                          dtype=int)
    for column in medians.columns:
        values = medians[column].to_numpy(dtype=float)
        if np.unique(values).size < 2:
            continue
        positive = values > two_level_threshold(values)
        called.loc[positive, column] = 1
        if positive.sum() > 1:
            cut = np.median(values[positive])
            called.loc[positive & (values > cut), column] = 2
    return called


def cluster_parents(
    medians: pd.DataFrame, cd3_cut: float, cd3_marker: str = "CD3"
) -> pd.Series:
    """Assign each cluster to the side of the CD3 gate its median falls on.

    The cut comes from the R half, which fits it at the density minimum between
    two separated modes and reports which rule placed it. Reusing it here means
    the two languages draw the same boundary, and it removes the largest source
    of cross talk in the labelling: a natural killer cell definition scoring
    against a T cell cluster.

    Args:
        medians: The output of `cluster_medians`.
        cd3_cut: The CD3 cut from the R half, on the same transformed scale.
        cd3_marker: The column holding CD3.

    Returns:
        A series indexed by cluster, holding `cd3_positive` or `cd3_negative`.

    Raises:
        KeyError: If the medians carry no CD3 column.
    """
    if cd3_marker not in medians.columns:
        raise KeyError(
            f"The cluster medians carry no '{cd3_marker}' column, so no cluster "
            f"can be placed against the CD3 cut."
        )
    return pd.Series(
        np.where(medians[cd3_marker] > cd3_cut, "cd3_positive", "cd3_negative"),
        index=medians.index,
        name="parent",
    )


def annotate(
    medians: pd.DataFrame,
    definitions: pd.DataFrame,
    parents: pd.Series | None = None,
    omit: list[str] | None = None,
    min_margin: float = 0.0,
    seed: int = 42,
) -> pd.DataFrame:
    """Label each cluster by scoring it against the cell type definitions.

    Every marker is first called positive or negative across the clusters by
    `binarise_medians`. A definition then scores the fraction of the markers it
    names that the cluster matches, so a definition naming two markers and one
    naming seven are on the same scale. A `high` expectation is met only by the
    brighter half of the positive clusters, which is what separates a CD56
    bright cluster from a CD56 positive one.

    When `parents` is given, a cluster is scored only against the definitions of
    its own parent. Without it the labelling is unstable: a filter that dropped
    the markers which fail to separate the clusters changed which markers
    survived from one clustering run to the next, and with them the labels. On
    two runs of this deposit that moved CD8 from a separation of 0.236 to 0.092
    and turned a set of CD4 T cell clusters into Vd1 gamma delta clusters
    covering 40 percent of the lymphocytes.

    Args:
        medians: The output of `cluster_medians`.
        definitions: The cell type table, with `cell_type`, an optional
            `parent`, one column per marker holding `pos`, `neg`, `high` or
            nothing, and a `note`.
        parents: The output of `cluster_parents`, or `None` to score every
            definition against every cluster.
        omit: Markers that carry no usable threshold, named rather than found by
            a rule. `gating/omip58_unusable_markers.csv` holds the list for this
            deposit with the measurement behind each entry. A computed filter
            was tried first and rejected: which markers it dropped changed from
            one clustering run to the next, and the labels changed with it.
            A definition is scored only when at least two of the markers it
            names survive. One surviving marker cannot separate a definition
            from the others of its parent, and scoring it anyway lets a
            population take clusters that belong to something else.
        min_margin: A cluster whose best score beats the runner up by less than
            this is left unlabelled.
        seed: Accepted for a stable signature. The scoring is deterministic.

    Returns:
        A table with `cluster`, `parent`, `cell_type`, `score`, `runner_up`,
        `margin` and `markers_named`.

    Raises:
        ValueError: If the two tables share no marker column.
    """
    omit = set(omit or ())
    marker_columns = [
        column
        for column in definitions.columns
        if column not in ("cell_type", "parent", "note")
        and column in medians.columns
    ]
    if not marker_columns:
        raise ValueError(
            "The definitions and the expression table share no marker column. "
            f"The definitions name {list(definitions.columns)} and the table holds "
            f"{list(medians.columns)}."
        )

    called = binarise_medians(medians[marker_columns], seed=seed)

    scores = pd.DataFrame(
        index=called.index, columns=definitions["cell_type"], dtype=float
    )
    named = {}
    for row in definitions.itertuples():
        matched = pd.Series(0.0, index=called.index)
        count = 0
        for column in marker_columns:
            rule = getattr(row, column, "")
            if not isinstance(rule, str) or rule not in ("pos", "neg", "high"):
                continue
            if column in omit:
                continue
            state = called[column]
            if rule == "pos":
                met = state >= 1
            elif rule == "high":
                met = state == 2
            else:
                met = state == 0
            matched = matched + met.astype(float)
            count += 1
        named[row.cell_type] = count
        if count < 2:
            scores[row.cell_type] = np.nan
            continue
        scores[row.cell_type] = matched / count

    if parents is not None and "parent" in definitions.columns:
        owner = dict(zip(definitions["cell_type"], definitions["parent"]))
        aligned = parents.reindex(scores.index)
        for cell_type in scores.columns:
            wrong_parent = aligned != owner[cell_type]
            scores.loc[wrong_parent, cell_type] = np.nan

    values = scores.to_numpy(dtype=float)
    usable = ~np.isnan(values)
    if not usable.any(axis=1).all():
        empty = scores.index[~usable.any(axis=1)]
        raise ValueError(
            f"No definition can be scored for cluster(s) {list(empty)}. Either "
            "no cell type belongs to their parent, or every cell type that does "
            "lost too many markers to `omit`."
        )
    filled = np.where(usable, values, -np.inf)
    order = np.argsort(-filled, axis=1)
    names = np.asarray(scores.columns)
    best = names[order[:, 0]]
    top = filled[np.arange(filled.shape[0]), order[:, 0]]
    second = filled[np.arange(filled.shape[0]), order[:, 1]]
    runner_up = np.where(np.isfinite(second), names[order[:, 1]], "none")
    second = np.where(np.isfinite(second), second, 0.0)
    margin = top - second
    labelled = np.where(margin >= min_margin, best, "unlabelled")

    return pd.DataFrame(
        {
            "cluster": scores.index,
            "parent": (parents.reindex(scores.index).to_numpy()
                       if parents is not None else "all"),
            "cell_type": labelled,
            "score": top,
            "runner_up": runner_up,
            "margin": margin,
            "markers_named": [named[name] for name in best],
        }
    ).reset_index(drop=True)


def population_frequencies(
    adata: ad.AnnData, labels: pd.DataFrame, key: str = "leiden"
) -> pd.DataFrame:
    """Count each labelled cell type per donor, as a percentage of the parent.

    Args:
        adata: A clustered object.
        labels: The output of `annotate`.
        key: The `.obs` column holding the cluster label.

    Returns:
        A table with `donor`, `cell_type`, `events`, `parent_events` and
        `percent_of_parent`.
    """
    mapping = dict(zip(labels["cluster"].astype(str), labels["cell_type"]))
    frame = pd.DataFrame(
        {
            "donor": adata.obs["donor"].astype(str).to_numpy(),
            "cell_type": [mapping.get(str(value), "unlabelled")
                          for value in adata.obs[key]],
        }
    )
    counted = (
        frame.groupby(["donor", "cell_type"], observed=True)
        .size()
        .reset_index(name="events")
    )
    totals = counted.groupby("donor", observed=True)["events"].transform("sum")
    counted["parent_events"] = totals
    counted["percent_of_parent"] = 100 * counted["events"] / totals
    return counted.sort_values(["donor", "cell_type"]).reset_index(drop=True)
