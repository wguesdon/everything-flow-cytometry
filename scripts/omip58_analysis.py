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
from sklearn.mixture import GaussianMixture

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


def marker_separation(medians: pd.DataFrame, seed: int = 42) -> pd.DataFrame:
    """Measure how cleanly each marker splits the clusters into two levels.

    A marker that takes one level across every cluster carries no information
    about which cluster is which, and letting it vote turns a label into a coin
    toss. The measure is the gap between the dimmest positive cluster and the
    brightest negative one, divided by the full range of the marker, so it does
    not depend on how bright the fluorochrome is.

    Args:
        medians: The output of `cluster_medians`.
        seed: The seed for the two group fit.

    Returns:
        A table with `marker`, `positive_clusters`, `gap` and `separation`,
        sorted from the best separated marker to the worst.
    """
    rows = []
    for column in medians.columns:
        values = medians[column].to_numpy(dtype=float)
        spread = values.max() - values.min()
        if np.unique(values).size < 2 or spread == 0:
            rows.append((column, 0, 0.0, 0.0))
            continue
        model = GaussianMixture(2, random_state=seed, n_init=5).fit(
            values.reshape(-1, 1)
        )
        assignment = model.predict(values.reshape(-1, 1))
        bright = int(np.argmax(model.means_.ravel()))
        positive = values[assignment == bright]
        negative = values[assignment != bright]
        gap = float(positive.min() - negative.max())
        rows.append((column, int(positive.size), gap, gap / spread))
    return (
        pd.DataFrame(rows,
                     columns=["marker", "positive_clusters", "gap", "separation"])
        .sort_values("separation", ascending=False)
        .reset_index(drop=True)
    )


def definition_coverage(
    medians: pd.DataFrame,
    definitions: pd.DataFrame,
    min_separation: float = 0.10,
    min_markers: float = 0.5,
    seed: int = 42,
) -> pd.DataFrame:
    """Say which cell type definitions still have enough markers to be scored.

    A claim about a population that no definition can score is unresolved, and
    it must not be read as a frequency of zero. This table is what tells the two
    cases apart.

    Args:
        medians: The output of `cluster_medians`.
        definitions: The cell type table.
        min_separation: The threshold passed to `marker_separation`.
        min_markers: The share of named markers that must survive.
        seed: The seed for the fits.

    Returns:
        A table with `cell_type`, `markers_named`, `markers_informative`,
        `lost` and `scorable`.
    """
    marker_columns = [
        column
        for column in definitions.columns
        if column not in ("cell_type", "note") and column in medians.columns
    ]
    separation = marker_separation(medians[marker_columns], seed=seed)
    informative = set(
        separation.loc[separation["separation"] >= min_separation, "marker"]
    )

    rows = []
    for row in definitions.itertuples():
        named = [
            column
            for column in marker_columns
            if isinstance(getattr(row, column, ""), str)
            and getattr(row, column) in ("pos", "neg", "high")
        ]
        kept = [column for column in named if column in informative]
        rows.append(
            {
                "cell_type": row.cell_type,
                "markers_named": len(named),
                "markers_informative": len(kept),
                "lost": " ".join(sorted(set(named) - set(kept))),
                "scorable": len(kept) >= 2 and len(kept) >= min_markers * len(named),
            }
        )
    return pd.DataFrame(rows)


def binarise_medians(medians: pd.DataFrame, seed: int = 42) -> pd.DataFrame:
    """Call each marker positive or negative in each cluster.

    A cluster median is an average over thousands of events, so the two levels a
    marker takes across clusters separate even where the single event
    distribution does not. Splitting the cluster medians of one marker into two
    groups is therefore a far easier problem than splitting its events, and it
    is the problem a reader solves by eye when looking at a heatmap.

    Args:
        medians: The output of `cluster_medians`.
        seed: The seed for the two group fit.

    Returns:
        A table of the same shape holding 0 for negative, 1 for positive and 2
        for the brighter half of the positive clusters, which is what a `high`
        expectation asks for.
    """
    called = pd.DataFrame(0, index=medians.index, columns=medians.columns,
                          dtype=int)
    for column in medians.columns:
        values = medians[column].to_numpy(dtype=float).reshape(-1, 1)
        if np.unique(values).size < 2:
            continue
        model = GaussianMixture(2, random_state=seed, n_init=3).fit(values)
        assignment = model.predict(values)
        bright = int(np.argmax(model.means_.ravel()))
        positive = assignment == bright
        called.loc[positive, column] = 1
        if positive.sum() > 1:
            cut = np.median(values[positive])
            called.loc[positive & (values.ravel() > cut), column] = 2
    return called


def annotate(
    medians: pd.DataFrame,
    definitions: pd.DataFrame,
    min_margin: float = 0.0,
    min_separation: float = 0.10,
    min_markers: float = 0.5,
    seed: int = 42,
) -> pd.DataFrame:
    """Label each cluster by scoring it against the cell type definitions.

    Every marker is first called positive or negative across the clusters by
    `binarise_medians`. A definition then scores the fraction of the markers it
    names that the cluster matches, so a definition naming two markers and one
    naming seven are on the same scale. A `high` expectation is met only by the
    brighter half of the positive clusters, which is what separates a CD56
    bright cluster from a CD56 positive one.

    Scoring against the scaled value rather than the call is what an earlier
    version did, and it fails twice. Min to max scaling puts a typical cluster
    near 0.5, so every definition scored close to every other and the margins
    ran from 0.013 to 0.27. It also labelled a CD4 T cell cluster as a Vd1
    gamma delta cluster, because a middling Vd1 median scored almost as well
    against `pos` as against `neg`.

    Args:
        medians: The output of `cluster_medians`.
        definitions: The cell type table, with a `cell_type` column, one column
            per marker holding `pos`, `neg`, `high` or nothing, and a `note`.
        min_margin: A cluster whose best score beats the runner up by less than
            this is left unlabelled.
        min_separation: A marker whose `marker_separation` falls below this does
            not vote. On this panel that removes the markers that also defeat a
            one dimensional cut.
        min_markers: The share of the markers a definition names that must
            survive the separation filter for the definition to be scored at
            all. A definition reduced to one generic marker would otherwise
            match every cluster carrying it.
        seed: The seed passed to `binarise_medians` and `marker_separation`.

    Returns:
        A table with `cluster`, `cell_type`, `score`, `runner_up`, `margin` and
        `markers_named`.

    Raises:
        ValueError: If the two tables share no marker column.
    """
    marker_columns = [
        column
        for column in definitions.columns
        if column not in ("cell_type", "note") and column in medians.columns
    ]
    if not marker_columns:
        raise ValueError(
            "The definitions and the expression table share no marker column. "
            f"The definitions name {list(definitions.columns)} and the table holds "
            f"{list(medians.columns)}."
        )

    separation = marker_separation(medians[marker_columns], seed=seed)
    informative = set(
        separation.loc[separation["separation"] >= min_separation, "marker"]
    )
    called = binarise_medians(medians[marker_columns], seed=seed)

    scores = pd.DataFrame(
        index=called.index, columns=definitions["cell_type"], dtype=float
    )
    named = {}
    for row in definitions.itertuples():
        matched = pd.Series(0.0, index=called.index)
        count = 0
        total_named = 0
        for column in marker_columns:
            rule = getattr(row, column, "")
            if not isinstance(rule, str) or rule not in ("pos", "neg", "high"):
                continue
            total_named += 1
            if column not in informative:
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
        # A definition that has lost most of its markers is no longer a
        # definition. Scoring it anyway lets a population survive on one generic
        # marker and take clusters that belong to something else.
        if count < 2 or count < min_markers * total_named:
            scores[row.cell_type] = np.nan
            continue
        scores[row.cell_type] = matched / count

    scores = scores.dropna(axis=1, how="all")
    if scores.shape[1] == 0:
        raise ValueError(
            "No definition kept enough markers to be scored. Lower "
            "min_separation or name more markers per cell type."
        )

    values = scores.to_numpy()
    order = np.argsort(-values, axis=1)
    names = np.asarray(scores.columns)
    best = names[order[:, 0]]
    runner_up = names[order[:, 1]] if values.shape[1] > 1 else best
    top = values[np.arange(values.shape[0]), order[:, 0]]
    second = (
        values[np.arange(values.shape[0]), order[:, 1]]
        if values.shape[1] > 1
        else np.zeros(values.shape[0])
    )
    margin = top - second
    labelled = np.where(margin >= min_margin, best, "unlabelled")

    return pd.DataFrame(
        {
            "cluster": scores.index,
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
