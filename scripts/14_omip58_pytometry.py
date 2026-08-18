"""The Python half of the OMIP-058 workflow.

`scripts/13_omip58_prepare.R` computes the spillover matrix from the deposited
single stains, applies it, removes the debris and the dead cells, splits the T
cells from the rest, and writes each population to `output/omip58/handoff` as
FCS. This script reads those files with pytometry and looks for the subsets that
the R half could not resolve.

The division is not arbitrary. Compensation and the scatter, viability and CD3
boundaries are one dimensional problems with two separated modes, and a fitted
cut settles each of them. Inside the T cells the panel's rare markers are
unimodal and the deposit carries no fluorescence minus one control, so a fitted
cut returns a number that is wrong rather than no number at all.
`output/omip58/one_dimensional_cuts.csv` records what those cuts do. Clustering
reads every marker at once, which is the part that does not reduce to a
threshold.

Run it through uv, never through pip:
    uv run --project python python scripts/14_omip58_pytometry.py
"""

from __future__ import annotations

import argparse
import sys
import warnings
from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import anndata as ad
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
from sklearn.mixture import GaussianMixture

sys.path.insert(0, str(Path(__file__).resolve().parent))

from omip58_analysis import (
    annotate,
    cluster,
    cluster_medians,
    cluster_parents,
    marker_separation,
    population_frequencies,
    read_handoff,
    subsample,
)

OUTPUT_DIR = Path("output") / "omip58"
HANDOFF_DIR = OUTPUT_DIR / "handoff"
GATING_DIR = Path("gating")
EVENTS_PER_DONOR = 100_000
COFACTOR = 150.0
SEED = 42
CLUSTERED_PATH = OUTPUT_DIR / "clustered.h5ad"


def write(frame: pd.DataFrame, name: str) -> pd.DataFrame:
    """Write a table into the output folder and return it unchanged.

    Args:
        frame: The table to write.
        name: The file name inside the output folder.

    Returns:
        The table, so a call can be chained.
    """
    frame.to_csv(OUTPUT_DIR / name, index=False)
    return frame


def bimodality(values: np.ndarray, seed: int = SEED) -> tuple[bool, float]:
    """Ask whether one Gaussian or two describe a set of values better.

    Args:
        values: The values to fit.
        seed: The seed for the fit.

    Returns:
        Whether two components win, and the drop in the Bayesian information
        criterion from one component to two. A positive drop favours two.
    """
    column = np.asarray(values, dtype=float).reshape(-1, 1)
    if column.shape[0] < 100:
        return False, float("nan")
    one = GaussianMixture(1, random_state=seed).fit(column).bic(column)
    two = GaussianMixture(2, random_state=seed).fit(column).bic(column)
    return bool(two < one), float(one - two)


def verdict(passed: bool | None) -> str:
    """Turn a test outcome into the word the claim table carries.

    Args:
        passed: The outcome, or `None` when the deposit cannot settle it.

    Returns:
        One of `reproduced`, `not reproduced` or `unresolved`.
    """
    if passed is None:
        return "unresolved"
    return "reproduced" if passed else "not reproduced"


def main() -> int:
    """Run the clustering, judge the claims and write every table.

    Returns:
        0 when the run completed.
    """
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--reuse-clustering", action="store_true",
        help="read the saved clustering instead of building the graph again"
    )
    arguments = parser.parse_args()

    warnings.filterwarnings("ignore", category=FutureWarning)
    sc.settings.verbosity = 1
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    print("Part 1: reading what the R half wrote")
    manifest = pd.read_csv(HANDOFF_DIR / "manifest.csv")
    print(manifest[["file", "donor", "population", "events"]].to_string(index=False))

    if arguments.reuse_clustering and CLUSTERED_PATH.exists():
        sampled = ad.read_h5ad(CLUSTERED_PATH)
        print(f"  reusing {CLUSTERED_PATH}, {sampled.n_obs:,} clustered events")
    else:
        adata = read_handoff(HANDOFF_DIR, cofactor=COFACTOR,
                             populations=["live_lymphocytes"])
        print(f"  live lymphocytes read: {adata.n_obs:,} events, "
              f"{adata.n_vars} markers")

        sampled = subsample(adata, per_donor=EVENTS_PER_DONOR, seed=SEED)
        print(f"  clustering on {sampled.n_obs:,} events, "
              f"{EVENTS_PER_DONOR:,} per donor")

        print("\nPart 2: clustering")
        sampled = cluster(sampled, resolution=1.0, n_neighbors=30, seed=SEED)
        sampled.write_h5ad(CLUSTERED_PATH)

    n_clusters = sampled.obs["leiden"].nunique()
    print(f"  leiden found {n_clusters} clusters")

    medians = cluster_medians(sampled)
    write(medians.reset_index(), "cluster_medians.csv")

    definitions = pd.read_csv(GATING_DIR / "omip58_cell_type_definitions.csv")
    marker_columns = [
        column
        for column in definitions.columns
        if column not in ("cell_type", "parent", "note")
        and column in medians.columns
    ]
    separation = marker_separation(medians[marker_columns])
    write(separation, "marker_separation.csv")
    print("\n  how well each marker separates the clusters:")
    print(separation.to_string(index=False, float_format=lambda x: f"{x:.3f}"))

    # The CD3 cut comes from the R half rather than being found again here, so
    # both languages draw one boundary and a natural killer cell definition
    # never scores against a T cell cluster.
    cuts = pd.read_csv(OUTPUT_DIR / "gate_cuts.csv")
    cd3_cuts = cuts.loc[cuts["marker"] == "CD3", "cut"]
    cd3_cut = float(cd3_cuts.mean())
    print(f"\n  CD3 cut from the R half: {cd3_cut:.3f} on the arcsinh scale, "
          f"the mean over {len(cd3_cuts)} donors")
    parents = cluster_parents(medians, cd3_cut)
    print(parents.value_counts().to_string())

    # The markers that carry no usable threshold are named rather than found by
    # a rule, with the measurement behind each one. Part 4 of the R half is the
    # evidence.
    unusable = pd.read_csv(GATING_DIR / "omip58_unusable_markers.csv")
    omit = list(unusable["marker"])
    print(f"\n  markers with no usable threshold: {', '.join(omit)}")

    labels = annotate(medians, definitions, parents=parents, omit=omit)
    scored = set(labels.loc[labels["markers_named"] > 0, "cell_type"])
    dropped = sorted(set(definitions["cell_type"]) - scored)
    if dropped:
        print(f"  populations no cluster can be scored against: "
              f"{', '.join(dropped)}")
    sizes = sampled.obs["leiden"].value_counts().rename("events")
    labels = labels.merge(sizes, left_on="cluster", right_index=True, how="left")
    labels["percent_of_sample"] = 100 * labels["events"] / sampled.n_obs
    write(labels, "cluster_labels.csv")
    print(labels.to_string(index=False, float_format=lambda x: f"{x:.3f}"))

    frequencies = population_frequencies(sampled, labels)
    write(frequencies, "cluster_frequencies.csv")
    print("\n" + frequencies.to_string(index=False,
                                       float_format=lambda x: f"{x:.3f}"))

    print("\nPart 3: figures")
    mapping = dict(zip(labels["cluster"].astype(str), labels["cell_type"]))
    sampled.obs["cell_type"] = [mapping.get(str(value), "unlabelled")
                                for value in sampled.obs["leiden"]]

    with plt.rc_context({"figure.figsize": (9, 7)}):
        sc.pl.umap(sampled, color="cell_type", show=False, legend_loc="right margin",
                   title="OMIP-058 live lymphocytes, labelled by cluster")
        plt.savefig(OUTPUT_DIR / "umap_cell_type.png", dpi=150, bbox_inches="tight")
        plt.close("all")

    with plt.rc_context({"figure.figsize": (9, 7)}):
        sc.pl.umap(sampled, color="donor", show=False,
                   title="OMIP-058 live lymphocytes, by donor")
        plt.savefig(OUTPUT_DIR / "umap_donor.png", dpi=150, bbox_inches="tight")
        plt.close("all")

    ordered = medians.loc[labels.sort_values("cell_type")["cluster"]]
    ordered.index = labels.sort_values("cell_type")["cell_type"].to_numpy()
    figure, axis = plt.subplots(figsize=(11, 6))
    image = axis.imshow(ordered.to_numpy(), aspect="auto", cmap="viridis")
    axis.set_xticks(range(ordered.shape[1]))
    axis.set_xticklabels(ordered.columns, rotation=90, fontsize=8)
    axis.set_yticks(range(ordered.shape[0]))
    axis.set_yticklabels(ordered.index, fontsize=8)
    axis.set_title("Median expression per cluster, arcsinh scale, cofactor 150")
    figure.colorbar(image, ax=axis, shrink=0.8)
    figure.tight_layout()
    figure.savefig(OUTPUT_DIR / "cluster_heatmap.png", dpi=150)
    plt.close("all")

    print("\nPart 4: the claims")
    claims = pd.read_csv(GATING_DIR / "omip58_paper_claims.csv")
    rows: list[tuple[str, str, str]] = []

    def observed_percent(cell_type: str, donor: str | None = None) -> float:
        subset = frequencies[frequencies["cell_type"] == cell_type]
        if donor is not None:
            subset = subset[subset["donor"] == donor]
        if subset.empty:
            return 0.0
        return float(subset["percent_of_parent"].mean())

    donors = sorted(frequencies["donor"].unique())
    found = set(labels["cell_type"])

    def unmeasurable(*names: str) -> str | None:
        """Name the populations a claim needs that no cluster carries."""
        missing = [name for name in names if name not in found]
        if not missing:
            return None
        return ", ".join(missing)

    reason = unmeasurable("iNKT cells")
    if reason:
        rows.append((
            "1",
            f"no cluster carries {reason}",
            "unresolved",
        ))
    else:
        inkt = [observed_percent("iNKT cells", donor) for donor in donors]
        rows.append(("1", f"{max(inkt):.3f} percent of live lymphocytes at most",
                     verdict(max(inkt) <= 1.0)))

    reason = unmeasurable("Vd2Vg9 gamma delta T cells", "Vd1 gamma delta T cells")
    vd2 = observed_percent("Vd2Vg9 gamma delta T cells")
    vd1 = observed_percent("Vd1 gamma delta T cells")
    if reason:
        rows.append((
            "2",
            f"no cluster carries {reason}",
            "unresolved",
        ))
    else:
        rows.append(("2", f"Vd2Vg9 {vd2:.3f} percent against Vd1 {vd1:.3f} percent",
                     verdict(vd2 > vd1)))

    gd_clusters = labels.loc[
        labels["cell_type"] == "Vd2Vg9 gamma delta T cells", "cluster"
    ].astype(str)
    if unmeasurable("Vd2Vg9 gamma delta T cells") or len(gd_clusters) == 0:
        rows.append(("3", "no Vd2Vg9 cluster to test", "unresolved"))
    else:
        inside = sampled[sampled.obs["leiden"].astype(str).isin(gd_clusters)]
        column = np.asarray(inside.layers["arcsinh"])[
            :, list(inside.var["short"]).index("Vd2")
        ]
        two_wins, drop = bimodality(column)
        rows.append((
            "3",
            (f"two components beat one by {drop:.0f} in the Bayesian "
             f"information criterion"),
            verdict(two_wins),
        ))

    rows.append(("4", f"MAIT cells carry a cluster: {'MAIT cells' in found}",
                 verdict("MAIT cells" in found)))

    unconventional = ("iNKT cells", "MAIT cells", "Vd1 gamma delta T cells",
                      "Vd2Vg9 gamma delta T cells")
    reason = unmeasurable(*unconventional)
    if reason:
        rows.append((
            "5",
            (f"the unconventional subsets are not all measurable, so the "
             f"share cannot be computed: no cluster carries {reason}"),
            "unresolved",
        ))
    else:
        conventional = observed_percent("CD4 T cells") + observed_percent("CD8 T cells")
        t_total = conventional + sum(observed_percent(name)
                                     for name in unconventional)
        share = 100 * conventional / t_total if t_total > 0 else 0.0
        rows.append(("5", f"{share:.1f} percent of the labelled T cells",
                     verdict(share > 50)))

    nk_names = ("Early NK cells", "Mature NK cells", "Terminal NK cells")
    nk_subsets = [name for name in nk_names if name in found]
    rows.append(("6", f"{len(nk_subsets)} of three NK subsets carry a cluster",
                 verdict(len(nk_subsets) == 3)))

    nk_values = {name: observed_percent(name) for name in nk_names}
    largest = max(nk_values, key=nk_values.get) if any(nk_values.values()) else None
    if len(nk_subsets) < 3:
        rows.append((
            "7",
            (f"only {len(nk_subsets)} of the three subsets carry a cluster, so "
             f"the largest cannot be named"),
            "unresolved",
        ))
    else:
        rows.append((
            "7",
            (f"largest NK subset is {largest}, "
             f"{nk_values.get(largest, 0):.3f} percent"),
            verdict(largest == "Mature NK cells"),
        ))

    gd_labels = ["Vd1 gamma delta T cells", "Vd2Vg9 gamma delta T cells"]
    gd_rows = labels[labels["cell_type"].isin(gd_labels)]
    conv_rows = labels[labels["cell_type"].isin(["CD4 T cells", "CD8 T cells"])]
    if unmeasurable(*gd_labels) or gd_rows.empty or conv_rows.empty:
        rows.append(("8", "no gamma delta cluster to test", "unresolved"))
    else:
        gd_cells = sampled[sampled.obs["leiden"].astype(str).isin(
            gd_rows["cluster"].astype(str))]
        conv_cells = sampled[sampled.obs["leiden"].astype(str).isin(
            conv_rows["cluster"].astype(str))]
        shorts = list(sampled.var["short"])
        gd_values = np.asarray(gd_cells.layers["arcsinh"])
        conv_values = np.asarray(conv_cells.layers["arcsinh"])
        fraction = float(np.mean(
            (gd_values[:, shorts.index("CD16")] >
             np.quantile(conv_values[:, shorts.index("CD16")], 0.999)) &
            (gd_values[:, shorts.index("CD56")] >
             np.quantile(conv_values[:, shorts.index("CD56")], 0.999))
        ))
        rows.append((
            "8",
            (f"{100 * fraction:.3f} percent of gamma delta cells are above "
             f"both thresholds"),
            verdict(fraction > 0),
        ))

    if conv_rows.empty:
        rows.append(("9", "no conventional T cell cluster to test", "unresolved"))
    else:
        conv_cells = sampled[sampled.obs["leiden"].astype(str).isin(
            conv_rows["cluster"].astype(str))]
        shorts = list(sampled.var["short"])
        conv_values = np.asarray(conv_cells.layers["arcsinh"])
        ccr7 = conv_values[:, shorts.index("CCR7")]
        cd45ra = conv_values[:, shorts.index("CD45RA")]
        ccr7_two, _ = bimodality(ccr7)
        cd45ra_two, _ = bimodality(cd45ra)
        quadrants = pd.crosstab(ccr7 > np.median(ccr7), cd45ra > np.median(cd45ra))
        rows.append((
            "9",
            (f"CCR7 bimodal {ccr7_two}, CD45RA bimodal {cd45ra_two}, "
             f"{int((quadrants > 0).to_numpy().sum())} quadrants occupied"),
            verdict(ccr7_two and cd45ra_two),
        ))

    expected_types = set(definitions["cell_type"]) - {"HLA-DR positive non T cells"}
    rows.append((
        "10",
        (f"{len(found & expected_types)} of {len(expected_types)} "
         f"populations carry a cluster"),
        verdict(found >= expected_types),
    ))

    verdicts = pd.DataFrame(rows, columns=["claim_id", "observed", "verdict"])
    verdicts["claim_id"] = verdicts["claim_id"].astype(int)
    verdicts = claims.merge(verdicts, on="claim_id", how="left")
    write(verdicts[["claim_id", "short_name", "expected", "observed", "verdict"]],
          "claim_verdicts.csv")
    print(verdicts[["claim_id", "short_name", "observed", "verdict"]].to_string(
        index=False))

    counts = verdicts["verdict"].value_counts().rename_axis("verdict")
    write(counts.reset_index(name="claims"), "verdict_counts.csv")
    print("\n" + counts.to_string())

    print(f"\nDone. Tables and figures are in {OUTPUT_DIR}.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
