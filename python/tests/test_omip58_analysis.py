"""Tests for scripts/omip58_analysis.py.

Every test builds its own small table. None of them reads the deposit or the
handoff, because a 77 MB FCS file has no place in a test suite.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest
from omip58_analysis import (
    annotate,
    binarise_medians,
    cluster_medians,
    definition_coverage,
    marker_separation,
    marker_tokens,
    population_frequencies,
    short_marker_names,
    subsample,
)


def build_medians() -> pd.DataFrame:
    """Return four clusters whose identity is obvious from three markers."""
    return pd.DataFrame(
        {
            "CD3": [4.0, 4.1, 0.1, 0.0],
            "CD4": [3.8, 0.1, 0.0, 0.1],
            "CD8": [0.0, 3.9, 0.1, 0.2],
            "CD56": [0.1, 0.0, 4.2, 4.0],
            "CD16": [0.0, 0.1, 0.2, 4.1],
            "noise": [1.0, 1.1, 0.9, 1.05],
        },
        index=["0", "1", "2", "3"],
    )


def build_definitions() -> pd.DataFrame:
    """Return definitions for the four clusters of `build_medians`."""
    return pd.DataFrame(
        {
            "cell_type": ["CD4 T cells", "CD8 T cells", "Early NK cells",
                          "Mature NK cells"],
            "CD3": ["pos", "pos", "neg", "neg"],
            "CD4": ["pos", "neg", "neg", "neg"],
            "CD8": ["neg", "pos", "neg", "neg"],
            "CD56": ["neg", "neg", "pos", "pos"],
            "CD16": ["neg", "neg", "neg", "pos"],
            "note": ["", "", "", ""],
        }
    )


def test_marker_tokens_splits_on_every_separator():
    assert marker_tokens("TCR Va7_2 BV711") == ["tcr", "va7", "2", "bv711"]
    assert marker_tokens("HLA-DR PE-Cy55") == ["hla", "dr", "pe", "cy55"]


def test_a_token_match_keeps_cd16_and_cd161_apart():
    assert "cd16" in marker_tokens("CD16 BUV496")
    assert "cd16" not in marker_tokens("CD161 PE-Cy5")
    assert "cd161" in marker_tokens("CD161 PE-Cy5")


def test_short_marker_names_rejects_a_panel_missing_a_marker():
    var = pd.DataFrame(index=["CD3 BV510", "CD4 BUV805"])
    with pytest.raises(ValueError, match="matches 0 markers"):
        short_marker_names(var)


FULL_PANEL = [
    "Live Dead UV Blue", "CD3 BV510", "CD4 BUV805", "CD8 BV570",
    "CD16 BUV496", "CD56 BUV563", "CD161 PE-Cy5", "HLA-DR PE-Cy55",
    "CCR7 BUV395", "CD45RA Ax700", "CD27 BV786", "CD28 BV650",
    "CD95 BUV737", "CD1d PBS57 tetramer APC", "TCR Va7_2 BV711",
    "TCR Vd1 FITC", "TCR Vd2 PE-CF594", "TCR Vg9 PE",
]


def test_short_marker_names_resolves_the_whole_panel():
    resolved = short_marker_names(pd.DataFrame(index=FULL_PANEL))
    assert resolved["CD3 BV510"] == "CD3"
    assert resolved["CD161 PE-Cy5"] == "CD161"
    assert resolved["CD16 BUV496"] == "CD16"
    assert len(resolved) == len(FULL_PANEL)


def test_short_marker_names_rejects_a_token_matching_twice():
    var = pd.DataFrame(index=[*FULL_PANEL, "CD3 BUV661"])
    with pytest.raises(ValueError, match="matches 2 markers"):
        short_marker_names(var)


def test_marker_separation_ranks_an_informative_marker_above_noise():
    separation = marker_separation(build_medians()).set_index("marker")
    assert separation.loc["CD3", "separation"] > 0.3
    assert separation.loc["noise", "separation"] == 0.0


def test_binarise_medians_calls_the_bright_clusters_positive():
    called = binarise_medians(build_medians())
    assert called.loc["0", "CD4"] >= 1
    assert called.loc["1", "CD4"] == 0
    assert called.loc["0", "CD3"] >= 1
    assert called.loc["2", "CD3"] == 0


def test_annotate_labels_each_cluster_with_its_own_definition():
    labels = annotate(build_medians(), build_definitions()).set_index("cluster")
    assert labels.loc["0", "cell_type"] == "CD4 T cells"
    assert labels.loc["1", "cell_type"] == "CD8 T cells"
    assert labels.loc["3", "cell_type"] == "Mature NK cells"
    assert (labels["score"] == 1.0).all()


def test_annotate_leaves_a_close_call_unlabelled():
    labels = annotate(build_medians(), build_definitions(), min_margin=0.9)
    assert set(labels["cell_type"]) == {"unlabelled"}


def test_annotate_rejects_tables_sharing_no_marker():
    medians = pd.DataFrame({"other": [1.0, 2.0]}, index=["0", "1"])
    with pytest.raises(ValueError, match="share no marker column"):
        annotate(medians, build_definitions())


def test_annotate_drops_a_definition_whose_markers_stopped_separating():
    # Every marker but CD3 is made uninformative, so no definition keeps the two
    # markers it needs and the call has nothing left to score.
    medians = build_medians()
    for column in ("CD4", "CD8", "CD56", "CD16"):
        medians[column] = [1.0, 1.01, 0.99, 1.0]
    with pytest.raises(ValueError, match="No definition kept enough markers"):
        annotate(medians, build_definitions(), min_separation=0.3)


def test_definition_coverage_names_the_markers_a_definition_lost():
    medians = build_medians()
    medians["CD16"] = [1.0, 1.01, 0.99, 1.0]
    coverage = definition_coverage(medians, build_definitions(),
                                   min_separation=0.3).set_index("cell_type")
    assert "CD16" in coverage.loc["Mature NK cells", "lost"]
    assert coverage.loc["CD4 T cells", "scorable"]


def test_subsample_takes_the_same_count_from_each_donor():
    import anndata as ad

    adata = ad.AnnData(
        np.zeros((300, 2), dtype=np.float32),
        obs=pd.DataFrame({"donor": ["a"] * 200 + ["b"] * 100}),
    )
    sampled = subsample(adata, per_donor=50, seed=1)
    assert sampled.n_obs == 100
    assert (sampled.obs["donor"].value_counts() == 50).all()


def test_subsample_takes_everything_when_a_donor_has_too_few():
    import anndata as ad

    adata = ad.AnnData(
        np.zeros((30, 2), dtype=np.float32),
        obs=pd.DataFrame({"donor": ["a"] * 20 + ["b"] * 10}),
    )
    sampled = subsample(adata, per_donor=50, seed=1)
    assert sampled.n_obs == 30


def test_subsample_is_repeatable():
    import anndata as ad

    adata = ad.AnnData(
        np.arange(200, dtype=np.float32).reshape(100, 2),
        obs=pd.DataFrame({"donor": ["a"] * 100}),
    )
    first = subsample(adata, per_donor=20, seed=7)
    second = subsample(adata, per_donor=20, seed=7)
    assert np.array_equal(first.X, second.X)


def test_cluster_medians_reads_the_named_layer():
    import anndata as ad

    adata = ad.AnnData(
        np.zeros((4, 2), dtype=np.float32),
        obs=pd.DataFrame({"leiden": ["0", "0", "1", "1"]}),
        var=pd.DataFrame({"short": ["CD3", "CD4"]}, index=["a", "b"]),
    )
    adata.layers["arcsinh"] = np.array(
        [[1.0, 5.0], [3.0, 7.0], [10.0, 0.0], [12.0, 2.0]], dtype=np.float32
    )
    medians = cluster_medians(adata)
    assert medians.loc["0", "CD3"] == pytest.approx(2.0)
    assert medians.loc["1", "CD4"] == pytest.approx(1.0)


def test_population_frequencies_sum_to_one_hundred_per_donor():
    import anndata as ad

    adata = ad.AnnData(
        np.zeros((10, 1), dtype=np.float32),
        obs=pd.DataFrame(
            {
                "donor": ["a"] * 6 + ["b"] * 4,
                "leiden": ["0", "0", "0", "1", "1", "1", "0", "0", "1", "1"],
            }
        ),
    )
    labels = pd.DataFrame({"cluster": ["0", "1"],
                           "cell_type": ["CD4 T cells", "CD8 T cells"]})
    frequencies = population_frequencies(adata, labels)
    totals = frequencies.groupby("donor")["percent_of_parent"].sum()
    assert totals.round(6).eq(100.0).all()
    assert frequencies.loc[
        (frequencies["donor"] == "a")
        & (frequencies["cell_type"] == "CD4 T cells"), "events"
    ].item() == 3
