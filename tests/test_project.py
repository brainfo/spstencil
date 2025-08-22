from __future__ import annotations

import numpy as np
import pandas as pd
import anndata as ad

import spstencil as sps


def make_adata(coords: np.ndarray, genes: list[str] = ["g1", "g2"]) -> ad.AnnData:
    n = coords.shape[0]
    X = np.vstack([
        np.linspace(0, 1, n),
        np.linspace(1, 0, n),
    ]).T
    adata = ad.AnnData(X=X, obs=pd.DataFrame({"imagecol": coords[:, 0], "imagerow": coords[:, 1]}), var=pd.DataFrame(index=genes))
    return adata


def test_projection_straight_line_horizontal():
    # Points along y=0 from x=0..10
    pts = np.stack([np.linspace(0, 10, 11), np.zeros(11)], axis=1)
    line = np.array([[0, 0], [10, 0]])
    res = sps.project_points_onto_polyline(pts, line)
    assert np.allclose(res.distances, 0)
    assert np.allclose(res.normalized_positions, np.linspace(0, 1, 11))


def test_projection_straight_line_vertical():
    pts = np.stack([np.zeros(11), np.linspace(0, 10, 11)], axis=1)
    line = np.array([[0, 0], [0, 10]])
    res = sps.project_points_onto_polyline(pts, line)
    assert np.allclose(res.distances, 0)
    assert np.allclose(res.normalized_positions, np.linspace(0, 1, 11))


def test_unroll_mean_aggregation_bins():
    # Create diagonal points projected on straight line
    coords = np.stack([np.linspace(0, 10, 20), np.linspace(0, 10, 20)], axis=1)
    adata = make_adata(coords)
    line = np.array([[0, 0], [10, 10]])
    agg_df, mapping_df = sps.unroll_adata_along_polyline(adata, line, n_bins=5, agg="mean")
    assert agg_df.shape[0] == 5
    # Ensure bin coverage
    assert mapping_df["bin"].min() >= 0 and mapping_df["bin"].max() <= 4


def test_polyline_from_roi_and_projection(tmp_path):
    import shutil
    import os
    roi_src = "/mnt/run/jh/my_github/spstencil/testdata/roi.jpg"
    roi_dst = tmp_path / "roi.jpg"
    shutil.copy(roi_src, roi_dst)
    poly = sps.extract_polyline_from_image(roi_dst, threshold=64, invert=False, step=5)
    assert poly.shape[1] == 2 and poly.shape[0] >= 2
    # Project a few random points near the image center to ensure it runs
    pts = np.array([[100, 100], [200, 150], [300, 250]], dtype=float)
    res = sps.project_points_onto_polyline(pts, poly)
    assert res.normalized_positions.shape == (3,)
    assert np.all(np.isfinite(res.distances))


