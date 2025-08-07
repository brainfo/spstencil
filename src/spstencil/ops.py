from __future__ import annotations

from pathlib import Path
from typing import Dict, Iterable, List, Mapping, Sequence, Tuple

import numpy as np
import pandas as pd
import anndata as ad


def subset_by_values(
    adata: ad.AnnData,
    key: str,
    include_values: Sequence[str] | None = None,
    exclude_values: Sequence[str] | None = None,
) -> ad.AnnData:
    """Return a new AnnData subset by obs[key] membership.

    At least one of include_values or exclude_values must be provided.
    """
    if include_values is None and exclude_values is None:
        raise ValueError("Provide include_values and/or exclude_values")

    if key not in adata.obs:
        raise KeyError(f"obs lacks column '{key}'")

    obs_series = adata.obs[key].astype("string")
    mask = pd.Series(True, index=obs_series.index)

    if include_values is not None:
        include_set = pd.Index(include_values).astype("string")
        mask &= obs_series.isin(include_set)

    if exclude_values is not None:
        exclude_set = pd.Index(exclude_values).astype("string")
        mask &= ~obs_series.isin(exclude_set)

    return adata[mask.values].copy()


def split_by_key(adata: ad.AnnData, key: str) -> Dict[str, ad.AnnData]:
    """Split AnnData into a dict keyed by unique values in obs[key]."""
    if key not in adata.obs:
        raise KeyError(f"obs lacks column '{key}'")
    groups = adata.obs[key].astype("string").fillna("NA").unique().tolist()
    return {val: adata[adata.obs[key].astype("string") == val].copy() for val in groups}


def aggregate_cell_counts(
    adata: ad.AnnData,
    cell_type_key: str,
    group_key: str | None = None,
    normalize: bool = True,
) -> pd.DataFrame:
    """Aggregate counts of cell types, optionally per group.

    Returns a DataFrame of percentages if normalize else counts.
    - If group_key is None, returns one-row DataFrame across all cells.
    - Expects categorical/string columns in obs.
    """
    if cell_type_key not in adata.obs:
        raise KeyError(f"obs lacks column '{cell_type_key}'")

    obs = adata.obs.copy()
    obs[cell_type_key] = obs[cell_type_key].astype("string")
    if group_key is not None:
        if group_key not in obs:
            raise KeyError(f"obs lacks column '{group_key}'")
        obs[group_key] = obs[group_key].astype("string")

    if group_key is None:
        counts = obs[cell_type_key].value_counts().rename_axis("cell_type").to_frame("count").T
        counts.index = ["all"]
    else:
        counts = (
            obs.groupby(group_key)[cell_type_key]
            .value_counts()
            .unstack(fill_value=0)
        )

    if normalize:
        counts = counts.div(counts.sum(axis=1), axis=0) * 100.0
    return counts


def save_split_by_key(
    adata: ad.AnnData,
    key: str,
    out_dir: str | Path,
    compression: str | None = "gzip",
) -> Dict[str, Path]:
    """Split by obs[key] and save each part to out_dir/<value>.h5ad. Returns paths."""
    from .io import save_anndata

    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    parts = split_by_key(adata, key)
    result: Dict[str, Path] = {}
    for value, sub in parts.items():
        safe_value = str(value).replace("/", "_")
        out_path = out_dir / f"{safe_value}.h5ad"
        save_anndata(sub, out_path, compression=compression)
        result[value] = out_path
    return result

