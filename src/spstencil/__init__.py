"""
spstencil: Utilities for subsetting, splitting, and aggregating h5ad files

Public API:
- load_anndata
- save_anndata
- subset_by_values
- split_by_key
- aggregate_cell_counts
"""

from .io import load_anndata, save_anndata
from .ops import subset_by_values, split_by_key, aggregate_cell_counts

__all__ = [
    "load_anndata",
    "save_anndata",
    "subset_by_values",
    "split_by_key",
    "aggregate_cell_counts",
]

