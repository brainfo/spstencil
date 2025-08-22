"""
spstencil
=========

Spatial transcriptomics utilities for subsetting, splitting, aggregating, and cropping
AnnData objects and related files.

Top-level conveniences mirror a Scanpy-style API for quick access, while full
functionality is organized under submodules `io`, `ops`, `st`, and `_utils`.

Example
-------

>>> import spstencil as sps
>>> adata = sps.load_anndata("sample.h5ad")
>>> subset = sps.subset_by_values(adata, key="tissue_type", include_values=["Epithelium"]) 
"""

from importlib.metadata import version as _pkg_version

# Re-export submodules for convenient access like `import spstencil as sps; sps.io.load_anndata`
from . import io as io  # noqa: F401
from . import ops as ops  # noqa: F401
from . import st as st  # noqa: F401
from ._utils import axes as axes  # noqa: F401
from ._utils import cyx2yxc as cyx2yxc  # noqa: F401
from . import project as project  # noqa: F401

# Common, stable top-level functions
from .io import load_anndata, save_anndata
from .ops import (
    subset_by_values,
    split_by_key,
    aggregate_cell_counts,
    save_split_by_key,
)
from .st import (
    get_spatial_coords,
    assign_tiles_from_coords,
    annotate_tissue_to_spots,
    subset_st_by_tissue,
    split_st_by_tissue,
    aggregate_st_by_cell_class,
    load_cell_predictions,
    validate_coordinate_compatibility,
    determine_crop_bounds,
    crop_spatial_data,
    crop_cell_data,
    crop_tissue_data,
    filter_cell_predictions_by_tissue,
    map_cells_to_tissue_classes,
    split_hdf5_by_tissue,
    aggregate_hdf5_by_cell_class,
    is_continuous,
    missing_ranges,
)
from .project import (
    project_points_onto_polyline,
    sample_curve,
    extract_polyline_from_image,
    unroll_adata_along_polyline,
    unroll_sdata,
    ProjectionResult,
)

__all__ = [
    # Submodules
    "io",
    "ops",
    "st",
    "axes",
    "cyx2yxc",
    "project",
    # Version
    "__version__",
    # IO
    "load_anndata",
    "save_anndata",
    # Generic ops
    "subset_by_values",
    "split_by_key",
    "aggregate_cell_counts",
    "save_split_by_key",
    # Spatial toolkit
    "get_spatial_coords",
    "assign_tiles_from_coords",
    "annotate_tissue_to_spots",
    "subset_st_by_tissue",
    "split_st_by_tissue",
    "aggregate_st_by_cell_class",
    "load_cell_predictions",
    "validate_coordinate_compatibility",
    "determine_crop_bounds",
    "crop_spatial_data",
    "crop_cell_data",
    "crop_tissue_data",
    "filter_cell_predictions_by_tissue",
    "map_cells_to_tissue_classes",
    "split_hdf5_by_tissue",
    "aggregate_hdf5_by_cell_class",
    "is_continuous",
    "missing_ranges",
    # Unroll utilities
    "project_points_onto_polyline",
    "sample_curve",
    "unroll_adata_along_polyline",
    "extract_polyline_from_image",
    "unroll_sdata",
    "ProjectionResult",
]

try:
    __version__ = _pkg_version("spstencil")
except Exception:  # pragma: no cover - during editable installs or local dev
    __version__ = "0.0.0"
