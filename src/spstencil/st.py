from __future__ import annotations

from pathlib import Path
from typing import Iterable, Sequence, Tuple, Dict, Optional, List

import h5py
import numpy as np
import pandas as pd
import anndata as ad
import logging
import sys


def load_cell_predictions(hdf5_path: str | Path) -> Tuple[np.ndarray, np.ndarray]:
    """Read predictions and coords arrays from an HDF5 produced by cell classifier.

    Expects datasets named 'predictions' and 'coords'. Returns (predictions, coords).
    """
    with h5py.File(str(hdf5_path), "r") as f:
        if "predictions" not in f or "coords" not in f:
            raise KeyError("HDF5 must contain 'predictions' and 'coords' datasets")
        preds = f["predictions"][()]
        coords = f["coords"][()]
    return preds, coords


def validate_coordinate_compatibility(
    h5ad_coords: np.ndarray, 
    hdf5_coords: Optional[np.ndarray] = None,
    tissue_coords: Optional[np.ndarray] = None,
    adata: Optional[ad.AnnData] = None,
    tolerance_factor: float = 10.0
) -> Dict[str, any]:
    """Validate coordinate system compatibility between files.
    
    Args:
        h5ad_coords: Coordinates from H5AD file [N, 2] in [x, y] format
        hdf5_coords: Optional coordinates from HDF5 file [M, 2] in [x, y] format  
        tissue_coords: Optional coordinates from tissue TSV [P, 2] in [x, y] format
        adata: Optional AnnData object for metadata extraction
        tolerance_factor: Maximum acceptable ratio difference between coordinate ranges
    
    Returns:
        Dict with validation results and warnings
    """
    logger = logging.getLogger(__name__)
    results = {"compatible": True, "warnings": [], "info": {}}
    
    # Extract coordinate ranges
    h5ad_x_range = (h5ad_coords[:, 0].min(), h5ad_coords[:, 0].max())
    h5ad_y_range = (h5ad_coords[:, 1].min(), h5ad_coords[:, 1].max())
    h5ad_x_span = h5ad_x_range[1] - h5ad_x_range[0]
    h5ad_y_span = h5ad_y_range[1] - h5ad_y_range[0]
    
    results["info"]["h5ad_ranges"] = {"x": h5ad_x_range, "y": h5ad_y_range, "x_span": h5ad_x_span, "y_span": h5ad_y_span}
    
    # Check H5AD metadata for scale factors
    if adata is not None:
        if "spatial" in adata.uns:
            spatial_info = adata.uns["spatial"]
            if isinstance(spatial_info, dict) and len(spatial_info) > 0:
                sample_key = list(spatial_info.keys())[0]
                sample_data = spatial_info[sample_key]
                if "scalefactors" in sample_data:
                    scale_factors = sample_data["scalefactors"]
                    results["info"]["scale_factors"] = scale_factors
                    logger.info(f"H5AD scale factors: {scale_factors}")
                if "metadata" in sample_data:
                    metadata = sample_data["metadata"]
                    if "source_image_height" in metadata and "source_image_width" in metadata:
                        source_dims = (metadata["source_image_width"], metadata["source_image_height"])
                        results["info"]["source_image_dims"] = source_dims
                        logger.info(f"Source image dimensions: {source_dims[0]}x{source_dims[1]}")
        
        if "resolution" in adata.uns:
            results["info"]["resolution"] = adata.uns["resolution"]
            logger.info(f"H5AD resolution: {adata.uns['resolution']}")
    
    # Validate against HDF5 coordinates
    if hdf5_coords is not None:
        hdf5_x_range = (hdf5_coords[:, 0].min(), hdf5_coords[:, 0].max())
        hdf5_y_range = (hdf5_coords[:, 1].min(), hdf5_coords[:, 1].max())
        hdf5_x_span = hdf5_x_range[1] - hdf5_x_range[0]
        hdf5_y_span = hdf5_y_range[1] - hdf5_y_range[0]
        
        results["info"]["hdf5_ranges"] = {"x": hdf5_x_range, "y": hdf5_y_range, "x_span": hdf5_x_span, "y_span": hdf5_y_span}
        
        # Check span ratios
        x_ratio = hdf5_x_span / h5ad_x_span if h5ad_x_span > 0 else float('inf')
        y_ratio = hdf5_y_span / h5ad_y_span if h5ad_y_span > 0 else float('inf')
        
        results["info"]["hdf5_h5ad_ratios"] = {"x": x_ratio, "y": y_ratio}
        
        if x_ratio > tolerance_factor or x_ratio < 1/tolerance_factor:
            warning = f"Large X coordinate scale difference: HDF5/H5AD ratio = {x_ratio:.2f}"
            results["warnings"].append(warning)
            logger.warning(warning)
            
        if y_ratio > tolerance_factor or y_ratio < 1/tolerance_factor:
            warning = f"Large Y coordinate scale difference: HDF5/H5AD ratio = {y_ratio:.2f}"
            results["warnings"].append(warning)  
            logger.warning(warning)
            
        if abs(x_ratio - y_ratio) > 0.5:  # Different aspect ratio scaling
            warning = f"Inconsistent aspect ratio scaling: X ratio={x_ratio:.2f}, Y ratio={y_ratio:.2f}"
            results["warnings"].append(warning)
            logger.warning(warning)
            
        logger.info(f"HDF5 vs H5AD coordinate spans - X: {hdf5_x_span} vs {h5ad_x_span} (ratio: {x_ratio:.2f}), Y: {hdf5_y_span} vs {h5ad_y_span} (ratio: {y_ratio:.2f})")
    
    # Validate against tissue coordinates
    if tissue_coords is not None:
        tissue_x_range = (tissue_coords[:, 0].min(), tissue_coords[:, 0].max())
        tissue_y_range = (tissue_coords[:, 1].min(), tissue_coords[:, 1].max())
        tissue_x_span = tissue_x_range[1] - tissue_x_range[0]
        tissue_y_span = tissue_y_range[1] - tissue_y_range[0]
        
        results["info"]["tissue_ranges"] = {"x": tissue_x_range, "y": tissue_y_range, "x_span": tissue_x_span, "y_span": tissue_y_span}
        
        # Check span ratios
        x_ratio = tissue_x_span / h5ad_x_span if h5ad_x_span > 0 else float('inf')
        y_ratio = tissue_y_span / h5ad_y_span if h5ad_y_span > 0 else float('inf')
        
        results["info"]["tissue_h5ad_ratios"] = {"x": x_ratio, "y": y_ratio}
        
        if x_ratio > tolerance_factor or x_ratio < 1/tolerance_factor:
            warning = f"Large X coordinate scale difference: Tissue/H5AD ratio = {x_ratio:.2f}"
            results["warnings"].append(warning)
            logger.warning(warning)
            
        if y_ratio > tolerance_factor or y_ratio < 1/tolerance_factor:
            warning = f"Large Y coordinate scale difference: Tissue/H5AD ratio = {y_ratio:.2f}"
            results["warnings"].append(warning)
            logger.warning(warning)
            
        logger.info(f"Tissue vs H5AD coordinate spans - X: {tissue_x_span} vs {h5ad_x_span} (ratio: {x_ratio:.2f}), Y: {tissue_y_span} vs {h5ad_y_span} (ratio: {y_ratio:.2f})")
        
        # Check if HDF5 and tissue coordinates match (they should for the same sample)
        if hdf5_coords is not None:
            hdf5_x_span = results["info"]["hdf5_ranges"]["x_span"]
            hdf5_y_span = results["info"]["hdf5_ranges"]["y_span"]
            hdf5_tissue_x_ratio = abs(hdf5_x_span - tissue_x_span) / max(hdf5_x_span, tissue_x_span)
            hdf5_tissue_y_ratio = abs(hdf5_y_span - tissue_y_span) / max(hdf5_y_span, tissue_y_span)
            
            if hdf5_tissue_x_ratio > 0.01 or hdf5_tissue_y_ratio > 0.01:  # >1% difference
                warning = f"HDF5 and tissue coordinate ranges don't match well - may be different samples"
                results["warnings"].append(warning)
                logger.warning(warning)
            else:
                logger.info("HDF5 and tissue coordinates match well - likely same sample")
    
    if len(results["warnings"]) > 0:
        results["compatible"] = False
        logger.warning(f"Coordinate compatibility check found {len(results['warnings'])} issues")
    else:
        logger.info("Coordinate systems appear compatible")
        
    return results


def load_tissue_predictions(tsv_path: str | Path) -> pd.DataFrame:
    """Load tissue tile predictions with columns [x, y, class]."""
    df = pd.read_csv(tsv_path, sep="\t")
    required = {"x", "y", "class"}
    if not required.issubset(df.columns):
        raise KeyError(f"tissue_preds.tsv missing required columns {required}")
    return df


def compute_keep_mask(
    cell_coords: np.ndarray,
    tissue_df: pd.DataFrame,
    exclude_classes: Sequence[str],
    swap_axes: bool = False,
) -> np.ndarray:
    """Compute a boolean mask of cells to keep by excluding tiles of given classes.

    - cell_coords: array shape (N, 2) with [x, y] pixel coordinates
    - tissue_df: DataFrame with columns [x, y, class]
    - exclude_classes: classes to exclude (e.g., ['Chorion'])
    - swap_axes: if True, swap (max_x, max_y) estimation to handle axis-swapped cases
    """
    if len(cell_coords) == 0:
        return np.array([], dtype=bool)

    max_x = int(tissue_df["x"].max())
    max_y = int(tissue_df["y"].max())
    if swap_axes:
        max_x, max_y = max_y, max_x

    # Image dimensions from cell coords (coords are [x, y])
    W = int(cell_coords[:, 0].max()) + 1  # x dimension
    H = int(cell_coords[:, 1].max()) + 1  # y dimension
    if max_y <= 0 or max_x <= 0 or H <= 0 or W <= 0:
        return np.ones(len(cell_coords), dtype=bool)

    tile_size_w = W / (max_x + 1)
    tile_size_h = H / (max_y + 1)

    # Exclusion set of tile (x, y) pairs
    excl = set(map(tuple, tissue_df[tissue_df["class"].isin(exclude_classes)][["x", "y"]].values))
    if not excl:
        return np.ones(len(cell_coords), dtype=bool)

    # Map cell pixel coords to tile indices (coords are [x, y])
    tile_coords_x = np.floor(cell_coords[:, 0] / tile_size_w).astype(int)
    tile_coords_y = np.floor(cell_coords[:, 1] / tile_size_h).astype(int)
    tile_pairs = np.vstack((tile_coords_x, tile_coords_y)).T

    keep = np.ones(len(cell_coords), dtype=bool)
    for i, pair in enumerate(tile_pairs):
        if (int(pair[0]), int(pair[1])) in excl:
            keep[i] = False
    return keep


def filter_cell_predictions_by_tissue(
    hdf5_path: str | Path,
    tissue_tsv: str | Path,
    exclude_classes: Sequence[str],
    swap_axes: bool = False,
) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Load HDF5 and tissue TSV, compute keep mask, and return (preds_kept, coords_kept, keep_mask)."""
    preds, coords = load_cell_predictions(hdf5_path)
    tissue_df = load_tissue_predictions(tissue_tsv)
    keep = compute_keep_mask(coords, tissue_df, exclude_classes=exclude_classes, swap_axes=swap_axes)
    return preds[keep], coords[keep], keep


# ---------- ST AnnData utilities ----------

def get_spatial_coords(adata: ad.AnnData) -> np.ndarray:
    """Return Nx2 pixel coordinates for ST spots/cells in [x, y] format.

    Prefers obsm['spatial'] (Scanpy/Visium convention). If absent, tries obs columns
    'x', 'y' or 'px', 'py' or 'array_row', 'array_col' (converted to pixels via max range).
    """
    if "spatial" in adata.obsm:
        coords = np.asarray(adata.obsm["spatial"])  # shape (N, 2), [x, y] format
        if coords.shape[1] != 2:
            raise ValueError("obsm['spatial'] must be of shape (N, 2)")
        return coords
    obs = adata.obs
    for xcol, ycol in (("x", "y"), ("px", "py")):
        if xcol in obs and ycol in obs:
            return np.vstack([obs[xcol].to_numpy(), obs[ycol].to_numpy()]).T
    # Fallback: Visium grid positions (array row/col), scaled to pixel-like space
    if "array_row" in obs and "array_col" in obs:
        rows = obs["array_row"].to_numpy()
        cols = obs["array_col"].to_numpy()
        # Scale grid to pseudo-pixels: treat max row/col as bounds
        x = cols.astype(float)
        y = rows.astype(float)
        return np.vstack([x, y]).T
    raise KeyError("Cannot find spatial coordinates. Provide obsm['spatial'] or obs['x','y'].")


def assign_tiles_from_coords(
    coords: np.ndarray,
    tissue_df: pd.DataFrame,
    swap_axes: bool = False,
) -> np.ndarray:
    """Map pixel coords to tile (x, y) indices defined by tissue_preds.tsv grid.
    
    coords should be in [x, y] format to match tissue_df columns.
    Returns an array of shape (N, 2) of int tile indices [tile_x, tile_y].
    """
    if len(coords) == 0:
        return np.zeros((0, 2), dtype=int)
    max_x = int(tissue_df["x"].max())
    max_y = int(tissue_df["y"].max())
    if swap_axes:
        max_x, max_y = max_y, max_x
    
    # coords are [x, y], so coords[:, 0] is x, coords[:, 1] is y
    W = int(coords[:, 0].max()) + 1  # x dimension
    H = int(coords[:, 1].max()) + 1  # y dimension
    tile_size_w = W / (max_x + 1) if max_x >= 0 else 1.0
    tile_size_h = H / (max_y + 1) if max_y >= 0 else 1.0
    tile_x = np.floor(coords[:, 0] / tile_size_w).astype(int)
    tile_y = np.floor(coords[:, 1] / tile_size_h).astype(int)
    return np.vstack([tile_x, tile_y]).T


def annotate_tissue_to_spots(
    adata: ad.AnnData,
    tissue_df: pd.DataFrame,
    *,
    swap_axes: bool = False,
    column_name: str = "tissue_type",
    fill_value: str = "Unlabelled",
) -> ad.AnnData:
    """Add a column in obs with tissue class by mapping each spot/cell to tile (x,y).

    The mapping uses obsm['spatial'] pixel coords when available.
    """
    coords = get_spatial_coords(adata)
    tiles = assign_tiles_from_coords(coords, tissue_df, swap_axes=swap_axes)
    # Build lookup from (x,y) -> class
    lut: Dict[Tuple[int, int], str] = {
        (int(row.x), int(row.y)): str(row["class"]) for _, row in tissue_df.iterrows()
    }
    classes = []
    for tx, ty in tiles:
        classes.append(lut.get((int(tx), int(ty)), fill_value))
    adata = adata.copy()
    adata.obs[column_name] = pd.Categorical(classes)
    return adata


def subset_st_by_tissue(
    adata: ad.AnnData,
    tissue_df: pd.DataFrame,
    include: Optional[Sequence[str]] = None,
    exclude: Optional[Sequence[str]] = None,
    swap_axes: bool = False,
    column_name: str = "tissue_type",
) -> ad.AnnData:
    """Return AnnData subset where mapped tissue_type is included/excluded.

    If include is provided, only those classes are kept. If exclude is provided, those are dropped.
    """
    adata_annot = annotate_tissue_to_spots(adata, tissue_df, swap_axes=swap_axes, column_name=column_name)
    series = adata_annot.obs[column_name].astype("string")
    mask = pd.Series(True, index=series.index)
    if include:
        mask &= series.isin(pd.Index(include).astype("string"))
    if exclude:
        mask &= ~series.isin(pd.Index(exclude).astype("string"))
    return adata_annot[mask.values].copy()


def split_st_by_tissue(
    adata: ad.AnnData,
    tissue_df: pd.DataFrame,
    swap_axes: bool = False,
    column_name: str = "tissue_type",
) -> Dict[str, ad.AnnData]:
    """Split AnnData into dict keyed by tissue class."""
    adata_annot = annotate_tissue_to_spots(adata, tissue_df, swap_axes=swap_axes, column_name=column_name)
    classes = adata_annot.obs[column_name].astype("string").fillna("NA").unique().tolist()
    return {c: adata_annot[adata_annot.obs[column_name].astype("string") == c].copy() for c in classes}


def aggregate_st_by_cell_class(
    adata: ad.AnnData,
    cell_hdf5_path: str | Path,
    *,
    normalize: bool = True,
    group_key: Optional[str] = None,
    tissue_df: Optional[pd.DataFrame] = None,
    swap_axes: bool = False,
) -> ad.AnnData:
    """Aggregate ST expression matrix by cell-class predictions from cell classifier HDF5.

    Strategy:
    - Map ST spots to tile indices (t_x, t_y) using their pixel coords.
    - Map cell predictions to tile indices using their pixel coords from HDF5.
    - Join on tile indices to assign cell classes to spots.
    - Aggregate gene expression by summing across spots with same cell class.
    - Optionally normalize by total expression per cell class.

    Returns AnnData with cell classes as obs and aggregated expression matrix.
    """
    preds, cell_coords = load_cell_predictions(cell_hdf5_path)
    # If a tissue grid is provided, use it for tile sizing; else infer from max cell tile grid via spot coords
    if tissue_df is None:
        # fabricate a grid using max tiles derived from spot coords to ensure consistent tiling
        # We set max_x,max_y by treating each spot tile as a grid from its coords
        # Using spot coords to compute tile indices relative grid of size based on max coords
        spot_coords = get_spatial_coords(adata)
        # Build a synthetic grid with max indices inferred from spot coords quantiles
        # Default to 100x100 grid if no tissue_df provided
        max_x = 100
        max_y = 100
        tissue_df = pd.DataFrame({
            "x": np.tile(np.arange(max_x + 1), max_y + 1), 
            "y": np.repeat(np.arange(max_y + 1), max_x + 1), 
            "class": "_"
        })
    
    # Assign tiles
    spot_coords = get_spatial_coords(adata)
    spot_tiles = assign_tiles_from_coords(spot_coords, tissue_df, swap_axes=swap_axes)
    cell_tiles = assign_tiles_from_coords(cell_coords, tissue_df, swap_axes=swap_axes)

    # Build mapping from tile -> cell classes in that tile
    tile_to_cell_classes: Dict[Tuple[int, int], List[int]] = {}
    for i, (tx, ty) in enumerate(cell_tiles):
        tile_to_cell_classes.setdefault((int(tx), int(ty)), []).append(int(preds[i]))

    # Assign cell classes to spots based on majority vote in each tile
    spot_cell_classes = []
    for tx, ty in spot_tiles:
        key = (int(tx), int(ty))
        cell_classes_in_tile = tile_to_cell_classes.get(key, [])
        if cell_classes_in_tile:
            # Use majority vote for cell class assignment
            from collections import Counter
            majority_class = Counter(cell_classes_in_tile).most_common(1)[0][0]
            spot_cell_classes.append(majority_class)
        else:
            # No cells mapped to this tile
            spot_cell_classes.append(-1)  # Use -1 for unmapped spots
    
    spot_cell_classes = np.array(spot_cell_classes)
    
    # Filter out unmapped spots
    mapped_mask = spot_cell_classes != -1
    mapped_spots = np.where(mapped_mask)[0]
    mapped_classes = spot_cell_classes[mapped_mask]
    
    if len(mapped_spots) == 0:
        # Return empty AnnData if no spots are mapped
        return ad.AnnData(np.zeros((0, adata.n_vars)), var=adata.var.copy())
    
    # Get unique cell classes and aggregate expression
    unique_classes = np.unique(mapped_classes)
    aggregated_X = []
    cell_class_names = []
    
    for cell_class in unique_classes:
        class_mask = mapped_classes == cell_class
        class_spots = mapped_spots[class_mask]
        
        # Sum expression across spots of this cell class
        if hasattr(adata.X, 'toarray'):  # Handle sparse matrices
            class_expression = adata.X[class_spots].sum(axis=0).A1
        else:
            class_expression = adata.X[class_spots].sum(axis=0)
        
        aggregated_X.append(class_expression)
        cell_class_names.append(f"cell_class_{cell_class}")
    
    aggregated_X = np.vstack(aggregated_X)
    
    # Create new AnnData with aggregated data
    adata_agg = ad.AnnData(
        X=aggregated_X,
        var=adata.var.copy(),
        obs=pd.DataFrame(
            index=cell_class_names,
            data={"cell_class": unique_classes, "n_spots_aggregated": [np.sum(mapped_classes == cc) for cc in unique_classes]}
        )
    )
    
    # Add grouping information if requested
    if group_key is not None:
        if group_key not in adata.obs:
            raise KeyError(f"obs lacks column '{group_key}'")
        # For aggregated data, we can't meaningfully preserve group info since we're aggregating across spots
        # Instead, we could create separate aggregations per group, but that's a different operation
        pass
    
    # Normalize if requested
    if normalize:
        # Normalize each cell class by total expression (TPM-like)
        totals = adata_agg.X.sum(axis=1, keepdims=True)
        totals[totals == 0] = 1  # Avoid division by zero
        adata_agg.X = (adata_agg.X / totals) * 1e6  # TPM normalization
    
    return adata_agg


# ---------- Cell-class HDF5 centric utilities ----------

def map_cells_to_tissue_classes(
    cell_coords: np.ndarray,
    tissue_df: pd.DataFrame,
    *,
    swap_axes: bool = False,
    fill_value: str = "Unlabelled",
) -> np.ndarray:
    """Return an array of tissue class labels per cell by tiling coords with tissue_preds grid."""
    cell_tiles = assign_tiles_from_coords(cell_coords, tissue_df, swap_axes=swap_axes)
    lut: Dict[Tuple[int, int], str] = {
        (int(row.x), int(row.y)): str(row["class"]) for _, row in tissue_df.iterrows()
    }
    classes = []
    for tx, ty in cell_tiles:
        classes.append(lut.get((int(tx), int(ty)), fill_value))
    return np.asarray(classes, dtype=object)


def split_hdf5_by_tissue(
    hdf5_path: str | Path,
    tissue_tsv: str | Path,
    *,
    swap_axes: bool = False,
) -> Dict[str, np.ndarray]:
    """Return dict mapping tissue class -> indices of cells in that tissue tile."""
    _, coords = load_cell_predictions(hdf5_path)
    tissue_df = load_tissue_predictions(tissue_tsv)
    classes = map_cells_to_tissue_classes(coords, tissue_df, swap_axes=swap_axes)
    result: Dict[str, list] = {}
    for i, c in enumerate(classes):
        result.setdefault(str(c), []).append(i)
    return {k: np.asarray(v, dtype=int) for k, v in result.items()}


def aggregate_hdf5_by_cell_class(
    hdf5_path: str | Path,
    *,
    tissue_tsv: Optional[str | Path] = None,
    group_by_tissue: bool = False,
    normalize: bool = True,
    swap_axes: bool = False,
) -> pd.DataFrame:
    """Aggregate counts/percentages of HDF5 cell classes overall or by tissue type.

    If group_by_tissue is True, requires tissue_tsv to map cells to tissue classes.
    """
    preds, coords = load_cell_predictions(hdf5_path)
    if group_by_tissue:
        if tissue_tsv is None:
            raise ValueError("tissue_tsv is required when group_by_tissue=True")
        tissue_df = load_tissue_predictions(tissue_tsv)
        classes = map_cells_to_tissue_classes(coords, tissue_df, swap_axes=swap_axes)
        df = pd.DataFrame({"group": classes, "cell_id": preds.astype(int)})
        counts = df.groupby("group")["cell_id"].value_counts().unstack(fill_value=0)
    else:
        counts = pd.Series(preds.astype(int)).value_counts().rename_axis("cell_id").to_frame("count").T
        counts.index = ["all"]
    if normalize:
        counts = counts.div(counts.sum(axis=1), axis=0) * 100.0
    return counts

