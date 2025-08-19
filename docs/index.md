# API

Import as:

```python
import spstencil as sps
```

## Submodules

You can import everything via the top-level too, e.g., `sps.load_anndata`, `sps.subset_by_values`, `sps.aggregate_st_by_cell_class`, etc.

### IO

::: spstencil.io

- **What**: Read/write AnnData `.h5ad`.
- **Key functions**:
  - `load_anndata(path, backed=None)`
  - `save_anndata(adata, path, compression="gzip")`
- **Example**:
```python
import spstencil as sps
adata = sps.load_anndata("data.h5ad")              # top-level alias
sps.io.save_anndata(adata, "out.h5ad")             # via submodule
```

### Generic ops

::: spstencil.ops

- **What**: Generic AnnData operations on `obs`.
- **Key functions**:
  - `subset_by_values(adata, key, include_values=None, exclude_values=None)`
  - `split_by_key(adata, key) -> dict[str, AnnData]`
  - `aggregate_cell_counts(adata, cell_type_key, group_key=None, normalize=True)`
  - `save_split_by_key(adata, key, out_dir, compression='gzip') -> dict[str, Path]`
- **Example**:
```python
subset = sps.ops.subset_by_values(adata, "tissue_type", include_values=["Epithelium"])
parts  = sps.ops.split_by_key(adata, "patient_id")
counts = sps.ops.aggregate_cell_counts(adata, "cell_type", group_key="patient_id")
written = sps.ops.save_split_by_key(adata, "patient_id", "split/")
```

### Spatial tools

::: spstencil.st

- **What**: Spatial tools: coordinate handling, tissue grid mapping, aggregation, cropping, HDF5 helpers.
- **Common workflow**:
  - Load inputs:
    - `load_tissue_predictions(tissue_tsv)`  # expects columns `x,y,class`
    - `load_cell_predictions(cell_hdf5) -> (preds, coords)`
    - `get_spatial_coords(adata) -> Nx2 [x,y]`
  - Map and subset/split:
    - `annota`te_tissue_to_spots(adata, tissue_df, swap_axes=False, column_name="tissue_type")`
    - `subset_st_by_tissue(adata, tissue_df, include=..., exclude=...)`
    - `split_st_by_tissue(adata, tissue_df)`
  - Aggregate by cell class:
    - `aggregate_st_by_cell_class(adata, cell_hdf5, normalize=True, group_key=None, tissue_df=None, swap_axes=False) -> AnnData`
  - Validate coordinates:
    - `validate_coordinate_compatibility(h5ad_coords, hdf5_coords=None, tissue_coords=None, adata=None, tolerance_factor=10.0)`
  - Cropping:
    - `determine_crop_bounds({"x": xs, "y": ys}, bin_size=100) -> bounds`
    - `crop_spatial_data(adata, bounds) -> AnnData`
    - `crop_cell_data(cell_hdf5, bounds, out_hdf5) -> int`
    - `crop_tissue_data(tissue_tsv, bounds, out_tsv, bin_size=100) -> int`
  - Utilities:
    - `map_cells_to_tissue_classes(cell_coords, tissue_df)`
    - `split_hdf5_by_tissue(cell_hdf5, tissue_tsv)`
    - `aggregate_hdf5_by_cell_class(cell_hdf5, tissue_tsv=None, group_by_tissue=False, normalize=True)`
    - `is_continuous(seq)`, `missing_ranges(seq)`
- **Example**:
```python
tissue = sps.st.load_tissue_predictions("tissue_preds.tsv")
adata_t = sps.st.subset_st_by_tissue(adata, tissue, include=["Epithelium"])

preds, coords = sps.st.load_cell_predictions("cells.hdf5")
adata_agg = sps.st.aggregate_st_by_cell_class(adata, "cells.hdf5", normalize=True)

xy = sps.st.get_spatial_coords(adata)
spatial = {"x": xy[:,0], "y": xy[:,1]}
bounds = sps.st.determine_crop_bounds(spatial, bin_size=100)
adata_crop = sps.st.crop_spatial_data(adata, bounds)
```

### Utilities

::: spstencil._utils.axes

- **What**: TIFF axis handling and conversions for arrays/images.
- **Key functions**:
  - `get_tiff_axes(path) -> Optional[str]`
  - `imread_yxc(path, force=False) -> np.ndarray`  # returns YXC
  - `convert_cyx_to_yxc_array(img) -> np.ndarray`
  - `get_yxc_dims(path, on_missing='infer', channel_upper=256, min_spatial=64) -> (Y, X, C)`
- **Example**:
```python
img_yxc = sps.axes.imread_yxc("he.tif", force=False)
y, x, c = sps.axes.get_yxc_dims("he.tif", on_missing="infer")
```

::: spstencil._utils.cyx2yxc

- **What**: Batch convert directories of TIFFs from CYX to YXC.
- **Key function**:
  - `convert_cyx_dir_to_yxc(input_dir, output_dir=None, workers=None, force=False) -> bool|int`
- **Example**:
```python
from pathlib import Path
num = sps.cyx2yxc.convert_cyx_dir_to_yxc(Path("raw/"), Path("yxc/"), workers=8, force=False)
```
- **CLI (equivalent)**:
```bash
python -m spstencil._utils.cyx2yxc /path/to/dir /path/to/out --force
```

## Top-level namespace

- Run: `spstencil --help`

::: spstencil
    options:
      members:
        - __version__
        - io
        - ops
        - st
        - axes
        - cyx2yxc
        - load_anndata
        - save_anndata
        - subset_by_values
        - split_by_key
        - aggregate_cell_counts
        - save_split_by_key
        - get_spatial_coords
        - assign_tiles_from_coords
        - annotate_tissue_to_spots
        - subset_st_by_tissue
        - split_st_by_tissue
        - aggregate_st_by_cell_class
        - load_cell_predictions
        - validate_coordinate_compatibility
        - determine_crop_bounds
        - crop_spatial_data
        - crop_cell_data
        - crop_tissue_data
        - filter_cell_predictions_by_tissue
        - map_cells_to_tissue_classes
        - split_hdf5_by_tissue
        - aggregate_hdf5_by_cell_class
        - is_continuous
        - missing_ranges
