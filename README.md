## spstencil

Spatial transcriptomics tools for subsetting, splitting, aggregating, and cropping data.

### Install

```bash
pip install -e .
```
### API docs

[docs](https://brainfo.github.io/spstencil/)

## Commands

### st-aggregate
Aggregate spatial transcriptomics data by cell classes.
```bash
spstencil st-aggregate data.h5ad cells.hdf5 --output aggregated.h5ad
```

### st-subset
Keep only specific tissue types from spatial data.
```bash
spstencil st-subset data.h5ad tissue.tsv --include Epithelium --output subset.h5ad
```

### st-split
Split spatial data into separate files by tissue type.
```bash
spstencil st-split data.h5ad tissue.tsv --out-dir split_data/
```

### st-crop
**Automatically crop spatial data based on detected gaps.**
```bash
# Basic cropping
spstencil st-crop data.h5ad --output cropped.h5ad

# Crop all related files together using same bounds
spstencil st-crop data.h5ad --cell-hdf5 cells.hdf5 --tissue-tsv tissue.tsv
```

**How it works:**
- Detects gaps in spatial coordinates (binned by `--bin-size`, default 100)
- When gaps are found near data edges, crops at the gap to remove empty regions
- Synchronously crops cell HDF5 and tissue TSV files using same coordinate bounds
- Example: gaps `[(33, 37), (38, 43)]` → crops at x=4300 (43×100)

### h5-aggregate  
Aggregate cell class counts from HDF5 files.
```bash
spstencil h5-aggregate cells.hdf5 tissue.tsv --output counts.csv
```

### h5-subset
Extract specific tissue classes from cell HDF5.
```bash
spstencil h5-subset cells.hdf5 tissue.tsv --include Epithelium --out-idx kept.npy
```

### h5-split
Split cell indices by tissue type.
```bash
spstencil h5-split cells.hdf5 tissue.tsv --out-dir indices/
```

### tif axes assert

- Convert a directory in-place with metadata checks:
    ```bash
    python -m spstencil._utils.cyx2yxc /path/to/dir
    ```

- Force convert (assume CYX) when metadata is absent/incorrect:
    ```bash
    python -m spstencil._utils.cyx2yxc /path/to/dir /path/to/out --force
    ```

## File Formats

- **ST .h5ad**: AnnData with spatial coordinates in `obsm['spatial']`
- **Cell .hdf5**: Contains `predictions` and `coords` datasets  
- **Tissue .tsv**: Tab-separated with columns `x`, `y`, `class`

## Options

- `--bin-size`: Coordinate binning for gap detection (default: 100)
- `--swap-axes`: Handle swapped coordinate systems
- `--no-normalize`: Skip expression normalization


