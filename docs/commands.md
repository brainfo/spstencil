# Commands

The `spstencil` CLI provides end-to-end utilities for spatial transcriptomics.

Show all commands:

```bash
spstencil --help
```

## st-aggregate

Aggregate an ST `.h5ad` to per-cell-class expression using a classifier HDF5.

```bash
spstencil st-aggregate ST.h5ad CELL.hdf5 \
  --tissue-tsv tissue_preds.tsv \
  --group-by patient_id \
  --no-normalize \
  --swap-axes \
  --output cell_bin.h5ad
```

- **st_h5ad**: Spatial AnnData with `obsm['spatial']`
- **cell_hdf5**: Classifier HDF5 with predictions and coordinates
- **--tissue-tsv**: Optional `x,y,class` TSV to filter/align spots
- **--group-by**: Group aggregation by `adata.obs` column
- **--no-normalize**: Keep counts instead of percentages
- **--swap-axes**: Handle samples with swapped tile axes
- **--output**: Output `.h5ad` (default: print summary)

## st-subset

Subset ST `.h5ad` by tissue classes from a TSV.

```bash
spstencil st-subset ST.h5ad tissue_preds.tsv \
  --include Epithelium Lamina \
  --exclude Artifact \
  --swap-axes \
  --output subset.h5ad
```

- **st_h5ad**: Spatial AnnData
- **tissue_tsv**: TSV with columns `x,y,class`
- **--include/--exclude**: Tissue classes to keep/drop
- **--swap-axes**: Handle swapped tile axes
- **--output**: Output subset `.h5ad`

## st-split

Split ST `.h5ad` into per-tissue `.h5ad` files.

```bash
spstencil st-split ST.h5ad tissue_preds.tsv \
  --swap-axes \
  --out-dir split/
```

- **st_h5ad**, **tissue_tsv** as above
- **--swap-axes**: Handle swapped tile axes
- **--out-dir**: Directory to write `<tissue>.h5ad`

## h5-split

Split classifier HDF5 indices by tissue TSV.

```bash
spstencil h5-split CELL.hdf5 tissue_preds.tsv \
  --swap-axes \
  --out-dir parts/
```

- **cell_hdf5**: Classifier HDF5
- **tissue_tsv**: TSV with columns `x,y,class`
- **--swap-axes**: Handle swapped tile axes
- **--out-dir**: Write `<tissue>.npy` of kept indices

## h5-aggregate

Aggregate classifier counts by tissue TSV.

```bash
spstencil h5-aggregate CELL.hdf5 tissue_preds.tsv \
  --by-tissue \
  --no-normalize \
  --swap-axes \
  --output counts.csv
```

- **cell_hdf5**: Classifier HDF5
- **tissue_tsv**: TSV with columns `x,y,class`
- **--by-tissue**: Group by tissue class
- **--no-normalize**: Keep counts instead of percentages
- **--swap-axes**: Handle swapped tile axes
- **--output**: Output CSV (default: stdout)

## h5-subset

Subset classifier HDF5 by tissue TSV; outputs indices.

```bash
spstencil h5-subset CELL.hdf5 tissue_preds.tsv \
  --include Epithelium \
  --exclude Artifact \
  --swap-axes \
  --out-idx kept.npy
```

- **cell_hdf5**, **tissue_tsv** as above
- **--include/--exclude**: Tissue classes to keep/drop
- **--swap-axes**: Handle swapped tile axes
- **--out-idx**: Output `.npy` of kept indices

## st-crop

Detect gaps along X/Y and crop spatial, cell, and tissue data in sync.

```bash
spstencil st-crop ST.h5ad \
  --bin-size 100 \
  --output ST_crop.h5ad \
  --cell-hdf5 CELL.hdf5 --cell-output CELL_crop.hdf5 \
  --tissue-tsv tissue_preds.tsv --tissue-output tissue_crop.tsv
```

- **st_h5ad**: Spatial AnnData
- **--bin-size**: Binning for gap detection
- **--output**: Output cropped ST `.h5ad`
- **--cell-hdf5/--cell-output**: Optional, crop classifier file
- **--tissue-tsv/--tissue-output**: Optional, crop tissue TSV

---

Tips:
- Run with `--help` on any subcommand for complete options.
- Coordinates must be in the same space; use `--swap-axes` for known swapped datasets.
