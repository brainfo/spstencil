from __future__ import annotations

import argparse
import logging
from pathlib import Path
from typing import List

from .io import load_anndata, save_anndata  # noqa: F401 (kept for future use)
from .st import (
    load_tissue_predictions,
    split_hdf5_by_tissue,
    aggregate_hdf5_by_cell_class,
    aggregate_st_by_cell_class,
    map_cells_to_tissue_classes,
    load_cell_predictions,
    subset_st_by_tissue,
    split_st_by_tissue,
    validate_coordinate_compatibility,
    get_spatial_coords,
)


def cmd_h5_subset(args: argparse.Namespace) -> None:
    import numpy as _np
    _, coords = load_cell_predictions(args.cell_hdf5)
    tissue = load_tissue_predictions(args.tissue_tsv)
    classes = map_cells_to_tissue_classes(coords, tissue, swap_axes=args.swap_axes)
    mask = _np.ones(len(classes), dtype=bool)
    if args.include:
        include_set = set(map(str, args.include))
        mask &= _np.array([c in include_set for c in classes], dtype=bool)
    if args.exclude:
        exclude_set = set(map(str, args.exclude))
        mask &= _np.array([c not in exclude_set for c in classes], dtype=bool)
    kept_idx = _np.nonzero(mask)[0]
    args.out_idx.parent.mkdir(parents=True, exist_ok=True)
    _np.save(args.out_idx, kept_idx)


def cmd_st_subset(args: argparse.Namespace) -> None:
    # Set up logging
    logging.basicConfig(level=logging.INFO, format='%(levelname)s: %(message)s')
    logger = logging.getLogger(__name__)
    
    # Load data
    adata = load_anndata(args.st_h5ad)
    tissue = load_tissue_predictions(args.tissue_tsv)
    
    # Validate coordinate compatibility
    h5ad_coords = get_spatial_coords(adata)
    tissue_coords = tissue[['x', 'y']].values
    
    validation = validate_coordinate_compatibility(
        h5ad_coords, tissue_coords=tissue_coords, adata=adata, tolerance_factor=3.0
    )
    
    if not validation["compatible"]:
        logger.warning("Coordinate compatibility issues detected:")
        for warning in validation["warnings"]:
            logger.warning(f"  - {warning}")
        logger.warning("Subset results may not be accurate due to coordinate system misalignment")
    
    # Proceed with subset
    sub = subset_st_by_tissue(
        adata,
        tissue,
        include=args.include,
        exclude=args.exclude,
        swap_axes=args.swap_axes,
    )
    save_anndata(sub, args.output)
    logger.info(f"Subset complete: {sub.n_obs} spots retained from {adata.n_obs} original spots")


def cmd_st_split(args: argparse.Namespace) -> None:
    adata = load_anndata(args.st_h5ad)
    tissue = load_tissue_predictions(args.tissue_tsv)
    parts = split_st_by_tissue(adata, tissue, swap_axes=args.swap_axes)
    args.out_dir.mkdir(parents=True, exist_ok=True)
    for label, sub in parts.items():
        safe = str(label).replace("/", "_")
        save_anndata(sub, args.out_dir / f"{safe}.h5ad")




def cmd_h5_split(args: argparse.Namespace) -> None:
    parts = split_hdf5_by_tissue(args.cell_hdf5, args.tissue_tsv, swap_axes=args.swap_axes)
    args.out_dir.mkdir(parents=True, exist_ok=True)
    import numpy as _np
    for label, idx in parts.items():
        safe = str(label).replace("/", "_")
        _np.save(args.out_dir / f"{safe}.npy", idx)


def cmd_h5_aggregate(args: argparse.Namespace) -> None:
    df = aggregate_hdf5_by_cell_class(
        args.cell_hdf5,
        tissue_tsv=args.tissue_tsv,
        group_by_tissue=args.by_tissue,
        normalize=not args.no_normalize,
        swap_axes=args.swap_axes,
    )
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        df.to_csv(args.output)
    else:
        import sys
        df.to_csv(sys.stdout)


def cmd_st_aggregate(args: argparse.Namespace) -> None:
    # Set up logging
    logging.basicConfig(level=logging.INFO, format='%(levelname)s: %(message)s')
    logger = logging.getLogger(__name__)
    
    # Load data
    adata = load_anndata(args.st_h5ad)
    tissue_df = load_tissue_predictions(args.tissue_tsv) if args.tissue_tsv else None
    
    # Validate coordinate compatibility
    h5ad_coords = get_spatial_coords(adata)
    hdf5_coords = None
    tissue_coords = None
    
    if args.cell_hdf5:
        _, hdf5_coords = load_cell_predictions(args.cell_hdf5)
    if tissue_df is not None:
        tissue_coords = tissue_df[['x', 'y']].values
    
    validation = validate_coordinate_compatibility(
        h5ad_coords, hdf5_coords, tissue_coords, adata, tolerance_factor=3.0
    )
    
    if not validation["compatible"]:
        logger.warning("Coordinate compatibility issues detected:")
        for warning in validation["warnings"]:
            logger.warning(f"  - {warning}")
        logger.warning("Results may not be accurate due to coordinate system misalignment")
    
    # Proceed with aggregation
    adata_agg = aggregate_st_by_cell_class(
        adata,
        args.cell_hdf5,
        normalize=not args.no_normalize,
        group_key=args.group_by if hasattr(args, 'group_by') and args.group_by else None,
        tissue_df=tissue_df,
        swap_axes=args.swap_axes,
    )
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        save_anndata(adata_agg, args.output)
        logger.info(f"Saved aggregated data to {args.output}")
    else:
        import sys
        print(f"Aggregated to {adata_agg.n_obs} cell classes with {adata_agg.n_vars} genes", file=sys.stderr)
        print("Use --output to save the aggregated AnnData object", file=sys.stderr)


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(prog="spstencil", description="Spatial transcriptomics tools for subsetting, splitting, and aggregating data")
    sub = p.add_subparsers(dest="command", required=True)

    p_st_aggr = sub.add_parser("st-aggregate", help="aggregate st.h5ad to cell_bin.h5ad by cell_class.hdf5")
    p_st_aggr.add_argument("st_h5ad", type=Path, help="ST AnnData .h5ad file")
    p_st_aggr.add_argument("cell_hdf5", type=Path, help="Cell classifier HDF5 with predictions and coords")
    p_st_aggr.add_argument("--tissue-tsv", type=Path, default=None, help="tissue_preds.tsv with columns x,y,class")
    p_st_aggr.add_argument("--group-by", type=str, default=None, help="Group aggregation by column in obs")
    p_st_aggr.add_argument("--no-normalize", action="store_true", help="Skip percentage normalization")
    p_st_aggr.add_argument("--swap-axes", action="store_true", help="Handle swapped tile axes for specific samples")
    p_st_aggr.add_argument("--output", type=Path, default=None, help="Output .h5ad file (default: print summary)")
    p_st_aggr.set_defaults(func=cmd_st_aggregate)

    p_st_subset = sub.add_parser("st-subset", help="Subset st.h5ad by tissue types from tissue_preds.tsv")
    p_st_subset.add_argument("st_h5ad", type=Path, help="ST AnnData .h5ad with obsm['spatial']")
    p_st_subset.add_argument("tissue_tsv", type=Path, help="tissue_preds.tsv with columns x,y,class")
    p_st_subset.add_argument("--include", nargs="*", default=None, help="Tissue classes to include")
    p_st_subset.add_argument("--exclude", nargs="*", default=None, help="Tissue classes to exclude")
    p_st_subset.add_argument("--swap-axes", action="store_true")
    p_st_subset.add_argument("--output", type=Path, required=True, help="Output subset .h5ad")
    p_st_subset.set_defaults(func=cmd_st_subset)

    p_st_split = sub.add_parser("st-split", help="Split st.h5ad into per-tissue .h5ad files")
    p_st_split.add_argument("st_h5ad", type=Path)
    p_st_split.add_argument("tissue_tsv", type=Path)
    p_st_split.add_argument("--swap-axes", action="store_true")
    p_st_split.add_argument("--out-dir", type=Path, required=True)
    p_st_split.set_defaults(func=cmd_st_split)

    p_h5_split = sub.add_parser("h5-split", help="Split cell_class.hdf5 indices by tissue_preds.tsv")
    p_h5_split.add_argument("cell_hdf5", type=Path)
    p_h5_split.add_argument("tissue_tsv", type=Path)
    p_h5_split.add_argument("--swap-axes", action="store_true")
    p_h5_split.add_argument("--out-dir", type=Path, required=True, help="Write <tissue>.npy of indices")
    p_h5_split.set_defaults(func=cmd_h5_split)

    p_h5_aggr = sub.add_parser("h5-aggregate", help="Aggregate cell_class.hdf5 counts by tissue_preds.tsv")
    p_h5_aggr.add_argument("cell_hdf5", type=Path)
    p_h5_aggr.add_argument("tissue_tsv", type=Path, help="tissue_preds.tsv with columns x,y,class")
    p_h5_aggr.add_argument("--by-tissue", action="store_true", help="Group by tissue class")
    p_h5_aggr.add_argument("--no-normalize", action="store_true")
    p_h5_aggr.add_argument("--output", type=Path, default=None)
    p_h5_aggr.add_argument("--swap-axes", action="store_true")
    p_h5_aggr.set_defaults(func=cmd_h5_aggregate)

    p_h5_subset = sub.add_parser("h5-subset", help="Subset cell_class.hdf5 by tissue types from tissue_preds.tsv")
    p_h5_subset.add_argument("cell_hdf5", type=Path, help="Cell classifier HDF5 with predictions and coords")
    p_h5_subset.add_argument("tissue_tsv", type=Path, help="tissue_preds.tsv with columns x,y,class")
    p_h5_subset.add_argument("--include", nargs="*", default=None, help="Tissue classes to include")
    p_h5_subset.add_argument("--exclude", nargs="*", default=None, help="Tissue classes to exclude")
    p_h5_subset.add_argument("--swap-axes", action="store_true", help="Handle swapped tile axes for specific samples")
    p_h5_subset.add_argument("--out-idx", type=Path, required=True, help="Output kept indices as .npy")
    p_h5_subset.set_defaults(func=cmd_h5_subset)

    return p


def main(argv: List[str] | None = None) -> None:
    parser = build_parser()
    args = parser.parse_args(argv)
    args.func(args)


if __name__ == "__main__":
    main()

