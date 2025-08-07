from __future__ import annotations

from pathlib import Path
from typing import Optional

import anndata as ad


def load_anndata(path: str | Path, backed: Optional[str] = None) -> ad.AnnData:
    """Load an AnnData object from an .h5ad file.

    - path: path to .h5ad
    - backed: pass 'r' for backed mode, or None for in-memory
    """
    file_path = Path(path)
    if not file_path.exists():
        raise FileNotFoundError(f"No such file: {file_path}")
    return ad.read_h5ad(str(file_path), backed=backed)


def save_anndata(adata: ad.AnnData, path: str | Path, compression: str | None = "gzip") -> None:
    """Save an AnnData object to .h5ad.

    - compression: None|'gzip'|'lzf' etc.
    """
    out = Path(path)
    out.parent.mkdir(parents=True, exist_ok=True)
    adata.write_h5ad(str(out), compression=compression)

