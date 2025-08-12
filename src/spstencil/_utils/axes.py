from __future__ import annotations

from pathlib import Path
from typing import Optional

import numpy as np
import tifffile as tfi


def get_tiff_axes(path: str | Path) -> Optional[str]:
    """Return axes string from first series of a TIFF, if available.

    Examples: "YXC", "CYX", "TCZYX". Returns None if not present.
    """
    p = Path(path)
    try:
        with tfi.TiffFile(str(p)) as tf:
            series = tf.series[0] if tf.series else None
            axes = getattr(series, "axes", None)
            return str(axes) if axes is not None else None
    except Exception:
        return None


def convert_cyx_to_yxc_array(img: np.ndarray) -> np.ndarray:
    """Convert a 3D CYX numpy array to YXC by transposing (1, 2, 0).

    Raises ValueError if the input is not 3D.
    """
    if img.ndim != 3:
        raise ValueError(f"Expected 3D array (CYX), got shape {img.shape}")
    return img.transpose(1, 2, 0)


def imread_yxc(path: str | Path, *, force: bool = False) -> np.ndarray:
    """Read a TIFF from disk and return a 3D YXC numpy array.

    Uses TIFF metadata to determine axes:
    - If axes == 'YXC': returns as-is
    - If axes == 'CYX': converts to YXC
    - If axes missing or different: raises unless force=True (assumes CYX)
    """
    p = Path(path)
    img = tfi.imread(str(p))
    if img.ndim != 3:
        raise ValueError(f"Expected 3D TIFF, got shape {img.shape} for {p}")
    axes = get_tiff_axes(p)
    if axes is None:
        if not force:
            raise ValueError(f"No TIFF axes metadata for {p}; pass force=True to assume CYX")
        return convert_cyx_to_yxc_array(img)
    axes_up = axes.upper()
    if axes_up == "YXC":
        return img
    if axes_up == "CYX":
        return convert_cyx_to_yxc_array(img)
    if not force:
        raise ValueError(f"Unsupported axes '{axes}' for {p}; expected CYX or YXC")
    return convert_cyx_to_yxc_array(img)


