from __future__ import annotations

from pathlib import Path
from typing import Optional
import warnings

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


def _get_samples_per_pixel(path: str | Path) -> Optional[int]:
    """Best-effort read of SamplesPerPixel tag from the first page.

    Returns an integer (e.g., 1, 3, 4) or None if unavailable.
    """
    p = Path(path)
    try:
        with tfi.TiffFile(str(p)) as tf:
            if not tf.pages:
                return None
            page = tf.pages[0]
            spp = getattr(page, "samplesperpixel", None)
            if spp is None:
                try:
                    spp = page.tags["SamplesPerPixel"].value  # type: ignore[index]
                except Exception:
                    return None
            try:
                return int(spp)
            except Exception:
                return None
    except Exception:
        return None


def get_yxc_dims(
    path: str | Path,
    *,
    on_missing: str = "infer",
    channel_upper: int = 256,
    min_spatial: int = 64,
) -> tuple[int, int, int]:
    """Return (Y, X, C) dimension sizes for a 3D TIFF.

    - Uses TIFF axes metadata when available (supports 'YXC' and 'CYX').
    - If metadata is missing or unsupported:
      - on_missing='error': raise ValueError
      - on_missing='force_cyx': assume CYX ordering
      - on_missing='infer':
        1) if TIFF SamplesPerPixel matches a dimension, use that as C
        2) else prefer a dimension <= channel_upper as C; others become Y/X
           (ideally both >= min_spatial)
    """
    p = Path(path)
    with tfi.TiffFile(str(p)) as tf:
        series = tf.series[0] if tf.series else None
        if series is None or series.shape is None:
            raise ValueError(f"Cannot determine shape for {p}")
        shape = tuple(int(x) for x in series.shape)
        axes = getattr(series, "axes", None)

    if len(shape) != 3:
        raise ValueError(f"Expected 3D TIFF, got shape {shape} for {p}")

    if axes:
        up = str(axes).upper()
        if up == "YXC":
            return shape[0], shape[1], shape[2]
        if up == "CYX":
            return shape[1], shape[2], shape[0]
        # fallthrough to missing policy for unsupported axes

    if on_missing == "error":
        raise ValueError(
            f"No usable axes metadata for {p} (axes={axes}); specify on_missing or convert first"
        )
    if on_missing == "force_cyx":
        return shape[1], shape[2], shape[0]
    if on_missing != "infer":
        raise ValueError("on_missing must be one of {'error','infer','force_cyx'}")

    # Heuristic inference
    d0, d1, d2 = shape
    warnings.warn(
        f"Inferring axes for {p}: shape={shape}, on_missing={on_missing}, "
        f"channel_upper={channel_upper}, min_spatial={min_spatial}",
        RuntimeWarning,
    )
    # 1) Prefer TIFF SamplesPerPixel if it matches a dimension
    spp = _get_samples_per_pixel(p)
    if spp is not None and spp in (d0, d1, d2):
        if spp == d0:
            return d1, d2, d0  # assume CYX
        if spp == d2:
            return d0, d1, d2  # assume YXC
        if spp == d1:
            # middle-as-channels; pick larger as X
            y, x = (d0, d2) if d2 >= d0 else (d2, d0)
            return y, x, d1

    # 2) If a boundary-small dim looks like channels
    if d0 <= channel_upper and d1 >= min_spatial and d2 >= min_spatial:
        return d1, d2, d0  # assume CYX
    if d2 <= channel_upper and d0 >= min_spatial and d1 >= min_spatial:
        return d0, d1, d2  # assume YXC
    # Generic: smallest is channels
    dims = [d0, d1, d2]
    c_idx = int(np.argmin(dims))
    c = dims[c_idx]
    spatial = [dims[i] for i in range(3) if i != c_idx]
    if not (spatial[0] >= min_spatial and spatial[1] >= min_spatial):
        # Prefer YXC if last dim is the smallest and reasonably small
        if c_idx == 2 and c <= channel_upper:
            return d0, d1, d2
        # otherwise assume CYX
        return d1, d2, d0
    # Map back to Y,X depending on which index is channels
    if c_idx == 0:  # CYX
        return d1, d2, d0
    if c_idx == 2:  # YXC
        return d0, d1, d2
    # Middle-as-channels is unusual; pick larger as X
    y, x = sorted(spatial)
    return y, x, c


