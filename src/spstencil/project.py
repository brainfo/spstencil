from __future__ import annotations

from dataclasses import dataclass
from typing import TYPE_CHECKING, Any, Callable, List, Optional, Sequence, Tuple
import os

import numpy as np
import pandas as pd
from PIL import Image

try:
    import anndata as ad  # type: ignore
except Exception:  # pragma: no cover - optional for non-AnnData usage
    ad = None  # type: ignore

if TYPE_CHECKING:  # pragma: no cover - only for type hints
    import anndata as anndata_typing


@dataclass
class ProjectionResult:
    """Result of projecting 2D points onto a polyline.

    Attributes
    - normalized_positions: Array of shape (N,) with values in [0, 1]
    - projected_points: Array of shape (N, 2) of closest points on the polyline
    - distances: Array of shape (N,) Euclidean distance from input points to projected_points
    - segment_index: Array of shape (N,) selected segment index for each point
    - t_along_segment: Array of shape (N,) param in [0,1] along the chosen segment
    """

    normalized_positions: np.ndarray
    projected_points: np.ndarray
    distances: np.ndarray
    segment_index: np.ndarray
    t_along_segment: np.ndarray


def _validate_polyline(polyline: np.ndarray) -> np.ndarray:
    arr = np.asarray(polyline, dtype=float)
    if arr.ndim != 2 or arr.shape[1] != 2:
        raise ValueError("polyline must be an array-like of shape (M, 2)")
    if arr.shape[0] < 2:
        raise ValueError("polyline must have at least two points")
    return arr


def _compute_segment_geometry(polyline: np.ndarray) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Return segment starts S, deltas D, lengths L, and cumulative length at segment start C.

    - S: (M-1, 2)
    - D: (M-1, 2)
    - L: (M-1,)
    - C: (M-1,) cumulative arclength at each segment start
    """
    S = polyline[:-1]
    E = polyline[1:]
    D = E - S
    L = np.linalg.norm(D, axis=1)
    # Remove zero-length segments if any
    nonzero = L > 0
    S = S[nonzero]
    D = D[nonzero]
    L = L[nonzero]
    if S.shape[0] == 0:
        raise ValueError("polyline contains only zero-length segments")
    C = np.concatenate([[0.0], np.cumsum(L[:-1])])
    return S, D, L, C


def project_points_onto_polyline(
    points: np.ndarray,
    polyline: np.ndarray,
    *,
    block_size: int = 10000,
    use_kdtree: bool = True,
    kdtree_k: int = 64,
    use_grid: bool = False,
    grid_cell_size: Optional[float] = None,
) -> ProjectionResult:
    """Project 2D points onto an arbitrary polyline and return normalized positions in [0,1].

    Parameters
    - points: array-like of shape (N, 2)
    - polyline: array-like of shape (M, 2). Any curve can be approximated by a dense polyline.
    - block_size: compute in blocks to cap memory; affects speed-memory tradeoff.

    Returns
    - ProjectionResult as defined above
    """
    P = np.asarray(points, dtype=float)
    if P.ndim != 2 or P.shape[1] != 2:
        raise ValueError("points must be an array-like of shape (N, 2)")
    poly = _validate_polyline(polyline)

    S, D, L, C = _compute_segment_geometry(poly)
    total_length = float(np.sum(L))
    if not np.isfinite(total_length) or total_length <= 0:
        raise ValueError("polyline total length must be positive and finite")

    num_points = P.shape[0]
    num_segments = S.shape[0]

    distances = np.empty(num_points, dtype=float)
    projected_points = np.empty((num_points, 2), dtype=float)
    segment_index = np.empty(num_points, dtype=np.int64)
    t_along_segment = np.empty(num_points, dtype=float)

    # Precompute segment squared lengths for t computation
    seg_len_sq = (L ** 2)

    # Optional KD-tree on segment midpoints
    kd = None
    mids = None
    if use_kdtree:
        mids = S + 0.5 * D
        try:
            from scipy.spatial import cKDTree  # type: ignore
            kd = cKDTree(mids)
        except Exception:
            kd = None  # fallback to numpy top-k on mids

    # Optional coarse grid over segment bounding boxes
    grid_index: Optional[dict[tuple[int, int], list[int]]] = None
    gx = gy = xmin = ymin = cell = None  # type: ignore
    if use_grid:
        seg_min = np.minimum(S, S + D)
        seg_max = np.maximum(S, S + D)
        xmin = float(np.min(seg_min[:, 0]))
        ymin = float(np.min(seg_min[:, 1]))
        xmax = float(np.max(seg_max[:, 0]))
        ymax = float(np.max(seg_max[:, 1]))
        if grid_cell_size is None:
            med_len = float(np.median(L))
            cell = max(med_len * 4.0, 1e-6)
        else:
            cell = float(grid_cell_size)
        nx = max(int(np.ceil((xmax - xmin) / cell)), 1)
        ny = max(int(np.ceil((ymax - ymin) / cell)), 1)
        gx, gy = nx, ny
        grid_index = {}
        # Map segment bounding boxes to cells
        for idx_seg in range(num_segments):
            minx, miny = seg_min[idx_seg]
            maxx, maxy = seg_max[idx_seg]
            cx0 = int(np.clip(np.floor((minx - xmin) / cell), 0, nx - 1))
            cy0 = int(np.clip(np.floor((miny - ymin) / cell), 0, ny - 1))
            cx1 = int(np.clip(np.floor((maxx - xmin) / cell), 0, nx - 1))
            cy1 = int(np.clip(np.floor((maxy - ymin) / cell), 0, ny - 1))
            for cx in range(cx0, cx1 + 1):
                for cy in range(cy0, cy1 + 1):
                    grid_index.setdefault((cx, cy), []).append(idx_seg)

    # Process in blocks to reduce peak memory
    for start in range(0, num_points, block_size):
        end = min(num_points, start + block_size)
        PB = P[start:end]  # (B, 2)

        # If neither KD-tree nor grid is used, keep the original fully vectorized path
        if not use_kdtree and not use_grid:
            SP = PB[:, None, :] - S[None, :, :]  # (B, S, 2)
            num = np.sum(SP * D[None, :, :], axis=2)
            t = num / seg_len_sq[None, :]
            t = np.clip(t, 0.0, 1.0)
            Q = S[None, :, :] + t[:, :, None] * D[None, :, :]
            diff = Q - PB[:, None, :]
            dist2 = np.sum(diff * diff, axis=2)
            idx_all = np.argmin(dist2, axis=1)
            best_t_all = t[np.arange(end - start), idx_all]
            best_Q_all = Q[np.arange(end - start), idx_all]
            best_dist_all = np.sqrt(dist2[np.arange(end - start), idx_all])
            projected_points[start:end, :] = best_Q_all
            distances[start:end] = best_dist_all
            segment_index[start:end] = idx_all
            t_along_segment[start:end] = best_t_all
            continue

        # Otherwise, for each point build a candidate set then compute exact best within it
        for i in range(end - start):
            p = PB[i]
            cand: set[int] = set()
            # Grid-based candidates (3x3 neighborhood)
            if grid_index is not None and gx is not None and gy is not None and xmin is not None and ymin is not None and cell is not None:
                cx = int(np.clip(np.floor((p[0] - xmin) / cell), 0, gx - 1))
                cy = int(np.clip(np.floor((p[1] - ymin) / cell), 0, gy - 1))
                for dx in (-1, 0, 1):
                    for dy in (-1, 0, 1):
                        key = (int(np.clip(cx + dx, 0, gx - 1)), int(np.clip(cy + dy, 0, gy - 1)))
                        if key in grid_index:
                            cand.update(grid_index[key])
            # KD-tree candidates
            if use_kdtree:
                if kd is not None:
                    k = min(int(kdtree_k), num_segments)
                    _, idxs = kd.query(p, k=k)
                    if np.ndim(idxs) == 0:
                        idxs = np.array([int(idxs)])
                    cand.update(map(int, np.asarray(idxs).tolist()))
                else:
                    # Fallback: nearest mids by Euclidean distance via argpartition
                    if mids is None:
                        mids = S + 0.5 * D
                    k = min(int(kdtree_k), num_segments)
                    d2m = np.sum((mids - p[None, :]) ** 2, axis=1)
                    idxs = np.argpartition(d2m, k - 1)[:k]
                    cand.update(map(int, np.asarray(idxs).tolist()))
            # Fallback to all segments if no candidates
            if not cand:
                cand_idx = np.arange(num_segments, dtype=int)
            else:
                cand_idx = np.fromiter(cand, dtype=int)

            # Compute exact projection to candidate segments
            SS = S[cand_idx]
            DD = D[cand_idx]
            LL = L[cand_idx]
            len_sq_sub = (LL ** 2)
            sp = p[None, :] - SS
            nump = np.sum(sp * DD, axis=1)
            tp = np.clip(nump / len_sq_sub, 0.0, 1.0)
            Qp = SS + tp[:, None] * DD
            d2 = np.sum((Qp - p[None, :]) ** 2, axis=1)
            j = int(np.argmin(d2))
            projected_points[start + i, :] = Qp[j]
            distances[start + i] = float(np.sqrt(d2[j]))
            segment_index[start + i] = int(cand_idx[j])
            t_along_segment[start + i] = float(tp[j])

    arc_positions = C[segment_index] + t_along_segment * L[segment_index]
    normalized_positions = arc_positions / total_length

    return ProjectionResult(
        normalized_positions=normalized_positions,
        projected_points=projected_points,
        distances=distances,
        segment_index=segment_index,
        t_along_segment=t_along_segment,
    )


def sample_curve(
    curve: Callable[[np.ndarray], np.ndarray],
    *,
    num_samples: int = 1000,
) -> np.ndarray:
    """Sample an arbitrary 2D parametric curve into a polyline.

    Parameters
    - curve: function mapping t in [0,1] to array-like [x, y]
    - num_samples: number of sample points along t=[0,1]
    """
    t = np.linspace(0.0, 1.0, num_samples)
    pts = np.asarray(curve(t), dtype=float)
    if pts.ndim != 2 or pts.shape[1] != 2:
        raise ValueError("curve(t) must return array of shape (num_samples, 2)")
    return pts


def extract_polyline_from_image(
    image_path: str | "os.PathLike[str]",
    *,
    threshold: int = 64,
    invert: bool = False,
    step: int = 1,
) -> np.ndarray:
    """Extract a polyline from a binary-like ROI image by tracing dark pixels.

    The image is converted to grayscale; pixels < threshold are treated as the curve.
    A single open curve is assumed; for loops, an arbitrary start is chosen and
    a traversal is performed along 8-connected neighbors.

    Parameters
    - image_path: path to an image (e.g., JPG/PNG) with a thin dark curve
    - threshold: grayscale threshold (0-255) for foreground detection
    - invert: interpret pixels >= threshold as foreground instead
    - step: subsampling step along the traced path (>=1)

    Returns
    - polyline as float array of shape (K, 2) with columns (x, y)
    """
    # Load and binarize
    im = Image.open(image_path).convert("L")
    arr = np.asarray(im)
    if invert:
        mask = arr >= threshold
    else:
        mask = arr < threshold

    coords = np.argwhere(mask)  # (y, x)
    if coords.size == 0:
        raise ValueError("No foreground pixels detected in ROI image")

    # Use 8-connectivity
    neighbor_offsets = [
        (-1, -1), (-1, 0), (-1, 1),
        (0, -1),           (0, 1),
        (1, -1),  (1, 0),  (1, 1),
    ]

    # Build set for O(1) membership
    coord_set = set(map(tuple, coords.tolist()))

    # Degree (number of neighbors) to choose endpoints
    def degree(y: int, x: int) -> int:
        deg = 0
        for dy, dx in neighbor_offsets:
            if (y + dy, x + dx) in coord_set:
                deg += 1
        return deg

    degrees = np.array([degree(int(y), int(x)) for y, x in coords], dtype=int)
    endpoints = coords[degrees == 1]

    # Start point selection: prefer endpoint; else minimal (y, x)
    if endpoints.shape[0] > 0:
        start_y, start_x = map(int, endpoints[0])
    else:
        start_y, start_x = map(int, coords[np.lexsort((coords[:, 1], coords[:, 0]))][0])

    # Greedy path tracing
    path: List[tuple[int, int]] = []
    visited: set[tuple[int, int]] = set()
    cur = (start_y, start_x)
    prev: Optional[tuple[int, int]] = None

    # Helper to choose next neighbor preferring straightness
    def choose_next(curr: tuple[int, int], prev_pt: Optional[tuple[int, int]]) -> Optional[tuple[int, int]]:
        cy, cx = curr
        candidates: List[tuple[int, int]] = []
        for dy, dx in neighbor_offsets:
            nb = (cy + dy, cx + dx)
            if nb in coord_set and nb not in visited:
                candidates.append(nb)
        if not candidates:
            return None
        if prev_pt is None:
            return candidates[0]
        # Prefer keeping direction (maximize dot product)
        vy, vx = cy - prev_pt[0], cx - prev_pt[1]
        best = None
        best_score = -1e9
        for ny, nx in candidates:
            wy, wx = ny - cy, nx - cx
            score = vy * wy + vx * wx
            if score > best_score:
                best_score = score
                best = (ny, nx)
        return best

    while cur not in visited:
        path.append(cur)
        visited.add(cur)
        nxt = choose_next(cur, prev)
        if nxt is None:
            break
        prev, cur = cur, nxt

    if len(path) < 2:
        raise ValueError("Failed to trace a continuous polyline from ROI image")

    # Convert to (x, y) and subsample
    path_xy = np.array([(x, y) for y, x in path], dtype=float)
    if step > 1 and path_xy.shape[0] > step:
        path_xy = path_xy[::step]

    # Remove zero-length duplicates if any
    diffs = np.linalg.norm(np.diff(path_xy, axis=0), axis=1)
    keep = np.concatenate([[True], diffs > 0])
    path_xy = path_xy[keep]
    if path_xy.shape[0] < 2:
        raise ValueError("Polyline after filtering is too short")
    return path_xy


def unroll_adata_along_polyline(
    adata: "anndata_typing.AnnData",
    polyline: np.ndarray,
    *,
    coord_keys: Tuple[str, str] = ("imagecol", "imagerow"),
    features: Optional[Sequence[str]] = None,
    layer: Optional[str] = None,
    n_bins: int = 50,
    agg: str = "mean",
    distance_max: Optional[float] = None,
    block_size: int = 10000,
) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """Project spots to a line and aggregate expression along normalized 0-1.

    Parameters
    - adata: AnnData with coordinates in .obs[coord_keys]
    - polyline: (M,2) array of points describing the line/curve in the same coordinate system
    - coord_keys: column names in obs for x, y
    - features: list of gene names to aggregate. If None, all variables will be used (may be large)
    - layer: AnnData layer to read expression from. If None, uses .X
    - n_bins: number of bins along [0, 1]
    - agg: one of {'mean', 'sum', 'median'}
    - distance_max: if provided, exclude spots farther than this distance from the line
    - block_size: passed to projector for memory-performance tuning

    Returns
    - (agg_df, mapping_df)
      - agg_df: rows=bins, columns: ['s_left','s_right','s_center', <features...>]
      - mapping_df: per-spot DataFrame with ['spot', 's', 'distance', 'bin']
    """
    if ad is None:
        raise ImportError("anndata is required for unroll_adata_along_polyline")

    for key in coord_keys:
        if key not in adata.obs:
            raise KeyError(f"obs lacks column '{key}'")

    coords = adata.obs.loc[:, list(coord_keys)].astype(float).to_numpy()
    proj = project_points_onto_polyline(coords, polyline, block_size=block_size)

    s = proj.normalized_positions
    d = proj.distances

    mapping_df = pd.DataFrame(
        {
            "spot": adata.obs_names.to_list(),
            "s": s,
            "distance": d,
        }
    )

    if distance_max is not None:
        keep_mask = mapping_df["distance"].to_numpy() <= float(distance_max)
        mapping_df = mapping_df.loc[keep_mask].reset_index(drop=True)
        adata_sub = adata[keep_mask].copy()
    else:
        adata_sub = adata

    # Compute bins
    if n_bins < 1:
        raise ValueError("n_bins must be >= 1")
    edges = np.linspace(0.0, 1.0, n_bins + 1)
    mids = 0.5 * (edges[:-1] + edges[1:])

    # Assign bins (right-closed, last bin includes 1.0)
    # Use searchsorted for performance consistency
    s_vals = mapping_df["s"].to_numpy()
    # For s==1.0 ensure it falls into last bin
    s_clipped = np.clip(s_vals, 0.0, np.nextafter(1.0, 0.0))
    bin_idx = np.searchsorted(edges, s_clipped, side="right") - 1
    bin_idx = np.clip(bin_idx, 0, n_bins - 1)
    mapping_df["bin"] = bin_idx

    # Prepare expression matrix for selected features and layer
    if features is None:
        feature_names: List[str] = adata_sub.var_names.to_list()
        if layer is None:
            X = adata_sub.X
        else:
            if layer not in adata_sub.layers:
                raise KeyError(f"Layer '{layer}' not found in adata.layers")
            X = adata_sub.layers[layer]
    else:
        feature_names = [f for f in features if f in adata_sub.var_names]
        if len(feature_names) == 0:
            raise ValueError("No requested features found in adata.var_names")
        if layer is None:
            X = adata_sub[:, feature_names].X
        else:
            if layer not in adata_sub.layers:
                raise KeyError(f"Layer '{layer}' not found in adata.layers")
            X = adata_sub[:, feature_names].layers[layer]

    # Convert to dense 2D numpy array for aggregation across few features
    if hasattr(X, "toarray"):
        X = X.toarray()
    else:
        X = np.asarray(X)
    if X.ndim != 2:
        X = np.atleast_2d(X)

    # Aggregation per bin
    agg_df = pd.DataFrame({
        "s_left": edges[:-1],
        "s_right": edges[1:],
        "s_center": mids,
    })

    # Initialize result storage
    values_by_bin = np.zeros((n_bins, X.shape[1]), dtype=float)
    counts_by_bin = np.bincount(bin_idx, minlength=n_bins).astype(float)

    # Compute aggregation
    for b in range(n_bins):
        mask = bin_idx == b
        if not np.any(mask):
            values = np.full((X.shape[1],), np.nan, dtype=float)
        else:
            XB = X[mask]
            if agg == "mean":
                values = np.nanmean(XB, axis=0)
            elif agg == "sum":
                values = np.nansum(XB, axis=0)
            elif agg == "median":
                values = np.nanmedian(XB, axis=0)
            else:
                raise ValueError("agg must be one of {'mean','sum','median'}")
        values_by_bin[b, :] = values

    # Create DataFrame columns for each feature
    for j, name in enumerate(feature_names):
        agg_df[name] = values_by_bin[:, j]

    # Attach counts for interpretability
    agg_df["n_spots"] = counts_by_bin

    return agg_df, mapping_df


__all__ = [
    "ProjectionResult",
    "project_points_onto_polyline",
    "sample_curve",
    "extract_polyline_from_image",
    "unroll_adata_along_polyline",
]


def _points_elem_to_coords(points_elem: Any) -> np.ndarray:
    """Convert a SpatialData Points element or GeoDataFrame to (N,2) coords (x,y)."""
    # Prefer method commonly present in SpatialData objects
    if hasattr(points_elem, "to_geodataframe"):
        gdf = points_elem.to_geodataframe()
    else:
        gdf = points_elem
    # GeoDataFrame-like with .geometry of shapely Points
    xs = np.asarray(gdf.geometry.x, dtype=float)
    ys = np.asarray(gdf.geometry.y, dtype=float)
    return np.stack([xs, ys], axis=1)


def _shapes_elem_to_polyline(shapes_elem: Any, *, pick: str = "longest") -> np.ndarray:
    """Convert SpatialData Shapes element or GeoDataFrame to a polyline array (M,2).

    - pick: 'longest' selects the longest LineString if multiple exist
    """
    if hasattr(shapes_elem, "to_geodataframe"):
        gdf = shapes_elem.to_geodataframe()
    else:
        gdf = shapes_elem

    # Extract candidate lines
    geoms = list(getattr(gdf, "geometry"))
    lines: List[Tuple[float, np.ndarray]] = []

    try:
        from shapely.geometry import LineString, MultiLineString  # type: ignore
    except Exception as e:  # pragma: no cover
        raise ImportError("shapely is required to parse line geometries") from e

    for geom in geoms:
        if isinstance(geom, LineString):
            coords = np.asarray(geom.coords, dtype=float)
            length = float(geom.length)
            lines.append((length, coords))
        elif hasattr(geom, "geoms") and isinstance(geom, MultiLineString):
            for sub in geom.geoms:
                coords = np.asarray(sub.coords, dtype=float)
                length = float(sub.length)
                lines.append((length, coords))

    if not lines:
        raise ValueError("No LineString geometries found in shapes element")

    # Pick longest
    lines.sort(key=lambda t: t[0], reverse=True)
    _, poly = lines[0]
    # Remove zero-length duplicates
    if poly.shape[0] < 2:
        raise ValueError("Polyline too short")
    diffs = np.linalg.norm(np.diff(poly, axis=0), axis=1)
    keep = np.concatenate([[True], diffs > 0])
    return poly[keep]


def unroll_sdata(
    sdata: Any,
    points_key: str,
    line_key: str,
    *,
    table_key: Optional[str] = None,
    features: Optional[Sequence[str]] = None,
    layer: Optional[str] = None,
    n_bins: int = 50,
    agg: str = "mean",
    distance_max: Optional[float] = None,
    block_size: int = 10000,
) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """Unroll along a line for a SpatialData-like object.

    - sdata.points[points_key]: Points element (or GeoDataFrame) with Point geometries
    - sdata.shapes[line_key]: Shapes element (or GeoDataFrame) with a LineString
    - If table_key is provided and exists in sdata.tables, aggregate features from that AnnData
      using the same spot order as the Points geometry index.
    """
    # Resolve points coords
    if not hasattr(sdata, "points") or points_key not in getattr(sdata, "points"):
        raise KeyError(f"sdata.points lacks key '{points_key}'")
    points_elem = sdata.points[points_key]
    coords = _points_elem_to_coords(points_elem)

    # Resolve polyline from shapes
    if not hasattr(sdata, "shapes") or line_key not in getattr(sdata, "shapes"):
        raise KeyError(f"sdata.shapes lacks key '{line_key}'")
    shapes_elem = sdata.shapes[line_key]
    polyline = _shapes_elem_to_polyline(shapes_elem)

    # If a table is provided, attempt to aggregate features from it
    if table_key is not None and hasattr(sdata, "tables") and table_key in getattr(sdata, "tables"):
        adata = sdata.tables[table_key].copy()
        # Store coords into obs with expected keys
        obs = adata.obs.copy()
        if coords.shape[0] != adata.n_obs:
            raise ValueError("Points count does not match table n_obs; provide a matching table")
        obs["imagecol"] = coords[:, 0]
        obs["imagerow"] = coords[:, 1]
        adata.obs = obs
        return unroll_adata_along_polyline(
            adata,
            polyline,
            coord_keys=("imagecol", "imagerow"),
            features=features,
            layer=layer,
            n_bins=n_bins,
            agg=agg,
            distance_max=distance_max,
            block_size=block_size,
        )

    # Otherwise, compute mapping and counts-only aggregation
    proj = project_points_onto_polyline(coords, polyline, block_size=block_size)
    s_vals = proj.normalized_positions
    d = proj.distances
    mapping_df = pd.DataFrame({"spot": np.arange(coords.shape[0]), "s": s_vals, "distance": d})

    if distance_max is not None:
        keep = mapping_df["distance"].to_numpy() <= float(distance_max)
        mapping_df = mapping_df.loc[keep].reset_index(drop=True)

    edges = np.linspace(0.0, 1.0, n_bins + 1)
    mids = 0.5 * (edges[:-1] + edges[1:])
    s_clipped = np.clip(mapping_df["s"].to_numpy(), 0.0, np.nextafter(1.0, 0.0))
    bin_idx = np.searchsorted(edges, s_clipped, side="right") - 1
    bin_idx = np.clip(bin_idx, 0, n_bins - 1)
    mapping_df["bin"] = bin_idx

    counts = np.bincount(bin_idx, minlength=n_bins).astype(float)
    agg_df = pd.DataFrame({
        "s_left": edges[:-1],
        "s_right": edges[1:],
        "s_center": mids,
        "n_spots": counts,
    })
    return agg_df, mapping_df


