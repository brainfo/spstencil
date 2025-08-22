from __future__ import annotations

import numpy as np
import spstencil as sps


def sine_polyline(num=400, xlim=(-500.0, 500.0)):
    xs = np.linspace(xlim[0], xlim[1], num)
    ys = np.sin(xs / 20.0)
    return np.stack([xs, ys], axis=1)


def random_points_in_bbox(n: int, bbox: tuple[float, float, float, float]) -> np.ndarray:
    xmin, ymin, xmax, ymax = bbox
    xs = np.random.uniform(xmin, xmax, size=n)
    ys = np.random.uniform(ymin, ymax, size=n)
    return np.stack([xs, ys], axis=1)


def spiral_polyline(num_turns=10, points_per_turn=500, radius_step=1.0):
    t = np.linspace(0, 2 * np.pi * num_turns, num_turns * points_per_turn)
    r = radius_step * t
    x = r * np.cos(t)
    y = r * np.sin(t)
    return np.stack([x, y], axis=1)


def random_walk_polyline(steps=10000, step_scale=1.0):
    steps_xy = np.random.normal(scale=step_scale, size=(steps, 2))
    pts = np.cumsum(steps_xy, axis=0)
    return pts


def grid_zigzag_polyline(nx=50, ny=50, cell=10.0):
    xs = []
    ys = []
    for j in range(ny):
        y = j * cell
        if j % 2 == 0:
            for i in range(nx):
                xs.append(i * cell)
                ys.append(y)
        else:
            for i in range(nx - 1, -1, -1):
                xs.append(i * cell)
                ys.append(y)
    return np.stack([np.array(xs), np.array(ys)], axis=1)


def test_kdtree_on_sine_curve():
    poly = sine_polyline(num=1000)
    bbox = (poly[:, 0].min(), poly[:, 1].min(), poly[:, 0].max(), poly[:, 1].max())
    pts = random_points_in_bbox(1000, bbox)
    res = sps.project_points_onto_polyline(pts, poly, use_kdtree=True, kdtree_k=64)
    assert res.normalized_positions.shape == (pts.shape[0],)
    assert np.all(np.isfinite(res.distances))


def test_kdtree_worst_case_spiral():
    poly = spiral_polyline(num_turns=8, points_per_turn=800, radius_step=0.5)
    # Points near center (many segments nearby)
    pts = np.random.normal(scale=5.0, size=(500, 2))
    res = sps.project_points_onto_polyline(pts, poly, use_kdtree=True, kdtree_k=128)
    assert res.projected_points.shape == (pts.shape[0], 2)
    assert np.all(res.distances >= 0)


def test_kdtree_on_random_walk():
    poly = random_walk_polyline(steps=10000, step_scale=1.0)
    xmin, ymin = poly.min(axis=0)
    xmax, ymax = poly.max(axis=0)
    pts = random_points_in_bbox(1000, (xmin, ymin, xmax, ymax))
    res = sps.project_points_onto_polyline(pts, poly, use_kdtree=True, kdtree_k=64)
    assert res.segment_index.shape == (pts.shape[0],)


def test_grid_pruning_on_zigzag():
    poly = grid_zigzag_polyline(nx=80, ny=40, cell=5.0)
    xmin, ymin = poly.min(axis=0)
    xmax, ymax = poly.max(axis=0)
    pts = random_points_in_bbox(2000, (xmin, ymin, xmax, ymax))
    res = sps.project_points_onto_polyline(
        pts, poly, use_kdtree=True, kdtree_k=64, use_grid=True, grid_cell_size=20.0
    )
    assert res.t_along_segment.shape == (pts.shape[0],)


