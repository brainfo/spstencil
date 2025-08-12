#!/usr/bin/env python3
"""
convert_cyx_to_yxc.py  ───────────────────────────────────────────────

Walk through a directory, find every *.tif(f) image assumed to be in
CYX (channels–first) order, transpose it to YXC (channels–last) order,
and save the result as a new TIFF.

The output keeps only the first three dimensions and adds an `axes`
entry to the TIFF metadata so downstream tools (e.g. `happy`) recognise
the layout.

Usage
-----
$ python -m spstencil._utils.cyx2yxc /path/to/input_dir [/path/to/output_dir]

• If *output_dir* is omitted, files are written next to the originals
  with the suffix ``_yxc.tif``.
• Recurses into sub‑directories.
"""

from __future__ import annotations

import argparse
import concurrent.futures as _fut
import multiprocessing as _mp
import pathlib
import sys
from typing import Sequence

import tifffile as tfi
from .axes import get_tiff_axes, convert_cyx_to_yxc_array


###############################################################################
# Helpers
###############################################################################

def _is_tiff(p: pathlib.Path) -> bool:
    return p.suffix.lower() in {".tif", ".tiff", ".tf8"}

def _find_tiffs(root: pathlib.Path) -> Sequence[pathlib.Path]:
    """Recursively yield TIFF files under *root*."""
    return [p for p in root.rglob("*") if _is_tiff(p)]

def _convert_one(src: pathlib.Path, out_dir: pathlib.Path | None, *, force: bool = False) -> bool:
    """Load *src*, validate axes metadata is CYX (unless forced), convert to YXC, and write.

    - If metadata axes is 'YXC', the file is skipped (already correct).
    - If metadata axes is present and not 'CYX', the file is skipped unless --force.
    - If metadata axes is missing, the file is skipped unless --force.
    """

    try:
        img = tfi.imread(src)
    except Exception as e:  # pragma: no cover  (informational)
        print(f"[ERROR] Failed to read {src}: {e}", file=sys.stderr)
        return False

    if img.ndim != 3:
        print(f"[SKIP ] {src.name}: expected 3‑D CYX, got shape {img.shape}")
        return False

    axes = get_tiff_axes(src)
    if axes is None and not force:
        print(f"[SKIP ] {src.name}: no TIFF axes metadata found; use --force to convert")
        return False
    if axes is not None:
        if axes.upper() == "YXC":
            print(f"[SKIP ] {src.name}: already in YXC order")
            return False
        if axes.upper() != "CYX" and not force:
            print(f"[SKIP ] {src.name}: axes '{axes}' not supported (expect CYX); use --force to convert")
            return False

    # CYX ➜ YXC
    yxc = convert_cyx_to_yxc_array(img)

    # Determine output path
    dst_dir = out_dir or src.parent
    dst_dir.mkdir(parents=True, exist_ok=True)
    dst_path = dst_dir / f"{src.stem}_yxc.tif"

    try:
        tfi.imwrite(dst_path, yxc, metadata={"axes": "YXC"})
        print(f"[OK   ] {src.name} → {dst_path.relative_to(dst_dir)}", flush=True)
        return True
    except Exception as e:  # pragma: no cover
        print(f"[ERROR] Failed to write {dst_path}: {e}", file=sys.stderr)
        return False


def convert_cyx_dir_to_yxc(
    input_dir: pathlib.Path,
    output_dir: pathlib.Path | None = None,
    *,
    workers: int | None = None,
    force: bool = False,
) -> int:
    """Convert all TIFFs under input_dir from CYX to YXC, respecting metadata.

    Returns the number of successfully converted files. Skips files that are already
    YXC or have incompatible/missing axes unless force=True.
    """
    if not input_dir.is_dir():
        raise NotADirectoryError(f"Input path '{input_dir}' is not a directory")

    tiffs = _find_tiffs(input_dir)
    if not tiffs:
        return 0

    max_workers = workers if workers is not None else max(_mp.cpu_count() // 2, 1)
    print(f"Found {len(tiffs)} TIFF(s), starting conversion…\n")
    converted = 0
    with _fut.ThreadPoolExecutor(max_workers=max_workers) as ex:
        futures = [ex.submit(_convert_one, p, output_dir, force=force) for p in tiffs]
        for f in futures:
            try:
                if f.result():
                    converted += 1
            except Exception as e:
                print(f"[ERROR] Worker failure: {e}", file=sys.stderr)
    print("\nDone.")
    return converted


###############################################################################
# Main CLI
###############################################################################

def _parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:  # noqa: D401
    """Return configured *argparse* namespace for given *argv*."""

    p = argparse.ArgumentParser(
        description="Convert CYX TIFFs to YXC order (channels last).",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    p.add_argument("input_dir", type=pathlib.Path, help="Directory containing TIFFs")
    p.add_argument(
        "output_dir",
        nargs="?",
        type=pathlib.Path,
        default=None,
        help="Write converted files here (defaults beside source with _yxc suffix)",
    )
    p.add_argument(
        "--force",
        action="store_true",
        help="Convert regardless of TIFF axes metadata (assume input is CYX)",
    )
    p.add_argument(
        "--workers",
        "-j",
        type=int,
        default=max(_mp.cpu_count() // 2, 1),
        help="How many parallel workers to use",
    )
    return p.parse_args(argv)


def main(argv: Sequence[str] | None = None) -> None:  # noqa: D401
    args = _parse_args(argv)

    try:
        converted = convert_cyx_dir_to_yxc(
            args.input_dir,
            args.output_dir,
            workers=args.workers,
            force=args.force,
        )
        if converted == 0:
            # Exit code 0; nothing converted isn't an error
            return
    except Exception as e:
        sys.exit(str(e))


###############################################################################
# Entry‑point guard
###############################################################################

if __name__ == "__main__":  # pragma: no cover
    main()
