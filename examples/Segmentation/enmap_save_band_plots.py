"""
Example: Save per-band plots for an EnMAP L2A dataset using viridis colorscale.

This script uses the utility function `save_enmap_band_plots` from
mineye.BayesianSegmentation.EnMap_feature_extraction to generate one PNG per
spectral band in an EnMAP product. Images are saved to:

  <enmap_folder>/intermediate_Data/plots_<dataset_name>/

Usage
-----
- Edit the `enmap_folder` variable below to point to your local dataset, or
  pass it via the --folder CLI argument.
- Optional CLI arguments are available for masking and display scaling.

Examples
--------
python examples/Segmentation/examples/enmap_save_band_plots.py \
  --folder "/path/to/ENMAP01-____L2A-..." \
  --apply-mask \
  --pmin 2 --pmax 98 \
  --dpi 150

Notes
-----
- The plotting uses a non-interactive Matplotlib backend (Agg), so it can run
  headless (e.g., on servers/CI) without opening windows.
- QA masks (if present) are combined to hide bad pixels when --apply-mask is set.
"""

from __future__ import annotations

import argparse
import os
import sys

from mineye.BayesianSegmentation.EnMap_feature_extraction import save_enmap_band_plots


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        description="Save per-band viridis plots for an EnMAP L2A dataset.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    # Default dataset path provided by the user/request
    parser.add_argument(
        "--folder",
        "-f",
        dest="folder_path",
        type=str,
        default="/Volumes/macminiextern/EnMap_ROI1/ENMAP01-____L2A-DT0000026661_20230712T114043Z_002_V010402_20240818T134102Z",
        help="Path to the EnMAP L2A product folder (contains SPECTRAL_IMAGE*.tif and METADATA*.xml)",
    )
    parser.add_argument(
        "--apply-mask",
        dest="apply_mask",
        action="store_true",
        help="Apply combined QA mask to hide bad pixels",
    )
    parser.add_argument(
        "--no-mask",
        dest="apply_mask",
        action="store_false",
        help="Do not apply QA mask",
    )
    parser.set_defaults(apply_mask=True)

    parser.add_argument("--pmin", type=float, default=2.0, help="Lower percentile for per-band stretch")
    parser.add_argument("--pmax", type=float, default=98.0, help="Upper percentile for per-band stretch")
    parser.add_argument("--vmin", type=float, default=None, help="Global lower bound (overrides pmin/pmax when used with --vmax)")
    parser.add_argument("--vmax", type=float, default=None, help="Global upper bound (overrides pmin/pmax when used with --vmin)")
    parser.add_argument("--dpi", type=int, default=150, help="Output image DPI")
    parser.add_argument("--figsize", type=float, nargs=2, metavar=("W", "H"), default=(6.0, 5.0), help="Figure size in inches (width height)")

    args = parser.parse_args(argv)

    folder = args.folder_path
    if not os.path.isdir(folder):
        print(f"[ERROR] Folder not found: {folder}")
        return 2

    print(f"[INFO] EnMAP folder: {folder}")

    save_enmap_band_plots(
        folder_path=folder,
        apply_mask=args.apply_mask,
        pmin=args.pmin,
        pmax=args.pmax,
        vmin=args.vmin,
        vmax=args.vmax,
        dpi=args.dpi,
        figsize=tuple(args.figsize),
    )

    print("[INFO] Done.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
