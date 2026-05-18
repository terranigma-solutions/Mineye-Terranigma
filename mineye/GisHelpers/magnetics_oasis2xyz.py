from pathlib import Path

import numpy as np
import pandas as pd
import rasterio
from rasterio.transform import xy
from rasterio.crs import CRS


# ============================================================
# CONFIG
# ============================================================

# Use paths relative to this script's location
SCRIPT_DIR = Path(__file__).resolve().parent
PROJECT_ROOT = SCRIPT_DIR.parents[1]

INPUT_DIR = PROJECT_ROOT / "examples/Data/General_Input_Data/Geophysical_Raw_Data/Ternove_Magnetics"
OUTPUT_DIR = PROJECT_ROOT / "examples/Data/General_Input_Data/Geophysical_Cleaned_Data/Ternove_Magnetics"

INPUT_FILES = [
    INPUT_DIR / "merged_B1B2_mean_median_201.ers",
    INPUT_DIR / "merged_masked_B1B2_mean_median_201_up3m.ers",
]

# Correct CRS for this dataset
# (the masked dataset metadata appears broken)
FORCED_CRS = CRS.from_epsg(32634)

# Additional nodata values sometimes present in these datasets
EXTRA_NODATA_VALUES = [-999999, -9999]


# ============================================================
# EXPORT FUNCTION
# ============================================================

def export_raster_to_xyz(input_path, output_dir):
    input_path = Path(input_path)
    output_dir = Path(output_dir)

    output_dir.mkdir(parents=True, exist_ok=True)

    print(f"\nProcessing: {input_path.name}")

    with rasterio.open(input_path) as src:

        # Read raster band
        data = src.read(1)

        # Get affine transform
        transform = src.transform

        # Get nodata value from metadata
        nodata = src.nodata

        # ----------------------------------------------------
        # Build valid-data mask
        # ----------------------------------------------------

        valid = np.isfinite(data)

        if nodata is not None:
            valid &= data != nodata

        for nd in EXTRA_NODATA_VALUES:
            valid &= data != nd

        # Get valid cell indices
        rows, cols = np.where(valid)

        # Extract values
        values = data[rows, cols]

        # Convert raster indices -> real-world coordinates
        xs, ys = xy(transform, rows, cols, offset="center")

        # ----------------------------------------------------
        # Create dataframe
        # ----------------------------------------------------

        df = pd.DataFrame({
            "x": np.array(xs),
            "y": np.array(ys),
            "value": values.astype(np.float32),
        })

        # Add CRS information
        df.attrs["crs"] = FORCED_CRS.to_string()

        # ----------------------------------------------------
        # Output paths
        # ----------------------------------------------------

        stem = input_path.stem.replace(".ers", "")

        csv_path = output_dir / f"{stem}.csv"
        parquet_path = output_dir / f"{stem}.parquet"

        # ----------------------------------------------------
        # Save CSV
        # ----------------------------------------------------

        df.to_csv(csv_path, index=False)

        # ----------------------------------------------------
        # Save Parquet
        # ----------------------------------------------------

        # Requires:
        # pip install pyarrow
        df.to_parquet(parquet_path, index=False)

        # ----------------------------------------------------
        # Print summary
        # ----------------------------------------------------

        print(f"  Valid points: {len(df):,}")
        print(f"  X range: {df.x.min():.2f} - {df.x.max():.2f}")
        print(f"  Y range: {df.y.min():.2f} - {df.y.max():.2f}")
        print(f"  Value range: {df.value.min():.2f} - {df.value.max():.2f}")

        print(f"  CSV saved: {csv_path}")
        print(f"  Parquet saved: {parquet_path}")


# ============================================================
# MAIN
# ============================================================

if __name__ == "__main__":

    for f in INPUT_FILES:
        export_raster_to_xyz(f, OUTPUT_DIR)

    print("\nDone.")