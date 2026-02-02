"""
CRS and bounds conversion helpers.

This module provides general-purpose functions to convert geographic bounds
between coordinate reference systems (CRS). It relies on rasterio/pyproj under
the hood and is safe for projected or geographic CRSs.

Typical use case
----------------
You have bounds defined in EPSG:32630 (WGS 84 / UTM zone 30N) but the raster
you work with is in EPSG:25830 (ETRS89 / UTM zone 30N). You can convert like:

>>> from mineye.GisHelpers.crs_utils import convert_bounds
>>> bounds_32630 = (xmin, ymin, xmax, ymax)
>>> bounds_25830 = convert_bounds(bounds_32630, src_crs="EPSG:32630", dst_crs="EPSG:25830")

Or, if you already have an open rasterio dataset:

>>> from mineye.GisHelpers.crs_utils import convert_bounds_to_dataset
>>> with rasterio.open(path) as src:
...     bounds_in_src_crs = convert_bounds_to_dataset(bounds_32630, "EPSG:32630", src)

Notes
-----
- The function accepts CRS in multiple forms: string (e.g., "EPSG:32630"),
  integer EPSG code (32630), or a rasterio.crs.CRS object.
- Bounds ordering (xmin, ymin, xmax, ymax) is ensured internally.
- Densification: By default, the bounds edges are densified with 21 points
  to reduce distortion during transformation (recommended for large extents
  or lat/long CRSs). You can change this via `densify_pts`.
"""
from __future__ import annotations

from typing import Iterable, Tuple, Union

import rasterio
from rasterio.crs import CRS
from rasterio.warp import transform_bounds, calculate_default_transform, reproject, Resampling

Bounds = Tuple[float, float, float, float]
CRSType = Union[str, int, CRS]


def _normalize_bounds(bounds: Iterable[float]) -> Bounds:
    """Ensure bounds are a tuple in (xmin, ymin, xmax, ymax) order.

    This function also swaps values if min/max were given in reverse order.
    """
    try:
        xmin, ymin, xmax, ymax = tuple(bounds)
    except Exception as e:
        raise TypeError("bounds must be an iterable of four numeric values (xmin, ymin, xmax, ymax)") from e

    if xmin > xmax:
        xmin, xmax = xmax, xmin
    if ymin > ymax:
        ymin, ymax = ymax, ymin
    return float(xmin), float(ymin), float(xmax), float(ymax)


def _crs_from_any(crs: CRSType) -> CRS:
    """Convert a CRS given as string/int/CRS into a rasterio CRS object."""
    try:
        return CRS.from_user_input(crs)
    except Exception as e:
        raise ValueError(f"Invalid CRS specification: {crs!r}") from e


def convert_bounds(
    bounds: Iterable[float],
    src_crs: CRSType,
    dst_crs: CRSType,
    *,
    densify_pts: int = 21,
    allow_invalid: bool = False,
) -> Bounds:
    """Convert bounds between CRSs using rasterio.warp.transform_bounds.

    Parameters
    ----------
    bounds : iterable of float
        Input bounds as (xmin, ymin, xmax, ymax) in the source CRS.
    src_crs : str | int | rasterio.crs.CRS
        Source CRS.
    dst_crs : str | int | rasterio.crs.CRS
        Destination CRS.
    densify_pts : int, default 21
        Number of points to densify each edge with during transformation
        to reduce edge distortion. Set to 0 to disable densification.
    allow_invalid : bool, default False
        If False, raises a ValueError when the transformed bounds are invalid.

    Returns
    -------
    tuple
        Converted bounds as (xmin, ymin, xmax, ymax) in the destination CRS.
    """
    src = _crs_from_any(src_crs)
    dst = _crs_from_any(dst_crs)
    xmin, ymin, xmax, ymax = _normalize_bounds(bounds)

    # Short-circuit if CRS are equivalent
    if src == dst:
        return xmin, ymin, xmax, ymax

    try:
        tb = transform_bounds(src, dst, xmin, ymin, xmax, ymax, densify_pts=densify_pts)
    except Exception as e:
        raise ValueError(
            f"Failed to transform bounds from {src.to_string()} to {dst.to_string()}: {e}"
        ) from e

    x0, y0, x1, y1 = tb
    if not allow_invalid:
        if not (isinstance(x0, (int, float)) and isinstance(x1, (int, float)) and isinstance(y0, (int, float)) and isinstance(y1, (int, float))):
            raise ValueError("Transformed bounds are not numeric.")
        if not (min(x0, x1) < max(x0, x1) and min(y0, y1) < max(y0, y1)):
            raise ValueError("Transformed bounds are invalid (no extent). Consider adjusting input or CRS.")

    # Normalize order before returning
    return _normalize_bounds((x0, y0, x1, y1))


def convert_bounds_to_dataset(
    bounds: Iterable[float],
    src_crs: CRSType,
    dataset: rasterio.io.DatasetReader,
    *,
    densify_pts: int = 21,
    allow_invalid: bool = False,
) -> Bounds:
    """Convert bounds from a source CRS to the CRS of an open rasterio dataset.

    Parameters
    ----------
    bounds : iterable of float
        Input bounds (xmin, ymin, xmax, ymax) in the source CRS.
    src_crs : str | int | rasterio.crs.CRS
        CRS of the input bounds.
    dataset : rasterio dataset
        Open rasterio dataset whose CRS the bounds should be transformed into.
    densify_pts : int
        See convert_bounds.
    allow_invalid : bool
        See convert_bounds.

    Returns
    -------
    tuple
        Converted bounds as (xmin, ymin, xmax, ymax) in the dataset CRS.
    """
    if not hasattr(dataset, "crs") or dataset.crs is None:
        raise ValueError("Dataset has no CRS; cannot transform bounds.")
    return convert_bounds(bounds, src_crs=src_crs, dst_crs=dataset.crs, densify_pts=densify_pts, allow_invalid=allow_invalid)


def reproject_geotiff(
    src_path: str,
    dst_path: str,
    dst_crs: CRSType,
    resampling: Resampling = Resampling.nearest,
) -> None:
    """Reproject a GeoTIFF file to a new coordinate reference system.

    Parameters
    ----------
    src_path : str
        Path to the source GeoTIFF file.
    dst_path : str
        Path where the reprojected GeoTIFF will be saved.
    dst_crs : str | int | rasterio.crs.CRS
        Destination CRS.
    resampling : rasterio.warp.Resampling, default Resampling.nearest
        Resampling method to use during reprojection.
    """
    dst_crs_obj = _crs_from_any(dst_crs)

    with rasterio.open(src_path) as src:
        transform, width, height = calculate_default_transform(
            src.crs, dst_crs_obj, src.width, src.height, *src.bounds
        )
        kwargs = src.meta.copy()
        kwargs.update({
            'crs': dst_crs_obj,
            'transform': transform,
            'width': width,
            'height': height
        })

        with rasterio.open(dst_path, 'w', **kwargs) as dst:
            for i in range(1, src.count + 1):
                reproject(
                    source=rasterio.band(src, i),
                    destination=rasterio.band(dst, i),
                    src_transform=src.transform,
                    src_crs=src.crs,
                    dst_transform=transform,
                    dst_crs=dst_crs_obj,
                    resampling=resampling
                )
