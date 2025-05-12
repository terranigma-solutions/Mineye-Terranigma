import geopandas as gpd
import numpy as np
import rasterio
from shapely.geometry import LineString

# ========== CONFIG ==========
POLYGON_PATH   = r"C:\Users\maxha\OneDrive\Desktop\clippedGeology.gpkg"
CODE_FIELD     = "CODE_UNIO"
DEM_PATH       = r"C:\Users\maxha\OneDrive\Desktop\Mineye\GIS\Tharsis\TerranigmaFiles\DEMs\AOI1_DEM_reprojected_30N.tif"
POINT_SPACING  = 50
OUTPUT_CSV     = r"C:\Users\maxha\OneDrive\Desktop\contact_points.csv"
# ============================

# Read the intended layer and show its columns
polygons = gpd.read_file(POLYGON_PATH).explode(index_parts=False).reset_index(drop=True)

# 2. Convert polygons to their boundary lines
polygons_boundary = polygons[[CODE_FIELD, 'geometry']].copy()
polygons_boundary['geometry'] = polygons_boundary.geometry.boundary
lines_gdf = (
    polygons_boundary
    .explode(index_parts=False)
    .reset_index(drop=True)
)

# 3. Function to sample points along a line
def sample_points(line: LineString, spacing: float):
    if not isinstance(line, LineString):
        return []
    length = line.length
    if length <= spacing:
        return [line.interpolate(0.5, normalized=True)]
    distances = np.arange(0, length, spacing)
    return [line.interpolate(d) for d in distances]

# 4. Generate points along each boundary line
point_records = []
for _, row in lines_gdf.iterrows():
    for pt in sample_points(row.geometry, POINT_SPACING):
        point_records.append({
            "geometry": pt,
            "X": pt.x,
            "Y": pt.y,
            "formation_code": row[CODE_FIELD],
        })

points_gdf = gpd.GeoDataFrame(point_records, geometry="geometry", crs=lines_gdf.crs)
points_gdf = points_gdf[~points_gdf["formation_code"].isin(["99", "5000"])].copy()


# 5. Add z values from the DEM
with rasterio.open(DEM_PATH) as dem:
    if points_gdf.crs != dem.crs:
        points_gdf = points_gdf.to_crs(dem.crs)

    coords = [(x, y) for x, y in zip(points_gdf.geometry.x, points_gdf.geometry.y)]
    zs = [val[0] if val is not None else None for val in dem.sample(coords)]
points_gdf["Z"] = zs

# 6. Find nearest neighbor with different code
def append_nearest_code(points, code_field, max_distance):
    """
    For each point in points_gdf, find the nearest neighbor within max_distance
    that has a different code_field value. Append two new columns:
      - f'nearest_{code_field}': the neighbor's code_field value (or None)
      - 'nearest_dist': the distance to that neighbor (or np.nan)
    """
    # Build spatial index
    idx = points.sindex

    nearest_code = []
    nearest_dist = []

    for i, point in enumerate(points.geometry):
        # Candidates within bounding box of buffer
        candidates = list(idx.intersection(point.buffer(max_distance).bounds))
        # Exclude self
        candidates = [j for j in candidates if j != i]
        if not candidates:
            nearest_code.append(None)
            nearest_dist.append(np.nan)
            continue

        # Filter candidates by different code
        my_code = points.at[i, code_field]
        cand = points.iloc[candidates]
        cand = cand[cand[code_field] != my_code]
        if cand.empty:
            nearest_code.append(None)
            nearest_dist.append(np.nan)
            continue

        # Compute distances and pick minimum
        dists = cand.geometry.distance(point)
        j = dists.idxmin()
        nearest_code.append(points.at[j, code_field])
        nearest_dist.append(dists[j])

    # Assign new columns
    points[f'nearest_{code_field}'] = nearest_code
    points['nearest_dist'] = nearest_dist
    return points

points_gdf = append_nearest_code(points_gdf, CODE_FIELD, POINT_SPACING)


# 7. Export to CSV
export_fields = [
    "X",
    "Y",
    "Z",
    CODE_FIELD,                     # the original polygon code at each point
    f"nearest_{CODE_FIELD}",        # the code of the nearest neighbor
    "nearest_dist"                  # the distance to that neighbor
]

points_gdf[export_fields].to_csv(OUTPUT_CSV, index=False)
print(f"✅ Exported {len(points_gdf)} points to {OUTPUT_CSV}")
