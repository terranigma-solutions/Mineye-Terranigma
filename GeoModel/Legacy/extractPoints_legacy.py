import geopandas as gpd
import numpy as np
import rasterio
from shapely.geometry import LineString

# ========== CONFIG ==========
POLYGON_PATH   = r"C:\Users\maxha\OneDrive\Desktop\clippedGeologyExport.gpkg"
CODE_FIELD     = "CODE_UNIO"
DEM_PATH       = r"C:\Users\maxha\OneDrive\Desktop\Mineye\GIS\Tharsis\TerranigmaFiles\DEMs\AOI1_DEM_reprojected_30N.tif"
POINT_SPACING  = 500
OUTPUT_CSV     = r"C:\Users\maxha\OneDrive\Desktop\contact_points.csv"
geologicalStack = [29, 40, 1, 51, 67, 76]
idToNameDict = {
    29: "Mid Devonian Siliciclastics",
    40: "Upper Devonian Siliciclastics",
    1: "Tournaisian Plutonites",
    51: "Upper Carboniferous Volcanics",
    67: "Visean Shales",
    76: "Mid Carboniferous Shales"
}
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
points_gdf['formation_code'] = points_gdf['formation_code'].astype(int)
points_gdf = points_gdf[~points_gdf['formation_code'].isin([99, 5000])].copy()



# 5. Add z values from the DEM
with rasterio.open(DEM_PATH) as dem:
    if points_gdf.crs != dem.crs:
        points_gdf = points_gdf.to_crs(dem.crs)

    coords = [(x, y) for x, y in zip(points_gdf.geometry.x, points_gdf.geometry.y)]
    zs = [val[0] if val is not None else None for val in dem.sample(coords)]
points_gdf["Z"] = zs

# 6. Find nearest neighbor with different code
def append_nearest_code(points_gdf, code_field, max_distance):
    idx = points_gdf.sindex
    nearest_codes, nearest_dists = [], []

    for label, pt in points_gdf.geometry.items():
        # get candidate *positions* from the spatial index
        candidate_positions = list(idx.intersection(pt.buffer(max_distance).bounds))
        candidate_positions = [pos for pos in candidate_positions if pos != label]
        if not candidate_positions:
            nearest_codes.append(None); nearest_dists.append(np.nan); continue

        # slice by position, not label
        cand = points_gdf.iloc[candidate_positions]
        my_code = points_gdf.at[label, code_field]
        cand = cand[cand[code_field] != my_code]
        if cand.empty:
            nearest_codes.append(None); nearest_dists.append(np.nan); continue

        # distances and pick
        dists = cand.geometry.distance(pt).values
        posmin = dists.argmin()
        nearest_codes.append(cand.iloc[posmin][code_field])
        nearest_dists.append(dists[posmin])

    points_gdf[f'nearest_{code_field}'] = nearest_codes
    points_gdf['nearest_dist']         = nearest_dists
    return points_gdf


points_gdf = append_nearest_code(points_gdf, "formation_code", POINT_SPACING)
# Filter out points with no nearest neighbor
points_gdf = points_gdf[points_gdf[f"nearest_formation_code"].notnull()].copy()

# 7. Filter out points which are closer than POINT_SPACING/2
def enforce_min_separation(points_gdf, min_dist):
    """
    Returns a new GeoDataFrame containing a subset of points_gdf
    such that no two points are closer than min_dist apart.
    Greedy algorithm: iterates in index order, keeps a point
    only if it’s ≥ min_dist from all previously kept points.
    """
    kept_idxs = []
    kept_geoms = []

    for row in points_gdf.itertuples():
        pt = row.geometry
        # check distance to all already kept points
        too_close = any(pt.distance(other) < min_dist for other in kept_geoms)
        if not too_close:
            kept_idxs.append(row.Index)
            kept_geoms.append(pt)

    return points_gdf.loc[kept_idxs].copy()

points_gdf = enforce_min_separation(points_gdf, POINT_SPACING / 2)

# 8. Define your stratigraphic stack (in order)
pos = {code: i for i, code in enumerate(geologicalStack)}

# Remove the valid_pairs filtering completely
# Instead, just filter to ensure both formations are in the stack
points_gdf = points_gdf[
    points_gdf.apply(
        lambda r: r["formation_code"] in pos and r["nearest_formation_code"] in pos,
        axis=1
    )
].copy()

def pick_younger(row):
    c1 = row["formation_code"]
    c2 = row["nearest_formation_code"]
    if c1 in pos and c2 in pos:
        # Compare positions but return the formation CODE, not the position
        return c1 if pos[c1] > pos[c2] else c2
    return None

# Apply the pick_younger function to create the formation column
points_gdf['formationId'] = points_gdf.apply(pick_younger, axis=1)
points_gdf['formationName'] = points_gdf['formationId'].map(idToNameDict)


# 8. Export to CSV
export_fields = [
    "X",
    "Y",
    "Z",
    "formationId",
    "formationName"
]

points_gdf[export_fields].to_csv(OUTPUT_CSV, index=False)
print(f"✅ Exported {len(points_gdf)} points to {OUTPUT_CSV}")
