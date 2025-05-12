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

# 1. load & explode your polygons
polygons = gpd.read_file(POLYGON_PATH) \
              .explode(index_parts=False) \
              .reset_index(drop=True)


# 2. find all touching pairs
adj = gpd.sjoin(
    polygons, polygons,
    how="inner",
    predicate="touches",
    lsuffix='',        # leave left column names unchanged
    rsuffix='_right'   # only suffix the right side
)


# keep each pair only once
adj = adj[adj.index < adj["index_right"]]

# 3. extract shared boundaries
contacts = []
for _, row in adj.iterrows():
    poly_l = row.geometry
    # pull the right polygon back from the original GeoDataFrame:
    poly_r = polygons.loc[row["index_right"], "geometry"]
    shared = poly_l.boundary.intersection(poly_r)
    if shared.is_empty:
        continue
    lines = [shared] if isinstance(shared, LineString) else list(shared.geoms)
    for ln in lines:
        contacts.append({
            "geometry": ln,
            "code_1": row[CODE_FIELD],
            "code_2": row[f"{CODE_FIELD}_right"]
        })

contacts_gdf = gpd.GeoDataFrame(contacts, crs=polygons.crs)


# 4. sample points along each segment
def sample_points(line, spacing):
    if line.length <= spacing:
        return [line.interpolate(0.5, normalized=True)]
    n = int(line.length // spacing)
    return [line.interpolate(spacing*i) for i in range(n+1)]

pts = []
for _, r in contacts_gdf.iterrows():
    for pt in sample_points(r.geometry, POINT_SPACING):
        pts.append((pt, r.code_1, r.code_2))

pts_gdf = gpd.GeoDataFrame(
    pts, columns=["geometry","formation_code_1","formation_code_2"], crs=contacts_gdf.crs
)
pts_gdf["X"] = pts_gdf.geometry.x
pts_gdf["Y"] = pts_gdf.geometry.y

# 5. sample elevation Z from DEM
with rasterio.open(DEM_PATH) as src:
    if pts_gdf.crs != src.crs:
        pts_gdf = pts_gdf.to_crs(src.crs)
    coords = [(x,y) for x,y in zip(pts_gdf.X, pts_gdf.Y)]
    pts_gdf["Z"] = [v[0] for v in src.sample(coords)]

# 6. export CSV
pts_gdf[["X","Y","Z","formation_code_1","formation_code_2"]].to_csv(OUTPUT_CSV, index=False)
print(f"✅ Exported {len(pts_gdf)} contact points to {OUTPUT_CSV}")
