import geopandas as gpd
import numpy as np
import pandas as pd
import rasterio
from shapely.geometry import LineString
import os

# ========== CONFIG ==========
script_dir = os.path.dirname(os.path.abspath(__file__))
project_root = os.path.join(script_dir, '..')
data_input_dir = os.path.join(project_root, 'Data', 'Input_Data')

POLYGON_PATH = os.path.join(data_input_dir, 'Stratigraphic_Data', 'QGIS', 'plutonite_outline.gpkg')
DEM_PATH = os.path.join(data_input_dir, 'Topographic_Data', 'topo_reduced_sf0.1.tif')
POINT_SPACING = 500
OUTPUT_CSV = os.path.join(project_root, 'Data', 'Output_Data', 'Simple-Models', 'contact_points.csv')

# Simple formation settings - all points get same formation
FORMATION_NAME = "Tournaisian Plutonites"
FORMATION_ID = 34
# ============================

# Read polygons and convert to boundary lines
polygons = gpd.read_file(POLYGON_PATH).explode(index_parts=False).reset_index(drop=True)
polygons_boundary = polygons[['geometry']].copy()
polygons_boundary['geometry'] = polygons_boundary.geometry.boundary
lines_gdf = polygons_boundary.explode(index_parts=False).reset_index(drop=True)

# Function to sample points along a line
def sample_points(line: LineString, spacing: float):
    if not isinstance(line, LineString):
        return []
    length = line.length
    if length <= spacing:
        return [line.interpolate(0.5, normalized=True)]
    distances = np.arange(0, length, spacing)
    return [line.interpolate(d) for d in distances]

# Generate points along each boundary line
point_records = []
for _, row in lines_gdf.iterrows():
    for pt in sample_points(row.geometry, POINT_SPACING):
        point_records.append({
            "geometry": pt,
            "X": pt.x,
            "Y": pt.y,
        })

points_gdf = gpd.GeoDataFrame(point_records, geometry="geometry", crs=lines_gdf.crs)

# Add elevation values from DEM
with rasterio.open(DEM_PATH) as dem:
    if points_gdf.crs != dem.crs:
        points_gdf = points_gdf.to_crs(dem.crs)

    coords = [(x, y) for x, y in zip(points_gdf.geometry.x, points_gdf.geometry.y)]
    zs = [val[0] if val is not None else 0 for val in dem.sample(coords)]
    points_gdf["Z"] = zs

# Create output DataFrame in plutonic_contact_points.csv format
output_df = pd.DataFrame({
    'X': points_gdf['X'],
    'Y': points_gdf['Y'],
    'Z': points_gdf['Z'],
    'formation': FORMATION_NAME,
    'formation_id': FORMATION_ID
})

# Save to CSV
output_df.to_csv(OUTPUT_CSV, index=False)
print(f"Generated {len(output_df)} contact points in {OUTPUT_CSV}")
print(f"Formation: {FORMATION_NAME} (ID: {FORMATION_ID})")
