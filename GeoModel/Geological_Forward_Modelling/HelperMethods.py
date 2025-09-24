import os
import gempy_viewer as gpv
import numpy as np
import rasterio
from skimage.transform import resize
from sklearn.cluster import KMeans
import pandas as pd

def simplify_formation_data(orientations, points, formation_id, simp_level):
    if simp_level == 0:
        return orientations, points

    orient_data = orientations[orientations['formation_id'] == formation_id]
    points_data = points[points['formation_id'] == formation_id]

    # Determine number of points to keep based on simplification level
    n_orient_keep = max(1, int(len(orient_data) * (1 - simp_level)))
    n_points_keep = max(3, int(len(points_data) * (1 - simp_level)))

    # Use KMeans to select representative orientations
    if len(orient_data) > n_orient_keep:
        coords = orient_data[['X', 'Y']].values
        kmeans = KMeans(n_clusters=n_orient_keep, random_state=42)
        clusters = kmeans.fit_predict(coords)
        orient_indices = []
        for i in range(n_orient_keep):
            cluster_points = orient_data.iloc[clusters == i]
            centroid_idx = cluster_points.index[0]  # Take first point in cluster
            orient_indices.append(centroid_idx)
        orient_data = orient_data.loc[orient_indices]

    # Select boundary points for contact points (keep convex hull-like points)
    if len(points_data) > n_points_keep:
        coords = points_data[['X', 'Y']].values
        center = coords.mean(axis=0)
        # Calculate angles from center to select points around perimeter
        angles = np.arctan2(coords[:, 1] - center[1], coords[:, 0] - center[0])
        sorted_indices = np.argsort(angles)
        # Select evenly spaced points around the perimeter
        step = len(sorted_indices) // n_points_keep
        selected_indices = sorted_indices[::step][:n_points_keep]
        points_data = points_data.iloc[selected_indices]

    return orient_data, points_data

def clean_topo_file(input_path: str, output_path: str, invalid_below: float = -100):
    # Load, clean, and save cleaned topography TIFF

    with rasterio.open(input_path) as src:
        profile = src.profile
        data = src.read(1)
        nodata = src.nodata

    data_cleaned = np.where(
        (data == nodata) | (data <= invalid_below) | np.isnan(data),1,data)

    profile.update(nodata=1)
    with rasterio.open(output_path, "w", **profile) as dst:
        dst.write(data_cleaned, 1)

def reduce_tif_resolution(input_path, output_path, scale_factor):
    # Add scale factor to output file name
    base, ext = os.path.splitext(output_path)
    output_path_sf = f"{base}_sf{scale_factor}{ext}"

    with rasterio.open(input_path) as src:
        data = src.read(1)
        profile = src.profile

        # Calculate new shape
        new_shape = (int(data.shape[0] * scale_factor), int(data.shape[1] * scale_factor))
        # Resample data
        data_resampled = resize(data, new_shape, preserve_range=True, anti_aliasing=True).astype(data.dtype)

        # Update profile
        profile.update(
            height=new_shape[0],
            width=new_shape[1],
            transform=src.transform * src.transform.scale(
                (src.width / new_shape[1]),
                (src.height / new_shape[0])
            )
        )

        with rasterio.open(output_path_sf, "w", **profile) as dst:
            dst.write(data_resampled, 1)


def push_model_to_local_server(gempy_model, le_api):
    unstructured_data = gempy_model.raw_arrays.meshes_to_subsurface()
    le_api.upload_mesh_to_new_space(
        space_name="[TEMP] Test python api",
        data=unstructured_data,
        file_name="ROI1",
        token=os.environ.get("APIKEY")
    )

def create_cross_section(geo_model, cross_section: int):
    i = 0
    while i < cross_section:
        if i % 2 == 0:
            gpv.plot_2d(geo_model, ve=6,
                        cell_number=i,
                        show_topography=True,
                        legend=True,
                        show_data=False,
                        direction="y")
        # Increment i to avoid infinite loop
        i += 1

    i = 0
    while i < cross_section:
        if i % 2 == 0:
            gpv.plot_2d(geo_model, ve=6,
                        cell_number=i,
                        show_topography=True,
                        legend=True,
                        show_data=False,
                        direction="x")
        # Increment i to avoid infinite loop
        i += 1

def drop_lithologies(points_df, orientations_df, lithologies_to_drop):
    points_df = points_df[~points_df['formation'].isin(lithologies_to_drop)]
    orientations_df = orientations_df[~orientations_df['formation'].isin(lithologies_to_drop)]
    points_df.to_csv('temp_points.csv', index=False)
    orientations_df.to_csv('temp_orientations.csv', index=False)

def plot_3d_async(geo_model):
    gpv.plot_3d(geo_model, show_lith=True, show_boundaries=True, ve=10, legend=False, show_data=True,
                show_topography=True, transformed_data=False)

def color_lithology(structural_elements):
        for element in structural_elements:
            if element.name == "Mid Carboniferous Shales":
                element.color = "#b2d9b2"  # light green
            elif element.name == "Upper Carboniferous Volcanics":
                element.color = "#ff8000"  # orange
            elif element.name == "Visean Shales":
                element.color = "#b200b2"  # purple
            elif element.name == "Tournaisian Plutonites":
                element.color = "#e37ecd"  # pink
            elif element.name == "Upper Devonian Siliciclastics":
                element.color = "#d9b280"  # light brown

def remove_boundary_artifacts(points_df, orientations_df, boundary_tolerance=500):
    """Remove points and orientations that are near mapping boundaries"""
    cleaned_points = []
    cleaned_orientations = []

    for fid in points_df['formation_id'].unique():
        formation_points = points_df[points_df['formation_id'] == fid].copy()
        formation_orientations = orientations_df[orientations_df['formation_id'] == fid].copy()

        if len(formation_points) == 0:
            cleaned_points.append(formation_points)
            cleaned_orientations.append(formation_orientations)
            continue

        # Get bounding box for both points and orientations combined
        all_x = pd.concat([formation_points['X'], formation_orientations['X']])
        all_y = pd.concat([formation_points['Y'], formation_orientations['Y']])
        x_min, x_max = all_x.min(), all_x.max()
        y_min, y_max = all_y.min(), all_y.max()

        # Simply remove all points/orientations within boundary tolerance
        points_mask = (
            (formation_points['X'] - x_min > boundary_tolerance) &
            (x_max - formation_points['X'] > boundary_tolerance) &
            (formation_points['Y'] - y_min > boundary_tolerance) &
            (y_max - formation_points['Y'] > boundary_tolerance)
        )

        orientations_mask = (
            (formation_orientations['X'] - x_min > boundary_tolerance) &
            (x_max - formation_orientations['X'] > boundary_tolerance) &
            (formation_orientations['Y'] - y_min > boundary_tolerance) &
            (y_max - formation_orientations['Y'] > boundary_tolerance)
        )

        # Keep only interior points and orientations
        formation_points_clean = formation_points[points_mask]
        formation_orientations_clean = formation_orientations[orientations_mask]

        cleaned_points.append(formation_points_clean)
        cleaned_orientations.append(formation_orientations_clean)

    result_orientations = pd.concat(cleaned_orientations, ignore_index=True) if cleaned_orientations else pd.DataFrame()
    result_points = pd.concat(cleaned_points, ignore_index=True) if cleaned_points else pd.DataFrame()

    return result_orientations, result_points
