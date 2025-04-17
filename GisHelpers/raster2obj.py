import numpy as np
from osgeo import gdal, osr
import os
import shutil


def get_raster_extent(dataset):
    """Get the extent of a raster in world coordinates."""
    geotransform = dataset.GetGeoTransform()
    cols = dataset.RasterXSize
    rows = dataset.RasterYSize

    # Calculate corner coordinates
    x_min = geotransform[0]
    y_max = geotransform[3]
    x_max = x_min + cols * geotransform[1]
    y_min = y_max + rows * geotransform[5]

    return (x_min, y_min, x_max, y_max)


def get_pixel_coords(x, y, geotransform):
    """Convert world coordinates to pixel coordinates."""
    x_origin = geotransform[0]
    y_origin = geotransform[3]
    pixel_width = geotransform[1]
    pixel_height = geotransform[5]

    col = int((x - x_origin) / pixel_width)
    row = int((y - y_origin) / pixel_height)

    return row, col


def reproject_extent(extent, src_srs, dst_srs):
    """Reproject an extent from one coordinate system to another."""
    # Create coordinate transformation
    transform = osr.CoordinateTransformation(src_srs, dst_srs)
    if transform is None:
        raise ValueError("Could not create coordinate transformation. Check if both coordinate systems are valid.")

    # Transform all four corners
    points = []
    for x, y in [(extent[0], extent[1]), (extent[2], extent[3])]:
        try:
            # Use Python API properly - it returns the transformed coordinates
            x_transformed, y_transformed, z_transformed = transform.TransformPoint(x, y, 0.0)
            points.append([x_transformed, y_transformed, z_transformed])
        except Exception as e:
            raise ValueError(f"Error transforming point ({x}, {y}): {str(e)}")

    # Return the new extent
    return (min(p[0] for p in points), min(p[1] for p in points),
            max(p[0] for p in points), max(p[1] for p in points))


def crop_dem_to_extent(dem_dataset, extent, raster_srs, scale_factor=1):
    """Crop DEM to match the given extent and optionally reduce resolution.
    
    Args:
        dem_dataset: GDAL dataset of the DEM
        extent: Tuple of (x_min, y_min, x_max, y_max) in world coordinates
        raster_srs: Spatial reference system of the raster
        scale_factor: Integer factor by which to reduce resolution (e.g., 2 means half resolution)
    """
    dem_geotransform = dem_dataset.GetGeoTransform()
    dem_srs = osr.SpatialReference()
    dem_srs.ImportFromWkt(dem_dataset.GetProjection())
    
    # Print coordinate systems and extents for debugging
    print(f"Raster SRS: {raster_srs.ExportToWkt()}")
    print(f"DEM SRS: {dem_srs.ExportToWkt()}")
    print(f"Original extent: {extent}")
    
    # If coordinate systems are different, reproject the extent
    if dem_srs.IsSame(raster_srs) == 0:
        print("Coordinate systems are different, attempting reprojection...")
        try:
            # Create a transformation between the coordinate systems
            transform = osr.CoordinateTransformation(raster_srs, dem_srs)
            if transform is None:
                raise ValueError("Could not create coordinate transformation")
            
            # Transform the extent coordinates
            points = []
            for x, y in [(extent[0], extent[1]), (extent[2], extent[3])]:
                # Create a point array with 3 elements
                point = [0.0, 0.0, 0.0]
                # Transform the point - note the order: result array, x, y, z
                transform.TransformPoint(point, x, y, 0.0)
                points.append(point)
                print(f"Transformed point ({x}, {y}) -> ({point[0]}, {point[1]})")
            
            # Create new extent from transformed points
            extent = (min(p[0] for p in points), min(p[1] for p in points),
                     max(p[0] for p in points), max(p[1] for p in points))
            print(f"Transformed extent: {extent}")
            print("Reprojection successful")
        except Exception as e:
            print(f"Warning: Could not reproject coordinates: {str(e)}")
            print("Attempting to use original extent...")
    
    # Get pixel coordinates for the extent
    row_start, col_start = get_pixel_coords(extent[0], extent[3], dem_geotransform)
    row_end, col_end = get_pixel_coords(extent[2], extent[1], dem_geotransform)
    
    print(f"Pixel coordinates - row: {row_start} to {row_end}, col: {col_start} to {col_end}")
    
    # Ensure we're within bounds and in the correct order
    row_start, row_end = min(row_start, row_end), max(row_start, row_end)
    col_start, col_end = min(col_start, col_end), max(col_start, col_end)
    
    # Ensure we're within the DEM bounds
    row_start = max(0, row_start)
    col_start = max(0, col_start)
    row_end = min(dem_dataset.RasterYSize, row_end)
    col_end = min(dem_dataset.RasterXSize, col_end)
    
    print(f"Adjusted pixel coordinates - row: {row_start} to {row_end}, col: {col_start} to {col_end}")
    
    # Check if we have valid dimensions
    if row_end <= row_start or col_end <= col_start:
        raise ValueError("Invalid cropping dimensions. The raster and DEM may not overlap or have different coordinate systems.")
    
    # Read the cropped area
    band = dem_dataset.GetRasterBand(1)
    cropped_data = band.ReadAsArray(col_start, row_start, 
                                   col_end - col_start, 
                                   row_end - row_start)
    
    print(f"Cropped data shape: {cropped_data.shape}")
    
    # Apply scaling if requested
    if scale_factor > 1:
        # Calculate new dimensions
        new_rows = cropped_data.shape[0] // scale_factor
        new_cols = cropped_data.shape[1] // scale_factor
        
        # Create a new array for the scaled data
        scaled_data = np.zeros((new_rows, new_cols), dtype=cropped_data.dtype)
        
        # Average blocks of pixels
        for i in range(new_rows):
            for j in range(new_cols):
                # Get the block of pixels to average
                block = cropped_data[i*scale_factor:(i+1)*scale_factor,
                                   j*scale_factor:(j+1)*scale_factor]
                # Calculate mean, ignoring NaN values
                scaled_data[i,j] = np.nanmean(block)
        
        cropped_data = scaled_data
        print(f"Scaled data shape: {cropped_data.shape}")
    
    # Create new geotransform for the cropped area
    new_geotransform = list(dem_geotransform)
    new_geotransform[0] = dem_geotransform[0] + col_start * dem_geotransform[1]  # x_min
    new_geotransform[3] = dem_geotransform[3] + row_start * dem_geotransform[5]  # y_max
    
    # Adjust geotransform for scaling
    if scale_factor > 1:
        new_geotransform[1] *= scale_factor  # pixel width
        new_geotransform[5] *= scale_factor  # pixel height
    
    return cropped_data, tuple(new_geotransform)


def read_raster(raster_path, raster_epsg=None):
    """Read a raster file and return its data and geotransform."""
    dataset = gdal.Open(raster_path)
    if dataset is None:
        raise ValueError(f"Could not open {raster_path}")

    band = dataset.GetRasterBand(1)
    data = band.ReadAsArray()
    geotransform = dataset.GetGeoTransform()
    projection = dataset.GetProjection()

    # Create spatial reference
    srs = osr.SpatialReference()
    if projection:
        srs.ImportFromWkt(projection)
    elif raster_epsg:
        srs.ImportFromEPSG(raster_epsg)
    else:
        raise ValueError("No projection defined in raster and no EPSG code provided")

    return data, geotransform, srs, dataset


def create_vertex_grid(dem_data, geotransform, base_height=0.0):
    """Create a grid of vertices from DEM data."""
    rows, cols = dem_data.shape
    vertices = []
    
    # Print debug information
    print(f"Creating vertex grid with {rows} rows and {cols} columns")
    print(f"Geotransform: {geotransform}")
    
    # Calculate the scale factors from the geotransform
    x_scale = geotransform[1]  # pixel width
    y_scale = geotransform[5]  # pixel height (negative for north-up images)
    
    # Calculate the origin point
    x_origin = geotransform[0]
    y_origin = geotransform[3]
    
    for row in range(rows):
        for col in range(cols):
            # Calculate world coordinates
            x = x_origin + col * x_scale
            y = y_origin + row * y_scale  # Note: y_scale is negative for north-up images
            z = float(dem_data[row, col])  # Ensure we get a float value
            if z <= 0.0 or np.isnan(z):
                z = base_height
            
            vertices.append((x, y, z))
    
    # Print some sample vertices for debugging
    print(f"First vertex: {vertices[0]}")
    print(f"Last vertex: {vertices[-1]}")
    print(f"Total vertices: {len(vertices)}")
    
    return vertices


def create_faces(rows, cols):
    """Create triangular faces for the mesh."""
    faces = []
    for row in range(rows - 1):
        for col in range(cols - 1):
            # Calculate vertex indices (1-based for OBJ format)
            v1 = row * cols + col + 1
            v2 = row * cols + col + 2
            v3 = (row + 1) * cols + col + 1
            v4 = (row + 1) * cols + col + 2
            
            # Create two triangles for each quad
            # First triangle
            faces.append((v1, v2, v3))
            # Second triangle
            faces.append((v2, v4, v3))
    
    print(f"Created {len(faces)} faces")
    return faces


def create_uv_coords(rows, cols):
    """Create UV coordinates for texture mapping."""
    uv_coords = []
    for row in range(rows):
        for col in range(cols):
            u = col / (cols - 1)
            v = 1 - (row / (rows - 1))  # Flip V coordinate
            uv_coords.append((u, v))
    return uv_coords


def write_mtl_file(texture_filename, output_path):
    """Write the MTL file for the texture."""
    mtl_path = os.path.splitext(output_path)[0] + '.mtl'
    with open(mtl_path, 'w') as f:
        f.write(f"newmtl material\n")
        f.write(f"map_Kd {texture_filename}\n")
    return mtl_path


def write_obj_file(vertices, faces, uv_coords, texture_path, output_path):
    """Write the OBJ file with vertices, faces, and texture coordinates."""
    # Get the base filename for the texture
    texture_filename = os.path.basename(texture_path)

    # Create the MTL file
    mtl_path = write_mtl_file(texture_filename, output_path)

    # Copy the texture to the output directory
    output_dir = os.path.dirname(output_path)
    if output_dir:
        texture_output_path = os.path.join(output_dir, texture_filename)
        shutil.copy2(texture_path, texture_output_path)

    with open(output_path, 'w') as f:
        # Write material library reference
        f.write(f"mtllib {os.path.basename(mtl_path)}\n")

        # Write vertices
        for v in vertices:
            f.write(f"v {v[0]} {v[1]} {v[2]}\n")

        # Write texture coordinates
        for uv in uv_coords:
            f.write(f"vt {uv[0]} {uv[1]}\n")

        # Write material reference
        f.write("usemtl material\n")

        # Write faces with texture coordinates
        for face in faces:
            f.write(f"f {face[0]}/{face[0]} {face[1]}/{face[1]} {face[2]}/{face[2]}\n")


def main(raster_path, dem_path, output_path, base_height=0.0, raster_epsg=None, dem_epsg=None, scale_factor=1):
    """Main function to convert raster and DEM to OBJ file."""
    # Read raster with specified EPSG
    raster_data, raster_geotransform, raster_srs, raster_dataset = read_raster(raster_path, raster_epsg)
    raster_extent = get_raster_extent(raster_dataset)
    
    # Read DEM with specified EPSG
    dem_dataset = gdal.Open(dem_path)
    if dem_dataset is None:
        raise ValueError(f"Could not open {dem_path}")
    
    # If DEM EPSG is specified, override the projection
    if dem_epsg:
        dem_srs = osr.SpatialReference()
        dem_srs.ImportFromEPSG(dem_epsg)
        dem_dataset.SetProjection(dem_srs.ExportToWkt())
    
    # Crop DEM to match raster extent with optional scaling
    dem_data, dem_geotransform = crop_dem_to_extent(dem_dataset, raster_extent, raster_srs, scale_factor)
    
    # Create vertices, faces, and UV coordinates
    vertices = create_vertex_grid(dem_data, dem_geotransform, base_height)
    faces = create_faces(dem_data.shape[0], dem_data.shape[1])
    uv_coords = create_uv_coords(dem_data.shape[0], dem_data.shape[1])
    
    # Write OBJ file with material and texture
    write_obj_file(vertices, faces, uv_coords, raster_path, output_path)
    
    print(f"OBJ file created successfully: {output_path}")
    print(f"MTL file created: {os.path.splitext(output_path)[0] + '.mtl'}")
    print(f"Texture copied to output directory")


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description='Convert raster and DEM to OBJ file')
    parser.add_argument('raster', help='Path to the PNG raster file')
    parser.add_argument('dem', help='Path to the DEM GeoTIFF file')
    parser.add_argument('output', help='Path to the output OBJ file')
    parser.add_argument('--base-height', type=float, default=0.0,
                        help='Base height for the flat plane before DEM displacement')
    parser.add_argument('--raster-epsg', type=int,
                        help='EPSG code for the raster coordinate system (e.g., 23030 for ED50/UTM zone 30N)')
    parser.add_argument('--dem-epsg', type=int,
                        help='EPSG code for the DEM coordinate system (e.g., 25829 for ETRS89/UTM zone 29N)')
    parser.add_argument('--scale-factor', type=int, default=1,
                        help='Factor by which to reduce DEM resolution (e.g., 2 means half resolution)')

    args = parser.parse_args()

    main(args.raster, args.dem, args.output, args.base_height, args.raster_epsg, args.dem_epsg, args.scale_factor)