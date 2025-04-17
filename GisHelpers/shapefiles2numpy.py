import geopandas as gpd
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path

def import_and_visualize_geojson(file_attributes):
    """
    Import multiple GeoJSON files containing point data, convert them to NumPy arrays, and visualize them together.
    
    Parameters:
    -----------
    file_attributes : dict
        Dictionary where keys are file paths and values are lists of attributes to import.
        Example: {
            "path/to/file1.geojson": ["attr1", "attr2"],
            "path/to/file2.geojson": ["attr3", "attr4"]
        }
    """
    if not file_attributes:
        return

    # Set the backend to TkAgg for PyCharm compatibility
    plt.switch_backend('TkAgg')
    
    # Create a figure and axis
    fig, ax = plt.subplots(figsize=(12, 10))
    
    # Set equal aspect ratio for proper spatial representation
    ax.set_aspect('equal')
    
    # Define different markers for each dataset
    markers = ['o', 's', '^', 'v', '<', '>', 'p', '*', 'h', 'H', 'D', 'd']
    
    # Process each file
    for idx, (file_path, attributes) in enumerate(file_attributes.items()):
        try:
            # Read the GeoJSON file
            gdf = gpd.read_file(file_path)
            
            # Print available attributes and their types
            print(f"\nAttributes in {Path(file_path).stem}:")
            for col in gdf.columns:
                if col != 'geometry':  # Skip the geometry column
                    print(f"- {col}: {gdf[col].dtype}")
            
            # Get the filename for the legend
            file_name = Path(file_path).stem
            
            # Extract all coordinates at once
            coords = np.array([(point.x, point.y) for point in gdf.geometry])
            
            # Get the first attribute for coloring
            if attributes and attributes[0] in gdf.columns:
                values = gdf[attributes[0]].values
                # Create scatter plot with color based on the first attribute
                scatter = ax.scatter(
                    coords[:, 0], 
                    coords[:, 1], 
                    c=values,
                    marker=markers[idx % len(markers)],
                    label=file_name,
                    s=10
                )
                # Add colorbar for this dataset
                cbar = plt.colorbar(scatter, ax=ax)
                cbar.set_label(attributes[0])
            else:
                # If no valid attribute, plot without color
                ax.scatter(
                    coords[:, 0], 
                    coords[:, 1], 
                    marker=markers[idx % len(markers)],
                    label=file_name,
                    s=10
                )
            
            # Print the selected attributes for this file
            print(f"\nSelected attributes for {file_name}:")
            for attr in attributes:
                if attr in gdf.columns:
                    print(f"- {attr}: {gdf[attr].dtype}")
                else:
                    print(f"- {attr}: Attribute not found")
                
        except Exception as e:
            print(f"Error processing file {file_path}: {str(e)}")
            continue
    
    # Set plot properties
    ax.set_title('Point Data Visualization')
    ax.set_xlabel('Longitude')
    ax.set_ylabel('Latitude')
    ax.legend()
    ax.grid(True)
    
    # Show the plot
    plt.tight_layout()
    plt.show(block=True)  # Use block=True to keep the plot window open

if __name__ == "__main__":
    # Example usage
    file_attributes = {
        "/Volumes/macminiextern/Tharsis/TerranigmaFiles /shapeFiles/AOI1_airborneMagnetics_EPSG23030.geojson": [
            "MAG",  # example attribute
            "POT"  # example attribute
        ],
        "/Volumes/macminiextern/Tharsis/TerranigmaFiles /shapeFiles/AOI1_GF_Gravity_cropped_EPSG23030.geojson": [
            "VALU_BOU267"
        ]
    }
    import_and_visualize_geojson(file_attributes)
