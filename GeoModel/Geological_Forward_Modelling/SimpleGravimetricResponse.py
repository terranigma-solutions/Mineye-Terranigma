import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import geopandas as gpd
import os
import pickle
import gempy as gp
import gempy_viewer as gpv

# ========== CONFIG ==========
BASE_DIR = os.getcwd()
data_dir = os.path.abspath(os.path.join(BASE_DIR, '..', '..', 'Data', 'Input_Data'))
geomodel_dir = os.path.abspath(os.path.join(BASE_DIR, '..', '..', 'Data', 'Output_Data'))
geophysical_dir = os.path.join(data_dir, 'Geophysical_Cleaned_Data')
forward_model_folder = os.path.abspath(os.path.join(BASE_DIR, '..', '..', 'GeoModel', 'Geological_Forward_Modelling'))

# Paths
pickle_model_path = os.path.join(forward_model_folder, 'temp_inputs', 'simple_geo_model.pkl')
gravity_data_path = os.path.join(geophysical_dir, 'cleaned_gravity_data.geojson')

# Model parameters
densities = np.array([2.67, 2.80])  # kg/m³ for different formations (host rock, plutonite)
gravity_resolution = 20  # Number of gravity measurement points per axis

# Triggers
trigger_load_gravity_data = True
trigger_forward_modeling = True
trigger_comparison_plots = True
trigger_save_results = True

print("=== Simple Gravimetric Forward Modeling ===")

# ========== LOAD PICKLED GEOLOGICAL MODEL ==========
if not os.path.exists(pickle_model_path):
    raise FileNotFoundError(f"Pickled model not found at {pickle_model_path}. Please run Simple_Model.py first with save_pickled_model=True")

with open(pickle_model_path, 'rb') as f:
    geo_model = pickle.load(f)

print(f"Loaded geological model from {pickle_model_path}")
print(f"Model extent: {geo_model.grid.regular_grid.extent}")

# ========== LOAD GRAVITY DATA ==========
if trigger_load_gravity_data and os.path.exists(gravity_data_path):
    print(f"Loading gravity data from {gravity_data_path}")


# ========== CREATE MEASUREMENT GRID ==========
# Extract model extent for creating measurement grid
extent = geo_model.grid.regular_grid.extent
min_x, max_x, min_y, max_y, min_z, max_z = extent

# Create regular measurement grid if no field data

# ========== FORWARD MODELING ==========
if trigger_forward_modeling:
    print("Computing forward gravity model...")


print("\n=== Forward Modeling Complete ===")
