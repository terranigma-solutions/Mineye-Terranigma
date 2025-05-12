import gempy as gp
import gempy_viewer as gpv
import pandas as pd

surface_points = r"C:\Users\maxha\OneDrive\Desktop\formationinputpoints.csv"

geo_model = gp.create_geomodel(project_name = 'AOI',
                               importer_helper = gp.data.ImporterHelper(
                                   path_to_surface_points=surface_points
                               ))

gpv.plot_2d(geo_model, show=True)