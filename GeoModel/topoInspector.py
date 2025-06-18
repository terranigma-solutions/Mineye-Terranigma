import rasterio
import numpy as np

path_to_topography = r"C:\Users\maxha\OneDrive\Desktop\topo.tif"

with rasterio.open(path_to_topography) as src:
    data = src.read(1)
    nodata = src.nodata

print('NoData:', nodata)
print('Min:', np.min(data))
print('Max:', np.max(data))

data_cleaned = np.where(data <= -100, np.nan, data)

print('Min (cleaned):', np.nanmin(data_cleaned))
print('Max (cleaned):', np.nanmax(data_cleaned))