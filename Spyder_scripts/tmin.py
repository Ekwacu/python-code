# -*- coding: utf-8 -*-
"""
Created on Sat Jun  3 20:16:37 2023

@author: 1
"""

import xarray as xr
import numpy as np
data_arrays = []
tanzania_lat_min = -11.75
tanzania_lat_max = -0.99
tanzania_lon_min = 29.34
tanzania_lon_max = 40.45
# Create an empty list to store the DataArray objects
data_arrays = []

# Iterate over the years
for year in range(1981, 2023):
    # Read the dataset for the current year
    file_path = f'tmin.{year}.nc'
    dataset = xr.open_dataset(file_path)
    # Extract latitude and longitude variables from the dataset
    latitudes = dataset.variables['lat'][:]
    longitudes = dataset.variables['lon'][:]

     # Find the indices of the latitude and longitude range for Tanzania
    tanzania_lat_indices = (latitudes >= tanzania_lat_min) & (latitudes <= tanzania_lat_max)
    tanzania_lon_indices = (longitudes >= tanzania_lon_min) & (longitudes <= tanzania_lon_max)

     # Extract the 'tmin' variable as a DataArray
    data_array = dataset['tmin'][:, tanzania_lat_indices, tanzania_lon_indices]
    
    # Append the DataArray to the list
    data_arrays.append(data_array)
# Concatenate the DataArray objects along the 'year' dimension
merged_data = xr.concat(data_arrays, dim='time')

# Save the merged dataset as a netCDF file
# output_file = 'tanzania_tmin_data.nc'
# merged_data.to_netcdf(output_file)

# print("Data saved successfully.")
# Concatenate the DataArray objects along the 'year' dimension
# merged_data = xr.concat(data_arrays, dim='time')
# # merged_data = pd.concat(frames, keys=["x", "y", "z"])
# # Save the merged dataset as a netCDF file
# output_file = 'merged_data.nc'
# merged_data.to_netcdf(output_file)

# # Print the dimensions of the merged dataset
# print(f"Merged data dimensions: {merged_data.dims}")