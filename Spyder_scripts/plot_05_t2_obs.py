# -*- coding: utf-8 -*-
"""
Created on Sun May 21 18:44:30 2023

@author: Joan
"""
import cdsapi
import xarray as xr
import numpy as np
import pandas as pd
from pathlib import Path

import matplotlib.pyplot as plt
import matplotlib.colors as colors
import matplotlib.cm as cm
import cartopy
import cartopy.crs as ccrs
import shapefile as shp
from descartes import PolygonPatch
#uploaded.keys()
home_dir = './' ## read in the data
data_dir = './'
data = xr.open_dataset(home_dir + data_dir + 't2_obs_days.nc')  ##input file                          
print(data.data_vars)

sf = shp.Reader(home_dir + data_dir + 'lv_basin.shp')  ## Used in ploting



data1= data.sel(time='2022-04-24')
t2m1  = data1['t2m']

data2= data.sel(time='2022-04-25')
t2m2  = data2['t2m']

data3= data.sel(time='2022-04-26')
t2m3  = data3['t2m']

##################### Set cartopy projection and plot map ######################
color ='jet'
cbrange = np.arange(10, 30, 1)


fig, ax = plt.subplots(1,3, subplot_kw=dict(projection=ccrs.PlateCarree()), figsize=(10, 5))



########################## 2022-04-24 ###############################
cb=ax[0].contourf(t2m1.longitude, t2m1.latitude, t2m1[0,0,:,:], levels=cbrange,extend='both',
                                cmap=color, transform=ccrs.PlateCarree())
## Plotting the figure
for poly in sf.shapes():
    poly_geo=poly.__geo_interface__
    ax[0].add_patch(PolygonPatch(poly_geo, fc='None', ec='black', linewidth= 2, alpha=1, fill='False', zorder=2))
ax[0].coastlines(resolution='50m', color='black')
ax[0].add_feature(cartopy.feature.BORDERS, linestyle='-',alpha=0.1)
ax[0].set_title ('(a)',loc='left')
###################### Adding the labels of grids ###############################
gl = ax[0].gridlines(draw_labels=True, linewidth=0.1, color='gray')
gl.top_labels = False  ###Changing these to true or false will add or remove the lon/lat labelling 
gl.right_labels = False
gl.bottom_labels = True
gl.left_labels = True


########################## 2022-04-25 ###############################
ax[1].contourf(t2m2.longitude, t2m2.latitude, t2m2[0,0,:,:], levels=cbrange,extend='both',
                                cmap=color, transform=ccrs.PlateCarree())
for poly in sf.shapes():
    poly_geo=poly.__geo_interface__
    ax[1].add_patch(PolygonPatch(poly_geo, fc='None', ec='black', linewidth= 2, alpha=1, fill='False', zorder=2))
ax[1].coastlines(resolution='50m', color='black')
ax[1].add_feature(cartopy.feature.BORDERS, linestyle='-', alpha=0.1)
ax[1].set_title ('(b)',loc='left')
###################### Adding the labels of grids ###############################
gl = ax[1].gridlines(draw_labels=True, linewidth=0.1, color='gray')
gl.top_labels = False  ###Changing these to true or false will add or remove the lon/lat labelling 
gl.right_labels = False
gl.bottom_labels = True
gl.left_labels = True


########################## 2022-04-26 ###############################
ax[2].contourf(t2m3.longitude, t2m3.latitude, t2m3[0,0,:,:], levels=cbrange,extend='both',
                                cmap=color, transform=ccrs.PlateCarree())
for poly in sf.shapes():
    poly_geo=poly.__geo_interface__
    ax[2].add_patch(PolygonPatch(poly_geo, fc='None', ec='black', linewidth= 2, alpha=1, fill='False', zorder=2))
ax[2].coastlines(resolution='50m', color='black')
ax[2].add_feature(cartopy.feature.BORDERS, linestyle='-',alpha=0.1)
ax[2].set_title ('(c)',loc='left')
###################### Adding the labels of grids ###############################
gl = ax[2].gridlines(draw_labels=True, linewidth=0.1, color='gray')
gl.top_labels = False  ###Changing these to true or false will add or remove the lon/lat labelling 
gl.right_labels = False
gl.bottom_labels = True
gl.left_labels = True


##########################  Setting titles and colorbar ##############################
#cax = plt.axes([0.85,0.1,0.075,0.8])
#fig.colorbar(cb,ax=ax,shrink = 0.5,location='bottom', orientation='horizontal')
plt.colorbar(cb,ax=ax,shrink = 0.5, label='W $m^{-2}$',orientation='horizontal')