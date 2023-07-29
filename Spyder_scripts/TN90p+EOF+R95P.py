#!/usr/bin/env python
# coding: utf-8

# In[12]:


#Import Library
import icclim
import datetime

#read infile and outfile 
files = "tanzania_tmin_data.nc"
out_f = "D:\mbawala\TN90P.nc"


# In[13]:


#Time range of the data
dt1 = datetime.datetime(1981, 1, 1)
dt2 = datetime.datetime(2022, 12, 31)

#Index calculation

icclim.index(
        index_name="TN90P",
        in_files= files,
        var_name=None,
        slice_mode=("year"),
        time_range=[dt1, dt2],
        out_file=out_f,
        base_period_time_range=[dt1, dt2],
        only_leap_years=False,
        ignore_Feb29th=False,
        interpolation= "median_unbiased",
        netcdf_version="NETCDF4",
        save_thresholds=False,
        logs_verbosity= "LOW",
        date_event=False
        
    )


# In[22]:


import cartopy.crs as ccrs
import cartopy.feature as cfeature
from netCDF4 import Dataset
import matplotlib.pyplot as plt
import numpy as np
import xarray as xr
import numpy as np
from eofs.standard import Eof


# In[23]:


data = xr.open_dataset('TNX.nc')


# In[24]:


#print(data)


# In[25]:


st = data['TNx'].values
lons = data['lon'].values
lats = data['lat'].values
#print(st)


# In[26]:


# Create an EOF solver to do the EOF analysis. Square-root of cosine of
# latitude weights are applied before the computation of EOFs.
coslat = np.cos(np.deg2rad(lats))
wgts = np.sqrt(coslat)[..., np.newaxis]
solver = Eof(st, weights=wgts)


# In[27]:


# Retrieve the leading EOF, expressed as the correlation between the leading
# PC time series and the input SST anomalies at each grid point, and the
# leading PC time series itself.
eof1 = solver.eofsAsCorrelation(neofs=1)
pc1 = solver.pcs(npcs=1, pcscaling=1)


# In[28]:



# Plot the leading EOF expressed as correlation in the Pacific domain.
clevs = np.linspace(-1, 1, 11)
ax = plt.axes(projection=ccrs.PlateCarree(central_longitude=34))
fill = ax.contourf(lons, lats, eof1.squeeze(), clevs,
                   transform=ccrs.PlateCarree(), cmap=plt.cm.RdBu_r)
ax.add_feature(cfeature.LAND, facecolor='w', edgecolor='k')
cb = plt.colorbar(fill, orientation='horizontal')
cb.set_label('TEMPERATURE', fontsize=12)
plt.title('EOF1: TN90p FOR MINIMUM TEMPERATURE', fontsize=12)
plt.savefig('D:\mbawala\TNx1.png', dpi = 300)

##### Plot the leading PC time series.

plt.figure()
years = range(1981, 2023)
var = solver.varianceFraction()
color=[]

plt.bar(years, pc1[:,0])
plt.axhline(0, color='k')
plt.title('PC1 TNX (1981-2022)')
plt.xlabel('Year')
plt.ylabel('TNx')
plt.xlim(1981, 2022)
plt.ylim(-3, 3)
print(years)
plt.show()
plt.savefig('D:\mbawala\TNN.png', dpi = 300)
var = solver.varianceFraction()
color=[]
#,color='b', linewidth=2


# In[ ]:


1. R95p INDEX Calculation
*********************************************************************************
#Import Library
import icclim
import datetime

#read infile and outfile 
files = ["chirps_daily_1981-2022.nc"]
out_f = "D:\IMSA\CDO\Kokotoa\R95Ptotal_1981-2022_NDJFMA.nc"

#Time range of the data
dt1 = datetime.datetime(1981, 1, 1)
dt2 = datetime.datetime(2022, 12, 31)

#Index calculation

icclim.index(
       index_name="R95PTOT",
       in_files= files,
       var_name='precip',
       slice_mode=("month", [11,12,1,2,3,4,]),
       time_range=[dt1, dt2],
       out_file=out_f,
       base_period_time_range=[dt1, dt2],
       only_leap_years=False,
       ignore_Feb29th=False,
       interpolation= "median_unbiased",
       netcdf_version="NETCDF4",
       save_thresholds=False,
       logs_verbosity= "LOW",
       date_event=False
       
   )
***************************************************************************************

