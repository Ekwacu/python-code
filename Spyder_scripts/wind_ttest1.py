# -*- coding: utf-8 -*-
"""
Created on Thu Mar 16 12:16:11 2023

@author: Jesse
"""

import xarray as xr
import numpy as np
from scipy import stats


def open_sst(file_path):
    
    ds = xr.open_dataset(file_path)
    ds.coords['lon']=(ds.coords['lon']+180)%360-180
    ds=ds.sortby(ds.lon) 
    
    ##-- group data into yearly means
    #ds1=ds.sel(levels=1)
    SST=ds.adia_expan.groupby('time.year').mean('time') 
        
    return SST
        
#open the files   
AA = open_sst("C:/Diabatic_data/extra3/advection/adiabtic_expansion2_1982-2000_mon_950hpa_Af.nc")
BB = open_sst("C:/Diabatic_data/extra3/advection/adiabtic_expansion2_2001-2021_mon_950hpa_Af.nc")


## perform a student's t-test using the (lat,lon)
diff_val = np.zeros((73, 121))
p_value = np.zeros((73, 121))

for i in np.arange(73):
    for j in np.arange(121):
        
        try:
            diff_val[i,j] = stats.ttest_ind(AA[:,i,j], BB[:,i,j], equal_var=True, alternative='two-sided').statistic  
            p_value[i,j]  = stats.ttest_ind(AA[:,i,j], BB[:,i,j], equal_var=True, alternative='two-sided').pvalue
        except:
            diff_val[i,j] = np.nan
            p_value[i,j]  = np.nan



## define data as xarray dataset and save as netcdf
lon=AA['lon'].values
lat=AA['lat'].values


var2=xr.DataArray(data=p_value, dims=('lat', 'lon'), coords={'lat':lat, 'lon':lon}, 
                     attrs=dict(description="significance of the difference b/n active and inactive based on t-test",),).rename('sig')

## save data as netcdf
var2.to_netcdf('C:/Diabatic_data/extra3/adiabtic_expansion_sig_pvalue_PHE_P2-P3_diff.nc', mode='w')
