# -*- coding: utf-8 -*-
"""
Created on Thu May 11 10:09:02 2023

@author: Samuel-NUIST
"""

import xarray as xr
import numpy as np
from scipy import stats
#from statsmodels.stats import weightstats as stests 

def open_sst(file_path):
    
    ds = xr.open_dataset(file_path)
    ds.coords['longitude']=(ds.coords['longitude']+180)%360-180
    ds=ds.sortby(ds.longitude) 
    
    ##-- group data into yearly means

    SST=ds.vws.groupby('time.year').mean('time')
    return SST
   
        
#open the files   
AA = open_sst("C:/Dinis_MOZ/Wind/vwsf_2000_2020.nc")
BB = open_sst("C:/Dinis_MOZ/Wind/VWSf_1980_2020.nc")


## perform a student's t-test using the 
diff_val = np.zeros((141, 321))
p_value = np.zeros((141, 321))

for i in np.arange(141):
    for j in np.arange(321):
        
        try:
            diff_val[i,j] = stats.ttest_ind(AA[:,i,j], BB[:,i,j], equal_var=True, alternative='two-sided').statistic  
            p_value[i,j]  = stats.ttest_ind(AA[:,i,j], BB[:,i,j], equal_var=True, alternative='two-sided').pvalue
        except:
            diff_val[i,j] = np.nan
            p_value[i,j]  = np.nan


#z_test ,p_val = stests.ztest(data1, x2 = data2, value = 0, alternative = 'two-sided')
## define data as xarray dataset and save as netcdf
lon=AA['longitude'].values
lat=AA['latitude'].values


var2=xr.DataArray(data=p_value, dims=('latitude', 'longitude'), coords={'latitude':lat, 'longitude':lon}, 
                     attrs=dict(description="significance of the difference b/n active and inactive based on t-test",),).rename('sig')

## save data as netcdf
var2.to_netcdf('C:/Dinis_MOZ/Wind/ttest_VWS_2000_2020.nc', mode='w')
