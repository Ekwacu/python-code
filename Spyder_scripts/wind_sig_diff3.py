# -*- coding: utf-8 -*-
"""
Created on Thu Mar 16 12:38:52 2023

@author: Jesse
"""


import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter
#from matplotlib.colors import ListedColormap
#import cmaps
import geocat.viz as gv



def open_sst(file_path):
    
    ds = xr.open_dataset(file_path)
    ds.coords['lon']=(ds.coords['lon']+180)%360-180
    ds=ds.sortby(ds.lon) 
    
    ##-- group data into yearly means
    #ds1=ds.sel(lev=95000)
    SST=ds.adia_expan.mean('time') 
        
    return SST
        
#open the files   
AA = open_sst("C:/Diabatic_data/extra3/advection/adiabtic_expansion2_1961-1981_mon_950hpa_Af.nc")
BB = open_sst("C:/Diabatic_data/extra3/advection/adiabtic_expansion2_2001-2021_mon_950hpa_Af.nc")




#compute for the difference in the variables
diff = AA - BB



def open_diff_sig(file_path):

    ds = xr.open_dataset(file_path)
    ds.coords['lon']=(ds.coords['lon']+180)%360-180
    ds=ds.sortby(ds.lon)

    diff_sig=ds.sig
        
    return diff_sig
        
#open the siginificance files    
diff_sig = open_diff_sig('C:/Diabatic_data/extra3/adiabtic_expansion_sig_pvalue_PHE_P1-P3_diff.nc')



## plotting
# Create a wider than normal figure to support our many plots
fig = plt.figure(figsize=(18,9))
plt.gcf().subplots_adjust(hspace=0.1, wspace=1)
plt.rcParams["font.family"] = "Arial"


def Plot(row, col, pos, diff, diff_sig, clevs, cbar_label, title):
    # Generate axes, using cartopy, drawing coastlines, and adding features
    projection = ccrs.PlateCarree()#(central_longitude=0)
    ax = fig.add_subplot(row, col, pos, projection=projection)
    plt.xlim([-20, 55])
    plt.ylim([-35, 40])
    plt.gca().set_yticks(np.arange(-30,50,10),crs=ccrs.PlateCarree())
    plt.gca().set_xticks(np.arange(-20,55,10),crs=ccrs.PlateCarree())
   # plt.gca().set_yticks(np.arange(-60,90,30),crs=ccrs.PlateCarree())
   # plt.gca().set_xticks(np.arange(-150,300,50),crs=ccrs.PlateCarree())
    lon_formatter=LongitudeFormatter(degree_symbol=''); lat_formatter=LatitudeFormatter(degree_symbol='')
    ax.xaxis.set_major_formatter(lon_formatter); ax.yaxis.set_major_formatter(lat_formatter);ax.tick_params(labelsize=20)
    xticks = ax.xaxis.get_major_ticks(); xticks[2].set_visible(True)
    
   # ax.coastlines(resolution='10m', color='black',linewidth=10)
    #ax.coastlines(resolution='10m', color='black', linewidth=0.7)
    #ax.add_feature(cfeature.BORDERS, linewidth=0.1)
    ax.set_xlabel("")
    ax.set_ylabel("")

    # Contourf-plot data
    temp = diff.plot.contourf(ax=ax,
                           transform=ccrs.PlateCarree(),
                           levels=clevs,
                           cmap=plt.cm.RdBu_r,
                           add_colorbar=False,
                           extend='both')
    
    #ax.coastlines(resolution='10m', color='black',linewidth=1)
    ax.coastlines(linewidth=1.8)
    ax.add_feature(cfeature.BORDERS, linewidth=1)
    
    # Plot Hatch
    pval = diff_sig
    cond = (pval <= 0.05)
    ## Mask out the areas that do not satisfy the conditions
    sig = pval.where(cond)
    
    ## make a hatch of significance
    plt.contourf(sig.lon,sig.lat,sig,hatches=['...'],alpha=0,
                 transform=ccrs.PlateCarree()) 
    
        # Add color bar
    cbar = plt.colorbar(temp,
                        orientation='vertical',
                        shrink=0.8,
                        extendfrac='auto', 
                        extendrect=False, 
                        drawedges=True)

    cbar.ax.tick_params(labelsize=18)
    cbar.set_ticks(clevs)
    cbar.ax.set_title(cbar_label, size=18)

    # Use geocat.viz.util convenience function to set titles and labels without calling several matplotlib functions
    gv.set_titles_and_labels(ax,
                            maintitle="",
                            lefttitle=title,
                            lefttitlefontsize=22,
                            righttitle="",
                            righttitlefontsize=22,
                            xlabel="",
                            ylabel="")


# define the levels for each variabledddddd
#clevs  = [-0.6 ,-0.5, -0.4, -0.3, -0.2, -0.1, 0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6]
clevs  = [-0.4 ,-0.3, -0.2, -0.1, 0, 0.1, 0.2, 0.3, 0.4, 0.5]
#clevs  = [-5, -3, -1, 0, 1, 3, 5]
#clevs  = [ -3, -2, -1, 0, 1, 2, 3,]
Plot(1, 1, 1, diff, diff_sig, clevs, '(degC day$^{-1}$)', "Adiabatic expansion/compression mean difference (P1-P3)")

fig.savefig('C:/Diabatic_data/extra3/adiabtic_expansion_sig-PHE_P1-P3_diff.png', bbox_inches='tight', pad_inches = 0.1, dpi=300)