"""
Created on Tue Apr 26 17:48:40 2022
Author: Jesse Kisembe


*Description*
This is a simple code that plots the spatial correlations between area averaged hotspell days and annual SSTs.

"""

import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER


#function to add plot_extras and reduce code lines
def plot_extras(ax):
    ax.coastlines()
    #ax.add_feature(cf.BORDERS)
    ax.set_extent([-90,180, -50, 50])
    gl=ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True, linewidth=1.0, color='gray', alpha=0.8, linestyle='--')
    gl.xlabels_top = False; gl.ylabels_left = True; gl.ylabels_right=False; gl.xlines = True
    gl.xformatter = LONGITUDE_FORMATTER; gl.yformatter = LATITUDE_FORMATTER



fig=plt.figure(figsize=(18,7))
fig.subplots_adjust(hspace=0.5, wspace=0.05,top=0.95, bottom=0.2, left=0.40, right=0.975)


nrows=3 #number of plot rows
ncols=2 #number of plot columns

#DATA_PATH = "E:/DGGCS/students/Samuel/data/cor_ssts_hotspells/"
DATA_PATH= "F:/Folder_2022/Research/Hotspells_spatial/cor_ssts_hotspells/Hotspell_tmin-12-cor/"
regions = ['EastAfrica','Burundi','Kenya','Rwanda','Uganda','Tanzania']


for index,region in enumerate(regions):
    
    '''
    open the data files
    '''
    #print(DATA_PATH +region +'_cor_sst_hotspell_tmax_15.nc')
    dset = xr.open_dataset(DATA_PATH +region +'_cor_sst_hotspell_tmin_12.nc')
    dset.cor.attrs = {'units':'','long_name': ''}
    dset.cor_sig.attrs = {'units':'','long_name': ''}
    
    '''
    setting the layout on which the stippling will be used.
    '''
    X = dset.cor.lon.values
    Y = dset.cor.lat.values
    lon,lat = np.meshgrid(X,Y)
    
    '''
    plotting out the correlation maps
    '''
    
    ax=fig.add_subplot(nrows,ncols,index+1, projection=ccrs.PlateCarree())
    dset.cor.plot(cmap='RdBu',vmin=-1.0, vmax=1.0,cbar_kwargs={'orientation':'vertical','label':''})
    plt.contourf(lon,lat,dset.cor_sig,hatches=['xxxx'],alpha=0)
    plot_extras(ax)
    ax.set_title(region,fontsize=12,fontweight='bold')
plt.savefig('F:/Folder_2022/Research/cor_plots/Tmin_12_cor_sst_hotspell.png',dpi=300, transparent=False, bbox_inches='tight', pad_inches=0.1)
#plt.savefig('E:/DGGCS/students/Samuel/plots/Tmax_15_cor_sst_hotspell.png',dpi=300,transparent=False,bbox_inches='tight',pad_inches=0.1)