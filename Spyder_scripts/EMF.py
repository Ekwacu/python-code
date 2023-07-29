# -*- coding: utf-8 -*-
"""
Created on Mon Jan 23 17:13:50 2023
Joseph Worwor
Bachelor of Atmospheric Science
Nanjing University of Information Science and Technology
@author: WORWOR
"""

import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cmaps
import geocat.viz as gv
from geocat.comp import eofunc_eofs, eofunc_pcs, month_to_season

###############################################################################
# User defined parameters and a convenience function:

# In order to specify region of the globe, time span, etc.
latS = 20.
latN = 40.
lonL = 100.
lonR = 125.

yearStart = 1979
yearEnd = 2021

neof = 3  # number of EOFs

###############################################################################
# Read in data:

# Open a netCDF data file using xarray default engine and load the data into xarrays
ds = xr.open_dataset('C:/ERSSTv5/SST_test/precip.mon.mean.nc')
ds["lon"]

###############################################################################
# Place latitudes in increasing order:

# To facilitate data subsetting

ds = ds.sortby("lat", ascending=True)

###############################################################################
# Limit data to the specified years:

startDate = f'{yearStart}-01-01'
endDate = f'{yearEnd}-12-31'

ds = ds.sel(time=slice(startDate, endDate))

###############################################################################
# Compute desired global seasonal mean using month_to_season()

# Choose the winter season (JJA)
season = "JJA"
precip = month_to_season(ds, season)

###############################################################################
# Create weights: sqrt(cos(lat))   [or sqrt(gw) ]
clat = precip['lat'].astype(np.float64)
clat = np.sqrt(np.cos(np.deg2rad(clat)))

###############################################################################
# Multiply precip by weights:

# Xarray will apply latitude-based weights to all longitudes and timesteps automatically.
# This is called "broadcasting".

wprecip = precip
wprecip['precip'] = precip['precip'] * clat

# For now, metadata for precip must be copied over explicitly; it is not preserved by binary operators like multiplication.
wprecip['precip'].attrs = ds['precip'].attrs
wprecip['precip'].attrs['long_name'] = 'Wgt: ' + wprecip['precip'].attrs['long_name']

###############################################################################
# Subset data to region:

xw = wprecip.sel(lat=slice(latS, latN), lon=slice(lonL, lonR))

###############################################################################
# Compute the EOFs:

# Transpose data to have 'time' in the first dimension
# as `eofunc` functions expects so for xarray inputs for now
xw_precip = xw["precip"].transpose('time', 'lat', 'lon')

eofs = eofunc_eofs(xw_precip, neofs=neof, meta=True)

pcs = eofunc_pcs(xw_precip, npcs=neof, meta=True)

# Change the sign of the second EOF and its time-series for
# consistent visualization purposes. See this explanation:
# https://www.ncl.ucar.edu/Support/talk_archives/2009/2015.html
# about that EOF signs are arbitrary and do not change the physical
# interpretation.
eofs[1, :, :] = eofs[1, :, :] * (-1)
pcs[1, :] = pcs[1, :] * (-1)

###############################################################################
# Normalize time series:

# Sum spatial weights over the area used.
nLon = xw.sizes["lon"]

# Bump the upper value of the slice, so that latitude values equal to latN are included.
clat_subset = clat.sel(lat=slice(latS, latN + 0.01))
weightTotal = clat_subset.sum() * nLon
pcs = pcs / weightTotal

###############################################################################
# Utility function:


# Define a utility function for creating a contour plot.
def make_contour_plot(ax, dataset):
    lat = dataset['lat']
    lon = dataset['lon']
    values = dataset.data

    # Import an NCL colormap
    cmap = cmaps.precip_diff_12lev

    cplot = ax.pcolormesh(lon,
                        lat,
                        values,
                        cmap = cmap,
                        vmin=-0.2, vmax=0.2,
                        transform=ccrs.PlateCarree())

    # Add coastlines
    ax.coastlines(linewidth=1.5)

    # Use geocat.viz.util convenience function to add minor and major tick lines
    gv.add_major_minor_ticks(ax,
                             x_minor_per_major=3,
                             y_minor_per_major=4,
                             labelsize=10)

    # Use geocat.viz.util convenience function to set axes tick values
    gv.set_axes_limits_and_ticks(ax,
                                 xticks=[100, 110, 120],
                                 yticks=[20, 30, 40])

    # Use geocat.viz.util convenience function to make plots look like NCL plots, using latitude & longitude tick labels
    gv.add_lat_lon_ticklabels(ax)

    return cplot, ax


###############################################################################
# Plot (1): Draw a contour plot for each EOF

# Generate figure and axes using Cartopy projection  and set figure size (width, height) in inches
fig, axs = plt.subplots(neof,
                        1,
                        subplot_kw={"projection": ccrs.PlateCarree()},
                        figsize=(6, 10.6))

# Add multiple axes to the figure as contour and contourf plots
for i in range(neof):
    eof_single = eofs.sel(eof=i)

    # Create contour plot for the current axes
    cplot, axs[i] = make_contour_plot(axs[i], eof_single)
    axs[i].set_extent([lonL, lonR, latS, latN], crs=ccrs.PlateCarree())
    # Use geocat.viz.util convenience function to add titles to left and right of the plot axis.
    pct = eofs.attrs['varianceFraction'].values[i] * 100
    gv.set_titles_and_labels(axs[i],
                             lefttitle=f'EOF {i + 1}',
                             lefttitlefontsize=10,
                             righttitle=f'{pct:.1f}%',
                             righttitlefontsize=10)

# Adjust subplot spacings and locations
plt.subplots_adjust(bottom=0.07, top=0.8, hspace=0.35)

# Add horizontal colorbar
cbar = plt.colorbar(cplot,
                    ax=axs,
                    orientation='horizontal',
                    shrink=0.9,
                    pad=0.05,
                    fraction=.02)
cbar.ax.tick_params(labelsize=10)

# Set a common title
axs[0].set_title(f'Precipitation: JJA: {yearStart}-{yearEnd}', fontsize=12, y=1.2)

# Show the plot
plt.savefig('EOF_pattern.png',bbox_inches='tight',dpi=300)
plt.show()

###############################################################################
# Utility function:

# Define a utility function for creating a bar plot.


def make_bar_plot(ax, dataset):
    years = list(dataset.time.dt.year)
    values = list(dataset.values)
    colors = ['lightblue' if val < 0 else 'coral' for val in values]

    ax.bar(years,
           values,
           color=colors,
           width=1.0,
           edgecolor='black',
           linewidth=0.5)
    ax.set_ylabel('mm/day')

    # Use geocat.viz.util convenience function to add minor and major tick lines
    gv.add_major_minor_ticks(ax,
                             x_minor_per_major=4,
                             y_minor_per_major=5,
                             labelsize=8)

    # Use geocat.viz.util convenience function to set axes tick values
    gv.set_axes_limits_and_ticks(ax,
                                 xticks=np.linspace(1980, 2020, 9),
                                 xlim=[1978.5, 2021.5])

    return ax


###############################################################################
# Plot (2): Produce a bar plot for each EOF.

# Generate figure and axes using Cartopy projection and set figure size (width, height) in inches
fig, axs = plt.subplots(neof, 1, constrained_layout=True, figsize=(6, 7.5))

# Add multiple axes to the figure as bar-plots
for i in range(neof):
    eof_single = pcs.sel(pc=i)

    axs[i] = make_bar_plot(axs[i], eof_single)
    pct = eofs.attrs['varianceFraction'].values[i] * 100
    gv.set_titles_and_labels(axs[i],
                             lefttitle=f'EOF {i + 1}',
                             lefttitlefontsize=10,
                             righttitle=f'{pct:.1f}%',
                             righttitlefontsize=10)

# Set a common title
axs[0].set_title(f'precip: JJA: {yearStart}-{yearEnd}', fontsize=12, y=1.12)

# Show the plot
plt.savefig('EOF_ts.png',bbox_inches='tight',dpi=300)
plt.show()

