"""
Filename:    plot_IVT.py
Author:      Deanna Nash, dnash@ucsd.edu
Description: Output figures for GEFS IVT for selected dates.
"""

import os, sys
import xarray as xr
import pandas as pd
import itertools ## need this for the cbarticks
import datetime

## import plotting modules
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
import cartopy.crs as ccrs
from matplotlib.colorbar import Colorbar # different way to handle colorbar

# import personal modules
sys.path.append('modules')
from read_preprocessed_GEFS import load_archive_GEFS_forecast
import domains
from plotter import draw_basemap
import cw3ecmaps as ccmaps

## Set domain name
domain_name = 'intwest'
varname = 'ivt'

## Set up data based on domain name
ext = domains.extent[domain_name]['ext']
dx = domains.extent[domain_name]['xticks']
dy = domains.extent[domain_name]['yticks']
figsize = domains.extent[domain_name]['figsize']

## load GEFS data
init_date_lst = ['20250213', '20250214', '20241218', '20241219']

for i, init_date in enumerate(init_date_lst):
    ds = load_archive_GEFS_forecast(date=init_date, varname='IVT', F=0)
    ## compute ensemble mean
    ds = ds.mean('number')

    ## Set domain name
    domain_name = 'intwest'
    varname = 'ivt'
    
    data = ds.ivt.values
    lats = ds.lat.values
    lons = ds.lon.values
    
    ## now let's plot and save figure with our domain
    
    fmt = 'png'
    nrows = 1
    ncols = 2
    
    # Set up projection
    datacrs = ccrs.PlateCarree()  ## the projection the data is in
    mapcrs = domains.extent[domain_name]['ccrs']
    
    fig = plt.figure(figsize=figsize)
    fig.dpi = 300
    fname = 'figs/IVT_{0}_{1}'.format(domain_name, init_date)
    
    ## Use gridspec to set up a plot with a series of subplots that is
    ## n-rows by n-columns
    gs = GridSpec(nrows, ncols, height_ratios=[1], width_ratios = [1, 0.05], wspace=0.025, hspace=0.05)
    ## use gs[rows index, columns index] to access grids
    
    ## add basemap
    ax = fig.add_subplot(gs[0, 0], projection=mapcrs)
    ax = draw_basemap(ax, extent=ext, xticks=dx, yticks=dy, grid= True, left_lats=True, right_lats=False)
    ax.set_extent(ext, datacrs)
    
    ## grab custom cmap based on domain and varname
    if (varname == 'ivt') & (domain_name == 'intwest'):
        cmap, norm, bnds, cbarticks, cbarlbl = ccmaps.cmap('ivt_intwest') # get cmap from our custom function
    else: 
        cmap, norm, bnds, cbarticks, cbarlbl = ccmaps.cmap(varname) # get cmap from our custom function
    ## plot random data
    cf = ax.contourf(lons, lats, data, transform=datacrs,
                     levels=bnds, cmap=cmap, norm=norm, alpha=0.9)
    
    # Add color bar
    cbax = plt.subplot(gs[0, -1]) # colorbar axis (first row, last column)
    cbarticks = list(itertools.compress(bnds, cbarticks)) ## this labels the cbarticks based on the cmap dictionary
    cb = Colorbar(ax = cbax, mappable = cf, orientation = 'vertical', ticklocation = 'right', ticks=cbarticks)
    cb.set_label(cbarlbl, fontsize=11)
    cb.ax.tick_params(labelsize=12)
    
    ## add fig label
    t = pd.to_datetime(ds.init_date.values).strftime("%Y-%m-%d %H:00")
    ax.set_title(t, loc='left')
    
    ## save figure
    fig.savefig('%s.%s' %(fname, fmt), bbox_inches='tight', dpi=fig.dpi)