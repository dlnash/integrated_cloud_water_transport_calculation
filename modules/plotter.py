"""
Filename:    plotter.py
Author:      Deanna Nash, dnash@ucsd.edu
Description: Functions for plotting maps
"""

# Import Python modules

import os, sys
import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
import matplotlib.ticker as mticker
import colorsys
import cmocean.cm as cmo
from matplotlib.colors import LinearSegmentedColormap # Linear interpolation for color maps
import matplotlib.patches as mpatches
from matplotlib import cm, colors as clr
from matplotlib.colorbar import Colorbar # different way to handle colorbar
import pandas as pd
import matplotlib.gridspec as gridspec
import cw3ecmaps as ccmaps
from PIL import Image
from matplotlib import font_manager as fm

def set_cw3e_font(current_dpi, scaling_factor):
    fm.fontManager.addfont('utils/fonts/helvetica.ttc')

    plt.rcParams.update({
                    'font.family' : 'Helvetica',
                    'figure.dpi': current_dpi,
                    'font.size': 8 * scaling_factor, #changes axes tick label
                    'axes.labelsize': 8 * scaling_factor,
                    'axes.titlesize': 8 * scaling_factor,
                    'xtick.labelsize': 8 * scaling_factor,#do nothing
                    'ytick.labelsize': 8 * scaling_factor, #do nothing
                    'legend.fontsize': 5 * scaling_factor,
                    'lines.linewidth': 0.7 * scaling_factor,
                    'axes.linewidth': 0.2 * scaling_factor,
                    'legend.fontsize': 12 * scaling_factor,
                    'xtick.major.width': 0.8 * scaling_factor,
                    'ytick.major.width': 0.8 * scaling_factor,
                    'xtick.minor.width': 0.6 * scaling_factor,
                    'ytick.minor.width': 0.6 * scaling_factor,
                    'lines.markersize': 6 * scaling_factor
                })
    
def plot_cw3e_logo(ax, orientation):
    ## location of CW3E logo
    if orientation == 'horizontal':
        im = '/common/CW3E_Logo_Suite/1-Horzontal-PRIMARY_LOGO/Digital/JPG-RGB/CW3E-Logo-Horizontal-FullColor-RGB.jpg'
    else:
        im = '/common/CW3E_Logo_Suite/5-Vertical-Acronym_Only/Digital/PNG/CW3E-Logo-Vertical-Acronym-FullColor.png'
    img = np.asarray(Image.open(im))
    ax.imshow(img)
    ax.axis('off')
    return ax

def plot_terrain(ax, ext):
    fname = '/work/bkawzenuk_work/Maps/data/ETOPO1_Bed_c_gmt4.grd'
    datacrs = ccrs.PlateCarree()
    grid = xr.open_dataset(fname)
    grid = grid.where(grid.z > 0) # mask below sea level
    grid = grid.sel(x=slice(ext[0], ext[1]), y=slice(ext[2], ext[3]))
    cs = ax.pcolormesh(grid.x, grid.y, grid.z,
                        cmap=cmo.gray_r, transform=datacrs, alpha=0.7)
    
    return ax

def draw_basemap(ax, datacrs=ccrs.PlateCarree(), extent=None, xticks=None, yticks=None, grid=False, left_lats=True, right_lats=False, bottom_lons=True, mask_ocean=False, coastline=True):
    """
    Creates and returns a background map on which to plot data. 
    
    Map features include continents and country borders.
    Option to set lat/lon tickmarks and draw gridlines.
    
    Parameters
    ----------
    ax : 
        plot Axes on which to draw the basemap
    
    datacrs : 
        crs that the data comes in (usually ccrs.PlateCarree())
        
    extent : float
        Set map extent to [lonmin, lonmax, latmin, latmax] 
        Default: None (uses global extent)
        
    grid : bool
        Whether to draw grid lines. Default: False
        
    xticks : float
        array of xtick locations (longitude tick marks)
    
    yticks : float
        array of ytick locations (latitude tick marks)
        
    left_lats : bool
        Whether to add latitude labels on the left side. Default: True
        
    right_lats : bool
        Whether to add latitude labels on the right side. Default: False
        
    Returns
    -------
    ax :
        plot Axes with Basemap
    
    Notes
    -----
    - Grayscale colors can be set using 0 (black) to 1 (white)
    - Alpha sets transparency (0 is transparent, 1 is solid)
    
    """
    ## some style dictionaries
    kw_ticklabels = {'color': 'black', 'weight': 'light'}
    kw_grid = {'linewidth': 0.6, 'color': 'k', 'linestyle': '--', 'alpha': 0.4}
    kw_ticks = {'length': 4, 'width': 0.5, 'pad': 2, 'color': 'black',
                             'labelcolor': 'dimgray'}

    # Use map projection (CRS) of the given Axes
    mapcrs = ax.projection    
    
    # Add map features (continents and country borders)
    ax.add_feature(cfeature.LAND, facecolor='0.9')      
    ax.add_feature(cfeature.BORDERS, edgecolor='0.4', linewidth=0.8)
    ax.add_feature(cfeature.STATES, edgecolor='0.2', linewidth=0.2)
    if coastline == True:
        ax.add_feature(cfeature.COASTLINE, edgecolor='0.4', linewidth=0.8)
    if mask_ocean == True:
        ax.add_feature(cfeature.OCEAN, edgecolor='0.4', zorder=12, facecolor='white') # mask ocean
        
    ## Tickmarks/Labels
    ## Add in meridian and parallels
    if mapcrs == ccrs.NorthPolarStereo(central_longitude=0.):
        gl = ax.gridlines(draw_labels=False,
                      linewidth=.5, color='black', alpha=0.5, linestyle='--')
    elif mapcrs == ccrs.SouthPolarStereo(central_longitude=0.):
        gl = ax.gridlines(draw_labels=True,
                      linewidth=.5, color='black', alpha=0.5, linestyle='--')
        
    else:
        gl = ax.gridlines(crs=datacrs, draw_labels=True, **kw_grid)
        gl.top_labels = False
        gl.left_labels = left_lats
        gl.right_labels = right_lats
        gl.bottom_labels = bottom_lons
        gl.xlocator = mticker.FixedLocator(xticks)
        gl.ylocator = mticker.FixedLocator(yticks)
        gl.xformatter = LONGITUDE_FORMATTER
        gl.yformatter = LATITUDE_FORMATTER
        gl.xlabel_style = kw_ticklabels
        gl.ylabel_style = kw_ticklabels
        
        # apply tick parameters
        ax.set_xticks(xticks, crs=datacrs)
        ax.set_yticks(yticks, crs=datacrs)
        plt.yticks(color='w', size=1) # hack: make the ytick labels white so the ticks show up but not the labels
        plt.xticks(color='w', size=1) # hack: make the ytick labels white so the ticks show up but not the labels
        ax.ticklabel_format(axis='both', style='plain')
    
    ## Gridlines
    # Draw gridlines if requested
    if (grid == True):
        gl.xlines = True
        gl.ylines = True
    if (grid == False):
        gl.xlines = False
        gl.ylines = False

    ## Map Extent
    # If no extent is given, use global extent
    if extent is None:        
        ax.set_global()
        extent = [-180., 180., -90., 90.]
    # If extent is given, set map extent to lat/lon bounding box
    else:
        ax.set_extent(extent, crs=datacrs)
    
    return ax

def add_subregion_boxes(ax, subregion_xy, width, height, ecolor, datacrs):
    '''This function will add subregion boxes to the given axes.
    subregion_xy 
    [[ymin, xmin], [ymin, xmin]]
    '''
    for i in range(len(subregion_xy)):
        ax.add_patch(mpatches.Rectangle(xy=subregion_xy[i], width=width[i], height=height[i],
                                        fill=False,
                                        edgecolor=ecolor,
                                        linewidth=1.0,
                                        transform=datacrs,
                                        zorder=100))
        
    return ax