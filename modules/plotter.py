"""
Filename:    plotter.py
Author:      Deanna Nash, dnash@ucsd.edu and Ricardo Vilela, rbatistavilela@ucsd.edu
Description: Functions for plotting maps
"""

import os, sys
import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter
import matplotlib.ticker as mticker
import colorsys
from datetime import datetime, timedelta
from matplotlib.colors import LinearSegmentedColormap # Linear interpolation for color maps
import matplotlib.patches as mpatches
from matplotlib import cm, colors as clr
from matplotlib.colorbar import Colorbar # different way to handle colorbar
import pandas as pd
from matplotlib.gridspec import GridSpec
import itertools
from PIL import Image
from matplotlib import font_manager as fm
from matplotlib.ticker import FuncFormatter
from scipy.ndimage import gaussian_filter
import copy
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import cartopy.feature as cf
import matplotlib.image as image


import matplotlib as mpl
mpl.use('agg')

import cw3ecmaps as ccmaps
import domains as domains
from plot_utils import crop_to_domain, plot_vectors, plot_contours, create_folder

current_dpi=300

def plot_cw3e_logo(ax, orientation):
    ## location of CW3E logo
    if orientation == 'horizontal':
        im = '/data/projects/operations/data/CW3E_Logo_Suite/1-Horzontal-PRIMARY_LOGO/Digital/JPG-RGB/CW3E-Logo-Horizontal-FullColor-RGB.jpg'
    else:
        im = '/data/projects/operations/data/CW3E_Logo_Suite/5-Vertical-Acronym_Only/Digital/PNG/CW3E-Logo-Vertical-Acronym-FullColor.png'
    img = np.asarray(Image.open(im))
    ax.imshow(img)
    ax.axis('off')
    return ax

def set_cw3e_font(current_dpi, scaling_factor, font_file):
    fm.fontManager.addfont(font_file)

    plt.rcParams.update({
                    'font.family' : 'Helvetica',
                    'figure.dpi': current_dpi,
                    'font.size': 6 * scaling_factor, #changes axes tick label
                    'axes.labelsize': 8 * scaling_factor,
                    'axes.titlesize': 6 * scaling_factor, #title fontsize
                    'xtick.labelsize': 8 * scaling_factor,#do nothing
                    'ytick.labelsize': 8 * scaling_factor, #do nothing
                    'legend.fontsize': 5 * scaling_factor,
                    'lines.linewidth': 0.4 * scaling_factor,
                    'axes.linewidth': 0.2 * scaling_factor,
                    'legend.fontsize': 12 * scaling_factor,
                    'xtick.major.width': 0.8 * scaling_factor,
                    'ytick.major.width': 0.8 * scaling_factor,
                    'xtick.minor.width': 0.6 * scaling_factor,
                    'ytick.minor.width': 0.6 * scaling_factor,
                    'lines.markersize': 6 * scaling_factor
                })



def draw_basemap(ax, domain_name, datacrs=ccrs.PlateCarree(), logofile=None, logopos='upper_left', extent=None, grid=True, mask_ocean=False, coastline=True):
    """
    Creates and returns a background map on which to plot data. 
    
    Map features include continents and country borders.
    Option to set lat/lon tickmarks and draw gridlines.
    
    Parameters
    ----------
    ax : 
        plot Axes on which to draw the basemap
    
    domain_name:
        Domain name from domains.py 

    datacrs: 
        Input data crs


    Returns
    -------
    ax :
        plot Axes with Basemap
    
    Notes
    -----
    - Grayscale colors can be set using 0 (black) to 1 (white)
    - Alpha sets transparency (0 is transparent, 1 is solid)
    
    """
    base_dpi=100


    scaling_factor = (current_dpi / base_dpi)**0.5

    
    ## Set up data based on domain name
    extent = domains.extent[domain_name]["ext"]

    # Set up projection
    mapcrs = domains.extent[domain_name]['ccrs']

    # Set tick/grid locations
    yticks = domains.extent[domain_name]["yticks"]
    xticks = domains.extent[domain_name]["xticks"]

 
    ## some style dictionaries
    kw_ticklabels = {'color': 'black', 'weight': 'light'}
    kw_grid = {'linewidth': 0.35*scaling_factor, 'color': 'k', 'linestyle': '--', 'alpha': 0.4}
    kw_ticks = {'length': 4, 'width': 0.5, 'pad': 2, 'color': 'black',
                             'labelcolor': 'dimgray'}

    # Use map projection (CRS) of the given Axes
    mapcrs = ax.projection    
    
    # Add map features (continents and country borders)
    ax.add_feature(cfeature.LAND.with_scale('10m'), facecolor='0.9')      
    ax.add_feature(cfeature.BORDERS.with_scale('10m'), edgecolor='0.4', linewidth=0.8)
    ax.add_feature(cfeature.STATES.with_scale('10m'), edgecolor='0.2', linewidth=0.2)
    ax.set_aspect('equal',adjustable='box')
    if coastline == True:
        ax.add_feature(cfeature.COASTLINE.with_scale('10m'), edgecolor='0.4', linewidth=0.8)
    if mask_ocean == True:
        ax.add_feature(cfeature.OCEAN, edgecolor='0.4', zorder=12, facecolor='white') # mask ocean
        
    ## Tickmarks/Labels
    ## Add in meridian and parallels

    if (mapcrs == ccrs.SouthPolarStereo(central_longitude=0.)) or (mapcrs == ccrs.NorthPolarStereo(central_longitude=0.)) :
        #gl = ax.gridlines(draw_labels=True,
        #                linewidth=.5, color='black', alpha=0.5, linestyle='--')
        gl = ax.gridlines(crs=datacrs, draw_labels=True, **kw_grid)
        gl.xlabel_style = kw_ticklabels
        gl.xlocator = mticker.FixedLocator(xticks)
        gl.ylocator = mticker.FixedLocator(yticks)
        gl.top_labels = True
        gl.left_labels = True #not working
        gl.right_labels = True #not working
        gl.bottom_labels = True
  
    else:
        gl = ax.gridlines(crs=datacrs, draw_labels=True, **kw_grid)
        gl.top_labels = False
        gl.left_labels = True
        gl.right_labels = False
        gl.bottom_labels = True
        gl.xlocator = mticker.FixedLocator(xticks)
        gl.ylocator = mticker.FixedLocator(yticks)
        gl.xformatter = LONGITUDE_FORMATTER
        gl.yformatter = LATITUDE_FORMATTER
        gl.xlabel_style = kw_ticklabels
        gl.ylabel_style = kw_ticklabels
        
        # apply tick parameters
        ax.set_xticks(xticks, crs=datacrs)
        ax.set_yticks(yticks, crs=datacrs)
        plt.yticks(color='w') # hack: make the ytick labels white so the ticks show up but not the labels
        plt.xticks(color='w') # hack: make the ytick labels white so the ticks show up but not the labels
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
    
    if (extent is None) or (domain_name=='Global'):        
        ax.set_global()
        extent = [-180., 180., -90., 90.]
    # If extent is given, set map extent to lat/lon bounding box
    else:
        extent2=copy.deepcopy(extent)
        #add some buffer so the tick labels are shown
        extent2[0]=extent[0]-0.001
        extent2[1]=extent[1]+0.001
        extent2[2]=extent[2]-0.001
        extent2[3]=extent[3]+0.001
        
        ax.set_extent(extent2, crs=datacrs)
    if logofile is not None:
        ## Add logo & attributions
        logo = Image.open(logofile)
        ax_width = np.abs(domains.extent[domain_name]["ext"][0] - domains.extent[domain_name]["ext"][1])
        ax_height = np.abs(domains.extent[domain_name]["ext"][2] - domains.extent[domain_name]["ext"][3])
        ax_width_dec = 1
        ax_height_dec = 1/(ax_width/ax_height)
        
        png_width, png_height = logo.size
        logo_png_proportion = png_width/png_height #width/height
        
        logo_width_prop = (9/100)/ax_width_dec #10% of it
        logo_height_prop = (9/100)/ax_height_dec/logo_png_proportion #10% of it
       
        
        logo_box_proportion = [a / b for a, b in zip([1,1], list(domains.extent[domain_name]["figsize"]))]
        
        if logopos =="upper_left":
            logo_box = [0.01, 0.99-logo_box_proportion[1],logo_box_proportion[0],logo_box_proportion[1]]
        elif logopos =="upper_right":
            logo_box = [0.99-logo_box_proportion[0],0.99-logo_box_proportion[1],logo_box_proportion[0],logo_box_proportion[1]]
        elif logopos == 'bottom_left':
            logo_box =[0.01,0.01,logo_box_proportion[0],logo_box_proportion[1]]
        elif logopos == 'bottom_right':
            logo_box = [0.99-logo_box_proportion[0],0.01,logo_box_proportion[0],logo_box_proportion[1]]

        logo_ax = ax.inset_axes(logo_box, zorder=1000)


        logo_ax.spines['bottom'].set_visible(False)
        logo_ax.spines['left'].set_visible(False)
        logo_ax.spines['right'].set_visible(False)
        logo_ax.set_xticks([])
        logo_ax.set_yticks([])
        logo_ax.imshow(logo,aspect='auto')
        logo_ax.axis('off')
        #logo_ax.set_facecolor('white')

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



def custom_tick_formatter(x, pos):

    if x.is_integer():
        return f'{int(x)}'
    else:
        return x

def plot_fields(model_data, config):
    output_fname = config["output_fname"]
    filled_var = config["filled_var"]
    contour_var = config.get("contour_var")
    vector_vars = config.get("vector_vars")
    domain_dash_line = config["domain_dash_line"]
    required_domains = config["required_domains"]
    ulc_title = config["ulc_title"]
    llc_title = config["llc_title"]
    lrc_title = config["lrc_title"]
    logo_file = config["logo_file"]

    base_dpi = 100
    scaling_factor = (current_dpi / base_dpi)**0.5
    set_cw3e_font(current_dpi, scaling_factor, font_file=config["font_file"])

    for domain_name in required_domains:
        cropped_ds = crop_to_domain(model_data, domain_name)

        gs = GridSpec(1, 2, width_ratios=[1, 0.05], wspace=0.01)
        fig = plt.figure(figsize=domains.extent[domain_name]["figsize"], dpi=current_dpi)

        ax = fig.add_subplot(gs[0, 0], projection=domains.extent[domain_name]['ccrs'])
        ax = draw_basemap(ax=ax, domain_name=domain_name,
                          logofile=logo_file, datacrs=cropped_ds.attrs["datacrs"],
                          mask_ocean=False, coastline=True)

        # Colormap
        colorbar_name = f"{filled_var}_{domain_name}"
        cmap, norm, bnds, cbarticks, cbarlbl = ccmaps.cmap(colorbar_name)

        # Filled contour
        cf = ax.contourf(cropped_ds[filled_var].longitude,
                         cropped_ds[filled_var].latitude,
                         cropped_ds[filled_var],
                         transform=cropped_ds.attrs["datacrs"],
                         levels=bnds, cmap=cmap, norm=norm, alpha=0.9)

        # Optional black contours of same field
        ax.contour(cropped_ds[filled_var].longitude,
                   cropped_ds[filled_var].latitude,
                   cropped_ds[filled_var],
                   transform=cropped_ds.attrs["datacrs"],
                   levels=bnds, colors='black', linewidths=0.3*scaling_factor, alpha=0.9)

        # Optional dashed domain
        if domain_dash_line:
            if cropped_ds[filled_var].latitude.ndim == 1:
                lon, lat = np.meshgrid(cropped_ds[filled_var].longitude,
                                       cropped_ds[filled_var].latitude)
            else:
                lat = cropped_ds[filled_var].latitude
                lon = cropped_ds[filled_var].longitude
            lat_perimeter = list(lat[:,0]) + list(lat[-1,:]) + list(lat[:,-1])[::-1] + list(lat[0,:])[::-1]
            lon_perimeter = list(lon[:,0]) + list(lon[-1,:]) + list(lon[:,-1])[::-1] + list(lon[0,:])[::-1]
            ax.plot(lon_perimeter, lat_perimeter, 'k--', transform=cropped_ds.attrs["datacrs"])

        # Optional vectors
        if vector_vars != [None, None]:
            uvar, vvar = vector_vars
            plot_vectors(ax, cropped_ds, uvar, vvar, filled_var, domain_name,
                         scaling_factor, bnds, config)

        # Optional contours
        if contour_var is not None:
            plot_contours(ax, cropped_ds, contour_var, scaling_factor)

        # Titles and colorbar
        pytimeinit = cropped_ds.attrs["init"].astype('datetime64[s]').astype(datetime)
        pytimevalidtime = cropped_ds.attrs["valid_time"].astype('datetime64[s]').astype(datetime)
        leadtime = int((pytimevalidtime - pytimeinit).total_seconds()/3600)

        cbax = plt.subplot(gs[0, 1])
        cb = Colorbar(ax=cbax, mappable=cf, orientation='vertical', ticks=cbarticks)
        cb.ax.yaxis.set_major_formatter(FuncFormatter(custom_tick_formatter))

        ax.set_title(ulc_title.format(model_name=cropped_ds.attrs["model"]) + "\n" +
                     llc_title.format(init_time=datetime.strftime(pytimeinit, '%H UTC %d %b %Y')),
                     loc="left", pad=0.5)
        ax.set_title(lrc_title.format(lead_time=leadtime,
                                      valid_time=datetime.strftime(pytimevalidtime, '%H UTC %d %b %Y')),
                     loc="right", pad=0.5)

        # Save
        fname = output_fname.format(model_name=cropped_ds.attrs["model"],
                                    domain_name=domain_name, lead_time=leadtime)
        create_folder(fname)
        fig.savefig(fname, bbox_inches='tight', dpi=fig.dpi)
        plt.close()
        print(f"[INFO] file saved: {fname}")
