import os
import copy
import numpy as np
import xarray as xr
from scipy.ndimage import gaussian_filter
import cartopy.crs as ccrs
import matplotlib.pyplot as plt
from matplotlib.colorbar import Colorbar
from matplotlib.ticker import FuncFormatter
import domains as domains

def create_folder(filename):
    '''
    :param path: path string
    :return:
    '''

    if not os.path.exists(os.path.dirname(filename)):
        os.makedirs(os.path.dirname(filename))

# ----------------------------------------------------------------------
# Domain cropping
# ----------------------------------------------------------------------
def crop_to_domain(model_data, domain_name):
    """
    Crop a dataset to the requested domain.
    """
    lat_min = domains.extent[domain_name]['ext'][2]
    lat_max = domains.extent[domain_name]['ext'][3]
    lon_min = domains.extent[domain_name]['ext'][0]
    lon_max = domains.extent[domain_name]['ext'][1]

    if domain_name == 'Antarctica':
        lat_max = 0
    if domain_name == 'Arctic':
        lat_min = 0

    if domain_name in ['NAtlantic', 'NH']:
        return copy.deepcopy(model_data)
    elif domain_name in ['Antarctica', 'Arctic']:
        if model_data.latitude[0] > model_data.latitude[-1]:
            return model_data.sel(latitude=slice(lat_max, lat_min))
        else:
            return model_data.sel(latitude=slice(lat_min, lat_max))
    else:
        if len(model_data.latitude.shape) == 2:
            return copy.deepcopy(model_data)
        else:
            if model_data.latitude[0] > model_data.latitude[-1]:
                return model_data.sel(latitude=slice(lat_max, lat_min),
                                      longitude=slice(lon_min, lon_max))
            else:
                return model_data.sel(latitude=slice(lat_min, lat_max),
                                      longitude=slice(lon_min, lon_max))

# ----------------------------------------------------------------------
# Vector plotting
# ----------------------------------------------------------------------
def plot_vectors(ax, cropped_ds, uvar, vvar, ref_var, domain_name, scaling_factor, bnds, config):
    """
    Plot quiver vectors on map axes with spacing rules.
    """
    u = cropped_ds[uvar].values.copy()
    v = cropped_ds[vvar].values.copy()

    image_dx = np.abs(domains.extent[domain_name]['ext'][1] -
                      domains.extent[domain_name]['ext'][0])
    datares = config.get('datares', np.abs(cropped_ds[ref_var].longitude.values.flatten()[1] -
                                           cropped_ds[ref_var].longitude.values.flatten()[0])/111)
    arrows_spacing_in_gridpoints = int((image_dx/datares)/35)

    # Mesh
    if cropped_ds[ref_var].latitude.values.ndim == 2:
        lats = cropped_ds[ref_var].latitude.values
        lons = cropped_ds[ref_var].longitude.values
    else:
        lons, lats = np.meshgrid(cropped_ds[ref_var].longitude.values,
                                 cropped_ds[ref_var].latitude.values)

    # Polar exception
    if isinstance(ax.projection, (ccrs.SouthPolarStereo, ccrs.NorthPolarStereo)):
        proj_geo = ccrs.PlateCarree()
        xy_proj = ax.projection.transform_points(proj_geo, lons, lats)
        x, y = xy_proj[..., 0], xy_proj[..., 1]

        target_spacing_m = 300 * 1000
        x_bins = np.arange(np.nanmin(x), np.nanmax(x)+target_spacing_m, target_spacing_m)
        y_bins = np.arange(np.nanmin(y), np.nanmax(y)+target_spacing_m, target_spacing_m)
        x_idx = np.digitize(x.ravel(), x_bins) - 1
        y_idx = np.digitize(y.ravel(), y_bins) - 1

        _, unique_indices = np.unique(np.stack([x_idx, y_idx], axis=1), axis=0, return_index=True)
        arrows_longitude = lons.ravel()[unique_indices]
        arrows_latitude  = lats.ravel()[unique_indices]
        arrows_u = u.ravel()[unique_indices]
        arrows_v = v.ravel()[unique_indices]
    else:
        if domain_name == 'NH':
            arrows_spacing_in_gridpoints = int((image_dx/datares)/45)
        arrows_latitude  = lats[::arrows_spacing_in_gridpoints, ::arrows_spacing_in_gridpoints]
        arrows_longitude = lons[::arrows_spacing_in_gridpoints, ::arrows_spacing_in_gridpoints]
        arrows_u = u[::arrows_spacing_in_gridpoints, ::arrows_spacing_in_gridpoints]
        arrows_v = v[::arrows_spacing_in_gridpoints, ::arrows_spacing_in_gridpoints]

    # Mask low values
    low_values = np.where(cropped_ds[ref_var].values < bnds[0])
    u[low_values] = np.nan
    v[low_values] = np.nan

    # Plot quiver
    q = ax.quiver(arrows_longitude, arrows_latitude, arrows_u, arrows_v,
                  scale=13/scaling_factor, scale_units='dots',
                  transform=ccrs.PlateCarree(),
                  width=0.0009*scaling_factor)

    # Quiver key inset
    quiver_box_proportion = [a / b for a, b in zip([1, 1], list(domains.extent[domain_name]["figsize"]))]
    if domain_name == 'NH':
        quiver_box = [1-quiver_box_proportion[0], 0,
                      quiver_box_proportion[0], quiver_box_proportion[1]]
    else:
        quiver_box = [1-quiver_box_proportion[0], 1-quiver_box_proportion[1],
                      quiver_box_proportion[0], quiver_box_proportion[1]]

    ax2 = ax.inset_axes(quiver_box)
    ax2.set_xticks([]); ax2.set_yticks([])
    ax2.quiverkey(Q=q,
                  X=quiver_box[0] + quiver_box_proportion[0]/2,
                  Y=quiver_box[1] + quiver_box_proportion[1]/1.3,
                  U=500, label='\n500\n[kg m$^{-1}$ s$^{-1}$]',
                  labelpos='S', color="k",
                  fontproperties={"size": 5*scaling_factor},
                  coordinates='axes')

# ----------------------------------------------------------------------
# Contour overlay
# ----------------------------------------------------------------------
def plot_contours(ax, cropped_ds, varname, scaling_factor, levels=np.arange(3000,3280,40)):
    """
    Plot smoothed contour overlays (e.g., MSLP).
    """
    smooth_field = gaussian_filter(cropped_ds[varname].values, sigma=1.0)
    contour_lines = ax.contour(cropped_ds[varname].longitude,
                               cropped_ds[varname].latitude,
                               smooth_field,
                               transform=cropped_ds.attrs["datacrs"],
                               colors='k',
                               levels=levels)

    kw_clabels = {'fontsize': 7*scaling_factor, 'inline': True,
                  'inline_spacing': 5, 'fmt': '%i',
                  'rightside_up': True, 'use_clabeltext': True}

    cl = ax.clabel(contour_lines, contour_lines.levels[::2], **kw_clabels)
    for txt in cl:
        txt.set_bbox(dict(facecolor='white', edgecolor='none', pad=0.5))
