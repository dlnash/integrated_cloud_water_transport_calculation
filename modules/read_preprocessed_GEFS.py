"""
Filename:    read_preprocessed_GEFS.py
Author:      Deanna Nash, dnash@ucsd.edu
Description: Functions for reading preprocessed GEFS data
"""

import xarray as xr
import numpy as np

def load_archive_GEFS_forecast(date, varname, F=None):
    fpath = '/expanse/nfs/cw3e/cwp186/data/preprocessed/GEFS/'

    ### load forecast from GEFS

    if F == None:
        fname_pattern = fpath + '{0}.t00z.0p50.f*.{1}'.format(date, varname)
        forecast = xr.open_mfdataset(fname_pattern, engine='netcdf4', concat_dim="step", combine='nested')
    else:
        F = str(F).zfill(3)
        fname = fpath + '{0}.t00z.0p50.f{2}.{1}'.format(date, varname, F)
        forecast = xr.open_dataset(fname)
    
    forecast = forecast.rename({'longitude': 'lon', 'latitude': 'lat', 
                                  "time": "init_date"}) # need to rename this to match GEFSv12 Reforecast

    return forecast