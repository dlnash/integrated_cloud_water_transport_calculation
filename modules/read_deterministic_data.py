#!/usr/bin/python3
"""
Filename:    read_deterministic_data.py
Author:      Deanna Nash, dnash@ucsd.edu
Description: functions to read deterministic data from GFS and ECMWF
"""

import sys
import os
import glob
import shutil
import subprocess
import re
import xarray as xr
import numpy as np
import pandas as pd
import datetime
import cartopy.crs as ccrs
import calc_funcs as cfuncs
import gc

def read_gfs_deterministic(filename, vardict, show_catalog=False):

    '''
    author: Ricardo Vilela
    email: rbatistavilela@ucsd.edu

    function usage:

    filename example:
    gfs_2024041512_f003.grb

    vardict example:
    vardict = {
            "u_wind":{"typeOfLevel":'isobaricInhPa',"shortName":"u"}, #U-component of wind
            "v_wind":{"typeOfLevel":'isobaricInhPa',"shortName":"v"}, #V-component of wind
            "iwv":{"typeOfLevel":'atmosphereSingleLayer',"shortName":"pwat"}, #Integrated precipitable water or IWV
            "temperature":{"typeOfLevel":'isobaricInhPa',"shortName":"t"}, #Temperature
            "rh":{"typeOfLevel":'isobaricInhPa',"shortName":"r"}, #Relative Humidity
            "sfc_pressure":{'name': 'Surface pressure', 'typeOfLevel': 'surface', 'level': 0, 'paramId': 134, 'shortName': 'sp'}, #surface pressure
            "sea_level_pressure":{'name': 'Pressure reduced to MSL', 'typeOfLevel': 'meanSea', 'level': 0, 'paramId': 260074, 'shortName': 'prmsl'} #mean sea level pressure
        }
    
    Output:
    Dictionary of objects for each variable set in the vardict argument. Ex.
    uwind = selected_vars["u_wind"].values
    uwind_latitude = selected_vars["u_wind"].latitude.values
    uwind_longitude = selected_vars["u_wind"].longitude.values

    for 3d variables:
    uwind_pressure_levels = selected_vars["u_wind"].isobaricInhPa.values

    '''
        
    selected_vars = {}
    #iterating over all variables in the vardict and storing each one in a new dictionary called selected_vars
    for var in vardict.keys():

        ds = xr.open_dataset(filename,decode_timedelta=True,
                        engine='cfgrib',filter_by_keys=vardict[var],backend_kwargs={"indexpath": ''})
        
        
        ## subset to specified pressure levels and points
        if vardict[var]["typeOfLevel"] =='isobaricInhPa':
            ds = ds.sel(isobaricInhPa=slice(1000, 200))
            
        else:
            ds = ds
            
        #exceptions when shortName does not match variable name in the grib    
        if var == 'u_wind_10m':
            selected_vars[var] = ds["u10"]
        elif var == 'v_wind_10m':
            selected_vars[var] = ds["v10"]
        else:    
            selected_vars[var] = ds[vardict[var]["shortName"]]    

    return selected_vars

def read_ecmwf_S2D(filename, vardict, show_catalog=False):
    
    '''
    author: Ricardo Vilela
    email: rbatistavilela@ucsd.edu

    function usage:
    
    filename example:
    S2D04151200041515001.grb

    vardict example:

    vardict = {
               
                "u_wind":{"shortName":'u'}, #U-component of wind
                "v_wind":{"shortName":'v'}, #V-component of wind
                "temperature":{"shortName":'t'}, #Temperature
                "specific_humidity":{"shortName":'q'}, #Specific humidity
                
                
                }
    Output:
    Dictionary of objects for each variable set in the vardict argument. Ex.
    uwind = selected_vars["u_wind"].values
    uwind_latitude = selected_vars["u_wind"].latitude.values
    uwind_longitude = selected_vars["u_wind"].longitude.values

    These files are in model levels instead of pressure levels, in other words, pressure is not a dimension but an extra variable built based on a lookup table of coefficients and the hybrid dimension.

    retrieving the pressure array
    pressure_array = selected_vars["pressure"].values


    '''
    
    print('[INFO] reading ECMWF file: '+filename)    


    selected_vars = {}
    
    # reads surface pressure first (this field is needed in order to calculate the pressure levels for each grid)
    sfc_pressure_ds = xr.open_dataset(filename,
                        engine='cfgrib',filter_by_keys={"shortName":'lnsp'},backend_kwargs={"indexpath": ''})
    
    sfc_pressure = 2.71828**sfc_pressure_ds.lnsp.values

    #read the coefficient lookup table
    coeff = pd.read_csv("/data/projects/operations/wvflux_meteograms/utils/ecmwf_coeffs.txt",names=["A","B"],sep=" ")
    coeffA = np.array(coeff["A"])[:, np.newaxis, np.newaxis]
    coeffB = np.array(coeff["B"])[:, np.newaxis, np.newaxis]

    #creating the pressure based on the coefficient table and the surface pressure
    pressure = (coeffA + (sfc_pressure * coeffB))/100

    #iterating over all variables in the vardict and storing each one in a new dictionary called selected_vars
    for var in vardict.keys():
    
        ds = xr.open_dataset(filename,
                        engine='cfgrib',filter_by_keys=vardict[var],backend_kwargs={"indexpath": ''})
        
        selected_vars[var] = ds[vardict[var]["shortName"]]
        
    ngrids = pressure.shape[0]*pressure.shape[1]*pressure.shape[2]
    
    #creating pressure 3d data array and its attributes
    pressure = xr.DataArray(pressure, name="pressure", dims=("hybrid", "latitude","longitude"), coords={"hybrid": coeff.index, "latitude":sfc_pressure_ds.latitude.values, "longitude":sfc_pressure_ds.longitude.values, "time":sfc_pressure_ds.time.values, "valid_time":sfc_pressure_ds.valid_time.values})
    pressure.attrs = {'GRIB_paramId': 0, 'GRIB_dataType': 'fc', 'GRIB_numberOfPoints': ngrids, 'GRIB_typeOfLevel': 'hybrid', 'GRIB_stepUnits': 1, 'GRIB_stepType': 'instant', 'GRIB_gridType': 'regular_ll', 'GRIB_NV': 276, 'GRIB_Nx': len(pressure.longitude.values), 'GRIB_Ny': len(pressure.latitude.values), 'GRIB_Nz': len(pressure.hybrid.values), 'GRIB_cfName': 'unknown', 'GRIB_cfVarName': 'pressure', 'GRIB_gridDefinitionDescription': 'Latitude/longitude', 'GRIB_iDirectionIncrementInDegrees': 0.1, 'GRIB_iScansNegatively': 0, 'GRIB_jDirectionIncrementInDegrees': 0.1, 'GRIB_jPointsAreConsecutive': 0, 'GRIB_jScansPositively': 0, 'GRIB_latitudeOfFirstGridPointInDegrees': np.nanmax(pressure.latitude.values), 'GRIB_latitudeOfLastGridPointInDegrees': np.nanmin(pressure.latitude.values), 'GRIB_longitudeOfFirstGridPointInDegrees': np.nanmin(pressure.longitude.values), 'GRIB_longitudeOfLastGridPointInDegrees': np.nanmax(pressure.longitude.values), 'GRIB_missingValue': 3.4028234663852886e+38, 'GRIB_name': 'Grid wise atmospheric pressure', 'GRIB_shortName': 'pressure', 'GRIB_units': 'Numeric', 'long_name': 'Atmospheric Pressure', 'units': 'hPa', 'standard_name': 'unknown'}
    
    
    #adding pressure variable to the object dictionary
    selected_vars["pressure"] = pressure
    

    #standardize latitude longitude time and level dimension name and position in the array

    return selected_vars


def read_ecmwf_S1D(filename, vardict, show_catalog=False):
    
    '''
    author: Ricardo Vilela
    email: rbatistavilela@ucsd.edu

    function usage:
    
    filename example:
    S1D04151200041515001

    vardict example:

    vardict = {
               
               "msl_pressure":{"shortName":'msl'}, #msl pressure
                
                
                }
    Output:
    Dictionary of objects for each variable set in the vardict argument. Ex.
    mslp = selected_vars["mslp"].values


    '''
    
    print('[INFO] reading ECMWF file: '+filename)
        
    selected_vars = {}

    #iterating over all variables in the vardict and storing each one in a new dictionary called selected_vars
    for var in vardict.keys():
    
        ds = xr.open_dataset(filename,
                        engine='cfgrib',filter_by_keys=vardict[var],backend_kwargs={"indexpath": ''})
        
        if var == 'u_wind_10m':
            selected_vars[var] = ds["u10"]
        elif var == 'v_wind_10m':
            selected_vars[var] = ds["v10"]
        else:    
            selected_vars[var] = ds[vardict[var]["shortName"]]   
        
    return selected_vars

def calc_gfs_data(F, fdate):
    #########################
    ### READ NEW GFS DATA ###
    #########################
    date_string = fdate
    fpath = '/data/projects/external_datasets/GFS/processed/{0}/'.format(date_string)

    fname = '{0}_F{1}.grb2'.format(date_string, str(F).zfill(3))

    fname = fpath+fname

    gfs_vardict = {

                    "rwmr":{'typeOfLevel': 'isobaricInhPa', 'shortName': 'rwmr'}, # rain mixing ratio (kg/kg)
                    "clwmr":{'typeOfLevel': 'isobaricInhPa', 'shortName': 'clwmr'}, # cloud water mixing ratio (kg/kg)
                    "icmr":{'typeOfLevel': 'isobaricInhPa', 'shortName': 'icmr'}, # ice mixing ratio (kg/kg)
                    "snmr":{'typeOfLevel': 'isobaricInhPa', 'shortName': 'snmr'},  # snow mixing ratio (kg/kg)
                    "grle":{'typeOfLevel': 'isobaricInhPa', 'shortName': 'grle'},  # graupel mixing ratio (kg/kg)
                    "u":{'typeOfLevel': 'isobaricInhPa', "shortName":'u'}, #U-component of wind
                    "v":{'typeOfLevel': 'isobaricInhPa', "shortName":'v'}, #V-component of wind
                    "t":{'typeOfLevel': 'isobaricInhPa', "shortName":'t'}, #Temperature
                    "q":{'typeOfLevel': 'isobaricInhPa', "shortName":'q'}, #Specific humidity
                    "sp": {'typeOfLevel': 'surface',  'level': 0, 'paramId': 134, 'shortName': 'sp'} # surface pressure
                    }

    gfs = read_gfs_deterministic(fname, gfs_vardict, show_catalog=False)

    # --- merge all selected vars to single dataset ---
    ds = xr.merge(gfs.values(), compat='no_conflicts', join="outer")

    # --- add condensates together as new var ---
    ds['twmr'] = ds['q'] + ds['rwmr'] + ds['clwmr'] + ds['icmr'] + ds['snmr'] + ds['grle']

    # --- Compute IVT and IWT ---
    ivt_ds = cfuncs.calc_transport_manual(ds, include_condensates=False)
    iwt_ds = cfuncs.calc_transport_manual(ds, include_condensates=True)
    ds = xr.merge([ivt_ds, iwt_ds], compat='no_conflicts')

    # --- trash collection ---
    del gfs, ivt_ds, iwt_ds
    gc.collect()

    # --- Compute ICT ---
    ds['ict'] = ds['iwt'] - ds['ivt']

    # --- Compute ICT/IVT ratio expressed as a percent (%) ---
    ds['ratio'] = (ds['ict']/ds['ivt']) * 100

    # --- Add attributes necessary for plotting function ---
    new_attrs = {
            "model": "GFS",
            "init": ds.time.values,
            "valid_time": ds.valid_time.values,
            "datacrs": ccrs.PlateCarree(central_longitude=0),
        }
    ds.attrs = new_attrs
    
    return ds