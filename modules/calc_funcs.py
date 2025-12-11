#!/usr/bin/python3
"""
Filename:    calc_funcs.py
Author:      Deanna Nash, dnash@ucsd.edu
Description: functions to calculate variables needed for the ivt cross sections
"""

import sys, os
import xarray as xr
import numpy as np
import pandas as pd

def format_timedelta_to_HHMMSS(td):
    td_in_seconds = td.total_seconds()
    hours, remainder = divmod(td_in_seconds, 3600)
    minutes, seconds = divmod(remainder, 60)
    hours = int(hours)
    minutes = int(minutes)
    seconds = int(seconds)
    if minutes < 10:
        minutes = "0{}".format(minutes)
    if seconds < 10:
        seconds = "0{}".format(seconds)
    return "{}:{}:{}".format(hours, minutes,seconds)

def calc_transport_manual(ds, include_condensates=False):
    '''
    Calculates vertically integrated vapor or water transport (IVT or IWT).
    
    Parameters:
    -----------
    ds : xarray.Dataset
        Must include variables: q, u, v
        Optional: rwmr, clwmr, icmr, snmr (if include_condensates=True)
    
    include_condensates : bool
        If True, includes condensate mixing ratios in the total mixing ratio.

    Returns:
    --------
    xarray.Dataset with ivtu, ivtv, and magnitude (ivt or iwt)
    '''

    pressure = ds.isobaricInhPa.values * 100  # convert hPa to Pa
    dp = np.diff(pressure)
    g = 9.81  # gravity
    ## choose the name for Q if include condensates is true
    Q_name = 'twmr' if include_condensates else 'q'

    qu_lst, qv_lst = [], []

    # enumerate through pressure levels so we select the layers
    for i, pres in enumerate(ds.isobaricInhPa.values[:-1]):
        pres2 = ds.isobaricInhPa.values[i+1]
        tmp = ds.sel(isobaricInhPa=[pres, pres2]) # select layer
        tmp = tmp.mean(dim='isobaricInhPa', skipna=True) # average q, u, v in layer
        # calculate ivtu in layer
        qu = ((tmp[Q_name]*tmp.u*dp[i])/g)*-1
        qu_lst.append(qu)
        # calculate ivtv in layer
        qv = ((tmp[Q_name]*tmp.v*dp[i])/g)*-1
        qv_lst.append(qv)

    ## add up u component of ivt from each layer
    qu = xr.concat(qu_lst, pd.Index(pressure[:-1], name="pres"))
    qu = qu.sum('pres')
    qu.name = 'ivtu' if not include_condensates else 'iwtu'
    
    # ## add up v component of ivt from each layer
    qv = xr.concat(qv_lst, pd.Index(pressure[:-1], name="pres"))
    qv = qv.sum('pres')
    qv.name = 'ivtv' if not include_condensates else 'iwtv'
    
    ## calculate IVT magnitude
    ivt = np.sqrt(qu**2 + qv**2)
    ivt.name = 'ivt' if not include_condensates else 'iwt'

    ds = xr.merge([qu, qv, ivt], compat='no_conflicts')
    
    return ds