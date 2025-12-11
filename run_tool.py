import time
start_time = time.time()  # Record the start time

import sys
import os
import glob
import re
import xarray as xr
import numpy as np
import pandas as pd
import gc
import cartopy.crs as ccrs
import multiprocessing as mp
import gc
import warnings
warnings.filterwarnings("ignore", message="Engine.*loading failed")

sys.path.append('modules')
from read_deterministic_data import calc_gfs_data
from plotter import plot_fields
from calc_funcs import format_timedelta_to_HHMMSS

model_name = sys.argv[1]
fdate = sys.argv[2]
F_lst = np.arange(0, 168+3, 3)

def plot_ICT_IWV(F, fdate):
    print('Reading {0} ...'.format(F))
    ds = calc_gfs_data(F, fdate)

    # Plot - use pconfig to change vars plotted
    outpath="/data/projects/website/mirror/htdocs/Projects/CO_landfalling_ARs/images/"
    pconfig = {
        'IVT': {
            'cfkey': 'ivt',
            'ckey': None,
            'ukey': 'ivtu',
            'vkey': 'ivtv',
            'title': "NCEP GFS IVT (kg m-1 s-1; shaded and vectors)"
        },
        'IWT': {
            'cfkey': 'iwt',
            'ckey': None,
            'ukey': 'ivtu',
            'vkey': 'ivtv',
            'title': "NCEP GFS IWT (kg m-1 s-1; shaded and vectors)"
        },
        'ICT': {
            'cfkey': 'ict',
            'ckey': None,
            'ukey': None,
            'vkey': None,
            'title': "NCEP GFS ICT (kg m-1 s-1; shaded)"
        },
        'RATIO': {
            'cfkey': 'ratio',
            'ckey': None,
            'ukey': None,
            'vkey': None,
            'title': "NCEP GFS ICT/IVT Ratio (%; shaded)"
        }
    }

    for varname in ['IVT', 'IWT', 'ICT', 'RATIO']:
        print(f'Plotting {varname}')

        config = {
                "model_name":"GFS",
                "output_fname" : outpath+'{model_name}_'+pconfig[varname]['cfkey']+'_{domain_name}_latest_F{lead_time}.png',
                "filled_var" : pconfig[varname]['cfkey'],
                "contour_var" : pconfig[varname]['ckey'],
                "vector_vars" : [pconfig[varname]['ukey'],pconfig[varname]['vkey']],
                "domain_dash_line" : False,
                "required_domains" :  ['IntWest'],
                "ulc_title":pconfig[varname]['title'],
                "llc_title":"Initialized: {init_time}",
                "lrc_title":"F-{lead_time:03d}  Valid: {valid_time}",
                "logo_file":'/home/dnash/repos/integrated_cloud_water_transport_calculation/CW3E-Logo-Vertical-Acronym-FullColor.png',
                "font_file":'/home/dnash/repos/integrated_cloud_water_transport_calculation/modules/helvetica.ttc'
            }


        plot_fields(ds, config)

    ds.close()
    del ds
    gc.collect()
    
def multiP_create_plots(F):
    plot_ICT_IWV(F, fdate)
    
if __name__ == '__main__':
    # sys.argv.append(None) ## add this in case date not specified in command line
    # fdate = sys.argv[1] ## set this to None to get most recently downloaded data    
    
    start_time = pd.Timestamp.today()
    print(f'Creating ICT plots for {model_name}')
    
    ###################################
    ### CREATE PLOTS USING POOL.MAP ###
    ###################################

    print('...create plots ...')

    with mp.get_context('spawn').Pool(processes=16) as pool:
        # debug_print("Via map with exception")
        # debug_print("\tKicking off pool via map with exception")
        pool.map(multiP_create_plots,F_lst)
        # debug_print("\tKicked off pool via map with exception")
        pool.close()
        pool.join()


    end_time = pd.Timestamp.today()
    td = end_time - start_time
    td = format_timedelta_to_HHMMSS(td)
    print(f'Plots for {model_name} took {td} to run')

