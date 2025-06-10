"""
Filename:    cw3ecmaps.py
Author:      Deanna Nash, dnash@ucsd.edu
Description: Functions for cmaps from the CW3E website adapted from from https://github.com/samwisehawkins/nclcmaps
"""

import numpy as np
import matplotlib.colors as mcolors
from matplotlib.colors import ListedColormap
from matplotlib.ticker import FuncFormatter

__all__ = ['cw3ecmaps']


cw3e_cmaps =    {
                "ivt" :{ 
                        "colors":[[255, 255, 3],
                                [255, 229, 3],
                                [255, 201, 2],
                                [255, 175, 2],
                                [255, 131, 1],
                                [255, 79, 1], 
                                [255, 24, 1],  
                                [235, 1, 7],  
                                [185, 0, 55], 
                                [134, 0, 99],
                                [86, 0, 137]],
                        "bounds":[250, 300, 400, 500, 600, 700, 800, 1000, 1200, 1400, 1600, 5000],
                        "ticks":[250, 300, 400, 500, 600, 700, 800, 1000, 1200, 1400, 1600],
                        "label": "IVT (kg m$^{-1}$ s$^{-1}$)",
                        },
                "ivt_intwest":{
                        "colors":[[12, 255, 255], # 100-150
                                [16, 185, 6], # 150-200
                                [140, 220, 2], # 200-250
                                [255, 255, 3], # 250-300
                                [255, 229, 3], # 300-400
                                [255, 201, 2], # 400-500
                                [255, 175, 2], # 500-600
                                [255, 131, 1], # 600-700
                                [255, 79, 1],  # 700-800
                                [255, 24, 1],  # 800-1000
                                [235, 1, 7],   # 1000-1200
                                [185, 0, 55],  # 1200-1400
                                [134, 0, 99], # 1400-1600
                                [86, 0, 137]],   # 1600+,
                        "bounds":[100, 150, 200, 250, 300, 400, 500, 600, 700, 800, 1000, 1200, 1400, 1600, 5000],
                        "ticks":[100, 150, 200, 250, 300, 400, 500, 600, 700, 800, 1000, 1200, 1400, 1600],
                        "label":'IVT (kg m$^{-1}$ s$^{-1}$)',
                        },
                }
                  
def cmap(cbarname):
        data = np.array(cw3e_cmaps[cbarname]["colors"])
        data = data / np.max(data)
        cmap = ListedColormap(data, name=cbarname)
        bnds = cw3e_cmaps[cbarname]["bounds"]
        norm = mcolors.BoundaryNorm(bnds, cmap.N)
        cbarticks = cw3e_cmaps[cbarname]["ticks"]
        cbarlbl = cw3e_cmaps[cbarname]["label"]
        return cmap, norm, bnds, cbarticks, cbarlbl