"""
Filename:    domains.py
Author:      Deanna Nash, dnash@ucsd.edu
Description: Dictionaries of domains for CW3E web maps
"""

## each domain can be used as extent argument for plotting: ext = [lonmin, lonmax, latmin, latmax] 
## Example: 
# "name_of_domain":  {'ext': [eastern most longitude, western most longitude, southern most latitude, northern-most latitude], # note that the extent should be in 0 to 360 longitude (unless it crosses the prime meridian, then it should be in -180 to 180 and -90 to 90 latitude
#                         'xticks': [locations of longitude label ticks]
#                         'yticks': [locations of latitude lable ticks]
#                         'center_longitude': the central longitude of the map
#                         'figsize': the figure size in inches, (width, height)
#                        }

import cartopy.crs as ccrs

extent = {"NPac":      {'ext': [115., 245., -10., 70.], # North Pacific (extent in 0-360)
                            'xticks': [125, 140, 155, 170, -175, -160, -145, -130, -115],
                            'yticks': [-10, 0, 10, 20, 30, 40, 50, 60, 70],
                            'ccrs': ccrs.PlateCarree(central_longitude=180),
                            'figsize': (10. ,5.5)
                           },
          "wwrf_npacific":  {'ext': [170., 260., 10., 65.], # North Pacific (extent in 0-360)
                            'xticks': [170, -175, -160, -145, -130, -115, -100],
                            'yticks': [10, 20, 30, 40, 50, 60, 70],
                            'ccrs': ccrs.PlateCarree(central_longitude=180),
                            'figsize': (10. ,5.5)
                           },
            
          
           "NEPac":    {'ext': [180., 245., 10., 70.], # Northeast Pacific (extent in 0-360)
                            'xticks': [-180, -170, -160, -150, -140, -130, -120],
                            'yticks': [10, 20, 30, 40, 50, 60, 70],
                            'ccrs': ccrs.PlateCarree(central_longitude=180),
                            'figsize': (10. ,8.)
                            },  
          
           "USWC":         {'ext': [220., 250., 20., 50.], # U.S. West Coast (extent in 0-360)
                            'xticks': [-140, -135, -130, -125, -120, -115, -110],
                            'yticks': [20, 25, 30, 35, 40, 45, 50],
                            'ccrs': ccrs.PlateCarree(central_longitude=180),
                            'figsize': (10. ,8.5)
                           },  
          
           "NAmerica":     {'ext': [220, 300, 10., 60.], # North America
                            'xticks': [-140, -130, -120, -110, -100, -90, -80, -70, -60],
                            'yticks': [10, 20, 30, 40, 50, 60],
                            'ccrs': ccrs.PlateCarree(central_longitude=0),
                            'figsize': (10. ,5.5)
                           }, 
          
           "NAtlantic":    {'ext': [-90., 15., 10., 70.], # North Atlantic (ext in -180 to 180)
                            'xticks': [-90, -75, -60, -45, -30, -15, 0, 15],
                            'yticks': [10, 20, 30, 40, 50, 60, 70],
                            'ccrs': ccrs.PlateCarree(central_longitude=0),
                            'figsize': (10. ,5.2)
                           }, 
          
           "IntWest":      {'ext': [230., 260., 20., 50.], # Interior West
                            'xticks': [-130, -125, -120, -115, -110, -105, -100],
                            'yticks': [20, 25, 30, 35, 40, 45, 50],
                            'ccrs': ccrs.PlateCarree(central_longitude=180),
                            'figsize': (10. ,8.5)
                           },
          
           "Antarctica":   {'ext': [-180., 180, -90, -50.], # Antarctica
                            'xticks': [-15, -30,-45, -60, -75, -90, -105, -120, -135, -150, -165, 180, 165, 150, 135, 120, 105, 90, 75, 60, 45, 30, 15, 0],
                            'yticks': [-50, -60, -70, -80, -90],
                            'ccrs': ccrs.SouthPolarStereo(central_longitude=0),
                            'figsize': (10. ,8.5)
                            }, 
            "Arctic":   {'ext': [-180., 180, 50, 90.], # Arctic
                            'xticks': [-15, -30,-45, -60, -75, -90, -105, -120, -135, -150, -165, 180, 165, 150, 135, 120, 105, 90, 75, 60, 45, 30, 15, 0],
                            'yticks': [50, 60, 70, 80, 90],
                            'ccrs': ccrs.NorthPolarStereo(central_longitude=0),
                            'figsize': (10. ,8.5)
                           }, 
           
            "EUS":         {'ext': [260, 295, 20., 50.], # southeast us
                            'xticks': [-100, -95, -90, -85, -80, -75, -70, -65],
                            'yticks': [20, 25, 30, 35, 40, 45, 50],
                            'ccrs': ccrs.PlateCarree(central_longitude=0),
                            'figsize': (10. ,7.5)},
            
            "WCQPF":        {'ext': [235, 255, 30., 50.], # west coast qpf
                            'xticks': [-125, -120, -115, -110, -105],
                            'yticks': [30, 35, 40, 45, 50],
                            'ccrs': ccrs.PlateCarree(central_longitude=0),
                            'figsize': (10. ,7.5)},
            "WCQPF2":       {'ext': [235, 245, 30., 50.], # west coast qpf 2
                            'xticks': [-125, -120, -115, -110, -105],
                            'yticks': [30, 35, 40, 45, 50],
                            'ccrs': ccrs.PlateCarree(central_longitude=0),
                            'figsize': (8. ,12)},
            "WCQPF3":        {'ext': [220, 255, 30., 60.], # west coast qpf + SE Alaska
                            'xticks': [-140, -135, -130, -125, -120, -115, -110, -105],
                            'yticks': [30, 35, 40, 45, 50, 55, 60],
                            'ccrs': ccrs.PlateCarree(central_longitude=0),
                            'figsize': (10. ,7)},
                            
            "SFBAY":        {'ext': [236, 239, 37., 40.], # west coast qpf 2
                            'xticks': [-124, -123, -122, -121],
                            'yticks': [37, 38, 39, 40],
                            'ccrs': ccrs.PlateCarree(central_longitude=0),
                            'figsize': (10. ,7.5)},
            "NH":           {'ext': [120, 380, 0., 80.], # northern hemisphere
                            'xticks': [120, 140, 160, 180, -160, -140, -120, -100, -80, -60, -40, -20, 0, 20],
                            'yticks': [0, 10, 20, 30, 40, 50, 60, 70, 80],
                            'ccrs': ccrs.PlateCarree(central_longitude=-110),
                            'figsize': (18., 5.)},
            "Global":       {'ext': [0, 360, -90., 90.], # northern hemisphere
                            'xticks': [120, 140, 160, 180, -160, -140, -120, -100, -80, -60, -40, -20, 0, 20],
                            'yticks': [-80,-60,-40,-20, 0, 20, 40, 60, 80],
                            'ccrs': ccrs.PlateCarree(central_longitude=180),
                            'figsize': (10., 5.)},
                            

                            
          }
           
           