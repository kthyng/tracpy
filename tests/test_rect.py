#!/usr/bin/env python

"""
Some basic tests of tracpy, using the simple rectangle input
"""

import tracpy

import os
import time
import datetime
import netCDF4

# some simple example data
currents_filename = os.path.join('input', 'ocean_his_0001.nc')
grid_filename = os.path.join('input', 'grid.nc')
time_units = 'seconds since 1970-01-01'
num_layers = 3


def test_run_2d():
    """
    can we initialize tracpy
    """

    name = 'test_run_2d'

    grid_nc = netCDF4.Dataset(currents_filename)
    
    loc = [currents_filename, grid_filename]
    
    start = time.time()
    grd = tracpy.inout.readgrid(loc,
                                nc=grid_nc)
    print "building grid took:", time.time() - start

    # Start date in date time formatting
    date = datetime.datetime(2013, 12, 19, 0)

    # Time between outputs
    tseas = 4*3600. # 4 hours between outputs, in seconds 

    # Number of days to run the drifters.
    ndays = tseas*9./(3600.*24)

    # Sets a smaller limit than between model outputs for when to force interpolation if hasn't already occurred.
    nsteps = 5

    # Controls the sampling frequency of the drifter tracks.
    N = 4

    # Use ff = 1 for forward in time and ff = -1 for backward in time.
    ff = 1

    ah = 0. # m^2/s
    av = 0. # m^2/s

    # turbulence/diffusion flag
    doturb = 0

    # two particles
    lon0 = [-123., -123.]
    lat0 = [48.55, 48.75]

    do3d = 0 # flag to set to 2-d

    z0 = 's' #'z' #'salt' #'s' 
    zpar = num_layers-1 # top layer

    lonp, latp, zp, t, grd = tracpy.run.run(loc,
                                            nsteps,
                                            ndays,
                                            ff,
                                            date,
                                            tseas,
                                            ah,
                                            av,
                                            lon0,
                                            lat0, 
                                            z0,
                                            zpar,
                                            do3d,
                                            doturb,
                                            name,
                                            grid=grd,
                                            dostream=0,
                                            N=N)
    print lonp
    print latp


    assert False






