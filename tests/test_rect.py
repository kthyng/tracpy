#!/usr/bin/env python

"""
Some basic tests of tracpy, using the simple rectangle input
"""

import tracpy
import tracpy.calcs
from tracpy.tracpy_class import Tracpy
import os
import time
import datetime
import numpy as np
import pyproj


def test_2dtransport():

    pass


def test_run_2d_ll():
    """
    Test initialization and running of tracpy with lat/lon coordinates.

    can we initialize and run tracpy (using rectangle example). Compare final
    location of drifters with known analytic answer. Using lon/lat coords.
    """

    # some simple example data
    currents_filename = os.path.join('input', 'ocean_his_0001.nc')
    grid_filename = os.path.join('input', 'grid.nc')
    time_units = 'seconds since 1970-01-01'
    num_layers = 3

    name = 'test_run_2d_ll'

    # Start date in date time formatting
    date = datetime.datetime(2013, 12, 19, 0)

    # Time between outputs
    tseas = 4*3600.  # 4 hours between outputs, in seconds

    # Number of days to run the drifters.
    ndays = tseas*9./(3600.*24)

    # Sets a smaller limit than between model outputs for when to force
    # interpolation if hasn't already occurred.
    nsteps = 5

    # Controls the sampling frequency of the drifter tracks.
    N = 4

    # This allows the user to call to TRACMASS for a different period of time
    # than between 2 model outputs
    # Just testing to try new loop, should have same behavior as before
    dtFromTracmass = tseas/2.

    # Use ff = 1 for forward in time and ff = -1 for backward in time.
    ff = 1  # will work for ff=1 or ff=-1 since checks by distance traveled

    ah = 0.  # m^2/s
    av = 0.  # m^2/s

    # turbulence/diffusion flag
    doturb = 0

    # two particles (starting positions)
    lon0 = [-123., -123.]
    lat0 = [48.55, 48.75]

    do3d = 0  # flag to set to 2-d

    z0 = 's'  # 'z' 'salt' 's'
    zpar = num_layers-1  # top layer

    usespherical = True

    # Get projection object
    proj = tracpy.tools.make_proj(setup='galveston', usebasemap=False)

    # Read in grid
    grid = tracpy.inout.readgrid(grid_filename, proj,
                                 vert_filename=currents_filename,
                                 usespherical=usespherical)

    # Initialize Tracpy class
    tp = Tracpy(currents_filename, grid, name=name, tseas=tseas,
                ndays=ndays, nsteps=nsteps, N=N, ff=ff, ah=ah, av=av,
                doturb=doturb, do3d=do3d, z0=z0, zpar=zpar,
                time_units=time_units, dtFromTracmass=dtFromTracmass,
                usespherical=usespherical)

    lonp, latp, zp, t, T0, U, V = tracpy.run.run(tp, date, lon0, lat0)

    ## check the results:
    print lonp.shape
    print lonp
    print latp

    # since eastward current, latitude should not change:
    assert np.allclose(lat0, latp.T)

    # current velocity -- 0.1 m/s
    # position
    distance = (ndays * 24 * 3600 * 0.1)*ff

    # better to use pyproj to compute the geodesic
    geod = pyproj.Geod(ellps='clrk66')
    end = geod.fwd(lon0, lat0, (90, 90), (distance, distance), radians=False)

    assert np.allclose(lonp[:, -1], end[0])


def test_run_2d_xy():
    """
    Test initialization and running of tracpy with projected coordinates.

    can we initialize and run tracpy (using rectangle example). Compare final
    location of drifters with known analytic answer. Using x/y coords for
    idealized type runs. Made this simulation information from the other test
    information using:
    """

    # some simple example data
    currents_filename = os.path.join('input', 'ocean_his_0001.nc')
    grid_filename = os.path.join('input', 'gridxy.nc')
    time_units = 'seconds since 1970-01-01'
    num_layers = 3

    name = 'test_run_2d_xy'

    start = time.time()

    # grd = tracpy.inout.readgrid(grid_filename, vert_filename=currents_filename)

    print "building grid took:", time.time() - start

    # Start date in date time formatting
    date = datetime.datetime(2013, 12, 19, 0)

    # Time between outputs
    tseas = 4*3600.  # 4 hours between outputs, in seconds

    # Number of days to run the drifters.
    ndays = tseas*9./(3600.*24)

    # Sets a smaller limit than between model outputs for when to force
    # interpolation if hasn't already occurred.
    nsteps = 5

    # Controls the sampling frequency of the drifter tracks.
    N = 4

    # This allows the user to call to TRACMASS for a different period of time
    # than between 2 model outputs
    # Just testing to try new loop, should have same behavior as before
    dtFromTracmass = tseas/2.

    # Use ff = 1 for forward in time and ff = -1 for backward in time.
    ff = 1  # will work for ff=1 or ff=-1 since checks by distance traveled

    ah = 0.  # m^2/s
    av = 0.  # m^2/s

    # turbulence/diffusion flag
    doturb = 0

    # two particles (starting positions)
    x0 = [22065., 22065.]
    y0 = [27777., 51587.]

    do3d = 0  # flag to set to 2-d

    z0 = 's'  # 'z' 'salt' 's'
    zpar = num_layers-1  # top layer

    usespherical = False

    # Get projection object
    proj = tracpy.tools.make_proj(setup='galveston', usebasemap=False)

    # Read in grid
    grid = tracpy.inout.readgrid(grid_filename, proj,
                                 vert_filename=currents_filename,
                                 usespherical=usespherical)

    # Initialize Tracpy class
    tp = Tracpy(currents_filename, grid, name=name, tseas=tseas,
                ndays=ndays, nsteps=nsteps, N=N, ff=ff, ah=ah, av=av,
                doturb=doturb, do3d=do3d, z0=z0, zpar=zpar,
                time_units=time_units, dtFromTracmass=dtFromTracmass,
                usespherical=usespherical)

    xp, yp, zp, t, T0, U, V = tracpy.run.run(tp, date, x0, y0)

    ## check the results:
    print xp.shape
    print xp
    print yp

    #eastward current, latitude should not change:
    assert np.allclose(y0, yp.T)

    # current velocity -- 0.1 m/s
    # position
    distance = (ndays * 24 * 3600 * 0.1)*ff

    print distance
    print xp[:, -1] - xp[:, 0]

    assert np.allclose(xp[:, -1] - xp[:, 0], distance)
