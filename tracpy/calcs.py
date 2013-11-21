"""
Functions to find helpful calculations based on drifter tracks.

Functions include:

* 
"""

import numpy as np
import pdb
from matplotlib.mlab import find
import netCDF4 as netCDF
from scipy import ndimage


def Var(xp, yp, tp, varin, nc):
    '''
    Calculate the given property, varin, along the input drifter tracks. This property can
    be changing in time and space.

    NOTE: Currently assuming surface tracks.

    Inputs:
        xp, yp      Horizontal drifter position in grid coordinates [ndrift, ntime]
        tp          Times for drifter [ntime]
        varin       Variable to calculate. Available options are: u, v, salt, temp, h, zeta
        nc          Netcdf file object where the model output can be accessed which includes 
                    all necessary times

    Outputs:
        varp        Variable along the drifter track

    Example call: 
    first re-save drifter tracks in grid space using tracpy.inout.save_ll2grid()
    then open that drifter track file calling the object d
        d = netCDF.Dataset([trackfile in grid coords])
    have the ocean simulation output available in nc
        loc = 'http://barataria.tamu.edu:8080/thredds/dodsC/NcML/txla_nesting6.nc'
        nc = netCDF.Dataset(loc)
    Call this function 
        tracpy.calcs.Var(d.variables['xg'][:], d.variables['yg'][:], d.variables['tp'][:], 'h', nc)
    '''

    # Time indices for the drifter track points
    units = 'seconds since 1970-01-01'
    t = nc.variables['ocean_time'][:] # model times
    istart = find(netCDF.num2date(t, units) <= netCDF.num2date(tp[0], units))[-1]
    iend = find(netCDF.num2date(t, units) >= netCDF.num2date(tp[-1], units))[0]
    tinds = np.arange(istart, iend)


    # Read in model information. Try reading it all in the for time, y, and x and then
    # interpolating from there.

    # 4D variables
    if varin in ('u', 'v', 'salt', 'temp'):
        var = nc.variables[varin][tinds, -1, :, :]

    # 3D variables
    elif varin in ('zeta'):
        var = nc.variables[varin][tinds, :, :]

    # 2D variables
    elif varin in ('h'):
        var = nc.variables[varin][:, :]


    # Grid location of var. xp and yp are on staggered grids, counting from the cell
    # edge and inward. Want to match the grid locations of these and var.
    # varin on ugrid
    if varin in ('u'):
        # xp correctly align with the u grid in the x direction, so that xp is at the cell center
        # this shifts the drifter y grid locations from the edge of the cell to the center
        yp = yp - 0.5 

    # varin on vgrid
    elif varin in ('v'):
        # this shifts the drifter x grid locations from the edge of the cell to the center
        xp = xp - 0.5 
        # yp correctly align with the v grid in the y direction, so that yp is at the cell center

    # varin on rho grid
    elif varin in ('salt', 'temp', 'zeta', 'h'):
        # shift both xp and yp to cell center from edge
        xp = xp - 0.5; yp = yp - 0.5;

    # Interpolate var to tracks. varp is the variable along the tracks
    # Need to know grid location of everything
    # h does not change in time so need to interpolate differently
    if varin == 'h':
        varp = ndimage.map_coordinates(var, np.array([yp.flatten(), \
                                    xp.flatten()]), \
                                    order=1, \
                                    mode='nearest').reshape(xp.shape)
    # these variables change in time
    else:
        tg = ((tp-tp[0])/(tp[-1]-tp[0]))*var.shape[0] # Make time into a "grid coordinate" that goes from 0 to number of time indices
        varp = ndimage.map_coordinates(var, np.array([tg.reshape((1,tg.shape[0])).repeat(xp.shape[0],axis=0).flatten(), \
                                    yp.flatten(), \
                                    xp.flatten()]), \
                                    order=1, \
                                    mode='nearest').reshape(xp.shape)

    return varp


