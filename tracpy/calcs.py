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
import time

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
        varp = tracpy.calcs.Var(d.variables['xg'][:], d.variables['yg'][:], d.variables['tp'][:], 'h', nc)
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


def get_dist(lon1, lons, lat1, lats): 
    '''
    Function to compute great circle distance between point lat1 and lon1 
    and arrays of points given by lons, lats or both same length arrays.
    Uses Haversine formula.
    '''

    lon1 = lon1*np.pi/180.
    lons = lons*np.pi/180.
    lat1 = lat1*np.pi/180.
    lats = lats*np.pi/180.

    earth_radius = 6373.
    distance = earth_radius*2.0*np.arcsin(np.sqrt(np.sin(0.50*(lat1-lats))**2 \
                                       + np.cos(lat1)*np.cos(lats) \
                                       * np.sin(0.50*(lon1-lons))**2))
    return distance


def rel_dispersion(lonp, latp, r=1, squared=True):
    '''
    Calculate the relative dispersion of a set of tracks. First, initial pairs
    of drifters are found, based on a maximum initial separation distance, then
    find the separation distance in time.

    Inputs:
        lonp, latp      Longitude/latitude of the drifter tracks [ndrifter,ntime]
        r               Initial separation distance (kilometers) defining pairs of drifters. 
                        Default is 1 kilometer.
        squared         Whether to present the results as separation distance squared or 
                        not squared. Squared by default.

    Outputs:
        D2              Relative dispersion (squared or not) averaged over drifter 
                        pairs [ntime].
        tp              Times along drifter track (input)
        nnans           Number of non-nan time steps in calculations for averaging properly.
                        Otherwise drifters that have exited the domain could affect calculations.

    To combine with other calculations of relative dispersion, first multiply by nnans, then
    combine with other relative dispersion calculations, then divide by the total number
    of nnans.
    '''

    # Find pairs of drifters based on initial position

    # Calculate the initial separation distances for each drifter from each other drifter
    tstart = time.time()
    dist = np.zeros((lonp.shape[0],lonp.shape[0]))*np.nan
    for idrifter in xrange(lonp.shape[0]):
        # dist contains all of the distances from other drifters for each drifter
        dist[idrifter, idrifter+1:] = get_dist(lonp[idrifter,0], lonp[idrifter+1:,0], 
                                    latp[idrifter,0], latp[idrifter+1:,0])
        # dist[idrifter,:] = get_dist(lonp[idrifter,0], lonp[:,0], latp[idrifter,0], latp[:,0])
    print 'time for initial particle separation: ', time.time()-tstart

    tstart = time.time()
    # let the index in axis 0 be the drifter id
    ID = np.arange(lonp.shape[0])

    # # save pairs to save time since they are always the same
    # if not os.path.exists('tracks/pairs.npz'):

    # Loop through all drifters and find initial separation distances smaller than r.
    # Then exclude repeated pairs.
    pairs = []
    for idrifter in xrange(lonp.shape[0]):
        ind = find(dist[idrifter,:]<=r)
        for i in ind:
            if ID[idrifter] != ID[i]:
                pairs.append([min(ID[idrifter], ID[i]), 
                                max(ID[idrifter], ID[i])])
    pairs_set = set(map(tuple,pairs))
    pairs = map(list,pairs_set)# now pairs has only unique pairs of drifters
    # pairs.sort() #unnecessary but handy for checking work
    #     np.savez('tracks/pairs.npz', pairs=pairs)
    # else:
    #     pairs = np.load('tracks/pairs.npz')['pairs']
    print 'time for finding pairs: ', time.time()-tstart
    # Calculate relative dispersion

    tstart = time.time()
    # Loop through pairs of drifters and calculate the separation distance in time
    D2 = np.ones(lonp.shape[1])*np.nan
    # to collect number of non-nans over all drifters and time steps
    nnans = np.zeros(lonp.shape[1]) 
    for ipair in xrange(len(pairs)):

        # calculate distance in time
        dist = get_dist(lonp[pairs[ipair][0],:], lonp[pairs[ipair][1],:], 
                    latp[pairs[ipair][0],:], latp[pairs[ipair][1],:])

        # dispersion can be presented as squared or not
        if squared:
            D2 = np.nansum(np.vstack([D2, dist**2]), axis=0)
        else:
            D2 = np.nansum(np.vstack([D2, dist]), axis=0)
        nnans = nnans + ~np.isnan(dist) # save these for averaging
    D2 = D2.squeeze()/nnans # average over all pairs
    print 'time for finding D: ', time.time()-tstart

    # # Distances squared, separately; times; number of non-nans for this set
    # np.savez(name[:-3] + 'D2.npz', D2=D2, t=t, nnans=nnans)
    # pdb.set_trace()
    return D2, nnans
