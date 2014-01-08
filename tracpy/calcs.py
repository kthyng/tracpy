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
    tstart = time.time()

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
    print 'time for finding ' + var + ': ', time.time()-tstart

    return varp


def get_dist(lon1, lons, lat1, lats): 
    '''
    Function to compute great circle distance between point lat1 and lon1 
    and arrays of points given by lons, lats or both same length arrays.
    Uses Haversine formula. Distance is in km.
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

    Example call:
    tracpy.calcs.rel_dispersion(dr.variables['lonp'][:], dr.variables['latp'][:], r=1, squared=True)
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


def rel_dispersion_comp(lonpc, latpc, tpc, lonp, latp, tp, r=1, squared=True):
    '''
    Calculate the relative dispersion of tracks lonp, latp as directly compared with
    the tracks described by lonpc, latpc. The two sets of tracks must start in the same 
    locations since this is assumed for making "pairs" of drifters for comparison (and 
    therefore pairs do not need to be found). The relative dispersion in this case is a
    measure of the difference between the two simulations, and is aimed at being used
    for examining differences in tracks due to changes in the numerical simulation.
    The tracks should also be coincident in time, but the script will find a way to match
    them up for the overlap periods.


    Inputs:
        lonpc, latpc    Longitude/latitude of the control drifter tracks [ndrifter,ntime]
        tpc             Time vector associated with lonpc, latpc
        lonp, latp      Longitude/latitude of the drifter tracks [ndrifter,ntime]
        tp              Time vector associated with lonp, latp
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

    Example call:
    dc = netCDF.Dataset('tracks/tseas_use300_nsteps1.nc') (5 min, control output)
    d = netCDF.Dataset('tracks/tseas_use1200_nsteps1.nc') (20 min, comparison output)
    tracpy.calcs.rel_dispersion_comp(dc.variables['lonp'][:], dc.variables['latp'][:], dc.variables['tp'][:],
                                     d.variables['lonp'][:], d.variables['latp'][:], d.variables['tp'][:],
                                     squared=True)
    '''


    # Account for time difference between two simulations
    # pdb.set_trace()
    tstart = time.time()

    # # Find overlapping period between the two sets of drifter tracks. For now, 
    # # assume that the "c" case is the one to match to start.
    # istart = find(tp==tpc[0]) # index to start secondary case at, in time
    # iend = find(tp==tpc[-1]) 
    # tp = tp[istart:iend+1]
    # # tp = tp[istart:]

    # Since drifters may be at different time steps, find the correct stride to
    # use to equalize, assuming one is a multiple of the other
    dt = (tp[1]-tp[0])
    dtc = (tpc[1]-tpc[0])

    if dt > dtc: # control case has higher temporal resolution
        tstride = int(dt/dtc)
    elif dtc > dt:
        tstride = int(dtc/dt)

    # Use tstride to equate arrays
    if dt > dtc:
        lonpc = lonpc[:,::tstride]
        latpc = latpc[:,::tstride]
        tpc = tpc[::tstride]
    elif dtc > dt:
        lonp = lonp[:,::tstride]
        latp = latp[:,::tstride]
        tp = tp[::tstride]

    print 'time for fixing timing: ', time.time()-tstart

    # Calculate relative dispersion
    tstart = time.time()
    # We know that drifters from the two sets have a one to one correspondence
    D2 = np.ones(lonp.shape[1])*np.nan
    nnans = np.zeros(lonp.shape[1]) # to collect number of non-nans over all drifters for a time
    for i in xrange(lonp.shape[0]): # loop through drifters, time is in array, axis=1
        dist = get_dist(lonpc[i,:], lonp[i,:], 
                    latpc[i,:], latp[i,:])
        if squared:
            D2 = np.nansum(np.vstack([D2, dist**2]), axis=0)
        else:
            D2 = np.nansum(np.vstack([D2, dist]), axis=0)
        nnans = nnans + ~np.isnan(dist)
    D2 = D2.squeeze()/nnans #len(pairs) # average over all pairs

    print 'time for finding numerical D: ', time.time()-tstart

    # # Distances squared, separately; times; number of non-nans for this set
    # np.savez(name[:-3] + 'D2.npz', D2=D2, t=t, nnans=nnans)
    # pdb.set_trace()
    return D2, nnans


def abs_dispersion(lonp, latp, squared=True):
    '''
    Calculate the absolute dispersion of a set of tracks. The distance between the position
    of a drifter and its original location are calculated in time, and averaged over all
    drifters.

    Inputs:
        lonp, latp      Longitude/latitude of the drifter tracks [ndrifter,ntime]
        squared         Whether to present the results as separation distance squared or 
                        not squared. Squared by default.

    Outputs:
        D2              Absolute dispersion (squared or not) averaged over drifter 
                        pairs [ntime].
        tp              Times along drifter track (input)
        nnans           Number of non-nan time steps in calculations for averaging properly.
                        Otherwise drifters that have exited the domain could affect calculations.

    To combine with other calculations of absolute dispersion, first multiply by nnans, then
    combine with other absolute dispersion calculations, then divide by the total number
    of nnans.

    Example call:
    tracpy.calcs.abs_dispersion(d.variables['lonp'][:], d.variables['latp'][:], squared=True)
    '''

    tstart = time.time()
    # Loop through pairs of drifters and calculate the separation distance in time
    D2 = np.ones(lonp.shape[1])*np.nan
    # to collect number of non-nans over all drifters and time steps
    nnans = np.zeros(lonp.shape[1]) 
    for idrifter in xrange(lonp.shape[0]):

        # calculate distance in time
        dist = get_dist(lonp[idrifter,0], lonp[idrifter,:], 
                    latp[idrifter,0], latp[idrifter,:])

        # dispersion can be presented as squared or not
        if squared:
            D2 = np.nansum(np.vstack([D2, dist**2]), axis=0)
        else:
            D2 = np.nansum(np.vstack([D2, dist]), axis=0)
        nnans = nnans + ~np.isnan(dist) # save these for averaging
    D2 = D2.squeeze()/nnans # average over all pairs
    print 'time for finding a: ', time.time()-tstart

    # # Distances squared, separately; times; number of non-nans for this set
    # np.savez(name[:-3] + 'D2.npz', D2=D2, t=t, nnans=nnans)
    # pdb.set_trace()
    return D2, nnans


def path(lonp, latp, squared=True):
    '''
    Calculate the path length in time of a set of tracks. The distance traveled or length
    of path from each drifter's initial location as it moves in time is calculated in time, 
    and the lengths are averaged over all drifters.

    Inputs:
        lonp, latp      Longitude/latitude of the drifter tracks [ndrifter,ntime]
        squared         Whether to present the results as separation distance squared or 
                        not squared. Squared by default.

    Outputs:
        D2              Path length (squared or not) averaged over drifter 
                        pairs [ntime].
        tp              Times along drifter track (input)
        nnans           Number of non-nan time steps in calculations for averaging properly.
                        Otherwise drifters that have exited the domain could affect calculations.

    To combine with other calculations of path length, first multiply by nnans, then
    combine with other path length calculations, then divide by the total number
    of nnans.

    Example call:
    tracpy.calcs.path(d.variables['lonp'][:], d.variables['latp'][:], squared=True)
    '''

    tstart = time.time()
    # Loop through pairs of drifters and calculate the separation distance in time
    D2 = np.ones(lonp.shape[1]-1)*np.nan
    # to collect number of non-nans over all drifters and time steps
    nnans = np.zeros(lonp.shape[1]-1) 
    for idrifter in xrange(lonp.shape[0]):

        # calculate distance in time
        dist = get_dist(lonp[idrifter,1:], lonp[idrifter,:-1], 
                    latp[idrifter,1:], latp[idrifter,:-1])
        dist = np.cumsum(dist)
        # pdb.set_trace()
        # dispersion can be presented as squared or not
        if squared:
            D2 = np.nansum(np.vstack([D2, dist**2]), axis=0)
        else:
            D2 = np.nansum(np.vstack([D2, dist]), axis=0)
        nnans = nnans + ~np.isnan(dist) # save these for averaging
    D2 = D2.squeeze()/nnans # average over all pairs
    print 'time for finding s: ', time.time()-tstart

    # # Distances squared, separately; times; number of non-nans for this set
    # np.savez(name[:-3] + 'D2.npz', D2=D2, t=t, nnans=nnans)
    # pdb.set_trace()
    return D2, nnans


def traj_ss(lon1, lat1, lon2, lat2):
    '''
    Trajectory skill score, from Liu and Weisberg, 2011
    '''

    # distance between drifters in time
    dist = get_dist(lon1, lon2, lat1, lat2) # in time

    # distance along path for control case, which is taken as lon1, lat1
    length = get_dist(lon1[:,:-1], lon1[:,1:], lat1[:,:-1], lat1[:,1:])

    # calculate s using cumulative sums
    # the first entry in time would be divided by zero, so this starts at the 2nd step
    s = np.cumsum(dist[:,1:], axis=1)/np.cumsum(length, axis=1)    

    # pdb.set_trace()
    # calculate skill score based on n=1
    ind = (s>1)
    ss = 1-s
    ss[ind] = 0.

    return ss


def moment1(xp):
    '''
    Calculate the 1st moment in a single direction of a set of tracks. in meters.

    Inputs:
        xp          x or y locations of the drifter tracks [ndrifter,ntime]

    Outputs:
        M               Relative dispersion (squared or not) averaged over drifter 
                        pairs [ntime].
        nnans           Number of non-nan time steps in calculations for averaging properly.
                        Otherwise drifters that have exited the domain could affect calculations.

    To combine with other calculations of relative dispersion, first multiply by nnans, then
    combine with other relative dispersion calculations, then divide by the total number
    of nnans.

    Example call:
    tracpy.calcs.moment1(xp)
    '''

    tstart = time.time()
    dists = (xp.T-xp[:,0]).T
    nnans = np.sum(~np.isnan(dists), axis=0)
    M = np.nansum(dists, axis=0)/nnans

    print 'time for finding M: ', time.time()-tstart

    # # Distances squared, separately; times; number of non-nans for this set
    # np.savez(name[:-3] + 'D2.npz', D2=D2, t=t, nnans=nnans)
    # pdb.set_trace()
    return M, nnans

def moment2(xp, M1):
    '''
    Calculate the 2nd moment in a single direction of a set of tracks. in meters.

    Inputs:
        xp          x or y locations of the drifter tracks [ndrifter,ntime]
        M1          First moment of tracks in a single direction

    Outputs:
        M2              Relative dispersion (squared or not) averaged over drifter 
                        pairs [ntime].
        nnans           Number of non-nan time steps in calculations for averaging properly.
                        Otherwise drifters that have exited the domain could affect calculations.

    To combine with other calculations of relative dispersion, first multiply by nnans, then
    combine with other relative dispersion calculations, then divide by the total number
    of nnans.

    Example call:
    tracpy.calcs.moment2(xp, M1x)
    '''

    tstart = time.time()
    dists = ((xp.T-xp[:,0]).T - M1)**2
    nnans = np.sum(~np.isnan(dists), axis=0) - 1
    M2 = np.nansum(dists, axis=0)/nnans

    print 'time for finding M: ', time.time()-tstart

    # # Distances squared, separately; times; number of non-nans for this set
    # np.savez(name[:-3] + 'D2.npz', D2=D2, t=t, nnans=nnans)
    # pdb.set_trace()
    return M2, nnans

def moment3(xp, M1):
    '''
    Calculate the 4th moment in a single direction of a set of tracks. in meters.

    Inputs:
        xp          x or y locations of the drifter tracks [ndrifter,ntime]
        M1          First moment of tracks in a single direction

    Outputs:
        M3              Relative dispersion (squared or not) averaged over drifter 
                        pairs [ntime].
        nnans           Number of non-nan time steps in calculations for averaging properly.
                        Otherwise drifters that have exited the domain could affect calculations.

    To combine with other calculations of relative dispersion, first multiply by nnans, then
    combine with other relative dispersion calculations, then divide by the total number
    of nnans.

    Example call:
    tracpy.calcs.moment3(xp, M1x)
    '''

    tstart = time.time()
    dist = xp - M1
    nnans = np.sum(~np.isnan(dist), axis=0)
    num = np.nansum(dist**3, axis=0)/nnans
    denom = (np.nansum(dist**2, axis=0)/nnans)**(3/2)
    M3 = num/denom

    print 'time for finding M: ', time.time()-tstart

    # # Distances squared, separately; times; number of non-nans for this set
    # np.savez(name[:-3] + 'D2.npz', D2=D2, t=t, nnans=nnans)
    # pdb.set_trace()
    return M3, nnans

def moment4(xp, M1):
    '''
    Calculate the 4th moment in a single direction of a set of tracks. in meters.

    Inputs:
        xp          x or y locations of the drifter tracks [ndrifter,ntime]
        M1          First moment of tracks in a single direction

    Outputs:
        M4              Relative dispersion (squared or not) averaged over drifter 
                        pairs [ntime].
        nnans           Number of non-nan time steps in calculations for averaging properly.
                        Otherwise drifters that have exited the domain could affect calculations.

    To combine with other calculations of relative dispersion, first multiply by nnans, then
    combine with other relative dispersion calculations, then divide by the total number
    of nnans.

    Example call:
    tracpy.calcs.moment4(xp, M1x)
    '''

    tstart = time.time()
    dist = xp - M1
    nnans = np.sum(~np.isnan(dist), axis=0)
    num = np.nansum(dist**4, axis=0)/nnans
    denom = (np.nansum(dist**2, axis=0)/nnans)**2
    # num = np.nansum(dist**4, axis=0)
    # denom = np.nansum(dist**2, axis=0)**2
    M4 = num/denom

    print 'time for finding M: ', time.time()-tstart

    # # Distances squared, separately; times; number of non-nans for this set
    # np.savez(name[:-3] + 'D2.npz', D2=D2, t=t, nnans=nnans)
    # pdb.set_trace()
    return M4, nnans