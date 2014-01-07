import numpy as np
import sys
import os
import op
import netCDF4 as netCDF
from mpl_toolkits.basemap import Basemap
import pdb
from matplotlib import delaunay
from matplotlib.pyplot import *
import glob
from datetime import datetime, timedelta
import time
from matplotlib.mlab import *
import inout
#import init
import plotting
import tools
from scipy import ndimage

def run(tp, date, lon0, lat0, T0=None, U=None, V=None):
# def run(loc, nsteps, ndays, ff, date, tseas, ah, av, lon0, lat0, z0, 
#         zpar, do3d, doturb, name, grid=None, dostream=0, N=1, 
#         T0=None, U=None, V=None, zparuv=None, tseas_use=None):
    '''
    FIX THIS FOR USING TRACPY CLASS
    To re-compile tracmass fortran code, type "make clean" and "make f2py", which will give 
    a file tracmass.so, which is the module we import above. Then in ipython, "run run.py"
    xend,yend,zend are particle locations at next step
    some variables are not specifically because f2py is hiding them from me:
     imt, jmt, km, ntractot
    Look at tracmass.step to see what it is doing and making optional at the end.
    Do this by importing tracmass and then tracmass.step?

    I am assuming here that the velocity field at two times are being input into tracmass
    such that the output is the position for the drifters at the time corresponding to the
    second velocity time. Each drifter may take some number of steps in between, but those
    are not saved.

    tp          TracPy object, from the Tracpy class.

    loc         Path to directory of grid and output files
    nsteps      Number of steps to do between model outputs (iter in tracmass) - sets the max 
                time step between drifter steps. Does not control the output sampling anymore.
    ndays       number of days to track the particles from start date
    ff          ff=1 to go forward in time and ff=-1 for backward in time
    date        Start date in datetime object
    tseas       Time between outputs in seconds
    ah          Horizontal diffusion in m^2/s. 
                See project values of 350, 100, 0, 2000. For -turb,-diffusion
    av          Vertical diffusion in m^2/s.
    do3d        for 3d flag, do3d=0 makes the run 2d and do3d=1 makes the run 3d
    doturb      turbulence/diffusion flag. 
                doturb=0 means no turb/diffusion,
                doturb=1 means adding parameterized turbulence
                doturb=2 means adding diffusion on a circle
                doturb=3 means adding diffusion on an ellipse (anisodiffusion)
    lon0        Drifter starting locations in x/zonal direction.
    lat0        Drifter starting locations in y/meridional direction.
    z0/zpar     For 3D drifter movement, turn off twodim flag in makefile.
                Then z0 should be an array of initial drifter depths. 
                The array should be the same size as lon0 and be negative
                for under water. Currently drifter depths need to be above 
                the seabed for every x,y particle location for the script to run.
                To do 3D but start at surface, use z0=zeros(ia.shape) and have
                 either zpar='fromMSL'
                choose fromMSL to have z0 starting depths be for that depth below the base 
                time-independent sea level (or mean sea level).
                choose 'fromZeta' to have z0 starting depths be for that depth below the
                time-dependent sea surface. Haven't quite finished the 'fromZeta' case.
                For 2D drifter movement, turn on twodim flag in makefile.
                Then: 
                set z0 to 's' for 2D along a terrain-following slice
                 and zpar to be the index of s level you want to use (0 to km-1)
                set z0 to 'rho' for 2D along a density surface
                 and zpar to be the density value you want to use
                 Can do the same thing with salinity ('salt') or temperature ('temp')
                 The model output doesn't currently have density though.
                set z0 to 'z' for 2D along a depth slice
                 and zpar to be the constant (negative) depth value you want to use
                To simulate drifters at the surface, set z0 to 's' 
                 and zpar = grid['km']-1 to put them in the upper s level
                 z0='s' is currently not working correctly!!!
                 In the meantime, do surface using the 3d set up option but with 2d flag set
    zparuv      (optional) Use this if the k index for the model output fields (e.g, u, v) is different
                 from the k index in the grid. This might happen if, for example, only the surface current
                 were saved, but the model run originally did have many layers. This parameter
                 represents the k index for the u and v output, not for the grid.
    tseas_use   (optional) Desired time between outputs in seconds, as opposed to the actual time between outputs
                 (tseas). Should be >= tseas since this is just an ability to use model output at less 
                 frequency than is available, probably just for testing purposes or matching other models.
                 Should to be a multiple of tseas (or will be rounded later).
    xp          x-locations in x,y coordinates for drifters
    yp          y-locations in x,y coordinates for drifters
    zp          z-locations (depths from mean sea level) for drifters
    t           time for drifter tracks
    name        Name of simulation to be used for netcdf file containing final tracks
    grid        (optional) Grid information, as read in by tracpy.inout.readgrid().
    N           Controls the output sampling. The length of time between model outputs is divided by N.
                Default is 1.

    The following inputs are for calculating Lagrangian stream functions
    dostream    Calculate streamfunctions (1) or not (0). Default is 0.
    U0, V0      (optional) Initial volume transports of drifters (m^3/s)
    U, V  (optional) Array aggregating volume transports as drifters move [imt-1,jmt], [imt,jmt-1]
    '''

    tic_start = time.time()
    tic_initial = time.time()

    # pdb.set_trace()
    # Convert date to number
    date = netCDF.date2num(date, tp.time_units)

    # Figure out what files will be used for this tracking
    nc, tinds = inout.setupROMSfiles(tp.currents_filename, date, tp.ff, tp.tout, tstride=tp.tstride)
    print tinds

    # Read in grid parameters into dictionary, grid, if haven't already
    if tp.grid is None:
        tp._readgrid()

    # Interpolate to get starting positions in grid space
    xstart0, ystart0, _ = tools.interpolate2d(lon0, lat0, tp.grid, 'd_ll2ij')
    # Do z a little lower down

    # Initialize seed locations 
    ia = np.ceil(xstart0) #[253]#,525]
    ja = np.ceil(ystart0) #[57]#,40]

    # don't use nan's
    # pdb.set_trace()
    ind2 = ~np.isnan(ia) * ~np.isnan(ja)
    ia = ia[ind2]
    ja = ja[ind2]
    xstart0 = xstart0[ind2]
    ystart0 = ystart0[ind2]

    dates = nc.variables['ocean_time'][:]   
    t0save = dates[tinds[0]] # time at start of drifter test from file in seconds since 1970-01-01, add this on at the end since it is big

    # Initialize drifter grid positions and indices
    xend = np.ones((ia.size,(len(tinds)-1)*tp.N+1))*np.nan
    yend = np.ones((ia.size,(len(tinds)-1)*tp.N+1))*np.nan
    zend = np.ones((ia.size,(len(tinds)-1)*tp.N+1))*np.nan
    zp = np.ones((ia.size,(len(tinds)-1)*tp.N+1))*np.nan
    # iend = np.ones((ia.size,(len(tinds)-1)*tp.N))*np.nan
    # jend = np.ones((ia.size,(len(tinds)-1)*tp.N))*np.nan
    # kend = np.ones((ia.size,(len(tinds)-1)*tp.N))*np.nan
    ttend = np.zeros((ia.size,(len(tinds)-1)*tp.N+1))
    t = np.zeros(((len(tinds)-1)*tp.N+1))
    flag = np.zeros((ia.size),dtype=np.int) # initialize all exit flags for in the domain

    # Initialize vertical stuff and fluxes
    # Read initial field in - to 'new' variable since will be moved
    # at the beginning of the time loop ahead
    if is_string_like(tp.z0): # isoslice case
        ufnew,vfnew,dztnew,zrtnew,zwtnew = inout.readfields(tinds[0],tp.grid,nc,tp.z0,tp.zpar,zparuv=tp.zparuv)
    else: # 3d case
        ufnew,vfnew,dztnew,zrtnew,zwtnew = inout.readfields(tinds[0],tp.grid,nc)
    # pdb.set_trace()
    ## Find zstart0 and ka
    # The k indices and z grid ratios should be on a wflux vertical grid,
    # which goes from 0 to km since the vertical velocities are defined
    # at the vertical cell edges. A drifter's grid cell is vertically bounded
    # above by the kth level and below by the (k-1)th level
    if is_string_like(tp.z0): # then doing a 2d isoslice
        # there is only one vertical grid cell, but with two vertically-
        # bounding edges, 0 and 1, so the initial ka value is 1 for all
        # isoslice drifters.
        ka = np.ones(ia.size) 

        # for s level isoslice, place drifters vertically at the center 
        # of the grid cell since that is where the u/v flux info is from.
        # For a rho/temp/density isoslice, we treat it the same way, such
        # that the u/v flux info taken at a specific rho/temp/density value
        # is treated as being at the center of the grid cells vertically.
        zstart0 = np.ones(ia.size)*0.5

    else:   # 3d case
        # Convert initial real space vertical locations to grid space
        # first find indices of grid cells vertically
        ka = np.ones(ia.size)*np.nan
        zstart0 = np.ones(ia.size)*np.nan

        if tp.zpar == 'fromMSL':
            for i in xrange(ia.size):
                # pdb.set_trace()
                ind = (tp.grid['zwt0'][ia[i],ja[i],:]<=tp.z0[i])
                # check to make sure there is at least one true value, so the z0 is shallower than the seabed
                if np.sum(ind): 
                    ka[i] = find(ind)[-1] # find value that is just shallower than starting vertical position
                # if the drifter starting vertical location is too deep for the x,y location, complain about it
                else:  # Maybe make this nan or something later
                    print 'drifter vertical starting location is too deep for its x,y location. Try again.'
                if (tp.z0[i] != tp.grid['zwt0'][ia[i],ja[i],ka[i]]) and (ka[i] != tp.grid['km']): # check this
                    ka[i] = ka[i]+1
                # Then find the vertical relative position in the grid cell by adding on the bit of grid cell
                zstart0[i] = ka[i] - abs(tp.z0[i]-tp.grid['zwt0'][ia[i],ja[i],ka[i]]) \
                                    /abs(tp.grid['zwt0'][ia[i],ja[i],ka[i]-1]-tp.grid['zwt0'][ia[i],ja[i],ka[i]])
        # elif zpar == 'fromZeta':
        #   for i in xrange(ia.size):
        #       pdb.set_trace()
        #       ind = (zwtnew[ia[i],ja[i],:]<=z0[i])
        #       ka[i] = find(ind)[-1] # find value that is just shallower than starting vertical position
        #       if (z0[i] != zwtnew[ia[i],ja[i],ka[i]]) and (ka[i] != grid['km']): # check this
        #           ka[i] = ka[i]+1
        #       # Then find the vertical relative position in the grid cell by adding on the bit of grid cell
        #       zstart0[i] = ka[i] - abs(z0[i]-zwtnew[ia[i],ja[i],ka[i]]) \
        #                           /abs(zwtnew[ia[i],ja[i],ka[i]-1]-zwtnew[ia[i],ja[i],ka[i]])

    # Find initial cell depths to concatenate to beginning of drifter tracks later
    zsave = tools.interpolate3d(xstart0, ystart0, zstart0, zwtnew)

    # Initialize x,y,z with initial seeded positions
    xend[:,0] = xstart0
    yend[:,0] = ystart0
    zend[:,0] = zstart0

    toc_initial = time.time()

    # j = 0 # index for number of saved steps for drifters
    tic_read = np.zeros(len(tinds))
    toc_read = np.zeros(len(tinds))
    tic_zinterp = np.zeros(len(tinds))
    toc_zinterp = np.zeros(len(tinds))
    tic_tracmass = np.zeros(len(tinds))
    toc_tracmass = np.zeros(len(tinds))
    # pdb.set_trace()
    xr3 = tp.grid['xr'].reshape((tp.grid['xr'].shape[0],tp.grid['xr'].shape[1],1)).repeat(zwtnew.shape[2],axis=2)
    yr3 = tp.grid['yr'].reshape((tp.grid['yr'].shape[0],tp.grid['yr'].shape[1],1)).repeat(zwtnew.shape[2],axis=2)

    # # Calculate subloop steps using input parameter dtFromTracmass. 
    # subloopsteps = 

    # Loop through model outputs. tinds is in proper order for moving forward
    # or backward in time, I think.
    for j,tind in enumerate(tinds[:-1]):

        print j

        # dtstep = 0.
        # while dtstep <= dtFromTracmass:

            # # interpolation constant. =1 if dtFromTracmass==tseas
            # r = dtFromTracmass/tseas

        #     ind = (flag[:] == 0) # indices where the drifters are still inside the domain
        #     xstart = xend[:,j*tp.N+i]
        #     ystart = yend[:,j*tp.N+i]
        #     zstart = zend[:,j*tp.N+i]



        #     dtstep = dtstep + dtFromTracmass



        # # # Loop over substeps between model outputs. This is for use with GNOME. Substeps will not necessarily divide
        # # # evenly into model output time, so there is a special statement for that. Also, if one doesn't need to 
        # # # access steps individually, such as in regular TracPy use, then this loop should collapse.
        # # for i in loopsteps:

        # # tic_read[j] = time.time()

        #     # mask out drifters that have exited the domain
        #     xstart = np.ma.masked_where(flag[:]==1,xstart)
        #     ystart = np.ma.masked_where(flag[:]==1,ystart)
        #     zstart = np.ma.masked_where(flag[:]==1,zstart)

        #     if not np.ma.compressed(xstart).any(): # exit if all of the drifters have exited the domain
        #         break

        #     # Do stepping in Tracpy class
        #     if tp.dostream:
        #         ufnew, vfnew, dztnew, zrtnew, zwtnew, xend[ind,j*tp.N+1:j*tp.N+tp.N+1],\
        #             yend[ind,j*tp.N+1:j*tp.N+tp.N+1],\
        #             zend[ind,j*tp.N+1:j*tp.N+tp.N+1],\
        #             zp[ind,j*tp.N+1:j*tp.N+tp.N+1],\
        #             flag[ind],\
        #             ttend[ind,j*tp.N+1:j*tp.N+tp.N+1], U, V = tp.step(tinds[j+1], nc, j, ttend[ind,j*tp.N], ufnew, vfnew, dztnew, zrtnew, zwtnew, 
        #                 xstart, ystart, zstart, T0[ind], U, V)
        #     else:
        #         ufnew, vfnew, dztnew, zrtnew, zwtnew, xend[ind,j*tp.N+1:j*tp.N+tp.N+1],\
        #             yend[ind,j*tp.N+1:j*tp.N+tp.N+1],\
        #             zend[ind,j*tp.N+1:j*tp.N+tp.N+1],\
        #             zp[ind,j*tp.N+1:j*tp.N+tp.N+1],\
        #             flag[ind],\
        #             ttend[ind,j*tp.N+1:j*tp.N+tp.N+1], U, V = tp.step(tinds[j+1], nc, j, ttend[ind,j*tp.N], ufnew, vfnew, dztnew, zrtnew, zwtnew, 
        #                 xstart, ystart, zstart)

        ind = (flag[:] == 0) # indices where the drifters are still inside the domain
        xstart = xend[:,j*tp.N]
        ystart = yend[:,j*tp.N]
        zstart = zend[:,j*tp.N]

        # mask out drifters that have exited the domain
        xstart = np.ma.masked_where(flag[:]==1,xstart)
        ystart = np.ma.masked_where(flag[:]==1,ystart)
        zstart = np.ma.masked_where(flag[:]==1,zstart)

        if not np.ma.compressed(xstart).any(): # exit if all of the drifters have exited the domain
            break

        # Do stepping in Tracpy class
        if tp.dostream:
            ufnew, vfnew, dztnew, zrtnew, zwtnew, xend[ind,j*tp.N+1:j*tp.N+tp.N+1],\
                yend[ind,j*tp.N+1:j*tp.N+tp.N+1],\
                zend[ind,j*tp.N+1:j*tp.N+tp.N+1],\
                zp[ind,j*tp.N+1:j*tp.N+tp.N+1],\
                flag[ind],\
                ttend[ind,j*tp.N+1:j*tp.N+tp.N+1], U, V = tp.step(tinds[j+1], nc, j, ttend[ind,j*tp.N], ufnew, vfnew, dztnew, zrtnew, zwtnew, 
                    xstart, ystart, zstart, T0[ind], U, V)
        else:
            ufnew, vfnew, dztnew, zrtnew, zwtnew, xend[ind,j*tp.N+1:j*tp.N+tp.N+1],\
                yend[ind,j*tp.N+1:j*tp.N+tp.N+1],\
                zend[ind,j*tp.N+1:j*tp.N+tp.N+1],\
                zp[ind,j*tp.N+1:j*tp.N+tp.N+1],\
                flag[ind],\
                ttend[ind,j*tp.N+1:j*tp.N+tp.N+1], U, V = tp.step(tinds[j+1], nc, j, ttend[ind,j*tp.N], ufnew, vfnew, dztnew, zrtnew, zwtnew, 
                    xstart, ystart, zstart)

    nc.close()
    # pdb.set_trace()
    ttend = ttend + t0save # add back in base time in seconds

    ## map coordinates interpolation
    # xp2, yp2, dt = tools.interpolate(xg,yg,grid,'m_ij2xy')
    # tic = time.time()
    lonp, latp, dt = tools.interpolate2d(xend,yend,tp.grid,'m_ij2ll',mode='constant',cval=np.nan)
    # print '2d interp time=', time.time()-tic

    # pdb.set_trace()

    runtime = time.time()-tic_start


    print "============================================="
    print ""
    print "Simulation name: ", tp.name
    print ""
    print "============================================="
    print "run time:\t\t\t", runtime
    print "---------------------------------------------"
    print "Time spent on:"

    initialtime = toc_initial-tic_initial
    print "\tInitial stuff: \t\t%4.2f (%4.2f%%)" % (initialtime, (initialtime/runtime)*100)

    readtime = np.sum(toc_read-tic_read)
    print "\tReading in fields: \t%4.2f (%4.2f%%)" % (readtime, (readtime/runtime)*100)

    zinterptime = np.sum(toc_zinterp-tic_zinterp)
    print "\tZ interpolation: \t%4.2f (%4.2f%%)" % (zinterptime, (zinterptime/runtime)*100)

    tractime = np.sum(toc_tracmass-tic_tracmass)
    print "\tTracmass: \t\t%4.2f (%4.2f%%)" % (tractime, (tractime/runtime)*100)
    print "============================================="

    # Save results to netcdf file
    if tp.dostream:
        inout.savetracks(lonp, latp, zp, ttend, tp.name, tp.nsteps, tp.N, tp.ff, tp.tseas_use, tp.ah, tp.av, \
                            tp.do3d, tp.doturb, tp.currents_filename, tp.T0, tp.U, tp.V)
        return lonp, latp, zp, ttend, tp.grid, T0, U, V
    else:
        inout.savetracks(lonp, latp, zp, ttend, tp.name, tp.nsteps, tp.N, tp.ff, tp.tseas_use, tp.ah, tp.av, \
                            tp.do3d, tp.doturb, tp.currents_filename)
        return lonp, latp, zp, ttend, tp.grid

# def start_run():
#     '''
#     Choose what initialization from above and then run.
#     '''

#     # Choose which initialization to use
#     loc,nsteps,ndays,ff,date,tseas,ah,av,lon0,lat0,z0,zpar,do3d,doturb,name = init.test1()

#     # Run tracmass!
#     lonp,latp,zp,t,grid = run(loc,nsteps,ndays,ff,date,tseas,ah,av,lon0,lat0,z0,zpar,do3d,doturb,name)

#     # pdb.set_trace()

#     # Plot tracks
#     plotting.tracks(lonp,latp,name,grid=grid)

#     # Plot final location (by time index) histogram
#     plotting.hist(lonp,latp,name,grid=grid,which='contour')
#     plotting.hist(lonp,latp,name,grid=grid,which='pcolor')  