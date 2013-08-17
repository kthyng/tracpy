import numpy as np
import sys
import os
import op
import tracmass
import netCDF4 as netCDF
from mpl_toolkits.basemap import Basemap
import pdb
from matplotlib import delaunay
from matplotlib.pyplot import *
import glob
from datetime import datetime, timedelta
from mpl_toolkits.basemap import Basemap
import time
from matplotlib.mlab import *
import inout
import init
import plotting
import tools
from scipy import ndimage

def run(loc, nsteps, ndays, ff, date, tseas, ah, av, lon0, lat0, z0, \
        zpar, do3d, doturb, name, grid=None, dostream=0, \
        T0=None, U=None, V=None):
    '''

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

    loc         Path to directory of grid and output files
    nsteps      Number of steps to do between model outputs (iter in tracmass)
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
    xp          x-locations in x,y coordinates for drifters
    yp          y-locations in x,y coordinates for drifters
    zp          z-locations (depths from mean sea level) for drifters
    t           time for drifter tracks
    name        Name of simulation to be used for netcdf file containing final tracks
    grid        (optional) Grid information, as read in by tracpy.inout.readgrid().

    The following inputs are for calculating Lagrangian stream functions
    dostream    Calculate streamfunctions (1) or not (0). Default is 0.
    U0, V0      (optional) Initial volume transports of drifters (m^3/s)
    U, V  (optional) Array aggregating volume transports as drifters move [imt-1,jmt], [imt,jmt-1]
    '''

    tic_start = time.time()
    tic_initial = time.time()

    # Units for time conversion with netCDF.num2date and .date2num
    units = 'seconds since 1970-01-01'

    # Number of model outputs to use
    # Adding one index so that all necessary indices are captured by this number.
    # Then the run loop uses only the indices determined by tout instead of needing
    # an extra one beyond
    tout = np.int((ndays*(24*3600))/tseas + 1)

    # Convert date to number
    date = netCDF.date2num(date, units)

    # Figure out what files will be used for this tracking
    nc, tinds = inout.setupROMSfiles(loc, date, ff, tout)

    # Read in grid parameters into dictionary, grid
    if grid is None:
        grid = inout.readgrid(loc, nc)
    else: # don't need to reread grid
        grid = grid

    # Interpolate to get starting positions in grid space
    xstart0, ystart0, _ = tools.interpolate2d(lon0, lat0, grid, 'd_ll2ij')
    # Do z a little lower down

    # Initialize seed locations 
    ia = np.ceil(xstart0) #[253]#,525]
    ja = np.ceil(ystart0) #[57]#,40]

    # don't use nan's
    ind2 = ~np.isnan(ia) * ~np.isnan(ja)
    ia = ia[ind2]
    ja = ja[ind2]
    xstart0 = xstart0[ind2]
    ystart0 = ystart0[ind2]

    dates = nc.variables['ocean_time'][:]   
    t0save = dates[tinds[0]] # time at start of drifter test from file in seconds since 1970-01-01, add this on at the end since it is big

    # Initialize drifter grid positions and indices
    xend = np.ones((ia.size,(len(tinds)-1)*nsteps))*np.nan
    yend = np.ones((ia.size,(len(tinds)-1)*nsteps))*np.nan
    zend = np.ones((ia.size,(len(tinds)-1)*nsteps))*np.nan
    zp = np.ones((ia.size,(len(tinds)-1)*nsteps))*np.nan
    iend = np.ones((ia.size,(len(tinds)-1)*nsteps))*np.nan
    jend = np.ones((ia.size,(len(tinds)-1)*nsteps))*np.nan
    kend = np.ones((ia.size,(len(tinds)-1)*nsteps))*np.nan
    ttend = np.ones((ia.size,(len(tinds)-1)*nsteps))*np.nan
    t = np.zeros(((len(tinds)-1)*nsteps+1))
    flag = np.zeros((ia.size),dtype=np.int) # initialize all exit flags for in the domain

    # Initialize vertical stuff and fluxes
    # Read initial field in - to 'new' variable since will be moved
    # at the beginning of the time loop ahead
    if is_string_like(z0): # isoslice case
        ufnew,vfnew,dztnew,zrtnew,zwtnew = inout.readfields(tinds[0],grid,nc,z0,zpar)
    else: # 3d case
        ufnew,vfnew,dztnew,zrtnew,zwtnew = inout.readfields(tinds[0],grid,nc)

    ## Find zstart0 and ka
    # The k indices and z grid ratios should be on a wflux vertical grid,
    # which goes from 0 to km since the vertical velocities are defined
    # at the vertical cell edges. A drifter's grid cell is vertically bounded
    # above by the kth level and below by the (k-1)th level
    if is_string_like(z0): # then doing a 2d isoslice
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

        if zpar == 'fromMSL':
            for i in xrange(ia.size):
                # pdb.set_trace()
                ind = (grid['zwt0'][ia[i],ja[i],:]<=z0[i])
                # check to make sure there is at least one true value, so the z0 is shallower than the seabed
                if np.sum(ind): 
                    ka[i] = find(ind)[-1] # find value that is just shallower than starting vertical position
                # if the drifter starting vertical location is too deep for the x,y location, complain about it
                else:  # Maybe make this nan or something later
                    print 'drifter vertical starting location is too deep for its x,y location. Try again.'
                if (z0[i] != grid['zwt0'][ia[i],ja[i],ka[i]]) and (ka[i] != grid['km']): # check this
                    ka[i] = ka[i]+1
                # Then find the vertical relative position in the grid cell by adding on the bit of grid cell
                zstart0[i] = ka[i] - abs(z0[i]-grid['zwt0'][ia[i],ja[i],ka[i]]) \
                                    /abs(grid['zwt0'][ia[i],ja[i],ka[i]-1]-grid['zwt0'][ia[i],ja[i],ka[i]])
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

    toc_initial = time.time()

    # j = 0 # index for number of saved steps for drifters
    tic_read = np.zeros(len(tinds))
    toc_read = np.zeros(len(tinds))
    tic_zinterp = np.zeros(len(tinds))
    toc_zinterp = np.zeros(len(tinds))
    tic_tracmass = np.zeros(len(tinds))
    toc_tracmass = np.zeros(len(tinds))
    # pdb.set_trace()
    xr3 = grid['xr'].reshape((grid['xr'].shape[0],grid['xr'].shape[1],1)).repeat(zwtnew.shape[2],axis=2)
    yr3 = grid['yr'].reshape((grid['yr'].shape[0],grid['yr'].shape[1],1)).repeat(zwtnew.shape[2],axis=2)
    # Loop through model outputs. tinds is in proper order for moving forward
    # or backward in time, I think.
    nmodel = 30 # number of extra model outputs to read in at a time
    for jj in xrange(0,len(tinds[:-1]),nmodel):
    # for j,tind in enumerate(tinds[:-1]):
        # pdb.set_trace()
        # # Move previous new time step to old time step info
        # ufold = ufnew
        # vfold = vfnew
        # dztold = dztnew
        # zrtold = zrtnew
        # zwtold = zwtnew

        tic_read[jj] = time.time()
        # pdb.set_trace()
        # Read stuff in for next time loop
        if is_string_like(z0): # isoslice case
            ufnew,vfnew,dztnew,zrtnew,zwtnew = inout.readfields(tinds[jj:jj+nmodel+1],grid,nc,z0,zpar)
        else: # 3d case
            ufnew,vfnew,dztnew,zrtnew,zwtnew = inout.readfields(tinds[jj:jj+nmodel+1],grid,nc)
        toc_read[jj] = time.time()
        # print "readfields run time:",toc_read-tic_read

        for j in xrange(nmodel): # loop between the model outputs

            print j+jj
            # pdb.set_trace()
            #  flux fields at starting time for this step
            if (j+jj) != 0:
                xstart = xend[:,(j+jj)*nsteps-1]
                ystart = yend[:,(j+jj)*nsteps-1]
                zstart = zend[:,(j+jj)*nsteps-1]
                # mask out drifters that have exited the domain
                xstart = np.ma.masked_where(flag[:]==1,xstart)
                ystart = np.ma.masked_where(flag[:]==1,ystart)
                zstart = np.ma.masked_where(flag[:]==1,zstart)
                ind = (flag[:] == 0) # indices where the drifters are still inside the domain
            else: # first loop, j==0
                xstart = xstart0
                ystart = ystart0
                zstart = zstart0
                # TODO: Do a check to make sure all drifter starting locations are within domain
                ind = (flag[:] == 0) # indices where the drifters are inside the domain to start

            # Find drifter locations
            # only send unmasked values to step
            if not np.ma.compressed(xstart).any(): # exit if all of the drifters have exited the domain
                break
            else:
                # Combine times for arrays for input to tracmass
                # from [ixjxk] to [ixjxkxt]
                # Change ordering for these three arrays here instead of in readfields since
                # concatenate does not seem to preserve ordering
                # uflux = np.asfortranarray(np.concatenate((ufold.reshape(np.append(ufold.shape,1)), \
                #                         ufnew.reshape(ufnew.shape)),axis=ufold.ndim))
                # vflux = np.asfortranarray(np.concatenate((vfold.reshape(np.append(vfold.shape,1)), \
                #                         vfnew.reshape(vfnew.shape)),axis=vfold.ndim))
                # dzt = np.asfortranarray(np.concatenate((dztold.reshape(np.append(dztold.shape,1)), \
                #                         dztnew.reshape(dztnew.shape)),axis=dztold.ndim))
                # # uflux = np.asfortranarray(np.concatenate((ufold.reshape(np.append(ufold.shape,1)), \
                # #                         ufnew.reshape(np.append(ufnew.shape,1))), \
                # #                         axis=ufold.ndim))
                # # vflux = np.asfortranarray(np.concatenate((vfold.reshape(np.append(vfold.shape,1)), \
                # #                         vfnew.reshape(np.append(vfnew.shape,1))), \
                # #                         axis=vfold.ndim))
                # # dzt = np.asfortranarray(np.concatenate((dztold.reshape(np.append(dztold.shape,1)), \
                # #                         dztnew.reshape(np.append(dztnew.shape,1))), \
                # #                         axis=dztold.ndim))

                # Change the horizontal indices from python to fortran indexing 
                # (vertical are zero-based in tracmass)
                xstart, ystart = tools.convert_indices('py2f',xstart,ystart)

                # km that is sent to tracmass is determined from uflux (see tracmass?)
                # so it will be the correct value for whether we are doing the 3D
                # or isoslice case.
                # vec = np.arange(j*nsteps,j*nsteps+nsteps) # indices for storing new track locations
                tic_tracmass[j+jj] = time.time()
                # pdb.set_trace()
                if dostream: # calculate Lagrangian stream functions
                    xend[ind,(j+jj)*nsteps:(j+jj)*nsteps+nsteps],\
                        yend[ind,(j+jj)*nsteps:(j+jj)*nsteps+nsteps],\
                        zend[ind,(j+jj)*nsteps:(j+jj)*nsteps+nsteps], \
                        iend[ind,(j+jj)*nsteps:(j+jj)*nsteps+nsteps],\
                        jend[ind,(j+jj)*nsteps:(j+jj)*nsteps+nsteps],\
                        kend[ind,(j+jj)*nsteps:(j+jj)*nsteps+nsteps],\
                        flag[ind],\
                        ttend[ind,(j+jj)*nsteps:(j+jj)*nsteps+nsteps], U, V = \
                            tracmass.step(np.ma.compressed(xstart),\
                                            np.ma.compressed(ystart),
                                            np.ma.compressed(zstart),
                                            tseas, ufnew[:,:,:,j:j+2], vfnew[:,:,:,j:j+2], \
                                            ff, grid['kmt'].astype(int), \
                                            dztnew[:,:,:,j:j+2], grid['dxdy'], grid['dxv'], \
                                            grid['dyu'], grid['h'], nsteps, \
                                            ah, av, do3d, doturb, dostream, \
                                            t0=T0[ind],
                                            ut=U, vt=V)
                    # xend[ind,j*nsteps:j*nsteps+nsteps],\
                    #     yend[ind,j*nsteps:j*nsteps+nsteps],\
                    #     zend[ind,j*nsteps:j*nsteps+nsteps], \
                    #     iend[ind,j*nsteps:j*nsteps+nsteps],\
                    #     jend[ind,j*nsteps:j*nsteps+nsteps],\
                    #     kend[ind,j*nsteps:j*nsteps+nsteps],\
                    #     flag[ind],\
                    #     ttend[ind,j*nsteps:j*nsteps+nsteps], U, V = \
                    #         tracmass.step(np.ma.compressed(xstart),\
                    #                         np.ma.compressed(ystart),
                    #                         np.ma.compressed(zstart),
                    #                         tseas, uflux, vflux, ff, \
                    #                         grid['kmt'].astype(int), \
                    #                         dzt, grid['dxdy'], grid['dxv'], \
                    #                         grid['dyu'], grid['h'], nsteps, \
                    #                         ah, av, do3d, doturb, dostream, \
                    #                         t0=T0[ind],
                    #                         ut=U, vt=V)
                else: # don't calculate Lagrangian stream functions
                    xend[ind,(j+jj)*nsteps:(j+jj)*nsteps+nsteps],\
                        yend[ind,(j+jj)*nsteps:(j+jj)*nsteps+nsteps],\
                        zend[ind,(j+jj)*nsteps:(j+jj)*nsteps+nsteps], \
                        iend[ind,(j+jj)*nsteps:(j+jj)*nsteps+nsteps],\
                        jend[ind,(j+jj)*nsteps:(j+jj)*nsteps+nsteps],\
                        kend[ind,(j+jj)*nsteps:(j+jj)*nsteps+nsteps],\
                        flag[ind],\
                        ttend[ind,(j+jj)*nsteps:(j+jj)*nsteps+nsteps], _, _ = \
                            tracmass.step(np.ma.compressed(xstart),\
                                            np.ma.compressed(ystart),
                                            np.ma.compressed(zstart),
                                            tseas, ufnew[:,:,:,j:j+2], vfnew[:,:,:,j:j+2], \
                                            ff, grid['kmt'].astype(int), \
                                            dztnew[:,:,:,j:j+2], grid['dxdy'], grid['dxv'], \
                                            grid['dyu'], grid['h'], nsteps, \
                                            ah, av, do3d, doturb, dostream)
                    # xend[ind,j*nsteps:j*nsteps+nsteps],\
                    #     yend[ind,j*nsteps:j*nsteps+nsteps],\
                    #     zend[ind,j*nsteps:j*nsteps+nsteps], \
                    #     iend[ind,j*nsteps:j*nsteps+nsteps],\
                    #     jend[ind,j*nsteps:j*nsteps+nsteps],\
                    #     kend[ind,j*nsteps:j*nsteps+nsteps],\
                    #     flag[ind],\
                    #     ttend[ind,j*nsteps:j*nsteps+nsteps], _, _ = \
                    #         tracmass.step(np.ma.compressed(xstart),\
                    #                         np.ma.compressed(ystart),
                    #                         np.ma.compressed(zstart),
                    #                         tseas, uflux, vflux, ff, \
                    #                         grid['kmt'].astype(int), \
                    #                         dzt, grid['dxdy'], grid['dxv'], \
                    #                         grid['dyu'], grid['h'], nsteps, \
                    #                         ah, av, do3d, doturb, dostream)
                toc_tracmass[j+jj] = time.time()
                # pdb.set_trace()

                # Change the horizontal indices from python to fortran indexing
                xend[ind,(j+jj)*nsteps:(j+jj)*nsteps+nsteps], \
                    yend[ind,(j+jj)*nsteps:(j+jj)*nsteps+nsteps] \
                                    = tools.convert_indices('f2py', \
                                        xend[ind,(j+jj)*nsteps:(j+jj)*nsteps+nsteps], \
                                        yend[ind,(j+jj)*nsteps:(j+jj)*nsteps+nsteps])

                # Calculate times for the output frequency
                if ff == 1:
                    t[(j+jj)*nsteps+1:(j+jj)*nsteps+nsteps+1] = t[(j+jj)*nsteps] + np.linspace(tseas/nsteps,tseas,nsteps) # update time in seconds to match drifters
                else:
                    t[(j+jj)*nsteps+1:(j+jj)*nsteps+nsteps+1] = t[(j+jj)*nsteps] - np.linspace(tseas/nsteps,tseas,nsteps) # update time in seconds to match drifters
                
                # Skip calculating real z position if we are doing surface-only drifters anyway
                if z0 != 's' and zpar != grid['km']-1:
                    tic_zinterp[j+jj] = time.time()
                    # Calculate real z position
                    r = np.linspace(1./nsteps,1,nsteps) # linear time interpolation constant that is used in tracmass

                    for n in xrange(nsteps): # loop through time steps
                        # interpolate to a specific output time
                        # pdb.set_trace()
                        zwt = (1.-r[n])*zwtold + r[n]*zwtnew
                        zp[ind,(j+jj)*nsteps:(j+jj)*nsteps+nsteps], dt = tools.interpolate3d(xend[ind,(j+jj)*nsteps:(j+jj)*nsteps+nsteps], \
                                                                yend[ind,(j+jj)*nsteps:(j+jj)*nsteps+nsteps], \
                                                                zend[ind,(j+jj)*nsteps:(j+jj)*nsteps+nsteps], \
                                                                zwt)
                    toc_zinterp[j+jj] = time.time()

    nc.close()
    t = t + t0save # add back in base time in seconds

    # pdb.set_trace()

    # Add on to front location for first time step
    xg=np.concatenate((xstart0.reshape(xstart0.size,1),xend),axis=1)
    yg=np.concatenate((ystart0.reshape(ystart0.size,1),yend),axis=1)
    # Concatenate zp with initial real space positions
    zp=np.concatenate((zsave[0].reshape(zstart0.size,1),zp),axis=1)

    # Delaunay interpolation
    # xp, yp, dt = tools.interpolate(xg,yg,grid,'d_ij2xy')
    # lonp, latp, dt = tools.interpolation(xg,yg,grid,'d_ij2ll')

    ## map coordinates interpolation
    # xp2, yp2, dt = tools.interpolate(xg,yg,grid,'m_ij2xy')
    # tic = time.time()
    lonp, latp, dt = tools.interpolate2d(xg,yg,grid,'m_ij2ll',mode='constant',cval=np.nan)
    # print '2d interp time=', time.time()-tic

    # pdb.set_trace()

    runtime = time.time()-tic_start

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

    # Save results to netcdf file
    if dostream:
        inout.savetracks(lonp, latp, zp, t, name, nsteps, ff, tseas, ah, av, \
                            do3d, doturb, loc, T0, U, V)
        return lonp, latp, zp, t, grid, T0, U, V
    else:
        inout.savetracks(lonp, latp, zp, t, name, nsteps, ff, tseas, ah, av, \
                            do3d, doturb, loc)
        return lonp, latp, zp, t, grid

def start_run():
    '''
    Choose what initialization from above and then run.
    '''

    # Choose which initialization to use
    loc,nsteps,ndays,ff,date,tseas,ah,av,lon0,lat0,z0,zpar,do3d,doturb,name = init.test1()

    # Run tracmass!
    lonp,latp,zp,t,grid = run(loc,nsteps,ndays,ff,date,tseas,ah,av,lon0,lat0,z0,zpar,do3d,doturb,name)

    # pdb.set_trace()

    # Plot tracks
    plotting.tracks(lonp,latp,name,grid=grid)

    # Plot final location (by time index) histogram
    plotting.hist(lonp,latp,name,grid=grid,which='contour')
    plotting.hist(lonp,latp,name,grid=grid,which='pcolor')  