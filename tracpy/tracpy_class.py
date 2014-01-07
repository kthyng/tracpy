#!/usr/bin/env python

'''
TracPy class
HAVE THINGS THAT ARE THE SAME FOR EACH STEP IN GNOME BE STORED IN 
THIS CLASS
* flush this out
* (keep adding tests)
* see where to put this into tracpy_mover.py
'''

import tracpy
import numpy as np
from matplotlib.pyplot import is_string_like
import pdb
import tracmass

class Tracpy(object):
    '''
    TracPy class.
    '''

    def __init__(self, currents_filename, grid_filename=None, nsteps=1, ndays=1, ff=1, tseas=3600.,
                ah=0., av=0., z0='s', zpar=1, do3d=0, doturb=0, name='test', dostream=0, N=1, 
                time_units='seconds since 1970-01-01', dtFromTracmass=None, zparuv=None, tseas_use=None):
        '''
        Initialize class.

        :param currents_filename: NetCDF file name (with extension) or OpenDAP url.
        :param grid_filename=None: NetCDF grid file name or OpenDAP url.
        :param nsteps=1: number of linearly interpolated steps between model outputs.
        :param ndays=1: number of run days
        :param ff=1: 1 is forward in time, -1 is backward
        :param tseas=3600.: number of seconds between model outputs
        :param ah=0.: horizontal diffusivity, in m^2/s
        :param av=0.: vertical diffusivity, in m^2/s
        :param z0='s': string flag in 2D case or array of initial z locations in 3D case
        :param zpar=1: isoslice value to in 2D case or string flag in 3D case
        :param do3d=0: 1 for 3D or 0 for 2D
        :param doturb=0: 0 for no added diffusion, 1 for diffusion vs velocity fluctuation, 2/3 for diffusion via random walk (3 for aligned with isobaths)
        :param name='test': name for output
        :param dostream=0: 1 to calculate transport for lagrangian stream functions, 0 to not
        :param N=None: number of steps between model outputs for outputting drifter locations. Defaults to output at nsteps. 
            If dtFromTracmass is being used, N is set by that.
        :param time_units='seconds since 1970-01-01': Reference for time, for changing between numerical times and datetime format
        :param dtFromTracmass=None: Time period for exiting from TRACMASS. If uninitialized, this is set to tseas 
            so that it only exits TRACMASS when it has gone through a full model output. If initialized by the user, 
            TRACMASS will run for 1 time step of length dtFromTracmass before exiting to the loop.
        :param zparuv=None: Defaults to zpar. Use this if the k index for the model output fields (e.g, u, v) is different from the k index in the grid
        :param tseas_use=None: Defaults to tseas. Desired time between outputs in seconds, as opposed to the actual time between outputs (tseas)
        '''

        self.currents_filename = currents_filename
        self.grid_filename = grid_filename
        self.grid = None

        # Initial parameters
        self.nsteps = nsteps
        self.ndays = ndays
        self.ff = ff
        self.tseas = tseas
        self.ah = ah
        self.av = av
        self.z0 = z0
        self.zpar = zpar
        self.do3d = do3d
        self.doturb = doturb
        self.name = name
        self.dostream = dostream
        self.N = N
        self.time_units = time_units

        # if loopsteps is None and nsteps is not None:
        #     # Use nsteps in TRACMASS and have inner loop collapse
        #     self.loopsteps = 1
        # elif loopsteps is not None and nsteps is None:
        #     # This means to use the inner loop (with loopsteps) and nsteps=1 to just do 1 step per call to TRACMASS
        #     self.nsteps = 1
        # elif loopsteps is None and nsteps is None:
        #     print 'need to input a value for nsteps or loopsteps.'
        #     break

        if dtFromTracmass is None:
            self.dtFromTracmass = tseas
        else:
            # If using dtFromTracmass, N=1, for steps between tracmass exits
            self.N = 1
            # # If using dtFromTracmass, N is set according to that.
            # self.N = (self.ndays*3600*24.)/self.tseas # this is the total number of steps
        if zparuv is None:
            self.zparuv = zpar
        if tseas_use is None:
            self.tseas_use = tseas

        # Calculate parameters that derive from other parameters

        # Number of model outputs to use (based on tseas, actual amount of model output)
        # This should not be updated with tstride since it represents the full amount of
        # indices in the original model output. tstride will be used separately to account
        # for the difference.
        # Adding one index so that all necessary indices are captured by this number.
        # Then the run loop uses only the indices determined by tout instead of needing
        # an extra one beyond
        self.tout = np.int((self.ndays*(24*3600.))/self.tseas + 1)

        # Calculate time outputs stride. Will be 1 if want to use all model output.
        self.tstride = int(self.tseas_use/self.tseas) # will round down

    def _readgrid(self):
        '''
        Read in horizontal and vertical grid.
        '''

        # BREAK UP READGRID INTO SMALLER FUNCTIONS LATER
        # LOOK AT TIMING AGAIN (why is grid so slow) LATER

        # if vertical grid information is not included in the grid file, or if all grid info
        # is not in output file, use two
        if self.grid_filename is not None:
            self.grid = tracpy.inout.readgrid(self.grid_filename, vert_filename=self.currents_filename)
        else:
            self.grid = tracpy.inout.readgrid(self.currents_filename)

    def step(self, tind, nc, j, tstart, ufnew, vfnew, dztnew, zrtnew, zwtnew, 
                xstart, ystart, zstart, T0=None, U=None, V=None):
        '''
        Take some number of steps between a start and end time.
        FIGURE OUT HOW TO KEEP TRACK OF TIME FOR EACH SET OF LINES

        :param tind: Time index to use for stepping
        FILL IN
        '''

        if self.grid is None:
            self._readgrid()

        # Figure out where in time we are 

        ## Load fields if not already loaded

        # Move previous new time step to old time step info
        ufold = ufnew.copy()
        vfold = vfnew.copy()
        dztold = dztnew.copy()
        zrtold = zrtnew.copy()
        zwtold = zwtnew.copy()

        # Read stuff in for next time loop
        if is_string_like(self.z0): # isoslice case
            ufnew,vfnew,dztnew,zrtnew,zwtnew = tracpy.inout.readfields(tind, self.grid, nc, self.z0, self.zpar, zparuv=self.zparuv)
        else: # 3d case
            ufnew,vfnew,dztnew,zrtnew,zwtnew = tracpy.inout.readfields(tind, self.grid, nc)

        #  flux fields at starting time for this step
        # ind = (flag[:] == 0) # indices where the drifters are still inside the domain

        # Combine times for arrays for input to tracmass
        # from [ixjxk] to [ixjxkxt]
        # Change ordering for these three arrays here instead of in readfields since
        # concatenate does not seem to preserve ordering
        uflux = np.asfortranarray(np.concatenate((ufold.reshape(np.append(ufold.shape,1)), \
                                ufnew.reshape(np.append(ufnew.shape,1))), \
                                axis=ufold.ndim))
        vflux = np.asfortranarray(np.concatenate((vfold.reshape(np.append(vfold.shape,1)), \
                                vfnew.reshape(np.append(vfnew.shape,1))), \
                                axis=vfold.ndim))
        dzt = np.asfortranarray(np.concatenate((dztold.reshape(np.append(dztold.shape,1)), \
                                dztnew.reshape(np.append(dztnew.shape,1))), \
                                axis=dztold.ndim))

        # Change the horizontal indices from python to fortran indexing 
        # (vertical are zero-based in tracmass)
        xstart, ystart = tracpy.tools.convert_indices('py2f',xstart,ystart)

        # km that is sent to tracmass is determined from uflux (see tracmass?)
        # so it will be the correct value for whether we are doing the 3D
        # or isoslice case.
        # vec = np.arange(j*N,j*N+N) # indices for storing new track locations
        # tic_tracmass[j] = time.time()
        # pdb.set_trace()
        if self.dostream: # calculate Lagrangian stream functions
            xend,\
                yend,\
                zend,\
                flag,\
                ttend, U, V = \
                    tracmass.step(np.ma.compressed(xstart),
                                    np.ma.compressed(ystart),
                                    np.ma.compressed(zstart),
                                    self.tseas_use, uflux, vflux, self.ff, 
                                    self.grid['kmt'].astype(int), 
                                    dzt, self.grid['dxdy'], self.grid['dxv'], 
                                    self.grid['dyu'], self.grid['h'], self.nsteps, 
                                    self.ah, self.av, self.do3d, self.doturb, self.dostream, self.N, 
                                    t0=T0,
                                    ut=U, vt=V)
        else: # don't calculate Lagrangian stream functions
            xend,\
                yend,\
                zend,\
                flag,\
                ttend, _, _ = \
                    tracmass.step(np.ma.compressed(xstart),
                                    np.ma.compressed(ystart),
                                    np.ma.compressed(zstart),
                                    self.tseas_use, uflux, vflux, self.ff, 
                                    self.grid['kmt'].astype(int), 
                                    dzt, self.grid['dxdy'], self.grid['dxv'], 
                                    self.grid['dyu'], self.grid['h'], self.nsteps, 
                                    self.ah, self.av, self.do3d, self.doturb, self.dostream, self.N)
        # toc_tracmass[j] = time.time()
        # pdb.set_trace()

        # Add initial step time to ttend
        ttend = (ttend.T + tstart).T

        # Change the horizontal indices from python to fortran indexing
        xend, \
            yend \
                            = tracpy.tools.convert_indices('f2py', \
                                xend, \
                                yend)

        
        # Skip calculating real z position if we are doing surface-only drifters anyway
        if self.z0 != 's' and self.zpar != self.grid['km']-1:
            # tic_zinterp[j] = time.time()
            # Calculate real z position
            r = np.linspace(1./self.N,1,self.N) # linear time interpolation constant that is used in tracmass

            for n in xrange(self.N): # loop through time steps
                # interpolate to a specific output time
                # pdb.set_trace()
                zwt = (1.-r[n])*zwtold + r[n]*zwtnew
                zp, dt = tools.interpolate3d(xend, \
                                                        yend, \
                                                        zend, \
                                                        zwt)
            # toc_zinterp[j] = time.time()
        else:
            zp = zend

        # return the new positions or the delta lat/lon
        return ufnew, vfnew, dztnew, zrtnew, zwtnew, xend, yend, zend, zp, flag, ttend, U, V