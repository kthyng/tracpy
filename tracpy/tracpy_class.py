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
import datetime
import netCDF4 as netCDF

class Tracpy(object):
    '''
    TracPy class.
    '''

    def __init__(self, currents_filename, grid_filename=None, nsteps=1, ndays=1, ff=1, tseas=3600.,
                ah=0., av=0., z0='s', zpar=1, do3d=0, doturb=0, name='test', dostream=0, N=1, 
                time_units='seconds since 1970-01-01', dtFromTracmass=None, zparuv=None, tseas_use=None,
                T0=None, U=None, V=None, usebasemap=False, savell=True, doperiodic=0):
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
        :param T0=None: Volume transport represented by each drifter. for use with dostream=1.
        :param U=None: east-west transport, is updated by TRACMASS. Only used if dostream=1.
        :param V=None: north-south transport, is updated by TRACMASS. Only used if dostream=1.
        :param usebasemap=False: whether to use basemap for projections in readgrid or not. Not is faster, but using basemap allows for plotting.
        :param savell=True: True to save drifter tracks in lon/lat and False to save them in grid coords
        :param doperiodic=0: Whether to use periodic boundary conditions for drifters and, if so, on which walls.
                0: do not use periodic boundary conditions
                1: use a periodic boundary condition in the east-west/x/i direction
                2: use a periodic boundary condition in the north-south/y/j direction
        '''

        self.currents_filename = currents_filename
        self.grid_filename = grid_filename
        self.grid = None

        # Initial parameters
        self.nsteps = nsteps
        self.ndays = ndays
        self.ff = ff
        self.tseas = float(tseas)
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
        self.T0 = T0
        self.U = U
        self.V = V
        self.usebasemap = usebasemap
        self.savell = savell
        self.doperiodic = doperiodic

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
            # self.N = (self.ndays*3600*24.)/self.tseas # this is the total number of model_step_is_done
            self.dtFromTracmass = dtFromTracmass

        # Find number of interior loop steps in case dtFromTracmass is not equal to tseas
        # NEEDS TO BE EVEN NUMBER FOR NOW: NEED TO GENERALIZE THIS LATER
        self.nsubsteps = int(self.tseas/self.dtFromTracmass)

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
        # now rounding up instead of down
        self.tout = np.int(np.ceil((ndays*(24*3600))/tseas + 1))

        # Calculate time outputs stride. Will be 1 if want to use all model output.
        self.tstride = int(self.tseas_use/self.tseas) # will round down

        # U,V currently initialized elsewhere if not input explicitly.
        if dostream:
            assert self.T0 is not None

        # For later use
        # fluxes
        self.uf = None
        self.vf = None
        self.dzt = None
        self.zrt = None
        self.zwt = None

    def _readgrid(self):
        '''
        Read in horizontal and vertical grid.
        '''

        # BREAK UP READGRID INTO SMALLER FUNCTIONS LATER

        # if vertical grid information is not included in the grid file, or if all grid info
        # is not in output file, use two
        if self.grid_filename is not None:
            self.grid = tracpy.inout.readgrid(self.grid_filename, vert_filename=self.currents_filename, usebasemap=self.usebasemap)
        else:
            self.grid = tracpy.inout.readgrid(self.currents_filename, usebasemap=self.usebasemap)

    def prepare_for_model_run(self, date, lon0, lat0):
        '''
        Get everything ready so that we can get to the simulation.

        FILL IN
        '''

    #     self.initialize_time()
    #     self.setup_initial_velocities()


        # Convert date to number
        date = netCDF.date2num(date, self.time_units)

        # Figure out what files will be used for this tracking
        nc, tinds = tracpy.inout.setupROMSfiles(self.currents_filename, date, self.ff, self.tout, tstride=self.tstride)

        # Read in grid parameters into dictionary, grid, if haven't already
        if self.grid is None:
            self._readgrid()

        # If dostream==1, do transport calculations and initialize to an empty array
        if self.U is None and self.dostream:
            self.U = np.ma.zeros(grid['xu'].shape, order='F')
            self.V = np.ma.zeros(grid['xv'].shape, order='F')

        # Interpolate to get starting positions in grid space
        xstart0, ystart0, _ = tracpy.tools.interpolate2d(lon0, lat0, self.grid, 'd_ll2ij')
        # Do z a little lower down

        # Initialize seed locations 
        ia = np.ceil(xstart0)
        ja = np.ceil(ystart0)

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
        xend = np.ones((ia.size,(len(tinds)-1)*self.N+1))*np.nan
        yend = np.ones((ia.size,(len(tinds)-1)*self.N+1))*np.nan
        zend = np.ones((ia.size,(len(tinds)-1)*self.N+1))*np.nan
        zp = np.ones((ia.size,(len(tinds)-1)*self.N+1))*np.nan
        # iend = np.ones((ia.size,(len(tinds)-1)*self.N))*np.nan
        # jend = np.ones((ia.size,(len(tinds)-1)*self.N))*np.nan
        # kend = np.ones((ia.size,(len(tinds)-1)*self.N))*np.nan
        ttend = np.zeros((ia.size,(len(tinds)-1)*self.N+1))
        flag = np.zeros((ia.size),dtype=np.int) # initialize all exit flags for in the domain

        # Initialize vertical stuff and fluxes
        # Read initial field in - to 'new' variable since will be moved
        # at the beginning of the time loop ahead
        if is_string_like(self.z0): # isoslice case
            # Now that we have the grid, initialize the info for the two bounding model 
            # steps using the grid size
            self.uf = np.asfortranarray(np.ones((self.grid['xu'].shape[0], 
                                                self.grid['xu'].shape[1], 1, 2)))*np.nan
            self.vf = np.asfortranarray(np.ones((self.grid['xv'].shape[0], 
                                                self.grid['xv'].shape[1], 1, 2)))*np.nan
            self.dzt = np.asfortranarray(np.ones((self.grid['xr'].shape[0], 
                                                self.grid['xr'].shape[1], 1, 2)))*np.nan
            self.zrt = np.asfortranarray(np.ones((self.grid['xr'].shape[0], 
                                                self.grid['xr'].shape[1], 1, 2)))*np.nan
            self.zwt = np.asfortranarray(np.ones((self.grid['xr'].shape[0], 
                                                self.grid['xr'].shape[1], self.grid['sc_r'].size, 2)))*np.nan
            self.uf[:,:,:,1], self.vf[:,:,:,1], self.dzt[:,:,:,1], self.zrt[:,:,:,1], self.zwt[:,:,:,1] = tracpy.inout.readfields(tinds[0],self.grid,nc,self.z0,self.zpar,zparuv=self.zparuv)
            # self.ufnew,self.vfnew,self.dztnew,self.zrtnew,self.zwtnew = tracpy.inout.readfields(tinds[0],self.grid,nc,self.z0,self.zpar,zparuv=self.zparuv)
        else: # 3d case
            # Now that we have the grid, initialize the info for the two bounding model 
            # steps using the grid size
            self.uf = np.asfortranarray(np.ones((self.grid['xu'].shape[0], 
                                                self.grid['xu'].shape[1],
                                                self.grid['sc_r'].size,
                                                2)))*np.nan
            self.vf = np.asfortranarray(np.ones((self.grid['xu'].shape[0], 
                                                self.grid['xu'].shape[1],
                                                self.grid['sc_r'].size,
                                                2)))*np.nan
            self.dzt = np.asfortranarray(np.ones((self.grid['xu'].shape[0], 
                                                self.grid['xu'].shape[1],
                                                self.grid['sc_r'].size,
                                                2)))*np.nan
            self.zrt = np.asfortranarray(np.ones((self.grid['xu'].shape[0], 
                                                self.grid['xu'].shape[1],
                                                self.grid['sc_r'].size,
                                                2)))*np.nan
            self.zwt = np.asfortranarray(np.ones((self.grid['xu'].shape[0], 
                                                self.grid['xu'].shape[1],
                                                self.grid['sc_r'].size,
                                                2)))*np.nan
            self.uf[:,:,:,1], self.vf[:,:,:,1], self.dzt[:,:,:,1], self.zrt[:,:,:,1], self.zwt[:,:,:,1] = tracpy.inout.readfields(tinds[0],self.grid,nc)
        # pdb.set_trace()
        ## Find zstart0 and ka
        # The k indices and z grid ratios should be on a wflux vertical grid,
        # which goes from 0 to km since the vertical velocities are defined
        # at the vertical cell edges. A drifter's grid cell is vertically bounded
        # above by the kth level and below by the (k-1)th level
        if is_string_like(self.z0): # then doing a 2d isoslice
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

            if self.zpar == 'fromMSL':
                for i in xrange(ia.size):
                    # pdb.set_trace()
                    ind = (self.grid['zwt0'][ia[i],ja[i],:]<=self.z0[i])
                    # check to make sure there is at least one true value, so the z0 is shallower than the seabed
                    if np.sum(ind): 
                        ka[i] = find(ind)[-1] # find value that is just shallower than starting vertical position
                    # if the drifter starting vertical location is too deep for the x,y location, complain about it
                    else:  # Maybe make this nan or something later
                        print 'drifter vertical starting location is too deep for its x,y location. Try again.'
                    if (self.z0[i] != self.grid['zwt0'][ia[i],ja[i],ka[i]]) and (ka[i] != self.grid['km']): # check this
                        ka[i] = ka[i]+1
                    # Then find the vertical relative position in the grid cell by adding on the bit of grid cell
                    zstart0[i] = ka[i] - abs(self.z0[i]-self.grid['zwt0'][ia[i],ja[i],ka[i]]) \
                                        /abs(self.grid['zwt0'][ia[i],ja[i],ka[i]-1]-self.grid['zwt0'][ia[i],ja[i],ka[i]])
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
        zsave = tracpy.tools.interpolate3d(xstart0, ystart0, zstart0, self.zwt[:,:,:,1])

        # Initialize x,y,z with initial seeded positions
        xend[:,0] = xstart0
        yend[:,0] = ystart0
        zend[:,0] = zstart0

        # toc_initial = time.time()

        # # j = 0 # index for number of saved steps for drifters
        # tic_read = np.zeros(len(tinds))
        # toc_read = np.zeros(len(tinds))
        # tic_zinterp = np.zeros(len(tinds))
        # toc_zinterp = np.zeros(len(tinds))
        # tic_tracmass = np.zeros(len(tinds))
        # toc_tracmass = np.zeros(len(tinds))
        # pdb.set_trace()
        # xr3 = self.grid['xr'].reshape((self.grid['xr'].shape[0],self.grid['xr'].shape[1],1)).repeat(zwtnew.shape[2],axis=2)
        # yr3 = self.grid['yr'].reshape((self.grid['yr'].shape[0],self.grid['yr'].shape[1],1)).repeat(zwtnew.shape[2],axis=2)

        # def initialize_time(self):

        return tinds, nc, t0save, xend, yend, zend, zp, ttend, flag

    def prepare_for_model_step(self, tind, nc, flag, xend, yend, zend, j, nsubstep):
        '''
        Already in a step, get ready to actually do step
        '''

        xstart = xend[:,j*self.N]
        ystart = yend[:,j*self.N]
        zstart = zend[:,j*self.N]

        # mask out drifters that have exited the domain
        xstart = np.ma.masked_where(flag[:]==1,xstart)
        ystart = np.ma.masked_where(flag[:]==1,ystart)
        zstart = np.ma.masked_where(flag[:]==1,zstart)

        # Move previous new time step to old time step info
        self.uf[:,:,:,0] = self.uf[:,:,:,1].copy()
        self.vf[:,:,:,0] = self.vf[:,:,:,1].copy()
        self.dzt[:,:,:,0] = self.dzt[:,:,:,1].copy()
        self.zrt[:,:,:,0] = self.zrt[:,:,:,1].copy()
        self.zwt[:,:,:,0] = self.zwt[:,:,:,1].copy()

        # Read stuff in for next time loop
        if is_string_like(self.z0): # isoslice case
            self.uf[:,:,:,1],self.vf[:,:,:,1],self.dzt[:,:,:,1],self.zrt[:,:,:,1],self.zwt[:,:,:,1] = tracpy.inout.readfields(tind, self.grid, nc, self.z0, self.zpar, zparuv=self.zparuv)
        else: # 3d case
            self.uf[:,:,:,1],self.vf[:,:,:,1],self.dzt[:,:,:,1],self.zrt[:,:,:,1],self.zwt[:,:,:,1] = tracpy.inout.readfields(tind, self.grid, nc)

        # Find the fluxes of the immediately bounding range for the desired time step, which can be less than 1 model output
        # SHOULD THIS BE PART OF SELF TOO? Leave uf and vf as is, though, because they may be used for interpolating the
        # input fluxes for substeps.
        ufsub = np.ones(self.uf.shape)*np.nan
        vfsub = np.ones(self.vf.shape)*np.nan
        # for earlier bounding flux info
        rp = nsubstep/self.nsubsteps # weighting for later time step
        rm = 1 - rp # timing for earlier time step
        ufsub[:,:,:,0] = rm*self.uf[:,:,:,0] + rp*self.uf[:,:,:,1]
        vfsub[:,:,:,0] = rm*self.vf[:,:,:,0] + rp*self.vf[:,:,:,1]
        # for later bounding flux info
        rp = (nsubstep+1)/self.nsubsteps # weighting for later time step
        rm = 1 - rp # timing for earlier time step
        ufsub[:,:,:,1] = rm*self.uf[:,:,:,0] + rp*self.uf[:,:,:,1]
        vfsub[:,:,:,1] = rm*self.vf[:,:,:,0] + rp*self.vf[:,:,:,1]

        # Change the horizontal indices from python to fortran indexing 
        # (vertical are zero-based in tracmass)
        xstart, ystart = tracpy.tools.convert_indices('py2f',xstart,ystart)

        return xstart, ystart, zstart, ufsub, vfsub

    def step(self, xstart, ystart, zstart, ufsub, vfsub):
        '''
        Take some number of steps between a start and end time.
        FIGURE OUT HOW TO KEEP TRACK OF TIME FOR EACH SET OF LINES

        :param tind: Time index to use for stepping
        FILL IN
        '''

        # Figure out where in time we are 

        # tic_tracmass[j] = time.time()
        xend, yend, zend, flag,\
            ttend, U, V = \
                tracmass.step(np.ma.compressed(xstart),
                                np.ma.compressed(ystart),
                                np.ma.compressed(zstart),
                                self.tseas_use, ufsub, vfsub, self.ff, 
                                self.grid['kmt'].astype(int), 
                                self.dzt, self.grid['dxdy'], self.grid['dxv'], 
                                self.grid['dyu'], self.grid['h'], self.nsteps, 
                                self.ah, self.av, self.do3d, self.doturb, 
                                self.doperiodic, self.dostream, self.N, 
                                t0=self.T0,
                                ut=self.U, vt=self.V)
        # toc_tracmass[j] = time.time()
        # pdb.set_trace()

        # return the new positions or the delta lat/lon
        return xend, yend, zend, flag, ttend, U, V

    def model_step_is_done(self, xend, yend, zend, ttend, tstart):
        '''
        Stuff to do after a call to TRACMASS
        '''

        # Add initial step time to ttend
        ttend = (ttend.T + tstart).T

        # Change the horizontal indices from python to fortran indexing
        xend, yend = tracpy.tools.convert_indices('f2py', xend, yend)

        # Skip calculating real z position if we are doing surface-only drifters anyway
        if self.z0 != 's' and self.zpar != self.grid['km']-1:
            # tic_zinterp[j] = time.time()
            # Calculate real z position
            r = np.linspace(1./self.N,1,self.N) # linear time interpolation constant that is used in tracmass

            for n in xrange(self.N): # loop through time steps
                # interpolate to a specific output time
                # pdb.set_trace()
                zwt = (1.-r[n])*self.zwt[:,:,:,0] + r[n]*self.zwt[:,:,:,1]
                zp, dt = tools.interpolate3d(xend, yend, zend, zwt)
            # toc_zinterp[j] = time.time()
        else:
            zp = zend

        # return the new positions or the delta lat/lon
        return xend, yend, zend, zp, ttend

    def finishSimulation(self, ttend, t0save, xend, yend, zp):
        '''
        Wrap up simulation.
        FILL IN
        NOT DOING TRANSPORT YET
        '''

        ttend = ttend + t0save # add back in base time in seconds

        ## map coordinates interpolation
        # xp2, yp2, dt = tools.interpolate(xg,yg,grid,'m_ij2xy')
        # tic = time.time()

        ## map coordinates interpolation if saving tracks as lon/lat
        if self.savell:
            lonp, latp, dt = tracpy.tools.interpolate2d(xend, yend, self.grid, 'm_ij2ll', mode='constant', cval=np.nan)
        else:
            # rename grid index locations as lon/lat to fit in with save syntax below
            lonp = xg; latp = yg;

        # print '2d interp time=', time.time()-tic
        # pdb.set_trace()

        # runtime = time.time()-tic_start


        # print "============================================="
        # print ""
        # print "Simulation name: ", self.name
        # print ""
        # print "============================================="
        # print "run time:\t\t\t", runtime
        # print "---------------------------------------------"
        # print "Time spent on:"

        # initialtime = toc_initial-tic_initial
        # print "\tInitial stuff: \t\t%4.2f (%4.2f%%)" % (initialtime, (initialtime/runtime)*100)

        # readtime = np.sum(toc_read-tic_read)
        # print "\tReading in fields: \t%4.2f (%4.2f%%)" % (readtime, (readtime/runtime)*100)

        # zinterptime = np.sum(toc_zinterp-tic_zinterp)
        # print "\tZ interpolation: \t%4.2f (%4.2f%%)" % (zinterptime, (zinterptime/runtime)*100)

        # tractime = np.sum(toc_tracmass-tic_tracmass)
        # print "\tTracmass: \t\t%4.2f (%4.2f%%)" % (tractime, (tractime/runtime)*100)
        # print "============================================="

        # Save results to netcdf file
        tracpy.inout.savetracks(lonp, latp, zp, ttend, self.name, self.nsteps, self.N, self.ff, self.tseas_use, self.ah, self.av, \
                            self.do3d, self.doturb, self.currents_filename, self.T0, self.U, self.V, savell=self.savell)
        return lonp, latp, zp, ttend, self.grid, self.T0, self.U, self.V

