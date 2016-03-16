#!/usr/bin/env python

'''
TracPy class
'''

import tracpy
import numpy as np
import tracmass
from matplotlib.mlab import find


class Tracpy(object):
    """TracPy class."""

    def __init__(self, currents_filename, grid, nsteps=1, ndays=1, ff=1,
                 tseas=3600., ah=0., av=0., z0='s', zpar=1, do3d=0, doturb=0,
                 name='test', dostream=0, N=1,
                 time_units='seconds since 1970-01-01', dtFromTracmass=None,
                 zparuv=None, tseas_use=None, savell=True, doperiodic=0,
                 usespherical=True, ellps='WGS84'):
        """Initialize class.

        Note: GCM==General Circulation Model, meaning the predicted u/v
        velocity fields that are input into TracPy to run the drifters.

        Args:
            currents_filename (str or List[str]): NetCDF file name
             (with extension), list of file names, or OpenDAP url to GCM
             output.
            grid (object): class containing all necessary grid information,
             as calculated from tracpy.inout.readgrid().
            nsteps (Optional[int]): sets the max time step between GCM model
             outputs between drifter steps. (iter in TRACMASS) Does not
             control the output sampling anymore. The velocity fields are
             assumed frozen while a drifter is stepped through a given grid
             cell. nsteps can force the reinterpolation of the fields by
             setting the max time before reinterpolation. Defaults to 1.
            ndays (Optional[float]): number of days to run for drifter tracks
             from start date. Defaults to 1.
            ff (Optional[int]): 1 is forward in time, -1 is backward.
             Defaults to 1.
            tseas (Optional[float]): number of seconds between GCM model
             outputs. Defaults to 3600.
            ah (Optional[float]): horizontal diffusivity, in m^2/s. Only used
             if doturb!=0. Defaults to 0.
            av (Optional[float]): vertical diffusivity, in m^2/s. Only used
             if doturb!=0 and do3d==1. Defaults to 0.
            z0 (Optional[str or array]): string flag in 2D case or array of
             initial z locations in 3D case. Defaults to 's'.
            zpar (Optional[int or str]): isoslice value to use in 2D case or
             string flag in 3D case. Default is 1.
             For 3D drifter movement, use do3d=1, and z0 should be an array
             of initial drifter depths. The array should be the same size as
             lon0 and be negative for under water. Currently drifter depths
             need to be above the seabed for every x,y particle location for
             the script to run.
             To do 3D but start at surface, use z0=zeros(ia.shape) and have
             either zpar = 'fromMSL'
             choose fromMSL to have z0 starting depths be for that depth
             below the base time-independent sea level (or mean sea level).
             choose 'fromZeta' to have z0 starting depths be for that depth
             below the time-dependent sea surface.
             Haven't quite finished the 'fromZeta' case.
             For 2D drifter movement, turn on twodim flag in makefile.
             Then:
             set z0 to 's' for 2D along a terrain-following slice
             and zpar to be the index of s level you want to use (0 to km-1)
             set z0 to 'rho' for 2D along a density surface
             and zpar to be the density value you want to use
             Can do the same thing with salinity ('salt') or temperature
             ('temp')
             The model output doesn't currently have density though.
             set z0 to 'z' for 2D along a depth slice
             and zpar to be the constant (negative) depth value you want to
             use
             To simulate drifters at the surface, set z0 to 's'
             and zpar = grid['km']-1 to put them in the upper s level
            do3d (Optional[int]): 1 for 3D or 0 for 2D. Default is 0.
            doturb (Optional[int]): 0 for no added diffusion, 1 for diffusion
             via velocity fluctuation, 2/3 for diffusion via random walk (3
             for aligned with isobaths). Default is 0.
            name (Optional[str]): name for output. Default is 'test'.
            dostream (Optional[int]): 1 to calculate transport for lagrangian
             stream functions, 0 to not. Default is 0.
            N (Optional[int]): number of steps between GCM model outputs for
             outputting drifter locations. Defaults to output at nsteps.
             If dtFromTracmass is being used, N is set by that.
            time_units (Optional[int]): Reference for time, for changing
             between numerical times and datetime format. Default is
             'seconds since 1970-01-01'.
            dtFromTracmass: Time period for exiting from TRACMASS. If
             uninitialized, this is set to tseas so that it only exits
             TRACMASS when it has gone through a full model output. If
             initialized by the user, TRACMASS will run for 1 time step of
             length dtFromTracmass before exiting to the loop.
            zparuv: Defaults to zpar. Use this if the k index for the model
             output fields (e.g, u, v) is different from the k index in the
             grid This might happen if, for example, only the surface current
             were saved, but the model run originally did have many layers.
             This parameter represents the k index for the u and v output,
             not for the grid.
            tseas_use: Defaults to tseas. Desired time between outputs in
             seconds, as opposed to the actual time between outputs (tseas).
             Should be >= tseas since this is just an ability to use model
             output at less frequency than is available, probably just for
             testing purposes or matching other models. Should be a multiple
             of tseas (or will be rounded later).
            savell (Optional[bool]): True to save drifter tracks in lon/lat
             and False to save them in grid coords. Default is True.
            doperiodic (Optional[int]): Whether to use periodic boundary
             conditions for drifters and, if so, on which walls. Default is 0.
             0: do not use periodic boundary conditions
             1: use a periodic boundary condition in the east-west/x/i
             direction
             2: use a periodic boundary condition in the north-south/y/j
             direction
            usespherical (Optional[bool]): True if want to use spherical
             (lon/lat) coordinates and False for idealized applications where
             it isn't necessary to project from spherical coordinates.
             Default is True.
        """

        self.currents_filename = currents_filename
        self.grid = grid

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
        self.savell = savell
        self.doperiodic = doperiodic
        self.usespherical = usespherical

        # if loopsteps is None and nsteps is not None:
        #     # Use nsteps in TRACMASS and have inner loop collapse
        #     self.loopsteps = 1
        # elif loopsteps is not None and nsteps is None:
        #     # This means to use the inner loop (with loopsteps) and nsteps=1
        #     # to just do 1 step per call to TRACMASS
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
            # # this is the total number of model_step_is_done
            # self.N = (self.ndays*3600*24.)/self.tseas
            self.dtFromTracmass = dtFromTracmass

        # Find number of interior loop steps in case dtFromTracmass is not
        # equal to tseas
        # NEEDS TO BE EVEN NUMBER FOR NOW: NEED TO GENERALIZE THIS LATER
        self.nsubsteps = int(self.tseas/self.dtFromTracmass)

        if zparuv is None:
            self.zparuv = zpar
        else:
            self.zparuv = zparuv

        if tseas_use is None:
            self.tseas_use = tseas

        # Calculate parameters that derive from other parameters

        # Number of model outputs to use (based on tseas, actual amount of
        # model output)
        # This should not be updated with tstride since it represents the full
        # amount of indices in the original model output. tstride will be used
        # separately to account for the difference.
        # Adding one index so that all necessary indices are captured by this
        # number.
        # Then the run loop uses only the indices determined by tout instead
        # of needing an extra one beyond now rounding up instead of down
        self.tout = np.int(np.ceil((ndays*(24*3600))/tseas + 1))

        # Calculate time outputs stride. Will be 1 if want to use all model
        # output.
        self.tstride = int(self.tseas_use/self.tseas)  # will round down

        # For later use
        # fluxes
        self.uf = None
        self.vf = None
        self.dzt = None
        self.zrt = None
        self.zwt = None

    def prepare_for_model_run(self, date, lon0, lat0):
        """Get everything ready so that we can get to the simulation."""

        # # Convert date to number
        # date = netCDF.date2num(date, self.time_units)

        # Figure out what files will be used for this tracking
        nc, tinds = tracpy.inout.setupROMSfiles(self.currents_filename, date,
                                                self.ff, self.tout,
                                                self.time_units,
                                                tstride=self.tstride)

        # Interpolate to get starting positions in grid space
        # convert from assumed input lon/lat coord locations to grid space
        if self.usespherical:
            xstart0, ystart0, _ = tracpy.tools.interpolate2d(lon0, lat0,
                                                             self.grid,
                                                             'd_ll2ij')
        # assume input seed locations are in projected/idealized space and
        # change to index space
        else:
            xstart0, ystart0, _ = tracpy.tools.interpolate2d(lon0, lat0,
                                                             self.grid,
                                                             'd_xy2ij')
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
        # time at start of drifter test from file in seconds since 1970-01-01
        # add this on at the end since it is big
        t0save = dates[tinds[0]]

        # Initialize drifter grid positions and indices
        xend = np.ones((ia.size, (len(tinds)-1)*self.N+1))*np.nan
        yend = np.ones((ia.size, (len(tinds)-1)*self.N+1))*np.nan
        zend = np.ones((ia.size, (len(tinds)-1)*self.N+1))*np.nan
        zp = np.ones((ia.size, (len(tinds)-1)*self.N+1))*np.nan
        ttend = np.ones((ia.size, (len(tinds)-1)*self.N+1))
        # initialize all exit flags for in the domain
        flag = np.zeros((ia.size), dtype=np.int)

        # Initialize vertical stuff and fluxes
        # Read initial field in - to 'new' variable since will be moved
        # at the beginning of the time loop ahead
        lx = self.grid.imt
        ly = self.grid.jmt
        lk = self.grid.sc_r.size

        # Now that we have the grid, initialize the info for the two
        # bounding model steps using the grid size
        self.uf = np.ones((2, lk-1, ly, lx-1))*np.nan
        self.vf = np.ones((2, lk-1, ly-1, lx))*np.nan
        self.dzt = np.ones((2, lk-1, ly, lx))*np.nan
        self.zrt = np.ones((2, lk-1, ly, lx))*np.nan
        self.zwt = np.ones((2, lk, ly, lx))*np.nan
        if isinstance(self.z0, str):  # isoslice case
            self.uf[1, :, :, :], self.vf[1, :, :, :], \
                self.dzt[1, :, :, :], self.zrt[1, :, :, :], \
                self.zwt[1, :, :, :] = \
                tracpy.inout.readfields(tinds[0], self.grid, nc, self.z0,
                                        self.zpar, zparuv=self.zparuv)

        else:  # 3d case
            self.uf[1, :, :, :], self.vf[1, :, :, :], \
                self.dzt[1, :, :, :], self.zrt[1, :, :, :], \
                self.zwt[1, :, :, :] = \
                tracpy.inout.readfields(tinds[0], self.grid, nc)

        ## Find zstart0 and ka
        # The k indices and z grid ratios should be on a wflux vertical grid,
        # which goes from 0 to km since the vertical velocities are defined
        # at the vertical cell edges. A drifter's grid cell is vertically
        # bounded above by the kth level and below by the (k-1)th level
        if isinstance(self.z0, str):  # then doing a 2d isoslice
            # there is only one vertical grid cell, but with two vertically-
            # bounding edges, 0 and 1, so the initial ka value is 1 for all
            # isoslice drifters.
            ka = np.ones(ia.size)

            # for s level isoslice, place drifters vertically at the center
            # of the grid cell since that is where the u/v flux info is from.
            # For a rho/temp/density isoslice, we treat it the same way, such
            # that the u/v flux info taken at a specific rho/temp/density
            # value is treated as being at the center of the grid cells
            # vertically.
            zstart0 = np.ones(ia.size)*0.5

        else:   # 3d case
            # Convert initial real space vertical locations to grid space
            # first find indices of grid cells vertically
            ka = np.ones(ia.size)*np.nan
            zstart0 = np.ones(ia.size)*np.nan

            if self.zpar == 'fromMSL':
                # print 'zpar==''fromMSL'' not implemented yet...'
                raise NotImplementedError("zpar==''fromMSL'' not implemented\
                                          yet...")
            #     for i in xrange(ia.size):
            #         # pdb.set_trace()
            #         ind = (self.grid['zwt0'][ia[i],ja[i],:]<=self.z0[i])
            #         # check to make sure there is at least one true value, so the z0 is shallower than the seabed
            #         if np.sum(ind):
            #             ka[i] = find(ind)[-1] # find value that is just shallower than starting vertical position
            #         # if the drifter starting vertical location is too deep for the x,y location, complain about it
            #         else:  # Maybe make this nan or something later
            #             print 'drifter vertical starting location is too deep for its x,y location. Try again.'
            #         if (self.z0[i] != self.grid['zwt0'][ia[i],ja[i],ka[i]]) and (ka[i] != self.grid['km']): # check this
            #             ka[i] = ka[i]+1
            #         # Then find the vertical relative position in the grid cell by adding on the bit of grid cell
            #         zstart0[i] = ka[i] - abs(self.z0[i]-self.grid['zwt0'][ia[i],ja[i],ka[i]]) \
            #                             /abs(self.grid['zwt0'][ia[i],ja[i],ka[i]-1]-self.grid['zwt0'][ia[i],ja[i],ka[i]])
            elif self.zpar == 'fromZeta':
                # In this case, the starting z values of the drifters are
                # found in grid space as z0 below the z surface for each
                # drifter
                for i in xrange(ia.size):
                    ind = (self.zwt[1, :, ja[i], ia[i]] <= self.z0[i])
                    # find value that is just shallower than starting vertical
                    # position
                    ka[i] = find(ind)[-1]
                    if (self.z0[i] != self.zwt[1, ka[i], ja[i], ia[i]]) and \
                       (ka[i] != self.grid.km):  # check this
                        ka[i] = ka[i]+1
                    # Then find the vertical relative position in the grid
                    # cell by adding on the bit of grid cell
                    zstart0[i] = ka[i] - \
                        abs(self.z0[i]-self.zwt[1, ka[i], ja[i], ia[i]]) \
                        / abs(self.zwt[1, ka[i]-1, ja[i], ia[i]] -
                              self.zwt[1, ka[i], ja[i], ia[i]])

        # Find initial cell depths to concatenate to beginning of drifter
        # tracks later
        # zsave = tracpy.tools.interpolate3d(xstart0, ystart0, zstart0,
        #                                    self.zwt[:, :, :, 1])

        # Initialize x,y,z with initial seeded positions
        xend[:, 0] = xstart0
        yend[:, 0] = ystart0
        zend[:, 0] = zstart0

        return tinds, nc, t0save, xend, yend, zend, zp, ttend, flag

    def prepare_for_model_step(self, tind, nc, flag, xend, yend, zend, j,
                               nsubstep, T0):
        """Already in a step, get ready to actually do step"""

        xstart = xend[:, j*self.N]
        ystart = yend[:, j*self.N]
        zstart = zend[:, j*self.N]

        # mask out drifters that have exited the domain
        xstart = np.ma.masked_where(flag[:] == 1, xstart)
        ystart = np.ma.masked_where(flag[:] == 1, ystart)
        zstart = np.ma.masked_where(flag[:] == 1, zstart)
        if T0 is not None:
            T0 = np.ma.masked_where(flag[:] == 1, T0)

        # Move previous new time step to old time step info
        self.uf[0, :, :, :] = self.uf[1, :, :, :].copy()
        self.vf[0, :, :, :] = self.vf[1, :, :, :].copy()
        self.dzt[0, :, :, :] = self.dzt[1, :, :, :].copy()
        self.zrt[0, :, :, :] = self.zrt[1, :, :, :].copy()
        self.zwt[0, :, :, :] = self.zwt[1, :, :, :].copy()

        # Read stuff in for next time loop
        if isinstance(self.z0, str):  # isoslice case
            self.uf[1, :, :, :], self.vf[1, :, :, :], self.dzt[1, :, :, :], \
                self.zrt[1, :, :, :], self.zwt[1, :, :, :] = \
                tracpy.inout.readfields(tind, self.grid, nc, self.z0,
                                        self.zpar, zparuv=self.zparuv)
        else:  # 3d case
            self.uf[1, :, :, :], self.vf[1, :, :, :], self.dzt[1, :, :, :], \
                self.zrt[1, :, :, :], self.zwt[1, :, :, :] = \
                tracpy.inout.readfields(tind, self.grid, nc)

        # Find the fluxes of the immediately bounding range for the desired
        # time step, which can be less than 1 model output
        # SHOULD THIS BE PART OF SELF TOO? Leave uf and vf as is, though,
        # because they may be used for interpolating the input fluxes for
        # substeps.
        ufsub = np.ones(self.uf.shape)*np.nan
        vfsub = np.ones(self.vf.shape)*np.nan
        # for earlier bounding flux info
        rp = nsubstep/self.nsubsteps  # weighting for later time step
        rm = 1 - rp  # timing for earlier time step
        ufsub[0, :, :, :] = rm*self.uf[0, :, :, :] + rp*self.uf[1, :, :, :]
        vfsub[0, :, :, :] = rm*self.vf[0, :, :, :] + rp*self.vf[1, :, :, :]
        # for later bounding flux info
        rp = (nsubstep+1)/self.nsubsteps  # weighting for later time step
        rm = 1 - rp  # timing for earlier time step
        ufsub[1, :, :, :] = rm*self.uf[0, :, :, :] + rp*self.uf[1, :, :, :]
        vfsub[1, :, :, :] = rm*self.vf[0, :, :, :] + rp*self.vf[1, :, :, :]

        # Change the horizontal indices from python to fortran indexing
        # (vertical are zero-based in tracmass)
        xstart, ystart = tracpy.tools.convert_indices('py2f', xstart, ystart)

        return xstart, ystart, zstart, ufsub, vfsub, T0

    def step(self, xstart, ystart, zstart, ufsub, vfsub, T0, U, V):
        """
        Take some number of steps between a start and end time.
        FIGURE OUT HOW TO KEEP TRACK OF TIME FOR EACH SET OF LINES

        FILL IN
        Transpose arrays when sent to Fortran.
        """

        if T0 is not None:
            xend, yend, zend, flag,\
                ttend, U, V = tracmass.step(np.ma.compressed(xstart),
                                            np.ma.compressed(ystart),
                                            np.ma.compressed(zstart),
                                            self.tseas_use, ufsub.T, vfsub.T,
                                            self.ff,
                                            self.grid.kmt.astype(int).T,
                                            self.dzt.T, self.grid.dxdy.T,
                                            self.grid.dxv.T,
                                            self.grid.dyu.T, self.grid.h.T,
                                            self.nsteps, self.ah, self.av,
                                            self.do3d, self.doturb,
                                            self.doperiodic, self.dostream,
                                            self.N, t0=np.ma.compressed(T0),
                                            ut=U.T, vt=V.T)
        else:
            xend, yend, zend, flag,\
                ttend, U, V = tracmass.step(np.ma.compressed(xstart),
                                            np.ma.compressed(ystart),
                                            np.ma.compressed(zstart),
                                            self.tseas_use, ufsub.T, vfsub.T,
                                            self.ff,
                                            self.grid.kmt.astype(int).T,
                                            self.dzt.T, self.grid.dxdy.T,
                                            self.grid.dxv.T,
                                            self.grid.dyu.T, self.grid.h.T,
                                            self.nsteps, self.ah, self.av,
                                            self.do3d, self.doturb,
                                            self.doperiodic, self.dostream,
                                            self.N)

        # return the new positions or the delta lat/lon
        return xend, yend, zend, flag, ttend, U, V

    def model_step_is_done(self, xend, yend, zend, ttend, tstart):
        """Stuff to do after a call to TRACMASS"""

        # Add initial step time to ttend
        ttend = (ttend.T + tstart).T

        # Change the horizontal indices from python to fortran indexing
        xend, yend = tracpy.tools.convert_indices('f2py', xend, yend)

        # Skip calculating real z position if we are doing surface-only
        # drifters anyway
        if self.z0 != 's' and self.zpar != self.grid.km-1:

            # Calculate real z position
            # linear time interpolation constant that is used in tracmass
            r = np.linspace(1./self.N, 1, self.N)

            for n in xrange(self.N):  # loop through time steps
                # interpolate to a specific output time
                # pdb.set_trace()
                zwt = (1.-r[n])*self.zwt[0, :, :, :] + \
                    r[n]*self.zwt[1, :, :, :]
                zp, dt = tracpy.tools.interpolate3d(xend, yend, zend, zwt)
        else:
            zp = zend

        # return the new positions or the delta lat/lon
        return xend, yend, zend, zp, ttend

    def finishSimulation(self, ttend, t0save, xend, yend, zp, T0, U, V):
        """
        Wrap up simulation.

        Note:
            Not doing transport yet.
        """

        ttend = ttend + t0save  # add back in base time in seconds

        ## map coordinates interpolation if saving tracks as lon/lat
        if self.savell:
            if self.usespherical:
                lonp, latp, dt = tracpy.tools.interpolate2d(xend, yend,
                                                            self.grid,
                                                            'm_ij2ll',
                                                            mode='constant',
                                                            cval=np.nan)
            else:
                lonp, latp, dt = tracpy.tools.interpolate2d(xend, yend,
                                                            self.grid,
                                                            'm_ij2xy',
                                                            mode='constant',
                                                            cval=np.nan)
        else:
            # rename grid index locations as lon/lat to fit in with save
            # syntax below
            lonp = xend
            latp = yend

        # Save results to netcdf file
        tracpy.inout.savetracks(lonp, latp, zp, ttend, self.name, self.nsteps,
                                self.N, self.ff, self.tseas_use, self.ah,
                                self.av, self.do3d, self.doturb,
                                self.currents_filename, self.doperiodic,
                                self.time_units, T0, U, V, savell=self.savell)

        return lonp, latp, zp, ttend, T0, U, V
