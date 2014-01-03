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

class Tracpy(object):
    '''
    TracPy class.
    '''

    def __init__(self, currents_filename, grid_filename=None, nsteps=1, ndays=1, ff=1, tseas=3600.,
                ah=0., av=0., z0='s', zpar=1, do3d=0, doturb=0, name='test', dostream=0, N=1, 
                time_units='seconds since 1970-01-01', zparuv=None, tseas_use=None):
        '''
        Initialize class.

        :param currents_filename: NetCDF file name (with extension) or OpenDAP url.
        :param grid_filename=None: NetCDF grid file name or OpenDAP url.
        :param nsteps: number of linearly interpolated steps between model outputs.
        :param ndays: number of run days
        :param ff: 1 is forward in time, -1 is backward
        :param tseas: number of seconds between model outputs
        :param ah: horizontal diffusivity, in m^2/s
        :param av: vertical diffusivity, in m^2/s
        :param z0: string flag in 2D case or array of initial z locations in 3D case
        :param zpar: isoslice value to in 2D case or string flag in 3D case
        :param do3d: 1 for 3D or 0 for 2D
        :param doturb: 0 for no added diffusion, 1 for diffusion vs velocity fluctuation, 2/3 for diffusion via random walk (3 for aligned with isobaths)
        :param name: name for output
        :param dostream: 1 to calculate transport for lagrangian stream functions, 0 to not
        :param N: number of steps between model outputs for outputting drifter locations
        :param time_units: Reference for time, for changing between numerical times and datetime format
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

    def step():
        '''
        Take some number of steps between a start and end time.
        '''

        if self.grid is None:
            self._readgrid()

        # Figure out where in time we are 

        # Load fields if not already loaded

        # Interpolate the velocity fields

        # Call TRACMASS

