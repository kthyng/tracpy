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

class Tracpy(object):
    '''
    TracPy class.
    '''

    def __init__(self, loc_source, grid_source=None):
        '''
        Initialize class.

        :param loc_source: NetCDF file name (with extension) or OpenDAP url.
        :param grid_source=None: NetCDF grid file name or OpenDAP url.
        '''

        self.loc_source = loc_source
        self.grid_source = grid_source
        self.grid = None

    def _readgrid(self):
        '''
        Read in horizontal and vertical grid.
        '''

        # BREAK UP READGRID INTO SMALLER FUNCTIONS LATER
        # LOOK AT TIMING AGAIN (why is grid so slow) LATER

        # if vertical grid information is not included in the grid file, or if all grid info
        # is not in output file, use two
        if self.grid_source is not None:
            self.grid = tracpy.inout.readgrid(self.grid_source, vert_source=self.loc_source)
        else:
            self.grid = tracpy.inout.readgrid(self.loc_source)

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

