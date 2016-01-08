'''
Testing TracPy class
Call with py.test test_tracpy_class.py
'''

import tracpy
from tracpy.tracpy_class import Tracpy
import os
import datetime
import numpy as np

# For niceties with file locations and such
here = os.path.dirname(__file__)


def test_init():
    '''
    Initialize class without vertical grid information.

    Make sure that we can initialize the class with a grid and without
    vertical grid information.
    '''

    # Get projection object
    proj = tracpy.tools.make_proj(setup='galveston', usebasemap=False)

    # Read in grid
    grid_filename = os.path.join(here, 'input', 'grid.nc')
    currents_filename = os.path.join(here, 'input', 'ocean_his_0001.nc')
    grid = tracpy.inout.readgrid(grid_filename, proj)

    tp = Tracpy(currents_filename, grid)

    assert True

    assert tp.grid

    # make sure there isn't a value for a vertical grid parameter
    assert not hasattr(grid, 'sc_r')


def test_readgridWithVertical():
    '''
    Test for initializing grid and saving vertical grid information.
    '''

    # Get projection object
    proj = tracpy.tools.make_proj(setup='galveston', usebasemap=False)

    # Read in grid
    grid_filename = os.path.join(here, 'input', 'grid.nc')
    currents_filename = os.path.join(here, 'input', 'ocean_his_0001.nc')
    grid = tracpy.inout.readgrid(grid_filename, proj,
                                 vert_filename=currents_filename)

    tp = Tracpy(currents_filename, grid)

    assert True

    assert tp.grid

    # make sure there is a value for a vertical grid parameter
    assert hasattr(grid, 'sc_r')


def test_prepareForSimulation():
    '''
    Test initial steps to get ready for a simulation.
    '''

    # Get projection object
    proj = tracpy.tools.make_proj(setup='galveston', usebasemap=False)

    # Read in grid
    grid_filename = os.path.join(here, 'input', 'grid.nc')
    currents_filename = os.path.join(here, 'input', 'ocean_his_0001.nc')
    grid = tracpy.inout.readgrid(grid_filename, proj,
                                 vert_filename=currents_filename)

    tp = Tracpy(currents_filename, grid)

    date = datetime.datetime(2013, 12, 19, 0)
    lon0 = [-123., -123.]
    lat0 = [48.55, 48.75]

    tinds, nc, t0save, xend, yend, zend, zp, ttend, flag \
        = tp.prepare_for_model_run(date, lon0, lat0)

    assert True

    assert tp.uf is not None

    assert np.sum(np.isnan(tp.uf[0, :, :, :])) == tp.uf[0, :, :, :].size
