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
    Make sure that we can initialize the class.
    '''
    
    tp = Tracpy(os.path.join(here, 'input', 'ocean_his_0001.nc'))

    assert True

def test_initWithGridSource():
    '''
    Make sure we can initialize with a grid source location
    '''
    
    tp = Tracpy(os.path.join(here, 'input', 'ocean_his_0001.nc'), grid_filename=os.path.join(here, 'input', 'grid.nc'))

    assert True

def test_readgrid():
    '''
    Test for initializing grid without vertical grid information.
    '''

    tp = Tracpy(os.path.join(here, 'input', 'grid.nc'))
    tp._readgrid()

    assert True

    assert tp.grid

def test_readgridWithVertical():
    '''
    Test for initializing grid and saving vertical grid information.
    '''

    tp = Tracpy(os.path.join(here, 'input', 'ocean_his_0001.nc'),
                os.path.join(here, 'input', 'grid.nc'))
    tp._readgrid()

    assert True

    assert tp.grid

    assert tp.grid['Vtransform'] # make sure vertical info is in there

def test_prepareForSimulation():
    '''
    Test initial steps to get ready for a simulation.
    '''

    date = datetime.datetime(2013, 12, 19, 0)
    lon0 = [-123., -123.]
    lat0 = [48.55, 48.75]
    tp = Tracpy(os.path.join(here, 'input', 'ocean_his_0001.nc'), grid_filename=os.path.join(here, 'input', 'grid.nc'))

    tinds, nc, t0save, xend, yend, zend, zp, ttend, flag = tp.prepare_for_model_run(date, lon0, lat0)

    assert True

    assert tp.uf is not None

    assert np.sum(np.isnan(tp.uf[:,:,:,0])) == tp.uf[:,:,:,0].size

def test_timestep():
    '''
    Test for moving between time indices and datetime.
    FILL IN
    '''

    pass