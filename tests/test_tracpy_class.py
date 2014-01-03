'''
Testing TracPy class
Call with py.test test_tracpy_class.py
'''

import tracpy
from tracpy.tracpy_class import Tracpy
import os

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
