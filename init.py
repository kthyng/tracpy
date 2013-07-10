'''
This file in tracpy is an example init file.

Functions to initialize various numerical experiments.
Contains:
	test1
	test2
	galveston
	hab1b

Make a new init_* for your application.

loc 	Path to directory of grid and output files
nsteps 	Number of steps to do between model outputs (iter in tracmass)
ndays 	number of days to track the particles from start date
ff 		ff=1 to go forward in time and ff=-1 for backward in time
date 	Start date in datetime object
tseas	Time between outputs in seconds
ah 		Horizontal diffusion in m^2/s. 
		See project values of 350, 100, 0, 2000. For -turb,-diffusion
av 		Vertical diffusion in m^2/s.
do3d 	for 3d flag, do3d=0 makes the run 2d and do3d=1 makes the run 3d
doturb	turbulence/diffusion flag. 
		doturb=0 means no turb/diffusion,
		doturb=1 means adding parameterized turbulence
		doturb=2 means adding diffusion on a circle
		doturb=3 means adding diffusion on an ellipse (anisodiffusion)
lon0 	Drifter starting locations in x/zonal direction.
lat0 	Drifter starting locations in y/meridional direction.
z0/zpar For 3D drifter movement, turn off twodim flag in makefile.
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
xp 		x-locations in x,y coordinates for drifters
yp 		y-locations in x,y coordinates for drifters
zp 		z-locations (depths from mean sea level) for drifters
t 		time for drifter tracks
name	Name of simulation to be used for netcdf file containing final tracks

'''

import numpy as np
import os
import netCDF4 as netCDF
import pdb
import glob
from datetime import datetime, timedelta
from matplotlib.mlab import *
import inout
import tools

def galveston():
	'''
	Start drifters outside Galveston Bay and see where they move backward in time.

	'''

	# Location of TXLA model output
	if 'rainier' in os.uname():
		loc = '/Users/kthyng/Documents/research/postdoc/' # for model outputs
	elif 'hafen.tamu.edu' in os.uname():
		loc = '/home/kthyng/shelf/' # for model outputs

	# Initialize parameters
	nsteps = 10
	ndays = 2
	ff = -1
	# Start date
	date = datetime(2009,11, 30, 0)
	# Time between outputs
	# Dt = 14400. # in seconds (4 hours), nc.variables['dt'][:] 
	tseas = 4*3600 # 4 hours between outputs, in seconds, time between model outputs 
	ah = 100.
	av = 1.e-5 # m^2/s, or try 5e-6

	## Input starting locations as real space lon,lat locations
	lon0,lat0 = np.meshgrid(np.linspace(-95.3,-94.3,10), 
							np.linspace(28.6,29.6,10))
	# pdb.set_trace()
	lon0 = lon0.flatten()
	lat0 = lat0.flatten()

	## Choose method for vertical placement of drifters
	# Also update makefile accordingly. Choose the twodim flag for isoslice.
	# See above for more notes, but do the following two lines for an isoslice
	z0 = 's' #'z' #'salt' #'s' 
	zpar = 29 #-10 #grid['km']-1 # 30 #grid['km']-1
	# Do the following two for a 3d simulation
	# z0 = np.ones(xstart0.shape)*-40 #  below the surface
	# zpar = 'fromMSL' 
	# pdb.set_trace()

	## Set flags
	
	# for 3d flag, do3d=0 makes the run 2d and do3d=1 makes the run 3d
	do3d = 0
	# turbulence/diffusion flag. doturb=0 means no turb/diffusion,
	# doturb=1 means adding parameterized turbulence
	# doturb=2 means adding diffusion on a circle
	# doturb=3 means adding diffusion on an ellipse (anisodiffusion)
	doturb = 3

	# simulation name, used for saving results into netcdf file
	name = 'galveston'

	return loc,nsteps,ndays,ff,date,tseas,ah,av,lon0,lat0,z0,zpar,do3d,doturb,name

def test1(loc=None, nsteps=None, ff=None, ah=None, grid=None, nlon=None, nlat=None, doturb=None, name=None):
	'''
	A drifter test using TXLA model output. 
	The comparison case for this simulation is 2D (do3d=0) 
	with no turbulence/diffusion (doturb=0).
	Drifters are started at the surface and run forward
	for ten days (ndays=10) from 11/25/09 (in date). Compare results with figure in examples/test1.png.

	Optional inputs for making tests easy to run:
		loc 			'thredds' or 'local', default = 'thredds'
		nsteps 			Number of particle steps to record between model outputs
						Default = 5
		ff 				Backward (-1) or forward (1) in time. Default is forward (1).
		ah 				Horizontal viscosity, default = 5
		grid 			If input, will not redo this step. Default is to load in grid.
		nlon, nlat 		Number of drifters to use in the lon/lat direction in seed array
						Default = 110, 98 (10 km spacing)
		doturb 			What, if any, subgrid parameterization to use. Default is 'none'
		name 			Specific name for track and figure files. Default is 'temp'
	'''

	# Location of TXLA model output
	# file and then grid. 
	# 0150 file goes from (2009, 11, 19, 12, 0) to (2009, 12, 6, 0, 0)
	if loc is None or loc == 'thredds':
		loc = 'http://barataria.tamu.edu:8080/thredds/dodsC/NcML/txla_nesting6.nc'		
	elif loc is 'local':
	# Location of TXLA model output
		if 'rainier' in os.uname():
			loc = '/Users/kthyng/Documents/research/postdoc/' # for model outputs
		elif 'hafen.tamu.edu' in os.uname():
			loc = '/home/kthyng/shelf/' # for model outputs

	# Initialize parameters
	if nsteps is None:
		nsteps = 5
	else:
		nsteps = nsteps

	ndays = .5 #1 #16
	if ff is None:
		ff = 1
	else:
		ff = ff
	# Start date
	date = datetime(2009,11, 25, 0)
	# date = datetime(2009,11, 20, 0)

	# Time between outputs
	# Dt = 14400. # in seconds (4 hours), nc.variables['dt'][:] 
	tseas = 4*3600 # 4 hours between outputs, in seconds, time between model outputs 
	if ah is None:
		ah = 5. #100.
	else:
		ah = ah
	
	av = 1.e-5 # m^2/s, or try 5e-6

	# grid = netCDF.Dataset(loc+'grid.nc')
	# lonr = grid.variables['lon_rho'][:]
	# latr = grid.variables['lat_rho'][:]
	if grid is None:
		grid = inout.readgrid(loc)
	else:
		grid = grid

	## Input starting locations as real space lon,lat locations
	# lon0,lat0 = np.meshgrid(-95.498218005315309,23.142258627126882) # [0,0] (SE) corner
	# lon0,lat0 = np.meshgrid(-97.748582291691989,23.000027311710628) # [-1,0] (SW) corner
	# lon0,lat0 = np.meshgrid(-87.757124031927574,29.235771320764623) # [0,-1] (NE) corner
	# lon0,lat0 = np.meshgrid(-88.3634073986196,30.388542615201313) # [-1,-1] (NW) corner
	# lon0,lat0 = np.meshgrid(np.linspace(-94,-93,10),np.linspace(28,29,10)) # grid outside Galveston Bay
	# lon0,lat0 = np.meshgrid(np.linspace(-95,-91,100),np.linspace(28,29,50)) # rectangle outside Galveston

	# lon0,lat0 = np.meshgrid(np.linspace(-98.5,-87.5,1100),np.linspace(22.5,31,980)) # whole domain, 1 km
	# lon0,lat0 = np.meshgrid(np.linspace(-98.5,-87.5,220),np.linspace(22.5,31,196)) # whole domain, 5 km
	# # FOR TEST1:
	# lon0,lat0 = np.meshgrid(np.linspace(-98.5,-87.5,110),np.linspace(22.5,31,98)) # whole domain, 10 km
	# lon0,lat0 = np.meshgrid(np.linspace(-98.5,-87.5,21),np.linspace(22.5,31,20)) # whole domain, 50 km

	if nlon is None:
		nlon = 110
	else:
		nlon = nlon
	if nlat is None:
		nlat = 98
	else:
		nlat = nlat
	lon0,lat0 = np.meshgrid(np.linspace(-98.5,-87.5,nlon),np.linspace(22.5,31,nlat)) # whole domain, 10 km

	# Eliminate points that are outside domain or in masked areas
	lon0,lat0 = tools.check_points(lon0,lat0,grid)

	## Choose method for vertical placement of drifters
	# Also update makefile accordingly. Choose the twodim flag for isoslice.
	# See above for more notes, but do the following two lines for an isoslice
	z0 = 's'  #'salt' #'s' #'z' #'salt' #'s' 
	zpar = 29 #30 #29 #-10 #grid['km']-1 # 30 #grid['km']-1
	# Do the following two for a 3d simulation
	# z0 = np.ones(xstart0.shape)*-40 #  below the surface
	# zpar = 'fromMSL' 
	# pdb.set_trace()

	# for 3d flag, do3d=0 makes the run 2d and do3d=1 makes the run 3d
	do3d = 0
	# turbulence/diffusion flag. doturb=0 means no turb/diffusion,
	# doturb=1 means adding parameterized turbulence
	# doturb=2 means adding diffusion on a circle
	# doturb=3 means adding diffusion on an ellipse (anisodiffusion)
	if doturb is None:
		doturb = 0
	else:
		doturb = doturb

	# simulation name, used for saving results into netcdf file
	if name is None:
		name = 'temp' #'5_5_D5_F'
	else:
		name = name

	return loc,nsteps,ndays,ff,date,tseas,ah,av,lon0,lat0,z0,zpar,do3d,doturb,name,grid

def test2():
	'''
	A drifter test using TXLA model output. 
	This simulation is 3D (do3d=1) with turbulence (doturb=1) added in.
	Drifters are started at 10 meters below the mean sea level and run backward (ff=-1)
	for five days from 11/25/09. Compare results with figure in examples/test2.png.
	'''

	# Location of TXLA model output
	# file and then grid
	loc = ['http://barataria.tamu.edu:8080/thredds/dodsC/txla_nesting6/ocean_his_0150.nc', \
			'http://barataria.tamu.edu:8080//thredds/dodsC/txla_nesting6_grid/txla_grd_v4_new.nc']

	# Initialize parameters
	nsteps = 10
	ndays = 5
	ff = -1
	# Start date
	date = datetime(2009,11, 25, 0)
	# Time between outputs
	# Dt = 14400. # in seconds (4 hours), nc.variables['dt'][:] 
	tseas = 4*3600 # 4 hours between outputs, in seconds, time between model outputs 
	ah = 100.
	av = 1.e-5 # m^2/s, or try 5e-6

	## Input starting locations as real space lon,lat locations
	lon0,lat0 = np.meshgrid(np.linspace(-94,-93,5), 
							np.linspace(28,29,5))
	lon0 = lon0.flatten()
	lat0 = lat0.flatten()

	## Choose method for vertical placement of drifters
	# # Also update makefile accordingly. Choose the twodim flag for isoslice.
	# # See above for more notes, but do the following two lines for an isoslice
	# z0 = 'z' #'salt' #'s' 
	# zpar = -10 #grid['km']-1 # 30 #grid['km']-1
	# Do the following two for a 3d simulation
	z0 = np.ones(lon0.shape)*-10 #  below the surface
	zpar = 'fromMSL' 

	# for 3d flag, do3d=0 makes the run 2d and do3d=1 makes the run 3d
	do3d = 1
	# turbulence/diffusion flag. doturb=0 means no turb/diffusion,
	# doturb=1 means adding parameterized turbulence
	# doturb=2 means adding diffusion on a circle
	# doturb=3 means adding diffusion on an ellipse (anisodiffusion)
	doturb = 1

	# simulation name, used for saving results into netcdf file
	name = 'test2'

	return loc,nsteps,ndays,ff,date,tseas,ah,av,lon0,lat0,z0,zpar,do3d,doturb,name

def hab1b():
	'''
	Initialize a drifter run using the starting locations from 
	HAB experiment 1b.
	'''

	if 'rainier' in os.uname():
		loc = '/Users/kthyng/Documents/research/postdoc/' # for model outputs
	elif 'hafen.tamu.edu' in os.uname():
		loc = '/home/kthyng/shelf/' # for model outputs

	# Initialize parameters
	nsteps = 10
	ndays = 10
	ff = 1
	# Start date
	date = datetime(2009,11, 30, 0)
	# Time between outputs
	# Dt = 14400. # in seconds (4 hours), nc.variables['dt'][:] 
	tseas = 4*3600 # 4 hours between outputs, in seconds, time between model outputs 
	ah = 100.
	av = 1.e-5 # m^2/s, or try 5e-6

	## Input starting locations as real space lon,lat locations
	# Read in starting locations from HAB experiment to test
	d = np.load(loc + 'hab/data/exp1b/starting_locations.npz')
	lon0 = d['lon0']
	lat0 = d['lat0']

	## Choose method for vertical placement of drifters
	# Also update makefile accordingly. Choose the twodim flag for isoslice.
	# See above for more notes, but do the following two lines for an isoslice
	z0 = 's' #'salt' #'s' 
	zpar = 29 #grid['km']-1 # 30 #grid['km']-1
	# Do the following two for a 3d simulation
	# z0 = np.ones(xstart0.shape)*-40 #  below the surface
	# zpar = 'fromMSL' 


	# for 3d flag, do3d=0 makes the run 2d and do3d=1 makes the run 3d
	do3d = 0
	# turbulence/diffusion flag. doturb=0 means no turb/diffusion,
	# doturb=1 means adding parameterized turbulence
	# doturb=2 means adding diffusion on a circle
	# doturb=3 means adding diffusion on an ellipse (anisodiffusion)
	doturb = 0

	# simulation name, used for saving results into netcdf file
	name = 'hab1b'

	return loc,nsteps,ndays,ff,date,tseas,ah,av,lon0,lat0,z0,zpar,do3d,doturb,name

