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
from readfields import *
from setupROMSfiles import *
from readgrid import *
from mpl_toolkits.basemap import Basemap
import time
from matplotlib.mlab import *
from convert_indices import *

'''

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

'''

def init_test1():
	'''
	A drifter test using TXLA model output. 
	This simulation is 2D (set twodim flag) with no turbulence (no turb flags).
	Drifters are started at 10 meters below the sea level and run forward
	for two days from 11/30/09. Compare results with figure in examples/test1.png.

	'''

	# Location of TXLA model output
	# file and then grid
	loc = ['http://barataria.tamu.edu:8080/thredds/dodsC/txla_nesting6/ocean_his_0150.nc', \
			'http://barataria.tamu.edu:8080//thredds/dodsC/txla_nesting6_grid/txla_grd_v4_new.nc']

	# Initialize parameters
	nsteps = 10
	ndays = 2
	ff = 1
	# Start date
	date = datetime(2009,11, 20, 0)
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
	# Also update makefile accordingly. Choose the twodim flag for isoslice.
	# See above for more notes, but do the following two lines for an isoslice
	z0 = 's' #'z' #'salt' #'s' 
	zpar = 29 #-10 #grid['km']-1 # 30 #grid['km']-1
	# Do the following two for a 3d simulation
	# z0 = np.ones(xstart0.shape)*-40 #  below the surface
	# zpar = 'fromMSL' 
	# pdb.set_trace()
	return loc,nsteps,ndays,ff,date,tseas,ah,av,lon0,lat0,z0,zpar

def init_test2():
	'''
	A drifter test using TXLA model output. 
	This simulation is 3D (comment out twodim flag) with turbulence added in.
	(uncomment turb flag).
	Drifters are started at 10 meters below the mean sea level and run forward
	for two days from 11/30/09. Compare results with figure in examples/test2.png.
	'''

	# Location of TXLA model output
	# file and then grid
	loc = ['http://barataria.tamu.edu:8080/thredds/dodsC/txla_nesting6/ocean_his_0150.nc', \
			'http://barataria.tamu.edu:8080//thredds/dodsC/txla_nesting6_grid/txla_grd_v4_new.nc']

	# Initialize parameters
	nsteps = 10
	ndays = 2
	ff = 1
	# Start date
	date = datetime(2009,11, 20, 0)
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
	z0 = np.ones(lon0.shape)*-1 #  below the surface
	zpar = 'fromMSL' 

	return loc,nsteps,ndays,ff,date,tseas,ah,av,lon0,lat0,z0,zpar

def init_hab1b():
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
	ndays = 2
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

	return loc,nsteps,ndays,ff,date,tseas,ah,av,lon0,lat0,z0,zpar

def run(loc,nsteps,ndays,ff,date,tseas,ah,av,lon0,lat0,z0,zpar):
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

	     hs 			Two time steps of ssh [imt,jmt,2] (m)

	Vertical location of drifters:
	z0 can be an array of starting depths (negative to be
	below the surface). 
	Or, if a 2D simulation is desired, make sure the 2dflag
	has been chosen for the simulation in the makefile and:
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
	To do 3D but start at surface, use z0=zeros(ia.shape) and have
	 either zpar='fromMSL'
	choose fromMSL to have z0 starting depths be for that depth below the base 
	time-independent sea level (or mean sea level).
	choose 'fromZeta' to have z0 starting depths be for that depth below the
	time-dependent sea surface. Haven't quite finished the 'fromZeta' case.

	'''

	tic = time.time()
	# Units for time conversion with netCDF.num2date and .date2num
	units = 'seconds since 1970-01-01'


	# ##### Necessary User Input #####

	# if 'rainier' in os.uname():
	# 	loc = '/Users/kthyng/Documents/research/postdoc/' # for model outputs
	# elif 'hafen.tamu.edu' in os.uname():
	# 	loc = '/home/kthyng/shelf/' # for model outputs

	# # Initialize parameters
	# # do2d = 1 # if do2d=1, 
	# nsteps = 10 # Number of steps to do between model outputs (iter in tracmass)
	# ndays = .75 # number of days to track the particles
	# ff = 1 # 1 forward, -1 backward
	# # Start date
	# date = datetime(2009,11, 30, 0)
	# # Time between outputs
	# Dt = 14400. # in seconds (4 hours), nc.variables['dt'][:] 
	# tseas = 4*3600 # 4 hours between outputs, in seconds, time between model outputs 
	# ah = 100. # horizontal diffusion in m^2/s. See project values of 350, 100, 0, 2000. For -turb,-diffusion
	# av = 1.e-5 # m^2/s, or try 5e-6

	# ## Input starting locations as real space x,y locations
	# # Read in starting locations from HAB experiment to test
	# d = np.load(loc + 'hab/data/exp1b/starting_locations.npz')
	# lon0 = d['lon0']
	# lat0 = d['lat0']

	# ## Choose method for vertical placement of drifters
	# # Also update makefile accordingly. Choose the twodim flag for isoslice.
	# # See above for more notes, but do the following two lines for an isoslice
	# z0 = 'salt' #'s' #'salt' #'s' 
	# zpar = 35 #grid['km']-1 # 30 #grid['km']-1
	# # Do the following two for a 3d simulation
	# # z0 = np.ones(xstart0.shape)*-40 #  below the surface
	# # zpar = 'fromMSL' 


	# ##### End of necessary user input #####

	# Number of model outputs to use
	tout = np.int((ndays*(24*3600))/tseas)

	# Convert date to number
	date = netCDF.date2num(date,units)

	# Figure out what files will be used for this tracking
	nc,tinds = setupROMSfiles(loc,date,ff,tout)

	# Read in grid parameters into dictionary, grid
	grid = readgrid(loc,nc)

	# Change input lat/lon into x,y then grid space
	# The basemap is set up for the northwestern Gulf of Mexico currently.
	x0,y0 = grid['basemap'](lon0,lat0)
	# Need to input x,y as relative to their grid box. Let's assume they are at 
	# some position relative to a grid node
	# Interpolate to get starting positions in grid space
	# tric is on a python index grid from 0 to imt-1 and 0 to jmt-1
	fX = grid['tric'].nn_interpolator(grid['X'].flatten())
	fY = grid['tric'].nn_interpolator(grid['Y'].flatten())
	xstart0 = fX(x0,y0)
	ystart0 = fY(x0,y0)
	# Do z a little lower down

	# Initialize seed locations 
	ia = np.ceil(xstart0) #[253]#,525]
	ja = np.ceil(ystart0) #[57]#,40]

	dates = nc.variables['ocean_time'][:]	
	t0save = dates[tinds[0]] # time at start of drifter test from file in seconds since 1970-01-01, add this on at the end since it is big

	# Initialize drifter grid positions and indices
	xend = np.ones(((len(tinds))*nsteps,ia.size))*np.nan
	yend = np.ones(((len(tinds))*nsteps,ia.size))*np.nan
	zend = np.ones(((len(tinds))*nsteps,ia.size))*np.nan
	zp = np.ones(((len(tinds))*nsteps,ia.size))*np.nan
	iend = np.ones(((len(tinds))*nsteps,ia.size))*np.nan
	jend = np.ones(((len(tinds))*nsteps,ia.size))*np.nan
	kend = np.ones(((len(tinds))*nsteps,ia.size))*np.nan
	ttend = np.ones(((len(tinds))*nsteps,ia.size))*np.nan
	t = np.zeros(((len(tinds))*nsteps+1))
	flag = np.zeros((ia.size),dtype=np.int) # initialize all exit flags for in the domain

	# Initialize vertical stuff and fluxes
	# Read initial field in - to 'new' variable since will be moved
	# at the beginning of the time loop ahead
	if is_string_like(z0): # isoslice case
		ufnew,vfnew,dztnew,zrtnew,zwtnew = readfields(tinds[0],grid,nc,z0,zpar)
	else: # 3d case
		ufnew,vfnew,dztnew,zrtnew,zwtnew = readfields(tinds[0],grid,nc)

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

	else:	# 3d case
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
				# Then find the vertical relative position in the grid cell	by adding on the bit of grid cell
				zstart0[i] = ka[i] - abs(z0[i]-grid['zwt0'][ia[i],ja[i],ka[i]]) \
									/abs(grid['zwt0'][ia[i],ja[i],ka[i]-1]-grid['zwt0'][ia[i],ja[i],ka[i]])
		# elif zpar == 'fromZeta':
		# 	for i in xrange(ia.size):
		# 		pdb.set_trace()
		# 		ind = (zwtnew[ia[i],ja[i],:]<=z0[i])
		# 		ka[i] = find(ind)[-1] # find value that is just shallower than starting vertical position
		# 		if (z0[i] != zwtnew[ia[i],ja[i],ka[i]]) and (ka[i] != grid['km']): # check this
		# 			ka[i] = ka[i]+1
		# 		# Then find the vertical relative position in the grid cell	by adding on the bit of grid cell
		# 		zstart0[i] = ka[i] - abs(z0[i]-zwtnew[ia[i],ja[i],ka[i]]) \
		# 							/abs(zwtnew[ia[i],ja[i],ka[i]-1]-zwtnew[ia[i],ja[i],ka[i]])

	# Find initial cell depths to concatenate to beginning of drifter tracks later
	# pdb.set_trace()
	zsave = zwtnew[ia.astype(int),ja.astype(int),ka.astype(int)]

	j = 0 # index for number of saved steps for drifters
	# Loop through model outputs. tinds is in proper order for moving forward
	# or backward in time, I think.
	for tind in tinds:
		# pdb.set_trace()
		# Move previous new time step to old time step info
		ufold = ufnew
		vfold = vfnew
		dztold = dztnew
		zrtold = zrtnew
		zwtold = zwtnew

		# Read stuff in for next time loop
		if is_string_like(z0): # isoslice case
			ufnew,vfnew,dztnew,zrtnew,zwtnew = readfields(tind+1,grid,nc,z0,zpar)
		else: # 3d case
			ufnew,vfnew,dztnew,zrtnew,zwtnew = readfields(tind+1,grid,nc)

		print j
		#  flux fields at starting time for this step
		if j != 0:
			xstart = xend[j*nsteps-1,:]
			ystart = yend[j*nsteps-1,:]
			zstart = zend[j*nsteps-1,:]
			ia = iend[j*nsteps-1,:]
			ja = jend[j*nsteps-1,:]
			ka = kend[j*nsteps-1,:]
			# mask out drifters that have exited the domain
			xstart = np.ma.masked_where(flag[:]==1,xstart)
			ystart = np.ma.masked_where(flag[:]==1,ystart)
			zstart = np.ma.masked_where(flag[:]==1,zstart)
			ia = np.ma.masked_where(flag[:]==1,ia)
			ja = np.ma.masked_where(flag[:]==1,ja)
			ka = np.ma.masked_where(flag[:]==1,ka)
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
			uflux = np.concatenate((ufold.reshape(np.append(ufold.shape,1)), \
									ufnew.reshape(np.append(ufnew.shape,1))), \
									axis=ufold.ndim)
			vflux = np.concatenate((vfold.reshape(np.append(vfold.shape,1)), \
									vfnew.reshape(np.append(vfnew.shape,1))), \
									axis=vfold.ndim)
			dzt = np.concatenate((dztold.reshape(np.append(dztold.shape,1)), \
									dztnew.reshape(np.append(dztnew.shape,1))), \
									axis=dztold.ndim)
			# pdb.set_trace()

			# Change the horizontal indices from python to fortran indexing 
			# (vertical are zero-based in tracmass)
			xstart,ystart,ia,ja = convert_indices('py2f',xstart,ystart,ia,ja)

			# km that is sent to tracmass is determined from uflux (see tracmass?)
			# so it will be the correct value for whether we are doing the 3D
			# or isoslice case.
			xend[j*nsteps:j*nsteps+nsteps,ind],\
				yend[j*nsteps:j*nsteps+nsteps,ind],\
				zend[j*nsteps:j*nsteps+nsteps,ind], \
				iend[j*nsteps:j*nsteps+nsteps,ind],\
				jend[j*nsteps:j*nsteps+nsteps,ind],\
				kend[j*nsteps:j*nsteps+nsteps,ind],\
				flag[ind],\
				ttend[j*nsteps:j*nsteps+nsteps,ind] = \
					tracmass.step(np.ma.compressed(xstart),np.ma.compressed(ystart),
						np.ma.compressed(zstart),np.ma.compressed(ia),np.ma.compressed(ja),
						np.ma.compressed(ka),tseas,uflux,vflux,ff,grid['kmt'].astype(int),
						dzt,grid['dxdy'],grid['dxv'],grid['dyu'],grid['h'],nsteps,ah,av)#dz.data,dxdy)
			# Change the horizontal indices from python to fortran indexing
			xend[j*nsteps:j*nsteps+nsteps,ind], \
				yend[j*nsteps:j*nsteps+nsteps,ind], \
				iend[j*nsteps:j*nsteps+nsteps,ind], \
				jend[j*nsteps:j*nsteps+nsteps,ind] = convert_indices('f2py',xend[j*nsteps:j*nsteps+nsteps,ind], \
																			yend[j*nsteps:j*nsteps+nsteps,ind], \
																			iend[j*nsteps:j*nsteps+nsteps,ind], \
																			jend[j*nsteps:j*nsteps+nsteps,ind])

			# Calculate times for the output frequency
			t[j*nsteps+1:j*nsteps+nsteps+1] = t[j*nsteps] + np.linspace(tseas/nsteps,tseas,nsteps) # update time in seconds to match drifters
			# Calculate real z position
			r = np.linspace(1./nsteps,1,nsteps) # linear time interpolation constant that is used in tracmass
			# Interpolate the drifter vertical depth in cell using the depths from the initial and later
			# time step.
			zp[j*nsteps:j*nsteps+nsteps,ind] = ((1.-r)*zwtold[iend[j*nsteps:j*nsteps+nsteps,ind].astype(int), \
															jend[j*nsteps:j*nsteps+nsteps,ind].astype(int), \
															kend[j*nsteps:j*nsteps+nsteps,ind].astype(int)].T \
												+ r*zwtnew[iend[j*nsteps:j*nsteps+nsteps,ind].astype(int), \
															jend[j*nsteps:j*nsteps+nsteps,ind].astype(int), \
															kend[j*nsteps:j*nsteps+nsteps,ind].astype(int)].T).T

		j = j + 1

	nc.close()

	t = t + t0save # add back in base time in seconds

	# Add on to front location for first time step
	xg=np.concatenate((xstart0.reshape(1,xstart0.size),xend),axis=0)
	yg=np.concatenate((ystart0.reshape(1,ystart0.size),yend),axis=0)
	# Concatenate zp with initial real space positions
	zp=np.concatenate((zsave.reshape(1,zstart0.size),zp),axis=0)

	# Recreate Cartesian particle locations from their index-relative locations
	# just by interpolating. These are in tracmass ordering
	fxr = grid['tri'].nn_interpolator(grid['xr'].flatten())
	fyr = grid['tri'].nn_interpolator(grid['yr'].flatten())
	xp = fxr(xg,yg)
	yp = fyr(xg,yg)

	# Plot tracks
	figure()
	grid['basemap'].drawcoastlines()
	grid['basemap'].fillcontinents('0.8')
	grid['basemap'].drawparallels(np.arange(18, 35), dashes=(1, 0), linewidth=0.15, labels=[1, 0, 0, 0])
	grid['basemap'].drawmeridians(np.arange(-100, -80), dashes=(1, 0), linewidth=0.15, labels=[0, 0, 0, 1])
	hold('on')
	contour(grid['xr'], grid['yr'], grid['h'], np.hstack(([10,20],np.arange(50,500,50))), colors='lightgrey', linewidths=0.15)
	plot(xp,yp,'g.')
	# Outline numerical domain
	plot(grid['xr'][0,:],grid['yr'][0,:],'k:')
	#ax1.plot(grid['xr'][-1,:],grid['yr'][-1,:],'k:')
	plot(grid['xr'][:,0],grid['yr'][:,0],'k:')
	plot(grid['xr'][:,-1],grid['yr'][:,-1],'k:')
	show()
	toc = time.time()
	print "run time:",toc-tic

	return xp,yp,zp,t

def start_run():
	'''
	Choose what initialization from above and then run.
	'''

	# Choose which initialization to use
	# loc,nsteps,ndays,ff,date,tseas,ah,av,lon0,lat0,z0,zpar = init_test1()
	loc,nsteps,ndays,ff,date,tseas,ah,av,lon0,lat0,z0,zpar = init_test2()
	# loc,nsteps,ndays,ff,date,tseas,ah,av,lon0,lat0,z0,zpar = init_hab1b()

	# Run tracmass!
	xp,yp,zp,t = run(loc,nsteps,ndays,ff,date,tseas,ah,av,lon0,lat0,z0,zpar)

	return xp,yp,zp,t