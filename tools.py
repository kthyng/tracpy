"""
Tools for dealing with drifter stuff.

Functions include:

* interpolate2d
* interpolate3d
* find_final
* convert_indices
* check_points
"""

import numpy as np
from matplotlib.mlab import *
import pdb
from scipy import ndimage
import time

def interpolate2d(x,y,grid,itype,xin=None,yin=None,order=1,mode='nearest',cval=0.):
	"""
	Horizontal interpolation to map between coordinate transformations.

	Inputs:
		x, y 	x, y
		grid 	grid as read in by inout.readgrid()
		itype 	'd_xy2ij' delaunay, from projected x, y to grid i, j
				'd_ij2xy' delaunay, from grid i, j to projected x, y
				'd_ll2ij' delaunay, from lon, lat to grid i, j
				'd_ij2ll' delaunay, from grid i, j to lon, lat
				'm_ij2xy' map_coordinates, from grid i, j to projected x, y
				  or if z, xin, and yin are also input, from grid i, j, k to 
				  projected x, y, z. Can use the 3d version of this for transforming
				  to lon/lat also if the xin/yin input are lon/lat arrays.
				'm_ij2ll' map_coordinates, from grid i, j to lon, lat
		xin 	3D array of x values that are mapped to the input x,y,z coordinates.
				This is only needed in the 3D mapping case. Normally, can just do this
				in 2D instead of 3D and get the same results.
		yin 	3D array of y values that are mapped to the input x,y,z coordinates.
				This is only needed in the 3D mapping case. Normally, can just do this
				in 2D instead of 3D and get the same results.
		order 	order of interpolation for map_coordinates. 1 for linear 
				and 3 for cubic. Default=1
		mode 	behavior for edge points. Default is 'nearest'.
				Notes on the map_coordinates function: 
				The "mode" kwarg here just controls how the boundaries are treated
				mode='nearest' is _not_ nearest neighbor interpolation, it just uses the
				value of the nearest cell if the point lies outside the grid.  The default is
				to treat the values outside the grid as zero, which can cause some edge
				effects if you're interpolating points near the edge.
				'constant', 'nearest', 'reflect' or 'wrap'
				The "order" kwarg controls the order of the splines used. The default is 
				cubic splines, order=3
		cval 	Constant value used in map_coordinates if mode='constant'


	Outputs:
		xi,yi 	Interpolated values
		dt 		Time required for interpolation
	"""

	tic = time.time()

	if itype == 'd_xy2ij':
		# Set up functions for interpolating 
		fx = grid['trir'].nn_interpolator(grid['X'].flatten())
		fy = grid['trir'].nn_interpolator(grid['Y'].flatten())
		# Need to shift indices to move from rho grid of interpolator to arakawa c grid
		xi = fx(x,y) - .5
		yi = fy(x,y) - .5

	elif itype == 'd_ij2xy':
		# Set up functions for interpolating 
		fx = grid['tri'].nn_interpolator(grid['xr'].flatten())
		fy = grid['tri'].nn_interpolator(grid['yr'].flatten())
		# Need to shift indices to move to rho grid of interpolator from arakawa c grid
		xi = fx(x+0.5, y+0.5)
		yi = fy(x+0.5, y+0.5)

	elif itype == 'd_ll2ij':
		# Set up functions for interpolating 
		fx = grid['trirllrho'].nn_interpolator(grid['X'].flatten())
		fy = grid['trirllrho'].nn_interpolator(grid['Y'].flatten())
		# Need to shift indices to move from rho grid of interpolator to arakawa c grid
		xi = fx(x,y) - .5
		yi = fy(x,y) - .5

	elif itype == 'd_ij2ll':
		# Set up functions for interpolating 
		fx = grid['tri'].nn_interpolator(grid['lonr'].flatten())
		fy = grid['tri'].nn_interpolator(grid['latr'].flatten())
		# Need to shift indices to move to rho grid of interpolator from arakawa c grid
		xi = fx(x+0.5, y+0.5)
		yi = fy(x+0.5, y+0.5)

	elif itype == 'm_ij2xy':
		xi = ndimage.map_coordinates(grid['xr'], np.array([x.flatten()+.5,\
										y.flatten()+.5]), \
										order=order,\
										mode=mode,\
										cval=cval).reshape(x.shape)
		yi = ndimage.map_coordinates(grid['yr'], np.array([x.flatten()+.5,\
										y.flatten()+.5]), \
										order=order,\
										mode=mode,\
										cval=cval).reshape(y.shape)

	elif itype == 'm_ij2ll':
		xi = ndimage.map_coordinates(grid['lonr'], np.array([x.flatten()+.5,\
										y.flatten()+.5]), \
										order=order, \
										mode=mode,\
										cval=cval).reshape(x.shape)
		yi = ndimage.map_coordinates(grid['latr'], np.array([x.flatten()+.5,\
										y.flatten()+.5]), \
										order=order, \
										mode=mode,\
										cval=cval).reshape(y.shape)

	# Need to retain nan's since are changed them to zeros here
	ind = np.isnan(x)
	xi[ind] = np.nan
	yi[ind] = np.nan

	dt = time.time() - tic

	return xi, yi, dt


def interpolate3d(x,y,z,zin,order=1,mode='nearest',cval=0.):
	"""
	3D interpolation for transforming from grid/index space to whatever space is
	input with zin.

	Inputs:
		x,y,z	x, y, z coordinates
		zin 	3D array of z values that are mapped to the input x,y,z coordinates.
		order 	order of interpolation for map_coordinates. 1 for linear 
				and 3 for cubic. Default=1
		mode 	behavior for edge points. Default is 'nearest'.
				Notes on the map_coordinates function: 
				The "mode" kwarg here just controls how the boundaries are treated
				mode='nearest' is _not_ nearest neighbor interpolation, it just uses the
				value of the nearest cell if the point lies outside the grid.  The default is
				to treat the values outside the grid as zero, which can cause some edge
				effects if you're interpolating points near the edge
				The "order" kwarg controls the order of the splines used. The default is 
				cubic splines, order=3


	Outputs:
		zi 	 	Interpolated values
		dt 		Time required for interpolation
	"""

	tic = time.time()

	zi = ndimage.map_coordinates(zin, np.array([x.flatten()+.5, \
								y.flatten()+.5, \
								z.flatten()]), \
								order=order, \
								mode=mode,cval=cval).reshape(z.shape)

	# Need to retain nan's since are changed them to zeros here
	ind = np.isnan(z)
	zi[ind] = np.nan

	dt = time.time() - tic

	# pdb.set_trace()
	return zi, dt

def find_final(xp,yp):
	"""
	Loop through drifters and find final location of drifters
	within the tracks arrays. This can be necessary because when
	drifters exit the numerical domain, they are nan'ed out.
	"""

	# pdb.set_trace()

	# Find final position for drifters (since they are nan'ed out after they hit the open boundary)
	# Make this a separate function later
	xpc = []
	ypc = []
	for idrift in xrange(xp.shape[1]-1):
		# pdb.set_trace()
		# print idrift
		# Find last non-nan and make sure it is in the desired month start time
		ind3 = ~np.isnan(xp[:,idrift])
		#pdb.set_trace()
		# only plot if last non-nan (back in time) is in 1 month period
		# in order to plot the tracks that "started" in the plotted month
		if np.sum(ind3) > 1: # don't want tracks that start on land
			# This is for if we care when the drifter stopped
			# if t[find(ind3)[-1]] >= datetime(year,startMonth,startDay,0) and \
			# 	t[find(ind3)[-1]] <= datetime(year,startMonth+1,startDay,0):
			# ind2 = ~np.isnan(xp[idrift,:])
			if np.sum(np.isnan(xp[:,idrift])) > 0 and np.sum(np.isnan(xp[:,idrift])) < xp.shape[0]: # if there is a nan
				# ax.plot(xp[idrift,find(ind2)[-1]].T,yp[idrift,find(ind2)[-1]].T,'o',color='orange',linewidth=.5,label='_nolegend_')
				xpc.append(xp[find(ind3)[-1],idrift])
				ypc.append(yp[find(ind3)[-1],idrift])
			else:
				# ax.plot(xp[idrift,-1].T,yp[idrift,-1].T,'o',color='orange',linewidth=.5,label='_nolegend_')
				xpc.append(xp[find(ind3)[-1],idrift])
				ypc.append(yp[find(ind3)[-1],idrift])

	return xpc,ypc


def convert_indices(direction,x,y):
	'''
	Converts indices between Python and Fortran indexing, assuming that
	Python indexing begins at 0 and Fortran (for x and y) begins at 1.
	In Tracmass, the vertical indexing does begin at zero so this script
	does nothing to vertical indexing.

	Usage:
		For before a call to tracmass:
			xstart,ystart = convert_indices('py2f',xstart,ystart)
		For after a call to tracmass:
			xend,yend = convert_indices('f2py',xend,yend)
	'''

	if direction == 'py2f':
		x = x+1
		y = y+1
	elif direction == 'f2py':
		x = x-1
		y = y-1

	return x, y

def check_points(lon0,lat0,grid):
	"""
	Eliminate starting locations for drifters that are outside numerical domain
	and that are masked out.

	Inputs:
		lon0,lat0 	Starting locations for drifters in lon/lat
		grid 		Grid made from readgrid.py

	Outputs:
		lon0,lat0 	Fixed lon0,lat0
	"""

	lonr = grid['lonr']
	latr = grid['latr']

	# If covering the whole domain, need to exclude points outside domain.
	# Use info just inside domain so points aren't right at the edge.
	xvert = np.hstack((np.flipud(lonr[1,:]),lonr[:,1],
		lonr[-2,:],np.flipud(lonr[:,-2])))
	yvert = np.hstack((np.flipud(latr[1,:]),latr[:,1],
		latr[-2,:],np.flipud(latr[:,-2])))
	verts = np.vstack((xvert,yvert))
	# Form path
	path = Path(verts.T)
	# Loop through particle tracks to eliminate drifters that start outside 
	# the domain
	# pdb.set_trace()
	for jd in range(lon0.shape[0]): # loop through drifters
		for it in range(lon0.shape[1]): 
			# if drifter is not inside path, nan out this and all 
			# subsequent points
			# print jd,it
			# if it > 20:
			# 	pdb.set_trace()
			if not path.contains_point(np.vstack((lon0[jd,it],lat0[jd,it]))):
				lon0[jd,it] = np.nan
				lat0[jd,it] = np.nan
				# break

	# pdb.set_trace()

	# Also nan out points that are masked
	fmask = grid['trirllrho'].nn_interpolator(grid['mask'].flatten())
	mask0 = fmask(lon0,lat0) # mask for lon0/lat0 points
	ind1 = (mask0==1.) # indices select out where points are masked

	ind2 = ~np.isnan(lon0)*ind1
	lon0 = lon0[ind2].flatten()
	lat0 = lat0[ind2].flatten()

	return lon0,lat0