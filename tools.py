"""
Tools for dealing with drifter stuff.
"""

import numpy as np
from matplotlib.mlab import *
import pdb

def interpolate(x,y,z=None,grid,order=None,mode=None,):
	"""


	Inputs:

	"""

	# Horizontal interpolation only
	if z is None:


	# 3D interpolation
	else:


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


def convert_indices(direction,x,y,i,j):
	'''
	Converts indices between Python and Fortran indexing, assuming that
	Python indexing begins at 0 and Fortran (for x and y) begins at 1.
	In Tracmass, the vertical indexing does begin at zero so this script
	does nothing to vertical indexing.

	Usage:
		For before a call to tracmass:
			xstart,ystart,ia,ja = convert_indices('py2f',xstart,ystart,ia,ja)
		For after a call to tracmass:
			xend,yend,iend,jend = convert_indices('f2py',xend,yend,iend,jend)
	'''

	if direction == 'py2f':
		x = x+1
		y = y+1
		i = i+1
		j = j+1
	elif direction == 'f2py':
		x = x-1
		y = y-1
		i = i-1
		j = j-1

	return x,y,i,j

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
	fmask = grid['tricllrho'].nn_interpolator(grid['mask'].flatten())
	mask0 = fmask(lon0,lat0) # mask for lon0/lat0 points
	ind1 = (mask0==1.) # indices select out where points are masked

	ind2 = ~np.isnan(lon0)*ind1
	lon0 = lon0[ind2].flatten()
	lat0 = lat0[ind2].flatten()

	return lon0,lat0