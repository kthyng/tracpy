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

'''

To re-compile tracmass fortran code, type "make f2py", which will give 
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
'''
tic = time.time()
# Units for time conversion with netCDF.num2date and .date2num
units = 'seconds since 1970-01-01'

if 'rainier' in os.uname():
	loc = '/Users/kthyng/Documents/research/postdoc/' # for model outputs
elif 'hafen.tamu.edu' in os.uname():
	loc = '/home/kthyng/shelf/' # for model outputs

# Initialize parameters
nsteps = 10 # Number of steps to do between model outputs (iter in tracmass)
ndays = .25 # number of days to track the particles
ff = -1 # 1 forward, -1 backward
# Start date
date = datetime(2009,11, 30, 0)
# Convert date to number
date = netCDF.date2num(date,units)
# Time between outputs
Dt = 14400. # in seconds (4 hours), nc.variables['dt'][:] 
# Number of model outputs to use
tout = np.int((ndays*(24*3600))/Dt)
tseas = 4*3600 # 4 hours between outputs, in seconds, time between model outputs 
ah = 100. # horizontal diffusion in m^2/s. See project values of 350, 100, 0, 2000. For -turb,-diffusion
av = 1.e-5 # m^2/s, or try 5e-6

# THINK I NEED THIS STARTED HERE AND THEN UPDATED IN STEP
# tt = ints*tseas

# Figure out what files will be used for this tracking
nc,tinds = setupROMSfiles(loc,date,ff,tout)

# Read in grid parameters into dictionary, grid
grid = readgrid(loc,nc)



# Need to input x,y as relative to their grid box. Let's assume they are at 
# some position relative to a grid node
# xstart0 = np.array([[525.2]])#525.2]])
# ystart0 = np.array([[40.3]])
# Read in starting locations from HAB experiment to test
d = np.load(loc + 'hab/data/exp1b/starting_locations.npz')
# pdb.set_trace()
x0,y0 = grid['basemap'](d['lon0'],d['lat0'])
# Interpolate to get starting positions in grid space
# tric is on a python index grid from 0 to imt-1 and 0 to jmt-1
fX = grid['tric'].nn_interpolator(grid['X'].flatten())
fY = grid['tric'].nn_interpolator(grid['Y'].flatten())
xstart0 = fX(x0,y0)
ystart0 = fY(x0,y0)
# Do z a little lower down

# Initialize seed locations # make these, e.g., ia=ceil(xstart0) (this was realized in cross.f95)
# # want ceil(xstart0) for drifter positions within the cell, but if xstart0 is on a grid cell border,
# # want ceil(xstart0+1)
# ia = np.ceil(xstart0)+(1-(np.ceil(xstart0)-np.floor(xstart0)))
# ja = np.ceil(ystart0)+(1-(np.ceil(ystart0)-np.floor(ystart0)))
# ka = np.ceil(zstart0)+(1-(np.ceil(zstart0)-np.floor(zstart0)))
ia = np.ceil(xstart0) #[253]#,525]
ja = np.ceil(ystart0) #[57]#,40]

# pdb.set_trace()
dates = nc.variables['ocean_time'][:]	
t0save = dates[0] # time at start of file in seconds since 1970-01-01, add this on at the end since it is big

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

# Initialize free surface and fluxes
uflux = np.ones((grid['imt']-1,grid['jmt'],grid['km'],2))*np.nan # uflux, 2 times
vflux = np.ones((grid['imt'],grid['jmt']-1,grid['km'],2))*np.nan # vflux, 2 times
dzt = np.ones((grid['imt'],grid['jmt'],grid['km'],2))*np.nan # 2 times
# Read initial field in - to 2nd time spot since will be moved
# at the beginning of the time loop ahead
uflux[:,:,:,1],vflux[:,:,:,1],dzt[:,:,:,1] = readfields(tinds[0],grid,nc)


# Input initial vertical locations in real space. 
# These are converted to grid space once dzt is calculated
# To place drifters on the surface, they should actually be just below the 
# surface because the vertical grid does not go to zero at the surface (just gets near it)
# For surface: Note that if surface is chosen here and the twodim flag is chosen in the 
# makefile, the drifters will stay at their relative position in the uppermost grid cells,
# that is, in the middle (vertically). Their actual, real space vertical position will 
# change in time.
z0 = 'surface'
# z0 = np.ones(xstart0.shape)*-.05 #  below the surface

# zstart0 = np.ones(xstart0.shape)*grid['km'] # The surface in tracmass is at k=KM

# For special case of the surface, take the surface grid cell index and the depth at 
# the middle of the surface layer
if z0 == 'surface':
	# subtract 1 from km since in python indexing here
	ka = np.ones(ia.size)*grid['km']-1
	# Find the corresponding depths at the rho grid in the topmost level
	zt = depths.get_zrho(1, 1, grid['km'], grid['theta_s'], grid['theta_b'], 
						grid['h'].T.copy(order='c'), 
						grid['hc'], zeta=hs[:,:,1].T.copy(order='c'), Hscale=3)
	zt = zt.T.copy(order='f')
	# change z0 to be what it is in the normal case: the real space depth for the drifters
	z0 = zt[ia.astype(int),ja.astype(int),ka.astype(int)]
	zstart0 = np.ones(ia.size)*(ka-.5) # start mid-way vertically through the uppermost grid cell
else:	
	# Convert initial real space vertical locations to grid space
	# first find indices of grid cells vertically
	ka = np.ones(ia.size)*np.nan
	zstart0 = np.ones(ia.size)*np.nan
	for i in xrange(ia.size):
		ind = (grid['dzt0'][ia[i],ja[i],:]<=z0[i])
		ka[i] = find(ind)[-1] # find value that is just shallower than starting vertical position
		if (z0[i] != grid['dzt0'][ia[i],ja[i],ka[i]]) and (ka[i] != grid['km']-1):
			ka[i] = ka[i]+1
		# Then find the vertical relative position in the grid cell	by adding on the bit of grid cell
		zstart0[i] = ka[i] - abs(z0[i]-grid['dzt0'][ia[i],ja[i],ka[i]])/abs(grid['dzt0'][ia[i],ja[i],ka[i]-1]-grid['dzt0'][ia[i],ja[i],ka[i]])

# Bump all grid-based fields up by one since in Fortran they are 1-based instead of 0-based
ia = ia + 1
ja = ja + 1
ka = ka + 1
xstart0 = xstart0 + 1
ystart0 = ystart0 + 1
zstart0 = zstart0 + 1

j = 0 # index for number of saved steps for drifters
# Loop through model outputs. tinds is in proper order for moving forward
# or backward in time, I think.
for tind in tinds:
	# pdb.set_trace()
	# Move previous new time step to old time step info
	uflux[:,:,:,0] = uflux[:,:,:,1]
	vflux[:,:,:,0] = vflux[:,:,:,1]
	dzt[:,:,:,0] = dzt[:,:,:,1]

	# Read stuff in for next time loop
	uflux[:,:,:,1],vflux[:,:,:,1],dzt[:,:,:,1] = readfields(tind+1,grid,nc)
	# pdb.set_trace()

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
		# t0[0] = tseas/nsteps # want to start time at the first moved step, will add on initial zero at the end

	# Find drifter locations
	# pdb.set_trace()
	# only send unmasked values to step
	if not np.ma.compressed(xstart).any(): # exit if all of the drifters have exited the domain
		break
	else:
		xend[j*nsteps:j*nsteps+nsteps,ind],\
			yend[j*nsteps:j*nsteps+nsteps,ind],\
			zend[j*nsteps:j*nsteps+nsteps,ind], \
			zp[j*nsteps:j*nsteps+nsteps,ind],\
			iend[j*nsteps:j*nsteps+nsteps,ind],\
			jend[j*nsteps:j*nsteps+nsteps,ind],\
			kend[j*nsteps:j*nsteps+nsteps,ind],\
			flag[ind],\
			ttend[j*nsteps:j*nsteps+nsteps,ind] = \
				tracmass.step(np.ma.compressed(xstart),np.ma.compressed(ystart),
					np.ma.compressed(zstart),np.ma.compressed(ia),np.ma.compressed(ja),
					np.ma.compressed(ka),tseas,uflux,vflux,ff,grid['kmt'].astype(int),
					dzt,grid['dxdy'],grid['dxv'],grid['dyu'],grid['h'],nsteps,ah,av)#dz.data,dxdy)
		# if np.sum(flag)>0:
		# pdb.set_trace()
		t[j*nsteps+1:j*nsteps+nsteps+1] = t[j*nsteps] + np.linspace(tseas/nsteps,tseas,nsteps)    #tseas/float(nsteps) # update time in seconds to match drifters

	j = j + 1

nc.close()
# grid.close()

# # Add on time for first time step
# t0=np.concatenate((0,t0),axis=0)

t = t + t0save # add back in base time in seconds

# Add on to front location for first time step
xg=np.concatenate((xstart0.reshape(1,xstart0.size),xend),axis=0)
yg=np.concatenate((ystart0.reshape(1,ystart0.size),yend),axis=0)
# Concatenate zp with initial real space positions
zp=np.concatenate((z0.reshape(1,zstart0.size),zp),axis=0)

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
# if __name__=='__main__':
# 	run()
