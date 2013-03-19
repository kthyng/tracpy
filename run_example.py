import numpy as np
import sys
import op
import tracmass
import netCDF4 as netCDF
from mpl_toolkits.basemap import Basemap
import pdb
from matplotlib import delaunay
from matplotlib.pyplot import *

# Use this to run simple examples before using the real model output.

# To re-compile tracmass fortran code, type "make f2py", which will give 
# a file tracmass.so, which is the module we import above. Then in ipython, "run run.py"

# Grid notes:
# With ghost cells, grid sizing is (in ROMS ordering, [j,i]):
#  psi grid: [JMT,IMT] (this doesn't change)
#  rho grid: [JMT+1,IMT+1]
#  u grid: [JMT+1,IMT]
#  v grid: [JMT,IMT+1]
# Without ghost cells, which is what I need to use to have full cells, grid sizing is:
#  psi grid: [JMT,IMT]
#  rho grid: [JMT-1,IMT-1]
#  u grid: [JMT-1,IMT]
#  v grid: [JMT,IMT-1]

# Initialize parameters
nsteps = 10 # Number of steps to do between model outputs
ff = 1 # forward
ex = 9           

grid = netCDF.Dataset('/Users/kthyng/Dropbox/python/tamu/hab/grid.nc')
maskr = grid.variables['mask_rho'][1:-1,1:-1]
grid.close()

# Need to input x,y as relative to their grid box. Let's assume they are at 
# some position relative to a grid node
xstart0 = np.array([[190.2,525.2]])
ystart0 = np.array([[100,40.3]])
# xstart, ystart = basemap(lonr[ja,ia],latr[ja,ia]) # locations in x,y coordiantes
zstart0 = np.array([[1,1]])
# Initialize seed locations # make these, e.g., ia=ceil(xstart0) (this was realized in cross.f95)
ia = np.ceil(xstart0) #[253]#,525]
ja = np.ceil(ystart0) #[57]#,40]
ka = np.ceil(zstart0) #[1]#,1]

# Load in example data
nc = netCDF.Dataset('example' + str(ex) + '.nc')
xpsi = nc.variables['xpsi'][:]
xu = nc.variables['xu'][:]
xv = nc.variables['xv'][:]
xr = nc.variables['xr'][:]
ypsi = nc.variables['ypsi'][:]
yu = nc.variables['yu'][:]
yv = nc.variables['yv'][:]
yr = nc.variables['yr'][:]
dxyz = nc.variables['dxyz'][:]
uflux = nc.variables['U'][:]
vflux = nc.variables['V'][:]
t = nc.variables['t'][:]
nc.close()

# Times
t0save = t[0] # time at start of file in seconds since 1970-01-01, add this on at the end since it is big
tseas = t[1]-t[0] #4*3600 # should be 4 hours between outputs 

# Just one time step for now
uflux = uflux[:,:,:,0:2]
vflux = vflux[:,:,:,0:2]
dxyz = dxyz[:,:,:,0:2]

# Loop over model outputs
xend = np.ones((nsteps,ia.size))*np.nan
yend = np.ones((nsteps,ia.size))*np.nan
zend = np.ones((nsteps,ia.size))*np.nan
ufluxinterp = uflux.copy()
vfluxinterp = vflux.copy()
# hsinterp = hs.copy()
dxyzinterp = dxyz.copy()
t0 = np.zeros((nsteps+1))
flag = np.zeros((ia.size),dtype=np.int) # initialize all exit flags for in the domain
# Loop over iterations between model outputs
for i in xrange(nsteps):
	#  flux fields at starting time for this step
	if i != 0:
		ufluxinterp[:,:,:,0] = ufluxinterp[:,:,:,1]
		vfluxinterp[:,:,:,0] = vfluxinterp[:,:,:,1]
		# hsinterp[:,:,0] = hsinterp[:,:,1]
		dxyzinterp[:,:,:,0] = dxyzinterp[:,:,:,1]
		xstart = xend[i-1,:]
		ystart = yend[i-1,:]
		zstart = zend[i-1,:]
		ia = xend[i-1,:].astype(int)
		ja = yend[i-1,:].astype(int)
		ka = zend[i-1,:].astype(int)
		# mask out drifters that have exited the domain
		xstart = np.ma.masked_where(flag[:]==1,xstart)
		ystart = np.ma.masked_where(flag[:]==1,ystart)
		zstart = np.ma.masked_where(flag[:]==1,zstart)
		ia = np.ma.masked_where(flag[:]==1,ia)
		ja = np.ma.masked_where(flag[:]==1,ja)
		ka = np.ma.masked_where(flag[:]==1,ka)
		ind = (flag[:] == 0) # indices where the drifters are still inside the domain
	else: # first loop
		xstart = xstart0
		ystart = ystart0
		zstart = zstart0
		# TODO: Do a check to make sure all drifter starting locations are within domain
		ind = (flag[:] == 0) # indices where the drifters are inside the domain to start
	# pdb.set_trace()
	# Linearly interpolate uflux and vflux to step between model outputs, field at ending time for this step
	ufluxinterp[:,:,:,1] = uflux[:,:,:,0] + (uflux[:,:,:,1]-uflux[:,:,:,0])*((i+1.)/nsteps)
	vfluxinterp[:,:,:,1] = vflux[:,:,:,0] + (vflux[:,:,:,1]-vflux[:,:,:,0])*((i+1.)/nsteps)
	# hsinterp[:,:,1] = hs[:,:,0] + (hs[:,:,1]-hs[:,:,0])*((i+1.)/nsteps)
	dxyzinterp[:,:,:,1] = dxyz[:,:,:,0] + (dxyz[:,:,:,1]-dxyz[:,:,:,0])*((i+1.)/nsteps)
	# only send unmasked values to step
	if not np.ma.compressed(xstart).any(): # exit if all of the drifters have exited the domain
		break
	else:
		xend[i,ind],yend[i,ind],zend[i,ind],flag[ind] = tracmass.step(np.ma.compressed(xstart),np.ma.compressed(ystart),
													np.ma.compressed(zstart),0.,np.ma.compressed(ia),np.ma.compressed(ja),
													np.ma.compressed(ka),tseas/float(nsteps),ufluxinterp,
													vfluxinterp,ff,dxyzinterp)#dz.data,dxdy)
		t0[i+1] = t0[i] + tseas/float(nsteps) # update time in seconds to match drifters

t0 = t0 + t0save # add back in base time in seconds

# Switch arrays to ROMS ordering from tracmass ordering
xr = xr.transpose(1,0)
yr = yr.transpose(1,0)

xg=np.concatenate((xstart0,xend))
yg=np.concatenate((ystart0,yend))
zg=np.concatenate((zstart0,zend))

# Recreate Cartesian particle locations from their index-relative locations
# just by interpolating
if ex != 1 and ex != 2 and ex !=3 and ex != 4 and ex !=5: # not on idealized grid only
	# # Angle on rho grid
	# grid = netCDF.Dataset('/Users/kthyng/Dropbox/python/tamu/hab/grid.nc')
	# # theta = grid.variables['angle'][1:-1,1:-1] #skip ghost cells
	# pm = grid.variables['pm'][1:-1,1:-1] #skip ghost cells
	# pn = grid.variables['pn'][1:-1,1:-1] #skip ghost cells
	# First interpolate grid angle to drifter locations
	X, Y = np.meshgrid(np.arange(xr.shape[1]),np.arange(yr.shape[0])) # grid in index coordinates, without ghost cells
	tri = delaunay.Triangulation(X.flatten(),Y.flatten())
	fxr = tri.nn_interpolator(xr.flatten())
	fyr = tri.nn_interpolator(yr.flatten())
	# ftheta = tri.nn_interpolator(theta.flatten())
	# fpm = tri.nn_interpolator(pm.flatten())
	# fpn = tri.nn_interpolator(pn.flatten())
	# thetap = ftheta(yg,xg)
	# pmp = fpm(xr[yg.astype(int),xg.astype(int)],yr[yg.astype(int),xg.astype(int)])
	# pnp = fpn(xr[yg.astype(int),xg.astype(int)],yr[yg.astype(int),xg.astype(int)])
	# xp = xr[yg.astype(int),xg.astype(int)]+(1/pmp)*(xg-np.floor(xg))
	# yp = yr[yg.astype(int),xg.astype(int)]+(1/pnp)*(yg-np.floor(yg))
	# xp = xr[yg.astype(int),xg.astype(int)]+pmp*(xg-np.floor(xg))
	# yp = yr[yg.astype(int),xg.astype(int)]+pnp*(yg-np.floor(yg))
	xp = fxr(xg,yg)
	yp = fyr(xg,yg)

# Plot tracks
figure()
colors = ['r','orange','yellow','green','b','purple','r','orange','yellow','green','b','purple','r','orange','yellow','green','b','purple',
	'r','orange','yellow','green','b','purple','r','orange','yellow','green','b','purple','r','orange','yellow','green','b','purple',
	'r','orange','yellow','green','b','purple','r','orange','yellow','green','b','purple','r','orange','yellow','green','b','purple',
	'r','orange','yellow','green','b','purple','r','orange','yellow','green','b','purple','r','orange','yellow','green','b','purple',
	'r','orange','yellow','green','b','purple','r','orange','yellow','green','b','purple','r','orange','yellow','green','b','purple',
	'r','orange','yellow','green','b','purple','r','orange','yellow','green','b','purple','r','orange','yellow','green','b','purple',
	'r','orange','yellow','green','b','purple','r','orange','yellow','green','b','purple','r','orange','yellow','green','b','purple',
	'r','orange','yellow','green','b','purple','r','orange','yellow','green','b','purple','r','orange','yellow','green','b','purple',
	'r','orange','yellow','green','b','purple','r','orange','yellow','green','b','purple','r','orange','yellow','green','b','purple',
	'r','orange','yellow','green','b','purple','r','orange','yellow','green','b','purple','r','orange','yellow','green','b','purple']
xr = np.ma.masked_where(maskr==0,xr)
yr = np.ma.masked_where(maskr==0,yr)
plot(xr,yr,'grey',xr.T,yr.T,'grey',alpha=.5)
# xlabel('y')
# ylabel('x')
hold('on')
for i in xrange(len(xp)):
	# have to plot these backwards because I am trying to stick with tracmass array ordering
	plot(xp[i],yp[i],'.',color=colors[i])
title('Example ' + str(ex))
show()
