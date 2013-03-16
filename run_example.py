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
ex = 6

# Initialize seed locations
ia = [51]#,525]
ja = [57]#,40]
ka = [1]#,1]
# Need to input x,y as relative to their grid box. Let's assume they are at 
# some position relative to a grid node
xstart0 = np.array([[50.5]])#,525.2]
ystart0 = np.array([[56.2]])#,40.3]
# xstart, ystart = basemap(lonr[ja,ia],latr[ja,ia]) # locations in x,y coordiantes
zstart0 = np.array([[1]])#,1]

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
xend = np.ones((nsteps,len(ia)))*np.nan
yend = np.ones((nsteps,len(ia)))*np.nan
zend = np.ones((nsteps,len(ia)))*np.nan
ufluxinterp = uflux.copy()
vfluxinterp = vflux.copy()
# hsinterp = hs.copy()
dxyzinterp = dxyz.copy()
t0 = np.zeros((nsteps+1,len(ia)))
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
	else:
		xstart = xstart0
		ystart = ystart0
		zstart = zstart0
	# pdb.set_trace()
	# Linearly interpolate uflux and vflux to step between model outputs, field at ending time for this step
	ufluxinterp[:,:,:,1] = uflux[:,:,:,0] + (uflux[:,:,:,1]-uflux[:,:,:,0])*((i+1.)/nsteps)
	vfluxinterp[:,:,:,1] = vflux[:,:,:,0] + (vflux[:,:,:,1]-vflux[:,:,:,0])*((i+1.)/nsteps)
	# hsinterp[:,:,1] = hs[:,:,0] + (hs[:,:,1]-hs[:,:,0])*((i+1.)/nsteps)
	dxyzinterp[:,:,:,1] = dxyz[:,:,:,0] + (dxyz[:,:,:,1]-dxyz[:,:,:,0])*((i+1.)/nsteps)
	# pdb.set_trace()
	# Find drifter locations
	# pdb.set_trace()
	xend[i,:],yend[i,:],zend[i,:],flag = tracmass.step(xstart,ystart,zstart,0.,ia,ja,ka,tseas/float(nsteps),ufluxinterp,vfluxinterp,ff,dxyzinterp)#dz.data,dxdy)
	t0[i+1] = t0[i] + tseas/float(nsteps) # update time in seconds to match drifters
	# FOR MORE THAN ONE DRIFTER, CHANGE SO THAT ONLY ACTIVE DRIFTERS ARE SENT TO STEP FUNCTION IN THE FIRST PLACE
	if flag == 1: # right now flag=1 if drifter has exited domain
		break

t0 = t0 + t0save # add back in base time in seconds

# # Recreate Cartesian particle locations from their index-relative locations
# # First interpolate grid angle to drifter locations
# ftheta = tri.nn_interpolator(theta.flatten())
# thetap = ftheta(xend,yend)
# # 
# xp = xpsi[yend.astype(int),xend.astype(int)]+(xend*np.cos(thetap) - yend*np.sin(thetap))
# yp = ypsi[yend.astype(int),xend.astype(int)]+(xend*np.sin(thetap) + yend*np.cos(thetap))
xp=np.concatenate((xstart0,xend))
yp=np.concatenate((ystart0,yend))
zp=np.concatenate((zstart0,zend))

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
plot(xr,yr,'grey',xr.T,yr.T,'grey',alpha=.5)
xlabel('y')
ylabel('x')
hold('on')
for i in xrange(len(xp)):
	# have to plot these backwards because I am trying to stick with tracmass array ordering
	plot(yp[i],xp[i],'.',color=colors[i])
title('Example ' + str(ex))
show()
