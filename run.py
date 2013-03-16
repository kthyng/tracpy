import numpy as np
import sys
import op
import tracmass
import netCDF4 as netCDF
from mpl_toolkits.basemap import Basemap
import pdb
from matplotlib import delaunay
from matplotlib.pyplot import *


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

# Basemap parameters. These match those in particleTracking_plot.py
llcrnrlon=-98; llcrnrlat=22.9; urcrnrlon=-87.5; urcrnrlat=30.5; projection='lcc'
lat_0=30; lon_0=-94; resolution='i'; area_thresh=0.
basemap = Basemap(llcrnrlon=llcrnrlon,
              llcrnrlat=llcrnrlat,
              urcrnrlon=urcrnrlon,
              urcrnrlat=urcrnrlat,
              projection=projection,
              lat_0=lat_0,
              lon_0=lon_0,
              resolution=resolution,
              area_thresh=area_thresh)

# def run():

ff = 1 # forward

# Loop through
it = 0
# Read in ROMS model output velocity fields
nc = netCDF.Dataset('/Users/kthyng/Documents/research/postdoc/ocean_his_0150.nc')
# nc = netCDF.Dataset('http://barataria.tamu.edu:8080/thredds/dodsC/txla_nesting6/ocean_his_0001.nc?u[0:1][0:1:29][1:1:189][0:1:669]')
# Read in at time index it and it+1 and all depth levels
# Also only read in actual cells not ghost cells. Have to use only actual cells
# in order to have all fully-formed cells (with all walls)
u = nc.variables['u'][it:it+2,:,1:-1,:] 
v = nc.variables['v'][it:it+2,:,:,1:-1]

# Read in grid parameters
# These include ghost cells so don't match u and v correctly
grid = netCDF.Dataset('/Users/kthyng/Dropbox/python/tamu/hab/grid.nc')
lonu = grid.variables['lon_u'][:]
latu = grid.variables['lat_u'][:]
xu, yu = basemap(lonu,latu)
lonv = grid.variables['lon_v'][:]
latv = grid.variables['lat_v'][:]
xv, yv = basemap(lonv,latv)
hfull = grid.variables['h'][:]
lonr = grid.variables['lon_rho'][:]
latr = grid.variables['lat_rho'][:]
xr, yr = basemap(lonr,latr)
lonpsi = grid.variables['lon_psi'][:]
latpsi = grid.variables['lat_psi'][:]
xpsi, ypsi = basemap(lonpsi,latpsi)
X, Y = np.meshgrid(np.arange(latpsi.shape[0]),np.arange(latpsi.shape[1]))
tri = delaunay.Triangulation(X.flatten(),Y.flatten())
# Angle on rho grid
theta = grid.variables['angle'][:]
# Interpolate theta to be on psi grid
theta = op.resize(op.resize(theta,0),1)

# WHY DO I HAVE TO CALCULATE THE FLUXES HERE WHEN THE BOX VOLUMES ARE ALSO CALCULATED IN THE CODE?
# Change velocity fields to fluxes. DO TWO TIME INDICES AT ONCE
t = nc.variables['ocean_time'][it:it+2]
# Some grid metrics
s = nc.variables['s_w'][:] # sigma coords, 31 layers
cs = nc.variables['Cs_w'][:] # stretching curve in sigma coords, 31 layers
### Set up z array on u grid for calculating u fluxes
# Want h only within domain in y direction (no ghost cells) and interpolated onto cell
# walls in x direction
h = op.resize(grid.variables['h'][1:-1,:],1)
# Get zeta. Same as with h (but there is a time dimension first)
zeta = op.resize(nc.variables['zeta'][it:it+2,1:-1,:],2) # [t,k,j,i]
# Create 3D grid information arrays	
# 4D zeta for creating full 3D grid in time # [t,k,j,i]
zeta4 = zeta.reshape(zeta.shape[0],1,zeta.shape[1],zeta.shape[2]) \
					.repeat(len(cs),axis=1)
cs4 = cs.reshape(1,len(cs),1,1) \
					.repeat(zeta.shape[0],axis=0) \
					.repeat(zeta.shape[1],axis=2) \
					.repeat(zeta.shape[2],axis=3)
h4 = h.reshape(1,1,h.shape[0],h.shape[1]) \
					.repeat(zeta.shape[0],axis=0) \
					.repeat(len(cs),axis=1)	
z4 = zeta4 + cs4*(zeta4+h4)
del(zeta4,cs4,h4)
dz4 = np.diff(z4,axis=1).reshape((u.shape)) # grid cell thickness in z coords
###
# Calculate uflux, dependent on time. Calculating for two time indices at once.
# Take y's on u grid and interpolate in y, then take difference in y direction
# WHAT IS THE BEST WAY TO DO THIS ON A CURVILINEAR GRID?
dyu = np.diff(op.resize(yu,0),axis=0) # [JMT-1,IMT]
# uflux is size [2,KM,JMT-1,IMT]
uflux = u*dyu*dz4
# pdb.set_trace()

# Do for vfluxes
### Set up z array on u grid for calculating v fluxes
# Want h only within domain in x direction (no ghost cells) and interpolated onto cell
# walls in y direction
# should be able to clean these up with array broadcasting
h = op.resize(grid.variables['h'][:,1:-1],0)
# Get zeta. Same as with h (but there is a time dimension first)
zeta = op.resize(nc.variables['zeta'][it:it+2,:,1:-1],1) # [t,k,j,i]
# Create 3D grid information arrays	
# 4D zeta for creating full 3D grid in time # [t,k,j,i]
zeta4 = zeta.reshape(zeta.shape[0],1,zeta.shape[1],zeta.shape[2]) \
					.repeat(len(cs),axis=1)
cs4 = cs.reshape(1,len(cs),1,1) \
					.repeat(zeta.shape[0],axis=0) \
					.repeat(zeta.shape[1],axis=2) \
					.repeat(zeta.shape[2],axis=3)
h4 = h.reshape(1,1,h.shape[0],h.shape[1]) \
					.repeat(zeta.shape[0],axis=0) \
					.repeat(len(cs),axis=1)	
z4 = zeta4 + cs4*(zeta4+h4)
del(zeta4,cs4,h4)
dz4 = np.diff(z4,axis=1).reshape((v.shape)) # grid cell thickness in z coords
###
# Calculate uflux, dependent on time. Calculating for two time indices at once.
# Take x's on v grid and interpolate in x, then take difference in x direction
dxv = np.diff(op.resize(xv,1),axis=1) # [JMT,IMT-1]
# vflux is size [2,KM,JMT,IMT-1]
vflux = v*dxv*dz4

nc.close()
grid.close()

# Initialize seed locations
ia = [132]#,525]
ja = [57]#,40]
ka = [1]#,1]
# Need to input x,y as relative to their grid box. Let's assume they are at 
# some position relative to a grid node
xstart = [132.1]#,525.2]
ystart = [57.2]#,40.3]
# xstart, ystart = basemap(lonr[ja,ia],latr[ja,ia]) # locations in x,y coordiantes
zstart = [1]#,1]

t0save = t[0] # time at start of file in seconds since 1970-01-01, add this on at the end since it is big
tseas = 4*3600 # 4 hours between outputs 
hs = op.resize(zeta,1) # sea surface height in meters
# dz = np.diff(z4,axis=1)[0,:,0,0] # just take one dz column for now

# Initialize parameters
nsteps = 10 # Number of steps to do between model outputs

# Loop through iterations between the two model outputs

# NEED TO CHECK/FIX THE GRID INPUTS LATER. THEY MAY ALL BE SHIFTED, ETC

# dxdy is the horizontal area of the cell walls, so it should be the size of the number of cells only
# I am calculating it by adding the area of the walls together for each cell
dxdy = ((dyu[:,1:]+dyu[:,0:-1])+(dxv[1:,:]+dxv[0:-1,:])) # [JMT-1,IMT-1]

# calculate dxyz here
dxyz = op.resize(dz4,2)*dxdy

# make arrays in same order as is expected in the fortran code
# ROMS expects [time x k x j x i] but tracmass is expecting [i x j x k x time]
uflux = uflux.T
vflux = vflux.T
# dxdy = dxdy.T
dxyz = dxyz.T
hs = hs.T
# pdb.set_trace()


# CHANGE ALL ARRAYS TO BE FORTRAN-DIRECTIONED INSTEAD OF PYTHON
uflux = uflux.copy(order='f') # [2,KM,JMT-1,IMT] (all of these are swapped now due to transposes above)
vflux = vflux.copy(order='f') #[2,KM,JMT,IMT-1]]
# dxdy = dxdy.copy(order='f') #[JMT,IMT]]
dxyz = dxyz.copy(order='f') #[2,KM,JMT-1,IMT-1]
hs = hs.copy(order='f') # [JMT-1,IMT-1]
# These are the size of the psi grid, at the corners of the grid boxes.
# This might need to change to the number of boxes in the future, which will
# change all of the minus one's. Also would change the default values in the step call.
IMT=u.shape[3] # 670
JMT=v.shape[2] # 190
# IMT=dxdy.shape[1] # 669
# JMT=dxdy.shape[0] # 189
KM=u.shape[1] # 30

# Call tracmass subroutine
# f2py_tracmass()
# xend,yend,zend are particle locations at next step
# some variables are not specifically because f2py is hiding them from me:
#  IMT, JMT, JM, ntractot
# Look at tracmass.step to see what it is doing and making optional at the end.
# Do this by importing tracmass and then tracmass.step?
# I am assuming here that the velocity field at two times are being input into tracmass
# such that the output is the position for the drifters at the time corresponding to the
# second velocity time. Each drifter may take some number of steps in between, but those
# are not saved.
# pdb.set_trace()

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

# Recreate Cartesian particle locations from their index-relative locations
# First interpolate grid angle to drifter locations
ftheta = tri.nn_interpolator(theta.flatten())
thetap = ftheta(xend,yend)
# 
xp = xpsi[yend.astype(int),xend.astype(int)]+(xend*np.cos(thetap) - yend*np.sin(thetap))
yp = ypsi[yend.astype(int),xend.astype(int)]+(xend*np.sin(thetap) + yend*np.cos(thetap))

# Plot tracks
basemap.drawcoastlines()
basemap.fillcontinents('0.8')
basemap.drawparallels(np.arange(18, 35), dashes=(1, 0), linewidth=0.15, labels=[1, 0, 0, 0])
basemap.drawmeridians(np.arange(-100, -80), dashes=(1, 0), linewidth=0.15, labels=[0, 0, 0, 1])
hold('on')
contour(xr, yr, hfull, np.hstack(([10,20],np.arange(50,500,50))), colors='lightgrey', linewidths=0.15)
plot(xp,yp,'.')
# Outline numerical domain
plot(xr[0,:],yr[0,:],'k:')
#ax1.plot(xr[-1,:],yr[-1,:],'k:')
plot(xr[:,0],yr[:,0],'k:')
plot(xr[:,-1],yr[:,-1],'k:')

# if __name__=='__main__':
# 	run()
