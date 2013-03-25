import numpy as np
import sys
import op
import tracmass
import netCDF4 as netCDF
from mpl_toolkits.basemap import Basemap
import pdb
from matplotlib import delaunay
from matplotlib.pyplot import *
import glob
from datetime import datetime, timedelta


# To re-compile tracmass fortran code, type "make f2py", which will give 
# a file tracmass.so, which is the module we import above. Then in ipython, "run run.py"
# xend,yend,zend are particle locations at next step
# some variables are not specifically because f2py is hiding them from me:
#  IMT, JMT, JM, ntractot
# Look at tracmass.step to see what it is doing and making optional at the end.
# Do this by importing tracmass and then tracmass.step?

# I am assuming here that the velocity field at two times are being input into tracmass
# such that the output is the position for the drifters at the time corresponding to the
# second velocity time. Each drifter may take some number of steps in between, but those
# are not saved.

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
# Units for time conversion with netCDF.num2date and .date2num
units = 'seconds since 1970-01-01'

loc = '/home/kthyng/shelf/' # for model outputs
# loc = '/Users/kthyng/Documents/research/postdoc/' # for model outputs

# Read in grid parameters
# These include ghost cells so don't match u and v correctly
grid = netCDF.Dataset(loc + 'grid.nc')
lonu = grid.variables['lon_u'][:]
latu = grid.variables['lat_u'][:]
xu, yu = basemap(lonu,latu)
lonv = grid.variables['lon_v'][:]
latv = grid.variables['lat_v'][:]
xv, yv = basemap(lonv,latv)
hfull = grid.variables['h'][:]
lonr = grid.variables['lon_rho'][:]#[1:-1,1:-1]
latr = grid.variables['lat_rho'][:]#[1:-1,1:-1]
xr, yr = basemap(lonr,latr)
lonpsi = grid.variables['lon_psi'][:]
latpsi = grid.variables['lat_psi'][:]
xpsi, ypsi = basemap(lonpsi,latpsi)
maskr = grid.variables['mask_rho'][:]#[1:-1,1:-1]
X, Y = np.meshgrid(np.arange(xr.shape[1]),np.arange(yr.shape[0])) # grid in index coordinates, without ghost cells
# Triangulation for grid space to curvilinear space
tri = delaunay.Triangulation(X.flatten(),Y.flatten())
# Triangulation for curvilinear space to grid space
tric = delaunay.Triangulation(xr.flatten(),yr.flatten())
# Angle on rho grid
theta = grid.variables['angle'][:]
# Interpolate theta to be on psi grid
theta = op.resize(op.resize(theta,0),1)
pm = grid.variables['pm'][:] # 1/dx
pn = grid.variables['pn'][:] # 1/dy
# The following are for setting up z array on u (v) grids for calculating u (v) fluxes
# Want h only within domain in y (x) direction (no ghost cells) and interpolated onto cell
# walls in x (y) direction
h = grid.variables['h'][:]
hu = op.resize(grid.variables['h'][1:-1,:],1)
hv = op.resize(grid.variables['h'][:,1:-1],0)

# Initialize parameters
nsteps = 10 # Number of steps to do between model outputs
ndays = .25 # number of days to track the particles
ff = 1 # forward
# Start date
date = datetime(2009, 9, 1, 0)
# Convert date to number
date = netCDF.date2num(date,units)
# Time between outputs
Dt = 14400. # in seconds (4 hours), nc.variables['dt'][:] 
# Number of model outputs to use
tout = np.int((ndays*(24*3600))/Dt)
tseas = 4*3600 # 4 hours between outputs 

# Need to input x,y as relative to their grid box. Let's assume they are at 
# some position relative to a grid node
# xstart0 = np.array([[525.2]])#525.2]])
# ystart0 = np.array([[40.3]])
# Read in starting locations from HAB experiment to test
d = np.load('/home/kthyng/hab/data/exp8/starting_locations.npz')
x0,y0 = basemap(d['lon0'],d['lat0'])
# Interpolate to get starting positions in grid space
fX = tric.nn_interpolator(X.flatten())
fY = tric.nn_interpolator(Y.flatten())
xstart0 = fX(x0,y0)
ystart0 = fY(x0,y0)
zstart0 = np.ones(xstart0.shape)*(29)
# Initialize seed locations # make these, e.g., ia=ceil(xstart0) (this was realized in cross.f95)
# # want ceil(xstart0) for drifter positions within the cell, but if xstart0 is on a grid cell border,
# # want ceil(xstart0+1)
# ia = np.ceil(xstart0)+(1-(np.ceil(xstart0)-np.floor(xstart0)))
# ja = np.ceil(ystart0)+(1-(np.ceil(ystart0)-np.floor(ystart0)))
# ka = np.ceil(zstart0)+(1-(np.ceil(zstart0)-np.floor(zstart0)))
ia = np.ceil(xstart0) #[253]#,525]
ja = np.ceil(ystart0) #[57]#,40]
ka = np.ceil(zstart0) #[1]#,1]

files = np.sort(glob.glob(loc + 'ocean_his_*.nc')) # sorted list of file names
filesfull = np.sort(glob.glob(loc + 'ocean_his_*.nc')) #full path of files
# Find the list of files that cover the desired time period
for i,name in enumerate(files): # Loop through files
	nctemp = netCDF.Dataset(name)
	ttemp = nctemp.variables['ocean_time'][:]
	nctemp.close()
	# If datenum_in is larger than the first time in the file but smaller
	# than the last time, then this is the correct file to use to start
	if date > ttemp[0] and date <= ttemp[-1]:
		ifile = i # this is the starting file identifier then
		break
# Since the number of indices per file can change, make the process
# of finding the necessary files a little more general
# Start by opening two files
i = 1
fname = [files[ifile]]


# pdb.set_trace()

# if ff: #forward - add 2nd file on end
# 	fname.append(files[ifile+i])
# else: #backward - add previous time file to beginning
# 	fname.insert(0,files[ifile-i])
nc = netCDF.MFDataset(fname) # files in fname are in chronological order
# number of indices included in opened files so far
ninds = nc.variables['ocean_time'][:].size 
# Find which output in ifile is closest to the user-input start time (choose lower index)
# Dates for drifters from this start date
dates = nc.variables['ocean_time'][:]	
ilow = date >= dates
# time index with time value just below date (relative to file ifile)
istart = dates[ilow].size - 1
nc.close()
# Select indices 
if ff:
	tinds = range(istart,istart+tout) # indices of model outputs desired
else: # backward in time
	# have to shift istart since there are now new indices behind since going backward
	tinds = range(istart,istart-tout,-1)
# If we need more indices than available in these files, add another

if ff:
	# if the final index we want is beyond the length of these files,
	# keep adding files on
	while tinds[-1] >= len(dates): 
		# if tdir: #forward - add 2nd file on end
		fname.append(files[ifile+i])
		nc = netCDF.MFDataset(fname) # files in fname are in chronological order
		dates = nc.variables['ocean_time'][:]	
		ilow = date >= dates
		# time index with time value just below datenum_in (relative to file ifile)
		istart = dates[ilow].size - 1
		tinds = range(istart,istart+tout)
		nc.close()
		i = i + 1
else: #backwards in time
	while tinds[-1] < 0:
		fname.insert(0,files[ifile-i])
		nc = netCDF.MFDataset(fname)
		dates = nc.variables['ocean_time'][:]	
		ilow = date >= dates
		# time index with time value just below datenum_in (relative to file ifile)
		istart = dates[ilow].size - 1
		tinds = range(istart,istart-tout,-1)
		nc.close()
		i = i + 1


# i = 2
# if ff:
# 	# if the final index we want is beyond the length of these files,
# 	# keep adding files on
# 	while tinds[-1] >= len(dates): 
# 		# if tdir: #forward - add 2nd file on end
# 		fname.append(files[ifile+i])
# 		nc = netCDF.MFDataset(fname) # files in fname are in chronological order
# 		dates = nc.variables['ocean_time'][:]	
# 		ilow = date >= dates
# 		# time index with time value just below datenum_in (relative to file ifile)
# 		istart = dates[ilow].size - 1
# 		tinds = range(istart,istart+tout)
# 		nc.close()
# 		i = i + 1
# else: #backwards in time
# 	while tinds[-1] < 0:
# 		fname.insert(0,files[ifile-i])
# 		nc = netCDF.MFDataset(fname)
# 		dates = nc.variables['ocean_time'][:]	
# 		ilow = date >= dates
# 		# time index with time value just below datenum_in (relative to file ifile)
# 		istart = dates[ilow].size - 1
# 		tinds = range(istart,istart-tout,-1)
# 		nc.close()
# 		i = i + 1
# model output files together containing all necessary model outputs
nc = netCDF.MFDataset(fname) # reopen since needed to close things in loop
# pdb.set_trace()
dates = nc.variables['ocean_time'][:]	
t0save = dates[0] # time at start of file in seconds since 1970-01-01, add this on at the end since it is big
# Some grid metrics
s = nc.variables['s_w'][:] # sigma coords, 31 layers
cs = nc.variables['Cs_w'][:] # stretching curve in sigma coords, 31 layers
# Basing this on setupgrid.f95 for rutgersNWA example project from Bror
IMT = h.shape[1] # 671
JMT = h.shape[0] # 191
KM = s.shape[0]-1 # 30
# ROMS ordering. dxv,dyu: [191,674]=[JMT,IMT+3]
dxv = np.ones((JMT,IMT+3))*np.nan
dxv[:,:-3] = xr
dxv[:,0:IMT-2] = dxv[:,1:IMT-1] - dxv[:,0:IMT-2]
dxv[:,IMT-1:IMT] = dxv[:,IMT-3:IMT-2]
dyu = np.ones((JMT,IMT+3))*np.nan
dyu[:,:-3] = yr
dyu[0:JMT-2,:] = dyu[1:JMT-1,:] - dyu[0:JMT-2,:]
dyu[JMT-1,:] = dyu[JMT-2,:]
# Adjust masking according to setupgrid.f95 for rutgersNWA example project from Bror
pdb.set_trace()
# my best guess is that kmt is an array containing the number of vertical
# grid cells in every (i,j) location. For sigma coords it doesn't change
# but in other codes it could. Also contains masking information.
kmt = np.ones((JMT,IMT))*KM
ind = (maskr[:,1:IMT-1]==1)
maskr[ind,0:IMT-2] = 1
ind = (maskr[1:IMT-1,:]==1)
maskr[0:JMT-2,:] = 1
ind = (maskr==0)
kmt[ind] = 0

# Initialize drifter grid positions and indices
xend = np.ones(((len(tinds))*nsteps,ia.size))*np.nan
yend = np.ones(((len(tinds))*nsteps,ia.size))*np.nan
zend = np.ones(((len(tinds))*nsteps,ia.size))*np.nan
iend = np.ones(((len(tinds))*nsteps,ia.size))*np.nan
jend = np.ones(((len(tinds))*nsteps,ia.size))*np.nan
kend = np.ones(((len(tinds))*nsteps,ia.size))*np.nan
t0 = np.zeros(((len(tinds))*nsteps+1))
flag = np.zeros((ia.size),dtype=np.int) # initialize all exit flags for in the domain

j = 0 # index for number of saved steps for drifters
# Loop through model outputs
for tind in tinds:
	# Read in ROMS model output velocity fields at time index tind one index
	# in front or behind depending on time direction
	# Only read in actual cells not ghost cells. Have to use only actual cells
	# in order to have all fully-formed cells (with all walls)
	if ff: # forward
		u = nc.variables['u'][tind:tind+2,:,1:-1,:] 
		v = nc.variables['v'][tind:tind+2,:,:,1:-1]
		zetau = op.resize(nc.variables['zeta'][tind:tind+2,1:-1,:],2) # [t,j,i]
		zetav = op.resize(nc.variables['zeta'][tind:tind+2,:,1:-1],1) # [t,j,i]
	else:
		u = nc.variables['u'][tind-1:tind-3,:,1:-1,:] 
		v = nc.variables['v'][tind-1:tind-3,:,:,1:-1]
		zetau = op.resize(nc.variables['zeta'][tind-1:tind-3,1:-1,:],2) # [t,j,i]
		zetav = op.resize(nc.variables['zeta'][tind-1:tind-3,:,1:-1],1) # [t,j,i]


	# WHY DO I HAVE TO CALCULATE THE FLUXES HERE WHEN THE BOX VOLUMES ARE ALSO CALCULATED IN THE CODE?
	# Creating full vertical grid in time # [t,k,j,i]
	zeta4u = zetau.reshape(zetau.shape[0],1,zetau.shape[1],zetau.shape[2]) \
						.repeat(len(cs),axis=1)
	cs4u = cs.reshape(1,len(cs),1,1) \
						.repeat(zetau.shape[0],axis=0) \
						.repeat(zetau.shape[1],axis=2) \
						.repeat(zetau.shape[2],axis=3)
	h4u = hu.reshape(1,1,hu.shape[0],hu.shape[1]) \
						.repeat(zetau.shape[0],axis=0) \
						.repeat(len(cs),axis=1)	
	z4u = zeta4u + cs4u*(zeta4u+h4u)
	del(zeta4u,cs4u,h4u)
	dz4u = np.diff(z4u,axis=1).reshape((u.shape)) # grid cell thickness in z coords
	# Change velocity fields to fluxes.
	# Calculate uflux, dependent on time. Calculating for two time indices at once.
	# # Use pn for dy
	# dyu = 1/op.resize(pn[1:-1,:],1) # [JMT-1,IMT]
	# uflux is size [2,KM,JMT-1,IMT]
	uflux = u*dyu*dz4u

	# Do for vfluxes
	zeta4v = zetav.reshape(zetav.shape[0],1,zetav.shape[1],zetav.shape[2]) \
						.repeat(len(cs),axis=1)
	cs4v = cs.reshape(1,len(cs),1,1) \
						.repeat(zetav.shape[0],axis=0) \
						.repeat(zetav.shape[1],axis=2) \
						.repeat(zetav.shape[2],axis=3)
	h4v = hv.reshape(1,1,hv.shape[0],hv.shape[1]) \
						.repeat(zetav.shape[0],axis=0) \
						.repeat(len(cs),axis=1)	
	z4v = zeta4v + cs4v*(zeta4v+h4v)
	del(zeta4v,cs4v,h4v)
	dz4v = np.diff(z4v,axis=1).reshape((v.shape)) # grid cell thickness in z coords
	# Calculate vflux, dependent on time. Calculating for two time indices at once.
	# # Use pm for dx
	# dxv = 1/op.resize(pm[:,1:-1],0) # [JMT,IMT-1]
	# vflux is size [2,KM,JMT,IMT-1]
	vflux = v*dxv*dz4v

	# dxdy is the horizontal area of the cell walls, so it should be the size of the number of cells only
	# I am calculating it by adding the area of the walls together for each cell
	dxdy = ((dyu[:,1:]+dyu[:,0:-1])+(dxv[1:,:]+dxv[0:-1,:])) # [JMT-1,IMT-1]

	# calculate dxyz here
	dxyz = op.resize(dz4v,2)*dxdy

	# hs = op.resize(zeta,1) # sea surface height in meters
	# dz = np.diff(z4,axis=1)[0,:,0,0] # just take one dz column for now

	# make arrays in same order as is expected in the fortran code
	# ROMS expects [time x k x j x i] but tracmass is expecting [i x j x k x time]
	uflux = uflux.T
	vflux = vflux.T
	# dxdy = dxdy.T
	dxyz = dxyz.T
	# hs = hs.T

	# change all arrays to be fortran-directioned instead of python
	uflux = uflux.copy(order='f') # [2,KM,JMT-1,IMT] (all of these are swapped now due to transposes above)
	vflux = vflux.copy(order='f') #[2,KM,JMT,IMT-1]]
	# dxdy = dxdy.copy(order='f') #[JMT,IMT]]
	dxyz = dxyz.copy(order='f') #[2,KM,JMT-1,IMT-1]
	# hs = hs.copy(order='f') # [JMT-1,IMT-1]
	# These are the size of the psi grid, at the corners of the grid boxes.
	# This might need to change to the number of boxes in the future, which will
	# change all of the minus one's. Also would change the default values in the step call.
	# IMT=u.shape[3] # 670
	# JMT=v.shape[2] # 190
	# KM=u.shape[1] # 30

	ufluxinterp = uflux.copy()
	vfluxinterp = vflux.copy()
	# hsinterp = hs.copy()
	dxyzinterp = dxyz.copy()
	# Loop over nsteps iterations between two model outputs
	for i in xrange(nsteps):
		print j
		#  flux fields at starting time for this step
		if j != 0:
			ufluxinterp[:,:,:,0] = ufluxinterp[:,:,:,1]
			vfluxinterp[:,:,:,0] = vfluxinterp[:,:,:,1]
			# hsinterp[:,:,0] = hsinterp[:,:,1]
			dxyzinterp[:,:,:,0] = dxyzinterp[:,:,:,1]
			xstart = xend[j-1,:]
			ystart = yend[j-1,:]
			zstart = zend[j-1,:]
			ia = iend[j-1,:]
			ja = jend[j-1,:]
			ka = kend[j-1,:]
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
		# Find drifter locations
		# pdb.set_trace()
		# only send unmasked values to step
		if not np.ma.compressed(xstart).any(): # exit if all of the drifters have exited the domain
			break
		else:
			xend[j,ind],yend[j,ind],zend[j,ind],iend[j,ind],jend[j,ind],kend[j,ind],flag[ind] = \
														tracmass.step(np.ma.compressed(xstart),np.ma.compressed(ystart),
														np.ma.compressed(zstart),0.,np.ma.compressed(ia),np.ma.compressed(ja),
														np.ma.compressed(ka),tseas/float(nsteps),ufluxinterp,
														vfluxinterp,ff,dxyzinterp)#dz.data,dxdy)
			# if np.sum(flag)>0:
			# pdb.set_trace()
			t0[j+1] = t0[j] + tseas/float(nsteps) # update time in seconds to match drifters
			# if i == 8:
			# 	pdb.set_trace()
		j = j + 1

nc.close()
grid.close()

t0 = t0 + t0save # add back in base time in seconds

xg=np.concatenate((xstart0.reshape(1,xstart0.size),xend),axis=0)
yg=np.concatenate((ystart0.reshape(1,ystart0.size),yend),axis=0)
zg=np.concatenate((zstart0.reshape(1,zstart0.size),zend),axis=0)

# Recreate Cartesian particle locations from their index-relative locations
# just by interpolating
fxr = tri.nn_interpolator(xr.flatten())
fyr = tri.nn_interpolator(yr.flatten())
xp = fxr(xg,yg)
yp = fyr(xg,yg)

# Plot tracks
figure()
basemap.drawcoastlines()
basemap.fillcontinents('0.8')
basemap.drawparallels(np.arange(18, 35), dashes=(1, 0), linewidth=0.15, labels=[1, 0, 0, 0])
basemap.drawmeridians(np.arange(-100, -80), dashes=(1, 0), linewidth=0.15, labels=[0, 0, 0, 1])
hold('on')
contour(xr, yr, hfull, np.hstack(([10,20],np.arange(50,500,50))), colors='lightgrey', linewidths=0.15)
plot(xp[:,0],yp[:,0],'g.')
# Outline numerical domain
plot(xr[0,:],yr[0,:],'k:')
#ax1.plot(xr[-1,:],yr[-1,:],'k:')
plot(xr[:,0],yr[:,0],'k:')
plot(xr[:,-1],yr[:,-1],'k:')
show()

# if __name__=='__main__':
# 	run()
