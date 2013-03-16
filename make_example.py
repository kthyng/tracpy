import numpy as np
import netCDF4 as netCDF
import pdb
from matplotlib.pyplot import *
from mpl_toolkits.basemap import Basemap
from matplotlib import delaunay
import op

# Make example situations to run with tracmass

### Which example to run ###
# example1: uniform cartesian grid with quiescent flows: runs correctly for ten steps between 
#  two model outputs. Single initialized drifter does not move.
# example2: uniform cartesian grid with uniform uflow in x: runs correctly for ten steps
#  between two model outputs. Single initialized drifter moves an equal amount east each step.
# example3: uniform cartesian grid with linearly increasing uflow in x: runs correctly for ten
#  steps between two model outputs. Single initialized drifter moves an increasing amount
#  east each step.
# example4: uniform cartesian grid with linearly increasing uflow in x and y: increasing
#  amount east and north each step
# example5: uniform cartesian grid with convergent horizontal flow in x: when x start position
#  is set to 51 (note the grid shift since in fortran indexing), drifter stays in place
# example6: curvilinear horizontal grid with quiescent flow
# later: changing z thicknesses
ex = 6

if ex == 1:
	rootgrp = netCDF.Dataset('example1.nc','w',format='NETCDF4')

	# Uniform Cartesian rectangular grid
	xl = 100 # size for psi grid
	yl = 200 # size for psi grid
	zl = 30
	tl = 3 # three time outputs
	dt = 4*3600 # 4 hours between outputs
	# x and y lengths are backward because of difference between python
	# and tracmass dimension ordering
	xpsi, ypsi = np.meshgrid(np.arange(0,yl),np.arange(0,xl))
	xr, yr = np.meshgrid(np.arange(0,yl-1),np.arange(0,xl-1))
	xu, yu = np.meshgrid(np.arange(0,yl-1),np.arange(0,xl))
	xv, yv = np.meshgrid(np.arange(0,yl),np.arange(0,xl-1))

	# Uniform vertical grid: [i,j,k,t]
	dxyz = np.ones((xl-1,yl-1,zl,tl))
	# dxyz = np.ones((zl,1)).reshape((1,1,zl,1)).repeat(xl,axis=0).repeat(yl,axis=1).repeat(tl,axis=3)

	# Quiescent flow in 4D: [i,j,k,t]
	U = np.zeros((xl,yl-1,zl,tl)) # use staggered grid still
	V = np.zeros((xl-1,yl,zl,tl)) # use staggered grid still

	# Time vector
	t = np.arange(0,3*dt,dt)

elif ex == 2:
	rootgrp = netCDF.Dataset('example2.nc','w',format='NETCDF4')

	# Uniform Cartesian rectangular grid
	xl = 100 # size for psi grid
	yl = 200 # size for psi grid
	zl = 30
	tl = 3 # three time outputs
	dt = (1/60.)*3600 # 1 minute between outputs
	xpsi, ypsi = np.meshgrid(np.arange(0,yl),np.arange(0,xl))
	xr, yr = np.meshgrid(np.arange(0,yl-1),np.arange(0,xl-1))
	xu, yu = np.meshgrid(np.arange(0,yl-1),np.arange(0,xl))
	xv, yv = np.meshgrid(np.arange(0,yl),np.arange(0,xl-1))

	# Uniform vertical grid: [i,j,k,t]
	# rho grid
	dz = np.ones((xl-1,yl-1,zl,tl))
	# u grid
	dzu = np.ones((xl,yl-1,zl,tl))
	# v grid
	dzv = np.ones((xl-1,yl,zl,tl))

	# Grid parameters
	dx = xr[0,1]-xr[0,0]
	dy = yr[1,0]-yr[0,0]
	dxv = np.ones(xv.shape)*dx
	dyu = np.ones(yu.shape)*dy
	dxdy = np.ones(xr.shape)*dx*dy
	dxyz = (dz.transpose(3,2,1,0)*dxdy.transpose(1,0)).transpose(3,2,1,0)

	# Linear flow in x in 4D: [i,j,k,t], use staggered grid still
	# velocities
	u = np.ones((xl,yl-1,zl,tl)) 
	v = np.zeros((xl-1,yl,zl,tl))
	# fluxes
	U = (u.transpose(3,2,1,0) * dyu.transpose(1,0) * dzu.transpose(3,2,1,0)).transpose(3,2,1,0)
	V = (v.transpose(3,2,1,0) * dxv.transpose(1,0) * dzv.transpose(3,2,1,0)).transpose(3,2,1,0)
	# pdb.set_trace()
	# Time vector
	t = np.arange(0,tl*dt,dt)

elif ex == 3:
	rootgrp = netCDF.Dataset('example3.nc','w',format='NETCDF4')

	# Uniform Cartesian rectangular grid
	xl = 100 # size for psi grid
	yl = 200 # size for psi grid
	zl = 30
	tl = 3 # three time outputs
	dt = (1/60.)*3600 # 1 minute between outputs
	xpsi, ypsi = np.meshgrid(np.arange(0,yl),np.arange(0,xl))
	xr, yr = np.meshgrid(np.arange(0,yl-1),np.arange(0,xl-1))
	xu, yu = np.meshgrid(np.arange(0,yl-1),np.arange(0,xl))
	xv, yv = np.meshgrid(np.arange(0,yl),np.arange(0,xl-1))

	# Uniform vertical grid: [i,j,k,t]
	# rho grid
	dz = np.ones((xl-1,yl-1,zl,tl))
	# u grid
	dzu = np.ones((xl,yl-1,zl,tl))
	# v grid
	dzv = np.ones((xl-1,yl,zl,tl))

	# Grid parameters
	dx = xr[0,1]-xr[0,0]
	dy = yr[1,0]-yr[0,0]
	dxv = np.ones(xv.shape)*dx
	dyu = np.ones(yu.shape)*dy
	dxdy = np.ones(xr.shape)*dx*dy
	dxyz = (dz.transpose(3,2,1,0)*dxdy.transpose(1,0)).transpose(3,2,1,0)

	# Linear flow in x in 4D: [i,j,k,t], use staggered grid still
	# velocities
	# u = np.ones((xl+1,yl,zl,tl)) 
	u = np.linspace(0,1,xl).reshape((xl,1,1,1)).repeat(zl,axis=2).repeat(yl-1,axis=1).repeat(tl,axis=3)
	v = np.zeros((xl-1,yl,zl,tl))
	# fluxes
	U = (u.transpose(3,2,1,0) * dyu.transpose(1,0) * dzu.transpose(3,2,1,0)).transpose(3,2,1,0)
	V = (v.transpose(3,2,1,0) * dxv.transpose(1,0) * dzv.transpose(3,2,1,0)).transpose(3,2,1,0)
	# pdb.set_trace()
	# Time vector
	t = np.arange(0,tl*dt,dt)

elif ex == 4:
	rootgrp = netCDF.Dataset('example4.nc','w',format='NETCDF4')

	# Uniform Cartesian rectangular grid
	xl = 100 # size for psi grid
	yl = 200 # size for psi grid
	zl = 30
	tl = 3 # three time outputs
	dt = (1/60.)*3600 # 1 minute between outputs
	xpsi, ypsi = np.meshgrid(np.arange(0,yl),np.arange(0,xl))
	xr, yr = np.meshgrid(np.arange(0,yl-1),np.arange(0,xl-1))
	xu, yu = np.meshgrid(np.arange(0,yl-1),np.arange(0,xl))
	xv, yv = np.meshgrid(np.arange(0,yl),np.arange(0,xl-1))

	# Uniform vertical grid: [i,j,k,t]
	# rho grid
	dz = np.ones((xl-1,yl-1,zl,tl))
	# u grid
	dzu = np.ones((xl,yl-1,zl,tl))
	# v grid
	dzv = np.ones((xl-1,yl,zl,tl))

	# Grid parameters
	dx = xr[0,1]-xr[0,0]
	dy = yr[1,0]-yr[0,0]
	dxv = np.ones(xv.shape)*dx
	dyu = np.ones(yu.shape)*dy
	dxdy = np.ones(xr.shape)*dx*dy
	dxyz = (dz.transpose(3,2,1,0)*dxdy.transpose(1,0)).transpose(3,2,1,0)

	# Linear flow in x in 4D: [i,j,k,t], use staggered grid still
	# velocities
	# u = np.ones((xl+1,yl,zl,tl)) 
	u = np.linspace(0,10,xl).reshape((xl,1,1,1)).repeat(zl,axis=2).repeat(yl-1,axis=1).repeat(tl,axis=3)
	v = np.linspace(0,10,yl).reshape((1,yl,1,1)).repeat(zl,axis=2).repeat(xl-1,axis=0).repeat(tl,axis=3)
	# fluxes
	U = (u.transpose(3,2,1,0) * dyu.transpose(1,0) * dzu.transpose(3,2,1,0)).transpose(3,2,1,0)
	V = (v.transpose(3,2,1,0) * dxv.transpose(1,0) * dzv.transpose(3,2,1,0)).transpose(3,2,1,0)
	# pdb.set_trace()
	# Time vector
	t = np.arange(0,tl*dt,dt)

elif ex == 5:
	rootgrp = netCDF.Dataset('example5.nc','w',format='NETCDF4')

	# Uniform Cartesian rectangular grid
	xl = 100 # size for psi grid
	yl = 200 # size for psi grid
	zl = 30
	tl = 3 # three time outputs
	dt = (1/60.)*3600 # 1 minute between outputs
	xpsi, ypsi = np.meshgrid(np.arange(0,yl),np.arange(0,xl))
	xr, yr = np.meshgrid(np.arange(0,yl-1),np.arange(0,xl-1))
	xu, yu = np.meshgrid(np.arange(0,yl-1),np.arange(0,xl))
	xv, yv = np.meshgrid(np.arange(0,yl),np.arange(0,xl-1))

	# Uniform vertical grid: [i,j,k,t]
	# rho grid
	dz = np.ones((xl-1,yl-1,zl,tl))
	# u grid
	dzu = np.ones((xl,yl-1,zl,tl))
	# v grid
	dzv = np.ones((xl-1,yl,zl,tl))

	# Grid parameters
	dx = xr[0,1]-xr[0,0]
	dy = yr[1,0]-yr[0,0]
	dxv = np.ones(xv.shape)*dx
	dyu = np.ones(yu.shape)*dy
	dxdy = np.ones(xr.shape)*dx*dy
	dxyz = (dz.transpose(3,2,1,0)*dxdy.transpose(1,0)).transpose(3,2,1,0)

	# Linear flow in x in 4D: [i,j,k,t], use staggered grid still
	# velocities
	# u = np.ones((xl+1,yl,zl,tl)) 
	u = np.linspace(-1,1,xl).reshape((xl,1,1,1)).repeat(zl,axis=2).repeat(yl-1,axis=1).repeat(tl,axis=3)
	v = np.zeros((xl-1,yl,zl,tl))
	# fluxes
	U = (u.transpose(3,2,1,0) * dyu.transpose(1,0) * dzu.transpose(3,2,1,0)).transpose(3,2,1,0)
	V = (v.transpose(3,2,1,0) * dxv.transpose(1,0) * dzv.transpose(3,2,1,0)).transpose(3,2,1,0)
	# pdb.set_trace()
	# Time vector
	t = np.arange(0,tl*dt,dt)

elif ex == 6:
	# Have to do real examples in ROMS ordering, then change for tracmass
	# using only non-ghost cells
	rootgrp = netCDF.Dataset('example6.nc','w',format='NETCDF4')

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

	# Shelf model curvilinear grid
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
	grid = netCDF.Dataset('/Users/kthyng/Dropbox/python/tamu/hab/grid.nc')
	lonu = grid.variables['lon_u'][1:-1,:]
	latu = grid.variables['lat_u'][1:-1,:]
	xu, yu = basemap(lonu,latu)
	lonv = grid.variables['lon_v'][:,1:-1]
	latv = grid.variables['lat_v'][:,1:-1]
	xv, yv = basemap(lonv,latv)
	hfull = grid.variables['h'][:]
	lonr = grid.variables['lon_rho'][1:-1,1:-1]
	latr = grid.variables['lat_rho'][1:-1,1:-1]
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
	grid.close()
	nc = netCDF.Dataset('/Users/kthyng/Documents/research/postdoc/ocean_his_0150.nc')
	t = nc.variables['ocean_time'][0:0+2]
	# Some grid metrics
	s = nc.variables['s_w'][:] # sigma coords, 31 layers
	cs = nc.variables['Cs_w'][:] # stretching curve in sigma coords, 31 layers
	nc.close()

	xl = xpsi.shape[1]
	yl = xpsi.shape[0]
	zl = len(cs)-1

	tl = 3 # three time outputs
	dt = 4*3600 # 4 hours between outputs
	# pdb.set_trace()

	# Change to tracmass ordering
	xr = xr.transpose(1,0)
	xu = xu.transpose(1,0)
	xv = xv.transpose(1,0)
	xpsi = xpsi.transpose(1,0)
	yr = yr.transpose(1,0)
	yu = yu.transpose(1,0)
	yv = yv.transpose(1,0)
	ypsi = ypsi.transpose(1,0)

	# Uniform vertical grid: [i,j,k,t]
	dxyz = np.ones((xl-1,yl-1,zl,tl))
	# dxyz = np.ones((zl,1)).reshape((1,1,zl,1)).repeat(xl,axis=0).repeat(yl,axis=1).repeat(tl,axis=3)

	# Quiescent flow in 4D: [i,j,k,t]
	
	U = np.zeros((xl,yl-1,zl,tl)) # use staggered grid still
	V = np.zeros((xl-1,yl,zl,tl)) # use staggered grid still

	# Time vector
	t = np.arange(0,3*dt,dt)

# Change all to Fortran ordering
xpsi = xpsi.copy(order='f')
xu = xu.copy(order='f')
xv = xv.copy(order='f')
xr = xr.copy(order='f')
ypsi = ypsi.copy(order='f')
yu = yu.copy(order='f')
yv = yv.copy(order='f')
yr = yr.copy(order='f')
dxyz = dxyz.copy(order='f')
U = U.copy(order='f')
V = V.copy(order='f')

## Make netCDF file
# Define dimensions
rootgrp.createDimension('xpsi',xl)
rootgrp.createDimension('xr',xl-1)
rootgrp.createDimension('ypsi',yl)
rootgrp.createDimension('yr',yl-1)
rootgrp.createDimension('zl',zl)
rootgrp.createDimension('tl',tl)
# Create variables
xpsis = rootgrp.createVariable('xpsi','f8',('xpsi','ypsi'))
xus = rootgrp.createVariable('xu','f8',('xpsi','yr'))
xvs = rootgrp.createVariable('xv','f8',('xr','ypsi'))
xrs = rootgrp.createVariable('xr','f8',('xr','yr'))
ypsis = rootgrp.createVariable('ypsi','f8',('xpsi','ypsi'))
yus = rootgrp.createVariable('yu','f8',('xpsi','yr'))
yvs = rootgrp.createVariable('yv','f8',('xr','ypsi'))
yrs = rootgrp.createVariable('yr','f8',('xr','yr'))
dxyzs = rootgrp.createVariable('dxyz','f8',('xr','yr','zl','tl'))
Us = rootgrp.createVariable('U','f8',('xpsi','yr','zl','tl'))
Vs = rootgrp.createVariable('V','f8',('xr','ypsi','zl','tl'))
ts = rootgrp.createVariable('t','f8',('tl',))
# Write data to netCDF variables
xpsis[:] = xpsi
xus[:] = xu
xvs[:] = xv
xrs[:] = xr
ypsis[:] = ypsi
yus[:] = yu
yvs[:] = yv
yrs[:] = yr
dxyzs[:] = dxyz
Us[:] = U
Vs[:] = V
ts[:] = t
figure()
pcolormesh(xu,yu,U[:,:,0,0])
title('Example ' + str(ex))
colorbar()
show()
