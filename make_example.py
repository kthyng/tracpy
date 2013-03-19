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
# example6: curvilinear horizontal grid with quiescent flow: works
# example7: curvilinear horizontal grid with constant flow in x (y) (x and y): drifter moves 
#  uniformly in along- (across-) (along- and across-) shore direction
# example8: curvilinear horizontal grid with linearly increasing flow in x (y) (x and y): drifter moves 
#  increasingly in the along- (across-) (along- and across-) shore direction
# example9: curvilinear horizontal grid with realistic z grid and:
#  * quiescent flow: drifter doesn't move
#  * uniform along-shore velocity: drifter moves with the flow seemingly appropriately, though it goes through land
#  * uniform across-shore velocity: drifter moves seemingly appropriately toward the shore. 
#    It went through the barrier island, but stopped at the edge of the domain.
#  * uniform in along- and across-shore: moves ok but stopped on land on a barrier island
#  * linear in along-shore: ok
#  * linear in across-shore: ok
#  * linear in both shores: ok
ex = 9

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
	dxyz = np.ones((xl-1,yl-1,zl,tl))*3000
	# dxyz = np.ones((zl,1)).reshape((1,1,zl,1)).repeat(xl,axis=0).repeat(yl,axis=1).repeat(tl,axis=3)

	# Quiescent flow in 4D: [i,j,k,t]
	
	U = np.zeros((xl,yl-1,zl,tl)) # use staggered grid still
	V = np.zeros((xl-1,yl,zl,tl)) # use staggered grid still

	# Time vector
	t = np.arange(0,3*dt,dt)

elif ex == 7:
	# Have to do real examples in ROMS ordering, then change for tracmass
	# using only non-ghost cells
	rootgrp = netCDF.Dataset('example7.nc','w',format='NETCDF4')

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
	dxyz = np.ones((xl-1,yl-1,zl,tl))*3000
	# dxyz = np.ones((zl,1)).reshape((1,1,zl,1)).repeat(xl,axis=0).repeat(yl,axis=1).repeat(tl,axis=3)

	# Linear along-shore flow in 4D: [i,j,k,t]
	U = np.ones((xl,yl-1,zl,tl)) # use staggered grid still
	V = np.ones((xl-1,yl,zl,tl)) # use staggered grid still

	# Time vector
	t = np.arange(0,3*dt,dt)

elif ex == 8:
	# Have to do real examples in ROMS ordering, then change for tracmass
	# using only non-ghost cells
	rootgrp = netCDF.Dataset('example8.nc','w',format='NETCDF4')

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
	dxyz = np.ones((xl-1,yl-1,zl,tl))*3000
	# dxyz = np.ones((zl,1)).reshape((1,1,zl,1)).repeat(xl,axis=0).repeat(yl,axis=1).repeat(tl,axis=3)

	# Linear flow in x in 4D: [i,j,k,t], use staggered grid still (treating velocity like flux for now)
	# velocities
	# u = np.ones((xl+1,yl,zl,tl)) 
	U = np.linspace(0,100,xl).reshape((xl,1,1,1)).repeat(zl,axis=2).repeat(yl-1,axis=1).repeat(tl,axis=3)
	# V = np.zeros((xl-1,yl,zl,tl))
	V = np.linspace(0,100,yl).reshape((yl,1,1,1)).repeat(zl,axis=2).repeat(xl-1,axis=0).repeat(tl,axis=3)
	# U = np.zeros((xl,yl-1,zl,tl))

elif ex == 9:
	# Have to do real examples in ROMS ordering, then change for tracmass
	# using only non-ghost cells
	rootgrp = netCDF.Dataset('example9.nc','w',format='NETCDF4')

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
	grid = netCDF.Dataset('/Users/kthyng/Dropbox/python/tamu/hab/grid.nc')
	lonu = grid.variables['lon_u'][:]#[1:-1,:] I NEED THE EXTRA INFORMATION FOR GRID DIFFERENCES
	latu = grid.variables['lat_u'][:]#[1:-1,:]
	xu, yu = basemap(lonu,latu)
	lonv = grid.variables['lon_v'][:]#[:,1:-1]
	latv = grid.variables['lat_v'][:]#[:,1:-1]
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
	pm = grid.variables['pm'][:] # 1/dx
	pn = grid.variables['pn'][:] # 1/dy

	nc = netCDF.Dataset('/Users/kthyng/Documents/research/postdoc/ocean_his_0150.nc')
	t = nc.variables['ocean_time'][0:0+2]
	# Some grid metrics
	s = nc.variables['s_w'][:] # sigma coords, 31 layers
	cs = nc.variables['Cs_w'][:] # stretching curve in sigma coords, 31 layers

	xl = xpsi.shape[1]
	yl = xpsi.shape[0]
	zl = len(cs)-1

	tl = 3 # three time outputs
	dt = 4*3600 # 4 hours between outputs

	# Real vertical grid
	# Want h only within domain in x direction (no ghost cells) and interpolated onto cell
	# walls in y direction
	h = op.resize(grid.variables['h'][:,1:-1],0)
	# Get zeta. Same as with h (but there is a time dimension first)
	zeta = op.resize(nc.variables['zeta'][0:0+tl,:,1:-1],1) # [t,k,j,i] JUST DOING FIRST tl TIME STEPS
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
	dz4v = np.diff(z4,axis=1).reshape((tl,zl,yl,xl-1)) # grid cell thickness in z coords

	# in u
	h = op.resize(grid.variables['h'][1:-1,:],1)
	# Get zeta. Same as with h (but there is a time dimension first)
	zeta = op.resize(nc.variables['zeta'][0:0+tl,1:-1,:],2) # [t,k,j,i] JUST DOING FIRST tl TIME STEPS
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
	dz4u = np.diff(z4,axis=1).reshape((tl,zl,yl-1,xl)) # grid cell thickness in z coords

	# Use pn for dy
	dxv = 1/op.resize(pn[:,1:-1],0) # [JMT,IMT-1]
	# Use pn for dy
	dyu = 1/op.resize(pn[1:-1,:],1) # [JMT-1,IMT]
	# dxdy is the horizontal area of the cell walls, so it should be the size of the number of cells only
	# I am calculating it by adding the area of the walls together for each cell
	dxdy = ((dyu[:,1:]+dyu[:,0:-1])+(dxv[1:,:]+dxv[0:-1,:])) # [JMT-1,IMT-1]
	# calculate dxyz here
	dxyz = op.resize(dz4v,2)*dxdy

	grid.close()
	nc.close()

	# Now eliminate ghost cells in xu/xv/yu/yv arrays
	xu = xu[1:-1,:]
	yu = yu[1:-1,:]
	xv = xv[:,1:-1]
	yv = yv[:,1:-1]

	# # Quiescent flow in 4D: [i,j,k,t]	
	# U = np.zeros((tl,zl,yl-1,xl)) # use staggered grid still
	# V = np.zeros((tl,zl,yl,xl-1)) # use staggered grid still
	## OR ##
	# # Uniform flow in x in 4D: [i,j,k,t], use staggered grid still velocities
	# u = np.ones((tl,zl,yl-1,xl))*.1
	# U = u*dyu*dz4u
	# V = np.zeros((tl,zl,yl,xl-1)) # use staggered grid still
	## OR ##
	# # Uniform flow in y in 4D: [i,j,k,t], use staggered grid still velocities
	# U = np.zeros((tl,zl,yl-1,xl)) # use staggered grid still
	# v = np.ones((tl,zl,yl,xl-1))*.1
	# V = v*dxv*dz4v
	## OR ##
	# # Uniform flow in x and y in 4D: [i,j,k,t], use staggered grid still velocities
	# u = np.ones((tl,zl,yl-1,xl))*.1
	# U = u*dyu*dz4u
	# v = np.ones((tl,zl,yl,xl-1))*.1
	# V = v*dxv*dz4v
	## OR ##
	# # Linear flow in x in 4D: [i,j,k,t], use staggered grid still velocities
	# u = np.linspace(0,.1,xl).reshape((1,1,1,xl)).repeat(zl,axis=1).repeat(yl-1,axis=2).repeat(tl,axis=0)
	# U = u*dyu*dz4u
	# V = np.zeros((tl,zl,yl,xl-1)) # use staggered grid still
	## OR ##
	# Linear flow in y in 4D: [i,j,k,t], use staggered grid still velocities
	U = np.zeros((tl,zl,yl-1,xl)) # use staggered grid still
	v = np.linspace(0,.1,yl).reshape((1,1,yl,1)).repeat(zl,axis=1).repeat(xl-1,axis=3).repeat(tl,axis=0)
	V = v*dxv*dz4v
	## OR ##
	# # Linear flow in x and y in 4D: [i,j,k,t], use staggered grid still velocities
	# u = np.linspace(0,.1,xl).reshape((1,1,1,xl)).repeat(zl,axis=1).repeat(yl-1,axis=2).repeat(tl,axis=0)
	# U = u*dyu*dz4u
	# v = np.linspace(0,.1,yl).reshape((1,1,yl,1)).repeat(zl,axis=1).repeat(xl-1,axis=3).repeat(tl,axis=0)
	# V = v*dxv*dz4v

	# Change to tracmass ordering
	xr = xr.transpose(1,0)
	xu = xu.transpose(1,0)
	xv = xv.transpose(1,0)
	xpsi = xpsi.transpose(1,0)
	yr = yr.transpose(1,0)
	yu = yu.transpose(1,0)
	yv = yv.transpose(1,0)
	ypsi = ypsi.transpose(1,0)
	dxyz = dxyz.transpose(3,2,1,0)
	U = U.transpose(3,2,1,0)
	V = V.transpose(3,2,1,0)

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
rootgrp.close()
figure()
pcolormesh(xu,yu,U[:,:,0,0])
title('Example ' + str(ex))
colorbar()
show()
