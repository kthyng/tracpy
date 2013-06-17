import netCDF4 as netCDF
import numpy as np
import pdb
import octant
from matplotlib.pyplot import *
import op
import time

def readfields(tind,grid,nc,z0=None,zpar=None):
	'''
	readfields()
	Kristen Thyng, March 2013

	Reads in model output in order to calculate fluxes and z grid 
	properties to send into step.f95.
	Should be called initially and then subsequently each time loop.
    All arrays are changed to Fortran ordering (from Python ordering)
    and to tracmass variables ordering from ROMS ordering 
    i.e. from [t,k,j,i] to [i,j,k,t]
    right away after reading in.

	Input:
	 tind 	Single time index for model output to read in
     grid   Dictionary containing all necessary time-independent grid fields
	 nc 	NetCDF object for relevant files
	 z0 	(optional) if doing 2d isoslice, z0 contains string saying which kind
	 zpar 	(optional) if doing 2d isoslice, zpar is the depth/level/density
	 		at which we are to get the level

    Output:
     uflux1 	Zonal (x) flux at tind
     vflux1 	Meriodional (y) flux at tind
     dzt 		Height of k-cells in 3 dim in meters on rho vertical grid. [imt,jmt,km]
     zrt 		Time-dependent depths of cells on vertical rho grid (meters).
     			For the isoslice case, zrt ends up with 1 vertical level which contains
     			the depths for the vertical center of the cell for that level.
     zwt 		Time-dependent depths of cells on vertical w grid (meters). zwt always
     			contains the depths at the vertical cell edges for the whole 3D grid
     			and the correct depths can be accessed using the drifter indices.

    Array descriptions:
     u,v 		Zonal (x) and meridional (y) velocities [imt,jmt,km] (m/s)
     ssh 		Free surface [imt,jmt] (m)
     dz   		Height of k-cells in 1 dim [km]
     			From coord.f95: z coordinates (z>0 going up) for layers in meters 
				bottom layer: k=0; surface layer: k=KM and zw=0
				dz = layer thickness
	 zt  		Depths (negative) in meters of w vertical grid [imt,jmt,km+1]
     dzt 		Height of k-cells in 3 dim in meters on rho vertical grid. [imt,jmt,km]
     dzt0 		Height of k-cells in 2 dim. [imt,jmt]
     dzu 	 	Height of each u grid cell [imt-1,jmt,km]
     dzv 	 	Height of each v grid cell [imt,jmt-1,km]
     uflux1 	Zonal (x) fluxes [imt-1,jmt,km] (m^3/s)?
     vflux1 	Meriodional (y) fluxes [imt,jmt-1,km] (m^3/s)?
	'''

	# tic_temp = time.time()
	# Read in model output for index tind
	if z0 == 's': # read in less model output to begin with, to save time
		u = nc.variables['u'][tind,zpar,:,:] 
		v = nc.variables['v'][tind,zpar,:,:]
		ssh = nc.variables['zeta'][tind,:,:] # [t,j,i], ssh in tracmass
	else:
		u = nc.variables['u'][tind,:,:,:] 
		v = nc.variables['v'][tind,:,:,:]
		ssh = nc.variables['zeta'][tind,:,:] # [t,j,i], ssh in tracmass
	# time_read = time.time()-tic_temp

	# tic_temp = time.time()
	# # make arrays in same order as is expected in the fortran code
	# # ROMS expects [time x k x j x i] but tracmass is expecting [i x j x k x time]
	# # change these arrays to be fortran-directioned instead of python
	# # u = u.T.copy(order='f')
	# # v = v.T.copy(order='f')
	# # ssh = ssh.T.copy(order='f')
	# # # Flip vertical dimension because in ROMS surface is at k=-1 
	# # # and in tracmass surface is at 1
	# # u = np.flipud(u)
	# # v = np.flipud(v)
	# time_flip1 = time.time()-tic_temp

	# This is code from tracmass itself, which is being replaced by Rob's octant code
	# # Copy calculations from rutgersNWA/readfields.f95
	# dzt = np.ones((grid['imt'],grid['jmt'],grid['km']))*np.nan
	# dzu = np.ones((grid['imt']-1,grid['jmt'],grid['km']))*np.nan
	# dzv = np.ones((grid['imt'],grid['jmt']-1,grid['km']))*np.nan
	# for k in xrange(grid['km']):
	# 	dzt0 = (grid['sc_r'][k]-grid['Cs_r'][k])*grid['hc'] \
	# 			+ grid['Cs_r'][k]*grid['h']
	# 	dzt[:,:,k] = dzt0 + ssh*(1.+dzt0/grid['h'])

	# dzt0 = dzt[:,:,grid['km']-1]
	# dzt[:,:,0:grid['km']-1] = dzt[:,:,1:grid['km']] - dzt[:,:,0:grid['km']-1]
	# dzt[:,:,grid['km']-1] = ssh - dzt0

	# tic_temp = time.time()
	h = grid['h'].T.copy(order='c')
	# Use octant to calculate depths for the appropriate vertical grid parameters
	# have to transform a few back to ROMS coordinates and python ordering for this
	zwt = octant.depths.get_zw(1, 1, grid['km']+1, grid['theta_s'], grid['theta_b'], 
					h, grid['hc'], zeta=ssh, Hscale=3)
	# Change dzt to tracmass/fortran ordering
	# zwt = zwt.T.copy(order='f')
	# dzt = zwt[:,:,1:] - zwt[:,:,:-1]
	dzt = zwt[1:,:,:] - zwt[:-1,:,:]
	# pdb.set_trace()
	# time_zw = time.time()-tic_temp

	# tic_temp = time.time()
	# also want depths on rho grid
	zrt = octant.depths.get_zrho(1, 1, grid['km'], grid['theta_s'], grid['theta_b'], 
					h, grid['hc'], zeta=ssh, Hscale=3)
	# Change dzt to tracmass/fortran ordering
	# zrt = zrt.T.copy(order='f')
	# time_zr = time.time()-tic_temp

	# tic_temp = time.time()
	dzu = .5*(dzt[:,:,0:grid['imt']-1] + dzt[:,:,1:grid['imt']])
	dzv = .5*(dzt[:,0:grid['jmt']-1,:] + dzt[:,1:grid['jmt'],:])
	# dzu = .5*(dzt[0:grid['imt']-1,:,:] + dzt[1:grid['imt'],:,:])
	# dzv = .5*(dzt[:,0:grid['jmt']-1,:] + dzt[:,1:grid['jmt'],:])
	# dzu = dzt[0:grid['imt']-1,:,:]*0.5 + dzt[1:grid['imt'],:,:]*0.5
	# dzv = dzt[:,0:grid['jmt']-1,:]*0.5 + dzt[:,1:grid['jmt'],:]*0.5
	# dzu[0:grid['imt']-1,:,:] = dzt[0:grid['imt']-1,:,:]*0.5 + dzt[1:grid['imt'],:,:]*0.5
	# dzv[:,0:grid['jmt']-1,:] = dzt[:,0:grid['jmt']-1,:]*0.5 + dzt[:,1:grid['jmt'],:]*0.5
	# pdb.set_trace()
	# time_interp = time.time()-tic_temp

	# tic_temp = time.time()
	# Change order back to ROMS/python for this calculation
	# u = u.T.copy(order='c')
	# v = v.T.copy(order='c')
	# dzu = dzu.T.copy(order='c')
	# dzv = dzv.T.copy(order='c')
	# dzt = dzt.T.copy(order='c')
	dyu = grid['dyu'].T.copy(order='c')
	dxv = grid['dxv'].T.copy(order='c')
	# zrt = zrt.T.copy(order='c')
	# time_flip2 = time.time()-tic_temp
	# pdb.set_trace()

	# tic_temp = time.time()
	# I think I can avoid this loop for the isoslice case
	if z0 == None: # 3d case
		uflux1 = u*dzu*dyu
		vflux1 = v*dzv*dxv
	elif z0 == 's': # want a specific s level zpar
		uflux1 = u*dzu[zpar,:,:]*dyu
		vflux1 = v*dzv[zpar,:,:]*dxv
		dzt = dzt[zpar,:,:]
		zrt = zrt[zpar,:,:]
	elif z0 == 'rho' or z0 == 'salt' or z0 == 'temp':
		# the vertical setup we're selecting an isovalue of
		vert = nc.variables[z0][tind,:,:,:]
		# pdb.set_trace()
		# Calculate flux and then take slice
		uflux1 = octant.tools.isoslice(u*dzu*dyu,op.resize(vert,2),zpar)
		vflux1 = octant.tools.isoslice(v*dzv*dxv,op.resize(vert,1),zpar)
		dzt = octant.tools.isoslice(dzt,vert,zpar)
		zrt = octant.tools.isoslice(zrt,vert,zpar)
		# pdb.set_trace()
	elif z0 == 'z':
		# Calculate flux and then take slice
		uflux1 = octant.tools.isoslice(u*dzu*dyu,op.resize(zrt,2),zpar)
		vflux1 = octant.tools.isoslice(v*dzv*dxv,op.resize(zrt,1),zpar)
		dzt = octant.tools.isoslice(dzt,zrt,zpar)
		zrt = np.ones(uflux1.shape)*zpar # array of the input desired depth
	# time_flux = time.time()-tic_temp

	# tic_temp = time.time()
	# Change all back to tracmass/fortran ordering if being used again
	uflux1 = uflux1.T.copy(order='f')
	vflux1 = vflux1.T.copy(order='f')
	dzt = dzt.T.copy(order='f')
	zrt = zrt.T.copy(order='f')
	ssh = ssh.T.copy(order='f')
	zwt = zwt.T.copy(order='f')
	# time_flip3 = time.time()-tic_temp


	# tic_temp = time.time()
	# make sure that all fluxes have a placeholder for depth
	if is_string_like(z0):
		uflux1 = uflux1.reshape(np.append(uflux1.shape,1))
		vflux1 = vflux1.reshape(np.append(vflux1.shape,1))
		dzt = dzt.reshape(np.append(dzt.shape,1))
		zrt = zrt.reshape(np.append(zrt.shape,1))
	# time_reshape = time.time()-tic_temp


	# # Flip vertical dimension because in ROMS surface is at k=-1 
	# # and in tracmass surface is at 1
	# # This is not true. In tracmass, surface is at k=KM
	# uflux1 = uflux1[:,:,::-1]
	# vflux1 = vflux1[:,:,::-1]
	# dzt = dzt[:,:,::-1]
	# uflux1 = np.flipud(uflux1)
	# vflux1 = np.flipud(vflux1)
	# dzt = np.flipud(dzt)


	# print "time to read=",time_read
	# # print "time to flip=",time_flip1
	# print "time to get zw=",time_zw
	# print "time to get zr=",time_zr
	# print "time to interp=",time_interp
	# print "time to flip2=",time_flip2
	# print "time to flux=",time_flux
	# print "time to flip3=",time_flip3
	# print "time to reshape=",time_reshape

	return uflux1, vflux1, dzt, zrt, zwt
