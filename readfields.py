import netCDF4 as netCDF
import numpy as np
import pdb
from octant import depths

def readfields(tind,grid,nc):
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

    Output:
     ssh 		Free surface at tind
     uflux1 	Zonal (x) flux at tind
     vflux1 	Meriodional (y) flux at tind
     dzt 		Height of k-cells in 3 dim in meters on rho vertical grid. [imt,jmt,km]

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

	# Read in model output for index tind
	u = nc.variables['u'][tind,:,:,:] 
	v = nc.variables['v'][tind,:,:,:]
	ssh = nc.variables['zeta'][tind,:,:] # [t,j,i], ssh in tracmass

	# make arrays in same order as is expected in the fortran code
	# ROMS expects [time x k x j x i] but tracmass is expecting [i x j x k x time]
	# change these arrays to be fortran-directioned instead of python
	u = u.T.copy(order='f')
	v = v.T.copy(order='f')
	ssh = ssh.T.copy(order='f')
	# # Flip vertical dimension because in ROMS surface is at k=-1 
	# # and in tracmass surface is at 1
	# u = np.flipud(u)
	# v = np.flipud(v)

	# This is code from tracmass itself, which is being replaced by Rob's octant code
	# # Copy calculations from rutgersNWA/readfields.f95
	# dzt = np.ones((grid['imt'],grid['jmt'],grid['km']))*np.nan
	dzu = np.ones((grid['imt']-1,grid['jmt'],grid['km']))*np.nan
	dzv = np.ones((grid['imt'],grid['jmt']-1,grid['km']))*np.nan
	# for k in xrange(grid['km']):
	# 	dzt0 = (grid['sc_r'][k]-grid['Cs_r'][k])*grid['hc'] \
	# 			+ grid['Cs_r'][k]*grid['h']
	# 	dzt[:,:,k] = dzt0 + ssh*(1.+dzt0/grid['h'])

	# dzt0 = dzt[:,:,grid['km']-1]
	# dzt[:,:,0:grid['km']-1] = dzt[:,:,1:grid['km']] - dzt[:,:,0:grid['km']-1]
	# dzt[:,:,grid['km']-1] = ssh - dzt0

	# Use octant to calculate depths for the appropriate vertical grid parameters
	# have to transform a few back to ROMS coordinates and python ordering for this
	zt = depths.get_zw(1, 1, grid['km']+1, grid['theta_s'], grid['theta_b'], 
					grid['h'].T.copy(order='c'), 
					grid['hc'], zeta=ssh.T.copy(order='c'), Hscale=3)
	# Change dzt to tracmass/fortran ordering
	zt = zt.T.copy(order='f')
	dzt = zt[:,:,1:] - zt[:,:,:-1]
	# pdb.set_trace()

	dzu[0:grid['imt']-1,:,:] = dzt[0:grid['imt']-1,:,:]*0.5 + dzt[1:grid['imt'],:,:]*0.5
	dzv[:,0:grid['jmt']-1,:] = dzt[:,0:grid['jmt']-1,:]*0.5 + dzt[:,1:grid['jmt'],:]*0.5
	# pdb.set_trace()
	uflux1 = np.ones((grid['imt']-1,grid['jmt'],grid['km']))*np.nan
	vflux1 = np.ones((grid['imt'],grid['jmt']-1,grid['km']))*np.nan
	for k in xrange(grid['km']):
		# pdb.set_trace()
		uflux1[:,:,k] = u[:,:,k]*dzu[:,:,k]*grid['dyu']
		vflux1[:,:,k] = v[:,:,k]*dzv[:,:,k]*grid['dxv']

	# # Flip vertical dimension because in ROMS surface is at k=-1 
	# # and in tracmass surface is at 1
	# # This is not true. In tracmass, surface is at k=KM
	# uflux1 = uflux1[:,:,::-1]
	# vflux1 = vflux1[:,:,::-1]
	# dzt = dzt[:,:,::-1]
	# uflux1 = np.flipud(uflux1)
	# vflux1 = np.flipud(vflux1)
	# dzt = np.flipud(dzt)

	return ssh, uflux1, vflux1, dzt