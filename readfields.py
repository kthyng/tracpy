import netCDF4 as netCDF
import numpy as np
import pdb

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

    Array descriptions:
     u,v 		 	Zonal (x) and meridional (y) velocities [imt,jmt,km] (m/s)
     ssh 			Free surface [imt,jmt] (m)
     dzt 			Height of k-cells in 3 dim. [imt,jmt,km]
     dzt0 			Height of k-cells in 2 dim. [imt,jmt]
     dzu 	 		Height of each u grid cell [imt-1,jmt,km]
     dzv 	 		Height of each v grid cell [imt,jmt-1,km]
     uflux1 	 	Zonal (x) fluxes [imt-1,jmt,km] (m^3/s)?
     vflux1 	 	Meriodional (y) fluxes [imt,jmt-1,km] (m^3/s)?
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

	# Copy calculations from rutgersNWA/readfields.f95
	dzt = np.ones((grid['imt'],grid['jmt'],grid['km']))*np.nan
	dzu = np.ones((grid['imt']-1,grid['jmt'],grid['km']))*np.nan
	dzv = np.ones((grid['imt'],grid['jmt']-1,grid['km']))*np.nan
	for k in xrange(grid['km']):
		dzt0 = (grid['sc_r'][k]-grid['Cs_r'][k])*grid['hc'] \
				+ grid['Cs_r'][k]*grid['h']
		dzt[:,:,k] = dzt0 + ssh*(1.+dzt0/grid['h'])

	dzt0 = dzt[:,:,grid['km']-1]
	dzt[:,:,0:grid['km']-1] = dzt[:,:,1:grid['km']] - dzt[:,:,0:grid['km']-1]
	dzt[:,:,grid['km']-1] = ssh - dzt0
	dzu[0:grid['imt']-1,:,:] = dzt[0:grid['imt']-1,:,:]*0.5 + dzt[1:grid['imt'],:,:]*0.5
	dzv[:,0:grid['jmt']-1,:] = dzt[:,0:grid['jmt']-1,:]*0.5 + dzt[:,1:grid['jmt'],:]*0.5
	# pdb.set_trace()
	uflux1 = np.ones((grid['imt']-1,grid['jmt'],grid['km']))*np.nan
	vflux1 = np.ones((grid['imt'],grid['jmt']-1,grid['km']))*np.nan
	for k in xrange(grid['km']):
		# pdb.set_trace()
		uflux1[:,:,k] = u[:,:,k]*dzu[:,:,k]*grid['dyu']
		vflux1[:,:,k] = v[:,:,k]*dzv[:,:,k]*grid['dxv']

	# Flip vertical dimension because in ROMS surface is at k=-1 
	# and in tracmass surface is at 1
	uflux1 = uflux1[:,:,::-1]
	vflux1 = vflux1[:,:,::-1]
	dzt = dzt[:,:,::-1]
	# uflux1 = np.flipud(uflux1)
	# vflux1 = np.flipud(vflux1)
	# dzt = np.flipud(dzt)

	return ssh, uflux1, vflux1, dzt


	# # Creating full vertical grid in time # [t,k,j,i]
	# zeta4u = zetau.reshape(zetau.shape[0],1,zetau.shape[1],zetau.shape[2]) \
	# 					.repeat(len(cs),axis=1)
	# cs4u = cs.reshape(1,len(cs),1,1) \
	# 					.repeat(zetau.shape[0],axis=0) \
	# 					.repeat(zetau.shape[1],axis=2) \
	# 					.repeat(zetau.shape[2],axis=3)
	# h4u = hu.reshape(1,1,hu.shape[0],hu.shape[1]) \
	# 					.repeat(zetau.shape[0],axis=0) \
	# 					.repeat(len(cs),axis=1)	
	# z4u = zeta4u + cs4u*(zeta4u+h4u)
	# del(zeta4u,cs4u,h4u)
	# dz4u = np.diff(z4u,axis=1).reshape((u.shape)) # grid cell thickness in z coords
	# # Change velocity fields to fluxes.
	# # Calculate uflux, dependent on time. Calculating for two time indices at once.
	# # # Use pn for dy
	# # dyu = 1/op.resize(pn[1:-1,:],1) # [jmt-1,imt]
	# # uflux is size [2,km,jmt-1,imt]
	# uflux = u*dyu*dz4u #UPDATE THIS

	# # Do for vfluxes
	# zeta4v = zetav.reshape(zetav.shape[0],1,zetav.shape[1],zetav.shape[2]) \
	# 					.repeat(len(cs),axis=1)
	# cs4v = cs.reshape(1,len(cs),1,1) \
	# 					.repeat(zetav.shape[0],axis=0) \
	# 					.repeat(zetav.shape[1],axis=2) \
	# 					.repeat(zetav.shape[2],axis=3)
	# h4v = hv.reshape(1,1,hv.shape[0],hv.shape[1]) \
	# 					.repeat(zetav.shape[0],axis=0) \
	# 					.repeat(len(cs),axis=1)	
	# z4v = zeta4v + cs4v*(zeta4v+h4v)
	# del(zeta4v,cs4v,h4v)
	# dz4v = np.diff(z4v,axis=1).reshape((v.shape)) # grid cell thickness in z coords
	# # Calculate vflux, dependent on time. Calculating for two time indices at once.
	# # # Use pm for dx
	# # dxv = 1/op.resize(pm[:,1:-1],0) # [jmt,imt-1]
	# # vflux is size [2,km,jmt,imt-1]
	# vflux = v*dxv*dz4v #UPDATE THIS

	# # # dxdy is the horizontal area of the cell walls, so it should be the size of the number of cells only
	# # # I am calculating it by adding the area of the walls together for each cell
	# # dxdy = ((dyu[:,1:]+dyu[:,0:-1])+(dxv[1:,:]+dxv[0:-1,:])) # [jmt-1,imt-1]

	# # calculate dxyz here
	# dxyz = op.resize(dz4v,2)*dxdy #UPDATE THIS

	# hs = op.resize(zeta,1) # sea surface height in meters
	# dz = np.diff(z4,axis=1)[0,:,0,0] # just take one dz column for now

	# # make arrays in same order as is expected in the fortran code
	# # ROMS expects [time x k x j x i] but tracmass is expecting [i x j x k x time]
	# uflux = uflux.T
	# vflux = vflux.T
	# # dxdy = dxdy.T
	# dxyz = dxyz.T
	# # hs = hs.T

	# # change all arrays to be fortran-directioned instead of python
	# uflux = uflux.copy(order='f') # [2,km,jmt-1,imt] (all of these are swapped now due to transposes above)
	# vflux = vflux.copy(order='f') #[2,km,jmt,imt-1]]
	# # dxdy = dxdy.copy(order='f') #[jmt,imt]]
	# dxyz = dxyz.copy(order='f') #[2,km,jmt-1,imt-1]
	# # hs = hs.copy(order='f') # [jmt-1,imt-1]
	# # These are the size of the psi grid, at the corners of the grid boxes.
	# # This might need to change to the number of boxes in the future, which will
	# # change all of the minus one's. Also would change the default values in the step call.
	# # imt=u.shape[3] # 670
	# # jmt=v.shape[2] # 190
	# # km=u.shape[1] # 30
