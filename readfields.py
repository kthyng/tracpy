'''
Read fields and get things ready for a call to step.f95
'''

def readfields(ff,tind,grid):




	if ff: # forward
		u = nc.variables['u'][tind:tind+2,:,:,:] 
		v = nc.variables['v'][tind:tind+2,:,:,:]
		# u = nc.variables['u'][tind:tind+2,:,1:-1,:] 
		# v = nc.variables['v'][tind:tind+2,:,:,1:-1]
		zeta = nc.variables['zeta'][tind:tind+2,:,:] # [t,j,i], ssh in tracmass
		# zetau = op.resize(nc.variables['zeta'][tind:tind+2,1:-1,:],2) # [t,j,i]
		# zetav = op.resize(nc.variables['zeta'][tind:tind+2,:,1:-1],1) # [t,j,i]
	else: # NEED TO ADJUST THESE STILL
		u = nc.variables['u'][tind-1:tind-3,:,1:-1,:] 
		v = nc.variables['v'][tind-1:tind-3,:,:,1:-1]
		zetau = op.resize(nc.variables['zeta'][tind-1:tind-3,1:-1,:],2) # [t,j,i]
		zetav = op.resize(nc.variables['zeta'][tind-1:tind-3,:,1:-1],1) # [t,j,i]


	# Copy calculations from rutgersNWA/readfields.f95
	dzt = np.ones((KM,JMT,IMT))*np.nan
	for k in xrange(KM):
		dzt0 = (sc_r[k]-Cs_r[k])*hc + Cs_r[k]*h
		dzt[k,:,:] = dzt0 + zeta*(1.+dzt0/h)
	dzt0 = dzt[KM,:,:]
	dzt[0:KM-2,:,:] = dzt[1:KM-1,:,:] - dzt[0:KM-2,:,:]
	dzt[KM-1,:,:] = zeta - dzt0
	dzu[:,:,0:IMT-2] = dzt[:,:,0:IMT-2]*0.5 + dzt[:,:,1:IMT-1]*0.5
	dzv[:,0:JMT-2,:] = dzt[:,0:JMT-2,:]*0.5 + dzt[:,1:JMT-1,:]*0.5
	pdb.set_trace()
	# for k in xrange(KM):
	# 	uflux[]



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
	uflux = u*dyu*dz4u #UPDATE THIS

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
	vflux = v*dxv*dz4v #UPDATE THIS

	# # dxdy is the horizontal area of the cell walls, so it should be the size of the number of cells only
	# # I am calculating it by adding the area of the walls together for each cell
	# dxdy = ((dyu[:,1:]+dyu[:,0:-1])+(dxv[1:,:]+dxv[0:-1,:])) # [JMT-1,IMT-1]

	# calculate dxyz here
	dxyz = op.resize(dz4v,2)*dxdy #UPDATE THIS

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
