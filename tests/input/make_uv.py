import numpy as np
import netCDF4 as netCDF
import datetime
# import pdb
# from matplotlib.pyplot import *
# from matplotlib import delaunay
# import op

# Make uniform eastward flow for rectangle simulation



# Read in grid for sizing info
grd = netCDF.Dataset('grid.nc')
xl = grd.variables['x_rho'][:].shape[1]
yl = grd.variables['x_rho'][:].shape[0]


# Add in vertical grid info
N = 3
s_rho = np.linspace(-0.975,-0.025,N)
s_w = np.linspace(-1,0,N+1)
hc = 0.
Cs_r = np.linspace(-0.975,-0.025,N)
Cs_w = np.linspace(-1,0,N+1)
theta_b = 0.
hmin = 100.
theta_s = 1e-4
tcline = 0.


# Time
units = 'seconds since 1970-01-01'
tl = 10 # three time outputs
dt = 4*3600 # 4 hours between outputs
startdate = netCDF.date2num(datetime.datetime(2013, 12, 19, 0), units)
t = np.arange(startdate, startdate + tl*dt, dt)


# Make velocity fields
# Linear flow in x in 4D: [i,j,k,t], use staggered grid still
u = 0.1*np.ones((tl,N,yl,xl-1)) 
v = np.zeros((tl,N,yl-1,xl))
# u = np.ones((xl,yl-1,N,tl)) 
# v = np.zeros((xl-1,yl,N,tl))


# Save file
rootgrp = netCDF.Dataset('ocean_his_0001.nc','w',format='NETCDF4')
# Define dimensions
rootgrp.createDimension('xpsi',xl-1)
rootgrp.createDimension('xr',xl)
rootgrp.createDimension('ypsi',yl-1)
rootgrp.createDimension('yr',yl)
rootgrp.createDimension('zl',N)
rootgrp.createDimension('zlp1',N+1)
rootgrp.createDimension('tl',tl)
# Create variables
# xpsis = rootgrp.createVariable('xpsi','f8',('xpsi','ypsi'))
# xus = rootgrp.createVariable('xu','f8',('xpsi','yr'))
# xvs = rootgrp.createVariable('xv','f8',('xr','ypsi'))
# xrs = rootgrp.createVariable('xr','f8',('xr','yr'))
# ypsis = rootgrp.createVariable('ypsi','f8',('xpsi','ypsi'))
# yus = rootgrp.createVariable('yu','f8',('xpsi','yr'))
# yvs = rootgrp.createVariable('yv','f8',('xr','ypsi'))
# yrs = rootgrp.createVariable('yr','f8',('xr','yr'))
# dxyzs = rootgrp.createVariable('dxyz','f8',('xr','yr','zl','tl'))
us = rootgrp.createVariable('u','f8',('tl','zl','yr','xpsi'))
vs = rootgrp.createVariable('v','f8',('tl','zl','ypsi','xr'))
ts = rootgrp.createVariable('ocean_time','f8',('tl',))
s_rhos = rootgrp.createVariable('s_rho','f8',('zl'))
s_ws = rootgrp.createVariable('s_w','f8',('zlp1'))
hcs = rootgrp.createVariable('hc','f8')
Cs_rs = rootgrp.createVariable('Cs_r','f8',('zl'))
Cs_ws = rootgrp.createVariable('Cs_w','f8',('zlp1'))
theta_bs = rootgrp.createVariable('theta_b','f8')
theta_ss = rootgrp.createVariable('theta_s','f8')
hmins = rootgrp.createVariable('hmin','f8')
tclines = rootgrp.createVariable('tcline','f8')
Vtransforms = rootgrp.createVariable('Vtransform','f8')
Vstretchings = rootgrp.createVariable('Vstretching','f8')

# Write data to netCDF variables
# xpsis[:] = xpsi
# xus[:] = xu
# xvs[:] = xv
# xrs[:] = xr
# ypsis[:] = ypsi
# yus[:] = yu
# yvs[:] = yv
# yrs[:] = yr
# dxyzs[:] = dxyz
us[:] = u
vs[:] = v
ts[:] = t
s_rhos[:] = s_rho
s_ws[:] = s_w
hcs[:] = hc
Cs_rs[:] = Cs_r
Cs_ws[:] = Cs_w
theta_bs[:] = theta_b
theta_ss[:] = theta_s
hmins[:] = hmin
tclines[:] = tcline
Vtransforms[:] = 1
Vstretchings[:] = 1
rootgrp.close()
