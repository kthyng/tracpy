import netCDF4 as netCDF
import numpy as np

def savetracks(xpin,ypin,zpin,tpin,name):
	"""
	Save tracks that have been calculated by tracmass into a netcdf file.

	Inputs:
		xpin,ypin,zpin 	Drifter track positions [time x drifter]
		tpin			Time vector for drifters [time]
		name 			Name of simulation, to use for saving file
	"""

	ntrac = xpin.shape[1] # number of drifters
	nt = xpin.shape[0] # number of time steps (with interpolation steps and starting point)

	# Open file for writing
	rootgrp = netCDF.Dataset(name + '.nc','w',format='NETCDF4')

	# Define dimensions
	rootgrp.createDimension('ntrac',ntrac)
	rootgrp.createDimension('nt',nt)

	# Create variables
	xp = rootgrp.createVariable('xp','f8',('nt','ntrac')) # 64-bit floating point
	yp = rootgrp.createVariable('yp','f8',('nt','ntrac')) # 64-bit floating point
	zp = rootgrp.createVariable('zp','f8',('nt','ntrac')) # 64-bit floating point
	tp = rootgrp.createVariable('tp','f8',('nt')) # 64-bit floating point

	# Write data to netCDF variables
	xp[:] = xpin
	yp[:] = ypin
	zp[:] = zpin
	tp[:] = tpin

	rootgrp.close()
