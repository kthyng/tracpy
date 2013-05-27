import netCDF4 as netCDF
import glob
import numpy as np
from datetime import datetime, timedelta
import pdb

def setupROMSfiles(loc,date,ff,tout):
	'''
	setupROMSfiles()
	Kristen Thyng, March 2013

	Figures out necessary files to read in for track times and what
	model output indices within those files to use.

	Input:
	 loc 	File location
	 date 	datetime format start date
	 ff 	Time direction. ff=1 forward, ff=-1 backward
	 tout 	Number of model outputs to use

	Output:
	 nc 	NetCDF object for relevant files
	 tinds 	Indices of outputs to use from fname files
	'''
	# pdb.set_trace()
	if len(loc) == 2: # contains specific file and grid file
	# if 'http' in loc: # use just input file
		fname = loc[0]
		nc = netCDF.Dataset(fname)
		if ff == 1: #forward in time
			dates = nc.variables['ocean_time'][:]	
			ilow = date >= dates
			# time index with time value just below datenum_in (relative to file ifile)
			istart = dates[ilow].size - 1
			tinds = range(istart,istart+tout)
		else: #backward
			dates = nc.variables['ocean_time'][:]	
			ilow = date >= dates
			# time index with time value just below datenum_in (relative to file ifile)
			istart = dates[ilow].size - 1
			tinds = range(istart,istart-tout,-1)

	# This is for the case when we have a bunch of files to sort through
	else:
		files = np.sort(glob.glob(loc + 'ocean_his_????.nc')) # sorted list of file names
		# files = np.sort(glob.glob(loc + 'ocean_his_*_tochange.nc')) # this is for idealized tests

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
		# pdb.set_trace()
		fname = [files[ifile]]

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
		if ff==1:
			tinds = range(istart,istart+tout) # indices of model outputs desired
		else: # backward in time
			# have to shift istart since there are now new indices behind since going backward
			tinds = range(istart,istart-tout,-1)
		# If we need more indices than available in these files, add another

		if ff==1:
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

		# model output files together containing all necessary model outputs
		nc = netCDF.MFDataset(fname) # reopen since needed to close things in loop
	# pdb.set_trace()
	return nc, tinds