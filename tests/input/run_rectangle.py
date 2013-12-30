import numpy as np
import netCDF4 as netCDF
import tracpy
import tracpy.plotting
import datetime
import matplotlib.pyplot as plt

units = 'seconds since 1970-01-01'

# Location of example model output file and grid
loc = ['ocean_his_0001.nc','grid.nc']

# Start date in date time formatting
date = datetime.datetime(2013, 12, 19, 0)

# Time between outputs
tseas = 4*3600. # 4 hours between outputs, in seconds 

# Number of days to run the drifters.
ndays = tseas*9./(3600.*24)

# Sets a smaller limit than between model outputs for when to force interpolation if hasn't already occurred.
nsteps = 5

# Controls the sampling frequency of the drifter tracks.
N = 4

# Use ff = 1 for forward in time and ff = -1 for backward in time.
ff = 1

ah = 0. # m^2/s
av = 0. # m^2/s

# turbulence/diffusion flag
doturb = 0

nc = netCDF.Dataset('ocean_his_0001.nc')
grd = tracpy.inout.readgrid(loc,nc=nc)

lon0 = [-123.,-123.]
lat0 = [48.55,48.75]

# # Input starting locations as real space lon,lat locations
# lon0, lat0 = np.meshgrid(np.linspace(-98.5,-87.5,55), \
#                             np.linspace(22.5,31,49)) # whole domain, 20 km

# # Eliminate points that are outside domain or in masked areas
# lon0, lat0 = tracpy.tools.check_points(lon0, lat0, grd)

# for 3d flag, do3d=0 makes the run 2d and do3d=1 makes the run 3d
do3d = 0

## Choose method for vertical placement of drifters
z0 = 's' #'z' #'salt' #'s' 
zpar = 2 #-10 #grid['km']-1 # 30 #grid['km']-1

# simulation name, used for saving results into netcdf file
name = 'rectangle'

lonp, latp, zp, t, grd = tracpy.run.run(loc, nsteps, ndays, ff, date, tseas, ah, av, lon0, lat0, 
                                            z0, zpar, do3d, doturb, name, grid=grd, dostream=0, N=N)

plt.figure(figsize=(14,10))
plt.plot(grd['lonr'], grd['latr'], 'lightgrey')
plt.plot(grd['lonr'].T,grd['latr'].T,'lightgrey')
plt.plot(lonp.T,latp.T,'k', linewidth=5)
plt.plot(lonp[:,0], latp[:,0], 'go', markersize=10)
plt.plot(lonp[:,-1], latp[:,-1], 'ro', markersize=10)
plt.xlim(grd['lonr'].min(), grd['lonr'].max()) 
plt.ylim(grd['latr'].min(), grd['latr'].max())
plt.show()