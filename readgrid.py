

# LABEL ALL VARIABLES

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

# Read in grid parameters
# These include ghost cells so don't match u and v correctly
grid = netCDF.Dataset(loc + 'grid.nc')
lonu = grid.variables['lon_u'][:]
latu = grid.variables['lat_u'][:]
xu, yu = basemap(lonu,latu)
lonv = grid.variables['lon_v'][:]
latv = grid.variables['lat_v'][:]
xv, yv = basemap(lonv,latv)
hfull = grid.variables['h'][:]
lonr = grid.variables['lon_rho'][:]#[1:-1,1:-1]
latr = grid.variables['lat_rho'][:]#[1:-1,1:-1]
xr, yr = basemap(lonr,latr)
lonpsi = grid.variables['lon_psi'][:]
latpsi = grid.variables['lat_psi'][:]
xpsi, ypsi = basemap(lonpsi,latpsi)
mask = grid.variables['mask_rho'][:]#[1:-1,1:-1]
X, Y = np.meshgrid(np.arange(xr.shape[1]),np.arange(yr.shape[0])) # grid in index coordinates, without ghost cells
# Triangulation for grid space to curvilinear space
tri = delaunay.Triangulation(X.flatten(),Y.flatten())
# Triangulation for curvilinear space to grid space
tric = delaunay.Triangulation(xr.flatten(),yr.flatten())
# # Angle on rho grid
# theta = grid.variables['angle'][:]
# # Interpolate theta to be on psi grid
# theta = op.resize(op.resize(theta,0),1)
# pm = grid.variables['pm'][:] # 1/dx
# pn = grid.variables['pn'][:] # 1/dy
# The following are for setting up z array on u (v) grids for calculating u (v) fluxes
# Want h only within domain in y (x) direction (no ghost cells) and interpolated onto cell
# walls in x (y) direction
h = grid.variables['h'][:]
# hu = op.resize(grid.variables['h'][1:-1,:],1)
# hv = op.resize(grid.variables['h'][:,1:-1],0)

# Some grid metrics
sc_r = nc.variables['s_w'][:] # sigma coords, 31 layers
Cs_r = nc.variables['Cs_w'][:] # stretching curve in sigma coords, 31 layers
hc = nc.variables['hc'][:]
# Basing this on setupgrid.f95 for rutgersNWA example project from Bror
imt = h.shape[1] # 671
jmt = h.shape[0] # 191
km = sc_r.shape[0]-1 # 30
# ROMS ordering.
dxv = np.ones((JMT,IMT))*np.nan
dxv = xr
dxv[:,0:IMT-2] = dxv[:,1:IMT-1] - dxv[:,0:IMT-2]
dxv[:,IMT-1:IMT] = dxv[:,IMT-3:IMT-2]
dyu = np.ones((JMT,IMT))*np.nan
dyu = yr
dyu[0:JMT-2,:] = dyu[1:JMT-1,:] - dyu[0:JMT-2,:]
dyu[JMT-1,:] = dyu[JMT-2,:]
dxdy = dyu*dxv
# Adjust masking according to setupgrid.f95 for rutgersNWA example project from Bror
# my best guess is that kmt is an array containing the number of vertical
# grid cells in every (i,j) location. For sigma coords it doesn't change
# but in other codes it could. Also contains masking information.
kmt = np.ones((JMT,IMT))*KM
ind = (maskr[:,1:IMT-1]==1)
maskr[ind,0:IMT-2] = 1
ind = (maskr[1:IMT-1,:]==1)
maskr[0:JMT-2,:] = 1
ind = (maskr==0)
kmt[ind] = 0
# make arrays in same order as is expected in the fortran code
# ROMS expects [time x k x j x i] but tracmass is expecting [i x j x k x time]
# change these arrays to be fortran-directioned instead of python
dxdy = dxdy.T.copy(order='f') # DO THIS NOW OR LATER?
kmt = kmt.T.copy(order='f') # DO THIS NOW OR LATER?

# Fill in grid
grid = {'imt':imt,'jmt':jmt,'km':km, # Grid index sizing constants: (imt,jmt,km) in (x,y,z), are for horizontal rho grid
		'dxv':dxv,'dyu':dyu,'dxdy':dxdy, # horizontal grid cell walls area parameters
		'mask':mask, # land/sea mask, on rho grid
		'xr':xr,'xu':xu,'xv':xv,'xpsi':xpsi, # zonal (x) coordinates on various grids
		'yr':yr,'yu':yu,'yv':yv,'ypsi':ypsi, # meriodional (y) coordinates on various grids
		'Cs_r':Cs_r,'sc_r':sc_r,'hc',hc,'h':h} # Vertical grid parameters