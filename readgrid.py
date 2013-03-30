import netCDF4 as netCDF
from mpl_toolkits.basemap import Basemap

def readgrid(loc):
    '''
    readgrid(loc)
    Kristen Thyng, March 2013
    This function should be read in at the beginnind of a run.py call.
    It reads in all necessary grid information that won't change in time
    and stores it in a dictionary called grid.
    All arrays are changed to Fortran ordering (from Python ordering)
    and to tracmass variables ordering from ROMS ordering 
    i.e. from [t,k,j,i] to [i,j,k,t]
    right away after reading in.

    Input:
     loc            File location

    Output:
     grid           Dictionary containing all necessary time-independent grid fields

    grid dictionary contains: (array sizing is for tracmass ordering)
     imt,jmt,km     Grid index sizing constants in (x,y,z), are for horizontal rho grid [scalar]
     dxv,dyu        Horizontal grid cell walls areas in x and y directions [imt,jmt]
     dxdy           Area of horizontal cell-walls [imt,jmt]
     mask           Land/sea mask [imt,jmt] 
     kmt            Number of vertical levels in horizontal space [imt,jmt]
     xr,yr          Rho grid zonal (x) and meriodional (y) coordinates [imt,jmt]
     xu,yu          U grid zonal (x) and meriodional (y) coordinates [imt,jmt]
     xv,yv          V grid zonal (x) and meriodional (y) coordinates [imt,jmt]
     xpsi,ypsi      Psi grid zonal (x) and meriodional (y) coordinates [imt,jmt]
     Cs_r,sc_r      Vertical grid streching paramters [km-1]
     hc             [scalar]
     h              Depths [imt,jmt]
     basemap        Basemap object
    '''

    # Basemap parameters.
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

    # Read in grid parameters and find x and y in domain on different grids
    grid = netCDF.Dataset(loc + 'grid.nc')
    lonu = grid.variables['lon_u'][:]
    latu = grid.variables['lat_u'][:]
    xu, yu = basemap(lonu,latu)
    lonv = grid.variables['lon_v'][:]
    latv = grid.variables['lat_v'][:]
    xv, yv = basemap(lonv,latv)
    lonr = grid.variables['lon_rho'][:]#[1:-1,1:-1]
    latr = grid.variables['lat_rho'][:]#[1:-1,1:-1]
    xr, yr = basemap(lonr,latr)
    lonpsi = grid.variables['lon_psi'][:]
    latpsi = grid.variables['lat_psi'][:]
    xpsi, ypsi = basemap(lonpsi,latpsi)
    mask = grid.variables['mask_rho'][:]#[1:-1,1:-1]

    # Vertical grid metrics
    sc_r = nc.variables['s_w'][:] # sigma coords, 31 layers
    Cs_r = nc.variables['Cs_w'][:] # stretching curve in sigma coords, 31 layers
    hc = nc.variables['hc'][:]
    h = grid.variables['h'][:]

    # make arrays in same order as is expected in the fortran code
    # ROMS expects [time x k x j x i] but tracmass is expecting [i x j x k x time]
    # change these arrays to be fortran-directioned instead of python
    dxv = dxv.T.copy(order='f')
    dyu = dyu.T.copy(order='f')
    dxdy = dxdy.T.copy(order='f')
    mask = mask.T.copy(order='f')
    kmt = kmt.T.copy(order='f')
    xr = xr.T.copy(order='f')
    yr = yr.T.copy(order='f')
    xu = xu.T.copy(order='f')
    yu = yu.T.copy(order='f')
    xv = xv.T.copy(order='f')
    yv = yv.T.copy(order='f')
    xpsi = xpsi.T.copy(order='f')
    ypsi = ypsi.T.copy(order='f')

    # Basing this on setupgrid.f95 for rutgersNWA example project from Bror
    # Grid sizes
    imt = h.shape[0] # 671
    jmt = h.shape[1] # 191
    km = sc_r.shape[0]-1 # 30


    # X, Y = np.meshgrid(np.arange(xr.shape[0]),np.arange(yr.shape[1])) # grid in index coordinates, without ghost cells
    # # Triangulation for grid space to curvilinear space
    # tri = delaunay.Triangulation(X.flatten(),Y.flatten())
    # # Triangulation for curvilinear space to grid space
    # tric = delaunay.Triangulation(xr.flatten(),yr.flatten())


    # tracmass ordering.
    dxv = xr
    dxv[0:imt-2,:] = dxv[1:imt-1,:] - dxv[0:imt-2,:]
    dxv[imt-1:imt,:] = dxv[imt-3:imt-2,:]

    dyu = yr
    dyu[:,0:jmt-2] = dyu[:,1:jmt-1] - dyu[:,0:jmt-2]
    dyu[:,jmt-1] = dyu[:,jmt-2]
    dxdy = dyu*dxv

    # Adjust masking according to setupgrid.f95 for rutgersNWA example project from Bror
    kmt = np.ones((imt,jmt))*km
    ind = (maskr[1:imt-1,:]==1)
    mask[0:imt-2,ind] = 1
    ind = (mask[:,1:imt-1]==1)
    mask[:,0:jmt-2] = 1
    ind = (mask==0)
    kmt[ind] = 0

    # Fill in grid structure
    grid = {'imt':imt,'jmt':jmt,'km':km, 
    	'dxv':dxv,'dyu':dyu,'dxdy':dxdy, 
    	'mask':mask,'kmt':kmt, 
    	'xr':xr,'xu':xu,'xv':xv,'xpsi':xpsi, 
    	'yr':yr,'yu':yu,'yv':yv,'ypsi':ypsi, 
    	'Cs_r':Cs_r,'sc_r':sc_r,'hc',hc,'h':h, 
        'basemap':basemap}

    return grid