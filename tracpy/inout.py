"""
Input/output routines for tracpy.

Contains:
    setupROMSfiles
    readgrid
    readfields
    savetracks
    loadtracks
    loadtransport
    save_ll2grid
"""

import netCDF4 as netCDF
import glob
import numpy as np
from datetime import datetime, timedelta
import pdb
from matplotlib import delaunay
import octant
import time
from matplotlib.pyplot import *
import op
import os
import tracpy
from matplotlib.mlab import find

def setupROMSfiles(loc,date,ff,tout, tstride=1):
    '''
    setupROMSfiles()
    Kristen Thyng, March 2013

    Figures out necessary files to read in for track times and what
    model output indices within those files to use.

    Input:
     loc        File location. loc can be:
                * a thredds server web address
                * a single string of a file location
                * a list of strings of multiple file locations to be searched
                through
     date       datetime format start date
     ff         Time direction. ff=1 forward, ff=-1 backward
     tout       Number of model outputs to use
     tstride    Stride in time, in case want to use less model output
                 than is available. Default is 1, using all output.

    Output:
     nc         NetCDF object for relevant files
     tinds      Indices of outputs to use from fname files
    '''
    # This addresses an issue in netCDF4 that was then fixed, but
    # this line makes updating unnecessary. Issue described here: 
    # http://code.google.com/p/netcdf4-python/issues/detail?id=170
    netCDF._set_default_format(format='NETCDF3_64BIT')

    # For thredds server where all information is available in one place
    # or for a single file
    if 'http' in loc or type(loc)==str:
        nc = netCDF.Dataset(loc)
        if ff == 1: #forward in time
            dates = nc.variables['ocean_time'][:] # don't stride here, need all times to make index determinations
            ilow = date >= dates
            # time index with time value just below datenum_in (relative to file ifile)
            istart = dates[ilow].size - 1
            tinds = range(istart,istart+tout, tstride) #use tstride here to get tinds correct
        else: #backward
            dates = nc.variables['ocean_time'][:]   
            ilow = date >= dates
            # time index with time value just below datenum_in (relative to file ifile)
            istart = dates[ilow].size - 1
            tinds = range(istart,istart-tout, -tstride) #use tstride here to get tinds correct

    # This is for the case when we have a bunch of files to sort through
    else:
        # the globbing should happen ahead of time so this case looks different than
        # the single file case
        # files = np.sort(glob.glob(loc)) # sorted list of file names
        files = loc # file list is now to be input
        # pdb.set_trace()
        # Find the list of files that cover the desired time period
        # First, check to see if there is one or more than one time index in each file because
        # if there is only one, we need to compare two files to see where we are in time
        # So, open the first file to check.
        nctemp = netCDF.Dataset(files[0])
        ttemp = nctemp.variables['ocean_time'][:]
        nctemp.close()
        if ttemp.size==1: # one output/time step per file. 
            noutputs = 1
            # Assume it is safe to open them all to search for starting place. Rest of indices
            # are checked below.
            nctemp = netCDF.MFDataset(files)
            ttemp = nctemp.variables['ocean_time'][:]
            nctemp.close()
            # find which time steps date is between and take earlier one as the starting file
            ifile = find(ttemp<=date)[0]

        else: # more than one output per file  
            noutputs = 0          
            for i,name in enumerate(files): # Loop through files
                nctemp = netCDF.Dataset(name)
                ttemp = nctemp.variables['ocean_time'][:]
                # pdb.set_trace()
                nctemp.close()
                # need to check subsequent file as well, in case starting time is between files
                if i<(len(files)-1):
                    nctemp2 = netCDF.Dataset(files[i+1])
                    ttemp2 = nctemp2.variables['ocean_time'][:]
                    nctemp2.close()
                # If datenum_in is larger than the first time in the file but smaller
                # than the last time, then this is the correct file to use to start
                # or, if the starting date is between this file and the next, start 
                # with this file (for interpolation between the outputs)
                if (date >= ttemp[0] and date <= ttemp[-1]) or (date>ttemp[-1] and date<ttemp2[0]):
                    ifile = i # this is the starting file identifier then
                    break

        # Since the number of indices per file can change, make the process
        # of finding the necessary files a little more general
        # Start by opening two files
        i = 1
        fname = [files[ifile]]

        nc = netCDF.MFDataset(fname) # files in fname are in chronological order
        # Find which output in ifile is closest to the user-input start time (choose lower index)
        # Dates for drifters from this start date
        dates = nc.variables['ocean_time'][:]   
        ilow = date >= dates
        # time index with time value just below date (relative to file ifile)
        istart = dates[ilow].size - 1
        nc.close()
        # Select indices 
        if ff==1:
            tinds = range(istart,istart+tout, tstride) # indices of model outputs desired
        else: # backward in time
            # have to shift istart since there are now new indices behind since going backward
            tinds = range(istart,istart-tout, -tstride)
        # If we need more indices than available in these files, add another

        if ff==1:
            # If there is one output per file, just grab the number of necessary files
            if noutputs:
                fname = files[ifile:ifile+tout]
                nc = netCDF.MFDataset(fname) # files in fname are in chronological order

            else: # multiple snapshots per file
                # pdb.set_trace()
                # if the final index we want is beyond the length of these files,
                # keep adding files on
                # need an extra index for interpolation
                while tinds[-1]+1 > len(dates):
                # while len(tinds)+1 > len(dates): 
                    # if tdir: #forward - add 2nd file on end
                    fname.append(files[ifile+i])
                    nc = netCDF.MFDataset(fname) # files in fname are in chronological order
                    dates = nc.variables['ocean_time'][:]   
                    # ilow = date >= dates
                    # # time index with time value just below datenum_in (relative to file ifile)
                    # istart = dates[ilow].size - 1
                    # tinds = range(istart,istart+tout, tstride)
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
                tinds = range(istart,istart-tout, -tstride)
                nc.close()
                i = i + 1
        
        # model output files together containing all necessary model outputs
        nc = netCDF.MFDataset(fname) # reopen since needed to close things in loop

    return nc, tinds

def readgrid(grid_filename, vert_filename=None, llcrnrlon=-98.5, llcrnrlat=22.5, 
            urcrnrlon=-87.5, urcrnrlat=31.0, lat_0=30, lon_0=-94, res='i', 
            usebasemap=False, usespherical=True):
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
     grid_filename    File name (with extension) where grid information is stored
     vert_filename    (optional) File name (with extension) where vertical grid information
                    is stored, if not in grid_loc. Can also skip this if don't need 
                    vertical grid info.
     also optional basemap box parameters. Default is for full shelf model.
     usebasemap          (False) Whether to use load basemap into grid (True) or pyproj (False).
                    Basemap is slower but can be used for plotting, and pyproj is the opposite.

    Output:
     grid           Dictionary containing all necessary time-independent grid fields

    grid dictionary contains: (array sizing is for tracmass ordering)
     imt,jmt,km     Grid index sizing constants in (x,y,z), are for horizontal rho grid [scalar]
     dxv            Horizontal grid cell walls areas in x direction [imt,jmt-1]
     dyu            Horizontal grid cell walls areas in y direction [imt-1,jmt]
     dxdy           Horizontal area of cells defined at cell centers [imt,jmt]
     mask           Land/sea mask [imt,jmt] 
     pm,pn          Difference in horizontal grid spacing in x and y [imt,jmt]
     kmt            Number of vertical levels in horizontal space [imt,jmt]
     dzt0           Thickness in meters of grid at each k-level with time-independent free surface. 
                    Surface is at km [imt,jmt,km].
     zrt0           Depth in meters of grid at each k-level on vertical rho grid with time-independent 
                    free surface. Surface is at km [imt,jmt,km]
     zwt0           Depth in meters of grid at each k-level on vertical w grid with time-independent 
                    free surface. Surface is at km [imt,jmt,km]
     xr,yr          Rho grid zonal (x) and meriodional (y) coordinates [imt,jmt]
     xu,yu          U grid zonal (x) and meriodional (y) coordinates [imt,jmt]
     xv,yv          V grid zonal (x) and meriodional (y) coordinates [imt,jmt]
     xpsi,ypsi      Psi grid zonal (x) and meriodional (y) coordinates [imt,jmt]
     X,Y            Grid index arrays
     tri,trir       Delaunay triangulations
     Cs_r,sc_r      Vertical grid streching paramters [km-1]
     hc             Critical depth [scalar]
     h              Depths [imt,jmt]
     theta_s        Vertical stretching parameter [scalar]. A parameter 
                    (typically 0.0 <= theta_s < 5.0) that defines the amount 
                    of grid focusing. A higher value for theta_s will focus the grid more.
     theta_b        Vertical stretching parameter [scalar]. A parameter (0.0 < theta_b < 1.0) 
                    that says whether the coordinate will be focused at the surface 
                    (theta_b -> 1.0) or split evenly between surface and bottom (theta_b -> 0)
     basemap        Basemap object

    Note: all are in fortran ordering and tracmass ordering except for X, Y, tri, and tric
    To test: [array].flags['F_CONTIGUOUS'] will return true if it is fortran ordering
    '''
    keeptime = 0 # do timing for readgrid
    if keeptime: starttime = time.time()

    # Read in grid parameters and find x and y in domain on different grids
    # use full dataset to get grid information
    # This addresses an issue in netCDF4 that was then fixed, but
    # this line makes updating unnecessary. Issue described here:
    # http://code.google.com/p/netcdf4-python/issues/detail?id=170
    netCDF._set_default_format(format='NETCDF3_64BIT')
    gridfile = netCDF.Dataset(grid_filename)

    # # Read in whether grid is spherical or not
    # try:
    #     usesphericaltemp = gridfile.variables['spherical'][:]
    #     if usesphericaltemp == 'T':
    #         usespherical = True
    #     else:
    #         usespherical = False
    # except KeyError: # Assume not lon/lat if spherical flag is not in grid file
    #     usespherical = False

    # Basemap parameters.
    if usespherical:
        llcrnrlon=llcrnrlon; llcrnrlat=llcrnrlat; 
        urcrnrlon=urcrnrlon; urcrnrlat=urcrnrlat; projection='lcc'
        lat_0=lat_0; lon_0=lon_0; resolution=res; area_thresh=0.
        # pdb.set_trace()
        if usebasemap:
            from mpl_toolkits.basemap import Basemap
            basemap = Basemap(llcrnrlon=llcrnrlon,
                         llcrnrlat=llcrnrlat,
                         urcrnrlon=urcrnrlon,
                         urcrnrlat=urcrnrlat,
                         projection=projection,
                         lat_0=lat_0,
                         lon_0=lon_0,
                         resolution=resolution,
                         area_thresh=area_thresh)
        else:
            from pyproj import Proj
            # this gives somewhat different differences between projected coordinates as compared with previous basemap
            # definition for the default values.
            basemap = Proj(proj='lcc', lat_1=llcrnrlat, lat_2=urcrnrlat, lat_0=lat_0, lon_0=lon_0, x_0=0, y_0=0,ellps='clrk66',datum='NAD27')
            # basemap = (proj='lcc',lat_1=44.33333333333334,lat_2=46,lat_0=43.66666666666666, lon_0=-120.5,x_0=609601.2192024384, y_0=0,ellps='clrk66',datum='NAD27')
            # basemap = Proj("+proj=lcc +lat_0=lat_0 +lon_0=lon_0")
                            # +x_0=1700000 \
                            # +y_0=300000 \
                            # +no_defs \
                            # +a=6378137 \
                            # +rf=298.257222101 \
                            # +to_meter=1")
        if keeptime: 
            basemaptime = time.time()
            print "basemap time ", basemaptime - starttime
    else:
        basemap = []

    if usespherical:
        lonu = gridfile.variables['lon_u'][:]
        latu = gridfile.variables['lat_u'][:]
        xu, yu = basemap(lonu,latu)
        lonv = gridfile.variables['lon_v'][:]
        latv = gridfile.variables['lat_v'][:]
        xv, yv = basemap(lonv,latv)
        lonr = gridfile.variables['lon_rho'][:]#[1:-1,1:-1]
        latr = gridfile.variables['lat_rho'][:]#[1:-1,1:-1]
        xr, yr = basemap(lonr,latr)
        lonpsi = gridfile.variables['lon_psi'][:]
        latpsi = gridfile.variables['lat_psi'][:]
        xpsi, ypsi = basemap(lonpsi,latpsi)
    else:
        # read cartesian data
        xu = gridfile.variables['x_u'][:]
        yu = gridfile.variables['y_u'][:]
        xv = gridfile.variables['x_v'][:]
        yv = gridfile.variables['y_v'][:]
        xr = gridfile.variables['x_rho'][:]
        yr = gridfile.variables['y_rho'][:]
        xpsi = gridfile.variables['x_psi'][:]
        ypsi = gridfile.variables['y_psi'][:]

        # assign to spherical arrays
        lonu, latu = xu, yu
        lonv, latv = xv, yv
        lonr, latr = xr, yr
        lonpsi, latpsi = xpsi, ypsi

    # Create mask of all active cells if there isn't one already
    try:
        mask = gridfile.variables['mask_rho'][:]#[1:-1,1:-1]
    except KeyError:
        mask = np.ones_like(xr)

    pm = gridfile.variables['pm'][:]
    pn = gridfile.variables['pn'][:]
    h = gridfile.variables['h'][:]
    if keeptime:
        hgridtime = time.time()
        print "horizontal grid time ", hgridtime - basemaptime

    # Vertical grid metrics
    if 's_w' in gridfile.variables: # then given grid file contains vertical grid info
        sc_r = gridfile.variables['s_w'][:] # sigma coords, 31 layers
        Cs_r = gridfile.variables['Cs_w'][:] # stretching curve in sigma coords, 31 layers
        hc = gridfile.variables['hc'][:]
        theta_s = gridfile.variables['theta_s'][:]
        theta_b = gridfile.variables['theta_b'][:]
        if 'Vtransform' in gridfile.variables:
            Vtransform = gridfile.variables['Vtransform'][:]
            Vstretching = gridfile.variables['Vstretching'][:]
        else:
            Vtransform = 1
            Vstretching = 1
    # Still want vertical grid metrics, but are in separate file
    elif vert_filename is not None:
        nc = netCDF.Dataset(vert_filename)
        if 's_w' in nc.variables:
            sc_r = nc.variables['s_w'][:] # sigma coords, 31 layers
        else:
            sc_r = nc.variables['sc_w'][:] # sigma coords, 31 layers
        Cs_r = nc.variables['Cs_w'][:] # stretching curve in sigma coords, 31 layers
        hc = nc.variables['hc'][:]
        theta_s = nc.variables['theta_s'][:]
        theta_b = nc.variables['theta_b'][:]
        if 'Vtransform' in nc.variables:
            Vtransform = nc.variables['Vtransform'][:]
            Vstretching = nc.variables['Vstretching'][:]
        else:
            Vtransform = 1
            Vstretching = 1

    if keeptime: 
        vgridtime = time.time()
        print "vertical grid time ", vgridtime - hgridtime

    # make arrays in same order as is expected in the fortran code
    # ROMS expects [time x k x j x i] but tracmass is expecting [i x j x k x time]
    # change these arrays to be fortran-directioned instead of python
    # This is faster than copying arrays. To test: .flags['F_CONTIGUOUS']
    mask = np.asfortranarray(mask.T)
    xr = np.asfortranarray(xr.T)
    yr = np.asfortranarray(yr.T)
    xu = np.asfortranarray(xu.T)
    yu = np.asfortranarray(yu.T)
    xv = np.asfortranarray(xv.T)
    yv = np.asfortranarray(yv.T)
    xpsi = np.asfortranarray(xpsi.T)
    ypsi = np.asfortranarray(ypsi.T)
    lonr = np.asfortranarray(lonr.T)

    latr = np.asfortranarray(latr.T)
    lonu = np.asfortranarray(lonu.T)
    latu = np.asfortranarray(latu.T)
    lonv = np.asfortranarray(lonv.T)
    latv = np.asfortranarray(latv.T)
    lonpsi = np.asfortranarray(lonpsi.T)
    latpsi = np.asfortranarray(latpsi.T)
    pm = np.asfortranarray(pm.T)
    pn = np.asfortranarray(pn.T)
    h = np.asfortranarray(h.T)

    if keeptime: 
        fortranflippingtime = time.time()
        print "fortran flipping time ", fortranflippingtime - vgridtime

    # Basing this on setupgrid.f95 for rutgersNWA example project from Bror
    # Grid sizes
    imt = h.shape[0] # 671
    jmt = h.shape[1] # 191
    if 'sc_r' in dir():
        km = sc_r.shape[0]-1 # 30 NOT SURE ON THIS ONE YET

    if keeptime: 
        gridsizetime = time.time()
        print "grid size time ", gridsizetime - fortranflippingtime

    # Index grid, for interpolation between real and grid space
    # this is for psi grid, so that middle of grid is min + .5 value
    # # X goes from 0 to imt-2 and Y goes from 0 to jmt-2
    # Y, X = np.meshgrid(np.arange(jmt-1),np.arange(imt-1)) # grid in index coordinates, without ghost cells
    # This is for rho
    # X goes from 0 to imt-1 and Y goes from 0 to jmt-1
    Y, X = np.meshgrid(np.arange(jmt),np.arange(imt)) # grid in index coordinates, without ghost cells
    # Triangulation for grid space to curvilinear space
    tri = delaunay.Triangulation(X.flatten(),Y.flatten())
    # Triangulation for curvilinear space to grid space
    # pdb.set_trace()
    trir = delaunay.Triangulation(xr.flatten(),yr.flatten())
    trirllrho = delaunay.Triangulation(lonr.flatten(),latr.flatten())

    if keeptime: 
        delaunaytime = time.time()
        print "delaunay time ", delaunaytime - gridsizetime

    # tracmass ordering.
    # Not sure how to convert this to pm, pn with appropriate shift
    dxv = 1/pm # pm is 1/\Delta x at cell centers
    dyu = 1/pn # pn is 1/\Delta y at cell centers

    dxdy = dyu*dxv

    # Change dxv,dyu to be correct u and v grid size after having 
    # them be too big for dxdy calculation. This is not in the 
    # rutgersNWA example and I am not sure why. [i,j]
    dxv = 0.5*(dxv[:,:-1]+dxv[:,1:])
    dyu = 0.5*(dyu[:-1,:]+dyu[1:,:])
    # # These should be interpolated
    # dxv = dxv[:,:-1]
    # dyu = dyu[:-1,:]

    if keeptime: 
        gridmetricstime = time.time()
        print "grid metrics time ", gridmetricstime - delaunaytime

    # Adjust masking according to setupgrid.f95 for rutgersNWA example project from Bror
    if 'sc_r' in dir():
        mask2 = mask.copy()
        kmt = np.ones((imt,jmt),order='f')*km
        ind = (mask2==1)
        ind[0:imt-1,:] = ind[1:imt,:]
        mask2[ind] = 1
        ind = (mask2==1)
        ind[:,0:jmt-1] = ind[:,1:jmt]
        mask2[ind] = 1
        # ind = (mask[1:imt-1,:]==1)
        # mask[0:imt-2,ind] = 1
        # ind = (mask[:,1:imt-1]==1)
        # mask[:,0:jmt-2] = 1
        ind = (mask2==0)
        kmt[ind] = 0

        # Use octant to calculate depths/thicknesses for the appropriate vertical grid parameters
        # have to transform a few back to ROMS coordinates and python ordering for this
        zwt0 = octant.depths.get_zw(Vtransform, Vstretching, km+1, theta_s, theta_b, 
                        h.T.copy(order='c'), 
                        hc, zeta=0, Hscale=3)
        zrt0 = octant.depths.get_zrho(Vtransform, Vstretching, km, theta_s, theta_b, 
                        h.T.copy(order='c'), 
                        hc, zeta=0, Hscale=3)
        # Change dzt to tracmass/fortran ordering
        zwt0 = zwt0.T.copy(order='f')
        zrt0 = zrt0.T.copy(order='f')
        # this should be the base grid layer thickness that doesn't change in time because it 
        # is for the reference vertical level
        dzt0 = zwt0[:,:,1:] - zwt0[:,:,:-1]

    if keeptime: 
        calculatingdepthstime = time.time()
        print "calculating depths time ", calculatingdepthstime - gridmetricstime

    # Fill in grid structure
    if 'sc_r' in dir():
        grid = {'imt':imt,'jmt':jmt,'km':km,#'angle':angle, 
            'dxv':dxv,'dyu':dyu,'dxdy':dxdy, 
            'mask':mask,'kmt':kmt,'dzt0':dzt0,
            'zrt0':zrt0,'zwt0':zwt0,
            'pm':pm,'pn':pn,'tri':tri,'trir':trir,'trirllrho':trirllrho,
            'xr':xr,'xu':xu,'xv':xv,'xpsi':xpsi,'X':X,
            'yr':yr,'yu':yu,'yv':yv,'ypsi':ypsi,'Y':Y,
            'lonr':lonr,'lonu':lonu,'lonv':lonv,'lonpsi':lonpsi,
            'latr':latr,'latu':latu,'latv':yv,'latpsi':latpsi,
            'Cs_r':Cs_r,'sc_r':sc_r,'hc':hc,'h':h, 
            'theta_s':theta_s,'theta_b':theta_b,
            'Vtransform':Vtransform, 'Vstretching':Vstretching,
            'basemap':basemap}
    else:
        grid = {'imt':imt,'jmt':jmt, #'angle':angle,
            'dxv':dxv,'dyu':dyu,'dxdy':dxdy, 
            'mask':mask,
            'pm':pm,'pn':pn,'tri':tri,'trir':trir,'trirllrho':trirllrho,
            'xr':xr,'xu':xu,'xv':xv,'xpsi':xpsi,'X':X,
            'yr':yr,'yu':yu,'yv':yv,'ypsi':ypsi,'Y':Y,
            'lonr':lonr,'lonu':lonu,'lonv':lonv,'lonpsi':lonpsi,
            'latr':latr,'latu':latu,'latv':yv,'latpsi':latpsi,
            'h':h, 
            'basemap':basemap}
 
    if keeptime: 
        griddicttime = time.time()
        print "saving grid dict time ", griddicttime - calculatingdepthstime

    gridfile.close()

    return grid


def readfields(tind,grid,nc,z0=None, zpar=None, zparuv=None):
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
     tind   Single time index for model output to read in
     grid   Dictionary containing all necessary time-independent grid fields
     nc     NetCDF object for relevant files
     z0     (optional) if doing 2d isoslice, z0 contains string saying which kind
     zpar   (optional) if doing 2d isoslice, zpar is the depth/level/density
            at which we are to get the level
     zparuv (optional) Use this if the k index for the model output fields (e.g, u, v) is different
            from the k index in the grid. This might happen if, for example, only the surface current
            were saved, but the model run originally did have many layers. This parameter
            represents the k index for the u and v output, not for the grid.

    Output:
     uflux1     Zonal (x) flux at tind
     vflux1     Meriodional (y) flux at tind
     dzt        Height of k-cells in 3 dim in meters on rho vertical grid. [imt,jmt,km]
     zrt        Time-dependent depths of cells on vertical rho grid (meters).
                For the isoslice case, zrt ends up with 1 vertical level which contains
                the depths for the vertical center of the cell for that level.
     zwt        Time-dependent depths of cells on vertical w grid (meters). zwt always
                contains the depths at the vertical cell edges for the whole 3D grid
                and the correct depths can be accessed using the drifter indices.

    Array descriptions:
     u,v        Zonal (x) and meridional (y) velocities [imt,jmt,km] (m/s)
     ssh        Free surface [imt,jmt] (m)
     dz         Height of k-cells in 1 dim [km]
                From coord.f95: z coordinates (z>0 going up) for layers in meters 
                bottom layer: k=0; surface layer: k=KM and zw=0
                dz = layer thickness
     zt         Depths (negative) in meters of w vertical grid [imt,jmt,km+1]
     dzt        Height of k-cells in 3 dim in meters on rho vertical grid. [imt,jmt,km]
     dzt0       Height of k-cells in 2 dim. [imt,jmt]
     dzu        Height of each u grid cell [imt-1,jmt,km]
     dzv        Height of each v grid cell [imt,jmt-1,km]
     uflux1     Zonal (x) fluxes [imt-1,jmt,km] (m^3/s)?
     vflux1     Meriodional (y) fluxes [imt,jmt-1,km] (m^3/s)?
    '''

    # this parameter is in case there is less model output available vertically than
    # was actually run on the grid
    # pdb.set_trace()
    if zparuv is None:
        zparuv = zpar

    # tic_temp = time.time()
    # Read in model output for index tind
    if z0 == 's': # read in less model output to begin with, to save time
        u = nc.variables['u'][tind,zparuv,:,:] 
        v = nc.variables['v'][tind,zparuv,:,:]
        if 'zeta' in nc.variables:
            ssh = nc.variables['zeta'][tind,:,:] # [t,j,i], ssh in tracmass
            sshread = True
        else:
            sshread = False
    else:
        u = nc.variables['u'][tind,:,:,:] 
        v = nc.variables['v'][tind,:,:,:]
        if 'zeta' in nc.variables:
            ssh = nc.variables['zeta'][tind,:,:] # [t,j,i], ssh in tracmass
            sshread = True
        else:
            sshread = False

    h = grid['h'].T.copy(order='c')
    # Use octant to calculate depths for the appropriate vertical grid parameters
    # have to transform a few back to ROMS coordinates and python ordering for this
    if sshread:
        zwt = octant.depths.get_zw(grid['Vtransform'], grid['Vstretching'], grid['km']+1, grid['theta_s'], grid['theta_b'], 
                        h, grid['hc'], zeta=ssh, Hscale=3)
    else: # if ssh isn't available, approximate as 0
        zwt = octant.depths.get_zw(grid['Vtransform'], grid['Vstretching'], grid['km']+1, grid['theta_s'], grid['theta_b'], 
                        h, grid['hc'], zeta=0, Hscale=3)
    # Change dzt to tracmass/fortran ordering
    dzt = zwt[1:,:,:] - zwt[:-1,:,:]

    # also want depths on rho grid
    if sshread:
        zrt = octant.depths.get_zrho(grid['Vtransform'], grid['Vstretching'], grid['km'], grid['theta_s'], grid['theta_b'], 
                        h, grid['hc'], zeta=ssh, Hscale=3)
    else:
        zrt = octant.depths.get_zrho(grid['Vtransform'], grid['Vstretching'], grid['km'], grid['theta_s'], grid['theta_b'], 
                        h, grid['hc'], zeta=0, Hscale=3)

    dzu = .5*(dzt[:,:,0:grid['imt']-1] + dzt[:,:,1:grid['imt']])
    dzv = .5*(dzt[:,0:grid['jmt']-1,:] + dzt[:,1:grid['jmt'],:])

    # Change order back to ROMS/python for this calculation
    dyu = grid['dyu'].T.copy(order='c')
    dxv = grid['dxv'].T.copy(order='c')

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
        # Calculate flux and then take slice
        uflux1 = octant.tools.isoslice(u*dzu*dyu,op.resize(vert,2),zpar)
        vflux1 = octant.tools.isoslice(v*dzv*dxv,op.resize(vert,1),zpar)
        dzt = octant.tools.isoslice(dzt,vert,zpar)
        zrt = octant.tools.isoslice(zrt,vert,zpar)
    elif z0 == 'z':
        # Calculate flux and then take slice
        uflux1 = octant.tools.isoslice(u*dzu*dyu,op.resize(zrt,2),zpar)
        vflux1 = octant.tools.isoslice(v*dzv*dxv,op.resize(zrt,1),zpar)
        dzt = octant.tools.isoslice(dzt,zrt,zpar)
        zrt = np.ones(uflux1.shape)*zpar # array of the input desired depth

    # Change all back to tracmass/fortran ordering if being used again
    # This is faster than copying arrays
    # Don't bother changing order of these arrays since I have to change it in 
    # run.py anyway (concatenate appears not to preserve ordering)
    uflux1 = uflux1.T
    vflux1 = vflux1.T
    dzt = np.asfortranarray(dzt.T)
    zrt = np.asfortranarray(zrt.T)
    if sshread:
        ssh = np.asfortranarray(ssh.T)
    zwt = np.asfortranarray(zwt.T)

    # make sure that all fluxes have a placeholder for depth
    if is_string_like(z0):
        uflux1 = uflux1.reshape(np.append(uflux1.shape,1))
        vflux1 = vflux1.reshape(np.append(vflux1.shape,1))
        dzt = dzt.reshape(np.append(dzt.shape,1))
        zrt = zrt.reshape(np.append(zrt.shape,1))

    return uflux1, vflux1, dzt, zrt, zwt

def savetracks(xin, yin ,zpin, tpin, name, nstepsin, Nin, ffin, tseasin,
                ahin, avin, do3din, doturbin, locin, 
                doperiodicin, time_unitsin, T0in=None, Uin=None, Vin=None,
                savell=True):
    """
    Save tracks that have been calculated by tracmass into a netcdf file.

    Inputs:
        xin,yin,zpin        Drifter track positions [drifter x time]
        tpin                Time vector for drifters [drifter x time]
        name                Name of simulation, to use for saving file
        savell              Whether saving in latlon (True) or grid coords (False). Default True.
    """

    # name for ll is basic, otherwise add 'gc' to indicate as grid indices
    if not savell:
        name += 'gc'

    ntrac = xin.shape[0] # number of drifters
    nt = xin.shape[1] # number of time steps (with interpolation steps and starting point)
    
    # # save hash for the particular commit version that is currently being used
    # tempfile = os.popen('git log -1 --format="%H"')
    # git_hash_in = tempfile.read()
    # tempfile.close()
    # # remove \n on the end that I can't get rid of
    # git_hash_in = git_hash_in[0:git_hash_in.find('\n')]

    # Save file into a local directory called tracks. Make directory if it doesn't exist.
    if 'tracks' not in name:
        if not os.path.exists('tracks'):
            os.makedirs('tracks')
        name = 'tracks/' + name

    # Open file for writing.
    # Using netCDF3-Classic because the ROMS output does and
    # MFDataset does not work with NetCDF4
    # Can't save variables over 2GB without special format:
    # http://www.ncl.ucar.edu/Support/talk_archives/2011/0599.html
    # rootgrp = netCDF.Dataset('tracks/' + name + '.nc','w',format='NETCDF3_CLASSIC')
    # # Hoping that this format will both allow large variables and aggregation
    # rootgrp = netCDF.Dataset('tracks/' + name + '.nc','w',format='NETCDF3_64BIT')
    # Really useful netCDF4 resource: 
    # http://www.unidata.ucar.edu/software/netcdf/workshops/2012/netcdf_python/netcdf4python.pdf
    # Looks like I might still be able to use MFDataset (with netCDF4_CLASSIC files)
    # Apply compression at the createVariable stage with zlib
    # Info about classic: http://www.unidata.ucar.edu/software/netcdf/docs/netcdf/NetCDF_002d4-Classic-Model-Format.html
    # Looks like I might be able to use this, still use MFDataset, have large variables, and compress too
    # 4-Classic can still only have 1 unlimited dimension
    rootgrp = netCDF.Dataset(name + '.nc', 'w', format='NETCDF4_CLASSIC')

    # Define dimensions
    rootgrp.createDimension('ntrac', ntrac)
    rootgrp.createDimension('nt', nt)
    if Uin is not None:
        xul = Uin.shape[0]
        yul = Uin.shape[1]
        rootgrp.createDimension('xul',xul)
        rootgrp.createDimension('yul',yul)
        xvl = Vin.shape[0]
        yvl = Vin.shape[1]
        rootgrp.createDimension('xvl',xvl)
        rootgrp.createDimension('yvl',yvl)

    # Do the rest of this by variable so they can be deleted as I go for memory.
    if savell: # if saving in latlon
        # Create variable
        lonp = rootgrp.createVariable('lonp','f8',('ntrac','nt'), zlib=True) # 64-bit floating point, with lossless compression
        # Set some attributes
        lonp.long_name = 'longitudinal position of drifter'
        lonp.units = 'degrees'
        lonp.time = 'tp'
        # Write data to netCDF variables
        lonp[:] = xin
        # Delete to save space
        del(xin)

        latp = rootgrp.createVariable('latp','f8',('ntrac','nt'), zlib=True) # 64-bit floating point, with lossless compression
        latp.long_name = 'latitudinal position of drifter'
        latp.units = 'degrees'
        latp.time = 'tp'
        latp[:] = yin
        del(yin)
    else: # then saving in grid coordinates
        # Create variable
        xg = rootgrp.createVariable('xg','f8',('ntrac','nt'), zlib=True) # 64-bit floating point, with lossless compression
        # Set some attributes
        xg.long_name = 'x grid position of drifter'
        xg.units = 'grid units'
        xg.time = 'tp'
        # Write data to netCDF variables
        xg[:] = xin
        # Delete to save space
        del(xin)

        yg = rootgrp.createVariable('yg','f8',('ntrac','nt'), zlib=True) # 64-bit floating point, with lossless compression
        yg.long_name = 'y grid position of drifter'
        yg.units = 'grid units'
        yg.time = 'tp'
        yg[:] = yin
        del(yin)


    if do3din:
        zp = rootgrp.createVariable('zp','f8',('ntrac','nt'), zlib=True) # 64-bit floating point, with lossless compression
        zp.long_name = 'vertical position of drifter (negative is downward from surface)'
        zp.units = 'meter'
        zp.time = 'tp'
        zp[:] = zpin
        del(zpin)
    else:
        del(zpin)

    tp = rootgrp.createVariable('tp','f8',('ntrac','nt'), zlib=True) # 64-bit floating point, with lossless compression
    tp.long_name = 'time at drifter locations'
    tp.units = time_unitsin
    tp[:] = tpin
    del(tpin)

    if Uin is not None:
        T0 = rootgrp.createVariable('T0','f8',('ntrac'), zlib=True) # 64-bit floating point, with lossless compression
        U = rootgrp.createVariable('U','f8',('xul','yul'), zlib=True) # 64-bit floating point, with lossless compression
        V = rootgrp.createVariable('V','f8',('xvl','yvl'), zlib=True) # 64-bit floating point, with lossless compression
        T0.long_name = 'Initial volume transport associated with each drifter'
        U.long_name = 'Aggregation of x volume transports of drifters'
        V.long_name = 'Aggregation of y volume transports of drifters'
        T0.units = 'meter3 second-1'
        U.units = 'meter3 second-1'
        V.units = 'meter3 second-1'
        T0[:] = T0in
        U[:] = Uin
        V[:] = Vin
        del(T0in,Uin,Vin)

    # Create variables
    # Main track information
    # Include other run details
    nsteps = rootgrp.createVariable('nsteps','i4')
    N = rootgrp.createVariable('N','i4')
    ff = rootgrp.createVariable('ff','i4')
    tseas = rootgrp.createVariable('tseas','f8')
    ah = rootgrp.createVariable('ah','f8')
    av = rootgrp.createVariable('av','f8')
    do3d = rootgrp.createVariable('do3d','i4')
    doturb = rootgrp.createVariable('doturb','i4')
    doperiodic = rootgrp.createVariable('doperiodic','i4')
    # pdb.set_trace()
    # loc = rootgrp.createVariable('loc','i4')
    # git_hash = rootgrp.createVariable('git_hash','i4')

    # Set some attributes
    nsteps.long_name = 'sets max time steps between time interpolations between model outputs'
    N.long_name = 'sets number of samplings of drifter track'
    ff.long_name = 'forward (1) or backward (-1) in time'
    tseas.long_name = 'time between model outputs'
    ah.long_name = 'horizontal diffusion'
    av.long_name = 'vertical diffusion'
    do3d.long_name = 'flag for running in 3d (1) or 2d (0)'
    doturb.long_name = 'flag for using no subgrid parameterization (0), added turbulent velocities (1), displacement to particle position on a circle (2), displacement to particle position on an ellipse (3)'
    doperiodic.long_name = 'flag for using periodic boundary conditions: none (0), in x-direction (1), in y-direction (2)'
    # if len(locin) == 2:
    #     loc.long_name = 'location of model output information used for drifter experiment\n' + locin[0]
    # else:
    #     loc.long_name = 'location of model output information used for drifter experiment\n' + locin
    # git_hash.long_name = 'unique identifier for commit version of tracpy\n' + git_hash_in

    tseas.units = 'second'
    ah.units = 'meter2 second-1'
    av.units = 'meter2 second-1'

    # Write data to netCDF variables
    nsteps[:] = nstepsin
    N[:] = Nin
    ff[:] = ffin
    tseas[:] = tseasin
    ah[:] = ahin
    av[:] = avin
    do3d[:] = do3din
    doturb[:] = doturbin
    doperiodic[:] = doperiodicin
    #loc[:] = ''
    #git_hash[:] = ''

    rootgrp.close()

def loadtracks(name,loc=None):
    """
    Load in track info from netcdf file.
    Inputs:
        name    Name of tracks file
        loc     (optional) Tracks file is assumed to be in local tracks directory.
                Use this to give location if it is not.
    """

    if loc is None:
        nc = netCDF.Dataset('tracks/' + name + '.nc')
    else:
        nc = netCDF.Dataset(loc + '/' + name + '.nc')

    lonp = nc.variables['lonp'][:]
    latp = nc.variables['latp'][:]
    zp = nc.variables['zp'][:]
    tp = nc.variables['tp'][:]

    return lonp,latp,zp,tp      

def loadtransport(name,fmod=None):
    '''

    Inputs:
        name    Name of project
        fmod    File modifier: a way to choose a subset of the file in 
                the project directory instead of all. Should be a string and
                can include asterisks as wildcards.

    Outputs:
        U, V    Transport of drifter volume in x and y directions over all
                used simulation files
        lon0    Initial lon location for drifters
        lat0    Initial lat location for drifters
        T0      Overall
    '''

    # Which files to read in.
    if fmod is None:
        Files = glob.glob('tracks/' + name + '/*.nc')
    elif type(fmod) == list and len(fmod)>1:
        Files = []
        for i in xrange(len(fmod)):
            Files = Files + glob.glob('tracks/' + fmod[i])
    else:
        Files = glob.glob('tracks/' + name + '/' + fmod + '.nc')

    Files.sort()

    # Load in U and V volume transports of drifters and add together for
    # all files
    for i, File in enumerate(Files):
        d = netCDF.Dataset(File)
        if i == 0: # initialize U and V transports from first file
            U = d.variables['U'][:]
            V = d.variables['V'][:]
            T0 = np.sum(d.variables['T0'][:])
        else: # add in transports from subsequent simulations
            U = U + d.variables['U'][:]
            V = V + d.variables['V'][:]
            T0 = T0 + np.sum(d.variables['T0'][:])

        # Add initial drifter location (all drifters start at the same location)
        lon0 = d.variables['lonp'][:,0]
        lat0 = d.variables['latp'][:,0]
        d.close()

    return U, V, lon0, lat0, T0

def save_ll2grid(name, grid, loc=None):
    '''
    Input drifter tracks from saved file in grid coordinates and save a new file with
    drifter tracks in lat/lon instead.

    Example usage:
     loc = 'http://barataria.tamu.edu:8080/thredds/dodsC/NcML/txla_nesting6.nc' # TXLA model/grid output location
     grid = tracpy.inout.readgrid(loc)
     tracpy.inout.save_ll2grid([trackfile], grid, loc=loc)

     Note that [trackfile] should be the name of the drifter tracks files, including .nc extension,
     and any location prefix after 'tracks/'
     Note: input a loc value if the drifter files do not have it saved (those run on hafen, for example)
    '''

    # load in tracks
    d = netCDF.Dataset(name)
    lonp = d.variables['lonp'][:]
    latp = d.variables['latp'][:]

    # Convert to grid coords
    x, y, dt = tracpy.tools.interpolate2d(lonp, latp, grid, 'd_ll2ij')
    del(lonp, latp, grid)
    print dt

    if 'loc' in d.variables:
        loc = d.variables['loc'][:]
    else:
        print 'will use input loc value for saving to file'


    # save new file
    # transport calculation included
    if 'U' in d.variables:
        if d.variables['do3d'][:]:
            savetracks(x, y, d.variables['zp'][:], d.variables['tp'][:], name.split('/')[-1][:-3], d.variables['nsteps'][:], d.variables['N'][:],
                        d.variables['ff'][:], d.variables['tseas'][:], d.variables['ah'][:], d.variables['av'][:], 
                        d.variables['do3d'][:], d.variables['doturb'][:], loc, 
                        d.variables['T0'][:], d.variables['U'][:], d.variables['V'][:], savell=False)
        else: # have to input something for z but it won't be saved
            savetracks(x, y, y, d.variables['tp'][:], name.split('/')[-1][:-3], d.variables['nsteps'][:], d.variables['N'][:],
                        d.variables['ff'][:], d.variables['tseas'][:], d.variables['ah'][:], d.variables['av'][:], 
                        d.variables['do3d'][:], d.variables['doturb'][:], loc, 
                        d.variables['T0'][:], d.variables['U'][:], d.variables['V'][:], savell=False)
    else:
        if d.variables['do3d'][:]:
            savetracks(x, y, d.variables['zp'][:], d.variables['tp'][:], name.split('/')[-1][:-3], d.variables['nsteps'][:], d.variables['N'][:],
                        d.variables['ff'][:], d.variables['tseas'][:], d.variables['ah'][:], d.variables['av'][:], 
                        d.variables['do3d'][:], d.variables['doturb'][:], loc, savell=False)
        else: # have to input something for z but it won't be saved
            savetracks(x, y, y, d.variables['tp'][:], name.split('/')[-1][:-3], d.variables['nsteps'][:], d.variables['N'][:],
                        d.variables['ff'][:], d.variables['tseas'][:], d.variables['ah'][:], d.variables['av'][:], 
                        d.variables['do3d'][:], d.variables['doturb'][:], loc, savell=False)

    d.close()
