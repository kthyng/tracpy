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
from scipy.spatial import Delaunay
import matplotlib.tri as mtri
import octant
import time
import op
import os
import tracpy
from matplotlib.mlab import find


def setupROMSfiles(loc, date, ff, tout, time_units, tstride=1):
    """
    setupROMSfiles()
    Kristen Thyng, March 2013

    Figures out necessary files to read in for track times and what
    model output indices within those files to use.

    Args:
        loc: File location. loc can be a thredds server web address, a single
         string of a file location, a list of strings of multiple file
         locations to be searched through.
        date: datetime format start date
        ff: Time direction. ff=1 forward, ff=-1 backward
        tout: Number of model outputs to use
        time_units: To convert to datetime
        tstride: Stride in time, in case want to use less model output than
         is available. Default is 1, using all output.

    Returns:
        * nc - NetCDF object for relevant files
        * tinds - Indices of outputs to use from fname files
    """

    # For thredds server where all information is available in one place
    # or for a single file
    if 'http' in loc or type(loc) == str:
        nc = netCDF.Dataset(loc)

    # This is for the case when we have a bunch of files to sort through
    else:
        # the globbing should happen ahead of time so this case looks
        # different than the single file case
        # files in fname are in chronological order
        nc = netCDF.MFDataset(loc)

    # Convert date to number
    # dates = netCDF.num2date(nc.variables['ocean_time'][:], time_units)
    # The calendar definition extends dates to before the year 1582 for use
    # with idealized simulations without meaningful dates.
    dates = netCDF.num2date(nc.variables['ocean_time'][:], time_units,
                            calendar='proleptic_gregorian')
    # time index with time value just below date (relative to file ifile)
    istart = find(dates <= date)[-1]

    # Select indices
    if ff == 1:
        # indices of model outputs desired
        tinds = range(istart, istart+tout, tstride)
    else:  # backward in time
        # have to shift istart since there are now new indices behind since
        # going backward
        tinds = range(istart, istart-tout, -tstride)

    return nc, tinds


def readgrid(grid_filename, proj, vert_filename=None, usespherical=True):
    """
    readgrid(loc)
    Kristen Thyng, March 2013
    This function should be read in at the beginnind of a run.py call.
    It reads in all necessary grid information that won't change in time
    and stores it in a dictionary called grid.
    All arrays are changed to Fortran ordering (from Python ordering)
    and to tracmass variables ordering from ROMS ordering
    i.e. from [t,k,j,i] to [i,j,k,t]
    right away after reading in.

    Args:
        grid_filename: File name (with extension) where grid information is
            stored
        vert_filename (optional): File name (with extension) where vertical
            grid information is stored, if not in grid_loc. Can also skip
            this if don't need vertical grid info. also optional prjection
            box parameters. Default is None.
        proj: Projection object.
        usespherical: Use spherical geometric coordinates (lat/lon) or not.

    Returns:
     * grid - Dictionary containing all necessary time-independent grid fields

    grid dictionary contains: (array sizing is for tracmass ordering)
     * imt,jmt,km: Grid index sizing constants in (x,y,z), are for horizontal
       rho grid [scalar]
     * dxv: Horizontal grid cell walls areas in x direction [imt,jmt-1]
     * dyu: Horizontal grid cell walls areas in y direction [imt-1,jmt]
     * dxdy: Horizontal area of cells defined at cell centers [imt,jmt]
     * mask: Land/sea mask [imt,jmt]
     * pm,pn: Difference in horizontal grid spacing in x and y [imt,jmt]
     * kmt: Number of vertical levels in horizontal space [imt,jmt]
     * dzt0: Thickness in meters of grid at each k-level with
       time-independent free surface. Surface is at km [imt,jmt,km].
     * zrt0: Depth in meters of grid at each k-level on vertical rho grid
       with time-independent free surface. Surface is at km [imt,jmt,km]
     * zwt0: Depth in meters of grid at each k-level on vertical w grid with
       time-independent free surface. Surface is at km [imt,jmt,km]
     * xr, yr: Rho grid zonal (x) and meriodional (y) coordinates [imt,jmt]
     * xu, yu: U grid zonal (x) and meriodional (y) coordinates [imt,jmt]
     * xv, yv: V grid zonal (x) and meriodional (y) coordinates [imt,jmt]
     * xpsi, ypsi: Psi grid zonal (x) and meriodional (y) coordinates
       [imt, jmt]
     * X, Y: Grid index arrays
     * tri, trir: Delaunay triangulations
     * Cs_r, sc_r: Vertical grid streching paramters [km-1]
     * hc: Critical depth [scalar]
     * h: Depths [imt,jmt]
     * theta_s: Vertical stretching parameter [scalar]. A parameter
       (typically 0.0 <= theta_s < 5.0) that defines the amount of grid
       focusing. A higher value for theta_s will focus the grid more.
     * theta_b: Vertical stretching parameter [scalar]. A parameter (0.0 <
       theta_b < 1.0) that says whether the coordinate will be focused at the
       surface (theta_b -> 1.0) or split evenly between surface and bottom
       (theta_b -> 0)
     * basemap: Basemap object

    Note: all are in fortran ordering and tracmass ordering except for X, Y,
    tri, and tric
    To test: [array].flags['F_CONTIGUOUS'] will return true if it is fortran
    ordering
    """

    # Read in grid parameters and find x and y in domain on different grids
    # use full dataset to get grid information
    gridfile = netCDF.Dataset(grid_filename)

    if usespherical:
        try:
            lon_vert = gridfile.variables['lon_vert'][:]
            lat_vert = gridfile.variables['lat_vert'][:]
        except:
            lon_rho = gridfile.variables['lon_rho'][:]
            lat_rho = gridfile.variables['lat_rho'][:]
            x_rho, y_rho = proj(lon_rho, lat_rho)

            # get vertex locations
            try:
                angle = gridfile.variables['angle'][:]
            except:
                angle = np.zeros(x_rho.shape)
            x_vert, y_vert = octant.grid.rho_to_vert(x_rho, y_rho,
                                                     gridfile.variables['pm'][:],
                                                     gridfile.variables['pn'][:],
                                                     angle)
            lon_vert, lat_vert = proj(x_vert, y_vert, inverse=True)

        try:
            mask_rho = gridfile.variables['mask'][:]
        except:
            mask_rho = gridfile.variables['mask_rho'][:]

        grid = octant.grid.CGrid_geo(lon_vert, lat_vert, proj)
        grid.mask_rho = mask_rho
    else:  # read cartesian data
        try:
            x_vert = gridfile.variables['x_vert'][:]
            y_vert = gridfile.variables['y_vert'][:]
        except:
            x_rho = gridfile.variables['x_rho'][:]
            y_rho = gridfile.variables['y_rho'][:]

            # get vertex locations
            try:
                angle = gridfile.variables['angle'][:]
            except:
                angle = np.zeros(x_rho.shape)
            x_vert, y_vert = octant.grid.rho_to_vert(x_rho, y_rho,
                                                     gridfile.variables['pm'][:],
                                                     gridfile.variables['pn'][:],
                                                     angle)

        grid = octant.grid.CGrid(x_vert, y_vert)

        try:
            mask_rho = gridfile.variables['mask'][:]
            grid.mask_rho = mask_rho
        # except KeyError as 'mask':
        #     mask_rho = gridfile.variables['mask_rho'][:]
        #     grid.mask_rho = mask_rho
        except KeyError:
            print('No mask.')

        # Add into grid spherical coord variables so they are avaiable as
        # expected for the code but set them equal to the projected coords.
        # Make this better in the future.
        grid.lon_rho = grid.x_rho
        grid.lat_rho = grid.y_rho
        grid.lon_psi = grid.x_psi
        grid.lat_psi = grid.y_psi
        grid.lon_u = grid.x_u
        grid.lat_u = grid.y_u
        grid.lon_v = grid.x_v
        grid.lat_v = grid.y_v

    # vertical grid info
    if (vert_filename is not None) or ('s_w' in gridfile.variables):
        if 's_w' in gridfile.variables:  # test for presence of vertical info
            nc = gridfile
        else:
            nc = netCDF.Dataset(vert_filename)

        if 's_w' in nc.variables:
            grid.sc_r = nc.variables['s_w'][:]  # sigma coords, 31 layers
        else:
            grid.c_r = nc.variables['sc_w'][:]  # sigma coords, 31 layers
        # stretching curve in sigma coords, 31 layers
        grid.Cs_r = nc.variables['Cs_w'][:]
        grid.hc = nc.variables['hc'][:]
        grid.theta_s = nc.variables['theta_s'][:]
        grid.theta_b = nc.variables['theta_b'][:]
        if 'Vtransform' in nc.variables:
            grid.Vtransform = nc.variables['Vtransform'][:]
            grid.Vstretching = nc.variables['Vstretching'][:]
        else:
            grid.Vtransform = 1
            grid.Vstretching = 1

    # Basing this on setupgrid.f95 for rutgersNWA example project from Bror
    grid.h = gridfile.variables['h'][:]
    # Grid sizes
    grid.imt = grid.h.shape[1]  # 191
    grid.jmt = grid.h.shape[0]  # 671
    if hasattr(grid, 'sc_r'):
        grid.km = grid.sc_r.shape[0]-1  # 30 NOT SURE ON THIS ONE YET

    # Index grid, for interpolation between real and grid space
    # This is for rho
    # X goes from 0 to imt-1 and Y goes from 0 to jmt-1
    # grid in index coordinates, without ghost cells
    grid.X, grid.Y = np.meshgrid(np.arange(grid.imt), np.arange(grid.jmt))
    # Triangulation for grid space to curvilinear space
    pts = np.column_stack((grid.X.flatten(), grid.Y.flatten()))
    tess = Delaunay(pts)
    grid.tri = mtri.Triangulation(grid.X.flatten(), grid.Y.flatten(),
                                  tess.simplices.copy())
    # Triangulation for curvilinear space to grid space
    # Have to use SciPy's Triangulation to be more robust.
    # http://matveichev.blogspot.com/2014/02/matplotlibs-tricontour-interesting.html
    pts = np.column_stack((grid.x_rho.flatten(), grid.y_rho.flatten()))
    tess = Delaunay(pts)
    grid.trir = mtri.Triangulation(grid.x_rho.flatten(),
                                   grid.y_rho.flatten(),
                                   tess.simplices.copy())
    # For the two triangulations that are not integer based, need to
    # preprocess the mask to get rid of potential flat triangles at the
    # boundaries
    # http://matplotlib.org/1.3.1/api/tri_api.html#matplotlib.tri.TriAnalyzer
    # Hopefully these will work for other cases too: for the xy spherical
    # unit test cases, I needed these both for the triangulation to be valid.
    mask = mtri.TriAnalyzer(grid.trir).get_flat_tri_mask(0.01, rescale=True)
    grid.trir.set_mask(mask)
    mask = mtri.TriAnalyzer(grid.trir).get_flat_tri_mask(0.01, rescale=False)
    grid.trir.set_mask(mask)

    pts = np.column_stack((grid.lon_rho.flatten(), grid.lat_rho.flatten()))
    tess = Delaunay(pts)
    grid.trirllrho = mtri.Triangulation(grid.lon_rho.flatten(),
                                        grid.lat_rho.flatten(),
                                        tess.simplices.copy())
    mask = mtri.TriAnalyzer(grid.trirllrho).get_flat_tri_mask(0.01, rescale=True)
    grid.trirllrho.set_mask(mask)
    mask = mtri.TriAnalyzer(grid.trirllrho).get_flat_tri_mask(0.01, rescale=False)
    grid.trirllrho.set_mask(mask)

    # tracmass ordering.
    # Not sure how to convert this to pm, pn with appropriate shift
    grid.dxv = 1/grid.pm  # pm is 1/\Delta x at cell centers
    grid.dyu = 1/grid.pn  # pn is 1/\Delta y at cell centers

    grid.dxdy = grid.dyu*grid.dxv

    # Change dxv,dyu to be correct u and v grid size after having
    # them be too big for dxdy calculation. This is not in the
    # rutgersNWA example and I am not sure why. [i,j]
    grid.dxv = 0.5*(grid.dxv[:-1, :] + grid.dxv[1:, :])
    grid.dyu = 0.5*(grid.dyu[:, :-1] + grid.dyu[:, 1:])

    # Adjust masking according to setupgrid.f95 for rutgersNWA example
    # project from Bror
    if hasattr(grid, 'sc_r'):
        mask2 = grid.mask.copy()
        grid.kmt = np.ones((grid.jmt, grid.imt))*grid.km
        ind = (mask2 == 1)
        ind[0:grid.jmt-1, :] = ind[1:grid.jmt, :]
        mask2[ind] = 1
        ind = (mask2 == 1)
        ind[:, 0:grid.imt-1] = ind[:, 1:grid.imt]
        mask2[ind] = 1
        ind = (mask2 == 0)
        grid.kmt[ind] = 0

        # Use octant to calculate depths/thicknesses for the appropriate
        # vertical grid parameters have to transform a few back to ROMS
        # coordinates and python ordering for this
        grid.zwt0 = octant.depths.get_zw(grid.Vtransform, grid.Vstretching,
                                         grid.km+1, grid.theta_s,
                                         grid.theta_b, grid.h, grid.hc,
                                         zeta=0, Hscale=3)
        grid.zrt0 = octant.depths.get_zrho(grid.Vtransform, grid.Vstretching,
                                           grid.km, grid.theta_s,
                                           grid.theta_b, grid.h, grid.hc,
                                           zeta=0, Hscale=3)

        # this should be the base grid layer thickness that doesn't change in
        # time because it is for the reference vertical level
        grid.dzt0 = grid.zwt0[1:, :, :] - grid.zwt0[:-1, :, :]

    gridfile.close()

    return grid


def readfields(tind, grid, nc, z0=None, zpar=None, zparuv=None):
    """
    readfields()
    Kristen Thyng, March 2013

    Reads in model output in order to calculate fluxes and z grid
    properties to send into step.f95.
    Should be called initially and then subsequently each time loop.
    All arrays are changed to Fortran ordering (from Python ordering)
    and to tracmass variables ordering from ROMS ordering
    i.e. from [t,k,j,i] to [i,j,k,t]
    right away after reading in.

    Args:
        tind: Single time index for model output to read in
        grid: Dictionary containing all necessary time-independent grid
         fields
        nc: NetCDF object for relevant files
        z0 (Optional): if doing 2d isoslice, z0 contains string saying which
         kind
        zpar (Optional): if doing 2d isoslice, zpar is the
         depth/level/density at which we are to get the level
        zparuv (Optional): Use this if the k index for the model output
         fields (e.g, u, v) is different from the k index in the grid. This
         might happen if, for example, only the surface current were saved,
         but the model run originally did have many layers. This parameter
         represents the k index for the u and v output, not for the grid.

    Returns:
        * uflux1 - Zonal (x) flux at tind
        * vflux1 - Meriodional (y) flux at tind
        * dzt - Height of k-cells in 3 dim in meters on rho vertical grid.
          [imt,jmt,km]
        * zrt - Time-dependent depths of cells on vertical rho grid (meters).
          For the isoslice case, zrt ends up with 1 vertical level which
          contains the depths for the vertical center of the cell for that
          level.
        * zwt - Time-dependent depths of cells on vertical w grid (meters).
          zwt always contains the depths at the vertical cell edges for the
          whole 3D grid and the correct depths can be accessed using the
          drifter indices.

    Array descriptions:
     * u,v - Zonal (x) and meridional (y) velocities [imt,jmt,km] (m/s)
     * ssh - Free surface [imt,jmt] (m)
     * dz - Height of k-cells in 1 dim [km]
       From coord.f95: z coordinates (z>0 going up) for layers in meters
       bottom layer: k=0; surface layer: k=KM and zw=0
       dz = layer thickness
     * zt - Depths (negative) in meters of w vertical grid [imt,jmt,km+1]
     * dzt - Height of k-cells in 3 dim in meters on rho vertical grid.
       [imt,jmt,km]
     * dzt0 - Height of k-cells in 2 dim. [imt,jmt]
     * dzu - Height of each u grid cell [imt-1,jmt,km]
     * dzv - Height of each v grid cell [imt,jmt-1,km]
     * uflux1 - Zonal (x) fluxes [imt-1,jmt,km] (m^3/s)?
     * vflux1 - Meriodional (y) fluxes [imt,jmt-1,km] (m^3/s)?

    """

    # this parameter is in case there is less model output available
    # vertically than was actually run on the grid
    if zparuv is None:
        zparuv = zpar

    # tic_temp = time.time()
    # Read in model output for index tind
    if z0 == 's':  # read in less model output to begin with, to save time
        u = nc.variables['u'][tind, zparuv, :, :]
        v = nc.variables['v'][tind, zparuv, :, :]
        if 'zeta' in nc.variables:
            # [t,j,i], ssh in tracmass
            ssh = nc.variables['zeta'][tind, :, :]
            sshread = True
        else:
            sshread = False
    else:
        u = nc.variables['u'][tind, :, :, :]
        v = nc.variables['v'][tind, :, :, :]
        if 'zeta' in nc.variables:
            # [t,j,i], ssh in tracmass
            ssh = nc.variables['zeta'][tind, :, :]
            sshread = True
        else:
            sshread = False

    # Use octant to calculate depths for the appropriate vertical grid
    # parameters have to transform a few back to ROMS coordinates and python
    # ordering for this
    if sshread:
        zwt = octant.depths.get_zw(grid.Vtransform, grid.Vstretching,
                                   grid.km+1, grid.theta_s,
                                   grid.theta_b, grid.h, grid.hc, zeta=ssh,
                                   Hscale=3)
    else:  # if ssh isn't available, approximate as 0
        zwt = octant.depths.get_zw(grid.Vtransform, grid.Vstretching,
                                   grid.km+1, grid.theta_s,
                                   grid.theta_b, grid.h, grid.hc, zeta=0,
                                   Hscale=3)

    # Change dzt to tracmass/fortran ordering
    dzt = zwt[1:, :, :] - zwt[:-1, :, :]

    # also want depths on rho grid
    if sshread:
        zrt = octant.depths.get_zrho(grid.Vtransform, grid.Vstretching,
                                     grid.km, grid.theta_s,
                                     grid.theta_b, grid.h, grid.hc,
                                     zeta=ssh, Hscale=3)
    else:
        zrt = octant.depths.get_zrho(grid.Vtransform, grid.Vstretching,
                                     grid.km, grid.theta_s,
                                     grid.theta_b, grid.h, grid.hc, zeta=0,
                                     Hscale=3)

    dzu = .5*(dzt[:, :, 0:grid.imt-1] + dzt[:, :, 1:grid.imt])
    dzv = .5*(dzt[:, 0:grid.jmt-1, :] + dzt[:, 1:grid.jmt, :])

    # I think I can avoid this loop for the isoslice case
    if z0 is None:  # 3d case
        uflux1 = u*dzu*grid.dyu
        vflux1 = v*dzv*grid.dxv
    elif z0 == 's':  # want a specific s level zpar
        uflux1 = u*dzu[zpar, :, :]*grid.dyu
        vflux1 = v*dzv[zpar, :, :]*grid.dxv
        dzt = dzt[zpar, :, :]
        zrt = zrt[zpar, :, :]
    elif z0 == 'rho' or z0 == 'salt' or z0 == 'temp':
        # the vertical setup we're selecting an isovalue of
        vert = nc.variables[z0][tind, :, :, :]
        # Calculate flux and then take slice
        uflux1 = octant.tools.isoslice(u*dzu*grid.dyu, op.resize(vert, 2), zpar)
        vflux1 = octant.tools.isoslice(v*dzv*grid.dxv, op.resize(vert, 1), zpar)
        dzt = octant.tools.isoslice(dzt, vert, zpar)
        zrt = octant.tools.isoslice(zrt, vert, zpar)
    elif z0 == 'z':
        # Calculate flux and then take slice
        uflux1 = octant.tools.isoslice(u*dzu*grid.dyu, op.resize(zrt, 2), zpar)
        vflux1 = octant.tools.isoslice(v*dzv*grid.dxv, op.resize(zrt, 1), zpar)
        dzt = octant.tools.isoslice(dzt, zrt, zpar)
        zrt = np.ones(uflux1.shape)*zpar  # array of the input desired depth

    # make sure that all fluxes have a placeholder for depth
    if isinstance(z0, str):
        uflux1 = uflux1.reshape(np.append(1, uflux1.shape))
        vflux1 = vflux1.reshape(np.append(1, vflux1.shape))
        dzt = dzt.reshape(np.append(1, dzt.shape))
        zrt = zrt.reshape(np.append(1, zrt.shape))

    return uflux1, vflux1, dzt, zrt, zwt


def savetracks(xin, yin, zpin, tpin, name, nstepsin, Nin, ffin, tseasin,
               ahin, avin, do3din, doturbin, locin, doperiodicin,
               time_unitsin, T0in=None, Uin=None, Vin=None, savell=True):
    """
    Save tracks that have been calculated by tracmass into a netcdf file.

    Args:
        xin,yin,zpin: Drifter track positions [drifter x time]
        tpin: Time vector for drifters [drifter x time]
        name: Name of simulation, to use for saving file
        savell: Whether saving in latlon (True) or grid coords (False).
         Default True.
    """

    # name for ll is basic, otherwise add 'gc' to indicate as grid indices
    if not savell:
        name += 'gc'

    ntrac = xin.shape[0]  # number of drifters
    # number of time steps (with interpolation steps and starting point)
    nt = xin.shape[1]

    # SAVE VERSION OF TRACPY USED

    # Save file into a local directory called tracks. Make directory if it
    # doesn't exist.
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
        rootgrp.createDimension('xul', xul)
        rootgrp.createDimension('yul', yul)
        xvl = Vin.shape[0]
        yvl = Vin.shape[1]
        rootgrp.createDimension('xvl', xvl)
        rootgrp.createDimension('yvl', yvl)

    # Do the rest of this by variable so they can be deleted as I go for memory.
    if savell:  # if saving in latlon
        # Create variable
        # 64-bit floating point, with lossless compression
        lonp = rootgrp.createVariable('lonp', 'f8', ('ntrac', 'nt'),
                                      zlib=True)
        # Set some attributes
        lonp.long_name = 'longitudinal position of drifter'
        lonp.units = 'degrees'
        lonp.time = 'tp'
        # Write data to netCDF variables
        lonp[:] = xin
        # Delete to save space
        del(xin)

        # 64-bit floating point, with lossless compression
        latp = rootgrp.createVariable('latp', 'f8', ('ntrac', 'nt'),
                                      zlib=True)
        latp.long_name = 'latitudinal position of drifter'
        latp.units = 'degrees'
        latp.time = 'tp'
        latp[:] = yin
        del(yin)
    else:  # then saving in grid coordinates
        # Create variable
        # 64-bit floating point, with lossless compression
        xg = rootgrp.createVariable('xg', 'f8', ('ntrac', 'nt'), zlib=True)
        # Set some attributes
        xg.long_name = 'x grid position of drifter'
        xg.units = 'grid units'
        xg.time = 'tp'
        # Write data to netCDF variables
        xg[:] = xin
        # Delete to save space
        del(xin)

        # 64-bit floating point, with lossless compression
        yg = rootgrp.createVariable('yg', 'f8', ('ntrac', 'nt'), zlib=True)
        yg.long_name = 'y grid position of drifter'
        yg.units = 'grid units'
        yg.time = 'tp'
        yg[:] = yin
        del(yin)

    if do3din:
        # 64-bit floating point, with lossless compression
        zp = rootgrp.createVariable('zp', 'f8', ('ntrac', 'nt'),
                                    zlib=True)
        zp.long_name = 'vertical position of drifter (negative is downward from surface)'
        zp.units = 'meter'
        zp.time = 'tp'
        zp[:] = zpin
        del(zpin)
    else:
        del(zpin)

    # 64-bit floating point, with lossless compression
    tp = rootgrp.createVariable('tp', 'f8', ('ntrac', 'nt'),
                                zlib=True)
    tp.long_name = 'time at drifter locations'
    tp.units = time_unitsin
    tp[:] = tpin
    del(tpin)

    if Uin is not None:
        # 64-bit floating point, with lossless compression
        T0 = rootgrp.createVariable('T0', 'f8', ('ntrac'), zlib=True)
        U = rootgrp.createVariable('U', 'f8', ('xul', 'yul'), zlib=True)
        V = rootgrp.createVariable('V', 'f8', ('xvl', 'yvl'), zlib=True)
        T0.long_name = 'Initial volume transport associated with each drifter'
        U.long_name = 'Aggregation of x volume transports of drifters'
        V.long_name = 'Aggregation of y volume transports of drifters'
        T0.units = 'meter3 second-1'
        U.units = 'meter3 second-1'
        V.units = 'meter3 second-1'
        T0[:] = T0in
        U[:] = Uin
        V[:] = Vin
        del(T0in, Uin, Vin)

    # Create variables
    # Main track information
    # Include other run details
    nsteps = rootgrp.createVariable('nsteps', 'i4')
    N = rootgrp.createVariable('N', 'i4')
    ff = rootgrp.createVariable('ff', 'i4')
    tseas = rootgrp.createVariable('tseas', 'f8')
    ah = rootgrp.createVariable('ah', 'f8')
    av = rootgrp.createVariable('av', 'f8')
    do3d = rootgrp.createVariable('do3d', 'i4')
    doturb = rootgrp.createVariable('doturb', 'i4')
    doperiodic = rootgrp.createVariable('doperiodic', 'i4')

    # Set some attributes
    nsteps.long_name = 'sets max time steps between time interpolations \
                        between model outputs'
    N.long_name = 'sets number of samplings of drifter track'
    ff.long_name = 'forward (1) or backward (-1) in time'
    tseas.long_name = 'time between model outputs'
    ah.long_name = 'horizontal diffusion'
    av.long_name = 'vertical diffusion'
    do3d.long_name = 'flag for running in 3d (1) or 2d (0)'
    doturb.long_name = 'flag for using no subgrid parameterization (0), \
                        added turbulent velocities (1), displacement to \
                        particle position on a circle (2), displacement to \
                        particle position on an ellipse (3)'
    doperiodic.long_name = 'flag for using periodic boundary conditions: \
                            none (0), in x-direction (1), in y-direction (2)'

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

    rootgrp.close()


def loadtracks(name, loc=None):
    """
    Load in track info from netcdf file.

    Args:
        name (str): Name of tracks file
        loc (Optional): Tracks file is assumed to be in local tracks
         directory. Use this to give location if it is not.
    """

    if loc is None:
        nc = netCDF.Dataset('tracks/' + name + '.nc')
    else:
        nc = netCDF.Dataset(loc + '/' + name + '.nc')

    lonp = nc.variables['lonp'][:]
    latp = nc.variables['latp'][:]
    zp = nc.variables['zp'][:]
    tp = nc.variables['tp'][:]

    return lonp, latp, zp, tp


def loadtransport(name, fmod=None):
    """
    Args:
        name: Name of project
        fmod: File modifier: a way to choose a subset of the file in the
         project directory instead of all. Should be a string and can include
         asterisks as wildcards.

    Returns:
        * U, V - Transport of drifter volume in x and y directions over all
          used simulation files
        * lon0 - Initial lon location for drifters
        * lat0 - Initial lat location for drifters
        * T0 - Overall
    """

    # Which files to read in.
    if fmod is None:
        Files = glob.glob('tracks/' + name + '/*.nc')
    elif type(fmod) == list and len(fmod) > 1:
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
        if i == 0:  # initialize U and V transports from first file
            U = d.variables['U'][:]
            V = d.variables['V'][:]
            T0 = np.sum(d.variables['T0'][:])
        else:  # add in transports from subsequent simulations
            U = U + d.variables['U'][:]
            V = V + d.variables['V'][:]
            T0 = T0 + np.sum(d.variables['T0'][:])

        # Add initial drifter location (all drifters start at the same
        # location)
        lon0 = d.variables['lonp'][:, 0]
        lat0 = d.variables['latp'][:, 0]
        d.close()

    return U, V, lon0, lat0, T0


def save_ll2grid(name, grid, loc=None):
    """
    Input drifter tracks from saved file in grid coordinates and save a new
    file with drifter tracks in lat/lon instead.

    Example:
        >>> loc = 'http://barataria.tamu.edu:8080/thredds/dodsC/NcML/txla_nesting6.nc' # TXLA model/grid output location
        >>> grid = tracpy.inout.readgrid(loc)
        >>> tracpy.inout.save_ll2grid([trackfile], grid, loc=loc)

    Note:
        [trackfile] should be the name of the drifter tracks files,
        including .nc extension, and any location prefix after 'tracks/'
    Note:
        input a loc value if the drifter files do not have it saved (those
        run on hafen, for example)
    """

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
            savetracks(x, y, d.variables['zp'][:], d.variables['tp'][:],
                       name.split('/')[-1][:-3], d.variables['nsteps'][:],
                       d.variables['N'][:], d.variables['ff'][:],
                       d.variables['tseas'][:], d.variables['ah'][:],
                       d.variables['av'][:], d.variables['do3d'][:],
                       d.variables['doturb'][:], loc, d.variables['T0'][:],
                       d.variables['U'][:], d.variables['V'][:], savell=False)
        else:  # have to input something for z but it won't be saved
            savetracks(x, y, y, d.variables['tp'][:],
                       name.split('/')[-1][:-3], d.variables['nsteps'][:],
                       d.variables['N'][:], d.variables['ff'][:],
                       d.variables['tseas'][:], d.variables['ah'][:],
                       d.variables['av'][:], d.variables['do3d'][:],
                       d.variables['doturb'][:], loc, d.variables['T0'][:],
                       d.variables['U'][:], d.variables['V'][:], savell=False)
    else:
        if d.variables['do3d'][:]:
            savetracks(x, y, d.variables['zp'][:], d.variables['tp'][:],
                       name.split('/')[-1][:-3], d.variables['nsteps'][:],
                       d.variables['N'][:], d.variables['ff'][:],
                       d.variables['tseas'][:], d.variables['ah'][:],
                       d.variables['av'][:], d.variables['do3d'][:],
                       d.variables['doturb'][:], loc, savell=False)
        else:  # have to input something for z but it won't be saved
            savetracks(x, y, y, d.variables['tp'][:],
                       name.split('/')[-1][:-3], d.variables['nsteps'][:],
                       d.variables['N'][:], d.variables['ff'][:],
                       d.variables['tseas'][:], d.variables['ah'][:],
                       d.variables['av'][:], d.variables['do3d'][:],
                       d.variables['doturb'][:], loc, savell=False)

    d.close()
