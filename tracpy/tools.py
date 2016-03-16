"""
Tools for dealing with drifter stuff.

Functions include:

* interpolate2d
* interpolate3d
* find_final
* convert_indices
* check_points
* seed
"""

import numpy as np
from scipy import ndimage
import time
import matplotlib.tri as mtri
from matplotlib.mlab import Path, find


def interpolate2d(x, y, grid, itype, xin=None, yin=None, order=1,
                  mode='nearest', cval=0.):
    """
    Horizontal interpolation to map between coordinate transformations.

    Args:
        x, y: x, y
        grid: grid as read in by inout.readgrid()
        itype:
         * 'd_xy2ij' delaunay, from projected x, y to grid i, j
         * 'd_ij2xy' delaunay, from grid i, j to projected x, y
         * 'd_ll2ij' delaunay, from lon, lat to grid i, j
         * 'd_ij2ll' delaunay, from grid i, j to lon, lat
         * 'm_ij2xy' map_coordinates, from grid i, j to projected x, y or if
           z, xin, and yin are also input, from grid i, j, k to projected x,
           y, z. Can use the 3d version of this for transforming to lon/lat
           also if the xin/yin input are lon/lat arrays.
         * 'm_ij2ll' map_coordinates, from grid i, j to lon, lat
        xin: 3D array of x values that are mapped to the input x,y,z
         coordinates. This is only needed in the 3D mapping case. Normally,
         can just do this in 2D instead of 3D and get the same results.
        yin: 3D array of y values that are mapped to the input x,y,z
         coordinates. This is only needed in the 3D mapping case. Normally,
         can just do this in 2D instead of 3D and get the same results.
        order: order of interpolation for map_coordinates. 1 for linear and 3
         for cubic. Default=1
        mode: behavior for edge points. Default is 'nearest'. Notes on the
         map_coordinates function: The "mode" kwarg here just controls how
         the boundaries are treated mode='nearest' is _not_ nearest neighbor
         interpolation, it just uses the value of the nearest cell if the
         point lies outside the grid. The default is to treat the values
         outside the grid as zero, which can cause some edge effects if
         you're interpolating points near the edge. 'constant', 'nearest',
         'reflect' or 'wrap'. The "order" kwarg controls the order of the
         splines used. The default is cubic splines, order=3
        cval: Constant value used in map_coordinates if mode='constant'

    Returns:
        * xi,yi - Interpolated values
        * dt - Time required for interpolation
    """

    tic = time.time()

    if itype == 'd_xy2ij':
        # Set up functions for interpolating
        # changing format to use more robust triangulation in grid set up
        fx = mtri.LinearTriInterpolator(grid.trir, grid.X.flatten())
        fy = mtri.LinearTriInterpolator(grid.trir, grid.Y.flatten())
        # Need to shift indices to move from rho grid of interpolator to
        # arakawa c grid
        xi = fx(x, y) - .5
        yi = fy(x, y) - .5

    elif itype == 'd_ij2xy':
        # Set up functions for interpolating
        fx = mtri.LinearTriInterpolator(grid.tri, grid.x_rho.flatten())
        fy = mtri.LinearTriInterpolator(grid.tri, grid.y_rho.flatten())
        # Need to shift indices to move to rho grid of interpolator from
        # arakawa c grid
        xi = fx(x+0.5, y+0.5)
        yi = fy(x+0.5, y+0.5)

    elif itype == 'd_ll2ij':
        # Set up functions for interpolating
        fx = mtri.LinearTriInterpolator(grid.trirllrho, grid.X.flatten())
        fy = mtri.LinearTriInterpolator(grid.trirllrho, grid.Y.flatten())
        # Need to shift indices to move from rho grid of interpolator to
        # arakawa c grid
        xi = fx(x, y) - .5
        yi = fy(x, y) - .5

    elif itype == 'd_ij2ll':
        # Set up functions for interpolating
        fx = mtri.LinearTriInterpolator(grid.tri, grid.lon_rho.flatten())
        fy = mtri.LinearTriInterpolator(grid.tri, grid.lat_rho.flatten())
        # Need to shift indices to move to rho grid of interpolator from
        # arakawa c grid
        xi = fx(x+0.5, y+0.5)
        yi = fy(x+0.5, y+0.5)

    elif itype == 'm_ij2xy':
        # .5's are to shift from u/v grid to rho grid for interpolator
        xi = ndimage.map_coordinates(grid.x_rho.T, np.array([x.flatten()+.5,
                                                             y.flatten()+.5]),
                                     order=order, mode=mode,
                                     cval=cval).reshape(x.shape)
        yi = ndimage.map_coordinates(grid.y_rho.T, np.array([x.flatten()+.5,
                                                             y.flatten()+.5]),
                                     order=order, mode=mode,
                                     cval=cval).reshape(y.shape)

    elif itype == 'm_ij2ll':
        xi = ndimage.map_coordinates(grid.lon_rho.T, np.array([x.flatten()+.5,
                                                               y.flatten()+.5]),
                                     order=order, mode=mode,
                                     cval=cval).reshape(x.shape)
        yi = ndimage.map_coordinates(grid.lat_rho.T, np.array([x.flatten()+.5,
                                                               y.flatten()+.5]),
                                     order=order, mode=mode,
                                     cval=cval).reshape(y.shape)

    # Need to retain nan's since are changed them to zeros here
    if xi.size > 1:
        ind = np.isnan(x)
        xi[ind] = np.nan
        yi[ind] = np.nan

    dt = time.time() - tic

    return xi, yi, dt


def interpolate3d(x, y, z, zin, order=1, mode='nearest', cval=0.):
    """
    3D interpolation for transforming from grid/index space to whatever space
    is input with zin.

    Args:
        x,y,z: x, y, z coordinates
        zin: 3D array of z values that are mapped to the input x,y,z
         coordinates.
        order: order of interpolation for map_coordinates. 1 for linear and 3
         for cubic. Default=1
        mode: behavior for edge points. Default is 'nearest'. Notes on the
         map_coordinates function: The "mode" kwarg here just controls how
         the boundaries are treated mode='nearest' is _not_ nearest neighbor
         interpolation, it just uses the value of the nearest cell if the
         point lies outside the grid. The default is to treat the values
         outside the grid as zero, which can cause some edge effects if
         you're interpolating points near the edge. The "order" kwarg
         controls the order of the splines used. The default is cubic
         splines, order=3

    Returns:
        * zi - Interpolated values
        * dt - Time required for interpolation
    """

    tic = time.time()

    # Shift of .5 is assuming that input x/y are on a staggered grid frame
    # (counting from the cell edge and across the cell) but that the z values
    # are at the cell center, or rho locations.
    zi = ndimage.map_coordinates(zin, np.array([x.flatten()+.5,
                                                y.flatten()+.5,
                                                z.flatten()]),
                                 order=order, mode=mode,
                                 cval=cval).reshape(z.shape)

    # Need to retain nan's since are changed them to zeros here
    ind = np.isnan(z)
    zi[ind] = np.nan

    dt = time.time() - tic

    return zi, dt


def find_final(xp, yp, ind=-1):
    """
    Loop through drifters and find final location of drifters within the
    tracks arrays. This can be necessary because when drifters exit the
    numerical domain, they are nan'ed out. default is to search for the final
    non-nan location (-1), but can search for others instead, for example,
    the first non-nan position, which is helpful if we are looking at the
    flipped output from a backward run.
    """

    # Find final position for drifters (since they are nan'ed out after they
    # hit the open boundary)
    # Make this a separate function later
    xpc = []
    ypc = []
    for idrift in xrange(xp.shape[0]):
        # Find last non-nan and make sure it is in the desired month start
        # time
        ind3 = ~np.isnan(xp[idrift, :])
        # only plot if last non-nan (back in time) is in 1 month period
        # in order to plot the tracks that "started" in the plotted month
        if np.sum(ind3) > 1:  # don't want tracks that start on land
            # This is for if we care when the drifter stopped
            # if t[find(ind3)[ind]] >= datetime(year,startMonth,startDay,0) and \
            #   t[find(ind3)[ind]] <= datetime(year,startMonth+1,startDay,0):
            # ind2 = ~np.isnan(xp[idrift,:])
            # if there is a nan
            if np.sum(np.isnan(xp[idrift, :])) > 0 and \
               np.sum(np.isnan(xp[idrift, :])) < xp.shape[1]:
                # ax.plot(xp[idrift,find(ind2)[ind]].T,yp[idrift,find(ind2)[ind]].T,'o',color='orange',linewidth=.5,label='_nolegend_')
                xpc.append(xp[idrift, find(ind3)[ind]])
                ypc.append(yp[idrift, find(ind3)[ind]])
            else:
                # ax.plot(xp[idrift,ind].T,yp[idrift,ind].T,'o',color='orange',linewidth=.5,label='_nolegend_')
                xpc.append(xp[idrift, find(ind3)[ind]])
                ypc.append(yp[idrift, find(ind3)[ind]])

    return xpc, ypc


def convert_indices(direction, x, y):
    """
    Converts indices between Python and Fortran indexing, assuming that
    Python indexing begins at 0 and Fortran (for x and y) begins at 1.
    In Tracmass, the vertical indexing does begin at zero so this script
    does nothing to vertical indexing.

    Examples:
        For before a call to tracmass:

        >>> xstart,ystart = convert_indices('py2f',xstart,ystart)

        For after a call to tracmass:

        >>> xend,yend = convert_indices('f2py',xend,yend)
    """

    if direction == 'py2f':
        x = x+1
        y = y+1
    elif direction == 'f2py':
        x = x-1
        y = y-1

    return x, y


def check_points(lon0, lat0, grid, z0=None, nobays=False):
    """
    Eliminate starting locations for drifters that are outside numerical
    domain and that are masked out. If provided an array of starting vertical
    locations in z0, it checks whether these points are at least 1m above the
    bottom.

    Args:
        lon0,lat0: Starting locations for drifters in lon/lat
        z0: Starting locations for drifters in z
        grid: Grid made from readgrid.py
        nobays: Whether to use points in bays or not. Default is False.

    Returns:
        * lon0,lat0 - Fixed lon0,lat0
        * z0 - Fixed z0
    """

    lonr = grid.lon_rho
    latr = grid.lat_rho

    # make interpolation function for water depth h.
    # Used to check if float is above the bottom
    if z0 is not None:
        from scipy.interpolate import interp2d
        h = grid.h
        hint = interp2d(lonr[:, 1], latr[1, :], h.T, fill_value=np.nan)

    # If covering the whole domain, need to exclude points outside domain.
    # Use info just inside domain so points aren't right at the edge.
    xvert = np.hstack((np.flipud(lonr[1, :]), lonr[:, 1], lonr[-2, :],
                       np.flipud(lonr[:, -2])))
    yvert = np.hstack((np.flipud(latr[1, :]), latr[:, 1], latr[-2, :],
                       np.flipud(latr[:, -2])))
    verts = np.vstack((xvert, yvert))
    # Form path
    path = Path(verts.T)
    # Loop through particle tracks to eliminate drifters that start outside
    # the domain
    if lon0.ndim == 2:
        for jd in range(lon0.shape[0]):  # loop through drifters
            for it in range(lon0.shape[1]):
                # if drifter is not inside path, nan out this and all
                # subsequent points
                if not path.contains_point(np.vstack((lon0[jd, it],
                                                      lat0[jd, it]))):
                    lon0[jd, it] = np.nan
                    lat0[jd, it] = np.nan
                    if z0 is not None:
                        z0[jd, it] = np.nan

                if z0 is not None:
                    # check that the drifter starts above the bottom
                    if z0[jd, it] <= -1*hint(lon0[jd, it], lat0[jd, it]):
                        lon0[jd, it] = np.nan
                        lat0[jd, it] = np.nan
                        if z0 is not None:
                            z0[jd, it] = np.nan

    elif lon0.ndim == 1:
        for jd in range(lon0.shape[0]):  # loop through drifters
            # if drifter is not inside path, nan out this and all
            # subsequent points
            if not path.contains_point(np.vstack((lon0[jd], lat0[jd]))):
                lon0[jd] = np.nan
                lat0[jd] = np.nan
                if z0 is not None:
                    z0[jd] = np.nan

            if z0 is not None:
                # check that the drifter starts above the bottom
                if z0[jd] <= -1*hint(lon0[jd], lat0[jd]):
                    lon0[jd] = np.nan
                    lat0[jd] = np.nan
                    z0[jd] = np.nan

    # Also nan out points that are masked
    fmask = mtri.LinearTriInterpolator(grid.trirllrho, grid.mask.flatten())
    mask0 = fmask(lon0, lat0)  # mask for lon0/lat0 points
    ind1 = (mask0 == 1.)  # indices select out where points are masked

    # If nobays, eliminate points with shallow bathymetry
    if nobays:
        fh = grid.trirllrho.nn_interpolator(grid.h.flatten())
        h0 = fh(lon0, lat0)
        ind2 = (h0 > 10.)
    else:
        ind2 = np.ones(ind1.shape).astype(bool)

    ind2 = ~np.isnan(lon0)*ind1*ind2

    L = len(ind2.ravel())
    Lnan = sum(ind2.ravel())
    print L - Lnan, '/', L, ' drifters NaN-ed out.'

    lon0 = lon0[ind2].flatten()
    lat0 = lat0[ind2].flatten()

    if z0 is not None:
        z0 = z0[ind2].flatten()

    if z0 is None:
        return lon0, lat0
    else:
        return lon0, lat0, z0


def seed(lon, lat, dlon=.5, dlat=.5, N=30):
    """
    Chose array of starting locations based on the start location. A Gaussian
    distribution is used to distribute points around the find location.

    Args:
        lon, lat: Start location
        dlon, dlat: Distance in degrees in which to seed drifters. Default is
         0.5 degrees.
        N: Number of drifters in x and y. Default is 30.

    Returns:
        * lon0, lat0 - Points in lon/lat at which to seed drifters
    """

    # Find 2D distribution of points around the package location
    # the center is indicated using (lon, lat)
    # The variance is given by [[.25,0],[0,.25]] indicates
    #  a standard deviation away from the center of .5 degrees
    #  or .5**2=.25
    # There are N points in both x and y
    dist = np.random.multivariate_normal((lon, lat),
                                         [[dlon**2, 0], [0, dlat**2]],
                                         [N, N])
    return dist[:, :, 0], dist[:, :, 1]


def make_proj(setup='nwgom', usebasemap=True, **kwargs):
    """Convenience function for setting up a projection object.

    This fills in some of the work of setting up a projection to use with a
    horizontal grid. User can choose between default setups and then change
    keyword arguments within the setups.

    Args:
        setup (str): Which default setup to use.

            'nwgom' - for NW Gulf of Mexico, for use with basemap
            'galveston' - for Galveston Bay, for use with pyproj
            'nwgom-pyproj' - for NW Gulf of Mexico, for use without basemap

        usebasemap (bool): True is use basemap for your projection, with
            specific keyword arguments. False is to use pyproj for your
            projection.

    Kwargs:
        usebasemap (bool): Whether to use load basemap into grid (True) or
            pyproj (False). Basemap is slower but can be used for plotting,
            and pyproj is the opposite. Default is False.
        llcrnrlon (float): Lower left corner longitude, in degrees. This and
            llcrnrlat, urcrnrlon, urcrnrlat, lat_0, lon_0, and res are only
            used if usebasemap=True.
        llcrnrlat (float): Lower left corner latitude.
        urcrnrlon (float): Upper right corner longitude.
        urcrnrlat (float): Upper right corner latitude.
        lat_0 (float): Latitude of center of projection.
        lon_0 (float): Longitude of center of projection.
        res (str): Resolution of coastal map used in basemap.
        zone (int): Zone for utm projection

    Returns:
        proj: Projection object for convering between geometric and projected
            coordinates.
    """

    # Default setups. Can use all of these parameters or can start with them
    # and replace some individually.
    if setup == 'nwgom':

        inputs = {'projection': 'lcc', 'llcrnrlon': -98.5, 'llcrnrlat': 22.5,
                  'urcrnrlon': -87.5, 'urcrnrlat': 31.0, 'lat_0': 30,
                  'lon_0': -94, 'resolution': 'i', 'area_thresh': 0.}

        usebasemap = True

    elif setup == 'galveston':

        inputs = {'proj': 'utm', 'zone': 15, 'ellps': 'clrk66',
                  'datum': 'NAD27'}

        usebasemap = False

    elif setup == 'nwgom-pyproj':

        inputs = {'proj': 'lcc', 'ellps': 'clrk66', 'datum': 'NAD27',
                  'lat_1': 22.5, 'lat_2': 31.0, 'lat_0': 30, 'lon_0': -94,
                  'x_0': 0, 'y_0': 0}

        usebasemap = False

    else:

        inputs = kwargs

    # If there are any input keyword arguments, use them to replace setup
    # arguments in inputs dictionary.
    if kwargs:
        for key, value in kwargs.items():
            if key in inputs:  # if key is in inputs, replace it
                inputs[key] = value
            else:  # otherwise, error
                print 'input key is not in inputs dict'

    # Set up projection using inputs dictionary
    if usebasemap:
        from mpl_toolkits.basemap import Basemap
        proj = Basemap(**inputs)
    else:
        from pyproj import Proj
        proj = Proj(**inputs)

    return proj
