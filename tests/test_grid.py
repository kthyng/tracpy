#!/usr/bin/env python

"""
Test projection and grid routines. Generally test all available pyproj
projection presets but not basemap since they are slower.
"""

import tracpy
import tracpy.calcs
import os
import numpy as np
import matplotlib.tri as mtri


# List projection setup for use in tests
# pyproj-based setups
projpyproj = ['galveston', 'nwgom-pyproj']
projbasemap = ['nwgom']  # don't usually test with this since slow


def test_proj_init():
    """Test initialization of preset pyproj projections."""

    # loop through projection presets
    for projsetup in projpyproj:

        # Get projection object
        proj = tracpy.tools.make_proj(setup=projsetup)

    assert proj


def test_grid_init():
    """Test initialization of grid."""

    # loop through projection presets
    for projsetup in projpyproj:

        # Get projection object
        proj = tracpy.tools.make_proj(setup=projsetup)

        grid_filename = os.path.join('input', 'grid.nc')

        # Read in grid
        grid = tracpy.inout.readgrid(grid_filename, proj, usespherical=True)

        assert grid


def test_proj_variant():
    """Test creating a projection with different than built-in variables."""

    pass


def test_proj_iteration():
    """Test projection conversion back and forth between spaces.

    Set up a projection, then convert between spaces and check that the
    result is close to the starting values.
    """

    # loop through projection presets
    for projsetup in projpyproj:

        # Get projection object. Can use either 'galveston' or 'nwgom-pyproj'
        # built in projection setups to test quickly ('nwgom' is for use with
        # usebasemap=True and thus is slow for testing).
        proj = tracpy.tools.make_proj(setup=projsetup)

        grid_filename = os.path.join('input', 'grid.nc')

        # Read in grid
        grid = tracpy.inout.readgrid(grid_filename, proj, usespherical=True)

        # convert back and forth
        lon_rho2, lat_rho2 = grid.proj(grid.x_rho, grid.y_rho, inverse=True)

        print(grid.lat_rho[0, :])
        print(lat_rho2[0, :])
        print(grid.lon_rho[0, :])
        print(lon_rho2[0, :])

        assert np.allclose(grid.lat_rho, lat_rho2)
        assert np.allclose(grid.lon_rho, lon_rho2)


def test_grid_triangulation_spherical():
    """Test that the grid triangulations are valid: spherical test cases."""

    # loop through projection presets
    for projsetup in projpyproj:

        # Get projection object
        proj = tracpy.tools.make_proj(setup=projsetup)

        grid_filename = os.path.join('input', 'grid.nc')

        # Read in grid
        grid = tracpy.inout.readgrid(grid_filename, proj, usespherical=True)

        assert mtri.LinearTriInterpolator(grid.trir, grid.x_rho.flatten())


def test_grid_triangulation_projected():
    """Test that the grid triangulations are valid: projected test cases."""

    # loop through projection presets
    for projsetup in projpyproj:

        # Get projection object
        proj = tracpy.tools.make_proj(setup=projsetup)

        grid_filename = os.path.join('input', 'gridxy.nc')

        # Read in grid
        grid = tracpy.inout.readgrid(grid_filename, proj, usespherical=False)

        assert mtri.LinearTriInterpolator(grid.trir, grid.x_rho.flatten())


def test_interpolation():
    """Test interpolation with grid space and projected grid the same.

    Create a test case with the 'projected' grid in grid space coordinates.
    When interpolating between them, there should be a shift because the
    rho points in projected space are not in the same setup as grid coords.
    """

    # Get projection object
    proj = tracpy.tools.make_proj(setup='nwgom-pyproj')

    grid_filename = os.path.join('input', 'gridij.nc')

    # Read in grid
    grid = tracpy.inout.readgrid(grid_filename, proj, usespherical=False)

    # Do some interpolating
    # projected grid to grid space, delaunay
    X, Y, _ = tracpy.tools.interpolate2d(grid.x_rho[2, 3], grid.y_rho[2, 3],
                                         grid, 'd_xy2ij')

    # There is a shift between the rho grid and the grid space grid because
    # of the staggered layout. Grid space counts from the u/v grid and
    # therefore is a little different from the rho grid.
    assert np.allclose(X, grid.x_rho[2, 3] - 0.5)
    assert np.allclose(Y, grid.y_rho[2, 3] - 0.5)

    # grid space to projected coordinates, delaunay
    x, y, _ = tracpy.tools.interpolate2d(grid.X[2, 3], grid.Y[2, 3], grid,
                                         'd_ij2xy')

    assert np.allclose(x, grid.X[2, 3] + 0.5)
    assert np.allclose(y, grid.Y[2, 3] + 0.5)

    # grid space to projected coordinates, map_coords
    x, y, _ = tracpy.tools.interpolate2d(grid.X[2, 3], grid.Y[2, 3], grid,
                                         'm_ij2xy')

    assert np.allclose(x, grid.X[2, 3] + 0.5)
    assert np.allclose(y, grid.Y[2, 3] + 0.5)


# def test_interpolation():
#     """Test that interpolation methods work.

#     Convert back and forth between spaces using multiple tools.interpolation2d
#     methods to make sure the values stay close.
#     """

#     # Get projection object
#     proj = tracpy.tools.make_proj(setup='nwgom-pyproj')

#     grid_filename = os.path.join('input', 'grid.nc')

#     # Read in grid
#     grid = tracpy.inout.readgrid(grid_filename, proj, usespherical=True)

#     # Do some interpolating
#     X_rho, Y_rho, _ = tracpy.tools.interpolate2d(grid.lon_rho[2, 3],
#                                                  grid.lat_rho[2,3], grid, 'd_ll2ij')

#     print X_rho
#     print grid.X[2,3]

#     assert np.allclose(X_rho, grid.X[2,3]-0.5)


def test_verts():
    """Test properties of vertices."""

    pass
