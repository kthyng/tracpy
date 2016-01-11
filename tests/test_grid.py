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

        print grid.lat_rho[0, :]
        print lat_rho2[0, :]
        print grid.lon_rho[0, :]
        print lon_rho2[0, :]

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

    # Get projection object
    proj = tracpy.tools.make_proj(setup='nwgom-pyproj')

    grid_filename = os.path.join('input', 'gridxy.nc')

    # Read in grid
    grid = tracpy.inout.readgrid(grid_filename, proj, usespherical=False)

    assert mtri.LinearTriInterpolator(grid.trir, grid.x_rho.flatten())


def test_interpolation():
    """Test that interpolation methods work.

    Convert back and forth between spaces using multiple tools.interpolation2d
    methods to make sure the values stay close.
    """

    pass
