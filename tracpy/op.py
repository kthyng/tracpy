# -*- coding: utf-8 -*-
"""
Created on Tue Sep  4 13:42:38 2012

@author: kthyng

Various operations
"""

import numpy as np


def resize(A, dim):
    """
    Average neighboring elements in an array A along a dimension dim.

    Args:
        A (array): array size, [m x n] in 2D case. Can be up to 3D.
        dim (int): dimension on which to act. Can be up to 2 (0, 1, 2).

    Returns:
        * A - array of size [(m-1) x n] if dim is 0
    """

    # B is A but rolled such that the dimension that is to be resized is in
    # the 0 position
    B = np.rollaxis(A, dim)

    # Do averaging
    B = 0.5*(B[0:-1]+B[1:])

    # Roll back to original
    return np.rollaxis(B, 0, dim+1)


def rotate(ucurv, vcurv, angle):
    """
    Rotate u and v on curvilinear grid about angle, in radians
    counterclockwise from east, to be on cartesian grid.

    Args:
        ucurv: u on a curvilinear grid, can be array
        vcurv: v on a curvilinear grid, can be array of same size as ucurv
        angle: angle between curvilinear grid and true east, in radians,
         can be array of same size as ucurv

    Returns:
        * u - true east u velocity
        * v - true north v velocity
    """

    u = ucurv*np.cos(angle) - vcurv*np.sin(angle)
    v = ucurv*np.sin(angle) + vcurv*np.cos(angle)
    return u, v


def find_nearest_index(xr, yr, x0, y0):
    """
    Find the index ind in vec of point pt.

    Args:
        xr (array): 2d array of x locations
        yr (array): 2d array of y locations
        x0: x-coord of Point that may be in Vector
        y0: y-coord of Point that may be in Vector

    Returns:
        * ind - index of pt in vec
    """

    d = np.sqrt((xr-x0)**2 + (yr-y0)**2)

    # index in flattened array of closest point/smallest distance
    ind = d.argmin()
    J, I = np.unravel_index(ind, d.shape)  # index in 2d array
    return J, I
