'''
Main run script for TracPy system.
'''


import numpy as np
from tracpy.time_class import Time


def run(tp, date, lon0, lat0, T0=None, U=None, V=None):
    """
    some variables are not specifically called because f2py is hides
    them like imt, jmt, km, ntractot.

    Look at tracmass.step to see what it is doing and making optional at the
    end.

    Args:
        tp: TracPy object, from the Tracpy class.
        date (datetime object): Start date.
        lon0 (array): Drifter starting locations in x/zonal direction.
        lat0 (array): Drifter starting locations in y/meridional direction.
        T0 (Optional[array]): Weighting of drifters for use with stream
         functions. Is not used if dostream=0.
        U,V (Optional[array]): For east-west/north-south transport, is
         updated by TRACMASS. Only used if dostream=1.

    Other variables:
     * xp: x-locations in x,y coordinates for drifters
     * yp: y-locations in x,y coordinates for drifters
     * zp: z-locations (depths from mean sea level) for drifters
     * t: time for drifter tracks
    """

    timer = Time()  # start timer for simulation

    # Initialize everything for a simulation
    tinds, nc, t0save, xend, yend, \
        zend, zp, ttend, flag = tp.prepare_for_model_run(date, lon0, lat0)

    timer.addtime('1: Preparing for simulation   ')

    # Loop through model outputs.
    for j, tind in enumerate(tinds[:-1]):

        print 'Using GCM model output index ', j

        # Loop through substeps in call to TRACMASS in case we want to add on
        # windage, etc, for each step
        for nsubstep in xrange(tp.nsubsteps):

            xstart, ystart, zstart, ufsub, vfsub, T0 = \
                tp.prepare_for_model_step(tinds[j+1], nc, flag, xend, yend,
                                          zend, j, nsubstep, T0)
            # indices where the drifters are still inside the domain
            ind = (flag[:] == 0)

            timer.addtime('2: Preparing for model step   ')

            # exit if all of the drifters have exited the domain
            if not np.ma.compressed(xstart).any():
                break

            # Do stepping in Tracpy class
            xend_temp,\
                yend_temp,\
                zend_temp,\
                flag[ind],\
                ttend_temp, U, V = tp.step(xstart, ystart, zstart, ufsub,
                                           vfsub, T0, U, V)

            timer.addtime('3: Stepping, using TRACMASS   ')

            xend[ind, j*tp.N+1:j*tp.N+tp.N+1], \
                yend[ind, j*tp.N+1:j*tp.N+tp.N+1], \
                zend[ind, j*tp.N+1:j*tp.N+tp.N+1], \
                zp[ind, j*tp.N+1:j*tp.N+tp.N+1], \
                ttend[ind, j*tp.N+1:j*tp.N+tp.N+1] = \
                tp.model_step_is_done(xend_temp, yend_temp, zend_temp,
                                      ttend_temp, ttend[ind, j*tp.N])

            timer.addtime('4: Processing after model step')

    nc.close()

    lonp, latp, zp, ttend, T0, U, V = tp.finishSimulation(ttend, t0save,
                                                          xend, yend, zp, T0,
                                                          U, V)

    timer.addtime('5: Processing after simulation')

    print "============================================="
    print ""
    print "Simulation name: ", tp.name
    print ""
    print "============================================="

    timer.write()

    return lonp, latp, zp, ttend, T0, U, V
