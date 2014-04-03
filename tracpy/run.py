import numpy as np
import op
import netCDF4 as netCDF
import pdb
from tracpy.time_class import Time

def run(tp, date, lon0, lat0):
    '''
    some variables are not specifically called because f2py is hides them
     like imt, jmt, km, ntractot
    Look at tracmass.step to see what it is doing and making optional at the end.

    Inputs:

    tp          TracPy object, from the Tracpy class.
    date        Start date in datetime object
    lon0        Drifter starting locations in x/zonal direction.
    lat0        Drifter starting locations in y/meridional direction.

    Other variables:

    xp          x-locations in x,y coordinates for drifters
    yp          y-locations in x,y coordinates for drifters
    zp          z-locations (depths from mean sea level) for drifters
    t           time for drifter tracks
    '''

    timer = Time() # start timer for simulation

    # Initialize everything for a simulation
    tinds, nc, t0save, xend, yend, zend, zp, ttend, flag = tp.prepare_for_model_run(date, lon0, lat0)

    timer.addtime('1: Preparing for simulation   ')

    # Loop through model outputs.
    for j,tind in enumerate(tinds[:-1]):

        print 'Using GCM model output index ', j

        # Loop through substeps in call to TRACMASS in case we want to add on windage, etc, for each step
        for nsubstep in xrange(tp.nsubsteps):

            xstart, ystart, zstart, ufsub, vfsub = tp.prepare_for_model_step(tinds[j+1], nc, flag, xend, yend, zend, j, nsubstep)
            ind = (flag[:] == 0) # indices where the drifters are still inside the domain

            timer.addtime('2: Preparing for model step   ')

            if not np.ma.compressed(xstart).any(): # exit if all of the drifters have exited the domain
                break

            # Do stepping in Tracpy class
            xend_temp,\
                yend_temp,\
                zend_temp,\
                flag[ind],\
                ttend_temp, U, V = tp.step(xstart, ystart, zstart, ufsub, vfsub)

            timer.addtime('3: Stepping, using TRACMASS   ')

            xend[ind,j*tp.N+1:j*tp.N+tp.N+1], \
                yend[ind,j*tp.N+1:j*tp.N+tp.N+1], \
                zend[ind,j*tp.N+1:j*tp.N+tp.N+1], \
                zp[ind,j*tp.N+1:j*tp.N+tp.N+1], \
                ttend[ind,j*tp.N+1:j*tp.N+tp.N+1] = tp.model_step_is_done(xend_temp, yend_temp, zend_temp, ttend_temp, ttend[ind,j*tp.N])

            timer.addtime('4: Processing after model step')

    nc.close()

    lonp, latp, zp, ttend, grid, T0, U, V = tp.finishSimulation(ttend, t0save, xend, yend, zp)

    timer.addtime('5: Processing after simulation')

    print "============================================="
    print ""
    print "Simulation name: ", tp.name
    print ""
    print "============================================="

    timer.write()

    return lonp, latp, zp, ttend, grid, T0, U, V
