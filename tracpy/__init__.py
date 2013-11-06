'''
tracpy is a wrapper written in Python for using Tracmass, a particle 
tracking system written in Fortran. Tracmass was written by other
researchers and is available online.

See README.md for an overview on instructions and the iPython 
notebook manual.ipynb for instructions on use. A user manual for 
Tracmass is available online.

Modules available in tracpy include:

* init.py
* inout.py
* manual.ipynb
* op.py
* plotting.py
* run.py
* tools.py

Modules available in tracmass include:

* calc_dxyz.f95
* calc_time.f95
* cross.f95
* diffusion.f95
* loop_pos.f95
* makefile
* pos.f95
* step.f95
* turb.f95
* vertvel.f95
'''

import init
import inout
import op
import plotting
import run
import tools

__authors__ = ['Kristen Thyng <kthyng@tamu.edu>']
