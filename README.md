tracpy
======

Fortran core of Tracmass + Python wrapping around the outside.

0. make clean
1. Compile tracmass code in Fortran: make f2py
2. Import the Python wrapper: import run
3. Run the Python wrapper in ipython: xp,yp,zp,t=run.start_run()
   Make sure to have the initialization you want to run set up to run in run.start_run()
