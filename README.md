#tracpy

Fortran core of Tracmass + Python wrapping around the outside.

To get the Fortran Tracmass code ready:

1. `make clean`
1. Compile tracmass code in Fortran: `make f2py`

To learn about the module Tracpy:

1. Open an iPython notebook server in the tracpy directory: `ipython notebook`
1. A webpage will open in your browser showing the available notebooks in the current directory. Open the notebook called manual
1. The cells of the notebook can be run in order by pushing "shift" and "enter" together or the whole notebook can be run by selecting Cell > Run all. The notebook demonstrates how to initialize and run a numerical drifter experiment using tracpy.
1. Alternatively, a static PDF version of the manual has been saved and can be viewed (manual.pdf) but not run.

For projects, it is suggested to start a separate directory with its own initialization and run file for the simulation(s). Then tracpy can be imported as a module and the run script can be run in the background from the terminal with standard output being redirected to a log file using `python run.py > log.txt &`