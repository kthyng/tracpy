#tracpy

Fortran core of Tracmass + Python wrapping around the outside.

To get the Fortran Tracmass code ready:

1. `make clean`
1. Compile tracmass code in Fortran: `make f2py`

Note: if the code will not compile, a first step could be to change the flag in the makefile from -m64 to -m32 if you are set up to use 32 bit instead of 64 bit.

To learn about the module Tracpy:

1. Open an iPython notebook server in the tracpy directory: `ipython notebook`
1. A webpage will open in your browser showing the available notebooks in the current directory. Open the notebook called manual
1. The cells of the notebook can be run in order by pushing "shift" and "enter" together or the whole notebook can be run by selecting Cell > Run all. The notebook demonstrates how to initialize and run a numerical drifter experiment using tracpy.
1. Alternatively, a static PDF version of the manual has been saved and can be viewed (manual.pdf) but not run.

For projects, it is suggested to start a separate directory with its own initialization and run file for the simulation(s). Then tracpy can be imported as a module and the run script can be run from the project. An example project can be found at https://github.com/kthyng/gisr.git.

Some more information about running on Linux machines with taskset to control what cores a given process is run on (the syntax of these commands depends on your run.py script):
Multiple instances of the simulation can be run on different cores if you have a multi-core Linux machine. To do this, use `taskset`:
`taskset 3 python2.7 run.py > temp_log.txt &`
This command runs the run.py script for a specific project in the background (due to the &), redirects the output from the screen to the text file 'temp_log.txt', and runs the process on core 3. This command can then be used for other instances by choosing a different core after the command `taskset`.

If a process is already running, you can check its current core using its PID:
`taskset -p [PID]`
Then you can move the process to a different core using
`taskset -p [CORE NUMBER] [PID]`