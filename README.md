# TracPy

[<img src="https://zenodo.org/badge/4563/kthyng/tracpy.png" class="picFloat">](https://zenodo.org/record/10433#.U6SWWBYxmd8)

Fortran core of TRACMASS + Python wrapping around the outside.


## To get the code

1. Make your new TracPy directory and change directories into it.
1. Clone the TracPy repository from GitHub. 
In the command prompt, type: 
`git clone https://github.com/kthyng/tracpy.git`
1. Install the package globally:
`pip install -e .`
This makes the package an editable install so that it can be updated with future additions to TracPy. Note that a required package is [octant](https://github.com/hetland/octant). To instead install the package locally: 
`pip install --user .`


## To update the code later

1. Move into your TracPy directory.
1. Update your GitHub repository.
`git pull`
1. Edit your install of TracPy.
`pip install -e .` 
or
`pip install --force-reinstall -e .`
or, for local installation: 
`pip install --ignore-installed --user .`


## To test the code

1. After making changes to the code, you can do some basic functionality testing with `py.test` or `nosetests` in the `tests` subdirectory.


## To learn more about the module TracPy

Learn more by running a small test case. Internet required.

1. Open an iPython notebook server in the TracPy `docs` directory.
`ipython notebook`
1. A webpage will open in your browser showing the available notebooks in the current directory. Open the notebook called manual.
1. The cells of the notebook can be run in order by pushing "shift" and "enter" together or the whole notebook can be run by selecting Cell > Run all. The notebook demonstrates how to initialize and run a numerical drifter experiment using TracPy.
1. Alternatively, a static PDF version of the manual can be viewed at `http://nbviewer.ipython.org/urls/raw.github.com/kthyng/tracpy/master/docs/manual.ipynb`.


## To run your own projects

For projects, it is suggested to start a separate directory with its own initialization and run file for the simulation(s). Then TracPy can be imported as a module and the run script can be run from the project. An example set of projects can be found at `https://github.com/kthyng/gisr.git`.


## To learn more about TRACMASS

* Döös, K., Kjellsson, J., & Jönsson, B. (2013). TRACMASS—A Lagrangian Trajectory Model. In Preventive Methods for Coastal Protection (pp. 225-249). Springer International Publishing.
* Döös, K., Rupolo, V., & Brodeau, L. (2011). Dispersion of surface drifters and model-simulated trajectories. Ocean Modelling, 39(3), 301-310.
* Döös, K., & Engqvist, A. (2007). Assessment of water exchange between a discharge region and the open sea–A comparison of different methodological concepts. Estuarine, Coastal and Shelf Science, 74(4), 709-721.
* TRACMASS on GitHub: https://github.com/TRACMASS
