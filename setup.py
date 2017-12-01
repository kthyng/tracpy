#!/usr/bin/env python

"""
setup.py for tracpy -- python wrappers around TRACMASS

"""
import shutil
from setuptools import setup # to support "develop" mode
from numpy.distutils.core import setup, Extension

# ## total kludge:
# ## this copies extension built with makefile to package
# ## better to have setup.py build TRACMASS

# try:
#     shutil.copyfile("src/tracmass.so", "tracpy/tracmass.so")

# except IOError:
#     raise Exception("TRAMASS extesnion not built -- built it with the makefile in src")

tracmass_mod = Extension(name = "tracmass",
                         sources=['src/calc_dxyz.f95',
                                  'src/cross.f95',
                                  'src/loop_pos.f95',
                                  'src/step.f95',
                                  'src/vertvel.f95',
                                  'src/calc_time.f95',
                                  'src/diffusion.f95',
                                  'src/pos.f95',
                                  'src/turb.f95',
                                  ],
                      #extra_f90_compile_args=["-ffixed-form"],
                      )

# modules = ["tracpy/inout.py",
#            "tracpy/op.py"
#            "tracpy/plotting.py",
#            "tracpy/run.py",
#            "tracpy/tools.py",
#            ]

print(tracmass_mod)

setup(
    name = "tracpy",
    version = "1.0",
    author = "Kristen Thyng",
    author_email = "",
    description = ("python wrappers around TRACMASS"),
    long_description=open('README.md').read(),
    classifiers=[
                 "Development Status :: 3 - Alpha",
    #             "Topic :: Utilities",
                 ],
    packages = ["tracpy"],
    # py_modules = modules,
    ext_package='tracpy',
    ext_modules = [tracmass_mod],
    scripts = [],
    )
