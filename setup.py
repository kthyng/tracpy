#!/usr/bin/env python

"""
setup.py for tracpy -- python wrappers around TRACMASS

"""
import shutil
from distutils.core import Extension
from setuptools import setup # to support "develop" mode

## total kludge:
## this copies extension built with makefile to package
## better to have setup.py build TRACMASS

try:
    shutil.copyfile("src/tracmass.so", "tracpy/tracmass.so")

except IOError:
    raise Exception("TRAMASS extesnion not built -- built it with the makefile in src")

tracmass = Extension("tracpy/tracmass",
                     sources=["src/tracmass.so"],
                     )

modules = ["tracpy/inout.py",
           "tracpy/op.py"
           "tracpy/plotting.py",
           "tracpy/run.py",
           "tracpy/tools.py",
           "tracpy/tracmass.so"
           ]

setup(
    name = "tracpy",
    version = "0.01",
    author = "Kristen Thyng",
    author_email = "",
    description = ("python wrappers around TRACMASS"),
    long_description=open('README.md').read(),
    classifiers=[
                 "Development Status :: 3 - Alpha",
    #             "Topic :: Utilities",
                 ],
    packages = ["tracpy"],
    modules = modules, 
    ext_modules = [tracmass],
    scripts = [],
    )
