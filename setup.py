# -*- coding: utf-8 -*-
"""
Created on Thu Aug  7 11:43:06 2014

@author: jmmauricio
"""

import os
from setuptools import setup, find_packages

# Utility function to read the README file.
# Used for the long_description.  It's nice, because now 1) we have a top level
# README file and 2) it's easier to type in the README file than to put a raw
# string in below ...
def read(fname):
    return open(os.path.join(os.path.dirname(__file__), fname)).read()

setup(
    name = "pysimu",
    version = "0.0.1c",
    author = "Juan Manuel Mauricio",
    author_email = "jmmauricio6@gmail.com",
    description = ("Simulations tools"),
    long_description = """\
This package provides a module pysimu, that can be used to simulate dynamical systems of the form
..
$\hat x = f(t,x)$

.. math::

   \frac{ \sum_{t=0}^{N}f(t,k) }{N}

Please see README.txt for examples.
""",
    license = "BSD",
    keywords = "dynamics model simulation",
    url = "http://packages.python.org/pysimu",
    packages=find_packages(where='.'),
    classifiers=[
        "Development Status :: 3 - Alpha",
        "Topic :: Utilities",
        "License :: OSI Approved :: BSD License",
    ],
)