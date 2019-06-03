#!/usr/bin/env python
from __future__ import absolute_import
#from distutils.core import setup,Extension
from setuptools import setup
import sys

setup(
  name='PolyXSim',
  version='1.1.1',
  description='Simulation of diffraction images',
  license='GPL', 
  maintainer='Henning Osholm Soerensen and Jette Oddershede',
  maintainer_email='henning.sorensen@risoe.dk',
  url='http://github.com/FABLE-3DXRD/PolyXSim',
  packages=['polyxsim'],
  package_dir={"polyxsim": "polyxsim"},
  scripts=["scripts/PolyXSim.py"]
)
