#!/usr/bin/env python
from distutils.core import setup,Extension
import sys

setup(
  name='simul_farfield',
  version='0.0.1',
  description='Image viewer for file series of diffraction images',
  license='GPL', maintainer='Henning Osholm Soerensen',
  maintainer_email='henning.sorensen@risoe.dk',
  url='http://fable.wiki.sourceforge.net',
  packages=['PolyXSim'],
  package_dir={"PolyXSim": "src"},
  scripts=["scripts/PolyXSim.py"]
)
