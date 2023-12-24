#!/usr/bin/env python

import os
import setuptools

here = os.path.abspath(os.path.dirname(__file__))
#exec(open(os.path.join(here)).read())

setuptools.setup(
    name        = 'cLCS',
    description = 'Climatological LCS code for python using OpenDrift as seen in Rodrigo Duran code for Matlab ',
    author      = 'Rodrigo Duran adapted to python by Mireya Montano',
    url         = 'https://github.com/MireyaMMO/cLCS',
    download_url = 'https://github.com/MireyaMMO/cLCS',
    version = 1,
    license = '',
    packages=['cLCS'],
#    long_description=open('README.rst').read(),
    install_requires=[
       'numpy',
       'scipy',
       'matplotlib',
       'netCDF4',
       'pyproj',
       'cartopy',
       'opendrift-landmask-data'
       'calendar'
       'opendrift'
       'pickle'
       'datetime'
    ],
    include_package_data = True
)
