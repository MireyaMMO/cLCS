#!/usr/bin/env python

import setuptools

setuptools.setup(
    name        = 'cLCS',
    version     = '1',
    author='Rodrigo Duran adapted to python by Mireya Montaño',
    packages=['cLCS'],
    license='',
    url='',
    description='Climatological LCS code for python using OpenDrift as seen in Rodrigo Duran code for Matlab',
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
)
