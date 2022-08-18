#!/usr/bin/env python

import os
import setuptools

here = os.path.abspath(os.path.dirname(__file__))
exec(open(os.path.join(here, 'opendrift/version.py')).read())

setuptools.setup(
    name        = 'cLCS',
    description = 'Climatological LCS code for python using OpenDrift as seen in Rodrigo Duran code for Matlab ',
    author      = 'Rodrigo Duran / Mireya Montaño',
    url         = 'https://github.com/MireyaMMO/cLCS',
    download_url = 'https://github.com/MireyaMMO/cLCS',
    version = __version__,
    license = '',
    install_requires = [
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
    packages = setuptools.find_packages(),
    include_package_data = True
)
