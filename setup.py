#!/bin/env python3
# -*- coding: utf-8 -*-

"""Setup.py for the SolTrack Python package."""

# Package version:
version='0.2.0'

# Get long description from README.md:
with open('README.md', 'r') as fh:
    long_description = fh.read()


# Set package properties:
from setuptools import setup
setup(
    name='soltrack',
    description='A free, fast and accurate Python package to compute the position of the Sun',
    author='Marc van der Sluys',
    url='http://soltrack.sf.net',
    
    packages=['soltrack'],
    install_requires=['astrotool','numpy','pandas'],
    long_description=long_description,
    long_description_content_type='text/markdown',
    
    version=version,
    license='EUPL 1.2',
    keywords=['astronomy','ephemeris','sun','solar','solar energy'],
    
    classifiers=[
        'Development Status :: 4 - Beta',
        'Intended Audience :: Developers',
        'Intended Audience :: Education',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: European Union Public Licence 1.2 (EUPL 1.2)',
        'Natural Language :: English',
        'Operating System :: OS Independent',
        'Programming Language :: Python :: 3',
        'Topic :: Scientific/Engineering',
        'Topic :: Scientific/Engineering :: Astronomy',
        'Topic :: Scientific/Engineering :: Physics',
    ]
)
