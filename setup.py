#!/bin/env python3

"""Setup.py for the SolTrack Python package."""

# Package version:
version="0.1.1"

# Get long description from README.md:
with open("README.md", "r") as fh:
    long_description = fh.read()


# Set package properties:
from setuptools import setup
setup(
    name='soltrack',
    description='A free, fast and accurate Python package to compute the position of the Sun',
    author='Marc van der Sluys',
    url='http://soltrack.sf.net',
    
    packages=['soltrack'],
    install_requires=['numpy','dataclasses','pytz'],
    long_description=long_description,
    long_description_content_type='text/markdown',
    
    version=version,
    license='GPLv3+',
    keywords=['astronomy','ephemeris','sun','solar','solar energy'],
    
    classifiers=[
        "Development Status :: 4 - Beta",
        "Intended Audience :: Developers",
        "Intended Audience :: Education",
        "Intended Audience :: Science/Research",
        "License :: OSI Approved :: GNU General Public License v3 or later (GPLv3+)",
        "Natural Language :: English",
        "Operating System :: OS Independent",
        "Programming Language :: Python :: 3",
        "Topic :: Scientific/Engineering",
        "Topic :: Scientific/Engineering :: Astronomy",
        "Topic :: Scientific/Engineering :: Physics",
    ]
)

