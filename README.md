# SolTrack #

![PyPI](https://img.shields.io/pypi/v/soltrack?color=%230A0)
![PyPI - Downloads](https://img.shields.io/pypi/dm/soltrack)
[![Documentation
Status](https://readthedocs.org/projects/soltrack/badge/?version=latest)](https://soltrack.readthedocs.io/en/latest/?badge=latest)
![PyPI - Licence](https://img.shields.io/pypi/l/soltrack?color=%230A0)

A free, fast and simple Python package to compute the position of the Sun, as well as its rise and set times.
SolTrack was originally written in C/C++ by [Marc van der Sluys](http://marc.vandersluys.nl) of the
department of astrophysics of the Radboud University Nijmegen, the Netherlands and the Sustainable energy
research group of the HAN University of Applied Sciences in Arnhem, the Netherlands, and Paul van Kan of the
Sustainable energy research group of the HAN University of Applied Sciences in Arnhem, the Netherlands.  The
code has now been translated to pure Python.

SolTrack performs about 12340 position calculations per second on a single 3.6 GHz core of my 2013 laptop,
with an accuracy of 0.0030 ± 0.0016°.  This makes it about 500x times faster than astropy, but around 50x
slower than pyEphem, which is written in C.  SolTrack has been used to control solar trackers, as well as for
modelling in solar energy.


## Installation ##

This package can be installed using `pip install soltrack`.  This should automatically install the dependency
packages `numpy`, `pytz` and `dataclasses` (if you're not using Python 3.7+) if they haven't been installed
already.  If you are installing by hand, ensure that these packages are installed as well.


## Example use ##

## Example for number of date/times ##

```python
"""Compute the position of the Sun and its rise and set times for a vector of instants."""

# Create a Pandas DatetimeIndex range every 20 days 1hour and 10 minutes, in my timezone:
import pandas as pd
dti = pd.date_range('2022-01-01 01:00', '2022-12-31 23:00', freq='20d1h10min', tz='Europe/Amsterdam')

# Set the geographic location to Arnhem, the Netherlands (we'll use degrees in SolTrack):
geo_lon =  5.950270  # Positive -> east of Greenwich (degrees)
geo_lat = 51.987380  # Positive -> northern hemisphere (degrees)


# Create a SolTrack instance and specify preferences:
from soltrack import SolTrack
st = SolTrack(geo_lon, geo_lat, use_degrees=True)  # Use default values for all but use_degrees
st.set_date_time(dti)  # Pass my dates and times to SolTrack
st.compute_position()  # Compute the Sun's position
st.compute_rise_set()  # Compute the rise and set times of the Sun


# Print some selected results as arrays and create chaos:
if st.lt is not None:  # If local time was used
    print('Local time:     ', *st.lt)  # Asterisk (*) unpacks the DTI
    
print('UTC:            ', *st.utc)
print('azimuth:        ', *st.azimuth)
print('altitude:       ', *st.altitude)
print('distance:       ', *st.distance)
print('riseTime:       ', *st.rise_time)
print('transTime:      ', *st.transit_time)
print('setTime:        ', *st.set_time)


# Store selected results in a Pandas DataFrame and print that in a more orderly fashion:
st.create_df(utc=True, jd=True, ecl=True, rts_pos=True)
with pd.option_context('display.max_columns',None, 'display.width',None):  # Want all columns
    print(st.df)
```

### Example use for a single instant ###
```python
```




## SolTrack pages ##

* [SourceForge](http://soltrack.sf.net): SolTrack homepage
* [Pypi](https://pypi.org/project/soltrack/): SolTrack Python package
* [GitHub](https://github.com/MarcvdSluys/SolTrack-Python): SolTrack Python source code
* [ReadTheDocs](https://soltrack.readthedocs.io/en/latest/): SolTrack Python documentation


## Author and licence ##

* Author: Marc van der Sluys
* Contact: http://marc.vandersluys.nl
* Licence: [EUPL 1.2](https://www.eupl.eu/1.2/en/)


## References ##

* [Celestial mechanics in a nutshell (CMiaNS)](https://cmians.sourceforge.io/)
* Meeus, [Astronomical algorithms](https://www.willbell.com/math/MC1.HTM), 2nd Ed.
* The C and Python codes have been adapted from the Fortran implementation in
  [libTheSky](http://libthesky.sourceforge.net/), which contains many references.
