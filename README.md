# SolTrack for Python

![PyPI](https://img.shields.io/pypi/v/soltrack?color=%230A0)
![PyPI - Downloads](https://img.shields.io/pypi/dm/soltrack)
[![Documentation
Status](https://readthedocs.org/projects/soltrack/badge/?version=latest)](https://soltrack.readthedocs.io/en/latest/?badge=latest)
![PyPI - Licence](https://img.shields.io/pypi/l/soltrack?color=%230A0)

A free, fast and simple Python package to compute the position of the Sun, as well as its rise and set times.
SolTrack was originally written in C/C++ by [Marc van der Sluys](http://marc.vandersluys.nl) of the department
of astrophysics of the Radboud University Nijmegen, the Netherlands and the Sustainable energy research group
of the HAN University of Applied Sciences in Arnhem, the Netherlands (now at the Netherlands Institute for
Nuclear and High-Energy Physics (Nikhef) in Amsterdam and the Institute for Gravitational and Subatomic
Physics (GRASP) at Utrecht University in the Netherlands), and Paul van Kan of the Sustainable energy research
group of the HAN University of Applied Sciences in Arnhem, the Netherlands.  The code has now been translated
to pure Python and can be used under the conditions of the EUPL 1.2 licence.

SolTrack can perform up to 1.5 million position calculations per second on a single 3.4 GHz core of my laptop
(the [C version of SolTrack](http://soltrack.sourceforge.net/) is about twice as fast) with an accuracy of
0.0030 ± 0.0016°.


## Installation

This package can be installed using `pip install soltrack`.  This should automatically install the dependency
packages `numpy`, `pandas`, `astroconst` and `astrotool` if they haven't been installed already (if you're not
using a Python version older than 3.7, you will need to install `dataclasses` as well).  If you are installing
by hand, ensure that these packages are installed as well.


## Example use

### Update from v0.1.x to v0.2.0

The update from SolTrack v0.1.4 to v0.3.0 is **not backwards compatible**.  Most public members/variables and
methods/functions, as well as dummy arguments have been renamed in order to comply better with the Python
standards.  I have tried to put as much nuisance into a single update as possible (as opposed to in multiple
updates, rather than as opposed to as little nuisance as possible), so that this will hopefully not happen
again in the future.  For public methods, the obsolescent old name is kept as an alias, and warnings are
issued when used, instructing how to adapt your code.  In many cases a variable or function name `likeThis`
will have been replaced to one `like_this`, i.e. upper case replaced with an underscore and lower case.  See
the [documentation](https://soltrack.readthedocs.io) for more details.


### Code example for a number or range of datetimes

```python """Compute the position of the Sun and its rise and set times for a vector of instances."""

# Create a Pandas DatetimeIndex range every 20 days 1 hour and 10 minutes, in my timezone:
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

### Code example for a single instant

Note that in most cases, you should use the [vector option](#code-example-for-a-number-or-range-of-datetimes)
instead, as it is faster starting from calculations for two instances.  See the section
[Performance](#performance) for more details.  The code listing below is provided for completeness only.


```python """Example Python script to compute the position of the Sun and its rise and set times for a single instant
and demonstrate some other features."""

from soltrack import SolTrack
import datetime as dt
import pytz as tz

# Set the geographic location to Arnhem, the Netherlands:
geo_lon =  5.950270  # Positive -> east of Greenwich (degrees)
geo_lat = 51.987380  # Positive -> northern hemisphere (degrees)

st = SolTrack(geo_lon, geo_lat, use_degrees=True)  # Same as above, using default values for all but use_degrees.

# Set SolTrack date and time using separate (UTC!) year, month, day, hour, minute and second variables:
st.set_date_and_time(2023, 7, 16,  6, 2, 49.217348)  # Date: 2023-07-16, time: 06:02:49.217348 UTC

# Alternatively, use a (localised) datetime object:
cet     = tz.timezone('Europe/Amsterdam')
my_date = dt.datetime(2023, 7, 16,  8, 2, 49, 217348)  # Same as above, in local time for TZ=+2 (08:02:49.217348 LT)
my_date = cet.localize(my_date)
st.set_date_time(my_date)  # Set SolTrack date and time using a Python datetime object.

# Compute the Sun's position:
st.compute_position()

# Compute the rise and set times of the Sun:
st.compute_rise_set()


# Write results to standard output:
print("Location:   %0.3lf E, %0.3lf N"  % (st.geo_longitude, st.geo_latitude))
# print("Date/time:  %4d-%02d-%02d %02d:%02d:%09.6lf" % (st.year, st.month, st.day,  st.hour, st.minute, st.second))
print("Date/time:  %s"                  % my_date)
print("JD:         %0.11lf"             % (st.julian_day))
print()

print("Ecliptic longitude, latitude:        %10.6lf° %10.6lf°"     % (st.longitude, 0.0))  # Note: latitude is always 0 in this model
print("Distance:                            %10.6lf°"              % (st.distance))
print("Right ascension, declination:        %10.6lf° %10.6lf°"     % (st._right_ascension_uncorr, st._declination_uncorr))
print("Uncorrected altitude:                            %10.6lf°"  % (st._altitude_uncorr))
print("Corrected azimuth, altitude:         %10.6lf° %10.6lf°"     % (st.azimuth, st.altitude))
print("Corrected hour angle, declination:   %10.6lf° %10.6lf°"     % (st.hour_angle, st.declination))
print()

print("Rise time:      %s,    azimuth:   %11.5lf" % (st.rise_time,     st.rise_azimuth))
print("Transit time:   %s,    altitude:  %11.5lf" % (st.transit_time,  st.transit_altitude))
print("Set time:       %s,    azimuth:   %11.5lf" % (st.set_time,      st.set_azimuth))
print()


# Change the location whilst keeping the same SolTrack object:
st.set_location(geo_lon, -geo_lat)

# Compute the current position of the Sun for the new location:
st.now()
st.compute_position()

print("Location:   %0.3lf E, %0.3lf N"                         % (st.geo_longitude, st.geo_latitude))
print("Date (UT):  %4d-%02d-%02d"                              % (st.year, st.month, st.day))
print("Time (UT):  %02d:%02d:%09.6lf"                          % (st.hour, st.minute, st.second))
print("Corrected azimuth, altitude:         %10.6lf° %10.6lf°" % (st.azimuth, st.altitude))
```


### Code example for SolTrack v0.1.x (old)
```python
"""Example Python script to compute the position of the Sun and its rise and set times for a single instant."""

import soltrack as st

# Set preferences (all are False by default):
useDegrees             = True   # Input (geographic position) and output are in degrees
useNorthEqualsZero     = True   # Azimuth: 0 = South, pi/2 (90deg) = West  ->  0 = North, pi/2 (90deg) = East
computeRefrEquatorial  = True   # Compure refraction-corrected equatorial coordinates (Hour angle, declination)
computeDistance        = True   # Compute the distance to the Sun

# Set up geographical location (in degrees, since useDegrees=True) in a SolTrack Location dataclass object:
loc = st.Location(5.950270, 51.987380)  # longitude (>0: east of Greenwich),  latitude (>0: northern hemisphere)

# Set (UT!) date and time in a SolTrack Time dataclass object:
time = st.Time(2045, 7, 16,  6, 2, 49.217348)  # Date: 2045-07-16, time: 06:02:49.217348


# Compute positions - returns a st.Position object:
pos = st.computeSunPosition(loc, time, useDegrees, useNorthEqualsZero, computeRefrEquatorial, computeDistance)

# Compute rise and set times - returns a st.RiseSet object:
riseSet = st.computeSunRiseSet(loc, time, 0.0, useDegrees, useNorthEqualsZero)


# Write results to standard output:
print("Location:  %0.3lf E, %0.3lf N" % (loc.longitude, loc.latitude))
print("Date:      %4d %2d %2d" % (time.year, time.month, time.day))
print("Time:      %2d %2d %9.6lf" % (time.hour, time.minute, time.second))
print("JD:        %0.11lf" % (pos.julianDay))
print()

print("Ecliptic longitude, latitude:        %10.6lf° %10.6lf°" % (pos.longitude, 0.0))
print("Right ascension, declination:        %10.6lf° %10.6lf°" % (pos.rightAscension, pos.declination))
print("Uncorrected altitude:                            %10.6lf°" % (pos.altitude))
print("Corrected azimuth, altitude:         %10.6lf° %10.6lf°" % (pos.azimuthRefract, pos.altitudeRefract))
print("Corrected hour angle, declination:   %10.6lf° %10.6lf°" % (pos.hourAngleRefract, pos.declinationRefract))
print()

print("Rise time:      %11.5lf,    azimuth:   %11.5lf" % (riseSet.riseTime, riseSet.riseAzimuth))
print("Transit time:   %11.5lf,    altitude:  %11.5lf" % (riseSet.transitTime, riseSet.transitAltitude))
print("Set time:       %11.5lf,    azimuth:   %11.5lf" % (riseSet.setTime, riseSet.setAzimuth))
print()
```

## Performance

### Performance of position calculations for an array/vector or range of instants

From v0.2.0 on, SolTrack uses Pandas to allow the use of arrays or vectors of datetimes (Series,
DatetimeIndex, ndarrays of datetime64, ...).  This causes some overhead, which is relatively more significant
for small numbers.  SolTrack performs fastest (per position calculation and on my laptop) for about
10<sup>5</sup>-10<sup>6</sup> calculations:

| N<sub>calc</sub> (-) | Time (s) | Speed (/s) | Speed (millions/s) |
|----------------------|----------|------------|--------------------|
| 1e0                  | 0.00081  | 1.23e+03   | 0.001              |
| 1e1                  | 0.00080  | 1.25e+04   | 0.013              |
| 1e2                  | 0.00105  | 9.52e+04   | 0.095              |
| 1e3                  | 0.00256  | 3.91e+05   | 0.39               |
| 1e4                  | 0.0132   | 7.58e+05   | 0.76               |
| 1e5                  | 0.0688   | 1.45e+06   | 1.45               |
| 3e5                  | 0.1960   | 1.53e+06   | 1.53               |
| 1e6                  | 0.6660   | 1.50e+06   | 1.50               |
| 3e6                  | 2.068    | 1.45e+06   | 1.45               |
| 1e7                  | 8.691    | 1.15e+06   | 1.15               |

These benchmarks were done on a single CPU core, capped at 3.4GHz.  The cpu and core were always the same, and
the minimum of 10 benchmarks is listed.  Timezone-naive datetimes (representing UTC) were used; using
timezone-aware datetimes slows down the code by a factor of ~4.6 compared to the numbers in the table.


### Performance of calculations of rise, transit and set times

Benchmarks for rise, transit and set computations were performed in the same manner as above, but using
`return_datetimes=False` in the `st.compute_rise_set()` call, which returns rise, transit and set times in
decimal hours.  If instead datetimes are desired, the code becomes a factor of about 3.1 times slower.  The
use of `utc=True` has little effect, and using timezone-aware datetimes instead of timezone-naive ones adds
6.1% to the computational times listed below.

| N<sub>calc</sub> (-) | Time (s) | Speed (/s) |
|----------------------|----------|------------|
| 1                    | 0.0020   | 500.0      |
| 10                   | 0.0134   | 746.3      |
| 100                  | 0.1247   | 801.9      |
| 1000                 | 1.207    | 828.5      |
| 1e4                  | 12.00    | 833.3      |


### Performance of position calculation for a single instant

Because of the overhead of Pandas and datetime-like objects, SolTrack has actually slowed down by a factor of
~10.4 between versions 0.1.4 and 0.2.0 when doing single calculations.  I don't consider that a big issue,
since a single position calculation still takes less than a millisecond, so that the effect is not noticable
for humans.  It does become noticable for large numbers, but then arrays can be used, which make the code
~1000x faster.  In fact, the array version is faster starting from two iterations, so there is usually very
little reason to use the single-datetime option.

Single calculations can be made about 19% faster by providing datetimes in UTC and specifying `utc=True` in
the `st.set_date_time()` method.  The use of timezone-naive versus timezone-aware datetimes has little
influence on the performance for single calculations, and calculation times scale linearly when doing multiple
calls (since the must then be done in a Python loop).

For single datetimes, 1000 calculations take about 0.68 seconds on the single CPU core of my laptop capped at
3.4 GHz, and about 0.55 seconds with `utc=True`.


## SolTrack pages

* [Van der Sluys & Van Kan (2022)](https://arxiv.org/abs/2209.01557): SolTrack: a free, fast and accurate routine to
  compute the position of the Sun.  Scientific paper with all technical details.

* [Pypi](https://pypi.org/project/soltrack/): SolTrack Python package
* [GitHub](https://github.com/MarcvdSluys/SolTrack-Python): SolTrack Python source code
* [ReadTheDocs](https://soltrack.readthedocs.io/en/latest/): SolTrack Python documentation

* [SolTrack for C/C++](http://soltrack.sourceforge.net/) on SourceForge
* [SolTrack for Arduino](https://github.com/MarcvdSluys/SolTrack-Arduino) on GitHub


## Author and licence

* Author: Marc van der Sluys
* Contact: http://marc.vandersluys.nl
* Licence: [EUPL 1.2](https://www.eupl.eu/1.2/en/)


## See also

* [AstroConst](https://pypi.org/project/astroconst/): a Python package that provides astronomical constants.
* [AstroTool](https://pypi.org/project/astrotool/): a Python package for astronomical calculations in Python
  or on the command line.
* [elp-mpp02](https://pypi.org/project/elp-mpp02/): accurate Moon positions using the lunar solution ELP/MPP02
  in Python.
* [libTheSky](http://libthesky.sourceforge.net/): a Fortran library to compute the positions of celestial
  bodies (Sun, Moon, planets, stars, asteroids, comets) and events (e.g. lunar phases) with great accuracy.
* [SolarEnergy](https://pypi.org/project/solarenergy/): A Python module to do simple modelling in the field of
  solar energy.


## References

* [Van der Sluys & Van Kan (2022)](https://arxiv.org/abs/2209.01557): SolTrack: a free, fast and accurate routine to
  compute the position of the Sun.  Scientific paper with all technical details.

* [Celestial mechanics in a nutshell (CMiaNS)](https://cmians.sourceforge.io/): online living document.
* Meeus, [Astronomical algorithms](https://www.willbell.com/math/MC1.HTM), 2nd Ed.
* The C/C++ and Python codes have been adapted from the Fortran implementation in
  [libTheSky](http://libthesky.sourceforge.net/), which contains many references.


<sub>Copyright (c) 2019-2022 Marc van der Sluys</sub>
