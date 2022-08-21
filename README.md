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

### Example use for SolTrack v0.1.x ###
```python
"""Example Python script to compute the position of the Sun and its rise and set times for a single instant
and demonstrate some other features."""

from soltrack import SolTrack
import datetime as dt

# Set the geographic location to Arnhem, the Netherlands:
geoLon =  5.950270  # Positive -> east of Greenwich (degrees)
geoLat = 51.987380  # Positive -> northern hemisphere (degrees)

# Create a SolTrack instance and specify preferences (the first two are False by default):
# - useDegrees: use degrees instead of radians for input (geographic location) and output (position).
#               Default: False.
# - useNorthEqualsZero: Set an azimuth of 0 to mean north instead of south.  Default: False.
# - computeRefrEquatorial: compute refraction-corrected equatorial coordinates.  Default: True.
# - computeDistance: compute the physical Eart-Sun distance in AU.  Default: True.
# 
# st = SolTrack(geoLon, geoLat, useDegrees=True, useNorthEqualsZero=False, computeRefrEquatorial=True,
#               computeDistance=True)

st = SolTrack(geoLon, geoLat, useDegrees=True)  # Same as above, using default values for all but useDegrees.

# Set SolTrack date and time using separate (UT!) year, month, day, hour, minute and second variables:
st.setDateAndTime(2045, 7, 16,  6, 2, 49.217348)  # Date: 2045-07-16, time: 06:02:49.217348 UTC

# Alternatively, use a datetime object:
myDate = dt.datetime(2045, 7, 16,  8, 2, 49, 217348)  # Same as above, in local time for TZ=+2 (08:02:49.217348 LT)
st.setDateTime(myDate)  # Set SolTrack date and time using a Python datetime object.


# Compute the Sun's position:
st.computePosition()

# Compute the rise and set times of the Sun:
st.computeRiseSet()


# Write results to standard output:
print("Location:   %0.3lf E, %0.3lf N"  % (st.geoLongitude, st.geoLatitude))
print("Date (UT):  %4d-%02d-%02d"       % (st.year, st.month, st.day))
print("Time (UT):  %02d:%02d:%09.6lf"   % (st.hour, st.minute, st.second))
print("JD:         %0.11lf"             % (st.julianDay))
print()

print("Ecliptic longitude, latitude:        %10.6lf° %10.6lf°"     % (st.longitude, 0.0))
print("Distance:                            %10.6lf°"              % (st.distance))
print("Right ascension, declination:        %10.6lf° %10.6lf°"     % (st._rightAscensionUncorr, st._declinationUncorr))
print("Uncorrected altitude:                            %10.6lf°"  % (st._altitudeUncorr))
print("Corrected azimuth, altitude:         %10.6lf° %10.6lf°"     % (st.azimuth, st.altitude))
print("Corrected hour angle, declination:   %10.6lf° %10.6lf°"     % (st.hourAngle, st.declination))
print()

print("Rise time:      %11.5lf,    azimuth:   %11.5lf" % (st.riseTime, st.riseAzimuth))
print("Transit time:   %11.5lf,    altitude:  %11.5lf" % (st.transitTime, st.transitAltitude))
print("Set time:       %11.5lf,    azimuth:   %11.5lf" % (st.setTime, st.setAzimuth))
print()


# Change the location whilst keeping the same SolTrack object:
st.setLocation(geoLon, -geoLat)

# Compute the current position of the Sun for the new location:
st.now()
st.computePosition()

print("Location:   %0.3lf E, %0.3lf N"                         % (st.geoLongitude, st.geoLatitude))
print("Date (UT):  %4d-%02d-%02d"                              % (st.year, st.month, st.day))
print("Time (UT):  %02d:%02d:%09.6lf"                          % (st.hour, st.minute, st.second))
print("Corrected azimuth, altitude:         %10.6lf° %10.6lf°" % (st.azimuth, st.altitude))
print()
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
