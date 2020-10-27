
"""SolTrack: a simple, free, fast and accurate C routine to compute the position of the Sun.
    
  Copyright (c) 2014-2020 Marc van der Sluys, Paul van Kan and Jurgen Reintjes, 
  Sustainable Energy research group, HAN University of applied sciences, Arnhem, The Netherlands
   
  This file is part of the SolTrack package, see: http://soltrack.sourceforge.net SolTrack is derived from
  libTheSky (http://libthesky.sourceforge.net) under the terms of the GPL v.3
  
  This is free software: you can redistribute it and/or modify it under the terms of the GNU Lesser General
  Public License as published by the Free Software Foundation, either version 3 of the License, or (at your
  option) any later version.
  
  This software is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the
  implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public
  License for more details.
  
  You should have received a copy of the GNU Lesser General Public License along with this code.  If not, see
  <http://www.gnu.org/licenses/>.

"""

import pytz as tz
import datetime as dt


class Location(object):
    """Class containing the geographic location to compute the Sun position for."""
    
    sinLat:      float = 0.0
    cosLat:      float = 0.0
    
    def __init__(self, longitude, latitude, pressure=101.0, temperature=283.0):
        self.longitude = longitude
        self.latitude  = latitude
        
        self.pressure    = pressure
        self.temperature = temperature
    

class Time(object):
    """Class containing the date and time (in UT) to compute the Sun position for."""
    
    def __init__(self, year=2000,month=1,day=1, hour=12,minute=0,second=0.0):
        self.year = year
        self.month = month
        self.day = day
        self.hour = hour
        self.minute = minute
        self.second = second

    
    def datetime2st(dtObj):
        """Convert a datetime object to a SolTrack time object.
        
           Args:
               dtObj (datetime):  Date and time in a Python datetime object.
               
           Returns:
               Time:  Date and time in a SolTrack time object.
        
        """
        
        time = Time()  # Create a SolTrack Time object
        
        utc         = dtObj.astimezone(tz.utc)  # Convert from local time to UTC
        
        time.year   = utc.year
        time.month  = utc.month
        time.day    = utc.day
        
        time.hour   = utc.hour
        time.minute = utc.minute
        time.second = utc.second + utc.microsecond/1e6
        
        return time

    
    def now():
        """Return the current system time as a SolTrack time object.
        
           Returns:
               Current system date and time in a SolTrack time object.
        """
        
        return Time.datetime2st(dt.datetime.now())
        
        
class Position(object):
    """Class containing the position of the Sun and related variables."""
    
    julianDay:           float = 0.0
    tJD:                 float = 0.0
    tJC:                 float = 0.0
    tJC2:                float = 0.0
    
    longitude:           float = 0.0
    distance:            float = 0.0
    
    obliquity:           float = 0.0
    cosObliquity:        float = 0.0
    nutationLon:         float = 0.0
    
    rightAscension:      float = 0.0
    declination:         float = 0.0
    agst:                float = 0.0
    
    altitude:            float = 0.0
    altitudeRefract:     float = 0.0
    azimuthRefract:      float = 0.0
    
    hourAngleRefract:    float = 0.0
    declinationRefract:  float = 0.0


class RiseSet(object):
    """Class containing rise,transit and set times of the Sun and their azimuths/altitudes."""
    
    riseTime:         float = 0.0
    transitTime:      float = 0.0
    setTime:          float = 0.0
    
    riseAzimuth:      float = 0.0
    transitAltitude:  float = 0.0
    setAzimuth:       float = 0.0
    
    
    
    


def copyObject(oldInst):
    """Deep copy an existing object by creating a new instance and copying its members.
    
    Note:
      - Simply copying an object copies it's *address*, and doesn't make a *data copy*.
    
    Parameters:
      oldInst (Class):  the existing object.
    
    Returns:
      (Class):  the new object.
    
    See:
      - https://stackoverflow.com/a/4794254/1386750
    
    """
    
    import copy
    
    newInst = copy.deepcopy(oldInst)
    
    return newInst
    
