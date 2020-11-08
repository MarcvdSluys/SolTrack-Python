#  Copyright (c) 2019-2020  Marc van der Sluys - marc.vandersluys.nl
#   
#  This file is part of the SolTrack Python package,
#  see: http://soltrack.sf.net
#   
#  This is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License
#  as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
#  
#  This software is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied
#  warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
#  
#  You should have received a copy of the GNU General Public License along with this code.  If not, see
#  <http://www.gnu.org/licenses/>.


from dataclasses import dataclass
import datetime as dt
import pytz as tz


@dataclass
class Time:
    """Class containing the date and time (in UT) to compute the Sun position for."""
    
    
    def setDateTime(self, year=2000,month=1,day=1, hour=12,minute=0,second=0.0):
        self.year   = year
        self.month  = month
        self.day    = day
        
        self.hour   = hour
        self.minute = minute
        self.second = second
        
        return
    
    
    def datetime2st(self, dtObj):
        """Convert a datetime object to a SolTrack time object.
        
           Args:
               dtObj (datetime):  Date and time in a Python datetime object.
               
           Returns:
               Time:  Date and time in a SolTrack time object.
        
        """
        
        utc         = dtObj.astimezone(tz.utc)  # Convert from local time to UTC
        
        self.year   = utc.year
        self.month  = utc.month
        self.day    = utc.day
        
        self.hour   = utc.hour
        self.minute = utc.minute
        self.second = utc.second + utc.microsecond/1e6
        
        return
    
    
    def now(self):
        """Return the current system time as a SolTrack time object.
        
           Returns:
               Current system date and time in a SolTrack time object.
        """
        self.datetime2st(dt.datetime.now())
        
        return
        
        
