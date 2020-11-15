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
import numpy as np
import datetime as dt
import pytz as tz


@dataclass
class Time:
    """Class containing the date and time (in UT) to compute the Sun position for."""
    
    
    def setDateAndTime(self, year=2000,month=1,day=1, hour=12,minute=0,second=0.0):
        """Set the SolTrack date and time using UTC year, month, day, hour, minute and second.
        
           Parameters:
             year   (int):    year of date.
             month  (int):    month of date.
             day    (int):    day of date.
             
             hour   (int):    hour of day  (default=0).
             minute (int):    minute of time  (default=0).
             second (float):  second of time  (default=0).
        
           Note:
             Use setDateTime() instead if you have a Python datetime object.
        
        """
        
        # If a Time object does not yet exist, create it:
        # if(not hasattr(self, "time")):
        #     self.time = Time()  # Create a SolTrack Time object
        
        self.year   = year
        self.month  = month
        self.day    = day
        
        self.hour   = hour
        self.minute = minute
        self.second = second
        
        return
    
    
    def setDateTime(self, dtObj):
        """Set the SolTrack date and time using a (local) Python datetime object.
        
           Parameters:
               dtObj (datetime):  Date and time in a Python datetime object.
               
           Returns:
               Time:  Date and time in a SolTrack time object.
        
           Note:
             Use setDateAndTime() instead if you have year, month,day, hour, minute and second as separate
             variables.

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
        self.setDateTime(dt.datetime.now())
        
        return
        
        
    def _computeJulianDay(self, year, month, day,  hour, minute, second):
        """Compute the Julian Day from the date and time.
        
        Note:
          - Gregorian calendar only (>~1582).
        
        Parameters:
          year    (int):    Year of date.
          month   (int):   Month of date.
          day     (int):     Day of date.
          hour    (int):    Hour of time.
          minute  (int):  Minute of time.
          second  (int):  Second of time.
        
        Returns:
          float:  Julian day for the given date and time.
        
        """
        
        if(month <= 2):  # Treat Jan, Feb as months 13, 14 of the previous year
            year  -= 1
            month += 12
            
        tmp1 = np.floor(year/100.0)
        tmp2 = 2 - tmp1 + np.floor(tmp1/4.0)
        
        dDay = day + hour/24.0 + minute/1440.0 + second/86400.0
        
        self.julianDay = np.floor(365.250*(year+4716)) + np.floor(30.60010*(month+1)) + dDay + tmp2 - 1524.5
        
        return
    
    
    
