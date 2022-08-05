# -*- coding: utf-8 -*-
# SPDX-License-Identifier: EUPL-1.2
#  
#  Copyright (c) 2019-2022  Marc van der Sluys - marc.vandersluys.nl
#  
#  This file is part of the SolTrack Python package,
#  see: http://soltrack.sf.net
#  
#  This is free software: you can redistribute it and/or modify it under the terms of the
#  European Union Public Licence 1.2 (EUPL 1.2).
#  
#  This software is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
#  without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
#  See the EU Public Licence for more details.
#  
#  You should have received a copy of the European Union Public Licence along with this code.
#  If not, see <https://www.eupl.eu/1.2/en/>.


from dataclasses import dataclass
import numpy as np
import pandas as pd
import datetime as dt


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
        
        # Combine the date/time values into a single "2D" array with a single row and the original values as
        # columns, and convert it to a Pandas df:
        self.df = pd.DataFrame(np.vstack([year,month,day, hour,minute,second]).transpose(),
                               columns=['year','month','day', 'hour','minute','second'])
        
        self.df = pd.DataFrame(data=pd.to_datetime(self.df), columns=['UTC'])  # Convert the date+time columns into a single datetime column
        
        self.year   = self.df.loc[:, 'UTC'].dt.year
        self.month  = self.df.loc[:, 'UTC'].dt.month
        self.day    = self.df.loc[:, 'UTC'].dt.day
        
        self.hour   = self.df.loc[:, 'UTC'].dt.hour
        self.minute = self.df.loc[:, 'UTC'].dt.minute
        self.second = self.df.loc[:, 'UTC'].dt.second + self.df.loc[:, 'UTC'].dt.microsecond/1e6
        
        return
    
    
    def setDateTime(self, dtObj):
        """Set the SolTrack date and time using a (local) Python datetime object.
        
           Parameters:
               dtObj (datetime(64)):  Date and time in a Python datetime object (UTC if timezone naive).
               
           Returns:
               Time:  Date and time in a SolTrack time object.
        
           Note:
             Use setDateAndTime() instead if you have year, month,day, hour, minute and second as separate
             variables.
        """
        
        if np.ndim(dtObj) == 0:  # Scalar, needs to be converted using [array]:
            self.df = pd.DataFrame(np.asarray([dtObj]), columns=['LT'])  # DF containing Series containing pd._libs.tslibs.timestamps.Timestamp == datetime64[ns]
        else:
            self.df = pd.DataFrame(np.asarray(dtObj), columns=['LT'])  # DF containing Series containing pd._libs.tslibs.timestamps.Timestamp == datetime64[ns]
        
        # Ensure timestamps are in UTC:
        self.df['UTC'] = pd.to_datetime(self.df['LT'], utc=True)  # utc=True: make timezone aware (if not already), and set TZ=UTC (converting if needed).
        self.df['UTC'] = self.df['UTC'].dt.tz_convert(None)       # Convert to UTC tz naive (i.e. convert to UTC and remove tz info)
        
        self.year   = self.df.loc[:, 'UTC'].dt.year
        self.month  = self.df.loc[:, 'UTC'].dt.month
        self.day    = self.df.loc[:, 'UTC'].dt.day
        
        self.hour   = self.df.loc[:, 'UTC'].dt.hour
        self.minute = self.df.loc[:, 'UTC'].dt.minute
        self.second = self.df.loc[:, 'UTC'].dt.second + self.df.loc[:, 'UTC'].dt.microsecond/1e6
        
        return
        
    
    def now(self):
        """Return the current system time as a SolTrack time object.
        
           Returns:
               Current system date and time in a SolTrack time object.
        """
        self.setDateTime(dt.datetime.now())
        
        return
        
        
    def _computeJulianDayScalar(self, year, month, day,  hour, minute, second):
        """Compute the Julian Day from the date and time.
        
        Parameters:
          year    (int):    Year of date.
          month   (int):   Month of date.
          day     (int):     Day of date.
          hour    (int):    Hour of time.
          minute  (int):  Minute of time.
          second  (int):  Second of time.
        
        Returns:
          float:  Julian day for the given date and time.
        
        Note:
          - Gregorian calendar only (>~1582).
        
        """
        
        if(month <= 2):  # Treat Jan, Feb as months 13, 14 of the previous year
            year  -= 1
            month += 12
            
        tmp1 = np.floor(year/100.0)
        tmp2 = 2 - tmp1 + np.floor(tmp1/4.0)
        
        dDay = day + hour/24.0 + minute/1440.0 + second/86400.0
        
        self.julianDay = np.floor(365.250*(year+4716)) + np.floor(30.60010*(month+1)) + dDay + tmp2 - 1524.5
        
        return
    
    
    
    def _computeJulianDay(self, year,month,day, hour,minute,second, jd_start_greg=2299160.5):
        """Compute the Julian Day from given year, month, day, hour, minute and second.
        
        Args:
          year (int):             Year CE.  Note that year=0 = 1 BCE, year=-1 = 2 BCE, etc.
          month (int):            Month number of year (1-12).
          day (float):            Day of month with fraction (1.0-31.999...).
        
          hour (int):             Hour of day (0-23).
          minute (int):           Minute of hour (0-59).
          second (float):         Second of minute (0.0-59.999...).
        
          jd_start_greg (float):  JD of start of Gregorian calendar (optional; default=2299160.5 = 1582-10-15.0).
        
        Returns:
          float:  jd: Julian day (days).
          
        Note:
          - The JD will be in the same timezone as the date and time (UT for the offical JD).
          - Decimals can be used in the day to take into account the time of day other than midnight, e.g. 1.5 for
            noon on the first day of the month.
        """
        
        # Copy and typecast input to numpy.ndarrays:
        year  = np.asarray(np.copy(year))
        month = np.asarray(np.copy(month))
        day   = np.asarray(np.copy(day) + hour/24 + minute/1440 + second/86400)
        
        # Jan/Feb are month 13/14 of the previous year:
        year[month  <= 2] -= 1
        month[month <= 2] += 12
        
        # JD for Julian calendar (ensure it always is an array):
        jd = np.asarray(np.floor(365.25*(year+4716)) + np.floor(30.6001*(month+1)) + day - 1524.5)
        
        # Apply correction for Gregorian calendar:
        sel      = jd >= jd_start_greg                # Select cases for Greg.cal.
        cent_1   = np.floor(year[sel]/100.0)          # Count: (century - 1)
        jd[sel] += 2 - cent_1 + np.floor(cent_1/4.0)  # Offset Julian-Gregorian
        
        self.julianDay = jd
        
        return
    
    
