# -*- coding: utf-8 -*-
# SPDX-License-Identifier: EUPL-1.2
#  
#  Copyright (c) 2019-2024  Marc van der Sluys - marc.vandersluys.nl
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
    """Class containing the date and time (in UTC) to compute the Sun position for."""
    
    
    def setDateAndTime(self, year=2000,month=1,day=1, hour=12,minute=0,second=0.0):
        """This function is obsolescent and will be removed in a future version.  Use set_date_and_time()
        instead."""
        _warn_obsolescent('setDateAndTime', 'set_date_and_time', rename=True)
        return self.set_date_and_time(year,month,day, hour,minute,second)
    
    
    def set_date_and_time(self, year=2000,month=1,day=1, hour=12,minute=0,second=0.0):
        """Set the SolTrack date and time using UTC year, month, day, hour, minute and second.
        
           Parameters:
             year   (int):    year of date.
             month  (int):    month of date.
             day    (int):    day of date.
        
             hour   (int):    hour of day  (default=0).
             minute (int):    minute of time  (default=0).
             second (float):  second of time  (default=0).
        
           Note:
             Use set_date_time() instead if you have a Python datetime object.
        """
        
        # Combine the date/time values into a single "2D" array with a single row and the original values as
        # columns, and convert it to a Pandas df:
        df = pd.DataFrame(np.vstack([year,month,day, hour,minute,second]).transpose(),
                          columns=['year','month','day', 'hour','minute','second'])
        
        df = pd.DataFrame(data=pd.to_datetime(df), columns=['UTC'])  # Convert the date+time columns into a single datetime column
        
        self.utc = df.loc[:, 'UTC']
        utc_dt = self.utc.dt
        
        self.year   = utc_dt.year.to_numpy()
        self.month  = utc_dt.month.to_numpy()
        self.day    = utc_dt.day.to_numpy()
        
        self.hour   = utc_dt.hour.to_numpy()
        self.minute = utc_dt.minute.to_numpy()
        self.second = utc_dt.second.to_numpy() + utc_dt.microsecond.to_numpy()/1e6
        
        return
    
    
    def setDateTime(self, dt_obj, utc=False):
        """This function is obsolescent and will be removed in a future version.  Use setDateTime()
        instead."""
        _warn_obsolescent('setDateTime', 'set_date_time', rename=True)
        return self.set_date_time(dt_obj, utc)
    
    
    def set_date_time(self, dt_obj, utc=False):
        """Set the SolTrack date and time using a (local) Python datetime object.
        
           Parameters:
               dt_obj (datetime(64)):  Date and time in a Python datetime object (UTC if timezone naive).
        
           Returns:
               Time:  Date and time in a SolTrack time object.
        
           Note:
             Use set_date_and_time() instead if you have year, month,day, hour, minute and second as separate
             variables.
        """
        
        dt_obj = np.asarray(np.copy(dt_obj))  # Copy and typecast to numpy.ndarray
        
        scalar_input = False
        if dt_obj.ndim == 0:
            dt_obj = dt_obj[None]  # Makes dt_obj a 1D array.  Comment: use np.newaxis instead?
            scalar_input = True
        
        if utc:  # up to ~29% faster if datetimes are known to be UTC.  Create DatetimeIndex with UTC times directly:
            if np.ndim(dt_obj) == 0:  # Scalar, needs to be converted using [array]:
                self.utc = pd.to_datetime(np.asarray([dt_obj]))  # DatetimeIndex
            else:
                self.utc = pd.to_datetime(np.asarray(dt_obj))    # DatetimeIndex
            
        else:  # Create DatetimeIndex with local times, then convert to UTC:
            if np.ndim(dt_obj) == 0:  # Scalar, needs to be converted using [array]:
                self.lt = pd.to_datetime(np.asarray([dt_obj]))   # DatetimeIndex
            else:
                self.lt = pd.to_datetime(np.asarray(dt_obj))     # DatetimeIndex
                
            # Ensure timestamps are in UTC:
            self.utc = pd.to_datetime(self.lt, utc=True)  # utc=True: make timezone aware (if not already), and set TZ=UTC (converting if needed).
            self.utc = self.utc.tz_convert(None)          # 8±2% faster if left out.  Convert to UTC tz naive (i.e. convert to UTC and remove tz info)
        
        
        # Note: the six lines below increase the cpu time of a compute_position() by a factor of 2.7!  The
        # .to_numpy() method parts have a negligable contribution to that.
        # CHECK1: using .utc.to_julian_date() would not require self.year-second here, but is slightly slower
        # and gives wrong results in computeRiseSet().  See CHECK1 in those places.
        self.year   = self.utc.year.to_numpy()
        self.month  = self.utc.month.to_numpy()
        self.day    = self.utc.day.to_numpy()
        
        self.hour   = self.utc.hour.to_numpy()
        self.minute = self.utc.minute.to_numpy()
        self.second = self.utc.second.to_numpy() + self.utc.microsecond.to_numpy()/1e6
        
        if scalar_input:
            self.year    = int(np.squeeze(self.year))      # Array -> scalar, int
            self.month   = int(np.squeeze(self.month))     # Array -> scalar, int
            self.day     = int(np.squeeze(self.day))       # Array -> scalar, int
            self.hour    = int(np.squeeze(self.hour))      # Array -> scalar, int
            self.minute  = int(np.squeeze(self.minute))    # Array -> scalar, int
            self.second  = float(np.squeeze(self.second))  # Array -> scalar, float
        
        return
    
    
    def now(self):
        """Return the current system time as a SolTrack time object.
        
           Returns:
               Current system date and time in a SolTrack time object.
        """
        
        self.set_date_time(dt.datetime.now())
        
        return
    
    
def _warn_obsolescent(old_name, new_name, rename=False, extra=False):
    """Warn that a function is obsolescent and will be removed.  Indicate whether this concerns a simple rename, possibly with extra features."""
    import sys
    sys.stderr.write('\nWarning: the SolTrack function '+old_name+'() is obsolescent and will be removed in a future version.')
    sys.stderr.write('  Use '+new_name+'() instead.')
    if rename:
        if extra:
            sys.stderr.write('  The interface has not changed much; a simple search and replace for the function names should suffice, but some dummy variables and class members have also be renamed, so please see the documentation  in case a simple rename does not work and for new features.\n\n')
        else:
            sys.stderr.write('  The interface has not changed; a simple search and replace for the function names may suffice, but some dummy variables and class members have also be renamed, so please see the documentation in case a simple rename does not work.\n\n')
    else:
        sys.stderr.write('  Please see the documentation for details.\n\n')
    return
