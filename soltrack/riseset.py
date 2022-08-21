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

from .data import Constants, Parameters
import pandas as pd
import astrotool as at


@dataclass
class RiseSet(Constants, Parameters):
    """Class concerning the rise, transit and set times and positions of the Sun and related attributes and methods."""
    
    
    def __init__(self):
        Parameters.__init__(self)
        
        # Rise, transit and set time:
        self.riseTime:               float = 0.0;      """Rise time of the Sun (hours LT or UTC)"""
        self.transitTime:            float = 0.0;      """Transit time of the Sun (hours LT or UTC)"""
        self.setTime:                float = 0.0;      """Set time of the Sun (hours LT or UTC)"""
        
        # Rise, transit and set position:
        self.riseAzimuth:            float = 0.0;      """Rise azimuth of the Sun (radians)"""
        self.transitAltitude:        float = 0.0;      """Transit altitude of the Sun (radians)"""
        self.setAzimuth:             float = 0.0;      """Set azimuth of the Sun (radians)"""
        
    
    
    def computeRiseSet(self, rsAlt=0.0, accur=1e-5, return_datetimes=True):
        """Compute rise, transit and set times for the Sun, as well as their azimuths/altitude.
        
        Parameters:
          rsAlt            (float):  Altitude to return rise/set data for (radians; optional, default=0.0 meaning actual rise/set).  Set rsAlt>pi/2 to compute transit only.
          accur:           (float):  Accuracy (rad).  Default: 1e-5 rad ~ 0.14s.  Don't make this smaller than 1e-16.
          return_datetimes (bool):   Return times as datetimes rather than decimal hours.  Defaults to True.
        
        Note:
          - rise/set/transit times are in the LOCAL timezone used for the input (hence UTC if UTC was used).
          - if rsAlt == 0.0, actual rise and set times are computed
          - if rsAlt != 0.0, the routine calculates when alt = rsAlt is reached
          - returns times, rise/set azimuth and transit altitude in the class riseSet
         
          See:
            - subroutine riset() in riset.f90 from libTheSky (libthesky.sf.net) for more info
        """
        
        trTime = np.empty(0);  riTime = np.empty(0);  seTime = np.empty(0)
        trAlt  = np.empty(0);  riAz   = np.empty(0);  seAz   = np.empty(0)
        
        # CHECK: we should use local time for local Sun rise/transit/set, but what do we do if UTC is available only?
        if self.lt is not None:
            origDT = self.lt
        else:
            origDT = self.utc  # Try to correct for tz using geoLon?  Important is that 'midnight' lies between sunset and sunrise (CHECK: true? sufficient?)
        
            
        for OrigDTi in origDT:  # Loop over dates in array
            tmHrs, azalt = self._computeRiseSetSingle(OrigDTi, rsAlt, accur, return_datetimes)  # Compute r/s/t for a single date
            
            # Returned datetimes are timezone naive.  Add the timezone of the original date, if any:
            if return_datetimes:
                tmHrs = tmHrs.dt.tz_localize(OrigDTi.tz)
            
            # Store intermediate results:
            trTime = np.append(trTime, tmHrs[0])  # Transit time - hours
            riTime = np.append(riTime, tmHrs[1])  # Rise time - hours
            seTime = np.append(seTime, tmHrs[2])  # Set time - hours
            
            trAlt = np.append(trAlt, azalt[0])  # Transit altitude (rad)
            riAz  = np.append(riAz,  azalt[1])  # Rise azimuth (rad)
            seAz  = np.append(seAz,  azalt[2])  # Set azimuth (rad)
        
        
        # Store final results:
        self.transitTime      = trTime
        self.riseTime         = riTime
        self.setTime          = seTime
        self.transitAltitude  = trAlt
        self.riseAzimuth      = riAz
        self.setAzimuth       = seAz
        
        return
    
    
    def _computeRiseSetSingle(self, origDT, rsAlt, accur, return_datetimes):
        """Compute rise, transit and set times for the Sun, as well as their azimuths/altitude for a single date.
        
        Parameters:
          origDt          (datetime):  Original datetime of the date to compute the rise, set and transit for..
          rsAlt              (float):  Altitude to return rise/set data for (radians).  Set rsAlt>pi/2 to compute transit only.
          accur:             (float):  Accuracy (rad).  Don't make this smaller than 1e-16.
          return_datetimes    (bool):  Return times as datetimes rather than decimal hours.
        
        See computeRiseSet() for more details.
        """
        
        rsa = -0.8333/self._R2D               # Standard altitude for the Sun in radians
        if(abs(rsAlt) > 1.e-9): rsa = rsAlt   # Use a user-specified altitude
        
        tmRad = np.zeros(3)
        azalt = np.zeros(3)
        alt=0.0;  ha=0.0; h0=0.0
        
        
        # We need a local SolTrack instance for the same location (but in radians!), but with different settings
        # (radians, south=0, need equatorial coordinates but not the distance), and independent times and
        # positions:
        if(self.param._useDegrees):
            st = self.__class__(self.geoLongitude/self._R2D, self.geoLatitude/self._R2D, useDegrees=False,
                                useNorthEqualsZero=False, computeRefrEquatorial=True, computeDistance=False)
        else:
            st = self.__class__(self.geoLongitude, self.geoLatitude, useDegrees=False, useNorthEqualsZero=False,
                                computeRefrEquatorial=True, computeDistance=False)
        
        # Set date and time to midnight of the desired date:
        midnight = origDT.normalize()  # Midnight for the date in OrigDT, keeping the timezone
        st.setDateTime(midnight)
        
        # Compute the Sun's position:
        st.computePosition()
        
        
        # Compute transit, rise and set times:
        agst0 = st._agst      # AGST for midnight
        
        evMax = 3                  # Compute transit, rise and set times by default (1-3)
        cosH0 = (np.sin(rsa) - np.sin(st.geoLatitude) * np.sin(st._declinationUncorr)) / (np.cos(st.geoLatitude) *
                                                                                          np.cos(st._declinationUncorr))
        
        if(abs(cosH0) > 1.0):      # Body never rises/sets
            evMax = 1              # Compute transit time and altitude only
        else:
            h0 = np.arccos(cosH0) % self._PI  # Should probably work without %
        
        
        tmRad[0] = (st._rightAscensionUncorr - st.geoLongitude - st._agst) % self._TWOPI  # Transit time in radians; lon0 > 0 for E
        if(evMax > 1):
            tmRad[1] = (tmRad[0] - h0) % self._TWOPI   # Rise time in radians
            tmRad[2] = (tmRad[0] + h0) % self._TWOPI   # Set time in radians
        
        
        for evi in range(evMax):  # Loop over transit, rise, set
            iter = 0
            dTmRad = np.inf
            
            while(abs(dTmRad) > accur):
                th0 = agst0 + 1.002737909350795*tmRad[evi]   # Solar day in sidereal days in 2000
                
                st.second = tmRad[evi]*self._R2H*3600.0       # Radians -> seconds - w.r.t. midnight (h=0,m=0)
                
                # CHECK1: Replacing the line above with the one below, and using utc.to_julian_date() in
                # computePosition() and removal of self.year-second in setDateTime() is more elegant, but
                # slightly slower.  See CHECK1 in thise places.  HOWEVER, this also gives different results
                # for the rise/set times...(?)
                # st.utc = st.utc.normalize() + dt.timedelta(hours=tmRad[evi]*self._R2H)
                
                st.computePosition()
                
                ha  = self._revPI(th0 + st.geoLongitude - st._rightAscensionUncorr)        # Hour angle: -PI - +PI
                alt = np.arcsin(np.sin(st.geoLatitude)*np.sin(st._declinationUncorr) +
                                np.cos(st.geoLatitude)*np.cos(st._declinationUncorr)*np.cos(ha))  # Altitude
                
                # Correction to transit/rise/set times:
                if(evi==0):           # Transit
                    dTmRad = -self._revPI(ha)  # -PI - +PI
                else:                 # Rise/set
                    dTmRad = (alt-rsa)/(np.cos(st._declinationUncorr)*np.cos(st.geoLatitude)*np.sin(ha))
                    
                tmRad[evi] = tmRad[evi] + dTmRad
                
                # Print debug output to stdOut:
                # print(" %4i %2i %2i  %2i %2i %9.3lf    " % (st.year,st.month,st.day, st.hour,st.minute,st.second))
                # print(" %3i %4i   %9.3lf %9.3lf %9.3lf \n" % (evi,iter, tmRad[evi]*24,abs(dTmRad)*24,accur*24))
                
                iter += 1
                if(iter > 30): break  # The while loop doesn't seem to converge
            # end while(abs(dTmRad) > accur)
            
            
            if(iter > 30):  # Convergence failed
                print('\n  *** WARNING:  riset():  Riset failed to converge: %i %9.3lf  ***\n' % (evi,rsAlt))
                tmRad[evi] = -np.inf
                azalt[evi] = -np.inf
            else:               # Result converged, store it
                if(evi == 0):
                    azalt[evi] = alt                                                                      # Transit altitude
                else:
                    azalt[evi] = np.arctan2( np.sin(ha), ( np.cos(ha) * np.sin(st.geoLatitude)  -
                                                           np.tan(st._declinationUncorr) * np.cos(st.geoLatitude) ) )   # Rise,set hour angle -> azimuth
            
            
            if(tmRad[evi] < 0.0 and abs(rsAlt) < 1.e-9):
                tmRad[evi] = -np.inf
                azalt[evi] = -np.inf
                
        # end for loop evi
        
        
        # Convert times from radians to hours:
        tmHrs = tmRad*self._R2H
        
        # Set north to zero radians for azimuth if desired (use the original parameters!):
        if(self.param._useNorthEqualsZero):
            azalt[1] = (azalt[1] + self._PI) % self._TWOPI  # Add PI and fold between 0 and 2pi
            azalt[2] = (azalt[2] + self._PI) % self._TWOPI  # Add PI and fold between 0 and 2pi
        
        # Convert resulting angles to degrees if desired (use the original parameters!):
        if(self.param._useDegrees):
            azalt[0] *= self._R2D   # Transit altitude
            azalt[1] *= self._R2D   # Rise azimuth
            azalt[2] *= self._R2D   # Set azimuth
        
        if return_datetimes:  # Return datetimes iso time in hours:
            # tmHrs now contains [transit, rise, set] times in the local timezone.  Convert to datetime Series.
            jds = at.jd_from_date(midnight.year, midnight.month, midnight.day) + tmHrs/24
            tmDTs = at.date_time_from_jd(jds)  # Array of 6 arrays containing N x year, N x month, NxD, NxHour,NxMin,NxSec.
            
            # Combine the six date/time values into a single "2D" array with a single row and the original values as
            # columns, and convert it to a Pandas df:
            df = pd.DataFrame(np.vstack(tmDTs).transpose(),
                              columns=['year','month','day', 'hour','minute','second'])
            
            # Convert the date+time columns into a Pandas Series containing datetimes and round off to nearest second.
            tmDTI = pd.Series.copy(pd.to_datetime(df).dt.round('S'))  # Without .copy(), None below throws an error
            
            if evMax == 1:  # No rise or set, transit only
                tmDTI[1:3] = None
                azalt[1:3] = None
                
            return tmDTI, azalt
        
        
        if evMax == 1:  # No rise or set, transit only
            tmHrs[1:3] = None
            azalt[1:3] = None
            
        return tmHrs, azalt
