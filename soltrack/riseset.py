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
import astrotool.date_time as at_dt
from astroconst import r2d as _R2D, pi as _PI, pi2 as _TWOPI, r2h as _R2H

from .data import Parameters


@dataclass
class RiseSet(Parameters):
    """Class concerning the rise, transit and set times and positions of the Sun and related attributes and methods."""
    
    
    def __init__(self):
        Parameters.__init__(self)
        
        # Rise, transit and set time:
        self.rise_time:               float = None;     """Rise time of the Sun (hours LT or UTC)"""
        self.transit_time:            float = None;     """Transit time of the Sun (hours LT or UTC)"""
        self.set_time:                float = None;     """Set time of the Sun (hours LT or UTC)"""
        
        # Rise, transit and set position:
        self.rise_azimuth:            float = None;     """Rise azimuth of the Sun (radians)"""
        self.transit_altitude:        float = None;     """Transit altitude of the Sun (radians)"""
        self.set_azimuth:             float = None;     """Set azimuth of the Sun (radians)"""
        
    
    
    def computeRiseSet(self, rs_alt=0.0, accur=1e-5, return_datetimes=True):
        """This function is obsolescent and will be removed in a future version.  Use compute_rise_set()
        instead."""
        _warn_obsolescent('computeRiseSet', 'compute_rise_set', rename=True)
        return self.compute_rise_set(rs_alt, accur, return_datetimes)
    
    
    def compute_rise_set(self, rs_alt=0.0, accur=1e-5, return_datetimes=True):
        """Compute rise, transit and set times for the Sun, as well as their azimuths/altitude.
        
        Parameters:
          rs_alt           (float):  Altitude to return rise/set data for (radians; optional, default=0.0 meaning actual rise/set).  Set rs_alt>pi/2 to compute transit only.
          accur:           (float):  Accuracy (rad).  Default: 1e-5 rad ~ 0.14s.  Don't make this smaller than 1e-16.
          return_datetimes (bool):   Return times as datetimes rather than decimal hours.  Defaults to True.
        
        Note:
          - rise/set/transit times are in the LOCAL timezone used for the input (hence UTC if UTC was used).
          - if rs_alt == 0.0, actual rise and set times are computed
          - if rs_alt != 0.0, the routine calculates when alt = rs_alt is reached
          - returns times, rise/set azimuth and transit altitude in the class riseSet
         
          See:
            - subroutine riset() in riset.f90 from libTheSky (libthesky.sf.net) for more info
        """
        
        trTime = np.empty(0);  riTime = np.empty(0);  seTime = np.empty(0)
        trAlt  = np.empty(0);  riAz   = np.empty(0);  seAz   = np.empty(0)
        
        # CHECK: we should use local time for local Sun rise/transit/set, but what do we do if UTC is available only?
        if self.lt is not None:
            orig_dt = self.lt
        else:
            orig_dt = self.utc  # Try to correct for tz using geoLon?  Important is that 'midnight' lies between sunset and sunrise (CHECK: true? sufficient?)
        
        
        for OrigDTi in orig_dt:  # Loop over dates in array
            tmHrs, azalt = self._compute_rise_set_single(OrigDTi, rs_alt, accur, return_datetimes)  # Compute r/s/t for a single date
            
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
        
        # If a single time was used, convert arrays to scalars:
        if orig_dt.size == 1:
            if return_datetimes:
                self.transit_time      = pd.Timestamp(tmHrs[0])  # Array -> scalar
                self.rise_time         = pd.Timestamp(tmHrs[1])  # Array -> scalar
                self.set_time          = pd.Timestamp(tmHrs[2])  # Array -> scalar
            else:
                self.transit_time      = tmHrs[0]                # Array -> scalar
                self.rise_time         = tmHrs[1]                # Array -> scalar
                self.set_time          = tmHrs[2]                # Array -> scalar
            
            self.transit_altitude  = azalt[0]                # Array -> scalar, float
            self.rise_azimuth      = azalt[1]                # Array -> scalar, float
            self.set_azimuth       = azalt[2]                # Array -> scalar, float
            
        else:
            
            # Store final arrays:
            self.transit_time      = trTime
            self.rise_time         = riTime
            self.set_time          = seTime
            
            self.transit_altitude  = trAlt
            self.rise_azimuth      = riAz
            self.set_azimuth       = seAz
        
        return
    
    
    def _computeRiseSetSingle(self, orig_dt, rs_alt, accur, return_datetimes):
        """This function is obsolescent and will be removed in a future version.  Use compute_rise_set()
        instead."""
        _warn_obsolescent('_computeRiseSetSingle', '_compute_rise_set_single', rename=True)
        return self._compute_rise_set_single(orig_dt, rs_alt, accur, return_datetimes)
    
        
    def _compute_rise_set_single(self, orig_dt, rs_alt, accur, return_datetimes):
        """Compute rise, transit and set times for the Sun, as well as their azimuths/altitude for a single date.
        
        Parameters:
          orig_dt         (datetime):  Original datetime of the date to compute the rise, set and transit for..
          rs_alt             (float):  Altitude to return rise/set data for (radians).  Set rs_alt>pi/2 to compute transit only.
          accur:             (float):  Accuracy (rad).  Don't make this smaller than 1e-16.
          return_datetimes    (bool):  Return times as datetimes rather than decimal hours.
        
        See computeRiseSet() for more details.
        """
        
        rsa = -0.8333/_R2D                    # Standard altitude for the Sun in radians
        if abs(rs_alt) > 1.e-9: rsa = rs_alt  # Use a user-specified altitude
        
        tmRad = np.zeros(3)
        azalt = np.zeros(3)
        alt=0.0;  ha=0.0; h0=0.0
        
        
        # We need a local SolTrack instance for the same location (but in radians!), but with different settings
        # (radians, south=0, need equatorial coordinates but not the distance), and independent times and
        # positions:
        if self.param._use_degrees:
            st = self.__class__(self.geo_longitude/_R2D, self.geo_latitude/_R2D, use_degrees=False,
                                use_north_equals_zero=False, compute_refr_equatorial=True, compute_distance=False)
        else:
            st = self.__class__(self.geo_longitude, self.geo_latitude, use_degrees=False, use_north_equals_zero=False,
                                compute_refr_equatorial=True, compute_distance=False)
        
        # Set date and time to midnight of the desired date:
        midnight = orig_dt.normalize()  # Midnight for the date in OrigDT, keeping the timezone
        st.set_date_time(midnight)
        
        # Compute the Sun's position:
        st.compute_position()
        
        
        # Compute transit, rise and set times:
        agst0 = st._agst      # AGST for midnight
        
        evMax = 3                  # Compute transit, rise and set times by default (1-3)
        cosH0 = (np.sin(rsa) - np.sin(st.geo_latitude) * np.sin(st._declination_uncorr)) / (np.cos(st.geo_latitude) *
                                                                                            np.cos(st._declination_uncorr))
        
        if abs(cosH0) > 1.0:       # Body never rises/sets
            evMax = 1              # Compute transit time and altitude only
        else:
            h0 = np.arccos(cosH0) % _PI  # Should probably work without %
        
        
        tmRad[0] = (st._right_ascension_uncorr - st.geo_longitude - st._agst) % _TWOPI  # Transit time in radians; lon0 > 0 for E
        if evMax > 1:
            tmRad[1] = (tmRad[0] - h0) % _TWOPI   # Rise time in radians
            tmRad[2] = (tmRad[0] + h0) % _TWOPI   # Set time in radians
        
        
        for evi in range(evMax):  # Loop over transit, rise, set
            iter = 0
            dTmRad = np.inf
            
            while abs(dTmRad) > accur:
                th0 = agst0 + 1.002737909350795*tmRad[evi]  # Solar day in sidereal days in 2000
                
                st.second = tmRad[evi]*_R2H*3600.0          # Radians -> seconds - w.r.t. midnight (h=0,m=0)
                
                # CHECK1: Replacing the line above with the one below, and using utc.to_julian_date() in
                # compute_position() and removal of self.year-second in set_date_time() is more elegant, but
                # slightly slower.  See CHECK1 in thise places.  HOWEVER, this also gives different results
                # for the rise/set times...(?)
                # st.utc = st.utc.normalize() + dt.timedelta(hours=tmRad[evi]*_R2H)
                
                st.compute_position()
                
                ha  = self._rev_pi(th0 + st.geo_longitude - st._right_ascension_uncorr)        # Hour angle: -PI - +PI
                alt = np.arcsin(np.sin(st.geo_latitude)*np.sin(st._declination_uncorr) +
                                np.cos(st.geo_latitude)*np.cos(st._declination_uncorr)*np.cos(ha))  # Altitude
                
                # Correction to transit/rise/set times:
                if evi==0:            # Transit
                    dTmRad = -self._rev_pi(ha)  # -PI - +PI
                else:                 # Rise/set
                    dTmRad = (alt-rsa)/(np.cos(st._declination_uncorr)*np.cos(st.geo_latitude)*np.sin(ha))
                    
                tmRad[evi] = tmRad[evi] + dTmRad
                
                # Print debug output to stdOut:
                # print(" %4i %2i %2i  %2i %2i %9.3lf    " % (st.year,st.month,st.day, st.hour,st.minute,st.second))
                # print(" %3i %4i   %9.3lf %9.3lf %9.3lf \n" % (evi,iter, tmRad[evi]*24,abs(dTmRad)*24,accur*24))
                
                iter += 1
                if iter > 30: break  # The while loop doesn't seem to converge
            # end while(abs(dTmRad) > accur)
            
            
            if iter > 30:  # Convergence failed
                print('\n  *** WARNING:  riset():  Riset failed to converge: %i %9.3lf  ***\n' % (evi,rs_alt))
                tmRad[evi] = -np.inf
                azalt[evi] = -np.inf
            else:               # Result converged, store it
                if evi == 0:
                    azalt[evi] = alt                                                                      # Transit altitude
                else:
                    azalt[evi] = np.arctan2( np.sin(ha), ( np.cos(ha) * np.sin(st.geo_latitude)  -
                                                           np.tan(st._declination_uncorr) * np.cos(st.geo_latitude) ) )   # Rise,set hour angle -> azimuth
            
            
            if tmRad[evi] < 0.0 and abs(rs_alt) < 1.e-9:
                tmRad[evi] = -np.inf
                azalt[evi] = -np.inf
                
        # end for loop evi
        
        
        # Convert times from radians to hours:
        tmHrs = tmRad*_R2H
        
        # Set north to zero radians for azimuth if desired (use the original parameters!):
        if self.param._use_north_equals_zero:
            azalt[1] = (azalt[1] + _PI) % _TWOPI  # Add PI and fold between 0 and 2pi
            azalt[2] = (azalt[2] + _PI) % _TWOPI  # Add PI and fold between 0 and 2pi
        
        # Convert resulting angles to degrees if desired (use the original parameters!):
        if self.param._use_degrees:
            azalt[0] *= _R2D   # Transit altitude
            azalt[1] *= _R2D   # Rise azimuth
            azalt[2] *= _R2D   # Set azimuth
        
        if return_datetimes:  # Return datetimes iso time in hours:
            # tmHrs now contains [transit, rise, set] times in the local timezone.  Convert to datetime Series.
            jds = at_dt.jd_from_date(midnight.year, midnight.month, midnight.day) + tmHrs/24
            tmDTs = at_dt.date_time_from_jd(jds)  # Array of 6 arrays containing N x year, N x month, NxD, NxHour,NxMin,NxSec.
            
            # Combine the six date/time values into a single "2D" array with a single row and the original values as
            # columns, and convert it to a Pandas df:
            df = pd.DataFrame(np.vstack(tmDTs).transpose(),
                              columns=['year','month','day', 'hour','minute','second'])
            
            # Convert the date+time columns into a Pandas Series containing datetimes and round off to nearest second.
            tmDTI = pd.Series.copy(pd.to_datetime(df).dt.round('s'))  # Without .copy(), None below throws an error
            
            if evMax == 1:  # No rise or set, transit only
                tmDTI[1:3] = None
                azalt[1:3] = None
                
            return tmDTI, azalt
        
        
        if evMax == 1:  # No rise or set, transit only
            tmHrs[1:3] = None
            azalt[1:3] = None
            
        return tmHrs, azalt
    
    
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
