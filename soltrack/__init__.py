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


"""SolTrack module

SolTrack is a simple, free, fast and accurate Python package to compute the position of the Sun, as well as
its rise and set times.  SolTrack can be used under the conditions of the GPLv3 licence.  These pages contain
the API documentation.  For more information on the Python package, licence, source code and data files, see
the [SolTrack homepage](http://soltrack.sf.net).

"""

name = "soltrack"

from dataclasses import dataclass
import numpy as np

from .data      import Parameters
from .location  import Location
from .time      import Time
from .position  import Position


@dataclass
class SolTrack(Location, Time, Position):
    
    
    def __init__(self, geoLongitude,geoLatitude, useDegrees=None, useNorthEqualsZero=None, computeRefrEquatorial=None, computeDistance=None):
        
        """Construct a SolTrack object with specified geographical location and parameters (settings).
        
        Parameters:
          geoLongitude (float):           Geographical longitude of the observer or installation (radians or degrees, depending on useDegrees).
          geoLatitude (float):            Geographical latitude of the observer or installation (radians or degrees, depending on useDegrees).
        
          useDegrees (bool):              Input (geographic position) and output are in degrees, rather than radians.
          useNorthEqualsZero (bool):      Azimuth: 0 = South, pi/2 (90deg) = West  ->  0 = North, pi/2 (90deg) = East.
          computeRefrEquatorial (bool):   Compute refraction-corrected equatorial coordinates (Hour angle, declination).
          computeDistance (bool):         Compute the distance to the Sun.
        
        """
        
        self.param     = Parameters()
        self.param.setParameters(useDegrees, useNorthEqualsZero, computeRefrEquatorial, computeDistance)
        
        Location.__init__(self, geoLongitude, geoLatitude)
        
        Time.__init__(self)
        
        Position.__init__(self)
        
        
        
        
    def computePosition(self):
        
        """ Method to compute the position of the Sun.
        """
                
        # If the user uses degrees, convert the geographic location to radians:
        if(self.param._useDegrees):
            self.geoLongitude /= self.R2D
            self.geoLatitude  /= self.R2D
        
        # Compute these once and reuse:
        self._sinLat = np.sin(self.geoLatitude)
        self._cosLat = np.sqrt(1.0 - self._sinLat**2)  # Cosine of a latitude is always positive or zero
        
        
        # Compute the Julian Day from the date and time:
        self._computeJulianDay(self.year, self.month, self.day, self.hour, self.minute, self.second)
        
        # Derived expressions for time, to be reused:
        self._tJD  = self.julianDay - 2451545.0                   # Time in Julian days since 2000.0
        self._tJC  = self._tJD/36525.0                             # Time in Julian centuries since 2000.0
        self._tJC2 = self._tJC**2                                  # T^2
        
        
        # Compute the ecliptic longitude of the Sun and the obliquity of the ecliptic:
        self._computeLongitude(self.param._computeDistance)
        
        # Convert ecliptic coordinates to geocentric equatorial coordinates:
        self._convertEclipticToEquatorial(self.longitude, self._cosObliquity)
        
        # Convert equatorial coordinates to horizontal coordinates, correcting for parallax and refraction:
        self._convertEquatorialToHorizontal()
        
        
        # Convert the corrected horizontal coordinates back to equatorial coordinates:
        if(self.param._computeRefrEquatorial):
            self._convertHorizontalToEquatorial(self._sinLat, self._cosLat, self.azimuthRefract,
                                                self.altitudeRefract)
            
        # Use the North=0 convention for azimuth and hour angle (default: South = 0) if desired:
        if(self.param._useNorthEqualsZero):
            self._setNorthToZero(self.azimuthRefract, self.hourAngleRefract)
            
        # If the user wants degrees, convert final results from radians to degrees:
        if(self.param._useDegrees):
            self.geoLongitude *= self.R2D  # Convert back to original
            self.geoLatitude  *= self.R2D  # Convert back to original
            self._convertRadiansToDegrees()    # Convert final results
        
        return
    
    
    
    def computeRiseSet(self, rsAlt=0.0):
        """Compute rise, transit and set times for the Sun, as well as their azimuths/altitude.
        
        Parameters:
          rsAlt              (float):     Altitude to return rise/set data for (radians; optional, default=0.0 meaning actual rise/set).  Set rsAlt>pi/2 to compute transit only.
        
        Note:
          - if rsAlt == 0.0, actual rise and set times are computed
          - if rsAlt != 0.0, the routine calculates when alt = rsAlt is reached
          - returns times, rise/set azimuth and transit altitude in the class riseSet
         
          See:
            - subroutine riset() in riset.f90 from libTheSky (libthesky.sf.net) for more info
        
        """
        
        rsa = -0.8333/self.R2D              # Standard altitude for the Sun in radians
        if(abs(rsAlt) > 1.e-9): rsa = rsAlt     # Use a user-specified altitude
        
        tmRad = np.zeros(3)
        azalt = np.zeros(3)
        alt=0.0;  ha=0.0; h0=0.0
        
        # We need a local SolTrack instance for the same location (but in radians!), but with different settings
        # (radians, south=0, need equatorial coordinates but not the distance), and independent times and
        # positions:
        if(self.param._useDegrees):
            st = SolTrack(self.geoLongitude/self.R2D, self.geoLatitude/self.R2D, useDegrees=False,
                          useNorthEqualsZero=False, computeRefrEquatorial=True, computeDistance=False)
        else:
            st = SolTrack(self.geoLongitude, self.geoLatitude, useDegrees=False, useNorthEqualsZero=False,
                          computeRefrEquatorial=True, computeDistance=False)
        
        # Set date and time to midnight of the desired date:
        st.setDateTime(self.year, self.month, self.day, 0,0,0.0)
        
        # Compute the Sun's position:
        st.computePosition()
        
        
        # Compute transit, rise and set times:
        agst0 = st._agst      # AGST for midnight
        
        evMax = 3                  # Compute transit, rise and set times by default (1-3)
        cosH0 = (np.sin(rsa) - np.sin(st.geoLatitude) * np.sin(st.declination)) / (np.cos(st.geoLatitude) *
                                                                                   np.cos(st.declination))
        
        if(abs(cosH0) > 1.0):      # Body never rises/sets
            evMax = 1              # Compute transit time and altitude only
        else:
            h0 = np.arccos(cosH0) % self.PI  # Should probably work without %
            
        
        tmRad[0] = (st.rightAscension - st.geoLongitude - st._agst) % self.TWO_PI  # Transit time in radians; lon0 > 0 for E
        if(evMax > 1):
            tmRad[1] = (tmRad[0] - h0) % self.TWO_PI   # Rise time in radians
            tmRad[2] = (tmRad[0] + h0) % self.TWO_PI   # Set time in radians
            
            
        accur = 1.0e-5            # Accuracy;  1e-5 rad ~ 0.14s. Don't make this smaller than 1e-16
        for evi in range(evMax):  # Loop over transit, rise, set
            iter = 0
            dTmRad = np.inf
            
            while(abs(dTmRad) > accur):
                th0 = agst0 + 1.002737909350795*tmRad[evi]   # Solar day in sidereal days in 2000
                
                st.second = tmRad[evi]*self.R2H*3600.0       # Radians -> seconds - w.r.t. midnight (h=0,m=0)
                st.computePosition()
                
                ha  = self._revPI(th0 + st.geoLongitude - st.rightAscension)        # Hour angle: -PI - +PI
                alt = np.arcsin(np.sin(st.geoLatitude)*np.sin(st.declination) +
                                np.cos(st.geoLatitude)*np.cos(st.declination)*np.cos(ha))  # Altitude
                
                # Correction to transit/rise/set times:
                if(evi==0):           # Transit
                    dTmRad = -self._revPI(ha)  # -PI - +PI
                else:                 # Rise/set
                    dTmRad = (alt-rsa)/(np.cos(st.declination)*np.cos(st.geoLatitude)*np.sin(ha))
                    
                tmRad[evi] = tmRad[evi] + dTmRad
                
                # Print debug output to stdOut:
                # print(" %4i %2i %2i  %2i %2i %9.3lf    " % (st.year,st.month,st.day, st.hour,st.minute,st.second))
                # print(" %3i %4i   %9.3lf %9.3lf %9.3lf \n" % (evi,iter, tmRad[evi]*24,abs(dTmRad)*24,accur*24))
                
                iter += 1
                if(iter > 30): break  # The while loop doesn't seem to converge
            # while(abs(dTmRad) > accur)
            
            
            if(iter > 30):  # Convergence failed
                print("\n  *** WARNING:  riset():  Riset failed to converge: %i %9.3lf  ***\n" % (evi,rsAlt))
                tmRad[evi] = -np.inf
                azalt[evi] = -np.inf
            else:               # Result converged, store it
                if(evi == 0):
                    azalt[evi] = alt                                                                      # Transit altitude
                else:
                    azalt[evi] = np.arctan2( np.sin(ha), ( np.cos(ha) * np.sin(st.geoLatitude)  -
                                                           np.tan(st.declination) * np.cos(st.geoLatitude) ) )   # Rise,set hour angle -> azimuth
            
            
            if(tmRad[evi] < 0.0 and abs(rsAlt) < 1.e-9):
                tmRad[evi] = -np.inf
                azalt[evi] = -np.inf
                
        # for-loop evi
        
        
        # Set north to zero radians for azimuth if desired (use the original parameters!):
        if(self.param._useNorthEqualsZero):
            azalt[1] = (azalt[1] + self.PI) % self.TWO_PI  # Add PI and fold between 0 and 2pi
            azalt[2] = (azalt[2] + self.PI) % self.TWO_PI  # Add PI and fold between 0 and 2pi
        
        
        # Convert resulting angles to degrees if desired (use the original parameters!):
        if(self.param._useDegrees):
            azalt[0] *= self.R2D   # Transit altitude
            azalt[1] *= self.R2D   # Rise azimuth
            azalt[2] *= self.R2D   # Set azimuth
            
            
        # Store results:
        self.transitTime     = tmRad[0]*self.R2H  # Transit time - radians -> hours
        self.riseTime        = tmRad[1]*self.R2H  # Rise time - radians -> hours
        self.setTime         = tmRad[2]*self.R2H  # Set time - radians -> hours
        
        self.transitAltitude = azalt[0]      # Transit altitude
        self.riseAzimuth     = azalt[1]      # Rise azimuth
        self.setAzimuth      = azalt[2]      # Set azimuth
        
        return
    
    
    
