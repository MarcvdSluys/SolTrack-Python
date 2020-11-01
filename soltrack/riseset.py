
"""SolTrack: a simple, free, fast and accurate C routine to compute the position of the Sun.
  
  Copyright (c) 2014-2020  Marc van der Sluys, Paul van Kan and Jurgen Reintjes,
  Sustainable Energy research group, HAN University of applied sciences, Arnhem, The Netherlands
   
  This file is part of the SolTrack package, see: http://soltrack.sourceforge.net
  SolTrack is derived from libTheSky (http://libthesky.sourceforge.net) under the terms of the GPL v.3
  
  This is free software: you can redistribute it and/or modify it under the terms of the GNU Lesser General
  Public License as published by the Free Software Foundation, either version 3 of the License, or (at your
  option) any later version.
  
  This software is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the
  implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public
  License for more details.
  
  You should have received a copy of the GNU Lesser General Public License along with this code.  If not, see
  <http://www.gnu.org/licenses/>.

"""


from dataclasses import dataclass
import numpy as np
import soltrack

from .data import Constants


@dataclass
class RiseSet(Constants):
    """Class containing attributes and methods to compute the rise,transit and set times of the Sun and their azimuths/altitudes."""
    
    def __init__(self, param):
        self.param = param
    
    
    def computeSunRiseSet(self, location, time, rsAlt=0.0):
        """Compute rise, transit and set times for the Sun, as well as their azimuths/altitude.
        
        Parameters:
          location           (Location):  Class containing the geographic location to compute the Sun's rise and set times for.
          time               (Time):      Class containing date and time to compute the position for, in UT.
          rsAlt              (float):     Altitude to return rise/set data for (radians; optional, default=0.0 meaning actual rise/set).  Set rsAlt>pi/2 to compute transit only.
        
        Returns:
          (RiseSet):   Class containing the Sun's rise, transit and set data.
        
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
        
        # Need a local SolTrack instance for the same location (but in radians!), but with different settings
        # (radians, south=0, need equatorial coordinates but not the distance), and independent times and
        # positions:
        if(self.param.useDegrees):
            st = soltrack.SolTrack(location.longitude/self.R2D, location.latitude/self.R2D, useDegrees=False,
                                   useNorthEqualsZero=False, computeRefrEquatorial=True,
                                   computeDistance=False)
        else:
            st = soltrack.SolTrack(location.longitude, location.latitude, useDegrees=False,
                                   useNorthEqualsZero=False, computeRefrEquatorial=True,
                                   computeDistance=False)
        
        # Set date and time to midnight of the desired date:
        st.setTime(time.year, time.month, time.day, 0,0,0.0)
        
        # Compute the Sun's position:
        st.pos.computeSunPosition(st, st)
        
        
        # Compute transit, rise and set times:
        agst0 = st.pos.agst      # AGST for midnight
        
        evMax = 3                  # Compute transit, rise and set times by default (1-3)
        cosH0 = (np.sin(rsa)-np.sin(st.latitude)*np.sin(st.pos.declination)) / (np.cos(st.latitude)*np.cos(st.pos.declination))
        
        if(abs(cosH0) > 1.0):      # Body never rises/sets
            evMax = 1              # Compute transit time and altitude only
        else:
            h0 = np.arccos(cosH0) % self.PI  # Should probably work without %
            
        
        tmRad[0] = (st.pos.rightAscension - st.longitude - st.pos.agst) % self.TWO_PI  # Transit time in radians; lon0 > 0 for E
        if(evMax > 1):
            tmRad[1] = (tmRad[0] - h0) % self.TWO_PI   # Rise time in radians
            tmRad[2] = (tmRad[0] + h0) % self.TWO_PI   # Set time in radians
            
            
        accur = 1.0e-5        # Accuracy;  1e-5 rad ~ 0.14s. Don't make this smaller than 1e-16
        for evi in range(evMax):  # Loop over transit, rise, set
            iter = 0
            dTmRad = np.inf
            
            while(abs(dTmRad) > accur):
                th0 = agst0 + 1.002737909350795*tmRad[evi]  # Solar day in sidereal days in 2000
                
                st.second = tmRad[evi]*self.R2H*3600.0       # Radians -> seconds - w.r.t. midnight (h=0,m=0)
                st.pos.computeSunPosition(st, st)
                
                ha  = self.revPI(th0 + st.longitude - st.pos.rightAscension)        # Hour angle: -PI - +PI
                alt = np.arcsin(np.sin(st.latitude)*np.sin(st.pos.declination) +
                                np.cos(st.latitude)*np.cos(st.pos.declination)*np.cos(ha))  # Altitude
                
                # Correction to transit/rise/set times:
                if(evi==0):           # Transit
                    dTmRad = -self.revPI(ha)  # -PI - +PI
                else:                 # Rise/set
                    dTmRad = (alt-rsa)/(np.cos(st.pos.declination)*np.cos(st.latitude)*np.sin(ha))
                    
                tmRad[evi] = tmRad[evi] + dTmRad
                
                # Print debug output to stdOut:
                # print(" %4i %2i %2i  %2i %2i %9.3lf    " % (st.year,st.month,st.day, st.hour,st.minute,st.second))
                # print(" %3i %4i   %9.3lf %9.3lf %9.3lf \n" % (evi,iter, tmRad[evi]*24,abs(dTmRad)*24,accur*24))
                
                iter += 1
                if(iter > 30): break  # while loop doesn't seem to converge
            # while(abs(dTmRad) > accur)
            
            
            if(iter > 30):  # Convergence failed
                print("\n  *** WARNING:  riset():  Riset failed to converge: %i %9.3lf  ***\n" % (evi,rsAlt))
                tmRad[evi] = -np.inf
                azalt[evi] = -np.inf
            else:               # Result converged, store it
                if(evi == 0):
                    azalt[evi] = alt                                                                      # Transit altitude
                else:
                    azalt[evi] = np.arctan2( np.sin(ha), ( np.cos(ha) * np.sin(st.latitude)  -
                                                           np.tan(st.pos.declination) * np.cos(st.latitude) ) )   # Rise,set hour angle -> azimuth
            
            
            if(tmRad[evi] < 0.0 and abs(rsAlt) < 1.e-9):
                tmRad[evi] = -np.inf
                azalt[evi] = -np.inf
                
        # for-loop evi
        
        
        # Set north to zero radians for azimuth if desired (use the original parameters!):
        if(self.param.useNorthEqualsZero):
            azalt[1] = (azalt[1] + self.PI) % self.TWO_PI  # Add PI and fold between 0 and 2pi
            azalt[2] = (azalt[2] + self.PI) % self.TWO_PI  # Add PI and fold between 0 and 2pi
        
        
        # Convert resulting angles to degrees if desired (use the original parameters!):
        if(self.param.useDegrees):
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
    
    
    
    def revPI(self, angle):
        """Fold an angle in radians to take a value between -PI and +PI.
        
        Parameters:
          angle (float):  Angle to fold (rad).
        
        """
        
        return ((angle + self.PI) % self.TWO_PI) - self.PI
    

