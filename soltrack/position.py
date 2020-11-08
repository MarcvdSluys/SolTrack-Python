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

from .data import Constants, Parameters


@dataclass
class Position(Constants, Parameters):
    """Class containing the position of the Sun and related attributes and methods."""

    
    def __init__(self):
        Parameters.__init__(self)
        
    
    
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
    
    
    
    def _computeLongitude(self, computeDistance=False):
        """Compute the ecliptic longitude of the Sun for a given instant.
        
        Note:
          - Also computes the obliquity of the ecliptic and nutation.
        
        Parameters:
          computeDistance (bool):  Compute distance to the Sun (optional, default=False).
        
        """
        
        l0 = 4.895063168 + 628.331966786 * self._tJC  +  5.291838e-6 * self._tJC2  # Mean longitude
        ma = 6.240060141 + 628.301955152 * self._tJC  -  2.682571e-6 * self._tJC2  # Mean anomaly
        
        sec = (3.34161088e-2 - 8.40725e-5* self._tJC - 2.443e-7*self._tJC2)*np.sin(ma) + \
            (3.489437e-4 - 1.76278e-6*self._tJC)*np.sin(2*ma)                            # Sun's equation of the centre
        odot = l0 + sec                                                                 # True longitude
        
        
        # Nutation, aberration:
        omg  = 2.1824390725 - 33.7570464271 * self._tJC  + 3.622256e-5 * self._tJC2       # Lon. of Moon's mean ascending node
        dpsi = -8.338601e-5*np.sin(omg)                                                 # Nutation in longitude
        dist = 1.0000010178                                                             # Mean distance to the Sun in AU
        if(computeDistance):
            ecc = 0.016708634 - 0.000042037   * self._tJC  -  0.0000001267 * self._tJC2   # Eccentricity of the Earth's orbit
            nu = ma + sec                                                               # True anomaly
            dist = dist*(1.0 - ecc**2)/(1.0 + ecc*np.cos(nu))                           # Geocentric distance of the Sun in AU
            
        aber = -9.93087e-5/dist                                                         # Aberration
        
        # Obliquity of the ecliptic and nutation - do this here, since we've already computed many of the ingredients:
        eps0 = 0.409092804222 - 2.26965525e-4*self._tJC - 2.86e-9*self._tJC2              # Mean obliquity of the ecliptic
        deps = 4.4615e-5*np.cos(omg)                                                    # Nutation in obliquity
        
        # Save position parameters:
        self.longitude = (odot + aber + dpsi) % self._TWOPI                             # Apparent geocentric longitude, referred to the true equinox of date
        
        self.distance = dist                                                            # Distance (AU)
        
        self._obliquity   = eps0 + deps                                                  # True obliquity of the ecliptic
        self._cosObliquity = np.cos(self._obliquity)                                     # Need the cosine later on
        self._nutationLon = dpsi                                                         # Nutation in longitude
        
        return
    
    
    
    def _convertEclipticToEquatorial(self, longitude, cosObliquity):
        """Convert ecliptic coordinates to equatorial coordinates.
        
        Note:
          - This function assumes that the ecliptic latitude = 0.
        
        Parameters:
          longitude    (float):  Ecliptic longitude of the Sun (rad).
          cosObliquity (float):  Cosine of the obliquity of the ecliptic.
        
        Returns:
          tuple (float,float):  Tuple containing (rightAscension, declination):
        
          - rightAscension (float):  Right ascension of the Sun (rad).
          - declination    (float):  Declination of the Sun (rad).
        
        """
        
        sinLon = np.sin(longitude)
        sinObl = np.sqrt(1.0 - cosObliquity**2)               # Sine of the obliquity of the ecliptic will be positive in the forseeable future
        
        self.rightAscension   = np.arctan2(cosObliquity*sinLon, np.cos(longitude)) % self._TWOPI  # 0 <= azimuth < 2pi
        self.declination      = np.arcsin(sinObl*sinLon)
        
        return
    
    
    def _convertEquatorialToHorizontal(self):
        """Convert equatorial to horizontal coordinates.
        
        Also corrects for parallax and atmospheric refraction.
        
        """
        
        gmst       = 4.89496121 + 6.300388098985*self._tJD + 6.77e-6*self._tJC2  # Greenwich mean sidereal time
        self._agst = gmst + self._nutationLon * self._cosObliquity              # Correction for equation of the equinoxes . apparent Greenwich sidereal time
        
        
        sinAlt=0.0
        # Azimuth does not need to be corrected for parallax or refraction, hence store the result in the 'azimuthRefract' variable directly:
        self.azimuthRefract, sinAlt = self._eq2horiz(self._sinLat,self._cosLat, self.geoLongitude,
                                                     self.rightAscension, self.declination, self._agst)
        
        alt = np.arcsin( sinAlt )                                  # Altitude of the Sun above the horizon (rad)
        cosAlt = np.sqrt(1.0 - sinAlt**2)                          # Cosine of the altitude is always positive or zero
        
        # Correct for parallax:
        alt -= 4.2635e-5 * cosAlt                                  # Horizontal parallax = 8.794" = 4.2635e-5 rad
        self._altitude = alt
        
        # Correct for atmospheric refraction:
        dalt = 2.967e-4 / np.tan(alt + 3.1376e-3/(alt + 8.92e-2))  # Refraction correction in altitude
        dalt *= self.pressure/101.0 * 283.0/self.temperature
        alt += dalt
        # to do: add pressure/temperature dependence
        self.altitudeRefract = alt
        
        return
    
    
    def _eq2horiz(self, sinLat, cosLat, longitude,  rightAscension, declination, agst):
        """Convert equatorial coordinates to horizontal coordinates.
        
        Parameters:
          sinLat         (float):  Sine of the geographic latitude of the observer.
          cosLat         (float):  Cosine of the geographic latitude of the observer.
          longitude      (float):  Geographic longitude of the observer (rad).
          rightAscension (float):  Right ascension of the Sun (rad).
          declination    (float):  Declination of the Sun (rad).
          agst           (float):  Apparent Greenwich sidereal time (Greenwich mean sidereal time, corrected for the equation of the equinoxes).
        
        Returns: 
          tuple (float,float):  Tuple containing (azimuth, sinAlt):
        
            - azimuth (float):  Azimuth ("wind direction") of the Sun (rad; 0=South).
            - sinAlt  (float):  Sine of the altitude of the Sun above the horizon.
        
        """
        
        ha  = agst + longitude - rightAscension                      # Local Hour Angle
        
        # Some preparation, saves ~29%:
        sinHa  = np.sin(ha)
        cosHa  = np.cos(ha)
        
        sinDec = np.sin(declination)
        cosDec = np.sqrt(1.0 - sinDec**2)                         # Cosine of a declination is always positive or zero
        tanDec = sinDec/cosDec
        
        azimuth = np.arctan2( sinHa,  cosHa  * sinLat - tanDec * cosLat )  # 0 <= azimuth < 2pi
        sinAlt = sinLat * sinDec + cosLat * cosDec * cosHa                 # Sine of the altitude above the horizon
        
        return azimuth, sinAlt
    
    
    def _convertHorizontalToEquatorial(self, sinLat, cosLat, azimuth, altitude):
        """Convert (refraction-corrected) horizontal coordinates to the equatorial coordinates hourAngle and
        declination.
        
        Parameters:
          sinLat   (float):  Sine of the geographic latitude of the observer.
          cosLat   (float):  Cosine of the geographic latitude of the observer.
          azimuth  (float):  Azimuth ("wind direction") of the Sun (rad; 0=South).
          altitude (float):  Altitude of the Sun above the horizon (rad).
        
        """
        
        # Multiply used variables:
        cosAz  = np.cos(azimuth)
        sinAz  = np.sin(azimuth)                                      # For symmetry
        
        sinAlt = np.sin(altitude)
        cosAlt = np.sqrt(1.0 - sinAlt**2)                             # Cosine of an altitude is always positive or zero
        tanAlt = sinAlt/cosAlt
        
        self.hourAngleRefract   = np.arctan2( sinAz,   cosAz  * sinLat + tanAlt * cosLat )      # Local Hour Angle:  0 <= hourAngle < 2pi
        self.declinationRefract = np.arcsin(  sinLat * sinAlt  -  cosLat * cosAlt * cosAz  )    # Declination
        
        return
    
    
    
    def _setNorthToZero(self, azimuthRefract, hourAngleRefract):
        """Convert the South=0 convention to North=0 convention for azimuth and hour angle.
        
        Note:
          - South=0 is the default in celestial astronomy.
          - This function makes the angles compatible with the compass/wind directions.
        
        Parameters:
          azimuth              (float):  Azimuth ("wind direction") of the Sun (rad).
          hourAngle            (float):  Hour angle of the Sun (rad).
        
        """
        
        self.azimuthRefract = (azimuthRefract + self._PI) % self._TWOPI                    # Add PI to set North=0
        
        if(self.param._computeRefrEquatorial):
            self.hourAngleRefract = (hourAngleRefract + self._PI) % self._TWOPI            # Add PI to set North=0
            
        return
    
    
    def _convertRadiansToDegrees(self):
        """Convert final results from radians to degrees.
        
        Note:
          - Does not touch intermediate results.
        
        """
        
        self.longitude *= self._R2D
        self.rightAscension *= self._R2D
        self.declination *= self._R2D
        
        self._altitude *= self._R2D
        self.azimuthRefract *= self._R2D
        self.altitudeRefract *= self._R2D
        
        if(self.param._computeRefrEquatorial):
            self.hourAngleRefract *= self._R2D
            self.declinationRefract *= self._R2D
            
        return
    
    
    def _revPI(self, angle):
        """Fold an angle in radians to take a value between -PI and +PI.
        
        Parameters:
          angle (float):  Angle to fold (radians).
        
        Returns:
          float: Angle between -PI and PI (radians).
        
        """
        
        return ((angle + self._PI) % self._TWOPI) - self._PI
    

