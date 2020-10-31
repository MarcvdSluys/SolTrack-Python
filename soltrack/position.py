
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

from dataclasses import dataclass
import numpy as np

from .data import Data
data = Data()


@dataclass
class Position:
    """Class containing the position of the Sun and related variables."""
    
    def computeSunPosition(self, location, time,  useDegrees=False, useNorthEqualsZero=False, computeRefrEquatorial=False, computeDistance=False):
        
        """
        Main function to compute the position of the Sun.
        
        Parameters:
          location          (Location):  Class containing the geographic location to compute the Sun's position for.
          time                  (Time):  Class containing date and time to compute the position for, in UT.
        
          useDegrees            (bool):  Use degrees for input and output angular variables, rather than radians (optional, default=False).
          useNorthEqualsZero    (bool):  Use the definition where azimuth=0 denotes north, rather than south (optional, default=False).
          computeRefrEquatorial (bool):  Compute refraction correction for equatorial coordinates (optional, default=False).
          computeDistance       (bool):  Compute distance to the Sun (optional, default=False).
        
        Returns:
          (Position):  Class containing the position of the Sun in horizontal (and equatorial if desired) coordinates (output).
        
        """
        
        # If the used uses degrees, convert the geographic location to radians:
        # In C, a local copy of location is made.  With Python objects, we need a deep copy.
        import copy
        loc = copy.deepcopy(location)  # Local instance of the Location class, so that it can be changed here
        if(useDegrees):
            loc.longitude /= data.R2D
            loc.latitude  /= data.R2D
        
        # Compute these once and reuse:
        loc.sinLat = np.sin(loc.latitude)
        loc.cosLat = np.sqrt(1.0 - loc.sinLat**2)  # Cosine of a latitude is always positive or zero
        
        
        # Compute the Julian Day from the date and time:
        self.computeJulianDay(time.year, time.month, time.day, time.hour, time.minute, time.second)
        
        
        # Derived expressions of time:
        self.tJD  = self.julianDay - 2451545.0                   # Time in Julian days since 2000.0
        self.tJC  = self.tJD/36525.0                             # Time in Julian centuries since 2000.0
        self.tJC2 = self.tJC**2                                  # T^2
        
        
        # Compute the ecliptic longitude of the Sun and the obliquity of the ecliptic:
        self.computeLongitude(computeDistance)
        
        # Convert ecliptic coordinates to geocentric equatorial coordinates:
        self.convertEclipticToEquatorial(self.longitude, self.cosObliquity)
        
        # Convert equatorial coordinates to horizontal coordinates, correcting for parallax and refraction:
        self.convertEquatorialToHorizontal(loc)
        
        
        # Convert the corrected horizontal coordinates back to equatorial coordinates:
        if(computeRefrEquatorial):
            self.convertHorizontalToEquatorial(loc.sinLat, loc.cosLat, self.azimuthRefract,
                                               self.altitudeRefract)
            
        # Use the North=0 convention for azimuth and hour angle (default: South = 0) if desired:
        if(useNorthEqualsZero):
            self.setNorthToZero(self.azimuthRefract, self.hourAngleRefract, computeRefrEquatorial)
            
        # If the user wants degrees, convert final results from radians to degrees:
        if(useDegrees):
            self.convertRadiansToDegrees(computeRefrEquatorial)
            
        return
    
    
    
    def computeJulianDay(self, year, month, day,  hour, minute, second):
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
    
    
    
    def computeLongitude(self, computeDistance=False):
        """Compute the ecliptic longitude of the Sun for a given instant.
        
        Note:
          - Also computes the obliquity of the ecliptic and nutation.
        
        Parameters:
          position    (Position):  Class containing the position of the Sun (I/O).
          computeDistance (bool):  Compute distance to the Sun (optional, default=False).
        
        """
        
        l0 = 4.895063168 + 628.331966786 * self.tJC  +  5.291838e-6 * self.tJC2  # Mean longitude
        ma = 6.240060141 + 628.301955152 * self.tJC  -  2.682571e-6 * self.tJC2  # Mean anomaly
        
        sec = (3.34161088e-2 - 8.40725e-5* self.tJC - 2.443e-7*self.tJC2)*np.sin(ma) + \
            (3.489437e-4 - 1.76278e-6*self.tJC)*np.sin(2*ma)                               # Sun's equation of the centre
        odot = l0 + sec                                                                    # True longitude
        
        
        # Nutation, aberration:
        omg  = 2.1824390725 - 33.7570464271 * self.tJC  + 3.622256e-5 * self.tJC2          # Lon. of Moon's mean ascending node
        dpsi = -8.338601e-5*np.sin(omg)                                                    # Nutation in longitude
        dist = 1.0000010178                                                                # Mean distance to the Sun in AU
        if(computeDistance):
            ecc = 0.016708634 - 0.000042037   * self.tJC  -  0.0000001267 * self.tJC2      # Eccentricity of the Earth's orbit
            nu = ma + sec                                                                  # True anomaly
            dist = dist*(1.0 - ecc**2)/(1.0 + ecc*np.cos(nu))                             # Geocentric distance of the Sun in AU
            
        aber = -9.93087e-5/dist                                                            # Aberration
        
        # Obliquity of the ecliptic and nutation - do this here, since we've already computed many of the ingredients:
        eps0 = 0.409092804222 - 2.26965525e-4*self.tJC - 2.86e-9*self.tJC2                 # Mean obliquity of the ecliptic
        deps = 4.4615e-5*np.cos(omg)                                                       # Nutation in obliquity
        
        # Save position parameters:
        self.longitude = (odot + aber + dpsi) % data.TWO_PI                                # Apparent geocentric longitude, referred to the true equinox of date
        
        self.distance = dist                                                               # Distance (AU)
        
        self.obliquity   = eps0 + deps                                                     # True obliquity of the ecliptic
        self.cosObliquity = np.cos(self.obliquity)                                         # Need the cosine later on
        self.nutationLon = dpsi                                                            # Nutation in longitude
        
        return
    
    
    
    def convertEclipticToEquatorial(self, longitude, cosObliquity):
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
        sinObl = np.sqrt(1.0-cosObliquity**2)               # Sine of the obliquity of the ecliptic will be positive in the forseeable future
        
        self.rightAscension   = np.arctan2(cosObliquity*sinLon, np.cos(longitude)) % data.TWO_PI  # 0 <= azimuth < 2pi
        self.declination      = np.arcsin(sinObl*sinLon)
        
        return
    
    
    def convertEquatorialToHorizontal(self, location):
        """Convert equatorial to horizontal coordinates.
        
        Also corrects for parallax and atmospheric refraction.
        
        Parameters:
          location (Location):  Class containing the geographic location of the observer (rad).
          position (Position):  Class containing the position of the Sun (rad, I/O).
        
        """
        
        gmst      = 4.89496121 + 6.300388098985*self.tJD + 6.77e-6*self.tJC2  # Greenwich mean sidereal time
        self.agst = gmst + self.nutationLon * self.cosObliquity               # Correction for equation of the equinoxes . apparent Greenwich sidereal time
        
        
        sinAlt=0.0
        # Azimuth does not need to be corrected for parallax or refraction, hence store the result in the 'azimuthRefract' variable directly:
        self.azimuthRefract, sinAlt = self.eq2horiz(location.sinLat,location.cosLat, location.longitude,
                                                    self.rightAscension, self.declination, self.agst)
        
        alt = np.arcsin( sinAlt )                                  # Altitude of the Sun above the horizon (rad)
        cosAlt = np.sqrt(1.0 - sinAlt**2)                          # Cosine of the altitude is always positive or zero
        
        # Correct for parallax:
        alt -= 4.2635e-5 * cosAlt                                  # Horizontal parallax = 8.794" = 4.2635e-5 rad
        self.altitude = alt
        
        # Correct for atmospheric refraction:
        dalt = 2.967e-4 / np.tan(alt + 3.1376e-3/(alt + 8.92e-2))  # Refraction correction in altitude
        dalt *= location.pressure/101.0 * 283.0/location.temperature
        alt += dalt
        # to do: add pressure/temperature dependence
        self.altitudeRefract = alt
        
        return
    
    
    def eq2horiz(self, sinLat, cosLat, longitude,  rightAscension, declination, agst):
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
    
    
    def convertHorizontalToEquatorial(self, sinLat, cosLat, azimuth, altitude):
        
        """Convert (refraction-corrected) horizontal coordinates to equatorial coordinates.
        
        Parameters:
          sinLat   (float):  Sine of the geographic latitude of the observer.
          cosLat   (float):  Cosine of the geographic latitude of the observer.
          azimuth  (float):  Azimuth ("wind direction") of the Sun (rad; 0=South).
          altitude (float):  Altitude of the Sun above the horizon (rad).
        
        Returns: 
          tuple (float,float):  Tuple containing (hourAngle, declination):
        
            - hourAngle   (float):  Hour angle of the Sun (rad; 0=South).
            - declination (float):  Declination of the Sun (rad).
        
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
    
    
    
    def setNorthToZero(self, azimuthRefract, hourAngleRefract, computeRefrEquatorial):
        """Convert the South=0 convention to North=0 convention for azimuth and hour angle.
        
        Note:
          - South=0 is the default in celestial astronomy.
          - This function makes the angles compatible with the compass/wind directions.
        
        Parameters:
          azimuth              (float):  Azimuth ("wind direction") of the Sun (rad).
          hourAngle            (float):  Hour angle of the Sun (rad).
          computeRefrEquatorial (bool):  Compute refraction correction for equatorial coordinates.
        
        Returns: 
          tuple (float,float):  Tuple containing (azimuth, hourAngle):
        
            - azimuth   (float):  Azimuth ("wind direction") of the Sun (rad).
            - hourAngle (float):  Hour angle of the Sun (rad; 0=South).
        
        """
        
        self.azimuthRefract = (azimuthRefract + data.PI) % data.TWO_PI                    # Add PI to set North=0
        
        if(computeRefrEquatorial):
            self.hourAngleRefract = (hourAngleRefract + data.PI) % data.TWO_PI            # Add PI to set North=0
            
        return
    
    
    def convertRadiansToDegrees(self, computeRefrEquatorial):
        """Convert final results from radians to degrees.
        
        Note:
          - Does not touch intermediate results.
        
        position          (Position):  Class containing Sun position (I/O).
        computeRefrEquatorial (bool):  Compute refraction correction for equatorial coordinates.
        
        """
        
        self.longitude *= data.R2D
        self.rightAscension *= data.R2D
        self.declination *= data.R2D
        
        self.altitude *= data.R2D
        self.azimuthRefract *= data.R2D
        self.altitudeRefract *= data.R2D
        
        if(computeRefrEquatorial):
            self.hourAngleRefract *= data.R2D
            self.declinationRefract *= data.R2D
            
        return

