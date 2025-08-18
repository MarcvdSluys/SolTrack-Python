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
from astroconst import r2d as _R2D, pi as _PI, pi2 as _TWOPI

from .data import Parameters


@dataclass
class Position(Parameters):
    """Class containing the position of the Sun and related attributes and methods."""

    
    def __init__(self):
        Parameters.__init__(self)
        
        # Time:
        self.lt:                       float = None;     """The local date/time for the desired instant, if any"""
        self.utc:                      float = 0.0;      """The universal date/time (UTC) for the desired instant"""
        self.julian_day:               float = 0.0;      """The Julian day for the desired instant"""
        
        self._tJD:                     float = 0.0;      """Time in Julian days since 2000.0"""
        self._tJC:                     float = 0.0;      """Time in Julian centuries since 2000.0"""
        self._tJC2:                    float = 0.0;      """Time in Julian centuries since 2000.0 squared"""
        
        # Ecliptical coordinates:
        self.longitude:                float = 0.0;      """Ecliptical longitude of the Sun (radians)"""
        self.distance:                 float = 0.0;      """Distance Earth-Sun (AU)"""
        
        # Obliquity of the ecliptic and nutation:
        self._obliquity:               float = 0.0;      """Obliquity of the ecliptic (radians)"""
        self._cos_obliquity:           float = 0.0;      """Cosine of the obliquity of the ecliptic"""
        self._nutation_lon:            float = 0.0;      """Nutation in longitude (radians)"""
        
        # Equatorial coordinates and sidereal time:
        self._right_ascension_uncorr:  float = 0.0;      """Right ascension of the Sun, UNCORRECTED for refraction (radians)"""
        self._declination_uncorr:      float = 0.0;      """Declination of the Sun, UNCORRECTED for refraction (radians)"""
        self.declination:              float = None;     """Declination of the Sun, corrected for refraction (radians)"""
        
        self._agst:                    float = 0.0;      """Apparent Greenwich sidereal time for the instant of interest (radians)"""
        self.hour_angle:               float = None;     """Hour angle of the Sun, corrected for refraction (radians)"""
        
        # Horizontal coordinates:
        self._altitude_uncorr:         float = 0.0;      """Altitude of the Sun, UNCORRECTED for refraction (radians)"""
        self.altitude:                 float = 0.0;      """Altitude of the Sun, corrected for refraction (radians)"""
        self.azimuth:                  float = 0.0;      """Azimuth of the Sun, corrected for refraction (radians)"""
        
        
    
    
    def computePosition(self):
        """This function is obsolescent and will be removed in a future version.  Use compute_position()
        instead."""
        _warn_obsolescent('computePosition', 'compute_position', rename=True)
        return self.compute_position()
    
    
    def compute_position(self):
        """ Method to compute the position of the Sun.
        """
        
        import astrotool.date_time as at_dt
        
        # If the user uses degrees, convert the geographic location to radians:
        if self.param._use_degrees:
            self.geo_longitude /= _R2D
            self.geo_latitude  /= _R2D
        
        # Compute these once and reuse:
        self._sin_lat = np.sin(self.geo_latitude)
        self._cos_lat = np.sqrt(1.0 - self._sin_lat**2)  # Cosine of a latitude is always positive or zero
        
        
        # Compute the Julian Day from the date and time:
        self.julian_day = at_dt.jd_from_date_time(self.year, self.month, self.day, self.hour, self.minute, self.second)
        
        # CHECK1: the line below instead of above alone makes the whole compute_position() call ~36% slower!
        # However, combining this with the removal of self.year-self.second in set_date_time() is only slightly
        # slower, but causes problems in computeRiseSet().  See CHECK1 in those places.
        # self.julian_day = self.utc.to_julian_date().to_numpy()
        
        # Derived expressions for time, to be reused:
        self._tJD  = self.julian_day - 2451545.0    # Time in Julian days since 2000.0
        self._tJC  = self._tJD/36525.0              # Time in Julian centuries since 2000.0
        self._tJC2 = self._tJC**2                   # T^2
        
        
        # Compute the ecliptic longitude of the Sun and the obliquity of the ecliptic and nutation:
        self._compute_longitude(self.param._compute_distance)
        
        # Convert ecliptic coordinates to geocentric equatorial coordinates:
        self._convert_ecliptic_to_equatorial(self.longitude, self._cos_obliquity)
        
        # Convert equatorial coordinates to horizontal coordinates, correcting for parallax and refraction:
        self._convert_equatorial_to_horizontal()
        
        
        # Convert the corrected horizontal coordinates back to equatorial coordinates:
        if self.param._compute_refr_equatorial:
            self._convert_horizontal_to_equatorial(self._sin_lat, self._cos_lat, self.azimuth,
                                                   self.altitude)
            
        # Use the North=0 convention for azimuth and hour angle (default: South = 0) if desired:
        if self.param._use_north_equals_zero:
            self._set_north_to_zero(self.azimuth, self.hour_angle)
            
        # If the user wants degrees, convert final results from radians to degrees:
        if self.param._use_degrees:
            self.geo_longitude *= _R2D          # Convert back to original
            self.geo_latitude  *= _R2D          # Convert back to original
            self._convert_radians_to_degrees()  # Convert final results
        
        return
    
    
    
    def _compute_longitude(self, compute_distance=True):
        """Compute the ecliptic longitude of the Sun for a given instant.
        
        Note:
          - Also computes and stores the obliquity of the ecliptic and nutation.
        
        Parameters:
          compute_distance (bool):  Compute distance to the Sun.  Note that this results in a marginally better
                                    longitude as well.  (optional, default=True).

        """
        
        l0 = 4.895063168 + 628.331966786 * self._tJC  +  5.291838e-6 * self._tJC2        # Mean longitude
        ma = 6.240060141 + 628.301955152 * self._tJC  -  2.682571e-6 * self._tJC2        # Mean anomaly
        
        sec = (3.34161088e-2 - 8.40725e-5* self._tJC - 2.443e-7*self._tJC2)*np.sin(ma) + \
            (3.489437e-4 - 1.76278e-6*self._tJC)*np.sin(2*ma)                            # Sun's equation of the centre
        odot = l0 + sec                                                                  # True longitude
        
        
        # Nutation, aberration:
        omg  = 2.1824390725 - 33.7570464271 * self._tJC  + 3.622256e-5 * self._tJC2      # Lon. of Moon's mean ascending node
        dpsi = -8.338601e-5*np.sin(omg)                                                  # Nutation in longitude
        dist = np.ones_like(l0)*1.0000010178                                             # Mean distance to the Sun in AU
        if compute_distance:
            ecc = 0.016708634 - 0.000042037   * self._tJC  -  0.0000001267 * self._tJC2  # Eccentricity of the Earth's orbit
            nu = ma + sec                                                                # True anomaly
            dist = dist*(1.0 - ecc**2)/(1.0 + ecc*np.cos(nu))                            # Geocentric distance of the Sun in AU
            
        aber = -9.93087e-5/dist                                                          # Aberration
        
        # Obliquity of the ecliptic and nutation - do this here, since we've already computed many of the ingredients:
        eps0 = 0.409092804222 - 2.26965525e-4*self._tJC - 2.86e-9*self._tJC2             # Mean obliquity of the ecliptic
        deps = 4.4615e-5*np.cos(omg)                                                     # Nutation in obliquity
        
        # Save position parameters:
        self.longitude = (odot + aber + dpsi) % _TWOPI                                   # Apparent geocentric longitude, referred to the true equinox of date
        
        self.distance = dist                                                             # Distance (AU)
        
        self._obliquity   = eps0 + deps                                                  # True obliquity of the ecliptic
        self._cos_obliquity = np.cos(self._obliquity)                                    # Need the cosine later on
        self._nutation_lon = dpsi                                                        # Nutation in longitude
        
        return
    
    
    
    def _convert_ecliptic_to_equatorial(self, longitude, cosObliquity):
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
        
        self._right_ascension_uncorr  = np.arctan2(cosObliquity*sinLon, np.cos(longitude)) % _TWOPI  # 0 <= azimuth < 2pi
        self._declination_uncorr     = np.arcsin(sinObl*sinLon)
        
        return
    
    
    def _convert_equatorial_to_horizontal(self):
        """Convert equatorial to horizontal coordinates.
        
        Also corrects for parallax and atmospheric refraction.
        
        """
        
        # We need the AGST for the coordinate transformation:
        gmst       = 4.89496121 + 6.300388098985*self._tJD + 6.77e-6*self._tJC2  # Greenwich mean sidereal time
        self._agst = gmst + self._nutation_lon * self._cos_obliquity             # Correction for equation of the equinoxes . apparent Greenwich sidereal time
        
        
        # Do the actual coordinate transformation:
        self.azimuth, sinAlt = self._eq2horiz_ct(self._sin_lat,self._cos_lat, self.geo_longitude,
                                                 self._right_ascension_uncorr, self._declination_uncorr,
                                                 self._agst)
        
        alt = np.arcsin( sinAlt )                                  # Altitude of the Sun above the horizon (rad)
        cosAlt = np.sqrt(1.0 - sinAlt**2)                          # Cosine of the altitude is always positive or zero
        
        # Correct for parallax:
        alt -= 4.2635e-5 * cosAlt                                  # Horizontal parallax = 8.794" = 4.2635e-5 rad
        self._altitude_uncorr = np.copy(alt)                       # Sun altitude, uncorrected for refraction.  If not copied, _altitude_uncorr and altitude will be converted to degrees TWICE!
        
        # Correct for atmospheric refraction:
        dalt = 2.967e-4 / np.tan(alt + 3.1376e-3/(alt + 8.92e-2))  # Refraction correction in altitude
        dalt *= self.pressure/101.0 * 283.0/self.temperature
        alt += dalt
        self.altitude = alt                                        # Sun altitude, corrected for atmospheric refraction
        
        return
    
    
    def _eq2horiz_ct(self, sinLat, cosLat, longitude,  rightAscension, declination, agst):
        """The actual coordinate transformation to convert equatorial to horizontal coordinates.
        
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
        cosDec = np.sqrt(1.0 - sinDec**2)                            # Cosine of a declination is always >=0
        tanDec = sinDec/cosDec
        
        azimuth = np.arctan2( sinHa,  cosHa  * sinLat - tanDec * cosLat )  # 0 <= azimuth < 2pi
        sinAlt = sinLat * sinDec + cosLat * cosDec * cosHa                 # Sine of the altitude above the horizon
        
        return azimuth, sinAlt
    
    
    def _convert_horizontal_to_equatorial(self, sinLat, cosLat, azimuth, altitude):
        """Convert (refraction-corrected) horizontal coordinates to the equatorial coordinates hour_angle and
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
        
        self.hour_angle   = np.arctan2( sinAz,   cosAz  * sinLat + tanAlt * cosLat )     # Local Hour Angle:  0 <= hour_angle < 2pi
        self.declination = np.arcsin(  sinLat * sinAlt  -  cosLat * cosAlt * cosAz  )    # Declination
        
        return
    
    
    
    def _set_north_to_zero(self, azimuth, hour_angle):
        """Convert the South=0 convention to North=0 convention for azimuth and hour angle.
        
        Note:
          - South=0 is the default in celestial astronomy.
          - This function makes the angles compatible with the compass/wind directions.
        
        Parameters:
          azimuth              (float):  Azimuth ("wind direction") of the Sun (rad).
          hour_angle            (float):  Hour angle of the Sun (rad).
        
        """
        
        self.azimuth = (azimuth + _PI) % _TWOPI                    # Add PI to set North=0
        
        if self.param._compute_refr_equatorial:
            self.hour_angle = (hour_angle + _PI) % _TWOPI          # Add PI to set North=0
            
        return
    
    
    def _convert_radians_to_degrees(self):
        """Convert final results from radians to degrees.
        
        Note:
          - Does not touch intermediate results.
        
        """
        
        self.longitude *= _R2D
        self._right_ascension_uncorr *= _R2D
        self._declination_uncorr *= _R2D
        
        self._altitude_uncorr *= _R2D
        self.azimuth *= _R2D
        self.altitude *= _R2D
        
        if self.param._compute_refr_equatorial:
            self.hour_angle *= _R2D
            self.declination *= _R2D
            
        return
    
    
    def _rev_pi(self, angle):
        """Fold an angle in radians to take a value between -PI and +PI.
        
        Parameters:
          angle (float):  Angle to fold (radians).
        
        Returns:
          float: Angle between -PI and PI (radians).
        
        """
        
        return ((angle + _PI) % _TWOPI) - _PI
    
    
    def create_df(self, utc=False, jd=False, ecl=False, eq=False, uncorr=False, rts_pos=False):
        """Create a Pandas DataFrame with the results of the Sun position and rise/set data.

        Parameters:
          utc     (bool):  Include utc, defaults to False.
          jd      (bool):  Include Julian day, defaults to False.
          ecl     (bool):  Include ecliptical coordinates, defaults to False.
          eq      (bool):  Include equatorial coordinates, defaults to False.
          uncorr  (bool):  Include coordinates uncorrected for refraction, defaults to False.
          rts_pos (bool):  Include the rise, transit and set positions, defaults to False.
        
        Note that if a desired variable is not available, the request will be silently ignored.
        """
        
        # Add date and time:
        if self.lt is not None:
            self.df = pd.DataFrame(data=self.lt, columns=['LT'])
            if utc: self.df['UTC'] = self.utc
        else:
            self.df = pd.DataFrame(data=self.utc, columns=['UTC'])  # New df with initial column 0-90 (91 ints)
        if jd: self.df['julianDay'] = self.julian_day
        
        # Add ecliptical coordinates:
        if ecl:
            self.df['longitude'] = self.longitude
            self.df['latitude']  = 0
        self.df['distance'] = self.distance
        
        # Add uncorrected coordinates:
        if uncorr:
            if eq:
                self.df['raUncorr']   = self._right_ascension_uncorr
                self.df['declUncorr'] = self._declination_uncorr
            self.df['altUncorr'] = self._altitude_uncorr
        
        # Add horizontal coordinates:
        self.df['azimuth']  = self.azimuth
        self.df['altitude'] = self.altitude
        
        # Add equatorial coordinates:
        if eq and self.hour_angle is not None:
            self.df['hourAngle'] = self.hour_angle
            self.df['declination'] = self.declination
        
        # Add rise, transit and set data if available:
        if self.transit_time is not None:
            self.df['riseTime']  = self.rise_time
            self.df['transTime'] = self.transit_time
            self.df['setTime']   = self.set_time
            
            # Add rise, transit and set positions if desired:
            if rts_pos:
                self.df['riseAzimuth']   = self.rise_azimuth
                self.df['transAltitude'] = self.transit_altitude
                self.df['setAzimuth']    = self.set_azimuth
        
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
