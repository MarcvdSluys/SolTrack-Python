
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


@dataclass
class Constants:
    """Class containing SolTrack constants."""
    
    PI: float         = 3.14159265358979323846;   """Pi"""
    TWO_PI: float     = 6.28318530717958647693;   """2 pi"""
    R2D: float        = 57.2957795130823208768;   """Radians to degrees conversion factor"""
    R2H: float        = 3.81971863420548805845;   """Radians to hours conversion factor"""



@dataclass
class Parameters:
    """Class containing SolTrack parameters/settings."""
    
    useDegrees: bool             = False;   """Input (geographic position) and output are in degrees"""
    useNorthEqualsZero: bool     = False;   """Azimuth: 0 = South, pi/2 (90deg) = West  ->  0 = North, pi/2 (90deg) = East"""
    computeRefrEquatorial: bool  = False;   """Compute refraction-corrected equatorial coordinates (Hour angle, declination)"""
    computeDistance: bool        = False;   """Compute the distance to the Sun"""
    
    
    def setParameters(self, useDegrees=False, useNorthEqualsZero=False, computeRefrEquatorial=False,
                      computeDistance=False):
        """Set the SolTrack parameters (settings).
        
        Parameters:
          useDegrees (bool):              Input (geographic position) and output are in degrees, rather than radians.
          useNorthEqualsZero (bool):      Azimuth: 0 = South, pi/2 (90deg) = West  ->  0 = North, pi/2 (90deg) = East.
          computeRefrEquatorial (bool):   Compute refraction-corrected equatorial coordinates (Hour angle, declination).
          computeDistance (bool):         Compute the distance to the Sun.
        """
        
        self.useDegrees             = useDegrees
        self.useNorthEqualsZero     = useNorthEqualsZero
        self.computeRefrEquatorial  = computeRefrEquatorial
        self.computeDistance        = computeDistance
        
        return
