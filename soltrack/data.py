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


@dataclass
class Constants:
    """Class containing SolTrack constants."""
    
    _PI:     float     = 3.14159265358979323846;   """Pi"""
    _TWOPI:  float     = 6.28318530717958647693;   """2 pi"""
    _R2D:    float     = 57.2957795130823208768;   """Radians to degrees conversion factor"""
    _R2H:    float     = 3.81971863420548805845;   """Radians to hours conversion factor"""
    _D2R:    float     = 0.01745329251994329576;   """Degrees to radians conversion factor"""



@dataclass
class Parameters:
    """Class containing SolTrack parameters/settings."""
    
    _useDegrees:             bool  = False;   """Input (geographic position) and output are in degrees."""
    _useNorthEqualsZero:     bool  = False;   """Azimuth: 0 = South, pi/2 (90deg) = West  ->  0 = North, pi/2 (90deg) = East."""
    _computeRefrEquatorial:  bool  = True;    """Compute refraction-corrected equatorial coordinates (Hour angle, declination)."""
    _computeDistance:        bool  = True;    """Compute the distance to the Sun."""
    
    
    def setParameters(self, useDegrees=None, useNorthEqualsZero=None, computeRefrEquatorial=None,
                      computeDistance=None):
        """Set the SolTrack parameters (settings).
        
        Parameters:
          useDegrees (bool):              Input (geographic position) and output are in degrees, rather than radians.
          useNorthEqualsZero (bool):      Azimuth: 0 = South, pi/2 (90deg) = West  ->  0 = North, pi/2 (90deg) = East.
          computeRefrEquatorial (bool):   Compute refraction-corrected equatorial coordinates (Hour angle, declination).
          computeDistance (bool):         Compute the distance to the Sun.
        """
        
        if(useDegrees is not None):             self._useDegrees             = useDegrees
        if(useNorthEqualsZero is not None):     self._useNorthEqualsZero     = useNorthEqualsZero
        if(computeRefrEquatorial is not None):  self._computeRefrEquatorial  = computeRefrEquatorial
        if(computeDistance is not None):        self._computeDistance        = computeDistance
        
        return
