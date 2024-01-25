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


@dataclass
class Parameters:
    """Class containing SolTrack parameters/settings."""
    
    _use_degrees:              bool  = False;   """Input (geographic position) and output are in degrees."""
    _use_north_equals_zero:    bool  = False;   """Azimuth: 0 = South, pi/2 (90deg) = West  ->  0 = North, pi/2 (90deg) = East."""
    _compute_refr_equatorial:  bool  = True;    """Compute refraction-corrected equatorial coordinates (Hour angle, declination)."""
    _compute_distance:         bool  = True;    """Compute the distance to the Sun."""
    
    
    def setParameters(self, use_degrees=None, use_north_equals_zero=None, compute_refr_equatorial=None,
                      compute_distance=None):
        """This function is obsolescent and will be removed in a future version.  Use set_parameters()
        instead."""
        _warn_obsolescent('setParameters', 'set_parameters', rename=True)
        return self.set_parameters(use_degrees, use_north_equals_zero, compute_refr_equatorial, compute_distance)
    
        
    def set_parameters(self, use_degrees=None, use_north_equals_zero=None, compute_refr_equatorial=None,
                       compute_distance=None):
        """Set the SolTrack parameters (settings).
        
        Parameters:
          use_degrees (bool):              Input (geographic position) and output are in degrees, rather than radians.
          use_north_equals_zero (bool):    Azimuth: 0 = South, pi/2 (90deg) = West  ->  0 = North, pi/2 (90deg) = East.
          compute_refr_equatorial (bool):  Compute refraction-corrected equatorial coordinates (Hour angle, declination).
          compute_distance (bool):         Compute the distance to the Sun.
        """
        
        if use_degrees is not None:              self._use_degrees              = use_degrees
        if use_north_equals_zero is not None:    self._use_north_equals_zero    = use_north_equals_zero
        if compute_refr_equatorial is not None:  self._compute_refr_equatorial  = compute_refr_equatorial
        if compute_distance is not None:         self._compute_distance         = compute_distance
        
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
