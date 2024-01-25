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
class Location:
    """Class containing the geographic location to compute the Sun position for."""
    
    geo_longitude: float = 0.0;     """Geographic longitude of the observer/site (>0 = east of Greenwich; radians or degrees)."""
    geo_latitude:  float = 0.0;     """Geographic latitude of the observer/site (>0 = northern hemisphere; radians or degrees)."""
    
    pressure:     float = 101.0;    """Air pressure at the site (kPa)."""
    temperature:  float = 283.0;    """Air temperature at the site (K)."""
    
    
    def setLocation(self, geo_longitude, geo_latitude, pressure=101.0, temperature=283.0):
        """This function is obsolescent and will be removed in a future version.  Use set_location()
        instead."""
        _warn_obsolescent('setLocation', 'set_location', rename=True)
        return self.set_location(geo_longitude, geo_latitude, pressure, temperature)
    
    
    def set_location(self, geo_longitude, geo_latitude, pressure=101.0, temperature=283.0):
        """Setter for the details of the observer/site location to compute the Sun position for.
        
        Parameters:
          geo_longitude (float):  Geographic longitude of the observer/site (>0 = east of Greenwich; radians or degrees).
          geo_latitude  (float):  Geographic latitude of the observer/site (>0 = northern hemisphere; radians or degrees).
          
          pressure      (float):  Air pressure at the site (kPa).
          temperature   (float):  Air temperature at the site (K).
        """
        
        self.geo_longitude = geo_longitude
        self.geo_latitude  = geo_latitude
        
        self.pressure      = pressure
        self.temperature   = temperature
        
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
    
