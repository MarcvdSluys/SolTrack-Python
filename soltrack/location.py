
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
class Location:
    """Class containing the geographic location to compute the Sun position for."""
    
    geoLongitude: float = 0.0;      """Geographic longitude of the observer/site (>0 = east of Greenwich; radians or degrees)."""
    geoLatitude:  float = 0.0;      """Geographic latitude of the observer/site (>0 = northern hemisphere; radians or degrees)."""
    
    pressure:     float = 101.0;    """Air pressure at the site (kPa)."""
    temperature:  float = 283.0;    """Air temperature at the site (K)."""
    
    
    def setLocation(self, geoLongitude, geoLatitude, pressure=101.0, temperature=283.0):
        """Setter for the details of the observer/site location to compute the Sun position for.
        
        Parameters:
          geoLongitude (float):  Geographic longitude of the observer/site (>0 = east of Greenwich; radians or degrees).
          geoLatitude  (float):  Geographic latitude of the observer/site (>0 = northern hemisphere; radians or degrees).
          
          pressure     (float):  Air pressure at the site (kPa).
          temperature  (float):  Air temperature at the site (K).
        
        """
        
        self.geoLongitude = geoLongitude
        self.geoLatitude  = geoLatitude
        
        self.pressure    = pressure
        self.temperature = temperature
        
        return

