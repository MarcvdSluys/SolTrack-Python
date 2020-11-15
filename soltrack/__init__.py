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

from .data      import Constants, Parameters
from .location  import Location
from .time      import Time
from .position  import Position
from .riseset   import RiseSet


@dataclass
class SolTrack(Location, Time, Position, RiseSet):
    
    
    def __init__(self, geoLongitude,geoLatitude, useDegrees=None, useNorthEqualsZero=None, computeRefrEquatorial=None, computeDistance=None):
        
        """Construct a SolTrack object with specified geographical location and parameters (settings).
        
        Parameters:
          geoLongitude (float):           Geographical longitude of the observer or installation (radians or degrees, depending on useDegrees).
          geoLatitude (float):            Geographical latitude of the observer or installation (radians or degrees, depending on useDegrees).
        
          useDegrees (bool):              Input (geographic position) and output are in degrees, rather than radians.
          useNorthEqualsZero (bool):      Azimuth: 0 = South, pi/2 (90deg) = West  ->  0 = North, pi/2 (90deg) = East.
          computeRefrEquatorial (bool):   Compute refraction-corrected equatorial coordinates (Hour angle, declination).
          computeDistance (bool):         Compute the distance to the Sun.
        
        Note:
          The SolTrack class is a composition of the Location, Time, Position and RiseSet classes.
        
        """
        
        # Create Constants and Parameters objects:
        self.cst       = Constants()
        self.param     = Parameters()
        self.param.setParameters(useDegrees, useNorthEqualsZero, computeRefrEquatorial, computeDistance)
        
        # Use composition to obtain the attributes from Location, Time, Position and RiseSet:
        Location.__init__(self, geoLongitude, geoLatitude)  # Set the geographic location
        Time.__init__(self)
        Position.__init__(self)
        RiseSet.__init__(self)
        
