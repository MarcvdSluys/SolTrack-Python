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


"""SolTrack module

SolTrack is a simple, free, fast and accurate Python package to compute the position of the Sun, as well as
its rise and set times.  SolTrack can be used under the conditions of the EUPL 1.2 licence.  These pages
contain the API documentation.  For more information on the Python package, licence, source code and data
files, see the [SolTrack homepage](http://soltrack.sf.net) and [Van der Sluys & Van Kan
(2022)](https://arxiv.org/abs/2209.01557) (open access scientific paper with all technical details).
"""

name = 'soltrack'

from dataclasses import dataclass

from .data      import Parameters
from .location  import Location
from .time      import Time
from .position  import Position
from .riseset   import RiseSet


@dataclass
class SolTrack(Location, Time, Position, RiseSet):
    
    
    def __init__(self, geo_longitude,geo_latitude, use_degrees=None, use_north_equals_zero=None, compute_refr_equatorial=None, compute_distance=None):
        
        """Construct a SolTrack object with specified geographical location and parameters (settings).
        
        Parameters:
          geo_longitude (float):           Geographical longitude of the observer or installation (radians or degrees, depending on use_degrees).
          geo_latitude (float):            Geographical latitude of the observer or installation (radians or degrees, depending on use_degrees).
        
          use_degrees (bool):              Input (geographic position) and output are in degrees, rather than radians.
          use_north_equals_zero (bool):    Azimuth: 0 = South, pi/2 (90deg) = West  ->  0 = North, pi/2 (90deg) = East.
          compute_refr_equatorial (bool):  Compute refraction-corrected equatorial coordinates (Hour angle, declination).
          compute_distance (bool):         Compute the distance to the Sun.
        
        Note:
          The SolTrack class is a composition of the Location, Time, Position and RiseSet classes.
        """
        
        # Create Parameters objects:
        self.param     = Parameters()
        self.param.set_parameters(use_degrees, use_north_equals_zero, compute_refr_equatorial, compute_distance)
        
        # Use composition to obtain the attributes from Location, Time, Position and RiseSet:
        Location.__init__(self, geo_longitude, geo_latitude)  # Set the geographic location
        Time.__init__(self)
        Position.__init__(self)
        RiseSet.__init__(self)
        
