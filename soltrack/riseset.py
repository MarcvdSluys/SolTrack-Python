
"""SolTrack: a simple, free, fast and accurate C routine to compute the position of the Sun.
  
  Copyright (c) 2014-2020  Marc van der Sluys, Paul van Kan and Jurgen Reintjes,
  Sustainable Energy research group, HAN University of applied sciences, Arnhem, The Netherlands
   
  This file is part of the SolTrack package, see: http://soltrack.sourceforge.net
  SolTrack is derived from libTheSky (http://libthesky.sourceforge.net) under the terms of the GPL v.3
  
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

from .data import Constants


@dataclass
class RiseSet(Constants):
    """Class containing attributes and methods to compute the rise,transit and set times of the Sun and their azimuths/altitudes."""
    
    def __init__(self, param):
        self.param = param
    
    
