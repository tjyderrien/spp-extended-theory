#!/usr/bin/env python2
#-*- coding: utf-8 -*-

# Copyright (C) 2013-2018 T. J.-Y. Derrien
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>

## @package Attosecond
# Module Attosecond explores ways to generate attosecond pulses using fs-laser pulses in solids. 

from libSPP import *

wavelength1 = 3100e-9

omega1 = 2*pi*c/wavelength1
omega2 = 2.18*omega1

wavelength2 = 2.*pi*c/omega2

print "Wavelength of pulse 1: "+str(wavelength1*1E9)+" nm."
print "Advised wavelength of pulse 2: "+str(wavelength2*1E9)+" nm."