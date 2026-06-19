#!/usr/bin/env python
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
from octopus_slabs.Libs.libAtomicUnits import *

wavelength1 = 1030e-9
wavelength2 = 800e-9
tau1=30E-15

omega1 = 2*pi*c/wavelength1
omega2 = 2*pi*c/wavelength2
#omega2 = 2.18*omega1

print(("Wavelength of pulse 1: "+str(wavelength1*1E9)+" nm."))
print(("Advised wavelength of pulse 2: "+str(wavelength2*1E9)+" nm."))

OneCycleMinDuration_SI = max(wavelength1,wavelength2)/c #3 cycles require this time in s

print(("Minimum required duration for a run: 3 cycles of "+str(round(1E15*OneCycleMinDuration_SI*3,0))+" fs each."))

RequiredThreeCycleTime_AU = Time_SI_to_AU(OneCycleMinDuration_SI*3)
print(("Total required duration for 3 cycles (a.u.): "+str(RequiredThreeCycleTime_AU)))
TimeStep_AU = 0.25
TimeStep_SI = Time_AU_to_SI(TimeStep_AU)
print(("TimeStep octopus = "+str(TimeStep_AU)))
print(("Required number of steps for octopus for having 3 cycles: "+str(RequiredThreeCycleTime_AU/TimeStep_AU)))
