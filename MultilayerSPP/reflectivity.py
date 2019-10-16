#!/usr/bin/env python2
#-*- coding: utf-8 -*-
## @package MultilayerReflectivity 
# Functions describing materials and their interaction with light. 

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

from libMaterials import *
import matplotlib.pyplot as plt

## EXAMPLE OF USAGE 
wavelength=1030e-9

## TODO: replace this by a function taking data in MaterialOpticalData.csv !
epsAir=1.;
if(wavelength == 1064e-9): 
    epsMo=-14.083065233570098+20.789041764340013j; epsSiO2=1.4496**2; epsSLG = 2.2889 + 0.000014899j
elif(wavelength == 800e-9): 
    epsMo=2.08+24.52j; epsSiO2=1.4533**2; epsSLG = 2.3018 + 0.0000075160j
elif(wavelength==1030e-9):
    epsSi   = 12.80259+0.0109j
    epsSiO2 = 2.1026565205

# Medium 1: thin film. 
eps2 = epsSiO2        #thin film
# Medium 2: substrate. 
eps3 = epsSi # | epsSiO2       #epsBK7 #environment | substrate
# Medium 3: environment
eps1 = epsAir       #environment | substrate
# Note: Inverting eps2 and eps3 should have no effect on the possible modes, but only on field amplification. 


#if(wavelength == 1064e-9): 
    #eps2=-14.083065233570098+20.789041764340013j; eps3=1.4496**2
#elif(wavelength == 800e-9): 
    #eps2=2.08+24.52j; eps3=1.4533**2
thickness_log = np.linspace(-9, np.log10(500e-9), 200)
thickness2 = np.power(10., thickness_log)

R = BiLayerReflectivity(wavelength, eps1, eps2, eps3, thickness2)

filename = "SiO2-Si-"+str(int(1E9*wavelength))+"-Reflectivity"

plt.figure()
plt.xlabel("Film thickness (nm)")
plt.ylabel("Reflectivity")
plt.plot(1E9*thickness2, R)
plt.grid()
plt.savefig(filename+".eps")
plt.show()
