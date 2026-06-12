#!/usr/bin/env python
#-*- coding: utf-8 -*-
## @package Reflectivity3layers
# Three-layer reflectivity, transmission and absorptivity

# Copyright (C) 2013-2024 T. J.-Y. Derrien
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

from spp_extended_theory.Libs.libMaterials import *
from spp_extended_theory.Libs.libImportOpticalData import *
import matplotlib.pyplot as plt

## EXAMPLE OF USAGE 
wavelengths = np.linspace(100e-9, 2e-6, 100)

## TODO: replace this by a function taking data in MaterialOpticalData.csv !
epsAir=1.;

epsSiO2 = ImportPalikDatabase_epsilon_fromNK(wavelengths, "SiO2-Palik")
epsAu   = ImportPalikDatabase_epsilon_fromNK(wavelengths, "Au-Palik")


"""if(wavelength == 1064e-9): 
    epsMo=-14.083065233570098+20.789041764340013j; epsSiO2=1.4496**2; epsSLG = 2.2889 + 0.000014899j
elif(wavelength == 800e-9): 
    epsMo=2.08+24.52j; epsSiO2=1.4533**2; epsSLG = 2.3018 + 0.0000075160j; epsAu = -26.154188586+1.8503881331j; epsBK7=(1.5108+9.2656e-9j)**2
elif(wavelength==1030e-9):
    epsSi   = 12.80259+0.0109j
    epsSiamorphous=13.0522658011
    epsSiO2 = 2.1026565205
    epsTi   = -4.2656+27.277j
    epsMo   = -11.6291789477+20.6107572133j
elif(wavelength==400e-9):
    epsMo=-1.1887745095+19.5149912779j; 
    epsSiO2=2.1614446988; #epsSLG = 2.2889 + 0.000014899j
"""

# Medium 1: External medium. 
eps1 = epsAir
# Medium 2: Thin film. 
eps2 = epsSiO2       #epsBK7 #environment | substrate
# Medium 3: Substrate
eps3 = MaxwellGarnett2(epsAu, epsAir, fraction)
# Note: Inverting eps2 and eps3 should have no effect on the possible modes, but only on field amplification. 
eps4 = epsSiO2

## TODO: need to interpolate eps1, 2, 3, 4 onto the same mesh!

R = TriLayerReflectivity(wavelength, eps1, eps2, eps3, eps4, thickness2, thickness3)
# T = ThreeLayerTransmission(wavelength, eps1, eps2, eps3, thickness2)
# A = 1.-R-T
filename = "Air-SiO2-Au-SiO2-Reflectivity"

plt.figure()
plt.xlabel("Wavelength (nm)")
plt.ylabel("Reflectivity")
plt.semilogx(1E9*wavelengths, R, label=r"$R(\lambda)$")
# plt.semilogx(1E9*thickness2, T, label="T")
# plt.semilogx(1E9*thickness2, A, label="A")
plt.legend(loc="best")
plt.grid()
plt.savefig(filename+".eps")
plt.savefig(filename+".png")
plt.tight_layout()
plt.show()
