#!/usr/bin/env python3
#-*- coding: utf-8 -*-

# Copyright (C) 2018 T.J.-Y. Derrien
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

# @package Multilayer_MirzaLevy
# Preparation of results for Prof. Bulgakova, Inam Mirza and Yoann Levy. 
# Explores the SPP theory at a thin film located 
# between two semi-infinite media. The formal model is presented in 
# T.J.-Y. Derrien et al, J. Appl. Phys. 116, 074902 (2014) and references 
# therein. 

from libMultilayerSPP import *
import matplotlib.pyplot as plt
from libMaterials import *
# === PRODUCTION OF SCIENTIFIC RESULTS ===

#precision
NumberOfPoints=30

#data
wavelength = 30e-9 #355e-9 #1030E-9 #1026
#neH20 = 1E28 #the most violent change
epsPd = 0.69885+0.45747j
epsFe = 0.94860+0.22374j
epsC  = 0.80621+0.12173j

#epsTibare   = -6.206969+25.2j #800 nm
#epsTiO2bare = 7.7841+0.j      #800 nm
#epsAir      = 1.+0.j          #air
#neTiO2      = np.power(10., np.linspace(np.log10(1E24), np.log10(1E29), NumberOfPoints))
#epsTiO2     = np.linspace(-2,1,NumberOfPoints)
#print("Shape of Ne array: ", np.shape(neTiO2)[0])

#t = 100E-9 #thickness of the layer in meters
#thickness = 10e-9
t_list = np.power(10., np.linspace(np.log10(1e-9), np.log10(30e-9), NumberOfPoints))
#meshes the initial guess area, all numbers are from the space of betas
x_min = -1E10
x_max = 1E10

y_min = -1E9
y_max = 1E9

x_steps = 40
y_steps = 40

# Medium 2: substrate. 
eps2 = epsPd      #epsBK7 #environment | substrate
# Medium 3: environment
eps3 = epsFe       #environment | substrate
# Note: Inverting eps2 and eps3 should have no effect on the possible modes, but only on field amplification. 
#for order in np.arange(0,np.shape(neTiO2)[0]-1):
#for order in np.arange(0,np.shape(epsTiO2)[0]-1):
# Medium 1: thin film. 
#print("Order: ", order)
#eps1 = Drude(wavelength, neTiO2[order], epsTiO2bare, 1.1e-15**-1, 0.18)
eps1=epsC #epsTiO2[order]
#thin film

for thickness in t_list:
    roots = findroots(eps1, eps2, eps3,
            wavelength, thickness,
            x_min, x_max,    
            y_min, y_max,    
            x_steps, y_steps)

    # Shaping the data to plot them with GNUplot
    roots_shape = np.shape(roots)
    #print(roots_shape)
    num_thickness = np.shape(thickness) #NOTE: is this used? 
    num_branches  = roots_shape[0]
    num_roots     = roots_shape[1]
    num_property  = roots_shape[2]

    for branch in np.arange(0,num_branches):
        for root_number in np.arange(0,num_roots): 
            #print(neTiO2[order], thickness, roots[branch][root_number][0], roots[branch][root_number][1], eps1.real, eps1.imag)
            print(eps1.real, thickness, roots[branch][root_number][0], roots[branch][root_number][1], eps1.real, eps1.imag, eps1.real*eps2.real+eps1.imag*eps2.imag, eps1.real*eps3.real+eps1.imag*eps3.imag)
            #print(neTiO2[order], thickness, roots[branch][root_number][0],
        #print("\n")

##TODO: from this, we would like to add a layer which will variate eps1 as function of oxide concentration
