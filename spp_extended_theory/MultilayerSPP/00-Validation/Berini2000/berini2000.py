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

# @package Multilayer_Dostovalov
# Preparation of results for Prof. Bulgakova and Sasha Dostovalov. 
# Explores the SPP theory at a thin film located 
# between two semi-infinite media. The formal model is presented in 
# T.J.-Y. Derrien et al, J. Appl. Phys. 116, 074902 (2014) and references 
# therein. 

from spp_extended_theory.MultilayerSPP.libMultilayerSPP import *

# === PRODUCTION OF SCIENTIFIC RESULTS ===

#precision
NumberOfPoints=50

#Charbonneau and Berini, Optics Letters Vol. 25, No. 11 (2000). 

#data
wavelength = 1550e-9
#epsAu       = -95.95924741872399+10.972582438155513j #1550 nm, Palik
epsAu       = -131.9475+12.65j
epsSiO2     = 2.085 #Berini #(1.4440+0j)**2 #Palik

# Medium 1: thin film. 
eps1 = epsAu        #thin film
# Medium 2: substrate. 
eps2 = epsSiO2       #epsBK7 #environment | substrate
# Medium 3: environment
eps3 = epsSiO2       #environment | substrate
# Note: Inverting eps2 and eps3 should have no effect on the possible modes, but only on field amplification. 

#t = 100E-9 #thickness of the layer in meters
#thickness = 100e-9
t_list = np.power(10., np.linspace(np.log10(1e-9), np.log10(300e-9), NumberOfPoints))
#meshes the initial guess area, all numbers are from the space of betas
x_min = -1E10
x_max = 1E10

y_min = -1E9
y_max = 1E9

x_steps = 40
y_steps = 40

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
        print((thickness, roots[branch][root_number][0], roots[branch][root_number][1]))
      #print("\n")

##TODO: from this, we would like to add a layer which will variate eps1 as function of oxide concentration
