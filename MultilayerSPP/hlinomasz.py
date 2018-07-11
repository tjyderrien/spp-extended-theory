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

# @package Multilayer_Hlinomasz
# Preparation of results for Prof. Bulgakova and Krystof Hlinomasz. 
# Explores the SPP theory at a thin film located 
# between two semi-infinite media. The formal model is presented in 
# T.J.-Y. Derrien et al, J. Appl. Phys. 116, 074902 (2014) and references 
# therein. 

from function import *
import matplotlib.pyplot as plt

# === PRODUCTION OF SCIENTIFIC RESULTS ===

#data
wavelength = 1064e-9
epsAir      = 1.+0.j          #air

if(wavelength == 1064e-9): 
    epsMo=-14.083065233570098+20.789041764340013j; epsSiO2=1.4496**2
elif(wavelength == 800e-9): 
    epsMo=2.08+24.52j; epsSiO2=1.4533**2

# Medium 1: thin film. 
eps1 = epsMo        #thin film
# Medium 2: substrate. 
eps2 = epsSiO2       #epsBK7 #environment | substrate
# Medium 3: environment
eps3 = epsAir       #environment | substrate
# Note: Inverting eps2 and eps3 should have no effect on the possible modes, but only on field amplification. 

#t = 100E-9 #thickness of the layer in meters
t_list = np.arange(10e-9, 100e-9, 10e-9)
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
    print("thickness, [[period, Lspp]]: ", thickness, roots, "\n")
    #print("")
