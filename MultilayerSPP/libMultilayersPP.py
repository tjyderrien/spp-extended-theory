#!/usr/bin/env python3
#-*- coding: utf-8 -*-

# Copyright (C) 2017 F. Preucil, T.J.-Y. Derrien
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

# @package libMultilayer
# Module libMultilayer explores the SPP theory at a thin film located 
# between two semi-infinite media. The formal model is presented in 
# T.J.-Y. Derrien et al, J. Appl. Phys. 116, 074902 (2014) and references 
# therein. 

import cmath
import numpy as np

from scipy.optimize import root
from numpy.linalg import norm
from itertools import product

## Function f(x,y)=0 to solve in the R^2 space. 
# f(x,y) = 0 is an equation defined in complex space. 
# All parameters can be complex-valued. 
# See: T.J.-Y. Derrien et al, J. Appl. Phys. 116, 074902 (2014) and references 
# therein. 
def func(betaR, eps1, eps2, eps3, k0, t, sgn1, sgn2, sgn3):
    beta = betaR[0] + betaR[1]*1j
    k1 = sgn1*cmath.sqrt(beta*beta - k0*k0*eps1)
    k2 = sgn2*cmath.sqrt(beta*beta - k0*k0*eps2)
    k3 = sgn3*cmath.sqrt(beta*beta - k0*k0*eps3)

    try:
        out = (k1/eps1 - k2/eps2) * (k1/eps1 - k3/eps3) * cmath.exp(-2*k1*t) - (k1/eps1 + k2/eps2) * (k1/eps1 + k3/eps3)
    except:
        out = 1e99
    return [out.real, out.imag]

## Finds many (all?) roots of the "func" function using two numerical solvers in RxR. 
# The solver is adapted to the function "func" and its mathematical dependencies. 
# WARNING: This routine cannot be used with another equation without modification. 
# In present state, the solver investigates the 8 possible branches of the thin film equation 
# for Surface Plasmon Polaritons. See T.J.-Y. Derrien et al, J. Appl. Phys. 116, 074902 (2014) and references 
# therein. 
# @param eps1: Dielectric permittivity of the thin film
# @param eps2: Dielectric permittivity of the half-plane below the thin film. 
# @param eps3: Dielectric permittivity of the half-plane above the thin film. Source light is supposed to come from this direction. 
# @param k0: wavenumber of the source light (SI units). 
# @param t: thickness of the film
# @param x_min: lower boundary of Re(roots)
# @param x_max: higher boundary of Re(roots)
# @param y_min: lower boundary of Im(roots)
# @param y_max: higher boundary of Im(roots)
# @param x_steps: number of steps used to mesh the Re(roots) space. 
# @param y_steps: number of steps used to mesh the Im(roots) space. 
# @param t_blur: tolerance to merge the identified solutions. 
def findroots(eps1, eps2, eps3, k0, t, x_min, x_max, y_min, y_max, x_steps, y_steps, t_blur):
    broots = []
    for sgn1, sgn2, sgn3 in product((-1,1), (-1,1), (-1,1)):
        roots = []
        unique = []
        for x in np.linspace(x_min, x_max, num=x_steps):
            for y in np.linspace(y_min, y_max, num=y_steps):                                                #take an initial guess [x, y] from within the specified grid, x_steps and y_steps determine the grid density
                nrt = root(func, [x, y], args=(eps1, eps2, eps3, k0, t, sgn1, sgn2, sgn3), method='hybr')   #uses the MINPACK method: http://www.netlib.org/minpack/
                if nrt.success:                                                                             #if it converges, add it to the list of roots
                    roots.append(nrt.x)

        ln = len(roots)
        while ln > 0:                                                                                       #this part should remove duplicate roots (considering tolerance t_blur)
            center = roots[0]                                                                               #take the first element of the list and put it in 'center'
            aux = [center]
            aux2 = []
            for rt in roots[1:]:                                                                            #go thru the rest of roots
                if norm(rt - center) < t_blur:                                                              #if the current root 'rt' is sufficiently near center, add it to the new list of similar roots 'aux'
                    aux.append(rt)
                    center = np.mean(aux, axis = 0)                                                         #update 'center' to include the new root, 'center' always lies in the middle
                else:
                    aux2.append(rt)
            if len(aux) > 2:
                unique.append(center)                                                                       #'unique' contains separated roots
            roots = list(aux2)                                                                              #continue the algorithm with the rest
            ln = len(roots)

        broots.append(unique)
    return broots
findroots = np.vectorize(findroots)

#Usage:
#
#   findroots(eps1, eps2, eps3,
#             k0, t,
#             x_min, x_max,       }
#             y_min, y_max,       } mesh parameters
#             x_steps, y_steps,   }
#             t_blur)
#
#returns a list of branches, each branch is a list of roots, each root is a 1D numpy array containing two values, first value = Re(beta), second value = Im(beta)

#Example:
#roots = findroots(-1+1j, 1, 1,
                  1, 1,
                  -10, 10,
                  -10, 10,
                  5, 5,
                  .0001)

#print(roots) #list of arrays
#print(roots[0])     #prints roots belonging to the zeroth branch (-, -, -)
#print(roots[7][0])  #prints the first root from the last branch (+, +, +)

#List of branch indexes:
#   0 (-, -, -)
#   1 (-, -, +)
#   2 (-, +, -)
#   3 (-, +, +)
#   4 (+, -, -)
#   5 (+, -, +)
#   6 (+, +, -)
#   7 (+, +, +)

