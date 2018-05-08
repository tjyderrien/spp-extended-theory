#!/usr/bin/env python3
#-*- coding: utf-8 -*-

# Copyright (C) 2018 F. Preucil, T.J.-Y. Derrien
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
import math, cmath
import numpy as np
from scipy.optimize import root
from itertools import product

def func(betaR, eps1, eps2, eps3, k0, t, sgn1, sgn2):
    beta = betaR[0] + betaR[1]*1.j
    kappa1 = cmath.sqrt(beta**2 - ke1)/eps1
    kappa2 = sgn1*cmath.sqrt(beta**2 - ke2)/eps2
    kappa3 = sgn2*cmath.sqrt(beta**2 - ke3)/eps3
    try:
        value = (kappa1-kappa2)*(kappa1-kappa3)*cmath.exp(-2*kappa1*eps1*t)-(kappa1+kappa2)*(kappa1+kappa3)
    except:
        value = 1E99
    return (value.real, value.imag)

def norm(vec):
    return math.sqrt(vec[0]*vec[0] + vec[1]*vec[1])

def cntr(inpt):
    px = 0
    py = 0
    ln = len(inpt)
    for pt in inpt:
        px += pt[0]
        py += pt[1]
    return (px/ln, py/ln)

# @param eps1: Dielectric permittivity of the thin film
# @param eps2: Dielectric permittivity of the half-plane below the thin film. 
# @param eps3: Dielectric permittivity of the half-plane above the thin film. Source light is supposed to come from this direction. 
# @param wavelength: wavelength of the source light (SI units). 
# @param t: thickness of the film
# @param x_min: lower boundary of Re(roots)
# @param x_max: higher boundary of Re(roots)
# @param y_min: lower boundary of Im(roots)
# @param y_max: higher boundary of Im(roots)
# @param x_steps: number of steps used to mesh the Re(roots) space. 
# @param y_steps: number of steps used to mesh the Im(roots) space. 
# @param tol_merge: tolerance to merge the identified solutions. 
def findroots(eps1, eps2, eps3, wavelength, t, x_min, x_max, y_min, y_max, x_steps, y_steps):
    global k0, ke1, ke2, ke3
    k0 = 2.*np.pi/wavelength
    ke1 = (k0**2)*eps1
    ke2 = (k0**2)*eps2
    ke3 = (k0**2)*eps3

    tol_merge = 1E3
    branches = []
    for sgn1, sgn2 in product((-1,1), (-1,1)):
        roots = []
        unique = []
        for x in np.linspace(x_min, x_max, num=x_steps):
            for y in np.linspace(y_min, y_max, num=y_steps):
                nrt = root(func, (x, y), args=(eps1, eps2, eps3, k0, t, sgn1, sgn2), method='hybr')
                if nrt.success:
                    roots.append(nrt.x)

        ln = len(roots)
        while ln > 0:
            center = roots[0]
            aux = [center]
            aux2 = []
            for rt in roots[1:]:
                if norm(rt-center) < tol_merge:
                    aux.append(rt)
                    center = cntr(aux)
                elif norm(rt+center) < tol_merge:
                    aux.append(-rt)
                    center = cntr(aux)
                else:
                    aux2.append(rt)
            if len(aux) > 4: #merging criterion
                if center[0] > 0:
                    unique.append(center)
                else:
                    unique.append([-center[0], -center[1]])
            roots = list(aux2)
            ln = len(roots)

##        valid = []
##        start = time()
##        for rt in unique:
##            beta = rt[0] + 1.j*rt[1]
##            k1 = cmath.sqrt(beta**2 - ke1)
##            k2 = sgn1*cmath.sqrt(beta**2 - ke2)
##            k3 = sgn2*cmath.sqrt(beta**2 - ke3)
##            C = cmath.exp((-k1-k3)*t/2)*(k1*eps3-k3*eps1)/(2*k1*eps3)
##            D = cmath.exp((k1-k3)*t/2)*(k1*eps3+k3*eps1)/(2*k1*eps3)
##            B1 = C*cmath.exp((k2-k1)*t/2) + D*cmath.exp((k2+k1)*t/2)
##            B2 = (C*cmath.exp((k2-k1)*t/2) - D*cmath.exp((k2+k1)*t/2))*(k1*eps2)/(k2*eps1)
##            if abs(1 - abs(B2/B1)) < tol_valid:
##                valid.append(rt)
##        prnt('Total (validated): %d' % len(valid))

##        print()
        branches.append(unique)
  # Selection of the maximum Lspp
    outs = []
    for branch in branches:
        maxy = 0
        maxx = 0
        for rt in branch:
            xi = 2.*np.pi/rt[0]
            yi = .5/rt[1]
            if yi > maxy: #select the branches with absolute maximum Lspp
                maxy = yi
                maxx = xi
        outs.append([maxx, maxy])
    return outs

#end of algorithm

#Usage:
#
#   findroots(eps1, eps2, eps3,
#             wavelength, t,
#             x_min, x_max,       }
#             y_min, y_max,       } mesh parameters
#             x_steps, y_steps)   }
#

#Example:
#roots = findroots(-1+1j, 1, 1,
                  #1, 1,
                  #-1E10, 1E10,
                  #-1E9, 1E9,
                  #30, 30)

# === PRODUCTION OF SCIENTIFIC RESULTS ===

#data
wavelength = 1026e-9 #355e-9 #1030E-9 #1026
epsTibare   = -6.206969+25.2j #800 nm
epsTiO2bare = 7.7841+0.j      #800 nm
epsCr       = -0.672122310000001+24.8657476j #-0.67+24.87j    #1026 nm
epsBK7      = 2.10277365777   #1026 nm
epsCr2O3    = 4.9713+0.1784j  #1 um [JDT Kruschwitz et al, Appl. Opt. 1997]
epsSi       = 12.8159503769+0.0114635303918j #1026 nm, Palik
epsAir      = 1.+0.j          #air
#epsCu       = -46.6046581932 + 4.7188669976j #1030 nm
epsCu       = -1.9937293241+4.9290716854j     #355  nm

# Medium 1: thin film. 
eps1 = epsCr        #thin film
# Medium 2: substrate. 
eps2 = epsBK7       #epsBK7 #environment | substrate
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
