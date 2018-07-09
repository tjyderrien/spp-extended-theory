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
# @param num_of_maxs: how many maxima to return
def findroots(eps1, eps2, eps3, wavelength, t, x_min, x_max, y_min, y_max, x_steps, y_steps, num_of_maxs):
    global k0, ke1, ke2, ke3
    k0 = 2.*np.pi/wavelength
    ke1 = (k0**2)*eps1
    ke2 = (k0**2)*eps2
    ke3 = (k0**2)*eps3

    tol_merge = 1E3
    merge_treshold = 4
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
            if len(aux) > merge_treshold: #merging criterion
                if center[0] > 0:
                    unique.append(center)
                else:
                    unique.append([-center[0], -center[1]])
            roots = list(aux2)
            ln = len(roots)
        branches.append(unique)

    #selection of maximal Lspps
    outs = []
    for branch in branches:
        brinv = [[2.*np.pi/rt[0], .5/rt[1]] for rt in branch]
        outs.append(sorted(brinv, key=lambda x:abs(x[1]), reverse=True)[:num_of_maxs])
    return outs

#end of algorithm

#Usage:
#
#   findroots(eps1, eps2, eps3,
#             wavelength, t,
#             x_min, x_max,       }
#             y_min, y_max,       } mesh parameters
#             x_steps, y_steps,   }
#             num_of_maxs)

#Example:
#roots = findroots(-1+1j, 1, 3,
                  #1, .4,
                  #-1E2, 1E2,
                  #-1E4, 1E4,
                  #30, 30,
                  #2)
