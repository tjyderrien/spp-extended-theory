#!/usr/bin/env python3
#-*- coding: utf-8 -*-

# Copyright (C) 2019 F. Preucil, T.J.-Y. Derrien
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
from time import time
import numpy.linalg as LA

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

## Plot the hyperbola
def plothyp(eps, col):
    radius = k0*k0*eps.imag/2.
    domain = np.linspace(xmi, min(xma, cmath.sqrt(eps*k0*k0).real), num=1000)
    plt.plot(domain, radius/domain, col, linewidth=.75)

def plotinvhyp(eps, col):
    radius = k0*k0*eps.imag/2.
    if radius != 0:
        domain = np.linspace(xmi, min(xma, cmath.sqrt(eps*k0*k0).real), num=1000)
        plt.plot(cx/domain, domain*cy/radius, col, linewidth=.75)

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

# Outputs: branch number correspond to the signs of (k2, k3). There is invariance per sign inversion of k1.  
# (0: --, 1: -+, 2:+-, 3:++). 

def findroots(eps1, eps2, eps3, wavelength, t, x_min, x_max, y_min, y_max, x_steps, y_steps, num_of_maxs=1):
    global k0, ke1, ke2, ke3
    k0 = 2.*np.pi/wavelength
    ke1 = (k0**2)*eps1
    ke2 = (k0**2)*eps2
    ke3 = (k0**2)*eps3

    tol_merge = 1E2 #1E3
    merge_treshold = 1 #4? 1: helps to not miss some modes
    branches = []
    for sgn1, sgn2 in product((-1,1), (-1,1)): 
    #for sgn1, sgn2 in product((1,1), (1,1)): #DEBUG LINE
    #for sgn1, sgn2 in ((-1, 1), (-1, 1), (1,-1), (1,1)): #DEBUG LINE
        roots = []
        unique = []
        start = time()
        for x in np.linspace(x_min, x_max, num=x_steps):
            for y in np.linspace(y_min, y_max, num=y_steps):
                nrt = root(func, (x, y), args=(eps1, eps2, eps3, k0, t, sgn1, sgn2), method='hybr')
                #checking = func(nrt.x, eps1, eps2, eps3, k0, t, sgn1, sgn2)
                if (nrt.success):
                    roots.append(nrt.x)
                #if (nrt.success and (checking.all() < 1E-10 ) ): #checking if solution verifies the equation
                    #nrt2 = root(func, nrt.x, args=(eps1, eps2, eps3, k0, t, sgn1, sgn2), method='lm')
                    #if nrt2.success:
                        #roots.append(nrt2.x)

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
        #outs.append(sorted(brinv, key=lambda x:x[0], reverse=True)[:num_of_maxs])
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

## We turn the nonlinear solver into a matrix form, and find out the eigen vealues directly, without going through a nonlinear solver. 
# @param k1: thin film wavenumber, 
# @param k2: environment wavenumber, 
# @param k3: substrate wavenumber
# @param eps1: thin film complex permittivity
# @param eps2: environment complex permittivity
# @param eps3: substrate complex permittivity
def ThreeLayerEigenSolver(t, k1, k2, k3, eps1, eps2, eps3):
    A11=np.exp(-k3*t/2.); A12=0.; A13=-np.exp(k1*t/2.); A14=-np.exp(-k1*t/2.);
    A21=k3/eps3*np.exp(-k3*t/2.); A22=0.; A23=k1/eps1*np.exp(k1*t/2.); A24=-k1/eps3*np.exp(-k1*t/2.)
    A31=0.; A32=-np.exp(-k2*t/2.); A33=np.exp(-k1*t/2.); A34=np.exp(k1*t/2.) 
    A41=0.; A42=k2/eps2*np.exp(-k2*t/2.); A43=-k1/eps1*np.exp(-k1*t/2.); A44=k1/eps1*np.exp(k1*t/2.)
    matrix=np.array([[A11, A12, A13, A14], [A21, A22, A23, A24], [A31, A32, A33, A34], [A41, A42, A43, A44]])
    del A11, A12, A13, A14, A21, A22, A23, A24, A31, A32, A33, A34, A41, A42, A43, A44
    w,v = LA.eig(matrix)
    return w

# Conversion from period to to beta/k0: 
def PeriodToBetaNorm(period, wavelength): 
    # period = 2.*np.pi / beta.real
    k0 = 2.*np.pi / wavelength
    beta_norm_re = np.divide(np.divide(2.*np.pi, period), k0)
    return beta_norm_re

# Conversion from SPP decay length to beta/k0:
def LsppToBetaNorm(Lspp): 
    # Lspp = 0.5/Im(beta) #Im(beta)=0.5/Lspp. 
    k0 = 2.*np.pi / wavelength
    beta_norm_im = np.divide(np.divide(0.5, Lspp), k0)
    return beta_norm_im

#wavelength=633e-9
#numberofroots = 10
## Medium 1: thin film.
#eps1 = -19.+0.53j     #Ag thin film
## Medium 2: substrate. 
#eps2 = 4. #3.999999+0.004j
## Medium 3: environment
#eps3 = 1.5**2 # eps2 #eps2: symmetric modes       #environment | substrate
## Note: Inverting eps2 and eps3 should have no effect on the possible modes, but only on field amplification. 
#k0 = 2.*np.pi/wavelength
#ke1 = (k0**2)*eps1
#ke2 = (k0**2)*eps2
#ke3 = (k0**2)*eps3

#sgn1 = 1.
#sgn2 = 1.

#beta = 1.+1j

#k1 = cmath.sqrt(beta**2 - ke1)/eps1
#k2 = sgn1*cmath.sqrt(beta**2 - ke2)/eps2
#k3 = sgn2*cmath.sqrt(beta**2 - ke3)/eps3
#sol = ThreeLayerEigenSolver(28e-9, k1, k2, k3, eps1, eps2, eps3)
#print(sol)

# The best would be to converge this towards finding beta such as matrix would be solved. 
# But still, for any beta, we obtain 4 eigen values. 
