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
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program. If not, see <http://www.gnu.org/licenses/>

# Package @libMultilayer explores the SPP theory at a thin film located 
# between two semi-infinite media. The formal model is presented in 
# T.J.-Y. Derrien et al, J. Appl. Phys. 116, 074902 (2014) and references 
# therein.
import math, cmath, pickle
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import interactive
from scipy.optimize import root
from itertools import product
from time import time

#to show info
showinfo = True

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

t = 25E-9 #thickness of the layer in meters

#branch indices
#0 (-, -, -) (+, -, -)
#1 (-, -, +) (+, -, +)
#2 (-, +, -) (+, +, -)
#3 (-, +, +) (+, +, +)

#meshes the initial guess area, all numbers are from the space of betas
x_min = -1E10
x_max = 1E10

y_min = -1E9
y_max = 1E9

x_steps = 40
y_steps = 40

#tolerances
tol_merge = 1E3   #absolute
tol_valid = 1E-1  #relative

#constants precache
k0 = 2.*np.pi/wavelength
ke1 = (k0**2)*eps1
ke2 = (k0**2)*eps2
ke3 = (k0**2)*eps3

def func(betaR):
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

def prnt(string):
    aux = '[%.1f s]' % (time() - start)
    print(string+' '*(35-len(string)-len(aux))+aux)

print('Guess rectangle:')
print('Re(beta) in [%.2e, %.2e]' % (x_min, x_max))
print('Im(beta) in [%.2e, %.2e]' % (y_min, y_max))
print('x_steps = %d' % x_steps)
print('y_steps = %d' % y_steps)
print()
print('Data:')
print('eps1:', eps1)
print('eps2:', eps2)
print('eps3:', eps3)
print('wavelength:', wavelength)
print('t:', t)
print()

#main algorithm
branches = []
for sgn1, sgn2 in product((-1,1), (-1,1)):
    roots = []
    unique = []
    print('Branch: (%s, %s)' % (('%+d' % sgn1)[0], ('%+d' % sgn2)[0]))
    start = time()
    for x in np.linspace(x_min, x_max, num=x_steps):
        for y in np.linspace(y_min, y_max, num=y_steps):
            nrt = root(func, (x, y), method='hybr')
            if nrt.success:
                roots.append(nrt.x)
    prnt('Total (converged): %d' % len(roots))

    ln = len(roots)
    start = time()
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
    prnt('Total (merged): %d' % len(unique))

    valid = []
    start = time()
    for rt in unique:
        beta = rt[0] + 1.j*rt[1]
        k1 = cmath.sqrt(beta**2 - ke1)
        k2 = sgn1*cmath.sqrt(beta**2 - ke2)
        k3 = sgn2*cmath.sqrt(beta**2 - ke3)
        C = cmath.exp((-k1-k3)*t/2)*(k1*eps3-k3*eps1)/(2*k1*eps3)
        D = cmath.exp((k1-k3)*t/2)*(k1*eps3+k3*eps1)/(2*k1*eps3)
        B1 = C*cmath.exp((k2-k1)*t/2) + D*cmath.exp((k2+k1)*t/2)
        B2 = (C*cmath.exp((k2-k1)*t/2) - D*cmath.exp((k2+k1)*t/2))*(k1*eps2)/(k2*eps1)
        if abs(1 - abs(B2/B1)) < tol_valid:
            valid.append(rt)
    prnt('Total (validated): %d' % len(valid))

    print()
    branches.append(valid)
#end of algorithm

#processing
cx = 1E6*2.*np.pi
cy = 1E9*.5

xmi, xma, ymi, yma = x_max, x_min, y_max, y_min
ixmi, ixma, iymi, iyma = cx/x_max, cx/x_min, cy/y_max, cy/y_min
ibranches = []
for branch in branches:
    ibranch = []
    for rt in branch:
        xi = cx/rt[0]
        yi = cy/rt[1]
        ibranch.append([xi, yi])
        if rt[0] < xmi:
            xmi = rt[0]
        if rt[0] > xma:
            xma = rt[0]
        if rt[1] < ymi:
            ymi = rt[1]
        if rt[1] > yma:
            yma = rt[1]
        if xi < ixmi:
            ixmi = xi
        if xi > ixma:
            ixma = xi
        if yi < iymi:
            iymi = yi
        if yi > iyma:
            iyma = yi
    ibranches.append(ibranch)
sx = (xma - xmi)/20
sy = (yma - ymi)/20
isx = (ixma - ixmi)/20
isy = (iyma - iymi)/20

#saves the roots and the parameters into a file
with open('sppdata.pkl', 'wb') as f:
    pickle.dump((branches, eps1, eps2, eps3, t, k0), f)

#plotting
colors = ('r', 'g', 'b', 'k')
plt.figure(1)
plt.gca().set_xlim((xmi-sx, xma+sx))
plt.gca().set_ylim((ymi-sy, yma+sy))

plothyp(eps1, 'r-')
plothyp(eps2, 'g-')
plothyp(eps3, 'b-')
for bind in range(len(branches)):
    for rind in range(len(branches[bind])):
        plt.scatter(*branches[bind][rind], c=colors[bind], marker='.')
        if showinfo:
            plt.gca().text(*branches[bind][rind], ' [%d, %d]' % (bind, rind))
plt.xlabel(r'Re $\beta$ [m$^{-1}$]')
plt.ylabel(r'Im $\beta$ [m$^{-1}$]')
plt.grid()
interactive(True)
plt.show()

plt.figure(2)
plt.gca().set_xlim((ixmi-isx, ixma+isx))
plt.gca().set_ylim((iymi-isy, iyma+isy))

plotinvhyp(eps1, 'r-')
plotinvhyp(eps2, 'g-')
plotinvhyp(eps3, 'b-')
for bind in range(len(ibranches)):
    for rind in range(len(ibranches[bind])):
        plt.scatter(*ibranches[bind][rind], c=colors[bind], marker='.')
        if showinfo:
            plt.gca().text(*ibranches[bind][rind], ' [%d, %d]' % (bind, rind))
plt.xlabel(r'Period [$\mu$m]')
plt.ylabel('Decay length [nm]')
plt.grid()
interactive(False)
plt.show()
