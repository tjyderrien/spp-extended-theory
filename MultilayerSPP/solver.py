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

# Module libMultilayer explores the SPP theory at a thin film located 
# between two semi-infinite media. The formal model is presented in 
# T.J.-Y. Derrien et al, J. Appl. Phys. 116, 074902 (2014) and references 
# therein. 

import math, cmath, pickle
import numpy as np

import matplotlib.pyplot as plt
from matplotlib import interactive
from scipy.optimize import root
from itertools import product

#to show info
showinfo = True

#data
wavelength = 1026E-9

epsTibare   = -6.206969+25.2j #800 nm
epsTiO2bare = 7.7841+0.j      #800 nm
epsCr       = -0.672122310000001+24.8657476j #-0.67+24.87j    #1026 nm
epsBK7      = 2.10277365777   #1026 nm
epsCr2O3    = 4.9713+0.1784j  #1 um [JDT Kruschwitz et al, Appl. Opt. 1997]
epsSi       = 12.8159503769+0.0114635303918j #1026 nm, Palik

eps1 = epsCr        #thin film
eps2 = 1.           #environment
eps3 = epsBK7       #epsBK7 #substrate

k0 = 2.*np.pi/wavelength
t = 100E-9 #thickness of the layer in meters

#branches
#0 (-, -, -) (+, -, -)
#1 (-, -, +) (+, -, +)
#2 (-, +, -) (+, +, -)
#3 (-, +, +) (+, +, +)

#Meshing the initial guess area
#Mesh enough to obtain a large enough number of solutions

x_min = -2.5E7
x_max = 2.5E7

y_min = -1E9
y_max = 1E9

x_steps = 10
y_steps =  50

#tolerances
tol_merge = 1E3 #from the space of betas

#constants

eps1s = eps1**2
eps2s = eps2**2
eps3s = eps3**2

ke1 = k0**2/eps1
ke2 = k0**2/eps2
ke3 = k0**2/eps3

def func(betaR):
    beta = betaR[0] + betaR[1]*1j
    kappa1 = cmath.sqrt(beta**2/eps1s - ke1)
    kappa2 = sgn1*cmath.sqrt(beta**2/eps2s - ke2)
    kappa3 = sgn2*cmath.sqrt(beta**2/eps3s - ke3)
    try:
        out = (kappa1-kappa2)*(kappa1-kappa3)*cmath.exp(-2*kappa1*eps1*t)-(kappa1+kappa2)*(kappa1+kappa3)
    except:
        out = 1e99
    return [out.real, out.imag]

def norm(vec):
    return math.sqrt(vec[0]*vec[0] + vec[1]*vec[1])

def cntr(inpt):
    px = 0
    py = 0
    ln = len(inpt)
    for pt in inpt:
        px += pt[0]
        py += pt[1]
    return [px/ln, py/ln]

def plothyp(eps, col):
    brpoint = cmath.sqrt(eps*k0*k0).real
    domain = np.linspace(-brpoint, -1E-9*brpoint, num=300)
    plt.plot(domain, (k0*k0*eps.imag/(2*domain)), col)
    domain = np.linspace(1E-9*brpoint, brpoint, num=300)
    plt.plot(domain, (k0*k0*eps.imag/(2*domain)), col)

def plot_segment(A, B):
    plt.plot((A[0], B[0]), (A[1], B[1]), 'gray')

def plot_hist(population, nbins):
    hist = np.histogram(population, nbins, density=True)
    bin_centers = [(hist[1][i] + hist[1][i+1])/2 for i in range(len(hist[1])-1)]
    bwdth = hist[1][2] - hist[1][1]
    plt.plot(bin_centers, hist[0])
    #plt.plot(bin_centers, np.cumsum(hist[0])*bwdth, label=str(nfibs[index]))
    plt.show()

def analyze(roots):
    dist = []
    for i in range(len(roots)):
        for j in range(i):
            dd = norm(roots[i] - roots[j])
            if dd < 1E5:
                dist.append(dd)
    plot_hist(dist, 100)

print('Guess area is the following rectangle:')
print('  Re(beta) in [%.2e, %.2e]' % (x_min, x_max))
print('  Im(beta) in [%.2e, %.2e]' % (y_min, y_max))
print('  x_steps = %d' % x_steps)
print('  y_steps = %d' % y_steps)
print()
print('Data:')
print('  eps1:', eps1)
print('  eps2:', eps2)
print('  eps3:', eps3)
print('  k0:', k0)
print('  t:', t)
print()

#main algorithm
branches = []
for sgn1, sgn2 in product((-1,1), (-1,1)):
    roots = []
    unique = []
    print('Branch: (%s, %s)' % (('%+d' % sgn1)[0], ('%+d' % sgn2)[0]))
    for x in np.linspace(x_min, x_max, num=x_steps):
        for y in np.linspace(y_min, y_max, num=y_steps):
            nrt = root(func, [x, y], method='hybr')
            if nrt.success:
                roots.append(nrt.x)
    print('Total (converged): %d   ' % len(roots))

    #analyze(roots)
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
        
    print('Total (merged): %d ' % len(unique))
    print()
    branches.append(unique)
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

sx = (xma - xmi)/10
sy = (yma - xmi)/10

isx = (ixma - ixmi)/10
isy = (iyma - ixmi)/10

with open('params.pkl', 'wb') as f:
    pickle.dump((eps1, eps2, eps3, t, k0), f)
with open('roots.pkl', 'wb') as f:
    pickle.dump(branches, f)

#plotting
plothyp(eps1, 'r-')
plothyp(eps2, 'g-')
plothyp(eps3, 'b-')

colors = ['r', 'g', 'b', 'k', 'm', 'c', 'lime', 'orangered']

plt.figure(1)
plt.gca().set_xlim([xmi-sx, xma+sx])
plt.gca().set_ylim([ymi-sy, yma+sy])
bnum = 0
for branch in branches:
    rnum = 0
    for rt in branch:
        plt.scatter(*rt, c=colors[bnum], marker='.')
        if showinfo:
            plt.gca().text(*rt, ' [%d, %d]' % (bnum, rnum))
        rnum += 1
    bnum += 1

plt.xlabel('Re ' + r'$\beta$' + ' (1/m)')
plt.ylabel('Im ' + r'$\beta$' + ' (1/m)')
plt.grid()

interactive(True)
plt.show()

plt.figure(2)
plt.gca().set_xlim([ixmi-isx, ixma+isx])
plt.gca().set_ylim([iymi-isy, iyma+isy])
bnum = 0
for branch in ibranches:
    rnum = 0
    for rt in branch:
        plt.scatter(*rt, c=colors[bnum], marker='.')
        if showinfo:
            plt.gca().text(*rt, ' [%d, %d]' % (bnum, rnum), size='10')
        rnum += 1
    bnum += 1

#plt.gca().set_xscale('log')
#plt.title(r'$\varepsilon_1=$'+str(eps1)+', '+r'$N_e=$'+str(ne)+', '+r'$\alpha=$'+str(alphainv))
plt.xlabel('Period (um)')
plt.ylabel('Decay length (nm)')
plt.grid()
interactive(False)
plt.show()
