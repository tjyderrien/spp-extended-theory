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

# Module libMultilayer explores the SPP theory at a thin film located 
# between two semi-infinite media. The formal model is presented in 
# T.J.-Y. Derrien et al, J. Appl. Phys. 116, 074902 (2014) and references 
# therein. 

import cmath
import numpy as np

import matplotlib.pyplot as plt
from scipy.optimize import root
from numpy.linalg import norm
from itertools import product
from libMaterials import Drude
#from libDatabase import ExportToTxt

prec = 8 #printing precision

#data
wavelength = 800E-9
ne = 5E27 #np.arange(1E25, 1E28, 10) #(m^-3) quantity of electrons in conduction band
nu = (1.1E-15)**-1 #collision time between conduction band electrons
meff = 0.18

eps2 = 1+0.j       #environment
eps3 = -6.206969+25.2j #substrate Si (no excitation)

#for ne in neList:
eps1 = Drude(wavelength, ne, 7.7841+0.j, nu, meff)

k0 = 2.*np.pi/wavelength
t = 150E-9 #thickness of the layer in meters

#branches
branches = [0, 1, 2, 3] #list of branches you want to use

#0 (-, -, -)
#1 (-, -, +)
#2 (-, +, -)
#3 (-, +, +)
#4 (+, -, -)
#5 (+, -, +)
#6 (+, +, -)
#7 (+, +, +)

#Meshing the initial guess area
#Mesh enough to obtain a large enough number of solutions

maxvalue = 1E8
maxsteps = 100

x_min   = 1E4
x_max   = 1E9

y_min   = -10E8
y_max   = 10E8

x_steps =  100
y_steps =  150

#tolerances
tol_merge = 1E5 #from the space of betas

#-----------------------------------------------------------------------------

def func(betaR):
    beta = betaR[0] + betaR[1]*1j
    k1 = sgn1*cmath.sqrt(beta*beta - k0*k0*eps1)
    k2 = sgn2*cmath.sqrt(beta*beta - k0*k0*eps2)
    k3 = sgn3*cmath.sqrt(beta*beta - k0*k0*eps3)
    try:
        out = (k1/eps1 - k2/eps2) * (k1/eps1 - k3/eps3) * cmath.exp(-2*k1*t) - (k1/eps1 + k2/eps2) * (k1/eps1 + k3/eps3)
    except:
        out = 1e99
    return [out.real, out.imag]

def plothyp(eps, col):
    brpoint = cmath.sqrt(eps*k0*k0).real
    domain = np.linspace(-brpoint, -1E-6*brpoint, num=2000)
    plt.plot(domain, k0*k0*eps.imag/(2*domain), col)
    domain = np.linspace(1E-6*brpoint, brpoint, num=2000)
    plt.plot(domain, k0*k0*eps.imag/(2*domain), col)

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
print('Selected branches:', branches)

#main algorithm
broots = []
for sgn1, sgn2, sgn3 in [list(product((-1,1), (-1,1), (-1,1)))[i] for i in branches]:
    roots = []
    unique = []
    num = 0
    print()
    print()
    print('Branch: (%s, %s, %s)' % (('%+d' % sgn1)[0], ('%+d' % sgn2)[0], ('%+d' % sgn3)[0]))
    for x in np.linspace(x_min, x_max, num=x_steps):
        print('Tracing the grid: %d%%' % (num/x_steps*100), end = '\r')
        for y in np.linspace(y_min, y_max, num=y_steps):                #take an initial guess [x, y] from within the specified grid, x_steps and y_steps determine the grid density
            nrt = root(func, [x, y], method='hybr')                     #uses the MINPACK method: http://www.netlib.org/minpack/
            if nrt.success:                                             #if it converges, add it to the list of roots
                roots.append(nrt.x)
        num += 1

    lns = len(roots)
    print('Total (converged): %d   ' % lns)

    ln = lns
    while ln > 0:                                                       #this part should remove duplicate roots (considering tolerance tol_merge)
        print('Merging roots: %d%%' % (100 - ln/lns*100), end = '\r')
        center = roots[0]                                               #take the first element of the list and put it in 'center'
        aux = [center]
        aux2 = []
        for rt in roots[1:]:                                            #go thru the rest of roots
            if norm(rt - center) < tol_merge:                           #if the current root 'rt' is sufficiently near center, add it to the new list of similar roots 'aux'
                aux.append(rt)
                center = np.mean(aux, axis = 0)                         #update 'center' to include the new root, 'center' always lies in the middle
            else:
                aux2.append(rt)
        if len(aux) > 2:
            unique.append(center)                                       #'unique' contains separated roots
        roots = list(aux2)                                              #continue the algorithm with the rest
        ln = len(roots)

    if len(unique) > 0:
        print(' '*18)
        print('Real%s     Imaginary%s│  Abs(F)  %s│  Levenberg–Marquardt' % tuple([' '*prec]*3))
        print('%s┼%s┼%s' % ('─'*(2*prec + 18), '─'*(prec + 10), '─'*(2*prec + 24)))
    for rt in unique:
        nrt = root(func, rt, method='lm')
        print(('%% .%de  %% .%de  │ %% .%de  │  %%s  %% .%de  %% .%de' % tuple([prec]*5)) % (rt[0], rt[1], norm(func(rt)), nrt.success, nrt.x[0], nrt.x[1])) #checking with Levenberg–Marquardt method
    if len(unique) > 0:
        print()
    print('Total (merged): %d ' % len(unique))
    broots.append(unique)

#plotting
bounds = [0, 0, 0, 0]
for rts in broots:
    for rt in rts:
        if rt[0] < bounds[0]:
            bounds[0] = rt[0]
        if rt[0] > bounds[1]:
            bounds[1] = rt[0]
        if rt[1] < bounds[2]:
            bounds[2] = rt[1]
        if rt[1] > bounds[3]:
            bounds[3] = rt[1]

#plothyp(eps1, 'r-')
#plothyp(eps2, 'g-')
#plothyp(eps3, 'b-')

broots2 = []
for branch in broots:
  aux = []
  for rt in branch:
    #print(rt, 2*np.pi/rt[0], .5/rt[1])
    aux.append([1E6*np.pi/rt[0], 1E9*.5/rt[1]])
  broots2.append(aux)
  
#print(broots2)
colors = ['r', 'g', 'b', 'k', 'm', 'c', 'lime', 'orangered']
for rts, col in zip(broots, colors):
    if len(rts) > 0:
        plt.scatter(*zip(*rts), c=col)

#ExportToTxt(broots, 'betaSolution-Ne'+str(ne)+'-t'+str(t)+'.log')

axes = plt.gca()
#axes.set_xlim([-1E7,1E7])
#axes.set_ylim([-1E7,1E7])
#axes.set_xscale('log')
alphainv = 1.E9/((4.*np.pi/wavelength * np.sqrt(eps1)).imag)
plt.title(r'$\varepsilon_1=$'+str(eps1)+', '+r'$N_e=$'+str(ne)+', '+r'$\alpha=$'+str(alphainv))
plt.xlabel('Period (um)')
plt.ylabel('Decay length (nm)')
#plt.axis([1.1*bnd for bnd in bounds])
plt.grid()
#plt.legend(loc=4)
#plt.savefig('betaSolution-wavelength'+str(wavelength*1E9)+'-thickness'+str(t)+'-density-'+str(ne)+'m-3.eps')
plt.show()
