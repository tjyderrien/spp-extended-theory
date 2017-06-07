#!/usr/bin/env python
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

## @package libMultilayer
# Module libMultilayer explores the SPP theory at a thin film located 
# between two semi-infinite media. The formal model is presented in 
# T.J.-Y. Derrien et al, J. Appl. Phys. 116, 074902 (2014) and references 
# therein. 

import cmath
import numpy as np

import matplotlib.pyplot as plt
from scipy.optimize import root
from scipy.linalg import norm
from itertools import product
from libMaterials import Drude
from libDatabase import ExportToTxt

pi = 2.*np.arcsin(1.0)

#data


def func(betaR): #{{{
    beta = betaR[0] + betaR[1]*1j
    k1 = sgn1*cmath.sqrt(beta*beta - k0*k0*eps1)
    k2 = sgn2*cmath.sqrt(beta*beta - k0*k0*eps2)
    k3 = sgn3*cmath.sqrt(beta*beta - k0*k0*eps3)

    try:
        out = (k1/eps1 - k2/eps2) * (k1/eps1 - k3/eps3) * cmath.exp(-2.*k1*t) - (k1/eps1 + k2/eps2) * (k1/eps1 + k3/eps3)
    except:
        out = 1e20
    return [out.real, out.imag]
#}}}

wavelength = 800e-9
ne=5e27 #np.arange(1E25, 1E28, 10) #(m^-3) quantity of electrons in conduction band
nu = (1.1E-15)**-1 #collision time between conduction band electrons
meff = 0.18

eps2 = 1.+0.j       #environment
eps3 = 13.64+0.048j #substrate Si (no excitation)

#for ne in neList:
eps1 = Drude(wavelength, ne, eps3, nu, meff)

print eps1

k0 = 2.*pi/wavelength
t = 100e-9 #Thickness of the layer in meters

#Meshing the initial guess area
#Mesh enough to obtain a large enough number of solutions

#delta = .01 #not used

maxvalue = 1E8
maxsteps = 100

x_min   = -maxvalue #real part of beta
x_max   =  maxvalue
x_steps =  maxsteps #round((x_max - x_min)/delta)

y_min   = -maxvalue
y_max   =  maxvalue
y_steps =  maxsteps #round((y_max - y_min)/delta)

#tolerances

t_root = 1E-30 #tolerance for the solution
t_sim  = 1E-30 #limits the radius for one converged solution (looks to be a relative quantity)

print('Guess area is a rectangle: [%s, %s]x[%si, %si], x_steps = %d, y_steps = %d' % (x_min, x_max, y_min, y_max, x_steps, y_steps))
print()
broots = []
for sgn1, sgn2, sgn3 in product((-1,1), (-1,1), (-1,1)):
    roots = []
    out = []
    print('Branch: [%s][%s][%s]' % (('%+d' % sgn1)[0], ('%+d' % sgn2)[0], ('%+d' % sgn3)[0]))
    for x in np.linspace(x_min, x_max, num=x_steps):
        for y in np.linspace(y_min, y_max, num=y_steps):

            nrt = root(func, [x, y], method='hybr', tol=t_root) #tolerance for the solution (relative ? absolute ?)
            if nrt.success:
                roots.append(nrt.x)

    while len(roots) > 0:
        cntr = roots[0]
        same = [cntr]
        nrest = []

        for rt in roots[1:]:
            if norm(rt - cntr, 1) < t_sim:
                same.append(rt)
                cntr = np.mean(same, axis = 0)
            else:
                nrest.append(rt)

        out.append(cntr)
        roots = list(nrest)

    for rt in out:
        print('%+.20g %+.20g    %+.20g %+.20g' % (rt[0], rt[1], func(rt)[0], func(rt)[1]))
    print('Total: %d' % len(out))
    print()
    broots.append(out)
    
#print("Number of roots"+str(broots))

##a = cmath.sqrt(eps1*k0*k0)
##plt.plot((a.real, -a.real), (a.imag, -a.imag), 'r')

colors = ['red', 'green', 'blue', 'black', 'magenta', 'cyan', 'lime', 'orangered']
for rts, col in zip(broots, colors):
  try:
    plt.scatter(*zip(*rts), c=col)
  except:
    print "Jumped"

ExportToTxt(broots, 'betaSolution-Ne'+str(ne)+'-t'+str(t)+'.log')

plt.xlabel('Re')
plt.ylabel('Im')
plt.grid()
plt.legend(loc=4)
plt.savefig('betaSolution-Ne'+str(ne)+'-t'+str(t)+'.eps')
plt.show()
