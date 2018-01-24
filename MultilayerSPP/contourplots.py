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

import matplotlib.pyplot as plt
import numpy as np
import cmath
from scipy.constants import c, epsilon_0

#which field you want to plot
switch = 2
#0 Hy
#1 Ex
#2 Ez

#data
wavelength = 1026E-9
ne = 1E16 #np.arange(1E25, 1E28, 10) #(m^-3) quantity of electrons in conduction band
nu = (1.1E-15)**-1 #collision time between conduction band electrons
meff = 0.18

epsTibare   = -6.206969+25.2j #800 nm
epsTiO2bare = 7.7841+0.j   #800 nm
epsCr       = -0.672122310000001+24.8657476j #-0.67+24.87j    #1026 nm
epsBK7      = 2.10277365777   #1026 nm
epsCr2O3    = 4.9713+0.1784j  #1 um [JDT Kruschwitz et al, Appl. Opt. 1997]
epsSi       = 12.8159503769+0.0114635303918j #1026 nm, Palik

eps1 = epsCr        #thin film
eps2 = 1+0.j        #environment
eps3 = epsBK7 # epsBK7       #substrate

#for ne in neList:
#eps1 = Drude(wavelength, ne, epsTiO2bare, nu, meff) #thin film

k0 = 2*np.pi/wavelength
t = 100E-9 #thickness of the layer in meters

omegaeps0 = k0*c*epsilon_0

#field amplitude
A = 1.

#beta
SPPperiod = 708E-9
SPPlength = 1.37E-6
beta = 2*np.pi/SPPperiod + .5/SPPlength*1j

#branch
sgn1 = 1
sgn2 = 1
sgn3 = 1

k1 = sgn1*cmath.sqrt(beta*beta - k0*k0*eps1)
k2 = sgn2*cmath.sqrt(beta*beta - k0*k0*eps2)
k3 = sgn3*cmath.sqrt(beta*beta - k0*k0*eps3)

#plotting
xrange = 2.*wavelength
zrange = wavelength/2.
steps = 200

x = np.linspace(0, xrange, steps)

#constants precache
C = A*np.exp((-k1-k3)*t/2)*(k1*eps3-k3*eps1)/(2*k1*eps3) #might produce error (dividing by 0)
D = A*np.exp((k1-k3)*t/2)*(k1*eps3+k3*eps1)/(2*k1*eps3)
B = C*np.exp((k2-k1)*t/2) + D*np.exp((k2+k1)*t/2)

#regions
#III
z = np.linspace(t/2, zrange, steps)
xx, zz = np.meshgrid(x, z, sparse=True)
Hy = A*np.exp(1j*beta*xx-k3*zz)
Ex = Hy*1j*k3/(omegaeps0*eps3)
Ez = -Hy*beta/(omegaeps0*eps3)

if switch == 0:
    field = Hy
elif switch == 1:
    field = Ex
elif switch == 2:
    field = Ez
plt.contourf(x*1E6, z*1E6, abs(field))

#II
z = np.linspace(-t/2, t/2, steps)
xx, zz = np.meshgrid(x, z, sparse=True)
Hy = C*np.exp(1j*beta*xx+k1*zz) + D*np.exp(1j*beta*xx-k1*zz)
Ex = C*np.exp(1j*beta*xx+k1*zz)*(-1j*k1)/(omegaeps0*eps1) + D*np.exp(1j*beta*xx-k1*zz)*(1j*k1)/(omegaeps0*eps1)
Ez = Hy*beta/(omegaeps0*eps1)

if switch == 0:
    field = Hy
    name = "Hy"
elif switch == 1:
    field = Ex
    name = "Ex"
elif switch == 2:
    field = Ez
    name = "Ez"
plt.contourf(x*1E6, z*1E6, abs(field))

#I
z = np.linspace(-zrange, -t/2, steps)
xx, zz = np.meshgrid(x, z, sparse=True)
Hy = B*np.exp(1j*beta*xx+k2*zz)
Ex = -Hy*1j*k2/(omegaeps0*eps2)
Ez = -Hy*beta/(omegaeps0*eps2)

if switch == 0:
    field = Hy
elif switch == 1:
    field = Ex
elif switch == 2:
    field = Ez
plt.contourf(x*1E6, z*1E6, abs(field))

plt.xlabel(r"X ($\mu$ m)")
plt.ylabel(r"Z ($\mu$ m)")
plt.colorbar()
filename=name+"-t"+str(t*1E9)+"nm"
plt.savefig(filename+".eps")
#plt.show()
