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
whichfield = 3
#0 Hy
#1 Ex
#2 Ez
#3 Sx
#4 Sz
#5 |E|
#6 |S|

#which part you want to plot (only for whichfield < 5)
whichpart = 0
#0 Re
#1 Im
#2 Abs

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

k0 = 2.*np.pi/wavelength
#t = 100e-9 #thickness of the layer in meters

omegaeps0 = k0*c*epsilon_0

#field amplitude
A = 1.

#beta

SPPperiod = 1030E-9; SPPlength = 5.36E-6; t = 42e-9   #SPP          42nm
#SPPperiod = 723E-9; SPPlength = 2.47E-6 ; t = 42e-9   #Hybride      42nm
#SPPperiod = 295E-9; SPPlength = 23E-9   ; t = 42e-9   #LambdaOverN  42nm

#SPPperiod = 1025E-9; SPPlength = 4.37E-6 ; t = 100E-9 #SPP         100nm
#SPPperiod = 708E-9;  SPPlength = 1.37E-6 ; t = 100E-9 #Hybride     100nm
#SPPperiod = 295E-9;  SPPlength = 23E-9   ; t = 100E-9 #LambdaOverN 100nm

beta = 2.*np.pi/SPPperiod + 1.j*.5/SPPlength

#branch
sgn1 = 1
sgn2 = 1
sgn3 = 1

k1 = sgn1*cmath.sqrt(beta**2 - k0**2*eps1)
k2 = sgn2*cmath.sqrt(beta**2 - k0**2*eps2)
k3 = sgn3*cmath.sqrt(beta**2 - k0**2*eps3)

#plotting
xrange = 10.*wavelength
zrange = wavelength/4.
steps = 500

x = np.linspace(0, xrange, steps)

#constants precache
C = A*np.exp((-k1-k3)*t/2)*(k1*eps3-k3*eps1)/(2*k1*eps3) #might produce error (dividing by 0)
D = A*np.exp((k1-k3)*t/2)*(k1*eps3+k3*eps1)/(2*k1*eps3)
B = C*np.exp((k2-k1)*t/2) + D*np.exp((k2+k1)*t/2)

def plotfield():
    global name, part
    toplot = 0*xx*zz
    if whichfield == 0:
        name = 'Hy'
        field = Hy
    elif whichfield == 1:
        name = 'Ex'
        field = Ex
    elif whichfield == 2:
        name = 'Ez'
        field = Ez
    elif whichfield == 3:
        name = 'Sx'
        field = Sx
    elif whichfield == 4:
        name = 'Sz'
        field = Sz
    elif whichfield == 5:
        name = ']E['
        toplot = np.sqrt(Ex*np.conj(Ex) + Ez*np.conj(Ez))
    elif whichfield == 6:
        name = ']S['
        toplot = np.sqrt(Sx*np.conj(Sx) + Sz*np.conj(Sz))

    if whichfield in (0, 1, 2, 3, 4):
        if whichpart == 0:
            part = 'RealPart'
            toplot = field.real
        elif whichpart == 1:
            part = 'ImaginaryPart'
            toplot = field.imag
        elif whichpart == 2:
            part = 'AbsoluteValue'
            toplot = abs(field)
    elif whichfield in (5, 6):
        part = 'TotalValue'

    plt.contourf(x*1E6, z*1E6, toplot)

#regions
#III
z = np.linspace(t/2, zrange, steps)
xx, zz = np.meshgrid(x, z, sparse=True)
Hy = A*np.exp(1j*beta*xx-k3*zz)
Ex = Hy*1j*k3/(omegaeps0*eps3)
Ez = -Hy*beta/(omegaeps0*eps3)
Sx = -Ez*Hy
Sz = Ex*Hy
plotfield()

#II
z = np.linspace(-t/2, t/2, steps)
xx, zz = np.meshgrid(x, z, sparse=True)
Hy = C*np.exp(1j*beta*xx+k1*zz) + D*np.exp(1j*beta*xx-k1*zz)
Ex = C*np.exp(1j*beta*xx+k1*zz)*(-1j*k1)/(omegaeps0*eps1) + D*np.exp(1j*beta*xx-k1*zz)*(1j*k1)/(omegaeps0*eps1)
Ez = Hy*beta/(omegaeps0*eps1)
Sx = -Ez*Hy
Sz = Ex*Hy
plotfield()

#I
z = np.linspace(-zrange, -t/2, steps)
xx, zz = np.meshgrid(x, z, sparse=True)
Hy = B*np.exp(1j*beta*xx+k2*zz)
Ex = -Hy*1j*k2/(omegaeps0*eps2)
Ez = -Hy*beta/(omegaeps0*eps2)
Sx = -Ez*Hy
Sz = Ex*Hy
plotfield()

plt.xlabel(r'x ($\mu$m)')
plt.ylabel(r'z ($\mu$m)')
plt.colorbar()
filename='Period'+str(SPPperiod*1E9)+'nm-Lspp'+str(SPPlength*1E6)+'um-'+name+'-'+part+'-t'+str(t*1E9)+'nm'
plt.savefig(filename+'.eps')
plt.show()
