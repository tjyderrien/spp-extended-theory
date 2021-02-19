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

## Package @MultilayerWithsources (experimental) explores the SPP theory at a thin film located 
# between two semi-infinite media. The formal model is presented in 
# T.J.-Y. Derrien et al, J. Appl. Phys. 116, 074902 (2014) and has been 
# enriched by F. Preucil (HiLASE Centre) to take into account more accurate
# source by adding a propagating wave in the ambient medium of irradiation.
# This module remains experimental and may not be valid. 

import matplotlib.pyplot as plt
import numpy as np
import cmath
from itertools import product
from scipy.constants import c, epsilon_0

branch_index = 1 #3

#two black lines to show the boundaries
showlines = True

#contourlevels
levels = 10

#branch indices
#0 (-, -, -) (+, -, -)
#1 (-, -, +) (+, -, +)
#2 (-, +, -) (+, +, -)
#3 (-, +, +) (+, +, +)

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
eps3 = 1.+0.j       #environment | substrate
# Note: Inverting eps2 and eps3 should have no effect on the possible modes, but only on field amplification. 

t = 30E-9 #thickness of the layer in meters
k0 = 2.*np.pi/wavelength

betaR = (500000, 600000) #What is this? 

#retrieves beta
beta = betaR[0] + 1.j*betaR[1]
SPPperiod = 2.*np.pi/betaR[0]
SPPlength = .5/betaR[1]

#which field you want to plot
whichfield = 0
#0 Hy
#1 Ex
#2 Ez
#3 Sx
#4 Sz
#5 |E|
#6 |Sinst|
#7 |Savg|

#which part you want to plot (only for whichfield < 5)
whichpart = 0
#0 Re
#1 Im
#2 Abs

#data
##wavelength = 1026E-9
##ne = 1E16 #np.arange(1E25, 1E28, 10) #(m^-3) quantity of electrons in conduction band
##nu = (1.1E-15)**-1 #collision time between conduction band electrons
##meff = 0.18
##
##epsTibare   = -6.206969+25.2j #800 nm
##epsTiO2bare = 7.7841+0.j   #800 nm
##epsCr       = -0.672122310000001+24.8657476j #-0.67+24.87j    #1026 nm
##epsBK7      = 2.10277365777   #1026 nm
##epsCr2O3    = 4.9713+0.1784j  #1 um [JDT Kruschwitz et al, Appl. Opt. 1997]
##epsSi       = 12.8159503769+0.0114635303918j #1026 nm, Palik
##
##eps1 = epsCr        #thin film
##eps2 = 1+0.j        #environment
##eps3 = epsBK7       #epsBK7       #substrate
##
##k0 = 2.*np.pi/wavelength
##t = 40E-9 #thickness of the layer in meters

##beta =  7.593326478142351e+06  + 1.j*6.229229632723185e+07
##SPPperiod = 1030E-9; SPPlength = 5.36E-6; t = 42e-9   #SPP          42nm
##SPPperiod = 723E-9; SPPlength = 2.47E-6 ; t = 42e-9   #Hybride      42nm
##SPPperiod = 295E-9; SPPlength = 23E-9   ; t = 42e-9   #LambdaOverN  42nm
##
##SPPperiod = 1025E-9; SPPlength = 4.37E-6 ; t = 100E-9 #SPP         100nm
##SPPperiod = 708E-9;  SPPlength = 1.37E-6 ; t = 100E-9 #Hybride     100nm
##SPPperiod = 295E-9;  SPPlength = 23E-9   ; t = 100E-9 #LambdaOverN 100nm
##
##beta = 2.*np.pi/SPPperiod + 1.j*.5/SPPlength

#constants precache
omegaeps0 = k0*c*epsilon_0

#setting the branch
sgn1, sgn2 = list(product((-1,1), (-1,1)))[branch_index]
k1 = cmath.sqrt(beta**2 - k0**2*eps1)
k2 = sgn1*cmath.sqrt(beta**2 - k0**2*eps2)
k3 = sgn2*cmath.sqrt(beta**2 - k0**2*eps3)

#plotting
xrange = wavelength*50
zrange = wavelength*.5
steps = 200

x = np.linspace(0., xrange, steps)
z = np.linspace(-zrange, zrange, steps)
x_mesh, z_mesh = np.meshgrid(x, z, sparse=True)
#x_mesh, z_mesh = x[None,:], z[:,None] #this might be faster than meshgrid but it seems it's not

#amplitudes precache
S = 1
matrix = [[cmath.exp(-t/2*k3), 0, -cmath.exp(t/2*k1), -cmath.exp(-t/2*k1)],
[1j*cmath.exp(-t/2*k3)*k3/(eps3*omegaeps0), 0, 1j*cmath.exp(t/2*k1)*k1/(eps1*omegaeps0), -1j*cmath.exp(-t/2*k1)*k1/(eps1*omegaeps0)],
[0, -cmath.exp(-t/2*k2), cmath.exp(-t/2*k1), cmath.exp(t/2*k1)],
[0, 1j*cmath.exp(-t/2*k2)*k2/(eps2*omegaeps0), -1j*cmath.exp(-t/2*k1)*k1/(eps1*omegaeps0), 1j*cmath.exp(t/2*k1)*k1/(eps1*omegaeps0)]]

GF = np.linalg.inv(matrix)

source = [-S*cmath.exp(t/2*k3), 1j*S*cmath.exp(t/2*k3)*k3/(eps3*omegaeps0), 0, 0]
result = np.matmul(GF, source)

A, B, C, D = result

print('Field amplitudes:')
print()
print(('A:', A))
print(('C:', C))
print(('D:', D))
print(('B:', B))

#field functions

def Hy(x, z):
    if z >= t/2:
        return A*np.exp(1.j*beta*x-k3*z) + S*np.exp(1.j*beta*x+k3*z)
    elif t/2 > z > -t/2:
        return C*np.exp(1.j*beta*x+k1*z) + D*np.exp(1.j*beta*x-k1*z)
    else:
        return B*np.exp(1.j*beta*x+k2*z)
def Ex(x, z):
    if z >= t/2:
        return A*np.exp(1.j*beta*x-k3*z)*1.j*k3/(omegaeps0*eps3) - S*np.exp(1.j*beta*x+k3*z)*1.j*k3/(omegaeps0*eps3)
    elif t/2 > z > -t/2:
        return (-C*np.exp(1.j*beta*x+k1*z) + D*np.exp(1.j*beta*x-k1*z))*1.j*k1/(omegaeps0*eps1)
    else:
        return -B*np.exp(1.j*beta*x+k2*z)*1.j*k2/(omegaeps0*eps2)
def Ez(x, z):
    if z >= t/2:
        return -A*np.exp(1.j*beta*x-k3*z)*beta/(omegaeps0*eps3) - S*np.exp(1.j*beta*x+k3*z)*beta/(omegaeps0*eps3)
    elif t/2 > z > -t/2:
        return -(C*np.exp(1.j*beta*x+k1*z) + D*np.exp(1.j*beta*x-k1*z))*beta/(omegaeps0*eps1)
    else:
        return -B*np.exp(1.j*beta*x+k2*z)*beta/(omegaeps0*eps2)
Hy_vect = np.vectorize(Hy)
Ex_vect = np.vectorize(Ex)
Ez_vect = np.vectorize(Ez)

if whichfield == 0:
    name = 'Hy'
    field = Hy_vect(x_mesh, z_mesh)
elif whichfield == 1:
    name = 'Ex'
    field = Ex_vect(x_mesh, z_mesh)
elif whichfield == 2:
    name = 'Ez'
    field = Ez_vect(x_mesh, z_mesh)
elif whichfield == 3:
    name = 'Sx_inst'
    Hy_mesh = Hy_vect(x_mesh, z_mesh)
    Ez_mesh = Ez_vect(x_mesh, z_mesh)
    field = -Ez_mesh.real*Hy_mesh.real
elif whichfield == 4:
    name = 'Sz_inst'
    Hy_mesh = Hy_vect(x_mesh, z_mesh)
    Ex_mesh = Ex_vect(x_mesh, z_mesh)
    field = Ex_mesh.real*Hy_mesh.real
elif whichfield == 5:
    name = ']E['
    Ex_mesh = Ex_vect(x_mesh, z_mesh)
    Ez_mesh = Ez_vect(x_mesh, z_mesh)
    toplot = np.sqrt(Ex_mesh.real**2 + Ez_mesh.real**2)
elif whichfield == 6:
    name = ']S_inst['
    Hy_mesh = Hy_vect(x_mesh, z_mesh)
    Ex_mesh = Ex_vect(x_mesh, z_mesh)
    Ez_mesh = Ez_vect(x_mesh, z_mesh)
    toplot = np.sqrt((Ez_mesh.real*Hy_mesh.real)**2 + (Ex_mesh.real*Hy_mesh.real)**2)
elif whichfield == 7:
    name = ']S_avg['
    Hy_mesh = Hy_vect(x_mesh, z_mesh)
    Ex_mesh = Ex_vect(x_mesh, z_mesh)
    Ez_mesh = Ez_vect(x_mesh, z_mesh)
    Sx_avg_mesh = -Ez_mesh*np.conj(Hy_mesh)
    Sz_avg_mesh = Ex_mesh*np.conj(Hy_mesh)
    toplot = np.sqrt(Sx_avg_mesh*np.conj(Sx_avg_mesh) + Sz_avg_mesh*np.conj(Sz_avg_mesh))
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
elif whichfield in (5, 6, 7):
    part = 'TotalValue'
plt.contourf(x*1E6, z*1E6, toplot, levels, cmap=plt.get_cmap('Blues'))

if showlines:
    plt.plot([x[0]*1E6, x[-1]*1E6], [t/2*1E6, t/2*1E6], 'k-', linewidth = .2)
    plt.plot([x[0]*1E6, x[-1]*1E6], [-t/2*1E6, -t/2*1E6], 'k-', linewidth = .2)

plt.xlabel(r'x [$\mu$m]')
plt.ylabel(r'z [$\mu$m]')
plt.colorbar()
filename = 'Period'+str(SPPperiod*1E9)+'nm-Lspp'+str(SPPlength*1E6)+'um-'+name+'-'+part+'-t'+str(t*1E9)+'nm'
plt.savefig(filename+'.eps')
plt.show()
