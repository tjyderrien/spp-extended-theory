#!/usr/bin/env python
#-*- coding: utf-8 -*-

# IMPORT LIBRARIES
from SimpleSPProutines import *

# Importing data from Palik book using graphs. 
# Optical data are given in csv files. 

folder = "Database/Palik/"
filename = "Ag-Palik"

nfile = folder+filename+"-n.csv"
kfile = folder+filename+"-k.csv"

narray = loadtxt(nfile, delimiter="\t", skiprows=1)
karray = loadtxt(kfile, delimiter="\t", skiprows=1)

# import wavelength, n and k from Palik
wavelength1 = narray[:,0]
wavelength2 = karray[:,0]
n = narray[:,1]
k = karray[:,1]
numrows = 10000
base = 10
# interpolate n and k on new wavelength mesh
order=1
#wavelengths = np.arange(np.amin(wavelength2),np.amax(wavelength2), precision) #regular mesh, AWFUL for memory
print "Generating new wavelength mesh: ("+str(np.amin(wavelength2))+", "+str(np.amax(wavelength2))+")"
wavelengths = np.logspace(np.amin(np.log10(wavelength2)), np.amax(np.log10(wavelength2)), num=numrows, base=base, endpoint = True)

print "New wavelength mesh has "+str(numrows)+" rows."
#print wavelengths

fni = InterpolatedUnivariateSpline(wavelength1, n, k=order)
fki = InterpolatedUnivariateSpline(wavelength2, k, k=order)

# defining the new n and k on a common mesh
ni = fni(wavelengths)
ki = fki(wavelengths)

plt.figure()
plt.xlabel(r'$\mathcal{R}e(\varepsilon)$ ($\mu$m)')
plt.ylabel('n, k')
plt.plot(wavelength1, n, 'bs', label='n Palik')
plt.plot(wavelength2, k, 'rs', label='k Palik')
plt.plot(wavelengths, ni, 'b-', label='n interp')
plt.plot(wavelengths, ki, 'r-', label='k interp')
plt.grid()
plt.legend(loc=1)
plt.savefig('PalikData.eps')
plt.show()

#TODO: def InterpolateDielectricPermittivity()