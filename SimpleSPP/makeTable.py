#!/usr/bin/env python
#-*- coding: utf-8 -*-

# IMPORT LIBRARIES
from SimpleSPProutines import *

# PHYSICAL INPUT
wavelength = 800e-9
epsAir=1e0
epsSi0=13.64+0.048j
meffe=0.18
nuSi=(1.1e-15)**-1

#example=period(betaSPP(wavelength, epsAir, Drude(wavelength, 1e28, epsSi0, nuSi)))

#print example

"""
We would like now to construct a database using available materials description with all possible interfaces
we will : 

1. For each material in database, select each material and verify, for each available wavelength, 
1.1: If SPP condition is verified, 
1.2. yes, then period can be calculated and shown;
1.4. SPP decay depth in medium 1
1.5. SPP decay depth in medium 2
2. Then extract a table which contains all possible scenarios
"""
# Select database
database="MaterialOpticalDatabaseForPlasmonics.csv"

# Build database array for choosing which material can be of interest to irradiate
dbarray = loadtxt(database, dtype='str', delimiter='\t')
#SPPactiveInterfaces(dbarray, 'new')

""" TODO: interface this with HTML for publication on the web. 
1. Put results into a NP.array.
2. Use a converter to HTML, CSV and PDF maybe. 
"""

# Now, we shall construct database for SPP lifetimes. Actually, SPP lifetime require the knowledge of all spectrum of response to be known. 
MaterialFolder="/usr/local/share/gsvit/data/spectra"
MaterialFile="Ag"

MaterialArray = loadtxt(MaterialFolder+'/'+MaterialFile, delimiter=' ', skiprows=4)

wavelengths=MaterialArray[:,0];
wavelengths=np.multiply(wavelengths,1e-6) #converting wavelength to emters

n=MaterialArray[:,1]; k=MaterialArray[:,2]

# Let's build epsilon
eps=np.add(n,np.multiply(1j, k))

plt.plot(wavelengths, eps)
#print MaterialArray
#print LifeTimeSpectrum(Material)