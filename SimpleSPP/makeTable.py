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

#f = open('SPPactiveInterfaces.dat', 'w+')

# To calculate symmetric SPP compatible interfaces, use the following line
SPPactiveInterfacesArray = SPPactiveInterfaces(dbarray, '')

# To calculate asymmetric POSITIVE SPP compatible interfaces, use the following line
#SPPactiveInterfacesArray = AsymmetricSPPposActiveInterfaces(dbarray, '')

# To calculate asymmetric NEGATIVE SPP compatible interfaces, use the following line
#SPPactiveInterfacesArray = AsymmetricSPPnegActiveInterfaces(dbarray, '')

#f.write(SPPactiveInterfacesArray)

""" TODO: interface this with HTML for publication on the web. 
1. Put results into a NP.array.
2. Use a converter to HTML, CSV and PDF maybe. 
"""
