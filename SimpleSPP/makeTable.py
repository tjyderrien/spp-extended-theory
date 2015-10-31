#!/usr/bin/env python
#-*- coding: utf-8 -*-

# IMPORT LIBRARIES
from SimpleSPProutines import *
#from plotGraph import *

# PHYSICAL INPUT
#wavelength = 800e-9
#epsAir=1e0
#epsSi0=13.64+0.048j
#meffe=0.18
#nuSi=(1.1e-15)**-1

#example=period(betaSPP(wavelength, epsAir, Drude(wavelength, 1e28, epsSi0, nuSi)))

#print example

SPPdb = GenerateDatabase()
SppOutput = 'SPPactiveInterfaces.dat'
ExportToTxt(SPPdb, SppOutput)

""" TODO: interface with HTML for publication on the web. 
1. Put results into a NP.array.
2. Use a converter to HTML, CSV and PDF maybe. 
"""

# Plotting the LIPSS period as function of materials

# Irradiation in air environment
# Concept: In the generated file with SPP active interfaces, plot the period as function of something

# SPPactiveInterfacesArray contains all data we need, just remove lines starting with #. 

#MaterialToSelect="Air"

#plt.figure()
#plt.xlabel('SPP decay length')
#plt.ylabel('SPP period (nm)')

#for i in SPPactiveInterfacesArray:
	#if (SPPactiveInterfacesArray[i,0]==MaterialToSelect): #select only entries which corresponds to interface of interest
		## Construct the interesting data point with x: decay length, y: period, txt: name of interface
		##TODO x = append()
		
##plt.plot(,y,label='')
		

