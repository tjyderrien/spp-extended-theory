#!/usr/bin/env python
#-*- coding: utf-8 -*-

# IMPORT LIBRARIES
from SimpleSPProutines import *
#from makeTable import *

def ConvertToFloat(ndarray):
  """ As our program is badly written, 
  we must make our own function for casting to float
  """
  ndarray2 = ndarray
  for i in range(0,ndarray.size):
    #print i
    try: 
      ndarray2[i] = float(ndarray[i])
    except: 
      ndarray2[i] = 0.
  return ndarray2


# Get the database
SPPdb = GenerateDatabase()
# TODO: Remove empty lines. Difficult in ndarrays. 

# Extract the data to plot
Material1 = SPPdb[:, 0]; Material2 = SPPdb[:,1]; 
Wavelength = SPPdb[:, 2]; 
OldSPPactiveBool = SPPdb[:,3]; NewSPPactiveBool = SPPdb[:,4]; 
SPPperiod = SPPdb[:,5]; SPPdecayDepth1 = SPPdb[:,6]; SPPdecayDepth2 = SPPdb[:,7]; 
Reflectivity = SPPdb[:, 8]; OpticalPenetration1 = SPPdb[:,9]; OpticalPenetration2 = SPPdb[:,10]; SPPdecayLength = SPPdb[:,11]

# TODO: Casting numbers to float
# No method is properly working... 

#Wavelength[:]=float(Wavelength[:])


#Wavelength = ConvertToFloat(Wavelength); SPPperiod = ConvertToFloat(SPPperiod); 
#SPPdecayDepth1 = ConvertToFloat(SPPdecayDepth1); SPPdecayDepth2 = ConvertToFloat(SPPdecayDepth2)
#Reflectivity = ConvertToFloat(Reflectivity); OpticalPenetration1 = ConvertToFloat(OpticalPenetration1); OpticalPenetration2 = ConvertToFloat(OpticalPenetration2)
#SPPdecayLength = ConvertToFloat(SPPdecayLength)

#print type(float(Wavelength))


# Plot period as function of materials
#plt.figure()
#plt.xlabel('Material 1')
#plt.ylabel('Period (nm)')
#plt.plot(Material1, SPPperiod, '-')
#plt.legend(loc=2)
#plt.title('SPP period')
#plt.savefig('MultiMaterial_PeriodSPP.png')

### This was another possibility
# Load the data file
#datafile = "SPPactiveInterfaces.dat"

# Build vectors from this datafile
#SPParray = loadtxt(datafile, delimiter='\t', skiprows=1)

# 