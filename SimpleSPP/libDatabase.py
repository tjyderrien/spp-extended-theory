#!/usr/bin/env python
#-*- coding: utf-8 -*-
## @package libDatabase
# Functions to manage the material databases

import numpy as np
from numpy import genfromtxt, loadtxt, chararray

def FilterDatabase(SPPdb, query, FieldIndex):
  """ Filter SPP database using query and returns a smaller database
  /!\ content of query cell should be exact
  """
  SPPdbFiltered = np.array(SPPdb[SPPdb[:,FieldIndex]==query,:]) #uses a table of booleans to select
  return SPPdbFiltered

def FilterDatabaseContains(SPPdb, query, FieldIndex):
  """ Filter SPP database using query and returns a smaller database
  """
  SPPdbFiltered = SPPdb[np.array(np.core.defchararray.find(SPPdb[:,FieldIndex], query)==0),:]
  return SPPdbFiltered

def FilterDatabaseLowerThan(SPPdb, query, FieldIndex):
  """ Filter SPP database using query and returns a smaller database
  /!\ content of query cell should be exact
  """
  SPPdbFiltered = np.array(SPPdb[SPPdb[:,FieldIndex]<query,:]) #uses a table of booleans to select
  return SPPdbFiltered

## Unfold data from database of materials (5 columns)
def ExtractMaterialData(Database): #TODO: Think how to take data from continuous database directly instead of the ponctual file.
  Material1 = Database[:, 0]; BandGap = Database[:,1]; 
  Wavelength = Database[:, 2]; 
  RealEps = Database[:,3]; ImagEps = Database[:,4]; 
  
  #Converts strings to floats
  #Material1 = np.asfarray(Material1)
  #BandGap = np.asfarray(BandGap)
  #Wavelength = np.asfarray(Wavelength)
  #RealEps = np.asfarray(RealEps)
  #ImagEps = np.asfarray(ImagEps)
  
  return Material1, BandGap, Wavelength, RealEps, ImagEps

## Defines the interface with SPPactiveInterfaces.dat
def ExtractDataDb(SPPdbFiltered):
  # Extract data from database
  Material1 = SPPdbFiltered[:, 0]; Material2 = SPPdbFiltered[:,1]; 
  Wavelength = SPPdbFiltered[:, 2]; 
  OldSPPactiveBool = SPPdbFiltered[:,3]; NewSPPactiveBool = SPPdbFiltered[:,4]; 
  RealEps = SPPdbFiltered[:,5]; RealEpsError = SPPdbFiltered[:,6];
  SPPdecayDepth1 = SPPdbFiltered[:,7]; SPPdecayDepth2 = SPPdbFiltered[:,8]; 
  Reflectivity = SPPdbFiltered[:, 9]; OpticalPenetration1 = SPPdbFiltered[:,10]; OpticalPenetration2 = SPPdbFiltered[:,10]; SPPdecayLength = SPPdbFiltered[:,12]
  eps1r = SPPdbFiltered[:, 13]; eps1c = SPPdbFiltered[:,14]; eps2r = SPPdbFiltered[:,15]; eps2c = SPPdbFiltered[:,16]; k1imag = SPPdbFiltered[:,17]; 
  k2imag = SPPdbFiltered[:,18]
  DeltaLsppValue = SPPdbFiltered[:,19]
  
  #Converts strings to floats
  Wavelength = np.asfarray(Wavelength)
  RealEps = np.asfarray(RealEps)
  RealEpsError = np.asfarray(RealEpsError)
  SPPdecayDepth1 = np.asfarray(SPPdecayDepth1)
  SPPdecayDepth2 = np.asfarray(SPPdecayDepth2)
  Reflectivity = np.asfarray(Reflectivity)
  OpticalPenetration1 = np.asfarray(OpticalPenetration1)
  OpticalPenetration2 = np.asfarray(OpticalPenetration2)
  SPPdecayLength = np.asfarray(SPPdecayLength)
  eps1r = np.asfarray(eps1r)
  eps1c = np.asfarray(eps1c)
  eps2r = np.asfarray(eps2r)
  eps2c = np.asfarray(eps2c)
  k1imag = np.asfarray(k1imag)
  k2imag = np.asfarray(k2imag)
  DeltaLsppValue = np.asfarray(DeltaLsppValue)

  return Material1, Material2, Wavelength, OldSPPactiveBool, NewSPPactiveBool, RealEps, RealEpsError, SPPdecayDepth1, SPPdecayDepth2, Reflectivity, OpticalPenetration1, OpticalPenetration2, SPPdecayLength, eps1r, eps1c, eps2r, eps2c, k1imag, k2imag, DeltaLsppValue


## Export an SPP array to a CSV file
# SPP array must be produced with one of the SPPactiveInterfaces functions
#
def ExportToTxt(dbarray, filename):
  
  try: 
    np.savetxt(filename, dbarray, fmt="%s", delimiter='\t', newline='\n',comments='#')
    out = 0
  except: 
    print "Could not output SPP database into a file"
    out = 1
  
  #counter=0
  
  #for i in dbarray:
    #for k in dbarray:
      #if (np.mod(counter, 20) == 0):
	#show the table line each 20 lines, but also put it in a table
	#if (comment):
	#print '{0:30s} {1:30s} {2:15s} {3:12s} {4:12s} {5:11s} {6:16s} {7:16s} {8:12s} {9:19s} {10:19s} {11:15s}'.format("# Substrate", "Layer", "Wavelength (nm)", "OldSPPactive", "NewSPPactive", "Period (nm)", "DecayDepth1 (nm)", "DecayDepth2 (nm)", "Reflectivity", "OpticalPenetration1", "OpticalPenetration2", "DecayLength")
      #print dbarray[counter,:]
      #Material1 = i
      #counter=counter+1
  #print Material1
      #print '{0:30s} {1:30s} {2:15f} {3:12s} {4:12s} {5:11f} {6:16f} {7:16f} {8:12f}  {9:19f} {10:19f} {11:15f}'.format(Material1, Material2, Wavelength, OldSPPactiveBool, NewSPPactiveBool, RealEps, SPPdecayDepth1, SPPdecayDepth2, Reflectivity, OpticalPenetration1, OpticalPenetration2, SPPdecayLength)
  return out