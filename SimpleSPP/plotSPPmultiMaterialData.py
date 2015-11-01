#!/usr/bin/env python
#-*- coding: utf-8 -*-

# IMPORT LIBRARIES
from SimpleSPProutines import *
#from makeTable import *

# Make the database
SPPdb = GenerateDatabase()
print "SPP database has "+str(len(SPPdb))+" entries."

# How to filter database ? 
## Remove lines from SPPdb where some word is found ? 
query = 'Air'
filtermap = np.char.count(SPPdb, query) #generate (int) matrix of identified pattern iterations
#print occurences[lines,cols]
#filtervector = np.sum(filtermap, axis=1, keepdims=False)
#TODO how to copy the highest value of filtermap to the whole row ? 
np.array(filtermap.shape, np.max(filtermap))

#exit()
# building mask
#mask = np.array(filtermap, dtype=np.int) #matching cells are not null
#print mask
SPPdbFiltered = np.ma.masked_where(filtermap == 0, SPPdb) #matching cells are true
print SPPdbFiltered
SPPdbFiltered2 = np.ma.mask_rowcols(SPPdbFiltered, axis=0)
print SPPdbFiltered2
exit()
# identify line number where value is not 1
#print "Filtered data using '"+query+"' returned "+str(len(SPPdbFiltered))+" entries."

# Extract the data to plot
Material1 = SPPdbFiltered[:, 0]; Material2 = SPPdbFiltered[:,1]; 
Wavelength = SPPdbFiltered[:, 2]; 
OldSPPactiveBool = SPPdbFiltered[:,3]; NewSPPactiveBool = SPPdbFiltered[:,4]; 
SPPperiod = SPPdbFiltered[:,5]; SPPdecayDepth1 = SPPdbFiltered[:,6]; SPPdecayDepth2 = SPPdbFiltered[:,7]; 
Reflectivity = SPPdbFiltered[:, 8]; OpticalPenetration1 = SPPdbFiltered[:,9]; OpticalPenetration2 = SPPdbFiltered[:,10]; SPPdecayLength = SPPdbFiltered[:,11]
eps1r = SPPdbFiltered[:, 12]; eps1c = SPPdbFiltered[:,13]; eps2r = SPPdbFiltered[:,14]; eps2c = SPPdbFiltered[:,15]

#print Material1.size
Wavelength = np.asfarray(Wavelength)
SPPperiod = np.asfarray(SPPperiod)
SPPdecayDepth1 = np.asfarray(SPPdecayDepth1)
SPPdecayDepth2 = np.asfarray(SPPdecayDepth2)
Reflectivity = np.asfarray(Reflectivity)
OpticalPenetration1 = np.asfarray(OpticalPenetration1)
OpticalPenetration2 = np.asfarray(OpticalPenetration2)
SPPdecayLength = np.asfarray(SPPdecayLength)
eps1r = np.asfarray(eps1r)
eps1c = np.asfarray(eps1c)
eps1r = np.asfarray(eps2r)
eps2c = np.asfarray(eps2c)


#print type(SPPdecayLength[4])


##Plot period as function of materials
plt.figure()
plt.xlabel('Material 1')
plt.ylabel('Period (nm)')
plt.plot(eps1r, SPPperiod, '+', label='')
plt.legend(loc=2)
plt.title('SPP period')
plt.savefig('MultiMaterial_PeriodSPP.png')
plt.show()

### This was another possibility
# Load the data file
#datafile = "SPPactiveInterfaces.dat"

# Build vectors from this datafile
#SPParray = loadtxt(datafile, delimiter='\t', skiprows=1)

# 