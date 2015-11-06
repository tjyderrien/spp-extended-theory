#!/usr/bin/env python
#-*- coding: utf-8 -*-

# IMPORT LIBRARIES
from SimpleSPProutines import *
from advancedPlotting import *

#from makeTable import *

def FilterDatabase(SPPdb, query, FieldIndex):
  """ Filter SPP database using query and returns a smaller database
  """
  SPPdbFiltered = np.array(SPPdb[SPPdb[:,FieldIndex]==query,:]) #uses a table of booleans to select
  return SPPdbFiltered


def CleanStrArray(Material2): 
  Material2clean = np.empty(Material2.shape, dtype='|S15')
  linenum=0
  for line in Material2: #for each line, replace Material2[line] with first word of Material2[line]
    fields = line.strip().split() #here is the first word, to replace the whole line. How to access id of line ?
    Material2clean[linenum] = fields[0]
    linenum = linenum + 1
  return(Material2clean)
#print Material2clean

def ExtractDataDb(SPPdbFiltered):
  # Extract data from database
  Material1 = SPPdbFiltered[:, 0]; Material2 = SPPdbFiltered[:,1]; 
  Wavelength = SPPdbFiltered[:, 2]; 
  OldSPPactiveBool = SPPdbFiltered[:,3]; NewSPPactiveBool = SPPdbFiltered[:,4]; 
  SPPperiod = SPPdbFiltered[:,5]; SPPdecayDepth1 = SPPdbFiltered[:,6]; SPPdecayDepth2 = SPPdbFiltered[:,7]; 
  Reflectivity = SPPdbFiltered[:, 8]; OpticalPenetration1 = SPPdbFiltered[:,9]; OpticalPenetration2 = SPPdbFiltered[:,10]; SPPdecayLength = SPPdbFiltered[:,11]
  eps1r = SPPdbFiltered[:, 12]; eps1c = SPPdbFiltered[:,13]; eps2r = SPPdbFiltered[:,14]; eps2c = SPPdbFiltered[:,15]
  
  #Converts strings to floats
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

  return Material1, Material2, Wavelength, OldSPPactiveBool, NewSPPactiveBool, SPPperiod, SPPdecayDepth1, SPPdecayDepth2, Reflectivity, OpticalPenetration1, OpticalPenetration2, SPPdecayLength, eps1r, eps1c, eps2r, eps2c


# Make the database
SPPdb = GenerateDatabase()
print "SPP database has "+str(len(SPPdb))+" entries."

# How to filter database ? 


## Choosing for which material 
#query = 'Air'
#query = 'Au (Johnson 1972)'
#query = 'Au (Palik)'
query = 'Ti (Palik)'
#query='SiC (Palik?)'
#query = 'TiO2 (Devore 1951, e)'
#query = 'SiO2 (Malitson 1965)'
SPPdb = FilterDatabase(SPPdb, query, 0)
SPPdb800 = FilterDatabase(SPPdb, '800.0', 2)
SPPdb400 = FilterDatabase(SPPdb, '400.0', 2)

print SPPdb800

# PLOT 1: 
# Extract the data for 800 nm
Material1, Material2, Wavelength, OldSPPactiveBool, NewSPPactiveBool, SPPperiod, SPPdecayDepth1, SPPdecayDepth2, Reflectivity, OpticalPenetration1, OpticalPenetration2, SPPdecayLength, eps1r, eps1c, eps2r, eps2c = ExtractDataDb(SPPdb800)

# CLean the first field
Material2clean = CleanStrArray(Material2)

# Prepare plot with arrows and text (but single wavelength)
makePlot(eps1r, SPPperiod, Material2clean, 'SPPperiodEnhanced800nm.eps', query, r'$Re(\varepsilon)$', 'Period (nm)', '800 nm', 'r')

# PLOT 2: 
# Extract data for 400 nm

# Prepare plot with arrow and text (for 400 nm wavelength)
Material1, Material2, Wavelength, OldSPPactiveBool, NewSPPactiveBool, SPPperiod, SPPdecayDepth1, SPPdecayDepth2, Reflectivity, OpticalPenetration1, OpticalPenetration2, SPPdecayLength, eps1r, eps1c, eps2r, eps2c = ExtractDataDb(SPPdb400)
# Clean the first field
Material2clean = CleanStrArray(Material2)
# make the plot2
makePlot(eps1r, SPPperiod, Material2clean, 'SPPperiodEnhanced400nm.eps', query, r'$Re(\varepsilon)$', 'Period (nm)', '400 nm', 'b')


###========================================================
##Plot period as function of materials for two wavelengths

# Extract (again) for 800 nm
Material1, Material2, Wavelength, OldSPPactiveBool, NewSPPactiveBool, SPPperiod, SPPdecayDepth1, SPPdecayDepth2, Reflectivity, OpticalPenetration1, OpticalPenetration2, SPPdecayLength, eps1r, eps1c, eps2r, eps2c = ExtractDataDb(SPPdb800)

# Calculation of effective refractive index: use vectorized function
epsilon1 = np.add(eps1r,np.multiply(1e0j, eps1c))
EpsilonToIndex = np.vectorize(EpsilonToIndex)
refractiveindex1 = EpsilonToIndex(epsilon1)
Radiation1=Wavelength/refractiveindex1.real
Radiation1=np.sort(Radiation1)

fig1=plt.figure()
plt.xlabel(r'$Re(\varepsilon)$')
plt.ylabel('Period (nm)')
plt.plot(eps1r, SPPperiod, 'or', label='800 nm', markersize=8)
plt.plot(np.sort(eps1r), Radiation1[::-1], 'r-', label=r'800 nm, $\lambda / n_1^{*}$')

# Extract (again) for 400 nm
Material1, Material2, Wavelength, OldSPPactiveBool, NewSPPactiveBool, SPPperiod, SPPdecayDepth1, SPPdecayDepth2, Reflectivity, OpticalPenetration1, OpticalPenetration2, SPPdecayLength, eps1r, eps1c, eps2r, eps2c = ExtractDataDb(SPPdb400)

# Calculation of effective refractive index: use vectorized function
epsilon1 = np.add(eps1r,np.multiply(1e0j, eps1c))
epsilon2 = np.add(eps2r,np.multiply(1e0j, eps2c))
#epsiloneff=np.divide(np.multiply(epsilon1, epsilon2), np.add(epsilon1, epsilon2))

EpsilonToIndex = np.vectorize(EpsilonToIndex)
#EffectiveIndex = np.vectorize(EffectiveIndex)

refractiveindex1 = EpsilonToIndex(epsilon1)
Radiation1=Wavelength/refractiveindex1
Radiation1=np.sort(Radiation1)

#plt.figure()
plt.plot(eps1r, SPPperiod, 'bs', label='400 nm', markersize=8)
plt.plot(np.sort(eps1r), Radiation1[::-1], 'b-', label=r'400 nm, $\lambda / n_1^{*}$')

#print eps1r

#plt.axis([-70,0,0,1000])
plt.axis([0,35,0,1000])
plt.legend(loc=1)
plt.title('SPP period ['+query+']')
plt.grid()
plt.savefig('MultiMaterial_PeriodSPP.eps')
plt.show()

### This was another possibility
# Load the data file
#datafile = "SPPactiveInterfaces.dat"

# Build vectors from this datafile
#SPParray = loadtxt(datafile, delimiter='\t', skiprows=1)

# 