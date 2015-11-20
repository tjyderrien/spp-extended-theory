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
  SPPperiod = SPPdbFiltered[:,5]; SPPperiodError = SPPdbFiltered[:,6];
  SPPdecayDepth1 = SPPdbFiltered[:,7]; SPPdecayDepth2 = SPPdbFiltered[:,8]; 
  Reflectivity = SPPdbFiltered[:, 9]; OpticalPenetration1 = SPPdbFiltered[:,10]; OpticalPenetration2 = SPPdbFiltered[:,10]; SPPdecayLength = SPPdbFiltered[:,12]
  eps1r = SPPdbFiltered[:, 13]; eps1c = SPPdbFiltered[:,14]; eps2r = SPPdbFiltered[:,15]; eps2c = SPPdbFiltered[:,16]; k1imag = SPPdbFiltered[:,17]; 
  k2imag = SPPdbFiltered[:,18]
  
  #Converts strings to floats
  Wavelength = np.asfarray(Wavelength)
  SPPperiod = np.asfarray(SPPperiod)
  SPPperiodError = np.asfarray(SPPperiodError)
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

  return Material1, Material2, Wavelength, OldSPPactiveBool, NewSPPactiveBool, SPPperiod, SPPperiodError, SPPdecayDepth1, SPPdecayDepth2, Reflectivity, OpticalPenetration1, OpticalPenetration2, SPPdecayLength, eps1r, eps1c, eps2r, eps2c, k1imag, k2imag

def plotDatabasePeriod(database, legend, outputfile): 
  """plot period of SPP at various interfaces contained in a database
  """
  # Extract the data for 800 nm
  Material1, Material2, Wavelength, OldSPPactiveBool, NewSPPactiveBool, SPPperiod, SPPperiodError, SPPdecayDepth1, SPPdecayDepth2, Reflectivity, OpticalPenetration1, OpticalPenetration2, SPPdecayLength, eps1r, eps1c, eps2r, eps2c, k1imag, k2imag = ExtractDataDb(database)
  
  # CLean the first field
  Material2clean = CleanStrArray(Material2)
  
  # Prepare plot with arrows and text (but single wavelength)
  makePlot(eps2r, SPPperiod, Material2clean, outputfile, query, r'$Re(\varepsilon)$', 'Period (nm)', legend, 'r')

  return 0

EpsilonToIndex = np.vectorize(EpsilonToIndex)
EffectiveIndex = np.vectorize(EffectiveIndex)

def Swap(eps1, eps2):
  eps3 = eps1
  eps1 = eps2
  eps2 = eps3
  del eps3
  return(eps1, eps2) 

Swap = np.vectorize(Swap)

def plotSeveralWavelengths(database1, database2, reverse, metal):
  """
  Plot period as function of materials for two wavelengths
  """
  # Extract data for 800 nm
  Material1, Material2, Wavelength, OldSPPactiveBool, NewSPPactiveBool, SPPperiod, SPPperiodError, SPPdecayDepth1, SPPdecayDepth2, Reflectivity, OpticalPenetration1, OpticalPenetration2, SPPdecayLength, eps1r, eps1c, eps2r, eps2c, k1imag, k2imag = ExtractDataDb(database1)
  
  if(reverse): #swap eps1 and eps2
    Material2, Material1, Wavelength, OldSPPactiveBool, NewSPPactiveBool, SPPperiod, SPPperiodError, SPPdecayDepth2, SPPdecayDepth1, Reflectivity, OpticalPenetration2, OpticalPenetration1, SPPdecayLength, eps2r, eps2c, eps1r, eps1c, k2imag, k1imag = ExtractDataDb(database1)
  else: 
    Material1, Material2, Wavelength, OldSPPactiveBool, NewSPPactiveBool, SPPperiod, SPPperiodError, SPPdecayDepth1, SPPdecayDepth2, Reflectivity, OpticalPenetration1, OpticalPenetration2, SPPdecayLength, eps1r, eps1c, eps2r, eps2c, k1imag, k2imag = ExtractDataDb(database1)
    
  # Calculation of refractive index
  eps1r=np.asfarray(eps1r)
  eps1c=np.asfarray(eps1c)
  eps2r=np.asfarray(eps2r)
  eps2c=np.asfarray(eps2c)

  epsilon1 = np.add(eps1r,np.multiply(1e0j, eps1c))
  epsilon2 = np.add(eps2r,np.multiply(1e0j, eps2c))
  
  if(not metal): 
    refractiveindex1 = EpsilonToIndex(epsilon1)
  else: 
    refractiveindex1 = EpsilonToIndex(epsilon2)
    
  Radiation1=Wavelength/refractiveindex1.real
  Radiation1=np.sort(Radiation1)

  fig1=plt.figure()
  plt.xlabel(r'Dielectric permittivity: $\mathcal{R}e(\varepsilon_2)$')
  plt.ylabel('SPP period $\Lambda$ (nm)')
  plt.plot(eps2r, SPPperiod, 'or', label='800 nm', markersize=8)
  plt.plot(np.sort(eps2r), Radiation1[::-1], 'r-')

  # Extract (again) for 400 nm
  
  if(reverse): #swap eps1 and eps2
    Material2, Material1, Wavelength, OldSPPactiveBool, NewSPPactiveBool, SPPperiod, SPPperiodError, SPPdecayDepth2, SPPdecayDepth1, Reflectivity, OpticalPenetration2, OpticalPenetration1, SPPdecayLength, eps2r, eps2c, eps1r, eps1c, k2imag, k1imag = ExtractDataDb(database2)
  else: 
    Material1, Material2, Wavelength, OldSPPactiveBool, NewSPPactiveBool, SPPperiod, SPPperiodError, SPPdecayDepth1, SPPdecayDepth2, Reflectivity, OpticalPenetration1, OpticalPenetration2, SPPdecayLength, eps1r, eps1c, eps2r, eps2c, k1imag, k2imag = ExtractDataDb(database2)
  
  # Calculation of refractive index
  eps1r=np.asfarray(eps1r)
  eps1c=np.asfarray(eps1c)
  eps2r=np.asfarray(eps2r)
  eps2c=np.asfarray(eps2c)

  epsilon1 = np.add(eps1r,np.multiply(1e0j, eps1c))
  epsilon2 = np.add(eps2r,np.multiply(1e0j, eps2c))

  if(not metal): 
    refractiveindex1 = EpsilonToIndex(epsilon1)
  else: 
    refractiveindex1 = EpsilonToIndex(epsilon2)
    
  Radiation1=Wavelength/refractiveindex1
  Radiation1=np.sort(Radiation1)

  #plt.figure()
  plt.plot(eps2r, SPPperiod, 'bs', label='400 nm', markersize=8)
  plt.plot(np.sort(eps2r), Radiation1[::-1], 'b-')

  #print eps1r
  if(metal):
    plt.axis([0,40,0,900]) ##KEEP 900 please
  else:
    plt.axis([-70,0,0,900]) ##KEEP 900 please
  plt.axis()
  plt.yticks([0,200,400,600,800])
  #plt.axis([0,35,0,1000])
  if(not metal): 
    plt.legend(loc=4)
  else: 
    plt.legend(loc=1)
  #plt.title(query)
  plt.grid()
  
  if(not metal): 
    #a = plt.axes([-70,300,-40,700], axisbg='g')
    #a = plt.axes([0.2,0.17,0.35,0.35], axisbg='w') #Good for SiO2
    a = plt.axes([0.2,0.2,0.35,0.35], axisbg='w') #Good for Air

    plt.xticks([-6,-4,-2,0])
    #plt.axis([-6,0,250,290]) #Good for SiO2
    plt.axis([-6,0,380,410]) #Good for Air
    #plt.yticks([250,270,290]) #Good for SiO2
    plt.yticks([380,390,400,410]) #Good for Air
    plt.grid()
    plt.plot(eps2r, SPPperiod, 'bs', markersize=8)
    plt.plot(np.sort(eps2r), Radiation1[::-1], 'b-')
    #plt.title('Zoom')
    #plt.xticks([])
    #plt.yticks([])

  plt.savefig('MultiMaterial_PeriodSPP.eps')
  plt.show()
  return 0

# =======================================================

## Choosing for which material 
query = 'Air'
#query = 'Au (Johnson 1972)'
#query = 'Au (Palik)'
#query = 'Ti (Palik)'
#query= 'SiC (Palik?)'
#query = 'TiO2 (Devore 1951, e)'
#query = 'SiO2 (Malitson 1965)'

# Make the database
SPPdb = GenerateDatabase()
print "SPP database has "+str(len(SPPdb))+" entries."

print "Full Database:"
print SPPdb

# Select the material of interface 1
SPPdb = FilterDatabase(SPPdb, query, 0)
print "Filter on materials: SPP database has now "+str(len(SPPdb))+" entries."
#print "Filtering Material 1"
#print SPPdb

SPPdb1030 = FilterDatabase(SPPdb, '1030.0', 2)
print "Filter on wavelength: SPP database 1030 nm has "+str(len(SPPdb1030))+" entries."

SPPdb800 = FilterDatabase(SPPdb, '800.0', 2)
print "Filter on wavelength: SPP database 800 nm has "+str(len(SPPdb800))+" entries."

SPPdb400 = FilterDatabase(SPPdb, '400.0', 2)
print "Filter on wavelength: SPP database 400 nm has "+str(len(SPPdb400))+" entries."

plotDatabasePeriod(SPPdb1030, '1030 nm', 'SPPperiodEnhanced1030nm.eps')
plotDatabasePeriod(SPPdb800, '800 nm', 'SPPperiodEnhanced800nm.eps')
plotDatabasePeriod(SPPdb400, '400 nm', 'SPPperiodEnhanced400nm.eps')

reverse = False #reverse eps1 and eps2 for plotting
metal = False
plotSeveralWavelengths(SPPdb800, SPPdb400, reverse, metal)


# ================== Check the k1imag, k2imag signs...
Material1, Material2, Wavelength, OldSPPactiveBool, NewSPPactiveBool, SPPperiod, SPPperiodError, SPPdecayDepth1, SPPdecayDepth2, Reflectivity, OpticalPenetration1, OpticalPenetration2, SPPdecayLength, eps1r, eps1c, eps2r, eps2c, k1imag, k2imag = ExtractDataDb(SPPdb800)

print k1imag, k2imag

## Plot the database materials of 800 nm

#print eps2r.shape, eps2c.shape

plt.figure()
plt.title('Materials of database at 800 nm')
plt.xlabel(r'$Re (\epsilon)$')
plt.ylabel(r'$Im (\epsilon)$')
plt.plot(eps2r, eps2c, 'rs')
plt.show()
plt.savefig('Database.eps')
