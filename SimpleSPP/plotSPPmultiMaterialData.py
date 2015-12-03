#!/usr/bin/env python
#-*- coding: utf-8 -*-

# IMPORT LIBRARIES
from SimpleSPProutines import *
from advancedPlotting import *

from matplotlib.ticker import MaxNLocator

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

def plotDatabaseMaterials(database, legend, outputfile, query): 
  """plot period of SPP at various interfaces contained in a database
  """
  # Extract the data for 800 nm
  Material1, Material2, Wavelength, OldSPPactiveBool, NewSPPactiveBool, SPPperiod, SPPperiodError, SPPdecayDepth1, SPPdecayDepth2, Reflectivity, OpticalPenetration1, OpticalPenetration2, SPPdecayLength, eps1r, eps1c, eps2r, eps2c, k1imag, k2imag = ExtractDataDb(database)
  
  # CLean the first field	
  Material2clean = CleanStrArray(Material2)
  
  # Prepare plot with arrows and text (but single wavelength)
  makePlot(eps2r, eps2c, Material2clean, outputfile, query, r'$Re(\varepsilon)$', r'$Im(\varepsilon)$', legend, 'r')
  #makePlot(eps2r, SPPdecayLength*1e-3, Material2clean, outputfile, query, r'$Re(\varepsilon)$', 'SPP decay length (um)', legend, 'r')

  return 0

def plotDatabasePeriod(database, legend, outputfile, query): 
  """plot period of SPP at various interfaces contained in a database
  """
  # Extract the data for 800 nm
  Material1, Material2, Wavelength, OldSPPactiveBool, NewSPPactiveBool, SPPperiod, SPPperiodError, SPPdecayDepth1, SPPdecayDepth2, Reflectivity, OpticalPenetration1, OpticalPenetration2, SPPdecayLength, eps1r, eps1c, eps2r, eps2c, k1imag, k2imag = ExtractDataDb(database)
  
  # CLean the first field	
  Material2clean = CleanStrArray(Material2)
  
  # Prepare plot with arrows and text (but single wavelength)
  makePlot(eps2r, SPPperiod, Material2clean, outputfile, query, r'$Re(\varepsilon)$', 'Period (nm)', legend, 'r')
  makePlot(eps2r, SPPdecayLength*1e-3, Material2clean, outputfile, query, r'$Re(\varepsilon)$', 'SPP decay length (um)', legend, 'r')

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
  
  fig1=plt.figure()
  plt.xlabel(r'Dielectric permittivity: $\mathcal{R}e(\varepsilon_2)$')
  plt.ylabel('SPP period $\Lambda$ (nm)')
  
  plt.plot(eps2r, SPPperiod, 'or', label='800 nm', markersize=8)
  
  eps2range = np.arange(0e0,40e0,0.1e0)
  
  if(not metal): 
    refractiveindex1 = EpsilonToIndex(epsilon1.real)
  else: 
    refractiveindex1 = EpsilonToIndex(epsilon2.real)
    #refractiveindex1 = EpsilonToIndex(eps2range)
    #Radiation1=Wavelength/refractiveindex1.real
    #Radiation1=np.sort(Radiation1)
    #plt.plot(eps2range, Radiation1, 'r-')
    
  Radiation1=Wavelength/refractiveindex1.real
  Radiation1=np.sort(Radiation1)
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
    plt.legend(loc=2)
  else: 
    plt.legend(loc=1)
  #plt.legend(handler_map={line1: HandlerLine2D(numpoints=1)})
  #plt.legend(handler_map={line2: HandlerLine2D(numpoints=1)})
  #plt.title(query)
  plt.grid()
  
  if(not metal): 
    #a = plt.axes([-70,300,-40,700], axisbg='g')
    #a = plt.axes([0.2,0.17,0.35,0.35], axisbg='w') #Good for SiO2
    #a = plt.axes([0.2,0.2,0.35,0.35], axisbg='w') #Good for Air

    plt.xticks([-6,-4,-2,0])
    plt.axis([-6,0,250,290]) #Good for SiO2
    #plt.axis([-6,0,380,410]) #Good for Air
    plt.yticks([250,270,290]) #Good for SiO2
    #plt.yticks([380,390,400,410]) #Good for Air
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

## Choose which material to select
query = 'Air'
#query = 'Au (Johnson 1972)'
#query = 'Au (Palik)'
#query = 'Ti (Palik)'
#query= 'SiC (Palik?)'
#query = 'TiO2 (Devore 1951, e)'
#query = 'SiO2 (Palik)'

## Generate the database
SPPdb = GenerateDatabase()
print "SPP database has "+str(len(SPPdb))+" entries."

print "Full Database:"
print SPPdb

# Select the material of interface 1
SPPdb = FilterDatabase(SPPdb, query, 0)
print "Filter on materials: SPP database has now "+str(len(SPPdb))+" entries."
#print "Filtering Material 1"
#print SPPdb


try: 
  title = '1064 nm'
  SPPdb1064 = FilterDatabase(SPPdb, '1064.0', 2)
  print "Filter on wavelength: SPP database "+title+" has "+str(len(SPPdb1064))+" entries."
  plotDatabasePeriod(SPPdb1064, title, 'SPPperiodEnhanced1064nm.eps', title)
except:
  print "Exception: no optical data is available for "+query+" at "+title+"."

try: 
  title = '1030 nm'
  SPPdb1030 = FilterDatabase(SPPdb, '1030.0', 2)
  print "Filter on wavelength: SPP database "+title+" has "+str(len(SPPdb1030))+" entries."
  plotDatabasePeriod(SPPdb1030, title, 'SPPperiodEnhanced1030nm.eps', title)
except:
  print "Exception: no optical data is available for "+query+" at "+title+"."

try: 
  title = '800 nm'
  SPPdb800 = FilterDatabase(SPPdb, '800.0', 2)
  print "Filter on wavelength: SPP database "+title+" has "+str(len(SPPdb800))+" entries."
  plotDatabasePeriod(SPPdb800, title, 'SPPperiodEnhanced800nm.eps', title)
except: 
  print "Exception: no optical data is available for "+query+" at "+title+"."

try: 
  title = '532 nm'
  SPPdb532 = FilterDatabase(SPPdb, '532.0', 2)
  print "Filter on wavelength: SPP database "+title+" has "+str(len(SPPdb532))+" entries."
  plotDatabasePeriod(SPPdb532, title, 'SPPperiodEnhanced532nm.eps', title)
except: 
  print "Exception: no optical data is available for "+query+" at "+title+"."

try:
  title = '400 nm'
  SPPdb400 = FilterDatabase(SPPdb, '400.0', 2)
  print "Filter on wavelength: SPP database "+title+" has "+str(len(SPPdb400))+" entries."
  plotDatabasePeriod(SPPdb400, title, 'SPPperiodEnhanced400nm.eps', title)
except: 
  print "Exception: no optical data is available for "+query+" at "+title+"."

reverse = False #reverse eps1 and eps2 for plotting
metal = False

plotSeveralWavelengths(SPPdb800, SPPdb400, reverse, metal)

# ===== PLOTTING the Lspp quantity as function of materials

## ================== PLOT DATABASE ...
Material1, Material2, Wavelength, OldSPPactiveBool, NewSPPactiveBool, SPPperiod, SPPperiodError, SPPdecayDepth1, SPPdecayDepth2, Reflectivity, OpticalPenetration1, OpticalPenetration2, SPPdecayLength, eps1r, eps1c, eps2r, eps2c, k1imag, k2imag = ExtractDataDb(SPPdb1030)

print k1imag, k2imag

## Plot the database materials of 1030 nm
title = '1030 nm'
plotDatabaseMaterials(SPPdb1030, title, 'Database.eps', title)
#print eps2r.shape, eps2c.shape

#plt.figure()
#plt.title('Materials of database at 1030 nm')
#plt.xlabel(r'$Re (\epsilon)$')
#plt.ylabel(r'$Im (\epsilon)$')
#plt.plot(eps2r, eps2c, 'rs')
#plt.show()
#plt.savefig('Database.eps')

# =============== Output 2D plot delta(periodSPP) [Re(eps), Im(eps)]

# Extract the material names from database
Material1, Material2, Wavelength, OldSPPactiveBool, NewSPPactiveBool, SPPperiod, SPPperiodError, SPPdecayDepth1, SPPdecayDepth2, Reflectivity, OpticalPenetration1, OpticalPenetration2, SPPdecayLength, eps1rM, eps1cM, eps2rM, eps2cM, k1imag, k2imag = ExtractDataDb(SPPdb1030)

# Calculation of refractive index
eps1rM=np.asfarray(eps1rM)
eps1cM=np.asfarray(eps1cM)
eps2rM=np.asfarray(eps2rM)
eps2cM=np.asfarray(eps2cM)

epsilon1 = np.add(eps1rM,np.multiply(1e0j, eps1cM))
epsilon2 = np.add(eps2rM,np.multiply(1e0j, eps2cM))

print "Plot the uncertainty on period as function of dielectric permittivity"

noise = 1e0

deltaBetaSPP = np.vectorize(deltaBetaSPP)
deltaPeriodSPP = np.vectorize(deltaPeriodSPP)

print deltaBetaSPP(800e-9, 1e0+0e0j, 1.1e0+1.1e0j, noise, noise, noise, noise)
print deltaPeriodSPP(800e-9, 1e0+0e0j, 1.1e0+1.1e0j, noise, noise, noise, noise)

precision = 1e-1

print "Mesh generation..."
epsr = np.arange(-5e0,5e0, precision)
epsc = np.arange(0e0,5e0, precision)

eps2r, eps2c = np.meshgrid(epsr, epsc)

print "Calculating uncertainty on SPP period..."
deltaPeriod = 1e9*(deltaPeriodSPP(1030e-9, 1e0, eps2r+eps2c*1e0j, 0e0, 0e0, noise, noise))
print "delta Period min = "+str(deltaPeriod.min())+", max = "+str(deltaPeriod.max())+"."
#levels = MaxNLocator(nbins=15).tick_values(0e0, deltaPeriod.max())
levels = [5, 10, 50, 100, 200, 300, 400, 500]

plt.figure()
CS = plt.contourf(eps2r, eps2c, deltaPeriod, levels=levels, cmap=plt.cm.Blues)
plt.xlabel(r'$Re(\varepsilon)$')
plt.ylabel(r'$Im(\varepsilon)$')

# adding the materials information !
plt.plot(eps2rM, eps2cM, 'or', label='1030 nm', markersize=8)
#plt.clabel(CS, levels=levels, inline=False, fontsize=20)

#plt.clabel(CS, inline=1, fontsize=20)
#plt.legend(pos=1)
plt.colorbar(CS)
plt.savefig('deltaPeriod.eps')
plt.savefig('deltaPeriod.png')
plt.show()