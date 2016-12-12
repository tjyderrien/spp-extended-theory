#!/usr/bin/env python
#-*- coding: utf-8 -*-

# IMPORT LIBRARIES
from SimpleSPProutines import *
from advancedPlotting import *

from matplotlib.ticker import MaxNLocator

#from makeTable import *

def CleanStrArray(Material2): #{{{
  Material2clean = np.empty(Material2.shape, dtype='|S15')
  linenum=0
  for line in Material2: #for each line, replace Material2[line] with first word of Material2[line]
    fields = line.strip().split() #here is the first word, to replace the whole line. How to access id of line ?
    Material2clean[linenum] = fields[0]
    linenum = linenum + 1
  return(Material2clean)
#}}}
#print Material2clean


def plotDatabaseMaterials(database, legend, outputfile, query,metal): #{{{
  """plot period of SPP at various interfaces contained in a database
  """
  # Extract the data for 800 nm
  if (not metal):
    Material1, Material2, Wavelength, OldSPPactiveBool, NewSPPactiveBool, SPPperiod, SPPperiodError, SPPdecayDepth1, SPPdecayDepth2, Reflectivity, OpticalPenetration1, OpticalPenetration2, SPPdecayLength, eps1r, eps1c, eps2r, eps2c, k1imag, k2imag, DeltaLsppValue = ExtractDataDb(database)
  else: 
    Material2, Material1, Wavelength, OldSPPactiveBool, NewSPPactiveBool, SPPperiod, SPPperiodError, SPPdecayDepth1, SPPdecayDepth2, Reflectivity, OpticalPenetration1, OpticalPenetration2, SPPdecayLength, eps2r, eps2c, eps1r, eps1c, k1imag, k2imag, DeltaLsppValue = ExtractDataDb(database)
  
  # CLean the first field	
  Material2clean = CleanStrArray(Material2)
  
  # Prepare plot with arrows and text (but single wavelength)
  makePlot(eps2r, eps2c, Material2clean, outputfile, query, r'$Re(\varepsilon)$', r'$Im(\varepsilon)$', legend, 'r')
  #makePlot(eps2r, SPPdecayLength*1e-3, Material2clean, outputfile, query, r'$Re(\varepsilon)$', 'SPP decay length (um)', legend, 'r')

  return 0
#}}}

def plotDatabasePeriod(database, legend, outputfile, query, metal): #{{{
  """plot period of SPP at various interfaces contained in a database
  """
  # Extract the data for 800 nm
  if (not metal):
    Material1, Material2, Wavelength, OldSPPactiveBool, NewSPPactiveBool, SPPperiod, SPPperiodError, SPPdecayDepth1, SPPdecayDepth2, Reflectivity, OpticalPenetration1, OpticalPenetration2, SPPdecayLength, eps1r, eps1c, eps2r, eps2c, k1imag, k2imag, DeltaLsppValue = ExtractDataDb(database)
  else: 
    Material2, Material1, Wavelength, OldSPPactiveBool, NewSPPactiveBool, SPPperiod, SPPperiodError, SPPdecayDepth1, SPPdecayDepth2, Reflectivity, OpticalPenetration1, OpticalPenetration2, SPPdecayLength, eps2r, eps2c, eps1r, eps1c, k1imag, k2imag, DeltaLsppValue = ExtractDataDb(database)
  
  # CLean the first field	
  Material2clean = CleanStrArray(Material2)
  
  # Prepare plot with arrows and text (but single wavelength)
  #makePlot(eps2r, eps2c, Material2clean, outputfile, query, r'$Re(\varepsilon)$', r'$Im(\varepsilon)$', legend, 'r')
  makePlot(eps2r, SPPperiod, Material2clean, outputfile, query, r'$Re(\varepsilon)$', 'Period (nm)', legend, 'r')
  #makePlot(eps2r, SPPdecayLength*1e-3, Material2clean, outputfile, query, r'$Re(\varepsilon)$', 'SPP decay length (um)', legend, 'r')

  return 0
#}}}

def plotDatabaseLspp(database, legend, outputfile, query, metal): #{{{
  """plot period of SPP at various interfaces contained in a database
  """
  # Extract the data for 800 nm
  if (not metal):
    Material1, Material2, Wavelength, OldSPPactiveBool, NewSPPactiveBool, SPPperiod, SPPperiodError, SPPdecayDepth1, SPPdecayDepth2, Reflectivity, OpticalPenetration1, OpticalPenetration2, SPPdecayLength, eps1r, eps1c, eps2r, eps2c, k1imag, k2imag, DeltaLsppValue = ExtractDataDb(database)
  else: 
    Material2, Material1, Wavelength, OldSPPactiveBool, NewSPPactiveBool, SPPperiod, SPPperiodError, SPPdecayDepth1, SPPdecayDepth2, Reflectivity, OpticalPenetration1, OpticalPenetration2, SPPdecayLength, eps2r, eps2c, eps1r, eps1c, k1imag, k2imag, DeltaLsppValue = ExtractDataDb(database)
  
  # CLean the first field	
  Material2clean = CleanStrArray(Material2)
  
  # Prepare plot with arrows and text (but single wavelength)
  #makePlot(eps2r, eps2c, Material2clean, outputfile, query, r'$Re(\varepsilon)$', r'$Im(\varepsilon)$', legend, 'r')
  makePlot(eps2r, SPPdecayLength*1e-3, Material2clean, outputfile, query, r'$Re(\varepsilon)$', 'SPP decay length (um)', legend, 'r')

  return 0
#}}}

def plotDatabaseDeltaLspp(database, legend, outputfile, query, metal): #{{{
  """plot deltaLspp at various interfaces contained in a database
  """
  if (not metal):
    Material1, Material2, Wavelength, OldSPPactiveBool, NewSPPactiveBool, SPPperiod, SPPperiodError, SPPdecayDepth1, SPPdecayDepth2, Reflectivity, OpticalPenetration1, OpticalPenetration2, SPPdecayLength, eps1r, eps1c, eps2r, eps2c, k1imag, k2imag, deltaLsppValue  = ExtractDataDb(database)
  else: 
    Material2, Material1, Wavelength, OldSPPactiveBool, NewSPPactiveBool, SPPperiod, SPPperiodError, SPPdecayDepth1, SPPdecayDepth2, Reflectivity, OpticalPenetration1, OpticalPenetration2, SPPdecayLength, eps2r, eps2c, eps1r, eps1c, k1imag, k2imag, deltaLsppValue = ExtractDataDb(database)
  # CLean the first field	
  Material2clean = CleanStrArray(Material2)
  
  # Prepare plot with arrows and text (but single wavelength)
  #makePlot(eps2r, eps2c, Material2clean, outputfile, query, r'$Re(\varepsilon)$', r'$Im(\varepsilon)$', legend, 'r')
  makePlot(eps2r, deltaLsppValue, Material2clean, outputfile, query, r'$Re(\varepsilon)$', r'$\delta L_{SPP}$ (nm)', legend, 'r')

  return 0
#}}}

#EpsilonToIndex = np.vectorize(EpsilonToIndex)

#Swap = np.vectorize(Swap)

def plotSeveralWavelengths(database1, database2, reverse, metal, query):#{{{
  """
  Plot period as function of materials for two wavelengths. 
  Two types of data are included: 
  - Points: 1 point per material. 
  - Lines : the continuum calculation.
  
  Metal mode: enable for materials found metallic on most of the wavelengths
  Reverse mode: invert the material indices from the SPP database. Useful some metallic materials. 
  
  """
  # Extract data for 800 nm
  Material1, Material2, Wavelength, OldSPPactiveBool, NewSPPactiveBool, SPPperiod, SPPperiodError, SPPdecayDepth1, SPPdecayDepth2, Reflectivity, OpticalPenetration1, OpticalPenetration2, SPPdecayLength, eps1r, eps1c, eps2r, eps2c, k1imag, k2imag, DeltaLsppValue = ExtractDataDb(database1)
  
  if(metal): #swap eps1 and eps2
    Material2, Material1, Wavelength, OldSPPactiveBool, NewSPPactiveBool, SPPperiod, SPPperiodError, SPPdecayDepth2, SPPdecayDepth1, Reflectivity, OpticalPenetration2, OpticalPenetration1, SPPdecayLength, eps2r, eps2c, eps1r, eps1c, k2imag, k1imag, DeltaLsppValue = ExtractDataDb(database1)
  else: 
    Material1, Material2, Wavelength, OldSPPactiveBool, NewSPPactiveBool, SPPperiod, SPPperiodError, SPPdecayDepth1, SPPdecayDepth2, Reflectivity, OpticalPenetration1, OpticalPenetration2, SPPdecayLength, eps1r, eps1c, eps2r, eps2c, k1imag, k2imag, DeltaLsppValue = ExtractDataDb(database1)
    
  # Calculation of refractive index array
  eps1r=np.asfarray(eps1r)
  eps1c=np.asfarray(eps1c)
  eps2r=np.asfarray(eps2r)
  eps2c=np.asfarray(eps2c)

  epsilon1 = np.add(eps1r,np.multiply(1e0j, eps1c))
  epsilon2 = np.add(eps2r,np.multiply(1e0j, eps2c))
  
  fig1=plt.figure()
  if (not metal):
    plt.xlabel(r'Dielectric permittivity: $\mathcal{R}e(\varepsilon_2)$')
  else:
    plt.xlabel(r'Dielectric permittivity: $\mathcal{R}e(\varepsilon_1)$')
  
  # Discontinuous plotting
  plt.ylabel('SPP period $\Lambda$ (nm)')
  plt.plot(eps2r, SPPperiod, 'or', label='800 nm', markersize=8)
  
  # Continuous plotting
  eps1range = epsilon1[0] #use external index
  eps2range = np.arange(-70,0e0,0.1e0)
  
  # Use only one object
  # TODO: reprogram the whole function with switches in function arguments
  if(not metal): 
    refractiveindex1 = EpsilonToIndex(epsilon1)
  else: 
    refractiveindex1 = EpsilonToIndex(epsilon2)
    
  Radiation1=Wavelength/refractiveindex1.real
  Radiation1=np.sort(Radiation1)
  LambdaPMA=Wavelength[1]/np.sqrt(np.multiply(eps1range, eps2range)/(np.add(eps1range, eps2range)))
  
  plt.plot(np.sort(eps2r), Radiation1[::-1], 'r--')
  plt.plot(eps2range, LambdaPMA, 'r-')
  
  # Extract (again) for 400 nm
  #TODO: use a function here! code is repeated! 
  if(metal): #swap eps1 and eps2
    Material2, Material1, Wavelength, OldSPPactiveBool, NewSPPactiveBool, SPPperiod, SPPperiodError, SPPdecayDepth2, SPPdecayDepth1, Reflectivity, OpticalPenetration2, OpticalPenetration1, SPPdecayLength, eps2r, eps2c, eps1r, eps1c, k2imag, k1imag, DeltaLsppValue = ExtractDataDb(database2)
  else: 
    Material1, Material2, Wavelength, OldSPPactiveBool, NewSPPactiveBool, SPPperiod, SPPperiodError, SPPdecayDepth1, SPPdecayDepth2, Reflectivity, OpticalPenetration1, OpticalPenetration2, SPPdecayLength, eps1r, eps1c, eps2r, eps2c, k1imag, k2imag, DeltaLsppValue = ExtractDataDb(database2)
  
  # Calculation of refractive index
  eps1r=np.asfarray(eps1r)
  eps1c=np.asfarray(eps1c)
  eps2r=np.asfarray(eps2r)
  eps2c=np.asfarray(eps2c)

  epsilon1 = np.add(eps1r,np.multiply(1e0j, eps1c))
  epsilon2 = np.add(eps2r,np.multiply(1e0j, eps2c))
  
  # Continuous plots
  #if(metal):
    #eps1range = np.arange(-70,0e0,0.1e0)
    #eps2range = epsilon1[0] #TODO: use external index
  #else:
  eps1range = epsilon1[0]
  eps2range = np.arange(-70,0e0,0.1e0)
  
  if(not metal): 
    refractiveindex1 = EpsilonToIndex(epsilon1)
  else: 
    refractiveindex1 = EpsilonToIndex(epsilon2)
  
  #TODO: make a function here! Some code is repeated! 
  #This is the lambda/n continuous curve - n: sqrt(epsilon) of material1
  Radiation1=Wavelength/refractiveindex1.real
  Radiation1=np.sort(Radiation1)
  
  #This is rather to plot lambda_PMA
  LambdaPMA=Wavelength[1]/np.sqrt(np.multiply(eps1range, eps2range)/(np.add(eps1range, eps2range)))

  #plt.figure()
  plt.plot(eps2r, SPPperiod, 'bs', label='400 nm', markersize=8)
  plt.plot(np.sort(eps2r), Radiation1[::-1], 'b--')
  plt.plot(eps2range, LambdaPMA, 'b-') #TODO: plot using a full range, not eps2r

  #print eps1r
  #TODO simplify + combine the 3 following tests
  if(metal):
    if(query=="Au (Palik)"):
      plt.axis([0,25,0,900]) ##KEEP 900 please #good for Au
    elif(query=="Ti (Johnson)"):
      plt.axis([0,20,0,900]) ##KEEP 900 please #good for Ti
    else:
      print "** Error: this query is not a planned case. Query="+query
      #exit()
  else: #non-metal
    plt.axis([-70,0,0,900]) ##KEEP 900 please
    
  plt.axis()
  plt.yticks([0,200,400,600,800])
  #plt.axis([0,35,0,1000])
  if(not metal): 
    if(query=="Air"):
      #plt.legend(loc=4) #Good for Air
      plt.legend(bbox_to_anchor=(0.95, 0.05), loc=4, borderaxespad=0.)
    elif(query=="SiO2 (Palik)"):
      plt.legend(loc=2) #Good for SiO_2
    else:
      #print "** Error: this query is not a planned case. Query="+query
      exit()
  else: #metal case
    plt.legend(loc=1)
  #plt.legend(handler_map={line1: HandlerLine2D(numpoints=1)})
  #plt.legend(handler_map={line2: HandlerLine2D(numpoints=1)})
  #plt.title(query)
  plt.grid()
  
  # Adding a sub-plot
  if(not metal): 
    #a = plt.axes([-70,300,-40,700], axisbg='g')
    if(query=="SiO2 (Palik)"):
      a = plt.axes([0.2,0.18,0.35,0.35], axisbg='w') #Good for SiO2
      plt.axis([-6,0,250,290]) #Good for SiO2
      plt.yticks([250,270,290]) #Good for SiO2
    elif (query=="Air"):
      # Position of the plot
      a = plt.axes([0.2,0.17,0.35,0.35], axisbg='w') #Good for Air
      plt.axis([-6,0,380,410]) #Good for Air
      plt.yticks([380,390,400,410])
    else: #TODO: set a general case
      a = plt.axes([0.2,0.2,0.35,0.35], axisbg='w')

    plt.xticks([-6,-4,-2,0])
    
    plt.grid()
    plt.plot(eps2r, SPPperiod, 'bs', markersize=8)
    plt.plot(np.sort(eps2r), Radiation1[::-1], 'b--')
    #plt.title('Zoom')
    #plt.xticks([])
    #plt.yticks([])

  plt.savefig('MultiMaterial_PeriodSPP.eps')
  plt.show()
  return 0
#}}}
# =======================================================
if(len(sys.argv)<=2):
  print "Usage: ./plotMultimaterials.py           \ "
  print "    <Name of the substrate (Air, Be, Au, ...)> \ "
  print "    <wavelength (nm)>                          \ "
  print "    <Source for data: Palik or name of the 1st author>"
  print "Example: ./plotMultimaterials.py Au 800 Palik"
  exit()
  
query = sys.argv[1]
wavelength = 1E-9*float(sys.argv[2])
try:
  source = " ("+sys.argv[3]+")"
except:
  source = ""

query = query+source

## Choose a wavelength
#wavelength = 1030e-9
print "Selected substrate = "+query+"."
print "Operating wavelength = "+str(wavelength*1E9)+"nm."

## Choose which material to select
#query = 'Air'
#query = 'Au (Palik)'
#query = 'Ti (Palik)'
#query = 'Ti (Johnson 1974)'
#query= 'SiC (Palik?)'
#query = 'TiO2 (Devore 1951, e)'
#query = 'SiO2 (Palik)'

# TODO: THIS VARIABLE MUST BE NOT DEFINED BY HAND !!! O_O
reverse = False #reverse eps1 and eps2 for plotting
metal = False

thickness=500e-9
print "Interface 1: "+query

## Generate the database
SPPdb = GenerateDatabase()
print "** Info: SPP database has "+str(len(SPPdb))+" entries."

#print "Full Database:"
#print SPPdb

# Select the material of interface 1 #TODO: This selector may be not clear for users. 
if (not metal):
  SPPdb = FilterDatabase(SPPdb, query, 0)
else:
  SPPdb = FilterDatabase(SPPdb, query, 1)

SPPdbSave = SPPdb
print "** Filtering materials: SPP database has now "+str(len(SPPdb))+" entries."
#print "Filtering Material 1"
#print SPPdb

niceWavelength = str(int(wavelength*1e9))
selectWavelength = str(niceWavelength)+'.0'
title = niceWavelength+' nm'

try: 
  SPPdb = FilterDatabase(SPPdb, selectWavelength, 2)
  print "Filter on wavelength: SPP database "+title+" has "+str(len(SPPdb))+" entries."
  plotDatabasePeriod(SPPdb, title, 'Period'+niceWavelength+'nm.eps', title, metal)
except:
  print "Warning: no optical data is available for "+query+" at "+title+"."

# Plot the double plot for publication
try: 
  localWavelength = '800.0'
  SPPdb1 = FilterDatabase(SPPdbSave, localWavelength, 2)
  #print "Filter on wavelength: SPP database "+title+" has "+str(len(SPPdb1))+" entries."
  #plotDatabasePeriod(SPPdb1, '800 nm', 'Period'+localWavelength+'nm.eps', title, metal)
except:
  print "** Warning: no optical data is available for "+query+" at "+title+"."
  
try: 
  localWavelength = '400.0'
  SPPdb2 = FilterDatabase(SPPdbSave, localWavelength, 2)
  #print "Filter on wavelength: SPP database "+title+" has "+str(len(SPPdb1))+" entries."
  #plotDatabasePeriod(SPPdb2, '400 nm', 'Period'+localWavelength+'nm.eps', title, metal)
except:
  print "** Warning: no optical data is available for "+query+" at "+title+"."
    
# This line is for the paper figure. 
#plotSeveralWavelengths(SPPdb1, SPPdb2, reverse, metal, query)

# ===== PLOTTING the Lspp quantity as function of materials

print "Plotting Lspp as function of materials."
## ================== PLOT DATABASE ...
Material1, Material2, Wavelength, OldSPPactiveBool, NewSPPactiveBool, SPPperiod, SPPperiodError, SPPdecayDepth1, SPPdecayDepth2, Reflectivity, OpticalPenetration1, OpticalPenetration2, SPPdecayLength, eps1r, eps1c, eps2r, eps2c, k1imag, k2imag, DeltaLsppValue = ExtractDataDb(SPPdb)

#print k1imag, k2imag

#print "** Plotting the database of materials..."
plotDatabaseMaterials(SPPdb, title, 'Database'+title+'.eps', title, metal)


# =============== Output 2D plots [Re(eps), Im(eps)]

# Extract the material names from database
Material1, Material2, Wavelength, OldSPPactiveBool, NewSPPactiveBool, SPPperiod, SPPperiodError, SPPdecayDepth1, SPPdecayDepth2, Reflectivity, OpticalPenetration1, OpticalPenetration2, SPPdecayLength, eps1rM, eps1cM, eps2rM, eps2cM, k1imag, k2imag, DeltaLsppValue = ExtractDataDb(SPPdb)

# Calculation of refractive index
eps1rM=np.asfarray(eps1rM)
eps1cM=np.asfarray(eps1cM)
eps2rM=np.asfarray(eps2rM)
eps2cM=np.asfarray(eps2cM)

epsilon1 = np.add(eps1rM,np.multiply(1e0j, eps1cM))
epsilon2 = np.add(eps2rM,np.multiply(1e0j, eps2cM))

print "Plot the period 2D as function of dielectric permittivity"

noise = eta #Lack of knowledge on dielectric permittivity
precision = 0.25e0 #Precision for meshing

period = np.vectorize(period)
betaSPP = np.vectorize(betaSPP)
DecayLengthSPP = np.vectorize(DecayLengthSPP)
deltaBetaSPP = np.vectorize(deltaBetaSPP)
deltaPeriodSPP = np.vectorize(deltaPeriodSPP)
deltaLspp = np.vectorize(deltaLspp)

print "== Testing the basic functions"
print "Period = "+str(period(betaSPP(wavelength, 1e0+0e0j, -49.5738812793e0+3.81282698968e0j)))+" m"
print "DeltaBetaSPP = "+str(deltaBetaSPP(wavelength, 1e0+0e0j, -49.5738812793e0+3.81282698968e0j, 0e0, 0e0, noise, noise))+" m-1"
print "Decay Length SPP = "+str(DecayLengthSPP(betaSPP(wavelength, 1e0+0e0j, -49.5738812793e0+3.81282698968e0j)))+" m"
print "DeltaPeriodSPP = "+str(deltaPeriodSPP(wavelength, 1e0+0e0j, -49.5738812793e0+3.81282698968e0j, 0e0, 0e0, noise, noise))+" m"
print "DeltaLspp = "+str(deltaLspp(wavelength, 1e0+0e0j, -49.5738812793e0+3.81282698968e0j, 0e0, 0e0, noise, noise))+" m"
#exit()
print "== Knowledge over dielectric permittivity: +/- "+str(noise)+"."

print "Mesh generation..."
epsr = np.arange(-120e0,15e0, precision)
epsc = np.arange(0e0,60e0, precision)

eps2r, eps2c = np.meshgrid(epsr, epsc)

print "== Multimaterial mode: calculating period of SPP and their uncertainty..."
Period = (period(betaSPP(wavelength, 1e0+0e0j, eps2r+eps2c*1e0j))/lengthunit)
deltaPeriod = 1e9*(deltaPeriodSPP(wavelength, 1e0, eps2r+eps2c*1e0j, 0e0, 0e0, noise, noise))

print "** Plotting 2D period(Re eps, Im eps)..."

plt.figure()
levels = MaxNLocator(nbins=15).tick_values(0e0, Period.max())
if(reverse): #{
	CS = plt.contourf(eps1r, eps1c, Period, levels=levels, cmap=plt.cm.Blues)
else:
	CS = plt.contourf(eps2r, eps2c, Period, levels=levels, cmap=plt.cm.Blues)
#}
	
plt.xlabel(r'$Re(\varepsilon)$')
plt.ylabel(r'$Im(\varepsilon)$')

# adding the dots for materials
if(reverse): #{
	plt.plot(eps1rM, eps1cM, 'or', label=str(wavelength)+'nm', markersize=8)
else:
	plt.plot(eps2rM, eps2cM, 'or', label=str(wavelength)+'nm', markersize=8)
#}
#plt.clabel(CS, levels=levels, inline=False, fontsize=20)

#plt.clabel(CS, inline=1, fontsize=20)
#plt.legend(pos=1)
plt.colorbar(CS)
plt.title(r"Period $\Lambda_{SPP}$, $\lambda =$"+niceWavelength+" nm.")
plt.savefig('Period2d'+niceWavelength+'.eps')
plt.savefig('Period2d'+niceWavelength+'.png')
#plt.show()

##### "** Plotting 2D DeltaPeriod(Re eps, Im eps)..."
#levels = MaxNLocator(nbins=15).tick_values(0e0, deltaPeriod.max())
levels = [1, 5, 10, 20, 30, 40, 50, 100]

plt.figure()

if(reverse): #{
	CS = plt.contourf(eps1r, eps1c, deltaPeriod, levels=levels, cmap=plt.cm.Blues)
else:
	CS = plt.contourf(eps2r, eps2c, deltaPeriod, levels=levels, cmap=plt.cm.Blues)
#}
	
plt.xlabel(r'$Re(\varepsilon)$')
plt.ylabel(r'$Im(\varepsilon)$')

# adding the dots for materials
if(reverse): #{
	plt.plot(eps1rM, eps1cM, 'or', label=str(wavelength)+'nm', markersize=8)
else:
	plt.plot(eps2rM, eps2cM, 'or', label=str(wavelength)+'nm', markersize=8)
#}
#plt.clabel(CS, levels=levels, inline=False, fontsize=20)

#plt.clabel(CS, inline=1, fontsize=20)
#plt.legend(pos=1)
plt.colorbar(CS)
plt.title("")
plt.title(r"Period fluctuations $\delta \Lambda_{SPP}$ (nm), $\lambda =$"+niceWavelength+" nm.")
plt.savefig(query+'deltaPeriod'+niceWavelength+'.eps')
plt.savefig(query+'deltaPeriod'+niceWavelength+'.png')
#plt.show()

####################################
print "Plotting the Lspp as function of Re eps, Im eps."
plotDatabaseMaterials(SPPdb, 'Database (1030 nm)', 'MaterialDatabase'+str(wavelength)+'nm.eps', 'Database (1030 nm)', metal)
plotDatabaseLspp(SPPdb, title, 'Database'+title+'.eps', title, metal)
LsppTable = 1e6*DecayLengthSPP(betaSPP(wavelength, 1e0, eps2r+eps2c*1e0j))
deltaLsppTable = 1e6*(deltaLspp(wavelength, 1e0, eps2r+eps2c*1e0j, 0e0, 0e0, noise, noise))

plt.figure()
levels = [5, 10, 20, 30, 40, 50, 100, 200] #um
#levels = MaxNLocator(nbins=15).tick_values(0e0, LsppTable.max())
if(reverse): #{
	CS = plt.contourf(eps1r, eps1c, LsppTable, levels=levels, cmap=plt.cm.Blues)
else:
	CS = plt.contourf(eps2r, eps2c, LsppTable, levels=levels, cmap=plt.cm.Blues)
#}
	
plt.xlabel(r'$Re(\varepsilon)$')
plt.ylabel(r'$Im(\varepsilon)$')

# adding the materials information !
if(reverse): #{
	plt.plot(eps1rM, eps1cM, 'or', label=str(wavelength)+'nm', markersize=8)
else:
	plt.plot(eps2rM, eps2cM, 'or', label=str(wavelength)+'nm', markersize=8)
#}
#plt.clabel(CS, levels=levels, inline=False, fontsize=20)

#plt.clabel(CS, inline=1, fontsize=20)
#plt.legend(pos=1)
plt.colorbar(CS)
plt.title(r'SPP mean-free-path $L_{SPP}$ ($\mu$m), $\lambda=$'+niceWavelength+' nm.')
plt.savefig('Lspp2d'+niceWavelength+'.eps')
plt.savefig('Lspp2d'+niceWavelength+'.png')
#plt.show()

####### "Plotting the Lspp uncertainty due to dispersion over dielectric permittivity."
#print "deltaLsppTable = "+str(deltaLsppTable)

#print "delta Lspp min = "+str(deltaLsppTable.min())+", max = "+str(deltaLsppTable.max())+"."
levels = [1, 5, 10, 20, 30, 40, 50, 100, 200] #um
#levels = MaxNLocator(nbins=15).tick_values(0e0, deltaLsppTable.max())

plt.figure()

if(reverse): #{
	CS = plt.contourf(eps1r, eps1c, deltaLsppTable, levels=levels, cmap=plt.cm.Blues)
else:
	CS = plt.contourf(eps2r, eps2c, deltaLsppTable, levels=levels, cmap=plt.cm.Blues)
#}
	
plt.xlabel(r'$Re(\varepsilon)$')
plt.ylabel(r'$Im(\varepsilon)$')

# adding the materials information !
if(reverse): #{
	plt.plot(eps1rM, eps1cM, 'or', label=str(wavelength)+'nm', markersize=8)
else:
	plt.plot(eps2rM, eps2cM, 'or', label=str(wavelength)+'nm', markersize=8)
#}
#plt.clabel(CS, levels=levels, inline=False, fontsize=20)

#plt.clabel(CS, inline=1, fontsize=20)
#plt.legend(pos=1)
plt.colorbar(CS)
plt.title(r'SPP mean-free-path fluctuations $\delta L_{SPP}$ ($\mu$m), $\lambda=$'+niceWavelength+' nm.')
plt.savefig('deltaLspp'+niceWavelength+'.eps')
plt.savefig('deltaLspp'+niceWavelength+'.png')
plt.show()
