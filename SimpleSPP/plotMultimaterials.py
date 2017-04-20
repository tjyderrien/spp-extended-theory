#!/usr/bin/env python
#-*- coding: utf-8 -*-

# Copyright (C) 2013-2017 T. J.-Y. Derrien
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>

## @package plotMultimaterials
## Generate plots for many materials contained into a given database.
# 
# This allows to perform comparisons of interesting quantities between many materials, for specific wavelengths. 

# IMPORT LIBRARIES
from libSPP import *
from libPlotting import *
from libDatabase import *
from matplotlib.ticker import MaxNLocator
#EpsilonToIndex = np.vectorize(EpsilonToIndex)

#Swap = np.vectorize(Swap)

# =======================================================
if(len(sys.argv)<=1):
  print "==== MULTI-MATERIAL SIMULATION ==== "
  print "MODE 1: Single substrate - multi film mode."
  print "Usage:  ./plotMultimaterials.py           \ "
  print "        <Name ONE substrate (Air, Be, Au, ...)> \ "
  print "        <ONE wavelength (nm)>                          \ "
  print "        <Select source for data: Palik | Name of the 1st author + Year>"
  print "Example: ./plotMultimaterials.py Au 800 \"Johnson 1974\""
  #TODO: DISABLED PART FOR NOW.
  #print ""
  #print "MODE 2:  Multi-substrate - multi film, based on a given list."
  #print "Usage:   ./plotMultimaterials.py <FileName.dat> <Wavelength (nm)>"
  #print "Example: ./plotMultimaterials.py OxideList.dat"
  exit()

# FROM THIS POINT, WE KNOW THAT USER USED A COMMAND LINE ARGUMENTS. 

if any("dat" in s for s in sys.argv):
  print "** Detected USAGE 2."
else: 
  #(len(sys.argv[3]) > 0):
  print "** Detected USAGE 1."
  query = sys.argv[1]
  try:
    source = " ("+sys.argv[3]+")"
  except:
    source = ""
  query = query+source
  print "Selected substrate = "+query+"."
#else:
  #print "** Error: not planned case."

# Valid for any case
try: 
  wavelength = 1E-9*float(sys.argv[2])
except: 
  print "** Error: Please indicate the light wavelength."
  exit()
  

## Choose a wavelength
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

# TODO: if file OxideList.dat is provided, then we can look for couples, instead of generating the list of materials by ourselves. 

# Select the material of interface 1 #TODO: This selector may be not clear for users. 
if (not metal):
  SPPdb = FilterDatabase(SPPdb, query, 0)
else:
  SPPdb = FilterDatabase(SPPdb, query, 1)

SPPdbSave = SPPdb
print "** Filtering materials: SPP database has now "+str(len(SPPdb))+" entries."

if(len(SPPdb)==0):
  print "** QUITTING..."
  exit()
  
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

print "** Plotting Lspp as function of materials."
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
precision = 1.0e-1 #Precision for meshing (space: dielectric permittivities)

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
epsr = np.arange(-110e0,10e0, precision)
epsc = np.arange(0e0,50e0, precision)

eps2r, eps2c = np.meshgrid(epsr, epsc)

print "== Multimaterial mode: calculating period of SPP and their uncertainty..."
Period = (period(betaSPP(wavelength, 1e0+0e0j, eps2r+eps2c*1e0j))/lengthunit)
deltaPeriod = 1e9*(deltaPeriodSPP(wavelength, 1e0, eps2r+eps2c*1e0j, 0e0, 0e0, noise, noise))

print "** Plotting 2D period(Re eps, Im eps)..."

plt.figure()
levels = MaxNLocator(nbins=10).tick_values(0e0, 1E9*2.0*wavelength) #
if(reverse): #{
	CS = plt.contourf(eps1r, eps1c, Period, levels=levels, cmap=plt.cm.RdBu_r)
else:
	CS = plt.contourf(eps2r, eps2c, Period, levels=levels, cmap=plt.cm.RdBu_r)
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

print "** Plotting 2D DeltaPeriod(Re eps, Im eps)..."
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
#levels = [1, 5, 10, 20, 30, 40, 50] #um
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
plt.title(r'Fluctuations $\delta L_{SPP}$ ($\mu$m), $\lambda=$'+niceWavelength+' nm.')
plt.savefig('deltaLspp'+niceWavelength+'.eps')
plt.savefig('deltaLspp'+niceWavelength+'.png')
plt.show()
