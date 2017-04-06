#!/usr/bin/env python
#-*- coding: utf-8 -*-

## Copyright (C) 2013-2017 T. J.-Y. Derrien
##
## This program is free software: you can redistribute it and/or modify
## it under the terms of the GNU General Public License as published by
## the Free Software Foundation, either version 3 of the License, or
## (at your option) any later version.
##
## This program is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.
##
## You should have received a copy of the GNU General Public License
## along with this program.  If not, see <http://www.gnu.org/licenses/>

## @package importPalikData
# Importing data from Palik book using digitized plots. 
# Optical data are input via using CSV-formatted input files, created using Engauge-digitizer software. 
# The lib generates (wavelength, ReEps, ImEps) tables to be used inside the program.

## Imports data (wavelength, n,k) from different wavelength meshes. 
#  Such data can be captured using a software like Engauge Digitizer. 
#  This leads to obtain (n,k) discretized on DIFFERENT MESHES. 
#  
#  NOTE: Warning: the data produced by this method are rather unprecise. SPP spectoscopy requires precision to 1E-3. 
#  This method gives a precision worst then 1E0. Use only in case no other data are available. 

# IMPORT LIBRARIES
from libSPP import *

def importFromNKtable(wavelength, folder, filename, plotting=1):#{{{
  nfile = folder+filename+"-n.csv"
  kfile = folder+filename+"-k.csv"

  narray = loadtxt(nfile, delimiter="\t", skiprows=1)
  karray = loadtxt(kfile, delimiter="\t", skiprows=1)
  
  unit = 1E-6
  
  # import wavelength, n and k from Palik
  wavelength1 = narray[:,0]*unit
  wavelength2 = karray[:,0]*unit
  n = narray[:,1]
  k = karray[:,1]
  numrows = 10000
  base = 10
  # interpolate n and k on new wavelength mesh
  order=1
  #wavelengths = np.arange(np.amin(wavelength2),np.amax(wavelength2), precision) #regular mesh, AWFUL for memory
  print "Generating new wavelength mesh: ("+str(np.amin(wavelength2))+", "+str(np.amax(wavelength2))+")"
  wavelengths = np.logspace(np.amin(np.log10(wavelength2)), np.amax(np.log10(wavelength2)), num=numrows, base=base, endpoint = True)

  print "New wavelength mesh has "+str(numrows)+" rows."
  #print wavelengths
  
  fni = InterpolatedUnivariateSpline(wavelength1, n, k=order)
  fki = InterpolatedUnivariateSpline(wavelength2, k, k=order)

  #Interpolated one optical constants
  try: 
	#wavelength = 800e-9
	ni = fni(wavelength); ki = fki(wavelength)
	epsilon = (ni+1j*ki)**2
	print "Interpolated permittivity at "+str(wavelength*1e9)+" nm = "+str(epsilon)
  except: 
	  print "Interpolation for "+str(wavelength*1E9)+" nm failed."
	  
  #try:
	#wavelength = 532e-9
	#ni = fni(wavelength); ki = fki(wavelength)
	#epsilon = (ni+1j*ki)**2
	#print "Interpolated permittivity at "+str(wavelength*1e9)+" nm = "+str(epsilon)
  #except: 
	  #print "Interpolation for "+str(wavelength*1E9)+" nm failed."
  #try:
	#wavelength = 400e-9
	#ni = fni(wavelength); ki = fki(wavelength)
	#epsilon = (ni+1j*ki)**2
	#print "Interpolated permittivity at "+str(wavelength*1e9)+" nm = "+str(epsilon)
  #except: 
	  #print "Interpolation for "+str(wavelength*1E9)+" nm failed."
  #try:
	#wavelength = 930e-9
	#ni = fni(wavelength); ki = fki(wavelength)
	#epsilon = (ni+1j*ki)**2
	#print "Interpolated permittivity at "+str(wavelength*1e9)+" nm = "+str(epsilon)
  #except: 
	  #print "Interpolation for "+str(wavelength*1E9)+" nm failed."

  # defining the new n and k on a common mesh
  ni = fni(wavelengths)
  ki = fki(wavelengths)

  plt.figure()
  plt.xlabel(r'$\mathcal{R}e(\varepsilon)$ (nm)')
  plt.ylabel('n, k')
  plt.semilogx(wavelength1*1e9, n, 'bs', label='n Palik')
  plt.semilogx(wavelength2*1e9, k, 'rs', label='k Palik')
  plt.semilogx(wavelengths*1e9, ni, 'b-', label='n interp')
  plt.semilogx(wavelengths*1e9, ki, 'r-', label='k interp')
  plt.grid()
  plt.legend(loc=1)
  plt.savefig('PalikData.eps')
  return epsilon
  #plt.show()
#}}}

## Imports (wavelength, ReEps, ImEps) data captured using a software like Engauge Digitizer. 
#   Input: (epsReal, epsImag) are discretized, also works on DIFFERENT MESHES. 
#   Output: (n,k, epsilon) discretized on the same mesh. 
#   
#   NOTE: Warning: the data produced by this method are rather unprecise. SPP spectoscopy requires precision to 1E-3. 
#   This method gives a precision worst then 1E0. Use only in case no other data are available. 
def importFromEpsilonTable(wavelength, folder, filename, plotting=True): #{{{
  nfile = folder+filename+"-epsR.csv" #TODO: rename nfile to ReEpsFile
  kfile = folder+filename+"-epsC.csv" #TODO: rename kfile to ImEpsFile
  
  print "** Info: opening "+nfile+" and "+kfile+"."
  narray = loadtxt(nfile, delimiter="\t", skiprows=1)
  karray = loadtxt(kfile, delimiter="\t", skiprows=1)
  
  unit = 1E-10 #Unit of the wavelength found in databases <nfile> and <kfile>. 
  
  # import wavelength, epsilonR, epsilonC
  wavelength1 = narray[:,0]*unit
  wavelength2 = karray[:,0]*unit
  n = narray[:,1]
  kk = karray[:,1]
  numrows = 10000 #Number of rows to interpolate the data on. 
  base = 10
  # interpolate on new wavelength mesh using 1st order
  order=1
  #wavelengths = np.arange(np.amin(wavelength2),np.amax(wavelength2), precision) #regular mesh, AWFUL for memory
  print "** Info: Generating the new wavelength mesh: ("+str(np.amin(wavelength2))+", "+str(np.amax(wavelength2))+")"
  wavelengths = np.logspace(np.amin(np.log10(wavelength2)), np.amax(np.log10(wavelength2)), num=numrows, base=base, endpoint = True)

  print "New wavelength mesh has "+str(numrows)+" rows."
  #print wavelengths
  
  print "** Info: definition of interpolation functions..."
  wavelength1=np.sort(wavelength1)
  wavelength2=np.sort(wavelength2)
    
  fni = InterpolatedUnivariateSpline(wavelength1, n, k=order)
  fki = InterpolatedUnivariateSpline(wavelength2, kk, k=order)
  
  #Interpolation for one optical constant
  try:
    ni = fni(wavelength); ki = fki(wavelength)
    #print wavelength, ni, ki
    epsilon = (ni+1j*ki) #NOTE: we are picking up the epsRe, and epsIm directly here
    print "" 
    print "Interpolated permittivity at "+str(wavelength*1E9)+" nm = "+str(epsilon)
  except: 
    print "Interpolation for "+str(wavelength*1E9)+" nm failed."
	  
  # defining the new epsR and epsC on a common mesh
  ni = fni(wavelengths)
  ki = fki(wavelengths)  

  if(plotting):
    plt.figure()
    plt.xlabel(r'$\mathcal{R}e(\varepsilon)$ (nm)')
    plt.ylabel('n, k')
    plt.semilogx(wavelength1*1e9, n, 'bs', label='Re(eps) data')
    plt.semilogx(wavelength2*1e9, kk, 'rs', label='Im(eps) data')
    plt.semilogx(wavelengths*1e9, ni, 'b-', label='Re(eps) interp')
    plt.semilogx(wavelengths*1e9, ki, 'r-', label='Im(eps) interp')
    plt.grid()
    plt.legend(loc=1)
    plt.savefig('GraphData.eps')
    plt.show()

  return epsilon
  
#}}}

## Mere function? Generates (ReEps, ImEps) from absorption data (given in m^{-1}). 
#  This routine is made to import data captured using sofware such as Engauge Digitized. 
#  Be very careful! The data produced by this method are very unprecise. SPP spectroscopy requires precision to 1E-3. 
#  This method gives a precision worst then 1E0. Then, it is only in case we have no other data. 
def importFromAbsorptionData(wavelength, folder, filename, plotting): #{{{
  absfile = folder+filename+".csv"
  narray = loadtxt(absfile, delimiter="\t", skiprows=1)
  unit = 1E9
  
  # import wavelength, alpha from spectroscopic data
  wavelength1 = h*c/(narray[::-1,0]*e)
  print wavelength1
  n = 1e2*narray[::-1,1]
  print n
  numrows = 10000
  base = 10
  # interpolate n and k on new wavelength mesh
  order=1
  #wavelengths = np.arange(np.amin(wavelength2),np.amax(wavelength2), precision) #regular mesh, AWFUL for memory
  print "Generating new wavelength mesh: ("+str(np.amin(wavelength1))+", "+str(np.amax(wavelength1))+")"
  wavelengths = np.logspace(np.amin(np.log10(wavelength1)), np.amax(np.log10(wavelength1)), num=numrows, base=base, endpoint = True)

  print "New wavelength mesh has "+str(numrows)+" rows."
  #print wavelengths
  
  fni = InterpolatedUnivariateSpline(wavelength1, n, k=order)

  #Interpolated one optical constants
  try:
    absnew = fni(wavelength)
    print "Interpolated absorptivity at "+str(wavelength*unit)+" nm = "+str(absnew)
  except:
    print "Interpolation for "+str(wavelength*unit)+" nm failed."

  # defining the new n and k on a common mesh
  ni = fni(wavelengths)

  if(plotting):
    plt.figure()
    plt.xlabel(r'Wavelength $\lambda$ (nm)')
    plt.ylabel(r'$\alpha$ (m$^{-1}$)')
    plt.loglog(wavelength1*unit, n, 'bs', label=r'$\alpha$ data')
    plt.loglog(wavelengths*unit, ni, 'b-', label=r'$\alpha$ interp')
    plt.grid()
    plt.legend(loc=1)
    plt.show()
    plt.savefig('SpectroscopicData.eps')
    
  return 0
  #plt.show()
#}}}

## Simply plots the optical data taken from a compatible database, 
# and interpolate palik data from tables of Palik at the given wavelength. 
# Useful to add one set of (ReEps, ImEps) for ONE wavelength in MaterialOpticalDatabaseForPlasmonics.csv. 
# INPUT
# @param wavelength: (float) value of desired output wavelength
# @param folder:     (str) name of the folder were database can be found
# @param filename:   (str) name of the material file
# @param plotting:   (boolean) plot the full data if True
def importFromTable(wavelength, folder, filename, plotting): #{{{
  # Look for Palik into the name
  if(filename.find("Palik") > 0):
    unit1 = 1E-10 #Palik data
  else:
    unit1 = 1E-6 #Other data
 
  # Fetch data
  DataFile = folder+filename
  try:
    DataArray = loadtxt(DataFile, delimiter="\t", skiprows=4)
    wavelengths = DataArray[:,0]*unit1; 
  except:
    DataArray = loadtxt(DataFile, delimiter=" ", skiprows=4)
    wavelengths = DataArray[:,0]*unit1
    
  n = DataArray[:,1]; kk = DataArray[:,2];

  # Interpolating using splines
  order = 1
  fni = InterpolatedUnivariateSpline(wavelengths, n, k=order)
  fki = InterpolatedUnivariateSpline(wavelengths, kk, k=order)

  #Interpolated one optical constants
  #wavelength = 1030e-9
  ni = fni(wavelength); ki = fki(wavelength)
  epsilon = (ni+1j*ki)**2
  print "Interpolated permittivity at "+str(wavelength*1e9)+" nm = "+str(epsilon)

  # Interpolate the full array and check it visually
  nimesh = fni(wavelengths); kimesh = fki(wavelengths)

  if (plotting): 
	plt.figure()
	plt.xlabel(r'$\lambda$ (nm)')
	plt.ylabel('n, k')
	plt.semilogx(1e9*wavelengths, n, 'bs', label='n Palik')
	plt.semilogx(1e9*wavelengths, kk, 'rs', label='k Palik')
	plt.semilogx(1e9*wavelengths, nimesh, 'b-', label='n interp')
	plt.semilogx(1e9*wavelengths, kimesh, 'r-', label='k interp')
	plt.grid()
	plt.legend(loc=2)
	plt.savefig('PalikData.eps')
	plt.show()
	
  return epsilon
#}}}

#==============================

if(len(sys.argv)<=2):
  print "Usage: ./importPalikData.py <Name of the material (Be, Au, ...)> <wavelength (nm)> <Source for data: Palik or name of the 1st author>"
  print "Example: ./importPalikData.py Au 800 Palik"
  exit()
 
material = sys.argv[1]
try:
  wavelength = float(sys.argv[2])
except: 
  print "** Error: command line must have the hape: <Material symbol> <Wavelength (nm)> <Publication author>"
  exit()
  
try:
  source = sys.argv[3]
except:
  print "** Warning: no data source given: using DEFAULT=Palik"
  source = "Palik"
  
filename = material+"-"+source
  
#filename = "Au-Palik"
#filename = "Be-Palik"
#filename = "Fe-Palik"
#filename = "Ni-Palik"
#filename = "a-Si-Palik"
#filename = "Si-Palik"
#filename = "SiC-Palik"
#filename = "Ti-Palik"
#filename = "Ag-Johnson"
#filename = "Ti-Johnson"
#filename = "SiO2-Palik"
#filename = "W-Palik"
#filename = "Cr-Johnson"
#filename = "BK7"
#filename = "Al-Palik"
plotting = True

#folder = "Database/PalikGraph/"
#filename = "c-Si-77K-Dash"
#filename = "c-Si-80K-Humlicek"
#filename = "Ti-Palik"
#filename="ZnO-Bond"

#=========================================

try: 
  if(source == "Palik"):
    folder = "Database/"
    print "Material: "+filename+"."
    print "Wavelength = "+str(wavelength)+" nm"
    epsilon = importFromTable(wavelength*1e-9, folder, filename, plotting)
    print epsilon
    print "You can add the following directly inside 'MaterialOpticalDatabaseForPlasmonics.csv'"
    print filename+"\t"+"?"+"\t"+str(int(wavelength))+"\t"+str(epsilon.real)+"\t"+str(epsilon.imag)+"\t?\t?\t?\t?\t?"
  elif(source == "Miller"):
    folder = "Database/LiquidMetals/"
    filename = "l-Cu-Miller"
    epsilon = importFromEpsilonTable(wavelength*1E-9, folder, filename, plotting)
    print epsilon
    print "You can add the following directly inside 'MaterialOpticalDatabaseForPlasmonics.csv'"
    print filename+"\t"+"?"+"\t"+str(int(wavelength))+"\t"+str(epsilon.real)+"\t"+str(epsilon.imag)+"\t?\t?\t?\t?\t?"
  else:
    folder = "Database/PalikGraph/"
    print "Material: "+filename+"."
    print "Wavelength = "+str(wavelength)+" nm"
    epsilon=importFromNKtable(wavelength*1e-9, folder, filename)
    print epsilon
    print "You can add the following directly inside 'MaterialOpticalDatabaseForPlasmonics.csv'"
    print filename+"\t"+"?"+"\t"+str(int(wavelength))+"\t"+str(epsilon.real)+"\t"+str(epsilon.imag)+"\t?\t?\t?\t?\t?"
except:
  print "Failed to import "+folder+filename+"!"
  print "Goto Database/importPalikData.sh for finding other sources."
  
#==========================================
#print "Lambda = 3000 nm"
#importFromTable(3000e-9, folder, filename, plotting=True)
#print "Lambda = 1064 nm"
#importFromTable(1064e-9, folder, filename, plotting=True)
#print "Lambda = 1060 nm"
#importFromTable(1060e-9, folder, filename, plotting)
#print "Lambda = 1030 nm"
#importFromTable(1030e-9, folder, filename, plotting)
#print "Lambda = 800 nm"
#importFromTable(800e-9, folder, filename, plotting)
#print "Lambda = 795 nm"
#importFromTable(795e-9, folder, filename, plotting)
#plotting = False
#print "Lambda = 625 nm"
#importFromTable(625e-9, folder, filename, plotting)
#print "Lambda = 532 nm"
#importFromTable(532e-9, folder, filename, plotting)
#print "Lambda = 515 nm"
#importFromTable(515e-9, folder, filename, plotting)
#print "Lambda = 400 nm"
#importFromTable(400e-9, folder, filename, plotting)

#===========================================
#importFromNKtable(folder, filename)
#print "Lambda = 3000 nm"
#importFromNKtable(3000e-9, folder, filename)

#importFromEpsilonTable(wavelength, folder, filename, plotting=True)
#===========================================
#print "Lambda = 515 nm"
#importFromAbsorptionData(515e-9, folder, filename, plotting=True)
