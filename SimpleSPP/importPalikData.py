#!/usr/bin/env python
#-*- coding: utf-8 -*-

# IMPORT LIBRARIES
from SimpleSPProutines import *

# Importing data from Palik book using graphs. 
# Optical data are given in csv files, created using Engauge-digitizer software. 

def importFromPalikGraph(folder, filename):
  """This routine is made to import data captured using Engauge Digitized. 
  Be very careful! The data produced by this method are very unprecise. SPP spectoscopy requires precision to 1E-3. 
  This method gives a precision worst then 1E0. Then, it is only in case we have no other data. 
  """
  nfile = folder+filename+"-n.csv"
  kfile = folder+filename+"-k.csv"

  narray = loadtxt(nfile, delimiter="\t", skiprows=1)
  karray = loadtxt(kfile, delimiter="\t", skiprows=1)

  # import wavelength, n and k from Palik
  wavelength1 = narray[:,0]*1e-6
  wavelength2 = karray[:,0]*1e-6
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
  wavelength = 800e-9
  ni = fni(wavelength); ki = fki(wavelength)
  epsilon = (ni+1j*ki)**2
  print "Interpolated permittivity at "+str(wavelength*1e9)+" nm = "+str(epsilon)

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
  return 0
  #plt.show()

def importFromTables(wavelength, folder, filename):
  # interpolate palik data from tables of Palik

  # Fetch data
  DataFile = folder+filename
  DataArray = loadtxt(DataFile, delimiter="\t", skiprows=4)
  wavelengths = DataArray[:,0]*1e-10; n = DataArray[:,1]; kk = DataArray[:,2];

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


folder = "Database/"
filename = "Pt-Palik"

#folder = "Database/PalikGraph/"
#filename = "Ti-Palik"
#print "Lambda = 1064 nm"
#print "Lambda = 1030 nm"
print "Lambda = 800 nm"
importFromTables(800e-9, folder, filename)
print "Lambda = 400 nm"
importFromTables(400e-9, folder, filename)
#importFromPalikGraph(folder, filename)