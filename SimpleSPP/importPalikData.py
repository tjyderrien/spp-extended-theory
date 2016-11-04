#!/usr/bin/env python
#-*- coding: utf-8 -*-

# IMPORT LIBRARIES
from SimpleSPProutines import *

# Importing data from Palik book using graphs. 
# Optical data are given in csv files, created using Engauge-digitizer software. 
# OUTPUTS: 

def importFromNKtable(wavelength, folder, filename, plotting=1):#{{{
  """This routine is made to import data captured using a software like Engauge Digitizer. 
  This leads to obtain (n,k) discretized on DIFFERENT MESHES. 
  
  NOTE: Warning: the data produced by this method are rather unprecise. SPP spectoscopy requires precision to 1E-3. 
  This method gives a precision worst then 1E0. Use only in case no other data are available. 
  """
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

def importFromEpsilonTable(wavelength, folder, filename, plotting): #{{{
  """This routine is made to import data captured using a software like Engauge Digitizer. 
  Input: (epsReal, epsImag) discretized on DIFFERENT MESHES. 
  Output: (n,k, epsilon) discretized on the same mesh. 
  
  NOTE: Warning: the data produced by this method are rather unprecise. SPP spectoscopy requires precision to 1E-3. 
  This method gives a precision worst then 1E0. Use only in case no other data are available. 
  """
  nfile = folder+filename+"-epsR.csv"
  kfile = folder+filename+"-epsC.csv"

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
	ni = fni(wavelength); ki = fki(wavelength)
	epsilon = (ni+1j*ki) #we are picking up the epsRe, and epsIm directly here
	print "Interpolated permittivity at "+str(wavelength*1e9)+" nm = "+str(epsilon)
  except: 
	  print "Interpolation for "+str(wavelength*1E9)+" nm failed."
	  
  #try:
	#wavelength = 532e-9
	#ni = fni(wavelength); ki = fki(wavelength)
	#epsilon = (ni+1j*ki)
	#print "Interpolated permittivity at "+str(wavelength*1e9)+" nm = "+str(epsilon)
  #except: 
	  #print "Interpolation for "+str(wavelength*1E9)+" nm failed."
  #try:
	#wavelength = 400e-9
	#ni = fni(wavelength); ki = fki(wavelength)
	#epsilon = (ni+1j*ki)
	#print "Interpolated permittivity at "+str(wavelength*1e9)+" nm = "+str(epsilon)
  #except: 
	  #print "Interpolation for "+str(wavelength*1E9)+" nm failed."
  #try:
	#wavelength = 930e-9
	#ni = fni(wavelength); ki = fki(wavelength)
	#epsilon = (ni+1j*ki)
	#print "Interpolated permittivity at "+str(wavelength*1e9)+" nm = "+str(epsilon)
  #except: 
	  #print "Interpolation for "+str(wavelength*1E9)+" nm failed."

  # defining the new n and k on a common mesh
  ni = fni(wavelengths)
  ki = fki(wavelengths)

  if(plotting):
	  plt.figure()
	  plt.xlabel(r'$\mathcal{R}e(\varepsilon)$ (nm)')
	  plt.ylabel('n, k')
	  plt.semilogx(wavelength1*1e9, n, 'bs', label='n Palik')
	  plt.semilogx(wavelength2*1e9, k, 'rs', label='k Palik')
	  plt.semilogx(wavelengths*1e9, ni, 'b-', label='n interp')
	  plt.semilogx(wavelengths*1e9, ki, 'r-', label='k interp')
	  plt.grid()
	  plt.legend(loc=1)
	  plt.savefig('GraphData.eps')
  return 0
  #plt.show()
#}}}

def importFromAbsorptionData(wavelength, folder, filename, plotting): #{{{
  """This routine is made to import data captured using sofware such as Engauge Digitized. 
  Be very careful! The data produced by this method are very unprecise. SPP spectroscopy requires precision to 1E-3. 
  This method gives a precision worst then 1E0. Then, it is only in case we have no other data. 
  """
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

def importFromTable(wavelength, folder, filename, plotting): #{{{
  """ Description: Simply plots the optical data taken from a prepared database, 
  and interpolate palik data from tables of Palik at the given wavelength. 
  Useful to add one set of (ReEps, ImEps) for ONE wavelength in MaterialOpticalDatabaseForPlasmonics.csv. 
  INPUT
  - wavelength: (float) value of desired output wavelength
  - folder: (str) name of the folder were database can be found
  - filename: (str) name of the material file
  - plotting: (boolean) plot the full data if True
  """
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
wavelength = float(sys.argv[2])
try:
  source = sys.argv[3]
except:
  print "No data source given: using DEFAULT=Palik"
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
  else:
    folder = "Database/PalikGraph/"
    print "Material: "+filename+"."
    print "Wavelength = "+str(wavelength)+" nm"
    epsilon=importFromNKtable(wavelength*1e-9, folder, filename)
    print epsilon
    print "You can add the following directly inside 'MaterialOpticalDatabaseForPlasmonics.csv'"
    print filename+"\t"+"?"+"\t"+str(int(wavelength))+"\t"+str(epsilon.real)+"\t"+str(epsilon.imag)+"\t?\t?\t?\t?\t?"
except:
  print "Failed to import "+filename+"!"
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

#importFromEpsilonTable(515e-9, folder, filename, plotting=True)
#===========================================
#print "Lambda = 515 nm"
#importFromAbsorptionData(515e-9, folder, filename, plotting=True)
