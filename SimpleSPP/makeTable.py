#!/usr/bin/env python
#-*- coding: utf-8 -*-

# IMPORT LIBRARIES
from SimpleSPProutines import *

# PHYSICAL INPUT
wavelength = 800e-9
epsAir=1e0
epsSi0=13.64+0.048j
meffe=0.18
nuSi=(1.1e-15)**-1

#example=period(betaSPP(wavelength, epsAir, Drude(wavelength, 1e28, epsSi0, nuSi)))

#print example

"""
We would like now to construct a database using available materials description with all possible interfaces
we will : 

1. For each material in database, select each material and verify, for each available wavelength, 
1.1: If SPP condition is verified, 
1.2. yes, then period can be calculated and shown;
1.4. SPP decay depth in medium 1
1.5. SPP decay depth in medium 2
2. Then extract a table which contains all possible scenarios
"""
# Select database
database="MaterialOpticalDatabaseForPlasmonics.csv"

# Build database array for choosing which material can be of interest to irradiate
dbarray = loadtxt(database, dtype='str', delimiter='\t')
#SPPactiveInterfaces(dbarray, 'new')

""" TODO: interface this with HTML for publication on the web. 
1. Put results into a NP.array.
2. Use a converter to HTML, CSV and PDF maybe. 
"""

# Now, we shall construct database for SPP lifetimes. Actually, SPP lifetime require the knowledge of all spectrum of response to be known. 
MaterialFolder="/usr/local/share/gsvit/data/spectra"
MaterialFile1="SiO2"
MaterialFile2="Ag"

# Loading Material dielectric complex permittivity into arrays
MaterialArray2 = loadtxt(MaterialFolder+'/'+MaterialFile2, delimiter=' ', skiprows=4)
wavelengths = MaterialArray2[:,0]
nx = MaterialArray2.size
try:
  MaterialArray1 = loadtxt(MaterialFolder+'/'+MaterialFile1, delimiter=' ', skiprows=4)
except:
  print "Material 1 ("+MaterialFile1+") was not found in "+MaterialFolder+"."
  #TODO: manage the exception
  #print "Material 1 was replaced by Air."
  
  #print np.array([wavelengths, np.ones(nx), np.zeros(nx)])
  #MaterialArray1 = np.zeros(nx, 3)
  #print MaterialArray1 = np.ones(nx)
  #print np.array(MaterialArray1, ndmin=2)
  #print MaterialArray1

# Let's build eps1 and eps2
wavelengths1=MaterialArray1[:,0]
wavelengths1=np.multiply(wavelengths1,1e-6) #converting wavelength to emters
n=MaterialArray1[:,1]; k=MaterialArray1[:,2]
eps1=np.power(np.add(n,np.multiply(1j, k)),2)
#print eps1

wavelengths2=MaterialArray2[:,0]
wavelengths2=np.multiply(wavelengths2,1e-6) #converting wavelength to emters
n=MaterialArray2[:,1]; k=MaterialArray2[:,2]
eps2=np.power(np.add(n,np.multiply(1j, k)),2)
#print eps2

# Before doing calculations, we shall interpolate the most dense mesh on the second, and take their intersection. 
print eps1.size, eps2.size

xnew = np.arange(100e-9,2e-6,10e-9)
order=1
feps1=InterpolatedUnivariateSpline(wavelengths1, eps1, k=order)
feps2=InterpolatedUnivariateSpline(wavelengths2, eps2, k=order)
eps1new=feps1(xnew)
eps2new=feps2(xnew)
# TODO: careful this does not work for complex-valued functions

def RealDerivativeByComplex(f,z):
  """Complex derivative a real-valued function f: z->f(z)
  """
  return 0.5*(np.diff(f,z.real) - 1j*np.diff(z.imag))

#def NewLifeTime():
  #"""
  #Input:
    #eps1:
    #eps2:
    #wavelength:
  #Output:
    #Lifetime array with wavelength
  #"""
  #beta=betaSPP(wavelengths, eps1, eps2)
  #SPPgroupVelocity=RealDerivativeByComplex(omega(wavelengths), beta)
  #SPPgroupVelocity=SPPgroupVelocity.real
  #LifeTime=0.5/betaSPP.real * (SPPgroupVelocity)**(-1)
  #return LifeTime

## Show dielectric function
plt.figure()
plt.xlabel('Wavelength (nm)')
plt.ylabel('\varepsilon')
plt.plot(1e9*wavelengths1, eps1.real, '-', label='Re('+MaterialFile1+')')
plt.plot(1e9*wavelengths1, eps1.imag, '--', label='Im('+MaterialFile1+')')
plt.plot(1e9*wavelengths2, eps2.real, '-', label='Re('+MaterialFile2+')')
plt.plot(1e9*wavelengths2, eps2.imag, '--', label='Im('+MaterialFile2+')')
plt.savefig('epsilon.png')
plt.show()

#print MaterialArray
#print LifeTimeSpectrum(Material)