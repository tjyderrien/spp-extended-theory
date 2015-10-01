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

#MaterialFolder="/usr/local/share/gsvit/data/spectra"
MaterialFolder="Database"

MaterialFile1="Air"
MaterialFile2='Si-Aspnes'

# Loading Material dielectric complex permittivity into arrays
MaterialArray2 = loadtxt(MaterialFolder+'/'+MaterialFile2, delimiter=' ', skiprows=4)
#wavelengths2 = MaterialArray2[:,0]
nlines, ncols = MaterialArray2.shape

#print nlines, ncols

try:
  MaterialArray1 = loadtxt(MaterialFolder+'/'+MaterialFile1, delimiter=' ', skiprows=4)
except:
  print "Material 1 ("+MaterialFile1+") was not found in "+MaterialFolder+"."
  print "Material 1 was replaced by Air."
  
  MaterialArray1 = np.zeros((nlines, 3))
  MaterialArray1[:,1] = np.ones(nlines) #Air index is 1. 
  MaterialArray1[:,0] = MaterialArray2[:,0] #Copy the table of wavelengths

#print "Shape of MaterialArray1 is "+str(MaterialArray1.shape)
#print "Shape of MaterialArray2 is "+str(MaterialArray2.shape)

# Let's build eps1 and eps2
wavelengths1=MaterialArray1[:,0]
wavelengths1=np.multiply(wavelengths1,1e-6) #converting wavelength to emters
n=MaterialArray1[:,1]; k=MaterialArray1[:,2]
eps1=np.power(np.add(n,np.multiply(1j, k)),2) #conversion to epsilon
del n, k
#print wavelengths1

wavelengths2=MaterialArray2[:,0]
wavelengths2=np.multiply(wavelengths2,1e-6) #converting wavelength to meters
n=MaterialArray2[:,1]; k=MaterialArray2[:,2]
eps2=np.power(np.add(n,np.multiply(1j, k)),2)
del n, k
#print eps2

# Before doing calculations, we shall interpolate the most dense mesh on the second, and take their intersection. 
print "Wavelength mesh size 1="+str(eps1.size)
print "Wavelength mesh size 2="+str(eps2.size)

# dense mesh generation
wavelengths = np.arange(np.amin(wavelengths2),np.amax(wavelengths2),1e-9)
order=3
feps1r=InterpolatedUnivariateSpline(wavelengths1, eps1.real, k=order)
feps1i=InterpolatedUnivariateSpline(wavelengths1, eps1.imag, k=order)
feps2r=InterpolatedUnivariateSpline(wavelengths2, eps2.real, k=order)
feps2i=InterpolatedUnivariateSpline(wavelengths2, eps2.imag, k=order)

print "Interpolating on Wavelength mesh size = "+str(wavelengths.size)

eps1new=np.add(feps1r(wavelengths),np.multiply(1.0j, feps1i(wavelengths)))
eps2new=np.add(feps2r(wavelengths),np.multiply(1.0j, feps2i(wavelengths)))

### Checking interpolation of the dielectric function
plt.figure()
plt.xlabel('Wavelength (nm)')
plt.ylabel('epsilon')
#plt.plot(1e9*wavelengths1, eps1.real, '-', label='Re('+MaterialFile1+')')
plt.plot(1e9*wavelengths, eps1new.real, '-', label='interp $Re('+MaterialFile1+')$')
#plt.plot(1e9*wavelengths1, eps1.imag, '-', label='Im('+MaterialFile1+')')
plt.plot(1e9*wavelengths, eps1new.imag, '--', label='interp $Im('+MaterialFile1+')$')
#plt.plot(1e9*wavelengths2, eps2.real, '-', label='Re('+MaterialFile2+')')
plt.plot(1e9*wavelengths, eps2new.real, '-', label='interp $Re('+MaterialFile2+')$')
#plt.plot(1e9*wavelengths2, eps2.imag, '-', label='Im('+MaterialFile2+')')
plt.plot(1e9*wavelengths, eps2new.imag, '--', label='interp $Im('+MaterialFile2+')$')
plt.legend(loc=2)
plt.title('Dielectric permittivity')
plt.savefig('epsilon.png')

## plot SPP dispersion relation

betaSPP = np.vectorize(betaSPP)
omega = np.vectorize(omega)

kspp = betaSPP(wavelengths, eps1new, eps2new)
omegaspp = omega(wavelengths)
#omegaspp = np.sort(omegaspp)

plt.figure()
plt.xlabel('k $(m^{-1})$')
plt.ylabel('$\omega$ ($s^{-1}$)')
plt.plot(kspp.real, omegaspp, label='SPP')
plt.plot(omega(wavelengths)/c, omega(wavelengths), label='Light line')
plt.plot(omega(wavelengths)/c, omega(np.add(np.multiply(wavelengths,0e0), 800e-9)), label='Laser 800 nm')
plt.title('Dispersion relation at $'+MaterialFile1+'$/$'+MaterialFile2+'$ interface')
plt.legend(loc=2)
plt.savefig('Dispersion.png')

## plot the lifetime with wavelength
#RealDerivativeByComplex = np.vectorize(RealDerivativeByComplex)
SPPgroupVelocity = RealDerivativeByComplex(omegaspp, kspp)
SPPphaseVelocity = np.divide(omegaspp,kspp)

print np.shape(omegaspp)
print np.shape(kspp)
print np.shape(SPPgroupVelocity)

#SPPgroupVelocityRe = SPPgroupVelocity.real
SPPgroupVelocityPlot = np.clip(SPPgroupVelocity.real, 0, 1000E8)

print SPPgroupVelocity.real

plt.figure()
plt.xlabel('Wavelength $(nm)$')
plt.ylabel('Velocity ($m/s$)')
plt.plot(1e9*2*pi*c/omegaspp[1:], SPPgroupVelocityPlot, label='$v_g$')
plt.plot(1e9*2*pi*c/omegaspp, SPPphaseVelocity, label='$v_{\phi}$')
#plt.plot(omega(wavelengths)/c, omega(wavelengths), label='Light line')
#plt.plot(omega(wavelengths)/c, omega(np.add(np.multiply(wavelengths,0e0), 800e-9)), label='Laser 800 nm')
plt.title('SPP velocities at $'+MaterialFile1+'$/$'+MaterialFile2+'$ interface')
plt.legend(loc=1)
#plt.axis([0,1000,0,4e8])
plt.savefig('Velocities.png')

## Now we can calculate SPP lifetime
LifeTimeOld = LifeTimeRaether(kspp[1:], eps1new[1:], eps2new[1:])
LifeTimeNew = LifeTimeDerrien(kspp[1:], SPPgroupVelocity.real)

plt.figure()
plt.xlabel('Wavelength $(nm)$')
plt.ylabel('Lifetime ($s$)')
plt.semilogy(1e9*2*pi*c/omegaspp[1:], LifeTimeOld, label='Raether')
plt.semilogy(1e9*2*pi*c/omegaspp[1:], LifeTimeNew, label='This work')
plt.title('SPP velocities at $'+MaterialFile1+'$/$'+MaterialFile2+'$ interface')
plt.legend(loc=4)
plt.savefig('Lifetime.png')
