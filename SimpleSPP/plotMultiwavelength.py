#!/usr/bin/env python
#-*- coding: utf-8 -*-

# IMPORT LIBRARIES
from SimpleSPProutines import *

# Now, we shall construct database for SPP lifetimes. Actually, SPP lifetime require the knowledge of all spectrum of response to be known. 

#MaterialFolder="/usr/local/share/gsvit/data/spectra"
MaterialFolder="Database"

MaterialFile1="Air"
#MaterialFile1="SiO2-Palik"
#MaterialFile1="Si-Palik"
#MaterialFile2="Ti-Palik"
#MaterialFile2="Ag-Palik"
# MaterialFile2="Au-Palik"
MaterialFile2="Au-Johnson"
#MaterialFile2="Ag-Johnson"
#MaterialFile2="Ti-Johnson"
#MaterialFile2="Mo-Palik"

UnitMat1=1e10
UnitMat2=1e6

# Loading Material dielectric complex permittivity into arrays
try: 
	MaterialArray2 = loadtxt(MaterialFolder+'/'+MaterialFile2, delimiter='\t', skiprows=4)
except: 
	print "Could not read "+MaterialFile2+" database."
	print "Attempting second method..."
	try:
		MaterialArray2 = loadtxt(MaterialFolder+'/'+MaterialFile2, delimiter=' ', skiprows=4)
		print "Success."
	except:
		print "Also failed reading of database... Exiting."
		exit()
		
#wavelengths2 = MaterialArray2[:,0]
nlines, ncols = MaterialArray2.shape

#print nlines, ncols

try:
  MaterialArray1 = loadtxt(MaterialFolder+'/'+MaterialFile1, delimiter='\t', skiprows=4)
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
wavelengths1=np.multiply(wavelengths1,UnitMat1**(-1)) #converting wavelength to emters
n=MaterialArray1[:,1]; k=MaterialArray1[:,2]
eps1=np.power(np.add(n,np.multiply(1j, k)),2) #conversion to epsilon
del n, k
#print wavelengths1

wavelengths2=MaterialArray2[:,0]
wavelengths2=np.multiply(wavelengths2,UnitMat2**(-1)) #converting wavelength to meters
n=MaterialArray2[:,1]; k=MaterialArray2[:,2]
eps2=np.power(np.add(n,np.multiply(1j, k)),2)
del n, k
#print eps2

# Before doing calculations, we shall interpolate the most dense mesh on the second, and take their intersection. 
print "Wavelength mesh size 1="+str(eps1.size)
print "Wavelength mesh size 2="+str(eps2.size)

# dense mesh generation
wavelengths = np.arange(np.amin(wavelengths2),np.amax(wavelengths2),1e-9)
order=1
feps1r=InterpolatedUnivariateSpline(wavelengths1, eps1.real, k=order)
feps1i=InterpolatedUnivariateSpline(wavelengths1, eps1.imag, k=order)
feps2r=InterpolatedUnivariateSpline(wavelengths2, eps2.real, k=order)
feps2i=InterpolatedUnivariateSpline(wavelengths2, eps2.imag, k=order)

print "Interpolating on Wavelength mesh size = "+str(wavelengths.size)

eps1new=np.add(feps1r(wavelengths),np.multiply(1.0j, feps1i(wavelengths)))
eps2new=np.add(feps2r(wavelengths),np.multiply(1.0j, feps2i(wavelengths)))

print "Checking quality of interpolation for the dielectric function..."
plt.figure()
plt.xlabel('Wavelength (nm)')
plt.ylabel('epsilon')
plt.semilogx(1e9*wavelengths1, eps1.real, '+', label='Re('+MaterialFile1+')')
plt.semilogx(1e9*wavelengths, eps1new.real, '-', label='interp $Re('+MaterialFile1+')$')
plt.semilogx(1e9*wavelengths1, eps1.imag, '+', label='Im('+MaterialFile1+')')
plt.semilogx(1e9*wavelengths, eps1new.imag, '-', label='interp $Im('+MaterialFile1+')$')
plt.semilogx(1e9*wavelengths2, eps2.real, '+', label='Re('+MaterialFile2+')')
plt.semilogx(1e9*wavelengths, eps2new.real, '-', label='interp $Re('+MaterialFile2+')$')
plt.semilogx(1e9*wavelengths2, eps2.imag, '+', label='Im('+MaterialFile2+')')
plt.semilogx(1e9*wavelengths, eps2new.imag, '-', label='interp $Im('+MaterialFile2+')$')
#plt.legend(loc=2)
plt.title('Dielectric permittivity')
plt.savefig('epsilon.png')
#plt.show()

print "Plot the SPP dispersion relation..."

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
plt.savefig('Dispersion.eps')
#plt.show()

print "Plot the SPP period with wavelength..."

plt.figure()
plt.xlabel('Wavelength $\lambda$ $(nm)$')
plt.ylabel('$\Lambda$ ($nm$)')
plt.plot(1e9*wavelengths, 1e9 * (2e0*pi/kspp.real), label='Near-field period')
plt.title('Period of field at $'+MaterialFile1+'$/$'+MaterialFile2+'$ interface')
plt.legend(loc=2)
plt.grid(True)
plt.savefig('Period.eps')
#plt.show()

print "Plot the lifetime with wavelength..."
#RealDerivativeByComplex = np.vectorize(RealDerivativeByComplex)
SPPgroupVelocity = RealDerivativeByComplex(omegaspp, kspp)
SPPphaseVelocity = np.divide(omegaspp,kspp)

print np.shape(omegaspp)
print np.shape(kspp)
print np.shape(SPPgroupVelocity)

#SPPgroupVelocityRe = SPPgroupVelocity.real
SPPgroupVelocityPlot = np.clip(SPPgroupVelocity.real, 0, 1000E8)

#print SPPgroupVelocity.real

plt.figure()
plt.xlabel(r'Wavelength $\lambda$ (nm)')
plt.ylabel(r'Velocity $v$ ($\mu$m/ps)')
plt.plot(1e9*2*pi*c/omegaspp[1:], 1E-6*SPPgroupVelocityPlot, 'b-', label=r'$v_g$')
plt.plot(1e9*2*pi*c/omegaspp, 1E-6*SPPphaseVelocity, 'r--', label=r'$v_{\phi}$')
#plt.plot(1e9*2*pi*c/omegaspp, c, label=r'$c$')
#plt.plot(omega(wavelengths)/c, omega(wavelengths), label='Light line')
#plt.plot(omega(wavelengths)/c, omega(np.add(np.multiply(wavelengths,0e0), 800e-9)), label=r'$c$')
#plt.title(MaterialFile1+'/'+MaterialFile2+' interface')
plt.legend(loc=4)
plt.xticks(np.arange(0, 2500, 500))
plt.axis([200,2000,0,300])
plt.savefig('Velocities.eps')
plt.savefig('Velocities.png')
#plt.show()

print "Plot SPP lifetime with wavelength..."
LifeTimeOld = LifeTimeRaether(kspp[1:], eps1new[1:], eps2new[1:])
LifeTimeNew = LifeTimeDerrien(kspp[1:], SPPgroupVelocity.real)
LifeTimePhase = LifeTimeDerrien(kspp[1:], SPPphaseVelocity[1:].real)
LifeTimeApprox = 2e0*(eps2new.real)**2 /( omegaspp * eps1new.real**2 * eps2new.imag)

print "Lifetime approx."
print LifeTimeApprox

# Let's clip Lifetime where they are negative. 
LifeTimeOld = np.clip(LifeTimeOld, 0, 1)
LifeTimeNew = np.clip(LifeTimeNew, 0, 1)
LifeTimePhase = np.clip(LifeTimePhase, 0, 1)
LifeTimeApprox = np.clip(LifeTimeApprox, 0, 1)

# Plotting the SPP lifetime

plt.figure()
plt.xlabel('Wavelength $\lambda$ (nm)')
plt.ylabel(r'SPP lifetime $\tau_{\mbox{SPP}}$ (ps)')
plt.plot(1e9*2*pi*c/omegaspp[1:], 1E12*LifeTimeOld, 'b-', label=r'Group, complex $\omega$, real $\beta$')
plt.plot(1e9*2*pi*c/omegaspp[1:], 1E12*LifeTimeNew, 'r-', label=r'Group, real $\omega$, complex $\beta$')
plt.plot(1e9*2*pi*c/omegaspp[1:], 1E12*LifeTimePhase, 'r--', label=r'Phase, real $\omega$, complex $\beta$')
plt.plot(1e9*2*pi*c/omegaspp, 1E12*LifeTimeApprox, 'b--', label=r'Phase, approx.')
#plt.title(MaterialFile1+'/'+MaterialFile2+' interface')
plt.legend(loc=2)
#plt.xticks(np.arange(0, 2500, 500))
plt.axis([400,1000,0,0.6])
plt.savefig('Lifetime.eps')
#plt.show()


# 
print "Plotting the SPP decay depth and optical penetration depth..."

SPPdecayDepth=DecayDepth(kzSPP(wavelengths, eps1new, eps2new))
SPPdecayDepth2=DecayDepth(kzSPP(wavelengths, eps2new, eps1new))
OpticalPenetrationDepth = 2e0*omega(wavelengths)/c*(eps2new**0.5e0)
OpticalPenetrationDepth = (OpticalPenetrationDepth.imag)**(-1e0)

plt.figure()
plt.xlabel('Wavelength $\lambda$ (nm)')
plt.ylabel(r'SPP decay depth $\delta_{\mbox{SPP}}$ ($\mu$lm)')
plt.semilogy(1e9*2*pi*c/omegaspp, 1E6*SPPdecayDepth, 'r-', label=r'$\delta_{SPP}$, medium 1')
plt.semilogy(1e9*2*pi*c/omegaspp, 1E6*SPPdecayDepth2, 'r--', label=r'$\delta_{SPP}$, medium 2')
plt.semilogy(1e9*2*pi*c/omegaspp, 1E6*OpticalPenetrationDepth, 'bx', label=r'$\delta_{OPD}$, medium 2')
plt.legend(loc=1)
#plt.axis([200,2000,0,10e0])
plt.savefig('SppDecayDepth.eps')
plt.savefig('SppDecayDepth.png')
plt.show()

