#!/usr/bin/env python
#-*- coding: utf-8 -*-
## @package Keldysh
# This module aims to calculate the density of excited electrons as function of laser parameters. 
# Two types of usage are planned :
# * Generating tables to use directly into simulation codes
# * Outputing density in certain conditions. 
from libKeldysh import *

print ""
print "** Welcome to SPP-extended-theory suite."
print "Author(s): T.J.-Y. Derrien"
print ""
print "** Loading Keldysh module [Keldysh, Sov. J. Exp. Th. Phys. 47, 5 (1964)]..."
print "** Loading Gruzdev formula [Gruzdev, Optical Engineering 53, 122515 (2014)]"

print "Defining material parameters..."
Egap = 2.58e0*e; #LDA band gap of Si: 2.58 eV. #1.12e0*e for indirect band gap; 
EgapEff = Egap; #starting before iterating to check convergence
meff=0.18e0; 

print "Meshing grids for intuitive calculations..."
Efield = 1E8*np.arange(1,100,0.25) #array for simple calculations

print "Defining the laser pulse..."
wavelength = 800e-9
tau=100e-15
PeakFluence = 100e-3*1E4 #J/cm2 * 1E4 = J/m2
PeakIntensity = PeakFluence/tau #scalar

tmin = -3.5*tau
tmax = 3.5*tau
dt = 1E-14
print "** Info: Number of time steps = "+str(int((tmax-tmin)/dt))+"."
instants=np.arange(tmin,tmax,dt)

PulseEnvelope=PulseGaussianTemporalShape(instants, tau, PeakIntensity)
print "** Info: Peak intensity = "+str(PulseEnvelope.max())+" W/m^2."
print "** Info: Peak field amplitude = "+str(IntensityToField(PulseEnvelope).max()/1E9)+" V/nm."

#print "** Starting the self-consistent loop..."

#for i in np.arange(1,50,1): #attempt of self consistent loop: divergent
#print "** ITERATION "+str(i)
print "Computing Adiabadicity coefficients..."

gamma = gammaKeldysh(Egap, meff, IntensityToField(PulseEnvelope), wavelength) #valid for scalar data
#gamma = gammaKeldysh(EgapEff, meff, IntensityToField(PulseEnvelope), wavelength) #self-consistent, divergent
print "** Info: Adiabadicity parameter = "+str(gamma.min())+"."

print "Computing Keldysh1, Keldysh2..."
k1 = Keldysh1(gamma); k2 = Keldysh2(gamma) #valid
EgapEff = EffectiveGap(Egap, k1, k2) # Original formula from Keldysh. Warning: scipy.special.ellipe (https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.ellipe.html#scipy.special.ellipe) uses a different convention than Maple, Wikipedia or mpmath.
#EgapEff = EffectiveGap(EgapEff, k1, k2) #This formula was made self-consistent, but divergent. 
order = 10

print ""
print "** Info: Egap = "+str(Egap/e)+" eV, Ueff = "+str(EgapEff.max()/e)+" eV."
print ""

KeldyshFunctionResult = KeldyshFunction( k1, k2, EgapEff, order, wavelength )
KeldyshFunctionResultG = KeldyshFunction_Gruzdev( k1, k2, EgapEff, order, wavelength )

#print "** End of self-consistent loop..."

wPI = IonizationRate(k1, k2, KeldyshFunctionResult, EgapEff, wavelength)
wPIg= IonizationRate_Gruzdev(k1, k2, KeldyshFunctionResultG, EgapEff, wavelength)
print "w_PI until order "+str(order)+" = ", wPI

# Current Keldysh model considers integration over time where w_PI is changing with instantaneous intensity. 
exit()

print ""
print "Temporal integration..."


print "Plotting..."
plt.figure()
plt.subplot(211)
plt.xlabel("Field (V/m)")
plt.ylabel("Excitation rate $w_{PI}$ (m$^{-3}$ s$^{-1}$)")
plt.loglog(Efield, wPI, linestyle="-", color="r", label=r"$w_{PI}$ original")
plt.loglog(Efield, wPIg, linestyle="-", color="b", label=r"$w_{PI}$ corrected")
plt.legend(loc=2)

plt.subplot(212)
plt.xlabel("Field (V/m)")
plt.ylabel("Adiabadicity parameter")
plt.loglog(Efield, gamma, color="k", linestyle="-", label=r"$\gamma$")
# plt.loglog(Efield, 0.1, label="Tunnelling limit")
# plt.loglog(Efield, 10.*np.ones(), label="MPI limit")
plt.legend(loc=3)

#plt.subplot(313)
#plt.xlabel("Intensity (W/m$^{2}$)") #Field (V/m)")
## plt.xlabel("Field (V/m)")
#plt.ylabel("Instant density (m$^{-3}$)")
#plt.loglog(FieldToIntensity(Efield), wPI*tau, label="$n_e$ estim.")
#plt.loglog(FieldToIntensity(Efield), wPIg*tau, label="$n_e$ estim. corrected")
#plt.legend(loc=2)
#plt.savefig("KeldyshAnalytic.eps")
#plt.show()
