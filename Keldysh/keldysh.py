#!/usr/bin/env python
#-*- coding: utf-8 -*-
## @package Keldysh
## Computes the Keldysh excitation rate of quasi-free electrons
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

ShortRefKeldysh = "[Keldysh (1965)]"
ShortRefGruzdev = "[Gruzdev (2014)]"

print "Defining material parameters..."
Egap = 2.58e0*e; #LDA band gap of Si: 2.58 eV. #1.12e0*e for indirect band gap; 
#EgapEff = Egap; #starting before iterating to check convergence
meff=0.18e0; 

print "Meshing grids for intuitive calculations..."
Efield = 1E8*np.arange(1,100,0.25) #array for simple calculations

print "Defining the laser pulse..."
wavelength = 800e-9
tau=10e-15
PeakFluence = 100e-3*1E4 #J/cm2 * 1E4 = J/m2
PeakIntensity = PeakFluence/tau #scalar

tmin = -3.5*tau
tmax = 3.5*tau
dt = 1E-17
print "** Info: Number of time steps = "+str(int((tmax-tmin)/dt))+"."
instants=np.arange(tmin,tmax,dt)

PulseEnvelope=PulseGaussianTemporalShape(instants, tau, PeakIntensity)
print "** Info: Peak intensity = "+str(PulseEnvelope.max())+" W/m^2."
print "** Info: Peak field amplitude = "+str(IntensityToField(PulseEnvelope).max()/1E9)+" V/nm."

#print "** Starting the self-consistent loop..."

#for i in np.arange(1,50,1): #attempt of self consistent loop: divergent
#print "** ITERATION "+str(i)
print "Computing Adiabadicity coefficients for the pulse envelope..."

gamma = gammaKeldysh(Egap, meff, IntensityToField(PulseEnvelope), wavelength) #valid for scalar data
#gamma = gammaKeldysh(EgapEff, meff, IntensityToField(PulseEnvelope), wavelength) #self-consistent, divergent
print "** Info: Adiabadicity parameter = "+str(gamma.min())+"."

#print "Computing Keldysh1, Keldysh2..."
k1 = Keldysh1(gamma); k2 = Keldysh2(gamma) #valid
EgapEff = EffectiveGap(Egap, k1, k2) # Original formula from Keldysh. Warning: scipy.special.ellipe (https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.ellipe.html#scipy.special.ellipe) uses a different convention than Maple, Wikipedia or mpmath.
#EgapEff = EffectiveGap(EgapEff, k1, k2) #This formula was made self-consistent, but divergent. 
order = 50

print ""
print "** Info: Egap = "+str(Egap/e)+" eV, max[Ueff] = "+str(EgapEff.max()/e)+" eV."
print ""

KeldyshFunctionResult = KeldyshFunction( k1, k2, EgapEff, order, wavelength )
KeldyshFunctionResultG = KeldyshFunction_Gruzdev( k1, k2, EgapEff, order, wavelength )

#print "** End of self-consistent loop..."
print "Computes w_PI (Keldysh), and w_PIg (Gruzdev) for the pulse envelope..."
wPI = IonizationRate(k1, k2, KeldyshFunctionResult, EgapEff, wavelength)
wPIg= IonizationRate_Gruzdev(k1, k2, KeldyshFunctionResultG, EgapEff, wavelength)
print "w_PI until order "+str(order)+" = ", wPI.max()

# Current Keldysh model considers integration over time where w_PI is changing with instantaneous intensity. 

print ""
print "Temporal integration..."
# Just multiply array of w_PI by dt, and calculate cumsum().

dN_excited_Keldysh = np.multiply(wPI, dt)
dN_excited_Gruzdev = np.multiply(wPIg, dt)

# Temporal integration
N_excited_Keldysh = dN_excited_Keldysh.cumsum()
N_excited_Gruzdev = dN_excited_Gruzdev.cumsum()

print "Plotting..."
print ""
print "** Warning: results may be not converged."
print "            Reduce dt, and increase order until convergence."
print ""
print "Maximum density N_ex "+ShortRefKeldysh+" = "+str(N_excited_Keldysh.max())+"."
print "Maximum density N_ex "+ShortRefGruzdev+" = "+str(N_excited_Gruzdev.max())+"."

xunit = 1E15
timeunit = "fs"

plt.figure()

plt.subplot(311)
plt.title(r"Gap = "+str(Egap/e)+" eV, $\lambda=$ "+str(wavelength*1E9)+r" nm, $\tau=$"+str(tau*xunit)+" "+timeunit+", "+r"$F_{max}=$"+str(PeakFluence/1E4)+" J/cm"+r"$^{2}$")
#plt.xlabel("Field (V/m)")
#plt.xlabel("Time (ps)")
plt.ylabel("Adiabadicity $\gamma$")
plt.semilogy(instants*xunit, gamma, color="k", linestyle="-", label=r"$\gamma$")
plt.grid()
# plt.loglog(Efield, 0.1, label="Tunnelling limit")
# plt.loglog(Efield, 10.*np.ones(), label="MPI limit")
plt.legend(loc=3)

plt.subplot(312)
#plt.xlabel("Field (V/m)")
#plt.xlabel("Time (ps)")
plt.ylabel("$w_{PI}$ (m$^{-3}$ s$^{-1}$)")
plt.plot(instants*xunit, wPI, linestyle="-", color="r", label=r"$w_{PI}$ "+ShortRefKeldysh)
plt.plot(instants*xunit, wPIg, linestyle="-", color="b", label=r"$w_{PI}$ "+ShortRefGruzdev)
plt.grid()
plt.legend(loc=2)

plt.subplot(313)
plt.xlabel("Intensity (W/m$^{2}$)")
plt.xlabel("Time ("+timeunit+")")
# plt.xlabel("Field (V/m)")
plt.ylabel("Density (m$^{-3}$)")
plt.plot(instants*xunit, N_excited_Keldysh, color="r", label="$n_e$ "+ShortRefKeldysh)
plt.plot(instants*xunit, N_excited_Gruzdev, color="b", label="$n_e$ "+ShortRefGruzdev)
plt.grid()
plt.legend(loc=2)
plt.savefig("KeldyshAnalytic.eps")
plt.show()

#TODO: print "Exporting density to a table for ZnO optical index calculations."