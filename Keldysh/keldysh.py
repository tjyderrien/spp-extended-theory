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
print "** TODO: Correct formula using [Gruzdev, Optical Engineering 53, 122515 (2014)"
print "** TODO: Change field amplitude to laser intensity."

Egap = 2.58e0*e; #LDA band gap of Si: 2.58 eV. #1.12e0*e; 
meff=0.18e0; 
Efield=1E9; 
wavelength = 800e-9
tau=100e-15


#}}}

print ""
print "** Generating mesh..."
Efield = 1E8*np.arange(1,100,0.25)

print ""
print "Computing gamma..."
gamma = gammaKeldysh(Egap, meff, Efield, wavelength) #valid
#print gamma

print "Computing Keldysh1, Keldysh2..."
k1 = Keldysh1(gamma); k2 = Keldysh2(gamma) #valid
EgapEff = EffectiveGap(Egap, k1, k2) #Warning: scipy.special.ellipe (https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.ellipe.html#scipy.special.ellipe) uses a different convention than Maple, Wikipedia or mpmath.
order = 10

print ""
print "** Info: "
print "Adiabadicity parameter = "+str(gamma)
print "Egap = "+str(Egap/e)+" eV, Ueff = "+str(EgapEff/e)+" eV."
print ""

KeldyshFunctionResult = KeldyshFunction( k1, k2, EgapEff, order, wavelength )
wPI = IonizationRate(k1, k2, KeldyshFunctionResult, EgapEff, wavelength)
print "w_PI until order "+str(order)+" = ", wPI

plt.figure()
plt.subplot(311)
plt.xlabel("Field (V/m)")
plt.ylabel("Excitation rate w_{PI} (m^{-3} s^{-1})")
plt.loglog(Efield, wPI, linestyle="-", color="k", label="w_{PI}")
plt.legend()

plt.subplot(312)
plt.xlabel("Field (V/m)")
plt.ylabel("Adiabadicity parameter")
plt.loglog(Efield, gamma, color="k", linestyle="-", label="gamma")
# plt.loglog(Efield, 0.1, label="Tunnelling limit")
# plt.loglog(Efield, 10.*np.ones(), label="MPI limit")
plt.legend()

plt.subplot(313)
plt.xlabel("Intensity (W/m^{2})") #Field (V/m)")
# plt.xlabel("Field (V/m)")
plt.ylabel("Density estimation (m^{-3})")
plt.loglog(FieldToIntensity(Efield), wPI*tau)
plt.savefig("KeldyshAnalytic.eps")
plt.show()
