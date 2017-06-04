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
print "** Author(s): T.J.-Y. Derrien"
print ""
print "** Loading Keldysh module [Keldysh, Sov. J. Exp. Th. Phys. 47, 5 (1964)]..."
print "** Loading Gruzdev formula [Gruzdev, Optical Engineering 53, 122515 (2014)]"

print "Defining material parameters..."

Egap = 2.56e0*e; #LDA band gap of Si: 2.58 eV. #1.12e0*e for indirect band gap; 
meff=0.2226e0; #Effective mass of Si
Ntotal=1.*5E28

wavelength = 800e-9; wavelength2 = 400e-9
tau=10e-15; dt = 1E-17; CEP=0e0
PeakFluence = 0.01*1E4 #J/cm2 * 1E4 = J/m2
PeakField   = 4655327068.03 # np.sqrt(2e0 * PeakFluence / (tau * c * epsilon_0))

t0=0. #defines the instant 0.
Delay = 0. #delay between maxima of the pulses
tmin=-1.*tau + t0; tmax=1.*tau + Delay + t0

instants = np.arange(tmin, tmax, dt)
#print "Time range: "+str(instants.min())+", "+str(instants.max())+"."

PeakField2  = PeakField #/ 2.
CEP2        = 0. #pi/3.
#wavelength2 = wavelength

# Test with a single pulse centered on 0
FieldEnvelope1, RealField1 = PulseSquaredSinTemporalShape(instants, tau, PeakField, wavelength, CEP, t0, 0.)

# We build a second pulse with a delay
FieldEnvelope2, RealField2 = PulseSquaredSinTemporalShape(instants, tau, PeakField2, wavelength2, CEP, t0, Delay)

print "# Test with a bicolor double pulse"

FieldEnvelopeTot, RealFieldTot = PulseSquaredSinTemporalShapeDoublePulse(instants, tau, tau, PeakField, PeakField, wavelength, wavelength2, CEP, CEP2, t0, Delay)

plt.plot(instants, RealField1.real, '-')
plt.plot(instants, FieldEnvelope1.real, '--')
plt.plot(instants, RealField2.real, '-')
plt.plot(instants, FieldEnvelope2.real, '--')
plt.plot(instants, RealFieldTot.real, '-')
plt.plot(instants, FieldEnvelopeTot.real, '--')
plt.xlabel('')
plt.savefig('PulseEnvelopes.eps')
plt.savefig('PulseEnvelopes.png')
#plt.show()

print "Computing the W_PI values..."
wPI_Zhukov = VZ_generateWpiTables(FieldEnvelope1, FieldEnvelope2, wavelength, wavelength2, CEP, CEP2, Egap, meff, tau, tau, Delay, dt, Ntotal, t0)

exit()

order = 50
ShowPlot = True

#timeKeldysh, N_Keldysh_SI, N_Gruzdev_SI = plotPulseToDensity(Egap, meff, wavelength, tau, FieldEnvelope.real, dt, order, ShowPlot, 0e0, Ntotal)

# print "Checking dt convergence..."
# plotPulseToDensity(Egap, meff, wavelength, tau, PeakFluence, 1E-17, order, ShowPlot)
# plotPulseToDensity(Egap, meff, wavelength, tau, PeakFluence, 5E-17, order, ShowPlot)
# plotPulseToDensity(Egap, meff, wavelength, tau, PeakFluence, 1E-16, order, ShowPlot)
# print ""
# print "Checking order convergence..."
# plotPulseToDensity(Egap, meff, wavelength, tau, PeakFluence, dt, 10, ShowPlot)
# plotPulseToDensity(Egap, meff, wavelength, tau, PeakFluence, dt, 20, ShowPlot)
# plotPulseToDensity(Egap, meff, wavelength, tau, PeakFluence, dt, 30, ShowPlot)
# plotPulseToDensity(Egap, meff, wavelength, tau, PeakFluence, dt, 40, ShowPlot)
# plotPulseToDensity(Egap, meff, wavelength, tau, PeakFluence, dt, 50, ShowPlot)

#TODO: print "Exporting density to a table for ZnO optical index calculations."
