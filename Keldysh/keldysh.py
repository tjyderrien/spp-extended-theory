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

wavelength = 800e-9
tau=10e-15; dt = 1E-17; CEP=0e0
PeakFluence = 0.01*1E4 #J/cm2 * 1E4 = J/m2
PeakField   = np.sqrt(2e0 * PeakFluence / (tau * c * epsilon_0))

tmin=-1.*tau; tmax=1.*tau
t0=0.

instants = np.arange(tmin, tmax, dt)
#print "Time range: "+str(instants.min())+", "+str(instants.max())+"."

# Test with a single pulse
FieldEnvelope, RealField = PulseSquaredSinTemporalShape(instants, tau, PeakField, wavelength, CEP, t0)

# Test with a double pulse
#Delay = 0.
#FieldEnvelope, RealField = PulseSquaredSinTemporalShapeDoublePulse(instants, tau, tau, PeakField, PeakField, wavelength, wavelength, CEP, CEP, t0, Delay)

#plt.plot(instants, RealField.real)
#plt.show()
#print FieldEnvelope
#exit()

order = 50
ShowPlot = True

timeKeldysh, N_Keldysh_SI, N_Gruzdev_SI = plotPulseToDensity(Egap, meff, wavelength, tau, FieldEnvelope.real, dt, order, ShowPlot, 0e0, Ntotal)

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
