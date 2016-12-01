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

Egap = 2.58e0*e; #LDA band gap of Si: 2.58 eV. #1.12e0*e for indirect band gap; 
meff=0.18e0; #Effective mass of Si
Ntotal=1.*5E28

wavelength = 800e-9
tau=10e-15; dt = 1E-18
PeakFluence = 100e-3*1E4 #J/cm2 * 1E4 = J/m2
order = 50
ShowPlot = True

plotPulseToDensity(Egap, meff, wavelength, tau, PeakFluence, dt, order, ShowPlot, 0e0, Ntotal)

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