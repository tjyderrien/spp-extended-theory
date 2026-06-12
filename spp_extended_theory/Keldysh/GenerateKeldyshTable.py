#!/usr/bin/env python
#-*- coding: utf-8 -*-
## @package GenerateTablesKeldysh
## Computes the Keldysh excitation rate of quasi-free electrons in a format compatible with interaction codes
# This module aims to calculate the density of excited electrons as function of laser parameters. 
# Two types of usage are planned :
# * Generating tables to use directly into simulation codes
import libKeldysh
from libKeldyshZhukov import IntensityToField
from spp_extended_theory.Libs.libDatabase import ExportToTxt
import numpy as np
from scipy.constants import pi, e

print("")
print("** Welcome to spp_extended_theory suite.")
print("** Author(s): T.J.-Y. Derrien")
print("")
print("** Loading Keldysh module [Keldysh, Sov. J. Exp. Th. Phys. 47, 5 (1964)]...")
print("** Loading Gruzdev formula [Gruzdev, Optical Engineering 53, 122515 (2014)]")

print("Defining material parameters...")

#Egap = 2.58e0*e; #LDA band gap of Si: 2.58 eV. #1.12e0*e for indirect band gap; 
meff=0.18e0; #Effective mass of Si
#Ntotal=1.*5E28

#wavelength = np.array([800e-9])
tau            = np.array([10e-15]);
PeakFluence    = 1E4*np.array([1e-3, 1e-2, 1e-1, 1E0]) #J/cm2 * 1E4 = J/m2
PeakIntensity  = np.divide(PeakFluence, np.divide(tau, np.sqrt(4e0*np.log(2e0)/pi)))
PeakField      = IntensityToField(PeakIntensity)
order = 50
ShowPlot = False

#print "E_peak = "+str(PeakField)+" V/m."

dEgap = 0.05e0 #eV
Egap = np.array([[np.multiply(e, np.arange(dEgap, 10e0, dEgap))]])
#Egap = np.array([e*1.12])

# Build the arrays
wavelength = np.array([1030e-9, 800e-9, 515e-9, 343e-9])
#PeakField = np.array([np.arange(1E4, 1E10, 1E9)]) #

# Building the tensors
wavelength = np.array([[wavelength]])
PeakField  = np.array([[PeakField]])

# Transposing for assisting the vectorization
wavelength = np.transpose( wavelength, (2,0,1) )
PeakField  = np.array(np.transpose( PeakField , (1,2,0) ))

# Checking before multiplexing
#print wavelength.shape
#print PeakField.shape
#print Egap.shape

# Computation
Database = libKeldysh.GenerateKeldyshDatabase(Egap, meff, wavelength, PeakField, order)

# Checking on stdout
print(Database)

# Output to a file
ExportToTxt(Database, "Keldysh.dat")
# TODO: normalize this format for C routines and Fortran routines, to plug into our codes. 