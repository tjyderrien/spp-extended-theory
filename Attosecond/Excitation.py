#!/usr/bin/env python
#-*- coding: utf-8 -*-

# Copyright (C) 2013-2020 T. J.-Y. Derrien
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>

## @package Excitation
## Excitation of materials is calculated for ZnO (for now) and later with any band gap material. 
# Basic inputs should be: 
# * Egap (J)
# * Wavelength
# * Fluence, pulse duration
# 
# Outputs will be:
# * (n,k) as function of time
# * max{(n,k)} as function of laser fluence. 

from libKeldysh import *

from libSPP import *

# Desired laser parameters
wavelength = 800e-9
tau=40e-15; dt = 1E-18
PeakFluence = 1.0*1E4  / 25. #np.multiply(np.arange(0.1,1.,0.1), 1E4) #J/cm2 * 1E4 = J/m2
order = 50
ShowPlot = True

## Extracting optical data for ZnO...
database = "MaterialOpticalDatabaseForPlasmonics.csv" #Import the database
MaterialDB = loadtxt(database, dtype='str', delimiter='\t')
query = "ZnO"
try:
  MaterialDB = FilterDatabaseContains(MaterialDB, query, 0) # Find ZnO for the given wavelength
  print("Refined database...")
  print(MaterialDB)
except: 
  print("** Error: Material is not contained into Materials database. ")
  exit()
  
query = str(int(wavelength*1E9))
try:
  MaterialDB = FilterDatabaseContains(MaterialDB, query, 2)
  print("Refined database...")
  print(MaterialDB)
except: 
  print("** Error: this material does not contain the required wavelength. ")

# Extract ZnO data from database
ZnO_name, ZnO_bandgap, ZnO_wavelength, ZnO_epsilonRe, ZnO_epsilonIm = ExtractMaterialData(MaterialDB)

Egap = float(ZnO_bandgap[0])*e; #band gap of ZnO: 3.3 eV. #Si: 1.12e0*e for indirect band gap; 

# Missing materials properties (TODO: to be added in database!)
collisionRate = 1E15 #TODO: arbitrary
meff=0.29e0; # Effective mass of ZnO
density_ZnO = 5.61E3 #[kg/m^3]
N_avogadro = 6.02E23 #[mol^-1]
MassMol = 81.408E-3 #[kg/mol]
N_total = density_ZnO / MassMol * N_avogadro # (Density [kg*m^-3] / MassMol [kg/mol] = mol / m3 ) * N_avogadro [mol^-1] = Density [m^-3]
print(("** Info: Maximum excitation density: "+str(N_total)))


print("** Generating tables of excited electron density...")
instants, N_excited_Keldysh, N_excited_Gruzdev = generateWpiTables(Egap, meff, wavelength, tau, PeakFluence, dt, order) #wPIg as function of Efield amplitudes (do we need to put the pulse there?)

ShortRefKeldysh = "[Keldysh (1964)]"
ShortRefGruzdev = "[Gruzdev (2014)]"

#TODO: shall we compute N_exc(t) for each pulse intensity? or make a simple law? Rather compute the whole thing. 
#N_exc = 
print("")
print("** Warning: results may be not converged.")
print("            Reduce dt, and increase order until convergence.")
print("")
#print "Maximum density N_ex "+ShortRefKeldysh+" = "+str(N_excited_Keldysh.max())+"."
print(("Maximum density N_ex "+ShortRefGruzdev+" = "+str(N_excited_Gruzdev.max())+"."))
print("")

print("** Pluging these excitations into dielectric permittivity change. ")


print("** Get the ZnO dielectric permittivity at equilibrium in database...")
#print ZnO_epsilonRe[0], ZnO_epsilonIm[0]
epsilon_inf = float(ZnO_epsilonRe[0]) + 1j * float(ZnO_epsilonIm[0])
print("** Excitation of the ZnO...")
epsilon_exc = Drude(wavelength, N_excited_Gruzdev.max(), epsilon_inf, collisionRate, meff)
print(("** Epsilon for excited ZnO (maximum value only): "+str(epsilon_exc)+" ."))
#plotPulseToDensity(Egap, meff, wavelength, tau, PeakFluence, dt, order, ShowPlot, 0e0, Ntotal)
