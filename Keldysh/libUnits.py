#!/usr/bin/env python2.7
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

import numpy as np
from scipy.constants import c, pi

## @package libUnits
# Provides routines to safely convert in SI units or CGS units. 

c_SI = c

# Convert length units from CGS to SI
# Validated on https://en.wikipedia.org/wiki/Centimetre%E2%80%93gram%E2%80%93second_system_of_units#Electromagnetic_units_in_various_CGS_systems
def Length_CGS_to_SI(CGS):
  return CGS / 1e2

## Convert length units from SI to CGS
# Validated on https://en.wikipedia.org/wiki/Centimetre%E2%80%93gram%E2%80%93second_system_of_units#Electromagnetic_units_in_various_CGS_systems
def Length_SI_to_CGS(SI):
  return SI * 1e2

## Converts a mass in g (CGS unit) to kg (SI)
# @param CGS: mass in g (CGS unit)
# Validated on https://en.wikipedia.org/wiki/Centimetre%E2%80%93gram%E2%80%93second_system_of_units#Electromagnetic_units_in_various_CGS_systems
def Mass_CGS_to_SI(CGS):
  return CGS * 1E-3

## Converts a mass from kg (SI) to g (CGS unit)
# @param SI: mass in kg (SI unit)
# Validated on https://en.wikipedia.org/wiki/Centimetre%E2%80%93gram%E2%80%93second_system_of_units#Electromagnetic_units_in_various_CGS_systems
def Mass_SI_to_CGS(SI):
  return SI * 1E3

## Converts velocity from cm/s (CGS unit) tp m/s (SI unit)
# @param CGS: velocity in cm/s
# Validated on https://en.wikipedia.org/wiki/Centimetre%E2%80%93gram%E2%80%93second_system_of_units#Electromagnetic_units_in_various_CGS_systems
def Velocity_CGS_to_SI(CGS):
  return CGS / 1E2

## Converts velocity from m/s (SI unit) to cm/s (CGS unit)
# @param SI: velocity in m/s
# Validated on https://en.wikipedia.org/wiki/Centimetre%E2%80%93gram%E2%80%93second_system_of_units#Electromagnetic_units_in_various_CGS_systems
def Velocity_SI_to_CGS(SI):
  return SI * 1E2

## Converts an energy in ergs (CGS unit) to Joules (SI unit)
# @param CGS: energy in ergs (CGS unit)
# Validated on https://en.wikipedia.org/wiki/Centimetre%E2%80%93gram%E2%80%93second_system_of_units#Electromagnetic_units_in_various_CGS_systems
def Energy_CGS_to_SI(CGS):
  return CGS / 1E7

## Converts an energy from Joules (SI unit) to ergs (CGS unit) 
# @param SI: energy in Joules (SI unit)
# Validated on https://en.wikipedia.org/wiki/Centimetre%E2%80%93gram%E2%80%93second_system_of_units#Electromagnetic_units_in_various_CGS_systems
def Energy_SI_to_CGS(SI):
  return SI * 1E7

## Converts electric charge in Coulomb (SI unit) to statC (CGS unit)
# Validated on Jackson book: 1 C ~ 3E9 statC 
def electric_charge_SI_to_CGS(SI):
  c_CGS = Velocity_SI_to_CGS(c_SI)
  conversion = c_CGS / 10.
  return SI * conversion

## Converts electric charge in statC (CGS unit) to Coulomb (SI unit)
# Validated on Jackson book: 1 C ~ 3E9 statC 
def electric_charge_CGS_to_SI(CGS):
  c_CGS = Velocity_SI_to_CGS(c_SI)
  conversion = c_CGS / 10.
  return CGS / conversion

## Converts field CGS units (statV/cm) in SI (V/m).
# @param CGS: input field in CGS units
# Retuns the field in SI units (V/m).
# Jackson: 1 V/m ~ 1 / 3 * 1E-4 
#                = 1E8 / c_SI * 1E-4 = 1E2 / c_SI
#                = 1E6 / c_CGS
def Field_CGS_to_SI(CGS):
  c_CGS      = Velocity_SI_to_CGS(c_SI)
  conversion = 1E6 / c_CGS
  return CGS/conversion

## Converts field SI units (V/m) to CGS units (statV/cm)
# @param SI: input field in SI units (V/m)
# Retuns the field in CGS units (statV/cm).
# Jackson: 1 V/m ~ 1 / 3 * 1E-4 
#                = 1E8 / c_SI * 1E-4 = 1E2 / c_SI
#                = 1E6 / c_CGS
def Field_SI_to_CGS(SI):
  c_CGS      = Velocity_SI_to_CGS(c_SI) #[cm/s]
  conversion = 1E6 / c_CGS
  return SI*conversion

def Energy_SI_to_Length(energy):
  wavelength = 2.*pi*c/energy
  return wavelength

def Length_SI_to_energy(wavelength):
  energy = 2.*pi*c/wavelength
  return energy

def Energy_Joules_to_eV(energy):
  return energy / e

def Energy_eV_to_Joules(energy):
  return energy * e

Energy_eV_to_Joules = np.vectorize(Energy_eV_to_Joules)
Energy_SI_to_Length = np.vectorize(Energy_SI_to_Length)
