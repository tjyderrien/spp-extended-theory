#!/usr/bin/env python
#-*- coding: utf-8 -*-

# Copyright (C) 2013-2017 T. J.-Y. Derrien
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

## @package libThermoCapillary
# Module libThermoCapillary defines the various regimes of thermo capillary instabilities found in the LIPSS formation. 
# We use as the following basics: 
# Chapter IV of PhD thesis: Derrien, T. J.-Y., Nanostructuring of solar cells by femtosecond laser irradiation. Theoretical study of the formation mechanisms. Université de la Méditerranée - Aix Marseille II, 2012. 
# Jean Berthier and Pascal Silberzan, "Microfluidics for Biotechnology", Artech House (2009).

from libThermalPropertiesMaterials import *

def Term1(T, k, depth, dynamic_viscosity, density, surface_tension):
  num1 = dynamic_viscosity * k**2
  num2 = np.sqrt( (np.exp(k*h)**2 + 1e0 ) * (np.exp(k*h)**2 - 1e0) )
  num3 = np.sqrt( density * np.sqrt( gravity * depth * density + surface_tension * k**3) ) 
  num  = - np.sqrt(2e0) * np.sqrt( num1 * num2 * num3  )
  denom = 2e0*np.sinh(2e0*k*h)
  return num / denom

def Term2(T, k, surface_tension, Absorptivity, laser_intensity, dynamic_viscosity):
  num1 = np.diff(surface_tension)/np.diff(T) * Absorptivity * laser_intensity / thermal_conductivity(T)
  num2 = k**2 * (gravity * thickness * density + surface_tension * k**3 )
  num3 = np.sqrt( thermal_conductivity(T) * k**2 * density * ( np.exp(k*depth)**2 + 1E0) / ( ( heat_capacity * density * np.sqrt( (  np.exp(k*thickness)**2 + 1) * np.exp(k*thickness)**2 - 1e0 )) * np.sqrt(density * (gravity*thickness*density+surface_tension*k**3)) ) )
  num4 = np.sqrt(2.) * (np.exp(k*thickness)**2+1e0)
  
  num = num1 * num2 * num3 * num4 
  
  denom1 = 4e0*density*np.sqrt( (np.exp(k*thickness)**2 + 1e0) * ( np.exp(k*thickness)**2 - 1e0 ) )
  denom2 = np.sqrt(density * (gravity*thickness*density+surface_tension*k**3))
  denom3 = ( (np.exp(k*thickness)**2-1e0) * (gravity*thickness*density+surface_tension*k**3) ) / (density * (np.exp(k*thickness)**2+1e0))
  denom4 = ( k**2 * np.abs(np.diff(surface_tension) / np.diff(T)) * Absorptivity * laser_intensity / thermal_conductivity ) / (density * (1e0 + np.sqrt( dynamic_viscosity * density * heat_capacity / (density * thermal_conductivity) ) ))
  denom = denom1 * denom2 * ( denom3 - denom4 )
  return num / denom

## Growth rate of thermoconvective instability 
# Reference is intentionally not given in order to avoid disclosing the model for now. 
def InstabilityGrowthRate(T, k, depth, dynamic_viscosity, density, surface_tension, Absorptivity, laser_intensity):
  return Term1(T, k, depth, dynamic_viscosity, density, surface_tension) / Term2(T, k,  surface_tension, Absorptivity, laser_intensity, dynamic_viscosity)
  
Term1                 = np.vectorize(Term1)
Term2                 = np.vectorize(Term2)
InstabilityGrowthRate = np.vectorize(InstabilityGrowthRate)

laser_wavelength = 1025e-9
laser_fluence = 4E4 #J/m2
laser_FWHM = 300e-15 #s

# BE VERY CAREFUL WITH UNITS !
T=np.arange(300.,400., 50.) #length should be greater than 1, strictly. 
print T
surface_tension = Silica_SurfaceTension(T)
print surface_tension
exit()


k = np.arange(0., 2.*pi/laser_wavelength / 2., 1E6)
thickness = 100e-9
dynamic_viscosity = Silica_DynamicViscosity(T)
density = Silica_Liquid_VolumicMass() #TODO: CHECK UNIT
Absorptivity = 1. #total
laser_intensity = laser_fluence / laser_FWHM

print Term1(T, k, thickness, dynamic_viscosity, density, surface_tension)
print Term2(T, k, surface_tension, Absorptivity, laser_intensity, dynamic_viscosity)