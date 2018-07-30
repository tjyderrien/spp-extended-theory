#!/usr/bin/env python2
#-*- coding: utf-8 -*-

# Copyright (C) 2013-2018 T. J.-Y. Derrien
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

## @package libThermalPropertiesMaterials
# Module libThermalPropertiesMaterials defines the thermodynamic properties of some materials. 

# TODO: this module could be helped by: 
# https://pypi.python.org/pypi/pycalphad/0.5.1
# https://github.com/guillemborrell/thermopy

from libSPP import *

from scipy.constants import h, hbar, e, gravitational_constant, Boltzmann
gravity = gravitational_constant
k_b     = Boltzmann

#print "Boltzmann constant: "+str(k_b)+" J.s."

Header = "[libThermalPropertiesMaterials] "

def Silica_Solid_VolumicMass():
  return 2.5e3 #kg/m3

def Silica_Liquid_VolumicMass():
  return 2.2e3 #kg/m3

def Silica_Solid_ThermalConductivity(T):
  if(T >= Silica_MeltingTemperature()): 
    print Header+"** Warning: Silica_Solid_ThermalConductivity() is used out validity range. "
  result = 0.14E2 #W/m/K [Bauerle data]
  return result

def Silica_Liquid_ThermalConductivity(T):
  if(T < Silica_MeltingTemperature()): 
    print Header+"** Warning: Silica_Liquid_ThermalConductivity() is used out of its validity range. "
  return 0.014E2 #W/m/K [Bauerle data]

def Silica_Solid_HeatDiffusivity(T):
  if(T >= Silica_MeltingTemperature()): 
    print Header+"** Warning: Silica_Solid_HeatDiffusivity() is used out of its validity range. "
  return 0.086e-4 #m2/s

def Silica_Liquid_HeatDiffusivity(T):
  if(T < Silica_MeltingTemperature()): 
    print Header+"** Warning: Silica_Liquid_HeatDiffusivity() is used out of its validity range. "
  return 0.009E-4 #m2/s

def Silica_Solid_HeatCapacity(T):
  if(T >= Silica_MeltingTemperature()): 
    print Header+"** Warning: Silica_Solid_HeatCapacity() is used out of its validity range. "
  # to be multiplied by density!
  return 0.74e3 #J / kg / K

def Silica_Liquid_HeatCapacity(T):
  # to be multiplied by density!
  if(T < Silica_MeltingTemperature()): 
    print Header+"** Warning: Silica_Liquid_HeatCapacity() is used out of its validity range. "
  return 0.72E3 #J / kg / K

## Thermal Conductivity of silica, for a wide range of temperatures. 
#Fitted on Wray, Kurt L. and Connolly, Thomas J., "Thermal Conductivity of Clear Fused Silica at High Temperatures", Journal of Applied Physics (1959), 1702--1705.
def Silica_HeatConductivity(T):
  if(T < 300. or T > 2000.):
    print Header+"** Warning: Silica_HeatConductivity() was used out of its validity range."
  a3 = 7.06418e-10
  a2 = -3.96976e-6
  a1 = 7.56664e-3
  a0 = 0.633319
  return a3*T**3 + a2*T**2 + a1*T + a0

## Surface tension of liquid silica (in N/m)
# @param T: temperature (K)
#Fitted on Boyd K et al., "Surface tension and viscosity measurement of optical glasses using a scanning CO 2 laser", Optical Materials Express (2012), 1101--1110.
#BUG: validity range? 
def Silica_SurfaceTension(T):
  a = 1.54E-5; b=0.267; 
  return a*T+b, a

## Silica melting temperature (in K)
# Source? 
def Silica_MeltingTemperature(): 
  return 1300. 

## Dynamic viscosity of fused silica (in Pa.s)
#Fitted on Urbain et al, "Viscosity of liquid silica, silicates and alumino-silicates", Geochimica et Cosmochimica Acta (1982), 1061--1072.
# @param T: temperature (K). Validity range: 1300-2000 K 
def Silica_DynamicViscosity(T):
  Tm_SiO2   = Silica_MeltingTemperature() #K
  Tmax_SiO2 = 2000e0 #K [arbitrary?]
  if(T < Tm_SiO2 or T>Tmax_SiO2):
    print Header+"** Warning: dynamic viscosity was taken out of range (SiO2 temperature must be liquid). "
    result = 0e0
  else: 
    a=6.23888379570223e0; b=-14.6668118241314e0; 
    result = 0.1*np.exp(1E4*a/T+b)
  return result

## Returns temperature-dependent band-gap energy of SiO2 (in eV).
# @param T: temperature (in K). Validity range: 300 K - 2000 K.
# Saito, K. & Ikushima, A. J. Absorption edge in silica glass Physical Review B, 2000, 62, 8584
def Silica_BandGapEnergy(T):
  Egap0  = 8.52 #eV
  L0     = 10.3e0
  omega0 = 0.079*e / hbar
  X      = 0.33e0
  return Egap0 - L0*( hbar * omega0 / e * (0.5e0 + 1E0/( np.exp(hbar*omega0 / kb / T) - 1E0 ) ) + 0.5E0 * X * hbar * omega0 )

Silica_HeatConductivity = np.vectorize(Silica_HeatConductivity)
Silica_SurfaceTension = np.vectorize(Silica_SurfaceTension)
Silica_DynamicViscosity = np.vectorize(Silica_DynamicViscosity)
Silica_BandGapEnergy = np.vectorize(Silica_BandGapEnergy)
