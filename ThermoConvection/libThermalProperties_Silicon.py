#!/usr/bin/env python2
#-*- coding: utf-8 -*-

# Copyright (C) 2013-2019 T. J.-Y. Derrien
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

## SILICON

def MeltingTemperature():
  return 1687. #K

def Liquid_Density():
    return 5.6368E28 #m-3

def Solid_VolumicMass():
  return 2.329e3 #kg/m3

def Liquid_MolarMass():
  return 

def Liquid_VolumicMass():
  return 2.553e3 #kg/m3

def Liquid_ThermalConductivity(T):
  Tm = MeltingTemperature()
  if(T<Tm):
      print "** Warning: K_lSi is used beyond its boundaries."
  return 1E2*(0.5+29.3E-5*(T-Tm))

def Liquid_HeatCapacity(T):
    Tm = MeltingTemperature()
    if(T<Tm):
        print Header+"** Warning: HeatDiffusivity is used beyond its boundaries."
    return 1045*Liquid_Density()

def Liquid_HeatDiffusivity(T):
    Tm = MeltingTemperature()
    if(T<Tm):
        print Header+"** Warning: HeatDiffusivity is used beyond its boundaries."
    return np.sqrt(np.divide(Liquid_ThermalConductivity(T), Liquid_HeatCapacity(T)))

def Solid_HeatConductivity(T):
    Tm = MeltingTemperature()
    if(T>Tm):
        print Header+"** Warning: HeatDiffusivity is used beyond its boundaries."
    return 1E2*(1585*np.power(T, -1.23))

def Solid_HeatCapacity(T):
    Tm = MeltingTemperature()
    if(T>Tm):
        print Header+"** Warning: HeatDiffusivity is used beyond its boundaries."
    return 1E6*(1.978+3.54E-4*T-3.68*T**-2)

def SolidHeat_Diffusivity(T):
    Tm = MeltingTemperature()
    if(T>Tm):
        print Header+"** Warning: HeatDiffusivity is used beyond its boundaries."
    return np.sqrt(np.divide(SolidHeatConductivity(T),SolidHeatCapacity(T)))

def Liquid_SurfaceTension(T):
    Tm = MeltingTemperature()
    if(T<Tm):
        print Header+"** Warning: HeatDiffusivity is used beyond its boundaries."
    return 0.885-0.28e-3*(T-Tm), 0. #N/m

## Dynamic viscosity of fused silica (in Pa.s)
#Fitted on Urbain et al, "Viscosity of liquid silica, silicates and alumino-silicates", Geochimica et Cosmochimica Acta (1982), 1061--1072.
# @param T: temperature (K). Validity range: 1300-2000 K 
def DynamicViscosity(T):
  Tm   = MeltingTemperature() #K
  Tmax = 2000e0 #K [arbitrary?]
  if(T < Tm or T>Tmax):
    print Header+"** Warning: dynamic viscosity was taken out of range (materials temperature must be liquid). "
    result = 0e0
  else: 
    Volume = np.power(Liquid_Density()*N_Av/MolarMass,-1)
    result = h * N_Av*np.exp(3.8*T_vap/T)/Volume
  return result


#Liquid_SurfaceTension = np.vectorize(Liquid_SurfaceTension)
