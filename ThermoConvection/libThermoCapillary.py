#!/usr/bin/env python2
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
from libNonDimensionalNumbers import *

#============== THESE ROUTINES ARE MESSY AND DONT EXACTLY FOLLOW THE SIMPLEST FORMULATION GIVEN IN ORIGINAL PAPER. 
def Term1(T, k, depth, dynamic_viscosity, density, surface_tension):
  num1 = dynamic_viscosity * k**2
  num2 = np.sqrt( (np.exp(k*h)**2 + 1e0 ) * (np.exp(k*h)**2 - 1e0) )
  num3 = np.sqrt( density * np.sqrt( gravity * depth * density + surface_tension * k**3) ) 
  num  = - np.sqrt(2e0) * np.sqrt( num1 * num2 * num3  )
  denom = 2e0*np.sinh(2e0*k*h)
  return num / denom

#TODO: Two ways can be used to estimate the derivative: numerically, or by taking analytical polynomial law for surface tension (but then it's exact).
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

#BUG: vectorization does not work yet. But let's validate the routine first.   
#Term1                 = np.vectorize(Term1)
#Term2                 = np.vectorize(Term2)
#InstabilityGrowthRate = np.vectorize(InstabilityGrowthRate)

#============== THESE ROUTINES DO FOLLOW THE ORIGINAL PAPER. 

## kinematic_viscosity [m2/s]: #TODO: check in fluid dynamics book!  
# @param dynamic_viscosity [Pa.s]
# @param mass_density [kg/m3]
def KinematicViscosity(dynamic_viscosity, mass_density):
  return dynamic_viscosity / mass_density

# Physical meaning? 
# Source: Levchenko paper
def omega0_Levchenko(density, gravity, thickness, surface_tension):
  result = np.sqrt(density * (gravity * thickness * density + surface_tension * k **3) ) / density
  return result

# Quantity omega (physical meaning?)
def omega_Levchenko(omega_0, thickness, k):
  # Formulation from Anisimov book
  formula1 = np.sqrt(np.exp(4.*k*thickness) - 1.) * omega_0 / (np.exp(2.*k*thickness)+1.) * omega_0
  # Formulation from original Levchenko article
  formula2 = np.sqrt(omega_0 ** 2 * np.tanh(k*thickness))
  return formula2

def SurfaceTensionDerivation(surface_tension_model=0):
  # We have two possile models to compute derivative of surface tension
  surface_tension, surface_tension_deriv = Silica_SurfaceTension(T) 
  if(surface_tension_model == 0): #analytic approach
    surface_tension_diff = surface_tension_deriv
  elif (surface_tension_model == 1):
    surface_tension_diff = surface_tension.diff() / T.diff() #numeric derivation
  else: 
    print Header+"** Error. Surface tension computation in SoundVelocity_Levchenko function."
  return surface_tension_diff

## Interfaces the laser with the model
# @param intensity: [in W/m2] laser intensity
# @param T [in K] temperature of the liquid layer
# @param thermal_conductivity [W/m/K] of the irradiated liquid
# @param absorption (ratio of absorbed optical energy [0-1]
# @param ModelIndex: 0 is for pulse>=ns duration, 1 is for pulse<=ps duration.
# @param e_ph_coupling_rate [W/m3/K] electron-phonon coupling rate within two-temperature-model
# @param Te [K] electron temperature. Taken equal to liquid temperature by default. 
def ThermalSourceModels(intensity, thermal_conductivity, absorption=1., ModelIndex=0,  e_ph_coupling_rate=0., Te=300., T=300.):
    if(ModelIndex==0): #assumes that all optical energy is transfered to heating. Neglects the potential energy. 
      print Header+"** Info: choosed ns+ absorption model in ThermalSourceModels."
      ThermalSource = absorption * intensity / thermal_conductivity
    elif(ModelIndex==1): #valid for fs and ps pulse durations, based on TTM
      print Header+"** Info: choosed fs-ps absorption model in ThermalSourceModels."
      ThermalSource = e_ph_coupling_rate * (Te - T)/thermal_conductivity
    else:
      print Header+"** Error: ThermalSource failed."
    return ThermalSource

## Compute the sound velocity from Levchenko paper #TODO: check physical meaning
# @param density [kg/m3] is the liquid density. 
# @param Peclet [no dimension] can be computer with PecletNumber()
# @param T [in K] temperature of the liquid
# @param surface_tension_diff must be computer with SurfaceTensionDerivation()
# @param ThermalSource_diff must be computed with function ThermalSourceModels()
def SoundVelocity_Levchenko(density, Peclet, T, surface_tension_diff, ThermalSource_diff):
  num   = np.sqrt( density * (1.+np.sqrt(Peclet))*surface_tension_diff * ThermalSource_diff)
  denom = density * (1. + np.sqrt( Peclet ) ) 
  return num / denom

## Growth rate of thermoconvective instability [Levchenko et al]
def InstabilityGrowthRate_Levchenko(kinematic_viscosity, k, omega, thickness, sound_velocity, omega_0, Peclet, diffusivity):
  Term1 = np.divide( np.sqrt(np.multiply(np.multiply(kinematic_viscosity, np.multiply(k,k)), omega) ), (2.**0.5 * np.sinh(2.*k*thickness) ) )
  Term2 = sound_velocity**2 * k**2 * omega_0 * (1. + np.sqrt(Peclet)) / (2.*np.sqrt(2.) * omega * (np.multiply(omega,omega) - sound_velocity**2 * np.multiply(k,k) ))
  Term3 = np.sqrt( np.divide(np.multiply(diffusivity, np.multiply(k,k)), omega) )
  gamma = - Term1 + np.multiply(Term2, Term3)
  return gamma

#============== HERE WE SIMPY USE THE ROUTINES

laser_wavelength = 1025e-9 #m
laser_fluence = 4E4 #J/m2
laser_FWHM = 300e-15 #s
thickness = 50e-9 #molten depth thickness [m]

#T=np.arange(1300.,2000., 10.) #length should be greater than 1, strictly. 
T = 2000. #K
print Header+"Temperature of the liquid [K]: "+str(T)

surface_tension, surface_tension_deriv = Silica_SurfaceTension(T)
print Header+"Surface tension [N.m2]: "+str(surface_tension)

k_laser = 2.*pi/laser_wavelength #m-1
k = np.arange(k_laser/10., 10.*k_laser, k_laser/100.) #NOTE: this is vector style. 
#k = k_laser #index-based programming style

print Header+"** Selected modes (1/m): "+str(k)+" 1/m, equiv. to "+str(1E9*2.*pi/k)+" nm."

#if( T < Silica_MeltingTemperature() ): 
  #print Header+"** Absurd: silica should reach melting temperature. "
  #exit()
  
Absorptivity = 1. #total
laser_intensity = laser_fluence / laser_FWHM

print Header+"Laser absorptivity: "+str(100.*Absorptivity)+" %"
print Header+"Laser intensity   : "+str(laser_intensity*1E-4)+ "W/cm2."

mass_density         = Silica_Liquid_VolumicMass()
dynamic_viscosity    = Silica_DynamicViscosity(T)
diffusivity          = Silica_Liquid_HeatDiffusivity(T)
thermal_conductivity = Silica_Liquid_ThermalConductivity(T)

print Header+"== Silica data =="
print Header+"Mass density: "+str(mass_density)+" kg/m3."
print Header+"Dynamic viscosity: "+str(dynamic_viscosity)+"Pa.s"
print Header+"Heat diffusivity: "+str(diffusivity)+" m2/s."
print Header+"Thermal conductivity: "+str(thermal_conductivity)+" W/m/K."

#print Header+"** PhD thesis complicated formulas"
#print Term1(T, k, thickness, dynamic_viscosity, mass_density, surface_tension) #TODO: compare with Maple
#print Term2(T, k, surface_tension, Absorptivity, laser_intensity, dynamic_viscosity) #TODO: compare with Maple<
#print InstabilityGrowthRate(T, k, depth, dynamic_viscosity, mass_density, surface_tension, Absorptivity, laser_intensity)

kinematic_viscosity  = KinematicViscosity(dynamic_viscosity, mass_density)
print Header+"Kinematic viscosity: "+str(kinematic_viscosity)+" m2/s."
print ""
print Header+"** Levchenko original paper formulas"
Peclet_Number        = PecletNumber(kinematic_viscosity, diffusivity)
omega_0              = omega0_Levchenko(mass_density, gravity, thickness, surface_tension)
omega                = omega_Levchenko(omega_0, thickness, k)

print Header+"Peclet number: "+str(Peclet_Number)+"."
print Header+"omega_0_Levchenko: "+str(omega_0)+" Hz."
print Header+"omega_Levchenko: "+str(omega)+" Hz."

surface_tension_diff = SurfaceTensionDerivation(0) #0: use analytical model, #1: use numerical model (use only if surfacetension(T) is non-linear)
#print Header+"** Checking diff(surface_tension): "+str(surface_tension_diff)
print ""

Te_Levy2017  = 1E6 #K
Nexc_Derrien = 1E26 #m-3
Ce = 1.5 * k_b * Nexc_Derrien
tau_e_ph_Levy2017 = 1e-12
ThermalSource_diff   = ThermalSourceModels(laser_intensity, thermal_conductivity, Absorptivity, 1, Ce/tau_e_ph_Levy2017, Te_Levy2017, T) #used without TTM for now
print Header+"** Checking diff(ThermalSource): "+str(ThermalSource_diff*1E-9)+" K/nm vs diff(T)/thickness: "+str(T/thickness*1E-9)+" K/nm."
print Header+"AI/kappa [K/nm]: "+str(1E-9*laser_intensity / thermal_conductivity)

print Header+"** Checking heat conductivity (m2/s) for solid silica at 300 K.: "
print Header+"T-dependent diffusivity (m2/s): "+str(Silica_HeatConductivity(300.) / ( Silica_Solid_HeatCapacity(300.) * mass_density)) #m2/s
print Header+"Constant data from Bauerle book: "+str(Silica_Solid_HeatDiffusivity(T))
print ""
print Header+"** Checking heat conductivity (m2/s) for liquid silica at 1300 K.: "
print Header+"T-dependent diffusivity (m2/s): "+str(Silica_HeatConductivity(1300.) / ( Silica_Liquid_HeatCapacity(300.) * mass_density)) #m2/s
print ""
print Header+"** Checking heat conductivity (m2/s) for liquid silica at 2000 K.: "
print Header+"T-dependent diffusivity (m2/s): "+str(Silica_HeatConductivity(2000.) / ( Silica_Liquid_HeatCapacity(300.) * mass_density)) #m2/s
print Header+"Constant data from Bauerle book: "+str(Silica_Liquid_HeatDiffusivity(T)) 

sound_velocity       = SoundVelocity_Levchenko(mass_density, Peclet_Number, T, surface_tension_diff, ThermalSource_diff)
gamma_Levchenko      = InstabilityGrowthRate_Levchenko(kinematic_viscosity, k, omega, thickness, sound_velocity, omega_0, Peclet_Number, diffusivity)

print Header+"Thermo-convective instability growth rate [1/s]: "+str(gamma_Levchenko)

plt.figure()
#plt.loglog(np.divide(2*pi,k),-gamma_Levchenko)
plt.loglog(np.divide(2*pi,k)/laser_wavelength,-gamma_Levchenko)
#plt.xlabel(r"$\Lambda$ (m)")
plt.xlabel(r"$\Lambda/\lambda$")
plt.ylabel(r"$\gamma$ (s$^-1$)")
filename = "T"+str(T)+"K-h"+str(thickness*1E9)+"nm"
plt.title(filename)
plt.savefig(filename+".png")
#plt.show()

