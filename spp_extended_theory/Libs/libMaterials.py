#!/usr/bin/env python
#-*- coding: utf-8 -*-

# Copyright (C) 2013-2023 T. J.-Y. Derrien
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

## @package libMaterials 
# Functions describing materials and their interaction with light. 

import cmath
from scipy.constants import Boltzmann, m_e, epsilon_0, e, pi, c
import spp_extended_theory.Libs.libLaser as libLaser
import numpy as np
import sys

k_b     = Boltzmann

# OPTICAL FUNCTIONS

## Return the value of dielectric function based on simplified Drude model
# Input:
# @param wavelength (float)
# @param ne (float)
# @param epsilon (complex): dielectric permittivity under wavelength, without 
# excitation
# @param nu (float): collision frequency
# Output: complex-valued dielectric permittivity
def Drude(wavelength, ne, epsilon, nu, meff=1.0):#{{{
  # from scipy.constants import m_e, epsilon_0, e
  omegap2 = ne * e**2 / (m_e * meff * epsilon_0)
  omega=2.0*pi*c/wavelength
  return epsilon - omegap2/(omega*omega) * 1e0/(1e0+1e0j*nu/omega)
#}}}

## Provides a collision frequency estimatead from electron (Te) and lattice (Ti) temperature for Cr.
# Model was taken from Sci. Rep. 7, 8485 (2017)
def CollisionFrequencyModelLevy_Cr(Te, Ti = 300):
    #AtomicDensity = 5.7E28
    A = 4.36E6; B = 9.10E13
    nu_e    = A*Te**2 + B*Ti
    E_Fermi = 8.*e
    v_Fermi = np.sqrt(2. * E_Fermi / m_e) #just converted Fermi energy to velocity
    #nu_c    = (4.*np.pi*AtomicDensity / 3.)**(1./3.)*np.sqrt(v_Fermi ** 2 + k_b * Te / m_e)
    nu_eff  = nu_e #np.min([nu_e, nu_c])
    return nu_eff

## Provides a collision frequency estimatead from electron (Te) and lattice (Ti) temperature for Cr.
# Model was taken from Sci. Rep. 7, 8485 (2017)
def CollisionFrequencyModelLevy_Ti(Te, Ti = 300):
    AtomicDensity = 5.7E28
    A = 2.2E6; B = 3.2E13 #NOTE: parameters for Ti
    nu_e    = A*Te**2 + B*Ti
    E_Fermi = 8.84*e #NOTE: parameter for Ti
    v_Fermi = np.sqrt(2. * E_Fermi / m_e) #just converted Fermi energy to velocity
    nu_c    = (4.*np.pi*AtomicDensity / 3.)**(1./3.)*np.sqrt(v_Fermi ** 2 + k_b * Te / m_e)
    nu_eff  = np.min([nu_e, nu_c])
    return nu_eff

## Return the value of dielectric function based on metallic Drude model
# Input:
# @param wavelength (m) (float)
# @param Te (K) (float)
# @param epsilon (complex): dielectric permittivity under wavelength, without 
# excitation
# Output: complex-valued dielectric permittivity
def Drude_Cr(wavelength, epsilon, Te):#{{{
  Header="[libMaterials] Drude_metal: "
  nu = CollisionFrequencyModelLevy_Cr(Te)
  omega=2.0*pi*c/wavelength
  #omegap2 = omega**2 * (1.-epsilon.real+epsilon.imag**2/(1.-epsilon.real)) #expression from Sci. Rep. #Valid
  omegap2 = omega**2 * (1.-2.*epsilon.real+epsilon.real**2+epsilon.imag**2)/(1.-epsilon.real)
  #print(Header+str(np.sqrt(omegap2)))
  #print(Header+str(nu))
  return 1. - omegap2/(omega*omega) * 1e0/(1e0+1e0j*nu/omega)
#}}}

## Return the value of dielectric function based on metallic Drude model
# Input:
# @param wavelength (m) (float)
# @param Te (K) (float)
# @param epsilon (complex): dielectric permittivity under wavelength, without 
# excitation
# Output: complex-valued dielectric permittivity
def Drude_Ti(wavelength, epsilon, Te):#{{{
  Header="[libMaterials] Drude_metal: "
  nu = CollisionFrequencyModelLevy_Ti(Te)
  omega=2.0*pi*c/wavelength
  #omegap2 = omega**2 * (1.-epsilon.real+epsilon.imag**2/(1.-epsilon.real)) #expression from Sci. Rep. #Valid
  omegap2 = omega**2 * (1.-2.*epsilon.real+epsilon.real**2+epsilon.imag**2)/(1.-epsilon.real)
  #print(Header+str(np.sqrt(omegap2)))
  #print(Header+str(nu))
  return 1. - omegap2/(omega*omega) * 1e0/(1e0+1e0j*nu/omega)
#}}}

## Fresnel reflectivity formula at single interface
# Input:
#   @param eps1: complex-valued permittivity 1+j0
#   @param eps2: idem, for medium2
#   @param angle: angle of incidence (deg or rad?)
# Output: 
#   interface reflectivity (float) R
def reflectivity(eps1, eps2, angle=0, pola="S"):#{{{
  #R=abs(((eps1**0.5e0-eps2**0.5e0)/(eps1**0.5e0+eps2**0.5e0))**2)
  if(pola=="S"): 
      term1 = eps1**0.5E0*np.cos(angle)
      term2 = (eps2-eps1*np.sin(angle)**2)**0.5E0
      R=abs(((term1-term2)/(term1+term2))**2)
  elif(pola=="P"):
      term1 = eps2*np.cos(angle)
      term2 = eps1**0.5E0*(eps2-eps1*np.sin(angle)**2)**0.5E0
      R=abs(((term1-term2)/(term1+term2))**2)
  else: 
      print("Error. Choose pola=S or P, nothing else.")
      sys.exit()
  return R
#}}}

## Complex transmission formula at single interface (normal incidence case)
# Input:
#   @param eps1: complex-valued permittivity 1+j0
#   @param eps2: idem, for medium2
# Output: 
#   interface complex transmission (complex)
def ComplexTransmission(eps1, eps2):
    return 2*np.sqrt(eps1)/(np.sqrt(eps2)+np.sqrt(eps1))

## Computes complex valued reflectivity with normal angle of incidence. Useful 
# for computing further media
# Input:
#   @param eps1: complex-valued permittivity 1+j0
#   @param eps2: idem, for medium2
# Output: 
#   Complex reflectivity. Must be |z|**2 to have the right reflectivity. 
def ComplexReflectivity(eps1, eps2): #{{{
  r = ((eps1**0.5e0-eps2**0.5e0)/(eps1**0.5e0+eps2**0.5e0))
  return r
#}}}

## Returns the complex refractive index
def EpsilonToIndex(eps):
  return cmath.sqrt(eps)

## Returns the complex permittivity from optical index
def IndexToEpsilon(n):
  return n*n

## Maxwell-Garnett method for mixing 2 materials together. 
# Returns the effective dielectric permittivity of 2-mixed materials. 
# @param eps1: dielectric permittivity (epsilon <complex>) of first medium
# @param eps2: dielectric permittivity (epsilon <complex>) of second medium
# @param fraction: fraction refers to epsilon2 which is mixed with (1.-fraction)*epsilon1 medium
# This function was validated by exact comparison with Sergei Lisunov. 
# Applicable for dielectric - metal mixtures. Maybe not applicable for 
# metal-metal mixtures. 
def MaxwellGarnett2(eps1, eps2, fraction):
  eps1r = eps1.real; eps1c = eps1.imag
  eps2r = eps2.real; eps2c = eps2.imag
  # Space for optimization is not so big: 
  # 23*storage+23*assignments+50*multiplications+37*additions+divisions. 
  # Whereas direct writing uses: 
  # 37*additions+104*multiplications+divisions+assignments
  epsilon_effective_real = \
    -(1.*(2.*eps1c**2*fraction**2-2.*fraction*eps2r**2-4.*eps1r+eps1r*fraction*eps2c
    **2-1.*eps1c**2*fraction*eps2r+4.*eps1r*fraction-1.*eps1r**2*fraction*eps2r+
    eps1r*fraction*eps2r**2+4.*eps1r*fraction*eps2r-1.*eps1r*eps2c**2-2.*eps1c**2*
    fraction+4.*eps1c*eps2c*fraction+2.*eps2c**2*fraction**2-2.*eps1r**2*fraction-1.
    *eps1r*eps2r**2-4.*eps1r*eps2r+2.*eps1r**2*fraction**2-4.*fraction*eps2r-2.*
    fraction*eps2c**2-4.*eps2c*fraction**2*eps1c+2.*fraction**2*eps2r**2-4.*fraction
    **2*eps2r*eps1r))/(2.*eps1r*fraction*eps2r+4.+4.*eps1r*fraction-2.*fraction*
    eps2r**2-4.*fraction*eps2r-2.*fraction*eps2c**2+eps2r**2+4.*eps2r+eps2c**2+2.*
    eps1c*eps2c*fraction+eps1r**2*fraction**2+fraction**2*eps2r**2-2.*fraction**2*
    eps2r*eps1r-2.*eps2c*fraction**2*eps1c+eps2c**2*fraction**2+eps1c**2*fraction**2
    )
  # This version uses: 26*additions+68*multiplications+divisions+assignments
  # Optimization would use: 
    #10*storage+10*assignments+46*multiplications+26*additions+divisions
  epsilon_effective_imag = \
    (4.*eps2c*fraction-1.*eps1c*fraction*eps2c**2+4.*eps1c*eps2r-4.*eps1c*fraction*
    eps2r+4.*eps2c*fraction*eps1r+eps2c*fraction*eps1c**2-4.*eps1c*fraction+eps2c*
    fraction*eps1r**2+4.*eps1c+eps1c*eps2c**2+eps1c*eps2r**2-1.*eps1c*fraction*eps2r
    **2)/(2.*eps1r*fraction*eps2r+4.+4.*eps1r*fraction-2.*fraction*eps2r**2-4.*
    fraction*eps2r-2.*fraction*eps2c**2+eps2r**2+4.*eps2r+eps2c**2+2.*eps1c*eps2c*
    fraction+eps1r**2*fraction**2+fraction**2*eps2r**2-2.*fraction**2*eps2r*eps1r-2.
    *eps2c*fraction**2*eps1c+eps2c**2*fraction**2+eps1c**2*fraction**2)
  
  epsilon_effective = epsilon_effective_real + 1.0j*epsilon_effective_imag
  return epsilon_effective

## Maxwell-Garnett method for mixing 3 materials together. 
# Returns the effective dielectric permittivity of 3-mixed materials. 
# Here we keep freedom in definition of fraction1-3 independent.
# It helps for plotting data more easily
# @param eps1: dielectric permittivity (epsilon <complex>) of first medium
# @param eps2: dielectric permittivity (epsilon <complex>) of second medium
# @param eps3: dielectric permittivity (epsilon <complex>) of third medium
# @param fraction1: fraction of epsilon1
# @param fraction2: fraction of epsilon2
# @param fraction3: fraction of epsilon3
# Applicable for dielectric - metal mixtures. Maybe not applicable for 
# metal-metal mixtures. 
def MaxwellGarnett3(eps1, eps2, eps3, fraction1, fraction2, fraction3):
  eps1r = eps1.real; eps1c = eps1.imag
  eps2r = eps2.real; eps2c = eps2.imag
  eps3r = eps3.real; eps3c = eps3.imag
  epsilon_effective_real = \
    -(-64+8*eps2c**2*fraction2**2*eps3c**2+32*fraction1+32*fraction2+32*fraction3-64
    *eps1r-64*eps2r-64*eps3r-16*eps1c**2-16*eps2c**2-16*eps3c**2-16*fraction3**2*
    eps1c**2*eps3r-16*fraction3**2*eps2r**2*eps3r+36*eps1c*fraction1*eps3c*fraction3
    *eps2r**2+8*eps2c**2*fraction2**2*eps3r**2+8*eps3c**2*fraction3**2*eps2r**2+8*
    fraction2**2*eps1c**2*eps2r**2-32*fraction2*eps2r*fraction1+16*fraction3*eps1r**
    2*fraction1*eps3r**2-32*fraction3*eps1r*eps2r**2*fraction2-32*fraction3*eps1r*
    eps2c**2*fraction2-8*fraction3*eps1c**2*eps2r**2*fraction2-16*fraction3**2*eps1c
    **2*eps2r*eps3r-4*fraction3**2*eps1c**2*eps2c**2*eps3r-8*fraction1*eps2r**2*
    eps3r**2*fraction2-8*fraction1*eps2r**2*eps3r**2*fraction3-32*fraction1*eps2r**2
    *eps3r*fraction2+16*fraction3*eps2r**2*eps3r*fraction2+8*eps3c**2*fraction3**2*
    eps1r*eps2r**2+32*eps2c**2*fraction2**2*eps1r*eps3r+8*eps2c**2*fraction2**2*
    eps1r**2*eps3r+32*eps3c**2*fraction3**2*eps1r*eps2r+16*fraction2*fraction1*eps1c
    **2*eps2r**2+16*fraction2*fraction1*eps1c**2*eps2c**2-8*fraction2*eps1c**2*eps2r
    *fraction3-32*fraction2*eps1r*eps3r**2*fraction3-32*fraction2*eps1r*eps3r*
    fraction3-32*fraction2*eps1r*eps3c**2*fraction3+2*eps2c**2*fraction2**2*eps1r**2
    *eps3c**2+16*fraction2*eps1c**2*eps2r*fraction1-8*fraction2*eps1r*eps3r**2*
    fraction1-32*fraction2*eps1r*eps3r*fraction1+8*eps3c**2*fraction3**2*eps1c**2+16
    *fraction3*eps2c**2*fraction2*eps3r**2-8*fraction2*eps1c**2*eps3r**2*fraction1-
    32*fraction2*eps1c**2*eps3r*fraction1-8*fraction2*eps1c**2*eps3c**2*fraction1+16
    *fraction3*eps2r**2*fraction2*eps3c**2+16*fraction3*eps2r**2*fraction2*eps3r**2-
    8*fraction1*eps2r**2*eps3c**2*fraction2-8*fraction1*eps2r**2*eps3c**2*fraction3-
    8*fraction1*eps2c**2*eps3r**2*fraction2-8*fraction1*eps2c**2*eps3r**2*fraction3-
    32*fraction1*eps2c**2*eps3r*fraction2-8*fraction1*eps2c**2*eps3c**2*fraction2-8*
    fraction1*eps2c**2*eps3c**2*fraction3-8*fraction3*eps1c**2*eps2c**2*fraction2-8*
    fraction2*eps1r**2*eps3r**2*fraction1-32*fraction2*eps1r**2*eps3r*fraction1-8*
    fraction2*eps1r**2*eps3c**2*fraction1+2*eps1c**2*fraction1**2*eps2r**2*eps3r**2+
    2*eps2c**2*fraction2**2*eps1r**2*eps3r**2+2*eps3c**2*fraction3**2*eps1r**2*eps2r
    **2+16*fraction3*fraction1*eps1c**2*eps3r**2+16*fraction3*fraction1*eps1c**2*
    eps3c**2+32*fraction3**2*eps1r*eps2r*eps3r**2+8*fraction3**2*eps1c**2*eps2r*
    eps3r**2+2*fraction3**2*eps1c**2*eps2c**2*eps3r**2-8*fraction3*eps1r**2*eps2r**2
    *fraction2-8*fraction3*eps1r**2*eps2c**2*fraction2-4*fraction3**2*eps1c**2*eps2r
    **2*eps3r-4*fraction3**2*eps1r**2*eps2r**2*eps3r-16*fraction3**2*eps1r**2*eps2r*
    eps3r+16*fraction3*eps1r**2*fraction1*eps3c**2+8*fraction2**2*eps3c**2*eps2r**2+
    8*fraction1**2*eps1r**2*eps2r*eps3c**2-32*fraction3*eps1r**2*eps2r*fraction1-8*
    fraction3*eps1r**2*eps2c**2*fraction1-8*fraction3*eps1r*eps2r**2*fraction1-32*
    fraction3*eps1r*eps2r*fraction1+144*eps1c*fraction1*eps3c*fraction3+144*eps2c*
    fraction2*eps3c*fraction3-8*fraction3*eps1r**2*eps2r**2*fraction1+144*eps1c*
    fraction1*eps2c*fraction2-4*fraction3**2*eps1r**2*eps2c**2*eps3r+8*fraction3**2*
    eps1r*eps2r**2*eps3r**2+16*fraction3*eps2c**2*fraction2*eps3c**2+32*eps1c**2*
    fraction1**2+32*fraction3**2*eps1r+32*fraction3**2*eps2r-64*fraction3**2*eps3r+8
    *fraction3**2*eps1r**2+8*fraction3**2*eps2r**2+8*fraction3**2*eps1c**2+8*
    fraction3**2*eps2c**2+32*fraction2**2*eps1r-64*fraction2**2*eps2r+32*fraction2**
    2*eps3r+8*fraction2**2*eps1r**2+8*fraction2**2*eps3r**2+8*fraction2**2*eps1c**2+
    8*fraction2**2*eps3c**2+32*eps3c**2*fraction3**2-4*eps1r**2*eps2r**2-16*eps1r**2
    *eps2r-4*eps1r**2*eps2c**2-16*eps1r*eps2r**2-64*eps1r*eps2r-16*eps1r*eps2c**2-4*
    eps1c**2*eps2r**2-16*eps1c**2*eps2r-16*eps1r**2-16*eps2r**2-16*eps3r**2+32*
    fraction1**2+32*fraction2**2+32*fraction3**2-fraction2*eps1r**2*eps2r**2*eps3c**
    2-4*fraction2*eps1r*eps2r**2*eps3r**2-fraction2*eps1r**2*eps2c**2*eps3c**2-4*
    fraction2*eps1r*eps2r*eps3c**2-fraction2*eps1r**2*eps2c**2*eps3r**2-16*fraction2
    *eps1r*eps2r*eps3r-fraction2*eps1r**2*eps2r**2*eps3r**2-16*fraction2*eps1r*eps2r
    **2*eps3r-fraction2*eps1c**2*eps2r**2*eps3r**2-16*fraction2*eps1r*eps2c**2*eps3r
    -fraction2*eps1r**2*eps2r*eps3r**2-4*fraction2*eps1r*eps2r**2*eps3c**2-4*
    fraction2*eps1c**2*eps2r**2*eps3r-4*fraction2*eps1r*eps2r*eps3r**2-4*fraction2*
    eps1r**2*eps2r**2*eps3r-4*fraction2*eps1r*eps2c**2*eps3c**2-fraction2*eps1c**2*
    eps2r**2*eps3c**2-fraction2*eps1c**2*eps2r*eps3r**2-4*fraction2*eps1c**2*eps2r*
    eps3r-fraction2*eps1c**2*eps2r*eps3c**2-fraction2*eps1c**2*eps2c**2*eps3r**2-4*
    fraction2*eps1c**2*eps2c**2*eps3r-fraction2*eps1c**2*eps2c**2*eps3c**2-4*
    fraction3*eps1r**2*eps2r*eps3r-4*fraction3*eps1r*eps2c**2*eps3r**2-4*fraction3*
    eps1r**2*eps2r*eps3c**2-fraction3*eps1r**2*eps2c**2*eps3r-fraction3*eps1r**2*
    eps2r**2*eps3c**2-4*fraction3*eps1r*eps2r**2*eps3r**2-16*fraction3*eps1r*eps2r*
    eps3c**2+36*eps1c*fraction1*eps3c*fraction3*eps2c**2+144*eps2c*fraction2*eps3c*
    fraction3*eps1r+36*eps2c*fraction2*eps3c*fraction3*eps1r**2+36*eps2c*fraction2*
    eps3c*fraction3*eps1c**2+4*fraction2*eps1r**2*eps2r*eps3r**2*fraction3+144*eps1c
    *fraction1*eps2c*fraction2*eps3r+36*eps1c*fraction1*eps2c*fraction2*eps3r**2+36*
    eps1c*fraction1*eps2c*fraction2*eps3c**2+144*eps1c*fraction1*eps3c*fraction3*
    eps2r+4*fraction2*eps1r**2*eps2r**2*eps3r**2*fraction1+16*fraction2*eps1r*eps2r*
    eps3c**2*fraction3-fraction3*eps1r**2*eps2c**2*eps3r**2-16*fraction3*eps1r*eps2r
    *eps3r-fraction3*eps1r**2*eps2r**2*eps3r**2-4*fraction3*eps1r*eps2r**2*eps3r-
    fraction3*eps1c**2*eps2r**2*eps3r**2-4*fraction3*eps1r*eps2c**2*eps3r-4*
    fraction3*eps1r**2*eps2r*eps3r**2-4*fraction3*eps1r*eps2r**2*eps3c**2-fraction3*
    eps1c**2*eps2r**2*eps3r-16*fraction3*eps1r*eps2r*eps3r**2-fraction3*eps1r**2*
    eps2r**2*eps3r-4*fraction3*eps1r*eps2c**2*eps3c**2-fraction3*eps1c**2*eps2r**2*
    eps3c**2-4*fraction3*eps1c**2*eps2r*eps3r**2-4*fraction3*eps1c**2*eps2r*eps3r-4*
    fraction3*eps1c**2*eps2r*eps3c**2-fraction3*eps1c**2*eps2c**2*eps3r**2-fraction3
    *eps1c**2*eps2c**2*eps3r-fraction3*eps1c**2*eps2c**2*eps3c**2+2*fraction2**2*
    eps1r**2*eps2r**2*eps3c**2+2*fraction2**2*eps1r**2*eps2r**2*eps3r**2+8*fraction2
    **2*eps1r*eps2r**2*eps3r**2+2*fraction2**2*eps1c**2*eps2r**2*eps3r**2+32*
    fraction2**2*eps1r*eps2r**2*eps3r-4*fraction2**2*eps1r**2*eps2r*eps3r**2-16*
    fraction2**2*eps1r*eps2r*eps3r**2+8*eps3c**2*fraction3**2*eps1r**2*eps2r+2*eps2c
    **2*fraction2**2*eps1c**2*eps3c**2+8*eps3c**2*fraction3**2*eps1r*eps2c**2-16*
    fraction2**2*eps1r*eps2r*eps3c**2-4*eps1c**2*eps2c**2-16*eps1r*eps3r**2-64*eps1r
    *eps3r-16*eps1r*eps3c**2-16*eps2r*eps3r**2-64*eps2r*eps3r-16*eps2r*eps3c**2-4*
    eps1r**2*eps3r**2-16*eps1r**2*eps3r-4*eps1r**2*eps3c**2-4*eps2r**2*eps3r**2-16*
    eps2r**2*eps3r-4*eps2r**2*eps3c**2-4*eps1c**2*eps3r**2-16*eps1c**2*eps3r-4*eps1c
    **2*eps3c**2-4*eps2c**2*eps3r**2-16*eps2c**2*eps3r-4*eps2c**2*eps3c**2-16*
    fraction1*eps1r+32*fraction1*eps2r+32*fraction1*eps3r-16*fraction1*eps1r**2+8*
    fraction1*eps2r**2+8*fraction1*eps3r**2-16*fraction1*eps1c**2+8*fraction1*eps2c**2
    +8*fraction1*eps3c**2+32*fraction2*eps1r-16*fraction2*eps2r+32*fraction2*eps3r
    +8*fraction2*eps1r**2-16*fraction2*eps2r**2+8*fraction2*eps3r**2+8*fraction2*
    eps1c**2-16*fraction2*eps2c**2+8*fraction2*eps3c**2+32*fraction3*eps1r+32*
    fraction3*eps2r-16*fraction3*eps3r+8*fraction3*eps1r**2+8*fraction3*eps2r**2-16*
    fraction3*eps3r**2+8*fraction3*eps1c**2+8*fraction3*eps2c**2-16*fraction3*eps3c**2
    +32*fraction1**2*eps1r**2+32*fraction2**2*eps2r**2+32*eps2c**2*fraction2**2+64
    *fraction2*fraction1+64*fraction3*fraction1+64*fraction2*fraction3+32*fraction3**2
    *eps3r**2-64*fraction2**2*eps1r*eps2r*eps3r+8*fraction2**2*eps1r**2*eps2r**2*
    eps3r+8*fraction2**2*eps2r**2*eps1r*eps3c**2+2*eps3c**2*fraction3**2*eps1c**2*
    eps2c**2+2*eps1c**2*fraction1**2*eps2c**2*eps3r**2+2*fraction3**2*eps1c**2*eps2r
    **2*eps3r**2+2*eps3c**2*fraction3**2*eps1r**2*eps2c**2+8*eps2c**2*fraction2**2*
    eps1r*eps3r**2+8*eps1c**2*fraction1**2*eps2r*eps3c**2+2*eps1c**2*fraction1**2*
    eps2r**2*eps3c**2+8*fraction2**2*eps1c**2*eps2r**2*eps3r+2*fraction2**2*eps1c**2
    *eps2r**2*eps3c**2+8*fraction3**2*eps1r*eps2c**2*eps3r**2+8*eps2c**2*fraction2**
    2*eps1c**2*eps3r+16*fraction1*eps1r*eps2r**2*fraction2+16*fraction1*eps1r*eps2c**2
    *fraction2+8*fraction3**2*eps3r**2*eps1r**2*eps2r+2*fraction3**2*eps1r**2*
    eps2c**2*eps3r**2+2*fraction3**2*eps2r**2*eps3r**2*eps1r**2-8*fraction3*eps1r*
    eps2c**2*fraction1+8*fraction1**2*eps1r**2*eps2r**2*eps3r+2*fraction1**2*eps1r**2
    *eps2r**2*eps3r**2+2*fraction1**2*eps1r**2*eps2r**2*eps3c**2+32*fraction1**2*
    eps1r**2*eps2r*eps3r+8*fraction1**2*eps1r**2*eps2r*eps3r**2+2*fraction1**2*eps1r
    **2*eps2c**2*eps3c**2+8*fraction1**2*eps1r**2*eps2c**2*eps3r+2*fraction1**2*
    eps1r**2*eps2c**2*eps3r**2+16*fraction2*eps1r*eps2r*eps3r*fraction3+4*fraction2*
    eps1r**2*eps2r**2*eps3c**2*fraction1+16*fraction2*eps1r**2*eps2r**2*eps3r*
    fraction1+4*fraction2*fraction3*eps1r**2*eps2c**2*eps3c**2+4*fraction2*eps1r**2*
    eps2r*eps3c**2*fraction3+4*fraction3*eps1r**2*eps2r**2*eps3r*fraction1+16*
    fraction3*eps1c**2*eps2r*eps3r*fraction1+4*fraction3*eps1c**2*eps2c**2*eps3r*
    fraction1+4*fraction3*eps1c**2*eps2r**2*eps3r*fraction1+4*fraction2*eps1r**2*
    eps2r*eps3r*fraction3+4*fraction3*eps1r*eps2c**2*eps3r*fraction1+4*fraction3*
    eps1r*eps2r**2*eps3r*fraction1+16*fraction3*eps1r*eps2r*eps3r*fraction1+4*
    fraction3*eps1r**2*eps2r**2*eps3r*fraction2+4*fraction3*eps1c**2*eps2c**2*eps3r*
    fraction2+4*fraction3*eps1r**2*eps2c**2*eps3r*fraction2+16*fraction1*eps1r*eps2r
    *eps3r**2*fraction3+4*fraction1*eps1r*eps2r**2*eps3c**2*fraction2+4*fraction1*
    eps1r*eps2r**2*eps3c**2*fraction3+16*fraction1*eps1r*eps2c**2*eps3r*fraction2+16
    *fraction1*eps1r*eps2r**2*eps3r*fraction2+16*fraction1*eps1r*eps2r*eps3c**2*
    fraction3+4*fraction1*eps1r*eps2r**2*eps3r**2*fraction2+4*fraction1*eps1r*eps2r**2
    *eps3r**2*fraction3+4*fraction2*eps3r**2*fraction1*eps1r**2*eps2c**2-16*
    fraction1**2*eps1r*eps2r**2-64*fraction1**2*eps1r*eps2r-16*eps1r**2*eps2r*eps3r-
    16*fraction1**2*eps1r*eps2c**2-16*fraction1**2*eps1r*eps3r**2-64*fraction1**2*
    eps1r*eps3r-16*fraction1**2*eps1r*eps3c**2+16*fraction2*eps3r**2*fraction1+32*
    eps2c**2*fraction2**2*eps1r+16*fraction3*eps2r**2*fraction1-32*fraction2*eps1r**
    2*fraction1+8*eps2c**2*fraction2**2*eps1c**2-32*fraction3*eps1r**2*fraction1+8*
    eps1c**2*fraction1**2*eps3r**2-4*eps1r*eps2c**2*eps3r**2-4*eps1r**2*eps2r*eps3c**2
    -4*eps1r**2*eps2c**2*eps3r-eps1r**2*eps2r**2*eps3c**2-4*eps1r*eps2r**2*eps3r**
    2-eps1r**2*eps2c**2*eps3c**2-16*eps1r*eps2r*eps3c**2-eps1r**2*eps2c**2*eps3r**2-
    64*eps1r*eps2r*eps3r-eps1r**2*eps2r**2*eps3r**2-16*eps1r*eps2r**2*eps3r-eps1c**2
    *eps2r**2*eps3r**2-16*eps1r*eps2c**2*eps3r-4*eps1r**2*eps2r*eps3r**2-4*eps1r*
    eps2r**2*eps3c**2-4*eps1c**2*eps2r**2*eps3r-16*eps1r*eps2r*eps3r**2-4*eps1r**2*
    eps2r**2*eps3r-4*eps1r*eps2c**2*eps3c**2-eps1c**2*eps2r**2*eps3c**2-4*eps1c**2*
    eps2r*eps3r**2-16*eps1c**2*eps2r*eps3r-4*eps1c**2*eps2r*eps3c**2-eps1c**2*eps2c**2
    *eps3r**2-4*eps1c**2*eps2c**2*eps3r-eps1c**2*eps2c**2*eps3c**2-4*fraction1*
    eps1r**2*eps2r**2-16*fraction1*eps1r**2*eps2r-4*fraction1*eps1r**2*eps2c**2-4*
    fraction1*eps1r*eps2r**2-16*fraction1*eps1r*eps2r-4*fraction1*eps1r*eps2c**2-4*
    fraction1*eps1c**2*eps2r**2-16*fraction1*eps1c**2*eps2r-4*fraction1*eps1c**2*
    eps2c**2-4*fraction1*eps1r*eps3r**2+4*fraction2*eps3r**2*fraction1*eps1c**2*
    eps2r**2+4*fraction2*eps3r**2*fraction1*eps1c**2*eps2c**2+4*fraction2*eps3c**2*
    fraction1*eps1r**2*eps2c**2+4*fraction2*eps3c**2*fraction1*eps1c**2*eps2r**2+4*
    fraction2*eps3c**2*fraction1*eps1c**2*eps2c**2+4*fraction1*eps1r*eps2c**2*eps3r**2
    *fraction2+4*fraction1*eps1r*eps2c**2*eps3r**2*fraction3+16*fraction2*eps3r*
    fraction1*eps1r**2*eps2c**2+16*fraction2*eps3r*fraction1*eps1c**2*eps2r**2+16*
    fraction2*eps3r*fraction1*eps1c**2*eps2c**2+4*fraction1*eps1r*eps2c**2*eps3c**2*
    fraction2+4*fraction1*eps1r*eps2c**2*eps3c**2*fraction3+4*fraction3*eps1r**2*
    eps2c**2*eps3r*fraction1+16*fraction3*eps1r**2*eps2r*eps3r*fraction1+4*fraction2
    *eps1c**2*eps2r*eps3c**2*fraction1+4*fraction2*eps1c**2*eps2r*eps3r**2*fraction1
    +16*fraction2*eps1c**2*eps2r*eps3r*fraction1+4*fraction2*eps1r*eps2r*eps3r**2*
    fraction1+4*fraction2*eps1r**2*eps2r*eps3r**2*fraction1+4*fraction2*eps1r*eps2r*
    eps3c**2*fraction1+16*fraction2*eps1r*eps2r*eps3r*fraction1+4*fraction2*eps1r**2
    *eps2r*eps3c**2*fraction1+16*fraction2*eps1r**2*eps2r*eps3r*fraction1+16*
    fraction3*fraction2*eps1r*eps2c**2*eps3r**2+4*fraction3*fraction2*eps1r**2*eps2r
    **2*eps3c**2-8*fraction2*eps2r*eps3c**2*fraction1-64*fraction1**2*eps1r+32*
    fraction1**2*eps2r+32*fraction1**2*eps3r+8*fraction1**2*eps2r**2+8*fraction1**2*
    eps3r**2+8*fraction1**2*eps2c**2+8*fraction1**2*eps3c**2+8*fraction2**2*eps3r**2
    *eps2r**2-16*fraction1*eps1r*eps3r-4*fraction1*eps1r*eps3c**2+8*fraction1*eps2r*
    eps3r**2+32*fraction1*eps2r*eps3r+8*fraction1*eps2r*eps3c**2-4*fraction1*eps1r**
    2*eps3r**2-16*fraction1*eps1r**2*eps3r-4*fraction1*eps1r**2*eps3c**2+2*fraction1
    *eps2r**2*eps3r**2+8*fraction1*eps2r**2*eps3r+2*fraction1*eps2r**2*eps3c**2-4*
    fraction1*eps1c**2*eps3r**2-16*fraction1*eps1c**2*eps3r-4*fraction1*eps1c**2*
    eps3c**2+2*fraction1*eps2c**2*eps3r**2+8*fraction1*eps2c**2*eps3r+2*fraction1*
    eps2c**2*eps3c**2-4*fraction2*eps1r**2*eps2r**2-4*fraction2*eps1r**2*eps2r-4*
    fraction2*eps1r**2*eps2c**2-16*fraction2*eps1r*eps2r**2-16*fraction2*eps1r*eps2r
    -16*fraction2*eps1r*eps2c**2-4*fraction2*eps1c**2*eps2r**2-4*fraction2*eps1c**2*
    eps2r-4*fraction2*eps1c**2*eps2c**2+8*fraction2*eps1r*eps3r**2+32*fraction2*
    eps1r*eps3r+8*fraction2*eps1r*eps3c**2-4*fraction2*eps2r*eps3r**2-16*fraction2*
    eps2r*eps3r-4*fraction2*eps2r*eps3c**2+2*fraction2*eps1r**2*eps3r**2+8*fraction2
    *eps1r**2*eps3r+2*fraction2*eps1r**2*eps3c**2-4*fraction2*eps2r**2*eps3r**2-16*
    fraction2*eps2r**2*eps3r-4*fraction2*eps2r**2*eps3c**2+2*fraction2*eps1c**2*
    eps3r**2+8*fraction2*eps1c**2*eps3r+2*fraction2*eps1c**2*eps3c**2-4*fraction2*
    eps2c**2*eps3r**2-16*fraction2*eps2c**2*eps3r-4*fraction2*eps2c**2*eps3c**2+2*
    fraction3*eps1r**2*eps2r**2+8*fraction3*eps1r**2*eps2r+2*fraction3*eps1r**2*
    eps2c**2+8*fraction3*eps1r*eps2r**2+32*fraction3*eps1r*eps2r+8*fraction3*eps1r*
    eps2c**2+2*fraction3*eps1c**2*eps2r**2+8*fraction3*eps1c**2*eps2r+2*fraction3*
    eps1c**2*eps2c**2-16*fraction3*eps1r*eps3r**2-16*fraction3*eps1r*eps3r-16*
    fraction3*eps1r*eps3c**2-16*fraction3*eps2r*eps3r**2-16*fraction3*eps2r*eps3r-16
    *fraction3*eps2r*eps3c**2-4*fraction3*eps1r**2*eps3r**2-4*fraction3*eps1r**2*
    eps3r-4*fraction3*eps1r**2*eps3c**2-4*fraction3*eps2r**2*eps3r**2-4*fraction3*
    eps2r**2*eps3r-4*fraction3*eps2r**2*eps3c**2-4*fraction3*eps1c**2*eps3r**2-4*
    fraction3*eps1c**2*eps3r-4*fraction3*eps1c**2*eps3c**2-4*fraction3*eps2c**2*
    eps3r**2-4*fraction3*eps2c**2*eps3r-4*fraction3*eps2c**2*eps3c**2+8*fraction3**2
    *eps2c**2*eps3r**2+8*fraction3**2*eps1r**2*eps3r**2+8*fraction3**2*eps2r**2*
    eps3r**2+8*fraction3**2*eps1c**2*eps3r**2+32*fraction3**2*eps1r*eps3r**2+32*
    fraction3**2*eps2r*eps3r**2-32*fraction3*eps3r*fraction1+16*fraction3*fraction2*
    eps1r*eps2r**2*eps3r**2+4*fraction3*fraction2*eps1r**2*eps2c**2*eps3r**2+4*
    fraction3*fraction2*eps1r**2*eps2r**2*eps3r**2+4*fraction3*fraction2*eps1c**2*
    eps2r**2*eps3r**2+16*fraction3*eps2r*fraction1*eps1c**2*eps3c**2+4*fraction3*
    eps2r**2*fraction1*eps1c**2*eps3c**2+4*fraction3*eps1c**2*fraction1*eps2c**2*
    eps3r**2+4*fraction3*eps1c**2*fraction1*eps2c**2*eps3c**2+16*fraction2*eps1r*
    eps2r*eps3r**2*fraction3+4*fraction2*eps1c**2*eps2r*eps3r**2*fraction3+4*
    fraction2*eps1c**2*eps2r*eps3r*fraction3+4*fraction2*eps1c**2*eps2r*eps3c**2*
    fraction3+16*fraction3*eps1r**2*fraction1*eps2r*eps3r**2+16*fraction3*eps1c**2*
    fraction1*eps2r*eps3r**2+16*fraction3*fraction2*eps1r*eps2c**2*eps3c**2+4*
    fraction3*fraction2*eps1c**2*eps2r**2*eps3c**2+4*fraction3*fraction2*eps1c**2*
    eps2c**2*eps3r**2+4*fraction3*fraction2*eps1c**2*eps2c**2*eps3c**2+16*fraction3*
    eps2r**2*fraction2*eps1r*eps3c**2+4*fraction3*eps1r**2*eps2r**2*fraction1*eps3r**2
    +4*fraction3*eps1r**2*eps2r**2*fraction1*eps3c**2+16*fraction3*eps1r**2*eps2r*
    fraction1*eps3c**2+4*fraction3*eps1r**2*eps2c**2*fraction1*eps3r**2+4*fraction3*
    eps1r**2*eps2c**2*fraction1*eps3c**2+4*fraction3*eps1c**2*eps2r**2*fraction1*
    eps3r**2+16*fraction3*eps1r*eps2r**2*eps3r*fraction2+16*fraction3*eps1r*eps2c**2
    *eps3r*fraction2+4*fraction3*eps1c**2*eps2r**2*eps3r*fraction2-32*fraction3*
    eps1r*fraction1-32*fraction2*eps1c**2*fraction1+16*fraction2*eps3c**2*fraction1-
    32*fraction3*eps1c**2*fraction1+16*fraction3*eps2c**2*fraction1-32*fraction2*
    eps1r*fraction1+8*eps2c**2*fraction2**2*eps1r**2+8*eps1c**2*fraction1**2*eps2c**
    2+32*eps3c**2*fraction3**2*eps1r+8*eps1c**2*fraction1**2*eps3c**2+2*fraction1**2
    *eps2c**2*eps3c**2+2*fraction1**2*eps2c**2*eps3r**2+8*fraction1**2*eps2c**2*
    eps3r+2*fraction1**2*eps2r**2*eps3c**2+8*fraction1**2*eps2r**2*eps3r+2*fraction1
    **2*eps2r**2*eps3r**2+32*fraction1**2*eps2r*eps3r+8*fraction1**2*eps2r*eps3c**2+
    8*fraction1**2*eps2r*eps3r**2+32*eps2c**2*fraction2**2*eps3r+32*eps1c**2*
    fraction1**2*eps3r-32*fraction2*eps2r*fraction3-32*fraction2*eps3r*fraction3+16*
    fraction2*eps1r**2*fraction3-32*fraction2*eps3r**2*fraction3+32*fraction1**2*
    eps1r**2*eps2r+8*fraction1**2*eps1r**2*eps2c**2+8*eps3c**2*fraction3**2*eps2c**2
    +8*fraction1**2*eps1r**2*eps2r**2-32*fraction1*eps2r**2*fraction2-32*fraction1*
    eps3r**2*fraction3-32*fraction1*eps2c**2*fraction2-32*fraction1*eps3c**2*
    fraction3+32*eps3c**2*fraction3**2*eps2r+64*fraction2*eps3r*fraction1+8*
    fraction2**2*eps1r**2*eps2r**2+8*fraction3**2*eps1r**2*eps2r+2*fraction3**2*
    eps1r**2*eps2c**2+8*fraction3**2*eps1r*eps2r**2-8*fraction2*eps1r*eps3c**2*
    fraction1-8*fraction2*eps2r*eps3r**2*fraction1-32*fraction2*eps2r*eps3r*
    fraction1-64*fraction3**2*eps1r*eps2r*eps3r-16*fraction3**2*eps1r*eps2r**2*eps3r
    -16*fraction3**2*eps1r*eps2c**2*eps3r-8*fraction2*eps1c**2*eps3r*fraction3-8*
    fraction2*eps1c**2*eps3c**2*fraction3+8*eps1c**2*fraction1**2*eps2r**2*eps3r-8*
    fraction2*eps1r**2*eps3r**2*fraction3-8*fraction2*eps1r**2*eps3r*fraction3-8*
    fraction2*eps1r**2*eps3c**2*fraction3-8*fraction2*eps1c**2*eps3r**2*fraction3+16
    *fraction2*eps1r**2*eps2r*fraction1+16*fraction2*eps1r*eps2r*fraction1+16*
    fraction2*eps2r*eps3r**2*fraction3+16*fraction2*eps2r*eps3r*fraction3+16*
    fraction2*eps2r*eps3c**2*fraction3-16*fraction1**2*eps1r*eps2r*eps3r**2-4*
    fraction1**2*eps1r*eps2c**2*eps3c**2+32*eps1c**2*fraction1**2*eps2r*eps3r-4*
    fraction1**2*eps1r*eps2r**2*eps3c**2-64*fraction1**2*eps1r*eps2r*eps3r-16*
    fraction1**2*eps1r*eps2r**2*eps3r-16*fraction1**2*eps1r*eps2c**2*eps3r-4*
    fraction1**2*eps1r*eps2c**2*eps3r**2-4*fraction1**2*eps1r*eps2r**2*eps3r**2-16*
    fraction1**2*eps1r*eps2r*eps3c**2-8*fraction3*eps2r**2*eps3r*fraction1+16*
    fraction3*eps1c**2*eps3r*fraction1-8*fraction3*eps2c**2*eps3r*fraction1-32*
    fraction3*eps1c**2*eps2r*fraction1-8*fraction3*eps1c**2*eps2c**2*fraction1+16*
    fraction3*eps1r*eps3r*fraction1-32*fraction3*eps2r*eps3r*fraction1+16*fraction3*
    eps1r**2*eps3r*fraction1-8*fraction3*eps1c**2*eps2r**2*fraction1+2*eps2c**2*
    fraction2**2*eps1c**2*eps3r**2-4*fraction2**2*eps1r**2*eps2r*eps3c**2-16*
    fraction2**2*eps1r**2*eps2r*eps3r+8*eps1c**2*fraction1**2*eps2r*eps3r**2-32*
    fraction1*eps2r*eps3r**2*fraction3-32*fraction1*eps2r*eps3c**2*fraction3+16*
    fraction3*eps2c**2*eps3r*fraction2+16*fraction1*eps1r*eps3r**2*fraction3+16*
    fraction1*eps1r*eps3c**2*fraction3+8*eps3c**2*fraction3**2*eps1c**2*eps2r+16*
    fraction2*fraction1*eps1r**2*eps2c**2+16*fraction2*fraction1*eps1r**2*eps2r**2+2
    *eps3c**2*fraction3**2*eps1c**2*eps2r**2+8*eps2c**2*fraction2**2*eps1r*eps3c**2-
    4*fraction2**2*eps1c**2*eps2r*eps3r**2-16*fraction2**2*eps1c**2*eps2r*eps3r-4*
    fraction2**2*eps1c**2*eps2r*eps3c**2+8*eps1c**2*fraction1**2*eps2c**2*eps3r+2*
    eps1c**2*fraction1**2*eps2c**2*eps3c**2-8*fraction2*eps1r**2*eps2r*fraction3-32*
    fraction2*eps1r*eps2r*fraction3-fraction3*eps1r**2*eps2c**2*eps3c**2-16*
    fraction1*eps1r**2*eps2r*eps3r+32*fraction3**2*eps1r*eps2r+8*fraction3**2*eps1r*
    eps2c**2+2*fraction3**2*eps1c**2*eps2r**2+8*fraction3**2*eps1c**2*eps2r+2*
    fraction3**2*eps1c**2*eps2c**2-64*fraction3**2*eps1r*eps3r-64*fraction3**2*eps2r
    *eps3r+2*fraction2**2*eps1r**2*eps3r**2+8*fraction2**2*eps1r**2*eps3r+2*
    fraction2**2*eps1r**2*eps3c**2+2*fraction2**2*eps1c**2*eps3r**2+8*fraction2**2*
    eps1c**2*eps3r+2*fraction2**2*eps1c**2*eps3c**2+2*fraction3**2*eps1r**2*eps2r**2
    -16*fraction2**2*eps1r**2*eps2r-64*fraction2**2*eps1r*eps2r-16*fraction2**2*
    eps1c**2*eps2r+8*fraction2**2*eps1r*eps3r**2+32*fraction2**2*eps1r*eps3r+8*
    fraction2**2*eps1r*eps3c**2-16*fraction2**2*eps2r*eps3r**2-64*fraction2**2*eps2r
    *eps3r-16*fraction2**2*eps2r*eps3c**2+32*eps1c**2*fraction1**2*eps2r-32*
    fraction3*eps2r**2*fraction2-32*fraction3*eps2c**2*fraction2+8*fraction1**2*
    eps1r**2*eps3r**2+8*eps3c**2*fraction3**2*eps1r**2+8*fraction1**2*eps1r**2*eps3c
    **2+32*fraction1**2*eps1r**2*eps3r+32*fraction2**2*eps1r*eps2r**2+8*eps1c**2*
    fraction1**2*eps2r**2+32*fraction2**2*eps3r*eps2r**2+64*fraction2*eps1r*
    fraction3-16*fraction3**2*eps2c**2*eps3r+16*fraction2*eps1c**2*fraction3-32*
    fraction2*eps3c**2*fraction3+64*fraction3*eps2r*fraction1-fraction1*eps1r*eps2c**2
    *eps3r**2-4*fraction1*eps1r**2*eps2r*eps3c**2-4*fraction1*eps1r**2*eps2c**2*
    eps3r-fraction1*eps1r**2*eps2r**2*eps3c**2-fraction1*eps1r*eps2r**2*eps3r**2-
    fraction1*eps1r**2*eps2c**2*eps3c**2-4*fraction1*eps1r*eps2r*eps3c**2-fraction1*
    eps1r**2*eps2c**2*eps3r**2-16*fraction1*eps1r*eps2r*eps3r-fraction1*eps1r**2*
    eps2r**2*eps3r**2-4*fraction1*eps1r*eps2r**2*eps3r-fraction1*eps1c**2*eps2r**2*
    eps3r**2-4*fraction1*eps1r*eps2c**2*eps3r-4*fraction1*eps1r**2*eps2r*eps3r**2-
    fraction1*eps1r*eps2r**2*eps3c**2-4*fraction1*eps1c**2*eps2r**2*eps3r-4*
    fraction1*eps1r*eps2r*eps3r**2-4*fraction1*eps1r**2*eps2r**2*eps3r-fraction1*
    eps1r*eps2c**2*eps3c**2-fraction1*eps1c**2*eps2r**2*eps3c**2-4*fraction1*eps1c**
    2*eps2r*eps3r**2-16*fraction1*eps1c**2*eps2r*eps3r-4*fraction1*eps1c**2*eps2r*
    eps3c**2-fraction1*eps1c**2*eps2c**2*eps3r**2-4*fraction1*eps1c**2*eps2c**2*
    eps3r-fraction1*eps1c**2*eps2c**2*eps3c**2-4*fraction2*eps1r**2*eps2r*eps3r-4*
    fraction2*eps1r*eps2c**2*eps3r**2-fraction2*eps1r**2*eps2r*eps3c**2-4*fraction2*
    eps1r**2*eps2c**2*eps3r-16*fraction3**2*eps1r**2*eps3r)/(64+4*eps2c**2*fraction2
    **2*eps3c**2+64*fraction1+64*fraction2+64*fraction3+64*eps1r+64*eps2r+64*eps3r+
    16*eps1c**2+16*eps2c**2+16*eps3c**2-8*fraction3**2*eps1c**2*eps3r-8*fraction3**2
    *eps2r**2*eps3r+18*eps1c*fraction1*eps3c*fraction3*eps2r**2+4*eps2c**2*fraction2
    **2*eps3r**2+4*eps3c**2*fraction3**2*eps2r**2+4*fraction2**2*eps1c**2*eps2r**2-
    16*fraction2*eps2r*fraction1+8*fraction3*eps1r**2*fraction1*eps3r**2-16*
    fraction3*eps1r*eps2r**2*fraction2-16*fraction3*eps1r*eps2c**2*fraction2-4*
    fraction3*eps1c**2*eps2r**2*fraction2-8*fraction3**2*eps1c**2*eps2r*eps3r-2*
    fraction3**2*eps1c**2*eps2c**2*eps3r-4*fraction1*eps2r**2*eps3r**2*fraction2-4*
    fraction1*eps2r**2*eps3r**2*fraction3-16*fraction1*eps2r**2*eps3r*fraction2+8*
    fraction3*eps2r**2*eps3r*fraction2+4*eps3c**2*fraction3**2*eps1r*eps2r**2+16*
    eps2c**2*fraction2**2*eps1r*eps3r+4*eps2c**2*fraction2**2*eps1r**2*eps3r+16*
    eps3c**2*fraction3**2*eps1r*eps2r+8*fraction2*fraction1*eps1c**2*eps2r**2+8*
    fraction2*fraction1*eps1c**2*eps2c**2-4*fraction2*eps1c**2*eps2r*fraction3-16*
    fraction2*eps1r*eps3r**2*fraction3-16*fraction2*eps1r*eps3r*fraction3-16*
    fraction2*eps1r*eps3c**2*fraction3+eps2c**2*fraction2**2*eps1r**2*eps3c**2+8*
    fraction2*eps1c**2*eps2r*fraction1-4*fraction2*eps1r*eps3r**2*fraction1-16*
    fraction2*eps1r*eps3r*fraction1+4*eps3c**2*fraction3**2*eps1c**2+8*fraction3*
    eps2c**2*fraction2*eps3r**2-4*fraction2*eps1c**2*eps3r**2*fraction1-16*fraction2
    *eps1c**2*eps3r*fraction1-4*fraction2*eps1c**2*eps3c**2*fraction1+8*fraction3*
    eps2r**2*fraction2*eps3c**2+8*fraction3*eps2r**2*fraction2*eps3r**2-4*fraction1*
    eps2r**2*eps3c**2*fraction2-4*fraction1*eps2r**2*eps3c**2*fraction3-4*fraction1*
    eps2c**2*eps3r**2*fraction2-4*fraction1*eps2c**2*eps3r**2*fraction3-16*fraction1
    *eps2c**2*eps3r*fraction2-4*fraction1*eps2c**2*eps3c**2*fraction2-4*fraction1*
    eps2c**2*eps3c**2*fraction3-4*fraction3*eps1c**2*eps2c**2*fraction2-4*fraction2*
    eps1r**2*eps3r**2*fraction1-16*fraction2*eps1r**2*eps3r*fraction1-4*fraction2*
    eps1r**2*eps3c**2*fraction1+eps1c**2*fraction1**2*eps2r**2*eps3r**2+eps2c**2*
    fraction2**2*eps1r**2*eps3r**2+eps3c**2*fraction3**2*eps1r**2*eps2r**2+8*
    fraction3*fraction1*eps1c**2*eps3r**2+8*fraction3*fraction1*eps1c**2*eps3c**2+16
    *fraction3**2*eps1r*eps2r*eps3r**2+4*fraction3**2*eps1c**2*eps2r*eps3r**2+
    fraction3**2*eps1c**2*eps2c**2*eps3r**2-4*fraction3*eps1r**2*eps2r**2*fraction2-
    4*fraction3*eps1r**2*eps2c**2*fraction2-2*fraction3**2*eps1c**2*eps2r**2*eps3r-2
    *fraction3**2*eps1r**2*eps2r**2*eps3r-8*fraction3**2*eps1r**2*eps2r*eps3r+8*
    fraction3*eps1r**2*fraction1*eps3c**2+4*fraction2**2*eps3c**2*eps2r**2+4*
    fraction1**2*eps1r**2*eps2r*eps3c**2-16*fraction3*eps1r**2*eps2r*fraction1-4*
    fraction3*eps1r**2*eps2c**2*fraction1-4*fraction3*eps1r*eps2r**2*fraction1-16*
    fraction3*eps1r*eps2r*fraction1+72*eps1c*fraction1*eps3c*fraction3+72*eps2c*
    fraction2*eps3c*fraction3-4*fraction3*eps1r**2*eps2r**2*fraction1+72*eps1c*
    fraction1*eps2c*fraction2-2*fraction3**2*eps1r**2*eps2c**2*eps3r+4*fraction3**2*
    eps1r*eps2r**2*eps3r**2+8*fraction3*eps2c**2*fraction2*eps3c**2+16*eps1c**2*
    fraction1**2+16*fraction3**2*eps1r+16*fraction3**2*eps2r-32*fraction3**2*eps3r+4
    *fraction3**2*eps1r**2+4*fraction3**2*eps2r**2+4*fraction3**2*eps1c**2+4*
    fraction3**2*eps2c**2+16*fraction2**2*eps1r-32*fraction2**2*eps2r+16*fraction2**
    2*eps3r+4*fraction2**2*eps1r**2+4*fraction2**2*eps3r**2+4*fraction2**2*eps1c**2+
    4*fraction2**2*eps3c**2+16*eps3c**2*fraction3**2+4*eps1r**2*eps2r**2+16*eps1r**2
    *eps2r+4*eps1r**2*eps2c**2+16*eps1r*eps2r**2+64*eps1r*eps2r+16*eps1r*eps2c**2+4*
    eps1c**2*eps2r**2+16*eps1c**2*eps2r+16*eps1r**2+16*eps2r**2+16*eps3r**2+16*
    fraction1**2+16*fraction2**2+16*fraction3**2-2*fraction2*eps1r**2*eps2r**2*eps3c
    **2-8*fraction2*eps1r*eps2r**2*eps3r**2-2*fraction2*eps1r**2*eps2c**2*eps3c**2-8
    *fraction2*eps1r*eps2r*eps3c**2-2*fraction2*eps1r**2*eps2c**2*eps3r**2-32*
    fraction2*eps1r*eps2r*eps3r-2*fraction2*eps1r**2*eps2r**2*eps3r**2-32*fraction2*
    eps1r*eps2r**2*eps3r-2*fraction2*eps1c**2*eps2r**2*eps3r**2-32*fraction2*eps1r*
    eps2c**2*eps3r-2*fraction2*eps1r**2*eps2r*eps3r**2-8*fraction2*eps1r*eps2r**2*
    eps3c**2-8*fraction2*eps1c**2*eps2r**2*eps3r-8*fraction2*eps1r*eps2r*eps3r**2-8*
    fraction2*eps1r**2*eps2r**2*eps3r-8*fraction2*eps1r*eps2c**2*eps3c**2-2*
    fraction2*eps1c**2*eps2r**2*eps3c**2-2*fraction2*eps1c**2*eps2r*eps3r**2-8*
    fraction2*eps1c**2*eps2r*eps3r-2*fraction2*eps1c**2*eps2r*eps3c**2-2*fraction2*
    eps1c**2*eps2c**2*eps3r**2-8*fraction2*eps1c**2*eps2c**2*eps3r-2*fraction2*eps1c
    **2*eps2c**2*eps3c**2-8*fraction3*eps1r**2*eps2r*eps3r-8*fraction3*eps1r*eps2c**
    2*eps3r**2-8*fraction3*eps1r**2*eps2r*eps3c**2-2*fraction3*eps1r**2*eps2c**2*
    eps3r-2*fraction3*eps1r**2*eps2r**2*eps3c**2-8*fraction3*eps1r*eps2r**2*eps3r**2
    -32*fraction3*eps1r*eps2r*eps3c**2+18*eps1c*fraction1*eps3c*fraction3*eps2c**2+
    72*eps2c*fraction2*eps3c*fraction3*eps1r+18*eps2c*fraction2*eps3c*fraction3*
    eps1r**2+18*eps2c*fraction2*eps3c*fraction3*eps1c**2+2*fraction2*eps1r**2*eps2r*
    eps3r**2*fraction3+72*eps1c*fraction1*eps2c*fraction2*eps3r+18*eps1c*fraction1*
    eps2c*fraction2*eps3r**2+18*eps1c*fraction1*eps2c*fraction2*eps3c**2+72*eps1c*
    fraction1*eps3c*fraction3*eps2r+2*fraction2*eps1r**2*eps2r**2*eps3r**2*fraction1
    +8*fraction2*eps1r*eps2r*eps3c**2*fraction3-2*fraction3*eps1r**2*eps2c**2*eps3r**2
    -32*fraction3*eps1r*eps2r*eps3r-2*fraction3*eps1r**2*eps2r**2*eps3r**2-8*
    fraction3*eps1r*eps2r**2*eps3r-2*fraction3*eps1c**2*eps2r**2*eps3r**2-8*
    fraction3*eps1r*eps2c**2*eps3r-8*fraction3*eps1r**2*eps2r*eps3r**2-8*fraction3*
    eps1r*eps2r**2*eps3c**2-2*fraction3*eps1c**2*eps2r**2*eps3r-32*fraction3*eps1r*
    eps2r*eps3r**2-2*fraction3*eps1r**2*eps2r**2*eps3r-8*fraction3*eps1r*eps2c**2*
    eps3c**2-2*fraction3*eps1c**2*eps2r**2*eps3c**2-8*fraction3*eps1c**2*eps2r*eps3r
    **2-8*fraction3*eps1c**2*eps2r*eps3r-8*fraction3*eps1c**2*eps2r*eps3c**2-2*
    fraction3*eps1c**2*eps2c**2*eps3r**2-2*fraction3*eps1c**2*eps2c**2*eps3r-2*
    fraction3*eps1c**2*eps2c**2*eps3c**2+fraction2**2*eps1r**2*eps2r**2*eps3c**2+
    fraction2**2*eps1r**2*eps2r**2*eps3r**2+4*fraction2**2*eps1r*eps2r**2*eps3r**2+
    fraction2**2*eps1c**2*eps2r**2*eps3r**2+16*fraction2**2*eps1r*eps2r**2*eps3r-2*
    fraction2**2*eps1r**2*eps2r*eps3r**2-8*fraction2**2*eps1r*eps2r*eps3r**2+4*eps3c
    **2*fraction3**2*eps1r**2*eps2r+eps2c**2*fraction2**2*eps1c**2*eps3c**2+4*eps3c**2
    *fraction3**2*eps1r*eps2c**2-8*fraction2**2*eps1r*eps2r*eps3c**2+4*eps1c**2*
    eps2c**2+16*eps1r*eps3r**2+64*eps1r*eps3r+16*eps1r*eps3c**2+16*eps2r*eps3r**2+64
    *eps2r*eps3r+16*eps2r*eps3c**2+4*eps1r**2*eps3r**2+16*eps1r**2*eps3r+4*eps1r**2*
    eps3c**2+4*eps2r**2*eps3r**2+16*eps2r**2*eps3r+4*eps2r**2*eps3c**2+4*eps1c**2*
    eps3r**2+16*eps1c**2*eps3r+4*eps1c**2*eps3c**2+4*eps2c**2*eps3r**2+16*eps2c**2*
    eps3r+4*eps2c**2*eps3c**2-32*fraction1*eps1r+64*fraction1*eps2r+64*fraction1*
    eps3r-32*fraction1*eps1r**2+16*fraction1*eps2r**2+16*fraction1*eps3r**2-32*
    fraction1*eps1c**2+16*fraction1*eps2c**2+16*fraction1*eps3c**2+64*fraction2*
    eps1r-32*fraction2*eps2r+64*fraction2*eps3r+16*fraction2*eps1r**2-32*fraction2*
    eps2r**2+16*fraction2*eps3r**2+16*fraction2*eps1c**2-32*fraction2*eps2c**2+16*
    fraction2*eps3c**2+64*fraction3*eps1r+64*fraction3*eps2r-32*fraction3*eps3r+16*
    fraction3*eps1r**2+16*fraction3*eps2r**2-32*fraction3*eps3r**2+16*fraction3*
    eps1c**2+16*fraction3*eps2c**2-32*fraction3*eps3c**2+16*fraction1**2*eps1r**2+16
    *fraction2**2*eps2r**2+16*eps2c**2*fraction2**2+32*fraction2*fraction1+32*
    fraction3*fraction1+32*fraction2*fraction3+16*fraction3**2*eps3r**2-32*fraction2
    **2*eps1r*eps2r*eps3r+4*fraction2**2*eps1r**2*eps2r**2*eps3r+4*fraction2**2*
    eps2r**2*eps1r*eps3c**2+eps3c**2*fraction3**2*eps1c**2*eps2c**2+eps1c**2*
    fraction1**2*eps2c**2*eps3r**2+fraction3**2*eps1c**2*eps2r**2*eps3r**2+eps3c**2*
    fraction3**2*eps1r**2*eps2c**2+4*eps2c**2*fraction2**2*eps1r*eps3r**2+4*eps1c**2
    *fraction1**2*eps2r*eps3c**2+eps1c**2*fraction1**2*eps2r**2*eps3c**2+4*fraction2
    **2*eps1c**2*eps2r**2*eps3r+fraction2**2*eps1c**2*eps2r**2*eps3c**2+4*fraction3**2
    *eps1r*eps2c**2*eps3r**2+4*eps2c**2*fraction2**2*eps1c**2*eps3r+8*fraction1*
    eps1r*eps2r**2*fraction2+8*fraction1*eps1r*eps2c**2*fraction2+4*fraction3**2*
    eps3r**2*eps1r**2*eps2r+fraction3**2*eps1r**2*eps2c**2*eps3r**2+fraction3**2*
    eps2r**2*eps3r**2*eps1r**2-4*fraction3*eps1r*eps2c**2*fraction1+4*fraction1**2*
    eps1r**2*eps2r**2*eps3r+fraction1**2*eps1r**2*eps2r**2*eps3r**2+fraction1**2*
    eps1r**2*eps2r**2*eps3c**2+16*fraction1**2*eps1r**2*eps2r*eps3r+4*fraction1**2*
    eps1r**2*eps2r*eps3r**2+fraction1**2*eps1r**2*eps2c**2*eps3c**2+4*fraction1**2*
    eps1r**2*eps2c**2*eps3r+fraction1**2*eps1r**2*eps2c**2*eps3r**2+8*fraction2*
    eps1r*eps2r*eps3r*fraction3+2*fraction2*eps1r**2*eps2r**2*eps3c**2*fraction1+8*
    fraction2*eps1r**2*eps2r**2*eps3r*fraction1+2*fraction2*fraction3*eps1r**2*eps2c
    **2*eps3c**2+2*fraction2*eps1r**2*eps2r*eps3c**2*fraction3+2*fraction3*eps1r**2*
    eps2r**2*eps3r*fraction1+8*fraction3*eps1c**2*eps2r*eps3r*fraction1+2*fraction3*
    eps1c**2*eps2c**2*eps3r*fraction1+2*fraction3*eps1c**2*eps2r**2*eps3r*fraction1+
    2*fraction2*eps1r**2*eps2r*eps3r*fraction3+2*fraction3*eps1r*eps2c**2*eps3r*
    fraction1+2*fraction3*eps1r*eps2r**2*eps3r*fraction1+8*fraction3*eps1r*eps2r*
    eps3r*fraction1+2*fraction3*eps1r**2*eps2r**2*eps3r*fraction2+2*fraction3*eps1c**2
    *eps2c**2*eps3r*fraction2+2*fraction3*eps1r**2*eps2c**2*eps3r*fraction2+8*
    fraction1*eps1r*eps2r*eps3r**2*fraction3+2*fraction1*eps1r*eps2r**2*eps3c**2*
    fraction2+2*fraction1*eps1r*eps2r**2*eps3c**2*fraction3+8*fraction1*eps1r*eps2c**2
    *eps3r*fraction2+8*fraction1*eps1r*eps2r**2*eps3r*fraction2+8*fraction1*eps1r*
    eps2r*eps3c**2*fraction3+2*fraction1*eps1r*eps2r**2*eps3r**2*fraction2+2*
    fraction1*eps1r*eps2r**2*eps3r**2*fraction3+2*fraction2*eps3r**2*fraction1*eps1r
    **2*eps2c**2-8*fraction1**2*eps1r*eps2r**2-32*fraction1**2*eps1r*eps2r+16*eps1r**2
    *eps2r*eps3r-8*fraction1**2*eps1r*eps2c**2-8*fraction1**2*eps1r*eps3r**2-32*
    fraction1**2*eps1r*eps3r-8*fraction1**2*eps1r*eps3c**2+8*fraction2*eps3r**2*
    fraction1+16*eps2c**2*fraction2**2*eps1r+8*fraction3*eps2r**2*fraction1-16*
    fraction2*eps1r**2*fraction1+4*eps2c**2*fraction2**2*eps1c**2-16*fraction3*eps1r
    **2*fraction1+4*eps1c**2*fraction1**2*eps3r**2+4*eps1r*eps2c**2*eps3r**2+4*eps1r
    **2*eps2r*eps3c**2+4*eps1r**2*eps2c**2*eps3r+eps1r**2*eps2r**2*eps3c**2+4*eps1r*
    eps2r**2*eps3r**2+eps1r**2*eps2c**2*eps3c**2+16*eps1r*eps2r*eps3c**2+eps1r**2*
    eps2c**2*eps3r**2+64*eps1r*eps2r*eps3r+eps1r**2*eps2r**2*eps3r**2+16*eps1r*eps2r
    **2*eps3r+eps1c**2*eps2r**2*eps3r**2+16*eps1r*eps2c**2*eps3r+4*eps1r**2*eps2r*
    eps3r**2+4*eps1r*eps2r**2*eps3c**2+4*eps1c**2*eps2r**2*eps3r+16*eps1r*eps2r*
    eps3r**2+4*eps1r**2*eps2r**2*eps3r+4*eps1r*eps2c**2*eps3c**2+eps1c**2*eps2r**2*
    eps3c**2+4*eps1c**2*eps2r*eps3r**2+16*eps1c**2*eps2r*eps3r+4*eps1c**2*eps2r*
    eps3c**2+eps1c**2*eps2c**2*eps3r**2+4*eps1c**2*eps2c**2*eps3r+eps1c**2*eps2c**2*
    eps3c**2-8*fraction1*eps1r**2*eps2r**2-32*fraction1*eps1r**2*eps2r-8*fraction1*
    eps1r**2*eps2c**2-8*fraction1*eps1r*eps2r**2-32*fraction1*eps1r*eps2r-8*
    fraction1*eps1r*eps2c**2-8*fraction1*eps1c**2*eps2r**2-32*fraction1*eps1c**2*
    eps2r-8*fraction1*eps1c**2*eps2c**2-8*fraction1*eps1r*eps3r**2+2*fraction2*eps3r
    **2*fraction1*eps1c**2*eps2r**2+2*fraction2*eps3r**2*fraction1*eps1c**2*eps2c**2
    +2*fraction2*eps3c**2*fraction1*eps1r**2*eps2c**2+2*fraction2*eps3c**2*fraction1
    *eps1c**2*eps2r**2+2*fraction2*eps3c**2*fraction1*eps1c**2*eps2c**2+2*fraction1*
    eps1r*eps2c**2*eps3r**2*fraction2+2*fraction1*eps1r*eps2c**2*eps3r**2*fraction3+
    8*fraction2*eps3r*fraction1*eps1r**2*eps2c**2+8*fraction2*eps3r*fraction1*eps1c**2
    *eps2r**2+8*fraction2*eps3r*fraction1*eps1c**2*eps2c**2+2*fraction1*eps1r*
    eps2c**2*eps3c**2*fraction2+2*fraction1*eps1r*eps2c**2*eps3c**2*fraction3+2*
    fraction3*eps1r**2*eps2c**2*eps3r*fraction1+8*fraction3*eps1r**2*eps2r*eps3r*
    fraction1+2*fraction2*eps1c**2*eps2r*eps3c**2*fraction1+2*fraction2*eps1c**2*
    eps2r*eps3r**2*fraction1+8*fraction2*eps1c**2*eps2r*eps3r*fraction1+2*fraction2*
    eps1r*eps2r*eps3r**2*fraction1+2*fraction2*eps1r**2*eps2r*eps3r**2*fraction1+2*
    fraction2*eps1r*eps2r*eps3c**2*fraction1+8*fraction2*eps1r*eps2r*eps3r*fraction1
    +2*fraction2*eps1r**2*eps2r*eps3c**2*fraction1+8*fraction2*eps1r**2*eps2r*eps3r*
    fraction1+8*fraction3*fraction2*eps1r*eps2c**2*eps3r**2+2*fraction3*fraction2*
    eps1r**2*eps2r**2*eps3c**2-4*fraction2*eps2r*eps3c**2*fraction1-32*fraction1**2*
    eps1r+16*fraction1**2*eps2r+16*fraction1**2*eps3r+4*fraction1**2*eps2r**2+4*
    fraction1**2*eps3r**2+4*fraction1**2*eps2c**2+4*fraction1**2*eps3c**2+4*
    fraction2**2*eps3r**2*eps2r**2-32*fraction1*eps1r*eps3r-8*fraction1*eps1r*eps3c**2
    +16*fraction1*eps2r*eps3r**2+64*fraction1*eps2r*eps3r+16*fraction1*eps2r*eps3c
    **2-8*fraction1*eps1r**2*eps3r**2-32*fraction1*eps1r**2*eps3r-8*fraction1*eps1r**2
    *eps3c**2+4*fraction1*eps2r**2*eps3r**2+16*fraction1*eps2r**2*eps3r+4*
    fraction1*eps2r**2*eps3c**2-8*fraction1*eps1c**2*eps3r**2-32*fraction1*eps1c**2*
    eps3r-8*fraction1*eps1c**2*eps3c**2+4*fraction1*eps2c**2*eps3r**2+16*fraction1*
    eps2c**2*eps3r+4*fraction1*eps2c**2*eps3c**2-8*fraction2*eps1r**2*eps2r**2-8*
    fraction2*eps1r**2*eps2r-8*fraction2*eps1r**2*eps2c**2-32*fraction2*eps1r*eps2r**2
    -32*fraction2*eps1r*eps2r-32*fraction2*eps1r*eps2c**2-8*fraction2*eps1c**2*
    eps2r**2-8*fraction2*eps1c**2*eps2r-8*fraction2*eps1c**2*eps2c**2+16*fraction2*
    eps1r*eps3r**2+64*fraction2*eps1r*eps3r+16*fraction2*eps1r*eps3c**2-8*fraction2*
    eps2r*eps3r**2-32*fraction2*eps2r*eps3r-8*fraction2*eps2r*eps3c**2+4*fraction2*
    eps1r**2*eps3r**2+16*fraction2*eps1r**2*eps3r+4*fraction2*eps1r**2*eps3c**2-8*
    fraction2*eps2r**2*eps3r**2-32*fraction2*eps2r**2*eps3r-8*fraction2*eps2r**2*
    eps3c**2+4*fraction2*eps1c**2*eps3r**2+16*fraction2*eps1c**2*eps3r+4*fraction2*
    eps1c**2*eps3c**2-8*fraction2*eps2c**2*eps3r**2-32*fraction2*eps2c**2*eps3r-8*
    fraction2*eps2c**2*eps3c**2+4*fraction3*eps1r**2*eps2r**2+16*fraction3*eps1r**2*
    eps2r+4*fraction3*eps1r**2*eps2c**2+16*fraction3*eps1r*eps2r**2+64*fraction3*
    eps1r*eps2r+16*fraction3*eps1r*eps2c**2+4*fraction3*eps1c**2*eps2r**2+16*
    fraction3*eps1c**2*eps2r+4*fraction3*eps1c**2*eps2c**2-32*fraction3*eps1r*eps3r**2
    -32*fraction3*eps1r*eps3r-32*fraction3*eps1r*eps3c**2-32*fraction3*eps2r*eps3r
    **2-32*fraction3*eps2r*eps3r-32*fraction3*eps2r*eps3c**2-8*fraction3*eps1r**2*
    eps3r**2-8*fraction3*eps1r**2*eps3r-8*fraction3*eps1r**2*eps3c**2-8*fraction3*
    eps2r**2*eps3r**2-8*fraction3*eps2r**2*eps3r-8*fraction3*eps2r**2*eps3c**2-8*
    fraction3*eps1c**2*eps3r**2-8*fraction3*eps1c**2*eps3r-8*fraction3*eps1c**2*
    eps3c**2-8*fraction3*eps2c**2*eps3r**2-8*fraction3*eps2c**2*eps3r-8*fraction3*
    eps2c**2*eps3c**2+4*fraction3**2*eps2c**2*eps3r**2+4*fraction3**2*eps1r**2*eps3r
    **2+4*fraction3**2*eps2r**2*eps3r**2+4*fraction3**2*eps1c**2*eps3r**2+16*
    fraction3**2*eps1r*eps3r**2+16*fraction3**2*eps2r*eps3r**2-16*fraction3*eps3r*
    fraction1+8*fraction3*fraction2*eps1r*eps2r**2*eps3r**2+2*fraction3*fraction2*
    eps1r**2*eps2c**2*eps3r**2+2*fraction3*fraction2*eps1r**2*eps2r**2*eps3r**2+2*
    fraction3*fraction2*eps1c**2*eps2r**2*eps3r**2+8*fraction3*eps2r*fraction1*eps1c
    **2*eps3c**2+2*fraction3*eps2r**2*fraction1*eps1c**2*eps3c**2+2*fraction3*eps1c**2
    *fraction1*eps2c**2*eps3r**2+2*fraction3*eps1c**2*fraction1*eps2c**2*eps3c**2+
    8*fraction2*eps1r*eps2r*eps3r**2*fraction3+2*fraction2*eps1c**2*eps2r*eps3r**2*
    fraction3+2*fraction2*eps1c**2*eps2r*eps3r*fraction3+2*fraction2*eps1c**2*eps2r*
    eps3c**2*fraction3+8*fraction3*eps1r**2*fraction1*eps2r*eps3r**2+8*fraction3*
    eps1c**2*fraction1*eps2r*eps3r**2+8*fraction3*fraction2*eps1r*eps2c**2*eps3c**2+
    2*fraction3*fraction2*eps1c**2*eps2r**2*eps3c**2+2*fraction3*fraction2*eps1c**2*
    eps2c**2*eps3r**2+2*fraction3*fraction2*eps1c**2*eps2c**2*eps3c**2+8*fraction3*
    eps2r**2*fraction2*eps1r*eps3c**2+2*fraction3*eps1r**2*eps2r**2*fraction1*eps3r**2
    +2*fraction3*eps1r**2*eps2r**2*fraction1*eps3c**2+8*fraction3*eps1r**2*eps2r*
    fraction1*eps3c**2+2*fraction3*eps1r**2*eps2c**2*fraction1*eps3r**2+2*fraction3*
    eps1r**2*eps2c**2*fraction1*eps3c**2+2*fraction3*eps1c**2*eps2r**2*fraction1*
    eps3r**2+8*fraction3*eps1r*eps2r**2*eps3r*fraction2+8*fraction3*eps1r*eps2c**2*
    eps3r*fraction2+2*fraction3*eps1c**2*eps2r**2*eps3r*fraction2-16*fraction3*eps1r
    *fraction1-16*fraction2*eps1c**2*fraction1+8*fraction2*eps3c**2*fraction1-16*
    fraction3*eps1c**2*fraction1+8*fraction3*eps2c**2*fraction1-16*fraction2*eps1r*
    fraction1+4*eps2c**2*fraction2**2*eps1r**2+4*eps1c**2*fraction1**2*eps2c**2+16*
    eps3c**2*fraction3**2*eps1r+4*eps1c**2*fraction1**2*eps3c**2+fraction1**2*eps2c**2
    *eps3c**2+fraction1**2*eps2c**2*eps3r**2+4*fraction1**2*eps2c**2*eps3r+
    fraction1**2*eps2r**2*eps3c**2+4*fraction1**2*eps2r**2*eps3r+fraction1**2*eps2r**2
    *eps3r**2+16*fraction1**2*eps2r*eps3r+4*fraction1**2*eps2r*eps3c**2+4*
    fraction1**2*eps2r*eps3r**2+16*eps2c**2*fraction2**2*eps3r+16*eps1c**2*fraction1
    **2*eps3r-16*fraction2*eps2r*fraction3-16*fraction2*eps3r*fraction3+8*fraction2*
    eps1r**2*fraction3-16*fraction2*eps3r**2*fraction3+16*fraction1**2*eps1r**2*
    eps2r+4*fraction1**2*eps1r**2*eps2c**2+4*eps3c**2*fraction3**2*eps2c**2+4*
    fraction1**2*eps1r**2*eps2r**2-16*fraction1*eps2r**2*fraction2-16*fraction1*
    eps3r**2*fraction3-16*fraction1*eps2c**2*fraction2-16*fraction1*eps3c**2*
    fraction3+16*eps3c**2*fraction3**2*eps2r+32*fraction2*eps3r*fraction1+4*
    fraction2**2*eps1r**2*eps2r**2+4*fraction3**2*eps1r**2*eps2r+fraction3**2*eps1r**2
    *eps2c**2+4*fraction3**2*eps1r*eps2r**2-4*fraction2*eps1r*eps3c**2*fraction1-4
    *fraction2*eps2r*eps3r**2*fraction1-16*fraction2*eps2r*eps3r*fraction1-32*
    fraction3**2*eps1r*eps2r*eps3r-8*fraction3**2*eps1r*eps2r**2*eps3r-8*fraction3**
    2*eps1r*eps2c**2*eps3r-4*fraction2*eps1c**2*eps3r*fraction3-4*fraction2*eps1c**2
    *eps3c**2*fraction3+4*eps1c**2*fraction1**2*eps2r**2*eps3r-4*fraction2*eps1r**2*
    eps3r**2*fraction3-4*fraction2*eps1r**2*eps3r*fraction3-4*fraction2*eps1r**2*
    eps3c**2*fraction3-4*fraction2*eps1c**2*eps3r**2*fraction3+8*fraction2*eps1r**2*
    eps2r*fraction1+8*fraction2*eps1r*eps2r*fraction1+8*fraction2*eps2r*eps3r**2*
    fraction3+8*fraction2*eps2r*eps3r*fraction3+8*fraction2*eps2r*eps3c**2*fraction3
    -8*fraction1**2*eps1r*eps2r*eps3r**2-2*fraction1**2*eps1r*eps2c**2*eps3c**2+16*
    eps1c**2*fraction1**2*eps2r*eps3r-2*fraction1**2*eps1r*eps2r**2*eps3c**2-32*
    fraction1**2*eps1r*eps2r*eps3r-8*fraction1**2*eps1r*eps2r**2*eps3r-8*fraction1**
    2*eps1r*eps2c**2*eps3r-2*fraction1**2*eps1r*eps2c**2*eps3r**2-2*fraction1**2*
    eps1r*eps2r**2*eps3r**2-8*fraction1**2*eps1r*eps2r*eps3c**2-4*fraction3*eps2r**2
    *eps3r*fraction1+8*fraction3*eps1c**2*eps3r*fraction1-4*fraction3*eps2c**2*eps3r
    *fraction1-16*fraction3*eps1c**2*eps2r*fraction1-4*fraction3*eps1c**2*eps2c**2*
    fraction1+8*fraction3*eps1r*eps3r*fraction1-16*fraction3*eps2r*eps3r*fraction1+8
    *fraction3*eps1r**2*eps3r*fraction1-4*fraction3*eps1c**2*eps2r**2*fraction1+
    eps2c**2*fraction2**2*eps1c**2*eps3r**2-2*fraction2**2*eps1r**2*eps2r*eps3c**2-8
    *fraction2**2*eps1r**2*eps2r*eps3r+4*eps1c**2*fraction1**2*eps2r*eps3r**2-16*
    fraction1*eps2r*eps3r**2*fraction3-16*fraction1*eps2r*eps3c**2*fraction3+8*
    fraction3*eps2c**2*eps3r*fraction2+8*fraction1*eps1r*eps3r**2*fraction3+8*
    fraction1*eps1r*eps3c**2*fraction3+4*eps3c**2*fraction3**2*eps1c**2*eps2r+8*
    fraction2*fraction1*eps1r**2*eps2c**2+8*fraction2*fraction1*eps1r**2*eps2r**2+
    eps3c**2*fraction3**2*eps1c**2*eps2r**2+4*eps2c**2*fraction2**2*eps1r*eps3c**2-2
    *fraction2**2*eps1c**2*eps2r*eps3r**2-8*fraction2**2*eps1c**2*eps2r*eps3r-2*
    fraction2**2*eps1c**2*eps2r*eps3c**2+4*eps1c**2*fraction1**2*eps2c**2*eps3r+
    eps1c**2*fraction1**2*eps2c**2*eps3c**2-4*fraction2*eps1r**2*eps2r*fraction3-16*
    fraction2*eps1r*eps2r*fraction3-2*fraction3*eps1r**2*eps2c**2*eps3c**2-32*
    fraction1*eps1r**2*eps2r*eps3r+16*fraction3**2*eps1r*eps2r+4*fraction3**2*eps1r*
    eps2c**2+fraction3**2*eps1c**2*eps2r**2+4*fraction3**2*eps1c**2*eps2r+fraction3**2
    *eps1c**2*eps2c**2-32*fraction3**2*eps1r*eps3r-32*fraction3**2*eps2r*eps3r+
    fraction2**2*eps1r**2*eps3r**2+4*fraction2**2*eps1r**2*eps3r+fraction2**2*eps1r**2
    *eps3c**2+fraction2**2*eps1c**2*eps3r**2+4*fraction2**2*eps1c**2*eps3r+
    fraction2**2*eps1c**2*eps3c**2+fraction3**2*eps1r**2*eps2r**2-8*fraction2**2*
    eps1r**2*eps2r-32*fraction2**2*eps1r*eps2r-8*fraction2**2*eps1c**2*eps2r+4*
    fraction2**2*eps1r*eps3r**2+16*fraction2**2*eps1r*eps3r+4*fraction2**2*eps1r*
    eps3c**2-8*fraction2**2*eps2r*eps3r**2-32*fraction2**2*eps2r*eps3r-8*fraction2**
    2*eps2r*eps3c**2+16*eps1c**2*fraction1**2*eps2r-16*fraction3*eps2r**2*fraction2-
    16*fraction3*eps2c**2*fraction2+4*fraction1**2*eps1r**2*eps3r**2+4*eps3c**2*
    fraction3**2*eps1r**2+4*fraction1**2*eps1r**2*eps3c**2+16*fraction1**2*eps1r**2*
    eps3r+16*fraction2**2*eps1r*eps2r**2+4*eps1c**2*fraction1**2*eps2r**2+16*
    fraction2**2*eps3r*eps2r**2+32*fraction2*eps1r*fraction3-8*fraction3**2*eps2c**2
    *eps3r+8*fraction2*eps1c**2*fraction3-16*fraction2*eps3c**2*fraction3+32*
    fraction3*eps2r*fraction1-2*fraction1*eps1r*eps2c**2*eps3r**2-8*fraction1*eps1r**2
    *eps2r*eps3c**2-8*fraction1*eps1r**2*eps2c**2*eps3r-2*fraction1*eps1r**2*eps2r
    **2*eps3c**2-2*fraction1*eps1r*eps2r**2*eps3r**2-2*fraction1*eps1r**2*eps2c**2*
    eps3c**2-8*fraction1*eps1r*eps2r*eps3c**2-2*fraction1*eps1r**2*eps2c**2*eps3r**2
    -32*fraction1*eps1r*eps2r*eps3r-2*fraction1*eps1r**2*eps2r**2*eps3r**2-8*
    fraction1*eps1r*eps2r**2*eps3r-2*fraction1*eps1c**2*eps2r**2*eps3r**2-8*
    fraction1*eps1r*eps2c**2*eps3r-8*fraction1*eps1r**2*eps2r*eps3r**2-2*fraction1*
    eps1r*eps2r**2*eps3c**2-8*fraction1*eps1c**2*eps2r**2*eps3r-8*fraction1*eps1r*
    eps2r*eps3r**2-8*fraction1*eps1r**2*eps2r**2*eps3r-2*fraction1*eps1r*eps2c**2*
    eps3c**2-2*fraction1*eps1c**2*eps2r**2*eps3c**2-8*fraction1*eps1c**2*eps2r*eps3r
    **2-32*fraction1*eps1c**2*eps2r*eps3r-8*fraction1*eps1c**2*eps2r*eps3c**2-2*
    fraction1*eps1c**2*eps2c**2*eps3r**2-8*fraction1*eps1c**2*eps2c**2*eps3r-2*
    fraction1*eps1c**2*eps2c**2*eps3c**2-8*fraction2*eps1r**2*eps2r*eps3r-8*
    fraction2*eps1r*eps2c**2*eps3r**2-2*fraction2*eps1r**2*eps2r*eps3c**2-8*
    fraction2*eps1r**2*eps2c**2*eps3r-8*fraction3**2*eps1r**2*eps3r)
  epsilon_effective_imag = \
    9*(16*eps1c*fraction1+16*eps2c*fraction2+16*eps3c*fraction3+16*eps1c*fraction1*
    eps2r+16*eps1c*fraction1*eps3r+4*eps1c*fraction1*eps2r**2+4*eps1c*fraction1*
    eps3r**2+4*eps1c*fraction1*eps2c**2+4*eps1c*fraction1*eps3c**2+16*eps2c*
    fraction2*eps1r+16*eps2c*fraction2*eps3r+4*eps2c*fraction2*eps1r**2+4*eps2c*
    fraction2*eps3r**2+4*eps2c*fraction2*eps1c**2+4*eps2c*fraction2*eps3c**2+16*
    eps3c*fraction3*eps1r+16*eps3c*fraction3*eps2r+4*eps3c*fraction3*eps1r**2+4*
    eps3c*fraction3*eps2r**2+4*eps3c*fraction3*eps1c**2+4*eps3c*fraction3*eps2c**2+4
    *eps1c*fraction1*eps2r*eps3r**2+16*eps1c*fraction1*eps2r*eps3r+4*eps1c*fraction1
    *eps2r*eps3c**2+eps1c*fraction1*eps2r**2*eps3r**2+4*eps1c*fraction1*eps2r**2*
    eps3r+eps1c*fraction1*eps2r**2*eps3c**2+eps1c*fraction1*eps2c**2*eps3r**2+4*
    eps1c*fraction1*eps2c**2*eps3r+eps1c*fraction1*eps2c**2*eps3c**2+4*eps2c*
    fraction2*eps1r*eps3r**2+16*eps2c*fraction2*eps1r*eps3r+4*eps2c*fraction2*eps1r*
    eps3c**2+eps2c*fraction2*eps1r**2*eps3r**2+4*eps2c*fraction2*eps1r**2*eps3r+
    eps2c*fraction2*eps1r**2*eps3c**2+eps2c*fraction2*eps1c**2*eps3r**2+4*eps2c*
    fraction2*eps1c**2*eps3r+eps2c*fraction2*eps1c**2*eps3c**2+eps3c*fraction3*eps1r
    **2*eps2r**2+4*eps3c*fraction3*eps1r**2*eps2r+eps3c*fraction3*eps1r**2*eps2c**2+
    4*eps3c*fraction3*eps1r*eps2r**2+16*eps3c*fraction3*eps1r*eps2r+4*eps3c*
    fraction3*eps1r*eps2c**2+eps3c*fraction3*eps1c**2*eps2r**2+4*eps3c*fraction3*
    eps1c**2*eps2r+eps3c*fraction3*eps1c**2*eps2c**2)/(64+4*eps2c**2*fraction2**2*
    eps3c**2+64*fraction1+64*fraction2+64*fraction3+64*eps1r+64*eps2r+64*eps3r+16*
    eps1c**2+16*eps2c**2+16*eps3c**2-8*fraction3**2*eps1c**2*eps3r-8*fraction3**2*
    eps2r**2*eps3r+18*eps1c*fraction1*eps3c*fraction3*eps2r**2+4*eps2c**2*fraction2**2
    *eps3r**2+4*eps3c**2*fraction3**2*eps2r**2+4*fraction2**2*eps1c**2*eps2r**2-16
    *fraction2*eps2r*fraction1+8*fraction3*eps1r**2*fraction1*eps3r**2-16*fraction3*
    eps1r*eps2r**2*fraction2-16*fraction3*eps1r*eps2c**2*fraction2-4*fraction3*eps1c
    **2*eps2r**2*fraction2-8*fraction3**2*eps1c**2*eps2r*eps3r-2*fraction3**2*eps1c**2
    *eps2c**2*eps3r-4*fraction1*eps2r**2*eps3r**2*fraction2-4*fraction1*eps2r**2*
    eps3r**2*fraction3-16*fraction1*eps2r**2*eps3r*fraction2+8*fraction3*eps2r**2*
    eps3r*fraction2+4*eps3c**2*fraction3**2*eps1r*eps2r**2+16*eps2c**2*fraction2**2*
    eps1r*eps3r+4*eps2c**2*fraction2**2*eps1r**2*eps3r+16*eps3c**2*fraction3**2*
    eps1r*eps2r+8*fraction2*fraction1*eps1c**2*eps2r**2+8*fraction2*fraction1*eps1c**2
    *eps2c**2-4*fraction2*eps1c**2*eps2r*fraction3-16*fraction2*eps1r*eps3r**2*
    fraction3-16*fraction2*eps1r*eps3r*fraction3-16*fraction2*eps1r*eps3c**2*
    fraction3+eps2c**2*fraction2**2*eps1r**2*eps3c**2+8*fraction2*eps1c**2*eps2r*
    fraction1-4*fraction2*eps1r*eps3r**2*fraction1-16*fraction2*eps1r*eps3r*
    fraction1+4*eps3c**2*fraction3**2*eps1c**2+8*fraction3*eps2c**2*fraction2*eps3r**2
    -4*fraction2*eps1c**2*eps3r**2*fraction1-16*fraction2*eps1c**2*eps3r*fraction1
    -4*fraction2*eps1c**2*eps3c**2*fraction1+8*fraction3*eps2r**2*fraction2*eps3c**2
    +8*fraction3*eps2r**2*fraction2*eps3r**2-4*fraction1*eps2r**2*eps3c**2*fraction2
    -4*fraction1*eps2r**2*eps3c**2*fraction3-4*fraction1*eps2c**2*eps3r**2*fraction2
    -4*fraction1*eps2c**2*eps3r**2*fraction3-16*fraction1*eps2c**2*eps3r*fraction2-4
    *fraction1*eps2c**2*eps3c**2*fraction2-4*fraction1*eps2c**2*eps3c**2*fraction3-4
    *fraction3*eps1c**2*eps2c**2*fraction2-4*fraction2*eps1r**2*eps3r**2*fraction1-
    16*fraction2*eps1r**2*eps3r*fraction1-4*fraction2*eps1r**2*eps3c**2*fraction1+
    eps1c**2*fraction1**2*eps2r**2*eps3r**2+eps2c**2*fraction2**2*eps1r**2*eps3r**2+
    eps3c**2*fraction3**2*eps1r**2*eps2r**2+8*fraction3*fraction1*eps1c**2*eps3r**2+
    8*fraction3*fraction1*eps1c**2*eps3c**2+16*fraction3**2*eps1r*eps2r*eps3r**2+4*
    fraction3**2*eps1c**2*eps2r*eps3r**2+fraction3**2*eps1c**2*eps2c**2*eps3r**2-4*
    fraction3*eps1r**2*eps2r**2*fraction2-4*fraction3*eps1r**2*eps2c**2*fraction2-2*
    fraction3**2*eps1c**2*eps2r**2*eps3r-2*fraction3**2*eps1r**2*eps2r**2*eps3r-8*
    fraction3**2*eps1r**2*eps2r*eps3r+8*fraction3*eps1r**2*fraction1*eps3c**2+4*
    fraction2**2*eps3c**2*eps2r**2+4*fraction1**2*eps1r**2*eps2r*eps3c**2-16*
    fraction3*eps1r**2*eps2r*fraction1-4*fraction3*eps1r**2*eps2c**2*fraction1-4*
    fraction3*eps1r*eps2r**2*fraction1-16*fraction3*eps1r*eps2r*fraction1+72*eps1c*
    fraction1*eps3c*fraction3+72*eps2c*fraction2*eps3c*fraction3-4*fraction3*eps1r**
    2*eps2r**2*fraction1+72*eps1c*fraction1*eps2c*fraction2-2*fraction3**2*eps1r**2*
    eps2c**2*eps3r+4*fraction3**2*eps1r*eps2r**2*eps3r**2+8*fraction3*eps2c**2*
    fraction2*eps3c**2+16*eps1c**2*fraction1**2+16*fraction3**2*eps1r+16*fraction3**
    2*eps2r-32*fraction3**2*eps3r+4*fraction3**2*eps1r**2+4*fraction3**2*eps2r**2+4*
    fraction3**2*eps1c**2+4*fraction3**2*eps2c**2+16*fraction2**2*eps1r-32*fraction2
    **2*eps2r+16*fraction2**2*eps3r+4*fraction2**2*eps1r**2+4*fraction2**2*eps3r**2+
    4*fraction2**2*eps1c**2+4*fraction2**2*eps3c**2+16*eps3c**2*fraction3**2+4*eps1r
    **2*eps2r**2+16*eps1r**2*eps2r+4*eps1r**2*eps2c**2+16*eps1r*eps2r**2+64*eps1r*
    eps2r+16*eps1r*eps2c**2+4*eps1c**2*eps2r**2+16*eps1c**2*eps2r+16*eps1r**2+16*
    eps2r**2+16*eps3r**2+16*fraction1**2+16*fraction2**2+16*fraction3**2-2*fraction2
    *eps1r**2*eps2r**2*eps3c**2-8*fraction2*eps1r*eps2r**2*eps3r**2-2*fraction2*
    eps1r**2*eps2c**2*eps3c**2-8*fraction2*eps1r*eps2r*eps3c**2-2*fraction2*eps1r**2
    *eps2c**2*eps3r**2-32*fraction2*eps1r*eps2r*eps3r-2*fraction2*eps1r**2*eps2r**2*
    eps3r**2-32*fraction2*eps1r*eps2r**2*eps3r-2*fraction2*eps1c**2*eps2r**2*eps3r**
    2-32*fraction2*eps1r*eps2c**2*eps3r-2*fraction2*eps1r**2*eps2r*eps3r**2-8*
    fraction2*eps1r*eps2r**2*eps3c**2-8*fraction2*eps1c**2*eps2r**2*eps3r-8*
    fraction2*eps1r*eps2r*eps3r**2-8*fraction2*eps1r**2*eps2r**2*eps3r-8*fraction2*
    eps1r*eps2c**2*eps3c**2-2*fraction2*eps1c**2*eps2r**2*eps3c**2-2*fraction2*eps1c
    **2*eps2r*eps3r**2-8*fraction2*eps1c**2*eps2r*eps3r-2*fraction2*eps1c**2*eps2r*
    eps3c**2-2*fraction2*eps1c**2*eps2c**2*eps3r**2-8*fraction2*eps1c**2*eps2c**2*
    eps3r-2*fraction2*eps1c**2*eps2c**2*eps3c**2-8*fraction3*eps1r**2*eps2r*eps3r-8*
    fraction3*eps1r*eps2c**2*eps3r**2-8*fraction3*eps1r**2*eps2r*eps3c**2-2*
    fraction3*eps1r**2*eps2c**2*eps3r-2*fraction3*eps1r**2*eps2r**2*eps3c**2-8*
    fraction3*eps1r*eps2r**2*eps3r**2-32*fraction3*eps1r*eps2r*eps3c**2+18*eps1c*
    fraction1*eps3c*fraction3*eps2c**2+72*eps2c*fraction2*eps3c*fraction3*eps1r+18*
    eps2c*fraction2*eps3c*fraction3*eps1r**2+18*eps2c*fraction2*eps3c*fraction3*
    eps1c**2+2*fraction2*eps1r**2*eps2r*eps3r**2*fraction3+72*eps1c*fraction1*eps2c*
    fraction2*eps3r+18*eps1c*fraction1*eps2c*fraction2*eps3r**2+18*eps1c*fraction1*
    eps2c*fraction2*eps3c**2+72*eps1c*fraction1*eps3c*fraction3*eps2r+2*fraction2*
    eps1r**2*eps2r**2*eps3r**2*fraction1+8*fraction2*eps1r*eps2r*eps3c**2*fraction3-
    2*fraction3*eps1r**2*eps2c**2*eps3r**2-32*fraction3*eps1r*eps2r*eps3r-2*
    fraction3*eps1r**2*eps2r**2*eps3r**2-8*fraction3*eps1r*eps2r**2*eps3r-2*
    fraction3*eps1c**2*eps2r**2*eps3r**2-8*fraction3*eps1r*eps2c**2*eps3r-8*
    fraction3*eps1r**2*eps2r*eps3r**2-8*fraction3*eps1r*eps2r**2*eps3c**2-2*
    fraction3*eps1c**2*eps2r**2*eps3r-32*fraction3*eps1r*eps2r*eps3r**2-2*fraction3*
    eps1r**2*eps2r**2*eps3r-8*fraction3*eps1r*eps2c**2*eps3c**2-2*fraction3*eps1c**2
    *eps2r**2*eps3c**2-8*fraction3*eps1c**2*eps2r*eps3r**2-8*fraction3*eps1c**2*
    eps2r*eps3r-8*fraction3*eps1c**2*eps2r*eps3c**2-2*fraction3*eps1c**2*eps2c**2*
    eps3r**2-2*fraction3*eps1c**2*eps2c**2*eps3r-2*fraction3*eps1c**2*eps2c**2*eps3c
    **2+fraction2**2*eps1r**2*eps2r**2*eps3c**2+fraction2**2*eps1r**2*eps2r**2*eps3r
    **2+4*fraction2**2*eps1r*eps2r**2*eps3r**2+fraction2**2*eps1c**2*eps2r**2*eps3r**2
    +16*fraction2**2*eps1r*eps2r**2*eps3r-2*fraction2**2*eps1r**2*eps2r*eps3r**2-8
    *fraction2**2*eps1r*eps2r*eps3r**2+4*eps3c**2*fraction3**2*eps1r**2*eps2r+eps2c**2
    *fraction2**2*eps1c**2*eps3c**2+4*eps3c**2*fraction3**2*eps1r*eps2c**2-8*
    fraction2**2*eps1r*eps2r*eps3c**2+4*eps1c**2*eps2c**2+16*eps1r*eps3r**2+64*eps1r
    *eps3r+16*eps1r*eps3c**2+16*eps2r*eps3r**2+64*eps2r*eps3r+16*eps2r*eps3c**2+4*
    eps1r**2*eps3r**2+16*eps1r**2*eps3r+4*eps1r**2*eps3c**2+4*eps2r**2*eps3r**2+16*
    eps2r**2*eps3r+4*eps2r**2*eps3c**2+4*eps1c**2*eps3r**2+16*eps1c**2*eps3r+4*eps1c
    **2*eps3c**2+4*eps2c**2*eps3r**2+16*eps2c**2*eps3r+4*eps2c**2*eps3c**2-32*
    fraction1*eps1r+64*fraction1*eps2r+64*fraction1*eps3r-32*fraction1*eps1r**2+16*
    fraction1*eps2r**2+16*fraction1*eps3r**2-32*fraction1*eps1c**2+16*fraction1*
    eps2c**2+16*fraction1*eps3c**2+64*fraction2*eps1r-32*fraction2*eps2r+64*
    fraction2*eps3r+16*fraction2*eps1r**2-32*fraction2*eps2r**2+16*fraction2*eps3r**
    2+16*fraction2*eps1c**2-32*fraction2*eps2c**2+16*fraction2*eps3c**2+64*fraction3
    *eps1r+64*fraction3*eps2r-32*fraction3*eps3r+16*fraction3*eps1r**2+16*fraction3*
    eps2r**2-32*fraction3*eps3r**2+16*fraction3*eps1c**2+16*fraction3*eps2c**2-32*
    fraction3*eps3c**2+16*fraction1**2*eps1r**2+16*fraction2**2*eps2r**2+16*eps2c**2
    *fraction2**2+32*fraction2*fraction1+32*fraction3*fraction1+32*fraction2*
    fraction3+16*fraction3**2*eps3r**2-32*fraction2**2*eps1r*eps2r*eps3r+4*fraction2
    **2*eps1r**2*eps2r**2*eps3r+4*fraction2**2*eps2r**2*eps1r*eps3c**2+eps3c**2*
    fraction3**2*eps1c**2*eps2c**2+eps1c**2*fraction1**2*eps2c**2*eps3r**2+fraction3
    **2*eps1c**2*eps2r**2*eps3r**2+eps3c**2*fraction3**2*eps1r**2*eps2c**2+4*eps2c**
    2*fraction2**2*eps1r*eps3r**2+4*eps1c**2*fraction1**2*eps2r*eps3c**2+eps1c**2*
    fraction1**2*eps2r**2*eps3c**2+4*fraction2**2*eps1c**2*eps2r**2*eps3r+fraction2**2
    *eps1c**2*eps2r**2*eps3c**2+4*fraction3**2*eps1r*eps2c**2*eps3r**2+4*eps2c**2*
    fraction2**2*eps1c**2*eps3r+8*fraction1*eps1r*eps2r**2*fraction2+8*fraction1*
    eps1r*eps2c**2*fraction2+4*fraction3**2*eps3r**2*eps1r**2*eps2r+fraction3**2*
    eps1r**2*eps2c**2*eps3r**2+fraction3**2*eps2r**2*eps3r**2*eps1r**2-4*fraction3*
    eps1r*eps2c**2*fraction1+4*fraction1**2*eps1r**2*eps2r**2*eps3r+fraction1**2*
    eps1r**2*eps2r**2*eps3r**2+fraction1**2*eps1r**2*eps2r**2*eps3c**2+16*fraction1**2
    *eps1r**2*eps2r*eps3r+4*fraction1**2*eps1r**2*eps2r*eps3r**2+fraction1**2*
    eps1r**2*eps2c**2*eps3c**2+4*fraction1**2*eps1r**2*eps2c**2*eps3r+fraction1**2*
    eps1r**2*eps2c**2*eps3r**2+8*fraction2*eps1r*eps2r*eps3r*fraction3+2*fraction2*
    eps1r**2*eps2r**2*eps3c**2*fraction1+8*fraction2*eps1r**2*eps2r**2*eps3r*
    fraction1+2*fraction2*fraction3*eps1r**2*eps2c**2*eps3c**2+2*fraction2*eps1r**2*
    eps2r*eps3c**2*fraction3+2*fraction3*eps1r**2*eps2r**2*eps3r*fraction1+8*
    fraction3*eps1c**2*eps2r*eps3r*fraction1+2*fraction3*eps1c**2*eps2c**2*eps3r*
    fraction1+2*fraction3*eps1c**2*eps2r**2*eps3r*fraction1+2*fraction2*eps1r**2*
    eps2r*eps3r*fraction3+2*fraction3*eps1r*eps2c**2*eps3r*fraction1+2*fraction3*
    eps1r*eps2r**2*eps3r*fraction1+8*fraction3*eps1r*eps2r*eps3r*fraction1+2*
    fraction3*eps1r**2*eps2r**2*eps3r*fraction2+2*fraction3*eps1c**2*eps2c**2*eps3r*
    fraction2+2*fraction3*eps1r**2*eps2c**2*eps3r*fraction2+8*fraction1*eps1r*eps2r*
    eps3r**2*fraction3+2*fraction1*eps1r*eps2r**2*eps3c**2*fraction2+2*fraction1*
    eps1r*eps2r**2*eps3c**2*fraction3+8*fraction1*eps1r*eps2c**2*eps3r*fraction2+8*
    fraction1*eps1r*eps2r**2*eps3r*fraction2+8*fraction1*eps1r*eps2r*eps3c**2*
    fraction3+2*fraction1*eps1r*eps2r**2*eps3r**2*fraction2+2*fraction1*eps1r*eps2r**2
    *eps3r**2*fraction3+2*fraction2*eps3r**2*fraction1*eps1r**2*eps2c**2-8*
    fraction1**2*eps1r*eps2r**2-32*fraction1**2*eps1r*eps2r+16*eps1r**2*eps2r*eps3r-
    8*fraction1**2*eps1r*eps2c**2-8*fraction1**2*eps1r*eps3r**2-32*fraction1**2*
    eps1r*eps3r-8*fraction1**2*eps1r*eps3c**2+8*fraction2*eps3r**2*fraction1+16*
    eps2c**2*fraction2**2*eps1r+8*fraction3*eps2r**2*fraction1-16*fraction2*eps1r**2
    *fraction1+4*eps2c**2*fraction2**2*eps1c**2-16*fraction3*eps1r**2*fraction1+4*
    eps1c**2*fraction1**2*eps3r**2+4*eps1r*eps2c**2*eps3r**2+4*eps1r**2*eps2r*eps3c**2
    +4*eps1r**2*eps2c**2*eps3r+eps1r**2*eps2r**2*eps3c**2+4*eps1r*eps2r**2*eps3r**
    2+eps1r**2*eps2c**2*eps3c**2+16*eps1r*eps2r*eps3c**2+eps1r**2*eps2c**2*eps3r**2+
    64*eps1r*eps2r*eps3r+eps1r**2*eps2r**2*eps3r**2+16*eps1r*eps2r**2*eps3r+eps1c**2
    *eps2r**2*eps3r**2+16*eps1r*eps2c**2*eps3r+4*eps1r**2*eps2r*eps3r**2+4*eps1r*
    eps2r**2*eps3c**2+4*eps1c**2*eps2r**2*eps3r+16*eps1r*eps2r*eps3r**2+4*eps1r**2*
    eps2r**2*eps3r+4*eps1r*eps2c**2*eps3c**2+eps1c**2*eps2r**2*eps3c**2+4*eps1c**2*
    eps2r*eps3r**2+16*eps1c**2*eps2r*eps3r+4*eps1c**2*eps2r*eps3c**2+eps1c**2*eps2c**2
    *eps3r**2+4*eps1c**2*eps2c**2*eps3r+eps1c**2*eps2c**2*eps3c**2-8*fraction1*
    eps1r**2*eps2r**2-32*fraction1*eps1r**2*eps2r-8*fraction1*eps1r**2*eps2c**2-8*
    fraction1*eps1r*eps2r**2-32*fraction1*eps1r*eps2r-8*fraction1*eps1r*eps2c**2-8*
    fraction1*eps1c**2*eps2r**2-32*fraction1*eps1c**2*eps2r-8*fraction1*eps1c**2*
    eps2c**2-8*fraction1*eps1r*eps3r**2+2*fraction2*eps3r**2*fraction1*eps1c**2*
    eps2r**2+2*fraction2*eps3r**2*fraction1*eps1c**2*eps2c**2+2*fraction2*eps3c**2*
    fraction1*eps1r**2*eps2c**2+2*fraction2*eps3c**2*fraction1*eps1c**2*eps2r**2+2*
    fraction2*eps3c**2*fraction1*eps1c**2*eps2c**2+2*fraction1*eps1r*eps2c**2*eps3r**2
    *fraction2+2*fraction1*eps1r*eps2c**2*eps3r**2*fraction3+8*fraction2*eps3r*
    fraction1*eps1r**2*eps2c**2+8*fraction2*eps3r*fraction1*eps1c**2*eps2r**2+8*
    fraction2*eps3r*fraction1*eps1c**2*eps2c**2+2*fraction1*eps1r*eps2c**2*eps3c**2*
    fraction2+2*fraction1*eps1r*eps2c**2*eps3c**2*fraction3+2*fraction3*eps1r**2*
    eps2c**2*eps3r*fraction1+8*fraction3*eps1r**2*eps2r*eps3r*fraction1+2*fraction2*
    eps1c**2*eps2r*eps3c**2*fraction1+2*fraction2*eps1c**2*eps2r*eps3r**2*fraction1+
    8*fraction2*eps1c**2*eps2r*eps3r*fraction1+2*fraction2*eps1r*eps2r*eps3r**2*
    fraction1+2*fraction2*eps1r**2*eps2r*eps3r**2*fraction1+2*fraction2*eps1r*eps2r*
    eps3c**2*fraction1+8*fraction2*eps1r*eps2r*eps3r*fraction1+2*fraction2*eps1r**2*
    eps2r*eps3c**2*fraction1+8*fraction2*eps1r**2*eps2r*eps3r*fraction1+8*fraction3*
    fraction2*eps1r*eps2c**2*eps3r**2+2*fraction3*fraction2*eps1r**2*eps2r**2*eps3c**2
    -4*fraction2*eps2r*eps3c**2*fraction1-32*fraction1**2*eps1r+16*fraction1**2*
    eps2r+16*fraction1**2*eps3r+4*fraction1**2*eps2r**2+4*fraction1**2*eps3r**2+4*
    fraction1**2*eps2c**2+4*fraction1**2*eps3c**2+4*fraction2**2*eps3r**2*eps2r**2-
    32*fraction1*eps1r*eps3r-8*fraction1*eps1r*eps3c**2+16*fraction1*eps2r*eps3r**2+
    64*fraction1*eps2r*eps3r+16*fraction1*eps2r*eps3c**2-8*fraction1*eps1r**2*eps3r**2
    -32*fraction1*eps1r**2*eps3r-8*fraction1*eps1r**2*eps3c**2+4*fraction1*eps2r**2
    *eps3r**2+16*fraction1*eps2r**2*eps3r+4*fraction1*eps2r**2*eps3c**2-8*fraction1
    *eps1c**2*eps3r**2-32*fraction1*eps1c**2*eps3r-8*fraction1*eps1c**2*eps3c**2+4*
    fraction1*eps2c**2*eps3r**2+16*fraction1*eps2c**2*eps3r+4*fraction1*eps2c**2*
    eps3c**2-8*fraction2*eps1r**2*eps2r**2-8*fraction2*eps1r**2*eps2r-8*fraction2*
    eps1r**2*eps2c**2-32*fraction2*eps1r*eps2r**2-32*fraction2*eps1r*eps2r-32*
    fraction2*eps1r*eps2c**2-8*fraction2*eps1c**2*eps2r**2-8*fraction2*eps1c**2*
    eps2r-8*fraction2*eps1c**2*eps2c**2+16*fraction2*eps1r*eps3r**2+64*fraction2*
    eps1r*eps3r+16*fraction2*eps1r*eps3c**2-8*fraction2*eps2r*eps3r**2-32*fraction2*
    eps2r*eps3r-8*fraction2*eps2r*eps3c**2+4*fraction2*eps1r**2*eps3r**2+16*
    fraction2*eps1r**2*eps3r+4*fraction2*eps1r**2*eps3c**2-8*fraction2*eps2r**2*
    eps3r**2-32*fraction2*eps2r**2*eps3r-8*fraction2*eps2r**2*eps3c**2+4*fraction2*
    eps1c**2*eps3r**2+16*fraction2*eps1c**2*eps3r+4*fraction2*eps1c**2*eps3c**2-8*
    fraction2*eps2c**2*eps3r**2-32*fraction2*eps2c**2*eps3r-8*fraction2*eps2c**2*
    eps3c**2+4*fraction3*eps1r**2*eps2r**2+16*fraction3*eps1r**2*eps2r+4*fraction3*
    eps1r**2*eps2c**2+16*fraction3*eps1r*eps2r**2+64*fraction3*eps1r*eps2r+16*
    fraction3*eps1r*eps2c**2+4*fraction3*eps1c**2*eps2r**2+16*fraction3*eps1c**2*
    eps2r+4*fraction3*eps1c**2*eps2c**2-32*fraction3*eps1r*eps3r**2-32*fraction3*
    eps1r*eps3r-32*fraction3*eps1r*eps3c**2-32*fraction3*eps2r*eps3r**2-32*fraction3
    *eps2r*eps3r-32*fraction3*eps2r*eps3c**2-8*fraction3*eps1r**2*eps3r**2-8*
    fraction3*eps1r**2*eps3r-8*fraction3*eps1r**2*eps3c**2-8*fraction3*eps2r**2*
    eps3r**2-8*fraction3*eps2r**2*eps3r-8*fraction3*eps2r**2*eps3c**2-8*fraction3*
    eps1c**2*eps3r**2-8*fraction3*eps1c**2*eps3r-8*fraction3*eps1c**2*eps3c**2-8*
    fraction3*eps2c**2*eps3r**2-8*fraction3*eps2c**2*eps3r-8*fraction3*eps2c**2*
    eps3c**2+4*fraction3**2*eps2c**2*eps3r**2+4*fraction3**2*eps1r**2*eps3r**2+4*
    fraction3**2*eps2r**2*eps3r**2+4*fraction3**2*eps1c**2*eps3r**2+16*fraction3**2*
    eps1r*eps3r**2+16*fraction3**2*eps2r*eps3r**2-16*fraction3*eps3r*fraction1+8*
    fraction3*fraction2*eps1r*eps2r**2*eps3r**2+2*fraction3*fraction2*eps1r**2*eps2c
    **2*eps3r**2+2*fraction3*fraction2*eps1r**2*eps2r**2*eps3r**2+2*fraction3*
    fraction2*eps1c**2*eps2r**2*eps3r**2+8*fraction3*eps2r*fraction1*eps1c**2*eps3c**2
    +2*fraction3*eps2r**2*fraction1*eps1c**2*eps3c**2+2*fraction3*eps1c**2*
    fraction1*eps2c**2*eps3r**2+2*fraction3*eps1c**2*fraction1*eps2c**2*eps3c**2+8*
    fraction2*eps1r*eps2r*eps3r**2*fraction3+2*fraction2*eps1c**2*eps2r*eps3r**2*
    fraction3+2*fraction2*eps1c**2*eps2r*eps3r*fraction3+2*fraction2*eps1c**2*eps2r*
    eps3c**2*fraction3+8*fraction3*eps1r**2*fraction1*eps2r*eps3r**2+8*fraction3*
    eps1c**2*fraction1*eps2r*eps3r**2+8*fraction3*fraction2*eps1r*eps2c**2*eps3c**2+
    2*fraction3*fraction2*eps1c**2*eps2r**2*eps3c**2+2*fraction3*fraction2*eps1c**2*
    eps2c**2*eps3r**2+2*fraction3*fraction2*eps1c**2*eps2c**2*eps3c**2+8*fraction3*
    eps2r**2*fraction2*eps1r*eps3c**2+2*fraction3*eps1r**2*eps2r**2*fraction1*eps3r**2
    +2*fraction3*eps1r**2*eps2r**2*fraction1*eps3c**2+8*fraction3*eps1r**2*eps2r*
    fraction1*eps3c**2+2*fraction3*eps1r**2*eps2c**2*fraction1*eps3r**2+2*fraction3*
    eps1r**2*eps2c**2*fraction1*eps3c**2+2*fraction3*eps1c**2*eps2r**2*fraction1*
    eps3r**2+8*fraction3*eps1r*eps2r**2*eps3r*fraction2+8*fraction3*eps1r*eps2c**2*
    eps3r*fraction2+2*fraction3*eps1c**2*eps2r**2*eps3r*fraction2-16*fraction3*eps1r
    *fraction1-16*fraction2*eps1c**2*fraction1+8*fraction2*eps3c**2*fraction1-16*
    fraction3*eps1c**2*fraction1+8*fraction3*eps2c**2*fraction1-16*fraction2*eps1r*
    fraction1+4*eps2c**2*fraction2**2*eps1r**2+4*eps1c**2*fraction1**2*eps2c**2+16*
    eps3c**2*fraction3**2*eps1r+4*eps1c**2*fraction1**2*eps3c**2+fraction1**2*eps2c**2
    *eps3c**2+fraction1**2*eps2c**2*eps3r**2+4*fraction1**2*eps2c**2*eps3r+
    fraction1**2*eps2r**2*eps3c**2+4*fraction1**2*eps2r**2*eps3r+fraction1**2*eps2r**2
    *eps3r**2+16*fraction1**2*eps2r*eps3r+4*fraction1**2*eps2r*eps3c**2+4*
    fraction1**2*eps2r*eps3r**2+16*eps2c**2*fraction2**2*eps3r+16*eps1c**2*fraction1
    **2*eps3r-16*fraction2*eps2r*fraction3-16*fraction2*eps3r*fraction3+8*fraction2*
    eps1r**2*fraction3-16*fraction2*eps3r**2*fraction3+16*fraction1**2*eps1r**2*
    eps2r+4*fraction1**2*eps1r**2*eps2c**2+4*eps3c**2*fraction3**2*eps2c**2+4*
    fraction1**2*eps1r**2*eps2r**2-16*fraction1*eps2r**2*fraction2-16*fraction1*
    eps3r**2*fraction3-16*fraction1*eps2c**2*fraction2-16*fraction1*eps3c**2*
    fraction3+16*eps3c**2*fraction3**2*eps2r+32*fraction2*eps3r*fraction1+4*
    fraction2**2*eps1r**2*eps2r**2+4*fraction3**2*eps1r**2*eps2r+fraction3**2*eps1r**2
    *eps2c**2+4*fraction3**2*eps1r*eps2r**2-4*fraction2*eps1r*eps3c**2*fraction1-4
    *fraction2*eps2r*eps3r**2*fraction1-16*fraction2*eps2r*eps3r*fraction1-32*
    fraction3**2*eps1r*eps2r*eps3r-8*fraction3**2*eps1r*eps2r**2*eps3r-8*fraction3**
    2*eps1r*eps2c**2*eps3r-4*fraction2*eps1c**2*eps3r*fraction3-4*fraction2*eps1c**2
    *eps3c**2*fraction3+4*eps1c**2*fraction1**2*eps2r**2*eps3r-4*fraction2*eps1r**2*
    eps3r**2*fraction3-4*fraction2*eps1r**2*eps3r*fraction3-4*fraction2*eps1r**2*
    eps3c**2*fraction3-4*fraction2*eps1c**2*eps3r**2*fraction3+8*fraction2*eps1r**2*
    eps2r*fraction1+8*fraction2*eps1r*eps2r*fraction1+8*fraction2*eps2r*eps3r**2*
    fraction3+8*fraction2*eps2r*eps3r*fraction3+8*fraction2*eps2r*eps3c**2*fraction3
    -8*fraction1**2*eps1r*eps2r*eps3r**2-2*fraction1**2*eps1r*eps2c**2*eps3c**2+16*
    eps1c**2*fraction1**2*eps2r*eps3r-2*fraction1**2*eps1r*eps2r**2*eps3c**2-32*
    fraction1**2*eps1r*eps2r*eps3r-8*fraction1**2*eps1r*eps2r**2*eps3r-8*fraction1**
    2*eps1r*eps2c**2*eps3r-2*fraction1**2*eps1r*eps2c**2*eps3r**2-2*fraction1**2*
    eps1r*eps2r**2*eps3r**2-8*fraction1**2*eps1r*eps2r*eps3c**2-4*fraction3*eps2r**2
    *eps3r*fraction1+8*fraction3*eps1c**2*eps3r*fraction1-4*fraction3*eps2c**2*eps3r
    *fraction1-16*fraction3*eps1c**2*eps2r*fraction1-4*fraction3*eps1c**2*eps2c**2*
    fraction1+8*fraction3*eps1r*eps3r*fraction1-16*fraction3*eps2r*eps3r*fraction1+8
    *fraction3*eps1r**2*eps3r*fraction1-4*fraction3*eps1c**2*eps2r**2*fraction1+
    eps2c**2*fraction2**2*eps1c**2*eps3r**2-2*fraction2**2*eps1r**2*eps2r*eps3c**2-8
    *fraction2**2*eps1r**2*eps2r*eps3r+4*eps1c**2*fraction1**2*eps2r*eps3r**2-16*
    fraction1*eps2r*eps3r**2*fraction3-16*fraction1*eps2r*eps3c**2*fraction3+8*
    fraction3*eps2c**2*eps3r*fraction2+8*fraction1*eps1r*eps3r**2*fraction3+8*
    fraction1*eps1r*eps3c**2*fraction3+4*eps3c**2*fraction3**2*eps1c**2*eps2r+8*
    fraction2*fraction1*eps1r**2*eps2c**2+8*fraction2*fraction1*eps1r**2*eps2r**2+
    eps3c**2*fraction3**2*eps1c**2*eps2r**2+4*eps2c**2*fraction2**2*eps1r*eps3c**2-2
    *fraction2**2*eps1c**2*eps2r*eps3r**2-8*fraction2**2*eps1c**2*eps2r*eps3r-2*
    fraction2**2*eps1c**2*eps2r*eps3c**2+4*eps1c**2*fraction1**2*eps2c**2*eps3r+
    eps1c**2*fraction1**2*eps2c**2*eps3c**2-4*fraction2*eps1r**2*eps2r*fraction3-16*
    fraction2*eps1r*eps2r*fraction3-2*fraction3*eps1r**2*eps2c**2*eps3c**2-32*
    fraction1*eps1r**2*eps2r*eps3r+16*fraction3**2*eps1r*eps2r+4*fraction3**2*eps1r*
    eps2c**2+fraction3**2*eps1c**2*eps2r**2+4*fraction3**2*eps1c**2*eps2r+fraction3**2
    *eps1c**2*eps2c**2-32*fraction3**2*eps1r*eps3r-32*fraction3**2*eps2r*eps3r+
    fraction2**2*eps1r**2*eps3r**2+4*fraction2**2*eps1r**2*eps3r+fraction2**2*eps1r**2
    *eps3c**2+fraction2**2*eps1c**2*eps3r**2+4*fraction2**2*eps1c**2*eps3r+
    fraction2**2*eps1c**2*eps3c**2+fraction3**2*eps1r**2*eps2r**2-8*fraction2**2*
    eps1r**2*eps2r-32*fraction2**2*eps1r*eps2r-8*fraction2**2*eps1c**2*eps2r+4*
    fraction2**2*eps1r*eps3r**2+16*fraction2**2*eps1r*eps3r+4*fraction2**2*eps1r*
    eps3c**2-8*fraction2**2*eps2r*eps3r**2-32*fraction2**2*eps2r*eps3r-8*fraction2**
    2*eps2r*eps3c**2+16*eps1c**2*fraction1**2*eps2r-16*fraction3*eps2r**2*fraction2-
    16*fraction3*eps2c**2*fraction2+4*fraction1**2*eps1r**2*eps3r**2+4*eps3c**2*
    fraction3**2*eps1r**2+4*fraction1**2*eps1r**2*eps3c**2+16*fraction1**2*eps1r**2*
    eps3r+16*fraction2**2*eps1r*eps2r**2+4*eps1c**2*fraction1**2*eps2r**2+16*
    fraction2**2*eps3r*eps2r**2+32*fraction2*eps1r*fraction3-8*fraction3**2*eps2c**2
    *eps3r+8*fraction2*eps1c**2*fraction3-16*fraction2*eps3c**2*fraction3+32*
    fraction3*eps2r*fraction1-2*fraction1*eps1r*eps2c**2*eps3r**2-8*fraction1*eps1r**2
    *eps2r*eps3c**2-8*fraction1*eps1r**2*eps2c**2*eps3r-2*fraction1*eps1r**2*eps2r
    **2*eps3c**2-2*fraction1*eps1r*eps2r**2*eps3r**2-2*fraction1*eps1r**2*eps2c**2*
    eps3c**2-8*fraction1*eps1r*eps2r*eps3c**2-2*fraction1*eps1r**2*eps2c**2*eps3r**2
    -32*fraction1*eps1r*eps2r*eps3r-2*fraction1*eps1r**2*eps2r**2*eps3r**2-8*
    fraction1*eps1r*eps2r**2*eps3r-2*fraction1*eps1c**2*eps2r**2*eps3r**2-8*
    fraction1*eps1r*eps2c**2*eps3r-8*fraction1*eps1r**2*eps2r*eps3r**2-2*fraction1*
    eps1r*eps2r**2*eps3c**2-8*fraction1*eps1c**2*eps2r**2*eps3r-8*fraction1*eps1r*
    eps2r*eps3r**2-8*fraction1*eps1r**2*eps2r**2*eps3r-2*fraction1*eps1r*eps2c**2*
    eps3c**2-2*fraction1*eps1c**2*eps2r**2*eps3c**2-8*fraction1*eps1c**2*eps2r*eps3r
    **2-32*fraction1*eps1c**2*eps2r*eps3r-8*fraction1*eps1c**2*eps2r*eps3c**2-2*
    fraction1*eps1c**2*eps2c**2*eps3r**2-8*fraction1*eps1c**2*eps2c**2*eps3r-2*
    fraction1*eps1c**2*eps2c**2*eps3c**2-8*fraction2*eps1r**2*eps2r*eps3r-8*
    fraction2*eps1r*eps2c**2*eps3r**2-2*fraction2*eps1r**2*eps2r*eps3c**2-8*
    fraction2*eps1r**2*eps2c**2*eps3r-8*fraction3**2*eps1r**2*eps3r)
  
  epsilon_effective = epsilon_effective_real + 1.0j*epsilon_effective_imag
  return epsilon_effective


## Lorentz-Lorenz method: see the Maxwell-Garnett model ( MaxwellGarnett2() )
def LorentzLorenz2(eps1, eps2, fraction):
  return MaxwellGarnett2(eps1, eps2, fraction)

## Lorentz-Lorenz method: see the Maxwell-Garnett model ( MaxwellGarnett3() )
def LorentzLorenz3(eps1, eps2, eps3, fraction1, fraction2):
  return MaxwellGarnett3(eps1, esp2, eps3, fraction1, fraction2)

## Bruggeman method to compute the effective dielectric permittivity
# @param eps1: dielectric permittivity of the ambient medium
# @param eps2: dielectric permittivity of the nanoparticle medium
# @param fraction2: volume fraction ratio of the eps2 medium inside the medium1. WARNING: if medium1 is a film, careful to compute the right geometry dependent volume fraction. 
def Bruggeman2(eps1, eps2, fraction2): #{{{
   # We use Bruggeman theory for 2 media, and provide both analytical solutions. 
   SquareRootContent=4*eps1**2 + 4*eps1*eps2 -12*fraction2*eps1**2 + 18*eps1*fraction2*eps2 + eps2**2 - 6*fraction2*eps2**2 + 9*fraction2**2*eps1**2 - 18*fraction2**2*eps1*eps2 + 9*fraction2**2*eps2**2
   FirstTerm = 0.5*eps1 - 0.25*eps2 - 0.75*fraction2*eps1 + 0.75*fraction2*eps2 
   sol1 = FirstTerm - 0.25 * np.sqrt(SquareRootContent)
   sol2 = FirstTerm + 0.25 * np.sqrt(SquareRootContent)
   return sol1, sol2
#}}}
         
Drude_Cr        = np.vectorize(Drude_Cr)
Drude_Ti        = np.vectorize(Drude_Ti)
EpsilonToIndex  = np.vectorize(EpsilonToIndex)
MaxwellGarnett2 = np.vectorize(MaxwellGarnett2)
MaxwellGarnett3 = np.vectorize(MaxwellGarnett3)
LorentzLorenz2  = np.vectorize(LorentzLorenz2)
LorentzLorenz3  = np.vectorize(LorentzLorenz3)
Bruggeman2      = np.vectorize(Bruggeman2)

#print "Attempt to use Maxwell-Garnett."

#print "Import Material 1"
#print "Import Material 2"
#print "Define a wavelength."
#print "Define a ratio of Material 1 (Material 2 = 1 - ratio)."

#MaxwellGarnett2(eps1, eps2, ratio)

## Computes Reflectivity in 3-material thin film configuration, where media 1 and 3 and half-infinite. 
# @param wavelength: wavelength (in meters) of the indicent photon
# @param eps123: complex dielectric permittivity of media 1 2 and 3
# @param thickness2: thickness of medium 2
def BiLayerReflectivity(wavelength, eps1, eps2, eps3, thickness2): #{{{
  phi2 = np.multiply(2.*np.pi*thickness2/wavelength, (np.sqrt(eps2))) 
  ##TODO: For the non-normal absorption, use Kovalenko formula of refraction. 
  r12  = ComplexReflectivity(eps1, eps2)
  r23  = ComplexReflectivity(eps2, eps3)
  r13  = (r12 + r23*np.exp(2.j*phi2)) / (1.+r12*r23*np.exp(2.j*phi2))
  return r13*np.conjugate(r13)
#}}}

## Compute transmission from a 3-material thin film configuration, where media 1 and 3 are semi-infinite. Formula originates from Stenzel book, page 108. Warning: formula from Born and Wolf may contain mistakes. It could not be validated for perfect dielectrics. 
# @param wavelength: wavelength (in meters) of the indicent photon
# @param eps123: complex dielectric permittivity of media 1 2 and 3
# @param thickness2: thickness of medium 2
def BiLayerTransmission(wavelength, eps1, eps2, eps3, thickness2): #{{{
   n3  = np.sqrt(eps3); n1 = np.sqrt(eps1)
   #r12  = ComplexReflectivity(eps1, eps2);  
   r21 = ComplexReflectivity(eps2, eps1)
   r23 = ComplexReflectivity(eps2, eps3); #r32 = ComplexReflectivity(eps3, eps2)                                                                            
   t12 = ComplexTransmission(eps1, eps2); #t21 = ComplexTransmission(eps2, eps1)
   t23 = ComplexTransmission(eps2, eps3); #t32 = ComplexTransmission(eps3, eps2)
   delta = 2*np.pi/wavelength * thickness2 * np.sqrt(eps2)
   t13 = t12 * t23 * np.exp(1j*delta)/(1.-r21*r23*np.exp(2j*delta))
   T13 = np.real(n3)/np.real(n1) * t13 * np.conjugate(t13) #total transmittivity
   return T13
#}}}                                                                           
BiLayerReflectivity = np.vectorize(BiLayerReflectivity)
BiLayerTransmission = np.vectorize(BiLayerTransmission)

