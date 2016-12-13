#!/usr/bin/env python
#-*- coding: utf-8 -*-
## @package libMaterials 
# Functions describing materials and their interaction with light. 

import numpy as np
import cmath
from libLaser import *

# OPTICAL FUNCTIONS

## Return the value of dielectric function based on simplified Drude model
# Input:
# @param wavelength (float)
# @param ne (float)
# @epsilon (complex): dielectric permittivity under wavelength, without excitation
# @nu (float): collision frequency
# Output: complex-valued dielectric permittivity
def Drude(wavelength, ne, epsilon, nu, meff=1.0):#{{{
  omegap2=ne * e**2 / (m_e * meff * epsilon_0)
  omega=2.0*pi*c/wavelength
  return epsilon - omegap2/(omega*omega) * 1e0/(1e0+1e0j*nu/omega)
#}}}

## Return Fresnel reflectivity 
# Input:
#   eps1: complex-valued permittivity 1+j0
#   eps2: idem, for medium2
# Output: 
#   interface reflectivity (float) R
def reflectivity(eps1, eps2):#{{{
  R=abs(((eps1**0.5e0-eps2**0.5e0)/(eps1**0.5e0+eps2**0.5e0))**2)
  return R
#}}}

## Returns the complex refractive index
def EpsilonToIndex(eps):
  return cmath.sqrt(eps)

## Returns the complex permittivity from optical index
def IndexToEpsilon(n):
  return n*n

## Calculate effective dielectric permittivity for 2-mediums using Maxwell-Garnett method
# @param eps1: dielectric permittivity (epsilon <complex>) of first medium
# @param eps2: dielectric permittivity (epsilon <complex>) of second medium
# @param fraction: fraction of epsilon2 mixed with (1.-fraction)*epsilon1 medium
def MaxwellGarnett2(eps1, eps2, fraction):
  eps1r = eps1.real; eps1c = eps1.imag
  eps2r = eps2.real; eps2c = eps2.imag
  # Space for optimization is not so big: 23*storage+23*assignments+50*multiplications+37*additions+divisions. 
  # Whereas direct writing uses: 37*additions+104*multiplications+divisions+assignments
  epsilon_effective_real = -(1.*(2.*eps1c**2*fraction**2-2.*fraction*eps2r**2-4.*eps1r+eps1r*fraction*eps2c**2-1.*eps1c**2*fraction*eps2r+4.*eps1r*fraction-1.*eps1r**2*fraction*eps2r+eps1r*fraction*eps2r**2+4.*eps1r*fraction*eps2r-1.*eps1r*eps2c**2-2.*eps1c**2*fraction+4.*eps1c*eps2c*fraction+2.*eps2c**2*fraction**2-2.*eps1r**2*fraction-1.*eps1r*eps2r**2-4.*eps1r*eps2r+2.*eps1r**2*fraction**2-4.*fraction*eps2r-2.*fraction*eps2c**2-4.*eps2c*fraction**2*eps1c+2.*fraction**2*eps2r**2-4.*fraction**2*eps2r*eps1r))/(2.*eps1r*fraction*eps2r+4.+4.*eps1r*fraction-2.*fraction*eps2r**2-4.*fraction*eps2r-2.*fraction*eps2c**2+eps2r**2+4.*eps2r+eps2c**2+2.*eps1c*eps2c*fraction+eps1r**2*fraction**2+fraction**2*eps2r**2-2.*fraction**2*eps2r*eps1r-2.*eps2c*fraction**2*eps1c+eps2c**2*fraction**2+eps1c**2*fraction**2)
  # This version uses: 26*additions+68*multiplications+divisions+assignments
  # Optimization would use: 10*storage+10*assignments+46*multiplications+26*additions+divisions
  epsilon_effective_imag = (4.*eps2c*fraction-1.*eps1c*fraction*eps2c**2+4.*eps1c*eps2r-4.*eps1c*fraction*eps2r+4.*eps2c*fraction*eps1r+eps2c*fraction*eps1c**2-4.*eps1c*fraction+eps2c*fraction*eps1r**2+4.*eps1c+eps1c*eps2c**2+eps1c*eps2r**2-1.*eps1c*fraction*eps2r**2)/(2.*eps1r*fraction*eps2r+4.+4.*eps1r*fraction-2.*fraction*eps2r**2-4.*fraction*eps2r-2.*fraction*eps2c**2+eps2r**2+4.*eps2r+eps2c**2+2.*eps1c*eps2c*fraction+eps1r**2*fraction**2+fraction**2*eps2r**2-2.*fraction**2*eps2r*eps1r-2.*eps2c*fraction**2*eps1c+eps2c**2*fraction**2+eps1c**2*fraction**2)
  
  epsilon_effective = epsilon_effective_real + 1.0j*epsilon_effective_imag
  return epsilon_effective

EpsilonToIndex = np.vectorize(EpsilonToIndex)
MaxwellGarnett2 = np.vectorize(MaxwellGarnett2)