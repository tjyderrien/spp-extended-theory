#!/usr/bin/env python
#-*- coding: utf-8 -*-
## @package libMaterials 
# Functions describing materials and their interaction with light. 

import numpy as np
import cmath

# OPTICAL FUNCTIONS
def Drude(wavelength, ne, epsilon, nu):#{{{
  """Return the value of dielectric function based on simplified Drude model
  Input:
    wavelength (float)
    ne (float)
    epsilon (complex): dielectric permittivity under wavelength, without excitation
    nu (float): collision frequency
  Output: complex-valued dielectric permittivity
  """
  omegap2=ne * e**2 / (m_e * meffe * epsilon_0)
  omega=2.0*pi*c/wavelength
  return epsilon - omegap2/(omega*omega) * 1e0/(1e0+1e0j*nu/omega)
#}}}

def reflectivity(eps1, eps2):#{{{
  """Return Fresnel reflectivity 
  Input:
    eps1: complex-valued permittivity 1+j0
    eps2: idem, for medium2
  Output: 
    interface reflectivity (float) R
  """
  R=abs(((eps1**0.5e0-eps2**0.5e0)/(eps1**0.5e0+eps2**0.5e0))**2)
  return R
#}}}

def EpsilonToIndex(eps):
  #returns the complex refractive index
  return cmath.sqrt(eps)

def IndexToEpsilon(n):
  #returns the complex permittivity from optical index
  return n*n

EpsilonToIndex = np.vectorize(EpsilonToIndex)