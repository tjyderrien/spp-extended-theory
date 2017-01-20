#!/usr/bin/env python
#-*- coding: utf-8 -*-
## @package libMath
# Defines some mathematical functions not yet defined in SciPy/NumPy. 

import numpy as np

## Define a complex sign function for Python
def csgn(x, y): #{{{
  # returns the complex sign, as in maple
  if(x.real > 0e0 or x.real == 0e0 and x.imag > 0e0):
    result = 1E0
  else:
    if(x.real < 0e0 or x.real == 0 and x.imag < 0e0):
      result=-1e0
    else: 
      result = -1E99
      print "csgn: Exception case, to be solved."
#}}}

##Complex derivative a real-valued function by a complex-number
#  Input:
#    f: z->f(z)
#    z: z complex-valued numbers
#  Output:
#    df/dz according to Wirtinger complex derivatives formula
#  Careful: Wirtinger formula contain a 1/2 usually, but it was removed here to obtain consistent results with Hohenau lifetime and Raether lifetime. 
def RealDerivativeByComplex(f,z):
  
  return (np.diff(f)/np.diff(z.real) - 1j*(np.diff(f)/np.diff(z.imag))) #2 x original

## Swap to variables
def Swap(eps1, eps2):#{{{
  eps3 = eps1
  eps1 = eps2
  eps2 = eps3
  del eps3
  return(eps1, eps2) 
#}}}

Swap = np.vectorize(Swap)