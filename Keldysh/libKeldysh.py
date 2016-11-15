#!/usr/bin/env python
#-*- coding: utf-8 -*-
## @package libKeldysh
# This module aims to calculate the density of excited electrons as function of laser parameters. 

# IMPORT LIBRARIES
import numpy as np
from numpy import genfromtxt, loadtxt, chararray
#from scipy.optimize import fsolve, root
from scipy.special import ellipk, ellipe, dawsn
#import cmath
import matplotlib as mp
import matplotlib.pyplot as plt
#from scipy.interpolate import InterpolatedUnivariateSpline
from matplotlib import rc
# from pylab import *
from scipy.constants import c, epsilon_0, mu_0, pi, e, m_e, h, hbar
#from matplotlib.legend_handler import HandlerLine2D
#import sys

## Computes the adiabadicity parameter
#@param gamma: Adiabadicity parameter (non-dimensional number)
# 
# * gammaKeldysh < 0.1: means tunneling effect is dominant, 
# * gammaKeldysh > 10: multi-photon excitation effect is dominant. 
def gammaKeldysh(Egap, meff, Efield, wavelength): #{{{
  
  #print Egap, meff, Efield	
  omegaLaser=2.*pi*c/wavelength
  if (Efield != 0e0):
    result = omegaLaser*np.sqrt(m_e*meff*Egap)/e/Efield
  else:
    print "gamma(): Divergence, as field equals = 0. Singular case of Keldysh functions. Should give w_PI = 0 then..."
  #print omegaLaser
  return result
# Numerically validated with comparison to Maple. 
#}}}

## Computes some intermediate quantity
def Keldysh1(gamma):
  return gamma/np.sqrt(1.0+gamma**2)

## Computes some intermediate quantity
def Keldysh2(gamma):
  return Keldysh1(gamma)/gamma

## Computes the effective gap for one material
# @param Egap: band gap of the transition (multi-photonic / tunnel transitions are DIRECT)
def EffectiveGap(Egap, k1, k2): #{{{
  if (k1 != 0):
    result = 2.0*Egap * ellipe(k2**2)/(pi*k1) #Warning: ellipe(x²) actually computes E(x). 
  else: 
    print "EffectiveGap(): Keldysh1 = 0."
    exit()
  return result
  # Validation: 
  # EffectiveGap(0, 0, 0) = Error. 
  # EffectiveGap(1., 1., 0) = 1
  # EffectiveGap(1.12*e, 1., 1.) = 0.71*e
  # EffectiveGap(Egap,k1,k2)/e = 1.12044049367 #Passed
#}}}

## Calculates the Dawson integral int(exp(y**2 - z**2), y=0..z)
# This is based on the Python library SciPy.special functions. 
# According to https://docs.scipy.org/doc/scipy-0.14.0/reference/generated/scipy.special.dawsn.html#scipy.special.dawsn
def DawsonIntegral(z): #{{{
  integral = dawsn(z)
  # Validation, compared with Maple. 
  # OK DawsonIntegral(0.)=0
  # OK DawsonIntegral(0.5) = 0.42443638350202229
  return integral
#}}}

def KeldyshFunction(Keldysh1, Keldysh2, Ueff, nmax, wavelength): #{{{
  omegaLaser = 2.*pi*c/wavelength #SI
  n=np.arange(0,nmax+1)
  #print n
  sumtable=np.exp(-pi*n*(ellipk(Keldysh1**2)-ellipe(Keldysh1**2))/ellipe(Keldysh2**2))*DawsonIntegral(pi*np.sqrt( ((2.0*np.trunc(Ueff/hbar/omegaLaser+1.))-2.0*Ueff/hbar/omegaLaser + n) / (2.0 * ellipk(Keldysh2**2)*ellipe(Keldysh2**2)) ) )
  #print "Effective gap: "+str(Ueff/e)+" eV."
  #print sumtable
  result = np.multiply(np.sqrt(pi/(2.*ellipk(Keldysh2**2))),np.sum(sumtable))
  return result
""" Validation
  KeldyshFunction(Keldysh1=0, Keldysh2=0, Ueff=0, nmax=0, wavelength=0): error. 
  KeldyshFunction(Keldysh1=0, Keldysh2=0, Ueff=0, nmax=0, wavelength=800e-9): 0.
  KeldyshFunction(Keldysh1=0, Keldysh2=0, Ueff=1.12*e, nmax=0, wavelength=800e-9): 0 
"""
#}}}

## Corrected Keldysh function according to Vitaly Gruzdev (see Ref in details). 
# Kane direct band gap structure
# This function is corrected according to Gruzdev, Opt. Eng. 53, 122515 (2014)
# Only factors of 2 are removed from inside the Keldysh integral, but inserted 
# as a proportion in IonizationRate_Gruzdev function. 
def KeldyshFunction_Gruzdev(Keldysh1, Keldysh2, Ueff, nmax, wavelength): #{{{
  omegaLaser = 2.*pi*c/wavelength #SI
  n=np.arange(0,nmax+1)
  #print n
  sumtable=np.exp(-pi*n*(ellipk(Keldysh1**2)-ellipe(Keldysh1**2))/ellipe(Keldysh2**2))*DawsonIntegral(pi*np.sqrt( ((np.trunc(Ueff/hbar/omegaLaser+1.))-Ueff/hbar/omegaLaser + n) / (2.0 * ellipk(Keldysh2**2)*ellipe(Keldysh2**2)) ) )
  #print "Effective gap: "+str(Ueff/e)+" eV."
  #print sumtable
  result = np.multiply(np.sqrt(pi/(2.*ellipk(Keldysh2**2))), np.sum(sumtable))
  return result

## The original Keldysh function for Kane direct band gap
# This function is known to contain mistakes. 
def IonizationRate(Keldysh1, Keldysh2, KeldyshFunctionResult, Ueff, wavelength):
  omegaLaser = 2.*pi*c/wavelength
  result = 2.*omegaLaser/(9.*pi)*((omegaLaser*m_e)/(hbar*Keldysh1))**(1.5)*KeldyshFunctionResult*np.exp(-pi*np.trunc(Ueff/hbar/omegaLaser+1)*((ellipk(Keldysh1**2)-ellipe(Keldysh1**2))/(ellipe(Keldysh2**2))))
  
  return result

## Corrected Keldysh photoionization probability according to Vitaly Gruzdev (see Ref in details). 
# Kane direct band gap structure
# This function is corrected according to Gruzdev, Opt. Eng. 53, 122515 (2014)
# Only factor of 2 was added outside the Keldysh integral, but removed
# inside. 
def IonizationRate_Gruzdev(Keldysh1, Keldysh2, KeldyshFunctionResult, Ueff, wavelength):
  omegaLaser = 2.*pi*c/wavelength
  result =2. * IonizationRate(Keldysh1, Keldysh2, KeldyshFunctionResult, Ueff, wavelength)
  #2.*omegaLaser/(9.*pi)*((omegaLaser*m_e)/(hbar*Keldysh1))**(1.5)*KeldyshFunctionResult*np.exp(-pi*np.trunc(Ueff/hbar/omegaLaser+1)*((ellipk(Keldysh1**2)-ellipe(Keldysh1**2))/(ellipe(Keldysh2**2))))
  
  return result

## Conversion between field and intensity (SI units)
# TODO: As it is absorbed field, it should multiplied by real(optical index)
# not exactly given by Keldysh theory, but rather by :
# 1. Value at rest (can be taken in "MaterialOpticalDatabase.dat" 
# 2. Value with excitation using a Drude model, associated with some effective mass and collision frequency. 
def FieldToIntensity(Field, permittivity=1):
  intensity = 0.5 * c * epsilon_0 * np.sqrt(permittivity) * Field**2
  return intensity.real

## This section can consume GB of RAM. Do not use. #{{{ 
#print ""
#print "** Vectorizing functions..."
FieldToIntensity = np.vectorize(FieldToIntensity)
gammaKeldysh = np.vectorize(gammaKeldysh)
Keldysh1 = np.vectorize(Keldysh1)
Keldysh2 = np.vectorize(Keldysh2)
EffectiveGap=np.vectorize(EffectiveGap)
KeldyshFunction=np.vectorize(KeldyshFunction)
KeldyshFunction_Gruzdev=np.vectorize(KeldyshFunction_Gruzdev)
IonizationRate_Gruzdev=np.vectorize(IonizationRate_Gruzdev)