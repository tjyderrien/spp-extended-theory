#!/usr/bin/env python
#-*- coding: utf-8 -*-
## @package Keldysh
# This module aims to calculate the density of excited electrons as function of laser parameters. 
# Two types of usage are planned :
# * Generating tables to use directly into simulation codes
# * Outputing density in certain conditions. 

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

def IonizationRate(Keldysh1, Keldysh2, KeldyshFunctionResult, Ueff, wavelength):
  omegaLaser = 2.*pi*c/wavelength
  IonizationRate=2.*omegaLaser/(9.*pi)*((omegaLaser*m_e)/(hbar*Keldysh1))**(1.5)*KeldyshFunctionResult*np.exp(-pi*np.trunc(Ueff/hbar/omegaLaser+1)*((ellipk(Keldysh1**2)-ellipe(Keldysh1**2))/(ellipe(Keldysh2**2))))
	  
  return IonizationRate

print ""
print "** Welcome to SPP-extended-theory suite."
print "Author(s): T.J.-Y. Derrien"
print ""
print "** Loading Keldysh module [Keldysh, Sov. J. Exp. Th. Phys. 47, 5 (1964)]..."
print "** TODO: Correct formula using [Gruzdev, Optical Engineering 53, 122515 (2014)"
Egap = 1.12e0*e; meff=0.18e0; Efield=1E9; wavelength = 800e-9

## This section can consume GB of RAM. Do not use. #{{{ 
#print ""
#print "** Vectorizing functions..."
gammaKeldysh = np.vectorize(gammaKeldysh)
Keldysh1 = np.vectorize(Keldysh1)
Keldysh2 = np.vectorize(Keldysh2)
EffectiveGap=np.vectorize(EffectiveGap)
KeldyshFunction=np.vectorize(KeldyshFunction)
IonizationRate=np.vectorize(IonizationRate)
#}}}

print ""
print "** Generating mesh..."
Efield = 1E8*np.arange(1,100,1)

print ""
print "Computing gamma..."
gamma = gammaKeldysh(Egap, meff, Efield, wavelength) #valid
#print gamma

print "Computing Keldysh1, Keldysh2..."
k1 = Keldysh1(gamma); k2 = Keldysh2(gamma) #valid
EgapEff = EffectiveGap(Egap, k1, k2) #Warning: scipy.special.ellipe (https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.ellipe.html#scipy.special.ellipe) uses a different convention than Maple, Wikipedia or mpmath.
order = 10

print ""
print "** Info: "
print "Adiabadicity parameter = "+str(gamma)
print "Egap = "+str(Egap/e)+" eV, Ueff = "+str(EgapEff/e)+" eV."
print ""

KeldyshFunctionResult = KeldyshFunction( k1, k2, EgapEff, order, wavelength )
wPI = IonizationRate(k1, k2, KeldyshFunctionResult, EgapEff, wavelength)
print "w_PI until order "+str(order)+" = ", wPI

plt.figure()
plt.xlabel("Field (V/m)")
plt.ylabel("Excitation rate w_{PI} (m^{-3} s^{-1})")
plt.plot(Efield, wPI, label="w_{PI}")
#plt.loglog(Efield, gamma, label="gamma")
plt.show()