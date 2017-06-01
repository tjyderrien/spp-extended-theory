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

## @package libKeldysh
# This module aims to calculate the density of excited electrons as function of laser parameters. 

# IMPORT LIBRARIES
import numpy as np
from numpy import genfromtxt, loadtxt, chararray
#from scipy.optimize import fsolve, root
from scipy.special import ellipk, ellipe, dawsn, factorial2, factorial, ellipkm1
#import cmath
import matplotlib as mp
import matplotlib.pyplot as plt
#from scipy.interpolate import InterpolatedUnivariateSpline
from matplotlib import rc
# from pylab import *
from scipy.constants import c, epsilon_0, mu_0, pi, e, m_e, h, hbar
#from matplotlib.legend_handler import HandlerLine2D
#import sys

from libDatabase import ExportToTxt
rc('font', **{'family':'serif', 'serif':['Georgia'], 'size':'16'})
rc('text', usetex=False)
mp.rcParams['legend.numpoints'] = 1

## Computes a step function
# Heaviside function
def step(x):
    return 1.0 * (x > 0.0)
step = np.vectorize(step)

## Computes the adiabadicity parameter
# @param gamma: Adiabadicity parameter (non-dimensional number)
# @param Egap: band gap energy (in Joules)
# @param meff: effective mass (in Arb. Units, as it is multiplied by electron mass)
# @param Efield: peak amplitude of the electric field (in V/m)
# @param wavelength: wavelength of the photon (in meters)
# 
# * gammaKeldysh < 0.1: means tunneling effect is dominant, 
# * gammaKeldysh > 10: multi-photon excitation effect is dominant. 
# 
# Warning: Don't use this function if Efield is too small (< 1 V/m), as it leads to divergence. 
#          Although problem should now be solved using higher precision numbers. 
def gammaKeldysh(Egap, meff, Efield, wavelength): #{{{
  #print Egap, meff, Efield
  omegaLaser=2.*pi*c/wavelength
  ErrorMessage=""
  # Validity limit
  if (Egap < hbar * omegaLaser ): 
    #TODO: this flow should be redirected to an error file. Stdout also goes into the variables. 
    ErrorMessage=ErrorMessage+"** Validity range error: the Keldysh model is not valid for linear absorption. INVALID RESULT...\n"
    ErrorMessage=ErrorMessage+"** Error details: "+str(int(wavelength*1E9))+" nm wavelength is too small for the gap "+str(float(Egap)/e)+".\n"
    #exit() #Avoid to quit, so that octopus still compare its results. 
  if (Efield > 1e-1): #if vectorial, then abs changed its meaning
    result = omegaLaser*np.sqrt(m_e*meff*Egap)/e/Efield
    #print Efield, result
  else:
    ErrorMessage=ErrorMessage+"gamma(): Divergence, as field equals = 0. Singular case of Keldysh functions. Should give w_PI = 0 then...\n"
    result = 1E9 #THIS VALUE IS ARBITRARY FOR A VERY SMALL FIELD. 
  #print omegaLaser
  return result
# Numerically validated with comparison to Maple. 
#}}}

## Computes some intermediate quantity, careful: long double precision.
def Keldysh1(gamma):
  value = np.float128(gamma) #K1(gamma) function has limit 1 when gamma > 5. Hence, we must compute k1(gamma) with a huge precision to stay out of unity. 
  return np.divide(value, np.sqrt( np.float128(1E0) + value * value) )

## Computes some intermediate quantity
def Keldysh2(gamma):
  value = np.float128(gamma) #idem about precision.
  try:
    result = np.divide( Keldysh1(value), value )
  except:
    result = 0e0
  return result

## Computes the effective gap for one material
# @param Egap: band gap of the transition (multi-photonic / tunnel transitions are DIRECT)
def EffectiveGap(Egap, k1, k2): #{{{
  # k1 = np.float64(k1); 
  k22 = np.float64(k2*k2) #reducing precision to call ellipe
  if (k1 != 0):
    result = 2.0*Egap * ellipe(k22)/(pi*k1) #Warning: ellipe(x²) actually computes E(x). 
  else: 
    print "EffectiveGap(): Singular error, Keldysh1 = 0."
    result = 0e0
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
  z2 = np.float64(z)
  integral = dawsn(z2)
  # Validation, compared with Maple. 
  # OK DawsonIntegral(0.)=0
  # OK DawsonIntegral(0.5) = 0.42443638350202229
  return integral
#}}}

def KeldyshFunction(Keldysh1, Keldysh2, Ueff, nmax, wavelength): #{{{
  omegaLaser = 2.*pi*c/wavelength #SI
  n_tab = np.arange(0,nmax+1)
 
  # print "Keldysh1 = "+str(Keldysh1)
  Keldysh11_128 = Keldysh1**2
  Keldysh11 = np.float64(Keldysh11_128)
  # print "--"
  # print "Keldysh11 = "+str(Keldysh11)
  Keldysh22 = np.float64(Keldysh2**2)
  # print "Keldysh11 = "+str(Keldysh11)
  distant_to_unity = 1E0 - Keldysh11_128
  if (distant_to_unity < 1E-320): #then it gonna crash for sure. 
    print "** Error on ellipk: argument 1 is singular. Distance to unit = "+str(distant_to_unity)+"Please increase precision on Keldysh1 or use ellipkm1 function (careful, argument IS not the same)."
  elif(distant_to_unity < 1E-10): 
    #threshold where functions ellipk and ellipkm1 give different values
    EllipticK1 = ellipkm1( np.float64(distant_to_unity) )
  else: #other cases, good for efficiency
    EllipticK1 = ellipk( Keldysh11 ) #inf if Keldysh11 = 1.  
  EllipticE1 = ellipe( Keldysh11 )
  EllipticE2 = ellipe( Keldysh22 )
  # print "EllipticK1 = "+str(EllipticK1)
  division = np.divide( EllipticK1 - EllipticE1 , EllipticE2 )
  # print "division= "+str(division)
  exponant = np.multiply( n_tab, division )
  sumtable = np.exp( - pi * exponant )
  sumtable = np.multiply(sumtable,DawsonIntegral(pi*np.sqrt( ((2.0*np.trunc(Ueff/hbar/omegaLaser+1.))-2.0*Ueff/hbar/omegaLaser + n_tab) / (2.0 * ellipk(Keldysh22)*ellipe(Keldysh22)) ) ) )
  #print "Effective gap: "+str(Ueff/e)+" eV."
  #print sumtable
  result = np.multiply(np.sqrt(pi/(2.*ellipk(Keldysh22))),np.sum(sumtable))
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
  Keldysh11_128 = Keldysh1**2
  Keldysh11 = np.float64(Keldysh11_128)
  # print "--"
  # print "Keldysh11 = "+str(Keldysh11)
  Keldysh22 = np.float64(Keldysh2**2)
  # print "Keldysh11 = "+str(Keldysh11)
  distant_to_unity = 1E0 - Keldysh11_128
  if (distant_to_unity < 1E-320): #then it gonna crash for sure. 
    print "** Error on ellipk: argument 1 is singular. Distance to unit = "+str(distant_to_unity)+"Please increase precision on Keldysh1 or use ellipkm1 function (careful, argument IS not the same)."
  elif(distant_to_unity < 1E-10):
    #threshold where functions ellipk and ellipkm1 give different values
    EllipticK1 = ellipkm1( np.float64(distant_to_unity) )
  else: #other cases, good for efficiency
    EllipticK1 = ellipk( Keldysh11 ) #inf if Keldysh11 = 1.i
  EllipticE1 = ellipe( Keldysh11 )
  EllipticE2 = ellipe( Keldysh22 )
  #print n
  sumtable=np.exp(-pi*n*(EllipticK1-EllipticE1)/EllipticE2)*DawsonIntegral(pi*np.sqrt( ((np.trunc(Ueff/hbar/omegaLaser+1.))-Ueff/hbar/omegaLaser + n) / (2.0 * ellipk(Keldysh22)*EllipticE2) ) )
  #print "Effective gap: "+str(Ueff/e)+" eV."
  #print sumtable
  result = np.multiply(np.sqrt(pi/(2.*ellipk(Keldysh22))), np.sum(sumtable))
  return result

## The original Keldysh function for Kane direct band gap
# This function is known to contain mistakes. 
def IonizationRate(Keldysh1, Keldysh2, KeldyshFunctionResult, Ueff, wavelength):
  omegaLaser = 2.*pi*c/wavelength
  Keldysh11 = np.float64(Keldysh1 * Keldysh1)
  Keldysh22 = np.float64(Keldysh2 * Keldysh2)
  try:
    result = 2.*omegaLaser/(9.*pi)*((omegaLaser*m_e)/(hbar*Keldysh1))**(1.5)*KeldyshFunctionResult*np.exp(-pi*np.trunc(Ueff/hbar/omegaLaser+1e0)*((ellipk(Keldysh11)-ellipe(Keldysh11))/(ellipe(Keldysh22))))
  except:
    result = 0e0
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
def FieldToIntensity(Field, permittivity=1e0):
  intensity = 0.5e0 * c * epsilon_0 * np.sqrt(permittivity) * Field * np.conjugate(Field)
  return intensity.real

## Converts intensity to electric field amplitude
def IntensityToField(intensity, permittivity=1.):
  Field = np.sqrt(2e0 * intensity / c / epsilon_0 / np.sqrt(permittivity))
  if (Field.imag > 1E-5): 
    print "** Error: unexpected imaginary part in the conversion from Intensity to Field!"
    exit(-1)
  return Field.real

## Builds Gaussian thickness from FWHM (in time or space)
def sigmaFWHM(FWHM):
  sigma = FWHM/(2.*np.sqrt(2.*np.log(2.)))
  return sigma

## Single pulse shape [table of peak_intensity(time)] evolution with time
#
# Defines the temporal shape of the laser pulse using a Gaussian law. 
# @param t: instant to output (can be a table)
# @param tau: pulse duration (s)
# @param PeakIntensity: peak intensity (W/m^2)
# @param t0: instant for the peak intensity (t0=0 by default)
def PulseGaussianTemporalShape(t, tau, PeakIntensity, t0=0.):
  sigmaTau = sigmaFWHM(tau)
  #PeakIntensity = fluence/tau 
  #TODO: Missing coefficient on peak intensity ? 
  intensity = PeakIntensity * np.exp(-0.5 * ((t-t0)/(sigmaTau))**2 )
  return intensity

## Single pulse shape [table of peak_intensity(time)] evolution with time
#
# Defines the temporal shape of the laser pulse using a squared sinus law. 
# Outputs: <Array of real-valued field envelope, array of complex electric field>
# @param t: instants to output (can be a table)
# @param tau: pulse duration FWHM (s)
# @param PeakField: peak of the electric field envelope (V/m) (scalar only)
# @param t0: central instant for the laser pulse (t0=0 by default)
# @param PulseDelay: temporal delay between 2 pulses (in seconds)
def PulseSquaredSinTemporalShape(t, tau, PeakField, wavelength, CEP=0., t0=0., PulseDelay=0.):
  t1 = t0 + PulseDelay
  omega = 2e0*pi*c/wavelength
  H1 = step(t - t1 + tau) #! theer could be a mistake in pulse duration here!
  H2 = step(t - t1 - tau)
  Envelope = PeakField*np.sin(pi*(t-t1-tau)/(2e0*tau))**2 * H1 * (1.-H2)
  Phase = np.exp(1e0j*omega*t+CEP)
  Field = Envelope * Phase
  return Envelope, Field

## Bi-color double pulse [table of TotalEnvelope(time), TotalField(time)] evolution with time (POLARIZATION IS FOR NOW NEGLECTED!)
# Output: Total envelope <array>, total field <array> at a given space point. 
# Construct the temporal shape of two-color laser pulses mixed together using a squared sinus law and a time delay. 
# Pulses CAN be of different wavelengths! 
# @param t: instants to output (can be a table)
# @param tau1: pulse 1 duration FWHM (s)
# @param tau2: pulse 2 duration FWHM (s)
# @param Efield1: pulse 1, (scalar) peak field of the envelope (V/m)
# @param Efield2: pulse 2, (scalar) peak field of the envelope (V/m)
# @param wavelength1: pulse 1, wavelength (meters)
# @param wavelength2: pulse 2, wavelength (meters)
# @param CEP1: pulse 1, carrier envelope phase
# @param CEP2: pulse 2, carrier envelope phase
# @param t1: instant for the peak field 1 (t1=0 by default)
# @param PulseDelay: delay between the amplitude maxima of pulse 1 and pulse 2 (seconds)
def PulseSquaredSinTemporalShapeDoublePulse(t, tau1, tau2, Efield1, Efield2, wavelength1, wavelength2, CEP1, CEP2, t1=0., PulseDelay=0.):
  omega1=2.*pi*c/wavelength1; omega2=2.*pi*c/wavelength2
  #sigmaTau1 = sigmaFWHM(tau1); sigmaTau2 = sigmaFWHM(tau2) #good for purely gaussian pulse, mmh? 
  t2 = t1 + PulseDelay
  H11        = step(t - t1 + tau1); H21 = step(t - t2 + tau2)
  H12        = step(t - t1 - tau1); H22 = step(t - t2 - tau2)
  FieldEnv1     = Efield1*np.sin(pi*(t-t1-tau1)/(2e0*tau1))**2 * H11 * (1.-H12) #could be bugged
  FieldEnv2     = Efield2*np.sin(pi*(t-t2-tau2)/(2e0*tau2))**2 * H21 * (1.-H22) #could be bugged
  Phase1 = np.exp(1e0j*omega1*t+CEP1)
  Phase2 = np.exp(1e0j*omega2*t+CEP2)
  #TotalEnvelope = np.sqrt( FieldEnv1*np.conj(FieldEnv1) + FieldEnv2*np.conj(FieldEnv2) + FieldEnv1*np.conj(FieldEnv2) * np.exp(1e0j*(omega1-omega2)*t) + np.conj(FieldEnv1)* FieldEnv2 * np.exp(1e0j*(omega2-omega1)*t) ) #complex square of the fields must provide the envelope
  TotalEnvelope = np.sqrt( FieldEnv1*np.conj(FieldEnv1) * np.exp(2.*CEP1) + FieldEnv2*np.conj(FieldEnv2) * np.exp(2.*CEP2) + FieldEnv1*np.conj(FieldEnv2) * np.exp(1e0j*(omega1-omega2)*t+CEP1+CEP2) + np.conj(FieldEnv1)* FieldEnv2 * np.exp(1e0j*(omega2-omega1)*t+CEP1+CEP2) ) #complex square of the fields must provide the envelope
  TotalField = FieldEnv1*Phase1 + FieldEnv2*Phase2
  return TotalEnvelope, TotalField

## Two-photon absorption probability from Bristow and Van Driel, 
# Applied Physics Letters 90, 191104 (2007)
def BristowLaw(wavelength, Egap):#{{{
  constant = 43E-11
  #Ep = 21. * e #plasmon peak for Si
  #Egap = 1.12 * e #gap of Si at 300 K
  
  ## Private function used in BristowLaw
  def F2(n, x): #{{{
    return ( pi*factorial2(2*n+1)/(2**(n+2)*factorial(n+2)) ) * (2.*x)**(-5) * (2.*x-1)**(n+2)
  #}}}
  beta = 2. * constant * F2(n, hbar*omega/Egap)
  return beta
#}}}

## Generate Keldysh tables for interfacing with codes
# Input: 
# @param Egap: scalar (J)
# @param meff: scalar (no unit)
# @param wavelength: laser wavelength (scalar, meters)
# @param PeakField: laser field amplitude (scalar, V/m)
# @param order: integration order for Keldysh model (integer, no unit)
def GenerateKeldyshDatabase(Egap, meff, wavelength, PeakField, order): #{{{
  ErrorMessage = ""
  gamma = gammaKeldysh(Egap, meff, PeakField, wavelength) #valid for scalar data
  k1 = Keldysh1(gamma); k2 = Keldysh2(gamma) #valid
  EgapEff = EffectiveGap(Egap, k1, k2) # Original formula from Keldysh. Warning: scipy.special.ellipe (https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.ellipe.html#scipy.special.ellipe) uses a different convention than Maple, Wikipedia or mpmath.
  KeldyshFunctionResultG = KeldyshFunction_Gruzdev( k1, k2, EgapEff, order, wavelength )
  wPIg = IonizationRate_Gruzdev(k1, k2, KeldyshFunctionResultG, EgapEff, wavelength)
  return wPIg
#}}}

#print "** Vectorizing functions..."
FieldToIntensity                        = np.vectorize(FieldToIntensity)
IntensityToField                        = np.vectorize(IntensityToField)
gammaKeldysh                            = np.vectorize(gammaKeldysh)
Keldysh1                                = np.vectorize(Keldysh1)
Keldysh2                                = np.vectorize(Keldysh2)
EffectiveGap                            = np.vectorize(EffectiveGap)
KeldyshFunction                         = np.vectorize(KeldyshFunction)
KeldyshFunction_Gruzdev                 = np.vectorize(KeldyshFunction_Gruzdev)
IonizationRate_Gruzdev                  = np.vectorize(IonizationRate_Gruzdev)
PulseGaussianTemporalShape              = np.vectorize(PulseGaussianTemporalShape)
PulseSquaredSinTemporalShape            = np.vectorize(PulseSquaredSinTemporalShape)
PulseSquaredSinTemporalShapeDoublePulse = np.vectorize(PulseSquaredSinTemporalShapeDoublePulse)
BristowLaw                              = np.vectorize(BristowLaw)
GenerateKeldyshDatabase                 = np.vectorize(GenerateKeldyshDatabase)

## Generate Keldysh tables for laser fields amplitudes, associated with a wavelength
# Input: 
# @param Egap: scalar (J)
# @param meff: scalar (no unit)
# @param wavelength: laser wavelength (scalar, meters)
# @param tau: total pulse duration (scalar, seconds)
# @param FieldEnvelope: envelope of the laser electric field (scalar|vector, in V/m)
# @param dt: precision (scalar, seconds)
# @param order: integration order for Keldysh model (integer, no unit)
# @N_total: limiter for the ionizable number of electrons (float, m^{-3})
def generateWpiTables(Egap = 2.56e0*e, meff = 0.2226e0, wavelength = 800e-9, tau = 10e-15, FieldEnvelope = 1e9, dt = 1E-17, order = 50, N_total=5E28, t0=0.0): #{{{
  Header="[libKeldysh] "
  ShortRefKeldysh = "[Keldysh (1964)]"
  ShortRefGruzdev = "[Gruzdev (2014)]"
  
  print Header+"Defining the laser pulse..."
  # t0 = 0e0
  tmin = -1.*tau+t0
  tmax =  1.*tau+t0
  #dt = 1E-17
  print Header+"** Info: Number of time steps = "+str(int((tmax-tmin)/dt))+"."
  instants=np.arange(tmin,tmax,dt)

  #Gaussian envelope
  #PulseEnvelope=PulseGaussianTemporalShape(instants, tau, Intensity, t0)
  #FieldEnvelope, PulseEnvelope = PulseSquaredSinTemporalShape(instants, tau, Intensity, wavelength, t0) #TODO: this should be generated outside this function
  IntensityEnvelope=FieldToIntensity(FieldEnvelope)
  print Header+"** Info: Peak intensity = "+str(np.max(IntensityEnvelope)/1E4)+" W/cm^2."
  print Header+"** Info: Peak field amplitude = "+str(np.max(FieldEnvelope)/1E9)+" V/nm."

  #print "** Starting the self-consistent loop..."

  #for i in np.arange(1,50,1): #attempt of self consistent loop: divergent
  #print "** ITERATION "+str(i)
  print Header+"Computing Adiabadicity coefficients for the pulse envelope..."
  
  gamma = gammaKeldysh(Egap, meff, FieldEnvelope, wavelength) #valid for scalar|vector data
  #gamma = gammaKeldysh(EgapEff, meff, IntensityToField(PulseEnvelope), wavelength) #self-consistent, divergent
  print Header+"** Info: Adiabadicity parameter = "+str(gamma.min())+"."

  #print "Computing Keldysh1, Keldysh2..."
  k1 = Keldysh1(gamma); k2 = Keldysh2(gamma) #valid
  EgapEff = EffectiveGap(Egap, k1, k2) # Original formula from Keldysh. Warning: scipy.special.ellipe (https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.ellipe.html#scipy.special.ellipe) uses a different convention than Maple, Wikipedia or mpmath.
  #EgapEff = EffectiveGap(EgapEff, k1, k2) #This formula was made self-consistent, but divergent. 
  #order = 50

  print ""
  print Header+"** Info: Egap = "+str(Egap/e)+" eV, max[Ueff] = "+str(EgapEff.max()/e)+" eV."
  print ""

  KeldyshFunctionResult  = KeldyshFunction( k1, k2, EgapEff, order, wavelength )
  KeldyshFunctionResultG = KeldyshFunction_Gruzdev( k1, k2, EgapEff, order, wavelength )

  #print "** End of self-consistent loop..."
  print Header+"Computes w_PI (Keldysh), and w_PIg (Gruzdev) for the pulse envelope..."
  wPI = IonizationRate(k1, k2, KeldyshFunctionResult, EgapEff, wavelength)
  wPIg = IonizationRate_Gruzdev(k1, k2, KeldyshFunctionResultG, EgapEff, wavelength)
  print Header+"w_PI until order "+str(order)+" = ", wPI.max()

  #print "Developing: exporting the table..."

  print ""
  print Header+"Temporal integration..."
  
  # Temporal integration without limiter
  
  # Just multiply array of w_PI by dt, with limited to Ntotal
  #N_excited_Keldysh = dN_excited_Keldysh.cumsum()
  #N_excited_Gruzdev = dN_excited_Gruzdev.cumsum()
  
  # Temporal integration with limiter
  dN_excited_Keldysh = np.multiply(wPI, dt)
  dN_excited_Gruzdev = np.multiply(wPIg, dt)
  #dN_excited_Bristow = np.multiply(BristowLaw(wavelength, Egap)*intensity**2/(2*hbar*omega)), dt)
  
  # Initial number of electrons in conduction band
  N_initial = np.zeros(wPI.shape)
  
  ExpArg_Keldysh = np.divide(dN_excited_Keldysh.cumsum(), N_total)
  ExpArg_Gruzdev = np.divide(dN_excited_Gruzdev.cumsum(), N_total)
  
  N_excited_Keldysh = np.multiply( np.exp(-ExpArg_Keldysh), N_total * np.exp(ExpArg_Keldysh) - N_total + N_initial)
  
  N_excited_Gruzdev = np.multiply( np.exp(-ExpArg_Gruzdev), N_total * np.exp(ExpArg_Gruzdev) - N_total + N_initial)
  return instants, N_excited_Keldysh, N_excited_Gruzdev, gamma, wPI, wPIg
#}}}

#generateWpiTables = np.vectorize(generateWpiTables)

## Compute and plot density evolution with time using the specific parameters.
# @param Egap (Joules): direct band gap of the modeled material
# @param meff (adim): effective mass of the conduction band
# @param wavelength (meters): photon wavelength of the excitation
# @param tau (seconds): duration of the Gaussian pulse (FWHM)
# @param FieldEnvelope (V/m): electric field envelope of the pulse (scalar|vector)
# @param dt (seconds): precision of the temporal envelope
# @param order (adim): order of the integration (default: 50).
def plotPulseToDensity(Egap = 2.56e0*e, meff = 0.2226e0, wavelength = 800e-9, tau = 10e-15, FieldEnvelope = 1e9, dt = 1E-17, order = 50, ShowPlot=False, t0=0., N_total=5E28): #{{{
  Header="[libKeldysh] "
  ShortRefKeldysh = "[Keldysh (1964)]"
  ShortRefGruzdev = "[Gruzdev (2014)]"
  #gamma = gammaKeldysh(Egap, meff, FieldEnvelope, wavelength)
  
  ### We shall generate the interesting pulse in the file from which we call the Keldysh generator
  instants, N_excited_Keldysh, N_excited_Gruzdev, gamma, wPI, wPIg = generateWpiTables(Egap, meff, wavelength, tau, FieldEnvelope, dt, order, N_total, t0)
  
  sizeGamma = len(gamma)
  
  print ""
  print Header+"** Warning: results may be not converged."
  print Header+"            Reduce dt, and increase order until convergence."
  print ""
  print Header+"Maximum density N_ex "+ShortRefKeldysh+" = "+str(N_excited_Keldysh.max())+"."
  print Header+"Maximum density N_ex "+ShortRefGruzdev+" = "+str(N_excited_Gruzdev.max())+"."
  print ""
#  return N_excited_Gruzdev 
  if(ShowPlot): 
    print Header+"Plotting..."

    xunit = 1E15
    timeunit = "fs"

    #plt.figure()
    plt.figure(figsize=(15,15))
    
    plt.subplot(411)
    plt.title(r"Gap = "+str(Egap/e)+" eV, $\lambda=$ "+str(wavelength*1E9)+r" nm, $\tau=$"+str(tau*xunit)+" "+timeunit+", "+r"$E_{max}=$"+str(FieldEnvelope.max()/1E9)+" V/nm")
    #plt.xlabel("Field (V/m)")
    #plt.xlabel("Time (ps)")
    plt.ylabel("Field envelope (V/m)")
    plt.plot(xunit*instants[0:sizeGamma], FieldEnvelope, color="k", linestyle="-", label=r"Field envelope")
    plt.grid()
    # plt.loglog(Efield, 0.1, label="Tunnelling limit")
    # plt.loglog(Efield, 10.*np.ones(), label="MPI limit")
    plt.legend(loc=3)
    
    plt.subplot(412)
    #plt.xlabel("Field (V/m)")
    #plt.xlabel("Time (ps)")
    plt.ylabel("Adiabadicity $\gamma$")
    plt.semilogy(xunit*instants[0:sizeGamma], gamma, color="k", linestyle="-", label=r"$\gamma$")
    plt.grid()
    # plt.loglog(Efield, 0.1, label="Tunnelling limit")
    # plt.loglog(Efield, 10.*np.ones(), label="MPI limit")
    plt.legend(loc=3)

    plt.subplot(413)
    #plt.xlabel("Field (V/m)")
    #plt.xlabel("Time (ps)")
    plt.ylabel("$w_{PI}$ (m$^{-3}$ s$^{-1}$)")
    plt.plot(instants[0:sizeGamma]*xunit, wPI, linestyle="-", color="r", label=r"$w_{PI}$ "+ShortRefKeldysh)
    plt.plot(instants[0:sizeGamma]*xunit, wPIg, linestyle="-", color="b", label=r"$w_{PI}$ "+ShortRefGruzdev)
    plt.grid()
    plt.legend(loc=2)

    plt.subplot(414)
    #plt.xlabel("Intensity (W/m$^{2}$)")
    plt.xlabel("Time ("+timeunit+")")
    # plt.xlabel("Field (V/m)")
    plt.ylabel("Density (m$^{-3}$)")
    plt.plot(instants[0:sizeGamma]*xunit, N_excited_Keldysh, color="r", label="$n_e$ "+ShortRefKeldysh)
    plt.plot(instants[0:sizeGamma]*xunit, N_excited_Gruzdev, color="b", label="$n_e$ "+ShortRefGruzdev)
    plt.grid()
    plt.legend(loc=2)
    plt.savefig("KeldyshSimple.eps") 
  
  return instants, N_excited_Keldysh, N_excited_Gruzdev
#}}}

##############################

## Converts field CGS units (statV/cm) in SI (V/m).
# @param CGS: input field in CGS units
# Retuns the field in SI units (V/m).
def Field_CGS_to_SI(CGS):
  c_CGS      = c*1e2
  conversion = 1E-6*c_CGS*1E2
  return CGS/conversion

def Field_SI_to_CGS(SI):
  c_CGS      = c*1e2
  conversion = 1E-6*c_CGS*1E2
  return SI*conversion

## Converts a mass in g (CGS unit) to kg (SI)
# @param CGS: mass in g (CGS unit)
def Mass_CGS_to_SI(CGS):
  return CGS * 1E-3

## Converts a mass from kg (SI) to g (CGS unit)
# @param SI: mass in kg (SI unit)
def Mass_SI_to_CGS(SI):
  return SI * 1E3

## Converts an energy in ergs (CGS unit) to Joules (SI unit)
# @param CGS: energy in ergs (CGS unit)
def Energy_CGS_to_SI(CGS):
  return CGS * 1E-7

## Converts an energy from Joules (SI unit) to ergs (CGS unit) 
# @param SI: energy in Joules (SI unit)
def Energy_SI_to_CGS(SI):
  return SI / 1E-7

## Converts velocity from m/s (SI unit) to cm/s (CGS unit)
# @param CGS: velocity in cm/s
def Velocity_CGS_to_SI(CGS):
  return CGS * 1E2

## Converts velocity from cm/s (CGS unit) to m/s (SI unit)
# @param SI: velocity in m/s
def Velocity_SI_to_CGS(SI):
  return SI / 1E-2

## Converts electric charge in statC (CGS unit) to Coulomb (SI unit)
def electric_charge_CGS_to_SI(CGS):
  c_CGS = Velocity_SI_to_CGS(c)
  return CGS / c_CGS #TODO: verify if consistent!

## Converts electric charge in Coulomb (SI unit) to statC (CGS unit)
def electric_charge_SI_to_CGS(SI):
  c_CGS = Velocity_SI_to_CGS(c)
  return SI * c_CGS

def Length_SI_to_CGS(SI):
  return SI * 1e2

def Length_CGS_to_SI(CGS):
  return CGS / 1e2

## Generates the normalization coefficients for electric field
# /!\ Vladimir uses adiabadicity coefficient for gas, which differs from a 1/sqrt(2) factor. 
# /!\ Vladimir also uses CGS units. 
# See the paper [Keldysh, L. Ionization in the field of a strong electromagnetic wave Journal of Experimental and Theoretical Physics, Lebedev Inst. of Physics, Moscow, 1964, 47, 5]
# @param Egap: band gap value (Joules, SI)
# @param meff: effective mass (no unit)
# @param wavelength: wavelength of the photons (in meters, SI)
# @param NormalizedPeakField: normalized field to be obtained in CGS
def VZ_FieldNormalization(Egap, meff, wavelength):
  # TODO: the function remains to be validated! 
  # field for which gamma_VZ = 1. 
  me_CGS = Mass_SI_to_CGS(m_e * meff) #[1 kg    (SI) = 1E3  g      (CGS) ]
  Eg_CGS = Energy_SI_to_CGS(Egap)     #[1 J     (SI) = 1E7  ergs   (CGS) ]
  c_CGS  = Velocity_SI_to_CGS(c)      #[1 [m/s] (SI) = 1E2 cm/s   (CGS) ]
  e_CGS  = electric_charge_SI_to_CGS(e) / c_CGS  #[1 C     (SI) = c_CGS \times statC (CGS) ]
  wavelength_CGS = Length_SI_to_CGS(wavelength)
  omega_CGS = 2.*pi*c_CGS/wavelength_CGS
  EfieldStar_VZ_CGS = np.sqrt(2.*me_CGS*omega_CGS**2 * Eg_CGS / e_CGS**2) #gas formula for Keldysh parameter
  EfieldStar_VZ_SI  = Field_CGS_to_SI(EfieldStar_VZ_CGS)
  return EfieldStar_VZ_SI, EfieldStar_VZ_CGS

VZ_FieldNormalization = np.vectorize(VZ_FieldNormalization)

## Interfaces bicolor tables of V. Zhukov bicolor Keldysh model with the laser parameters
# Returns the W_PI coefficients from V. Zhukov model to be integrated in time for bi-color laser pulses
def VZ_generateWpiTables(FieldEnvelope1, FieldEnvelope2, wavelength1 = 800e-9, wavelength2 = 800e-9, Egap=2.56*e, meff=0.2226, tau1=10e-15, tau2=10e-15, Delay=0., dt = 1E-17, N_total=5E28, t0=0.0): #{{{
  Header="[libKeldysh] "
  
  print Header+"Defining the laser pulse..."
  # t0 = 0e0
  tmin = -1.*tau1+t0
  tmax =  1.*tau2+t0 + Delay
  
  print Header+"** Info: Number of time steps = "+str(int((tmax-tmin)/dt))+"."
  instants=np.arange(tmin,tmax,dt)
  
  # Printing info on the pulses
  IntensityEnvelope1=FieldToIntensity(FieldEnvelope1)
  print Header+"** Info: Peak intensity 1= "+str(np.max(IntensityEnvelope1)/1E4)+" W/cm^2."
  print Header+"** Info: Peak field amplitude 1= "+str(np.max(FieldEnvelope1)/1E9)+" V/nm."
  
  IntensityEnvelope2=FieldToIntensity(FieldEnvelope2)
  print Header+"** Info: Peak intensity 2= "+str(np.max(IntensityEnvelope2)/1E4)+" W/cm^2."
  print Header+"** Info: Peak field amplitude 2= "+str(np.max(FieldEnvelope2)/1E9)+" V/nm."

  # 0. Expressing FieldEnvelopes in CGS
  FieldNormalizationCoeff1_SI, FieldNormalizationCoeff1_CGS = VZ_FieldNormalization(Egap, meff, wavelength1)
  FieldNormalizationCoeff2_SI, FieldNormalizationCoeff2_CGS = VZ_FieldNormalization(Egap, meff, wavelength2)
  
  print Header+"** Normalization coefficient Field 1 [SI]: "+str(FieldNormalizationCoeff1_SI)
  print Header+"** Normalization coefficient Field 1 [CGS]: "+str(FieldNormalizationCoeff1_CGS)
  print Header+"** Normalization coefficient Field 2 [SI]: "+str(FieldNormalizationCoeff2_SI)
  print Header+"** Normalization coefficient Field 2 [CGS]: "+str(FieldNormalizationCoeff2_CGS)
  
  # Normalize FieldEnvelope 1,2 in CGS. Vladimir requires normalized field1 and normalized field2 to deliver a W_PI. Note: his formula for gamma is the one for gas and does not account for optical Stark effect (increase of gap with field strength). #TODO: Why ? Stark effect also happens in gas. 
  
  FieldEnvelope1_CGS = Field_SI_to_CGS(FieldEnvelope1)
  FieldEnvelope2_CGS = Field_SI_to_CGS(FieldEnvelope2)
  
  print Header+"Field1 [CGS] = "+str(FieldEnvelope1_CGS.max())
  print Header+"Field2 [CGS] = "+str(FieldEnvelope2_CGS.max())
  
  FieldEnvelopeNormalized1_CGS = FieldEnvelope1_CGS / FieldNormalizationCoeff1_CGS
  FieldEnvelopeNormalized2_CGS = FieldEnvelope2_CGS / FieldNormalizationCoeff2_CGS
  
  FieldEnvelopeNormalized1_SI = FieldEnvelope1 / FieldNormalizationCoeff1_SI
  FieldEnvelopeNormalized2_SI = FieldEnvelope2 / FieldNormalizationCoeff2_SI
  
  print Header+"Normalized Field1 is now [CGS]: "+str(FieldEnvelopeNormalized1_CGS.max())
  print Header+"Normalized Field2 is now [CGS]: "+str(FieldEnvelopeNormalized2_CGS.max())
  print Header+"Normalized Field1 is now [SI]: "+str(FieldEnvelopeNormalized1_SI.max())
  print Header+"Normalized Field2 is now [SI]: "+str(FieldEnvelopeNormalized2_SI.max())
  print ""
  print Header+"** Info: Egap = "+str(Egap/e)+" eV"
  print ""

  #print "** End of self-consistent loop..."
  print Header+"Computes w_PI (Keldysh) for the bicolor pulse envelope..."
  
  print Header+"max(w_PI) = ", w_PI.max()

  #print "Developing: exporting the table..."

  print ""
  print Header+"Temporal integration..."
  
  # Temporal integration without limiter
  
  # Just multiply array of w_PI by dt, with limited to Ntotal
  #N_excited_Keldysh = dN_excited_Keldysh.cumsum()
  #N_excited_Gruzdev = dN_excited_Gruzdev.cumsum()
  
  # Temporal integration with limiter
  dN_excited_Zhukov = np.multiply(wPI, dt)
  
  # Initial number of electrons in conduction band
  N_initial = np.zeros(wPI.shape)
  
  ExpArg_Zhukov = np.divide(dN_excited_Zhukov.cumsum(), N_total)
  
  N_excited_Zhukov = np.multiply( np.exp(-ExpArg_Zhukov), N_total * np.exp(ExpArg_Zhukov) - N_total + N_initial)
  
  #N_excited_Gruzdev = np.multiply( np.exp(-ExpArg_Gruzdev), N_total * np.exp(ExpArg_Gruzdev) - N_total + N_initial)
  return instants, N_excited_Zhukov, wPI, wPIg
#}}}
  
  
  
  
  