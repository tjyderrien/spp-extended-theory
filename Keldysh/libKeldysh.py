#!/usr/bin/env python
#-*- coding: utf-8 -*-
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
rc('font', **{'family':'serif', 'serif':['Palatino'], 'size':'18'})
rc('text', usetex=True)
mp.rcParams['legend.numpoints'] = 1

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
  if (abs(Efield) > 1e0):
    result = omegaLaser*np.sqrt(m_e*meff*Egap)/e/Efield
  else:
    ErrorMessage=ErrorMessage+"gamma(): Divergence, as field equals = 0. Singular case of Keldysh functions. Should give w_PI = 0 then...\n"
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
def FieldToIntensity(Field, permittivity=1):
  intensity = 0.5 * c * epsilon_0 * np.sqrt(permittivity) * Field**2
  return intensity.real

## Converts intensity to electric field amplitude
def IntensityToField(intensity, permittivity=1):
  Field = np.sqrt(2.0 * intensity / c / epsilon_0 / np.sqrt(permittivity))
  if (Field.imag != 0): 
    print "** Error: unexpected imaginary part in the conversion from Intensity to Field!"
    exit(-1)
  return Field.real

## Builds Gaussian thickness from FWHM (in time or space)
def sigmaFWHM(FWHM):
  sigma = FWHM/(2.*np.sqrt(2.*np.log(2.)))
  return sigma

## Pulse shape [table of peak_intensity(time)] evolution with time
#
# Defines the temporal shape of the laser pulse using a Gaussian law. 
# @param t: instant to output (can be a table)
# @param tau: pulse duration (s)
# @param intensity: peak intensity (W/m^2)
# @param t0: instant for the peak intensity (t0=0 by default)
def PulseGaussianTemporalShape(t, tau, PeakIntensity, t0=0.):
  sigmaTau = sigmaFWHM(tau)
  #PeakIntensity = fluence/tau 
  #TODO: Missing coefficient on peak intensity ? 
  intensity = PeakIntensity * np.exp(-0.5 * ((t-t0)/(sigmaTau))**2 )
  return intensity

## Calculate the time-dependent excited electron density for a given pulse shape
# @param IntensityShape: function describing the temporal enveloppe of the pulse
# @param 
#def TimeDependentDensityKeldysh(IntensityShape, tau, intensity):
  #%IntensityShape(

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
FieldToIntensity           = np.vectorize(FieldToIntensity)
IntensityToField           = np.vectorize(IntensityToField)
gammaKeldysh               = np.vectorize(gammaKeldysh)
Keldysh1                   = np.vectorize(Keldysh1)
Keldysh2                   = np.vectorize(Keldysh2)
EffectiveGap               = np.vectorize(EffectiveGap)
KeldyshFunction            = np.vectorize(KeldyshFunction)
KeldyshFunction_Gruzdev    = np.vectorize(KeldyshFunction_Gruzdev)
IonizationRate_Gruzdev     = np.vectorize(IonizationRate_Gruzdev)
PulseGaussianTemporalShape = np.vectorize(PulseGaussianTemporalShape)
BristowLaw                 = np.vectorize(BristowLaw)
GenerateKeldyshDatabase    = np.vectorize(GenerateKeldyshDatabase)

## Generate Keldysh tables for a range of laser intensity, associated with a wavelength
# Input: 
# @param Egap: scalar (J)
# @param meff: scalar (no unit)
# @param wavelength: laser wavelength (scalar, meters)
# @param tau: pulse duration (scalar, seconds)
# @param PeakFluence: peak fluence of the laser (scalar, J/m2)
# @param dt: precision (scalar, seconds)
# @param order: integration order for Keldysh model (integer, no unit)
# @N_total: limiter for the ionizable number of electrons (float, m^{-3})
def generateWpiTables(Egap = 2.58e0*e, meff = 0.18e0, wavelength = 800e-9, tau = 10e-15, PeakFluence = 100e-3*1E4, dt = 1E-17, order = 50, N_total=5E28, t0=0.0): #{{{
  ShortRefKeldysh = "[Keldysh (1964)]"
  ShortRefGruzdev = "[Gruzdev (2014)]"
  
  print "Defining the laser pulse..."
  PeakIntensity = PeakFluence/tau #scalar, TODO: isn't it multiplied by sqrt(4 ln 2 / Pi) ? 
  # t0 = 0e0
  tmin = -3.5*tau+t0
  tmax = 3.5*tau+t0
  #dt = 1E-17
  print "** Info: Number of time steps = "+str(int((tmax-tmin)/dt))+"."
  instants=np.arange(tmin,tmax,dt)

  PulseEnvelope=PulseGaussianTemporalShape(instants, tau, PeakIntensity, t0)
  print "** Info: Peak intensity = "+str(PulseEnvelope.max()/1E4)+" W/cm^2."
  print "** Info: Peak field amplitude = "+str(IntensityToField(PulseEnvelope).max()/1E9)+" V/nm."

  #print "** Starting the self-consistent loop..."

  #for i in np.arange(1,50,1): #attempt of self consistent loop: divergent
  #print "** ITERATION "+str(i)
  print "Computing Adiabadicity coefficients for the pulse envelope..."

  gamma = gammaKeldysh(Egap, meff, IntensityToField(PulseEnvelope), wavelength) #valid for scalar data
  #gamma = gammaKeldysh(EgapEff, meff, IntensityToField(PulseEnvelope), wavelength) #self-consistent, divergent
  print "** Info: Adiabadicity parameter = "+str(gamma.min())+"."

  #print "Computing Keldysh1, Keldysh2..."
  k1 = Keldysh1(gamma); k2 = Keldysh2(gamma) #valid
  EgapEff = EffectiveGap(Egap, k1, k2) # Original formula from Keldysh. Warning: scipy.special.ellipe (https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.ellipe.html#scipy.special.ellipe) uses a different convention than Maple, Wikipedia or mpmath.
  #EgapEff = EffectiveGap(EgapEff, k1, k2) #This formula was made self-consistent, but divergent. 
  #order = 50

  print ""
  print "** Info: Egap = "+str(Egap/e)+" eV, max[Ueff] = "+str(EgapEff.max()/e)+" eV."
  print ""

  KeldyshFunctionResult = KeldyshFunction( k1, k2, EgapEff, order, wavelength )
  KeldyshFunctionResultG = KeldyshFunction_Gruzdev( k1, k2, EgapEff, order, wavelength )

  #print "** End of self-consistent loop..."
  print "Computes w_PI (Keldysh), and w_PIg (Gruzdev) for the pulse envelope..."
  wPI = IonizationRate(k1, k2, KeldyshFunctionResult, EgapEff, wavelength)
  wPIg = IonizationRate_Gruzdev(k1, k2, KeldyshFunctionResultG, EgapEff, wavelength)
  print "w_PI until order "+str(order)+" = ", wPI.max()

  #print "Developing: exporting the table..."

  print ""
  print "Temporal integration..."
  
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
  return instants, N_excited_Keldysh, N_excited_Gruzdev
#}}}

#generateWpiTables = np.vectorize(generateWpiTables)

## Compute and plot density evolution with time using the specific parameters.
# @param Egap (Joules): direct band gap of the modeled material
# @param meff (adim): effective mass of the conduction band
# @param wavelength (meters): photon wavelength of the excitation
# @param tau (seconds): duration of the Gaussian pulse (FWHM)
# @param PeakFluence (J/m2): maximum fluence of the pulse
# @param dt (seconds): precision of the temporal envelope
# @param order (adim): order of the integration (default: 50).
def plotPulseToDensity(Egap = 2.56e0*e, meff = 0.18e0, wavelength = 800e-9, tau = 10e-15, PeakFluence = 100e-3*1E4, dt = 1E-17, order = 50, ShowPlot=False, t0=0., N_total=5E28): #{{{

  ShortRefKeldysh = "[Keldysh (1964)]"
  ShortRefGruzdev = "[Gruzdev (2014)]"
  
  instants, N_excited_Keldysh, N_excited_Gruzdev = generateWpiTables(Egap, meff, wavelength, tau, PeakFluence, dt, order, N_total, t0)
  
  print ""
  print "** Warning: results may be not converged."
  print "            Reduce dt, and increase order until convergence."
  print ""
  print "Maximum density N_ex "+ShortRefKeldysh+" = "+str(N_excited_Keldysh.max())+"."
  print "Maximum density N_ex "+ShortRefGruzdev+" = "+str(N_excited_Gruzdev.max())+"."
  print ""
#  return N_excited_Gruzdev 
  if(ShowPlot): 
    print "Plotting..."

    xunit = 1E15
    timeunit = "fs"

    #plt.figure()
    plt.figure(figsize=(15,15))
    plt.subplot(311)
    plt.title(r"Gap = "+str(Egap/e)+" eV, $\lambda=$ "+str(wavelength*1E9)+r" nm, $\tau=$"+str(tau*xunit)+" "+timeunit+", "+r"$F_{max}=$"+str(PeakFluence/1E4)+" J/cm"+r"$^{2}$")
    #plt.xlabel("Field (V/m)")
    #plt.xlabel("Time (ps)")
    plt.ylabel("Adiabadicity $\gamma$")
    plt.semilogy(instants*xunit, gamma, color="k", linestyle="-", label=r"$\gamma$")
    plt.grid()
    # plt.loglog(Efield, 0.1, label="Tunnelling limit")
    # plt.loglog(Efield, 10.*np.ones(), label="MPI limit")
    plt.legend(loc=3)

    plt.subplot(312)
    #plt.xlabel("Field (V/m)")
    #plt.xlabel("Time (ps)")
    plt.ylabel("$w_{PI}$ (m$^{-3}$ s$^{-1}$)")
    plt.plot(instants*xunit, wPI, linestyle="-", color="r", label=r"$w_{PI}$ "+ShortRefKeldysh)
    plt.plot(instants*xunit, wPIg, linestyle="-", color="b", label=r"$w_{PI}$ "+ShortRefGruzdev)
    plt.grid()
    plt.legend(loc=2)

    plt.subplot(313)
    plt.xlabel("Intensity (W/m$^{2}$)")
    plt.xlabel("Time ("+timeunit+")")
    # plt.xlabel("Field (V/m)")
    plt.ylabel("Density (m$^{-3}$)")
    plt.plot(instants*xunit, N_excited_Keldysh, color="r", label="$n_e$ "+ShortRefKeldysh)
    plt.plot(instants*xunit, N_excited_Gruzdev, color="b", label="$n_e$ "+ShortRefGruzdev)
    plt.grid()
    plt.legend(loc=2)
    plt.savefig("KeldyshAnalytic.eps") 
  
  return instants, N_excited_Keldysh, N_excited_Gruzdev
#}}}
