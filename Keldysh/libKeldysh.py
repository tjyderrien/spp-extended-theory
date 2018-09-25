#!/usr/bin/env python2
#-*- coding: utf-8 -*-

# Copyright (C) 2013-2018 T. J.-Y. Derrien
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
# Several flavors of the Keldysh theory are available: 
# - Keldysh original paper in solid, for Kane band structure [compared with td-dft]
# - Keldysh paper with few terms corrected by Gruzdev [compared with td-dft]
# - Keldysh-Zhukov tables, where Keldysh theory was computed numerically without using the saddle point method
# - Keldysh-Shcheblanov model, improving rigor on the analytical integration [https://arxiv.org/abs/1706.07303]
# - Keldysh-Corkum model, allowing for analytical treatment of mulltiwavelength fields [Physical Review Letters, 2017, 118, 173601]

from libKeldyshZhukov import *
#from libKeldyshUlrich import * #NOT READY YET
from libKeldyshPulses import *

rc('font', **{'family':'serif', 'serif':['Helvetica'], 'size':'16'})
rc('text', usetex=False)
mp.rcParams['legend.numpoints'] = 1

Header="[libKeldysh] "

ShortRefKeldysh = "[Keldysh (1964)]"
ShortRefGruzdev = "[Gruzdev (2014)]"
ShortRefGulley  = "[Gulley (2012)]"

## Computes the adiabadicity parameter
# @param gamma: Adiabadicity parameter (non-dimensional number)
# @param Egap: band gap energy (in Joules)
# @param meff: effective mass (in Arb. Units, as it is multiplied by electron mass INSIDE the function)
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
    result = 1E99 #THIS VALUE IS ARBITRARY FOR A VERY SMALL FIELD. 
  #print(ErrorMessage)
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
  print "** Debug info: Keldysh1 = "+str(k1)+", Keldysh2 = "+str(k2)
  
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

## Formula for tunnel ionization in semiconductors [Keldysh 1964, Eq. (60)]
def KeldyshTunnelingLimit(Egap, meff, wavelength, Efield): #{{{
    
    gamma = gammaKeldysh(Egap, meff, Efield, wavelength) #valid for scalar data
    k1 = Keldysh1(gamma); k2 = Keldysh2(gamma) #valid
    EgapEff = EffectiveGap(Egap, k1, k2)
    
    omegaLaser=2.*pi*c/wavelength
    
    w_tunnel = 2./9./np.pi**2 * EgapEff / hbar * (m_e*meff*EgapEff/hbar**2)**1.5 * (e*hbar*Efield/(m_e*meff)**0.5/EgapEff**1.5)**2.5*np.exp(-0.5*np.pi*(m_e*meff)**0.5*EgapEff**1.5/e/hbar/Efield * (1.-1./8.*(m_e*meff)*omegaLaser**2*EgapEff/e**2/Efield**2))
    return w_tunnel
#}}}
    
## Compute and plot density evolution with time using the specific parameters.
# @param Egap (Joules): direct band gap of the modeled material
# @param meff (adim): effective mass of the conduction band
# @param wavelength (meters): photon wavelength of the excitation
# @param tau (seconds): duration of the Gaussian pulse (FWHM)
# @param FieldEnvelope (V/m): electric field envelope of the pulse (scalar|vector)
# @param dt (seconds): precision of the temporal envelope
# @param order (adim): order of the integration (default: 50).
def plotPulseToDensity(Egap = 2.56e0*e, meff = 0.2226e0, wavelength = 800e-9, tau = 10e-15, FieldEnvelope = 1e9, dt = 1E-17, order = 50, ShowPlot=False, t0=0., N_total=5E28): #{{{
  #Header="[libKeldysh] "
  #gamma = gammaKeldysh(Egap, meff, FieldEnvelope, wavelength)
  
  ### We shall generate the interesting pulse in the file from which we call the Keldysh generator
  instants, N_excited_Keldysh, N_excited_Gruzdev, gamma, wPI, wPIg = generateWpiTables(Egap, meff, wavelength, tau, FieldEnvelope, dt, order, N_total, t0)
  
  try: 
    sizeGamma = len(gamma)
  except: 
    print Header+"Error: problem on gamma in plotPulseToDensity():"+str(gamma)
  
  print ""
  print Header+"** Warning: results may be not converged."
  print Header+"            Reduce dt, and increase order until convergence."
  print ""
  print Header+"Maximum density N_ex "+ShortRefKeldysh+" = "+str(N_excited_Keldysh.max())+"."
  print Header+"Maximum density N_ex "+ShortRefGruzdev+" = "+str(N_excited_Gruzdev.max())+"."
  print ""
#  return N_excited_Gruzdev 
  if(ShowPlot): 
    print Header+"Plotting as function of time..."

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
    plt.legend(loc=3)
    
    plt.subplot(412)
    #plt.xlabel("Field (V/m)")
    #plt.xlabel("Time (ps)")
    plt.ylabel("Adiabadicity $\gamma$")
    plt.semilogy(xunit*instants[0:sizeGamma], gamma, color="k", linestyle="-", label=r"$\gamma$")
    plt.grid()
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
    
    plt.xlabel("Time ("+timeunit+")")
    # plt.xlabel("Field (V/m)")
    plt.ylabel("Density (m$^{-3}$)")
    plt.plot(instants[0:sizeGamma]*xunit, N_excited_Keldysh, color="r", label="$n_e$ "+ShortRefKeldysh)
    plt.plot(instants[0:sizeGamma]*xunit, N_excited_Gruzdev, color="b", label="$n_e$ "+ShortRefGruzdev)
    plt.grid()
    plt.legend(loc=2)
    plt.savefig("KeldyshSimple.eps") 
    
    ### Second plot
    print Header+"Importing Gulley [2012] data..."
    try:
      Gulley2012=np.loadtxt("Gulley-Fig2.csv", dtype='float', delimiter='\t')
      #print Gulley2012[:,0]
    except: 
      print Header+"** Warning: failed to import Gulley2012 data table..."
      
    
    print Header+"Plotting as function of laser field envelope..."
    plt.figure()
    plt.xlabel("Intensity (W/m$^{2}$)")
    plt.ylabel("$w_{PI}$ (m$^{-3}$ s$^{-1}$)")
    #plt.loglog(FieldToIntensity(FieldEnvelope.real), wPI,  linestyle="-", color="r", label=r"$w_{PI}$ "+ShortRefKeldysh)
    plt.loglog(FieldToIntensity(FieldEnvelope.real), wPIg, linestyle="-", color="b", label=r"$w_{PI}$ "+ShortRefGruzdev)
    #plt.loglog(Gulley2012[:,0], Gulley2012[:,1], linestyle="-", color="k", label="Data from "+ShortRefGulley) #JUST FOR VALIDATION. 
    plt.grid()
    plt.legend(loc=2)
    plt.xlim((1E11*1E4, 1E13*1E4))
    plt.ylim((1E20*1E6,1E40*1E6))
    plt.savefig("Keldysh-Field-Wpi.eps") 
    
  return instants, N_excited_Keldysh, N_excited_Gruzdev
#}}}


## Extends the Keldysh-Gruzdev models to parameters required by Stephane Gräf to analyze nanostructure formation in SiO2. 
def SilicaGraef2017(PeakFluence): #{{{
  print "Defining SiO2 material parameters from [Gräf2017]..."

  Egap = 8.024234328*e; #band gap of SiO2
  meff = 1e0; #Effective mass of SiO2
  Ntotal=10.*5E28; #valence band electron density #to avoid limitation

  wavelength = 1025e-9; #wavelength2 = 800e-9
  tau=300e-15; dt = 1E-17; CEP=0e0
  #PeakFluence = 5E4
  PeakField   = np.sqrt(2e0 * PeakFluence / (tau * c * epsilon_0))

  t0=0. #defines the instant 0.
  Delay = 0. #delay between maxima of the pulses
  tmin=-1.*tau + t0; tmax=1.*tau + Delay + t0

  instants = np.arange(tmin, tmax, dt)
  #print "Time range: "+str(instants.min())+", "+str(instants.max())+"."

  PeakField2  = 0. #/ 2.
  CEP2        = 0. #pi/3.
  #wavelength2 = wavelength
  
  print Header+"** Test: building single pulse centered on 0..."
  FieldEnvelope1, RealField1 = PulseSquaredSinTemporalShape(instants, tau, PeakField, wavelength, CEP, t0, 0.)

  #print Header+"** Test: We build a second pulse with a delay..."
  #FieldEnvelope2, RealField2 = PulseSquaredSinTemporalShape(instants, tau, PeakField2, wavelength2, CEP, t0, Delay)

  #print Header+"** Test: building a bicolor double pulse"

  #FieldEnvelopeTot, RealFieldTot = PulseSquaredSinTemporalShapeDoublePulse(instants, tau, tau, PeakField, PeakField, wavelength, wavelength2, CEP, CEP2, t0, Delay)

  plt.plot(instants, RealField1.real, '-')
  plt.plot(instants, FieldEnvelope1.real, '--')
  #plt.plot(instants, RealField2.real, '-')
  #plt.plot(instants, FieldEnvelope2.real, '--')
  #plt.plot(instants, RealFieldTot.real, '-')
  #plt.plot(instants, FieldEnvelopeTot.real, '--')
  plt.xlabel('')
  plt.savefig('PulseEnvelopes.eps')
  plt.savefig('PulseEnvelopes.png')
  #plt.show()

  print Header+"** Info: PulseEnvelope.EPS and PNG were written in the current folder. "

  order = 50
  ShowPlot = True

  print Header+"** Test 1: computing the W_PI values from self-coded and validated Gruzdev theory..."
  timeKeldysh, N_Keldysh_SI, N_Gruzdev_SI = plotPulseToDensity(Egap, meff, wavelength, tau, FieldEnvelope1.real, dt, order, ShowPlot, 0e0, Ntotal)
  
  #print Header+"** Test 0: Convergence test using the Keldysh-Gruzdev formulas..."
  #plotPulseToDensity(Egap, meff, wavelength, tau, PeakFluence, 1E-17, order, ShowPlot)
  #plotPulseToDensity(Egap, meff, wavelength, tau, PeakFluence, 5E-17, order, ShowPlot)
  #plotPulseToDensity(Egap, meff, wavelength, tau, PeakFluence, 1E-16, order, ShowPlot)
  
  #print "Checking dt convergence..."
  #print ""
  #print "Checking order convergence..."
  
  #plotPulseToDensity(Egap, meff, wavelength, tau, PeakFluence, dt, 10, ShowPlot)
  #plotPulseToDensity(Egap, meff, wavelength, tau, PeakFluence, dt, 20, ShowPlot)
  #plotPulseToDensity(Egap, meff, wavelength, tau, PeakFluence, dt, 30, ShowPlot)
  #plotPulseToDensity(Egap, meff, wavelength, tau, PeakFluence, dt, 40, ShowPlot)
  #plotPulseToDensity(Egap, meff, wavelength, tau, PeakFluence, dt, 50, ShowPlot)

  #print Header+"** Test 1: computing the W_PI values from self-coded and validated Gruzdev theory..."
  #timeKeldysh, N_Keldysh_SI, N_Gruzdev_SI = plotPulseToDensity(Egap, meff, wavelength, tau, FieldEnvelope1.real, dt, order, ShowPlot, 0e0, Ntotal)

  #print Header+"** Test 2: computing the W_PI values from Vladimir Zhukov tables..."
  #print Header+"           WE DONT HAVE THEM FOR THIS BAND GAP. Contact zukov@ict.nsc.ru."
  #wPI_Zhukov = VZ_generateWpiTables(FieldEnvelope1, FieldEnvelope2, wavelength, wavelength2, CEP, CEP2, Egap, meff, tau, tau, Delay, dt, Ntotal, t0)
  
  return N_Gruzdev_SI.max()
#}}}

# Extends the Keldysh-Gruzdev model to parameters required by Hamed Merdji group (CEA > LYDIL France) for the nanocone irradiation in ZnO
def ZnOMerdji2017(intensity):#{{{
  print "Defining ZnO material parameters from [Huang2014]..."

  VolumicMass = 5.606e3 #kg/m-3
  MolarMass   = 81.38e-3 #kg/mol
 
  Egap = 3.42*e; #band gap of ZnO [Tsoi2006: Tsoi, S. and Lu, X. and Ramdas, A. K. and Alawadhi, H. and Grimsditch, M. and Cardona, M. and Lauck, R., "Isotopic-mass dependence of the A, B, and C excitonic band gaps in ZnO at low temperatures", Physical Review B (2006).]
  meff = 0.19e0; #Effective mass of ZnO [Huang2014] #TODO: not so serious paper on ZnO! Find a pump probe of ZnO to be more sure. 
  Ntotal=VolumicMass * Avogadro / MolarMass #10.*5E28; #valence band electron density #to avoid limitation
  print Header+"Limiting the excitation degree to Z*=1. Density: "+str(Ntotal*1E-6)+" cm-3"

  wavelength = 3200e-9;
  tau=100e-15; dt = 1E-17; CEP=0e0
  #intensity = 1E12*1E4
  PeakFluence = intensity * tau
  PeakField   = np.sqrt(2e0 * PeakFluence / (tau * c * epsilon_0))

  t0=0. #defines the instant 0.
  Delay = 0. #delay between maxima of the pulses
  tmin=-1.*tau + t0; tmax=1.*tau + Delay + t0

  instants = np.arange(tmin, tmax, dt)
  #print "Time range: "+str(instants.min())+", "+str(instants.max())+"."

  PeakField2  = 0. #/ 2.
  CEP2        = 0. #pi/3.
  #wavelength2 = wavelength
  
  print Header+"** Test: building single pulse centered on 0..."
  FieldEnvelope1, RealField1 = PulseSquaredSinTemporalShape(instants, tau, PeakField, wavelength, CEP, t0, 0.)

  #print Header+"** Test: We build a second pulse with a delay..."
  #FieldEnvelope2, RealField2 = PulseSquaredSinTemporalShape(instants, tau, PeakField2, wavelength2, CEP, t0, Delay)

  #print Header+"** Test: building a bicolor double pulse"

  #FieldEnvelopeTot, RealFieldTot = PulseSquaredSinTemporalShapeDoublePulse(instants, tau, tau, PeakField, PeakField, wavelength, wavelength2, CEP, CEP2, t0, Delay)

  plt.plot(instants, RealField1.real, '-')
  plt.plot(instants, FieldEnvelope1.real, '--')
  #plt.plot(instants, RealField2.real, '-')
  #plt.plot(instants, FieldEnvelope2.real, '--')
  #plt.plot(instants, RealFieldTot.real, '-')
  #plt.plot(instants, FieldEnvelopeTot.real, '--')
  plt.xlabel('')
  plt.savefig('PulseEnvelopes.eps')
  plt.savefig('PulseEnvelopes.png')
  #plt.show()

  print Header+"** Info: PulseEnvelope.EPS and PNG were written in the current folder. "

  order = 50
  ShowPlot = True

  print Header+"** Test 1: computing the W_PI values from self-coded and validated Gruzdev theory..."
  timeKeldysh, N_Keldysh_SI, N_Gruzdev_SI = plotPulseToDensity(Egap, meff, wavelength, tau, FieldEnvelope1.real, dt, order, ShowPlot, 0e0, Ntotal)  

  return N_Gruzdev_SI.max()
#}}}  

