#!/usr/bin/env python2
#-*- coding: utf-8 -*-

# Copyright (C) 2013-2019 T. J.-Y. Derrien
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
ShortRefShcheblanov = "[Shcheblanov (2017)]"

## Computes the adiabadicity parameter
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
def gammaKeldysh(Egap, meff, Efield, wavelength, RefractiveIndex=1.): #{{{
  RefractiveIndex_real = RefractiveIndex.real
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
    result = omegaLaser*np.sqrt(m_e*meff*Egap)/e/Efield #Efield must not be ponderated by the optical index in adiabadicity coefficient. Everywhere else, it must be. 
    #print Efield, result
  else:
    ErrorMessage=ErrorMessage+"gamma(): Divergence, as field equals = 0. Singular case of Keldysh functions. Should give w_PI = 0 then...\n"
    result = 1E9 #THIS VALUE IS ARBITRARY FOR A VERY SMALL FIELD. 
  #print(ErrorMessage)
  #print omegaLaser
  return np.float64(result)
# Numerically validated with comparison to Maple. 
#}}}


## Computes the adiabadicity parameter
# Input are same for gammaKeldysh(). 
# Only difference is coefficient 2 in gamma for atoms. For solids, it does not appear. 
def gammaKeldysh_Atoms(Egap, meff, Efield, wavelength, RefractiveIndex=1.): #{{{
  RefractiveIndex_real = RefractiveIndex.real
  #print Egap, meff, Efield
  RecoverSolidGamma=1. #1: no action. 0.5: recovers the solid gamma value
  omegaLaser=2.*pi*c/wavelength
  ErrorMessage=""
  # Validity limit
  if ( Egap < hbar * omegaLaser ): 
    #TODO: this flow should be redirected to an error file. Stdout also goes into the variables. 
    ErrorMessage=ErrorMessage+"** Validity range error: the Keldysh model is not valid for linear absorption. INVALID RESULT...\n"
    ErrorMessage=ErrorMessage+"** Error details: "+str(int(wavelength*1E9))+" nm wavelength is too small for the gap "+str(float(Egap)/e)+".\n"
    #exit() #Avoid to quit, so that octopus still compare its results. 
  if (Efield > 1e-1): #if vectorial, then abs changed its meaning
    result = omegaLaser*np.sqrt(RecoverSolidGamma*2.*m_e*meff*Egap)/e/(Efield) #NOTE: no dependence to optical index for adiabadicity coefficient
    #print Efield, result
  else:
    ErrorMessage=ErrorMessage+"gamma(): Divergence, as field equals = 0. Singular case of Keldysh functions. Should give w_PI = 0 then...\n"
    result = 1E9 #THIS VALUE IS ARBITRARY FOR A VERY SMALL FIELD. 
  #print(ErrorMessage)
  #print omegaLaser
  return np.float64(result)
# Numerically validated with comparison to Maple. 
#}}}

## Computes some intermediate quantity, careful: long double precision.
def Keldysh1phi(gamma): #phi() in Gruzdev2014
  value = np.float128(gamma) #K1(gamma) function has limit 1 when gamma > 5. Hence, we must compute k1(gamma) with a huge precision to stay out of unity. 
  return np.divide(value, np.sqrt( np.float128(1E0) + np.power(value, 2)) )

## Computes some intermediate quantity for Gulley2012, careful: long double precision.
def Keldysh1phiGulley(gamma): #phi() in Gruzdev2014
  value = np.float128(gamma) #K1(gamma) function has limit 1 when gamma > 5. Hence, we must compute k1(gamma) with a huge precision to stay out of unity. 
  return np.divide(np.power(value, 2), np.float128(1E0) + np.power(value, 2))

## Computes some intermediate quantity for Keldysh and Gruzdev models
def Keldysh2theta(gamma): #theta() in Gruzdev2014
  value = np.float128(gamma) #idem about precision.
  #try:
  #except:
    #result = np.float128(0e0)
  #return result
  return np.divide(1E0, np.sqrt( np.float128(1E0) + np.power(value, 2)) )

## Computes intermediate quantity for Gulley model
def Keldysh2thetaGulley(k1): #theta() in Gulley2012
  value = np.float128(k1) #idem about precision.
  #try:
  #except:
    #result = np.float128(0e0)
  return 1. - value
  #return np.add(1E0, -np.sqrt( np.float128(1E0) + np.power(value, 2)) )

## Computes the effective gap for one material
# @param Egap: band gap of the transition (multi-photonic transitions are DIRECT. Tunnel transitions can be INDIRECT)
def EffectiveGap(Egap, k1, k2): #{{{
  k11 = np.float64(k1); #Keldysh1phi()
  k22 = np.float64(k2*k2) #Keldysh2theta() #reducing precision to call ellipe
  if (k11 != 0):
    result = 2.0 / pi * Egap * ellipe(k22)/(k11) #Warning: ellipe(x²) actually computes E(x). 
  else: 
    print "EffectiveGap(): Singular error, Keldysh1phi = 0."
    result = 0e0
  return result
  # Validation: 
  # EffectiveGap(0, 0, 0) = Error. 
  # EffectiveGap(1., 1., 0) = 1
  # EffectiveGap(1.12*e, 1., 1.) = 0.71*e
  # EffectiveGap(Egap,k1,k2)/e = 1.12044049367 #Passed
#}}}

## Taylor development of the effective band gap in Gulley2012 theory.
def EffectiveGapGulley(Egap, Efield, meff, wavelength): #{{{
    omegaLaser = 2.*np.pi*c/wavelength
    return Egap + 0.25 * (e*Efield)**2/(m_e*meff*omegaLaser**2)
#}}}

## Effective potential for excited atom
# @param Potential: in Coulomb, please
# @param wavelength: in meter, please
def EffectiveIonizationPotentialAtom(Potential, meff, wavelength, Field, RefractiveIndex=1.): #{{{
    omegaLaser=2.*np.pi*c/wavelength
    EffectivePotential=Potential+e**2*(Field*np.sqrt(RefractiveIndex.real))**2/(4.*m_e*meff*omegaLaser**2) #this is for Stark effect, where optical refracitve index DOES count. 
    return EffectivePotential
#}}}

## Calculates the Dawson integral int(exp(y**2 - z**2), y=0..z)
# This is based on the Python library SciPy.special functions. 
# According to https://docs.scipy.org/doc/scipy-0.14.0/reference/generated/scipy.special.dawsn.html#scipy.special.dawsn
def DawsonIntegral(z): #{{{
  z2 = np.float64(z)
  integral = dawsn(z2)
  # Validation, compared with Maple. 
  # OK DawsnoIntegral(0.)=0
  # OK DawsonIntegral(0.5) = 0.42443638350202229
  return integral
#}}}

# Generalized definition of hyperbolic arcsinus. 
def arcsinhln(z): #{{{
    return np.log(z+np.sqrt(1.+z**2))
#}}}

## Function S(gamma, x) in Keldysh (1965), employed for atoms
def KeldyshFunction_Atoms(gamma, x, nmax=10): #{{{
    n_tab = np.arange(0,nmax+1)
    expTermL = np.trunc(x+1.)-x+n_tab
    expTermR = arcsinhln(gamma) - gamma/(np.sqrt(1.+np.multiply(gamma, gamma)))
    expTerm  = -2. * expTermL * expTermR
    sumtable = np.exp(expTerm) * DawsonIntegral(
        np.sqrt(2.*gamma/np.sqrt(1.+gamma*gamma) * (np.trunc(x+1)-x+n_tab))
        )
    result   = np.sum(sumtable)
    return result
#}}}

## Function Q(phi, theta) in Keldysh (1965), employed for solids. 
def KeldyshFunction(Keldysh1phi, Keldysh2theta, Ueff, nmax, wavelength): #{{{
  omegaLaser = 2.*pi*c/wavelength #SI
  n_tab = np.arange(0,nmax+1)
 
  # print "Keldysh1phi = "+str(Keldysh1phi)
  Keldysh1phi2_128 = np.power(Keldysh1phi,2)
  Keldysh1phi2 = np.float64(Keldysh1phi2_128)
  # print "--"
  # print "Keldysh1phi2 = "+str(Keldysh1phi2)
  Keldysh2theta2 = np.float64(Keldysh2theta**2)
  # print "Keldysh1phi2 = "+str(Keldysh1phi2)
  distant_to_unity = np.float128(1E0) - Keldysh1phi2_128
  if (distant_to_unity < 1E-320): #then it gonna crash for sure. 
    print "** Error on ellipk: argument 1 is singular. Distance to unit = "+str(distant_to_unity)+"Please increase precision on Keldysh1phi or use ellipkm1 function (careful, argument IS not the same)."
  elif(distant_to_unity < 1E-10): 
    #threshold where functions ellipk and ellipkm1 give different values
    EllipticK1_phi = ellipkm1( np.float64(distant_to_unity) )
  else: #other cases, good for efficiency
    EllipticK1_phi = ellipk( Keldysh1phi2 ) #inf if Keldysh1phi1 = 1.  
  EllipticE1_phi = ellipe( Keldysh1phi2 )
  EllipticE2_theta = ellipe( Keldysh2theta2 )
  # print "EllipticK1_phi = "+str(EllipticK1_phi)
  division = np.divide( EllipticK1_phi - EllipticE1_phi , EllipticE2_theta )
  # print "division= "+str(division)
  exponant = np.multiply( n_tab, division )
  sumtable = np.exp( - pi * exponant )
  sumtable = np.multiply(sumtable,DawsonIntegral(pi*np.sqrt( ((2.0*np.trunc(Ueff/hbar/omegaLaser+1.))-2.0*Ueff/hbar/omegaLaser + n_tab) / (2.0 * ellipk(Keldysh2theta2)*ellipe(Keldysh2theta2)) ) ) ) #/4 K E correction, according to Zhukov
  #print "Effective gap: "+str(Ueff/e)+" eV."
  #print sumtable
  result = np.multiply(np.sqrt(pi/(2.*ellipk(Keldysh2theta2))),np.sum(sumtable))
  return result
""" Validation
  KeldyshFunction(Keldysh1phi=0, Keldysh2theta=0, Ueff=0, nmax=0, wavelength=0): error. 
  KeldyshFunction(Keldysh1phi=0, Keldysh2theta=0, Ueff=0, nmax=0, wavelength=800e-9): 0.
  KeldyshFunction(Keldysh1phi=0, Keldysh2theta=0, Ueff=1.12*e, nmax=0, wavelength=800e-9): 0 
"""
#}}}

## Strange function used by Gulley, not by others. 
def GulleyX(Egap, gamma, Keldysh2theta, wavelength): #{{{
    #Keldysh2theta2 = np.float64(Keldysh2theta**2)
    #EllipticE2_theta = ellipe( Keldysh2theta2 )
    EllipticE2_theta = ellipe( np.float64(Keldysh2theta) ) #NOTE, bug is in Gulley2012: sqrt is necessary if one wants to recover the definition of Keldysh. Gulley2012 made a mistake in the definition of x (only there!). But E(...) implementation requires E(x^2) to work. 
    omegaLaser = 2.*pi*c/wavelength #SI
    #xGulley     = 2.*Egap/(np.pi * omegaLaser) * np.sqrt(1.-gamma**2)/(gamma) * EllipticE2_theta #BUG 1: this generates complex numbers in the MPI regime. It must be sqrt(1+gamma**2), like in Keldysh paper. 
    xGulley     = 2.*Egap/(np.pi * hbar * omegaLaser) * np.sqrt(1.+gamma**2)/(gamma) * EllipticE2_theta #BUG 2: Gulley must have forgotten a hbar. It is needed for dimensional consistency. 
    return xGulley
#}}}

def Gulley_Compute_Elliptics(Keldysh1phi, Keldysh2theta): #{{{
## This section computes K(...) in a safe manner. 
  # print "Keldysh1phi = "+str(Keldysh1phi)
  Keldysh1phi2_128 = np.power(Keldysh1phi,2)
  Keldysh1phi2 = np.float64(Keldysh1phi2_128)
  # print "--"
  # print "Keldysh1phi2 = "+str(Keldysh1phi2)
  Keldysh2theta2_128 = np.power(Keldysh2theta,2)
  Keldysh2theta2 = np.float64(Keldysh2theta**2)
  # print "Keldysh1phi2 = "+str(Keldysh1phi2)
  distant_to_unity = np.float128(1E0) - Keldysh1phi2_128
  distant_to_unity2= np.float128(1E0) - Keldysh2theta2_128
  if (distant_to_unity.all() < 1E-320 or distant_to_unity2.all() < 1E-320): #then it gonna crash for sure. 
    print "** Error on ellipk: argument 1 is singular. Distance to unit = "+str(distant_to_unity)+"Please increase precision on Keldysh1phi or use ellipkm1 function (careful, argument IS not the same)."
  elif(distant_to_unity.all() < 1E-10 or distant_to_unity2.all() < 1E-10): 
    #threshold where functions ellipk and ellipkm1 give different values
    EllipticK1_phi   = ellipkm1( np.float64(distant_to_unity ) )
    EllipticK2_theta = ellipkm1( np.float64(distant_to_unity2) )
  else: #other cases, good for efficiency
    EllipticK1_phi   = ellipk ( Keldysh1phi2   ) #inf if Keldysh1phi1 = 1.  
    EllipticK2_theta = ellipk ( Keldysh2theta2 )
 
  ## This section computes E(...) that are less problematic. 
  EllipticE1_phi   = ellipe( Keldysh1phi2   ) #not used
  EllipticE2_theta = ellipe( Keldysh2theta2 )
  
  return EllipticK1_phi, EllipticK2_theta, EllipticE1_phi, EllipticE2_theta
#}}}

## Keldysh function from Gulley
# NOTE: what is a Keldysh function? Change the name. 
def KeldyshFunctionGulley(Keldysh1phi, Keldysh2theta, xGulley, gamma, nmax, wavelength): #{{{
  omegaLaser = 2.*pi*c/wavelength #SI
  n_tab = np.arange(0,nmax+1)
  
  EllipticK1_phi, EllipticK2_theta, EllipticE1_phi, EllipticE2_theta = Gulley_Compute_Elliptics(Keldysh1phi, Keldysh2theta) #validated by comparison with Maple
  
  division = np.divide( EllipticK1_phi - EllipticE2_theta , EllipticE2_theta ) #should be fine
  
  omegaGulley = np.pi * division
  sumtable1   = np.exp( - n_tab * omegaGulley )
  thetaGulley = np.pi**2/(4. * EllipticK2_theta * EllipticE2_theta) #may be sensitive here
  #xGulley     = 2.*Egap/(np.pi * omegaLaser) * np.sqrt(1.-gamma**2)/(gamma) * EllipticE2_theta
  nu          = np.trunc(xGulley+1.)-xGulley
  sumtable = np.multiply(sumtable1, DawsonIntegral(np.sqrt(thetaGulley * (n_tab + 2. * nu))))
  #pi*np.sqrt( ((2.0*np.trunc(Ueff/hbar/omegaLaser+1.))-2.0*Ueff/hbar/omegaLaser + n_tab) / (4.0 * ellipk(Keldysh2theta2)*ellipe(Keldysh2theta2)) ) ) ) #/4 K E correction, according to Zhukov
  #print "Effective gap: "+str(Ueff/e)+" eV."
  #print sumtable
  result = np.multiply(np.sqrt(pi/(2.*EllipticK2_theta)),np.sum(sumtable))
  return result #validated by exact numerical comparison with Maple
#}}}

## Corrected Keldysh function according to Vitaly Gruzdev (see Ref in details). 
# Kane direct band gap structure
# This function is corrected according to Gruzdev, Opt. Eng. 53, 122515 (2014)
# Only factors of 2 are removed from inside the Keldysh integral, but inserted 
# as a proportion in IonizationRate_Gruzdev function. 
def KeldyshFunction_Gruzdev(Keldysh1phi, Keldysh2theta, Ueff, nmax, wavelength): #{{{
  omegaLaser = 2.*pi*c/wavelength #SI
  n=np.arange(0,nmax)
  Keldysh1phi2_128 = Keldysh1phi**2
  Keldysh1phi2 = np.float64(Keldysh1phi2_128)
  # print "--"
  # print "Keldysh1phi2 = "+str(Keldysh1phi2)
  Keldysh2theta2 = np.float64(Keldysh2theta**2)
  # print "Keldysh1phi2 = "+str(Keldysh1phi2)
  distant_to_unity = 1E0 - Keldysh1phi2_128
  if (distant_to_unity < 1E-320): #then it gonna crash for sure. 
    print "** Error on ellipk: argument 1 is singular. Distance to unit = "+str(distant_to_unity)+"Please increase precision on Keldysh1phi or use ellipkm1 function (careful, argument IS not the same)."
  elif(distant_to_unity < 1E-10):
    #threshold where functions ellipk and ellipkm1 give different values
    EllipticK1_phi = ellipkm1( np.float64(distant_to_unity) )
  else: #other cases, good for efficiency
    EllipticK1_phi = ellipk( Keldysh1phi2 ) #inf if Keldysh1phi2 = 1.i
  EllipticE1_phi   = ellipe( Keldysh1phi2 )  
  EllipticE2_theta = ellipe( Keldysh2theta2 ) 
  #print n
  sumtable=np.exp(-pi*n*(EllipticK1_phi-EllipticE1_phi)/EllipticE2_theta)*DawsonIntegral(pi*np.sqrt( ((np.trunc(Ueff/hbar/omegaLaser+1.))-Ueff/hbar/omegaLaser + n) / (2.0 * ellipk(Keldysh2theta2)*EllipticE2_theta) ) )
  
  # If including Gulley correction...
  #sumtable=np.exp(-pi*n*(EllipticK1_phi-EllipticE1_phi)/EllipticE2_theta)*DawsonIntegral(pi*np.sqrt( ((np.trunc(Ueff/hbar/omegaLaser+1.))-Ueff/hbar/omegaLaser + n) / (4.0 * ellipk(Keldysh2theta2)*EllipticE2_theta) ) )
  #print "Effective gap: "+str(Ueff/e)+" eV."
  #print sumtable
  result = np.multiply(np.sqrt(pi/(2.*ellipk(Keldysh2theta2))), np.sum(sumtable))
  return result
#}}}

## Atomic Keldysh ionization rate [Keldysh 1965, Eq. (1)]
# @param Acoeff: adjustable coefficient, no unit
# @param wavelength: wavelength of the photons in meters, please
# @param atomic_potential: atomic potential in Coulomb, please
# @param field_strength: field amplitude, in V/m
def IonizationRateAtomsEq1(Acoeff, wavelength, meff, atomic_potential, field_strength, nmax=50, RefractiveIndex=1): #{{{
  omegaLaser = 2.*np.pi*c/wavelength
  gamma      = gammaKeldysh_Atoms(atomic_potential, meff, field_strength, wavelength, RefractiveIndex=1)
  EffAtomPotential = EffectiveIonizationPotentialAtom(atomic_potential, meff, wavelength, field_strength, RefractiveIndex.real)
  expTerm    = -2.*EffAtomPotential/hbar/omegaLaser * (
      arcsinhln(gamma)-gamma*np.sqrt(1.+gamma**2)/(1.+2.*gamma**2))
  wAtom      = Acoeff * omegaLaser * np.power(atomic_potential/hbar/omegaLaser, 1.5) * np.power(gamma / np.sqrt(1.+gamma**2),5./2.) * KeldyshFunction_Atoms(gamma, EffAtomPotential/hbar/omegaLaser, nmax) * np.exp(expTerm)
  return wAtom
#}}}


## Atomic Keldysh ionization rate [Keldysh 1965, 2 times the Eq. (16)]
# @param Acoeff: adjustable coefficient, no unit
# @param wavelength: wavelength of the photons in meters, please
# @param atomic_potential: atomic potential in Coulomb, please
# @param field_strength: field amplitude, in V/m
def IonizationRateAtoms(wavelength, meff, atomic_potential, field_strength, nmax=50, RefractiveIndex=1, n0=5E28): #{{{
  omegaLaser = 2.*np.pi*c/wavelength
  gamma      = gammaKeldysh_Atoms(atomic_potential, meff, field_strength, wavelength, RefractiveIndex=1)
  EffAtomPotential = EffectiveIonizationPotentialAtom(atomic_potential, meff, wavelength, field_strength, 1.) #NOTE: should depend on optical refractive index in solids? A priori, no.
  expTerm    = -2.*EffAtomPotential/hbar/omegaLaser * (
      arcsinhln(gamma)-gamma*np.sqrt(1.+gamma**2)/(1.+2.*gamma**2))
  wAtom      = 2.*omegaLaser * np.sqrt(2*atomic_potential/hbar/omegaLaser) * np.power(gamma/np.sqrt(1.+gamma**2), 1.5) * KeldyshFunction_Atoms(gamma, EffAtomPotential/hbar/omegaLaser, nmax) * np.exp(expTerm)
  
  ## DEBUG IonizationRateAtoms
  #Acoeff=1.; wavelength=800e-9; meff=0.2226; atomic_potential=2.56*e; field_strength=10E9; nmax=500; 
  #gamma=gammaKeldysh_Atoms(atomic_potential, meff, field_strength, wavelength, 1)
  #print "Gamma: "+str(gamma)
  #x=EffectiveIonizationPotentialAtom(atomic_potential, meff, wavelength, field_strength) #J, works. 
  #print "Effective Ionization potential (eV): "+str(x/e)
  #print("KeldyshFunction_Atoms: %10.3E"% KeldyshFunction_Atoms(gamma, x, nmax)) 

  #print("KeldyshFunction_Atoms (debug): %10.3E"% KeldyshFunction_Atoms(1., 0., nmax)) 

  #n0=5E28
  #print("Ionization rate (s-1): %10.3E"% IonizationRateAtoms(wavelength, meff, atomic_potential, field_strength, nmax, 1))
  #print "Ionization rate (m^-3 s-1): "+str(n0*IonizationRateAtoms(wavelength, meff, atomic_potential, field_strength, nmax, 1))
  
  return wAtom*n0
#}}}

## The original Keldysh function for Kane direct band gap (solids)
# This function is known to contain mistakes. 
def IonizationRate(Keldysh1phi, Keldysh2theta, KeldyshFunctionResult, Ueff, wavelength, meff=1.0):
  omegaLaser = 2.*np.pi*c/wavelength
  Keldysh1phi2 = np.float64(Keldysh1phi * Keldysh1phi)
  Keldysh2theta2 = np.float64(Keldysh2theta * Keldysh2theta)
  #try:
  result = 2.*omegaLaser/(9.*pi)*((omegaLaser*m_e*meff)/(hbar*Keldysh1phi))**(1.5)*KeldyshFunctionResult*np.exp(-pi*np.trunc(Ueff/hbar/omegaLaser+1e0)*((ellipk(Keldysh1phi2)-ellipe(Keldysh1phi2))/(ellipe(Keldysh2theta2)))) #NOTE: the effective mass in the first term changes everything. It may force us to consider the correction of Gulley for the KeldyshFunction
  #except:
    #print Header+"IonizationRate: ** Error in computation of IonizationRate."
    #result = 0e0
  return result

## Corrected Keldysh photoionization probability according to Vitaly Gruzdev (see Ref in details). 
# Kane direct band gap structure
# This function is corrected according to Gruzdev, Opt. Eng. 53, 122515 (2014)
# Only factor of 2 was added outside the Keldysh integral, but removed
# inside. 
def IonizationRate_Gruzdev(Keldysh1phi, Keldysh2theta, KeldyshFunctionResult, Ueff, wavelength, meff=1.0):
  omegaLaser = 2.*np.pi*c/wavelength
  result =2. * IonizationRate(Keldysh1phi, Keldysh2theta, KeldyshFunctionResult, Ueff, wavelength, meff)
  #2.*omegaLaser/(9.*pi)*((omegaLaser*m_e)/(hbar*Keldysh1phi))**(1.5)*KeldyshFunctionResult*np.exp(-pi*np.trunc(Ueff/hbar/omegaLaser+1)*((ellipk(Keldysh1phi**2)-ellipe(Keldysh1phi**2))/(ellipe(Keldysh2theta**2))))
  
  return result

## The Gulley function for Kane direct band gap structure
# Compared to Keldysh 1965 and Gruzdev 2014, the effective band gap is not exactly apparent. 
def IonizationRate_Gulley(Keldysh1phi, Keldysh2theta, GulleyFunctionResult, xGulley, wavelength, meff=1.0):
  omegaLaser = 2.*np.pi*c/wavelength
  #Keldysh1phi2 = np.float64(Keldysh1phi * Keldysh1phi)
  #Keldysh2theta2 = np.float64(Keldysh2theta * Keldysh2theta)
  #try:
  
  EllipticK1_phi, EllipticK2_theta, EllipticE1_phi, EllipticE2_theta = Gulley_Compute_Elliptics(Keldysh1phi, Keldysh2theta) #OK
  
  division = np.divide( EllipticK1_phi - EllipticE2_theta , EllipticE2_theta ) #OK
  omegaGulley = np.pi * division
  result = 2.*omegaLaser/(9.*pi) * ((omegaLaser*m_e*meff)/(hbar*np.sqrt(Keldysh1phi)))**(1.5)*GulleyFunctionResult*np.exp(-omegaGulley*np.trunc(xGulley+1e0))                                                                                                               
                                                                                                               #*((ellipk(Keldysh1phi2)-ellipe(Keldysh1phi2))/(ellipe(Keldysh2theta2)))) #NOTE: the effective mass in the first term changes everything. It may force us to consider the correction of Gulley for the KeldyshFunction
  #except:
    #print Header+"IonizationRate: ** Error in computation of IonizationRate."
    #result = 0e0
  return result

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

## Generate Keldysh-Gruzdev tables for interfacing with codes
# Input: 
# @param Egap: scalar (J)
# @param meff: scalar (no unit, electron mass is accounted directly in the routine)
# @param wavelength: laser wavelength (scalar, meters)
# @param PeakField: laser field amplitude (scalar, V/m)
# @param order: integration order for Keldysh model (integer, no unit)
# @param RefractiveIndex: in Gruzdev2014, Epeak must be multiplied by sqrt(RefractiveIndex) to EXACTLY repeat his results. This originates that pulse duration is shortened in matter. 
def GenerateKeldyshDatabase(Egap, meff, wavelength, PeakField, order, RefractiveIndex=1): #{{{
  KillStarkEffect=False #WARNING: just for DEBUG purposes: this disables the effective gap for computation of Wpi. 
  ErrorMessage = ""
  # I = 0.5 c epsilon_0 n0 E**2
  # E = np.sqrt(2 I / c / epsilon_0 / n0)
  gamma = gammaKeldysh(Egap, meff, PeakField, wavelength, RefractiveIndex.real) #valid for scalar data
  k1 = Keldysh1phi(gamma); k2 = Keldysh2theta(gamma) #valid
  EgapEff = EffectiveGap(Egap, k1, k2) # Original formula from Keldysh. Warning: scipy.special.ellipe (https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.ellipe.html#scipy.special.ellipe) uses a different convention than Maple, Wikipedia or mpmath.
  if(KillStarkEffect):
      EgapEff=Egap
  KeldyshFunctionResultG = KeldyshFunction_Gruzdev( k1, k2, EgapEff, order, wavelength )
  wPIg = IonizationRate_Gruzdev(k1, k2, KeldyshFunctionResultG, EgapEff, wavelength, meff)
  return wPIg
#}}}

## Generate Keldysh-Gulley tables for interfacing with codes
# Input: 
# @param Egap: scalar (J)
# @param meff: scalar (no unit, electron mass is accounted directly in the routine)
# @param wavelength: laser wavelength (scalar, meters)
# @param PeakField: laser field amplitude (scalar, V/m)
# @param order: integration order for Keldysh model (integer, no unit)
# @param RefractiveIndex: in Gruzdev2014, Epeak must be multiplied by sqrt(RefractiveIndex) to EXACTLY repeat his results. This originates that pulse duration is shortened in matter. 
def GenerateKeldyshGulleyDatabase(Egap, meff, wavelength, PeakField, order, RefractiveIndex=1): #{{{
  ErrorMessage = ""
  # I = 0.5 c epsilon_0 n0 E**2
  # E = np.sqrt(2 I / c / epsilon_0 / n0)
  gamma = gammaKeldysh(Egap, meff, PeakField, wavelength, 1.) #same formula for each, although Gulley does NOT use optical refractive index. 
  k1 = Keldysh1phiGulley(gamma); 
  k2 = Keldysh2thetaGulley(k1)
  EgapEff = EffectiveGapGulley(Egap, PeakField, meff, wavelength) #NOTE: the effective gap of Gulley seems to be only representative. 
  
  xGulley = GulleyX(Egap, gamma, k2, wavelength) #generating the X parameter of Gulley
  
  KeldyshFunctionResultGulley = KeldyshFunctionGulley( k1, k2, xGulley, gamma, order, wavelength )
  wPI = IonizationRate_Gulley(k1, k2,  KeldyshFunctionResultGulley, xGulley, wavelength, meff)
  return wPI, KeldyshFunctionResultGulley, xGulley
#}}}

#print "** Vectorizing functions..."
FieldToIntensity                        = np.vectorize(FieldToIntensity)
IntensityToField                        = np.vectorize(IntensityToField)
gammaKeldysh                            = np.vectorize(gammaKeldysh)
gammaKeldysh_Atoms                      = np.vectorize(gammaKeldysh_Atoms)
Keldysh1phi                             = np.vectorize(Keldysh1phi)
Keldysh2theta                           = np.vectorize(Keldysh2theta)
EffectiveGap                            = np.vectorize(EffectiveGap)
KeldyshFunction                         = np.vectorize(KeldyshFunction)
KeldyshFunction_Gruzdev                 = np.vectorize(KeldyshFunction_Gruzdev)
IonizationRate_Gruzdev                  = np.vectorize(IonizationRate_Gruzdev)
Keldysh1phiGulley                       = np.vectorize(Keldysh1phiGulley)
Keldysh2thetaGulley                     = np.vectorize(Keldysh2thetaGulley)

EffectiveIonizationPotentialAtom        = np.vectorize(EffectiveIonizationPotentialAtom)

EffectiveGapGulley                      = np.vectorize(EffectiveGapGulley)
KeldyshFunctionGulley                   = np.vectorize(KeldyshFunctionGulley)

KeldyshFunction_Atoms                   = np.vectorize(KeldyshFunction_Atoms)

BristowLaw                              = np.vectorize(BristowLaw)
GenerateKeldyshDatabase                 = np.vectorize(GenerateKeldyshDatabase)

## Generate Keldysh excitation rate tables for solids, for a batch of laser fields amplitudes, associated with a wavelength
# Input: 
# @param Egap: scalar (J)
# @param meff: scalar (no unit)
# @param wavelength: laser wavelength (scalar, meters)
# @param tau: total pulse duration (scalar, seconds)
# @param FieldEnvelope: envelope of the laser electric field (scalar|vector, in V/m)
# @param dt: precision (scalar, seconds)
# @param order: integration order for Keldysh model (integer, no unit)
# @N_total: limiter for the ionizable number of electrons (float, m^{-3})
def generateWpiTables(Egap = 2.56e0*e, meff = 0.2226e0, wavelength = 800e-9, tau = 10e-15, FieldEnvelope = 1e9, dt = 1E-17, order = 50, N_total=4*5E28, t0=0.0, TemporalIntegration=False, OpticalIndex=1.): #{{{
                                                                                                                                                                                                                     
  KillStarkEffect=False #WARNING: just for DEBUG purposes: this disables the effective gap for computation of Wpi. 
  
  Header="[libKeldysh] generateWpiTables: "
  print Header+"Defining the laser pulse..."
  # t0 = 0e0
  tmin = -2.*tau+t0
  tmax =  2.*tau+t0
  #dt = 1E-17
  print Header+"** Info: Number of time steps = "+str(int((tmax-tmin)/dt))+"."
  #instants=np.arange(tmin,tmax,dt)
  instants=np.linspace(tmin,tmax,len(FieldEnvelope))

  #Gaussian envelope
  #PulseEnvelope=PulseGaussianTemporalShape(instants, tau, Intensity, t0)
  #FieldEnvelope, PulseEnvelope = PulseSquaredSinTemporalShape(instants, tau, Intensity, wavelength, t0) #TODO: this should be generated outside this function
  IntensityEnvelope=FieldToIntensity(FieldEnvelope)
  print Header+"** Info: Peak intensity = ["+str(np.min(IntensityEnvelope)/1E4)+", "+str(np.max(IntensityEnvelope)/1E4)+"] W/cm^2."
  print Header+"** Info: Peak field amplitude = ["+str(np.min(FieldEnvelope)/1E9)+", "+str(np.max(FieldEnvelope)/1E9)+"] V/nm."

  #print "** Starting the self-consistent loop..."

  #for i in np.arange(1,50,1): #attempt of self consistent loop: divergent
  #print "** ITERATION "+str(i)
  print Header+"Computing Adiabadicity coefficients for the pulse envelope..."
  
  gamma = gammaKeldysh(Egap, meff, FieldEnvelope, wavelength, 1.) #valid for scalar|vector data
  #gamma = gammaKeldysh(EgapEff, meff, IntensityToField(PulseEnvelope), wavelength) #self-consistent, divergent
  print Header+"** Info: Adiabadicity parameter: (min,max) = ("+str(gamma.min())+", "+str(gamma.max())+")."

  #print "Computing Keldysh1phi, Keldysh2theta..."
  k1 = Keldysh1phi(gamma); k2 = Keldysh2theta(gamma) #valid
  EgapEff = EffectiveGap(Egap, k1, k2) # Original formula from Keldysh. Warning: scipy.special.ellipe (https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.ellipe.html#scipy.special.ellipe) uses a different convention than Maple, Wikipedia or mpmath.
  #EgapEff = EffectiveGap(EgapEff, k1, k2) #This formula was made self-consistent, but divergent. 
  #order = 50
  
  EffAtomPotential = EffectiveIonizationPotentialAtom(Egap, meff, wavelength, FieldEnvelope)

  print ""
  print Header+"** Info: Egap = "+str(Egap/e)+" eV." 
  print Header+"** Info: effective Egap: ["+str(np.min(EgapEff/e))+", "+str(np.max(EgapEff)/e)+"] eV."
  print Header+"** Info: effective Atomic potential: ["+str(np.min(EffAtomPotential/e))+", "+str(np.max(EffAtomPotential)/e)+"] eV."
  #print Header+"** Debug info: Keldysh1phi: min,max = ( "+str(np.min(k1))+", "+str(np.max(k1))+"), Keldysh2theta = "+str(k2)
  if(KillStarkEffect):
      EgapEff=Egap
  KeldyshFunctionResult  = KeldyshFunction( k1, k2, EgapEff, order, wavelength )
  KeldyshFunctionResultG = KeldyshFunction_Gruzdev( k1, k2, EgapEff, order, wavelength )
  
  #print Header+"** Debug info: KeldyshFunctionResult: "+str(np.min(KeldyshFunctionResult))+", "+str(np.max(KeldyshFunctionResult))
  
  #print Header+"** Debug info: KeldyshFunctionResultG: "+str(np.min(KeldyshFunctionResultG))+", "+str(np.max(KeldyshFunctionResultG))

  #print "** End of self-consistent loop..."
  print Header+"Computes w_PI^a (Keldysh atoms), w_PI^s (Keldysh solid), and w_PIg (Gruzdev) for the pulse envelope..."
  wPI = IonizationRate(k1, k2, KeldyshFunctionResult, EgapEff, wavelength, meff)
  wPIg = IonizationRate_Gruzdev(k1, k2, KeldyshFunctionResultG, EgapEff, wavelength, meff)
  
  Acoeff = 1. #for now
  #wPIa = IonizationRateAtoms(Acoeff, wavelength, meff, Egap, FieldEnvelope, order, OpticalIndex.real)
  wPIa = IonizationRateAtoms(wavelength, meff, Egap, FieldEnvelope, order, OpticalIndex.real)
  
  print Header+"w_PI (Keldysh atomic) until order "+str(order)+" (min,max) = ", np.min(wPIa), np.max(wPIa)
  print Header+"w_PI (Keldysh solid) until order "+str(order)+" (min,max) = ", np.min(wPI), np.max(wPI)
  print Header+"w_PI (Gruzdev) until order "+str(order)+" (min,max) = ", np.min(wPIg), np.max(wPIg)

  #print "Developing: exporting the table..."
  #exit()
  print ""
  
  if(TemporalIntegration): #{{{
    print Header+"Temporal integration..."
    
    # Temporal integration without limiter
    dN_excited_Keldysh = np.multiply(wPI, dt)
    dN_excited_Gruzdev = np.multiply(wPIg, dt)
    
    # Just multiply array of w_PI by dt, with limited to Ntotal
    N_excited_Keldysh = np.zeros(np.shape(dN_excited_Keldysh))
    N_excited_Gruzdev = np.zeros(np.shape(dN_excited_Gruzdev))
    
    ## Temporal integration using trapeze method to plot Ne(t) 
    for i in np.arange(1,len(dN_excited_Keldysh),1):
        N_excited_Keldysh[i] = N_excited_Keldysh[i-1] + 0.5*(dN_excited_Keldysh[i-1] + dN_excited_Keldysh[i])
        N_excited_Gruzdev[i] = N_excited_Gruzdev[i-1] + 0.5*(dN_excited_Gruzdev[i-1] + dN_excited_Gruzdev[i])
        
    # Trapeze integration method
    N_excited_Keldysh_trapz = np.trapz(wPI, instants)
    N_excited_Gruzdev_trapz = np.trapz(wPIg, instants)
  #}}}
  else: 
    N_excited_Keldysh = 0e0; N_excited_Gruzdev=0e0; N_excited_Keldysh_trapz=0e0; N_excited_Gruzdev_trapz=0e0
  # Temporal integration with limiter
  ##dN_excited_Bristow = np.multiply(BristowLaw(wavelength, Egap)*intensity**2/(2*hbar*omega)), dt)
  
  ## Initial number of electrons in conduction band
  #N_initial = np.zeros(wPI.shape)
  
  #ExpArg_Keldysh = np.divide(dN_excited_Keldysh.cumsum(), N_total)
  #ExpArg_Gruzdev = np.divide(dN_excited_Gruzdev.cumsum(), N_total)
  
  #N_excited_Keldysh = np.multiply( np.exp(-ExpArg_Keldysh), N_total * np.exp(ExpArg_Keldysh) - N_total + N_initial)
  
  #N_excited_Gruzdev = np.multiply( np.exp(-ExpArg_Gruzdev), N_total * np.exp(ExpArg_Gruzdev) - N_total + N_initial)
  return instants, N_excited_Keldysh, N_excited_Gruzdev, gamma, wPI, wPIg, N_excited_Keldysh_trapz, N_excited_Gruzdev_trapz  
#}}}

#generateWpiTables = np.vectorize(generateWpiTables)


## Formula for tunnel ionization in semiconductors [Keldysh 1964, Eq. (60)]
# @param Egap (Joules): direct band gap of the modeled material
# @param meff (adim): effective mass of the conduction band
# @param wavelength (meters): photon wavelength of the excitation
# @param Efield (V/m): electric field envelope of the pulse (scalar|vector)
def KeldyshTunnelingLimit(Egap, meff, wavelength, Efield, EnableStarkEffect=False, OpticalIndex=1.): #{{{
    
    gamma = gammaKeldysh(Egap, meff, Efield, wavelength, 1.) #NOTE: optical index must not affect the tunneling time
    
    k1 = Keldysh1phi(gamma); k2 = Keldysh2theta(gamma) #valid
    
    EgapEff = EffectiveGap(Egap, k1, k2)
    omegaLaser=2.*pi*c/wavelength
    
    #EgapEff = np.clip(EgapEff0, 2*Egap, 10*Egap)
    if(not EnableStarkEffect):
        EgapEff = Egap #For Keldysh1phi965 and Kaiser2000, this line MUST NOT be commented. 
    
    # Direct from Keldysh paper #BUG: overflow in the exp() !!!
    #w_tunnel = 2./9./np.pi**2 * EgapEff / hbar * (m_e*meff*EgapEff/hbar**2)**1.5 * (e*hbar*Efield/(m_e*meff)**0.5/EgapEff**1.5)**2.5  *np.exp(-0.5*np.pi*(m_e*meff)**0.5*EgapEff**1.5/e/hbar/Efield * (1.-1./8.*(m_e*meff)*omegaLaser**2*EgapEff/e**2/Efield**2))
    
    # Taken from Kaiser Phys Rev B, 2000 [not completely validated - gives same result as Keldysh, but still could not obtain Fig. 2 of Kaiser 2000.]
    
    w_tunnel = 2./9./np.pi**2 * EgapEff / hbar * (m_e*meff*EgapEff/hbar**2)**1.5 * (hbar*omegaLaser/EgapEff/gamma)**2.5 * np.exp(-0.5*np.pi*EgapEff*gamma/hbar/omegaLaser * (1.-1./8.*gamma**2))
    
    return w_tunnel
#}}}

## Provides the intensity for which gamma has the given value.
# Useful to normalize the peak field. 
def IntensityAtGamma(gamma, Egap_SI, meff, wavelength, RefractiveIndex=1):
  omega = 2.*np.pi * c / wavelength
  RefractiveIndex_real = 1. #NOTE: tunneling time should not be affected by RefractiveIndex
  #gamma = 1 <=> omega * sqrt(m_e*meff*Egap)/e/Efield/sqrt(RefractiveIndex_real) = 1 <=> m_e*meff*Egap/e**2 = 0.5 * Efield**2*(RefractiveIndex_real)
  Ipeak = meff*m_e*Egap_SI*omega**2*c*epsilon_0 / (2. * e**2 * gamma ** 2) / RefractiveIndex_real #adiab. coeff. does not take optical refractive index
  return Ipeak
    
## Compute and plot density evolution with time using the specific parameters.
# @param Egap (Joules): direct band gap of the modeled material
# @param meff (adim): effective mass of the conduction band
# @param wavelength (meters): photon wavelength of the excitation
# @param tau (seconds): duration of the Gaussian pulse (FWHM)
# @param FieldEnvelope (V/m): electric field envelope of the pulse (scalar|vector)
# @param dt (seconds): precision of the temporal envelope
# @param order (adim): order of the integration (default: 50).
def plotPulseToDensity(Egap = 2.56e0*e, meff = 0.2226e0, wavelength = 800e-9, tau = 10e-15, FieldEnvelope = 1e9, dt = 1E-17, order = 50, ShowPlot=False, t0=0., N_total=4.*5E28): #{{{
  #Header="[libKeldysh] "
  #gamma = gammaKeldysh(Egap, meff, FieldEnvelope, wavelength)
  
  print Header+"======== WARNING: convergence of Ne(t) can be hard to reach. ==========\n"
  
  ### We shall generate the interesting pulse in the file from which we call the Keldysh generator
  instants, N_excited_Keldysh, N_excited_Gruzdev, gamma, wPI, wPIg, N_excited_Keldysh_trapz, N_excited_Gruzdev_trapz = generateWpiTables(Egap, meff, wavelength, tau, FieldEnvelope, dt, order, N_total, t0)
  
  try: 
    sizeGamma = len(gamma)
  except: 
    print Header+"Error: problem on gamma in plotPulseToDensity():"+str(gamma)
  
  print ""
  print Header+"** Warning: results may be not converged."
  print Header+"            Reduce dt, and increase order until convergence."
  print ""
  
  # Direct sum
  print Header+"Maximum density N_ex "+ShortRefKeldysh+" = "+str(np.max(N_excited_Keldysh))+"."
  print Header+"Maximum density N_ex "+ShortRefGruzdev+" = "+str(np.max(N_excited_Gruzdev))+"."
  
  # Trapeze integration rule
  print Header+"Maximum density N_ex "+ShortRefKeldysh+" = "+str(np.max(N_excited_Keldysh_trapz))+"."
  print Header+"Maximum density N_ex "+ShortRefGruzdev+" = "+str(np.max(N_excited_Gruzdev_trapz))+"."
  print ""
#  return N_excited_Gruzdev 
  if(ShowPlot): 
    print Header+"Plotting as function of time..."

    xunit = 1E15
    timeunit = "fs"

    #plt.figure()
    plt.figure(figsize=(15,15))
    
    plt.subplot(311)
    plt.title(r"Gap = "+str(Egap/e)+" eV, $\lambda=$ "+str(wavelength*1E9)+r" nm, $\tau=$"+str(tau*xunit)+" "+timeunit+", "+r"$E_{max}=$"+str(np.max(FieldEnvelope)/1E9)+" V/nm")
    #plt.xlabel("Field (V/m)")
    #plt.xlabel("Time (ps)")
    plt.ylabel("Field envelope (V/m)")
    plt.plot(xunit*instants[0:sizeGamma], FieldEnvelope, color="k", linestyle="-", label=r"Field envelope")
    plt.grid()
    plt.legend(loc=3)
    
    plt.subplot(312)
    #plt.xlabel("Field (V/m)")
    #plt.xlabel("Time (ps)")
    plt.ylabel("Adiabadicity $\gamma$")
    plt.semilogy(xunit*instants[0:sizeGamma], gamma, color="k", linestyle="-", label=r"$\gamma$")
    plt.grid()
    plt.legend(loc=3)

    plt.subplot(313)
    #plt.xlabel("Field (V/m)")
    #plt.xlabel("Time (ps)")
    plt.ylabel("$w_{PI}$ (m$^{-3}$ s$^{-1}$)")
    plt.plot(instants[0:sizeGamma]*xunit, wPI, linestyle="-", color="r", label=r"$w_{PI}$ "+ShortRefKeldysh)
    plt.plot(instants[0:sizeGamma]*xunit, wPIg, linestyle="-", color="b", label=r"$w_{PI}$ "+ShortRefGruzdev)
    plt.grid()
    plt.legend(loc=2)

    #plt.subplot(414)
    
    #plt.xlabel("Time ("+timeunit+")")
    ## plt.xlabel("Field (V/m)")
    #plt.ylabel("Density (m$^{-3}$)")
    #plt.plot(instants[0:sizeGamma]*xunit, N_excited_Keldysh, color="r", label="$n_e$ "+ShortRefKeldysh)
    #plt.plot(instants[0:sizeGamma]*xunit, N_excited_Gruzdev, color="b", label="$n_e$ "+ShortRefGruzdev)
    #plt.grid()
    
    plt.legend(loc=2)
    plt.tight_layout()
    plt.savefig("KeldyshSimple.eps") 
    
    ### Second plot
    #print Header+"Importing Gulley [2012] data..."
    #try:
      ##Gulley2012=np.loadtxt("Gulley-Fig2.csv", dtype='float', delimiter='\t')
      ##Gulley2012=np.loadtxt("Gruzdev2014-Fig1.csv", dtype='float', delimiter=',')
      ##print Gulley2012[:,0]
    #except: 
      #print Header+"** Warning: failed to import Gulley2012 data table..."
      
    
    #print Header+"Plotting as function of laser field intensity ..."
    #plt.figure()
    #plt.xlabel("Intensity (W/cm$^{2}$)")
    #plt.ylabel("$w_{PI}$ (cm$^{-3}$ fs$^{-1}$)")
    #plt.loglog(1e-4*FieldToIntensity(FieldEnvelope.real), 1e-6*1e-15*wPI,  linestyle="-", color="r", label=r"$w_{PI}$ "+ShortRefKeldysh)
    #plt.loglog(1e-4*FieldToIntensity(FieldEnvelope.real), 1e-6*1e-15*wPIg, linestyle="-", color="b", label=r"$w_{PI}$ "+ShortRefGruzdev)
    #plt.loglog(Gulley2012[:,0], Gulley2012[:,1], linestyle="-", color="k", label="Data from "+ShortRefGulley) #JUST FOR VALIDATION. 
    #plt.grid()
    #plt.legend(loc=2)
    #plt.xlim((1E10, 1E14))
    ##plt.ylim((1E20*1E6,1E40*1E6))
    #plt.tight_layout()
    #plt.savefig("Keldysh-Field-Wpi.eps")
    #plt.show()
    
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
  tmin=-4.*tau + t0; tmax=4.*tau + Delay + t0

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

