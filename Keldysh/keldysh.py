#!/usr/bin/env python
#-*- coding: utf-8 -*-

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

# Compute the Keldysh excitation

def gammaKeldysh(Egap, meff, Efield, wavelength): #{{{
  #> Adiabadicity parameter
  #> gammaKeldysh < 0.1: means tunneling effect is dominant, 
  #> gammaKeldysh > 10: multi-photon excitation effect is dominant. 
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

def Keldysh1(gamma):
  return gamma/np.sqrt(1.0+gamma**2)

def Keldysh2(gamma):
  return Keldysh1(gamma)/gamma

def EffectiveGap(Egap, k1, k2): #{{{
  if (k1 != 0):
    result = 2.0*Egap * ellipe(k2)/(pi*k1) #BUG in ellipe ? 
  else: 
    print "EffectiveGap(): Keldysh1 = 0."
    exit()
  return result
  # Validation: 
  # EffectiveGap(0, 0, 0) = Error. 
  # EffectiveGap(1., 1., 0) = 1
  # EffectiveGap(1.12*e, 1., 1.) = 0.71*e
#}}}
  
def DawsonIntegral(z): #{{{
  #According to https://docs.scipy.org/doc/scipy-0.14.0/reference/generated/scipy.special.dawsn.html#scipy.special.dawsn
  #Compute int(exp(y**2 - z**2), y=0..z)
  integral = dawsn(z)
  # Validation, compared with Maple. 
  # OK DawsonIntegral(0.)=0
  # OK DawsonIntegral(0.5) = 0.42443638350202229
  return integral
#}}}

def KeldyshFunction(Keldysh1, Keldysh2, Ueff, nmax, wavelength): #{{{
  omegaLaser = 2.*pi*c/wavelength #SI
  n=np.arange(0,nmax+1)
  print n
  sumtable=np.exp(-pi*n*(ellipk(Keldysh1)-ellipe(Keldysh1))/ellipe(Keldysh2))*DawsonIntegral(pi*np.sqrt( (2.0*np.trunc(Ueff/hbar/omegaLaser+1.))-2.0*Ueff/hbar/omegaLaser + n) / (2.0 * ellipk(Keldysh2)*ellipe(Keldysh2)) )
  #print "Effective gap: "+str(Ueff/e)+" eV."
  print sumtable
  result = np.multiply(np.sqrt(pi/(2.*ellipk(Keldysh2))),np.sum(sumtable))
  return result
""" Validation
  KeldyshFunction(Keldysh1=0, Keldysh2=0, Ueff=0, nmax=0, wavelength=0): error. 
  KeldyshFunction(Keldysh1=0, Keldysh2=0, Ueff=0, nmax=0, wavelength=800e-9): 0.
  KeldyshFunction(Keldysh1=0, Keldysh2=0, Ueff=1.12*e, nmax=0, wavelength=800e-9): 0
  
"""
#}}}
Egap = 1.12e0*e; meff=0.18e0; Efield=1E8; wavelength = 800e-9

gamma = gammaKeldysh(Egap, meff, Efield, wavelength) #valid
k1 = Keldysh1(gamma); k2 = Keldysh2(gamma) #valid
EgapEff = EffectiveGap(Egap, k1, k2) #BUG: scipy.special.ellipe (https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.ellipe.html#scipy.special.ellipe) and sympy.mpmath.ellipe (http://docs.sympy.org/0.7.1/modules/mpmath/functions/elliptic.html#legendre-elliptic-integrals) may not give the same result. USE mpmath module instead (http://mpmath.org/). 
order = 10

print "** Info: "
print "gamma = "+str(gamma)
print "k1 = "+str(k1)+", k2 = "+str(k2)+"." #VALID
print "Egap = "+str(Egap/e)+" eV, Ueff = "+str(EgapEff/e)+" eV."
print "PROBLEM ON EFFECTIVE GAP. Maple said: increase. Python said: decrease."

#print "Test on Effective gap function: "+str(EffectiveGap(Egap,k1,k2)/e)
print "Test on direct calculation: "+str(2e0*Egap/(pi*k1) * ellipe(k2)/e) #BUG
print "Ellipe("+str(k2)+")="+str(ellipe(k2))

#print "order = "+str(order)+"."

#wPI = KeldyshFunction( k1, k2, EgapEff, order, wavelength )

#print "Keldysh w_PI until order "+str(order)+" = ", wPI

#def IonizationRate(omegaLaser, Keldysh1, ):
  
  #IonizationRate=(2.*omegaLaser/(9.*pi)*((omegaLaser*m_e)/(hbar*Keldysh1))**(1.5)*KeldyshFunction()*np.exp(-pi*truncate(Ueff/hbar/omegaLaser+1)*((ellipk(Keldysh1)-ellipe(Keldysh1))/(ellipe(Keldysh2))))
	  
  #return IonizationRate
   