#!/usr/bin/env python
#-*- coding: utf-8 -*-

# IMPORT LIBRARIES
import numpy as np
from numpy import genfromtxt, loadtxt, chararray
#from scipy.optimize import fsolve, root
from scipy.special import ellipk, ellipe
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

def gamma(Egap, meff, Efield, omegaLaser): #{{{
  #> Adiabadicity parameter
  #> gamma < 0.1: means tunneling effect is dominant, 
  #> gamma > 10: multi-photon excitation effect is dominant. 
  return(omegaLaser*sqrt(m_e*Egap*meff)/e/Efield)
#}}}

def Keldysh1(gamma):
  return gamma/sqrt(1.0+gamma**2)

def Keldysh2(gamma):
  return Keldysh1(gamma)/gamma

def EffectiveGap(Egap):
  return 2.0*Egap/(pi*Keldysh1)*ellipe(Keldysh2)
  
def KeldyshFunction2(y,z): #{{{
  #TODO: VALIDATE
  #Compute int(y**2 - z**2, y=0..z)
  samples=10000
  dy=z/float(samples)
  integral=0.
  for i in np.arange(0,samples):
    integral=integral+np.exp(y**2-z**2)*dy #TODO: error in treating y
  return integral
#}}}

def KeldyshFunction():
  return sqrt(pi/(2.*ellipk(Keldysh2)))*np.sum(np.exp(-pi*n*(ellipk(Keldysh1)-ellipe(Keldysh1))/ellipe(Keldysh2))*KeldyshFunction2(pi*sqrt( (2.0*truncate(Ueff/hbar/omegaLaser+1.))-2.0*Ueff/hbar/omegaLaser + n) / (2.0 * ellipk(Keldysh2)*ellipe(Keldysh2)) ) )

def IonizationRate(omegaLaser, Keldysh1, ):
  
  IonizationRate=(2.*omegaLaser/(9.*pi)*((omegaLaser*m_e)/(hbar*Keldysh1))**(1.5)*KeldyshFunction()*np.exp(-pi*truncate(Ueff/hbar/omegaLaser+1)*((ellipk(Keldysh1)-ellipe(Keldysh1))/(ellipe(Keldysh2))))
	  
  return IonizationRate
   