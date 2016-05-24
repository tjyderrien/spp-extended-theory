
#!/usr/bin/env python
#-*- coding: utf-8 -*-

# IMPORT LIBRARIES
import numpy as np
from numpy import genfromtxt, loadtxt, chararray
from scipy.optimize import fsolve, root
import cmath
import matplotlib as mp
import matplotlib.pyplot as plt
from scipy.interpolate import InterpolatedUnivariateSpline
from matplotlib import rc
# from pylab import *
from scipy.constants import c, epsilon_0, mu_0, pi, e, m_e, h
from matplotlib.legend_handler import HandlerLine2D

lengthunit = 1e-9
eta = 5e0 #assumed precision error on the dielectric permittivity

# Settings for matplotlib
#rc('font',**{'family':'sans-serif','sans-serif':['Helvetica'], 'size':'16'})
## for Palatino and other serif fonts use:
rc('font', **{'family':'serif', 'serif':['Palatino'], 'size':'20'})
rc('text', usetex=True)
mp.rcParams['legend.numpoints'] = 1

# basic wave function
def omega(wavelength):#{{{
    return 2.0*pi*c/wavelength
#}}}

# SPP BASIC FUNCTIONS
def betaSPP(wavelength, eps1, eps2):#{{{
    """calculate the SPP wave number on a flat interface
    input: wavelength (float), eps1 (complex), eps2(complex)
    """
    omega = 2.0*pi*c/wavelength
    try:
        value = omega/c * cmath.sqrt(eps1 * eps2 / (eps1 + eps2))
    except: 
        print "betaSPP: singular case, error code: -1"
        value = -1e0+0e0j
    return value
#}}}

def AsymmetricSPPconditionPos(eps1, eps2):#{{{
    """ Assume that Re(k1).Re(k2) > 0 and verify the subsequent consequences alltogether.
    It exists then two sub-modes, let's say a positive one (Im k1<0, Im k2>0), and a negative one (Im k1>0, Im k2<0)
    This routine is about: Im(k1)>0, Im(k2)<0
    """
    #value = eps1/eps2
    value = eps1.imag * eps2.real - eps1.real * eps2.imag
    condition = (value.imag < 0e0)
    #condition = (eps1.imag / eps2.imag * eps2.real < eps1.real)
    return condition
#}}}

def AsymmetricSPPconditionNeg(eps1, eps2):#{{{
    """ Assume that Re(k1).Re(k2) > 0 and verify the subsequent consequences alltogether.
    It exists then two sub-modes, let's say a positive one (Im k1<0, Im k2>0), and a negative one (Im k1>0, Im k2<0)
    This routine is about: Im(k1)<0, Im(k2)>0
    """
    #value = eps1/eps2
    value = eps1.imag * eps2.real - eps1.real * eps2.imag
    condition = (value.imag > 0e0)
    #condition = (eps1.imag / eps2.imag * eps2.real < eps1.real)
    return condition
#}}}

def SPPconditionValue(eps1, eps2):#{{{
  """SPPconditionValue() returns the value of condition for SPP. If its negative, then SPP can be excited at a flat interface. 
  /!\ This condition is restricted to checking the real part of the dispersion relation for symmetric SPP only. 
    Input: eps1, eps2: complex-valued quantities
    Output: float
  """
  condition=eps1.real*eps2.real+eps1.imag*eps2.imag
  return condition
#}}}

def SPPcondition(eps1, eps2):#{{{
  """ Returns a boolean claiming if SPP are excitable on an interface
  """
  if (SPPconditionValue(eps1, eps2) < 0.0):
           output=True
  else:
    output=False
  return output
#}}}

def OldSPPcondition(eps1, eps2):#{{{
  """SPPconditionValue() returns the value of condition for SPP IN PERFECT MATERIALS (Im(eps)<<|Re(eps)). If its negative, then SPP can be excited at a flat interface. 
    Input: eps1, eps2: complex-valued quantities
    Output: float
  """
  condition1=(eps1.real*eps2.real<0.0)
  #condition2= eps2.real < abs(eps1.real) #this version is not symmetric, hence strange
  # Let's use its generalization which is actually symmetric. 
  condition2 = (eps1.real * eps2.real / (eps1.real + eps2.real) > 0e0)
  return (condition1 and condition2)
#}}}

def period(betaSPP):#{{{
  """ Returns the period of the light-SPP field at a given interface
  """
  try:
    result = 2.0*pi/betaSPP.real
  except: 
    result = -1e0
  return result
#}}}

# Precision over knowledge of period
def deltaBetaSPP(wavelength, eps1, eps2, deps1r, deps1c, deps2r, deps2c):#{{{
  """Calculates the precision over Re(beta) using uncertainty calculations
  """
  eps1r = eps1.real; eps1c = eps1.imag
  eps2r = eps2.real; eps2c = eps2.imag
  
  try:
    term1 = 1e0/2e0*pi*abs(1e0/(2e0*(((eps1r*eps2r-eps1c*eps2c)*(eps1r+eps2r)/((eps1r+eps2r)**2+
            (eps1c+eps2c)**2)+(eps1c*eps2r+eps1r*eps2c)*(eps1c+eps2c)/((eps1r+eps2r)**2+
            (eps1c+eps2c)**2))**2+((eps1c*eps2r+eps1r*eps2c)*(eps1r+eps2r)/((eps1r+eps2r)**2+
            (eps1c+eps2c)**2)-(eps1r*eps2r-eps1c*eps2c)*(eps1c+eps2c)/((eps1r+eps2r)**2+
            (eps1c+eps2c)**2))**2)**(1e0/2e0)+(2e0*eps1r*eps2r-2e0*eps1c*eps2c)*
            (eps1r+eps2r)/((eps1r+eps2r)**2+(eps1c+eps2c)**2)+(2e0*eps1c*eps2r+2e0*eps1r*eps2c)*
            (eps1c+eps2c)/((eps1r+eps2r)**2+(eps1c+eps2c)**2))**(1e0/2e0)/wavelength*(1e0/
            (((eps1r*eps2r-eps1c*eps2c)*(eps1r+eps2r)/((eps1r+eps2r)**2+(eps1c+eps2c)**2)+
            (eps1c*eps2r+eps1r*eps2c)*(eps1c+eps2c)/((eps1r+eps2r)**2+(eps1c+eps2c)**2))**2+
            ((eps1c*eps2r+eps1r*eps2c)*(eps1r+eps2r)/((eps1r+eps2r)**2+(eps1c+eps2c)**2)-
             (eps1r*eps2r-eps1c*eps2c)*(eps1c+eps2c)/((eps1r+eps2r)**2+(eps1c+eps2c)**2))**2)**(1e0/2e0)*
            (2e0*((eps1r*eps2r-eps1c*eps2c)*(eps1r+eps2r)/((eps1r+eps2r)**2+(eps1c+eps2c)**2)+
            (eps1c*eps2r+eps1r*eps2c)*(eps1c+eps2c)/((eps1r+eps2r)**2+(eps1c+eps2c)**2))*
            (eps2r*(eps1r+eps2r)/((eps1r+eps2r)**2+(eps1c+eps2c)**2)+
            (eps1r*eps2r-eps1c*eps2c)/((eps1r+eps2r)**2+(eps1c+eps2c)**2)
            -(eps1r*eps2r-eps1c*eps2c)*(eps1r+eps2r)/((eps1r+eps2r)**2
            +(eps1c+eps2c)**2)**2*(2*eps1r+2*eps2r)+eps2c*(eps1c+eps2c)/((eps1r+eps2r)**2
            +(eps1c+eps2c)**2)-(eps1c*eps2r+eps1r*eps2c)*(eps1c+eps2c)/((eps1r+eps2r)**2
            +(eps1c+eps2c)**2)**2*(2*eps1r+2*eps2r))+2e0*((eps1c*eps2r+eps1r*eps2c)
            *(eps1r+eps2r)/((eps1r+eps2r)**2+(eps1c+eps2c)**2)-(eps1r*eps2r-eps1c*eps2c)
            *(eps1c+eps2c)/((eps1r+eps2r)**2+(eps1c+eps2c)**2))
            *(eps2c*(eps1r+eps2r)/((eps1r+eps2r)**2+(eps1c+eps2c)**2)
            +(eps1c*eps2r+eps1r*eps2c)/((eps1r+eps2r)**2+(eps1c+eps2c)**2)
            -(eps1c*eps2r+eps1r*eps2c)*(eps1r+eps2r)/((eps1r+eps2r)**2
            +(eps1c+eps2c)**2)**2*(2e0*eps1r+2e0*eps2r)-eps2r*(eps1c+eps2c)/((eps1r+eps2r)**2
            +(eps1c+eps2c)**2)+(eps1r*eps2r-eps1c*eps2c)*(eps1c+eps2c)/((eps1r+eps2r)**2
            +(eps1c+eps2c)**2)**2*(2*eps1r+2*eps2r)))+2*eps2r*(eps1r+eps2r)/((eps1r+eps2r)**2
            +(eps1c+eps2c)**2)+(2*eps1r*eps2r-2*eps1c*eps2c)/((eps1r+eps2r)**2+(eps1c+eps2c)**2)
            -(2*eps1r*eps2r-2*eps1c*eps2c)*(eps1r+eps2r)
            /((eps1r+eps2r)**2+(eps1c+eps2c)**2)**2*(2*eps1r+2*eps2r)+2*eps2c*(eps1c+eps2c)
            /((eps1r+eps2r)**2+(eps1c+eps2c)**2)-(2*eps1c*eps2r+2*eps1r*eps2c)*(eps1c+eps2c)
            /((eps1r+eps2r)**2+(eps1c+eps2c)**2)**2*(2*eps1r+2*eps2r)))*deps1r
            
    term2 = (1e0/2e0)*pi*abs((((2e0*((eps1r*eps2r-eps1c*eps2c)*(eps1r+eps2r)/((eps1r+eps2r)**2+(eps1c+eps2c)**2)
            +(eps1c*eps2r+eps1r*eps2c)*(eps1c+eps2c)/((eps1r+eps2r)**2+(eps1c+eps2c)**2)))*(-eps2c*(eps1r+eps2r)/
            ((eps1r+eps2r)**2+(eps1c+eps2c)**2)-(eps1r*eps2r-eps1c*eps2c)*(eps1r+eps2r)*
            (2*eps1c+2*eps2c)/((eps1r+eps2r)**2+(eps1c+eps2c)**2)**2+eps2r*(eps1c+eps2c)/((eps1r+eps2r)**2
            +(eps1c+eps2c)**2)+(eps1c*eps2r+eps1r*eps2c)/((eps1r+eps2r)**2+(eps1c+eps2c)**2)-(eps1c*eps2r+eps1r*eps2c)
            *(eps1c+eps2c)*(2*eps1c+2*eps2c)/((eps1r+eps2r)**2+(eps1c+eps2c)**2)**2)
            +(2e0*((eps1c*eps2r+eps1r*eps2c)*(eps1r+eps2r)/((eps1r+eps2r)**2+(eps1c+eps2c)**2)
            -(eps1r*eps2r-eps1c*eps2c)*(eps1c+eps2c)/((eps1r+eps2r)**2+(eps1c+eps2c)**2)))
            *(eps2r*(eps1r+eps2r)/((eps1r+eps2r)**2+(eps1c+eps2c)**2)-(eps1c*eps2r+eps1r*eps2c)*(eps1r+eps2r)
            *(2*eps1c+2*eps2c)/((eps1r+eps2r)**2+(eps1c+eps2c)**2)**2+eps2c*(eps1c+eps2c)/((eps1r+eps2r)**2
            +(eps1c+eps2c)**2)-(eps1r*eps2r-eps1c*eps2c)/((eps1r+eps2r)**2+(eps1c+eps2c)**2)+(eps1r*eps2r-eps1c*eps2c)
            *(eps1c+eps2c)*(2*eps1c+2*eps2c)/((eps1r+eps2r)**2+(eps1c+eps2c)**2)**2))/np.sqrt(((eps1r*eps2r-eps1c*eps2c)
            *(eps1r+eps2r)/((eps1r+eps2r)**2+(eps1c+eps2c)**2)+(eps1c*eps2r+eps1r*eps2c)*(eps1c+eps2c)/((eps1r+eps2r)**2
            +(eps1c+eps2c)**2))**2+((eps1c*eps2r+eps1r*eps2c)*(eps1r+eps2r)/((eps1r+eps2r)**2+(eps1c+eps2c)**2)
            -(eps1r*eps2r-eps1c*eps2c)*(eps1c+eps2c)/((eps1r+eps2r)**2+(eps1c+eps2c)**2))**2)-2*eps2c*(eps1r+eps2r)/
            ((eps1r+eps2r)**2+(eps1c+eps2c)**2)-(2*eps1r*eps2r-2*eps1c*eps2c)*(eps1r+eps2r)*(2*eps1c+2*eps2c)/
            ((eps1r+eps2r)**2+(eps1c+eps2c)**2)**2+2*eps2r*(eps1c+eps2c)/((eps1r+eps2r)**2+(eps1c+eps2c)**2)+
            (2*eps1c*eps2r+2*eps1r*eps2c)/((eps1r+eps2r)**2+(eps1c+eps2c)**2)-(2*eps1c*eps2r+2*eps1r*eps2c)*
            (eps1c+eps2c)*(2*eps1c+2*eps2c)/((eps1r+eps2r)**2+(eps1c+eps2c)**2)**2)/(np.sqrt(2e0*np.sqrt(((eps1r*eps2r
            -eps1c*eps2c)*(eps1r+eps2r)/((eps1r+eps2r)**2+(eps1c+eps2c)**2)+(eps1c*eps2r+eps1r*eps2c)*(eps1c+eps2c)/
            ((eps1r+eps2r)**2+(eps1c+eps2c)**2))**2+((eps1c*eps2r+eps1r*eps2c)*(eps1r+eps2r)/((eps1r+eps2r)**2+
            (eps1c+eps2c)**2)-(eps1r*eps2r-eps1c*eps2c)*(eps1c+eps2c)/((eps1r+eps2r)**2+(eps1c+eps2c)**2))**2)+
            (2*eps1r*eps2r-2*eps1c*eps2c)*(eps1r+eps2r)/((eps1r+eps2r)**2+(eps1c+eps2c)**2)+(2*eps1c*eps2r+2*eps1r*eps2c)
            *(eps1c+eps2c)/((eps1r+eps2r)**2+(eps1c+eps2c)**2))*wavelength))*deps1c
            
    term3 = 1e0/2e0*pi*abs(1e0/(2e0*(((eps1r*eps2r-eps1c*eps2c)*(eps1r+eps2r)/((eps1r+eps2r)**2+(eps1c+eps2c)**2)+
            (eps1c*eps2r+eps1r*eps2c)*(eps1c+eps2c)/((eps1r+eps2r)**2+(eps1c+eps2c)**2))**2+((eps1c*eps2r+eps1r*eps2c)
            *(eps1r+eps2r)/((eps1r+eps2r)**2+(eps1c+eps2c)**2)-(eps1r*eps2r-eps1c*eps2c)*(eps1c+eps2c)/((eps1r+eps2r)**2
            +(eps1c+eps2c)**2))**2)**(1e0/2e0)+(2e0*eps1r*eps2r-2e0*eps1c*eps2c)*(eps1r+eps2r)/((eps1r+eps2r)**2+(eps1c+eps2c)**2)
            +(2e0*eps1c*eps2r+2e0*eps1r*eps2c)*(eps1c+eps2c)/((eps1r+eps2r)**2+(eps1c+eps2c)**2))**(1e0/2e0)/wavelength
            *(1e0/(((eps1r*eps2r-eps1c*eps2c)*(eps1r+eps2r)/((eps1r+eps2r)**2+(eps1c+eps2c)**2)+(eps1c*eps2r+eps1r*eps2c)
            *(eps1c+eps2c)/((eps1r+eps2r)**2+(eps1c+eps2c)**2))**2+((eps1c*eps2r+eps1r*eps2c)*(eps1r+eps2r)/((eps1r+eps2r)**2
            +(eps1c+eps2c)**2)-(eps1r*eps2r-eps1c*eps2c)*(eps1c+eps2c)/((eps1r+eps2r)**2+(eps1c+eps2c)**2))**2)**(1e0/2e0)
            *(2e0*((eps1r*eps2r-eps1c*eps2c)*(eps1r+eps2r)/((eps1r+eps2r)**2+(eps1c+eps2c)**2)+(eps1c*eps2r+eps1r*eps2c)
            *(eps1c+eps2c)/((eps1r+eps2r)**2+(eps1c+eps2c)**2))*(eps1r*(eps1r+eps2r)/((eps1r+eps2r)**2+(eps1c+eps2c)**2)
            +(eps1r*eps2r-eps1c*eps2c)/((eps1r+eps2r)**2+(eps1c+eps2c)**2)-(eps1r*eps2r-eps1c*eps2c)*(eps1r+eps2r)/
            ((eps1r+eps2r)**2+(eps1c+eps2c)**2)**2*(2e0*eps1r+2e0*eps2r)+eps1c*(eps1c+eps2c)/((eps1r+eps2r)**2
            +(eps1c+eps2c)**2)-(eps1c*eps2r+eps1r*eps2c)*(eps1c+eps2c)/((eps1r+eps2r)**2+(eps1c+eps2c)**2)**2*
            (2*eps1r+2*eps2r))+2*((eps1c*eps2r+eps1r*eps2c)*(eps1r+eps2r)/((eps1r+eps2r)**2+(eps1c+eps2c)**2)
            -(eps1r*eps2r-eps1c*eps2c)*(eps1c+eps2c)/((eps1r+eps2r)**2+(eps1c+eps2c)**2))*(eps1c*(eps1r+eps2r)/((eps1r+eps2r)**2
            +(eps1c+eps2c)**2)+(eps1c*eps2r+eps1r*eps2c)/((eps1r+eps2r)**2+(eps1c+eps2c)**2)-(eps1c*eps2r+eps1r*eps2c)*
            (eps1r+eps2r)/((eps1r+eps2r)**2+(eps1c+eps2c)**2)**2*(2*eps1r+2*eps2r)-eps1r*(eps1c+eps2c)/((eps1r+eps2r)**2
            +(eps1c+eps2c)**2)+(eps1r*eps2r-eps1c*eps2c)*(eps1c+eps2c)/((eps1r+eps2r)**2+(eps1c+eps2c)**2)**2
            *(2*eps1r+2*eps2r)))+2*eps1r*(eps1r+eps2r)/((eps1r+eps2r)**2+(eps1c+eps2c)**2)+(2*eps1r*eps2r-2*eps1c*eps2c)
            /((eps1r+eps2r)**2+(eps1c+eps2c)**2)-(2*eps1r*eps2r-2*eps1c*eps2c)*(eps1r+eps2r)/((eps1r+eps2r)**2
            +(eps1c+eps2c)**2)**2*(2*eps1r+2*eps2r)+2*eps1c*(eps1c+eps2c)/((eps1r+eps2r)**2+(eps1c+eps2c)**2)
            -(2*eps1c*eps2r+2*eps1r*eps2c)*(eps1c+eps2c)/((eps1r+eps2r)**2+(eps1c+eps2c)**2)**2*(2*eps1r+2*eps2r)))*deps2r
            
    term4 = 1e0/2e0*pi*abs(1e0/(2e0*(((eps1r*eps2r-eps1c*eps2c)*(eps1r+eps2r)/((eps1r+eps2r)**2+(eps1c+eps2c)**2)
            +(eps1c*eps2r+eps1r*eps2c)*(eps1c+eps2c)/((eps1r+eps2r)**2+(eps1c+eps2c)**2))**2+((eps1c*eps2r+eps1r*eps2c)
            *(eps1r+eps2r)/((eps1r+eps2r)**2+(eps1c+eps2c)**2)-(eps1r*eps2r-eps1c*eps2c)*(eps1c+eps2c)/((eps1r+eps2r)**2
            +(eps1c+eps2c)**2))**2)**(1e0/2e0)+(2*eps1r*eps2r-2*eps1c*eps2c)*(eps1r+eps2r)/((eps1r+eps2r)**2+(eps1c+eps2c)**2)
            +(2*eps1c*eps2r+2*eps1r*eps2c)*(eps1c+eps2c)/((eps1r+eps2r)**2+(eps1c+eps2c)**2))**(1e0/2e0)/wavelength*(1e0/
            (((eps1r*eps2r-eps1c*eps2c)*(eps1r+eps2r)/((eps1r+eps2r)**2+(eps1c+eps2c)**2)+(eps1c*eps2r+eps1r*eps2c)
            *(eps1c+eps2c)/((eps1r+eps2r)**2+(eps1c+eps2c)**2))**2+((eps1c*eps2r+eps1r*eps2c)*(eps1r+eps2r)/((eps1r+eps2r)**2
            +(eps1c+eps2c)**2)-(eps1r*eps2r-eps1c*eps2c)*(eps1c+eps2c)/((eps1r+eps2r)**2+(eps1c+eps2c)**2))**2)**(1e0/2e0)
            *(2*((eps1r*eps2r-eps1c*eps2c)*(eps1r+eps2r)/((eps1r+eps2r)**2+(eps1c+eps2c)**2)+(eps1c*eps2r+eps1r*eps2c)*
            (eps1c+eps2c)/((eps1r+eps2r)**2+(eps1c+eps2c)**2))*(-eps1c*(eps1r+eps2r)/((eps1r+eps2r)**2+(eps1c+eps2c)**2)
            -(eps1r*eps2r-eps1c*eps2c)*(eps1r+eps2r)/((eps1r+eps2r)**2+(eps1c+eps2c)**2)**2*(2*eps1c+2*eps2c)+eps1r*
            (eps1c+eps2c)/((eps1r+eps2r)**2+(eps1c+eps2c)**2)+(eps1c*eps2r+eps1r*eps2c)/((eps1r+eps2r)**2+(eps1c+eps2c)**2)
            -(eps1c*eps2r+eps1r*eps2c)*(eps1c+eps2c)/((eps1r+eps2r)**2+(eps1c+eps2c)**2)**2*(2*eps1c+2*eps2c))
            +2*((eps1c*eps2r+eps1r*eps2c)*(eps1r+eps2r)/((eps1r+eps2r)**2+(eps1c+eps2c)**2)-(eps1r*eps2r-eps1c*eps2c)*
            (eps1c+eps2c)/((eps1r+eps2r)**2+(eps1c+eps2c)**2))*(eps1r*(eps1r+eps2r)/((eps1r+eps2r)**2+(eps1c+eps2c)**2)
            -(eps1c*eps2r+eps1r*eps2c)*(eps1r+eps2r)/((eps1r+eps2r)**2+(eps1c+eps2c)**2)**2*(2*eps1c+2*eps2c)+eps1c*
            (eps1c+eps2c)/((eps1r+eps2r)**2+(eps1c+eps2c)**2)-(eps1r*eps2r-eps1c*eps2c)/((eps1r+eps2r)**2+(eps1c+eps2c)**2)
            +(eps1r*eps2r-eps1c*eps2c)*(eps1c+eps2c)/((eps1r+eps2r)**2+(eps1c+eps2c)**2)**2*(2*eps1c+2*eps2c)))-2*eps1c*
            (eps1r+eps2r)/((eps1r+eps2r)**2+(eps1c+eps2c)**2)-(2*eps1r*eps2r-2*eps1c*eps2c)*(eps1r+eps2r)/((eps1r+eps2r)**2
            +(eps1c+eps2c)**2)**2*(2*eps1c+2*eps2c)+2*eps1r*(eps1c+eps2c)/((eps1r+eps2r)**2+(eps1c+eps2c)**2)+
            (2*eps1c*eps2r+2*eps1r*eps2c)/((eps1r+eps2r)**2+(eps1c+eps2c)**2)-(2*eps1c*eps2r+2*eps1r*eps2c)*(eps1c+eps2c)/
            ((eps1r+eps2r)**2+(eps1c+eps2c)**2)**2*(2*eps1c+2*eps2c)))*deps2c
            
    deltaReBeta = abs(term1) + abs(term2) + abs(term3) + abs(term4)
  except: 
    deltaReBeta = -1e0
  
  return deltaReBeta
#}}}

def deltaPeriodSPP(wavelength, eps1, eps2, deps1r, deps1c, deps2r, deps2c): #{{{
  ## returns the absolute uncertainty on the SPP period.
  ## this function was validated on one typical value where function is close to singularity. 
  deltaBeta = deltaBetaSPP(wavelength, eps1, eps2, deps1r, deps1c, deps2r, deps2c)
  beta = betaSPP(wavelength, eps1, eps2)
  beta = beta.real
  try: 
    value = 2e0*pi*deltaBeta / (beta**2)
  except:
    value = -1e0
  return value.real
#}}}

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

# Precision over knowledge of period
def deltaImBetaSPP(wavelength, eps1, eps2, deps1r, deps1c, deps2r, deps2c):#{{{
  """Calculates the precision over Im(beta) using uncertainty calculations
  """
  Er1 = eps1.real; Ec1 = eps1.imag
  Er2 = eps2.real; Ec2 = eps2.imag
  
  try:
    term1 = Pi*csgn(1e0,(Ec1*Er2**2+Er1**2*Ec2+Ec1**2*Ec2+Ec1*Ec2**2-1E0j*Er1**2*Er2-1E0j*Er1*Er2**2-1E0j*Ec1**2*Er2-1E0j*Er1*Ec2**2)*(Er1**2+2e0*Er1*Er2+Er2**2+Ec1**2+2e0*Ec1*Ec2+Ec2**2))*(2e0*(((Er1*Er2-Ec1*Ec2)*(Er1+Er2)/((Er1+Er2)**2+(Ec1+Ec2)**2)+(Ec1*Er2+Er1*Ec2)*(Ec1+Ec2)/((Er1+Er2)**2+(Ec1+Ec2)**2))**2+((Ec1*Er2+Er1*Ec2)*(Er1+Er2)/((Er1+Er2)**2+(Ec1+Ec2)**2)-(Er1*Er2-Ec1*Ec2)*(Ec1+Ec2)/((Er1+Er2)**2+(Ec1+Ec2)**2))**2)**(0.5e0)-(2e0*Er1*Er2-2e0*Ec1*Ec2)*(Er1+Er2)/((Er1+Er2)**2+(Ec1+Ec2)**2)-(2e0*Ec1*Er2+2e0*Er1*Ec2)*(Ec1+Ec2)/((Er1+Er2)**2+(Ec1+Ec2)**2))**(0.5e0)/wavelength+0.5e0*Pi*csgn((Ec1*Er2**2+Er1**2e0*Ec2+Ec1**2*Ec2+Ec1*Ec2**2-1E0j*Er1**2*Er2-1E0j*Er1*Er2**2-1E0j*Ec1**2*Er2-1E0j*Er1*Ec2**2)*(Er1**2+2e0*Er1*Er2+Er2**2+Ec1**2+2e0*Ec1*Ec2+Ec2**2))/(2e0*(((Er1*Er2-Ec1*Ec2)*(Er1+Er2)/((Er1+Er2)**2+(Ec1+Ec2)**2)+(Ec1*Er2+Er1*Ec2)*(Ec1+Ec2)/((Er1+Er2)**2+(Ec1+Ec2)**2))**2+((Ec1*Er2+Er1*Ec2)*(Er1+Er2)/((Er1+Er2)**2+(Ec1+Ec2)**2)-(Er1*Er2-Ec1*Ec2)*(Ec1+Ec2)/((Er1+Er2)**2+(Ec1+Ec2)**2))**2)**(0.5e0)-(2e0*Er1*Er2-2e0*Ec1*Ec2)*(Er1+Er2)/((Er1+Er2)**2+(Ec1+Ec2)**2)-(2e0*Ec1*Er2+2e0*Er1*Ec2)*(Ec1+Ec2)/((Er1+Er2)**2+(Ec1+Ec2)**2))**(0.5e0)/wavelength*(1E0/(((Er1*Er2-Ec1*Ec2)*(Er1+Er2)/((Er1+Er2)**2+(Ec1+Ec2)**2)+(Ec1*Er2+Er1*Ec2)*(Ec1+Ec2)/((Er1+Er2)**2+(Ec1+Ec2)**2))**2+((Ec1*Er2+Er1*Ec2)*(Er1+Er2)/((Er1+Er2)**2+(Ec1+Ec2)**2)-(Er1*Er2-Ec1*Ec2)*(Ec1+Ec2)/((Er1+Er2)**2+(Ec1+Ec2)**2))**2)**(0.5e0)*(2e0*((Er1*Er2-Ec1*Ec2)*(Er1+Er2)/((Er1+Er2)**2+(Ec1+Ec2)**2)+(Ec1*Er2+Er1*Ec2)*(Ec1+Ec2)/((Er1+Er2)**2+(Ec1+Ec2)**2))*(Er2*(Er1+Er2)/((Er1+Er2)**2+(Ec1+Ec2)**2)+(Er1*Er2-Ec1*Ec2)/((Er1+Er2)**2+(Ec1+Ec2)**2)-(Er1*Er2-Ec1*Ec2)*(Er1+Er2)/((Er1+Er2)**2+(Ec1+Ec2)**2)**2*(2e0*Er1+2e0*Er2)+Ec2*(Ec1+Ec2)/((Er1+Er2)**2+(Ec1+Ec2)**2)-(Ec1*Er2+Er1*Ec2)*(Ec1+Ec2)/((Er1+Er2)**2+(Ec1+Ec2)**2)**2*(2e0*Er1+2e0*Er2))+2e0*((Ec1*Er2+Er1*Ec2)*(Er1+Er2)/((Er1+Er2)**2+(Ec1+Ec2)**2)-(Er1*Er2-Ec1*Ec2)*(Ec1+Ec2)/((Er1+Er2)**2+(Ec1+Ec2)**2))*(Ec2*(Er1+Er2)/((Er1+Er2)**2+(Ec1+Ec2)**2)+(Ec1*Er2+Er1*Ec2)/((Er1+Er2)**2+(Ec1+Ec2)**2)-(Ec1*Er2+Er1*Ec2)*(Er1+Er2)/((Er1+Er2)**2+(Ec1+Ec2)**2)**2*(2e0*Er1+2e0*Er2)-Er2*(Ec1+Ec2)/((Er1+Er2)**2+(Ec1+Ec2)**2)+(Er1*Er2-Ec1*Ec2)*(Ec1+Ec2)/((Er1+Er2)**2+(Ec1+Ec2)**2)**2*(2e0*Er1+2e0*Er2)))-2e0*Er2*(Er1+Er2)/((Er1+Er2)**2+(Ec1+Ec2)**2)-(2e0*Er1*Er2-2e0*Ec1*Ec2)/((Er1+Er2)**2+(Ec1+Ec2)**2)+(2e0*Er1*Er2-2e0*Ec1*Ec2)*(Er1+Er2)/((Er1+Er2)**2+(Ec1+Ec2)**2)**2*(2e0*Er1+2e0*Er2)-2e0*Ec2*(Ec1+Ec2)/((Er1+Er2)**2+(Ec1+Ec2)**2)+(2e0*Ec1*Er2+2e0*Er1*Ec2)*(Ec1+Ec2)/((Er1+Er2)**2+(Ec1+Ec2)**2)**2*(2e0*Er1+2e0*Er2))
            
    term2 = Pi*csgn(1,(Ec1*Er2**2+Er1**2*Ec2+Ec1**2e0*Ec2+Ec1*Ec2**2-1E0j*Er1**2*Er2-1E0j*Er1*Er2**2-1E0j*Ec1**2*Er2-1E0j*Er1*Ec2**2)*(Er1**2+2e0*Er1*Er2+Er2**2+Ec1**2+2e0*Ec1*Ec2+Ec2**2))*(2e0*(((Er1*Er2-Ec1*Ec2)*(Er1+Er2)/((Er1+Er2)**2+(Ec1+Ec2)**2)+(Ec1*Er2+Er1*Ec2)*(Ec1+Ec2)/((Er1+Er2)**2+(Ec1+Ec2)**2))**2+((Ec1*Er2+Er1*Ec2)*(Er1+Er2)/((Er1+Er2)**2+(Ec1+Ec2)**2)-(Er1*Er2-Ec1*Ec2)*(Ec1+Ec2)/((Er1+Er2)**2+(Ec1+Ec2)**2))**2)**(0.5e0)-(2e0*Er1*Er2-2e0*Ec1*Ec2)*(Er1+Er2)/((Er1+Er2)**2+(Ec1+Ec2)**2)-(2e0*Ec1*Er2+2e0*Er1*Ec2)*(Ec1+Ec2)/((Er1+Er2)**2+(Ec1+Ec2)**2))**(0.5e0)/wavelength+0.5e0*Pi*csgn((Ec1*Er2**2+Er1**2*Ec2+Ec1**2*Ec2+Ec1*Ec2**2-1E0j*Er1**2*Er2-1E0j*Er1*Er2**2-1E0j*Ec1**2*Er2-1E0j*Er1*Ec2**2)*(Er1**2+2e0*Er1*Er2+Er2**2+Ec1**2+2e0*Ec1*Ec2+Ec2**2))/(2e0*(((Er1*Er2-Ec1*Ec2)*(Er1+Er2)/((Er1+Er2)**2+(Ec1+Ec2)**2)+(Ec1*Er2+Er1*Ec2)*(Ec1+Ec2)/((Er1+Er2)**2+(Ec1+Ec2)**2))**2+((Ec1*Er2+Er1*Ec2)*(Er1+Er2)/((Er1+Er2)**2+(Ec1+Ec2)**2)-(Er1*Er2-Ec1*Ec2)*(Ec1+Ec2)/((Er1+Er2)**2+(Ec1+Ec2)**2))**2)**(0.5e0)-(2e0*Er1*Er2-2e0*Ec1*Ec2)*(Er1+Er2)/((Er1+Er2)**2+(Ec1+Ec2)**2)-(2e0*Ec1*Er2+2e0*Er1*Ec2)*(Ec1+Ec2)/((Er1+Er2)**2+(Ec1+Ec2)**2))**(0.5e0)/wavelength*(1E0/(((Er1*Er2-Ec1*Ec2)*(Er1+Er2)/((Er1+Er2)**2+(Ec1+Ec2)**2)+(Ec1*Er2+Er1*Ec2)*(Ec1+Ec2)/((Er1+Er2)**2+(Ec1+Ec2)**2))**2+((Ec1*Er2+Er1*Ec2)*(Er1+Er2)/((Er1+Er2)**2+(Ec1+Ec2)**2)-(Er1*Er2-Ec1*Ec2)*(Ec1+Ec2)/((Er1+Er2)**2+(Ec1+Ec2)**2))**2)**(0.5e0)*(2e0*((Er1*Er2-Ec1*Ec2)*(Er1+Er2)/((Er1+Er2)**2+(Ec1+Ec2)**2)+(Ec1*Er2+Er1*Ec2)*(Ec1+Ec2)/((Er1+Er2)**2+(Ec1+Ec2)**2))*(-Ec2*(Er1+Er2)/((Er1+Er2)**2+(Ec1+Ec2)**2)-(Er1*Er2-Ec1*Ec2)*(Er1+Er2)/((Er1+Er2)**2+(Ec1+Ec2)**2)**2*(2e0*Ec1+2e0*Ec2)+Er2*(Ec1+Ec2)/((Er1+Er2)**2+(Ec1+Ec2)**2)+(Ec1*Er2+Er1*Ec2)/((Er1+Er2)**2+(Ec1+Ec2)**2)-(Ec1*Er2+Er1*Ec2)*(Ec1+Ec2)/((Er1+Er2)**2+(Ec1+Ec2)**2)**2*(2e0*Ec1+2e0*Ec2))+2e0*((Ec1*Er2+Er1*Ec2)*(Er1+Er2)/((Er1+Er2)**2+(Ec1+Ec2)**2)-(Er1*Er2-Ec1*Ec2)*(Ec1+Ec2)/((Er1+Er2)**2+(Ec1+Ec2)**2))*(Er2*(Er1+Er2)/((Er1+Er2)**2+(Ec1+Ec2)**2)-(Ec1*Er2+Er1*Ec2)*(Er1+Er2)/((Er1+Er2)**2+(Ec1+Ec2)**2)**2*(2e0*Ec1+2e0*Ec2)+Ec2*(Ec1+Ec2)/((Er1+Er2)**2+(Ec1+Ec2)**2)-(Er1*Er2-Ec1*Ec2)/((Er1+Er2)**2+(Ec1+Ec2)**2)+(Er1*Er2-Ec1*Ec2)*(Ec1+Ec2)/((Er1+Er2)**2+(Ec1+Ec2)**2)**2*(2e0*Ec1+2e0*Ec2)))+2e0*Ec2*(Er1+Er2)/((Er1+Er2)**2+(Ec1+Ec2)**2)+(2e0*Er1*Er2-2e0*Ec1*Ec2)*(Er1+Er2)/((Er1+Er2)**2+(Ec1+Ec2)**2)**2*(2e0*Ec1+2e0*Ec2)-2e0*Er2*(Ec1+Ec2)/((Er1+Er2)**2+(Ec1+Ec2)**2)-(2e0*Ec1*Er2+2e0*Er1*Ec2)/((Er1+Er2)**2+(Ec1+Ec2)**2)+(2e0*Ec1*Er2+2e0*Er1*Ec2)*(Ec1+Ec2)/((Er1+Er2)**2+(Ec1+Ec2)**2)**2*(2e0*Ec1+2e0*Ec2))
            
    term3 = Pi*csgn(1,(Ec1*Er2**2+Er1**2*Ec2+Ec1**2*Ec2+Ec1*Ec2**2-1E0j*Er1**2*Er2-1E0j*Er1*Er2**2-1E0j*Ec1**2*Er2-1E0j*Er1*Ec2**2)*(Er1**2+2e0*Er1*Er2+Er2**2+Ec1**2+2e0*Ec1*Ec2+Ec2**2))*(2e0*(((Er1*Er2-Ec1*Ec2)*(Er1+Er2)/((Er1+Er2)**2+(Ec1+Ec2)**2)+(Ec1*Er2+Er1*Ec2)*(Ec1+Ec2)/((Er1+Er2)**2+(Ec1+Ec2)**2))**2+((Ec1*Er2+Er1*Ec2)*(Er1+Er2)/((Er1+Er2)**2+(Ec1+Ec2)**2)-(Er1*Er2-Ec1*Ec2)*(Ec1+Ec2)/((Er1+Er2)**2+(Ec1+Ec2)**2))**2)**(0.5e0)-(2e0*Er1*Er2-2e0*Ec1*Ec2)*(Er1+Er2)/((Er1+Er2)**2+(Ec1+Ec2)**2)-(2e0*Ec1*Er2+2e0*Er1*Ec2)*(Ec1+Ec2)/((Er1+Er2)**2+(Ec1+Ec2)**2))**(0.5e0)/wavelength+0.5e0*Pi*csgn((Ec1*Er2**2+Er1**2*Ec2+Ec1**2*Ec2+Ec1*Ec2**2-1E0j*Er1**2*Er2-1E0j*Er1*Er2**2-1E0j*Ec1**2*Er2-1E0j*Er1*Ec2**2)*(Er1**2+2e0*Er1*Er2+Er2**2+Ec1**2+2e0*Ec1*Ec2+Ec2**2))/(2e0*(((Er1*Er2-Ec1*Ec2)*(Er1+Er2)/((Er1+Er2)**2+(Ec1+Ec2)**2)+(Ec1*Er2+Er1*Ec2)*(Ec1+Ec2)/((Er1+Er2)**2+(Ec1+Ec2)**2))**2+((Ec1*Er2+Er1*Ec2)*(Er1+Er2)/((Er1+Er2)**2+(Ec1+Ec2)**2)-(Er1*Er2-Ec1*Ec2)*(Ec1+Ec2)/((Er1+Er2)**2+(Ec1+Ec2)**2))**2)**(0.5e0)-(2e0*Er1*Er2-2e0*Ec1*Ec2)*(Er1+Er2)/((Er1+Er2)**2+(Ec1+Ec2)**2)-(2e0*Ec1*Er2+2e0*Er1*Ec2)*(Ec1+Ec2)/((Er1+Er2)**2+(Ec1+Ec2)**2))**(0.5e0)/wavelength*(1E0/(((Er1*Er2-Ec1*Ec2)*(Er1+Er2)/((Er1+Er2)**2+(Ec1+Ec2)**2)+(Ec1*Er2+Er1*Ec2)*(Ec1+Ec2)/((Er1+Er2)**2+(Ec1+Ec2)**2))**2+((Ec1*Er2+Er1*Ec2)*(Er1+Er2)/((Er1+Er2)**2+(Ec1+Ec2)**2)-(Er1*Er2-Ec1*Ec2)*(Ec1+Ec2)/((Er1+Er2)**2+(Ec1+Ec2)**2))**2)**(0.5e0)*(2e0*((Er1*Er2-Ec1*Ec2)*(Er1+Er2)/((Er1+Er2)**2+(Ec1+Ec2)**2)+(Ec1*Er2+Er1*Ec2)*(Ec1+Ec2)/((Er1+Er2)**2+(Ec1+Ec2)**2))*(Er1*(Er1+Er2)/((Er1+Er2)**2+(Ec1+Ec2)**2)+(Er1*Er2-Ec1*Ec2)/((Er1+Er2)**2+(Ec1+Ec2)**2)-(Er1*Er2-Ec1*Ec2)*(Er1+Er2)/((Er1+Er2)**2+(Ec1+Ec2)**2)**2*(2e0*Er1+2e0*Er2)+Ec1*(Ec1+Ec2)/((Er1+Er2)**2+(Ec1+Ec2)**2)-(Ec1*Er2+Er1*Ec2)*(Ec1+Ec2)/((Er1+Er2)**2+(Ec1+Ec2)**2)**2*(2e0*Er1+2e0*Er2))+2e0*((Ec1*Er2+Er1*Ec2)*(Er1+Er2)/((Er1+Er2)**2+(Ec1+Ec2)**2)-(Er1*Er2-Ec1*Ec2)*(Ec1+Ec2)/((Er1+Er2)**2+(Ec1+Ec2)**2))*(Ec1*(Er1+Er2)/((Er1+Er2)**2+(Ec1+Ec2)**2)+(Ec1*Er2+Er1*Ec2)/((Er1+Er2)**2+(Ec1+Ec2)**2)-(Ec1*Er2+Er1*Ec2)*(Er1+Er2)/((Er1+Er2)**2+(Ec1+Ec2)**2)**2*(2e0*Er1+2e0*Er2)-Er1*(Ec1+Ec2)/((Er1+Er2)**2+(Ec1+Ec2)**2)+(Er1*Er2-Ec1*Ec2)*(Ec1+Ec2)/((Er1+Er2)**2+(Ec1+Ec2)**2)**2*(2e0*Er1+2e0*Er2)))-2e0*Er1*(Er1+Er2)/((Er1+Er2)**2+(Ec1+Ec2)**2)-(2e0*Er1*Er2-2e0*Ec1*Ec2)/((Er1+Er2)**2+(Ec1+Ec2)**2)+(2e0*Er1*Er2-2e0*Ec1*Ec2)*(Er1+Er2)/((Er1+Er2)**2+(Ec1+Ec2)**2)**2*(2e0*Er1+2e0*Er2)-2e0*Ec1*(Ec1+Ec2)/((Er1+Er2)**2+(Ec1+Ec2)**2)+(2e0*Ec1*Er2+2e0*Er1*Ec2)*(Ec1+Ec2)/((Er1+Er2)**2+(Ec1+Ec2)**2)**2*(2e0*Er1+2e0*Er2))
            
    term4 = Pi*csgn(1,(Ec1*Er2**2+Er1**2*Ec2+Ec1**2*Ec2+Ec1*Ec2**2-1E0j*Er1**2*Er2-1E0j*Er1*Er2**2-1E0j*Ec1**2*Er2-1E0j*Er1*Ec2**2)*(Er1**2+2e0*Er1*Er2+Er2**2+Ec1**2+2e0*Ec1*Ec2+Ec2**2))*(2e0*(((Er1*Er2-Ec1*Ec2)*(Er1+Er2)/((Er1+Er2)**2+(Ec1+Ec2)**2)+(Ec1*Er2+Er1*Ec2)*(Ec1+Ec2)/((Er1+Er2)**2+(Ec1+Ec2)**2))**2+((Ec1*Er2+Er1*Ec2)*(Er1+Er2)/((Er1+Er2)**2+(Ec1+Ec2)**2)-(Er1*Er2-Ec1*Ec2)*(Ec1+Ec2)/((Er1+Er2)**2+(Ec1+Ec2)**2))**2)**(0.5e0)-(2e0*Er1*Er2-2e0*Ec1*Ec2)*(Er1+Er2)/((Er1+Er2)**2+(Ec1+Ec2)**2)-(2e0*Ec1*Er2+2e0*Er1*Ec2)*(Ec1+Ec2)/((Er1+Er2)**2+(Ec1+Ec2)**2))**(0.5e0)/wavelength+0.5e0*Pi*csgn((Ec1*Er2**2+Er1**2*Ec2+Ec1**2*Ec2+Ec1*Ec2**2-1E0j*Er1**2*Er2-1E0j*Er1*Er2**2-1E0j*Ec1**2*Er2-1E0j*Er1*Ec2**2)*(Er1**2+2e0*Er1*Er2+Er2**2+Ec1**2+2e0*Ec1*Ec2+Ec2**2))/(2e0*(((Er1*Er2-Ec1*Ec2)*(Er1+Er2)/((Er1+Er2)**2+(Ec1+Ec2)**2)+(Ec1*Er2+Er1*Ec2)*(Ec1+Ec2)/((Er1+Er2)**2+(Ec1+Ec2)**2))**2+((Ec1*Er2+Er1*Ec2)*(Er1+Er2)/((Er1+Er2)**2+(Ec1+Ec2)**2)-(Er1*Er2-Ec1*Ec2)*(Ec1+Ec2)/((Er1+Er2)**2+(Ec1+Ec2)**2))**2)**(0.5e0)-(2e0*Er1*Er2-2e0*Ec1*Ec2)*(Er1+Er2)/((Er1+Er2)**2+(Ec1+Ec2)**2)-(2e0*Ec1*Er2+2e0*Er1*Ec2)*(Ec1+Ec2)/((Er1+Er2)**2+(Ec1+Ec2)**2))**(0.5e0)/wavelength*(1E0/(((Er1*Er2-Ec1*Ec2)*(Er1+Er2)/((Er1+Er2)**2+(Ec1+Ec2)**2)+(Ec1*Er2+Er1*Ec2)*(Ec1+Ec2)/((Er1+Er2)**2+(Ec1+Ec2)**2))**2+((Ec1*Er2+Er1*Ec2)*(Er1+Er2)/((Er1+Er2)**2+(Ec1+Ec2)**2)-(Er1*Er2-Ec1*Ec2)*(Ec1+Ec2)/((Er1+Er2)**2+(Ec1+Ec2)**2))**2)**(0.5e0)*(2e0*((Er1*Er2-Ec1*Ec2)*(Er1+Er2)/((Er1+Er2)**2+(Ec1+Ec2)**2)+(Ec1*Er2+Er1*Ec2)*(Ec1+Ec2)/((Er1+Er2)**2+(Ec1+Ec2)**2))*(-Ec1*(Er1+Er2)/((Er1+Er2)**2+(Ec1+Ec2)**2)-(Er1*Er2-Ec1*Ec2)*(Er1+Er2)/((Er1+Er2)**2+(Ec1+Ec2)**2)**2*(2e0*Ec1+2e0*Ec2)+Er1*(Ec1+Ec2)/((Er1+Er2)**2+(Ec1+Ec2)**2)+(Ec1*Er2+Er1*Ec2)/((Er1+Er2)**2+(Ec1+Ec2)**2)-(Ec1*Er2+Er1*Ec2)*(Ec1+Ec2)/((Er1+Er2)**2+(Ec1+Ec2)**2)**2*(2e0*Ec1+2e0*Ec2))+2e0*((Ec1*Er2+Er1*Ec2)*(Er1+Er2)/((Er1+Er2)**2+(Ec1+Ec2)**2)-(Er1*Er2-Ec1*Ec2)*(Ec1+Ec2)/((Er1+Er2)**2+(Ec1+Ec2)**2))*(Er1*(Er1+Er2)/((Er1+Er2)**2+(Ec1+Ec2)**2)-(Ec1*Er2+Er1*Ec2)*(Er1+Er2)/((Er1+Er2)**2+(Ec1+Ec2)**2)**2*(2e0*Ec1+2e0*Ec2)+Ec1*(Ec1+Ec2)/((Er1+Er2)**2+(Ec1+Ec2)**2)-(Er1*Er2-Ec1*Ec2)/((Er1+Er2)**2+(Ec1+Ec2)**2)+(Er1*Er2-Ec1*Ec2)*(Ec1+Ec2)/((Er1+Er2)**2+(Ec1+Ec2)**2)**2*(2e0*Ec1+2e0*Ec2)))+2e0*Ec1*(Er1+Er2)/((Er1+Er2)**2+(Ec1+Ec2)**2)+(2e0*Er1*Er2-2e0*Ec1*Ec2)*(Er1+Er2)/((Er1+Er2)**2+(Ec1+Ec2)**2)**2*(2e0*Ec1+2e0*Ec2)-2e0*Er1*(Ec1+Ec2)/((Er1+Er2)**2+(Ec1+Ec2)**2)-(2e0*Ec1*Er2+2e0*Er1*Ec2)/((Er1+Er2)**2+(Ec1+Ec2)**2)+(2e0*Ec1*Er2+2e0*Er1*Ec2)*(Ec1+Ec2)/((Er1+Er2)**2+(Ec1+Ec2)**2)**2*(2e0*Ec1+2e0*Ec2))
            
    deltaImBeta = abs(term1*deps1r) + abs(term2*deps1c) + abs(term3*deps2r) + abs(term4*deps2c)
  except: 
    deltaImBeta = -1e0
  
  return deltaImBeta
#}}}

def deltaLspp(wavelength, eps1, eps2, deps1r, deps1c, deps2r, deps2c): #{{{
  ## returns the absolute uncertainty on the SPP mean free path.
  deltaBeta = deltaImBetaSPP(wavelength, eps1, eps2, deps1r, deps1c, deps2r, deps2c)
  beta = betaSPP(wavelength, eps1, eps2)
  beta = beta.imag
  try: 
    value = 0.5e0*abs(deltaBeta)/abs(beta)**2
  except:
    value = -1e0
  return value.real


### SPP decay depth
def DecayDepth(kzSPP):#{{{
  return 2e0*pi/kzSPP.real
#}}}

def kzSPP(wavelength,eps1,eps2):#{{{
  return cmath.sqrt(betaSPP(wavelength,eps1,eps2)**2-eps1*(omega(wavelength)**2/c**2))
#}}}

def DecayLengthSPP(beta):#{{{
    """Return the coherent length of SPPs
    """
    try:
      result = 1e0/(2e0*beta.imag)
    except: 
      result = -1
    return result
#}}}
  
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
  
## More elaborated functions

def ExperimentallyAchievable(OpticalPenetrationDepth, DecayDepth):#{{{
  """ Print only the experimentally possible cases: SPP active depth must be smaller than absorption depth. 
    # if Opd1 > 0, then: 
    #   return (Opd1 > SPPdecayDepth1)
    # else: 
    #   return true
  """
  if (OpticalPenetrationDepth != -1):
    #test if OPD > SPPdecayDepth
    return (OpticalPenetrationDepth > DecayDepth)
  else:
    return True
#}}}

def SPPactiveInterfaces(dbarray, comment):#{{{
  """Print all the SPP-active interfaces available in database
  CONSIDERS ONLY SYMMETRIC CASES
  If comment=="new", old SPP-active interfaces are removed from the table
  """
  counter=0
  #print len(dbarray)
  #sizeDatabase = len(dbarray)
  #sizeOfArray = sizeDatabase**2
  #print sizeOfArray	

  SPParray = np.empty((0,20)) #, dtype='|S30')
  
  # double loop to test all configurations (brute-forcing...)
  for i in dbarray:
    for k in dbarray:
      #extracting info on medium1 and medium2 in array; building eps1 and eps2
      name1=i[0]; name2=k[0]; wavelength1=1e-9*float(i[2]); wavelength2=1e-9*float(k[2])
      eps1=float(i[3])+1j*float(i[4]); eps2=float(k[3])+1j*float(k[4]) 
      try:
        gap1=float(i[1]); 
      except:
        gap1=10; 
      try: 
        gap2=float(k[1]); 
      except: 
        gap2=10;
      RegularLIPSScondition = False #(eps2.imag > 1e0*abs(eps1.real))  
      #ConditionOnGap=(gap1<0.1)
      ConditionOnGap=True #always true to avoid selection
      
      #calculate SPP condition
      if ((wavelength1 == wavelength2) and ConditionOnGap): #we must consider same wavelength, otherwise there is no meaning, but we also select only metallic substrates
        #if (gap2<0.1): #we select only metallic materials for interface 2 = substrate
        #print wavelength1, wavelength2
        #if (SPPcondition(eps1, eps2)): #the interface is SPP active
          
        # Build the table of SPP active interfaces
        Material1=name1
        Material2=name2
        Wavelength=(wavelength1/lengthunit)
        
        Absorption1 = (2e0*omega(wavelength1) / c) * (eps1)**0.5
        Absorption2 = (2e0*omega(wavelength2) / c) * (eps2)**0.5
        # Calculate optical penetration depth. -1 means infinite. 			
        if (Absorption1.imag == 0e0): 
	        OpticalPenetration1 = -1
        else: 
	        OpticalPenetration1 = 1e9 * 1e0/Absorption1.imag
	        
        if (Absorption2.imag == 0e0): 
	        OpticalPenetration2 = -1
        else:
	        OpticalPenetration2 = 1e9 * 1e0/Absorption2.imag
        # SPP activitivity condition with Perfect Medium Approximation ? 
        OldSPPactiveBool=''; 
        if (OldSPPcondition(eps1, eps2)): 
	        OldSPPactiveBool='Yes'
        else: 
	        OldSPPactiveBool='No'
        if (SPPcondition(eps1, eps2)): 
	        NewSPPactiveBool='Yes'
        else: 
	        NewSPPactiveBool='No'

        # If new or old SPP active condition is true, then show	
        if ((SPPcondition(eps1,eps2)) or (OldSPPcondition(eps1,eps2)) or RegularLIPSScondition):
	        Period=(period(betaSPP(wavelength1,eps1, eps2))/lengthunit)
	        SPPdecayDepth1=(DecayDepth(kzSPP(wavelength1, eps1, eps2))/lengthunit)
	        SPPdecayDepth2=(DecayDepth(kzSPP(wavelength2, eps2, eps1))/lengthunit)
	        SPPdecayLength=DecayLengthSPP(betaSPP(wavelength1,eps1, eps2))/lengthunit

	        SPPdepthImagk1 = kzSPP(wavelength1, eps1, eps2).imag
	        SPPdepthImagk2 = kzSPP(wavelength2, eps2, eps1).imag

	        PeriodError=deltaPeriodSPP(wavelength1, eps1, eps2, eta, eta, eta, eta)/lengthunit

        else: 
	        Period=0
	        SPPdecayDepth1=0
	        SPPdecayDepth2=0
	        SPPdecayLength=0
	        SPPdepthImagk1=0
	        SPPdepthImagk2=0
	        PeriodError=0
        
        Reflectivity=(reflectivity(eps1, eps2))
        deltaLsppValues=deltaLspp(wavelength1,eps1,eps2,eta,eta,eta,eta)/lengthunit
        
        # Print only the experimentally possible cases: SPP active depth must be smaller than absorption depth. 
        # ensure that SPPdecayDepth is smaller than layer thickness, to avoid shift of dispersion relation
        ExperimentalAchievable = True #ExperimentallyAchievable(OpticalPenetration1, SPPdecayDepth1)
        
        # Print the table of active SPP interfaces for all cases or only new SPP interfaces					
        #if( not (comment=="new")): 
        Condition = ExperimentalAchievable and (Period!=0) 
        #and (SPPdecayLength < 20000e0) and (abs(eps2.real) < eps2.imag)
        if(Condition):
          counter=counter+1
          #print SPParray.shape
          SPParray = np.vstack((SPParray, [Material1, Material2, Wavelength, OldSPPactiveBool, NewSPPactiveBool, 
	Period, PeriodError, SPPdecayDepth1, SPPdecayDepth2, Reflectivity, #10 
	OpticalPenetration1, OpticalPenetration2, SPPdecayLength, eps1.real, eps1.imag, #15
	eps2.real, eps2.imag, SPPdepthImagk1, SPPdepthImagk2, deltaLsppValues.real]))
          
  return SPParray
#}}}

#def AsymmetricSPPposActiveInterfaces(dbarray, comment):
  #"""Print all the SPP-active interfaces available in database
  #If comment=="new", old SPP-active interfaces are removed from the table
  #"""
  #counter=0
  
  ## double loop to test all configurations (brute-forcing...)
  #for i in dbarray:
    #for k in dbarray:
      ##extracting info on medium1 and medium2 in array; building eps1 and eps2
      #name1=i[0]; name2=k[0]; wavelength1=1e-9*float(i[2]); wavelength2=1e-9*float(k[2])
      #eps1=float(i[3])+1j*float(i[4]); eps2=float(k[3])+1j*float(k[4]) 
      #try:
        #gap1=float(i[1]); 
      #except:
        #gap1=10; 
      #try: 
        #gap2=float(k[1]); 
      #except: 
        #gap2=10;
        
      ##ConditionOnGap=(gap1<0.1)
      #ConditionOnGap=True #always true to avoid selection
      
      ##calculate SPP condition
      #if ((wavelength1 == wavelength2) and ConditionOnGap): #we must consider same wavelength, otherwise there is no meaning, but we also select only metallic substrates
        ##if (gap2<0.1): #we select only metallic materials for interface 2 = substrate
        ##print wavelength1, wavelength2
        ##if (SPPcondition(eps1, eps2)): #the interface is SPP active
          
        ## Build the table of SPP active interfaces
        #Material1=name1
        #Material2=name2
        #Wavelength=(wavelength1/lengthunit)
        
        #Absorption1 = (2e0*omega(wavelength1) / c) * (eps1)**0.5
        #Absorption2 = (2e0*omega(wavelength2) / c) * (eps2)**0.5
        ## Calculate optical penetration depth. -1 means infinite. 			
        #if (Absorption1.imag == 0e0): 
	        #OpticalPenetration1 = -1
        #else: 
	        #OpticalPenetration1 = 1e9 * 1e0/Absorption1.imag
	        
        #if (Absorption2.imag == 0e0): 
	        #OpticalPenetration2 = -1
        #else:
	        #OpticalPenetration2 = 1e9 * 1e0/Absorption2.imag
        ## SPP activitivity condition with Perfect Medium Approximation ? 
        #OldSPPactiveBool=''; 
        #if (OldSPPcondition(eps1, eps2)): 
	        #OldSPPactiveBool='Yes'
        #else: 
	        #OldSPPactiveBool='No'
        #if (AsymmetricSPPconditionPos(eps1, eps2)): 
	        #NewSPPactiveBool='Yes'
        #else: 
	        #NewSPPactiveBool='No'
        ## If new or old SPP active condition is true, then show	
        #if ((AsymmetricSPPconditionPos(eps1,eps2)) or (OldSPPcondition(eps1,eps2))):
	        #RealEps=(period(betaSPP(wavelength1,eps1, eps2))/lengthunit)
	        #SPPdecayDepth1=(DecayDepth(kzSPP(wavelength1, eps1, eps2))/lengthunit)
	        #SPPdecayDepth2=(DecayDepth(kzSPP(wavelength2, eps2, eps1))/lengthunit)
        #else: 
	        #RealEps=0
	        #SPPdecayDepth1=0
	        #SPPdecayDepth2=0
        
        #Reflectivity=(reflectivity(eps1, eps2))
        
        ## Print only the experimentally possible cases: SPP active depth must be smaller than absorption depth. 
        ## ensure that SPPdecayDepth is smaller than layer thickness, to avoid shift of dispersion relation
        #ExperimentalAchievable = ExperimentallyAchievable(OpticalPenetration1, SPPdecayDepth1)
        
        ## Print the table of active SPP interfaces for all cases or only new SPP interfaces					
        ##if( not (comment=="new")): 
        #if(ExperimentalAchievable and (RealEps!=0)):
          #if (np.mod(counter, 20) == 0): 
            ##show the table line each 20 lines
            #print '{0:30s} {1:30s} {2:15s} {3:12s} {4:12s} {5:11s} {6:16s} {7:16s} {8:12s} {9:19s} {10:19s}'.format("# Substrate", "Layer", "Wavelength (nm)", "OldSPPactive", "NewSPPactive", "Period (nm)", "DecayDepth1 (nm)", "DecayDepth2 (nm)", "Reflectivity", "OpticalPenetration1", "OpticalPenetration2")
          
          #counter=counter+1
          #print '{0:30s} {1:30s} {2:15f} {3:12s} {4:12s} {5:11f} {6:16f} {7:16f} {8:12f}  {9:19f} {10:19f}'.format(Material1, Material2, Wavelength, OldSPPactiveBool, NewSPPactiveBool, RealEps, SPPdecayDepth1, SPPdecayDepth2, Reflectivity, OpticalPenetration1, OpticalPenetration2)
        
  #return 0

#def AsymmetricSPPnegActiveInterfaces(dbarray, comment):
  #"""Print all the SPP-active interfaces available in database
  #If comment=="new", old SPP-active interfaces are removed from the table
  #"""
  #counter=0
  
  ## double loop to test all configurations (brute-forcing...)
  #for i in dbarray:
    #for k in dbarray:
      ##extracting info on medium1 and medium2 in array; building eps1 and eps2
      #name1=i[0]; name2=k[0]; wavelength1=1e-9*float(i[2]); wavelength2=1e-9*float(k[2])
      #eps1=float(i[3])+1j*float(i[4]); eps2=float(k[3])+1j*float(k[4]) 
      #try:
        #gap1=float(i[1]); 
      #except:
        #gap1=10; 
      #try: 
        #gap2=float(k[1]); 
      #except: 
        #gap2=10;
        
      ##ConditionOnGap=(gap1<0.1)
      #ConditionOnGap=True #always true to avoid selection
      
      ##calculate SPP condition
      #if ((wavelength1 == wavelength2) and ConditionOnGap): #we must consider same wavelength, otherwise there is no meaning, but we also select only metallic substrates
        ##if (gap2<0.1): #we select only metallic materials for interface 2 = substrate
        ##print wavelength1, wavelength2
        ##if (SPPcondition(eps1, eps2)): #the interface is SPP active
          
        ## Build the table of SPP active interfaces
        #Material1=name1
        #Material2=name2
        #Wavelength=(wavelength1/lengthunit)
        
        #Absorption1 = (2e0*omega(wavelength1) / c) * (eps1)**0.5
        #Absorption2 = (2e0*omega(wavelength2) / c) * (eps2)**0.5
        ## Calculate optical penetration depth. -1 means infinite. 			
        #if (Absorption1.imag == 0e0): 
	        #OpticalPenetration1 = -1
        #else: 
	        #OpticalPenetration1 = 1e9 * 1e0/Absorption1.imag
	        
        #if (Absorption2.imag == 0e0): 
	        #OpticalPenetration2 = -1
        #else:
	        #OpticalPenetration2 = 1e9 * 1e0/Absorption2.imag
        ## SPP activitivity condition with Perfect Medium Approximation ? 
        #OldSPPactiveBool=''; 
        #if (OldSPPcondition(eps1, eps2)): 
	        #OldSPPactiveBool='Yes'
        #else: 
	        #OldSPPactiveBool='No'
        #if (AsymmetricSPPconditionNeg(eps1, eps2)): 
	        #NewSPPactiveBool='Yes'
        #else: 
	        #NewSPPactiveBool='No'
        ## If new or old SPP active condition is true, then show	
        #if ((AsymmetricSPPconditionNeg(eps1,eps2)) or (OldSPPcondition(eps1,eps2))):
	        #RealEps=(period(betaSPP(wavelength1,eps1, eps2))/lengthunit)
	        #SPPdecayDepth1=(DecayDepth(kzSPP(wavelength1, eps1, eps2))/lengthunit)
	        #SPPdecayDepth2=(DecayDepth(kzSPP(wavelength2, eps2, eps1))/lengthunit)
        #else: 
	        #RealEps=0
	        #SPPdecayDepth1=0
	        #SPPdecayDepth2=0
        
        #Reflectivity=(reflectivity(eps1, eps2))
        
        ## Print only the experimentally possible cases: SPP active depth must be smaller than absorption depth. 
        ## ensure that SPPdecayDepth is smaller than layer thickness, to avoid shift of dispersion relation
        #ExperimentalAchievable = True #ExperimentallyAchievable(OpticalPenetration1, SPPdecayDepth1)
        
        ## Print the table of active SPP interfaces for all cases or only new SPP interfaces					
        ##if( not (comment=="new")): 
        #if(ExperimentalAchievable and (RealEps!=0)):
          #if (np.mod(counter, 20) == 0): 
            ##show the table line each 20 lines
            #print '{0:30s} {1:30s} {2:15s} {3:12s} {4:12s} {5:11s} {6:16s} {7:16s} {8:12s} {9:19s} {10:19s}'.format("# Substrate", "Layer", "Wavelength (nm)", "OldSPPactive", "NewSPPactive", "Period (nm)", "DecayDepth1 (nm)", "DecayDepth2 (nm)", "Reflectivity", "OpticalPenetration1", "OpticalPenetration2")
          
          #counter=counter+1
          #print '{0:30s} {1:30s} {2:15f} {3:12s} {4:12s} {5:11f} {6:16f} {7:16f} {8:12f}  {9:19f} {10:19f}'.format(Material1, Material2, Wavelength, OldSPPactiveBool, NewSPPactiveBool, RealEps, SPPdecayDepth1, SPPdecayDepth2, Reflectivity, OpticalPenetration1, OpticalPenetration2)
        
  #return 0

def GenerateDatabase():
  """
  We would like now to construct a database using available materials description with all possible interfaces
  we will : 

  1. For each material in database, select each material and verify, for each available wavelength, 
  1.1: If SPP condition is verified, 
  1.2. yes, then period can be calculated and shown;
  1.4. SPP decay depth in medium 1
  1.5. SPP decay depth in medium 2
  2. Then extract a table which contains all possible scenarios
  """
  # Select database
  database="MaterialOpticalDatabaseForPlasmonics.csv"

  # Build database array for choosing which material can be of interest to irradiate
  dbarray = loadtxt(database, dtype='str', delimiter='\t')

  # To calculate symmetric SPP compatible interfaces, use the following line
  SPPactiveInterfacesArray = SPPactiveInterfaces(dbarray, '')

  # To calculate asymmetric POSITIVE SPP compatible interfaces, use the following line
  #SPPactiveInterfacesArray = AsymmetricSPPposActiveInterfaces(dbarray, '')

  # To calculate asymmetric NEGATIVE SPP compatible interfaces, use the following line
  #SPPactiveInterfacesArray = AsymmetricSPPnegActiveInterfaces(dbarray, '')
  return SPPactiveInterfacesArray


def ExportToTxt(dbarray, filename):
  """
  Export an SPP array to a CSV file
  SPP array must be produced with one of the SPPactiveInterfaces functions
  """
  try: 
    np.savetxt(filename, dbarray, fmt="%s", delimiter='\t', newline='\n',comments='#')
    out = 0
  except: 
    print "Could not output SPP database into a file"
    out = 1
  
  #counter=0
  
  #for i in dbarray:
    #for k in dbarray:
      #if (np.mod(counter, 20) == 0):
	#show the table line each 20 lines, but also put it in a table
	#if (comment):
	#print '{0:30s} {1:30s} {2:15s} {3:12s} {4:12s} {5:11s} {6:16s} {7:16s} {8:12s} {9:19s} {10:19s} {11:15s}'.format("# Substrate", "Layer", "Wavelength (nm)", "OldSPPactive", "NewSPPactive", "Period (nm)", "DecayDepth1 (nm)", "DecayDepth2 (nm)", "Reflectivity", "OpticalPenetration1", "OpticalPenetration2", "DecayLength")
      #print dbarray[counter,:]
      #Material1 = i
      #counter=counter+1
  #print Material1
      #print '{0:30s} {1:30s} {2:15f} {3:12s} {4:12s} {5:11f} {6:16f} {7:16f} {8:12f}  {9:19f} {10:19f} {11:15f}'.format(Material1, Material2, Wavelength, OldSPPactiveBool, NewSPPactiveBool, RealEps, SPPdecayDepth1, SPPdecayDepth2, Reflectivity, OpticalPenetration1, OpticalPenetration2, SPPdecayLength)
  return out

def RealDerivativeByComplex(f,z):
  """Complex derivative a real-valued function by a complex-number
  Input:
    f: z->f(z)
    z: z complex-valued numbers
  Output:
    df/dz according to Wirtinger complex derivatives formula
  """
  return 0.5e0 * (np.diff(f)/np.diff(z.real) - 1j*(np.diff(f)/np.diff(z.imag))) #original

def LifeTimeRaether(beta, eps2, eps1):
	omegasppimag=beta.real * c * eps1.imag/(2e0*eps1.real**2) * (eps1.real * eps2.real)/(eps1.real + eps2.real)
	# lifetime=2e0*pi/omegasppimag #Raether formula
	lifetime=1e0/omegasppimag #modified Raether formula to match with complex group velocity approach
	return lifetime

def SPPlength(beta): #{{{
  length = 1e0/(2e0 * beta.imag) #Maier formula
  return length
#}}}

def LifeTimeDerrien(beta, vg):
  length = SPPlength(beta)
  lifetime = length * (vg)**(-1e0)
  return lifetime
              
def EpsilonToIndex(eps):
  #returns the complex refractive index
  return cmath.sqrt(eps)

def IndexToEpsilon(n):
  #returns the complex permittivity from optical index
  return n*n

def EffectiveIndex(eps1, eps2): 
  #returns effective optical index of SPP
  #neff = 
  return cmath.sqrt( eps1*eps2 / (eps1+eps2) )


def FilterDatabase(SPPdb, query, FieldIndex):
  """ Filter SPP database using query and returns a smaller database
  /!\ content of query cell should be exact
  """
  SPPdbFiltered = np.array(SPPdb[SPPdb[:,FieldIndex]==query,:]) #uses a table of booleans to select
  return SPPdbFiltered

def FilterDatabaseContains(SPPdb, query, FieldIndex):
  """ Filter SPP database using query and returns a smaller database
  """
  SPPdbFiltered = SPPdb[np.array(np.core.defchararray.find(SPPdb[:,FieldIndex], query)==0),:]
  return SPPdbFiltered

def FilterDatabaseLowerThan(SPPdb, query, FieldIndex):
  """ Filter SPP database using query and returns a smaller database
  /!\ content of query cell should be exact
  """
  SPPdbFiltered = np.array(SPPdb[SPPdb[:,FieldIndex]<query,:]) #uses a table of booleans to select
  return SPPdbFiltered

def ExtractMaterialData(Database): #TODO: Think how to take data from continuous database directly instead of the ponctual file.
  # Extract data from database of materials (5 columns)
  Material1 = Database[:, 0]; BandGap = Database[:,1]; 
  Wavelength = Database[:, 2]; 
  RealEps = Database[:,3]; ImagEps = Database[:,4]; 
  
  #Converts strings to floats
  #Material1 = np.asfarray(Material1)
  #BandGap = np.asfarray(BandGap)
  #Wavelength = np.asfarray(Wavelength)
  #RealEps = np.asfarray(RealEps)
  #ImagEps = np.asfarray(ImagEps)
  
  return Material1, BandGap, Wavelength, RealEps, ImagEps

def ExtractDataDb(SPPdbFiltered):
  # Extract data from database
  Material1 = SPPdbFiltered[:, 0]; Material2 = SPPdbFiltered[:,1]; 
  Wavelength = SPPdbFiltered[:, 2]; 
  OldSPPactiveBool = SPPdbFiltered[:,3]; NewSPPactiveBool = SPPdbFiltered[:,4]; 
  RealEps = SPPdbFiltered[:,5]; RealEpsError = SPPdbFiltered[:,6];
  SPPdecayDepth1 = SPPdbFiltered[:,7]; SPPdecayDepth2 = SPPdbFiltered[:,8]; 
  Reflectivity = SPPdbFiltered[:, 9]; OpticalPenetration1 = SPPdbFiltered[:,10]; OpticalPenetration2 = SPPdbFiltered[:,10]; SPPdecayLength = SPPdbFiltered[:,12]
  eps1r = SPPdbFiltered[:, 13]; eps1c = SPPdbFiltered[:,14]; eps2r = SPPdbFiltered[:,15]; eps2c = SPPdbFiltered[:,16]; k1imag = SPPdbFiltered[:,17]; 
  k2imag = SPPdbFiltered[:,18]
  DeltaLsppValue = SPPdbFiltered[:,19]
  
  #Converts strings to floats
  Wavelength = np.asfarray(Wavelength)
  RealEps = np.asfarray(RealEps)
  RealEpsError = np.asfarray(RealEpsError)
  SPPdecayDepth1 = np.asfarray(SPPdecayDepth1)
  SPPdecayDepth2 = np.asfarray(SPPdecayDepth2)
  Reflectivity = np.asfarray(Reflectivity)
  OpticalPenetration1 = np.asfarray(OpticalPenetration1)
  OpticalPenetration2 = np.asfarray(OpticalPenetration2)
  SPPdecayLength = np.asfarray(SPPdecayLength)
  eps1r = np.asfarray(eps1r)
  eps1c = np.asfarray(eps1c)
  eps2r = np.asfarray(eps2r)
  eps2c = np.asfarray(eps2c)
  k1imag = np.asfarray(k1imag)
  k2imag = np.asfarray(k2imag)
  DeltaLsppValue = np.asfarray(DeltaLsppValue)

  return Material1, Material2, Wavelength, OldSPPactiveBool, NewSPPactiveBool, RealEps, RealEpsError, SPPdecayDepth1, SPPdecayDepth2, Reflectivity, OpticalPenetration1, OpticalPenetration2, SPPdecayLength, eps1r, eps1c, eps2r, eps2c, k1imag, k2imag, DeltaLsppValue

def Swap(eps1, eps2):#{{{
  eps3 = eps1
  eps1 = eps2
  eps2 = eps3
  del eps3
  return(eps1, eps2) 
#}}}

kzSPP = np.vectorize(kzSPP)
DecayDepth = np.vectorize(DecayDepth)
omega = np.vectorize(omega) 