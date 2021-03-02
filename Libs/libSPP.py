#!/usr/bin/env python2
#-*- coding: utf-8 -*-

# Copyright (C) 2013-2020 T. J.-Y. Derrien
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

## @package libSPP
# Module libSPP explores the SPP theory at a single interface between 
# two semi-infinite media. The formal model is presented in 
# T.J.-Y. Derrien et al, Journal of Optics 18, 115007 (2016)

# IMPORT PYTHON LIBRARIES
import numpy as np
from numpy import genfromtxt, loadtxt, chararray
from scipy.optimize import fsolve, root
import cmath
import matplotlib as mp
import matplotlib.pyplot as plt
from scipy.interpolate import InterpolatedUnivariateSpline
from matplotlib import rc, font_manager
# from pylab import *
from scipy.constants import c, epsilon_0, mu_0, pi, e, m_e, h
from matplotlib.legend_handler import HandlerLine2D
import sys
from colorama import Fore
from colorama import Style

# IMPORT CUSTOM LIBRARIES
# from libKeldysh import *
from libDatabase import *
from libLaser import *
from libMaterials import *
from libMath import *

lengthunit = 1e-9
eta = 0.1 #assumed precision error on the dielectric permittivity
UsingTeX=True #TODO: set to False for Windows or linux cluster users

## 0: all permisive, no verification on SPP excitation condition
## 1: use the RegularLIPSScondition, softer than pure SPP excitation condition
## 2: Period != 0 is necessary for a material to be listed in results
## 3: Extreme level: use ExperimentallyAchievable() to verify possibility of decay depth > optical penetration depth
LevelOfSPPaccuracy=0

# Settings for matplotlib: taken from https://stackoverflow.com/questions/12322738/how-do-i-change-the-axis-tick-font-in-a-matplotlib-plot-when-rendering-using-lat
sizeOfFont = 14
FontName='cm' #Helvetica
fontProperties = {'family':'sans-serif','sans-serif':[FontName],
    'weight' : 'normal', 'size' : sizeOfFont}
ticks_font = font_manager.FontProperties(family=FontName, style='italic',
    size=sizeOfFont, weight='normal', stretch='normal')	
rc('font',**fontProperties)
#rc('text.latex', preamble=r'\usepackage{cmbright}')
#rc('font',**{'family':'sans-serif','sans-serif':['Arial'], 'size':'18'})
## for Palatino and other serif fonts use:
#rc('font', **{'family':'serif', 'serif':['Palatino'], 'size':'18'})
rc('text', usetex=UsingTeX)
mp.rcParams['legend.numpoints'] = 1

## basic wave function
#def omega(wavelength):#{{{
    #return 2.0*pi*c/wavelength
##}}}

# SPP BASIC FUNCTIONS

## Computes the SPP wave number on a flat interface
# @param wavelength (float), 
# @param eps1 (complex), 
# @param eps2 (complex)
def betaSPP(wavelength, eps1, eps2):#{{{

    omega = 2.0*pi*c/wavelength
    try:
        value = omega/c * cmath.sqrt(eps1 * eps2 / (eps1 + eps2))
    except: 
        print("**Info: betaSPP: singular case")
        value = -1e0+0e0j
    return value
#}}}

## OBSELETE. Assume that Re(k1).Re(k2) > 0 and verify the subsequent consequences altogether.
#  It exists then two sub-modes, let's say a positive one (Im k1<0, Im k2>0), 
#  and a negative one (Im k1>0, Im k2<0).
#  This routine is about: Im(k1)>0, Im(k2)<0
def AsymmetricSPPconditionPos(eps1, eps2):#{{{
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

## SPPconditionValue() returns the value of condition for SPP. If its negative, then SPP can be excited at a flat interface. 
#  /!\ This condition is restricted to checking the real part of the dispersion relation for symmetric SPP only. 
#  Input: 
#  @param eps1: complex-valued permittivity of medium 1
#  @param eps2: complex-valued permittivity of medium 2
#  Output: float
def SPPconditionValue(eps1, eps2):#{{{
  condition=eps1.real*eps2.real+eps1.imag*eps2.imag
  return condition
#}}}

## Returns a boolean stating if SPP are excitable at the interface
# defined by eps1 | eps2 joint media. 
def SPPcondition(eps1, eps2):#{{{
  if (SPPconditionValue(eps1, eps2) < 0.0):
           output=True
  else:
    output=False
  return output
#}}}

## SPPconditionValue() returns the value of condition for SPP IN PERFECT MATERIALS (Im(eps)<<|Re(eps)). If its negative, then SPP can be excited at a flat interface. 
#  Input: eps1, eps2: complex-valued quantities
#  Output: float
def OldSPPcondition(eps1, eps2):#{{{
  condition1=(eps1.real*eps2.real<0.0)
  #condition2= eps2.real < abs(eps1.real) #this version is not symmetric, hence strange
  # Let's use its generalization which is actually symmetric. 
  condition2 = (eps1.real * eps2.real / (eps1.real + eps2.real) > 0e0)
  return (condition1 and condition2)
#}}}

## Returns the period of the light-SPP field at a given interface
def period(betaSPP):#{{{
  try:
    result = 2.0*pi/betaSPP.real
  except: 
    result = -1e0
  return result
#}}}

## Precision over knowledge of the Re(complex-valued SPP wavenumber)
def deltaBetaSPP(wavelength, eps1, eps2, deps1r, deps1c, deps2r, deps2c):#{{{
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

## returns the absolute uncertainty on the SPP period.
#  this function was validated on one typical value where function is close to singularity.  
def deltaPeriodSPP(wavelength, eps1, eps2, deps1r, deps1c, deps2r, deps2c): #{{{
  deltaBeta = deltaBetaSPP(wavelength, eps1, eps2, deps1r, deps1c, deps2r, deps2c)
  beta = betaSPP(wavelength, eps1, eps2)
  beta = beta.real
  try: 
    value = 2e0*pi*deltaBeta / (beta**2)
  except:
    value = -1e0
  return value.real
#}}}

# Precision over knowledge of Im(complex SPP wavenumber)
def deltaImBetaSPP(wavelength, eps1, eps2, deps1r, deps1c, deps2r, deps2c):#{{{
  Er1 = eps1.real; Ec1 = eps1.imag
  Er2 = eps2.real; Ec2 = eps2.imag
  
  try:
  #term1 = pi*csgn(1e0,(Ec1*Er2**2+Er1**2*Ec2+Ec1**2*Ec2+Ec1*Ec2**2-1E0j*Er1**2*Er2-1E0j*Er1*Er2**2-1E0j*Ec1**2*Er2-1E0j*Er1*Ec2**2)*(Er1**2+2e0*Er1*Er2+Er2**2+Ec1**2+2e0*Ec1*Ec2+Ec2**2))*(2e0*(((Er1*Er2-Ec1*Ec2)*(Er1+Er2)/((Er1+Er2)**2+(Ec1+Ec2)**2)+(Ec1*Er2+Er1*Ec2)*(Ec1+Ec2)/((Er1+Er2)**2+(Ec1+Ec2)**2))**2+((Ec1*Er2+Er1*Ec2)*(Er1+Er2)/((Er1+Er2)**2+(Ec1+Ec2)**2)-(Er1*Er2-Ec1*Ec2)*(Ec1+Ec2)/((Er1+Er2)**2+(Ec1+Ec2)**2))**2)**(0.5e0)-(2e0*Er1*Er2-2e0*Ec1*Ec2)*(Er1+Er2)/((Er1+Er2)**2+(Ec1+Ec2)**2)-(2e0*Ec1*Er2+2e0*Er1*Ec2)*(Ec1+Ec2)/((Er1+Er2)**2+(Ec1+Ec2)**2))**(0.5e0)/wavelength+0.5e0*pi*csgn((Ec1*Er2**2+Er1**2e0*Ec2+Ec1**2*Ec2+Ec1*Ec2**2-1E0j*Er1**2*Er2-1E0j*Er1*Er2**2-1E0j*Ec1**2*Er2-1E0j*Er1*Ec2**2)*(Er1**2+2e0*Er1*Er2+Er2**2+Ec1**2+2e0*Ec1*Ec2+Ec2**2))/(2e0*(((Er1*Er2-Ec1*Ec2)*(Er1+Er2)/((Er1+Er2)**2+(Ec1+Ec2)**2)+(Ec1*Er2+Er1*Ec2)*(Ec1+Ec2)/((Er1+Er2)**2+(Ec1+Ec2)**2))**2+((Ec1*Er2+Er1*Ec2)*(Er1+Er2)/((Er1+Er2)**2+(Ec1+Ec2)**2)-(Er1*Er2-Ec1*Ec2)*(Ec1+Ec2)/((Er1+Er2)**2+(Ec1+Ec2)**2))**2)**(0.5e0)-(2e0*Er1*Er2-2e0*Ec1*Ec2)*(Er1+Er2)/((Er1+Er2)**2+(Ec1+Ec2)**2)-(2e0*Ec1*Er2+2e0*Er1*Ec2)*(Ec1+Ec2)/((Er1+Er2)**2+(Ec1+Ec2)**2))**(0.5e0)/wavelength*(1E0/(((Er1*Er2-Ec1*Ec2)*(Er1+Er2)/((Er1+Er2)**2+(Ec1+Ec2)**2)+(Ec1*Er2+Er1*Ec2)*(Ec1+Ec2)/((Er1+Er2)**2+(Ec1+Ec2)**2))**2+((Ec1*Er2+Er1*Ec2)*(Er1+Er2)/((Er1+Er2)**2+(Ec1+Ec2)**2)-(Er1*Er2-Ec1*Ec2)*(Ec1+Ec2)/((Er1+Er2)**2+(Ec1+Ec2)**2))**2)**(0.5e0)*(2e0*((Er1*Er2-Ec1*Ec2)*(Er1+Er2)/((Er1+Er2)**2+(Ec1+Ec2)**2)+(Ec1*Er2+Er1*Ec2)*(Ec1+Ec2)/((Er1+Er2)**2+(Ec1+Ec2)**2))*(Er2*(Er1+Er2)/((Er1+Er2)**2+(Ec1+Ec2)**2)+(Er1*Er2-Ec1*Ec2)/((Er1+Er2)**2+(Ec1+Ec2)**2)-(Er1*Er2-Ec1*Ec2)*(Er1+Er2)/((Er1+Er2)**2+(Ec1+Ec2)**2)**2*(2e0*Er1+2e0*Er2)+Ec2*(Ec1+Ec2)/((Er1+Er2)**2+(Ec1+Ec2)**2)-(Ec1*Er2+Er1*Ec2)*(Ec1+Ec2)/((Er1+Er2)**2+(Ec1+Ec2)**2)**2*(2e0*Er1+2e0*Er2))+2e0*((Ec1*Er2+Er1*Ec2)*(Er1+Er2)/((Er1+Er2)**2+(Ec1+Ec2)**2)-(Er1*Er2-Ec1*Ec2)*(Ec1+Ec2)/((Er1+Er2)**2+(Ec1+Ec2)**2))*(Ec2*(Er1+Er2)/((Er1+Er2)**2+(Ec1+Ec2)**2)+(Ec1*Er2+Er1*Ec2)/((Er1+Er2)**2+(Ec1+Ec2)**2)-(Ec1*Er2+Er1*Ec2)*(Er1+Er2)/((Er1+Er2)**2+(Ec1+Ec2)**2)**2*(2e0*Er1+2e0*Er2)-Er2*(Ec1+Ec2)/((Er1+Er2)**2+(Ec1+Ec2)**2)+(Er1*Er2-Ec1*Ec2)*(Ec1+Ec2)/((Er1+Er2)**2+(Ec1+Ec2)**2)**2*(2e0*Er1+2e0*Er2)))-2e0*Er2*(Er1+Er2)/((Er1+Er2)**2+(Ec1+Ec2)**2)-(2e0*Er1*Er2-2e0*Ec1*Ec2)/((Er1+Er2)**2+(Ec1+Ec2)**2)+(2e0*Er1*Er2-2e0*Ec1*Ec2)*(Er1+Er2)/((Er1+Er2)**2+(Ec1+Ec2)**2)**2*(2e0*Er1+2e0*Er2)-2e0*Ec2*(Ec1+Ec2)/((Er1+Er2)**2+(Ec1+Ec2)**2)+(2e0*Ec1*Er2+2e0*Er1*Ec2)*(Ec1+Ec2)/((Er1+Er2)**2+(Ec1+Ec2)**2)**2*(2e0*Er1+2e0*Er2))
    term1 = pi*(2e0*(((Er1*Er2-Ec1*Ec2)*(Er1+Er2)/((Er1+Er2)**2+(Ec1+Ec2)**2)+(Ec1*Er2+Er1*Ec2)*(Ec1+Ec2)/((Er1+Er2)**2+(Ec1+Ec2)**2))**2+((Ec1*Er2+Er1*Ec2)*(Er1+Er2)/((Er1+Er2)**2+(Ec1+Ec2)**2)-(Er1*Er2-Ec1*Ec2)*(Ec1+Ec2)/((Er1+Er2)**2+(Ec1+Ec2)**2))**2)**(0.5e0)-(2e0*Er1*Er2-2e0*Ec1*Ec2)*(Er1+Er2)/((Er1+Er2)**2+(Ec1+Ec2)**2)-(2e0*Ec1*Er2+2e0*Er1*Ec2)*(Ec1+Ec2)/((Er1+Er2)**2+(Ec1+Ec2)**2))**(0.5e0)/wavelength+0.5e0*pi/(2e0*(((Er1*Er2-Ec1*Ec2)*(Er1+Er2)/((Er1+Er2)**2+(Ec1+Ec2)**2)+(Ec1*Er2+Er1*Ec2)*(Ec1+Ec2)/((Er1+Er2)**2+(Ec1+Ec2)**2))**2+((Ec1*Er2+Er1*Ec2)*(Er1+Er2)/((Er1+Er2)**2+(Ec1+Ec2)**2)-(Er1*Er2-Ec1*Ec2)*(Ec1+Ec2)/((Er1+Er2)**2+(Ec1+Ec2)**2))**2)**(0.5e0)-(2e0*Er1*Er2-2e0*Ec1*Ec2)*(Er1+Er2)/((Er1+Er2)**2+(Ec1+Ec2)**2)-(2e0*Ec1*Er2+2e0*Er1*Ec2)*(Ec1+Ec2)/((Er1+Er2)**2+(Ec1+Ec2)**2))**(0.5e0)/wavelength*(1E0/(((Er1*Er2-Ec1*Ec2)*(Er1+Er2)/((Er1+Er2)**2+(Ec1+Ec2)**2)+(Ec1*Er2+Er1*Ec2)*(Ec1+Ec2)/((Er1+Er2)**2+(Ec1+Ec2)**2))**2+((Ec1*Er2+Er1*Ec2)*(Er1+Er2)/((Er1+Er2)**2+(Ec1+Ec2)**2)-(Er1*Er2-Ec1*Ec2)*(Ec1+Ec2)/((Er1+Er2)**2+(Ec1+Ec2)**2))**2)**(0.5e0)*(2e0*((Er1*Er2-Ec1*Ec2)*(Er1+Er2)/((Er1+Er2)**2+(Ec1+Ec2)**2)+(Ec1*Er2+Er1*Ec2)*(Ec1+Ec2)/((Er1+Er2)**2+(Ec1+Ec2)**2))*(Er2*(Er1+Er2)/((Er1+Er2)**2+(Ec1+Ec2)**2)+(Er1*Er2-Ec1*Ec2)/((Er1+Er2)**2+(Ec1+Ec2)**2)-(Er1*Er2-Ec1*Ec2)*(Er1+Er2)/((Er1+Er2)**2+(Ec1+Ec2)**2)**2*(2e0*Er1+2e0*Er2)+Ec2*(Ec1+Ec2)/((Er1+Er2)**2+(Ec1+Ec2)**2)-(Ec1*Er2+Er1*Ec2)*(Ec1+Ec2)/((Er1+Er2)**2+(Ec1+Ec2)**2)**2*(2e0*Er1+2e0*Er2))+2e0*((Ec1*Er2+Er1*Ec2)*(Er1+Er2)/((Er1+Er2)**2+(Ec1+Ec2)**2)-(Er1*Er2-Ec1*Ec2)*(Ec1+Ec2)/((Er1+Er2)**2+(Ec1+Ec2)**2))*(Ec2*(Er1+Er2)/((Er1+Er2)**2+(Ec1+Ec2)**2)+(Ec1*Er2+Er1*Ec2)/((Er1+Er2)**2+(Ec1+Ec2)**2)-(Ec1*Er2+Er1*Ec2)*(Er1+Er2)/((Er1+Er2)**2+(Ec1+Ec2)**2)**2*(2e0*Er1+2e0*Er2)-Er2*(Ec1+Ec2)/((Er1+Er2)**2+(Ec1+Ec2)**2)+(Er1*Er2-Ec1*Ec2)*(Ec1+Ec2)/((Er1+Er2)**2+(Ec1+Ec2)**2)**2*(2e0*Er1+2e0*Er2)))-2e0*Er2*(Er1+Er2)/((Er1+Er2)**2+(Ec1+Ec2)**2)-(2e0*Er1*Er2-2e0*Ec1*Ec2)/((Er1+Er2)**2+(Ec1+Ec2)**2)+(2e0*Er1*Er2-2e0*Ec1*Ec2)*(Er1+Er2)/((Er1+Er2)**2+(Ec1+Ec2)**2)**2*(2e0*Er1+2e0*Er2)-2e0*Ec2*(Ec1+Ec2)/((Er1+Er2)**2+(Ec1+Ec2)**2)+(2e0*Ec1*Er2+2e0*Er1*Ec2)*(Ec1+Ec2)/((Er1+Er2)**2+(Ec1+Ec2)**2)**2*(2e0*Er1+2e0*Er2))
          
  #term2 = pi*csgn(1,(Ec1*Er2**2+Er1**2*Ec2+Ec1**2e0*Ec2+Ec1*Ec2**2-1E0j*Er1**2*Er2-1E0j*Er1*Er2**2-1E0j*Ec1**2*Er2-1E0j*Er1*Ec2**2)*(Er1**2+2e0*Er1*Er2+Er2**2+Ec1**2+2e0*Ec1*Ec2+Ec2**2))*(2e0*(((Er1*Er2-Ec1*Ec2)*(Er1+Er2)/((Er1+Er2)**2+(Ec1+Ec2)**2)+(Ec1*Er2+Er1*Ec2)*(Ec1+Ec2)/((Er1+Er2)**2+(Ec1+Ec2)**2))**2+((Ec1*Er2+Er1*Ec2)*(Er1+Er2)/((Er1+Er2)**2+(Ec1+Ec2)**2)-(Er1*Er2-Ec1*Ec2)*(Ec1+Ec2)/((Er1+Er2)**2+(Ec1+Ec2)**2))**2)**(0.5e0)-(2e0*Er1*Er2-2e0*Ec1*Ec2)*(Er1+Er2)/((Er1+Er2)**2+(Ec1+Ec2)**2)-(2e0*Ec1*Er2+2e0*Er1*Ec2)*(Ec1+Ec2)/((Er1+Er2)**2+(Ec1+Ec2)**2))**(0.5e0)/wavelength+0.5e0*pi*csgn((Ec1*Er2**2+Er1**2*Ec2+Ec1**2*Ec2+Ec1*Ec2**2-1E0j*Er1**2*Er2-1E0j*Er1*Er2**2-1E0j*Ec1**2*Er2-1E0j*Er1*Ec2**2)*(Er1**2+2e0*Er1*Er2+Er2**2+Ec1**2+2e0*Ec1*Ec2+Ec2**2))/(2e0*(((Er1*Er2-Ec1*Ec2)*(Er1+Er2)/((Er1+Er2)**2+(Ec1+Ec2)**2)+(Ec1*Er2+Er1*Ec2)*(Ec1+Ec2)/((Er1+Er2)**2+(Ec1+Ec2)**2))**2+((Ec1*Er2+Er1*Ec2)*(Er1+Er2)/((Er1+Er2)**2+(Ec1+Ec2)**2)-(Er1*Er2-Ec1*Ec2)*(Ec1+Ec2)/((Er1+Er2)**2+(Ec1+Ec2)**2))**2)**(0.5e0)-(2e0*Er1*Er2-2e0*Ec1*Ec2)*(Er1+Er2)/((Er1+Er2)**2+(Ec1+Ec2)**2)-(2e0*Ec1*Er2+2e0*Er1*Ec2)*(Ec1+Ec2)/((Er1+Er2)**2+(Ec1+Ec2)**2))**(0.5e0)/wavelength*(1E0/(((Er1*Er2-Ec1*Ec2)*(Er1+Er2)/((Er1+Er2)**2+(Ec1+Ec2)**2)+(Ec1*Er2+Er1*Ec2)*(Ec1+Ec2)/((Er1+Er2)**2+(Ec1+Ec2)**2))**2+((Ec1*Er2+Er1*Ec2)*(Er1+Er2)/((Er1+Er2)**2+(Ec1+Ec2)**2)-(Er1*Er2-Ec1*Ec2)*(Ec1+Ec2)/((Er1+Er2)**2+(Ec1+Ec2)**2))**2)**(0.5e0)*(2e0*((Er1*Er2-Ec1*Ec2)*(Er1+Er2)/((Er1+Er2)**2+(Ec1+Ec2)**2)+(Ec1*Er2+Er1*Ec2)*(Ec1+Ec2)/((Er1+Er2)**2+(Ec1+Ec2)**2))*(-Ec2*(Er1+Er2)/((Er1+Er2)**2+(Ec1+Ec2)**2)-(Er1*Er2-Ec1*Ec2)*(Er1+Er2)/((Er1+Er2)**2+(Ec1+Ec2)**2)**2*(2e0*Ec1+2e0*Ec2)+Er2*(Ec1+Ec2)/((Er1+Er2)**2+(Ec1+Ec2)**2)+(Ec1*Er2+Er1*Ec2)/((Er1+Er2)**2+(Ec1+Ec2)**2)-(Ec1*Er2+Er1*Ec2)*(Ec1+Ec2)/((Er1+Er2)**2+(Ec1+Ec2)**2)**2*(2e0*Ec1+2e0*Ec2))+2e0*((Ec1*Er2+Er1*Ec2)*(Er1+Er2)/((Er1+Er2)**2+(Ec1+Ec2)**2)-(Er1*Er2-Ec1*Ec2)*(Ec1+Ec2)/((Er1+Er2)**2+(Ec1+Ec2)**2))*(Er2*(Er1+Er2)/((Er1+Er2)**2+(Ec1+Ec2)**2)-(Ec1*Er2+Er1*Ec2)*(Er1+Er2)/((Er1+Er2)**2+(Ec1+Ec2)**2)**2*(2e0*Ec1+2e0*Ec2)+Ec2*(Ec1+Ec2)/((Er1+Er2)**2+(Ec1+Ec2)**2)-(Er1*Er2-Ec1*Ec2)/((Er1+Er2)**2+(Ec1+Ec2)**2)+(Er1*Er2-Ec1*Ec2)*(Ec1+Ec2)/((Er1+Er2)**2+(Ec1+Ec2)**2)**2*(2e0*Ec1+2e0*Ec2)))+2e0*Ec2*(Er1+Er2)/((Er1+Er2)**2+(Ec1+Ec2)**2)+(2e0*Er1*Er2-2e0*Ec1*Ec2)*(Er1+Er2)/((Er1+Er2)**2+(Ec1+Ec2)**2)**2*(2e0*Ec1+2e0*Ec2)-2e0*Er2*(Ec1+Ec2)/((Er1+Er2)**2+(Ec1+Ec2)**2)-(2e0*Ec1*Er2+2e0*Er1*Ec2)/((Er1+Er2)**2+(Ec1+Ec2)**2)+(2e0*Ec1*Er2+2e0*Er1*Ec2)*(Ec1+Ec2)/((Er1+Er2)**2+(Ec1+Ec2)**2)**2*(2e0*Ec1+2e0*Ec2))
    term2 = pi*(2e0*(((Er1*Er2-Ec1*Ec2)*(Er1+Er2)/((Er1+Er2)**2+(Ec1+Ec2)**2)+(Ec1*Er2+Er1*Ec2)*(Ec1+Ec2)/((Er1+Er2)**2+(Ec1+Ec2)**2))**2+((Ec1*Er2+Er1*Ec2)*(Er1+Er2)/((Er1+Er2)**2+(Ec1+Ec2)**2)-(Er1*Er2-Ec1*Ec2)*(Ec1+Ec2)/((Er1+Er2)**2+(Ec1+Ec2)**2))**2)**(0.5e0)-(2e0*Er1*Er2-2e0*Ec1*Ec2)*(Er1+Er2)/((Er1+Er2)**2+(Ec1+Ec2)**2)-(2e0*Ec1*Er2+2e0*Er1*Ec2)*(Ec1+Ec2)/((Er1+Er2)**2+(Ec1+Ec2)**2))**(0.5e0)/wavelength+0.5e0*pi/(2e0*(((Er1*Er2-Ec1*Ec2)*(Er1+Er2)/((Er1+Er2)**2+(Ec1+Ec2)**2)+(Ec1*Er2+Er1*Ec2)*(Ec1+Ec2)/((Er1+Er2)**2+(Ec1+Ec2)**2))**2+((Ec1*Er2+Er1*Ec2)*(Er1+Er2)/((Er1+Er2)**2+(Ec1+Ec2)**2)-(Er1*Er2-Ec1*Ec2)*(Ec1+Ec2)/((Er1+Er2)**2+(Ec1+Ec2)**2))**2)**(0.5e0)-(2e0*Er1*Er2-2e0*Ec1*Ec2)*(Er1+Er2)/((Er1+Er2)**2+(Ec1+Ec2)**2)-(2e0*Ec1*Er2+2e0*Er1*Ec2)*(Ec1+Ec2)/((Er1+Er2)**2+(Ec1+Ec2)**2))**(0.5e0)/wavelength*(1E0/(((Er1*Er2-Ec1*Ec2)*(Er1+Er2)/((Er1+Er2)**2+(Ec1+Ec2)**2)+(Ec1*Er2+Er1*Ec2)*(Ec1+Ec2)/((Er1+Er2)**2+(Ec1+Ec2)**2))**2+((Ec1*Er2+Er1*Ec2)*(Er1+Er2)/((Er1+Er2)**2+(Ec1+Ec2)**2)-(Er1*Er2-Ec1*Ec2)*(Ec1+Ec2)/((Er1+Er2)**2+(Ec1+Ec2)**2))**2)**(0.5e0)*(2e0*((Er1*Er2-Ec1*Ec2)*(Er1+Er2)/((Er1+Er2)**2+(Ec1+Ec2)**2)+(Ec1*Er2+Er1*Ec2)*(Ec1+Ec2)/((Er1+Er2)**2+(Ec1+Ec2)**2))*(-Ec2*(Er1+Er2)/((Er1+Er2)**2+(Ec1+Ec2)**2)-(Er1*Er2-Ec1*Ec2)*(Er1+Er2)/((Er1+Er2)**2+(Ec1+Ec2)**2)**2*(2e0*Ec1+2e0*Ec2)+Er2*(Ec1+Ec2)/((Er1+Er2)**2+(Ec1+Ec2)**2)+(Ec1*Er2+Er1*Ec2)/((Er1+Er2)**2+(Ec1+Ec2)**2)-(Ec1*Er2+Er1*Ec2)*(Ec1+Ec2)/((Er1+Er2)**2+(Ec1+Ec2)**2)**2*(2e0*Ec1+2e0*Ec2))+2e0*((Ec1*Er2+Er1*Ec2)*(Er1+Er2)/((Er1+Er2)**2+(Ec1+Ec2)**2)-(Er1*Er2-Ec1*Ec2)*(Ec1+Ec2)/((Er1+Er2)**2+(Ec1+Ec2)**2))*(Er2*(Er1+Er2)/((Er1+Er2)**2+(Ec1+Ec2)**2)-(Ec1*Er2+Er1*Ec2)*(Er1+Er2)/((Er1+Er2)**2+(Ec1+Ec2)**2)**2*(2e0*Ec1+2e0*Ec2)+Ec2*(Ec1+Ec2)/((Er1+Er2)**2+(Ec1+Ec2)**2)-(Er1*Er2-Ec1*Ec2)/((Er1+Er2)**2+(Ec1+Ec2)**2)+(Er1*Er2-Ec1*Ec2)*(Ec1+Ec2)/((Er1+Er2)**2+(Ec1+Ec2)**2)**2*(2e0*Ec1+2e0*Ec2)))+2e0*Ec2*(Er1+Er2)/((Er1+Er2)**2+(Ec1+Ec2)**2)+(2e0*Er1*Er2-2e0*Ec1*Ec2)*(Er1+Er2)/((Er1+Er2)**2+(Ec1+Ec2)**2)**2*(2e0*Ec1+2e0*Ec2)-2e0*Er2*(Ec1+Ec2)/((Er1+Er2)**2+(Ec1+Ec2)**2)-(2e0*Ec1*Er2+2e0*Er1*Ec2)/((Er1+Er2)**2+(Ec1+Ec2)**2)+(2e0*Ec1*Er2+2e0*Er1*Ec2)*(Ec1+Ec2)/((Er1+Er2)**2+(Ec1+Ec2)**2)**2*(2e0*Ec1+2e0*Ec2))
  
  #term3 = pi*csgn(1,(Ec1*Er2**2+Er1**2*Ec2+Ec1**2*Ec2+Ec1*Ec2**2-1E0j*Er1**2*Er2-1E0j*Er1*Er2**2-1E0j*Ec1**2*Er2-1E0j*Er1*Ec2**2)*(Er1**2+2e0*Er1*Er2+Er2**2+Ec1**2+2e0*Ec1*Ec2+Ec2**2))*(2e0*(((Er1*Er2-Ec1*Ec2)*(Er1+Er2)/((Er1+Er2)**2+(Ec1+Ec2)**2)+(Ec1*Er2+Er1*Ec2)*(Ec1+Ec2)/((Er1+Er2)**2+(Ec1+Ec2)**2))**2+((Ec1*Er2+Er1*Ec2)*(Er1+Er2)/((Er1+Er2)**2+(Ec1+Ec2)**2)-(Er1*Er2-Ec1*Ec2)*(Ec1+Ec2)/((Er1+Er2)**2+(Ec1+Ec2)**2))**2)**(0.5e0)-(2e0*Er1*Er2-2e0*Ec1*Ec2)*(Er1+Er2)/((Er1+Er2)**2+(Ec1+Ec2)**2)-(2e0*Ec1*Er2+2e0*Er1*Ec2)*(Ec1+Ec2)/((Er1+Er2)**2+(Ec1+Ec2)**2))**(0.5e0)/wavelength+0.5e0*pi*csgn((Ec1*Er2**2+Er1**2*Ec2+Ec1**2*Ec2+Ec1*Ec2**2-1E0j*Er1**2*Er2-1E0j*Er1*Er2**2-1E0j*Ec1**2*Er2-1E0j*Er1*Ec2**2)*(Er1**2+2e0*Er1*Er2+Er2**2+Ec1**2+2e0*Ec1*Ec2+Ec2**2))/(2e0*(((Er1*Er2-Ec1*Ec2)*(Er1+Er2)/((Er1+Er2)**2+(Ec1+Ec2)**2)+(Ec1*Er2+Er1*Ec2)*(Ec1+Ec2)/((Er1+Er2)**2+(Ec1+Ec2)**2))**2+((Ec1*Er2+Er1*Ec2)*(Er1+Er2)/((Er1+Er2)**2+(Ec1+Ec2)**2)-(Er1*Er2-Ec1*Ec2)*(Ec1+Ec2)/((Er1+Er2)**2+(Ec1+Ec2)**2))**2)**(0.5e0)-(2e0*Er1*Er2-2e0*Ec1*Ec2)*(Er1+Er2)/((Er1+Er2)**2+(Ec1+Ec2)**2)-(2e0*Ec1*Er2+2e0*Er1*Ec2)*(Ec1+Ec2)/((Er1+Er2)**2+(Ec1+Ec2)**2))**(0.5e0)/wavelength*(1E0/(((Er1*Er2-Ec1*Ec2)*(Er1+Er2)/((Er1+Er2)**2+(Ec1+Ec2)**2)+(Ec1*Er2+Er1*Ec2)*(Ec1+Ec2)/((Er1+Er2)**2+(Ec1+Ec2)**2))**2+((Ec1*Er2+Er1*Ec2)*(Er1+Er2)/((Er1+Er2)**2+(Ec1+Ec2)**2)-(Er1*Er2-Ec1*Ec2)*(Ec1+Ec2)/((Er1+Er2)**2+(Ec1+Ec2)**2))**2)**(0.5e0)*(2e0*((Er1*Er2-Ec1*Ec2)*(Er1+Er2)/((Er1+Er2)**2+(Ec1+Ec2)**2)+(Ec1*Er2+Er1*Ec2)*(Ec1+Ec2)/((Er1+Er2)**2+(Ec1+Ec2)**2))*(Er1*(Er1+Er2)/((Er1+Er2)**2+(Ec1+Ec2)**2)+(Er1*Er2-Ec1*Ec2)/((Er1+Er2)**2+(Ec1+Ec2)**2)-(Er1*Er2-Ec1*Ec2)*(Er1+Er2)/((Er1+Er2)**2+(Ec1+Ec2)**2)**2*(2e0*Er1+2e0*Er2)+Ec1*(Ec1+Ec2)/((Er1+Er2)**2+(Ec1+Ec2)**2)-(Ec1*Er2+Er1*Ec2)*(Ec1+Ec2)/((Er1+Er2)**2+(Ec1+Ec2)**2)**2*(2e0*Er1+2e0*Er2))+2e0*((Ec1*Er2+Er1*Ec2)*(Er1+Er2)/((Er1+Er2)**2+(Ec1+Ec2)**2)-(Er1*Er2-Ec1*Ec2)*(Ec1+Ec2)/((Er1+Er2)**2+(Ec1+Ec2)**2))*(Ec1*(Er1+Er2)/((Er1+Er2)**2+(Ec1+Ec2)**2)+(Ec1*Er2+Er1*Ec2)/((Er1+Er2)**2+(Ec1+Ec2)**2)-(Ec1*Er2+Er1*Ec2)*(Er1+Er2)/((Er1+Er2)**2+(Ec1+Ec2)**2)**2*(2e0*Er1+2e0*Er2)-Er1*(Ec1+Ec2)/((Er1+Er2)**2+(Ec1+Ec2)**2)+(Er1*Er2-Ec1*Ec2)*(Ec1+Ec2)/((Er1+Er2)**2+(Ec1+Ec2)**2)**2*(2e0*Er1+2e0*Er2)))-2e0*Er1*(Er1+Er2)/((Er1+Er2)**2+(Ec1+Ec2)**2)-(2e0*Er1*Er2-2e0*Ec1*Ec2)/((Er1+Er2)**2+(Ec1+Ec2)**2)+(2e0*Er1*Er2-2e0*Ec1*Ec2)*(Er1+Er2)/((Er1+Er2)**2+(Ec1+Ec2)**2)**2*(2e0*Er1+2e0*Er2)-2e0*Ec1*(Ec1+Ec2)/((Er1+Er2)**2+(Ec1+Ec2)**2)+(2e0*Ec1*Er2+2e0*Er1*Ec2)*(Ec1+Ec2)/((Er1+Er2)**2+(Ec1+Ec2)**2)**2*(2e0*Er1+2e0*Er2))
    term3 = pi*(2e0*(((Er1*Er2-Ec1*Ec2)*(Er1+Er2)/((Er1+Er2)**2+(Ec1+Ec2)**2)+(Ec1*Er2+Er1*Ec2)*(Ec1+Ec2)/((Er1+Er2)**2+(Ec1+Ec2)**2))**2+((Ec1*Er2+Er1*Ec2)*(Er1+Er2)/((Er1+Er2)**2+(Ec1+Ec2)**2)-(Er1*Er2-Ec1*Ec2)*(Ec1+Ec2)/((Er1+Er2)**2+(Ec1+Ec2)**2))**2)**(0.5e0)-(2e0*Er1*Er2-2e0*Ec1*Ec2)*(Er1+Er2)/((Er1+Er2)**2+(Ec1+Ec2)**2)-(2e0*Ec1*Er2+2e0*Er1*Ec2)*(Ec1+Ec2)/((Er1+Er2)**2+(Ec1+Ec2)**2))**(0.5e0)/wavelength+0.5e0*pi/(2e0*(((Er1*Er2-Ec1*Ec2)*(Er1+Er2)/((Er1+Er2)**2+(Ec1+Ec2)**2)+(Ec1*Er2+Er1*Ec2)*(Ec1+Ec2)/((Er1+Er2)**2+(Ec1+Ec2)**2))**2+((Ec1*Er2+Er1*Ec2)*(Er1+Er2)/((Er1+Er2)**2+(Ec1+Ec2)**2)-(Er1*Er2-Ec1*Ec2)*(Ec1+Ec2)/((Er1+Er2)**2+(Ec1+Ec2)**2))**2)**(0.5e0)-(2e0*Er1*Er2-2e0*Ec1*Ec2)*(Er1+Er2)/((Er1+Er2)**2+(Ec1+Ec2)**2)-(2e0*Ec1*Er2+2e0*Er1*Ec2)*(Ec1+Ec2)/((Er1+Er2)**2+(Ec1+Ec2)**2))**(0.5e0)/wavelength*(1E0/(((Er1*Er2-Ec1*Ec2)*(Er1+Er2)/((Er1+Er2)**2+(Ec1+Ec2)**2)+(Ec1*Er2+Er1*Ec2)*(Ec1+Ec2)/((Er1+Er2)**2+(Ec1+Ec2)**2))**2+((Ec1*Er2+Er1*Ec2)*(Er1+Er2)/((Er1+Er2)**2+(Ec1+Ec2)**2)-(Er1*Er2-Ec1*Ec2)*(Ec1+Ec2)/((Er1+Er2)**2+(Ec1+Ec2)**2))**2)**(0.5e0)*(2e0*((Er1*Er2-Ec1*Ec2)*(Er1+Er2)/((Er1+Er2)**2+(Ec1+Ec2)**2)+(Ec1*Er2+Er1*Ec2)*(Ec1+Ec2)/((Er1+Er2)**2+(Ec1+Ec2)**2))*(Er1*(Er1+Er2)/((Er1+Er2)**2+(Ec1+Ec2)**2)+(Er1*Er2-Ec1*Ec2)/((Er1+Er2)**2+(Ec1+Ec2)**2)-(Er1*Er2-Ec1*Ec2)*(Er1+Er2)/((Er1+Er2)**2+(Ec1+Ec2)**2)**2*(2e0*Er1+2e0*Er2)+Ec1*(Ec1+Ec2)/((Er1+Er2)**2+(Ec1+Ec2)**2)-(Ec1*Er2+Er1*Ec2)*(Ec1+Ec2)/((Er1+Er2)**2+(Ec1+Ec2)**2)**2*(2e0*Er1+2e0*Er2))+2e0*((Ec1*Er2+Er1*Ec2)*(Er1+Er2)/((Er1+Er2)**2+(Ec1+Ec2)**2)-(Er1*Er2-Ec1*Ec2)*(Ec1+Ec2)/((Er1+Er2)**2+(Ec1+Ec2)**2))*(Ec1*(Er1+Er2)/((Er1+Er2)**2+(Ec1+Ec2)**2)+(Ec1*Er2+Er1*Ec2)/((Er1+Er2)**2+(Ec1+Ec2)**2)-(Ec1*Er2+Er1*Ec2)*(Er1+Er2)/((Er1+Er2)**2+(Ec1+Ec2)**2)**2*(2e0*Er1+2e0*Er2)-Er1*(Ec1+Ec2)/((Er1+Er2)**2+(Ec1+Ec2)**2)+(Er1*Er2-Ec1*Ec2)*(Ec1+Ec2)/((Er1+Er2)**2+(Ec1+Ec2)**2)**2*(2e0*Er1+2e0*Er2)))-2e0*Er1*(Er1+Er2)/((Er1+Er2)**2+(Ec1+Ec2)**2)-(2e0*Er1*Er2-2e0*Ec1*Ec2)/((Er1+Er2)**2+(Ec1+Ec2)**2)+(2e0*Er1*Er2-2e0*Ec1*Ec2)*(Er1+Er2)/((Er1+Er2)**2+(Ec1+Ec2)**2)**2*(2e0*Er1+2e0*Er2)-2e0*Ec1*(Ec1+Ec2)/((Er1+Er2)**2+(Ec1+Ec2)**2)+(2e0*Ec1*Er2+2e0*Er1*Ec2)*(Ec1+Ec2)/((Er1+Er2)**2+(Ec1+Ec2)**2)**2*(2e0*Er1+2e0*Er2))
          
  #term4 = pi*csgn(1,(Ec1*Er2**2+Er1**2*Ec2+Ec1**2*Ec2+Ec1*Ec2**2-1E0j*Er1**2*Er2-1E0j*Er1*Er2**2-1E0j*Ec1**2*Er2-1E0j*Er1*Ec2**2)*(Er1**2+2e0*Er1*Er2+Er2**2+Ec1**2+2e0*Ec1*Ec2+Ec2**2))*(2e0*(((Er1*Er2-Ec1*Ec2)*(Er1+Er2)/((Er1+Er2)**2+(Ec1+Ec2)**2)+(Ec1*Er2+Er1*Ec2)*(Ec1+Ec2)/((Er1+Er2)**2+(Ec1+Ec2)**2))**2+((Ec1*Er2+Er1*Ec2)*(Er1+Er2)/((Er1+Er2)**2+(Ec1+Ec2)**2)-(Er1*Er2-Ec1*Ec2)*(Ec1+Ec2)/((Er1+Er2)**2+(Ec1+Ec2)**2))**2)**(0.5e0)-(2e0*Er1*Er2-2e0*Ec1*Ec2)*(Er1+Er2)/((Er1+Er2)**2+(Ec1+Ec2)**2)-(2e0*Ec1*Er2+2e0*Er1*Ec2)*(Ec1+Ec2)/((Er1+Er2)**2+(Ec1+Ec2)**2))**(0.5e0)/wavelength+0.5e0*pi*csgn((Ec1*Er2**2+Er1**2*Ec2+Ec1**2*Ec2+Ec1*Ec2**2-1E0j*Er1**2*Er2-1E0j*Er1*Er2**2-1E0j*Ec1**2*Er2-1E0j*Er1*Ec2**2)*(Er1**2+2e0*Er1*Er2+Er2**2+Ec1**2+2e0*Ec1*Ec2+Ec2**2))/(2e0*(((Er1*Er2-Ec1*Ec2)*(Er1+Er2)/((Er1+Er2)**2+(Ec1+Ec2)**2)+(Ec1*Er2+Er1*Ec2)*(Ec1+Ec2)/((Er1+Er2)**2+(Ec1+Ec2)**2))**2+((Ec1*Er2+Er1*Ec2)*(Er1+Er2)/((Er1+Er2)**2+(Ec1+Ec2)**2)-(Er1*Er2-Ec1*Ec2)*(Ec1+Ec2)/((Er1+Er2)**2+(Ec1+Ec2)**2))**2)**(0.5e0)-(2e0*Er1*Er2-2e0*Ec1*Ec2)*(Er1+Er2)/((Er1+Er2)**2+(Ec1+Ec2)**2)-(2e0*Ec1*Er2+2e0*Er1*Ec2)*(Ec1+Ec2)/((Er1+Er2)**2+(Ec1+Ec2)**2))**(0.5e0)/wavelength*(1E0/(((Er1*Er2-Ec1*Ec2)*(Er1+Er2)/((Er1+Er2)**2+(Ec1+Ec2)**2)+(Ec1*Er2+Er1*Ec2)*(Ec1+Ec2)/((Er1+Er2)**2+(Ec1+Ec2)**2))**2+((Ec1*Er2+Er1*Ec2)*(Er1+Er2)/((Er1+Er2)**2+(Ec1+Ec2)**2)-(Er1*Er2-Ec1*Ec2)*(Ec1+Ec2)/((Er1+Er2)**2+(Ec1+Ec2)**2))**2)**(0.5e0)*(2e0*((Er1*Er2-Ec1*Ec2)*(Er1+Er2)/((Er1+Er2)**2+(Ec1+Ec2)**2)+(Ec1*Er2+Er1*Ec2)*(Ec1+Ec2)/((Er1+Er2)**2+(Ec1+Ec2)**2))*(-Ec1*(Er1+Er2)/((Er1+Er2)**2+(Ec1+Ec2)**2)-(Er1*Er2-Ec1*Ec2)*(Er1+Er2)/((Er1+Er2)**2+(Ec1+Ec2)**2)**2*(2e0*Ec1+2e0*Ec2)+Er1*(Ec1+Ec2)/((Er1+Er2)**2+(Ec1+Ec2)**2)+(Ec1*Er2+Er1*Ec2)/((Er1+Er2)**2+(Ec1+Ec2)**2)-(Ec1*Er2+Er1*Ec2)*(Ec1+Ec2)/((Er1+Er2)**2+(Ec1+Ec2)**2)**2*(2e0*Ec1+2e0*Ec2))+2e0*((Ec1*Er2+Er1*Ec2)*(Er1+Er2)/((Er1+Er2)**2+(Ec1+Ec2)**2)-(Er1*Er2-Ec1*Ec2)*(Ec1+Ec2)/((Er1+Er2)**2+(Ec1+Ec2)**2))*(Er1*(Er1+Er2)/((Er1+Er2)**2+(Ec1+Ec2)**2)-(Ec1*Er2+Er1*Ec2)*(Er1+Er2)/((Er1+Er2)**2+(Ec1+Ec2)**2)**2*(2e0*Ec1+2e0*Ec2)+Ec1*(Ec1+Ec2)/((Er1+Er2)**2+(Ec1+Ec2)**2)-(Er1*Er2-Ec1*Ec2)/((Er1+Er2)**2+(Ec1+Ec2)**2)+(Er1*Er2-Ec1*Ec2)*(Ec1+Ec2)/((Er1+Er2)**2+(Ec1+Ec2)**2)**2*(2e0*Ec1+2e0*Ec2)))+2e0*Ec1*(Er1+Er2)/((Er1+Er2)**2+(Ec1+Ec2)**2)+(2e0*Er1*Er2-2e0*Ec1*Ec2)*(Er1+Er2)/((Er1+Er2)**2+(Ec1+Ec2)**2)**2*(2e0*Ec1+2e0*Ec2)-2e0*Er1*(Ec1+Ec2)/((Er1+Er2)**2+(Ec1+Ec2)**2)-(2e0*Ec1*Er2+2e0*Er1*Ec2)/((Er1+Er2)**2+(Ec1+Ec2)**2)+(2e0*Ec1*Er2+2e0*Er1*Ec2)*(Ec1+Ec2)/((Er1+Er2)**2+(Ec1+Ec2)**2)**2*(2e0*Ec1+2e0*Ec2))
    term4 = pi*(2e0*(((Er1*Er2-Ec1*Ec2)*(Er1+Er2)/((Er1+Er2)**2+(Ec1+Ec2)**2)+(Ec1*Er2+Er1*Ec2)*(Ec1+Ec2)/((Er1+Er2)**2+(Ec1+Ec2)**2))**2+((Ec1*Er2+Er1*Ec2)*(Er1+Er2)/((Er1+Er2)**2+(Ec1+Ec2)**2)-(Er1*Er2-Ec1*Ec2)*(Ec1+Ec2)/((Er1+Er2)**2+(Ec1+Ec2)**2))**2)**(0.5e0)-(2e0*Er1*Er2-2e0*Ec1*Ec2)*(Er1+Er2)/((Er1+Er2)**2+(Ec1+Ec2)**2)-(2e0*Ec1*Er2+2e0*Er1*Ec2)*(Ec1+Ec2)/((Er1+Er2)**2+(Ec1+Ec2)**2))**(0.5e0)/wavelength+0.5e0*pi/(2e0*(((Er1*Er2-Ec1*Ec2)*(Er1+Er2)/((Er1+Er2)**2+(Ec1+Ec2)**2)+(Ec1*Er2+Er1*Ec2)*(Ec1+Ec2)/((Er1+Er2)**2+(Ec1+Ec2)**2))**2+((Ec1*Er2+Er1*Ec2)*(Er1+Er2)/((Er1+Er2)**2+(Ec1+Ec2)**2)-(Er1*Er2-Ec1*Ec2)*(Ec1+Ec2)/((Er1+Er2)**2+(Ec1+Ec2)**2))**2)**(0.5e0)-(2e0*Er1*Er2-2e0*Ec1*Ec2)*(Er1+Er2)/((Er1+Er2)**2+(Ec1+Ec2)**2)-(2e0*Ec1*Er2+2e0*Er1*Ec2)*(Ec1+Ec2)/((Er1+Er2)**2+(Ec1+Ec2)**2))**(0.5e0)/wavelength*(1E0/(((Er1*Er2-Ec1*Ec2)*(Er1+Er2)/((Er1+Er2)**2+(Ec1+Ec2)**2)+(Ec1*Er2+Er1*Ec2)*(Ec1+Ec2)/((Er1+Er2)**2+(Ec1+Ec2)**2))**2+((Ec1*Er2+Er1*Ec2)*(Er1+Er2)/((Er1+Er2)**2+(Ec1+Ec2)**2)-(Er1*Er2-Ec1*Ec2)*(Ec1+Ec2)/((Er1+Er2)**2+(Ec1+Ec2)**2))**2)**(0.5e0)*(2e0*((Er1*Er2-Ec1*Ec2)*(Er1+Er2)/((Er1+Er2)**2+(Ec1+Ec2)**2)+(Ec1*Er2+Er1*Ec2)*(Ec1+Ec2)/((Er1+Er2)**2+(Ec1+Ec2)**2))*(-Ec1*(Er1+Er2)/((Er1+Er2)**2+(Ec1+Ec2)**2)-(Er1*Er2-Ec1*Ec2)*(Er1+Er2)/((Er1+Er2)**2+(Ec1+Ec2)**2)**2*(2e0*Ec1+2e0*Ec2)+Er1*(Ec1+Ec2)/((Er1+Er2)**2+(Ec1+Ec2)**2)+(Ec1*Er2+Er1*Ec2)/((Er1+Er2)**2+(Ec1+Ec2)**2)-(Ec1*Er2+Er1*Ec2)*(Ec1+Ec2)/((Er1+Er2)**2+(Ec1+Ec2)**2)**2*(2e0*Ec1+2e0*Ec2))+2e0*((Ec1*Er2+Er1*Ec2)*(Er1+Er2)/((Er1+Er2)**2+(Ec1+Ec2)**2)-(Er1*Er2-Ec1*Ec2)*(Ec1+Ec2)/((Er1+Er2)**2+(Ec1+Ec2)**2))*(Er1*(Er1+Er2)/((Er1+Er2)**2+(Ec1+Ec2)**2)-(Ec1*Er2+Er1*Ec2)*(Er1+Er2)/((Er1+Er2)**2+(Ec1+Ec2)**2)**2*(2e0*Ec1+2e0*Ec2)+Ec1*(Ec1+Ec2)/((Er1+Er2)**2+(Ec1+Ec2)**2)-(Er1*Er2-Ec1*Ec2)/((Er1+Er2)**2+(Ec1+Ec2)**2)+(Er1*Er2-Ec1*Ec2)*(Ec1+Ec2)/((Er1+Er2)**2+(Ec1+Ec2)**2)**2*(2e0*Ec1+2e0*Ec2)))+2e0*Ec1*(Er1+Er2)/((Er1+Er2)**2+(Ec1+Ec2)**2)+(2e0*Er1*Er2-2e0*Ec1*Ec2)*(Er1+Er2)/((Er1+Er2)**2+(Ec1+Ec2)**2)**2*(2e0*Ec1+2e0*Ec2)-2e0*Er1*(Ec1+Ec2)/((Er1+Er2)**2+(Ec1+Ec2)**2)-(2e0*Ec1*Er2+2e0*Er1*Ec2)/((Er1+Er2)**2+(Ec1+Ec2)**2)+(2e0*Ec1*Er2+2e0*Er1*Ec2)*(Ec1+Ec2)/((Er1+Er2)**2+(Ec1+Ec2)**2)**2*(2e0*Ec1+2e0*Ec2))

    deltaImBeta = abs(term1*deps1r) + abs(term2*deps1c) + abs(term3*deps2r) + abs(term4*deps2c)
    
  except: 
    deltaImBeta = -1e0
  
  return deltaImBeta
#}}}

## returns the absolute uncertainty on the SPP mean free path L_SPP.
def deltaLspp(wavelength, eps1, eps2, deps1r, deps1c, deps2r, deps2c): #{{{
  deltaBeta = deltaImBetaSPP(wavelength, eps1, eps2, deps1r, deps1c, deps2r, deps2c)
  beta = betaSPP(wavelength, eps1, eps2)
  betaImag = beta.imag
  try: 
    value = 0.5e0*abs(deltaBeta)/(betaImag**2)
  except:
    value = -1e0
  return value.real


## Computes the SPP decay depth in one slab
def DecayDepth(kzSPP):#{{{ #trying with module instead of real part
  if(kzSPP.real != 0): 
      kzSPPnorm=np.sqrt(kzSPP.real**2 + kzSPP.imag**2)
  else:
      kzSPPnorm = -1
  return 1e0/kzSPPnorm
  #return 2e0*pi/kzSPP.real
#}}}

## Computes the complex wavenumber in direction of incident laser, perp. to SPP propagation. 
def kzSPP(wavelength,eps1,eps2):#{{{
  return cmath.sqrt(betaSPP(wavelength,eps1,eps2)**2-eps1*(omega(wavelength)**2/c**2))
#}}}

## Return the coherent length of SPPs
def DecayLengthSPP(beta):#{{{
    try:
      result = 1e0/(2e0*beta.imag)
    except: 
      result = -1
    return result
#}}}
  
## More elaborated functions

## Print only the experimentally possible cases: SPP active depth must be smaller than absorption depth. 
# if Opd1 > 0, then: 
#   return (Opd1 > SPPdecayDepth1)
# else: 
#   return true

## Tests if optical penetration depth is larger than SPP decay depth
# TODO: this intution was not disproved experimentally. Invent experiment to solve this. 
def ExperimentallyAchievable(OpticalPenetrationDepth, DecayDepth):#{{{
  if (OpticalPenetrationDepth != -1):
    #test if OPD > SPPdecayDepth
    return (OpticalPenetrationDepth > DecayDepth)
  else:
    return True
#}}}

## Print all the SPP-active interfaces available in database
#  CONSIDERS ONLY CASES where group velocity v_g > 0
#  TODO: rewrite this function to control better the conditional parameters (ExperimentallyAchievable, Period!=0, and SPPcondition.). 
#  If comment=="new", old SPP-active interfaces are removed from the table
def SPPactiveInterfaces(dbarray, comment):#{{{
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
        # TODO: this condition should consider several levels of accuracy. 
        # 0: compute for every materials, blindly
        # 1: compute when SPPcondition() or OldSPPcondition() return True. 
        # 2: compute only using the RegularLIPSScondition (more permissive than SPP excitation conditions)
        if(LevelOfSPPaccuracy==0): 
            Condition=True
        elif(LevelOfSPPaccuracy==1): 
            Condition=((SPPcondition(eps1,eps2)) or (OldSPPcondition(eps1,eps2)))
        elif(LevelOfSPPaccuracy==2):
            Condition=RegularLIPSScondition
        else: 
            print(Header+"** Error: level of tolerance over SPP conditions is not well indicated. See libSPP.py: LevelOfSPPaccuracy.")
            exit()
        if(Condition):
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
        if(LevelOfSPPaccuracy < 3):
          ExperimentalAchievable = True #ExperimentallyAchievable(OpticalPenetration1, SPPdecayDepth1)
        else:
          ExperimentalAchievable = ExperimentallyAchievable(OpticalPenetration1, SPPdecayDepth1)
        
        # Print the table of active SPP interfaces for all cases or only new SPP interfaces					
        #if( not (comment=="new")):
        ##TODO: introduce RigorLevels: 
        ## 0: all permisive, no verification on SPP excitation condition
        ## 1: use the RegularLIPSScondition, softer than pure SPP excitation condition
        ## 2: use the generalized SPP excitation condition
        ## 3: Extreme level: use ExperimentallyAchievable() to verify possibility of decay depth > optical penetration depth
        if(LevelOfSPPaccuracy == 1):
          Condition = (SPPdecayLength < 20000e0) #and (abs(eps2.real) < eps2.imag)
        elif (LevelOfSPPaccuracy == 2):
          Condition = (Period!=0)
        elif (LevelOfSPPaccuracy >= 3):
          Condition = ExperimentalAchievable and (Period!=0)
        else: #super permissive case
          Condition = True
		
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

## Builds a database with SPP active interfaces
# We would like now to construct a database using available materials description with all possible interfaces
# we will : 
# 
# 1. For each material in database, select each material and verify, for each available wavelength, 
# 1.1: If SPP condition is verified, 
# 1.2. yes, then period can be calculated and shown;
# 1.4. SPP decay depth in medium 1
# 1.5. SPP decay depth in medium 2
# 2. Then extract a table which contains all possible scenarios
def GenerateDatabase():
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

## calculates group velocity for a set of dispersion curves neglecting complex space
#  applies to:
#  (<array of omegas>, array of k) -> <array of group velocities>
#  
def HohenauGroupVelocity(omega, k): 
  #vg = c/(nspp - wavelength * dnspp / dwavelength)
  vg = np.diff(omega) / np.diff(k.real)
  return vg

## Complex derivative a real-valued function by a complex-number
# Input:
#   f: z->f(z)
#   z: z complex-valued numbers
# Output:
#   df/dz according to Wirtinger complex derivatives formula
# Careful: Wirtinger formula contain a 1/2 usually, but it was removed here to obtain consistent results with Hohenau lifetime and Raether lifetime. 
def RealDerivativeByComplex(f,z):
  return (np.diff(f)/np.diff(z.real) - 1j*(np.diff(f)/np.diff(z.imag))) #2 x original

## SPP lifetime as defined by H. Raether, Surface Plasmons on Smooth and Rough Surfaces and on Gratings Springer-Verlag (1986)
def LifeTimeRaether(beta, eps2, eps1):
	omegasppimag=beta.real * c * eps1.imag/(2e0*eps1.real**2) * (eps1.real * eps2.real)/(eps1.real + eps2.real)
	#lifetime=2e0*pi/omegasppimag #Raether original formula
	#lifetime=2e0/omegasppimag #Raether modified formula, we just dropped the pi. 
	lifetime=0.5e0/omegasppimag #modified Raether formula to match with complex group velocity approach
	return lifetime

## Returns the SPP decay length (in meters)
# Maier, S. A. Science, S. (Ed.) Plasmonics, Fundamentals and Applications Springer, 2007
#According to [R. Krenn et al, PRB 78, 155405 (2008)], this formula considers the 1/e decay of I_spp. 
def SPPlength(beta): #{{{
  length = 1e0/(2e0 * beta.imag) #Maier formula
  return length
#}}}

## Returns a lifetime calculated using group velocity @param vg
# Uses a lifetime based on group velocity
def LifeTimeVg(beta, vg): 
  length = SPPlength(beta)
  lifetime = length * (vg)**(-1e0)
  return lifetime
              
## Returns effective optical index of SPP
def EffectiveIndex(eps1, eps2): 
  #neff = 
  return cmath.sqrt( eps1*eps2 / (eps1+eps2) )

kzSPP = np.vectorize(kzSPP)
DecayDepth = np.vectorize(DecayDepth)
EffectiveIndex = np.vectorize(EffectiveIndex)
#omega = np.vectorize(omega) 
