#!/usr/bin/env python2
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

## @package libKeldysh-Zhukov
# This module aims to calculate the density of excited electrons as function of laser parameters. 
# Several flavors of the Keldysh theory are available: 
# - Keldysh original paper in solid, for Kane band structure [compared with td-dft]
# - Keldysh paper with few terms corrected by Gruzdev [compared with td-dft]
# - Keldysh-Zhukov tables, where Keldysh theory was computed numerically without using the saddle point method
# - Keldysh-Shcheblanov model, improving rigor on the analytical integration [https://arxiv.org/abs/1706.07303]
# - Keldysh-Corkum model, allowing for analytical treatment of mulltiwavelength fields [Physical Review Letters, 2017, 118, 173601]

# IMPORT LIBRARIES
import numpy as np
from numpy import genfromtxt, loadtxt, chararray
#from scipy.optimize import fsolve, root
from scipy.special import ellipk, ellipe, dawsn, factorial2, factorial, ellipkm1
#import cmath
import matplotlib as mp
import matplotlib.pyplot as plt
from scipy.interpolate import interp2d, InterpolatedUnivariateSpline
from matplotlib import rc
# from pylab import *
from scipy.constants import c, epsilon_0, mu_0, pi, e, m_e, h, hbar, Avogadro
#from matplotlib.legend_handler import HandlerLine2D
#import sys

from libUnits import *
from libDatabase import *

## Generates the normalization coefficients for electric field
# /!\ Vladimir uses adiabadicity coefficient for gas, which differs from a 1/sqrt(2) factor. 
# /!\ Vladimir also uses CGS units. 
# /!\ Vladimir neglects the effective mass to compute normalization of the field! #TODO Check with him!
# See the paper [Keldysh, L. Ionization in the field of a strong electromagnetic wave Journal of Experimental and Theoretical Physics, Lebedev Inst. of Physics, Moscow, 1964, 47, 5]
# @param Egap: band gap value (Joules, SI)
# @param meff: effective mass (no unit)
# @param wavelength: wavelength of the photons (in meters, SI)
# @param NormalizedPeakField: normalized field to be obtained in CGS
def VZ_FieldNormalization(Egap, meff, wavelength):
  # field for which gamma_VZ = 1.  
  me_CGS = Mass_SI_to_CGS(m_e) * meff #[1 kg    (SI) = 1E3  g      (CGS) ]
  Eg_CGS = Energy_SI_to_CGS(Egap)     #[1 J     (SI) = 1E7  ergs   (CGS) ]
  c_CGS  = Velocity_SI_to_CGS(c)      #[1 [m/s] (SI) = 1E2 cm/s   (CGS) ]
  e_CGS  = electric_charge_SI_to_CGS(e)  #[1 C     (SI) = c_CGS \times statC (CGS) ]
  wavelength_CGS = Length_SI_to_CGS(wavelength)
  omega_CGS = 2.*pi*c_CGS/wavelength_CGS #should be equal to SI
  omega_SI  = 2.*pi*c    /wavelength #SI
  if(omega_CGS != omega_SI):
    print "Error on omega_CGS"
    exit()
  EfieldStar_VZ_CGS = np.sqrt(2. * omega_CGS**2 / e_CGS**2 * me_CGS * Eg_CGS) #gas formula for Keldysh parameter
  EfieldStar_VZ_SI  = np.sqrt(2. * omega_SI**2 / e**2 * m_e * meff * Egap)
  #print EfieldStar_VZ_SI
  #Field_CGS_to_SI(EfieldStar_VZ_CGS)
  return EfieldStar_VZ_SI, EfieldStar_VZ_CGS

VZ_FieldNormalization = np.vectorize(VZ_FieldNormalization)

## Interfaces bicolor tables of V. Zhukov bicolor Keldysh model with the laser parameters
# Returns the W_PI coefficients from V. Zhukov model to be integrated in time for bi-color laser pulses
def VZ_generateWpiTables(FieldEnvelope1, FieldEnvelope2, wavelength1 = 800e-9, wavelength2 = 800e-9, CEP1=0., CEP2=0., Egap=2.56*e, meff=0.2226, tau1=10e-15, tau2=10e-15, Delay=0., dt = 1E-17, N_total=5E28, t0=0.0): #{{{
  Header="[libKeldysh] "
  
  print Header+"Defining the laser pulse..."
  # t0 = 0e0
  tmin = -1.*tau1+t0
  tmax =  1.*tau2+t0 + Delay
  
  print Header+"** Info: Number of time steps = "+str(int((tmax-tmin)/dt))+"."
  instants=np.arange(tmin,tmax,dt)
  
  # Printing info on the pulses
  IntensityEnvelope1=FieldToIntensity(FieldEnvelope1)
  print Header+"** Info: wavelength 1=", wavelength1
  print Header+"** Info: Peak intensity 1= "+str(np.max(IntensityEnvelope1)/1E4)+" W/cm^2."
  print Header+"** Info: Peak field amplitude 1= "+str(np.max(FieldEnvelope1)/1E9)+" V/nm."
  
  IntensityEnvelope2=FieldToIntensity(FieldEnvelope2)
  print Header+"** Info: wavelength 2=", wavelength2
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
  
  FieldEnvelopeNormalized2_CGSMax = FieldEnvelopeNormalized2_CGS.max()
  
  print Header+"Normalized Field1 is now [CGS]: "+str(FieldEnvelopeNormalized1_CGS.max())
  print Header+"Normalized Field2 is now [CGS]: "+str(FieldEnvelopeNormalized2_CGS.max())
  print Header+"Normalized Field1 is now [SI]: "+str(FieldEnvelopeNormalized1_SI.max())
  print Header+"Normalized Field2 is now [SI]: "+str(FieldEnvelopeNormalized2_SI.max())
  print ""
  print Header+"** Info: Egap = "+str(Egap/e)+" eV"
  print ""

  #Choosing the right column index in the data files
  print Header+"** Selecting the right headers..."
  if( FieldEnvelopeNormalized2_CGS.max() == 0. ): # was wavelength1 == wavelength2.  
    Dictionnary = {'FieldSquaredLog10': 0, 'log10wpi': 1, 'photons': 2, 'energy': 3, 'wpi': 4, 'FieldSquared1': 5} #Monochromatic case
    DataFolder  = 'Zhukov/Monochrome/'
    DataFileName={'1030': DataFolder+'DLG1030mono.dat', '800': DataFolder+'DLG800mono.dat', '400': DataFolder+'DLG400mono.dat'}
    print Header+"Choosing the right database..."
    if(wavelength1   == 800e-9):
      VZ_basename = DataFileName['800']
    elif(wavelength1 == 400e-9):
      VZ_basename = DataFileName['400']
    elif(wavelength1 == 1030e-9):
      VZ_basename = DataFileName['1030']
    else: 
      print Header+"** Warning: single color general Keldysh model is available in this library. "
  elif((wavelength1 == 800e-9 and wavelength2 == 1030e-9) or (wavelength1 == 1030e-9 and wavelength2 == 800e-9)): 
    #TODO: the two sets could be inverted! Therefore data should be swept.
    print "THIS SET IS BROKEN. Waiting for the input of Vladimir Zhukov."
    exit()
    
  elif((wavelength1 == 400e-9 and wavelength2 == 2*wavelength1) or (wavelength1 == 800e-9 and wavelength2 == wavelength1/2.)):
    #TODO: the two sets could be inverted! 
    Dictionnary={'FieldSquared1': 0, 'FieldSquared2': 1, 'wpi': 2 } #Bichromatic case
    DataFolder  = 'Zhukov/800x400nm/'
    DataFileName={'phi=0': 'W400x800fi=0.dat', 'phi=pi/4': 'W400x800fi=pina4.dat'}
    if(CEP2==0.):
      VZ_basename = DataFolder+DataFileName['phi=0']
    elif(CEP2 == pi/4.): 
      VZ_basename = DataFolder+DataFileName['phi=pi/4']
    elif(CEP2 == pi/2.):
      VZ_basename = DataFolder+DataFileName['phi=0']
    else: 
      print Header+"Fields value are not available for the specified particular case of 400x800 nm."
      exit()
  elif((wavelength1 == 800e-9 and wavelength2 == 2.*wavelength1) or (wavelength1 == 800e-9 and wavelength2 == wavelength1*2.)):
    #TODO: the two sets could be inverted! 
    Dictionnary={'FieldSquared1': 0, 'FieldSquared2': 1, 'wpi': 2} #Bichromatic case
    DataFolder  = 'Zhukov/800x1600/'
    
    DataFileName={'phi=0': 'Wpi800x1600fi=0.dat', 'phi=pi/2': 'Wpi800x1600fi=0.dat', 'phi=pi/3': 'Wpi800x1600fi=pi_over_3.dat', 'phi=pi/4': 'Wpi800x1600fi=pi_over_4.dat'}
    
    if(CEP2==0. or CEP==pi/2.):
      VZ_basename = DataFolder+DataFileName['phi=0']
    elif(CEP2 == pi/3.):
      VZ_basename = DataFolder+DataFileName['phi=pi/3']
    elif(CEP2 == pi/4.): 
      VZ_basename = DataFolder+DataFileName['phi=pi/4']
    else: 
      print Header+"Fields value are not available for 800x1600 nm."
      exit()
  else:
    print Header+"THIS COMBINATION OF WAVES IS NOT AVAILABLE. Please kindly ask the corresponding data to Prof. Vladimir Zhukov, zukov@ict.nsc.ru." 
    exit()
  print Header+"Path: "+VZ_basename
  IndexWpi           = Dictionnary['wpi']
  IndexFieldSquared1 = Dictionnary['FieldSquared1']
  if (FieldEnvelopeNormalized2_CGS.max() > 0.):
    IndexFieldSquared2 = Dictionnary['FieldSquared2']
    NumberOfColors = 2
  else:
    NumberOfColors = 1
  
  # Fetching content of the files
  databasecontents   = loadtxt(VZ_basename, skiprows=2)
  
  deltaField_CGS = 0.0025 #TODO: automatic step from the database file? Isnt it a bit small ?!
  deltaField_SI = Field_CGS_to_SI(deltaField_CGS)
  
  # Time to filter the entries with the required normaliezd field in the relevant database
  #databasecontentsfilter=FilterDatabaseLowerThan(databasecontentsfilter,FieldEnvelopeNormalized1_CGS+deltaField_CGS,IndexFieldSquared)
  #databasecontentsfilter=FilterDatabaseGreaterThan(databasecontentsfilter,FieldEnvelopeNormalized1_CGS-deltaField_CGS,IndexFieldSquared)
  DB_FieldSquaredNorm1 = databasecontents[:,IndexFieldSquared1] #this is a mapping
  DB_Wpi               = databasecontents[:,IndexWpi]           #this is a mapping
  if(NumberOfColors == 2):
    DB_FieldSquaredNorm2 = databasecontents[:,IndexFieldSquared2] #this is a mapping
  
    #print DB_FieldSquaredNorm1.shape, DB_FieldSquaredNorm2.shape, DB_Wpi.shape
    print Header+"Data well imported from DB."
    
    # This works for 2D matrix interpolation
    print Header+"** Preparing interpolation of w_pi(E1,E2)..."
    xdim = int(np.sqrt(len(DB_FieldSquaredNorm1)))
    print xdim
    Wpi_2D_X  = DB_FieldSquaredNorm1.reshape(xdim, xdim)
    Wpi_2D_Y  = DB_FieldSquaredNorm2.reshape(xdim, xdim)
    Wpi_2D    = DB_Wpi.reshape((xdim, xdim))
    Wpi_X     = Wpi_2D_X[:,0]
    Wpi_Y     = Wpi_2D_Y[0,:]
  
    print Header+"** Interpolating the w_PI..."
    #InterpolationOrder=1
    WPI_func = interp2d(Wpi_X, Wpi_Y, Wpi_2D) #interfaces w_PI with E1, E2 values, as a function in RxR. 
  
    print Header+"Range of the interpolant: "
    print FieldEnvelopeNormalized1_CGS.min(), FieldEnvelopeNormalized1_CGS.max()
    
    # Interpolating the right W_PI [CGS unit!]
    # WPI [SI] = m^-3 s^-1
    # WPI [CGS]= cm^-3.s^-1
    w_PI_CGS_matrix = WPI_func(FieldEnvelopeNormalized1_CGS**2, FieldEnvelopeNormalized2_CGS**2) #hence w_PI_CGS: This function is giving a full matrix in RxR. I just want a list indexed on instants. Maybe capture the diagonal of this matrix? This would correspond to the same instants.
    w_PI_CGS = w_PI_CGS_matrix.diagonal()
    # w_PI_CGS is in particles per cm^-3. 
  elif(NumberOfColors == 1): 
    print Header+"** Info: selected the single color Keldysh tables."
    print Header+"Linear interpolation from the single-color W_PI(E1) table... "
    InterpolationOrder=1
    print Header+"Size of field envelope: "+str(DB_FieldSquaredNorm1.shape)
    print Header+"Size of the WPI database: "+str(DB_Wpi.shape)
    WPI_func_1d = InterpolatedUnivariateSpline(DB_FieldSquaredNorm1, DB_Wpi, k=InterpolationOrder)
    w_PI_CGS = WPI_func_1d(FieldEnvelopeNormalized1_CGS**2)
    w_PI_CGS = np.clip(w_PI_CGS, 0., None)
    
  else: 
    print Header+"Number of colors is too high. Keldysh-Zhukov model is made for 2 colors. "
    exit()
  
  print Header+"range(w_PI_CGS) = ", w_PI_CGS.min(), w_PI_CGS.max()
  print Header+"dimension(w_PI_CGS) = ", w_PI_CGS.shape
  
  print Header+"** Conversion to W_PI (SI)..."
  # Normalization coefficient given by Vladimir. 
  OverallCoefficient = 128.*Egap/(pi*hbar)*N_total #CGS unit? #What is this 128 ? 
  w_PI_SI = OverallCoefficient * (Length_CGS_to_SI(1.))**-3 * w_PI_CGS 
  
  print Header+"range(w_PI_SI) = ", w_PI_SI.min(), w_PI_SI.max()
  print Header+"dimension(w_PI_SI) = ", w_PI_SI.shape
  
  #print "Waiting for Vladimir's response on the W_PI CGS unit."
  
  # Still we cannot use this data, as they don't share the same timeline. 
  # Did we introduce E1(t) and E2(t) in each dimension? 
  # Or did we introduce (E1, E2) as a mapping? 
  # Vladimir provided a mapping. We want to use it as a 
  # We shall call w_PI_SI for each instantaneous field1 and field2. 
  # 
  
  print ""
  print Header+"Temporal integration..."
  
  # Temporal integration without limiter
  
  # Just multiply array of w_PI by dt, with limited to Ntotal
  #N_excited_Keldysh = dN_excited_Keldysh.cumsum()
  #N_excited_Gruzdev = dN_excited_Gruzdev.cumsum()
  
  # Temporal integration with limiter
  dN_excited_Zhukov = np.multiply(w_PI_SI, dt)
  print dN_excited_Zhukov.cumsum()
  # Initial number of electrons in conduction band
  N_initial = np.zeros(w_PI_SI.shape)
  
  ExpArg_Zhukov = np.divide(dN_excited_Zhukov.cumsum(), N_total)
  
  N_excited_Zhukov = np.multiply( np.exp(-ExpArg_Zhukov), N_total * np.exp(ExpArg_Zhukov) - N_total + N_initial)
  
  #N_excited_Gruzdev = np.multiply( np.exp(-ExpArg_Gruzdev), N_total * np.exp(ExpArg_Gruzdev) - N_total + N_initial)
  return instants, N_excited_Zhukov, w_PI_SI, w_PI_CGS
#}}}