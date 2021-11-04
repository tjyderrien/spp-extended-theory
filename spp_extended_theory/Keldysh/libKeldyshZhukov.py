#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Copyright (C) 2013-2021 T. J.-Y. Derrien
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
# from numpy import loadtxt
#from scipy.optimize import fsolve, root
#from scipy.misc import factorial2, factorial
#import cmath

#from matplotlib.legend_handler import HandlerLine2D
#import sys

# from scipy.optimize import fsolve, root
# from scipy.misc import factorial2, factorial
# import cmath
import matplotlib as mp
import matplotlib.pyplot as plt
from matplotlib import rc
# IMPORT LIBRARIES
from numpy import loadtxt
# from pylab import *
from scipy.constants import epsilon_0, e, m_e
from scipy.interpolate import interp2d, InterpolatedUnivariateSpline

from spp_extended_theory.Keldysh.libUnits import *
from spp_extended_theory.Libs.libDatabase import *
from   octopus_slabs.Libs.libLogging         import init_logger
logger = init_logger(__name__, verbose=True) #"plotFinalQuantities")

# from matplotlib.legend_handler import HandlerLine2D
import os

OctopusData=os.environ["QuantumLaPruns"]
OctopusSources=os.environ["QuantumLaPsources"]
sppextendedtheory=os.environ["spp_extended_theory"]

PathPrefix = sppextendedtheory+"/spp_extended_theory/Keldysh/" #PathPrefix no longer needs to be defined in octopus-slabs, but is inherent to libKeldyshZhukov.
logger.debug("Path prefix="+PathPrefix)

rc('font',**{'family':'sans-serif','sans-serif':['Helvetica'], 'size':'14'})
## for Palatino and other serif fonts use:
#rc('font', **{'family':'serif', 'serif':['Palatino'], 'size':'16'})
rc('text', usetex=True) #True: Does not work on Draco. Just put False, then, and it will be smooth. 
mp.rcParams['legend.numpoints'] = 1

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
    print("** Error: unexpected imaginary part in the conversion from Intensity to Field!")
    exit(-1)
  return Field.real


## Generates the normalization coefficient for electric field E*. 
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
    print("Error on omega_CGS")
    exit()
  EfieldStar_VZ_CGS = np.sqrt(2. * omega_CGS**2 / e_CGS**2 * me_CGS * Eg_CGS) #gas formula for Keldysh parameter
  EfieldStar_VZ_SI  = np.sqrt(2. * omega_SI**2 / e**2 * m_e * meff * Egap)
  #print EfieldStar_VZ_SI
  #Field_CGS_to_SI(EfieldStar_VZ_CGS)
  return EfieldStar_VZ_SI, EfieldStar_VZ_CGS

VZ_FieldNormalization = np.vectorize(VZ_FieldNormalization)

## Defines the interface to the simulation data provided by VP Zhukov. Format of the data was not systematic, hence we had to define a dictionnary for reading each produced files. In bicolor datasets, fields E1 and E2 were sometimes inverted. 
# @param FieldEnvelope1: Electric field 1 (SI)
# @param FieldEnvelope2: Electric field 2 (SI)
# @param wavelength1: wavelength pulse 1 (SI)
# @param wavelength2: wavelength pulse 2 (SI)
# @param CEP1: phase pulse 1 (SI)
# @param CEP2: phase pulse 2 (SI)
# @param Egap: band gap (SI) in Joules. 
# @param meff: effective mass (no unit)
# @param PathPrefix: string to be added before all employed machine paths
def VP_ChooseLibrary(FieldEnvelope1, FieldEnvelope2, wavelength1, wavelength2, CEP1, CEP2, Egap=2.56*e, meff=0.2226, PathPrefix=""):
    #print "Path prefix: "+PathPrefix
    Header="[libKeldyshZhukov: VP_ChooseLibrary: ]"
    print("** Info: selected wavelength: "+str(wavelength1*1E9)+" nm.")
    InvertedFields=True
    VZ_basename = ''
    Header = "[libKeldyshZhukov] VP_ChooseLibrary: "
    if( np.max(FieldEnvelope2) < 1e-3 ): #single pulse mode
        
        DataFolder  = PathPrefix+'Zhukov/Monochrome/'
        DataFileName={'1030': DataFolder+'2_56ev1030mkm.dat', 
                      '800': DataFolder+'2_56ev800tot.dat',
                      #'800': DataFolder+'DLG800mono.dat', 
                      '400': DataFolder+'DLG400mono.dat', 
                      '1600': DataFolder+'2_56ev1600mkm-A.dat',
                      '3200': DataFolder+'2_56_3200new.dat'}
        print(Header+"Choosing the right database...")
        if(wavelength1   == 800e-9):
            VZ_basename = DataFileName['800']
            Dictionnary = {'FieldSquaredLog10': 0, 'log10wpi': 1}
            #Dictionnary = {'FieldSquaredLog10': 0, 'log10wpi': 1, 'photons': 2, 'energy': 3, 'wpi': 4, 'FieldSquared1': 5} #Monochromatic case
        elif(wavelength1 == 400e-9):
            VZ_basename = DataFileName['400']
            Dictionnary = {'FieldSquaredLog10': 0, 'log10wpi': 1, 'photons': 2, 'energy': 3, 'wpi': 4, 'FieldSquared1': 5} #Monochromatic case
        elif(wavelength1 == 1030e-9):
            VZ_basename = DataFileName['1030']
            #Dictionnary = {'FieldSquaredLog10': 0, 'log10wpi': 1, 'photons': 2, 'energy': 3, 'wpi': 4, 'FieldSquared1': 5} #Monochromatic case
            Dictionnary = {'FieldSquaredLog10': 0, 'log10wpi': 1}
        elif(wavelength1 == 1600e-9): 
            VZ_basename = DataFileName['1600']
            Dictionnary = {'FieldSquaredLog10': 0, 'log10wpi': 1}
        elif(wavelength1 == 3200e-9): 
            VZ_basename = DataFileName['3200']
            Dictionnary = {'FieldSquaredLog10': 0, 'log10wpi': 1}
        else: 
            print(Header+"** Warning: numerical integration from VP Zhukov is not available. Please use single color Keldysh-Gruzdev model, that is available in this library. ")
            Dictionnary = {}
    elif((wavelength1 == 800e-9 and wavelength2 == 1030e-9) or (wavelength1 == 1030e-9 and wavelength2 == 800e-9)): #BicolorCase.
        #TODO: the two sets could be inverted! Therefore data should be swept.
        InvertedFields=True
        Dictionnary={'FieldSquared1': 0, 'FieldSquared2': 1, 'wpi': 2} #Bichromatic case
        DataFolder  = PathPrefix+'Zhukov/800x1030/'
        
        DataFileName={'phi=0': 'Wpi800x1600fi=0.dat', 'phi=pi/2': 'Wpi800x1600fi=0.dat', 'phi=pi/3': 'Wpi800x1600fi=pi_over_3.dat', 'phi=pi/4': 'Wpi800x1600fi=pi_over_4.dat'}
        
        if(CEP2==0. or CEP2==pi/2.):
            VZ_basename = DataFolder+DataFileName['phi=0']
        elif(CEP2 == pi/3.):
            VZ_basename = DataFolder+DataFileName['phi=pi/3']
        elif(CEP2 == pi/4.): 
            VZ_basename = DataFolder+DataFileName['phi=pi/4']
        else: 
            print(Header+"Fields value are not available for 800x1600 nm.")
            exit()
        print("THIS SET IS BROKEN. Waiting for the input of Vladimir Zhukov.")
        exit()
        
    elif((wavelength1 == 400e-9 and wavelength2 == 2*wavelength1) or (wavelength1 == 800e-9 and wavelength2 == wavelength1/2.)):
        #TODO: the two sets could be inverted! 
        InvertedFields=True
        Dictionnary={'FieldSquared1': 0, 'FieldSquared2': 1, 'wpi': 2 } #Bichromatic case
        DataFolder  = PathPrefix+'Zhukov/800x400nm/'
        DataFileName={'phi=0': 'W400x800fi=0.dat', 'phi=pi/4': 'W400x800fi=pina4.dat'}
        if(CEP2==0.):
            VZ_basename = DataFolder+DataFileName['phi=0']
        elif(CEP2 == pi/4.): 
            VZ_basename = DataFolder+DataFileName['phi=pi/4']
        elif(CEP2 == pi/2.):
            VZ_basename = DataFolder+DataFileName['phi=0']
        else: 
            print(Header+"Fields value are not available for the specified particular case of 400x800 nm.")
            exit()
    elif((wavelength1 == 800e-9 and wavelength2 == 2.*wavelength1) or (wavelength1 == 800e-9 and wavelength2 == wavelength1*2.)):
        #TODO: the two sets could be inverted! 
        InvertedFields=False
        Dictionnary={'FieldSquared1': 0, 'FieldSquared2': 1, 'wpi': 2} #Bichromatic case
        DataFolder  = PathPrefix+'Zhukov/800x1600/'
        
        DataFileName={'phi=0': 'Wpi800x1600fi=0.dat', 'phi=pi/2': 'Wpi800x1600fi=0.dat', 'phi=pi/3': 'Wpi800x1600fi=pi_over_3.dat', 'phi=pi/4': 'Wpi800x1600fi=pi_over_4.dat'}
        
        if(CEP2==0. or CEP2==pi/2.):
            VZ_basename = DataFolder+DataFileName['phi=0']
        elif(CEP2 == pi/3.):
            VZ_basename = DataFolder+DataFileName['phi=pi/3']
        elif(CEP2 == pi/4.): 
            VZ_basename = DataFolder+DataFileName['phi=pi/4']
        else: 
            print(Header+"Fields value are not available for 800x1600 nm.")
            exit()
    else:
        print(Header+"THIS COMBINATION OF WAVES IS NOT AVAILABLE. Please kindly ask the corresponding data to Prof. Vladimir Zhukov, zukov@ict.nsc.ru.")
        Dictionnary={}
        #exit() 
    return VZ_basename, Dictionnary, InvertedFields

## Provide bicolor tables of V. Zhukov bicolor Keldysh model for the selected wavelengths
# Returns the W_PI coefficients from V. Zhukov model to be integrated in time for bi-color laser pulses
# @param FieldEnvelope1: field1 (SI units) where Wpi will be interpolated at.
# @param FieldEnvelope2: field2 (SI units) where Wpi will be interpolated at.
# @param wavelength1: 1st color. Not order sensitive. 
# @param wavelength2: 2nd color. Not order sensitive. 
# @param CEP1: should stay to 0
# @param CEP2: can change to simple values (pi/2, pi/3, pi/4)
# @param Egap: value in Joules
# @param meff: effective mass (no dimension)
def VZ_generateWpiTables(FieldEnvelope1, FieldEnvelope2, wavelength1 = 800e-9, wavelength2 = 800e-9, CEP1=0., CEP2=0.,
                         Egap=2.56*e, meff=0.2226, crystal_density=5E28, PathPrefix=PathPrefix): #{{{
  Header="[libKeldyshZhukov] VZ_generateWpiTables: "
  
  # Printing info on the pulses
  IntensityEnvelope1=FieldToIntensity(FieldEnvelope1)
  print(Header+"** Info: wavelength 1=", wavelength1)
  print(Header+"** Info: Peak intensity 1= "+str(np.max(IntensityEnvelope1)/1E4)+" W/cm^2.")
  print(Header+"** Info: Peak field amplitude 1= "+str(np.max(FieldEnvelope1)/1E9)+" V/nm.")
  
  IntensityEnvelope2=FieldToIntensity(FieldEnvelope2)
  print(Header+"** Info: wavelength 2=", wavelength2)
  print(Header+"** Info: Peak intensity 2= "+str(np.max(IntensityEnvelope2)/1E4)+" W/cm^2.")
  print(Header+"** Info: Peak field amplitude 2= "+str(np.max(FieldEnvelope2)/1E9)+" V/nm.")

  # 0. Expressing Fields in CGS
  FieldNormalizationCoeff1_SI, FieldNormalizationCoeff1_CGS = VZ_FieldNormalization(Egap, meff, wavelength1)
  FieldNormalizationCoeff2_SI, FieldNormalizationCoeff2_CGS = VZ_FieldNormalization(Egap, meff, wavelength2)
  
  print(Header+"** Normalization coefficient Field 1 [SI]: "+str(FieldNormalizationCoeff1_SI))
  print(Header+"** Normalization coefficient Field 1 [CGS]: "+str(FieldNormalizationCoeff1_CGS))
  print(Header+"** Normalization coefficient Field 2 [SI]: "+str(FieldNormalizationCoeff2_SI))
  print(Header+"** Normalization coefficient Field 2 [CGS]: "+str(FieldNormalizationCoeff2_CGS))
  
  # Normalize FieldEnvelope 1,2 in CGS. Vladimir requires normalized field1 and normalized field2 to deliver a W_PI. Note: his formula for gamma is the one for gas and does not account for optical Stark effect (increase of gap with field strength). #TODO: Why ? Stark effect also happens in gas. 
  
  FieldEnvelope1_CGS = Field_SI_to_CGS(FieldEnvelope1)
  FieldEnvelope2_CGS = Field_SI_to_CGS(FieldEnvelope2)
  
  print(Header+"Field1 [CGS] = "+str(np.max(FieldEnvelope1_CGS)))
  print(Header+"Field2 [CGS] = "+str(np.max(FieldEnvelope2_CGS)))
  
  FieldEnvelopeNormalized1_CGS = FieldEnvelope1_CGS / FieldNormalizationCoeff1_CGS
  FieldEnvelopeNormalized2_CGS = FieldEnvelope2_CGS / FieldNormalizationCoeff2_CGS
  
  FieldEnvelopeNormalized1_SI = FieldEnvelope1 / FieldNormalizationCoeff1_SI
  FieldEnvelopeNormalized2_SI = FieldEnvelope2 / FieldNormalizationCoeff2_SI
  
  FieldEnvelopeNormalized2_CGSMax = FieldEnvelopeNormalized2_CGS.max()
  
  print(Header+"Normalized Field1 is now [CGS]: "+str(FieldEnvelopeNormalized1_CGS.max()))
  print(Header+"Normalized Field2 is now [CGS]: "+str(FieldEnvelopeNormalized2_CGS.max()))
  print(Header+"Normalized Field1 is now [SI]: "+str(FieldEnvelopeNormalized1_SI.max()))
  print(Header+"Normalized Field2 is now [SI]: "+str(FieldEnvelopeNormalized2_SI.max()))
  print("")
  print(Header+"** Info: Egap = "+str(Egap/e)+" eV")
  print("")

  # 1. Choosing the right data file
  print(Header+"** Selecting the right VP Zhukov datafile...")
  VZ_basename, Dictionnary, InvertedFields=VP_ChooseLibrary(FieldEnvelope1, FieldEnvelope2, wavelength1, wavelength2,
                                                            CEP1, CEP2, 2.56*e, 0.2226, PathPrefix)
  
  logger.debug("Path: "+VZ_basename)
  #IndexWpi           = Dictionnary['wpi']
  #IndexLogWpi        = Dictionnary['log10wpi']
  #IndexFieldSquared1 = Dictionnary['FieldSquared1'] 
  if (FieldEnvelopeNormalized2_CGS.max() > 0.):
    #IndexFieldSquared2 = Dictionnary['FieldSquared2']
    NumberOfColors = 2
  else:
    NumberOfColors = 1
  
  # Fetching content of the files
  if(VZ_basename != ''):
    databasecontents   = loadtxt(VZ_basename, skiprows=2)
  else: 
    print(Header+"Path was empty.")
    return 0 #we leave the function
  
  #deltaField_CGS = 0.0025 #TODO: automatic step from the database file? 
  #deltaField_SI  = Field_CGS_to_SI(deltaField_CGS)
  
  # Time to filter the entries with the required normaliezd field in the relevant database
  #databasecontentsfilter=FilterDatabaseLowerThan(databasecontentsfilter,FieldEnvelopeNormalized1_CGS+deltaField_CGS,IndexFieldSquared)
  #databasecontentsfilter=FilterDatabaseGreaterThan(databasecontentsfilter,FieldEnvelopeNormalized1_CGS-deltaField_CGS,IndexFieldSquared)
  
  ## VP Zhukov provided non-systematic files. We prioritize the information from Wpi. 
  #if(wavelength1 in [1600e-9, 3200e-9]):
  DB_Wpi           = np.power(10.,databasecontents[:,Dictionnary['log10wpi']])
  #else: 
      #DB_Wpi           = databasecontents[:,Dictionnary['wpi']]           #this is a mapping
      
  
  if(NumberOfColors == 2): #{{{
    DB_FieldSquaredNorm1 = databasecontents[:,Dictionnary['FieldSquared1']] #this is a mapping
    DB_FieldSquaredNorm2 = databasecontents[:,Dictionnary['FieldSquared2']] #this is a mapping
  
    #print DB_FieldSquaredNorm1.shape, DB_FieldSquaredNorm2.shape, DB_Wpi.shape
    print(Header+"Data well imported from DB.")
    
    # Do we need a 2D matrix interpolation? 
    # This works for 2D matrix interpolation
    print(Header+"** Preparing interpolation of w_pi(E1,E2)...")
    xdim = int(np.sqrt(len(DB_FieldSquaredNorm1))) #NOTE: assumes data are square. True for 800+1600nm, CEP=0 at least.
    print(xdim)
    Wpi_2D_X  = DB_FieldSquaredNorm1.reshape(xdim, xdim)
    #print Wpi_2D_X
    Wpi_2D_Y  = DB_FieldSquaredNorm2.reshape(xdim, xdim)
    #print Wpi_2D_Y
    Wpi_2D    = DB_Wpi.reshape((xdim, xdim))
    #print Wpi_2D
    Wpi_X     = Wpi_2D_X[:,0]
    #print Wpi_X
    Wpi_Y     = Wpi_2D_Y[0,:]
    #print Wpi_Y
    Wpi_2D_t = np.transpose(Wpi_2D)
    print(Header+"** Interpolating the w_PI...")
    #InterpolationOrder=1
    # WARNING: the 800x1600 files have #1: E2**1 !! and #2: E1**2 ! Therefore we had to transpose the matrix. 
    if(InvertedFields):
        WPI_func = interp2d(Wpi_X, Wpi_Y, Wpi_2D_t, fill_value=0) #interfaces w_PI with E1, E2 values, as a function in RxR. 
    else:
        WPI_func = interp2d(Wpi_X, Wpi_Y, Wpi_2D, fill_value=0)
    print(Header+"Range of the interpolant: ")
    print(FieldEnvelopeNormalized1_CGS.min(), FieldEnvelopeNormalized1_CGS.max())
    print(FieldEnvelopeNormalized2_CGS.min(), FieldEnvelopeNormalized2_CGS.max())
    
    print(Header+"Trying to 2D-interpolate at one (E1,E2) point: ")
    print(WPI_func(FieldEnvelopeNormalized1_CGS**2, FieldEnvelopeNormalized2_CGS**2))
    # NOTE: The code works until here. 
    
    # Interpolating the right W_PI [CGS unit!]
    # WPI [SI] = m^-3 s^-1
    # WPI [CGS]= cm^-3.s^-1
    w_PI_CGS = WPI_func(FieldEnvelopeNormalized1_CGS**2, FieldEnvelopeNormalized2_CGS**2) 
    # w_PI_CGS is in particles per cm^-3. 
  elif(NumberOfColors == 1): 
    if(wavelength1 in [800e-9, 1030e-9, 1030e-9, 1600E-9, 3200E-9]):
        DB_FieldSquaredNorm1 = np.power(10.,databasecontents[:,Dictionnary['FieldSquaredLog10']]) #this is a mapping
    else:
        DB_FieldSquaredNorm1 = databasecontents[:,Dictionnary['FieldSquared1']] #this is a mapping
    print(Header+"** Info: selected the single color Keldysh tables.")
    print(Header+"Linear interpolation from the single-color W_PI(E1) table... ")
    InterpolationOrder=1
    print(Header+"Size of field envelope: "+str(DB_FieldSquaredNorm1.shape))
    print(Header+"Size of the WPI database: "+str(DB_Wpi.shape))
    WPI_func_1d = InterpolatedUnivariateSpline(DB_FieldSquaredNorm1, DB_Wpi, k=InterpolationOrder, ext='zeros')
    w_PI_CGS = WPI_func_1d(FieldEnvelopeNormalized1_CGS**2)
    #w_PI_CGS = np.clip(w_PI_CGS, 0., None)
    
  else: 
    print(Header+"Number of colors is too high. Keldysh-Zhukov model is made for 2 colors. ")
    exit()
  #}}}
  print(Header+"range(w_PI_CGS) = ", w_PI_CGS.min(), w_PI_CGS.max())
  print(Header+"dimension(w_PI_CGS) = ", w_PI_CGS.shape)
  
  print(Header+"** Conversion to W_PI (SI)...")
  # Normalization coefficient given by Vladimir. 
  # OverallCoefficient = 128.*Egap/(pi*hbar)*N_total #CGS unit?
  w_PI_SI = crystal_density * w_PI_CGS * 1E15 #as simple as this, according to VP Zhukov. 
  
  print(Header+"range(w_PI_SI) = ", w_PI_SI.min(), w_PI_SI.max())
  print(Header+"dimension(w_PI_SI) = ", w_PI_SI.shape)
  
  #print "Waiting for Vladimir's response on the W_PI CGS unit."
  
  # Still we cannot use this data, as they don't share the same timeline. 
  # Did we introduce E1(t) and E2(t) in each dimension? 
  # Or did we introduce (E1, E2) as a mapping? 
  # Vladimir provided a mapping. We want to use it as a 
  # We shall call w_PI_SI for each instantaneous field1 and field2. 
  # 
  
  return w_PI_SI
#}}}

## Perform temporal integration on the VZ_Bicolor tables for a given pulse. 
def VP_BicolorNexc(tau1, tau2, t0, Delay):
  print(Header+"Defining the laser pulse...")
  # t0 = 0e0
  tmin = -1.*tau1+t0
  tmax =  1.*tau2+t0 + Delay
  
  print(Header+"** Info: Number of time steps = "+str(int((tmax-tmin)/dt))+".")
  instants=np.arange(tmin,tmax,dt)
  
  print(Header+"Temporal integration...")
  
  
  
  # Just multiply array of w_PI by dt, with limited to Ntotal
  #N_excited_Keldysh = dN_excited_Keldysh.cumsum()
  #N_excited_Gruzdev = dN_excited_Gruzdev.cumsum()
  
  # Temporal integration with limiter
  dN_excited_Zhukov = np.multiply(w_PI_SI, dt)
  print(dN_excited_Zhukov.cumsum())
  # Initial number of electrons in conduction band
  N_initial = np.zeros(w_PI_SI.shape)
  
  ExpArg_Zhukov = np.divide(dN_excited_Zhukov.cumsum(), N_total)
  
  N_excited_Zhukov = np.multiply( np.exp(-ExpArg_Zhukov), N_total * np.exp(ExpArg_Zhukov) - N_total + N_initial)
  
  #N_excited_Gruzdev = np.multiply( np.exp(-ExpArg_Gruzdev), N_total * np.exp(ExpArg_Gruzdev) - N_total + N_initial)
  return instants, N_excited_Zhukov

## Test function for fetching a value from files provided by VP Zhukov. 
def VPZ_Wpi0D(): #{{{
# Trying to extract one single value from Zhukov files. 
    Efield1=1E9; Efield2=0; wavelength1=1600e-9; wavelength2=1600e-9; CEP1=0; CEP2=0; 
    Egap = 2.56*e; meff=0.2226
    VZ_generateWpiTables(Efield1, Efield2, wavelength1, wavelength2, CEP1, CEP2, Egap, meff)
#}}}

## Test function for fetching a batch of values from files provided by VP Zhukov. 
def VPZ_Wpi1D(): #{{{
    # Now trying to extract a sequence of values from Zhukov files. 
    Egap = 2.56*e; meff=0.2226
    Efield2=0.0E0; 
    
    wavelength1=3200e-9; CEP1=0; 
    wavelength2=1600e-9; CEP2=0; 
    
    Efield1_log = np.linspace(8,10,100)
    Efield1 = np.power(10.,Efield1_log)
    
    Wpi = VZ_generateWpiTables(Efield1, Efield2, wavelength1, wavelength2, CEP1, CEP2, Egap, meff)
    # NOTE: add a warning when interpolation occurs OUT of boundaries!

    plt.figure()
    plt.loglog(Efield1, Wpi, 'o-')
    plt.xlabel(r"$E_1$ (V/m)")
    plt.ylabel(r"$w_{PI}$ (m$^{-3}$.s$^{-1}$)")
    plt.show()
#}}}


## Plot the (E1, E2) map of W_PI according to VP Zhukov theories.
def VPZ_Wpi2D(wavelength1, wavelength2, Efield1_max, Efield2_max, CEP2=0., Egap=2.56*e, meff=0.2226):
    # Now trying to extract a sequence of values from Zhukov files. 
    #Egap = 2.56*e; meff=0.2226
    #wavelength1=800e-9; wavelength2=1600e-9; CEP1=0; CEP2=0; 
    CEP1=0.
    #Efield1_max = 0.8E10
    #Efield2_max = 0.4E10
    Efield1_log = np.linspace(8,np.log10(Efield1_max),800)
    Efield1 = np.power(10.,Efield1_log)
    Efield2_log = np.linspace(8,np.log10(Efield2_max),800)
    Efield2 = np.power(10.,Efield2_log)

    Wpi = VZ_generateWpiTables(Efield1, Efield2, wavelength1, wavelength2, CEP1, CEP2, Egap, meff)
    # NOTE: add a warning when interpolation occurs OUT of boundaries!

    print(np.shape(Wpi))

    QuantityTitle = r"$w_{PI}$ (Zhukov)"
    filename      = "Zhukov_"+str(wavelength1*1E9)+"nm-"+str(wavelength2*1E9)+"nm-CEP2-"+str(round(CEP2/pi,2))
    
    #heatmap, xedges, yedges = np.histogram2d(Efield1, Efield2, bins=(np.size(Efield1), np.size(Efield2)))
    #extent = [xedges[0], xedges[-1], yedges[0], yedges[-1]]
    #print extent
    #plt.imshow(Wpi, extent=extent, interpolation='nearest')
    #plot2dHeatMap(Efield1, Efield2, Wpi, QuantityTitle, filename)
    #plt.plot(Efield2, Wpi)
    
    plt.figure()
    cmap = mp.cm.get_cmap(name='Blues', lut=None)
    plt.xlabel(r"$E_1$ (V/nm)")
    plt.ylabel(r"$E_2$ (V/nm)")
    plt.title(r"$w_{PI}(\lambda_1=$"+str(wavelength1*1E9)+r"$,\lambda_2=$"+str(wavelength2*1E9)+r"$)$, [m$^{-3}$.s$^{-1}$], $\phi=$"+str(round(CEP2/pi,2))+r"$\pi$")
    plt.contourf(Efield1, Efield2, Wpi, cmap=cmap)
    plt.colorbar()
    plt.tight_layout()
    plt.savefig(filename+".eps")
    plt.savefig(filename+".png")
    plt.show()

#VPZ_Wpi0D()
#VPZ_Wpi1D()
# VPZ_Wpi2D(800e-9, 400e-9, 9E9, 0E0, 0.)
#VPZ_Wpi2D(800e-9, 400e-9, 9E9, 1.6E10, 0.)
#VPZ_Wpi2D(800e-9, 400e-9, 9E9, 1.6E10, pi/4.)
#VPZ_Wpi2D(800e-9, 400e-9, 9E9, 1.6E10, pi/3.)
#VPZ_Wpi2D(800e-9, 1600e-9, 4E9, 2E9, 0.)
#VPZ_Wpi2D(800e-9, 1600e-9, 4E9, 2E9, pi/4.)
#VPZ_Wpi2D(800e-9, 1600e-9, 4E9, 2E9, pi/3.)
#VPZ_Wpi2D(800e-9, 1030e-9, 0.8E10, 1E10)
