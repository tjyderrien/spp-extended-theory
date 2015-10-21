#!/usr/bin/env python
#-*- coding: utf-8 -*-

# IMPORT LIBRARIES
import numpy as np
from numpy import genfromtxt, loadtxt
from scipy.optimize import fsolve, root
import cmath
import matplotlib as mp
import matplotlib.pyplot as plt
from scipy.interpolate import InterpolatedUnivariateSpline
from matplotlib import rc
# from pylab import *
from scipy.constants import c, epsilon_0, mu_0, pi, e, m_e

lengthunit=1e-9

# Settings for matplotlib
rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
## for Palatino and other serif fonts use:
#rc('font',**{'family':'serif','serif':['Palatino']})
rc('text', usetex=True)

# basic wave function
def omega(wavelength):
  return 2.0*pi*c/wavelength

# SPP BASIC FUNCTIONS
def betaSPP(wavelength, eps1, eps2):
  """calculate the SPP wave number on a flat interface
    input: wavelength (float), eps1 (complex), eps2(complex)
  """
  omega=2.0*pi*c/wavelength
  return omega/c * cmath.sqrt(eps1 * eps2 / (eps1 + eps2))
  
def AsymmetricSPPconditionPos(eps1, eps2):
	""" Assume that Re(k1).Re(k2) < 0 and verify the subsequent consequences alltogether.
	It exists then two sub-modes, let's say a positive one (k1<0, k2>0), and a negative one (k1>0, k2<0)
	"""
	value = eps1/eps2
	condition = (value.imag < 0e0)
	#condition = (eps1.imag / eps2.imag * eps2.real < eps1.real)
	return condition

def AsymmetricSPPconditionNeg(eps1, eps2):
	""" Assume that Re(k1).Re(k2) < 0 and verify the subsequent consequences alltogether.
	It exists then two sub-modes, let's say a positive one (k1<0, k2>0), and a negative one (k1>0, k2<0)
	"""
	value = eps1/eps2
	condition = (value.imag > 0e0)
	#condition = (eps1.imag / eps2.imag * eps2.real < eps1.real)
	return condition

def SPPconditionValue(eps1, eps2):
  """SPPconditionValue() returns the value of condition for SPP. If its negative, then SPP can be excited at a flat interface. 
  /!\ This condition is restricted to checking the real part of the dispersion relation for symmetric SPP only. 
    Input: eps1, eps2: complex-valued quantities
    Output: float
  """
  condition=eps1.real*eps2.real+eps1.imag*eps2.imag
  return condition

def SPPcondition(eps1, eps2):
  """ Returns a boolean claiming if SPP are excitable on an interface
  """
  if (SPPconditionValue(eps1, eps2) < 0.0):
           output=True
  else:
    output=False
  return output

def OldSPPcondition(eps1, eps2):
  """SPPconditionValue() returns the value of condition for SPP IN PERFECT MATERIALS (Im(eps)<<|Re(eps)). If its negative, then SPP can be excited at a flat interface. 
    Input: eps1, eps2: complex-valued quantities
    Output: float
  """
  condition1=(eps1.real*eps2.real<0.0)
  #condition2= eps2.real < abs(eps1.real) #this version is not symmetric, hence strange
  # Let's use its generalization which is actually symmetric. 
  condition2 = (eps1.real * eps2.real / (eps1.real + eps2.real) > 0e0)
  return (condition1 and condition2)

def period(betaSPP):
  """ Returns the period of the light-SPP field at a given interface
  """
  return 2.0*pi/betaSPP.real

### SPP decay depth
def DecayDepth(kzSPP):
  return 2e0*pi/kzSPP.real

def kzSPP(wavelength,eps1,eps2):
  return cmath.sqrt(betaSPP(wavelength,eps1,eps2)**2-eps1*(omega(wavelength)**2/c**2))

def DecayLengthSPP(beta):
	"""Return the coherent length of SPPs
	"""
	return 2e0/beta.imag
  
# OPTICAL FUNCTIONS
def Drude(wavelength, ne, epsilon, nu):
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
  return epsilon - omegap2/(omega*omega) * 1/(1+1j*nu/omega)

def reflectivity(eps1, eps2):
  """Return Fresnel reflectivity 
  Input:
    eps1: complex-valued permittivity 1+j0
    eps2: idem, for medium2
  Output: 
    interface reflectivity (float) R
  """
  R=abs(((eps1**0.5e0-eps2**0.5e0)/(eps1**0.5e0+eps2**0.5e0))**2)
  return R
  
## More elaborated functions

def ExperimentallyAchievable(OpticalPenetrationDepth, DecayDepth):
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

def SPPactiveInterfaces(dbarray, comment):
  """Print all the SPP-active interfaces available in database
  If comment=="new", old SPP-active interfaces are removed from the table
  """
  counter=0
  
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
        if ((SPPcondition(eps1,eps2)) or (OldSPPcondition(eps1,eps2))):
	        SPPperiod=(period(betaSPP(wavelength1,eps1, eps2))/lengthunit)
	        SPPdecayDepth1=(DecayDepth(kzSPP(wavelength1, eps1, eps2))/lengthunit)
	        SPPdecayDepth2=(DecayDepth(kzSPP(wavelength2, eps2, eps1))/lengthunit)
	        SPPdecayLength=DecayLengthSPP(betaSPP(wavelength1,eps1, eps2))/lengthunit
        else: 
	        SPPperiod=0
	        SPPdecayDepth1=0
	        SPPdecayDepth2=0
	        SPPdecayLength=0
        
        Reflectivity=(reflectivity(eps1, eps2))
        
        # Print only the experimentally possible cases: SPP active depth must be smaller than absorption depth. 
        # ensure that SPPdecayDepth is smaller than layer thickness, to avoid shift of dispersion relation
        ExperimentalAchievable = True #ExperimentallyAchievable(OpticalPenetration1, SPPdecayDepth1)
        
        # Print the table of active SPP interfaces for all cases or only new SPP interfaces					
        #if( not (comment=="new")): 
        if(ExperimentalAchievable and (SPPperiod!=0)):
          if (np.mod(counter, 20) == 0): 
            #show the table line each 20 lines
            if (comment):
							print '{0:12s} {1:12s} {2:15s} {3:12s} {4:12s} {5:11s} {6:16s} {7:16s} {8:12s} {9:19s} {10:19s} {11:15s}'.format("# Substrate", "Layer", "Wavelength (nm)", "OldSPPactive", "NewSPPactive", "Period (nm)", "DecayDepth1 (nm)", "DecayDepth2 (nm)", "Reflectivity", "OpticalPenetration1", "OpticalPenetration2", "DecayLength")

          counter=counter+1
          print '{0:12s} {1:12s} {2:15f} {3:12s} {4:12s} {5:11f} {6:16f} {7:16f} {8:12f}  {9:19f} {10:19f} {11:15f}'.format(Material1, Material2, Wavelength, OldSPPactiveBool, NewSPPactiveBool, SPPperiod, SPPdecayDepth1, SPPdecayDepth2, Reflectivity, OpticalPenetration1, OpticalPenetration2, SPPdecayLength)
        
  return 0

def AsymmetricSPPposActiveInterfaces(dbarray, comment):
  """Print all the SPP-active interfaces available in database
  If comment=="new", old SPP-active interfaces are removed from the table
  """
  counter=0
  
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
        if (AsymmetricSPPconditionPos(eps1, eps2)): 
	        NewSPPactiveBool='Yes'
        else: 
	        NewSPPactiveBool='No'
        # If new or old SPP active condition is true, then show	
        if ((AsymmetricSPPconditionPos(eps1,eps2)) or (OldSPPcondition(eps1,eps2))):
	        SPPperiod=(period(betaSPP(wavelength1,eps1, eps2))/lengthunit)
	        SPPdecayDepth1=(DecayDepth(kzSPP(wavelength1, eps1, eps2))/lengthunit)
	        SPPdecayDepth2=(DecayDepth(kzSPP(wavelength2, eps2, eps1))/lengthunit)
        else: 
	        SPPperiod=0
	        SPPdecayDepth1=0
	        SPPdecayDepth2=0
        
        Reflectivity=(reflectivity(eps1, eps2))
        
        # Print only the experimentally possible cases: SPP active depth must be smaller than absorption depth. 
        # ensure that SPPdecayDepth is smaller than layer thickness, to avoid shift of dispersion relation
        ExperimentalAchievable = True #ExperimentallyAchievable(OpticalPenetration1, SPPdecayDepth1)
        
        # Print the table of active SPP interfaces for all cases or only new SPP interfaces					
        #if( not (comment=="new")): 
        if(ExperimentalAchievable and (SPPperiod!=0)):
          if (np.mod(counter, 20) == 0): 
            #show the table line each 20 lines
            print '{0:12s} {1:12s} {2:15s} {3:12s} {4:12s} {5:11s} {6:16s} {7:16s} {8:12s} {9:19s} {10:19s}'.format("# Substrate", "Layer", "Wavelength (nm)", "OldSPPactive", "NewSPPactive", "Period (nm)", "DecayDepth1 (nm)", "DecayDepth2 (nm)", "Reflectivity", "OpticalPenetration1", "OpticalPenetration2")
          
          counter=counter+1
          print '{0:12s} {1:12s} {2:15f} {3:12s} {4:12s} {5:11f} {6:16f} {7:16f} {8:12f}  {9:19f} {10:19f}'.format(Material1, Material2, Wavelength, OldSPPactiveBool, NewSPPactiveBool, SPPperiod, SPPdecayDepth1, SPPdecayDepth2, Reflectivity, OpticalPenetration1, OpticalPenetration2)
        
  return 0


def AsymmetricSPPnegActiveInterfaces(dbarray, comment):
  """Print all the SPP-active interfaces available in database
  If comment=="new", old SPP-active interfaces are removed from the table
  """
  counter=0
  
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
        if (AsymmetricSPPconditionNeg(eps1, eps2)): 
	        NewSPPactiveBool='Yes'
        else: 
	        NewSPPactiveBool='No'
        # If new or old SPP active condition is true, then show	
        if ((AsymmetricSPPconditionNeg(eps1,eps2)) or (OldSPPcondition(eps1,eps2))):
	        SPPperiod=(period(betaSPP(wavelength1,eps1, eps2))/lengthunit)
	        SPPdecayDepth1=(DecayDepth(kzSPP(wavelength1, eps1, eps2))/lengthunit)
	        SPPdecayDepth2=(DecayDepth(kzSPP(wavelength2, eps2, eps1))/lengthunit)
        else: 
	        SPPperiod=0
	        SPPdecayDepth1=0
	        SPPdecayDepth2=0
        
        Reflectivity=(reflectivity(eps1, eps2))
        
        # Print only the experimentally possible cases: SPP active depth must be smaller than absorption depth. 
        # ensure that SPPdecayDepth is smaller than layer thickness, to avoid shift of dispersion relation
        ExperimentalAchievable = True #ExperimentallyAchievable(OpticalPenetration1, SPPdecayDepth1)
        
        # Print the table of active SPP interfaces for all cases or only new SPP interfaces					
        #if( not (comment=="new")): 
        if(ExperimentalAchievable and (SPPperiod!=0)):
          if (np.mod(counter, 20) == 0): 
            #show the table line each 20 lines
            print '{0:12s} {1:12s} {2:15s} {3:12s} {4:12s} {5:11s} {6:16s} {7:16s} {8:12s} {9:19s} {10:19s}'.format("# Substrate", "Layer", "Wavelength (nm)", "OldSPPactive", "NewSPPactive", "Period (nm)", "DecayDepth1 (nm)", "DecayDepth2 (nm)", "Reflectivity", "OpticalPenetration1", "OpticalPenetration2")
          
          counter=counter+1
          print '{0:12s} {1:12s} {2:15f} {3:12s} {4:12s} {5:11f} {6:16f} {7:16f} {8:12f}  {9:19f} {10:19f}'.format(Material1, Material2, Wavelength, OldSPPactiveBool, NewSPPactiveBool, SPPperiod, SPPdecayDepth1, SPPdecayDepth2, Reflectivity, OpticalPenetration1, OpticalPenetration2)
        
  return 0

def RealDerivativeByComplex(f,z):
  """Complex derivative a real-valued function by a complex-number
  Input:
    f: z->f(z)
    z: z complex-valued numbers
  Output:
    df/dz according to complex derivatives formula
  """
  return (0.5e0+0j) * (np.diff(f)/np.diff(z.real) - 1j*(np.diff(f)/np.diff(z.imag))) #original
  #return 0.5 * np.add(np.divide(np.diff(f),np.diff(z.real)), - 1j*np.divide(np.diff(f),np.diff(z.imag))) #original

def LifeTimeRaether(beta, eps2, eps1):
	omegasppimag=beta.real * c * eps1.imag/(2.*eps1.real**2) * (eps1.real * eps2)/(eps1.real + eps2)
	lifetime=1e0/(2e0*omegasppimag)
	return lifetime

def SPPlength(beta):
	return 1e0/(2e0 * beta.imag)

def LifeTimeDerrien(beta, vg):
	length = SPPlength(beta)
	lifetime = 0.5 * length * (vg)**(-1e0)
	return lifetime