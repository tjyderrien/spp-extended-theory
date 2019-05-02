#!/usr/bin/env python2.7
#-*- coding: utf-8 -*-
## @package libStark
## Computes the Stark effect on a band structure at gamma points
# This module aims at computing the band gap energy as function of the average laser fieled induced by the Stark effect

import numpy as np
import numpy.linalg as LA
import matplotlib.pyplot as plt
from scipy.constants import c, epsilon_0, mu_0, pi, e, m_e, h, hbar

import math, cmath
from scipy.optimize import root
from time import time

from libAtomicUnits import *

## Provides the change of the branch electronic levels
# From simple Floquet Hamiltonian on constant pulse of frequency omega, the shift of 6 bands with the electric field is given. The eigen values have been computed from the Hamiltonian given in the Nano Letters. 
def Stark2bands1photon_EnergyShift_modified_exact(Efield_AU, omega_AU, E_gap_AU, DME_AU=1): #{{{                                                                                            
    M         = DME_AU
    Mbar      = np.conj(DME_AU)
    
    E1 =  0.5*np.sqrt(2.*Efield_AU**2*M*Mbar+E_gap_AU**2)
    
    E2 = -E1
    
    Tmp = np.sqrt(
            Mbar**2*Efield_AU**4*M**2 + 16.*Efield_AU**2*M*Mbar*omega_AU**2+64.*omega_AU**4+16.*omega_AU**2*E_gap_AU**2
            )
    
    E3 = 0.5 * np.sqrt(Tmp
        + Efield_AU**2*M*Mbar+8.*omega_AU**2+E_gap_AU**2
        )
        
    E4 = - E3
    E5 = 0.5 * np.sqrt(-Tmp
        + Efield_AU**2*M*Mbar+8.*omega_AU**2+E_gap_AU**2
        )
    
    E6 = -E5
    
    return E1, E2, E3, E4, E5, E6
#}}}

### Solves the characteristic polynom for the eigen values of band gap modification    
def Stark2bands1photon_EnergyShift_polynom_numerical(eigen, Efield_AU, omega_AU, E_gap_AU, DME_AU=1):
    M         = DME_AU
    Mbar      = np.conj(DME_AU)
    
    A6 = 1
    
    A4 = -Efield_AU**2*M*Mbar-0.75*E_gap_AU**2+3*omega_AU**2
    
    A3 = 2*omega_AU**3
    
    A2 = 0.5*Efield_AU**2*Mbar*M*E_gap_AU**2 + 0.5*omega_AU**2*E_gap_AU**2 + 0.25*Efield_AU**4*Mbar**2*M**2+3./16.*E_gap_AU**4+1.5*Efield_AU**2*M*Mbar*omega_AU**2
    
    A1 = -0.5*Efield_AU**2+M*Mbar*omega_AU**3-0.5*omega_AU**3*E_gap_AU**2
    
    A0 = - 1./64.*E_gap_AU**6     -1./16.*Efield_AU**2*M*Mbar*E_gap_AU**4     +1./16.*omega_AU**2*E_gap_AU**4     +1./8.*Efield_AU**2*M*Mbar*omega_AU**2*E_gap_AU**2     - 1./16. * E_gap_AU**2*Efield_AU**4*Mbar**2*M**2     - 1./16. * omega_AU**2*Efield_AU**4*Mbar**2*M**2
    
    return A6 * eigen**6 + A4 * eigen**4 + A3 * eigen**3 + A2*eigen**2 + A1 * eigen + A0

def Stark2bands1photon_EnergyShift_notcorrected_numerical(Efield_AU, omega_AU, E_gap_AU, DME_AU=1, x_steps=50):
    # Solver parameters
    x_min = -1E15       #-1E10
    x_max = 1E15       #1E10
    #x_steps = 30   #70
    
    # Invoking a multiple root solver similar to the one developed by F. Preucil
    # findroots(eps1, eps2, eps3, wavelength, t, x_min, x_max, y_min, y_max, x_steps, y_steps, num_of_maxs=1):
    tol_merge = 1E-10 #1E3
    merge_treshold = 1 #4? 1: helps to not miss some modes
    roots = []
    branches = []
    unique = []
    start = time()
    for x in np.linspace(x_min, x_max, num=x_steps):
        nrt = root(Stark2bands1photon_EnergyShift_polynom_numerical, (x), args=(Efield_AU, omega_AU, E_gap_AU, DME_AU), method='hybr')
        #print nrt
        #checking = func(nrt.x, eps1, eps2, eps3, k0, t, sgn1, sgn2)
        if (nrt.success):
            roots.append(nrt.x)
        
    return np.unique(roots)
    
def Stark2bands1photon_Cropped_EnergyShift_notcorrected_exact(Efield_AU, omega_AU, E_gap_AU, DME_AU=1): #{{{
    M         = DME_AU
    Mbar      = np.conj(DME_AU)
    
    Tmp = omega_AU**2+Efield_AU**2*M*Mbar+E_gap_AU**2
    E1  = 0.5*omega_AU + 0.5*np.sqrt(Tmp+2*E_gap_AU*omega_AU)
    E2  = 0.5*omega_AU - 0.5*np.sqrt(Tmp+2*E_gap_AU*omega_AU)
    E3  = 0.5*omega_AU + 0.5*np.sqrt(Tmp-2*E_gap_AU*omega_AU)
    E4  = 0.5*omega_AU - 0.5*np.sqrt(Tmp-2*E_gap_AU*omega_AU)
    
    EgapShift_AU = E4 - E2
    
    return E1, E2, E3, E4, EgapShift_AU
#}}}   


### Solves the characteristic polynom for the eigen values of band gap modification    
def Stark2bands1photon_Cropped_EnergyShift_polynom_numerical(eigen, Efield_AU, omega_AU, E_gap_AU, DME_AU=1):
    M         = DME_AU
    Mbar      = np.conj(DME_AU)
    
    A4 = 1.
    
    A3 = -2*omega_AU
    
    A2 = -0.5*Efield_AU**2*M*Mbar-0.5*E_gap_AU**2+omega_AU**2
    
    A1 = 0.5*Efield_AU**2*M*Mbar*omega_AU+0.5*E_gap_AU**2*omega_AU
    
    A0 = 1./8. * Efield_AU**2*Mbar*E_gap_AU**2*M + 1./16.*Efield_AU**4*M**2*Mbar**2 + 1./16.*E_gap_AU**4 - 0.25*E_gap_AU**2*omega_AU**2
    
    return A4 * eigen**4 + A3 * eigen**3 + A2*eigen**2 + A1 * eigen + A0

def Stark2bands1photon_EnergyShift_notcorrected_numerical(Efield_AU, omega_AU, E_gap_AU, DME_AU=1, x_steps=100):
    # Solver parameters
    x_min = -1       #-1E10
    x_max = 1       #1E10
    #x_steps = 30   #70
    
    tol_merge = 1E-10 #1E3
    merge_treshold = 1 #4? 1: helps to not miss some modes
    roots = []
    branches = []
    unique = []
    start = time()
    for x in np.linspace(x_min, x_max, num=x_steps):
        nrt = root(Stark2bands1photon_Cropped_EnergyShift_polynom_numerical, (x), args=(Efield_AU, omega_AU, E_gap_AU, DME_AU), method='hybr')
        #print nrt
        #checking = func(nrt.x, eps1, eps2, eps3, k0, t, sgn1, sgn2)
        if (nrt.success):
            roots.append(nrt.x)
        
    return np.unique(roots)

### Solves the characteristic polynom for the eigen values of band gap modification    
def Stark4bands1photon_EnergyShift_polynom_numerical(eigen, Efield_AU, omega_AU, E_gap_AU, DME_AU=1):
    M         = DME_AU
    Mbar      = np.conj(DME_AU)
    
    
    A12 = 1.
    
    A11 = -1./2.*E_gap_AU
    
    A10 = -6.*Efield_AU**2*M*Mbar-16.*omega_AU**2-5./4.*E_gap_AU**2
    
    A9 = 23./8.*Efield_AU**2*Mbar*M*E_gap_AU+8.*E_gap_AU*omega_AU**2+5./8.*E_gap_AU**3
    
    A8 = 12.*omega_AU**2*E_gap_AU**2 - 1./4.*Efield_AU**2*Mbar*M*omega_AU*E_gap_AU + 5./8.*E_gap_AU**4 - 8.*Efield_AU**2*M**2*omega_AU**2 - 0.5*Efield_AU**4 * Mbar**3*M + 17./2.*Efield_AU**4*Mbar**2*M**2 - 0.5*Efield_AU**4*Mbar*M**3 - 8.*Efield_AU**2*Mbar**2*omega_AU**2 + 83./16.*Efield_AU**2*Mbar*M*E_gap_AU**2 + 0.25*Efield_AU**2*Mbar**2*E_gap_AU*omega_AU+40.*Efield_AU**2*M*Mbar*omega_AU**2
    
    A7 = -6.*E_gap_AU**3*omega_AU**2-17.*E_gap_AU*omega_AU**2*Efield_AU**2*M*Mbar + 0.25*Efield_AU**2*M*Mbar*E_gap_AU**2*omega_AU - 5./16.*E_gap_AU**5 + 3.*E_gap_AU*Efield_AU**2*Mbar**2*omega_AU**2 - 5./2.*Efield_AU**2*M*Mbar*E_gap_AU**3-0.25*Efield_AU**2*Mbar**2*E_gap_AU**2*omega_AU+3./16.*Efield_AU**4*M*Mbar**3*E_gap_AU-59./16.*Efield_AU**4*M**2*Mbar**2.*E_gap_AU+7./16.*Efield_AU**4*M**3*Mbar*E_gap_AU+4.*E_gap_AU*Efield_AU**2*M**2*omega_AU**2
    
    A6 = -16.*Efield_AU**2*M*Mbar*omega_AU**2*E_gap_AU**2-3*omega_AU**2*E_gap_AU**4+5./32.*Efield_AU**4*M**3*Mbar*E_gap_AU**2-133./32.*Efield_AU**4*M**2*Mbar**2*E_gap_AU**2-5./32.*E_gap_AU**6+3./4.*Efield_AU**4*M**2*Mbar**2*E_gap_AU*omega_AU+1./8.*Efield_AU**2*M*Mbar*E_gap_AU**3*omega_AU-1./2.*Efield_AU**4*Mbar**3*M*E_gap_AU*omega_AU+5*E_gap_AU**2*Efield_AU**2*Mbar**2*omega_AU**2-1./8.*Efield_AU**2*Mbar**2*E_gap_AU**3*omega_AU-5./2.*Efield_AU**6*Mbar**3*M**3-1./2.*Efield_AU**6*Mbar**4*M**2-1./2.*Efield_AU**6*Mbar**2*M**4+13./32.*Efield_AU**4*Mbar**3*M*E_gap_AU**2-13./8.*Efield_AU**2*M*Mbar*E_gap_AU**4-20.*Efield_AU**4*M**2*Mbar**2*omega_AU**2+4*omega_AU**2*Efield_AU**4*M**3*Mbar+4.*omega_AU**2*Efield_AU**4*Mbar**3*M-1./4.*E_gap_AU*Efield_AU**4*M**3*Mbar*omega_AU+4.*E_gap_AU**2*Efield_AU**2*M**2*omega_AU**2
    
    A5 = 3./2.*E_gap_AU**5*omega_AU**2+3./32.*Efield_AU**6*M**4*Mbar**2*E_gap_AU + 51./64.*Efield_AU**2*M*Mbar*E_gap_AU**5-13./64.*Efield_AU**4*M**3*Mbar*E_gap_AU**3-9./64.*Efield_AU**4*Mbar**3*M*E_gap_AU**3+5./64.*E_gap_AU**7-3./2.*E_gap_AU*omega_AU**2*Efield_AU**4*M**3*Mbar+3./8.*Efield_AU**4*M**3*Mbar*omega_AU*E_gap_AU**2-3./16.*Efield_AU**2*M*Mbar*E_gap_AU**4*omega_AU+3./8.*Efield_AU**4*M*Mbar**3*E_gap_AU**2*omega_AU+27./4.*E_gap_AU**3*Efield_AU**2*M*Mbar*omega_AU**2+13./2*E_gap_AU*omega_AU**2*Efield_AU**4*M**2*Mbar**2-E_gap_AU*omega_AU**2*Efield_AU**4*Mbar**3*M-3./4.*Efield_AU**4*M**2*Mbar**2*E_gap_AU**2*omega_AU+3./16.*Efield_AU**2*Mbar**2*E_gap_AU**4*omega_AU-7./4.*E_gap_AU**3*Efield_AU**2*Mbar**2*omega_AU**2+7./32.*Efield_AU**6*Mbar**4*M**2*E_gap_AU+117./64.*Efield_AU**4*M**2*Mbar**2*E_gap_AU**3+27./32.*Efield_AU**6*M**3*Mbar**3*E_gap_AU-2*E_gap_AU**3*Efield_AU**2*M**2*omega_AU**2
    
    A4 = 1./4.*omega_AU**2*E_gap_AU**6-3./2.*E_gap_AU**2*Efield_AU**4*M**3*Mbar*omega_AU**2-3./2.*E_gap_AU**2*Efield_AU**4*M*Mbar**3*omega_AU**2+3/2*Efield_AU**2*M*Mbar*omega_AU**2*E_gap_AU**4+4*Efield_AU**4*M**2*Mbar**2*omega_AU**2*E_gap_AU**2+5./256.*E_gap_AU**8+1./16.*Efield_AU**6*Mbar**2*M**4*E_gap_AU*omega_AU+2*Efield_AU**6*M**3*Mbar**3*omega_AU**2+1./8.*Efield_AU**8*Mbar**5*M**3+3./16.*Efield_AU**8*Mbar**4*M**4+1./16.*Efield_AU**8*Mbar**6*M**2+1./8.*Efield_AU**8*Mbar**3*M**5+1./16.*Efield_AU**8*Mbar**2*M**6+27./128.*Efield_AU**2*M*Mbar*E_gap_AU**6+55./64*Efield_AU**6*M**3*Mbar**3*E_gap_AU**2+7./64.*Efield_AU**6*Mbar**4*M**2*E_gap_AU**2-1./8.*E_gap_AU*Efield_AU**6*M**3*Mbar**3*omega_AU+1./16.*Efield_AU**6*Mbar**5*M*E_gap_AU*omega_AU+75./128.*Efield_AU**4*M**2*Mbar**2*E_gap_AU**4+1./128.*Efield_AU**4*M**3*Mbar*E_gap_AU**4-15./128.*Efield_AU**4*Mbar**3*M*E_gap_AU**4-1./16.*Efield_AU**4*M**2*Mbar**2*E_gap_AU**3*omega_AU-1./8.*Efield_AU**4*M**3*Mbar*omega_AU*E_gap_AU**3+3./16.*Efield_AU**4*M*Mbar**3*E_gap_AU**3*omega_AU+7./64.*Efield_AU**6*M**4*Mbar**2*E_gap_AU**2-Efield_AU**2*Mbar**2*E_gap_AU**4*omega_AU**2-1./2.*Efield_AU**2*M**2*E_gap_AU**4*omega_AU**2
    
    A3 = -3./64.*omega_AU*Efield_AU**2*Mbar**2*E_gap_AU**6+5./16.*Efield_AU**2*Mbar**2*E_gap_AU**5*omega_AU**2+1./4.*Efield_AU**2*M**2*omega_AU**2*E_gap_AU**5-11./16.*Efield_AU**2*M*Mbar*E_gap_AU**5*omega_AU**2-1./64.*Efield_AU**8*M**2*Mbar**6*E_gap_AU-7./64.*Efield_AU**2*M*Mbar*E_gap_AU**7+9./256.*Efield_AU**4*M*Mbar**3*E_gap_AU**5-73./256.*Efield_AU**4*M**2*Mbar**2*E_gap_AU**5-15./64.*Efield_AU**6*M**3*Mbar**3*E_gap_AU**3-3./64.*Efield_AU**8*M**4*Mbar**4*E_gap_AU+1./64.*Efield_AU**6*M**4*Mbar**2*E_gap_AU**3-1./32.*Efield_AU**8*M**3*Mbar**5*E_gap_AU-3./64.*Efield_AU**6*M**2*Mbar**4*E_gap_AU**3-1./32.*Efield_AU**8*M**5*Mbar**3*E_gap_AU-1./64.*Efield_AU**8*M**6*Mbar**2*E_gap_AU+5./256.*Efield_AU**4*M**3*Mbar*E_gap_AU**5+1./4.*omega_AU*Efield_AU**4*M**2*Mbar**2*E_gap_AU**4-1./16.*omega_AU*Efield_AU**4*M**3*Mbar*E_gap_AU**4-3./16.*omega_AU*Efield_AU**4*Mbar**3*M*E_gap_AU**4-3./4.*Efield_AU**4*M**2*Mbar**2*E_gap_AU**3*omega_AU**2+1./2.*Efield_AU**4*M**3*Mbar*omega_AU**2*E_gap_AU**3+1./4.*Efield_AU**4*M*Mbar**3*E_gap_AU**3*omega_AU**2+1./16.*Efield_AU**6*M**3*Mbar**3*E_gap_AU**2*omega_AU+3./64.*Efield_AU**2*M*Mbar*E_gap_AU**6*omega_AU-1./16.*Efield_AU**6*M**4*Mbar**2*omega_AU*E_gap_AU**2-1./2.*E_gap_AU*omega_AU**2*Efield_AU**6*M**3*Mbar**3-1./8.*E_gap_AU**7*omega_AU**2-5./512.*E_gap_AU**9
    
    A2 = -1./1024.*E_gap_AU**10+1./128.*Efield_AU**6*M**2*Mbar**4*E_gap_AU**4-1./128.*Efield_AU**2*M*Mbar*E_gap_AU**8-7./128.*Efield_AU**6*M**3*Mbar**3*E_gap_AU**4+1./128.*Efield_AU**6*M**4*Mbar**2*E_gap_AU**4-7./512.*Efield_AU**4*M**2*Mbar**2*E_gap_AU**6-3./64.*E_gap_AU**2*Efield_AU**8*Mbar**3*M**5-9./128.*Efield_AU**8*M**4*Mbar**4*E_gap_AU**2-3./64.*E_gap_AU**2*Efield_AU**8*Mbar**5*M**3-3./128.*E_gap_AU**2*Efield_AU**8*Mbar**2*M**6-1./512.*Efield_AU**4*M**3*Mbar*E_gap_AU**6+1./16.*Efield_AU**2*Mbar**2*E_gap_AU**6*omega_AU**2+1./128.*omega_AU*Efield_AU**2*Mbar**2*E_gap_AU**7+7./512.*Efield_AU**4*M*Mbar**3*E_gap_AU**6-3./128.*E_gap_AU**2*Efield_AU**8*Mbar**6*M**2+3./64.*omega_AU*Efield_AU**4*M**3*Mbar*E_gap_AU**5-3./64.*omega_AU*Efield_AU**4*M**2*Mbar**2*E_gap_AU**5-1./128.*Efield_AU**2*M*Mbar*E_gap_AU**7*omega_AU+1./8.*Efield_AU**4*M*Mbar**3*E_gap_AU**4*omega_AU**2-1./4.*Efield_AU**4*M**2*Mbar**2*E_gap_AU**4*omega_AU**2-1./4.*Efield_AU**6*M**3*Mbar**3*omega_AU**2*E_gap_AU**2+1./32.*Efield_AU**6*M**3*Mbar**3*E_gap_AU**3*omega_AU+1./8.*Efield_AU**4*M**3*Mbar*omega_AU**2*E_gap_AU**4-1./32.*Efield_AU**6*M*Mbar**5*E_gap_AU**3*omega_AU
    
    A1 = 1./64.*Efield_AU**6*M**4*Mbar**2*omega_AU*E_gap_AU**4+1./256.*E_gap_AU**3*Efield_AU**8*Mbar**2*M**6-5./512.*Efield_AU**6*M**4*Mbar**2*E_gap_AU**5+1./128.*E_gap_AU**3*Efield_AU**8*Mbar**5*M**3+1./128.*E_gap_AU**3*Efield_AU**8*Mbar**3*M**5+3./512.*Efield_AU**6*M**3*Mbar**3*E_gap_AU**5+3./256.*Efield_AU**8*M**4*Mbar**4*E_gap_AU**3+15./1024.*Efield_AU**4*M**2*Mbar**2*E_gap_AU**7-1./512.*Efield_AU**6*M**2*Mbar**4*E_gap_AU**5+11./2048.*Efield_AU**2*M*Mbar*E_gap_AU**9-1./64.*Efield_AU**2*Mbar**2*omega_AU**2*E_gap_AU**7+1./256.*Efield_AU**2*Mbar**2*omega_AU*E_gap_AU**8-3./1024.*Efield_AU**4*M*Mbar**3*E_gap_AU**7+1./256.*E_gap_AU**3*Efield_AU**8*Mbar**6*M**2+1./1024.*Efield_AU**4*M**3*Mbar*E_gap_AU**7-1./256.*Efield_AU**2*M*Mbar*E_gap_AU**8*omega_AU-1./128.*omega_AU*Efield_AU**4*M**3*Mbar*E_gap_AU**6-1./64.*Efield_AU**6*M**3*Mbar**3*E_gap_AU**4*omega_AU-1./64.*omega_AU*Efield_AU**4*M**2*Mbar**2*E_gap_AU**6+3./128.*omega_AU*Efield_AU**4*Mbar**3*M*E_gap_AU**6-1./32.*Efield_AU**4*M**3*Mbar*omega_AU**2*E_gap_AU**5+1./64.*Efield_AU**2*M*Mbar*E_gap_AU**7*omega_AU**2+1./32.*Efield_AU**4*M**2*Mbar**2*E_gap_AU**5*omega_AU**2+1./2048.*E_gap_AU**11
    
    A0 = -3./2048.*Efield_AU**4*Mbar**2*M**2*E_gap_AU**8+3./512.*Efield_AU**8*M**4*Mbar**4*E_gap_AU**4+1./512.*Efield_AU**8*M**6*Mbar**2*E_gap_AU**4+1./256.*Efield_AU**8*M**5*Mbar**3*E_gap_AU**4+1./256.*Efield_AU**8*M**3*Mbar**5*E_gap_AU**4+1./512.*E_gap_AU**4*Efield_AU**8*Mbar**6*M**2-1./4096.*Efield_AU**2*M*Mbar*E_gap_AU**10-1./2048.*Efield_AU**4*M*Mbar**3*E_gap_AU**8-1./1024.*Efield_AU**6*M**3*Mbar**3*E_gap_AU**6-1./1024.*Efield_AU**6*M**4*Mbar**2*E_gap_AU**6-1./1024.*Efield_AU**6*M**2*Mbar**4*E_gap_AU**6-1./2048.*Efield_AU**4*M**3*Mbar*E_gap_AU**8-1./1024.*Efield_AU**2*Mbar**2*omega_AU*E_gap_AU**9+1./256.*Efield_AU**6*M*Mbar**5*E_gap_AU**5*omega_AU+1./1024.*Efield_AU**2*M*Mbar*E_gap_AU**9*omega_AU-1./256.*Efield_AU**6*M**4*Mbar**2*omega_AU*E_gap_AU**5+1./256.*omega_AU*Efield_AU**4*M**2*Mbar**2*E_gap_AU**7-1./256.*omega_AU*Efield_AU**4*Mbar**3*M*E_gap_AU**7
    
    return A12 * eigen**12 + A11*eigen**11 + A10*eigen**10 + A9*eigen**9 + A8*eigen**8 + A7*eigen**7 + A6*eigen**6 + A5*eigen**5 + A4 * eigen**4 + A3 * eigen**3 + A2*eigen**2 + A1 * eigen + A0

def Stark4bands1photon_EnergyShift_notcorrected_numerical(Efield_AU, omega_AU, E_gap_AU, DME_AU=1, x_steps=100):
    # Solver parameters
    x_min = -1       #-1E10
    x_max = 1       #1E10
    #x_steps = 30   #70
    
    tol_merge = 1E-10 #1E3
    merge_treshold = 1 #4? 1: helps to not miss some modes
    roots = []
    branches = []
    unique = []
    start = time()
    for x in np.linspace(x_min, x_max, num=x_steps):
        nrt = root(Stark4bands1photon_EnergyShift_polynom_numerical, (x), args=(Efield_AU, omega_AU, E_gap_AU, DME_AU), method='hybr')
        #print nrt
        #checking = func(nrt.x, eps1, eps2, eps3, k0, t, sgn1, sgn2)
        if (nrt.success):
            roots.append(nrt.x)
        
    return np.unique(roots)

# We call a linear algebra library instead of using Filip solver.
def Stark4bands1photon_EnergyShift_eigen(Efield_AU, omega_AU, E_gap_AU, DME_AU=1):
    M         = DME_AU
    Mbar      = np.conj(DME_AU)
    AMover2=Efield_AU*M/2.
    AMover2c=Efield_AU*Mbar/2.
    matrix = np.array(
        [[E_gap_AU/2.-omega_AU, -omega_AU, -omega_AU, -omega_AU, 0, AMover2, AMover2,AMover2,0,0,0,0], 
         [-omega_AU, E_gap_AU/2.-omega_AU,-omega_AU, -omega_AU, AMover2c, 0, AMover2, AMover2, 0, 0, 0, 0], 
         [-omega_AU, -omega_AU, -E_gap_AU/2.-omega_AU, -omega_AU, AMover2c, AMover2c, 0.,AMover2, 0, 0, 0, 0], 
         [-omega_AU, -omega_AU, -0.5*E_gap_AU-omega_AU, -omega_AU, AMover2c, AMover2c, AMover2c, 0, 0, 0, 0, 0], 
         [0, AMover2, AMover2, AMover2, E_gap_AU/2., 0, 0, 0, 0, AMover2, AMover2, AMover2],
         [AMover2c, 0, AMover2, AMover2, 0, 0.5*E_gap_AU, 0, 0, AMover2c, 0, AMover2, AMover2], 
         [AMover2c, AMover2c, 0, AMover2, 0, 0, -0.5*E_gap_AU, 0, AMover2c, AMover2c, 0, AMover2],
         [AMover2c, AMover2c, AMover2c, 0, 0, 0, 0, -0.5*E_gap_AU, AMover2c, AMover2c, AMover2c, 0],
         [0, 0, 0, 0, 0, AMover2, AMover2, AMover2, 0.5*E_gap_AU+omega_AU, omega_AU, omega_AU, omega_AU], 
         [0, 0, 0, 0, AMover2c, 0, AMover2, AMover2, omega_AU, 0.5*E_gap_AU+omega_AU, omega_AU, omega_AU], 
         [0, 0, 0, 0, AMover2c, AMover2c, 0, AMover2, omega_AU, omega_AU, -0.5*E_gap_AU+omega_AU, omega_AU], 
         [0, 0, 0, 0, AMover2c, AMover2c, AMover2c, 0, omega_AU, omega_AU, omega_AU, -0.5*E_gap_AU+omega_AU]])
    
    #print np.size(matrix)
    
    w, v = LA.eig(matrix)
    return w

# We call a linear algebra library instead of using Filip solver.
def Stark2bands2photons_EnergyShift_eigen(Efield_AU, omega_AU, E_gap_AU, DME_AU=1):
    M         = DME_AU
    Mbar      = np.conj(DME_AU)
    AMover2=Efield_AU*M/2.
    AMover2c=Efield_AU*Mbar/2.
    matrix = np.array(
        [[E_gap_AU/2.-2.*omega_AU, -2.*omega_AU, 0, AMover2, 0, 0, 0, 0, 0, 0], 
         [-2.*omega_AU, -E_gap_AU/2.-2.*omega_AU, AMover2c, 0, 0, 0, 0, 0, 0, 0],
         
         [0, AMover2, E_gap_AU/2.-omega_AU, -omega_AU, 0, AMover2, 0, 0, 0, 0], 
         [AMover2c, 0, -omega_AU, -0.5*E_gap_AU-omega_AU, AMover2c, 0, 0, 0, 0, 0], 
         
         [0, 0, 0, AMover2, E_gap_AU/2., 0, 0, AMover2, 0, 0],
         [0, 0, AMover2c, 0, 0, -0.5*E_gap_AU, AMover2c, 0, 0, 0], 
         
         [0, 0, 0, 0, 0, AMover2, 0.5*E_gap_AU+omega_AU, 0, 0, AMover2],
         [0, 0, 0, 0, AMover2c, 0, 0, -0.5*E_gap_AU+omega_AU, AMover2c, 0],
         
         [0, 0, 0, 0, 0, 0, 0, AMover2, 0.5*E_gap_AU+2.*omega_AU, 0], 
         [0, 0, 0, 0, 0, 0, AMover2c, 0, 0, -0.5*E_gap_AU+2*omega_AU]])
    #print np.size(matrix)
    
    w, v = LA.eig(matrix)
    return w

def TestingNumericalSolver4():
    samples = 10
    DME = 1. #-1+2j #arbitrary!
    Efield_AU = np.linspace(0,1, samples)

    #E_gap_SI      = 2.56*e
    #wavelength    = 800e-9
    #omega_SI      = 2.*np.pi * c / wavelength

    E_gap_AU  = 1 #Energy_eV_to_Hartree(E_gap_SI/e)
    omega_AU  = E_gap_AU #Energy_eV_to_Hartree(omega_SI*hbar/e)

    print("Attempting numerical solution ...")
    roots=[]
    #Eexact=np.zeros((samples,20))
    #index=0
    for index in np.arange(0,len(Efield_AU)):
        Eexact= Stark2bands1photon_EnergyShift_notcorrected_numerical(Efield_AU[index], omega_AU, E_gap_AU, DME, 20)
        roots.append(np.unique(Eexact)) #removes numerical degeneracies

    print np.shape(Efield_AU)
    print np.shape(roots)
    print roots[0][:]
    return roots




#exit()

#E1, E2, E3, E4, E5, E6 = Stark2bands1photon_EnergyShift_modified_exact(Efield_AU, omega_AU, E_gap_AU, DME)
#ENC1, ENC2, ENC3, ENC4, EgapShift_AU = Stark2bands1photon_Cropped_EnergyShift_notcorrected_exact(Efield_AU, omega_AU, E_gap_AU, DME)

#plt.figure()
##plt.scatter(Efield_AU, roots, 'g+', label="6x6 (original)")
#plt.plot(Efield_AU, E1, 'k-', label="6x6 (deduced)")
#plt.plot(Efield_AU, E2, 'k-')
#plt.plot(Efield_AU, E3, 'k-')
#plt.plot(Efield_AU, E4, 'k-')
#plt.plot(Efield_AU, E5, 'k-')
#plt.plot(Efield_AU, E6, 'k-')
#plt.plot(Efield_AU, ENC1, 'b--', label="4x4 (original)")
#plt.plot(Efield_AU, ENC2, 'r--')
#plt.plot(Efield_AU, ENC3, 'b--')
#plt.plot(Efield_AU, ENC4, 'r--')
##plt.plot(Efield_AU, EgapShift_AU, 'g--', label=r"$E_g$ (Stark)")

#plt.xlabel("Field amplitude (at.u.)")
#plt.ylabel("Energy shifts (Hartree)")
#plt.legend(loc="best")
#plt.show()
