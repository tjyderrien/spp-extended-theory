#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Copyright (C) 2016-2021 T. J.-Y. Derrien, K. Gazdova
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
# along with this program. If not, see <htt

## @package libStark
## Computes the Stark effect for a energy band structure at the gamma points
# This module aims at computing the band gap energy as function of the average laser field induced by the Stark effect

from time import time

import numpy as np
import numpy.linalg as LA
from scipy.constants import c, pi, e, h, hbar, h
from scipy.optimize import root
from scipy.special import jn_zeros as BesselJzeros
import octopus_slabs.Libs.libAtomicUnits as au

# import logging
from octopus_slabs.Libs.libLogging import init_logger

#Enable_A2=True

logger = init_logger(__name__, verbose=False) #"plotFinalQuantities")
# logger.basicConfig(format='%(asctime)s %(message)s', datefmt='%m/%d/%Y %I:%M:%S %p')
# logger.basicConfig(filename='libStark.log', level=logger.INFO)

## Provides the shift of the quasi electronic levels
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
    
    A8 = 12.*omega_AU**2*E_gap_AU**2 - 1./4.*Efield_AU**2*Mbar*M*omega_AU*E_gap_AU + 5./8.*E_gap_AU**4 \
         - 8.*Efield_AU**2*M**2*omega_AU**2 - 0.5*Efield_AU**4 * Mbar**3*M + 17./2.*Efield_AU**4*Mbar**2*M**2 \
         - 0.5*Efield_AU**4*Mbar*M**3 - 8.*Efield_AU**2*Mbar**2*omega_AU**2 + 83./16.*Efield_AU**2*Mbar*M*E_gap_AU**2 \
         + 0.25*Efield_AU**2*Mbar**2*E_gap_AU*omega_AU+40.*Efield_AU**2*M*Mbar*omega_AU**2
    
    A7 = -6.*E_gap_AU**3*omega_AU**2-17.*E_gap_AU*omega_AU**2*Efield_AU**2*M*Mbar \
         + 0.25*Efield_AU**2*M*Mbar*E_gap_AU**2*omega_AU - 5./16.*E_gap_AU**5 \
         + 3.*E_gap_AU*Efield_AU**2*Mbar**2*omega_AU**2 - 5./2.*Efield_AU**2*M*Mbar*E_gap_AU**3\
         -0.25*Efield_AU**2*Mbar**2*E_gap_AU**2*omega_AU+3./16.*Efield_AU**4*M*Mbar**3*E_gap_AU\
         -59./16.*Efield_AU**4*M**2*Mbar**2.*E_gap_AU+7./16.*Efield_AU**4*M**3*Mbar*E_gap_AU\
         +4.*E_gap_AU*Efield_AU**2*M**2*omega_AU**2
    
    A6 = -16.*Efield_AU**2*M*Mbar*omega_AU**2*E_gap_AU**2-3*omega_AU**2*E_gap_AU**4+5./32.*Efield_AU**4*M**3*Mbar*E_gap_AU**2-133./32.*Efield_AU**4*M**2*Mbar**2*E_gap_AU**2-5./32.*E_gap_AU**6+3./4.*Efield_AU**4*M**2*Mbar**2*E_gap_AU*omega_AU+1./8.*Efield_AU**2*M*Mbar*E_gap_AU**3*omega_AU-1./2.*Efield_AU**4*Mbar**3*M*E_gap_AU*omega_AU+5*E_gap_AU**2*Efield_AU**2*Mbar**2*omega_AU**2-1./8.*Efield_AU**2*Mbar**2*E_gap_AU**3*omega_AU-5./2.*Efield_AU**6*Mbar**3*M**3-1./2.*Efield_AU**6*Mbar**4*M**2-1./2.*Efield_AU**6*Mbar**2*M**4+13./32.*Efield_AU**4*Mbar**3*M*E_gap_AU**2-13./8.*Efield_AU**2*M*Mbar*E_gap_AU**4-20.*Efield_AU**4*M**2*Mbar**2*omega_AU**2+4*omega_AU**2*Efield_AU**4*M**3*Mbar+4.*omega_AU**2*Efield_AU**4*Mbar**3*M-1./4.*E_gap_AU*Efield_AU**4*M**3*Mbar*omega_AU+4.*E_gap_AU**2*Efield_AU**2*M**2*omega_AU**2
    
    A5 = 3./2.*E_gap_AU**5*omega_AU**2+3./32.*Efield_AU**6*M**4*Mbar**2*E_gap_AU\
         + 51./64.*Efield_AU**2*M*Mbar*E_gap_AU**5-13./64.*Efield_AU**4*M**3*Mbar*E_gap_AU**3\
         -9./64.*Efield_AU**4*Mbar**3*M*E_gap_AU**3+5./64.*E_gap_AU**7\
         -3./2.*E_gap_AU*omega_AU**2*Efield_AU**4*M**3*Mbar+3./8.*Efield_AU**4*M**3*Mbar*omega_AU*E_gap_AU**2\
         -3./16.*Efield_AU**2*M*Mbar*E_gap_AU**4*omega_AU+3./8.*Efield_AU**4*M*Mbar**3*E_gap_AU**2*omega_AU\
         +27./4.*E_gap_AU**3*Efield_AU**2*M*Mbar*omega_AU**2+13./2*E_gap_AU*omega_AU**2*Efield_AU**4*M**2*Mbar**2\
         -E_gap_AU*omega_AU**2*Efield_AU**4*Mbar**3*M\
         -3./4.*Efield_AU**4*M**2*Mbar**2*E_gap_AU**2*omega_AU+3./16.*Efield_AU**2*Mbar**2*E_gap_AU**4*omega_AU\
         -7./4.*E_gap_AU**3*Efield_AU**2*Mbar**2*omega_AU**2\
         +7./32.*Efield_AU**6*Mbar**4*M**2*E_gap_AU+117./64.*Efield_AU**4*M**2*Mbar**2*E_gap_AU**3\
         +27./32.*Efield_AU**6*M**3*Mbar**3*E_gap_AU-2*E_gap_AU**3*Efield_AU**2*M**2*omega_AU**2
    
    A4 = 1./4.*omega_AU**2*E_gap_AU**6-3./2.*E_gap_AU**2*Efield_AU**4*M**3*Mbar*omega_AU**2\
         -3./2.*E_gap_AU**2*Efield_AU**4*M*Mbar**3*omega_AU**2+3/2*Efield_AU**2*M*Mbar*omega_AU**2*E_gap_AU**4\
         +4*Efield_AU**4*M**2*Mbar**2*omega_AU**2*E_gap_AU**2+5./256.*E_gap_AU**8\
         +1./16.*Efield_AU**6*Mbar**2*M**4*E_gap_AU*omega_AU+2*Efield_AU**6*M**3*Mbar**3*omega_AU**2\
         +1./8.*Efield_AU**8*Mbar**5*M**3+3./16.*Efield_AU**8*Mbar**4*M**4+1./16.*Efield_AU**8*Mbar**6*M**2\
         +1./8.*Efield_AU**8*Mbar**3*M**5+1./16.*Efield_AU**8*Mbar**2*M**6+27./128.*Efield_AU**2*M*Mbar*E_gap_AU**6\
         +55./64*Efield_AU**6*M**3*Mbar**3*E_gap_AU**2+7./64.*Efield_AU**6*Mbar**4*M**2*E_gap_AU**2\
         -1./8.*E_gap_AU*Efield_AU**6*M**3*Mbar**3*omega_AU+1./16.*Efield_AU**6*Mbar**5*M*E_gap_AU*omega_AU\
         +75./128.*Efield_AU**4*M**2*Mbar**2*E_gap_AU**4+1./128.*Efield_AU**4*M**3*Mbar*E_gap_AU**4\
         -15./128.*Efield_AU**4*Mbar**3*M*E_gap_AU**4-1./16.*Efield_AU**4*M**2*Mbar**2*E_gap_AU**3*omega_AU\
         -1./8.*Efield_AU**4*M**3*Mbar*omega_AU*E_gap_AU**3+3./16.*Efield_AU**4*M*Mbar**3*E_gap_AU**3*omega_AU\
         +7./64.*Efield_AU**6*M**4*Mbar**2*E_gap_AU**2-Efield_AU**2*Mbar**2*E_gap_AU**4*omega_AU**2\
         -1./2.*Efield_AU**2*M**2*E_gap_AU**4*omega_AU**2
    
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
def Stark4bands1photon_EnergyShift_eigen(Efield_AU, omega_AU, E_gap_AU, DME_AU=1, Enable_A2=False):
    M         = DME_AU
    Mbar      = np.conj(DME_AU)
    AMover2=Efield_AU*M/2.
    AMover2c=Efield_AU*Mbar/2.
    if(Enable_A2):
        A2over4 = Efield_AU ** 2 * omega_AU / 4 / pi
    else:
        A2over4 = 0.
    matrix = np.array(
        [[E_gap_AU / 2. - omega_AU + A2over4, -omega_AU + A2over4, -omega_AU + A2over4, -omega_AU + A2over4, 0, AMover2, AMover2, AMover2, A2over4, A2over4, A2over4, A2over4], 
         [-omega_AU + A2over4, E_gap_AU/2.-omega_AU + A2over4,-omega_AU + A2over4, -omega_AU + A2over4, AMover2c, 0, AMover2, AMover2, A2over4, A2over4, A2over4, A2over4], 
         [-omega_AU + A2over4, -omega_AU + A2over4, -E_gap_AU/2.-omega_AU + A2over4, -omega_AU + A2over4, AMover2c, AMover2c, 0.,AMover2, A2over4, A2over4, A2over4, A2over4], 
         [-omega_AU + A2over4, -omega_AU + A2over4, -0.5*E_gap_AU-omega_AU + A2over4, -omega_AU + A2over4, AMover2c, AMover2c, AMover2c, 0, A2over4, A2over4, A2over4, A2over4], 
         [0, AMover2, AMover2, AMover2, E_gap_AU/2. + A2over4, A2over4, A2over4, A2over4, 0, AMover2, AMover2, AMover2],
         [AMover2c, 0, AMover2, AMover2, A2over4, 0.5*E_gap_AU + A2over4, A2over4, A2over4, AMover2c, 0, AMover2, AMover2], 
         [AMover2c, AMover2c, 0, AMover2, A2over4, A2over4, -0.5*E_gap_AU + A2over4, A2over4, AMover2c, AMover2c, 0, AMover2],
         [AMover2c, AMover2c, AMover2c, 0, A2over4, A2over4, A2over4, -0.5*E_gap_AU + A2over4, AMover2c, AMover2c, AMover2c, 0],
         [A2over4, A2over4, A2over4, A2over4, 0, AMover2, AMover2, AMover2, 0.5*E_gap_AU + omega_AU + A2over4, omega_AU + A2over4, omega_AU + A2over4, omega_AU + A2over4], 
         [A2over4, A2over4, A2over4, A2over4, AMover2c, 0, AMover2, AMover2, omega_AU + A2over4, 0.5*E_gap_AU + omega_AU + A2over4, omega_AU + A2over4, omega_AU + A2over4], 
         [A2over4, A2over4, A2over4, A2over4, AMover2c, AMover2c, 0, AMover2, omega_AU + A2over4, omega_AU + A2over4, -0.5*E_gap_AU+omega_AU + A2over4, omega_AU + A2over4], 
         [A2over4, A2over4, A2over4, A2over4, AMover2c, AMover2c, AMover2c, 0, omega_AU + A2over4, omega_AU + A2over4, omega_AU + A2over4, -0.5*E_gap_AU+omega_AU + A2over4]])
    # TODO: validation by - determinant of diagonalized matrix should be 0
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


# Floquet two-band model with scalar A2 terms included. Terms calculated by K. Gazdova.
def Stark2bands1photons_EnergyShift_eigen_A2(Efield_AU, omega_AU, E_gap_AU, DME_AU=1):
    M = DME_AU
    Mbar = np.conj(DME_AU)
    AMover2 = Efield_AU * M / 2.
    AMover2c = Efield_AU * Mbar / 2.
    A2over4 = Efield_AU ** 2 * omega_AU / 4 / pi
    matrix = np.array(
        [[E_gap_AU / 2. - omega_AU + A2over4, -omega_AU + A2over4, 0, AMover2, Efield_AU**2/4., Efield_AU**2/4.],
         [- omega_AU + A2over4, -E_gap_AU / 2. - omega_AU + A2over4, AMover2c, 0, Efield_AU**2/4., Efield_AU**2/4.],

         [0, AMover2, E_gap_AU / 2. + A2over4, A2over4, 0, AMover2],
         [AMover2c, 0, A2over4, -0.5 * E_gap_AU + A2over4, AMover2c, 0],

         [Efield_AU**2/4., Efield_AU**2/4., 0, AMover2, 0.5 * E_gap_AU + omega_AU + A2over4, omega_AU + A2over4],
         [Efield_AU**2/4., Efield_AU**2/4., AMover2c, 0, omega_AU+A2over4, -0.5 * E_gap_AU + omega_AU + A2over4]]
    )
    # print np.size(matrix)

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

    logger.info("Attempting numerical solution ...")
    roots=[]
    #Eexact=np.zeros((samples,20))
    #index=0
    for index in np.arange(0,len(Efield_AU)):
        Eexact= Stark2bands1photon_EnergyShift_notcorrected_numerical(Efield_AU[index], omega_AU, E_gap_AU, DME, 20)
        roots.append(np.unique(Eexact)) #removes numerical degeneracies

    logger.info(np.shape(Efield_AU))
    logger.info(np.shape(roots))
    logger.info(roots[0][:])
    return roots

## Kronecker of two numbers
def Kronecker(i,j):
  sol=0
  if(i==j):
      sol=1
  else:
      sol=0
  return sol


## Returns the Floquet band structure at a given k-point assuming <MPI_number> photons transitions. 
def ComputeFloquetBandStructure(filename, nb_atoms, Z_electrons, kpoints, unocc_states, MPI_number, Efield_SI,
                                wavelength_SI, num_time_steps, eigenvalues, matrixelements,
                                Complex_valued_ME=False, field_polarization_dir=0,
                                spin_occupation=2, verbose=1, Enable_A2=False): #{{{
    Header="[libStark: ComputeFloquetBandStructure(): ] "
    logger.info("Floquet for wavelength: "+str(wavelength_SI*1E9)+" nm, E = "+str(Efield_SI*1E-9)+" V/nm.")
    matrix_side=Z_electrons+unocc_states
    ## 6. Build then diagonalize the Floquet Hamiltonian. 
    # 6.1:Build the eigenvalued matrix
    H_GS=np.zeros((matrix_side, matrix_side))*1j
    if(verbose==1):
        logger.info("** Info: side of the matrix")
        logger.info("         "+str(np.size(H_GS[0])))
        logger.info("Matrix_side: "+str(matrix_side)+" Z_el: "+str(Z_electrons)+" unocc: "+str(unocc_states))
    for i in np.arange(0,np.size(H_GS[0])):
        H_GS[i,i]=eigenvalues[i]
    if(verbose==1):
        logger.debug("Info: matrix of eigenvalues for the selected k-point. Matrix shape is "+str(np.shape(H_GS)))

    # 6.2: Add the dipolar matrix elements to the eigenvalues

    ## Function that returns the light-perturbed GS Hamiltonian for a specific k-point
    def H_perturb(t, omega_AU, Efield_AU, H_GS, matrixelements, Enable_A2=False):
        integral1_sum=0
        integral2_sum=0
        if Enable_A2:
            for t in np.arange(tmin,tmax,dt):
                integral1=Efield_AU**2*np.cos(omega_AU*t)**2*dt
                integral1_sum=np.add(integral1_sum, integral1)
#           for t in np.arange(tmin,tmax,dt):
#                integral2=integral1_sum*dt
#                integral2_sum=np.add(integral2_sum, integral2)
        else:
            integral1_sum=0
#            integral2_sum=0
        return np.add(np.add(H_GS, Efield_AU*np.cos(omega_AU*t)*matrixelements), integral1_sum)
        
    ## Returns the Rabi frequency for each possible dipolar transition (atomic units)
    # WARNING: 1/c error may be found (when employed convention is E=-1/c dA/dt). 
    def RabiMatrix_AU(t, omega_AU, Efield_AU, matrixelements):
#        c_AU = 1./137.
        #Formula: Efield_AU(t) * d / omega_AU
        RabiMatrix_AU = Efield_AU*np.cos(omega_AU*t)*matrixelements/omega_AU 
        return RabiMatrix_AU

    # 6.3: introducing Fourier transform, multiphotonic levels and temporal averaging over one period.

    omega_SI  = 2.*np.pi*c / wavelength_SI #this is in SI...

    Efield_AU = au.Field_SI_to_AU(Efield_SI)
    omega_AU  = au.Energy_eV_to_Hartree(omega_SI*hbar/e)
    
    tmin=0.; tmax=2.*pi/omega_AU #NOTE: changing this induces a shift to higher energies. Find out why. It should be phase dependent (2 pi / omega). 
    dt=(tmax-tmin)/num_time_steps

    # Initialization for t=0
    if(verbose==1):
        logger.info("Info: H(t)=H_GS + dipolar_matrix * Efield(t)")
        logger.info("Info: shape of H_GS: "+str(np.shape(H_GS)))
        logger.info("Info: shape of matrixelements: "+str(np.shape(matrixelements)))
    
    # Only cosmetic!
    H0=H_perturb(tmin, omega_AU, Efield_AU, H_GS, matrixelements, Enable_A2) #this gives max of RabiFreq. Spanning fields will give the rest.
    Rabi0=RabiMatrix_AU(tmin, omega_AU, Efield_AU, matrixelements)
    if(verbose==1):
        logger.info("Info: shape of H(t): "+str(np.shape(H0)))
        logger.info("Info: shape of RabiFreq(t): "+str(np.shape(Rabi0)))

    ## Temporal integration of H_Floquet^{m,n}
    def H_Floquet_mn_func(tmin, tmax, dt, omega_AU, MPI_number, n, Efield_AU, H_GS, matrixelements, Enable_A2=False): #{{{
        H0_init = H_perturb(tmin, omega_AU, Efield_AU, H_GS, matrixelements, Enable_A2) #just to initialize
        H0_sum = np.zeros(np.shape(H0_init))*0j #initialization
        for t in np.arange(tmin,tmax,dt): #integral
            H0=H_perturb(t, omega_AU, Efield_AU, H_GS, matrixelements, Enable_A2) #H(t)
            H0_t1=omega_AU/2./np.pi * np.multiply(np.exp(1j*(MPI_number-n)*omega_AU*t), H0)*dt #e( i(m-n) omega t ) H(t)
            H0_sum = np.add(H0_sum, H0_t1) #have to integrate this
           
        H0_sum = H0_sum + Kronecker(MPI_number,n)*MPI_number*omega_AU
        return H0_sum
    #}}}
    
    ## Rabi frequeqncy is time dependent, and is a large matrix. We mainly need the extrema for now. 
    # These correspond to 
    # d e f Rabi_extrema(tmin, tmax, dt, omega_AU, Efield_AU, matrixelements): #{{{
        #Rabi_init=RabiMatrix_AU(tmin, omega_AU, Efield_AU, matrixelements)
        #Rabi_sum = RabiMatrix_AU(tmax, omega_AU, Efield_AU, matrixelements)     
        #Rabi_t = [] # Rabi_sum = np.zeros(np.shape(Rabi_init))*0j #initialization
        #for t in np.arange(tmin,tmax,dt): #integral
            #Rabi0 = RabiMatrix_AU(t, omega_AU, Efield_AU, matrixelements) #Rabi(t)
            #Rabi_t.append(Rabi0)
        #return np.min(np.min(Rabi_t)), np.max(np.max(Rabi_t))
            

    H_Floquet_mn=H_Floquet_mn_func(tmin, tmax, dt, omega_AU, MPI_number, MPI_number, Efield_AU, H_GS, matrixelements, Enable_A2)
    
    #Rabi0=RabiMatrix_AU(tmin, omega_AU, Efield_AU, matrixelements) #Rabi(t) #This gives the maximum Rabi frequencies already. Spanning in field will give the right set of values. No need to span in time. 
    #RabiFloquet_AU_min, RabiFloquet_AU_max = Rabi_extrema(tmin, tmax, dt, omega_AU, Efield_AU, matrixelements)
    RabiFloquet_AU = Rabi0 #just for name consistency
    
    RabiFloquet_AU_min = np.min(np.min(np.abs(RabiFloquet_AU)))
    RabiFloquet_AU_max = np.max(np.max(np.abs(RabiFloquet_AU)))
    RabiFloquet_SI_min = au.Energy_Hartree_to_eV(RabiFloquet_AU_min)
    RabiFloquet_SI_max = au.Energy_Hartree_to_eV(RabiFloquet_AU_max)
    
    #print RabiFloquet_AU
    #print "Rabi energy (Ha): ["+ str(RabiFloquet_AU_min) +", "+str(RabiFloquet_AU_max)+"]"
    if(verbose==1):
        logger.info("Rabi energy (eV): ["+ str(RabiFloquet_SI_min) +", "+str(RabiFloquet_SI_max)+"]")
    
    ## Print here the ratio Rabi/Laser and the value of the BesselFunction.
    ArgForBesselJ=np.abs(RabiFloquet_AU)/omega_AU
    
    #print "x for BesselJ(x): "+str(np.shape(ArgForBesselJ))
    #print "Rabi/laser (a.u.): ["+str(np.min(ArgForBesselJ))+", "+str(np.max(ArgForBesselJ))+"]"
    #print "BesselJ(Rabi/laser): "+str(BesselJ(0, RabiFloquet_AU_max/omega_AU))
    
    ## Returns the optimal omega_AU for disabling tunneling
    # In principle....
    #d e f OmegaCutoff(Efield_AU, dipole, MPI_number, rootnum=1):
        #return np.sqrt( Efield_AU*np.divide( dipole, BesselJzeros( MPI_number, rootnum ) ) )
    
    #print "Test: BesselJzeros(MPI_number, 1): "+str(BesselJzeros(MPI_number, 1))
    
    #root_orders = np.arange(1,5)
    #OmegaCutoff=np.vectorize(OmegaCutoff)
    OmegaCutOff_AU=[]
    
    for root_orders in np.arange(1,5):
        sol=np.sqrt(Efield_AU*matrixelements/BesselJzeros(MPI_number, root_orders)[0]) #cutoff energies (a.u.)
        #print sol
        OmegaCutOff_AU.append(sol)
    
    OmegaCutOff_eV          = au.Energy_Hartree_to_eV(OmegaCutOff_AU)
    OmegaCutOff_m           = h*c/au.Energy_Hartree_to_eV(OmegaCutOff_AU)/e #WARNING: divide by 0 if empty
    OmegaCutOff_eV_filtered = (list(set(OmegaCutOff_eV[OmegaCutOff_eV>0.1]))) #0.1 eV
    OmegaCutOff_m_filtered  = (list(set(OmegaCutOff_m[OmegaCutOff_m<5e-6]))) #5 um
    OmegaCutOff_eV_sorted   = np.sort(OmegaCutOff_eV_filtered)
    OmegaCutOff_m_sorted    = np.sort(OmegaCutOff_m_filtered)
    
    if(verbose==1):
        logger.info("Ideal photon energies for "+str(Efield_SI*1E-9)+" V/nm: \n"+str(list(set(OmegaCutOff_eV_sorted)))+" (eV)")
        logger.info("Ideal photon wavelengths for "+str(Efield_SI*1E-9)+" V/nm: \n"+str(list(set(OmegaCutOff_m_sorted)))+" (m)")

    #  wavelength (m) = h * c / (E(eV) * e)
    #print H_Floquet_mn 
    
    if(verbose==1):
        logger.info("Info: element (m,n) of Floquet matrix constructed with success.")
        logger.info("Info: Now, one has to write the full matrix.")

    if(verbose==1):
        logger.info("Info: shape of element [m,n] of the matrix we want to write: "+str(np.shape(H_Floquet_mn)))
    H_Floquet_shape = (MPI_number*2+1) * (matrix_side)
    H_Floquet=np.zeros((H_Floquet_shape, H_Floquet_shape))*1j
    if(verbose==1): 
        logger.info("Info: shape of the complete Floquet matrix")
        logger.info(H_Floquet_shape, H_Floquet_shape)
        logger.info(np.shape(H_Floquet))
    for m in np.arange(-MPI_number,MPI_number+1,1):
        for n in np.arange(-MPI_number,MPI_number+1,1):
            H_Floquet_mn=H_Floquet_mn_func(tmin, tmax, dt, omega_AU, m, n, Efield_AU, H_GS, matrixelements, Enable_A2) #tensor of 4th order...
            #print "m="+str(m)
            #print "m+MPI_number="+str(m+MPI_number)
            #print "n="+str(n)
            #print "n+MPI_number="+str(n+MPI_number)
            for p in np.arange(0,matrix_side,1):
                for q in np.arange(0,matrix_side,1):
                    #print p,q #OK
                    #print (m+MPI_number)*matrix_side+p, (n+MPI_number)*matrix_side+q
                    H_Floquet[(m+MPI_number)*matrix_side+p, (n+MPI_number)*matrix_side+q]=H_Floquet_mn[p,q]
    #print H_Floquet

    # 6.4: Diagonalize the Floquet matrix
    replicas, v = LA.eig(H_Floquet)
    #print replicas
    
    return eigenvalues, replicas, matrixelements
#}}}



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
