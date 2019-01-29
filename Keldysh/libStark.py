#!/usr/bin/env python2.7
#-*- coding: utf-8 -*-
## @package libStark
## Computes the Stark effect on a band structure at gamma points
# This module aims at computing the band gap energy as function of the average laser fieled induced by the Stark effect

import numpy as np
import matplotlib.pyplot as plt
from scipy.constants import c, epsilon_0, mu_0, pi, e, m_e, h, hbar

import math, cmath
from scipy.optimize import root
from time import time

from libAtomicUnits import *

## Provides the change of the branch electronic levels
# From simple Floquet Hamiltonian on constant pulse of frequency omega, the shift of 6 bands with the electric field is given. The eigen values have been computed from the Hamiltonian given in the Nano Letters. 
def Stark6bandsEnergyShift_modified(Efield_AU, omega_AU, E_gap_AU, DME_AU=1): #{{{
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
def Stark6bandsEnergyShift_polynom(eigen, Efield_AU, omega_AU, E_gap_AU, DME_AU=1):
    M         = DME_AU
    Mbar      = np.conj(DME_AU)
    
    A6 = 1
    
    A4 = -Efield_AU**2*M*Mbar-0.75*E_gap_AU**2-3*omega_AU**2
    
    A3 = 2*omega_AU**3
    
    A2 = 0.5*Efield_AU**2*Mbar*M*E_gap_AU**2 + 0.5*omega_AU**2*E_gap_AU**2 + 0.5*Efield_AU**4+Mbar**2*M**2+3./16.*E_gap_AU**4+1.5*Efield_AU**2*M*Mbar*omega_AU**2
    
    A1 = -0.5*Efield_AU**2+M*Mbar*omega_AU**3-0.5*omega_AU**3*E_gap_AU**2
    
    A0 = -1./16.*omega_AU**2*Efield_AU**4*Mbar**2*M**2
    -1./16.*Efield_AU**2*M*Mbar*E_gap_AU**4+1./8.*Efield_AU**2*M*Mbar*omega_AU**2*E_gap_AU**2 - 1./64.*E_gap_AU**6 - 1./16. * E_gap_AU**2*Efield_AU**4*Mbar**2*M**2+1./16.*omega_AU**2*E_gap_AU**4
    
    return A6 * eigen**6 + A4 * eigen**4 + A3 * eigen**3 + A2*eigen**2 + A1 * eigen + A0

def Stark6bandsEnergyShift_notcorrected(Efield_AU, omega_AU, E_gap_AU, DME_AU=1, x_steps=30):
    # Solver parameters
    x_min = -1E10       #-1E10
    x_max = 1E10       #1E10
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
        nrt = root(Stark6bandsEnergyShift_polynom, (x), args=(Efield_AU, omega_AU, E_gap_AU, DME_AU), method='hybr')
        #print nrt
        #checking = func(nrt.x, eps1, eps2, eps3, k0, t, sgn1, sgn2)
        if (nrt.success):
            roots.append(nrt.x)
        
        ## Filtering the degeneracies
        #ln = len(roots)
        #while ln > 0:
            #center = roots[0]
            #aux = [center]
            #aux2 = []
            #for rt in roots[1:]:
                #if norm(rt-center) < tol_merge:
                    #aux.append(rt)
                    #center = cntr(aux)
                #elif norm(rt+center) < tol_merge:
                    #aux.append(-rt)
                    #center = cntr(aux)
                #else:
                    #aux2.append(rt)
                    
            #if len(aux) > merge_treshold: #merging criterion
                #if center[0] > 0:
                    #unique.append(center)
                #else:
                    #unique.append([-center[0], -center[1]])
            #roots = list(aux2)
            #ln = len(roots)
            
        #branches.append(unique)

    #outs = branches
    
    #selection of maximal Lspps
    #for branch in branches:
        #brinv = [[rt[0]] for rt in branch]
        #outs.append(sorted(brinv, key=lambda x:abs(x[1]), reverse=True)[:num_of_maxs])
        ##outs.append(sorted(brinv, key=lambda x:x[0], reverse=True)[:num_of_maxs])
    #return branches
    #single_roots = set(roots)
    return np.unique(roots)
    
def Stark4bandsEnergyShift_notcorrected(Efield_AU, omega_AU, E_gap_AU, DME_AU=1): #{{{
    M         = DME_AU
    Mbar      = np.conj(DME_AU)
    
    Tmp = omega_AU**2+Efield_AU**2*M*Mbar+E_gap_AU**2
    E1  = 0.5*omega_AU + 0.5*np.sqrt(Tmp+2*E_gap_AU*omega_AU)
    E2  = 0.5*omega_AU - 0.5*np.sqrt(Tmp+2*E_gap_AU*omega_AU)
    E3  = 0.5*omega_AU + 0.5*np.sqrt(Tmp-2*E_gap_AU*omega_AU)
        
    E4  = 0.5*omega_AU - 0.5*np.sqrt(Tmp-2*E_gap_AU*omega_AU)
    
    return E1, E2, E3, E4
#}}}   


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
    Eexact= Stark6bandsEnergyShift_notcorrected(Efield_AU[index], omega_AU, E_gap_AU, DME, 20)
    roots.append(np.unique(Eexact)) #removes numerical degeneracies

print np.shape(Efield_AU)
print np.shape(roots)
print roots[0][:]

#exit()

E1, E2, E3, E4, E5, E6 = Stark6bandsEnergyShift_modified(Efield_AU, omega_AU, E_gap_AU, DME)
ENC1, ENC2, ENC3, ENC4 = Stark4bandsEnergyShift_notcorrected(Efield_AU, omega_AU, E_gap_AU, DME)

plt.figure()
#plt.scatter(Efield_AU, roots, 'g+', label="6x6 (original)")
plt.plot(Efield_AU, E1, 'k-', label="6x6 (deduced)")
plt.plot(Efield_AU, E2, 'k-')
plt.plot(Efield_AU, E3, 'k-')
plt.plot(Efield_AU, E4, 'k-')
plt.plot(Efield_AU, E5, 'k-')
plt.plot(Efield_AU, E6, 'k-')
plt.plot(Efield_AU, ENC1, 'b--', label="4x4 (original)")
plt.plot(Efield_AU, ENC2, 'b--')
plt.plot(Efield_AU, ENC3, 'b--')
plt.plot(Efield_AU, ENC4, 'b--')

plt.xlabel("Field amplitude (at.u.)")
plt.ylabel("Energy shifts (Hartree)")
plt.legend(loc="best")
plt.show()
