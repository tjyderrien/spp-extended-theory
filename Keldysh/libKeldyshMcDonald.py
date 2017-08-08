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

## @package libKeldysh-McDonald
# This module aims to calculate the density of excited electrons as function of laser parameters. 
# Several flavors of the Keldysh theory are available: 
# - Keldysh original paper in solid, for Kane band structure [compared with td-dft]
# - Keldysh paper with few terms corrected by Gruzdev [compared with td-dft]
# - Keldysh-Zhukov tables, where Keldysh theory was computed numerically without using the saddle point method
# - Keldysh-Shcheblanov model, improving rigor on the analytical integration [https://arxiv.org/abs/1706.07303]
# - Keldysh-McDonald-Corkum model, allowing for analytical treatment of mulltiwavelength fields [Physical Review Letters, 2017, 118, 173601]

# IMPORT LIBRARIES
import numpy as np
import numpy.linalg as npl
from numpy import genfromtxt, loadtxt, chararray
#from scipy.optimize import fsolve, root
from scipy.special import ellipk, ellipe, dawsn, factorial2, factorial, ellipkm1
#import cmath
import matplotlib as mp
import matplotlib.pyplot as plt
from scipy.interpolate import interp2d, InterpolatedUnivariateSpline
from matplotlib import rc
# from pylab import *
from scipy.constants import c, epsilon_0, mu_0, pi, e, m_e, h, hbar
#from matplotlib.legend_handler import HandlerLine2D
#import sys

from libUnits import *       #part of the spp-extended-theory/Keldysh
from libAtomicUnits import * #part of the octopus-slabs repository. 
from libDatabase import *    #part of the spp-extended-theory
from libKeldyshPulses import * #part ot the spp-etended-theory/Keldysh

Header="[libKeldyshMcDonald] "

## Builds a conduction band model epsilon(k) as function of longitunal k
# mesh1D k -> mesh1D epsilon (Bohr)
# @param A_projected_AU: A projected 
# @param EgapDirect: value of the direct band gap energy (eV)
# @param kpoints: number of k-points in 1D direction
# @param AtomicDistance: distance between lattice sites (in Bohr)
def SiBandStructureLongitudinal(A_projected_AU, EgapDirect=2.65, kpoints=64, AtomicDistance=5.32): #{{{
  print Header+"Defining the band structure."
  #AtomicDistance=5.32 #interatomic distance in Si (at. u.) 

  #kN = 64 #resolution of k space

  # Building the objects
  epsilonLong = [ Energy_eV_to_Hartree(EgapDirect) ] #LDA band gap of Si #TODO: this must be field depedent to account for stark effect

  Delta = np.max(epsilonLong)
  alpha = np.array([Delta/2., -Delta/2.]) #|alpha> # First, we use the coeffs of the McDonald paper. 
  index = np.arange(0,len(alpha)) # |j>
  #print alpha, index

  # Give values to alpha.
  kmin = 0.
  kmax = 2.*pi/AtomicDistance
  k = np.arange(0., kmax, (kmax-kmin) / kpoints ) #meshing |k> space until 2pi/a
  
  # We shall apply the time dependence of k HERE.
  print np.shape(A_projected_AU)
  k_td = np.broadcast_to(k,(len(A_projected_AU),len(k))) #extending size of k into time
  #print np.shape(k_td)
  A_projected_AU_td = np.broadcast_to(A_projected_AU, (len(k), len(A_projected_AU) ) )
  A_projected_AU_td = np.transpose(A_projected_AU_td)
  #print Header+"Boundaries of A_projected_AU (Hartree/Bohr)."
  print Header+"Boundaries of A_projected_SI (V/m)."
  print Field_AU_to_SI(A_projected_AU.max())
  k_td = k + A_projected_AU_td #shifts k|| by A|| in a time-dependent fashion. 
  
  print Header+"Dimension of k_td="+str(np.shape(k_td)) #GOOD. 
  
  #print index, k
  #print Header+"|j><k|"
  #print np.outer(index, k) #generates |k><index|. How to compute |index><k| ?
  #print Header+"|k><j|"
  #print np.outer(k, index)

  # How to create a new dimension with len(time)? Introduce it before giving k to index. 
  cosTerm = np.cos(np.outer(k_td, index)*AtomicDistance)

  print Header+"cosTerm"
  print np.shape(cosTerm)
  print cosTerm

  termsM1  = cosTerm*alpha #What * is exactly doing here? 

  #print Header+"cos(:,:) * <alpha|"
  ##print termsM1

  ## Indicial writing: very clear, very intuitive, computationally expensive
  #termsL = np.zeros((len(k), len(alpha)))
  #for i in np.arange(0,kN):
    #for j in np.arange(0,len(alpha)):
      #termsL[i,j] = alpha[j]*np.cos(j*k[i]*AtomicDistance)

  #print Header+"** Comparison Matrix vs indice method"
  #print Header+"== VECTOR =="
  #print termsL
  #print Header+"== MATRIX =="
  #print termsM1
  #print Header+"Computing band gap using alpha_j"
  epsilonLong_AU = np.sum(termsM1, 1)
  epsilonLong_SI = Energy_Hartree_to_eV(epsilonLong_AU)

  print Header+"Dimension of k-space..."
  print len(k)
  
  print Header+"Dimension of epsilon_Long_SI"
  print len(epsilonLong_SI)

  plt.figure()
  plt.xlabel('k')
  plt.ylabel(r'$\varepsilon_{||}$')
  plt.plot(k,epsilonLong_SI)
  plt.savefig("BandStructure1D.eps")
  plt.savefig("BandStructure1D.png")
  #plt.show()
          
  return epsilonLong_AU
#}}}



## Build the dipolar transition matrix
# @param Egap: k-dependent band gap energy (at.u.)
# McDonald et al, Phys Rev A 92, 033845 (2015)
# Returns a tuple of 4th order. 
def DipolarTransition(Eg):
  Epx = 0.302 
  Epy = 0.302
  Epz = 0.375
  dx = np.sqrt(Epx / (2e0*Eg**2))
  dy = np.sqrt(Epy / (2e0*Eg**2))
  dz = np.sqrt(Epz / (2e0*Eg**2))
  return np.array([dx, dy, dz])

# Laser pulse is already generated via keldysh.py yes? 
wavelength = 800e-9; PolarizationAngle = 0. 
tau=10e-15; dt = 1E-17; CEP=0e0
PeakFluence = 1.*1E4 #J/cm2 * 1E4 = J/m2
PeakField   = np.sqrt(2e0 * PeakFluence / (tau * c * epsilon_0))

t0=0. #defines the instant 0.
tmin=-1.*tau + t0; tmax=1.*tau + t0

instants = np.arange(tmin, tmax, dt)
#print "Time range: "+str(instants.min())+", "+str(instants.max())+"."

print Header+"** Test: building single pulse centered on 0..."
RealField1x, RealField1y, RealField1z = PulseSquaredSinTemporalShape_vectorial_linear(instants, tau, PeakField, wavelength, PolarizationAngle, CEP, t0, 0.)

## Converts electric field to vector potential |Ax,Ay,Az>(t)
# @param ElectricField in Hartree atomic units. 
# @param instants contains the instants carrying the pulse electric field (any temporal interval is allowed)
def LaserFieldAU_to_VectorPotential(ElectricField, instants):
  dA = ElectricField[1:] * np.diff(instants)
  A = - dA.cumsum() #computes the integral of int(E.dt)
  return A #is A(t) in Hartree atomic units
   
# Converts the electric field (SI) to Atomic Units first. 
ElectricField_x_AU = Field_SI_to_AU(RealField1x)
ElectricField_y_AU = Field_SI_to_AU(RealField1y)
ElectricField_z_AU = Field_SI_to_AU(RealField1z)

# Each component should be independent with time derivation. 
Ax_AU = LaserFieldAU_to_VectorPotential(ElectricField_x_AU, instants)
Ay_AU = LaserFieldAU_to_VectorPotential(ElectricField_y_AU, instants)
Az_AU = LaserFieldAU_to_VectorPotential(ElectricField_z_AU, instants)

A_AU = np.array([Ax_AU, Ay_AU, Az_AU ])
A_norm = np.sqrt(Ax_AU.max()**2+Ay_AU.max()**2+Az_AU.max()**2)
print Header+"Building the vector potential..."
print np.shape(A_AU)
#print A_AU

print Header+"Build the vectorial electric field..."
ElectricField_AU = np.array([ ElectricField_x_AU, ElectricField_y_AU, ElectricField_z_AU ])

# Ok, first we have to construct time-dependent dipolar momentum.
# Question: should we change the way we construct epsilon(k) to epsilon(K+A(t)). 
# But then, we have a time-dependent band gap in one dimension? 

# We build a 3d-TD band structure (TD k-space)
EpsilonX = SiBandStructureLongitudinal(Ax_AU) #A_|| should be here
exit()
EpsilonY = SiBandStructureLongitudinal(Ay_AU) #A_perp should be here
EpsilonZ = SiBandStructureLongitudinal(Az_AU) #A_perp should be here
print np.shape(EpsilonX), np.shape(EpsilonY), np.shape(EpsilonZ)

Epsilon = np.einsum('i,j,k->ijk', EpsilonX, EpsilonY, EpsilonZ)
print np.shape(Epsilon)
print Epsilon

#plt.matshow(Epsilon[32,:,:]) #apparently band structure is well defined
#plt.show()

d = DipolarTransition(Epsilon)
print np.shape(d)

## Have to compute the inner product <d|F(t)>
print Header+"Computing the dipolar coupling strength..."
print np.shape(d), np.shape(ElectricField_AU)
#CouplingStrengh = np.inner(d, ElectricField_AU)
AtomicDistance_AU = 5.32 #Interatomic distances in Bohrs for Silicon
print A_norm, 2.*pi/AtomicDistance_AU