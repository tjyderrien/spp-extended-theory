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

Header="[libKeldyshMcDonald] "

## Builds the band structure 1D with band structure from McDonald paper
# @param EgapDirect: value of the direct band gap energy (eV)
# @param kpoints: number of k-points in 1D direction
# @param AtomicDistance: distance between lattice sites (in Bohr)
def BuildSiBandStructure(EgapDirect=2.65, kpoints=64, AtomicDistance=5.32):
  print Header+"Defining the band structure."
  #AtomicDistance=5.32 #interatomic distance in Si (at. u.) 

  #kN = 64 #resolution of k space

  # Building the objects
  epsilonLong = [ Energy_eV_to_Hartree(EgapDirect) ] #LDA band gap of Si #TODO: his must be field depedent to account for stark effect

  Delta = np.max(epsilonLong)
  alpha = np.array([Delta/2., -Delta/2.]) #|alpha> # First, we use the coeffs of the McDonald paper. 
  index = np.arange(0,len(alpha)) # |j>
  #print alpha, index

  # Give values to alpha.
  kmin = 0.
  kmax = 2.*pi/AtomicDistance
  k = np.arange(0., kmax, (kmax-kmin) / kN ) #meshing |k> space until 2pi/a

  #print index, k
  #print Header+"|j><k|"
  #print np.outer(index, k) #generates |k><index|. How to compute |index><k| ?
  #print Header+"|k><j|"
  #print np.outer(k, index)

  cosTerm = np.cos(np.outer(k, index)*AtomicDistance)

  #print Header+"cosTerm"
  #print cosTerm

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

  plt.figure()
  plt.xlabel('k')
  plt.ylabel(r'$\varepsilon$_{||}')
  plt.plot(k,epsilonLong)
  plt.savefig("BandStructure.eps")
  plt.savefig("BandStructure.png")
  #plt.show()
          
  return epsilonLong_AU