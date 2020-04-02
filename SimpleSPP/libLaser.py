#!/usr/bin/env python
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

## @package libLaser
# Functions to describe the laser pulse

import numpy as np
from scipy.constants import c, epsilon_0, mu_0, pi, e, m_e, h, e

## calculate laser frequency (Hz) from wavelength (m)
def omega(wavelength):#{{{
    return 2.0*pi*c/wavelength
#}}}

## Compute photon wavelength (m) from a band gap energy (eV)
# @param Energy (eV)
# Output: corresponding photon wavelength (m)
def Energy_to_Wavelength(energy):
  wavelength = h * c / Energy / e
  return wavelength

## Converts photon wavelength (m) to a band gap energy (eV)
# @param Wavelength (m)
# Output: corresponding photon wavelength (m)
def Wavelength_to_Energy(wavelength):
  Energy = h * c / wavelength / e
  return Energy
  
omega                = np.vectorize(omega) 
Wavelength_to_Energy = np.vectorize(Wavelength_to_Energy)
