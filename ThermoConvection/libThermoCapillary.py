#!/usr/bin/env python
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

## @package libThermoCapillary
# Module libThermoCapillary defines the various regimes of thermo capillary instabilities found in the LIPSS formation. 
# We use as the following basics: 
# Chapter IV of PhD thesis: Derrien, T. J.-Y., Nanostructuring of solar cells by femtosecond laser irradiation. Theoretical study of the formation mechanisms. Université de la Méditerranée - Aix Marseille II, 2012. 
# Jean Berthier and Pascal Silberzan, "Microfluidics for Biotechnology", Artech House (2009).

def Term1(T, k, depth, dynamic_viscosity, density, surface_tension):
  num1 = dynamic_viscosity * k**2
  num2 = np.sqrt( (np.exp(k*h)**2 + 1e0 ) * (np.exp(k*h)**2 - 1e0) )
  num3 = np.sqrt( density * np.sqrt( gravity * depth * density + surface_tension * k**3) ) 
  num  = - np.sqrt(2e0) * np.sqrt( num1 * num2 * num3  )
  denom = 2e0*np.sinh(2e0*k*h)
  return num / denom

def Term2():