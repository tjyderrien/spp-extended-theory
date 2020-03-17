#!/usr/bin/env python2
#-*- coding: utf-8 -*-

# Copyright (C) 2013-2018 T. J.-Y. Derrien
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

## @package libNonDimensionalNumbers
# Module libNonDimensionalNumbers simply computes the non-dimensional numbers
# of the Navier-Stokes equation to define the regime of matter flowing. 
# We use as the following basics: 
# Chapter IV of PhD thesis: Derrien, T. J.-Y., Nanostructuring of solar cells by femtosecond laser irradiation. Theoretical study of the formation mechanisms. Université de la Méditerranée - Aix Marseille II, 2012. 
# Jean Berthier and Pascal Silberzan, "Microfluidics for Biotechnology", Artech House (2009).

from libThermalProperties_Silicon import *
#from libThermalProperties_Silica import *


Header="libNonDimensionalNumbers: "

## Collision time for the matter momentum damping
# @param density: density of the liquid
# @param viscosity: DYNAMIC viscosity of the liquid
# @param length: caracteristic length of the problem
def MomentumDampingTime(density, dynamic_viscosity, length):
  tau_damping = density * length**2 / dynamic_viscosity
  return tau_damping

def PressureInducedVelocity(pressure, dynamic_viscosity, length):
  velocity = pressure * length / dynamic_viscosity
  duration = length / velocity
  return velocity, duration

def SurfaceTensionInducedVelocity(SurfaceTension, depth, length, dynamic_viscosity):
  velocity = 2e0 * SurfaceTension * depth / (length * dynamic_viscosity)
  duration = length / velocity
  return velocity, duration

def PressureInducedNonlinearVelocity(pressure, density):
  velocity = np.sqrt( pressure / density )
  duration = length / velocity
  return velocity, duration

def SurfaceTensionInducedNonlinearVelocity(SurfaceTension, depth, density, length):
  velocity = np.sqrt( 2e0 * SurfaceTension * depth / (density * length**2) )
  duration = length / velocity
  return velocity, duration

def ReynoldsNumber(density, velocity, length, dynamic_viscosity):
  return density * velocity * length / dynamic_viscosity

def CapillaryNumber(dynamic_viscosity, velocity, SurfaceTension):
  return dynamic_viscosity * velocity / SurfaceTension

def WebberNumber(density, velocity, length, SurfaceTension): 
  return density * velocity**2 * length / SurfaceTension

## Peclet number
# @param kinematic_viscosity [m2/s]: 
# @param diffusivity [m2/s]: 
def PecletNumber(kinematic_viscosity, diffusivity):
  Pe = kinematic_viscosity / diffusivity
  return Pe

MomentumDampingTime                     = np.vectorize(MomentumDampingTime                       )
PressureInducedVelocity                 = np.vectorize(PressureInducedVelocity                   )
SurfaceTensionInducedVelocity           = np.vectorize(SurfaceTensionInducedVelocity             )
PressureInducedNonlinearVelocity        = np.vectorize(PressureInducedNonlinearVelocity          )
SurfaceTensionInducedNonlinearVelocity  = np.vectorize(SurfaceTensionInducedNonlinearVelocity    )
ReynoldsNumber                          = np.vectorize(ReynoldsNumber                            )
CapillaryNumber                         = np.vectorize(CapillaryNumber                           )
WebberNumber                            = np.vectorize(WebberNumber                              )
PecletNumber                            = np.vectorize(PecletNumber)
