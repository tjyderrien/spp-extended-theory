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

## Preparation of the multilayer SPP model. 
# This script aims to calculate the dispersion relation of 
# multilayer SPPs following the paper: 
# TJY Derrien et al, Journal of Applied Physics 116, 074902 (2014)

# DEFINITION OF BASIC FUNCTIONS

import cmath

## Just a simple test function. 
def simplereal(x):
    return x-1

## Two test equations, coupled
def doublereal(x):
    out = [x[0]*cos(x[1])-4, x[0]*x[1]-x[1]-5]
    # out.append(x[0]*x[1]-x[1]-5)
    return out
            
## The non-linear dispersion relation for multilayer SPP
# @param eps1 is the complex-valued dielectric permittivity composing the layer
    # @param eps2: see eps3.
    # @param eps3 are the complex-valued dielectric permittivities composing the bulk media 2,3 sandwitching the medium 1. 
    # @param t is the thickness of medium 1
    # @param k0 is the laser frequency
    # @param beta is the SPP wavenumber in layer 1 and is the unknown of the system, so just indicate a guess here
def dispersionSPP(eps1,eps2,eps3,t,k0,beta):    
    k1=cmath.sqrt(beta*beta-k0*k0 * eps1)
    k2=cmath.sqrt(beta*beta-k0*k0 * eps2)
    k3=cmath.sqrt(beta*beta-k0*k0 * eps3)
    
    eq = cmath.exp(-2.*k1*t) - ( (k1/eps1 + k2/eps2) / (k1/eps1 - k2/eps2) ) * ( (k1/eps1 + k3/eps3) / (k1/eps1 - k3/eps3) )
    return eq;

## Testing non-linear solver in R^2 space.
def f(z):
    # z^2+1=0, directly detailed in R^2 space
    temp=z[0]**2-z[1]**2+2j*z[0]*z[1]+1;
    return [temp.real, temp.imag]

## Testing non-linear solver in complex space. 
def f2(z):
    # z^2+1=0, without explicitly detailing real and imaginary part
    temp=cmath.exp(2. * cmath.log( z ))+1
    return temp;

## Interfacing complex space to R^2 space.
def c_to_r2(z):
    return [z.real, z.imag]

## Interfacing R^2 space to complex space
def r2_to_c(z):
    return z[0]+z[1]*1j
    
# we first define the dispersion relation for any epsilon1, 2, 3
# by finding root of the complex-valued function


# SOLVING THE EQUATION

import numpy as np
from scipy.optimize import fsolve, root
import cmath


# try to find root of simple real-valued function
# VALID
# print simplereal(1)
# x0=fsolve(simplereal, 10)
# print x0

# try to find root of a R^2->R valued function
# VALID
# guess=np.array([0., 10.])
# z=f(guess);
# x0=fsolve(f, guess)
# print 'z0=',x0


# try to find root of complex polynom function
# Real and imaginary parts are EXPLICIT
# let's solve a simple complex-valued equation
# VALID
# guess=np.array([1,2])
# z0=np.array([1,5])
# print f(z0)
# z1=root(f,z0)
# print 'Solution is', z1

# try to find root of a R^2->R valued function
# f2 is treated like an IMPLICIT expression into R2
# f2: C -> R
guessC=1-1j
guessR=np.array([guessC.real, guessC.imag]) #in R^2

print r2_to_c(guessR)

f2GuessC=f2(guessC);
f2GuessR=f2(r2_to_c(guessR))
print 'f2(', guessR,')=', f2GuessR
print 'f2(', guessC, ')=', f2GuessC

z1=fsolve([np.real(f2), np.imag(f2)], guessR)
# print c_to_r2(fz0)

# generates a mistake since the function must be real valued, but still defined as a function
# hence, real and imaginary part must be separated ? 
# z0=root(c_to_r2(f2), guess) 
# print z0




# PLOTTING THE RESULTS

import matplotlib as mp
from pylab import *

# let's plot a simple law
xmin=0.0
xmax=2.0
Nsteps=100
# x = linspace(xmin,xmax,Nsteps)
x=z0
y=z1
plot(x, y)

xlabel('X (m)')
ylabel('Y (m)')
title('Complex plane plot')
grid(True)
savefig("test.png")
show()
