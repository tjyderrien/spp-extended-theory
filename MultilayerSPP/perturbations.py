#!/usr/bin/env python3
#-*- coding: utf-8 -*-

# Copyright (C) 2018 F. Preucil, T.J.-Y. Derrien
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

# @package libMultilayer
# Module libMultilayer explores the SPP theory at a thin film located 
# between two semi-infinite media. The formal model is presented in 
# T.J.-Y. Derrien et al, J. Appl. Phys. 116, 074902 (2014) and references 
# therein. 
import cmath
import numpy as np

#parameters
wavelength = 1
eps1 = 1
eps2 = 2
eps3 = 3

dt = 1 #dt << wavelength

#signs:
sbeta0 = 1
sgn1 = 1
sgn2 = -1
sgn3 = 1

k0 = 2.*np.pi/wavelength
beta0sq = k0**2*eps2*eps3/(eps2+eps3) #2-layer solution squared
beta0 = sbeta0*cmath.sqrt(beta0sq)    #square root including sign

k1 = sgn1*cmath.sqrt(beta0sq-k0**2*eps1)
k2 = sgn2*cmath.sqrt(beta0sq-k0**2*eps2)
k3 = sgn3*cmath.sqrt(beta0sq-k0**2*eps3)
dk1 = beta0/k1
dk2 = beta0/k2
dk3 = beta0/k3
dFdt = -2*k1*(k1/eps1-k2/eps2)*(k1/eps1-k3/eps3)
dFdbeta = (dk1/eps1-dk2/eps2)*(k1/eps1-k3/eps3)+(k1/eps1-k2/eps2)*(dk1/eps1-dk3/eps3)-(dk1/eps1+dk2/eps2)*(k1/eps1+k3/eps3)-(k1/eps1+k2/eps2)*(dk1/eps1+dk3/eps3)

dbeta = (-dFdt/dFdbeta)*dt
beta = beta0 + dbeta

print('beta0:', beta0)
print('dbeta:', dbeta)
print()
print('beta:', beta)
