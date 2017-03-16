#!/usr/bin/env python
#-*- coding: utf-8 -*-

## Copyright (C) 2013-2017 T. J.-Y. Derrien
##
## This program is free software: you can redistribute it and/or modify
## it under the terms of the GNU General Public License as published by
## the Free Software Foundation, either version 3 of the License, or
## (at your option) any later version.
##
## This program is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.
##
## You should have received a copy of the GNU General Public License
## along with this program.  If not, see <http://www.gnu.org/licenses/>

## @package libLaser
# Functions to describe the laser pulse

import numpy as np
from scipy.constants import c, epsilon_0, mu_0, pi, e, m_e, h

## calculate laser frequency (Hz) from wavelength (m)
def omega(wavelength):#{{{
    return 2.0*pi*c/wavelength
#}}}

omega = np.vectorize(omega) 