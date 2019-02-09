#!/usr/bin/env python3
#-*- coding: utf-8 -*-

# Copyright (C) 2019 T.J.-Y. Derrien
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

## @package MultilayerVanadium
# Preparation of results for F. Giorgianni (SwissFEL)
# Explores the SPP theory at a thin film located 
# between two semi-infinite media. The formal model is presented in 
# T.J.-Y. Derrien et al, J. Appl. Phys. 116, 074902 (2014) and references 
# therein. 

from libMultilayerSPP import *
from libMaterials     import *
import libUnits
from importPalikData import importFromEpsilonTable_batch
# import libDatabase
#from joblib import Parallel, delayed
#import multiprocessing

#num_cores = multiprocessing.cpu_count()

import matplotlib.pyplot as plt
from scipy.interpolate import InterpolatedUnivariateSpline
Header="# [vanadium.py]: "
# === PRODUCTION OF SCIENTIFIC RESULTS ===

#precision
NumberOfPoints=50

#data
wavelength = 85e-6
epsAir      = 1.+0.j #air
epsAl2O3    = 3.3635560000000115+0j #Palik
thickness = 80e-9
temperature = 300 #K

# optical data parameters
folder = "VanadiumOxides"

# Optical data for VO2
VO2={'300K_o': 'VO2-300K-EperpAaxis', 
     '355K_o': 'VO2-355K-EperpAaxis', 
     '300K_e': 'VO2-300K-EparallelAaxis', 
     '355K_e': 'VO2-355K-EparallelAaxis'}

select = str(temperature)+"K_"

# Building the optical data with same wavelength space, using a unit enabling to convert eV to wavelength
unit = libAtomicUnits.Energy_SI_to_Length(libAtomicUnits.Energy_eV_to_Joules(1.)) #
epsVO2_o_r = importFromEpsilonTable_batch(folder, VO2[select+"_o"+"_r"], False, unit)
epsVO2_e = importFromEpsilonTable_batch(folder, VO2[select+"_o"+"_r"], False, unit)


# Verification of the optical data validity by reflectivity computation

# mixing the materials that can coexist as amorphous Vanadium oxides. 


# computing the possible multilayer modes
