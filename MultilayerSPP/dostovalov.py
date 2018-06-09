#!/usr/bin/env python3
#-*- coding: utf-8 -*-

# Copyright (C) 2018 T.J.-Y. Derrien
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

# @package Multilayer_Dostovalov
# Preparation of results for Prof. Bulgakova and Sasha Dostovalov. 
# Explores the SPP theory at a thin film located 
# between two semi-infinite media. The formal model is presented in 
# T.J.-Y. Derrien et al, J. Appl. Phys. 116, 074902 (2014) and references 
# therein. 

from libMultilayerSPP import *
from libMaterials     import *
from joblib import Parallel, delayed
import multiprocessing

num_cores = multiprocessing.cpu_count()

import matplotlib.pyplot as plt
Header="# [dostovalov.py]: "
# === PRODUCTION OF SCIENTIFIC RESULTS ===

#precision
NumberOfPoints=50

#data
wavelength = 1026e-9 #355e-9 #1030E-9 #1026

epsCr2O3   = 3.8273816+0.0483803j
#epsCr2O3    = 4.9713+0.1784j  #1 um [JDT Kruschwitz et al, Appl. Opt. 1997]
epsCr      = -0.6721223+24.8657476j

epsBK7      = 2.10277365777   #1026 nm
epsSi       = 12.8159503769+0.0114635303918j #1026 nm, Palik
epsAir      = 1.+0.j          #air
#epsCu       = -46.6046581932 + 4.7188669976j #1030 nm
epsCu       = -1.9937293241+4.9290716854j     #355  nm

#meshes the initial guess area, all numbers are from the space of betas
x_min = -1E10
x_max = 1E10

y_min = -1E9
y_max = 1E9

x_steps = 40
y_steps = 40

# Scenario proposed by Thibault: an oxide layer grows at the top of the Cr sample, reducing progressively the periodicity by lambda/n. 
def ScenarioOfOxidePrecipitation(): 
    # Medium 1: thin film.
    eps1 = epsCr     #thin film
    # Medium 2: substrate. 
    eps2 = epsCr2O3       #epsBK7 #environment | substrate
    # Medium 3: environment
    eps3 = epsBK7       #environment | substrate
    # Note: Inverting eps2 and eps3 should have no effect on the possible modes, but only on field amplification. 
    
    thickness_size = 20
    thickness_min  = 10e-9
    thickness_max  = 300e-9
    
    t_list = np.linspace(thickness_min, thickness_max, thickness_size, endpoint=True)
    for thickness in t_list:
        roots = findroots(eps1, eps2, eps3,
                wavelength, thickness,
                x_min, x_max,    
                y_min, y_max,    
                x_steps, y_steps)

        # Shaping the data to plot them with GNUplot
        roots_shape = np.shape(roots)
        #print(roots_shape)
        num_thickness= np.shape(thickness)
        num_branches = roots_shape[0]
        num_property = roots_shape[1]

        for branch in np.arange(0,num_branches-1):
            print("1.", thickness, roots[branch][0], roots[branch][1], eps1.real, eps1.imag)

# Scenario proposed by Nadya
def ScenarioOfCrOxideMixture():
    print(Header, "# Info: Considering a mixed fraction of Cr with Cr2O3 with several thicknesses.")

    fraction_size = 20
    fraction_min = 0.5
    fraction_max = 1.0

    thickness_size = 1
    thickness_min  = 28e-9
    thickness_max  = 28e-9

    ## Running 

    fraction = np.arange(fraction_min, fraction_max, (fraction_max-fraction_min)/float(fraction_size)) #fraction of Cr
    epsCrCr2O3_list = MaxwellGarnett2(epsCr, epsCr2O3, 1.-fraction)
    print(Header, "# Info: size of the fraction matrix: ", fraction_size)

    print(Header, "# Info: preparation of the root finder.")
    
    #results = Parallel(n_jobs=num_cores)(delayed(processInput)(i) for i in inputs)

    #t = 100E-9 #thickness of the layer in meters
    #thickness = 10e-9

    print("# Fraction of Cr: ", fraction)
    summary = np.zeros((0, 7))
    ## Preparation of the thin film modeling for various compositions
    for fraction_index in np.arange(0,fraction_size): #arange excludes the last one
        # Medium 1: thin film. 
        eps1 = epsCrCr2O3_list[fraction_index]        #thin film
        # Medium 2: substrate. 
        eps2 = epsBK7       #epsBK7 #environment | substrate
        # Medium 3: environment
        eps3 = epsAir       #environment | substrate
        # Note: Inverting eps2 and eps3 should have no effect on the possible modes, but only on field amplification. 
        t_list = np.linspace(thickness_min, thickness_max, thickness_size, endpoint=True)
        
        
        for thickness in t_list:
            roots = findroots(eps1, eps2, eps3,
                    wavelength, thickness,
                    x_min, x_max,    
                    y_min, y_max,    
                    x_steps, y_steps)

        # Shaping the data to plot them with GNUplot
        roots_shape = np.shape(roots)
        #print(roots_shape)
        num_thickness = np.shape(thickness) #NOTE: is this used? 
        num_branches  = roots_shape[0]
        num_roots     = roots_shape[1]
        num_property  = roots_shape[2]

        for branch in np.arange(0,num_branches-1):
                for root_number in np.arange(0,num_roots): 
                    ToBeAdded = [thickness, branch, root_number, roots[branch][root_number][0], roots[branch][root_number][1], fraction[fraction_index], eps1]
                    print(ToBeAdded)
                    if(roots[branch][root_number][0] != 0e0): 
                        summary = np.vstack((summary, ToBeAdded ))

#== Extract the constructed table
    thickness     = summary[:, 0] 
    branch        = summary[:, 1]
    root_number   = summary[:, 2]
    period        = summary[:, 3]
    lspp          = summary[:, 4]
    fractionOxide = summary[:, 5]
    epsilonFilm   = summary[:, 6]

    print(summary)

    # Now, we shall sort out the data
    lists = sorted(zip(*[root_number, branch, thickness, period, lspp, fractionOxide, epsilonFilm]))
    root_number_s, branch_s, thickness_s, period_s, lspp_s, fractionOxide_s, epsilonFilm_s = list(zip(*lists))

    #ReBeta = np.divide(2.*np.pi,period_s)
    #ImBeta = np.divide(0.5,lspp_s)

    plt.figure()
    ax1 = plt.subplot(111)
    plt.xlabel(r'Fraction of Cr$_2$O$_3$ (perc.)')
    plt.ylabel('SPP period (m)')
    ax12 = ax1.twinx()
    plt.ylabel(r'SPP mean free path $L_{SPP}$ (m)')
    # Then we could plot them in the right order
    plot11, = ax1.loglog(fractionOxide, period_s, 'r+', label=r'Period $\Lambda$')
    plot12, = ax1.loglog(fractionOxide, wavelength*np.ones(np.shape(thickness_s)), '-', label=r'$\lambda$')
    plot13, = ax12.loglog(fractionOxide, lspp_s, 'b+', label=r'$L_{SPP}$')
    plt.tight_layout()

    plt.legend(loc='best')
    plt.tight_layout()
    plt.savefig("Dostovalov_Cr2O3mixedWithCr_Thickness.eps")
    plt.show()
    
#}}}

# =====================
# ** Info: computing 3-layer reflectivity..."
#R = BiLayerReflectivity(epsAir, epsCrCr2O3_list, epsBK7, t_list) #dimension is good for a HeatMap picture

ScenarioOfOxidePrecipitation()
