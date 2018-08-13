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

## @package MultilayerDostovalov
# Preparation of results for Prof. Bulgakova and Sasha Dostovalov. 
# Explores the SPP theory at a thin film located 
# between two semi-infinite media. The formal model is presented in 
# T.J.-Y. Derrien et al, J. Appl. Phys. 116, 074902 (2014) and references 
# therein. 

from libMultilayerSPP import *
from libMaterials     import *
#from joblib import Parallel, delayed
#import multiprocessing

#num_cores = multiprocessing.cpu_count()

import matplotlib.pyplot as plt
Header="# [dostovalov.py]: "
# === PRODUCTION OF SCIENTIFIC RESULTS ===

#precision
NumberOfPoints=50

#data
wavelength = 1026e-9 #355e-9 #1030E-9 #1026

epsCr2O3   = 3.8273816+0.0483803j #Al-Kuhaili, M. & Durrani, S. Optical properties of chromium oxide thin films deposited by electron-beam evaporation Optical Materials, 2007, 29, 709-713
#epsCr2O3    = 4.9713+0.1784j  #1 um [JDT Kruschwitz et al, Appl. Opt. 1997]
epsCr      = -0.6721223+24.8657476j
epsCrO2   = 1.3587463082734004+9.00595525243578j #Chase, L. L. Optical properties of Cr O 2 and Mo O 2 from 0.1 to 6 eV Physical Review B, 1974, 10, 2226-2231


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
    eps2 = epsCrO2 #2O3       #epsBK7 #environment | substrate
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

## Scenario proposed by Nadya
# Pulse by pulse, oxygen from ambient air diffuses into the Cr, and leads to formation of Cr2O3. 
# By making use of a Lorenz-Lorentz model, one can construct a dielectric permittivity of the oxide
# and therefore observe the transition between SPP and waveguiding modes. 
def ScenarioOfCrOxideMixture():
    print(Header, "# Info: Considering a mixed fraction of Cr with Cr2O3 with several thicknesses.")
    wavelength = 1026e-9
    fraction_size = 60
    fraction_min = 0.
    fraction_max = 1.

    thickness_size = 1
    thickness_min  = 28e-9
    thickness_max  = 28e-9

    ## Running 

    fraction = np.linspace(fraction_min, fraction_max, fraction_size) #fraction of Cr
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

    ExperimentalData_velocity   = np.array([1e-6, 10e-16, 50e-6, 100e-6, 200e-6, 300e-6]) #m/s
    ExperimentalData_LSFL       = np.array([696e-9, 704e-9, 816e-9, 858e-9, -100e0, -100e-9]) #m #better observed for low velocities, i.e., high number of pulses, i.e., largest amounts of oxide
    ExperimentalData_LSFL_error = np.array([78e-9,71e-9,139e-9,140e-9,0e-9,0e-9]) #m
    ExperimentalData_HSFL       = np.array([170e-9, 159e-9, 217e-9, 244e-9, 249e-9, 238e-9]) #m
    ExperimentalData_HSFL_error = np.array([64e-9,38e-9,101e-9,110e-9,128e-9,72.5e-9]) #m
    CrFraction_Fitted           = np.array([0.6e0, 0.7e0, 0.8e0, 0.86e0, 0.90e0, 0.95e0]) #hand fitted to match period with existing modes [on request of Nadya]

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
    ax1 = plt.subplot(211)
    plt.title(r'Film thickness $t=$'+str(round(np.real(thickness[0])*1e9))+' nm')
    #plt.xlabel(r'Fraction of Cr (perc.)')
    plt.ylabel(r'SPP period $\Lambda$ (nm)') 
    #ax1y = ax1.twiny()
    plot1y1 = ax1.errorbar(CrFraction_Fitted, 1e9*ExperimentalData_LSFL, yerr=1e9*ExperimentalData_LSFL_error, fmt='ro', label=r'Period LSFL')
    plot1y2 = ax1.errorbar(CrFraction_Fitted, 1e9*ExperimentalData_HSFL, yerr=1e9*ExperimentalData_HSFL_error, fmt='r^', label=r'Period HSFL')
    #ax1.set_yscale('log')
    plt.ylim((0.,1.1e9*wavelength))
    #ax12 = ax1.twinx()
    #plt.ylabel(r'SPP mean free path $L_{SPP}$ (m)')
    #plt.ylabel(r'Re($\varepsilon$), Im($\varepsilon$)')
    # Then we could plot them in the right order
    plot11, = ax1.plot(fractionOxide_s, np.multiply(1e9,period_s), 'r+', label=r'Period $\Lambda$')
    plot12, = ax1.plot(fractionOxide_s, np.multiply(1e9,wavelength*np.ones(np.shape(thickness_s))), 'k--', label=r'$\lambda$')
    
    #plot14, = ax12.plot(fractionOxide_s, np.real(epsilonFilm_s), 'b+', label=r'Re($\varepsilon$)')
    #plot15, = ax12.plot(fractionOxide_s, np.imag(epsilonFilm_s), 'b^', label=r'Im($\varepsilon$)')
    
    ax1.yaxis.label.set_color(plot11.get_color()) #colorizes the label
    ax1.spines["left"].set_edgecolor(plot11.get_color()) #colorizes the axis
    ax1.tick_params(axis='y', colors=plot11.get_color()) #colorizes the tics and numbers
    
    #ax12.yaxis.label.set_color(plot14.get_color()) #colorizes the label
    #ax12.spines["right"].set_edgecolor(plot14.get_color()) #colorizes the axis
    #ax12.tick_params(axis='y', colors=plot14.get_color()) #colorizes the tics and numbers
    
    plt.tight_layout()
    
    #plot1   = [plot11, plot12, plot14, plot15]
    #labels1 = [l.get_label() for l in plot1]
    #ax1.legend(plot1, labels1, loc='best')
    
    ax2 = plt.subplot(212)
    plt.ylabel(r'$L_{SPP}$ decay length (m)')
    plt.xlabel(r'Fraction of Cr (perc.)')
    plot21, = ax2.semilogy(fractionOxide_s, lspp_s,   'r^', label=r'SPP decay length $L_{SPP}$')
    ax22 = ax2.twinx()    
    plt.ylabel(r'Re($\varepsilon$), Im($\varepsilon$)')
    plot22, = ax22.plot(fractionOxide_s, np.real(epsilonFilm_s), 'b+', label=r'Re($\varepsilon$)')
    plot23, = ax22.plot(fractionOxide_s, np.imag(epsilonFilm_s), 'b^', label=r'Im($\varepsilon$)')
    plot24, = ax22.plot(fractionOxide_s, np.multiply(epsBK7, np.ones(np.shape(fractionOxide_s))), 'b--', label=r'$Re[\varepsilon$(BK7)] ')
    
    ax2.yaxis.label.set_color(plot21.get_color()) #colorizes the label
    ax2.spines["left"].set_edgecolor(plot21.get_color()) #colorizes the axis
    ax2.tick_params(axis='y', colors=plot21.get_color()) #colorizes the tics and numbers
    
    ax22.yaxis.label.set_color(plot22.get_color()) #colorizes the label
    ax22.spines["right"].set_edgecolor(plot22.get_color()) #colorizes the axis
    ax22.tick_params(axis='y', colors=plot22.get_color()) #colorizes the tics and numbers
    
    plot2   = [plot21, plot22, plot23, plot24]
    plotComb= [plot11, plot12, plot21, plot22, plot23, plot24]
    #labels2 = [l.get_label() for l in plot2]
    labelsComb = [l.get_label() for l in plotComb]
    #ax2.legend(plot2, labels2, loc='best')
    ax2.legend(plotComb, labelsComb, loc='best')
    
    plt.tight_layout()
    plt.savefig("Dostovalov_Cr2O3mixedWithCr_Thickness.eps")
    plt.show()
    
#}}}

# =====================
# ** Info: computing 3-layer reflectivity..."
#R = BiLayerReflectivity(epsAir, epsCrCr2O3_list, epsBK7, t_list) #dimension is good for a HeatMap picture

#ScenarioOfOxidePrecipitation()
ScenarioOfCrOxideMixture()
