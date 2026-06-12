#!/usr/bin/env python3
#-*- coding: utf-8 -*-

# Copyright (C) 2018-2022 T.J.-Y. Derrien
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

## @package MultilayerMo
# Preparation of results for .Mo film on SiO2 substrate
# Explores the SPP theory at a thin film located 
# between two semi-infinite media. The formal model is presented in 
# T.J.-Y. Derrien et al, J. Appl. Phys. 116, 074902 (2014) and references 
# therein.
# Dostovalov, A. V.; Derrien, T. J.-Y.; Appl. Surf. Sci., 2019, 491, 650-658 and references therein.

from spp_extended_theory.MultilayerSPP import libMultilayerSPP as lml
from spp_extended_theory.Libs import libMaterials  as lmat
# import libDatabase
#from joblib import Parallel, delayed
#import multiprocessing

#num_cores = multiprocessing.cpu_count()

import numpy as np
import sys
import matplotlib.pyplot as plt
from spp_extended_theory.MultilayerSPP import libMultilayerSPP as lml
from scipy.interpolate import InterpolatedUnivariateSpline
Header="# [dostovalov.py]: "
# === PRODUCTION OF SCIENTIFIC RESULTS ===

#precision
NumberOfPoints=50

#data
            
#epsBK7      = 2.10277365777   #1030 nm
epsAir      = 1.+0.j          #air
#epsCu       = -46.6046581932 + 4.7188669976j #1030 nm

#meshes the initial guess area, all numbers are from the space of betas
x_min = -6e7       #-1E10
x_max = 6e7       #1E10
               #
y_min = -1E9
y_max = 1E9
        
x_steps =  30    #70
y_steps =  30    #70

## Reorganizes the order of fields and output only necessary information
def SplitSummaryTable(summary): #{{{
    if(len(summary[:,0])==0): #if table is empty, avoids the crash
        thickness_s = 0e0; branch_s=-1; root_number_s=-1; period_s=-1; lspp_s=-1; fractionOxide_s=-1; epsilonFilm_s=1E99
    else:
        thickness     = summary[:, 0] 
        branch        = summary[:, 1]
        root_number   = summary[:, 2]
        period        = summary[:, 3]
        lspp          = summary[:, 4]
        fractionOxide = summary[:, 5]
        epsilonFilm   = summary[:, 6]
        
        lists = sorted(zip(*[thickness, fractionOxide, epsilonFilm, branch, root_number, period, lspp]))
        thickness_s, fractionOxide_s, epsilonFilm_s, branch_s, root_number_s, period_s, lspp_s  = list(zip(*lists))
    return thickness_s, fractionOxide_s, epsilonFilm_s, branch_s, root_number_s, period_s, lspp_s
#}}}

# Thin film of metal is deposited at the surface of a dielectric substrate in air atmosphere
# and irradiated by wavelength. 
def Derrien_HRLIPSSonMoFilms(wavelength, thickness_size):
    Header="[Mo Film: Derrien_HRLIPSSonMoFilms: ]"
    if( wavelength==1030E-9):
        epsCu       = -46.6046581932+4.7188669976j #Palik
        epsSi       = 12.80259+0.0109j #Palik
        epsSiO2     = 2.1026565205
        epsMo       = -11.6291789477+20.6107572133j
    elif(wavelength==800E-9):
        epsCu       = -25.2739024871+2.5157955079j #Palik
        epsSi       = 13.6338991279 + 0.0479330857j  #Palik
        epsSiO2     = 2.1124739397
    elif(wavelength==515E-9):
        epsCu       = -5.4822221975+5.821281567j
        epsSi       = 17.8251990451 + 0.5066899508j
        epsSiO2     = 2.1361639897
    elif(wavelength==400E-9):
        epsCu       = -3.4917+5.2156j
        epsSi       = 30.8542847158+4.300769121j
        epsSiO2     = 2.1614446988
    else: 
        print((Header,"No optical data were provided for this input."))
        sys.exit()
 
    numberofroots = 10
    # Medium 1: thin film.
    eps1 = epsMo     #Au thin film
    # Medium 2: substrate. 
    eps2 = epsSiO2                 
    # Medium 3: environment
    eps3 = 1.0 #.5**2 # eps2 #eps2: symmetric modes       #environment | substrate
    # Note: Inverting eps2 and eps3 should have no effect on the possible modes, but only on field amplification. 
    
    #thickness_size = 20
    thickness_min  = 0.1e-9
    thickness_max  = 300e-9
    t_list_log = np.linspace(np.log10(thickness_min), np.log10(thickness_max), thickness_size, endpoint=True)
    t_list = np.power(10., t_list_log)
    
    summary = np.zeros((0, 7))
    for thickness in t_list:
        roots = lml.findroots(eps1, eps2, eps3,
                wavelength, thickness,
                x_min, x_max,    
                y_min, y_max,    
                x_steps, y_steps, numberofroots)

        ## Shaping the data to plot them with GNUplot
        #roots_shape = np.shape(roots)
        ##print(roots_shape)
        #num_thickness= np.shape(thickness)
        #num_branches = roots_shape[0]
        #num_property = roots_shape[1]

        #for branch in np.arange(0,num_branches-1):
            #print("1.", thickness, roots[branch][0], roots[branch][1], eps1.real, eps1.imag)
        
        num_branches = len(roots) 
        
        num_thickness = np.shape(t_list) #NOTE: number of tested thicknesses
        for branch in np.arange(0,num_branches):
            roots_in_branch = roots[branch]
            #print("\n")
            print(("Roots in branch #"+str(branch)))
            for order in np.arange(0,len(roots_in_branch)):
                roots_in_branch_order = roots_in_branch[order]
                #print("\n")
                fraction = 0e0 #irrelevant in this context
                print(("Thickness:"+str(thickness)+", Order #"+str(order)+": Period="+str(roots_in_branch_order[0])+" Lspp="+str(roots_in_branch_order[1])))
                print(("beta/k0="+str(lml.PeriodToBetaNorm(roots_in_branch_order[0], wavelength))))
                ToBeAdded = [thickness, branch, order, roots_in_branch_order[0], roots_in_branch_order[1], fraction, eps1]
                if(abs(roots_in_branch[order][0]) > 1E-15 and abs(roots_in_branch[order][1]) > 1E-10): 
                    # We remove modes were |Lspp| < 0.1 nm or |period| < 0. 
                    summary = np.vstack((summary, ToBeAdded ))
    
    print(summary)
    ## Exporting the results 
    summary_branch0 = np.array(summary[summary[:,1]==0,:]) #-,-
    summary_branch1 = np.array(summary[summary[:,1]==1,:]) #-,+
    summary_branch2 = np.array(summary[summary[:,1]==2,:]) #+,-
    summary_branch3 = np.array(summary[summary[:,1]==3,:]) #+,+
    
    #ReBeta = np.divide(2.*np.pi,period_s)
    #ImBeta = np.divide(0.5,lspp_s)
    
    # Then we could plot them in the right order
    #thickness, branch, root_number, period, lspp, fractionOxide, epsilonFilm = SplitSummaryTable(summary)
    
    thickness, fractionOxide, epsilonFilm, branch, root_number, period, lspp = SplitSummaryTable(summary)
    
    thickness0, fractionOxide0, epsilonFilm0, branch0, root_number0, period0, lspp0 = SplitSummaryTable(summary_branch0)
    thickness1, fractionOxide1, epsilonFilm1, branch1, root_number1, period1, lspp1 = SplitSummaryTable(summary_branch1)
    thickness2, fractionOxide2, epsilonFilm2, branch2, root_number2, period2, lspp2 = SplitSummaryTable(summary_branch2)
    thickness3, fractionOxide3, epsilonFilm3, branch3, root_number3, period3, lspp3 = SplitSummaryTable(summary_branch3)
    
    summary_branch0_export = np.array([np.real(thickness0), np.real(fractionOxide0), np.real(epsilonFilm0), np.imag(epsilonFilm0), np.real(branch0), np.real(root_number0), np.real(period0), np.real(lspp0)])
    summary_branch1_export = np.array([np.real(thickness1), np.real(fractionOxide1), np.real(epsilonFilm1), np.imag(epsilonFilm1), np.real(branch1), np.real(root_number1), np.real(period1), np.real(lspp1)])
    summary_branch2_export = np.array([np.real(thickness2), np.real(fractionOxide2), np.real(epsilonFilm2), np.imag(epsilonFilm2), np.real(branch2), np.real(root_number2), np.real(period2), np.real(lspp2)])
    summary_branch3_export = np.array([np.real(thickness3), np.real(fractionOxide3), np.real(epsilonFilm3), np.imag(epsilonFilm3), np.real(branch3), np.real(root_number3), np.real(period3), np.real(lspp3)])
    
    summary_branch0_export_t = np.transpose(summary_branch0_export)
    summary_branch1_export_t = np.transpose(summary_branch1_export)
    summary_branch2_export_t = np.transpose(summary_branch2_export)
    summary_branch3_export_t = np.transpose(summary_branch3_export)
        
    print("SPP branches are ready. Exporting to CSV...")
    filename = "Derrien2019-SPPmodes-branch"
    np.savetxt(filename+"0"+".csv", summary_branch0_export_t)
    np.savetxt(filename+"1"+".csv", summary_branch1_export_t)
    np.savetxt(filename+"2"+".csv", summary_branch2_export_t)
    np.savetxt(filename+"3"+".csv", summary_branch3_export_t)
    
    # plt.figure()
    fig, (ax1, ax2) = plt.subplots(2, sharex=True)
    plt.xlabel(r'Thickness (nm)')
    ax1.set_ylabel(r'SPP period $\Lambda$ (nm)')
    # ax1.set_yscale('log')
    plot110, = ax1.plot(np.multiply(1e9,thickness0), np.multiply(1e9,period0), 'r+', label=r'SPP period $\Lambda$, branch (-,-)')
    plot111, = ax1.plot(np.multiply(1e9,thickness1), np.multiply(1e9,period1), 'k+', label=r'SPP period $\Lambda$, branch (-,+)')
    plot112, = ax1.plot(np.multiply(1e9,thickness2), np.multiply(1e9,period2), 'b+', label=r'SPP period $\Lambda$, branch (+,-)')
    plot113, = ax1.plot(np.multiply(1e9,thickness3), np.multiply(1e9,period3), 'g+', label=r'SPP period $\Lambda$, branch ( +,+)')
    ax1.set_ylim((0.,1.1e9*wavelength))
    
    # plt.figure()
    # ax2 = plt.subplot(212)
    ax2.set_xlabel(r'Thickness (nm)')
    ax2.set_ylabel(r'SPP mean free path $L_{SPP}$ (m)')
    plot110, = ax2.plot(np.multiply(1e9,thickness0), np.multiply(1e0,np.abs(lspp0)), 'r+', label=r'SPP period $\Lambda$, branch (-,-)')
    plot111, = ax2.plot(np.multiply(1e9,thickness1), np.multiply(1e0,np.abs(lspp1)), 'k+', label=r'SPP period $\Lambda$, branch (-,+)')
    plot112, = ax2.plot(np.multiply(1e9,thickness2), np.multiply(1e0,np.abs(lspp2)), 'b+', label=r'SPP period $\Lambda$, branch (+,-)')
    plot113, = ax2.plot(np.multiply(1e9,thickness3), np.multiply(1e0,np.abs(lspp3)), 'g+', label=r'SPP period $\Lambda$, branch ( +,+)')
    ax2.set_ylim((20E-9,1E-3))
    ax2.set_xscale('log')
    ax2.set_yscale('log')
    
    #plot12,  = ax1.plot(np.multiply(1e9, thickness), np.multiply(1e9,wavelength*np.ones(np.shape(fractionOxide))), 'k-', linewidth=0.5, label=r'Laser wavelength $\lambda$')
    
    #ax12 = ax1.twinx()
    #plot14, = ax12.plot(fractionOxide_s, np.real(epsilonFilm_s), 'b+', label=r'Re($\varepsilon$)')
    #plot15, = ax12.plot(fractionOxide_s, np.imag(epsilonFilm_s), 'b^', label=r'Im($\varepsilon$)')
    #plot14,  = ax1.plot(np.multiply(1e9,thickness0), np.multiply(1e0,period0), 'r^', label=r'Period $\Lambda$, (-,-)')
    #plot15,  = ax1.plot(np.multiply(1e9,thickness1), np.multiply(1e0,period1), 'k^', label=r'Period $\Lambda$, (-,+)')
    #plot16,  = ax1.plot(np.multiply(1e9,thickness2), np.multiply(1e0,period2), 'b^', label=r'Period $\Lambda$, (+,-)')
    #plot17,  = ax1.plot(np.multiply(1e9,thickness3), np.multiply(1e0,period3), 'g^', label=r'Period $\Lambda$, (+,+)')
    
    #plt.ylabel(r'Re($\varepsilon$), Im($\varepsilon$)')
    #plt.ylabel(r'$\beta/k_0$')
    #ax1.yaxis.label.set_color(plot110.get_color()) #colorizes the label
    #ax1.spines["left"].set_edgecolor(plot110.get_color()) #colorizes the axis
    #ax1.tick_params(axis='y', colors=plot110.get_color()) #colorizes the tics and numbers
    plotComb1 = []
    plotComb1+= [plot110, plot111, plot112, plot113]; 
    #plotComb1+=[plot14, plot15, plot16, plot17]
    
    #ax12.yaxis.label.set_color(plot14.get_color()) #colorizes the label
    #ax12.spines["right"].set_edgecolor(plot14.get_color()) #colorizes the axis
    #ax12.tick_params(axis='y', colors=plot14.get_color()) #colorizes the tics and numbers
    plt.tight_layout()
    labelsComb1 = [l.get_label() for l in plotComb1]
    #ax1.legend(plotComb1, labelsComb1, loc="best")
    fig.savefig("Air-MoFilm-SiO2-PeriodLsppWithThickness.eps")
    plt.show()



# =====================
#ScenarioOfOxidePrecipitation()
# ** Info: computing 3-layer reflectivity..."
#R = BiLayerReflectivity(epsAir, epsCrCr2O3_list, epsBK7, t_list) #dimension is good for a HeatMap picture

## Validation cases in Python. 
wavelength = 1030E-9 # 800e-9 # 1030e-9
Fraction_size = 30 #number of samples
thickness_size = 60 #Fraction_size
#Burke_SymmetricModes(thickness_size)
#Derrien_HRLIPSSonAuFilms(1030e-9, thickness_size)
#Derrien_HRLIPSSonAuFilms(800e-9, thickness_size)
#Derrien_HRLIPSSonAuFilms(515e-9, thickness_size)
#Derrien_HRLIPSSonAuFilms(400e-9, thickness_size)

#Derrien_HRLIPSSonAlFilms(1030e-9, thickness_size)
#Derrien_HRLIPSSonAlFilms(800e-9, thickness_size)
#Derrien_HRLIPSSonAlFilms(515e-9, thickness_size)
#Derrien_HRLIPSSonAlFilms(400e-9, thickness_size)

Derrien_HRLIPSSonMoFilms(1030e-9, thickness_size)
#Derrien_HRLIPSSonMoFilms(800e-9, thickness_size)
#Derrien_HRLIPSSonMoFilms(515e-9, thickness_size)
#Derrien_HRLIPSSonMoFilms(400e-9, thickness_size)

# Derrien_HRLIPSSonCuFilms(wavelength, thickness_size)

sys.exit()

