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
# import libDatabase
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
def ScenarioOfOxidePrecipitation(epsMedium, epsSubstrate, epsEnvironment=1.): 
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

## Compute the spatial period and mean-free path of Surface Plasmon Polaritons 
# at the interface between three materials. 
# Here we investigate the role of the materials mixing happening pulse by pulse, oxygen from ambient air diffuses into the Cr, and leads to formation of Cr2O3. 
# By making use of a Lorenz-Lorentz model, one can construct a dielectric permittivity of the oxide
# and therefore observe the transition between SPP and waveguiding modes. 
# Note: Some results exist showing that L_spp < 0. These modes may be physical, and should be discussed 
# in the frame of the works of P. Berini. 
def ScenarioOfCrOxideMixture(epsSample, epsOxide=1., epsSubstrate=1., fraction_size=60, SampleName='Cr', OxideName='Cr2O3', SubstrateName='BK7', PlotLspp=False, FilterNegativeLspp=False):
    print(Header, "# Info: Considering a mixed fraction of Cr with Cr2O3 with several thicknesses.")
    wavelength = 1026e-9
    #fraction_size = 60
    numberofroots = 10 #BUG: does not work if numberofroots > 1. 
    
    fraction_min = 0.
    fraction_max = 1.

    thickness_size = 1
    thickness_min  = 28e-9
    thickness_max  = 28e-9
    
    #PlotLspp = False
    if(PlotLspp): 
        plotA=211; plotB=212
    else: 
        plotA=111; 
    ## Running 

    fraction = np.linspace(fraction_min, fraction_max, fraction_size) #fraction of Cr (includes the final value)
    epsCrCr2O3_list = MaxwellGarnett2(epsSample, epsOxide, 1.-fraction)
    print(Header, "# Info: size of the fraction matrix: ", fraction_size)

    print(Header, "# Info: preparation of the root finder.")
    
    #results = Parallel(n_jobs=num_cores)(delayed(processInput)(i) for i in inputs)

    #t = 100E-9 #thickness of the layer in meters
    #thickness = 10e-9
    print("\n")
    print("# Fraction of Cr: ", fraction)
    summary = np.zeros((0, 7))
    ## Preparation of the thin film modeling for various compositions
    for fraction_index in np.arange(0,fraction_size): #arange excludes the last one, linspace includes it
        # Medium 1: thin film. 
        eps1 = epsCrCr2O3_list[fraction_index]        #thin film
        # Medium 2: substrate. 
        eps2 = epsBK7       #epsBK7 #environment | substrate
        # Medium 3: environment
        eps3 = epsAir       #environment | substrate
        # Note: Inverting eps2 and eps3 should have no effect on the possible modes, but only on field amplification. 
        t_list = np.linspace(thickness_min, thickness_max, thickness_size, endpoint=True)
        
        #for thickness in t_list: #BUG: works only if using 1 thickness
        thickness = t_list[0] 
        roots = findroots(eps1, eps2, eps3,
                    wavelength, thickness,
                    x_min, x_max,    
                    y_min, y_max,    
                    x_steps, y_steps, numberofroots)

        # Preparation of the data for plotting
        print("\nOxide fraction #"+str(fraction_index)+"="+str(fraction[fraction_index]))
        num_branches = len(roots) 
        print("Number of branches: "+str(num_branches)) #number of SPP branches for this sample. 
        #print("ndimn: "+str(np.ndim(roots)))
        
        # This means we have tested <roots_oxidationDegree> samples with different oxidations. 
        
        #roots_t = np.ndarray(roots)
        
        #print(Header+"** Summary of the roots. Number of oxidation degrees: "+str(roots_oxidationDegree))
        #print("Roots")
        #print(roots) #For each branch, each sample, we have two data Period, and L_spp
        
        ### This version is not general enough for number_of_roots > 1. 
        num_thickness = np.shape(t_list) #NOTE: number of tested thicknesses
        for branch in np.arange(0,num_branches):
            roots_in_branch = roots[branch]
            #print("\n")
            print("Roots in branch #"+str(branch))
            for order in np.arange(0,len(roots_in_branch)):
                roots_in_branch_order = roots_in_branch[order]
                #print("\n")
                print("Order #"+str(order)+": Period="+str(roots_in_branch_order[0])+" Lspp="+str(roots_in_branch_order[1]))
                ToBeAdded = [thickness, branch, order, roots_in_branch_order[0], roots_in_branch_order[1], fraction[fraction_index], eps1]
                if(abs(roots_in_branch[order][0]) > 1E-15 and abs(roots_in_branch[order][1]) > 1E-10): 
                    # We remove modes were Lspp < 0.1 nm or period < 0. 
                    summary = np.vstack((summary, ToBeAdded ))
                    
        #num_branches  = roots_oxidationDegree[0]
        #return(roots, num_thickness, num_branches)
        #num_roots     = roots_oxidationDegree[1] #BUG: this crashes since MultilayerSPP repository was merged in develop... 
        #num_property  = roots_oxidationDegree[2]

        #for branch in np.arange(0,num_branches-1):
                #for root_number in np.arange(0,num_roots): 
                    #ToBeAdded = [thickness, branch, root_number, roots[branch][root_number][0], roots[branch][root_number][1], fraction[fraction_index], eps1]
                    #print(ToBeAdded)
                    #if(roots[branch][root_number][0] != 0e0): 
                        #
    #print(summary)
    #exit()
    #return(summary)
    
    print(Header+"** Preparation of the plots")

    ExperimentalData_velocity   = np.array([1e-6, 10e-16, 50e-6, 100e-6, 200e-6, 300e-6]) #m/s
    ExperimentalData_LSFL       = np.array([696e-9, 704e-9, 816e-9, 858e-9, -100e0, -100e-9]) #m #better observed for low velocities, i.e., high number of pulses, i.e., largest amounts of oxide
    ExperimentalData_LSFL_error = np.array([78e-9,71e-9,139e-9,140e-9,0e-9,0e-9]) #m
    ExperimentalData_HSFL       = np.array([170e-9, 159e-9, 217e-9, 244e-9, 249e-9, 238e-9]) #m
    ExperimentalData_HSFL_error = np.array([64e-9,38e-9,101e-9,110e-9,128e-9,72.5e-9]) #m
    CrFraction_Fitted           = np.array([0.6e0, 0.7e0, 0.8e0, 0.86e0, 0.90e0, 0.95e0]) #hand fitted to match period with existing modes [on request of Nadya]
    
    print(summary)

#== Extract the constructed table
    
    if(FilterNegativeLspp): 
        print("Filtering the negative Lspp / Imag(beta) < 0 ...")
        summary_filtered = np.array(summary[summary[:,4]>0,:])
        del summary
        summary = summary_filtered
        del summary_filtered
    
    ## We shall split the results by branch number. 
    # 1. Conditional filtering of the table for branch_s == 0, 1, 2 or 3. Similar filtering is available in libDatabase.py. 
    summary_branch0 = np.array(summary[summary[:,1]==0,:])
    summary_branch1 = np.array(summary[summary[:,1]==1,:])
    summary_branch2 = np.array(summary[summary[:,1]==2,:])
    summary_branch3 = np.array(summary[summary[:,1]==3,:])
    
    
    def SplitSummaryTable(summary):
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
            
            lists = sorted(zip(*[root_number, branch, thickness, period, lspp, fractionOxide, epsilonFilm]))
            root_number_s, branch_s, thickness_s, period_s, lspp_s, fractionOxide_s, epsilonFilm_s = list(zip(*lists))
        return root_number_s, branch_s, thickness_s, period_s, lspp_s, fractionOxide_s, epsilonFilm_s


    #ReBeta = np.divide(2.*np.pi,period_s)
    #ImBeta = np.divide(0.5,lspp_s)
    
    # Then we could plot them in the right order
    thickness, branch, root_number, period, lspp, fractionOxide, epsilonFilm = SplitSummaryTable(summary)
    
    thickness, branch0, root_number, period0, lspp0, fractionOxide0, epsilonFilm = SplitSummaryTable(summary_branch0)
    thickness, branch1, root_number, period1, lspp1, fractionOxide1, epsilonFilm = SplitSummaryTable(summary_branch1)
    thickness, branch2, root_number, period2, lspp2, fractionOxide2, epsilonFilm = SplitSummaryTable(summary_branch2)
    thickness, branch3, root_number, period3, lspp3, fractionOxide3, epsilonFilm = SplitSummaryTable(summary_branch3)
        
    print("SPP branches are ready. Exporting to CSV...")
    filename = "Dostovalov-SPPmodes-branch"
    np.savetxt(filename+"0"+".csv", summary_branch0)
    np.savetxt(filename+"1"+".csv", summary_branch1)
    np.savetxt(filename+"2"+".csv", summary_branch2)
    np.savetxt(filename+"3"+".csv", summary_branch3)
    
    plt.figure()
    ax1 = plt.subplot(plotA)
    #plt.title(r'Film thickness $t=$'+str(round(np.real(thickness[0])*1e9))+' nm')
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
    
    plot110, = ax1.plot(fractionOxide0, np.multiply(1e9,period0), 'r+', label=r'SPP period $\Lambda$, branch (-,-)')
    plot111, = ax1.plot(fractionOxide1, np.multiply(1e9,period1), 'k+', label=r'SPP period $\Lambda$, branch (-,+)')
    plot112, = ax1.plot(fractionOxide2, np.multiply(1e9,period2), 'b+', label=r'SPP period $\Lambda$, branch (+,-)')
    plot113, = ax1.plot(fractionOxide3, np.multiply(1e9,period3), 'go', label=r'SPP period $\Lambda$, branch ( +,+)')
    
    plot12, = ax1.plot(fractionOxide, np.multiply(1e9,wavelength*np.ones(np.shape(fractionOxide))), 'k-', linewidth=0.5, label=r'Laser wavelength $\lambda$')
    
    #plot14, = ax12.plot(fractionOxide_s, np.real(epsilonFilm_s), 'b+', label=r'Re($\varepsilon$)')
    #plot15, = ax12.plot(fractionOxide_s, np.imag(epsilonFilm_s), 'b^', label=r'Im($\varepsilon$)')
    
    #ax1.yaxis.label.set_color(plot110.get_color()) #colorizes the label
    #ax1.spines["left"].set_edgecolor(plot110.get_color()) #colorizes the axis
    #ax1.tick_params(axis='y', colors=plot110.get_color()) #colorizes the tics and numbers
    plotComb= [plot110, plot111, plot112, plot113, plot12]; 
    
    #ax12.yaxis.label.set_color(plot14.get_color()) #colorizes the label
    #ax12.spines["right"].set_edgecolor(plot14.get_color()) #colorizes the axis
    #ax12.tick_params(axis='y', colors=plot14.get_color()) #colorizes the tics and numbers
    
    plt.tight_layout()
    
    #plot1   = [plot11, plot12, plot14, plot15]
    #labels1 = [l.get_label() for l in plot1]
    #ax1.legend(plot1, labels1, loc='best')
    if(PlotLspp):
        ax2 = plt.subplot(plotB)
        plt.ylabel(r'$L_{SPP}$ decay length (m)')
        plt.xlabel(r'Fraction of Cr (perc.)')
        #plot21, = ax2.semilogy(fractionOxide_s, lspp_s,   'r^', label=r'SPP decay length $L_{SPP}$')
        plot211, = ax2.semilogy(fractionOxide0, np.abs(lspp0),   'r^', label=r'SPP decay length $L_{SPP}$, --')
        plot212, = ax2.semilogy(fractionOxide1, np.abs(lspp1),   'k^', label=r'SPP decay length $L_{SPP}$, -+')
        plot213, = ax2.semilogy(fractionOxide2, np.abs(lspp2),   'b^', label=r'SPP decay length $L_{SPP}$, +-')
        plot214, = ax2.semilogy(fractionOxide3, np.abs(lspp3),   'g^', label=r'SPP decay length $L_{SPP}$, ++')
        
        ax22 = ax2.twinx()    
        plt.ylabel(r'Re($\varepsilon$), Im($\varepsilon$)')
        #plot22, = ax22.plot(fractionOxide_s, np.real(epsilonFilm_s), 'b+', label=r'Re($\varepsilon$)')
        #plot23, = ax22.plot(fractionOxide_s, np.imag(epsilonFilm_s), 'b^', label=r'Im($\varepsilon$)')
        #plot24, = ax22.plot(fractionOxide, np.multiply(epsBK7, np.ones(np.shape(fractionOxide))), 'b--', label=r'$Re[\varepsilon$(BK7)] ')
        
        #ax2.yaxis.label.set_color(plot21.get_color()) #colorizes the label
        #ax2.spines["left"].set_edgecolor(plot21.get_color()) #colorizes the axis
        #ax2.tick_params(axis='y', colors=plot21.get_color()) #colorizes the tics and numbers
        
        #ax22.yaxis.label.set_color(plot22.get_color()) #colorizes the label
        #ax22.spines["right"].set_edgecolor(plot22.get_color()) #colorizes the axis
        #ax22.tick_params(axis='y', colors=plot22.get_color()) #colorizes the tics and numbers
        
        #plot2   = [plot21, plot22, plot23, plot24]
        plot2   = [plot211, plot212, plot213, plot214] #, plot24]
        plotComb+=plot2
        #labels2 = [l.get_label() for l in plot2]
        #ax2.legend(plot2, labels2, loc='best')
    labelsComb = [l.get_label() for l in plotComb]
    ax1.legend(plotComb, labelsComb, loc='upper left')
    
    
    plt.xlim((0,1))
    plt.tight_layout()
    filename="Dostovalov_"+OxideName+"-mixedWith-"+SampleName+"-Thickness-"+str(1E9*thickness_max)+"nm"
    plt.savefig(filename+".eps")
    plt.savefig(filename+".png")
    plt.show()
    
#}}}

# =====================
# ** Info: computing 3-layer reflectivity..."
#R = BiLayerReflectivity(epsAir, epsCrCr2O3_list, epsBK7, t_list) #dimension is good for a HeatMap picture

Fraction_size = 100 #samples

#ScenarioOfOxidePrecipitation()
#roots, num_thickness, num_branches = 
ScenarioOfCrOxideMixture(epsCr, epsCr2O3, epsBK7, Fraction_size, 'Cr', 'Cr2O3')
#ScenarioOfCrOxideMixture(epsCr, epsCrO2, epsBK7, Fraction_size,  'Cr', 'CrO2')

# Estimation of the sample heating

# n_opt = cmath.sqrt(epsilon)
# alpha = 2*omega/c * n_opt.imag
# S = alpha * I
# Equation: C_l * dT / dt = nabla( kappa \nabla ( T ) ) + S
