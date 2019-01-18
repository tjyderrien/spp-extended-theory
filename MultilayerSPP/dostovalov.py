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
from scipy.interpolate import InterpolatedUnivariateSpline
Header="# [dostovalov.py]: "
# === PRODUCTION OF SCIENTIFIC RESULTS ===

#precision
NumberOfPoints=50

#data
wavelength = 1026e-9 #355e-9 #1030E-9 #1026

epsCr2O3   = 3.82738158014083 + 0.0483802637311967j  #Al-Kuhaili, M. & Durrani, S. Optical properties of chromium oxide thin films deposited by electron-beam evaporation Optical Materials, 2007, 29, 709-713
#epsCr2O3    = 4.9713+0.1784j  #1 um [JDT Kruschwitz et al, Appl. Opt. 1997]
epsCr      = -0.6721223+24.8657476j
epsCrO2_e    = 1.3587463082734004+9.00595525243578j #Chase, L. L. Optical properties of Cr O 2 and Mo O 2 from 0.1 to 6 eV Physical Review B, 1974, 10, 2226-2231 (E || C mode. C: axis for extraordinary mode). 
#epsCrO2_e = epsCrO2 #(E || c)
epsCrO2_o  = 0.5784455474251002+6.477144040000001j #Chase, L. L. Optical properties of Cr O 2 and Mo O 2 from 0.1 to 6 eV Physical Review B, 1974, 10, 2226-2231 (E perpendicular to C).

epsTi = -4.289599704142011+27.217715606508875j

epsBK7      = 2.10277365777   #1026 nm
epsSi       = 12.8159503769+0.0114635303918j #1026 nm, Palik
epsAir      = 1.+0.j          #air
#epsCu       = -46.6046581932 + 4.7188669976j #1030 nm
epsCu       = -1.9937293241+4.9290716854j     #355  nm

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

# Scenario proposed by Thibault: an oxide layer grows at the top of the Cr sample, reducing progressively the periodicity by lambda/n. 
# This function can also be used to generate simple results for verification
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

## Preparation of the excited SPP as function of the electron temperature (estimation)
def PeriodsAsFunctionOfTemperature(epsSample, epsSubstrate=1., epsEnvironment=1., Te_size=10, SampleName='Cr', SubstrateName='BK7', PlotLspp=False, FilterNegativeLspp=False, PlotEpsilons=True, Te_max=1E6):
    
    ## In this function, Te should be interpreted as state of matter (in term of electron temperature)
    
    print(Header, "# Info: Considering temperature Te of electron in Cr with Cr2O3 with several thicknesses.")
    wavelength = 1026e-9
    #Te_size = 60
    numberofroots = 10 #per branch. 10 exceeds the final number of roots per branch
    
    Te_min = 300e0
    #Te_max = 1E6

    thickness_size = 1
    thickness_min  = 50e-9 #28e-9 #28e-9
    thickness_max  = thickness_min
    
    if(PlotLspp and PlotEpsilons):
        plotA = 311; plotB=312; plotC=313
    elif((PlotLspp and not PlotEpsilons) or (not PlotLspp and PlotEpsilons)): 
        plotA=211; plotB=212
    else: 
        plotA=111; 
        
    ## Running 
    logTe = np.linspace(np.log10(Te_min), np.log10(Te_max), Te_size)
    Te = np.power(10, logTe)
    epsCr_list = Drude_Ti(wavelength, epsSample, Te) #Ti opt data
    #epsCr_list = Drude_Cr(wavelength, epsSample, Te)
    # MaxwellGarnett2(epsSample, epsOxide, Te) #Te refers to the electron temperature Te here! 
    #NOTE: epsOxide here is Cr oxide, not environment. 
    print(Header, "# Info: size of the Te matrix: ", Te_size)
    print(Header, "# Info: preparation of the root finder.")
    
    #results = Parallel(n_jobs=num_cores)(delayed(processInput)(i) for i in inputs)

    #t = 100E-9 #thickness of the layer in meters
    #thickness = 10e-9
    print("\n")
    print("# Te of electrons in Cr: ", Te)
    summary = np.zeros((0, 7))
    ## Preparation of the thin film modeling for various compositions
    for Te_index in np.arange(0,Te_size): #arange excludes the last one, linspace includes it
        # Medium 1: thin film. 
        eps1 = epsCr_list[Te_index]        #thin film
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
        print("\nTe #"+str(Te_index)+"="+str(Te[Te_index]))
        num_branches = len(roots) 
        print("Number of branches: "+str(num_branches)) #number of SPP branches for this sample. 
        #print("ndimn: "+str(np.ndim(roots)))
        
        # This means we have tested <roots_oxidationDegree> samples with different oxidations. 
        
        #roots_t = np.ndarray(roots)
        
        #print(Header+"** Summary of the roots. Number of oxidation degrees: "+str(roots_oxidationDegree))
        #print("Roots")
        #print(roots) #For each branch, each sample, we have two data Period, and L_spp
        
        ### This version is general enough for number_of_roots > 1. 
        num_thickness = np.shape(t_list) #NOTE: number of tested thicknesses
        for branch in np.arange(0,num_branches):
            roots_in_branch = roots[branch]
            #print("\n")
            print("Roots in branch #"+str(branch))
            for order in np.arange(0,len(roots_in_branch)):
                roots_in_branch_order = roots_in_branch[order]
                #print("\n")
                print("Order #"+str(order)+": Period="+str(roots_in_branch_order[0])+" Lspp="+str(roots_in_branch_order[1]))
                ToBeAdded = [thickness, branch, order, roots_in_branch_order[0], roots_in_branch_order[1], Te[Te_index], eps1]
                if(roots_in_branch[order][0] > 1E-15 and abs(roots_in_branch[order][1]) > 1E-10): 
                    # We remove modes were Lspp < 0.1 nm or period < 0. 
                    summary = np.vstack((summary, ToBeAdded ))
                    
        #num_branches  = roots_oxidationDegree[0]
        #return(roots, num_thickness, num_branches)
        #num_roots     = roots_oxidationDegree[1] #BUG: this crashes since MultilayerSPP repository was merged in develop... 
        #num_property  = roots_oxidationDegree[2]

        #for branch in np.arange(0,num_branches-1):
                #for root_number in np.arange(0,num_roots): 
                    #ToBeAdded = [thickness, branch, root_number, roots[branch][root_number][0], roots[branch][root_number][1], Te[Te_index], eps1]
                    #print(ToBeAdded)
                    #if(roots[branch][root_number][0] != 0e0): 
                        #
    #print(summary)
    #exit()
    #return(summary)
    
    print(Header+"** Preparation of the plots as function of oxide ratio")

    ExperimentalData_velocity   = np.array([1e-6, 10e-16, 50e-6, 100e-6, 200e-6, 300e-6]) #m/s
    ExperimentalData_LSFL       = np.array([696e-9, 704e-9, 816e-9, 858e-9, -100e0, -100e-9]) #m #better observed for low velocities, i.e., high number of pulses, i.e., largest amounts of oxide
    ExperimentalData_LSFL_error = np.array([78e-9,71e-9,139e-9,140e-9,0e-9,0e-9]) #m
    ExperimentalData_HSFL       = np.array([170e-9, 159e-9, 217e-9, 244e-9, 249e-9, 238e-9]) #m
    ExperimentalData_HSFL_error = np.array([64e-9,38e-9,101e-9,110e-9,128e-9,72.5e-9]) #m
    CrFraction_Fitted           = np.array([0.6e0, 0.7e0, 0.8e0, 0.86e0, 0.90e0, 0.95e0]) #hand fitted to match period with existing modes [on request of Nadya]
    CrO2Fraction_Fitted         = np.add(1, - np.array([0.6e0, 0.7e0, 0.8e0, 0.86e0, 0.90e0, 0.95e0])) #we express it in term of oxidized
    
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
    
    #ReBeta = np.divide(2.*np.pi,period_s)
    #ImBeta = np.divide(0.5,lspp_s)
    
    # Then we could plot them in the right order
    thickness, Te, epsilonFilm, branch, root_number, period, lspp = SplitSummaryTable(summary)
    
    thickness0, Te0, epsilonFilm0, branch0, root_number0, period0, lspp0 = SplitSummaryTable(summary_branch0)
    thickness1, Te1, epsilonFilm1, branch1, root_number1, period1, lspp1 = SplitSummaryTable(summary_branch1)
    thickness2, Te2, epsilonFilm2, branch2, root_number2, period2, lspp2 = SplitSummaryTable(summary_branch2)
    thickness3, Te3, epsilonFilm3, branch3, root_number3, period3, lspp3 = SplitSummaryTable(summary_branch3)
    
    summary_export         = np.array([np.real(Te) , np.real(epsilonFilm) , np.imag(epsilonFilm) , np.real(branch) , np.real(root_number) , np.real(period) ])
    summary_branch0_export = np.array([np.real(Te0), np.real(epsilonFilm0), np.imag(epsilonFilm0), np.real(branch0), np.real(root_number0), np.real(period0)])
    summary_branch1_export = np.array([np.real(Te1), np.real(epsilonFilm1), np.imag(epsilonFilm1), np.real(branch1), np.real(root_number1), np.real(period1)])
    summary_branch2_export = np.array([np.real(Te2), np.real(epsilonFilm2), np.imag(epsilonFilm2), np.real(branch2), np.real(root_number2), np.real(period2)])
    summary_branch3_export = np.array([np.real(Te3), np.real(epsilonFilm3), np.imag(epsilonFilm3), np.real(branch3), np.real(root_number3), np.real(period3)])
    
    summary_export_t         = np.transpose(summary_export)
    summary_branch0_export_t = np.transpose(summary_branch0_export)
    summary_branch1_export_t = np.transpose(summary_branch1_export)
    summary_branch2_export_t = np.transpose(summary_branch2_export)
    summary_branch3_export_t = np.transpose(summary_branch3_export)

    # NOTE: I would like to get clean numbers, not the content of summary_branch0
    # Is there a problem with root_number0 for example? 
    print("SPP branches are ready. Exporting to CSV...")
    filename = "Dostovalov-SPPmodes-Cr-CrO2"
    np.savetxt(filename+".csv", summary_export_t)
    np.savetxt(filename+"-branch0"+".csv", summary_branch0_export_t)
    np.savetxt(filename+"-branch1"+".csv", summary_branch1_export_t)
    np.savetxt(filename+"-branch2"+".csv", summary_branch2_export_t)
    np.savetxt(filename+"-branch3"+".csv", summary_branch3_export_t)
    
    plt.figure()
    ax1 = plt.subplot(plotA)
    #plt.title(r'Film thickness $t=$'+str(round(np.real(thickness[0])*1e9))+' nm')
    #plt.xlabel(r'Fraction of Cr oxide  (perc.)')
    plt.ylabel(r'SPP period $\Lambda$ (nm)') 
    #ax1y = ax1.twiny()
    #plot1y1 = ax1.errorbar(np.multiply(CrFraction_Fitted,100), 1e9*ExperimentalData_LSFL, yerr=1e9*ExperimentalData_LSFL_error, fmt='ro', label=r'Period LSFL')
    #plot1y2 = ax1.errorbar(np.multiply(CrFraction_Fitted,100), 1e9*ExperimentalData_HSFL, yerr=1e9*ExperimentalData_HSFL_error, fmt='r^', label=r'Period HSFL')
    #ax1.set_yscale('log')
    plt.ylim((0.,1.1e9*wavelength))
    #ax12 = ax1.twinx()
    #plt.ylabel(r'SPP mean free path $L_{SPP}$ (m)')
    #plt.ylabel(r'Re($\varepsilon$), Im($\varepsilon$)')
    
    plot110, = ax1.semilogx(Te0, np.multiply(1e9,period0), 'r+', label=r'SPP period $\Lambda$, branch (-,-)')
    plot111, = ax1.semilogx(Te1, np.multiply(1e9,period1), 'k+', label=r'SPP period $\Lambda$, branch (-,+)')
    plot112, = ax1.semilogx(Te2, np.multiply(1e9,period2), 'b+', label=r'SPP period $\Lambda$, branch (+,-)')
    plot113, = ax1.semilogx(Te3, np.multiply(1e9,period3), 'go', label=r'SPP period $\Lambda$, branch ( +,+)')
    
    plot12, = ax1.semilogx(Te, np.multiply(1e9,wavelength*np.ones(np.shape(Te))), 'k-', linewidth=0.5, label=r'Laser wavelength $\lambda$')
    
    #plot14, = ax12.plot(Te_s, np.real(epsilonFilm_s), 'b+', label=r'Re($\varepsilon$)')
    #plot15, = ax12.plot(Te_s, np.imag(epsilonFilm_s), 'b^', label=r'Im($\varepsilon$)')
    
    #ax1.yaxis.label.set_color(plot110.get_color()) #colorizes the label
    #ax1.spines["left"].set_edgecolor(plot110.get_color()) #colorizes the axis
    #ax1.tick_params(axis='y', colors=plot110.get_color()) #colorizes the tics and numbers
    plotComb1= [plot110, plot111, plot112, plot113, plot12]; 
    
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
        ax2.set_xlabel(r'Te (K)')
        #plot21, = ax2.semilogy(Te_s, lspp_s,   'r^', label=r'SPP decay length $L_{SPP}$')
        plot21, = ax2.loglog(Te0, np.abs(lspp0),   'r^', label=r'SPP decay length $L_{SPP}$, --')
        plot22, = ax2.loglog(Te1, np.abs(lspp1),   'k^', label=r'SPP decay length $L_{SPP}$, -+')
        plot23, = ax2.loglog(Te2, np.abs(lspp2),   'b^', label=r'SPP decay length $L_{SPP}$, +-')
        plot24, = ax2.loglog(Te3, np.abs(lspp3),   'g^', label=r'SPP decay length $L_{SPP}$, ++')
        plot2   = [plot21, plot22, plot23, plot24] #, plot24]
        plotComb2 = plot2
        plt.tight_layout()
        
        #ax22 = ax2.twinx()
    if(PlotEpsilons and not PlotLspp): #plot in plotB window, keep ax3 name. 
        ax3 = plt.subplot(plotB)
    elif(PlotEpsilons and PlotLspp):
        ax3 = plt.subplot(plotC)
    if(PlotEpsilons):
        plt.ylabel(r'Re($\varepsilon$), Im($\varepsilon$)')
        plot31, = ax3.semilogx(Te, np.real(epsilonFilm), 'b+', label=r'Re$(\varepsilon(T_e))$')
        plot32, = ax3.semilogx(Te, np.imag(epsilonFilm), 'b^', label=r'Im$(\varepsilon(T_e))$')
        plot33, = ax3.semilogx(Te, np.multiply(epsBK7, np.ones(np.shape(Te))), 'k-', label=r'$Re[\varepsilon$(BK7)] ')
        plot3   = [plot31, plot32, plot33]
        plotComb3 = plot3
        plt.tight_layout()
        
        #ax2.yaxis.label.set_color(plot21.get_color()) #colorizes the label
        #ax2.spines["left"].set_edgecolor(plot21.get_color()) #colorizes the axis
        #ax2.tick_params(axis='y', colors=plot21.get_color()) #colorizes the tics and numbers
        
        #ax22.yaxis.label.set_color(plot22.get_color()) #colorizes the label
        #ax22.spines["right"].set_edgecolor(plot22.get_color()) #colorizes the axis
        #ax22.tick_params(axis='y', colors=plot22.get_color()) #colorizes the tics and numbers
        #labels2 = [l.get_label() for l in plot2]
        #ax2.legend(plot2, labels2, loc='best')
        
    if(PlotLspp and PlotEpsilons):
        ax3.set_xlabel('Te (K)')
        labelsComb2 = [l.get_label() for l in plotComb2]
        labelsComb3 = [l.get_label() for l in plotComb3]
        ax2.legend(plotComb2, labelsComb2, bbox_to_anchor=(1.04,1), loc="upper left", mode="expand")
        ax3.legend(plotComb3, labelsComb3, bbox_to_anchor=(1.04,1), loc="upper left", mode="expand")
    elif(PlotLspp and not PlotEpsilons):
        ax2.set_xlabel('Te (K)')
        labelsComb2 = [l.get_label() for l in plotComb2]
        ax2.legend(plotComb2, labelsComb2, bbox_to_anchor=(1.04,1), loc="upper left", mode="expand")
    elif(not PlotLspp and PlotEpsilons): #ax2 does not exist
        ax3.set_xlabel('Te (K)')
        labelsComb3 = [l.get_label() for l in plotComb3]
        ax3.legend(plotComb3, labelsComb3, bbox_to_anchor=(1.04,1), loc="upper left", mode="expand")
    else:
        ax1.set_xlabel('Te (K)')
    
    labelsComb1 = [l.get_label() for l in plotComb1]
    ax1.legend(plotComb1, labelsComb1, bbox_to_anchor=(1.04,1), loc="upper left", mode="expand")
    
    
    #plt.xlim((0,100))
    plt.tight_layout()
    filename="Dostovalov-Cr-Thickness-"+str(1E9*thickness_max)+"nm-TemperaturePeriod"
    plt.savefig(filename+".eps")
    plt.savefig(filename+".png")
    plt.show()
#}}}

    
# Scenario proposed by Thibault: an oxide layer grows at the top of the Cr sample, reducing progressively the periodicity by lambda/n. 
# This function can also be used to generate simple results for verification
def Burke_SymmetricModes(thickness_size): 
    wavelength=633e-9
    numberofroots = 10
    # Medium 1: thin film.
    eps1 = -19.+0.53j     #thin film
    # Medium 2: substrate. 
    eps2 = 4. #3.999999+0.004j
    # Medium 3: environment
    eps3 = 1.5**2 # eps2 #eps2: symmetric modes       #environment | substrate
    # Note: Inverting eps2 and eps3 should have no effect on the possible modes, but only on field amplification. 
    
    # Conversion to beta/k0: 
    def PeriodToBetaNorm(period): 
        # period = 2.*np.pi / beta.real
        k0 = 2.*np.pi / wavelength
        beta_norm_re = np.divide(np.divide(2.*np.pi, period), k0)
        return beta_norm_re
    
    #thickness_size = 20
    thickness_min  = 10e-9
    thickness_max  = 100e-9
    t_list = np.linspace(thickness_min, thickness_max, thickness_size, endpoint=True)
    
    summary = np.zeros((0, 7))
    for thickness in t_list:
        roots = findroots(eps1, eps2, eps3,
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
            print("Roots in branch #"+str(branch))
            for order in np.arange(0,len(roots_in_branch)):
                roots_in_branch_order = roots_in_branch[order]
                #print("\n")
                fraction = 0e0 #irrelevant in this context
                print("Thickness:"+str(thickness)+", Order #"+str(order)+": Period="+str(roots_in_branch_order[0])+" Lspp="+str(roots_in_branch_order[1]))
                print("beta/k0="+str(PeriodToBetaNorm(roots_in_branch_order[0])))
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
    
    summary_branch0_export = np.array([np.real(thickness0), np.real(fractionOxide0), np.real(epsilonFilm0), np.imag(epsilonFilm0), np.real(branch0), np.real(root_number0), np.real(period0)])
    summary_branch1_export = np.array([np.real(thickness1), np.real(fractionOxide1), np.real(epsilonFilm1), np.imag(epsilonFilm1), np.real(branch1), np.real(root_number1), np.real(period1)])
    summary_branch2_export = np.array([np.real(thickness2), np.real(fractionOxide2), np.real(epsilonFilm2), np.imag(epsilonFilm2), np.real(branch2), np.real(root_number2), np.real(period2)])
    summary_branch3_export = np.array([np.real(thickness3), np.real(fractionOxide3), np.real(epsilonFilm3), np.imag(epsilonFilm3), np.real(branch3), np.real(root_number3), np.real(period3)])
    
    summary_branch0_export_t = np.transpose(summary_branch0_export)
    summary_branch1_export_t = np.transpose(summary_branch1_export)
    summary_branch2_export_t = np.transpose(summary_branch2_export)
    summary_branch3_export_t = np.transpose(summary_branch3_export)
        
    print("SPP branches are ready. Exporting to CSV...")
    filename = "Burke-SPPmodes-branch"
    np.savetxt(filename+"0"+".csv", summary_branch0_export_t)
    np.savetxt(filename+"1"+".csv", summary_branch1_export_t)
    np.savetxt(filename+"2"+".csv", summary_branch2_export_t)
    np.savetxt(filename+"3"+".csv", summary_branch3_export_t)
    
    
    plt.figure()
    ax1 = plt.subplot(111)
    plt.xlabel(r'Thickness (nm)')
    #plt.ylabel(r'SPP period $\Lambda$ (nm)') 
    #ax1.set_yscale('log')
    #plt.ylim((0.,1.1e9*wavelength))
    #plot110, = ax1.plot(np.multiply(1e9,thickness0), np.multiply(1e9,period0), 'r+', label=r'SPP period $\Lambda$, branch (-,-)')
    #plot111, = ax1.plot(np.multiply(1e9,thickness1), np.multiply(1e9,period1), 'k+', label=r'SPP period $\Lambda$, branch (-,+)')
    #plot112, = ax1.plot(np.multiply(1e9,thickness2), np.multiply(1e9,period2), 'b+', label=r'SPP period $\Lambda$, branch (+,-)')
    #plot113, = ax1.plot(np.multiply(1e9,thickness3), np.multiply(1e9,period3), 'g+', label=r'SPP period $\Lambda$, branch ( +,+)')
    
    #plot12,  = ax1.plot(np.multiply(1e9, thickness), np.multiply(1e9,wavelength*np.ones(np.shape(fractionOxide))), 'k-', linewidth=0.5, label=r'Laser wavelength $\lambda$')
    
    #ax12 = ax1.twinx()
    #plot14, = ax12.plot(fractionOxide_s, np.real(epsilonFilm_s), 'b+', label=r'Re($\varepsilon$)')
    #plot15, = ax12.plot(fractionOxide_s, np.imag(epsilonFilm_s), 'b^', label=r'Im($\varepsilon$)')
    plot14,  = ax1.plot(np.multiply(1e9,thickness0), np.multiply(1e0,PeriodToBetaNorm(period0)), 'r^', label=r'$\beta/k_0$, (-,-)')
    plot15,  = ax1.plot(np.multiply(1e9,thickness1), np.multiply(1e0,PeriodToBetaNorm(period1)), 'k^', label=r'$\beta/k_0$, (-,+)')
    plot16,  = ax1.plot(np.multiply(1e9,thickness2), np.multiply(1e0,PeriodToBetaNorm(period2)), 'b^', label=r'$\beta/k_0$, (+,-)')
    plot17,  = ax1.plot(np.multiply(1e9,thickness3), np.multiply(1e0,PeriodToBetaNorm(period3)), 'g^', label=r'$\beta/k_0$, (+,+)')
    
    #plt.ylabel(r'SPP mean free path $L_{SPP}$ (m)')
    #plt.ylabel(r'Re($\varepsilon$), Im($\varepsilon$)')
    plt.ylabel(r'$\beta/k_0$')
    #ax1.yaxis.label.set_color(plot110.get_color()) #colorizes the label
    #ax1.spines["left"].set_edgecolor(plot110.get_color()) #colorizes the axis
    #ax1.tick_params(axis='y', colors=plot110.get_color()) #colorizes the tics and numbers
    plotComb1 = []
    #plotComb1+= [plot110, plot111, plot112, plot113, plot12]; 
    plotComb1+=[plot14, plot15, plot16, plot17]
    
    #ax12.yaxis.label.set_color(plot14.get_color()) #colorizes the label
    #ax12.spines["right"].set_edgecolor(plot14.get_color()) #colorizes the axis
    #ax12.tick_params(axis='y', colors=plot14.get_color()) #colorizes the tics and numbers
    plt.tight_layout()
    labelsComb1 = [l.get_label() for l in plotComb1]
    ax1.legend(plotComb1, labelsComb1, loc="upper right")
    plt.show()


## Compute the spatial period and mean-free path of Surface Plasmon Polaritons 
# at the interface between three materials. 
# Here we investigate the role of the materials mixing happening pulse by pulse, oxygen from ambient air diffuses into the Cr, and leads to formation of Cr2O3. 
# By making use of a Lorenz-Lorentz model, one can construct a dielectric permittivity of the oxide
# and therefore observe the transition between SPP and waveguiding modes. 
# Note: Some results exist showing that L_spp < 0. These modes may be physical, and should be discussed 
# in the frame of the works of P. Berini. 
def ScenarioOfCrOxideMixture(epsSample, epsOxide=1., epsSubstrate=1., epsEnvironment=1., fraction_size=10, SampleName='Cr', OxideName='Cr2O3', SubstrateName='BK7', PlotLspp=False, FilterNegativeLspp=False, PlotEpsilons=True):
    print(Header, "# Info: Considering a mixed fraction of Cr with Cr2O3 with several thicknesses.")
    wavelength = 1026e-9
    #fraction_size = 60
    numberofroots = 10 #per branch. 10 exceeds the final number of roots per branch
    
    fraction_min = 0.
    fraction_max = 1.

    thickness_size = 1
    thickness_min  = 70e-9 #28e-9
    thickness_max  = thickness_min
    
    if(PlotLspp and PlotEpsilons):
        plotA = 311; plotB=312; plotC=313
    elif((PlotLspp and not PlotEpsilons) or (not PlotLspp and PlotEpsilons)): 
        plotA=211; plotB=212
    else: 
        plotA=111; 
        
    ## Running 
    fraction = np.linspace(fraction_min, fraction_max, fraction_size) #fraction of CrO2 (includes the final value)
    epsCrCr2O3_list = MaxwellGarnett2(epsSample, epsOxide, fraction) #fraction refers to the fraction oxide here! 
    #NOTE: epsOxide here is Cr oxide, not environment. 
    print(Header, "# Info: size of the fraction matrix: ", fraction_size)

    print(Header, "# Info: preparation of the root finder.")
    
    #results = Parallel(n_jobs=num_cores)(delayed(processInput)(i) for i in inputs)

    #t = 100E-9 #thickness of the layer in meters
    #thickness = 10e-9
    print("\n")
    print("# Fraction of CrO2: ", fraction)
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
        
        ### This version is general enough for number_of_roots > 1. 
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
                if(roots_in_branch[order][0] > 1E-15 and abs(roots_in_branch[order][1]) > 1E-10): 
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
    
    print(Header+"** Preparation of the plots as function of oxide ratio")

    ExperimentalData_velocity   = np.array([1e-6, 10e-16, 50e-6, 100e-6, 200e-6, 300e-6]) #m/s
    ExperimentalData_LSFL       = np.array([696e-9, 704e-9, 816e-9, 858e-9, -100e0, -100e-9]) #m #better observed for low velocities, i.e., high number of pulses, i.e., largest amounts of oxide
    ExperimentalData_LSFL_error = np.array([78e-9,71e-9,139e-9,140e-9,0e-9,0e-9]) #m
    ExperimentalData_HSFL       = np.array([170e-9, 159e-9, 217e-9, 244e-9, 249e-9, 238e-9]) #m
    ExperimentalData_HSFL_error = np.array([64e-9,38e-9,101e-9,110e-9,128e-9,72.5e-9]) #m
    CrFraction_Fitted           = np.array([0.6e0, 0.7e0, 0.8e0, 0.86e0, 0.90e0, 0.95e0]) #hand fitted to match period with existing modes [on request of Nadya]
    CrO2Fraction_Fitted         = np.add(1, - np.array([0.6e0, 0.7e0, 0.8e0, 0.86e0, 0.90e0, 0.95e0])) #we express it in term of oxidized
    
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
    
    #ReBeta = np.divide(2.*np.pi,period_s)
    #ImBeta = np.divide(0.5,lspp_s)
    
    # Then we could plot them in the right order
    thickness, fractionOxide, epsilonFilm, branch, root_number, period, lspp = SplitSummaryTable(summary)
    
    thickness0, fractionOxide0, epsilonFilm0, branch0, root_number0, period0, lspp0 = SplitSummaryTable(summary_branch0)
    thickness1, fractionOxide1, epsilonFilm1, branch1, root_number1, period1, lspp1 = SplitSummaryTable(summary_branch1)
    thickness2, fractionOxide2, epsilonFilm2, branch2, root_number2, period2, lspp2 = SplitSummaryTable(summary_branch2)
    thickness3, fractionOxide3, epsilonFilm3, branch3, root_number3, period3, lspp3 = SplitSummaryTable(summary_branch3)
    
    summary_export         = np.array([np.real(fractionOxide) , np.real(epsilonFilm) , np.imag(epsilonFilm) , np.real(branch) , np.real(root_number) , np.real(period) ])
    summary_branch0_export = np.array([np.real(fractionOxide0), np.real(epsilonFilm0), np.imag(epsilonFilm0), np.real(branch0), np.real(root_number0), np.real(period0)])
    summary_branch1_export = np.array([np.real(fractionOxide1), np.real(epsilonFilm1), np.imag(epsilonFilm1), np.real(branch1), np.real(root_number1), np.real(period1)])
    summary_branch2_export = np.array([np.real(fractionOxide2), np.real(epsilonFilm2), np.imag(epsilonFilm2), np.real(branch2), np.real(root_number2), np.real(period2)])
    summary_branch3_export = np.array([np.real(fractionOxide3), np.real(epsilonFilm3), np.imag(epsilonFilm3), np.real(branch3), np.real(root_number3), np.real(period3)])
    
    summary_export_t         = np.transpose(summary_export)
    summary_branch0_export_t = np.transpose(summary_branch0_export)
    summary_branch1_export_t = np.transpose(summary_branch1_export)
    summary_branch2_export_t = np.transpose(summary_branch2_export)
    summary_branch3_export_t = np.transpose(summary_branch3_export)

    # NOTE: I would like to get clean numbers, not the content of summary_branch0
    # Is there a problem with root_number0 for example? 
    print("SPP branches are ready. Exporting to CSV...")
    filename = "Dostovalov-SPPmodes-Cr-CrO2"
    np.savetxt(filename+".csv", summary_export_t)
    np.savetxt(filename+"-branch0"+".csv", summary_branch0_export_t)
    np.savetxt(filename+"-branch1"+".csv", summary_branch1_export_t)
    np.savetxt(filename+"-branch2"+".csv", summary_branch2_export_t)
    np.savetxt(filename+"-branch3"+".csv", summary_branch3_export_t)
    
    plt.figure()
    ax1 = plt.subplot(plotA)
    #plt.title(r'Film thickness $t=$'+str(round(np.real(thickness[0])*1e9))+' nm')
    #plt.xlabel(r'Fraction of Cr oxide  (perc.)')
    plt.ylabel(r'SPP period $\Lambda$ (nm)') 
    #ax1y = ax1.twiny()
    #plot1y1 = ax1.errorbar(np.multiply(CrFraction_Fitted,100), 1e9*ExperimentalData_LSFL, yerr=1e9*ExperimentalData_LSFL_error, fmt='ro', label=r'Period LSFL')
    #plot1y2 = ax1.errorbar(np.multiply(CrFraction_Fitted,100), 1e9*ExperimentalData_HSFL, yerr=1e9*ExperimentalData_HSFL_error, fmt='r^', label=r'Period HSFL')
    #ax1.set_yscale('log')
    plt.ylim((0.,1.1e9*wavelength))
    #ax12 = ax1.twinx()
    #plt.ylabel(r'SPP mean free path $L_{SPP}$ (m)')
    #plt.ylabel(r'Re($\varepsilon$), Im($\varepsilon$)')
    
    plot110, = ax1.plot(np.multiply(100,fractionOxide0), np.multiply(1e9,period0), 'r+', label=r'SPP period $\Lambda$, branch (-,-)')
    plot111, = ax1.plot(np.multiply(100,fractionOxide1), np.multiply(1e9,period1), 'k+', label=r'SPP period $\Lambda$, branch (-,+)')
    plot112, = ax1.plot(np.multiply(100,fractionOxide2), np.multiply(1e9,period2), 'b+', label=r'SPP period $\Lambda$, branch (+,-)')
    plot113, = ax1.plot(np.multiply(100,fractionOxide3), np.multiply(1e9,period3), 'go', label=r'SPP period $\Lambda$, branch ( +,+)')
    
    plot12, = ax1.plot(np.multiply(100, fractionOxide), np.multiply(1e9,wavelength*np.ones(np.shape(fractionOxide))), 'k-', linewidth=0.5, label=r'Laser wavelength $\lambda$')
    
    #plot14, = ax12.plot(fractionOxide_s, np.real(epsilonFilm_s), 'b+', label=r'Re($\varepsilon$)')
    #plot15, = ax12.plot(fractionOxide_s, np.imag(epsilonFilm_s), 'b^', label=r'Im($\varepsilon$)')
    
    #ax1.yaxis.label.set_color(plot110.get_color()) #colorizes the label
    #ax1.spines["left"].set_edgecolor(plot110.get_color()) #colorizes the axis
    #ax1.tick_params(axis='y', colors=plot110.get_color()) #colorizes the tics and numbers
    plotComb1= [plot110, plot111, plot112, plot113, plot12]; 
    
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
        ax2.set_xlabel(r'Fraction of Cr (%)')
        #plot21, = ax2.semilogy(fractionOxide_s, lspp_s,   'r^', label=r'SPP decay length $L_{SPP}$')
        plot21, = ax2.semilogy(np.multiply(1E2,fractionOxide0), np.abs(lspp0),   'r^', label=r'SPP decay length $L_{SPP}$, --')
        plot22, = ax2.semilogy(np.multiply(1E2,fractionOxide1), np.abs(lspp1),   'k^', label=r'SPP decay length $L_{SPP}$, -+')
        plot23, = ax2.semilogy(np.multiply(1E2,fractionOxide2), np.abs(lspp2),   'b^', label=r'SPP decay length $L_{SPP}$, +-')
        plot24, = ax2.semilogy(np.multiply(1E2,fractionOxide3), np.abs(lspp3),   'g^', label=r'SPP decay length $L_{SPP}$, ++')
        plot2   = [plot21, plot22, plot23, plot24] #, plot24]
        plotComb2 = plot2
        plt.tight_layout()
        
        #ax22 = ax2.twinx()
    if(PlotEpsilons and not PlotLspp): #plot in plotB window, keep ax3 name. 
        ax3 = plt.subplot(plotB)
    if(PlotEpsilons):
        plt.ylabel(r'Re($\varepsilon$), Im($\varepsilon$)')
        plot31, = ax3.plot(np.multiply(1E2,fractionOxide), np.real(epsilonFilm), 'b+', label=r'Re$(\varepsilon)$ (Cr + oxide)')
        plot32, = ax3.plot(np.multiply(1E2,fractionOxide), np.imag(epsilonFilm), 'b^', label=r'Im$(\varepsilon)$ (Cr + oxide)')
        plot33, = ax3.plot(np.multiply(1E2,fractionOxide), np.multiply(epsBK7, np.ones(np.shape(fractionOxide))), 'k-', label=r'$Re[\varepsilon$(BK7)] ')
        plot3   = [plot31, plot32, plot33]
        plotComb3 = plot3
        plt.tight_layout()
        
        #ax2.yaxis.label.set_color(plot21.get_color()) #colorizes the label
        #ax2.spines["left"].set_edgecolor(plot21.get_color()) #colorizes the axis
        #ax2.tick_params(axis='y', colors=plot21.get_color()) #colorizes the tics and numbers
        
        #ax22.yaxis.label.set_color(plot22.get_color()) #colorizes the label
        #ax22.spines["right"].set_edgecolor(plot22.get_color()) #colorizes the axis
        #ax22.tick_params(axis='y', colors=plot22.get_color()) #colorizes the tics and numbers
        #labels2 = [l.get_label() for l in plot2]
        #ax2.legend(plot2, labels2, loc='best')
        
    if(PlotLspp and PlotEpsilons):
        ax3.set_xlabel('Fraction of Cr oxide (%)')
        labelsComb2 = [l.get_label() for l in plotComb2]
        labelsComb3 = [l.get_label() for l in plotComb3]
        ax2.legend(plotComb2, labelsComb2, bbox_to_anchor=(1.04,1), loc="upper left", mode="expand")
        ax3.legend(plotComb3, labelsComb3, bbox_to_anchor=(1.04,1), loc="upper left", mode="expand")
    elif(PlotLspp and not PlotEpsilons):
        ax2.set_xlabel('Fraction of Cr oxide (%)')
        labelsComb2 = [l.get_label() for l in plotComb2]
        ax2.legend(plotComb2, labelsComb2, bbox_to_anchor=(1.04,1), loc="upper left", mode="expand")
    elif(not PlotLspp and PlotEpsilons): #ax2 does not exist
        ax3.set_xlabel('Fraction of Cr oxide (%)')
        labelsComb3 = [l.get_label() for l in plotComb3]
        ax3.legend(plotComb3, labelsComb3, bbox_to_anchor=(1.04,1), loc="upper left", mode="expand")
    else:
        ax1.set_xlabel('Fraction of Cr oxide (%)')
    
    labelsComb1 = [l.get_label() for l in plotComb1]
    ax1.legend(plotComb1, labelsComb1, bbox_to_anchor=(1.04,1), loc="upper left", mode="expand")
    
    
    #plt.xlim((0,100))
    plt.tight_layout()
    filename="Dostovalov-"+OxideName+"-mixedWith-"+SampleName+"-Thickness-"+str(1E9*thickness_max)+"nm"
    plt.savefig(filename+".eps")
    plt.savefig(filename+".png")
    plt.show()
#}}}

## Compute the spatial period and mean-free path of Surface Plasmon Polaritons 
# at the interface between three materials. 
# Here we investigate the role of the materials mixing happening pulse by pulse, oxygen from ambient air diffuses into the Cr, and leads to formation of Cr2O3/CrO2. 
# By making use of a Lorenz-Lorentz model, one can construct a dielectric permittivity of the oxide
# and therefore observe the transition between SPP and waveguiding modes. 
# Note: Some results exist showing that L_spp < 0. These modes may be physical, and should be discussed 
# in the frame of the works of P. Berini. 
def ScenarioOfCrOxideMixture3(epsSample, epsOxide1=1., epsOxide2=1., epsSubstrate=1., fraction_size=60, SampleName='Cr', OxideName1='Cr2O3', OxideName2='CrO2', SubstrateName='BK7', PlotLspp=False, FilterNegativeLspp=False, PlotEpsilons=True, PlotExperimentalData=True):
    print(Header, "# Info: Considering a mixed fraction of Cr with Cr2O3 and CrO2 with several thicknesses.")
    wavelength = 1026e-9
    #fraction_size = 60
    numberofroots = 10
    
    number_of_ticks = 5
    
    fraction_min = 0.
    fraction_max = 1.

    scanning_velocity_min = 0e0
    scanning_velocity_max = 350e-6 #m/s
    
    thickness_size = 1
    thickness_min  = 28e-9
    thickness_max  = 28e-9
    
    #PlotLspp = False
    
    if(PlotLspp and PlotEpsilons):
        plotA = 311; plotB=312; plotC=313
    elif((PlotLspp and not PlotEpsilons) or (not PlotLspp and PlotEpsilons)): 
        plotA=211; plotB=212
    else: 
        plotA=111; 

    print(Header, "The fraction of Cr2O3 over CrO2 is fixed by Raman measurements.")
    # TODO: we hereby prepare the fraction of Cr2O3 / CrO2 from experimental data
    # FractionA is Cr
    # FractionB is Cr2O3
    # FractionC is CrO2
    # A+B+C = 1
    # B/C = f(scanning_velocity) [experimental data provided by A. Dostovalov]
    # Solution: A = 1-B-C. B = f(v)*C, C = parameter to plot with. 
    
    # Generating the final velocity mesh
    scanning_velocities        = np.linspace(scanning_velocity_min, scanning_velocity_max, fraction_size)
    scanning_velocities_xticks = np.linspace(scanning_velocity_min, scanning_velocity_max, number_of_ticks)
    print(Header, "Acquisition of the f(scanning_velocity)")
    Cr2O3overCrO2_ratio_min = np.loadtxt("Cr-Cr2O3-CrO2/RatioCr2O3overCrO2-min.csv", delimiter="\t", skiprows=1)
    Cr2O3overCrO2_ratio_max = np.loadtxt("Cr-Cr2O3-CrO2/RatioCr2O3overCrO2-min.csv", delimiter="\t", skiprows=1)
    
    v_exp_min     = Cr2O3overCrO2_ratio_min[:,0] * 1E-6 #m/s
    v_exp_max     = Cr2O3overCrO2_ratio_max[:,0] * 1E-6 #m/s
    ratio_exp_min = Cr2O3overCrO2_ratio_min[:,1]
    ratio_exp_max = Cr2O3overCrO2_ratio_max[:,1]
    
    print(Header, "Generating the function interpolation of the f(v) data.") #f_Cr2O3overCrO2_ratio_min are functions
    f_Cr2O3overCrO2_ratio_min = InterpolatedUnivariateSpline(v_exp_min, ratio_exp_min, k=1) #other orders dont work
    f_Cr2O3overCrO2_ratio_max = InterpolatedUnivariateSpline(v_exp_max, ratio_exp_max, k=1)
    
    v_exp_av = 0.5 * (v_exp_min + v_exp_max)
    ratio_exp_av = 0.5*(ratio_exp_min + ratio_exp_max)
    ratio_exp_err = 0.5*(ratio_exp_max - ratio_exp_min)
    print(Header, "Build the average of interpolated functions") #(f_Cr2O3overCrO2_ratio_av is an array)
    f_Cr2O3overCrO2_ratio_av = 0.5*(f_Cr2O3overCrO2_ratio_min(scanning_velocities) + f_Cr2O3overCrO2_ratio_max(scanning_velocities))
    
    print(Header, "Conclusion of the materials fractions.")
    fractionC = np.linspace(fraction_min, fraction_max, fraction_size) #fraction of CrO2 drags all the [0:1] interval
    fractionB = f_Cr2O3overCrO2_ratio_av * fractionC #associating fraction of B with f(v)
    fractionA = 1E0 - f_Cr2O3overCrO2_ratio_av * fractionC - fractionC #Cr fraction
    
    fraction_total          = fractionA+fractionB+fractionC
    fraction_total_expected = np.ones(np.shape(fractionA))
    
    print(Header, "Test: fractionA+fractionB+fractionC: min, max", np.min(fraction_total), np.max(fraction_total))
    print(Header, "Test: is it equal to expected value everywhere?", np.array_equal(fraction_total, fraction_total_expected)) 
    print(Header, "Test: max of the differences: ", np.max(np.abs(np.add(fraction_total, np.multiply(-1,fraction_total_expected)))))
    
    
    
    print(Header, "Applying the MaxwellGarnett3 routine...")
    epsCrCr2O3CrO2_list = MaxwellGarnett3(epsSample, epsOxide1, epsOxide2, fractionA, fractionB, fractionC)
    
    print(Header, "# Info: size of the fraction matrix: ", fraction_size)

    print(Header, "# Info: preparation of the root finder.")
    fraction = fractionC #we will plot using fractionCrO2 now. 
    
    #results = Parallel(n_jobs=num_cores)(delayed(processInput)(i) for i in inputs)

    #t = 100E-9 #thickness of the layer in meters
    #thickness = 10e-9
    print("\n")
    print("# Fraction of CrO2: ", fraction)
    summary = np.zeros((0, 7))
    ## Preparation of the thin film modeling for various compositions
    for fraction_index in np.arange(0,fraction_size): #arange excludes the last one, linspace includes it
        # Medium 1: thin film. 
        eps1 = epsCrCr2O3CrO2_list[fraction_index]        #thin film
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
        print("\nCrO2 oxide fraction #"+str(fraction_index)+"="+str(fraction[fraction_index]))
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
    
    print(Header+"** Preparation of the plots as function of oxide ratio")

    ExperimentalData_velocity   = np.array([1e-6, 10e-16, 50e-6, 100e-6, 200e-6, 300e-6]) #m/s
    #NOTE: with increasing the scanning velocity, the period of LIPSS increases
    #NOTE: with increasing the scanning velocity, the ratio Cr2O3/CrO2 increases, Cr2O3 still dominates
    #NOTE: hence, logically, the period of LIPSS should increase with increasing Cr2O3 ratio (or equivalently decreasing CrO2 ratio), i.e., Period decreases if increasing CrO2 ratio.
    # Therefore, d(period)/d(CrO2 ratio) < 0. <CrO2Fraction_Fitted> should match with this. 
    #NOTE: the content of <scanning_velocities> variable should be much more precise. 
    ExperimentalData_LSFL       = np.array([696e-9, 704e-9, 816e-9, 858e-9, -100e0, -100e-9]) #m #better observed for low velocities, i.e., high number of pulses, i.e., largest amounts of oxide
    ExperimentalData_LSFL_error = np.array([78e-9,71e-9,139e-9,140e-9,0e-9,0e-9]) #m
    ExperimentalData_HSFL       = np.array([170e-9, 159e-9, 217e-9, 244e-9, 249e-9, 238e-9]) #m
    ExperimentalData_HSFL_error = np.array([64e-9,38e-9,101e-9,110e-9,128e-9,72.5e-9]) #m
    #CrO2Fraction_Fitted         = np.array([0.6e0, 0.7e0, 0.8e0, 0.86e0, 0.90e0, 0.95e0]) #there is now more Cr2O3 when increasing scanning speed 
    CrO2Fraction_Fitted         = 0.5*np.array([0.6, 0.5, 0.4, 0.3, 0.2, 0.1]) #there is now more Cr2O3 when increasing scanning speed
    # To compute the resulting fraction of Cr2O3, just pass CrO2Fraction_Fitted to f(scanning_velocity). 
    Cr2O3overCrO2_av     = 0.5*(f_Cr2O3overCrO2_ratio_min(v_exp_min)+f_Cr2O3overCrO2_ratio_max(v_exp_max)) #experimental ratio f(v) = Cr2O3/CrO2
    #print(Header, "Test: Cr2O3overCrO2_av: ", np.min(Cr2O3overCrO2_av), np.max(Cr2O3overCrO2_av)) #Passed
    #Cr2O3Fraction_Fitted  = Cr2O3overCrO2_av * CrO2Fraction_Fitted #This is the ratio of Cr2O3 corresponding to ratio of CrO2 that is plotting on the other X axis. WARNING: nothing tells that Python plots the X2 axis at the position that corresponds!!!
    
    #x1ticks  = np.linspace(0, 1, number_of_ticks) #graduation for C (CrO2 concentration)
    #generated a mesh of velocities with same dimension as x1ticks: scanning_velocities_xticks
    # We build X2 axis using definition of B = C * f(v). 
    #x12ticks = np.round(x1ticks * 0.5*(f_Cr2O3overCrO2_ratio_min(scanning_velocities_xticks)+f_Cr2O3overCrO2_ratio_max(scanning_velocities_xticks)), decimals=0) #BUG: does not work for some reason
    
    
    
    #print(summary)

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
    
    #ReBeta = np.divide(2.*np.pi,period_s)
    #ImBeta = np.divide(0.5,lspp_s)
    
    # Then we could plot them in the right order
    #thickness, branch, root_number, period, lspp, fractionOxide, epsilonFilm = SplitSummaryTable(summary)
    
    thickness, fractionOxide, epsilonFilm, branch, root_number, period, lspp = SplitSummaryTable(summary)
    
    thickness0, fractionOxide0, epsilonFilm0, branch0, root_number0, period0, lspp0 = SplitSummaryTable(summary_branch0)
    thickness1, fractionOxide1, epsilonFilm1, branch1, root_number1, period1, lspp1 = SplitSummaryTable(summary_branch1)
    thickness2, fractionOxide2, epsilonFilm2, branch2, root_number2, period2, lspp2 = SplitSummaryTable(summary_branch2)
    thickness3, fractionOxide3, epsilonFilm3, branch3, root_number3, period3, lspp3 = SplitSummaryTable(summary_branch3)
    
    summary_branch0_export = np.array([np.real(fractionOxide0), np.real(epsilonFilm0), np.imag(epsilonFilm0), np.real(branch0), np.real(root_number0), np.real(period0)])
    summary_branch1_export = np.array([np.real(fractionOxide1), np.real(epsilonFilm1), np.imag(epsilonFilm1), np.real(branch1), np.real(root_number1), np.real(period1)])
    summary_branch2_export = np.array([np.real(fractionOxide2), np.real(epsilonFilm2), np.imag(epsilonFilm2), np.real(branch2), np.real(root_number2), np.real(period2)])
    summary_branch3_export = np.array([np.real(fractionOxide3), np.real(epsilonFilm3), np.imag(epsilonFilm3), np.real(branch3), np.real(root_number3), np.real(period3)])
    
    summary_branch0_export_t = np.transpose(summary_branch0_export)
    summary_branch1_export_t = np.transpose(summary_branch1_export)
    summary_branch2_export_t = np.transpose(summary_branch2_export)
    summary_branch3_export_t = np.transpose(summary_branch3_export)
        
    print("SPP branches are ready. Exporting to CSV...")
    filename = "Dostovalov-SPPmodes-branch"
    np.savetxt(filename+"0"+".csv", summary_branch0)
    np.savetxt(filename+"1"+".csv", summary_branch1)
    np.savetxt(filename+"2"+".csv", summary_branch2)
    np.savetxt(filename+"3"+".csv", summary_branch3)
    
    plt.figure()
    ax1 = plt.subplot(plotA)
    #ax1.set_xticks(np.multiply(100,x1ticks))
    #ax1.set_xticklabels(np.multiply(100,x1ticks))
    #plt.title(r'Film thickness $t=$'+str(round(np.real(thickness[0])*1e9))+' nm')
    #plt.xlabel(r'Fraction of Cr (perc.)')
    plt.ylabel(r'SPP period $\Lambda$ (nm)') 
    #ax1y = ax1.twiny()
    if(PlotExperimentalData):
        ## Manual entry of oxidation ratio
        plot1y1 = ax1.errorbar(np.multiply(CrO2Fraction_Fitted,100), 1e9*ExperimentalData_LSFL, yerr=1e9*ExperimentalData_LSFL_error, fmt='ro', label=r'Period LSFL')
        plot1y2 = ax1.errorbar(np.multiply(CrO2Fraction_Fitted,100), 1e9*ExperimentalData_HSFL, yerr=1e9*ExperimentalData_HSFL_error, fmt='r^', label=r'Period HSFL')
        # We can also plot the data using Fraction of Cr2O3 in top axis !
        #ax12 = ax1.twiny() #BUG: did not work as we don't have the mapping between scanning_velocities and FractionC !
        #ax12.set_xticks(np.multiply(100,x1ticks)) #position indicated as function of X-axis
        #ax12.set_xticklabels(np.multiply(100,x12ticks))
        ##ax12.errorbar(np.multiply(Cr2O3Fraction_Fitted,100), 1e9*ExperimentalData_LSFL, yerr=1e9*ExperimentalData_LSFL_error, fmt='go')
        #ax12.errorbar(np.multiply(Cr2O3Fraction_Fitted,100), 1e9*ExperimentalData_HSFL, yerr=1e9*ExperimentalData_HSFL_error, fmt='g^')
        #ax12.set_xlabel(r'Fraction of Cr$_2$O$_3$')
        # Automatic adjustement is not available, as FractionC is a PARAMETER of all other fractions. 
        #ax1.set_yscale('log')
    plt.ylim((0.,1.1e9*wavelength))
    #ax12 = ax1.twinx()
    #plt.ylabel(r'SPP mean free path $L_{SPP}$ (m)')
    #plt.ylabel(r'Re($\varepsilon$), Im($\varepsilon$)')
    
    plot110, = ax1.plot(np.multiply(100,fractionOxide0), np.multiply(1e9,period0), 'r+', label=r'SPP period $\Lambda$, branch (-,-)')
    plot111, = ax1.plot(np.multiply(100,fractionOxide1), np.multiply(1e9,period1), 'k+', label=r'SPP period $\Lambda$, branch (-,+)')
    plot112, = ax1.plot(np.multiply(100,fractionOxide2), np.multiply(1e9,period2), 'b+', label=r'SPP period $\Lambda$, branch (+,-)')
    plot113, = ax1.plot(np.multiply(100,fractionOxide3), np.multiply(1e9,period3), 'go', label=r'SPP period $\Lambda$, branch ( +,+)')
    
    plot12, = ax1.plot(np.multiply(100, fractionOxide), np.multiply(1e9,wavelength*np.ones(np.shape(fractionOxide))), 'k-', linewidth=0.5, label=r'Laser wavelength $\lambda$')
    
    #plot14, = ax12.plot(fractionOxide_s, np.real(epsilonFilm_s), 'b+', label=r'Re($\varepsilon$)')
    #plot15, = ax12.plot(fractionOxide_s, np.imag(epsilonFilm_s), 'b^', label=r'Im($\varepsilon$)')
    
    #ax1.yaxis.label.set_color(plot110.get_color()) #colorizes the label
    #ax1.spines["left"].set_edgecolor(plot110.get_color()) #colorizes the axis
    #ax1.tick_params(axis='y', colors=plot110.get_color()) #colorizes the tics and numbers
    plotComb1= [plot110, plot111, plot112, plot113, plot12]; 
    
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
        ax2.set_xlabel(r'Fraction of Cr (%)')
        #plot21, = ax2.semilogy(fractionOxide_s, lspp_s,   'r^', label=r'SPP decay length $L_{SPP}$')
        plot21, = ax2.semilogy(np.multiply(1E2,fractionOxide0), np.abs(lspp0),   'r^', label=r'SPP decay length $L_{SPP}$, --')
        plot22, = ax2.semilogy(np.multiply(1E2,fractionOxide1), np.abs(lspp1),   'k^', label=r'SPP decay length $L_{SPP}$, -+')
        plot23, = ax2.semilogy(np.multiply(1E2,fractionOxide2), np.abs(lspp2),   'b^', label=r'SPP decay length $L_{SPP}$, +-')
        plot24, = ax2.semilogy(np.multiply(1E2,fractionOxide3), np.abs(lspp3),   'g^', label=r'SPP decay length $L_{SPP}$, ++')
        plot2   = [plot21, plot22, plot23, plot24] #, plot24]
        plotComb2 = plot2
        plt.tight_layout()
        
        #ax22 = ax2.twinx()
    if(PlotEpsilons and not PlotLspp): #plot in plotB window, keep ax3 name. 
        ax3 = plt.subplot(plotB)
    if(PlotEpsilons):
        plt.ylabel(r'Re($\varepsilon$), Im($\varepsilon$)')
        plot31, = ax3.plot(np.multiply(1E2,fractionOxide), np.real(epsilonFilm), 'b+', label=r'Re$(\varepsilon)$ (Cr + oxides)')
        plot32, = ax3.plot(np.multiply(1E2,fractionOxide), np.imag(epsilonFilm), 'b^', label=r'Im$(\varepsilon)$ (Cr + oxides)')
        plot33, = ax3.plot(np.multiply(1E2,fractionOxide), np.multiply(epsBK7, np.ones(np.shape(fractionOxide))), 'k-', label=r'$Re[\varepsilon$(BK7)] ')
        plot3   = [plot31, plot32, plot33]
        plotComb3 = plot3
        plt.tight_layout()
        
        #ax2.yaxis.label.set_color(plot21.get_color()) #colorizes the label
        #ax2.spines["left"].set_edgecolor(plot21.get_color()) #colorizes the axis
        #ax2.tick_params(axis='y', colors=plot21.get_color()) #colorizes the tics and numbers
        
        #ax22.yaxis.label.set_color(plot22.get_color()) #colorizes the label
        #ax22.spines["right"].set_edgecolor(plot22.get_color()) #colorizes the axis
        #ax22.tick_params(axis='y', colors=plot22.get_color()) #colorizes the tics and numbers
        #labels2 = [l.get_label() for l in plot2]
        #ax2.legend(plot2, labels2, loc='best')
        
    if(PlotLspp and PlotEpsilons):
        ax3.set_xlabel('Fraction of CrO$_2$ (%)')
        labelsComb2 = [l.get_label() for l in plotComb2]
        labelsComb3 = [l.get_label() for l in plotComb3]
        ax2.legend(plotComb2, labelsComb2, loc='upper left')
        ax3.legend(plotComb3, labelsComb3, loc='upper left')
    elif(PlotLspp and not PlotEpsilons):
        ax2.set_xlabel('Fraction of CrO$_2$ (%)')
        labelsComb2 = [l.get_label() for l in plotComb2]
        ax2.legend(plotComb2, labelsComb2, loc='upper right')
    elif(not PlotLspp and PlotEpsilons): #ax2 does not exist
        ax3.set_xlabel('Fraction of CrO$_2$ (%)')
        labelsComb3 = [l.get_label() for l in plotComb3]
        ax3.legend(plotComb3, labelsComb3, loc='upper right')
    else:
        ax1.set_xlabel(r'Fraction of CrO$_2$ (%)')
    
    labelsComb1 = [l.get_label() for l in plotComb1]
    ax1.legend(plotComb1, labelsComb1, loc='upper right')
    
    
    plt.xlim((0,100))
    plt.tight_layout()
    filename="Dostovalov_"+OxideName1+"-mixedWith-"+SampleName+"-and-"+OxideName2+"-Thickness-"+str(1E9*thickness_max)+"nm"
    plt.savefig(filename+".eps")
    plt.savefig(filename+".png")
    plt.show()
#}}}



## Computes the heating of the sample
# @param wavelength: wavelength of the irradiating photons
# @param epsCr: complex dielectric permittivity of the material
def ThinFilmHeating(wavelength, epsCr): #{{{
    # Estimation of the sample heating
    omega = 2e0*np.pi*c/wavelength
    n_opt = cmath.sqrt(epsCr)
    alpha = 2E0*omega/c * n_opt.imag
    tau = 232e-15
    spot_diam = 15e-6 #defined at 1/e2 
    experimental_energy = 100e-9 #[J]

    print("INPUT: experimental energy ", experimental_energy, "J")
    print("INPUT: spot size at 1/e2 ", spot_diam, "m")
    r_size = 5000
    rmin = 0; rmax = 10*spot_diam; dr = ( rmax - rmin ) / r_size
    r  = np.arange(rmin, rmax, dr)
    intensity_peak = 97623195435.4045E0
    print("ESTIMATION: peak intensity: %5.2e" % intensity_peak, "W/m2")

    intensity_r = intensity_peak * np.exp(-2e0*(r**2/(0.5*spot_diam)**2))

    total_power = np.trapz(intensity_r, r, dx=dr)
    print("total_power: ", total_power, "W")

    total_energy = tau * total_power * np.sqrt(4e0 * np.log(2E0) / np.pi) #normalization due to Gaussian pulse temporal envelope
    print("Total energy: ", total_energy, "J")
    intensity_peak_norm = (total_energy / experimental_energy)**-1
    print("Normalized intensity peak to: ", intensity_peak_norm)
    peak_fluence = intensity_peak * tau * np.sqrt(4e0 * np.log(2E0) / np.pi)

    print("Resulting peak fluence for 1/e2 spot size convention: ", peak_fluence, "J/m2")

    S = alpha * intensity_peak
    C_l = 0.46E3 * 7.2E3 #Bauerle, Edition 4. 
    kappa = 0.97E2 #Bauerle, Edition 4. 
    dt = tau
    Length = 30e-9
    T_surf = 500e0
    T0 = 300E0 #K
    dT_x = T_surf-T0 #temperature gradient over the film depth
    # Equation: C_l * dT / dt = nabla( kappa \nabla ( T ) ) + S
    dT_t = (1E0 / Length**2 * kappa * dT_x + S) / C_l * dt #increase per dt

    print("Increase of lattice temperature:", dT_t, "K")
#}}}

## Compute the SPP modes for mixed oxide ratio using an external set of optical data
# @param fraction: array of oxide fraction
# Size of epsSample and fraction should be of the same dimension
def ScenarioOfCrOxideMixture_ext(epsSample, epsEnvironment=1., epsSubstrate=1., fraction=1, SampleName='CrCompoundOxide', EnvironmentName='Air', SubstrateName='BK7', PlotLspp=False, FilterNegativeLspp=False, PlotEpsilons=False):
    print(Header, "# Info: Considering a mixed fraction of Cr with Cr2O3 and CrO2 (using external data) with several thicknesses.")
    wavelength = 1026e-9
    
    fraction_size = len(fraction)
    numberofroots = 10
    
    #fraction_min = 0.
    #fraction_max = 1.

    thickness_size = 1
    thickness_min  = 100e-9 #35e-9 #50e-9 #50e-9 #70e-9 # 28e-9
    thickness_max  = thickness_min
    
    #PlotLspp = False
    
    if(PlotLspp and PlotEpsilons):
        plotA = 311; plotB=312; plotC=313
    elif((PlotLspp and not PlotEpsilons) or (not PlotLspp and PlotEpsilons)): 
        plotA=211; plotB=212
    else: 
        plotA=111; 
    ## Running 

    #fraction = np.linspace(fraction_min, fraction_max, fraction_size) #fraction of Cr (includes the final value)
    
    print(Header, "# Info: size of the external fraction matrix: ", fraction_size)
    print(Header, "# Info: preparation of the root finder.")
    
    #results = Parallel(n_jobs=num_cores)(delayed(processInput)(i) for i in inputs)

    #t = 100E-9 #thickness of the layer in meters
    #thickness = 10e-9
    print("\n")
    print("# Fraction of Cr oxide: ", fraction)
    summary = np.zeros((0, 7))
    ## Preparation of the thin film modeling for various compositions
    for fraction_index in np.arange(0,fraction_size): #arange excludes the last one, linspace includes it
        # Medium 1: thin film. 
        eps1 = epsSample[fraction_index]        #thin film
        # Medium 2: substrate. 
        eps2 = epsBK7       #epsBK7 #environment | substrate
        # Medium 3: environment
        eps3 = epsEnvironment       #environment | substrate
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
        for branch in np.arange(0, num_branches):
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
    
    print(Header+"** Preparation of the plots as function of oxide ratio")

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
    
    #ReBeta = np.divide(2.*np.pi,period_s)
    #ImBeta = np.divide(0.5,lspp_s)
    
    # Then we could plot them in the right order
    thickness, fractionOxide, epsilonFilm, branch, root_number, period, lspp = SplitSummaryTable(summary)
    
    thickness0, fractionOxide0, epsilonFilm0, branch0, root_number0, period0, lspp0 = SplitSummaryTable(summary_branch0)
    thickness1, fractionOxide1, epsilonFilm1, branch1, root_number1, period1, lspp1 = SplitSummaryTable(summary_branch1)
    thickness2, fractionOxide2, epsilonFilm2, branch2, root_number2, period2, lspp2 = SplitSummaryTable(summary_branch2)
    thickness3, fractionOxide3, epsilonFilm3, branch3, root_number3, period3, lspp3 = SplitSummaryTable(summary_branch3)
    
    summary_export         = np.array([np.real(fractionOxide) , np.real(epsilonFilm) , np.imag(epsilonFilm) , np.real(branch) , np.real(root_number) , np.real(period) ])
    summary_branch0_export = np.array([np.real(fractionOxide0), np.real(epsilonFilm0), np.imag(epsilonFilm0), np.real(branch0), np.real(root_number0), np.real(period0)])
    summary_branch1_export = np.array([np.real(fractionOxide1), np.real(epsilonFilm1), np.imag(epsilonFilm1), np.real(branch1), np.real(root_number1), np.real(period1)])
    summary_branch2_export = np.array([np.real(fractionOxide2), np.real(epsilonFilm2), np.imag(epsilonFilm2), np.real(branch2), np.real(root_number2), np.real(period2)])
    summary_branch3_export = np.array([np.real(fractionOxide3), np.real(epsilonFilm3), np.imag(epsilonFilm3), np.real(branch3), np.real(root_number3), np.real(period3)])
    
    # Adding manual analysis of the obtained results
    # Plotting the lambda/n analysis. Branch number 4. 
    # Plotting the lambda/2n analysis. Branch number 5. 
    PeriodWaveGuideMode = np.divide(wavelength, np.real(np.power(epsilonFilm, 0.5)))
              
    summary_analysis_export = np.array([np.real(fractionOxide) , np.real(epsilonFilm) , np.imag(epsilonFilm), 4*np.ones(np.shape(branch)) , np.zeros(np.shape(root_number)), PeriodWaveGuideMode ])
    
    summary_export_t         = np.transpose(summary_export)
    summary_branch0_export_t = np.transpose(summary_branch0_export)
    summary_branch1_export_t = np.transpose(summary_branch1_export)
    summary_branch2_export_t = np.transpose(summary_branch2_export)
    summary_branch3_export_t = np.transpose(summary_branch3_export)
    summary_analysis_export_t = np.transpose(summary_analysis_export)
    
    # NOTE: I would like to get clean numbers, not the content of summary_branch0
    # Is there a problem with root_number0 for example? 
    print("SPP branches are ready. Exporting to CSV...")
    filename = "Dostovalov-SPPmodes-CrOxideCompound"
    np.savetxt(filename+".csv", summary_export_t)
    np.savetxt(filename+"-Cr-branch0"+".csv", summary_branch0_export_t)
    np.savetxt(filename+"-Cr-branch1"+".csv", summary_branch1_export_t)
    np.savetxt(filename+"-Cr-branch2"+".csv", summary_branch2_export_t)
    np.savetxt(filename+"-Cr-branch3"+".csv", summary_branch3_export_t)
    np.savetxt(filename+"-Cr-branch3"+".csv", summary_branch3_export_t)
    np.savetxt(filename+"-Cr-WaveGuideAnalysis"+".csv", summary_analysis_export_t)
    
    plt.figure()
    ax1 = plt.subplot(plotA)
    #plt.title(r'Film thickness $t=$'+str(round(np.real(thickness[0])*1e9))+' nm')
    #plt.xlabel(r'Fraction of Cr (perc.)')
    plt.ylabel(r'SPP period $\Lambda$ (nm)') 
    #ax1y = ax1.twiny()
    #plot1y1 = ax1.errorbar(np.multiply(CrFraction_Fitted,100), 1e9*ExperimentalData_LSFL, yerr=1e9*ExperimentalData_LSFL_error, fmt='ro', label=r'Period LSFL')
    #plot1y2 = ax1.errorbar(np.multiply(CrFraction_Fitted,100), 1e9*ExperimentalData_HSFL, yerr=1e9*ExperimentalData_HSFL_error, fmt='r^', label=r'Period HSFL')
    #ax1.set_yscale('log')
    plt.ylim((0.,1.1e9*wavelength))
    #ax12 = ax1.twinx()
    #plt.ylabel(r'SPP mean free path $L_{SPP}$ (m)')
    #plt.ylabel(r'Re($\varepsilon$), Im($\varepsilon$)')
    
    plot110, = ax1.plot(np.multiply(100,fractionOxide0), np.multiply(1e9,period0), 'r+', label=r'SPP period $\Lambda$, branch (-,-)')
    plot111, = ax1.plot(np.multiply(100,fractionOxide1), np.multiply(1e9,period1), 'k+', label=r'SPP period $\Lambda$, branch (-,+)')
    plot112, = ax1.plot(np.multiply(100,fractionOxide2), np.multiply(1e9,period2), 'b+', label=r'SPP period $\Lambda$, branch (+,-)')
    plot113, = ax1.plot(np.multiply(100,fractionOxide3), np.multiply(1e9,period3), 'go', label=r'SPP period $\Lambda$, branch (+,+)')
    
    plot12, = ax1.plot(np.multiply(100, fractionOxide), np.multiply(1e9,wavelength*np.ones(np.shape(fractionOxide))), 'k-', linewidth=0.5, label=r'Laser wavelength $\lambda$')
    plot13, = ax1.plot(np.multiply(100, fractionOxide), np.multiply(1E9,PeriodWaveGuideMode), 'k--', linewidth=0.5, label=r'Waveguide modes $\lambda/n$')
    
    #plot14, = ax12.plot(fractionOxide_s, np.real(epsilonFilm_s), 'b+', label=r'Re($\varepsilon$)')
    #plot15, = ax12.plot(fractionOxide_s, np.imag(epsilonFilm_s), 'b^', label=r'Im($\varepsilon$)')
    
    #ax1.yaxis.label.set_color(plot110.get_color()) #colorizes the label
    #ax1.spines["left"].set_edgecolor(plot110.get_color()) #colorizes the axis
    #ax1.tick_params(axis='y', colors=plot110.get_color()) #colorizes the tics and numbers
    plotComb1= [plot110, plot111, plot112, plot113, plot12, plot13]; 
    
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
        ax2.set_xlabel(r'Fraction of Cr oxide (%)')
        #plot21, = ax2.semilogy(fractionOxide_s, lspp_s,   'r^', label=r'SPP decay length $L_{SPP}$')
        plot21, = ax2.semilogy(np.multiply(1E2,fractionOxide0), np.abs(lspp0),   'r^', label=r'SPP decay length $L_{SPP}$, --')
        plot22, = ax2.semilogy(np.multiply(1E2,fractionOxide1), np.abs(lspp1),   'k^', label=r'SPP decay length $L_{SPP}$, -+')
        plot23, = ax2.semilogy(np.multiply(1E2,fractionOxide2), np.abs(lspp2),   'b^', label=r'SPP decay length $L_{SPP}$, +-')
        plot24, = ax2.semilogy(np.multiply(1E2,fractionOxide3), np.abs(lspp3),   'g^', label=r'SPP decay length $L_{SPP}$, ++')
        plot2   = [plot21, plot22, plot23, plot24] #, plot24]
        plotComb2 = plot2
        plt.tight_layout()
        
        #ax22 = ax2.twinx()
    if(PlotEpsilons and not PlotLspp): #plot in plotB window, keep ax3 name. 
        ax3 = plt.subplot(plotB)
    if(PlotEpsilons):
        plt.ylabel(r'Re($\varepsilon$), Im($\varepsilon$)')
        plot31, = ax3.plot(np.multiply(1E2,fractionOxide), np.real(epsilonFilm), 'b+', label=r'Re$(\varepsilon)$ (Cr + oxide)')
        plot32, = ax3.plot(np.multiply(1E2,fractionOxide), np.imag(epsilonFilm), 'b^', label=r'Im$(\varepsilon)$ (Cr + oxide)')
        plot33, = ax3.plot(np.multiply(1E2,fractionOxide), np.multiply(epsBK7, np.ones(np.shape(fractionOxide))), 'k-', label=r'$Re[\varepsilon$(BK7)] ')
        plot3   = [plot31, plot32, plot33]
        plotComb3 = plot3
        plt.tight_layout()
        
        #ax2.yaxis.label.set_color(plot21.get_color()) #colorizes the label
        #ax2.spines["left"].set_edgecolor(plot21.get_color()) #colorizes the axis
        #ax2.tick_params(axis='y', colors=plot21.get_color()) #colorizes the tics and numbers
        
        #ax22.yaxis.label.set_color(plot22.get_color()) #colorizes the label
        #ax22.spines["right"].set_edgecolor(plot22.get_color()) #colorizes the axis
        #ax22.tick_params(axis='y', colors=plot22.get_color()) #colorizes the tics and numbers
        #labels2 = [l.get_label() for l in plot2]
        #ax2.legend(plot2, labels2, loc='best')
        
    if(PlotLspp and PlotEpsilons):
        ax3.set_xlabel('Fraction of Cr oxide (%)')
        labelsComb2 = [l.get_label() for l in plotComb2]
        labelsComb3 = [l.get_label() for l in plotComb3]
        ax2.legend(plotComb2, labelsComb2, bbox_to_anchor=(1.04,1), loc="upper left", mode="expand")
        ax3.legend(plotComb3, labelsComb3, bbox_to_anchor=(1.04,1), loc="upper left", mode="expand")
    elif(PlotLspp and not PlotEpsilons):
        ax2.set_xlabel('Fraction of Cr oxide (%)')
        labelsComb2 = [l.get_label() for l in plotComb2]
        ax2.legend(plotComb2, labelsComb2, bbox_to_anchor=(1.04,1), loc="upper left", mode="expand")
    elif(not PlotLspp and PlotEpsilons): #ax2 does not exist
        ax3.set_xlabel('Fraction of Cr oxide (%)')
        labelsComb3 = [l.get_label() for l in plotComb3]
        ax3.legend(plotComb3, labelsComb3, bbox_to_anchor=(1.04,1), loc="upper left", mode="expand")
    else:
        ax1.set_xlabel('Fraction of Cr oxide (%)')
    
    labelsComb1 = [l.get_label() for l in plotComb1]
    ax1.legend(plotComb1, labelsComb1, bbox_to_anchor=(1.04,1), loc="upper left", mode="expand")
    
    
    #plt.xlim((0,100))
    plt.tight_layout()
    filename="Dostovalov-"+EnvironmentName+"-"+SampleName+"-Thickness-"+str(1E9*thickness_max)+"nm"
    plt.savefig(filename+".eps")
    plt.savefig(filename+".png")
    plt.show()
#}}}

## Trying to repeat optical data provided by Sergei Lisunov. 
def RepeatLisunovMixtureOfOxides(): #{{{
    ## Method 1: using Maxwell-Garnett2 twice. 

    # (CrO2_o = 2, CrO2_e = 1)
    # Means that Fraction(CrO2_o) = 2/3, Fraction(CrO2_e) = 1/3. 
    epsCrO2_2o1e = MaxwellGarnett2(epsCrO2_o, epsCrO2_e, 1./3.) #This does not strongly affect the optical index of CrO2. 
    print("Mixing CrO2 (2o+e): ", epsCrO2_2o1e)
    # Now constructing the mixed Maxwell Garnett2 data for Cr2O3+CrO2. 
    # Careful! epsCr2O3 ratio is not 0.35. 
    # ratio(Cr2O3)/ratio(CrO2) = 0.35, and ratio(Cr2O3)+ratio(CrO2) = 1. 
    # Therefore, ratio(Cr2O3) = 0.25925 and ratio(CrO2) = 0.74074. 
    #/!\ Optics Express from Dostovalov 2018 shows we have MORE CrO2 than Cr2O3. 
    epsCrXOY = MaxwellGarnett2(epsCrO2_2o1e, epsCr2O3, 1.-0.25925) #fraction here refers to medium2
    #epsCrXOY = MaxwellGarnett2(epsCrO2_2o1e, epsCr2O3, 1.-0.35) #fraction here refers to medium2
    print("Mixing CrO2 (mixed) + Cr2O3 gives", epsCrXOY)
    fractionOxide    = np.linspace(0, 1, Fraction_size) #0: 100% Cr, 1: 100% oxide
    CrMixedWithCrXOY = MaxwellGarnett2(epsCr, epsCrXOY, fractionOxide)

    #CrCrXOY_Lisunov  = np.loadtxt("Dostovalov-Cr/Cr-Cr2O3-CrO2/Sergei_Lisunov/OptProperties_Cr_with_oxides.csv", skiprows=2)
    CrCrXOY_Lisunov  = np.loadtxt("Dostovalov-Cr/Cr-Cr2O3-CrO2/Sergei_Lisunov/OptProperties_Cr_with_oxides_corrected.csv", skiprows=2)
    limiter = 2 #limit the number of cells to get, then we can update the plot without recomputing the whole thing.
    CrCrXOY_fraction   = CrCrXOY_Lisunov[:,0]
    epsR_CrCrXOY_L     = CrCrXOY_Lisunov[:,1]
    epsC_CrCrXOY_L     = CrCrXOY_Lisunov[:,2]
    eps_CrCrXOY_L = np.add(epsR_CrCrXOY_L, np.multiply(1.j, epsC_CrCrXOY_L))

    plt.figure()
    plt.plot(100*fractionOxide, np.ones(np.shape(fractionOxide))*epsCrO2_2o1e.real, 'k-', label=r'Cr$_2$O$_3$ (2o+e), Re$(\varepsilon)$')
    plt.plot(100*fractionOxide, np.ones(np.shape(fractionOxide))*epsCrO2_2o1e.imag, 'k--', label=r'Cr$_2$O$_3$ (2o+e), Im$(\varepsilon)$')
    plt.plot(100*fractionOxide, np.ones(np.shape(fractionOxide))*epsCrXOY.real, 'g-', label=r'Cr$_2$O$_3$+CrO$_2$, Re$(\varepsilon)$')
    plt.plot(100*fractionOxide, np.ones(np.shape(fractionOxide))*epsCrXOY.imag, 'g--', label=r'Cr$_2$O$_3$+CrO$_2$, Im$(\varepsilon)$')
    plt.plot(100*fractionOxide, CrMixedWithCrXOY.real, 'r-', label=r'Cr+Cr$_x$O$_y$, Re$(\varepsilon)$')
    plt.plot(100*fractionOxide, CrMixedWithCrXOY.imag, 'r--', label=r'Cr+Cr$_x$O$_y$, Im$(\varepsilon)$')
    plt.plot(100*CrCrXOY_fraction, epsR_CrCrXOY_L, 'b+', label='Re(eps), Sergei data using inverted o/e modes')
    plt.plot(100*CrCrXOY_fraction, epsC_CrCrXOY_L, 'b.', label='Im(eps), Sergei data using inverted o/e modes')
    plt.xlabel('Fraction of oxide (%)')
    plt.ylabel(r'Re$(\varepsilon)$, Im$(\varepsilon)$')
    plt.legend(loc='best')
    plt.show()

    # This does not lead to repeat data from Sergei. 
    # Email mentioned they have used a mixture of (o) mode and (e) mode. 
    # Impossible to find his results as well. I conclude that it originates from lack of precision in sampling CrO2(o) and CrO2(e) in the paper of Chase. Error should come from Sergey as I did carefully dit it with Digitizer. He might have done it with a ruler, maybe?  

    ### Method 2: using MaxwellGarnett3
    #ratio_Cr2O3_CrO2 = 0.35
    #FractionCrO2  = np.linspace(0.,1.,Fraction_size)
    #FractionCr2O3 = ratio_Cr2O3_CrO2 * FractionCrO2
    #FractionCr    = 1-ratio_Cr2O3_CrO2 * FractionCrO2 - FractionCrO2

    #CrMixedWithCrXOY_3 = MaxwellGarnett3(epsCr, epsCr2O3, epsCrO2, FractionCr, FractionCr2O3, FractionCr)

    #print(CrMixedWithCrXOY_3)

    #plt.figure()
    #plt.plot(100*FractionCrO2, CrMixedWithCrXOY_3.real, label='Method 2, Re(eps)')
    #plt.plot(100*FractionCrO2, CrMixedWithCrXOY_3.imag, label='Method 2, Im(eps)')
    #plt.xlabel('Fraction of oxide (%)')
    #plt.ylabel(r'Re$(\varepsilon)$, Im$(\varepsilon)$')
    #plt.legend(loc='best')
    #plt.show()
    ##ThinFilmHeating(wavelength, epsCr)

    ## Result is extremely different using this approach than using MaxwellGarnett2(Cr, MaxwellGarnett2(Cr2O3, MaxwellGarnett2(CrO2_o, CrO2_e))). 
    ## I believe the right approach would be to use MaxwellGarnett4(Cr, Cr2O3, CrO2_o, CrO2_e). 
#}}}

def PreparePublicationFigure_OxideFraction():
    ## Takes ~ 30 min run
    ## Preparing SPP period using an external file
    #CrCrXOY_Lisunov  = np.loadtxt("Dostovalov-Cr/Cr-Cr2O3-CrO2/Sergei_Lisunov/OptProperties_Cr_with_oxides.csv", skiprows=2)
    CrCrXOY_Lisunov  = np.loadtxt("Dostovalov-Cr/Cr-Cr2O3-CrO2/Sergei_Lisunov/OptProperties_Cr_with_oxides_corrected.csv", skiprows=0)
    limiter = 2 #limit the number of cells to get, then we can update the plot without recomputing the whole thing.
    CrCrXOY_fraction   = CrCrXOY_Lisunov[:,0]
    epsR_CrCrXOY_L     = CrCrXOY_Lisunov[:,1]
    epsC_CrCrXOY_L     = CrCrXOY_Lisunov[:,2]
    eps_CrCrXOY_L = np.add(epsR_CrCrXOY_L, np.multiply(1.j, epsC_CrCrXOY_L))
    NumberOfSuperImposedPlots=1
    Every = 20*NumberOfSuperImposedPlots
    Shift = int(0*Every/NumberOfSuperImposedPlots) #enable to plot shifted plots to avoid superimposition
    ScenarioOfCrOxideMixture_ext(eps_CrCrXOY_L[Shift::Every], epsAir, epsBK7, CrCrXOY_fraction[Shift::Every], 'Cr_compounds_oxide', 'Air', 'BK7')

# =====================
#ScenarioOfOxidePrecipitation()
# ** Info: computing 3-layer reflectivity..."
#R = BiLayerReflectivity(epsAir, epsCrCr2O3_list, epsBK7, t_list) #dimension is good for a HeatMap picture

## Validation cases in Python. 
Fraction_size = 30 #number of samples
thickness_size = 60 #Fraction_size
#Burke_SymmetricModes(thickness_size)

#ScenarioOfCrOxideMixture(epsCr, epsCr2O3, epsBK7, epsAir, Fraction_size, 'Cr', 'Cr2O3')
#ScenarioOfCrOxideMixture(epsCr, epsCrO2, epsBK7, epsAir, Fraction_size,  'Cr', 'CrO2')
#ScenarioOfCrOxideMixture3(epsCr, epsCr2O3, epsCrO2, epsBK7, Fraction_size, 'Cr', 'Cr2O3', 'CrO2')

RepeatLisunovMixtureOfOxides()

#PreparePublicationFigure_OxideFraction()
#ScenarioOfCrOxideMixture_ext(eps_CrCrXOY_L[Shift::Every], epsAir, epsBK7, CrCrXOY_fraction[Shift::Every], 'Cr_compounds_oxide', 'Air', 'BK7')
Te_max = 1E10
#PeriodsAsFunctionOfTemperature(epsTi, epsBK7, epsAir, Fraction_size, 'Cr', 'BK7', True, False, True, Te_max)
