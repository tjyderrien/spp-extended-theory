#!/usr/bin/env python3
#-*- coding: utf-8 -*-

# Copyright (C) 2018-2020 T.J.-Y. Derrien
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

## @package MultilayerKovaricek
# Preparation of results for Prof. Bulgakova and Sasha Dostovalov. 
# Explores the SPP theory at a thin film located 
# between two semi-infinite media. The formal model is presented in 
# T.J.-Y. Derrien et al, J. Appl. Phys. 116, 074902 (2014) and references 
# therein. 

from libMultilayerSPP import *
import sys
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
NumberOfPoints=5

#data
wavelength = 1030e-9 #1064e-9

#epsCr2O3   = 3.82738158014083 + 0.0483802637311967j  #Al-Kuhaili, M. & Durrani, S. Optical properties of chromium oxide thin films deposited by electron-beam evaporation Optical Materials, 2007, 29, 709-713
#epsCr2O3    = 4.9713+0.1784j  #1 um [JDT Kruschwitz et al, Appl. Opt. 1997]
#epsCr      = -0.6721223+24.8657476j
#epsCrO2_e    = 1.3587463082734004+9.00595525243578j #Chase, L. L. Optical properties of Cr O 2 and Mo O 2 from 0.1 to 6 eV Physical Review B, 1974, 10, 2226-2231 (E || C mode. C: axis for extraordinary mode). 
#epsCrO2_e = epsCrO2 #(E || c)
#epsCrO2_o  = 0.5784455474251002+6.477144040000001j #Chase, L. L. Optical properties of Cr O 2 and Mo O 2 from 0.1 to 6 eV Physical Review B, 1974, 10, 2226-2231 (E perpendicular to C).

#epsTi = -4.289599704142011+27.217715606508875j

## TODO: replace this by a function taking data in MaterialOpticalData.csv !
if(wavelength == 1064e-9): 
    epsSiO2=1.4496**2; 
    epsSi=12.6893281+0.0067935529j
    epsSiLiquid=-16.92943996+62.91836213j
elif(wavelength==1030e-9):
    epsSiO2 = 2.1026565205
    epsSi   = 12.80259+0.0109j
    epsSiLiquid=-16.90553111+60.8265925j

epsAir      = 1.+0.j          #air

neSiLog     = np.linspace(25,np.log10(4*5.E28),NumberOfPoints)
neSi        = np.power(10., neSiLog)

numberofroots_per_branch=10

t_list = [300e-9] #500e-9]

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
        fractionOxide = summary[:, 1] #can play the role of excited density
        branch        = summary[:, 2]
        root_number   = summary[:, 3]
        period        = summary[:, 4]
        lspp          = summary[:, 5]
        epsilonFilm   = summary[:, 6]
        
        lists = sorted(zip(*[thickness, fractionOxide, epsilonFilm, branch, root_number, period, lspp]))
        thickness_s, fractionOxide_s, epsilonFilm_s, branch_s, root_number_s, period_s, lspp_s  = list(zip(*lists))
    return thickness_s, fractionOxide_s, epsilonFilm_s, branch_s, root_number_s, period_s, lspp_s
#}}}


eps1 = epsSiO2
eps2 = epsAir

#print(Header+"** Debug: neSi")
#print(neSi)

def ExcitationOfSiLayer():
    summary = np.zeros((0, 7))
    for excitation_index in np.arange(0,np.shape(neSi)[0]-1):
    #for excitation_index in np.arange(0,np.shape(epsTiO2)[0]-1):
        # Medium 1: thin film. 
        #print("Order: ", excitation_index)
        eps3 = Drude(wavelength, neSi[excitation_index], epsSi, 1.1e-15**-1, 0.18)
        #eps3 = epsSiLiquid
        #print(excitation_index)
        #print(neSi[excitation_index])
        #eps1=epsSiO2 #epsTiO2[excitation_index]
        #thin film

        for thickness in t_list:
            roots = findroots(eps1, eps2, eps3,
                    wavelength, thickness,
                    x_min, x_max,    
                    y_min, y_max,    
                    x_steps, y_steps, numberofroots_per_branch)

            # Shaping the data to plot them with GNUplot
            roots_shape = np.shape(roots)
            #print(roots_shape)
            num_thickness = np.shape(t_list) #NOTE: is this used? 
            
            num_branches = len(roots) 
            print("Number of branches: "+str(num_branches)) #number of SPP branches for this sample. 
            
            for branch in np.arange(0,num_branches):
                roots_in_branch = roots[branch]
                print("Roots in branch #"+str(branch))
                for order in np.arange(0,len(roots_in_branch)): 
                    roots_in_branch_order = roots_in_branch[order]
                    print("Order #"+str(order)+": Period="+str(roots_in_branch_order[0])+" Lspp="+str(roots_in_branch_order[1]))
                    #TODO: #signs are necessary here... 
                    decaydepth1=0 
                    decaydepth2=0
                    decaydepth3=0
                    ToBeAdded = [thickness, neSi[excitation_index], branch, order, roots_in_branch_order[0], roots_in_branch_order[1], eps1] #, decaydepth1, decaydepth2, decaydepth3]
                    if(roots_in_branch[order][0] > 1E-15 and abs(roots_in_branch[order][1]) > 1E-10): 
                        # We remove modes were Lspp < 0.1 nm or period < 0. 
                        summary = np.vstack((summary, ToBeAdded ))
                    #print(neSi[excitation_index], thickness, roots[branch][order][0], roots[branch][order][1], eps1.real, eps1.imag)
                    #print(eps1.real, thickness, roots[branch][order][0], roots[branch][order][1], eps1.real, eps1.imag, eps1.real*eps2.real+eps1.imag*eps2.imag, eps1.real*eps3.real+eps1.imag*eps3.imag)
                    #print(neSi[excitation_index], thickness, roots[branch][order][0],
                #print("\n")

    ## Exporting results per branch
    summary_branch0 = np.array(summary[summary[:,2]==0,:])
    summary_branch1 = np.array(summary[summary[:,2]==1,:])
    summary_branch2 = np.array(summary[summary[:,2]==2,:])
    summary_branch3 = np.array(summary[summary[:,2]==3,:])

    #ReBeta = np.divide(2.*np.pi,period_s)
    #ImBeta = np.divide(0.5,lspp_s)

    #print(summary)

    # Then we could plot them in the right order
    thickness, n_exc, epsilonFilm, branch, root_number, period, lspp = SplitSummaryTable(summary)

    thickness0, n_exc_0, epsilonFilm0, branch0, root_number0, period0, lspp0 = SplitSummaryTable(summary_branch0)
    thickness1, n_exc_1, epsilonFilm1, branch1, root_number1, period1, lspp1 = SplitSummaryTable(summary_branch1)
    thickness2, n_exc_2, epsilonFilm2, branch2, root_number2, period2, lspp2 = SplitSummaryTable(summary_branch2)
    thickness3, n_exc_3, epsilonFilm3, branch3, root_number3, period3, lspp3 = SplitSummaryTable(summary_branch3)

    print(n_exc)
    #print(n_exc_0)

    summary_export         = np.array([n_exc,   np.real(epsilonFilm ), np.imag(epsilonFilm) , np.real(branch) , np.real(root_number) , np.real(period) ])
    summary_branch0_export = np.array([n_exc_0, np.real(epsilonFilm0), np.imag(epsilonFilm0), np.real(branch0), np.real(root_number0), np.real(period0)])
    summary_branch1_export = np.array([n_exc_1, np.real(epsilonFilm1), np.imag(epsilonFilm1), np.real(branch1), np.real(root_number1), np.real(period1)])
    summary_branch2_export = np.array([n_exc_2, np.real(epsilonFilm2), np.imag(epsilonFilm2), np.real(branch2), np.real(root_number2), np.real(period2)])
    summary_branch3_export = np.array([n_exc_3, np.real(epsilonFilm3), np.imag(epsilonFilm3), np.real(branch3), np.real(root_number3), np.real(period3)])

    summary_export_t         = np.transpose(summary_export)
    summary_branch0_export_t = np.transpose(summary_branch0_export)
    summary_branch1_export_t = np.transpose(summary_branch1_export)
    summary_branch2_export_t = np.transpose(summary_branch2_export)
    summary_branch3_export_t = np.transpose(summary_branch3_export)

    # NOTE: I would like to get clean numbers, not the content of summary_branch0
    # Is there a problem with root_number0 for example? 
    print("SPP branches are ready. Exporting to CSV...")
    filename = "SiO2Si_exc-SPPmodes-"+str(int(1E9*wavelength))
    np.savetxt(filename+".csv", summary_export_t)
    np.savetxt(filename+"-branch0"+".csv", summary_branch0_export_t)
    np.savetxt(filename+"-branch1"+".csv", summary_branch1_export_t)
    np.savetxt(filename+"-branch2"+".csv", summary_branch2_export_t)
    np.savetxt(filename+"-branch3"+".csv", summary_branch3_export_t)

    # Best would be to output confinement k1 directly from the solver. 

    plt.figure()
    ax1 = plt.subplot(111)
    plt.xlabel(r'Excited electron density in Si (m$^{-3}$)') 
    plt.ylabel(r'SPP period $\Lambda$ (nm)') 
    #ax1y = ax1.twiny()
    #plot1y1 = ax1.errorbar(np.multiply(CrFraction_Fitted,100), 1e9*ExperimentalData_LSFL, yerr=1e9*ExperimentalData_LSFL_error, fmt='ro', label=r'Period LSFL')
    #plot1y2 = ax1.errorbar(np.multiply(CrFraction_Fitted,100), 1e9*ExperimentalData_HSFL, yerr=1e9*ExperimentalData_HSFL_error, fmt='r^', label=r'Period HSFL')
    #ax1.set_yscale('log')
    plt.ylim((0.,2e9*wavelength))
    #ax12 = ax1.twinx()
    #plt.ylabel(r'SPP mean free path $L_{SPP}$ (m)')
    #plt.ylabel(r'Re($\varepsilon$), Im($\varepsilon$)')

    plot110, = ax1.semilogx(n_exc_0, np.multiply(1e9,period0), 'r+', label=r'SPP period $\Lambda$, branch (-,-)')
    plot111, = ax1.semilogx(n_exc_1, np.multiply(1e9,period1), 'k+', label=r'SPP period $\Lambda$, branch (-,+)')
    plot112, = ax1.semilogx(n_exc_2, np.multiply(1e9,period2), 'b+', label=r'SPP period $\Lambda$, branch (+,-)')
    plot113, = ax1.semilogx(n_exc_3, np.multiply(1e9,period3), 'go', label=r'SPP period $\Lambda$, branch (+,+)')

    RefractiveIndex = np.power(epsilonFilm,0.5)
    plot12, = ax1.semilogx(n_exc,    np.multiply(1e9,wavelength*np.ones(np.shape(n_exc))), 'k-', linewidth=0.5, label=r'Laser wavelength $\lambda$')
    plot13, = ax1.semilogx(n_exc,    np.multiply(1e9,np.divide(wavelength, RefractiveIndex.real)), 'k--', linewidth=0.5, label=r'Waveguiding modes $\lambda/n$')

    #plot14, = ax12.plot(fractionOxide_s, np.real(epsilonFilm_s), 'b+', label=r'Re($\varepsilon$)')
    #plot15, = ax12.plot(fractionOxide_s, np.imag(epsilonFilm_s), 'b^', label=r'Im($\varepsilon$)')

    #ax1.yaxis.label.set_color(plot110.get_color()) #colorizes the label
    #ax1.spines["left"].set_edgecolor(plot110.get_color()) #colorizes the axis
    #ax1.tick_params(axis='y', colors=plot110.get_color()) #colorizes the tics and numbers
    plotComb1= [plot110, plot111, plot112, plot113, plot12, plot13]; 
    labelsComb1 = [l.get_label() for l in plotComb1]
    ax1.legend(plotComb1, labelsComb1, loc="best")
        
    thickness_max = np.max(t_list)

    #ax12.yaxis.label.set_color(plot14.get_color()) #colorizes the label
    #ax12.spines["right"].set_edgecolor(plot14.get_color()) #colorizes the axis
    #ax12.tick_params(axis='y', colors=plot14.get_color()) #colorizes the tics and numbers

    plt.tight_layout()
    filename="SiO2Si_exc-wavelength-"+str(int(1E9*wavelength))+"-Thickness-"+str(1E9*thickness_max)+"nm"
    plt.savefig(filename+".eps")
    plt.savefig(filename+".png")
    plt.show()


def ExcitationOfSiLayer(wavelength, thickness_size): 
    #wavelength=800e-9 #1030e-9
    numberofroots = 10
    # Medium 1: thin film.
    eps1 = epsSiO2  
    # Medium 2: substrate. 
    eps2 = epsSiLiquid #4. #3.999999+0.004j
    # Medium 3: environment
    eps3 = epsAir #.5**2 # eps2 #eps2: symmetric modes       #environment | substrate
    # Note: Inverting eps2 and eps3 should have no effect on the possible modes, but only on field amplification. 
    
    #thickness_size = 20
    thickness_min  = 1e-9
    thickness_max  = 1000e-9
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
                print("beta/k0="+str(PeriodToBetaNorm(roots_in_branch_order[0], wavelength)))
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
    
    plt.figure()
    ax1 = plt.subplot(111)
    plt.xlabel(r'Thickness (nm)')
    plt.ylabel(r'SPP period $\Lambda$ (nm)') 
    #ax1.set_yscale('log')
    plt.ylim((0.,2*1e9*wavelength))
    plot110, = ax1.plot(np.multiply(1e9,thickness0), np.multiply(1e9,period0), 'r+', label=r'SPP period $\Lambda$, branch (-,-)')
    plot111, = ax1.plot(np.multiply(1e9,thickness1), np.multiply(1e9,period1), 'k+', label=r'SPP period $\Lambda$, branch (-,+)')
    plot112, = ax1.plot(np.multiply(1e9,thickness2), np.multiply(1e9,period2), 'b+', label=r'SPP period $\Lambda$, branch (+,-)')
    plot113, = ax1.plot(np.multiply(1e9,thickness3), np.multiply(1e9,period3), 'g+', label=r'SPP period $\Lambda$, branch ( +,+)')
    
    plot12,  = ax1.plot(np.multiply(1e9, thickness), np.multiply(1e9,wavelength*np.ones(np.shape(fractionOxide))), 'k-', linewidth=0.5, label=r'Laser wavelength $\lambda$')
    
    #ax12 = ax1.twinx()
    #plot14, = ax12.plot(fractionOxide_s, np.real(epsilonFilm_s), 'b+', label=r'Re($\varepsilon$)')
    #plot15, = ax12.plot(fractionOxide_s, np.imag(epsilonFilm_s), 'b^', label=r'Im($\varepsilon$)')
    #plot14,  = ax1.plot(np.multiply(1e9,thickness0), np.multiply(1e9,np.abs(period0)), 'r^', label=r'$\Lambda_{\mathrm{SPP}}$, (-,-)')
    #plot15,  = ax1.plot(np.multiply(1e9,thickness1), np.multiply(1e9,np.abs(period1)), 'k^', label=r'$\Lambda_{\mathrm{SPP}}$, (-,+)')
    #plot16,  = ax1.plot(np.multiply(1e9,thickness2), np.multiply(1e9,np.abs(period2)), 'b^', label=r'$\Lambda_{\mathrm{SPP}}$, (+,-)')
    #plot17,  = ax1.plot(np.multiply(1e9,thickness3), np.multiply(1e9,np.abs(period3)), 'g^', label=r'$\Lambda_{\mathrm{SPP}}$, (+,+)')
    
    #plt.ylabel(r'SPP mean free path $L_{SPP}$ (m)')
    #plt.ylabel(r'Re($\varepsilon$), Im($\varepsilon$)')
    #plt.ylabel(r'$\beta/k_0$')
    #ax1.yaxis.label.set_color(plot110.get_color()) #colorizes the label
    #ax1.spines["left"].set_edgecolor(plot110.get_color()) #colorizes the axis
    #ax1.tick_params(axis='y', colors=plot110.get_color()) #colorizes the tics and numbers
    plotComb1 = []
    plotComb1+= [plot110, plot111, plot112, plot113, plot12]; 
    #plotComb1+=[plot14, plot15, plot16, plot17]
    
    #ax12.yaxis.label.set_color(plot14.get_color()) #colorizes the label
    #ax12.spines["right"].set_edgecolor(plot14.get_color()) #colorizes the axis
    #ax12.tick_params(axis='y', colors=plot14.get_color()) #colorizes the tics and numbers
    plt.tight_layout()
    labelsComb1 = [l.get_label() for l in plotComb1]
    ax1.legend(plotComb1, labelsComb1, loc="upper right")
    plt.show()


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

# =====================
#ScenarioOfOxidePrecipitation()
# ** Info: computing 3-layer reflectivity..."
#R = BiLayerReflectivity(epsAir, epsCrCr2O3_list, epsBK7, t_list) #dimension is good for a HeatMap picture

## Validation cases in Python. 
Fraction_size = 30 #number of samples
thickness_size = 10 #Fraction_size

ExcitationOfSiLayer(wavelength, thickness_size) #good, but it may occur that Si does not get so excited during irradiation of Si in ps regime. 

#EffectofSiMoltenThickness()

sys.exit()
#ScenarioOfCrOxideMixture(epsCr, epsCr2O3, epsBK7, epsAir, Fraction_size, 'Cr', 'Cr2O3')
#ScenarioOfCrOxideMixture(epsCr, epsCrO2, epsBK7, epsAir, Fraction_size,  'Cr', 'CrO2')
#ScenarioOfSimultaneousMixingMG3(epsCr, epsCr2O3, epsCrO2, epsBK7, Fraction_size, 'Cr', 'Cr2O3', 'CrO2')
