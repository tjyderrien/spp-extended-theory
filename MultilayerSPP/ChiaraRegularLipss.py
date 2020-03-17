#!/usr/bin/env python2
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

## @module plotTemperature.py presents properties of SPP as function of temperature of a material. 
# Of course, temperature-dependent optical data of the material are necessary. 

# IMPORT LIBRARIES
from libSPP import *
from libMaterials import *
from scipy.constants import c, epsilon_0, mu_0, pi, e, m_e, h, hbar
from libMultilayerSPP import *
from plotSipe import *

## Returns the dielectric permittivity of liquid Si at 1687 K: 
# Gellison Applied Physics Letters, 51, 352 (1987)
def DielectricPermittivity_MoltenSilicon(wavelength):
  beta = 8.2e0*e/hbar
  omegac = 34.1e0*e/hbar
  omega = 2e0*pi*c/wavelength
  dielectricR = 1e0 - (omegac * omegac + beta * beta)/(omega * omega + beta * beta);
  dielectricI = beta / (2e0 * omega) * (omegac * omegac + beta * beta) / (omega * omega + beta * beta)
  eps = dielectricR + dielectricI * 1j
  return eps

def SPPonFullyMoltenSi(wavelength): 
    print "SIMPLE SPP MODEL with MOLTEN Si data"
    eps_MoltenSi = DielectricPermittivity_MoltenSilicon(wavelength)
    #print eps_MoltenSi
    print "Dielectric permittivity of molten Si at "+str(1E9*wavelength)+" nm = "+str(eps_MoltenSi)
    beta = betaSPP(wavelength, 1., eps_MoltenSi) 
    print "SPP period (nm): "+str(period(beta)*1E9)
    print "SPP mean-free-path (um): "+str(DecayLengthSPP(beta)*1E6)


def MoltenSiLayerOnSi(wavelength): #{{{ 
    print "MultilayerSPP: computing the possible SPP periods for various thicknesses of molten Si."

    # Medium 1: thin film. 
    eps1 = eps_MoltenSi        #thin film
    # Medium 2: substrate. 
    eps2 = epsAir       #epsBK7 #environment | substrate
    # Medium 3: environment
    eps3 = eps_Si_1030nm       #environment | substrate
    # Note: Inverting eps2 and eps3 should have no effect on the possible modes, but only on field amplification. 

    #t = 100E-9 #thickness of the layer in meters
    #thickness = 100e-9
    t_list = np.power(10., np.linspace(np.log10(1e-9), np.log10(50e-9), NumberOfPoints))
    #meshes the initial guess area, all numbers are from the space of betas
    x_max = 1E8 # 6e7       #1E10
    y_max = 1E9
                #
    x_min = -x_max # -6e7       #-1E10
    y_min = - y_max
            
    x_steps =  50    #70
    y_steps =  50    #70

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


    summary = np.zeros((0, 7))

    for thickness in t_list:
        roots = findroots(eps1, eps2, eps3,
                wavelength, thickness,
                x_min, x_max,    
                y_min, y_max,    
                x_steps, y_steps)
                
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


    ## Shaping the data to plot them with GNUplot
    #roots_shape = np.shape(roots)
    ##print(roots_shape)
    #num_thickness = np.shape(thickness) #NOTE: is this used? 
    #num_branches  = roots_shape[0]
    #num_roots     = roots_shape[1]
    #num_property  = roots_shape[2]

    #for branch in np.arange(0,num_branches):
        #for root_number in np.arange(0,num_roots): 
            #print(thickness, roots[branch][root_number][0], roots[branch][root_number][1])
        ##print("\n")
        
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
    plt.ylabel(r'SPP period $\Lambda$ (nm)') 
    #ax1.set_yscale('log')
    #plt.ylim((0.,1.1e9*wavelength))
    plot110, = ax1.plot(np.multiply(1e9,thickness0), np.multiply(1e9,period0), 'r+', label=r'SPP period $\Lambda$, branch (-,-)')
    plot111, = ax1.plot(np.multiply(1e9,thickness1), np.multiply(1e9,period1), 'k+', label=r'SPP period $\Lambda$, branch (-,+)')
    plot112, = ax1.plot(np.multiply(1e9,thickness2), np.multiply(1e9,period2), 'b+', label=r'SPP period $\Lambda$, branch (+,-)')
    plot113, = ax1.plot(np.multiply(1e9,thickness3), np.multiply(1e9,period3), 'g+', label=r'SPP period $\Lambda$, branch ( +,+)')

    plot12,  = ax1.plot(np.multiply(1e9, thickness), np.multiply(1e9,wavelength*np.ones(np.shape(fractionOxide))), 'k-', linewidth=0.5, label=r'Laser wavelength $\lambda$')

    #ax12 = ax1.twinx()
    #plot14, = ax12.plot(fractionOxide_s, np.real(epsilonFilm_s), 'b+', label=r'Re($\varepsilon$)')
    #plot15, = ax12.plot(fractionOxide_s, np.imag(epsilonFilm_s), 'b^', label=r'Im($\varepsilon$)')
    #plot14,  = ax1.plot(np.multiply(1e9,thickness0), np.multiply(1e0,PeriodToBetaNorm(period0)), 'r^', label=r'$\beta/k_0$, (-,-)')
    #plot15,  = ax1.plot(np.multiply(1e9,thickness1), np.multiply(1e0,PeriodToBetaNorm(period1)), 'k^', label=r'$\beta/k_0$, (-,+)')
    #plot16,  = ax1.plot(np.multiply(1e9,thickness2), np.multiply(1e0,PeriodToBetaNorm(period2)), 'b^', label=r'$\beta/k_0$, (+,-)')
    #plot17,  = ax1.plot(np.multiply(1e9,thickness3), np.multiply(1e0,PeriodToBetaNorm(period3)), 'g^', label=r'$\beta/k_0$, (+,+)')

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
    plt.savefig("MoltenSionBareSi-1030nm.png")
#}}}
    
wavelength = 1030e-9
NumberOfPoints=50

epsAir        = 1.
eps_MoltenSi  = DielectricPermittivity_MoltenSilicon(wavelength)
eps_Si_1030nm = 12.80259 + 0.0109j

# SPPonFullyMoltenSi(wavelength) # not so interesting
# MoltenSiLayerOnSi(wavelength)   # there are smaller modes, but nothing thickness-stable around 700 nm. 

## DO NOT FORGET: try Sipe model on molten Si
plotSipe1D_sectionX(1030e-9, eps_MoltenSi)
