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

## @package MultilayerMirzaLevy
# Preparation of results for Prof. Bulgakova, Inam Mirza and Yoann Levy. 
# Explores the SPP theory at a thin film located 
# between two semi-infinite media. The formal model is presented in 
# T.J.-Y. Derrien et al, J. Appl. Phys. 116, 074902 (2014) and references 
# therein. 

from libMultilayerSPP import *
import matplotlib.pyplot as plt
from libMaterials import *

#precision
NumberOfPoints           = 200
numberofroots_per_branch = 10
#data
wavelength = 800e-9 #355e-9 #1030E-9 #1026

#neH20 = 1E28 #the most violent change
epsTibare   = -6.206969+25.2j #800 nm
epsTiO2bare = 7.7841+0.j      #800 nm
epsAir      = 1.+0.j          #air

neTiO2      = np.power(10., np.linspace(np.log10(1E24), np.log10(1E29), NumberOfPoints))
#epsTiO2     = np.linspace(-2,1,NumberOfPoints)
#print("Shape of Ne array: ", np.shape(neTiO2)[0])

#t = 100E-9 #thickness of the layer in meters
t_list = [151e-9] #thickness = 10e-9
#t_list = np.power(10., np.linspace(np.log10(1e-9), np.log10(300e-9), NumberOfPoints))
#meshes the initial guess area, all numbers are from the space of betas
x_min = -6e7 #-1E10
x_max = 6e7 #1E10

y_min = -1E9
y_max = 1E9

x_steps = 30 #40
y_steps = 30 #40

# Medium 2: substrate. 
eps2 = epsTibare      #epsBK7 #environment | substrate
# Medium 3: environment
eps3 = epsAir       #environment | substrate
# Note: Inverting eps2 and eps3 should have no effect on the possible modes, but only on field amplification. 

## Reorganizes the order of fields and output only necessary information
def SplitSummaryTable(summary): #{{{
    if(len(summary[:,0])==0): #if table is empty, avoids the crash
        thickness_s = 0e0; n_exc_s = -1e0; branch_s=-1; root_number_s=-1; period_s=-1; lspp_s=-1; epsilonFilm_s=1E99
    else:
        thickness     = summary[:, 0] 
        n_exc         = summary[:, 1]
        branch        = summary[:, 2]
        root_number   = summary[:, 3]
        period        = summary[:, 4]
        lspp          = summary[:, 5]
        epsilonFilm   = summary[:, 6]
        
        lists = sorted(zip(*[thickness, n_exc, epsilonFilm, branch, root_number, period, lspp]))
        thickness_s, n_exc_s, epsilonFilm_s, branch_s, root_number_s, period_s, lspp_s  = list(zip(*lists))
    return thickness_s, n_exc_s, epsilonFilm_s, branch_s, root_number_s, period_s, lspp_s
#}}}

summary = np.zeros((0, 7))
for excitation_index in np.arange(0,np.shape(neTiO2)[0]-1):
#for excitation_index in np.arange(0,np.shape(epsTiO2)[0]-1):
    # Medium 1: thin film. 
    #print("Order: ", excitation_index)
    eps1 = Drude(wavelength, neTiO2[excitation_index], epsTiO2bare, 1.1e-15**-1, 0.18)
    #eps1=epsTiO2bare #epsTiO2[excitation_index]
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
                ToBeAdded = [thickness, neTiO2[excitation_index], branch, order, roots_in_branch_order[0], roots_in_branch_order[1], eps1]
                if(roots_in_branch[order][0] > 1E-15 and abs(roots_in_branch[order][1]) > 1E-10): 
                    # We remove modes were Lspp < 0.1 nm or period < 0. 
                    summary = np.vstack((summary, ToBeAdded ))
                #print(neTiO2[excitation_index], thickness, roots[branch][order][0], roots[branch][order][1], eps1.real, eps1.imag)
                #print(eps1.real, thickness, roots[branch][order][0], roots[branch][order][1], eps1.real, eps1.imag, eps1.real*eps2.real+eps1.imag*eps2.imag, eps1.real*eps3.real+eps1.imag*eps3.imag)
                #print(neTiO2[excitation_index], thickness, roots[branch][order][0],
            #print("\n")

## Exporting results per branch
summary_branch0 = np.array(summary[summary[:,2]==0,:])
summary_branch1 = np.array(summary[summary[:,2]==1,:])
summary_branch2 = np.array(summary[summary[:,2]==2,:])
summary_branch3 = np.array(summary[summary[:,2]==3,:])

#ReBeta = np.divide(2.*np.pi,period_s)
#ImBeta = np.divide(0.5,lspp_s)

# Then we could plot them in the right order
thickness, n_exc, epsilonFilm, branch, root_number, period, lspp = SplitSummaryTable(summary)

thickness0, n_exc_0, epsilonFilm0, branch0, root_number0, period0, lspp0 = SplitSummaryTable(summary_branch0)
thickness1, n_exc_1, epsilonFilm1, branch1, root_number1, period1, lspp1 = SplitSummaryTable(summary_branch1)
thickness2, n_exc_2, epsilonFilm2, branch2, root_number2, period2, lspp2 = SplitSummaryTable(summary_branch2)
thickness3, n_exc_3, epsilonFilm3, branch3, root_number3, period3, lspp3 = SplitSummaryTable(summary_branch3)

#print(n_exc)
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
filename = "TiO2Ti-SPPmodes"
np.savetxt(filename+".csv", summary_export_t)
np.savetxt(filename+"-branch0"+".csv", summary_branch0_export_t)
np.savetxt(filename+"-branch1"+".csv", summary_branch1_export_t)
np.savetxt(filename+"-branch2"+".csv", summary_branch2_export_t)
np.savetxt(filename+"-branch3"+".csv", summary_branch3_export_t)

plt.figure()
ax1 = plt.subplot(111)
#plt.title(r'Film thickness $t=$'+str(round(np.real(thickness[0])*1e9))+' nm')
#plt.xlabel(r'Fraction of Cr oxide  (perc.)')
plt.xlabel(r'Excited electron density (m$^{-3}$)') 
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
plot113, = ax1.semilogx(n_exc_3, np.multiply(1e9,period3), 'go', label=r'SPP period $\Lambda$, branch ( +,+)')

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
filename="TiO2Ti-Thickness-"+str(1E9*thickness_max)+"nm"
plt.savefig(filename+".eps")
plt.savefig(filename+".png")
plt.show()
