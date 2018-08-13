#!/usr/bin/env python3
#-*- coding: utf-8 -*-

# Copyright (C) 2018 F. Preucil, T.J.-Y. Derrien
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

## @package libMultilayerSPP
# Module libMultilayerSPP explores the SPP theory at a thin film located 
# between two semi-infinite media. The formal model is presented in 
# T.J.-Y. Derrien et al, J. Appl. Phys. 116, 074902 (2014) and references 
# therein. 

from libMaterials import *

import math, cmath
import numpy             as np
from   scipy.optimize    import root
from   itertools         import product
import matplotlib.pyplot as plt

def func(betaR, eps1, eps2, eps3, k0, t, sgn1, sgn2):
    beta = betaR[0] + betaR[1]*1.j
    kappa1 = cmath.sqrt(beta**2 - ke1)/eps1
    kappa2 = sgn1*cmath.sqrt(beta**2 - ke2)/eps2
    kappa3 = sgn2*cmath.sqrt(beta**2 - ke3)/eps3
    try:
        value = (kappa1-kappa2)*(kappa1-kappa3)*cmath.exp(-2*kappa1*eps1*t)-(kappa1+kappa2)*(kappa1+kappa3)
    except:
        value = 1E99
    return (value.real, value.imag)

def norm(vec):
    return math.sqrt(vec[0]*vec[0] + vec[1]*vec[1])

def cntr(inpt):
    px = 0
    py = 0
    ln = len(inpt)
    for pt in inpt:
        px += pt[0]
        py += pt[1]
    return (px/ln, py/ln)

## Computes the roots of SPP period and propagation length in a 3-layer system for many sample thicknesses
# @param eps1: Dielectric permittivity of the thin film
# @param eps2: Dielectric permittivity of the half-plane below the thin film. 
# @param eps3: Dielectric permittivity of the half-plane above the thin film. Source light is supposed to come from this direction. 
# @param wavelength: wavelength of the source light (SI units). 
# @param t: thickness of the film
# @param x_min: lower boundary of Re(roots)
# @param x_max: higher boundary of Re(roots)
# @param y_min: lower boundary of Im(roots)
# @param y_max: higher boundary of Im(roots)
# @param x_steps: number of steps used to mesh the Re(roots) space. 
# @param y_steps: number of steps used to mesh the Im(roots) space. 
# @param tol_merge: tolerance to merge the identified solutions. 
#Usage:
#
#   findroots(eps1, eps2, eps3,
#             wavelength, t,
#             x_min, x_max,       }
#             y_min, y_max,       } mesh parameters
#             x_steps, y_steps)   }
#
#
#Example:
#roots = findroots(-1+1j, 1, 1,
#                  #1, 1,
#                  #-10, 10,
#                  #-10, 10,
#                  #5, 5,
#                  #.0001)
#
#print(roots) #list of arrays
#print(roots[0])     #prints roots belonging to the zeroth branch (-, -, -)
#print(roots[7][0])  #prints the first root from the last branch (+, +, +)
#
#List of branch indices:
#   0 (-, -, -)
#   1 (-, -, +)
#   2 (-, +, -)
#   3 (-, +, +)
#   4 (+, -, -)
#   5 (+, -, +)
#   6 (+, +, -)
#   7 (+, +, +)
def findroots(eps1, eps2, eps3, wavelength, t, x_min, x_max, y_min, y_max, x_steps, y_steps):
    global k0, ke1, ke2, ke3
    k0 = 2.*np.pi/wavelength
    ke1 = (k0**2)*eps1
    ke2 = (k0**2)*eps2
    ke3 = (k0**2)*eps3

    tol_merge = 1E3
    branches = []
    for sgn1, sgn2 in product((-1,1), (-1,1)):
        roots = []
        unique = []
        for x in np.linspace(x_min, x_max, num=x_steps):
            for y in np.linspace(y_min, y_max, num=y_steps):
                nrt = root(func, (x, y), args=(eps1, eps2, eps3, k0, t, sgn1, sgn2), method='hybr')
                if nrt.success:
                    roots.append(nrt.x)

        ln = len(roots)
        while ln > 0:
            center = roots[0]
            aux = [center]
            aux2 = []
            for rt in roots[1:]:
                if norm(rt-center) < tol_merge:
                    aux.append(rt)
                    center = cntr(aux)
                elif norm(rt+center) < tol_merge:
                    aux.append(-rt)
                    center = cntr(aux)
                else:
                    aux2.append(rt)
            if len(aux) > 4: #merging criterion
                if center[0] > 0:
                    unique.append(center)
                else:
                    unique.append([-center[0], -center[1]])
            roots = list(aux2)
            ln = len(roots)

##        valid = []
##        start = time()
##        for rt in unique:
##            beta = rt[0] + 1.j*rt[1]
##            k1 = cmath.sqrt(beta**2 - ke1)
##            k2 = sgn1*cmath.sqrt(beta**2 - ke2)
##            k3 = sgn2*cmath.sqrt(beta**2 - ke3)
##            C = cmath.exp((-k1-k3)*t/2)*(k1*eps3-k3*eps1)/(2*k1*eps3)
##            D = cmath.exp((k1-k3)*t/2)*(k1*eps3+k3*eps1)/(2*k1*eps3)
##            B1 = C*cmath.exp((k2-k1)*t/2) + D*cmath.exp((k2+k1)*t/2)
##            B2 = (C*cmath.exp((k2-k1)*t/2) - D*cmath.exp((k2+k1)*t/2))*(k1*eps2)/(k2*eps1)
##            if abs(1 - abs(B2/B1)) < tol_valid:
##                valid.append(rt)
##        prnt('Total (validated): %d' % len(valid))
##        print()
        branches.append(unique)

    #selection of the maximum Lspp
    outs = []
    for branch in branches:
        maxy, maxx, maxy2, maxx2 = float('-inf'), float('-inf'), float('-inf'), float('-inf')
        for rt in branch:
            xi = 2.*np.pi/rt[0]
            yi = .5/rt[1]
            if yi > maxy:
                maxy, maxy2 = yi, maxy
                maxx, maxx2 = xi, maxx
            elif yi > maxy2:
                maxy2 = yi
                maxx2 = xi
        if len(branch) == 0:
            outs.append([[0, 0], [0, 0]])
        elif len(branch) == 1:
            outs.append([[maxx, maxy], [0, 0]])
        else:
            outs.append([[maxx, maxy], [maxx2, maxy2]])
    return outs

#end of algorithm

## Repeat results from Charbonneau and Berini, Optics Letters Vol. 25, No. 11 (2000). 
def Berini2000(NumberOfPoints=5): #{{{
    
    #data
    wavelength = 1550e-9
    #epsAu       = -95.95924741872399+10.972582438155513j #1550 nm, Palik
    epsAu       = -131.9475+12.65j
    epsSiO2     = 2.085 #Berini #(1.4440+0j)**2 #Palik

    # Medium 1: thin film. 
    eps1 = epsAu        #thin film
    # Medium 2: substrate. 
    eps2 = epsSiO2       #epsBK7 #environment | substrate
    # Medium 3: environment
    eps3 = epsSiO2       #environment | substrate
    # Note: Inverting eps2 and eps3 should have no effect on the possible modes, but only on field amplification. 

    #t = 100E-9 #thickness of the layer in meters
    #thickness = 100e-9
    t_list = np.power(10., np.linspace(np.log10(1e-9), np.log10(300e-9), NumberOfPoints))
    #meshes the initial guess area, all numbers are from the space of betas
    x_min = -1E10
    x_max = 1E10

    y_min = -1E9
    y_max = 1E9

    x_steps = 40
    y_steps = 40

    summary = np.zeros((0, 5))

    for thickness in t_list:
        roots = findroots(eps1, eps2, eps3,
                wavelength, thickness,
                x_min, x_max,    
                y_min, y_max,    
                x_steps, y_steps)

        # Shaping the data to plot them with GNUplot
        roots_shape = np.shape(roots)
        #print(roots_shape)
        num_thickness = np.shape(thickness) 
        num_branches  = roots_shape[0] #We have 2 or 4 branches
        num_roots     = roots_shape[1] #In each branch, we select 2-3 roots. But there are infinitly of them. 
        num_property  = roots_shape[2] #Period and Lspp

        for branch in np.arange(0,num_branches):
            #thickness_t  = np.zeros((1,0))    
            #spp_period_t = np.zeros((1,0))
            #l_spp_t      = np.zeros((1,0))
            for root_number in np.arange(0,num_roots): 
                ToBeAdded = [thickness, branch, root_number, roots[branch][root_number][0], roots[branch][root_number][1]]
                print(ToBeAdded)
                if(roots[branch][root_number][0] != 0e0): 
                    summary = np.vstack((summary, ToBeAdded ))

    ## === PLOTTING THE RESULTS ==

    # Extract the constructed table
    thickness   = summary[:, 0] 
    branch      = summary[:, 1]
    root_number = summary[:, 2]
    period      = summary[:, 3]
    lspp        = summary[:, 4]

    print(summary)

    # Now, we shall sort out the data
    lists = sorted(zip(*[root_number, branch, thickness, period, lspp]))
    root_number_s, branch_s, thickness_s, period_s, lspp_s = list(zip(*lists))

    ReBeta = np.divide(2.*np.pi,period_s)
    ImBeta = np.divide(0.5,lspp_s)

    plt.figure()
    ax1 = plt.subplot(121)
    plt.xlabel('Film thickness (m)')
    plt.ylabel('SPP period (m)')
    ax12 = ax1.twinx()
    plt.ylabel(r'SPP mean free path $L_{SPP}$ (m)')
    # Then we could plot them in the right order
    plot11,  = ax1.loglog(thickness_s, period_s, 'r+', label=r'Period $\Lambda$')
    plot12,  = ax1.loglog(thickness_s, wavelength*np.ones(np.shape(thickness_s)), '-', label=r'$\lambda$')
    plot13, = ax12.loglog(thickness, lspp_s, 'b+', label=r'$L_{SPP}$')
    plt.tight_layout()

    ax2 = plt.subplot(122)
    plt.xlabel('Film thickness (m)')
    plt.ylabel(r'Re(k) (m$^{-1}$)')
    ax22 = ax2.twinx()
    plt.ylabel(r'Im(k) (m$^{-1}$)')
    # Then we could plot them in the right order
    plot21, = ax2.loglog(thickness_s, ReBeta, '+', label=r'Re($\beta$)')
    #ax2.loglog(thickness_s, wavelength*np.ones(np.shape(thickness_s)), '-', label=r'$\lambda$')
    plot22, = ax22.loglog(thickness, ImBeta, 'b+', label=r'Im($\beta$)')
    plt.legend(loc='best')
    plt.tight_layout()
    plt.savefig("Berini2000.eps")
    plt.show()
#}}}

## Repeats results from Derrien, Bonse et al J. Appl. Phys. (2014), Fig. 6. 
# Air / Si* film / Si
# TODO: add more number of branches to obtain the points that are missing from the article. 
# NOTE: for now, only 2-3 solutions are selected with the following criterion max{L_spp}
def DerrienBonse2014_Air_SiExcited_Si(NumberOfPoints=20): #{{{
    #data
    wavelength = 790e-9 #355e-9 #1030E-9 #1026
    neSi = np.power(10.,np.linspace(27.,29.,NumberOfPoints)) #the most violent change
    #epsTibare   = -6.206969+25.2j #800 nm
    #epsTiO2bare = 7.7841+0.j      #800 nm
    #epsCr       = -0.672122310000001+24.8657476j #-0.67+24.87j    #1026 nm
    #epsBK7      = 2.10277365777   #1026 nm
    #epsCr2O3    = 4.9713+0.1784j  #1 um [JDT Kruschwitz et al, Appl. Opt. 1997]
    epsSi       = (3.693+0.006j)**2
    epsAir      = 1.+0.j          #air
    #epsCu       = -46.6046581932 + 4.7188669976j #1030 nm
    #epsCu       = -1.9937293241+4.9290716854j     #355  nm
    epsH2O      = 1.326**2
    #epsSiO2     = 1.453**2
    
    # Medium 1: thin film. 
    # Varying. See below. 
    # Medium 2: substrate. 
    eps2 = epsSi
    # Medium 3: environment
    eps3 = epsAir       #environment | substrate
    # Note: Inverting eps2 and eps3 should have no effect on the possible modes, but only on field amplification. 

    #t = 100E-9 #thickness of the layer in meters
    #thickness = 10e-9
    t_list = np.array([100e-9, 10e-9]) #, 10e-9, 20e-9])
    #meshes the initial guess area, all numbers are from the space of betas
    x_min = -1E10
    x_max = 1E10

    y_min = -1E9
    y_max = 1E9

    x_steps = 40
    y_steps = 40
    
    summary = np.zeros((0, 7))
    for neSiLocal in neSi:
        eps1 = Drude(wavelength, neSiLocal, epsSi, 1.1e-15**-1, 0.18)     #environment | substrate
        for thickness in t_list:
            roots = findroots(eps1, eps2, eps3,
                    wavelength, thickness,
                    x_min, x_max,    
                    y_min, y_max,    
                    x_steps, y_steps)

            # Shaping the data to plot them with GNUplot
            roots_shape = np.shape(roots)
            #print(roots_shape)
            num_thickness = np.shape(thickness)
            num_branches  = roots_shape[0]
            num_roots     = roots_shape[1] #In each branch, we select 2-3 roots. But there are infinitly of them. 
            num_property  = roots_shape[2] #Period and Lspp
            #num_epsilons  = np.shape(eps3)

            for branch in np.arange(0,num_branches-1):
                for root_number in np.arange(0,num_roots): 
                    ToBeAdded = [thickness, branch, root_number, roots[branch][root_number][0], roots[branch][root_number][1], neSiLocal, eps1]
                    print(ToBeAdded)
                    if(roots[branch][root_number][0] != 0e0): 
                        summary = np.vstack((summary, ToBeAdded ))
                #print("\n")
            
#== Extract the constructed table
    thickness    = summary[:, 0] 
    branch       = summary[:, 1]
    root_number  = summary[:, 2]
    period       = summary[:, 3]
    lspp         = summary[:, 4]
    densitySi    = summary[:, 5]
    epsilonFilm  = summary[:, 6]

    print(summary)

    # Now, we shall sort out the data
    lists = sorted(zip(*[root_number, branch, thickness, period, lspp, densitySi, epsilonFilm]))
    root_number_s, branch_s, thickness_s, period_s, lspp_s, densitySi_s, epsilonFilm_s = list(zip(*lists))

    #ReBeta = np.divide(2.*np.pi,period_s)
    #ImBeta = np.divide(0.5,lspp_s)

    plt.figure()
    ax1 = plt.subplot(111)
    plt.xlabel(r'Density of excited electrons in Si film (m$^{-3}$)')
    plt.ylabel('SPP period (m)')
    ax12 = ax1.twinx()
    plt.ylabel(r'SPP mean free path $L_{SPP}$ (m)')
    # Then we could plot them in the right order
    plot11, = ax1.loglog(densitySi_s, period_s, 'r+', label=r'Period $\Lambda$')
    plot12, = ax1.loglog(densitySi_s, wavelength*np.ones(np.shape(thickness_s)), '-', label=r'$\lambda$')
    plot13, =ax12.loglog(densitySi_s, lspp_s, 'b+', label=r'$L_{SPP}$')
    plt.tight_layout()

    plt.legend(loc='best')
    plt.tight_layout()
    plt.savefig("Derrien2014_Fig5.eps")
    plt.show()
    
#}}}



## Repeats results from Derrien, Bonse et al J. Appl. Phys. (2014), Fig. 6. 
# H20 / Si* film / Si
# TODO: add more number of branches to obtain the points that are missing from the article. 
# NOTE: for now, only 2-3 solutions are selected with the following criterion max{L_spp}
def DerrienBonse2014_H20_SiExcited_Si(NumberOfPoints=20): #{{{
    #data
    wavelength = 790e-9 #355e-9 #1030E-9 #1026
    neSi = np.power(10.,np.linspace(27.,29.,NumberOfPoints)) #the most violent change
    #epsTibare   = -6.206969+25.2j #800 nm
    #epsTiO2bare = 7.7841+0.j      #800 nm
    #epsCr       = -0.672122310000001+24.8657476j #-0.67+24.87j    #1026 nm
    #epsBK7      = 2.10277365777   #1026 nm
    #epsCr2O3    = 4.9713+0.1784j  #1 um [JDT Kruschwitz et al, Appl. Opt. 1997]
    epsSi       = (3.693+0.006j)**2
    epsAir      = 1.+0.j          #air
    #epsCu       = -46.6046581932 + 4.7188669976j #1030 nm
    #epsCu       = -1.9937293241+4.9290716854j     #355  nm
    epsH2O      = 1.326**2
    #epsSiO2     = 1.453**2
    
    # Medium 1: thin film. 
    # Varying. See below. 
    # Medium 2: substrate. 
    eps2 = epsSi
    # Medium 3: environment
    eps3 = epsH2O       #environment | substrate
    # Note: Inverting eps2 and eps3 should have no effect on the possible modes, but only on field amplification. 

    #t = 100E-9 #thickness of the layer in meters
    #thickness = 10e-9
    t_list = np.array([100e-9, 10e-9]) #, 10e-9, 20e-9])
    #meshes the initial guess area, all numbers are from the space of betas
    x_min = -1E10
    x_max = 1E10

    y_min = -1E9
    y_max = 1E9

    x_steps = 40
    y_steps = 40
    
    summary = np.zeros((0, 7))
    for neSiLocal in neSi:
        eps1 = Drude(wavelength, neSiLocal, epsSi, 1.1e-15**-1, 0.18)     #environment | substrate
        for thickness in t_list:
            roots = findroots(eps1, eps2, eps3,
                    wavelength, thickness,
                    x_min, x_max,    
                    y_min, y_max,    
                    x_steps, y_steps)

            # Shaping the data to plot them with GNUplot
            roots_shape = np.shape(roots)
            #print(roots_shape)
            num_thickness = np.shape(thickness)
            num_branches  = roots_shape[0]
            num_roots     = roots_shape[1] #In each branch, we select 2-3 roots. But there are infinitly of them. 
            num_property  = roots_shape[2] #Period and Lspp
            #num_epsilons  = np.shape(eps3)

            for branch in np.arange(0,num_branches-1):
                for root_number in np.arange(0,num_roots): 
                    ToBeAdded = [thickness, branch, root_number, roots[branch][root_number][0], roots[branch][root_number][1], neSiLocal, eps1]
                    print(ToBeAdded)
                    if(roots[branch][root_number][0] != 0e0): 
                        summary = np.vstack((summary, ToBeAdded ))
                #print("\n")
            
#== Extract the constructed table
    thickness    = summary[:, 0] 
    branch       = summary[:, 1]
    root_number  = summary[:, 2]
    period       = summary[:, 3]
    lspp         = summary[:, 4]
    densitySi    = summary[:, 5]
    epsilonFilm  = summary[:, 6]

    print(summary)

    # Now, we shall sort out the data
    lists = sorted(zip(*[root_number, branch, thickness, period, lspp, densitySi, epsilonFilm]))
    root_number_s, branch_s, thickness_s, period_s, lspp_s, densitySi_s, epsilonFilm_s = list(zip(*lists))

    #ReBeta = np.divide(2.*np.pi,period_s)
    #ImBeta = np.divide(0.5,lspp_s)

    plt.figure()
    ax1 = plt.subplot(111)
    plt.xlabel(r'Density of excited electrons in Si film (m$^{-3}$)')
    plt.ylabel('SPP period (m)')
    ax12 = ax1.twinx()
    plt.ylabel(r'SPP mean free path $L_{SPP}$ (m)')
    # Then we could plot them in the right order
    plot11, = ax1.loglog(densitySi_s, period_s, 'r+', label=r'Period $\Lambda$')
    plot12, = ax1.loglog(densitySi_s, wavelength*np.ones(np.shape(thickness_s)), '-', label=r'$\lambda$')
    plot13, =ax12.loglog(densitySi_s, lspp_s, 'b+', label=r'$L_{SPP}$')
    plt.tight_layout()

    plt.legend(loc='best')
    plt.tight_layout()
    plt.savefig("Derrien2014_Fig6.eps")
    plt.show()
    
#}}}


## Repeats results from Derrien, Bonse et al J. Appl. Phys. (2014), Fig. 7. 
# H20* (ne) / SiO2 / Si*
# TODO: add more number of branches to obtain the points that are missing from the article. 
# NOTE: for now, only 2-3 solutions are selected with the following criterion max{L_spp}
def DerrienBonse2014_H20_SiO2_SiExcited(NumberOfPoints=20): #{{{
    #data
    wavelength = 790e-9 #355e-9 #1030E-9 #1026
    neH20 = np.power(10.,np.linspace(27.,29.,NumberOfPoints)) #the most violent change
    #epsTibare   = -6.206969+25.2j #800 nm
    #epsTiO2bare = 7.7841+0.j      #800 nm
    #epsCr       = -0.672122310000001+24.8657476j #-0.67+24.87j    #1026 nm
    #epsBK7      = 2.10277365777   #1026 nm
    #epsCr2O3    = 4.9713+0.1784j  #1 um [JDT Kruschwitz et al, Appl. Opt. 1997]
    epsSi       = (3.693+0.006j)**2
    epsAir      = 1.+0.j          #air
    #epsCu       = -46.6046581932 + 4.7188669976j #1030 nm
    #epsCu       = -1.9937293241+4.9290716854j     #355  nm
    epsH2O      = 1.326**2
    epsSiO2     = 1.453**2
    # Medium 1: thin film. 
    eps1 = epsSiO2
    # Medium 2: substrate. 
    eps2 = Drude(wavelength, 5.1E27, epsSi, 1.1e-15**-1, 0.18)        
    # Medium 3: environment
    #eps3 = epsAir       #environment | substrate
    # Note: Inverting eps2 and eps3 should have no effect on the possible modes, but only on field amplification. 

    #t = 100E-9 #thickness of the layer in meters
    #thickness = 10e-9
    t_list = np.array([5e-9]) #, 10e-9, 20e-9])
    #meshes the initial guess area, all numbers are from the space of betas
    x_min = -1E10
    x_max = 1E10

    y_min = -1E9
    y_max = 1E9

    x_steps = 40
    y_steps = 40
    
    summary = np.zeros((0, 7))
    for neH20Local in neH20:
        eps3 = Drude(wavelength, neH20Local, epsH2O, 1.7e-15**-1, 0.5)      #epsBK7 #environment | substrate
        for thickness in t_list:
            roots = findroots(eps1, eps2, eps3,
                    wavelength, thickness,
                    x_min, x_max,    
                    y_min, y_max,    
                    x_steps, y_steps)

            # Shaping the data to plot them with GNUplot
            roots_shape = np.shape(roots)
            #print(roots_shape)
            num_thickness = np.shape(thickness)
            num_branches  = roots_shape[0]
            num_roots     = roots_shape[1] #In each branch, we select 2-3 roots. But there are infinitly of them. 
            num_property  = roots_shape[2] #Period and Lspp
            #num_epsilons  = np.shape(eps3)

            for branch in np.arange(0,num_branches-1):
                for root_number in np.arange(0,num_roots): 
                    ToBeAdded = [thickness, branch, root_number, roots[branch][root_number][0], roots[branch][root_number][1], neH20Local, eps3]
                    print(ToBeAdded)
                    if(roots[branch][root_number][0] != 0e0): 
                        summary = np.vstack((summary, ToBeAdded ))
                #print("\n")
            
#== Extract the constructed table
    thickness    = summary[:, 0] 
    branch       = summary[:, 1]
    root_number  = summary[:, 2]
    period       = summary[:, 3]
    lspp         = summary[:, 4]
    densityWater = summary[:, 5]
    epsilonWater = summary[:, 6]

    print(summary)

    # Now, we shall sort out the data
    lists = sorted(zip(*[root_number, branch, thickness, period, lspp, densityWater, epsilonWater]))
    root_number_s, branch_s, thickness_s, period_s, lspp_s, densityWater_s, epsilonWater_s = list(zip(*lists))

    #ReBeta = np.divide(2.*np.pi,period_s)
    #ImBeta = np.divide(0.5,lspp_s)

    plt.figure()
    ax1 = plt.subplot(111)
    plt.xlabel(r'Density of excited electrons in H$_2$O (m$^{-3}$)')
    plt.ylabel('SPP period (m)')
    ax12 = ax1.twinx()
    plt.ylabel(r'SPP mean free path $L_{SPP}$ (m)')
    # Then we could plot them in the right order
    plot11, = ax1.loglog(densityWater_s, period_s, 'r+', label=r'Period $\Lambda$')
    plot12, = ax1.loglog(densityWater_s, wavelength*np.ones(np.shape(thickness_s)), '-', label=r'$\lambda$')
    plot13, = ax12.loglog(densityWater_s, lspp_s, 'b+', label=r'$L_{SPP}$')
    plt.tight_layout()

    plt.legend(loc='best')
    plt.tight_layout()
    plt.savefig("Derrien2014_Fig7.eps")
    plt.show()
    
#}}}

## Repeats results from Derrien, Bonse et al J. Appl. Phys. (2014), Fig. 7. 
# H20* (ne) / SiO2 / Si*
# TODO: add more number of branches to obtain the points that are missing from the article. 
# NOTE: for now, only 2-3 solutions are selected with the following criterion max{L_spp}
def DostovalovMetanano2018(NumberOfPoints=20): #{{{
    #data
    wavelength = 1026e-9 #355e-9 #1030E-9 #1026
    epsCr       = -0.672122310000001+24.8657476j #-0.67+24.87j    #1026 nm
    epsBK7      = 2.10277365777   #1026 nm
    #epsCr2O3    = 4.9713+0.1784j  #1 um [JDT Kruschwitz et al, Appl. Opt. 1997]
    #epsSi       = (3.693+0.006j)**2
    epsAir      = 1.+0.j          #air
    # Medium 1: thin film. 
    eps1 = epsSiO2
    # Medium 2: substrate. 
    eps2 = Drude(wavelength, 5.1E27, epsSi, 1.1e-15**-1, 0.18)        
    # Medium 3: environment
    #eps3 = epsAir       #environment | substrate
    # Note: Inverting eps2 and eps3 should have no effect on the possible modes, but only on field amplification. 

    #t = 100E-9 #thickness of the layer in meters
    #thickness = 10e-9
    t_list = np.array([5e-9]) #, 10e-9, 20e-9])
    #meshes the initial guess area, all numbers are from the space of betas
    x_min = -1E10
    x_max = 1E10

    y_min = -1E9
    y_max = 1E9

    x_steps = 40
    y_steps = 40
    
    summary = np.zeros((0, 7))
    for neH20Local in neH20:
        eps3 = Drude(wavelength, neH20Local, epsH2O, 1.7e-15**-1, 0.5)      #epsBK7 #environment | substrate
        for thickness in t_list:
            roots = findroots(eps1, eps2, eps3,
                    wavelength, thickness,
                    x_min, x_max,    
                    y_min, y_max,    
                    x_steps, y_steps)

            # Shaping the data to plot them with GNUplot
            roots_shape = np.shape(roots)
            #print(roots_shape)
            num_thickness = np.shape(thickness)
            num_branches  = roots_shape[0]
            num_roots     = roots_shape[1] #In each branch, we select 2-3 roots. But there are infinitly of them. 
            num_property  = roots_shape[2] #Period and Lspp
            #num_epsilons  = np.shape(eps3)

            for branch in np.arange(0,num_branches-1):
                for root_number in np.arange(0,num_roots): 
                    ToBeAdded = [thickness, branch, root_number, roots[branch][root_number][0], roots[branch][root_number][1], neH20Local, eps3]
                    print(ToBeAdded)
                    if(roots[branch][root_number][0] != 0e0): 
                        summary = np.vstack((summary, ToBeAdded ))
                #print("\n")
            
#== Extract the constructed table
    thickness    = summary[:, 0] 
    branch       = summary[:, 1]
    root_number  = summary[:, 2]
    period       = summary[:, 3]
    lspp         = summary[:, 4]
    densityWater = summary[:, 5]
    epsilonWater = summary[:, 6]

    print(summary)

    # Now, we shall sort out the data
    lists = sorted(zip(*[root_number, branch, thickness, period, lspp, densityWater, epsilonWater]))
    root_number_s, branch_s, thickness_s, period_s, lspp_s, densityWater_s, epsilonWater_s = list(zip(*lists))

    #ReBeta = np.divide(2.*np.pi,period_s)
    #ImBeta = np.divide(0.5,lspp_s)

    plt.figure()
    ax1 = plt.subplot(111)
    plt.xlabel(r'Density of excited electrons in H$_2$O (m$^{-3}$)')
    plt.ylabel('SPP period (m)')
    ax12 = ax1.twinx()
    plt.ylabel(r'SPP mean free path $L_{SPP}$ (m)')
    # Then we could plot them in the right order
    plot11, = ax1.loglog(densityWater_s, period_s, 'r+', label=r'Period $\Lambda$')
    plot12, = ax1.loglog(densityWater_s, wavelength*np.ones(np.shape(thickness_s)), '-', label=r'$\lambda$')
    plot13, = ax12.loglog(densityWater_s, lspp_s, 'b+', label=r'$L_{SPP}$')
    plt.tight_layout()

    plt.legend(loc='best')
    plt.tight_layout()
    plt.savefig("Derrien2014_Fig7.eps")
    plt.show()
    
#}}}

#===============  VALIDATION ON PREVIOUS WORKS ==========
#berini2000() #VALID
#DerrienBonse2014_H20_SiO2_SiExcited() #Incomplete: further roots should be captured to make the complete figure. 
#DerrienBonse2014_H20_SiExcited_Si()   #More than complete. NOTE: Novel results are available, we can already find transition betweel plasmon polaritons and phonons polaritons by using discontinuity of L_spp. 
#DerrienBonse2014_Air_SiExcited_Si()   # More than complete as well. TODO: show to Jörn !!!!

# =============== NOVEL CONCLUSIONS ==========
#DostovalovMetanano2018()

    ##TODO: from this, we would like to add a layer which will variate eps1 as function of oxide concentration
