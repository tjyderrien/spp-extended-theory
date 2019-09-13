#!/usr/bin/env python2
#-*- coding: utf-8 -*-

# Copyright (C) 2013-2019 T. J.-Y. Derrien
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

## @package Starinskiy
# We explore here the transmission spectrum of: 
# - a thin film with nanoparticles included using effective medium theory, 
# - and compare with the bulk model of Mie scattering in a dielectric matrix. 
# We have 6 systems: 
# 1-3: vacuum | nanoparticles (7 nm, 10 nm, 12 nm) in vacuum | Si substrates
# 4-6: vacuum | SiOx (with 6%, 8%, 10% volume of nanoparticles (7, 10, 12 nm diameter) | Si substrates

import numpy as np
from scipy.interpolate import InterpolatedUnivariateSpline
from libMaterials import MaxwellGarnett2
from libDatabase import ExportToTxt

epsAir=1.

## Compute the dielectric permittivity of the film containing the set of nanoparticles. 
# This will be introduced into the Mie scattering theory, as a background, to somehow reconstruct the effect of the periodic boundary conditions. 
def BackgroundForMie(order=1):
    # Import optical data of the SiOx.
    filename="Starinskiy_Properties_SiOx.csv"
    SiOx_nk = np.loadtxt(filename)
    wavelength_t = SiOx_nk[:, 0]*1E-9 #in m
    n_t          = SiOx_nk[:, 1]
    k_t          = SiOx_nk[:, 2]
    # Interpolate between the points
    n_func = InterpolatedUnivariateSpline(wavelength_t, n_t, k=order) #wavelength in m
    k_func = InterpolatedUnivariateSpline(wavelength_t, k_t, k=order)
    #print n_func(500e-9) #WORKS.
    #print k_func(500e-9) #WORKS.
    
    # Import (n,k)[wavelength] spectrum of bulk Ag from Palik data
    filename = "/home/hilase/Documents/spp-extended-theory/SimpleSPP/Database/Ag-Johnson"
    nk_NP_t  = np.loadtxt(filename, skiprows=4)
    wavelength_NP_t = nk_NP_t[:,0]*1E-6
    n_NP_t   = nk_NP_t[:,1]
    k_NP_t   = nk_NP_t[:,2]
    n_NP_func = InterpolatedUnivariateSpline(wavelength_NP_t, n_NP_t, k=order)
    k_NP_func = InterpolatedUnivariateSpline(wavelength_NP_t, k_NP_t, k=order)
    
    #print n_NP_t 
    #print n_NP_func(wavelength_NP_t) #WORKS
    
    # We interpolate the freshly captured data from Johnson/Palik/other on the wavelenth of interest
    
    n_NP_t_new = n_NP_func(wavelength_t) #WORKS
    k_NP_t_new = k_NP_func(wavelength_t) #WORKS
    
    #print n_NP_t, n_NP_t_new #WORKS
    #print k_NP_t, k_NP_t_new #WORKS
    
    eps_NP_t = np.power(np.add(n_NP_t_new, np.multiply(1.j, k_NP_t_new)), 2) #(n+ik)**2
    #print n_NP_t+1j*k_NP_t, np.sqrt(eps_NP_t) #WORKS
        
    # We mix nanoparticles with the air, first. 
    # BUG: did Sergey Starinskiy renormalized the volume fraction before deducing it? After discussion, it should clearly be recomputed. 
    
    eps_AirNP_7nm  = MaxwellGarnett2(epsAir, eps_NP_t, 0.06)
    eps_AirNP_10nm = MaxwellGarnett2(epsAir, eps_NP_t, 0.08)
    eps_AirNP_12nm = MaxwellGarnett2(epsAir, eps_NP_t, 0.10)
    
    # Converting optical index to dielectric permittivity
    eps_SiOx_t = np.power(np.add(n_t,np.multiply(1j,k_t)), 2)
    
    # We mix nanoparticles with SiOx
    # BUG: did Sergey Starinskiy renormalized the volume fraction before deducing it? 
    eps_SiOxNP_7nm  = MaxwellGarnett2(eps_SiOx_t, eps_NP_t, 0.06) #BUG: we generate negative imaginary parts ?! Error in MaxwellGarnett2 ?
    eps_SiOxNP_10nm = MaxwellGarnett2(eps_SiOx_t, eps_NP_t, 0.08)
    eps_SiOxNP_12nm = MaxwellGarnett2(eps_SiOx_t, eps_NP_t, 0.10)
    
    # Converting the dielectric permittivites back to N, K for Sergey
    nk_AirNP_7nm   = np.sqrt(eps_AirNP_7nm)
    nk_AirNP_10nm  = np.sqrt(eps_AirNP_10nm)
    nk_AirNP_12nm  = np.sqrt(eps_AirNP_12nm)
    
    nk_SiOxNP_7nm   = np.sqrt(eps_SiOxNP_7nm)
    nk_SiOxNP_10nm  = np.sqrt(eps_SiOxNP_10nm)
    nk_SiOxNP_12nm  = np.sqrt(eps_SiOxNP_12nm)
    
    # We write the final contents into CSV files. 
    Sample1="nk_AirNP_7nm.csv"
    Sample2="nk_AirNP_10nm.csv"
    Sample3="nk_AirNP_12nm.csv"
             
    Sample4="nk_SiOxNP_7nm.csv"
    Sample5="nk_SiOxNP_10nm.csv"
    Sample6="nk_SiOxNP_12nm.csv"
    
    ExportToTxt([1E9*wavelength_NP_t, nk_AirNP_7nm.real,   nk_AirNP_7nm.imag  ] , Sample1)
    ExportToTxt([1E9*wavelength_NP_t, nk_AirNP_10nm.real , nk_AirNP_10nm.imag ], Sample2)
    ExportToTxt([1E9*wavelength_NP_t, nk_AirNP_12nm.real , nk_AirNP_12nm.imag ], Sample3)
    ExportToTxt([1E9*wavelength_NP_t, nk_SiOxNP_7nm.real , nk_SiOxNP_7nm.imag ], Sample4)
    ExportToTxt([1E9*wavelength_NP_t, nk_SiOxNP_10nm.real, nk_SiOxNP_10nm.imag], Sample5)
    ExportToTxt([1E9*wavelength_NP_t, nk_SiOxNP_12nm.real, nk_SiOxNP_12nm.imag], Sample6)
    
    return 0

BackgroundForMie()
    
