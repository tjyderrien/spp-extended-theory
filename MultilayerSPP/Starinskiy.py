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
import matplotlib.pyplot as plt
from scipy.interpolate import InterpolatedUnivariateSpline

from libMaterials import MaxwellGarnett2, Bruggeman2, BiLayerTransmission
from libDatabase import ExportToTxt

epsAir=1.

## Compute the volume fraction of nanoparticles in a film assuming measurement was performed on a square area
# @param t: film thickness (meters)
# @param r: radius of nanoparticles
# @param d_ref: side of the square area employed for nanoparticle count
# @param N: number of nanoparticles on the area of measurement
def VolumeFraction_SquareArea(thickness, radius, d_ref, N):
    diameter=2.*radius
    volume_fraction=4./3.*N*np.pi*(diameter/2.)**3/(d_ref**2 * thickness)
    return volume_fraction

VolumeFraction_SquareArea=np.vectorize(VolumeFraction_SquareArea)

## Import the optical data from the set provided by Sergey
def ImportOpticalDataForSiOxFilm(order=2): 
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
    #filename = "/home/hilase/Documents/spp-extended-theory/SimpleSPP/Database/Au-Johnson" #just changing 1 single character and you get Gold instead of Ag... 
    nk_NP_t  = np.loadtxt(filename, skiprows=4)
    wavelength_NP_t = nk_NP_t[:,0]*1E-6
    n_NP_t   = nk_NP_t[:,1]
    k_NP_t   = nk_NP_t[:,2]
    n_NP_func = InterpolatedUnivariateSpline(wavelength_NP_t, n_NP_t, k=order)
    k_NP_func = InterpolatedUnivariateSpline(wavelength_NP_t, k_NP_t, k=order)
    
    # Import(n,k)[wavelength] spectrum of bulk SiO2
    filename = "/home/hilase/Documents/spp-extended-theory/SimpleSPP/Database/SiO2-Palik"
    nk_SiO2_t  = np.loadtxt(filename, skiprows=4)
    wavelength_SiO2_t = nk_SiO2_t[:,0]*1E-10
    n_SiO2_t   = nk_SiO2_t[:,1]
    k_SiO2_t   = nk_SiO2_t[:,2]
    n_SiO2_func = InterpolatedUnivariateSpline(wavelength_SiO2_t, n_SiO2_t, k=order)
    k_SiO2_func = InterpolatedUnivariateSpline(wavelength_SiO2_t, k_SiO2_t, k=order)
    
    #print n_NP_t 
    #print n_NP_func(wavelength_NP_t) #WORKS
    
    # New wavelength mesh will be more precise, but is bound to minimum overlapping wavelength interval
    # NewSet = Set1 inter Set2
    wavelength_t_min = np.max([wavelength_NP_t.min(), wavelength_t.min(), wavelength_SiO2_t.min()])
    wavelength_t_max = np.min([wavelength_NP_t.max(), wavelength_t.max(), wavelength_SiO2_t.max()])
    wavelength_new_t = np.linspace(wavelength_t_min, wavelength_t_max, 2000) #OK
    #print(wavelength_t, wavelength_new_t) #OK
    
    n_SiOx_new=n_func(wavelength_new_t) #OK
    k_SiOx_new=k_func(wavelength_new_t) #OK
    
    #print(n_t, n_SiOx_new) #OK
    #print(k_t, k_SiOx_new) #OK
    
    # We interpolate the freshly captured data from Johnson/Palik/other on the wavelenth of interest
    n_NP_t_new = n_NP_func(wavelength_new_t)
    k_NP_t_new = k_NP_func(wavelength_new_t)
    
    # Interpolating the SiO2 data from Palik
    n_SiO2_t_new = n_SiO2_func(wavelength_new_t)
    k_SiO2_t_new = k_SiO2_func(wavelength_new_t)
    
    #print n_NP_t, n_NP_t_new #OK
    #print k_NP_t, k_NP_t_new #OK
    # Converting optical index to dielectric permittivity
    #eps_SiOx_t = np.power(np.add(n_t,np.multiply(1j,k_t)), 2)
    eps_SiOx_t = np.power(np.add(n_SiOx_new,np.multiply(1j,k_SiOx_new)), 2)
    eps_NP_t   = np.power(np.add(n_NP_t_new, np.multiply(1.j, k_NP_t_new)), 2) #(n+ik)**2
    eps_SiO2_t = np.power(np.add(n_SiO2_t_new, np.multiply(1.j, k_SiO2_t_new)), 2)
    #print n_NP_t+1j*k_NP_t, np.sqrt(eps_NP_t) #WORKS
    return wavelength_new_t, eps_NP_t, eps_SiOx_t, eps_SiO2_t


## Compute the dielectric permittivity of the film containing the set of nanoparticles. 
# This will be introduced into the Mie scattering theory, as a background, to somehow reconstruct the effect of the periodic boundary conditions. 
def BackgroundForMie(): #{{{
    wavelength_new_t, eps_NP_t, eps_SiOx_t, eps_SiO2_t = ImportOpticalDataForSiOxFilm()
    
    # Recovering n, k from eps_NP_t
    nk_NP_t_new = np.sqrt(eps_NP_t)
    n_NP_t_new = nk_NP_t_new.real
    k_NP_t_new = nk_NP_t_new.imag
    
    nk_SiOx_new = np.sqrt(eps_SiOx_t)
    n_SiOx_new  = nk_SiOx_new.real
    k_SiOx_new  = nk_SiOx_new.imag
    
    # We mix nanoparticles with the air, first. 
    # BUG: did Sergey Starinskiy renormalized the volume fraction before deducing it? After discussion, it should clearly be recomputed. 
    
    ## Now we compute spectrum of several configurations
    ## CONFIG 1: air  |(NP mixed with air ) | SiO2 substrate
    ## CONFIG 2: SiOx |(NP mixed with SiOx) | SiO2 substrate (denom. "strict")
    ## CONFIG 3: air  |(NP mixed with SiOx) | SiO2 substrate
    ## 
    ## We also consider configurations A, B, C, as no
    
    ## OPTICAL DATA FOR CONFIG 1. 
    ## NOTE: we can use MaxwellGarnett2 or Bruggeman2 theory. 
    eps_AirNP_7nm  = MaxwellGarnett2(epsAir, eps_NP_t, 0.1795)
    eps_AirNP_10nm = MaxwellGarnett2(epsAir, eps_NP_t, 0.2994)
    eps_AirNP_12nm = MaxwellGarnett2(epsAir, eps_NP_t, 0.3533)
    
    epsAir_t = np.ones(np.shape(eps_NP_t)) #dirty fix to make zip function happy. 
    #print epsAir_t
    
    ## NOTE: we now try to build Set#3 of data with Bruggeman theory. 
    #eps_AirNP_7nm_Bruggeman_1, eps_AirNP_7nm_Bruggeman_2 = Bruggeman2(epsAir, eps_NP_t, 0.1795)
    
    #print "== ATTEMPT WITH BRUGGERMAN =="
    #print eps_AirNP_7nm_Bruggeman_1, eps_AirNP_7nm_Bruggeman_2
    
    ## PREPARING OPTICAL DATA FOR CONFIG 2. 
    eps_SiOxNP_7nm_strict  = MaxwellGarnett2(eps_SiOx_t, eps_NP_t, 0.1795)
    eps_SiOxNP_10nm_strict = MaxwellGarnett2(eps_SiOx_t, eps_NP_t, 0.2994)
    eps_SiOxNP_12nm_strict = MaxwellGarnett2(eps_SiOx_t, eps_NP_t, 0.3533)
    
    ## PREPARING OPTICAL DATA FOR CONFIG 3. 
    # We mix nanoparticles with SiOx, with complete SiOx thickness
    eps_SiOxNP_7nm  = MaxwellGarnett2(eps_SiOx_t, eps_NP_t, 0.06)
    eps_SiOxNP_10nm = MaxwellGarnett2(eps_SiOx_t, eps_NP_t, 0.08)
    eps_SiOxNP_12nm = MaxwellGarnett2(eps_SiOx_t, eps_NP_t, 0.10)
    
    ## EXPORTING THE DATA FOR SERGEY
    # Converting the dielectric permittivites back to N, K for Sergey
    nk_AirNP_7nm   = np.sqrt(eps_AirNP_7nm)
    nk_AirNP_10nm  = np.sqrt(eps_AirNP_10nm)
    nk_AirNP_12nm  = np.sqrt(eps_AirNP_12nm)
    
    nk_SiOxNP_7nm_strict  = np.sqrt(eps_SiOxNP_7nm_strict )
    nk_SiOxNP_10nm_strict = np.sqrt(eps_SiOxNP_10nm_strict)
    nk_SiOxNP_12nm_strict = np.sqrt(eps_SiOxNP_12nm_strict)
    
    nk_SiOxNP_7nm   = np.sqrt(eps_SiOxNP_7nm )
    nk_SiOxNP_10nm  = np.sqrt(eps_SiOxNP_10nm)
    nk_SiOxNP_12nm  = np.sqrt(eps_SiOxNP_12nm)
    
    ## OPTICAL DATA FOR SUBSTRATE
    #eps_SiO2_t
    
    ## COMPUTING TRANMISSION SPECTRA
    ## 
    ## CONFIG 1, 3 NP sizes
    ## air | air + Ag NP (t=t_SiOx (exp.)) | SiO2
    
    film_thickness = 30E-9 #experimental value given by Sergey
    
    Tomega1_7nm  = BiLayerTransmission(wavelength_new_t, epsAir_t, eps_AirNP_7nm,  eps_SiO2_t, film_thickness) #Sample1
    Tomega1_10nm = BiLayerTransmission(wavelength_new_t, epsAir_t, eps_AirNP_10nm, eps_SiO2_t, film_thickness) #Sample2
    Tomega1_12nm = BiLayerTransmission(wavelength_new_t, epsAir_t, eps_AirNP_12nm, eps_SiO2_t, film_thickness) #Sample3
    
    ## CONFIG 2, 3 NP sizes
    ## SiOx | SiOx + Ag NP (t=t_NP strict) | SiO2
    Tomega2_7nm  = BiLayerTransmission(wavelength_new_t, eps_SiOx_t, eps_SiOxNP_7nm_strict,  eps_SiO2_t, 7e-9) #Sample7
    Tomega2_10nm = BiLayerTransmission(wavelength_new_t, eps_SiOx_t, eps_SiOxNP_10nm_strict, eps_SiO2_t, 10E-9)#Sample8
    Tomega2_12nm = BiLayerTransmission(wavelength_new_t, eps_SiOx_t, eps_SiOxNP_12nm_strict, eps_SiO2_t, 12E-9)#Sample9
    
    ## CONFIG 3, 3 NP sizes
    ## air | SiOx + Ag NP (t=t_SiOx experimental) | SiO2 
    Tomega3_7nm  = BiLayerTransmission(wavelength_new_t, epsAir_t, eps_SiOxNP_7nm,  eps_SiO2_t, film_thickness) #Sample4
    Tomega3_10nm = BiLayerTransmission(wavelength_new_t, epsAir_t, eps_SiOxNP_10nm, eps_SiO2_t, film_thickness) #Sample5
    Tomega3_12nm = BiLayerTransmission(wavelength_new_t, epsAir_t, eps_SiOxNP_12nm, eps_SiO2_t, film_thickness) #Sample6
    
    print "** DEBUG **"
    print Tomega1_7nm
    
    ## Plotting the final results. 
    plt.figure()
    ax1=plt.subplot(311)
    #ax1.set_xlabel("Light wavelength (nm)")
    #ax1.set_ylabel("Transmission")
    
    ax2=plt.subplot(312)
    #ax2.set_xlabel("Light wavelength (nm)")
    #ax2.set_ylabel("Transmission")
    
    ax3=plt.subplot(313)
    ax3.set_xlabel("Light wavelength (nm)")
    ax3.set_ylabel("Transmission")
    
    wavelength_new_t_show = np.multiply(1E9, wavelength_new_t)
    
    ax1.plot(wavelength_new_t_show, Tomega1_7nm , "-" ,  label=r"air | 30 nm [air + (NP 7 nm)] | SiO$_2$")
    ax1.plot(wavelength_new_t_show, Tomega2_7nm , "--",  label=r"SiO$_x$ | (SiO$_x$ + NP) 7 nm | SiO$_2$")
    ax1.plot(wavelength_new_t_show, Tomega3_7nm , "-.",  label=r"air | 30 nm [SiO$_x$ + (NP 7 nm)] | SiO$_2$")
    ax1.legend(loc="best")
        
    ax2.plot(wavelength_new_t_show, Tomega1_10nm, "-" , label=r"air | 30 nm [air + (NP 10 nm)] | SiO$_2$")
    ax2.plot(wavelength_new_t_show, Tomega2_10nm, "--", label=r"SiO$_x$ | (SiO$_x$ + NP) 10 nm | SiO$_2$")
    ax2.plot(wavelength_new_t_show, Tomega3_10nm, "-.", label=r"air | 30 nm [SiO$_x$ + (NP 10 nm)] | SiO$_2$")
    ax2.legend(loc="best")
    
    ax3.plot(wavelength_new_t_show, Tomega1_12nm, "-" , label=r"air | 30 nm [air + (NP 12 nm)] | SiO$_2$")
    ax3.plot(wavelength_new_t_show, Tomega2_12nm, "--", label=r"SiO$_x$ | (SiO$_x$ + NP) 12 nm | SiO$_2$")
    ax3.plot(wavelength_new_t_show, Tomega3_12nm, "-.", label=r"air | 30 nm [SiO$_x$ + (NP 12 nm)] | SiO$_2$")
    ax3.legend(loc="best")
    
    filename="Starinskiy_AgNP"
    
    plt.tight_layout()
    
    plt.savefig(filename+".eps")
    plt.savefig(filename+".png")
    plt.show()
    
    ## We export write the final contents into CSV files, for Sergey. 
    Sample1="Air-AirAndAgNP7nm-30nmThick-SiO2.csv"
    Sample2="Air-AirAndAgNP10nm-30nmThick-SiO2.csv"
    Sample3="Air-AirAndAgNP12nm-30nmThick-SiO2.csv"
             
    Sample4="Air-SiOxAndAgNP7nm-30nmThick-SiO2.csv"
    Sample5="Air-SiOxAndAgNP10nm-30nmThick-SiO2.csv"
    Sample6="Air-SiOxAndAgNP12nm-30nmThick-SiO2.csv"
    
    Sample7="SiOx-SiOxAndAgNP7nm-7nmThick-SiO2.csv"
    Sample8="SiOx-SiOxAndAgNP10nm-10nmThick-SiO2.csv"
    Sample9="SiOx-SiOxAndAgNP12nm-12nmThick-SiO2.csv"
    
    header="wavelength (nm)\tn(matrix)\tk(matrix)\tn(NP)\tk(NP)\tn(eff)\tk(eff)\tTransmission"
    
    Sample1_Array=list(zip(*[np.real(1E9*wavelength_new_t), np.real(epsAir_t),   np.imag(epsAir_t)  , np.real(n_NP_t_new), np.real(k_NP_t_new), np.real(nk_AirNP_7nm  ), np.imag(nk_AirNP_7nm  ), np.real(Tomega1_7nm )]))
    Sample2_Array=list(zip(*[np.real(1E9*wavelength_new_t), np.real(epsAir_t),   np.imag(epsAir_t)  , np.real(n_NP_t_new), np.real(k_NP_t_new), np.real(nk_AirNP_10nm ), np.imag(nk_AirNP_10nm ), np.real(Tomega1_10nm)]))
    Sample3_Array=list(zip(*[np.real(1E9*wavelength_new_t), np.real(epsAir_t),   np.imag(epsAir_t)  , np.real(n_NP_t_new), np.real(k_NP_t_new), np.real(nk_AirNP_12nm ), np.imag(nk_AirNP_12nm ), np.real(Tomega1_12nm)]))
    Sample4_Array=list(zip(*[np.real(1E9*wavelength_new_t), n_SiOx_new,          k_SiOx_new,          np.real(n_NP_t_new), np.real(k_NP_t_new), np.real(nk_SiOxNP_7nm ), np.imag(nk_SiOxNP_7nm ), np.real(Tomega3_7nm )]))
    Sample5_Array=list(zip(*[np.real(1E9*wavelength_new_t), n_SiOx_new,          k_SiOx_new,          np.real(n_NP_t_new), np.real(k_NP_t_new), np.real(nk_SiOxNP_10nm), np.imag(nk_SiOxNP_10nm), np.real(Tomega3_10nm)]))
    Sample6_Array=list(zip(*[np.real(1E9*wavelength_new_t), n_SiOx_new,          k_SiOx_new,          np.real(n_NP_t_new), np.real(k_NP_t_new), np.real(nk_SiOxNP_12nm), np.imag(nk_SiOxNP_12nm), np.real(Tomega3_12nm)]))
    Sample7_Array=list(zip(*[np.real(1E9*wavelength_new_t), n_SiOx_new,          k_SiOx_new,          np.real(n_NP_t_new), np.real(k_NP_t_new), np.real(nk_SiOxNP_7nm_strict),  np.imag(nk_SiOxNP_7nm_strict),  np.real(Tomega2_7nm )]))
    Sample8_Array=list(zip(*[np.real(1E9*wavelength_new_t), n_SiOx_new,          k_SiOx_new,          np.real(n_NP_t_new), np.real(k_NP_t_new), np.real(nk_SiOxNP_10nm_strict), np.imag(nk_SiOxNP_10nm_strict), np.real(Tomega2_10nm)]))
    Sample9_Array=list(zip(*[np.real(1E9*wavelength_new_t), n_SiOx_new,          k_SiOx_new,          np.real(n_NP_t_new), np.real(k_NP_t_new), np.real(nk_SiOxNP_12nm_strict), np.imag(nk_SiOxNP_12nm_strict), np.real(Tomega2_12nm)]))
    
    np.savetxt(Sample1, Sample1_Array, fmt='%1.4e', header=header, delimiter=",")
    np.savetxt(Sample2, Sample2_Array, fmt='%1.4e', header=header, delimiter=",")
    np.savetxt(Sample3, Sample3_Array, fmt='%1.4e', header=header, delimiter=",")
    np.savetxt(Sample4, Sample4_Array, fmt='%1.4e', header=header, delimiter=",")
    np.savetxt(Sample5, Sample5_Array, fmt='%1.4e', header=header, delimiter=",")
    np.savetxt(Sample6, Sample6_Array, fmt='%1.4e', header=header, delimiter=",")
    np.savetxt(Sample7, Sample7_Array, fmt='%1.4e', header=header, delimiter=",")
    np.savetxt(Sample8, Sample8_Array, fmt='%1.4e', header=header, delimiter=",")
    np.savetxt(Sample9, Sample9_Array, fmt='%1.4e', header=header, delimiter=",")
    
    return 0
#}}}

#Data provided by Sergey
diameters=np.array([7e-9,10e-9,12e-9])
radius=0.5*diameters #m
area_measurement=1E-12 #m² (1um x 1um)
Number_of_particles_per_area=[7000,5718,4686]

# Example
film_thickness = 30e-9 #Experiment was made with 30 nm. 

# Logical
d_ref = np.sqrt(area_measurement)
print "== VOLUME FRACTION =="
print "Nanoparticles in air, no film"
# VolumeFraction_SquareArea(thickness, radius, d_ref, N)
print VolumeFraction_SquareArea(2.*radius, radius, d_ref, Number_of_particles_per_area)

print "Nanoparticles in SiOx film (warning: volume fraction will vary with film thickness)"
print VolumeFraction_SquareArea(film_thickness, radius, d_ref, Number_of_particles_per_area)
print "Warning: these values don't take into account the particle distribution."

print "Nanoparticles in a SiOx film of thickness settled by the nanoparticle average diameter."
print VolumeFraction_SquareArea(2*radius, radius, d_ref, Number_of_particles_per_area)
print "Warning: these values don't take into account the particle distribution."

print "== EFFECTIVE OPTICAL DATA =="
BackgroundForMie()
print "See .csv files."
