#!/usr/bin/env python2
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

## @package midIRxp
# Preparation of results for Prof. Bulgakova mid-IR experiments. 

import numpy as np
from libKeldysh import *

wavelength = [1600e-9, 1700e-9, 2200e-9, 2600e-9, 3200e-9] #available wavelengths
tau = 2e-12 #mid-IR laser HiLASE
w0 = 10e-6 #diameter (m)
EnergyPerPulse = 50e-6 #maximum

PhotonEnergy=np.divide(h*c,wavelength)
E = PhotonEnergy / e

def BeamProperties(): 
    print("Wavelengths (m)")
    print(wavelength)
    print("Photon Energies (eV)")
    print(E)
    print("Duration for 1 cycle (s)")
    omega = np.divide(2.*np.pi*c,wavelength)
    T = np.divide(2.*np.pi, omega)
    print(T)
    print("Photon frequency (Hz)")
    print(np.divide(2.*np.pi,T))
    print("Number of cycles in each pulse")
    print(np.divide(tau, T))

def Si():
    Egap_d = 3.4*e
    Egap_i = 1.12*e
    meff = 1.0 #which effective mass should we take ?
    material = "ZnO"
    return Egap_d, Egap_i, meff
    
def Ge(): 
    Egap_d = 0.8*e
    Egap_i = 0.66*e
    meff = 1.0
    material = "ZnO"
    return Egap_d, Egap_i, meff, material

def ZnO(): 
    Egap_d = 3.6*e
    Egap_i = 3.6*e
    meff = 1.0
    material = "ZnO"
    return Egap_d, Egap_i, meff, material
    
BeamProperties()
#Egap_d, Egap_i, meff, material = Si()
#Egap_d, Egap_i, meff, material = Ge()
Egap_d, Egap_i, meff, material = ZnO()

Nphotons_d=np.ceil(Egap_d / PhotonEnergy)
Nphotons_i=np.ceil(Egap_i / PhotonEnergy)
wavelength_index = 0 #1, 2
numpoints = 10
    
def BeamIntensity(): 
    print("Peak fluence (J/cm2)")
    PeakFluence = 2.*EnergyPerPulse / np.pi / w0**2
    print(PeakFluence*1e-4)
    print("Peak intensity (J/m2)")
    Ipeak = PeakFluence / tau * np.sqrt(4.*np.log(2.)/np.pi)
    print(Ipeak)
    print("Peak field (V/nm)")
    Efield = np.sqrt(2.*Ipeak / c / epsilon_0) #V/m
    print(Efield/1e9)
    return PeakFluence, Ipeak, Efield

PeakFluence, Ipeak, Efield = BeamIntensity()

print("Keldysh Adiabadicity coefficients")
print(gammaKeldysh(Egap_d, meff, Efield, wavelength))

print("=== Intensity dependent results ===")
Ipeak_max = 14.95E16 #W/m2
Efield_max = np.sqrt(2.*Ipeak / c / epsilon_0)
Efields_log = np.linspace(6, np.log10(Efield_max), numpoints)
Efields     = np.power(10., Efields_log)
Ipeak_list = 0.5*c*epsilon_0*np.power(Efields,2) * 1E-4

#colors={'r', 'b', 'g'}

SizeX = 8
SizeY = SizeX / ((1.+np.sqrt(5.))/2.)

plt.figure(figsize=(SizeX,SizeY))
#plt.title(r"$\lambda$="+str(wavelength[wavelength_index]*1E9)+" nm")
plt.xlabel(r"Peak intensity (W/cm$^2$)")
plt.ylabel(r"Keldysh parameter $\gamma$")
for wavelength_index in [0, 1, 2, 3, 4]: 
    gammas = gammaKeldysh(Egap_d, meff, Efields, wavelength[wavelength_index])
    plt.loglog(Ipeak_list, gammas, label=r"$\gamma(\lambda=$"+str(wavelength[wavelength_index]*1E9)+" nm)")
    
plt.loglog(Ipeak_list, np.ones(np.shape(Ipeak_list))*0.1, 'k-' ) #, label=r"$\gamma=0.1$")
plt.loglog(Ipeak_list, np.ones(np.shape(Ipeak_list))*1,   'k--') #, label=r"$\gamma=1$")
plt.loglog(Ipeak_list, np.ones(np.shape(Ipeak_list))*10,  'k-.') #, label=r"$\gamma=10$")
plt.legend(loc="upper left")
plt.tight_layout()
plt.savefig(material+"-mid-IR-KeldyshParameters.png")
plt.show()
