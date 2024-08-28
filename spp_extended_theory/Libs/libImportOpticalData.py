#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Copyright (C) 2013-2024 T. J.-Y. Derrien
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

## @package importPalikData
# Importing data from Palik book using digitized plots. 
# Optical data are input via using CSV-formatted input files, created using Engauge-digitizer software. 
# The lib generates (wavelength, ReEps, ImEps) values to be used inside the program.
# Two possibilities are available: 
# - Interpolating one single value at a given wavelength (importFromTable). 
# - Interpolating on a wide spectrum to combine several sources (other functions). 

# IMPORT LIBRARIES
from spp_extended_theory.Libs.libSPP import *
# from spp_extended_theory.Keldysh.libUnits import *
import spp_extended_theory.Keldysh.libUnits as libUnits
from spp_extended_theory.Libs.libDatabase import *
import numpy as np
from scipy.interpolate import InterpolatedUnivariateSpline
import matplotlib.pyplot as plt
from octopus_slabs.Libs import libLogging

logger = libLogging.init_logger(__name__, verbose=True)

# from libImportOpticalData import *

Header = "[importOpticalData] "

## Imports 2 files: one for N, one for K. Does the interpolation of wavelength mesh, then returns table of complex epsilon.
# @param wavelength: wavelength in meters
# @param folder: location of optical data
# @param filename: filename before -n and -k.csv.
# @param plotting=1: show plot?
# @param unit=1E-6: basic unit is in um.
def importFromNKtable(wavelength, folder, filename, plotting=1, unit=1E-6):  # {{{
    nfile = folder + filename + "-n.csv"
    kfile = folder + filename + "-k.csv"

    narray = np.loadtxt(nfile, delimiter="\t", skiprows=1)
    karray = np.loadtxt(kfile, delimiter="\t", skiprows=1)

    # unit = 1E-6

    # import wavelength, n and k from Palik
    wavelength1 = narray[:, 0] * unit
    wavelength2 = karray[:, 0] * unit
    n = narray[:, 1]
    k = karray[:, 1]
    numrows = 10000
    base = 10
    # interpolate n and k on new wavelength mesh
    order = 1
    # wavelengths = np.arange(np.amin(wavelength2),np.amax(wavelength2), precision) #regular mesh, AWFUL for memory
    print(("importFromNKtable: Generating new wavelength mesh: (" + str(np.amin(wavelength2)) + ", " + str(
        np.amax(wavelength2)) + ")"))
    wavelengths = np.logspace(np.amin(np.log10(wavelength2)), np.amax(np.log10(wavelength2)), num=numrows, base=base,
                              endpoint=True)

    print(("New wavelength mesh has " + str(numrows) + " rows."))
    # print wavelengths

    fni = InterpolatedUnivariateSpline(wavelength1, n, k=order)
    fki = InterpolatedUnivariateSpline(wavelength2, k, k=order)

    # Interpolated one optical constants
    try:
        # wavelength = 800e-9
        ni = fni(wavelength);
        ki = fki(wavelength)
        epsilon = (ni + 1j * ki) ** 2
        print(("importFromNKtable: Interpolated permittivity at " + str(wavelength * 1e9) + " nm = " + str(epsilon)))
    except:
        print(("importFromNKtable: Interpolation for " + str(wavelength * 1E9) + " nm failed."))

    # try:
    # wavelength = 532e-9
    # ni = fni(wavelength); ki = fki(wavelength)
    # epsilon = (ni+1j*ki)**2
    # print "Interpolated permittivity at "+str(wavelength*1e9)+" nm = "+str(epsilon)
    # except:
    # print "Interpolation for "+str(wavelength*1E9)+" nm failed."
    # try:
    # wavelength = 400e-9
    # ni = fni(wavelength); ki = fki(wavelength)
    # epsilon = (ni+1j*ki)**2
    # print "Interpolated permittivity at "+str(wavelength*1e9)+" nm = "+str(epsilon)
    # except:
    # print "Interpolation for "+str(wavelength*1E9)+" nm failed."
    # try:
    # wavelength = 930e-9
    # ni = fni(wavelength); ki = fki(wavelength)
    # epsilon = (ni+1j*ki)**2
    # print "Interpolated permittivity at "+str(wavelength*1e9)+" nm = "+str(epsilon)
    # except:
    # print "Interpolation for "+str(wavelength*1E9)+" nm failed."

    # defining the new n and k on a common mesh
    ni = fni(wavelengths)
    ki = fki(wavelengths)

    plt.figure()
    plt.xlabel(r'$\mathcal{R}e(\varepsilon)$ (nm)')
    plt.ylabel('n, k')
    plt.semilogx(wavelength1 * 1e9, n, 'bs', label='n Palik')
    plt.semilogx(wavelength2 * 1e9, k, 'rs', label='k Palik')
    plt.semilogx(wavelengths * 1e9, ni, 'b-', label='n interp')
    plt.semilogx(wavelengths * 1e9, ki, 'r-', label='k interp')
    plt.grid()
    plt.legend(loc=1)
    plt.savefig('PalikData.eps')
    return epsilon
    # plt.show()


# }}}

## Imports optical data of type (wavelength, n,k) from different wavelength meshes. Return (n+k*1j)
#  Such data can be captured using a software like Engauge Digitizer. 
#  This leads to obtain (n,k) discretized on DIFFERENT MESHES. 
#  
#  NOTE: Warning: the data produced by this method are rather unprecise. SPP spectoscopy requires precision to 1E-3. 
#  This method gives a precision worst then 1E0. Use only in case no other data are available. 
def importFromNKtable_batch(folder, filename, plotting=1, unit=1E-6):  # {{{
    nfile = folder + filename + "-n.csv"
    kfile = folder + filename + "-k.csv"

    narray = np.loadtxt(nfile, delimiter="\t", skiprows=1)
    karray = np.loadtxt(kfile, delimiter="\t", skiprows=1)

    # unit = 1E-6

    # import wavelength, n and k from Palik
    wavelength1 = narray[:, 0] * unit
    wavelength2 = karray[:, 0] * unit
    n = narray[:, 1]
    k = karray[:, 1]
    numrows = 10000
    base = 10
    # interpolate n and k on new wavelength mesh
    order = 1
    # wavelengths = np.arange(np.amin(wavelength2),np.amax(wavelength2), precision) #regular mesh, AWFUL for memory
    print(("Generating new wavelength mesh: (" + str(np.amin(wavelength2)) + ", " + str(np.amax(wavelength2)) + ")"))
    wavelengths = np.logspace(np.amin(np.log10(wavelength2)), np.amax(np.log10(wavelength2)), num=numrows, base=base,
                              endpoint=True)

    print(("New wavelength mesh has " + str(numrows) + " rows."))
    # print wavelengths

    fni = InterpolatedUnivariateSpline(wavelength1, n, k=order)
    fki = InterpolatedUnivariateSpline(wavelength2, k, k=order)

    # Interpolated one optical constants
    # try:
    ##wavelength = 800e-9
    # ni = fni(wavelength); ki = fki(wavelength)
    # epsilon = (ni+1j*ki)**2
    # print "Interpolated permittivity at "+str(wavelength*1e9)+" nm = "+str(epsilon)
    # except:
    # print "Interpolation for "+str(wavelength*1E9)+" nm failed."

    # defining the new n and k on a common mesh
    ni = fni(wavelengths)
    ki = fki(wavelengths)

    if (plotting == 1):
        plt.figure()
        plt.xlabel(r'$\lambda$ (nm)')
        plt.ylabel('n, k')
        plt.semilogx(wavelength1 * 1e9, n, 'bs', label='n Palik')
        plt.semilogx(wavelength2 * 1e9, k, 'rs', label='k Palik')
        plt.semilogx(wavelengths * 1e9, ni, 'b-', label='n interp')
        plt.semilogx(wavelengths * 1e9, ki, 'r-', label='k interp')
        plt.grid()
        plt.legend(loc=1)
        plt.savefig('PalikData.eps')

    FinalNKarray = np.transpose([wavelengths, ni, ki])
    np.savetxt(folder + filename + ".txt", FinalNKarray, fmt="%s", delimiter="\t", header="wavelength (nm), \t n \t k",
               newline="\n", comments="# ")
    print("Merged optical data into a shared file...")

    return wavelengths, ni + 1.j * ki
    # plt.show()


# }}}

## Imports (wavelength, ReEps, ImEps) data captured using a software like Engauge Digitizer. 
#   Input: (epsReal, epsImag) are discretized, also works on DIFFERENT MESHES. 
#   Output: (epsilon) interpolated at a given value
#           and plots the whole spectrum for verification. 
#   
#   NOTE: Warning: the data produced by this method are rather unprecise. SPP spectoscopy requires precision to 1E-3. 
#   This method gives a precision worst then 1E0. Use only in case no other data are available. 
def importFromEpsilonTable(wavelength, folder, filename, plotting=True, unit=1E-10):  # {{{
    nfile = folder + filename + "-epsR.csv"  # TODO: rename nfile to ReEpsFile
    kfile = folder + filename + "-epsC.csv"  # TODO: rename kfile to ImEpsFile

    print(("** Info: opening " + nfile + " and " + kfile + "."))
    narray = np.loadtxt(nfile, delimiter="\t", skiprows=1)
    karray = np.loadtxt(kfile, delimiter="\t", skiprows=1)

    # unit = 1E-10 #Unit of the wavelength found in databases <nfile> and <kfile>.

    # import wavelength, epsilonR, epsilonC
    wavelength1 = narray[:, 0] * unit
    wavelength2 = karray[:, 0] * unit
    n = narray[:, 1]
    kk = karray[:, 1]
    numrows = 10000  # Number of rows to interpolate the data on.
    base = 10
    # interpolate on new wavelength mesh using 1st order
    order = 1
    # wavelengths = np.arange(np.amin(wavelength2),np.amax(wavelength2), precision) #regular mesh, AWFUL for memory
    print(("** Info: Generating the new wavelength mesh: (" + str(np.amin(wavelength2)) + ", " + str(
        np.amax(wavelength2)) + ")"))
    wavelengths = np.logspace(np.amin(np.log10(wavelength2)), np.amax(np.log10(wavelength2)), num=numrows, base=base,
                              endpoint=True)

    print(("New wavelength mesh has " + str(numrows) + " rows."))
    # print wavelengths

    print("** Info: definition of interpolation functions...")
    wavelength1 = np.sort(wavelength1)
    wavelength2 = np.sort(wavelength2)

    fni = InterpolatedUnivariateSpline(wavelength1, n, k=order)
    fki = InterpolatedUnivariateSpline(wavelength2, kk, k=order)

    # Interpolation for one optical constant
    try:
        ni = fni(wavelength);
        ki = fki(wavelength)
        # print wavelength, ni, ki
        epsilon = (ni + 1j * ki)  # NOTE: we are picking up the epsRe, and epsIm directly here
        print("")
        print(
            ("importFromEpsilonTable: Interpolated permittivity at " + str(wavelength * 1E9) + " nm = " + str(epsilon)))
    except:
        print(("importFromEpsilonTable: Interpolation for " + str(wavelength * 1E9) + " nm failed."))

    # defining the new epsR and epsC on a common mesh
    ni = fni(wavelengths)
    ki = fki(wavelengths)
    epsilons = (ni + 1j * ki)  ##!! Names are misleading here: we actually work dielectric permittivities!

    if (plotting):
        plt.figure()
        plt.xlabel(r'$\mathcal{R}e(\varepsilon)$ (nm)')
        plt.ylabel(r'Re$(\varepsilon)$, Im($\varepsilon$)')
        plt.semilogx(wavelength1 * 1e9, n, 'bs', label=r'Re$(\varepsilon)$ data')
        plt.semilogx(wavelength2 * 1e9, kk, 'rs', label=r'Im$(\varepsilon)$ data')
        plt.semilogx(wavelengths * 1e9, ni, 'b-', label=r'Re$(\varepsilon)$ interp')
        plt.semilogx(wavelengths * 1e9, ki, 'r-', label=r'Im$(\varepsilon)$ interp')
        plt.grid()
        plt.legend(loc=1)
        plt.savefig('GraphData.eps')
        plt.show()

    return epsilon


# }}}

## Imports (wavelength, ReEps, ImEps) data captured using a software like Engauge Digitizer. 
#   Input: file with (wavelengths, epsReal, epsImag) are discretized, also works on DIFFERENT MESHES. 
#   Output: epsilons interpolated on a whole mesh
#   
#   NOTE: Warning: the data produced by this method are rather unprecise. SPP spectoscopy requires precision to 1E-3. 
#   This method gives a precision worst then 1E0. Use only in case no other data are available. 
def importFromEpsilonTable_batch(folder, filename, plotting=True, unit=1E-10):  # {{{
    nfile = folder + filename + "-epsR.csv"  # TODO: rename nfile to ReEpsFile
    kfile = folder + filename + "-epsC.csv"  # TODO: rename kfile to ImEpsFile

    print(("** Info: opening " + nfile + " and " + kfile + "."))
    narray = np.loadtxt(nfile, delimiter="\t", skiprows=1)
    karray = np.loadtxt(kfile, delimiter="\t", skiprows=1)

    # unit = 1E-10 #Unit of the wavelength found in databases <nfile> and <kfile>.

    # import wavelength, epsilonR, epsilonC
    wavelength1 = narray[:, 0] * unit
    wavelength2 = karray[:, 0] * unit
    n = narray[:, 1]
    kk = karray[:, 1]
    numrows = 10000  # Number of rows to interpolate the data on.
    base = 10
    # interpolate on new wavelength mesh using 1st order
    order = 1
    # wavelengths = np.arange(np.amin(wavelength2),np.amax(wavelength2), precision) #regular mesh, AWFUL for memory
    print(("** Info: Generating the new wavelength mesh: (" + str(np.amin(wavelength2)) + ", " + str(
        np.amax(wavelength2)) + ")"))
    wavelengths = np.logspace(np.amin(np.log10(wavelength2)), np.amax(np.log10(wavelength2)), num=numrows, base=base,
                              endpoint=True)

    print(("New wavelength mesh has " + str(numrows) + " rows."))
    # print wavelengths

    print("** Info: definition of interpolation functions...")
    wavelength1 = np.sort(wavelength1)
    wavelength2 = np.sort(wavelength2)

    fni = InterpolatedUnivariateSpline(wavelength1, n, k=order)
    fki = InterpolatedUnivariateSpline(wavelength2, kk, k=order)

    ##Interpolation for one optical constant
    # try:
    # ni = fni(wavelength); ki = fki(wavelength)
    ##print wavelength, ni, ki
    # epsilon = (ni+1j*ki) #NOTE: we are picking up the epsRe, and epsIm directly here
    # print ""
    # print "Interpolated permittivity at "+str(wavelength*1E9)+" nm = "+str(epsilon)
    # except:
    # print "Interpolation for "+str(wavelength*1E9)+" nm failed."

    # defining the new epsR and epsC on a common mesh
    ni = fni(wavelengths)
    ki = fki(wavelengths)
    epsilons = (ni + 1.j * ki)  ##WARNING: Names are misleading here: we actually work dielectric permittivities!

    if (plotting):
        plt.figure()
        plt.xlabel(r'$\mathcal{R}e(\varepsilon)$ (nm)')
        plt.ylabel(r'$Re(\varepsilon)$, $Im(\varepsilon)$')
        plt.semilogx(wavelength1 * 1e9, n, 'bs', label=r'Re$(\varepsilon)$ data')
        plt.semilogx(wavelength2 * 1e9, kk, 'rs', label=r'Im$(\varepsilon)$ data')
        plt.semilogx(wavelengths * 1e9, ni, 'b-', label=r'Re$(\varepsilon)$ interp')
        plt.semilogx(wavelengths * 1e9, ki, 'r-', label=r'Im$(\varepsilon)$ interp')
        plt.grid()
        plt.legend(loc=1)
        plt.savefig('GraphData.eps')
        # plt.show()

    return wavelengths, epsilons


# }}}

## Mere function? Generates (ReEps, ImEps) from absorption data (given in m^{-1}). HOW?
#  This routine is made to import data captured using sofware such as Engauge Digitized. 
#  Be very careful! The data produced by this method are very unprecise. SPP spectroscopy requires precision to 1E-3. 
#  This method gives a precision worst then 1E0. Then, it is only in case we have no other data.
# @param wavelength: range to express data onto
# @param folder: location of the data
# @param filename: ".CSV" will be added to the name
# @param plotting: True|False
# @param unit: unit of the wavelength (nm by default)
def importFromAbsorptionData(folder, filename, energyAxis=False, plotting=True, unit1=1E9, unit2=1E-2):  # {{{
    absfile = folder + filename + ".csv"
    narray = np.loadtxt(absfile, delimiter="\t", skiprows=1)
    #unit = 1E9

    # import wavelength, alpha from spectroscopic data
    if(energyAxis):
        wavelength1 = h * c / (narray[::-1, 0] * e)
    else:
        wavelength1 = narray[::, 0] / unit1
    print(wavelength1)
    absorptivity = narray[::, 1] / unit2
    print(absorptivity)
    numrows = 10000
    base = 10
    # interpolate absorptivity k on new wavelength mesh
    order = 1
    # wavelengths = np.arange(np.amin(wavelength2),np.amax(wavelength2), precision) #regular mesh, AWFUL for memory
    print(("Generating new wavelength mesh: (" + str(np.amin(wavelength1)) + ", " + str(np.amax(wavelength1)) + ")"))
    wavelengths = np.logspace(np.amin(np.log10(wavelength1)), np.amax(np.log10(wavelength1)), num=numrows, base=base,
                              endpoint=True)

    print(("New wavelength mesh has " + str(numrows) + " rows."))
    # print wavelengths

    fni = InterpolatedUnivariateSpline(wavelength1, absorptivity, k=order)

    # Interpolated one optical constants
    try:
        absnew = fni(wavelengths)
        print(("Interpolated absorptivity at " + str(wavelengths * unit1) + " nm = " + str(absnew)))
    except:
        print(("Interpolation for " + str(wavelengths * unit1) + " nm failed."))

    # defining k on a common mesh
    absorptivity_new = fni(wavelengths)

    if (plotting):
        plt.figure()
        plt.xlabel(r'Wavelength $\lambda$ (nm)')
        plt.ylabel(r'$\alpha$ (m$^{-1}$)')
        plt.loglog(wavelength1 * unit1, absorptivity, 'bs', label=r'$\alpha$ data')
        plt.loglog(wavelengths * unit1, absorptivity_new, 'b-', label=r'$\alpha$ interp')
        plt.grid()
        plt.legend(loc=1)
        plt.show()
        plt.savefig('Absorptivity.eps')

    return wavelengths, absorptivity_new * wavelengths / 4. / np.pi
    # plt.show()


# }}}

## Simply plots the optical data taken from a compatible database, 
# and interpolate palik data from tables of Palik at the given wavelength. 
# Useful to add one set of (ReEps, ImEps) for ONE wavelength in MaterialOpticalDatabaseForPlasmonics.csv. 
# INPUT
# @param wavelength: (float) value of desired output wavelength
# @param folder:     (str) name of the folder were database can be found
# @param filename:   (str) name of the material file
# @param plotting:   (boolean) plot the full data if True
def importFromTable(wavelength, folder, filename, plotting):  # {{{
    # Look for Palik into the name
    if (filename.find("Palik") > 0):
        print("Palik data identified.")
        unit1 = 1E-10  # Palik data
    elif (filename.find("Gori") > 0):
        print("Gori data identified.")
        unit1 = 1E-9
    else:
        print("** Warning: Default case: choosing um for the input.")
        unit1 = 1E-6  # Other data

    # Fetch data
    DataFile = folder + filename
    try:
        DataArray = np.loadtxt(DataFile, delimiter="\t", skiprows=4)
        wavelengths = DataArray[:, 0] * unit1;
    except:
        DataArray = np.loadtxt(DataFile, delimiter=" ", skiprows=4)
        wavelengths = DataArray[:, 0] * unit1

    n = DataArray[:, 1];
    kk = DataArray[:, 2];

    print(n)

    # Interpolating using splines
    order = 1
    fni = InterpolatedUnivariateSpline(wavelengths, n, k=order)
    fki = InterpolatedUnivariateSpline(wavelengths, kk, k=order)

    # Interpolated one optical constants
    # wavelength = 1030e-9
    ni = fni(wavelength);
    ki = fki(wavelength)
    print((wavelength, ni, ki))
    epsilon = (ni + 1j * ki) ** 2
    print(("importFromTable: Interpolated permittivity at " + str(wavelength * 1e9) + " nm = " + str(epsilon)))

    # Interpolate the full array and check it visually
    nimesh = fni(wavelengths);
    kimesh = fki(wavelengths)

    if (plotting):
        plt.figure()
        plt.xlabel(r'$\lambda$ (nm)')
        plt.ylabel('n, k')
        plt.semilogx(1e9 * wavelengths, n, 'bs', label='n Palik')
        plt.semilogx(1e9 * wavelengths, kk, 'rs', label='k Palik')
        plt.semilogx(1e9 * wavelengths, nimesh, 'b-', label='n interp')
        plt.semilogx(1e9 * wavelengths, kimesh, 'r-', label='k interp')
        plt.grid()
        plt.legend(loc=2)
        plt.savefig('PalikData.eps')
        plt.show()

    return epsilon


# }}}

## Imports the optical data from Palik database
# Returns the spectrum (in Re(eps), Im(eps)) interpolate of a material interpolated on the given grid of wavelengths
# @param wavelengths: vector of wavelengths
# @param DataFile="Si-Palik": filename from Palik book, where data files are in Angstroms, n, k form.
def ImportPalikDatabase_epsilon_fromNK(wavelengths, DataFile="Si-Palik"):  # {{{
    import os
    from scipy.interpolate import InterpolatedUnivariateSpline
    PathFile=os.environ["spp_extended_theory"]+"/spp_extended_theory/SimpleSPP/Database/"
    DataFile=PathFile+DataFile
    DataArray = np.loadtxt(DataFile, delimiter="\t", skiprows=4)
    wavelengths_source = DataArray[:, 0]
    if (np.min(np.diff(wavelengths_source)) == 0e0):
        logger.warning(
            Header + "Error in wavelengths_source. Careful! Palik book combines different optical data, hence there can be overlap, sometimes! It must be removed to be imported properly here. ")
        # print np.diff(wavelengths_source)
        sys.exit()
    OpticalN = DataArray[:, 1]
    OpticalK = DataArray[:, 2]
    order = 1
    fni = InterpolatedUnivariateSpline(wavelengths_source * 1E-10, OpticalN, k=order)
    fki = InterpolatedUnivariateSpline(wavelengths_source * 1E-10, OpticalK, k=order)
    OpticalN_new = fni(wavelengths)
    OpticalK_new = fki(wavelengths)
    epsilon = np.power(OpticalN_new + np.multiply(1j, OpticalK_new), 2)
    # epsR = np.real(epsilon); epsC = np.imag(epsilon)
    return epsilon


# }}}

## Combines 2 sets of optical data using harmonic interpolation. 
# @param wavelength1: mesh set of the 1st optical data
# @param wavelength2: mesh set of the 2nd optical data
# @param eps1: complex scalar field of 1st optical data
# @param eps2: complex scalar field of 2nd optical data
# @param unit: base unit for the file
def interpolateTwoSetsOfOpticalData(wavelength1, wavelength2, eps1, eps2):
    # import wavelength, n and k from Palik
    numrows = 10000
    base = 10
    # interpolate n and k on new wavelength mesh
    order = 1
    # wavelengths = np.arange(np.amin(wavelength2),np.amax(wavelength2), precision) #regular mesh, AWFUL for memory
    wavelength_min = np.amin([wavelength1.min(), wavelength2.min()])
    wavelength_max = np.amax([wavelength1.max(), wavelength2.max()])
    print(("Generating new wavelength mesh: (" + str(np.amin(wavelength_min)) + ", " + str(
        np.amax(wavelength_max)) + ")"))
    wavelengths = np.logspace(np.amin(np.log10(wavelength_min)), np.amax(np.log10(wavelength_max)), num=numrows,
                              base=base, endpoint=True)

    print(("New wavelength mesh has " + str(numrows) + " rows."))
    # print wavelengths
    eps1r = eps1.real;
    eps1c = eps1.imag
    eps2r = eps2.real;
    eps2c = eps2.imag
    feps1r = InterpolatedUnivariateSpline(wavelength1, eps1r, k=order, ext=1)
    feps1c = InterpolatedUnivariateSpline(wavelength1, eps1c, k=order, ext=1)
    feps2r = InterpolatedUnivariateSpline(wavelength2, eps2r, k=order, ext=1)
    feps2c = InterpolatedUnivariateSpline(wavelength2, eps2c, k=order, ext=1)
    print((Header + "Interpolation functions are ready."))

    ## METHOD ?
    ## Harmonic average? Does not work with bounded-energy permittivities
    ## Algebraic average? Have no physical meaning.
    ## Adding them? If they are with separate support, then yes.
    ## Maxwell Garnett? With which fraction then?

    print("Combining sets of optical data via ADDING them [!they must have different support!]...")
    # print feps2r(wavelengths)+1j*feps2c(wavelengths)
    # epsilon_final = feps1r(wavelengths)+1.j*feps1c(wavelengths) + feps2r(wavelengths)+1.j*feps2c(wavelengths) #they have different support, hence it should be fine]
    epsilon_final = feps1r(wavelengths)  # + feps2r(wavelengths)
    # print np.shape(wavelengths), np.shape(epsilon_final)
    # return wavelengths, epsilonR_final+1j*epsilonC_final
    return wavelengths, epsilon_final

# ==============================
