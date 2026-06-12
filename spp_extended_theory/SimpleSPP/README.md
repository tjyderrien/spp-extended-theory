# README

This file details the spp_extended_theory code. 

# Concept

This set of Python functions is aiming at performing calculations to predict spectroscopical features of laser-matter interaction. In particular, emphasis is put on prediction of Surface Plasmon Polaritons on various materials, and for various wavelengths. 

# Note

SimpleSPP folder is limited to two semi-infinite materials sharing one interface. 

# Description of executable files

## plotMultimaterials.py: 

Plots the SPP period as function of material's Re(epsilon). 

This file uses self-made database Material.csv, where different data were collected by hand. 

## plotMultiwavelength.py: 

Generates several plots as function of irradiation wavelength: 

* Dielectric permittivity as function of wavelength

* Dispersion relation w(k)

* SPP period as a function of wavelength

* SPP life time as function of wavelength, and compared with Raether formula. 

This file uses database with the shape (lambda, n, k) that can be found on some websites: 

* refractiveindex.info

* SOPRA database (given with GsVit software)

/!\ Refractive index databases can completely modify results. In particular, SOPRA database is very bad for SPP group velocities (> c !!). 

## makeTable.py

* Generate a table of SPP-active materials. Modify the source to get list of active interfaces for:  symmetric SPP, asymmetric +/- SPP 


# Libraries: 

## libSPP.py: 

This file contains many basic routines supporting the SPP theory: 

### Basics for lasers

* omega(wavelength): just gives laser frequency from laser wavelength

* reflectivity(eps1, eps2): calculate the normal incidence reflectivity

* EpsilonToIndex(eps): generate n and k value for one complex epsilon value

### Basic routines of the SPP model

* betaSPP(wavelength, eps1, eps2): calculate SPP wavenumber from dielectric permittivities

* AsymmetricSPPconditionPos(eps1, eps2): + mode of asymmetric SPP condition

* AsymmetricSPPconditionNeg(eps1, eps2): - mode of asymmetric SPP condition

* SPPconditionValue(eps1, eps2): if <0, then symmetric SPP can be excited

* SPPcondition(eps1, eps2): tests if SPPconditionValue < 0 for one material combination

* OldSPPcondition(eps1, eps2): tests the SPP condition of non-absorbing materials (this theory is still widely used though not sufficient for most materials)

### Properties of SPPs

* period(beta): returns the period out of a (complex-valued) SPP wave number

* DecayDepth(kzSPP): Calculate decay depth with kzSPP (complex)

* kzSPP(wavelength, eps1, eps2): calculate wavenumber in direction of propagation

* DecayLengthSPP(beta): gives the decay distance of SPPs on surface

* SPPlength(beta): SAME FUNCTION. FIX THIS. 

* LifeTimeRaether(beta, eps2, eps1): provides SPP lifetime using Raether formula

* LifetimeDerrien(beta, vg): lifetime based on group velocity

* EffectiveIndex(eps1, eps2): generates effective optical index of an SPP (complex-valued quantity)

### SPP database generation

* Drude(wavelength, ne, epsilon, nu): provides excitation of carriers - this function is incomplete for now

* ExperimentallyAchievable(OPD, DecayDepth): some function which checks if SPP decay depth is larger than optical penetration depth - (high risk experimental criterion)

* SPPactiveInterfaces(dbarray, comment): this produces a list of SPP active interfaces based on a database

* AsymmetricSPPposActiveInterfaces(dbarray, comment): same but for Asymmetric (+) SPPs

* AsymmetricSPPnegActiveInterfaces: same but for Asymmetric (-) SPPs. 

* GenerateDatabase(): build SPP database array from MaterialOpticalDatabaseForPlasmonics.csv

* ExportToTxt(dbarray, filename): exports an array to a textfile

### Mathematical functions

* RealDerivativeByComplex(f,z): complex derivative 
