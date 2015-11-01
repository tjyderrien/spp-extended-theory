# README #

# Concept

This set of Python functions allow to predict possibilities of SPP on various materials, and for various wavelengths.

The model is limited to two semi-infinite materials sharing one interface. 

# List of executable files 

## plotSPPmultiMaterialData.py: 

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
