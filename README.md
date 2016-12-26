# README #

# Structure of the code

The code is now made of 2 branches. 
* Branch master: contains the public code for the community. 
* Branch develop: contains the private code, with unreleased features. 

# Main developers

* T.J.-Y. Derrien: developing theory, programming. 

* Y. Levy: LIPSS database contributor.

* J. Bonse: Tabulated Palik data.

# Current collaborators

* Inam Mirza, Hilase/FZU, Prague, Czech Republic

* Iaroslav Gnilitskyi, UNIMORE, University of Modena, Italy

# Financial supports

* Adolf Martens Fellowship: BAM Federal Institute for Materials Research and Testing. 

* HiLASE projects (State Budget of Czech Republic). 

* Marie Sklodowska-Curie Actions: project QuantumLaP. 

# Purpose of this project

You may find a collection of programs aiming to develop control over the transient excitation of Surface Plasmon Polaritons (SPPs) upon laser irradiation. 
SPPs can be used with light to enhance local power density, but also to induce Laser-Induced Periodic Surface Structures (LIPSS). 

## Concept

This code contains optical data for many materials on a large range of wavelengths. These databasis were tabulated from Palik (Handbook of optical constants, Academic Press 1985), and sometimes refractiveindex.info databases. There exists possibility to import optical data using importPalikData.py. 

## Contents

For each of these data, the following theories can be applied. 

* SPP theory for single interface (coded, validated, multi-material, published), 

* SPP theory for multilayer systems: 3 layer model (available with Maple: to be coded in Python), 

* SPP theory for rough interfaces: Sipe model (coded, validated), Sipe-Drude model (coded in Python).

# Install

* For beginners: To simply download programs, you can click here: https://bitbucket.org/tjyderrien/spp-extended-theory/downloads

* For more advanced users: git clone git@bitbucket.org:tjyderrien/spp-extended-theory.git

# How to use ?

* To use the code, install Anaconda on your system (https://www.continuum.io/downloads), and then execute: "git clone git@bitbucket.org:tjyderrien/spp-extended-theory.git" in a terminal. 

# How to contribute ?

* Write message in the forum

* Develop / correct a part of the program. Please document your contributions using dOyxgen style (see https://www.stack.nl/~dimitri/doxygen/manual/docblocks.html#pythonblocks).

* Point out some mistakes. Suggest corrections. 

* Cite this work if it was useful in your publications.  

## Use git to publish your changes. Here are basic commands

* To download the simulation code, type: git clone git@bitbucket.org:tjyderrien/spp-extended-theory.git

* git pull: download latest changes to your version of the code (no worries, it is non-destructive to your changes!)

* git add [File1 [File2 [...]]]: use this command to add the files you would like to update/add. 

* git commit -m 'LIPSS database: update of Aug 7th 2015': Just give the description of changes you have applied. 

* git push: Upload your changes to the central server on Bitbucket.org.
