[![Homepage](https://img.shields.io/badge/Home-quantumlap.eu-green.svg)](http://www.quantumlap.eu)
[![License: GPL v3](https://img.shields.io/badge/License-GPL%20v3-blue.svg)](http://www.gnu.org/licenses/gpl-3.0)
![Bitbucket open issues](https://img.shields.io/bitbucket/issues/tjyderrien/spp_extended_theory?style=plastic)
![Twitter Follow](https://img.shields.io/twitter/follow/tjyderrien?style=social)
![Website](https://img.shields.io/website?down_color=red&down_message=quantumlap.eu%20%5Btemporarily%20down%21%5D&up_color=green&up_message=quantumlap.eu&url=http%3A%2F%2Fwww.quantumlap.eu)
![Maintenance](https://img.shields.io/maintenance/yes/2021)


# README #

# Structure of the code

The code is now made of 3 branches. 
* Branch release: contains the seemingly stable public code for the community.
* Branch master:  contains a more advanced version of code for the community. 
* Branch develop: contains the latest code, with non-tested features and possible bugs. 

# Main developers

* T.J.-Y. Derrien: developing SPP and keldysh theory, programming. 

* F. Preucil: programming the multilayer solver. Started improvement of multilayer SPP model. 

* K. Hlinomaz: a bit of Python maintenance at the moment

* J. Bonse: Provided the tabulated Palik data.

# Users of the code

* Inam Mirza (Hilase Centre, Prague): user for compound materials

* Aleksander Kovacevic (Bratislava, Serbia): multilayered materials

* Stephan Gräf (Jena University): user for Sipe theory.

# Financial supports

* Adolf Martens Fellowship: BAM Federal Institute for Materials Research and Testing. 

* HiLASE projects (State Budget of Czech Republic). 

* Marie Sklodowska-Curie Actions: funding from the European Commission for the Marie Sklodowska-Curie Individual Fellowship under QuantumLaP project No. 657424. 

# Purpose of this project

You may find a collection of programs aiming to develop control over the transient excitation of Surface Plasmon Polaritons (SPPs) upon laser irradiation. 
The project has diverged into a number of programs. Refactoring is planned. 
SPPs can be used with light to enhance local power density, but also to induce Laser-Induced Periodic Surface Structures (LIPSS). 

## Concept

This code contains optical data for many materials on a large range of wavelengths. These databasis were tabulated from Palik (Handbook of optical constants, Academic Press 1985), and sometimes refractiveindex.info databases. There exists possibility to import optical data using importPalikData.py. An automatic documentation can be found in the folder ./doc. If not up to date, run <bash generate.sh> (Linux). 

## Contents

For each (material, wavelength), the following theories can be applied. 

* SPP theory for single interface (coded in Python 2, validated, multi-material, published), 

* SPP theory for multilayer systems: 3 layer model (available with Maple: under developement in Python 3), published.

* SPP theory for rough interfaces: Sipe model (coded in Python 2, validated), Sipe-Drude model (coded in Python 2, validated, published by Jörn Bonse in 2005, 2009).

# Installation 

## Simple install for beginners and developpers

pip install -e git+https://bitbucket.org/tjyderrien/spp-extended-theory/

or

1. git clone https://bitbucket.org/tjyderrien/spp-extended-theory/
2. git checkout Modules
3. conda env create -f environment.yaml #if you want to install a separate conda environment for this code

If this fails, look for your install of Conda ("whereis conda" or "locate conda"). Usually it is practical to introduce this into .bashrc: export PATH=/opt/miniconda3/bin:$PATH

4. conda activate spp-extended-theory #this name should be present after using Step 3.
5. pip install -e <folder_with_setup.py> (can be "." or "spp-extended-theory")

## Simple uninstall for beginners

pip uninstall spp-extended-theory

## Note to Windows users

The code works with Windows: 
* Install Anaconda (Python) for Windows (during install: no need for the symbolic linking compatibility).
* Download the present repository. 

Most of routines were not tested with Windows, so feel free to address any problem to <derrien@fzu.cz>.


# How to contribute ?

* Write message about the encountered problems in the forum

* Develop / correct a part of the program. Please document your contributions using dOyxgen style (see https://www.stack.nl/~dimitri/doxygen/manual/docblocks.html#pythonblocks).

* Cite this work if it was useful in your publications.  

## Use git to publish your changes. Here are basic commands

* To download the simulation code, type: git clone git@bitbucket.org:tjyderrien/spp_extended_theory.git

* git pull: download latest changes to your version of the code (no worries, it is non-destructive to your changes!)

* git add [File1 [File2 [...]]]: use this command to add the files you would like to update/add. 

* git commit -m 'LIPSS database: update of Aug 7th 2015': Just give the description of changes you have applied. 

* git push: Upload your changes to the central server on Bitbucket.org.
