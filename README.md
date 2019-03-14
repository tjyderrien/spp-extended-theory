# README #

# Null modification

# Structure of the code

The code is now made of 3 branches. 
* Branch release: contains the seemingly stable public code for the community.
* Branch master:  contains a more advanced version of code for the community. 
* Branch develop: contains the latest code, with non-tested features and possible bugs. 

# Main developers

* T.J.-Y. Derrien: developing theory, programming. 

* Y. Levy: LIPSS database contributor.

* J. Bonse: Provided the tabulated Palik data.

# Current collaborators for this code

* Stephan Gräf (Jena University): user for Sipe theory.

* Inam Mirza (Hilase Centre, Prague): user for compound materials

* Aleksander Kovacevic (Bratislava, Serbia): multilayered materials

# Financial supports

* Adolf Martens Fellowship: BAM Federal Institute for Materials Research and Testing. 

* HiLASE projects (State Budget of Czech Republic). 

* Marie Sklodowska-Curie Actions: funding from the European Commission for the Marie Sklodowska-Curie Individual Fellowship under QuantumLaP project No. 657424. 

# Purpose of this project

You may find a collection of programs aiming to develop control over the transient excitation of Surface Plasmon Polaritons (SPPs) upon laser irradiation. 
SPPs can be used with light to enhance local power density, but also to induce Laser-Induced Periodic Surface Structures (LIPSS). 

## Concept

This code contains optical data for many materials on a large range of wavelengths. These databasis were tabulated from Palik (Handbook of optical constants, Academic Press 1985), and sometimes refractiveindex.info databases. There exists possibility to import optical data using importPalikData.py. 

## Contents

For each of these data, the following theories can be applied. 

* SPP theory for single interface (coded, validated, multi-material, published), 

* SPP theory for multilayer systems: 3 layer model (available with Maple: to be coded in Python), published.

* SPP theory for rough interfaces: Sipe model (coded, validated), Sipe-Drude model (coded in Python, validated, published by Jörn Bonse in 2005, 2009).

# Installation 

## Simple install for beginners

* Simply click "download" in the repository bar. You can also click here: https://bitbucket.org/tjyderrien/spp-extended-theory/downloads. 

* To get the regular updates without loosing your modifications, it is recommended to learn basics of Git. 

## For users familiar with Git. 

* For educated users: it is advised to download in HTTP mode using git clone https://tjyderrien@bitbucket.org/tjyderrien/spp-extended-theory.git

* For hardcore users, you may need SSH protocol (key is required, contact me via email) to donwload: 
 git clone git@bitbucket.org:tjyderrien/spp-extended-theory.git

## Note to Windows users

The code works with Windows: 
* Install Anaconda (Python) for Windows (during install: no need for the symbolic linking compatibility). 
* Download the present repository,  
* Open/Edit the .py files and enjoy. 

Most of routines were not tested with Windows, so feel free to address any problem to <derrien@fzu.cz.>

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
