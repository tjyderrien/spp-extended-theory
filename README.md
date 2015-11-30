# README #

# Main contributors

* T.J.-Y. Derrien : scientist, Hilase Centre FZU, AV CR

* J. Bonse, BAM Berlin: scientist, BAM Federal Institute of Materials Research and Testing, Berlin

# Technical support

* R. Garai, Hilase Centre FZU, AV CR : helped for some LIPSS databasing

# Current users

* Yoann Levy, FZU

* Inam Mirza, FZU

* Iaroslav Gnilitskyi

# Financial supports

* Adolf Martens Fellowship: BAM Federal Institute for Materials Research and Testing. 

* HiLASE project. 

* Marie Sklodowska-Curie Actions: project QuantumLaP. 

# Purpose of this project

You may find a collection of programs aiming to develop control over the excitation of Surface Plasmon Polaritons (SPPs). SPPs can be used with light to enhance local power density, but also to induce Laser-Induced Periodic Surface Structures (LIPSS). 

## Concept

This code contains optical data for many materials on large range of wavelength. These databasis were tabulated from Palik (Handbook of optical constants, Academic Press 1985), and refractiveindex.info databasis. 

## Contents

For each of these data, the following theories can be applied. 

* SPP theory for single interface (coded, validated, multi-material, ready for publication), 

* SPP theory for multilayer systems: 3 layer model (to be coded), 

* SPP theory for rough interfaces: Sipe model (coded, validated), Sipe-Drude model (coded).

# How to use ?

* To simply download programs, you can click here: https://bitbucket.org/tjyderrien/spp-extended-theory/downloads

* To use the code, install Anaconda on your system (https://www.continuum.io/downloads), and then execute: "git clone git@bitbucket.org:tjyderrien/spp-extended-theory.git" in a terminal. 

# How to contribute ?

* Write message in the forum

* Develop / correct a part of the program

* Point out some mistakes. Suggest corrections. 

* Cite this work if it was useful in your publications.  

## Use git to publish your changes. Here are basic commands

* To download the simulation code, type: git clone git@bitbucket.org:tjyderrien/spp-extended-theory.git

* git pull: download latest changes to your version of the code (no worries, it is non-destructive to your changes!)

* git add [File1 [File2 [...]]]: use this command to add the files you would like to update/add. 

* git commit -m 'LIPSS database: update of Aug 7th 2015': Just give the description of changes you have applied. 

* git push: Upload your changes to the central server on Bitbucket.org. 
