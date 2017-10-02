#!/usr/bin/env python
#-*- coding: utf-8 -*-

# Copyright (C) 2013-2017 T. J.-Y. Derrien
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

## @package plotSipe
# Module plotSipe provides the coupling efficiency factor as function of wavenumber kappa in 2D. 
# The system is a semi-infinite medium described by an homogeneous complex-valued dielectric permittivity. 
# The model was strictly taken from [Bonse, J. et al, J. Appl. Phys. 97, 013538 (2005)], which is clever summary of 
# the Sipe model given in [Sipe, J. E. et al. Phys. Rev. B 27, 1141-1154 (1983)]
# NOTE: module was only validated for normal incidence. 

#TODO: To verify coupling efficiency factor from Sipe theory: plot the efficiency factor maximum as function of the laser wavelength, and correlate with papers such as Endriz and Spicer, PRB 4, 4144 (1971); Benneth and Porteus, JOSA 51, 123 (1961)

## IMPORT LIBRARIES
# IMPORT PYTHON LIBRARIES


import numpy as np
from numpy import genfromtxt, loadtxt, chararray
from scipy.optimize import fsolve, root
import cmath
import matplotlib as mp
import matplotlib.pyplot as plt
from scipy.interpolate import InterpolatedUnivariateSpline
from matplotlib import rc, font_manager
# from pylab import *
from scipy.constants import c, epsilon_0, mu_0, pi, e, m_e, h
from matplotlib.legend_handler import HandlerLine2D
import sys

# IMPORT CUSTOM LIBRARIES

# from libKeldysh import *
from libDatabase import *
from libLaser import *
from libMaterials import *
from libMath import *
from libSPP import *
from libSipe import *
#from plotGraph import *

Header="[plotSipe.py] "



#==== Attempting a 1D plot
# Known quantities
theta = 0e0 #Single value here, but we can vectorize functions easily later. 
f = 0.1e0 #Filling factor: taken from Bonse et al, JAP (2009)
s = 0.4e0 #Shape factor: taken from Bonse et al, JAP (2009)
#wavelength = 800e-9
#epsilon = 12.80259+0.00109j
#epsilon = -97.593456+25.2698472743j

# Meshes for solution
#kappax = np.arange(0, 4, 0.1)
#kappay = np.arange(0, 4, 0.1)
#kappa = np.array([wavelength * 1, wavelength * 0]); #test values
kappax = 4e0 ; kappay = 0e0*kappax; #test values
numberlevels = 8 #for the final 2D plot
#for kappax in meshkappa:

#ftab = np.arange(0, 1, 0.1)
print "Info: Generating the mesh..."
kapparange = np.arange(0.1,4,0.1)
#for wavelength in wavelengths
#for f in ftab:
#idtab = 0
#etaresult = np.zeros(kapparange.shape)
#kappax = 0e0
#for kappay in kapparange:
  
  ### Defining simple quantities for Sipe model
  #kappa = np.array([kappax, kappay])
  #kappai = np.array([-cmath.sin(theta), 0])
  #kappap = kappai + kappa; kappam = kappai - kappa

  ##print "kappax = "+str(kappax)
  ##print idtab
  #etaresult[idtab] = etas(theta, f, s, epsilon, kappa, kappap, kappam)
  ##print etaresult[idtab]
  ##print "eta = "+str(etaresult)
  #idtab = idtab+1

##=========== Make a 1D plot

#plt.figure()
#plt.xlabel(r'$\kappa_x$')
#plt.ylabel(r'$\eta$')
#print kapparange.shape, etaresult.shape
#plt.plot(kapparange, etaresult, '-', label='Sipe')
#print etaresult
#plt.savefig('SipeEtaKappaX.eps')
##exit()

#=========== 2D plot

query = 'Air'
#query2= 'InP (Bonse 2005)'
#query2= 'Mo (Ordal 1988)'
#query2= 'Cu (Palik)'
query2='SiO2 (Palik)'
print "Caution: the expression must be exactly the one of MaterialDatabase.csv."

wavelength = 1025
select = str(wavelength)
unit = 1E-9
wavelength = wavelength * unit
print "Wavelength = "+str(wavelength/unit)+" nm."

## Prepares the usual SipeEfficiencyFactor(kx, ky) for a specific material query2 immersed in Air. 
# @param wavelength: photon energy given in SI (meters)
# @param query2: <string> linking to a material given in ../MaterialDatabase.csv.
# NOTE: Sipe model is limited to air-material interface. Cannot be used with water-material for example. 
#       For more advanced combinations of materials, see [T.J.-Y. Derrien et al, Journal of Optics 18, 115007 (2016)]
def plotSipeFromDatabase(wavelength, query2, query='Air'): #{{{
  ## Generate the database
  SPPdb = GenerateDatabase() #Generate from MaterialDatabase.csv
  print "SPP database has "+str(len(SPPdb))+" entries."

  #print "Full Database:"
  #print SPPdb

  # Select the material of interface 1
  SPPdb = FilterDatabase(SPPdb, query, 0)
  print "Filter on materials: SPP database has now "+str(len(SPPdb))+" entries."

  #print SPPdb #works well

  # Filter database on wavelength
  try: 
	  title = select+' nm'
	  SPPdb = FilterDatabase(SPPdb, select+".0", 2)
	  print "Filter on wavelength: SPP database "+title+" has "+str(len(SPPdb))+" entries."
  except:
	  print "Exception: no optical data is available for "+query+" at "+title+"."
	  print SPPdb
	  exit()
    


  # Filter database on materials

  try: 
	  title = query2
	  SPPdb = FilterDatabase(SPPdb, query2, 1)
	  print "Filter on material: SPP database "+title+" has "+str(len(SPPdb))+" entries."
	  #print SPPdb
  except:
	  print "Exception: no optical data is available for "+query+" at "+title+"."
	  print SPPdb
	  exit()

  if(len(SPPdb)==0):
    print "SPP database returned 0 matching result."
    exit()
    
  # Extract materials from database
  Material1, Material2, Wavelength, OldSPPactiveBool, NewSPPactiveBool, SPPperiod, SPPperiodError, SPPdecayDepth1, SPPdecayDepth2, Reflectivity, OpticalPenetration1, OpticalPenetration2, SPPdecayLength, eps1rM, eps1cM, eps2rM, eps2cM, k1imag, k2imag, DeltaLsppValue = ExtractDataDb(SPPdb)

  # Calculation of refractive index
  eps1rM=np.asfarray(eps1rM)
  eps1cM=np.asfarray(eps1cM)
  eps2rM=np.asfarray(eps2rM)
  eps2cM=np.asfarray(eps2cM)

  epsilon1 = np.add(eps1rM,np.multiply(1e0j, eps1cM))
  epsilon2 = np.add(eps2rM,np.multiply(1e0j, eps2cM))
    
  print "Mesh generation..."
  k_precision = 5e-2
  SipeRanges = 2e0
  kx = np.arange(-SipeRanges,SipeRanges,k_precision)
  ky = np.arange(-SipeRanges,SipeRanges,k_precision)
  kxx, kyy = np.meshgrid(ky, kx)

  # calculating Sipe efficiency for many materials
  
  print "Calculating efficiency for all (kx, ky) values at wavelength "+title+"."

  print "kxx shape = "+str(kxx.shape)+"."
  etaSipe = np.zeros(kxx.shape)

  materialIndex = 0
  print "Preparing 2D figure for material "+str(Material2[materialIndex])
  print epsilon2[materialIndex]

  for m in np.arange(0,(kx.size),1):
	  #idy=0
	  for n in np.arange(0,ky.size,1):
		  #print "[Debug]"+str(m)+", "+str(n)
		  kappa = np.array([kx[m], ky[n]])
		  kappai = np.array([-cmath.sin(theta), 0])
		  kappap = kappai + kappa; kappam = kappai - kappa
		  etaSipe[m,n] = etap(theta, f, s, epsilon2[materialIndex], kappa, kappap, kappam)
		  #idy=idy+1
	  #idx=idx+1	

  print etaSipe

  maximum = np.amax(etaSipe)
  print maximum
  print Header+"Plot the graph for one given dielectric permittivity."
  plt.figure()
  levels = np.arange(0,maximum,maximum/numberlevels)
  CS = plt.contourf(kxx, kyy, etaSipe, levels=levels, cmap=plt.cm.Blues)
  plt.title(query+"/"+query2+r": $\lambda=$"+str(int(wavelength/unit))+" nm	")
  plt.xlabel(r'$\kappa_x$')
  plt.ylabel(r'$\kappa_y$')
  plt.colorbar(CS)
  filename='Sipe2d'+str(wavelength/unit)+'nm-'+query2
  plt.savefig(filename+'.eps')
  plt.savefig(filename+'.png')
  plt.show()
#}}}

#====================== GENERIC PLOTTING of the Sipe model =================
print Header+"Plot the kappaX for which maximum efficiency is found as function of dielectric permittivity. "
# This could help to localize problems and limitations of the Sipe theory. 

# Generating kappaX, kappaY meshes. 
print "Mesh generation..."
k_precision = 0.5
SipeRanges = 2e0
#kx = np.arange(0.,SipeRanges,k_precision)
#ky = np.arange(0.,SipeRanges,k_precision)
kx = [1.5e0]; ky= [0.0e0] #single value of kx,ky
#kxx, kyy = np.meshgrid(ky, kx)

# calculating Sipe efficiency for many materials
title = select+' nm'

# Generating mapping of dielectric permittivities
epsR = np.arange(-20, 2., 0.05) #eta: precision on epsilon inherited from libSPP.py
epsI = np.arange(  0., 10., 0.05) #eta: precision on epsilon inherited from libSPP.py

epsI2, epsR2 = np.meshgrid(epsI, epsR)
epsilon2 = np.add(epsR2,np.multiply(1e0j, epsI2)) #map of all possible dielectric permittivities

print "Calculating efficiency for all (kx, ky) values at wavelength "+title+"."

#print "kxx shape = "+str(kxx.shape)+"."
etaSipe = np.zeros((len(kx), len(ky), len(epsR), len(epsI)))
print Header+"Memory usage: "+str(len(kx)*len(ky)*len(epsR)*len(epsI)*64./8./1024./1024.)+" MB."
#print Header+"Memory usage: "+str(len(etaSipe)*64./8./1024.)+" kB."

print Header+"Shape (kx,ky,epsR,epsI)="+str(np.shape(etaSipe))

# Building the etaSipe(kx,ky) distribution
for j in np.arange(0,len(epsR),1):
  for k in np.arange(0,len(epsI),1):
    for m in np.arange(0,len(kx),1):
      for n in np.arange(0,len(ky),1):
        #print "[Debug] kx["+str(m)+"], ky["+str(n)+"], epsR["+str(j)+"], epsI["+str(k)+"]."
        kappa = np.array([kx[m], ky[n]])
        kappai = np.array([-cmath.sin(theta), 0])
        kappap = kappai + kappa; kappam = kappai - kappa
        try:
          etaSipe[m,n,j,k] = etap(theta, f, s, epsilon2[j,k], kappa, kappap, kappam)
        except:
          etaSipe[m,n,j,k] = 0.

# Ok, it's time to get a picture mapping of the efficiency.
#print etaSipe
maximum = np.amax(etaSipe) #finds maximum value of efficiency
print etaSipe[0,0,:,:]
print Header+"Maximum efficiency = "+str(maximum)
# TODO: find the corresponding value of kappaX, kappaY. 
plt.figure()
plt.matshow(etaSipe[:,:,0,0]) #eta(kx,ky;epsR=-2, epsI=0)
plt.colorbar()
plt.figure()
levels = [-2, -1, 0, 1, 2]
CS=plt.contourf(epsR2,epsI2,np.log10(etaSipe[0,0,:,:]),levels=levels, cmap=plt.cm.RdBu_r)
plt.colorbar(CS)
plt.xlabel(r'Re($\varepsilon$)')
plt.ylabel(r'Im($\varepsilon$)')
plt.show()

# Or we can integrate on a certain range of kx in 1. to 1.10. 
# TODO: About LIPSS regularity: Efficiency factor is maybe not the best quantity to look at, as highest factor is obtained for materials where -Re(eps) ~ Im(eps). If efficacy factor correspond to field enhancement, why do we find a low coupling with Au and Ag? And high coupling with W and Ti ? 
# TODO: About LIPSS period: Shall we automatically capture the (kx,ky) where efficiency is the highest in the Sipe(kx,ky)? Then we could plot Most_probable_period ( Re(eps) , Im(eps) ). 

# TODO 
# - Automatize the inverse Fourier transform to check regularity and pattern shape: see formula in my thesis. 
# - Automatic calculation of orientation angle precision
# - Can be great to plot directly precision angle as a function of materials. 