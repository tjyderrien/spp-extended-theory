#!/usr/bin/env python2
#-*- coding: utf-8 -*-

# Copyright (C) 2013-2018 T. J.-Y. Derrien
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
from matplotlib.ticker import MaxNLocator
import sys

# IMPORT CUSTOM LIBRARIES

# from libKeldysh import *
from libDatabase import *
from libLaser import *
from libMaterials import *
from libMath import *
from libSPP import *
from libSipe import *
from libPlotting import *
#from plotGraph import *

Header="[plotSipe.py] "
numberlevels = 100 #12 #for the final 2D plot

# Size of main canvas
SizeX = 8.
SizeY = 6.

#plt.figure(figsize=(SizeX,SizeY))

## Compute efficacy factor for a table of filling factors. 
#@param theta: angle of incidence (rad) 
#@param f: Filling factor
#@param s: Shape factor
def plotSipe1D_sectionX(wavelength, epsilon, f=0.1e0, s=0.4e0, theta=0e0): #{{{
    # Meshes for solution
    #kappax = np.arange(0, 4, 0.1)
    #kappay = np.arange(0, 4, 0.1)
    #kappa = np.array([wavelength * 1, wavelength * 0]); #test values
    #kappax = 4e0 ; kappay = 0e0*kappax; #test values
    kappa_max = 4E0
    kappa_min = 0.1E0
    numberofkpoints = 1000 #NOTE: be generous here, otherwise peaks will not be well resolved. 
    #for kappax in meshkappa:
    #ftab = np.arange(0, 1, 0.1)
    print "Info: Generating the mesh..."
    kapparange = np.arange(kappa_min,kappa_max,(kappa_max-kappa_min)/numberofkpoints)
    #for wavelength in wavelengths:
    #for f in ftab:
    idtab = 0
    etaSresult = np.zeros(kapparange.shape)
    etaPresult = np.zeros(kapparange.shape)
    kappax = 0e0
    for kappay in kapparange:
        ## Defining simple quantities for Sipe model
        kappa = np.array([kappax, kappay])
        kappai = np.array([-cmath.sin(theta), 0.])
        kappap = kappai + kappa; kappam = kappai - kappa

        #print "kappax = "+str(kappax)
        #print idtab
        etaPresult[idtab] = etap(theta, f, s, epsilon, kappa, kappap, kappam)
        etaSresult[idtab] = etas(theta, f, s, epsilon, kappa, kappap, kappam)
        #print etaresult[idtab]
        #print "eta = "+str(etaresult)
        idtab = idtab+1

    #=========== Make a 1D plot

    plt.figure(figsize=(SizeX,SizeY))
    plt.xlabel(r'$\kappa_x$')
    plt.ylabel(r'$\eta$')
    print kapparange.shape, etaSresult.shape
    plt.plot(kapparange, etaPresult, '-', label=r'$\eta_P$')
    plt.plot(kapparange, etaSresult, '-', label=r'$\eta_S$')
    print etaSresult
    plt.grid()
    plt.legend()
    plt.savefig('SipeEtaKappaX.eps')
    plt.show()
    #exit()
#}}}

#=========== 2D plot

#query2= 'InP (Bonse 2005)'
#query2= 'Mo (Ordal 1988)'
#query2= 'Cu (Palik)'
#query2= 'SiO2 (Palik)'
print "Caution: the values of queries must be exactly the one of MaterialDatabase.csv."

## Prepares the classical SipeEfficiencyFactor(kx, ky) for a specific material query2 immersed in Air. 
# @param wavelength: photon energy given in SI (m)
# @param query2: <string> linking to a material given in ../MaterialDatabase.csv.
# @param numberofkpoints: be generous here. At least 100 points are required per dimension. 
# NOTE: Sipe model is limited to air-material interface. Cannot be used with water-material for example. 
#       For more advanced combinations of materials, see [T.J.-Y. Derrien et al, Journal of Optics 18, 115007 (2016)]
def plotSipeFromDatabase(select, query2, wavelength, numberofkpoints=100, query='Air', theta=0, f=0.1, s=0.4): #{{{
  ## Generate the database
  SPPdb = GenerateDatabase() #Generate from MaterialDatabase.csv
  print "SPP database has "+str(len(SPPdb))+" entries."

  #print "Full Database:"
  #print SPPdb

  # Select the material of interface 1
  SPPdb = FilterDatabase(SPPdb, query, 0)
  print "Filter on materials: SPP database has now "+str(len(SPPdb))+" entries."

  # Filter database on materials
  try: 
    title = query2
    SPPdb = FilterDatabase(SPPdb, query2, 1)
    print "Filter on material: SPP database "+title+" has "+str(len(SPPdb))+" entries."
    print Header+"Check the level of tolerance in libSPP.py."
    #print SPPdb
  except:
    print "Exception: no optical data is available for "+query+" at "+title+"."
    print SPPdb
    exit()

  print SPPdb #works well
  # Filter database on wavelength
  title = select+' nm'
  try: 
    SPPdb = FilterDatabase(SPPdb, select, 2)
    print "Filter on wavelength: SPP database "+title+" has "+str(len(SPPdb))+" entries."
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
  
  SipeRanges = 4e0
  kx = np.arange(-SipeRanges,SipeRanges,(2.*SipeRanges)/numberofkpoints)
  ky = np.arange(-SipeRanges,SipeRanges,(2.*SipeRanges)/numberofkpoints)
  kxx, kyy = np.meshgrid(ky, kx)

  # calculating Sipe efficiency for many materials
  
  print "Calculating efficiency for all (kx, ky) values at wavelength "+title+"."

  print "kxx shape = "+str(kxx.shape)+"."
  etaSipe = np.zeros(kxx.shape)

  materialIndex = 0
  print "Preparing 2D figure for material "+str(Material2[materialIndex])
  print epsilon2[materialIndex]
  matrixsize = kx.size*ky.size
  print Header+"Size: "+str(matrixsize)
  index = 0
  for m in np.arange(0,(kx.size),1):
	  #idy=0
	  for n in np.arange(0,ky.size,1):
		  #print "[Debug]"+str(m)+", "+str(n)
		  kappa = np.array([kx[m], ky[n]])
		  kappai = np.array([-cmath.sin(theta), 0])
		  kappap = kappai + kappa; kappam = kappai - kappa
		  etaSipe[m,n] = etap(theta, f, s, epsilon2[materialIndex], kappa, kappap, kappam)
		  #idy=idy+1
		  index = index + 1
                  print Header+"** Progress: "+str(round(float(index)/float(matrixsize)*100.))+" percents."
  print etaSipe

  maximum = np.amax(etaSipe)
  print maximum
  print Header+"Plot the graph for one given dielectric permittivity."
  plt.figure(figsize=(SizeX,SizeY))
  levels = np.arange(0,maximum,maximum/numberlevels)
  CS = plt.contourf(kxx, kyy, etaSipe, levels=levels, cmap=plt.cm.Greys) #plt.cm.Blues #plt.cm.binary
  plt.title(query+"/"+query2+r": $\lambda=$"+str(int(wavelength/unit))+" nm	")
  plt.xlabel(r'$\kappa_x$')
  plt.ylabel(r'$\kappa_y$')
  plt.colorbar(CS)
  filename='Sipe2d'+str(wavelength/unit)+'nm-'+query2
  plt.savefig(filename+'.eps')
  plt.savefig(filename+'.png')
  plt.show()
#}}}

## Prepare a generalized plot2d of the Sipe model for any material at equilibrium. 
# This could help to localize problems and limitations of the Sipe theory. 
def plotGenericSipeMaps(kx_value, ky_value, epsilon_precision = 0.05, theta=0., f=0.1, s=0.4, PlotMaterials=False, PlotDrude=False, PlotLegend=False): #{{{
  print Header+"Plot the kappaX for which maximum efficiency is found as function of dielectric permittivity. "

  # Generating kappaX, kappaY meshes. 
  print "Mesh generation..."
  #k_precision = 0.5
  
  SipeRanges = 2e0
  #kx = np.arange(0.,SipeRanges,k_precision)
  #ky = np.arange(0.,SipeRanges,k_precision)

  # Manual definition of kappa_x, kappa_y. 
  #kx = [0.0e0, 0.8e0, 0.9e0, 1.0e0, 1.1e0, 1.2e0]; 
  kx = [kx_value]
  ky = [ky_value]
  #ky = [0.8e0, 0.9e0, 1.0e0, 1.1e0, 1.2e0] #single value of kx,ky
  #kxx, kyy = np.meshgrid(ky, kx)
  epsilon_real_min = -20. #-60 is more realistic
  epsilon_real_max = 10.   
  epsilon_imag_min = 0.
  epsilon_imag_max = 10.  #30 is more realistic
  # calculating Sipe efficiency for many materials
  title = select+' nm'

  # Generating mapping of dielectric permittivities
  epsR = np.arange(epsilon_real_min, epsilon_real_max, epsilon_precision) #eta: precision on epsilon could be inherited from libSPP.py by using the variable "precision". 
  epsI = np.arange(  epsilon_imag_min, epsilon_imag_max, epsilon_precision) #eta: precision on epsilon could be inherited from libSPP.py by using the variable "precision"
  
  
  # According to the paper: 
  wavelength = 1025e-9; epsilonSiO2 = (1.4504+0.j)**2; nuSiO2 = 0.4E-15**-1; meffSiO2 = 0.49; 
  neSiO2   = np.arange(0.,1E28,1E27)
  epsDrude_1 = Drude(wavelength, neSiO2, epsilonSiO2, 0.4e-15**-1, meffSiO2)
  epsDrude_2 = Drude(wavelength, neSiO2, epsilonSiO2, 0.5e-15**-1, meffSiO2)
  epsDrude_3 = Drude(wavelength, neSiO2, epsilonSiO2, 1.0e-15**-1, meffSiO2)
  
  epsI2, epsR2 = np.meshgrid(epsI, epsR)
  epsilon2 = np.add(epsR2,np.multiply(1e0j, epsI2)) #map of all possible dielectric permittivities

  print "Calculating efficiency for all (kx, ky) values at wavelength "+title+"."

  #print "kxx shape = "+str(kxx.shape)+"."
  number_iterations = len(kx)*len(ky)*len(epsR)*len(epsI)
  etaSipe = np.zeros((len(kx), len(ky), len(epsR), len(epsI)))
  print Header+"Memory usage: "+str(number_iterations*64./8./1024./1024.)+" MB."
  #print Header+"Memory usage: "+str(len(etaSipe)*64./8./1024.)+" kB."
  print Header+"Number of iterations: "+str(number_iterations)
  print Header+"Shape (kx,ky,epsR,epsI)="+str(np.shape(etaSipe))
  iteration_number=0
  # Building the etaSipe(kx,ky) distribution for all values of permittivities and all kappa_x, kappa_y.
  for j in np.arange(0,len(epsR),1):
    for k in np.arange(0,len(epsI),1):
      for m in np.arange(0,len(kx),1):
        for n in np.arange(0,len(ky),1):
          #print "[Debug] kx["+str(m)+"], ky["+str(n)+"], epsR["+str(j)+"], epsI["+str(k)+"]."
          kappa = np.array([kx[m], ky[n]])
          kappai = np.array([-cmath.sin(theta), 0.])
          kappap = kappai + kappa; kappam = kappai - kappa
          try:
            etaSipe[m,n,j,k] = etap(theta, f, s, epsilon2[j,k], kappa, kappap, kappam)
          except:
            etaSipe[m,n,j,k] = 0.
            print Header+"** Exception case was met."
          iteration_number = iteration_number + 1
          progress = 100.*iteration_number / number_iterations
    print Header+"** Progress: "+str(round(progress))+" percents."

  #print etaSipe
  maximum = np.amax(etaSipe) #finds maximum value of efficiency

  print Header+"** Info: Maximum value of efficacy: "+str(maximum)
  #print etaSipe[0,0,:,:]
  # TODO: find the value of kappaX, kappaY and epsilon for which efficiency is maximum. 

  # Ok, it's time to get a picture mapping of the efficiency.

  #plt.figure()
  #plt.contourf(etaSipe[:,:,0,0]) #eta(kx,ky;epsR=-2, epsI=0)
  #plt.xlabel(r'$\kappa_x$')
  #plt.ylabel(r'$\kappa_y$')
  #plt.title(r"$\eta (\kappa_x, \kappa_y; \varepsilon_r=$"+str(epsR[0])+r", $\varepsilon_i=$"+str(epsI[0])+")")
  #plt.colorbar()
  # ============================== ADDING THE CORRESPONDING MATERIALS TO THE PLOT ====================
  ## Generate the database
  SPPdb = GenerateDatabase()
  print "** Info: SPP database has "+str(len(SPPdb))+" entries."
  
  #print "Full Database:"
  #print SPPdb
  
  # TODO: if file OxideList.dat is provided, then we can look for couples, instead of generating the list of materials by ourselves. 
  
  # Select the material of interface 1 #TODO: This selector may be not clear for users. 
  query = "Air"
  SPPdb = FilterDatabase(SPPdb, query, 0)
  niceWavelength = str(int(wavelength*1e9))
  selectWavelength = str(niceWavelength)+'.0'
  title = niceWavelength+' nm'
  if(PlotMaterials==True): 
    try: 
        SPPdb = FilterDatabase(SPPdb, selectWavelength, 2) #BUG: crashes when treating Air
        print "Filter on wavelength: SPP database "+title+" has "+str(np.shape(SPPdb))+" entries."
        #plotDatabasePeriod(SPPdb, title, 'Period'+niceWavelength+'nm.eps', title, metal)
    except:
        print "Warning: no optical data is available for "+query+" at "+title+"."
    
    SPPdbSave = SPPdb
    print "** Limiting the Re(epsilon) space maximum ..."
    try: 
        SPPdb_filtered = FilterDatabaseLowerThan(SPPdb, epsilon_real_max, 15) #Field 16: eps.real
        SPPdb = SPPdb_filtered[0][:][:] #shape is getting one extra dimension for nothing! 
        print "Filter on wavelength: SPP database "+title+" has "+str(np.shape(SPPdb))+" entries."
        print SPPdb
        #plotDatabasePeriod(SPPdb, title, 'Period'+niceWavelength+'nm.eps', title, metal)
    except:
        print "Warning: no optical data is available for "+query+" at "+title+"."
        exit()
    
    #print "** Limiting the Re(epsilon) space minimum ..."
    #try: 
        #SPPdb_filtered = FilterDatabaseGreaterThan(SPPdb, epsilon_real_min, 15)
        #SPPdb = SPPdb_filtered[0][:][:] #shape is getting one extra dimension for nothing! 
        #print "Filter on wavelength: SPP database "+title+" has "+str(len(SPPdb))+" entries."
        #print SPPdb
        ##plotDatabasePeriod(SPPdb, title, 'Period'+niceWavelength+'nm.eps', title, metal)
    #except:
        #print "Warning: no optical data is available for "+query+" at "+title+"."
        #exit()
    print "** Limiting the Im(epsilon) space maximum ..."
    try: 
        SPPdb_filtered = FilterDatabaseGreaterThan(SPPdb, epsilon_imag_max, 16)
        SPPdb = SPPdb_filtered[0][:][:] #shape is getting one extra dimension for nothing! 
        print "Filter on wavelength: SPP database "+title+" has "+str(len(SPPdb))+" entries."
        #plotDatabasePeriod(SPPdb, title, 'Period'+niceWavelength+'nm.eps', title, metal)
    except:
        print "Warning: no optical data is available for "+query+" at "+title+"."
    
    print "** Limiting the Im(epsilon) space minimum ..."
    try: 
        SPPdb_filtered = FilterDatabaseGreaterThan(SPPdb, epsilon_imag_min, 16)
        SPPdb = SPPdb_filtered[0][:][:] #shape is getting one extra dimension for nothing! 
        print "Filter on wavelength: SPP database "+title+" has "+str(len(SPPdb))+" entries."
        #plotDatabasePeriod(SPPdb, title, 'Period'+niceWavelength+'nm.eps', title, metal)
    except:
        print "Warning: no optical data is available for "+query+" at "+title+"."
    
  
  #SPPdb = FilterDatabase(SPPdb, selectWavelength, 2)
  
  print "** Filtering materials: SPP database has now "+str(np.shape(SPPdb))+" entries."
  
  if(len(SPPdb)==0):
    print "** QUITTING..."
    exit()
  #print Header+"** Saving the materials database. "
  
  print Header+"** Preparing the list of materials"
  #print SPPdb
  
  Material1, Material2, Wavelength, OldSPPactiveBool, NewSPPactiveBool, SPPperiod, SPPperiodError, SPPdecayDepth1, SPPdecayDepth2, Reflectivity, OpticalPenetration1, OpticalPenetration2, SPPdecayLength, eps1rM, eps1cM, eps2rM, eps2cM, k1imag, k2imag, DeltaLsppValue = ExtractDataDb(SPPdb)
  
  # Calculation of refractive index
  eps1rM=np.asfarray(eps1rM)
  eps1cM=np.asfarray(eps1cM)
  eps2rM=np.asfarray(eps2rM)
  eps2cM=np.asfarray(eps2cM)

  plt.figure(figsize=(SizeX,SizeY))
  #levels = np.arange(0,10,1) 
  levels = np.arange(-2,4,1)
  CS=plt.contourf(epsR2,epsI2,np.log10(etaSipe[0,0,:,:]),levels=levels, cmap=plt.cm.RdYlBu_r) #RdBu_r: this one has white in center. Uncool. #Greys
  plt.xlabel(r'Re($\varepsilon$)')
  plt.ylabel(r'Im($\varepsilon$)')
  
  # adding the dots for materials
  #if(reverse): #{
  # plt.plot(eps1rM, eps1cM, 'or', label=str(wavelength)+'nm', markersize=8)
  #else:
  if(PlotMaterials==True): 
    plt.plot(eps2rM, eps2cM, 'or', label=str(wavelength)+'nm', markersize=8)
  #}
  if(PlotDrude==True): 
    plt.plot(epsDrude_1.real, epsDrude_1.imag, 'r--', label=r'Drude, $\tau_{D}=0.4$ fs', markersize=8)
    plt.plot(epsDrude_2.real, epsDrude_2.imag, 'g--', label=r'Drude, $\tau_{D}=0.5$ fs', markersize=8)
    plt.plot(epsDrude_3.real, epsDrude_3.imag, 'b--', label=r'Drude, $\tau_{D}=1.0$ fs', markersize=8)
  #}
  if(PlotLegend): 
    plt.legend(loc=2)
  #CS=plt.contourf(epsR2,epsI2,(etaSipe[0,0,:,:]),levels=levels, cmap=plt.cm.Greys) #RdBu_r
  plt.colorbar(CS)
  plt.title(r"$\eta (\kappa_x=$"+str(kx[0])+r"$, \kappa_y=$"+str(ky[0])+r"$;\varepsilon_r, \varepsilon_i$)")
  plt.xlabel(r'Re($\varepsilon$)')
  plt.ylabel(r'Im($\varepsilon$)')
  filename = "GrafDerrien-GeneralizedSipe-kx-"+str(kx[0])+"-ky-"+str(ky[0])
  print Header+"** Writing file: "+filename 
  plt.savefig(filename+".eps")
  plt.savefig(filename+".png")
  #plt.show()
  # Or we can integrate on a certain range of kx in 1. to 1.10. 
  # TODO: About LIPSS regularity: Efficiency factor is maybe not the best quantity to look at, as highest factor is obtained for materials where -Re(eps) ~ Im(eps). If efficacy factor correspond to field enhancement, why do we find a low coupling with Au and Ag? And high coupling with W and Ti ? 
  # TODO: About LIPSS period: Shall we automatically capture the (kx,ky) where efficiency is the highest in the Sipe(kx,ky)? Then we could plot Most_probable_period ( Re(eps) , Im(eps) ). 

  # TODO 
  # - Automatize the inverse Fourier transform to check regularity and pattern shape: see formula in my thesis. 
  # - Automatic calculation of orientation angle precision
  # - Can be great to plot directly precision angle as a function of materials. 
  
  return 0
#}}}

wavelength = 1026.0 #1064.0
select = str(int(wavelength))
request = select+".0"
unit = 1E-9
#wavelength = wavelength * unit
#print "Wavelength = "+str(wavelength/unit)+" nm."
kpointnumber = 500

# ======================= VALIDATION CASES ======================

# Bonse 2005 provides a 1D plot that can be compared quantitatively. 
# Bonse, J.; Munz, M. & Sturm, H. Structure formation on the surface of indium phosphide irradiated by femtosecond laser pulses J. Appl. Phys., 2005, 97, 013538
#plotSipe1D_sectionX(800e-9, 11.9296534741+1.4568728233j) #c-InP at 800 nm [Bonse2005]
#plotSipe1D_sectionX(800e-9, 12.21+1.4j) #c-InP at 800 nm [Bonse2005] NOTE: optical refractive index given in caption is not accurate, although a reference to Palik has been indicated. Data from pure Palik look to match better with the results provided in the article. 
#plotSipe1D_sectionX(800e-9, 14.4+1.52j)  #a-InP at 800 nm [Bonse2005]

# Generalization of this result 
# plotSipeFromDatabase(request, "InP (Palik)", 800e-9, kpointnumber)

#executing typical Sipe figure (like in [Bonse et al, Journal of Applied Physics (2009)]
#plotSipeFromDatabase(request, "Si (Palik)", 800e-9, kpointnumber)
# NOTE: Results are qualitatively okay', but it remains difficult to be sure of the complete repetition of the obtained ones. It might originate from the lack of decimals in the approximation of optical refractive index. 2 decimals are not sufficient for plasmonics. 
# NOTE: Colormap is also not provided in the manuscript. It is difficult to re-use the same mapping. 


# This paper provides Sipe maps computed by FDTD, but no amplitude is given on the efficacy factor. 
# Zhang, H.; Colombier, J.-P.; Li, C.; Faure, N.; Cheng, G. & Stoian, R. Coherence in ultrafast laser-induced periodic surface structures Physical Review B, 2015, 92. 

#TODO: Repeat figures from Colombier et al

# ======================= SCIENTIFIC PRODUCTION DATA =====================
filling = 0.1
shape = 0.4
angle = 0.
#plotSipeFromDatabase(request, "Cr (Johnson 1974)", 1026e-9, kpointnumber, 'Air', angle, filling, shape)
epsCr = -0.6721223+24.8657476j
nCr = np.sqrt(epsCr)-3.5j
print nCr
plotSipe1D_sectionX(1026e-9, nCr**2, filling, shape)

def plotStephanGraf_Materials2018(): #{{{
    #maps prepared for Stephane Gräf on generalized Sipe model (2018)
    k_precision    = 0.05
    filling_factor = 0.1
    shape_factor   = 0.4
    PlotMaterials  = False
    PlotDrude      = True

    plotGenericSipeMaps(0.8, 0.0, k_precision, 0., 0.1, 0.4, PlotMaterials, PlotDrude, True)
    plotGenericSipeMaps(0.9, 0.0, k_precision, 0., 0.1, 0.4, PlotMaterials, PlotDrude, False)
    plotGenericSipeMaps(1.0, 0.0, k_precision, 0., 0.1, 0.4, PlotMaterials, PlotDrude, False)
    plotGenericSipeMaps(1.1, 0.0, k_precision, 0., 0.1, 0.4, PlotMaterials, PlotDrude, False)
    plotGenericSipeMaps(1.2, 0.0, k_precision, 0., 0.1, 0.4, PlotMaterials, PlotDrude, False)
    plotGenericSipeMaps(0.0, 0.8, k_precision, 0., 0.1, 0.4, PlotMaterials, PlotDrude, False)
    plotGenericSipeMaps(0.0, 0.9, k_precision, 0., 0.1, 0.4, PlotMaterials, PlotDrude, False)
    plotGenericSipeMaps(0.0, 1.0, k_precision, 0., 0.1, 0.4, PlotMaterials, PlotDrude, False)
    plotGenericSipeMaps(0.0, 1.1, k_precision, 0., 0.1, 0.4, PlotMaterials, PlotDrude, False)
    plotGenericSipeMaps(0.0, 1.2, k_precision, 0., 0.1, 0.4, PlotMaterials, PlotDrude, False)
#}}}

#plotSipe1D_sectionX(1030e-9, -0.6721223+24.8657476j)
