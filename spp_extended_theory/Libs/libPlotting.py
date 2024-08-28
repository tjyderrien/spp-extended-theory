#!/usr/bin/env python
#-*- coding: utf-8 -*-

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

## @package libPlotting
## Provides several advanced plotting features.
# Library to plot various results from the simulations. 
# Also, plotting using annotations on each point is also available. 
# Prepared from source: http://stackoverflow.com/questions/8850142/matplotlib-overlapping-annotations

import matplotlib.pyplot as plt
import matplotlib.backends.backend_ps
import matplotlib.backends.backend_svg
import numpy as np
import spp_extended_theory.Libs.libDatabase as libDatabase
import spp_extended_theory.Libs.libMaterials as libMaterials
import sys


## Defines text positions for plotting
#  Original author and descriptinons are given here:
#  http://stackoverflow.com/questions/8850142/matplotlib-overlapping-annotations
def get_text_positions(x_data, y_data, txt_width, txt_height):

  a = list(zip(y_data, x_data))
  text_positions = y_data.copy()
  for index, (y, x) in enumerate(a):
    local_text_positions = [i for i in a if i[0] > (y - txt_height) 
	        and (abs(i[1] - x) < txt_width * 2) and i != (y,x)]
    if local_text_positions:
      sorted_ltp = sorted(local_text_positions)
      if abs(sorted_ltp[0][0] - y) < txt_height: #True == collision
        differ = np.diff(sorted_ltp, axis=0)
        a[index] = (sorted_ltp[-1][0] + txt_height, a[index][1])
        text_positions[index] = sorted_ltp[-1][0] + txt_height
        for k, (j, m) in enumerate(differ):
          #j is the vertical distance between words
          if j > txt_height * 2: #if True then room to fit a word in
            a[index] = (sorted_ltp[k][0] + txt_height, a[index][1])
            text_positions[index] = sorted_ltp[k][0] + txt_height
            break
  return text_positions

## Prepares the arrows for plotting using overlapping annotations
def text_plotter(x_data, y_data, text_content, text_positions, axis,txt_width,txt_height, color):
    for x,y,s,t in zip(x_data, y_data, text_content, text_positions):
        axis.text(x - txt_width, 1.01*t, s, rotation=0, color=color)
        if y != t:
            axis.arrow(x, t, 0, y-t, color='grey', alpha=0.3, width=5*0.01,
                       head_width=3*0.2, head_length=0.05*txt_width*0.5,
                       zorder=0, length_includes_head=True)

## Generate a plot with text labels on each point. 
# Useful to address many materials in the same figure. 
def makePlot(x_data, y_data, tags, filename, plottitle, labelx, labely, functionlabel, textcolor):
  #random test data:
  #x_data = random_sample(100)
  #y_data = random_integers(10,50,(100))

  #GOOD PLOT:
  fig2 = plt.figure()
  ax2 = fig2.add_subplot(111)
  ax2.plot(x_data, y_data, 's'+textcolor, markersize=8, label=functionlabel)
  plt.xlabel(labelx)
  plt.ylabel(labely)
  ax2.grid()
  plt.title(plottitle)
  plt.legend()
  #set the bbox for the text. Increase txt_width for wider text.
  txt_height = 0.08*(plt.ylim()[1] - plt.ylim()[0])
  txt_width = 0.01*(plt.xlim()[1] - plt.xlim()[0])
  #Get the corrected text positions, then write the text.
  text_positions = get_text_positions(x_data, y_data, txt_width, txt_height)
  text_plotter(x_data, y_data, tags, text_positions, ax2, txt_width, txt_height, textcolor)

  # plt.ylim(0.,max(text_positions)+2*txt_height)
  #plt.xlim(-0.1,1.1)
  
  plt.savefig(filename)
  plt.show()
  return 0
  
## Plots the visual list (ReEps, ImEps) of materials present in a database, for a given wavelength. 
#  @param database:   the database with data to plot
#  @param legend:     legend string to put onto the plot
#  @param outputfile: name of the file to output
#  @param query:      
def plotDatabaseMaterials(database, legend, outputfile, query, metal): #{{{

  # Material1, Material2, Wavelength, OldSPPactiveBool, NewSPPactiveBool,  # 5
  # SPPperiod, SPPperiodError, SPPdecayDepth1, SPPdecayDepth2, Reflectivity,  # 10
  # OpticalPenetration1, OpticalPenetration2, SPPdecayLength, eps1r, eps1c,  # 15
  # eps2r, eps2c, k1imag, k2imag, deltaLsppValues,  # 20
  # Fa, Jo, LifeTime

  # Unfold the data from database
  if (not metal):
    Material1, Material2, Wavelength, OldSPPactiveBool, NewSPPactiveBool, SPPperiod, SPPperiodError, SPPdecayDepth1, \
    SPPdecayDepth2, Reflectivity, OpticalPenetration1, OpticalPenetration2, SPPdecayLength, eps1r, eps1c, eps2r, eps2c, \
    k1imag, k2imag, DeltaLsppValue, Fa, Jo, LifeTime = libDatabase.ExtractDataDb(database)
  else: 
    Material2, Material1, Wavelength, OldSPPactiveBool, NewSPPactiveBool, SPPperiod, SPPperiodError, SPPdecayDepth1, \
    SPPdecayDepth2, Reflectivity, OpticalPenetration1, OpticalPenetration2, SPPdecayLength, eps2r, eps2c, eps1r, eps1c, \
    k1imag, k2imag, DeltaLsppValue, Fa, Jo, LifeTime = libDatabase.ExtractDataDb(database)
  
  # Clean the first field
  Material2clean = libDatabase.CleanStrArray(Material2)
  
  # Prepare plot with arrows and text (but single wavelength)
  makePlot(eps2r, eps2c, Material2clean, outputfile, query, r'$Re(\varepsilon)$', r'$Im(\varepsilon)$', legend, 'r')
  #makePlot(eps2r, SPPdecayLength*1e-3, Material2clean, outputfile, query, r'$Re(\varepsilon)$', 'SPP decay length (um)', legend, 'r')

  return 0
#}}}

## Plots the SPP-period for interfaces given in the database, for a given wavelength. 
def plotDatabasePeriod(database, legend, outputfile, query, metal): #{{{
  # Unfold data from database

  # Material1, Material2, Wavelength, OldSPPactiveBool, NewSPPactiveBool,  # 5
  # SPPperiod, SPPperiodError, SPPdecayDepth1, SPPdecayDepth2, Reflectivity,  # 10
  # OpticalPenetration1, OpticalPenetration2, SPPdecayLength, eps1r, eps1c,  # 15
  # eps2r, eps2c, k1imag, k2imag, deltaLsppValues,  # 20
  # Fa, Jo, LifeTime

  if (not metal):
    Material1, Material2, Wavelength, OldSPPactiveBool, NewSPPactiveBool, SPPperiod, SPPperiodError, SPPdecayDepth1, \
    SPPdecayDepth2, Reflectivity, OpticalPenetration1, OpticalPenetration2, SPPdecayLength, eps1r, eps1c, eps2r, eps2c, \
    k1imag, k2imag, DeltaLsppValue, Fa, Jo, LifeTime = libDatabase.ExtractDataDb(database)
  else: 
    Material2, Material1, Wavelength, OldSPPactiveBool, NewSPPactiveBool, SPPperiod, SPPperiodError, SPPdecayDepth1, \
    SPPdecayDepth2, Reflectivity, OpticalPenetration1, OpticalPenetration2, SPPdecayLength, eps2r, eps2c, eps1r, eps1c, \
    k1imag, k2imag, DeltaLsppValue, Fa, Jo, LifeTime = libDatabase.ExtractDataDb(database)
  
  # Clean the first field	
  Material2clean = libDatabase.CleanStrArray(Material2)
  
  # Prepare plot with arrows and text (but single wavelength)
  makePlot(eps2r, SPPperiod, Material2clean, outputfile, query, r'$Re(\varepsilon)$', 'Period (nm)', legend, 'r')
  
  return 0
#}}}

## Plots the SPP-period for interfaces given in the database, for a given wavelength.
def plotDatabaseLifetimeRaether(database, legend, outputfile, query, metal):  # {{{
  # Unfold data from database
  if (not metal):

    #Material1, Material2, Wavelength, OldSPPactiveBool, NewSPPactiveBool,  # 5
    #SPPperiod, SPPperiodError, SPPdecayDepth1, SPPdecayDepth2, Reflectivity,  # 10
    #OpticalPenetration1, OpticalPenetration2, SPPdecayLength, eps1r, eps1c,  # 15
    #eps2r, eps2c, k1imag, k2imag, deltaLsppValues,  # 20
    #Fa, Jo, LifeTime

    Material1, Material2, Wavelength, OldSPPactiveBool, NewSPPactiveBool, \
    SPPperiod, SPPperiodError, SPPdecayDepth1, SPPdecayDepth2, Reflectivity, \
    OpticalPenetration1, OpticalPenetration2, SPPdecayLength, eps1r, eps1c, \
    eps2r, eps2c, k1imag, k2imag, DeltaLsppValue, \
    Fa, Jo, LifeTime = libDatabase.ExtractDataDb(database)

  else:
    Material2, Material1, Wavelength, OldSPPactiveBool, NewSPPactiveBool, \
    SPPperiod, SPPperiodError, SPPdecayDepth1, SPPdecayDepth2, Reflectivity, \
    OpticalPenetration1, OpticalPenetration2, SPPdecayLength, eps2r, eps2c, \
    eps1r, eps1c, k1imag, k2imag, DeltaLsppValue, \
    Fa, Jo, LifeTime = libDatabase.ExtractDataDb(database)

  # Clean the first field
  Material2clean = libDatabase.CleanStrArray(Material2)

  # Prepare plot with arrows and text (but single wavelength)
  makePlot(eps2r, LifeTime, Material2clean, outputfile, query, r'$Re(\varepsilon)$', 'Lifetime Raether (s)', legend, 'r')

  return 0


# }}}

## Plots the SPP decay length L_{SPP} at various interfaces contained in a database
# This should be restricted to a single wavelength for clarity. 
def plotDatabaseLspp(database, legend, outputfile, query, metal): #{{{
  
  # Unfolding data from database
  if (not metal):
    Material1, Material2, Wavelength, OldSPPactiveBool, NewSPPactiveBool, SPPperiod, SPPperiodError, SPPdecayDepth1, \
    SPPdecayDepth2, Reflectivity, OpticalPenetration1, OpticalPenetration2, SPPdecayLength, eps1r, eps1c, eps2r, eps2c, \
    k1imag, k2imag, DeltaLsppValue, Fa, Jo, LifeTime = libDatabase.ExtractDataDb(database)
  else: 
    Material2, Material1, Wavelength, OldSPPactiveBool, NewSPPactiveBool, SPPperiod, SPPperiodError, SPPdecayDepth1, \
    SPPdecayDepth2, Reflectivity, OpticalPenetration1, OpticalPenetration2, SPPdecayLength, eps2r, eps2c, eps1r, eps1c, \
    k1imag, k2imag, DeltaLsppValue, Fa, Jo, LifeTime = libDatabase.ExtractDataDb(database)
  
  # CLean the first field	
  Material2clean = libDatabase.CleanStrArray(Material2)
  
  # Prepare plot with arrows and text (but single wavelength)
  makePlot(eps2r, SPPdecayLength*1e-3, Material2clean, outputfile, query, r'$Re(\varepsilon)$', 'SPP decay length (um)', legend, 'r')

  return 0
#}}}

## Plots uncertainty <deltaL_{spp}> on SPP decay length (L_{SPP}) at interfaces contained in the given database.
# For clarity, it is advised to reduce database to single wavelength. 
def plotDatabaseDeltaLspp(database, legend, outputfile, query, metal): #{{{
  if (not metal):
    Material1, Material2, Wavelength, OldSPPactiveBool, NewSPPactiveBool, SPPperiod, SPPperiodError, SPPdecayDepth1, \
    SPPdecayDepth2, Reflectivity, OpticalPenetration1, OpticalPenetration2, SPPdecayLength, eps1r, eps1c, eps2r, eps2c, \
    k1imag, k2imag, deltaLsppValue, Fa, Jo, LifeTime = libDatabase.ExtractDataDb(database)
  else: 
    Material2, Material1, Wavelength, OldSPPactiveBool, NewSPPactiveBool, SPPperiod, SPPperiodError, SPPdecayDepth1, \
    SPPdecayDepth2, Reflectivity, OpticalPenetration1, OpticalPenetration2, SPPdecayLength, eps2r, eps2c, eps1r, eps1c, \
    k1imag, k2imag, deltaLsppValue, Fa, Jo, LifeTime  = libDatabase.ExtractDataDb(database)
  # CLean the first field	
  Material2clean = libDatabase.CleanStrArray(Material2)
  
  # Prepare plot with arrows and text (but single wavelength)
  makePlot(eps2r, deltaLsppValue, Material2clean, outputfile, query, r'$Re(\varepsilon)$', r'$\delta L_{SPP}$ (nm)', legend, 'r')

  return 0
#}}}  
  
## Plot period as function of materials for two wavelengths. 
# Two types of data are included: 
# - Points: 1 point per material. 
# - Lines : the continuum calculation.
# 
# Metal mode: enable for materials found metallic on most of the wavelengths
# Reverse mode: invert the material indices from the SPP database. Useful some metallic materials. 
def plotSeveralWavelengths(database1, database2, reverse, metal, query):#{{{
  # Extract data for 800 nm
  Material1, Material2, Wavelength, OldSPPactiveBool, NewSPPactiveBool, SPPperiod, SPPperiodError, SPPdecayDepth1, \
  SPPdecayDepth2, Reflectivity, OpticalPenetration1, OpticalPenetration2, SPPdecayLength, eps1r, eps1c, eps2r, eps2c, \
  k1imag, k2imag, DeltaLsppValue, Fa, Jo, LifeTime = libDatabase.ExtractDataDb(database1)
  
  if(metal): #swap eps1 and eps2
    Material2, Material1, Wavelength, OldSPPactiveBool, NewSPPactiveBool, SPPperiod, SPPperiodError, SPPdecayDepth2, \
    SPPdecayDepth1, Reflectivity, OpticalPenetration2, OpticalPenetration1, SPPdecayLength, eps2r, eps2c, eps1r, eps1c, \
    k2imag, k1imag, DeltaLsppValue, Fa, Jo, LifeTime  = libDatabase.ExtractDataDb(database1)
  else: 
    Material1, Material2, Wavelength, OldSPPactiveBool, NewSPPactiveBool, SPPperiod, SPPperiodError, SPPdecayDepth1, \
    SPPdecayDepth2, Reflectivity, OpticalPenetration1, OpticalPenetration2, SPPdecayLength, eps1r, eps1c, eps2r, eps2c, \
    k1imag, k2imag, DeltaLsppValue, Fa, Jo, LifeTime  = libDatabase.ExtractDataDb(database1)
    
  # Calculation of refractive index array
  eps1r=np.asfarray(eps1r)
  eps1c=np.asfarray(eps1c)
  eps2r=np.asfarray(eps2r)
  eps2c=np.asfarray(eps2c)

  epsilon1 = np.add(eps1r,np.multiply(1e0j, eps1c))
  epsilon2 = np.add(eps2r,np.multiply(1e0j, eps2c))
  
  fig1=plt.figure()
  if (not metal):
    plt.xlabel(r'Dielectric permittivity: $\mathcal{R}e(\varepsilon_2)$')
  else:
    plt.xlabel(r'Dielectric permittivity: $\mathcal{R}e(\varepsilon_1)$')
  
  # Discontinuous plotting
  plt.ylabel('SPP period $\Lambda$ (nm)')
  plt.plot(eps2r, SPPperiod, 'or', label='800 nm', markersize=8)
  
  # Continuous plotting
  eps1range = epsilon1[0] #use external index
  eps2range = np.arange(-70,0e0,0.1e0)
  
  # Use only one object
  # TODO: reprogram the whole function with switches in function arguments
  if(not metal): 
    refractiveindex1 = libMaterials.EpsilonToIndex(epsilon1)
  else: 
    refractiveindex1 = libMaterials.EpsilonToIndex(epsilon2)
    
  Radiation1=Wavelength/refractiveindex1.real
  Radiation1=np.sort(Radiation1)
  LambdaPMA=Wavelength[1]/np.sqrt(np.multiply(eps1range, eps2range)/(np.add(eps1range, eps2range)))
  
  plt.plot(np.sort(eps2r), Radiation1[::-1], 'r--')
  plt.plot(eps2range, LambdaPMA, 'r-')
  
  # Extract (again) for 400 nm
  #TODO: use a function here! code is repeated! 
  if(metal): #swap eps1 and eps2
    Material2, Material1, Wavelength, OldSPPactiveBool, NewSPPactiveBool, SPPperiod, SPPperiodError, SPPdecayDepth2, \
    SPPdecayDepth1, Reflectivity, OpticalPenetration2, OpticalPenetration1, SPPdecayLength, eps2r, eps2c, eps1r, eps1c, \
    k2imag, k1imag, DeltaLsppValue, Fa, Jo, LifeTime  = libDatabase.ExtractDataDb(database2)
  else: 
    Material1, Material2, Wavelength, OldSPPactiveBool, NewSPPactiveBool, SPPperiod, SPPperiodError, SPPdecayDepth1, \
    SPPdecayDepth2, Reflectivity, OpticalPenetration1, OpticalPenetration2, SPPdecayLength, eps1r, eps1c, eps2r, eps2c, \
    k1imag, k2imag, DeltaLsppValue, Fa, Jo, LifeTime  = libDatabase.ExtractDataDb(database2)
  
  # Calculation of refractive index
  eps1r=np.asfarray(eps1r)
  eps1c=np.asfarray(eps1c)
  eps2r=np.asfarray(eps2r)
  eps2c=np.asfarray(eps2c)

  epsilon1 = np.add(eps1r,np.multiply(1e0j, eps1c))
  epsilon2 = np.add(eps2r,np.multiply(1e0j, eps2c))
  
  # Continuous plots
  #if(metal):
    #eps1range = np.arange(-70,0e0,0.1e0)
    #eps2range = epsilon1[0] #TODO: use external index
  #else:
  eps1range = epsilon1[0]
  eps2range = np.arange(-70,0e0,0.1e0)
  
  if(not metal): 
    refractiveindex1 = libMaterials.EpsilonToIndex(epsilon1)
  else: 
    refractiveindex1 = libMaterials.EpsilonToIndex(epsilon2)
  
  #TODO: make a function here! Some code is repeated! 
  #This is the lambda/n continuous curve - n: sqrt(epsilon) of material1
  Radiation1=Wavelength/refractiveindex1.real
  Radiation1=np.sort(Radiation1)
  
  #This is rather to plot lambda_PMA
  LambdaPMA=Wavelength[1]/np.sqrt(np.multiply(eps1range, eps2range)/(np.add(eps1range, eps2range)))

  #plt.figure()
  plt.plot(eps2r, SPPperiod, 'bs', label='400 nm', markersize=8)
  plt.plot(np.sort(eps2r), Radiation1[::-1], 'b--')
  plt.plot(eps2range, LambdaPMA, 'b-') #TODO: plot using a full range, not eps2r

  #print eps1r
  #TODO simplify + combine the 3 following tests
  if(metal):
    if(query=="Au (Palik)"):
      plt.axis([0,25,0,900]) ##KEEP 900 please #good for Au
    elif(query=="Ti (Johnson)"):
      plt.axis([0,20,0,900]) ##KEEP 900 please #good for Ti
    else:
      print("** Error: this query is not a planned case. Query="+query)
      #sys.exit()
  else: #non-metal
    plt.axis([-70,0,0,900]) ##KEEP 900 please
    
  plt.axis()
  plt.yticks([0,200,400,600,800])
  #plt.axis([0,35,0,1000])
  if(not metal): 
    if(query=="Air"):
      #plt.legend(loc=4) #Good for Air
      plt.legend(bbox_to_anchor=(0.95, 0.05), loc=4, borderaxespad=0.)
    elif(query=="SiO2 (Palik)"):
      plt.legend(loc=2) #Good for SiO_2
    else:
      #print "** Error: this query is not a planned case. Query="+query
      sys.exit()
  else: #metal case
    plt.legend(loc=1)
  #plt.legend(handler_map={line1: HandlerLine2D(numpoints=1)})
  #plt.legend(handler_map={line2: HandlerLine2D(numpoints=1)})
  #plt.title(query)
  plt.grid()
  
  # Adding a sub-plot
  if(not metal): 
    #a = plt.axes([-70,300,-40,700], axisbg='g')
    if(query=="SiO2 (Palik)"):
      a = plt.axes([0.2,0.18,0.35,0.35], axisbg='w') #Good for SiO2
      plt.axis([-6,0,250,290]) #Good for SiO2
      plt.yticks([250,270,290]) #Good for SiO2
    elif (query=="Air"):
      # Position of the plot
      a = plt.axes([0.2,0.17,0.35,0.35], axisbg='w') #Good for Air
      plt.axis([-6,0,380,410]) #Good for Air
      plt.yticks([380,390,400,410])
    else: #TODO: set a general case
      a = plt.axes([0.2,0.2,0.35,0.35], axisbg='w')

    plt.xticks([-6,-4,-2,0])
    
    plt.grid()
    plt.plot(eps2r, SPPperiod, 'bs', markersize=8)
    plt.plot(np.sort(eps2r), Radiation1[::-1], 'b--')
    #plt.title('Zoom')
    #plt.xticks([])
    #plt.yticks([])

  plt.savefig('MultiMaterial_PeriodSPP.eps')
  plt.show()
  return 0
#}}}
