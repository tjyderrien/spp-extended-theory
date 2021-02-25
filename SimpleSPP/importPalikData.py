#!/usr/bin/env python2
#-*- coding: utf-8 -*-

# Copyright (C) 2013-2020 T. J.-Y. Derrien
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
#from libSPP import *
#from libUnits import *
#from libDatabase import *
from libImportOpticalData import *

#==============================

if(len(sys.argv)<=2):
  print "Usage: ./importPalikData.py <Name of the material (Be, Au, ...)> <wavelength (nm)> <Source for data: Palik or name of the 1st author>"
  print "Example: ./importPalikData.py Au 800 Palik"
  exit()
 
material = sys.argv[1]
try:
  wavelength = float(sys.argv[2])
except: 
  print "** Error: command line must have the hape: <Material symbol> <Wavelength (nm)> <Publication author>"
  exit()
  
try:
  source = sys.argv[3]
except:
  print "** Warning: no data source given: using DEFAULT=Palik"
  source = "Palik"
  
filename = material+"-"+source
folder = ""
  
#filename = "Au-Palik"
#filename = "Be-Palik"
#filename = "Fe-Palik"
#filename = "Ni-Palik"
#filename = "a-Si-Palik"
#filename = "Si-Palik"
#filename = "SiC-Palik"
#filename = "Ti-Palik"
#filename = "Ag-Johnson"
#filename = "Ti-Johnson"
#filename = "SiO2-Palik"
#filename = "W-Palik"
#filename = "Cr-Johnson"
#filename = "BK7"
#filename = "Al-Palik"
plotting = True

#folder = "Database/PalikGraph/"
#filename = "c-Si-77K-Dash"
#filename = "c-Si-80K-Humlicek"
#filename = "Ti-Palik"
#filename="ZnO-Bond"

#=========================================

try: #TODO: should we select by author? Or by units? 
# NOTE: if Palik, then wavelength is given in nm. 
# NOTE: if other, then wavelengths are usually given in um. 
# BUG: access to databases are treated differently between importPalikData.py and plotMultiwavelength.py. 
  if(source == "Palik"):
    folder = "Database/"
    print "Material: "+filename+"."
    print "Wavelength = "+str(wavelength)+" nm"
    epsilon = importFromTable(wavelength*1e-9, folder, filename, plotting)
    print epsilon
    print "You can add the following directly inside 'MaterialOpticalDatabaseForPlasmonics.csv'"
    print filename+"\t"+"?"+"\t"+str(int(wavelength))+"\t"+str(epsilon.real)+"\t"+str(epsilon.imag)+"\t?\t?\t?\t?\t?"
  elif(source == "Miller"):
    folder = "Database/LiquidMetals/"
    filename = "l-Cu-Miller"
    epsilon = importFromEpsilonTable(wavelength*1E-9, folder, filename, plotting)
    print epsilon
    print "You can add the following directly inside 'MaterialOpticalDatabaseForPlasmonics.csv'"
    print filename+"\t"+"?"+"\t"+str(int(wavelength))+"\t"+str(epsilon.real)+"\t"+str(epsilon.imag)+"\t?\t?\t?\t?\t?"
  elif(source == "GoriAndBond"):
    folder = "Database/"
    filename = "ZnO-GoriAndBond"
    print "Material: "+filename+"."
    print "Wavelength = "+str(wavelength)+" nm"
    epsilon = importFromTable(wavelength*1e-9, folder, filename, plotting) #Careful: misleading neames. These data have been stored as dielectric permittivities, but this function inteprets it as nk. 
    #epsilon = nk**2
    #print nk
    print "You can add the following directly inside 'MaterialOpticalDatabaseForPlasmonics.csv'"
    print filename+"\t"+"?"+"\t"+str(int(wavelength))+"\t"+str(epsilon.real)+"\t"+str(epsilon.imag)+"\t?\t?\t?\t?\t?"
  elif(source == "Gori"):
    # We plot Gori data along with Bond data. 
    # Gori is meshed on (energy (eV), epsilon)
    # Bond is meshed on (wavelength (, n, k)
    print "** Warning: ZnO-Gori data must be completed around 4 eV using Bond data."
    print "WE NOW GENERATE THE INTERPOLATED FILES TO ALLOW FOR MERGING."
    print "To use the merged Gori-Bond files, type GoriAndBond source instead of Gori."
    print ""
    print "Importing ZnO-Gori data..."
    plotting = True
    folderGori   = "Database/"
    filenameGori = "ZnO-Gori"
    
    energiesGori, epsilonGori  = importFromEpsilonTable_batch(folderGori, filenameGori, plotting, 1E0) #unit in nm
    #print energiesGori, epsilonGori #output format: eV, epsilon
    # BUG: impossible to call Energy_SI_to_Length ?
    #wavelengthsGori = Energy_SI_to_Length(Energy_eV_to_Joules(energiesGori)) #BUG HERE
    energiesGori_J   = e*energiesGori
    #print energiesGori_J, epsilonGori #output format: J, epsilon
    wavelengthsGori  = h*c/(energiesGori_J) #output format: meters
    #print wavelengthsGori, epsilonGori
    nGori=np.sqrt(epsilonGori)
    GoriFile = np.array([wavelengthsGori*1E9, nGori.real, nGori.imag])
    ExportToTxt(np.flipud(np.transpose(GoriFile)), "ZnO-Gori.csv")
    print Header+"** Exported Gori file. "
    
    print "Importing more detailed ZnO-Gori data..."
    filenameGori2 = "ZnO-Gori2"
    
    energiesGori2, epsilonGori2  = importFromEpsilonTable_batch(folderGori, filenameGori2, plotting, 1E0) #unit in nm
    energiesGori_J2   = e*energiesGori2
    wavelengthsGori2  = h*c/(energiesGori_J2) #output format: meters
    
    nGori2 = np.sqrt(epsilonGori2)
    
    Gori2File = np.array([wavelengthsGori2*1E9, nGori2.real, nGori2.imag])
    ExportToTxt(np.flipud(np.transpose(Gori2File)), "ZnO-Gori2.csv")
    print Header+"** Exported Gori2 file. "
    #print "Combining the two Gori sets of data..."
    #wavelengths_Gori_final, epsilonGori_final = interpolateTwoSetsOfOpticalData(wavelengthsGori, wavelengthsGori2, epsilonGori, epsilonGori2)
    #print wavelengths_Gori_final, epsilonGori_final
    
    print ""
    print Header+"Info: Successfully imported ZnO-Gori data."
    print ""
    print Header+"Completing with ZnO-Bond data..."
    folderBond   = "Database/PalikGraph/"
    filenameBond = "ZnO-Bond"
    wavelengthsBond, nkBond       = importFromNKtable_batch(folderBond, filenameBond, plotting, 1E-6) #unit in um
    print wavelengthsBond, nkBond
    print Header+"Info: Imported ZnO-Bond data."
    epsilonBond  =  np.multiply(nkBond, nkBond) #converting (n,k) to (epsR, epsC)
    print Header+"Info: Converted ZnO-Bond to dielectric permittivity."
    BondFile = np.array([wavelengthsBond*1E9, nkBond.real, nkBond.imag])
    ExportToTxt(np.transpose(BondFile), "ZnO-Bond.csv")
    print Header+"** Exported Bond file. "
    
    plt.figure()
    plt.semilogx(wavelengthsGori, epsilonGori.real, "r-", label='Gori')
    plt.semilogx(wavelengthsGori2, epsilonGori2.real, "b-", label='Gori-2')
    plt.semilogx(wavelengthsBond, epsilonBond.real, "g-", label='Bond')
    plt.semilogx(wavelengthsGori, epsilonGori.imag, "r--", label='')
    plt.semilogx(wavelengthsGori2, epsilonGori2.imag, "b--", label='')
    plt.semilogx(wavelengthsBond, epsilonBond.imag, "g--", label='')
    plt.xlabel("Wavelength (nm)")
    plt.ylabel(r"$\varepsilon$")
    plt.legend()
    plt.grid()
    plt.savefig("ZnO-reconstructed.eps")
    plt.show()
    
    print "DONT FORGET TO ACCOLATE THE OPTICAL DATA TO A FILE ZnO-GariAndBond."
  elif(source == "Chase"):
    print "** Info: branching with CrO2 Chase optical data..."
    folder = "Database/"
    filename = "CrO2-Chase"
    importFromEpsilonTable(wavelength, folder, filename, True, 1E-6)
  elif(source == "Chase-X"):
    print "** Info: branching with CrO2 Chase optical data..."
    folder = "Database/"
    filename = "CrO2-Chase-X"
    importFromEpsilonTable(wavelength, folder, filename, True, 1E-6)
  else: #TODO: revise the design here. 
    #folder = "Database/PalikGraph/" #TODO: Ag-Johnson and BK7-Maliton are in ./Database actually.
    folder = "Database/"
    print "Material: "+filename+"."
    print "Wavelength = "+str(wavelength)+" nm"
    try: 
      epsilon=importFromTable(wavelength*1e-9, folder, filename, True) #NOTE: sometimes have to change the unit here. 
    except:
      print "Import failed using importFromTable(). Trying with importFromNKtable()."
      try: 
        importFromNKtable(wavelength*1e-6, folder, filename) #BUG: error in 
      except:
        print "Failed even using importFromNKtable(). Call that damn developer. "
        exit()
    print epsilon
    print "You can add the following directly inside 'MaterialOpticalDatabaseForPlasmonics.csv'"
    print filename+"\t"+"?"+"\t"+str(int(wavelength))+"\t"+str(epsilon.real)+"\t"+str(epsilon.imag)+"\t?\t?\t?\t?\t?"
except:
  print "Failed to import "+folder+filename+"!"
  print "Use SimpleSPP/Database/importPalikData.sh for finding other sources."
  
#==========================================
#print "Lambda = 3000 nm"
#importFromTable(3000e-9, folder, filename, plotting=True)
#print "Lambda = 1064 nm"
#importFromTable(1064e-9, folder, filename, plotting=True)
#print "Lambda = 1060 nm"
#importFromTable(1060e-9, folder, filename, plotting)
#print "Lambda = 1030 nm"
#importFromTable(1030e-9, folder, filename, plotting)
#print "Lambda = 800 nm"
#importFromTable(800e-9, folder, filename, plotting)
#print "Lambda = 795 nm"
#importFromTable(795e-9, folder, filename, plotting)
#plotting = False
#print "Lambda = 625 nm"
#importFromTable(625e-9, folder, filename, plotting)
#print "Lambda = 532 nm"
#importFromTable(532e-9, folder, filename, plotting)
#print "Lambda = 515 nm"
#importFromTable(515e-9, folder, filename, plotting)
#print "Lambda = 400 nm"
#importFromTable(400e-9, folder, filename, plotting)

#===========================================
#importFromNKtable(folder, filename)
#print "Lambda = 3000 nm"
#importFromNKtable(3000e-9, folder, filename)

#importFromEpsilonTable(wavelength, folder, filename, plotting=True)
#===========================================
#print "Lambda = 515 nm"
#importFromAbsorptionData(515e-9, folder, filename, plotting=True)
