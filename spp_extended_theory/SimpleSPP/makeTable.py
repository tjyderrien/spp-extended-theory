#!/usr/bin/env python
#-*- coding: utf-8 -*-
# Copyright (C) 2013-2025 T. J.-Y. Derrien
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

## @package makeTable
# Generates the table of SPP-active interfaces. 

# IMPORT LIBRARIES
from spp_extended_theory.Libs import libSPP #libSPP import *
from spp_extended_theory.Libs import libDatabase
import numpy as np
import pandas as pd
#from plotGraph import *

# PHYSICAL INPUT
#wavelength = 800e-9
#epsAir=1e0
#epsSi0=13.64+0.048j
#meffe=0.18
#nuSi=(1.1e-15)**-1

#example=period(betaSPP(wavelength, epsAir, Drude(wavelength, 1e28, epsSi0, nuSi)))

#print example

print("Generating SPP database...")
SPPdb = libSPP.GenerateDatabase()

SppOutput = 'SPPactiveInterfaces.dat'
print("Exporting to "+SppOutput+"...")

print()
SPPdb_export=np.vectorize(lambda x: float(x) if x.replace(".", "", 1).isdigit() else np.nan)(SPPdb)
libDatabase.ExportToTxt(SPPdb_export, SppOutput)

# Writing table caption
f=open(SppOutput, "a")

#Material1, Material2, Wavelength, OldSPPactiveBool, NewSPPactiveBool, #5
#Period, PeriodError, SPPdecayDepth1, SPPdecayDepth2, Reflectivity, #10
#OpticalPenetration1, OpticalPenetration2, SPPdecayLength, eps1.real, eps1.imag, #15
#eps2.real, eps2.imag, SPPdepthImagk1, SPPdepthImagk2, deltaLsppValues.real, #20
#Fa, Jo, LifeTime

f.write("#Material1\tMaterial2\tWavelength\tOldSPPactiveBool\tNewSPPactiveBool"
        "\tSPPperiod\tSPPperiodError\tSPPdecayDepth1\tSPPdecayDepth2\tReflectivity"
        "\tOpticalPenetration1\tOpticalPenetration2\tSPPdecayLength\teps1real\teps1imag"
        "\teps2real\teps2imag\tSPPdepthImagk1\tSPPdepthImagk2\tDeltaLsppValues"
        "\tFaradayNumber\tJouleNumber\tLifetimeRaether")
f.close()
print() 

# Writing table caption
# f=open(SppOutput, "a")

data=pd.DataFrame(SPPdb)
print("Exported. Please open file "+SppOutput+".")
data.to_csv(SppOutput)


deltaBetaSPP = np.vectorize(libSPP.deltaBetaSPP)
deltaPeriodSPP = np.vectorize(libSPP.deltaPeriodSPP)
deltaLspp = np.vectorize(libSPP.deltaLspp)

db = pd.DataFrame(SPPdb)
db.describe(include="all")
db.dtypes


# Writing table caption
print("Exported. Please open file "+SppOutput+".")


##  plot precision of Lambda over precision of epsilon

#eta = np.arange(1e-3, 1e0, 1e-3)

#wavelength = 1030e-9
#epsAir = 1+0j
#epsAu = -26.154188586+1.8503881331j
#epsTi = -4.2656+27.277j
#epsMo = -11.728+20.297j
#epsAl = -90.19+27.70j

#eps1 = epsAir
#eps2 = epsAl

#DeltaPeriod = deltaPeriodSPP(wavelength, eps1, eps2, eta, eta, eta, eta)
#SPPperiod = period(betaSPP(wavelength, eps1, eps2))
#DeltaPeriodRel = DeltaPeriod / SPPperiod

#print "Wavelength (nm)"
#print "Precision on period: "+str(1e9*DeltaPeriod)+" nm"
#print "Relative precision on period "+str(100*DeltaPeriodRel)+"%"

#fig = plt.figure()
#plt.title(r'Period uncertainty at $\lambda=$'+str(wavelength*1E9)+' nm')
#ax1 = fig.add_subplot(111)
#ax1.semilogx(eta, 1e2*DeltaPeriodRel, '-r', label='abs')
#ax1.set_xlabel(r'$\eta$')
#ax1.set_ylabel('Relative uncertainty (\%)', color='r')
#for tl in ax1.get_yticklabels():
    #tl.set_color('r')
#ax2 = ax1.twinx()
#ax2.semilogx(eta, 1e9*DeltaPeriod, '-b', label='rel')
#ax2.set_ylabel(r'$\delta \Lambda_{SPP}$ (nm)', color='b')
#for tl in ax2.get_yticklabels():
    #tl.set_color('b')
#plt.grid()
#plt.show()
#plt.savefig('PeriodError.eps')

#print 1e9*deltaPeriodSPP(800e-9, 1+0j, -26.154188586+1.8503881331j, 1e-4, 1e-4, 1e-4, 1e-4)
#print deltaBetaSPP(800e-9, 1e0+0e0j, -26.154188586e0+1j*1.8503881331e0, 1e-4, 1e-4, 1e-4, 1e-4)


""" TODO: interface with HTML for publication on the web. 
1. Put results into a NP.array.
2. Use a converter to HTML, CSV and PDF maybe. 
"""

# SPPactiveInterfacesArray contains all data we need, just remove lines starting with #. 

