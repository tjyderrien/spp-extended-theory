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

# IMPORT LIBRARIES
from libSPP import *

#folder = "Database/"

Material1 = 'ZnO (Bond 1965, o)'
Material2 = 'Au (Palik)'
wavelength = 3000e-9
unit = 1E9

# Select database
database="MaterialOpticalDatabaseForPlasmonics.csv"

# Build database array for choosing which material can be of interest to irradiate
dbarray = loadtxt(database, dtype='str', delimiter='\t')

# Select the material of interface 1
DataMaterial1 = FilterDatabase(dbarray, Material1, 0)
DataMaterial2 = FilterDatabase(dbarray, Material2, 0)

# Select the wavelength
DataMaterial1 = FilterDatabase(DataMaterial1, wavelength*unit, 2)
DataMaterial2 = FilterDatabase(DataMaterial2, wavelength*unit, 2)
#print DataMaterial1, DataMaterial2
#print DataMaterial2[:,0:4]
Material1loc, BandGap, wavelengthLoc, RealEps, ImagEps = ExtractMaterialData(DataMaterial1)

#print Material1loc, BandGap, wavelengthLoc, RealEps, ImagEps 

#print "Calculating SPP condition with excitation level..."
#SPPconditionValue(eps1, Drude(eps2, Ne, CollisionRate, OpticalMass))
## plot SPPconditionNew(Drude(Material, CarrierDensity, CollisionRate, OpticalMass)

#plt.figure()
#plt.xlabel('N_{e-h} $(m^{-3})$')
#plt.ylabel('F_{spp}')
#plt.plot(Ne, SPPcondition, label='Ext-SPP cond.')
#plt.title('SPP condition $'+MaterialFile1+'$/$'+MaterialFile2+'$ interface')
#plt.legend(loc=2)
#plt.savefig('Dispersion.eps')
##plt.show()
