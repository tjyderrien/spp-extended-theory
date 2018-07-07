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

# IMPORT LIBRARIES
from libSPP import *

#folder = "Database/"

Material1 = 'Air' #'ZnO (Bond 1965, o)'
Material2 = 'Si (Palik)' #'Au (Palik)'
wavelength = 800e-9
unit = 1E9

# Select database
database="MaterialOpticalDatabaseForPlasmonics.csv"

# Build database array for choosing which material can be of interest to irradiate
dbarray = loadtxt(database, dtype='str', delimiter='\t')

# Select the material of interface 1
DataMaterial1 = FilterDatabase(dbarray, Material1, 0)
DataMaterial2 = FilterDatabase(dbarray, Material2, 0)
#print DataMaterial1, DataMaterial2

# Select the wavelength
DataMaterial1_filtered = FilterDatabase(DataMaterial1, str(int(wavelength*unit)), 2)
DataMaterial2_filtered = FilterDatabase(DataMaterial2, str(int(wavelength*unit)), 2)
#print DataMaterial1_filtered, DataMaterial2_filtered


##print DataMaterial2[:,0:4]
Material1loc, BandGap, wavelengthLoc, RealEps, ImagEps = ExtractMaterialData(DataMaterial1_filtered)
RealEps = np.asfarray(RealEps); ImagEps = np.asfarray(ImagEps)
eps1 = RealEps + 1j * ImagEps
del RealEps, ImagEps

Material2loc, BandGap, wavelengthLoc, RealEps, ImagEps = ExtractMaterialData(DataMaterial2_filtered)
RealEps = np.asfarray(RealEps); ImagEps = np.asfarray(ImagEps)
eps2 = RealEps + 1j * ImagEps

print "** Plotting SPP condition as function of excitation"

LogNe = np.linspace(16.,np.log10(4*5E28),500)
Ne    = np.power(10., LogNe)
CollisionRate = 1.1E-15**-1
OpticalMass   = 0.18e0
print "** Calculating SPP condition with excitation level..."
SPPcondition = SPPconditionValue(eps1, Drude(wavelength, Ne, eps2, CollisionRate, OpticalMass))

epsAir = eps1
epsSi  = Drude(wavelength, Ne, eps2, CollisionRate, OpticalMass)

plt.figure()
plt.xlabel(r'$N_{e-h}$ (m$^{-3}$)')
plt.ylabel(r'$F_{spp}$')
plt.loglog(Ne, SPPcondition, 'r', label='SPP forbidden')
plt.loglog(Ne, -SPPcondition, 'b', label='SPP allowed')
plt.title('SPP condition '+str(Material1loc[0])+'/'+str(Material2loc[0]))
plt.legend(loc='best')
plt.tight_layout()
plt.savefig('Dispersion.eps')
#plt.show()

print "** Plotting the L_spp as function of excitation"

betaSPP        = np.vectorize(betaSPP)
DecayLengthSPP = np.vectorize(DecayLengthSPP)

print epsAir.shape, epsSi.shape
beta = betaSPP(wavelength, epsAir, epsSi)
Lspp = DecayLengthSPP(beta)

print np.shape(Lspp)

plt.figure()
plt.ylabel(r'$\beta$ (m$^{-1}$)')
plt.semilogx(Ne,  beta.real, label = r'$\beta^{+}$')
#plt.plot(Ne, -beta.real, label = r'$\beta^{-}$')
plt.legend(loc='best')
plt.tight_layout()
plt.savefig('Beta-ExcitationSi.eps')
#plt.show()

plt.figure()
plt.ylabel(r'$L_{SPP}$ (m)')
plt.loglog(Ne, Lspp,  label = r'$L_{SPP}^{+}$')
#plt.plot(Ne, -Lspp, label = r'$L_{SPP}^{-}$')
plt.legend(loc='best')
plt.tight_layout()
plt.grid()
plt.savefig('Lspp-ExcitationSi.eps')
plt.show()
