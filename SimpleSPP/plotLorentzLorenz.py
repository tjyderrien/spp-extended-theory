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

## @package LorentzLorenz
# Computes the dielectric permittivity as function of fraction of one material at a single wavelength. 
# Computes the SPP properties as function of the materials fraction

# IMPORT LIBRARIES
from libSPP import *
from libMaterials import *
from libMultilayerSPP import *

Header="[plotLorentzLorenz] "

wavelength = 1026e-9
epsCr2O3   = 3.8273816+ 0.0483803j
epsCr      = -0.6721223+24.8657476j
epsAg = -48.121511807025215 + 3.1028275101335927
epsAu = -49.150229501143876 + 3.7710634965016174j

eps1 = epsAg
eps2 = epsAu #epsCr

MaterialFile1 = "Ag"
#MaterialFile2 = r"Cr$_2$O$_3$"
MaterialFile2 = "Au"

fraction = np.arange(0., 1., 1e-2)
EffectivePermittivity = MaxwellGarnett2(eps1, eps2, 1.-fraction)
Reflectivity = reflectivity(1., EffectivePermittivity, 0E0, 'S')

print "Plotting Maxwell-Garnett 2-material mixing."
plt.figure()
plt.title('Mixture: '+MaterialFile1+'/'+MaterialFile2+r" ($\lambda=$"+str(int(1E9*wavelength))+" nm)")

ax1 = plt.subplot(111)
          
ax1.set_xlabel(r'Fraction of '+MaterialFile1+' ($\%$)')
plot11, = ax1.plot(fraction*100., EffectivePermittivity.real, 'b-', label=r'Re$(\varepsilon_{eff})$')
plot12, = ax1.plot(fraction*100., EffectivePermittivity.imag, 'b--', label=r'Im$(\varepsilon_{eff})$')

ax12 = ax1.twinx()
plot2, = ax12.plot(fraction*100., Reflectivity, 'k-', label=r'Re$(\varepsilon_{eff})$')

plotComb1 = [plot11, plot12, plot2]

labelsComb1 = [l.get_label() for l in plotComb1]
ax1.legend(plotComb1, labelsComb1, loc='best')

#plt.legend(loc='best')
#plt.xlim((200.,2000.))
filename=MaterialFile1+"-"+MaterialFile2+'-wavelength'+str(int(wavelength))
plt.grid()
plt.tight_layout()
plt.savefig(filename+'.eps')
plt.savefig(filename+'.png')
plt.show()

#print Header+"** Preparing multilayer SPP computation..."




