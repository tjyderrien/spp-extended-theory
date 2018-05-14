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

# IMPORT LIBRARIES
from libSPP import *
from libMaterials import *

wavelength = 1026e-9
epsCr2O3   = 3.8273816+ 0.0483803j
epsCr      = -0.6721223+24.8657476j

eps1 = epsCr
eps2 = epsCr2O3

MaterialFile1 = "Cr"
MaterialFile2 = r"Cr$_2$O$_3$"

fraction = np.arange(0., 1., 1e-2)
EffectivePermittivity = MaxwellGarnett2(eps1, eps2, 1.-fraction)
print "Plotting Maxwell-Garnett 2-material mixing."
plt.figure()
plt.title('Mixture: '+MaterialFile1+'/'+MaterialFile2+r" ($\lambda=$"+str(int(1E9*wavelength))+" nm)")
plt.xlabel(r'Fraction')
plt.plot(fraction*100., EffectivePermittivity.real, 'b-', label=r'$Re(\varepsilon_{eff})$')
plt.plot(fraction*100., EffectivePermittivity.imag, 'r-', label=r'$Im(\varepsilon_{eff})$')
plt.legend(loc=1)
#plt.xlim((200.,2000.))
filename=MaterialFile1+"-"+MaterialFile2+'-wavelength'+str(int(wavelength))
plt.grid()
plt.savefig(filename+'.eps')
plt.savefig(filename+'.png')
plt.show()
