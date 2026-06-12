#!/usr/bin/env python2
#-*- coding: utf-8 -*-

import matplotlib.pyplot as plt
import numpy as np
from scipy.constants import e

from libKeldysh import *

Efield=np.linspace(1,300E8, 100)

# Parameters to repeat paper of Kaiser, Phys. Rev. B 2000. 
Egap = 9e0*e
meff = 1E0
wavelength=500E-9

W_tunnel = KeldyshTunnelingLimit(Egap, meff, wavelength, Efield)

plt.figure()
plt.semilogy(Efield, np.clip(W_tunnel,1E30, 1E45))
plt.xlabel('Electric field amplitude (V/m)')
plt.ylabel(r'$w_{TI}$')
plt.show()
