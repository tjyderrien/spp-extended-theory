#!/usr/bin/env python2
#-*- coding: utf-8 -*-

# Copyright (C) 2017 F. Preucil, T.J.-Y. Derrien
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

import math, cmath
from scipy.integrate import ode

#physical constants
atomic_dists = (5.32, 6.14, 9.83)                                    #atomic distances (units = bohrs)
band_coefs = ((0.1826, -0.1519, -0.0224, -0.0036, -0.0032, -0.0015), #coefficients (units = hartrees)
              (0.1454, -0.1454, 0, 0, 0, 0),
              (0.0494, -0.0494, 0, 0, 0, 0))

E_g = .1213 #gap energy (units = hartrees)
d_0x = 3.46 #x-component of the dipole moment (units = e*bohr)

#dispersion curve
def Epsilon(k):
    total = E_g
    for i in range(0, 3):       #iterating over directions
        for j in range(0, 6):   #cosine sum
            total += (band_coefs[i][j])*math.cos(j*k[i]*atomic_dists[i])/3
    return total

#we only have 1D field
#x-component of the laser field (squared sine converted to cosine)
def F_x(t):
    if t1 < t < t2:
        t_shifted = t-t0
        return peak2*(1+math.cos(piFWHM*t_shifted))*math.cos(omega*t_shifted+phi)
    else:
        return 0

#x-component of the vector potential (analytic integration)
def A_x(t):
    if t >= t2:
        return cnstvalA
    elif t1 < t < t2:
        t_shifted = t-t0
        return (c1*math.sin((omega-piFWHM)*t_shifted+phi)
                + c2*math.sin((omega+piFWHM)*t_shifted+phi)
                + c3*math.sin(omega*t_shifted+phi)
                + c4)
    else:
        return 0

#contains the ODE system RHS
def RHS(t, X, K, T_2):
    pi = X[0]
    n_v = X[1]
    n_c = X[2]
    S = X[3]
    pi_dot = -pi/T_2 -1j*d_0x*F_x(t)*(n_v-n_c)*cmath.exp(-1j*S)
    n_v_dot = 2*d_0x*F_x(t)*(pi*cmath.exp(1j*S)).imag
    n_c_dot = -n_v_dot
    S_dot = Epsilon([K[0]+A_x(t), K[1], K[2]])
    return [pi_dot, n_v_dot, n_c_dot, S_dot]

#main function
#https://docs.scipy.org/doc/scipy/reference/generated/scipy.integrate.ode.html

def solve(t, field_params, t_init, init_conds, K, T_2):
    global peak, omega, FWHM, t0, phi
    peak, omega, FWHM, t0, phi = field_params

    #precaching field values
    global pi2, t1, t2, FWHMomg, piFWHM, FWHM2mpi2omg, c1, c2, c3, c4, peak2, cnstvalA
    pi2 = math.pi**2
    t1 = t0-FWHM
    t2 = t0+FWHM

    FWHMomg = FWHM*omega
    piFWHM = math.pi/FWHM
    FWHM2mpi2omg = (FWHMomg**2-pi2)*omega

    c1 = -peak*FWHM/(4*(FWHMomg-math.pi))
    c2 = -peak*FWHM/(4*(FWHMomg+math.pi))
    c3 = -peak/(2*omega)
    c4 = peak*pi2/(2*FWHM2mpi2omg)*math.sin(FWHMomg-phi)

    peak2 = peak/2
    cnstvalA = peak*pi2*math.sin(FWHMomg)*math.cos(phi)/FWHM2mpi2omg
    #end of precaching

    intg = ode(RHS).set_integrator('zvode', method='bdf', nsteps=50000)
    intg.set_initial_value(init_conds, t_init).set_f_params(K, T_2)
    intg.integrate(t)
    return intg.y[1].real, intg.y[2].real
