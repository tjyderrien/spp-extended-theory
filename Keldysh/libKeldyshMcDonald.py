import math, cmath
import numpy as np
from scipy.integrate import ode
import matplotlib.pyplot as plt

#some constants
F_0 = 1
d_0x = 1
omega_0 = 1
t_0 = 4
T_2 = 1
log2 = math.log(2)

#external functions
def Omega(t):
    return F_0*d_0x*math.cos(omega_0*t)*math.exp(-log2*((t - t_0)/(5*omega_0))**2) #no K-dependency?
def S(t):
    return t #this is a placeholder, it should contain the action, which is K-dependent

#initial condition
t_start = 0
t_stop = 8
points = 200

pi_start = 3
n_v_start = 2
n_c_start = 1

#contains the ODE system RHS
def RHS(t, X):
    pi = X[0]
    n_v = X[1]
    n_c = X[2]
    pi_dot = -pi/T_2 - 1j*Omega(t)*(n_v - n_c)*cmath.exp(-1j*S(t))
    n_v_dot = 2*Omega(t)*(pi*cmath.exp(1j*S(t))).imag
    n_c_dot = -n_v_dot
    return [pi_dot, n_v_dot, n_c_dot]

#initializes the integrator
#https://docs.scipy.org/doc/scipy/reference/generated/scipy.integrate.ode.html
r = ode(RHS).set_integrator('zvode', method='bdf')

#integrates the problem to get results
times = np.linspace(t_start, t_stop, points)
pi_Re = []
pi_Im = []
n_v = []
n_c = []
for t in times:
    r.set_initial_value([pi_start, n_v_start, n_c_start], t_start)
    r.integrate(t)
    pi_Re.append(r.y[0].real)
    pi_Im.append(r.y[0].imag)
    n_v.append(r.y[1].real)
    n_c.append(r.y[2].real)

plt.plot(times, pi_Re, 'r', label='Re(pi)')
plt.plot(times, pi_Im, 'r', linestyle='--', label='Im(pi)')
plt.plot(times, n_v, 'g', label='n_v')
plt.plot(times, n_c, 'b', label='n_c')
plt.legend(loc='upper right')
plt.show()
