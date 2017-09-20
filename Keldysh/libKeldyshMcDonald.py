import math, cmath
import numpy as np
from scipy.integrate import ode
from scipy.integrate import quad
import matplotlib.pyplot as plt

#some constants
log2 = math.log(2)
band_freq = (5.32, 6.14, 9.83)                                          #band frequencies
band_amps_v = ((-0.0928, 0.0705, 0.0200, -0.0012, 0.0029, 0.0006),      #band amplitudes
               (-0.0307, 0.0307, 0, 0, 0, 0),
               (-0.0059, 0.0059, 0, 0, 0, 0))
band_amps_c = ((0.0898, -0.0814, -0.0024, -0.0048, -0.0003, -0.0009),
               (0.1147, -0.1147, 0, 0, 0, 0),
               (0.0435, -0.0435, 0, 0, 0, 0))
E_g = 3.3                                                               #gap energy
d_0x = 3.46                                                             #x-component of the dipole moment

F_0 = 1                                                                 #laser field parameters
omega_0 = 2
t_0 = 20

T_2 = 1                                                                 #what is this?

#dispersion curve
def Epsilon(k):
    total = E_g
    for i in range(0, 3):       #iterating over directions
        total_in_dir = 0
        for j in range(0, 6):   #cosine sum
            total_in_dir += (band_amps_c[i][j]-band_amps_v[i][j])*math.cos(j*k[i]*band_freq[i])/3
        total += total_in_dir
    return total

#we only have 1D field
#x-component of the laser field
def F_x(t):
    return F_0*math.cos(omega_0*t)*math.exp(-log2*((t-t_0)/(5*omega_0))**2)
#x-component of the vector potential, temporal integral of the laser field
def A_x(t):
    return quad(F_x, -np.inf, t, limit=100)[0]

#action integral
def S(K, t):
    return quad(lambda t, K: Epsilon([K[0] + A_x(t), K[1], K[2]]), -np.inf, t, args=K, limit=100)[0]

#initial condition
t_start = -10
t_stop = 18
samples = 600

pi_start = 0
n_v_start = 2
n_c_start = 1
S_start = 0
A_start = 0

#contains the ODE system RHS
def RHS(t, X):
    pi = X[0]
    n_v = X[1]
    n_c = X[2]
    S = X[3]
    A = X[4]
    pi_dot = -pi/T_2-1j*d_0x*F_x(t)*(n_v-n_c)*cmath.exp(-1j*S)
    n_v_dot = 2*d_0x*F_x(t)*(pi*cmath.exp(1j*S)).imag
    n_c_dot = -n_v_dot
    S_dot = Epsilon([A.real, 0, 0])
    A_dot = F_x(t)
    return [pi_dot, n_v_dot, n_c_dot, S_dot, A_dot]

#initializes the integrator
#https://docs.scipy.org/doc/scipy/reference/generated/scipy.integrate.ode.html
r = ode(RHS).set_integrator('zvode', method='bdf', nsteps=5000)

#integrates the problem to get results
times = np.linspace(t_start, t_stop, samples)
pi_Re = []
pi_Im = []
n_v = []
n_c = []
Ss = []
As = []
for t in times:
    r.set_initial_value([pi_start, n_v_start, n_c_start, S_start, A_start], t_start)
    r.integrate(t)
    pi_Re.append(r.y[0].real)
    pi_Im.append(r.y[0].imag)
    n_v.append(r.y[1].real)
    n_c.append(r.y[2].real)
    Ss.append(r.y[3].real)
    As.append(r.y[4].real)

plt.plot(times, pi_Re, 'r', label='Re(pi)')
plt.plot(times, pi_Im, 'r', linestyle='--', label='Im(pi)')
plt.plot(times, n_v, 'g', label='n_v')
plt.plot(times, n_c, 'b', label='n_c')
#plt.plot(times, Ss, 'yellow', label='S')
plt.plot(times, As, 'magenta', label='A')
plt.legend(loc='upper right')
plt.show()
