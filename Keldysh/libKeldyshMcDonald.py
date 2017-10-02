import math, cmath
import numpy as np
from scipy.integrate import ode
from scipy.integrate import quad
from scipy.special import erf #for analytic
from scipy.special import wofz
import matplotlib.pyplot as plt

#some constants
log2 = math.log(2)
slog2 = math.sqrt(math.log(2))
spilog2 = math.sqrt(math.pi/math.log(2))
atomic_dist = (5.32, 6.14, 9.83)                                    #atomic distances
band_amps_v = ((-0.0928, 0.0705, 0.0200, -0.0012, 0.0029, 0.0006),  #band amplitudes
               (-0.0307, 0.0307, 0, 0, 0, 0),
               (-0.0059, 0.0059, 0, 0, 0, 0))
band_amps_c = ((0.0898, -0.0814, -0.0024, -0.0048, -0.0003, -0.0009),
               (0.1147, -0.1147, 0, 0, 0, 0),
               (0.0435, -0.0435, 0, 0, 0, 0))
E_g = 0 #gap energy
d_0x = 3.46 #x-component of the dipole moment

F_0 = 1 #laser field parameters
omega_0 = .28
t_0 = 5

T_2 = 10 #damping factor

#some precached numbers for analytic vector potential
prefactor = -2.5*F_0*omega_0*spilog2
a_anl = 1j*slog2/(5.*omega_0)
b_anl = (-1j*t_0*2*log2-25*omega_0*omega_0*omega_0)/(10*omega_0*slog2)
shift = 2*math.exp(-25*omega_0*omega_0*omega_0*omega_0/(4*log2))*math.cos(t_0*omega_0)

#dispersion curve
def Epsilon(k):
    total = E_g
    for i in range(0, 3):       #iterating over directions
        for j in range(0, 6):   #cosine sum
            total += (band_amps_c[i][j]-band_amps_v[i][j])*math.cos(j*k[i]*atomic_dist[i])/3
    return total

#we only have 1D field
#x-component of the laser field
def F_x(t):
    return F_0*math.cos(omega_0*t)*math.exp(-log2*(t-t_0)*(t-t_0)/(25*omega_0*omega_0))
#x-component of the vector potential, temporal integral of the laser field
def A_x(t):
    #return analytic_A_x(t)
    return -quad(F_x, -np.inf, t, limit=500)[0]

#analytic vector potential
def analytic_A_x(t):
    aux = (wofz(a_anl*t+b_anl)*cmath.exp(-1j*omega_0*t)).real
    return prefactor*(shift - math.exp(-log2*(t-t_0)*(t-t_0)/(25*omega_0*omega_0))*aux)

#action integral
def S(K, t):
    return quad(lambda t, K: Epsilon([K[0]+A_x(t), K[1], K[2]]), -np.inf, t, args=K, limit=500)[0]

#--------------------------------------------------------------------------------------------------

K = [0, 0, 0]
#initial condition
t_init = 0
t_end = 30
samples = 500

pi_init = 5 + 1j
n_v_init = 3
n_c_init = 0
A_x_init = A_x(t_init)
S_init = S(K, t_init)

print('t_init = %f' % t_init)
print('t_end = %f' % t_end)
print('A_x_init = %.15f' % A_x_init)
print('S_init = %.15f' % S_init)

#contains the ODE system RHS
def RHS(t, X):
    pi = X[0]
    n_v = X[1]
    n_c = X[2]
    A_x = X[3]
    S = X[3]
    pi_dot = -pi/T_2 -1j*d_0x*F_x(t)*(n_v-n_c)*cmath.exp(-1j*S)
    n_v_dot = 2*d_0x*F_x(t)*(pi*cmath.exp(1j*S)).imag
    n_c_dot = -n_v_dot
    A_x_dot = -F_x(t)
    S_dot = Epsilon([K[0]+A_x.real, K[1], K[2]])
    return [pi_dot, n_v_dot, n_c_dot, A_x_dot, S_dot]

#initializes the integrator
#https://docs.scipy.org/doc/scipy/reference/generated/scipy.integrate.ode.html
intg = ode(RHS).set_integrator('zvode', method='bdf', nsteps=5000)

#integrates the problem to get results
times = np.linspace(t_init, t_end, samples)
pi_Re = []
pi_Im = []
n_v = []
n_c = []
Ss = []
A_xs = []
F_xs = []
aA_xs = []

tlen = len(times)
ind = 0
for t in times:
    ind += 1
    intg.set_initial_value([pi_init, n_v_init, n_c_init, A_x_init, S_init], t_init)
    intg.integrate(t)
    pi_Re.append(intg.y[0].real)
    pi_Im.append(intg.y[0].imag)
    n_v.append(intg.y[1].real)
    n_c.append(intg.y[2].real)
    A_xs.append(intg.y[3].real)
    Ss.append(intg.y[4].real)
    F_xs.append(F_x(t))
    aA_xs.append(analytic_A_x(t))
    print('Progress: %d%%' % (ind/tlen * 100), end = '\r')

plt.plot(times, pi_Re, 'r', label='Re(pi)')
plt.plot(times, pi_Im, 'r', linestyle='--', label='Im(pi)')
plt.plot(times, n_v, 'g', label='n_v')
plt.plot(times, n_c, 'b', label='n_c')
plt.plot(times, A_xs, 'm', label='A_x')
plt.plot(times, Ss, 'c', label='S')
plt.plot(times, F_xs, 'gray', label='F_x')
#plt.plot(times, aA_xs, 'lime', label='analytic_A_x')
plt.legend(loc='upper right')
plt.show()
