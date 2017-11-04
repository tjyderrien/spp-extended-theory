import math, cmath
import numpy as np
from scipy.integrate import ode
import matplotlib.pyplot as plt

#physical constants
atomic_dist = (5.32, 6.14, 9.83)                                    #atomic distances
band_amps_v = ((-0.0928, 0.0705, 0.0200, -0.0012, 0.0029, 0.0006),  #band amplitudes
               (-0.0307, 0.0307, 0, 0, 0, 0),
               (-0.0059, 0.0059, 0, 0, 0, 0))
band_amps_c = ((0.0898, -0.0814, -0.0024, -0.0048, -0.0003, -0.0009),
               (0.1147, -0.1147, 0, 0, 0, 0),
               (0.0435, -0.0435, 0, 0, 0, 0))
E_g = 3.3 #gap energy
d_0x = 3.46 #x-component of the dipole moment
T_2 = 10 #damping factor

#laser field parameters
peak = 1
omega = 20
FWHM = 1
t0 = 1
phi = 0

#precaching values
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

#dispersion curve
def Epsilon(k):
    total = E_g
    for i in range(0, 3):       #iterating over directions
        for j in range(0, 6):   #cosine sum
            total += (band_amps_c[i][j]-band_amps_v[i][j])*math.cos(j*k[i]*atomic_dist[i])/3
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

#--------------------------------------------------------------------------------------------------

K = [0, 0, 0]
#initial condition
t_init = 0
t_end = 3
samples = 300

pi_init = 2 + 1j
n_v_init = 3
n_c_init = 1
S_init = 0

print('Two bands model')
print('-'*15)
print()
print('E_g = %.3e' % E_g)
print('d_0x = %.3e' % d_0x)
print()
print('peak = %.3e' % peak)
print('omega = %.3e' % omega)
print('FWHM = %.3e' % FWHM)
print('t0 = %.3e' % t0)
print('phi = %.3e' % phi)
print()
print('T_2 = %.3e' % T_2)
print()
print('K =', K)
print('t_init = %.3e' % t_init)
print('t_end = %.3e' % t_end)
print('samples = %d' % samples)
print()
print('pi_init = %.3e + %.3ej' % (pi_init.real, pi_init.imag))
print('n_v_init = %.3e' % n_v_init)
print('n_c_init = %.3e' % n_c_init)
print('S_init = %.3e' % S_init)
print()

#contains the ODE system RHS
def RHS(t, X):
    pi = X[0]
    n_v = X[1]
    n_c = X[2]
    S = X[3]
    pi_dot = -pi/T_2 -1j*d_0x*F_x(t)*(n_v-n_c)*cmath.exp(-1j*S)
    n_v_dot = 2*d_0x*F_x(t)*(pi*cmath.exp(1j*S)).imag
    n_c_dot = -n_v_dot
    S_dot = Epsilon([K[0]+A_x(t), K[1], K[2]])
    return [pi_dot, n_v_dot, n_c_dot, S_dot]

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
F_xs = []
A_xs = []

tlen = len(times)
ind = 0
for t in times:
    ind += 1
    intg.set_initial_value([pi_init, n_v_init, n_c_init, S_init], t_init)
    intg.integrate(t)
    pi_Re.append(intg.y[0].real)
    pi_Im.append(intg.y[0].imag)
    n_v.append(intg.y[1].real)
    n_c.append(intg.y[2].real)
    Ss.append(intg.y[3].real)
    F_xs.append(F_x(t))
    A_xs.append(A_x(t))
    print('Progress: %d%%' % (ind/tlen * 100), end = '\r')

plt.plot(times, pi_Re, 'r', label='Re(pi)')
plt.plot(times, pi_Im, 'r', linestyle='--', label='Im(pi)')
plt.plot(times, n_v, 'g', label='n_v')
plt.plot(times, n_c, 'b', label='n_c')
#plt.plot(times, Ss, 'c', label='S')
plt.plot(times, F_xs, 'k', label='F_x')
plt.plot(times, A_xs, 'gray', label='A_x')
plt.legend(loc='upper right')
plt.grid(True)
plt.show()
