import math, cmath
from scipy.integrate import ode
import matplotlib.pyplot as plt

#some constants
F_0 = 1
d_0x = 1
omega_0 = 1
t_0 = 5
T_2 = 1
log2 = math.log(2)

#external functions
def Omega(t):
    return F_0*d_0x*math.cos(omega_0*t)*math.exp(-log2*((t - t_0)/(5*omega_0))**2)
def S(t):
    return t

#initial condition, etc...
t_start = 0
t_stop = 10
dt = .1
X_start = [2, 2, 1]

#contains the ODE system matrix
def ODE_matrix(t, X):
    pi = X[0]
    n_v = X[1]
    n_c = X[2]
    pi_dot = -pi/T_2 - 1j*Omega(t)*(n_v - n_c)*cmath.exp(-1j*S(t))
    n_v_dot = 2*Omega(t)*(pi*cmath.exp(1j*S(t))).imag
    n_c_dot = -n_v_dot
    return [pi_dot, n_v_dot, n_c_dot]

#initializes the integrator
r = ode(ODE_matrix).set_integrator('zvode', method='bdf')
r.set_initial_value(X_start, t_start)

#integrates the problem to get results
times = []
pi_Re = []
pi_Im = []
n_v = []
n_c = []
while r.successful() and r.t < t_stop:
    print(r.t, r.y)
    times.append(r.t)
    pi_Re.append(r.y[0].real)
    pi_Im.append(r.y[0].imag)
    n_v.append(r.y[1].real)
    n_c.append(r.y[2].real)
    r.integrate(r.t + dt)

plt.plot(times, pi_Re, 'r')
plt.plot(times, pi_Im, 'r', linestyle='--')
plt.plot(times, n_v, 'g')
plt.plot(times, n_c, 'b')
plt.show()
