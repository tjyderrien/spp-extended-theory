import math, cmath
from scipy.integrate import ode
import matplotlib.pyplot as plt

#some constants
T_2 = 1
F_0 = 1
d_0x = 1

#external functions
def f(t):
    return math.sin(t)
def S(t):
    return 1

#initial condition, etc...
t_start = 0.0
t_stop = 10.0
dt = .05
X_start = [1.0, 2.0, 3.0]

#contains the ODE system matrix
def ODE_matrix(t, X):
    pi = X[0]
    n_v = X[1]
    n_c = X[2]
    pi_dot = -pi/T_2 - 1j*F_0*f(t)*d_0x*(n_v - n_c)*cmath.exp(-1j*S(t))
    n_v_dot = 2*F_0*f(t)*d_0x*(pi*cmath.exp(1j*S(t))).imag
    n_c_dot = -n_v_dot
    return [pi_dot, n_v_dot, n_c_dot]

#initializing the integrator
r = ode(ODE_matrix).set_integrator('zvode', method='bdf')
r.set_initial_value(X_start, t_start)

#integrate the problem to get results
while r.successful() and r.t < t_stop:
    r.integrate(r.t + dt)
    print(r.t, r.y)
    plt.scatter(r.t, r.y[0].real, c='r')
    plt.scatter(r.t, r.y[0].imag, c='g')
    plt.scatter(r.t, r.y[1].real, c='b')
    plt.scatter(r.t, r.y[1].imag, c='k')
    plt.scatter(r.t, r.y[2].real, c='cyan')
    plt.scatter(r.t, r.y[2].imag, c='magenta')

plt.show()
