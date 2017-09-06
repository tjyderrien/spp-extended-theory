import math, cmath
from scipy.integrate import ode
import matplotlib.pyplot as plt

T_2 = 1
F_0 = 1
d_0x = 1

def f(t):
    return math.sin(t)
def S(t):
    return t

t_start = 0.0
t_stop = 10.0
dt = .05
X_start = [1.0, 2.0j, -1.0j]

def odeRHS(t, X):
    pi = X[0]
    nv = X[1]
    nc = X[2]
    pi_dot = -pi/T_2 - 1j*F_0*f(t)*d_0x*(nv - nc)*cmath.exp(-1j*S(t))
    nv_dot = 2*F_0*f(t)*d_0x*(pi*cmath.exp(1j*S(t))).imag
    nc_dot = -nv_dot
    return [pi_dot, nv_dot, nc_dot]

r = ode(odeRHS).set_integrator('zvode', method='bdf')
r.set_initial_value(X_start, t_start)

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
