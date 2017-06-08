import cmath
import numpy as np

import matplotlib.pyplot as plt
from scipy.optimize import root
from scipy.linalg import norm
from itertools import product
from libMaterials import Drude

pi = np.pi

#data
wavelength = 800e-9
ne=1E28 #(m^-3) quantity of electrons in conduction band
nu = (1.1e-15)**-1 #collision time between conduction band electrons
meff = 0.18

epsilon=13.64+0.048j
eps2 = 1.+0.j       #environment
eps3 = 13.64+0.048j #substrate Si (no excitation)
eps1 = Drude(wavelength, ne, eps3, nu, meff)

k0 = 2*pi/wavelength
t = 10e-9 #in meters

#area of initial guesses

#x = real part, y = imaginary part
x_min = -1e8
x_max = 1e8
x_steps = 20

y_min = -1e9
y_max = 1e9
y_steps = 20

#tolerances

t_zero = 100000
t_smear = 10

def func(betaR):
    beta = betaR[0] + betaR[1]*1j
    k1 = sgn1*cmath.sqrt(beta*beta - k0*k0*eps1)
    k2 = sgn2*cmath.sqrt(beta*beta - k0*k0*eps2)
    k3 = sgn3*cmath.sqrt(beta*beta - k0*k0*eps3)

    try:
        out = (k1/eps1 - k2/eps2) * (k1/eps1 - k3/eps3) * cmath.exp(-2*k1*t) - (k1/eps1 + k2/eps2) * (k1/eps1 + k3/eps3)
    except:
        out = 1e99
    return [out.real, out.imag]

print('Guess area is a rectangle: [%s, %s]x[%si, %si], x_steps = %d, y_steps = %d' % (x_min, x_max, y_min, y_max, x_steps, y_steps))
print('Initial data:')
print('    Eps1:', eps1)
print('    Eps2:', eps2)
print('    Eps3:', eps3)
print('    k0:', k0)
print('    t:', t)

colors = ['red', 'green', 'blue', 'black', 'magenta', 'cyan', 'lime', 'orangered']
data = []
i = 0
for sgn1, sgn2, sgn3 in product((-1,1), (-1,1), (-1,1)):
    for x in np.linspace(x_min, x_max, num=x_steps):
        for y in np.linspace(y_min, y_max, num=y_steps):
            nrt = root(func, [x, y], method='hybr')
            if nrt.success and norm(nrt.fun) < t_zero:
                data.append((x, nrt.x[0]))
                data.append((y, nrt.x[1]))
                data.append(colors[i])
    i += 1

plt.plot(*data)

plt.xlabel('Re')
plt.ylabel('Im')
plt.grid()
#plt.savefig('betaSolution.eps')
plt.show()
