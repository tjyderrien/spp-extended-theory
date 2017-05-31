import cmath
import numpy as np

import matplotlib.pyplot as plt
from scipy.optimize import root, fsolve


eps1 = 1+1j
eps2 = -2+1j
eps3 = 1+3j
t = 1
k0 = 1

def disp(beta):
    b = beta[0] + 1j * beta[1]
    k1 = cmath.sqrt(b * b - k0 * k0 * eps1)
    k2 = cmath.sqrt(b * b - k0 * k0 * eps2)
    k3 = cmath.sqrt(b * b - k0 * k0 * eps3)

    eq = (k1/eps1 - k2/eps2)*(k1/eps1 - k3/eps3)*cmath.exp(-2*k1*t) - (k1/eps1 + k2/eps2)*(k1/eps1 + k3/eps3)
    return [eq.real, eq.imag];

def func(w):
    z = w[0] + 1j * w[1]
    f  = -cmath.exp(-2*z) + z
    return [f.real, f.imag]

i = 0
data = []
for x in np.linspace(-10, 10, num=20):
    for y in np.linspace(-10, 10, num=20):
        i += 1
        F = root(disp, [x, y])
        print(x, y, F)
        data.append((x, F.x[0]))
        data.append((y, F.x[1]))

        if (i % 5 == 0):
            data.append('r')
        elif (i % 5 == 1):
            data.append('g')
        elif (i % 5 == 2):
            data.append('b')
        elif (i % 5 == 3):
            data.append('m')
        else:
            data.append('black')

plt.plot(*data)

ax = plt.gca()
ax.set_aspect('equal')
plt.show()
