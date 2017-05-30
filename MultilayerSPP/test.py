import cmath
import numpy as np

import matplotlib.pyplot as plt
from scipy.optimize import root, fsolve

def func(w):
    z = w[0] + 1j*w[1]
    f  = z**5 + 1
    return [f.real, f.imag]

i = 0
data = []
for x in np.linspace(-1.2, 1.2, num=20):
    for y in np.linspace(-1.2, 1.2, num=20):
        i += 1
        F = root(func, [x, y])
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
