import cmath
import numpy as np

import matplotlib.pyplot as plt
from scipy.optimize import root, fsolve
from scipy.linalg import norm

eps1 = -1+1j
eps2 = 2-3j
eps3 = 1-7j
t = 1
k0 = 2

def disp(betaR, denom):
    betaC = betaR[0] + 1j * betaR[1]
    k1 = cmath.sqrt(betaC*betaC - k0*k0*eps1)
    k2 = cmath.sqrt(betaC*betaC - k0*k0*eps2)
    k3 = cmath.sqrt(betaC*betaC - k0*k0*eps3)

    if denom:
        out = (k1/eps1 - k2/eps2) * (k1/eps1 - k3/eps3)
    else:
        out = (k1/eps1 - k2/eps2) * (k1/eps1 - k3/eps3) * cmath.exp(-2*k1*t) - (k1/eps1 + k2/eps2) * (k1/eps1 + k3/eps3)      
    return [out.real, out.imag]

def func(w):
    z = w[0] + 1j * w[1]
    f  = z**10 + 1
    return [f.real, f.imag]

roots = []

for x in np.linspace(-1.5, 1.5, num=40):
    for y in np.linspace(-120, 120, num=2400):
        nrt = root(disp, [x, y], args=(False,))
        if nrt.success:
            nrt2 = root(disp, nrt.x, args=(False,))
            if nrt2.success:
                new = True
                if norm(disp(nrt2.x, True)) < 1e-10:
                    new = False
                for rt in roots:
                    if norm(rt - nrt2.x) < 1e-5:
                        new = False
                if new:
                    roots.append(nrt2.x)
                    print('%+.9f' % nrt2.x[0], '%+.9fi' % nrt2.x[1], '      ', '%+.9f' % norm(disp(nrt2.x, False)), '%+.9f' % norm(disp(nrt2.x, True)))



print(len(roots))

plt.scatter(*zip(*roots))

plt.show()
##
##i = 0
##data = []
##for x in np.linspace(-100, 100, num=20):
##    for y in np.linspace(-100, 100, num=20):
##        i += 1
##        F = root(disp, [x, y])
##        print(x, y, F)
##        data.append((x, F.x[0]))
##        data.append((y, F.x[1]))
##
##        if (i % 5 == 0):
##            data.append('r')
##        elif (i % 5 == 1):
##            data.append('g')
##        elif (i % 5 == 2):
##            data.append('b')
##        elif (i % 5 == 3):
##            data.append('m')
##        else:
##            data.append('black')
##
##plt.plot(*data)
##
##ax = plt.gca()
##ax.set_aspect('equal')
##plt.show()
