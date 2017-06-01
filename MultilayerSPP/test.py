import cmath
import numpy as np

import matplotlib.pyplot as plt
from scipy.optimize import root, fsolve
from scipy.linalg import norm
from itertools import product

eps1 = 1 + 1j
eps2 = 3 + 1j
eps3 = 1 + 3j
k0 = 1
t = 2

def disp(betaR):
    betaC = betaR[0] + betaR[1]*1j
    k1 = sgn1*cmath.sqrt(betaC*betaC - k0*k0*eps1)
    k2 = sgn2*cmath.sqrt(betaC*betaC - k0*k0*eps2)
    k3 = sgn3*cmath.sqrt(betaC*betaC - k0*k0*eps3)

    out = (k1/eps1 - k2/eps2) * (k1/eps1 - k3/eps3) * cmath.exp(-2*k1*t) - (k1/eps1 + k2/eps2) * (k1/eps1 + k3/eps3)      
    return [out.real, out.imag]

mroots = []
for sgn1, sgn2, sgn3 in product((-1,1), (-1,1), (-1,1)):
    roots = []
    for x in np.linspace(-2, 2, num=500):
        for y in np.linspace(-10, 10, num=500):
            nrt = root(disp, [x, y], method='hybr')
            if nrt.success:
                new = True
                if norm(disp(nrt.x)) > .01:
                    new = False
                for rt in roots:
                    if norm(rt - nrt.x) < 1e-6:
                        new = False
                if new:
                    roots.append(nrt.x)
                    print('%+.9f' % nrt.x[0], '%+.9fi' % nrt.x[1], '      ', '%+.9f' % norm(disp(nrt.x)), root(disp, nrt.x, method='lm').success)
    print(len(roots))
    mroots.append(roots)

a = cmath.sqrt(eps1*k0*k0)
plt.plot((a.real, -a.real), (a.imag, -a.imag), 'r')
plt.scatter(*zip(*mroots[0]), c='red')
plt.scatter(*zip(*mroots[1]), c='green')
plt.scatter(*zip(*mroots[2]), c='blue')
plt.scatter(*zip(*mroots[3]), c='black')
plt.scatter(*zip(*mroots[4]), c='magenta')
plt.scatter(*zip(*mroots[5]), c='cyan')
plt.scatter(*zip(*mroots[6]), c='lime')
plt.scatter(*zip(*mroots[7]), c='orangered')

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
