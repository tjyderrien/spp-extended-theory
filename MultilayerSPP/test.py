import cmath
import numpy as np

import matplotlib.pyplot as plt
from scipy.optimize import root
from scipy.linalg import norm
from itertools import product
from libMaterials import Drude

pi = 2.*np.arcsin(1.0)

#data
wavelength = 800e-9
ne=5E27 #(m^-3) quantity of electrons in conduction band
nu = (1.1E-15)**-1 #collision time between conduction band electrons
meff = 0.18

epsilon=13.64+0.048j
eps2 = 1.+0.j       #environment
eps3 = 13.64+0.048j #substrate Si (no excitation)
eps1 = Drude(wavelength, ne, eps3, nu, meff)

k0 = 2*pi/wavelength
t = 10e-9 #in meters

#guess area
step = .01

x_min = -1E9 #real part of beta
x_max = 1E9
x_steps = 60  #round((x_max - x_min)/step)

y_min = -1E9
y_max = 1E9
y_steps = 300#round((y_max - y_min)/step)

#tolerances
t_sim = 1e-1

def func(betaR):
    beta = betaR[0] + betaR[1]*1j
    k1 = sgn1*cmath.sqrt(beta*beta - k0*k0*eps1)
    k2 = sgn2*cmath.sqrt(beta*beta - k0*k0*eps2)
    k3 = sgn3*cmath.sqrt(beta*beta - k0*k0*eps3)

    try:
        out = (k1/eps1 - k2/eps2) * (k1/eps1 - k3/eps3) * cmath.exp(-2*k1*t) - (k1/eps1 + k2/eps2) * (k1/eps1 + k3/eps3)
    except:
        out = 1e20
    return [out.real, out.imag]

##    try:
##        outj = - (k1/eps1 + k3/eps3)*(beta/(k1*eps1) + beta/(k2*eps2)) - (k1/eps1 + k2/eps2)*(beta/(k1*eps1) + beta/(k3*eps3))
##        + cmath.exp(-2*k1*t)*(
##            (k1/eps1 - k3/eps3)*(beta/(k1*eps1) - beta/(k2*eps2))
##            + (k1/eps1 - k2/eps2)*(beta/(k1*eps1) - beta/(k3*eps3))
##            - 2*beta*t*(k1/eps1 - k2/eps2)*(k1/eps1 - k3/eps3)/k1
##            )
##        
##    except:
##        outj = 1e20

     #, np.array([[outj.real, (1j*outj).real], [outj.imag, (1j*outj).imag]])]

print('Guess area is a rectangle: [%f, %f]x[%fi, %fi], x_steps = %d, y_steps = %d' % (x_min, x_max, y_min, y_max, x_steps, y_steps))
print()
mroots = []
for sgn1, sgn2, sgn3 in product((-1,1), (-1,1), (-1,1)):
    roots = []
    out = []
    print('Branch: [%s][%s][%s]' % (('%+d' % sgn1)[0], ('%+d' % sgn2)[0], ('%+d' % sgn3)[0]))
    for x in np.linspace(x_min, x_max, num=x_steps):
        for y in np.linspace(y_min, y_max, num=y_steps):
            nrt = root(func, [x, y], method='hybr', tol=1e-14)
            if nrt.success:
                roots.append(nrt.x)

    while len(roots) > 0:
        cntr = roots[0]
        same = [cntr]
        nrest = []

        for rt in roots[1:]:
            if norm(rt - cntr, 1) < t_sim:
                same.append(rt)
                cntr = np.mean(same, axis = 0)
            else:
                nrest.append(rt)

        out.append(cntr)
        roots = list(nrest)

    for rt in out:
        print('%+2.15f %+2.15f    %+2.15f %+2.15f' % (rt[0], rt[1], func(rt)[0], func(rt)[1]))
    print('Total: %d' % len(out))
    print()
    mroots.append(out)

##a = cmath.sqrt(eps1*k0*k0)
##plt.plot((a.real, -a.real), (a.imag, -a.imag), 'r')
plt.scatter(*zip(*mroots[0]), c='red')
plt.scatter(*zip(*mroots[1]), c='green')
plt.scatter(*zip(*mroots[2]), c='blue')
plt.scatter(*zip(*mroots[3]), c='black')
plt.scatter(*zip(*mroots[4]), c='magenta')
plt.scatter(*zip(*mroots[5]), c='cyan')
plt.scatter(*zip(*mroots[6]), c='lime')
plt.scatter(*zip(*mroots[7]), c='orangered')
plt.xlabel('Re')
plt.ylabel('Im')
plt.grid()
plt.savefig('betaSolution.eps')
plt.show()
