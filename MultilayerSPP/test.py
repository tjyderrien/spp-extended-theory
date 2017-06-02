import cmath
import numpy as np

import matplotlib.pyplot as plt
from scipy.optimize import root, fsolve
from scipy.linalg import norm
from itertools import product

#data
eps1 = 1 + 1j
eps2 = 3 + 1j
eps3 = 1 + 3j
k0 = 1
t = 1

#guess area
x_start = -5
x_stop = 5
x_steps = 50

y_start = -5
y_stop = 5
y_steps = 50

#tolerances
t_true = .01
t_sim = 1e-3

def disp(betaR, new):
    betaC = betaR[0] + betaR[1]*1j
    k1 = sgn1*cmath.sqrt(betaC*betaC - k0*k0*eps1)
    k2 = sgn2*cmath.sqrt(betaC*betaC - k0*k0*eps2)
    k3 = sgn3*cmath.sqrt(betaC*betaC - k0*k0*eps3)

    if new:
        try:
            out = (k1/eps1 - k2/eps2) * (k1/eps1 - k3/eps3) * cmath.exp(-2*k1*t) - (k1/eps1 + k2/eps2) * (k1/eps1 + k3/eps3)
        except:
            out = 1e20
    else:
        try:
            out = cmath.exp(-2*k1*t) - ((k1/eps1 + k2/eps2) * (k1/eps1 + k3/eps3))/((k1/eps1 - k2/eps2) * (k1/eps1 - k3/eps3))
        except:
            out = 1e20        
    return [out.real, out.imag]

print('Guess area is a rectangle in the complex plane: [%s, %s]x[%si, %si], x_steps = %s, y_steps = %s' % (x_start, x_stop, y_start, y_stop, x_steps, y_steps))
print()
mroots = []
for sgn1, sgn2, sgn3 in product((-1,1), (-1,1), (-1,1)):
    roots = []
    print('Branch:', ('%+d' % sgn1)[0]+('%+d' % sgn2)[0]+('%+d' % sgn3)[0])
    print('Re(beta)        Im(beta)             Abs(new form)   Abs(old form)')
    for x in np.linspace(x_start, x_stop, num=x_steps):
        for y in np.linspace(y_start, y_stop, num=y_steps):
            nrt = root(disp, [x, y], args=(True,), method='hybr')
            if nrt.success:
                new = True
                if norm(disp(nrt.x, True)) > t_true:
                    new = False
                for rt in roots:
                    if norm(rt - nrt.x) < t_sim:
                        new = False
                if new:
                    roots.append(nrt.x)
                    print('%+.12f' % nrt.x[0], '%+.12fi' % nrt.x[1], '   ', '%+.12f' % norm(disp(nrt.x, True)), '%+.12f' % norm(disp(nrt.x, False)))
    print('Total:', len(roots))
    print()
    mroots.append(roots)

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

plt.show()
