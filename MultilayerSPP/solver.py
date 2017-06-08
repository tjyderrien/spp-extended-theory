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
x_steps = 30

y_min = -1e9
y_max = 1e9
y_steps = 30

totalg = x_steps*y_steps

#tolerances

t_zero = 100000
t_blur = 10

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
print('    eps1:', eps1)
print('    eps2:', eps2)
print('    eps3:', eps3)
print('    k0:', k0)
print('    t:', t)

broots = []
for sgn1, sgn2, sgn3 in product((-1,1), (-1,1), (-1,1)): #choose a branch
    roots = []
    separed = []
    num = 0
    print()
    print('Branch: (%s, %s, %s)' % (('%+d' % sgn1)[0], ('%+d' % sgn2)[0], ('%+d' % sgn3)[0]))
    for x in np.linspace(x_min, x_max, num=x_steps):
        num += 1
        print('Tracing the grid: %d%%' % (num/x_steps*100), end = '\r')
        for y in np.linspace(y_min, y_max, num=y_steps): #take an initial guess [x, y] from within the specified grid, x_steps and y_steps determine the grid density
            nrt = root(func, [x, y], method='hybr')
            if nrt.success and norm(nrt.fun) < t_zero:   #if it converges and the value is less than tolerance t_zero, add it to the list of roots
                roots.append(nrt.x)

    ln = len(roots)
    print('Total (converged): %d   ' % ln)

    for i in range(len(roots)):                          #this part should remove duplicate roots (considering tolerance t_blur)
        print('Similarizing roots: %.2f%%' % (i/ln*100), end = '\r')
        for j in range(i + 1, len(roots)):
            if norm(roots[i] - roots[j], 1) < t_blur:
                switch = True
                for sep in separed:
                    if i in sep and j not in sep:
                        separed[separed.index(sep)].append(j)
                        switch = False
                    elif j in sep and i not in sep:
                        separed[separed.index(sep)].append(i)
                        switch = False
                    elif i in sep and j in sep:
                        switch = False
                if switch:
                    separed.append([i, j])

    unique = [np.mean([roots[i] for i in sep], axis = 0) for sep in separed]

    for rt in unique:
        print('% .10e % .10e    % .10e % .10e' % (rt[0], rt[1], func(rt)[0], func(rt)[1]))
    print('Total (similarized): %d' % len(unique))
    broots.append(unique)

#a = cmath.sqrt(eps1*k0*k0)
#plt.plot((a.real, -a.real), (a.imag, -a.imag), 'red')
colors = ['red', 'green', 'blue', 'black', 'magenta', 'cyan', 'lime', 'orangered']
for rts, col in zip(broots, colors):
    if len(rts) > 0:
        plt.scatter(*zip(*rts), c=col)

plt.xlabel('Re')
plt.ylabel('Im')
plt.grid()
plt.savefig('betaSolution.eps')
plt.show()
