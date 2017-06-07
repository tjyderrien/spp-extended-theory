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
x_steps = 50

y_min = -1e9
y_max = 1e9
y_steps = 50

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

broots = []
for sgn1, sgn2, sgn3 in product((-1,1), (-1,1), (-1,1)): #choose a branch
    roots = []
    unique = []
    print()
    print('Branch: [%s][%s][%s]' % (('%+d' % sgn1)[0], ('%+d' % sgn2)[0], ('%+d' % sgn3)[0]))
    for x in np.linspace(x_min, x_max, num=x_steps):
        for y in np.linspace(y_min, y_max, num=y_steps):
            #take an initial guess [x, y] from within the specified grid (x_min, x_max, x_steps, y_min, y_max, y_steps), x_steps and y_steps determine the grid density
            nrt = root(func, [x, y], method='hybr') #try to calculate a root for this initial guess
            if nrt.success and norm(nrt.fun) < t_zero: #if it converges and the value is less than tolerance t_zero, add it to the list of roots
                roots.append(nrt.x)

    while len(roots) > 0: #now go thru the list of roots
        center = roots[0] #take the first element of the list and put it in 'center'
        same = [center]
        nrest = []

        for rt in roots[1:]: #go thru the rest of roots
            if norm(rt - center) < t_smear: #if the current root 'rt' is sufficiently near center, add it to the new list of similar roots 'same'
                same.append(rt)
                center = np.mean(same, axis = 0) #update 'center' to include the new root, 'center' always lies in the middle
            else:
                nrest.append(rt)

        unique.append(center) #'unique' contains separated roots
        roots = list(nrest) #continue the algorithm with the rest

    #while loop ended, now all the roots are separated and their centers are inside 'unique'

    for rt in unique:
        print('% .10e  % .10e       % .10e  % .10e' % (rt[0], rt[1], func(rt)[0], func(rt)[1]))
    print('Total: %d' % len(unique))
    broots.append(unique)

colors = ['red', 'green', 'blue', 'black', 'magenta', 'cyan', 'lime', 'orangered']
for rts, col in zip(broots, colors):
    if len(rts) > 0:
        plt.scatter(*zip(*rts), c=col)

plt.xlabel('Re')
plt.ylabel('Im')
plt.grid()
plt.savefig('betaSolution.eps')
plt.show()
