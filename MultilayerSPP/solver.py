import cmath
import numpy as np

import matplotlib.pyplot as plt
from scipy.optimize import root
from numpy.linalg import norm
from itertools import product
from libMaterials import Drude

pi = np.pi
prec = 8 #printing precision

#data
wavelength = 800e-9
ne=1e28 #(m^-3) quantity of electrons in conduction band
nu = (1.1e-15)**-1 #collision time between conduction band electrons
meff = 0.18

epsilon = 13.64+0.048j
eps2 = 1.+0.j       #environment
eps3 = 13.64+0.048j #substrate Si (no excitation)
eps1 = Drude(wavelength, ne, eps3, nu, meff)

k0 = 2*pi/wavelength
t = 10e-9 #in meters

#branches
branches = [0, 1, 2, 3, 4, 5, 6, 7] #list of branches you want to use

#0 (-, -, -)
#1 (-, -, +)
#2 (-, +, -)
#3 (-, +, +)
#4 (+, -, -)
#5 (+, -, +)
#6 (+, +, -)
#7 (+, +, +)

#area of initial guesses
#x = real part, y = imaginary part
x_min = -5e3
x_max = 5e3

y_min = 2e9
y_max = 3e9

x_steps = 20
y_steps = 5000

#tolerances
t_blur = 100000

#-----------------------------------------------------------------------------

a = k0*cmath.sqrt(eps1*eps2*(eps1-eps2)/(eps1*eps1 - eps2*eps2))
print('% .5e, % .5e, % .5e, % .5e' % (a.real, a.imag, -a.real, -a.imag))

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

print('Guess area is a rectangle:')
print('[%.2e, %.2e]x[%.2ei, %.2ei], x_steps = %d, y_steps = %d' % (x_min, x_max, y_min, y_max, x_steps, y_steps))
print()
print('Initial data:')
print('    eps1:', eps1)
print('    eps2:', eps2)
print('    eps3:', eps3)
print('    k0:', k0)
print('    t:', t)

broots = []
for sgn1, sgn2, sgn3 in [list(product((-1,1), (-1,1), (-1,1)))[i] for i in branches]:
    roots = []
    unique = []
    num = 0
    print()
    print()
    print('Branch: (%s, %s, %s)' % (('%+d' % sgn1)[0], ('%+d' % sgn2)[0], ('%+d' % sgn3)[0]))
    for x in np.linspace(x_min, x_max, num=x_steps):
        print('Tracing the grid: %d%%' % (num/x_steps*100), end = '\r') 
        for y in np.linspace(y_min, y_max, num=y_steps): #take an initial guess [x, y] from within the specified grid, x_steps and y_steps determine the grid density
            nrt = root(func, [x, y], method='hybr') #uses the MINPACK method: http://www.netlib.org/minpack/
            if nrt.success: #if it converges, add it to the list of roots
                roots.append(nrt.x)
        num += 1

    lns = len(roots)
    print('Total (converged): %d   ' % lns)
    
##    for i in range(len(roots)): #this part should remove duplicate roots (considering tolerance t_blur)
##        print('Assimilating roots: %d%%' % (i/ln*100), end = '\r')
##        for j in range(i + 1, len(roots)):
##            if norm(roots[i] - roots[j]) < t_blur:
##                switch = True
##                for sep in separed:
##                    if i in sep:
##                        if j in sep:
##                            switch = False
##                            break
##                        else:
##                            separed[separed.index(sep)].append(j)
##                            switch = False
##                            break
##                    else:
##                        if j in sep:
##                            separed[separed.index(sep)].append(i)
##                            switch = False
##                            break
##                if switch:
##                    separed.append([i, j])
##
##    unique = [np.mean([roots[i] for i in sep], axis = 0) for sep in separed]

    ln = lns
    while ln > 0: #this part should remove duplicate roots (considering tolerance t_blur)
        print('Assimilating roots: %d%%' % (100 - ln/lns*100), end = '\r')
        center = roots[0] #take the first element of the list and put it in 'center'
        aux = [center]
        aux2 = []
        for rt in roots[1:]: #go thru the rest of roots
            if norm(rt - center) < t_blur: #if the current root 'rt' is sufficiently near center, add it to the new list of similar roots 'aux'
                aux.append(rt)
                center = np.mean(aux, axis = 0) #update 'center' to include the new root, 'center' always lies in the middle
            else:
                aux2.append(rt)
        if len(aux) > 2:
            unique.append(center) #'unique' contains separated roots
        roots = list(aux2) #continue the algorithm with the rest
        ln = len(roots)

    if len(unique) > 0:
        print(' '*23)
        print('Real%s     Imaginary%s│  Abs(F)  %s│  Levenberg–Marquardt' % tuple([' '*prec]*3))
        print('%s┼%s┼%s' % ('─'*(2*prec + 18), '─'*(prec + 10), '─'*(2*prec + 24)))
    for rt in unique:
        nrt = root(func, rt, method='lm')
        print(('%% .%de  %% .%de  │ %% .%de  │  %%s  %% .%de  %% .%de' % tuple([prec]*5)) % (rt[0], rt[1], norm(func(rt)), nrt.success, nrt.x[0], nrt.x[1])) #checking with Levenberg–Marquardt method
    if len(unique) > 0:
        print()
    print('Total (assimilated): %d ' % len(unique))
    #a = k0*cmath.sqrt(eps1*eps2*(eps1-eps2)/(eps1*eps1 - eps2*eps2))
    #print('% .5e, % .5e, % .5e, % .5e' % (a.real, a.imag, -a.real, -a.imag))

    broots.append(unique)

##a = cmath.sqrt(eps1*k0*k0)
##plt.plot((a.real, -a.real), (a.imag, -a.imag), 'red')
colors = ['red', 'green', 'blue', 'black', 'magenta', 'cyan', 'lime', 'orangered']
for rts, col in zip(broots, colors):
    if len(rts) > 0:
        plt.scatter(*zip(*rts), c=col)

plt.xlabel('Re')
plt.ylabel('Im')
plt.grid()
plt.savefig('betaSolution.eps')
plt.show()
