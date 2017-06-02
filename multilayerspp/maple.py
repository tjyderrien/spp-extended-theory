import cmath
import numpy as np

import matplotlib.pyplot as plt
from scipy.optimize import root, fsolve
from scipy.linalg import norm
from itertools import product

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

eps1 = 1 + 1j
eps2 = 3 + 1j
eps3 = 1 + 3j
k0 = 1
t = 2

sgn1 = 1
sgn2 = 1
sgn3 = -1

pnt = [.78068843594832496398, .41062177797714372379]
print('Starting point:', pnt)
print('Branch:', ('%+d' % sgn1)[0]+('%+d' % sgn2)[0]+('%+d' % sgn3)[0])
print()
print(root(disp, pnt, args=(True,)))
print()
print(root(disp, pnt, args=(False,)))
